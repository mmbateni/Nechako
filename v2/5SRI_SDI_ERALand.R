##############################################################################
# SRI & SDI CALCULATION - PARALLEL VERSION
#   SRI : Standardized Runoff Index (Shukla & Wood, 2008)
#   SDI : Streamflow Drought Index  (Nalbantis & Tsakiris, 2009)
# Input : ERA5-Land monthly runoff products (00 UTC 24-hour accumulation,
#         monthly mean of daily accumulations), as written by the project
#         downloader into monthly_data_direct/.
# Shared-path version aligned with 1SPI_ERALand.R
##############################################################################

# ---- Libraries ----
rm(list = ls())

# ---- Working directory (framework rule) -------------------------------------
# Run every script from the project root: set the working directory there
# BEFORE running (setwd(), or open the RStudio project). All paths below are
# derived from getwd(); no path is hard-coded.
WD_PATH    <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
UTILS_FILE <- file.path(WD_PATH, "DROUGHT_ANALYSIS_utils.R")
if (!file.exists(UTILS_FILE))
  stop("DROUGHT_ANALYSIS_utils.R not found in the working directory: ", WD_PATH,
       "\nSet the working directory to the project root before running this script.")

library(terra)
library(ncdf4)
library(zoo)
library(writexl)
library(lmomco)
library(parallel)

# ---- Paths / shared utilities ----
source(UTILS_FILE)

data_dir <- file.path(WD_PATH, "monthly_data_direct")
sri_root <- file.path(WD_PATH, "sri_results_seasonal")
sdi_root <- file.path(WD_PATH, "sdi_results")

# ==============================================================================
# USER SETTINGS
# ==============================================================================
# Total runoff (ro = sro + ssro) is the standard SRI/SDI input. Surface and
# sub-surface runoff are optional component diagnostics; remove if not needed.
runoff_inputs <- c(
  runoff             = "runoff_monthly.nc",
  surface_runoff     = "surface_runoff_monthly.nc",
  sub_surface_runoff = "sub_surface_runoff_monthly.nc"
)

sri_timescales    <- c(1, 3, 6, 12, 24, 36)
# SRI distribution: zero-inflated gamma fitted by L-moments (lmomco), with
# Stagge et al. (2015) centre-of-mass zero handling -- identical to SPI and SSPI
# (zero_inflated_gamma_prob() in DROUGHT_ANALYSIS_utils.R).
sri_distribution  <- "gamma (L-moments)"
zero_prob_method  <- "stagge2015"
low_var_threshold <- 0.01              # mm^2, same rule as SPI
min_positive_fit  <- 10                # min. positive values for a parametric fit
min_sample_size   <- 5                 # min. years per calendar month / SDI period
z_clip            <- 4.75
eps_clip          <- 1e-6

sdi_start_month   <- 10                # hydrological year start (10 = October water year)
sdi_periods       <- c(3, 6, 9, 12)    # cumulative months -> k1..k4
sdi_transform     <- "log"             # "log" (log-normal SDI, recommended) or "none"

negative_policy   <- "zero"            # ERA5-Land runoff is signed: "zero" or "na"
use_parallel      <- TRUE

stopifnot(sdi_transform %in% c("log", "none"),
          negative_policy %in% c("zero", "na"),
          sdi_start_month %in% 1:12,
          all(sdi_periods %in% 1:12))

cfg <- list(low_var = low_var_threshold, min_pos = min_positive_fit,
            min_n = min_sample_size, z_clip = z_clip, eps = eps_clip,
            timescales = sri_timescales, sdi_periods = sdi_periods,
            sdi_transform = sdi_transform)

# ==============================================================================
# I/O HELPERS
# ==============================================================================
nc_time_to_date <- function(vals, units) {
  m <- regmatches(units, regexec("^\\s*(seconds|minutes|hours|days)\\s+since\\s+(.+)$", units))[[1]]
  if (length(m) != 3L) stop("Unrecognised NetCDF time units: ", units)
  mult <- c(seconds = 1, minutes = 60, hours = 3600, days = 86400)[[m[2]]]
  origin <- as.POSIXct(substr(trimws(m[3]), 1, 19), tz = "UTC",
                       tryFormats = c("%Y-%m-%d %H:%M:%S", "%Y-%m-%d"))
  as.Date(origin + vals * mult, tz = "UTC")
}

read_runoff_product <- function(path) {
  if (!file.exists(path)) stop("Input file not found: ", path)
  nc <- ncdf4::nc_open(path)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  
  data_vars <- names(nc$var)[vapply(nc$var, function(v) v$ndims == 3L, logical(1))]
  if (length(data_vars) != 1L)
    stop("Expected exactly one 3-D data variable in ", path, "; found: ",
         paste(data_vars, collapse = ", "))
  vname <- data_vars[[1L]]
  
  get_att <- function(var, att) {
    a <- ncdf4::ncatt_get(nc, var, att)
    if (isTRUE(a$hasatt)) trimws(as.character(a$value)) else NA_character_
  }
  
  tdim <- intersect(c("valid_time", "time"), names(nc$dim))
  if (!length(tdim)) stop("No time dimension found in ", path)
  tdim <- tdim[[1L]]
  dates_nc <- nc_time_to_date(nc$dim[[tdim]]$vals, nc$dim[[tdim]]$units)
  
  r <- tryCatch(terra::rast(path, subds = vname), error = function(e) terra::rast(path))
  if (terra::nlyr(r) != length(dates_nc)) stop("Layer/time mismatch in ", path)
  
  list(rast = r, var = vname, units = get_att(vname, "units"),
       long_name = get_att(vname, "long_name"),
       global = ncdf4::ncatt_get(nc, 0), dates = dates_nc)
}

normalise_monthly_dates <- function(d) {
  d <- as.Date(format(d, "%Y-%m-01"))
  mi <- as.integer(format(d, "%Y")) * 12L + as.integer(format(d, "%m"))
  if (anyDuplicated(mi) || any(diff(mi) != 1L))
    stop("Time axis is not a continuous monthly sequence.")
  d
}

# Do not infer units from magnitude: only explicit unit strings are accepted.
runoff_unit_mode <- function(u) {
  if (is.na(u)) return(NA_character_)
  u0 <- gsub("\\s+", " ", tolower(trimws(u)))
  if (u0 %in% c("mm/day", "mm/d", "mm day-1", "mm d-1", "mm per day")) return("mm_per_day")
  if (u0 %in% c("mm/month", "mm month-1", "mm per month"))            return("mm_per_month")
  if (u0 %in% c("m/day", "m/d", "m day-1", "m d-1", "m per day") ||
      grepl("^m (of water equivalent )?per 24[- ]?h(our)?s? accumulation$", u0))
    return("m_per_day")
  NA_character_
}

days_in_month_vec <- function(d) {
  nxt <- seq(d[1L], by = "month", length.out = length(d) + 1L)[-1L]
  as.integer(nxt - d)
}

# Matrix is cells x time; conversion is applied column-wise (unambiguous).
to_mm_month <- function(mat, mode, dim_vec) {
  switch(mode,
         mm_per_day   = sweep(mat, 2, dim_vec, `*`),
         m_per_day    = sweep(mat * 1000, 2, dim_vec, `*`),
         mm_per_month = mat)
}

# ==============================================================================
# STATISTICAL CORE (exported to workers)
# ==============================================================================
roll_sum_right <- function(x, k) {
  if (k == 1) return(x)
  if (length(x) < k) return(rep(NA_real_, length(x)))
  as.numeric(zoo::rollapply(x, width = k, FUN = sum, align = "right",
                            fill = NA_real_, partial = FALSE))
}

clip_prob  <- function(p, eps = 1e-6) pmax(pmin(p, 1 - eps), eps)
gringorten <- function(x) (rank(x, ties.method = "average") - 0.44) / (length(x) + 0.12)

# Variance-aware standardisation by calendar month (same logic as SPI).
# Method codes: 1 = zero-inflated gamma (L-moments, Stagge 2015), 2 = empirical, 3 = zero variance
standardize_monthly <- function(x, month_vec, cfg) {
  n <- length(x); z <- rep(NA_real_, n); meth <- rep(NA_integer_, n)
  for (m in 1:12) {
    idx <- which(month_vec == m & is.finite(x))
    if (length(idx) < cfg$min_n) next
    samp <- x[idx]; sv <- stats::var(samp)
    
    if (!is.finite(sv) || sv < .Machine$double.eps) { z[idx] <- 0; meth[idx] <- 3L; next }
    if (sv < cfg$low_var) {
      z[idx] <- stats::qnorm(clip_prob(gringorten(samp), cfg$eps)); meth[idx] <- 2L; next
    }
    
    p <- zero_inflated_gamma_prob(samp, samp, min_pos = cfg$min_pos)
    if (!is.null(p) && all(is.finite(p))) {
      z[idx] <- stats::qnorm(clip_prob(p, cfg$eps)); meth[idx] <- 1L; next
    }
    z[idx] <- stats::qnorm(clip_prob(gringorten(samp), cfg$eps)); meth[idx] <- 2L
  }
  f <- is.finite(z)
  z[f] <- pmin(pmax(z[f], -cfg$z_clip), cfg$z_clip)
  list(z = z, method = meth)
}

# Hydrological-year index shared by all pixels.
build_hydro_index <- function(dates, start_month, periods) {
  yr  <- as.integer(format(dates, "%Y")); mo <- as.integer(format(dates, "%m"))
  hy  <- yr - as.integer(mo < start_month)
  pos <- ((mo - start_month) %% 12L) + 1L
  hy_all <- sort(unique(hy))
  avail <- matrix(FALSE, length(hy_all), 12L)
  avail[cbind(match(hy, hy_all), pos)] <- TRUE
  keep <- vapply(seq_along(hy_all),
                 function(r) all(avail[r, seq_len(min(periods))]), logical(1))
  years <- hy_all[keep]
  list(row = match(hy, years), pos = pos, years = years, n = length(years),
       start_month = start_month, periods = periods)
}

# SDI: cumulative volumes from the start of the hydrological year.
# Method codes: 1 = parametric (log-normal or plain), 2 = empirical, 3 = zero variance
compute_sdi <- function(x, hyi, cfg) {
  H <- matrix(NA_real_, hyi$n, 12L)
  ok <- !is.na(hyi$row)
  H[cbind(hyi$row[ok], hyi$pos[ok])] <- x[ok]
  np <- length(hyi$periods)
  V <- matrix(NA_real_, hyi$n, np); Z <- V; M <- matrix(NA_integer_, hyi$n, np)
  for (k in seq_len(np)) {
    V[, k] <- rowSums(H[, seq_len(hyi$periods[k]), drop = FALSE])  # NA if any month missing
    idx <- which(is.finite(V[, k]))
    if (length(idx) < cfg$min_n) next
    s <- V[idx, k]; sv <- stats::var(s)
    if (!is.finite(sv) || sv < .Machine$double.eps) { Z[idx, k] <- 0; M[idx, k] <- 3L; next }
    if (cfg$sdi_transform == "log" && all(s > 0)) {
      y <- log(s); Z[idx, k] <- (y - mean(y)) / stats::sd(y); M[idx, k] <- 1L
    } else if (cfg$sdi_transform == "none") {
      Z[idx, k] <- (s - mean(s)) / stats::sd(s); M[idx, k] <- 1L
    } else {
      Z[idx, k] <- stats::qnorm(clip_prob(gringorten(s), cfg$eps)); M[idx, k] <- 2L
    }
  }
  list(V = V, z = Z, method = M)
}

process_pixel <- function(x, month_vec, hyi, cfg) {
  sri <- lapply(cfg$timescales, function(sc)
    standardize_monthly(roll_sum_right(x, sc), month_vec, cfg))
  list(sri = sri, sdi = compute_sdi(x, hyi, cfg))
}

# Worker reads data exported once to each worker's global environment.
worker_fun <- function(i) {
  tryCatch(process_pixel(.wk_mat[i, ], .wk_month, .wk_hyi, .wk_cfg),
           error = function(e) list(error = conditionMessage(e)))
}

# ==============================================================================
# OUTPUT HELPERS
# ==============================================================================
write_index_nc <- function(sub_mat, layer_dates, template, valid_cells, file, varname, longname) {
  r <- terra::rast(template, nlyrs = ncol(sub_mat))
  full <- matrix(NA_real_, terra::ncell(template), ncol(sub_mat))
  full[valid_cells, ] <- sub_mat
  terra::values(r) <- full
  terra::time(r) <- layer_dates
  terra::writeCDF(r, file, varname = varname, longname = longname,
                  unit = "standardized_index", missval = -9999, overwrite = TRUE)
}

plot_method_map <- function(meth, template, valid_cells, basin, file, title, empirical_pct) {
  dom <- apply(meth, 1, function(x) {
    x <- x[is.finite(x)]
    if (!length(x)) NA_real_ else as.numeric(names(which.max(table(x))))
  })
  full <- rep(NA_real_, terra::ncell(template)); full[valid_cells] <- dom
  mr <- terra::rast(template); terra::values(mr) <- full
  labs <- c("Gamma (Stagge zeros)", "Empirical", "Zero-Variance")
  cols <- c("#4575b4", "#d73027", "#91bfdb")
  png(file, width = 1200, height = 800, res = 150); on.exit(dev.off(), add = TRUE)
  plot(mr, col = cols, breaks = c(0.5, 1.5, 2.5, 3.5), legend = FALSE,
       main = title, axes = FALSE, box = FALSE)
  plot(basin, add = TRUE, border = "darkgray", lwd = 1.5)
  cnt <- table(factor(dom, levels = 1:3, labels = labs))
  legend("bottomright", legend = sprintf("%s (%d pixels)", names(cnt), as.integer(cnt)),
         fill = cols, title = "Method", cex = 0.85, bg = "white", bty = "o", border = "gray50")
  mtext(sprintf("Valid pixels: %d | Empirical: %.1f%%", length(valid_cells), empirical_pct),
        side = 1, line = 0.5, cex = 0.75, col = "gray30")
}

classify_sdi <- function(z) {
  as.character(cut(z, breaks = c(-Inf, -2, -1.5, -1, 0, Inf), right = FALSE,
                   labels = c("Extreme", "Severe", "Moderate", "Mild", "Non-drought")))
}

period_end_dates <- function(years, start_month, P) {
  em <- start_month + P - 1L
  ey <- years + (em - 1L) %/% 12L
  em <- ((em - 1L) %% 12L) + 1L
  as.Date(sprintf("%04d-%02d-01", ey, em))
}

pct <- function(m, code) 100 * sum(m == code, na.rm = TRUE) / max(1, sum(!is.na(m)))

# ==============================================================================
# CANONICAL GRID / BASIN MASK (identical to the SPI mask)
# ==============================================================================
cat("\n============================================================\n")
cat("SRI / SDI CALCULATION (PARALLEL MODE)\n")
cat("============================================================\n")

# Whole basin = dissolved union of the 7 authoritative subbasins
# (identical construction to 3SPEI_ERALand.R).
basin <- nechako_subbasin_union(load_nechako_subbasins())
cat("Basin boundary loaded (union of 7 subbasins)\n")

# Rebuild the canonical mask from the same precipitation grid and source label
# used by 1SPI_ERALand.R so that SPI, SRI and SDI share the same basin pixels.
# The runoff grids must match this geometry exactly; no resampling is done.
precip_ref_file <- file.path(data_dir, "total_precipitation_monthly.nc")
if (file.exists(precip_ref_file)) {
  grid_ref   <- terra::rast(precip_ref_file)[[1]]
  grid_label <- "1SPI_ERALand.R canonical native ERA5-Land precipitation grid"
} else {
  warning("Precipitation reference grid not found; using the first runoff grid as canonical grid.")
  grid_ref   <- terra::rast(file.path(data_dir, runoff_inputs[[1]]))[[1]]
  grid_label <- "5SRI_SDI_ERALand.R native ERA5-Land runoff grid"
}
basin <- terra::project(basin, terra::crs(grid_ref))
common_grid <- save_common_analysis_grid(grid_ref, basin = basin, source = grid_label)
common_mask <- terra::unwrap(common_grid$basin_mask)
mask_vals   <- terra::values(common_mask, mat = FALSE)
basin_cells <- which(!is.na(mask_vals) & mask_vals != 0)
if (!length(basin_cells)) stop("Canonical basin mask contains zero cells.")

# Fractional-coverage basin weights (shared framework rule): w = cell area x
# fractional coverage; basin means need >= BASIN_MIN_VALID_AREA_COVERAGE (0.20)
# of the total basin weight to be reported.
sri_bw <- basin_fractional_weights(grid_ref, basin)
if (!identical(which(sri_bw$mask), basin_cells))
  stop("Fractional-weight basin mask differs from the canonical SRI/SDI basin mask.")
cat(sprintf("Basin weights: %d cells (%d partial) | weighted area %.1f km2 | min valid coverage %.2f\n",
            sri_bw$n_cells, sri_bw$n_partial, sri_bw$area_m2 / 1e6, BASIN_MIN_VALID_AREA_COVERAGE))
cat(sprintf("Canonical basin cells: %d\n", length(basin_cells)))

# ==============================================================================
# MAIN PROCESSING FUNCTION (one runoff product)
# ==============================================================================
process_runoff_product <- function(label, file_name, cl) {
  cat("\n============================================================\n")
  cat(sprintf("RUNOFF PRODUCT: %s (%s)\n", label, file_name))
  cat("============================================================\n")
  
  sri_dir <- file.path(sri_root, label); dir.create(sri_dir, recursive = TRUE, showWarnings = FALSE)
  sdi_dir <- file.path(sdi_root, label); dir.create(sdi_dir, recursive = TRUE, showWarnings = FALSE)
  
  prod      <- read_runoff_product(file.path(data_dir, file_name))
  dates     <- normalise_monthly_dates(prod$dates)
  month_vec <- as.integer(format(dates, "%m"))
  template  <- prod$rast[[1]]
  
  # ---- Unit contract check ----
  mode <- runoff_unit_mode(prod$units)
  cat("\n===== RUNOFF UNIT CONTRACT CHECK =====\n")
  cat(sprintf("  NetCDF variable : %s\n", prod$var))
  cat(sprintf("  Long name       : %s\n", prod$long_name))
  cat(sprintf("  Stored units    : %s\n", ifelse(is.na(prod$units), "<missing>", prod$units)))
  for (a in c("accumulation_interpretation", "downstream_conversion"))
    if (!is.null(prod$global[[a]])) cat(sprintf("  %s : %s\n", a, prod$global[[a]]))
  if (is.na(mode))
    stop("Unsupported or missing runoff units in ", file_name, ": ", prod$units)
  cat(sprintf("  Conversion mode : %s -> mm/month (actual calendar days)\n", mode))
  
  if (!terra::compareGeom(prod$rast, common_mask, lyrs = FALSE, stopOnError = FALSE))
    stop("Grid of ", file_name, " does not match the canonical SPI grid. ",
         "Check extent/resolution of the downloaded product.")
  
  # ---- Extract basin matrix, convert units ----
  vals <- terra::values(prod$rast, mat = TRUE)
  vals[!is.finite(vals) | vals < -1e20] <- NA_real_        # fill-value guard
  vals <- vals[basin_cells, , drop = FALSE]
  vals <- to_mm_month(vals, mode, days_in_month_vec(dates))
  
  # ERA5-Land runoff is signed; tiny negatives are numerical artefacts.
  neg <- is.finite(vals) & vals < 0
  n_neg <- sum(neg); min_neg <- if (n_neg) min(vals[neg]) else 0
  if (n_neg) vals[neg] <- if (negative_policy == "zero") 0 else NA_real_
  cat(sprintf("  Negative values : %d (min %.4f mm/month) -> set to %s\n",
              n_neg, min_neg, negative_policy))
  if (min_neg < -1) warning(sprintf("[%s] negative runoff below -1 mm/month found (%.3f).", label, min_neg))
  
  keep        <- rowSums(is.finite(vals)) > 0
  valid_cells <- basin_cells[keep]
  runoff_mat  <- vals[keep, , drop = FALSE]
  n_pix <- nrow(runoff_mat)
  if (!n_pix) stop("No valid basin pixels for ", label)
  cat(sprintf("  Valid basin pixels: %d of %d canonical cells (%d all-NA dropped)\n",
              n_pix, length(basin_cells), sum(!keep)))
  
  # ---- Basin series: fractional-coverage area weights, 0.20 minimum coverage ----
  # Coverage is measured against the whole-basin weight, so cells dropped above
  # (all-NA) count as missing area rather than shrinking the basin.
  basin_depth_mm <- basin_weighted_mean(runoff_mat, sri_bw, cells = valid_cells)
  basin_cov_frac <- basin_valid_coverage(runoff_mat, sri_bw, cells = valid_cells)
  basin_area_m2  <- sri_bw$area_m2          # coverage-weighted basin area
  basin_vol_hm3  <- basin_depth_mm / 1000 * basin_area_m2 / 1e6
  clim <- tapply(basin_depth_mm, month_vec, mean, na.rm = TRUE)
  cat(sprintf("  Basin area: %.1f km2 | mean annual runoff: %.1f mm\n",
              basin_area_m2 / 1e6, sum(clim)))
  
  # ---- Hydrological-year index ----
  hyi <- build_hydro_index(dates, sdi_start_month, sdi_periods)
  hy_lab <- if (sdi_start_month == 1L) as.character(hyi$years) else paste0(hyi$years, "-", hyi$years + 1L)
  period_lab <- vapply(sdi_periods, function(P)
    paste0(month.abb[sdi_start_month], "-", month.abb[((sdi_start_month - 1L + P - 1L) %% 12L) + 1L]), "")
  
  # ---- Pixel-wise SRI + SDI ----
  assign(".wk_mat", runoff_mat, envir = globalenv())
  assign(".wk_month", month_vec, envir = globalenv())
  assign(".wk_hyi", hyi, envir = globalenv())
  assign(".wk_cfg", cfg, envir = globalenv())
  
  cat(sprintf("\nProcessing %d pixels (SRI-%s + SDI k1-k%d)... ",
              n_pix, paste(sri_timescales, collapse = "/"), length(sdi_periods)))
  t0 <- Sys.time()
  if (!is.null(cl)) {
    clusterExport(cl, c(".wk_mat", ".wk_month", ".wk_hyi", ".wk_cfg"), envir = globalenv())
    res <- parLapply(cl, seq_len(n_pix), worker_fun)
  } else {
    res <- lapply(seq_len(n_pix), worker_fun)
  }
  cat(sprintf("done (%.1f min)\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  
  errs <- vapply(res, function(p) !is.null(p$error), logical(1))
  if (any(errs))
    stop(sprintf("[%s] %d/%d pixels failed. First worker error: %s",
                 label, sum(errs), n_pix, res[[which(errs)[1L]]]$error))
  
  sri_z <- lapply(seq_along(sri_timescales), function(s) do.call(rbind, lapply(res, function(p) p$sri[[s]]$z)))
  sri_m <- lapply(seq_along(sri_timescales), function(s) do.call(rbind, lapply(res, function(p) p$sri[[s]]$method)))
  sdi_z <- lapply(seq_along(sdi_periods), function(k) do.call(rbind, lapply(res, function(p) p$sdi$z[, k])))
  sdi_m <- lapply(seq_along(sdi_periods), function(k) do.call(rbind, lapply(res, function(p) p$sdi$method[, k])))
  rm(res); gc()
  
  basin_sri <- lapply(sri_timescales, function(sc)
    standardize_monthly(roll_sum_right(basin_depth_mm, sc), month_vec, cfg))
  basin_sdi <- compute_sdi(basin_vol_hm3, hyi, cfg)
  
  coords <- terra::xyFromCell(template, valid_cells)
  
  # ---- SRI outputs ----
  for (s in seq_along(sri_timescales)) {
    sc <- sri_timescales[s]; z <- sri_z[[s]]; meth <- sri_m[[s]]
    cat(sprintf("  -> SRI-%d outputs (CSV, Excel, NetCDF, map)\n", sc))
    excel_data <- list()
    for (m in 1:12) {
      idx_m <- which(month_vec == m); if (!length(idx_m)) next
      df <- data.frame(lon = coords[, 1], lat = coords[, 2], z[, idx_m, drop = FALSE])
      names(df)[-(1:2)] <- format(dates[idx_m], "%Y")
      write.csv(df, file.path(sri_dir, sprintf("sri_%02d_month%02d_%s.csv", sc, m, month.abb[m])),
                row.names = FALSE)
      excel_data[[month.abb[m]]] <- df
      write_index_nc(z[, idx_m, drop = FALSE], dates[idx_m], template, valid_cells,
                     file.path(sri_dir, sprintf("sri_%02d_month%02d_%s.nc", sc, m, month.abb[m])),
                     "sri", sprintf("Standardized Runoff Index (SRI-%d), ERA5-Land %s", sc, label))
    }
    xlsx_file <- file.path(sri_dir, sprintf("sri_%02d_all_months.xlsx", sc))
    if (file.exists(xlsx_file)) file.remove(xlsx_file)
    write_xlsx(excel_data, xlsx_file)
    plot_method_map(meth, template, valid_cells, basin,
                    file.path(sri_dir, sprintf("sri_%02d_distribution_map.png", sc)),
                    sprintf("SRI-%d (%s)", sc, label), pct(meth, 2))
  }
  
  basin_sri_df <- data.frame(date = dates, runoff_mm_month = basin_depth_mm, volume_hm3 = basin_vol_hm3,
                             valid_area_coverage_fraction = basin_cov_frac)
  for (s in seq_along(sri_timescales))
    basin_sri_df[[sprintf("SRI_%d", sri_timescales[s])]] <- basin_sri[[s]]$z
  write.csv(basin_sri_df, file.path(sri_dir, "sri_basin_timeseries.csv"), row.names = FALSE)
  
  png(file.path(sri_dir, "sri_basin_timeseries.png"), width = 1800, height = 1400, res = 150)
  op <- par(mfrow = c(ceiling(length(sri_timescales) / 2), 2), mar = c(3, 4, 2.5, 1))
  for (s in seq_along(sri_timescales)) {
    z <- basin_sri[[s]]$z
    plot(dates, z, type = "h", col = ifelse(!is.na(z) & z < 0, "#d73027", "#4575b4"),
         ylim = range(c(-3, 3, z), na.rm = TRUE), xlab = "", ylab = "SRI",
         main = sprintf("Basin SRI-%d (%s)", sri_timescales[s], label))
    abline(h = c(-1, -1.5, -2), lty = 2, col = "gray50"); abline(h = 0)
  }
  par(op); dev.off()
  
  # ---- SDI outputs ----
  cat("  -> SDI outputs (CSV, Excel, NetCDF, basin plot)\n")
  sdi_excel <- list()
  for (k in seq_along(sdi_periods)) {
    z <- sdi_z[[k]]
    df <- data.frame(lon = coords[, 1], lat = coords[, 2], z)
    names(df)[-(1:2)] <- hy_lab
    write.csv(df, file.path(sdi_dir, sprintf("sdi_k%d_%s.csv", k, period_lab[k])), row.names = FALSE)
    sdi_excel[[sprintf("k%d_%s", k, period_lab[k])]] <- df
    has <- colSums(is.finite(z)) > 0
    if (any(has))
      write_index_nc(z[, has, drop = FALSE],
                     period_end_dates(hyi$years, sdi_start_month, sdi_periods[k])[has],
                     template, valid_cells,
                     file.path(sdi_dir, sprintf("sdi_k%d_%s.nc", k, period_lab[k])), "sdi",
                     sprintf("Streamflow Drought Index k%d (%s), ERA5-Land %s", k, period_lab[k], label))
  }
  
  basin_sdi_df <- data.frame(hydro_year = hy_lab, start_year = hyi$years)
  for (k in seq_along(sdi_periods)) {
    tag <- sprintf("k%d_%s", k, period_lab[k])
    basin_sdi_df[[paste0("volume_hm3_", tag)]] <- basin_sdi$V[, k]
    basin_sdi_df[[paste0("SDI_", tag)]]        <- basin_sdi$z[, k]
    basin_sdi_df[[paste0("class_", tag)]]      <- classify_sdi(basin_sdi$z[, k])
  }
  write.csv(basin_sdi_df, file.path(sdi_dir, "sdi_basin_timeseries.csv"), row.names = FALSE)
  sdi_excel[["basin"]] <- basin_sdi_df
  xlsx_sdi <- file.path(sdi_dir, "sdi_all_periods.xlsx")
  if (file.exists(xlsx_sdi)) file.remove(xlsx_sdi)
  write_xlsx(sdi_excel, xlsx_sdi)
  
  png(file.path(sdi_dir, "sdi_basin_timeseries.png"), width = 1800, height = 1200, res = 150)
  op <- par(mfrow = c(2, 2), mar = c(4.5, 4, 2.5, 1))
  for (k in seq_along(sdi_periods)) {
    z <- basin_sdi$z[, k]
    barplot(z, names.arg = hyi$years, col = ifelse(is.na(z) | z >= 0, "#4575b4", "#d73027"),
            border = NA, las = 2, cex.names = 0.5, ylab = "SDI",
            main = sprintf("Basin SDI k%d (%s) - %s", k, period_lab[k], label))
    abline(h = c(-1, -1.5, -2), lty = 2, col = "gray50"); abline(h = 0)
  }
  par(op); dev.off()
  
  # ---- Summary file ----
  S <- character(0)
  add <- function(fmt, ...) S <<- c(S, sprintf(fmt, ...))
  add("============================================================")
  add("SRI / SDI SUMMARY - ERA5-Land %s", label)
  add("============================================================")
  add("Calculation date    : %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  add("Input file          : %s (variable '%s')", file_name, prod$var)
  add("Stored units        : %s -> %s -> mm/month", prod$units, mode)
  add("Negative values     : %d (min %.4f mm/month), policy = %s", n_neg, min_neg, negative_policy)
  add("Time period         : %s to %s (%d months)", format(min(dates)), format(max(dates)), length(dates))
  add("Grid                : %d x %d | valid basin pixels %d", ncol(template), nrow(template), n_pix)
  add("Basin area          : %.1f km2 (cell area x fractional coverage, 7-subbasin union)", basin_area_m2 / 1e6)
  add("Basin mean rule     : fractional-coverage area weights, min valid coverage %.2f", BASIN_MIN_VALID_AREA_COVERAGE)
  add("SRI distribution    : %s | zero handling: %s", sri_distribution, zero_prob_method)
  add("Basin runoff climatology (mm/month): %s",
      paste(sprintf("%s=%.1f", month.abb, clim), collapse = " "))
  
  for (s in seq_along(sri_timescales)) {
    z <- sri_z[[s]]; m <- sri_m[[s]]
    add(""); add("---------- SRI-%d ----------", sri_timescales[s])
    add("NA rate (basin pixels) : %.3f%%", 100 * mean(is.na(z)))
    add("Methods: parametric=%.1f%%  empirical=%.1f%%  zero-var=%.1f%%", pct(m, 1), pct(m, 2), pct(m, 3))
    add("Pixel SRI mean / sd    : %.3f / %.3f", mean(z, na.rm = TRUE), stats::sd(as.vector(z), na.rm = TRUE))
    add("Clipped at +/-%.2f     : %d dry, %d wet", z_clip,
        sum(z <= -z_clip, na.rm = TRUE), sum(z >= z_clip, na.rm = TRUE))
    add("  %-5s %9s %7s %9s %7s %9s %7s", "Month", "Case1(n)", "C1(%)", "Case2(n)", "C2(%)", "Case3(n)", "C3(%)")
    for (mm in 1:12) {
      v <- as.integer(m[, month_vec == mm]); v <- v[!is.na(v)]; nt <- length(v)
      if (!nt) next
      add("  %-5s %9d %6.1f%% %9d %6.1f%% %9d %6.1f%%", month.abb[mm],
          sum(v == 1L), 100 * sum(v == 1L) / nt, sum(v == 2L), 100 * sum(v == 2L) / nt,
          sum(v == 3L), 100 * sum(v == 3L) / nt)
    }
    add("Class frequency (all pixels, full record):")
    add("  W3 extremely wet  (SRI >= 2.0)        : %.1f%%", 100 * mean(z >= 2, na.rm = TRUE))
    add("  W2 severely wet   (1.5 <= SRI < 2.0)  : %.1f%%", 100 * mean(z >= 1.5 & z < 2, na.rm = TRUE))
    add("  W1 moderately wet (1.0 <= SRI < 1.5)  : %.1f%%", 100 * mean(z >= 1 & z < 1.5, na.rm = TRUE))
    add("  N0 normal         (-1 < SRI < 1)      : %.1f%%", 100 * mean(z > -1 & z < 1, na.rm = TRUE))
    add("  D1 moderately dry (-1.5 < SRI <= -1)  : %.1f%%", 100 * mean(z <= -1 & z > -1.5, na.rm = TRUE))
    add("  D2 severely dry   (-2 < SRI <= -1.5)  : %.1f%%", 100 * mean(z <= -1.5 & z > -2, na.rm = TRUE))
    add("  D3 extremely dry  (SRI <= -2.0)       : %.1f%%", 100 * mean(z <= -2, na.rm = TRUE))
  }
  
  add(""); add("---------- SDI (hydrological year starts %s, transform = %s) ----------",
               month.abb[sdi_start_month], sdi_transform)
  for (k in seq_along(sdi_periods)) {
    z <- sdi_z[[k]]; m <- sdi_m[[k]]
    add("SDI k%d (%s): hydro-years=%d | parametric=%.1f%% empirical=%.1f%% zero-var=%.1f%%",
        k, period_lab[k], sum(colSums(is.finite(z)) > 0), pct(m, 1), pct(m, 2), pct(m, 3))
    tab <- table(factor(classify_sdi(basin_sdi$z[, k]),
                        levels = c("Non-drought", "Mild", "Moderate", "Severe", "Extreme")))
    add("  Basin classes: %s", paste(sprintf("%s=%d", names(tab), as.integer(tab)), collapse = ", "))
    w <- which.min(basin_sdi$z[, k])
    if (length(w)) add("  Driest basin hydro-year: %s (SDI = %.2f, volume = %.0f hm3)",
                       hy_lab[w], basin_sdi$z[w, k], basin_sdi$V[w, k])
  }
  
  summ_name <- sprintf("sri_sdi_%s_summary.txt", label)
  writeLines(S, file.path(sri_dir, summ_name))
  writeLines(S, file.path(sdi_dir, summ_name))
  cat(sprintf("  Summary written: %s\n", summ_name))
  invisible(TRUE)
}

# ==============================================================================
# RUN
# ==============================================================================
cl <- NULL
if (use_parallel) {
  n_cores <- max(1L, detectCores() - 1L)
  if (n_cores > 1L) {
    cat(sprintf("\nStarting parallel cluster with %d cores\n", n_cores))
    dir.create(tempdir(), recursive = TRUE, showWarnings = FALSE)  # Windows temp-dir fix
    cl <- makeCluster(n_cores)
    clusterEvalQ(cl, { library(zoo); library(lmomco); NULL })
    clusterExport(cl, c("roll_sum_right", "clip_prob", "gringorten",
                        "zero_inflated_gamma_prob", "stagge_zero_probabilities",
                        "standardize_monthly", "compute_sdi", "process_pixel", "worker_fun"),
                  envir = globalenv())
  }
}

tryCatch({
  for (lab in names(runoff_inputs)) process_runoff_product(lab, runoff_inputs[[lab]], cl)
}, finally = {
  if (!is.null(cl)) { stopCluster(cl); cat("\nParallel cluster stopped\n") }
  rm(list = intersect(c(".wk_mat", ".wk_month", ".wk_hyi", ".wk_cfg"), ls(globalenv(), all.names = TRUE)),
     envir = globalenv())
})

cat("\n============================================================\n")
cat("SRI / SDI CALCULATION COMPLETE\n")
cat("============================================================\n")
cat(sprintf("SRI output root: %s\n", normalizePath(sri_root)))
cat(sprintf("SDI output root: %s\n", normalizePath(sdi_root)))