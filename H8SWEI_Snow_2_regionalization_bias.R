##############################################
# 02_regionalization_bias.R
# STAGE 2 of 4 — Function definitions, RFA regionalization machinery,
#                 and basin-wide bias correction
#
# WHAT THIS SCRIPT DOES:
#   - Defines every helper/RFA/bias-correction function used downstream
#   - Runs the basin-wide bias correction against station SWE obs (if
#     BIAS_CORRECT is enabled) — this includes the slow LOSO-CV step
#   - Runs the snow-cover-timing validation and other QA checks
#   - Runs the RFA construction block (regionalization, discordancy,
#     heterogeneity, pooled regional gamma fits) if REGIONALIZE_SSPI is TRUE
#
# This is often the SLOWEST stage before the main parallel loop (LOSO-CV and
# the heterogeneity Monte Carlo are themselves parallelized). Isolating it
# means a bug/crash here never forces you to reload and re-crop the rasters
# from script 01.
#
# INPUT:  stage_checkpoints/stage1.rds  (from 01_setup_and_data.R)
# OUTPUT: stage_checkpoints/stage2.rds  (everything script 03 needs)
##############################################

STAGE_DIR <- "stage_checkpoints"
stage1_path <- file.path(STAGE_DIR, "stage1.rds")
if (!file.exists(stage1_path)) {
  stop(sprintf(
    "Cannot find %s. Run 01_setup_and_data.R first (it must complete and write this checkpoint).",
    stage1_path
  ), call. = FALSE)
}

cat("============================================================\n")
cat("STAGE 2 - loading checkpoint from script 01\n")
cat("============================================================\n")
stage1_data <- readRDS(stage1_path)
list2env(stage1_data, envir = environment())   # dates_swe, dates, etc. now exist
rm(stage1_data)

# ---- Rehydrate SpatRaster/SpatVector objects from terra's own files ----
# Done AFTER list2env() so these are the last word - nothing here gets
# clobbered by the (deliberately swe/scf_stack/basin-free) RDS payload.
swe <- terra::rast(file.path(STAGE_DIR, "stage1_swe.tif"))
terra::time(swe) <- dates_swe          # GeoTIFF doesn't store the time dimension
scf_stack_path <- file.path(STAGE_DIR, "stage1_scf_stack.tif")
scf_stack <- if (file.exists(scf_stack_path)) {
  s <- terra::rast(scf_stack_path)
  terra::time(s) <- dates_scf
  s
} else NULL
basin <- terra::vect(file.path(STAGE_DIR, "stage1_basin.gpkg"))

cat(sprintf("[OK] Restored %d object(s) from %s\n", length(ls()), stage1_path))
# ==============================================================================
#   Helper Functions
# ==============================================================================
clip_prob <- function(p, eps = 1e-6) pmax(pmin(p, 1 - eps), eps)
## Converts H0 into a continuous pooling weight instead of a hard cliff at
## H0_REGION_FLAG_THRESHOLD. weight=1 -> fully pooled (regional gamma only),
## weight=0 -> fully local. Only trustworthy once HETEROGENEITY_NSIM is large
## (Fix #1) - blending a noisy H0 is still noisy, just less abruptly so.
h0_shrink_weight <- function(H0, threshold = H0_REGION_FLAG_THRESHOLD, band = 0.2) {
  if (!is.finite(H0)) return(0)
  lo <- threshold - band; hi <- threshold + band
  if (H0 <= lo) return(1)
  if (H0 >= hi) return(0)
  1 - (H0 - lo) / (hi - lo)
}
# ---- Progress-reporting helper ----
# Formats a duration in seconds as H:MM:SS (or M:SS under an hour) for the
# elapsed-time / ETA messages added throughout the slower sections of this
# script (LOSO-CV, regionalized bias-correction fitting, the k-fold/temporal-
# transfer validation protocol, etc.) so a long-running step reports how long
# it's been going and roughly how much longer it will take, instead of
# leaving the console silent.
format_hms <- function(seconds) {
  if (!is.finite(seconds) || seconds < 0) return("unknown")
  seconds <- round(seconds)
  h <- seconds %/% 3600
  m <- (seconds %% 3600) %/% 60
  s <- seconds %% 60
  if (h > 0) sprintf("%d:%02d:%02d", h, m, s) else sprintf("%d:%02d", m, s)
}
## ---- Effective sample size for spatially autocorrelated pixel pooling ----
## Pixels in a region aren't independent "sites" like gauged catchments -
## they're adjacent grid cells sharing weather. Raw pooled pixel-month count
## overstates the effective degrees of freedom behind MIN_REGIONAL_GAMMA_SAMPLE
## and the GoF/H0 diagnostics. Approximated via the average pairwise
## correlation of DE-SEASONALIZED monthly SWE anomalies within the region -
## one cor() call per region (17 regions), not per pixel-month.
compute_regional_effective_n <- function(swe_matrix, region_id, r, dates_vec, ref_idx) {
  pix_r <- which(region_id == r)
  if (length(pix_r) < 2) return(list(avg_rho = 0, n_pix = length(pix_r), n_eff_multiplier = 1))
  
  mon <- as.integer(format(dates_vec, "%m"))
  sub <- swe_matrix[pix_r, ref_idx, drop = FALSE]
  
  anom <- sub
  for (m in 1:12) {
    cc <- which(mon[ref_idx] == m)
    if (length(cc) < 2) next
    anom[, cc] <- sub[, cc] - rowMeans(sub[, cc, drop = FALSE], na.rm = TRUE)
  }
  
  keep <- rowSums(is.finite(anom)) >= 10
  anom <- anom[keep, , drop = FALSE]
  n_pix_eff <- nrow(anom)
  if (n_pix_eff < 2) return(list(avg_rho = 0, n_pix = length(pix_r), n_eff_multiplier = 1))
  
  cm <- suppressWarnings(stats::cor(t(anom), use = "pairwise.complete.obs"))
  cm[!is.finite(cm)] <- NA_real_
  avg_rho <- mean(cm[upper.tri(cm)], na.rm = TRUE)
  if (!is.finite(avg_rho)) avg_rho <- 0
  avg_rho <- max(0, avg_rho)
  
  # Design-effect deflation: N_eff = K / (1 + (K-1)*rho)
  list(avg_rho = avg_rho, n_pix = n_pix_eff,
       n_eff_multiplier = n_pix_eff / (1 + (n_pix_eff - 1) * avg_rho))
}
# ============================================================================
#   OPTIONAL NETWORK / SAR HELPERS
# ============================================================================
canon_name <- function(x) gsub("[^a-z0-9]+", "", tolower(x))

pick_first_col <- function(df, candidates) {
  cn <- canon_name(names(df))
  for (cand in candidates) {
    idx <- match(canon_name(cand), cn)
    if (!is.na(idx)) return(names(df)[idx])
  }
  NULL
}

standardize_station_obs <- function(df, source_label = "station") {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  c_station <- pick_first_col(df, c("station_id", "station", "site_id", "id"))
  c_lon     <- pick_first_col(df, c("lon", "longitude", "x"))
  c_lat     <- pick_first_col(df, c("lat", "latitude", "y"))
  c_date    <- pick_first_col(df, c("date", "datetime", "timestamp", "start_of_interval_utc", "startofintervalutc"))
  c_swe     <- pick_first_col(df, c("swe_mm", "swe", "snow_water_equivalent_mm", "sw_daily_average_mm", "swdailyaveragemm"))
  if (any(vapply(list(c_station, c_lon, c_lat, c_date, c_swe), is.null, logical(1)))) {
    return(NULL)
  }
  
  out <- data.frame(
    station_id = as.character(df[[c_station]]),
    lon = suppressWarnings(as.numeric(df[[c_lon]])),
    lat = suppressWarnings(as.numeric(df[[c_lat]])),
    date = as.character(df[[c_date]]),
    swe_mm = suppressWarnings(as.numeric(df[[c_swe]])),
    stringsAsFactors = FALSE
  )
  
  if (!inherits(out$date, "Date")) {
    n_date_rows <- length(out$date)
    cat(sprintf("    Parsing %d date value(s) for '%s'... ", n_date_rows, source_label))
    t0_date <- Sys.time()
    
    # Try each format on a SMALL SAMPLE first (cheap) to pick the right one,
    # instead of running as.Date() on the full ~1.4M-row vector up to 4 times
    # in sequence. Falls back to the original full-vector behaviour if none
    # of the sampled formats look promising.
    sample_idx <- sample.int(n_date_rows, min(200, n_date_rows))
    sample_vals <- out$date[sample_idx]
    
    formats_to_try <- list(
      list(fmt = NULL,          label = "ISO/default"),
      list(fmt = "%m/%d/%Y",    label = "%m/%d/%Y"),
      list(fmt = "%Y-%m-%d",    label = "%Y-%m-%d"),
      list(fmt = "%d/%m/%Y",    label = "%d/%m/%Y")
    )
    
    chosen_fmt <- NULL
    for (f in formats_to_try) {
      d_sample <- if (is.null(f$fmt)) as.Date(sample_vals) else as.Date(sample_vals, format = f$fmt)
      if (!all(is.na(d_sample))) {
        chosen_fmt <- f
        break
      }
    }
    
    if (is.null(chosen_fmt)) {
      # Nothing worked even on the sample - fall back to the original
      # exhaustive full-vector attempts so behaviour is unchanged, just slow.
      d <- as.Date(out$date)
      if (all(is.na(d))) d <- as.Date(out$date, format = "%m/%d/%Y")
      if (all(is.na(d))) d <- as.Date(out$date, format = "%Y-%m-%d")
      if (all(is.na(d))) d <- as.Date(out$date, format = "%d/%m/%Y")
      out$date <- d
    } else {
      d <- if (is.null(chosen_fmt$fmt)) as.Date(out$date) else as.Date(out$date, format = chosen_fmt$fmt)
      out$date <- d
    }
    
    elapsed_date <- as.numeric(difftime(Sys.time(), t0_date, units = "secs"))
    n_parsed <- sum(!is.na(out$date))
    cat(sprintf("done in %.1fs (%s format, %d/%d parsed)\n",
                elapsed_date,
                if (is.null(chosen_fmt)) "none matched" else chosen_fmt$label,
                n_parsed, n_date_rows))
  }
  
  out <- out[is.finite(out$lon) & is.finite(out$lat) & is.finite(out$swe_mm) & !is.na(out$date), ]
  if (nrow(out) == 0) return(NULL)
  out$station_id <- ifelse(is.na(out$station_id) | out$station_id == "", source_label, out$station_id)
  out
}

merge_station_networks <- function(asws_path, canswe_path = NULL, swe_template = NULL,
                                   require_canswe = FALSE) {
  files <- c()
  file_sources <- c()
  if (!is.null(asws_path) && file.exists(asws_path)) {
    files <- c(files, asws_path); file_sources <- c(file_sources, "primary")
  }
  if (!is.null(canswe_path) && file.exists(canswe_path)) {
    files <- c(files, canswe_path); file_sources <- c(file_sources, "canswe_direct")
  }
  if (length(files) == 0) return(NULL)
  
  read_one <- function(path, tag) {
    file_mb <- file.info(path)$size / 1e6
    cat(sprintf("  Reading %s (%.1f MB)... ", path, file_mb))
    t0 <- Sys.time()
    
    raw <- tryCatch({
      if (requireNamespace("data.table", quietly = TRUE)) {
        as.data.frame(data.table::fread(path, showProgress = FALSE))
      } else {
        cat("\n    [NOTE] data.table not installed - falling back to read.csv() (much slower for large files).\n    ")
        read.csv(path, stringsAsFactors = FALSE)
      }
    }, error = function(e) {
      cat(sprintf("\n    [ERROR] Failed to read '%s': %s\n", path, conditionMessage(e)))
      NULL
    })
    
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (is.null(raw)) {
      cat(sprintf("failed after %.1fs\n", elapsed))
      return(NULL)
    }
    cat(sprintf("done in %.1fs (%d rows)\n", elapsed, nrow(raw)))
    
    std <- standardize_station_obs(raw, source_label = basename(path))
    if (is.null(std)) {
      log_pipeline_fallback(sprintf("Could not standardize station file '%s' into station_id/lon/lat/date/swe_mm.", path))
      return(NULL)
    }
    # Tag provenance explicitly rather than relying on value coincidence.
    # "primary" file may itself carry CanSWE_* rows (merged in by the download
    # script); flag those on content, not on which file they came from, so the
    # is_canswe flag survives whichever file happens to contain them.
    std$is_canswe <- grepl("^CanSWE_", std$station_id) | tag == "canswe_direct"
    std
  }  
  obs_list <- Map(read_one, files, file_sources)
  obs_list <- obs_list[!vapply(obs_list, is.null, logical(1))]
  if (length(obs_list) == 0) return(NULL)
  
  out <- do.call(rbind, obs_list)
  out <- out[complete.cases(out[c("station_id", "lon", "lat", "date", "swe_mm")]), ]
  
  # ---- Fix #1 (corrected): drop stations that can NEVER contribute a ----
  # ---- matched pair anywhere downstream, using the SAME criterion every ----
  # ---- fit function already uses - not "inside the basin polygon". ----
  # Tier 0/1/2/3, LOSO-CV, and the validation protocol all discover a
  # station is unusable the same way: cells(swe_template, pts) returns NA
  # (the station falls outside the raster's extent), or the cell it DOES
  # resolve to is NA in the raster (inside the extent's bounding box but on
  # a pixel masked out when the raster was cropped/masked to the basin in
  # script 01). That discovery currently happens implicitly and repeatedly -
  # once per row, inside every fit, every LOSO fold, every validation fold.
  # Filtering ONCE here, on that exact criterion, removes the truly-useless
  # rows without over-restricting to "inside the drawn basin boundary": a
  # station just outside the polygon but still landing on a valid non-NA
  # edge pixel (grid/polygon misalignment, a slightly buffered catchment,
  # etc.) is KEPT here, exactly as it already would be by every downstream
  # fit - this is a superset of "inside the basin", not a subset. Run before
  # the CanSWE dedup/aggregation below purely so those steps don't spend
  # time on rows that are discarded anyway; it does not change which rows
  # survive dedup/aggregation, since a row dropped here would never have
  # matched anything downstream regardless of when it's removed.
  if (!is.null(swe_template)) {
    n_before_grid <- length(unique(out$station_id))
    pts <- terra::vect(out, geom = c("lon", "lat"), crs = terra::crs(swe_template))
    cellnum <- terra::cells(swe_template, pts)[, "cell"]
    template_vals <- as.numeric(terra::values(swe_template))
    station_val <- rep(NA_real_, length(cellnum))
    ok_cell <- is.finite(cellnum) & cellnum >= 1 & cellnum <= length(template_vals)
    station_val[ok_cell] <- template_vals[cellnum[ok_cell]]
    keep_grid <- is.finite(station_val)
    out <- out[keep_grid, , drop = FALSE]
    cat(sprintf("[OK] Valid-grid-cell pre-filter: kept %d/%d station(s) that resolve to a non-NA SWE pixel (not restricted to inside the basin polygon).\n",
                length(unique(out$station_id)), n_before_grid))
  }
  
  # Dedup non-CanSWE rows the old way (exact-value uniqueness is fine here -
  # there's only ever one source for ASWS/primary rows).
  non_canswe <- out[!out$is_canswe, , drop = FALSE]
  non_canswe <- unique(non_canswe)
  
  # Dedup CanSWE rows by IDENTITY (station_id + date), not by value match.
  # If daily.csv's CanSWE snapshot and CanSWE_daily.csv's snapshot disagree
  # (different swe_mm for the same station/date, e.g. after a CanSWE revision
  # or an out-of-sync re-run), this keeps exactly one row per station/date
  # instead of silently keeping both because the value tuples differ.
  canswe <- out[out$is_canswe, , drop = FALSE]
  if (nrow(canswe) > 0) {
    key <- paste(canswe$station_id, canswe$date)
    dup_keys <- key[duplicated(key)]
    if (length(unique(dup_keys)) > 0) {
      # Check whether the "duplicates" actually agree in value - if they
      # don't, this IS the stale-snapshot scenario the naming-coincidence
      # dedup used to miss, so it's worth a loud, explicit flag.
      # (Vectorized via ave() instead of a per-key vapply/rows-scan, which
      # was O(D x N) and effectively hung on large CanSWE-merged inputs.)
      n_unique_per_key <- ave(canswe$swe_mm, key, FUN = function(x) length(unique(x)))
      disagreeing_keys <- unique(key[n_unique_per_key > 1])
      if (length(disagreeing_keys) > 0) {
        log_pipeline_fallback(sprintf(
          "%d CanSWE station/date pairs had DIFFERING swe_mm values across '%s' and '%s' (likely from an out-of-sync re-run of the download script's Part A/B). Kept one value per station/date (first occurrence) rather than double-counting or guessing which is correct.",
          length(disagreeing_keys), asws_path, canswe_path))
      }
    }
    canswe <- canswe[!duplicated(key), , drop = FALSE]
  }
  
  out <- rbind(non_canswe, canswe)
  out$is_canswe <- NULL
  
  # ---- Fix #2: aggregate to one row per station per calendar month ----
  # Every downstream fit (Tier 0 QA screen, Tier 1/2/3 bias-correction fits,
  # LOSO-CV, and the k-fold/temporal-transfer validation protocol) matches a
  # station observation against a SINGLE monthly ERA5-Land grid value (one
  # column per calendar month x year - see matched_pairs()/build_row()/etc.
  # throughout this script). Any station network supplied at DAILY
  # resolution (CanSWE's native format; potentially the "primary"/ASWS file
  # too, since its underlying columns are daily readings - see
  # prep_station_file()) is therefore pseudo-replicated once loaded: every
  # day within the same calendar month gets compared against that SAME
  # monthly ERA5-Land value, which is not new information - it just inflates
  # n_ratios/leverage for whichever stations happen to report more days that
  # month, and biases sample-size-based decisions
  # (MIN_STATION_MONTHS_TIER2/3, the LOSO_CV_MIN_STATIONS/auto-tier-mode
  # gate, etc.) upward without adding real degrees of freedom.
  #
  # Collapsing here - once, for BOTH networks together, after the identity-
  # based CanSWE dedup above - to one row per station per calendar month
  # (median of that month's daily values) removes the redundancy without
  # discarding information: the median IS the station's best single
  # estimate of that month's typical SWE, the same kind of monthly
  # central-tendency statistic compute_index_values()/build_regions() use on
  # the ERA5-Land side. This step is a NO-OP for any source that is already
  # monthly (e.g. a network pre-aggregated upstream, or prep_station_file()
  # output fed in separately): grouping by station + calendar-month then
  # yields exactly one row per group already present, and the median of a
  # single value is that value.
  n_before_agg <- nrow(out)
  out$year_month <- format(out$date, "%Y-%m")
  agg_key <- paste(out$station_id, out$year_month)
  monthly_val <- ave(out$swe_mm, agg_key, FUN = function(x) stats::median(x, na.rm = TRUE))
  keep_first <- !duplicated(agg_key)
  out <- out[keep_first, , drop = FALSE]
  out$swe_mm <- monthly_val[keep_first]
  out$date <- as.Date(paste0(out$year_month, "-01"))
  out$year_month <- NULL
  cat(sprintf("[OK] Aggregated station observations to monthly resolution: %d row(s) -> %d station-month row(s) (median; no-op for sources already monthly).\n",
              n_before_agg, nrow(out)))
  
  out
}

load_optional_sar_proxy <- function(swe_template, basin, sar_path, n_time = NULL, hard_stop = FALSE) {
  if (is.null(sar_path) || !file.exists(sar_path)) {
    log_pipeline_fallback(
      sprintf("SAR_DEPTH_PATH ('%s') not found - SAR-augmented correction tiers will be skipped or fall back to non-SAR predictors.", sar_path),
      hard_stop = hard_stop
    )
    return(NULL)
  }
  
  sar <- tryCatch(rast(sar_path), error = function(e) NULL)
  if (is.null(sar)) {
    log_pipeline_fallback(sprintf("SAR proxy raster could not be read from '%s' - skipping SAR-augmented tier.", sar_path))
    return(NULL)
  }
  if (!is.null(basin)) sar <- crop(sar, project(basin, crs(sar)))
  if (!compareGeom(sar[[1]], swe_template[[1]], stopOnError = FALSE)) {
    sar <- tryCatch(resample(sar, swe_template, method = "bilinear"), error = function(e) sar)
  }
  
  if (nlyr(sar) == 1 && !is.null(n_time) && n_time > 1) {
    sar_mat <- values(sar, mat = TRUE)
    sar_mat <- sar_mat[, rep(1, n_time), drop = FALSE]
  } else {
    if (!is.null(n_time) && nlyr(sar) != n_time) {
      log_pipeline_fallback(sprintf(
        "SAR proxy layer count (%d) does not match SWE time dimension (%d) - skipping SAR-augmented tier.",
        nlyr(sar), n_time
      ))
      return(NULL)
    }
    sar_mat <- values(sar, mat = TRUE)
  }
  
  if (nrow(sar_mat) != ncell(swe_template)) {
    log_pipeline_fallback("SAR proxy raster did not align cleanly to the SWE template - skipping SAR-augmented tier.")
    return(NULL)
  }
  
  sar_max <- suppressWarnings(max(sar_mat, na.rm = TRUE))
  if (is.finite(sar_max) && sar_max < 10) {
    sar_mat <- sar_mat * 1000
  }
  sar_mat
}

# ==============================================================================
#   REGIONAL FREQUENCY ANALYSIS (RFA) HELPER FUNCTIONS  [SSPI-1 / method 2]
# ==============================================================================
# Implements the index-flood regionalization workflow described in
# Fontana et al. (2026), adapted from snow-depth extremes (block maxima, GEV)
# to monthly SWE values (full distribution, zero-inflated gamma). The unit of
# pooling here is "pixel-month within a homogeneous region" rather than
# "station-year within a homogeneous region", but the L-moment machinery
# (index value normalization, discordancy Di, heterogeneity H0, regional
# fitting, goodness-of-fit) is otherwise a direct analogue.

## ---- 1. Regionalization: elevation x aspect bands (DEM) or spatial fallback ----
build_regions <- function(swe_template, basin, dem_path,
                          n_elev_bands = 4, use_aspect = TRUE,
                          n_aspect_classes = 4, n_regions_fallback = 6,
                          min_region_pixels = 8, require_dem = REQUIRE_DEM) {
  
  n_pix <- ncell(swe_template)
  coords <- xyFromCell(swe_template, 1:n_pix)
  valid_mask <- !is.na(values(swe_template))
  
  region_id  <- rep(NA_integer_, n_pix)
  elevation  <- rep(NA_real_, n_pix)
  aspect_deg <- rep(NA_real_, n_pix)
  method_used <- "none"
  
  if (!is.null(dem_path) && file.exists(dem_path)) {
    cat(sprintf("[OK] Building regions from DEM: %s\n", dem_path))
    dem <- rast(dem_path)
    if (!is.null(basin)) dem <- crop(dem, project(basin, crs(dem)))
    if (!compareGeom(dem, swe_template, stopOnError = FALSE)) {
      dem <- resample(dem, swe_template, method = "bilinear")
    }
    elev_vals <- values(dem)[, 1]
    elevation[valid_mask] <- elev_vals[valid_mask]
    
    # Elevation quantile bands (so each band gets a comparable sample size)
    elev_breaks <- quantile(elevation[valid_mask], probs = seq(0, 1, length.out = n_elev_bands + 1),
                            na.rm = TRUE, type = 7)
    elev_breaks[1] <- -Inf; elev_breaks[length(elev_breaks)] <- Inf
    elev_band <- cut(elevation, breaks = unique(elev_breaks), labels = FALSE, include.lowest = TRUE)
    
    if (use_aspect) {
      asp <- terrain(dem, v = "aspect", unit = "degrees")
      asp_vals <- values(asp)[, 1]
      aspect_deg[valid_mask] <- asp_vals[valid_mask]
      # 4-class compass binning (N/E/S/W), flat (-1 aspect from terra) -> class 0
      aspect_class <- rep(NA_integer_, n_pix)
      aspect_class[valid_mask] <- {
        a <- aspect_deg[valid_mask]
        cls <- floor(((a + 45) %% 360) / 90) + 1   # 1=N,2=E,3=S,4=W
        cls[is.na(a) | a < 0] <- 0                  # flat / undefined
        as.integer(cls)
      }
      region_id <- as.integer(elev_band) * 10L + aspect_class
    } else {
      region_id <- as.integer(elev_band)
    }
    method_used <- "dem_elevation_aspect"
    
  } else {
    log_pipeline_fallback(
      sprintf(paste0(
        "DEM_PATH ('%s') not found - falling back to UNWEIGHTED SPATIAL K-MEANS ",
        "regionalization (no physical elevation/aspect basis). Regional pooling ",
        "(SSPI RFA and/or Tier 2/3 bias correction) will use spatially-clustered ",
        "regions instead of elevation x aspect bands. Set REQUIRE_DEM <- TRUE to ",
        "make this a hard stop instead, or fix DEM_PATH to point at the real DEM."),
        dem_path),
      hard_stop = isTRUE(require_dem)
    )
    valid_idx <- which(valid_mask)
    if (length(valid_idx) >= n_regions_fallback) {
      set.seed(42)
      km <- kmeans(coords[valid_idx, , drop = FALSE], centers = n_regions_fallback, nstart = 10)
      region_id[valid_idx] <- km$cluster
    } else {
      region_id[valid_idx] <- 1L
    }
    method_used <- "kmeans_spatial_fallback"
  }
  
  # ---- Merge undersized regions into their nearest (by centroid) neighbour ----
  # FIX: the previous version computed `others <- setdiff(all_regions,
  # small_regions)` ONCE per region, i.e. "every OTHER region not itself
  # flagged small". When every region in the domain happened to be below
  # min_region_pixels (e.g. a fine elevation x aspect grid on a small basin,
  # or a small n_regions_fallback with a lot of pixels split many ways),
  # `others` was empty for every single region and NOTHING merged - the
  # undersized regions silently survived and were pooled anyway downstream,
  # exactly the small-sample instability min_region_pixels exists to avoid.
  # Fixed by merging iteratively: at each step, merge the CURRENT smallest
  # region into its nearest still-distinct neighbour (which always exists as
  # long as more than one region remains), then recompute sizes/centroids and
  # repeat. This guarantees no region remains below min_region_pixels unless
  # exactly one region is left (nothing left to merge into).
  region_id[!valid_mask] <- NA_integer_
  n_merged_total <- 0L
  repeat {
    tab <- table(region_id[valid_mask])
    if (length(tab) <= 1) break
    sizes <- as.integer(tab)
    if (all(sizes >= min_region_pixels)) break
    
    r_small <- as.integer(names(tab)[which.min(sizes)])  # smallest region, merge first
    all_ids <- as.integer(names(tab))
    others  <- setdiff(all_ids, r_small)
    
    centroids <- sapply(all_ids, function(r) {
      idx <- which(region_id == r)
      colMeans(coords[idx, , drop = FALSE])
    })
    centroids <- t(centroids)
    rownames(centroids) <- as.character(all_ids)
    
    d <- sqrt(rowSums(sweep(centroids[as.character(others), , drop = FALSE], 2,
                            centroids[as.character(r_small), ], "-")^2))
    nearest <- others[which.min(d)]
    region_id[region_id == r_small & !is.na(region_id)] <- nearest
    n_merged_total <- n_merged_total + 1L
  }
  if (n_merged_total > 0) {
    cat(sprintf("[OK] Merged %d undersized region(s) (< %d pixels) into their nearest neighbour (iterative smallest-first).\n",
                n_merged_total, min_region_pixels))
  }
  
  region_id <- match(region_id, sort(unique(region_id[valid_mask])))  # renumber 1..K, NA stays NA
  n_regions <- length(unique(region_id[valid_mask]))
  cat(sprintf("[OK] Regionalization complete (%s): %d homogeneous sub-regions, %d valid pixels.\n",
              method_used, n_regions, sum(valid_mask)))
  
  list(region_id = region_id, elevation = elevation, aspect_deg = aspect_deg,
       n_regions = n_regions, method = method_used)
}

## ---- 2. Index value (Eq. 1-2 analogue): per-pixel, per-calendar-month scale ----
compute_index_values <- function(swe_matrix, dates_vec, ref_idx, stat = "median") {
  mon <- as.integer(format(dates_vec, "%m"))
  n_pix <- nrow(swe_matrix)
  idx_val <- matrix(NA_real_, nrow = n_pix, ncol = 12)
  fn <- if (stat == "mean") mean else median
  for (m in 1:12) {
    cols <- which(mon == m & seq_along(dates_vec) %in% ref_idx)
    if (length(cols) == 0) next
    sub <- swe_matrix[, cols, drop = FALSE]
    # Index value computed from the FULL reference-period record for that
    # calendar month (mirrors mu_i = mean of yearly maxima in Eq. 1; here the
    # "extreme" is replaced by the central tendency of the monthly record).
    idx_val[, m] <- apply(sub, 1, function(x) {
      x <- x[is.finite(x)]
      if (length(x) < 3) return(NA_real_)
      fn(x)
    })
  }
  idx_val[idx_val <= 0 | !is.finite(idx_val)] <- NA_real_
  idx_val
}

## ---- 3. L-moment ratios per pixel-month (for discordancy / heterogeneity) ----
lmom_ratios_site <- function(x) {
  x <- x[is.finite(x) & x > 0]
  if (length(x) < 5) return(c(t2 = NA_real_, t3 = NA_real_, t4 = NA_real_, n = length(x)))
  lm <- try(lmomco::lmoms(x), silent = TRUE)
  if (inherits(lm, "try-error") || is.null(lm) || any(!is.finite(lm$ratios[2:4]))) {
    return(c(t2 = NA_real_, t3 = NA_real_, t4 = NA_real_, n = length(x)))
  }
  c(t2 = lm$ratios[2], t3 = lm$ratios[3], t4 = lm$ratios[4], n = length(x))
}

## ---- 4. Discordancy measure Di (Hosking & Wallis 1993, Eq. 3-5) ----
# Standard critical values of Di (Hosking & Wallis 1997, Table 3)
.Dcr_table <- c(`5` = 1.333, `6` = 1.648, `7` = 1.917, `8` = 2.140, `9` = 2.329,
                `10` = 2.491, `11` = 2.632, `12` = 2.757, `13` = 2.869, `14` = 2.971,
                `15` = 3.000)
get_Dcr <- function(N) {
  if (N < 5) return(Inf)   # too few sites for a meaningful discordancy test - never flag
  key <- as.character(min(N, 15))
  unname(.Dcr_table[key])
}

discordancy_test <- function(lmom_mat) {
  # lmom_mat: matrix with columns t2, t3, t4 (one row per site/pixel)
  ok <- which(rowSums(is.finite(lmom_mat[, c("t2", "t3", "t4"), drop = FALSE])) == 3)
  N <- length(ok)
  Di <- rep(NA_real_, nrow(lmom_mat))
  if (N < 4) {
    # S becomes singular for N=3 and Di is non-informative for N<=4
    # (exactly as noted in the paper); skip the test, keep all sites.
    return(list(Di = Di, Dcr = get_Dcr(N), discordant = rep(FALSE, nrow(lmom_mat)), N = N))
  }
  U <- lmom_mat[ok, c("t2", "t3", "t4"), drop = FALSE]
  ubar <- colMeans(U)
  # Hosking & Wallis (1993) Eq. (5): S = (N-1)^-1 * sum_i (u_i-ubar)(u_i-ubar)^T
  # which is exactly the sample covariance matrix.
  S <- stats::cov(U)
  Sinv <- try(solve(S), silent = TRUE)
  if (inherits(Sinv, "try-error")) {
    return(list(Di = Di, Dcr = get_Dcr(N), discordant = rep(FALSE, nrow(lmom_mat)), N = N))
  }
  Dcr <- get_Dcr(N)
  for (k in seq_along(ok)) {
    d <- as.numeric(U[k, ] - ubar)
    Di[ok[k]] <- (1 / 3) * as.numeric(t(d) %*% Sinv %*% d)
  }
  list(Di = Di, Dcr = Dcr, discordant = !is.na(Di) & Di >= Dcr, N = N)
}

## ---- 5. Heterogeneity measure H0 (Hosking & Wallis 1993/1997, Eq. 6) ----
# Monte-Carlo simulation from a regional 4-parameter kappa distribution fitted
# to the (discordancy-screened) pooled normalized sample. Falls back to a
# gamma-based simulation if the kappa fit fails (common with heavily
# zero-inflated / short pooled samples).
heterogeneity_test <- function(site_samples, Nsim = 500, seed = 123) {
  # site_samples: list of numeric vectors (normalized, nonzero, finite values),
  # one per site/pixel within the candidate region for one calendar month.
  site_samples <- site_samples[sapply(site_samples, length) >= 5]
  Nsites <- length(site_samples)
  if (Nsites < 2) return(list(H0 = NA_real_, V0 = NA_real_, status = "too_few_sites"))
  
  pooled <- unlist(site_samples)
  reg_lm <- try(lmomco::lmoms(pooled), silent = TRUE)
  if (inherits(reg_lm, "try-error")) return(list(H0 = NA_real_, V0 = NA_real_, status = "lmom_failed"))
  
  par_obj <- try(lmomco::parkap(reg_lm), silent = TRUE)
  sim_fn <- NULL
  if (!inherits(par_obj, "try-error") && !is.null(par_obj) && isTRUE(par_obj$ifail == 0)) {
    sim_fn <- function(n) lmomco::quakap(stats::runif(n), par_obj)
    dist_used <- "kappa"
  } else {
    gam_par <- try(lmomco::pargam(reg_lm), silent = TRUE)
    if (inherits(gam_par, "try-error") || is.null(gam_par)) {
      return(list(H0 = NA_real_, V0 = NA_real_, status = "fit_failed"))
    }
    sim_fn <- function(n) lmomco::quagam(stats::runif(n), gam_par)
    dist_used <- "gamma"
  }
  
  site_n <- sapply(site_samples, length)
  lcv_at_site <- function(samples_list) {
    sapply(samples_list, function(x) {
      lm <- try(lmomco::lmoms(x), silent = TRUE)
      if (inherits(lm, "try-error") || is.null(lm)) return(NA_real_)
      lm$ratios[2]
    })
  }
  V0 <- function(lcv, n) {
    w <- n / sum(n)
    sqrt(sum(w * (lcv - sum(w * lcv, na.rm = TRUE))^2, na.rm = TRUE))
  }
  
  lcv_obs <- lcv_at_site(site_samples)
  V0_obs <- V0(lcv_obs, site_n)
  
  set.seed(seed)
  V0_sim <- numeric(Nsim)
  for (s in seq_len(Nsim)) {
    sim_samples <- lapply(site_n, function(n) pmax(sim_fn(n), 1e-6))
    lcv_sim <- lcv_at_site(sim_samples)
    V0_sim[s] <- V0(lcv_sim, site_n)
  }
  V0_sim <- V0_sim[is.finite(V0_sim)]
  if (length(V0_sim) < 10) return(list(H0 = NA_real_, V0 = V0_obs, status = "sim_failed"))
  
  H0 <- (V0_obs - mean(V0_sim)) / stats::sd(V0_sim)
  # First-order Monte-Carlo SE of standardized H0 from uncertainty in the
  # simulated mean. This intentionally ignores additional uncertainty in the
  # simulated SD, so it is a diagnostic rather than a formal HW confidence CI.
  H0_mc_se <- 1 / sqrt(length(V0_sim))
  z_mc <- stats::qnorm(0.5 + H0_MC_CONF_LEVEL / 2)
  H0_upper <- H0 + z_mc * H0_mc_se
  list(H0 = H0, H0_mc_se = H0_mc_se, H0_upper = H0_upper, V0 = V0_obs,
       status = paste0("ok_", dist_used))
}

## ---- 6. Pooled regional gamma fit + goodness-of-fit (KS test) ----
fit_regional_gamma <- function(pooled_sample, effective_n = NULL) {
  pooled_sample <- pooled_sample[is.finite(pooled_sample) & pooled_sample > 0]
  n_check <- if (!is.null(effective_n)) effective_n else length(pooled_sample)
  if (n_check < MIN_REGIONAL_GAMMA_SAMPLE) return(list(par = NULL, status = "insufficient_data"))
  
  lm_obj <- try(lmomco::lmoms(pooled_sample), silent = TRUE)
  if (inherits(lm_obj, "try-error") || is.null(lm_obj)) {
    return(list(par = NULL, status = "lmoms_failed"))
  }
  par_obj <- try(lmomco::pargam(lm_obj), silent = TRUE)
  if (inherits(par_obj, "try-error") || is.null(par_obj)) {
    return(list(par = NULL, status = "pargam_failed"))
  }
  list(par = par_obj, status = "ok")
}

gof_test_gamma <- function(pooled_sample, par_obj, alpha = 0.05, effective_n = NULL) {
  pooled_sample <- pooled_sample[is.finite(pooled_sample) & pooled_sample > 0]
  n_check <- if (!is.null(effective_n)) effective_n else length(pooled_sample)
  if (is.null(par_obj) || n_check < MIN_REGIONAL_GAMMA_SAMPLE) {
    return(list(D = NA_real_, p_value = NA_real_, pass = FALSE))
  }
  cdf_fun <- function(x) as.numeric(lmomco::cdfgam(x, par_obj))
  has_ties <- anyDuplicated(pooled_sample) > 0L
  ks <- try(suppressWarnings(stats::ks.test(pooled_sample, cdf_fun)), silent = TRUE)
  if (inherits(ks, "try-error"))
    return(list(D = NA_real_, p_value = NA_real_, pass = FALSE, ties = has_ties))
  list(D = unname(ks$statistic), p_value = ks$p.value, pass = ks$p.value >= alpha,
       ties = has_ties)
}

## ---- 7. Bias correction vs. independent station SWE obs (Eq. 9-10 analogue) ----
# Simple parametric scaling: x*_corr = C * x*_E5L, with C estimated as the
# ratio of station-to-ERA5-Land dimensionless quantiles, pooled across a
# probability grid (more stable than a single quantile). One C per
# region x calendar-month, applied to raw pixel SWE before normalization.
#
# NOTE: kept for reference / future use with a denser station network, but
# NOT called anywhere in this script anymore - see fit_basin_bias_correction()
# below and Section F of the RFA construction block for why it was replaced
# by a single basin-wide monthly correction.
fit_bias_correction <- function(station_obs, swe, region_id, swe_matrix, dates_vec,
                                ref_idx, index_values, clip_range = c(0.5, 2.0)) {
  mon <- as.integer(format(dates_vec, "%m"))
  n_regions <- max(region_id, na.rm = TRUE)
  C <- matrix(1, nrow = n_regions, ncol = 12)
  details <- data.frame()
  
  if (is.null(station_obs) || nrow(station_obs) == 0) {
    cat("[OK] No station observations supplied: bias correction coefficients set to 1 (no correction).\n")
    return(list(C = C, details = details))
  }
  
  station_obs$date <- as.Date(station_obs$date)
  pts <- vect(station_obs, geom = c("lon", "lat"), crs = crs(swe))
  cellnum <- cells(swe[[1]], pts)[, "cell"]
  station_obs$cell <- cellnum
  station_obs$region <- region_id[cellnum]
  station_obs <- station_obs[is.finite(station_obs$region), ]
  
  if (nrow(station_obs) == 0) {
    cat("[OK] Station observations did not intersect any valid basin pixel: skipping bias correction.\n")
    return(list(C = C, details = details))
  }
  
  probs <- seq(0.1, 0.9, by = 0.1)
  for (r in sort(unique(station_obs$region))) {
    for (m in 1:12) {
      sub <- station_obs[station_obs$region == r &
                           as.integer(format(station_obs$date, "%m")) == m &
                           station_obs$date >= dates_vec[ref_idx[1]] &
                           station_obs$date <= dates_vec[ref_idx[length(ref_idx)]], ]
      if (nrow(sub) < 10) next  # not enough station-months to estimate a quantile ratio
      
      # station dimensionless values: normalize each station by ITS OWN
      # reference-period median/mean for that month (paper's index-value logic
      # applied to the station record itself).
      station_idx_val <- stats::median(sub$swe_mm, na.rm = TRUE)
      if (!is.finite(station_idx_val) || station_idx_val <= 0) next
      station_star <- sub$swe_mm / station_idx_val
      q_station <- stats::quantile(station_star, probs, na.rm = TRUE, type = 7)
      
      # ERA5-Land dimensionless values pooled over all pixels in the region for
      # this month, reference period only.
      pix_in_region <- which(region_id == r)
      cols <- which(mon == m & seq_along(dates_vec) %in% ref_idx)
      if (length(pix_in_region) == 0 || length(cols) == 0) next
      era_vals <- swe_matrix[pix_in_region, cols, drop = FALSE]
      era_idx  <- index_values[pix_in_region, m]
      era_star <- sweep(era_vals, 1, era_idx, "/")
      era_star <- era_star[is.finite(era_star) & era_star > 0]
      if (length(era_star) < 30) next
      q_era <- stats::quantile(era_star, probs, na.rm = TRUE, type = 7)
      
      ratio <- q_station / q_era
      ratio <- ratio[is.finite(ratio) & ratio > 0]
      if (length(ratio) < 3) next
      c_val <- stats::median(ratio)
      c_val <- max(min(c_val, clip_range[2]), clip_range[1])
      C[r, m] <- c_val
      details <- rbind(details, data.frame(region = r, month = m, n_station = nrow(sub),
                                           n_era5land = length(era_star), C = c_val))
    }
  }
  cat(sprintf("[OK] Bias correction: %d region/month combinations corrected from station data.\n",
              nrow(details)))
  list(C = C, details = details)
}

## ---- 7a-QA. Tier 0: outlier screening on station observations ----
## Kanda & Fletcher (2025) removed |z|>3 outliers (~5% of records) before any
## fitting step. Screening is done on the STATION/ERA5-Land RATIO at the
## matched pixel/date, not on raw SWE, so legitimate high-SWE events (the
## peaks a drought index most needs) are never discarded just for being
## large - only station-ERA5-Land pairs whose *mismatch* is itself an outlier
## are flagged. See the bias-correction framework document, Tier 0.
qa_screen_station_obs <- function(station_obs, swe_template, swe_matrix, dates_vec,
                                  z_thresh = 3, out_dir = NULL) {
  if (is.null(station_obs) || nrow(station_obs) == 0) return(station_obs)
  station_obs$date <- as.Date(station_obs$date)
  pts <- vect(station_obs, geom = c("lon", "lat"), crs = crs(swe_template))
  cellnum <- cells(swe_template, pts)[, "cell"]
  station_obs$cell <- cellnum
  
  mon <- as.integer(format(dates_vec, "%m")); yr <- as.integer(format(dates_vec, "%Y"))
  s_mon <- as.integer(format(station_obs$date, "%m")); s_yr <- as.integer(format(station_obs$date, "%Y"))
  
  era_val <- rep(NA_real_, nrow(station_obs))
  n_qa_rows <- nrow(station_obs)
  qa_report_every <- max(1, round(n_qa_rows / 5))  # ~5 progress lines for large inputs
  if (n_qa_rows > 2000) cat(sprintf("Tier 0 (QA): matching %d station-observation row(s) to grid cells/dates...\n", n_qa_rows))
  for (i in seq_len(n_qa_rows)) {
    cix <- station_obs$cell[i]
    if (is.na(cix) || cix < 1 || cix > nrow(swe_matrix)) next
    col <- which(mon == s_mon[i] & yr == s_yr[i])
    if (length(col) == 1) era_val[i] <- swe_matrix[cix, col]
    if (n_qa_rows > 2000 && i %% qa_report_every == 0) {
      cat(sprintf("  ... %d/%d rows matched (%.0f%%)\n", i, n_qa_rows, 100 * i / n_qa_rows))
    }
  }
  
  ok <- is.finite(era_val) & is.finite(station_obs$swe_mm) & era_val > 0 & station_obs$swe_mm > 0
  ratio <- rep(NA_real_, nrow(station_obs))
  ratio[ok] <- station_obs$swe_mm[ok] / era_val[ok]
  z <- rep(NA_real_, nrow(station_obs))
  if (sum(ok) >= 3) {
    z[ok] <- (ratio[ok] - mean(ratio[ok], na.rm = TRUE)) / stats::sd(ratio[ok], na.rm = TRUE)
  }
  flagged <- ok & is.finite(z) & abs(z) > z_thresh
  n_flagged <- sum(flagged)
  
  cat(sprintf("Tier 0 (QA): flagged %d/%d station-months as |z|>%.1f outliers on the station/ERA5-Land ratio (%.1f%%).\n",
              n_flagged, sum(ok), z_thresh, 100 * n_flagged / max(1, sum(ok))))
  if (!is.null(out_dir) && n_flagged > 0) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    write.csv(station_obs[flagged, setdiff(names(station_obs), "cell"), drop = FALSE],
              file.path(out_dir, "qa_outliers_removed.csv"), row.names = FALSE)
  }
  station_obs[!flagged, setdiff(names(station_obs), "cell"), drop = FALSE]
}

## ---- 7b. Basin-wide monthly bias correction against station SWE obs ----
prep_station_file <- function(path, station_id, lon, lat, good_grades = c(0, 1)) {
  raw <- read.csv(path, stringsAsFactors = FALSE)
  raw$date <- as.Date(raw$`Start.of.Interval..UTC.`, format = "%m/%d/%Y")  # confirm actual format
  raw <- raw[raw$Grade.Code %in% good_grades & is.finite(raw$SW.Daily.Average..mm.), ]
  
  raw$year_month <- format(raw$date, "%Y-%m")
  monthly <- aggregate(SW.Daily.Average..mm. ~ year_month, data = raw, FUN = mean)
  monthly$date <- as.Date(paste0(monthly$year_month, "-01"))
  
  data.frame(station_id = station_id, lon = lon, lat = lat,
             date = monthly$date, swe_mm = monthly$SW.Daily.Average..mm.)
}
# Single multiplicative scaling factor per calendar month (<=12 numbers
# total), pooled across ALL stations and applied uniformly to every basin
# pixel for that month. For each station, uses the nearest ERA5-Land pixel;
# for each calendar month, pools the ratio (station SWE mm) / (ERA5-Land SWE
# mm at that pixel) across all stations and overlapping years, and takes the
# median as the correction factor. Deliberately NOT quantile mapping (too few
# points per probability bin) and NOT spatially interpolated (a handful of
# stations can't support kriging/IDW without inventing spatial structure the
# data can't back up) - this is the standard low-variance "delta method"
# correction appropriate for a scarcely-observed basin.
## Extended per the bias-correction framework document (Tier 1, Section 5):
##  - fits BOTH a multiplicative factor and an additive offset per calendar
##    month, and keeps whichever has the lower in-sample RMSE for that month
##    (SWE bias is not always proportional, especially near zero);
##  - replaces the old hard min_ratios_per_month on/off cutoff with
##    empirical-Bayes-style shrinkage of sparse months toward a basin-wide,
##    all-months prior, weighted by sample size (n / min_ratios_per_month,
##    capped at 1) - this degrades gracefully instead of forcing factor=1.
## v2 (this version): adds a per-pixel-scaling SLR (simple linear regression,
## station ~ a + b*ERA5) form alongside mult/add, following King et al. (2020)
## and Kanda & Fletcher (2025), and FIXES a real bug in v1: the additive
## offset was clipped only implicitly via clip_range on the ratio, so an
## unbounded d_val could and did pass straight through (observed on the real
## Nechako run as a +400 to +800 mm offset applied to every pixel in the
## basin, including low-elevation pixels far from any station). d_val (and
## the SLR intercept b0) are now bounded by the SAME clip_range, expressed
## relative to that month's median matched ERA5-Land magnitude, so the
## "sanity clip" actually protects the additive/SLR forms and not just the
## multiplicative one. See test_bias_functions.R for the validation harness
## this was developed and checked against (synthetic basin with a
## non-monotonic, elevation-dependent bias that a flat correction cannot
## fix), and the Test A/B/C checks confirming the clip engages and SLR is
## selected when the station-ERA5 relationship has both slope and intercept.
fit_basin_bias_correction <- function(station_obs, swe_template, swe_matrix,
                                      dates_vec, clip_range = c(0.5, 2.0),
                                      min_ratios_per_month = 5) {
  factor_by_month <- rep(1, 12)
  offset_by_month <- rep(0, 12)
  b0_by_month      <- rep(NA_real_, 12)
  b1_by_month      <- rep(NA_real_, 12)
  form_by_month   <- rep("mult", 12)
  details <- data.frame()
  
  empty_out <- list(factor = factor_by_month, offset = offset_by_month,
                    b0 = b0_by_month, b1 = b1_by_month,
                    form = form_by_month, details = details)
  
  if (is.null(station_obs) || nrow(station_obs) == 0) {
    cat("No station observations supplied: basin-wide bias correction set to identity (no correction).\n")
    return(empty_out)
  }
  
  station_obs$date <- as.Date(station_obs$date)
  pts <- vect(station_obs, geom = c("lon", "lat"), crs = crs(swe_template))
  
  # Nearest cell to each station (works whether or not the station falls
  # exactly inside a basin-masked NA-free pixel).
  cellnum <- cells(swe_template, pts)[, "cell"]
  station_obs$cell <- cellnum
  station_obs <- station_obs[is.finite(station_obs$cell), ]
  if (nrow(station_obs) == 0) {
    cat("Stations did not resolve to any raster cell: skipping basin-wide bias correction.\n")
    return(empty_out)
  }
  
  mon <- as.integer(format(dates_vec, "%m"))
  yr  <- as.integer(format(dates_vec, "%Y"))
  station_obs$s_mon <- as.integer(format(station_obs$date, "%m"))
  station_obs$s_yr  <- as.integer(format(station_obs$date, "%Y"))
  
  # Helper: matched (era5, station) pairs for an arbitrary subset of rows.
  matched_pairs <- function(sub) {
    era <- rep(NA_real_, nrow(sub)); sta <- rep(NA_real_, nrow(sub))
    for (i in seq_len(nrow(sub))) {
      cix <- sub$cell[i]
      if (is.na(cix) || cix < 1 || cix > nrow(swe_matrix)) next
      col <- which(mon == sub$s_mon[i] & yr == sub$s_yr[i])
      if (length(col) != 1) next
      era[i] <- swe_matrix[cix, col]; sta[i] <- sub$swe_mm[i]
    }
    ok <- is.finite(era) & is.finite(sta)
    list(era = era[ok], sta = sta[ok])
  }
  
  # Basin-wide, all-months-pooled prior for shrinkage of sparse months.
  all_pairs <- matched_pairs(station_obs)
  all_ratio <- all_pairs$sta[all_pairs$era > 0 & all_pairs$sta > 0] /
    all_pairs$era[all_pairs$era > 0 & all_pairs$sta > 0]
  all_diff  <- all_pairs$sta - all_pairs$era
  prior_ratio <- if (length(all_ratio) >= 3) stats::median(all_ratio) else 1
  prior_diff  <- if (length(all_diff)  >= 3) stats::median(all_diff)  else 0
  
  for (m in 1:12) {
    sub_m <- station_obs[station_obs$s_mon == m, ]
    pr <- matched_pairs(sub_m)
    n_m <- length(pr$era)
    
    if (n_m == 0) {
      details <- rbind(details, data.frame(month = m, n_ratios = 0, form = "mult",
                                           factor = 1, offset = 0, b0 = NA_real_, b1 = NA_real_,
                                           applied = FALSE, shrunk = FALSE))
      next
    }
    
    ratio_ok <- pr$era > 0 & pr$sta > 0
    ratios_m <- pr$sta[ratio_ok] / pr$era[ratio_ok]
    diffs_m  <- pr$sta - pr$era
    
    # Empirical-Bayes-style shrinkage toward the basin-wide prior, weighted
    # by sample size relative to min_ratios_per_month (framework doc, Section 5).
    w <- min(1, n_m / min_ratios_per_month)
    shrunk <- w < 1
    c_val_local <- if (length(ratios_m) > 0) stats::median(ratios_m) else prior_ratio
    d_val_local <- if (length(diffs_m)  > 0) stats::median(diffs_m)  else prior_diff
    c_val <- w * c_val_local + (1 - w) * prior_ratio
    d_val <- w * d_val_local + (1 - w) * prior_diff
    c_val <- max(min(c_val, clip_range[2]), clip_range[1])
    
    # FIX: bound the additive offset by the SAME clip_range, expressed
    # relative to this month's typical (median) ERA5-Land magnitude, so the
    # "sanity clip" protects the additive form too, not just the
    # multiplicative one. Without this, an unbounded additive offset (as
    # happened in the real Nechako run: +400 to +800 mm applied to every
    # pixel in the basin) passes through unchecked.
    med_era_m <- if (length(pr$era) > 0) stats::median(pr$era, na.rm = TRUE) else NA_real_
    d_bound <- c(-Inf, Inf)
    if (is.finite(med_era_m) && med_era_m > 0) {
      d_bound <- med_era_m * (clip_range - 1)
      d_val <- max(min(d_val, max(d_bound)), min(d_bound))
    }
    
    # SLR form (King et al. 2020; Kanda & Fletcher 2025) - station ~ a +
    # b*ERA5, fit by OLS. Unlike mult/add (one number applied identically to
    # every pixel), the slope+intercept scale WITH each pixel's own raw ERA5
    # value, so low-SWE and high-SWE pixels in the same calendar month are
    # not forced to receive the same absolute correction. Falls back to
    # mult/add when too few points or a degenerate (zero-variance) fit.
    slr_b0 <- NA_real_; slr_b1 <- NA_real_; rmse_slr <- Inf
    if (n_m >= 3 && stats::sd(pr$era) > 1e-8) {
      fit_slr <- tryCatch(
        stats::lm(sta ~ era, data = data.frame(sta = pr$sta, era = pr$era)),
        error = function(e) NULL)
      if (!is.null(fit_slr) && all(is.finite(stats::coef(fit_slr)))) {
        b0 <- unname(stats::coef(fit_slr)[1]); b1 <- unname(stats::coef(fit_slr)[2])
        b1 <- max(min(b1, clip_range[2]), clip_range[1])
        if (is.finite(med_era_m) && med_era_m > 0) {
          b0 <- max(min(b0, max(d_bound)), min(d_bound))
        }
        slr_b0 <- b0; slr_b1 <- b1
        rmse_slr <- sqrt(mean((slr_b0 + slr_b1 * pr$era - pr$sta)^2, na.rm = TRUE))
      }
    }
    
    # In-sample RMSE: keep whichever form (additive, multiplicative, or SLR)
    # fits this month's matched station-ERA5-Land pairs better.
    if (n_m >= MIN_BIAS_FORM_SELECTION_N) {
      rmse_mult <- sqrt(mean((pr$era * c_val - pr$sta)^2, na.rm = TRUE))
      rmse_add  <- sqrt(mean((pr$era + d_val - pr$sta)^2, na.rm = TRUE))
      rmses <- c(mult = rmse_mult, add = rmse_add, slr = rmse_slr)
      form_m <- names(which.min(rmses))
    } else {
      form_m <- "mult"
    }
    
    factor_by_month[m] <- c_val
    offset_by_month[m] <- d_val
    b0_by_month[m] <- slr_b0
    b1_by_month[m] <- slr_b1
    form_by_month[m]   <- form_m
    details <- rbind(details, data.frame(month = m, n_ratios = n_m, form = form_m,
                                         factor = c_val, offset = d_val,
                                         b0 = slr_b0, b1 = slr_b1,
                                         applied = TRUE, shrunk = shrunk))
  }
  
  n_applied <- sum(details$applied, na.rm = TRUE)
  n_shrunk  <- sum(details$shrunk, na.rm = TRUE)
  cat(sprintf("Basin-wide bias correction (Tier 1): %d/12 months corrected from %d station(s); %d partially shrunk toward the basin-wide prior (n < %d ratios/month).\n",
              n_applied, length(unique(station_obs$station_id)), n_shrunk, min_ratios_per_month))
  if (n_applied > 0) {
    applied_rows <- details[details$applied, ]
    print(applied_rows[, c("month", "form", "n_ratios", "factor", "offset", "b0", "b1", "shrunk")])
  }
  
  list(factor = factor_by_month, offset = offset_by_month, b0 = b0_by_month, b1 = b1_by_month,
       form = form_by_month, details = details)
}

## ==========================================================================
##   TIERED BIAS CORRECTION: shared predictors, Tier 2 (regionalized CDF-
##   matching), Tier 3 (MLR), and LOSO-CV tier selection
## ==========================================================================
## See the Nechako SWE bias-correction framework document for the rationale:
## Tier 1 (fit_basin_bias_correction() above) is the floor/fallback; Tier 2
## adds an elevation x aspect regional basis (reusing build_regions(), the
## same function SSPI pooling uses); Tier 3 adds a multiple linear
## regression on physical predictors (Kanda & Fletcher 2025 Table 1). A
## locally-trained Random Forest tier is deliberately NOT included - with
## only a handful of Nechako stations, training RF on them risks the same
## poor spatial-transfer behaviour both source papers document. The
## framework's recommended RF path is to train on the public CanSWE/NorSWE
## network for this ecozone offline and feed the result in like a station
## file; that step is out of scope for this script.

## ---- Shared predictor extraction for Tier 2 / Tier 3 ----
## Reuses build_regions() (defined above) so bias correction and SSPI
## regionalization draw on the same DEM-derived elevation/aspect fields, and
## adds slope + latitude/longitude. Runs once, independent of
## REGIONALIZE_SSPI (that flag controls SSPI *pooling* only; regionalizing
## the bias correction itself is a separate, always-available step).
get_bias_predictors <- function(swe_template, basin, dem_path,
                                n_elev_bands, use_aspect, n_aspect_classes,
                                n_regions_fallback, min_region_pixels,
                                n_time, sar_path = SAR_DEPTH_PATH) {
  regions_built <- build_regions(swe_template, basin, dem_path,
                                 n_elev_bands = n_elev_bands, use_aspect = use_aspect,
                                 n_aspect_classes = n_aspect_classes,
                                 n_regions_fallback = n_regions_fallback,
                                 min_region_pixels = min_region_pixels)
  n_pix  <- ncell(swe_template)
  coords <- xyFromCell(swe_template, 1:n_pix)
  slope_vals <- rep(NA_real_, n_pix)
  if (!is.null(dem_path) && file.exists(dem_path)) {
    dem <- rast(dem_path)
    if (!is.null(basin)) dem <- crop(dem, project(basin, crs(dem)))
    if (!compareGeom(dem, swe_template, stopOnError = FALSE)) {
      dem <- resample(dem, swe_template, method = "bilinear")
    }
    slp <- terrain(dem, v = "slope", unit = "degrees")
    slope_vals <- values(slp)[, 1]
  }
  sar_proxy <- load_optional_sar_proxy(swe_template, basin, sar_path, n_time = n_time)
  
  list(region_id = regions_built$region_id,
       elevation = regions_built$elevation,
       aspect_deg = regions_built$aspect_deg,
       slope_deg = slope_vals,
       lon = coords[, 1], lat = coords[, 2],
       sar_proxy = sar_proxy,
       n_regions = regions_built$n_regions,
       method = regions_built$method)
}

## ---- Degree-adaptive, monotonicity-checked polynomial CDF-matching ----
## Zahmatkesh et al. (2019) style: BS = p1*S^3 + p2*S^2 + p3*S + p4, used as
## the fitting engine for Tier 2 (fit_regionalized_qm(), immediately below).
## Degree adapts to how many matched (station, ERA5) pairs are actually
## available in a region-month, and a fitted curve that is non-monotonic (an
## overfitting artifact with few points) is rejected in favour of a lower
## degree. Validated standalone in test_poly_qm.R (linear-ish data recovers a
## monotonic high-degree fit; n=3 degrades to a 2-point line; n=2 returns
## NULL cleanly; a sparse/gapped n=8 sample correctly rejects an unstable
## non-monotonic cubic; predictions outside the training range stay
## non-negative and bounded rather than exploding).
##
## sta, era: numeric vectors of matched (station, ERA5-Land) pairs (raw mm,
## NOT dimensionless - simpler and more robust than ratio-of-quantiles when
## n is small, and it's what Zahmatkesh et al. (2019) actually fit).
## Returns NULL if not even a 2-point (degree-1) fit is possible.
fit_poly_qm <- function(sta, era, max_degree = 3, eval_grid_n = 50) {
  ok <- is.finite(sta) & is.finite(era) & era >= 0
  sta <- sta[ok]; era <- era[ok]
  n <- length(sta)
  if (n < 3) return(NULL)
  
  degree_cap <- min(max_degree, max(1, floor((n - 1) / 3)))
  grid_x <- seq(min(era), max(era), length.out = eval_grid_n)
  
  ## Leave-one-out RMSE for a given FITTED model, via the exact hat-matrix/
  ## PRESS-statistic shortcut: for ordinary least squares, the leave-one-out
  ## residual at point i equals resid_i / (1 - h_ii), where h_ii is point i's
  ## leverage (diagonal of the hat matrix) from the model fit on the FULL
  ## sample - no refitting required. This is mathematically identical to
  ## explicitly dropping each point, refitting, and predicting (verified
  ## numerically to ~1e-14), but costs one lm() fit total instead of n+1, and
  ## reuses the SAME fit object already computed below for the monotonicity/
  ## negativity screen - no extra fitting at all.
  ##
  ## This replaces an earlier version that explicitly refit
  ## `lm(sta ~ poly(era, deg), data = data.frame(sta = sta[-i], era = era[-i]))`
  ## inside a for-loop over every point, for every candidate degree, for every
  ## region-month cell in fit_regionalized_qm()'s 17 x 12 = 204-cell loop. At
  ## the real Nechako station-network scale (2958 stations feeding well-
  ## observed accumulation-season region-months with pools in the hundreds to
  ## low thousands), that cost seconds to tens of seconds PER cell with zero
  ## progress output in between - indistinguishable from a hang on a long
  ## run. The hat-matrix shortcut is 100-2500x+ faster at those sample sizes
  ## (benchmarked) and returns the identical RMSE value, so degree selection
  ## is completely unchanged - only the cost of computing it is.
  loo_rmse_from_fit <- function(fit_obj, deg) {
    if (n < deg + 2) return(Inf)
    h <- tryCatch(stats::lm.influence(fit_obj, do.coef = FALSE)$hat, error = function(e) NULL)
    if (is.null(h) || length(h) != n) return(Inf)
    denom <- 1 - h
    # Points with leverage ~1 (e.g. an isolated/extreme era value with a
    # high-degree fit) would divide by ~0 and blow up the LOO residual out
    # of proportion to any real predictive error - treat as "this degree is
    # too flexible for this sample" rather than propagate Inf/NaN silently.
    if (any(denom < 1e-6)) return(Inf)
    loo_resid <- stats::residuals(fit_obj) / denom
    if (any(!is.finite(loo_resid))) return(Inf)
    sqrt(mean(loo_resid^2))
  }
  
  best <- NULL; best_loo <- Inf
  train_df_full <- data.frame(sta = sta, era = era)
  for (deg in degree_cap:1) {
    if (n < deg + 2) next
    fit <- tryCatch(
      stats::lm(sta ~ stats::poly(era, degree = deg, raw = TRUE), data = train_df_full),
      error = function(e) NULL)
    if (is.null(fit) || !all(is.finite(stats::coef(fit)))) next
    
    pred_grid <- tryCatch(stats::predict(fit, newdata = data.frame(era = grid_x)),
                          error = function(e) rep(NA_real_, length(grid_x)))
    if (any(!is.finite(pred_grid))) next
    if (any(diff(pred_grid) < -1e-6)) next   # monotonicity screen (unchanged)
    if (any(pred_grid < -1e-6)) next          # negativity screen (unchanged)
    
    # Monotonic + non-negative is necessary but not sufficient - a
    # monotonic-yet-wild cubic on ~11 points can still overfit smoothly.
    # Only keep this degree if its LOO RMSE actually beats the best
    # lower-degree candidate found so far.
    this_loo <- loo_rmse_from_fit(fit, deg)
    if (this_loo < best_loo) {
      best <- list(degree = deg, coef = stats::coef(fit), n = n,
                   era_range = range(era), fit = fit, loo_rmse = this_loo)
      best_loo <- this_loo
    }
  }
  best
}

predict_poly_qm <- function(pqm, era_vals) {
  if (is.null(pqm)) return(era_vals)
  # A degree>=2 polynomial fit on a handful of points is only trustworthy
  # WITHIN the span of the training data; evaluating it far outside that
  # span can explode (a real risk here, since most grid pixels were never
  # within reach of any of the 3 stations). Clamp the INPUT to the fitted
  # era range before evaluating the curve (constant extrapolation at the
  # boundary), then apply the boundary-to-input ratio linearly beyond that
  # so pixels far outside the sampled range still move smoothly away from
  # 1:1, rather than being frozen at the boundary value or blowing up.
  lo <- pqm$era_range[1]; hi <- pqm$era_range[2]
  clamped <- pmin(pmax(era_vals, lo), hi)
  pred_at_clamped <- tryCatch(
    stats::predict(pqm$fit, newdata = data.frame(era = clamped)),
    error = function(e) rep(NA_real_, length(era_vals)))
  bad <- !is.finite(pred_at_clamped)
  pred_at_clamped[bad] <- clamped[bad]
  
  # boundary correction ratio (bounded, same spirit as BIAS_CORRECT_CLIP)
  ratio_at_hi <- if (hi > 0) pred_at_clamped[which.max(clamped)] / hi else 1
  ratio_at_lo <- if (lo > 0) pred_at_clamped[which.min(clamped)] / lo else 1
  ratio_at_hi <- max(min(ratio_at_hi, 2.0), 0.5)
  ratio_at_lo <- max(min(ratio_at_lo, 2.0), 0.5)
  
  pred <- pred_at_clamped
  above <- era_vals > hi
  below <- era_vals < lo
  pred[above] <- era_vals[above] * ratio_at_hi
  pred[below] <- era_vals[below] * ratio_at_lo
  pred[pred < 0] <- 0             # SWE cannot be negative
  pred
}

## ---- Tier 2: regionalized, degree-adaptive polynomial CDF-matching ----
## v2 (this version): replaces the v1 ratio-of-quantiles fit (a single
## scalar C per region-month, i.e. still a flat multiplicative correction,
## just region-specific) with a degree-adaptive, monotonicity-checked
## polynomial fit directly on matched (station, ERA5-Land) RAW mm pairs -
## BS = p1*S^3 + p2*S^2 + p3*S + p4, following Zahmatkesh et al. (2019).
## This is a genuine curve, not a constant ratio: it can express a bias that
## GROWS or SHRINKS with SWE magnitude within a single region-month (e.g. a
## small overestimate at low SWE and a larger underestimate at high SWE),
## which a single scalar ratio cannot. Fitting on raw mm rather than
## dimensionless quantile ratios is simpler and more robust when n is small
## (no probability-grid interpolation needed), and it's what Zahmatkesh et
## al. (2019) actually fit.
##
## Degree adapts to how many matched pairs are available in a region-month
## (roughly >=3 points per coefficient, capped at max_degree=3/cubic - the
## paper only uses cubic with "a large number of gridded data"), and a fit
## that is non-monotonic on its own training range (a classic overfitting
## artifact with few points) is rejected in favour of a lower degree; if even
## a degree-1 (2-point) fit isn't possible the region-month falls back to
## pooling across all snow-season months within the region, then further to
## NULL (the caller, apply_tiered_correction(), falls back to the Tier-1
## basin-wide correction when the region-month has no usable poly fit).
## See test_poly_qm.R for the standalone validation of fit_poly_qm()/
## predict_poly_qm() (degree selection, monotonicity rejection, and safe
## behaviour both interpolating and extrapolating beyond the training range).
fit_regionalized_qm <- function(station_obs, swe_template, swe_matrix, dates_vec,
                                region_id, index_values, min_n = 8,
                                clip_range = c(0.5, 2.0), max_degree = 3) {
  mon <- as.integer(format(dates_vec, "%m"))
  n_regions <- suppressWarnings(max(region_id, na.rm = TRUE))
  if (!is.finite(n_regions) || n_regions < 1) {
    return(list(C = matrix(list(NULL), 0, 12), details = data.frame()))
  }
  # C is a region x month LIST-matrix of fitted poly models (or NULL where
  # unfit), not a scalar-ratio matrix - apply_tiered_correction() calls
  # predict_poly_qm() on whichever entry is non-NULL for a pixel's region/month.
  C <- matrix(list(NULL), nrow = n_regions, ncol = 12)
  details <- data.frame()
  
  if (is.null(station_obs) || nrow(station_obs) == 0) return(list(C = C, details = details))
  
  station_obs$date <- as.Date(station_obs$date)
  pts <- vect(station_obs, geom = c("lon", "lat"), crs = crs(swe_template))
  cellnum <- cells(swe_template, pts)[, "cell"]
  station_obs$cell   <- cellnum
  station_obs$region <- region_id[cellnum]
  station_obs$mon    <- as.integer(format(station_obs$date, "%m"))
  station_obs <- station_obs[is.finite(station_obs$region) & is.finite(station_obs$cell), ]
  if (nrow(station_obs) == 0) return(list(C = C, details = details))
  
  # Matched (station, ERA5-Land) RAW mm pairs for an arbitrary pool of
  # station rows (sub) against an arbitrary pool of ERA5-Land pixels/columns
  # (pix_pool/cols_pool). No normalization needed - fit_poly_qm() fits on
  # raw mm directly (Zahmatkesh et al. 2019 style).
  matched_pairs_region <- function(sub, pix_pool, cols_pool) {
    if (nrow(sub) < min_n || length(pix_pool) == 0 || length(cols_pool) == 0) {
      return(list(sta = numeric(0), era = numeric(0)))
    }
    # Nearest-in-region-and-time pairing: for each station row, use its own
    # matched pixel/date if that pixel falls in this region's pool; this
    # keeps the pairing consistent with how Tier 1/Tier 3 match station rows
    # to specific pixels rather than pooling every pixel against every
    # station row indiscriminately.
    sta_out <- numeric(0); era_out <- numeric(0)
    for (i in seq_len(nrow(sub))) {
      cix <- sub$cell[i]
      if (!(cix %in% pix_pool)) next
      col <- cols_pool[mon[cols_pool] == sub$mon[i]]
      # sub$mon[i] already reflects the station's own calendar month; when
      # cols_pool spans multiple months (the region-pooled-months fallback),
      # this still only matches columns of the SAME month as the station row.
      if (length(col) == 0) next
      # Match by year too, mirroring fit_basin_bias_correction()'s pairing.
      yrs_col <- as.integer(format(dates_vec[col], "%Y"))
      s_yr <- as.integer(format(sub$date[i], "%Y"))
      col <- col[yrs_col == s_yr]
      if (length(col) != 1) next
      era_val <- swe_matrix[cix, col]
      if (is.finite(era_val)) { sta_out <- c(sta_out, sub$swe_mm[i]); era_out <- c(era_out, era_val) }
    }
    list(sta = sta_out, era = era_out)
  }
  
  regions_here <- sort(unique(station_obs$region))
  n_cells_total <- length(regions_here) * 12L
  cat(sprintf("Tier 2 (regionalized polynomial CDF-matching): fitting %d region-month cell(s)\n",
              n_cells_total))
  cat("          (each cell's fit_poly_qm() runs a leave-one-out CV over its matched\n")
  cat("          station-ERA5 pairs to pick a polynomial degree - this can take from\n")
  cat("          well under a second up to tens of seconds for a well-observed\n")
  cat("          region/month, so progress is reported below rather than only after\n")
  cat("          the full loop finishes.)\n")
  qm_t0 <- Sys.time()
  cell_i <- 0L
  
  for (r in regions_here) {
    pix_r <- which(region_id == r)
    sub_r <- station_obs[station_obs$region == r, ]
    for (m in 1:12) {
      cell_i <- cell_i + 1L
      sub_rm <- sub_r[sub_r$mon == m, ]
      cell_t0 <- Sys.time()
      pr <- matched_pairs_region(sub_rm, pix_r, which(mon == m))
      used <- "region_month"
      fit  <- if (length(pr$sta) >= 3) fit_poly_qm(pr$sta, pr$era, max_degree = max_degree) else NULL
      if (is.null(fit)) {
        pr2 <- matched_pairs_region(sub_r, pix_r, which(mon %in% unique(sub_r$mon)))
        fit <- if (length(pr2$sta) >= 3) fit_poly_qm(pr2$sta, pr2$era, max_degree = max_degree) else NULL
        used <- "region_pooled_months"
        if (!is.null(fit)) pr <- pr2
      }
      if (!is.null(fit)) {
        C[[r, m]] <- fit
        details <- rbind(details, data.frame(region = r, month = m, n_station = length(pr$sta),
                                             degree = fit$degree, basis = used))
      }
      cell_elapsed <- as.numeric(difftime(Sys.time(), cell_t0, units = "secs"))
      # Only print a per-cell line when it's slow enough to matter (avoids
      # flooding the console with 204 near-instant lines for sparse
      # region-months), plus a running summary every 20 cells regardless.
      if (cell_elapsed > 1) {
        cat(sprintf("  [Tier 2] region %2d, month %2d (%s): %s, n=%d, %.1fs\n",
                    r, m, month.abb[m],
                    if (!is.null(fit)) sprintf("degree-%d fit (%s)", fit$degree, used) else "no usable fit",
                    if (!is.null(fit)) length(pr$sta) else length(sub_rm$swe_mm),
                    cell_elapsed))
      }
      if (cell_i %% 20 == 0 || cell_i == n_cells_total) {
        total_elapsed <- as.numeric(difftime(Sys.time(), qm_t0, units = "secs"))
        avg_per_cell <- total_elapsed / cell_i
        eta_remaining <- avg_per_cell * (n_cells_total - cell_i)
        cat(sprintf("  [Tier 2] %d/%d region-month cell(s) done | elapsed %s | ETA %s remaining\n",
                    cell_i, n_cells_total, format_hms(total_elapsed), format_hms(eta_remaining)))
      }
    }
  }
  cat(sprintf("Tier 2 (regionalized polynomial CDF-matching): %d/%d region-months corrected from station data (rest fall back to Tier 1).\n",
              nrow(details), n_regions * 12))
  if (nrow(details) > 0) print(details)
  list(C = C, details = details)
}

## Applies the (form-aware) Tier-1 basin-wide monthly correction to a vector
## of raw ERA5-Land SWE values: additive offset, multiplicative factor, or
## SLR (station ~ b0 + b1*ERA5), whichever fit_basin_bias_correction()
## selected for that calendar month. Floored at 0 since SWE cannot be
## negative (an add/SLR form can otherwise push a small ERA5 value below 0).
apply_tier1_correction <- function(era_vals, month, bbc) {
  form <- bbc$form[month]
  if (identical(form, "add")) {
    out <- era_vals + bbc$offset[month]
  } else if (identical(form, "slr") && is.finite(bbc$b0[month]) && is.finite(bbc$b1[month])) {
    out <- bbc$b0[month] + bbc$b1[month] * era_vals
  } else {
    out <- era_vals * bbc$factor[month]
  }
  pmax(out, 0)
}

## Per-pixel CORRECTED VALUE for a given month: region-specific degree-
## adaptive polynomial CDF-match (Zahmatkesh-style, via fit_poly_qm()/
## predict_poly_qm()) where a region/month has a usable fit, else the
## form-aware Tier-1 basin-wide correction as fallback. C is now a region x
## month LIST-matrix of fitted poly models (see fit_regionalized_qm()), not
## a scalar-ratio matrix, since the correction curve is generally nonlinear.
apply_tiered_correction <- function(era_vals, region_id, C, bbc, month) {
  out <- apply_tier1_correction(era_vals, month, bbc)
  if (is.null(C) || nrow(C) == 0) return(out)
  has_region <- !is.na(region_id) & region_id >= 1 & region_id <= nrow(C)
  idx <- which(has_region)
  if (length(idx) == 0) return(out)
  # Group by region so predict_poly_qm() is called once per distinct
  # (region present among idx) rather than once per pixel.
  regs_here <- unique(region_id[idx])
  for (r in regs_here) {
    fit <- C[[r, month]]
    if (is.null(fit)) next
    sub_idx <- idx[region_id[idx] == r]
    out[sub_idx] <- predict_poly_qm(fit, era_vals[sub_idx])
  }
  out
}

## ---- Tier 3: multiple linear regression bias correction ----
## station_SWE ~ ERA5_SWE + elevation + slope + sin(aspect) + cos(aspect) +
## longitude + latitude + sin(2*pi*doy/365) + cos(2*pi*doy/365), following
## the Kanda & Fletcher (2025) predictor set. Their station-specific
## "elevation difference" term is dropped here since it has no well-defined
## value away from a station (see the bias-correction framework document);
## elevation/slope/aspect/lat/lon instead come from the same DEM-derived
## fields used for SSPI regionalization, so the gridded prediction stays
## well-posed everywhere ERA5-Land itself has a value.
fit_mlr_correction <- function(station_obs, swe_template, swe_matrix, dates_vec,
                               elevation, aspect_deg, slope_deg, lon, lat,
                               min_n = 15, sar_proxy = NULL) {
  station_obs$date <- as.Date(station_obs$date)
  pts <- vect(station_obs, geom = c("lon", "lat"), crs = crs(swe_template))
  cellnum <- cells(swe_template, pts)[, "cell"]
  station_obs$cell <- cellnum
  station_obs <- station_obs[is.finite(station_obs$cell), ]
  if (nrow(station_obs) < min_n) {
    cat(sprintf("Tier 3 (MLR): only %d station-obs (< %d) - skipping.\n", nrow(station_obs), min_n))
    return(NULL)
  }
  
  mon_all <- as.integer(format(dates_vec, "%m"))
  yr_all  <- as.integer(format(dates_vec, "%Y"))
  station_obs$s_mon <- as.integer(format(station_obs$date, "%m"))
  station_obs$s_yr  <- as.integer(format(station_obs$date, "%Y"))
  doy <- as.integer(format(station_obs$date, "%j"))
  station_obs$sin_doy <- sin(2 * pi * doy / 365)
  station_obs$cos_doy <- cos(2 * pi * doy / 365)
  asp_rad <- aspect_deg * pi / 180
  
  build_row <- function(i) {
    cix <- station_obs$cell[i]
    col <- which(mon_all == station_obs$s_mon[i] & yr_all == station_obs$s_yr[i])
    if (length(col) != 1 || is.na(cix)) return(NULL)
    era_val <- swe_matrix[cix, col]
    if (!is.finite(era_val)) return(NULL)
    row <- data.frame(station_swe = station_obs$swe_mm[i], era5_swe = era_val,
                      elevation = elevation[cix], slope = slope_deg[cix],
                      aspect_sin = sin(asp_rad[cix]), aspect_cos = cos(asp_rad[cix]),
                      lon = lon[cix], lat = lat[cix],
                      sin_doy = station_obs$sin_doy[i], cos_doy = station_obs$cos_doy[i])
    if (!is.null(sar_proxy)) {
      sar_val <- sar_proxy[cix, col]
      if (is.finite(sar_val)) row$sar_proxy <- sar_val
    }
    row
  }
  train_rows <- lapply(seq_len(nrow(station_obs)), build_row)
  train_df <- do.call(rbind, train_rows[!sapply(train_rows, is.null)])
  if (is.null(train_df) || nrow(train_df) == 0) {
    cat("Tier 3 (MLR): no rows survived pixel/date matching - skipping.\n")
    return(NULL)
  }
  train_df <- train_df[stats::complete.cases(train_df) &
                         is.finite(train_df$station_swe) & is.finite(train_df$era5_swe), ]
  if (nrow(train_df) < min_n) {
    cat(sprintf("Tier 3 (MLR): only %d complete rows after matching (< %d) - skipping.\n",
                nrow(train_df), min_n))
    return(NULL)
  }
  
  formula_terms <- "station_swe ~ era5_swe + elevation + slope + aspect_sin + aspect_cos + lon + lat + sin_doy + cos_doy"
  if (!is.null(sar_proxy) && "sar_proxy" %in% names(train_df)) {
    formula_terms <- paste0(formula_terms, " + sar_proxy")
  }
  fit <- tryCatch(
    stats::lm(stats::as.formula(formula_terms), data = train_df),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    cat("Tier 3 (MLR): model fit failed - skipping.\n")
    return(NULL)
  }
  cat(sprintf("Tier 3 (MLR): fit on %d station-months, R^2 = %.2f\n",
              nrow(train_df), summary(fit)$r.squared))
  fit
}

## Applies a fitted Tier-3 MLR model to a full grid column (one date at a
## time). Falls back to the raw ERA5-Land value wherever prediction fails
## or would be negative (SWE cannot be negative).
predict_mlr_grid <- function(fit, era5_vec, elevation, slope_deg, aspect_deg, lon, lat, date, sar_proxy = NULL) {
  if (is.null(fit)) return(era5_vec)
  asp_rad <- aspect_deg * pi / 180
  doy <- as.integer(format(date, "%j"))
  newdata <- data.frame(era5_swe = era5_vec, elevation = elevation, slope = slope_deg,
                        aspect_sin = sin(asp_rad), aspect_cos = cos(asp_rad),
                        lon = lon, lat = lat,
                        sin_doy = sin(2 * pi * doy / 365), cos_doy = cos(2 * pi * doy / 365))
  if (!is.null(sar_proxy)) newdata$sar_proxy <- sar_proxy
  pred <- tryCatch(stats::predict(fit, newdata = newdata), error = function(e) rep(NA_real_, length(era5_vec)))
  bad <- !is.finite(pred)
  pred[bad] <- era5_vec[bad]
  pred[is.finite(pred) & pred < 0] <- 0
  pred
}
## ---- Optional Tier 2 SAR-augmented candidate (selector hook) ----
## When SAR proxy data are available, this branch reuses the MLR engine with
## SAR included as an additional predictor. When SAR is absent, callers should
## treat the candidate as unavailable and fall back to the non-SAR tiers.
fit_tier2_sar_correction <- function(station_obs, swe_template, swe_matrix, dates_vec,
                                     elevation, aspect_deg, slope_deg, lon, lat,
                                     min_n = MIN_STATION_MONTHS_TIER3,
                                     sar_proxy = NULL) {
  fit_mlr_correction(station_obs, swe_template, swe_matrix, dates_vec,
                     elevation, aspect_deg, slope_deg, lon, lat,
                     min_n = min_n, sar_proxy = sar_proxy)
}



## ---- Reduce the LOSO-CV held-out station set: dedupe by pixel, then ----
## ---- stratify by homogeneous region, capped at LOSO_CV_MAX_STATIONS   ----
## Addresses the 2958-distinct-station / 790-pixel mismatch described in the
## LOSO_CV_MAX_STATIONS comment above. Two-stage reduction:
##
##  Stage 1 (pixel dedup): station_obs can carry many station_id values that
##  resolve to the SAME ERA5-Land grid cell (dense CanSWE points, nearby ASWS
##  gauges, etc). For LOSO-CV purposes, holding out a second station in a cell
##  already represented tests no new spatial location, so only ONE station
##  per occupied cell is kept as a holdout CANDIDATE (the one with the most
##  observation rows, for a slightly more informative held-out test).
##
##  Stage 2 (region-stratified sampling): the surviving one-per-cell
##  candidates are then sampled WITHOUT replacement, allocated proportionally
##  to how many candidate stations fall in each build_regions() region, with
##  a floor of `min_per_region` per region (so small regions are not silently
##  dropped from validation) and a total cap of `max_stations`.
##
## Returns the character vector of station_id values to actually hold out one
## at a time in loso_cv_bias_tiers(); Tier 1/2/3 are still FIT on the full,
## un-subsampled station_obs in every fold - only which stations get pulled
## out and scored is reduced.
build_loso_holdout_stations <- function(station_obs, swe_template, region_id,
                                        max_stations = LOSO_CV_MAX_STATIONS,
                                        min_per_region = LOSO_CV_MIN_PER_REGION,
                                        seed = BASIN_SWEI_RANDOM_SEED) {
  all_stations <- unique(station_obs$station_id)
  n_all <- length(all_stations)
  if (!is.finite(max_stations) || n_all <= max_stations) {
    cat(sprintf("[LOSO-CV setup] %d distinct station(s) is at/under LOSO_CV_MAX_STATIONS (%s) - using all of them, no subsampling needed.\n",
                n_all, if (is.finite(max_stations)) max_stations else "Inf"))
    return(all_stations)
  }
  
  cat(sprintf("\n[LOSO-CV setup] %d distinct stations exceeds LOSO_CV_MAX_STATIONS (%d).\n", n_all, max_stations))
  cat("[LOSO-CV setup] Reducing the HELD-OUT set (Tier 1/2/3 still fit on the FULL network every fold):\n")
  cat("[LOSO-CV setup]   Stage 1: one representative station per occupied grid cell (co-located\n")
  cat("[LOSO-CV setup]            stations are redundant for spatial-transfer testing).\n")
  cat("[LOSO-CV setup]   Stage 2: stratified sampling across build_regions() regions, so every\n")
  cat("[LOSO-CV setup]            region keeps validation coverage rather than being crowded out.\n")
  
  ## ---- Stage 1: map each station to a cell, keep the best-observed one per cell ----
  station_obs$date <- as.Date(station_obs$date)
  pts <- vect(station_obs, geom = c("lon", "lat"), crs = crs(swe_template))
  cellnum <- cells(swe_template, pts)[, "cell"]
  station_obs$cell <- cellnum
  station_obs <- station_obs[is.finite(station_obs$cell), ]
  
  n_rows_tab <- table(station_obs$station_id)
  cell_of_station <- tapply(station_obs$cell, station_obs$station_id, function(x) x[1])
  
  cand_df <- data.frame(
    station_id = names(cell_of_station),
    cell       = as.integer(cell_of_station),
    n_rows     = as.integer(n_rows_tab[names(cell_of_station)]),
    stringsAsFactors = FALSE
  )
  cand_df <- cand_df[is.finite(cand_df$cell), ]
  cand_df <- cand_df[order(cand_df$cell, -cand_df$n_rows), ]
  cand_df <- cand_df[!duplicated(cand_df$cell), ]  # best (most-observed) station per cell
  
  cat(sprintf("[LOSO-CV setup]   Stage 1 result: %d occupied cell(s) -> %d candidate station(s) (from %d total).\n",
              length(unique(cand_df$cell)), nrow(cand_df), n_all))
  
  ## ---- Stage 2: stratify candidates by region, proportional allocation + floor ----
  cand_df$region <- region_id[cand_df$cell]
  cand_df <- cand_df[is.finite(cand_df$region), ]
  
  if (nrow(cand_df) <= max_stations) {
    cat(sprintf("[LOSO-CV setup]   Stage 2: %d one-per-cell candidate(s) already <= cap (%d) - keeping all of them, no further sampling needed.\n",
                nrow(cand_df), max_stations))
    return(cand_df$station_id)
  }
  
  set.seed(seed)
  regions_present <- sort(unique(cand_df$region))
  n_per_region_avail <- table(cand_df$region)
  
  # Proportional allocation of the total cap across regions, by candidate
  # count, then enforce the per-region floor (capped at that region's own
  # availability) and re-normalize if the floors overran the total budget.
  prop_alloc <- floor(max_stations * n_per_region_avail / sum(n_per_region_avail))
  alloc <- pmax(prop_alloc, pmin(min_per_region, n_per_region_avail))
  names(alloc) <- names(n_per_region_avail)
  
  # If floors pushed total allocation over the cap, shave down the largest
  # (proportionally-allocated) regions first rather than breaking the floor.
  while (sum(alloc) > max_stations) {
    shaveable <- which(alloc > pmin(min_per_region, n_per_region_avail))
    if (length(shaveable) == 0) break  # every region already at its floor; accept slight overrun
    r_shave <- shaveable[which.max(alloc[shaveable])]
    alloc[r_shave] <- alloc[r_shave] - 1
  }
  # If proportional+floor allocation left room under the cap, hand the
  # leftover to regions that still have unused candidates, largest region first.
  leftover <- max_stations - sum(alloc)
  if (leftover > 0) {
    growable_order <- order(-n_per_region_avail)
    gi <- 1
    while (leftover > 0 && any(alloc[growable_order] < n_per_region_avail[growable_order])) {
      r <- growable_order[gi]
      if (alloc[r] < n_per_region_avail[r]) {
        alloc[r] <- alloc[r] + 1
        leftover <- leftover - 1
      }
      gi <- gi %% length(growable_order) + 1
    }
  }
  
  holdout <- character(0)
  alloc_report <- data.frame(region = integer(0), n_candidates = integer(0), n_selected = integer(0))
  for (r in regions_present) {
    idx_r <- which(cand_df$region == r)
    n_take <- min(alloc[as.character(r)], length(idx_r))
    sel <- if (length(idx_r) > 0) sample(cand_df$station_id[idx_r], n_take) else character(0)
    holdout <- c(holdout, sel)
    alloc_report <- rbind(alloc_report, data.frame(region = r, n_candidates = length(idx_r), n_selected = n_take))
  }
  
  cat(sprintf("[LOSO-CV setup]   Stage 2 result: %d station(s) selected across %d region(s) (target cap: %d).\n",
              length(holdout), length(regions_present), max_stations))
  print(alloc_report)
  cat(sprintf("[LOSO-CV setup] -> LOSO-CV will hold out %d station(s) instead of %d (reduction: %.1fx). Tier 1/2/3 fits\n",
              length(holdout), n_all, n_all / length(holdout)))
  cat("[LOSO-CV setup]    themselves still use the FULL station network in every fold - only the\n")
  cat("[LOSO-CV setup]    held-out/scored set is reduced.\n")
  
  holdout
}

## ---- Leave-one-station-out (LOSO) CV: pick the best-performing tier ----
## Refits each tier with one station withheld at a time and predicts at the
## withheld station's matched pixel/date, pooling squared errors across all
## held-out predictions - the "spatial-block" style validation Kanda &
## Fletcher (2025) show matters most (their random-CV R^2 of 0.94 fell to
## 0.68 under spatial transfer). Requires at least LOSO_CV_MIN_STATIONS
## distinct stations; the caller checks that before invoking this.
##
## `holdout_stations`: optional character vector restricting which distinct
## station_id values are actually pulled out and scored one at a time (see
## build_loso_holdout_stations() above). Defaults to NULL, which reproduces
## the original exhaustive behaviour (every distinct station in station_obs
## is held out in turn) - training data for Tier 1/2/3 is ALWAYS the full
## station_obs minus whichever single station is currently held out,
## regardless of this argument; only the loop's iteration list is affected.
## ---- LOSO-CV per-station worker ----
## Defined at top level (not nested inside loso_cv_bias_tiers()) so it
## serializes cleanly against each cluster node's .GlobalEnv - the same
## mechanism the .het_worker() pattern above relies on. This function's body
## is IDENTICAL in substance to the inside of the old serial "for (s_i in
## seq_along(stations))" loop it replaces: same train/test split, same four
## fits, same scoring logic, same accumulation into a per-station results
## data.frame. Nothing about what gets computed has changed - only that a
## worker process computes one station's fold instead of the master doing
## every station one after another. Free variables referenced here
## (station_obs, swe_template, swe_matrix, dates_vec, region_id, elevation,
## aspect_deg, slope_deg, lon, lat, clip_range, min_n_tier2, min_n_tier3,
## sar_proxy, mon_all, yr_all, index_values_cv) are clusterExport()-ed into
## each worker's .GlobalEnv by loso_cv_bias_tiers() before any task is sent.
.loso_worker <- function(s) {
  station_t0 <- Sys.time()
  empty_results <- data.frame(tier = character(), sq_err = numeric(), pred = numeric(), obs = numeric())
  train_obs <- station_obs[station_obs$station_id != s, ]
  test_obs  <- station_obs[station_obs$station_id == s, ]
  if (nrow(train_obs) == 0 || nrow(test_obs) == 0) {
    return(list(station = s, results = empty_results,
                elapsed = as.numeric(difftime(Sys.time(), station_t0, units = "secs")),
                skipped = TRUE))
  }
  
  bbc <- fit_basin_bias_correction(train_obs, swe_template, swe_matrix, dates_vec,
                                   clip_range = clip_range)
  qm  <- fit_regionalized_qm(train_obs, swe_template, swe_matrix, dates_vec,
                             region_id, index_values_cv, min_n = min_n_tier2,
                             clip_range = clip_range)
  mlr <- fit_mlr_correction(train_obs, swe_template, swe_matrix, dates_vec,
                            elevation, aspect_deg, slope_deg, lon, lat, min_n = min_n_tier3)
  mlr_sar <- if (!is.null(sar_proxy)) fit_tier2_sar_correction(train_obs, swe_template, swe_matrix, dates_vec,
                                                               elevation, aspect_deg, slope_deg, lon, lat, min_n = min_n_tier3,
                                                               sar_proxy = sar_proxy) else NULL
  
  test_obs$date <- as.Date(test_obs$date)
  pts <- vect(test_obs, geom = c("lon", "lat"), crs = crs(swe_template))
  test_obs$cell <- cells(swe_template, pts)[, "cell"]
  
  station_results <- empty_results
  for (i in seq_len(nrow(test_obs))) {
    cix <- test_obs$cell[i]
    if (is.na(cix)) next
    d <- test_obs$date[i]
    col <- which(mon_all == as.integer(format(d, "%m")) & yr_all == as.integer(format(d, "%Y")))
    if (length(col) != 1) next
    era_val <- swe_matrix[cix, col]
    obs_val <- test_obs$swe_mm[i]
    if (!is.finite(era_val) || !is.finite(obs_val)) next
    m <- as.integer(format(d, "%m"))
    
    pred1 <- apply_tier1_correction(era_val, m, bbc)
    pred2 <- apply_tiered_correction(era_val, region_id[cix], qm$C, bbc, m)
    station_results <- rbind(station_results,
                             data.frame(tier = "tier1", sq_err = (pred1 - obs_val)^2, pred = pred1, obs = obs_val),
                             data.frame(tier = "tier2", sq_err = (pred2 - obs_val)^2, pred = pred2, obs = obs_val))
    if (!is.null(mlr)) {
      pred3 <- predict_mlr_grid(mlr, era_val, elevation[cix], slope_deg[cix],
                                aspect_deg[cix], lon[cix], lat[cix], d)
      if (is.finite(pred3)) {
        station_results <- rbind(station_results, data.frame(tier = "tier3", sq_err = (pred3 - obs_val)^2, pred = pred3, obs = obs_val))
      }
    }
    if (!is.null(mlr_sar)) {
      sar_col <- if (!is.null(sar_proxy)) sar_proxy[cix, which(mon_all == as.integer(format(d, "%m")) & yr_all == as.integer(format(d, "%Y")))] else NULL
      pred4 <- predict_mlr_grid(mlr_sar, era_val, elevation[cix], slope_deg[cix],
                                aspect_deg[cix], lon[cix], lat[cix], d,
                                sar_proxy = sar_col)
      if (is.finite(pred4)) {
        station_results <- rbind(station_results, data.frame(tier = "tier2_sar", sq_err = (pred4 - obs_val)^2, pred = pred4, obs = obs_val))
      }
    }
  }
  
  list(station = s, results = station_results,
       elapsed = as.numeric(difftime(Sys.time(), station_t0, units = "secs")),
       skipped = FALSE)
}

loso_cv_bias_tiers <- function(station_obs, swe_template, swe_matrix, dates_vec,
                               region_id, elevation, aspect_deg, slope_deg, lon, lat,
                               clip_range, min_n_tier2, min_n_tier3, sar_proxy = NULL,
                               holdout_stations = NULL) {
  stations <- if (!is.null(holdout_stations)) holdout_stations else unique(station_obs$station_id)
  n_stations_loso <- length(stations)
  index_values_cv <- compute_index_values(swe_matrix, dates_vec, seq_along(dates_vec),
                                          stat = INDEX_VALUE_STAT)
  mon_all <- as.integer(format(dates_vec, "%m"))
  yr_all  <- as.integer(format(dates_vec, "%Y"))
  
  cat(sprintf("\n[LOSO-CV] Starting leave-one-station-out cross-validation across %d distinct station(s).\n",
              n_stations_loso))
  cat("[LOSO-CV] Each station refits Tier 1 (basin-wide), Tier 2 (regionalized QM),\n")
  cat("          Tier 3 (MLR), and Tier 2-SAR (if available) with that station held out -\n")
  cat("          this is the most time-consuming step in bias correction. Every station's\n")
  cat("          fold only depends on its own train/test split, so folds are farmed out\n")
  cat("          across parallel workers (the exact same per-station computation as\n")
  cat("          before, just running concurrently instead of one after another).\n")
  
  if (n_stations_loso == 0) return(list(best = "tier1", table = data.frame()))
  
  ## ---- Decide serial vs. parallel execution ----
  ## PSOCK cluster startup (makeCluster()) spawns n_cores_loso brand-new R
  ## processes, each of which must open a TCP socket back to this session. On
  ## a machine where that handshake is blocked or slowed (Windows Defender/
  ## antivirus process-spawn scanning, a corporate firewall, VPN-intercepted
  ## loopback traffic - all invisible from inside R, since PSOCK workers
  ## default to outfile = NULL, which silently discards anything a stuck
  ## worker would have printed), makeCluster() can sit for HOURS with zero
  ## CPU use and no error - indistinguishable, from the console alone, from
  ## real computation. That fixed startup risk/cost is the same whether
  ## there are 15 folds or 500; it only pays off once the per-fold work (a
  ## full Tier 1/2/3 refit) is large enough to dwarf it. Below
  ## LOSO_CV_SEQUENTIAL_MAX distinct held-out stations, .loso_worker() runs
  ## directly in a plain loop instead of via a cluster - IDENTICAL numerical
  ## results to the parallel path (same function, same inputs; nothing in
  ## fit_basin_bias_correction()/fit_regionalized_qm()/fit_poly_qm()/
  ## fit_mlr_correction() sets a seed or otherwise depends on execution
  ## order), just without ever calling makeCluster(). This is a pure
  ## execution-strategy switch, not a methodology change.
  if (!exists("LOSO_CV_SEQUENTIAL_MAX", inherits = TRUE)) LOSO_CV_SEQUENTIAL_MAX <- 30
  
  ## .loso_worker() is defined at top level, so (like every cluster worker)
  ## it looks up its free variables in .GlobalEnv, not in loso_cv_bias_tiers()'s
  ## local frame. The parallel path satisfies that via clusterExport() into
  ## each WORKER's own .GlobalEnv; running serially in THIS process means
  ## temporarily mirroring the identical values into the actual .GlobalEnv
  ## instead - saving/restoring whatever those names were bound to
  ## beforehand (on.exit) so this can never leak state into any other part
  ## of the script that happens to reuse a generic name like `region_id` or
  ## `elevation` for something else.
  .loso_worker_vars <- c("station_obs", "swe_template", "swe_matrix", "dates_vec",
                         "region_id", "elevation", "aspect_deg", "slope_deg",
                         "lon", "lat", "clip_range", "min_n_tier2", "min_n_tier3",
                         "sar_proxy", "mon_all", "yr_all", "index_values_cv")
  
  run_loso_serial <- function() {
    prior_existed <- setNames(logical(length(.loso_worker_vars)), .loso_worker_vars)
    prior_vals    <- setNames(vector("list", length(.loso_worker_vars)), .loso_worker_vars)
    for (nm in .loso_worker_vars) {
      prior_existed[nm] <- exists(nm, envir = .GlobalEnv, inherits = FALSE)
      if (prior_existed[nm]) prior_vals[[nm]] <- get(nm, envir = .GlobalEnv, inherits = FALSE)
    }
    on.exit({
      for (nm in .loso_worker_vars) {
        if (prior_existed[nm]) assign(nm, prior_vals[[nm]], envir = .GlobalEnv)
        else if (exists(nm, envir = .GlobalEnv, inherits = FALSE)) rm(list = nm, envir = .GlobalEnv)
      }
    }, add = TRUE)
    live_vals <- mget(.loso_worker_vars, inherits = TRUE)
    for (nm in .loso_worker_vars) assign(nm, live_vals[[nm]], envir = .GlobalEnv)
    
    n_tasks <- n_stations_loso
    out <- vector("list", n_tasks)
    t0 <- Sys.time()
    for (i in seq_len(n_tasks)) {
      out[[i]] <- .loso_worker(stations[i])
      s_val <- out[[i]]
      elapsed_total <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      rate    <- i / elapsed_total
      eta_sec <- if (is.finite(rate) && rate > 0) (n_tasks - i) / rate else 0
      if (isTRUE(s_val$skipped)) {
        cat(sprintf("[LOSO-CV] %4d/%4d done | station '%s' skipped (no train or test rows) | elapsed %s | ETA %s remaining\n",
                    i, n_tasks, s_val$station, format_hms(elapsed_total), format_hms(eta_sec)))
      } else {
        cat(sprintf("[LOSO-CV] %4d/%4d done | station '%s' took %.1fs | elapsed %s | ETA %s remaining\n",
                    i, n_tasks, s_val$station, s_val$elapsed, format_hms(elapsed_total), format_hms(eta_sec)))
      }
    }
    out
  }
  
  if (!isTRUE(ENABLE_PSOCK_PARALLEL) || n_stations_loso <= LOSO_CV_SEQUENTIAL_MAX) {
    cat(sprintf("[LOSO-CV] %d station(s) is at/under LOSO_CV_SEQUENTIAL_MAX (%d): running folds\n",
                n_stations_loso, LOSO_CV_SEQUENTIAL_MAX))
    cat("          SEQUENTIALLY rather than starting a parallel cluster - for this few folds,\n")
    cat("          PSOCK cluster startup overhead (and its risk of hanging silently on some\n")
    cat("          networks/firewalls) outweighs any speedup. Identical .loso_worker()\n")
    cat("          computation as the parallel path, just one fold at a time.\n")
    station_out <- run_loso_serial()
    
  } else {
    # ---- DIAGNOSTIC: object sizes being handed to the cluster ----
    # SpatRaster objects (terra) hold an external C++ pointer that does NOT
    # serialize meaningfully across a PSOCK socket connection - handing one to
    # clusterExport() is the classic cause of a hang that shows near-zero CPU
    # in Task Manager (the workers are stuck on a bad/blocking serialization,
    # not computing). swe_template here is `swe[[1]]`, a live SpatRaster, so it
    # must NEVER be included in a clusterExport() varlist - see "raster recipe"
    # rebuild below, which exports only plain, fully-serializable geometry
    # (crs string + extent + nrow/ncol) instead of the raster object itself.
    # station_obs / swe_matrix are large-but-ordinary R objects (data.frame /
    # matrix) and DO serialize correctly over sockets - they're just slow to
    # copy if large, which is why their sizes are reported below.
    cat(sprintf("[LOSO-CV setup] Object sizes for cluster transfer: station_obs = %.1f MB | swe_matrix = %.1f MB\n",
                as.numeric(object.size(station_obs)) / 1e6,
                as.numeric(object.size(swe_matrix)) / 1e6))
    
    # ---- Build a lightweight, fully-serializable "raster recipe" ----
    # Captures only what .loso_worker() / fit_basin_bias_correction() /
    # fit_regionalized_qm() / fit_mlr_correction() actually need from
    # swe_template: crs(swe_template) and cells(swe_template, pts) for
    # lon/lat -> cell-index lookup. Neither of those touches pixel VALUES, only
    # geometry (extent, resolution, dimensions) and the CRS string - so a
    # brand-new empty SpatRaster built from this recipe on each worker is
    # functionally identical to the real swe_template for every call made
    # inside this cluster, without ever serializing the real object itself.
    swe_template_recipe <- list(
      crs_wkt = terra::crs(swe_template),
      extent  = as.vector(terra::ext(swe_template)),
      nrow    = terra::nrow(swe_template),
      ncol    = terra::ncol(swe_template)
    )
    cat(sprintf("[LOSO-CV setup] swe_template raster recipe captured (crs/extent/dims only, no pixel data) - this is what gets exported, not the SpatRaster itself.\n"))
    
    n_cores_loso <- min(MAX_PSOCK_WORKERS, max(1L, detectCores() - 1L))
    cat(sprintf("[LOSO-CV] Running %d station fold(s) in parallel (%d cores)\n", n_stations_loso, n_cores_loso))
    cluster_t0 <- Sys.time()
    cluster_log <- file.path(tempdir(), "loso_cv_cluster_worker_log.txt")
    cat(sprintf("[%s] [LOSO-CV setup] makeCluster(%d) starting (60s per-worker connect timeout;\n            worker stdout/stderr -> %s;\n            falls back to sequential execution automatically if this fails)...\n",
                format(Sys.time(), "%H:%M:%S"), n_cores_loso, cluster_log))
    flush.console()
    
    # ---- Robust cluster startup ----
    # timeout = 60 bounds the per-worker socket handshake (the default is 30
    # DAYS - effectively "hang forever", which is almost certainly what a
    # multi-hour stall right at this line means: some worker's TCP connect-
    # back to this session is being blocked - firewall, antivirus process-
    # spawn scanning, or a corporate proxy intercepting loopback traffic -
    # and the silence is itself expected: PSOCK workers default to
    # outfile = NULL, which DISCARDS anything a stuck worker would have
    # printed, so a real connection failure looks identical to "just slow"
    # from the console. outfile = cluster_log routes worker stdout/stderr to
    # a real file instead, so a future large-N run that hits this same
    # problem leaves a diagnosable trace rather than silence. If
    # makeCluster() still fails or times out, fall back to the identical
    # serial loop used above rather than losing the whole bias-correction
    # run to a cluster that will never come up.
    cl_loso <- tryCatch(
      parallel::makeCluster(n_cores_loso, timeout = 60, outfile = cluster_log),
      error = function(e) {
        cat(sprintf("[LOSO-CV setup] makeCluster() failed or timed out (%s).\n            See worker log: %s\n",
                    conditionMessage(e), cluster_log))
        NULL
      }
    )
    
    if (is.null(cl_loso)) {
      cat("[LOSO-CV] Falling back to SEQUENTIAL execution for this run (identical per-fold\n")
      cat("          computation, just without a cluster). Check the worker log above if\n")
      cat("          this keeps happening on larger station networks.\n")
      station_out <- run_loso_serial()
      
    } else {
      cat(sprintf("[%s] [LOSO-CV setup] makeCluster() done (%.1fs). Exporting objects to workers...\n",
                  format(Sys.time(), "%H:%M:%S"),
                  as.numeric(difftime(Sys.time(), cluster_t0, units = "secs"))))
      flush.console()
      
      export_t0 <- Sys.time()
      clusterExport(cl_loso, varlist = c(".loso_worker", "fit_basin_bias_correction",
                                         "fit_regionalized_qm", "fit_poly_qm", "predict_poly_qm",
                                         "fit_mlr_correction", "fit_tier2_sar_correction",
                                         "predict_mlr_grid", "apply_tier1_correction",
                                         "apply_tiered_correction",
                                         "station_obs", "swe_template_recipe", "swe_matrix", "dates_vec",
                                         "region_id", "elevation", "aspect_deg", "slope_deg",
                                         "lon", "lat", "clip_range", "min_n_tier2", "min_n_tier3",
                                         "sar_proxy", "mon_all", "yr_all", "index_values_cv"),
                    envir = environment())
      cat(sprintf("[%s] [LOSO-CV setup] clusterExport() done (%.1fs). Loading terra + rebuilding template raster on each worker...\n",
                  format(Sys.time(), "%H:%M:%S"),
                  as.numeric(difftime(Sys.time(), export_t0, units = "secs"))))
      flush.console()
      
      eval_t0 <- Sys.time()
      clusterEvalQ(cl_loso, library(terra))
      # Rebuild a real (but pixel-empty) SpatRaster on EACH WORKER from the
      # plain-data recipe exported above, and bind it to the name `swe_template`
      # in that worker's .GlobalEnv - this is the object .loso_worker() and the
      # fit_*() helper functions actually reference, so no other code in this
      # script needs to change. Built once per worker (not once per task), so
      # this cost is paid n_cores_loso times total, not n_stations_loso times.
      clusterEvalQ(cl_loso, {
        swe_template <- terra::rast(
          xmin = swe_template_recipe$extent["xmin"], xmax = swe_template_recipe$extent["xmax"],
          ymin = swe_template_recipe$extent["ymin"], ymax = swe_template_recipe$extent["ymax"],
          nrows = swe_template_recipe$nrow, ncols = swe_template_recipe$ncol,
          crs = swe_template_recipe$crs_wkt
        )
        NULL
      })
      cat(sprintf("[%s] [LOSO-CV setup] Worker setup complete (%.1fs total: cluster+export+rebuild). Dispatching tasks...\n",
                  format(Sys.time(), "%H:%M:%S"),
                  as.numeric(difftime(Sys.time(), cluster_t0, units = "secs"))))
      flush.console()
      
      # ---- Manual load-balanced dispatch with a live progress readout ----
      # Same low-level sendCall()/recvOneResult() pattern used for the
      # heterogeneity tests earlier in the script, so a station taking several
      # minutes doesn't look like a hang - the master prints a line every time
      # one station's fold finishes rather than waiting on the whole batch.
      n_tasks   <- n_stations_loso
      n_workers <- length(cl_loso)
      station_out <- vector("list", n_tasks)
      loso_t0   <- Sys.time()
      
      n_initial <- min(n_workers, n_tasks)
      for (i in seq_len(n_initial)) {
        parallel:::sendCall(cl_loso[[i]], .loso_worker, list(stations[i]), tag = i)
      }
      next_task <- n_initial + 1L
      completed <- 0L
      
      tryCatch({
        while (completed < n_tasks) {
          res <- parallel:::recvOneResult(cl_loso)
          station_out[[res$tag]] <- res$value
          completed <- completed + 1L
          
          if (next_task <= n_tasks) {
            parallel:::sendCall(cl_loso[[res$node]], .loso_worker, list(stations[next_task]), tag = next_task)
            next_task <- next_task + 1L
          }
          
          s_val         <- station_out[[res$tag]]
          elapsed_total <- as.numeric(difftime(Sys.time(), loso_t0, units = "secs"))
          rate          <- completed / elapsed_total
          eta_sec       <- if (is.finite(rate) && rate > 0) (n_tasks - completed) / rate else 0
          if (isTRUE(s_val$skipped)) {
            cat(sprintf("[LOSO-CV] %4d/%4d done | station '%s' skipped (no train or test rows) | elapsed %s | ETA %s remaining\n",
                        completed, n_tasks, s_val$station, format_hms(elapsed_total), format_hms(eta_sec)))
          } else {
            cat(sprintf("[LOSO-CV] %4d/%4d done | station '%s' took %.1fs | elapsed %s | ETA %s remaining\n",
                        completed, n_tasks, s_val$station, s_val$elapsed, format_hms(elapsed_total), format_hms(eta_sec)))
          }
        }
      }, error = function(e) { try(stopCluster(cl_loso), silent = TRUE); stop(e) })
      
      stopCluster(cl_loso)
    }
  }
  
  results <- do.call(rbind, lapply(station_out, `[[`, "results"))
  if (is.null(results) || nrow(results) == 0) return(list(best = "tier1", table = data.frame()))
  rmse_tab <- aggregate(sq_err ~ tier, data = results,
                        FUN = function(x) sqrt(mean(x, na.rm = TRUE)))
  names(rmse_tab)[2] <- "rmse"
  n_tab <- aggregate(sq_err ~ tier, data = results, FUN = length)
  names(n_tab)[2] <- "n_holdout"
  rmse_tab <- merge(rmse_tab, n_tab, by = "tier")
  
  # Rank-fidelity metrics logged alongside RMSE (Fix #7). Tier SELECTION
  # below still uses RMSE only - this is diagnostic, not a logic change.
  rank_tab <- do.call(rbind, lapply(split(results, results$tier), function(df) {
    rho <- if (nrow(df) >= 4 && stats::sd(df$pred) > 0 && stats::sd(df$obs) > 0)
      suppressWarnings(stats::cor(df$pred, df$obs, method = "spearman")) else NA_real_
    tau <- if (nrow(df) >= 4 && stats::sd(df$pred) > 0 && stats::sd(df$obs) > 0)
      suppressWarnings(stats::cor(df$pred, df$obs, method = "kendall")) else NA_real_
    data.frame(tier = df$tier[1], spearman_rho = rho, kendall_tau = tau)
  }))
  rmse_tab <- merge(rmse_tab, rank_tab, by = "tier")
  
  best <- rmse_tab$tier[which.min(rmse_tab$rmse)]
  
  cat("LOSO-CV RMSE by tier (lower is better) - selection stays RMSE-based; rank\n")
  cat("metrics are logged alongside since SSPI-1 downstream is rank-based:\n")
  for (i in seq_len(nrow(rmse_tab))) {
    cat(sprintf("  %s: RMSE = %.2f mm | Spearman rho = %.3f | Kendall tau = %.3f (n = %d held-out obs)\n",
                rmse_tab$tier[i], rmse_tab$rmse[i], rmse_tab$spearman_rho[i], rmse_tab$kendall_tau[i], rmse_tab$n_holdout[i]))
  }
  cat(sprintf("-> Selected tier: %s\n", best))
  # station_results (one row per held-out obs x tier) is kept so
  # loso_cv_stability_check() below can bootstrap-resample it without
  # needing to re-fit any tier.
  list(best = best, table = rmse_tab, station_results = results)
}

## ==========================================================================
##   LOSO-CV BOUNDARY STABILITY CHECK
## ==========================================================================
## Addresses: "worth stress-testing LOSO_CV_MIN_STATIONS by running LOSO
## manually at fewer stations to see how unstable the tier ranking is right
## at the boundary". loso_cv_bias_tiers() already computed one held-out
## squared error per (station-observation, tier); this re-derives per-tier
## RMSE from BOOTSTRAP RESAMPLES of that same result set (resampling WHICH
## held-out rows are included, with replacement) rather than re-fitting
## anything - a cheap, honest way to ask "if the sample of stations near the
## LOSO_CV_MIN_STATIONS boundary had looked slightly different, would the
## SAME tier still win?", without the cost of re-running
## fit_regionalized_qm()/fit_mlr_correction() many times over.
## Only meaningful when n_stations is close to LOSO_CV_MIN_STATIONS - well
## above that threshold, tier selection is already reasonably well-supported
## and this check is skipped to avoid pointless computation.
loso_cv_stability_check <- function(loso_result, n_stations, min_stations_threshold,
                                    near_boundary_margin = 3, n_boot = 500, seed = 40) {
  if (is.null(loso_result) || is.null(loso_result$station_results) ||
      nrow(loso_result$station_results) == 0) {
    return(NULL)
  }
  if (n_stations > min_stations_threshold + near_boundary_margin) {
    cat(sprintf(
      "LOSO-CV stability check skipped: %d stations is comfortably above LOSO_CV_MIN_STATIONS (%d) + margin (%d).\n",
      n_stations, min_stations_threshold, near_boundary_margin))
    return(NULL)
  }
  
  cat(sprintf(
    "\n----- LOSO-CV BOUNDARY STABILITY CHECK (%d stations, near LOSO_CV_MIN_STATIONS = %d) -----\n",
    n_stations, min_stations_threshold))
  cat("Bootstrap-resampling the held-out station errors to see how often the SAME tier wins;\n")
  cat("a low win-rate for the nominally 'best' tier means the ranking is not well-supported at\n")
  cat("this station count, even though a single LOSO-CV pass picked a winner.\n")
  
  res <- loso_result$station_results
  set.seed(seed)
  
  n_rows <- nrow(res)
  boot_winners <- character(n_boot)
  for (b in seq_len(n_boot)) {
    idx <- sample.int(n_rows, n_rows, replace = TRUE)
    boot_res <- res[idx, ]
    boot_rmse <- aggregate(sq_err ~ tier, data = boot_res,
                           FUN = function(x) sqrt(mean(x, na.rm = TRUE)))
    boot_winners[b] <- if (nrow(boot_rmse) == 0) NA_character_ else boot_rmse$tier[which.min(boot_rmse$rmse)]
  }
  
  win_tab <- as.data.frame(table(boot_winners, useNA = "ifany"))
  names(win_tab) <- c("tier", "n_bootstrap_wins")
  win_tab$win_fraction <- win_tab$n_bootstrap_wins / n_boot
  win_tab <- win_tab[order(-win_tab$win_fraction), ]
  
  nominal_best <- loso_result$best
  nominal_best_win_frac <- win_tab$win_fraction[win_tab$tier == nominal_best]
  if (length(nominal_best_win_frac) == 0) nominal_best_win_frac <- 0
  
  for (i in seq_len(nrow(win_tab))) {
    cat(sprintf("  %s: won %d/%d bootstrap resamples (%.0f%%)\n",
                win_tab$tier[i], win_tab$n_bootstrap_wins[i], n_boot, 100 * win_tab$win_fraction[i]))
  }
  if (nominal_best_win_frac < 0.60) {
    cat(sprintf(paste0(
      "[CAUTION] The single LOSO-CV pass picked '%s' as best, but it only wins %.0f%% of bootstrap ",
      "resamples of the held-out errors - the tier ranking is NOT well-supported at this station count ",
      "(n=%d, near LOSO_CV_MIN_STATIONS=%d). Consider treating 'tier1_only' as the safer choice until ",
      "more stations are available, even though 'auto' nominally selected '%s' this run.\n"),
      nominal_best, 100 * nominal_best_win_frac, n_stations, min_stations_threshold, nominal_best))
  } else {
    cat(sprintf("[OK] '%s' wins %.0f%% of bootstrap resamples - reasonably stable given the available stations.\n",
                nominal_best, 100 * nominal_best_win_frac))
  }
  
  list(win_table = win_tab, nominal_best = nominal_best, nominal_best_win_fraction = nominal_best_win_frac,
       n_stations = n_stations, n_boot = n_boot)
}

## ==========================================================================
##   VALIDATION PROTOCOL (framework doc, Section 5)
## ==========================================================================
## Kanda & Fletcher (2025) recommend running ALL THREE transferability tests,
## since each exposes a different failure mode:
##   1. Random k-fold CV        - optimistic upper bound; compares tiers only.
##   2. Spatial-block CV        - leave-one-station-out; this IS
##                                 loso_cv_bias_tiers() above, reused here.
##   3. Temporal transfer       - train on early years / test on held-out
##                                 later years, and vice versa.
## Reports bias, RMSE, and R^2 per tier x scheme, so every run documents not
## just what correction was applied but how much it can be trusted - written
## to bias_correction_validation.csv alongside basin_wide_bias_correction_by_month.csv.
run_bias_correction_validation <- function(station_obs, swe_template, swe_matrix_uncorrected,
                                           dates_vec, region_id, bias_predictors,
                                           clip_range, min_n_tier2, min_n_tier3,
                                           loso_table = NULL, k_folds = 5) {
  out <- data.frame()
  if (!is.null(loso_table) && nrow(loso_table) > 0) {
    out <- rbind(out, data.frame(scheme = "spatial_block_loso", tier = loso_table$tier,
                                 rmse = loso_table$rmse, bias = NA_real_, r2 = NA_real_,
                                 spearman_rho = loso_table$spearman_rho, kendall_tau = loso_table$kendall_tau,
                                 n_holdout = loso_table$n_holdout))
  }
  if (is.null(station_obs) || nrow(station_obs) == 0) return(out)
  
  station_obs$date <- as.Date(station_obs$date)
  index_values_v <- compute_index_values(swe_matrix_uncorrected, dates_vec,
                                         seq_along(dates_vec), stat = INDEX_VALUE_STAT)
  mon_all <- as.integer(format(dates_vec, "%m")); yr_all <- as.integer(format(dates_vec, "%Y"))
  
  # ---- Progress reporting ----
  # Each call to eval_split() below refits Tier 1/2/3 (and Tier 2-SAR) from
  # scratch, exactly like one LOSO-CV station iteration - the k-fold loop
  # alone does this k_folds (default 5) times, plus 2 more for the temporal-
  # transfer split, so this validation protocol can itself take as long as
  # the LOSO-CV step. Report which split is running and how long it took.
  cat("\n[VALIDATION] Starting random k-fold + temporal-transfer validation protocol\n")
  cat("             (each split below refits Tier 1/2/3 from scratch).\n")
  validation_t0 <- Sys.time()
  
  ## Fits all three tiers on train_obs and scores them against test_obs.
  eval_split <- function(train_obs, test_obs, scheme_name) {
    split_t0 <- Sys.time()
    if (nrow(train_obs) == 0 || nrow(test_obs) == 0) return(data.frame())
    cat(sprintf("[VALIDATION] Fitting split '%s' (%d train / %d test obs) ... ",
                scheme_name, nrow(train_obs), nrow(test_obs)))
    bbc_v <- fit_basin_bias_correction(train_obs, swe_template, swe_matrix_uncorrected,
                                       dates_vec, clip_range = clip_range)
    qm_v  <- fit_regionalized_qm(train_obs, swe_template, swe_matrix_uncorrected, dates_vec,
                                 region_id, index_values_v, min_n = min_n_tier2,
                                 clip_range = clip_range)
    mlr_v <- fit_mlr_correction(train_obs, swe_template, swe_matrix_uncorrected, dates_vec,
                                bias_predictors$elevation, bias_predictors$aspect_deg,
                                bias_predictors$slope_deg, bias_predictors$lon, bias_predictors$lat,
                                min_n = min_n_tier3)
    mlr_sar_v <- if (!is.null(bias_predictors$sar_proxy)) fit_tier2_sar_correction(
      train_obs, swe_template, swe_matrix_uncorrected, dates_vec,
      bias_predictors$elevation, bias_predictors$aspect_deg, bias_predictors$slope_deg,
      bias_predictors$lon, bias_predictors$lat, min_n = min_n_tier3,
      sar_proxy = bias_predictors$sar_proxy
    ) else NULL
    pts <- vect(test_obs, geom = c("lon", "lat"), crs = crs(swe_template))
    test_obs$cell <- cells(swe_template, pts)[, "cell"]
    
    preds <- data.frame()
    for (i in seq_len(nrow(test_obs))) {
      cix <- test_obs$cell[i]; if (is.na(cix)) next
      d <- test_obs$date[i]
      col <- which(mon_all == as.integer(format(d, "%m")) & yr_all == as.integer(format(d, "%Y")))
      if (length(col) != 1) next
      era_val <- swe_matrix_uncorrected[cix, col]; obs_val <- test_obs$swe_mm[i]
      if (!is.finite(era_val) || !is.finite(obs_val)) next
      m <- as.integer(format(d, "%m"))
      
      p1 <- apply_tier1_correction(era_val, m, bbc_v)
      p2 <- apply_tiered_correction(era_val, region_id[cix], qm_v$C, bbc_v, m)
      preds <- rbind(preds, data.frame(tier = "tier1", pred = p1, obs = obs_val),
                     data.frame(tier = "tier2", pred = p2, obs = obs_val))
      if (!is.null(mlr_v)) {
        p3 <- predict_mlr_grid(mlr_v, era_val, bias_predictors$elevation[cix], bias_predictors$slope_deg[cix],
                               bias_predictors$aspect_deg[cix], bias_predictors$lon[cix], bias_predictors$lat[cix], d)
        if (is.finite(p3)) preds <- rbind(preds, data.frame(tier = "tier3", pred = p3, obs = obs_val))
      }
      if (!is.null(mlr_sar_v)) {
        sar_col <- if (!is.null(bias_predictors$sar_proxy)) bias_predictors$sar_proxy[cix, which(mon_all == as.integer(format(d, "%m")) & yr_all == as.integer(format(d, "%Y")))] else NULL
        p4 <- predict_mlr_grid(mlr_sar_v, era_val, bias_predictors$elevation[cix], bias_predictors$slope_deg[cix],
                               bias_predictors$aspect_deg[cix], bias_predictors$lon[cix], bias_predictors$lat[cix], d,
                               sar_proxy = sar_col)
        if (is.finite(p4)) preds <- rbind(preds, data.frame(tier = "tier2_sar", pred = p4, obs = obs_val))
      }
    }
    if (nrow(preds) == 0) {
      cat(sprintf("no usable predictions (%.1fs)\n", as.numeric(difftime(Sys.time(), split_t0, units = "secs"))))
      return(data.frame())
    }
    cat(sprintf("done in %.1fs\n", as.numeric(difftime(Sys.time(), split_t0, units = "secs"))))
    do.call(rbind, lapply(split(preds, preds$tier), function(df) {
      err <- df$pred - df$obs
      r2 <- if (stats::sd(df$obs) > 0) 1 - sum(err^2) / sum((df$obs - mean(df$obs))^2) else NA_real_
      # Rank-fidelity metrics (Fix #7): SSPI-1 is a rank-based transform
      # (Gringorten / gamma CDF), so tier selection by RMSE alone can favor
      # a "smoother" tier that quietly compresses the distribution's tails.
      # Logged alongside RMSE using the SAME held-out predictions - no
      # extra refits, just two more cor() calls per fold.
      rho <- if (nrow(df) >= 4 && stats::sd(df$pred) > 0 && stats::sd(df$obs) > 0)
        suppressWarnings(stats::cor(df$pred, df$obs, method = "spearman")) else NA_real_
      tau <- if (nrow(df) >= 4 && stats::sd(df$pred) > 0 && stats::sd(df$obs) > 0)
        suppressWarnings(stats::cor(df$pred, df$obs, method = "kendall")) else NA_real_
      data.frame(scheme = scheme_name, tier = df$tier[1], rmse = sqrt(mean(err^2)),
                 bias = mean(err), r2 = r2, spearman_rho = rho, kendall_tau = tau,
                 n_holdout = nrow(df))
    }))
  }
  
  ## ---- 1. Random k-fold CV (row-wise, on station-month observations) ----
  set.seed(42)
  n_obs <- nrow(station_obs)
  if (n_obs >= k_folds) {
    folds <- sample(rep(1:k_folds, length.out = n_obs))
    for (k in 1:k_folds) {
      cat(sprintf("[VALIDATION] Random k-fold CV: fold %d/%d\n", k, k_folds))
      out <- rbind(out, eval_split(station_obs[folds != k, ], station_obs[folds == k, ],
                                   sprintf("random_kfold_%d", k)))
    }
  }
  
  ## ---- 3. Temporal transfer: early-half train / late-half test, and vice versa ----
  yrs <- sort(unique(as.integer(format(station_obs$date, "%Y"))))
  if (length(yrs) >= 4) {
    mid <- yrs[ceiling(length(yrs) / 2)]
    early <- station_obs[as.integer(format(station_obs$date, "%Y")) <= mid, ]
    late  <- station_obs[as.integer(format(station_obs$date, "%Y")) >  mid, ]
    cat("[VALIDATION] Temporal-transfer CV: early-to-late\n")
    out <- rbind(out, eval_split(early, late, "temporal_transfer_early_to_late"))
    cat("[VALIDATION] Temporal-transfer CV: late-to-early\n")
    out <- rbind(out, eval_split(late, early, "temporal_transfer_late_to_early"))
  } else {
    cat("Temporal-transfer CV: fewer than 4 distinct station years - skipping (need an early/half split).\n")
  }
  
  cat(sprintf("[VALIDATION] Validation protocol complete in %s.\n",
              format_hms(as.numeric(difftime(Sys.time(), validation_t0, units = "secs")))))
  out
}

gringorten_swei_seasonal <- function(v_smoothed, dates_vec, scf_mask_list,
                                     cell_index = NULL,
                                     basin_scf_mask = NULL,
                                     eps = 1e-6) {
  if (length(v_smoothed) == 0 || all(is.na(v_smoothed)) || is.null(v_smoothed)) {
    return(list(swei = rep(NA_real_, length(dates_vec)), method = rep(NA_integer_, length(dates_vec))))
  }
  v_clean <- v_smoothed
  v_clean[!is.finite(v_clean)] <- NA_real_
  n <- length(v_clean)
  z <- rep(NA_real_, n)
  method_used <- rep(NA_integer_, n) 
  mon <- as.integer(format(dates_vec, "%m"))
  
  for (m in 1:12) {
    all_m_idx <- which(mon == m)
    if (length(all_m_idx) == 0) next
    nonzero_frac <- sum(v_clean[all_m_idx] > 0, na.rm = TRUE) / length(all_m_idx)
    if (nonzero_frac < 0.75) next
    
    idx <- which(mon == m & is.finite(v_clean))
    if (length(idx) == 0) next
    
    if (is.null(scf_mask_list) || length(scf_mask_list) < m) next
    if (!is.null(cell_index)) {
      mask_vec <- scf_mask_list[[m]]
      if (length(mask_vec) < cell_index || !mask_vec[cell_index]) next
    } else if (!is.null(basin_scf_mask)) {
      if (length(basin_scf_mask) < m || !basin_scf_mask[m]) next
    }
    
    samp <- v_clean[idx]
    zero_idx <- which(samp == 0 & !is.na(samp))
    if (length(zero_idx) > 0) {
      nonzero_vals <- samp[samp > 0 & !is.na(samp)]
      if (length(nonzero_vals) > 0) {
        min_nonzero <- min(nonzero_vals)
        samp[zero_idx] <- runif(length(zero_idx), min = 0, max = min_nonzero)
      }
    }
    valid_idx <- which(!is.na(samp))
    if (length(valid_idx) < 3) next
    
    samp_valid <- samp[valid_idx]
    idx_valid <- idx[valid_idx]
    
    r <- rank(samp_valid, ties.method = "average")
    N <- length(samp_valid)
    p_val <- (r - 0.44) / (N + 0.12)
    p_val <- clip_prob(p_val, eps = eps)
    z_val <- qnorm(p_val)
    
    z[idx_valid] <- z_val
    method_used[idx_valid] <- 1
  }
  list(swei = z, method = method_used)
}

monthly_sspi <- function(v, dates_vec, ref_idx, eps = 1e-6) {
  n <- length(v)
  if (n == 0 || all(is.na(v)) || is.null(v)) {
    return(list(sspi = rep(NA_real_, n), method = rep(NA_integer_, n)))
  }
  v_clean <- v
  v_clean[!is.finite(v_clean)] <- NA_real_
  z <- rep(NA_real_, n)
  method_used <- rep(NA_integer_, n) 
  mon <- as.integer(format(dates_vec, "%m"))
  
  for (m in 1:12) {
    ref_mon_idx <- which(mon == m & seq_along(dates_vec) %in% ref_idx)
    ref_samp <- v_clean[ref_mon_idx]
    ref_samp <- ref_samp[is.finite(ref_samp)]
    obs_mon_idx <- which(mon == m)
    
    if (length(ref_samp) < MIN_PIXEL_GAMMA_SAMPLE) {
      if (length(obs_mon_idx) >= 3) {
        obs_vals <- v_clean[obs_mon_idx]
        fin_idx  <- is.finite(obs_vals)
        n_fin    <- sum(fin_idx)
        if (n_fin >= 3) {
          p <- (rank(obs_vals[fin_idx], ties.method = "average") - 0.44) / (n_fin + 0.12)
          p <- clip_prob(p, eps = eps)
          z[obs_mon_idx[fin_idx]] <- qnorm(p)
          method_used[obs_mon_idx[fin_idx]] <- 2
        }
      }
      next
    }
    ref_var <- var(ref_samp, na.rm = TRUE)
    if (!is.finite(ref_var) || ref_var < .Machine$double.eps) {
      z[obs_mon_idx] <- 0
      method_used[obs_mon_idx] <- 3
      next
    }
    if (ref_var < 0.01) {
      obs_vals <- v_clean[obs_mon_idx]
      fin_idx <- is.finite(obs_vals)
      if (sum(fin_idx) >= 3) {
        p <- (rank(obs_vals[fin_idx], ties.method = "average") - 0.44) /
          (sum(fin_idx) + 0.12)
        p <- clip_prob(p, eps = eps)
        z[obs_mon_idx[fin_idx]] <- qnorm(p)
        method_used[obs_mon_idx[fin_idx]] <- 2
      }
      next
    }
    p0 <- mean(ref_samp <= 0, na.rm = TRUE)
    ref_pos <- ref_samp[ref_samp > 0]
    if (length(ref_pos) >= 10) {
      lm_obj <- try(lmomco::lmoms(ref_pos), silent = TRUE)
      if (!inherits(lm_obj, "try-error") && !is.null(lm_obj)) {
        par_obj <- try(lmomco::pargam(lm_obj), silent = TRUE)
        if (!inherits(par_obj, "try-error") && !is.null(par_obj)) {
          obs_vals <- v_clean[obs_mon_idx]
          p_m <- rep(NA_real_, length(obs_vals))
          pos_idx <- which(is.finite(obs_vals) & obs_vals > 0)
          if (length(pos_idx) > 0) {
            Fg <- try(lmomco::cdfgam(obs_vals[pos_idx], par_obj), silent = TRUE)
            if (!inherits(Fg, "try-error")) p_m[pos_idx] <- p0 + (1 - p0) * as.numeric(Fg)
          }
          zero_idx2 <- which(is.finite(obs_vals) & obs_vals <= 0)
          if (length(zero_idx2) > 0) p_m[zero_idx2] <- p0
          p_m <- clip_prob(p_m, eps = eps)
          z[obs_mon_idx] <- qnorm(p_m)
          method_used[obs_mon_idx] <- 1
          next
        }
      }
    }
    all_samp <- c(ref_samp, v_clean[obs_mon_idx])
    all_samp <- all_samp[is.finite(all_samp)]
    obs_vals <- v_clean[obs_mon_idx]
    fin_obs <- is.finite(obs_vals)
    if (length(all_samp) >= 10 && sum(fin_obs) > 0) {
      combined <- c(all_samp, obs_vals[fin_obs])
      r_all <- rank(combined, ties.method = "average")
      r_obs <- r_all[(length(all_samp) + 1):(length(all_samp) + sum(fin_obs))]
      p <- (r_obs - 0.44) / (length(combined) + 0.12)
      p <- clip_prob(p, eps = eps)
      z[obs_mon_idx[fin_obs]] <- qnorm(p)
      method_used[obs_mon_idx[fin_obs]] <- 2
    }
  }
  fin_z <- is.finite(z)
  z[fin_z & z < -5] <- -5
  z[fin_z & z >  5] <-  5
  list(sspi = z, method = method_used)
}

## ---- 8. Regionalized SSPI-1 (uses pooled regional gamma / empirical fit) ----
# Method codes: 1=regional gamma (GoF pass), 2=regional pooled empirical
# (GoF fail, still regional), 3=local gamma fallback (no usable region),
# 4=local empirical fallback, 5=local zero-variance fallback.
# Falls back lazily to the original monthly_sspi() (per-pixel) whenever the
# regional fit is unusable for a given pixel/month (no region assigned,
# missing/zero index value, discordant pixel, or insufficient pooled sample).
monthly_sspi_regional <- function(v, dates_vec, ref_idx, pixel_idx,
                                  region_vec, index_values, C, region_par,
                                  discordant_mat, eps = 1e-6) {
  n <- length(v)
  z <- rep(NA_real_, n)
  method_used <- rep(NA_integer_, n)
  mon <- as.integer(format(dates_vec, "%m"))
  v_clean <- v
  v_clean[!is.finite(v_clean)] <- NA_real_
  
  region_here <- if (pixel_idx <= length(region_vec)) region_vec[pixel_idx] else NA_integer_
  old_result <- NULL  # lazily computed pixel-local fallback (only if actually needed)
  
  for (m in 1:12) {
    obs_mon_idx <- which(mon == m)
    if (length(obs_mon_idx) == 0) next
    
    mu <- if (!is.na(region_here)) index_values[pixel_idx, m] else NA_real_
    disc_flag <- if (!is.null(discordant_mat)) isTRUE(discordant_mat[pixel_idx, m]) else FALSE
    par_info <- if (!is.na(region_here) && !is.null(region_par[[region_here]])) {
      region_par[[region_here]][[m]]
    } else NULL
    # pool_ok now additionally requires gof_pass == TRUE. Region/months that
    # failed the KS goodness-of-fit, OR were excluded from pooling upstream
    # for failing the H0 heterogeneity screen (region_par[[r]][[m]]$status ==
    # "heterogeneous_skip", set in the RFA construction block), have
    # par = NULL / gof_pass = FALSE and so never qualify here. There is
    # deliberately no regional-empirical fallback: a region/month either uses
    # the pooled regional gamma, or falls straight through to per-pixel local
    # fitting below (monthly_sspi) - see prior discussion on why a pooled
    # empirical CDF was dropped as an intermediate option.
    pool_ok <- !is.null(par_info) && isTRUE(par_info$gof_pass) && !is.null(par_info$par) &&
      !is.null(par_info$pooled_sample) && length(par_info$pooled_sample) >= 15
    use_regional <- !is.na(region_here) && is.finite(mu) && mu > 0 && !disc_flag && pool_ok
    
    if (use_regional) {
      c_val <- C[region_here, m]
      mu_corr <- mu * c_val
      x_star <- (v_clean[obs_mon_idx] * c_val) / mu_corr
      
      ref_mon_idx <- which(mon == m & seq_along(dates_vec) %in% ref_idx)
      ref_samp <- v_clean[ref_mon_idx]
      ref_samp <- ref_samp[is.finite(ref_samp)]
      p0 <- if (length(ref_samp) >= 5) mean(ref_samp <= 0, na.rm = TRUE) else mean(v_clean[obs_mon_idx] <= 0, na.rm = TRUE)
      if (!is.finite(p0)) p0 <- 0
      
      p_reg <- rep(NA_real_, length(obs_mon_idx))
      pos_idx <- which(is.finite(x_star) & x_star > 0)
      if (length(pos_idx) > 0) {
        Fg <- try(lmomco::cdfgam(x_star[pos_idx], par_info$par), silent = TRUE)
        if (!inherits(Fg, "try-error")) p_reg[pos_idx] <- p0 + (1 - p0) * as.numeric(Fg)
      }
      zero_idx2 <- which(is.finite(x_star) & x_star <= 0)
      if (length(zero_idx2) > 0) p_reg[zero_idx2] <- p0
      
      # Soft H0 shrinkage (Fix #5): blend with the local per-pixel fit,
      # weighted by shrink_weight in [0,1]. Replaces the old hard cliff at
      # H0_REGION_FLAG_THRESHOLD. Blend happens in PROBABILITY space (bounded
      # [0,1]) rather than z-space, so it can't produce an out-of-range value.
      w <- if (!is.null(par_info$shrink_weight)) par_info$shrink_weight else 1
      
      if (w >= 1 - 1e-9) {
        p_m <- clip_prob(p_reg, eps = eps)
        z[obs_mon_idx] <- qnorm(p_m)
        method_used[obs_mon_idx] <- 1L
      } else {
        if (is.null(old_result)) {
          old_result <- monthly_sspi(v, dates_vec, ref_idx, eps = eps)
        }
        p_local <- pnorm(old_result$sspi[obs_mon_idx])
        p_reg_filled <- ifelse(is.finite(p_reg), p_reg, p_local)
        p_blend <- clip_prob(w * p_reg_filled + (1 - w) * p_local, eps = eps)
        z[obs_mon_idx] <- qnorm(p_blend)
        method_used[obs_mon_idx] <- 6L   # NEW: blended regional/local
      }
      next
    }
    
    # ---- Fallback: no usable regional fit for this pixel/month ----
    if (is.null(old_result)) {
      old_result <- monthly_sspi(v, dates_vec, ref_idx, eps = eps)
    }
    z[obs_mon_idx] <- old_result$sspi[obs_mon_idx]
    old_codes <- old_result$method[obs_mon_idx]   # 1=gamma,2=empirical,3=zero-var (local)
    method_used[obs_mon_idx] <- ifelse(is.na(old_codes), NA_integer_, old_codes + 2L)
  }
  
  fin_z <- is.finite(z)
  z[fin_z & z < -5] <- -5
  z[fin_z & z >  5] <-  5
  list(sspi = z, method = method_used)
}

# ==============================================================================
#   BASIN-WIDE BIAS CORRECTION AGAINST STATION SWE OBSERVATIONS  [method 2 only]
# ==============================================================================
# Runs ONCE, here - after fit_basin_bias_correction() has been defined above,
# and BEFORE the RFA construction block / main calculation loop use swe_matrix
# for anything. Applying it at this point (rather than deep inside the
# REGIONALIZE_SSPI==TRUE branch) means BOTH the regionalized path
# (monthly_sspi_regional, via build_regions/compute_index_values/the pooled
# gamma fits) and the pixel-by-pixel path (monthly_sspi) see the same
# corrected swe_matrix, regardless of REGIONALIZE_SSPI.
#
# This is the ONLY place station correction is applied - see Section F inside
# the RFA construction block below, which intentionally no-ops the old
# per-region fit_bias_correction() to avoid double-correcting.
if (scf_method == "2") {
  cat("\n===== BASIN-WIDE BIAS CORRECTION vs STATION OBSERVATIONS =====\n")
  cat("This step can take a while (station QA -> regionalization -> Tier 1/2/3 fits ->\n")
  cat("LOSO-CV -> validation protocol). Progress and elapsed/ETA times are reported\n")
  cat("below at each stage so a long run doesn't look stalled.\n")
  bc_step_t0 <- Sys.time()
  bc_stage <- function(msg) {
    cat(sprintf("\n[BIAS-CORRECT %s elapsed] %s\n",
                format_hms(as.numeric(difftime(Sys.time(), bc_step_t0, units = "secs"))), msg))
  }
  
  bc_diag_dir <- file.path(out_dir, "rfa_diagnostics")
  
  if (BIAS_CORRECT && (file.exists(STATION_OBS_PATH) || file.exists(CANSWE_STATION_OBS_PATH))) {
    if (!dir.exists(bc_diag_dir)) dir.create(bc_diag_dir, recursive = TRUE)
    
    bc_stage("Loading and merging station SWE observation network(s)...")
    station_obs_raw <- merge_station_networks(STATION_OBS_PATH, CANSWE_STATION_OBS_PATH,
                                              swe_template = swe[[1]])
    if (is.null(station_obs_raw) || nrow(station_obs_raw) == 0) {
      log_pipeline_fallback(
        sprintf("No usable station SWE observations found at '%s' or '%s' - bias correction will be skipped.",
                STATION_OBS_PATH, CANSWE_STATION_OBS_PATH),
        hard_stop = isTRUE(REQUIRE_STATION_OBS)
      )
      station_obs_raw <- NULL
    } else {
      cat(sprintf("  -> Loaded %d station-observation row(s) from %d distinct station(s).\n",
                  nrow(station_obs_raw), length(unique(station_obs_raw$station_id))))
    }
    
    ## ---- Tier 0: QA screen (Z-score outlier removal on the station/ERA5-
    ## Land ratio, not on raw SWE - see framework doc, Tier 0) ----
    bc_stage("Tier 0: QA-screening station observations for station/ERA5-Land outliers...")
    station_obs_raw <- qa_screen_station_obs(station_obs_raw, swe[[1]], swe_matrix, dates,
                                             z_thresh = 3, out_dir = bc_diag_dir)
    
    ## ---- Predictors shared by Tier 2 (regionalized QM) and Tier 3 (MLR) -
    ## reuses build_regions() (via get_bias_predictors()) so bias correction
    ## draws on the SAME elevation x aspect regions as SSPI pooling, per the
    ## framework doc's Tier 2 design, independent of REGIONALIZE_SSPI. ----
    bc_stage("Building elevation x aspect regions and bias-correction predictor fields\n            (DEM crop/resample/terrain + optional SAR proxy load - can take a minute)...")
    bias_predictors <- get_bias_predictors(swe[[1]], basin, DEM_PATH,
                                           n_elev_bands = N_ELEV_BANDS, use_aspect = USE_ASPECT,
                                           n_aspect_classes = N_ASPECT_CLASSES,
                                           n_regions_fallback = N_REGIONS_FALLBACK,
                                           min_region_pixels = MIN_REGION_PIXELS,
                                           n_time = ncol(swe_matrix), sar_path = SAR_DEPTH_PATH)
    region_id_bc <- bias_predictors$region_id
    
    ## ---- Tier 1 (always fit - every other tier falls back to it) ----
    bc_stage("Tier 1: fitting basin-wide monthly bias correction (mult/add/SLR)...")
    bbc <- fit_basin_bias_correction(station_obs_raw, swe[[1]], swe_matrix, dates,
                                     clip_range = BIAS_CORRECT_CLIP)
    index_values_bc <- compute_index_values(swe_matrix, dates, seq_along(dates), stat = INDEX_VALUE_STAT)
    
    n_stations_bc <- if (is.null(station_obs_raw)) 0 else length(unique(station_obs_raw$station_id))
    tier_mode <- BIAS_TIER_MODE
    if (is.null(station_obs_raw) || n_stations_bc == 0) {
      tier_mode <- "tier1_only"
    } else if (tier_mode == "auto" && n_stations_bc < LOSO_CV_MIN_STATIONS) {
      cat(sprintf("Only %d distinct station(s) (< LOSO_CV_MIN_STATIONS = %d): forcing BIAS_TIER_MODE to 'tier1_only' for this run - too few stations for a meaningful held-out test.\n",
                  n_stations_bc, LOSO_CV_MIN_STATIONS))
      tier_mode <- "tier1_only"
    }
    
    qm <- NULL; mlr <- NULL; mlr_sar <- NULL; loso <- NULL
    if (tier_mode %in% c("tier2", "tier3", "tier2_sar", "auto")) {
      bc_stage("Tier 2: fitting regionalized polynomial CDF-matching (per region x month)...")
      qm <- fit_regionalized_qm(station_obs_raw, swe[[1]], swe_matrix, dates,
                                region_id_bc, index_values_bc, min_n = MIN_STATION_MONTHS_TIER2,
                                clip_range = BIAS_CORRECT_CLIP)
    }
    if (tier_mode %in% c("tier3", "tier2_sar", "auto")) {
      bc_stage("Tier 3: fitting multiple linear regression (MLR) bias correction...")
      mlr <- fit_mlr_correction(station_obs_raw, swe[[1]], swe_matrix, dates,
                                bias_predictors$elevation, bias_predictors$aspect_deg,
                                bias_predictors$slope_deg, bias_predictors$lon, bias_predictors$lat,
                                min_n = MIN_STATION_MONTHS_TIER3)
      if (!is.null(bias_predictors$sar_proxy)) {
        bc_stage("Tier 2-SAR: fitting SAR-augmented MLR bias correction...")
      }
      mlr_sar <- if (!is.null(bias_predictors$sar_proxy)) fit_tier2_sar_correction(
        station_obs_raw, swe[[1]], swe_matrix, dates,
        bias_predictors$elevation, bias_predictors$aspect_deg,
        bias_predictors$slope_deg, bias_predictors$lon, bias_predictors$lat,
        min_n = MIN_STATION_MONTHS_TIER3, sar_proxy = bias_predictors$sar_proxy
      ) else NULL
    }
    if (tier_mode == "auto") {
      ## ---- Reduce the LOSO-CV held-out station set before running it ----
      ## See LOSO_CV_MAX_STATIONS config comment and build_loso_holdout_stations()
      ## for the full rationale: with 790 pixels / 17 regions but potentially
      ## thousands of distinct station_id values, exhaustively holding out
      ## every station individually is a validation-design mismatch, not just
      ## slow. Tier 1/2/3 are still fit on the FULL station network in every
      ## fold; only the set of stations actually pulled out and scored shrinks.
      bc_stage("Selecting LOSO-CV held-out station set (pixel dedup + region-stratified sampling)...")
      loso_holdout_stations <- build_loso_holdout_stations(
        station_obs_raw, swe[[1]], region_id_bc,
        max_stations = LOSO_CV_MAX_STATIONS, min_per_region = LOSO_CV_MIN_PER_REGION)
      n_loso_holdout <- length(loso_holdout_stations)
      write.csv(data.frame(station_id = loso_holdout_stations),
                file.path(bc_diag_dir, "loso_cv_holdout_stations_used.csv"), row.names = FALSE)
      
      bc_stage(sprintf(
        "Running LOSO-CV across %d held-out station(s) (of %d total; see loso_cv_holdout_stations_used.csv)\n            to select the best tier - this refits Tier 1/2/3 once PER held-out station and is\n            typically the slowest step - see per-station progress below)...",
        n_loso_holdout, n_stations_bc))
      loso <- loso_cv_bias_tiers(station_obs_raw, swe[[1]], swe_matrix, dates,
                                 region_id_bc, bias_predictors$elevation, bias_predictors$aspect_deg,
                                 bias_predictors$slope_deg, bias_predictors$lon, bias_predictors$lat,
                                 clip_range = BIAS_CORRECT_CLIP,
                                 min_n_tier2 = MIN_STATION_MONTHS_TIER2, min_n_tier3 = MIN_STATION_MONTHS_TIER3,
                                 sar_proxy = bias_predictors$sar_proxy,
                                 holdout_stations = loso_holdout_stations)
      chosen_tier <- loso$best
      if (!is.null(loso$table) && nrow(loso$table) > 0) {
        write.csv(loso$table, file.path(bc_diag_dir, "bias_correction_loso_cv_by_tier.csv"), row.names = FALSE)
      }
      
      # Stress-test the tier ranking near the LOSO_CV_MIN_STATIONS boundary
      # (skips automatically, with a message, once well above threshold).
      # NOTE: uses n_loso_holdout (the ACTUAL number of folds LOSO-CV just
      # ran), not n_stations_bc (the full, pre-subsampling station count) -
      # the boundary-stability question is about how many held-out folds the
      # tier ranking is resting on, which is n_loso_holdout after subsampling.
      bc_stage("Bootstrap-resampling LOSO-CV results to stress-test tier ranking stability...")
      loso_stability <- loso_cv_stability_check(loso, n_loso_holdout, LOSO_CV_MIN_STATIONS)
      if (!is.null(loso_stability)) {
        write.csv(loso_stability$win_table,
                  file.path(bc_diag_dir, "bias_correction_loso_cv_stability.csv"), row.names = FALSE)
      }
    } else {
      chosen_tier <- switch(tier_mode, tier1_only = "tier1", tier2 = "tier2", tier2_sar = "tier2_sar", tier3 = "tier3", "tier1")
    }
    
    # Fall back gracefully if the selected/forced tier didn't actually fit
    # (mirrors the "fall back to simpler fit when the fancier one fails"
    # logic already used for the regional-gamma vs. local-gamma choice).
    if (chosen_tier == "tier3" && is.null(mlr)) chosen_tier <- "tier2"
    if (chosen_tier == "tier2_sar" && is.null(mlr_sar)) chosen_tier <- "tier2"
    if (chosen_tier == "tier2" && (is.null(qm) || nrow(qm$C) == 0)) chosen_tier <- "tier1"
    cat(sprintf("==> Bias-correction tier applied to the full grid: %s\n", chosen_tier))
    
    ## ---- Apply the chosen tier basin-wide, BEFORE the RFA/SSPI blocks use
    ## swe_matrix for anything (same placement as before). Keep an
    ## uncorrected copy for the validation protocol below. ----
    bc_stage(sprintf("Applying chosen tier ('%s') to every pixel/month across the full grid...", chosen_tier))
    apply_t0 <- Sys.time()
    swe_matrix_pre_bc <- swe_matrix
    months_vec_bc <- as.integer(format(dates, "%m"))
    for (m in 1:12) {
      cols_m <- which(months_vec_bc == m)
      if (length(cols_m) == 0) next
      if (chosen_tier %in% c("tier2", "tier3", "tier2_sar")) {
        cat(sprintf("  Applying %s: month %2d/12 (%s, %d time-step column(s))...\n",
                    chosen_tier, m, month.abb[m], length(cols_m)))
      }
      if (chosen_tier == "tier1") {
        swe_matrix[, cols_m] <- apply_tier1_correction(swe_matrix[, cols_m], m, bbc)
      } else if (chosen_tier == "tier2") {
        for (cc in cols_m) {
          swe_matrix[, cc] <- apply_tiered_correction(swe_matrix[, cc], region_id_bc, qm$C, bbc, m)
        }
      } else if (chosen_tier == "tier3") {
        for (cc in cols_m) {
          pred <- predict_mlr_grid(mlr, swe_matrix[, cc], bias_predictors$elevation,
                                   bias_predictors$slope_deg, bias_predictors$aspect_deg,
                                   bias_predictors$lon, bias_predictors$lat, dates[cc])
          bad <- !is.finite(pred)
          if (any(bad)) {
            pred[bad] <- apply_tiered_correction(swe_matrix[bad, cc], region_id_bc[bad], qm$C, bbc, m)
          }
          swe_matrix[, cc] <- pred
        }
      } else if (chosen_tier == "tier2_sar" && !is.null(mlr_sar)) {
        for (cc in cols_m) {
          sar_col <- if (!is.null(bias_predictors$sar_proxy)) bias_predictors$sar_proxy[, cc] else NULL
          pred <- predict_mlr_grid(mlr_sar, swe_matrix[, cc], bias_predictors$elevation,
                                   bias_predictors$slope_deg, bias_predictors$aspect_deg,
                                   bias_predictors$lon, bias_predictors$lat, dates[cc],
                                   sar_proxy = sar_col)
          bad <- !is.finite(pred)
          if (any(bad)) {
            pred[bad] <- apply_tiered_correction(swe_matrix[bad, cc], region_id_bc[bad], qm$C, bbc, m)
          }
          swe_matrix[, cc] <- pred
        }
      }
    }
    cat(sprintf("  -> Tier '%s' applied to full grid in %s.\n",
                chosen_tier, format_hms(as.numeric(difftime(Sys.time(), apply_t0, units = "secs")))))
    
    basin_bias_factor <- bbc$factor   
    basin_bias_log     <- bbc$details
    if (nrow(basin_bias_log) > 0) {
      # Dry-season coefficients frequently hit the hard bounds because SWE is
      # near zero and scalar correction is ill-conditioned. Flag these rows
      # explicitly so clipped values are not interpreted as estimated effects.
      tol_clip <- sqrt(.Machine$double.eps)
      basin_bias_log$factor_at_clip <- abs(basin_bias_log$factor - BIAS_CORRECT_CLIP[1]) <= tol_clip |
        abs(basin_bias_log$factor - BIAS_CORRECT_CLIP[2]) <= tol_clip
      basin_bias_log$b1_at_clip <- is.finite(basin_bias_log$b1) &
        (abs(basin_bias_log$b1 - BIAS_CORRECT_CLIP[1]) <= tol_clip |
           abs(basin_bias_log$b1 - BIAS_CORRECT_CLIP[2]) <= tol_clip)
      dry_clip <- basin_bias_log$month %in% 7:10 & (basin_bias_log$factor_at_clip | basin_bias_log$b1_at_clip)
      if (any(dry_clip, na.rm = TRUE)) {
        warning("Bias-correction coefficients hit hard clip bounds in Jul-Oct; treat those dry-season coefficients as constraint artifacts, not well-identified bias estimates.", call. = FALSE)
      }
      write.csv(basin_bias_log, file.path(bc_diag_dir, "basin_wide_bias_correction_by_month.csv"),
                row.names = FALSE)
    }
    if (!is.null(qm) && nrow(qm$details) > 0) {
      write.csv(qm$details, file.path(bc_diag_dir, "regionalized_qm_bias_correction_by_region_month.csv"),
                row.names = FALSE)
    }
    
    ## ---- CROSS-RUN TIER-USAGE HISTORY (addresses "does Tier 2/3 actually
    ## trigger in practice, or is BIAS_TIER_MODE == 'auto' silently falling
    ## back to Tier 1 every run because the station network never clears
    ## LOSO_CV_MIN_STATIONS / MIN_STATION_MONTHS_TIER2 / _TIER3?") ----
    ## A single run's console output only shows "Bias-correction tier applied:
    ## <x>" once and is easy to lose track of across many runs over a season.
    ## This appends one row per run to a persistent CSV (never overwritten) so
    ## that after N runs you can literally count what fraction resolved to
    ## Tier 1 vs Tier 2/2_sar/3, rather than assuming the tiered framework is
    ## doing something it isn't. Written even when chosen_tier == "tier1" so
    ## an all-Tier-1 history is visible, not just absent.
    tier_history_path <- file.path(bc_diag_dir, "bias_correction_tier_history.csv")
    loso_rmse_summary <- if (!is.null(loso) && !is.null(loso$table) && nrow(loso$table) > 0) {
      paste(sprintf("%s=%.2f(n=%d)", loso$table$tier, loso$table$rmse, loso$table$n_holdout), collapse = "; ")
    } else {
      NA_character_
    }
    loso_stability_summary <- if (exists("loso_stability", inherits = FALSE) && !is.null(loso_stability)) {
      sprintf("%s won %.0f%% of %d bootstrap resamples", loso_stability$nominal_best,
              100 * loso_stability$nominal_best_win_fraction, loso_stability$n_boot)
    } else {
      NA_character_
    }
    tier_history_row <- data.frame(
      run_timestamp        = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      n_stations           = n_stations_bc,
      n_loso_holdout       = if (exists("n_loso_holdout", inherits = FALSE)) n_loso_holdout else NA_integer_,
      loso_cv_max_stations_config = LOSO_CV_MAX_STATIONS,
      bias_tier_mode_config = BIAS_TIER_MODE,
      tier_mode_effective  = tier_mode,
      forced_tier1_only    = (BIAS_TIER_MODE == "auto" && tier_mode == "tier1_only" && n_stations_bc > 0),
      chosen_tier          = chosen_tier,
      tier2_fit_available  = !is.null(qm) && nrow(qm$C) > 0,
      tier3_fit_available  = !is.null(mlr),
      tier2_sar_available  = !is.null(mlr_sar),
      loso_cv_ran          = tier_mode == "auto",
      loso_rmse_by_tier    = loso_rmse_summary,
      loso_stability_check = loso_stability_summary,
      stringsAsFactors = FALSE
    )
    write.table(tier_history_row, tier_history_path, sep = ",", row.names = FALSE,
                col.names = !file.exists(tier_history_path), append = file.exists(tier_history_path))
    
    # Warn (not just log to CSV) if this run's history shows Tier 1 has been
    # the ONLY tier ever chosen across a meaningful number of past runs - that
    # is the concrete, checkable version of "Tier 2/3 exist but may never
    # actually be triggering for this station network."
    hist_check <- tryCatch(read.csv(tier_history_path, stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(hist_check) && nrow(hist_check) >= 5) {
      frac_tier1 <- mean(hist_check$chosen_tier == "tier1", na.rm = TRUE)
      if (frac_tier1 == 1) {
        cat(sprintf(paste0(
          "[NOTE] Across all %d logged run(s) in %s, the chosen bias-correction tier has ALWAYS been ",
          "'tier1' (basin-wide). If BIAS_TIER_MODE = 'auto' is set expecting Tier 2/3 to sometimes win, ",
          "this basin's station network may not be clearing the thresholds needed to even ATTEMPT them ",
          "(LOSO_CV_MIN_STATIONS = %d distinct stations; MIN_STATION_MONTHS_TIER2 = %d station-months per ",
          "region x month; MIN_STATION_MONTHS_TIER3 = %d station-months pooled). Check the ",
          "'tier2_fit_available'/'tier3_fit_available' columns in %s to see whether Tier 2/3 are even ",
          "being fit, versus fit-but-losing on LOSO-CV RMSE.\n"),
          nrow(hist_check), tier_history_path, LOSO_CV_MIN_STATIONS, MIN_STATION_MONTHS_TIER2,
          MIN_STATION_MONTHS_TIER3, tier_history_path))
      }
    }
    
    ## ---- Validation protocol (framework doc, Section 5): random k-fold,
    ## spatial-block (LOSO, already computed above when tier_mode=="auto"),
    ## and temporal-transfer CV - bias/RMSE/R^2 per tier x scheme. ----
    bc_stage("Running full validation protocol (random k-fold + temporal-transfer CV)...")
    validation_tab <- run_bias_correction_validation(
      station_obs_raw, swe[[1]], swe_matrix_uncorrected = swe_matrix_pre_bc,
      dates_vec = dates, region_id = region_id_bc, bias_predictors = bias_predictors,
      clip_range = BIAS_CORRECT_CLIP, min_n_tier2 = MIN_STATION_MONTHS_TIER2,
      min_n_tier3 = MIN_STATION_MONTHS_TIER3, loso_table = loso$table)
    if (!is.null(validation_tab) && nrow(validation_tab) > 0) {
      write.csv(validation_tab, file.path(bc_diag_dir, "bias_correction_validation.csv"), row.names = FALSE)
      cat(sprintf("Bias-correction validation written to: %s\n",
                  file.path(bc_diag_dir, "bias_correction_validation.csv")))
    }
    
    bc_stage("Basin-wide bias correction complete.")
    cat(sprintf("===== END BASIN-WIDE BIAS CORRECTION (total time: %s) =====\n",
                format_hms(as.numeric(difftime(Sys.time(), bc_step_t0, units = "secs")))))
    
  } else {
    basin_bias_factor <- rep(1, 12)
    basin_bias_log <- data.frame()
    chosen_tier <- NA_character_
    if (isTRUE(BIAS_CORRECT)) {
      log_pipeline_fallback(
        sprintf(paste0(
          "BIAS_CORRECT = TRUE but no station file found at '%s' - proceeding ",
          "WITHOUT bias correction (factor = 1 for every month). ERA5-Land snow ",
          "variables carry known systematic biases (Kouki et al. 2023; Blau et al. ",
          "2024). Supply station SWE observations at STATION_OBS_PATH (columns: ",
          "station_id, lon, lat, date, swe_mm) to enable this correction, or set ",
          "REQUIRE_STATION_OBS <- TRUE to make this a hard stop instead."),
          STATION_OBS_PATH),
        hard_stop = isTRUE(REQUIRE_STATION_OBS)
      )
    } else {
      cat("BIAS_CORRECT = FALSE: skipping basin-wide bias correction.\n")
    }
  }
}

# ==============================================================================
#   SNOW-COVER-TIMING VALIDATION AGAINST OPENLANDMAP CLIMATOLOGY  [method 2]
# ==============================================================================
# Framework doc, Section 4: the bias correction above targets SWE MAGNITUDE;
# it says nothing about whether the corrected snow season's TIMING (onset,
# melt-out) agrees with an independent, long (2000-2020) satellite record -
# Bresson et al. (2026) show quantile-mapping-style corrections routinely get
# magnitude right while getting season length wrong. This closes the
# SCF_CLIMATOLOGY_DIR / SCF_CLIMATOLOGY_VALIDATE stub declared in the
# CONFIGURATION block: those files are written by Section 4 of
# Nechako_Snow_Data_download.R (OpenLandMap "Global MODIS-based snow cover
# monthly long-term", Hengl 2021, https://doi.org/10.5281/zenodo.5774953),
# one P05/P95 snow-cover-fraction raster PER CALENDAR MONTH, each already
# averaged across 2000-2020 and cropped to the Nechako extent.
validate_scf_timing <- function(swe_matrix, dates_vec, basin_pix_idx,
                                scf_dir, swe_threshold_mm, out_dir) {
  if (!dir.exists(scf_dir) || length(list.files(scf_dir, pattern = "\\.tif$")) == 0) {
    cat(sprintf("SCF_CLIMATOLOGY_VALIDATE = TRUE but no climatology rasters found in '%s' - skipping.\n", scf_dir))
    cat("  (Run Section 4 of Nechako_Snow_Data_download.R to populate it from the\n")
    cat("  OpenLandMap climatology - see framework doc, Section 4.)\n")
    return(invisible(NULL))
  }
  mon_lab <- tolower(format(as.Date(paste0("2020-", 1:12, "-01")), "%b"))
  mon <- as.integer(format(dates_vec, "%m"))
  
  # Basin-wide monthly snow-cover FREQUENCY from the corrected ERA5-Land SWE,
  # using the Bresson et al. (2026) SWE threshold for "snow covered".
  scf_era <- sapply(1:12, function(m) {
    cols <- which(mon == m)
    if (length(cols) == 0 || length(basin_pix_idx) == 0) return(NA_real_)
    mean(swe_matrix[basin_pix_idx, cols, drop = FALSE] >= swe_threshold_mm, na.rm = TRUE)
  })
  
  results <- data.frame()
  for (m in 1:12) {
    p05_f <- list.files(scf_dir, pattern = sprintf("p0?5.*%s", mon_lab[m]), full.names = TRUE, ignore.case = TRUE)
    p95_f <- list.files(scf_dir, pattern = sprintf("p95.*%s", mon_lab[m]), full.names = TRUE, ignore.case = TRUE)
    if (length(p05_f) == 0 || length(p95_f) == 0 || !is.finite(scf_era[m])) next
    p05_r <- try(rast(p05_f[1]), silent = TRUE)
    p95_r <- try(rast(p95_f[1]), silent = TRUE)
    if (inherits(p05_r, "try-error") || inherits(p95_r, "try-error")) next
    p05_v <- suppressWarnings(mean(values(p05_r), na.rm = TRUE))
    p95_v <- suppressWarnings(mean(values(p95_r), na.rm = TRUE))
    if (is.finite(p05_v) && p05_v > 1.5) p05_v <- p05_v / 100  # percent -> fraction
    if (is.finite(p95_v) && p95_v > 1.5) p95_v <- p95_v / 100
    in_band <- is.finite(p05_v) && is.finite(p95_v) && scf_era[m] >= p05_v && scf_era[m] <= p95_v
    results <- rbind(results, data.frame(month = m, month_label = month.abb[m],
                                         era5_scf_corrected = scf_era[m],
                                         climatology_p05 = p05_v, climatology_p95 = p95_v,
                                         within_p05_p95 = in_band))
  }
  if (nrow(results) == 0) {
    cat("SCF climatology timing check: no matching climatology rasters found for any month - skipping.\n")
    return(invisible(NULL))
  }
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  write.csv(results, file.path(out_dir, "scf_climatology_validation.csv"), row.names = FALSE)
  n_out <- sum(!results$within_p05_p95, na.rm = TRUE)
  cat(sprintf("SCF climatology timing check: %d/%d month(s) with corrected-ERA5-Land snow-cover frequency outside the 2000-2020 P05-P95 band (see scf_climatology_validation.csv).\n",
              n_out, nrow(results)))
  if (n_out > 0) {
    cat("  A systematic mismatch here means the magnitude correction above fixed the\n")
    cat("  mean/CDF but not the snow-season TIMING (Bresson et al. 2026) - consider a\n")
    cat("  phenological shift/stretch correction, or excluding the affected weeks from\n")
    cat("  Method-1 SCF-mask output until timing is reconciled (framework doc, Section 4).\n")
  }
  results
}

if (scf_method == "2" && isTRUE(SCF_CLIMATOLOGY_VALIDATE)) {
  cat("\n===== SNOW-COVER-TIMING VALIDATION (OpenLandMap climatology) =====\n")
  scf_timing_log <- validate_scf_timing(swe_matrix, dates, basin_rows_prescreen,
                                        SCF_CLIMATOLOGY_DIR, SNOW_COVER_SWE_THRESHOLD_MM,
                                        file.path(out_dir, "rfa_diagnostics"))
}
# ==============================================================================
#This part reads the basin-mean monthly CSVs produced by Nechako_Snow_Bridge.R 
# (snow_data/cross_validation/) and reports
# correlation + bias of the bias-corrected ERA5-Land basin-mean series
# against each independent product that's actually available, WITHOUT
# touching swe_matrix itself - this is diagnostic-only, same spirit as the
# existing OpenLandMap check.
## Uses the SAME basin_pix_idx / swe_matrix / dates already in scope at the
# point in H8SWEI_Snow.R - no new inputs needed.
#
# COVERAGE NOTE: as of ... only
# download/produce the calendar months listed in that script's SEASON_MONTHS
# config (default: every month except Aug/Sep) rather than all 12, since this
# cross-check is the only thing ... feeds - it does not touch
# SWEI/SSPI-1 itself. That's fine for this block (still comfortably above the
# 6-overlapping-month minimum below with the default 10/12 months), but if
# SEASON_MONTHS is narrowed further on the download side, some products may
# start reporting "fewer than 6 overlapping months" below more often - that's
# expected in that case, not a bug here.
# ==============================================================================

CROSS_VALIDATION_DIR     <- "snow_data/cross_validation"  # written by Part D of Nechako_Snow_Data_download.R
CROSS_VALIDATE_EXTRA_PRODUCTS <- TRUE

validate_against_cross_products <- function(swe_matrix, dates_vec, basin_pix_idx, cv_dir, out_dir) {
  if (!dir.exists(cv_dir)) {
    cat(sprintf("Cross-validation dir '%s' not found - skipping run Nechako_Snow_Bridge.R first(between the download script and H8)).\n", cv_dir))
    return(invisible(NULL))
  }
  csvs <- list.files(cv_dir, pattern = "_monthly_basinmean\\.csv$", full.names = TRUE)
  if (length(csvs) == 0) {
    cat(sprintf("No *_monthly_basinmean.csv files found in '%s' - skipping.\n", cv_dir))
    return(invisible(NULL))
  }
  
  # H8's own basin-mean monthly SWE, for comparison against every available series.
  if (length(basin_pix_idx) == 0) {
    cat("No basin pixels available - skipping cross-product validation.\n")
    return(invisible(NULL))
  }
  era5_basin_mean <- colMeans(swe_matrix[basin_pix_idx, , drop = FALSE], na.rm = TRUE)
  era5_df <- data.frame(date = as.Date(dates_vec), era5land_swe_mm = era5_basin_mean)
  
  # ----------------------------------------------------------------------
  # CMC (NSIDC-0447) MONTHLY-ONLY CADENCE
  # ----------------------------------------------------------------------
  # CMC has no daily product at all - it is documented (see the download
  # script's header and build_cmc_basinmean() in Nechako_Snow_Bridge.R) as
  # producing SWE estimates at MONTHLY resolution, for Oct-Jun only, each
  # water year. The merge below already only compares on a shared
  # calendar-month key (both sides truncated to year-month), so it never
  # fabricates a sub-monthly CMC value - but nothing previously recorded
  # that fact for a reader of the OUTPUT table, and nothing checked whether
  # CMC's own within-season monthly coverage (Oct-Jun) was actually complete
  # for the months it did claim to cover. Both are addressed below: any
  # product whose file matches CMC's naming convention is tagged with its
  # native cadence and season restriction in the results table, and its
  # Oct-Jun completeness is checked and logged (informational only - this
  # does not change what gets compared, since the year-month merge already
  # protects against inventing resolution CMC does not have).
  is_cmc_product <- function(varname, fpath) {
    grepl("^cmc", varname, ignore.case = TRUE) || grepl("cmc", basename(fpath), ignore.case = TRUE)
  }
  CMC_SEASON_MONTHS <- c(10, 11, 12, 1, 2, 3, 4, 5, 6)  # Oct-Jun, per NSIDC-0447 documentation
  
  results <- data.frame()
  for (f in csvs) {
    varname <- names(read.csv(f, nrows = 1))[2]
    cv <- tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)
    if (is.null(cv) || !("date" %in% names(cv))) next
    cv$date <- as.Date(cv$date)
    
    is_cmc <- is_cmc_product(varname, f)
    if (is_cmc) {
      # Informational completeness check only - CMC's own within-season
      # (Oct-Jun) monthly coverage across the years present in this file.
      # A gap here reflects a real hole in CMC's record for that water year,
      # not a bug in this script; it's logged so a thin comparison isn't
      # mistaken for a full one.
      cmc_years <- unique(as.integer(format(cv$date, "%Y")))
      cv_months <- as.integer(format(cv$date, "%m"))
      n_present <- sum(cv_months %in% CMC_SEASON_MONTHS)
      n_expected <- length(cmc_years) * length(CMC_SEASON_MONTHS)
      cat(sprintf(
        "  [%s] CMC (NSIDC-0447): monthly-only, Oct-Jun each water year. %d/%d Oct-Jun month-slots present across %d year(s) in this file (gaps reflect CMC's own record, not this comparison).\n",
        varname, n_present, n_expected, length(cmc_years)))
    }
    
    # Compare on a common calendar-month basis (cv series may be daily-derived
    # monthly means with day-of-month != 1, or SWE-vs-SCF units - this is a
    # RELATIVE-AGREEMENT check, not a unit-matched validation, except where
    # units already coincide, e.g. AU_DySno SWE mm vs ERA5-Land SWE mm).
    # For CMC specifically, this year-month merge is also what keeps the
    # comparison honestly monthly - it truncates BOTH series to whole
    # calendar months rather than ever upsampling CMC to a finer timescale.
    cv$yearmon <- format(cv$date, "%Y-%m")
    era5_df$yearmon <- format(era5_df$date, "%Y-%m")
    merged <- merge(era5_df[, c("yearmon", "era5land_swe_mm")], cv[, c("yearmon", varname)], by = "yearmon")
    if (nrow(merged) < 6) {
      cat(sprintf("  [%s] fewer than 6 overlapping months with ERA5-Land - skipping correlation.\n", varname))
      next
    }
    is_swe_like <- grepl("swe", varname, ignore.case = TRUE)
    r_val <- suppressWarnings(cor(merged$era5land_swe_mm, merged[[varname]], use = "complete.obs"))
    bias_val <- if (is_swe_like) {
      mean(merged$era5land_swe_mm - merged[[varname]], na.rm = TRUE)
    } else {
      NA_real_  # bias in mm doesn't make sense against a fraction/NDSI series
    }
    results <- rbind(results, data.frame(
      product = varname, n_months = nrow(merged), correlation_vs_era5land = r_val,
      mean_bias_mm_vs_era5land = bias_val,
      native_cadence = if (is_cmc) "monthly_only_Oct-Jun" else "see_source_documentation"
    ))
    cat(sprintf("  [%s] n=%d months, r=%.2f%s\n", varname, nrow(merged), r_val,
                if (is_swe_like) sprintf(", mean bias=%.1f mm (ERA5-Land minus %s)", bias_val, varname) else ""))
  }
  
  if (nrow(results) == 0) {
    cat("No cross-validation product had enough overlapping months - nothing written.\n")
    return(invisible(NULL))
  }
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  write.csv(results, file.path(out_dir, "cross_validation_extra_products.csv"), row.names = FALSE)
  cat(sprintf("Cross-validation vs. AU_DySno/MODIS-VIIRS/HLS/Crocus-ERA5 written to: %s\n",
              file.path(out_dir, "cross_validation_extra_products.csv")))
  results
}

if (scf_method == "2" && isTRUE(CROSS_VALIDATE_EXTRA_PRODUCTS)) {
  cat("\n===== CROSS-VALIDATION vs. AU_DySno / MODIS-VIIRS / HLS / Crocus-ERA5 =====\n")
  extra_cv_log <- validate_against_cross_products(
    swe_matrix, dates, basin_rows_prescreen,
    CROSS_VALIDATION_DIR, file.path(out_dir, "rfa_diagnostics")
  )
}

robustness_diag_dir <- file.path(out_dir, "robustness_diagnostics")
if (!dir.exists(robustness_diag_dir)) dir.create(robustness_diag_dir, recursive = TRUE)

# ==============================================================================
# SHARED HELPER: monthly basin-mean series from a *_monthly_basinmean.csv
# (same file convention validate_against_cross_products() already reads) - kept
# separate/local rather than editing that existing, tested function.
# ==============================================================================
.load_cv_basinmean <- function(cv_dir, filename, varname) {
  f <- file.path(cv_dir, filename)
  if (!file.exists(f)) return(NULL)
  df <- tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(df) || !all(c("date", varname) %in% names(df))) return(NULL)
  df$date <- as.Date(df$date)
  df[, c("date", varname)]
}

# Station observations, collapsed to a monthly basin-mean series (simple mean
# across all stations reporting that month - matches the "pooled across ALL
# stations" philosophy already used in fit_basin_bias_correction()).
.station_basin_mean <- function(station_obs) {
  if (is.null(station_obs) || nrow(station_obs) == 0) return(NULL)
  d <- as.Date(station_obs$date)
  ym <- format(d, "%Y-%m-01")
  agg <- stats::aggregate(swe_mm ~ ym, data = data.frame(ym = ym, swe_mm = station_obs$swe_mm),
                          FUN = function(x) mean(x, na.rm = TRUE))
  data.frame(date = as.Date(agg$ym), station_swe_mm = agg$swe_mm)
}

# Raw (pre-correction) ERA5-Land basin-mean series, with a graceful fallback
# to the current swe_matrix if bias correction never ran this session (in
# which case swe_matrix IS the raw series, so it's still a fair comparison).
.era5_raw_basin_mean <- function() {
  src <- if (exists("swe_matrix_pre_bc", inherits = FALSE)) swe_matrix_pre_bc else swe_matrix
  data.frame(date = as.Date(dates),
             era5land_swe_mm = colMeans(src[basin_rows_prescreen, , drop = FALSE], na.rm = TRUE))
}


# ==============================================================================
# 1. TRIPLE COLLOCATION - ERA5-Land (raw) x CMC x station obs
# ==============================================================================
# Classic covariance-based TC error-variance estimator (Stoffelen 1998;
# McColl et al. 2014). Assumes a linear error model per source
# (x = alpha_x + beta_x*truth + eps_x), errors mutually uncorrelated and
# uncorrelated with truth - it does NOT require treating any of the three as
# ground truth, which is the whole point: it lets you check whether the
# "station obs = truth" assumption baked into fit_basin_bias_correction() /
# fit_regionalized_qm() / fit_mlr_correction() is actually reasonable.
#
#   sigma2_ex = var(x) - cov(x,y)*cov(x,z)/cov(y,z)
#   sigma2_ey = var(y) - cov(x,y)*cov(y,z)/cov(x,z)
#   sigma2_ez = var(z) - cov(x,z)*cov(y,z)/cov(x,y)
#
# Optionally also builds an extended-TC rescaled, minimum-variance merged
# series purely as a DIAGNOSTIC comparison against the Tier 1-3 corrected
# series - it is NOT written back into swe_matrix. Wiring it in for real
# would be a separate, deliberate follow-up once you've reviewed a few runs
# of this diagnostic and trust the triplet's overlap size.
triple_collocation_swe <- function(era5_df, cmc_df, station_df, min_overlap_months = 24) {
  if (is.null(cmc_df) || is.null(station_df)) {
    cat("[TC] Need CMC and station basin-mean series both present - skipping Triple Collocation.\n")
    return(invisible(NULL))
  }
  m <- merge(era5_df, cmc_df, by = "date")
  m <- merge(m, station_df, by = "date")
  m <- m[stats::complete.cases(m), ]
  if (nrow(m) < min_overlap_months) {
    cat(sprintf("[TC] Only %d overlapping month(s) across ERA5-Land/CMC/station (need >= %d) - skipping.\n",
                nrow(m), min_overlap_months))
    return(invisible(NULL))
  }
  
  x <- m$era5land_swe_mm; y <- m$cmc_swe_mm; z <- m$station_swe_mm
  cxy <- stats::cov(x, y); cxz <- stats::cov(x, z); cyz <- stats::cov(y, z)
  
  sigma2_ex <- stats::var(x) - (cxy * cxz) / cyz
  sigma2_ey <- stats::var(y) - (cxy * cyz) / cxz
  sigma2_ez <- stats::var(z) - (cxz * cyz) / cxy
  
  res <- data.frame(
    source = c("era5land_raw", "cmc", "station_obs"),
    error_variance_mm2 = c(sigma2_ex, sigma2_ey, sigma2_ez),
    error_sd_mm = sqrt(pmax(c(sigma2_ex, sigma2_ey, sigma2_ez), 0)),
    n_overlap_months = nrow(m)
  )
  
  # Negative error variance is a known TC failure mode (poor overlap, strongly
  # correlated errors between two of the three sources, or a near-singular
  # denominator) - flag rather than silently clip.
  if (any(res$error_variance_mm2 < 0)) {
    cat("[TC] WARNING: at least one source produced a negative error variance estimate -\n")
    cat("     treat this run's TC output as unreliable (too little/too noisy overlap data)\n")
    cat("     rather than as evidence that source has zero error.\n")
  }
  
  write.csv(res, file.path(robustness_diag_dir, "triple_collocation_error_variances.csv"), row.names = FALSE)
  cat("\n[TC] Triple Collocation error-variance estimates (mm SWE, RMS):\n")
  print(res)
  
  # Sanity check against the "station = truth" assumption used everywhere
  # else in the bias-correction block: if the station error SD comes out
  # comparable to or larger than ERA5-Land's own, the tiered correction may
  # be pulling the grid toward a noisier reference than intended.
  tc_valid <- all(is.finite(res$error_variance_mm2)) && all(res$error_variance_mm2 > 0) &&
    all(is.finite(c(cxy, cxz, cyz))) && all(abs(c(cxy, cxz, cyz)) > .Machine$double.eps)
  if (tc_valid && res$error_sd_mm[res$source == "station_obs"] >=
      0.75 * res$error_sd_mm[res$source == "era5land_raw"]) {
    log_pipeline_fallback(sprintf(
      "Triple Collocation estimates station-obs error SD (%.1f mm) at >=75%% of ERA5-Land's own raw error SD (%.1f mm) over %d overlap months. The Tier 1-3 bias correction treats station data as ground truth - review triple_collocation_error_variances.csv before trusting the corrected grid.",
      res$error_sd_mm[res$source == "station_obs"], res$error_sd_mm[res$source == "era5land_raw"], nrow(m)))
  }
  
  # --- Optional: extended-TC minimum-variance merged series (diagnostic only) ---
  beta_y <- cxz / cyz   # rescales CMC onto ERA5-Land's units
  beta_z <- cxy / cyz   # rescales station onto ERA5-Land's units
  y_resc <- mean(x) + beta_y * (y - mean(y))
  z_resc <- mean(x) + beta_z * (z - mean(z))
  w <- if (tc_valid) 1 / c(sigma2_ex, beta_y^2 * sigma2_ey, beta_z^2 * sigma2_ez) else rep(NA_real_, 3)
  if (tc_valid && all(is.finite(w)) && all(w > 0)) {
    w <- w / sum(w)
    merged <- w[1] * x + w[2] * y_resc + w[3] * z_resc
    merge_diag <- data.frame(date = m$date, era5land_raw = x, cmc = y, station = z,
                             tc_merged = merged, tier_corrected =
                               if (exists("swe_matrix", inherits = FALSE))
                                 colMeans(swe_matrix[basin_rows_prescreen, match(m$date, as.Date(dates)), drop = FALSE], na.rm = TRUE)
                             else NA_real_)
    write.csv(merge_diag, file.path(robustness_diag_dir, "triple_collocation_merged_vs_tiered.csv"), row.names = FALSE)
    cat(sprintf("[TC] Diagnostic-only merged series written (weights: ERA5=%.2f, CMC=%.2f, station=%.2f) -\n     compare tc_merged vs tier_corrected columns before deciding whether to adopt TC weighting for real.\n",
                w[1], w[2], w[3]))
  }
  
  invisible(res)
}

era5_raw_bm    <- .era5_raw_basin_mean()
cmc_bm         <- .load_cv_basinmean(CROSS_VALIDATION_DIR, "cmc_swe_monthly_basinmean.csv", "cmc_swe_mm")
station_bm     <- if (exists("station_obs_raw", inherits = FALSE)) .station_basin_mean(station_obs_raw) else NULL

cat("\n===== TRIPLE COLLOCATION (ERA5-Land x CMC x station obs) =====\n")
tc_result <- triple_collocation_swe(era5_raw_bm, cmc_bm, station_bm)


# ==============================================================================
# 2. ERA5-LAND HOMOGENEITY CHECK (double-mass curve + Pettitt test vs. CMC)
# ==============================================================================
# Tests the RAW ERA5-Land record itself for an internal structural break,
# using CMC as the independent reference over their shared window. This is
# necessarily limited to whatever period CMC actually covers - it cannot see
# a break that predates CMC's own record. Say so explicitly in the output
# rather than implying full-record coverage.
check_era5_homogeneity <- function(era5_df, cmc_df, out_dir) {
  if (is.null(cmc_df)) {
    cat("[Homogeneity] No CMC basin-mean series available - skipping.\n")
    return(invisible(NULL))
  }
  m <- merge(era5_df, cmc_df, by = "date")
  m <- m[stats::complete.cases(m), ]
  m <- m[order(m$date), ]
  if (nrow(m) < 60) {
    cat(sprintf("[Homogeneity] Only %d overlapping month(s) with CMC (need >= 60, i.e. ~5 years) - skipping.\n", nrow(m)))
    return(invisible(NULL))
  }
  
  cat(sprintf("[Homogeneity] NOTE: this test can only detect a break WITHIN the ERA5-Land/CMC overlap window (%s to %s). It cannot see a break earlier in the 1950-present ERA5-Land record, since no independent reference that old is available in this pipeline.\n",
              min(m$date), max(m$date)))
  
  # Double-mass curve: cumulative ERA5-Land vs cumulative CMC. A stable
  # relationship is a straight line; a slope change marks a candidate break.
  m$cum_era5 <- cumsum(m$era5land_swe_mm)
  m$cum_cmc  <- cumsum(m$cmc_swe_mm)
  write.csv(m[, c("date", "cum_era5", "cum_cmc")],
            file.path(out_dir, "era5land_double_mass_curve.csv"), row.names = FALSE)
  
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    p <- ggplot2::ggplot(m, ggplot2::aes(x = cum_cmc, y = cum_era5)) +
      ggplot2::geom_point(size = 0.8) +
      ggplot2::geom_smooth(method = "lm", se = FALSE, linewidth = 0.5) +
      ggplot2::labs(title = "ERA5-Land vs CMC double-mass curve",
                    x = "Cumulative CMC SWE (mm)", y = "Cumulative ERA5-Land SWE (mm)")
    ggplot2::ggsave(file.path(out_dir, "era5land_double_mass_curve.png"), p, width = 6, height = 5)
  }
  
  # Pettitt's test (trend::pettitt.test) on the ERA5-CMC ratio series: a
  # single most-likely change-point plus a significance test, simpler and
  # better-supported than a hand-rolled SNHT critical-value table.
  if (!requireNamespace("trend", quietly = TRUE)) {
    cat("[Homogeneity] Package 'trend' not installed (install.packages(\"trend\")) - skipping Pettitt test, double-mass CSV/plot written above only.\n")
    return(invisible(m))
  }
  ratio <- m$era5land_swe_mm / m$cmc_swe_mm
  ratio[!is.finite(ratio)] <- NA
  ratio <- ratio[!is.na(ratio)]
  pt <- tryCatch(trend::pettitt.test(ratio), error = function(e) NULL)
  if (!is.null(pt)) {
    pt <- trend::pettitt.test(ratio)
    
    break_idx <- as.integer(unname(pt$estimate))[1L]
    
    stopifnot(
      length(break_idx) == 1L,
      break_idx >= 1L,
      break_idx <= nrow(m)
    )
    
    pettitt_tbl <- data.frame(
      candidate_break_date = as.Date(m$date[break_idx]),
      statistic            = as.numeric(unname(pt$statistic)),
      p_value              = as.numeric(unname(pt$p.value)),
      stringsAsFactors     = FALSE,
      row.names            = NULL
    )
    
    write.csv(
      pettitt_tbl,
      file.path(out_dir, "era5land_homogeneity_pettitt_test.csv"),
      row.names = FALSE
    )
    break_date <-as.Date( m$date[break_idx])
    cat(sprintf("[Homogeneity] Pettitt test on ERA5-Land/CMC ratio: candidate break at %s (K=%.1f, p=%.3f).\n",
                break_date, pt$statistic, pt$p.value))
    if (pt$p.value < 0.05) {
      log_pipeline_fallback(sprintf(
        "Pettitt homogeneity test flags a likely structural break in ERA5-Land relative to CMC near %s (p=%.3f). REFERENCE_START/REFERENCE_END (%s to %s) and any long-term trend statements should be checked against this before being reported as a physical signal.",
        break_date, pt$p.value, REFERENCE_START, REFERENCE_END))
    }
    pettitt_tbl <- data.frame(
      candidate_break_date = break_date,
      statistic = as.numeric(unname(pt$statistic)),
      p_value = as.numeric(unname(pt$p.value)),
      stringsAsFactors = FALSE
    )
    rownames(pettitt_tbl) <- NULL
    write.csv(
      pettitt_tbl,
      file.path(out_dir,
                "era5land_homogeneity_pettitt_test.csv"),
      row.names = FALSE
    )
  }
  invisible(m)
}

cat("\n===== ERA5-LAND HOMOGENEITY CHECK (vs. CMC) =====\n")
homog_result <- check_era5_homogeneity(era5_raw_bm, cmc_bm, robustness_diag_dir)


# ==============================================================================
# 3. TREND-PRESERVATION CHECK (raw vs. bias-corrected ERA5-Land)
# ==============================================================================
# Quantile mapping (Tier 2) is not trend-preserving by construction - this
# checks whether the tier actually applied this run measurably changed the
# long-term Sen's-slope trend of basin-mean SWE relative to the raw record.
check_trend_preservation <- function(out_dir) {
  if (!exists("swe_matrix_pre_bc", inherits = TRUE)) {
    cat("[Trend] No pre-correction snapshot available this run (bias correction did not run) - skipping.\n")
    return(invisible(NULL))
  }
  if (!requireNamespace("trend", quietly = TRUE)) {
    cat("[Trend] Package 'trend' not installed (install.packages(\"trend\")) - skipping.\n")
    return(invisible(NULL))
  }
  raw_bm  <- colMeans(swe_matrix_pre_bc[basin_rows_prescreen, , drop = FALSE], na.rm = TRUE)
  corr_bm <- colMeans(swe_matrix[basin_rows_prescreen, , drop = FALSE], na.rm = TRUE)
  
  # Annual peak (max monthly) basin SWE is the standard snow-drought trend
  # metric - less sensitive to within-year noise than the monthly series.
  yr <- as.integer(format(dates, "%Y"))
  ann_raw  <- tapply(raw_bm, yr, max, na.rm = TRUE)
  ann_corr <- tapply(corr_bm, yr, max, na.rm = TRUE)
  ann_raw[!is.finite(ann_raw)] <- NA; ann_corr[!is.finite(ann_corr)] <- NA
  
  ts_raw  <- ts(ann_raw,  start = min(yr))
  ts_corr <- ts(ann_corr, start = min(yr))
  
  sen_raw  <- tryCatch(trend::sens.slope(ts_raw[!is.na(ts_raw)]),  error = function(e) NULL)
  sen_corr <- tryCatch(trend::sens.slope(ts_corr[!is.na(ts_corr)]), error = function(e) NULL)
  mk_raw   <- tryCatch(trend::mk.test(ts_raw[!is.na(ts_raw)]),   error = function(e) NULL)
  mk_corr  <- tryCatch(trend::mk.test(ts_corr[!is.na(ts_corr)]), error = function(e) NULL)
  
  if (is.null(sen_raw) || is.null(sen_corr)) {
    cat("[Trend] Sen's slope could not be computed (too few annual values) - skipping.\n")
    return(invisible(NULL))
  }
  
  res <- data.frame(
    series = c("raw_era5land", "bias_corrected"),
    tier_applied = c(NA, if (exists("chosen_tier", inherits = TRUE)) chosen_tier else NA),
    sen_slope_mm_per_year = c(sen_raw$estimates, sen_corr$estimates),
    mk_p_value = c(mk_raw$p.value, mk_corr$p.value)
  )
  write.csv(res, file.path(out_dir, "trend_preservation_check.csv"), row.names = FALSE)
  cat("\n[Trend] Sen's slope, annual peak basin SWE, raw vs. corrected:\n")
  print(res)
  
  sign_flip <- sign(sen_raw$estimates) != sign(sen_corr$estimates) &&
    mk_raw$p.value < 0.1 && mk_corr$p.value < 0.1
  pct_change <- if (sen_raw$estimates != 0) abs((sen_corr$estimates - sen_raw$estimates) / sen_raw$estimates) * 100 else NA
  if (isTRUE(sign_flip) || (is.finite(pct_change) && pct_change > 50)) {
    log_pipeline_fallback(sprintf(
      "Bias correction (tier: %s) changed the Sen's-slope trend in annual peak basin SWE from %.2f to %.2f mm/yr (%s). If this run used Tier 2 quantile mapping, this is the known non-trend-preserving failure mode - see trend_preservation_check.csv.",
      if (exists("chosen_tier", inherits = TRUE)) chosen_tier else "unknown",
      sen_raw$estimates, sen_corr$estimates,
      if (isTRUE(sign_flip)) "SIGN FLIPPED" else sprintf("%.0f%% change", pct_change)))
  }
  invisible(res)
}

cat("\n===== TREND-PRESERVATION CHECK =====\n")
trend_result <- check_trend_preservation(robustness_diag_dir)


# ==============================================================================
# 4. EXTREME / RETURN-PERIOD CHECK (raw vs. bias-corrected, low-snow years)
# ==============================================================================
# Fits a GEV to CORE-WINTER (Jan-Mar) seasonal-minimum basin SWE rather than
# calendar-year minima. Calendar-year minima are structurally ~0 in a seasonal
# snow basin and produce a degenerate lower tail, negative GEV return levels,
# and meaningless percentage differences. The check is skipped when the winter
# minima are themselves near-degenerate.
check_extreme_return_periods <- function(out_dir, return_periods = c(20, 50)) {
  if (!exists("swe_matrix_pre_bc", inherits = TRUE)) {
    cat("[Extremes] No pre-correction snapshot available this run - skipping.\n")
    return(invisible(NULL))
  }
  raw_bm  <- colMeans(swe_matrix_pre_bc[basin_rows_prescreen, , drop = FALSE], na.rm = TRUE)
  corr_bm <- colMeans(swe_matrix[basin_rows_prescreen, , drop = FALSE], na.rm = TRUE)
  yr <- as.integer(format(dates, "%Y"))
  mo <- as.integer(format(dates, "%m"))
  winter <- mo %in% c(1L, 2L, 3L)
  winter_min_raw  <- tapply(raw_bm[winter],  yr[winter], min, na.rm = TRUE)
  winter_min_corr <- tapply(corr_bm[winter], yr[winter], min, na.rm = TRUE)
  winter_min_raw  <- winter_min_raw[is.finite(winter_min_raw)]
  winter_min_corr <- winter_min_corr[is.finite(winter_min_corr)]
  
  if (length(winter_min_raw) < 20 || length(winter_min_corr) < 20) {
    cat(sprintf("[Extremes] Only %d/%d Jan-Mar winter-minimum value(s) (need >= 20) - skipping.\n",
                length(winter_min_raw), length(winter_min_corr)))
    return(invisible(NULL))
  }
  raw_scale <- stats::IQR(winter_min_raw, na.rm = TRUE)
  corr_scale <- stats::IQR(winter_min_corr, na.rm = TRUE)
  if (!is.finite(raw_scale) || !is.finite(corr_scale) || raw_scale <= 1e-6 || corr_scale <= 1e-6) {
    cat("[Extremes] Jan-Mar winter minima are near-degenerate (IQR ~ 0); skipping GEV check.\n")
    return(invisible(NULL))
  }
  
  fit_return_levels <- function(x, periods) {
    lm3 <- lmomco::lmoms(x)
    par_gev <- tryCatch(lmomco::pargev(lm3), error = function(e) NULL)
    if (is.null(par_gev)) return(rep(NA_real_, length(periods)))
    sapply(periods, function(rp) lmomco::quagev(1 / rp, par_gev))
  }
  
  
  rl_raw  <- fit_return_levels(winter_min_raw,  return_periods)
  rl_corr <- fit_return_levels(winter_min_corr, return_periods)
  
  res <- data.frame(
    return_period_years = return_periods,
    raw_return_level_mm = rl_raw,
    corrected_return_level_mm = rl_corr,
    pct_difference = ifelse(abs(rl_raw) > 1e-6, 100 * (rl_corr - rl_raw) / rl_raw, NA_real_)
  )
  write.csv(res, file.path(out_dir, "extreme_return_period_check.csv"), row.names = FALSE)
  cat("\n[Extremes] GEV-implied low-snow return levels, Jan-Mar winter-minimum basin SWE, raw vs. corrected:\n")
  print(res)
  
  if (any(is.finite(res$pct_difference) & abs(res$pct_difference) > 25)) {
    log_pipeline_fallback(sprintf(
      "Bias correction (tier: %s) changes the GEV-implied low-snow return level by >25%% at at least one return period (see extreme_return_period_check.csv) - drought-frequency statements in the summary output (Exceptional/Extreme deficit %%) should be read alongside this, since they inherit whichever series (corrected) was used to compute them.",
      if (exists("chosen_tier", inherits = TRUE)) chosen_tier else "unknown"))
  }
  invisible(res)
}

cat("\n===== EXTREME / RETURN-PERIOD CHECK =====\n")
extreme_result <- check_extreme_return_periods(robustness_diag_dir)


# ==============================================================================
# 5. CHEAP BASIN-MEAN UNCERTAINTY BAND ON SSPI (reuses monthly_sspi() as-is)
# ==============================================================================
# NOT a per-pixel uncertainty map (that would mean re-running the parallel
# per-pixel loop N times - a separate, bigger change). This perturbs the
# basin-mean SWE series by +/- the tier's own LOSO-CV RMSE (already computed,
# already written to bias_correction_validation.csv) and recomputes SSPI on
# just that one vector, three times, to get a defensible basin-mean band.
build_basin_sspi_uncertainty_band <- function(out_dir) {
  if (!exists("validation_tab", inherits = TRUE) || is.null(validation_tab) || nrow(validation_tab) == 0) {
    cat("[Uncertainty] No bias-correction validation table available this run - skipping.\n")
    return(invisible(NULL))
  }
  if (!exists("chosen_tier", inherits = TRUE)) {
    cat("[Uncertainty] No chosen_tier in scope - skipping.\n")
    return(invisible(NULL))
  }
  rmse_row <- validation_tab[validation_tab$tier == chosen_tier & validation_tab$scheme == "spatial_block_loso", ]
  if (nrow(rmse_row) == 0) rmse_row <- validation_tab[validation_tab$tier == chosen_tier, ]
  if (nrow(rmse_row) == 0 || !is.finite(rmse_row$rmse[1])) {
    cat("[Uncertainty] No usable RMSE found for the applied tier - skipping.\n")
    return(invisible(NULL))
  }
  rmse_mm <- rmse_row$rmse[1]
  
  corr_bm <- colMeans(swe_matrix[basin_rows_prescreen, , drop = FALSE], na.rm = TRUE)
  central <- monthly_sspi(corr_bm, dates, ref_idx)
  plus    <- monthly_sspi(pmax(corr_bm + rmse_mm, 0), dates, ref_idx)
  minus   <- monthly_sspi(pmax(corr_bm - rmse_mm, 0), dates, ref_idx)
  
  res <- data.frame(date = dates, sspi = central$sspi,
                    sspi_lower = pmin(plus$sspi, minus$sspi, na.rm = FALSE),
                    sspi_upper = pmax(plus$sspi, minus$sspi, na.rm = FALSE),
                    rmse_mm_used = rmse_mm, tier = chosen_tier)
  write.csv(res, file.path(out_dir, "basin_mean_sspi_uncertainty_band.csv"), row.names = FALSE)
  cat(sprintf("[Uncertainty] Basin-mean SSPI +/- band written using spatial-block-LOSO RMSE = %.1f mm (tier %s).\n",
              rmse_mm, chosen_tier))
  invisible(res)
}

cat("\n===== BASIN-MEAN SSPI UNCERTAINTY BAND =====\n")
uncertainty_result <- build_basin_sspi_uncertainty_band(robustness_diag_dir)

cat("\n[OK] Robustness diagnostics complete - see: ", robustness_diag_dir, "\n", sep = "")
# ==============================================================================
#   REGIONAL FREQUENCY ANALYSIS (RFA) CONSTRUCTION  [SSPI-1 / method 2 only]
# ==============================================================================
# Builds everything monthly_sspi_regional() needs: homogeneous regions, per-
# pixel index values, discordancy screening, heterogeneity diagnostics, and
# pooled regional gamma fits with goodness-of-fit screening. All of this runs
# ONCE (not per-pixel-in-parallel) and the resulting small objects are
# exported to the cluster workers. Station bias correction has already been
# applied basin-wide, upstream, directly to swe_matrix (see above) - so this
# block operates on already-corrected SWE.
if (scf_method == "2" && REGIONALIZE_SSPI) {
  
  cat("\n============================================================\n")
  cat("REGIONAL FREQUENCY ANALYSIS (RFA) CONSTRUCTION FOR SSPI-1\n")
  cat("============================================================\n")
  
  ## ---- A. Regionalization (elevation x aspect, or spatial k-means fallback) ----
  regions_built <- build_regions(swe[[1]], basin, DEM_PATH,
                                 n_elev_bands = N_ELEV_BANDS, use_aspect = USE_ASPECT,
                                 n_aspect_classes = N_ASPECT_CLASSES,
                                 n_regions_fallback = N_REGIONS_FALLBACK,
                                 min_region_pixels = MIN_REGION_PIXELS)
  region_id <- regions_built$region_id
  region_id[invalid_pix_sspi] <- NA_integer_   # zero-inflation exclusions stay excluded
  n_regions <- length(unique(region_id[!is.na(region_id)]))
  
  ## ---- B. Per-pixel, per-calendar-month index values (Eq. 1-2 analogue) ----
  cat("\n===== COMPUTING PER-PIXEL INDEX VALUES =====\n")
  index_values <- compute_index_values(swe_matrix, dates, ref_idx, stat = INDEX_VALUE_STAT)
  cat(sprintf("[OK] Index values (%s of reference-period monthly SWE): %.1f%% pixel-months valid\n",
              INDEX_VALUE_STAT, 100 * mean(is.finite(index_values))))
  
  ## ---- C. L-moment ratios per pixel-month, on normalized reference samples ----
  cat("\n===== COMPUTING L-MOMENT RATIOS PER PIXEL-MONTH =====\n")
  mon_ref <- as.integer(format(dates[ref_idx], "%m"))
  lmom_array <- array(NA_real_, dim = c(n_pixels, 12, 3),
                      dimnames = list(NULL, NULL, c("t2", "t3", "t4")))
  for (m in 1:12) {
    cols <- ref_idx[mon_ref == m]
    if (length(cols) == 0) next
    sub  <- swe_matrix[, cols, drop = FALSE]
    idxv <- index_values[, m]
    elig <- which(!is.na(region_id) & is.finite(idxv) & idxv > 0)
    for (i in elig) {
      x_star <- sub[i, ]
      x_star <- x_star[is.finite(x_star)] / idxv[i]
      lr <- lmom_ratios_site(x_star)
      lmom_array[i, m, ] <- lr[c("t2", "t3", "t4")]
    }
  }
  cat("[OK] L-moment ratios computed for all eligible pixel-months.\n")
  
  ## ---- D. Discordancy screening (Hosking & Wallis Di, Eq. 3-5) ----
  cat("\n===== DISCORDANCY SCREENING (Hosking & Wallis D_i) =====\n")
  discordant_mat <- matrix(FALSE, n_pixels, 12)
  discordancy_log <- data.frame()
  if (DISCORDANCY_SCREENING) {
    for (r in seq_len(n_regions)) {
      pix_r <- which(region_id == r)
      if (length(pix_r) < 4) next  # Di non-informative for N<=4 (S singular at N=3)
      for (m in 1:12) {
        lm_mat <- cbind(t2 = lmom_array[pix_r, m, "t2"],
                        t3 = lmom_array[pix_r, m, "t3"],
                        t4 = lmom_array[pix_r, m, "t4"])
        dt <- discordancy_test(lm_mat)
        if (dt$N >= 4 && any(dt$discordant)) {
          flagged <- pix_r[which(dt$discordant)]
          discordant_mat[flagged, m] <- TRUE
          discordancy_log <- rbind(discordancy_log,
                                   data.frame(region = r, month = m, pixel = flagged,
                                              Di = dt$Di[dt$discordant], Dcr = dt$Dcr))
        }
      }
    }
  }
  cat(sprintf("[OK] Discordancy screening: %d pixel-months flagged discordant (excluded from pooling).\n",
              sum(discordant_mat)))
  
  ## ---- E. Heterogeneity test (Hosking & Wallis H0, Monte Carlo) ----
  cat("\n===== HETEROGENEITY TEST (Hosking & Wallis H0) =====\n")
  
  # Each region/month H0 test is completely independent of every other one
  # (own pooled sample, own kappa/gamma fit, own Nsim simulations), so build
  # the task list up front and farm it out across cores rather than looping
  # serially. This is the dominant cost in the RFA construction block.
  het_tasks <- list()
  for (r in seq_len(n_regions)) {
    pix_r <- which(region_id == r)
    if (length(pix_r) < 2) next
    for (m in 1:12) {
      keep <- pix_r[!discordant_mat[pix_r, m]]
      if (length(keep) < 2) next
      het_tasks[[length(het_tasks) + 1]] <- list(region = r, month = m, keep = keep)
    }
  }
  cat(sprintf("[OK] %d region-month heterogeneity tests queued (Nsim=%d each)\n",
              length(het_tasks), HETEROGENEITY_NSIM))
  
  if (length(het_tasks) > 0) {
    n_tasks <- length(het_tasks)
    
    # Worker defined at top level so it serializes against each cluster node's
    # own .GlobalEnv (same mechanism clusterExport() relies on) AND, when run
    # serially below, resolves its free variables (swe_matrix, index_values,
    # ref_idx, mon_ref, HETEROGENEITY_NSIM) directly via normal lexical
    # scoping to this script's actual .GlobalEnv - no export/mirroring step
    # needed either way, since this whole RFA construction block already runs
    # at top level (unlike loso_cv_bias_tiers(), which is a function and
    # needed to mirror values into .GlobalEnv explicitly for its own serial
    # fallback).
    .het_worker <- function(task) {
      cols <- ref_idx[mon_ref == task$month]
      site_samples <- lapply(task$keep, function(i) {
        x <- swe_matrix[i, cols] / index_values[i, task$month]
        x[is.finite(x) & x > 0]
      })
      ht <- tryCatch(
        heterogeneity_test(site_samples, Nsim = HETEROGENEITY_NSIM),
        error = function(e) list(H0 = NA_real_, V0 = NA_real_, status = "error")
      )
      data.frame(region = task$region, month = task$month, H0 = ht$H0,
                 H0_mc_se = if (is.null(ht$H0_mc_se)) NA_real_ else ht$H0_mc_se,
                 H0_upper = if (is.null(ht$H0_upper)) NA_real_ else ht$H0_upper,
                 status = ht$status, n_sites = length(task$keep))
    }
    
    # ---- Serial execution path (used directly for small task counts, and
    # as the automatic fallback if the parallel cluster can't start) ----
    # Reproducibility note: the parallel path uses clusterSetRNGStream()
    # (L'Ecuyer-CMRG substreams, one per worker) so every worker's Monte
    # Carlo draws are independent and reproducible together. A single-
    # process set.seed() below is the serial equivalent for reproducibility
    # of THIS run, but draws from a different RNG stream than the parallel
    # path would - so a serial run and a parallel run of the same data will
    # each be internally reproducible but won't produce bit-identical H0
    # values against each other. Both are valid Monte Carlo estimates at the configured Nsim
    # estimates of the same H0 statistic; only the specific random draws
    # differ, not the methodology.
    run_het_serial <- function() {
      cat("[OK] Running heterogeneity tests SEQUENTIALLY (identical .het_worker()\n")
      cat("     computation, just one region-month at a time).\n")
      set.seed(BASIN_SWEI_RANDOM_SEED)
      out <- vector("list", n_tasks)
      t0 <- Sys.time()
      for (i in seq_len(n_tasks)) {
        out[[i]] <- .het_worker(het_tasks[[i]])
        if (i %% 20 == 0 || i == n_tasks) {
          elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
          rate    <- i / elapsed
          eta_sec <- if (is.finite(rate) && rate > 0) (n_tasks - i) / rate else NA_real_
          cat(sprintf("\r  Progress: %4d/%4d (%5.1f%%) | elapsed %5.0fs | est. remaining %5.0fs   ",
                      i, n_tasks, 100 * i / n_tasks, elapsed, ifelse(is.na(eta_sec), 0, eta_sec)))
          flush.console()
        }
      }
      cat("\n")
      out
    }
    
    ## ---- Decide serial vs. parallel execution ----
    ## Same rationale as loso_cv_bias_tiers(): PSOCK makeCluster() spawns
    ## brand-new R processes that each open a TCP socket back to this
    ## session, and on a machine where that handshake is blocked/slowed
    ## (firewall, antivirus process-spawn scanning, corporate loopback
    ## interception - all invisible from inside R, since PSOCK workers
    ## default to outfile = NULL, silently discarding anything a stuck
    ## worker would have printed), it can hang for HOURS with zero CPU and
    ## no error. Below HETEROGENEITY_SEQUENTIAL_MAX tasks, skip the cluster
    ## entirely; above it, still attempt parallel but with a bounded connect
    ## timeout, a real worker-output log, and an automatic serial fallback
    ## if it fails - the same two-layer defense already applied to LOSO-CV.
    if (!exists("HETEROGENEITY_SEQUENTIAL_MAX", inherits = TRUE)) HETEROGENEITY_SEQUENTIAL_MAX <- 20
    
    if (!isTRUE(ENABLE_PSOCK_PARALLEL) || n_tasks <= HETEROGENEITY_SEQUENTIAL_MAX) {
      cat(sprintf("[OK] %d task(s) is at/under HETEROGENEITY_SEQUENTIAL_MAX (%d): running\n",
                  n_tasks, HETEROGENEITY_SEQUENTIAL_MAX))
      cat("     sequentially rather than starting a parallel cluster.\n")
      results <- run_het_serial()
      
    } else {
      n_cores_het <- min(MAX_PSOCK_WORKERS, max(1L, detectCores() - 1L))
      cat(sprintf("[OK] Running heterogeneity tests in parallel (%d cores)\n", n_cores_het))
      cluster_log_het <- file.path(tempdir(), "heterogeneity_cluster_worker_log.txt")
      cat(sprintf("[OK] makeCluster(%d) starting (60s per-worker connect timeout; worker\n     stdout/stderr -> %s;\n     falls back to sequential execution automatically if this fails)...\n",
                  n_cores_het, cluster_log_het))
      flush.console()
      
      # ---- Robust cluster startup (same reasoning as loso_cv_bias_tiers()) ----
      # timeout = 60 bounds the per-worker socket handshake (default is 30
      # DAYS - effectively "hang forever", which is almost certainly what a
      # multi-hour stall right here means). outfile routes worker stdout/
      # stderr to a real file so a genuine connection failure leaves a
      # diagnosable trace instead of silence. If makeCluster() still fails
      # or times out, fall back to the identical serial loop above rather
      # than losing this diagnostic step to a cluster that will never come up.
      cl_het <- tryCatch(
        parallel::makeCluster(n_cores_het, timeout = 60, outfile = cluster_log_het),
        error = function(e) {
          cat(sprintf("[OK] makeCluster() failed or timed out (%s).\n     See worker log: %s\n",
                      conditionMessage(e), cluster_log_het))
          NULL
        }
      )
      
      if (is.null(cl_het)) {
        cat("[OK] Falling back to SEQUENTIAL execution for this run (identical per-task\n")
        cat("     computation, just without a cluster). Check the worker log above if\n")
        cat("     this keeps happening.\n")
        results <- run_het_serial()
        
      } else {
        clusterSetRNGStream(cl_het, iseed = BASIN_SWEI_RANDOM_SEED)
        clusterExport(cl_het, varlist = c("heterogeneity_test", "swe_matrix", "index_values",
                                          "ref_idx", "mon_ref", "HETEROGENEITY_NSIM", "H0_MC_CONF_LEVEL"),
                      envir = environment())
        clusterEvalQ(cl_het, library(lmomco))
        
        # ---- Manual load-balanced dispatch with a live progress readout ----
        # parLapply()/clusterApplyLB() block silently until the WHOLE batch is
        # done, which is exactly what produced the opaque stall. This reproduces
        # the same load-balanced scheduling using the low-level sendCall() /
        # recvOneResult() primitives that those functions are themselves built
        # on, so the master can print a line every time one (region, month) test
        # finishes instead of waiting on all 240 at once.
        n_workers  <- length(cl_het)
        results    <- vector("list", n_tasks)
        start_time <- Sys.time()
        
        n_initial <- min(n_workers, n_tasks)
        for (i in seq_len(n_initial)) {
          parallel:::sendCall(cl_het[[i]], .het_worker, list(het_tasks[[i]]), tag = i)
        }
        next_task <- n_initial + 1L
        completed <- 0L
        
        tryCatch({
          while (completed < n_tasks) {
            res <- parallel:::recvOneResult(cl_het)
            results[[res$tag]] <- res$value
            completed <- completed + 1L
            
            if (next_task <= n_tasks) {
              parallel:::sendCall(cl_het[[res$node]], .het_worker, list(het_tasks[[next_task]]), tag = next_task)
              next_task <- next_task + 1L
            }
            
            elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
            rate    <- completed / elapsed
            eta_sec <- if (is.finite(rate) && rate > 0) (n_tasks - completed) / rate else NA_real_
            cat(sprintf("\r  Progress: %4d/%4d (%5.1f%%) | elapsed %5.0fs | est. remaining %5.0fs   ",
                        completed, n_tasks, 100 * completed / n_tasks, elapsed,
                        ifelse(is.na(eta_sec), 0, eta_sec)))
            flush.console()
          }
        }, error = function(e) { try(stopCluster(cl_het), silent = TRUE); stop(e) })
        cat("\n")
        
        stopCluster(cl_het)
      }
    }
    
    heterogeneity_log <- do.call(rbind, results)
  } else {
    heterogeneity_log <- data.frame()
  }
  
  
  # ---- E.1 Turn the H0 flag into an actual pooling decision ----
  # Previously heterogeneity_log was written to disk and printed as a
  # diagnostic but never actually consulted when fitting the regional gamma
  # (every region/month pooled regardless of H0). Region/months flagged
  # H0 >= H0_REGION_FLAG_THRESHOLD are now excluded from regional pooling
  # entirely and fall through to per-pixel local fitting in
  # monthly_sspi_regional() instead - see heterogeneous_mat below and its use
  # in section G. NA/failed H0 tests (lmom_failed, fit_failed, sim_failed,
  # too_few_sites) are treated as heterogeneous too, since "unknown" is not
  # the same as "known homogeneous" and pooling on an untested sample
  # defeats the purpose of the screen.
  heterogeneous_mat <- matrix(TRUE, n_regions, 12)  # default: don't pool
  if (nrow(heterogeneity_log) > 0) {
    for (k in seq_len(nrow(heterogeneity_log))) {
      r <- heterogeneity_log$region[k]; m <- heterogeneity_log$month[k]
      h <- heterogeneity_log$H0[k]
      h_upper <- if ("H0_upper" %in% names(heterogeneity_log)) heterogeneity_log$H0_upper[k] else NA_real_
      h_decision <- if (is.finite(h_upper)) h_upper else h
      heterogeneous_mat[r, m] <- !(is.finite(h_decision) && h_decision < H0_REGION_FLAG_THRESHOLD)
    }
    n_flag <- sum(heterogeneity_log$H0 >= H0_REGION_FLAG_THRESHOLD, na.rm = TRUE)
    n_untested <- sum(!is.finite(heterogeneity_log$H0))
    cat(sprintf("[OK] Heterogeneity test: %d/%d region-months flagged H0 >= %.0f (possibly/definitely heterogeneous).\n",
                n_flag, nrow(heterogeneity_log), H0_REGION_FLAG_THRESHOLD))
    cat(sprintf("  %d additional region-months had no valid H0 (untested) and are also excluded from pooling.\n",
                n_untested))
    cat("  Flagged/untested region-months are EXCLUDED from regional pooling and fall back\n")
    cat("  to per-pixel local fitting (monthly_sspi) for that region/month - no regional-\n")
    cat("  empirical fallback is used. See rfa_diagnostics/ for the full H0 breakdown.\n")
  }
  
  ## ---- F. Bias correction vs. independent station SWE observations ----
  # Station correction already applied basin-wide (once) upstream, directly
  # to swe_matrix, before this RFA block runs - see the "BASIN-WIDE BIAS
  # CORRECTION" section above. Re-running fit_bias_correction() here would
  # double-correct every pixel, so it is intentionally NOT called. C_matrix
  # is fixed at identity so monthly_sspi_regional()'s `mu_corr <- mu * c_val`
  # step remains a no-op. fit_bias_correction() is left defined earlier in
  # the script (unused) in case a denser station network later makes
  # per-region quantile mapping viable.
  C_matrix <- matrix(1, nrow = max(n_regions, 1), ncol = 12)
  bias_correction_log <- data.frame()
  cat("[OK] Per-region station bias correction: skipped (basin-wide correction already applied upstream; see basin_bias_log).\n")
  
  ## ---- G. Pooled regional gamma fit + KS goodness-of-fit per region/month ----
  cat("\n===== FITTING POOLED REGIONAL GAMMA DISTRIBUTIONS =====\n")
  region_par <- vector("list", n_regions)
  gof_log <- data.frame()
  for (r in seq_len(n_regions)) {
    region_par[[r]] <- vector("list", 12)
    pix_r <- which(region_id == r)
    for (m in 1:12) {
      keep <- pix_r[!discordant_mat[pix_r, m]]
      cols <- ref_idx[mon_ref == m]
      if (length(keep) == 0 || length(cols) == 0) {
        region_par[[r]][[m]] <- list(par = NULL, pooled_sample = numeric(0),
                                     gof_pass = FALSE, status = "no_data")
        # Previously silent (no gof_log row) - logging it so gof_log always
        # has n_regions*12 rows and "missing" region-months don't have to be
        # inferred by diffing row counts against n_regions*12.
        gof_log <- rbind(gof_log, data.frame(region = r, month = m, n_pool = 0L,
                                             fit_status = "no_data", ks_D = NA_real_,
                                             ks_p = NA_real_, ks_ties = NA, gamma_used = FALSE))
        next
      }
      if (isTRUE(heterogeneous_mat[r, m])) {
        # Region/month failed (or never passed) the H0 heterogeneity screen -
        # do NOT pool. Leaving par = NULL and pooled_sample empty means
        # pool_ok is FALSE in monthly_sspi_regional(), so every pixel in this
        # region/month falls through to per-pixel local fitting instead of a
        # regional gamma OR a regional-empirical CDF.
        region_par[[r]][[m]] <- list(par = NULL, pooled_sample = numeric(0),
                                     gof_pass = FALSE, status = "heterogeneous_skip")
        # n_pool here is the pooled SAMPLE size the region/month WOULD have
        # had (matching the units of every other row's n_pool), not the site
        # count - a prior version of this log logged length(keep) (site
        # count, e.g. 36) here instead of the actual pooled sample size (e.g.
        # ~450-3000), which made heterogeneous_skip rows look like they had
        # far less data than they actually did.
        c_val_diag <- C_matrix[r, m]
        sub_diag  <- swe_matrix[keep, cols, drop = FALSE] * c_val_diag
        idxv_diag <- index_values[keep, m] * c_val_diag
        x_star_diag <- sweep(sub_diag, 1, idxv_diag, "/")
        n_pool_diag <- sum(is.finite(x_star_diag) & x_star_diag > 0)
        gof_log <- rbind(gof_log, data.frame(region = r, month = m, n_pool = n_pool_diag,
                                             fit_status = "heterogeneous_skip", ks_D = NA_real_,
                                             ks_p = NA_real_, ks_ties = NA, gamma_used = FALSE))
        next
      }
      c_val <- C_matrix[r, m]
      sub  <- swe_matrix[keep, cols, drop = FALSE] * c_val
      idxv <- index_values[keep, m] * c_val
      x_star <- sweep(sub, 1, idxv, "/")
      x_star <- x_star[is.finite(x_star) & x_star > 0]
      
      fit <- fit_regional_gamma(x_star)
      gof <- if (!is.null(fit$par)) {
        gof_test_gamma(x_star, fit$par, alpha = GOF_ALPHA)
      } else list(D = NA_real_, p_value = NA_real_, pass = FALSE)
      
      region_par[[r]][[m]] <- list(par = fit$par, pooled_sample = x_star,
                                   gof_pass = isTRUE(gof$pass), status = fit$status,
                                   n_pool = length(x_star))
      gof_log <- rbind(gof_log, data.frame(region = r, month = m, n_pool = length(x_star),
                                           fit_status = fit$status, ks_D = gof$D,
                                           ks_p = gof$p_value, ks_ties = isTRUE(gof$ties),
                                           gamma_used = isTRUE(gof$pass)))
    }
  }
  if (nrow(gof_log) > 0) {
    n_gamma <- sum(gof_log$gamma_used, na.rm = TRUE)
    n_het_skip <- sum(gof_log$fit_status == "heterogeneous_skip", na.rm = TRUE)
    cat(sprintf("[OK] Regional gamma fitting: %d/%d region-months pass KS goodness-of-fit (alpha=%.2f)\n",
                n_gamma, nrow(gof_log), GOF_ALPHA))
    cat(sprintf("  Region-months excluded from pooling for failing the H0 heterogeneity screen: %d\n",
                n_het_skip))
    # No regional-empirical fallback: every region-month that fails GoF or the
    # H0 screen falls straight through to per-pixel local fitting instead
    # (see monthly_sspi_regional / "Regional Empirical" is always 0% below).
    cat(sprintf("  Region-months failing GoF (but passing H0) that fall back to LOCAL fitting: %d\n",
                sum(!gof_log$gamma_used & gof_log$fit_status != "heterogeneous_skip" &
                      gof_log$n_pool >= 15, na.rm = TRUE)))
  }
  
  ## ---- H. Diagnostic outputs (mirrors paper's Table 3/6 + Fig. 6) ----
  cat("\n===== WRITING RFA DIAGNOSTIC OUTPUTS =====\n")
  rfa_dir <- file.path(out_dir, "rfa_diagnostics")
  if (!dir.exists(rfa_dir)) dir.create(rfa_dir, recursive = TRUE)
  if (nrow(discordancy_log) > 0)
    write.csv(discordancy_log, file.path(rfa_dir, "discordancy_flagged_pixels.csv"), row.names = FALSE)
  if (nrow(heterogeneity_log) > 0)
    write.csv(heterogeneity_log, file.path(rfa_dir, "heterogeneity_H0_by_region_month.csv"), row.names = FALSE)
  if (nrow(gof_log) > 0)
    write.csv(gof_log, file.path(rfa_dir, "gamma_goodness_of_fit_by_region_month.csv"), row.names = FALSE)
  if (nrow(bias_correction_log) > 0)
    write.csv(bias_correction_log, file.path(rfa_dir, "bias_correction_coefficients.csv"), row.names = FALSE)
  
  region_raster <- rast(swe[[1]])
  values(region_raster) <- region_id
  png(file.path(rfa_dir, "regions_map.png"), width = 1200, height = 900, res = 150)
  plot(region_raster, main = sprintf("SSPI-1 homogeneous regions (%s)", regions_built$method),
       col = grDevices::hcl.colors(max(n_regions, 1), "Dark 3"), axes = FALSE, box = FALSE)
  if (!is.null(basin)) plot(basin, add = TRUE, border = "black", lwd = 1.5)
  dev.off()
  # region_raster is only needed for this PNG - it's a fresh terra SpatRaster
  # (external C++ pointer), and unlike swe/scf_stack/basin it is never
  # rehydrated from disk on reload. Left in the environment, it gets swept
  # into the end-of-script saveRDS() checkpoint, where its pointer is either
  # already stale or unserializable -> "external pointer is not valid" at
  # save time. Drop it now; nothing downstream needs it.
  rm(region_raster)
  cat(sprintf("[OK] RFA diagnostics written to: %s\n", normalizePath(rfa_dir)))
  cat("============================================================\n")
  
} else if (scf_method == "2" && !REGIONALIZE_SSPI) {
  cat("\n- REGIONALIZE_SSPI = FALSE: using original pixel-by-pixel SSPI-1 fitting.\n")
  region_id <- NULL; index_values <- NULL; C_matrix <- NULL
  region_par <- NULL; discordant_mat <- NULL
}


##############################################
# ---- STAGE 2 CHECKPOINT: save everything for script 03 ----
##############################################
cat("\n============================================================\n")
cat("STAGE 2 COMPLETE - saving workspace checkpoint for script 03\n")
cat("============================================================\n")

# ---- Persist SpatRaster/SpatVector objects via terra's own I/O ----
# Same reason as script 01's checkpoint: saveRDS() cannot round-trip terra's
# external C++ pointers ("external pointer is not valid" at save/load time).
# swe/scf_stack/basin were rehydrated from stage1's .tif/.gpkg files at the
# top of this script (lines 42-50) and are unchanged here, so simply re-copy
# those same files forward as this stage's checkpoint rather than writing
# from the (possibly stale, e.g. after any crop/mask) in-memory objects -
# consistent with how they were produced. EXCLUDE them from the RDS payload
# below so nothing overwrites the rehydrated versions when script 03 unpacks
# the RDS, and so the RDS write itself doesn't touch a live pointer.
file.copy(file.path(STAGE_DIR, "stage1_swe.tif"),
          file.path(STAGE_DIR, "stage2_swe.tif"), overwrite = TRUE)
if (!is.null(scf_stack)) {
  file.copy(file.path(STAGE_DIR, "stage1_scf_stack.tif"),
            file.path(STAGE_DIR, "stage2_scf_stack.tif"), overwrite = TRUE)
}
if (!is.null(basin)) {
  file.copy(file.path(STAGE_DIR, "stage1_basin.gpkg"),
            file.path(STAGE_DIR, "stage2_basin.gpkg"), overwrite = TRUE)
}

# Defensive terra-pointer sweep. Do not rely on manually rm()'ing every
# temporary raster/vector: automatically exclude every live Spat* object from
# the RDS payload, then remove the nonessential ones before RStudio can inspect
# or autosave them. Core swe/scf_stack/basin remain available until the sweep.
# Recursive (not just a top-level inherits() check) so a Spat* object buried
# inside a list or data.frame column - which the old top-level-only check
# would silently miss and let flow straight into saveRDS() - gets caught too.
has_nested_spat <- function(x, depth = 0L, max_depth = 6L) {
  if (depth > max_depth) return(FALSE)
  if (inherits(x, c("SpatRaster", "SpatVector",
                    "SpatRasterDataset", "SpatRasterCollection"))) return(TRUE)
  if (is.list(x)) return(any(vapply(x, has_nested_spat, logical(1),
                                    depth = depth + 1L, max_depth = max_depth)))
  FALSE
}
all_stage2_names <- ls(all.names = TRUE, envir = environment())
spat_leftover <- Filter(function(nm) {
  obj <- tryCatch(get(nm, envir = environment(), inherits = FALSE), error = function(e) NULL)
  !is.null(obj) && has_nested_spat(obj)
}, all_stage2_names)
if (length(spat_leftover) > 0L) {
  cat(sprintf("[CHECKPOINT] Excluding %d live terra Spat* object(s) from stage2.rds: %s\n",
              length(spat_leftover), paste(spat_leftover, collapse = ", ")))
}

stage2_objects <- setdiff(all_stage2_names,
                          unique(c("STAGE_DIR", "stage1_path", "stage2_path",
                                   "has_nested_spat", "all_stage2_names",
                                   "spat_leftover", spat_leftover)))
stage2_path <- file.path(STAGE_DIR, "stage2.rds")
saveRDS(mget(stage2_objects, envir = environment(), inherits = FALSE), stage2_path, compress = "xz")

# Verify the checkpoint immediately in an isolated object before releasing
# terra pointers. This catches a corrupt/incomplete RDS while the stage is
# still running, rather than discovering it at script 03 startup.
stage2_verify <- readRDS(stage2_path)
if (!is.list(stage2_verify) || length(stage2_verify) != length(stage2_objects)) {
  stop("stage2.rds verification failed: object count/type mismatch after saveRDS().", call. = FALSE)
}
rm(stage2_verify)

# Release all live terra objects before the script returns to the RStudio top
# level. This prevents Environment-pane refresh or .RData autosave from
# dereferencing stale external pointers after the checkpoint has succeeded.
if (length(spat_leftover) > 0L) rm(list = spat_leftover, envir = environment())
gc(verbose = FALSE)
gc(verbose = FALSE)  # second pass: catches wrapper objects whose finalizer
# only actually freed the underlying pointer during the
# first pass, before RStudio's Environment pane refresh
# or .RData autosave gets a chance to touch them

cat(sprintf("[OK] Saved %d object(s) to: %s\n", length(stage2_objects), normalizePath(stage2_path)))
cat(sprintf("[OK] File size: %.1f MB\n", file.size(stage2_path) / 1024^2))
cat("\nNext: run 03_main_parallel_calc.R\n")
cat("(Once 03 has run successfully, you can delete BOTH stage1.rds and\n")
cat(" stage2.rds to save space - 03's checkpoint (stage3.rds) is small and\n")
cat(" is all that scripts 03->04 need going forward.)\n")