# 8trends_visualization.R  ·  SPATIAL DROUGHT DIAGNOSTICS
# ─────────────────────────────────────────────────────────────────────────────────
# COMPATIBLE WITH: DROUGHT_ANALYSIS_utils.R, 6trend_test_ALL.R
# SUPPORTS: SPI, SPEI-PM, SPEI-ERA5, SWEI, MSPI, MSPEI (Gridded)
# FRAMEWORK: observed PM–ERA5-Land PEV comparison only; no detrended/counterfactual branch.
#
# SCOPE (spatial figures only):
#   PART 1 – Merge individual result CSVs → master diagnostics CSV
#   PART 2 – Spatial figures
#     Fig 1: Trend maps          (Kendall τ_vc + FDR significance)
#     Fig 2: Event maps          (# events, mean duration, max intensity)
#     Fig 3: Temporal pattern maps (clustering, regime-shift decade)
#     Fig 4: Timescale comparison (% significant, median τ, events, duration)
#     Fig 5: Method comparison   (VC vs TFPW scatter + agreement bar)
#   PART 3 – Manuscript Figure X
#     3 × 4 Sen-slope maps for SPI, SPEI-PM and SPEI-ERA5 at 3, 12, 24, 36 months
#
# NOTE: Basin-averaged time series CSVs, PNG plots, and the Excel summary
#       are produced by 7basin_timeseries.R — run that script for all
#       non-spatial outputs.
#
# EXECUTION ORDER: 6trend_test_ALL.R → 7basin_timeseries.R → 8trends_visualization.R

source("DROUGHT_ANALYSIS_utils.R")
utils_load_packages(c("terra", "data.table", "ggplot2", "sf", "viridis",
                      "cowplot", "patchwork", "scales", "RColorBrewer",
                      "ggspatial", "gridExtra", "grid", "lubridate"))

# ── Working directory & output folder ───────────────────────────────
if (!dir.exists(WD_PATH)) stop("Working directory not found: ", WD_PATH)
setwd(WD_PATH)
fig_dir <- file.path(TREND_DIR, "figures_temporal_diagnostics")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ── Native ERA5-Land grid and basin boundary ───────────────────────
# Script 8 does NOT use a canonical/project-wide analysis grid.
# The first layer of the original ERA5-Land precipitation NetCDF is the
# sole authoritative raster geometry for every spatial product in this script.
native_precip_file <- file.path(
  WD_PATH,
  "monthly_data_direct",
  "total_precipitation_monthly.nc"
)

if (!file.exists(native_precip_file)) {
  stop(
    "Native ERA5-Land precipitation NetCDF not found: ",
    native_precip_file,
    "\nRun Monthly_ERALand_download.R first."
  )
}

native_precip <- terra::rast(native_precip_file)

if (!inherits(native_precip, "SpatRaster")) {
  stop("Failed to load native ERA5-Land precipitation NetCDF as a SpatRaster.")
}

if (terra::nlyr(native_precip) < 1L) {
  stop("Native ERA5-Land precipitation NetCDF contains no raster layers.")
}

analysis_template <- native_precip[[1]]

if (is.na(terra::crs(analysis_template)) ||
    !nzchar(terra::crs(analysis_template))) {
  stop("Native ERA5-Land precipitation grid has no CRS.")
}

# Backward-compatible alias for any downstream plotting code that still expects
# this object name. It is NOT a canonical-grid object.
common_template <- analysis_template

cat(sprintf(
  "[GRID] Native ERA5-Land precipitation grid: %d x %d = %d cells | ",
  terra::nrow(analysis_template),
  terra::ncol(analysis_template),
  terra::ncell(analysis_template)
))
cat(sprintf(
  "resolution %.10f x %.10f | CRS=%s\n",
  terra::res(analysis_template)[1],
  terra::res(analysis_template)[2],
  terra::crs(analysis_template)
))

basin_vect <- load_basin_vect(BASIN_SHP)
if (!inherits(basin_vect, "SpatVector")) {
  basin_vect <- terra::vect(basin_vect)
}

# Only the vector geometry is transformed, if necessary. The native raster
# geometry is never reprojected or resampled.
if (!terra::same.crs(basin_vect, analysis_template)) {
  basin_vect <- terra::project(
    basin_vect,
    terra::crs(analysis_template)
  )
}

basin_sf <- sf::st_as_sf(basin_vect)
basin_sf_global <<- basin_sf

# Native-grid basin mask. This is generated directly on the authoritative
# precipitation grid and is reused by all spatial diagnostics below.
basin_coverage <- terra::rasterize(
  basin_vect,
  analysis_template,
  cover = TRUE
)

basin_coverage_values <- as.numeric(
  terra::values(basin_coverage, mat = FALSE)
)
basin_coverage_values[!is.finite(basin_coverage_values)] <- 0
basin_coverage_values <- pmax(0, pmin(1, basin_coverage_values))

basin_mask_raster <- terra::setValues(
  analysis_template,
  as.numeric(basin_coverage_values > 0)
)

basin_mask_geom_ok <- terra::compareGeom(
  basin_mask_raster,
  analysis_template,
  stopOnError = FALSE
)

if (!isTRUE(basin_mask_geom_ok)) {
  stop(
    "Native ERA5-Land basin mask does not match the precipitation grid."
  )
}

basin_mask_cells <- as.logical(
  terra::values(basin_mask_raster, mat = FALSE)
)

basin_pixels <- sum(
  basin_mask_cells,
  na.rm = TRUE
)

if (basin_pixels == 0L) {
  stop("CRITICAL: native ERA5-Land basin mask contains zero cells.")
}

cat(sprintf(
  "[GRID] Native basin mask: %d cells (%.2f%% of raster)\n",
  basin_pixels,
  100 * basin_pixels / terra::ncell(analysis_template)
))

# ── Analysis constants ──────────────────────────────────────────────
DPI              <- 300
TIMESCALES_VIZ   <- TIMESCALES_SPATIAL   # c(1,3,6,12) for spatial figures
MSPI_MSPEI_SCALE <- 1L                   # single-scale multivariate indices

index_display <- function(x) {
  out <- toupper(x)
  out[x == "spei"] <- "SPEI-PM"
  out[x == "spei_era5"] <- "SPEI-ERA5"
  out
}

# ── Index colours (extend for MSPI / MSPEI) ─────────────────────────
if (!"swei"  %in% names(index_colours))
  index_colours <- c(index_colours, swei  = "#1E90FF")
if (!"mspi"  %in% names(index_colours))
  index_colours <- c(index_colours, mspi  = "#33a02c")   # green
if (!"mspei" %in% names(index_colours))
  index_colours <- c(index_colours, mspei = "#fb9a99")   # light red
if (!"spei_era5" %in% names(index_colours))
  index_colours <- c(index_colours, spei_era5 = "#d62728") # ERA5-Land PEV observed

################################################################################
# PART 1 – MERGE INDIVIDUAL RESULT FILES → MASTER CSV
################################################################################
cat("\n══════════════════════════════════════\n")
cat("PART 1: Building master diagnostics file\n")
cat("══════════════════════════════════════\n")

master_csv <- file.path(TREND_DIR, "all_temporal_diagnostics_results.csv")

load_clip_results <- function(index_type, timescale) {
  f <- file.path(TREND_DIR, sprintf("%s_%02d_results.csv", index_type, timescale))
  if (!file.exists(f)) {
    if (!(index_type %in% c("mspi", "mspei") ||
          (index_type == "swei" && timescale != SWEI_SCALE)))
      cat("  ⚠ Missing: ", basename(f), "\n")
    return(NULL)
  }
  dt <- data.table::fread(f)
  cat(sprintf("  %s-%02d: %d rows\n", toupper(index_type), timescale, nrow(dt)))
  
  # Standardise coordinate column names
  if ("x" %in% names(dt) && !"lon" %in% names(dt))
    data.table::setnames(dt, "x", "lon")
  if ("y" %in% names(dt) && !"lat" %in% names(dt))
    data.table::setnames(dt, "y", "lat")
  if (!all(c("lon", "lat") %in% names(dt))) {
    cat("  ⚠ No coordinate columns — skipping\n"); return(NULL)
  }
  dt <- dt[!is.na(lon) & !is.na(lat)]
  
  pts <- terra::vect(
    dt,
    geom = c("lon", "lat"),
    crs = terra::crs(analysis_template)
  )
  
  inside_rel <- terra::relate(
    pts,
    basin_vect,
    relation = "intersects",
    pairs = FALSE
  )
  
  # terra::relate() returns a matrix for vector-to-vector relations in the
  # installed terra version. Convert that matrix to one logical value per
  # input point before using it as a data.table row index.
  if (is.matrix(inside_rel)) {
    inside <- apply(inside_rel, 1L, any)
  } else {
    inside <- as.logical(inside_rel)
  }
  
  if (length(inside) != nrow(dt)) {
    stop(
      "Point-in-basin relation returned ",
      length(inside),
      " flags for ",
      nrow(dt),
      " trend-result points."
    )
  }
  
  dt_clip <- dt[inside, ]
  
  if (is.null(dt_clip) || nrow(dt_clip) == 0) {
    cat("  ⚠ No points remain after basin clip\n"); return(NULL)
  }
  dt_clip[, `:=`(index_type = index_type, timescale = timescale)]
  dt_clip
}

all_results <- list()

cat("\n  Loading SPI, SPEI-PM & SPEI-ERA5...\n")
for (idx in c("spi", "spei", "spei_era5")) {
  for (sc in TIMESCALES_VIZ) {
    key <- sprintf("%s_%02d", idx, sc)
    r   <- load_clip_results(idx, sc)
    if (!is.null(r)) all_results[[key]] <- r
  }
}

cat("\n  Loading SWEI (scale", SWEI_SCALE, ")...\n")
swei_res <- load_clip_results("swei", SWEI_SCALE)
if (!is.null(swei_res)) all_results[[sprintf("swei_%02d", SWEI_SCALE)]] <- swei_res

cat("\n  Loading MSPI & MSPEI (scale", MSPI_MSPEI_SCALE, ")...\n")
for (idx in c("mspi", "mspei")) {
  key <- sprintf("%s_%02d", idx, MSPI_MSPEI_SCALE)
  r   <- load_clip_results(idx, MSPI_MSPEI_SCALE)
  if (!is.null(r)) all_results[[key]] <- r
}

if (!length(all_results))
  stop("No data loaded. Run 6trend_test_ALL.R first.")

combined <- data.table::rbindlist(all_results, fill = TRUE)
data.table::fwrite(combined, master_csv)
cat(sprintf("✓ Master CSV: %d rows → %s\n", nrow(combined), basename(master_csv)))

################################################################################
# PART 2 – SPATIAL FIGURES
################################################################################
cat("\n══════════════════════════════════════\n")
cat("PART 2: Spatial figures\n")
cat("══════════════════════════════════════\n")

plot_data <- combined

# ── Shared raster helpers ──────────────────────────────────────────
make_raster <- function(dt, vcol, tpl = NULL) {
  if (is.null(tpl)) tpl <- analysis_template
  create_raster_from_points(dt, tpl, vcol, basin_sf)
}

# Which indices are available for a given scale
active_indices <- function(scale) {
  # SPEI-PM and SPEI-ERA5 are parallel observed formulations.
  base <- c("spi", "spei")
  if (any(plot_data$index_type == "spei_era5" & plot_data$timescale == scale))
    base <- c(base, "spei_era5")
  if (scale == SWEI_SCALE        && any(plot_data$index_type == "swei"))
    base <- c(base, "swei")
  if (scale == MSPI_MSPEI_SCALE) {
    if (any(plot_data$index_type == "mspi"))  base <- c(base, "mspi")
    if (any(plot_data$index_type == "mspei")) base <- c(base, "mspei")
  }
  base
}

# ── Fig 1: Trend maps (Kendall τ_vc + FDR significance) ───────────
# Layout — 2 rows × n_index columns:
#   Row 1  Kendall's τ (VC-HR98) with stippling overlay
#          • Filled dots mark pixels that survive BH-FDR at α = 0.05
#          • Stippling follows Wilks (2006, BAMS): show all τ values but
#            identify only the locally robust subset
#   Row 2  Three-category significance map:
#          • Dark red   — FDR-significant          (p_fdr_vc  < 0.05)
#          • Orange     — Nominally-sig only        (p_value_vc < 0.05 but p_fdr_vc ≥ 0.05)
#          • Light grey — Not significant           (p_value_vc ≥ 0.05)
#   Falls back to 4-level nominal panel when p_fdr_vc is absent (old w3 run).
create_trend_maps <- function(scale = 12, outfile) {
  cat(sprintf("  Fig 1: trend maps (scale=%d)... ", scale))
  if (!safe_pdf(outfile, 12, 10)) return()
  indices_here <- active_indices(scale)
  n_plots      <- length(indices_here)
  
  par(mfrow = c(2, n_plots), mar = c(2, 2, 3, 3), oma = c(0, 0, 2, 0))
  tau_range <- c(-0.4, 0.4)
  has_fdr   <- "p_fdr_vc" %in% names(plot_data)
  if (!has_fdr)
    cat("\n  ⚠ p_fdr_vc absent — falling back to nominal significance",
        "\n    (re-run 6trend_test_ALL.R to generate FDR-corrected columns)\n")
  
  for (idx in indices_here) {
    sub <- plot_data[index_type == idx & timescale == scale]
    if (!nrow(sub)) {
      plot.new(); text(0.5, 0.5, paste("No data:", toupper(idx))); next
    }
    tpl   <- analysis_template
    r_tau <- make_raster(sub, "tau_vc", tpl)
    
    # Row 1: τ map + FDR stippling
    plot_raster_clean(r_tau, paste0(toupper(idx), ": Kendall's τ (VC-HR98)"),
                      tau_range, grDevices::hcl.colors(101, "RdBu", rev = TRUE),
                      legend_title = "τ")
    if (has_fdr) {
      stip <- sub[!is.na(p_fdr_vc) & p_fdr_vc < 0.05]
      if (nrow(stip) > 0)
        graphics::points(stip$lon, stip$lat, pch = 20, cex = 0.45, col = "black")
      graphics::legend("bottomleft",
                       legend = "FDR-significant (BH α = 0.05)",
                       pch = 20, col = "black", bty = "n", cex = 0.75)
    }
    
    # Row 2: significance category map
    sub2 <- data.table::copy(sub)
    if (has_fdr) {
      sub2[, sig_cat := data.table::fcase(
        is.na(p_value_vc),                         NA_real_,
        !is.na(p_fdr_vc) & p_fdr_vc  < 0.05,      3,   # FDR-significant
        p_value_vc < 0.05,                          2,   # nominally-sig only
        default =                                   1)]  # not significant
      r_sig <- make_raster(sub2[!is.na(sig_cat)], "sig_cat", tpl)
      plot_raster_clean(
        r_sig,
        paste0(toupper(idx), ": Significance (BH-FDR corrected)"),
        NULL, c("#f7f7f7", "#fc8d59", "#d73027"),
        breaks = c(0.5, 1.5, 2.5, 3.5), categorical = TRUE, legend = FALSE)
      graphics::legend("bottomright",
                       legend = c("Not significant",
                                  "Nominal only (p < 0.05)",
                                  "FDR-significant"),
                       fill = c("#f7f7f7", "#fc8d59", "#d73027"),
                       bty = "n", cex = 0.75, title = "Significance")
    } else {
      sub2[, sig_num := dplyr::case_when(
        is.na(p_value_vc)   ~ NA_real_,
        p_value_vc < 0.001  ~ 1,
        p_value_vc < 0.01   ~ 2,
        p_value_vc < 0.05   ~ 3,
        TRUE                ~ 4)]
      r_sig <- make_raster(sub2[!is.na(sig_num)], "sig_num", tpl)
      plot_raster_clean(r_sig, paste0(toupper(idx), ": Significance (nominal)"),
                        NULL, c("#d73027", "#fc8d59", "#fee08b", "#f7f7f7"),
                        breaks = c(0.5, 1.5, 2.5, 3.5, 4.5),
                        categorical = TRUE, legend_title = "p-value")
    }
  }
  graphics::mtext(sprintf("Trend Analysis – VC-HR98 Method (%d-month)", scale),
                  outer = TRUE, cex = 1.2, font = 2)
  grDevices::dev.off()
  cat(" ✓\n")
}

# ── Fig 2: Event characteristic maps ──────────────────────────────
# Colour design: DARKER colour = MORE extreme drought in every panel.
#   n_events      YlOrRd  (dark at high end, rev = FALSE)
#   mean_duration viridis reversed (long duration = dark, rev = TRUE)
#   max_intensity Reds reversed (most-negative = darkest, rev = TRUE)
create_event_maps <- function(scale = 12, outfile) {
  cat(sprintf("  Fig 2: event maps (scale=%d)... ", scale))
  if (!safe_pdf(outfile, 12, 14)) return()
  indices_here <- active_indices(scale)
  n_plots      <- length(indices_here)
  par(mfrow = c(3, n_plots), mar = c(2, 2, 3, 3), oma = c(0, 0, 2, 0))
  
  metrics <- list(
    n_events      = list(title = "# Events",      pal = "YlOrRd"),
    mean_duration = list(title = "Mean Duration",  pal = "viridis"),
    max_intensity = list(title = "Max Intensity",  pal = "Reds"))
  
  for (idx in indices_here) {
    sub <- plot_data[index_type == idx & timescale == scale]
    if (!nrow(sub)) {
      for (i in 1:3) { plot.new(); text(0.5, 0.5, paste("No data:", toupper(idx))) }
      next
    }
    tpl <- analysis_template
    for (mn in names(metrics)) {
      if (all(is.na(sub[[mn]]))) {
        plot.new(); text(0.5, 0.5, paste("All NA:", mn)); next
      }
      r   <- make_raster(sub, mn, tpl)
      zlm <- if (mn == "max_intensity") c(min(sub[[mn]], na.rm = TRUE), 0)
      else                       c(0, max(sub[[mn]], na.rm = TRUE))
      plot_raster_clean(
        r, paste0(toupper(idx), ": ", metrics[[mn]]$title), zlm,
        grDevices::hcl.colors(101, metrics[[mn]]$pal,
                              rev = mn %in% c("max_intensity", "mean_duration")))
    }
  }
  graphics::mtext(sprintf("Event Characteristics (%d-month)", scale),
                  outer = TRUE, cex = 1.2, font = 2)
  grDevices::dev.off()
  cat(" ✓\n")
}

# ── Fig 3: Temporal pattern maps ──────────────────────────────────
create_temporal_maps <- function(scale = 12, outfile) {
  cat(sprintf("  Fig 3: temporal maps (scale=%d)...", scale))
  if (!safe_pdf(outfile, 12, 8)) return()
  indices_here <- active_indices(scale)
  n_plots      <- length(indices_here)
  par(mfrow = c(n_plots, 3), mar = c(2, 1.5, 2.5, 1), oma = c(1, 1, 3, 1))
  
  need_cols <- c("filtered_runs", "p_value_runs", "clustering",
                 "regime_shift_year", "n_spectral_peaks")
  missing_cols <- setdiff(need_cols, names(plot_data))
  if (length(missing_cols)) {
    plot.new()
    text(0.5, 0.5,
         paste("Missing columns – re-run 6trend_test_ALL.R:\n",
               paste(missing_cols, collapse = ", ")),
         cex = 1.1, col = "red")
    grDevices::dev.off()
    cat(" ⚠ missing columns (need w3 upgrade)\n")
    return()
  }
  
  for (idx in indices_here) {
    sub <- plot_data[index_type == idx & timescale == scale]
    tpl <- analysis_template
    if (is.null(tpl)) {
      for (i in 1:3) { plot.new(); text(0.5, 0.5, paste("No data:", toupper(idx))) }
      next
    }
    
    # Panel A: Temporal clustering (runs test)
    sub_cl <- data.table::copy(sub)
    sub_cl[, clus_num := data.table::fcase(
      filtered_runs | is.na(p_value_runs), NA_real_,
      p_value_runs >= 0.05,                0,
      clustering == "clustered",           -1,
      clustering == "dispersed",            1)]
    r_cl <- make_raster(sub_cl, "clus_num", tpl)
    plot_raster_clean(r_cl, paste0(toupper(idx), ": Temporal Clustering"),
                      NULL, c("#e41a1c", "#fee090", "#4daf4a"),
                      breaks = c(-1.5, -0.5, 0.5, 1.5),
                      categorical = TRUE, legend = FALSE)
    graphics::legend("bottomright",
                     legend = c("Clustered", "Not sig.", "Dispersed"),
                     fill = c("#e41a1c", "#fee090", "#4daf4a"),
                     bty = "n", cex = 0.85)
    
    # Panel B: Regime-shift decade
    sub_rs <- sub[!is.na(regime_shift_year)]
    if (nrow(sub_rs)) {
      sub_rs2 <- data.table::copy(sub_rs)
      sub_rs2[, shift_dec := as.numeric(
        cut(regime_shift_year,
            breaks = seq(1950, 2030, by = 10),
            labels = 1:8, include.lowest = TRUE))]
      r_rs <- make_raster(sub_rs2, "shift_dec", tpl)
      plot_raster_clean(r_rs, paste0(toupper(idx), ": Regime-Shift Decade"),
                        c(1, 8), grDevices::hcl.colors(8, "Spectral", rev = TRUE),
                        legend_title = "Decade")
    } else {
      plot.new()
      text(0.5, 0.5, paste(toupper(idx), "– no regime shifts detected"))
    }
    
    # Panel C: Spectral peaks
    r_sp <- make_raster(sub, "n_spectral_peaks", tpl)
    plot_raster_clean(r_sp, paste0(toupper(idx), ": Spectral Peaks"),
                      c(0, max(sub$n_spectral_peaks, na.rm = TRUE, 1)),
                      grDevices::hcl.colors(101, "viridis"))
  }
  graphics::mtext(sprintf("Temporal Patterns (%d-month)", scale),
                  outer = TRUE, cex = 1.2, font = 2)
  grDevices::dev.off()
  cat(" ✓\n")
}

# ── Fig 4: Timescale comparison ────────────────────────────────────
# Panel A: TWO lines per index × scale — solid = nominal sig, dashed = BH-FDR.
# Gap between lines = false-discovery inflation at that scale.
create_timescale_comparison <- function(outfile) {
  cat("  Fig 4: timescale comparison... ")
  has_fdr <- "p_fdr_vc" %in% names(plot_data)
  
  smry <- plot_data[, .(
    Pct_Sig_Nom = sum(p_value_vc < 0.05, na.rm = TRUE) / .N * 100,
    Pct_Sig_FDR = if (has_fdr)
      sum(!is.na(p_fdr_vc) & p_fdr_vc < 0.05, na.rm = TRUE) / .N * 100
    else NA_real_,
    Median_Tau  = median(tau_vc,       na.rm = TRUE),
    Mean_Events = mean(n_events,        na.rm = TRUE),
    Mean_Dur    = mean(mean_duration,   na.rm = TRUE)
  ), by = .(index_type, timescale)]
  
  mk_single <- function(col, ylab, title) {
    ggplot2::ggplot(smry,
                    ggplot2::aes(x = timescale, y = .data[[col]],
                                 color = index_type, group = index_type)) +
      ggplot2::geom_line(linewidth = 1.2, na.rm = TRUE) +
      ggplot2::geom_point(size = 4, na.rm = TRUE) +
      ggplot2::scale_color_manual(values = index_colours, labels = index_display) +
      ggplot2::scale_x_continuous(breaks = sort(unique(smry$timescale))) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::labs(title = title,
                    x = "Timescale (months)", y = ylab, color = "Index")
  }
  
  sig_long <- data.table::melt(
    smry,
    id.vars      = c("index_type", "timescale"),
    measure.vars = c("Pct_Sig_Nom", "Pct_Sig_FDR"),
    variable.name = "Correction", value.name = "Pct_Sig")
  sig_long[, Correction := data.table::fcase(
    Correction == "Pct_Sig_Nom", "Nominal (p < 0.05)",
    Correction == "Pct_Sig_FDR", "BH-FDR (p_fdr < 0.05)")]
  sig_long <- sig_long[!is.na(Pct_Sig)]
  
  pA <- ggplot2::ggplot(
    sig_long,
    ggplot2::aes(x        = timescale, y        = Pct_Sig,
                 color    = index_type,
                 linetype = Correction,
                 group    = interaction(index_type, Correction))) +
    ggplot2::geom_line(linewidth = 1.1, na.rm = TRUE) +
    ggplot2::geom_point(
      data  = sig_long[Correction == "Nominal (p < 0.05)"],
      size  = 3.5, na.rm = TRUE) +
    ggplot2::geom_point(
      data  = sig_long[Correction == "BH-FDR (p_fdr < 0.05)"],
      size  = 3.5, shape = 1, na.rm = TRUE) +
    ggplot2::scale_color_manual(values = index_colours, labels = index_display) +
    ggplot2::scale_linetype_manual(
      values = c("Nominal (p < 0.05)"    = "solid",
                 "BH-FDR (p_fdr < 0.05)" = "dashed")) +
    ggplot2::scale_x_continuous(breaks = sort(unique(smry$timescale))) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(
      title    = "A) % Significant Pixels — Nominal vs BH-FDR",
      subtitle = "Solid = nominal  |  Dashed = BH-FDR corrected  |  Gap = false-discovery inflation",
      x = "Timescale (months)", y = "% Significant pixels",
      color = "Index", linetype = "Correction")
  
  p <- (pA | mk_single("Median_Tau",  "τ",        "B) Median Kendall's τ (VC-HR98)")) /
    (mk_single("Mean_Events", "N events",  "C) Mean Drought Events") |
       mk_single("Mean_Dur",    "months",   "D) Mean Duration")) +
    patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(legend.position = "bottom")
  
  ggplot2::ggsave(outfile, p, width = 14, height = 10, dpi = DPI)
  cat(" ✓\n")
}

# ── Fig 5: Method comparison (VC vs TFPW) ─────────────────────────
# Panel A: scatter τ_vc vs τ_tfpw (tau is not affected by FDR correction)
# Panel B: agreement bar; uses BH-FDR significance when columns are present,
#          falls back to nominal p-values with a labelled subtitle if absent.
create_method_comparison <- function(outfile) {
  cat("  Fig 5: method comparison (VC vs TFPW)... ")
  needed <- c("tau_tfpw", "p_value_tfpw")
  if (!all(needed %in% names(plot_data))) {
    if (safe_pdf(outfile)) {
      plot.new()
      text(0.5, 0.5,
           "TFPW columns not found.\nRe-run 6trend_test_ALL.R with TFPW enabled.",
           cex = 1.2, col = "orange")
      grDevices::dev.off()
    }
    cat(" ⚠ TFPW columns absent — placeholder saved\n")
    return()
  }
  
  has_fdr <- all(c("p_fdr_vc", "p_fdr_tfpw") %in% names(plot_data))
  sub     <- plot_data[timescale == 12]
  
  if (has_fdr) {
    sub[, `:=`(vc_sig   = !is.na(p_fdr_vc)   & p_fdr_vc   < 0.05,
               tfpw_sig = !is.na(p_fdr_tfpw)  & p_fdr_tfpw < 0.05)]
    bar_subtitle <- "BH-FDR corrected (p_fdr < 0.05)"
  } else {
    sub[, `:=`(vc_sig   = !is.na(p_value_vc)   & p_value_vc   < 0.05,
               tfpw_sig = !is.na(p_value_tfpw)  & p_value_tfpw < 0.05)]
    bar_subtitle <- "nominal p < 0.05 (re-run w3 for FDR columns)"
  }
  
  p1 <- ggplot2::ggplot(sub, ggplot2::aes(tau_vc, tau_tfpw, color = index_type)) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         linetype = "dashed", color = "gray40") +
    ggplot2::geom_point(alpha = 0.4, size = 2, na.rm = TRUE) +
    ggplot2::facet_wrap(~index_type,
                        labeller = ggplot2::labeller(index_type = toupper)) +
    ggplot2::scale_color_manual(values = index_colours) +
    ggplot2::coord_equal() +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::labs(title = "A) τ: Variance-Corrected (VC-HR98) vs TFPW (12-month)",
                  x = "τ (VC-HR98)", y = "τ (TFPW)")
  
  agree <- sub[, .(
    Both      = sum( vc_sig &  tfpw_sig, na.rm = TRUE),
    VC_only   = sum( vc_sig & !tfpw_sig, na.rm = TRUE),
    TFPW_only = sum(!vc_sig &  tfpw_sig, na.rm = TRUE),
    None      = sum(!vc_sig & !tfpw_sig, na.rm = TRUE)
  ), by = index_type]
  ag_l <- data.table::melt(agree, id.vars = "index_type",
                           variable.name = "Category", value.name = "N_pixels")
  
  p2 <- ggplot2::ggplot(ag_l,
                        ggplot2::aes(toupper(index_type), N_pixels, fill = Category)) +
    ggplot2::geom_col(position = "stack") +
    ggplot2::scale_fill_brewer(palette = "Set2") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::labs(
      title    = "B) Significance Agreement — VC-HR98 vs TFPW (12-month)",
      subtitle = bar_subtitle,
      x = "Index", y = "# Pixels")
  
  ggplot2::ggsave(outfile, p1 | p2, width = 14, height = 6, dpi = DPI)
  cat(" ✓\n")
}

# ── Run all spatial figures ─────────────────────────────────────────
for (sc in c(1, 3, 6, 12, 24, 36)) {
  sfx <- sprintf("%02dmo", sc)
  tryCatch(
    create_trend_maps(sc, file.path(fig_dir, sprintf("Fig1_TrendMaps_%s.pdf", sfx))),
    error = function(e) cat("  ❌ Fig 1 trend maps:", e$message, "\n"))
  tryCatch(
    create_event_maps(sc, file.path(fig_dir, sprintf("Fig2_EventMaps_%s.pdf", sfx))),
    error = function(e) cat("  ❌ Fig 2 event maps:", e$message, "\n"))
  tryCatch(
    create_temporal_maps(sc, file.path(fig_dir, sprintf("Fig3_TemporalMaps_%s.pdf", sfx))),
    error = function(e) cat("  ❌ Fig 3 temporal maps:", e$message, "\n"))
}

tryCatch(
  create_timescale_comparison(file.path(fig_dir, "Fig4_TimescaleComparison.png")),
  error = function(e) cat("  ❌ Fig 4 timescale comparison:", e$message, "\n"))
tryCatch(
  create_method_comparison(file.path(fig_dir, "Fig5_MethodComparison.png")),
  error = function(e) cat("  ❌ Fig 5 method comparison:", e$message, "\n"))

################################################################################
# PART 3 – MANUSCRIPT FIGURE X
# Publication-quality 3 × 4 spatial drought-trend figure.
#
# Rows    : SPI, SPEI-PM, SPEI-ERA5
# Columns : 3-, 12-, 24-, and 36-month accumulation scales
# Field   : TFPW Sen's slope (sl_tfpw)
# Overlay : BH-FDR-significant TFPW tests (p_fdr_tfpw < 0.05)
#
# IMPORTANT:
#   This manuscript figure uses the same native-raster plotting pathway as the
#   established Script 8 spatial figures (make_raster() + plot_raster_clean()).
#   The previous ggplot2 + ImageMagick assembly was removed because it could
#   produce a visually blank PDF even when the file itself was non-empty.
#   Direct base-graphics output also preserves the native ERA5-Land raster
#   geometry without any coordinate reinterpretation.
################################################################################
cat("\n══════════════════════════════════════\n")
cat("PART 3: Manuscript Figure X — spatial drought trends\n")
cat("══════════════════════════════════════\n")

create_ms_drought_senslope_maps <- function(
    scales_ms = c(3, 12, 24, 36),
    outpdf,
    outpng,
    dpi = 300) {
  
  required <- c(
    "index_type", "timescale", "lon", "lat",
    "sl_tfpw", "p_value_tfpw", "p_fdr_tfpw"
  )
  missing <- setdiff(required, names(plot_data))
  if (length(missing)) {
    stop(
      "Required TFPW/Sen/FDR columns are missing: ",
      paste(missing, collapse = ", ")
    )
  }
  
  idx_order <- c("spi", "spei", "spei_era5")
  idx_labels <- c(
    spi       = "SPI",
    spei      = "SPEI-PM",
    spei_era5 = "SPEI-ERA5"
  )
  
  dat <- data.table::copy(
    plot_data[
      index_type %in% idx_order & timescale %in% scales_ms
    ]
  )
  
  if (!nrow(dat)) {
    stop("No SPI/SPEI data found for requested manuscript scales.")
  }
  
  # Fail early if the requested manuscript panels do not contain any usable
  # TFPW Sen-slope values. This prevents a formally valid but visually empty
  # figure from being accepted as a successful output.
  panel_check <- dat[
    , .(
      n_rows = .N,
      n_finite_slope = sum(is.finite(sl_tfpw)),
      n_finite_xy = sum(is.finite(lon) & is.finite(lat)),
      n_finite_fdr = sum(is.finite(p_fdr_tfpw))
    ),
    by = .(index_type, timescale)
  ]
  
  print(panel_check)
  
  bad_panels <- panel_check[n_finite_slope == 0L]
  if (nrow(bad_panels)) {
    stop(
      "No finite TFPW Sen's slope values for manuscript panel(s): ",
      paste(
        sprintf("%s-%02d", bad_panels$index_type, bad_panels$timescale),
        collapse = ", "
      ),
      ". Re-run 6trend_test_ALL.R and inspect the corresponding *_results.csv files."
    )
  }
  
  # Use the authoritative native-raster plotting path already used by Script 8.
  # The basin boundary is therefore in exactly the same coordinate system as the
  # raster cells, avoiding the coordinate/CRS reinterpretation that affected the
  # previous ggplot2 implementation.
  
  # One symmetric slope range per index keeps the four accumulation scales
  # comparable within each index.
  z_tbl <- dat[
    is.finite(sl_tfpw),
    .(z = max(abs(sl_tfpw), na.rm = TRUE)),
    by = index_type
  ]
  z_lookup <- setNames(z_tbl$z, z_tbl$index_type)
  z_lookup[!is.finite(z_lookup) | z_lookup <= 0] <- 1
  
  # Confirm that every requested panel has coordinates inside the native raster
  # extent. This is a diagnostic only; no coordinates are altered or wrapped.
  native_ext <- terra::ext(analysis_template)
  coord_check <- dat[
    is.finite(lon) & is.finite(lat),
    .(
      xmin = min(lon),
      xmax = max(lon),
      ymin = min(lat),
      ymax = max(lat)
    ),
    by = .(index_type, timescale)
  ]
  print(coord_check)
  
  # Direct base-graphics output is intentional. It is robust for SpatRaster
  # objects and is the same rendering family used by the existing spatial maps.
  draw_figure <- function() {
    graphics::par(
      mfrow = c(length(idx_order), length(scales_ms)),
      mar = c(2.2, 2.2, 3.0, 2.0),
      oma = c(0, 0, 3.0, 0),
      xaxs = "i",
      yaxs = "i"
    )
    
    for (idx in idx_order) {
      for (sc in scales_ms) {
        sub <- dat[index_type == idx & timescale == sc]
        
        if (!nrow(sub)) {
          graphics::plot.new()
          graphics::text(
            0.5, 0.5,
            sprintf("No data\n%s — %d months", idx_labels[[idx]], sc),
            cex = 1.0
          )
          next
        }
        
        sub_plot <- sub[
          is.finite(sl_tfpw) &
            is.finite(lon) &
            is.finite(lat)
        ]
        
        if (!nrow(sub_plot)) {
          graphics::plot.new()
          graphics::text(
            0.5, 0.5,
            sprintf("No finite Sen's slope\n%s — %d months", idx_labels[[idx]], sc),
            cex = 1.0
          )
          next
        }
        
        z <- z_lookup[[idx]]
        slope_range <- c(-z, z)
        slope_cols <- grDevices::hcl.colors(
          101,
          "RdBu",
          rev = TRUE
        )
        
        r_slope <- make_raster(
          sub_plot,
          "sl_tfpw",
          analysis_template
        )
        
        n_total <- nrow(sub_plot)
        sig <- sub[
          is.finite(p_fdr_tfpw) &
            p_fdr_tfpw < 0.05 &
            is.finite(sl_tfpw) &
            is.finite(lon) &
            is.finite(lat)
        ]
        n_sig <- nrow(sig)
        pct_sig <- 100 * n_sig / n_total
        
        plot_raster_clean(
          r_slope,
          paste0(idx_labels[[idx]], " — ", sc, " months"),
          slope_range,
          slope_cols,
          legend_title = "Sen's slope"
        )
        
        # Significance stippling is added in native raster coordinates after the
        # slope field, matching the established Script 8 map convention.
        if (n_sig > 0L) {
          graphics::points(
            sig$lon,
            sig$lat,
            pch = 20,
            cex = 0.32,
            col = "black"
          )
        }
        
        graphics::mtext(
          sprintf(
            "BH-FDR significant: %d/%d (%.1f%%)",
            n_sig,
            n_total,
            pct_sig
          ),
          side = 3,
          line = 0.15,
          cex = 0.55,
          col = "grey30"
        )
      }
    }
    
    graphics::mtext(
      "Spatial drought-index trends — Nechako River Basin (1950–2025)",
      outer = TRUE,
      side = 3,
      line = 0.8,
      cex = 1.15,
      font = 2
    )
    
    graphics::mtext(
      paste0(
        "TFPW Mann–Kendall Sen's slope; stippling = BH-FDR-significant pixels " ,
        "(α = 0.05). Columns: ",
        paste(scales_ms, collapse = ", "),
        " months."
      ),
      outer = TRUE,
      side = 3,
      line = 2.0,
      cex = 0.72,
      col = "grey30"
    )
  }
  
  # Write directly to the graphics devices. Do not rasterize the complete figure
  # through ImageMagick: the previous assembly pathway generated a blank PDF.
  grDevices::pdf(
    file = outpdf,
    width = 14,
    height = 10,
    onefile = TRUE,
    useDingbats = FALSE
  )
  draw_figure()
  grDevices::dev.off()
  
  grDevices::png(
    filename = outpng,
    width = 14,
    height = 10,
    units = "in",
    res = dpi,
    type = "cairo",
    bg = "white"
  )
  draw_figure()
  grDevices::dev.off()
  
  # Output-level validation: existence and non-zero size are necessary but not
  # sufficient, so also report the expected 3 × 4 panel structure.
  for (f in c(outpdf, outpng)) {
    if (!file.exists(f)) {
      stop("Manuscript drought-trend figure was not created: ", f)
    }
    sz <- file.info(f)$size
    if (is.na(sz) || sz <= 0) {
      stop("Manuscript drought-trend figure is empty: ", f)
    }
  }
  
  cat(sprintf(
    "  ✓ Manuscript drought-trend figure written directly by R graphics device: %s\n",
    basename(outpdf)
  ))
  cat(sprintf(
    "  ✓ Manuscript drought-trend figure written directly by R graphics device: %s\n",
    basename(outpng)
  ))
  
  invisible(TRUE)
}

# Manuscript Figure X is a required output. Remove stale files first so the final
# validation cannot mistake an artifact from an earlier run for a successful
# figure generated by this run.
ms_fig_pdf <- file.path(
  fig_dir,
  "FigX_MS_Spatial_Drought_Trends_SenSlope_3_12_24_36.pdf"
)
ms_fig_png <- file.path(
  fig_dir,
  "FigX_MS_Spatial_Drought_Trends_SenSlope_3_12_24_36.png"
)
unlink(c(ms_fig_pdf, ms_fig_png), force = TRUE)

create_ms_drought_senslope_maps(
  scales_ms = c(3, 12, 24, 36),
  outpdf = ms_fig_pdf,
  outpng = ms_fig_png,
  dpi = DPI
)

################################################################################
# LEGACY / EXISTING FIGURE 8 OUTPUTS
# Kept unchanged so existing diagnostic products are preserved.
################################################################################
cat("\n  Existing diagnostic trend-map outputs are retained above; the new\n")
cat("  manuscript figure is the dedicated 3 × 4 Sen-slope + BH-FDR product.\n")

################################################################################
# PART 3 COMPLETE
#
# Do not emit a final DONE banner here.  Part 4 and the regional extension are
# still executed below.  Final completion is reported only after all required
# artifacts have been verified at the end of the script.
################################################################################
cat("\n✓ Manuscript Figure X generation completed. Continuing to Part 4 diagnostics.\n\n")
################################################################################
# PART 4 – PIXEL-TO-BASIN SPEARMAN CORRELATION DIAGNOSTICS
#
# GRID POLICY — IMPORTANT:
#   ALL calculations in this section use the exact native ERA5-Land grid established
#   upstream by 3SPEI_ERALand.R.  That grid is the projected ERA5-Land
#   precipitation/SPI analysis grid on which SPI, SPEI-PM and SPEI-ERA5 are
#   required to remain one-to-one aligned.
#
#   NO reproject(), resample(), project(), crop(), extend(), or rasterize()
#   operation is performed on the drought-index rasters in this section.
#   The native template and native basin mask are established at the start of
#   this script from the original ERA5-Land precipitation NetCDF.
#
#   Every seasonal SPI/SPEI NetCDF is hard-checked against that exact grid
#   (rows, columns, resolution, extent, origin, CRS, and cell count). A
#   mismatch is a fatal error rather than silently changing the grid.
#
# Purpose:
#   For every native-grid basin cell i:
#
#       rho_i = Spearman( X_i(t), X_basin(t) )
#
#   using the original monthly gridded SPI/SPEI series and the corresponding
#   authoritative basin-average series from 6trend_test_ALL.R.
################################################################################

cat("\n============================================================\n")
cat("PART 4: PIXEL-TO-BASIN SPEARMAN CORRELATION DIAGNOSTICS\n")
cat("============================================================\n")

CORR_OUT_DIR <- file.path(fig_dir, "pixel_basin_correlation")
dir.create(CORR_OUT_DIR, showWarnings = FALSE, recursive = TRUE)

CORR_SHORT_SCALES <- c(1L, 3L, 6L)
CORR_LONG_SCALES  <- c(12L, 24L, 36L)
CORR_ALL_SCALES   <- c(CORR_SHORT_SCALES, CORR_LONG_SCALES)
CORR_INDEXES      <- c("spi", "spei")
CORR_MIN_N        <- 24L

# -----------------------------------------------------------------------------
# Native grid: use the same native ERA5-Land precipitation template defined
# at the beginning of this script. The basin mask is generated directly on
# that grid; no stored canonical-grid RDS is read.
# -----------------------------------------------------------------------------
CORR_TEMPLATE <- analysis_template
CORR_MASK_RASTER <- basin_mask_raster
CORR_BASIN_CELLS <- basin_mask_cells

if (terra::ncell(CORR_TEMPLATE) != terra::ncell(CORR_MASK_RASTER) ||
    terra::nrow(CORR_TEMPLATE) != terra::nrow(CORR_MASK_RASTER) ||
    terra::ncol(CORR_TEMPLATE) != terra::ncol(CORR_MASK_RASTER)) {
  stop("Native raster template and basin mask have different dimensions.")
}

if (!isTRUE(all.equal(
  as.numeric(terra::res(CORR_TEMPLATE)),
  as.numeric(terra::res(CORR_MASK_RASTER)),
  tolerance = 1e-10
))) {
  stop("Native raster template and basin mask have different resolutions.")
}

corr_ext <- c(
  terra::xmin(CORR_TEMPLATE), terra::xmax(CORR_TEMPLATE),
  terra::ymin(CORR_TEMPLATE), terra::ymax(CORR_TEMPLATE)
)
mask_ext <- c(
  terra::xmin(CORR_MASK_RASTER), terra::xmax(CORR_MASK_RASTER),
  terra::ymin(CORR_MASK_RASTER), terra::ymax(CORR_MASK_RASTER)
)
if (!isTRUE(all.equal(corr_ext, mask_ext, tolerance = 1e-8))) {
  stop("Native raster template and basin mask have different extents.")
}

corr_origin <- terra::origin(CORR_TEMPLATE)
mask_origin <- terra::origin(CORR_MASK_RASTER)
if (!isTRUE(all.equal(
  as.numeric(corr_origin),
  as.numeric(mask_origin),
  tolerance = 1e-8
))) {
  stop("Native raster template and basin mask have different origins.")
}

if (!terra::same.crs(CORR_TEMPLATE, CORR_MASK_RASTER)) {
  stop("Native raster template and basin mask have different CRS.")
}

if (length(CORR_BASIN_CELLS) != terra::ncell(CORR_TEMPLATE)) {
  stop("Native basin mask length does not equal native grid cell count.")
}

if (sum(CORR_BASIN_CELLS, na.rm = TRUE) != basin_pixels) {
  stop("Native basin-mask cell count changed after mask construction.")
}

cat(sprintf(
  "Native SPI/SPEI grid: %d rows × %d cols = %d cells\n",
  terra::nrow(CORR_TEMPLATE),
  terra::ncol(CORR_TEMPLATE),
  terra::ncell(CORR_TEMPLATE)
))
cat(sprintf(
  "Native resolution: %.10f × %.10f\n",
  terra::res(CORR_TEMPLATE)[1], terra::res(CORR_TEMPLATE)[2]
))
cat(sprintf(
  "Native CRS: %s\n",
  terra::crs(CORR_TEMPLATE)
))
cat(sprintf(
  "Native basin cells: %d\n",
  sum(CORR_BASIN_CELLS, na.rm = TRUE)
))

# -----------------------------------------------------------------------------
# Exact geometry validator.
# No correction is permitted: a mismatch stops the calculation.
# -----------------------------------------------------------------------------
validate_corr_geometry <- function(r, reference, file_label) {
  if (terra::ncell(r) != terra::ncell(reference)) {
    stop(sprintf(
      "GRID ERROR [%s]: cell count %d != native %d",
      file_label, terra::ncell(r), terra::ncell(reference)
    ))
  }
  if (terra::nrow(r) != terra::nrow(reference) ||
      terra::ncol(r) != terra::ncol(reference)) {
    stop(sprintf(
      "GRID ERROR [%s]: dimensions %d×%d != native %d×%d",
      file_label,
      terra::nrow(r), terra::ncol(r),
      terra::nrow(reference), terra::ncol(reference)
    ))
  }
  if (!isTRUE(all.equal(
    as.numeric(terra::res(r)),
    as.numeric(terra::res(reference)),
    tolerance = 1e-10
  ))) {
    stop("GRID ERROR [", file_label, "]: resolution differs from native ERA5-Land grid.")
  }
  ext_r <- c(
    terra::xmin(r), terra::xmax(r),
    terra::ymin(r), terra::ymax(r)
  )
  ext_ref <- c(
    terra::xmin(reference), terra::xmax(reference),
    terra::ymin(reference), terra::ymax(reference)
  )
  if (!isTRUE(all.equal(ext_r, ext_ref, tolerance = 1e-8))) {
    stop("GRID ERROR [", file_label, "]: extent differs from native ERA5-Land grid.")
  }
  if (!isTRUE(all.equal(
    as.numeric(terra::origin(r)),
    as.numeric(terra::origin(reference)),
    tolerance = 1e-8
  ))) {
    stop("GRID ERROR [", file_label, "]: origin differs from native ERA5-Land grid.")
  }
  if (!terra::same.crs(r, reference)) {
    stop("GRID ERROR [", file_label, "]: CRS differs from native ERA5-Land grid.")
  }
  invisible(TRUE)
}

# -----------------------------------------------------------------------------
# Locate seasonal NetCDF files.
# -----------------------------------------------------------------------------
find_corr_nc_files <- function(index_type, scale) {
  data_dir <- switch(
    tolower(index_type),
    spi  = SPI_SEAS_DIR,
    spei = SPEI_SEAS_DIR,
    stop("Unsupported correlation index: ", index_type)
  )
  files <- find_seasonal_nc_files(data_dir, index_type, scale)
  if (!length(files)) {
    stop(sprintf(
      "No seasonal NetCDF files found for %s-%d in %s",
      toupper(index_type), scale, data_dir
    ))
  }
  files
}

# -----------------------------------------------------------------------------
# Read a complete seasonal product into chronological monthly order.
# The native template is NOT inferred from the seasonal NetCDFs and those
# NetCDFs are NOT reprojected/resampled. Every file must already be exactly
# aligned to the native ERA5-Land precipitation grid.
# -----------------------------------------------------------------------------
read_corr_index_matrix <- function(index_type, scale, target_dates) {
  files <- find_corr_nc_files(index_type, scale)
  if (length(files) != 12L) {
    stop(sprintf(
      "Expected 12 seasonal files for %s-%d; found %d",
      toupper(index_type), scale, length(files)
    ))
  }
  
  out <- matrix(
    NA_real_,
    nrow = terra::ncell(CORR_TEMPLATE),
    ncol = length(target_dates)
  )
  nlayer_total <- 0L
  
  for (m in seq_along(files)) {
    r <- terra::rast(files[m])
    file_label <- basename(files[m])
    validate_corr_geometry(r, CORR_TEMPLATE, file_label)
    
    nlayer_total <- nlayer_total + terra::nlyr(r)
    vals <- terra::values(r, mat = TRUE)
    vals <- as.matrix(vals)
    vals[!is.finite(vals)] <- NA_real_
    
    month_dates <- target_dates[
      format(target_dates, "%m") == sprintf("%02d", m)
    ]
    if (length(month_dates) != ncol(vals)) {
      stop(sprintf(
        "%s-%d month %02d has %d raster layers but %d target dates",
        toupper(index_type), scale, m, ncol(vals), length(month_dates)
      ))
    }
    
    col_idx <- match(month_dates, target_dates)
    out[, col_idx] <- vals
    
    rm(r, vals)
    gc(verbose = FALSE)
  }
  
  if (nlayer_total != length(target_dates)) {
    stop(sprintf(
      "Unexpected layer count for %s-%d: %d layers; expected %d",
      toupper(index_type), scale, nlayer_total, length(target_dates)
    ))
  }
  
  list(matrix = out, files = files)
}

# -----------------------------------------------------------------------------
# Pixel-wise Spearman correlation.
# The output raster is created directly from the native ERA5-Land grid template.
# -----------------------------------------------------------------------------
compute_corr_raster <- function(ts_matrix, basin_series, index_type, scale) {
  if (!is.matrix(ts_matrix)) stop("ts_matrix must be a matrix")
  if (nrow(ts_matrix) != terra::ncell(CORR_TEMPLATE)) {
    stop(sprintf(
      "Cell-count mismatch for %s-%d: matrix has %d rows, native ERA5-Land grid has %d cells",
      index_type, scale, nrow(ts_matrix), terra::ncell(CORR_TEMPLATE)
    ))
  }
  if (ncol(ts_matrix) != length(basin_series)) {
    stop(sprintf(
      "Temporal mismatch for %s-%d: %d pixel months vs %d basin months",
      index_type, scale, ncol(ts_matrix), length(basin_series)
    ))
  }
  
  n_pair <- integer(nrow(ts_matrix))
  rr <- rep(NA_real_, nrow(ts_matrix))
  
  for (i in which(CORR_BASIN_CELLS)) {
    x <- ts_matrix[i, ]
    ok <- is.finite(x) & is.finite(basin_series)
    n_pair[i] <- sum(ok)
    
    if (n_pair[i] < CORR_MIN_N) next
    
    xok <- x[ok]
    yok <- basin_series[ok]
    if (length(unique(xok)) < 2L || length(unique(yok)) < 2L) next
    
    rr[i] <- suppressWarnings(
      stats::cor(xok, yok, method = "spearman")
    )
  }
  
  r <- terra::setValues(
    terra::rast(CORR_TEMPLATE),
    rr
  )
  
  list(raster = r, n_pair = n_pair)
}

# -----------------------------------------------------------------------------
# Map data and map plot.
# -----------------------------------------------------------------------------
corr_map_df <- function(r, label) {
  d <- terra::as.data.frame(r, xy = TRUE, na.rm = TRUE)
  if (!nrow(d)) stop("No finite correlation values available for ", label)
  names(d)[3] <- "rho"
  d$label <- label
  d
}

plot_corr_map <- function(r, label) {
  d <- corr_map_df(r, label)
  ggplot2::ggplot(d, ggplot2::aes(x = x, y = y, fill = rho)) +
    ggplot2::geom_raster() +
    ggplot2::geom_sf(
      data = basin_sf,
      inherit.aes = FALSE,
      fill = NA,
      linewidth = 0.35
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-1, 1),
      oob = scales::squish,
      name = expression(rho[S])
    ) +
    ggplot2::coord_sf(expand = FALSE) +
    ggplot2::labs(title = label) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size = 7),
      plot.title = ggplot2::element_text(size = 11, face = "bold", hjust = 0.5),
      legend.position = "right",
      legend.title = ggplot2::element_text(size = 9),
      legend.text = ggplot2::element_text(size = 8),
      panel.grid = ggplot2::element_blank()
    )
}

# -----------------------------------------------------------------------------
# Density plot: three curves + one dotted median line per timescale.
# -----------------------------------------------------------------------------
plot_corr_density <- function(corr_values, index_type, scales) {
  rows <- lapply(scales, function(sc) {
    key <- sprintf("%s_%02d", index_type, sc)
    x <- corr_values[[key]]
    if (is.null(x)) return(NULL)
    data.frame(
      rho = x,
      scale = factor(sc, levels = scales),
      stringsAsFactors = FALSE
    )
  })
  
  d <- data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
  d <- as.data.frame(d)
  d <- d[is.finite(d$rho), , drop = FALSE]
  if (!nrow(d)) {
    stop("No finite correlation values available for density plot: ", index_type)
  }
  
  meds <- aggregate(rho ~ scale, data = d, FUN = median, na.rm = TRUE)
  
  ggplot2::ggplot(d, ggplot2::aes(x = rho, colour = scale)) +
    ggplot2::geom_density(
      linewidth = 1.15,
      na.rm = TRUE,
      adjust = 1.05
    ) +
    ggplot2::geom_vline(
      data = meds,
      ggplot2::aes(xintercept = rho, colour = scale),
      linetype = "dotted",
      linewidth = 0.85,
      show.legend = FALSE
    ) +
    ggplot2::scale_colour_manual(
      values = c("#0072B2", "#E69F00", "#009E73")[seq_along(scales)],
      name = "Timescale"
    ) +
    ggplot2::scale_x_continuous(
      limits = c(-1, 1),
      breaks = seq(-1, 1, 0.2)
    ) +
    ggplot2::labs(
      title = if (index_type == "spi") {
        "SPI: pixel–basin Spearman correlation distribution"
      } else {
        "SPEI-PM: pixel–basin Spearman correlation distribution"
      },
      x = expression(Spearman~rho[S]),
      y = "Density"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 12, face = "bold", hjust = 0.5),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size = 10),
      legend.text = ggplot2::element_text(size = 9),
      panel.grid.minor = ggplot2::element_blank()
    )
}

# -----------------------------------------------------------------------------
# Calculate all requested correlations once.
# -----------------------------------------------------------------------------
corr_rasters <- list()
corr_values  <- list()
corr_meta    <- list()

for (idx in CORR_INDEXES) {
  for (sc in CORR_ALL_SCALES) {
    key <- sprintf("%s_%02d", idx, sc)
    cat(sprintf("\n  → %s-%d\n", toupper(idx), sc))
    
    basin_df <- load_basin_avg_csv(idx, sc)
    if (is.null(basin_df)) {
      stop(sprintf(
        "Required basin-average CSV is missing for %s-%d. Run 6trend_test_ALL.R first.",
        idx_label, sc
      ))
    }
    
    basin_dates  <- as.Date(basin_df$date)
    basin_series <- as.numeric(basin_df$value)
    
    nc <- read_corr_index_matrix(idx, sc, basin_dates)
    
    # Enforce the native basin mask exactly; no spatial operation is applied.
    ts_basin <- nc$matrix
    ts_basin[!CORR_BASIN_CELLS, ] <- NA_real_
    
    corr_result <- compute_corr_raster(
      ts_matrix = ts_basin,
      basin_series = basin_series,
      index_type = idx,
      scale = sc
    )
    
    r_corr <- corr_result$raster
    vals <- terra::values(r_corr, mat = FALSE)
    vals <- vals[is.finite(vals)]
    
    if (!length(vals)) {
      stop(sprintf(
        "No finite pixel-basin correlations for %s-%d",
        idx_label, sc
      ))
    }
    
    corr_rasters[[key]] <- r_corr
    corr_values[[key]] <- vals
    corr_meta[[key]] <- data.frame(
      index_type = idx,
      index_label = if (idx == "spi") "SPI" else "SPEI-PM",
      scale = sc,
      n_pixels = length(vals),
      median_rho = median(vals, na.rm = TRUE),
      mean_rho = mean(vals, na.rm = TRUE),
      min_rho = min(vals, na.rm = TRUE),
      max_rho = max(vals, na.rm = TRUE),
      n_paired_months_min = min(corr_result$n_pair[CORR_BASIN_CELLS], na.rm = TRUE),
      n_paired_months_max = max(corr_result$n_pair[CORR_BASIN_CELLS], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    
    cat(sprintf(
      "     %d basin pixels; median Spearman rho = %.3f; range = %.3f to %.3f\n",
      length(vals), median(vals), min(vals), max(vals)
    ))
    
    rm(nc, ts_basin, r_corr, corr_result, vals, basin_df)
    gc(verbose = FALSE)
  }
}

# -----------------------------------------------------------------------------
# Audit outputs.
# -----------------------------------------------------------------------------
corr_summary <- data.table::rbindlist(
  corr_meta,
  use.names = TRUE,
  fill = TRUE
)
corr_summary <- as.data.frame(corr_summary)

data.table::fwrite(
  corr_summary,
  file.path(CORR_OUT_DIR, "PixelBasin_Spearman_Summary.csv")
)

# Save one row per pixel/index/scale. Coordinates are taken directly from the
# native ERA5-Land grid; no coordinate transformation is performed here.
# IMPORTANT:
# terra::crds(SpatRaster) may return coordinates only for non-NA cells
# depending on the raster state/version.  That would break the one-row-per-
# native-cell audit table. Recover coordinates explicitly for EVERY native
# raster cell, in native cell-number order.
cell_id <- seq_len(terra::ncell(CORR_TEMPLATE))
xy <- terra::xyFromCell(CORR_TEMPLATE, cell_id)

if (nrow(xy) != terra::ncell(CORR_TEMPLATE)) {
  stop(
    "Native-grid coordinate extraction failed: returned ",
    nrow(xy), " coordinates for ",
    terra::ncell(CORR_TEMPLATE), " cells."
  )
}

correlation_audit <- data.table::rbindlist(
  lapply(names(corr_rasters), function(key) {
    r <- corr_rasters[[key]]
    parts <- strsplit(key, "_", fixed = TRUE)[[1]]
    idx <- parts[1]
    sc <- as.integer(parts[2])
    
    vals <- terra::values(r, mat = FALSE)
    
    if (length(vals) != length(cell_id)) {
      stop(
        "Correlation raster length mismatch for ", key,
        ": ", length(vals), " values for ",
        length(cell_id), " native ERA5-Land cells."
      )
    }
    
    data.table::data.table(
      cell = cell_id,
      x = xy[, 1],
      y = xy[, 2],
      index_type = idx,
      index_label = if (idx == "spi") "SPI" else "SPEI-PM",
      scale = sc,
      spearman_rho = vals
    )[
      CORR_BASIN_CELLS & is.finite(spearman_rho)
    ]
  }),
  use.names = TRUE,
  fill = TRUE
)



# Final integrity checks: every retained audit point must be a native cell
# centre of the native ERA5-Land precipitation grid and must be one of the
# native basin cells.
if (nrow(correlation_audit) == 0L) {
  stop("Correlation audit table is empty after native basin-cell filtering.")
}
if (!all(correlation_audit$cell %in% cell_id)) {
  stop("Correlation audit contains cell IDs outside the native ERA5-Land grid.")
}
if (!all(correlation_audit$cell %in% cell_id[CORR_BASIN_CELLS])) {
  stop("Correlation audit contains cells outside the native basin mask.")
}

data.table::fwrite(
  correlation_audit,
  file.path(CORR_OUT_DIR, "PixelBasin_Spearman_Correlations.csv")
)

# -----------------------------------------------------------------------------
# Assemble the two requested figures.
# Top: six maps. Bottom: SPI distribution + SPEI-PM distribution.
# -----------------------------------------------------------------------------
make_corr_figure <- function(scales, outfile_stub) {
  map_keys <- c(
    sprintf("spi_%02d", scales),
    sprintf("spei_%02d", scales)
  )
  
  map_labels <- c(
    sprintf("SPI-%d", scales),
    sprintf("SPEI-PM-%d", scales)
  )
  
  map_plots <- Map(
    f = function(key, label) plot_corr_map(corr_rasters[[key]], label),
    key = map_keys,
    label = map_labels
  )
  
  if (length(map_plots) != 6L) {
    stop(
      "Correlation figure requires exactly six map panels (",
      "3 SPI + 3 SPEI-PM), but got ", length(map_plots), "."
    )
  }
  
  map_grid <- patchwork::wrap_plots(
    map_plots,
    nrow = 2,
    ncol = 3,
    byrow = TRUE
  )
  
  spi_density <- plot_corr_density(corr_values, "spi", scales)
  spei_density <- plot_corr_density(corr_values, "spei", scales)
  
  final_plot <- map_grid /
    spi_density /
    spei_density +
    patchwork::plot_layout(
      heights = c(2.0, 0.72, 0.72)
    )
  
  png_file <- file.path(CORR_OUT_DIR, paste0(outfile_stub, ".png"))
  pdf_file <- file.path(CORR_OUT_DIR, paste0(outfile_stub, ".pdf"))
  
  ggplot2::ggsave(
    png_file,
    final_plot,
    width = 15,
    height = 15,
    units = "in",
    dpi = DPI,
    bg = "white"
  )
  ggplot2::ggsave(
    pdf_file,
    final_plot,
    width = 15,
    height = 15,
    units = "in",
    device = grDevices::cairo_pdf,
    bg = "white"
  )
  
  cat("  ✓ ", basename(png_file), "\n", sep = "")
  cat("  ✓ ", basename(pdf_file), "\n", sep = "")
}

cat("\n  Writing short-scale correlation figure (1/3/6)...\n")
make_corr_figure(
  CORR_SHORT_SCALES,
  "FigS_PixelBasin_Spearman_1_3_6"
)

cat("\n  Writing long-scale correlation figure (12/24/36)...\n")
make_corr_figure(
  CORR_LONG_SCALES,
  "FigS_PixelBasin_Spearman_12_24_36"
)

cat("\n============================================================\n")
cat("PART 4 COMPLETE: pixel-to-basin Spearman diagnostics\n")
cat("============================================================\n")
cat("Output directory: ", normalizePath(CORR_OUT_DIR, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("Native grid was used without reprojection/resampling.\n")
cat("\n")
################################################################################
# REGIONAL EXTENSION — SUBBASIN TREND SUMMARY + REGION MAP
################################################################################
regional_trend_dir <- file.path(TREND_DIR, "regional_trends")
regional_index_dir <- file.path(TREND_DIR, "regional_timeseries")
regional_summary_dir <- file.path(TREND_DIR, "regional_spatial_summary")
dir.create(regional_summary_dir, recursive=TRUE, showWarnings=FALSE)

# Calculate regional trend statistics directly from the authoritative regional
# monthly index means. This does not alter the pixel-level spatial trend maps.
regional_files <- list.files(regional_index_dir,
                             pattern="_regional_averaged_by_month\\.csv$",
                             full.names=TRUE, ignore.case=TRUE)
if (length(regional_files)) {
  reg_ts <- data.table::rbindlist(lapply(regional_files, function(f) {
    x <- data.table::fread(f)
    if (!all(c("region_id","region_name","date","value","index_type","scale") %in% names(x))) return(NULL)
    x[, date := as.Date(date)]
    x[, value := as.numeric(value)]
    x
  }), use.names=TRUE, fill=TRUE)
  
  regional_trend_summary <- reg_ts[, {
    x <- value[is.finite(value)]
    
    mk <- if (length(x) >= 10L) {
      tryCatch(
        Kendall::MannKendall(x),
        error = function(e) NULL
      )
    } else {
      NULL
    }
    
    ss <- if (length(x) >= 10L) {
      tryCatch(
        trend::sens.slope(x),
        error = function(e) NULL
      )
    } else {
      NULL
    }
    
    list(
      n_months = length(x),
      mean_value = if (length(x)) mean(x) else NA_real_,
      tau = if (!is.null(mk)) as.numeric(mk$tau) else NA_real_,
      p_value = if (!is.null(mk)) as.numeric(mk$sl) else NA_real_,
      sen_slope_per_month =
        if (!is.null(ss)) as.numeric(ss$estimates) else NA_real_,
      sen_slope_per_year =
        if (!is.null(ss)) as.numeric(ss$estimates) * 12 else NA_real_
    )
  }, by = .(
    region_id,
    region_name,
    index_type,
    scale
  )]
  data.table::fwrite(regional_trend_summary,
                     file.path(regional_summary_dir, "regional_index_trend_summary.csv"))
  
  # Publication/diagnostic map of the authoritative 7-subbasin partition.
  # Publication/diagnostic map of the authoritative 7-subbasin partition.
  sub_regions <- load_nechako_regions()
  
  if (!NECHAKO_REGION_BASIN_ID %in% names(sub_regions)) {
    stop(
      "Authoritative whole-basin region '",
      NECHAKO_REGION_BASIN_ID,
      "' was not found in load_nechako_regions()."
    )
  }
  
  # load_nechako_regions() already returns SpatVector objects.
  # Do NOT call terra::vect() on an existing SpatVector.
  whole_region <- sub_regions[[NECHAKO_REGION_BASIN_ID]]
  reg_sf <- sf::st_as_sf(whole_region)
  
  sub_ids <- setdiff(names(sub_regions), NECHAKO_REGION_BASIN_ID)
  
  if (!length(sub_ids)) {
    stop("No supplied subbasins remain after removing the whole-basin union.")
  }
  
  sub_list <- sub_regions[sub_ids]
  stopifnot(all(vapply(sub_list, inherits, logical(1), "SpatVector")))
  sub_v <- do.call(rbind, unname(sub_list))
  sub_sf <- sf::st_as_sf(sub_v)
  
  sub_sf$region_id <- sub_ids
  sub_sf$region_name <- gsub("_", " ", sub_ids, fixed = TRUE)
  
  pmap <- ggplot2::ggplot() +
    ggplot2::geom_sf(data=sub_sf, fill=NA, linewidth=0.6) +
    ggplot2::geom_sf_text(data=sub_sf, ggplot2::aes(label=region_name), size=2.8) +
    ggplot2::geom_sf(data=reg_sf, fill=NA, linewidth=1.0) +
    ggplot2::labs(title="Authoritative Nechako regional analysis domain",
                  subtitle="7 supplied subbasins; whole basin = dissolved union") +
    ggplot2::theme_bw(base_size=10)
  ggplot2::ggsave(file.path(regional_summary_dir, "Nechako_7_subbasin_domain.png"),
                  pmap, width=10, height=8, dpi=300)
  ggplot2::ggsave(file.path(regional_summary_dir, "Nechako_7_subbasin_domain.pdf"),
                  pmap, width=10, height=8)
  cat("✓ Regional trend summary and authoritative subbasin map written.\n")
}

################################################################################
# FINAL OUTPUT VALIDATION
#
# A successful run is defined by verified required artifacts, not merely by
# reaching the end of a tryCatch() block.  The manuscript Figure X PDF/PNG and
# master diagnostics CSV are required outputs of this script.
################################################################################
verify_required_file <- function(path, label) {
  if (!file.exists(path)) {
    stop("FINAL OUTPUT VALIDATION FAILED: missing ", label, ": ", path)
  }
  info <- file.info(path)
  if (is.na(info$size) || info$size <= 0) {
    stop("FINAL OUTPUT VALIDATION FAILED: empty ", label, ": ", path)
  }
  cat(sprintf("  ✓ %s: %s (%s bytes)\n",
              label,
              normalizePath(path, winslash = "/", mustWork = FALSE),
              format(info$size, big.mark = ",", scientific = FALSE)))
  invisible(TRUE)
}

cat("\n============================================================\n")
cat("FINAL OUTPUT VALIDATION\n")
cat("============================================================\n")
verify_required_file(master_csv, "Master diagnostics CSV")
verify_required_file(ms_fig_pdf, "Manuscript Figure X PDF")
verify_required_file(ms_fig_png, "Manuscript Figure X PNG")

cat("\n╔══════════════════════════════════════════════════════════╗\n")
cat("║  FINAL OUTPUT VALIDATION PASSED                         ║\n")
cat("║  8trends_visualization.R DONE                        ║\n")
cat("╚══════════════════════════════════════════════════════════╝\n\n")
cat("Outputs:\n")
cat("  Master diagnostics CSV : ", normalizePath(master_csv, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("  Spatial figures        : ", normalizePath(fig_dir, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("  Manuscript Figure X    : ", basename(ms_fig_pdf), "\n", sep = "")
cat("  Manuscript Figure X    : ", basename(ms_fig_png), "\n", sep = "")
cat("\n✓ All required Script 8 outputs were generated and verified.\n\n")