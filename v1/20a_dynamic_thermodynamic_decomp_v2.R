# ==============================================================================
#   10a_dynamic_thermodynamic_decomp_v2.R (observed-framework version)
# ==============================================================================
# *** LEGACY / SUPPLEMENTARY SCRIPT — NOT THE CURRENT MANUSCRIPT FRAMING ***
#
# This script is retained as a legacy/supplementary decomposition-style diagnostic.
# The current manuscript does not treat SPEI_PM - SPEI_Thw as an additive physical decomposition of drought severity. The current
# manuscript (Data and Methods 3e, "PET-Formulation Sensitivity") explicitly
# rejects this interpretation:
#
#   "Rather than interpreting the difference between the two SPEI
#    formulations as a decomposition of physical drought drivers, we use
#    their agreement and divergence to quantify PET-formulation
#    sensitivity ... subtraction of separately standardized SPEI series
#    does not constitute an additive physical decomposition of drought
#    forcing."
#
# 11nechako_drought_manuscript_v2.R accordingly gates all outputs that used
# this script's "thermodynamic fraction" framing behind
# RUN_LEGACY_OUTPUTS <- FALSE. The CURRENT Results section is produced by
# 12PET_formulation_sensitivity_v2.R, which recasts this exact PM/Thw
# comparison as ΔSPEI_PET = SPEI_PM - SPEI_THW (a formulation-sensitivity
# residual, not a "dynamic" component) and reports the four Methods-3e
# agreement metrics (Pearson r, RMSE, drought-state agreement, mean |ΔSPEI|)
# instead of a thermodynamic/dynamic percentage split.
#
# This script is RETAINED as a supplementary intermediate for two reasons:
#   1. It produces decomp_results/decomp_full_timeseries.csv, which
#      10b_PET_bias_nonstationarity.R depends on for its (Methods-3b-aligned,
#      non-decompositional) ΔPET bias/non-stationarity diagnostics.
#   2. Its "thermodynamic fraction" figures may still be useful as a
#      SUPPLEMENTARY / robustness appendix showing how the framing evolved,
#      but they should NOT be captioned or cited as primary Results figures
#      without first relabelling SPEI_Thw / SPEI_Dynamic away from
#      "thermodynamic" / "dynamic" language to avoid contradicting Methods 3e.
#
# ORIGINAL DESCRIPTION (kept for provenance; read the disclaimer above first):
# OBSERVED PM vs. THORNTHWAITE DROUGHT DECOMPOSITION
# Nechako River Basin, BC — 2022-2025 drought study
# Loads pre-computed SPEI from 3SPEI_ERALand_v2.R outputs
# - SPEI_PM  (spei_results_seasonal/)      = full observed drought
# - SPEI_Thw (spei_results_seasonal_thw/)  = Thornthwaite observed formulation
# - Difference = PET-formulation sensitivity residual (PM − Thw)
# Needs the SPEI output directories from 3SPEI_ERALand_v2.R
# ==============================================================================
library(terra)
library(lubridate)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggtext)
library(patchwork)
setwd("D:/Nechako_Drought/Nechako/")
# Source shared utilities — constants (DROUGHT_ONSET, SEVERE_THRESHOLD, EQUAL_AREA_CRS,
# SPEI_SEAS_DIR, etc.) and statistical functions must be available before any
# analysis code runs.  All other scripts source this file; 10a was the only one
# that did not, which would silently use wrong values if any shared constant was
# referenced in future edits.
source("DROUGHT_ANALYSIS_utils_v2.R")
out_dir <- "decomp_results"
if (!dir.exists(out_dir)) dir.create(out_dir)

# ------------------------------------------------------------------------------
# GLOBAL THRESHOLDS (used throughout all steps)
# ------------------------------------------------------------------------------
THR_EVENT <- -0.5   # event-identification threshold (Section 3.2 of manuscript)
THR_DROUGHT <- -0.5 # retained drought gate for legacy/seasonal diagnostics
THR_FTHM <- THR_DROUGHT  # backward-compatible alias; no physical fraction is computed

# ==============================================================================
#   STEP 1: LOAD PRE-COMPUTED SPEI RESULTS
# ==============================================================================
cat("\n===== STEP 1: LOAD PRE-COMPUTED SPEI RESULTS =====\n")
# Directories from 3SPEI_ERALand_v2.R
spei_pm_dir   <- "spei_results_seasonal"
spei_thw_dir  <- "spei_results_seasonal_thw"
if (!dir.exists(spei_pm_dir))  stop("SPEI_PM directory not found. Run 3SPEI_ERALand_v2.R first.")
if (!dir.exists(spei_thw_dir)) stop("SPEI_Thw directory not found. Run 3SPEI_ERALand_v2.R first.")
# Scales available from 3SPEI_ERALand_v2.R
decomp_scales <- c(1, 3, 6, 12, 24, 36)
# Helper function to load SPEI NetCDF files for one scale.
# Each monthly NetCDF is expected to contain 76 annual spatial layers, one
# for each year from 1950 through 2025. The loader reads all layers and
# interleaves the 12 calendar-month products into a chronological
# [canonical SPEI grid cells x 912 months] matrix.
load_spei_scale <- function(scale, pet_type, target_dir) {
  cat(sprintf("  Loading SPEI_%s-%d...\n", pet_type, scale))
  
  month_names <- c(
    "Jan", "Feb", "Mar", "Apr", "May", "Jun",
    "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"
  )
  
  expected_years <- 76L
  expected_months <- expected_years * 12L
  expected_year_vector <- 1950L:2025L
  
  # --------------------------------------------------------------------------
  # Read one calendar-month product.
  #
  # Production SPEI NetCDF files contain one spatial raster layer for each
  # year of the 1950-2025 record (76 layers). The previous version of this
  # loader used r[[1]], which silently retained only the first year and then
  # cbind()'d the 12 January-December files. That produced a cells x 12
  # matrix and caused the downstream 912-month QA failure.
  #
  # The loader now explicitly reads all annual layers and reconstructs the
  # chronological monthly matrix:
  #   Jan-1950, Feb-1950, ..., Dec-1950, Jan-1951, ..., Dec-2025.
  # --------------------------------------------------------------------------
  read_month_product <- function(month_index) {
    nc_file <- file.path(
      target_dir,
      sprintf(
        "spei_%02d_month%02d_%s.nc",
        scale,
        month_index,
        month_names[month_index]
      )
    )
    
    csv_file <- file.path(
      target_dir,
      sprintf(
        "spei_%02d_month%02d_%s.csv",
        scale,
        month_index,
        month_names[month_index]
      )
    )
    
    if (file.exists(nc_file)) {
      r <- terra::rast(nc_file)
      n_layers <- terra::nlyr(r)
      
      if (n_layers != expected_years) {
        stop(
          "Unexpected layer count in ", basename(nc_file), ": ",
          n_layers, ". Expected ", expected_years,
          " annual layers for 1950-2025. The previous loader read only ",
          "layer 1 and therefore could not reconstruct the 912-month record."
        )
      }
      
      v <- terra::values(r, mat = TRUE)
      
      if (!is.matrix(v) || nrow(v) <= 0L || ncol(v) != expected_years) {
        stop(
          "Unexpected NetCDF value dimensions in ", basename(nc_file), ": ",
          paste(dim(v), collapse = " x "),
          ". Expected cells x 76."
        )
      }
      
      if (!is.numeric(v)) {
        stop("Non-numeric SPEI values in: ", nc_file)
      }
      
      # If the NetCDF carries a usable time axis, verify that its layer order
      # corresponds to the intended 1950-2025 annual sequence.
      tt <- tryCatch(terra::time(r), error = function(e) NULL)
      if (!is.null(tt) && length(tt) == expected_years && !all(is.na(tt))) {
        tt_date <- as.Date(tt)
        tt_year <- as.integer(format(tt_date, "%Y"))
        
        if (!identical(tt_year, expected_year_vector)) {
          stop(
            "Unexpected annual layer order in ", basename(nc_file),
            ". NetCDF years are: ", paste(tt_year, collapse = ", "),
            "; expected 1950-2025."
          )
        }
      }
      
      return(v)
    }
    
    if (file.exists(csv_file)) {
      df <- read.csv(
        csv_file,
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
      
      if (ncol(df) < 3L) {
        stop(
          "SPEI CSV has fewer than 3 columns: ",
          csv_file
        )
      }
      
      # Preferred CSV structure: lon, lat, 1950, 1951, ..., 2025.
      year_cols <- names(df)[
        grepl("^(1950|195[1-9]|19[6-9][0-9]|20[0-1][0-9]|202[0-5])$", names(df))
      ]
      
      if (length(year_cols) == expected_years &&
          identical(as.integer(year_cols), expected_year_vector)) {
        v <- as.matrix(df[, year_cols, drop = FALSE])
        storage.mode(v) <- "double"
        return(v)
      }
      
      # Compatibility fallback for a simple lon/lat/value CSV. This fallback
      # can represent only one observation and is therefore not valid for the
      # 1950-2025 monthly production record.
      stop(
        "CSV fallback does not contain the expected 76 annual columns in: ",
        csv_file,
        ". Expected columns named 1950 through 2025."
      )
    }
    
    stop(
      "Missing SPEI-", scale,
      " month ", month_index,
      " file in ", target_dir,
      ". Expected either:\n  ", nc_file,
      "\nor\n  ", csv_file
    )
  }
  
  month_data <- lapply(seq_len(12L), read_month_product)
  
  n_cells_each <- vapply(
    month_data,
    nrow,
    integer(1)
  )
  
  if (length(unique(n_cells_each)) != 1L) {
    stop(
      "Inconsistent spatial cell counts among SPEI-", scale, " ", pet_type,
      " monthly files: ",
      paste(n_cells_each, collapse = ", ")
    )
  }
  
  if (any(vapply(month_data, ncol, integer(1)) != expected_years)) {
    stop(
      "Each calendar-month SPEI product must contain 76 annual observations."
    )
  }
  
  n_cells <- n_cells_each[1L]
  
  # --------------------------------------------------------------------------
  # Interleave the 12 calendar-month matrices by year.
  # --------------------------------------------------------------------------
  spei_mat <- matrix(
    NA_real_,
    nrow = n_cells,
    ncol = expected_months
  )
  
  out_col <- 1L
  
  for (year_index in seq_len(expected_years)) {
    for (month_index in seq_len(12L)) {
      spei_mat[, out_col] <- month_data[[month_index]][, year_index]
      out_col <- out_col + 1L
    }
  }
  
  if (ncol(spei_mat) != expected_months) {
    stop(
      "Failed to reconstruct SPEI-", scale, " ", pet_type,
      " monthly matrix: obtained ", ncol(spei_mat),
      " columns; expected ", expected_months, "."
    )
  }
  
  cat(sprintf(
    "    Reconstructed %d cells x %d months (1950-01 to 2025-12)\n",
    nrow(spei_mat),
    ncol(spei_mat)
  ))
  
  spei_mat
}

# Load all scales for both PET types in one pass (avoids two separate loops)
scale_data  <- lapply(setNames(decomp_scales, paste0("scale_", decomp_scales)), function(sc) {
  list(
    pm  = load_spei_scale(sc, "PM",  spei_pm_dir),
    thw = load_spei_scale(sc, "Thw", spei_thw_dir)
  )
})
spei_pm_all   <- lapply(scale_data, `[[`, "pm")
spei_thw_all  <- lapply(scale_data, `[[`, "thw")
cat("✓ SPEI results loaded from both directories\n")

# ==============================================================================
#   STEP 2: LOAD METADATA AND THE CANONICAL SPEI BASIN MASK
# ==============================================================================
cat("\n===== STEP 2: LOAD METADATA AND CANONICAL SPEI MASK =====\n")

precip <- terra::rast("monthly_data_direct/total_precipitation_monthly.nc")

dates <- as.Date(terra::time(precip))
if (is.null(dates) || length(dates) != terra::nlyr(precip) || all(is.na(dates))) {
  cat("Reconstructing dates from 1950-01...\n")
  dates <- seq(
    as.Date("1950-01-01"),
    by = "month",
    length.out = terra::nlyr(precip)
  )
}

month_nums <- as.integer(format(dates, "%m"))
year_nums  <- as.integer(format(dates, "%Y"))

if (length(dates) != 912L) {
  stop(
    "Expected 912 monthly dates (1950-01 through 2025-12), found ",
    length(dates), "."
  )
}

# The authoritative spatial geometry for the SPEI products is stored by
# 3SPEI_ERALand_v2.R in spei_analysis_grid.rds. Do NOT infer the basin from
# NA values in an arbitrary SPEI layer and do NOT use a separate canonical
# mask file that is not part of this production workflow.
spei_grid_file <- file.path(
  spei_pm_dir,
  "spei_analysis_grid.rds"
)

if (!file.exists(spei_grid_file)) {
  stop(
    "Canonical SPEI analysis grid not found: ",
    spei_grid_file,
    "\nRun 3SPEI_ERALand_v2.R first."
  )
}

spei_grid <- readRDS(spei_grid_file)

if (is.null(spei_grid$raster_template) ||
    is.null(spei_grid$basin_mask)) {
  stop(
    "Invalid spei_analysis_grid.rds: raster_template and/or basin_mask missing."
  )
}

# PackedSpatRaster objects may require unwrap(); ordinary SpatRaster objects
# can be used directly.
spei_template <- tryCatch(
  {
    if (inherits(spei_grid$raster_template, "PackedSpatRaster")) {
      terra::unwrap(spei_grid$raster_template)
    } else {
      spei_grid$raster_template
    }
  },
  error = function(e) {
    stop(
      "Could not unwrap raster_template from spei_analysis_grid.rds: ",
      e$message
    )
  }
)

if (!inherits(spei_template, "SpatRaster")) {
  stop("spei_analysis_grid.rds$raster_template is not a SpatRaster.")
}

raw_basin_mask <- spei_grid$basin_mask

basin_mask_rast <- tryCatch(
  {
    if (inherits(raw_basin_mask, "PackedSpatRaster")) {
      terra::unwrap(raw_basin_mask)
    } else if (inherits(raw_basin_mask, "SpatRaster")) {
      raw_basin_mask
    } else {
      NULL
    }
  },
  error = function(e) {
    stop(
      "Could not unwrap basin_mask from spei_analysis_grid.rds: ",
      e$message
    )
  }
)

if (!is.null(basin_mask_rast)) {
  if (!terra::compareGeom(
    basin_mask_rast,
    spei_template,
    stopOnError = FALSE,
    crs = TRUE,
    ext = TRUE,
    rowcol = TRUE,
    res = TRUE
  )) {
    stop("Canonical basin mask geometry does not match the canonical SPEI grid.")
  }
  
  basin_mask <- as.logical(
    terra::values(basin_mask_rast, mat = FALSE)
  )
} else {
  basin_mask <- as.logical(raw_basin_mask)
}

n_cells <- terra::ncell(spei_template)
n_pixels <- sum(basin_mask, na.rm = TRUE)

if (length(basin_mask) != n_cells) {
  stop(
    "Canonical basin mask length (", length(basin_mask),
    ") does not equal canonical grid cell count (", n_cells, ")."
  )
}

if (n_pixels <= 0L) {
  stop("Canonical basin mask contains zero valid basin cells.")
}

if (!is.null(spei_grid$basin_pixels) &&
    is.finite(spei_grid$basin_pixels) &&
    n_pixels != as.integer(spei_grid$basin_pixels)) {
  stop(
    "Canonical basin-pixel count mismatch: mask contains ",
    n_pixels,
    " cells but spei_analysis_grid.rds records ",
    spei_grid$basin_pixels,
    "."
  )
}

# Verify that the production SPEI raster geometry is the same canonical grid.
spei_sample <- terra::rast(
  file.path(spei_pm_dir, "spei_03_month01_Jan.nc")
)

if (!terra::compareGeom(
  spei_sample[[1]],
  spei_template,
  stopOnError = FALSE,
  crs = TRUE,
  ext = TRUE,
  rowcol = TRUE,
  res = TRUE
)) {
  stop(
    "SPEI NetCDF geometry does not match spei_analysis_grid.rds. ",
    "Do not proceed with basin averaging until the upstream SPEI grid is consistent."
  )
}

if (terra::ncell(spei_sample) != n_cells) {
  stop(
    "SPEI raster cell count (", terra::ncell(spei_sample),
    ") does not match canonical grid cell count (", n_cells, ")."
  )
}

cat(sprintf(
  "[OK] Canonical SPEI grid: %d cells; basin mask: %d cells; %d months\n",
  n_cells, n_pixels, length(dates)
))
cat(sprintf(
  "     CRS: %s\n",
  terra::crs(spei_template)
))
cat(sprintf(
  "     Resolution: %.6f x %.6f map units\n",
  terra::res(spei_template)[1],
  terra::res(spei_template)[2]
))

# EPSG:3005 is an equal-area production grid. A constant positive weight for
# every basin cell is therefore exactly equivalent to physical cell-area
# weighting; the common cell area cancels from the weighted mean.
area_weights_all <- numeric(n_cells)
area_weights_all[basin_mask] <- 1

if (length(area_weights_all) != n_cells ||
    any(!is.finite(area_weights_all)) ||
    any(area_weights_all < 0) ||
    sum(area_weights_all) <= 0) {
  stop("Failed to construct valid full-grid basin weights.")
}

cat("  Weighting: unit weights inside canonical basin; zero outside.\n")
cat("  Interpretation: exact area weighting on the equal-area SPEI grid.\n")

# ==============================================================================
#   STEP 3: COMPUTE PET-FORMULATION SENSITIVITY DIAGNOSTIC
# ==============================================================================
cat("\n===== STEP 3: COMPUTE PET-FORMULATION SENSITIVITY =====\n")

# The PM–Thornthwaite difference is retained for compatibility with 10b, but
# it is explicitly a formulation-sensitivity residual:
#
#   DeltaSPEI_PET = SPEI_PM - SPEI_Thw
#
# It is NOT an additive thermodynamic/dynamic physical decomposition.
compute_legacy_decomposition <- function(
    spei_pm,
    spei_thw,
    dates,
    month_nums,
    year_nums,
    scale,
    weights
) {
  if (!is.matrix(spei_pm) || !is.matrix(spei_thw)) {
    stop("SPEI inputs must be matrices.")
  }
  
  if (nrow(spei_pm) != nrow(spei_thw)) {
    stop(
      "PM/Thw pixel-count mismatch for SPEI-", scale, ": ",
      nrow(spei_pm), " vs ", nrow(spei_thw), "."
    )
  }
  
  if (ncol(spei_pm) != length(dates) ||
      ncol(spei_thw) != length(dates)) {
    stop(
      "Time dimension mismatch for SPEI-", scale,
      ". Matrix columns must equal supplied dates."
    )
  }
  
  if (length(weights) != nrow(spei_pm)) {
    stop(
      "Area-weight length mismatch for SPEI-", scale,
      ": ", length(weights), " weights for ", nrow(spei_pm),
      " raster cells."
    )
  }
  
  if (any(!is.finite(weights)) || any(weights < 0) ||
      sum(weights) <= 0) {
    stop("Invalid full-grid basin weights for SPEI-", scale, ".")
  }
  
  spei_pm_mean <- area_weighted_colmeans(spei_pm, weights)
  spei_thw_mean <- area_weighted_colmeans(spei_thw, weights)
  
  delta_spei_pet <- spei_pm_mean - spei_thw_mean
  
  data.frame(
    date          = dates,
    year          = year_nums,
    month         = month_nums,
    scale         = scale,
    SPEI_PM       = spei_pm_mean,
    SPEI_Thw      = spei_thw_mean,
    DeltaSPEI_PET = delta_spei_pet,
    # Backward-compatible column required by 10b. This is an alias only.
    SPEI_Dynamic  = delta_spei_pet,
    # Retained as NA compatibility fields; no physical fractions are computed.
    thm_frac      = NA_real_,
    dyn_frac      = NA_real_,
    stringsAsFactors = FALSE
  )
}

decomp_all <- dplyr::bind_rows(lapply(decomp_scales, function(sc) {
  
  pm_mat <- spei_pm_all[[paste0("scale_", sc)]]
  thw_mat <- spei_thw_all[[paste0("scale_", sc)]]
  
  if (is.null(pm_mat) || is.null(thw_mat)) {
    stop("Missing PM or Thornthwaite SPEI matrix for scale ", sc, ".")
  }
  
  if (nrow(pm_mat) != n_cells || nrow(thw_mat) != n_cells) {
    stop(
      "SPEI-", sc,
      " matrix has ", nrow(pm_mat),
      " cells but canonical grid has ", n_cells, "."
    )
  }
  
  if (ncol(pm_mat) != 912L || ncol(thw_mat) != 912L) {
    stop(
      "SPEI-", sc,
      " matrices must contain 912 monthly observations."
    )
  }
  
  result <- compute_legacy_decomposition(
    pm_mat,
    thw_mat,
    dates,
    month_nums,
    year_nums,
    sc,
    weights = area_weights_all
  )
  
  cat(sprintf("  Scale %d decomposition complete\n", sc))
  result
}))

write.csv(
  decomp_all,
  file.path(out_dir, "decomp_full_timeseries.csv"),
  row.names = FALSE
)

cat("[OK] Decomposition timeseries saved\n")

# ==============================================================================
#   STEP 4: 2022-2025 PET-FORMULATION SENSITIVITY SUMMARY
# ==============================================================================
cat("\n===== STEP 4: 2022-2025 PET-FORMULATION SENSITIVITY =====\n")

drought_decomp <- decomp_all %>%
  dplyr::filter(
    year >= 2022,
    year <= 2025,
    SPEI_PM <= THR_DROUGHT
  ) %>%
  dplyr::group_by(scale) %>%
  dplyr::summarise(
    SPEI_PM_mean          = mean(SPEI_PM, na.rm = TRUE),
    SPEI_Thw_mean         = mean(SPEI_Thw, na.rm = TRUE),
    DeltaSPEI_PET_mean    = mean(DeltaSPEI_PET, na.rm = TRUE),
    mean_abs_SPEI_PM      = mean(abs(SPEI_PM), na.rm = TRUE),
    mean_abs_SPEI_Thw     = mean(abs(SPEI_Thw), na.rm = TRUE),
    mean_abs_DeltaSPEI_PET = mean(abs(DeltaSPEI_PET), na.rm = TRUE),
    n_obs                 = dplyr::n(),
    .groups = "drop"
  )

# Legacy fraction columns are retained as NA so downstream legacy readers do
# not silently interpret them as valid physical contribution estimates.
drought_decomp$thm_frac_mean <- NA_real_
drought_decomp$dyn_frac_mean <- NA_real_
drought_decomp$thm_pct <- NA_real_
drought_decomp$dyn_pct <- NA_real_

cat("\n2022-2025 PET-formulation sensitivity summary:\n")
print(drought_decomp)

write.csv(
  drought_decomp,
  file.path(out_dir, "decomp_2022_2025_summary.csv"),
  row.names = FALSE
)

# ==============================================================================
#   STEP 5: LONG-TERM TREND IN PET-FORMULATION SENSITIVITY
# ==============================================================================
cat("\n===== STEP 5: TREND IN PET-FORMULATION SENSITIVITY =====\n")

# JJA SPEI-3 diagnostic. The response variable is the magnitude of the
# formulation-sensitivity residual, not a thermodynamic fraction.
trend_df <- decomp_all %>%
  dplyr::filter(
    scale == 3L,
    month %in% c(6L, 7L, 8L),
    SPEI_PM <= THR_DROUGHT
  ) %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(
    SPEI_PM = mean(SPEI_PM, na.rm = TRUE),
    SPEI_Thw = mean(SPEI_Thw, na.rm = TRUE),
    DeltaSPEI_PET = mean(DeltaSPEI_PET, na.rm = TRUE),
    mean_abs_DeltaSPEI_PET = mean(abs(DeltaSPEI_PET), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::filter(is.finite(mean_abs_DeltaSPEI_PET))

if (nrow(trend_df) < 10L) {
  stop("Insufficient annual JJA observations for the SPEI-3 sensitivity trend.")
}

delta_trend <- lm(mean_abs_DeltaSPEI_PET ~ year, data = trend_df)
delta_slope <- unname(coef(delta_trend)["year"])
delta_pval  <- summary(delta_trend)$coefficients["year", "Pr(>|t|)"]
delta_ci95  <- confint(delta_trend, "year", level = 0.95)

trend_post90 <- dplyr::filter(trend_df, year >= 1990L)

if (nrow(trend_post90) >= 10L) {
  delta_trend_post90 <- lm(mean_abs_DeltaSPEI_PET ~ year, data = trend_post90)
  delta_slope_post90 <- unname(coef(delta_trend_post90)["year"])
  delta_pval_post90  <- summary(delta_trend_post90)$coefficients[
    "year", "Pr(>|t|)"
  ]
} else {
  delta_trend_post90 <- NULL
  delta_slope_post90 <- NA_real_
  delta_pval_post90  <- NA_real_
}

cat(sprintf(
  "SPEI-3 JJA mean |DeltaSPEI_PET|, 1950-2025:\n  Slope = %+.6f yr^-1; p = %.4f\n",
  delta_slope,
  delta_pval
))

if (is.finite(delta_slope_post90)) {
  cat(sprintf(
    "SPEI-3 JJA mean |DeltaSPEI_PET|, 1990-2025:\n  Slope = %+.6f yr^-1; p = %.4f\n",
    delta_slope_post90,
    delta_pval_post90
  ))
}

# ==============================================================================
#   STEP 6: MANUSCRIPT / SUPPLEMENTARY FIGURES
# ==============================================================================
cat("\n===== STEP 6: GENERATE FIGURES =====\n")

theme_ms <- ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(size = 10, face = "bold", hjust = 0),
    plot.subtitle = ggplot2::element_text(
      size = 8.5, colour = "grey40", hjust = 0,
      margin = ggplot2::margin(b = 4)
    ),
    legend.position = "bottom",
    legend.key.size = ggplot2::unit(0.4, "cm"),
    legend.text = ggplot2::element_text(size = 8.5),
    legend.title = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_line(
      colour = "grey94", linewidth = 0.3
    ),
    panel.grid.minor = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(4, 8, 4, 4)
  )

COL_PM <- "#1a1a1a"
COL_THW <- "#e31a1c"
COL_DELTA <- "#1f78b4"

HIGHLIGHT_START <- as.Date("2022-01-01")
HIGHLIGHT_END <- as.Date("2025-12-31")

# ------------------------------------------------------------------------------
# Supplementary S1: SPEI-3 observed formulations and PET-sensitivity residual
# ------------------------------------------------------------------------------
df_s1 <- decomp_all %>%
  dplyr::filter(scale == 3L) %>%
  dplyr::select(date, SPEI_PM, SPEI_Thw, DeltaSPEI_PET) %>%
  tidyr::pivot_longer(
    cols = -date,
    names_to = "component",
    values_to = "SPEI"
  ) %>%
  dplyr::mutate(
    component = dplyr::recode(
      component,
      SPEI_PM = "SPEI_PM",
      SPEI_Thw = "SPEI_Thw",
      DeltaSPEI_PET = "DeltaSPEI_PET"
    )
  )

p_s1 <- ggplot2::ggplot(
  df_s1,
  ggplot2::aes(date, SPEI, colour = component)
) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey50") +
  ggplot2::geom_hline(
    yintercept = THR_EVENT,
    linetype = "dashed",
    colour = "grey40"
  ) +
  ggplot2::geom_rect(
    xmin = HIGHLIGHT_START,
    xmax = HIGHLIGHT_END,
    ymin = -Inf,
    ymax = Inf,
    inherit.aes = FALSE,
    fill = "grey85",
    alpha = 0.45
  ) +
  ggplot2::geom_line(linewidth = 0.6, alpha = 0.85) +
  ggplot2::scale_colour_manual(
    values = c(
      SPEI_PM = COL_PM,
      SPEI_Thw = COL_THW,
      DeltaSPEI_PET = COL_DELTA
    )
  ) +
  ggplot2::labs(
    title = "Nechako Basin: SPEI-3 formulation sensitivity",
    subtitle = "DeltaSPEI_PET = SPEI_PM - SPEI_Thw; grey band = 2022-2025",
    x = NULL,
    y = "Index value",
    colour = NULL
  ) +
  theme_ms

ggplot2::ggsave(
  file.path(out_dir, "fig_decomp_spei3_timeseries.pdf"),
  p_s1,
  width = 10,
  height = 6,
  units = "in",
  device = "pdf"
)

# ------------------------------------------------------------------------------
# Manuscript Fig. 5a: SPEI-1 and SPEI-3 time series
# ------------------------------------------------------------------------------
ms_scales_ts <- MANUSCRIPT_SCALES

panels_ts <- lapply(seq_along(ms_scales_ts), function(i) {
  
  sc <- ms_scales_ts[i]
  
  df_sc <- decomp_all %>%
    dplyr::filter(
      scale == sc,
      date >= as.Date("2020-01-01"),
      date <= as.Date("2025-12-31")
    ) %>%
    dplyr::select(date, SPEI_PM, SPEI_Thw, DeltaSPEI_PET) %>%
    tidyr::pivot_longer(
      cols = -date,
      names_to = "component",
      values_to = "value"
    )
  
  ggplot2::ggplot(
    df_sc,
    ggplot2::aes(
      x = date,
      y = value,
      colour = component,
      linetype = component
    )
  ) +
    ggplot2::annotate(
      "rect",
      xmin = HIGHLIGHT_START,
      xmax = HIGHLIGHT_END,
      ymin = -Inf,
      ymax = Inf,
      fill = "grey85",
      alpha = 0.45
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      colour = "grey50",
      linewidth = 0.35
    ) +
    ggplot2::geom_hline(
      yintercept = THR_EVENT,
      colour = "grey35",
      linewidth = 0.5,
      linetype = "dashed"
    ) +
    ggplot2::geom_line(linewidth = 0.65, alpha = 0.9) +
    ggplot2::scale_colour_manual(
      values = c(
        SPEI_PM = COL_PM,
        SPEI_Thw = COL_THW,
        DeltaSPEI_PET = COL_DELTA
      ),
      labels = c(
        SPEI_PM = "SPEI_PM",
        SPEI_Thw = "SPEI_Thw",
        DeltaSPEI_PET = "DeltaSPEI_PET"
      )
    ) +
    ggplot2::scale_linetype_manual(
      values = c(
        SPEI_PM = "solid",
        SPEI_Thw = "solid",
        DeltaSPEI_PET = "dashed"
      ),
      labels = c(
        SPEI_PM = "SPEI_PM",
        SPEI_Thw = "SPEI_Thw",
        DeltaSPEI_PET = "DeltaSPEI_PET"
      )
    ) +
    ggplot2::scale_x_date(
      date_breaks = "1 year",
      date_labels = "%Y",
      expand = ggplot2::expansion(add = c(15, 15))
    ) +
    ggplot2::labs(
      title = sprintf("(%s) SPEI-%d", letters[i], sc),
      subtitle = "DeltaSPEI_PET is a PET-formulation sensitivity residual, not a physical component.",
      x = if (i == length(ms_scales_ts)) "Year" else NULL,
      y = "Index value",
      colour = NULL,
      linetype = NULL
    ) +
    theme_ms +
    ggplot2::theme(
      axis.text.x = if (i == length(ms_scales_ts)) {
        ggplot2::element_text(size = 7.5)
      } else {
        ggplot2::element_blank()
      },
      axis.ticks.x = if (i == length(ms_scales_ts)) {
        ggplot2::element_line()
      } else {
        ggplot2::element_blank()
      }
    )
})

fig5a <- patchwork::wrap_plots(panels_ts, ncol = 1) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(
    title = "Observed SPEI formulations and PET-formulation sensitivity",
    subtitle = paste0(
      "Nechako River Basin, 2020-2025. ",
      "Black = SPEI_PM; red = SPEI_Thw; blue dashed = DeltaSPEI_PET."
    ),
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = 10, face = "bold", hjust = 0.5
      ),
      plot.subtitle = ggplot2::element_text(
        size = 7.5, colour = "grey35", hjust = 0
      ),
      legend.position = "bottom"
    )
  )

ggplot2::ggsave(
  file.path(out_dir, "Fig5a_decomp_timeseries_2020_2025.pdf"),
  fig5a,
  width = 7,
  height = 7.5,
  units = "in",
  device = "pdf"
)

ggplot2::ggsave(
  file.path(out_dir, "Fig5a_decomp_timeseries_2020_2025.png"),
  fig5a,
  width = 7,
  height = 7.5,
  units = "in",
  dpi = 300,
  device = "png"
)

# ------------------------------------------------------------------------------
# Manuscript Fig. 5b: non-additive magnitude comparison
# ------------------------------------------------------------------------------
decomp_ms_bar <- drought_decomp %>%
  dplyr::filter(scale %in% MANUSCRIPT_SCALES) %>%
  dplyr::select(
    scale,
    mean_abs_SPEI_Thw,
    mean_abs_DeltaSPEI_PET
  ) %>%
  tidyr::pivot_longer(
    cols = c(mean_abs_SPEI_Thw, mean_abs_DeltaSPEI_PET),
    names_to = "metric",
    values_to = "magnitude"
  ) %>%
  dplyr::mutate(
    component = dplyr::recode(
      metric,
      mean_abs_SPEI_Thw = "mean |SPEI_Thw|",
      mean_abs_DeltaSPEI_PET = "mean |DeltaSPEI_PET|"
    ),
    scale_lab = factor(
      paste0("SPEI-", scale),
      levels = paste0("SPEI-", MANUSCRIPT_SCALES)
    )
  )

fig5b <- ggplot2::ggplot(
  decomp_ms_bar,
  ggplot2::aes(
    x = scale_lab,
    y = magnitude,
    fill = component
  )
) +
  ggplot2::geom_col(
    position = ggplot2::position_dodge(width = 0.7),
    width = 0.6
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "mean |SPEI_Thw|" = COL_THW,
      "mean |DeltaSPEI_PET|" = COL_DELTA
    )
  ) +
  ggplot2::labs(
    title = "Magnitude of PET-formulation sensitivity during 2022-2025 drought",
    subtitle = paste0(
      "Means are calculated only for months with SPEI_PM <= ",
      THR_DROUGHT,
      ". Values are not additive contributions."
    ),
    x = "Accumulation scale",
    y = "Mean absolute index magnitude",
    fill = NULL
  ) +
  theme_ms

ggplot2::ggsave(
  file.path(out_dir, "Fig5b_scale_comparison.pdf"),
  fig5b,
  width = 5.5,
  height = 4.8,
  units = "in",
  device = "pdf"
)

ggplot2::ggsave(
  file.path(out_dir, "Fig5b_scale_comparison.png"),
  fig5b,
  width = 5.5,
  height = 4.8,
  units = "in",
  dpi = 300,
  device = "png"
)

# ------------------------------------------------------------------------------
# Manuscript Fig. 6: JJA trend in mean |DeltaSPEI_PET|
# ------------------------------------------------------------------------------
pred_full <- predict(
  delta_trend,
  newdata = data.frame(year = trend_df$year),
  interval = "confidence",
  level = 0.95
)

trend_df$fit_full <- pred_full[, "fit"]
trend_df$lwr95 <- pred_full[, "lwr"]
trend_df$upr95 <- pred_full[, "upr"]

fig6 <- ggplot2::ggplot() +
  ggplot2::annotate(
    "rect",
    xmin = 2022,
    xmax = 2026,
    ymin = -Inf,
    ymax = Inf,
    fill = "grey85",
    alpha = 0.5
  ) +
  ggplot2::geom_ribbon(
    data = trend_df,
    ggplot2::aes(
      x = year,
      ymin = lwr95,
      ymax = upr95
    ),
    fill = COL_DELTA,
    alpha = 0.15
  ) +
  ggplot2::geom_line(
    data = trend_df,
    ggplot2::aes(x = year, y = fit_full),
    colour = COL_DELTA,
    linewidth = 0.9
  ) +
  ggplot2::geom_point(
    data = trend_df,
    ggplot2::aes(x = year, y = mean_abs_DeltaSPEI_PET),
    size = 2,
    colour = COL_DELTA
  ) +
  ggplot2::geom_vline(
    xintercept = 2022,
    linetype = "dashed",
    colour = "grey50"
  ) +
  ggplot2::scale_x_continuous(
    breaks = seq(1950, 2025, by = 10),
    expand = ggplot2::expansion(add = c(1, 2))
  ) +
  ggplot2::labs(
    title = "Temporal trend in PET-formulation sensitivity of JJA SPEI-3",
    subtitle = paste0(
      "Annual JJA mean |DeltaSPEI_PET|, 1950-2025. ",
      "OLS slope = ", sprintf("%+.5f yr^-1", delta_slope),
      "; p = ", format.pval(delta_pval, digits = 3, eps = 0.001),
      "."
    ),
    x = "Year",
    y = "Mean |DeltaSPEI_PET|"
  ) +
  theme_ms

ggplot2::ggsave(
  file.path(out_dir, "Fig6_PET_formulation_sensitivity_trend.pdf"),
  fig6,
  width = 7.5,
  height = 5,
  units = "in",
  device = "pdf"
)

ggplot2::ggsave(
  file.path(out_dir, "Fig6_PET_formulation_sensitivity_trend.png"),
  fig6,
  width = 7.5,
  height = 5,
  units = "in",
  dpi = 300,
  device = "png"
)

# ==============================================================================
#   STEP 7: SEASONAL 2022-2025 PET-FORMULATION SENSITIVITY
# ==============================================================================
# IMPORTANT: the requested seasonal drought filtering is intentionally retained:
# only observations with SPEI_PM <= THR_DROUGHT are included.
cat("\n===== STEP 7: SEASONAL SENSITIVITY 2022-2025 =====\n")

MONTH_ABBR <- c(
  "Jan", "Feb", "Mar", "Apr", "May", "Jun",
  "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"
)

SEAS_SCALES <- MANUSCRIPT_SCALES

seas_stats <- decomp_all %>%
  dplyr::filter(
    year >= 2022,
    year <= 2025,
    scale %in% SEAS_SCALES,
    SPEI_PM <= THR_DROUGHT
  ) %>%
  dplyr::mutate(
    month_abb = factor(
      MONTH_ABBR[month],
      levels = MONTH_ABBR
    )
  ) %>%
  dplyr::group_by(scale, month, month_abb) %>%
  dplyr::summarise(
    pm_mean = mean(SPEI_PM, na.rm = TRUE),
    thw_mean = mean(SPEI_Thw, na.rm = TRUE),
    delta_mean = mean(DeltaSPEI_PET, na.rm = TRUE),
    
    pm_sd = sd(SPEI_PM, na.rm = TRUE),
    thw_sd = sd(SPEI_Thw, na.rm = TRUE),
    delta_sd = sd(DeltaSPEI_PET, na.rm = TRUE),
    
    pm_abs = mean(abs(SPEI_PM), na.rm = TRUE),
    thw_abs = mean(abs(SPEI_Thw), na.rm = TRUE),
    delta_abs = mean(abs(DeltaSPEI_PET), na.rm = TRUE),
    
    n_obs = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    pm_lo = pm_mean - pm_sd,
    pm_hi = pm_mean + pm_sd,
    thw_lo = thw_mean - thw_sd,
    thw_hi = thw_mean + thw_sd,
    delta_lo = delta_mean - delta_sd,
    delta_hi = delta_mean + delta_sd
  )

cat(sprintf(
  "  Seasonal statistics: %d rows across SPEI-%s.\n",
  nrow(seas_stats),
  paste(SEAS_SCALES, collapse = ", ")
))

theme_seas <- ggplot2::theme_classic(base_size = 9) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      size = 9, face = "bold", hjust = 0
    ),
    plot.subtitle = ggplot2::element_text(
      size = 7.5, colour = "grey35", hjust = 0,
      margin = ggplot2::margin(b = 3)
    ),
    axis.title.x = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(size = 7.5),
    axis.text.y = ggplot2::element_text(size = 7.5),
    axis.title.y = ggplot2::element_text(size = 8),
    legend.position = "bottom",
    legend.key.size = ggplot2::unit(0.38, "cm"),
    legend.text = ggplot2::element_text(size = 7.5),
    legend.title = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_line(
      colour = "grey92", linewidth = 0.3
    ),
    panel.grid.major.x = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(3, 6, 2, 4)
  )

season_panels <- lapply(seq_along(SEAS_SCALES), function(i) {
  
  sc <- SEAS_SCALES[i]
  
  d <- seas_stats %>%
    dplyr::filter(scale == sc)
  
  ggplot2::ggplot(
    d,
    ggplot2::aes(x = month_abb)
  ) +
    ggplot2::geom_hline(
      yintercept = 0,
      colour = "grey65",
      linewidth = 0.3
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = pm_lo, ymax = pm_hi, group = 1),
      fill = COL_PM,
      alpha = 0.08
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = thw_lo, ymax = thw_hi, group = 1),
      fill = COL_THW,
      alpha = 0.10
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = delta_lo, ymax = delta_hi, group = 1),
      fill = COL_DELTA,
      alpha = 0.10
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = pm_mean, colour = "SPEI_PM", group = 1),
      linewidth = 0.65
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = thw_mean, colour = "SPEI_Thw", group = 1),
      linewidth = 0.65
    ) +
    ggplot2::geom_line(
      ggplot2::aes(
        y = delta_mean,
        colour = "DeltaSPEI_PET",
        group = 1
      ),
      linewidth = 0.65,
      linetype = "dashed"
    ) +
    ggplot2::scale_colour_manual(
      values = c(
        SPEI_PM = COL_PM,
        SPEI_Thw = COL_THW,
        DeltaSPEI_PET = COL_DELTA
      )
    ) +
    ggplot2::labs(
      title = sprintf("SPEI-%d", sc),
      subtitle = "Means ± 1 SD; drought gate SPEI_PM <= -0.5 retained.",
      y = "Index value",
      colour = NULL
    ) +
    theme_seas
  
})

fig5c <- patchwork::wrap_plots(season_panels, ncol = 1) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(
    title = "Seasonal PET-formulation sensitivity, 2022-2025",
    subtitle = paste0(
      "DeltaSPEI_PET = SPEI_PM - SPEI_Thw. ",
      "This residual is not interpreted as a physical thermodynamic/dynamic decomposition."
    )
  )

ggplot2::ggsave(
  file.path(out_dir, "Fig5c_seasonal_sensitivity_2022_2025.pdf"),
  fig5c,
  width = 7,
  height = 7.5,
  units = "in",
  device = "pdf"
)

ggplot2::ggsave(
  file.path(out_dir, "Fig5c_seasonal_sensitivity_2022_2025.png"),
  fig5c,
  width = 7,
  height = 7.5,
  units = "in",
  dpi = 300,
  device = "png"
)

# Seasonal magnitude comparison; retained seasonal drought filter.
season_bar <- seas_stats %>%
  dplyr::select(
    scale,
    month_abb,
    thw_abs,
    delta_abs
  ) %>%
  tidyr::pivot_longer(
    cols = c(thw_abs, delta_abs),
    names_to = "metric",
    values_to = "magnitude"
  ) %>%
  dplyr::mutate(
    component = dplyr::recode(
      metric,
      thw_abs = "mean |SPEI_Thw|",
      delta_abs = "mean |DeltaSPEI_PET|"
    )
  )

fig5d <- ggplot2::ggplot(
  season_bar,
  ggplot2::aes(
    x = month_abb,
    y = magnitude,
    fill = component
  )
) +
  ggplot2::geom_col(
    position = ggplot2::position_dodge(width = 0.75),
    width = 0.65
  ) +
  ggplot2::facet_wrap(
    ~ scale,
    ncol = 1,
    labeller = ggplot2::labeller(
      scale = function(x) paste0("SPEI-", x)
    )
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "mean |SPEI_Thw|" = COL_THW,
      "mean |DeltaSPEI_PET|" = COL_DELTA
    )
  ) +
  ggplot2::labs(
    title = "Seasonal magnitude of PET-formulation sensitivity",
    subtitle = "2022-2025 drought-gated observations; magnitudes are non-additive.",
    x = "Calendar month",
    y = "Mean absolute index magnitude",
    fill = NULL
  ) +
  theme_seas

ggplot2::ggsave(
  file.path(out_dir, "Fig5d_seasonal_sensitivity_2022_2025.pdf"),
  fig5d,
  width = 7,
  height = 8,
  units = "in",
  device = "pdf"
)

ggplot2::ggsave(
  file.path(out_dir, "Fig5d_seasonal_sensitivity_2022_2025.png"),
  fig5d,
  width = 7,
  height = 8,
  units = "in",
  dpi = 300,
  device = "png"
)

# ==============================================================================
#   FINAL QA / OUTPUT SUMMARY
# ==============================================================================
cat("\n===== FINAL QA =====\n")

stopifnot(
  nrow(decomp_all) == length(decomp_scales) * length(dates),
  all(decomp_all$scale %in% decomp_scales),
  all(is.finite(area_weights_all)),
  all(area_weights_all >= 0),
  sum(area_weights_all) > 0
)

cat(sprintf(
  "[OK] %d scales x %d months = %d decomposition records.\n",
  length(decomp_scales),
  length(dates),
  nrow(decomp_all)
))

cat(sprintf(
  "[OK] Canonical grid: %d cells; basin cells: %d.\n",
  n_cells,
  n_pixels
))

cat("[OK] No unweighted colMeans() fallback is permitted.\n")
cat("[OK] No SPEI-2 is requested; manuscript scales are SPEI-1 and SPEI-3.\n")
cat("[OK] DeltaSPEI_PET is treated as a PET-formulation sensitivity residual.\n")
cat("[OK] Seasonal drought filtering was intentionally retained.\n")
cat(sprintf(
  "[OK] Output directory: %s\n",
  normalizePath(out_dir, winslash = "/", mustWork = FALSE)
))

cat("\n10a PET-formulation sensitivity diagnostic complete.\n")
