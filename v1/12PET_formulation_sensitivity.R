################################################################################
# 12PET_source_sensitivity_v2.R
# PET-FORMULATION SENSITIVITY — CURRENT MANUSCRIPT RESULTS (Methods/Results 3e/4e)
#
# Purpose
#   Compare independently standardized SPEI calculated with FAO-56
#   Penman-Monteith PET (SPEI_PM) and ERA5-Land native PEV PET (SPEI_ERA5).
#
# IMPORTANT INTERPRETATION
#   delta_spei_pet = SPEI_PM - SPEI_ERA5 is a PET-source sensitivity residual.
#   It is NOT an additive decomposition of physical drought drivers.
#
# Primary metrics:
#   1. Pearson r
#   2. RMSE
#   3. drought-state agreement (%; threshold = DROUGHT_ONSET)
#   4. mean absolute delta SPEI
#
# Metrics are reported for:
#   full record; 2022-2025; DJF/MAM/JJA/SON; warm season (MJJAS);
#   and drought-active months.
#
# Inputs
#   Authoritative basin-average CSVs loaded by load_basin_avg_csv():
#     spei_results_seasonal/spei_XX_basin_averaged_by_month.csv
#     spei_results_seasonal_era5land/spei_era5_XX_basin_averaged_by_month.csv
#
# Outputs -> pet_sensitivity/
#   PET_sensitivity_monthly.csv
#   PET_agreement_by_scale_period.csv
#   PET_divergence_extremes.csv
#   PET_agreement_bootstrap_CI.csv
#   PET_divergence_by_drought_intensity.csv
#   PET_divergence_drought_intensity_association.csv
#   Fig_PET_source_sensitivity_timeseries.pdf/.png
#   Fig_PET_source_agreement_metrics.pdf/.png
################################################################################

setwd(Sys.getenv("NECHAKO_WD", "D:/Nechako_Drought/Nechako/"))
source("DROUGHT_ANALYSIS_utils_v2.R")
utils_load_packages(c("dplyr", "tidyr", "ggplot2", "data.table", "lubridate",
                      "patchwork", "scales"))

if (!dir.exists(WD_PATH)) stop("Working directory not found: ", WD_PATH)
setwd(WD_PATH)

OUT_DIR <- file.path(WD_PATH, "pet_sensitivity")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

ANALYSIS_START <- as.Date("1950-01-01")
ANALYSIS_END   <- as.Date("2025-12-01")
RECENT_START   <- as.Date("2022-01-01")
RECENT_END     <- as.Date("2025-12-01")
DROUGHT_THR    <- DROUGHT_ONSET
SCALES         <- SPEI_SCALES
BOOT_R         <- 1000L
BOOT_BLOCK     <- 12L
BOOT_SEED      <- 12013L

season_of <- function(m) {
  dplyr::case_when(
    m %in% c(12L, 1L, 2L) ~ "DJF",
    m %in% 3:5             ~ "MAM",
    m %in% 6:8             ~ "JJA",
    m %in% 9:11            ~ "SON",
    TRUE                   ~ NA_character_
  )
}

load_pair <- function(scale) {
  pm  <- load_basin_avg_csv("spei", scale)
  era5 <- load_basin_avg_csv("spei_era5", scale)
  if (is.null(pm) || is.null(era5))
    stop("Missing basin-average SPEI input for scale ", scale,
         ". Run 3SPEI_ERALand_v2.R and 6trend_test_ALL_v2.R first.")
  
  validate_monthly_record(pm$date,  label = paste0("SPEI_PM-", scale))
  validate_monthly_record(era5$date, label = paste0("SPEI_ERA5-", scale))
  
  z <- merge(pm, era5, by = "date", all = TRUE, suffixes = c("_pm", "_era5"))
  z <- z[z$date >= ANALYSIS_START & z$date <= ANALYSIS_END, , drop = FALSE]
  validate_monthly_record(z$date, ANALYSIS_START, ANALYSIS_END,
                          label = paste0("paired SPEI-", scale),
                          require_complete = TRUE)
  
  z$scale <- as.integer(scale)
  z$year  <- as.integer(format(z$date, "%Y"))
  z$month <- as.integer(format(z$date, "%m"))
  z$season <- season_of(z$month)
  # Year/season-year cluster used by the dependence-aware bootstrap.
  # DJF months are deliberately kept within the same calendar-year cluster
  # (December, January, February) so every complete analysis year contributes
  # one seasonal cluster without creating artificial edge years.
  z$bootstrap_year <- z$year
  z$SPEI_PM  <- as.numeric(z$value_pm)
  z$SPEI_ERA5 <- as.numeric(z$value_era5)
  z$delta_spei_pet <- z$SPEI_PM - z$SPEI_ERA5
  z$abs_delta_spei_pet <- abs(z$delta_spei_pet)
  z$drought_pm  <- is.finite(z$SPEI_PM)  & z$SPEI_PM  < DROUGHT_THR
  z$drought_era5 <- is.finite(z$SPEI_ERA5) & z$SPEI_ERA5 < DROUGHT_THR
  z$drought_active <- z$drought_pm | z$drought_era5
  z
}

all_ts <- dplyr::bind_rows(lapply(SCALES, load_pair))
validate_scales(unique(all_ts$scale), expected = SCALES, allow_extra = FALSE,
                label = "PET-sensitivity scales")

write.csv(all_ts, file.path(OUT_DIR, "PET_sensitivity_monthly.csv"),
          row.names = FALSE)

metric_row <- function(df, scale, period) {
  ok <- is.finite(df$SPEI_PM) & is.finite(df$SPEI_ERA5)
  x <- df$SPEI_PM[ok]
  y <- df$SPEI_ERA5[ok]
  n <- length(x)
  
  data.frame(
    scale = as.integer(scale),
    period = period,
    n_months = n,
    pearson_r = if (n >= 3L && stats::sd(x) > 0 && stats::sd(y) > 0)
      stats::cor(x, y, method = "pearson") else NA_real_,
    rmse = if (n) sqrt(mean((x - y)^2)) else NA_real_,
    drought_state_agreement_pct = if (n)
      100 * mean((x < DROUGHT_THR) == (y < DROUGHT_THR)) else NA_real_,
    mean_abs_delta_spei = if (n) mean(abs(x - y)) else NA_real_,
    mean_delta_spei = if (n) mean(x - y) else NA_real_,
    stringsAsFactors = FALSE
  )
}

period_subsets <- function(df) {
  list(
    "Full record" = df,
    "2022-2025" = df[df$date >= RECENT_START & df$date <= RECENT_END, ],
    "DJF" = df[df$season == "DJF", ],
    "MAM" = df[df$season == "MAM", ],
    "JJA" = df[df$season == "JJA", ],
    "SON" = df[df$season == "SON", ],
    "Warm season (MJJAS)" = df[df$month %in% 5:9, ],
    "Drought-active" = df[df$drought_active %in% TRUE, ],
    "2022-2025 drought-active" =
      df[df$date >= RECENT_START & df$date <= RECENT_END &
           df$drought_active %in% TRUE, ]
  )
}

metrics <- dplyr::bind_rows(lapply(split(all_ts, all_ts$scale), function(df) {
  ss <- period_subsets(df)
  dplyr::bind_rows(Map(function(x, nm) metric_row(x, unique(df$scale), nm),
                       ss, names(ss)))
}))

# ------------------------------------------------------------------------------
# DEPENDENCE-AWARE BOOTSTRAP CONFIDENCE INTERVALS FOR AGREEMENT METRICS
#
# Bootstrap design:
#   * Full record and contiguous 2022-2025 record:
#       circular monthly block bootstrap with 12-month blocks.
#   * Seasonal, warm-season, and drought-active subsets:
#       year-cluster bootstrap. Whole calendar years are resampled together,
#       avoiding the invalid assumption that filtered seasonal observations
#       form a contiguous monthly time series.
#
# PERFORMANCE / PROGRESS:
#   The bootstrap operates on numeric vectors and integer row indices only.
#   No dplyr::bind_rows() or data-frame construction occurs inside the
#   replicate loop.
#
#   Console progress is reported:
#       - once when each scale/period begins;
#       - every BOOT_PROGRESS_STEP replicates;
#       - at completion of every scale/period;
#       - cumulative elapsed time is reported for the complete bootstrap.
#
# IMPORTANT:
#   For seasonal/irregular subsets with fewer than 8 independent calendar
#   years, a dependence-aware CI is not reported.
# ------------------------------------------------------------------------------

set.seed(BOOT_SEED)

BOOT_PROGRESS_STEP <- 100L

bootstrap_time_start <- Sys.time()

bootstrap_elapsed <- function(t0) {
  as.numeric(Sys.time() - t0, units = "secs")
}

format_elapsed <- function(seconds) {
  seconds <- as.numeric(seconds)
  
  if (!is.finite(seconds)) {
    return("NA")
  }
  
  if (seconds < 60) {
    sprintf("%.1f s", seconds)
  } else {
    sprintf("%dm %.1fs", floor(seconds / 60), seconds %% 60)
  }
}

bootstrap_progress <- function(
    scale,
    period,
    b,
    B,
    t0) {
  
  pct <- 100 * b / B
  
  cat(
    sprintf(
      "[%s] Bootstrap scale %d | %-28s | %4d/%d (%5.1f%%) | elapsed %s\n",
      format(Sys.time(), "%H:%M:%S"),
      as.integer(scale),
      period,
      as.integer(b),
      as.integer(B),
      pct,
      format_elapsed(bootstrap_elapsed(t0))
    )
  )
  
  flush.console()
}

circular_block_indices <- function(
    n,
    block_length = BOOT_BLOCK) {
  
  n <- as.integer(n)
  block_length <- as.integer(block_length)
  
  if (n <= 1L) {
    return(seq_len(n))
  }
  
  if (block_length < 1L) {
    stop("block_length must be >= 1.")
  }
  
  block_length <- min(block_length, n)
  
  starts <- sample.int(
    n,
    ceiling(n / block_length),
    replace = TRUE
  )
  
  idx <- unlist(
    lapply(
      starts,
      function(st) {
        ((st - 1L + seq_len(block_length) - 1L) %% n) + 1L
      }
    ),
    use.names = FALSE
  )
  
  idx[seq_len(n)]
}

bootstrap_metric_ci <- function(
    df,
    scale,
    period,
    B = BOOT_R,
    block_length = BOOT_BLOCK) {
  
  ok <- is.finite(df$SPEI_PM) &
    is.finite(df$SPEI_ERA5)
  
  d <- df[ok, , drop = FALSE]
  
  x <- as.numeric(d$SPEI_PM)
  y <- as.numeric(d$SPEI_ERA5)
  
  n <- length(x)
  
  empty_result <- data.frame(
    scale = as.integer(scale),
    period = period,
    n_months = n,
    bootstrap_reps = 0L,
    bootstrap_method = NA_character_,
    block_length = NA_integer_,
    pearson_r_ci_low = NA_real_,
    pearson_r_ci_high = NA_real_,
    rmse_ci_low = NA_real_,
    rmse_ci_high = NA_real_,
    drought_state_agreement_pct_ci_low = NA_real_,
    drought_state_agreement_pct_ci_high = NA_real_,
    mean_abs_delta_spei_ci_low = NA_real_,
    mean_abs_delta_spei_ci_high = NA_real_,
    stringsAsFactors = FALSE
  )
  
  if (n < 12L) {
    empty_result$bootstrap_method <- "insufficient_months"
    return(empty_result)
  }
  
  n_clusters <- if ("bootstrap_year" %in% names(d)) {
    length(
      unique(
        d$bootstrap_year[
          is.finite(d$bootstrap_year)
        ]
      )
    )
  } else {
    NA_integer_
  }
  
  monthly_block_period <- period %in% c(
    "Full record",
    "2022-2025"
  )
  
  if (
    !monthly_block_period &&
    (
      !is.finite(n_clusters) ||
      n_clusters < 8L
    )
  ) {
    empty_result$bootstrap_method <-
      "year_cluster_insufficient_years"
    return(empty_result)
  }
  
  if (monthly_block_period) {
    
    method <- "circular_monthly_block"
    effective_block <- min(
      as.integer(block_length),
      n
    )
    
    year_values <- NULL
    year_groups <- NULL
    
  } else {
    
    method <- "year_cluster"
    effective_block <- NA_integer_
    
    year_values <- sort(
      unique(
        as.integer(d$bootstrap_year)
      )
    )
    
    year_groups <- lapply(
      year_values,
      function(y) which(d$bootstrap_year == y)
    )
    
    names(year_groups) <- as.character(year_values)
  }
  
  vals <- matrix(
    NA_real_,
    nrow = B,
    ncol = 4L,
    dimnames = list(
      NULL,
      c(
        "pearson_r",
        "rmse",
        "drought_state_agreement_pct",
        "mean_abs_delta_spei"
      )
    )
  )
  
  t0 <- Sys.time()
  
  cat(
    sprintf(
      "[%s] START bootstrap | scale %d | %-28s | n=%d | B=%d | method=%s\n",
      format(Sys.time(), "%H:%M:%S"),
      as.integer(scale),
      period,
      n,
      as.integer(B),
      method
    )
  )
  flush.console()
  
  for (b in seq_len(B)) {
    
    if (monthly_block_period) {
      
      ii <- circular_block_indices(
        n,
        block_length = block_length
      )
      
    } else {
      
      sampled_years <- sample(
        year_values,
        size = length(year_values),
        replace = TRUE
      )
      
      ii <- unlist(
        year_groups[
          as.character(sampled_years)
        ],
        use.names = FALSE
      )
    }
    
    xb <- x[ii]
    yb <- y[ii]
    
    sx <- stats::sd(xb)
    sy <- stats::sd(yb)
    
    vals[b, 1L] <- if (
      length(xb) >= 3L &&
      is.finite(sx) &&
      is.finite(sy) &&
      sx > 0 &&
      sy > 0
    ) {
      stats::cor(xb, yb, method = "pearson")
    } else {
      NA_real_
    }
    
    vals[b, 2L] <- sqrt(
      mean((xb - yb)^2)
    )
    
    vals[b, 3L] <- 100 * mean(
      (xb < DROUGHT_THR) ==
        (yb < DROUGHT_THR)
    )
    
    vals[b, 4L] <- mean(
      abs(xb - yb)
    )
    
    if (
      b == 1L ||
      b == B ||
      b %% BOOT_PROGRESS_STEP == 0L
    ) {
      bootstrap_progress(
        scale = scale,
        period = period,
        b = b,
        B = B,
        t0 = t0
      )
    }
  }
  
  ci <- function(z) {
    
    z <- z[is.finite(z)]
    
    if (!length(z)) {
      return(c(NA_real_, NA_real_))
    }
    
    as.numeric(
      stats::quantile(
        z,
        probs = c(0.025, 0.975),
        names = FALSE,
        type = 7
      )
    )
  }
  
  c1 <- ci(vals[, 1L])
  c2 <- ci(vals[, 2L])
  c3 <- ci(vals[, 3L])
  c4 <- ci(vals[, 4L])
  
  elapsed <- bootstrap_elapsed(t0)
  
  cat(
    sprintf(
      "[%s] DONE bootstrap | scale %d | %-28s | %d reps | %s\n",
      format(Sys.time(), "%H:%M:%S"),
      as.integer(scale),
      period,
      as.integer(B),
      format_elapsed(elapsed)
    )
  )
  flush.console()
  
  data.frame(
    scale = as.integer(scale),
    period = period,
    n_months = n,
    bootstrap_reps = B,
    bootstrap_method = method,
    block_length = effective_block,
    pearson_r_ci_low = c1[1L],
    pearson_r_ci_high = c1[2L],
    rmse_ci_low = c2[1L],
    rmse_ci_high = c2[2L],
    drought_state_agreement_pct_ci_low = c3[1L],
    drought_state_agreement_pct_ci_high = c3[2L],
    mean_abs_delta_spei_ci_low = c4[1L],
    mean_abs_delta_spei_ci_high = c4[2L],
    stringsAsFactors = FALSE
  )
}

cat(
  sprintf(
    "\n[%s] ============================================================\n",
    format(Sys.time(), "%H:%M:%S")
  )
)
cat(
  sprintf(
    "[%s] STARTING DEPENDENCE-AWARE BOOTSTRAP | scales=%d | periods/scale=%d | B=%d\n",
    format(Sys.time(), "%H:%M:%S"),
    length(SCALES),
    9L,
    BOOT_R
  )
)
cat(
  sprintf(
    "[%s] Progress interval: every %d replicates\n",
    format(Sys.time(), "%H:%M:%S"),
    BOOT_PROGRESS_STEP
  )
)
cat(
  sprintf(
    "[%s] ============================================================\n\n",
    format(Sys.time(), "%H:%M:%S")
  )
)
flush.console()

bootstrap_ci <- dplyr::bind_rows(
  lapply(
    split(all_ts, all_ts$scale),
    function(df) {
      
      scale_id <- unique(df$scale)
      ss <- period_subsets(df)
      
      dplyr::bind_rows(
        lapply(
          names(ss),
          function(nm) {
            bootstrap_metric_ci(
              df = ss[[nm]],
              scale = scale_id,
              period = nm
            )
          }
        )
      )
    }
  )
)

bootstrap_total_elapsed <- bootstrap_elapsed(
  bootstrap_time_start
)

cat(
  sprintf(
    "\n[%s] ============================================================\n",
    format(Sys.time(), "%H:%M:%S")
  )
)
cat(
  sprintf(
    "[%s] BOOTSTRAP COMPLETE | total elapsed: %s\n",
    format(Sys.time(), "%H:%M:%S"),
    format_elapsed(bootstrap_total_elapsed)
  )
)
cat(
  sprintf(
    "[%s] ============================================================\n\n",
    format(Sys.time(), "%H:%M:%S")
  )
)
flush.console()

metrics <- dplyr::left_join(
  metrics,
  bootstrap_ci,
  by = c("scale", "period", "n_months")
)

write.csv(
  metrics,
  file.path(OUT_DIR, "PET_agreement_by_scale_period.csv"),
  row.names = FALSE
)

write.csv(
  bootstrap_ci,
  file.path(OUT_DIR, "PET_agreement_bootstrap_CI.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# PET-FORMULATION DIVERGENCE AS A FUNCTION OF DROUGHT INTENSITY
#
# This addresses whether |Delta SPEI_PET| becomes systematically larger as the
# diagnosed drought becomes more severe. Severity is represented conservatively
# by the mean of the two independently standardized SPEI formulations; the
# resulting statistic remains a PET-source sensitivity diagnostic, not a causal
# drought decomposition.
# ------------------------------------------------------------------------------
intensity_tbl <- all_ts |>
  dplyr::mutate(
    mean_spei = (SPEI_PM + SPEI_ERA5) / 2,
    drought_intensity = pmax(0, -mean_spei),
    drought_class = dplyr::case_when(
      mean_spei >= DROUGHT_THR ~ "Non-drought",
      mean_spei < DROUGHT_THR & mean_spei >= -1.0 ~ "Mild",
      mean_spei < -1.0 & mean_spei >= -1.5 ~ "Moderate",
      mean_spei < -1.5 & mean_spei >= -2.0 ~ "Severe",
      mean_spei < -2.0 ~ "Extreme",
      TRUE ~ NA_character_
    )
  )

intensity_summary <- intensity_tbl |>
  dplyr::filter(is.finite(abs_delta_spei_pet), is.finite(mean_spei)) |>
  dplyr::group_by(scale, drought_class) |>
  dplyr::summarise(
    n_months = dplyr::n(),
    mean_abs_delta_spei = mean(abs_delta_spei_pet),
    median_abs_delta_spei = median(abs_delta_spei_pet),
    q25_abs_delta_spei = quantile(abs_delta_spei_pet, 0.25, names=FALSE),
    q75_abs_delta_spei = quantile(abs_delta_spei_pet, 0.75, names=FALSE),
    mean_drought_intensity = mean(drought_intensity),
    .groups = "drop"
  ) |>
  dplyr::mutate(drought_class = factor(drought_class,
                                       levels=c("Non-drought","Mild","Moderate","Severe","Extreme"))) |>
  dplyr::arrange(scale, drought_class)

intensity_association <- dplyr::bind_rows(lapply(split(intensity_tbl, intensity_tbl$scale), function(df) {
  ok <- is.finite(df$drought_intensity) & is.finite(df$abs_delta_spei_pet)
  data.frame(
    scale = unique(df$scale), n_months = sum(ok),
    pearson_r = if(sum(ok)>=3) cor(df$drought_intensity[ok], df$abs_delta_spei_pet[ok]) else NA_real_,
    spearman_rho = if(sum(ok)>=3) suppressWarnings(cor(df$drought_intensity[ok], df$abs_delta_spei_pet[ok], method="spearman")) else NA_real_
  )
}))

write.csv(intensity_summary, file.path(OUT_DIR, "PET_divergence_by_drought_intensity.csv"), row.names=FALSE)
write.csv(intensity_association, file.path(OUT_DIR, "PET_divergence_drought_intensity_association.csv"), row.names=FALSE)

extremes <- all_ts |>
  dplyr::filter(is.finite(abs_delta_spei_pet)) |>
  dplyr::group_by(scale) |>
  dplyr::arrange(dplyr::desc(abs_delta_spei_pet), .by_group = TRUE) |>
  dplyr::slice_head(n = 20L) |>
  dplyr::ungroup()

write.csv(extremes, file.path(OUT_DIR, "PET_divergence_extremes.csv"),
          row.names = FALSE)

plot_scales <- intersect(c(1L, 3L, 6L, 12L), unique(all_ts$scale))
pdat <- all_ts |>
  dplyr::filter(scale %in% plot_scales) |>
  tidyr::pivot_longer(c(SPEI_PM, SPEI_ERA5), names_to = "formulation",
                      values_to = "SPEI") |>
  dplyr::mutate(
    formulation = factor(formulation, c("SPEI_PM", "SPEI_ERA5"),
                         c("Penman-Monteith", "ERA5-Land native PEV")),
    scale_lab = factor(scale, levels = plot_scales,
                       labels = paste0("SPEI-", plot_scales))
  )

p1 <- ggplot2::ggplot(pdat, ggplot2::aes(date, SPEI, linetype = formulation)) +
  ggplot2::geom_hline(yintercept = DROUGHT_THR, linewidth = 0.35) +
  ggplot2::geom_line(linewidth = 0.35, na.rm = TRUE) +
  ggplot2::facet_wrap(~scale_lab, ncol = 1, scales = "free_x") +
  ggplot2::labs(x = NULL, y = "Basin-mean SPEI", linetype = "PET source/product",
                title = "SPEI sensitivity to PET source/product",
                subtitle = "PM and ERA5-Land native PEV are independently standardized") +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(legend.position = "bottom")

ggplot2::ggsave(file.path(OUT_DIR, "Fig_PET_source_sensitivity_timeseries.png"),
                p1, width = 11, height = 9, dpi = FIG_DPI)
ggplot2::ggsave(file.path(OUT_DIR, "Fig_PET_source_sensitivity_timeseries.pdf"),
                p1, width = 11, height = 9)

mplot <- metrics |>
  dplyr::filter(period %in% c("Full record", "2022-2025",
                              "Warm season (MJJAS)", "Drought-active")) |>
  dplyr::select(scale, period, pearson_r, rmse,
                drought_state_agreement_pct, mean_abs_delta_spei) |>
  tidyr::pivot_longer(-c(scale, period), names_to = "metric", values_to = "value") |>
  dplyr::mutate(
    metric = factor(metric,
                    levels = c("pearson_r", "rmse", "drought_state_agreement_pct",
                               "mean_abs_delta_spei"),
                    labels = c("Pearson r", "RMSE", "Drought-state agreement (%)",
                               "Mean |Delta SPEI|"))
  )

p2 <- ggplot2::ggplot(mplot, ggplot2::aes(factor(scale), value, group = period,
                                          linetype = period)) +
  ggplot2::geom_point() +
  ggplot2::geom_line() +
  ggplot2::facet_wrap(~metric, scales = "free_y", ncol = 2) +
  ggplot2::labs(x = "Accumulation scale (months)", y = NULL,
                linetype = "Period",
                title = "PET-source/product agreement metrics") +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(legend.position = "bottom")

ggplot2::ggsave(file.path(OUT_DIR, "Fig_PET_source_agreement_metrics.png"),
                p2, width = 10, height = 7, dpi = FIG_DPI)
ggplot2::ggsave(file.path(OUT_DIR, "Fig_PET_source_agreement_metrics.pdf"),
                p2, width = 10, height = 7)

cat("\n12PET_source_sensitivity_v2.R completed.\n")
cat("Outputs:", OUT_DIR, "\n")
