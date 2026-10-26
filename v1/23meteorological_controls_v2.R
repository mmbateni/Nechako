################################################################################
# 13meteorological_controls_v2.R
# METEOROLOGICAL CONTROLS OF PET-FORMULATION SENSITIVITY
# CURRENT MANUSCRIPT RESULTS (Methods/Results 3f/4f)
#
# Purpose
#   Diagnose meteorological conditions associated with large |delta SPEI_PET|.
#   This is association/context, NOT additive physical attribution.
#
# Inputs
#   monthly_data_direct/meteorological_controls_monthly_summary.csv
#     Required: date, year, month, t2m_degC, rn_MJ_m2_d, vpd_kPa, wind2m_ms
#     Optional: pressure_kPa
#   pet_sensitivity/PET_sensitivity_monthly.csv
#
# Standardization
#   Each meteorological control is standardized relative to its own calendar-
#   month climatology over 1950-2025, preventing the seasonal cycle from
#   dominating the association analysis.
#
# Outputs -> met_controls/
#   meteorological_controls_standardized.csv
#   met_anomaly_correlations.csv
#   met_anomaly_composites_warm_season.csv
#   met_large_divergence_months.csv
#   met_multivariable_conditioning.csv
#   Fig_met_controls_large_divergence.pdf/.png
#
# MODIFICATION
#   pressure_kPa is OPTIONAL. The four manuscript-facing controls remain the
#   primary association/composite variables: temperature, net radiation, VPD,
#   and 2-m wind speed. Surface pressure is included in the multivariable
#   conditioning model only when the upstream PET summary provides it.
#
# INTERPRETATION
#   The analysis does not decompose delta SPEI_PET into physical drivers.
#   Predictors are statistical conditioning variables for formulation
#   divergence. Correlation and regression coefficients are associations.
################################################################################

setwd(Sys.getenv("NECHAKO_WD", "D:/Nechako_Drought/Nechako/"))
source("DROUGHT_ANALYSIS_utils_v2.R")
utils_load_packages(c("dplyr", "tidyr", "ggplot2", "data.table", "lubridate"))

if (!dir.exists(WD_PATH)) stop("Working directory not found: ", WD_PATH)
setwd(WD_PATH)

OUT_DIR <- file.path(WD_PATH, "met_controls")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

MET_FILE <- file.path(
  WD_PATH, "monthly_data_direct",
  "meteorological_controls_monthly_summary.csv"
)
SENS_FILE <- file.path(
  WD_PATH, "pet_sensitivity",
  "PET_sensitivity_monthly.csv"
)

if (!file.exists(MET_FILE))
  stop("Missing ", MET_FILE, ". Re-run 2b_PET_ERALand_v2.R.")
if (!file.exists(SENS_FILE))
  stop("Missing ", SENS_FILE, ". Run 12PET_formulation_sensitivity_v2.R first.")

met <- data.table::fread(MET_FILE)
sens <- data.table::fread(SENS_FILE)

################################################################################
# 1. INPUT VALIDATION
################################################################################

required_met <- c(
  "date", "year", "month", "t2m_degC", "rn_MJ_m2_d",
  "vpd_kPa", "wind2m_ms"
)

miss <- setdiff(required_met, names(met))
if (length(miss)) {
  stop(
    "Meteorological summary missing required core column(s): ",
    paste(miss, collapse = ", ")
  )
}

required_sens <- c(
  "date", "scale", "SPEI_PM", "SPEI_THW",
  "delta_spei_pet", "abs_delta_spei_pet", "drought_active"
)

miss_sens <- setdiff(required_sens, names(sens))
if (length(miss_sens)) {
  stop(
    "PET-sensitivity file missing required column(s): ",
    paste(miss_sens, collapse = ", "),
    ". Re-run 12PET_formulation_sensitivity_v2.R."
  )
}

normalise_month_date <- function(x, label = "date") {
  x <- trimws(as.character(x))
  if (!length(x)) stop(label, " is empty.")
  
  is_ym  <- grepl("^\\d{4}-\\d{2}$", x)
  is_ymd <- grepl("^\\d{4}-\\d{2}-\\d{2}$", x)
  unsupported <- !(is_ym | is_ymd)
  
  if (any(unsupported)) {
    stop(
      label, " contains unsupported date format(s): ",
      paste(head(unique(x[unsupported]), 10L), collapse = ", "),
      ". Expected YYYY-MM or YYYY-MM-DD."
    )
  }
  
  x_canonical <- ifelse(is_ym, paste0(x, "-01"), x)
  out <- as.Date(x_canonical, format = "%Y-%m-%d")
  
  if (anyNA(out)) {
    stop(
      label, " contains invalid calendar date(s): ",
      paste(head(unique(x[is.na(out)]), 10L), collapse = ", ")
    )
  }
  
  as.Date(format(out, "%Y-%m-01"), format = "%Y-%m-%d")
}

met$date  <- normalise_month_date(met$date, "meteorological-controls date")
sens$date <- normalise_month_date(sens$date, "PET-sensitivity date")
met$month <- as.integer(met$month)
met$year  <- as.integer(met$year)
sens$scale <- as.integer(sens$scale)

validate_monthly_record(
  met$date,
  as.Date("1950-01-01"),
  as.Date("2025-12-01"),
  label = "meteorological-controls record",
  require_complete = TRUE
)

if (any(!is.finite(met$month)) || any(!met$month %in% 1:12))
  stop("Meteorological summary contains invalid calendar-month values.")

################################################################################
# 2. CORE METEOROLOGICAL CONTROLS
################################################################################

# These four variables are the manuscript-facing PM controls. Surface pressure
# is deliberately NOT added to this primary association/composite set because
# the current Results-3f framework identifies temperature, radiation, humidity/
# VPD and wind as the meteorological controls represented explicitly in the PM
# formulation beyond Thornthwaite's temperature-only structure.
driver_map <- c(
  t2m_degC   = "Temperature",
  rn_MJ_m2_d = "Net radiation",
  vpd_kPa    = "VPD",
  wind2m_ms  = "Wind speed"
)

z_by_calendar_month <- function(x, month) {
  x <- as.numeric(x)
  month <- as.integer(month)
  out <- rep(NA_real_, length(x))
  
  for (m in 1:12) {
    ii <- which(month == m & is.finite(x))
    if (!length(ii)) next
    
    mu <- mean(x[ii])
    sig <- stats::sd(x[ii])
    
    if (!is.finite(sig) || sig <= 0) {
      out[ii] <- NA_real_
    } else {
      out[ii] <- (x[ii] - mu) / sig
    }
  }
  
  out
}

for (v in names(driver_map)) {
  met[[paste0(v, "_z")]] <- z_by_calendar_month(
    met[[v]],
    met$month
  )
}

# Optional pressure branch. It is retained for statistical conditioning but is
# not promoted to the primary four-driver composite analysis.
has_pressure <- "pressure_kPa" %in% names(met)

if (has_pressure) {
  met$pressure_kPa_z <- z_by_calendar_month(
    met$pressure_kPa,
    met$month
  )
  cat("  Optional pressure predictor detected: pressure_kPa\n")
} else {
  cat("  Optional pressure predictor not present; continuing with core four controls.\n")
}

write.csv(
  met,
  file.path(OUT_DIR, "meteorological_controls_standardized.csv"),
  row.names = FALSE
)

################################################################################
# 3. JOIN PET-FORMULATION SENSITIVITY AND METEOROLOGICAL CONTROLS
################################################################################

sens_keep <- sens[, required_sens, with = FALSE]
met_join <- met

joined <- merge(
  sens_keep,
  met_join,
  by = "date",
  all.x = TRUE,
  sort = TRUE
)

if (anyNA(joined$month)) {
  stop("PET-sensitivity dates failed to match meteorological-controls dates.")
}

if (anyNA(joined$scale)) {
  stop("PET-sensitivity records contain missing scale values.")
}

joined$warm_season <- joined$month %in% 5:9
joined$recent <- joined$date >= as.Date("2022-01-01") &
  joined$date <= as.Date("2025-12-01")

# Primary association predictors.
driver_z <- paste0(names(driver_map), "_z")

################################################################################
# 4. BIVARIATE ASSOCIATION ANALYSIS
################################################################################

cor_one <- function(df, zcol, scale, period) {
  ok <- is.finite(df[[zcol]]) & is.finite(df$abs_delta_spei_pet)
  x <- df[[zcol]][ok]
  y <- df$abs_delta_spei_pet[ok]
  n <- length(x)
  
  pear <- if (
    n >= 3L &&
    is.finite(stats::sd(x)) &&
    is.finite(stats::sd(y)) &&
    stats::sd(x) > 0 &&
    stats::sd(y) > 0
  ) {
    tryCatch(
      stats::cor.test(x, y, method = "pearson"),
      error = function(e) NULL
    )
  } else {
    NULL
  }
  
  spear <- if (n >= 3L) {
    tryCatch(
      suppressWarnings(
        stats::cor.test(x, y, method = "spearman", exact = FALSE)
      ),
      error = function(e) NULL
    )
  } else {
    NULL
  }
  
  data.frame(
    scale = as.integer(scale),
    period = period,
    driver = unname(driver_map[sub("_z$", "", zcol)]),
    n = n,
    pearson_r = if (is.null(pear)) NA_real_ else unname(pear$estimate),
    pearson_p = if (is.null(pear)) NA_real_ else pear$p.value,
    spearman_rho = if (is.null(spear)) NA_real_ else unname(spear$estimate),
    spearman_p = if (is.null(spear)) NA_real_ else spear$p.value,
    stringsAsFactors = FALSE
  )
}

correlations <- dplyr::bind_rows(
  lapply(split(joined, joined$scale), function(df) {
    periods <- list(
      "Full record" = df,
      "Warm season (MJJAS)" = df[df$warm_season, , drop = FALSE],
      "2022-2025" = df[df$recent, , drop = FALSE],
      "2022-2025 warm season" =
        df[df$recent & df$warm_season, , drop = FALSE],
      "Drought-active" =
        df[df$drought_active %in% TRUE, , drop = FALSE]
    )
    
    dplyr::bind_rows(
      lapply(names(periods), function(pn) {
        dplyr::bind_rows(
          lapply(
            driver_z,
            function(zc) cor_one(periods[[pn]], zc, unique(df$scale), pn)
          )
        )
      })
    )
  })
)

correlations <- correlations |>
  dplyr::group_by(scale, period) |>
  dplyr::mutate(
    pearson_p_fdr = stats::p.adjust(pearson_p, method = "BH"),
    spearman_p_fdr = stats::p.adjust(spearman_p, method = "BH")
  ) |>
  dplyr::ungroup()

write.csv(
  correlations,
  file.path(OUT_DIR, "met_anomaly_correlations.csv"),
  row.names = FALSE
)

################################################################################
# 5. MULTIVARIABLE STATISTICAL CONDITIONING
################################################################################

# This is a conditioning analysis only. It does not assign physical portions
# of PM PET or SPEI divergence to individual predictors.
BOOT_R <- 1000L
BOOT_BLOCK <- 12L
BOOT_SEED <- 13013L
set.seed(BOOT_SEED)

model_predictors <- driver_z

if (has_pressure) {
  model_predictors <- c(model_predictors, "pressure_kPa_z")
}

predictor_labels <- c(
  t2m_degC_z = "Temperature",
  rn_MJ_m2_d_z = "Net radiation",
  vpd_kPa_z = "VPD",
  wind2m_ms_z = "Wind speed",
  pressure_kPa_z = "Surface pressure"
)

block_bootstrap_rows <- function(n, block_length = BOOT_BLOCK) {
  n <- as.integer(n)
  block_length <- as.integer(block_length)
  
  if (n <= 1L) return(seq_len(n))
  if (block_length < 1L) stop("block_length must be >= 1.")
  
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

vif_from_matrix <- function(x) {
  x <- as.matrix(x)
  if (ncol(x) < 2L) {
    return(setNames(rep(NA_real_, ncol(x)), colnames(x)))
  }
  
  out <- setNames(rep(NA_real_, ncol(x)), colnames(x))
  
  for (j in seq_len(ncol(x))) {
    others <- x[, -j, drop = FALSE]
    
    fit <- tryCatch(
      stats::lm(x[, j] ~ others),
      error = function(e) NULL
    )
    
    r2 <- if (is.null(fit)) NA_real_ else summary(fit)$r.squared
    
    out[j] <- if (
      is.finite(r2) && r2 < 1
    ) {
      1 / (1 - r2)
    } else if (isTRUE(all.equal(r2, 1))) {
      Inf
    } else {
      NA_real_
    }
  }
  
  out
}

fit_multivariable_one <- function(
    df,
    scale,
    period,
    B = BOOT_R,
    block_length = BOOT_BLOCK) {
  
  # IMPORTANT:
  # split(joined, joined$scale) returns data.table objects because joined is
  # created with data.table::fread()/merge(). Base-style row indexing such as
  # d[ii, , drop = FALSE] can then be interpreted by [.data.table as a join,
  # producing:
  #   Incompatible join types: x.date (...) and i.V1 (character)
  # Convert the per-period object to a plain data.frame before any positional
  # bootstrap indexing. This is a type-safety fix; it does not change the model.
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  
  keep <- is.finite(df$abs_delta_spei_pet) &
    Reduce(`&`, lapply(df[model_predictors], is.finite))
  
  d <- df[
    keep,
    c("date", "abs_delta_spei_pet", model_predictors),
    drop = FALSE
  ]
  d <- as.data.frame(d, stringsAsFactors = FALSE)
  
  n <- nrow(d)
  p <- length(model_predictors)
  
  empty_result <- function(reps = 0L) {
    data.frame(
      scale = as.integer(scale),
      period = period,
      n = n,
      predictor = model_predictors,
      beta = NA_real_,
      beta_ci_low = NA_real_,
      beta_ci_high = NA_real_,
      model_r2 = NA_real_,
      model_adj_r2 = NA_real_,
      vif = NA_real_,
      bootstrap_reps = as.integer(reps),
      block_length = as.integer(min(block_length, max(1L, n))),
      stringsAsFactors = FALSE
    )
  }
  
  if (n < max(20L, p + 5L)) {
    return(empty_result(0L))
  }
  
  form <- stats::as.formula(
    paste(
      "abs_delta_spei_pet ~",
      paste(model_predictors, collapse = "+")
    )
  )
  
  fit <- tryCatch(
    stats::lm(form, data = d),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(empty_result(0L))
  }
  
  co <- stats::coef(fit)
  sm <- summary(fit)
  vif <- vif_from_matrix(as.matrix(d[, model_predictors, drop = FALSE]))
  
  bb <- matrix(
    NA_real_,
    nrow = B,
    ncol = p,
    dimnames = list(NULL, model_predictors)
  )
  
  for (b in seq_len(B)) {
    ii <- block_bootstrap_rows(n, block_length)
    
    fb <- tryCatch(
      stats::lm(form, data = d[ii, , drop = FALSE]),
      error = function(e) NULL
    )
    
    if (!is.null(fb)) {
      bc <- stats::coef(fb)
      bb[b, ] <- unname(bc[model_predictors])
    }
  }
  
  ci <- lapply(seq_len(p), function(j) {
    z <- bb[, j]
    z <- z[is.finite(z)]
    
    if (!length(z)) {
      c(NA_real_, NA_real_)
    } else {
      as.numeric(
        stats::quantile(
          z,
          probs = c(0.025, 0.975),
          names = FALSE,
          na.rm = TRUE
        )
      )
    }
  })
  
  data.frame(
    scale = as.integer(scale),
    period = period,
    n = n,
    predictor = model_predictors,
    beta = unname(co[model_predictors]),
    beta_ci_low = vapply(ci, `[[`, numeric(1), 1L),
    beta_ci_high = vapply(ci, `[[`, numeric(1), 2L),
    model_r2 = sm$r.squared,
    model_adj_r2 = sm$adj.r.squared,
    vif = unname(vif[model_predictors]),
    bootstrap_reps = as.integer(B),
    block_length = as.integer(min(block_length, n)),
    stringsAsFactors = FALSE
  )
}

multivariable_results <- dplyr::bind_rows(
  lapply(split(joined, joined$scale), function(df) {
    periods <- list(
      "Full record" = df,
      "Warm season (MJJAS)" = df[df$warm_season, , drop = FALSE],
      "2022-2025" = df[df$recent, , drop = FALSE],
      "2022-2025 warm season" =
        df[df$recent & df$warm_season, , drop = FALSE],
      "Drought-active" =
        df[df$drought_active %in% TRUE, , drop = FALSE]
    )
    
    dplyr::bind_rows(
      lapply(
        names(periods),
        function(pn) {
          fit_multivariable_one(
            periods[[pn]],
            unique(df$scale),
            pn
          )
        }
      )
    )
  })
)

multivariable_results$predictor <- factor(
  multivariable_results$predictor,
  levels = model_predictors,
  labels = unname(predictor_labels[model_predictors])
)

write.csv(
  multivariable_results,
  file.path(OUT_DIR, "met_multivariable_conditioning.csv"),
  row.names = FALSE
)

################################################################################
# 6. SCALE-SPECIFIC LARGE-DIVERGENCE THRESHOLD
################################################################################

# Define large divergence independently within each scale from the upper
# quartile of |delta SPEI| over the full record. This avoids imposing an
# arbitrary standardized-index difference threshold.
thresholds <- joined |>
  dplyr::group_by(scale) |>
  dplyr::summarise(
    abs_delta_q75 = stats::quantile(
      abs_delta_spei_pet,
      0.75,
      na.rm = TRUE,
      names = FALSE
    ),
    .groups = "drop"
  )

joined <- dplyr::left_join(joined, thresholds, by = "scale")
joined$large_divergence <- is.finite(joined$abs_delta_spei_pet) &
  is.finite(joined$abs_delta_q75) &
  joined$abs_delta_spei_pet >= joined$abs_delta_q75

large_recent <- joined |>
  dplyr::filter(recent, large_divergence) |>
  dplyr::arrange(scale, dplyr::desc(abs_delta_spei_pet))

write.csv(
  large_recent,
  file.path(OUT_DIR, "met_large_divergence_months.csv"),
  row.names = FALSE
)

################################################################################
# 7. WARM-SEASON COMPOSITES
################################################################################

long <- joined |>
  dplyr::filter(warm_season) |>
  tidyr::pivot_longer(
    dplyr::all_of(driver_z),
    names_to = "driver_code",
    values_to = "anomaly_z"
  ) |>
  dplyr::mutate(
    driver = unname(driver_map[sub("_z$", "", driver_code)]),
    divergence_class = ifelse(
      large_divergence,
      "Large PET-formulation divergence",
      "Other months"
    )
  )

composites <- long |>
  dplyr::group_by(scale, driver, divergence_class) |>
  dplyr::summarise(
    n = sum(is.finite(anomaly_z)),
    mean_anomaly_z = if (any(is.finite(anomaly_z)))
      mean(anomaly_z, na.rm = TRUE) else NA_real_,
    median_anomaly_z = if (any(is.finite(anomaly_z)))
      stats::median(anomaly_z, na.rm = TRUE) else NA_real_,
    sd_anomaly_z = if (sum(is.finite(anomaly_z)) >= 2L)
      stats::sd(anomaly_z, na.rm = TRUE) else NA_real_,
    .groups = "drop"
  )

write.csv(
  composites,
  file.path(OUT_DIR, "met_anomaly_composites_warm_season.csv"),
  row.names = FALSE
)

################################################################################
# 8. MANUSCRIPT-STYLE COMPOSITE FIGURE
#
# IMPORTANT:
#   This figure intentionally uses the FULL 1950-2025 record.
#   The large-divergence threshold is defined from the full-record upper
#   quartile of |Delta SPEI_PET| within each scale.
#
#   The plotted observations are warm-season (MJJAS) months classified as
#   large-divergence months. This is distinct from the separate 2022-2025
#   recent-period association analysis.
################################################################################

plot_scales <- intersect(
  c(1L, 3L, 6L, 12L, 24L, 36L),
  sort(unique(as.integer(long$scale)))
)

plot_dat <- long |>
  dplyr::filter(
    scale %in% plot_scales,
    divergence_class == "Large PET-formulation divergence"
  ) |>
  dplyr::group_by(scale, driver) |>
  dplyr::summarise(
    n = sum(is.finite(anomaly_z)),
    mean_z = if (n > 0L) {
      mean(anomaly_z, na.rm = TRUE)
    } else {
      NA_real_
    },
    se_z = if (n >= 2L) {
      stats::sd(anomaly_z, na.rm = TRUE) / sqrt(n)
    } else {
      NA_real_
    },
    .groups = "drop"
  ) |>
  dplyr::mutate(
    scale_label = factor(
      paste0("SPEI-", scale),
      levels = paste0("SPEI-", plot_scales)
    ),
    driver = factor(
      driver,
      levels = c(
        "Temperature",
        "Net radiation",
        "VPD",
        "Wind speed"
      )
    )
  )

# Colorblind-friendly Okabe-Ito-inspired palette.
driver_cols <- c(
  "Temperature"   = "#D55E00",
  "Net radiation" = "#0072B2",
  "VPD"            = "#009E73",
  "Wind speed"    = "#CC79A7"
)

driver_linetypes <- c(
  "Temperature"   = "solid",
  "Net radiation" = "dashed",
  "VPD"           = "dotdash",
  "Wind speed"    = "dotted"
)

driver_shapes <- c(
  "Temperature"   = 16,
  "Net radiation" = 17,
  "VPD"           = 15,
  "Wind speed"    = 18
)

if (nrow(plot_dat) > 0L) {
  
  p <- ggplot2::ggplot(
    plot_dat,
    ggplot2::aes(
      x = scale_label,
      y = mean_z,
      group = driver,
      colour = driver,
      linetype = driver
    )
  ) +
    
    # Reference climatological state.
    ggplot2::geom_hline(
      yintercept = 0,
      linewidth = 0.55,
      colour = "grey45"
    ) +
    
    # Standard-error intervals.
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = mean_z - se_z,
        ymax = mean_z + se_z
      ),
      width = 0.10,
      linewidth = 0.65,
      alpha = 0.75,
      na.rm = TRUE
    ) +
    
    # Connecting lines.
    ggplot2::geom_line(
      linewidth = 1.00,
      lineend = "round",
      na.rm = TRUE
    ) +
    
    # Point estimates.
    ggplot2::geom_point(
      ggplot2::aes(shape = driver),
      size = 3.1,
      stroke = 0.25,
      na.rm = TRUE
    ) +
    
    ggplot2::scale_colour_manual(
      values = driver_cols,
      name = "Meteorological control"
    ) +
    
    ggplot2::scale_linetype_manual(
      values = driver_linetypes,
      name = "Meteorological control"
    ) +
    
    ggplot2::scale_shape_manual(
      values = driver_shapes,
      name = "Meteorological control"
    ) +
    
    ggplot2::scale_y_continuous(
      name = "Calendar-month standardized anomaly (z)",
      breaks = scales::pretty_breaks(n = 6),
      expand = ggplot2::expansion(mult = c(0.05, 0.08))
    ) +
    
    ggplot2::scale_x_discrete(
      name = "SPEI accumulation scale (months)"
    ) +
    
    ggplot2::labs(
      title = "Meteorological conditions during large PET-formulation divergence",
      subtitle = paste0(
        "Warm-season (MJJAS) large-divergence months, 1950–2025; ",
        "threshold = scale-specific upper quartile of |ΔSPEI\u209a\u2091\u209c|"
      ),
      caption = paste0(
        "Points show mean calendar-month standardized anomalies; ",
        "error bars denote ±1 SE. Anomalies are standardized separately ",
        "for each calendar month over 1950–2025."
      )
    ) +
    
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        order = 1,
        override.aes = list(
          linewidth = 1.1,
          size = 3
        )
      ),
      linetype = ggplot2::guide_legend(
        order = 1
      ),
      shape = ggplot2::guide_legend(
        order = 1
      )
    ) +
    
    ggplot2::theme_minimal(
      base_size = 11,
      base_family = "sans"
    ) +
    
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 14,
        colour = "#202020",
        margin = ggplot2::margin(b = 5)
      ),
      
      plot.subtitle = ggplot2::element_text(
        size = 10.5,
        colour = "#505050",
        margin = ggplot2::margin(b = 10)
      ),
      
      plot.caption = ggplot2::element_text(
        size = 8.5,
        colour = "#606060",
        hjust = 0,
        margin = ggplot2::margin(t = 10)
      ),
      
      axis.title.x = ggplot2::element_text(
        face = "bold",
        size = 11,
        margin = ggplot2::margin(t = 8)
      ),
      
      axis.title.y = ggplot2::element_text(
        face = "bold",
        size = 11,
        margin = ggplot2::margin(r = 8)
      ),
      
      axis.text.x = ggplot2::element_text(
        size = 10.5,
        face = "bold",
        colour = "#303030"
      ),
      
      axis.text.y = ggplot2::element_text(
        size = 10,
        colour = "#303030"
      ),
      
      panel.grid.major.x = ggplot2::element_blank(),
      
      panel.grid.major.y = ggplot2::element_line(
        colour = "#D9D9D9",
        linewidth = 0.35
      ),
      
      panel.grid.minor = ggplot2::element_blank(),
      
      panel.border = ggplot2::element_rect(
        colour = "#B5B5B5",
        fill = NA,
        linewidth = 0.55
      ),
      
      legend.position = "bottom",
      
      legend.title = ggplot2::element_text(
        face = "bold",
        size = 10
      ),
      
      legend.text = ggplot2::element_text(
        size = 9.5
      ),
      
      legend.key.width = grid::unit(1.4, "cm"),
      
      legend.spacing.x = grid::unit(0.35, "cm"),
      
      plot.margin = ggplot2::margin(
        t = 8,
        r = 10,
        b = 8,
        l = 10
      )
    )
  
  ggplot2::ggsave(
    file.path(
      OUT_DIR,
      "Fig_met_controls_large_divergence.png"
    ),
    p,
    width = 9.2,
    height = 6.2,
    dpi = 600,
    bg = "white"
  )
  
  ggplot2::ggsave(
    file.path(
      OUT_DIR,
      "Fig_met_controls_large_divergence.pdf"
    ),
    p,
    width = 9.2,
    height = 6.2,
    device = cairo_pdf,
    bg = "white"
  )
  
} else {
  
  warning(
    paste(
      "No warm-season large-divergence observations available",
      "for the requested figure scales."
    ),
    call. = FALSE
  )
}
################################################################################
# 9. FINAL SUMMARY
################################################################################

cat("\n13meteorological_controls_v2.R completed.\n")
cat("Outputs: ", normalizePath(OUT_DIR, mustWork = FALSE), "\n", sep = "")
cat("Core association predictors: Temperature, Net radiation, VPD, Wind speed\n")
cat(
  "Multivariable conditioning predictors: ",
  paste(
    unname(predictor_labels[model_predictors]),
    collapse = ", "
  ),
  "\n",
  sep = ""
)
if (!has_pressure) {
  cat("Surface pressure conditioning: unavailable in meteorological summary; core analysis retained.\n")
} else {
  cat("Surface pressure conditioning: included as an optional statistical predictor.\n")
}
cat("Interpretation: association/context only; no additive physical attribution.\n")
