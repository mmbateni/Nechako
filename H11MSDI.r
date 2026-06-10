# ==============================================================================
# NECHAKO BASIN — VINE COPULA MULTIVARIATE STANDARDISED DROUGHT INDEX (MSDI)
# Enhanced version — R-vine & D-vine with seasonal stratification,
# reference-period-calibrated margins, proper Kendall distribution via
# Monte-Carlo, tail-dependence diagnostics, vine truncation comparison,
# drought-event ranking, component contribution, and publication plots.
#
# COMPONENT INDICES (5)
# ─────────────────────────────────────────────────────────────────────────────
#  Variable    Index      Script source              Hydrological role
#  ─────────── ─────────  ─────────────────────────  ───────────────────────
#  SWEI        Snow       8SWEI_Snow.R               SWE / snowpack state
#  SSI_soil    Soil       3SSI_SSMI_ERALand.R        Soil-moisture deficit
#  SGI         Ground     10SGI.R                    Groundwater level
#  SSI_wsc     Stream     7WSC_StreamFlow.R          Streamflow (WSC gauges)
#  NTSDI       Reservoir  9NTSDI_GeoLakes_Extended.r Nechako Reservoir storage
#
# PHYSICAL CASCADE — PARTIAL (see MANAGED RESERVOIR NOTE below)
# ─────────────────────────────────────────────────────────────────────────────
#  NATURAL sub-cascade (unambiguous, used as D-vine spine):
#    SWEI → SSI_soil → SGI
#    Snowpack accumulation → spring/summer soil recharge from melt →
#    percolation to groundwater storage
#
#  MANAGED variables (NO imposed cascade order — R-vine governs):
#    NTSDI  — Nechako Reservoir (Kenney Dam, 1952; operated by BC Hydro)
#             Level = f(natural inflow) MINUS f(Kemano diversion, spill rules,
#             power-demand scheduling).  Operationally decoupled from SGI.
#    SSI_wsc — WSC gauges DOWNSTREAM of Kenney Dam (e.g. 08JC001 Nechako
#             at Vanderhoof) measure REGULATED RELEASES, not natural runoff.
#             Release timing is driven by BC Hydro power demand — NOT by
#             reservoir storage in a simple physical sense.
#
#  ⚠  VINE STRUCTURE RECOMMENDATION
#     D-vine: spine restricted to the 3 natural variables; NTSDI and SSI_wsc
#             appended at the tail with order chosen by BIC (data-driven).
#     R-vine: PRIMARY model — data-driven MST structure imposes no physical
#             cascade assumption and is robust to operational management.
#
#  ⚠ SEASONAL STRATIFICATION (OPTION A — EXCLUSIVE)
#     Snow season (Oct–Apr):  all 5 variables — fit vines_snow5
#     Warm season (May–Sep):  4 variables — SWEI excluded, fit vines_warm4
#     Full-year vine REMOVED — structurally incompatible with seasonal SWEI
#
# METHOD ENHANCEMENTS OVER DRAFT
# ─────────────────────────────────────────────────────────────────────────────
#  • Reference-period-calibrated empirical marginal CDFs (per calendar month)
#  • Kendall distribution estimated via MC simulation from the fitted vine
#    (vs. simple Blom ranking on the observed sample)
#  • Vine truncation analysis (Tree 1 / 2 / full) with BIC comparison
#  • Empirical tail-dependence coefficients for all pairs (λ_L, λ_U)
#  • Conditional contribution scores — which variable drives drought
#  • Bootstrap 90 % confidence bands on MSDI (optional, toggle in cfg)
#  • Drought event table: duration, severity (deficit sum), peak MSDI
#  • Extended plot suite: heatmap, contribution chart, bivariate contours
#
# REFERENCES
#   Sklar (1959); Joe (2015) Dependence Modelling with Copulas
#   Aas et al. (2009) Ins. Math. Econ. 44, 182–198 — pair-copula construction
#   Bedford & Cooke (2002) Ann. Stat. 30, 1031–1068 — R-vine
#   Genest & Favre (2007) J. Hydrol. Eng. 12, 347–368
#   Farahmand & AghaKouchak (2015) Adv. Water Resour. 76, 140–145
#   Bloomfield & Marchant (2013) HESS 17, 4769–4787 (SGI)
#   Huning & AghaKouchak (2020) Nat. Commun. 11, 965 (SWEI)
# ==============================================================================

rm(list = ls())

# ==============================================================================
# 1.  PACKAGES
# ==============================================================================
required_pkgs <- c(
  "ggplot2","dplyr",
  "tidyverse", "lubridate",
  "rvinecopulib",   # vine fitting, CDF, simulation
  "VineCopula"      # GOF, Rosenblatt transform
)

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cran.rstudio.com/")
  library(pkg, character.only = TRUE)
}

# ==============================================================================
# 2.  USER CONFIGURATION
# ==============================================================================
cfg <- list(
  
  # ── Input file paths (match output conventions of upstream scripts) ─────────
  
  # SWEI  — wide CSV:  Year column + 3-letter month columns (Jan…Dec)
  swei_file     = "swei_results_seasonal/swei_01_basin_averaged_by_month.csv",
  
  # SSI soil  — per-month CSVs assembled by assemble_ssi_from_monthly_csvs().
  # Set folder and prefix to match 3SSI_SSMI_ERALand.R output naming exactly
  # (prefix is case-sensitive on Linux/macOS — verify against actual filenames).
  ssi_soil_folder = "ssi_results_seasonal",
  ssi_soil_prefix = "ssi_L1_3_01",
  
  # SGI  — long CSV:  columns 'month_label' (YYYY-MM), 'sgi_median', 'date'
  sgi_file = "sgi_nechako_multisource.csv",
  
  # WSC streamflow SSI — RDS files saved by 7WSC_StreamFlow.R
  # Each RDS is a data frame with columns: date (YYYY-MM-15), ssi, ...
  ssi_wsc_rds = c(
    "streamflow_results/ssi/monthly_data/08JE001_monthly_ssi.rds",
    "streamflow_results/ssi/monthly_data/08JC002_monthly_ssi.rds"
  ),
  # NTSDI — long CSV from 9NTSDI_GeoLakes_Extended.r Pipeline 2
  # Columns: date, NTSDI_norm  (Blom-standardised reservoir index)
  ntsdi_file    = "Nechako_NTSDI_Pipeline2_Output.csv",
  
  # ── Reference period (calibrate pseudo-obs marginals here only) ────────────
  ref_start     = as.Date("1981-01-01"),
  ref_end       = as.Date("2020-12-31"),
  
  # ── Seasonal stratification ─────────────────────────────────────────────────
  # snow_months  — months with reliable SWE signal; SWEI is valid here
  # warm_months  — SWEI excluded (or replaced by SSPI-1 if snow_index = "sspi1")
  snow_months   = c(10, 11, 12, 1, 2, 3, 4),   # Oct–Apr
  warm_months   = c(5, 6, 7, 8, 9),             # May–Sep
  
  # ── Snow index choice ──────────────────────────────────────────────────────
  # "swei"  — use SWEI (NAs in warm season → seasonal vines)
  # "sspi1" — use SSPI-1 year-round (from sspi_results_monthly/ folder)
  #            set sspi1_file path below if choosing "sspi1"
  snow_index    = "swei",
  sspi1_file    = "sspi_results_monthly/sspi_01_basin_averaged_by_month.csv",
  
  # ── Vine copula settings ────────────────────────────────────────────────────
  family_set    = "parametric",   # "tll" for nonparametric (slow, flexible)
  selcrit       = "bic",          # "aic" for larger models
  par_method    = "mle",
  
  # D-vine ORDER RATIONALE
  # ─────────────────────────────────────────────────────────────────────────
  # Positions 1-3: NATURAL cascade (physically unambiguous)
  #   Snow(1) → SSI_soil(2) → SGI(3)
  #
  # Positions 4-5: MANAGED variables — no physical cascade can be assumed.
  #   The Nechako Reservoir is operated for hydropower (Kemano diversion).
  #   Downstream WSC flow is a regulated release, not a natural runoff signal.
  #   TWO D-vine tail orderings are tried; BIC selects between them:
  #     Variant A: … → NTSDI(4) → SSI_wsc(5)
  #     Variant B: … → SSI_wsc(4) → NTSDI(5)
  #
  # var_list column index mapping (used in dvine_order_* below):
  #   Snow=1  SSI_soil=2  SGI=3  NTSDI=4  SSI_wsc=5
  dvine_order_5A = c(1L, 2L, 3L, 4L, 5L),   # Variant A: natural → NTSDI → SSI_wsc
  dvine_order_5B = c(1L, 2L, 3L, 5L, 4L),   # Variant B: natural → SSI_wsc → NTSDI
  
  # Warm-season 4-variable vine (Snow excluded); same two tail variants:
  #   SSI_soil=1  SGI=2  NTSDI=3  SSI_wsc=4   (positions within the 4-var submatrix)
  dvine_order_4A = c(1L, 2L, 3L, 4L),        # SSI_soil → SGI → NTSDI → SSI_wsc
  dvine_order_4B = c(1L, 2L, 4L, 3L),        # SSI_soil → SGI → SSI_wsc → NTSDI
  
  # Vine truncation levels to compare (NULL = fit full vine only)
  truncation_levels = c(1L, 2L),
  
  # Monte-Carlo samples for vine CDF integration (pvinecop) and Kendall dist.
  n_mc          = 2e4L,           # increase to 1e5 for publication
  n_kc_sim      = 1e5L,           # MC draws for Kendall distribution (1e5 publication; 5e4 for speed)
  
  # ── Bootstrap CI (set to TRUE for publication-grade output) ────────────────
  run_bootstrap = FALSE,
  n_boot        = 200L,
  
  # ── Output ──────────────────────────────────────────────────────────────────
  out_dir       = "vine_msdi_output"
)

dir.create(cfg$out_dir, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 3.  HELPER FUNCTIONS
# ==============================================================================
plot_standard_basin_ts <- function(df, index_name, color = "steelblue", drought_threshold = -0.5) {
  if (!all(c("date", "value") %in% names(df))) stop("df must have 'date' and 'value' columns.")
  df <- df[order(df$date), ]
  df <- df[!is.na(df$value), ]
  
  # Drought event count
  in_drought <- df$value < drought_threshold
  rle_res <- rle(in_drought)
  n_ev <- sum(rle_res$values & rle_res$lengths >= 1)
  
  # Dynamic y-axis
  y_range <- max(df$value, na.rm = TRUE) - min(df$value, na.rm = TRUE)
  y_pad   <- y_range * 0.10
  y_lo    <- min(-3.5, min(df$value, na.rm = TRUE) - y_pad)
  y_hi    <- max( 3.5, max(df$value, na.rm = TRUE) + y_pad)
  
  # Drought bands (matching 7basin_timeseries.R)
  drought_bands <- list(
    ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -2.0, fill = "#7B0025", alpha = 0.15),
    ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = -2.0, ymax = -1.5, fill = "#D73027", alpha = 0.15),
    ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = -1.5, ymax = -1.0, fill = "#F46D43", alpha = 0.15),
    ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = -1.0, ymax =  1.0, fill = "#2E7D32", alpha = 0.15),
    ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin =  1.0, ymax =  1.5, fill = "#90EE90", alpha = 0.15),
    ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin =  1.5, ymax =  2.0, fill = "#66BB6A", alpha = 0.15),
    ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin =  2.0, ymax =  Inf, fill = "#4575B4", alpha = 0.15)
  )
  
  ggplot2::ggplot(df, ggplot2::aes(x = date, y = value)) +
    drought_bands +
    ggplot2::geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
    ggplot2::geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "gray50", linewidth = 0.3) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = ifelse(value < drought_threshold, value, 0), ymax = 0), 
                         fill = "#c0392b", alpha = 0.35) +
    ggplot2::geom_line(color = color, linewidth = 0.6) +
    ggplot2::scale_x_date(date_breaks = "5 years", date_labels = "%Y") +
    ggplot2::scale_y_continuous(limits = c(y_lo, y_hi)) +
    ggplot2::labs(
      title    = sprintf("%s Basin-Averaged Time Series", toupper(index_name)),
      subtitle = sprintf("%d drought events detected (onset < %.1f)", n_ev, drought_threshold),
      x = "Year", y = sprintf("%s Index Value", toupper(index_name))
    ) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "gray30"),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom"
    )
}
## 3a. Pivot wide month-CSV (Year + 3-letter month cols) → long data frame ----
pivot_wide_monthly <- function(path, value_name) {
  df <- read_csv(path, show_col_types = FALSE)
  result <- df %>%
    pivot_longer(-Year, names_to = "month_abbr", values_to = value_name) %>%
    mutate(date = as.Date(paste0(Year, "-", month_abbr, "-01"), "%Y-%b-%d")) %>%
    filter(!is.na(date), !is.na(.data[[value_name]]))
  result[, c("date", value_name), drop = FALSE]}

## 3b. Reference-period empirical CDF transform (per calendar month) ----------
calibrate_marginals <- function(df, value_col, ref_start, ref_end) {
  ref <- df |>
    filter(date >= ref_start, date <= ref_end) |>
    mutate(mon = month(date))
  
  lapply(1:12, function(m) {
    vals <- ref[[value_col]][ref$mon == m]
    vals <- vals[!is.na(vals)]
    if (length(vals) < 5L) return(NULL)
    list(fn = ecdf(vals), n = length(vals))   # store n explicitly alongside ecdf
  })
}

## 3c. Apply calibrated ECDFs to full series → U[0,1] -------------------------
apply_marginals <- function(df, value_col, ecdfs) {
  u    <- rep(NA_real_, nrow(df))
  mons <- month(df$date)
  for (m in 1:12) {
    idx <- which(mons == m & !is.na(df[[value_col]]))
    if (length(idx) == 0L || is.null(ecdfs[[m]])) next
    n_ref  <- ecdfs[[m]]$n                          # safe: stored explicitly
    raw    <- ecdfs[[m]]$fn(df[[value_col]][idx])   # apply the ecdf
    raw    <- pmax(1 / (n_ref + 1), pmin(n_ref / (n_ref + 1), raw))
    u[idx] <- raw
  }
  u
}

## 3d. Kendall distribution via Monte-Carlo from a fitted vine -----------------
# More rigorous than simple Blom ranking on the observed sample:
# simulate n_sim obs from the vine, evaluate the vine CDF at each, and the
# empirical distribution of those CDF values IS the estimated KC(t).
estimate_kendall_cdf <- function(vine_fit, n_sim, n_mc) {
  sim_u  <- rvinecop(n_sim, vine = vine_fit)        # d-dimensional uniform
  p_sim  <- pvinecop(sim_u, vine = vine_fit, n_mc = n_mc)
  ecdf(p_sim)                                        # return KC as an ecdf object
}

## 3e. Transform observed vine CDF values → MSDI using estimated KC ----------
vine_to_msdi <- function(p_vine, KC_ecdf) {
  kc_vals <- KC_ecdf(p_vine)
  # clamp to open interval before qnorm
  kc_vals <- pmax(1e-6, pmin(1 - 1e-6, kc_vals))
  qnorm(kc_vals)
}

## 3f. Empirical tail-dependence coefficients for one pair --------------------
# λ_L = P(V ≤ q | U ≤ q) as q → 0
# λ_U = P(V > 1-q | U > 1-q) as q → 1
# Estimated at q = 0.10 and 0.05; averaged.
tail_dep <- function(u, v, q = c(0.05, 0.10)) {
  n    <- length(u)
  lL   <- mean(sapply(q, function(qi)
    sum(u <= qi & v <= qi) / sum(u <= qi)))
  lU   <- mean(sapply(q, function(qi)
    sum(u >= 1 - qi & v >= 1 - qi) / sum(u >= 1 - qi)))
  c(lambda_L = lL, lambda_U = lU)
}

## 3g. Drought event extraction ------------------------------------------------
extract_drought_events <- function(date_vec, msdi_vec,
                                   threshold = -1.0,
                                   min_dur   = 3L) {
  df   <- tibble(date = date_vec, msdi = msdi_vec) |>
    arrange(date) |>
    mutate(in_drought = msdi < threshold)
  
  rle_res <- rle(df$in_drought)
  ends    <- cumsum(rle_res$lengths)
  starts  <- ends - rle_res$lengths + 1L
  
  events <- purrr::map_dfr(seq_along(rle_res$values), function(i) {
    if (!rle_res$values[i]) return(NULL)
    idx <- starts[i]:ends[i]
    if (length(idx) < min_dur) return(NULL)
    tibble(
      start_date      = df$date[idx[1]],
      end_date        = df$date[idx[length(idx)]],
      duration_months = length(idx),
      peak_msdi       = min(df$msdi[idx]),
      severity        = sum(abs(df$msdi[idx][df$msdi[idx] < threshold]),
                            na.rm = TRUE),
      drought_class   = case_when(
        min(df$msdi[idx]) < -2.0 ~ "Extreme",
        min(df$msdi[idx]) < -1.5 ~ "Severe",
        TRUE                     ~ "Moderate"
      )
    )
  })
  
  if (nrow(events) == 0L)
    return(tibble(start_date = as.Date(NA), end_date = as.Date(NA),
                  duration_months = integer(0), peak_msdi = numeric(0),
                  severity = numeric(0), drought_class = character(0)))
  events
}
## 3h. Assemble SSI soil from per-month CSVs -----------------------------------
# The upstream script saved one CSV per calendar month, not a wide basin file.
# This function reads all 12, extracts a basin representative value (median
# across grid cells), and returns a long data frame identical in structure to
# pivot_wide_monthly() output.
#
# Assumptions about the per-month CSV structure (adjust col name if needed):
#   - Contains a numeric column named "ssi" (or similar — detected automatically)
#   - Each row is a grid cell or pixel for that calendar month
#   - Files are named  ssi_L1_3_01_month01_Jan.csv … _month12_Dec.csv
assemble_ssi_from_monthly_csvs <- function(folder, prefix, value_name,
                                           years_range = NULL) {
  month_files <- tibble(
    mon   = 1:12,
    abbr  = c("Jan","Feb","Mar","Apr","May","Jun",
              "Jul","Aug","Sep","Oct","Nov","Dec"),
    path  = file.path(folder,
                      sprintf("%s_month%02d_%s.csv",
                              prefix,
                              1:12,
                              c("Jan","Feb","Mar","Apr","May","Jun",
                                "Jul","Aug","Sep","Oct","Nov","Dec")))
  )
  
  purrr::map_dfr(seq_len(nrow(month_files)), function(i) {
    p <- month_files$path[i]
    if (!file.exists(p)) { warning("Missing: ", p); return(NULL) }
    
    raw <- read_csv(p, show_col_types = FALSE)
    
    # Auto-detect the SSI value column (first numeric col that isn't Year/year/pixel)
    num_cols <- names(raw)[sapply(raw, is.numeric)]
    id_cols  <- grep("^(year|Year|pixel|id|lon|lat|x|y)$",
                     num_cols, value = TRUE, ignore.case = TRUE)
    val_col  <- setdiff(num_cols, id_cols)[1]
    if (is.na(val_col)) { warning("No numeric value col in ", p); return(NULL) }
    
    # Detect year column
    yr_col <- grep("^year$", names(raw), value = TRUE, ignore.case = TRUE)[1]
    
    if (!is.na(yr_col)) {
      # Long format: multiple years, one row per (year × pixel)
      raw %>%
        group_by(Year = .data[[yr_col]]) %>%
        summarise(!!value_name := median(.data[[val_col]], na.rm = TRUE),
                  .groups = "drop") %>%
        mutate(date = as.Date(paste0(Year, "-", month_files$mon[i], "-01"))) %>%
        select(date, !!value_name)
    } else {
      # Wide format: one row per pixel, columns = years
      yr_cols <- names(raw)[grepl("^\\d{4}$", names(raw))]
      if (length(yr_cols) == 0) { warning("Can't detect year cols in ", p); return(NULL) }
      raw %>%
        select(all_of(yr_cols)) %>%
        summarise(across(everything(), ~ median(.x, na.rm = TRUE))) %>%
        pivot_longer(everything(), names_to = "Year", values_to = value_name) %>%
        mutate(Year = as.integer(Year),
               date = as.Date(paste0(Year, "-", month_files$mon[i], "-01"))) %>%
        select(date, !!value_name)
    }
  }) %>%
    filter(!is.na(date), !is.na(.data[[value_name]])) %>%
    { if (!is.null(years_range))
      filter(., between(year(date), years_range[1], years_range[2]))
      else . } %>%
    arrange(date)
}
## 3i. drop_na with row-count diagnostics ─────────────────────────────────
# Wraps drop_na() and emits a warning whenever rows are silently removed,
# so date-alignment gaps from the full_join merge are never invisible.
drop_na_logged <- function(df, label) {
  n_in  <- nrow(df)
  out   <- tidyr::drop_na(df)
  n_out <- nrow(out)
  if (n_out < n_in)
    message(sprintf(
      "  \u26A0 %s: %d of %d rows retained after drop_na (%d removed — check date alignment).",
      label, n_out, n_in, n_in - n_out))
  else
    message(sprintf("  \u2714 %s: all %d rows complete (no NAs dropped).", label, n_out))
  out
}

# ==============================================================================
# 4.  LOAD & MERGE INDICES
# ==============================================================================
cat("\n", strrep("=", 60), "\n  LOADING INDEX TIME SERIES\n",
    strrep("=", 60), "\n\n", sep = "")

# ── 4a. Snow index: SWEI or SSPI-1 ────────────────────────────────────────
snow_col <- if (cfg$snow_index == "sspi1") "SSPI1" else "SWEI"

if (cfg$snow_index == "sspi1" && file.exists(cfg$sspi1_file)) {
  snow_long <- pivot_wide_monthly(cfg$sspi1_file, "SSPI1")
  cat(sprintf("  Snow (SSPI-1):  %d records (%s – %s)\n",
              nrow(snow_long), min(snow_long$date), max(snow_long$date)))
} else {
  snow_long <- pivot_wide_monthly(cfg$swei_file, "SWEI")
  cat(sprintf("  Snow (SWEI):    %d records (%s – %s)\n",
              nrow(snow_long), min(snow_long$date), max(snow_long$date)))
  if (cfg$snow_index == "sspi1")
    message("  ⚠ SSPI-1 file not found — falling back to SWEI.")
}

# ── 4b. SSI soil moisture ─────────────────────────────────────────────────
ssi_soil_long <- assemble_ssi_from_monthly_csvs(
  folder     = cfg$ssi_soil_folder,
  prefix     = cfg$ssi_soil_prefix,
  value_name = "SSI_soil"
)
cat(sprintf("  SSI soil:       %d records (%s – %s)\n",                      
            nrow(ssi_soil_long), min(ssi_soil_long$date), max(ssi_soil_long$date)))
# ── 4c. SGI groundwater ───────────────────────────────────────────────────
# 10SGI.R saves: month_label (YYYY-MM), sgi_median, date
sgi_raw  <- read_csv(cfg$sgi_file, show_col_types = FALSE)
sgi_long <- sgi_raw |>                                                       
  mutate(date = if ("date" %in% names(sgi_raw))
    floor_date(as.Date(date), "month")
    else
      as.Date(paste0(month_label, "-01"))) |>
  select(date, SGI = sgi_median) |>
  filter(!is.na(SGI))
cat(sprintf("  SGI:            %d records (%s – %s)\n",
            nrow(sgi_long), min(sgi_long$date), max(sgi_long$date)))

# ── 4d. WSC streamflow SSI ────────────────────────────────────────────────
# REGULATORY NOTE:
#   08JC001  Nechako River at Vanderhoof — DOWNSTREAM of Kenney Dam.
#            Records regulated RELEASES governed by BC Hydro power scheduling.
#   08JE001  Stuart River near Fort St. James — natural tributary (unregulated).
#            Drains Stuart Lake into the Nechako below the reservoir.
# Basin-averaging these two stations produces a MIXED signal: partly natural
# (Stuart) and partly operational (Nechako main stem).  The vine copula will
# capture the net statistical association; causal interpretation of the
# SSI_wsc node must account for this operational confounding.
#
# 7WSC_StreamFlow.R saves date as the 15th of each month; snap to month-start.
ssi_wsc_list <- lapply(cfg$ssi_wsc_rds, function(f) {
  if (!file.exists(f)) { warning("Not found: ", f); return(NULL) }
  readRDS(f) |>
    mutate(date = floor_date(as.Date(date), "month")) |>
    select(date, ssi) |>
    filter(!is.na(ssi))
})
ssi_wsc_list <- Filter(Negate(is.null), ssi_wsc_list)
if (length(ssi_wsc_list) == 0) stop("No WSC SSI RDS found — check cfg$ssi_wsc_rds.")

ssi_wsc_long <- bind_rows(ssi_wsc_list) |>
  group_by(date) |>
  summarise(SSI_wsc = mean(ssi, na.rm = TRUE), .groups = "drop") |>
  filter(!is.na(SSI_wsc))
cat(sprintf("  SSI stream:     %d records (%s – %s)\n",
            nrow(ssi_wsc_long), min(ssi_wsc_long$date), max(ssi_wsc_long$date)))

# ── 4e. NTSDI Nechako Reservoir ───────────────────────────────────────────
# MANAGED RESERVOIR NOTE:
#   The Nechako Reservoir (created by Kenney Dam, 1952; ~890 km², ~850–854 m ASL)
#   is operated by BC Hydro for hydropower.  Water is diverted westward through
#   the Nechako Tunnel to the Kemano generating station, NOT released downstream.
#   Reservoir level therefore reflects:
#     (a) Natural inflow from basin snowmelt / rainfall (coupled to SWEI/SSI/SGI)
#     (b) BC Hydro's diversion schedule (power demand — operationally decoupled)
#     (c) Minimum downstream flow obligations (environmental flow orders)
#   NTSDI is a valid DROUGHT INDICATOR (low storage = stressed system) but its
#   statistical coupling to upstream natural variables is confounded by operations.
#   The R-vine structure is preferred because it makes no cascade assumption.
#
# 9NTSDI writes: date, NTSDI_norm  (column name = paste0("NTSDI","_norm"))
ntsdi_raw  <- read_csv(cfg$ntsdi_file, show_col_types = FALSE)
ntsdi_long <- ntsdi_raw |>                                                       
  mutate(date = floor_date(as.Date(date), "month")) |>
  select(date, NTSDI = NTSDI_norm) |>
  filter(!is.na(NTSDI))
cat(sprintf("  NTSDI:          %d records (%s – %s)\n",
            nrow(ntsdi_long), min(ntsdi_long$date), max(ntsdi_long$date)))

# ── 4f. Merge all five ────────────────────────────────────────────────────
all_raw <- snow_long |>
  rename(Snow = all_of(snow_col)) |>
  full_join(ssi_soil_long, by = "date") |>
  full_join(sgi_long,      by = "date") |>
  full_join(ssi_wsc_long,  by = "date") |>
  full_join(ntsdi_long,    by = "date") |>
  arrange(date)

cat(sprintf("\n  ➤ Combined record: %d months (%s – %s)\n",
            nrow(all_raw), min(all_raw$date), max(all_raw$date)))

# ── Date-alignment diagnostic: flag any column with unexpected NAs ─────────
na_counts <- colSums(is.na(all_raw[, c("Snow", "SSI_soil", "SGI", "SSI_wsc", "NTSDI")]))
if (any(na_counts > 0)) {
  message("  \u26A0 NA counts in merged table (may indicate date-format misalignment):")
  for (nm in names(na_counts)[na_counts > 0])
    message(sprintf("      %-12s  %d NAs  (%.1f%% of %d months)",
                    nm, na_counts[nm], 100 * na_counts[nm] / nrow(all_raw), nrow(all_raw)))
} else {
  message("  \u2714 All five index columns fully populated — no date-alignment gaps detected.")
}
cat("\n")

# ==============================================================================
# 5.  REFERENCE-PERIOD CALIBRATED PSEUDO-OBSERVATIONS
# ==============================================================================
# Per calendar-month ECDFs are fitted on the reference period only, then
# applied to the full series.  This prevents out-of-reference anomalies from
# shifting the marginal rank structure of the vine.

cat(strrep("=", 60), "\n  CALIBRATING MARGINAL ECDF TRANSFORMS\n",
    strrep("=", 60), "\n\n", sep = "")

# Variable order is fixed here and must match dvine_order_* indices in cfg:
#   Snow=1  SSI_soil=2  SGI=3  NTSDI=4  SSI_wsc=5
# NTSDI and SSI_wsc are placed last because they are managed/regulated;
# their mutual ordering is resolved by BIC comparison of D-vine variants A/B.
var_list   <- c("Snow", "SSI_soil", "SGI", "NTSDI", "SSI_wsc")
ecdf_cache <- list()

for (v in var_list) {
  ecdf_cache[[v]] <- calibrate_marginals(all_raw, v, cfg$ref_start, cfg$ref_end)
  u_col <- paste0("u_", v)
  all_raw[[u_col]] <- apply_marginals(all_raw, v, ecdf_cache[[v]])
  n_u <- sum(!is.na(all_raw[[u_col]]))
  cat(sprintf("  %-12s  %d uniform values (%.1f%% coverage)\n",
              v, n_u, 100 * n_u / nrow(all_raw)))
}

# ==============================================================================
# 6.  TAIL-DEPENDENCE DIAGNOSTICS (all 10 pairs)
# ==============================================================================
cat("\n", strrep("=", 60), "\n  TAIL-DEPENDENCE MATRIX (λ_L / λ_U)\n",
    strrep("=", 60), "\n\n", sep = "")

u_cols <- paste0("u_", var_list)

# Compute pairwise Kendall tau and tail coefficients using reference period
ref_u <- all_raw |>
  filter(date >= cfg$ref_start, date <= cfg$ref_end) |>
  select(all_of(u_cols)) |>
  drop_na_logged("Tail-dependence reference matrix (all months)")

td_df <- purrr::map_dfr(
  seq_along(var_list),
  function(i) purrr::map_dfr(seq_along(var_list), function(j) {
    if (j <= i) return(NULL)
    u <- ref_u[[paste0("u_", var_list[i])]]
    v <- ref_u[[paste0("u_", var_list[j])]]
    td <- tail_dep(u, v)
    tibble(
      Var1    = var_list[i], Var2 = var_list[j],
      tau     = cor(u, v, method = "kendall"),
      lambda_L = td["lambda_L"],
      lambda_U = td["lambda_U"]
    )
  })
)

print(td_df |> mutate(across(where(is.numeric), \(x) round(x, 3))),
      n = nrow(td_df))
cat("\n  Note: λ_L > 0.2 → lower-tail (drought co-occurrence) dependence\n")
cat("        λ_U > 0.2 → upper-tail (wet co-occurrence) dependence\n\n")

# Save tail-dep table
write_csv(td_df, file.path(cfg$out_dir, "tail_dependence_matrix.csv"))

# ==============================================================================
# 7.  PAIRWISE KENDALL TAU MATRIX
# ==============================================================================
tau_mat <- cor(ref_u, method = "kendall")
rownames(tau_mat) <- colnames(tau_mat) <- var_list
cat(strrep("=", 60), "\n  PAIRWISE KENDALL TAU (reference period)\n",
    strrep("=", 60), "\n\n", sep = "")
print(round(tau_mat, 3))

# ==============================================================================
# 8.  VINE COPULA FITTING — SEASONAL STRATIFICATION (OPTION A ONLY)
# ==============================================================================
# Two vine regimes are fitted (Option A — exclusive seasonal stratification):
#   (B) Snow-season 5-variable vine (D + R)  — Oct–Apr
#   (C) Warm-season 4-variable vine (D + R)  — May–Sep (Snow excluded)
#
# Full-year 5-variable vine REMOVED — structurally incompatible with seasonal SWEI
#
# D-vine cascade logic:
#   • Positions 1–3 (natural, unambiguous):  Snow → SSI_soil → SGI
#   • Positions 4–5 (managed/regulated):     NTSDI and SSI_wsc
#     The mutual ordering of NTSDI and SSI_wsc is NOT physically determined
#     because the Nechako Reservoir is operated for BC Hydro hydropower.
#     Both tail variants are fitted; the lower-BIC variant is retained.
#     The R-vine is the PRIMARY recommended model because it imposes no
#     cascade structure on the operationally managed variables.

fit_vine_pair <- function(u_mat_fit, dvine_ord_A, dvine_ord_B, label) {
  
  # ── D-vine variant A ──────────────────────────────────────────────────────
  cat(sprintf("\n  Fitting D-vine-A [%s] …  ", label))
  vine_dA <- vinecop(
    u_mat_fit,
    structure  = dvine_structure(dvine_ord_A),
    family_set = cfg$family_set,
    par_method = cfg$par_method,
    selcrit    = cfg$selcrit,
    cores      = max(1L, parallel::detectCores() - 1L)
  )
  cat(sprintf("logLik=%.1f  BIC=%.1f\n", logLik(vine_dA), BIC(vine_dA)))
  
  # ── D-vine variant B ──────────────────────────────────────────────────────
  cat(sprintf("  Fitting D-vine-B [%s] …  ", label))
  vine_dB <- vinecop(
    u_mat_fit,
    structure  = dvine_structure(dvine_ord_B),
    family_set = cfg$family_set,
    par_method = cfg$par_method,
    selcrit    = cfg$selcrit,
    cores      = max(1L, parallel::detectCores() - 1L)
  )
  cat(sprintf("logLik=%.1f  BIC=%.1f\n", logLik(vine_dB), BIC(vine_dB)))
  
  # Select the better D-vine tail ordering by BIC
  vine_d       <- if (BIC(vine_dA) <= BIC(vine_dB)) vine_dA else vine_dB
  dvine_winner <- if (BIC(vine_dA) <= BIC(vine_dB)) "A" else "B"
  cat(sprintf("  D-vine tail ordering selected: Variant %s (ΔBIC vs other = %.2f)\n",
              dvine_winner, abs(BIC(vine_dA) - BIC(vine_dB))))
  
  cat("  NOTE: Variant A tail = [NTSDI → SSI_wsc]; Variant B tail = [SSI_wsc → NTSDI]\n")
  cat("        Neither tail ordering has a causal physical basis (managed reservoir); BIC selection is purely statistical.\n")
  # ── R-vine (data-driven — PRIMARY for managed system) ────────────────────
  cat(sprintf("  Fitting R-vine  [%s] (PRIMARY) …  ", label))
  vine_r <- vinecop(
    u_mat_fit,
    family_set = cfg$family_set,
    par_method = cfg$par_method,
    selcrit    = cfg$selcrit,
    cores      = max(1L, parallel::detectCores() - 1L)
  )
  cat(sprintf("logLik=%.1f  BIC=%.1f\n", logLik(vine_r), BIC(vine_r)))
  
  # ── Vine truncation comparison (applied to R-vine) ────────────────────────
  trunc_bic <- tibble(truncation = "full", BIC = BIC(vine_r))
  for (k in setdiff(cfg$truncation_levels, list(NULL))) {
    vt <- tryCatch(truncate_model(vine_r, trunc_lvl = k), error = function(e) NULL)
    if (!is.null(vt))
      trunc_bic <- bind_rows(trunc_bic,
                             tibble(truncation = paste0("tree", k), BIC = BIC(vt)))
  }
  cat("  R-vine truncation BIC:\n")
  print(trunc_bic)
  
  # ── Best overall (D-vine vs R-vine) ──────────────────────────────────────
  # For reporting: select best by BIC, but note the R-vine recommendation
  which_best <- if (BIC(vine_d) <= BIC(vine_r)) "D-vine" else "R-vine"
  vine_best  <- if (which_best == "D-vine") vine_d else vine_r
  cat(sprintf("  BIC-selected model: %s",  which_best))
  if (which_best == "D-vine")
    cat("  ← note: R-vine recommended for managed variables despite BIC\n")
  else
    cat("  ✓ (consistent with R-vine recommendation for managed system)\n")
  
  list(
    dA         = vine_dA,
    dB         = vine_dB,
    d          = vine_d,         # best-BIC D-vine variant
    dvine_tail = dvine_winner,   # "A" or "B"
    r          = vine_r,         # PRIMARY recommended model
    best       = vine_best,
    which      = which_best,
    trunc      = trunc_bic
  )
}

# ── 8B. Snow season (Oct–Apr, 5 variables) ────────────────────────────────
cat("\n", strrep("=", 60),
    "\n  VINE FITTING — SNOW SEASON Oct–Apr (5 variables)\n",
    strrep("=", 60), "\n", sep = "")

snow5_ref <- all_raw |>
  filter(date >= cfg$ref_start, date <= cfg$ref_end,
         month(date) %in% cfg$snow_months) |>
  select(all_of(u_cols)) |>
  drop_na_logged("Snow-season reference period (5-var)")

if (nrow(snow5_ref) < 40L)
  warning("Fewer than 40 complete snow-season months — estimates unreliable.")

u_snow5 <- as.matrix(snow5_ref)
colnames(u_snow5) <- var_list

vines_snow5 <- fit_vine_pair(u_snow5, cfg$dvine_order_5A, cfg$dvine_order_5B, "Snow-season")
saveRDS(vines_snow5, file.path(cfg$out_dir, "vines_snow5.rds"))

# ── 8C. Warm season (May–Sep, 4 variables — Snow excluded) ───────────────
cat("\n", strrep("=", 60),
    "\n  VINE FITTING — WARM SEASON May–Sep (4 variables)\n",
    strrep("=", 60), "\n", sep = "")

# Drop Snow column; variable order becomes: SSI_soil=1, SGI=2, NTSDI=3, SSI_wsc=4
u_cols_warm <- paste0("u_", c("SSI_soil", "SGI", "NTSDI", "SSI_wsc"))
var_list_warm <- c("SSI_soil", "SGI", "NTSDI", "SSI_wsc")

warm4_ref <- all_raw |>
  filter(date >= cfg$ref_start, date <= cfg$ref_end,
         month(date) %in% cfg$warm_months) |>
  select(all_of(u_cols_warm)) |>
  drop_na_logged("Warm-season reference period (4-var)")

if (nrow(warm4_ref) < 30L)
  warning("Fewer than 30 complete warm-season months — estimates unreliable.")

u_warm4 <- as.matrix(warm4_ref)
colnames(u_warm4) <- var_list_warm

vines_warm4 <- fit_vine_pair(u_warm4, cfg$dvine_order_4A, cfg$dvine_order_4B, "Warm-season")
saveRDS(vines_warm4, file.path(cfg$out_dir, "vines_warm4.rds"))

# ==============================================================================
# ==== SECTION 9: VINE CDF EVALUATION ====
# ==============================================================================
set.seed(40L)   # ← random seeds pvinecop() MC integration
cat("\n", strrep("=", 60), "\n  EVALUATING VINE CDF (full record)\n",
    strrep("=", 60), "\n\n", sep = "")

# Seasonal evaluation only (Option A): snow season (5 vars) + warm season (4 vars)
all_snow <- all_raw |>
  filter(month(date) %in% cfg$snow_months) |>
  select(date, all_of(u_cols)) |>
  drop_na_logged("Snow-season full record (CDF evaluation)")

all_warm <- all_raw |>
  filter(month(date) %in% cfg$warm_months) |>
  select(date, all_of(u_cols_warm)) |>
  drop_na_logged("Warm-season full record (CDF evaluation)")

# for snow season
p_snow_best <- pvinecop(
  `colnames<-`(as.matrix(all_snow[, u_cols]), var_list),
  vines_snow5$best, n_mc = cfg$n_mc
)

# for warm season
p_warm_best <- pvinecop(
  `colnames<-`(as.matrix(all_warm[, u_cols_warm]), var_list_warm),
  vines_warm4$best, n_mc = cfg$n_mc
)

# ==============================================================================
# ==== SECTION 10: KENDALL DISTRIBUTION ESTIMATION VIA MONTE-CARLO ====
# ==============================================================================
set.seed(40L)   # ← random seeds rvinecop() inside estimate_kendall_cdf()
cat("\n  Estimating Kendall distributions (MC, n =",
    format(cfg$n_kc_sim, big.mark = ","), ") …\n")

# Seasonal Kendall distributions only (Option A)
KC_snow5   <- estimate_kendall_cdf(vines_snow5$best, cfg$n_kc_sim, cfg$n_mc)
KC_warm4   <- estimate_kendall_cdf(vines_warm4$best, cfg$n_kc_sim, cfg$n_mc)
cat("  Kendall distributions estimated.\n\n")

# ==============================================================================
# 11. MSDI COMPUTATION — SEASONAL ONLY (OPTION A)
# ==============================================================================
# Seasonal MSDI (best-structure vine per season)
MSDI_snow <- vine_to_msdi(p_snow_best, KC_snow5)
MSDI_warm <- vine_to_msdi(p_warm_best, KC_warm4)

# Combine seasonal MSDI into a single continuous time series (primary output)
seasonal_out <- bind_rows(
  tibble(date = all_snow$date, MSDI_seasonal = MSDI_snow,
         season = "Snow (Oct-Apr)"),
  tibble(date = all_warm$date, MSDI_seasonal = MSDI_warm,
         season = "Warm (May-Sep)")
) |> arrange(date)

# ==============================================================================
# ==== SECTION 12: COMPONENT CONTRIBUTIONS - SEASONAL  ==== 
# ==============================================================================
# For each observation and variable, estimate how much that variable's anomaly
# drives the joint drought.  Method: compute vine CDF with that variable's
# pseudo-obs replaced by the neutral level (0.50) and take the difference
# from the original CDF.  Larger positive difference = larger contribution to
# drought (variable pulled the joint probability down).
set.seed(40L)   # ← random seeds pvinecop() inside compute_contributions()

cat(strrep("=", 60),
    "\n  COMPUTING COMPONENT CONTRIBUTIONS\n",
    strrep("=", 60), "\n\n", sep = "")

compute_contributions <- function(u_full_mat, vine_obj, var_names, n_mc) {
  p_base <- pvinecop(u_full_mat, vine_obj, n_mc = n_mc)
  contrib <- matrix(NA_real_, nrow = nrow(u_full_mat), ncol = length(var_names),
                    dimnames = list(NULL, paste0("contrib_", var_names)))
  
  for (k in seq_along(var_names)) {
    u_neu      <- u_full_mat
    u_neu[, k] <- 0.50       # replace with neutral (50th percentile)
    p_neu      <- pvinecop(u_neu, vine_obj, n_mc = n_mc)
    contrib[, k] <- p_neu - p_base   # >0 if variable k deepens drought
  }
  
  # Normalise to proportion of total deviation from 0.5
  row_sum <- rowSums(abs(contrib))
  contrib_norm <- sweep(contrib, 1, ifelse(row_sum == 0, 1, row_sum), "/")
  as.data.frame(contrib_norm)
}

# Run on seasonal best models separately, then merge chronologically
contrib_snow <- compute_contributions(
  as.matrix(all_snow[, u_cols]), vines_snow5$best, var_list,
  n_mc = min(5000L, cfg$n_mc)
)
contrib_warm <- compute_contributions(
  as.matrix(all_warm[, u_cols_warm]), vines_warm4$best, var_list_warm,
  n_mc = min(5000L, cfg$n_mc)
)

# Merge contributions chronologically with seasonal labels
contrib_df <- bind_rows(
  tibble(date = all_snow$date, contrib_snow, season = "Snow"),
  tibble(date = all_warm$date, contrib_warm, season = "Warm")
) |> arrange(date)

# ==============================================================================
# 13. DROUGHT EVENT IDENTIFICATION & RANKING
# ==============================================================================
cat(strrep("=", 60),
    "\n  DROUGHT EVENT IDENTIFICATION\n",
    strrep("=", 60), "\n\n", sep = "")

# Use the seasonal MSDI as the primary composite
drought_minor    <- extract_drought_events(
  seasonal_out$date, seasonal_out$MSDI_seasonal,
  threshold = -1.0, min_dur = 3L
)
drought_moderate <- extract_drought_events(
  seasonal_out$date, seasonal_out$MSDI_seasonal,
  threshold = -1.5, min_dur = 2L
)
drought_severe   <- extract_drought_events(
  seasonal_out$date, seasonal_out$MSDI_seasonal,
  threshold = -2.0, min_dur = 1L
)

cat(sprintf("  Minor  events (MSDI < −1.0, ≥ 3 months): %d\n",
            nrow(drought_minor)))
cat(sprintf("  Moderate events (MSDI < −1.5, ≥ 2 months): %d\n",
            nrow(drought_moderate)))
cat(sprintf("  Severe events  (MSDI < −2.0): %d\n\n",
            nrow(drought_severe)))

if (nrow(drought_moderate) > 0) {
  cat("  Top moderate/severe events by severity:\n")
  print(drought_moderate |>
          arrange(desc(severity)) |>
          head(10) |>
          mutate(across(where(is.numeric), \(x) round(x, 2))))
}

write_csv(drought_minor,    file.path(cfg$out_dir, "drought_events_minor.csv"))
write_csv(drought_moderate, file.path(cfg$out_dir, "drought_events_moderate.csv"))
write_csv(drought_severe,   file.path(cfg$out_dir, "drought_events_severe.csv"))

# ==============================================================================
# 14. DIAGNOSTICS
# ==============================================================================
cat("\n", strrep("=", 60), "\n  DIAGNOSTICS\n", strrep("=", 60), "\n\n", sep = "")

# ── 14a. KS normality test — seasonal MSDI only ───────────────────────────
for (nm in c("MSDI_snow", "MSDI_warm")) {
  x   <- get(nm)
  ks  <- ks.test(x, "pnorm")
  sw  <- if (length(x) >= 7L) shapiro.test(sample(x, min(5000L, length(x))))
  else list(statistic = NA, p.value = NA)
  cat(sprintf("  %-14s  KS-p=%.4f (%s)  SW-p=%.4f  mean=%.3f  sd=%.3f\n",
              nm, ks$p.value,
              ifelse(ks$p.value > 0.05, "✓ N(0,1)", "✗ non-N"),
              sw$p.value, mean(x, na.rm = TRUE), sd(x, na.rm = TRUE)))
}

# ── 14b. D-vine vs R-vine model comparison table — seasonal only ──────────
model_cmp <- tibble(
  Season     = c("Snow (Oct-Apr)", "Snow (Oct-Apr)",
                 "Warm (May-Sep)", "Warm (May-Sep)"),
  Model      = rep(c("D-vine", "R-vine"), 2),
  d          = c(5L, 5L, 4L, 4L),
  logLik     = c(logLik(vines_snow5$d),  logLik(vines_snow5$r),
                 logLik(vines_warm4$d),  logLik(vines_warm4$r)),
  AIC        = c(AIC(vines_snow5$d),  AIC(vines_snow5$r),
                 AIC(vines_warm4$d),  AIC(vines_warm4$r)),
  BIC        = c(BIC(vines_snow5$d),  BIC(vines_snow5$r),
                 BIC(vines_warm4$d),  BIC(vines_warm4$r))
)
model_cmp <- model_cmp |>
  group_by(Season) |>
  mutate(ΔBIC = BIC - min(BIC)) |>
  ungroup()

cat("\n  Model comparison table:\n")
print(model_cmp |> mutate(across(where(is.numeric), \(x) round(x, 2))),
      n = nrow(model_cmp))
write_csv(model_cmp, file.path(cfg$out_dir, "model_comparison.csv"))

# ── 14c. Rosenblatt GOF (VineCopula package) — seasonal only ─────────────
cat("\n  Rosenblatt GOF (White test, B=100):\n")
# Skip if as_VineCopula not available; seasonal models use rvinecopulib directly
cat("    (skipped — rvinecopulib does not export as_VineCopula; use VineCopula::RVineGOFTest with manual conversion if needed)\n")

# ==============================================================================
# 15. ASSEMBLE MAIN OUTPUT TABLE — SEASONAL MSDI PRIMARY
# ==============================================================================
msdi_out <- seasonal_out |>
  left_join(contrib_df, by = c("date", "season")) |>
  left_join(all_raw |> select(date, all_of(var_list)), by = "date") |>
  mutate(
    # Drought classification for primary (seasonal) MSDI
    drought_cat = case_when(
      is.na(MSDI_seasonal)    ~ NA_character_,
      MSDI_seasonal >   0.00  ~ "No drought",
      MSDI_seasonal >= -1.00  ~ "Minor",
      MSDI_seasonal >= -1.50  ~ "Moderate",
      MSDI_seasonal >= -2.00  ~ "Severe",
      TRUE                    ~ "Extreme"
    )
  )

write_csv(msdi_out, file.path(cfg$out_dir, "Nechako_MSDI_VineCopula.csv"))
cat("\n  Main CSV written: Nechako_MSDI_VineCopula.csv\n")                 

# ==============================================================================
# 16. OPTIONAL BOOTSTRAP CONFIDENCE BANDS — SEASONAL (OPTIONAL)
# ==============================================================================
if (cfg$run_bootstrap) {
  cat("\n", strrep("=", 60), "\n  BOOTSTRAP (n = ", cfg$n_boot, ")\n",
      strrep("=", 60), "\n\n", sep = "")
  
  set.seed(40L)
  # Bootstrap on seasonal data separately, then combine
  boot_msdi_snow <- matrix(NA_real_, nrow = nrow(all_snow), ncol = cfg$n_boot)
  boot_msdi_warm <- matrix(NA_real_, nrow = nrow(all_warm), ncol = cfg$n_boot)
  
  # Snow season bootstrap
  if (nrow(all_snow) > 10) {
    pb <- txtProgressBar(min = 0, max = cfg$n_boot, style = 3)
    for (b in seq_len(cfg$n_boot)) {
      idx_b <- sample(nrow(all_snow), nrow(all_snow), replace = TRUE)
      u_b   <- as.matrix(all_snow[, u_cols])[idx_b, , drop = FALSE]
      tryCatch({
        vb   <- vinecop(u_b, family_set = cfg$family_set,
                        par_method = cfg$par_method,
                        selcrit    = cfg$selcrit, cores = 1L)
        KC_b <- estimate_kendall_cdf(vb, n_sim = 10000L, n_mc = 2000L)
        p_b  <- pvinecop(as.matrix(all_snow[, u_cols]), vb, n_mc = 2000L)
        boot_msdi_snow[, b] <- vine_to_msdi(p_b, KC_b)
      }, error = function(e) NULL)
      setTxtProgressBar(pb, b)
    }
    close(pb)
  }
  
  # Warm season bootstrap
  if (nrow(all_warm) > 10) {
    pb <- txtProgressBar(min = 0, max = cfg$n_boot, style = 3)
    for (b in seq_len(cfg$n_boot)) {
      idx_b <- sample(nrow(all_warm), nrow(all_warm), replace = TRUE)
      u_b   <- as.matrix(all_warm[, u_cols_warm])[idx_b, , drop = FALSE]
      tryCatch({
        vb   <- vinecop(u_b, family_set = cfg$family_set,
                        par_method = cfg$par_method,
                        selcrit    = cfg$selcrit, cores = 1L)
        KC_b <- estimate_kendall_cdf(vb, n_sim = 10000L, n_mc = 2000L)
        p_b  <- pvinecop(as.matrix(all_warm[, u_cols_warm]), vb, n_mc = 2000L)
        boot_msdi_warm[, b] <- vine_to_msdi(p_b, KC_b)
      }, error = function(e) NULL)
      setTxtProgressBar(pb, b)
    }
    close(pb)
  }
  
  # Combine bootstrap CIs chronologically
  msdi_boot_ci <- bind_rows(
    tibble(date = all_snow$date,
           MSDI_lo = apply(boot_msdi_snow, 1, quantile, 0.05, na.rm = TRUE),
           MSDI_hi = apply(boot_msdi_snow, 1, quantile, 0.95, na.rm = TRUE),
           season = "Snow"),
    tibble(date = all_warm$date,
           MSDI_lo = apply(boot_msdi_warm, 1, quantile, 0.05, na.rm = TRUE),
           MSDI_hi = apply(boot_msdi_warm, 1, quantile, 0.95, na.rm = TRUE),
           season = "Warm")
  ) |> arrange(date)
  
  write_csv(msdi_boot_ci, file.path(cfg$out_dir, "MSDI_bootstrap_CI.csv"))
  cat("  Bootstrap CI saved.\n")
} else {
  msdi_boot_ci <- NULL
}

# ==============================================================================
# 17. PLOTS
# ==============================================================================
cat("\n", strrep("=", 60), "\n  GENERATING PLOTS\n",
    strrep("=", 60), "\n\n", sep = "")

drought_pal <- c(
  "No drought" = "#4575b4",
  "Minor"      = "#fee090",
  "Moderate"   = "#fdae61",
  "Severe"     = "#d73027",
  "Extreme"    = "#a50026"
)

# ── Plot 1: Seasonal MSDI time series (primary output) ───────────────────
df_msdi <- seasonal_out %>% dplyr::select(date, value = MSDI_seasonal, season)

# Plot with seasonal facets and save
p1 <- plot_standard_basin_ts(df_msdi, index_name = "MSDI", color = "#d73027") +
  ggplot2::facet_wrap(~season, ncol = 1, scales = "free_x")

ggsave(file.path(cfg$out_dir, "plot1_MSDI_seasonal.png"),
       p1, width = 13, height = 7, dpi = 300)
cat("  Plot 1: MSDI seasonal time series\n")

# ── Plot 2 REMOVED (Option A: no full-year comparison) ───────────────────
# Plot 2 (D-vine vs R-vine full-year) deleted per Option A

# ── Plot 3: MSDI heat map (Year × Month) ─────────────────────────────────
hm_data <- seasonal_out |>
  mutate(Year = year(date), Month = month(date, label = TRUE, abbr = TRUE)) |>
  filter(!is.na(MSDI_seasonal))

p3 <- ggplot(hm_data, aes(Month, factor(Year, levels = rev(sort(unique(Year)))),
                          fill = MSDI_seasonal)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_gradientn(
    colours = c("#a50026", "#d73027", "#fdae61", "grey90", "#74add1", "#313695"),
    limits  = c(-3, 3), oob = scales::squish,
    name    = "MSDI"
  ) +
  labs(title   = "Nechako Basin — MSDI Heat Map",
       x = NULL, y = NULL) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.y = element_text(size = 7))

ggsave(file.path(cfg$out_dir, "plot3_MSDI_heatmap.png"),
       p3, width = 9, height = max(4, nrow(distinct(hm_data, Year)) * 0.22 + 1.5),
       dpi = 300)
cat("  Plot 3: MSDI heat map\n")

# ── Plot 4: Component contribution stacked area chart ─────────────────────
# Use seasonal contributions merged chronologically
if (any(grepl("^contrib_", names(contrib_df)))) {
  contrib_long <- contrib_df |>
    select(date, season, starts_with("contrib_")) |>
    pivot_longer(-c(date, season),
                 names_to  = "variable",
                 values_to = "contribution") |>
    mutate(
      variable  = gsub("contrib_", "", variable),
      direction = ifelse(contribution >= 0, "drought", "wet"),
      contrib_abs = abs(contribution)
    )
  
  p4 <- ggplot(
    contrib_long |> filter(direction == "drought"),
    aes(date, contrib_abs, fill = variable)
  ) +
    geom_area(position = "stack", alpha = 0.8) +
    scale_fill_brewer(palette = "Set1", name = "Component") +
    labs(
      title    = "Component Drought Contribution",
      subtitle = "Proportion of joint drought driven by each water-cycle variable",
      x = NULL, y = "Normalised contribution"
    ) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())
  
  ggsave(file.path(cfg$out_dir, "plot4_component_contribution.png"),
         p4, width = 13, height = 4.5, dpi = 300)
  cat("  Plot 4: Component contribution chart\n")
}

# ── Plot 5: Pair scatter (pseudo-obs, coloured by seasonal MSDI) ──────────
# Uses the snow-season observations only (all 5 vars available)

# Get the dates corresponding to u_snow5 (reference period, snow months, complete cases)
# FIX: Include 'date' in select() before pull()
snow5_dates <- all_raw |>
  filter(date >= cfg$ref_start, date <= cfg$ref_end,
         month(date) %in% cfg$snow_months) |>
  select(date, all_of(u_cols)) |>
  drop_na_logged("Plot 5 snow-season date index") |>
  pull(date)

# Subset MSDI_snow to match the same dates
MSDI_snow_plot <- seasonal_out |>
  filter(date %in% snow5_dates) |>
  pull(MSDI_seasonal)

pd5 <- as_tibble(u_snow5) |>
  setNames(var_list) |>
  mutate(MSDI = MSDI_snow_plot,
         Drought = cut(MSDI, c(-Inf, -2, -1.5, -1, 0, Inf),
                       labels = c("Extreme", "Severe", "Moderate", "Minor", "Normal")))

if (requireNamespace("GGally", quietly = TRUE)) {
  library(GGally)
  p5  <- ggpairs(
    pd5, columns = 1:5,
    aes(colour = Drought, alpha = 0.5),
    upper = list(continuous = wrap("cor", method = "kendall", size = 3.2)),
    lower = list(continuous = wrap("points", size = 0.7)),
    diag  = list(continuous = wrap("densityDiag", alpha = 0.5)),
    legend = 1
  ) +
    scale_colour_manual(
      name   = "Drought class",
      values = c(Extreme = "#a50026", Severe = "#d73027",
                 Moderate = "#fdae61", Minor = "#fee090", Normal = "#4575b4")
    ) +
    labs(title = "Snow-season Pseudo-obs Pair Plot (Kendall τ in upper panel)") +
    theme_bw(base_size = 9)
  ggsave(file.path(cfg$out_dir, "plot5_pair_plot_snow.png"),
         p5, width = 10, height = 9, dpi = 250)
  cat("  Plot 5: Pair scatter matrix (snow season)\n")
} else {
  cat("  Plot 5 skipped (GGally not installed)\n")
}
# ── Plot 6: Drought event severity vs duration scatter ───────────────────
if (nrow(drought_minor) > 0) {
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    # Fallback build if ggrepel is not installed
    p6 <- ggplot(drought_minor,
                 aes(duration_months, severity, colour = drought_class,
                     size = abs(peak_msdi), label = format(start_date, "%Y-%m"))) +
      geom_point(alpha = 0.8) +
      geom_text(size = 3, vjust = -0.8, show.legend = FALSE) +
      scale_colour_manual(values = c(Moderate = "#fdae61", Severe = "#d73027",
                                     Extreme  = "#a50026"),
                          name = "Peak class") +
      scale_size_continuous(range = c(2, 8), name = "|Peak MSDI|") +
      labs(title    = "Nechako Drought Event Catalogue (MSDI < −1.0)",
           x = "Duration (months)", y = "Severity (deficit sum)") +
      theme_bw(base_size = 11) +
      theme(plot.title = element_text(face = "bold"))
  } else {
    # Standard build with ggrepel
    p6 <- ggplot(drought_minor,
                 aes(duration_months, severity, colour = drought_class,
                     size = abs(peak_msdi), label = format(start_date, "%Y-%m"))) +
      geom_point(alpha = 0.8) +
      ggrepel::geom_text_repel(size = 3, max.overlaps = 15, show.legend = FALSE) +
      scale_colour_manual(values = c(Moderate = "#fdae61", Severe = "#d73027",
                                     Extreme  = "#a50026"),
                          name = "Peak class") +
      scale_size_continuous(range = c(2, 8), name = "|Peak MSDI|") +
      labs(title    = "Nechako Drought Event Catalogue (MSDI < −1.0)",
           x = "Duration (months)", y = "Severity (deficit sum)") +
      theme_bw(base_size = 11) +
      theme(plot.title = element_text(face = "bold"))
  }
  
  ggsave(file.path(cfg$out_dir, "plot6_drought_event_scatter.png"),
         p6, width = 8, height = 6, dpi = 250)
  cat("  Plot 6: Drought event severity vs duration\n")
}

# ── Plot 7: Tail-dependence bar chart ─────────────────────────────────────
p7  <- td_df |>
  mutate(pair = paste0(Var1, "\n", Var2)) |>
  pivot_longer(c(lambda_L, lambda_U), names_to = "type", values_to = "lambda") |>
  mutate(type = recode(type, lambda_L = "Lower tail (λ_L)",
                       lambda_U = "Upper tail (λ_U)")) |>
  ggplot(aes(pair, lambda, fill = type)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_hline(yintercept = 0.2, linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values = c("Lower tail (λ_L)" = "#d73027",
                               "Upper tail (λ_U)" = "#4575b4"),
                    name = NULL) +
  labs(title   = "Empirical Tail-Dependence Coefficients",
       subtitle = "Dashed line: λ = 0.2 (weak dependence threshold)",
       x = NULL, y = "Tail-dep. coefficient λ") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(size = 7),
        plot.title = element_text(face = "bold"),
        legend.position = "top")

ggsave(file.path(cfg$out_dir, "plot7_tail_dependence.png"),
       p7, width = 12, height = 4, dpi = 300)
cat("  Plot 7: Tail-dependence bar chart\n")

# ==============================================================================
# 18. FINAL SUMMARY — SEASONAL ONLY (OPTION A)
# ==============================================================================
cat("\n", strrep("=", 60),
    "\n  VINE COPULA MSDI — COMPLETE (Option A: Seasonal Stratification Only)\n",
    strrep("=", 60), "\n\n", sep = "")

cat(sprintf("  Records (seasonal MSDI):  %d months (%s – %s)\n",
            nrow(seasonal_out), min(seasonal_out$date), max(seasonal_out$date)))

best_snow_label <- vines_snow5$which
best_warm_label <- vines_warm4$which

cat(sprintf("  Snow-season best model:    %s (BIC=%.1f)\n",
            best_snow_label,
            min(BIC(vines_snow5$d), BIC(vines_snow5$r))))
cat(sprintf("  Warm-season best model:    %s (BIC=%.1f)\n",
            best_warm_label,
            min(BIC(vines_warm4$d), BIC(vines_warm4$r))))
cat(sprintf("  Drought events catalogued: %d moderate/severe\n",
            nrow(drought_moderate)))
cat(sprintf("  Outputs in:                %s/\n\n", cfg$out_dir))

cat("  Output files:\n")
for (f in list.files(cfg$out_dir, full.names = FALSE))
  cat("    •", f, "\n")

# ==============================================================================
# END
# ==============================================================================