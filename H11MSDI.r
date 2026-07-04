# ==============================================================================
# NECHAKO BASIN — VINE COPULA MULTIVARIATE STANDARDISED DROUGHT INDEX (MSDI)
# v2 — REVISED per basin-critique discussion:
#   (1) D-vine spine reordered to match empirical Kendall-tau coupling
#       strength instead of an assumed slow storage-cascade (mountainous,
#       snow-dominated basin: fast snowmelt->streamflow response is
#       structurally weighted, not buried in a conditional/high-tree
#       position). A fully data-driven (greedy tau-chaining) order is also
#       fitted and BIC-compared automatically each run (D-vine Variant C).
#   (2) SSI_wsc split into TWO physically distinct nodes instead of one
#       basin-averaged column, using H7WSC_StreamFlow.R's own dual-pipeline
#       output:
#         SSI_wsc_natural   = Pipeline B (naturalized/de-regulated 08JC001
#                              reconstruction + observed 08JE001 tributary)
#         SSI_wsc_regulated = Pipeline C (observed 08JC001, POST-REGULATION
#                              YEARS ONLY — the true BC Hydro operational
#                              signal)
#       This replaces the previous single `mean(08JC001, 08JE001)` column,
#       which silently re-averaged the two pipelines' outputs back together
#       and discarded exactly the natural/regulated distinction H7 computed.
#
# COMPONENT INDICES (6, was 5)
# ─────────────────────────────────────────────────────────────────────────────
#  Variable            Index      Script source                 Hydrological role
#  ──────────────────  ─────────  ────────────────────────────  ───────────────────────
#  SWEI                Snow       8SWEI_Snow.R                  SWE / snowpack state
#  SSI_soil            Soil       3SSI_SSMI_ERALand.R            Soil-moisture deficit
#  SGI                 Ground     10SGI.R                        Groundwater level
#  SSI_wsc_natural     Stream(N)  7WSC_StreamFlow.R Pipeline B   Naturalized/unregulated streamflow
#  SSI_wsc_regulated   Stream(R)  7WSC_StreamFlow.R Pipeline C   Regulated (post-dam) streamflow
#  NTSDI               Reservoir  9NTSDI_GeoLakes_Extended.r     Nechako Reservoir storage
#
# PHYSICAL CASCADE — REVISED
# ─────────────────────────────────────────────────────────────────────────────
#  NATURAL sub-cascade (mountainous, snow-dominated basin — fast response):
#    SWEI -> SSI_wsc_natural -> SSI_soil -> SGI
#    Snowmelt reaches a naturalized streamflow signal quickly (direct
#    runoff / fast subsurface flow typical of steep terrain and thin
#    soils); slower soil-moisture recharge and deep groundwater storage
#    are placed AFTER streamflow, not before it, because the basin's own
#    empirical Kendall-tau matrix (Section 6/7) shows Snow-SSI_soil and
#    SSI_soil-SGI are near-zero while Snow-streamflow and
#    SSI_soil-streamflow are the strongest pairs in the dataset. Putting
#    weakly-correlated pairs adjacent (Tree 1) and strongly-correlated
#    pairs in conditional, higher-tree positions was the flaw in the
#    original cascade order.
#
#  MANAGED tail (no imposed physical order — BIC + R-vine govern):
#    SSI_wsc_regulated — observed post-regulation 08JC001 release signal,
#                        driven by BC Hydro power-demand scheduling.
#    NTSDI             — Nechako Reservoir storage, decoupled from natural
#                        inflow by Kemano diversion / spill rules.
#    Two tail orderings are still fitted and BIC-selects between them
#    (Variant A/B), exactly as before — only the tail members changed
#    (SSI_wsc_regulated instead of the old single averaged SSI_wsc).
#
#  A THIRD D-vine variant (C) is now fitted automatically: a fully
#  data-driven order obtained by greedily chaining variables from the
#  empirical |tau| matrix (see order_dvine_by_tau()). This is compared by
#  BIC alongside variants A and B, so the "assumed cascade" is no longer
#  taken on faith — it is checked against the data every run.
#
#  ⚠  VINE STRUCTURE RECOMMENDATION (unchanged)
#     R-vine: PRIMARY model — data-driven MST structure imposes no physical
#             cascade assumption and is robust to operational management.
#     D-vine: reported for interpretability / sensitivity comparison only.
#
#  ⚠ SEASONAL STRATIFICATION (OPTION A — EXCLUSIVE, unchanged)
#     Snow season (Oct–Apr):  all 6 variables — fit vines_snow
#     Warm season (May–Sep):  5 variables — SWEI excluded — fit vines_warm
#     Full-year vine REMOVED — structurally incompatible with seasonal SWEI
#
# METHOD ENHANCEMENTS OVER DRAFT (unchanged from v1, retained)
# ─────────────────────────────────────────────────────────────────────────────
#  • Reference-period-calibrated empirical marginal CDFs (per calendar month)
#  • Kendall distribution estimated via MC simulation from the fitted vine
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
  ssi_soil_folder = "ssi_results_seasonal",
  ssi_soil_prefix = "ssi_L1_3_01",
  
  # SGI  — long CSV:  columns 'month_label' (YYYY-MM), 'sgi_median', 'date'
  sgi_file =  path.expand("./sgi_results/sgi_nechako_multisource_basin.csv"),
  
  # ── WSC streamflow SSI — SPLIT into natural vs regulated (v2) ──────────────
  # Uses H7WSC_StreamFlow.R's own dual-pipeline output instead of averaging
  # both stations into one column:
  #   Pipeline B ("ssi/monthly_data") — 08JC001 NATURALIZED (de-regulated)
  #     reconstruction + 08JE001 observed (always-unregulated tributary).
  #     Basin-averaged here into a single "natural-equivalent" signal.
  #   Pipeline C ("ssi_C/monthly_data") — 08JC001 ORIGINAL/observed record,
  #     restricted to POST-REGULATION years only. The true regulated signal.
  # NOTE: 08JE001 is intentionally NOT included in the regulated set — it is
  # the same unregulated series in both pipelines and would dilute the
  # regulated signal with unregulated data if included there too.
  ssi_wsc_natural_rds = c(
    "streamflow_results/ssi/monthly_data/08JC001_monthly_ssi.rds",   # Pipeline B: naturalized 08JC001
    "streamflow_results/ssi/monthly_data/08JE001_monthly_ssi.rds"    # Pipeline B: observed 08JE001 (unregulated)
  ),
  ssi_wsc_regulated_rds = "streamflow_results/ssi_C/monthly_data/08JC001_monthly_ssi.rds",  # Pipeline C: observed, post-regulation only
  
  # Optional diagnostic: H7's Pettitt-test homogeneity flag for the
  # naturalized 08JC001 splice (observed pre-reg + naturalized post-reg,
  # merged with no marker of the transition point). A significant result
  # means SSI_wsc_natural should not be treated as a clean single-regime
  # series without caveat.
  wsc_homogeneity_rds = "streamflow_results/ssi/naturalized_homogeneity_check.rds",
  
  # NTSDI — long CSV from 9NTSDI_GeoLakes_Extended.r Pipeline 2
  # Columns: date, NTSDI_norm  (Blom-standardised reservoir index)
  ntsdi_file    =  "ntsdi_results/Nechako_NTSDI_Pipeline2_Output.csv",
  
  # ── Reference period (calibrate pseudo-obs marginals here only) ────────────
  ref_start     = as.Date("1981-01-01"),
  ref_end       = as.Date("2020-12-31"),
  
  # ── Seasonal stratification ─────────────────────────────────────────────────
  snow_months   = c(10, 11, 12, 1, 2, 3, 4),   # Oct–Apr
  warm_months   = c(5, 6, 7, 8, 9),             # May–Sep
  
  # ── Snow index choice ──────────────────────────────────────────────────────
  snow_index    = "swei",
  sspi1_file    = "sspi_results_monthly/sspi_01_basin_averaged_by_month.csv",
  
  # ── Vine copula settings ────────────────────────────────────────────────────
  family_set    = "parametric",   # "tll" for nonparametric (slow, flexible)
  selcrit       = "bic",          # "aic" for larger models
  par_method    = "mle",
  
  # D-VINE ORDER RATIONALE (REVISED — v2)
  # ─────────────────────────────────────────────────────────────────────────
  # Positions reflect empirical coupling strength (Kendall tau), not an
  # assumed slow storage cascade. Natural streamflow is placed adjacent to
  # Snow/soil (fast mountain-basin response); the managed pair
  # (SSI_wsc_regulated, NTSDI) is placed at the tail with BOTH orderings
  # fitted and the lower-BIC one retained, exactly as in v1.
  #
  # A THIRD, fully data-driven order (Variant C) is computed automatically
  # inside fit_vine_pair() via order_dvine_by_tau() and BIC-compared against
  # A and B every run — so this hand-specified rationale is checked against
  # the data rather than assumed.
  #
  # var_list column index mapping — SNOW SEASON (6 vars):
  #   Snow=1  SSI_soil=2  SGI=3  SSI_wsc_natural=4  SSI_wsc_regulated=5  NTSDI=6
  dvine_order_snow_A = c(1L, 4L, 2L, 3L, 5L, 6L),   # Snow->natural flow->soil->GW->regulated->NTSDI
  dvine_order_snow_B = c(1L, 4L, 2L, 3L, 6L, 5L),   # Snow->natural flow->soil->GW->NTSDI->regulated
  
  # var_list_warm column index mapping — WARM SEASON (5 vars, Snow excluded):
  #   SSI_soil=1  SGI=2  SSI_wsc_natural=3  SSI_wsc_regulated=4  NTSDI=5
  dvine_order_warm_A = c(1L, 3L, 2L, 4L, 5L),       # SSI_soil->natural flow->SGI->regulated->NTSDI
  dvine_order_warm_B = c(1L, 3L, 2L, 5L, 4L),       # SSI_soil->natural flow->SGI->NTSDI->regulated
  
  # Vine truncation levels to compare (NULL = fit full vine only)
  truncation_levels = c(1L, 2L),
  
  # Monte-Carlo samples for vine CDF integration (pvinecop) and Kendall dist.
  n_mc          = 2e4L,           # increase to 1e5 for publication
  n_kc_sim      = 1e5L,           # MC draws for Kendall distribution
  
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
  
  in_drought <- df$value < drought_threshold
  rle_res <- rle(in_drought)
  n_ev <- sum(rle_res$values & rle_res$lengths >= 1)
  
  y_range <- max(df$value, na.rm = TRUE) - min(df$value, na.rm = TRUE)
  y_pad   <- y_range * 0.10
  y_lo    <- min(-3.5, min(df$value, na.rm = TRUE) - y_pad)
  y_hi    <- max( 3.5, max(df$value, na.rm = TRUE) + y_pad)
  
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
  result[, c("date", value_name), drop = FALSE]
}

## 3b. Reference-period empirical CDF transform (per calendar month) ----------
calibrate_marginals <- function(df, value_col, ref_start, ref_end) {
  ref <- df |>
    filter(date >= ref_start, date <= ref_end) |>
    mutate(mon = month(date))
  
  lapply(1:12, function(m) {
    vals <- ref[[value_col]][ref$mon == m]
    vals <- vals[!is.na(vals)]
    if (length(vals) < 5L) return(NULL)
    list(fn = ecdf(vals), n = length(vals))
  })
}

## 3c. Apply calibrated ECDFs to full series → U[0,1] -------------------------
apply_marginals <- function(df, value_col, ecdfs) {
  u    <- rep(NA_real_, nrow(df))
  mons <- month(df$date)
  for (m in 1:12) {
    idx <- which(mons == m & !is.na(df[[value_col]]))
    if (length(idx) == 0L || is.null(ecdfs[[m]])) next
    n_ref  <- ecdfs[[m]]$n
    raw    <- ecdfs[[m]]$fn(df[[value_col]][idx])
    raw    <- pmax(1 / (n_ref + 1), pmin(n_ref / (n_ref + 1), raw))
    u[idx] <- raw
  }
  u
}

## 3d. Kendall distribution via Monte-Carlo from a fitted vine -----------------
estimate_kendall_cdf <- function(vine_fit, n_sim, n_mc) {
  sim_u  <- rvinecop(n_sim, vine = vine_fit)
  p_sim  <- pvinecop(sim_u, vine = vine_fit, n_mc = n_mc)
  ecdf(p_sim)
}

## 3e. Transform observed vine CDF values → MSDI using estimated KC ----------
vine_to_msdi <- function(p_vine, KC_ecdf) {
  kc_vals <- KC_ecdf(p_vine)
  kc_vals <- pmax(1e-6, pmin(1 - 1e-6, kc_vals))
  qnorm(kc_vals)
}

## 3f. Empirical tail-dependence coefficients for one pair --------------------
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
    
    num_cols <- names(raw)[sapply(raw, is.numeric)]
    id_cols  <- grep("^(year|Year|pixel|id|lon|lat|x|y)$",
                     num_cols, value = TRUE, ignore.case = TRUE)
    val_col  <- setdiff(num_cols, id_cols)[1]
    if (is.na(val_col)) { warning("No numeric value col in ", p); return(NULL) }
    
    yr_col <- grep("^year$", names(raw), value = TRUE, ignore.case = TRUE)[1]
    
    if (!is.na(yr_col)) {
      raw %>%
        group_by(Year = .data[[yr_col]]) %>%
        summarise(!!value_name := median(.data[[val_col]], na.rm = TRUE),
                  .groups = "drop") %>%
        mutate(date = as.Date(paste0(Year, "-", month_files$mon[i], "-01"))) %>%
        dplyr::select(date, !!value_name)
    } else {
      yr_cols <- names(raw)[grepl("^\\d{4}$", names(raw))]
      if (length(yr_cols) == 0) { warning("Can't detect year cols in ", p); return(NULL) }
      raw %>%
        dplyr::select(all_of(yr_cols)) %>%
        summarise(across(everything(), ~ median(.x, na.rm = TRUE))) %>%
        pivot_longer(everything(), names_to = "Year", values_to = value_name) %>%
        mutate(Year = as.integer(Year),
               date = as.Date(paste0(Year, "-", month_files$mon[i], "-01"))) %>%
        dplyr::select(date, !!value_name)
    }
  }) %>%
    filter(!is.na(date), !is.na(.data[[value_name]])) %>%
    { if (!is.null(years_range))
      filter(., between(year(date), years_range[1], years_range[2]))
      else . } %>%
    arrange(date)
}

## 3i. drop_na with row-count diagnostics ─────────────────────────────────
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

## 3j. Data-driven D-vine spine order via greedy Kendall-tau chaining --------
# Builds a variable ordering for the D-vine spine directly from the
# empirical Kendall tau matrix, rather than assuming a physical cascade.
# Starts from the most strongly correlated pair (max |tau|), then greedily
# extends the chain by attaching whichever remaining variable has the
# strongest |tau| to either current end of the chain. This exists because,
# in this mountainous snow-dominated basin, an assumed physical cascade can
# place weakly-correlated variables adjacent in the vine (Tree 1, fully
# flexible) while strongly-correlated pairs are pushed into conditional,
# higher-tree positions — exactly backwards. This order is fitted as
# D-vine Variant C and BIC-compared against the hand-specified variants.
order_dvine_by_tau <- function(u_mat) {
  var_names <- colnames(u_mat)
  d <- ncol(u_mat)
  if (d < 2L) return(seq_len(d))
  
  tau <- suppressWarnings(cor(u_mat, method = "kendall", use = "pairwise.complete.obs"))
  diag(tau) <- 0
  
  ut <- upper.tri(tau)
  best_flat <- which(abs(tau) == max(abs(tau[ut]), na.rm = TRUE) & ut, arr.ind = TRUE)[1, ]
  chain <- var_names[best_flat]
  remaining <- setdiff(var_names, chain)
  
  while (length(remaining) > 0) {
    left_end  <- chain[1]
    right_end <- chain[length(chain)]
    scores <- sapply(remaining, function(v)
      max(abs(tau[v, left_end]), abs(tau[v, right_end]), na.rm = TRUE))
    nxt <- remaining[which.max(scores)]
    if (abs(tau[nxt, left_end]) >= abs(tau[nxt, right_end])) {
      chain <- c(nxt, chain)
    } else {
      chain <- c(chain, nxt)
    }
    remaining <- setdiff(remaining, nxt)
  }
  match(chain, var_names)
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
sgi_raw  <- read_csv(cfg$sgi_file, show_col_types = FALSE)
sgi_long <- sgi_raw |>
  mutate(date = if ("date" %in% names(sgi_raw))
    floor_date(as.Date(date), "month")
    else
      as.Date(paste0(month_label, "-01"))) |>
  dplyr::select(date, SGI = sgi_median) |>
  filter(!is.na(SGI))
cat(sprintf("  SGI:            %d records (%s – %s)\n",
            nrow(sgi_long), min(sgi_long$date), max(sgi_long$date)))

# ── 4d. WSC streamflow SSI — NATURAL vs REGULATED (v2 split) ──────────────
# REGULATORY NOTE:
#   08JC001  Nechako River at Vanderhoof — DOWNSTREAM of Kenney Dam.
#   08JE001  Stuart River near Fort St. James — natural tributary (unregulated).
#
# H7WSC_StreamFlow.R produces TWO relevant pipelines for 08JC001:
#   Pipeline B ("ssi/monthly_data")   — NATURALIZED (de-regulated) flow:
#       pre-regulation observed + post-regulation naturalized reconstruction,
#       pre-merged upstream, Pettitt-tested for an undocumented splice.
#   Pipeline C ("ssi_C/monthly_data") — ORIGINAL/observed flow, restricted to
#       POST-REGULATION years only — the real, as-measured regulated signal.
# 08JE001 (Stuart River) is unregulated and identical in both pipelines, so
# it is only used once, folded into the "natural" node alongside the
# naturalized 08JC001 reconstruction.
#
# This produces two SEPARATE nodes instead of one basin-averaged SSI_wsc,
# preserving the natural/regulated distinction the previous v1 script
# discarded at the merge step.

## -- Natural streamflow signal (Pipeline B: naturalized 08JC001 + 08JE001) --
ssi_wsc_natural_list <- lapply(cfg$ssi_wsc_natural_rds, function(f) {
  if (!file.exists(f)) { warning("Not found: ", f); return(NULL) }
  readRDS(f) |>
    mutate(date = floor_date(as.Date(date), "month")) |>
    dplyr::select(date, ssi) |>
    filter(!is.na(ssi))
})
ssi_wsc_natural_list <- Filter(Negate(is.null), ssi_wsc_natural_list)
if (length(ssi_wsc_natural_list) == 0) stop("No natural WSC SSI RDS found — check cfg$ssi_wsc_natural_rds.")

ssi_wsc_natural_long <- bind_rows(ssi_wsc_natural_list) |>
  group_by(date) |>
  summarise(SSI_wsc_natural = mean(ssi, na.rm = TRUE), .groups = "drop") |>
  filter(!is.na(SSI_wsc_natural))

cat(sprintf("  SSI stream (natural):    %d records (%s – %s)\n",
            nrow(ssi_wsc_natural_long), min(ssi_wsc_natural_long$date), max(ssi_wsc_natural_long$date)))

## -- Regulated streamflow signal (Pipeline C: observed 08JC001, post-reg) --
if (!file.exists(cfg$ssi_wsc_regulated_rds))
  stop("Regulated WSC SSI RDS not found: ", cfg$ssi_wsc_regulated_rds,
       "  (check cfg$ssi_wsc_regulated_rds against H7's ssi_C/monthly_data output).")

ssi_wsc_regulated_long <- readRDS(cfg$ssi_wsc_regulated_rds) |>
  mutate(date = floor_date(as.Date(date), "month")) |>
  dplyr::select(date, SSI_wsc_regulated = ssi) |>
  filter(!is.na(SSI_wsc_regulated))

cat(sprintf("  SSI stream (regulated):  %d records (%s – %s)\n",
            nrow(ssi_wsc_regulated_long), min(ssi_wsc_regulated_long$date), max(ssi_wsc_regulated_long$date)))

## -- Optional diagnostic: carry over H7's naturalized-splice homogeneity flag --
if (file.exists(cfg$wsc_homogeneity_rds)) {
  wsc_homog <- readRDS(cfg$wsc_homogeneity_rds)
  hg_08jc001 <- tryCatch(wsc_homog[["08JC001"]], error = function(e) NULL)
  if (!is.null(hg_08jc001) && isTRUE(hg_08jc001$tested)) {
    if (isTRUE(hg_08jc001$significant)) {
      message(sprintf(
        "  \u26A0 SSI_wsc_natural (naturalized 08JC001) flagged INHOMOGENEOUS by Pettitt test (change year ~%s, p=%.4f).",
        hg_08jc001$change_year, hg_08jc001$p_value))
      message("    -> Treat SSI_wsc_natural as a potentially spliced/non-stationary series; interpret its vine coupling with caution.")
    } else {
      message("  \u2714 SSI_wsc_natural (naturalized 08JC001) homogeneity check: no significant splice detected.")
    }
  }
} else {
  message("  (H7 homogeneity check RDS not found — skipping naturalized-splice diagnostic.)")
}

# ── 4e. NTSDI Nechako Reservoir ───────────────────────────────────────────
ntsdi_raw  <- read_csv(cfg$ntsdi_file, show_col_types = FALSE)
ntsdi_long <- ntsdi_raw |>
  mutate(date = floor_date(as.Date(date), "month")) |>
  dplyr::select(date, NTSDI = NTSDI_norm) |>
  filter(!is.na(NTSDI))
cat(sprintf("  NTSDI:          %d records (%s – %s)\n",
            nrow(ntsdi_long), min(ntsdi_long$date), max(ntsdi_long$date)))

# ── 4f. Merge all six ─────────────────────────────────────────────────────
all_raw <- snow_long |>
  rename(Snow = all_of(snow_col)) |>
  full_join(ssi_soil_long,          by = "date") |>
  full_join(sgi_long,               by = "date") |>
  full_join(ssi_wsc_natural_long,   by = "date") |>
  full_join(ssi_wsc_regulated_long, by = "date") |>
  full_join(ntsdi_long,             by = "date") |>
  arrange(date)

cat(sprintf("\n  ➤ Combined record: %d months (%s – %s)\n",
            nrow(all_raw), min(all_raw$date), max(all_raw$date)))

# ── Date-alignment diagnostic: flag any column with unexpected NAs ─────────
na_counts <- colSums(is.na(all_raw[, c("Snow", "SSI_soil", "SGI",
                                       "SSI_wsc_natural", "SSI_wsc_regulated", "NTSDI")]))
if (any(na_counts > 0)) {
  message("  \u26A0 NA counts in merged table (may indicate date-format misalignment):")
  for (nm in names(na_counts)[na_counts > 0])
    message(sprintf("      %-18s  %d NAs  (%.1f%% of %d months)",
                    nm, na_counts[nm], 100 * na_counts[nm] / nrow(all_raw), nrow(all_raw)))
} else {
  message("  \u2714 All six index columns fully populated — no date-alignment gaps detected.")
}
cat("\n")

# ==============================================================================
# 5.  REFERENCE-PERIOD CALIBRATED PSEUDO-OBSERVATIONS
# ==============================================================================
cat(strrep("=", 60), "\n  CALIBRATING MARGINAL ECDF TRANSFORMS\n",
    strrep("=", 60), "\n\n", sep = "")

# Variable order is fixed here and must match dvine_order_* indices in cfg:
#   Snow=1  SSI_soil=2  SGI=3  SSI_wsc_natural=4  SSI_wsc_regulated=5  NTSDI=6
var_list   <- c("Snow", "SSI_soil", "SGI", "SSI_wsc_natural", "SSI_wsc_regulated", "NTSDI")
ecdf_cache <- list()

for (v in var_list) {
  ecdf_cache[[v]] <- calibrate_marginals(all_raw, v, cfg$ref_start, cfg$ref_end)
  u_col <- paste0("u_", v)
  all_raw[[u_col]] <- apply_marginals(all_raw, v, ecdf_cache[[v]])
  n_u <- sum(!is.na(all_raw[[u_col]]))
  cat(sprintf("  %-18s  %d uniform values (%.1f%% coverage)\n",
              v, n_u, 100 * n_u / nrow(all_raw)))
}

# ==============================================================================
# 6.  TAIL-DEPENDENCE DIAGNOSTICS (all 15 pairs — 6 variables)
# ==============================================================================
cat("\n", strrep("=", 60), "\n  TAIL-DEPENDENCE MATRIX (λ_L / λ_U)\n",
    strrep("=", 60), "\n\n", sep = "")

u_cols <- paste0("u_", var_list)

ref_u <- all_raw |>
  filter(date >= cfg$ref_start, date <= cfg$ref_end) |>
  dplyr::select(all_of(u_cols)) |>
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
fit_vine_pair <- function(u_mat_fit, dvine_ord_A, dvine_ord_B, label) {
  
  # ── Data-driven D-vine order (Variant C) — empirical tau greedy chaining ──
  dvine_ord_C <- order_dvine_by_tau(u_mat_fit)
  cat(sprintf("\n  Data-driven D-vine order [%s]: %s\n",
              label, paste(colnames(u_mat_fit)[dvine_ord_C], collapse = " -> ")))
  
  # ── D-vine variant A (revised physical cascade) ────────────────────────────
  cat(sprintf("  Fitting D-vine-A (revised physical cascade) [%s] …  ", label))
  vine_dA <- vinecop(
    u_mat_fit,
    structure  = dvine_structure(dvine_ord_A),
    family_set = cfg$family_set,
    par_method = cfg$par_method,
    selcrit    = cfg$selcrit,
    cores      = max(1L, parallel::detectCores() - 1L)
  )
  cat(sprintf("logLik=%.1f  BIC=%.1f\n", logLik(vine_dA), BIC(vine_dA)))
  
  # ── D-vine variant B (revised physical cascade, alt managed tail) ─────────
  cat(sprintf("  Fitting D-vine-B (revised physical cascade, alt tail) [%s] …  ", label))
  vine_dB <- vinecop(
    u_mat_fit,
    structure  = dvine_structure(dvine_ord_B),
    family_set = cfg$family_set,
    par_method = cfg$par_method,
    selcrit    = cfg$selcrit,
    cores      = max(1L, parallel::detectCores() - 1L)
  )
  cat(sprintf("logLik=%.1f  BIC=%.1f\n", logLik(vine_dB), BIC(vine_dB)))
  
  # ── D-vine variant C (fully data-driven tau order) ─────────────────────────
  cat(sprintf("  Fitting D-vine-C (data-driven tau order) [%s] …  ", label))
  vine_dC <- vinecop(
    u_mat_fit,
    structure  = dvine_structure(dvine_ord_C),
    family_set = cfg$family_set,
    par_method = cfg$par_method,
    selcrit    = cfg$selcrit,
    cores      = max(1L, parallel::detectCores() - 1L)
  )
  cat(sprintf("logLik=%.1f  BIC=%.1f\n", logLik(vine_dC), BIC(vine_dC)))
  
  # Select the best D-vine ordering among A / B / C by BIC
  bics <- c(A = BIC(vine_dA), B = BIC(vine_dB), C = BIC(vine_dC))
  dvine_winner <- names(which.min(bics))
  vine_d <- switch(dvine_winner, A = vine_dA, B = vine_dB, C = vine_dC)
  
  cat(sprintf("  D-vine order selected: Variant %s (BIC=%.1f; others: %s)\n",
              dvine_winner, min(bics),
              paste(sprintf("%s=%.1f", names(bics)[names(bics) != dvine_winner],
                            bics[names(bics) != dvine_winner]), collapse = ", ")))
  cat("  NOTE: Variants A/B = revised physical cascade (natural streamflow adjacent\n")
  cat("        to Snow/soil; managed tail order BIC-selected). Variant C = fully\n")
  cat("        empirical tau-greedy order, independent of physical assumptions.\n")
  
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
  which_best <- if (BIC(vine_d) <= BIC(vine_r)) "D-vine" else "R-vine"
  vine_best  <- if (which_best == "D-vine") vine_d else vine_r
  cat(sprintf("  BIC-selected model: %s",  which_best))
  if (which_best == "D-vine")
    cat("  <- note: R-vine recommended for managed variables despite BIC\n")
  else
    cat("  OK (consistent with R-vine recommendation for managed system)\n")
  
  list(
    dA            = vine_dA,
    dB            = vine_dB,
    dC            = vine_dC,
    dvine_order_C = dvine_ord_C,
    d             = vine_d,         # best-BIC D-vine variant among A/B/C
    dvine_tail    = dvine_winner,   # "A", "B", or "C"
    r             = vine_r,         # PRIMARY recommended model
    best          = vine_best,
    which         = which_best,
    trunc         = trunc_bic
  )
}

# ── 8B. Snow season (Oct–Apr, 6 variables) ────────────────────────────────
cat("\n", strrep("=", 60),
    "\n  VINE FITTING — SNOW SEASON Oct–Apr (6 variables)\n",
    strrep("=", 60), "\n", sep = "")

snow_ref <- all_raw |>
  filter(date >= cfg$ref_start, date <= cfg$ref_end,
         month(date) %in% cfg$snow_months) |>
  dplyr::select(all_of(u_cols)) |>
  drop_na_logged("Snow-season reference period (6-var)")

if (nrow(snow_ref) < 40L)
  warning("Fewer than 40 complete snow-season months — estimates unreliable.")

u_snow <- as.matrix(snow_ref)
colnames(u_snow) <- var_list

vines_snow <- fit_vine_pair(u_snow, cfg$dvine_order_snow_A, cfg$dvine_order_snow_B, "Snow-season")
saveRDS(vines_snow, file.path(cfg$out_dir, "vines_snow.rds"))

# ── 8C. Warm season (May–Sep, 5 variables — Snow excluded) ───────────────
cat("\n", strrep("=", 60),
    "\n  VINE FITTING — WARM SEASON May–Sep (5 variables)\n",
    strrep("=", 60), "\n", sep = "")

u_cols_warm    <- paste0("u_", c("SSI_soil", "SGI", "SSI_wsc_natural", "SSI_wsc_regulated", "NTSDI"))
var_list_warm  <- c("SSI_soil", "SGI", "SSI_wsc_natural", "SSI_wsc_regulated", "NTSDI")

warm_ref <- all_raw |>
  filter(date >= cfg$ref_start, date <= cfg$ref_end,
         month(date) %in% cfg$warm_months) |>
  dplyr::select(all_of(u_cols_warm)) |>
  drop_na_logged("Warm-season reference period (5-var)")

if (nrow(warm_ref) < 30L)
  warning("Fewer than 30 complete warm-season months — estimates unreliable.")

u_warm <- as.matrix(warm_ref)
colnames(u_warm) <- var_list_warm

vines_warm <- fit_vine_pair(u_warm, cfg$dvine_order_warm_A, cfg$dvine_order_warm_B, "Warm-season")
saveRDS(vines_warm, file.path(cfg$out_dir, "vines_warm.rds"))

# ==============================================================================
# ==== SECTION 9: VINE CDF EVALUATION ====
# ==============================================================================
set.seed(40L)
cat("\n", strrep("=", 60), "\n  EVALUATING VINE CDF (full record)\n",
    strrep("=", 60), "\n\n", sep = "")

all_snow <- all_raw |>
  filter(month(date) %in% cfg$snow_months) |>
  dplyr::select(date, all_of(u_cols)) |>
  drop_na_logged("Snow-season full record (CDF evaluation)")

all_warm <- all_raw |>
  filter(month(date) %in% cfg$warm_months) |>
  dplyr::select(date, all_of(u_cols_warm)) |>
  drop_na_logged("Warm-season full record (CDF evaluation)")

p_snow_best <- pvinecop(
  `colnames<-`(as.matrix(all_snow[, u_cols]), var_list),
  vines_snow$best, n_mc = cfg$n_mc
)

p_warm_best <- pvinecop(
  `colnames<-`(as.matrix(all_warm[, u_cols_warm]), var_list_warm),
  vines_warm$best, n_mc = cfg$n_mc
)

# ==============================================================================
# ==== SECTION 10: KENDALL DISTRIBUTION ESTIMATION VIA MONTE-CARLO ====
# ==============================================================================
set.seed(40L)
cat("\n  Estimating Kendall distributions (MC, n =",
    format(cfg$n_kc_sim, big.mark = ","), ") …\n")

KC_snow   <- estimate_kendall_cdf(vines_snow$best, cfg$n_kc_sim, cfg$n_mc)
KC_warm   <- estimate_kendall_cdf(vines_warm$best, cfg$n_kc_sim, cfg$n_mc)
cat("  Kendall distributions estimated.\n\n")

# ==============================================================================
# 11. MSDI COMPUTATION — SEASONAL ONLY (OPTION A)
# ==============================================================================
MSDI_snow <- vine_to_msdi(p_snow_best, KC_snow)
MSDI_warm <- vine_to_msdi(p_warm_best, KC_warm)

seasonal_out <- bind_rows(
  tibble(date = all_snow$date, MSDI_seasonal = MSDI_snow,
         season = "Snow (Oct-Apr)"),
  tibble(date = all_warm$date, MSDI_seasonal = MSDI_warm,
         season = "Warm (May-Sep)")
) |> arrange(date)

# ==============================================================================
# ==== SECTION 12: COMPONENT CONTRIBUTIONS - SEASONAL  ====
# ==============================================================================
set.seed(40L)

cat(strrep("=", 60),
    "\n  COMPUTING COMPONENT CONTRIBUTIONS\n",
    strrep("=", 60), "\n\n", sep = "")

compute_contributions <- function(u_full_mat, vine_obj, var_names, n_mc) {
  p_base <- pvinecop(u_full_mat, vine_obj, n_mc = n_mc)
  contrib <- matrix(NA_real_, nrow = nrow(u_full_mat), ncol = length(var_names),
                    dimnames = list(NULL, paste0("contrib_", var_names)))
  
  for (k in seq_along(var_names)) {
    u_neu      <- u_full_mat
    u_neu[, k] <- 0.50
    p_neu      <- pvinecop(u_neu, vine_obj, n_mc = n_mc)
    contrib[, k] <- p_neu - p_base
  }
  
  row_sum <- rowSums(abs(contrib))
  contrib_norm <- sweep(contrib, 1, ifelse(row_sum == 0, 1, row_sum), "/")
  as.data.frame(contrib_norm)
}

contrib_snow <- compute_contributions(
  as.matrix(all_snow[, u_cols]), vines_snow$best, var_list,
  n_mc = min(5000L, cfg$n_mc)
)
contrib_warm <- compute_contributions(
  as.matrix(all_warm[, u_cols_warm]), vines_warm$best, var_list_warm,
  n_mc = min(5000L, cfg$n_mc)
)

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
  d          = c(6L, 6L, 5L, 5L),
  logLik     = c(logLik(vines_snow$d),  logLik(vines_snow$r),
                 logLik(vines_warm$d),  logLik(vines_warm$r)),
  AIC        = c(AIC(vines_snow$d),  AIC(vines_snow$r),
                 AIC(vines_warm$d),  AIC(vines_warm$r)),
  BIC        = c(BIC(vines_snow$d),  BIC(vines_snow$r),
                 BIC(vines_warm$d),  BIC(vines_warm$r))
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
cat("    (skipped — rvinecopulib does not export as_VineCopula; use VineCopula::RVineGOFTest with manual conversion if needed)\n")

# ── 14d. D-vine order diagnostic — which variant won, per season ─────────
cat("\n  D-vine order selection summary:\n")
cat(sprintf("    Snow season: Variant %s selected  (data-driven order = %s)\n",
            vines_snow$dvine_tail,
            paste(var_list[vines_snow$dvine_order_C], collapse = " -> ")))
cat(sprintf("    Warm season: Variant %s selected  (data-driven order = %s)\n",
            vines_warm$dvine_tail,
            paste(var_list_warm[vines_warm$dvine_order_C], collapse = " -> ")))

# ==============================================================================
# 15. ASSEMBLE MAIN OUTPUT TABLE — SEASONAL MSDI PRIMARY
# ==============================================================================
msdi_out <- seasonal_out |>
  left_join(contrib_df, by = c("date", "season")) |>
  left_join(all_raw |> dplyr::select(date, all_of(var_list)), by = "date") |>
  mutate(
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
  boot_msdi_snow <- matrix(NA_real_, nrow = nrow(all_snow), ncol = cfg$n_boot)
  boot_msdi_warm <- matrix(NA_real_, nrow = nrow(all_warm), ncol = cfg$n_boot)
  
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

p1 <- plot_standard_basin_ts(df_msdi, index_name = "MSDI", color = "#d73027") +
  ggplot2::facet_wrap(~season, ncol = 1, scales = "free_x")

ggsave(file.path(cfg$out_dir, "plot1_MSDI_seasonal.png"),
       p1, width = 13, height = 7, dpi = 300)
cat("  Plot 1: MSDI seasonal time series\n")

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
if (any(grepl("^contrib_", names(contrib_df)))) {
  contrib_long <- contrib_df |>
    dplyr::select(date, season, starts_with("contrib_")) |>
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
snow_dates <- all_raw |>
  filter(date >= cfg$ref_start, date <= cfg$ref_end,
         month(date) %in% cfg$snow_months) |>
  dplyr::select(date, all_of(u_cols)) |>
  drop_na_logged("Plot 5 snow-season date index") |>
  pull(date)

MSDI_snow_plot <- seasonal_out |>
  filter(date %in% snow_dates) |>
  pull(MSDI_seasonal)

pd5 <- as_tibble(u_snow) |>
  setNames(var_list) |>
  mutate(MSDI = MSDI_snow_plot,
         Drought = cut(MSDI, c(-Inf, -2, -1.5, -1, 0, Inf),
                       labels = c("Extreme", "Severe", "Moderate", "Minor", "Normal")))

if (requireNamespace("GGally", quietly = TRUE)) {
  library(GGally)
  p5  <- ggpairs(
    pd5, columns = seq_along(var_list),
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
         p5, width = 11, height = 10, dpi = 250)
  cat("  Plot 5: Pair scatter matrix (snow season)\n")
} else {
  cat("  Plot 5 skipped (GGally not installed)\n")
}

# ── Plot 6: Drought event severity vs duration scatter ───────────────────
if (nrow(drought_minor) > 0) {
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
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
  theme(axis.text.x = element_text(size = 6.5),
        plot.title = element_text(face = "bold"),
        legend.position = "top")

ggsave(file.path(cfg$out_dir, "plot7_tail_dependence.png"),
       p7, width = 14, height = 4.5, dpi = 300)
cat("  Plot 7: Tail-dependence bar chart\n")

# ==============================================================================
# 18. FINAL SUMMARY — SEASONAL ONLY (OPTION A)
# ==============================================================================
cat("\n", strrep("=", 60),
    "\n  VINE COPULA MSDI — COMPLETE (v2: natural/regulated WSC split +\n",
    "  empirical tau-informed D-vine ordering)\n",
    strrep("=", 60), "\n\n", sep = "")

cat(sprintf("  Records (seasonal MSDI):  %d months (%s – %s)\n",
            nrow(seasonal_out), min(seasonal_out$date), max(seasonal_out$date)))

best_snow_label <- vines_snow$which
best_warm_label <- vines_warm$which

cat(sprintf("  Snow-season best model:    %s (BIC=%.1f)\n",
            best_snow_label,
            min(BIC(vines_snow$d), BIC(vines_snow$r))))
cat(sprintf("  Warm-season best model:    %s (BIC=%.1f)\n",
            best_warm_label,
            min(BIC(vines_warm$d), BIC(vines_warm$r))))
cat(sprintf("  Drought events catalogued: %d moderate/severe\n",
            nrow(drought_moderate)))
cat(sprintf("  Outputs in:                %s/\n\n", cfg$out_dir))

cat("  Output files:\n")
for (f in list.files(cfg$out_dir, full.names = FALSE))
  cat("    •", f, "\n")

# ==============================================================================
# END
# ==============================================================================