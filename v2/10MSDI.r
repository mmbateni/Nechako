# ==============================================================================
# 10MSDI_pipeline.R - NECHAKO BASIN - MSDI-N index, relationship analysis and
#                     README diagrams in one script
#
# Merged from (run order kept):
#   Part A  10MSDI.r         vine-copula MSDI-N (v5)       -> msdi_results_<VERSION_TAG>/
#   Part B  diags_10.R       relationship analysis (A2 v5)  -> msdi_results_<VERSION_TAG>/relationship_v5/
#   Part C  make_diagrams.R  README diagrams (Graphviz)     -> <DIAG_README_DIR>/figures/
#
# Part B reuses Part A's inputs, reference period and Gringorten pseudo-
# observations, so the two can no longer drift apart. Its consistency check
# compares against the relationship_by_regime values computed in this run.
# Part C edge labels remain manually curated (some source screening tables are
# external). A source-hash audit and explicit warning prevent silent staleness;
# the diagram labels are NOT automatically certified against numerical results.
#
# Run from the project root (the folder holding DROUGHT_ANALYSIS_utils.R).
# R >= 4.1 (uses \(x) lambdas and r"(...)" raw strings).
# Sections 1-6 always run (Part B needs them); switch parts on/off in RUN.
# ==============================================================================
rm(list = ls())
t_start <- Sys.time()

RUN <- list(
  msdi         = TRUE,   # Part A, sections 7-13: vines, MSDI-N, events, plots, summary
  relationship = TRUE,   # Part B: tau by month, partial correlations, cascades, mediation
  diagrams     = TRUE    # Part C: README Graphviz diagrams
)
DIAG_README_DIR <- getwd()   # EDIT: folder holding the README (Part C writes figures/ there)

# ==============================================================================
# PART A - VINE-COPULA MULTIVARIATE STANDARDISED DROUGHT INDEX (10MSDI.r)
#
# v5 - VERIFIED NATURALIZED FLOW TO 2025
#   (1) Naturalized 08JC002 SSI (Script 6, Pipeline B) is verified to 2025-12.
#       The Flow mask and the SRI-1 provisional splice are removed: MSDI-N uses
#       naturalized flow for every month 1956-01..2025-12. 
#       However, since WSC not approved the last tree monthes yet, those data is still provisional.
#   (2) Vine fit AND marginal reference extended to 1956-2025 (70 years).
#       Set cfg$fit_end <- as.Date("2021-12-01") to reproduce the v4 reference
#       (2022-2025 are then evaluated out of sample).
#   (3) Structure criterion: effective-sample-size BIC (cfg$bic_mode = "ess"),
#         -2 logL (n_years / n_rows) + k log(n_years)
#       "rows" = standard BIC. v4 combined the row log-likelihood with a years
#       penalty, which was less conservative than either.
#   (4) Monte-Carlo resolution: a month whose C(u) is 0 at tail-sample
#       resolution, while the Kendall sample also has ties at 0, is now flagged
#       (mc_floor). C(u) is written (C_u) and the lowest attainable MSDI per
#       window is reported (msdi_resolution_floor, model_comparison.csv).
#   (5) SRI-1 (ERA5-Land) is kept only as an optional near-real-time variant /
#       sensitivity run (cfg$run_sri_variant): separate outputs, agreement with
#       the main index overall and by regime. Never spliced into MSDI-N.
#   (6) Input_status from Script 6 v5.4 is carried through; months with
#       provisional/persisted daily inputs are flagged provisional = TRUE
#       (none expected with the verified files).
#   (7) Bivariate (recession) windows fit one copula (no candidate contest);
#       known-cosmetic GGally warnings in plot 5 are muffled.
#
# v4 (retained): natural-state scope (no storage / regulated-flow / SGI nodes);
#   nodes Supply SPEI-PM-3, Snow SSPI-1, Soil SSI-SM-1 L1_3, Flow naturalized
#   SSI 08JC002; regimes and design D-vine orders
#     Accumulation Nov-Apr : Snow - Supply - Flow - Soil
#     Melt         May-Jun : Supply - Soil - Flow - Snow
#     Recession    Jul-Oct : Soil - Flow
#   one vine per calendar month on a +/-1-month window within the regime;
#   Gringorten margins per calendar month; Kendall transform by Monte Carlo;
#   stationarity tests at 1977 / 1989; MSDI_seasonal / season names kept for
#   Script 14. Design rationale: README_Multivariate_Drought_Framework(H_Copula).md
#
# v3 (retained): framework contract (project-root WD, DROUGHT_ANALYSIS_utils.R,
#   Gringorten rule, identify_drought_events_ts() events, McKee classes,
#   input manifest, Year x Jan..Dec wide table).
#
# COMPONENT INDICES
# ------------------------------------------------------------------------------
#  Node    Index              Script source                         Compartment
#  ------  -----------------  ------------------------------------  -------------
#  Supply  SPEI-PM-3          SPI/SPEI script (basin mean)          P - PET supply
#  Snow    SSPI-1 / SWEI-3    4snow_drought_indices_ERALand.R       Snowpack
#  Soil    SSI-SM-1 L1_3      8SSI_SSMI_ERALand.R                   Root-zone soil
#  Flow    SSI 08JC002 (B)    6StreamFlow_Index_WSCLocal.R Pipe B   Naturalized flow
#  (SRI1   SRI-1              Script 5, ERA5-Land runoff            optional variant)
#
# MSDI = PHI^-1( K_C( C(u) ) )  (Kao & Govindaraju 2010), per calendar month.
#
# OUTPUTS (msdi_results_<VERSION_TAG>/):
#   Nechako_MSDI_VineCopula.csv     monthly MSDI + regime + flags + C_u + contributions
#   msdi_01_by_month.csv            Year x Jan..Dec
#   msdi_drought_events.csv         framework events + primary driver
#   msdi_input_manifest.csv         inputs, sources, records
#   relationship_by_regime.csv      tau and tail co-occurrence per regime pair
#   stationarity_tests.csv          tree-1 pairs, 1977 / 1989 splits
#   model_comparison.csv            per-month BIC, selected structure, MC resolution
#   vine_edges.csv                  every pair-copula (family, rotation, tau)
#   msdi_normality_by_regime.csv
#   msdi_sri1_variant.csv, sri1_variant_agreement.csv   (if run_sri_variant)
#   MSDI_bootstrap_CI.csv (optional), msdi_summary.txt, figures/*.png, plot3..plot7_*.png
#
# EXECUTION ORDER: 4, 5, 6, 8 (and the SPI/SPEI script) -> this script -> 14
#
# REFERENCES
#   Sklar (1959); Joe (2015) Dependence Modeling with Copulas
#   Aas et al. (2009) Insur. Math. Econ. 44, 182-198 - pair-copula construction
#   Bedford & Cooke (2002) Ann. Stat. 30, 1031-1068 - R-vine
#   Dissmann et al. (2013) Comput. Stat. Data Anal. 59, 52-69 - vine selection
#   Genest & Favre (2007) J. Hydrol. Eng. 12, 347-368
#   Kao & Govindaraju (2010) J. Hydrol. 380, 121-134 - Kendall-distribution JDI
#   Hao & AghaKouchak (2013) Adv. Water Resour. 57, 12-18 - MSDI
#   Vicente-Serrano et al. (2010) J. Climate 23, 1696-1718 - SPEI
#   Huning & AghaKouchak (2020) Nat. Commun. 11, 965 - SWEI
#   Pettitt (1979) Appl. Stat. 28, 126-135
#   Kass & Raftery (1995) J. Am. Stat. Assoc. 90, 773-795 - BIC evidence scale
# ==============================================================================

# ---- Working directory (framework rule) -------------------------------------
WD_PATH    <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
UTILS_FILE <- file.path(WD_PATH, "DROUGHT_ANALYSIS_utils.R")
if (!file.exists(UTILS_FILE))
  stop("DROUGHT_ANALYSIS_utils.R not found in the working directory: ", WD_PATH,
       "\nSet the working directory to the project root before running this script.")
source(UTILS_FILE)

# ==============================================================================
# 1.  PACKAGES
# ==============================================================================
utils_load_packages(c("ggplot2", "dplyr", "tidyr", "readr", "purrr", "tibble",
                      "lubridate", "rvinecopulib","DiagrammeR","DiagrammeRsvg","rsvg","ggrepel","GGally"))

MSDI_EVENT_SCALE <- 1L
DATE_NEG_INF <- structure(-Inf, class = "Date")
DATE_POS_INF <- structure( Inf, class = "Date")
ALL_NODES    <- c("Supply", "Snow", "Soil", "Flow")

# ==============================================================================
# 2.  USER CONFIGURATION
# ==============================================================================
SWEI_DIR <- if (exists("SWEI_SEAS_DIR")) SWEI_SEAS_DIR else file.path(WD_PATH, "swei_results_seasonal")
SSPI_DIR <- file.path(WD_PATH, "sspi_results_seasonal")
VERSION_TAG <- "N5_regime"

cfg <- list(
  # -- Periods ---------------------------------------------------------------
  analysis_start = as.Date("1956-01-01"),   # first full year of naturalized flow
  analysis_end   = as.Date("2025-12-01"),
  fit_start      = as.Date("1956-01-01"),   # vine fit AND marginal reference (also Part B period)
  fit_end        = as.Date("2025-12-01"),   # v5: full verified record (v4: 2021-12-01)
  min_ref_n      = 5L,
  
  # -- Supply node: SPEI-PM-3 basin mean --------------------------------------
  # EDIT to the SPEI-PM output of your SPI/SPEI script. Layouts accepted:
  # wide (Year + Jan..Dec) or long (date + value column named in supply_cols).
  # NOTE: diags_10.R pointed to spei_results/ instead of spei_results_seasonal/;
  # the merged script uses this one path for Parts A and B.
  supply_file = file.path(WD_PATH, "spei_results_seasonal", "spei_03_basin_averaged_by_month.csv"),
  supply_cols = c("spei_03", "SPEI_03", "SPEI_3", "spei_pm_03", "spei", "SPEI"),
  
  # -- Snow node: "sspi" (SSPI-1, default) or "swei" (SWEI-3, sensitivity) ----
  snow_index = "sspi",
  swei_file  = file.path(SWEI_DIR, "swei_03_basin_averaged_by_month.csv"),
  sspi_file  = file.path(SSPI_DIR, "sspi_01_basin_averaged_by_month.csv"),
  
  # -- Soil node: SSI-SM (Beta) L1_3, 1-month ---------------------------------
  soil_index = "ssi_sm",
  soil_layer = "L1_3",
  soil_scale = 1L,
  
  # -- Flow node: naturalized SSI, Pipeline B (verified to 2025-12) -----------
  flow_station = "08JC002",
  flow_rds     = file.path(WD_PATH, "streamflow_results", "ssi", "monthly_data",
                           "08JC002_monthly_ssi.rds"),
  provisional_status = c("provisional", "persisted"),   # Script 6 Input_status values
  provisional_manual = as.Date(c("2025-10-01", "2025-11-01", "2025-12-01")),  # README 2.3; empty once approved
  
  # -- Optional near-real-time variant: ERA5-Land SRI-1 in the Flow slot ------
  run_sri_variant = TRUE,
  sri_file = file.path(WD_PATH, "sri_results_seasonal", "runoff", "sri_01_basin_averaged_by_month.csv"),
  sri_cols = c("SRI_1", "sri_01", "SRI_01", "SRI1"),
  
  # -- Fallback: the index panel used for the relationship analysis -----------
  panel_file = file.path(WD_PATH, "index_panel_monthly.csv"),
  allow_panel_fallback = FALSE,   # single authoritative panel fallback switch
  # -- Hydrological regimes (node set + design D-vine order) ------------------
  regimes = list(
    acc  = list(label = "Accumulation (Nov-Apr)", months = c(11L, 12L, 1L, 2L, 3L, 4L),
                order = c("Snow", "Supply", "Flow", "Soil")),
    melt = list(label = "Melt (May-Jun)", months = c(5L, 6L),
                order = c("Supply", "Soil", "Flow", "Snow")),
    rec  = list(label = "Recession (Jul-Oct)", months = 7:10,
                order = c("Soil", "Flow"))
  ),
  window_within_regime = TRUE,   # Part A AND Part B +/-1-month windows respect regime
  
  # -- Vine settings ----------------------------------------------------------
  family_set = "parametric",     # rotations allowed (negative Snow-Soil edge)
  selcrit    = "bic",
  par_method = "mle",
  bic_mode   = "ess",            # "ess" (one observation per year) or "rows" (standard)
  bic_margin = 10,               # BIC gain needed to leave the design D-vine
  
  # -- Monte-Carlo Kendall engine ---------------------------------------------
  n_ref      = 1e5L,             # reference sample for C(u)
  n_ref_tail = 1e6L,             # dry-tail reference sample
  n_eval     = 2e4L,             # sample for the Kendall CDF
  tail_count = 20L,              # points with <= this many dominated -> tail sample
  prob_eps   = 1e-6,
  
  # -- Diagnostics ------------------------------------------------------------
  tail_q            = c(0.10, 0.20),
  stationarity_cuts = c(1977L, 1989L),
  n_perm            = 999L,
  
  # -- Bootstrap (year blocks, same structure) --------------------------------
  run_bootstrap   = FALSE,
  n_boot          = 100L,
  boot_n_ref      = 2e4L,
  boot_n_ref_tail = 2e5L,
  boot_n_eval     = 5e3L,
  # -- Plots ------------------------------------------------------------------
  plot_label_top   = 10L,                    # event scatter: label the n most severe
  plot_label_extra = c("2000-12", "2025-06"),# ... plus these event starts (YYYY-MM)
  plot_zoom_n      = 1L,                     # contribution close-ups for the n most severe events
  plot_zoom_pad    = 2L,                     # months shown either side of a close-up event
  heatmap_binned   = FALSE,                  # TRUE: McKee class colours instead of a gradient
  pairs_regimes    = c("acc", "melt"),       # regimes with a pairs plot (>= 3 nodes)
  
  seed    = 40L,
  out_dir = file.path(WD_PATH, paste0("msdi_results_", VERSION_TAG))
)

reg_months <- unlist(lapply(cfg$regimes, `[[`, "months"))
stopifnot(cfg$snow_index %in% c("swei", "sspi"),
          cfg$soil_index %in% c("ssi_sm", "ssmi"),
          cfg$bic_mode %in% c("ess", "rows"),
          setequal(reg_months, 1:12), !anyDuplicated(reg_months),
          all(unlist(lapply(cfg$regimes, `[[`, "order")) %in% ALL_NODES),
          cfg$fit_start >= cfg$analysis_start, cfg$fit_end <= cfg$analysis_end)

cfg$soil_file <- file.path(WD_PATH, sprintf("%s_results_seasonal", cfg$soil_index), cfg$soil_layer,
                           sprintf("%s_%02d_basin_averaged_by_month.csv", cfg$soil_index, cfg$soil_scale))
FIG_DIR <- file.path(cfg$out_dir, "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
REF_LABEL <- sprintf("%s-%s", format(cfg$fit_start, "%Y"), format(cfg$fit_end, "%Y"))

# ==============================================================================
# 3.  HELPER FUNCTIONS
# ==============================================================================
month_start <- function(d) as.Date(format(as.Date(d), "%Y-%m-01"))

## 3a. Wide basin CSV (Year + Jan..Dec) -> long -------------------------------
pivot_wide_monthly <- function(path, value_name) {
  if (!file.exists(path)) stop("Input not found: ", path)
  df <- readr::read_csv(path, show_col_types = FALSE)
  yr_col <- intersect(c("Year", "year", "YEAR"), names(df))[1]
  mcols  <- intersect(month.abb, names(df))
  if (is.na(yr_col) || length(mcols) != 12L)
    stop("Expected a Year column and Jan..Dec columns in ", path)
  out <- df |>
    dplyr::select(Year = dplyr::all_of(yr_col), dplyr::all_of(mcols)) |>
    tidyr::pivot_longer(-Year, names_to = "month_abbr", values_to = value_name) |>
    dplyr::mutate(date = as.Date(sprintf("%d-%02d-01", as.integer(Year), match(month_abbr, month.abb))),
                  !!value_name := suppressWarnings(as.numeric(.data[[value_name]]))) |>
    dplyr::filter(!is.na(date), is.finite(.data[[value_name]]))
  out[, c("date", value_name), drop = FALSE]
}

## 3b. Any layout (wide Year x Jan..Dec, or long date + value) -> long --------
read_index_any <- function(path, value_name, cols = NULL) {
  df <- readr::read_csv(path, show_col_types = FALSE)
  yr_col <- intersect(c("Year", "year", "YEAR"), names(df))[1]
  if (!is.na(yr_col) && all(month.abb %in% names(df)))
    return(pivot_wide_monthly(path, value_name))
  dcol <- intersect(c("date", "Date", "DATE", "month_label", "time"), names(df))[1]
  if (is.na(dcol)) stop("No date column (or Year + Jan..Dec layout) in ", path)
  vcol <- if (length(cols)) intersect(cols, names(df))[1] else NA_character_
  if (is.na(vcol)) {
    num <- setdiff(names(df)[vapply(df, is.numeric, logical(1))], c(dcol, yr_col))
    if (length(num) != 1L)
      stop(sprintf("Cannot identify the value column in %s (numeric columns: %s). Set the *_cols option.",
                   path, paste(num, collapse = ", ")))
    vcol <- num
  }
  dv <- df[[dcol]]
  dd <- if (inherits(dv, c("Date", "POSIXt"))) as.Date(dv)
  else as.Date(paste0(substr(as.character(dv), 1, 7), "-01"))
  out <- data.frame(date = month_start(dd), v = suppressWarnings(as.numeric(df[[vcol]])))
  names(out)[2] <- value_name
  out[!is.na(out$date) & is.finite(out[[value_name]]), , drop = FALSE]
}

## 3c. Streamflow SSI RDS (Script 6) -> long, keeping Input_status ------------
read_ssi_rds <- function(f) {
  x <- readRDS(f); nm <- names(x)
  if (!all(c("date", "ssi") %in% nm)) stop("Expected date and ssi columns in ", f)
  scol <- nm[tolower(nm) == "input_status"][1]          # Input_status or input_status
  if (is.na(scol)) message("  Flow RDS has no input_status column; provisional flags come from cfg$provisional_manual only.")
  data.frame(date = month_start(x$date), ssi = as.numeric(x$ssi),
             flow_status = if (!is.na(scol)) tolower(as.character(x[[scol]])) else NA_character_,
             stringsAsFactors = FALSE) |>
    dplyr::filter(is.finite(ssi))
}
## 3d. Input manifest ------------------------------------------------------------
manifest <- list()
log_input <- function(node, index, file, script, df, col) {
  ok <- df[is.finite(df[[col]]), , drop = FALSE]
  cat(sprintf("  %-7s %-30s %5d records (%s - %s)\n", node, index, nrow(ok),
              if (nrow(ok)) format(min(ok$date), "%Y-%m") else "NA",
              if (nrow(ok)) format(max(ok$date), "%Y-%m") else "NA"))
  manifest[[length(manifest) + 1L]] <<- data.frame(
    node = node, index = index, script = script, file = paste(file, collapse = "; "),
    n_records = nrow(ok),
    start = if (nrow(ok)) format(min(ok$date), "%Y-%m") else NA_character_,
    end   = if (nrow(ok)) format(max(ok$date), "%Y-%m") else NA_character_,
    stringsAsFactors = FALSE)
}

## 3e. Load one node: RDS -> file -> panel fallback -----------------------------
load_node <- function(node, label, script, file = NULL, rds = NULL, cols = NULL, panel_col = NULL) {
  if (!is.null(rds) && file.exists(rds)) {
    df <- read_ssi_rds(rds); names(df)[2] <- node; src <- rds
  } else if (!is.null(file) && file.exists(file)) {
    df <- read_index_any(file, node, cols); src <- file
  } else if (isTRUE(cfg$allow_panel_fallback) && !is.null(panel_col) && file.exists(cfg$panel_file)) {    pn <- readr::read_csv(cfg$panel_file, show_col_types = FALSE)
  if (!panel_col %in% names(pn)) stop(sprintf("Column '%s' not in %s", panel_col, cfg$panel_file))
  warning(sprintf("%s: upstream file not found - using column '%s' of %s.",
                  node, panel_col, basename(cfg$panel_file)), call. = FALSE)
  df <- data.frame(date = month_start(pn$date), v = suppressWarnings(as.numeric(pn[[panel_col]])))
  names(df)[2] <- node
  df  <- df[is.finite(df[[node]]), , drop = FALSE]
  src <- paste0(cfg$panel_file, "::", panel_col)
  } else {
    stop(sprintf("No input for node %s (looked for: %s; panel: %s).", node,
                 paste(c(rds, file), collapse = " | "), cfg$panel_file))
  }
  df <- dplyr::distinct(df, date, .keep_all = TRUE)
  log_input(node, label, src, script, df, node)
  df
}

## 3f. Per-calendar-month reference + Gringorten pseudo-observations -------------
calibrate_marginals <- function(df, value_col, ref_start, ref_end, min_n = cfg$min_ref_n) {
  ref <- df[df$date >= ref_start & df$date <= ref_end, , drop = FALSE]
  mon <- lubridate::month(ref$date)
  lapply(1:12, function(m) {
    vals <- sort(ref[[value_col]][mon == m & is.finite(ref[[value_col]])])
    if (length(vals) < min_n) NULL else vals
  })
}
apply_marginals <- function(df, value_col, refs) {
  u    <- rep(NA_real_, nrow(df))
  mons <- lubridate::month(df$date)
  for (m in 1:12) {
    idx <- which(mons == m & is.finite(df[[value_col]]))
    if (!length(idx) || is.null(refs[[m]])) next
    ref <- refs[[m]]; n <- length(ref); x <- df[[value_col]][idx]
    n_lt <- findInterval(x, ref, left.open = TRUE)
    n_le <- findInterval(x, ref)
    r    <- n_lt + (n_le - n_lt + 1) / 2
    u[idx] <- (r - 0.44) / (n + 0.12)
  }
  u
}

## 3g. Regime of a month and its fit window ------------------------------------
regime_of <- function(m)
  names(cfg$regimes)[vapply(cfg$regimes, function(r) m %in% r$months, logical(1))]
fit_window <- function(m) {
  w <- ((m + (-1:1) - 1) %% 12) + 1
  if (cfg$window_within_regime) w <- intersect(w, cfg$regimes[[regime_of(m)]]$months)
  sort(w)
}

# Safe pairwise correlation for descriptive diagnostics: only finite paired data.
# Undefined for fewer than 3 pairs or constant samples (return NA, never NaN).
safe_cor <- function(x, y, method = "pearson") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L) return(NA_real_)
  x <- x[ok]; y <- y[ok]
  if (!is.finite(stats::var(x)) || !is.finite(stats::var(y)) ||
      stats::var(x) <= 0 || stats::var(y) <= 0) return(NA_real_)
  z <- suppressWarnings(stats::cor(x, y, method = method))
  if (is.finite(z)) z else NA_real_
}

## 3h. Vine fitting and autocorrelation-aware BIC ------------------------------
n_cores <- max(1L, parallel::detectCores() - 1L)
fit_one <- function(u, structure = NULL) {
  args <- list(data = u, family_set = cfg$family_set, par_method = cfg$par_method,
               selcrit = cfg$selcrit, cores = n_cores)
  if (!is.null(structure)) args$structure <- structure
  do.call(rvinecopulib::vinecop, args)
}
# "ess" : each year's 2-3 window rows count as one observation (conservative)
# "rows": standard BIC on all rows
bic_years <- function(fit, n_rows, n_years) {
  L <- stats::logLik(fit)
  k <- attr(L, "df"); if (is.null(k)) k <- fit$npars
  L <- as.numeric(L)
  if (identical(cfg$bic_mode, "rows")) -2 * L + k * log(n_rows)
  else -2 * L * n_years / n_rows + k * log(n_years)
}

## 3i. Data-driven D-vine spine via greedy Kendall-tau chaining --------------
order_dvine_by_tau <- function(u_mat) {
  var_names <- colnames(u_mat); d <- ncol(u_mat)
  if (d < 3L) return(seq_len(d))
  tau <- suppressWarnings(stats::cor(u_mat, method = "kendall", use = "pairwise.complete.obs"))
  diag(tau) <- 0
  ut <- upper.tri(tau)
  best <- which(abs(tau) == max(abs(tau[ut]), na.rm = TRUE) & ut, arr.ind = TRUE)[1, ]
  chain <- var_names[best]; remaining <- setdiff(var_names, chain)
  while (length(remaining)) {
    le <- chain[1]; re <- chain[length(chain)]
    sc <- sapply(remaining, function(v) max(abs(tau[v, le]), abs(tau[v, re]), na.rm = TRUE))
    nxt <- remaining[which.max(sc)]
    chain <- if (abs(tau[nxt, le]) >= abs(tau[nxt, re])) c(nxt, chain) else c(chain, nxt)
    remaining <- setdiff(remaining, nxt)
  }
  match(chain, var_names)
}

## 3j. Pair-copula table and independence-edge count ---------------------------
to_idx <- function(i) {
  i <- unlist(i)
  if (is.character(i)) i <- unlist(strsplit(i, "[^0-9]+"))
  i <- suppressWarnings(as.integer(i)); i[!is.na(i)]
}
edge_table <- function(fit, node_names) {
  s <- tryCatch(as.data.frame(summary(fit)), error = function(e) NULL)
  if (is.null(s) || !nrow(s) || !all(c("tree", "conditioned", "family") %in% names(s))) return(NULL)
  lab <- function(x) vapply(x, function(i) paste(node_names[to_idx(i)], collapse = "-"), character(1))
  data.frame(tree = s$tree, edge = s$edge,
             conditioned  = lab(s$conditioned),
             conditioning = if ("conditioning" %in% names(s)) lab(s$conditioning) else "",
             family   = s$family,
             rotation = if ("rotation" %in% names(s)) s$rotation else NA,
             tau      = if ("tau" %in% names(s)) round(as.numeric(s$tau), 3) else NA_real_,
             stringsAsFactors = FALSE)
}
n_indep <- function(fit) {
  s <- tryCatch(as.data.frame(summary(fit)), error = function(e) NULL)
  if (is.null(s)) NA_integer_ else sum(s$family == "indep")
}

## 3k. Monte-Carlo Kendall engine ------------------------------------------------
# C(x) is estimated as the share of a simulated reference sample dominated by x.
# For each query x only the sorted prefix R[, j] <= x_j of its smallest
# coordinate j is scanned, so dry-tail queries are cheap. Queries with few
# dominated points are re-evaluated against a 10x larger tail sample.
make_dom_index <- function(R) {
  R <- as.matrix(R); d <- ncol(R)
  mats <- lapply(seq_len(d), function(j) R[order(R[, j]), , drop = FALSE])
  list(mats = mats, keys = lapply(seq_len(d), function(j) mats[[j]][, j]), n = nrow(R), d = d)
}
dom_count <- function(S, idx) {
  S <- as.matrix(S); out <- integer(nrow(S)); d <- idx$d
  for (i in seq_len(nrow(S))) {
    x <- S[i, ]; j <- which.min(x)
    k <- findInterval(x[j], idx$keys[[j]])
    if (k == 0L) next
    if (d == 1L) { out[i] <- k; next }
    cand <- idx$mats[[j]][seq_len(k), -j, drop = FALSE]
    out[i] <- sum(rowSums(cand <= rep(x[-j], each = k)) == d - 1L)
  }
  out
}
t_hat <- function(S, eng) {
  S  <- as.matrix(S)
  c1 <- dom_count(S, eng$R1)
  t  <- c1 / eng$R1$n
  tl <- which(c1 <= cfg$tail_count)
  if (length(tl)) t[tl] <- dom_count(S[tl, , drop = FALSE], eng$R2) / eng$R2$n
  t
}
kendall_engine <- function(fit, seed, n_ref, n_ref_tail, n_eval) {
  set.seed(seed)
  eng <- list(R1 = make_dom_index(rvinecopulib::rvinecop(n_ref, fit)),
              R2 = make_dom_index(rvinecopulib::rvinecop(n_ref_tail, fit)))
  eng$t_sim <- sort(t_hat(rvinecopulib::rvinecop(n_eval, fit), eng))
  eng
}
K_eval <- function(t, t_sorted) {   # mid-rank Kendall CDF
  (findInterval(t, t_sorted, left.open = TRUE) + findInterval(t, t_sorted)) / (2 * length(t_sorted))
}
k_lo <- function(eng) max(cfg$prob_eps, 0.5 / length(eng$t_sim))
# Lowest attainable MSDI of an engine: a query with C = 0 gets K = n0 / (2 n_eval),
# n0 = simulated points that also have C = 0 at tail-sample resolution.
resolution_floor <- function(eng) {
  n0 <- sum(eng$t_sim <= 0)
  stats::qnorm(max(n0 / (2 * length(eng$t_sim)), k_lo(eng)))
}
msdi_from_u <- function(U, eng) {
  t  <- t_hat(U, eng)
  k  <- K_eval(t, eng$t_sim)
  lo <- k_lo(eng)
  list(t = t, k = k,
       floor = k <= 0 | k >= 1 | (t <= 0 & any(eng$t_sim <= 0)),
       msdi = stats::qnorm(pmin(pmax(k, lo), 1 - lo)))
}

## 3l. Empirical tail co-occurrence: mean over q of P(V <= q | U <= q) --------
tail_dep <- function(u, v, q = cfg$tail_q) {
  lL <- mean(sapply(q, function(qi) { a <- sum(u <= qi); if (a) sum(u <= qi & v <= qi) / a else NA_real_ }), na.rm = TRUE)
  lU <- mean(sapply(q, function(qi) { a <- sum(u >= 1 - qi); if (a) sum(u >= 1 - qi & v >= 1 - qi) / a else NA_real_ }), na.rm = TRUE)
  c(lambda_L = lL, lambda_U = lU)
}

## 3m. Split-sample dependence test (year-block permutation) -------------------
split_test <- function(x, y, yr, cut, B) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]; yr <- yr[ok]
  yrs <- unique(yr); early <- yrs[yrs < cut]
  if (length(early) < 5L || length(yrs) - length(early) < 5L)
    return(c(tau_early = NA, tau_late = NA, diff = NA, p_perm = NA))
  tk  <- function(s) stats::cor(x[s], y[s], method = "kendall")
  g   <- function(e) { s <- yr %in% e; tk(s) - tk(!s) }
  obs <- g(early)
  nul <- replicate(B, g(sample(yrs, length(early))))
  c(tau_early = tk(yr %in% early), tau_late = tk(!(yr %in% early)),
    diff = obs, p_perm = mean(abs(nul) >= abs(obs)))
}

## 3n. McKee et al. (1993) classes ---------------------------------------------
msdi_class <- function(x) {
  dplyr::case_when(
    !is.finite(x) ~ NA_character_,
    x >=  2.0     ~ "W3 Extremely wet",
    x >=  1.5     ~ "W2 Severely wet",
    x >=  1.0     ~ "W1 Moderately wet",
    x >  -1.0     ~ "N0 Near normal",
    x >  -1.5     ~ "D1 Moderately dry",
    x >  -2.0     ~ "D2 Severely dry",
    TRUE          ~ "D3 Extremely dry")
}
CLASS_LEVELS <- c("D3 Extremely dry", "D2 Severely dry", "D1 Moderately dry", "N0 Near normal",
                  "W1 Moderately wet", "W2 Severely wet", "W3 Extremely wet")
CLASS_PAL <- c("D3 Extremely dry" = "#a50026", "D2 Severely dry" = "#d73027",
               "D1 Moderately dry" = "#fdae61", "N0 Near normal" = "grey70",
               "W1 Moderately wet" = "#abd9e9", "W2 Severely wet" = "#74add1",
               "W3 Extremely wet" = "#313695")

## 3o. Full monthly calendar ----------------------------------------------------
complete_monthly <- function(df, from = min(df$date), to = max(df$date)) {
  cal <- data.frame(date = seq(month_start(from), month_start(to), by = "month"))
  dplyr::left_join(cal, df, by = "date")
}

## 3p. Framework drought events + primary driver -------------------------------
msdi_event_table <- function(dates, values, season, contrib_mat, provisional) {
  ev <- identify_drought_events_ts(dates, values, scale = MSDI_EVENT_SCALE)
  in_ev <- attr(ev, "in_event")
  if (is.null(in_ev) || length(in_ev) != length(values)) in_ev <- rep(FALSE, length(values))
  in_ev[is.na(in_ev)] <- FALSE
  r <- rle(in_ev); ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1L
  runs <- which(r$values)
  if (length(runs) != nrow(ev))
    warning(sprintf("MSDI: %d event runs vs %d events from identify_drought_events_ts().",
                    length(runs), nrow(ev)), call. = FALSE)
  tab <- dplyr::bind_rows(lapply(seq_along(runs), function(k) {
    idx <- starts[runs[k]]:ends[runs[k]]; v <- values[idx]
    cm  <- suppressWarnings(colMeans(pmax(contrib_mat[idx, , drop = FALSE], 0), na.rm = TRUE))
    cm[!is.finite(cm)] <- NA_real_
    drv <- if (all(is.na(cm))) NA_character_ else sub("^contrib_", "", names(which.max(cm)))
    inten <- mean(DROUGHT_ONSET - v, na.rm = TRUE)
    data.frame(index = "MSDI", event_id = k,
               start_date = dates[min(idx)], end_date = dates[max(idx)],
               duration_months = length(idx),
               msdi_min = min(v, na.rm = TRUE), msdi_mean = mean(v, na.rm = TRUE),
               intensity = inten, severity = inten * length(idx),
               peak_class = msdi_class(min(v, na.rm = TRUE)),
               onset_season = season[min(idx)],
               seasons_spanned = paste(unique(stats::na.omit(season[idx])), collapse = " + "),
               primary_driver = drv,
               left_censored = min(idx) == 1L,
               provisional_months = sum(provisional[idx], na.rm = TRUE),
               stringsAsFactors = FALSE)
  }))
  list(events = tab, in_event = in_ev, min_duration = attr(ev, "min_duration"))
}

## 3q. Framework-style time-series PNG (same as Scripts 7 and 9) -------------
#   threshold_line : optional event threshold drawn as a red dashed line
#   flag_dates     : optional dates marked with open circles (e.g. provisional input)
plot_index_timeseries <- function(dates_vec, values_vec, index_label, title_label, file,
                                  subtitle_extra = NULL, event_scale = MSDI_EVENT_SCALE,
                                  threshold_line = NULL, flag_dates = NULL) {
  df <- data.frame(date = as.Date(dates_vec), value = values_vec)
  df <- df[order(df$date), ]
  if (!any(is.finite(df$value))) return(invisible(NULL))
  ev       <- identify_drought_events_ts(df$date, df$value, scale = event_scale)
  n_events <- nrow(ev); min_dur <- attr(ev, "min_duration")
  df$below <- attr(ev, "in_event")
  df <- df[is.finite(df$value), ]
  month_id <- as.integer(format(df$date, "%Y")) * 12L + as.integer(format(df$date, "%m"))
  df$grp   <- cumsum(c(1L, diff(as.integer(df$below)) != 0L | diff(month_id) != 1L))
  fl <- if (length(flag_dates)) df[df$date %in% as.Date(flag_dates), , drop = FALSE] else df[0, ]
  y_lo <- min(-2.2, floor(min(df$value) * 10) / 10)
  y_hi <- max( 2.2, ceiling(max(df$value) * 10) / 10)
  bnd <- function(ymin, ymax, fill, alpha)
    ggplot2::annotate("rect", xmin = DATE_NEG_INF, xmax = DATE_POS_INF,
                      ymin = ymin, ymax = ymax, fill = fill, alpha = alpha)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = date, y = value)) +
    bnd( 1, Inf, "#4472C4", 0.10) + bnd(0, 1, "#70AD47", 0.08) +
    bnd(-1, 0,   "#70AD47", 0.05) + bnd(-Inf, -1, "#C00000", 0.08) +
    ggplot2::geom_hline(yintercept = 0, color = "black", linewidth = 0.6) +
    ggplot2::geom_hline(yintercept = c(-1, 1), color = "grey30", linetype = "dashed", linewidth = 0.4) +
    ggplot2::geom_hline(yintercept = c(-1.5, 1.5, -2, 2), color = "grey55", linetype = "dotted", linewidth = 0.4)
  if (!is.null(threshold_line))
    p <- p + ggplot2::geom_hline(yintercept = threshold_line, color = "#C0504D",
                                 linetype = "longdash", linewidth = 0.55)
  if (any(df$below))
    p <- p + ggplot2::geom_ribbon(data = subset(df, below),
                                  ggplot2::aes(ymin = value, ymax = 0, group = grp),
                                  fill = "#C0504D", alpha = 0.45)
  p <- p + ggplot2::geom_line(color = "#1F77B4", linewidth = 0.45)
  if (nrow(fl))
    p <- p + ggplot2::geom_point(data = fl, shape = 21, fill = "white", color = "black",
                                 size = 1.8, stroke = 0.6)
  cap <- c(if (!is.null(threshold_line)) sprintf("Red dashed: event threshold (%.1f)", threshold_line),
           "dashed/dotted: McKee class bounds",
           if (nrow(fl)) "open circles: provisional input")
  p <- p +
    ggplot2::scale_x_date(date_breaks = "5 years", date_labels = "%Y",
                          expand = ggplot2::expansion(mult = c(0.01, 0.015))) +
    ggplot2::scale_y_continuous(breaks = seq(-5, 5, 1)) +
    ggplot2::coord_cartesian(ylim = c(y_lo, y_hi)) +
    ggplot2::labs(
      title    = sprintf("%s Time Series  |  NECHAKO BASIN", title_label),
      subtitle = sprintf("%d drought events  (onset < %.2f, end >= %.2f, min %d months)%s",
                         n_events, DROUGHT_ONSET, DROUGHT_END, min_dur,
                         if (!is.null(subtitle_extra)) paste0("  |  ", subtitle_extra) else ""),
      caption  = paste(cap, collapse = "; "),
      x = "Year", y = sprintf("%s Index Value", index_label)) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_line(color = "grey85"),
                   panel.grid.major.y = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(face = "bold", size = 15),
                   plot.subtitle = ggplot2::element_text(color = "grey40", size = 11),
                   plot.caption = ggplot2::element_text(color = "grey40", size = 9),
                   axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
                   plot.margin = ggplot2::margin(10, 16, 10, 10))
  ggplot2::ggsave(file, p, width = 11, height = 5.8, dpi = 200, bg = "white")
  cat(sprintf("  [PLOT] %s saved: %s  (%d events)\n", index_label, basename(file), n_events))
  invisible(p)
}

## 3r. Muffle only named cosmetic warnings ------------------------------------
muffle_cosmetic <- function(expr, patterns) {
  withCallingHandlers(expr, warning = function(w) {
    if (any(vapply(patterns, grepl, logical(1), x = conditionMessage(w))))
      invokeRestart("muffleWarning")
  })
}

# ==============================================================================
# 4.  LOAD & MERGE INDICES
# ==============================================================================
cat("\n", strrep("=", 60), "\n  LOADING INDEX TIME SERIES\n", strrep("=", 60), "\n\n", sep = "")

supply_long <- load_node("Supply", "SPEI-PM-3", "SPI/SPEI script (basin mean)",
                         file = cfg$supply_file, cols = cfg$supply_cols, panel_col = "spei_03")

snow_label <- if (cfg$snow_index == "sspi") "SSPI-1" else "SWEI-3"
snow_long  <- load_node("Snow", snow_label, "4snow_drought_indices_ERALand.R",
                        file = if (cfg$snow_index == "sspi") cfg$sspi_file else cfg$swei_file,
                        panel_col = if (cfg$snow_index == "sspi") "sspi_01" else "swei_03")

soil_label <- sprintf("%s-%d %s", if (cfg$soil_index == "ssmi") "SSMI" else "SSI-SM",
                      cfg$soil_scale, cfg$soil_layer)
soil_panel <- if (cfg$soil_index == "ssi_sm" && cfg$soil_layer == "L1_3" && cfg$soil_scale == 1L)
  "ssism_L13_01" else NULL
soil_long  <- load_node("Soil", soil_label, "8SSI_SSMI_ERALand.R",
                        file = cfg$soil_file, panel_col = soil_panel)

flow_label <- sprintf("SSI %s naturalized (B)", cfg$flow_station)
flow_long  <- load_node("Flow", flow_label, "6StreamFlow_Index_WSCLocal.R (Pipeline B)",
                        rds = cfg$flow_rds, panel_col = paste0("ssi_", cfg$flow_station))

sri_long <- if (cfg$run_sri_variant)
  load_node("SRI1", "SRI-1 (ERA5-Land, variant)", "SRI script (Script 5)",
            file = cfg$sri_file, cols = cfg$sri_cols, panel_col = "SRI_1") else NULL

# Record the inputs actually used in this analysis
manifest_df <- dplyr::bind_rows(manifest)

readr::write_csv(
  manifest_df,
  file.path(cfg$out_dir, "msdi_input_manifest.csv"),
  na = ""
)

cal <- data.frame(date = seq(cfg$analysis_start, cfg$analysis_end, by = "month"))
dat <- Reduce(function(a, b) dplyr::left_join(a, b, by = "date"),
              c(list(cal), purrr::compact(list(supply_long, snow_long, soil_long, flow_long, sri_long))))
if (!"flow_status" %in% names(dat)) dat$flow_status <- NA_character_
dat$year  <- as.integer(format(dat$date, "%Y"))
dat$month <- as.integer(format(dat$date, "%m"))
NODE_COLS <- c(ALL_NODES, if (cfg$run_sri_variant) "SRI1")

n_flow_prov <- sum(dat$flow_status %in% cfg$provisional_status)
cat(sprintf("\n  Calendar: %d months (%s - %s); reference/fit %s; Flow months with provisional input: %d\n",
            nrow(dat), format(min(dat$date), "%Y-%m"), format(max(dat$date), "%Y-%m"),
            REF_LABEL, n_flow_prov))
if (sum(is.finite(dat$Flow)) < nrow(dat))
  warning(sprintf("Flow node missing in %d of %d calendar months.",
                  nrow(dat) - sum(is.finite(dat$Flow)), nrow(dat)), call. = FALSE)

# ==============================================================================
# 5.  GRINGORTEN PSEUDO-OBSERVATIONS (per calendar month, reference period)
# ==============================================================================
cat("\n", strrep("=", 60), "\n  MARGINAL PSEUDO-OBSERVATIONS (reference ", REF_LABEL, ")\n",
    strrep("=", 60), "\n\n", sep = "")
for (v in NODE_COLS) {
  refs <- calibrate_marginals(dat, v, cfg$fit_start, cfg$fit_end)
  dat[[paste0("u_", v)]] <- apply_marginals(dat, v, refs)
  n_u <- sum(is.finite(dat[[paste0("u_", v)]]))
  cat(sprintf("  %-7s %d uniform values (%.1f%% of calendar)\n", v, n_u, 100 * n_u / nrow(dat)))
}

# ==============================================================================
# 6.  RELATIONSHIP & STATIONARITY DIAGNOSTICS PER REGIME
# ==============================================================================
cat("\n", strrep("=", 60), "\n  RELATIONSHIP DIAGNOSTICS BY REGIME (", REF_LABEL, ")\n",
    strrep("=", 60), "\n", sep = "")
in_fit <- dat$date >= cfg$fit_start & dat$date <= cfg$fit_end
td_rows <- list(); st_rows <- list()
for (rg in names(cfg$regimes)) {
  RG  <- cfg$regimes[[rg]]; ord <- RG$order
  sub <- dat[in_fit & dat$month %in% RG$months, c("year", paste0("u_", ord)), drop = FALSE]
  sub <- sub[stats::complete.cases(sub), , drop = FALSE]
  tm  <- stats::cor(sub[, paste0("u_", ord)], method = "kendall")
  dimnames(tm) <- list(ord, ord)
  cat(sprintf("\n  %s  (n = %d months, %d years)\n", RG$label, nrow(sub), length(unique(sub$year))))
  print(round(tm, 3))
  prs <- utils::combn(ord, 2)
  for (k in seq_len(ncol(prs))) {
    a <- paste0("u_", prs[1, k]); b <- paste0("u_", prs[2, k])
    td <- tail_dep(sub[[a]], sub[[b]])
    td_rows[[length(td_rows) + 1L]] <- data.frame(
      regime = RG$label, Var1 = prs[1, k], Var2 = prs[2, k],
      tree1_design = abs(match(prs[1, k], ord) - match(prs[2, k], ord)) == 1L,
      tau = stats::cor(sub[[a]], sub[[b]], method = "kendall"),
      lambda_L = td[["lambda_L"]], lambda_U = td[["lambda_U"]],
      independence_level = mean(cfg$tail_q), n = nrow(sub))
  }
  for (k in seq_len(length(ord) - 1L)) {
    a <- paste0("u_", ord[k]); b <- paste0("u_", ord[k + 1L])
    for (ct in cfg$stationarity_cuts) {
      set.seed(cfg$seed + ct + k)
      st <- split_test(sub[[a]], sub[[b]], sub$year, ct, cfg$n_perm)
      st_rows[[length(st_rows) + 1L]] <- data.frame(
        regime = RG$label, pair = paste(ord[k], ord[k + 1L], sep = "-"), split_year = ct,
        tau_early = st[["tau_early"]], tau_late = st[["tau_late"]],
        diff_early_minus_late = st[["diff"]], p_perm = st[["p_perm"]])
    }
  }
}
td_df <- dplyr::bind_rows(td_rows)
st_df <- dplyr::bind_rows(st_rows)
readr::write_csv(td_df, file.path(cfg$out_dir, "relationship_by_regime.csv"))
readr::write_csv(st_df, file.path(cfg$out_dir, "stationarity_tests.csv"))
cat("\n  Stationarity of tree-1 design pairs (year-block permutation):\n")
print(st_df |> dplyr::mutate(dplyr::across(dplyr::where(is.numeric), \(x) round(x, 3))), row.names = FALSE)
if (any(st_df$p_perm < 0.05, na.rm = TRUE))
  warning("At least one tree-1 pair shows a significant split-sample change in dependence - see stationarity_tests.csv.",
          call. = FALSE)

# ==============================================================================
# Sections 7-13 run only when RUN$msdi is TRUE. The block is left unindented
# so the code matches 10MSDI.r line for line.
# ==============================================================================
if (RUN$msdi) {
  
  # ==============================================================================
  # 7.  PER-MONTH VINES, KENDALL TRANSFORM, CONTRIBUTIONS
  # ==============================================================================
  fit_variant <- function(flow_col, tag, eval_end) {
    cat("\n", strrep("=", 60), sprintf("\n  VINE FITTING - %s (Flow slot = %s; BIC mode = %s)\n",
                                       tag, flow_col, cfg$bic_mode),
        strrep("=", 60), "\n", sep = "")
    cache <- list(); idx_rows <- list(); mod_rows <- list(); edge_rows <- list(); keep <- list()
    for (m in 1:12) {
      rg    <- regime_of(m); ord <- cfg$regimes[[rg]]$order
      dcols <- ifelse(ord == "Flow", flow_col, ord)
      ucols <- paste0("u_", dcols)
      w     <- fit_window(m); key <- paste(w, collapse = "-")
      
      if (is.null(cache[[key]])) {
        fd <- dat[in_fit & dat$month %in% w, c("date", "year", ucols), drop = FALSE]
        fd <- fd[stats::complete.cases(fd), , drop = FALSE]
        if (nrow(fd) < 40L) warning(sprintf("[%s] window %s: only %d complete months.", tag, key, nrow(fd)), call. = FALSE)
        u  <- as.matrix(fd[, ucols, drop = FALSE]); colnames(u) <- ord
        ny <- length(unique(fd$year))
        if (length(ord) < 3L) {
          # Bivariate: only one structure exists
          f0 <- fit_one(u, rvinecopulib::dvine_structure(seq_along(ord)))
          cand <- list(D_design = f0, D_greedy = f0, R_vine = f0)
          greedy <- seq_along(ord)
        } else {
          greedy <- order_dvine_by_tau(u)
          cand <- list(D_design = fit_one(u, rvinecopulib::dvine_structure(seq_along(ord))),
                       D_greedy = fit_one(u, rvinecopulib::dvine_structure(greedy)),
                       R_vine   = fit_one(u))
        }
        bic  <- vapply(cand, bic_years, numeric(1), n_rows = nrow(fd), n_years = ny)
        alt  <- names(bic)[-1][which.min(bic[-1])]
        pick <- if (bic[[alt]] < bic[["D_design"]] - cfg$bic_margin) alt else "D_design"
        cache[[key]] <- list(fit = cand[[pick]], u = u, year = fd$year,
                             hyear = fd$year + (as.integer(format(fd$date, "%m")) >= 11L),
                             pick = pick, bic = bic,
                             greedy = paste(ord[greedy], collapse = "-"), n = nrow(fd), ny = ny)
        cat(sprintf("  window %-8s %-24s n=%3d (%2d yr)  BIC design=%8.1f greedy=%8.1f R=%8.1f -> %s\n",
                    key, cfg$regimes[[rg]]$label, nrow(fd), ny,
                    bic[["D_design"]], bic[["D_greedy"]], bic[["R_vine"]], pick))
        e <- edge_table(cand[[pick]], ord)
        if (!is.null(e)) edge_rows[[key]] <- cbind(variant = tag, window = key, regime = rg, e)
      }
      cc <- cache[[key]]
      
      ed <- dat[dat$month == m & dat$date >= cfg$analysis_start & dat$date <= eval_end,
                c("date", ucols), drop = FALSE]
      ed <- ed[stats::complete.cases(ed), , drop = FALSE]
      if (!nrow(ed)) next
      U <- as.matrix(ed[, ucols, drop = FALSE]); colnames(U) <- ord
      
      eng <- kendall_engine(cc$fit, cfg$seed + 100L * m, cfg$n_ref, cfg$n_ref_tail, cfg$n_eval)
      res <- msdi_from_u(U, eng)
      n0     <- sum(eng$t_sim <= 0)
      floorv <- resolution_floor(eng)
      
      # Contribution: change in the vine CDF when a node is set to its median;
      # positive = the node was drier than normal. Normalised by the row |sum|.
      contrib <- matrix(NA_real_, nrow(U), length(ALL_NODES),
                        dimnames = list(NULL, paste0("contrib_", ALL_NODES)))
      for (j in seq_along(ord)) {
        Uj <- U; Uj[, j] <- 0.5
        contrib[, paste0("contrib_", ord[j])] <- t_hat(Uj, eng) - res$t
      }
      rs <- rowSums(abs(contrib), na.rm = TRUE)
      contrib <- contrib / ifelse(rs == 0, 1, rs)
      rm(eng); invisible(gc(verbose = FALSE))
      
      idx_rows[[m]] <- data.frame(date = ed$date, MSDI = res$msdi, mc_floor = res$floor, C_u = res$t,
                                  regime = rg, season = cfg$regimes[[rg]]$label,
                                  contrib, check.names = FALSE)
      mod_rows[[m]] <- data.frame(variant = tag, month = m, month_abb = month.abb[m], regime = rg,
                                  window = key, n_fit = cc$n, n_years = cc$ny, bic_mode = cfg$bic_mode,
                                  BIC_D_design = cc$bic[["D_design"]],
                                  BIC_D_greedy = cc$bic[["D_greedy"]],
                                  BIC_R_vine   = cc$bic[["R_vine"]],
                                  selected = cc$pick, design_order = paste(ord, collapse = "-"),
                                  greedy_order = cc$greedy, n_indep_edges = n_indep(cc$fit),
                                  n_eval = nrow(ed), n_mc_floor = sum(res$floor),
                                  n_tsim_zero = n0, msdi_resolution_floor = floorv,
                                  msdi_min_month = min(res$msdi))
      keep[[m]] <- list(key = key, U = U, dates = ed$date)
      cat(sprintf("    %s: %d months evaluated, %d at MC resolution limit (floor MSDI %.2f, min %.2f)\n",
                  month.abb[m], nrow(ed), sum(res$floor), floorv, min(res$msdi)))
    }
    list(index  = dplyr::bind_rows(idx_rows) |> dplyr::arrange(date),
         models = dplyr::bind_rows(mod_rows),
         edges  = dplyr::bind_rows(edge_rows),
         cache  = cache, keep = keep)
  }
  
  main    <- fit_variant("Flow", "main", cfg$analysis_end)
  sri_var <- if (cfg$run_sri_variant) fit_variant("SRI1", "sri1_variant", cfg$analysis_end) else NULL
  saveRDS(list(main = main$cache, sri1_variant = if (!is.null(sri_var)) sri_var$cache),
          file.path(cfg$out_dir, "vines_by_window.rds"))
  
  # ==============================================================================
  # 8.  ASSEMBLE MAIN SERIES; SRI-1 VARIANT AND AGREEMENT
  # ==============================================================================
  seasonal_out <- main$index |>
    dplyr::left_join(dat[, c("date", "flow_status")], by = "date") |>
    dplyr::mutate(provisional = flow_status %in% cfg$provisional_status | date %in% cfg$provisional_manual,
                  flow_source = flow_label) |>
    dplyr::select(-flow_status) |>
    dplyr::rename(MSDI_seasonal = MSDI) |>
    dplyr::arrange(date)
  contrib_cols <- paste0("contrib_", ALL_NODES)
  if (any(seasonal_out$provisional))
    warning(sprintf("%d MSDI month(s) rest on provisional/persisted flow input (flagged provisional = TRUE).",
                    sum(seasonal_out$provisional)), call. = FALSE)
  
  sri_agreement <- NULL
  if (!is.null(sri_var)) {
    ag <- dplyr::inner_join(main$index[, c("date", "regime", "MSDI")],
                            sri_var$index[, c("date", "MSDI")],
                            by = "date", suffix = c("_main", "_sri1"))
    agree <- function(g, label) {
      # Use identical finite pairs for all agreement metrics and guard empty D1 sets.
      g <- g[is.finite(g$MSDI_main) & is.finite(g$MSDI_sri1), , drop = FALSE]
      main_d1 <- g$MSDI_main <= -1
      sri_d1  <- g$MSDI_sri1 <= -1
      n_main_d1 <- sum(main_d1)
      n_sri_d1  <- sum(sri_d1)
      data.frame(
        regime = label, n_months = nrow(g),
        kendall_tau = safe_cor(g$MSDI_main, g$MSDI_sri1, method = "kendall"),
        pearson_r   = safe_cor(g$MSDI_main, g$MSDI_sri1),
        rmse = if (nrow(g)) sqrt(mean((g$MSDI_main - g$MSDI_sri1)^2)) else NA_real_,
        hit_rate_D1 = if (n_main_d1) sum(main_d1 & sri_d1) / n_main_d1 else NA_real_,
        false_alarm_D1 = if (n_sri_d1) sum(sri_d1 & !main_d1) / n_sri_d1 else NA_real_,
        n_main_D1 = n_main_d1, n_sri_D1 = n_sri_d1)
    }
    sri_agreement <- do.call(rbind, c(
      list(agree(ag, "all")),
      unname(lapply(split(ag, ag$regime), function(g) agree(g, g$regime[1])))))
    readr::write_csv(sri_agreement, file.path(cfg$out_dir, "sri1_variant_agreement.csv"))
    readr::write_csv(dplyr::rename(sri_var$index, MSDI_sri1 = MSDI, contrib_SRI1 = contrib_Flow),
                     file.path(cfg$out_dir, "msdi_sri1_variant.csv"), na = "NA")
    cat("\n  SRI-1 variant vs main MSDI-N (both fitted on ", REF_LABEL, "; in-sample):\n", sep = "")
    print(sri_agreement |> dplyr::mutate(dplyr::across(dplyr::where(is.numeric), \(x) round(x, 3))),
          row.names = FALSE)
  }
  
  # ==============================================================================
  # 9.  DROUGHT EVENTS & MAIN TABLE
  # ==============================================================================
  cat("\n", strrep("=", 60), "\n  DROUGHT EVENT IDENTIFICATION (framework)\n", strrep("=", 60), "\n\n", sep = "")
  msdi_cal <- complete_monthly(seasonal_out, from = cfg$analysis_start, to = cfg$analysis_end) |>
    dplyr::left_join(dat[, c("date", NODE_COLS)], by = "date")
  msdi_cal$provisional[is.na(msdi_cal$provisional)] <- FALSE
  
  ev_res <- msdi_event_table(msdi_cal$date, msdi_cal$MSDI_seasonal, msdi_cal$season,
                             as.matrix(msdi_cal[, contrib_cols]), msdi_cal$provisional)
  drought_events <- ev_res$events
  cat(sprintf("  Events (onset < %.2f, end >= %.2f, min %d months): %d\n",
              DROUGHT_ONSET, DROUGHT_END, ev_res$min_duration, nrow(drought_events)))
  if (nrow(drought_events)) {
    ev_floor <- vapply(seq_len(nrow(drought_events)), function(k) {
      s <- msdi_cal$date >= drought_events$start_date[k] & msdi_cal$date <= drought_events$end_date[k]
      sum(msdi_cal$mc_floor[s], na.rm = TRUE)
    }, numeric(1))
    drought_events$months_at_mc_floor <- ev_floor
    cat("  Top events by severity:\n")
    print(drought_events |> dplyr::arrange(dplyr::desc(severity)) |> utils::head(10) |>
            dplyr::select(start_date, end_date, duration_months, msdi_min, severity,
                          peak_class, primary_driver, left_censored, months_at_mc_floor) |>
            dplyr::mutate(dplyr::across(dplyr::where(is.numeric), \(x) round(x, 2))))
  }
  readr::write_csv(drought_events, file.path(cfg$out_dir, "msdi_drought_events.csv"), na = "")
  
  msdi_out <- msdi_cal |>
    dplyr::mutate(drought_cat = msdi_class(MSDI_seasonal), in_drought_event = ev_res$in_event) |>
    dplyr::select(date, MSDI_seasonal, season, regime, drought_cat, in_drought_event,
                  provisional, flow_source, mc_floor, C_u,
                  dplyr::all_of(contrib_cols), dplyr::all_of(NODE_COLS))
  readr::write_csv(msdi_out, file.path(cfg$out_dir, "Nechako_MSDI_VineCopula.csv"), na = "NA")
  
  wide <- data.frame(Year = sort(unique(as.integer(format(msdi_out$date, "%Y")))))
  for (m in 1:12) {
    s <- msdi_out[as.integer(format(msdi_out$date, "%m")) == m, ]
    wide[[month.abb[m]]] <- s$MSDI_seasonal[match(wide$Year, as.integer(format(s$date, "%Y")))]
  }
  utils::write.csv(wide, file.path(cfg$out_dir, "msdi_01_by_month.csv"), row.names = FALSE, na = "")
  cat("\n  Main CSV written: Nechako_MSDI_VineCopula.csv\n")
  
  # ==============================================================================
  # 10. DIAGNOSTICS
  # ==============================================================================
  cat("\n", strrep("=", 60), "\n  DIAGNOSTICS\n", strrep("=", 60), "\n\n", sep = "")
  ref_out <- seasonal_out |>
    dplyr::filter(date >= cfg$fit_start, date <= cfg$fit_end, !provisional, is.finite(MSDI_seasonal))
  norm_tab <- ref_out |>
    dplyr::group_by(regime) |>
    dplyr::summarise(n = dplyr::n(), mean = mean(MSDI_seasonal), sd = stats::sd(MSDI_seasonal),
                     ks_p = suppressWarnings(stats::ks.test(MSDI_seasonal, "pnorm")$p.value),
                     n_mc_floor = sum(mc_floor), .groups = "drop")
  norm_tab <- dplyr::bind_rows(norm_tab, tibble::tibble(
    regime = "all", n = nrow(ref_out), mean = mean(ref_out$MSDI_seasonal), sd = stats::sd(ref_out$MSDI_seasonal),
    ks_p = suppressWarnings(stats::ks.test(ref_out$MSDI_seasonal, "pnorm")$p.value),
    n_mc_floor = sum(ref_out$mc_floor)))
  print(norm_tab |> dplyr::mutate(dplyr::across(dplyr::where(is.numeric), \(x) round(x, 4))))
  readr::write_csv(norm_tab, file.path(cfg$out_dir, "msdi_normality_by_regime.csv"))
  
  model_cmp <- dplyr::bind_rows(main$models, if (!is.null(sri_var)) sri_var$models)
  readr::write_csv(model_cmp, file.path(cfg$out_dir, "model_comparison.csv"))
  cat("\n  Structure selected per month (main):\n")
  print(main$models[, c("month_abb", "regime", "window", "selected", "design_order",
                        "greedy_order", "n_indep_edges", "msdi_resolution_floor")], row.names = FALSE)
  edges_df <- dplyr::bind_rows(main$edges, if (!is.null(sri_var)) sri_var$edges)
  readr::write_csv(edges_df, file.path(cfg$out_dir, "vine_edges.csv"))
  cat("\n  Expected checks: accumulation tree-3 Snow-Soil edge tau < 0;",
      "melt tree-2 Supply-Flow | Soil edge 'indep'. See vine_edges.csv.\n")
  
  # ==============================================================================
  # 11. OPTIONAL BOOTSTRAP (Nov-Oct hydrological-year blocks, same structure as reported)
  # ==============================================================================
  if (cfg$run_bootstrap) {
    cat("\n", strrep("=", 60), "\n  BOOTSTRAP (n = ", cfg$n_boot, ")\n", strrep("=", 60), "\n\n", sep = "")
    boot_out <- list()
    for (m in 1:12) {
      kp <- main$keep[[m]]; if (is.null(kp)) next
      cc <- main$cache[[kp$key]]
      blocks <- split(seq_len(nrow(cc$u)), cc$hyear)
      mat <- matrix(NA_real_, nrow(kp$U), cfg$n_boot)
      for (b in seq_len(cfg$n_boot)) {
        set.seed(cfg$seed + 1000L * m + b)
        ub <- cc$u[unlist(blocks[sample(length(blocks), replace = TRUE)]), , drop = FALSE]
        vb <- tryCatch(rvinecopulib::vinecop(ub, structure = cc$fit$structure, family_set = cfg$family_set,
                                             par_method = cfg$par_method, selcrit = cfg$selcrit),
                       error = function(e) NULL)
        if (is.null(vb)) next
        eng <- kendall_engine(vb, cfg$seed + b, cfg$boot_n_ref, cfg$boot_n_ref_tail, cfg$boot_n_eval)
        mat[, b] <- msdi_from_u(kp$U, eng)$msdi
      }
      boot_out[[m]] <- data.frame(date = kp$dates,
                                  MSDI_lo = apply(mat, 1, stats::quantile, 0.05, na.rm = TRUE),
                                  MSDI_hi = apply(mat, 1, stats::quantile, 0.95, na.rm = TRUE))
      cat(sprintf("  %s done\n", month.abb[m]))
    }
    readr::write_csv(dplyr::bind_rows(boot_out) |> dplyr::arrange(date),
                     file.path(cfg$out_dir, "MSDI_bootstrap_CI.csv"))
  }
  
  # ==============================================================================
  # 12. PLOTS
  # ==============================================================================
  cat("\n", strrep("=", 60), "\n  GENERATING PLOTS\n", strrep("=", 60), "\n\n", sep = "")
  NODE_PAL <- c(Flow = "#E41A1C", Snow = "#377EB8", Soil = "#4DAF4A", Supply = "#984EA3")
  plot_theme <- ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                   plot.caption = ggplot2::element_text(colour = "grey40", size = 8),
                   panel.grid.minor = ggplot2::element_blank())
  dec_of <- function(d) paste0(10L * (lubridate::year(d) %/% 10L), "s")
  
  ## 12a. Time series ------------------------------------------------------------
  plot_index_timeseries(
    msdi_cal$date, msdi_cal$MSDI_seasonal, "MSDI-N", "MSDI-N (Vine Copula, Regime-based)",
    file.path(FIG_DIR, "msdi_basin_timeseries.png"),
    subtitle_extra = sprintf("Nov-Apr 4 | May-Jun 4 | Jul-Oct 2 nodes; reference %s", REF_LABEL),
    threshold_line = DROUGHT_ONSET,
    flag_dates     = msdi_cal$date[msdi_cal$provisional])
  
  ## 12b. Heat map ---------------------------------------------------------------
  hm_data <- msdi_out |>
    dplyr::filter(is.finite(MSDI_seasonal)) |>
    dplyr::mutate(Year  = factor(lubridate::year(date), levels = rev(sort(unique(lubridate::year(date))))),
                  Month = lubridate::month(date, label = TRUE, abbr = TRUE),
                  Class = factor(drought_cat, levels = CLASS_LEVELS))
  HM_PAL <- CLASS_PAL; HM_PAL["N0 Near normal"] <- "grey92"
  fill_layers <- if (isTRUE(cfg$heatmap_binned)) list(
    ggplot2::geom_tile(ggplot2::aes(fill = Class), colour = "white", linewidth = 0.3),
    ggplot2::scale_fill_manual(values = HM_PAL, name = "MSDI-N class", drop = FALSE)
  ) else list(
    ggplot2::geom_tile(ggplot2::aes(fill = MSDI_seasonal), colour = "white", linewidth = 0.3),
    ggplot2::scale_fill_gradientn(colours = c("#a50026", "#d73027", "#fdae61", "grey90", "#74add1", "#313695"),
                                  limits = c(-3, 3), oob = scales::squish, name = "MSDI-N"))
  p3 <- ggplot2::ggplot(hm_data, ggplot2::aes(Month, Year)) +
    fill_layers +
    ggplot2::geom_tile(data = hm_data[hm_data$provisional, ], fill = NA, colour = "black", linewidth = 0.7) +
    ggplot2::labs(title = "Nechako Basin - MSDI-N Heat Map",
                  subtitle = sprintf("Reference %s; naturalized SSI %s in the Flow node throughout",
                                     REF_LABEL, cfg$flow_station),
                  caption = if (any(hm_data$provisional)) "Outlined tiles: provisional flow input" else NULL,
                  x = NULL, y = NULL) +
    plot_theme + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 7))
  ggplot2::ggsave(file.path(cfg$out_dir, "plot3_MSDI_heatmap.png"), p3, width = 9,
                  height = max(4, dplyr::n_distinct(hm_data$Year) * 0.22 + 1.5), dpi = 300, bg = "white")
  
  ## 12c. Component contributions -----------------------------------------------
  # Complete month x node grid: nodes absent from a regime (Supply, Snow in Jul-Oct)
  # contribute 0, so every month has exactly one value per node and bars never interpolate.
  contrib_long <- msdi_out |>
    dplyr::select(date, dplyr::all_of(contrib_cols)) |>
    tidyr::pivot_longer(-date, names_to = "node", values_to = "contribution") |>
    dplyr::mutate(node  = factor(sub("^contrib_", "", node), levels = names(NODE_PAL)),
                  share = dplyr::if_else(is.finite(contribution), pmax(contribution, 0), 0),
                  decade = dec_of(date))
  stopifnot(all(table(contrib_long$date) == length(ALL_NODES)))
  stack_tot <- tapply(contrib_long$share, contrib_long$date, sum)
  if (any(stack_tot > 1 + 1e-9)) warning("Stacked positive contribution exceeds 1 - check normalisation.", call. = FALSE)
  cat(sprintf("  Contributions: max stacked positive share = %.3f (<= 1 by construction)\n", max(stack_tot)))
  
  # Full record, one panel per decade; every panel spans a full decade
  d0 <- as.integer(sub("s$", "", sort(unique(contrib_long$decade))))
  dec_span <- data.frame(decade = rep(paste0(d0, "s"), 2),
                         date   = c(as.Date(sprintf("%d-01-01", d0)), as.Date(sprintf("%d-01-01", d0 + 10L))))
  ev_shade <- msdi_out |>
    dplyr::filter(in_drought_event %in% TRUE) |>
    dplyr::transmute(xmin = date - 15, xmax = date + 16, decade = dec_of(date))
  p4 <- ggplot2::ggplot(contrib_long, ggplot2::aes(date, share, fill = node)) +
    ggplot2::geom_rect(data = ev_shade, ggplot2::aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 1),
                       inherit.aes = FALSE, fill = "grey88") +
    ggplot2::geom_col(width = 30) +
    ggplot2::geom_blank(data = dec_span, ggplot2::aes(x = date), inherit.aes = FALSE) +
    ggplot2::facet_wrap(~ decade, ncol = 1, scales = "free_x") +
    ggplot2::scale_fill_manual(values = NODE_PAL, name = "Node") +
    ggplot2::scale_x_date(date_breaks = "1 year", date_labels = "%Y", expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), expand = c(0, 0)) +
    ggplot2::labs(title = "Component Drought Contribution",
                  subtitle = "Positive share of the joint-CDF change when each node is set to its median; grey = MSDI-N event months",
                  caption = "Jul-Oct: only Soil and Flow are nodes. Gap below 1 = share of nodes wetter than their median.",
                  x = NULL, y = "Normalised positive contribution") +
    plot_theme + ggplot2::theme(legend.position = "bottom",
                                strip.text = ggplot2::element_text(face = "bold"))
  ggplot2::ggsave(file.path(cfg$out_dir, "plot4_component_contribution.png"), p4,
                  width = 11, height = 1.6 * length(d0) + 1.5, dpi = 300, bg = "white")
  
  # Close-ups of the most severe events
  ev_rank <- drought_events |> dplyr::arrange(dplyr::desc(severity))
  for (k in seq_len(min(cfg$plot_zoom_n, nrow(ev_rank)))) {
    e   <- ev_rank[k, ]; pad <- cfg$plot_zoom_pad + 1L
    x0  <- seq(e$start_date, by = "-1 month", length.out = pad)[pad]
    x1  <- seq(e$end_date,   by =  "1 month", length.out = pad)[pad]
    zd  <- dplyr::filter(contrib_long, date >= x0, date <= x1)
    zl  <- msdi_out |>
      dplyr::filter(date >= x0, date <= x1, is.finite(MSDI_seasonal)) |>
      dplyr::mutate(lab = paste0(sprintf("%.1f", MSDI_seasonal), ifelse(provisional, "*", "")))
    pz <- ggplot2::ggplot(zd, ggplot2::aes(date, share, fill = node)) +
      ggplot2::annotate("rect", xmin = e$start_date - 15, xmax = e$end_date + 15,
                        ymin = 0, ymax = 1, fill = "grey92") +
      ggplot2::geom_col(width = 27, colour = "white", linewidth = 0.2) +
      ggplot2::geom_text(data = zl, ggplot2::aes(date, 1.03, label = lab), inherit.aes = FALSE,
                         size = 2.6, vjust = 0) +
      ggplot2::scale_fill_manual(values = NODE_PAL, name = "Node", drop = FALSE) +
      ggplot2::scale_x_date(date_breaks = "2 months", date_labels = "%b\n%Y",
                            expand = ggplot2::expansion(add = 20)) +
      ggplot2::scale_y_continuous(limits = c(0, 1.12), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
      ggplot2::labs(title = sprintf("Node contributions - event rank %d (%s to %s, severity %.1f)", k,
                                    format(e$start_date, "%Y-%m"), format(e$end_date, "%Y-%m"), e$severity),
                    subtitle = "Positive share of the joint-CDF change when each node is set to its median; grey band = event; numbers = MSDI-N",
                    caption = paste("Gap below 1 = share of nodes wetter than their median. Jul-Oct: only Soil and Flow are nodes.",
                                    if (any(zl$provisional)) "* provisional input." else ""),
                    x = NULL, y = "Normalised positive contribution") +
      plot_theme + ggplot2::theme(legend.position = "bottom")
    ggplot2::ggsave(file.path(cfg$out_dir, sprintf("plot4_component_contribution_event%d.png", k)), pz,
                    width = 9, height = 4.8, dpi = 300, bg = "white")
  }
  
  ## 12d. Pairs plots (pooled tau only; class shown in points and histograms) ---
  if (requireNamespace("GGally", quietly = TRUE)) {
    upper_tau <- function(data, mapping, ...) {
      x <- GGally::eval_data_col(data, mapping$x); y <- GGally::eval_data_col(data, mapping$y)
      tau <- stats::cor(x, y, method = "kendall", use = "complete.obs")
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5, label = sprintf("\u03c4 = %.3f", tau), size = 5,
                          colour = if (tau < 0) "#b2182b" else "grey15") +
        ggplot2::scale_x_continuous(limits = c(0, 1)) + ggplot2::scale_y_continuous(limits = c(0, 1))
    }
    for (rg in intersect(cfg$pairs_regimes, names(cfg$regimes))) {
      RG <- cfg$regimes[[rg]]; ord <- RG$order
      if (length(ord) < 3L) next
      muffle_cosmetic({
        pd <- dat[in_fit & dat$month %in% RG$months, c("date", paste0("u_", ord))]
        pd <- pd[stats::complete.cases(pd), ]
        names(pd) <- c("date", ord)
        pd$Class <- factor(msdi_class(seasonal_out$MSDI_seasonal[match(pd$date, seasonal_out$date)]),
                           levels = CLASS_LEVELS)
        p5 <- GGally::ggpairs(pd, columns = 2:(1 + length(ord)),
                              mapping = ggplot2::aes(colour = Class, fill = Class),
                              upper = list(continuous = upper_tau),
                              lower = list(continuous = GGally::wrap("points", size = 0.7, alpha = 0.6)),
                              diag  = list(continuous = GGally::wrap("barDiag", bins = 10, colour = NA)),
                              legend = 1) +
          ggplot2::scale_colour_manual(values = CLASS_PAL, name = "MSDI-N class", drop = TRUE) +
          ggplot2::scale_fill_manual(values = CLASS_PAL, name = "MSDI-N class", drop = TRUE) +
          ggplot2::labs(title = sprintf("%s pseudo-observations, %s (n = %d months)", RG$label, REF_LABEL, nrow(pd)),
                        caption = "Upper: pooled Kendall tau. Lower: pseudo-observations by MSDI-N class. Diagonal: class composition per decile.") +
          ggplot2::theme_bw(base_size = 9) + ggplot2::theme(legend.position = "bottom")
        ggplot2::ggsave(file.path(cfg$out_dir, sprintf("plot5_pair_plot_%s.png",
                                                       c(acc = "accumulation", melt = "melt")[[rg]])),
                        p5, width = 10, height = 10, dpi = 250, bg = "white")
      }, patterns = c("No shared levels found", "exact p-value with ties", "in `summarise"))
    }
  } else cat("  Plot 5 skipped (GGally not installed)\n")
  
  ## 12e. Event scatter (selective labels) --------------------------------------
  if (nrow(drought_events)) {
    ev_plot <- drought_events |>
      dplyr::mutate(rank  = rank(-severity, ties.method = "first"),
                    start = format(start_date, "%Y-%m"),
                    label = ifelse(rank <= cfg$plot_label_top | start %in% cfg$plot_label_extra,
                                   paste0(start, ifelse(provisional_months > 0, "*", "")), ""),
                    peak_class = factor(peak_class, levels = CLASS_LEVELS))
    lab_layer <- if (requireNamespace("ggrepel", quietly = TRUE))
      ggrepel::geom_text_repel(size = 3, max.overlaps = Inf, min.segment.length = 0,
                               box.padding = 0.45, seed = cfg$seed, show.legend = FALSE) else
                                 ggplot2::geom_text(data = function(x) x[x$label != "", ], size = 3, vjust = -0.9, show.legend = FALSE)
    p6 <- ggplot2::ggplot(ev_plot, ggplot2::aes(duration_months, severity, colour = peak_class,
                                                size = abs(msdi_min), label = label)) +
      ggplot2::geom_point(alpha = 0.8) + lab_layer +
      ggplot2::scale_colour_manual(values = CLASS_PAL, name = "Peak class", drop = TRUE) +
      ggplot2::scale_size_continuous(range = c(2, 8), name = "|Peak MSDI-N|") +
      ggplot2::labs(title = "Nechako MSDI-N Drought Event Catalogue",
                    subtitle = sprintf("Framework events: onset < %.2f, end >= %.2f, min %d months; S = I x D",
                                       DROUGHT_ONSET, DROUGHT_END, ev_res$min_duration),
                    caption = sprintf("Labelled: %d most severe events%s; * event includes provisional months.",
                                      cfg$plot_label_top,
                                      if (length(cfg$plot_label_extra))
                                        paste0(" and ", paste(cfg$plot_label_extra, collapse = ", ")) else ""),
                    x = "Duration (months)", y = "Severity (I x D)") +
      plot_theme
    ggplot2::ggsave(file.path(cfg$out_dir, "plot6_drought_event_scatter.png"), p6,
                    width = 8, height = 6, dpi = 250, bg = "white")
  }
  
  ## 12f. Tail co-occurrence ------------------------------------------------------
  p7 <- td_df |>
    dplyr::mutate(pair = paste0(Var1, "-", Var2),
                  regime = factor(regime, levels = vapply(cfg$regimes, `[[`, character(1), "label"))) |>
    tidyr::pivot_longer(c(lambda_L, lambda_U), names_to = "type", values_to = "lambda") |>
    dplyr::mutate(type = dplyr::recode(type, lambda_L = "Lower tail (dry)", lambda_U = "Upper tail (wet)")) |>
    ggplot2::ggplot(ggplot2::aes(pair, lambda, fill = type)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.75), width = 0.7) +
    ggplot2::geom_hline(yintercept = mean(cfg$tail_q), linetype = "dashed", colour = "grey40") +
    ggplot2::facet_grid(. ~ regime, scales = "free_x", space = "free_x") +
    ggplot2::scale_fill_manual(values = c("Lower tail (dry)" = "#d73027", "Upper tail (wet)" = "#4575b4"), name = NULL) +
    ggplot2::labs(title = "Empirical tail co-occurrence by regime",
                  subtitle = sprintf("Mean over q = %s of P(V <= q | U <= q); dashed = independence; %s",
                                     paste(cfg$tail_q, collapse = ", "), REF_LABEL),
                  caption = "Accumulation Snow-Soil is negatively dependent: both concordant tails at or below independence, as expected under substitution.",
                  x = NULL, y = "Conditional co-occurrence") +
    plot_theme +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1), legend.position = "top")
  ggplot2::ggsave(file.path(cfg$out_dir, "plot7_tail_dependence.png"), p7,
                  width = 12, height = 4.8, dpi = 300, bg = "white")
  cat("  Plots written.\n")
  
  # ==============================================================================
  # 13. SUMMARY FILE
  # ==============================================================================
  z   <- msdi_out$MSDI_seasonal
  fin <- msdi_out[is.finite(z), ]
  sel_tab <- table(main$models$selected)
  sri_all <- if (!is.null(sri_agreement)) sri_agreement[sri_agreement$regime == "all", ] else NULL
  lines <- c(
    "============================================================",
    "MSDI-N SUMMARY - Nechako River Basin (vine copula, regime-based)",
    "============================================================",
    sprintf("Run date          : %s  (version %s)", format(Sys.time(), "%Y-%m-%d %H:%M"), VERSION_TAG),
    sprintf("Nodes             : Supply SPEI-PM-3 | Snow %s | Soil %s | Flow %s", snow_label, soil_label, flow_label),
    sprintf("Regimes           : %s", paste(vapply(cfg$regimes, function(r)
      sprintf("%s [%s]", r$label, paste(r$order, collapse = "-")), character(1)), collapse = " ; ")),
    sprintf("Fit / reference   : %s to %s; Gringorten per calendar month; +/-1-month window within regime",
            format(cfg$fit_start, "%Y-%m"), format(cfg$fit_end, "%Y-%m")),
    sprintf("Structure choice  : %s (BIC mode %s, margin %g)",
            paste(sprintf("%s=%d", names(sel_tab), as.integer(sel_tab)), collapse = ", "),
            cfg$bic_mode, cfg$bic_margin),
    sprintf("MC engine         : n_ref=%d, n_ref_tail=%d, n_eval=%d; months at resolution limit: %d; floor MSDI range %.2f to %.2f",
            as.integer(cfg$n_ref), as.integer(cfg$n_ref_tail), as.integer(cfg$n_eval),
            sum(msdi_out$mc_floor, na.rm = TRUE),
            min(main$models$msdi_resolution_floor), max(main$models$msdi_resolution_floor)),
    sprintf("MSDI months       : %d, %s to %s (%d flagged provisional flow input)", nrow(fin),
            format(min(fin$date), "%Y-%m"), format(max(fin$date), "%Y-%m"), sum(fin$provisional)),
    if (!is.null(sri_all))
      sprintf("SRI-1 variant     : tau=%.3f r=%.3f RMSE=%.3f hit(D1)=%.2f false-alarm(D1)=%.2f (all regimes, in-sample)",
              sri_all$kendall_tau, sri_all$pearson_r, sri_all$rmse,
              sri_all$hit_rate_D1, sri_all$false_alarm_D1) else "SRI-1 variant     : not run",
    sprintf("Stationarity      : min p_perm over %d tree-1 tests = %.3f",
            sum(is.finite(st_df$p_perm)), suppressWarnings(min(st_df$p_perm, na.rm = TRUE))),
    sprintf("Normality (KS p)  : %s", paste(sprintf("%s=%.3f", norm_tab$regime, norm_tab$ks_p), collapse = ", ")),
    sprintf("Drought events    : %d (identify_drought_events_ts, scale %d; onset < %.2f, end >= %.2f, min %d months)",
            nrow(drought_events), MSDI_EVENT_SCALE, DROUGHT_ONSET, DROUGHT_END, ev_res$min_duration),
    "",
    "Class frequency (McKee et al. 1993 bins):",
    sprintf("  W3 Extremely wet  (>= 2.0)    : %.1f%%", 100 * mean(z >= 2, na.rm = TRUE)),
    sprintf("  W2 Severely wet   [1.5, 2)    : %.1f%%", 100 * mean(z >= 1.5 & z < 2, na.rm = TRUE)),
    sprintf("  W1 Moderately wet [1, 1.5)    : %.1f%%", 100 * mean(z >= 1 & z < 1.5, na.rm = TRUE)),
    sprintf("  N0 Near normal    (-1, 1)     : %.1f%%", 100 * mean(z > -1 & z < 1, na.rm = TRUE)),
    sprintf("  D1 Moderately dry (-1.5, -1]  : %.1f%%", 100 * mean(z <= -1 & z > -1.5, na.rm = TRUE)),
    sprintf("  D2 Severely dry   (-2, -1.5]  : %.1f%%", 100 * mean(z <= -1.5 & z > -2, na.rm = TRUE)),
    sprintf("  D3 Extremely dry  (<= -2.0)   : %.1f%%", 100 * mean(z <= -2, na.rm = TRUE)),
    "",
    "Inputs:",
    sprintf("  %-7s %-30s %5d  %s-%s  %s", manifest_df$node, manifest_df$index, manifest_df$n_records,
            manifest_df$start, manifest_df$end, manifest_df$file)
  )
  writeLines(lines, file.path(cfg$out_dir, "msdi_summary.txt"))
  
  cat("\n", strrep("=", 60), "\n  VINE COPULA MSDI-N - COMPLETE\n", strrep("=", 60), "\n\n", sep = "")
  cat(sprintf("  MSDI months: %d (%s - %s), provisional flow input: %d\n", nrow(fin),
              format(min(fin$date), "%Y-%m"), format(max(fin$date), "%Y-%m"), sum(fin$provisional)))
  cat(sprintf("  Drought events: %d\n", nrow(drought_events)))
  cat(sprintf("  Outputs in: %s\n", cfg$out_dir))
  
}  # end if (RUN$msdi)

# ==============================================================================
# PART B - RELATIONSHIP ANALYSIS (formerly diags_10.R, A2_relationship_v5)
#   Nodes : Supply SPEI-PM-3 | Snow SSPI-1 | Soil SSI-SM-1 L1_3 | Flow SSI 08JC002
#   Period: cfg$fit_start .. cfg$fit_end, Gringorten margins from Section 5
#   Bootstrap: year blocks = Nov-Oct hydrological years (a winter stays in one block)
#   Outputs (<cfg$out_dir>/relationship_v5/):
#     A2_tau_by_month.csv, A2_unexplained_pcor.csv, A2_cascades.csv, A2_mediation.csv
# ==============================================================================
if (RUN$relationship) {
  cat("\n", strrep("=", 60), "\n  RELATIONSHIP ANALYSIS (A2 v5)\n", strrep("=", 60), "\n", sep = "")
  t_rel <- Sys.time()
  rcfg <- list(
    B       = 500L,   # bootstrap replicates, Tables 1 and 3
    B_med   = 200L,   # bootstrap replicates, Table 4 (copula refits)
    seed    = 1L,
    out_dir = file.path(cfg$out_dir, "relationship_v5"))
  dir.create(rcfg$out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # ---- B0: fit-period panel, hydrological year, normal scores ---------------
  rel <- dat[in_fit, c("date", "year", "month", paste0("u_", ALL_NODES))]
  rel <- rel[order(rel$date), ]
  rel$hyear <- rel$year + (rel$month >= 11L)        # Nov-Oct hydrological year: bootstrap block
  for (v in ALL_NODES) rel[[paste0("z_", v)]] <- stats::qnorm(rel[[paste0("u_", v)]])
  if (any(!is.finite(as.matrix(rel[, paste0("u_", ALL_NODES)]))))
    warning("Some pseudo-observations are missing; check the inputs.", call. = FALSE)
  
  # Consistency check against Section 6 (v5 values: 0.435 / 0.323 / 0.448)
  acc_rows <- rel[rel$month %in% cfg$regimes$acc$months, ]
  acc_rows <- acc_rows[stats::complete.cases(acc_rows[, paste0("u_", ALL_NODES)]), ]
  chk <- td_df[td_df$regime == cfg$regimes$acc$label, ]
  for (p in list(c("Snow", "Supply"), c("Supply", "Flow"), c("Flow", "Soil"))) {
    here <- stats::cor(acc_rows[[paste0("u_", p[1])]], acc_rows[[paste0("u_", p[2])]], method = "kendall")
    ref  <- chk$tau[(chk$Var1 == p[1] & chk$Var2 == p[2]) | (chk$Var1 == p[2] & chk$Var2 == p[1])]
    cat(sprintf("  Check %s-%s: %.3f (relationship_by_regime %.3f)\n", p[1], p[2], here, ref))
    if (length(ref) != 1L || abs(here - ref) > 1e-8)
      warning(sprintf("Part B tau for %s-%s differs from Section 6.", p[1], p[2]), call. = FALSE)
  }
  
  # The Part B target-month windows use the *same* regime boundary rule as
  # Part A's vine fits. In particular, Apr cannot include May; Jun cannot
  # include Jul; Oct cannot include Nov when window_within_regime = TRUE.
  win3 <- fit_window
  boot_ci <- function(s, fun, B = rcfg$B) {      # hydrological-year block bootstrap, 95% CI
    if (nrow(s) < 3L || length(unique(s$hyear)) < 2L)
      return(c(NA_real_, NA_real_))
    blocks <- split(seq_len(nrow(s)), s$hyear)
    bs <- replicate(B, tryCatch(fun(s[unlist(blocks[sample(length(blocks), replace = TRUE)]), , drop = FALSE]),
                                error = function(e) NA_real_))
    bs <- bs[is.finite(bs)]
    if (length(bs) < 2L) return(c(NA_real_, NA_real_))
    unname(stats::quantile(bs, c(0.025, 0.975), na.rm = TRUE))
  }
  set.seed(rcfg$seed)
  
  # ---- B1: Kendall tau per target month (regime-clipped window), with CI ----
  cat("\nPart B 1/4: tau by month\n")
  prs <- utils::combn(ALL_NODES, 2)
  t1 <- do.call(rbind, lapply(1:12, function(m) {
    window <- rel[rel$month %in% win3(m), , drop = FALSE]
    do.call(rbind, lapply(seq_len(ncol(prs)), function(k) {
      a <- paste0("u_", prs[1, k]); b <- paste0("u_", prs[2, k])
      s <- window[is.finite(window[[a]]) & is.finite(window[[b]]), , drop = FALSE]
      f <- function(q) safe_cor(q[[a]], q[[b]], method = "kendall")
      ci <- boot_ci(s, f)
      data.frame(month = month.abb[m], pair = paste(prs[, k], collapse = "-"),
                 tau = f(s), lo = ci[1], hi = ci[2], n = nrow(s))
    }))
  }))
  write.csv(t1, file.path(rcfg$out_dir, "A2_tau_by_month.csv"), row.names = FALSE)
  
  # ---- B2: share unexplained (1 - R^2) and partial correlations -------------
  cat("Part B 2/4: unexplained share and partial correlations\n")
  zc <- paste0("z_", ALL_NODES)
  upper_ij <- which(upper.tri(matrix(0, length(zc), length(zc))), arr.ind = TRUE)
  unexpl_names <- paste0("unexpl_", ALL_NODES)
  pcor_names <- paste0("pcor_", ALL_NODES[upper_ij[, 1]], "-", ALL_NODES[upper_ij[, 2]])
  t2 <- do.call(rbind, lapply(1:12, function(m) {
    s <- rel[rel$month %in% win3(m), zc, drop = FALSE]
    s <- s[stats::complete.cases(s) & rowSums(!is.finite(as.matrix(s))) == 0L, , drop = FALSE]
    vals <- stats::setNames(rep(NA_real_, length(unexpl_names) + length(pcor_names)),
                            c(unexpl_names, pcor_names))
    if (nrow(s) >= ncol(s) + 1L) {
      R <- suppressWarnings(stats::cor(s))
      Ri <- tryCatch(if (all(is.finite(R))) solve(R) else NULL,
                     error = function(e) NULL)   # singular covariance: undefined partial r
      if (!is.null(Ri) && all(is.finite(Ri)) && all(diag(Ri) > 0)) {
        P <- -stats::cov2cor(Ri)
        vals[unexpl_names] <- 1 / diag(Ri)
        vals[pcor_names] <- P[upper.tri(P)]
      }
    }
    data.frame(month = month.abb[m], as.list(vals), n_complete = nrow(s), check.names = FALSE)
  }))
  write.csv(t2, file.path(rcfg$out_dir, "A2_unexplained_pcor.csv"), row.names = FALSE)
  
  # ---- B3: cascades, partial r(X(t-k), Y(t) | Y(t-1)), lags 1-6 -------------
  cat("Part B 3/4: cascades\n")
  lagv <- function(x, k) c(rep(NA, k), head(x, -k))
  mset <- list(NovMar = c(11, 12, 1, 2, 3), AprJun = 4:6, MayJul = 5:7, JulOct = 7:10)
  links <- list(c("Supply", "Snow", "NovMar"), c("Snow", "Soil", "AprJun"), c("Snow", "Flow", "MayJul"),
                c("Soil", "Flow", "AprJun"), c("Soil", "Flow", "JulOct"),
                c("Supply", "Soil", "JulOct"), c("Supply", "Flow", "JulOct"))
  t3 <- do.call(rbind, lapply(links, function(L) do.call(rbind, lapply(c("forward", "reverse"), function(dir) {
    X <- if (dir == "forward") L[1] else L[2]; Y <- if (dir == "forward") L[2] else L[1]
    do.call(rbind, lapply(1:6, function(k) {
      s <- data.frame(hyear = rel$hyear, month = rel$month, x = lagv(rel[[paste0("z_", X)]], k),
                      y = rel[[paste0("z_", Y)]], ylag1 = lagv(rel[[paste0("z_", Y)]], 1))
      s <- s[s$month %in% mset[[L[3]]] & stats::complete.cases(s) &
               is.finite(s$x) & is.finite(s$y) & is.finite(s$ylag1), , drop = FALSE]
      f <- function(q) {
        if (nrow(q) < 4L) return(NA_real_)
        tryCatch(safe_cor(stats::resid(stats::lm(y ~ ylag1, q)),
                          stats::resid(stats::lm(x ~ ylag1, q))),
                 error = function(e) NA_real_)
      }
      ci <- boot_ci(s, f)
      data.frame(link = paste(L[1], "->", L[2]), months = L[3], direction = dir, lag = k,
                 partial_r = f(s), lo = ci[1], hi = ci[2], n = nrow(s))
    }))
  }))))
  write.csv(t3, file.path(rcfg$out_dir, "A2_cascades.csv"), row.names = FALSE)
  
  # ---- B4: mediation, conditional tau of (A, C) given B ---------------------
  cat("Part B 4/4: mediation (slowest part)\n")
  clamp <- function(h) pmin(pmax(h, 1e-10), 1 - 1e-10)
  h_given_2 <- function(f, u) {                 # P(U1 <= u1 | U2 = u2)
    h <- tryCatch(rvinecopulib::hbicop(u, cond_var = 2, family = f), error = function(e) NULL)
    if (is.null(h)) h <- predict(f, newdata = u, what = "hfunc2")
    clamp(as.numeric(h))
  }
  hfun <- function(s, X, Bn) {
    # Keep original row alignment while filtering incomplete/non-finite pairs
    # *before* bicop() and h-function evaluation.
    uX <- as.matrix(s[, paste0("u_", c(X, Bn)), drop = FALSE])
    ok <- stats::complete.cases(uX) & rowSums(!is.finite(uX)) == 0L
    out <- rep(NA_real_, nrow(s))
    if (sum(ok) < 5L) return(out)
    u_ok <- uX[ok, , drop = FALSE]
    f <- tryCatch(rvinecopulib::bicop(u_ok, family_set = "parametric", selcrit = "bic"),
                  error = function(e) NULL)
    if (is.null(f)) return(out)
    h <- tryCatch(h_given_2(f, u_ok), error = function(e) NULL)
    if (!is.null(h) && length(h) == nrow(u_ok)) out[ok] <- h
    out
  }
  med <- list(c("Snow", "Flow", "Soil", "melt"),   c("Supply", "Flow", "Soil", "melt"),
              c("Snow", "Soil", "Supply", "melt"), c("Soil", "Flow", "Snow", "melt"),
              c("Snow", "Soil", "Flow", "melt"),   c("Supply", "Flow", "Snow", "acc"),
              c("Supply", "Flow", "Soil", "acc"),  c("Snow", "Soil", "Supply", "acc"),
              c("Supply", "Flow", "Soil", "rec"))
  t4 <- do.call(rbind, lapply(med, function(M) {
    cat(sprintf("  %s: %s-%s | %s\n", M[4], M[1], M[2], M[3]))
    cols <- paste0("u_", M[1:3])
    s <- rel[rel$month %in% cfg$regimes[[M[4]]]$months, , drop = FALSE]
    s <- s[stats::complete.cases(s[, cols, drop = FALSE]) &
             rowSums(!is.finite(as.matrix(s[, cols, drop = FALSE]))) == 0L, , drop = FALSE]
    n <- nrow(s)
    ny <- length(unique(s$hyear))
    # Stable one-row schema even when a sample or fit is unavailable.
    empty <- function(n, ny) data.frame(
      regime = M[4], pair = paste(M[1], M[2], sep = "-"), given = M[3],
      tau_marginal = NA_real_, tau_conditional = NA_real_, lo = NA_real_, hi = NA_real_,
      ratio = NA_real_, dBIC_years = NA_real_, family = NA_character_,
      rotation = NA_integer_, n = n, n_hyears = ny,
      BIC_fitted_years = NA_real_, BIC_independence_years = 0)
    if (n < 5L || ny < 2L) return(empty(n, ny))
    
    ha <- hfun(s, M[1], M[3]); hc <- hfun(s, M[2], M[3])
    ok <- is.finite(ha) & is.finite(hc)
    s <- s[ok, , drop = FALSE]
    hmat <- cbind(ha[ok], hc[ok])
    n <- nrow(s); ny <- length(unique(s$hyear))
    if (n < 5L || ny < 2L) return(empty(n, ny))
    fc <- tryCatch(rvinecopulib::bicop(hmat, family_set = "parametric", selcrit = "bic"),
                   error = function(e) NULL)
    if (is.null(fc)) return(empty(n, ny))
    L <- tryCatch(stats::logLik(fc), error = function(e) NULL)
    if (is.null(L)) return(empty(n, ny))
    k <- attr(L, "df"); if (is.null(k)) k <- fc$npars
    if (length(k) != 1L || !is.finite(k) || !is.finite(as.numeric(L))) return(empty(n, ny))
    
    # Part A's effective-years BIC: BIC_fit = -2*logLik*n_years/n_rows + k*log(n_years).
    # The independence copula has logLik=0 and k=0, thus BIC_independence=0.
    # dBIC_years = BIC_independence - BIC_fit: positive favors dependence.
    bic_fit_years <- -2 * as.numeric(L) * ny / n + k * log(ny)
    bic_indep_years <- 0
    bic_gain_years <- bic_indep_years - bic_fit_years
    marg_tau <- safe_cor(s[[cols[1]]], s[[cols[2]]], method = "kendall")
    cond_tau <- safe_cor(hmat[, 1], hmat[, 2], method = "kendall")
    ratio <- if (is.finite(marg_tau) && abs(marg_tau) > 1e-10) cond_tau / marg_tau else NA_real_
    ci <- boot_ci(s, function(q)
      safe_cor(hfun(q, M[1], M[3]), hfun(q, M[2], M[3]), method = "kendall"),
      B = rcfg$B_med)
    data.frame(regime = M[4], pair = paste(M[1], M[2], sep = "-"), given = M[3],
               tau_marginal = marg_tau, tau_conditional = cond_tau, lo = ci[1], hi = ci[2],
               ratio = ratio, dBIC_years = bic_gain_years,
               family = fc$family, rotation = fc$rotation, n = n, n_hyears = ny,
               BIC_fitted_years = bic_fit_years, BIC_independence_years = bic_indep_years)
  }))
  write.csv(t4, file.path(rcfg$out_dir, "A2_mediation.csv"), row.names = FALSE)
  
  cat(sprintf("\nPart B done in %.1f min. Files in %s:\n",
              as.numeric(difftime(Sys.time(), t_rel, units = "mins")), rcfg$out_dir))
  print(list.files(rcfg$out_dir, pattern = "^A2_.*\\.csv$"))
}

# ==============================================================================
# PART C - README DIAGRAMS (Graphviz via DiagrammeR)
#   Writes <DIAG_README_DIR>/figures/<name>.png (for the README), .svg and the
#   DOT sources in figures/dot/<name>.gv. Numeric labels below are MANUAL,
#   including screening comparisons that this script does not recalculate.
#   The source manifest records the current output-table hashes and warns on
#   missing/changed tables; human comparison of labels is still required.
# ==============================================================================
if (RUN$diagrams) {
  cat("\n", strrep("=", 60), "\n  README DIAGRAMS\n", strrep("=", 60), "\n", sep = "")
  
  DIAG_DIR <- file.path(DIAG_README_DIR, "figures")
  dir.create(file.path(DIAG_DIR, "dot"), showWarnings = FALSE, recursive = TRUE)
  
  # There is no table in this script for every hard-coded screening value.
  # Keep traceability to result tables and warn instead of silently presenting
  # manually curated numbers as automatically validated ones.
  diag_source_files <- c(
    model_comparison = file.path(cfg$out_dir, "model_comparison.csv"),
    vine_edges = file.path(cfg$out_dir, "vine_edges.csv"),
    relationship_by_regime = file.path(cfg$out_dir, "relationship_by_regime.csv"),
    A2_tau_by_month = file.path(cfg$out_dir, "relationship_v5", "A2_tau_by_month.csv"),
    A2_cascades = file.path(cfg$out_dir, "relationship_v5", "A2_cascades.csv"),
    A2_mediation = file.path(cfg$out_dir, "relationship_v5", "A2_mediation.csv"))
  source_exists <- file.exists(diag_source_files)
  source_md5 <- rep(NA_character_, length(diag_source_files))
  source_md5[source_exists] <- as.character(tools::md5sum(diag_source_files[source_exists]))
  diagram_manifest <- data.frame(
    source = names(diag_source_files), file = unname(diag_source_files),
    exists = unname(source_exists), md5 = source_md5,
    labels_automatically_verified = FALSE, stringsAsFactors = FALSE)
  diagram_manifest_file <- file.path(DIAG_DIR, "diagram_source_manifest.csv")
  if (file.exists(diagram_manifest_file)) {
    old_manifest <- tryCatch(utils::read.csv(diagram_manifest_file, stringsAsFactors = FALSE),
                             error = function(e) NULL)
    if (!is.null(old_manifest) && all(c("source", "md5") %in% names(old_manifest))) {
      prev <- old_manifest$md5[match(diagram_manifest$source, old_manifest$source)]
      changed <- (is.na(prev) != is.na(diagram_manifest$md5)) |
        (!is.na(prev) & !is.na(diagram_manifest$md5) & prev != diagram_manifest$md5)
      changed[is.na(changed)] <- FALSE
      if (any(changed))
        warning("Diagram input tables changed since the last rendering: ",
                paste(diagram_manifest$source[changed], collapse = ", "),
                ". Review manual diagram labels before publication.", call. = FALSE)
    }
  }
  utils::write.csv(diagram_manifest, diagram_manifest_file, row.names = FALSE)
  if (any(!source_exists))
    warning("Part C source tables missing: ",
            paste(names(diag_source_files)[!source_exists], collapse = ", "), call. = FALSE)
  warning("Part C numeric edge labels are manually curated, NOT validated against current tables. ",
          "Review all labels (especially external screening results) and see diagram_source_manifest.csv.",
          call. = FALSE)
  
  render_dot <- function(dot, name, scale = 3) {          # scale 3 = ~216 dpi
    writeLines(dot, file.path(DIAG_DIR, "dot", paste0(name, ".gv")), useBytes = TRUE)
    svg <- DiagrammeRsvg::export_svg(DiagrammeR::grViz(dot))
    writeLines(svg, file.path(DIAG_DIR, paste0(name, ".svg")), useBytes = TRUE)
    m <- regmatches(svg, regexpr('width="[0-9.]+pt"', svg))
    w <- if (length(m)) as.numeric(gsub("[^0-9.]", "", m)) else 1000
    rsvg::rsvg_png(charToRaw(svg), file.path(DIAG_DIR, paste0(name, ".png")), width = round(w * scale))
    message("  wrote ", file.path(DIAG_DIR, paste0(name, ".png")))
  }
  
  # ---- C1. Candidate screening and node selection -------------------------------
  dot_screening <- r"---(
digraph screening {
  graph [rankdir=LR, newrank=true, fontname="Arial", fontsize=13, labelloc=t,
         nodesep=0.22, ranksep=0.85, pad=0.2,
         label=<<b>Candidate screening and node selection</b><br/><font point-size="10">Edge labels: Kendall &tau; range over calendar months (&plusmn;1-month windows, 1956&ndash;2025; Table 4.1)</font>>]
  node [fontname="Arial", fontsize=10, shape=box, style="rounded,filled", fillcolor=white, color=grey40, margin="0.12,0.06"]
  edge [fontname="Arial", fontsize=9, arrowhead=none, color=grey40, fontcolor=grey20]

  subgraph cluster_supply {
    label="Supply"; style=rounded; color="#984EA3"; fontcolor="#984EA3"; fontsize=11
    spi3    [label="SPI-3", fillcolor=grey93, fontcolor=grey45, color=grey65]
    speiera [label="SPEI-ERA5-3\n(ERA5-Land PEV)", style="rounded,dashed"]
    speipm  [label="SPEI-PM-3", fillcolor="#e5d8ee", color="#984EA3", penwidth=2.2]
    spi3    -> speipm [label="redundant\n0.83&ndash;0.97", color="grey35:invis:grey35"]
    speiera -> speipm [label="different PET\n0.67&ndash;0.83"]
  }
  subgraph cluster_snow {
    label="Snow"; style=rounded; color="#377EB8"; fontcolor="#377EB8"; fontsize=11
    swei [label="SWEI-3", style="rounded,dashed"]
    sspi [label="SSPI-1", fillcolor="#deebf7", color="#377EB8", penwidth=2.2]
    swei -> sspi [label="interchangeable Oct&ndash;Apr\n0.50&ndash;0.90"]
  }
  subgraph cluster_soil {
    label="Soil"; style=rounded; color="#4DAF4A"; fontcolor="#4DAF4A"; fontsize=11
    ssmi [label="SSMI-1 L1_3", fillcolor=grey93, fontcolor=grey45, color=grey65]
    l14  [label="SSI-SM-1 L1_4\n(0&ndash;289 cm)", style="rounded,dashed"]
    l13  [label="SSI-SM-1 L1_3\n(0&ndash;100 cm)", fillcolor="#e5f5e0", color="#4DAF4A", penwidth=2.2]
    ssmi -> l13 [label="duplicate\n0.96&ndash;0.99", color="grey35:invis:grey35"]
    l14  -> l13 [label="distinct Dec&ndash;Jul\n0.50&ndash;0.81"]
  }
  subgraph cluster_flow {
    label="Flow"; style=rounded; color="#E41A1C"; fontcolor="#E41A1C"; fontsize=11
    ja017   [label="SSI 08JA017", fillcolor=grey93, fontcolor=grey45, color=grey65]
    jc001   [label="SSI 08JC001\n(Vanderhoof)", style="rounded,dashed"]
    srilong [label="SRI-3 / -6 / -12", fillcolor=grey93, fontcolor=grey45, color=grey65]
    sri1    [label="SRI-1\n(ERA5-Land runoff)", style="rounded,dashed"]
    jc002   [label="SSI 08JC002\nnaturalized (Isle Pierre)", fillcolor="#fde0dd", color="#E41A1C", penwidth=2.2]
    ja017   -> jc001 [label="redundant\n0.88&ndash;0.94", color="grey35:invis:grey35"]
    jc001   -> jc002 [label="redundant\n0.75&ndash;0.87", color="grey35:invis:grey35"]
    sri1    -> jc002 [label="best substitute\n0.49&ndash;0.73"]
    srilong -> jc002 [label="weaker\n0.15&ndash;0.62"]
  }

  msdi [label=<<b>MSDI-N</b><br/><font point-size="9">4 complementary nodes,<br/>regime-based vine</font>>,
        fillcolor="#fff7bc", color=grey20, penwidth=2.5, fontsize=12]
  speipm -> msdi [arrowhead=normal, color=grey20, penwidth=1.4]
  sspi   -> msdi [arrowhead=normal, color=grey20, penwidth=1.4]
  l13    -> msdi [arrowhead=normal, color=grey20, penwidth=1.4]
  jc002  -> msdi [arrowhead=normal, color=grey20, penwidth=1.4]

  nrt [label="Near-real-time variant\n(SRI-1 in the Flow slot; Section 3.5)", style="rounded,dashed"]
  sri1 -> nrt [arrowhead=normal, style=dashed]

  subgraph cluster_out {
    label="Outside the index"; style="rounded,dashed"; color=grey60; fontcolor=grey40; fontsize=11
    derived [label="H-SWSI, O-SWSI, SDI, NTSDI\nderived / composite: excluded", fillcolor=grey93, fontcolor=grey45, color=grey65]
    valid   [label="SGI (wells), SSI 08JE001 (Stuart)\nvalidation only (Section 4.14)", style="rounded,dashed"]
    storage [label="Reservoir storage,\nregulated flow", fillcolor="#fff2cc", color="#bf9000"]
    managed [label="Managed-supply layer, coupled to\nMSDI-N by a bivariate copula (Section 9)", fillcolor="#fff2cc", color="#bf9000", penwidth=1.8]
    storage -> managed [arrowhead=normal, color="#bf9000", label="managed"]
  }
}
)---"

# ---- C2. Relationship map by regime --------------------------------------------
dot_relmap <- r"---(
digraph relmap {
  graph [rankdir=TB, newrank=true, fontname="Arial", fontsize=16, labelloc=t,
         nodesep=0.4, ranksep=0.8, pad=0.2,
         label=<<b>Relationship map by regime</b><br/><font point-size="12">Supply (input) &rarr; Snow, Soil (stores) &rarr; Flow (output). Labels: &tau; from the fitted vines (Table 9b), lag-1 partial r (Table 4.4) or marginal &rarr; conditional &tau; (Table 4.5)</font>>]
  node [fontname="Arial", fontsize=14, shape=box, style="rounded,filled", width=1.3, color=grey35]
  edge [fontname="Arial", fontsize=12, dir=none, fontcolor=grey15]

  subgraph cluster_acc {
    label=<<b>Accumulation (Nov&ndash;Apr)</b>>; style=rounded; color=grey60; fontsize=15
    a_sup  [label="Supply", fillcolor="#e5d8ee"]
    a_snow [label="Snow",   fillcolor="#deebf7"]
    a_soil [label="Soil",   fillcolor="#e5f5e0"]
    a_flow [label="Flow",   fillcolor="#fde0dd"]
    a_sup  -> a_snow [label="direct\n0.40&ndash;0.48", penwidth=3.6]
    a_sup  -> a_soil [label="via Flow\n(indep. Nov&ndash;Feb)", style=dotted, color=grey55, penwidth=1.5]
    a_sup  -> a_flow [label="direct (rain,\nrain-on-snow)\n0.20&ndash;0.39", penwidth=2.8]
    a_snow -> a_flow [label="substitution | Supply\n&minus;0.08 to &minus;0.27", style=dashed, color="#b2182b", fontcolor="#b2182b", penwidth=2.2]
    a_soil -> a_flow [label="direct (baseflow)\n0.37&ndash;0.60", penwidth=4]
    a_snow -> a_soil [label="substitution | Supply, Flow\n&minus;0.32 to &minus;0.36", style=dashed, color="#b2182b", fontcolor="#b2182b", penwidth=3, constraint=false]
  }

  subgraph cluster_melt {
    label=<<b>Melt (May&ndash;Jun)</b>>; style=rounded; color=grey60; fontsize=15
    m_sup  [label="Supply", fillcolor="#e5d8ee"]
    m_snow [label="Snow",   fillcolor="#deebf7"]
    m_soil [label="Soil",   fillcolor="#e5f5e0"]
    m_flow [label="Flow",   fillcolor="#fde0dd"]
    m_sup  -> m_soil [label="direct\n0.40", penwidth=3.4]
    m_sup  -> m_snow [label="none given\nSoil, Flow", style=dotted, color=grey70, fontcolor=grey45]
    m_sup  -> m_flow [label="via Soil\n0.16 &rarr; 0.01", style=dotted, color=grey55, penwidth=1.5]
    m_snow -> m_flow [label="cascade, lag-1 r 0.71\n(contemp. 0.21)", dir=forward, color="#2166ac", fontcolor="#2166ac", penwidth=2.4]
    m_soil -> m_flow [label="common melt response\n+ weak direct 0.25", penwidth=2.2]
    m_snow -> m_soil [label="direct, lower tail\n0.43 | Flow (lag-1 r 0.33)", penwidth=3.6, constraint=false]
  }

  subgraph cluster_rec {
    label=<<b>Recession (Jul&ndash;Oct)</b>>; style=rounded; color=grey60; fontsize=15
    r_sup  [label="Supply\n(dropped)",  style="rounded,filled,dashed", fillcolor=grey93, fontcolor=grey45, color=grey65]
    r_snow [label="Snow\n(excluded)",   style="rounded,filled,dashed", fillcolor=grey93, fontcolor=grey45, color=grey65]
    r_soil [label="Soil", fillcolor="#e5f5e0"]
    r_flow [label="Flow", fillcolor="#fde0dd"]
    r_sup  -> r_soil [label="near-redundant\n0.60&ndash;0.74", color="grey25:invis:grey25", penwidth=1.3]
    r_sup  -> r_flow [label="via Soil\n0.40 &rarr; 0.13", style=dotted, color=grey55, penwidth=1.5]
    r_soil -> r_flow [label="cascade, lag-1 r 0.44\n0.26 (Jul) &rarr; 0.50 (Oct)", dir=forward, color="#2166ac", fontcolor="#2166ac", penwidth=3.2]
    r_snow -> r_flow [label="spring memory\ncarried by Flow", style=dotted, color=grey70, fontcolor=grey45]
  }

  subgraph cluster_key {
    label=<<b>Edge types</b>>; style=rounded; color=grey70; fontsize=14
    node [shape=point, width=0.05, color=grey45]
    k1a -> k1b [label="direct\n(width ~ |&tau;|)", penwidth=3]
    k2a -> k2b [label="cascade\n(lagged)", dir=forward, color="#2166ac", fontcolor="#2166ac", penwidth=2.4]
    k3a -> k3b [label="substitution\n(negative)", style=dashed, color="#b2182b", fontcolor="#b2182b", penwidth=2.4]
    k4a -> k4b [label="mediated /\ncond. indep.", style=dotted, color=grey55, penwidth=1.5]
    k5a -> k5b [label="near-\nredundant", color="grey25:invis:grey25", penwidth=1.3]
  }
     // 2 x 2 layout: Recession below Accumulation, key below Melt
     a_flow -> r_sup [style=invis]
     m_flow -> k1a   [style=invis]
}
)---"

# ---- C3. Fitted vine structures -------------------------------------------------
dot_vine_acc <- r"---(
graph vine_acc {
  graph [rankdir=TB, newrank=true, fontname="Arial", fontsize=13, labelloc=t,
         nodesep=0.45, ranksep=0.75, pad=0.2,
         label=<<b>(a) Accumulation (Nov&ndash;Apr)</b>: design D-vine Snow &ndash; Supply &ndash; Flow &ndash; Soil<br/><font point-size="10">Edge labels: Kendall &tau; range over the six monthly windows (families vary by window; Table 9b)</font>>]
  node [fontname="Arial", fontsize=11, shape=box, style="rounded,filled", fillcolor=white, color=grey35]
  edge [fontname="Arial", fontsize=10, fontcolor=grey15]

  T1 [label="Tree 1", shape=plaintext, style="", fontcolor=grey40, fontsize=12]
  T2 [label="Tree 2", shape=plaintext, style="", fontcolor=grey40, fontsize=12]
  T3 [label="Tree 3", shape=plaintext, style="", fontcolor=grey40, fontsize=12]
  T1 -- T2 -- T3 [style=invis]

  snow [label="Snow",   fillcolor="#deebf7"]
  sup  [label="Supply", fillcolor="#e5d8ee"]
  flow [label="Flow",   fillcolor="#fde0dd"]
  soil [label="Soil",   fillcolor="#e5f5e0"]
  { rank=same; T1; snow; sup; flow; soil }
  T1 -- snow [style=invis]
  snow -- sup  [label="0.40&ndash;0.48", penwidth=3.6]
  sup  -- flow [label="0.20&ndash;0.39", penwidth=2.8]
  flow -- soil [label="0.37&ndash;0.60", penwidth=4]

  e12 [label="Snow, Supply", fillcolor=grey96]
  e23 [label="Supply, Flow", fillcolor=grey96]
  e34 [label="Flow, Soil",   fillcolor=grey96]
  { rank=same; T2; e12; e23; e34 }
  T2 -- e12 [style=invis]
  e12 -- e23 [label="Snow, Flow | Supply\n&minus;0.08 to &minus;0.27", style=dashed, color="#b2182b", fontcolor="#b2182b", penwidth=2.4]
  e23 -- e34 [label="Supply, Soil | Flow\nindep. Nov&ndash;Feb;\n&minus;0.17, &minus;0.21 Mar&ndash;Apr", style=dotted, color=grey50, penwidth=1.6]

  e13 [label="Snow, Flow | Supply", fillcolor=grey96]
  e24 [label="Supply, Soil | Flow", fillcolor=grey96]
  { rank=same; T3; e13; e24 }
  T3 -- e13 [style=invis]
  e13 -- e24 [label="Snow, Soil | Supply, Flow\nGauss &minus;0.32 to &minus;0.36", style=dashed, color="#b2182b", fontcolor="#b2182b", penwidth=3]

  sup -- e12 [style=invis]; flow -- e23 [style=invis]; soil -- e34 [style=invis]
  e23 -- e13 [style=invis]; e34 -- e24 [style=invis]
}
)---"

dot_vine_melt <- r"---(
graph vine_melt {
  graph [rankdir=TB, newrank=true, fontname="Arial", fontsize=13, labelloc=t,
         nodesep=0.45, ranksep=0.75, pad=0.2,
         label=<<b>(b) Melt (May&ndash;Jun)</b>: design D-vine Supply &ndash; Soil &ndash; Flow &ndash; Snow<br/><font point-size="10">One fit window (May&ndash;Jun); family and Kendall &tau; per edge</font>>]
  node [fontname="Arial", fontsize=11, shape=box, style="rounded,filled", fillcolor=white, color=grey35]
  edge [fontname="Arial", fontsize=10, fontcolor=grey15]

  T1 [label="Tree 1", shape=plaintext, style="", fontcolor=grey40, fontsize=12]
  T2 [label="Tree 2", shape=plaintext, style="", fontcolor=grey40, fontsize=12]
  T3 [label="Tree 3", shape=plaintext, style="", fontcolor=grey40, fontsize=12]
  T1 -- T2 -- T3 [style=invis]

  sup  [label="Supply", fillcolor="#e5d8ee"]
  soil [label="Soil",   fillcolor="#e5f5e0"]
  flow [label="Flow",   fillcolor="#fde0dd"]
  snow [label="Snow",   fillcolor="#deebf7"]
  { rank=same; T1; sup; soil; flow; snow }
  T1 -- sup [style=invis]
  sup  -- soil [label="Gauss 0.40", penwidth=3.4]
  soil -- flow [label="Frank 0.25", penwidth=2.5]
  flow -- snow [label="Gauss 0.21", penwidth=2.3]

  e12 [label="Supply, Soil", fillcolor=grey96]
  e23 [label="Soil, Flow",   fillcolor=grey96]
  e34 [label="Flow, Snow",   fillcolor=grey96]
  { rank=same; T2; e12; e23; e34 }
  T2 -- e12 [style=invis]
  e12 -- e23 [label="Supply, Flow | Soil\nindependence", style=dotted, color=grey50, penwidth=1.6]
  e23 -- e34 [label="Soil, Snow | Flow\nsurvival Gumbel 0.43\n(lower tail)", penwidth=3.6]

  e13 [label="Supply, Flow | Soil", fillcolor=grey96]
  e24 [label="Soil, Snow | Flow",   fillcolor=grey96]
  { rank=same; T3; e13; e24 }
  T3 -- e13 [style=invis]
  e13 -- e24 [label="Supply, Snow | Soil, Flow\nindependence", style=dotted, color=grey50, penwidth=1.6]

  soil -- e12 [style=invis]; flow -- e23 [style=invis]; snow -- e34 [style=invis]
  e23 -- e13 [style=invis]; e34 -- e24 [style=invis]
}
)---"

dot_vine_rec <- r"---(
graph vine_rec {
  graph [rankdir=LR, fontname="Arial", fontsize=13, labelloc=t, nodesep=0.4, ranksep=1.8, pad=0.2,
         label=<<b>(c) Recession (Jul&ndash;Oct)</b>: bivariate copula Soil &ndash; Flow, one fit per month>]
  node [fontname="Arial", fontsize=11, shape=box, style="rounded,filled", fillcolor=white, color=grey35]
  edge [fontname="Arial", fontsize=10, fontcolor=grey15]
  soil [label="Soil", fillcolor="#e5f5e0"]
  flow [label="Flow", fillcolor="#fde0dd"]
  soil -- flow [label="Jul  Gumbel 0.26\lAug  Gumbel 0.32\lSep  Gauss 0.42\lOct  survival Gumbel 0.50\l", penwidth=3.2]
  sup  [label="Supply: dropped\n(near-redundant with Soil)", style="rounded,filled,dashed", fillcolor=grey93, fontcolor=grey45, color=grey65]
  snow [label="Snow: excluded\n(no current storage)",       style="rounded,filled,dashed", fillcolor=grey93, fontcolor=grey45, color=grey65]
}
)---"

render_dot(dot_screening, "fig_screening")
render_dot(dot_relmap,    "fig_relationship_map")
render_dot(dot_vine_acc,  "fig_vine_accumulation")
render_dot(dot_vine_melt, "fig_vine_melt")
render_dot(dot_vine_rec,  "fig_vine_recession")
}

cat(sprintf("\nAll requested parts done in %.1f min.\n",
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
# ==============================================================================
# END
# ==============================================================================