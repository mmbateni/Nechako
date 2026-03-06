####################################################################################
# utils_teleconnection_addon.R  ·  TELECONNECTION DATA UTILITIES
# Source this file AFTER sourcing DROUGHT_ANALYSIS_utils.R
####################################################################################
# Provides:
#   A.  Source URL constants for ONI, PDO, PNA (NOAA CPC / PSL)
#   B.  Parsers for each index's raw ASCII format
#   C.  load_teleconnection()  — unified loader with caching
#   D.  build_analysis_dataframe() — aligns drought index + teleconnections at lag τ
#   E.  assumption_checks()  — VIF, normality, homoscedasticity before lavaan
#
# USAGE:
#   source("DROUGHT_ANALYSIS_utils.R")
#   source("utils_teleconnection_addon.R")
#   tc <- load_teleconnection("oni")      # returns data.frame(date, value)
#   df <- build_analysis_dataframe(...)   # returns wide analysis frame
####################################################################################

## ── A. CONFIGURATION ─────────────────────────────────────────────────────────

# Output directory for downloaded raw files and processed CSV
TELE_DIR <- file.path(WD_PATH, "teleconnections")

# NOAA source URLs (publicly accessible, no authentication needed)
TELE_URLS <- list(
  oni = "https://www.cpc.ncep.noaa.gov/data/indices/oni.ascii.txt",
  pdo = "https://psl.noaa.gov/data/correlation/pdo.data",
  pna = "https://www.cpc.ncep.noaa.gov/products/precip/CWlink/pna/norm.pna.monthly.b5001.current.ascii.table"
)

# Standard data start years for each index
TELE_START_YEARS <- list(oni = 1950, pdo = 1900, pna = 1950)

## ── B. RAW FORMAT PARSERS ────────────────────────────────────────────────────

#' Parse ONI ASCII (NOAA CPC format)
#' Format: Season  YR  TOTAL  ANOM — 3-month running means
#' We extract the centre month of each season to get monthly values.
#'
#' @param raw_lines  Character vector from readLines()
#' @return data.frame(date, oni) monthly, one row per centre month
parse_oni_raw <- function(raw_lines) {
  # Drop header lines (any line that doesn't start with 3 letters)
  data_lines <- raw_lines[grepl("^[A-Z]{3}\\s+\\d{4}", raw_lines)]
  if (!length(data_lines)) stop("ONI: no data lines matched expected pattern")
  
  # Season abbreviation → centre month number
  season_centre <- c(
    DJF = 1, JFM = 2, FMA = 3, MAM = 4, AMJ = 5, MJJ = 6,
    JJA = 7, JAS = 8, ASO = 9, SON = 10, OND = 11, NDJ = 12
  )
  
  rows <- do.call(rbind, lapply(data_lines, function(ln) {
    parts <- strsplit(trimws(ln), "\\s+")[[1]]
    if (length(parts) < 4) return(NULL)
    season <- parts[1]; yr <- as.integer(parts[2])
    anom   <- suppressWarnings(as.numeric(parts[4]))
    cm     <- season_centre[season]
    if (is.na(cm) || is.na(yr) || is.na(anom)) return(NULL)
    data.frame(date  = as.Date(sprintf("%d-%02d-01", yr, cm)),
               oni   = anom, stringsAsFactors = FALSE)
  }))
  
  if (is.null(rows) || nrow(rows) == 0) stop("ONI parsing produced no rows")
  rows <- rows[order(rows$date), ]
  # Remove duplicates (NDJ wraps across years — keep first occurrence)
  rows <- rows[!duplicated(rows$date), ]
  rows
}

#' Parse PDO ASCII (NOAA PSL format)
#' Format: year  Jan  Feb  ...  Dec  — one year per row, missing = -9.99
#'
#' @param raw_lines Character vector from readLines()
#' @return data.frame(date, pdo)
parse_pdo_raw <- function(raw_lines) {
  # Header line contains "PDO" or starts with "Year"
  # Data lines start with a 4-digit year
  data_lines <- raw_lines[grepl("^\\s*\\d{4}\\s+", raw_lines)]
  if (!length(data_lines)) stop("PDO: no data lines matched expected pattern")
  
  month_names_short <- c("Jan","Feb","Mar","Apr","May","Jun",
                         "Jul","Aug","Sep","Oct","Nov","Dec")
  rows <- do.call(rbind, lapply(data_lines, function(ln) {
    parts <- suppressWarnings(as.numeric(strsplit(trimws(ln), "\\s+")[[1]]))
    if (length(parts) < 13 || is.na(parts[1])) return(NULL)
    yr     <- as.integer(parts[1])
    values <- parts[2:13]
    values[abs(values) >= 9.99] <- NA_real_   # PSL missing-value code
    do.call(rbind, lapply(1:12, function(m) {
      data.frame(date = as.Date(sprintf("%d-%02d-01", yr, m)),
                 pdo  = values[m], stringsAsFactors = FALSE)
    }))
  }))
  
  if (is.null(rows) || nrow(rows) == 0) stop("PDO parsing produced no rows")
  rows <- rows[order(rows$date), ]
  rows <- rows[!duplicated(rows$date), ]
  rows
}

#' Parse PNA ASCII (NOAA CPC format)
#' Format: header row + year rows, 12 monthly values per row
#' Missing = -99.90
#'
#' @param raw_lines Character vector from readLines()
#' @return data.frame(date, pna)
parse_pna_raw <- function(raw_lines) {
  data_lines <- raw_lines[grepl("^\\s*\\d{4}\\s+", raw_lines)]
  if (!length(data_lines)) stop("PNA: no data lines matched expected pattern")
  
  rows <- do.call(rbind, lapply(data_lines, function(ln) {
    parts <- suppressWarnings(as.numeric(strsplit(trimws(ln), "\\s+")[[1]]))
    if (length(parts) < 13 || is.na(parts[1])) return(NULL)
    yr     <- as.integer(parts[1])
    values <- parts[2:13]
    values[abs(values) >= 99.0] <- NA_real_
    do.call(rbind, lapply(1:12, function(m) {
      data.frame(date = as.Date(sprintf("%d-%02d-01", yr, m)),
                 pna  = values[m], stringsAsFactors = FALSE)
    }))
  }))
  
  if (is.null(rows) || nrow(rows) == 0) stop("PNA parsing produced no rows")
  rows <- rows[order(rows$date), ]
  rows <- rows[!duplicated(rows$date), ]
  rows
}

## ── C. UNIFIED LOADER WITH CACHING ──────────────────────────────────────────

#' Load a teleconnection index with local file caching.
#'
#' On first call, downloads the raw ASCII file from NOAA, parses it, and
#' saves a clean CSV to TELE_DIR.  Subsequent calls read from the cache.
#'
#' @param index_name   "oni", "pdo", or "pna"
#' @param start_year   Integer; earliest year to retain (default from TELE_START_YEARS)
#' @param end_year     Integer; latest year to retain (default: current year)
#' @param force_download Logical; if TRUE re-download even if cache exists
#' @return data.frame(date [Date], value [numeric])
load_teleconnection <- function(index_name,
                                start_year   = NULL,
                                end_year     = as.integer(format(Sys.Date(), "%Y")),
                                force_download = FALSE) {
  idx <- tolower(index_name)
  if (!idx %in% names(TELE_URLS))
    stop(sprintf("Unknown teleconnection index '%s'. Choose: %s",
                 index_name, paste(names(TELE_URLS), collapse = ", ")))
  
  dir.create(TELE_DIR, showWarnings = FALSE, recursive = TRUE)
  cache_file <- file.path(TELE_DIR, sprintf("%s_monthly.csv", idx))
  
  # ── Use cache if present and not forced ─────────────────────────────────
  if (file.exists(cache_file) && !force_download) {
    cat(sprintf("  ✓ Loading %s from cache: %s\n", toupper(idx), basename(cache_file)))
    df <- utils::read.csv(cache_file, stringsAsFactors = FALSE)
    df$date <- as.Date(df$date)
    colnames(df)[2] <- "value"
  } else {
    # ── Download ─────────────────────────────────────────────────────────
    raw_file <- file.path(TELE_DIR, sprintf("%s_raw.txt", idx))
    cat(sprintf("  ↓ Downloading %s from NOAA...\n", toupper(idx)))
    tryCatch(
      utils::download.file(TELE_URLS[[idx]], destfile = raw_file, quiet = TRUE),
      error = function(e) stop(sprintf(
        "Download failed for %s.\n  URL: %s\n  Error: %s\n  → Download manually and place at: %s",
        toupper(idx), TELE_URLS[[idx]], e$message, raw_file))
    )
    
    # ── Parse ─────────────────────────────────────────────────────────────
    raw_lines <- readLines(raw_file, warn = FALSE)
    parsed <- switch(idx,
                     oni = parse_oni_raw(raw_lines),
                     pdo = parse_pdo_raw(raw_lines),
                     pna = parse_pna_raw(raw_lines),
                     stop("No parser for: ", idx)
    )
    
    df <- parsed
    colnames(df)[2] <- "value"
    
    # ── Cache ─────────────────────────────────────────────────────────────
    utils::write.csv(df, cache_file, row.names = FALSE)
    cat(sprintf("  ✓ Parsed & cached: %d rows (%s to %s)\n",
                nrow(df), min(df$date), max(df$date)))
  }
  
  # ── Filter by year range ─────────────────────────────────────────────────
  sy <- if (!is.null(start_year)) start_year else TELE_START_YEARS[[idx]]
  df <- df[as.integer(format(df$date, "%Y")) >= sy &
             as.integer(format(df$date, "%Y")) <= end_year, ]
  df <- df[!is.na(df$value), ]
  df <- df[order(df$date), ]
  rownames(df) <- NULL
  cat(sprintf("  ✓ %s ready: %d months (%s – %s)\n",
              toupper(idx), nrow(df), min(df$date), max(df$date)))
  df
}

## ── D. ANALYSIS DATAFRAME BUILDER ───────────────────────────────────────────

#' Align drought index with lagged teleconnection predictors.
#'
#' For moderated mediation the outcome Y is the drought index measured τ months
#' AFTER the teleconnection observations.  This function:
#'   1. Loads the basin-averaged drought index (via load_basin_avg_csv)
#'   2. Loads ONI, PDO, PNA (via load_teleconnection)
#'   3. Merges on date with the drought index shifted forward by `lag_months`
#'   4. Adds centred/scaled versions (suffix _c) and interaction term (ONI_PDO_int,
#'      PNA_PDO_int)
#'   5. Returns a complete-cases data.frame ready for lavaan
#'
#' @param index_type   "spi", "spei", or "swei"
#' @param scale        Integer timescale (e.g. 6 for SPEI-6)
#' @param lag_months   Integer; drought index is measured lag_months after
#'                     the teleconnection month (0 = contemporaneous)
#' @param start_year   Integer; first year to include (default 1950)
#' @param end_year     Integer; last year to include (default 2024)
#' @return data.frame with columns:
#'         date, SPI_outcome, ONI, PDO, PNA,
#'         ONI_c, PDO_c, PNA_c, ONI_PDO_int, PNA_PDO_int,
#'         pdv_phase (factor: Cool / Neutral / Warm)
build_analysis_dataframe <- function(index_type,
                                     scale,
                                     lag_months  = 0L,
                                     start_year  = 1950L,
                                     end_year    = 2024L) {
  
  # ── Load drought index ───────────────────────────────────────────────────
  di <- load_basin_avg_csv(index_type, scale)
  if (is.null(di) || nrow(di) == 0)
    stop(sprintf("Drought index %s-%d not found. Run w1_basin_timeseries.R first.",
                 toupper(index_type), scale))
  colnames(di) <- c("date", "drought_value")
  di$date <- as.Date(di$date)
  
  # ── Load teleconnection indices ──────────────────────────────────────────
  oni_df <- load_teleconnection("oni", start_year = start_year, end_year = end_year + 1)
  pdo_df <- load_teleconnection("pdo", start_year = start_year, end_year = end_year + 1)
  pna_df <- load_teleconnection("pna", start_year = start_year, end_year = end_year + 1)
  
  colnames(oni_df) <- c("date", "ONI")
  colnames(pdo_df) <- c("date", "PDO")
  colnames(pna_df) <- c("date", "PNA")
  
  # ── Merge teleconnections (same month t) ─────────────────────────────────
  tele <- merge(oni_df, pdo_df, by = "date", all = FALSE)
  tele <- merge(tele,   pna_df, by = "date", all = FALSE)
  
  # ── Lag: drought index at t + lag_months ─────────────────────────────────
  # Strategy: shift drought dates BACK by lag_months so they align with the
  # teleconnection date when we merge on date.
  # i.e., drought measured at month (t + lag) ↔ tele measured at month t
  if (lag_months > 0L) {
    di$date_tele <- di$date - months(lag_months)   # align key for merge
  } else {
    di$date_tele <- di$date
  }
  
  # ── Merge ─────────────────────────────────────────────────────────────────
  df <- merge(tele,
              di[, c("date_tele", "drought_value")],
              by.x = "date", by.y = "date_tele",
              all.x = FALSE)
  
  df <- df[as.integer(format(df$date, "%Y")) >= start_year &
             as.integer(format(df$date, "%Y")) <= end_year, ]
  df <- df[stats::complete.cases(df), ]
  df <- df[order(df$date), ]
  
  if (nrow(df) < 30)
    warning(sprintf(
      "Only %d complete rows after alignment (lag=%d). Results may be unreliable.",
      nrow(df), lag_months))
  
  # ── Standardise all continuous variables ─────────────────────────────────
  # Centring and scaling is required for lavaan moderated mediation because:
  #  (a) interaction terms (product terms) are collinear with main effects
  #      unless centred; (b) standardised coefficients are directly comparable
  #      across indices and lags.
  df$ONI_c <- as.numeric(scale(df$ONI))
  df$PDO_c <- as.numeric(scale(df$PDO))
  df$PNA_c <- as.numeric(scale(df$PNA))
  df$Y_c   <- as.numeric(scale(df$drought_value))   # outcome
  
  # ── Interaction (product) terms ──────────────────────────────────────────
  # These are the W×X and W×M terms in moderated mediation.
  # Because all components are centred, multicollinearity is minimised.
  df$ONI_PDO_int <- df$ONI_c * df$PDO_c   # moderates the a-path (ONI → PNA)
  df$PNA_PDO_int <- df$PNA_c * df$PDO_c   # moderates the b-path (PNA → Y)
  
  # ── PDO phase as categorical variable ────────────────────────────────────
  df$pdv_phase <- cut(df$PDO,
                      breaks = c(-Inf, -0.5, 0.5, Inf),
                      labels = c("Cool", "Neutral", "Warm"),
                      right  = FALSE)
  
  # ── Attach metadata ───────────────────────────────────────────────────────
  attr(df, "index_label") <- sprintf("%s-%d", toupper(index_type), scale)
  attr(df, "lag_months")  <- lag_months
  cat(sprintf("  ✓ Analysis frame: %d rows | index: %s | lag: %d months\n",
              nrow(df), attr(df, "index_label"), lag_months))
  df
}

## ── E. ASSUMPTION CHECKS ─────────────────────────────────────────────────────

#' Run pre-modelling assumption checks and print a diagnostic summary.
#'
#' Checks: (1) VIF on main effects, (2) Shapiro-Wilk normality on all
#' continuous variables, (3) Breusch-Pagan homoscedasticity in the a-path
#' and b-path preliminary OLS regressions, (4) Durbin-Watson autocorrelation.
#'
#' @param df   Output of build_analysis_dataframe()
#' @return Invisible list of test results; prints a formatted summary
check_mm_assumptions <- function(df) {
  utils_load_packages(c("car", "lmtest", "nortest"))
  
  cat("\n══════════════════════════════════════════════════════════\n")
  cat("  PRE-MODELLING ASSUMPTION CHECKS\n")
  cat("══════════════════════════════════════════════════════════\n")
  
  results <- list()
  
  # 1. Variance Inflation Factor ───────────────────────────────────────────
  cat("\n── 1. VIF (multicollinearity check on main effects) ──\n")
  # Fit OLS with main effects + interaction terms
  ols_y <- lm(Y_c ~ ONI_c + PDO_c + PNA_c + ONI_PDO_int + PNA_PDO_int, data = df)
  tryCatch({
    vif_vals <- car::vif(ols_y)
    results$vif <- vif_vals
    print(round(vif_vals, 3))
    if (any(vif_vals > 10))
      cat("  ⚠ VIF > 10 detected — check centring of interaction terms\n")
    else if (any(vif_vals > 5))
      cat("  ⚠ VIF > 5 for some terms — moderate multicollinearity\n")
    else
      cat("  ✓ All VIF < 5 — acceptable multicollinearity\n")
  }, error = function(e) cat("  VIF check failed:", e$message, "\n"))
  
  # 2. Normality of key variables ──────────────────────────────────────────
  cat("\n── 2. Shapiro-Wilk normality (n ≤ 5000 subsample) ──\n")
  vars_to_check <- c("ONI_c", "PDO_c", "PNA_c", "Y_c")
  norm_results <- data.frame(Variable = vars_to_check,
                             W = NA_real_, p_value = NA_real_,
                             Normal = NA_character_,
                             stringsAsFactors = FALSE)
  for (i in seq_along(vars_to_check)) {
    v   <- df[[vars_to_check[i]]]
    v   <- v[is.finite(v)]
    sub <- if (length(v) > 5000) sample(v, 5000) else v
    sw  <- tryCatch(shapiro.test(sub), error = function(e) NULL)
    if (!is.null(sw)) {
      norm_results$W[i]      <- round(sw$statistic, 4)
      norm_results$p_value[i] <- round(sw$p.value, 4)
      norm_results$Normal[i]  <- if (sw$p.value > 0.05) "Yes (p>0.05)" else "No (p≤0.05)"
    }
  }
  print(norm_results)
  cat("  Note: lavaan bootstrap CIs are robust to non-normality.\n")
  results$normality <- norm_results
  
  # 3. Homoscedasticity (a-path and b-path OLS) ────────────────────────────
  cat("\n── 3. Breusch-Pagan homoscedasticity ──\n")
  ols_a <- lm(PNA_c ~ ONI_c + PDO_c + ONI_PDO_int, data = df)
  ols_b <- lm(Y_c   ~ PNA_c + ONI_c + PDO_c + PNA_PDO_int, data = df)
  
  for (nm in c("a-path (ONI→PNA)", "b-path (PNA→Y)")) {
    mdl <- if (nm == "a-path (ONI→PNA)") ols_a else ols_b
    bp  <- tryCatch(lmtest::bptest(mdl), error = function(e) NULL)
    if (!is.null(bp)) {
      cat(sprintf("  %s:  BP = %.3f, p = %.4f  %s\n",
                  nm, bp$statistic, bp$p.value,
                  if (bp$p.value > 0.05) "✓ Homoscedastic" else "⚠ Heteroscedastic → use HC3 SEs"))
    }
  }
  
  # 4. Durbin-Watson autocorrelation ───────────────────────────────────────
  cat("\n── 4. Durbin-Watson autocorrelation ──\n")
  for (nm in c("a-path", "b-path")) {
    mdl <- if (nm == "a-path") ols_a else ols_b
    dw  <- tryCatch(lmtest::dwtest(mdl), error = function(e) NULL)
    if (!is.null(dw)) {
      cat(sprintf("  %s:  DW = %.3f, p = %.4f  %s\n",
                  nm, dw$statistic, dw$p.value,
                  if (dw$p.value > 0.10) "✓ No autocorrelation" else "⚠ Serial correlation present"))
    }
  }
  
  cat("\n══════════════════════════════════════════════════════════\n")
  cat("  Note: Serial correlation is EXPECTED in monthly climate\n")
  cat("  data. lavaan bootstrap (5000 reps) handles this via\n")
  cat("  resampling, not distributional assumptions.\n")
  cat("══════════════════════════════════════════════════════════\n\n")
  
  invisible(results)
}

cat("✓ utils_teleconnection_addon.R loaded\n")