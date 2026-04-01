####################################################################################
# utils_teleconnection_addon.R  ·  TELECONNECTION DATA UTILITIES + GAM TOOLS
# Source this file AFTER sourcing DROUGHT_ANALYSIS_utils.R
####################################################################################
# Provides:
#   A.  URL constants for ONI, PDO, PNA (NOAA CPC / PSL)
#   B.  ASCII parsers for each index's raw format
#   C.  load_teleconnection()  — unified loader with caching
#   D.  build_analysis_dataframe() — aligns drought index + teleconnections at lag τ
#       NOTE: Returns raw ONI/PDO/PNA (for GAM) AND centred/scaled _c versions
#             (for lavaan MM) in a single data.frame
#       PERF: Accepts optional pre-merged `tele_df` to avoid re-loading from disk
#             on every call. Pass the TELE_MASTER object built in w8.
#   E.  check_mm_assumptions()  — VIF, normality, BP, DW before lavaan
#   F.  GAM UTILITIES:
#         fit_gam_full()           — full model with main smooths + ti() interactions
#         gam_variance_partition() — deviance explained per smooth term
#                                    PERF: accepts precomputed_s to skip summary()
#         gam_concurvity_check()   — concurvity (GAM analog of VIF)
#         predict_gam_surface()    — prediction grid for a 2-D response surface
#         plot_gam_surfaces()      — 3-panel ONI×PDO conditioned on PNA (key figure)
#         plot_gam_partials()      — partial effect plots for all smooths
#         plot_gam_residuals()     — residual diagnostics (4 panels)
#         gam_lag_summary()        — one-row summary for lag-profile table
#                                    PERF: accepts precomputed_s to skip summary()
####################################################################################

## ── A. CONFIGURATION ─────────────────────────────────────────────────────────

WD_PATH <- "D:/Nechako_Drought/Nechako/"
TELE_DIR <- file.path(WD_PATH, "teleconnections")

TELE_URLS <- list(
  oni    = "https://www.cpc.ncep.noaa.gov/data/indices/oni.ascii.txt",
  pdo    = "https://psl.noaa.gov/data/correlation/pdo.data",
  pna    = "https://www.cpc.ncep.noaa.gov/products/precip/CWlink/pna/norm.pna.monthly.b5001.current.ascii.table",
  amo    = "https://psl.noaa.gov/data/correlation/amo.data",
  nino34 = "https://psl.noaa.gov/data/correlation/nina34.data",
  ao     = "https://www.cpc.ncep.noaa.gov/products/precip/CWlink/daily_ao_index/monthly.ao.index.b50.current.ascii.table"
)

TELE_START_YEARS <- list(oni = 1950, pdo = 1900, pna = 1950,
                         amo = 1948, nino34 = 1950, ao = 1950)

## ── B. RAW FORMAT PARSERS ────────────────────────────────────────────────────

parse_oni_raw <- function(raw_lines) {
  data_lines <- raw_lines[grepl("^[A-Z]{3}\\s+\\d{4}", raw_lines)]
  if (!length(data_lines)) stop("ONI: no data lines matched expected pattern")
  season_centre <- c(DJF=1,JFM=2,FMA=3,MAM=4,AMJ=5,MJJ=6,
                     JJA=7,JAS=8,ASO=9,SON=10,OND=11,NDJ=12)
  rows <- dplyr::bind_rows(lapply(data_lines, function(ln) {
    parts <- strsplit(trimws(ln), "\\s+")[[1]]
    if (length(parts) < 4) return(NULL)
    season <- parts[1]; yr <- as.integer(parts[2])
    anom   <- suppressWarnings(as.numeric(parts[4]))
    cm     <- season_centre[season]
    if (is.na(cm) || is.na(yr) || is.na(anom)) return(NULL)
    data.frame(date=as.Date(sprintf("%d-%02d-01", yr, cm)),
               oni=anom, stringsAsFactors=FALSE)
  }))
  if (is.null(rows) || nrow(rows)==0) stop("ONI parsing produced no rows")
  rows <- rows[order(rows$date),]
  rows[!duplicated(rows$date),]
}

parse_pdo_raw <- function(raw_lines) {
  data_lines <- raw_lines[grepl("^\\s*\\d{4}\\s+", raw_lines)]
  if (!length(data_lines)) stop("PDO: no data lines matched")
  rows <- dplyr::bind_rows(lapply(data_lines, function(ln) {
    parts <- suppressWarnings(as.numeric(strsplit(trimws(ln), "\\s+")[[1]]))
    if (length(parts) < 13 || is.na(parts[1])) return(NULL)
    yr <- as.integer(parts[1]); values <- parts[2:13]
    values[abs(values) >= 9.99] <- NA_real_
    data.frame(date=as.Date(paste(yr, sprintf("%02d", 1:12), "01", sep="-")),
               pdo=values, stringsAsFactors=FALSE)
  }))
  if (is.null(rows) || nrow(rows)==0) stop("PDO parsing produced no rows")
  rows <- rows[order(rows$date),]
  rows[!duplicated(rows$date),]
}

parse_pna_raw <- function(raw_lines) {
  data_lines <- raw_lines[grepl("^\\s*\\d{4}\\s+", raw_lines)]
  if (!length(data_lines)) stop("PNA: no data lines matched")
  rows <- dplyr::bind_rows(lapply(data_lines, function(ln) {
    parts <- suppressWarnings(as.numeric(strsplit(trimws(ln), "\\s+")[[1]]))
    if (length(parts) < 13 || is.na(parts[1])) return(NULL)
    yr <- as.integer(parts[1]); values <- parts[2:13]
    values[abs(values) >= 99.0] <- NA_real_
    data.frame(date=as.Date(paste(yr, sprintf("%02d", 1:12), "01", sep="-")),
               pna=values, stringsAsFactors=FALSE)
  }))
  if (is.null(rows) || nrow(rows)==0) stop("PNA parsing produced no rows")
  rows <- rows[order(rows$date),]
  rows[!duplicated(rows$date),]
}

## ── Internal helper: year × 12-month table format (PSL/CPC standard) ────────
## Used by AMO, Nino34, AO parsers.  PDO and PNA keep their own parsers for
## backward-compatibility.
.parse_year_month_table <- function(raw_lines, na_threshold = 99.0) {
  data_lines <- raw_lines[grepl("^\\s*\\d{4}\\s+", raw_lines)]
  if (!length(data_lines)) stop("No data lines matched year × month pattern")
  rows <- dplyr::bind_rows(lapply(data_lines, function(ln) {
    parts <- suppressWarnings(as.numeric(strsplit(trimws(ln), "\\s+")[[1]]))
    if (length(parts) < 13 || is.na(parts[1])) return(NULL)
    yr     <- as.integer(parts[1])
    values <- parts[2:13]
    values[abs(values) >= na_threshold] <- NA_real_
    data.frame(date  = as.Date(paste(yr, sprintf("%02d", 1:12), "01", sep = "-")),
               value = values, stringsAsFactors = FALSE)
  }))
  if (is.null(rows) || nrow(rows) == 0) stop("Parsing produced no rows")
  rows <- rows[order(rows$date), ]
  rows[!duplicated(rows$date), ]
}

parse_amo_raw    <- function(raw_lines) .parse_year_month_table(raw_lines, na_threshold = 9.99)
parse_nino34_raw <- function(raw_lines) .parse_year_month_table(raw_lines, na_threshold = 9.99)
parse_ao_raw     <- function(raw_lines) .parse_year_month_table(raw_lines, na_threshold = 99.0)

## ── C. UNIFIED LOADER WITH CACHING ──────────────────────────────────────────

load_teleconnection <- function(index_name,
                                start_year     = NULL,
                                end_year       = as.integer(format(Sys.Date(),"%Y")),
                                force_download = FALSE) {
  idx <- tolower(index_name)
  if (!idx %in% names(TELE_URLS))
    stop(sprintf("Unknown index '%s'. Choose: %s",
                 index_name, paste(names(TELE_URLS), collapse=", ")))
  dir.create(TELE_DIR, showWarnings=FALSE, recursive=TRUE)
  cache_file <- file.path(TELE_DIR, sprintf("%s_monthly.csv", idx))
  
  if (file.exists(cache_file) && !force_download) {
    cat(sprintf("  ✓ %s from cache\n", toupper(idx)))
    df <- utils::read.csv(cache_file, stringsAsFactors=FALSE)
    df$date <- as.Date(df$date); colnames(df)[2] <- "value"
  } else {
    raw_file <- file.path(TELE_DIR, sprintf("%s_raw.txt", idx))
    cat(sprintf("  ↓ Downloading %s...\n", toupper(idx)))
    tryCatch(
      utils::download.file(TELE_URLS[[idx]], destfile=raw_file, quiet=TRUE),
      error = function(e) stop(sprintf(
        "Download failed for %s.\n  URL: %s\n  Error: %s\n  Place file at: %s",
        toupper(idx), TELE_URLS[[idx]], e$message, raw_file)))
    raw_lines <- readLines(raw_file, warn=FALSE)
    parsed <- switch(idx,
                     oni    = parse_oni_raw(raw_lines),
                     pdo    = parse_pdo_raw(raw_lines),
                     pna    = parse_pna_raw(raw_lines),
                     amo    = parse_amo_raw(raw_lines),
                     nino34 = parse_nino34_raw(raw_lines),
                     ao     = parse_ao_raw(raw_lines),
                     stop("No parser for: ", idx))
    df <- parsed; colnames(df)[2] <- "value"
    utils::write.csv(df, cache_file, row.names=FALSE)
    cat(sprintf("  ✓ Parsed & cached: %d rows\n", nrow(df)))
  }
  sy <- if (!is.null(start_year)) start_year else TELE_START_YEARS[[idx]]
  df <- df[as.integer(format(df$date,"%Y")) >= sy &
             as.integer(format(df$date,"%Y")) <= end_year,]
  df <- df[!is.na(df$value),]; df <- df[order(df$date),]
  rownames(df) <- NULL
  cat(sprintf("  ✓ %s: %d months (%s – %s)\n",
              toupper(idx), nrow(df), min(df$date), max(df$date)))
  df
}

## ── D. ANALYSIS DATAFRAME BUILDER ───────────────────────────────────────────

#' Build a single data.frame suitable for BOTH lavaan (MM) and mgcv (GAM).
#'
#' PERFORMANCE: Pass `tele_df` (a pre-merged ONI/PDO/PNA data frame built once
#' in w8 as TELE_MASTER) to skip the three load_teleconnection() calls on every
#' invocation.  When `tele_df=NULL` (default) the original download/cache path
#' is used — full backward-compatibility is preserved.
#'
#' @param tele_df  Optional pre-merged teleconnection frame with columns
#'                 date, ONI, PDO, PNA covering at least start_year:end_year+1.
#'                 Create with:
#'                   oni <- load_teleconnection("oni"); colnames(oni)[2]<-"ONI"
#'                   pdo <- load_teleconnection("pdo"); colnames(pdo)[2]<-"PDO"
#'                   pna <- load_teleconnection("pna"); colnames(pna)[2]<-"PNA"
#'                   TELE_MASTER <- Reduce(function(a,b) merge(a,b,by="date"),
#'                                         list(oni, pdo, pna))
build_analysis_dataframe <- function(index_type,
                                     scale,
                                     lag_months = 0L,
                                     start_year = 1950L,
                                     end_year   = 2024L,
                                     tele_df    = NULL) {
  di <- load_basin_avg_csv(index_type, scale)
  if (is.null(di) || nrow(di)==0)
    stop(sprintf("Index %s-%d not found. Run w1_basin_timeseries.R first.",
                 toupper(index_type), scale))
  colnames(di) <- c("date","drought_value"); di$date <- as.Date(di$date)
  
  # ── Teleconnection loading: use pre-built master if supplied ──────────────
  if (!is.null(tele_df)) {
    # Subset the pre-loaded master to the requested year window (+1 year buffer)
    tele_years <- as.integer(format(tele_df$date, "%Y"))
    tele <- tele_df[tele_years >= start_year & tele_years <= (end_year + 1L), ]
    if (nrow(tele) == 0)
      stop("tele_df has no rows in the requested year range — check TELE_MASTER.")
  } else {
    # Original path: load & merge from cache / network
    oni_df <- load_teleconnection("oni", start_year=start_year, end_year=end_year+1)
    pdo_df <- load_teleconnection("pdo", start_year=start_year, end_year=end_year+1)
    pna_df <- load_teleconnection("pna", start_year=start_year, end_year=end_year+1)
    colnames(oni_df) <- c("date","ONI")
    colnames(pdo_df) <- c("date","PDO")
    colnames(pna_df) <- c("date","PNA")
    tele <- merge(oni_df, pdo_df, by="date", all=FALSE)
    tele <- merge(tele,   pna_df, by="date", all=FALSE)
  }
  
  # BUG-FIX: lubridate::months() is not namespace-exported in lubridate >= 1.9
  # when the package is loaded but not attached (common inside future workers).
  # Use base-R POSIXlt month arithmetic instead — no external dependencies.
  di$date_tele <- if (lag_months > 0L) {
    lt      <- as.POSIXlt(di$date)
    lt$mon  <- lt$mon - as.integer(lag_months)
    as.Date(lt)
  } else {
    di$date
  }
  
  df <- merge(tele, di[, c("date_tele","drought_value")],
              by.x="date", by.y="date_tele", all.x=FALSE)
  df <- df[as.integer(format(df$date,"%Y")) >= start_year &
             as.integer(format(df$date,"%Y")) <= end_year,]
  df <- df[stats::complete.cases(df),]; df <- df[order(df$date),]
  if (nrow(df) < 30)
    warning(sprintf("Only %d complete rows (lag=%d).", nrow(df), lag_months))
  
  # Centred/scaled variables — for lavaan MM
  df$ONI_c <- as.numeric(scale(df$ONI))
  df$PDO_c <- as.numeric(scale(df$PDO))
  df$PNA_c <- as.numeric(scale(df$PNA))
  df$Y_c   <- as.numeric(scale(df$drought_value))
  
  # Product interaction terms — for lavaan MM only
  df$ONI_PDO_int <- df$ONI_c * df$PDO_c
  df$PNA_PDO_int <- df$PNA_c * df$PDO_c
  
  # Categorical PDO phase — used in both MM and GAM
  df$pdo_phase <- cut(df$PDO, breaks=c(-Inf,-0.5,0.5,Inf),
                      labels=c("Cool","Neutral","Warm"), right=FALSE)
  
  # Calendar month + year — for GAM seasonality smooth and time trends
  df$month <- as.integer(format(df$date,"%m"))
  df$year  <- as.integer(format(df$date,"%Y"))
  
  attr(df,"index_label") <- sprintf("%s-%d", toupper(index_type), scale)
  attr(df,"lag_months")  <- lag_months
  cat(sprintf("  ✓ Analysis frame: %d rows | %s | lag %d months\n",
              nrow(df), attr(df,"index_label"), lag_months))
  df
}

## ── E. LAVAAN ASSUMPTION CHECKS ──────────────────────────────────────────────

check_mm_assumptions <- function(df) {
  utils_load_packages(c("car","lmtest"))
  cat("\n══════════════════════════════════════════════════════════\n")
  cat("  PRE-MODELLING ASSUMPTION CHECKS (lavaan MM)\n")
  cat("══════════════════════════════════════════════════════════\n")
  results <- list()
  
  cat("\n── 1. VIF ──\n")
  ols_y <- lm(Y_c ~ ONI_c + PDO_c + PNA_c + ONI_PDO_int + PNA_PDO_int, data=df)
  tryCatch({
    vif_vals <- car::vif(ols_y); results$vif <- vif_vals
    print(round(vif_vals,3))
    if (any(vif_vals>10)) cat("  ⚠ VIF>10 — check centring\n")
    else if (any(vif_vals>5)) cat("  ⚠ VIF 5-10 — moderate\n")
    else cat("  ✓ All VIF<5\n")
  }, error=function(e) cat("  VIF failed:", e$message, "\n"))
  
  cat("\n── 2. Shapiro-Wilk normality ──\n")
  vcols <- c("ONI_c","PDO_c","PNA_c","Y_c")
  norm_res <- data.frame(Variable=vcols, W=NA_real_, p=NA_real_,
                         Normal=NA_character_, stringsAsFactors=FALSE)
  for (i in seq_along(vcols)) {
    v <- df[[vcols[i]]]; v <- v[is.finite(v)]
    sub <- if (length(v)>5000) sample(v,5000) else v
    sw <- tryCatch(shapiro.test(sub), error=function(e) NULL)
    if (!is.null(sw)) {
      norm_res$W[i] <- round(sw$statistic,4); norm_res$p[i] <- round(sw$p.value,4)
      norm_res$Normal[i] <- if (sw$p.value>0.05) "Yes" else "No (→ bootstrap)"
    }
  }
  print(norm_res); results$normality <- norm_res
  
  cat("\n── 3. Breusch-Pagan & Durbin-Watson ──\n")
  ols_a <- lm(PNA_c ~ ONI_c + PDO_c + ONI_PDO_int, data=df)
  ols_b <- lm(Y_c   ~ PNA_c + ONI_c + PDO_c + PNA_PDO_int, data=df)
  for (nm in c("a-path","b-path")) {
    mdl <- if (nm=="a-path") ols_a else ols_b
    bp  <- tryCatch(lmtest::bptest(mdl), error=function(e) NULL)
    dw  <- tryCatch(lmtest::dwtest(mdl), error=function(e) NULL)
    if (!is.null(bp)) cat(sprintf("  %s BP: p=%.4f %s\n", nm, bp$p.value,
                                  if (bp$p.value>0.05) "✓" else "⚠ heteroscedastic"))
    if (!is.null(dw)) cat(sprintf("  %s DW: p=%.4f %s\n", nm, dw$p.value,
                                  if (dw$p.value>0.10) "✓" else "⚠ autocorrelated (expected for monthly data)"))
  }
  cat("══════════════════════════════════════════════════════════\n\n")
  invisible(results)
}

## ── F. GAM UTILITIES ─────────────────────────────────────────────────────────

#' Fit full GAM: main smooths + all three tensor-product interactions + seasonality.
fit_gam_full <- function(df,
                         response_var        = "drought_value",
                         k_main              = 8L,
                         k_tensor            = 6L,
                         include_seasonality = TRUE,
                         select              = TRUE) {
  utils_load_packages(c("mgcv","parallel"))
  if (!response_var %in% names(df))
    stop(sprintf("Response '%s' not in data frame.", response_var))
  
  seas_term <- if (include_seasonality && "month" %in% names(df))
    "+ s(month, bs='cc', k=6)" else ""
  
  frm_str <- paste0(
    response_var,
    " ~ s(ONI,k=", k_main, ",bs='cr')",
    " + s(PDO,k=", k_main, ",bs='cr')",
    " + s(PNA,k=", k_main, ",bs='cr')",
    " + ti(ONI,PDO,k=c(", k_tensor, ",", k_tensor, "))",
    " + ti(ONI,PNA,k=c(", k_tensor, ",", k_tensor, "))",
    " + ti(PDO,PNA,k=c(", k_tensor, ",", k_tensor, "))",
    " ", seas_term
  )
  frm <- as.formula(frm_str)
  
  cat(sprintf("  Fitting GAM: %s\n", deparse(frm, width.cutoff=120)[1]))
  fit <- tryCatch(
    mgcv::gam(frm, data=df, method="REML", select=select,
              control=mgcv::gam.control(
                nthreads=max(1L, parallel::detectCores()-1L))),
    error=function(e) { cat("  ❌ GAM failed:", e$message, "\n"); NULL })
  if (is.null(fit)) return(NULL)
  
  s <- summary(fit)
  cat(sprintf("  ✓ GAM fitted | Dev.expl=%.1f%% | adj-R²=%.3f | EDF=%.1f\n",
              s$dev.expl*100, s$r.sq, sum(s$edf)))
  fit
}

#' Partition variance across all GAM smooth terms.
#'
#' PERF: Pass `precomputed_s = summary(gam_fit)` to avoid a redundant
#' summary() call when the caller has already computed it.
gam_variance_partition <- function(fit, precomputed_s = NULL) {
  if (is.null(fit)) return(NULL)
  s   <- if (!is.null(precomputed_s)) precomputed_s else summary(fit)
  smt <- s$s.table
  vp  <- data.frame(
    Term        = rownames(smt),
    EDF         = round(smt[,"edf"], 2),
    F_stat      = round(smt[,"F"], 3),
    p_value     = round(smt[,"p-value"], 4),
    Significant = smt[,"p-value"] < 0.05,
    stringsAsFactors=FALSE)
  vp$Type <- dplyr::case_when(
    grepl("^s\\(ONI\\)",   vp$Term) ~ "Main: ONI",
    grepl("^s\\(PDO\\)",   vp$Term) ~ "Main: PDO",
    grepl("^s\\(PNA\\)",   vp$Term) ~ "Main: PNA",
    grepl("^s\\(month\\)", vp$Term) ~ "Seasonality",
    grepl("ti\\(ONI,PDO",  vp$Term) ~ "Interaction: ONI×PDO  ← key term",
    grepl("ti\\(ONI,PNA",  vp$Term) ~ "Interaction: ONI×PNA",
    grepl("ti\\(PDO,PNA",  vp$Term) ~ "Interaction: PDO×PNA",
    TRUE                            ~ vp$Term)
  cat("\n  ── GAM Variance Partition ──\n")
  print(vp[, c("Type","EDF","F_stat","p_value","Significant")], row.names=FALSE)
  cat(sprintf("  Total deviance explained: %.1f%% | adj-R²: %.3f\n",
              s$dev.expl*100, s$r.sq))
  invisible(vp)
}

#' Concurvity check — GAM analog of multicollinearity/VIF.
gam_concurvity_check <- function(fit) {
  if (is.null(fit)) return(NULL)
  utils_load_packages(c("mgcv"))
  cc  <- mgcv::concurvity(fit, full=TRUE)
  worst_vec <- as.numeric(cc["worst",    ])
  obs_vec   <- as.numeric(cc["observed", ])
  term_nms  <- colnames(cc)
  res <- data.frame(Term=term_nms,
                    Worst=round(worst_vec, 3),
                    Observed=round(obs_vec, 3),
                    Problem=worst_vec > 0.8,
                    stringsAsFactors=FALSE)
  cat("\n  ── GAM Concurvity (analog of VIF) ──\n")
  print(res, row.names=FALSE)
  if (any(res$Problem, na.rm=TRUE))
    cat(sprintf("  ⚠ Concurvity>0.8: %s — term may be redundant\n",
                paste(res$Term[res$Problem & !is.na(res$Problem)], collapse=", ")))
  else cat("  ✓ All concurvity<0.8\n")
  invisible(res)
}

#' Generate a prediction grid on a fine mesh over (x_var, y_var).
predict_gam_surface <- function(fit, x_var, y_var,
                                cond_list=list(), n_grid=80L) {
  if (is.null(fit)) return(NULL)
  utils_load_packages(c("mgcv"))
  x_seq <- seq(range(fit$model[[x_var]], na.rm=TRUE)[1],
               range(fit$model[[x_var]], na.rm=TRUE)[2], length.out=n_grid)
  y_seq <- seq(range(fit$model[[y_var]], na.rm=TRUE)[1],
               range(fit$model[[y_var]], na.rm=TRUE)[2], length.out=n_grid)
  grid  <- expand.grid(setNames(list(x_seq, y_seq), c(x_var, y_var)))
  for (nm in names(cond_list)) grid[[nm]] <- cond_list[[nm]]
  resp_col <- names(fit$model)[1]
  for (mv in names(fit$model)) {
    if (mv == resp_col || mv %in% names(grid)) next
    col_d <- fit$model[[mv]]
    grid[[mv]] <- if (is.numeric(col_d)) median(col_d, na.rm=TRUE) else levels(col_d)[1]
  }
  pred        <- mgcv::predict.gam(fit, newdata=grid, type="response", se.fit=TRUE)
  grid$fit    <- as.numeric(pred$fit)
  grid$se     <- as.numeric(pred$se.fit)
  grid$ci_lo  <- grid$fit - 1.96*grid$se
  grid$ci_hi  <- grid$fit + 1.96*grid$se
  grid
}

#' Primary publication figure: 3-panel ONI×PDO response surfaces conditioned on PNA.
plot_gam_surfaces <- function(fit,
                              index_label = "SPEI-6",
                              lag         = 0L,
                              drought_thr = -1.0,
                              pna_vals    = c("Low PNA (−1.5)"=-1.5,
                                              "Neutral PNA (0)"=0,
                                              "High PNA (+1.5)"=1.5)) {
  if (is.null(fit)) return(NULL)
  utils_load_packages(c("mgcv","ggplot2","dplyr"))
  
  month_cond <- if ("month" %in% names(fit$model)) list(month=1L) else list()
  
  all_grids <- do.call(rbind, lapply(seq_along(pna_vals), function(i) {
    g <- predict_gam_surface(fit, x_var="ONI", y_var="PDO",
                             cond_list=c(list(PNA=pna_vals[i]), month_cond))
    g$PNA_label <- factor(names(pna_vals)[i], levels=names(pna_vals))
    g
  }))
  
  zlim <- max(abs(all_grids$fit), na.rm=TRUE) * 1.05
  
  ggplot2::ggplot(all_grids, ggplot2::aes(x=ONI, y=PDO, fill=fit)) +
    ggplot2::geom_tile() +
    ggplot2::geom_contour(ggplot2::aes(z=fit), breaks=drought_thr,
                          colour="#C0392B", linewidth=1.1) +
    ggplot2::geom_contour(ggplot2::aes(z=fit), breaks=1.0,
                          colour="#1A5276", linewidth=0.7, linetype="dashed") +
    ggplot2::geom_contour(ggplot2::aes(z=fit), breaks=0,
                          colour="grey30", linewidth=0.4, linetype="dotted") +
    ggplot2::geom_hline(yintercept=c(-0.5,0.5), linetype="dotted",
                        colour="grey60", linewidth=0.4) +
    ggplot2::geom_vline(xintercept=c(-0.5,0.5), linetype="dotted",
                        colour="grey60", linewidth=0.4) +
    ggplot2::scale_fill_gradientn(
      colours=c("#1A5276","#2E86C1","#AED6F1","#FDFEFE",
                "#FAD7A0","#E59866","#C0392B"),
      limits=c(-zlim, zlim),
      name=sprintf("E[%s]", index_label),
      guide=ggplot2::guide_colourbar(barwidth=1, barheight=12)) +
    ggplot2::facet_wrap(~PNA_label, ncol=3) +
    shared_ts_theme(11) +
    ggplot2::labs(
      title    = sprintf("GAM Response Surface: E[%s | ONI, PDO] | PNA conditioned", index_label),
      subtitle = sprintf("Lag=%d mo | Red=drought threshold(%.1f) | REML, ti() pure interactions, select=TRUE",
                         lag, drought_thr),
      x="ONI (El Nino -> La Nina axis)", y="PDO Index",
      caption="Blue dashed=wet(+1.0) | Dotted=zero | Quadrant lines=phase boundaries") +
    ggplot2::theme(legend.position="right",
                   strip.text=ggplot2::element_text(face="bold",size=10),
                   panel.spacing=ggplot2::unit(0.8,"lines"))
}

#' Partial effect plots for all smooth terms (gratia if available, else base graphics).
plot_gam_partials <- function(fit, index_label="SPEI-6", lag=0L) {
  if (is.null(fit)) return(NULL)
  if (requireNamespace("gratia", quietly=TRUE)) {
    utils_load_packages(c("gratia","ggplot2","patchwork"))
    p <- gratia::draw(fit, residuals=TRUE, rug=TRUE, scales="fixed", n_contour=8L) +
      patchwork::plot_annotation(
        title=sprintf("GAM Partial Effects: %s | Lag=%d months", index_label, lag),
        subtitle="Each panel: marginal smooth | Points=partial residuals | Shading=95% CI")
    return(p)
  }
  cat("  (gratia not installed — using base mgcv::plot.gam)\n")
  mgcv::plot.gam(fit, pages=1, residuals=TRUE, rug=TRUE, shade=TRUE, shade.col="lightblue",
                 main=sprintf("GAM Partial Effects: %s | Lag=%d", index_label, lag))
  invisible(NULL)
}

#' Four-panel GAM residual diagnostics.
plot_gam_residuals <- function(fit, index_label="SPEI-6") {
  if (is.null(fit)) return(NULL)
  utils_load_packages(c("ggplot2","patchwork"))
  rd <- data.frame(fitted=fitted(fit), residual=residuals(fit, type="deviance"))
  
  p1 <- ggplot2::ggplot(rd, ggplot2::aes(x=fitted, y=residual)) +
    ggplot2::geom_point(alpha=0.3, size=0.8) +
    ggplot2::geom_hline(yintercept=0, linetype="dashed", colour="red") +
    ggplot2::geom_smooth(method="loess", se=FALSE, colour="#2980B9", linewidth=0.8) +
    shared_ts_theme(10) +
    ggplot2::labs(title="Residuals vs Fitted", x="Fitted", y="Deviance residuals")
  
  p2 <- ggplot2::ggplot(rd, ggplot2::aes(sample=residual)) +
    ggplot2::stat_qq(alpha=0.4, size=0.8) + ggplot2::stat_qq_line(colour="red") +
    shared_ts_theme(10) + ggplot2::labs(title="Normal Q-Q", x="Theoretical", y="Sample")
  
  p3 <- ggplot2::ggplot(rd, ggplot2::aes(x=residual)) +
    ggplot2::geom_histogram(bins=35, fill="#2980B9", alpha=0.7, colour="white") +
    shared_ts_theme(10) + ggplot2::labs(title="Residual histogram", x="Residual", y="Count")
  
  acf_obj <- acf(rd$residual, plot=FALSE, lag.max=24)
  acf_df  <- data.frame(lag=as.integer(acf_obj$lag[-1]), acf=as.numeric(acf_obj$acf[-1]))
  ci_val  <- qnorm(0.975) / sqrt(nrow(rd))
  p4 <- ggplot2::ggplot(acf_df, ggplot2::aes(x=lag, y=acf)) +
    ggplot2::geom_hline(yintercept=0) +
    ggplot2::geom_hline(yintercept=c(-ci_val,ci_val),
                        linetype="dashed", colour="blue", linewidth=0.6) +
    ggplot2::geom_segment(ggplot2::aes(xend=lag, yend=0), colour="#E74C3C") +
    ggplot2::geom_point(size=1.5) +
    shared_ts_theme(10) + ggplot2::labs(title="ACF of residuals", x="Lag", y="ACF")
  
  (p1+p2)/(p3+p4) +
    patchwork::plot_annotation(
      title=sprintf("GAM Residual Diagnostics: %s", index_label),
      subtitle="Residual autocorrelation is expected for monthly climate data — monitor but not disqualifying")
}

#' One-row summary of a fitted GAM for the lag-profile comparison table.
#'
#' PERF: Pass `precomputed_s = summary(gam_fit)` to avoid a redundant
#' summary() call when the caller has already computed it.
gam_lag_summary <- function(fit, lag, index_label, precomputed_s = NULL) {
  empty <- data.frame(Index=index_label, Lag=lag, Dev_Expl_pct=NA, AdjR2=NA,
                      EDF_total=NA, p_s_ONI=NA, p_s_PDO=NA, p_s_PNA=NA,
                      p_ti_ONI_PDO=NA, p_ti_ONI_PNA=NA, p_ti_PDO_PNA=NA,
                      GCV=NA, stringsAsFactors=FALSE)
  if (is.null(fit)) return(empty)
  s   <- if (!is.null(precomputed_s)) precomputed_s else summary(fit)
  smt <- s$s.table
  gp  <- function(pat) {
    i <- grep(pat, rownames(smt))
    if (!length(i)) NA_real_ else round(smt[i[1],"p-value"], 4)
  }
  data.frame(Index=index_label, Lag=lag,
             Dev_Expl_pct = round(s$dev.expl*100, 2),
             AdjR2        = round(s$r.sq, 3),
             EDF_total    = round(sum(s$edf), 1),
             p_s_ONI      = gp("^s\\(ONI"),
             p_s_PDO      = gp("^s\\(PDO"),
             p_s_PNA      = gp("^s\\(PNA"),
             p_ti_ONI_PDO = gp("ti\\(ONI,PDO"),
             p_ti_ONI_PNA = gp("ti\\(ONI,PNA"),
             p_ti_PDO_PNA = gp("ti\\(PDO,PNA"),
             GCV          = round(fit$gcv.ubre, 4),
             stringsAsFactors=FALSE)
}

cat("✓ utils_teleconnection_addon.R loaded (MM + GAM utilities)\n")