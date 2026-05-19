####################################################################################
# w8_6_regime_causal_test.R
# REGIME-DEPENDENT CAUSAL TEST — Structural break in ENSO drought coefficient
# Kretschmer et al. (2021, Part V) / Saggioro et al. (2020)
#
# Tests whether the causal effect of ENSO on Nechako drought has shifted between
# the early and late periods of the record, linking to:
#   - PELT changepoints from w1_trend_test.R
#   - F_thm thermodynamic fraction from w11
#
# Two formal tests are applied:
#   1. OLS-CUSUM structural break test (strucchange::efp)
#   2. Chow test at the PELT-detected changepoint year
#
# Usage:
#   source("w8_0_setup_moderated_mediation.R")
#   source("w8_6_regime_causal_test.R")
#
# Outputs written to MM_DIR/regime_causal_tests/:
#   regime_causal_<index>.pdf      -- split-period ONI coefficient plot
#   structural_break_summary.csv   -- CUSUM + Chow results, all indices
#   regime_causal_all.xlsx         -- Excel workbook
####################################################################################

source("w8_0_setup_moderated_mediation.R")
utils_load_packages(c("strucchange","ggplot2","patchwork","openxlsx","dplyr"))

RC_START <- proc.time()
cat("\n================================================================\n")
cat("  w8_6  REGIME-DEPENDENT CAUSAL TEST\n")
cat("================================================================\n\n")

REGIME_DIR <- file.path(MM_DIR, "regime_causal_tests")
dir.create(REGIME_DIR, showWarnings=FALSE, recursive=TRUE)

####################################################################################
# SPLIT YEAR — read from w1 PELT output or fall back to 1990
####################################################################################

pelt_path <- file.path(WD_PATH, "trend_analysis", "pelt_changepoints.csv")
SPLIT_YEAR <- if (file.exists(pelt_path)) {
  cp_df <- tryCatch(read.csv(pelt_path, stringsAsFactors=FALSE), error=function(e) NULL)
  if (!is.null(cp_df) && "changepoint_year" %in% names(cp_df) && nrow(cp_df) > 0) {
    yr <- as.integer(cp_df$changepoint_year[1])
    cat(sprintf("  PELT changepoint loaded: split year = %d\n", yr)); yr
  } else {
    cat("  PELT file empty or malformed; using default split year = 1990\n"); 1990L
  }
} else {
  cat("  PELT file not found (expected: trend_analysis/pelt_changepoints.csv)\n")
  cat("  Using default split year = 1990\n"); 1990L
}
cat(sprintf("  Early period: 1950 – %d  |  Late period: %d – 2024\n\n",
            SPLIT_YEAR - 1L, SPLIT_YEAR))

####################################################################################
# HELPER — fit and summarise the causal regression for one period
####################################################################################

.fit_period <- function(df_period, period_label) {
  if (nrow(df_period) < 10) {
    cat(sprintf("    ! %s: too few obs (%d) to fit\n", period_label, nrow(df_period)))
    return(NULL)
  }
  tryCatch(
    lm(Y_c ~ ONI_c + PDO_c + F_thm_std, data=df_period),
    error=function(e) {
      cat(sprintf("    ! %s lm failed: %s\n", period_label, e$message)); NULL
    })
}

####################################################################################
# MAIN LOOP OVER ALL DROUGHT INDICES
####################################################################################

all_break_rows <- list()

for (idx_spec in ALL_INDICES) {

  lbl <- sprintf("%s-%d", toupper(idx_spec$index), idx_spec$scale)
  cat(sprintf("  -- %s --\n", lbl))

  # Use optimal lag from GAM lag profile if available; else default to 2
  gam_lp_path <- file.path(MM_DIR, lbl, "gam_lag_profile.csv")
  opt_lag <- if (file.exists(gam_lp_path)) {
    glp <- tryCatch(read.csv(gam_lp_path, stringsAsFactors=FALSE), error=function(e) NULL)
    if (!is.null(glp) && "Dev_Expl_pct" %in% names(glp))
      glp$Lag[which.max(glp$Dev_Expl_pct)] else 2L
  } else 2L

  df_full <- tryCatch(
    build_analysis_dataframe(idx_spec$index, idx_spec$scale,
                             lag_months=opt_lag, start_year=1950L, end_year=2024L,
                             tele_df=TELE_MASTER),
    error=function(e) { cat(sprintf("    XX Data: %s\n", e$message)); NULL })
  if (is.null(df_full)) next

  # Merge F_thm thermodynamic confounder
  if (!is.null(FTHM_SERIES) && "F_thm" %in% names(FTHM_SERIES)) {
    df_full <- merge(df_full, FTHM_SERIES, by="date", all.x=TRUE)
    df_full$F_thm_std <- as.numeric(scale(df_full$F_thm))
  } else {
    df_full$F_thm_std <- 0   # zeroes out F_thm if w11 output is unavailable
    cat("    Note: F_thm unavailable; F_thm_std set to 0\n")
  }

  df_full$year <- as.integer(format(df_full$date, "%Y"))
  df_early     <- df_full[df_full$year <  SPLIT_YEAR, ]
  df_late      <- df_full[df_full$year >= SPLIT_YEAR, ]

  # ── Regression for each period ─────────────────────────────────────────────
  fit_full  <- .fit_period(df_full,  "Full")
  fit_early <- .fit_period(df_early, sprintf("Early (< %d)",  SPLIT_YEAR))
  fit_late  <- .fit_period(df_late,  sprintf("Late  (>=%d)",  SPLIT_YEAR))
  if (is.null(fit_full)) next

  .oni_stats <- function(fit) {
    if (is.null(fit)) return(rep(NA_real_, 4L))
    s <- tryCatch(coef(summary(fit))["ONI_c", ], error=function(e) rep(NA_real_,4))
    s
  }
  full_oni  <- .oni_stats(fit_full)
  early_oni <- .oni_stats(fit_early)
  late_oni  <- .oni_stats(fit_late)

  cat(sprintf("    Full  (n=%d):  ONI = %6.4f  p = %.4f\n",
              nrow(df_full),  full_oni[1],  full_oni[4]))
  cat(sprintf("    Early (n=%d):  ONI = %6.4f  p = %.4f\n",
              nrow(df_early), early_oni[1], early_oni[4]))
  cat(sprintf("    Late  (n=%d):  ONI = %6.4f  p = %.4f\n\n",
              nrow(df_late),  late_oni[1],  late_oni[4]))

  # ── CUSUM structural break test ────────────────────────────────────────────
  cusum_stat <- NA_real_; cusum_p <- NA_real_
  tryCatch({
    efp_test   <- strucchange::efp(Y_c ~ ONI_c + PDO_c + F_thm_std,
                                    data=df_full, type="OLS-CUSUM")
    sc_test    <- strucchange::sctest(efp_test)
    cusum_stat <- round(sc_test$statistic, 4)
    cusum_p    <- round(sc_test$p.value,   4)
    cat(sprintf("    CUSUM:  stat = %.4f  p = %.4f%s\n",
                cusum_stat, cusum_p, if (cusum_p < 0.05) "  ** structural break" else ""))
  }, error=function(e) cat(sprintf("    XX CUSUM: %s\n", e$message)))

  # ── Chow test at SPLIT_YEAR ────────────────────────────────────────────────
  chow_f <- NA_real_; chow_p <- NA_real_
  tryCatch({
    k <- 4L  # intercept + ONI + PDO + F_thm
    n_e <- nrow(df_early); n_l <- nrow(df_late)
    if (!is.null(fit_early) && !is.null(fit_late) &&
        n_e >= k + 2L && n_l >= k + 2L) {
      rss_r <- sum(resid(fit_full)^2)
      rss_u <- sum(resid(fit_early)^2) + sum(resid(fit_late)^2)
      chow_f <- round(((rss_r - rss_u) / k) / (rss_u / (n_e + n_l - 2*k)), 4)
      chow_p <- round(pf(chow_f, k, n_e + n_l - 2*k, lower.tail=FALSE), 4)
      cat(sprintf("    Chow:   F = %.4f  p = %.4f%s\n",
                  chow_f, chow_p, if (chow_p < 0.05) "  ** structural break" else ""))
    } else cat("    Chow: insufficient obs in one sub-period; skipped\n")
  }, error=function(e) cat(sprintf("    XX Chow: %s\n", e$message)))
  cat("\n")

  # ── Coefficient comparison plot ────────────────────────────────────────────
  coef_df <- data.frame(
    Period  = factor(
      c(sprintf("Early\n(< %d)",  SPLIT_YEAR),
        sprintf("Late\n(>= %d)", SPLIT_YEAR),
        "Full\nperiod"),
      levels=c(sprintf("Early\n(< %d)", SPLIT_YEAR),
               sprintf("Late\n(>= %d)", SPLIT_YEAR),
               "Full\nperiod")),
    Coeff   = c(early_oni[1], late_oni[1], full_oni[1]),
    SE      = c(early_oni[2], late_oni[2], full_oni[2]),
    stringsAsFactors=FALSE)
  coef_df$lower <- coef_df$Coeff - 1.96 * coef_df$SE
  coef_df$upper <- coef_df$Coeff + 1.96 * coef_df$SE
  coef_df$col   <- c("#2980B9","#C0392B","#27AE60")  # early/late/full colours

  p_coef <- ggplot2::ggplot(coef_df,
      ggplot2::aes(x=Period, y=Coeff, colour=Period)) +
    ggplot2::geom_hline(yintercept=0, linetype="dashed", colour="grey50") +
    ggplot2::geom_errorbar(ggplot2::aes(ymin=lower, ymax=upper),
                           width=0.15, linewidth=1.1) +
    ggplot2::geom_point(size=4) +
    ggplot2::scale_colour_manual(values=setNames(coef_df$col, coef_df$Period)) +
    ggplot2::guides(colour="none") +
    shared_ts_theme(12) +
    ggplot2::labs(
      title    = sprintf("Regime-dependent ENSO causal coefficient | %s", lbl),
      subtitle = sprintf(
        "CUSUM p=%.3f | Chow p=%.3f | split=%d | lag=%d months | 95%% CI shown",
        cusum_p, if (is.na(chow_p)) NaN else chow_p,
        SPLIT_YEAR, opt_lag),
      x="Period",
      y="ONI coefficient (Y_c ~ ONI + PDO + F_thm)")

  reg_pdf <- file.path(REGIME_DIR, sprintf("regime_causal_%s.pdf", lbl))
  if (safe_pdf(reg_pdf, width=8, height=6)) {
    print(p_coef); grDevices::dev.off()
    cat(sprintf("  Saved: regime_causal_%s.pdf\n", lbl))
  }

  # ── Collect summary row ────────────────────────────────────────────────────
  all_break_rows[[lbl]] <- data.frame(
    Index           = lbl,
    OptimalLag      = opt_lag,
    Split_Year      = SPLIT_YEAR,
    N_full          = nrow(df_full),
    N_early         = nrow(df_early),
    N_late          = nrow(df_late),
    Full_ONI_coeff  = round(full_oni[1],  4),
    Full_ONI_p      = round(full_oni[4],  4),
    Early_ONI_coeff = round(early_oni[1], 4),
    Early_ONI_p     = round(early_oni[4], 4),
    Late_ONI_coeff  = round(late_oni[1],  4),
    Late_ONI_p      = round(late_oni[4],  4),
    CUSUM_stat      = cusum_stat,
    CUSUM_p         = cusum_p,
    Chow_F          = chow_f,
    Chow_p          = chow_p,
    Structural_break = (!is.na(cusum_p) && cusum_p < 0.05) |
                       (!is.na(chow_p)  && chow_p  < 0.05),
    stringsAsFactors=FALSE)
}

####################################################################################
# COLLATE AND SAVE
####################################################################################

if (length(all_break_rows) == 0) {
  cat("  No results produced.  Check data availability.\n\n")
  stop("w8_6: no output produced.", call.=FALSE)
}

break_df <- do.call(rbind, all_break_rows)
rownames(break_df) <- NULL

out_csv <- file.path(MM_DIR, "structural_break_summary.csv")
write.csv(break_df, out_csv, row.names=FALSE)
cat(sprintf("\n  Saved: %s\n\n", out_csv))

# Print summary table
cat("  Structural break summary:\n")
print(break_df[, c("Index","OptimalLag","Full_ONI_coeff","Early_ONI_coeff",
                   "Late_ONI_coeff","CUSUM_p","Chow_p","Structural_break")],
      row.names=FALSE)

# Excel workbook
wb  <- openxlsx::createWorkbook()
hdr <- openxlsx::createStyle(textDecoration="bold", fgFill="#1B3A6B", fontColour="#FFFFFF")
.add <- function(wb, name, df) {
  name <- substr(name, 1, 31)
  openxlsx::addWorksheet(wb, name)
  openxlsx::writeData(wb, name, df)
  openxlsx::addStyle(wb, name, hdr, rows=1, cols=seq_len(ncol(df)))
}
.add(wb, "Structural_Break_Summary", break_df)
.add(wb, "Metadata", data.frame(
  Parameter=c("Split_year","PELT_file_found","F_thm_used","Test_CUSUM","Test_Chow",
               "Regression_formula","Run_date"),
  Value=c(as.character(SPLIT_YEAR), as.character(file.exists(pelt_path)),
          as.character(!is.null(FTHM_SERIES)),
          "OLS-CUSUM via strucchange::efp()",
          "Chow F-test at SPLIT_YEAR",
          "Y_c ~ ONI_c + PDO_c + F_thm_std",
          format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  stringsAsFactors=FALSE))
out_xlsx <- file.path(MM_DIR, "regime_causal_all.xlsx")
openxlsx::saveWorkbook(wb, out_xlsx, overwrite=TRUE)
cat(sprintf("  Saved: %s\n", out_xlsx))

cat("\n================================================================\n")
cat("  w8_6 REGIME CAUSAL TEST COMPLETE\n")
cat("================================================================\n")
cat(sprintf("  Time    : %s\n",  fmt_dur(elapsed_sec(RC_START))))
cat(sprintf("  Split   : %d\n",  SPLIT_YEAR))
sig_n <- sum(break_df$Structural_break, na.rm=TRUE)
cat(sprintf("  Indices with structural break (p<0.05): %d / %d\n\n",
            sig_n, nrow(break_df)))
cat("  Interpretation:\n")
cat("    Structural break = ENSO causal effect has shifted between early/late periods.\n")
cat("    Cross-reference with PDO phase shift and F_thm trend from w11.\n\n")
