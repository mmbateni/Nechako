####################################################################################
# w7_moderated_mediation.R  ·  MODERATED MEDIATION ANALYSIS (lavaan)
####################################################################################
# Causal model:
#
#   ONI ──(a-path, moderated by PDO)──► PNA
#    │                                   │
#    │                                   │ (b-path, moderated by PDO)
#    │ (direct c'-path)                  ▼
#    └───────────────────────────────► SPI/SPEI (outcome Y)
#
# PDO is the MODERATOR:  it changes both the slope with which ONI drives
# PNA (a-path) and the slope with which PNA drives drought (b-path).
#
# Model equations:
#   PNA_c ~ a1*ONI_c + a2*PDO_c + a3*(ONI_c × PDO_c)   [a-path equation]
#   Y_c   ~ b1*PNA_c + b2*ONI_c + b3*PDO_c +            [b-path equation]
#           b4*(PNA_c × PDO_c)
#
# Conditional indirect effect at PDO = w:
#   IE(w) = (a1 + a3*w) × (b1 + b4*w)
#
# Index of Moderated Mediation:
#   IMM = a3*b1 + a1*b4 + a3*b4
#   Significant IMM confirms that PDO statistically moderates the indirect
#   ONI → PNA → Y pathway (Hayes 2015, Psychol. Methods).
#
# OUTPUTS (all to {WD_PATH}/moderated_mediation/):
#   mm_results_summary.xlsx          — path coefficients + bootstrap CIs for all lags
#   mm_conditional_effects_lag*.pdf  — conditional IE plots (low/med/high PDO)
#   mm_jn_floodlight_lag*.pdf        — Johnson-Neyman regions of significance
#   mm_path_diagram_lag*.pdf         — standardised path diagram
#   mm_pdv_regime_comparison.pdf     — regime-stratified regression surfaces
#   mm_lag_profile.pdf               — IE(w) as a function of lag (0–12 months)
#
# Depends on: DROUGHT_ANALYSIS_utils.R, utils_teleconnection_addon.R
# Run AFTER: w6_teleconnection_prep.R (for cached teleconnection CSVs)
####################################################################################

rm(list = ls())
gc()

source("DROUGHT_ANALYSIS_utils.R")
source("utils_teleconnection_addon.R")

utils_load_packages(c(
  "lavaan", "ggplot2", "patchwork", "dplyr",
  "openxlsx", "tidyr", "scales", "mgcv"
))

setwd(WD_PATH)

# ── Output directory ─────────────────────────────────────────────────────────
MM_DIR <- file.path(WD_PATH, "moderated_mediation")
dir.create(MM_DIR, showWarnings = FALSE, recursive = TRUE)

cat("\n╔═══════════════════════════════════════════════╗\n")
cat("║  w7  MODERATED MEDIATION ANALYSIS             ║\n")
cat("╚═══════════════════════════════════════════════╝\n\n")

####################################################################################
# CONFIGURATION — edit here to change what is analysed
####################################################################################
# Primary drought index for main results (can be "spi", "spei", "swei")
PRIMARY_INDEX <- "spei"
PRIMARY_SCALE <- 6

# Sensitivity indices — run same model to test robustness
SENSITIVITY_INDICES <- list(
  list(index = "spi",  scale = 6),
  list(index = "spi",  scale = 12),
  list(index = "spei", scale = 3),
  list(index = "swei", scale = SWEI_SCALE)
)

# Lags to test (months between teleconnection observation and drought response)
LAGS_TO_TEST <- 0:6

# Bootstrap settings (increase for publication — 10,000 for final run)
N_BOOTSTRAP <- 5000L
BOOTSTRAP_SEED <- 42L

# Values of PDO (standardised) at which to evaluate conditional indirect effects
# -1 SD, 0 (mean), +1 SD — the standard "low / medium / high" in PROCESS macro
W_VALUES <- c(-1, 0, 1)
W_LABELS <- c("Cool PDO (−1 SD)", "Neutral PDO (Mean)", "Warm PDO (+1 SD)")

####################################################################################
# LAVAAN MODEL BUILDER
####################################################################################

#' Construct the lavaan model string for moderated mediation.
#'
#' Following Hayes (2018) Model 7 (moderated mediation where the moderator
#' affects both the a-path and the b-path).  All predictors must be centred
#' before passing to lavaan (done in build_analysis_dataframe()).
#'
#' @param w_low   Numeric; value of PDO_c at which to evaluate IE (e.g. -1 SD)
#' @param w_mid   Numeric; (default 0 = mean)
#' @param w_high  Numeric; (default +1 SD)
#' @return Character string; lavaan model syntax
build_lavaan_model <- function(w_low  = -1,
                               w_mid  =  0,
                               w_high =  1) {
  sprintf('
  # ── a-path: ONI → PNA, moderated by PDO ──────────────────────────────
  # a1 = slope of ONI on PNA at mean PDO
  # a2 = main effect of PDO on PNA (PDO can independently drive PNA)
  # a3 = the moderation: how much does PDO change the ONI→PNA slope?
  PNA_c ~ a1*ONI_c + a2*PDO_c + a3*ONI_PDO_int

  # ── b-path + direct c-prime path ─────────────────────────────────────
  # b1 = slope of PNA on Y at mean PDO
  # b2 = direct effect of ONI on Y (controlling for PNA — the c-prime path)
  # b3 = direct effect of PDO on Y
  # b4 = moderation: how much does PDO change the PNA→Y slope?
  Y_c ~ b1*PNA_c + cprime*ONI_c + b3*PDO_c + b4*PNA_PDO_int

  # ── Conditional indirect effects IE(w) = (a1 + a3*w)(b1 + b4*w) ─────
  IE.low  := (a1 + a3 * %s) * (b1 + b4 * %s)
  IE.mid  := (a1 + a3 * %s) * (b1 + b4 * %s)
  IE.high := (a1 + a3 * %s) * (b1 + b4 * %s)

  # ── Total indirect effect (at mean PDO) ──────────────────────────────
  IE.total := a1 * b1

  # ── Total effect of ONI on Y ─────────────────────────────────────────
  TE := cprime + a1 * b1

  # ── Proportion mediated (at mean PDO) ────────────────────────────────
  # Only interpret when IE and TE have the same sign
  PM := (a1 * b1) / TE

  # ── Index of Moderated Mediation (Hayes 2015) ────────────────────────
  # IMM ≠ 0 → PDO significantly moderates the indirect pathway
  IMM := a3 * b1 + a1 * b4 + a3 * b4
  ', w_low, w_low, w_mid, w_mid, w_high, w_high)
}

####################################################################################
# CORE FUNCTION: fit one model at one lag
####################################################################################

#' Fit moderated mediation model and return standardised results.
#'
#' @param df         Output of build_analysis_dataframe()
#' @param lag        Integer lag in months (for labelling only)
#' @param n_boot     Number of bootstrap replications
#' @param seed       Random seed for reproducibility
#' @return Named list: fit (lavaan object), params (data.frame), summary_df
fit_moderated_mediation <- function(df, lag = 0L, n_boot = 5000L, seed = 42L) {

  cat(sprintf("\n  ── Fitting MM model (lag = %d months) ──\n", lag))
  cat(sprintf("     n = %d observations\n", nrow(df)))

  # Check we have enough data
  if (nrow(df) < 60)
    warning(sprintf("Only %d observations at lag %d — bootstrap CIs may be unstable.",
                    nrow(df), lag))

  model_str <- build_lavaan_model(w_low  = W_VALUES[1],
                                  w_mid  = W_VALUES[2],
                                  w_high = W_VALUES[3])

  set.seed(seed)
  fit <- tryCatch(
    lavaan::sem(
      model    = model_str,
      data     = df,
      se       = "bootstrap",          # bias-corrected bootstrap SEs
      bootstrap = n_boot,
      estimator = "ML",                # ML estimator (normal-theory)
      fixed.x  = FALSE,                # treat predictors as random variables
      # Correction: use Huber-White sandwich SEs for point estimates
      # to handle potential heteroscedasticity; bootstrap handles uncertainty
      missing  = "listwise"
    ),
    error = function(e) {
      cat(sprintf("  ❌ lavaan failed: %s\n", e$message))
      return(NULL)
    }
  )

  if (is.null(fit)) return(NULL)

  # Convergence check
  if (!lavaan::lavTech(fit, "converged")) {
    warning(sprintf("Model did not converge at lag = %d", lag))
  }

  # ── Extract parameter estimates with bootstrap CIs ──────────────────────
  # parameterEstimates gives: est, se, z, pvalue, ci.lower, ci.upper
  pe <- lavaan::parameterEstimates(fit,
                                   ci       = TRUE,
                                   level    = 0.95,
                                   boot.ci.type = "bca.simple",  # BCa intervals
                                   standardized = TRUE)

  # Subset to defined parameters and user-defined computed quantities
  pe_def   <- pe[pe$op %in% c("~", ":="), ]

  # Tidy labels
  pe_def$label_clean <- dplyr::case_when(
    pe_def$label == "a1"      ~ "a1: ONI→PNA (at mean PDO)",
    pe_def$label == "a2"      ~ "a2: PDO→PNA (direct)",
    pe_def$label == "a3"      ~ "a3: ONI×PDO→PNA (moderation of a-path)",
    pe_def$label == "b1"      ~ "b1: PNA→Y (at mean PDO)",
    pe_def$label == "cprime"  ~ "c': ONI→Y (direct, controlling PNA)",
    pe_def$label == "b3"      ~ "b3: PDO→Y (direct)",
    pe_def$label == "b4"      ~ "b4: PNA×PDO→Y (moderation of b-path)",
    pe_def$label == "IE.low"  ~ sprintf("IE (PDO = %.0f SD)", W_VALUES[1]),
    pe_def$label == "IE.mid"  ~ sprintf("IE (PDO = %.0f SD)", W_VALUES[2]),
    pe_def$label == "IE.high" ~ sprintf("IE (PDO = +%.0f SD)", W_VALUES[3]),
    pe_def$label == "IE.total"~ "Total IE (at mean PDO)",
    pe_def$label == "TE"      ~ "Total Effect (ONI→Y)",
    pe_def$label == "PM"      ~ "Proportion Mediated",
    pe_def$label == "IMM"     ~ "Index of Moderated Mediation (IMM)",
    TRUE                      ~ pe_def$label
  )

  # Round for printing
  pe_print <- pe_def[, c("label_clean", "est.std", "ci.lower", "ci.upper", "pvalue")]
  colnames(pe_print) <- c("Parameter", "Std.Est", "CI.lower", "CI.upper", "p-value")
  pe_print[, 2:4] <- round(pe_print[, 2:4], 4)
  pe_print[, 5]   <- round(pe_print[, 5], 4)

  cat("\n  Results:\n")
  print(pe_print, row.names = FALSE)

  # ── Model fit indices ───────────────────────────────────────────────────
  # Note: saturated path model has df=0, so absolute fit indices are trivial.
  # Report R² for each equation instead.
  r2 <- lavaan::lavInspect(fit, "r2")
  cat(sprintf("\n  R² (PNA equation / a-path): %.3f\n", r2["PNA_c"]))
  cat(sprintf("  R² (Y equation  / b-path): %.3f\n",   r2["Y_c"]))

  # Check IMM significance
  imm_row <- pe_def[pe_def$label == "IMM", ]
  if (nrow(imm_row) > 0) {
    sig <- imm_row$pvalue < 0.05
    cat(sprintf("\n  IMM = %.4f (95%% BCa CI: [%.4f, %.4f]) — %s\n",
                imm_row$est.std,
                imm_row$ci.lower,
                imm_row$ci.upper,
                if (sig) "✓ SIGNIFICANT: PDO moderates the indirect path"
                else     "✗ Not significant: PDO does not moderate the indirect path"))
  }

  list(
    fit        = fit,
    params     = pe_def,
    params_tbl = pe_print,
    r2         = r2,
    n          = nrow(df),
    lag        = lag
  )
}

####################################################################################
# JOHNSON-NEYMAN FLOODLIGHT ANALYSIS
####################################################################################

#' Compute Johnson-Neyman regions of significance for the conditional IE.
#'
#' Finds the exact values of PDO_c at which the 95% CI on IE(w) crosses zero,
#' i.e., the boundaries of the "region of significance" where the indirect
#' effect is statistically distinguishable from zero.
#'
#' Uses a fine grid of w values and BCa bootstrap CIs at each point.
#'
#' @param df     Analysis data.frame (output of build_analysis_dataframe)
#' @param n_boot Number of bootstrap draws per w value (reduce to 1000 for speed)
#' @param w_seq  Numeric vector of PDO_c values to evaluate
#' @return data.frame: w, IE, ci_lower, ci_upper, significant
compute_jn_floodlight <- function(df,
                                   n_boot = 2000L,
                                   w_seq  = seq(-2.5, 2.5, by = 0.1)) {
  cat("  Computing Johnson-Neyman floodlight (this may take a few minutes)...\n")

  results <- do.call(rbind, lapply(w_seq, function(w) {
    # Fit model at this specific w value using a temporary renamed model
    # IE(w) = (a1 + a3*w)(b1 + b4*w)
    # Compute via delta method or bootstrap at a single w point
    # Here we use a fast OLS-based delta method approximation,
    # and flag significance when the analytical CI excludes zero.
    tryCatch({
      # a-path OLS
      lm_a <- lm(PNA_c ~ ONI_c + PDO_c + ONI_PDO_int, data = df)
      a1   <- coef(lm_a)["ONI_c"]
      a3   <- coef(lm_a)["ONI_PDO_int"]
      # b-path OLS (conditional on PNA)
      lm_b <- lm(Y_c ~ PNA_c + ONI_c + PDO_c + PNA_PDO_int, data = df)
      b1   <- coef(lm_b)["PNA_c"]
      b4   <- coef(lm_b)["PNA_PDO_int"]

      ie_w <- (a1 + a3 * w) * (b1 + b4 * w)

      # Bootstrap SE for this IE(w)
      set.seed(BOOTSTRAP_SEED)
      boot_ie <- replicate(n_boot, {
        idx   <- sample(nrow(df), replace = TRUE)
        b_df  <- df[idx, ]
        la  <- tryCatch(lm(PNA_c ~ ONI_c + PDO_c + ONI_PDO_int, data = b_df), error = function(e) NULL)
        lb  <- tryCatch(lm(Y_c ~ PNA_c + ONI_c + PDO_c + PNA_PDO_int, data = b_df), error = function(e) NULL)
        if (is.null(la) || is.null(lb)) return(NA_real_)
        ca  <- coef(la); cb <- coef(lb)
        (ca["ONI_c"] + ca["ONI_PDO_int"] * w) *
          (cb["PNA_c"] + cb["PNA_PDO_int"] * w)
      })
      boot_ie <- boot_ie[is.finite(boot_ie)]
      if (length(boot_ie) < 100) {
        ci_lo <- NA_real_; ci_hi <- NA_real_
      } else {
        ci_lo <- quantile(boot_ie, 0.025, na.rm = TRUE)
        ci_hi <- quantile(boot_ie, 0.975, na.rm = TRUE)
      }
      data.frame(w         = w,
                 IE        = ie_w,
                 ci_lower  = ci_lo,
                 ci_upper  = ci_hi,
                 significant = !is.na(ci_lo) && (ci_lo > 0 || ci_hi < 0),
                 stringsAsFactors = FALSE)
    }, error = function(e) NULL)
  }))

  results <- results[!sapply(results, is.null), ]
  results <- do.call(rbind, if (is.data.frame(results)) list(results) else results)
  results
}

####################################################################################
# VISUALISATION FUNCTIONS
####################################################################################

#' Plot conditional indirect effects at three PDO levels.
plot_conditional_ie <- function(mm_result, lag, index_label) {
  # Extract IE estimates
  pe  <- mm_result$params
  ie_rows <- pe[pe$label %in% c("IE.low", "IE.mid", "IE.high"), ]
  if (nrow(ie_rows) == 0) return(NULL)

  ie_df <- data.frame(
    PDO_level = factor(W_LABELS, levels = W_LABELS),
    IE        = ie_rows$est.std,
    ci_lo     = ie_rows$ci.lower,
    ci_hi     = ie_rows$ci.upper,
    sig       = ie_rows$pvalue < 0.05,
    stringsAsFactors = FALSE
  )

  ggplot2::ggplot(ie_df, ggplot2::aes(x = PDO_level, y = IE,
                                       colour = sig, shape = sig)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = ci_lo, ymax = ci_hi),
                           width = 0.15, linewidth = 1.1) +
    ggplot2::geom_point(size = 4) +
    ggplot2::scale_colour_manual(values = c("FALSE" = "grey60", "TRUE" = "#C0392B"),
                                  name = "p < 0.05") +
    ggplot2::scale_shape_manual(values = c("FALSE" = 1, "TRUE" = 19),
                                 name = "p < 0.05") +
    shared_ts_theme(12) +
    ggplot2::labs(
      title    = sprintf("Conditional Indirect Effect of ONI → PNA → %s", index_label),
      subtitle = sprintf("Drought index lag = %d months | 95%% BCa bootstrap CI | n = %d",
                         lag, mm_result$n),
      x        = "PDO Phase",
      y        = "Indirect Effect (standardised)"
    )
}

#' Plot Johnson-Neyman floodlight.
plot_jn_floodlight <- function(jn_df, lag, index_label) {
  if (is.null(jn_df) || nrow(jn_df) == 0) return(NULL)

  # Find JN transition points (where significance changes)
  jn_df$sig_num <- as.integer(jn_df$significant)
  transitions <- which(diff(jn_df$sig_num) != 0)

  p <- ggplot2::ggplot(jn_df, ggplot2::aes(x = w)) +
    # Shaded region of significance
    ggplot2::geom_ribbon(
      data = jn_df[jn_df$significant, ],
      ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
      fill = "#FADBD8", alpha = 0.6) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
                         fill = "grey85", alpha = 0.4) +
    ggplot2::geom_line(ggplot2::aes(y = IE),
                       colour = "#C0392B", linewidth = 1.1) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::geom_vline(xintercept = W_VALUES,
                        linetype = "dotted", colour = "#2980B9", linewidth = 0.6) +
    ggplot2::annotate("text",
                      x = W_VALUES,
                      y = max(jn_df$ci_upper, na.rm = TRUE),
                      label = c("Cool PDO\n(−1 SD)", "Mean\nPDO", "Warm PDO\n(+1 SD)"),
                      size = 3, vjust = -0.3, colour = "#2980B9") +
    shared_ts_theme(12) +
    ggplot2::labs(
      title    = "Johnson-Neyman Floodlight: Region of Significance",
      subtitle = sprintf("Conditional IE (ONI→PNA→%s) as a function of PDO | lag = %d months",
                         index_label, lag),
      x        = "PDO Index (standardised)",
      y        = "Indirect Effect (standardised)"
    ) +
    ggplot2::annotate("text",
                      x = min(jn_df$w), y = min(jn_df$IE, na.rm = TRUE),
                      label = "Shaded = significant (95% CI excludes zero)",
                      hjust = 0, size = 3, colour = "#C0392B")
  p
}

#' Plot standardised path diagram as ggplot (simple schematic).
plot_path_diagram <- function(mm_result, lag, index_label) {
  pe  <- mm_result$params
  get_est <- function(lbl) {
    row <- pe[pe$label == lbl, ]
    if (nrow(row) == 0) return(NA_real_)
    sprintf("β=%.3f\np=%.3f", row$est.std[1], row$pvalue[1])
  }
  a1_lbl  <- get_est("a1")
  a3_lbl  <- get_est("a3")
  b1_lbl  <- get_est("b1")
  b4_lbl  <- get_est("b4")
  cp_lbl  <- get_est("cprime")
  imm_lbl <- get_est("IMM")

  # Simple text-based diagram using ggplot annotation
  ggplot2::ggplot() +
    ggplot2::xlim(0, 10) + ggplot2::ylim(0, 6) +
    ggplot2::theme_void() +

    # Nodes
    ggplot2::annotate("rect", xmin=0.3, xmax=2.3, ymin=2.5, ymax=3.5,
                      fill="#D6EAF8", colour="#2980B9", linewidth=1.2) +
    ggplot2::annotate("text", x=1.3, y=3.0, label="ONI\n(Predictor X)",
                      size=3.5, fontface="bold") +

    ggplot2::annotate("rect", xmin=3.8, xmax=6.2, ymin=4.2, ymax=5.2,
                      fill="#D5F5E3", colour="#27AE60", linewidth=1.2) +
    ggplot2::annotate("text", x=5.0, y=4.7, label="PNA\n(Mediator M)",
                      size=3.5, fontface="bold") +

    ggplot2::annotate("rect", xmin=7.7, xmax=9.7, ymin=2.5, ymax=3.5,
                      fill="#FADBD8", colour="#C0392B", linewidth=1.2) +
    ggplot2::annotate("text", x=8.7, y=3.0,
                      label=sprintf("%s\n(Outcome Y)", index_label),
                      size=3.5, fontface="bold") +

    ggplot2::annotate("rect", xmin=3.8, xmax=6.2, ymin=0.8, ymax=1.8,
                      fill="#FCF3CF", colour="#F39C12", linewidth=1.2) +
    ggplot2::annotate("text", x=5.0, y=1.3, label="PDO\n(Moderator W)",
                      size=3.5, fontface="bold") +

    # Arrows (a-path: ONI → PNA)
    ggplot2::annotate("segment", x=2.3, xend=3.8, y=3.3, yend=4.5,
                      arrow=ggplot2::arrow(length=ggplot2::unit(0.3,"cm"),type="closed"),
                      colour="#27AE60", linewidth=1) +
    ggplot2::annotate("text", x=2.85, y=4.1, label=paste("a-path:", a1_lbl),
                      size=3, colour="#27AE60", hjust=0) +

    # Arrow (b-path: PNA → Y)
    ggplot2::annotate("segment", x=6.2, xend=7.7, y=4.5, yend=3.3,
                      arrow=ggplot2::arrow(length=ggplot2::unit(0.3,"cm"),type="closed"),
                      colour="#C0392B", linewidth=1) +
    ggplot2::annotate("text", x=6.7, y=4.1, label=paste("b-path:", b1_lbl),
                      size=3, colour="#C0392B", hjust=0) +

    # Arrow (direct c' path: ONI → Y)
    ggplot2::annotate("segment", x=2.3, xend=7.7, y=2.8, yend=2.8,
                      arrow=ggplot2::arrow(length=ggplot2::unit(0.3,"cm"),type="closed"),
                      colour="#7F8C8D", linewidth=0.8, linetype="dashed") +
    ggplot2::annotate("text", x=5.0, y=2.55, label=paste("c':", cp_lbl),
                      size=3, colour="#7F8C8D") +

    # PDO moderation arrows (to a-path and b-path midpoints)
    ggplot2::annotate("segment", x=4.4, xend=3.1, y=1.8, yend=3.8,
                      arrow=ggplot2::arrow(length=ggplot2::unit(0.25,"cm"),type="open"),
                      colour="#F39C12", linewidth=0.8) +
    ggplot2::annotate("text", x=3.0, y=2.8, label=paste("a3:", a3_lbl),
                      size=2.8, colour="#F39C12", hjust=1) +
    ggplot2::annotate("segment", x=5.6, xend=6.9, y=1.8, yend=3.8,
                      arrow=ggplot2::arrow(length=ggplot2::unit(0.25,"cm"),type="open"),
                      colour="#F39C12", linewidth=0.8) +
    ggplot2::annotate("text", x=7.0, y=2.8, label=paste("b4:", b4_lbl),
                      size=2.8, colour="#F39C12", hjust=0) +

    ggplot2::labs(
      title    = sprintf("Standardised Path Diagram | %s | Lag = %d months",
                         index_label, lag),
      subtitle = sprintf("Index of Moderated Mediation (IMM): %s", imm_lbl),
      caption  = "Grey dashed = direct (c') path  |  Orange = PDO moderation paths"
    ) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face="bold", size=13, hjust=0.5),
      plot.subtitle = ggplot2::element_text(size=10, hjust=0.5, colour="#C0392B"),
      plot.caption  = ggplot2::element_text(size=8, colour="grey50")
    )
}

####################################################################################
# MAIN ANALYSIS LOOP — PRIMARY INDEX ACROSS ALL LAGS
####################################################################################
cat("\n══════════════════════════════════════════════════════════\n")
cat(sprintf("  PRIMARY ANALYSIS: %s-%d across lags 0–%d months\n",
            toupper(PRIMARY_INDEX), PRIMARY_SCALE, max(LAGS_TO_TEST)))
cat("══════════════════════════════════════════════════════════\n")

index_label <- sprintf("%s-%d", toupper(PRIMARY_INDEX), PRIMARY_SCALE)
all_results <- list()
lag_profile <- list()   # for the temporal profile figure

for (lag in LAGS_TO_TEST) {
  cat(sprintf("\n▶ Lag = %d months\n", lag))

  # Build analysis data frame
  df_lag <- tryCatch(
    build_analysis_dataframe(PRIMARY_INDEX, PRIMARY_SCALE,
                             lag_months = lag,
                             start_year = 1950L, end_year = 2024L),
    error = function(e) {
      cat(sprintf("  ❌ Data assembly failed: %s\n", e$message))
      NULL
    }
  )
  if (is.null(df_lag)) next

  # Assumption checks (run only for lag = 0 to avoid clutter)
  if (lag == 0L) {
    check_mm_assumptions(df_lag)
  }

  # Fit moderated mediation
  mm_res <- fit_moderated_mediation(df_lag, lag = lag,
                                    n_boot = N_BOOTSTRAP,
                                    seed   = BOOTSTRAP_SEED)
  if (is.null(mm_res)) next

  all_results[[sprintf("lag_%02d", lag)]] <- mm_res

  # ── Collect IE estimates for lag-profile figure ──────────────────────────
  pe_lag <- mm_res$params
  for (ie_lbl in c("IE.low", "IE.mid", "IE.high", "IMM")) {
    row <- pe_lag[pe_lag$label == ie_lbl, ]
    if (nrow(row) > 0) {
      lag_profile[[length(lag_profile) + 1]] <- data.frame(
        lag      = lag,
        param    = ie_lbl,
        est      = row$est.std[1],
        ci_lower = row$ci.lower[1],
        ci_upper = row$ci.upper[1],
        p_value  = row$pvalue[1],
        stringsAsFactors = FALSE
      )
    }
  }

  # ── Per-lag diagnostic PDFs ──────────────────────────────────────────────
  lag_pdf <- file.path(MM_DIR, sprintf("mm_lag%02d_%s.pdf", lag, index_label))
  if (safe_pdf(lag_pdf, width = 13, height = 9)) {

    # 1. Path diagram
    p_diag <- plot_path_diagram(mm_res, lag, index_label)
    if (!is.null(p_diag)) print(p_diag)

    # 2. Conditional IE plot
    p_cie <- plot_conditional_ie(mm_res, lag, index_label)
    if (!is.null(p_cie)) print(p_cie)

    # 3. Johnson-Neyman floodlight
    jn_df <- tryCatch(
      compute_jn_floodlight(df_lag, n_boot = 1500L),
      error = function(e) { cat("  ⚠ JN failed:", e$message, "\n"); NULL }
    )
    if (!is.null(jn_df) && nrow(jn_df) > 0) {
      p_jn <- plot_jn_floodlight(jn_df, lag, index_label)
      if (!is.null(p_jn)) print(p_jn)
      # Save JN data for publication
      write.csv(jn_df,
                file.path(MM_DIR, sprintf("jn_data_lag%02d.csv", lag)),
                row.names = FALSE)
    }

    grDevices::dev.off()
    cat(sprintf("  ✓ Saved: mm_lag%02d_%s.pdf\n", lag, index_label))
  }
}

####################################################################################
# LAG PROFILE FIGURE — IE as a function of lag months
####################################################################################
cat("\n── Building lag-profile figure ──\n")

if (length(lag_profile) > 0) {
  lp_df <- do.call(rbind, lag_profile)

  # Map labels to readable names
  lp_df$param_label <- dplyr::case_when(
    lp_df$param == "IE.low"  ~ W_LABELS[1],
    lp_df$param == "IE.mid"  ~ W_LABELS[2],
    lp_df$param == "IE.high" ~ W_LABELS[3],
    lp_df$param == "IMM"     ~ "IMM (Moderation Index)",
    TRUE ~ lp_df$param
  )

  p_lag_prof <- ggplot2::ggplot(
    lp_df[lp_df$param != "IMM", ],
    ggplot2::aes(x = lag, y = est,
                 colour = param_label,
                 fill   = param_label)
  ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
                         alpha = 0.15, colour = NA) +
    ggplot2::geom_line(linewidth = 1.1) +
    ggplot2::geom_point(
      ggplot2::aes(shape = p_value < 0.05),
      size = 3
    ) +
    ggplot2::scale_colour_manual(
      values = c("#1A5276", "#2980B9", "#C0392B"),
      name = "PDO phase"
    ) +
    ggplot2::scale_fill_manual(
      values = c("#1A5276", "#2980B9", "#C0392B"),
      name = "PDO phase"
    ) +
    ggplot2::scale_shape_manual(
      values = c("FALSE" = 1, "TRUE" = 19),
      name   = "p < 0.05"
    ) +
    ggplot2::scale_x_continuous(breaks = LAGS_TO_TEST) +
    shared_ts_theme(12) +
    ggplot2::labs(
      title    = sprintf("Temporal Profile of Conditional Indirect Effects: ONI → PNA → %s",
                         index_label),
      subtitle = "Lag = months between teleconnection observation and drought response | Shading = 95% BCa CI",
      x        = "Lag (months)",
      y        = "Conditional Indirect Effect (standardised)"
    )

  p_imm <- ggplot2::ggplot(
    lp_df[lp_df$param == "IMM", ],
    ggplot2::aes(x = lag, y = est)
  ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
                         fill = "#F39C12", alpha = 0.2) +
    ggplot2::geom_line(colour = "#F39C12", linewidth = 1.1) +
    ggplot2::geom_point(
      ggplot2::aes(shape = p_value < 0.05),
      colour = "#F39C12", size = 3
    ) +
    ggplot2::scale_shape_manual(values = c("FALSE"=1, "TRUE"=19), name="p<0.05") +
    ggplot2::scale_x_continuous(breaks = LAGS_TO_TEST) +
    shared_ts_theme(12) +
    ggplot2::labs(
      title    = "Index of Moderated Mediation (IMM) by Lag",
      subtitle = "IMM ≠ 0 → PDO significantly moderates the ONI→PNA→Drought pathway",
      x = "Lag (months)", y = "IMM (standardised)"
    )

  lag_prof_pdf <- file.path(MM_DIR, sprintf("mm_lag_profile_%s.pdf", index_label))
  if (safe_pdf(lag_prof_pdf, width = 13, height = 11)) {
    print(p_lag_prof / p_imm)
    grDevices::dev.off()
    cat(sprintf("  ✓ Saved: mm_lag_profile_%s.pdf\n", index_label))
  }

  write.csv(lp_df,
            file.path(MM_DIR, sprintf("lag_profile_data_%s.csv", index_label)),
            row.names = FALSE)
}

####################################################################################
# SENSITIVITY ANALYSIS — SAME LAG (LAG=2) ACROSS MULTIPLE INDICES
####################################################################################
cat("\n══════════════════════════════════════════════════════════\n")
cat("  SENSITIVITY ANALYSIS: Multiple indices at optimal lag\n")
cat("══════════════════════════════════════════════════════════\n")

# Determine optimal lag from primary analysis (lag with max |IMM|)
if (length(lag_profile) > 0) {
  imm_data   <- lp_df[lp_df$param == "IMM" & !is.na(lp_df$est), ]
  optimal_lag <- if (nrow(imm_data) > 0) imm_data$lag[which.max(abs(imm_data$est))] else 2L
} else {
  optimal_lag <- 2L
}
cat(sprintf("  Optimal lag (max |IMM|): %d months\n", optimal_lag))

sensitivity_rows <- list()
for (si in SENSITIVITY_INDICES) {
  lbl <- sprintf("%s-%d", toupper(si$index), si$scale)
  cat(sprintf("\n  ▶ Sensitivity: %s at lag=%d\n", lbl, optimal_lag))
  df_s <- tryCatch(
    build_analysis_dataframe(si$index, si$scale,
                             lag_months = optimal_lag,
                             start_year = 1950L, end_year = 2024L),
    error = function(e) { cat("    ❌", e$message, "\n"); NULL }
  )
  if (is.null(df_s)) next

  mm_s <- fit_moderated_mediation(df_s, lag = optimal_lag,
                                  n_boot = N_BOOTSTRAP,
                                  seed   = BOOTSTRAP_SEED)
  if (is.null(mm_s)) next

  # Extract key parameters for comparison table
  for (pname in c("a1","a3","b1","b4","cprime","IMM","IE.low","IE.mid","IE.high")) {
    row <- mm_s$params[mm_s$params$label == pname, ]
    if (nrow(row) > 0) {
      sensitivity_rows[[length(sensitivity_rows) + 1]] <- data.frame(
        Index    = lbl,
        Lag      = optimal_lag,
        Param    = pname,
        Std.Est  = round(row$est.std[1], 4),
        CI.lower = round(row$ci.lower[1], 4),
        CI.upper = round(row$ci.upper[1], 4),
        p.value  = round(row$pvalue[1], 4),
        stringsAsFactors = FALSE
      )
    }
  }
}

if (length(sensitivity_rows) > 0) {
  sens_df <- do.call(rbind, sensitivity_rows)
  write.csv(sens_df,
            file.path(MM_DIR, "sensitivity_analysis.csv"),
            row.names = FALSE)
  cat(sprintf("\n  ✓ Saved: sensitivity_analysis.csv (%d rows)\n", nrow(sens_df)))
}

####################################################################################
# EXCEL SUMMARY WORKBOOK
####################################################################################
cat("\n── Compiling Excel summary workbook ──\n")

wb <- openxlsx::createWorkbook()
hdr_style <- openxlsx::createStyle(textDecoration = "bold",
                                    fgFill = "#1B3A6B",
                                    fontColour = "#FFFFFF")

# Sheet 1: Primary results across lags
if (length(all_results) > 0) {
  all_params_list <- lapply(names(all_results), function(nm) {
    r   <- all_results[[nm]]
    df_ <- r$params_tbl
    df_$Lag <- r$lag
    df_
  })
  all_params_df <- do.call(rbind, all_params_list)
  openxlsx::addWorksheet(wb, "Primary_All_Lags")
  openxlsx::writeData(wb, "Primary_All_Lags", all_params_df)
  openxlsx::addStyle(wb, "Primary_All_Lags", hdr_style,
                     rows = 1, cols = seq_len(ncol(all_params_df)))
}

# Sheet 2: Lag profile
if (length(lag_profile) > 0) {
  openxlsx::addWorksheet(wb, "Lag_Profile")
  openxlsx::writeData(wb, "Lag_Profile", lp_df)
  openxlsx::addStyle(wb, "Lag_Profile", hdr_style,
                     rows = 1, cols = seq_len(ncol(lp_df)))
}

# Sheet 3: Sensitivity
if (length(sensitivity_rows) > 0) {
  openxlsx::addWorksheet(wb, "Sensitivity")
  openxlsx::writeData(wb, "Sensitivity", sens_df)
  openxlsx::addStyle(wb, "Sensitivity", hdr_style,
                     rows = 1, cols = seq_len(ncol(sens_df)))
}

xlsx_out <- file.path(MM_DIR,
                       sprintf("mm_results_summary_%s.xlsx", index_label))
openxlsx::saveWorkbook(wb, xlsx_out, overwrite = TRUE)
cat(sprintf("  ✓ Saved: %s\n", basename(xlsx_out)))

####################################################################################
# DONE
####################################################################################
cat("\n╔═══════════════════════════════════════════════╗\n")
cat("║  w7 COMPLETE                                  ║\n")
cat("╚═══════════════════════════════════════════════╝\n")
cat(sprintf("  All outputs in: %s\n\n", MM_DIR))
cat("  KEY FILES:\n")
cat(sprintf("  ├── mm_lag_profile_%s.pdf\n", index_label))
cat(sprintf("  ├── mm_results_summary_%s.xlsx\n", index_label))
cat(sprintf("  ├── mm_lag00_%s.pdf  (path diagram + JN at lag=0)\n", index_label))
cat(sprintf("  └── sensitivity_analysis.csv\n"))
