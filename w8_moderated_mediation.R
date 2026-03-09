####################################################################################
# w8_moderated_mediation_gam.R
# DUAL ANALYSIS: Moderated Mediation (lavaan) + GAM with Tensor Smooths (mgcv)
####################################################################################
#
# ANALYTICAL STRATEGY (three-method architecture):
#
#  TRACK 1 — GAM tensor smooths (mgcv)   ← run first: no causal assumptions
#    "What is the full nonlinear response surface of drought to any combination
#     of ONI, PDO, PNA?"
#    Model: SPEI-6 ~ s(ONI) + s(PDO) + s(PNA) + ti(ONI,PDO) + ti(ONI,PNA)
#                  + ti(PDO,PNA) + s(month,bs="cc")
#    → 3-panel response surface figure (publication Fig 2)
#    → Variance partition: how much does each term explain?
#    → Concurvity check (GAM analog of VIF)
#
#  TRACK 2 — Moderated Mediation (lavaan)  ← causal decomposition
#    "How much of the ONI drought signal is transmitted through PNA,
#     and does PDO change that fraction?"
#    Model: Hayes (2018) Model 7 with BCa bootstrap (5,000 reps)
#    → Conditional indirect effects IE(w) at Cool/Neutral/Warm PDO
#    → Index of Moderated Mediation (IMM) — the key inferential test
#    → Johnson-Neyman floodlight: exact PDO threshold for significance
#
#  TRACK 3 — Cross-validation of results
#    → Compare GAM ti(ONI,PDO) significance against lavaan IMM significance
#    → Compare GAM response surfaces against MM conditional IE direction/magnitude
#    → If both agree: strong convergent evidence for publication
#
# OUTPUTS (all to {WD_PATH}/moderated_mediation/):
#
#  GAM outputs:
#    gam_surfaces_lag*.pdf           — 3-panel response surface (key figure)
#    gam_partials_lag*.pdf           — partial effect plots for all smooths
#    gam_residuals_lag*.pdf          — residual diagnostics
#    gam_lag_profile.csv/.pdf        — Dev.expl, p-values, EDF across lags 0-6
#    gam_variance_partition_lag*.csv — per-term deviance contribution
#    gam_sensitivity.csv             — same GAM across SPI-6/12, SPEI-3, SWEI
#
#  MM outputs (unchanged from previous version):
#    mm_lag*_SPEI-6.pdf              — path diagram + conditional IE + JN floodlight
#    mm_lag_profile_SPEI-6.pdf       — IE(w) and IMM across lags 0-6
#    mm_results_summary_SPEI-6.xlsx  — all path coefficients + BCa CIs
#    sensitivity_analysis.csv
#
#  Joint output:
#    convergence_table_SPEI-6.csv    — side-by-side GAM vs MM results per lag
#
# Depends on: DROUGHT_ANALYSIS_utils.R, utils_teleconnection_addon.R
# Run AFTER: w6_teleconnection_prep.R
####################################################################################

rm(list=ls()); gc()
source("DROUGHT_ANALYSIS_utils.R")
source("utils_teleconnection_addon.R")

utils_load_packages(c(
  "lavaan","mgcv","ggplot2","patchwork","dplyr","openxlsx","tidyr","scales"
))

setwd(WD_PATH)
MM_DIR <- file.path(WD_PATH, "moderated_mediation")
dir.create(MM_DIR, showWarnings=FALSE, recursive=TRUE)

cat("\n╔════════════════════════════════════════════════════════╗\n")
cat("║  w7  MODERATED MEDIATION + GAM TENSOR SMOOTHS          ║\n")
cat("╚════════════════════════════════════════════════════════╝\n\n")

####################################################################################
# CONFIGURATION
####################################################################################

PRIMARY_INDEX <- "spei";  PRIMARY_SCALE  <- 6
LAGS_TO_TEST  <- 0:6
N_BOOTSTRAP   <- 5000L;   BOOTSTRAP_SEED <- 42L
W_VALUES      <- c(-1, 0, 1)
W_LABELS      <- c("Cool PDO (−1 SD)","Neutral PDO (Mean)","Warm PDO (+1 SD)")

SENSITIVITY_INDICES <- list(
  list(index="spi",  scale=6),
  list(index="spi",  scale=12),
  list(index="spei", scale=3),
  list(index="swei", scale=SWEI_SCALE)
)

# GAM configuration
GAM_K_MAIN   <- 8L    # basis dimension for univariate smooths
GAM_K_TENSOR <- 6L    # basis dimension per dimension for ti()
GAM_RESPONSE <- "drought_value"   # raw values; use "Y_c" for standardised comparison

####################################################################################
# LAVAAN MODEL BUILDER (unchanged)
####################################################################################

build_lavaan_model <- function(w_low=-1, w_mid=0, w_high=1) {
  sprintf('
  PNA_c ~ a1*ONI_c + a2*PDO_c + a3*ONI_PDO_int
  Y_c   ~ b1*PNA_c + cprime*ONI_c + b3*PDO_c + b4*PNA_PDO_int
  IE.low  := (a1 + a3*%s) * (b1 + b4*%s)
  IE.mid  := (a1 + a3*%s) * (b1 + b4*%s)
  IE.high := (a1 + a3*%s) * (b1 + b4*%s)
  IE.total := a1 * b1
  TE       := cprime + a1*b1
  PM       := (a1*b1) / TE
  IMM      := a3*b1 + a1*b4 + a3*b4
  ', w_low,w_low, w_mid,w_mid, w_high,w_high)
}

####################################################################################
# FIT MODERATED MEDIATION
####################################################################################

fit_moderated_mediation <- function(df, lag=0L, n_boot=5000L, seed=42L) {
  cat(sprintf("\n  [MM] Fitting lag=%d | n=%d\n", lag, nrow(df)))
  if (nrow(df) < 60) warning(sprintf("Only %d obs at lag=%d", nrow(df), lag))
  
  model_str <- build_lavaan_model(W_VALUES[1], W_VALUES[2], W_VALUES[3])
  set.seed(seed)
  fit <- tryCatch(
    lavaan::sem(model=model_str, data=df, se="bootstrap", bootstrap=n_boot,
                estimator="ML", fixed.x=FALSE, missing="listwise"),
    error=function(e) { cat("  ❌ lavaan:", e$message, "\n"); NULL })
  if (is.null(fit)) return(NULL)
  if (!lavaan::lavTech(fit,"converged")) warning(sprintf("MM not converged lag=%d",lag))
  
  pe <- lavaan::parameterEstimates(fit, ci=TRUE, level=0.95,
                                   boot.ci.type="bca.simple", standardized=TRUE)
  pe_def <- pe[pe$op %in% c("~",":="),]
  pe_def$label_clean <- tryCatch(
    dplyr::case_when(
      pe_def$label=="a1"      ~ "a1: ONI->PNA (mean PDO)",
      pe_def$label=="a2"      ~ "a2: PDO->PNA direct",
      pe_def$label=="a3"      ~ "a3: ONIxPDO->PNA (mod a-path)",
      pe_def$label=="b1"      ~ "b1: PNA->Y (mean PDO)",
      pe_def$label=="cprime"  ~ "c': ONI->Y direct",
      pe_def$label=="b3"      ~ "b3: PDO->Y direct",
      pe_def$label=="b4"      ~ "b4: PNAxPDO->Y (mod b-path)",
      pe_def$label=="IE.low"  ~ sprintf("IE (PDO=%.0f SD)",W_VALUES[1]),
      pe_def$label=="IE.mid"  ~ sprintf("IE (PDO=%.0f SD)",W_VALUES[2]),
      pe_def$label=="IE.high" ~ sprintf("IE (PDO=+%.0f SD)",W_VALUES[3]),
      pe_def$label=="IE.total"~ "Total IE (mean PDO)",
      pe_def$label=="TE"      ~ "Total Effect ONI->Y",
      pe_def$label=="PM"      ~ "Proportion mediated",
      pe_def$label=="IMM"     ~ "IMM (Index of Mod. Mediation)",
      TRUE ~ pe_def$label),
    error = function(e) { cat("  warning: label_clean fallback\n"); pe_def$label })
  
  # Defensive: ensure all needed columns exist before subsetting
  if (!"label_clean" %in% names(pe_def)) pe_def$label_clean <- pe_def$label
  for (nc in c("ci.lower","ci.upper","pvalue")) {
    if (!nc %in% names(pe_def)) pe_def[[nc]] <- NA_real_
  }
  est_col <- if ("est.std" %in% names(pe_def)) "est.std" else "est"
  pe_print <- pe_def[, c("label_clean", est_col, "ci.lower","ci.upper","pvalue")]
  colnames(pe_print) <- c("Parameter","Std.Est","CI.lower","CI.upper","p-value")
  pe_print[,2:4] <- round(pe_print[,2:4],4); pe_print[,5] <- round(pe_print[,5],4)
  print(pe_print, row.names=FALSE)
  
  r2 <- lavaan::lavInspect(fit,"r2")
  cat(sprintf("  R²(PNA)=%.3f | R²(Y)=%.3f\n", r2["PNA_c"], r2["Y_c"]))
  
  imm_row <- pe_def[pe_def$label=="IMM",]
  if (nrow(imm_row))
    cat(sprintf("  IMM=%.4f [%.4f,%.4f] %s\n",
                .get_est(imm_row), imm_row$ci.lower, imm_row$ci.upper,
                if (imm_row$pvalue<0.05) "✓ PDO moderates indirect path"
                else "✗ PDO moderation not significant"))
  
  list(fit=fit, params=pe_def, params_tbl=pe_print,
       r2=r2, n=nrow(df), lag=lag)
}


# Helper: safely extract est.std (or fall back to est) from a params row
.get_est <- function(row) {
  if ("est.std" %in% names(row) && length(row$est.std) > 0 && !is.na(row$est.std[1]))
    row$est.std[1]
  else if ("est" %in% names(row) && length(row$est) > 0)
    row$est[1]
  else NA_real_
}

####################################################################################
# JOHNSON-NEYMAN FLOODLIGHT
####################################################################################

compute_jn_floodlight <- function(df, n_boot=1500L,
                                  w_seq=seq(-2.5,2.5,by=0.1)) {
  cat("  [JN] Floodlight analysis...\n")
  # Build list first, filter NULLs, then bind — avoids sapply-on-dataframe bug
  res_list <- lapply(w_seq, function(w) {
    tryCatch({
      la <- lm(PNA_c ~ ONI_c + PDO_c + ONI_PDO_int, data=df)
      lb <- lm(Y_c   ~ PNA_c + ONI_c + PDO_c + PNA_PDO_int, data=df)
      ca <- coef(la); cb <- coef(lb)
      ie_w <- (ca["ONI_c"]+ca["ONI_PDO_int"]*w)*(cb["PNA_c"]+cb["PNA_PDO_int"]*w)
      set.seed(BOOTSTRAP_SEED)
      boot_ie <- replicate(n_boot, {
        idx <- sample(nrow(df),replace=TRUE)
        la2 <- tryCatch(lm(PNA_c~ONI_c+PDO_c+ONI_PDO_int,data=df[idx,]),error=function(e)NULL)
        lb2 <- tryCatch(lm(Y_c~PNA_c+ONI_c+PDO_c+PNA_PDO_int,data=df[idx,]),error=function(e)NULL)
        if (is.null(la2)||is.null(lb2)) return(NA_real_)
        ca2<-coef(la2);cb2<-coef(lb2)
        (ca2["ONI_c"]+ca2["ONI_PDO_int"]*w)*(cb2["PNA_c"]+cb2["PNA_PDO_int"]*w)
      })
      boot_ie <- boot_ie[is.finite(boot_ie)]
      # Use unname() to strip quantile name attributes so data.frame() cols stay clean
      ci_lo <- if (length(boot_ie)<100) NA_real_ else unname(quantile(boot_ie,0.025))
      ci_hi <- if (length(boot_ie)<100) NA_real_ else unname(quantile(boot_ie,0.975))
      data.frame(w=w, IE=unname(ie_w), ci_lower=ci_lo, ci_upper=ci_hi,
                 significant=!is.na(ci_lo)&&(ci_lo>0||ci_hi<0),
                 stringsAsFactors=FALSE)
    }, error=function(e) NULL)
  })
  res_list <- res_list[!vapply(res_list, is.null, logical(1))]
  if (length(res_list) == 0) return(NULL)
  results <- do.call(rbind, res_list)
  results
}

####################################################################################
# MM VISUALISATION
####################################################################################

plot_conditional_ie <- function(mm_result, lag, index_label) {
  pe <- mm_result$params
  ie_rows <- pe[pe$label %in% c("IE.low","IE.mid","IE.high"),]
  if (!nrow(ie_rows)) return(NULL)
  ie_df <- data.frame(PDO_level=factor(W_LABELS,levels=W_LABELS),
                      IE=sapply(seq_len(nrow(ie_rows)), function(i) .get_est(ie_rows[i,])),
                      ci_lo=ie_rows$ci.lower, ci_hi=ie_rows$ci.upper,
                      sig=ie_rows$pvalue<0.05)
  ggplot2::ggplot(ie_df, ggplot2::aes(x=PDO_level,y=IE,colour=sig,shape=sig)) +
    ggplot2::geom_hline(yintercept=0, linetype="dashed", colour="grey50") +
    ggplot2::geom_errorbar(ggplot2::aes(ymin=ci_lo,ymax=ci_hi),
                           width=0.15, linewidth=1.1) +
    ggplot2::geom_point(size=4) +
    ggplot2::scale_colour_manual(values=c("FALSE"="grey60","TRUE"="#C0392B"),
                                 name="p<0.05") +
    ggplot2::scale_shape_manual(values=c("FALSE"=1,"TRUE"=19), name="p<0.05") +
    shared_ts_theme(12) +
    ggplot2::labs(title=sprintf("Conditional Indirect Effect: ONI→PNA→%s",index_label),
                  subtitle=sprintf("Lag=%d months | 95%% BCa bootstrap CI | n=%d",
                                   lag, mm_result$n),
                  x="PDO Phase", y="IE (standardised)")
}

plot_jn_floodlight <- function(jn_df, lag, index_label) {
  if (is.null(jn_df)||!nrow(jn_df)) return(NULL)
  ggplot2::ggplot(jn_df, ggplot2::aes(x=w)) +
    ggplot2::geom_ribbon(data=jn_df[jn_df$significant,],
                         ggplot2::aes(ymin=ci_lower,ymax=ci_upper),
                         fill="#FADBD8", alpha=0.6) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin=ci_lower,ymax=ci_upper),
                         fill="grey85", alpha=0.4) +
    ggplot2::geom_line(ggplot2::aes(y=IE), colour="#C0392B", linewidth=1.1) +
    ggplot2::geom_hline(yintercept=0, linetype="dashed") +
    ggplot2::geom_vline(xintercept=W_VALUES, linetype="dotted",
                        colour="#2980B9", linewidth=0.6) +
    shared_ts_theme(12) +
    ggplot2::labs(title="Johnson-Neyman: Region of significance",
                  subtitle=sprintf("IE(ONI→PNA→%s) as f(PDO) | lag=%d | Pink=significant",
                                   index_label, lag),
                  x="PDO (standardised)", y="IE (standardised)")
}

plot_path_diagram <- function(mm_result, lag, index_label) {
  pe <- mm_result$params
  ge <- function(lbl) {
    row <- pe[pe$label==lbl,]
    if (!nrow(row)) return("n/a")
    sprintf("b=%.3f\np=%.3f", .get_est(row), row$pvalue[1])
  }
  ggplot2::ggplot() + ggplot2::xlim(0,10) + ggplot2::ylim(0,6) +
    ggplot2::theme_void() +
    ggplot2::annotate("rect",xmin=0.3,xmax=2.3,ymin=2.5,ymax=3.5,
                      fill="#D6EAF8",colour="#2980B9",linewidth=1.2) +
    ggplot2::annotate("text",x=1.3,y=3.0,label="ONI\n(X)",size=3.5,fontface="bold") +
    ggplot2::annotate("rect",xmin=3.8,xmax=6.2,ymin=4.2,ymax=5.2,
                      fill="#D5F5E3",colour="#27AE60",linewidth=1.2) +
    ggplot2::annotate("text",x=5.0,y=4.7,label="PNA\n(Mediator M)",size=3.5,fontface="bold") +
    ggplot2::annotate("rect",xmin=7.7,xmax=9.7,ymin=2.5,ymax=3.5,
                      fill="#FADBD8",colour="#C0392B",linewidth=1.2) +
    ggplot2::annotate("text",x=8.7,y=3.0,label=sprintf("%s\n(Y)",index_label),
                      size=3.5,fontface="bold") +
    ggplot2::annotate("rect",xmin=3.8,xmax=6.2,ymin=0.8,ymax=1.8,
                      fill="#FCF3CF",colour="#F39C12",linewidth=1.2) +
    ggplot2::annotate("text",x=5.0,y=1.3,label="PDO\n(Moderator W)",size=3.5,fontface="bold") +
    ggplot2::annotate("segment",x=2.3,xend=3.8,y=3.3,yend=4.5,
                      arrow=ggplot2::arrow(length=ggplot2::unit(0.3,"cm"),type="closed"),
                      colour="#27AE60",linewidth=1) +
    ggplot2::annotate("text",x=2.9,y=4.1,label=paste("a1:",ge("a1")),
                      size=3,colour="#27AE60",hjust=0) +
    ggplot2::annotate("segment",x=6.2,xend=7.7,y=4.5,yend=3.3,
                      arrow=ggplot2::arrow(length=ggplot2::unit(0.3,"cm"),type="closed"),
                      colour="#C0392B",linewidth=1) +
    ggplot2::annotate("text",x=6.7,y=4.1,label=paste("b1:",ge("b1")),
                      size=3,colour="#C0392B",hjust=0) +
    ggplot2::annotate("segment",x=2.3,xend=7.7,y=2.8,yend=2.8,
                      arrow=ggplot2::arrow(length=ggplot2::unit(0.3,"cm"),type="closed"),
                      colour="#7F8C8D",linewidth=0.8,linetype="dashed") +
    ggplot2::annotate("text",x=5.0,y=2.55,label=paste("c':",ge("cprime")),
                      size=3,colour="#7F8C8D") +
    ggplot2::annotate("segment",x=4.4,xend=3.1,y=1.8,yend=3.8,
                      arrow=ggplot2::arrow(length=ggplot2::unit(0.25,"cm"),type="open"),
                      colour="#F39C12",linewidth=0.8) +
    ggplot2::annotate("text",x=3.0,y=2.8,label=paste("a3:",ge("a3")),
                      size=2.8,colour="#F39C12",hjust=1) +
    ggplot2::annotate("segment",x=5.6,xend=6.9,y=1.8,yend=3.8,
                      arrow=ggplot2::arrow(length=ggplot2::unit(0.25,"cm"),type="open"),
                      colour="#F39C12",linewidth=0.8) +
    ggplot2::annotate("text",x=7.0,y=2.8,label=paste("b4:",ge("b4")),
                      size=2.8,colour="#F39C12",hjust=0) +
    ggplot2::labs(
      title=sprintf("Path Diagram: %s | Lag=%d months", index_label, lag),
      subtitle=sprintf("IMM: %s", ge("IMM")),
      caption="Dashed=direct c' path | Orange=PDO moderation") +
    ggplot2::theme(plot.title=ggplot2::element_text(face="bold",size=13,hjust=0.5),
                   plot.subtitle=ggplot2::element_text(size=10,hjust=0.5,colour="#C0392B"))
}

####################################################################################
# MAIN LOOP — BOTH TRACKS ACROSS ALL LAGS
####################################################################################

index_label    <- sprintf("%s-%d", toupper(PRIMARY_INDEX), PRIMARY_SCALE)
all_mm_results <- list()
all_gam_fits   <- list()
mm_lag_profile <- list()
gam_lag_prof   <- list()
convergence_rows <- list()

cat(sprintf("\n══════════════════════════════════════════════════════════\n"))
cat(sprintf("  PRIMARY: %s | Lags 0–%d\n", index_label, max(LAGS_TO_TEST)))
cat(sprintf("══════════════════════════════════════════════════════════\n"))

for (lag in LAGS_TO_TEST) {
  cat(sprintf("\n▶▶ LAG = %d months ─────────────────────────────────────\n", lag))
  
  df_lag <- tryCatch(
    build_analysis_dataframe(PRIMARY_INDEX, PRIMARY_SCALE,
                             lag_months=lag, start_year=1950L, end_year=2024L),
    error=function(e) { cat("  ❌ Data:", e$message, "\n"); NULL })
  if (is.null(df_lag)) next
  
  if (lag == 0L) check_mm_assumptions(df_lag)
  
  # ══════════════════════════════════════════════════════════
  # TRACK 1: GAM
  # ══════════════════════════════════════════════════════════
  cat(sprintf("\n  ─── TRACK 1: GAM (lag=%d) ───\n", lag))
  
  gam_fit <- fit_gam_full(df_lag,
                          response_var        = GAM_RESPONSE,
                          k_main              = GAM_K_MAIN,
                          k_tensor            = GAM_K_TENSOR,
                          include_seasonality = TRUE,
                          select              = TRUE)
  
  all_gam_fits[[sprintf("lag_%02d", lag)]] <- gam_fit
  
  if (!is.null(gam_fit)) {
    # Variance partition
    vp <- gam_variance_partition(gam_fit)
    if (!is.null(vp)) {
      vp$lag <- lag; vp$index <- index_label
      write.csv(vp, file.path(MM_DIR, sprintf("gam_variance_partition_lag%02d.csv",lag)),
                row.names=FALSE)
    }
    
    # Concurvity
    gam_concurvity_check(gam_fit)
    
    # Lag summary for profile
    gam_lag_prof[[length(gam_lag_prof)+1]] <- gam_lag_summary(gam_fit, lag, index_label)
    
    # Per-lag PDF
    gam_pdf <- file.path(MM_DIR, sprintf("gam_lag%02d_%s.pdf", lag, index_label))
    if (safe_pdf(gam_pdf, width=14, height=10)) {
      
      # Page 1: 3-panel response surface (THE key publication figure)
      p_surf <- plot_gam_surfaces(gam_fit, index_label=index_label, lag=lag,
                                  drought_thr=DROUGHT_THRESHOLD)
      if (!is.null(p_surf)) print(p_surf)
      
      # Page 2: Partial effects
      p_part <- plot_gam_partials(gam_fit, index_label=index_label, lag=lag)
      if (!is.null(p_part)) print(p_part)
      
      # Page 3: Residual diagnostics
      p_resid <- plot_gam_residuals(gam_fit, index_label=index_label)
      if (!is.null(p_resid)) print(p_resid)
      
      grDevices::dev.off()
      cat(sprintf("  ✓ Saved: gam_lag%02d_%s.pdf (3 pages)\n", lag, index_label))
    }
  }
  
  # ══════════════════════════════════════════════════════════
  # TRACK 2: MODERATED MEDIATION
  # ══════════════════════════════════════════════════════════
  cat(sprintf("\n  ─── TRACK 2: Moderated Mediation (lag=%d) ───\n", lag))
  
  mm_res <- fit_moderated_mediation(df_lag, lag=lag,
                                    n_boot=N_BOOTSTRAP, seed=BOOTSTRAP_SEED)
  all_mm_results[[sprintf("lag_%02d",lag)]] <- mm_res
  
  if (!is.null(mm_res)) {
    # Collect IE profile
    pe_lag <- mm_res$params
    for (ie_lbl in c("IE.low","IE.mid","IE.high","IMM")) {
      row <- pe_lag[pe_lag$label==ie_lbl,]
      if (nrow(row)) {
        # Use est.std if present (standardised), fall back to est (unstandardised)
        est_val <- if ("est.std" %in% names(row) && !is.na(row$est.std[1]))
          row$est.std[1] else if ("est" %in% names(row)) row$est[1] else NA_real_
        ci_lo_val <- if ("ci.lower" %in% names(row)) row$ci.lower[1] else NA_real_
        ci_hi_val <- if ("ci.upper" %in% names(row)) row$ci.upper[1] else NA_real_
        pv_val    <- if ("pvalue"   %in% names(row)) row$pvalue[1]   else NA_real_
        mm_lag_profile[[length(mm_lag_profile)+1]] <- data.frame(
          lag=lag, param=ie_lbl,
          est=est_val, ci_lower=ci_lo_val,
          ci_upper=ci_hi_val, p_value=pv_val,
          stringsAsFactors=FALSE)
      }
    }
    
    # Per-lag MM PDF
    mm_pdf <- file.path(MM_DIR, sprintf("mm_lag%02d_%s.pdf", lag, index_label))
    if (safe_pdf(mm_pdf, width=13, height=9)) {
      print(plot_path_diagram(mm_res, lag, index_label))
      print(plot_conditional_ie(mm_res, lag, index_label))
      jn_df <- tryCatch(compute_jn_floodlight(df_lag, n_boot=1200L),
                        error=function(e) NULL)
      if (!is.null(jn_df) && nrow(jn_df)) {
        print(plot_jn_floodlight(jn_df, lag, index_label))
        write.csv(jn_df, file.path(MM_DIR,sprintf("jn_data_lag%02d.csv",lag)),
                  row.names=FALSE)
      }
      grDevices::dev.off()
      cat(sprintf("  ✓ Saved: mm_lag%02d_%s.pdf\n", lag, index_label))
    }
  }
  
  # ══════════════════════════════════════════════════════════
  # TRACK 3: CONVERGENCE CHECK
  # ══════════════════════════════════════════════════════════
  gam_ti_p <- NA_real_; mm_imm_p <- NA_real_
  gam_dev  <- NA_real_; mm_imm_e <- NA_real_
  
  if (!is.null(gam_fit)) {
    s <- summary(gam_fit); smt <- s$s.table
    ti_idx <- grep("ti\\(ONI,PDO", rownames(smt))
    if (length(ti_idx)) gam_ti_p <- round(smt[ti_idx[1],"p-value"],4)
    gam_dev <- round(s$dev.expl*100, 2)
  }
  if (!is.null(mm_res)) {
    imm_row <- mm_res$params[mm_res$params$label=="IMM",]
    if (nrow(imm_row)) {
      mm_imm_p <- round(imm_row$pvalue[1],4)
      mm_imm_e <- round(.get_est(imm_row),4)
    }
  }
  
  convergence_rows[[length(convergence_rows)+1]] <- data.frame(
    Lag              = lag,
    GAM_Dev_Expl_pct = gam_dev,
    GAM_ti_ONI_PDO_p = gam_ti_p,
    GAM_ti_sig       = !is.na(gam_ti_p) && gam_ti_p < 0.05,
    MM_IMM_est       = mm_imm_e,
    MM_IMM_p         = mm_imm_p,
    MM_IMM_sig       = !is.na(mm_imm_p) && mm_imm_p < 0.05,
    Both_agree       = (!is.na(gam_ti_p) && gam_ti_p<0.05) &&
      (!is.na(mm_imm_p) && mm_imm_p<0.05),
    stringsAsFactors=FALSE)
}

####################################################################################
# GAM LAG PROFILE FIGURE + TABLE
####################################################################################
cat("\n── GAM lag profile ──\n")

if (length(gam_lag_prof) > 0) {
  gam_lp_df <- do.call(rbind, gam_lag_prof)
  write.csv(gam_lp_df, file.path(MM_DIR,"gam_lag_profile.csv"), row.names=FALSE)
  cat(sprintf("  Saved: gam_lag_profile.csv\n"))
  
  # Build a multi-panel profile plot
  # Panel A: Deviance explained across lags
  # Panel B: p-values for the three tensor interaction terms
  lp_long <- tidyr::pivot_longer(
    gam_lp_df,
    cols=c("p_ti_ONI_PDO","p_ti_ONI_PNA","p_ti_PDO_PNA"),
    names_to="Term", values_to="p_value")
  lp_long$Term <- dplyr::recode(lp_long$Term,
                                p_ti_ONI_PDO="ti(ONI,PDO)  ← key",
                                p_ti_ONI_PNA="ti(ONI,PNA)",
                                p_ti_PDO_PNA="ti(PDO,PNA)")
  ti_clrs <- c("ti(ONI,PDO)  ← key"="#C0392B",
               "ti(ONI,PNA)"="#2980B9",
               "ti(PDO,PNA)"="#27AE60")
  
  pA <- ggplot2::ggplot(gam_lp_df, ggplot2::aes(x=Lag, y=Dev_Expl_pct)) +
    ggplot2::geom_line(colour="#1A5276", linewidth=1.2) +
    ggplot2::geom_point(size=3, colour="#1A5276") +
    ggplot2::scale_x_continuous(breaks=LAGS_TO_TEST) +
    shared_ts_theme(11) +
    ggplot2::labs(title="GAM: Deviance explained across lags",
                  subtitle=sprintf("%s | Main+tensor smooths + seasonality | select=TRUE",
                                   index_label),
                  x="Lag (months)", y="Deviance explained (%)")
  
  pB <- ggplot2::ggplot(lp_long[!is.na(lp_long$p_value),],
                        ggplot2::aes(x=Lag, y=p_value, colour=Term)) +
    ggplot2::geom_hline(yintercept=0.05, linetype="dashed", colour="grey30") +
    ggplot2::annotate("text",x=max(LAGS_TO_TEST),y=0.05,label="p=0.05",
                      hjust=1,vjust=-0.4,size=3,colour="grey30") +
    ggplot2::geom_line(linewidth=1.0) +
    ggplot2::geom_point(size=2.5) +
    ggplot2::scale_colour_manual(values=ti_clrs) +
    ggplot2::scale_x_continuous(breaks=LAGS_TO_TEST) +
    ggplot2::scale_y_log10(breaks=c(0.001,0.01,0.05,0.1,0.5,1)) +
    shared_ts_theme(11) +
    ggplot2::labs(title="GAM tensor interaction p-values across lags",
                  subtitle="Log10 scale | Points below dashed line are significant (p<0.05)",
                  x="Lag (months)", y="p-value (log10)")
  
  gam_prof_pdf <- file.path(MM_DIR, sprintf("gam_lag_profile_%s.pdf", index_label))
  if (safe_pdf(gam_prof_pdf, width=13, height=11)) {
    print(pA / pB)
    grDevices::dev.off()
    cat(sprintf("  Saved: gam_lag_profile_%s.pdf\n", index_label))
  }
}

####################################################################################
# MM LAG PROFILE FIGURE
####################################################################################
cat("\n── MM lag profile ──\n")

if (length(mm_lag_profile) > 0) {
  lp_df <- do.call(rbind, mm_lag_profile)
  lp_df$param_label <- dplyr::case_when(
    lp_df$param=="IE.low"  ~ W_LABELS[1],
    lp_df$param=="IE.mid"  ~ W_LABELS[2],
    lp_df$param=="IE.high" ~ W_LABELS[3],
    lp_df$param=="IMM"     ~ "IMM",
    TRUE ~ lp_df$param)
  
  p_ie <- ggplot2::ggplot(
    lp_df[lp_df$param!="IMM",],
    ggplot2::aes(x=lag,y=est,colour=param_label,fill=param_label)) +
    ggplot2::geom_hline(yintercept=0, linetype="dashed", colour="grey40") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin=ci_lower,ymax=ci_upper),alpha=0.15,colour=NA) +
    ggplot2::geom_line(linewidth=1.1) +
    ggplot2::geom_point(ggplot2::aes(shape=p_value<0.05),size=3) +
    ggplot2::scale_colour_manual(values=c("#1A5276","#2980B9","#C0392B"),name="PDO phase") +
    ggplot2::scale_fill_manual(values=c("#1A5276","#2980B9","#C0392B"),name="PDO phase") +
    ggplot2::scale_shape_manual(values=c("FALSE"=1,"TRUE"=19),name="p<0.05") +
    ggplot2::scale_x_continuous(breaks=LAGS_TO_TEST) +
    shared_ts_theme(12) +
    ggplot2::labs(title=sprintf("MM: Conditional Indirect Effects across lags | %s",index_label),
                  subtitle="Shading=95% BCa CI | Filled point=p<0.05",
                  x="Lag (months)",y="Conditional IE (standardised)")
  
  p_imm <- ggplot2::ggplot(lp_df[lp_df$param=="IMM",],
                           ggplot2::aes(x=lag,y=est)) +
    ggplot2::geom_hline(yintercept=0, linetype="dashed") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin=ci_lower,ymax=ci_upper),
                         fill="#F39C12", alpha=0.25) +
    ggplot2::geom_line(colour="#F39C12", linewidth=1.1) +
    ggplot2::geom_point(ggplot2::aes(shape=p_value<0.05),colour="#F39C12",size=3) +
    ggplot2::scale_shape_manual(values=c("FALSE"=1,"TRUE"=19),name="p<0.05") +
    ggplot2::scale_x_continuous(breaks=LAGS_TO_TEST) +
    shared_ts_theme(12) +
    ggplot2::labs(title="MM: Index of Moderated Mediation (IMM) across lags",
                  subtitle="IMM≠0 → PDO significantly moderates the ONI→PNA→Drought indirect pathway",
                  x="Lag (months)",y="IMM (standardised)")
  
  mm_prof_pdf <- file.path(MM_DIR, sprintf("mm_lag_profile_%s.pdf", index_label))
  if (safe_pdf(mm_prof_pdf, width=13, height=11)) {
    print(p_ie / p_imm)
    grDevices::dev.off()
    cat(sprintf("  Saved: mm_lag_profile_%s.pdf\n", index_label))
  }
  write.csv(lp_df, file.path(MM_DIR,sprintf("lag_profile_data_%s.csv",index_label)),
            row.names=FALSE)
}

####################################################################################
# CONVERGENCE TABLE — GAM vs MM side by side
####################################################################################
cat("\n── Convergence table (GAM vs MM) ──\n")

if (length(convergence_rows) > 0) {
  conv_df <- do.call(rbind, convergence_rows)
  conv_df$Interpretation <- dplyr::case_when(
    conv_df$Both_agree  ~ "✓ CONVERGENT: both methods agree",
    conv_df$GAM_ti_sig & !conv_df$MM_IMM_sig
    ~ "GAM detects interaction; MM pathway not significant",
    !conv_df$GAM_ti_sig & conv_df$MM_IMM_sig
    ~ "MM detects mediation; GAM interaction marginal",
    TRUE                ~ "Neither significant at this lag")
  cat("\n")
  print(conv_df, row.names=FALSE)
  write.csv(conv_df, file.path(MM_DIR,sprintf("convergence_table_%s.csv",index_label)),
            row.names=FALSE)
  cat(sprintf("  Saved: convergence_table_%s.csv\n", index_label))
}

####################################################################################
# SENSITIVITY ANALYSIS (optimal lag, multiple drought indices)
####################################################################################
cat("\n── Sensitivity analysis ──\n")

if (length(gam_lag_prof) > 0) {
  gam_lp_df2 <- do.call(rbind, gam_lag_prof)
  optimal_lag <- gam_lp_df2$Lag[which.max(gam_lp_df2$Dev_Expl_pct)]
} else optimal_lag <- 2L
cat(sprintf("  Optimal lag (max GAM dev.expl): %d months\n", optimal_lag))

sens_mm  <- list(); sens_gam <- list()

for (si in SENSITIVITY_INDICES) {
  lbl <- sprintf("%s-%d", toupper(si$index), si$scale)
  cat(sprintf("\n  Sensitivity: %s at lag=%d\n", lbl, optimal_lag))
  df_s <- tryCatch(
    build_analysis_dataframe(si$index, si$scale, lag_months=optimal_lag,
                             start_year=1950L, end_year=2024L),
    error=function(e) { cat("  ❌", e$message, "\n"); NULL })
  if (is.null(df_s)) next
  
  # GAM sensitivity
  gf <- fit_gam_full(df_s, response_var=GAM_RESPONSE,
                     k_main=GAM_K_MAIN, k_tensor=GAM_K_TENSOR,
                     include_seasonality=TRUE, select=TRUE)
  sens_gam[[length(sens_gam)+1]] <- gam_lag_summary(gf, optimal_lag, lbl)
  
  # MM sensitivity
  mm_s <- fit_moderated_mediation(df_s, lag=optimal_lag,
                                  n_boot=N_BOOTSTRAP, seed=BOOTSTRAP_SEED)
  if (!is.null(mm_s)) {
    for (pn in c("a1","a3","b1","b4","cprime","IMM","IE.low","IE.mid","IE.high")) {
      row <- mm_s$params[mm_s$params$label==pn,]
      if (nrow(row))
        sens_mm[[length(sens_mm)+1]] <- data.frame(
          Index=lbl, Lag=optimal_lag, Param=pn,
          Std.Est=round(.get_est(row),4),
          CI.lo=round(row$ci.lower[1],4), CI.hi=round(row$ci.upper[1],4),
          p=round(row$pvalue[1],4), stringsAsFactors=FALSE)
    }
  }
}

if (length(sens_gam)) {
  sg <- do.call(rbind, sens_gam)
  write.csv(sg, file.path(MM_DIR,"gam_sensitivity.csv"), row.names=FALSE)
  cat("\n  GAM sensitivity summary:\n"); print(sg, row.names=FALSE)
}
if (length(sens_mm)) {
  sm <- do.call(rbind, sens_mm)
  write.csv(sm, file.path(MM_DIR,"mm_sensitivity.csv"), row.names=FALSE)
}

####################################################################################
# EXCEL WORKBOOK — all results combined
####################################################################################
cat("\n── Building Excel workbook ──\n")

wb  <- openxlsx::createWorkbook()
hdr <- openxlsx::createStyle(textDecoration="bold", fgFill="#1B3A6B",
                             fontColour="#FFFFFF")
add_sheet <- function(wb, name, df) {
  openxlsx::addWorksheet(wb, name)
  openxlsx::writeData(wb, name, df)
  openxlsx::addStyle(wb, name, hdr, rows=1, cols=seq_len(ncol(df)))
}

# MM primary results
if (length(all_mm_results)) {
  mm_all_df <- do.call(rbind, lapply(names(all_mm_results), function(nm) {
    r <- all_mm_results[[nm]]; if (is.null(r)) return(NULL)
    d <- r$params_tbl; d$Lag <- r$lag; d
  }))
  if (!is.null(mm_all_df)) add_sheet(wb, "MM_All_Lags", mm_all_df)
}

# MM lag profile
if (length(mm_lag_profile)) add_sheet(wb, "MM_Lag_Profile", do.call(rbind,mm_lag_profile))

# MM sensitivity
if (length(sens_mm))  add_sheet(wb, "MM_Sensitivity",  do.call(rbind,sens_mm))

# GAM lag profile
if (length(gam_lag_prof)) add_sheet(wb, "GAM_Lag_Profile", do.call(rbind,gam_lag_prof))

# GAM sensitivity
if (length(sens_gam)) add_sheet(wb, "GAM_Sensitivity", do.call(rbind,sens_gam))

# Convergence
if (length(convergence_rows)) add_sheet(wb, "Convergence_GAM_vs_MM",
                                        do.call(rbind,convergence_rows))

xlsx_out <- file.path(MM_DIR, sprintf("results_summary_%s.xlsx", index_label))
openxlsx::saveWorkbook(wb, xlsx_out, overwrite=TRUE)
cat(sprintf("  Saved: %s\n", basename(xlsx_out)))

####################################################################################
# DONE
####################################################################################
cat("\n╔════════════════════════════════════════════════════════╗\n")
cat("║  w7 COMPLETE                                           ║\n")
cat("╚════════════════════════════════════════════════════════╝\n")
cat(sprintf("  All outputs: %s\n\n", MM_DIR))
cat("  KEY OUTPUTS:\n")
cat(sprintf("  ├── gam_lag0*_%s.pdf          (surface + partials + residuals)\n",index_label))
cat(sprintf("  ├── gam_lag_profile_%s.pdf    (dev.expl + tensor p-values)\n",index_label))
cat(sprintf("  ├── mm_lag_profile_%s.pdf     (conditional IE + IMM)\n",index_label))
cat(sprintf("  ├── convergence_table_%s.csv  (GAM vs MM agreement per lag)\n",index_label))
cat(sprintf("  └── results_summary_%s.xlsx  (all results in one workbook)\n",index_label))