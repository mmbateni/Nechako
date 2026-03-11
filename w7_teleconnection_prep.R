####################################################################################
# w7_teleconnection_prep.R  ·  TELECONNECTION DATA PREPARATION + GAM EXPLORATION
####################################################################################
# Downloads and parses ONI, PDO, PNA. Produces QC plots and descriptive stats.
# Section 6 runs predictor-level GAMs characterising pairwise nonlinear
# relationships — informs k choices and concurvity diagnostics used.
#
# OUTPUTS (all to {WD_PATH}/teleconnections/):
#   oni_monthly.csv / pdo_monthly.csv / pna_monthly.csv
#   teleconnections_merged.csv          — wide: date, ONI, PDO, PNA, pdo_phase
#   teleconnections_stats.xlsx          — descriptive + Pearson/Spearman cors
#   pdo_phase_summary.csv
#   teleconnections_QC_plots.pdf        — 5 pages: TS / scatter / density / clim / ACF
#   teleconnections_predictor_GAM.pdf   — NEW: pairwise predictor GAM smooths
#
# Run BEFORE: w8_moderated_mediation_gam.R
# Depends on: DROUGHT_ANALYSIS_utils.R, utils_teleconnection_addon.R
####################################################################################

rm(list=ls()); gc()
source("DROUGHT_ANALYSIS_utils.R")
source("utils_teleconnection_addon.R")

utils_load_packages(c(
  "ggplot2","patchwork","dplyr","lubridate",
  "openxlsx","GGally","mgcv","tidyr"
))

setwd(WD_PATH)
dir.create(TELE_DIR, showWarnings=FALSE, recursive=TRUE)

cat("\n╔═══════════════════════════════════════════════════╗\n")
cat("║  w6  TELECONNECTION PREP + PREDICTOR GAM          ║\n")
cat("╚═══════════════════════════════════════════════════╝\n\n")

START_YEAR <- 1950L
END_YEAR   <- 2024L
tele_clrs  <- c(ONI="#E74C3C", PDO="#2980B9", PNA="#27AE60")

####################################################################################
# STEP 1: DOWNLOAD / LOAD
####################################################################################
cat("── STEP 1: Load indices ──\n\n")

oni_df <- load_teleconnection("oni", start_year=START_YEAR, end_year=END_YEAR)
pdo_df <- load_teleconnection("pdo", start_year=START_YEAR, end_year=END_YEAR)
pna_df <- load_teleconnection("pna", start_year=START_YEAR, end_year=END_YEAR)
colnames(oni_df) <- c("date","ONI")
colnames(pdo_df) <- c("date","PDO")
colnames(pna_df) <- c("date","PNA")

####################################################################################
# STEP 2: MERGE
####################################################################################
cat("\n── STEP 2: Merge ──\n")

tc_wide <- Reduce(function(a,b) merge(a,b,by="date",all=FALSE),
                  list(oni_df, pdo_df, pna_df))
tc_wide <- tc_wide[order(tc_wide$date), ]
tc_wide <- tc_wide[stats::complete.cases(tc_wide), ]
tc_wide$month <- as.integer(format(tc_wide$date, "%m"))
tc_wide$year  <- as.integer(format(tc_wide$date, "%Y"))

cat(sprintf("  Merged: %d months (%s – %s)\n",
            nrow(tc_wide), min(tc_wide$date), max(tc_wide$date)))

write.csv(tc_wide, file.path(TELE_DIR,"teleconnections_merged.csv"), row.names=FALSE)
cat("  Saved: teleconnections_merged.csv\n")

####################################################################################
# STEP 3: DESCRIPTIVE STATISTICS
####################################################################################
cat("\n── STEP 3: Descriptive statistics ──\n")

desc_stats <- data.frame(
  Index    = c("ONI","PDO","PNA"),
  N        = sapply(c("ONI","PDO","PNA"), function(v) sum(!is.na(tc_wide[[v]]))),
  Mean     = round(sapply(c("ONI","PDO","PNA"), function(v) mean(tc_wide[[v]],na.rm=TRUE)),4),
  SD       = round(sapply(c("ONI","PDO","PNA"), function(v) sd(tc_wide[[v]],na.rm=TRUE)),4),
  Min      = round(sapply(c("ONI","PDO","PNA"), function(v) min(tc_wide[[v]],na.rm=TRUE)),4),
  Max      = round(sapply(c("ONI","PDO","PNA"), function(v) max(tc_wide[[v]],na.rm=TRUE)),4),
  Skewness = round(sapply(c("ONI","PDO","PNA"), function(v) {
    x <- tc_wide[[v]][!is.na(tc_wide[[v]])]
    (sum((x-mean(x))^3)/length(x)) / sd(x)^3
  }), 4),
  stringsAsFactors=FALSE
)
print(desc_stats)

cor_pear <- round(cor(tc_wide[,c("ONI","PDO","PNA")], method="pearson",  use="complete.obs"), 3)
cor_spear<- round(cor(tc_wide[,c("ONI","PDO","PNA")], method="spearman", use="complete.obs"), 3)
cat("  Pearson:\n");  print(cor_pear)
cat("  Spearman:\n"); print(cor_spear)

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb,"Descriptive");  openxlsx::writeData(wb,"Descriptive",desc_stats)
openxlsx::addWorksheet(wb,"Pearson");      openxlsx::writeData(wb,"Pearson",as.data.frame(cor_pear))
openxlsx::addWorksheet(wb,"Spearman");     openxlsx::writeData(wb,"Spearman",as.data.frame(cor_spear))
openxlsx::saveWorkbook(wb, file.path(TELE_DIR,"teleconnections_stats.xlsx"), overwrite=TRUE)
cat("  Saved: teleconnections_stats.xlsx\n")

####################################################################################
# STEP 4: QC DIAGNOSTIC PLOTS
####################################################################################
cat("\n── STEP 4: QC plots ──\n")

tc_long <- do.call(rbind, lapply(c("ONI","PDO","PNA"), function(idx)
  data.frame(date=tc_wide$date, Index=idx, Value=tc_wide[[idx]], stringsAsFactors=FALSE)))

# Fig 1: Full time series
p_ts <- ggplot2::ggplot(tc_long, ggplot2::aes(x=date, y=Value, colour=Index)) +
  ggplot2::geom_line(linewidth=0.55, alpha=0.85) +
  ggplot2::geom_hline(yintercept=0, linetype="dashed", colour="grey40") +
  ggplot2::facet_wrap(~Index, ncol=1, scales="free_y") +
  ggplot2::scale_colour_manual(values=tele_clrs) +
  ggplot2::scale_x_date(date_breaks="10 years", date_labels="%Y") +
  shared_ts_theme(11) +
  ggplot2::labs(title="Teleconnection Indices — Nechako Watershed Attribution Study",
                subtitle=sprintf("Monthly values %d–%d | NOAA CPC / PSL", START_YEAR, END_YEAR),
                x=NULL, y="Index value") +
  ggplot2::theme(legend.position="none", strip.text=ggplot2::element_text(face="bold"))

# Fig 2: Scatter matrix
p_scatter <- GGally::ggpairs(
  tc_wide[,c("ONI","PDO","PNA")],
  upper=list(continuous=GGally::wrap("cor", method="spearman", size=4.5)),
  lower=list(continuous=GGally::wrap("points", alpha=0.2, size=0.7)),
  diag =list(continuous=GGally::wrap("densityDiag", fill="steelblue", alpha=0.4)),
  title="Pairwise relationships among teleconnection predictors (Spearman ρ)")

# Fig 3: ONI × PDO joint density
p_oniXpdo <- ggplot2::ggplot(tc_wide, ggplot2::aes(x=ONI, y=PDO)) +
  ggplot2::geom_bin2d(bins=35) +
  ggplot2::scale_fill_viridis_c(option="plasma", name="Count") +
  ggplot2::geom_hline(yintercept=c(-0.5,0.5), linetype="dashed") +
  ggplot2::geom_vline(xintercept=c(-0.5,0.5), linetype="dashed") +
  ggplot2::annotate("text",x=1.8,y=1.2,
                    label="El Niño\nWarm PDO\n(drought risk)",
                    size=3.2, colour="#C0392B", fontface="bold") +
  ggplot2::annotate("text",x=-1.8,y=-1.2,
                    label="La Niña\nCool PDO\n(wet conditions)",
                    size=3.2, colour="#1A5276", fontface="bold") +
  shared_ts_theme(11) +
  ggplot2::labs(title="Joint distribution of ONI and PDO",
                subtitle="Quadrant lines at ±0.5 define PDO and ENSO phase boundaries",
                x="ONI", y="PDO")

# Fig 4: Monthly climatology
tc_clim <- tc_long %>%
  dplyr::mutate(Month=as.integer(format(date,"%m"))) %>%
  dplyr::group_by(Index, Month) %>%
  dplyr::summarise(Mean=mean(Value,na.rm=TRUE), SD=sd(Value,na.rm=TRUE), .groups="drop")

p_clim <- ggplot2::ggplot(tc_clim,
                           ggplot2::aes(x=Month, y=Mean, colour=Index, fill=Index)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin=Mean-SD, ymax=Mean+SD), alpha=0.15, colour=NA) +
  ggplot2::geom_line(linewidth=1) +
  ggplot2::geom_hline(yintercept=0, linetype="dashed", colour="grey40") +
  ggplot2::scale_colour_manual(values=tele_clrs) +
  ggplot2::scale_fill_manual(values=tele_clrs) +
  ggplot2::scale_x_continuous(breaks=1:12, labels=month.abb) +
  shared_ts_theme(11) +
  ggplot2::labs(title="Seasonal climatology of teleconnection indices",
                subtitle="Mean ± 1 SD by calendar month", x="Month", y="Mean value")

# Fig 5: ACF for each index
acf_df <- do.call(rbind, lapply(c("ONI","PDO","PNA"), function(idx) {
  ac <- acf(tc_wide[[idx]], lag.max=36, plot=FALSE)
  data.frame(Index=idx, Lag=as.integer(ac$lag[-1]),
             ACF=as.numeric(ac$acf[-1]), stringsAsFactors=FALSE)
}))
ci_acf <- qnorm(0.975)/sqrt(nrow(tc_wide))
p_acf <- ggplot2::ggplot(acf_df, ggplot2::aes(x=Lag, y=ACF, colour=Index)) +
  ggplot2::geom_hline(yintercept=0) +
  ggplot2::geom_hline(yintercept=c(-ci_acf,ci_acf), linetype="dashed",
                      colour="blue", linewidth=0.5) +
  ggplot2::geom_segment(ggplot2::aes(xend=Lag,yend=0,colour=Index)) +
  ggplot2::geom_point(size=1.2) +
  ggplot2::facet_wrap(~Index, ncol=1) +
  ggplot2::scale_colour_manual(values=tele_clrs) +
  shared_ts_theme(10) +
  ggplot2::labs(title="Autocorrelation functions of teleconnection indices",
                subtitle="Blue dashed = 95% CI under white noise | Persistence informs effective N",
                x="Lag (months)", y="ACF") +
  ggplot2::theme(legend.position="none")

pdf_qc <- file.path(TELE_DIR,"teleconnections_QC_plots.pdf")
if (safe_pdf(pdf_qc, width=13, height=9)) {
  print(p_ts); print(p_scatter); print(p_oniXpdo); print(p_clim); print(p_acf)
  grDevices::dev.off()
  cat("  Saved: teleconnections_QC_plots.pdf (5 pages)\n")
}

####################################################################################
# STEP 5: PDO PHASE CLASSIFICATION
####################################################################################
cat("\n── STEP 5: PDO phase classification ──\n")

tc_wide$pdo_phase <- cut(tc_wide$PDO, breaks=c(-Inf,-0.5,0.5,Inf),
                          labels=c("Cool","Neutral","Warm"), right=FALSE)

phase_summary <- tc_wide %>%
  dplyr::group_by(pdo_phase) %>%
  dplyr::summarise(
    N_months  = dplyr::n(),
    Pct       = round(100*dplyr::n()/nrow(tc_wide),1),
    Mean_ONI  = round(mean(ONI,na.rm=TRUE),3),
    Mean_PNA  = round(mean(PNA,na.rm=TRUE),3),
    Mean_PDO  = round(mean(PDO,na.rm=TRUE),3),
    .groups="drop")
print(phase_summary)
write.csv(phase_summary, file.path(TELE_DIR,"pdo_phase_summary.csv"), row.names=FALSE)
write.csv(tc_wide, file.path(TELE_DIR,"teleconnections_merged.csv"), row.names=FALSE)
cat("  Saved: pdo_phase_summary.csv + teleconnections_merged.csv (with pdo_phase)\n")

####################################################################################
# STEP 6 (NEW): PREDICTOR-LEVEL GAMs
# Purpose: characterise the nonlinear pairwise relationships BETWEEN teleconnection
# predictors before fitting the full drought~f(ONI,PDO,PNA) model in w7.
# This step:
#   (a) Fits three bivariate GAMs: PDO~s(ONI), PNA~s(ONI), PNA~s(PDO)
#   (b) Plots the fitted smooths with partial residuals
#   (c) Reports deviance explained and effective df — if EDF≈1 the relationship is
#       well-approximated by a linear term, confirming that the full tensor model
#       in w7 is not adding spurious nonlinearity
#   (d) Fits a trivariate GAM: PNA ~ s(ONI) + s(PDO) + ti(ONI,PDO) to check whether
#       the ONI→PNA relationship is meaningfully modulated by PDO state
#       (empirical confirmation that the moderated mediation causal graph is
#       supported by the teleconnection data alone)
####################################################################################
cat("\n── STEP 6 (NEW): Predictor-level GAMs ──\n")
cat("  Purpose: characterise nonlinear predictor relationships before w7\n\n")

# ── 6a: Bivariate GAMs ────────────────────────────────────────────────────────
bivar_models <- list(
  list(y="PDO", x="ONI", label="PDO ~ s(ONI)"),
  list(y="PNA", x="ONI", label="PNA ~ s(ONI)  [a-path check]"),
  list(y="PNA", x="PDO", label="PNA ~ s(PDO)")
)

bivar_results <- list()
bivar_plots   <- list()

for (bm in bivar_models) {
  cat(sprintf("  Fitting: %s\n", bm$label))
  frm <- as.formula(sprintf("%s ~ s(%s, k=8, bs='cr')", bm$y, bm$x))
  fit <- tryCatch(
    mgcv::gam(frm, data=tc_wide, method="REML"),
    error=function(e) { cat("  ❌", e$message, "\n"); NULL })
  if (is.null(fit)) next

  s <- summary(fit)
  edf <- round(s$s.table[1,"edf"], 2)
  pv  <- round(s$s.table[1,"p-value"], 4)
  dev <- round(s$dev.expl*100, 1)
  cat(sprintf("    EDF=%.2f | p=%.4f | Dev.expl=%.1f%%\n", edf, pv, dev))
  cat(sprintf("    Interpretation: %s\n",
              if (edf < 1.5) "near-linear relationship — modest nonlinearity"
              else if (edf < 3) "moderate nonlinearity"
              else "substantial nonlinearity — tensor interactions in w7 are justified"))

  bivar_results[[bm$label]] <- list(label=bm$label, edf=edf, p=pv, dev_expl=dev)

  # Build partial effect plot via predict grid
  x_seq <- seq(range(tc_wide[[bm$x]], na.rm=TRUE)[1],
               range(tc_wide[[bm$x]], na.rm=TRUE)[2], length.out=200)
  nd    <- data.frame(x_seq); colnames(nd) <- bm$x
  pr    <- mgcv::predict.gam(fit, newdata=nd, type="link", se.fit=TRUE)
  pred_df <- data.frame(x=x_seq, fit=pr$fit,
                         lo=pr$fit-1.96*pr$se.fit,
                         hi=pr$fit+1.96*pr$se.fit)

  # Partial residuals
  pr_resid <- data.frame(
    x    = tc_wide[[bm$x]],
    resid= residuals(fit, type="working") + fitted(fit) - mean(fitted(fit))
  )

  bivar_plots[[bm$label]] <- ggplot2::ggplot() +
    ggplot2::geom_point(data=pr_resid,
                        ggplot2::aes(x=x, y=resid), alpha=0.2, size=0.7, colour="grey50") +
    ggplot2::geom_ribbon(data=pred_df,
                         ggplot2::aes(x=x, ymin=lo, ymax=hi),
                         fill="#AED6F1", alpha=0.5) +
    ggplot2::geom_line(data=pred_df,
                       ggplot2::aes(x=x, y=fit), colour="#1A5276", linewidth=1.1) +
    ggplot2::geom_hline(yintercept=0, linetype="dashed", colour="grey40") +
    shared_ts_theme(10) +
    ggplot2::labs(
      title    = bm$label,
      subtitle = sprintf("EDF=%.2f | p=%.4f | Dev.expl=%.1f%%", edf, pv, dev),
      x        = bm$x, y=bm$y)
}

# ── 6b: Trivariate GAM — PNA ~ s(ONI) + s(PDO) + ti(ONI,PDO) ─────────────────
# This is the empirical test of the a-path from the MM model:
# Does PDO state change the ONI→PNA relationship in the data?
cat("\n  Fitting trivariate predictor GAM: PNA ~ s(ONI) + s(PDO) + ti(ONI,PDO)\n")
fit_pna_trivar <- tryCatch(
  mgcv::gam(PNA ~ s(ONI,k=8,bs="cr") + s(PDO,k=8,bs="cr") + ti(ONI,PDO,k=c(6,6)),
            data=tc_wide, method="REML", select=TRUE),
  error=function(e) { cat("  ❌", e$message, "\n"); NULL })

p_pna_surface <- NULL
if (!is.null(fit_pna_trivar)) {
  s3 <- summary(fit_pna_trivar)
  cat(sprintf("  PNA trivar GAM: Dev.expl=%.1f%% | adj-R²=%.3f\n",
              s3$dev.expl*100, s3$r.sq))
  cat("  Term p-values:\n")
  print(round(s3$s.table[,c("edf","p-value")], 4))

  # Is ti(ONI,PDO) significant? If yes, the tensor interaction in w7 is empirically warranted
  ti_row <- grep("ti\\(ONI,PDO", rownames(s3$s.table))
  if (length(ti_row)) {
    ti_p <- s3$s.table[ti_row,"p-value"]
    cat(sprintf("\n  ti(ONI,PDO) in PNA model: p=%.4f  %s\n", ti_p,
                if (ti_p < 0.05) "✓ Significant — PDO modulates ONI→PNA (supports MM a-path)"
                else "Not significant — ONI→PNA relationship is PDO-phase independent"))
  }

  # Response surface: PNA as a function of (ONI, PDO)
  grid3 <- predict_gam_surface(fit_pna_trivar, x_var="ONI", y_var="PDO",
                                cond_list=list(), n_grid=80L)
  if (!is.null(grid3)) {
    p_pna_surface <- ggplot2::ggplot(grid3, ggplot2::aes(x=ONI, y=PDO, fill=fit)) +
      ggplot2::geom_tile() +
      ggplot2::geom_contour(ggplot2::aes(z=fit), colour="white",
                            linewidth=0.5, breaks=seq(-2,2,by=0.5)) +
      ggplot2::scale_fill_gradientn(
        colours=c("#1A5276","#AED6F1","white","#FAD7A0","#C0392B"),
        name="E[PNA]") +
      ggplot2::geom_hline(yintercept=c(-0.5,0.5), linetype="dashed",
                          colour="grey80", linewidth=0.5) +
      ggplot2::geom_vline(xintercept=c(-0.5,0.5), linetype="dashed",
                          colour="grey80", linewidth=0.5) +
      shared_ts_theme(11) +
      ggplot2::labs(
        title    = "GAM surface: E[PNA | ONI, PDO]",
        subtitle = sprintf("Empirical test of MM a-path: does PDO modulate ONI→PNA? | Dev.expl=%.1f%%",
                           s3$dev.expl*100),
        x="ONI", y="PDO",
        caption="If colour gradient shifts meaningfully across PDO axis, tensor interaction is warranted in w7")
  }
}

# ── Save predictor GAM PDF ────────────────────────────────────────────────────
pdf_gam <- file.path(TELE_DIR,"teleconnections_predictor_GAM.pdf")
if (safe_pdf(pdf_gam, width=13, height=9)) {

  # Page 1: three bivariate smooth plots
  if (length(bivar_plots) >= 3) {
    print(bivar_plots[[1]] + bivar_plots[[2]] + bivar_plots[[3]] +
            patchwork::plot_layout(ncol=3) +
            patchwork::plot_annotation(
              title="Bivariate GAM smooths among teleconnection predictors",
              subtitle="Points = partial residuals | Shading = 95% CI | EDF > 1.5 justifies nonlinear terms in w7"))
  } else {
    for (pp in bivar_plots) print(pp)
  }

  # Page 2: PNA surface
  if (!is.null(p_pna_surface)) print(p_pna_surface)

  grDevices::dev.off()
  cat(sprintf("\n  Saved: teleconnections_predictor_GAM.pdf\n"))
}

# ── Summary table for bivariate GAMs ─────────────────────────────────────────
if (length(bivar_results) > 0) {
  bv_df <- do.call(rbind, lapply(bivar_results, function(r)
    data.frame(Model=r$label, EDF=r$edf, p_value=r$p, Dev_Expl_pct=r$dev_expl,
               Justifies_Nonlinear=r$edf>1.5, stringsAsFactors=FALSE)))
  cat("\n  Bivariate predictor GAM summary:\n")
  print(bv_df, row.names=FALSE)
  write.csv(bv_df, file.path(TELE_DIR,"predictor_GAM_summary.csv"), row.names=FALSE)
  cat("  Saved: predictor_GAM_summary.csv\n")
}

####################################################################################
# DONE
####################################################################################
cat("\n╔═══════════════════════════════════════════════════╗\n")
cat("║  w7 COMPLETE                                      ║\n")
cat("╚═══════════════════════════════════════════════════╝\n")
cat("  Key files written to:", TELE_DIR, "\n")
cat("  Next: run w8_moderated_mediation_gam.R\n\n")
