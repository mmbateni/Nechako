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

# ── Windows temp-directory guard ─────────────────────────────────────────────
# Windows periodically cleans the R session temp directory while R is running.
# openxlsx::saveWorkbook() creates a subdirectory under tempdir(); if the parent
# is gone, the call fails with "Failed to create temporary directory …".
# Recreating tempdir() here (and before each saveWorkbook) is sufficient.
# If recreate still fails, redirect to a writable location on the project drive.
.fix_tempdir <- function() {
  td <- tempdir()
  if (!dir.exists(td)) {
    ok <- dir.create(td, recursive = TRUE, showWarnings = FALSE)
    if (!ok) {
      # Hard fallback: use a subdirectory of the working directory
      td_alt <- file.path(getwd(), ".rtmp_wb")
      dir.create(td_alt, recursive = TRUE, showWarnings = FALSE)
      Sys.setenv(TMPDIR = td_alt, TMP = td_alt, TEMP = td_alt)
      cat(sprintf("  \u26a0 R tempdir() recreated at: %s\n", td_alt))
    } else {
      cat(sprintf("  \u26a0 R tempdir() was missing; recreated at: %s\n", td))
    }
  }
}

# Wrapper: ensures tempdir() exists before every saveWorkbook call
.safe_save_workbook <- function(wb, path, overwrite = TRUE) {
  .fix_tempdir()
  tryCatch(
    openxlsx::saveWorkbook(wb, path, overwrite = overwrite),
    error = function(e) {
      cat(sprintf("  \u26a0 saveWorkbook failed (%s): %s\n",
                  basename(path), conditionMessage(e)))
      cat("    Retrying after temp-dir reset...\n")
      td_alt <- file.path(normalizePath("."), ".rtmp_wb")
      dir.create(td_alt, recursive = TRUE, showWarnings = FALSE)
      Sys.setenv(TMPDIR = td_alt, TMP = td_alt, TEMP = td_alt)
      .fix_tempdir()
      tryCatch(
        openxlsx::saveWorkbook(wb, path, overwrite = overwrite),
        error = function(e2) cat(sprintf(
          "  \u274c saveWorkbook still failed: %s\n", conditionMessage(e2)))
      )
    }
  )
}

.fix_tempdir()   # run once at startup
# ─────────────────────────────────────────────────────────────────────────────

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

####################################################################################
# STEP 3: DESCRIPTIVE STATISTICS
####################################################################################
cat("\n── STEP 3: Descriptive statistics ──\n")

# Compute all descriptive stats in one pass per variable
desc_stats <- dplyr::bind_rows(lapply(c("ONI","PDO","PNA"), function(v) {
  x <- tc_wide[[v]][!is.na(tc_wide[[v]])]
  data.frame(
    Index    = v,
    N        = length(x),
    Mean     = round(mean(x), 4),
    SD       = round(sd(x),   4),
    Min      = round(min(x),  4),
    Max      = round(max(x),  4),
    Skewness = round((sum((x - mean(x))^3) / length(x)) / sd(x)^3, 4),
    stringsAsFactors = FALSE
  )
}))
print(desc_stats)

cor_pear <- round(cor(tc_wide[,c("ONI","PDO","PNA")], method="pearson",  use="complete.obs"), 3)
cor_spear<- round(cor(tc_wide[,c("ONI","PDO","PNA")], method="spearman", use="complete.obs"), 3)
cat("  Pearson:\n");  print(cor_pear)
cat("  Spearman:\n"); print(cor_spear)

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb,"Descriptive");  openxlsx::writeData(wb,"Descriptive",desc_stats)
openxlsx::addWorksheet(wb,"Pearson");      openxlsx::writeData(wb,"Pearson",as.data.frame(cor_pear))
openxlsx::addWorksheet(wb,"Spearman");     openxlsx::writeData(wb,"Spearman",as.data.frame(cor_spear))
.safe_save_workbook(wb, file.path(TELE_DIR,"teleconnections_stats.xlsx"))
cat("  Saved: teleconnections_stats.xlsx\n")

####################################################################################
# STEP 4: QC DIAGNOSTIC PLOTS
####################################################################################
cat("\n── STEP 4: QC plots ──\n")

tc_long <- dplyr::bind_rows(lapply(c("ONI","PDO","PNA"), function(idx)
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
# STEP 7: MANUSCRIPT FIGURE 9 — Teleconnection context for drought record
#
# Requires the basin SPEI-3 time series produced by w1/w2.  Loaded via the same
# load_basin_avg_csv() helper used throughout the pipeline.  If the file is
# absent the figure is skipped gracefully with a clear warning.
#
# Two complementary panels saved as a single figure:
#
#   Fig9a  —  Three-panel time series
#     Row 1: ONI (monthly, ±0.5 phase boundary, ENSO event annotations)
#     Row 2: PDO (monthly, ±0.5 phase boundary, warm/cool shading)
#     Row 3: Basin SPEI-3 with colour fill by drought category
#     All panels share a common x-axis 1950–2025; 2022–2025 grey band on all.
#
#   Fig9b  —  PDO-phase composite bar chart
#     Mean SPEI-3 ± bootstrapped 95 % CI stratified by PDO phase
#     (Cool / Neutral / Warm, threshold ±0.5).  A horizontal reference at 0
#     and at −0.5 (event threshold) is shown.  Each bar labelled with N.
#     Inset violin shows the full distribution within each phase.
#
# Outputs (TELE_DIR):
#   Fig9a_teleconnection_timeseries.pdf/.png    — 7.5 × 8.5 in, 300 DPI
#   Fig9b_PDO_phase_SPEI3_composite.pdf/.png    — 6.0 × 5.0 in, 300 DPI
####################################################################################
cat("\n── STEP 7: Manuscript Figure 9 (teleconnection context) ──\n")

tryCatch({
  
  utils_load_packages(c("scales"))
  
  # ── shared constants ────────────────────────────────────────────────────────
  HIGHLIGHT_START <- as.Date("2022-01-01")
  HIGHLIGHT_END   <- as.Date("2025-12-31")
  THR_EVENT       <- -0.5    # event threshold from w5
  
  COL_WARM  <- "#FADBD8"   # pale red  — warm PDO / El Niño fill
  COL_COOL  <- "#D6EAF8"   # pale blue — cool PDO / La Niña fill
  COL_ONI   <- tele_clrs["ONI"]
  COL_PDO   <- tele_clrs["PDO"]
  COL_SPEI  <- "#2c7bb6"   # blue line for SPEI-3
  
  theme_fig9 <- ggplot2::theme_classic(base_size = 9) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(size = 9,  face = "bold", hjust = 0),
      plot.subtitle = ggplot2::element_text(size = 7.5, colour = "grey40",
                                            hjust = 0,
                                            margin = ggplot2::margin(b = 3)),
      axis.text     = ggplot2::element_text(size = 7.5),
      axis.title    = ggplot2::element_text(size = 8.0),
      panel.grid.major = ggplot2::element_line(colour = "grey94", linewidth = 0.25),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin   = ggplot2::margin(3, 6, 2, 4)
    )
  
  # ── load SPEI-3 basin time series ───────────────────────────────────────────
  spei3_raw <- tryCatch(
    load_basin_avg_csv("spei", 3),
    error = function(e) {
      cat(sprintf("  \u26a0 SPEI-3 not found: %s\n  Trying load_any_ts()...\n",
                  e$message))
      NULL
    })
  
  if (is.null(spei3_raw)) {
    spei3_raw <- tryCatch(
      {
        f <- file.path(WD_PATH, "temporal_drought", "basin_averaged_timeseries",
                       "spei_03_basin_average.csv")
        if (!file.exists(f)) stop("CSV not found")
        df <- data.table::fread(f)
        df$date  <- as.Date(df$date)
        df$value <- as.numeric(df$value)
        as.data.frame(df)
      },
      error = function(e) {
        cat(sprintf("  \u26a0 SPEI-3 basin CSV also missing: %s\n", e$message))
        NULL
      })
  }
  
  if (is.null(spei3_raw)) {
    cat("  \u26a0 Cannot build Fig9: SPEI-3 data unavailable. Run w1/w2 first.\n")
    stop("SPEI-3 missing")
  }
  
  spei3 <- as.data.frame(spei3_raw)
  names(spei3)[which(names(spei3) %in% c("date","Date"))[1]] <- "date"
  names(spei3)[which(names(spei3) %in% c("value","Value","Mean_Value"))[1]] <- "value"
  spei3$date  <- as.Date(spei3$date)
  spei3$value <- as.numeric(spei3$value)
  
  # align to teleconnection date range
  spei3 <- spei3[spei3$date >= min(tc_wide$date) &
                   spei3$date <= max(tc_wide$date), ]
  
  # merge tc_wide + spei3 for PDO-phase composite
  spei3$year  <- as.integer(format(spei3$date, "%Y"))
  spei3$month <- as.integer(format(spei3$date, "%m"))
  tc_spei <- merge(tc_wide[, c("date","ONI","PDO","PNA","pdo_phase")],
                   spei3[, c("date","value")],
                   by = "date", all = FALSE)
  names(tc_spei)[names(tc_spei) == "value"] <- "SPEI3"
  
  cat(sprintf("  Merged teleconnection + SPEI-3: %d months\n", nrow(tc_spei)))
  
  # ── drought-category fill for SPEI-3 panel ──────────────────────────────────
  spei3$fill_cat <- dplyr::case_when(
    spei3$value >  0     ~ "Wet/normal",
    spei3$value > -0.5   ~ "Mild deficit",
    spei3$value > -1.0   ~ "Moderate drought",
    spei3$value > -1.5   ~ "Severe drought",
    TRUE                 ~ "Extreme drought"
  )
  spei3$fill_cat <- factor(spei3$fill_cat,
                           levels = c("Wet/normal","Mild deficit",
                                      "Moderate drought","Severe drought","Extreme drought"))
  spei3$date_end <- spei3$date + 31
  
  fill_pal_spei <- c(
    "Wet/normal"       = "#4393c3",
    "Mild deficit"     = "#ffffbf",
    "Moderate drought" = "#fdae61",
    "Severe drought"   = "#f46d43",
    "Extreme drought"  = "#d73027"
  )
  
  # ══════════════════════════════════════════════════════════════════════════════
  # PANEL 1: ONI time series
  # ══════════════════════════════════════════════════════════════════════════════
  p_oni <- ggplot2::ggplot() +
    
    # 2022-2025 band
    ggplot2::annotate("rect",
                      xmin = HIGHLIGHT_START, xmax = HIGHLIGHT_END,
                      ymin = -Inf, ymax = Inf, fill = "grey82", alpha = 0.55) +
    
    # warm ENSO fill (ONI > 0.5)
    ggplot2::geom_rect(
      data = dplyr::filter(tc_wide, ONI >= 0.5),
      ggplot2::aes(xmin = date, xmax = date + 31,
                   ymin = 0.5, ymax = pmax(ONI, 0.5)),
      fill = COL_WARM, alpha = 0.55, inherit.aes = FALSE) +
    
    # cool ENSO fill (ONI < -0.5)
    ggplot2::geom_rect(
      data = dplyr::filter(tc_wide, ONI <= -0.5),
      ggplot2::aes(xmin = date, xmax = date + 31,
                   ymin = pmin(ONI, -0.5), ymax = -0.5),
      fill = COL_COOL, alpha = 0.55, inherit.aes = FALSE) +
    
    ggplot2::geom_hline(yintercept = 0,
                        colour = "grey45", linewidth = 0.35) +
    ggplot2::geom_hline(yintercept = c(-0.5, 0.5),
                        colour = "grey35", linewidth = 0.4, linetype = "dashed") +
    ggplot2::geom_line(
      data = tc_wide,
      ggplot2::aes(x = date, y = ONI),
      colour = COL_ONI, linewidth = 0.55, inherit.aes = FALSE) +
    
    ggplot2::annotate("text", x = as.Date("1951-01-01"), y = 2.5,
                      label = "(a)  ONI", size = 3.0, hjust = 0,
                      fontface = "bold", colour = "grey10") +
    ggplot2::annotate("text", x = as.Date("1997-01-01"), y = 2.8,
                      label = "El Ni\u00f1o\n(ONI \u2265 +0.5)", size = 2.3,
                      colour = "#c0392b", hjust = 0.5) +
    ggplot2::annotate("text", x = as.Date("1973-01-01"), y = -2.3,
                      label = "La Ni\u00f1a\n(ONI \u2264 \u22120.5)", size = 2.3,
                      colour = "#1a5276", hjust = 0.5) +
    
    ggplot2::scale_x_date(date_breaks = "10 years", date_labels = "%Y",
                          expand = ggplot2::expansion(add = c(0, 0))) +
    ggplot2::scale_y_continuous(
      breaks = c(-2, -1, 0, 1, 2),
      expand = ggplot2::expansion(mult = c(0.02, 0.12))) +
    ggplot2::labs(x = NULL, y = "ONI") +
    theme_fig9 +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank())
  
  # ══════════════════════════════════════════════════════════════════════════════
  # PANEL 2: PDO time series
  # ══════════════════════════════════════════════════════════════════════════════
  p_pdo <- ggplot2::ggplot() +
    
    ggplot2::annotate("rect",
                      xmin = HIGHLIGHT_START, xmax = HIGHLIGHT_END,
                      ymin = -Inf, ymax = Inf, fill = "grey82", alpha = 0.55) +
    
    # warm PDO fill
    ggplot2::geom_rect(
      data = dplyr::filter(tc_wide, PDO >= 0.5),
      ggplot2::aes(xmin = date, xmax = date + 31,
                   ymin = 0.5, ymax = pmax(PDO, 0.5)),
      fill = COL_WARM, alpha = 0.55, inherit.aes = FALSE) +
    
    # cool PDO fill
    ggplot2::geom_rect(
      data = dplyr::filter(tc_wide, PDO <= -0.5),
      ggplot2::aes(xmin = date, xmax = date + 31,
                   ymin = pmin(PDO, -0.5), ymax = -0.5),
      fill = COL_COOL, alpha = 0.55, inherit.aes = FALSE) +
    
    ggplot2::geom_hline(yintercept = 0,
                        colour = "grey45", linewidth = 0.35) +
    ggplot2::geom_hline(yintercept = c(-0.5, 0.5),
                        colour = "grey35", linewidth = 0.4, linetype = "dashed") +
    ggplot2::geom_line(
      data = tc_wide,
      ggplot2::aes(x = date, y = PDO),
      colour = COL_PDO, linewidth = 0.55, inherit.aes = FALSE) +
    
    # annotate warm/cool PDO regime shifts
    ggplot2::annotate("segment",
                      x = as.Date("1977-01-01"), xend = as.Date("1998-01-01"),
                      y = 1.6, yend = 1.6, arrow = NULL,
                      colour = "#922b21", linewidth = 0.5) +
    ggplot2::annotate("text",
                      x = as.Date("1987-06-01"), y = 1.75,
                      label = "Warm PDO regime 1977\u20131998",
                      size = 2.1, colour = "#922b21", hjust = 0.5) +
    
    ggplot2::annotate("text",
                      x = as.Date("1951-01-01"), y = -1.5,
                      label = "(b)  PDO", size = 3.0, hjust = 0,
                      fontface = "bold", colour = "grey10") +
    
    ggplot2::scale_x_date(date_breaks = "10 years", date_labels = "%Y",
                          expand = ggplot2::expansion(add = c(0, 0))) +
    ggplot2::scale_y_continuous(
      breaks = c(-2, -1, 0, 1, 2),
      expand = ggplot2::expansion(mult = c(0.05, 0.20))) +
    ggplot2::labs(x = NULL, y = "PDO") +
    theme_fig9 +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank())
  
  # ══════════════════════════════════════════════════════════════════════════════
  # PANEL 3: SPEI-3 with drought-category fill
  # ══════════════════════════════════════════════════════════════════════════════
  ylo_s <- min(-2.6, min(spei3$value, na.rm = TRUE) * 1.08)
  yhi_s <- max( 2.0, max(spei3$value, na.rm = TRUE) * 1.08)
  
  p_spei3 <- ggplot2::ggplot() +
    
    ggplot2::annotate("rect",
                      xmin = HIGHLIGHT_START, xmax = HIGHLIGHT_END,
                      ymin = -Inf, ymax = Inf, fill = "grey82", alpha = 0.55) +
    ggplot2::annotate("text",
                      x = HIGHLIGHT_START + 180, y = yhi_s * 0.92,
                      label = "2022\u20132025",
                      size = 2.3, hjust = 0, colour = "grey35", fontface = "italic") +
    
    # drought-category fill rectangles
    ggplot2::geom_rect(
      data = spei3,
      ggplot2::aes(xmin = date, xmax = date_end,
                   ymin = -Inf, ymax = Inf, fill = fill_cat),
      alpha = 0.28, inherit.aes = FALSE) +
    ggplot2::scale_fill_manual(values = fill_pal_spei,
                               name = "Drought category",
                               guide = ggplot2::guide_legend(
                                 nrow = 1,
                                 title.position = "top")) +
    
    ggplot2::geom_hline(yintercept = 0,
                        colour = "grey45", linewidth = 0.35) +
    ggplot2::geom_hline(yintercept = THR_EVENT,
                        colour = "grey30", linewidth = 0.45, linetype = "dashed") +
    ggplot2::annotate("text",
                      x = as.Date("1951-01-01"), y = THR_EVENT + 0.09,
                      label = "x\u2080 = \u22120.5", size = 2.1,
                      hjust = 0, colour = "grey30") +
    ggplot2::geom_line(
      data = spei3,
      ggplot2::aes(x = date, y = value),
      colour = COL_SPEI, linewidth = 0.55, inherit.aes = FALSE) +
    
    ggplot2::annotate("text", x = as.Date("1951-01-01"), y = yhi_s * 0.92,
                      label = "(c)  SPEI-3", size = 3.0, hjust = 0,
                      fontface = "bold", colour = "grey10") +
    
    ggplot2::scale_x_date(date_breaks = "10 years", date_labels = "%Y",
                          expand = ggplot2::expansion(add = c(0, 0))) +
    ggplot2::scale_y_continuous(
      breaks = c(-2, -1.5, -1, -0.5, 0, 1, 2),
      limits = c(ylo_s, yhi_s),
      expand = ggplot2::expansion(mult = c(0, 0))) +
    ggplot2::labs(x = "Year", y = "SPEI-3") +
    theme_fig9
  
  # ── assemble Fig9a ─────────────────────────────────────────────────────────
  fig9a <- p_oni / p_pdo / p_spei3 +
    patchwork::plot_layout(heights = c(1, 1, 1.2), guides = "collect") +
    patchwork::plot_annotation(
      title    = "Teleconnection indices and basin SPEI-3 \u2014 Nechako River Basin (1950\u20132025)",
      subtitle = paste0(
        "Red/blue shading: warm (ONI/PDO \u2265 +0.5) / cool (ONI/PDO \u2264 \u22120.5) phases.  ",
        "Grey band: 2022\u20132025 focus period.  Dashed: phase boundaries (\u00b10.5)."),
      theme = ggplot2::theme(
        plot.title    = ggplot2::element_text(size = 10, face = "bold", hjust = 0.5),
        plot.subtitle = ggplot2::element_text(size = 7.5, colour = "grey40", hjust = 0),
        legend.position = "bottom",
        legend.text  = ggplot2::element_text(size = 7),
        legend.key.size = ggplot2::unit(0.35, "cm")
      )
    )
  
  fig9a_pdf <- file.path(TELE_DIR, "Fig9a_teleconnection_timeseries.pdf")
  fig9a_png <- file.path(TELE_DIR, "Fig9a_teleconnection_timeseries.png")
  tryCatch({
    ggplot2::ggsave(fig9a_pdf, fig9a, width = 7.5, height = 8.5,
                    units = "in", device = "pdf")
    cat(sprintf("  \u2713 Fig9a (PDF): %s\n", basename(fig9a_pdf)))
  }, error = function(e) cat(sprintf("  \u26a0 Fig9a PDF: %s\n", e$message)))
  tryCatch({
    ggplot2::ggsave(fig9a_png, fig9a, width = 7.5, height = 8.5,
                    units = "in", dpi = 300, device = "png")
    cat(sprintf("  \u2713 Fig9a (PNG): %s\n", basename(fig9a_png)))
  }, error = function(e) cat(sprintf("  \u26a0 Fig9a PNG: %s\n", e$message)))
  
  
  # ══════════════════════════════════════════════════════════════════════════════
  # FIGURE 9b — PDO-phase composite (3 panels)  [SECTION 5.3 SUPPORT]
  #
  # Root causes of previous blank figure:
  #   (1) SPEI-3 load failures silently caught → tc_spei had 0 rows
  #   (2) phase_stats = NULL → ggplot(NULL) = blank device
  #   (3) annotate(x = 0.55) on discrete axis → annotation outside plot area
  #
  # Fix strategy:
  #   • 3-tier SPEI-3 loader with direct seasonal CSV fallback
  #   • Hard validation of tc_spei before any plotting
  #   • All x-references use factor levels, not numeric coordinates
  #
  # Three panels for Section 5.3:
  #   (a) Mean SPEI-3 ± boot CI per PDO phase  — magnitude of drought signal
  #   (b) Drought frequency (% months SPEI-3 < -0.5) per phase — occurrence rate
  #   (c) Full SPEI-3 distribution per phase (violin + box) — shape of signal
  # ══════════════════════════════════════════════════════════════════════════════
  
  # ── Tier-3 SPEI-3 loader: direct per-month seasonal CSVs ────────────────────
  # Called only if load_basin_avg_csv() and the basin_average.csv flat file
  # both fail (e.g. w1 has not been run but 3SPEI_ERALand.R has).
  load_spei3_from_seasonal_csvs <- function(seas_dir) {
    MONTH_ABB <- c("Jan","Feb","Mar","Apr","May","Jun",
                   "Jul","Aug","Sep","Oct","Nov","Dec")
    rows <- vector("list", 12L)
    for (m in seq_len(12L)) {
      csv_f <- file.path(seas_dir,
                         sprintf("spei_03_month%02d_%s.csv", m, MONTH_ABB[m]))
      if (!file.exists(csv_f)) next
      df <- tryCatch(data.table::fread(csv_f, data.table = FALSE),
                     error = function(e) NULL)
      if (is.null(df)) next
      yr_cols <- setdiff(names(df), c("lon","lat"))
      yr_cols <- yr_cols[grepl("^[0-9]{4}$", yr_cols)]
      if (!length(yr_cols)) next
      yrs  <- as.integer(yr_cols)
      vals <- colMeans(df[, yr_cols, drop = FALSE], na.rm = TRUE)
      rows[[m]] <- data.frame(
        date  = as.Date(paste(yrs, m, "01", sep = "-")),
        value = as.numeric(vals),
        stringsAsFactors = FALSE)
    }
    valid <- Filter(Negate(is.null), rows)
    if (!length(valid)) return(NULL)
    out <- do.call(rbind, valid)
    out <- out[order(out$date), ]
    out[!is.na(out$value), ]
  }
  
  # ── Load SPEI-3 with 3 fallback tiers ────────────────────────────────────────
  spei3_raw <- tryCatch(load_basin_avg_csv("spei", 3), error = function(e) NULL)
  
  if (is.null(spei3_raw) || !nrow(spei3_raw)) {
    cat("  Tier 1 failed → trying flat basin-average CSV...\n")
    spei3_raw <- tryCatch({
      f <- file.path(WD_PATH, "temporal_drought", "basin_averaged_timeseries",
                     "spei_03_basin_average.csv")
      if (!file.exists(f)) stop("not found")
      df <- data.table::fread(f)
      df$date <- as.Date(df$date); df$value <- as.numeric(df$value)
      as.data.frame(df)
    }, error = function(e) NULL)
  }
  
  if (is.null(spei3_raw) || !nrow(spei3_raw)) {
    cat("  Tier 2 failed → reading per-month seasonal CSVs directly...\n")
    seas_dir_pm <- file.path(WD_PATH, "spei_results_seasonal")
    spei3_raw   <- load_spei3_from_seasonal_csvs(seas_dir_pm)
    if (!is.null(spei3_raw) && nrow(spei3_raw))
      cat(sprintf("  Tier 3: loaded %d months from seasonal CSVs\n", nrow(spei3_raw)))
  }
  
  if (is.null(spei3_raw) || !nrow(spei3_raw))
    stop("SPEI-3 unavailable after 3 loading tiers. Run 3SPEI_ERALand.R first.")
  
  spei3 <- as.data.frame(spei3_raw)
  nm <- names(spei3)
  names(spei3)[which(nm %in% c("date","Date"))[1]]             <- "date"
  names(spei3)[which(nm %in% c("value","Value","Mean_Value"))[1]] <- "value"
  spei3$date  <- as.Date(spei3$date)
  spei3$value <- as.numeric(spei3$value)
  spei3 <- spei3[!is.na(spei3$date) & !is.na(spei3$value), ]
  
  # ── Merge with teleconnection data ───────────────────────────────────────────
  tc_sub <- tc_wide[, c("date","ONI","PDO","PNA","pdo_phase")]
  tc_sub$date <- as.Date(tc_sub$date)
  tc_spei <- merge(tc_sub, spei3[, c("date","value")], by = "date", all = FALSE)
  names(tc_spei)[names(tc_spei) == "value"] <- "SPEI3"
  tc_spei <- tc_spei[!is.na(tc_spei$SPEI3) & !is.na(tc_spei$pdo_phase), ]
  cat(sprintf("  Merged teleconnection + SPEI-3: %d months\n", nrow(tc_spei)))
  
  if (!nrow(tc_spei))
    stop("Merge produced 0 rows: check date formats in tc_wide and SPEI-3 loader.")
  
  # ── Bootstrap statistics ──────────────────────────────────────────────────────
  set.seed(42L)
  N_BOOT       <- 2000L
  PHASE_LEVELS <- c("Cool", "Neutral", "Warm")
  PHASE_COLS   <- c(Cool = "#2980B9", Neutral = "#888888", Warm = "#C0392B")
  
  tc_spei$pdo_phase_f <- factor(tc_spei$pdo_phase, levels = PHASE_LEVELS)
  
  phase_stats_list <- lapply(PHASE_LEVELS, function(ph) {
    x <- tc_spei$SPEI3[tc_spei$pdo_phase_f == ph]
    x <- x[!is.na(x)]
    if (length(x) < 5L) {
      cat(sprintf("  ⚠ Phase '%s' has only %d observations — skipping\n",
                  ph, length(x)))
      return(NULL)
    }
    boot_means <- replicate(N_BOOT, mean(sample(x, length(x), replace = TRUE)))
    n_drought  <- sum(x < -0.5)
    data.frame(
      pdo_phase   = ph,
      N           = length(x),
      mean_spei   = mean(x),
      ci_lo       = quantile(boot_means, 0.025),
      ci_hi       = quantile(boot_means, 0.975),
      pct_drought = 100 * n_drought / length(x),
      stringsAsFactors = FALSE)
  })
  
  # Hard stop if all phases failed
  phase_stats <- do.call(rbind, Filter(Negate(is.null), phase_stats_list))
  if (is.null(phase_stats) || !nrow(phase_stats))
    stop("phase_stats is empty — all PDO phases have < 5 observations.")
  
  phase_stats$pdo_phase_f <- factor(phase_stats$pdo_phase, levels = PHASE_LEVELS)
  phase_stats$n_label     <- sprintf("n=%d", phase_stats$N)
  phase_stats$mean_label  <- sprintf("%.3f", phase_stats$mean_spei)
  phase_stats$pct_label   <- sprintf("%.0f%%", phase_stats$pct_drought)
  
  cat(sprintf("  Phase statistics computed for %d phases\n", nrow(phase_stats)))
  print(phase_stats[, c("pdo_phase","N","mean_spei","ci_lo","ci_hi","pct_drought")])
  
  # ── Shared theme ──────────────────────────────────────────────────────────────
  theme_9b <- ggplot2::theme_classic(base_size = 9) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(size = 9,  face = "bold", hjust = 0),
      plot.subtitle    = ggplot2::element_text(size = 7.5, colour = "grey35",
                                               hjust = 0, margin = ggplot2::margin(b=3)),
      axis.text        = ggplot2::element_text(size = 8),
      axis.title       = ggplot2::element_text(size = 8.5),
      panel.grid.major.y = ggplot2::element_line(colour = "grey93", linewidth = 0.3),
      panel.grid.major.x = ggplot2::element_blank(),
      plot.margin      = ggplot2::margin(3, 6, 3, 4))
  
  # ── Panel (a): Mean SPEI-3 ± boot CI ─────────────────────────────────────────
  pa_9b <- ggplot2::ggplot(phase_stats,
                           ggplot2::aes(x = pdo_phase_f, y = mean_spei,
                                        fill = pdo_phase_f)) +
    ggplot2::geom_hline(yintercept = 0,    colour = "grey40", linewidth = 0.45) +
    ggplot2::geom_hline(yintercept = -0.5, colour = "grey40", linewidth = 0.45,
                        linetype = "dashed") +
    # violin behind bars (full distribution)
    ggplot2::geom_violin(
      data   = tc_spei,
      ggplot2::aes(x = pdo_phase_f, y = SPEI3, fill = pdo_phase_f),
      alpha  = 0.10, colour = NA, width = 0.75,
      inherit.aes = FALSE) +
    ggplot2::geom_col(width = 0.42, colour = "white", linewidth = 0.25,
                      alpha = 0.90) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = ci_lo, ymax = ci_hi),
                           width = 0.16, linewidth = 0.65, colour = "grey15") +
    # n= label above CI bar — uses phase name on x (discrete), NOT 0.55
    ggplot2::geom_text(ggplot2::aes(y = ci_hi + 0.015, label = n_label),
                       size = 2.7, colour = "grey20", vjust = 0) +
    # mean value inside bar
    ggplot2::geom_text(ggplot2::aes(y = mean_spei / 2, label = mean_label),
                       size = 2.9, fontface = "bold", colour = "white") +
    # threshold label anchored to the "Cool" bar (first factor level)
    ggplot2::annotate("text",
                      x     = "Cool",
                      y     = -0.5 + 0.02,
                      label = "x\u2080 = \u22120.5",
                      size  = 2.4, hjust = 0.1, vjust = 0, colour = "grey30") +
    ggplot2::scale_fill_manual(values = PHASE_COLS, guide = "none") +
    ggplot2::scale_y_continuous(
      breaks = seq(-0.8, 0.5, by = 0.1),
      expand = ggplot2::expansion(mult = c(0.12, 0.20))) +
    ggplot2::labs(
      title    = "(a)  Mean SPEI-3 by PDO phase",
      subtitle = paste0("Bars = mean \u00b1 bootstrapped 95% CI (", N_BOOT, " reps).  ",
                        "Violin = full distribution.  Dashed = event threshold."),
      x = "PDO phase", y = "Mean basin SPEI-3") +
    theme_9b
  
  # ── Panel (b): Drought frequency (% months < -0.5) per PDO phase ─────────────
  pb_9b <- ggplot2::ggplot(phase_stats,
                           ggplot2::aes(x = pdo_phase_f, y = pct_drought,
                                        fill = pdo_phase_f)) +
    ggplot2::geom_col(width = 0.52, colour = "white", linewidth = 0.25,
                      alpha = 0.88) +
    ggplot2::geom_text(ggplot2::aes(y = pct_drought + 0.8, label = pct_label),
                       size = 3.0, fontface = "bold", colour = "grey20", vjust = 0) +
    ggplot2::geom_text(ggplot2::aes(y = pct_drought / 2, label = n_label),
                       size = 2.6, colour = "white", vjust = 0.5) +
    ggplot2::scale_fill_manual(values = PHASE_COLS, guide = "none") +
    ggplot2::scale_y_continuous(
      limits = c(0, max(phase_stats$pct_drought, na.rm = TRUE) * 1.25),
      expand = ggplot2::expansion(mult = c(0, 0.05)),
      labels = function(x) paste0(x, "%")) +
    ggplot2::labs(
      title    = "(b)  Drought frequency by PDO phase",
      subtitle = "% of months with SPEI-3 < \u22120.5 within each PDO phase.",
      x = "PDO phase", y = "Drought months (%)") +
    theme_9b
  
  # ── Panel (c): Full SPEI-3 distribution per phase ────────────────────────────
  pc_9b <- ggplot2::ggplot(tc_spei,
                           ggplot2::aes(x = pdo_phase_f, y = SPEI3,
                                        fill = pdo_phase_f, colour = pdo_phase_f)) +
    ggplot2::geom_hline(yintercept = 0,    colour = "grey40", linewidth = 0.4) +
    ggplot2::geom_hline(yintercept = -0.5, colour = "grey40", linewidth = 0.4,
                        linetype = "dashed") +
    ggplot2::geom_violin(alpha = 0.18, colour = NA, width = 0.75, trim = TRUE) +
    ggplot2::geom_boxplot(width = 0.22, outlier.size = 0.7, outlier.alpha = 0.4,
                          fill = "white", linewidth = 0.5) +
    ggplot2::scale_fill_manual(values   = PHASE_COLS, guide = "none") +
    ggplot2::scale_colour_manual(values = PHASE_COLS, guide = "none") +
    ggplot2::scale_y_continuous(breaks = seq(-3, 3, by = 0.5),
                                expand = ggplot2::expansion(mult = c(0.05, 0.05))) +
    ggplot2::labs(
      title    = "(c)  SPEI-3 distribution by PDO phase",
      subtitle = "Violin + box-and-whisker; dashed = event threshold (x\u2080 = \u22120.5).",
      x = "PDO phase", y = "SPEI-3") +
    theme_9b
  
  # ── Assemble 3-panel figure ───────────────────────────────────────────────────
  period_str <- sprintf("%d\u2013%d",
                        as.integer(format(min(tc_spei$date),"%Y")),
                        as.integer(format(max(tc_spei$date),"%Y")))
  
  fig9b <- pa_9b | pb_9b | pc_9b +
    patchwork::plot_annotation(
      title    = paste0("Basin SPEI-3 conditioned on PDO phase \u2014",
                        " Nechako River Basin  (", period_str, ")"),
      subtitle = paste0(
        "PDO phase thresholds: Cool \u2264 \u22120.5 < Neutral \u2264 +0.5 < Warm.",
        "  Supports Section 5.3: Warm PDO amplifies drought frequency and severity."),
      theme = ggplot2::theme(
        plot.title    = ggplot2::element_text(size = 10, face = "bold", hjust = 0.5),
        plot.subtitle = ggplot2::element_text(size = 8, colour = "grey35", hjust = 0,
                                              margin = ggplot2::margin(b = 5))))
  
  fig9b_pdf <- file.path(TELE_DIR, "Fig9b_PDO_phase_SPEI3_composite.pdf")
  fig9b_png <- file.path(TELE_DIR, "Fig9b_PDO_phase_SPEI3_composite.png")
  tryCatch({
    ggplot2::ggsave(fig9b_pdf, fig9b, width = 10.0, height = 5.0,
                    units = "in", device = "pdf")
    cat(sprintf("  \u2713 Fig9b (PDF): %s\n", basename(fig9b_pdf)))
  }, error = function(e) cat(sprintf("  \u26a0 Fig9b PDF: %s\n", e$message)))
  tryCatch({
    ggplot2::ggsave(fig9b_png, fig9b, width = 10.0, height = 5.0,
                    units = "in", dpi = 300, device = "png")
    cat(sprintf("  \u2713 Fig9b (PNG): %s\n", basename(fig9b_png)))
  }, error = function(e) cat(sprintf("  \u26a0 Fig9b PNG: %s\n", e$message)))
  
  # save phase statistics CSV for Table in manuscript
  write.csv(
    phase_stats[, c("pdo_phase","N","mean_spei","ci_lo","ci_hi","pct_drought")],
    file.path(TELE_DIR, "PDO_phase_SPEI3_composite_stats.csv"),
    row.names = FALSE)
  cat(sprintf("  \u2713 Phase stats: PDO_phase_SPEI3_composite_stats.csv\n"))
  
  
}, error = function(e) {
  cat(sprintf("  \u26a0 Step 7 failed: %s\n  Check that w1/w2 have been run first.\n",
              e$message))
})

####################################################################################
# DONE
####################################################################################
cat("\n╔═══════════════════════════════════════════════════╗\n")
cat("║  w7 COMPLETE                                      ║\n")
cat("╚═══════════════════════════════════════════════════╝\n")
cat("  Key files written to:", TELE_DIR, "\n")
cat("  Existing outputs: QC plots, GAM PDF, CSVs, XLSX\n")
cat("  NEW manuscript figures:\n")
cat("    Fig9a_teleconnection_timeseries.pdf/.png\n")
cat("    Fig9b_PDO_phase_SPEI3_composite.pdf/.png\n")
cat("    PDO_phase_SPEI3_composite_stats.csv\n")
cat("  Next: run w8_moderated_mediation_gam.R\n\n")