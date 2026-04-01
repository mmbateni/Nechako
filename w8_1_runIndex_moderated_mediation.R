####################################################################################
# w8_01_run_index.R
# SINGLE-INDEX WORKER: all lags for ONE drought index (GAM + MM + JN)
#
# USAGE — two modes:
#
#   A) Interactive / debug a single index:
#        source("w8_00_setup.R")
#        SCRIPT_START <- proc.time()
#        idx_spec <- list(index="spei", scale=6)
#        source("w8_01_run_index.R")
#
#   B) Called automatically by w8_02_batch.R (parallel worker).
#        Do NOT source directly when running the full batch.
#
# INPUTS  (must exist in calling environment):
#   idx_spec, SCRIPT_START, and all globals from w8_00_setup.R
#
# RETURN VALUE (used by w8_02_batch.R):
#   .idx_result  list(index, status, optimal_lag, elapsed_sec)
####################################################################################

if (exists("WD_PATH")) setwd(WD_PATH)
idx_spec <- list(index='spei', scale=6)
if (!exists("TELE_MASTER")) source("DROUGHT_ANALYSIS_utils.R")
suppressPackageStartupMessages({ library(boot); library(parallel) })
if (!exists("idx_spec"))
  stop("idx_spec not found. Set it before sourcing: idx_spec <- list(index='spei', scale=6)",
       call.=FALSE)
if (!exists("SCRIPT_START"))
  SCRIPT_START <- proc.time()   # fallback: measure from this script's own start
index_label   <- sprintf("%s-%d", toupper(idx_spec$index), idx_spec$scale)
IDX_DIR       <- file.path(MM_DIR, index_label)
dir.create(IDX_DIR, showWarnings=FALSE, recursive=TRUE)

idx_t0        <- proc.time()
lag_sec_vec__ <- numeric(0L)

all_mm_results <- list(); all_gam_fits <- list()
mm_lag_profile <- list(); gam_lag_prof  <- list(); convergence_rows <- list()

cat(sprintf("\n%s\n  INDEX: %s | Lags 0-%d | script elapsed %s\n%s\n\n",
            strrep("=",62), index_label, max(LAGS_TO_TEST),
            fmt_dur(elapsed_sec(SCRIPT_START)), strrep("=",62)))

####################################################################################
# INNER LAG LOOP
####################################################################################

for (lag in LAGS_TO_TEST) {
  
  lag_t0__  <- proc.time()
  lag_pos__ <- which(LAGS_TO_TEST == lag)
  cat(sprintf("\n>> LAG = %d months  [%d/%d | index elapsed %s | script elapsed %s]\n",
              lag, lag_pos__, length(LAGS_TO_TEST),
              fmt_dur(elapsed_sec(idx_t0)), fmt_dur(elapsed_sec(SCRIPT_START))))
  
  df_lag <- tryCatch(
    build_analysis_dataframe(idx_spec$index, idx_spec$scale,
                             lag_months=lag, start_year=1950L, end_year=2024L,
                             tele_df=TELE_MASTER),
    error=function(e) { cat("  XX Data:", e$message, "\n"); NULL })
  if (is.null(df_lag)) next
  if (lag == 0L) check_mm_assumptions(df_lag)
  
  # ── TRACK 1: GAM ───────────────────────────────────────────────────────────────
  cat(sprintf("\n  --- TRACK 1: GAM (lag=%d) ---\n", lag))
  gam_fit <- fit_gam_full(df_lag, response_var=GAM_RESPONSE,
                          k_main=GAM_K_MAIN, k_tensor=GAM_K_TENSOR,
                          include_seasonality=TRUE, select=TRUE)
  all_gam_fits[[sprintf("lag_%02d",lag)]] <- gam_fit
  gam_s <- if (!is.null(gam_fit)) summary(gam_fit) else NULL   # PERF-3
  
  if (!is.null(gam_fit)) {
    vp <- gam_variance_partition(gam_fit, precomputed_s=gam_s)
    if (!is.null(vp)) {
      vp$lag <- lag; vp$index <- index_label
      write.csv(vp, file.path(IDX_DIR,sprintf("gam_variance_partition_lag%02d.csv",lag)),
                row.names=FALSE)
    }
    gam_concurvity_check(gam_fit)
    gam_lag_prof[[length(gam_lag_prof)+1]] <-
      gam_lag_summary(gam_fit, lag, index_label, precomputed_s=gam_s)
    
    gam_pdf <- file.path(IDX_DIR, sprintf("gam_lag%02d_%s.pdf",lag,index_label))
    if (safe_pdf(gam_pdf, width=14, height=10)) {
      p_surf  <- plot_gam_surfaces(gam_fit, index_label=index_label, lag=lag,
                                   drought_thr=DROUGHT_THRESHOLD)
      p_part  <- plot_gam_partials(gam_fit,  index_label=index_label, lag=lag)
      p_resid <- plot_gam_residuals(gam_fit, index_label=index_label)
      if (!is.null(p_surf))  print(p_surf)
      if (!is.null(p_part))  print(p_part)
      if (!is.null(p_resid)) print(p_resid)
      grDevices::dev.off()
      cat(sprintf("  OK Saved: gam_lag%02d_%s.pdf\n",lag,index_label))
    }
  }
  
  # ── TRACK 2: MODERATED MEDIATION  (PERF-4 screening) ──────────────────────────
  cat(sprintf("\n  --- TRACK 2: Moderated Mediation (lag=%d) ---\n", lag))
  gam_ti_p_screen <- NA_real_
  if (!is.null(gam_s)) {
    smt  <- gam_s$s.table; ti_i <- grep("ti\\(ONI,PDO",rownames(smt))
    if (length(ti_i)) gam_ti_p_screen <- smt[ti_i[1],"p-value"]
  }
  n_boot_this <- if (is.na(gam_ti_p_screen) || gam_ti_p_screen <= GAM_SCREEN_P) {
    cat(sprintf("  [PERF-4] ti(ONI,PDO) p=%.3f -> full %d-rep bootstrap\n",
                ifelse(is.na(gam_ti_p_screen),NaN,gam_ti_p_screen), N_BOOTSTRAP))
    N_BOOTSTRAP
  } else {
    cat(sprintf("  [PERF-4] ti(ONI,PDO) p=%.3f > %.2f -> reduced %d-rep bootstrap\n",
                gam_ti_p_screen, GAM_SCREEN_P, N_BOOT_SCREEN))
    N_BOOT_SCREEN
  }
  
  mm_res <- fit_moderated_mediation(df_lag, lag=lag, n_boot=n_boot_this,
                                    seed=BOOTSTRAP_SEED, block_length=BLOCK_LENGTH)
  all_mm_results[[sprintf("lag_%02d",lag)]] <- mm_res
  
  if (!is.null(mm_res)) {
    pe_lag <- mm_res$params
    for (ie_lbl in c("IE.low","IE.mid","IE.high","IMM")) {
      row <- pe_lag[pe_lag$label==ie_lbl,]
      if (nrow(row))
        mm_lag_profile[[length(mm_lag_profile)+1]] <- data.frame(
          lag=lag, param=ie_lbl, est=.get_est(row),
          ci_lower=if("ci.lower"%in%names(row))row$ci.lower[1] else NA_real_,
          ci_upper=if("ci.upper"%in%names(row))row$ci.upper[1] else NA_real_,
          p_value =if("pvalue"  %in%names(row))row$pvalue[1]   else NA_real_,
          stringsAsFactors=FALSE)
    }
    mm_pdf <- file.path(IDX_DIR, sprintf("mm_lag%02d_%s.pdf",lag,index_label))
    if (safe_pdf(mm_pdf, width=13, height=9)) {
      print(plot_path_diagram(mm_res, lag, index_label))
      print(plot_conditional_ie(mm_res, lag, index_label))
      jn_df <- tryCatch(compute_jn_floodlight(df_lag, n_boot=JN_BOOT_N),error=function(e)NULL)
      if (!is.null(jn_df)&&nrow(jn_df)) {
        print(plot_jn_floodlight(jn_df, lag, index_label))
        write.csv(jn_df, file.path(IDX_DIR,sprintf("jn_data_lag%02d.csv",lag)),row.names=FALSE)
      }
      grDevices::dev.off()
      cat(sprintf("  OK Saved: mm_lag%02d_%s.pdf\n",lag,index_label))
    }
  }
  
  # ── TRACK 3: CONVERGENCE CHECK  (reuses gam_s — PERF-3) ───────────────────────
  gam_ti_p<-NA_real_; mm_imm_p<-NA_real_; gam_dev<-NA_real_; mm_imm_e<-NA_real_
  if (!is.null(gam_s)) {
    smt<-gam_s$s.table; ti_idx<-grep("ti\\(ONI,PDO",rownames(smt))
    if (length(ti_idx)) gam_ti_p <- round(smt[ti_idx[1],"p-value"],4)
    gam_dev <- round(gam_s$dev.expl*100,2)
  }
  if (!is.null(mm_res)) {
    imm_row <- mm_res$params[mm_res$params$label=="IMM",]
    if (nrow(imm_row)) { mm_imm_p<-round(imm_row$pvalue[1],4); mm_imm_e<-round(.get_est(imm_row),4) }
  }
  convergence_rows[[length(convergence_rows)+1]] <- data.frame(
    Lag=lag, GAM_Dev_Expl_pct=gam_dev, GAM_ti_ONI_PDO_p=gam_ti_p,
    GAM_ti_sig=!is.na(gam_ti_p)&&gam_ti_p<0.05,
    MM_IMM_est=mm_imm_e, MM_IMM_p=mm_imm_p,
    MM_IMM_sig=!is.na(mm_imm_p)&&mm_imm_p<0.05,
    Both_agree=(!is.na(gam_ti_p)&&gam_ti_p<0.05)&&(!is.na(mm_imm_p)&&mm_imm_p<0.05),
    Boot_reps=n_boot_this, Block_length=BLOCK_LENGTH, stringsAsFactors=FALSE)
  
  # ── Lag timing & progress bar ──────────────────────────────────────────────────
  lag_sec_vec__ <- c(lag_sec_vec__, elapsed_sec(lag_t0__))
  lags_done__   <- length(lag_sec_vec__)
  eta_lag__     <- mean(lag_sec_vec__) * (length(LAGS_TO_TEST) - lags_done__)
  cat(fmt_progress(lags_done__, length(LAGS_TO_TEST), elapsed_sec(idx_t0), eta_lag__), "\n")
  
}  # end lag loop

####################################################################################
# GAM LAG PROFILE
####################################################################################
cat(sprintf("\n-- GAM lag profile  [index elapsed %s] --\n", fmt_dur(elapsed_sec(idx_t0))))

if (length(gam_lag_prof) > 0) {
  gam_lp_df <- do.call(rbind, gam_lag_prof)
  write.csv(gam_lp_df, file.path(IDX_DIR,"gam_lag_profile.csv"), row.names=FALSE)
  cat("  Saved: gam_lag_profile.csv\n")
  
  lp_long <- tidyr::pivot_longer(gam_lp_df,cols=c("p_ti_ONI_PDO","p_ti_ONI_PNA","p_ti_PDO_PNA"),
                                 names_to="Term",values_to="p_value")
  lp_long$Term <- dplyr::recode(lp_long$Term,
                                p_ti_ONI_PDO="ti(ONI,PDO)  <- key",
                                p_ti_ONI_PNA="ti(ONI,PNA)", p_ti_PDO_PNA="ti(PDO,PNA)")
  ti_clrs <- c("ti(ONI,PDO)  <- key"="#C0392B","ti(ONI,PNA)"="#2980B9","ti(PDO,PNA)"="#27AE60")
  
  pA <- ggplot2::ggplot(gam_lp_df,ggplot2::aes(x=Lag,y=Dev_Expl_pct))+
    ggplot2::geom_line(colour="#1A5276",linewidth=1.2)+ggplot2::geom_point(size=3,colour="#1A5276")+
    ggplot2::scale_x_continuous(breaks=LAGS_TO_TEST)+shared_ts_theme(11)+
    ggplot2::labs(title="GAM: Deviance explained across lags",
                  subtitle=sprintf("%s | Main+tensor smooths + seasonality | select=TRUE",index_label),
                  x="Lag (months)",y="Deviance explained (%)")
  pB <- ggplot2::ggplot(lp_long[!is.na(lp_long$p_value),],ggplot2::aes(x=Lag,y=p_value,colour=Term))+
    ggplot2::geom_hline(yintercept=0.05,linetype="dashed",colour="grey30")+
    ggplot2::annotate("text",x=max(LAGS_TO_TEST),y=0.05,label="p=0.05",hjust=1,vjust=-0.4,size=3,colour="grey30")+
    ggplot2::geom_line(linewidth=1.0)+ggplot2::geom_point(size=2.5)+
    ggplot2::scale_colour_manual(values=ti_clrs)+ggplot2::scale_x_continuous(breaks=LAGS_TO_TEST)+
    ggplot2::scale_y_log10(breaks=c(0.001,0.01,0.05,0.1,0.5,1))+shared_ts_theme(11)+
    ggplot2::labs(title="GAM tensor interaction p-values across lags",
                  subtitle="Log10 scale | Points below dashed line are significant (p<0.05)",
                  x="Lag (months)",y="p-value (log10)")
  
  gam_prof_pdf <- file.path(IDX_DIR,sprintf("gam_lag_profile_%s.pdf",index_label))
  if (safe_pdf(gam_prof_pdf,width=13,height=11)) {
    print(pA/pB); grDevices::dev.off()
    cat(sprintf("  Saved: gam_lag_profile_%s.pdf\n",index_label))
  }
}

####################################################################################
# MM LAG PROFILE
####################################################################################
cat("\n-- MM lag profile --\n")

if (length(mm_lag_profile) > 0) {
  lp_df <- do.call(rbind, mm_lag_profile)
  lp_df$param_label <- dplyr::case_when(
    lp_df$param=="IE.low"  ~ W_LABELS[1], lp_df$param=="IE.mid"  ~ W_LABELS[2],
    lp_df$param=="IE.high" ~ W_LABELS[3], lp_df$param=="IMM"     ~ "IMM",
    TRUE ~ lp_df$param)
  
  p_ie <- ggplot2::ggplot(lp_df[lp_df$param!="IMM",],
                          ggplot2::aes(x=lag,y=est,colour=param_label,fill=param_label))+
    ggplot2::geom_hline(yintercept=0,linetype="dashed",colour="grey40")+
    ggplot2::geom_ribbon(ggplot2::aes(ymin=ci_lower,ymax=ci_upper),alpha=0.15,colour=NA)+
    ggplot2::geom_line(linewidth=1.1)+
    ggplot2::geom_point(ggplot2::aes(shape=p_value<0.05),size=3)+
    ggplot2::scale_colour_manual(values=c("#1A5276","#2980B9","#C0392B"),name="PDO phase")+
    ggplot2::scale_fill_manual(values=c("#1A5276","#2980B9","#C0392B"),name="PDO phase")+
    ggplot2::scale_shape_manual(values=c("FALSE"=1,"TRUE"=19),name="p<0.05")+
    ggplot2::scale_x_continuous(breaks=LAGS_TO_TEST)+shared_ts_theme(12)+
    ggplot2::labs(title=sprintf("MM: Conditional Indirect Effects across lags | %s",index_label),
                  subtitle=sprintf("Shading=95%% Block-bootstrap CI (l=%d) | Filled=p<0.05",BLOCK_LENGTH),
                  x="Lag (months)",y="Conditional IE (standardised)")
  p_imm <- ggplot2::ggplot(lp_df[lp_df$param=="IMM",],ggplot2::aes(x=lag,y=est))+
    ggplot2::geom_hline(yintercept=0,linetype="dashed")+
    ggplot2::geom_ribbon(ggplot2::aes(ymin=ci_lower,ymax=ci_upper),fill="#F39C12",alpha=0.25)+
    ggplot2::geom_line(colour="#F39C12",linewidth=1.1)+
    ggplot2::geom_point(ggplot2::aes(shape=p_value<0.05),colour="#F39C12",size=3)+
    ggplot2::scale_shape_manual(values=c("FALSE"=1,"TRUE"=19),name="p<0.05")+
    ggplot2::scale_x_continuous(breaks=LAGS_TO_TEST)+shared_ts_theme(12)+
    ggplot2::labs(title="MM: Index of Moderated Mediation (IMM) across lags",
                  subtitle="IMM!=0 -> PDO significantly moderates the ONI->PNA->Drought indirect pathway",
                  x="Lag (months)",y="IMM (standardised)")
  
  mm_prof_pdf <- file.path(IDX_DIR,sprintf("mm_lag_profile_%s.pdf",index_label))
  if (safe_pdf(mm_prof_pdf,width=13,height=11)) {
    print(p_ie/p_imm); grDevices::dev.off()
    cat(sprintf("  Saved: mm_lag_profile_%s.pdf\n",index_label))
  }
  write.csv(lp_df, file.path(IDX_DIR,sprintf("lag_profile_data_%s.csv",index_label)),row.names=FALSE)
}

####################################################################################
# CONVERGENCE TABLE
####################################################################################
cat("\n-- Convergence table (GAM vs MM) --\n")

if (length(convergence_rows) > 0) {
  conv_df <- do.call(rbind, convergence_rows)
  conv_df$Interpretation <- dplyr::case_when(
    conv_df$Both_agree                       ~ "OK CONVERGENT: both methods agree",
    conv_df$GAM_ti_sig & !conv_df$MM_IMM_sig ~ "GAM detects interaction; MM pathway not significant",
    !conv_df$GAM_ti_sig & conv_df$MM_IMM_sig ~ "MM detects mediation; GAM interaction marginal",
    TRUE                                     ~ "Neither significant at this lag")
  cat("\n"); print(conv_df, row.names=FALSE)
  write.csv(conv_df, file.path(IDX_DIR,sprintf("convergence_table_%s.csv",index_label)),row.names=FALSE)
  cat(sprintf("  Saved: convergence_table_%s.csv\n",index_label))
}

####################################################################################
# SENSITIVITY ANALYSIS
####################################################################################
cat("\n-- Sensitivity analysis --\n")

if (length(gam_lag_prof) > 0) {
  gam_lp_df2 <- do.call(rbind, gam_lag_prof)
  optimal_lag <- gam_lp_df2$Lag[which.max(gam_lp_df2$Dev_Expl_pct)]
} else optimal_lag <- 2L
cat(sprintf("  Optimal lag (max GAM dev.expl): %d months\n", optimal_lag))

sens_mm <- list(); sens_gam <- list()
for (si in SENSITIVITY_INDICES) {
  lbl <- sprintf("%s-%d", toupper(si$index), si$scale)
  cat(sprintf("\n  Sensitivity: %s at lag=%d\n", lbl, optimal_lag))
  df_s <- tryCatch(
    build_analysis_dataframe(si$index, si$scale, lag_months=optimal_lag,
                             start_year=1950L, end_year=2024L, tele_df=TELE_MASTER),
    error=function(e) { cat("  XX", e$message, "\n"); NULL })
  if (is.null(df_s)) next
  gf   <- fit_gam_full(df_s, response_var=GAM_RESPONSE, k_main=GAM_K_MAIN,
                       k_tensor=GAM_K_TENSOR, include_seasonality=TRUE, select=TRUE)
  gf_s <- if (!is.null(gf)) summary(gf) else NULL
  sens_gam[[length(sens_gam)+1]] <- gam_lag_summary(gf, optimal_lag, lbl, precomputed_s=gf_s)
  mm_s <- fit_moderated_mediation(df_s, lag=optimal_lag, n_boot=N_BOOTSTRAP,
                                  seed=BOOTSTRAP_SEED, block_length=BLOCK_LENGTH)
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
  write.csv(sg, file.path(MM_DIR,sprintf("gam_sensitivity_%s.csv",index_label)),row.names=FALSE)
  cat("\n  GAM sensitivity summary:\n"); print(sg, row.names=FALSE)
}
if (length(sens_mm)) {
  sm <- do.call(rbind, sens_mm)
  write.csv(sm, file.path(MM_DIR,sprintf("mm_sensitivity_%s.csv",index_label)),row.names=FALSE)
}

####################################################################################
# EXCEL WORKBOOK
####################################################################################
cat("\n-- Building Excel workbook --\n")

wb  <- openxlsx::createWorkbook()
hdr <- openxlsx::createStyle(textDecoration="bold",fgFill="#1B3A6B",fontColour="#FFFFFF")
add_sheet <- function(wb, name, df) {
  openxlsx::addWorksheet(wb, name); openxlsx::writeData(wb, name, df)
  openxlsx::addStyle(wb, name, hdr, rows=1, cols=seq_len(ncol(df)))
}

if (length(all_mm_results)) {
  mm_all_df <- do.call(rbind, lapply(names(all_mm_results), function(nm) {
    r<-all_mm_results[[nm]]; if (is.null(r)) return(NULL); d<-r$params_tbl; d$Lag<-r$lag; d }))
  if (!is.null(mm_all_df)) add_sheet(wb,"MM_All_Lags",mm_all_df)
}
if (length(mm_lag_profile))  add_sheet(wb,"MM_Lag_Profile", do.call(rbind,mm_lag_profile))
if (length(sens_mm))         add_sheet(wb,"MM_Sensitivity",  do.call(rbind,sens_mm))
if (length(gam_lag_prof))    add_sheet(wb,"GAM_Lag_Profile", do.call(rbind,gam_lag_prof))
if (length(sens_gam))        add_sheet(wb,"GAM_Sensitivity", do.call(rbind,sens_gam))
if (length(convergence_rows))
  add_sheet(wb,"Convergence_GAM_vs_MM", do.call(rbind,convergence_rows))

boot_meta <- data.frame(
  Parameter=c("Bootstrap_type","Block_length_months","Full_R","Screened_R",
              "Screen_threshold_ti_ONI_PDO_p","CI_type","P_value_method"),
  Value=c("Stationary block bootstrap (boot::tsboot, sim='geom')",
          as.character(BLOCK_LENGTH), as.character(N_BOOTSTRAP),
          as.character(N_BOOT_SCREEN), as.character(GAM_SCREEN_P),
          "Percentile (2.5%, 97.5% of bootstrap distribution)",
          "Proportion method: 2*min(P(theta*<0), P(theta*>0))"),
  stringsAsFactors=FALSE)
add_sheet(wb,"Bootstrap_Metadata",boot_meta)

xlsx_out <- file.path(IDX_DIR,sprintf("results_summary_%s.xlsx",index_label))
openxlsx::saveWorkbook(wb, xlsx_out, overwrite=TRUE)
cat(sprintf("  Saved: %s\n", basename(xlsx_out)))

idx_elapsed__ <- elapsed_sec(idx_t0)
cat(sprintf("\n  OK %s complete in %s | script elapsed %s\n  Outputs: %s\n",
            index_label, fmt_dur(idx_elapsed__),
            fmt_dur(elapsed_sec(SCRIPT_START)), IDX_DIR))

.idx_result <- list(index=index_label, status="complete",
                    optimal_lag=optimal_lag, elapsed_sec=idx_elapsed__)