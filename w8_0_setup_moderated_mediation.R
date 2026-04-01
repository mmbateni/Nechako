####################################################################################
# w8_00_setup.R
# SHARED SETUP: configuration, timing helpers, all analysis functions, TELE_MASTER
#
# SOURCE THIS FILE FIRST.  Every other w8_*.R script begins with:
#   source("w8_00_setup.R")
#
# Depends on: DROUGHT_ANALYSIS_utils.R, utils_teleconnection_addon.R
####################################################################################

source("DROUGHT_ANALYSIS_utils.R")
source("utils_teleconnection_addon.R")

utils_load_packages(c(
  "lavaan","mgcv","ggplot2","patchwork","dplyr","openxlsx","tidyr","scales",
  "boot","parallel","future","future.apply","lubridate"
))

setwd(WD_PATH)
MM_DIR <- file.path(WD_PATH, "moderated_mediation")
dir.create(MM_DIR, showWarnings=FALSE, recursive=TRUE)

####################################################################################
# TIMING HELPERS
####################################################################################

fmt_dur <- function(secs) {
  secs <- max(0, round(secs))
  h <- secs %/% 3600L; m <- (secs %% 3600L) %/% 60L; s <- secs %% 60L
  if (h > 0) sprintf("%dh %02dm %02ds", h, m, s)
  else if (m > 0) sprintf("%dm %02ds", m, s)
  else sprintf("%ds", s)
}

elapsed_sec <- function(t0) (proc.time() - t0)[["elapsed"]]

fmt_progress <- function(done, total, elapsed_s, eta_s, width=40L) {
  pct    <- if (total > 0) done / total else 0
  filled <- round(pct * width)
  bar    <- paste0(strrep("=", max(0, filled - 1L)),
                   if (filled > 0 && done < total) ">" else if (filled > 0) "=" else "",
                   strrep(" ", max(0, width - filled)))
  sprintf("  [%s] %3.0f%% | %d/%d | elapsed %s | ETA %s",
          bar, pct * 100, done, total,
          fmt_dur(elapsed_s),
          if (done >= total) "done" else fmt_dur(eta_s))
}

####################################################################################
# CONFIGURATION
####################################################################################

ALL_INDICES <- list(
  list(index="spi",  scale=1),  list(index="spi",  scale=3),
  list(index="spi",  scale=6),  list(index="spi",  scale=9),
  list(index="spi",  scale=12), list(index="spei", scale=1),
  list(index="spei", scale=3),  list(index="spei", scale=6),
  list(index="spei", scale=9),  list(index="spei", scale=12)
)

LAGS_TO_TEST  <- 0:6
N_BOOTSTRAP   <- 5000L;  BOOTSTRAP_SEED <- 42L
W_VALUES      <- c(-1, 0, 1)
W_LABELS      <- c("Cool PDO (-1 SD)", "Neutral PDO (Mean)", "Warm PDO (+1 SD)")
BLOCK_LENGTH  <- 12L     # stationary bootstrap block length (months)
GAM_SCREEN_P  <- 0.20    # PERF-4: if ti(ONI,PDO) p > this, use reduced bootstrap
N_BOOT_SCREEN <- 500L    # reduced bootstrap reps when screened
JN_BOOT_N     <- 1200L   # JN floodlight reps per w-value
N_WORKERS     <- min(length(ALL_INDICES), max(1L, parallel::detectCores() - 1L))
SENSITIVITY_INDICES <- list(list(index="swei", scale=SWEI_SCALE))
GAM_K_MAIN   <- 8L;  GAM_K_TENSOR <- 6L;  GAM_RESPONSE <- "drought_value"

####################################################################################
# PERF-1: LOAD TELECONNECTIONS ONCE
####################################################################################

.load_tele_master <- function() {
  cat("-- Pre-loading teleconnection masters (PERF-1) --\n")
  oni <- load_teleconnection("oni", start_year=1950L, end_year=2025L)
  pdo <- load_teleconnection("pdo", start_year=1950L, end_year=2025L)
  pna <- load_teleconnection("pna", start_year=1950L, end_year=2025L)
  colnames(oni)[2] <- "ONI"; colnames(pdo)[2] <- "PDO"; colnames(pna)[2] <- "PNA"
  tm  <- Reduce(function(a,b) merge(a, b, by="date", all=FALSE), list(oni, pdo, pna))
  cat(sprintf("  OK TELE_MASTER: %d months (%s - %s)\n\n",
              nrow(tm), min(tm$date), max(tm$date)))
  tm
}
if (!exists("TELE_MASTER")) TELE_MASTER <- .load_tele_master()

####################################################################################
# LAVAAN MODEL BUILDER
####################################################################################

build_lavaan_model <- function(w_low=-1, w_mid=0, w_high=1) {
  sprintf('
  PNA_c ~ a1*ONI_c + a2*PDO_c + a3*ONI_PDO_int
  Y_c   ~ b1*PNA_c + cprime*ONI_c + b3*PDO_c + b4*PNA_PDO_int
  IE.low   := (a1 + a3*%s) * (b1 + b4*%s)
  IE.mid   := (a1 + a3*%s) * (b1 + b4*%s)
  IE.high  := (a1 + a3*%s) * (b1 + b4*%s)
  IE.total := a1 * b1
  TE       := cprime + a1*b1
  PM       := (a1*b1) / TE
  IMM      := a3*b1 + a1*b4 + a3*b4
  ', w_low,w_low, w_mid,w_mid, w_high,w_high)
}

.get_est <- function(row) {
  if ("est.std" %in% names(row) && length(row$est.std) > 0 && !is.na(row$est.std[1]))
    row$est.std[1]
  else if ("est" %in% names(row) && length(row$est) > 0) row$est[1]
  else NA_real_
}

####################################################################################
# FIT MODERATED MEDIATION — STATIONARY BLOCK BOOTSTRAP (PERF-5)
####################################################################################

fit_moderated_mediation <- function(df, lag=0L, n_boot=N_BOOTSTRAP,
                                    seed=BOOTSTRAP_SEED, block_length=BLOCK_LENGTH) {
  cat(sprintf("\n  [MM] Fitting lag=%d | n=%d | block_boot R=%d l=%d\n",
              lag, nrow(df), n_boot, block_length))
  if (nrow(df) < 60) warning(sprintf("Only %d obs at lag=%d", nrow(df), lag))
  model_str <- build_lavaan_model(W_VALUES[1], W_VALUES[2], W_VALUES[3])
  
  set.seed(seed)
  fit0 <- tryCatch(
    lavaan::sem(model=model_str, data=df, se="none",
                estimator="ML", fixed.x=FALSE, missing="listwise"),
    error=function(e) { cat("  XX lavaan:", e$message, "\n"); NULL })
  if (is.null(fit0)) return(NULL)
  if (!lavaan::lavTech(fit0, "converged")) warning(sprintf("MM not converged lag=%d", lag))
  pe0 <- lavaan::parameterEstimates(fit0, ci=FALSE, standardized=TRUE)
  r2  <- lavaan::lavInspect(fit0, "r2")
  cat(sprintf("  R2(PNA)=%.3f | R2(Y)=%.3f\n", r2["PNA_c"], r2["Y_c"]))
  
  PARAM_LBLS <- c("a1","a2","a3","b1","cprime","b3","b4",
                  "IE.low","IE.mid","IE.high","IE.total","TE","PM","IMM")
  .boot_stat <- function(d) {
    fb <- tryCatch(lavaan::sem(model=model_str, data=d, se="none",
                               estimator="ML", fixed.x=FALSE, missing="listwise"),
                   error=function(e) NULL)
    if (is.null(fb) || !isTRUE(lavaan::lavTech(fb,"converged")))
      return(rep(NA_real_, length(PARAM_LBLS)))
    pb <- lavaan::parameterEstimates(fb, ci=FALSE, standardized=TRUE)
    pb_def <- pb[pb$op %in% c("~",":="), ]
    sapply(PARAM_LBLS, function(lbl) {
      r <- pb_def[pb_def$label==lbl, ]
      if (!nrow(r)) return(NA_real_)
      if ("est.std" %in% names(r) && !is.na(r$est.std[1])) r$est.std[1]
      else if ("est" %in% names(r)) r$est[1] else NA_real_
    })
  }
  
  set.seed(seed)
  boot_out <- tryCatch(
    boot::tsboot(tseries=df, statistic=.boot_stat, R=n_boot,
                 l=block_length, sim="geom"),
    error=function(e) { cat("  XX block bootstrap:", e$message, "\n"); NULL })
  
  pe_def <- pe0[pe0$op %in% c("~",":="), ]
  pe_def$ci.lower <- NA_real_; pe_def$ci.upper <- NA_real_; pe_def$pvalue <- NA_real_
  
  if (!is.null(boot_out)) {
    boot_mat <- boot_out$t; colnames(boot_mat) <- PARAM_LBLS
    n_ok <- sum(apply(boot_mat, 1, function(r) all(is.finite(r))))
    cat(sprintf("  Block bootstrap: %d/%d replicates converged\n", n_ok, n_boot))
    for (lbl in PARAM_LBLS) {
      ri <- which(pe_def$label==lbl); if (!length(ri)) next
      b  <- boot_mat[,lbl]; b <- b[is.finite(b)]; if (length(b)<100) next
      ci <- quantile(b, c(0.025,0.975))
      pe_def$ci.lower[ri] <- unname(ci[1]); pe_def$ci.upper[ri] <- unname(ci[2])
      pt <- if ("est.std" %in% names(pe_def) && !is.na(pe_def$est.std[ri]))
        pe_def$est.std[ri] else pe_def$est[ri]
      if (!is.na(pt))
        pe_def$pvalue[ri] <- min(1.0, 2*min(mean(b<0,na.rm=TRUE),mean(b>0,na.rm=TRUE)))
    }
  } else cat("  ! Block bootstrap failed -- CIs/p-values will be NA\n")
  
  pe_def$label_clean <- tryCatch(dplyr::case_when(
    pe_def$label=="a1"       ~ "a1: ONI->PNA (mean PDO)",
    pe_def$label=="a2"       ~ "a2: PDO->PNA direct",
    pe_def$label=="a3"       ~ "a3: ONIxPDO->PNA (mod a-path)",
    pe_def$label=="b1"       ~ "b1: PNA->Y (mean PDO)",
    pe_def$label=="cprime"   ~ "c': ONI->Y direct",
    pe_def$label=="b3"       ~ "b3: PDO->Y direct",
    pe_def$label=="b4"       ~ "b4: PNAxPDO->Y (mod b-path)",
    pe_def$label=="IE.low"   ~ sprintf("IE (PDO=%.0f SD)",  W_VALUES[1]),
    pe_def$label=="IE.mid"   ~ sprintf("IE (PDO=%.0f SD)",  W_VALUES[2]),
    pe_def$label=="IE.high"  ~ sprintf("IE (PDO=+%.0f SD)", W_VALUES[3]),
    pe_def$label=="IE.total" ~ "Total IE (mean PDO)",
    pe_def$label=="TE"       ~ "Total Effect ONI->Y",
    pe_def$label=="PM"       ~ "Proportion mediated",
    pe_def$label=="IMM"      ~ "IMM (Index of Mod. Mediation)",
    TRUE ~ pe_def$label),
    error=function(e) { cat("  warning: label_clean fallback\n"); pe_def$label })
  if (!"label_clean" %in% names(pe_def)) pe_def$label_clean <- pe_def$label
  for (nc in c("ci.lower","ci.upper","pvalue"))
    if (!nc %in% names(pe_def)) pe_def[[nc]] <- NA_real_
  
  est_col  <- if ("est.std" %in% names(pe_def)) "est.std" else "est"
  pe_print <- pe_def[, c("label_clean",est_col,"ci.lower","ci.upper","pvalue")]
  colnames(pe_print) <- c("Parameter","Std.Est","CI.lower","CI.upper","p-value")
  pe_print[,2:4] <- round(pe_print[,2:4],4); pe_print[,5] <- round(pe_print[,5],4)
  print(pe_print, row.names=FALSE)
  
  imm_row <- pe_def[pe_def$label=="IMM",]
  if (nrow(imm_row))
    cat(sprintf("  IMM=%.4f [%.4f,%.4f] %s\n",
                .get_est(imm_row), imm_row$ci.lower, imm_row$ci.upper,
                if (!is.na(imm_row$pvalue) && imm_row$pvalue<0.05)
                  "OK PDO moderates indirect path" else "-- PDO moderation not significant"))
  
  list(fit=fit0, params=pe_def, params_tbl=pe_print, r2=r2, n=nrow(df), lag=lag,
       boot_out=boot_out, block_length=block_length)
}

####################################################################################
# JOHNSON-NEYMAN FLOODLIGHT  (sequential — no nested PSOCK)
####################################################################################

compute_jn_floodlight <- function(df, n_boot=JN_BOOT_N, w_seq=seq(-2.5,2.5,by=0.1)) {
  cat("  [JN] Floodlight analysis...\n")
  .jn_one <- function(w, df_inner, n_boot_inner, boot_seed) {
    tryCatch({
      la   <- lm(PNA_c~ONI_c+PDO_c+ONI_PDO_int,         data=df_inner)
      lb   <- lm(Y_c~PNA_c+ONI_c+PDO_c+PNA_PDO_int,     data=df_inner)
      ca   <- coef(la); cb <- coef(lb)
      ie_w <- (ca["ONI_c"]+ca["ONI_PDO_int"]*w)*(cb["PNA_c"]+cb["PNA_PDO_int"]*w)
      set.seed(boot_seed)
      boot_ie <- replicate(n_boot_inner, {
        idx <- sample(nrow(df_inner),replace=TRUE)
        la2 <- tryCatch(lm(PNA_c~ONI_c+PDO_c+ONI_PDO_int,     data=df_inner[idx,]),error=function(e)NULL)
        lb2 <- tryCatch(lm(Y_c~PNA_c+ONI_c+PDO_c+PNA_PDO_int, data=df_inner[idx,]),error=function(e)NULL)
        if (is.null(la2)||is.null(lb2)) return(NA_real_)
        ca2<-coef(la2); cb2<-coef(lb2)
        (ca2["ONI_c"]+ca2["ONI_PDO_int"]*w)*(cb2["PNA_c"]+cb2["PNA_PDO_int"]*w)
      })
      boot_ie <- boot_ie[is.finite(boot_ie)]
      ci_lo <- if (length(boot_ie)<100) NA_real_ else unname(quantile(boot_ie,0.025))
      ci_hi <- if (length(boot_ie)<100) NA_real_ else unname(quantile(boot_ie,0.975))
      data.frame(w=w,IE=unname(ie_w),ci_lower=ci_lo,ci_upper=ci_hi,
                 significant=!is.na(ci_lo)&&(ci_lo>0||ci_hi<0),stringsAsFactors=FALSE)
    }, error=function(e) NULL)
  }
  cat("  [JN] Running sequentially (safe inside future worker)\n")
  res_list <- lapply(w_seq, function(w) .jn_one(w, df, n_boot, BOOTSTRAP_SEED))
  res_list <- res_list[!vapply(res_list, is.null, logical(1))]
  if (length(res_list)==0) return(NULL)
  do.call(rbind, res_list)
}

####################################################################################
# MM VISUALISATION FUNCTIONS
####################################################################################

plot_conditional_ie <- function(mm_result, lag, index_label) {
  pe <- mm_result$params
  ie_rows <- pe[pe$label %in% c("IE.low","IE.mid","IE.high"),]
  if (!nrow(ie_rows)) return(NULL)
  ie_df <- data.frame(PDO_level=factor(W_LABELS,levels=W_LABELS),
                      IE=sapply(seq_len(nrow(ie_rows)),function(i).get_est(ie_rows[i,])),
                      ci_lo=ie_rows$ci.lower, ci_hi=ie_rows$ci.upper,
                      sig=!is.na(ie_rows$pvalue)&ie_rows$pvalue<0.05)
  ggplot2::ggplot(ie_df,ggplot2::aes(x=PDO_level,y=IE,colour=sig,shape=sig)) +
    ggplot2::geom_hline(yintercept=0,linetype="dashed",colour="grey50") +
    ggplot2::geom_errorbar(ggplot2::aes(ymin=ci_lo,ymax=ci_hi),width=0.15,linewidth=1.1) +
    ggplot2::geom_point(size=4) +
    ggplot2::scale_colour_manual(values=c("FALSE"="grey60","TRUE"="#C0392B"),name="p<0.05") +
    ggplot2::scale_shape_manual(values=c("FALSE"=1,"TRUE"=19),name="p<0.05") +
    shared_ts_theme(12) +
    ggplot2::labs(title=sprintf("Conditional Indirect Effect: ONI->PNA->%s",index_label),
                  subtitle=sprintf("Lag=%d months | 95%% Block-bootstrap CI (l=%d) | n=%d",
                                   lag,mm_result$block_length,mm_result$n),
                  x="PDO Phase",y="IE (standardised)")
}

plot_jn_floodlight <- function(jn_df, lag, index_label) {
  if (is.null(jn_df)||!nrow(jn_df)) return(NULL)
  ggplot2::ggplot(jn_df,ggplot2::aes(x=w)) +
    ggplot2::geom_ribbon(data=jn_df[jn_df$significant,],
                         ggplot2::aes(ymin=ci_lower,ymax=ci_upper),fill="#FADBD8",alpha=0.6) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin=ci_lower,ymax=ci_upper),fill="grey85",alpha=0.4) +
    ggplot2::geom_line(ggplot2::aes(y=IE),colour="#C0392B",linewidth=1.1) +
    ggplot2::geom_hline(yintercept=0,linetype="dashed") +
    ggplot2::geom_vline(xintercept=W_VALUES,linetype="dotted",colour="#2980B9",linewidth=0.6) +
    shared_ts_theme(12) +
    ggplot2::labs(title="Johnson-Neyman: Region of significance",
                  subtitle=sprintf("IE(ONI->PNA->%s) as f(PDO) | lag=%d | Pink=significant",
                                   index_label,lag),
                  x="PDO (standardised)",y="IE (standardised)")
}

plot_path_diagram <- function(mm_result, lag, index_label) {
  pe <- mm_result$params
  ge <- function(lbl) {
    row <- pe[pe$label==lbl,]; if (!nrow(row)) return("n/a")
    p_str <- if (!is.na(row$pvalue[1])) sprintf("p=%.3f",row$pvalue[1]) else "p=NA"
    sprintf("b=%.3f\n%s",.get_est(row),p_str)
  }
  ggplot2::ggplot()+ggplot2::xlim(0,10)+ggplot2::ylim(0,6)+ggplot2::theme_void()+
    ggplot2::annotate("rect",xmin=0.3,xmax=2.3,ymin=2.5,ymax=3.5,fill="#D6EAF8",colour="#2980B9",linewidth=1.2)+
    ggplot2::annotate("text",x=1.3,y=3.0,label="ONI\n(X)",size=3.5,fontface="bold")+
    ggplot2::annotate("rect",xmin=3.8,xmax=6.2,ymin=4.2,ymax=5.2,fill="#D5F5E3",colour="#27AE60",linewidth=1.2)+
    ggplot2::annotate("text",x=5.0,y=4.7,label="PNA\n(Mediator M)",size=3.5,fontface="bold")+
    ggplot2::annotate("rect",xmin=7.7,xmax=9.7,ymin=2.5,ymax=3.5,fill="#FADBD8",colour="#C0392B",linewidth=1.2)+
    ggplot2::annotate("text",x=8.7,y=3.0,label=sprintf("%s\n(Y)",index_label),size=3.5,fontface="bold")+
    ggplot2::annotate("rect",xmin=3.8,xmax=6.2,ymin=0.8,ymax=1.8,fill="#FCF3CF",colour="#F39C12",linewidth=1.2)+
    ggplot2::annotate("text",x=5.0,y=1.3,label="PDO\n(Moderator W)",size=3.5,fontface="bold")+
    ggplot2::annotate("segment",x=2.3,xend=3.8,y=3.3,yend=4.5,arrow=ggplot2::arrow(length=ggplot2::unit(0.3,"cm"),type="closed"),colour="#27AE60",linewidth=1)+
    ggplot2::annotate("text",x=2.9,y=4.1,label=paste("a1:",ge("a1")),size=3,colour="#27AE60",hjust=0)+
    ggplot2::annotate("segment",x=6.2,xend=7.7,y=4.5,yend=3.3,arrow=ggplot2::arrow(length=ggplot2::unit(0.3,"cm"),type="closed"),colour="#C0392B",linewidth=1)+
    ggplot2::annotate("text",x=6.7,y=4.1,label=paste("b1:",ge("b1")),size=3,colour="#C0392B",hjust=0)+
    ggplot2::annotate("segment",x=2.3,xend=7.7,y=2.8,yend=2.8,arrow=ggplot2::arrow(length=ggplot2::unit(0.3,"cm"),type="closed"),colour="#7F8C8D",linewidth=0.8,linetype="dashed")+
    ggplot2::annotate("text",x=5.0,y=2.55,label=paste("c':",ge("cprime")),size=3,colour="#7F8C8D")+
    ggplot2::annotate("segment",x=4.4,xend=3.1,y=1.8,yend=3.8,arrow=ggplot2::arrow(length=ggplot2::unit(0.25,"cm"),type="open"),colour="#F39C12",linewidth=0.8)+
    ggplot2::annotate("text",x=3.0,y=2.8,label=paste("a3:",ge("a3")),size=2.8,colour="#F39C12",hjust=1)+
    ggplot2::annotate("segment",x=5.6,xend=6.9,y=1.8,yend=3.8,arrow=ggplot2::arrow(length=ggplot2::unit(0.25,"cm"),type="open"),colour="#F39C12",linewidth=0.8)+
    ggplot2::annotate("text",x=7.0,y=2.8,label=paste("b4:",ge("b4")),size=2.8,colour="#F39C12",hjust=0)+
    ggplot2::labs(title=sprintf("Path Diagram: %s | Lag=%d months",index_label,lag),
                  subtitle=sprintf("IMM: %s",ge("IMM")),
                  caption=sprintf("Dashed=direct c' path | Orange=PDO moderation | Block bootstrap l=%d",mm_result$block_length))+
    ggplot2::theme(plot.title=ggplot2::element_text(face="bold",size=13,hjust=0.5),
                   plot.subtitle=ggplot2::element_text(size=10,hjust=0.5,colour="#C0392B"))
}

####################################################################################
# GLOBALS EXPORT LIST FOR future WORKERS
####################################################################################

.w8_globals <- c(
  "TELE_MASTER",
  "build_lavaan_model","fit_moderated_mediation",".get_est",
  "compute_jn_floodlight","plot_conditional_ie","plot_jn_floodlight","plot_path_diagram",
  "fmt_dur","elapsed_sec","fmt_progress","SCRIPT_START",
  "ALL_INDICES","LAGS_TO_TEST","N_BOOTSTRAP","BOOTSTRAP_SEED",
  "BLOCK_LENGTH","GAM_SCREEN_P","N_BOOT_SCREEN","JN_BOOT_N",
  "W_VALUES","W_LABELS","MM_DIR","GAM_K_MAIN","GAM_K_TENSOR","GAM_RESPONSE",
  "SENSITIVITY_INDICES","N_WORKERS","SWEI_SCALE","WD_PATH"
)

cat("-- w8_00_setup.R loaded: config, functions, TELE_MASTER ready --\n\n")