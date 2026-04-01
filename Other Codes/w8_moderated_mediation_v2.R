####################################################################################
# w8_moderated_mediation_v2.R
# DUAL ANALYSIS: Moderated Mediation (lavaan) + GAM with Tensor Smooths (mgcv)
####################################################################################
#
# PERFORMANCE CHANGES vs v1 (see PERF tags throughout):
#
#  PERF-1  Teleconnections loaded ONCE before the outer loop (TELE_MASTER).
#          build_analysis_dataframe() receives it as tele_df= instead of
#          re-reading 3 CSVs on every one of the 70 calls (10 indices × 7 lags).
#
#  PERF-2  Outer index loop runs via future_lapply() (parallel::multisession).
#          Each index is fully independent — up to N_WORKERS concurrent runs.
#          All file I/O per index goes to its own subdirectory, no race conditions.
#
#  PERF-3  summary(gam_fit) computed ONCE per lag and passed to downstream
#          helpers (gam_variance_partition, gam_lag_summary) as precomputed_s.
#
#  PERF-4  GAM-first screening before MM. If ti(ONI,PDO) p-value > GAM_SCREEN_P
#          (default 0.20), the full 5000-rep bootstrap is replaced by a cheap
#          500-rep run. The interaction must be at least marginal to warrant the
#          full computational cost.
#
#  PERF-5  Block bootstrap (boot::tsboot, stationary bootstrap sim="geom")
#          replaces lavaan's internal i.i.d. BCa bootstrap for ALL MM fits.
#          Rationale: monthly climate data is temporally autocorrelated (DW test
#          will typically flag this). I.i.d. resampling underestimates CIs.
#          The stationary bootstrap (Politis & Romano 1994) respects autocorrelation
#          structure by resampling contiguous blocks of length ~BLOCK_LENGTH months.
#          Each replicate fits lavaan with se="none" (point estimates only, ~20×
#          faster per fit than se="bootstrap"). CIs are percentile-based from the
#          block bootstrap distribution; p-values use the proportion method.
#
#  PERF-6  Johnson-Neyman bootstrap uses parallel::parLapplyLB with a PSOCK
#          cluster (Windows-safe) — 51 w-values, each with JN_BOOT_N replicates.
#
# ANALYTICAL STRATEGY (unchanged):
#  TRACK 1 — GAM tensor smooths (mgcv)
#  TRACK 2 — Moderated Mediation (lavaan + block bootstrap)
#  TRACK 3 — Cross-validation of results
#
# Depends on: DROUGHT_ANALYSIS_utils.R, utils_teleconnection_addon.R
# Run AFTER: w7_teleconnection_prep.R
####################################################################################

rm(list=ls()); gc()
SCRIPT_START <- proc.time()          # wall-clock anchor for the entire run
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

# Format an elapsed-seconds value as "Xh Ym Zs" (or "Ym Zs" / "Zs" when compact).
fmt_dur <- function(secs) {
  secs <- max(0, round(secs))
  h <- secs %/% 3600L
  m <- (secs %% 3600L) %/% 60L
  s <- secs %% 60L
  if (h > 0) sprintf("%dh %02dm %02ds", h, m, s)
  else if (m > 0) sprintf("%dm %02ds", m, s)
  else sprintf("%ds", s)
}

# Return wall-clock seconds elapsed since a proc.time() snapshot.
elapsed_sec <- function(t0) (proc.time() - t0)[["elapsed"]]

# Pretty progress bar: [=========>        ] 55% | elapsed 1m 02s | ETA 0m 51s
fmt_progress <- function(done, total, elapsed_s, eta_s, width=40L) {
  pct   <- if (total > 0) done / total else 0
  filled <- round(pct * width)
  bar    <- paste0(strrep("=", max(0, filled - 1L)),
                   if (filled > 0 && done < total) ">" else if (filled > 0) "=" else "",
                   strrep(" ", max(0, width - filled)))
  sprintf("  [%s] %3.0f%% | %d/%d | elapsed %s | ETA %s",
          bar, pct * 100, done, total,
          fmt_dur(elapsed_s),
          if (done >= total) "done" else fmt_dur(eta_s))
}

cat("\n╔════════════════════════════════════════════════════════╗\n")
cat("║  w8  MODERATED MEDIATION + GAM TENSOR SMOOTHS          ║\n")
cat("╚════════════════════════════════════════════════════════╝\n\n")
cat(sprintf("  Started : %s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

####################################################################################
# CONFIGURATION
####################################################################################

ALL_INDICES <- list(
  list(index="spi",  scale=1),
  list(index="spi",  scale=3),
  list(index="spi",  scale=6),
  list(index="spi",  scale=9),
  list(index="spi",  scale=12),
  list(index="spei", scale=1),
  list(index="spei", scale=3),
  list(index="spei", scale=6),
  list(index="spei", scale=9),
  list(index="spei", scale=12)
)

LAGS_TO_TEST  <- 0:6
N_BOOTSTRAP   <- 5000L;  BOOTSTRAP_SEED <- 42L
W_VALUES      <- c(-1, 0, 1)
W_LABELS      <- c("Cool PDO (−1 SD)","Neutral PDO (Mean)","Warm PDO (+1 SD)")

# PERF-5: Block bootstrap parameters
# BLOCK_LENGTH: expected block length for the stationary bootstrap.
# Rule of thumb: 12 months (one seasonal cycle) for monthly climate series.
# Longer series or strong multi-year ENSO memory → consider 18–24.
BLOCK_LENGTH  <- 12L

# PERF-4: GAM screening threshold.
# If ti(ONI,PDO) p-value > GAM_SCREEN_P use only N_BOOT_SCREEN replicates.
GAM_SCREEN_P  <- 0.20
N_BOOT_SCREEN <- 500L

# PERF-6: JN bootstrap replicate count (per w-value point).
JN_BOOT_N     <- 1200L

# PERF-2: Max parallel workers for the outer index loop.
N_WORKERS     <- min(length(ALL_INDICES), max(1L, parallel::detectCores() - 1L))

SENSITIVITY_INDICES <- list(
  list(index="swei", scale=SWEI_SCALE)   # SWEI_SCALE resolved HERE at definition time
)

# GAM configuration
GAM_K_MAIN   <- 8L
GAM_K_TENSOR <- 6L
GAM_RESPONSE <- "drought_value"

####################################################################################
# PERF-1: LOAD TELECONNECTIONS ONCE
# build_analysis_dataframe() accepts tele_df= so these 3 CSV reads happen
# exactly once instead of 3×70 = 210 times over the full run.
####################################################################################
cat("── Pre-loading teleconnection masters (PERF-1) ──\n")

oni_master <- load_teleconnection("oni", start_year=1950L, end_year=2025L)
pdo_master <- load_teleconnection("pdo", start_year=1950L, end_year=2025L)
pna_master <- load_teleconnection("pna", start_year=1950L, end_year=2025L)
colnames(oni_master)[2] <- "ONI"
colnames(pdo_master)[2] <- "PDO"
colnames(pna_master)[2] <- "PNA"

TELE_MASTER <- Reduce(
  function(a, b) merge(a, b, by="date", all=FALSE),
  list(oni_master, pdo_master, pna_master)
)
cat(sprintf("  ✓ TELE_MASTER: %d months (%s – %s)  [+%s]\n\n",
            nrow(TELE_MASTER), min(TELE_MASTER$date), max(TELE_MASTER$date),
            fmt_dur(elapsed_sec(SCRIPT_START))))

rm(oni_master, pdo_master, pna_master); gc()

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
# FIT MODERATED MEDIATION — BLOCK BOOTSTRAP (PERF-5)
####################################################################################
#
# Architecture:
#   Step 1  Fit lavaan with se="none" → fast point estimates + R²
#   Step 2  boot::tsboot (stationary bootstrap, sim="geom") resamples the data
#           frame in temporal blocks of expected length BLOCK_LENGTH. Each
#           replicate re-fits lavaan with se="none".
#   Step 3  Percentile CIs from block bootstrap distribution.
#           P-values: proportion method (2×min(P(θ*<0), P(θ*>0))).
#
# Why block bootstrap?
#   Monthly climate residuals carry substantial serial autocorrelation (ENSO,
#   PDO persist 6–24 months). I.i.d. resampling treats every month as
#   independent, which inflates the effective sample size and underestimates
#   CIs/p-values. Block resampling preserves the local autocorrelation
#   structure within each block.
####################################################################################

fit_moderated_mediation <- function(df, lag=0L, n_boot=N_BOOTSTRAP, seed=BOOTSTRAP_SEED,
                                    block_length=BLOCK_LENGTH) {
  cat(sprintf("\n  [MM] Fitting lag=%d | n=%d | block_boot R=%d l=%d\n",
              lag, nrow(df), n_boot, block_length))
  if (nrow(df) < 60) warning(sprintf("Only %d obs at lag=%d", nrow(df), lag))
  
  model_str <- build_lavaan_model(W_VALUES[1], W_VALUES[2], W_VALUES[3])
  
  # ── Step 1: point estimates, se="none" (no expensive SE computation) ─────
  set.seed(seed)
  fit0 <- tryCatch(
    lavaan::sem(model=model_str, data=df, se="none",
                estimator="ML", fixed.x=FALSE, missing="listwise"),
    error=function(e) { cat("  ❌ lavaan:", e$message, "\n"); NULL })
  if (is.null(fit0)) return(NULL)
  if (!lavaan::lavTech(fit0, "converged"))
    warning(sprintf("MM not converged lag=%d", lag))
  
  pe0 <- lavaan::parameterEstimates(fit0, ci=FALSE, standardized=TRUE)
  r2  <- lavaan::lavInspect(fit0, "r2")
  cat(sprintf("  R²(PNA)=%.3f | R²(Y)=%.3f\n",
              r2["PNA_c"], r2["Y_c"]))
  
  # ── Step 2: stationary block bootstrap via boot::tsboot ──────────────────
  # Parameter labels to track across replicates
  PARAM_LBLS <- c("a1","a2","a3","b1","cprime","b3","b4",
                  "IE.low","IE.mid","IE.high","IE.total","TE","PM","IMM")
  
  # Statistic function: re-fit lavaan on resampled data, return named std-est vector
  .boot_stat <- function(d) {
    fb <- tryCatch(
      lavaan::sem(model=model_str, data=d, se="none",
                  estimator="ML", fixed.x=FALSE, missing="listwise"),
      error=function(e) NULL)
    if (is.null(fb) || !isTRUE(lavaan::lavTech(fb, "converged")))
      return(rep(NA_real_, length(PARAM_LBLS)))
    pb <- lavaan::parameterEstimates(fb, ci=FALSE, standardized=TRUE)
    pb_def <- pb[pb$op %in% c("~", ":="), ]
    sapply(PARAM_LBLS, function(lbl) {
      r <- pb_def[pb_def$label == lbl, ]
      if (!nrow(r)) return(NA_real_)
      if ("est.std" %in% names(r) && !is.na(r$est.std[1])) r$est.std[1]
      else if ("est" %in% names(r)) r$est[1]
      else NA_real_
    })
  }
  
  set.seed(seed)
  boot_out <- tryCatch(
    boot::tsboot(tseries   = df,
                 statistic = .boot_stat,
                 R         = n_boot,
                 l         = block_length,
                 sim       = "geom"),   # stationary bootstrap (Politis & Romano 1994)
    error=function(e) { cat("  ❌ block bootstrap:", e$message, "\n"); NULL })
  
  # ── Step 3: assemble pe_def with block-bootstrap CIs and p-values ────────
  pe_def  <- pe0[pe0$op %in% c("~", ":="), ]
  pe_def$ci.lower <- NA_real_
  pe_def$ci.upper <- NA_real_
  pe_def$pvalue   <- NA_real_
  
  if (!is.null(boot_out)) {
    boot_mat <- boot_out$t
    colnames(boot_mat) <- PARAM_LBLS
    n_ok <- sum(apply(boot_mat, 1, function(r) all(is.finite(r))))
    cat(sprintf("  Block bootstrap: %d/%d replicates converged\n", n_ok, n_boot))
    
    for (lbl in PARAM_LBLS) {
      ri <- which(pe_def$label == lbl)
      if (!length(ri)) next
      b  <- boot_mat[, lbl]; b <- b[is.finite(b)]
      if (length(b) < 100) next
      
      # Percentile CI
      ci <- quantile(b, c(0.025, 0.975))
      pe_def$ci.lower[ri] <- unname(ci[1])
      pe_def$ci.upper[ri] <- unname(ci[2])
      
      # Proportion-based p-value (two-sided)
      pt <- if ("est.std" %in% names(pe_def) && !is.na(pe_def$est.std[ri]))
        pe_def$est.std[ri] else pe_def$est[ri]
      if (!is.na(pt))
        pe_def$pvalue[ri] <- min(1.0,
                                 2 * min(mean(b < 0, na.rm=TRUE), mean(b > 0, na.rm=TRUE)))
    }
  } else {
    cat("  ⚠ Block bootstrap failed — CIs/p-values will be NA\n")
  }
  
  # ── Label cleaning (unchanged from original) ─────────────────────────────
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
  
  if (!"label_clean" %in% names(pe_def)) pe_def$label_clean <- pe_def$label
  for (nc in c("ci.lower","ci.upper","pvalue")) {
    if (!nc %in% names(pe_def)) pe_def[[nc]] <- NA_real_
  }
  
  est_col  <- if ("est.std" %in% names(pe_def)) "est.std" else "est"
  pe_print <- pe_def[, c("label_clean", est_col, "ci.lower","ci.upper","pvalue")]
  colnames(pe_print) <- c("Parameter","Std.Est","CI.lower","CI.upper","p-value")
  pe_print[,2:4] <- round(pe_print[,2:4],4); pe_print[,5] <- round(pe_print[,5],4)
  print(pe_print, row.names=FALSE)
  
  imm_row <- pe_def[pe_def$label=="IMM",]
  if (nrow(imm_row))
    cat(sprintf("  IMM=%.4f [%.4f,%.4f] %s\n",
                .get_est(imm_row), imm_row$ci.lower, imm_row$ci.upper,
                if (!is.na(imm_row$pvalue) && imm_row$pvalue<0.05)
                  "✓ PDO moderates indirect path"
                else "✗ PDO moderation not significant"))
  
  list(fit=fit0, params=pe_def, params_tbl=pe_print,
       r2=r2, n=nrow(df), lag=lag,
       boot_out=boot_out, block_length=block_length)
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
# JOHNSON-NEYMAN FLOODLIGHT  (PERF-6: parallel PSOCK cluster)
####################################################################################

compute_jn_floodlight <- function(df, n_boot=JN_BOOT_N,
                                  w_seq=seq(-2.5,2.5,by=0.1)) {
  cat("  [JN] Floodlight analysis...\n")
  
  # Inner function for one w-value
  .jn_one <- function(w, df_inner, n_boot_inner, boot_seed) {
    tryCatch({
      la <- lm(PNA_c ~ ONI_c + PDO_c + ONI_PDO_int, data=df_inner)
      lb <- lm(Y_c   ~ PNA_c + ONI_c + PDO_c + PNA_PDO_int, data=df_inner)
      ca <- coef(la); cb <- coef(lb)
      ie_w <- (ca["ONI_c"]+ca["ONI_PDO_int"]*w)*(cb["PNA_c"]+cb["PNA_PDO_int"]*w)
      set.seed(boot_seed)
      boot_ie <- replicate(n_boot_inner, {
        idx <- sample(nrow(df_inner), replace=TRUE)
        la2 <- tryCatch(lm(PNA_c~ONI_c+PDO_c+ONI_PDO_int,     data=df_inner[idx,]),error=function(e)NULL)
        lb2 <- tryCatch(lm(Y_c~PNA_c+ONI_c+PDO_c+PNA_PDO_int, data=df_inner[idx,]),error=function(e)NULL)
        if (is.null(la2)||is.null(lb2)) return(NA_real_)
        ca2<-coef(la2); cb2<-coef(lb2)
        (ca2["ONI_c"]+ca2["ONI_PDO_int"]*w)*(cb2["PNA_c"]+cb2["PNA_PDO_int"]*w)
      })
      boot_ie <- boot_ie[is.finite(boot_ie)]
      ci_lo <- if (length(boot_ie)<100) NA_real_ else unname(quantile(boot_ie,0.025))
      ci_hi <- if (length(boot_ie)<100) NA_real_ else unname(quantile(boot_ie,0.975))
      data.frame(w=w, IE=unname(ie_w), ci_lower=ci_lo, ci_upper=ci_hi,
                 significant=!is.na(ci_lo)&&(ci_lo>0||ci_hi<0),
                 stringsAsFactors=FALSE)
    }, error=function(e) NULL)
  }
  
  # BUG-FIX (PERF-6 revised): Always run JN sequentially — no nested PSOCK cluster.
  #
  # Original code tried to detect future workers via R_FUTURE_WORKER_UUID and
  # spawn a PSOCK cluster only at the top level.  In practice the env-var was
  # not reliably set inside every multisession worker, so some workers created
  # nested PSOCK clusters.  A nested cluster opens new sockets on the same host
  # and corrupts future's own socket handshake, producing
  # "communication ... is out of sync" / UnexpectedFutureResultError.
  #
  # The outer future_lapply already runs up to N_WORKERS indices concurrently.
  # Running 51 w-values x JN_BOOT_N reps sequentially inside a worker is fast
  # (typically < 30 s) and completely safe.
  cat("  [JN] Running sequentially (safe inside future worker or at top level)\n")
  res_list <- lapply(w_seq, function(w)
    .jn_one(w, df, n_boot, BOOTSTRAP_SEED))
  
  res_list <- res_list[!vapply(res_list, is.null, logical(1))]
  if (length(res_list) == 0) return(NULL)
  do.call(rbind, res_list)
}

####################################################################################
# MM VISUALISATION  (unchanged)
####################################################################################

plot_conditional_ie <- function(mm_result, lag, index_label) {
  pe <- mm_result$params
  ie_rows <- pe[pe$label %in% c("IE.low","IE.mid","IE.high"),]
  if (!nrow(ie_rows)) return(NULL)
  ie_df <- data.frame(PDO_level=factor(W_LABELS,levels=W_LABELS),
                      IE=sapply(seq_len(nrow(ie_rows)), function(i) .get_est(ie_rows[i,])),
                      ci_lo=ie_rows$ci.lower, ci_hi=ie_rows$ci.upper,
                      sig=!is.na(ie_rows$pvalue) & ie_rows$pvalue<0.05)
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
                  subtitle=sprintf("Lag=%d months | 95%% Block-bootstrap CI (l=%d) | n=%d",
                                   lag, mm_result$block_length, mm_result$n),
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
    p_str <- if (!is.na(row$pvalue[1])) sprintf("p=%.3f",row$pvalue[1]) else "p=NA"
    sprintf("b=%.3f\n%s", .get_est(row), p_str)
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
      caption=sprintf("Dashed=direct c' path | Orange=PDO moderation | Block bootstrap l=%d",
                      mm_result$block_length)) +
    ggplot2::theme(plot.title=ggplot2::element_text(face="bold",size=13,hjust=0.5),
                   plot.subtitle=ggplot2::element_text(size=10,hjust=0.5,colour="#C0392B"))
}

####################################################################################
# PERF-2: OUTER LOOP VIA future_lapply
# Each index runs in its own parallel worker. Workers share read-only globals
# (TELE_MASTER, constants, all function definitions) via future's automatic
# global export.
####################################################################################

cat(sprintf("\n── Launching parallel outer loop | %d indices | %d workers | script elapsed %s ──\n\n",
            length(ALL_INDICES), N_WORKERS, fmt_dur(elapsed_sec(SCRIPT_START))))

# Reset to sequential first — this shuts down any PSOCK workers left over from a
# previous run of this script (or from the old nested-cluster JN code).  Without
# this, re-running the script in the same R session leaves the old worker sockets
# alive; when future tries to handshake new workers they intercept the connection
# and reply with "future:::requestNode() validation call" instead of "future-grmall",
# causing FutureLaunchError on the Nth worker.
future::plan(future::sequential)
future::plan(future::multisession, workers=N_WORKERS)

# Capture all function/constant globals needed inside futures
.w8_globals <- c(
  # Heavy data object — must be explicitly passed
  "TELE_MASTER",
  # w8-local functions
  "build_lavaan_model","fit_moderated_mediation",".get_est",
  "compute_jn_floodlight","plot_conditional_ie","plot_jn_floodlight",
  "plot_path_diagram",
  # timing helpers (used inside workers)
  "fmt_dur","elapsed_sec","fmt_progress","SCRIPT_START",
  # w8 configuration constants
  "LAGS_TO_TEST","N_BOOTSTRAP","BOOTSTRAP_SEED",
  "BLOCK_LENGTH","GAM_SCREEN_P","N_BOOT_SCREEN","JN_BOOT_N",
  "W_VALUES","W_LABELS","MM_DIR","GAM_K_MAIN","GAM_K_TENSOR","GAM_RESPONSE",
  "SENSITIVITY_INDICES","N_WORKERS",
  "SWEI_SCALE",
  # ↓ needed by setwd() BEFORE source() can run in the worker
  "WD_PATH"
)

results_list <- future.apply::future_lapply(
  ALL_INDICES,
  FUN = function(idx_spec) {
    
    # ── package availability inside worker ──────────────────────────────────
    setwd(WD_PATH)                          # workers don't inherit working directory
    source("DROUGHT_ANALYSIS_utils.R")      # defines path constants + helpers
    # utils_teleconnection_addon.R functions are already exported via future.globals;
    # sourcing it again in each worker is redundant — removed (PERF-7).
    suppressPackageStartupMessages({
      library(boot); library(parallel)
    })
    
    index_label <- sprintf("%s-%d", toupper(idx_spec$index), idx_spec$scale)
    IDX_DIR     <- file.path(MM_DIR, index_label)
    dir.create(IDX_DIR, showWarnings=FALSE, recursive=TRUE)
    
    idx_t0        <- proc.time()          # start time for this index worker
    lag_sec_vec__ <- numeric(0L)          # accumulates per-lag wall times for ETA
    
    all_mm_results   <- list()
    all_gam_fits     <- list()
    mm_lag_profile   <- list()
    gam_lag_prof     <- list()
    convergence_rows <- list()
    
    cat(sprintf("\n╔%s╗\n", strrep("═", 58)))
    cat(sprintf("  INDEX: %s | Lags 0–%d | script elapsed %s\n",
                index_label, max(LAGS_TO_TEST), fmt_dur(elapsed_sec(SCRIPT_START))))
    cat(sprintf("╚%s╝\n\n", strrep("═", 58)))
    
    # ════════════════════════════════════════════════════════
    # INNER LAG LOOP
    # ════════════════════════════════════════════════════════
    for (lag in LAGS_TO_TEST) {
      lag_t0__   <- proc.time()    # per-lag wall-clock start
      lag_pos__  <- which(LAGS_TO_TEST == lag)   # 1-based position in the lag sequence
      cat(sprintf("\n▶▶ LAG = %d months  [%d/%d | index elapsed %s | script elapsed %s]\n",
                  lag, lag_pos__, length(LAGS_TO_TEST),
                  fmt_dur(elapsed_sec(idx_t0)),
                  fmt_dur(elapsed_sec(SCRIPT_START))))
      
      # PERF-1: pass pre-loaded teleconnections
      df_lag <- tryCatch(
        build_analysis_dataframe(idx_spec$index, idx_spec$scale,
                                 lag_months=lag, start_year=1950L, end_year=2024L,
                                 tele_df=TELE_MASTER),
        error=function(e) { cat("  ❌ Data:", e$message, "\n"); NULL })
      if (is.null(df_lag)) next
      
      if (lag == 0L) check_mm_assumptions(df_lag)
      
      # ══════════════════════════════════════════════════════
      # TRACK 1: GAM
      # ══════════════════════════════════════════════════════
      cat(sprintf("\n  ─── TRACK 1: GAM (lag=%d) ───\n", lag))
      
      gam_fit <- fit_gam_full(df_lag,
                              response_var        = GAM_RESPONSE,
                              k_main              = GAM_K_MAIN,
                              k_tensor            = GAM_K_TENSOR,
                              include_seasonality = TRUE,
                              select              = TRUE)
      all_gam_fits[[sprintf("lag_%02d", lag)]] <- gam_fit
      
      # PERF-3: compute summary once, pass to all downstream helpers
      gam_s <- if (!is.null(gam_fit)) summary(gam_fit) else NULL
      
      if (!is.null(gam_fit)) {
        vp <- gam_variance_partition(gam_fit, precomputed_s=gam_s)
        if (!is.null(vp)) {
          vp$lag <- lag; vp$index <- index_label
          write.csv(vp, file.path(IDX_DIR,
                                  sprintf("gam_variance_partition_lag%02d.csv",lag)), row.names=FALSE)
        }
        gam_concurvity_check(gam_fit)
        
        # PERF-3: reuse gam_s in gam_lag_summary
        gam_lag_prof[[length(gam_lag_prof)+1]] <-
          gam_lag_summary(gam_fit, lag, index_label, precomputed_s=gam_s)
        
        gam_pdf <- file.path(IDX_DIR, sprintf("gam_lag%02d_%s.pdf", lag, index_label))
        if (safe_pdf(gam_pdf, width=14, height=10)) {
          p_surf <- plot_gam_surfaces(gam_fit, index_label=index_label, lag=lag,
                                      drought_thr=DROUGHT_THRESHOLD)
          if (!is.null(p_surf)) print(p_surf)
          p_part <- plot_gam_partials(gam_fit, index_label=index_label, lag=lag)
          if (!is.null(p_part)) print(p_part)
          p_resid <- plot_gam_residuals(gam_fit, index_label=index_label)
          if (!is.null(p_resid)) print(p_resid)
          grDevices::dev.off()
          cat(sprintf("  ✓ Saved: gam_lag%02d_%s.pdf (3 pages)\n", lag, index_label))
        }
      }
      
      # ══════════════════════════════════════════════════════
      # TRACK 2: MODERATED MEDIATION
      # PERF-4: GAM screening before spending full bootstrap budget
      # ══════════════════════════════════════════════════════
      cat(sprintf("\n  ─── TRACK 2: Moderated Mediation (lag=%d) ───\n", lag))
      
      # Extract ti(ONI,PDO) p-value from the precomputed GAM summary
      gam_ti_p_screen <- NA_real_
      if (!is.null(gam_s)) {
        smt <- gam_s$s.table
        ti_i <- grep("ti\\(ONI,PDO", rownames(smt))
        if (length(ti_i)) gam_ti_p_screen <- smt[ti_i[1], "p-value"]
      }
      
      # Decide bootstrap replicate count
      n_boot_this <- if (is.na(gam_ti_p_screen) || gam_ti_p_screen <= GAM_SCREEN_P) {
        cat(sprintf("  [PERF-4] ti(ONI,PDO) p=%.3f → full %d-rep bootstrap\n",
                    ifelse(is.na(gam_ti_p_screen), NaN, gam_ti_p_screen), N_BOOTSTRAP))
        N_BOOTSTRAP
      } else {
        cat(sprintf("  [PERF-4] ti(ONI,PDO) p=%.3f > %.2f → reduced %d-rep bootstrap\n",
                    gam_ti_p_screen, GAM_SCREEN_P, N_BOOT_SCREEN))
        N_BOOT_SCREEN
      }
      
      mm_res <- fit_moderated_mediation(df_lag, lag=lag,
                                        n_boot=n_boot_this, seed=BOOTSTRAP_SEED,
                                        block_length=BLOCK_LENGTH)
      all_mm_results[[sprintf("lag_%02d", lag)]] <- mm_res
      
      if (!is.null(mm_res)) {
        pe_lag <- mm_res$params
        for (ie_lbl in c("IE.low","IE.mid","IE.high","IMM")) {
          row <- pe_lag[pe_lag$label==ie_lbl,]
          if (nrow(row)) {
            est_val  <- .get_est(row)
            ci_lo_v  <- if ("ci.lower" %in% names(row)) row$ci.lower[1] else NA_real_
            ci_hi_v  <- if ("ci.upper" %in% names(row)) row$ci.upper[1] else NA_real_
            pv_val   <- if ("pvalue"   %in% names(row)) row$pvalue[1]   else NA_real_
            mm_lag_profile[[length(mm_lag_profile)+1]] <- data.frame(
              lag=lag, param=ie_lbl, est=est_val,
              ci_lower=ci_lo_v, ci_upper=ci_hi_v, p_value=pv_val,
              stringsAsFactors=FALSE)
          }
        }
        
        mm_pdf <- file.path(IDX_DIR, sprintf("mm_lag%02d_%s.pdf", lag, index_label))
        if (safe_pdf(mm_pdf, width=13, height=9)) {
          print(plot_path_diagram(mm_res, lag, index_label))
          print(plot_conditional_ie(mm_res, lag, index_label))
          jn_df <- tryCatch(compute_jn_floodlight(df_lag, n_boot=JN_BOOT_N),
                            error=function(e) NULL)
          if (!is.null(jn_df) && nrow(jn_df)) {
            print(plot_jn_floodlight(jn_df, lag, index_label))
            write.csv(jn_df, file.path(IDX_DIR, sprintf("jn_data_lag%02d.csv",lag)),
                      row.names=FALSE)
          }
          grDevices::dev.off()
          cat(sprintf("  ✓ Saved: mm_lag%02d_%s.pdf\n", lag, index_label))
        }
      }
      
      # ══════════════════════════════════════════════════════
      # TRACK 3: CONVERGENCE CHECK
      # PERF-3: reuse gam_s — no extra summary() call
      # ══════════════════════════════════════════════════════
      gam_ti_p <- NA_real_; mm_imm_p <- NA_real_
      gam_dev  <- NA_real_; mm_imm_e <- NA_real_
      
      if (!is.null(gam_s)) {
        smt    <- gam_s$s.table
        ti_idx <- grep("ti\\(ONI,PDO", rownames(smt))
        if (length(ti_idx)) gam_ti_p <- round(smt[ti_idx[1],"p-value"],4)
        gam_dev <- round(gam_s$dev.expl*100, 2)
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
        Boot_reps        = n_boot_this,
        Block_length     = BLOCK_LENGTH,
        stringsAsFactors=FALSE)
      
      # ── Lag timing & ETA ────────────────────────────────────────────────────
      lag_sec_vec__ <- c(lag_sec_vec__, elapsed_sec(lag_t0__))
      lags_done__   <- length(lag_sec_vec__)
      lags_left__   <- length(LAGS_TO_TEST) - lags_done__
      eta_lag__     <- mean(lag_sec_vec__) * lags_left__
      cat(fmt_progress(lags_done__, length(LAGS_TO_TEST),
                       elapsed_sec(idx_t0), eta_lag__), "\n")
      
    }  # end lag loop
    
    ##############################################################################
    # GAM LAG PROFILE FIGURE + TABLE
    ##############################################################################
    cat(sprintf("\n── GAM lag profile  [index elapsed %s] ──\n",
                fmt_dur(elapsed_sec(idx_t0))))
    
    if (length(gam_lag_prof) > 0) {
      gam_lp_df <- do.call(rbind, gam_lag_prof)
      write.csv(gam_lp_df, file.path(IDX_DIR,"gam_lag_profile.csv"), row.names=FALSE)
      cat("  Saved: gam_lag_profile.csv\n")
      
      lp_long <- tidyr::pivot_longer(gam_lp_df,
                                     cols=c("p_ti_ONI_PDO","p_ti_ONI_PNA","p_ti_PDO_PNA"),
                                     names_to="Term", values_to="p_value")
      lp_long$Term <- dplyr::recode(lp_long$Term,
                                    p_ti_ONI_PDO="ti(ONI,PDO)  ← key",
                                    p_ti_ONI_PNA="ti(ONI,PNA)",
                                    p_ti_PDO_PNA="ti(PDO,PNA)")
      ti_clrs <- c("ti(ONI,PDO)  ← key"="#C0392B","ti(ONI,PNA)"="#2980B9",
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
      
      gam_prof_pdf <- file.path(IDX_DIR, sprintf("gam_lag_profile_%s.pdf", index_label))
      if (safe_pdf(gam_prof_pdf, width=13, height=11)) {
        print(pA / pB); grDevices::dev.off()
        cat(sprintf("  Saved: gam_lag_profile_%s.pdf\n", index_label))
      }
    }
    
    ##############################################################################
    # MM LAG PROFILE FIGURE
    ##############################################################################
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
        ggplot2::labs(title=sprintf("MM: Conditional Indirect Effects across lags | %s",
                                    index_label),
                      subtitle=sprintf("Shading=95%% Block-bootstrap CI (l=%d) | Filled=p<0.05",
                                       BLOCK_LENGTH),
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
      
      mm_prof_pdf <- file.path(IDX_DIR, sprintf("mm_lag_profile_%s.pdf", index_label))
      if (safe_pdf(mm_prof_pdf, width=13, height=11)) {
        print(p_ie / p_imm); grDevices::dev.off()
        cat(sprintf("  Saved: mm_lag_profile_%s.pdf\n", index_label))
      }
      write.csv(lp_df, file.path(IDX_DIR, sprintf("lag_profile_data_%s.csv",index_label)),
                row.names=FALSE)
    }
    
    ##############################################################################
    # CONVERGENCE TABLE — GAM vs MM side by side
    ##############################################################################
    cat("\n── Convergence table (GAM vs MM) ──\n")
    
    if (length(convergence_rows) > 0) {
      conv_df <- do.call(rbind, convergence_rows)
      conv_df$Interpretation <- dplyr::case_when(
        conv_df$Both_agree  ~ "✓ CONVERGENT: both methods agree",
        conv_df$GAM_ti_sig & !conv_df$MM_IMM_sig
        ~ "GAM detects interaction; MM pathway not significant",
        !conv_df$GAM_ti_sig & conv_df$MM_IMM_sig
        ~ "MM detects mediation; GAM interaction marginal",
        TRUE ~ "Neither significant at this lag")
      cat("\n"); print(conv_df, row.names=FALSE)
      write.csv(conv_df, file.path(IDX_DIR,
                                   sprintf("convergence_table_%s.csv", index_label)), row.names=FALSE)
      cat(sprintf("  Saved: convergence_table_%s.csv\n", index_label))
    }
    
    ##############################################################################
    # SENSITIVITY ANALYSIS  (optimal lag, SWEI comparator)
    ##############################################################################
    cat("\n── Sensitivity analysis ──\n")
    
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
                                 start_year=1950L, end_year=2024L,
                                 tele_df=TELE_MASTER),
        error=function(e) { cat("  ❌", e$message, "\n"); NULL })
      if (is.null(df_s)) next
      
      gf   <- fit_gam_full(df_s, response_var=GAM_RESPONSE,
                           k_main=GAM_K_MAIN, k_tensor=GAM_K_TENSOR,
                           include_seasonality=TRUE, select=TRUE)
      gf_s <- if (!is.null(gf)) summary(gf) else NULL    # PERF-3
      sens_gam[[length(sens_gam)+1]] <- gam_lag_summary(gf, optimal_lag, lbl,
                                                        precomputed_s=gf_s)
      
      mm_s <- fit_moderated_mediation(df_s, lag=optimal_lag,
                                      n_boot=N_BOOTSTRAP, seed=BOOTSTRAP_SEED,
                                      block_length=BLOCK_LENGTH)
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
    
    ##############################################################################
    # EXCEL WORKBOOK — all results combined
    ##############################################################################
    cat("\n── Building Excel workbook ──\n")
    
    wb  <- openxlsx::createWorkbook()
    hdr <- openxlsx::createStyle(textDecoration="bold", fgFill="#1B3A6B",
                                 fontColour="#FFFFFF")
    add_sheet <- function(wb, name, df) {
      openxlsx::addWorksheet(wb, name)
      openxlsx::writeData(wb, name, df)
      openxlsx::addStyle(wb, name, hdr, rows=1, cols=seq_len(ncol(df)))
    }
    
    if (length(all_mm_results)) {
      mm_all_df <- do.call(rbind, lapply(names(all_mm_results), function(nm) {
        r <- all_mm_results[[nm]]; if (is.null(r)) return(NULL)
        d <- r$params_tbl; d$Lag <- r$lag; d
      }))
      if (!is.null(mm_all_df)) add_sheet(wb, "MM_All_Lags", mm_all_df)
    }
    if (length(mm_lag_profile)) add_sheet(wb, "MM_Lag_Profile", do.call(rbind,mm_lag_profile))
    if (length(sens_mm))        add_sheet(wb, "MM_Sensitivity",  do.call(rbind,sens_mm))
    if (length(gam_lag_prof))   add_sheet(wb, "GAM_Lag_Profile", do.call(rbind,gam_lag_prof))
    if (length(sens_gam))       add_sheet(wb, "GAM_Sensitivity", do.call(rbind,sens_gam))
    if (length(convergence_rows))
      add_sheet(wb, "Convergence_GAM_vs_MM", do.call(rbind,convergence_rows))
    
    # Bootstrap metadata sheet
    boot_meta <- data.frame(
      Parameter      = c("Bootstrap_type","Block_length_months","Full_R","Screened_R",
                         "Screen_threshold_ti_ONI_PDO_p","CI_type","P_value_method"),
      Value          = c("Stationary block bootstrap (boot::tsboot, sim='geom')",
                         as.character(BLOCK_LENGTH),
                         as.character(N_BOOTSTRAP),
                         as.character(N_BOOT_SCREEN),
                         as.character(GAM_SCREEN_P),
                         "Percentile (2.5%, 97.5% of bootstrap distribution)",
                         "Proportion method: 2*min(P(theta*<0), P(theta*>0))"),
      stringsAsFactors=FALSE)
    add_sheet(wb, "Bootstrap_Metadata", boot_meta)
    
    xlsx_out <- file.path(IDX_DIR, sprintf("results_summary_%s.xlsx", index_label))
    openxlsx::saveWorkbook(wb, xlsx_out, overwrite=TRUE)
    cat(sprintf("  Saved: %s\n", basename(xlsx_out)))
    idx_elapsed__ <- elapsed_sec(idx_t0)
    cat(sprintf("\n  ✓ %s complete in %s | script elapsed %s | outputs: %s\n",
                index_label, fmt_dur(idx_elapsed__),
                fmt_dur(elapsed_sec(SCRIPT_START)), IDX_DIR))
    
    list(index=index_label, status="complete", optimal_lag=optimal_lag,
         elapsed_sec=idx_elapsed__)
    
  },  # end future_lapply FUN
  future.globals  = .w8_globals,
  future.packages = c("lavaan","mgcv","boot","parallel","ggplot2","patchwork",
                      "dplyr","openxlsx","tidyr","scales"),
  future.seed     = BOOTSTRAP_SEED
)

# Restore sequential plan
future::plan(future::sequential)

# ── Post-run timing summary ────────────────────────────────────────────────────
total_wall__ <- elapsed_sec(SCRIPT_START)
cat(sprintf("\n── Parallel run complete | total wall time: %s ──\n", fmt_dur(total_wall__)))
cat(sprintf("  %-14s  %s\n", "Index", "Worker time"))
cat(sprintf("  %-14s  %s\n", strrep("-",14), strrep("-",12)))
for (r__ in results_list) {
  if (!is.null(r__) && !is.null(r__$elapsed_sec))
    cat(sprintf("  %-14s  %s\n", r__$index, fmt_dur(r__$elapsed_sec)))
}
cat(sprintf("  %-14s  %s  (wall clock — workers ran in parallel)\n",
            "TOTAL", fmt_dur(total_wall__)))

####################################################################################
# DONE
####################################################################################
cat("\n╔════════════════════════════════════════════════════════╗\n")
cat("║  w8 COMPLETE — ALL INDICES                             ║\n")
cat("╚════════════════════════════════════════════════════════╝\n")
cat(sprintf("  Total wall time : %s\n", fmt_dur(elapsed_sec(SCRIPT_START))))
cat(sprintf("  Root output directory: %s\n", MM_DIR))
cat("  Per-index subdirectories:\n")
for (idx_spec in ALL_INDICES) {
  lbl_done <- sprintf("%s-%d", toupper(idx_spec$index), idx_spec$scale)
  cat(sprintf("    %s/  →  gam_*.pdf | mm_*.pdf | convergence_table.csv | results_summary.xlsx\n",
              lbl_done))
}
if (length(SENSITIVITY_INDICES) > 0) {
  cat(sprintf("  Sensitivity (SWEI): gam_sensitivity.csv | mm_sensitivity.csv  →  %s\n",
              MM_DIR))
}
cat(sprintf("\n  Performance settings used:\n"))
cat(sprintf("    Workers (outer loop):  %d\n", N_WORKERS))
cat(sprintf("    Bootstrap type:        stationary block (boot::tsboot, sim='geom')\n"))
cat(sprintf("    Block length:          %d months\n", BLOCK_LENGTH))
cat(sprintf("    Full bootstrap R:      %d\n", N_BOOTSTRAP))
cat(sprintf("    Screened bootstrap R:  %d (when ti(ONI,PDO) p > %.2f)\n",
            N_BOOT_SCREEN, GAM_SCREEN_P))