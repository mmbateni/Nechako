####################################################################################
# w8_02_batch.R
# RESUME-AWARE PARALLEL ORCHESTRATOR
#
# USAGE:
#   source("w8_02_batch.R")
#
# Scans MM_DIR for already-completed indices (results_summary_*.xlsx present),
# submits only PENDING indices to future_lapply, then prints a timing table.
#
# RESUME: if interrupted, re-run this script — completed indices are skipped.
# FORCE FULL RE-RUN: set FORCE_RERUN <- TRUE before sourcing.
####################################################################################

rm(list=ls()); gc()
SCRIPT_START <- proc.time()
source("w8_0_setup_moderated_mediation.R")

FORCE_RERUN <- FALSE

# ── Resume check ─────────────────────────────────────────────────────────────────
.is_complete <- function(idx_spec) {
  lbl <- sprintf("%s-%d", toupper(idx_spec$index), idx_spec$scale)
  file.exists(file.path(MM_DIR, lbl, sprintf("results_summary_%s.xlsx", lbl)))
}

pending <- if (FORCE_RERUN) ALL_INDICES else Filter(Negate(.is_complete), ALL_INDICES)
skipped <- if (FORCE_RERUN) list()     else Filter(.is_complete,         ALL_INDICES)

n_total <- length(ALL_INDICES); n_done <- length(skipped); n_pending <- length(pending)

cat(sprintf("\n╔════════════════════════════════════════════════════════╗\n"))
cat(sprintf("║  w8  MODERATED MEDIATION + GAM TENSOR SMOOTHS          ║\n"))
cat(sprintf("╚════════════════════════════════════════════════════════╝\n\n"))
cat(sprintf("  Started  : %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat(sprintf("  Total    : %d indices\n", n_total))
cat(sprintf("  Complete : %d  (will be skipped)\n", n_done))
cat(sprintf("  Pending  : %d  (will be run now)\n\n", n_pending))

if (n_done > 0) {
  cat("  Already complete:\n")
  for (s in skipped) cat(sprintf("    [done]  %s-%d\n", toupper(s$index), s$scale))
  cat("\n")
}

if (n_pending == 0) {
  cat("  All indices already complete.  Run w8_03_collate.R to build combined outputs.\n\n")
  stop("All done — no pending indices.", call.=FALSE)
}

cat("  Pending:\n")
for (p in pending) cat(sprintf("    [ run]  %s-%d\n", toupper(p$index), p$scale))
cat("\n")

# ── Parallel execution ────────────────────────────────────────────────────────────
n_workers_actual <- min(n_pending, N_WORKERS)
cat(sprintf("-- Launching %d workers for %d pending indices | script elapsed %s --\n\n",
            n_workers_actual, n_pending, fmt_dur(elapsed_sec(SCRIPT_START))))

future::plan(future::sequential)
future::plan(future::multisession, workers=n_workers_actual)

results_list <- future.apply::future_lapply(
  pending,
  FUN = function(idx_spec) {
    setwd(WD_PATH)
    source("DROUGHT_ANALYSIS_utils.R")
    suppressPackageStartupMessages({ library(boot); library(parallel) })
    assign("idx_spec", idx_spec, envir=.GlobalEnv)
    source("w8_1_runIndex_moderated_mediation.R")
    .idx_result
  },
  future.globals  = .w8_globals,
  future.packages = c("lavaan","mgcv","boot","parallel","ggplot2","patchwork",
                      "dplyr","openxlsx","tidyr","scales"),
  future.seed     = BOOTSTRAP_SEED
)

future::plan(future::sequential)

# ── Post-run timing table ─────────────────────────────────────────────────────────
total_wall__ <- elapsed_sec(SCRIPT_START)
cat(sprintf("\n-- Parallel run complete | total wall time: %s --\n", fmt_dur(total_wall__)))
cat(sprintf("  %-14s  %-10s  %s\n","Index","Status","Worker time"))
cat(sprintf("  %-14s  %-10s  %s\n",strrep("-",14),strrep("-",10),strrep("-",12)))
for (s in skipped)
  cat(sprintf("  %-14s  %-10s  %s\n",
              sprintf("%s-%d",toupper(s$index),s$scale), "skipped", "--"))
for (r in results_list)
  if (!is.null(r))
    cat(sprintf("  %-14s  %-10s  %s\n", r$index,
                ifelse(isTRUE(r$status=="complete"),"complete","ERROR"),
                if (!is.null(r$elapsed_sec)) fmt_dur(r$elapsed_sec) else "??"))
cat(sprintf("  %-14s  %-10s  %s  (wall clock)\n","TOTAL","",fmt_dur(total_wall__)))

# ── Done banner ───────────────────────────────────────────────────────────────────
n_newly_done <- sum(vapply(results_list,function(r) isTRUE(!is.null(r)&&r$status=="complete"),logical(1)))
n_errors     <- n_pending - n_newly_done

cat(sprintf("\n╔════════════════════════════════════════════════════════╗\n"))
cat(sprintf("║  w8_02 BATCH COMPLETE                                   ║\n"))
cat(sprintf("╚════════════════════════════════════════════════════════╝\n"))
cat(sprintf("  Total wall time  : %s\n",  fmt_dur(total_wall__)))
cat(sprintf("  Indices complete : %d / %d\n", n_done+n_newly_done, n_total))
if (n_errors > 0) cat(sprintf("  !! Errors        : %d\n", n_errors))
cat(sprintf("  Output root      : %s\n\n", MM_DIR))
if (n_done+n_newly_done == n_total) {
  cat("  All complete. Run w8_03_collate.R to build combined outputs.\n\n")
} else {
  cat("  Some indices still pending.  Re-run w8_02_batch.R to resume.\n\n")
}