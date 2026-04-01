####################################################################################
# w8_03_collate.R
# POST-RUN COLLATION: merges all per-index outputs into one combined workbook
# and produces cross-index comparison figures.
#
# Run AFTER w8_02_batch.R has finished (fully or partially).
# Safe to run on partial results -- skips indices whose output folders are missing.
#
# Usage:
#   source("w8_03_collate.R")
#   # or: Rscript w8_03_collate.R
#
# Outputs written to MM_DIR:
#   combined_results_all_indices.xlsx   -- all per-index sheets in one workbook
#   cross_index_gam_deviance.pdf        -- GAM deviance across lags, all indices
#   cross_index_gam_ti_pvalues.pdf      -- ti(ONI,PDO) p-value heatmap
#   cross_index_mm_imm.pdf              -- IMM estimates across lags, all indices
#   cross_index_convergence_summary.csv -- convergence table, all indices combined
####################################################################################

COLLATE_START <- proc.time()
source("w8_0_setup_moderated_mediation.R")   # config + helpers (TELE_MASTER not needed here)

utils_load_packages(c("openxlsx","ggplot2","patchwork","dplyr","tidyr","scales"))

cat("\n================================================================\n")
cat("  w8  COLLATION -- CROSS-INDEX AGGREGATION\n")
cat("================================================================\n")
cat(sprintf("  Started : %s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

####################################################################################
# HELPER: read one CSV from an index output directory (returns NULL silently)
####################################################################################

.read_idx_csv <- function(index_label, filename) {
  path <- file.path(MM_DIR, index_label, filename)
  if (!file.exists(path)) return(NULL)
  tryCatch(read.csv(path, stringsAsFactors=FALSE),
           error=function(e) { cat(sprintf("  ! Could not read %s\n", path)); NULL })
}

####################################################################################
# COLLECT PER-INDEX DATA
####################################################################################

all_labels <- vapply(ALL_INDICES, function(x)
  sprintf("%s-%d", toupper(x$index), x$scale), character(1))

gam_lp_all   <- list()
mm_lp_all    <- list()
conv_all     <- list()
sens_gam_all <- list()
sens_mm_all  <- list()
found        <- character(0)
missing      <- character(0)

for (lbl in all_labels) {
  xlsx_path <- file.path(MM_DIR, lbl, sprintf("results_summary_%s.xlsx", lbl))
  if (!file.exists(xlsx_path)) {
    missing <- c(missing, lbl)
    cat(sprintf("  ! Missing: %s (run w8_2_batch_moderated_mediation.R first)\n", lbl))
    next
  }
  found <- c(found, lbl)
  
  g <- .read_idx_csv(lbl, "gam_lag_profile.csv")
  if (!is.null(g)) { g$Index <- lbl; gam_lp_all[[lbl]] <- g }
  
  m <- .read_idx_csv(lbl, sprintf("lag_profile_data_%s.csv", lbl))
  if (!is.null(m)) { m$Index <- lbl; mm_lp_all[[lbl]] <- m }
  
  cv <- .read_idx_csv(lbl, sprintf("convergence_table_%s.csv", lbl))
  if (!is.null(cv)) { cv$Index <- lbl; conv_all[[lbl]] <- cv }
}

# Sensitivity CSVs live in MM_DIR root
sg_path <- file.path(MM_DIR, "gam_sensitivity.csv")
sm_path <- file.path(MM_DIR, "mm_sensitivity.csv")
if (file.exists(sg_path)) sens_gam_all <- list(read.csv(sg_path, stringsAsFactors=FALSE))
if (file.exists(sm_path)) sens_mm_all  <- list(read.csv(sm_path, stringsAsFactors=FALSE))

cat(sprintf("\n  Found %d / %d indices.\n", length(found), length(all_labels)))
if (length(missing) > 0)
  cat(sprintf("  Missing: %s\n", paste(missing, collapse=", ")))
if (length(found) == 0) stop("No completed index outputs found.", call.=FALSE)

gam_lp_df <- if (length(gam_lp_all))  do.call(rbind, gam_lp_all)  else NULL
mm_lp_df  <- if (length(mm_lp_all))   do.call(rbind, mm_lp_all)   else NULL
conv_df   <- if (length(conv_all))    do.call(rbind, conv_all)     else NULL
sens_gam  <- if (length(sens_gam_all)) sens_gam_all[[1]]           else NULL
sens_mm   <- if (length(sens_mm_all))  sens_mm_all[[1]]            else NULL

####################################################################################
# CROSS-INDEX FIGURES
####################################################################################

cat("\n-- Building cross-index figures --\n")

# 1. GAM deviance explained across lags, all indices
if (!is.null(gam_lp_df) && "Dev_Expl_pct" %in% names(gam_lp_df)) {
  gam_lp_df$IndexType <- ifelse(grepl("^SPI-",  gam_lp_df$Index), "SPI",
                                ifelse(grepl("^SPEI-", gam_lp_df$Index), "SPEI", "Other"))
  p_dev <- ggplot2::ggplot(gam_lp_df,
                           ggplot2::aes(x=Lag, y=Dev_Expl_pct, colour=Index, group=Index)) +
    ggplot2::geom_line(linewidth=1.1) +
    ggplot2::geom_point(size=2.5) +
    ggplot2::facet_wrap(~IndexType, ncol=2, scales="free_y") +
    ggplot2::scale_x_continuous(breaks=LAGS_TO_TEST) +
    shared_ts_theme(11) +
    ggplot2::labs(
      title    = "GAM Deviance Explained (%) across Lags -- All Indices",
      subtitle = "Main + tensor smooths + seasonality | select=TRUE",
      x="Lag (months)", y="Deviance explained (%)", colour="Index")
  dev_pdf <- file.path(MM_DIR, "cross_index_gam_deviance.pdf")
  if (safe_pdf(dev_pdf, width=14, height=7)) {
    print(p_dev); grDevices::dev.off()
    cat("  Saved: cross_index_gam_deviance.pdf\n")
  }
}

# 2. ti(ONI,PDO) p-value heatmap
if (!is.null(gam_lp_df) && "p_ti_ONI_PDO" %in% names(gam_lp_df)) {
  p_heat <- ggplot2::ggplot(gam_lp_df,
                            ggplot2::aes(x=factor(Lag), y=Index,
                                         fill=pmin(p_ti_ONI_PDO, 0.5))) +
    ggplot2::geom_tile(colour="white", linewidth=0.4) +
    ggplot2::geom_text(
      ggplot2::aes(label=ifelse(!is.na(p_ti_ONI_PDO),
                                formatC(p_ti_ONI_PDO, digits=3, format="f"), ""),
                   colour=p_ti_ONI_PDO < 0.05),
      size=2.8, fontface="bold") +
    ggplot2::scale_fill_gradient(low="#C0392B", high="white",
                                 limits=c(0, 0.5), name="p-value") +
    ggplot2::scale_colour_manual(values=c("FALSE"="grey40","TRUE"="white"),
                                 guide="none") +
    shared_ts_theme(11) +
    ggplot2::labs(
      title    = "ti(ONI,PDO) p-value heatmap across Lags and Indices",
      subtitle = "Red = low p-value | white text = significant (p<0.05)",
      x="Lag (months)", y="Index")
  heat_pdf <- file.path(MM_DIR, "cross_index_gam_ti_pvalues.pdf")
  if (safe_pdf(heat_pdf, width=13, height=7)) {
    print(p_heat); grDevices::dev.off()
    cat("  Saved: cross_index_gam_ti_pvalues.pdf\n")
  }
}

# 3. IMM estimates across lags, all indices
if (!is.null(mm_lp_df) && "param" %in% names(mm_lp_df)) {
  imm_df <- mm_lp_df[mm_lp_df$param=="IMM",]
  if (nrow(imm_df) > 0) {
    p_imm <- ggplot2::ggplot(imm_df,
                             ggplot2::aes(x=lag, y=est, colour=Index, group=Index)) +
      ggplot2::geom_hline(yintercept=0, linetype="dashed", colour="grey40") +
      ggplot2::geom_ribbon(ggplot2::aes(ymin=ci_lower, ymax=ci_upper, fill=Index),
                           alpha=0.10, colour=NA) +
      ggplot2::geom_line(linewidth=1.0) +
      ggplot2::geom_point(ggplot2::aes(shape=p_value<0.05), size=2.5) +
      ggplot2::scale_shape_manual(values=c("FALSE"=1,"TRUE"=19), name="p<0.05") +
      ggplot2::scale_x_continuous(breaks=LAGS_TO_TEST) +
      ggplot2::facet_wrap(~Index, ncol=5) +
      shared_ts_theme(10) +
      ggplot2::labs(
        title    = "Index of Moderated Mediation (IMM) -- All Indices",
        subtitle = sprintf("Block-bootstrap CI (l=%d) | Filled point = p<0.05",
                           BLOCK_LENGTH),
        x="Lag (months)", y="IMM (standardised)") +
      ggplot2::theme(legend.position="none")
    imm_pdf <- file.path(MM_DIR, "cross_index_mm_imm.pdf")
    if (safe_pdf(imm_pdf, width=16, height=8)) {
      print(p_imm); grDevices::dev.off()
      cat("  Saved: cross_index_mm_imm.pdf\n")
    }
  }
}

# 4. Convergence summary CSV
if (!is.null(conv_df)) {
  write.csv(conv_df, file.path(MM_DIR, "cross_index_convergence_summary.csv"),
            row.names=FALSE)
  cat("  Saved: cross_index_convergence_summary.csv\n")
  cat(sprintf("\n  Convergence across %d index-lag combinations:\n", nrow(conv_df)))
  tbl <- table(conv_df$Interpretation)
  for (nm in names(tbl)) cat(sprintf("    [%3d]  %s\n", tbl[[nm]], nm))
}

####################################################################################
# COMBINED EXCEL WORKBOOK
####################################################################################

cat("\n-- Building combined Excel workbook --\n")

wb  <- openxlsx::createWorkbook()
hdr <- openxlsx::createStyle(textDecoration="bold", fgFill="#1B3A6B",
                             fontColour="#FFFFFF")
.add <- function(wb, name, df) {
  name <- substr(name, 1, 31)  # Excel sheet name limit
  openxlsx::addWorksheet(wb, name)
  openxlsx::writeData(wb, name, df)
  openxlsx::addStyle(wb, name, hdr, rows=1, cols=seq_len(ncol(df)))
}

if (!is.null(gam_lp_df)) .add(wb, "GAM_Lag_Profile_ALL",  gam_lp_df)
if (!is.null(mm_lp_df))  .add(wb, "MM_Lag_Profile_ALL",   mm_lp_df)
if (!is.null(conv_df))   .add(wb, "Convergence_ALL",       conv_df)
if (!is.null(sens_gam))  .add(wb, "GAM_Sensitivity",       sens_gam)
if (!is.null(sens_mm))   .add(wb, "MM_Sensitivity",        sens_mm)

# Per-index MM_All_Lags sheets pulled from individual workbooks
for (lbl in found) {
  xlsx_path <- file.path(MM_DIR, lbl, sprintf("results_summary_%s.xlsx", lbl))
  tryCatch({
    mm_sheet       <- openxlsx::read.xlsx(xlsx_path, sheet="MM_All_Lags")
    mm_sheet$Index <- lbl
    .add(wb, sprintf("MM_%s", lbl), mm_sheet)
  }, error=function(e) NULL)
}

# Run status and metadata sheets
.add(wb, "Run_Status", data.frame(
  Index  = all_labels,
  Status = ifelse(all_labels %in% found, "complete", "missing"),
  stringsAsFactors=FALSE))

.add(wb, "Collation_Info", data.frame(
  Parameter = c("Collated_at","Indices_found","Indices_missing",
                "Block_length","Full_Bootstrap_R","Screened_Bootstrap_R"),
  Value     = c(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                paste(found,   collapse=", "),
                paste(missing, collapse=", "),
                as.character(BLOCK_LENGTH),
                as.character(N_BOOTSTRAP),
                as.character(N_BOOT_SCREEN)),
  stringsAsFactors=FALSE))

openxlsx::saveWorkbook(wb,
                       file.path(MM_DIR, "combined_results_all_indices.xlsx"),
                       overwrite=TRUE)
cat("  Saved: combined_results_all_indices.xlsx\n")

####################################################################################
# DONE
####################################################################################

cat("\n================================================================\n")
cat("  w8 COLLATION COMPLETE\n")
cat("================================================================\n")
cat(sprintf("  Collation time : %s\n", fmt_dur(elapsed_sec(COLLATE_START))))
cat(sprintf("  Output root    : %s\n\n", MM_DIR))
cat("  Key outputs:\n")
cat("    combined_results_all_indices.xlsx\n")
cat("    cross_index_gam_deviance.pdf\n")
cat("    cross_index_gam_ti_pvalues.pdf\n")
cat("    cross_index_mm_imm.pdf\n")
cat("    cross_index_convergence_summary.csv\n\n")