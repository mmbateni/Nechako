####################################################################################
# w8_4_bayes_factor_table.R
# BAYES FACTOR TABLE — Cross-index nonparametric causal quantification
# Kretschmer et al. (2021) Part IV / Example 5
#
# Reads bf_optimal_lag_<index>.csv files written by w8_1_runIndex, and builds:
#   - A combined conditional probability table P(drought | ENSO, PDO)
#   - Bayes factors isolating the PDO amplification *given* ENSO phase
#   - A cross-index heatmap figure
#   - A summary Excel workbook
#
# Run AFTER w8_2_batch_moderated_mediation.R has finished.
#
# Usage:
#   source("w8_4_bayes_factor_table.R")
#   # or: Rscript w8_4_bayes_factor_table.R
#
# Outputs written to MM_DIR:
#   bf_cross_index_summary.xlsx     -- BF tables for all indices
#   bf_heatmap_enso_pdo.pdf         -- BF heatmap (ENSO x PDO phase) for each index
#   bf_delta_p_la_nina.pdf          -- delta_P (La Nina vs neutral) across indices
####################################################################################

BF_START <- proc.time()
source("w8_0_setup_moderated_mediation.R")
utils_load_packages(c("openxlsx","ggplot2","patchwork","dplyr","tidyr","scales"))

cat("\n================================================================\n")
cat("  w8_4  BAYES FACTOR TABLE — CROSS-INDEX COLLATION\n")
cat("================================================================\n")
cat(sprintf("  Started : %s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

all_labels <- vapply(ALL_INDICES, function(x)
  sprintf("%s-%d", toupper(x$index), x$scale), character(1))

bf_all  <- list()
found   <- character(0)
missing <- character(0)

####################################################################################
# COLLECT PER-INDEX BAYES FACTOR FILES
####################################################################################

for (lbl in all_labels) {
  bf_path <- file.path(MM_DIR, lbl, sprintf("bf_optimal_lag_%s.csv", lbl))
  if (!file.exists(bf_path)) {
    missing <- c(missing, lbl)
    cat(sprintf("  ! Missing BF file: %s  (run w8_2_batch first)\n", lbl)); next
  }
  bf_df <- tryCatch(read.csv(bf_path, stringsAsFactors=FALSE),
                    error=function(e) { cat(sprintf("  ! Could not read %s\n", bf_path)); NULL })
  if (!is.null(bf_df)) {
    bf_df$Index <- lbl
    bf_all[[lbl]] <- bf_df
    found <- c(found, lbl)
  }
}

cat(sprintf("\n  Found %d / %d indices.\n", length(found), length(all_labels)))
if (length(missing) > 0)
  cat(sprintf("  Missing: %s\n", paste(missing, collapse=", ")))
if (length(found) == 0)
  stop("No BF output files found.  Run w8_2_batch_moderated_mediation.R first.",
       call.=FALSE)

bf_combined <- do.call(rbind, bf_all)
rownames(bf_combined) <- NULL

####################################################################################
# FIGURE 1: Bayes factor heatmap — La Niña row (PDO amplification)
####################################################################################

cat("\n-- Building Bayes factor figures --\n")

ln_bf <- bf_combined[bf_combined$ENSO_cat == "LaNina" & !is.na(bf_combined$BayesFactor), ]

if (nrow(ln_bf) > 0) {
  p_heatmap <- ggplot2::ggplot(ln_bf,
      ggplot2::aes(x=PDO_phase, y=Index, fill=BayesFactor)) +
    ggplot2::geom_tile(colour="white", linewidth=0.5) +
    ggplot2::geom_text(ggplot2::aes(label=round(BayesFactor, 2)),
                       size=3.5, fontface="bold") +
    ggplot2::scale_fill_gradient2(low="#2980B9", mid="white", high="#C0392B",
                                  midpoint=1, name="Bayes\nFactor") +
    ggplot2::geom_vline(xintercept=1.5, colour="grey30", linewidth=0.3) +
    shared_ts_theme(12) +
    ggplot2::labs(
      title    = "Bayes Factors: P(drought | La Niña, PDO) / P(drought | La Niña)",
      subtitle = paste("Values > 1: PDO phase amplifies La Niña drought probability",
                       "| Drought threshold: SPEI-3 < -1.0"),
      x="PDO Phase", y="Drought Index")
  bf_heat_pdf <- file.path(MM_DIR, "bf_heatmap_enso_pdo.pdf")
  if (safe_pdf(bf_heat_pdf, width=9, height=7)) {
    print(p_heatmap); grDevices::dev.off()
    cat("  Saved: bf_heatmap_enso_pdo.pdf\n")
  }
}

####################################################################################
# FIGURE 2: Delta-P (La Niña vs ENSO-neutral) across all indices
# This is the nonparametric do(ENSO=LaNina) intervention estimate
####################################################################################

delta_p_rows <- lapply(found, function(lbl) {
  idx_dir <- file.path(MM_DIR, lbl)
  bf_path <- file.path(idx_dir, sprintf("bf_optimal_lag_%s.csv", lbl))
  if (!file.exists(bf_path)) return(NULL)
  bf_df <- read.csv(bf_path, stringsAsFactors=FALSE)
  p_ln   <- bf_df$ProbDrought[bf_df$ENSO_cat == "LaNina"][1]
  # P(drought | neutral) approximated from the Neutral row average
  p_neu  <- mean(bf_df$ProbDrought[bf_df$ENSO_cat == "Neutral"], na.rm=TRUE)
  data.frame(Index=lbl, P_LaNina=p_ln, P_Neutral=p_neu,
             DeltaP=p_ln - p_neu, stringsAsFactors=FALSE)
})
delta_p_df <- do.call(rbind, delta_p_rows[!vapply(delta_p_rows, is.null, logical(1))])

if (!is.null(delta_p_df) && nrow(delta_p_df) > 0) {
  delta_p_df$IndexType <- ifelse(grepl("^SPI-",  delta_p_df$Index), "SPI",
                                 ifelse(grepl("^SPEI-", delta_p_df$Index), "SPEI", "Other"))
  p_delta <- ggplot2::ggplot(delta_p_df,
      ggplot2::aes(x=reorder(Index, DeltaP), y=DeltaP, fill=IndexType)) +
    ggplot2::geom_hline(yintercept=0, linetype="dashed", colour="grey40") +
    ggplot2::geom_col(width=0.65) +
    ggplot2::scale_fill_manual(values=c("SPI"="#2980B9","SPEI"="#C0392B","Other"="#7F8C8D")) +
    ggplot2::coord_flip() +
    shared_ts_theme(12) +
    ggplot2::labs(
      title    = "Delta-P (La Niña vs ENSO-neutral) — causal intervention estimate",
      subtitle = "Nonparametric equivalent of do(ENSO = La Niña) across all drought indices",
      x="Index", y="P(drought|La Niña) - P(drought|neutral)")
  dp_pdf <- file.path(MM_DIR, "bf_delta_p_la_nina.pdf")
  if (safe_pdf(dp_pdf, width=10, height=7)) {
    print(p_delta); grDevices::dev.off()
    cat("  Saved: bf_delta_p_la_nina.pdf\n")
  }
  write.csv(delta_p_df, file.path(MM_DIR, "bf_delta_p_summary.csv"), row.names=FALSE)
  cat("  Saved: bf_delta_p_summary.csv\n")
}

####################################################################################
# COMBINED EXCEL WORKBOOK
####################################################################################

cat("\n-- Building Bayes factor Excel workbook --\n")

wb  <- openxlsx::createWorkbook()
hdr <- openxlsx::createStyle(textDecoration="bold", fgFill="#1B3A6B", fontColour="#FFFFFF")
.add <- function(wb, name, df) {
  name <- substr(name, 1, 31)
  openxlsx::addWorksheet(wb, name)
  openxlsx::writeData(wb, name, df)
  openxlsx::addStyle(wb, name, hdr, rows=1, cols=seq_len(ncol(df)))
}

.add(wb, "BF_All_Indices",       bf_combined)
if (!is.null(ln_bf) && nrow(ln_bf) > 0)
  .add(wb, "BF_LaNina_Focus",    ln_bf)
if (!is.null(delta_p_df) && nrow(delta_p_df) > 0)
  .add(wb, "BF_DeltaP_LaNina",   delta_p_df)

.add(wb, "BF_Metadata", data.frame(
  Parameter=c("DroughtThreshold","PDO_phase_rule","ENSO_tercile_method",
               "BayesFactor_formula","Indices_found","Collated_at"),
  Value=c("SPEI/SPI < -1.0 SD",
          "Standardised PDO >= 0 = Warm; < 0 = Cool",
          "Tercile quantile split of standardised ONI",
          "P(drought|ENSO,PDO) / P(drought|ENSO alone)",
          paste(found, collapse=", "),
          format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  stringsAsFactors=FALSE))

openxlsx::saveWorkbook(wb, file.path(MM_DIR, "bf_cross_index_summary.xlsx"),
                       overwrite=TRUE)
cat("  Saved: bf_cross_index_summary.xlsx\n")

####################################################################################
# DONE
####################################################################################

cat("\n================================================================\n")
cat("  w8_4 BAYES FACTOR COLLATION COMPLETE\n")
cat("================================================================\n")
cat(sprintf("  Time    : %s\n",  fmt_dur(elapsed_sec(BF_START))))
cat(sprintf("  Indices : %d / %d found\n\n", length(found), length(all_labels)))
cat("  Key outputs:\n")
cat("    bf_cross_index_summary.xlsx\n")
cat("    bf_heatmap_enso_pdo.pdf\n")
cat("    bf_delta_p_la_nina.pdf\n")
cat("    bf_delta_p_summary.csv\n\n")
cat("  Interpretation:\n")
cat("    BF > 1  = PDO phase amplifies drought probability given ENSO state\n")
cat("    delta_P = causal contribution of La Nina to severe drought probability\n\n")
