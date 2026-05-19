#!/usr/bin/env Rscript
# w8_5_pcmci_lag_validation.R
# PCMCI CAUSAL DISCOVERY — Optimal lag detection for Nechako teleconnection DAG
# Translated from Python (w8_5_pcmci_lag_validation.py)
# Requires: R (>=4.0), dplyr, argparse, reticulate
# Python env must have: pip install tigramite matplotlib numpy pandas

library(argparse)
library(dplyr)
library(tidyr)
library(reticulate)

# ─────────────────────────────────────────────────────────────────────────────
# CLI Argument Parsing
# ─────────────────────────────────────────────────────────────────────────────
parser <- ArgumentParser(
  description = "PCMCI causal discovery for Nechako teleconnections"
)
parser$add_argument(
  "--data_path", default = NULL,
  help = "Path to pcmci_input_data.csv. Defaults to ./moderated_mediation/pcmci_input_data.csv"
)
parser$add_argument(
  "--tau_max", type = "integer", default = 12L,
  help = "Maximum lag in months to test (default: 12)"
)
parser$add_argument(
  "--alpha", type = "double", default = 0.05,
  help = "PC-step significance threshold (default: 0.05)"
)
parser$add_argument(
  "--out_dir", default = NULL,
  help = "Output directory. Defaults to the data file's parent directory."
)
args <- parser$parse_args()

# ─────────────────────────────────────────────────────────────────────────────
# Locate Input & Setup Output Directory
# ─────────────────────────────────────────────────────────────────────────────
if (is.null(args$data_path)) {
  args$data_path <- file.path("moderated_mediation", "pcmci_input_data.csv")
}

if (!file.exists(args$data_path)) {
  stop(
    "ERROR: Input data not found: ", args$data_path, "\n",
    "Export from R first:\n",
    "  write.csv(pcmci_df, file.path(MM_DIR, 'pcmci_input_data.csv'), row.names=FALSE)"
  )
}

out_dir <- args$out_dir %||% dirname(normalizePath(args$data_path))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# Load and Validate Data
# ─────────────────────────────────────────────────────────────────────────────
df <- read.csv(args$data_path, stringsAsFactors = FALSE)
if ("date" %in% names(df)) df$date <- as.Date(df$date)

var_names <- c("ONI", "PDO", "PNA", "Z500_ridge", "SPEI3", "F_thm")
missing_cols <- setdiff(var_names, names(df))
if (length(missing_cols) > 0) {
  stop("ERROR: Missing columns in input data: ", paste(missing_cols, collapse = ", "), "\nCheck the R export step.")
}

# Use only complete cases across all variables
data_arr <- as.matrix(df[, var_names])
data_arr <- data_arr[complete.cases(data_arr), , drop = FALSE]
n_obs <- nrow(data_arr)
n_vars <- ncol(data_arr)
cat(sprintf("PCMCI input: %d obs x %d variables.\n", n_obs, n_vars))
cat(sprintf("Variables:   %s\n\n", paste(var_names, collapse = ", ")))

# ─────────────────────────────────────────────────────────────────────────────
# Import tigramite via reticulate
# ─────────────────────────────────────────────────────────────────────────────
tryCatch({
  # Ensure Python environment has tigramite
  use_python(reticulate::py_config()$python, required = TRUE)
  pp       <- import("tigramite.data_processing")
  PCMCI    <- import("tigramite.pcmci")$PCMCI
  ParCorr  <- import("tigramite.independence_tests.parcorr")$ParCorr
}, error = function(e) {
  stop("ERROR: tigramite not installed in Python environment.\nInstall with: pip install tigramite")
})

has_plot <- tryCatch({
  plt <- import("matplotlib")
  plt$use("Agg")
  import("matplotlib.pyplot")
  import("tigramite.plotting")
  TRUE
}, error = function(e) {
  message("Warning: matplotlib not available; skipping causal graph figure.\n")
  FALSE
})

# ─────────────────────────────────────────────────────────────────────────────
# Run PCMCI
# ─────────────────────────────────────────────────────────────────────────────
cat(sprintf("Running PCMCI (tau_max=%d, alpha=%.2f) ...\n", args$tau_max, args$alpha))
tg_df   <- pp$DataFrame(data_arr, var_names = var_names)
pcmci   <- PCMCI(dataframe = tg_df, cond_ind_test = ParCorr(), verbosity = 1L)
results <- pcmci$run_pcmci(tau_max = args$tau_max, pc_alpha = args$alpha)

# ─────────────────────────────────────────────────────────────────────────────
# FDR-corrected p-values & Extract Significant Links
# ─────────────────────────────────────────────────────────────────────────────
# Convert numpy arrays to R arrays for safe 1-based indexing
p_matrix  <- reticulate::py_to_r(results$p_matrix)
val_matrix <- reticulate::py_to_r(results$val_matrix)

q_matrix <- pcmci$get_corrected_pvalues(
  p_matrix = p_matrix,
  tau_max = args$tau_max,
  fdr_method = "fdr_bh"
)
q_matrix <- reticulate::py_to_r(q_matrix)

# Extract links
sig_rows <- list()
for (j in seq_along(var_names)) {
  target <- var_names[j]
  for (i in seq_along(var_names)) {
    source <- var_names[i]
    if (i == j) next  # skip auto-links
    for (tau in seq_len(args$tau_max)) {
      q_val <- q_matrix[i, j, tau]
      coef  <- val_matrix[i, j, tau]
      if (q_val < args$alpha) {
        sig_rows[[length(sig_rows) + 1]] <- data.frame(
          Source  = source,
          Target  = target,
          Lag     = tau,
          Coeff   = round(coef, 4),
          q_value = round(q_val, 4),
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

sig_df <- if (length(sig_rows) == 0) {
  tibble()
} else {
  bind_rows(sig_rows) %>% arrange(Target, Source, Lag)
}

sig_path <- file.path(out_dir, "pcmci_results_summary.csv")
write.csv(sig_df, sig_path, row.names = FALSE)
cat(sprintf("\nSignificant links found: %d\n  Saved: %s\n", nrow(sig_df), sig_path))

# ─────────────────────────────────────────────────────────────────────────────
# Optimal Lag per Link (max |coeff|)
# ─────────────────────────────────────────────────────────────────────────────
if (nrow(sig_df) > 0) {
  opt_df <- sig_df %>%
    group_by(Source, Target) %>%
    arrange(desc(abs(Coeff))) %>%
    slice(1) %>%
    ungroup() %>%
    rename(OptimalLag = Lag, Coeff_at_OptLag = Coeff) %>%
    select(Source, Target, OptimalLag, Coeff_at_OptLag, q_value)
} else {
  opt_df <- tibble()
}

opt_path <- file.path(out_dir, "pcmci_optimal_lags.csv")
write.csv(opt_df, opt_path, row.names = FALSE)
cat(sprintf("  Saved: %s\n", opt_path))

# ─────────────────────────────────────────────────────────────────────────────
# Validate against proposed Nechako DAG
# ─────────────────────────────────────────────────────────────────────────────
dag_links <- data.frame(
  Source = c("ENSO", "PDO", "PDO", "PNA", "Z500_ridge", "F_thm"),
  Target = c("PNA", "PNA", "SPEI3", "Z500_ridge", "SPEI3", "SPEI3"),
  Label  = c("a-path", "PDO->PNA modulation", "PDO direct path",
             "PNA->Z500 chain", "proximate atm driver", "thermodynamic confounder"),
  stringsAsFactors = FALSE
)

dag_links_adapted <- dag_links %>% mutate(Source = ifelse(Source == "ENSO", "ONI", Source))

val_rows <- list()
cat("\n  DAG link validation (PCMCI vs proposed DAG): \n")
cat(sprintf("  %-30s %-6s %-8s %-8s %-8s\n", "Link", "DAG?", "PCMCI?", "OptLag", "Coeff"))
cat(paste0("  ", strrep("-", 62), "\n"))

for (k in seq_len(nrow(dag_links_adapted))) {
  row_dag <- dag_links_adapted[k, ]
  in_dag <- TRUE
  
  match_df <- if (nrow(opt_df) > 0) {
    opt_df %>% filter(Source == row_dag$Source, Target == row_dag$Target)
  } else tibble()
  
  if (nrow(match_df) > 0) {
    opt_lag  <- match_df$OptimalLag[1]
    coeff    <- match_df$Coeff_at_OptLag[1]
    detected <- TRUE
    status   <- "OK  "
  } else {
    opt_lag  <- NA
    coeff    <- NA
    detected <- FALSE
    status   <- "--  "
  }
  
  link_str <- sprintf("%s -> %s", row_dag$Source, row_dag$Target)
  cat(sprintf("  %s %-28s %-6s %-8s %-8s %s\n",
              status, link_str, "Yes",
              ifelse(detected, "Yes", "No"),
              ifelse(detected, as.character(opt_lag), "n/a"),
              ifelse(detected, sprintf("%.3f", coeff), "n/a")))
  
  val_rows[[k]] <- tibble(
    Source    = row_dag$Source,
    Target    = row_dag$Target,
    Label     = row_dag$Label,
    In_DAG    = in_dag,
    PCMCI_sig = detected,
    OptLag    = opt_lag,
    Coeff     = coeff
  )
}

val_df <- bind_rows(val_rows)

# Report unexpected links
if (nrow(opt_df) > 0) {
  dag_pairs <- paste(dag_links_adapted$Source, dag_links_adapted$Target, sep = "->")
  opt_pairs <- paste(opt_df$Source, opt_df$Target, sep = "->")
  unexpected <- opt_df[!opt_pairs %in% dag_pairs, ]
} else {
  unexpected <- tibble()
}

if (nrow(unexpected) > 0) {
  cat("\n  Unexpected significant links (not in proposed DAG): \n")
  for (i in seq_len(nrow(unexpected))) {
    cat(sprintf("    %s -> %s  lag=%d coeff=%.3f\n",
                unexpected$Source[i], unexpected$Target[i],
                unexpected$OptimalLag[i], unexpected$Coeff_at_OptLag[i]))
  }
  cat("  -> Consider revising the DAG in w8_0b_dag_specification.R\n")
} else {
  cat("\n  No unexpected significant links outside the proposed DAG.\n")
}

val_path <- file.path(out_dir, "pcmci_dag_validation.csv")
write.csv(val_df, val_path, row.names = FALSE)
cat(sprintf("\n  Saved: %s\n", val_path))

# ─────────────────────────────────────────────────────────────────────────────
# Key Lag Recommendations
# ─────────────────────────────────────────────────────────────────────────────
cat("\n  Key lag recommendations for updating the Nechako DAG: \n")
key_pairs <- list(
  list(src = "ONI", tgt = "PNA",     label = "ONI -> PNA (expected: 0-1 months)"),
  list(src = "PNA", tgt = "SPEI3",   label = "PNA -> SPEI3 (expected: 1-3 months)"),
  list(src = "PDO", tgt = "SPEI3",   label = "PDO -> SPEI3 direct (check necessity)")
)

for (kp in key_pairs) {
  if (nrow(opt_df) > 0) {
    m <- opt_df %>% filter(Source == kp$src, Target == kp$tgt)
    lag_str <- if (nrow(m) > 0) sprintf("optimal lag = %d months", m$OptimalLag[1]) else "not significant"
  } else lag_str <- "not significant"
  cat(sprintf("    %s\n      PCMCI result: %s\n", kp$label, lag_str))
}

# ─────────────────────────────────────────────────────────────────────────────
# Causal Graph Figure (via reticulate)
# ─────────────────────────────────────────────────────────────────────────────
if (has_plot && !is.null(results$graph)) {
  tryCatch({
    fig_path <- file.path(out_dir, "pcmci_graph.png")
    plt      <- import("matplotlib.pyplot")
    tp       <- import("tigramite.plotting")
    
    fig <- plt$figure(figsize = c(10, 8))
    tp$plot_graph(
      val_matrix          = val_matrix,
      graph               = results$graph,
      var_names           = var_names,
      link_colorbar_label = "MCI coefficient",
      node_colorbar_label = "Auto-MCI",
      show_colorbar       = TRUE
    )
    plt$suptitle(
      sprintf("PCMCI Causal Graph | tau_max=%d | alpha=%.2f", args$tau_max, args$alpha),
      fontsize = 12
    )
    plt$savefig(fig_path, dpi = 150, bbox_inches = "tight")
    plt$close(fig)
    cat(sprintf("  Saved: %s\n", fig_path))
  }, error = function(e) {
    cat(sprintf("  Warning: could not save causal graph (%s)\n", e$message))
  })
}

# ─────────────────────────────────────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────────────────────────────────────
cat("\n", strrep("= ", 30), "\n")
cat("  PCMCI validation complete.\n")
cat(strrep("= ", 30), "\n")
cat(sprintf("  Outputs in: %s\n", out_dir))
cat("    pcmci_results_summary.csv\n")
cat("    pcmci_optimal_lags.csv\n")
cat("    pcmci_dag_validation.csv\n")
if (has_plot) cat("    pcmci_graph.png\n")
cat("\n  Next step: Update lag structure in w8_0b_dag_specification.R\n")
cat("  if PCMCI optimal lags differ from defaults (0-3 months).\n")