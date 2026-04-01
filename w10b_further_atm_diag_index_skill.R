# ==============================================================================
# w10b_index_skill.R
# B2: CLIMATE INDEX CORRELATIONS & PREDICTABILITY SKILL  (SPI & SPEI)
#
# Split from w10_further_atmospheric_diagnostics.R.
# This module has NO dependency on ERA5 rasters — it only needs:
#   • pre-computed basin-averaged SPI/SPEI CSVs  (from w1/w3 pipeline)
#   • NOAA climate indices — downloaded & cached automatically via
#     load_teleconnection() in utils_teleconnection_addon.R
#     (PDO, PNA from w7 cache; AMO, Nino34, AO downloaded on first run)
#
# PERF: SPI/SPEI basin means are loaded from the pre-computed CSVs via
#       load_basin_avg_csv() (in DROUGHT_ANALYSIS_utils.R) instead of reading
#       12 per-month NetCDF files and computing basin means from scratch.
#       This cuts the I/O for this module from ~60 NC reads to ~10 CSV reads.
#
# Outputs:
#   skill_results/
#     lag_correlation_matrix_SPI_SPEI.csv
#     best_lag_summary_SPI_SPEI.csv
#     loyo_skill_scores_SPI_SPEI.csv
#     fig_lag_correlation_profiles_SPI_SPEI.pdf
#     fig_skill_summary_SPI_SPEI.pdf
#     fig_best_lag_comparison.pdf
#     fig_hit_rate_comparison.pdf
#     summary_SPI_SPEI_skill.csv
#
# Run AFTER: w1/w3 (SPI/SPEI basin CSVs must exist)
# Independent of w10a — can run concurrently. Needs only the SPI/SPEI basin-averaged
# CSVs from w1/w2, and downloads teleconnection indices automatically
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(lubridate)
})

setwd("D:/Nechako_Drought/Nechako/")
source("DROUGHT_ANALYSIS_utils.R")          # provides load_basin_avg_csv, path constants
source("utils_teleconnection_addon.R")      # provides load_teleconnection(), TELE_DIR

out_b2 <- "skill_results"
if (!dir.exists(out_b2)) dir.create(out_b2)

cat("\n========== w10b: CLIMATE INDEX SKILL ANALYSIS ==========\n")

# ------------------------------------------------------------------------------
# LOAD SPI/SPEI BASIN MEANS via pre-computed CSVs
# load_basin_avg_csv() reads the wide Year × Month CSV and returns data.frame(date, value)
# This replaces the expensive per-month NC load used in the original w10.
# ------------------------------------------------------------------------------
cat("\nLoading SPI & SPEI basin means from pre-computed CSVs...\n")

scales      <- c(1, 3, 6, 12, 24)
index_types <- c("SPI", "SPEI")

spi_spei_df <- dplyr::bind_rows(lapply(index_types, function(it) {
  dplyr::bind_rows(lapply(scales, function(sc) {
    df <- tryCatch(
      load_basin_avg_csv(tolower(it), sc),
      error = function(e) NULL)
    if (is.null(df) || nrow(df) == 0) {
      cat(sprintf("  ⚠ %s-%d: not found, skipped\n", it, sc))
      return(NULL)
    }
    df %>%
      rename(Index_Value = value) %>%
      mutate(scale = sc, year = year(date), month_num = month(date),
             index_type = it,
             index = sprintf("%s%d", it, sc))
  }))
}))

cat(sprintf("✓ SPI & SPEI loaded: %d records (%d index × scale combinations)\n",
            nrow(spi_spei_df),
            length(unique(paste(spi_spei_df$index_type, spi_spei_df$scale)))))

# ------------------------------------------------------------------------------
# LOAD CLIMATE INDICES via load_teleconnection()
# PDO and PNA reuse the cache written by w7 (teleconnections/).
# AMO, Nino34, and AO are downloaded and cached on first run.
# All six indices share the same TELE_DIR cache — no manual file placement needed.
# ------------------------------------------------------------------------------
cat("\nLoading climate indices via load_teleconnection()...\n")

.load_idx <- function(name, col) {
  df <- tryCatch(
    load_teleconnection(name, start_year = 1950, end_year = 2025),
    error = function(e) { cat(sprintf("  ⚠ %s: %s\n", toupper(name), e$message)); NULL })
  if (is.null(df)) return(NULL)
  df %>%
    rename(!!col := value) %>%
    mutate(year = as.integer(format(date, "%Y")),
           month_num = as.integer(format(date, "%m")))
}

idx_list <- Filter(Negate(is.null), list(
  .load_idx("pdo",    "PDO"),
  .load_idx("pna",    "PNA"),
  .load_idx("amo",    "AMO"),
  .load_idx("nino34", "Nino34"),
  .load_idx("ao",     "AO")
))

if (!length(idx_list)) {
  cat("WARNING: No climate indices could be loaded – B2 skipped.\n")
  quit(save = "no")
}

clim_idx_df <- Reduce(function(a, b) full_join(a, b, by = c("date","year","month_num")),
                      idx_list)
index_names <- setdiff(colnames(clim_idx_df), c("date","year","month_num"))
cat(sprintf("✓ Climate indices: %s\n", paste(index_names, collapse = ", ")))

# ------------------------------------------------------------------------------
# B2-A: LAG CORRELATIONS
# ------------------------------------------------------------------------------
cat("\n── B2-A: Lag correlations ──\n")
max_lag <- 12

lagged_cor <- function(index_vals, idx_vals, max_lag = 12) {
  n <- length(index_vals)
  dplyr::bind_rows(lapply(0:max_lag, function(lag) {
    if (lag == 0) { s <- index_vals; x <- idx_vals
    } else        { s <- index_vals[(lag+1):n]; x <- idx_vals[1:(n-lag)] }
    ok <- !is.na(s) & !is.na(x)
    if (sum(ok) < 20) return(data.frame(lag_k=lag, r=NA_real_, p=NA_real_, n=sum(ok)))
    ct <- cor.test(s[ok], x[ok])
    data.frame(lag_k=lag, r=ct$estimate, p=ct$p.value, n=sum(ok))
  }))
}

lag_cor_results <- list()
for (it in index_types) {
  for (sc in scales) {
    index_sc <- spi_spei_df %>% filter(index_type == it, scale == sc) %>% arrange(date)
    if (!nrow(index_sc)) next
    for (idx_nm in index_names) {
      merged <- index_sc %>%
        left_join(clim_idx_df %>% select(date, all_of(idx_nm)), by = "date")
      lc <- lagged_cor(merged$Index_Value, merged[[idx_nm]], max_lag)
      lc$index      <- idx_nm
      lc$scale      <- sc
      lc$index_type <- it
      lag_cor_results[[paste(it, sc, idx_nm)]] <- lc
    }
  }
}
lag_cor_df <- dplyr::bind_rows(lag_cor_results)
write_csv(lag_cor_df, file.path(out_b2,"lag_correlation_matrix_SPI_SPEI.csv"))

best_lag_summary <- lag_cor_df %>%
  group_by(index_type, scale, index) %>%
  slice_max(abs(r), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(significant = p < 0.05)
write_csv(best_lag_summary, file.path(out_b2,"best_lag_summary_SPI_SPEI.csv"))

# ------------------------------------------------------------------------------
# B2-B: LOYO CROSS-VALIDATION SKILL
# ------------------------------------------------------------------------------
cat("── B2-B: LOYO cross-validation ──\n")

loyo_skill <- function(index_vals, idx_vals, dates, lag_k) {
  n <- length(index_vals)
  if (lag_k > 0) {
    target    <- index_vals[(lag_k+1):n]
    predictor <- idx_vals[1:(n-lag_k)]
    yr_target <- year(dates[(lag_k+1):n])
  } else {
    target <- index_vals; predictor <- idx_vals; yr_target <- year(dates)
  }
  ok <- !is.na(target) & !is.na(predictor)
  target <- target[ok]; predictor <- predictor[ok]; yr_target <- yr_target[ok]
  years_unique <- unique(yr_target)
  if (length(years_unique) < 10) return(NULL)
  
  predictions <- numeric(length(target))
  for (yr in years_unique) {
    tr  <- yr_target != yr; te <- yr_target == yr
    if (sum(tr) < 10 || !any(te)) { predictions[te] <- NA; next }
    m   <- lm(target[tr] ~ predictor[tr])
    predictions[te] <- coef(m)[1] + coef(m)[2] * predictor[te]
  }
  ok2 <- !is.na(predictions)
  if (sum(ok2) < 5) return(NULL)
  r_cv  <- cor(target[ok2], predictions[ok2])
  rmse  <- sqrt(mean((target[ok2] - predictions[ok2])^2))
  d_obs <- target[ok2] < -1.282; d_prd <- predictions[ok2] < -1.282
  data.frame(r_cv = r_cv, rmse = rmse,
             hit_rate       = sum(d_obs & d_prd,  na.rm=TRUE) / max(sum(d_obs), 1),
             false_alarm_rate = sum(!d_obs & d_prd, na.rm=TRUE) / max(sum(!d_obs), 1))
}

target_predictors <- best_lag_summary %>%
  group_by(index_type, scale) %>%
  slice_max(abs(r), n = 2, with_ties = FALSE) %>%
  ungroup()

skill_results <- dplyr::bind_rows(lapply(seq_len(nrow(target_predictors)), function(i) {
  row <- target_predictors[i, ]
  merged <- spi_spei_df %>%
    filter(index_type == row$index_type, scale == row$scale) %>%
    arrange(date) %>%
    left_join(clim_idx_df %>% select(date, all_of(row$index)), by = "date")
  sk <- loyo_skill(merged$Index_Value, merged[[row$index]], merged$date, row$lag_k)
  if (is.null(sk)) return(NULL)
  sk$index_type <- row$index_type; sk$scale <- row$scale
  sk$index <- row$index; sk$lag_k <- row$lag_k
  sk
}))
write_csv(skill_results, file.path(out_b2,"loyo_skill_scores_SPI_SPEI.csv"))

# ------------------------------------------------------------------------------
# B2-C: FIGURES
# ------------------------------------------------------------------------------
cat("── B2-C: Figures ──\n")

pdf(file.path(out_b2,"fig_lag_correlation_profiles_SPI_SPEI.pdf"), width=14, height=10)
print(
  lag_cor_df %>%
    filter(scale %in% c(3,6,12), !is.na(r)) %>%
    mutate(sig = ifelse(p < 0.05,"p<0.05","n.s."),
           label = sprintf("%s-%d", index_type, scale)) %>%
    ggplot(aes(x=lag_k, y=r, colour=index, linetype=sig)) +
    geom_hline(yintercept=0, colour="gray50") +
    geom_line() + geom_point(aes(shape=sig), size=2) +
    facet_grid(index_type ~ label, scales="free_y") +
    scale_colour_brewer(palette="Dark2") +
    labs(title="Lag Correlation: SPI & SPEI vs. Climate Indices",
         x="Lag (months)", y="Pearson r",
         colour="Climate Index", linetype="Significance") +
    theme_minimal(base_size=11) + theme(legend.position="bottom")
)
dev.off()

pdf(file.path(out_b2,"fig_skill_summary_SPI_SPEI.pdf"), width=10, height=7)
print(
  skill_results %>%
    mutate(predictor_label = sprintf("%s @ lag-%d", index, lag_k),
           index_type = factor(index_type, levels=c("SPI","SPEI"))) %>%
    ggplot(aes(x=predictor_label, y=r_cv, fill=index_type)) +
    geom_bar(stat="identity", position="dodge") +
    geom_hline(yintercept=0) +
    geom_hline(yintercept=0.3, linetype="dashed", colour="red", alpha=0.5) +
    coord_flip() + scale_fill_manual(values=c(SPI="#2E86AB",SPEI="#A23B72")) +
    labs(title="LOYO Cross-Validation Skill (r_cv): SPI vs SPEI",
         x="Predictor (Climate Index @ Lag)", y="Cross-validated correlation r_cv") +
    theme_minimal(base_size=11) + theme(legend.position="top")
)
dev.off()

pdf(file.path(out_b2,"fig_best_lag_comparison.pdf"), width=10, height=6)
print(
  best_lag_summary %>%
    filter(scale %in% c(6,12)) %>%
    mutate(index_type=factor(index_type,levels=c("SPI","SPEI"))) %>%
    ggplot(aes(x=index, y=abs(r), fill=index_type)) +
    geom_bar(stat="identity", position="dodge") +
    geom_text(aes(label=paste0("lag=",lag_k)),
              position=position_dodge(width=0.9), vjust=-0.5, size=3) +
    facet_wrap(~paste0("Scale-",scale), ncol=2) +
    scale_fill_manual(values=c(SPI="#2E86AB",SPEI="#A23B72")) +
    labs(title="Best Lag Correlation Comparison: SPI vs SPEI",
         x="Climate Index", y="|Correlation| (best lag)") +
    theme_minimal(base_size=11) + theme(legend.position="bottom")
)
dev.off()

pdf(file.path(out_b2,"fig_hit_rate_comparison.pdf"), width=9, height=6)
print(
  skill_results %>%
    mutate(index_type=factor(index_type,levels=c("SPI","SPEI"))) %>%
    ggplot(aes(x=paste0(index,"-lag",lag_k), y=hit_rate, fill=index_type)) +
    geom_bar(stat="identity", position="dodge") +
    geom_hline(yintercept=0.5, linetype="dashed", colour="red", alpha=0.5) +
    coord_flip() + scale_fill_manual(values=c(SPI="#2E86AB",SPEI="#A23B72")) +
    labs(title="Drought Hit Rate: SPI vs SPEI Prediction Skill",
         x="Predictor", y="Hit Rate") +
    theme_minimal(base_size=11) + theme(legend.position="top")
)
dev.off()

write_csv(
  skill_results %>%
    select(index_type, scale, index, lag_k, r_cv, rmse, hit_rate, false_alarm_rate) %>%
    arrange(index_type, scale, desc(r_cv)),
  file.path(out_b2,"summary_SPI_SPEI_skill.csv")
)

cat("✓ B2 complete.\n")
cat("\n========== w10b COMPLETE ==========\n")
cat("  skill_results/ — lag correlations, LOYO skill, figures\n\n")