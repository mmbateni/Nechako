# ==============================================================================
# w11b_PET_bias_nonstationarity.R
# THORNTHWAITE vs. PENMAN-MONTEITH PET BIAS NON-STATIONARITY ANALYSIS
# Nechako River Basin, BC — 2022-2025 Drought Study
#
# PREREQUISITE:
#   Run 3SPEI_ERALand.R with the ADD_TO_3SPEI_ERALand.R block inserted
#   (after line 191). This writes:
#     spei_results_seasonal/WB_basin_average_monthly.csv
#   containing monthly basin-average WB_PM and WB_Thw in mm/month.
#
# UNITS — VERIFIED:
#   In 3SPEI_ERALand.R:
#     precip (ERA5-Land) starts as m/day mean rate
#     × 1000            -> mm/day
#     × days_in_month   -> mm/month  [this is wb_basin]
#   PET CSVs (from 2preq_PET_ERALand.R):
#     mean_pet           -> mm/day
#   Therefore:
#     bias = PET_Thw - PET_PM                   [mm/day]
#     bias_mm_month = bias × days_in_month      [mm/month]
#   Both WB_PM and bias are converted to mm/month BEFORE k-month rolling sums.
#   This ensures all quantities in the correction formula share the same units.
#
# SCALE-SPECIFIC CORRECTION (for SPEI-k, ending calendar month m):
#   WB_k(t)        = k-month rolling sum of WB_PM(t)      [mm/month accumulated]
#   bias_k(t)      = k-month rolling sum of bias_mm_month(t)
#   SD_WB_k(m)     = SD of WB_k over P1 years with ending month m
#   delta_bias_k(m)= mean(bias_k over P3, ending month m)
#                  - mean(bias_k over P1, ending month m)
#   delta_SPEI_k(m)= delta_bias_k(m) / SD_WB_k(m)         [dimensionless]
#   SPEI_Thw_corr  = SPEI_Thw + delta_SPEI_k(m)
#   F_thm_corr     = |SPEI_Thw_corr| / (|SPEI_PM| + eps)
#
# PERIODS:
#   P1  1950-1990  pre-warming reference (41 years)
#   P3  2022-2025  event window (4 years)
#
# INPUTS:
#   spei_results_seasonal/WB_basin_average_monthly.csv   (from 3SPEI_ERALand.R)
#   monthly_data_direct/ERA5Land_Nechako_PET_monthly_summary.csv
#   monthly_data_direct/ERA5Land_Nechako_PET_Thornthwaite_monthly_summary.csv
#   decomp_results/decomp_full_timeseries.csv
#
# OUTPUTS -> decomp_results/bias_nonstationarity/
#   bias_monthly_timeseries.csv          (monthly bias, mm/day and mm/month)
#   correction_params_by_scale.csv       (delta_bias_k, SD_WB_k, delta_SPEI_k)
#   bias_period_seasonal_summary.csv
#   significance_tests.csv
#   SPEI_Thw_corrected_P3.csv
#   bias_manuscript_tables.xlsx
#   bias_nonstationarity.pdf
# ==============================================================================

library(zoo)
library(lubridate)
library(dplyr)
library(tidyr)
library(ggplot2)
library(writexl)

# ── Paths ─────────────────────────────────────────────────────────────────────
setwd("D:/Nechako_Drought/Nechako/")
PET_DIR    <- "monthly_data_direct"
WB_DIR     <- "spei_results_seasonal"          # written by 3SPEI_ERALand.R
DECOMP_DIR <- "decomp_results"
OUT_DIR    <- file.path(DECOMP_DIR, "bias_nonstationarity")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Settings ──────────────────────────────────────────────────────────────────
P1_START <- 1950L; P1_END <- 1990L
P3_START <- 2022L; P3_END <- 2025L
SCALES   <- c(1L, 2L, 3L)
EPS      <- 1e-6

season_of <- function(m) dplyr::case_when(
  m %in% c(12L,1L,2L)  ~ "DJF",
  m %in% c(3L,4L,5L)   ~ "MAM",
  m %in% c(6L,7L,8L)   ~ "JJA",
  m %in% c(9L,10L,11L) ~ "SON"
)

days_in_month_lut <- c(31,28.25,31,30,31,30,31,31,30,31,30,31)  # avg; leap-year aware below

cat("==============================================================\n")
cat("  w11b: PET BIAS NON-STATIONARITY ANALYSIS (P1 vs P3)\n")
cat("==============================================================\n\n")

# ==============================================================================
# SECTION 1: LOAD ALL INPUTS
# ==============================================================================
cat("--- SECTION 1: Loading inputs ---\n")

wb_csv <- file.path(WB_DIR, "WB_basin_average_monthly.csv")
for (f in c(wb_csv,
            file.path(PET_DIR, "ERA5Land_Nechako_PET_monthly_summary.csv"),
            file.path(PET_DIR, "ERA5Land_Nechako_PET_Thornthwaite_monthly_summary.csv"),
            file.path(DECOMP_DIR, "decomp_full_timeseries.csv"))) {
  if (!file.exists(f)) stop("File not found — check prerequisites: ", f)
}

wb_raw  <- read.csv(wb_csv)                  # mm/month, from 3SPEI_ERALand.R
pm_raw  <- read.csv(file.path(PET_DIR, "ERA5Land_Nechako_PET_monthly_summary.csv"))
thw_raw <- read.csv(file.path(PET_DIR, "ERA5Land_Nechako_PET_Thornthwaite_monthly_summary.csv"))
decomp  <- read.csv(file.path(DECOMP_DIR, "decomp_full_timeseries.csv"))

cat(sprintf("  WB CSV:  %d rows, cols: %s\n", nrow(wb_raw),
            paste(names(wb_raw), collapse = ", ")))
cat(sprintf("  PM PET:  %d rows (%s to %s)\n",
            nrow(pm_raw), pm_raw$date[1], tail(pm_raw$date, 1)))
cat(sprintf("  Thw PET: %d rows\n", nrow(thw_raw)))
cat(sprintf("  Decomp:  %d rows, scales: %s\n",
            nrow(decomp), paste(sort(unique(decomp$scale)), collapse = ",")))

# ==============================================================================
# SECTION 2: BUILD MONTHLY BASE TABLE
#   - WB_PM in mm/month (direct from CSV — units already verified)
#   - bias in mm/day (from PET CSVs)
#   - bias_mm_month = bias × actual days_in_month (leap-year aware)
# ==============================================================================
cat("\n--- SECTION 2: Building monthly base table ---\n")

# Merge WB and PET CSVs on date/year/month
base <- inner_join(
  wb_raw  %>% select(date, year, month, days_in_month, wb_pm_mm_month),
  pm_raw  %>% select(date, year, month, mean_pet) %>% rename(pet_pm_mm_day  = mean_pet),
  by = c("date", "year", "month")
) %>%
  inner_join(
    thw_raw %>% select(date, year, month, mean_pet) %>% rename(pet_thw_mm_day = mean_pet),
    by = c("date", "year", "month")
  ) %>%
  arrange(year, month) %>%
  mutate(
    bias_mm_day   = pet_thw_mm_day - pet_pm_mm_day,          # mm/day
    bias_mm_month = bias_mm_day * days_in_month,              # mm/month (same units as WB)
    season        = season_of(month),
    in_P1         = year >= P1_START & year <= P1_END,
    in_P3         = year >= P3_START & year <= P3_END
  )

cat(sprintf("  Merged: %d months\n", nrow(base)))
cat(sprintf("  WB_PM:        %.2f to %.2f mm/month\n",
            min(base$wb_pm_mm_month,  na.rm = TRUE),
            max(base$wb_pm_mm_month,  na.rm = TRUE)))
cat(sprintf("  bias (mm/day):  %.4f to %.4f\n",
            min(base$bias_mm_day,     na.rm = TRUE),
            max(base$bias_mm_day,     na.rm = TRUE)))
cat(sprintf("  bias (mm/month): %.3f to %.3f\n",
            min(base$bias_mm_month,   na.rm = TRUE),
            max(base$bias_mm_month,   na.rm = TRUE)))

write.csv(base %>% select(date, year, month, season, days_in_month,
                           wb_pm_mm_month, pet_pm_mm_day, pet_thw_mm_day,
                           bias_mm_day, bias_mm_month, in_P1, in_P3),
          file.path(OUT_DIR, "bias_monthly_timeseries.csv"), row.names = FALSE)
cat("  Saved: bias_monthly_timeseries.csv\n")

# ==============================================================================
# SECTION 3: SCALE-SPECIFIC CORRECTION PARAMETERS
#
# For each scale k in {1, 2, 3} and each ending calendar month m:
#   WB_k(t)         = zoo::rollsum(WB_PM, k, align="right")   [mm/month ×k]
#   bias_k(t)       = zoo::rollsum(bias_mm_month, k, align="right")
#   SD_WB_k(m)      = SD of WB_k over P1 years with ending month == m
#   mean_bias_k_P1(m)= mean(bias_k) over P1, ending month == m
#   mean_bias_k_P3(m)= mean(bias_k) over P3, ending month == m
#   delta_bias_k(m)  = mean_bias_k_P3(m) - mean_bias_k_P1(m)
#   delta_SPEI_k(m)  = delta_bias_k(m) / SD_WB_k(m)
#
# All quantities are in mm/month (k-month accumulated), so the ratio is
# dimensionless — correct SPEI units. No unit mismatch possible.
# ==============================================================================
cat("\n--- SECTION 3: Scale-specific correction parameters ---\n")

correction_all <- list()

for (k in SCALES) {
  cat(sprintf("\n  Scale k=%d:\n", k))

  # k-month rolling sums (align="right": window ENDS at position i)
  wb_k    <- rollsum(base$wb_pm_mm_month, k, align = "right", fill = NA)
  bias_k  <- rollsum(base$bias_mm_month,  k, align = "right", fill = NA)

  tmp <- base %>%
    mutate(wb_k   = wb_k,
           bias_k = bias_k)

  # P1 statistics per ending calendar month
  p1_stats <- tmp %>%
    filter(in_P1, !is.na(wb_k), !is.na(bias_k)) %>%
    group_by(month) %>%
    summarise(
      n_P1              = n(),
      mean_bias_k_P1    = mean(bias_k, na.rm = TRUE),
      sd_wb_k_P1        = sd(wb_k,    na.rm = TRUE),   # CORRECT denominator
      mean_wb_k_P1      = mean(wb_k,   na.rm = TRUE),
      .groups = "drop"
    )

  # P3 statistics per ending calendar month
  p3_stats <- tmp %>%
    filter(in_P3, !is.na(bias_k)) %>%
    group_by(month) %>%
    summarise(
      n_P3              = n(),
      mean_bias_k_P3    = mean(bias_k, na.rm = TRUE),
      .groups = "drop"
    )

  corr_k <- inner_join(p1_stats, p3_stats, by = "month") %>%
    mutate(
      scale            = k,
      month_name       = month.abb[month],
      season           = season_of(month),
      # Change in k-month accumulated bias [mm/month accumulated]
      delta_bias_k     = mean_bias_k_P3 - mean_bias_k_P1,
      # Convert to SPEI units: divide by SD of k-month accumulated WB_PM
      # Both numerator and denominator in same units (mm/month accumulated)
      delta_SPEI_k     = delta_bias_k / sd_wb_k_P1
    ) %>%
    select(scale, month, month_name, season,
           n_P1, n_P3,
           mean_bias_k_P1, mean_bias_k_P3, delta_bias_k,
           sd_wb_k_P1, mean_wb_k_P1,
           delta_SPEI_k)

  correction_all[[as.character(k)]] <- corr_k

  cat(sprintf("  %-5s %-6s %-14s %-14s %-14s %-14s %-12s\n",
              "Month","Season","bias_k_P1","bias_k_P3",
              "delta_bias_k","SD_WB_k_P1","delta_SPEI_k"))
  for (i in seq_len(nrow(corr_k))) {
    r <- corr_k[i, ]
    cat(sprintf("  %-5s %-6s %+13.3f  %+13.3f  %+13.3f  %13.3f  %+11.4f\n",
                r$month_name, r$season,
                r$mean_bias_k_P1, r$mean_bias_k_P3, r$delta_bias_k,
                r$sd_wb_k_P1,     r$delta_SPEI_k))
  }
}

correction_df <- bind_rows(correction_all)
write.csv(correction_df,
          file.path(OUT_DIR, "correction_params_by_scale.csv"),
          row.names = FALSE)
cat("\n  Saved: correction_params_by_scale.csv\n")

# ==============================================================================
# SECTION 4: SEASONAL AND ANNUAL SUMMARIES (aggregated AFTER Section 3)
# ==============================================================================
cat("\n--- SECTION 4: Seasonal and annual summaries ---\n")

seas_summary <- correction_df %>%
  group_by(scale, season) %>%
  summarise(
    mean_delta_bias_k  = mean(delta_bias_k,  na.rm = TRUE),
    mean_delta_SPEI_k  = mean(delta_SPEI_k,  na.rm = TRUE),
    mean_sd_wb_k_P1    = mean(sd_wb_k_P1,     na.rm = TRUE),
    .groups = "drop"
  )

ann_summary <- correction_df %>%
  group_by(scale) %>%
  summarise(
    mean_delta_bias_k  = mean(delta_bias_k,  na.rm = TRUE),
    mean_delta_SPEI_k  = mean(delta_SPEI_k,  na.rm = TRUE),
    mean_sd_wb_k_P1    = mean(sd_wb_k_P1,     na.rm = TRUE),
    .groups = "drop"
  ) %>% mutate(season = "Annual")

cat("\n  Annual delta_SPEI per scale:\n")
print(ann_summary %>% select(scale, mean_delta_bias_k, mean_sd_wb_k_P1, mean_delta_SPEI_k))

write.csv(bind_rows(ann_summary, seas_summary),
          file.path(OUT_DIR, "bias_period_seasonal_summary.csv"),
          row.names = FALSE)
cat("  Saved: bias_period_seasonal_summary.csv\n")

# ==============================================================================
# SECTION 5: SIGNIFICANCE TESTS (annual-aggregated to reduce autocorrelation)
# ==============================================================================
cat("\n--- SECTION 5: Significance tests ---\n")

sig_list <- list()

# Use scale=1 monthly bias in mm/month for annual-aggregated test
# (scale=1: bias_k = bias_mm_month; most straightforward interpretation)
p1_ann <- base %>%
  filter(in_P1) %>%
  group_by(year) %>%
  summarise(bias_ann = mean(bias_mm_month, na.rm = TRUE), .groups = "drop") %>%
  pull(bias_ann)

p3_ann <- base %>%
  filter(in_P3) %>%
  group_by(year) %>%
  summarise(bias_ann = mean(bias_mm_month, na.rm = TRUE), .groups = "drop") %>%
  pull(bias_ann)

tt_ann <- t.test(p3_ann, p1_ann, alternative = "two.sided", var.equal = FALSE)
sig_list[["Annual_k1"]] <- data.frame(
  comparison    = "P3 vs P1 — Annual, bias in mm/month (k=1 monthly bias)",
  n_P1          = length(p1_ann),
  n_P3          = length(p3_ann),
  mean_bias_P1  = round(mean(p1_ann), 4),
  mean_bias_P3  = round(mean(p3_ann), 4),
  delta_bias    = round(mean(p3_ann) - mean(p1_ann), 4),
  t_stat        = round(tt_ann$statistic, 3),
  p_value       = round(tt_ann$p.value, 4),
  note          = "n_P3=4 years; p-value indicative only"
)

# Seasonal tests
for (s in c("DJF","MAM","JJA","SON")) {
  v1 <- base %>% filter(in_P1, season == s) %>% pull(bias_mm_month)
  v3 <- base %>% filter(in_P3, season == s) %>% pull(bias_mm_month)
  if (length(v3) < 3) next
  tt <- t.test(v3, v1, alternative = "two.sided", var.equal = FALSE)
  wt <- wilcox.test(v3, v1, alternative = "two.sided", exact = FALSE)
  sig_list[[s]] <- data.frame(
    comparison   = paste0("P3 vs P1 — ", s, " (monthly obs, mm/month)"),
    n_P1 = length(v1), n_P3 = length(v3),
    mean_bias_P1 = round(mean(v1), 4),
    mean_bias_P3 = round(mean(v3), 4),
    delta_bias   = round(mean(v3) - mean(v1), 4),
    t_stat       = round(tt$statistic, 3),
    p_value      = round(tt$p.value, 4),
    note         = "monthly autocorrelation present"
  )
}

sig_df <- bind_rows(sig_list)
cat("\n  Results:\n"); print(sig_df %>% select(-note))
write.csv(sig_df, file.path(OUT_DIR, "significance_tests.csv"), row.names = FALSE)
cat("  Saved: significance_tests.csv\n")

# ==============================================================================
# SECTION 6: APPLY SCALE-SPECIFIC CORRECTION TO SPEI_Thw (P3 only)
#
# For each P3 month, look up delta_SPEI_k(month) from the scale-k correction
# table and add it to SPEI_Thw. Recompute f_thm from corrected value.
# ==============================================================================
cat("\n--- SECTION 6: Applying scale-specific correction to SPEI_Thw (P3) ---\n")

corrected_list <- list()

for (k in SCALES) {
  corr_k  <- correction_all[[as.character(k)]]
  delta_lk <- setNames(corr_k$delta_SPEI_k, as.character(corr_k$month))

  d_sc <- decomp %>%
    filter(scale == k, year >= P3_START, year <= P3_END) %>%
    mutate(
      delta_SPEI_k   = delta_lk[as.character(month)],
      SPEI_Thw_corr  = SPEI_Thw + delta_SPEI_k,
      f_thm_orig     = abs(SPEI_Thw)      / (abs(SPEI_PM) + EPS),
      f_thm_corr     = abs(SPEI_Thw_corr) / (abs(SPEI_PM) + EPS),
      delta_f_pp     = 100 * (f_thm_corr - f_thm_orig),
      season         = season_of(month)
    ) %>%
    select(date, year, month, season, scale,
           SPEI_PM, SPEI_Thw, delta_SPEI_k,
           SPEI_Thw_corr, f_thm_orig, f_thm_corr, delta_f_pp)

  corrected_list[[as.character(k)]] <- d_sc

  seas_sc <- d_sc %>%
    group_by(season) %>%
    summarise(
      f_orig_pct  = round(100 * mean(f_thm_orig, na.rm = TRUE), 2),
      f_corr_pct  = round(100 * mean(f_thm_corr, na.rm = TRUE), 2),
      delta_f_pp  = round(mean(delta_f_pp, na.rm = TRUE), 2),
      .groups = "drop"
    )

  ann_sc <- d_sc %>%
    summarise(
      f_orig_pct = round(100 * mean(f_thm_orig, na.rm = TRUE), 2),
      f_corr_pct = round(100 * mean(f_thm_corr, na.rm = TRUE), 2),
      delta_f_pp = round(mean(delta_f_pp, na.rm = TRUE), 2)
    )

  cat(sprintf("\n  SPEI-%d (48-month mean): f_thm %.2f%% → %.2f%% (Δ = %+.2f pp)\n",
              k, ann_sc$f_orig_pct, ann_sc$f_corr_pct, ann_sc$delta_f_pp))
  cat("    Seasonal breakdown:\n")
  print(seas_sc)
}

corrected_all_df <- bind_rows(corrected_list)
write.csv(corrected_all_df,
          file.path(OUT_DIR, "SPEI_Thw_corrected_P3.csv"),
          row.names = FALSE)
cat("\n  Saved: SPEI_Thw_corrected_P3.csv\n")

# ==============================================================================
# SECTION 7: FIGURES
# ==============================================================================
cat("\n--- SECTION 7: Figures ---\n")

pdf(file.path(OUT_DIR, "bias_nonstationarity.pdf"), width = 14, height = 9)

# Fig 1: Full monthly bias time series (mm/month, scale=1)
base_plt <- base %>%
  mutate(date_fmt  = as.Date(paste0(date, "-01")),
         bias_12ma = stats::filter(bias_mm_month, rep(1/12,12), sides=2) %>% as.numeric())
p1_m <- mean(base$bias_mm_month[base$in_P1], na.rm = TRUE)
p3_m <- mean(base$bias_mm_month[base$in_P3], na.rm = TRUE)

f1 <- ggplot(base_plt, aes(x = date_fmt, y = bias_mm_month)) +
  annotate("rect", xmin=as.Date("1950-01-01"), xmax=as.Date("1990-12-31"),
           ymin=-Inf, ymax=Inf, alpha=0.07, fill="#4393c3") +
  annotate("rect", xmin=as.Date("2022-01-01"), xmax=as.Date("2025-12-31"),
           ymin=-Inf, ymax=Inf, alpha=0.15, fill="#762a83") +
  geom_line(colour="grey70", linewidth=0.35, alpha=0.8) +
  geom_line(aes(y=bias_12ma), colour="#e31a1c", linewidth=1.0, na.rm=TRUE) +
  geom_hline(yintercept=0,   linetype="dashed", colour="black") +
  geom_hline(yintercept=p1_m, colour="#4393c3", linewidth=0.9) +
  geom_hline(yintercept=p3_m, colour="#762a83", linewidth=0.9) +
  annotate("text", x=as.Date("1970-01-01"), y=p1_m+0.5,
           label=sprintf("P1 mean = %+.3f mm/month", p1_m),
           colour="#4393c3", size=3.0) +
  annotate("text", x=as.Date("2023-07-01"), y=p3_m+0.5,
           label=sprintf("P3 mean = %+.3f mm/month", p3_m),
           colour="#762a83", size=3.0) +
  scale_x_date(date_breaks="10 years", date_labels="%Y") +
  labs(title="Monthly bias: (PET_Thw - PET_PM) × days_in_month [mm/month]",
       subtitle="Blue = P1 (1950-1990); Purple = P3 (2022-2025); Red = 12-month MA",
       x="Date", y="Bias (mm/month)") +
  theme_bw(base_size=11)
print(f1)

# Fig 2: delta_SPEI_k by month, faceted by scale
f2 <- ggplot(correction_df,
             aes(x=factor(month_name, levels=month.abb),
                 y=delta_SPEI_k, fill=season)) +
  geom_col(alpha=0.85) +
  geom_hline(yintercept=0, linewidth=0.5) +
  facet_wrap(~paste0("SPEI-", scale), ncol=3) +
  scale_fill_manual(values=c(DJF="#4393c3",MAM="#74c476",JJA="#fd8d3c",SON="#9e9ac8"),
                    name="Season") +
  labs(title="Scale-specific delta_SPEI per calendar month (P3 minus P1)",
       subtitle=paste0("delta_SPEI_k(m) = delta_bias_k(m) / SD_WB_k_P1(m). ",
                       "Units: mm/month (k-month accumulated) / mm/month = dimensionless."),
       x="Month", y="delta_SPEI (SPEI units)") +
  theme_bw(base_size=11)
print(f2)

# Fig 3: SD_WB_k(m) by scale — shows why scale matters for denominator
sd_plot <- correction_df %>%
  mutate(month_name = factor(month_name, levels=month.abb),
         scale_lbl  = paste0("SPEI-", scale, " (k=", scale, " month SD_WB)"))

f3 <- ggplot(sd_plot, aes(x=month_name, y=sd_wb_k_P1,
                            colour=scale_lbl, group=scale_lbl)) +
  geom_line(linewidth=1.1) +
  geom_point(size=2.2) +
  labs(title="SD of k-month accumulated WB_PM per calendar month (P1: 1950-1990)",
       subtitle=paste0("SD grows with k because it accumulates variance over more months. ",
                       "Using k=1 SD for k=3 SPEI would severely understate the correction."),
       x="Ending calendar month", y="SD_WB_k (mm/month accumulated)",
       colour=NULL) +
  theme_bw(base_size=11) + theme(legend.position="bottom")
print(f3)

# Fig 4: Corrected vs original f_thm by season and scale
corr_seas_plot <- bind_rows(lapply(SCALES, function(k) {
  d <- corrected_list[[as.character(k)]]
  bind_rows(
    d %>% group_by(season) %>%
      summarise(f_pct=100*mean(f_thm_orig,na.rm=TRUE),.groups="drop") %>%
      mutate(type="Original"),
    d %>% group_by(season) %>%
      summarise(f_pct=100*mean(f_thm_corr,na.rm=TRUE),.groups="drop") %>%
      mutate(type="Bias-corrected")
  ) %>% mutate(scale_lbl=paste0("SPEI-",k))
})) %>%
  mutate(type=factor(type,levels=c("Original","Bias-corrected")),
         season=factor(season,levels=c("DJF","MAM","JJA","SON")))

f4 <- ggplot(corr_seas_plot, aes(x=season, y=f_pct, fill=type)) +
  geom_col(position=position_dodge(0.7), width=0.6, alpha=0.85) +
  geom_hline(yintercept=100, linetype="dashed") +
  facet_wrap(~scale_lbl, ncol=3) +
  scale_fill_manual(values=c("Original"="#d6604d","Bias-corrected"="#4393c3"),
                    name=NULL) +
  labs(title="Thermodynamic fraction: original vs bias-corrected (P3: 2022-2025)",
       subtitle="Each scale uses its own k-month accumulated delta_bias and SD_WB",
       x="Season", y="Thermodynamic fraction (%)") +
  theme_bw(base_size=11) + theme(legend.position="bottom")
print(f4)

dev.off()
cat(sprintf("  Saved: %s\n", file.path(OUT_DIR, "bias_nonstationarity.pdf")))

# ==============================================================================
# SECTION 8: EXCEL WORKBOOK
# ==============================================================================
cat("\n--- SECTION 8: Excel workbook ---\n")
write_xlsx(
  list(
    "Correction_params_by_scale" = correction_df,
    "Seasonal_annual_summary"    = bind_rows(ann_summary, seas_summary),
    "Significance_tests"         = sig_df,
    "SPEI_Thw_corrected_P3"      = corrected_all_df,
    "Base_monthly_timeseries"    = base %>%
      select(date, year, month, season, days_in_month,
             wb_pm_mm_month, bias_mm_day, bias_mm_month, in_P1, in_P3)
  ),
  path = file.path(OUT_DIR, "bias_manuscript_tables.xlsx")
)
cat("  Saved: bias_manuscript_tables.xlsx\n")

# ==============================================================================
# SECTION 9: MANUSCRIPT SUMMARY
# ==============================================================================
cat("\n")
cat("==============================================================\n")
cat("  MANUSCRIPT REPORTING SUMMARY\n")
cat("==============================================================\n")
cat(sprintf("
  Bias (PET_Thw - PET_PM) converted to mm/month (× days_in_month):
    P1 annual mean (1950-1990):  %+.3f mm/month
    P3 annual mean (2022-2025):  %+.3f mm/month
    Annual delta bias:           %+.3f mm/month
    Welch t-test (annual agg.):  t=%.3f, p=%.4f (n_P3=%d — indicative)\n\n",
  mean(p1_ann, na.rm=TRUE), mean(p3_ann, na.rm=TRUE),
  mean(p3_ann, na.rm=TRUE) - mean(p1_ann, na.rm=TRUE),
  tt_ann$statistic, tt_ann$p.value, length(p3_ann)
))

for (k in SCALES) {
  d <- corrected_list[[as.character(k)]]
  fo <- round(100 * mean(d$f_thm_orig, na.rm=TRUE), 2)
  fc <- round(100 * mean(d$f_thm_corr, na.rm=TRUE), 2)
  df <- round(fc - fo, 2)
  jja <- d %>% filter(season=="JJA") %>%
    summarise(o=100*mean(f_thm_orig,na.rm=TRUE), c=100*mean(f_thm_corr,na.rm=TRUE))
  cat(sprintf("  SPEI-%d: f_thm %.2f%% → %.2f%% (Δ = %+.2f pp); JJA: %.2f%% → %.2f%%\n",
              k, fo, fc, df, jja$o, jja$c))
}

cat(sprintf("
  Interpretation guide:
    |Δ f_thm| < 1 pp  → bias non-stationarity negligible; report original.
    |Δ f_thm| >= 1 pp → report both values and note that warming-amplified
                         Thornthwaite divergence from PM contributes to the
                         thermodynamic fraction at that time scale.\n"))

cat("==============================================================\n")
cat("  w11b COMPLETE. Outputs: decomp_results/bias_nonstationarity/\n")
cat("==============================================================\n")
