# ============================================================================
# NECHAKO BASIN — TOTAL STORAGE DEFICIT INDEX (TSDI)
#   + PIT normalization to standard Normal (TSDI_norm / Standardized TSDI)
#
# Based on:
#   Awange et al. (2016) Adv. Water Resources 94, 45–59
#   Yirdaw et al. (2008) J. Hydrology 356, 84–92
#   Farahmand & AghaKouchak (2015) Adv. Water Resources 76, 140–145
#     (nonparametric Blom plotting-position PIT — same principle as SPI)
#
# INPUT:
#   BEST_ESTIMATE_basin_total_Mm3.csv
#     Required columns:  date, total_storage_Mm3
#
# OUTPUT:
#   Nechako_TSDI_Output.csv   — full time series with TSDI and TSDI_norm
#   Four diagnostic plots
# ============================================================================

rm(list = ls())

# ============================================================================
# PACKAGES
# ============================================================================

packages <- c("tidyverse", "lubridate", "zoo")

for (p in packages) {
  if (!require(p, character.only = TRUE)) {
    install.packages(p)
    library(p, character.only = TRUE)
  }
}

# ============================================================================
# USER SETTINGS
# ============================================================================

INPUT_FILE <- "BEST_ESTIMATE_basin_total_Mm3.csv"

# Palmer drought severity class used to anchor p and q
# –1 mild | –2 moderate | –3 severe | –4 extreme
DROUGHT_CLASS_C <- -3

# Optional rolling smoothing of raw storage before computing TSD
USE_SMOOTHING  <- FALSE
SMOOTH_WINDOW  <- 3   # months

# Blom plotting-position constant (standard for nonparametric SPI-style PIT)
BLOM_A <- 0.44

# ============================================================================
# LOAD DATA
# ============================================================================

df <- read_csv(INPUT_FILE, show_col_types = FALSE)

cat("Columns found:", paste(names(df), collapse = ", "), "\n")

stopifnot(
  "date"              %in% names(df),
  "total_storage_Mm3" %in% names(df)
)

# ============================================================================
# DATE HANDLING
# ============================================================================

df <- df %>%
  mutate(
    date  = as.Date(date),
    year  = year(date),
    month = month(date)
  ) %>%
  arrange(date)

# ============================================================================
# OPTIONAL SMOOTHING
# ============================================================================

if (USE_SMOOTHING) {
  df <- df %>%
    mutate(
      storage_Mm3 = zoo::rollmean(
        total_storage_Mm3,
        k     = SMOOTH_WINDOW,
        fill  = NA,
        align = "center"
      )
    )
} else {
  df <- df %>%
    mutate(storage_Mm3 = total_storage_Mm3)
}

# ============================================================================
# MONTHLY CLIMATOLOGY  (mean / max / min per calendar month)
# Follows Awange et al. (2016) Eq. 1
# ============================================================================

monthly_stats <- df %>%
  group_by(month) %>%
  summarise(
    clim_mean = mean(storage_Mm3, na.rm = TRUE),
    clim_max  = max( storage_Mm3, na.rm = TRUE),
    clim_min  = min( storage_Mm3, na.rm = TRUE),
    .groups   = "drop"
  )

df <- left_join(df, monthly_stats, by = "month")

# ============================================================================
# TOTAL STORAGE DEFICIT  (TSD %)
# Eq. 1 of Awange et al. (2016)
# ============================================================================

df <- df %>%
  mutate(
    TSD = (storage_Mm3 - clim_mean) /
      (clim_max    - clim_min  ) * 100
  )

# ============================================================================
# CUMULATIVE TSD  (used to fit the drought monograph)
# ============================================================================

df <- df %>%
  mutate(cumulative_TSD = cumsum(replace_na(TSD, 0)))

# ============================================================================
# IDENTIFY THE DOMINANT (LONGEST) CONTINUOUSLY NEGATIVE TSD PERIOD
# ============================================================================

is_dry <- df$TSD < 0
r      <- rle(is_dry)

ends   <- cumsum(r$lengths)
starts <- ends - r$lengths + 1

dry_runs <- which(r$values == TRUE)

if (length(dry_runs) == 0) {
  stop("No drought periods detected (no consecutive negative TSD months).")
}

longest_run <- dry_runs[which.max(r$lengths[dry_runs])]
dry_start   <- starts[longest_run]
dry_end     <- ends[longest_run]

drought_df <- df[dry_start:dry_end, ] %>%
  mutate(drought_time = seq_len(n()))

cat(
  "\nDominant drought episode:\n",
  format(min(drought_df$date)), "to", format(max(drought_df$date)),
  "(", nrow(drought_df), "months )\n"
)

# ============================================================================
# LINEAR FIT TO CUMULATIVE TSD DURING DROUGHT  (the "drought monograph")
# Follows Fig. 2 of Awange et al. (2016)
# ============================================================================

fit <- lm(cumulative_TSD ~ drought_time, data = drought_df)
m   <- coef(fit)[2]   # slope
b   <- coef(fit)[1]   # intercept

cat(
  "\nDrought monograph regression:\n",
  "  slope (m) =", round(m, 4), "\n",
  "  intercept (b) =", round(b, 4), "\n"
)

# Guard against degenerate regression (m + b ≈ 0)
denom <- m + b
if (abs(denom) < 1e-6) {
  stop(
    "m + b is near zero (", round(denom, 6), ").",
    " Cannot compute p and q. ",
    "Check that the dominant dry period has a clearly negative cumulative trend."
  )
}

# ============================================================================
# COMPUTE p AND q   (Eq. 3 of Awange et al. 2016)
# ============================================================================

p <- 1 - (m / denom)
q <- -DROUGHT_CLASS_C / denom   # C is negative, so -C is positive → q > 0

cat(
  "\nDerived coefficients:\n",
  "  p =", round(p, 6), "\n",
  "  q =", round(q, 6), "\n"
)

# Soft warning if p is outside the unit interval (memory decay factor should
# ideally be 0 < p < 1 for stability; values slightly outside are possible)
if (p <= 0 || p >= 1) {
  warning(
    "p = ", round(p, 4), " is outside (0, 1). ",
    "Check the regression period and drought class setting."
  )
}

# ============================================================================
# RECURSIVE TSDI   (Eq. 2 of Awange et al. 2016)
# ============================================================================

n    <- nrow(df)
TSDI <- rep(NA_real_, n)

TSDI[1] <- 0.02 * df$TSD[1]   # initial condition per Yirdaw et al. (2008)

for (i in 2:n) {
  if (!is.na(df$TSD[i])) {
    TSDI[i] <- p * TSDI[i - 1] + q * df$TSD[i]
  }
  # if TSD is NA, TSDI stays NA (propagation stops)
}

df$TSDI <- TSDI

# ============================================================================
# PIT NORMALIZATION  →  TSDI_norm  (approximately N(0,1))
#
# Principle — identical to nonparametric SPI (Farahmand & AghaKouchak 2015):
#
#   Step 1 — empirical CDF via Blom plotting position:
#               u_i = (rank_i − 0.44) / (n + 0.12)
#             This maps the ranked TSDI values to U[0,1].
#
#   Step 2 — inverse Normal transform:
#               z_i = Φ⁻¹(u_i)
#             z_i is approximately N(0,1).
#
# Note: this is a univariate transform, so the simple PIT (without the
# Kendall distribution correction) is theoretically correct here — the
# Kendall device is only needed when aggregating multiple variables through
# their joint CDF (see De Michele et al. 2026, Water Resources Research).
# ============================================================================

tsdi_valid <- df %>%
  filter(!is.na(TSDI)) %>%
  mutate(
    rank_tsdi  = rank(TSDI, ties.method = "average"),
    n_valid    = n(),
    u_blom     = (rank_tsdi - BLOM_A) / (n_valid + 1 - 2 * BLOM_A),
    # equivalent to (rank - 0.44) / (n + 0.12)
    TSDI_norm  = qnorm(u_blom)
  ) %>%
  select(date, rank_tsdi, n_valid, u_blom, TSDI_norm)

df <- left_join(df, tsdi_valid, by = "date")

# ============================================================================
# SAVE OUTPUT
# ============================================================================

write_csv(df, "Nechako_TSDI_Output.csv")
cat("\nOutput written: Nechako_TSDI_Output.csv\n")

# ============================================================================
# PLOTS
# ============================================================================

# --- 1. Raw storage -----------------------------------------------------------
p_storage <- ggplot(df, aes(date, storage_Mm3)) +
  geom_line(linewidth = 0.8) +
  labs(
    title    = "Nechako Basin Total Lake Storage",
    subtitle = "GloLakes basin-total storage",
    x        = NULL,
    y        = "Storage (Mm³)"
  ) +
  theme_bw()

print(p_storage)

# --- 2. TSD (%) ---------------------------------------------------------------
p_tsd <- ggplot(df, aes(date, TSD)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_line(linewidth = 0.8) +
  labs(
    title    = "Total Storage Deficit (TSD %)",
    subtitle = "Nechako Basin",
    x        = NULL,
    y        = "TSD (%)"
  ) +
  theme_bw()

print(p_tsd)

# --- 3. Raw TSDI (Palmer-style thresholds) ------------------------------------
p_tsdi <- ggplot(df, aes(date, TSDI)) +
  
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  geom_line(linewidth = 0.8) +
  
  labs(
    title    = "Nechako Basin TSDI",
    subtitle = paste0("p = ", round(p, 4), "  |  q = ", round(q, 4),
                      "  |  C = ", DROUGHT_CLASS_C),
    x        = NULL,
    y        = "TSDI"
  ) +
  theme_bw()

print(p_tsdi)

# --- 4. TSDI_norm (standard-Normal / SPI-style thresholds) -------------------
p_tsdi_norm <- ggplot(df, aes(date, TSDI_norm)) +
  
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  geom_line(linewidth = 0.8) +
  
  labs(
    title    = "Nechako Basin Standardized TSDI  (TSDI_norm)",
    subtitle = "Nonparametric PIT via Blom plotting position → N(0,1)",
    x        = NULL,
    y        = "TSDI_norm  [standard normal units]"
  ) +
  theme_bw()

print(p_tsdi_norm)

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n============================================\n")
cat("TSDI ANALYSIS COMPLETE\n")
cat("============================================\n")

cat("\nTime range:\n")
cat(" ", format(min(df$date)), "to", format(max(df$date)), "\n")

cat("\nRaw TSDI range:\n")
cat(" ", round(min(df$TSDI, na.rm = TRUE), 3), "to",
    round(max(df$TSDI, na.rm = TRUE), 3), "\n")

cat("\nStandardized TSDI_norm range:\n")
cat(" ", round(min(df$TSDI_norm, na.rm = TRUE), 3), "to",
    round(max(df$TSDI_norm, na.rm = TRUE), 3), "\n")

cat("\nColumns in output CSV:\n")
cat(" ", paste(names(df), collapse = ", "), "\n")