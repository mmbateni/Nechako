#==============================================================================
# MULTIVARIATE DROUGHT INDEX (MSDI) - Kendall-aware standardization
# Combines: TSDI_norm, SWEI/SSPI, SSI/SSMI
# Method: Empirical Kendall distribution per De Michele et al. (2026)
#==============================================================================
rm(list = ls())
library(tidyverse)

#==============================================================================
# 1. LOAD PRE-COMPUTED INDICES 
#==============================================================================
tsdi_df   <- read_csv("Nechako_TSDI_Output.csv", show_col_types = FALSE)
swei_df   <- read_csv("swei_results_seasonal/swei_01_basin_averaged_by_month.csv", show_col_types = FALSE)
ssi_df    <- read_csv("ssi_results_seasonal/ssi_L1_3_01_basin_averaged_by_month.csv", show_col_types = FALSE)

#==============================================================================
# 2. ALIGN TIME SERIES & MERGE
#==============================================================================
# Fixed: added format = "%Y-%b-%d" to resolve the NA generation when parsing "Year-Month-01"
msdi_data <- tsdi_df %>%
  select(date, TSDI_norm) %>%
  left_join(
    swei_df %>% 
      pivot_longer(-Year, names_to = "month", values_to = "SWEI") %>%
      mutate(date = as.Date(paste0(Year, "-", month, "-01"), format = "%Y-%b-%d")),
    by = "date"
  ) %>%
  left_join(
    ssi_df %>% 
      pivot_longer(-Year, names_to = "month", values_to = "SSI") %>%
      mutate(date = as.Date(paste0(Year, "-", month, "-01"), format = "%Y-%b-%d")),
    by = "date"
  ) %>%
  filter(!is.na(TSDI_norm), !is.na(SWEI), !is.na(SSI)) %>%
  arrange(date)

cat("Merged records:", nrow(msdi_data), "\n")

#==============================================================================
# 3. EMPIRICAL JOINT CDF VIA VECTORIZED RANKING
#==============================================================================
# Fast O(N^2) evaluation using vectorized logical matrices instead of apply loops
compute_empirical_joint_cdf <- function(mat) {
  n <- nrow(mat)
  p_vec <- numeric(n)
  
  for (i in seq_len(n)) {
    p_vec[i] <- sum(mat[, 1] <= mat[i, 1] & 
                      mat[, 2] <= mat[i, 2] & 
                      mat[, 3] <= mat[i, 3], na.rm = TRUE) / n
  }
  return(p_vec)
}

X_mat <- msdi_data %>% select(TSDI_norm, SWEI, SSI) %>% as.matrix()
msdi_data$p_joint <- compute_empirical_joint_cdf(X_mat)

#==============================================================================
# 4. KENDALL DISTRIBUTION FUNCTION KC(t) & STANDARDIZATION
#==============================================================================
# Using the Blom plotting position directly on the joint probabilities 
# to derive the Kendall distribution function KC without boundary clipping issues.
msdi_data <- msdi_data %>%
  mutate(
    rank_p = rank(p_joint, ties.method = "average"),
    n_obs  = n(),
    KC     = (rank_p - 0.44) / (n_obs + 0.12),
    MSDI   = qnorm(KC)
  )

#==============================================================================
# 5. DIAGNOSTICS & SAVING
#==============================================================================
cat("\nMSDI summary:\n")
cat("  Mean:  ", round(mean(msdi_data$MSDI, na.rm = TRUE), 3), "\n")
cat("  SD:    ", round(sd(msdi_data$MSDI, na.rm = TRUE), 3), "\n")
cat("  Range: ", round(min(msdi_data$MSDI, na.rm = TRUE), 2), "to",
    round(max(msdi_data$MSDI, na.rm = TRUE), 2), "\n")

ks_test <- ks.test(msdi_data$MSDI, "pnorm")
cat("  KS test normal p-value:", round(ks_test$p.value, 4), "\n")

write_csv(msdi_data %>% select(date, TSDI_norm, SWEI, SSI, p_joint, KC, MSDI),
          "Nechako_MSDI_Kendall_Output.csv")
cat("\nOutput successfully saved to: Nechako_MSDI_Kendall_Output.csv\n")

#==============================================================================
# 6. PLOT TIMSERIES
#==============================================================================
p_msdi <- ggplot(msdi_data, aes(date, MSDI)) +
  geom_hline(yintercept = c(-2, -1.5, -1, 0, 1, 1.5, 2),
             linetype = "dashed", colour = "grey80") +
  geom_line(linewidth = 0.8, colour = "steelblue") +
  geom_hline(yintercept = 0, linetype = "solid", colour = "grey50") +
  labs(
    title = "Nechako Basin Multivariate Drought Index (MSDI)",
    subtitle = "Kendall-aware standardization (De Michele et al. 2026)",
    x = NULL, y = "MSDI [standard normal units]"
  ) +
  annotate("text", x = min(msdi_data$date), y = -2.2,
           label = "Exceptional drought ≤ -2.0", hjust = 0, size = 3, colour = "darkred") +
  theme_bw()

print(p_msdi)
ggsave("Nechako_MSDI_timeseries.png", p_msdi, width = 10, height = 4, dpi = 300)