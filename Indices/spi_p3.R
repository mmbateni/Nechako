
# spi_p3.R
# Compute Standardized Precipitation Index (SPI) using Pearson Type III (PE3)
# with L-moment parameter estimation, per Guttman (1999).

# ---- Dependencies ----
# install.packages(c("zoo", "lmomco"))
library(zoo)
library(lmomco) # pelpe3 (fit) and cdfpe3 (CDF)

accumulate_series <- function(x, scale) {
  rollapply(x, width = scale, by = 1, align = "right",
            FUN = sum, na.rm = FALSE, partial = FALSE)
}

month_factor <- function(dates) as.integer(format(dates, "%m"))

lmoms_safe <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 10) return(NULL)
  try(lmoments(x), silent = TRUE)
}

spi_p3 <- function(prcp_monthly, dates,
                   scales = c(1, 2, 3, 6, 9, 12, 24, 36, 48, 60, 72),
                   baseline_start = NULL, baseline_end = NULL,
                   min_samples_per_month = 30) {

  stopifnot(length(prcp_monthly) == length(dates))
  z <- zoo(prcp_monthly, order.by = as.Date(dates))

  if (!is.null(baseline_start) && !is.null(baseline_end)) {
    baseline_idx <- index(z) >= as.Date(baseline_start) & index(z) <= as.Date(baseline_end)
  } else {
    baseline_idx <- rep(TRUE, length(z))
  }

  out <- list()

  for (k in scales) {
    ak <- accumulate_series(z, k)
    idx <- index(ak)
    mon <- month_factor(idx)

    spi_vals <- rep(NA_real_, length(ak))

    for (m in 1:12) {
      sel_all  <- which(mon == m)
      if (!is.null(baseline_start) && !is.null(baseline_end)) {
        sel_base <- sel_all[idx[sel_all] >= as.Date(baseline_start) & idx[sel_all] <= as.Date(baseline_end)]
      } else {
        sel_base <- sel_all
      }

      x_base <- coredata(ak)[sel_base]
      x_obs  <- coredata(ak)[sel_all]

      # Enough baseline samples per month
      if (sum(!is.na(x_base)) < min_samples_per_month) next

      # Fraction of zeros (q)
      n_base <- sum(!is.na(x_base))
      m_zero <- sum(x_base == 0, na.rm = TRUE)
      q_zero <- m_zero / n_base

      # PE3 fit on positive baseline values
      x_pos <- x_base[x_base > 0 & !is.na(x_base)]
      if (length(x_pos) < 10) next

      lm <- lmoms_safe(x_pos)
      if (is.null(lm) || inherits(lm, "try-error")) next

      para <- pelpe3(lm) # PE3 parameters (XI, ALPHA, KAPPA)

      # Mixed CDF: H(x) = q + (1-q) * F_PE3(x)
      Fx <- rep(NA_real_, length(x_obs))
      pos_mask <- !is.na(x_obs)
      x_clipped <- pmax(x_obs[pos_mask], 0)
      Fx[pos_mask] <- cdfpe3(q = x_clipped, para = para)
      Hx <- q_zero + (1 - q_zero) * Fx

      # SPI via inverse-normal
      spi_vals[sel_all] <- qnorm(Hx)
    }

    out[[paste0("SPI_", k, "mo_P3")]] <- zoo(spi_vals, order.by = idx)
  }

  return(out)
}

# ---- Example ----
# df <- data.frame(date = seq(as.Date("1980-01-01"), as.Date("2024-12-01"), by = "month"),
#                  prcp = pmax(rnorm(541, mean = 60, sd = 25), 0))
# res_p3 <- spi_p3(prcp_monthly = df$prcp, dates = df$date,
#                  baseline_start = "1991-01-01", baseline_end = "2020-12-01")
# plot(res_p3$SPI_36mo_P3, main = "SPI-36 (Pearson Type III)", ylab = "SPI", xlab = "Time")
``
