#============================================================================
# METHODOLOGY GUIDE: SPI/SPEI DROUGHT SEVERITY-AREA-FREQUENCY ANALYSIS
#============================================================================
library(terra); library(ggplot2); library(dplyr); library(copula); library(gridExtra)
DROUGHT_THRESHOLD <- -0.5

analyze_multiple_thresholds <- function(index_rast, cell_area, dates, thresholds=c(-0.5, -1.0, -1.5), index_name="SPI") {
  results <- list()
  for(thresh in thresholds) {
    n_time <- terra::nlyr(index_rast); sev <- numeric(n_time); area <- numeric(n_time)
    for(t in 1:n_time) {
      v <- terra::values(index_rast[[t]]); valid <- is.finite(v); dry <- valid & v < thresh
      if(sum(dry)>0 & sum(valid)>0) { area[t] <- 100*sum(dry)/sum(valid); sev[t] <- mean(thresh - v[dry]) }
    }
    results[[as.character(thresh)]] <- data.frame(date=dates, severity=sev, area_pct=area)
  }
  df <- dplyr::bind_rows(lapply(names(results), function(th) data.frame(date=results[[th]]$date, area_pct=results[[th]]$area_pct, threshold=th)))
  p <- ggplot(df, aes(x=as.Date(date), y=area_pct, colour=factor(threshold))) + geom_line() + theme_bw()
  pdf(sprintf("%s_threshold_sensitivity.pdf", index_name), width=12, height=5); print(p); dev.off()
  invisible(results)
}

seasonal_drought_analysis <- function(drought_chars, index_name, output_dir=".") {
  stats <- drought_chars %>% group_by(month) %>% summarise(mean_severity=mean(severity[severity>0],na.rm=TRUE), mean_area=mean(area_pct[area_pct>0],na.rm=TRUE), .groups="drop")
  p1 <- ggplot(stats, aes(x=month, y=mean_severity)) + geom_col(fill="#f46d43") + scale_x_continuous(breaks=1:12, labels=month.abb) + theme_bw()
  p2 <- ggplot(stats, aes(x=month, y=mean_area)) + geom_col(fill="#4575b4") + scale_x_continuous(breaks=1:12, labels=month.abb) + theme_bw()
  pdf(file.path(output_dir, sprintf("%s_seasonal.pdf", index_name)), width=12, height=6); grid.arrange(p1, p2, ncol=2); dev.off()
  stats
}

map_drought_frequency <- function(index_rast, threshold=DROUGHT_THRESHOLD, output_dir="drought_maps", index_name="SPI") {
  dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
  freq <- terra::app(index_rast, function(x) 100*mean(x<threshold, na.rm=TRUE))
  sev  <- terra::app(index_rast, function(x) { d <- x[is.finite(x) & x<threshold]; if(length(d)>0) mean(threshold-d) else 0 })
  terra::writeCDF(freq, file.path(output_dir, sprintf("%s_freq.nc", index_name))); terra::writeCDF(sev, file.path(output_dir, sprintf("%s_sev.nc", index_name)))
  invisible(list(frequency=freq, severity=sev))
}

compare_return_period_methods <- function(saf_cond, saf_kend, index_name, output_file=NULL) {
  if(is.null(saf_cond)||is.null(saf_kend)) return(invisible(NULL))
  if(is.null(output_file)) output_file <- sprintf("%s_method_comparison.pdf", index_name)
  comb <- rbind(saf_cond, saf_kend)
  p <- ggplot(comb, aes(x=Area_pct, y=Severity, colour=factor(ReturnPeriod_years), linetype=Method, shape=Method)) +
    geom_line(linewidth=1.2) + geom_point(size=3) + scale_color_viridis_d(name="Return Period\n(years)") +
    scale_linetype_manual(values=c(Conditional="solid", Kendall="dashed")) + scale_shape_manual(values=c(Conditional=19, Kendall=17)) +
    labs(x="Percent of Area Under Drought (%)", y="Drought Severity", title=sprintf("%s: Method Comparison", index_name)) + theme_bw()
  pdf(output_file, width=10, height=8); print(p); dev.off(); p
}

cross_correlate_indices <- function(dc1, dc3, output_dir=".") {
  m <- merge(dc1, dc3, by="date", suffixes=c("_1","_3"))
  r <- cor(m$severity_1, m$severity_3)
  pdf(file.path(output_dir, "SPI_cross_correlation.pdf"), width=8, height=5)
  plot(m$severity_1, m$severity_3, main=sprintf("Severity r=%.3f", r)); abline(0,1,lty=2,col="grey"); dev.off()
  m
}

cat("METHODOLOGY GUIDE LOADED (SPI/SPEI ONLY)\n")