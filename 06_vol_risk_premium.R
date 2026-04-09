# Volatility risk premium and bivariate VAR(r, dVIX)

library(xts)
library(zoo)
library(vars)
library(tseries)

full_data  <- readRDS("data/full_data.rds")
rv_22d     <- readRDS("data/rv_22d.rds")
vix_close  <- readRDS("data/vix_close.rds")

returns_all <- full_data[, "r"]

# VRP = VIX (annualised %) - realised 22-day vol (annualised %)
vix_aligned <- merge(returns_all, vix_close, join = "left")[, "VIX.Close"]
vix_aligned <- na.locf(vix_aligned)

ann_rv_22d <- sqrt((rv_22d / 22) * 252) * 100
colnames(ann_rv_22d) <- "ann_rv_22d"

vrp_xts <- merge(vix_aligned, ann_rv_22d, join = "inner")
vrp_xts  <- na.omit(vrp_xts)
vrp_xts$VRP <- vrp_xts[, 1] - vrp_xts[, 2]

cat(sprintf("VRP sample: %s to %s\n", start(vrp_xts), end(vrp_xts)))
cat(sprintf("  Mean VRP: %.2f pp\n", mean(coredata(vrp_xts$VRP))))
cat(sprintf("  Std  VRP: %.2f pp\n", sd(coredata(vrp_xts$VRP))))

vrp_df <- data.frame(
  date       = index(vrp_xts),
  VIX        = coredata(vrp_xts)[, 1],
  ann_rv_22d = coredata(vrp_xts)[, 2],
  VRP        = coredata(vrp_xts$VRP)
)
saveRDS(vrp_df, "outputs/vrp_series.rds")

# Bivariate VAR in (r_t, dVIX_t)
delta_vix <- diff(vix_aligned)
delta_vix <- na.omit(delta_vix)
colnames(delta_vix) <- "dVIX"

# Align returns and dVIX
var_data <- merge(returns_all, delta_vix, join = "inner")
var_data <- na.omit(var_data)
colnames(var_data) <- c("r", "dVIX")

# ADF stationarity check
adf_r    <- adf.test(coredata(var_data[, "r"]))
adf_dvix <- adf.test(coredata(var_data[, "dVIX"]))

cat(sprintf("\nADF log-returns: stat=%.3f  p=%.4f\n", adf_r$statistic, adf_r$p.value))
cat(sprintf("ADF dVIX:        stat=%.3f  p=%.4f\n", adf_dvix$statistic, adf_dvix$p.value))

var_mat <- as.matrix(var_data)
lag_sel <- VARselect(var_mat, lag.max = 10, type = "const")
opt_lag <- lag_sel$selection["SC(n)"]

cat(sprintf("\nVAR lag: BIC=%d\n", opt_lag))

# Estimate VAR
var_model <- VAR(var_mat, p = opt_lag, type = "const")

gc_r_causes_dvix <- causality(var_model, cause = "r")
gc_dvix_causes_r <- causality(var_model, cause = "dVIX")

granger_results <- list(
  r_causes_dVIX = list(
    F_stat = gc_r_causes_dvix$Granger$statistic,
    pval   = gc_r_causes_dvix$Granger$p.value
  ),
  dVIX_causes_r = list(
    F_stat = gc_dvix_causes_r$Granger$statistic,
    pval   = gc_dvix_causes_r$Granger$p.value
  )
)

cat(sprintf("Returns → dVIX: F=%.3f  p=%.4f\n",
            granger_results$r_causes_dVIX$F_stat,
            granger_results$r_causes_dVIX$pval))
cat(sprintf("dVIX → Returns: F=%.3f  p=%.4f\n",
            granger_results$dVIX_causes_r$F_stat,
            granger_results$dVIX_causes_r$pval))

irf_r_to_dvix <- irf(var_model, impulse = "r", response = "dVIX",
                     n.ahead = 20, boot = TRUE, ci = 0.95, runs = 200)

irf_dvix_to_r <- irf(var_model, impulse = "dVIX", response = "r",
                     n.ahead = 20, boot = TRUE, ci = 0.95, runs = 200)

saveRDS(var_model,       "outputs/var_model.rds")
saveRDS(granger_results, "outputs/granger_results.rds")
saveRDS(irf_r_to_dvix,  "outputs/irf_r_to_dvix.rds")
saveRDS(irf_dvix_to_r,  "outputs/irf_dvix_to_r.rds")
