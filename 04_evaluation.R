# forecast evaluation: RMSE, Mincer-Zarnowitz regressions, encompassing tests

library(xts)
library(lmtest)
library(car)
library(dplyr)

test_data       <- readRDS("data/test_data.rds")
ewma_forecasts  <- readRDS("data/ewma_forecasts.rds")
garch_forecasts <- readRDS("data/garch_forecasts.rds")
full_data       <- readRDS("data/full_data.rds")
returns_all     <- full_data[, "r"]
HOLDOUT_START   <- start(test_data)

holdout <- data.frame(
  date        = index(test_data),
  r           = coredata(test_data[, "r"]),
  rv_1d       = coredata(test_data[, "rv_1d"]),
  rv_22d      = coredata(test_data[, "rv_22d"]),
  vix_var_1d  = coredata(test_data[, "vix_var_1d"]),
  vix_var_22d = coredata(test_data[, "vix_var_22d"])
) |>
  merge(ewma_forecasts, by = "date") |>
  merge(garch_forecasts[, c("date", "garch_var_1d", "garch_var_22d")], by = "date")

holdout <- holdout[complete.cases(holdout) & holdout$rv_22d > 0 & holdout$garch_var_22d > 0, ]

# RMSE
rmse <- function(f, a) sqrt(mean((f - a)^2))

rmse_table <- data.frame(
  model    = c("EWMA", "GARCH", "VIX"),
  rmse_1d  = c(rmse(holdout$ewma_var_1d,  holdout$rv_1d),
               rmse(holdout$garch_var_1d, holdout$rv_1d),
               rmse(holdout$vix_var_1d,   holdout$rv_1d)),
  rmse_22d = c(rmse(holdout$ewma_var_22d,  holdout$rv_22d),
               rmse(holdout$garch_var_22d, holdout$rv_22d),
               rmse(holdout$vix_var_22d,   holdout$rv_22d))
)
rmse_table$pct_impr_over_ewma_1d  <- with(rmse_table,
  round(100 * (rmse_1d[model == "EWMA"]  - rmse_1d)  / rmse_1d[model == "EWMA"],  2))
rmse_table$pct_impr_over_ewma_22d <- with(rmse_table,
  round(100 * (rmse_22d[model == "EWMA"] - rmse_22d) / rmse_22d[model == "EWMA"], 2))

print(rmse_table)

# Mincer-Zarnowitz: RV = a + b*forecast  H0: a=0, b=1
mz_reg <- function(rv, fc, label) {
  fit <- lm(rv ~ fc, data = data.frame(rv = rv, fc = fc))
  cf  <- coef(fit)
  jt  <- tryCatch(
    linearHypothesis(fit, c("(Intercept) = 0", "fc = 1"), test = "F"),
    error = function(e) NULL
  )
  data.frame(
    label  = label,
    alpha  = round(cf[1], 8),
    beta   = round(cf[2], 4),
    R2     = round(summary(fit)$r.squared, 4),
    F_stat = if (!is.null(jt)) round(jt$F[2], 3)        else NA,
    F_pval = if (!is.null(jt)) round(jt$`Pr(>F)`[2], 4) else NA
  )
}

mz_results <- rbind(
  mz_reg(holdout$rv_1d,  holdout$ewma_var_1d,   "EWMA_1d"),
  mz_reg(holdout$rv_1d,  holdout$garch_var_1d,  "GARCH_1d"),
  mz_reg(holdout$rv_1d,  holdout$vix_var_1d,    "VIX_1d"),
  mz_reg(holdout$rv_22d, holdout$ewma_var_22d,  "EWMA_22d"),
  mz_reg(holdout$rv_22d, holdout$garch_var_22d, "GARCH_22d"),
  mz_reg(holdout$rv_22d, holdout$vix_var_22d,   "VIX_22d"))
print(mz_results)

# encompassing: RV = a + b1*GARCH + b2*VIX
encompassing <- function(rv, fc_garch, fc_vix, label) {
  fit <- lm(rv ~ garch + vix,
            data = data.frame(rv = rv, garch = fc_garch, vix = fc_vix))
  cf  <- coef(summary(fit))
  data.frame(
    label       = label,
    beta1_garch = round(cf["garch", "Estimate"],  4),
    t_garch     = round(cf["garch", "t value"],   3),
    p_garch     = round(cf["garch", "Pr(>|t|)"],  4),
    beta2_vix   = round(cf["vix",   "Estimate"],  4),
    t_vix       = round(cf["vix",   "t value"],   3),
    p_vix       = round(cf["vix",   "Pr(>|t|)"],  4),
    R2          = round(summary(fit)$r.squared,   4)
  )
}

enc_results <- rbind(
  encompassing(holdout$rv_1d,  holdout$garch_var_1d,  holdout$vix_var_1d,  "1d"),
  encompassing(holdout$rv_22d, holdout$garch_var_22d, holdout$vix_var_22d, "22d"))
print(enc_results)

# encompassing at h = 1, 5, 10, 22 days
# intermediate GARCH horizons approximated by linear interpolation of cumulative var
r_sq_all    <- coredata(returns_all)^2
dates       <- index(returns_all)
n_all       <- length(r_sq_all)
holdout_pos <- which(dates >= HOLDOUT_START)

enc_by_h <- lapply(c(1, 5, 10, 22), function(h) {
  rv_h  <- sapply(holdout_pos, function(i)
    if (i + h <= n_all) sum(r_sq_all[(i + 1):(i + h)]) else NA)
  vix_h <- coredata(test_data[, "vix_var_1d"]) * h
  g1    <- garch_forecasts$garch_var_1d
  g22   <- garch_forecasts$garch_var_22d
  gh    <- pmax(g1 * h + (g22 - g1 * 22) * (h - 1) / 21, .Machine$double.eps)
  keep  <- !is.na(rv_h) & !is.na(vix_h) & rv_h > 0
  fit   <- lm(rv ~ garch + vix,
              data = data.frame(rv    = rv_h[keep],
                                garch = gh[keep],
                                vix   = as.numeric(vix_h)[keep]))
  cf <- coef(summary(fit))
  data.frame(h       = h,
             t_garch = round(cf["garch", "t value"],  3),
             p_garch = round(cf["garch", "Pr(>|t|)"], 4),
             t_vix   = round(cf["vix",   "t value"],  3),
             p_vix   = round(cf["vix",   "Pr(>|t|)"], 4))
})

enc_by_horizon <- do.call(rbind, enc_by_h)
print(enc_by_horizon)

saveRDS(rmse_table,     "outputs/table1_rmse.rds")
saveRDS(mz_results,     "outputs/table2_mz.rds")
saveRDS(enc_results,    "outputs/table2_encompassing.rds")
saveRDS(enc_by_horizon, "outputs/encompassing_by_horizon.rds")
saveRDS(holdout,        "data/holdout_merged.rds")

# --------------------------------------------------------------------------- #
# Load data
# --------------------------------------------------------------------------- #
test_data        <- readRDS("data/test_data.rds")
ewma_forecasts   <- readRDS("data/ewma_forecasts.rds")
garch_forecasts  <- readRDS("data/garch_forecasts.rds")
full_data        <- readRDS("data/full_data.rds")
returns_all      <- full_data[, "r"]

HOLDOUT_START <- start(test_data)

# --------------------------------------------------------------------------- #
# Build master holdout data frame
# --------------------------------------------------------------------------- #
holdout <- data.frame(
  date          = index(test_data),
  r             = coredata(test_data[, "r"]),
  rv_1d         = coredata(test_data[, "rv_1d"]),
  rv_22d        = coredata(test_data[, "rv_22d"]),
  vix_var_1d    = coredata(test_data[, "vix_var_1d"]),
  vix_var_22d   = coredata(test_data[, "vix_var_22d"]),
  stringsAsFactors = FALSE
) |>
  merge(ewma_forecasts,  by = "date") |>
  merge(garch_forecasts[, c("date", "garch_var_1d", "garch_var_22d")], by = "date")

# Drop rows where any forecast or realized vol is NA/zero
holdout <- holdout[complete.cases(holdout) &
                   holdout$rv_22d > 0 &
                   holdout$garch_var_22d > 0, ]

N <- nrow(holdout)
message(sprintf("Holdout rows used for evaluation: %d", N))

# --------------------------------------------------------------------------- #
# Step 15: RMSE by model × horizon
# --------------------------------------------------------------------------- #
rmse <- function(forecast, actual) sqrt(mean((forecast - actual)^2))

rmse_table <- data.frame(
  model    = c("EWMA", "GARCH", "VIX"),
  rmse_1d  = c(
    rmse(holdout$ewma_var_1d,  holdout$rv_1d),
    rmse(holdout$garch_var_1d, holdout$rv_1d),
    rmse(holdout$vix_var_1d,   holdout$rv_1d)
  ),
  rmse_22d = c(
    rmse(holdout$ewma_var_22d,  holdout$rv_22d),
    rmse(holdout$garch_var_22d, holdout$rv_22d),
    rmse(holdout$vix_var_22d,   holdout$rv_22d)
  )
)

rmse_table$pct_impr_over_ewma_1d <- with(rmse_table,
  round(100 * (rmse_1d[model == "EWMA"] - rmse_1d) /
              rmse_1d[model == "EWMA"], 2))

rmse_table$pct_impr_over_ewma_22d <- with(rmse_table,
  round(100 * (rmse_22d[model == "EWMA"] - rmse_22d) /
              rmse_22d[model == "EWMA"], 2))

message("\n--- Table 1: RMSE by model × horizon ---")
print(rmse_table)

# --------------------------------------------------------------------------- #
# Step 16: Mincer-Zarnowitz regressions
# RV_t = α + β * Forecast_t + ε_t   (H₀: α=0, β=1)
# --------------------------------------------------------------------------- #
mz_reg <- function(rv, fc, label) {
  dat <- data.frame(rv = rv, fc = fc)
  fit <- lm(rv ~ fc, data = dat)
  cf  <- coef(fit)
  r2  <- summary(fit)$r.squared

  # Joint F-test H₀: α=0, β=1
  jt  <- tryCatch(
    linearHypothesis(fit, c("(Intercept) = 0", "fc = 1"), test = "F"),
    error = function(e) NULL
  )
  f_stat <- if (!is.null(jt)) jt$F[2]    else NA
  f_pval <- if (!is.null(jt)) jt$`Pr(>F)`[2] else NA

  data.frame(
    label    = label,
    alpha    = round(cf[1], 8),
    beta     = round(cf[2], 4),
    R2       = round(r2, 4),
    F_stat   = round(f_stat, 3),
    F_pval   = round(f_pval, 4),
    stringsAsFactors = FALSE
  )
}

mz_results <- rbind(
  mz_reg(holdout$rv_1d,  holdout$ewma_var_1d,  "EWMA_1d"),
  mz_reg(holdout$rv_1d,  holdout$garch_var_1d, "GARCH_1d"),
  mz_reg(holdout$rv_1d,  holdout$vix_var_1d,   "VIX_1d"),
  mz_reg(holdout$rv_22d, holdout$ewma_var_22d,  "EWMA_22d"),
  mz_reg(holdout$rv_22d, holdout$garch_var_22d, "GARCH_22d"),
  mz_reg(holdout$rv_22d, holdout$vix_var_22d,   "VIX_22d")
)

# --------------------------------------------------------------------------- #
# Step 17: Encompassing tests
# RV_t = α + β₁ * GARCH_t + β₂ * VIX_t + ε_t
# --------------------------------------------------------------------------- #
encompassing <- function(rv, fc_garch, fc_vix, label) {
  dat <- data.frame(rv = rv, garch = fc_garch, vix = fc_vix)
  fit <- lm(rv ~ garch + vix, data = dat)
  cf  <- coef(summary(fit))
  data.frame(
    label          = label,
    beta1_garch    = round(cf["garch", "Estimate"], 4),
    t_garch        = round(cf["garch", "t value"], 3),
    p_garch        = round(cf["garch", "Pr(>|t|)"], 4),
    beta2_vix      = round(cf["vix", "Estimate"], 4),
    t_vix          = round(cf["vix", "t value"], 3),
    p_vix          = round(cf["vix", "Pr(>|t|)"], 4),
    R2             = round(summary(fit)$r.squared, 4),
    stringsAsFactors = FALSE
  )
}

enc_results <- rbind(
  encompassing(holdout$rv_1d,  holdout$garch_var_1d,  holdout$vix_var_1d,  "1d"),
  encompassing(holdout$rv_22d, holdout$garch_var_22d, holdout$vix_var_22d, "22d")
)

message("\n--- Table 2: Mincer-Zarnowitz ---")
print(mz_results)
message("\n--- Encompassing tests (GARCH vs VIX) ---")
print(enc_results)

# --------------------------------------------------------------------------- #
# Step 18: Encompassing at multiple horizons h = 1, 5, 10, 22
# For each h, build forward-looking RV_h and a GARCH h-day variance forecast.
# We reuse the step-h GARCH sigmas stored per row. Since we only stored
# garch_var_1d and garch_var_22d (not intermediate horizons), we approximate
# intermediate horizons by linear interpolation of cumulative variance.
# For a more exact version one would run ugarchforecast for each h.
# --------------------------------------------------------------------------- #
r_sq_all <- coredata(returns_all)^2
r_vec    <- coredata(returns_all)[, 1]
dates    <- index(returns_all)
n_all    <- length(r_vec)
holdout_pos <- which(dates >= HOLDOUT_START)

horizons <- c(1, 5, 10, 22)

enc_by_h <- lapply(horizons, function(h) {
  # Forward RV_h for each holdout day: strictly forward-looking sum of r²
  # indices (i+1):(i+h) matches the rv_22d convention in 01_data.R
  rv_h <- sapply(holdout_pos, function(i) {
    if (i + h <= n_all) sum(r_sq_all[(i + 1):(i + h)]) else NA
  })

  # VIX-implied variance at horizon h
  vix_h <- coredata(test_data[, "vix_var_1d"]) * h   # scale 1d VIX var to h days

  # GARCH h-day: interpolate between garch_var_1d (h=1) and garch_var_22d (h=22)
  g1  <- garch_forecasts$garch_var_1d
  g22 <- garch_forecasts$garch_var_22d
  garch_h <- g1 * h + (g22 - g1 * 22) * (h - 1) / (22 - 1)
  garch_h <- pmax(garch_h, .Machine$double.eps)   # guard against negatives

  keep <- !is.na(rv_h) & !is.na(vix_h) & !is.na(garch_h) & rv_h > 0
  dat  <- data.frame(rv = rv_h[keep],
                     garch = garch_h[keep],
                     vix   = as.numeric(vix_h)[keep])
  fit  <- lm(rv ~ garch + vix, data = dat)
  cf   <- coef(summary(fit))
  data.frame(
    h       = h,
    t_garch = round(cf["garch", "t value"], 3),
    p_garch = round(cf["garch", "Pr(>|t|)"], 4),
    t_vix   = round(cf["vix",   "t value"], 3),
    p_vix   = round(cf["vix",   "Pr(>|t|)"], 4)
  )
})

enc_by_horizon <- do.call(rbind, enc_by_h)
message("\n--- Encompassing t-stats by horizon ---")
print(enc_by_horizon)

# --------------------------------------------------------------------------- #
# Save outputs
# --------------------------------------------------------------------------- #
saveRDS(rmse_table,        "outputs/table1_rmse.rds")
saveRDS(mz_results,        "outputs/table2_mz.rds")
saveRDS(enc_results,       "outputs/table2_encompassing.rds")
saveRDS(enc_by_horizon,    "outputs/encompassing_by_horizon.rds")
saveRDS(holdout,           "data/holdout_merged.rds")   # convenient for later phases

message("\nPhase 3 complete. Evaluation outputs saved to outputs/")
