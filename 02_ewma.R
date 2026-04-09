# EWMA(lambda=0.94) baseline forecasts
# sigma2_t = lam*sigma2_{t-1} + (1-lam)*r_{t-1}^2, initialised on first 60 obs
# 22-day forecast = 22 * sigma2_t  (no mean reversion)

library(xts)

full_data  <- readRDS("data/full_data.rds")
test_data  <- readRDS("data/test_data.rds")

returns_all   <- full_data[, "r"]
HOLDOUT_START <- start(test_data)

lam    <- 0.94
r_vec  <- coredata(returns_all)[, 1]
n      <- length(r_vec)
sigma2 <- numeric(n)
sigma2[1] <- var(r_vec[1:60])

for (t in 2:n)
  sigma2[t] <- lam * sigma2[t - 1] + (1 - lam) * r_vec[t - 1]^2

ewma_var <- xts(sigma2, order.by = index(returns_all))
colnames(ewma_var) <- "ewma_var_1d"

ewma_var_22d <- 22 * ewma_var
colnames(ewma_var_22d) <- "ewma_var_22d"

ewma_hold <- merge(ewma_var, ewma_var_22d)[index(ewma_var) >= HOLDOUT_START]

ewma_forecasts <- data.frame(
  date         = index(ewma_hold),
  ewma_var_1d  = coredata(ewma_hold[, "ewma_var_1d"]),
  ewma_var_22d = coredata(ewma_hold[, "ewma_var_22d"])
)
rownames(ewma_forecasts) <- NULL

saveRDS(ewma_forecasts, "data/ewma_forecasts.rds")
