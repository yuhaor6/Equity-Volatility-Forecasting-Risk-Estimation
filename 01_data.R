# data download and preprocessing
# S&P 500 + VIX from Yahoo Finance (1990-2024), log returns, realized variance,
# VIX-implied variance. Train/test split at 2019-01-01.

library(quantmod)
library(xts)
library(zoo)

if (!dir.exists("data"))    dir.create("data")
if (!dir.exists("outputs")) dir.create("outputs")
if (!dir.exists("figures")) dir.create("figures")

getSymbols("^GSPC", src = "yahoo", from = "1990-01-01", to = "2024-12-31",
           auto.assign = TRUE, warnings = FALSE)
getSymbols("^VIX",  src = "yahoo", from = "1990-01-01", to = "2024-12-31",
           auto.assign = TRUE, warnings = FALSE)

spx_close <- Ad(GSPC)
vix_close  <- Cl(VIX)

returns_xts <- diff(log(spx_close))
returns_xts <- na.omit(returns_xts)
colnames(returns_xts) <- "r"

# 1-day realized variance
rv_1d <- returns_xts^2
colnames(rv_1d) <- "rv_1d"

# cumulative 22-day forward variance: sum r_{t+1}^2 ... r_{t+22}^2
# same scale as GARCH cumulative 22-step forecast and VIX-implied 22d var
r_sq <- coredata(returns_xts)^2
n    <- length(r_sq)
rv_22d_vec <- rep(NA_real_, n)
for (i in seq_len(n - 22))
  rv_22d_vec[i] <- sum(r_sq[(i + 1):(i + 22)])

rv_22d <- xts(rv_22d_vec, order.by = index(returns_xts))
colnames(rv_22d) <- "rv_22d"

# align VIX to return dates; fill gaps forward
vix_aligned <- merge(returns_xts, vix_close, join = "left")[, "VIX.Close"]
vix_aligned <- na.locf(vix_aligned)

vix_var_1d  <- (vix_aligned / 100)^2 / 252
vix_var_22d <- (vix_aligned / 100)^2 * (22 / 252)
colnames(vix_var_1d)  <- "vix_var_1d"
colnames(vix_var_22d) <- "vix_var_22d"

full_data     <- na.omit(merge(returns_xts, rv_1d, rv_22d, vix_var_1d, vix_var_22d))
HOLDOUT_START <- as.Date("2019-01-01")
train_data    <- full_data[index(full_data) <  HOLDOUT_START, ]
test_data     <- full_data[index(full_data) >= HOLDOUT_START, ]

cat(sprintf("train: %d obs  (%s to %s)\n", nrow(train_data), start(train_data), end(train_data)))
cat(sprintf("test:  %d obs  (%s to %s)\n", nrow(test_data),  start(test_data),  end(test_data)))

saveRDS(returns_xts, "data/returns_xts.rds")
saveRDS(vix_close,   "data/vix_close.rds")
saveRDS(rv_1d,       "data/rv_1d.rds")
saveRDS(rv_22d,      "data/rv_22d.rds")
saveRDS(vix_var_1d,  "data/vix_var_1d.rds")
saveRDS(vix_var_22d, "data/vix_var_22d.rds")
saveRDS(full_data,   "data/full_data.rds")
saveRDS(train_data,  "data/train_data.rds")
saveRDS(test_data,   "data/test_data.rds")
