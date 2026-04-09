# rolling GARCH(2,1)-sstd out-of-sample forecasts
# refits every REFIT_EVERY days; ugarchforecast() filters through current params
# on non-refit days. stores 1-day and cumulative 22-day variance, plus nu and xi.

library(rugarch)
library(xts)

full_data     <- readRDS("data/full_data.rds")
test_data     <- readRDS("data/test_data.rds")

returns_all   <- full_data[, "r"]
HOLDOUT_START <- start(test_data)
holdout_idx   <- which(index(returns_all) >= HOLDOUT_START)
n_holdout     <- length(holdout_idx)

garch_spec <- ugarchspec(
  variance.model = list(
    model      = "sGARCH",
    garchOrder = c(2, 1)
  ),
  mean.model = list(
    armaOrder    = c(0, 0),
    include.mean = TRUE
  ),
  distribution.model = "sstd"
)

REFIT_EVERY <- 20
MAX_HORIZON <- 22

garch_var_1d  <- numeric(n_holdout)
garch_var_22d <- numeric(n_holdout)
nu_vec        <- numeric(n_holdout)
xi_vec        <- numeric(n_holdout)

current_fit <- NULL
set.seed(2019)

cat(sprintf("rolling GARCH: %d holdout days, refit every %d\n", n_holdout, REFIT_EVERY))

for (k in seq_len(n_holdout)) {
  t        <- holdout_idx[k]
  r_window <- returns_all[1:(t - 1), ]

  if (is.null(current_fit) || k %% REFIT_EVERY == 1) {
    tryCatch({
      current_fit <- ugarchfit(spec = garch_spec, data = r_window,
                               solver = "hybrid", solver.control = list(trace = 0))
    }, error = function(e) NULL)
  }

  if (!is.null(current_fit)) {
    fc <- tryCatch(
      ugarchforecast(current_fit, data = r_window, n.ahead = MAX_HORIZON),
      error = function(e) NULL
    )
    if (!is.null(fc)) {
      var_fc           <- as.numeric(sigma(fc))^2
      garch_var_1d[k]  <- var_fc[1]
      garch_var_22d[k] <- sum(var_fc)
      cf               <- coef(current_fit)
      nu_vec[k]        <- cf["shape"]
      xi_vec[k]        <- cf["skew"]
    } else {
      garch_var_1d[k]  <- as.numeric(tail(r_window, 1))^2
      garch_var_22d[k] <- 22 * garch_var_1d[k]
      nu_vec[k]        <- if (k > 1) nu_vec[k - 1] else 8
      xi_vec[k]        <- if (k > 1) xi_vec[k - 1] else 1
    }
  }

  if (k %% 100 == 0) cat(sprintf("  %d / %d\n", k, n_holdout))
}

garch_forecasts <- data.frame(
  date          = index(test_data),
  garch_var_1d  = garch_var_1d,
  garch_var_22d = garch_var_22d,
  nu            = nu_vec,
  xi            = xi_vec
)

saveRDS(garch_forecasts, "data/garch_forecasts.rds")
cat(sprintf("done. mean ann vol = %.1f%%\n",
            sqrt(mean(garch_forecasts$garch_var_1d) * 252) * 100))
