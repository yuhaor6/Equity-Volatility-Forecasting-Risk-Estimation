# VaR estimation and backtesting: Kupiec + Christoffersen tests

library(rugarch)
library(xts)

holdout         <- readRDS("data/holdout_merged.rds")
garch_forecasts <- readRDS("data/garch_forecasts.rds")

g <- merge(holdout[, c("date", "r")], garch_forecasts, by = "date")

mu_hat <- mean(holdout$r)

var_normal <- mu_hat + sqrt(g$garch_var_1d) * qnorm(0.01)

var_sstd <- mapply(
  function(sigma, nu, xi) {
    mu_hat + sigma * tryCatch(
      qdist("sstd", p = 0.01, mu = 0, sigma = 1, shape = nu, skew = xi),
      error = function(e) qnorm(0.01)
    )
  },
  sigma = sqrt(g$garch_var_1d),
  nu    = g$nu,
  xi    = g$xi
)

excess_normal <- as.integer(g$r < var_normal)
excess_sstd   <- as.integer(g$r < var_sstd)

rate_normal <- mean(excess_normal)
rate_sstd   <- mean(excess_sstd)
n           <- nrow(g)

cat(sprintf("Normal VaR  exceedance rate: %.4f (target 0.01, n=%d)\n", rate_normal, n))
cat(sprintf("Skewed-t VaR exceedance rate: %.4f (target 0.01, n=%d)\n", rate_sstd,   n))

# --------------------------------------------------------------------------- #
kupiec <- function(exceedances, n, p0 = 0.01) {
  T1   <- sum(exceedances)
  T0   <- n - T1
  p_hat <- T1 / n
  if (p_hat == 0 || p_hat == 1) {
    return(data.frame(T1 = T1, p_hat = p_hat, LR = NA, pval = NA))
  }
  LR <- -2 * (T1 * log(p0) + T0 * log(1 - p0) -
              T1 * log(p_hat) - T0 * log(1 - p_hat))
  pval <- pchisq(LR, df = 1, lower.tail = FALSE)
  data.frame(T1 = T1, p_hat = round(p_hat, 4), LR = round(LR, 3), pval = round(pval, 4))
}

kup_normal <- kupiec(excess_normal, n)
kup_sstd   <- kupiec(excess_sstd,   n)

christoffersen <- function(exceedances) {
  n00 <- sum(exceedances[-length(exceedances)] == 0 & exceedances[-1] == 0)
  n01 <- sum(exceedances[-length(exceedances)] == 0 & exceedances[-1] == 1)
  n10 <- sum(exceedances[-length(exceedances)] == 1 & exceedances[-1] == 0)
  n11 <- sum(exceedances[-length(exceedances)] == 1 & exceedances[-1] == 1)

  pi_hat  <- (n01 + n11) / (n00 + n01 + n10 + n11)
  pi01    <- n01 / (n00 + n01)
  pi11    <- if ((n10 + n11) > 0) n11 / (n10 + n11) else 0

  if (pi_hat == 0 || pi_hat == 1 || pi01 == 0 || pi01 == 1) {
    return(data.frame(LR_ind = NA, pval_ind = NA))
  }

  L_pi  <- (1 - pi_hat)^(n00 + n10) * pi_hat^(n01 + n11)
  L_01  <- tryCatch(
    (1 - pi01)^n00 * pi01^n01 * (1 - pi11)^n10 * pi11^n11,
    warning = function(w) NA
  )
  if (is.na(L_01) || L_01 == 0) return(data.frame(LR_ind = NA, pval_ind = NA))

  LR_ind  <- -2 * (log(L_pi) - log(L_01))
  pval_ind <- pchisq(LR_ind, df = 1, lower.tail = FALSE)
  data.frame(LR_ind = round(LR_ind, 3), pval_ind = round(pval_ind, 4))
}

chr_normal <- christoffersen(excess_normal)
chr_sstd   <- christoffersen(excess_sstd)

table3 <- data.frame(
  model        = c("GARCH-Normal", "GARCH-sstd"),
  exceedances  = c(kup_normal$T1,    kup_sstd$T1),
  exc_rate     = c(kup_normal$p_hat, kup_sstd$p_hat),
  kupiec_LR    = c(kup_normal$LR,    kup_sstd$LR),
  kupiec_pval  = c(kup_normal$pval,  kup_sstd$pval),
  chr_LR_ind   = c(chr_normal$LR_ind,  chr_sstd$LR_ind),
  chr_pval_ind = c(chr_normal$pval_ind, chr_sstd$pval_ind)
)

print(table3)

var_series <- data.frame(
  date       = g$date,
  r          = g$r,
  var_normal = var_normal,
  var_sstd   = var_sstd,
  exc_normal = excess_normal,
  exc_sstd   = excess_sstd
)

saveRDS(table3,     "outputs/table3_var.rds")
saveRDS(var_series, "outputs/var_series.rds")
