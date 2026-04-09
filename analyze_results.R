setwd("C:/Users/yuhao/Documents/Equity Volatility Forecasting & Risk Estimation")

cat("=== TABLE 1: RMSE ===\n")
t1 <- read.csv("outputs/table1_rmse.csv")
print(t1)

cat("\n=== TABLE 2: Mincer-Zarnowitz ===\n")
t2 <- read.csv("outputs/table2_mz.csv")
print(t2)

cat("\n=== TABLE 2b: Encompassing ===\n")
t2e <- read.csv("outputs/table2_encompassing.csv")
print(t2e)

cat("\n=== TABLE 3: VaR Backtesting ===\n")
t3 <- read.csv("outputs/table3_var.csv")
print(t3)

cat("\n=== Encompassing by Horizon ===\n")
enc <- readRDS("outputs/encompassing_by_horizon.rds")
print(enc)

cat("\n=== Granger Causality ===\n")
gr <- readRDS("outputs/granger_results.rds")
cat("Returns -> dVIX  F=", round(gr$r_causes_dVIX$F_stat, 3),
    " p=", round(gr$r_causes_dVIX$pval, 4), "\n")
cat("dVIX -> Returns  F=", round(gr$dVIX_causes_r$F_stat, 3),
    " p=", round(gr$dVIX_causes_r$pval, 4), "\n")

cat("\n=== VRP Summary ===\n")
vrp <- readRDS("outputs/vrp_series.rds")
cat("Mean VRP:", round(mean(vrp$VRP, na.rm=TRUE), 2), "pp\n")
cat("Std  VRP:", round(sd(vrp$VRP,   na.rm=TRUE), 2), "pp\n")
cat("Min  VRP:", round(min(vrp$VRP,  na.rm=TRUE), 2), "pp\n")
cat("Max  VRP:", round(max(vrp$VRP,  na.rm=TRUE), 2), "pp\n")
cat("Pct positive:", round(100 * mean(vrp$VRP > 0, na.rm=TRUE), 1), "%\n")
