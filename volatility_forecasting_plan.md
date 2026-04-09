# Equity Volatility Forecasting & Risk Estimation — Implementation Plan

## Phase 0: Environment & Data

1. **Set up R with packages**: `rugarch`, `rmgarch`, `vars`, `forecast`, `xts`, `quantmod`, `PerformanceAnalytics`, `lmtest`, `car`, `sgt` (skewed generalized t), `ggplot2`.

2. **Download data**: Pull daily S&P 500 adjusted close prices from Yahoo Finance via `quantmod::getSymbols("^GSPC", from="1990-01-01", to="2024-12-31")`. Separately pull daily VIX close (`^VIX`).

3. **Compute returns**: Log returns `r_t = ln(P_t / P_{t-1})`. Drop the first NA row.

4. **Compute realized variance proxies**:
   - 1-day: squared daily return `r_t^2` (or use absolute return — document the choice).
   - 22-day: rolling 22-day realized variance = sum of squared returns over the next 22 trading days, divided by 22. This is the **forward-looking** target for the 22-day horizon forecast.

5. **Align VIX to variance scale**: VIX is annualized vol in percentage points. Convert to daily variance: `VIX_var_1d = (VIX/100)^2 / 252`. For 22-day: `VIX_var_22d = (VIX/100)^2 * (22/252)`.

6. **Train/test split**: Use the last ~6 years (roughly 2019–2024, ~1,510 obs) as holdout. Everything before is the estimation window.

---

## Phase 1: EWMA Baseline

7. **Implement EWMA(λ=0.94)**: Recursively compute `σ²_t = λ * σ²_{t-1} + (1-λ) * r²_{t-1}`, initialized at the sample variance of the first 60 returns.

8. **Generate out-of-sample 1-day forecasts**: On each holdout day, the EWMA forecast is just `σ²_t` (no re-estimation needed — it's fully recursive).

9. **For 22-day horizon**: Since EWMA has no mean-reversion, the h-step-ahead forecast is just `h * σ²_t`. Compute `22 * σ²_t` for each holdout day.

10. **Store forecasts** in a data frame aligned by date.

---

## Phase 2: GARCH Estimation & Forecasting

11. **Specify the model**: Use `rugarch::ugarchspec()` with `variance.model = list(model="sGARCH", garchOrder=c(2,1))` and `distribution.model = "sstd"` (skewed Student-t). Mean model: `armaOrder=c(0,0)` with `include.mean=TRUE`.

12. **Rolling out-of-sample forecasts**: Use `ugarchroll()` or a manual loop — on each holdout day t, estimate GARCH on all data up to t-1, then forecast:
    - 1-day ahead: `ugarchforecast(..., n.ahead=1)`.
    - 22-day ahead: `ugarchforecast(..., n.ahead=22)`, then sum the 22 conditional variance forecasts to get cumulative 22-day variance.
    - **Practical note**: Re-estimating daily on 34 years of data is expensive. A common shortcut is to re-estimate every 20–50 days and use `ugarchfilter` in between. Document whichever choice you make.

13. **Extract fitted distribution parameters** (shape ν, skew ξ) from each estimation window — you'll need these for VaR.

14. **Store GARCH 1-day and 22-day variance forecasts** alongside EWMA and VIX-implied.

---

## Phase 3: Forecast Evaluation

15. **RMSE comparison**: For each horizon (1-day, 22-day), compute RMSE of each model's variance forecast against realized variance. Report percentage improvement of GARCH over EWMA.

16. **Mincer-Zarnowitz regressions**: For each model m and horizon h, run `RV_t = α + β * Forecast_m_t + ε_t`. Test H₀: α=0, β=1 (joint F-test). An ideal forecast has R² near 1, α≈0, β≈1.

17. **Encompassing tests**: Run `RV_t = α + β₁ * GARCH_t + β₂ * VIX_t + ε_t`. If β₁ is insignificant while β₂ is significant, VIX "encompasses" GARCH at that horizon. Do this for both horizons to identify where GARCH adds incremental content.

18. **Vary the horizon**: Repeat the encompassing regression at h = 1, 5, 10, 22 days. Plot the t-statistic of β₁ (GARCH) as a function of horizon to show it fades beyond ~5 days.

---

## Phase 4: VaR Estimation & Backtesting

19. **Compute 99% 1-day VaR under two assumptions**:
    - **Normal**: `VaR_t = μ̂ + σ̂_t * qnorm(0.01)` where σ̂_t is the GARCH conditional SD.
    - **Skewed-t**: `VaR_t = μ̂ + σ̂_t * qdist("sstd", 0.01, shape=ν̂, skew=ξ̂)` using the fitted skew-t quantile from `rugarch`.

20. **Exceedance series**: Flag days where actual return < VaR. Compute empirical exceedance rate = count / holdout length. Target: 1%.

21. **Kupiec unconditional coverage test**: Likelihood ratio test comparing observed exceedance rate to 1%. Report the test statistic and p-value. A p-value > 0.05 means you cannot reject correct coverage.

22. **Christoffersen independence test** (bonus): Tests whether exceedances cluster. Combined conditional coverage test = Kupiec + independence.

---

## Phase 5: Volatility Risk Premium & Lead-Lag Dynamics

23. **Volatility risk premium**: Compute `VRP_t = VIX_t (annualized) − RV_t (annualized)` where RV is the subsequent 22-day realized vol (annualized). Report the time-series mean (~3–4 pts typically), plot it over time, note it widens in crises.

24. **Bivariate VAR**: Estimate a VAR in `(r_t, ΔVIX_t)` where `ΔVIX_t = VIX_t − VIX_{t-1}`.
    - Use `vars::VARselect()` to pick lag order by AIC/BIC.
    - Estimate with `vars::VAR()`.

25. **Granger causality tests**: `vars::causality(var_model, cause="returns")` and vice versa. You expect returns Granger-cause VIX changes (leverage effect / asymmetric volatility), and possibly weaker reverse causality.

26. **Impulse response functions**: Plot IRFs from the VAR — a negative return shock should produce a positive VIX response (asymmetry). Use `vars::irf()` with bootstrapped confidence bands.

---

## Phase 6: Outputs & Visualization

27. **Table 1**: Out-of-sample RMSE by model × horizon.

28. **Table 2**: Mincer-Zarnowitz regression coefficients + encompassing test results.

29. **Table 3**: VaR exceedance rates and Kupiec test p-values.

30. **Figure 1**: Time series of conditional volatility (GARCH σ_t vs. EWMA σ_t vs. VIX-implied σ_t) over the holdout period, with realized vol overlaid.

31. **Figure 2**: VaR exceedances — plot returns with VaR thresholds, flagging breach days.

32. **Figure 3**: Volatility risk premium over time.

33. **Figure 4**: IRF plots from the bivariate VAR.

34. **Figure 5**: GARCH encompassing t-statistic as a function of forecast horizon (1–22 days).

---

## Implementation Notes

- The most computationally expensive step is the rolling GARCH re-estimation (Phase 2, step 12). If the agent hits time constraints, re-estimating every 20 trading days and filtering in between is the standard compromise.
- The 22-day realized variance target must be strictly forward-looking; make sure the alignment doesn't accidentally include contemporaneous information.
- When computing VRP, ensure VIX and realized vol are on the same annualized scale before differencing.
- For Granger causality, stationarity of both series in the VAR is required — log returns are stationary by construction; first-differenced VIX should be as well, but run an ADF test to confirm.
