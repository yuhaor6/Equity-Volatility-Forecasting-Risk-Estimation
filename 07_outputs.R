# Write CSV summary tables and produce all figures

library(ggplot2)
library(patchwork)
library(xts)
library(vars)
library(dplyr)

# --------------------------------------------------------------------------- #
# Load all output objects
# --------------------------------------------------------------------------- #
table1          <- readRDS("outputs/table1_rmse.rds")
table2_mz       <- readRDS("outputs/table2_mz.rds")
table2_enc      <- readRDS("outputs/table2_encompassing.rds")
table3          <- readRDS("outputs/table3_var.rds")
holdout         <- readRDS("data/holdout_merged.rds")
ewma_f          <- readRDS("data/ewma_forecasts.rds")
garch_f         <- readRDS("data/garch_forecasts.rds")
vrp_df          <- readRDS("outputs/vrp_series.rds")
var_series      <- readRDS("outputs/var_series.rds")
irf_r_to_dvix   <- readRDS("outputs/irf_r_to_dvix.rds")
enc_by_horizon  <- readRDS("outputs/encompassing_by_horizon.rds")

write.csv(table1,     "outputs/table1_rmse.csv",         row.names = FALSE)
write.csv(table2_mz,  "outputs/table2_mz.csv",           row.names = FALSE)
write.csv(table2_enc, "outputs/table2_encompassing.csv",  row.names = FALSE)
write.csv(table3,     "outputs/table3_var.csv",           row.names = FALSE)

HOLDOUT_START <- min(holdout$date)

save_png <- function(p, path, w = 12, h = 6, dpi = 180) {
  ggsave(path, plot = p, width = w, height = h, dpi = dpi)
}

vol_df <- holdout |>
  select(date, rv_1d) |>
  mutate(realized = sqrt(rv_1d * 252) * 100) |>
  left_join(ewma_f  |> transmute(date, ewma  = sqrt(ewma_var_1d  * 252) * 100), by = "date") |>
  left_join(garch_f |> transmute(date, garch = sqrt(garch_var_1d * 252) * 100), by = "date") |>
  mutate(vix = sqrt(holdout$vix_var_1d * 252) * 100) |>
  tidyr::pivot_longer(cols = c(realized, ewma, garch, vix),
                      names_to = "series", values_to = "ann_vol")

fig1 <- ggplot(vol_df, aes(x = date, y = ann_vol, colour = series, linewidth = series)) +
  geom_line(alpha = 0.85) +
  scale_colour_manual(values = c(realized = "black", ewma = "#E69F00",
                                 garch = "#0072B2", vix = "#D55E00"),
                      labels = c(realized = "Realized", ewma = "EWMA",
                                 garch = "GARCH(2,1)", vix = "VIX-implied")) +
  scale_linewidth_manual(values = c(realized = 0.5, ewma = 0.7, garch = 0.7, vix = 0.7),
                         guide = "none") +
  labs(title = "Conditional Volatility — Holdout Period",
       x = NULL, y = "Annualised Volatility (%)", colour = "Model") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

save_png(fig1, "figures/fig1_conditional_vol.png")

var_df <- var_series |>
  mutate(date = as.Date(date))

fig2_normal <- ggplot(var_df, aes(x = date)) +
  geom_col(aes(y = r * 100), fill = "steelblue", alpha = 0.5, width = 1) +
  geom_line(aes(y = var_normal * 100), colour = "red", linewidth = 0.6) +
  geom_point(data = var_df[var_df$exc_normal == 1, ],
             aes(y = r * 100), colour = "red", size = 1.2) +
  labs(title = "Normal VaR (99%)", x = NULL, y = "Return (%)") +
  theme_bw(base_size = 10)

fig2_sstd <- ggplot(var_df, aes(x = date)) +
  geom_col(aes(y = r * 100), fill = "steelblue", alpha = 0.5, width = 1) +
  geom_line(aes(y = var_sstd * 100), colour = "darkred", linewidth = 0.6) +
  geom_point(data = var_df[var_df$exc_sstd == 1, ],
             aes(y = r * 100), colour = "darkred", size = 1.2) +
  labs(title = "Skewed-t VaR (99%)", x = NULL, y = "Return (%)") +
  theme_bw(base_size = 10)

fig2 <- fig2_normal / fig2_sstd +
  plot_annotation(title = "99% 1-Day VaR — Returns with Breach Days (red)")

save_png(fig2, "figures/fig2_var_exceedances.png", h = 8)

vrp_plot_df <- vrp_df |>
  mutate(date = as.Date(date))

fig3 <- ggplot(vrp_plot_df, aes(x = date, y = VRP)) +
  geom_line(colour = "#0072B2", alpha = 0.8, linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = mean(vrp_plot_df$VRP, na.rm = TRUE),
             linetype = "dotdash", colour = "red", linewidth = 0.7) +
  annotate("text", x = max(vrp_plot_df$date, na.rm = TRUE),
           y = mean(vrp_plot_df$VRP, na.rm = TRUE) + 0.5,
           label = sprintf("Mean = %.2f pp", mean(vrp_plot_df$VRP, na.rm = TRUE)),
           hjust = 1, size = 3.5, colour = "red") +
  labs(title = "Volatility Risk Premium (VIX − Realised 22-day Vol)",
       x = NULL, y = "VRP (annualised percentage points)") +
  theme_bw(base_size = 11)

save_png(fig3, "figures/fig3_vrp.png")

irf_df <- data.frame(
  h        = 0:(length(irf_r_to_dvix$irf$r) - 1),
  irf      = as.numeric(irf_r_to_dvix$irf$r),
  irf_lo   = as.numeric(irf_r_to_dvix$Lower$r),
  irf_hi   = as.numeric(irf_r_to_dvix$Upper$r)
)

fig4 <- ggplot(irf_df, aes(x = h)) +
  geom_ribbon(aes(ymin = irf_lo, ymax = irf_hi), fill = "steelblue", alpha = 0.25) +
  geom_line(aes(y = irf), colour = "steelblue", linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "IRF — Return Shock → ΔVIX Response (95% bootstrap CI)",
       x = "Horizon (days)", y = "Response of ΔVIX") +
  theme_bw(base_size = 11)

save_png(fig4, "figures/fig4_irf.png", w = 8, h = 5)

fig5 <- ggplot(enc_by_horizon, aes(x = h, y = t_garch)) +
  geom_line(colour = "#0072B2", linewidth = 1) +
  geom_point(colour = "#0072B2", size = 3) +
  geom_hline(yintercept =  1.96, linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = -1.96, linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = 0,     linetype = "solid",  colour = "grey30") +
  annotate("text", x = max(enc_by_horizon$h), y = 1.96 + 0.15,
           label = "±1.96", hjust = 1, size = 3.5, colour = "grey50") +
  scale_x_continuous(breaks = enc_by_horizon$h) +
  labs(title = "GARCH Encompassing t-Statistic by Forecast Horizon",
       subtitle = "Dashed lines = 5% significance thresholds",
       x = "Horizon (days)", y = "t-statistic (GARCH coefficient)") +
  theme_bw(base_size = 11)

save_png(fig5, "figures/fig5_garch_tstat_horizon.png", w = 8, h = 5)
