# run all scripts in order
# open in RStudio and press Ctrl+Shift+S, or:
#   setwd("C:/Users/yuhao/Documents/Equity Volatility Forecasting & Risk Estimation")
#   source("main.R")

if (interactive())
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))

source("00_setup.R")   # packages
source("01_data.R")    # data download (~1 min)
source("02_ewma.R")    # EWMA forecasts
source("03_garch.R")   # rolling GARCH (~20-40 min)
source("04_evaluation.R")
source("05_var_backtesting.R")
source("06_vol_risk_premium.R")
source("07_outputs.R")
