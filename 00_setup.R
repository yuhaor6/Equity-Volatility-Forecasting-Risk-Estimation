pkgs <- c(
  "rugarch", "rmgarch", "vars", "forecast",
  "xts", "zoo", "quantmod", "PerformanceAnalytics",
  "lmtest", "car", "sgt", "ggplot2", "patchwork",
  "dplyr", "tidyr", "tseries"
)

missing <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(missing))
  install.packages(missing, repos = "https://cloud.r-project.org")

invisible(lapply(pkgs, library, character.only = TRUE))
