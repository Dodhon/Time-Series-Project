library(tidyverse)
library(lubridate)
library(forecast)
library(tseries)

flu_df <- read.csv("/Users/thuptenwangpo/Downloads/VIW_FNT.csv", stringsAsFactors = FALSE)
flu_df$ISO_WEEKSTARTDATE <- ymd(flu_df$ISO_WEEKSTARTDATE)
flu_df <- arrange(flu_df, ISO_WEEKSTARTDATE)

flu_global <- flu_df %>% 
  group_by(ISO_WEEKSTARTDATE) %>% 
  summarise(Total_Flu_Cases = sum(INF_ALL, na.rm = TRUE)) %>% 
  ungroup()

ts_flu <- ts(flu_global$Total_Flu_Cases, frequency = 52, start = c(year(min(flu_global$ISO_WEEKSTARTDATE)), week(min(flu_global$ISO_WEEKSTARTDATE))))

autoplot(ts_flu) + ggtitle("Global Weekly Influenza Cases") + xlab("Time") + ylab("Cases")

adf_result <- adf.test(ts_flu)
print(adf_result)

if(adf_result$p.value > 0.05){
  ts_log <- log(ts_flu + 1) 
  ts_diff <- diff(ts_log)
  autoplot(ts_diff) + ggtitle("Differenced Log-Transformed Series") + xlab("Time") + ylab("Diff(log(Cases))")
  ts_preprocessed <- ts_diff
} else {
  ts_preprocessed <- ts_flu
}

# Fit models:
fit_ar <- arima(ts_preprocessed, order = c(3, 0, 0))
fit_ma <- arima(ts_preprocessed, order = c(0, 0, 3))
fit_arma <- arima(ts_preprocessed, order = c(2, 0, 2))
fit_arima <- auto.arima(ts_flu)


model_list <- list(AR = fit_ar, MA = fit_ma, ARMA = fit_arma, ARIMA = fit_arima)
results <- data.frame(Model = names(model_list),
                      AIC = sapply(model_list, AIC),
                      BIC = sapply(model_list, BIC))
print(results)

best_model <- fit_arima
fc <- forecast(best_model, h = 12)
autoplot(fc) + ggtitle("Forecast of Weekly Influenza Cases") + xlab("Time") + ylab("Cases")
