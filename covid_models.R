library(tidyverse)
library(lubridate)
library(forecast)
library(tseries)
library(scales)

df <- read.csv("C:/Users/visma/Downloads/time_series_covid19_confirmed_global.csv")
meta_cols <- names(df)[1:4]
raw_date_cols <- names(df)[5:ncol(df)]
cleaned_dates <- gsub("\\.", "/", gsub("^X", "", raw_date_cols))
names(df)[5:ncol(df)] <- make.unique(map_chr(cleaned_dates, ~ as.character(mdy(.x))))

df_long <- df %>%
  pivot_longer(cols = -all_of(meta_cols), names_to = "Date", values_to = "Confirmed") %>%
  mutate(Date = as.Date(Date))

# globally total confirmed cases per day
global_cases <- df_long %>%
  group_by(Date) %>%
  summarise(Total_Confirmed = sum(Confirmed)) %>%
  arrange(Date)

# Convert to monthly new cases
monthly_cases <- global_cases %>%
  mutate(Month = floor_date(Date, "month")) %>%
  group_by(Month) %>%
  summarise(Monthly_Cases = max(Total_Confirmed, na.rm = TRUE)) %>%
  mutate(New_Cases = Monthly_Cases - lag(Monthly_Cases)) %>%
  drop_na()

# Convert to time series
ts_covid <- ts(monthly_cases$New_Cases,
               frequency = 12,
               start = c(year(min(monthly_cases$Month)), month(min(monthly_cases$Month))))

autoplot(ts_covid) +
  ggtitle("Global Monthly New COVID-19 Cases") +
  xlab("Time") + ylab("New Cases")


# ADF Test for stationarity
adf_result <- adf.test(ts_covid)
print(adf_result)

if(adf_result$p.value > 0.05){
  ts_log <- log(ts_covid + 1)
  ts_diff <- diff(ts_log)
  autoplot(ts_diff) +
    ggtitle("Differenced Log-Transformed COVID-19 Cases") +
    xlab("Time") + ylab("Diff(log(Cases))")
  ts_preprocessed <- ts_diff
} else {
  ts_preprocessed <- ts_covid
}

# ACF and PACF for Model Diagnostics
acf(ts_preprocessed, main = "ACF - COVID-19 Monthly New Cases")
pacf(ts_preprocessed, main = "PACF - COVID-19 Monthly New Cases")

# Fit models
fit_ar <- arima(ts_preprocessed, order = c(3, 0, 0))
fit_ma <- arima(ts_preprocessed, order = c(0, 0, 3))
fit_arma <- arima(ts_preprocessed, order = c(2, 0, 2))
fit_arima <- auto.arima(ts_covid)

# Evaluate with AIC/BIC
model_list <- list(AR = fit_ar, MA = fit_ma, ARMA = fit_arma, ARIMA = fit_arima)
results <- data.frame(Model = names(model_list),
                      AIC = sapply(model_list, AIC),
                      BIC = sapply(model_list, BIC))
print(results)

# Assign the best model (if it's lowest AIC)
best_model <- fit_ma  # Change this if another model scores better

# Forecast next 6 months
fc <- forecast(best_model, h = 6)

autoplot(fc) +
  ggtitle("6-Month Forecast of Global Monthly COVID-19 Cases") +
  xlab("Time") + ylab("Forecasted Cases")


# Extract residuals
resid_covid <-checkresiduals(fit_ma)


# Residuals over time
autoplot(resid_covid) +
  ggtitle("Residuals over Time - COVID-19 ARIMA(0,0,3)") +
  xlab("Time") + ylab("Residuals")

# Histogram and QQ plot
ggplot(data.frame(resid = resid_covid), aes(x = resid)) +
  geom_histogram(bins = 15, fill = "skyblue", color = "black") +
  ggtitle("Histogram of Residuals")

qqnorm(resid_covid)
qqline(resid_covid, col = "red")

# ACF of residuals
acf(resid_covid, main = "ACF of Residuals - COVID-19")

# Ljung-Box Test (to check white noise)
Box.test(resid_covid, lag = 10, type = "Ljung-Box")

# Shapiro-Wilk Normality Test
shapiro.test(resid_covid)



