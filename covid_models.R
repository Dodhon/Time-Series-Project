# Load the required libraries
library(tidyverse)
library(lubridate)
library(forecast)
library(tseries)
library(scales)
library(plotly)

# Load and preprocess the COVID-19 dataset
df <- read.csv("time series project/time_series_covid19_confirmed_global.csv")

meta_cols <- names(df)[1:4]
raw_date_cols <- names(df)[5:ncol(df)]
cleaned_dates <- gsub("\\.", "/", gsub("^X", "", raw_date_cols))
names(df)[5:ncol(df)] <- make.unique(map_chr(cleaned_dates, ~ as.character(mdy(.x))))

df_long <- df %>%
  pivot_longer(cols = -all_of(meta_cols), names_to = "Date", values_to = "Confirmed") %>%
  mutate(Date = as.Date(Date))

# EDA: Total global confirmed cases per day
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
ts_covid <- ts(monthly_cases$New_Cases, frequency = 12,
               start = c(year(min(monthly_cases$Month)), month(min(monthly_cases$Month))))

ts_covis_plot <- autoplot(ts_covid) +
  ggtitle("Global Monthly New COVID-19 Cases") +
  xlab("Time") + ylab("New Cases")

# ADF test for checking stationarity
adf_result <- adf.test(ts_covid)
print(adf_result)

# Transformation if needed
if(adf_result$p.value > 0.05){
  ts_preprocessed <- diff(log(ts_covid + 1))
} else {
  ts_preprocessed <- ts_covid
}

# Check autocorrelation with ACF and PACF plots
acf(ts_preprocessed, main = "ACF Plot - COVID-19 Monthly New Cases")
pacf(ts_preprocessed, main = "PACF Plot - COVID-19 Monthly New Cases")

# Fit and compare multiple models including SARIMA
fit_ar <- arima(ts_preprocessed, order = c(3, 0, 0))
fit_ma <- arima(ts_preprocessed, order = c(0, 0, 3))
fit_arma <- arima(ts_preprocessed, order = c(2, 0, 2))
fit_arima <- auto.arima(ts_covid, seasonal = FALSE)
fit_sarima <- auto.arima(ts_covid, seasonal = TRUE)

# Evaluate using AIC/BIC
model_list <- list(AR = fit_ar, MA = fit_ma, ARMA = fit_arma, ARIMA = fit_arima, SARIMA = fit_sarima)
results <- data.frame(Model = names(model_list),
                      AIC = sapply(model_list, AIC),
                      BIC = sapply(model_list, BIC))
print(results)

# Automatically picks best model with lowest AIC
best_model_name <- results$Model[which.min(results$AIC)]
best_model <- model_list[[best_model_name]]
cat("Best model selected:", best_model_name, "\n")

# Forecasting with the selected model for the next 6 months
fc <- forecast(best_model, h = 6)

# Interactive Forecast Plot (Zoomable with plotly)
ggplotly(
  autoplot(fc) +
    labs(title = paste("Interactive Forecast Plot -", best_model_name, "Model"),
         x = "Time", y = "Forecasted Cases")
)

# Residual diagnostics
resid_covid <- residuals(best_model)

# Residuals over time plot
autoplot(resid_covid) +
  ggtitle(paste("Residuals Over Time -", best_model_name, "Model")) +
  xlab("Time") + ylab("Residuals")

# Histogram of residuals
ggplot(data.frame(resid = resid_covid), aes(x = resid)) +
  geom_histogram(bins = 15, fill = "skyblue", color = "black") +
  ggtitle("Histogram of Residuals")

# QQ plot of residuals for normality
qqnorm(resid_covid)
qqline(resid_covid, col = "red")

# ACF plot of residuals (check for remaining autocorrelation)
acf(resid_covid, main="ACF Plot of Residuals for COVID-19")

# Ljung-Box Test for residual autocorrelation
lb_covid <- Box.test(resid_covid, lag = 10, type = "Ljung-Box")
print(lb_covid)

# Shapiro-Wilk Test for residual normality
sw_covid <- shapiro.test(resid_covid)
print(sw_covid)

# Evaluate and enhance forecast accuracy using tsCV since no external covariates
  # Use only errors for the horizon
for (h in 1:6) {
  cv_errors <- tsCV(ts_covid, 
                    forecastfunction = function(y, h) {
                      fit <- auto.arima(y)
                      forecast(fit, h = h)
                    },
                    h = h)

  valid_errors <- cv_errors[!is.na(cv_errors)]
  cv_rmse <- if (length(valid_errors) > 0) sqrt(mean(valid_errors^2)) else NA # Remove NAs and calculate RMSE
  cat(paste("h =", h, "Cross-Validation RMSE:"), round(cv_rmse, 2), "\n")
}

# Inspecting the error matrix 
summary(cv_errors)  # See min, max, NA count
sum(is.na(cv_errors))  # If it's all NA, you'll get NaN for RMSE

# Plot cross-validation errors
autoplot(cv_errors) +
  ggtitle("Cross-Validation Errors Over Time") +
  xlab("Time") + ylab("Forecast Error")

# Cross-validation errors histogram (missing values handled)
cv_errors_clean <- na.omit(as.numeric(cv_errors))
ggplot(data.frame(errors = cv_errors_clean), aes(x = errors)) +
  geom_histogram(bins = 15, fill = "orange", color = "black") +
  ggtitle("Histogram of Cross-Validation Errors") +
  xlab("Forecast Errors") + 
  ylab("Frequency")
