# Load the required libraries
library(tidyverse)
library(lubridate)
library(scales)
library(forecast)
library(tseries)
library(dplyr)

# Load and prepare the Influenza dataset
flu_df <- read.csv("time series project/VIW_FNT.csv", stringsAsFactors = FALSE)
flu_df$ISO_WEEKSTARTDATE <- ymd(flu_df$ISO_WEEKSTARTDATE)
flu_df <- arrange(flu_df, ISO_WEEKSTARTDATE)

# EDA: Data structure and summary
str(flu_df)
summary(flu_df)
range(flu_df$ISO_WEEKSTARTDATE, na.rm = TRUE)
unique(flu_df$WHOREGION)
unique(flu_df$HEMISPHERE)
colSums(is.na(flu_df)) %>% sort(decreasing = TRUE) %>% head(10)

# Top 10 countries by total flu cases
top_countries <- flu_df %>%
  group_by(COUNTRY_AREA_TERRITORY) %>%
  summarise(Total_Flu_Cases = sum(INF_ALL, na.rm = TRUE)) %>%
  arrange(desc(Total_Flu_Cases)) %>%
  slice_head(n = 10)

ggplot(top_countries, aes(x = reorder(COUNTRY_AREA_TERRITORY, Total_Flu_Cases),
                          y = Total_Flu_Cases,
                          fill = COUNTRY_AREA_TERRITORY)) +
  geom_col() +
  coord_flip() +
  labs(title = "Top 10 Countries by Total Influenza Cases",
       x = "Country",
       y = "Total Flu Cases") +
  scale_y_continuous(labels = comma) +
  theme_minimal()

# Weekly global flu cases
flu_global <- flu_df %>%
  group_by(ISO_WEEKSTARTDATE) %>%
  summarise(Total_Flu_Cases = sum(INF_ALL, na.rm = TRUE))

ggplot(flu_global, aes(x = ISO_WEEKSTARTDATE, y = Total_Flu_Cases)) +
  geom_line(color = "blue") +
  labs(title = "Global Weekly Influenza Cases",
       x = "Week", y = "Cases") +
  scale_y_continuous(labels = comma) +
  theme_minimal()

# Weekly flu cases in one country (e.g., China)
country_name <- "China"
country_cases <- flu_df %>%
  filter(COUNTRY_AREA_TERRITORY == country_name) %>%
  group_by(ISO_WEEKSTARTDATE) %>%
  summarise(Total_Flu_Cases = sum(INF_ALL, na.rm = TRUE))

ggplot(country_cases, aes(x = ISO_WEEKSTARTDATE, y = Total_Flu_Cases)) +
  geom_line(color = "darkgreen") +
  labs(title = paste("Weekly Influenza Cases in", country_name),
       x = "Week", y = "Influenza Cases") +
  scale_y_continuous(labels = comma) +
  theme_minimal()

# Outlier detection with Z-scores
mean_cases <- mean(flu_global$Total_Flu_Cases, na.rm = TRUE)
sd_cases <- sd(flu_global$Total_Flu_Cases, na.rm = TRUE)

flu_global <- flu_global %>%
  mutate(z_score = (Total_Flu_Cases - mean_cases) / sd_cases,
         outlier = abs(z_score) > 2)

ggplot(flu_global, aes(x = ISO_WEEKSTARTDATE, y = Total_Flu_Cases)) +
  geom_line(color = "blue") +
  geom_point(aes(color = outlier), size = 2) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  labs(title = "Global Weekly Influenza Cases (Outliers Highlighted)",
       x = "Week", y = "Cases", color = "Outlier") +
  scale_y_continuous(labels = comma) +
  theme_minimal()

print(filter(flu_global, outlier))

# Combine weekly global influenza cases and the processed specimen covariate
flu_model_data <- flu_df %>%
  group_by(ISO_WEEKSTARTDATE) %>%
  summarise(
    total_flu_cases = sum(INF_ALL, na.rm = TRUE),
    specimens_processed = sum(SPEC_PROCESSED_NB, na.rm = TRUE)
  ) %>%
  arrange(ISO_WEEKSTARTDATE)

# Convert to time series
ts_flu <- ts(flu_model_data$total_flu_cases,
             frequency = 52,
             start = c(year(min(flu_model_data$ISO_WEEKSTARTDATE)),
                       week(min(flu_model_data$ISO_WEEKSTARTDATE))))
# Covariate in time series
ts_covariate <- ts(flu_model_data$specimens_processed,
                   start = start(ts_flu),
                   frequency = frequency(ts_flu))

# ADF test for checking stationarity
adf_result <- adf.test(ts_flu)
print(adf_result)

# Transformation if required
if (adf_result$p.value > 0.05) {
  ts_flu_preprocessed <- diff(log(ts_flu + 1))
  covariate_specimens <- ts_covariate[-1]  # remove first to align with diff()
  message("Applied logarithm-differencing for stationarity.")
} else {
  ts_flu_preprocessed <- ts_flu
  covariate_specimens <- ts_covariate  # same length, no diff needed
  message("Data is already stationary.")
}

# ACF and PACF plots
acf(ts_flu_preprocessed, main = "ACF - Weekly Global Influenza Cases with Processed Specimen Covariate")
pacf(ts_flu_preprocessed, main = "PACF - Weekly Global Influenza Cases with Processed Specimen Covariate")

# Fit ARIMA without covariate data
fit_no_xreg <- auto.arima(ts_flu_preprocessed)

# Fit ARIMA with covariate data (SPEC_PROCESSED_NB) using same preprocessed data
fit_xreg <- auto.arima(ts_flu_preprocessed, xreg = ts_covariate) # The xreg function is used for external covariates to improve the accuracy of the forecast
summary(fit_xreg)

# Evaluate using AIC/BIC
comparison <- data.frame(
  Model = c("ARIMA with No Covariate", "ARIMA with Covariate"),
  AIC = c(AIC(fit_no_xreg), AIC(fit_xreg)),
  BIC = c(BIC(fit_no_xreg), BIC(fit_xreg))
)
print(comparison)

# Forecast the next 12 weeks using historical mean of the covariate
future_covariate <- rep(mean(ts_covariate, na.rm = TRUE), 12) # Assuming the covariate stays constant
forecast_xreg <- forecast(fit_xreg, xreg = future_covariate, h = 12)

# Forecast plot
autoplot(forecast_xreg) +
  ggtitle("12-Week Influenza Forecast with Processed Specimen Covariate") +
  xlab("Week") + ylab("Forecasted Influenza Cases (log-diff scale)")

# Residual diagnostics
resid_flu <- residuals(fit_xreg)

# Residuals over time plot
autoplot(resid_flu) +
  ggtitle("Residuals Over Time - ARIMA with Covariate") +
  xlab("Time") + ylab("Residuals")

# Histogram and density of residuals
ggplot(data.frame(resid = resid_flu), aes(x = resid)) +
  geom_histogram(aes(y = ..density..), bins = 15, fill = "lightcoral", color = "black") +
  geom_density(color = "blue") +
  ggtitle("Histogram and Density of Residuals") # Density is used for better visualization

# QQ plot of residuals for normality
qqnorm(resid_flu)
qqline(resid_flu, col = "red")

# ACF plot of residuals (check for remaining autocorrelation)
acf(resid_flu, main = "ACF Plot of Residuals for Influenza")

# Ljung-Box Test for residual autocorrelation
lb_flu <- Box.test(resid_flu, lag=10, type="Ljung-Box")
print(lb_flu)

# Shapiro-Wilk Test for residual normality
sw_flu <- shapiro.test(resid_flu)
print(sw_flu)

# Forecast accuracy (RMSE) using train-test split
n <- length(ts_flu_preprocessed)
train_size <- floor(0.8 * n)

# Split data into training and testing sets to help calculate RMSE
ts_train <- ts_flu_preprocessed[1:train_size] # Data without the covariate
ts_test <- ts_flu_preprocessed[(train_size + 1):n]

xreg_train <- covariate_specimens[1:train_size] # Data with the covariate
xreg_test <- covariate_specimens[(train_size + 1):n]

# Fit model on training set
fit_train <- auto.arima(ts_train, xreg = xreg_train)

# Forecast on test set
h <- length(ts_test)
fc_test <- forecast(fit_train, xreg = xreg_test, h = h)

# Compute RMSE (still on log-diff scale)
rmse_flu <- sqrt(mean((fc_test$mean - ts_test)^2, na.rm = TRUE))
cat("Influenza Model RMSE (Transformed):", round(rmse_flu, 2), "\n")

# Get the last actual value before the test period
last_actual_val <- ts_flu[train_size]

# Reconstruct and back-transform the forecast and actual values in original scale
forecast_original <- last_actual_val * exp(cumsum(fc_test$mean))
actual_original <- ts_flu[(train_size + 1):n]

# Sanity checks
# Check for negative values
any(forecast_original < 0)
any(actual_original < 0)

# Check lengths match
length(forecast_original) == length(actual_original)

# Create proper date sequence
date_seq <- seq.Date(
  as.Date(date_decimal(time(ts_flu)[train_size + 1])),
  by = "week",
  length.out = length(actual_original)
)

# Create plot data frame
df_plot <- data.frame(
  Date = date_seq,
  Actual = actual_original,
  Forecast = forecast_original
)

# Forecast vs actual plot (in original scale)
ggplot(plot_data, aes(x = Date)) +
  geom_line(aes(y = Actual, color = "Actual"), size = 1) +
  geom_line(aes(y = Forecast, color = "Forecast"), linetype = "dashed", size = 1) +
  scale_color_manual(name = "", values = c("Actual" = "darkgreen", "Forecast" = "red")) +
  labs(title = "Influenza Cases: Forecast vs Actual", y = "Weekly Flu Cases (Back-Transformed)", x = "Date") +
  theme_minimal() +
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_y_continuous(labels = scales::comma) # Format y-axis numbers

# Baseline model for comparison -----
# Fit model without covariate
fit_no_covariate <- auto.arima(ts_train)

# Forecast without covariate
fc_no_covariate <- forecast(fit_no_covariate, h = length(ts_test))

# Compute RMSE without covariate
rmse_no_covariate <- sqrt(mean((fc_no_covariate$mean - ts_test)^2, na.rm = TRUE))
cat("Influenza Model RMSE without Covariate:", round(rmse_no_covariate, 2), "\n")

# Compare to the one of with covariate
cat("Influenza Model RMSE with Covariate:", round(rmse_flu, 2), "\n")

# Residual diagnostics for no-covariate model
resid_nocov <- residuals(fit_no_covariate)

# Residuals over time plot without covariate
autoplot(resid_nocov) +
  ggtitle("Residuals Over Time - ARIMA without Covariate") +
  xlab("Time") + ylab("Residuals")

# Histogram and density of residuals without covariate
ggplot(data.frame(resid = resid_nocov), aes(x = resid)) +
  geom_histogram(aes(y = ..density..), bins = 15, fill = "lightblue", color = "black") +
  geom_density(color = "blue") +
  ggtitle("Histogram and Density of Residuals without Covariate")

# QQ plot of residuals without covariate
qqnorm(resid_nocov)
qqline(resid_nocov, col = "red")

# ACF plot of residuals without covariate
acf(resid_nocov, main = "ACF Plot of Residuals without Covariate")

# Ljung-Box Test without covariate
lb_nocov <- Box.test(resid_nocov, lag=10, type="Ljung-Box")
print(lb_nocov)

# Shapiro-Wilk Test without covariate
sw_nocov <- shapiro.test(resid_nocov)
print(sw_nocov)