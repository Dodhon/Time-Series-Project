# Load libraries
library(tidyverse)
library(lubridate)
library(scales)
library(forecast)
library(tseries)
library(gridExtra)
library(dplyr)

# Load data
flu_df <- read.csv("C:/Users/visma/Downloads/VIW_FNT.csv", stringsAsFactors = FALSE)
flu_df$ISO_WEEKSTARTDATE <- ymd(flu_df$ISO_WEEKSTARTDATE)
flu_df <- arrange(flu_df, ISO_WEEKSTARTDATE)

# Data structure and summary
str(flu_df)
summary(flu_df)
range(flu_df$ISO_WEEKSTARTDATE, na.rm = TRUE)
unique(flu_df$WHOREGION)
unique(flu_df$HEMISPHERE)
colSums(is.na(flu_df)) %>% sort(decreasing = TRUE) %>% head(10)

# Top 10 countries by cases
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

# Global weekly influenza cases
flu_global <- flu_df %>%
  group_by(ISO_WEEKSTARTDATE) %>%
  summarise(Total_Flu_Cases = sum(INF_ALL, na.rm = TRUE))

ggplot(flu_global, aes(x = ISO_WEEKSTARTDATE, y = Total_Flu_Cases)) +
  geom_line(color = "blue") +
  labs(title = "Global Weekly Influenza Cases",
       x = "Week", y = "Cases") +
  scale_y_continuous(labels = comma) +
  theme_minimal()

# Cases in one country (e.g., China)
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

# Stationarity check for ARIMA
ts_flu <- ts(flu_global$Total_Flu_Cases,
             frequency = 52,
             start = c(year(min(flu_global$ISO_WEEKSTARTDATE)),
                       week(min(flu_global$ISO_WEEKSTARTDATE))))

adf_result <- adf.test(ts_flu)
print(adf_result)

if (adf_result$p.value > 0.05) {
  ts_log <- log(ts_flu + 1)
  ts_diff <- diff(ts_log)
  ts_preprocessed <- ts_diff
  message("Used log-differencing for stationarity.")
} else {
  ts_preprocessed <- ts_flu
  message("Data is already stationary.")
}

acf(ts_preprocessed, main = "ACF - Weekly Global Influenza Cases")
pacf(ts_preprocessed, main = "PACF - Weekly Global Influenza Cases")

# Forecasting across different cutoffs
dec_to_yearweek <- function(dec) {
  yr <- floor(dec)
  wk <- round((dec - yr) * 52) + 1
  sprintf("%d-W%02d", yr, wk)
}

cutoffs <- c("2024-01-01", "2024-06-01", "2025-01-01", "2025-06-01")
plot_list <- list()

for(cut_off in cutoffs) {
  df_sub <- flu_df %>% filter(ISO_WEEKSTARTDATE < ymd(cut_off))
  agg <- df_sub %>%
    group_by(ISO_WEEKSTARTDATE) %>%
    summarise(Total_Flu_Cases = sum(INF_ALL, na.rm = TRUE), .groups="drop")
  
  ts_flu <- ts(agg$Total_Flu_Cases,
               frequency = 52,
               start = c(year(min(agg$ISO_WEEKSTARTDATE)),
                         week(min(agg$ISO_WEEKSTARTDATE))))
  
  fit <- auto.arima(ts_flu)
  fc  <- forecast(fit, h = 12)
  
  n <- length(fc$mean)
  start_dec <- max(time(ts_flu)) + 1/52
  new_dec <- seq(from = start_dec, by = 1/52, length.out = n)
  
  df_fc <- tibble(
    decimal_time   = new_dec,
    Point_Forecast = as.numeric(fc$mean),
    Lo95           = fc$lower[,"95%"],
    Hi95           = fc$upper[,"95%"]
  ) %>%
    mutate(YearWeek = map_chr(decimal_time, dec_to_yearweek))
  
  p <- ggplot(df_fc, aes(x = YearWeek, y = Point_Forecast, group = 1)) +
    geom_ribbon(aes(ymin = Lo95, ymax = Hi95), fill = "gray80", alpha = .4) +
    geom_line(color = "steelblue", size = 1) +
    labs(title = paste("Forecast up to", cut_off),
         x = "Year-Week", y = "Weekly Flu Cases") +
    theme_minimal(base_size = 14) +
    theme(
      plot.margin = margin(1,1,1,1, "cm"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank()
    ) +
    scale_x_discrete(expand = expansion(mult = c(0.05, 0.05)))
  
  plot_list[[cut_off]] <- p
}

# Show all forecasts in a column
grid.arrange(grobs = plot_list, ncol = 1)



resid_flu <- residuals(fit_arima_flu)

# Time plot of residuals
autoplot(resid_flu) +
  ggtitle("Residuals over Time - Global Influenza ARIMA Model") +
  xlab("Time") + ylab("Residuals")

# Histogram
ggplot(data.frame(resid = resid_flu), aes(x = resid)) +
  geom_histogram(bins = 15, fill = "lightcoral", color = "black") +
  ggtitle("Histogram of Residuals - Flu")

# QQ plot
qqnorm(resid_flu)
qqline(resid_flu, col = "red")
  
# ACF
acf(resid_flu, main = "ACF of Residuals - Flu")

# Ljung-Box
Box.test(resid_flu, lag = 10, type = "Ljung-Box")x

# Shapiro-Wilk
shapiro.test(resid_flu)

