library(tidyverse)
library(lubridate)
library(forecast)
library(tseries)
library(gridExtra)

# Helper: decimal year → "YYYY-Www"
dec_to_yearweek <- function(dec) {
  yr <- floor(dec)
  wk <- round((dec - yr) * 52) + 1
  sprintf("%d-W%02d", yr, wk)
}

flu_df <- read.csv("/Users/thuptenwangpo/Downloads/VIW_FNT.csv", stringsAsFactors = FALSE)
flu_df$ISO_WEEKSTARTDATE <- ymd(flu_df$ISO_WEEKSTARTDATE)
flu_df <- arrange(flu_df, ISO_WEEKSTARTDATE)

cutoffs <- c("2024-01-01", "2024-06-01",
             "2025-01-01", "2025-06-01")

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

# Display all plots in one column
grid.arrange(grobs = plot_list, ncol = 1)
