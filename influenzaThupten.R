library(tidyverse)
library(lubridate)
library(scales)
library(dplyr)

flu_df <- read.csv("/Users/thuptenwangpo/Downloads/VIW_FNT.csv")
flu_df$ISO_WEEKSTARTDATE <- ymd(flu_df$ISO_WEEKSTARTDATE)
str(flu_df)
summary(flu_df)
range(flu_df$ISO_WEEKSTARTDATE, na.rm = TRUE)
unique(flu_df$WHOREGION)
unique(flu_df$HEMISPHERE)
colSums(is.na(flu_df)) %>% sort(decreasing = TRUE) %>% head(10)

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

flu_global <- flu_df %>%
  group_by(ISO_WEEKSTARTDATE) %>%
  summarise(Total_Flu_Cases = sum(INF_ALL, na.rm = TRUE))

ggplot(flu_global, aes(x = ISO_WEEKSTARTDATE, y = Total_Flu_Cases)) +
  geom_line(color = "blue") +
  labs(title = "Global Weekly Influenza Cases",
       x = "Week", y = "Cases") +
  scale_y_continuous(labels = comma) +
  theme_minimal()

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

# Outlier detection Z-scores
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

# Print outlier observations
print(filter(flu_global, outlier))
