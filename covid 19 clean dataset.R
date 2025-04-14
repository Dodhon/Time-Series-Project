
## Covid-19 

library(tidyverse)
library(lubridate)
library(scales)
df<- read.csv("C:/Users/visma/Downloads/time_series_covid19_confirmed_global.csv")

str(df)


meta_cols <- names(df)[1:4]
raw_date_cols <- names(df)[5:ncol(df)]


cleaned_raw <- gsub("\\.", "/", gsub("^X", "", raw_date_cols)) 


cleaned_dates <- map_chr(cleaned_raw, ~ as.character(suppressWarnings(mdy(.x))))


sum(is.na(cleaned_dates))       
anyDuplicated(cleaned_dates)    


cleaned_dates_unique <- make.unique(cleaned_dates)
names(df)[5:ncol(df)] <- cleaned_dates_unique


df_long <- df %>%
  pivot_longer(cols = -all_of(meta_cols),
               names_to = "Date",
               values_to = "Confirmed") %>%
  mutate(Date = as.Date(Date))


head(df_long)


global_cases <- df_long %>%
  group_by(Date) %>%
  summarise(Total_Confirmed = sum(Confirmed))


ggplot(global_cases, aes(x = Date, y = Total_Confirmed)) +
  geom_line(color = "steelblue") +
  labs(title = "Global Cumulative Confirmed COVID-19 Cases",
       x = "Date", y = "Total Confirmed Cases") +
  scale_y_continuous(labels = comma) +
  theme_minimal()


# top 10 countries
latest_date <- max(df_long$Date)

top_countries <- df_long %>%
  filter(Date == latest_date) %>%
  group_by(`Country.Region`) %>%
  summarise(Total = sum(Confirmed)) %>%
  arrange(desc(Total)) %>%
  slice(1:10)

ggplot(top_countries, aes(x = reorder(`Country.Region`, Total), y = Total)) +
  geom_col(fill = "red") +
  coord_flip() +
  labs(title = "Top 10 Countries by Confirmed COVID-19 Cases",
       x = "Country", y = "Total Confirmed Cases") +
  scale_y_continuous(labels = comma) +
  theme_minimal()



# cases per country

country_name <- "US"

country_cases <- df_long %>%
  filter(`Country.Region` == country_name) %>%
  group_by(Date) %>%
  summarise(Confirmed = sum(Confirmed))

ggplot(country_cases, aes(x = Date, y = Confirmed)) +
  geom_line(color = "purple") +
  labs(title = paste("Confirmed Cases in", country_name),
       x = "Date", y = "Confirmed Cases") +
  scale_y_continuous(labels = comma)+
  theme_minimal()


