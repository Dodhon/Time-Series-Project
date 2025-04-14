library(tidyverse)
library(lubridate)
library(scales)

# only used countries, dates and inf_all(infected with both A and B types of influenza)


flu_df <- read.csv("C:/Users/visma/Downloads/VIW_FNT.csv")

str(flu_df)
summary(flu_df)


flu_df$ISO_WEEKSTARTDATE <- ymd(flu_df$ISO_WEEKSTARTDATE)


range(flu_df$ISO_WEEKSTARTDATE, na.rm = TRUE)


unique(flu_df$WHOREGION)
unique(flu_df$HEMISPHERE)

colSums(is.na(flu_df)) %>%
  sort(decreasing = TRUE) %>%
  head(10)



# top 10 countries infected
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
  scale_y_continuous(labels = comma)+
  theme_minimal() 



# global weekly influenza cases
flu_global <- flu_df %>%
  group_by(ISO_WEEKSTARTDATE) %>%
  summarise(Total_Flu_Cases = sum(INF_ALL, na.rm = TRUE))

ggplot(flu_global, aes(x = ISO_WEEKSTARTDATE, y = Total_Flu_Cases)) +
  geom_line(color = "blue") +
  labs(title = "Global Weekly Influenza Cases",
       x = "Week", y = "Cases") +
  scale_y_continuous(labels = comma)+
  theme_minimal()




# Choose the country
country_name <- "China"  


country_cases <- flu_df %>%
  filter(COUNTRY_AREA_TERRITORY == country_name) %>%
  group_by(ISO_WEEKSTARTDATE) %>%
  summarise(Total_Flu_Cases = sum(INF_ALL, na.rm = TRUE))


ggplot(country_cases, aes(x = ISO_WEEKSTARTDATE, y = Total_Flu_Cases)) +
  geom_line(color = "darkgreen") +
  labs(title = paste("Weekly Influenza Cases in", country_name),
       x = "Week",
       y = "Influenza Cases") +
  scale_y_continuous(labels = comma) +
  theme_minimal()

