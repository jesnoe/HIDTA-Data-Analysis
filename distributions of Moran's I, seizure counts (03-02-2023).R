# setwd("/Users/euseongjang/Documents/R")
# setwd("C:/Users/gkfrj/Documents/R")
library(fpp2)
library(spdep)
library(scales)
library(readxl)
library(urbnmapr)
library(tidyverse)
library(gridExtra)
library(lubridate)

crack.rel <- read.csv("crack counts KNN5 R codes 999 two-sided relative (02-21-2023).csv") %>% as_tibble
cocaine.rel <- read.csv("cocaine other counts KNN5 R codes 999 two-sided relative (02-21-2023).csv") %>% as_tibble
crack.abs.1 <- read.csv("crack counts KNN5 R codes 999 two-sided absolute base_1 (02-21-2023).csv") %>% as_tibble
crack.abs.t_1 <- read.csv("crack counts KNN5 R codes 999 two-sided absolute base_t1 (02-21-2023).csv") %>% as_tibble

crack.rel %>% ggplot(aes(LISA_IJan8)) +
  geom_histogram(bins=100) +
  labs(x="Moran's I", title="Dist. of Moran's I in Jan 2020")
crack.rel %>% ggplot(aes(LISA_IJan8)) +
  geom_histogram(bins=100) +
  labs(x="Moran's I", title="Dist. of Moran's I in Jan 2020") +
  scale_y_continuous(limits=c(0,10), oob = rescale_none)

crack.rel %>% ggplot(aes(Jan_2020)) +
  geom_density() +
  labs(x="Seizure Count", title="Dist. of Seizure Counts I in Jan 2020")
crack.rel %>% ggplot(aes(Jan_2020)) +
  geom_histogram(bins=100) +
  labs(x="Seizure Count", title="Dist. of Seizure Counts I in Jan 2020") +
  scale_y_continuous(limits=c(0,10), oob = rescale_none)

crack.rel %>% ggplot(aes(Feb_2020)) +
  geom_histogram(bins=100)
crack.rel %>% ggplot(aes(Feb_2020)) +
  geom_histogram(bins=100) +
  scale_y_continuous(limits=c(0,10), oob = rescale_none)

LISA_C_to_Month <- function(col_names) {
  col_names <- gsub("0", "2020", col_names)
  col_names <- gsub("1", "2021", col_names)
  col_names <- gsub("8", "2018", col_names)
  col_names <- gsub("9", "2019", col_names)
  col_names <- parse_date(col_names, "LISA_C%b%Y") %>% substr(1,7)
  return(col_names)
}

LISA_I_to_Month <- function(col_names) {
  col_names <- gsub("0", "2020", col_names)
  col_names <- gsub("1", "2021", col_names)
  col_names <- gsub("8", "2018", col_names)
  col_names <- gsub("9", "2019", col_names)
  col_names <- parse_date(col_names, "LISA_I%b%Y") %>% substr(1,7)
  return(col_names)
}

# names(crack.rel)[5:52] <- LISA_C_to_Month(names(crack.rel)[5:52])
names(crack.rel)[5:52] <- parse_date(names(crack.rel)[5:52], "%b_%Y") %>% substr(1,7)

crack.rel.I <- crack.rel %>% select(state:HIDTA, contains("LISA_I"))
names(crack.rel.I)[5:52] <- LISA_I_to_Month(names(crack.rel.I)[5:52])

crack.rel %>%
  select(state:HIDTA, `2020-01`:`2020-12`) %>% 
  pivot_longer(c(-state, -county, -GEOID, -HIDTA), names_to="year_month", values_to="seizure_counts") %>% 
  ggplot(aes(seizure_counts)) +
  facet_wrap(. ~ year_month) +
  geom_histogram(bins=100)

crack.rel.I %>%
  select(state:HIDTA, `2020-01`:`2020-12`) %>% 
  pivot_longer(c(-state, -county, -GEOID, -HIDTA), names_to="year_month", values_to="LISA_I") %>% 
  ggplot(aes(LISA_I)) +
  facet_wrap(. ~ year_month) +
  geom_histogram(bins=100)

crack.rel %>%
  select(state:HIDTA, `2020-01`:`2020-12`) %>% 
  pivot_longer(c(-state, -county, -GEOID, -HIDTA), names_to="year_month", values_to="seizure_counts") %>% 
  ggplot(aes(seizure_counts)) +
  facet_wrap(. ~ year_month) +
  geom_histogram(bins=100) +
  scale_y_continuous(limits=c(0,10), oob = rescale_none)

crack.rel.I %>%
  select(state:HIDTA, `2020-01`:`2020-12`) %>% 
  pivot_longer(c(-state, -county, -GEOID, -HIDTA), names_to="year_month", values_to="LISA_I") %>% 
  ggplot(aes(LISA_I)) +
  facet_wrap(. ~ year_month) +
  geom_histogram(bins=100) +
  scale_y_continuous(limits=c(0,10), oob = rescale_none)

crack.rel$`2018-01` %>% table
