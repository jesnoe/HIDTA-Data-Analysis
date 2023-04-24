library(readxl)
library(urbnmapr)
library(tidyverse)
library(gridExtra)
library(lubridate)

# urbnmapr data
unique(countydata$county_fips) # some ips have 5 digits. need to be converted into integer.
HIDTA.dist <- read.csv("HIDTA Regions.csv", header=T) %>% as_tibble
coordinate.HIDTA <- left_join(counties, HIDTA.dist, by=c("state_name", "county_name")) %>% filter(!(state_name %in% c("Alaska", "Hawaii", "Puerto Rico")))
names(coordinate.HIDTA)[7] <- "GEOID"
coordinate.HIDTA$county <- str_split(coordinate.HIDTA$county_name, " ") %>%
  lapply(function(x) str_c(c(x[(1:(length(x)-1))]), collapse=" ") ) %>%
  unlist
HIDTA.regions <- coordinate.HIDTA %>% filter(!is.na(HIDTA)) %>% select(GEOID, HIDTA) %>% unique

non_US_territory <- read_xlsx("Drug Seizures All HIDTAs All Drugs 2018-2021 Original.xlsx") %>% 
  filter(!(State %in% coordinate.HIDTA$state_name)) %>% 
  pull(State) %>% unique
non_US_territory

seizures <- read_xlsx("Drug Seizures All HIDTAs All Drugs 2018-2021 Original.xlsx") %>% filter(State %in% coordinate.HIDTA$state_name)
# seizures <- seizures %>% filter(County != "Michigan")
seizures$County <- substring(seizures$County, 1, str_locate(seizures$County, ",")[,1]-1)
seizures$Year <- substring(seizures$SeizureDate, 1, 4)
seizures$Month <- substring(seizures$SeizureDate, 6, 7)
seizures$Year <- as.numeric(seizures$Year)
seizures$Month <- as.numeric(seizures$Month)
names(seizures)[1:2] <- c("state_name", "county")
all_counts <- seizures %>% 
  group_by(Drug) %>%
  summarise(count=sum(Quantity > 0)) %>% 
  arrange(desc(count))
all_counts
seizures %>% 
  group_by(Drug) %>%
  summarise(count=n()) %>% 
  arrange(desc(count))

seizures_renamed <- read_xlsx("Drug Seizures All HIDTAs All Drugs 2018-2021 Combined.xlsx") %>% filter(State %in% coordinate.HIDTA$state_name)
# seizures_renamed <- seizures_renamed %>% filter(County != "Michigan")
seizures_renamed$County <- substring(seizures_renamed$County, 1, str_locate(seizures_renamed$County, ",")[,1]-1)
seizures_renamed$Year <- substring(seizures_renamed$SeizureDate, 1, 4)
seizures_renamed$Month <- substring(seizures_renamed$SeizureDate, 6, 7)
seizures_renamed$Year <- as.numeric(seizures_renamed$Year)
seizures_renamed$Month <- as.numeric(seizures_renamed$Month)
names(seizures_renamed)[1:2] <- c("state_name", "county")
all_counts_renamed <- seizures_renamed %>% 
  group_by(Drug) %>%
  summarise(count=sum(Quantity > 0)) %>% 
  arrange(desc(count))
all_counts_renamed
seizures_renamed %>% 
  group_by(Drug) %>%
  summarise(count=n()) %>% 
  arrange(desc(count))

seizures %>% filter(Quantity == 0)
seizures_renamed %>% filter(Quantity == 0)

seizures %>% filter(!(county %in% unique(seizures_renamed$county))) %>% select(state_name, county) %>% unique
seizures_renamed %>% filter(!(county %in% unique(seizures$county))) %>% select(state_name, county) %>% unique

seizures %>% select(state_name, county) %>% unique %>% filter(state_name=="Maryland")
seizures_renamed %>% select(state_name, county) %>% unique %>% filter(state_name=="Maryland")
