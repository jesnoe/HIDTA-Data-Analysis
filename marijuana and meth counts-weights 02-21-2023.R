library(tidyverse)
library(fpp2)
library(spdep)
library(readxl)
library(urbnmapr)
library(gridExtra)

overdose <- read.csv("VSRR_Provisional_County-Level_Drug_Overdose_Death_Counts.csv") %>% as_tibble

# urbnmapr data
unique(countydata$county_fips) # some ips have 5 digits. need to be converted into integer.
HIDTA.dist <- read.csv("HIDTA Regions.csv", header=T) %>% as_tibble
coordinate.HIDTA <- left_join(counties, HIDTA.dist, by=c("state_name", "county_name")) %>% filter(!(state_name %in% c("Alaska", "Hawaii", "Puerto Rico")))
names(coordinate.HIDTA)[7] <- "GEOID"
coordinate.HIDTA$county <- str_split(coordinate.HIDTA$county_name, " ") %>%
  lapply(function(x) str_c(c(x[(1:(length(x)-1))]), collapse=" ") ) %>%
  unlist
HIDTA.regions <- coordinate.HIDTA %>% filter(!is.na(HIDTA)) %>% select(GEOID, HIDTA) %>% unique

seizures <- read_xlsx("Drug Seizures All HIDTAs All Drugs 2018-2021 Combined.xlsx") %>% filter(State %in% coordinate.HIDTA$state_name)
seizures <- seizures %>% filter(County != "Michigan")
seizures$County <- substring(seizures$County, 1, str_locate(seizures$County, ",")[,1]-1)
# seizures$County <- sapply(seizures$County, function(x) paste(x[1], str_to_lower(substring(x, 2, length(x)-1)), sep=""))
seizures$Year <- substring(seizures$SeizureDate, 1, 4)
seizures$Month <- substring(seizures$SeizureDate, 6, 7)
seizures$Year <- as.numeric(seizures$Year)
seizures$Month <- as.numeric(seizures$Month)
seizures$County[seizures$State == "District of Columbia"] <- "District of"

names(seizures)[1:2] <- c("state_name", "county")

all_counts <- seizures %>% left_join(unique(coordinate.HIDTA[,c(7, 12, 14, 15)]), by=c("state_name", "county")) %>% 
  group_by(state_name, county, GEOID, Year, Month) %>%
  summarise(count=n()) %>% 
  ungroup %>%
  select(GEOID:count) %>%
  complete(GEOID, Year, Month)

all_counts <- all_counts %>% mutate(Month2=month(Month, label=T))
all_counts_missingness <- all_counts %>% select(-Month) %>% 
  pivot_wider(names_from=c(Month2, Year), values_from=count) %>% 
  arrange(GEOID)
all_counts_missingness$missing <- all_counts_missingness %>%
  apply(1, function(x) ifelse(sum(x, na.rm=T) == 0, 1, 0)) # 1 if missing (no observation for 48 months)
non_missing_GEOID <- all_counts_missingness$GEOID[all_counts_missingness$missing == 0]


# meth count
{
  meth <- seizures %>% filter(Drug %in% c("Methamphetamine", "Ice", "Amphetamine", "Methamphetamine in solution",
                                          "Meth precursor: Ephedrine hydrochloride", "Methamphetamine Oil", "Meth Precursor Chemicals",
                                          "Methamphetamine Tablet", "Ephedrine", "Pseudoephedrine"), Unit == "Kg")
  meth <- meth %>% mutate(Month2=month(Month, label=T))
  meth <- meth %>% group_by(state_name, county, Year, Month2) %>%
    summarise(count=sum(Quantity > 0)) %>%
    complete(state_name, county, Year, Month2) %>% 
    pivot_wider(names_from=c(Month2, Year), names_sep="_", values_from=count)
  cal_order <- c(str_c(calender$Month2, 2018, sep="_"),
                 str_c(calender$Month2, 2019, sep="_"),
                 str_c(calender$Month2, 2020, sep="_"),
                 str_c(calender$Month2, 2021, sep="_"))
  meth <- meth[, c(1:2, match(cal_order, names(meth)))]
  
  meth <- left_join(unique(coordinate.HIDTA[,c(7, 12, 14, 15)]), meth, by=c("state_name", "county")) %>%
    select(state_name, county, GEOID, HIDTA, Jan_2018:Dec_2021) %>% arrange(GEOID)
  names(meth)[1] <- "state"
}

meth %>% filter(is.na(county))
grep("Jan_2018", names(meth)) # 5
grep("Dec_2021", names(meth)) # 52
meth <- meth %>% filter(GEOID %in% non_missing_GEOID | !is.na(HIDTA))
meth[,5:52][is.na(meth[,5:52])] <- 0
meth %>% filter(is.na(GEOID))
meth # 1518 counties in total.
# write.csv(meth, "meth count HIDTA (02-21-2023).csv", row.names=F)

# meth weight
{
  meth <- seizures %>% filter(Drug %in% c("Methamphetamine", "Ice", "Amphetamine", "Methamphetamine in solution",
                                          "Meth precursor: Ephedrine hydrochloride", "Methamphetamine Oil", "Meth Precursor Chemicals",
                                          "Methamphetamine Tablet", "Ephedrine", "Pseudoephedrine"), Unit == "Kg")
  meth <- meth %>% mutate(Month2=month(Month, label=T))
  meth <- meth %>% group_by(state_name, county, Year, Month2) %>%
    summarise(weight=sum(Quantity)) %>%
    complete(state_name, county, Year, Month2) %>% 
    pivot_wider(names_from=c(Month2, Year), names_sep="_", values_from=weight)
  meth <- meth[, c(1:2, match(cal_order, names(meth)))]
  
  meth <- left_join(unique(coordinate.HIDTA[,c(7, 12, 14, 15)]), meth, by=c("state_name", "county")) %>%
    select(state_name, county, GEOID, HIDTA, Jan_2018:Dec_2021) %>% arrange(GEOID)
  names(meth)[1] <- "state"
}

meth %>% filter(is.na(county))
grep("Jan_2018", names(meth)) # 5
grep("Dec_2021", names(meth)) # 52
meth <- meth %>% filter(GEOID %in% non_missing_GEOID | !is.na(HIDTA))
meth[,5:52][is.na(meth[,5:52])] <- 0
meth %>% filter(is.na(GEOID))
meth # 1518 counties in total.
# write.csv(meth, "meth weight HIDTA (02-21-2023).csv", row.names=F)

# marijuana count
{
  marijuana <- seizures %>% filter(grepl("Marijuana", Drug), Unit == "Kg")
  marijuana <- marijuana %>% mutate(Month2=month(Month, label=T))
  marijuana <- marijuana %>% group_by(state_name, county, Year, Month2) %>%
    summarise(count=sum(Quantity > 0)) %>%
    complete(state_name, county, Year, Month2) %>% 
    pivot_wider(names_from=c(Month2, Year), names_sep="_", values_from=count)
  cal_order <- c(str_c(calender$Month2, 2018, sep="_"),
                 str_c(calender$Month2, 2019, sep="_"),
                 str_c(calender$Month2, 2020, sep="_"),
                 str_c(calender$Month2, 2021, sep="_"))
  marijuana <- marijuana[, c(1:2, match(cal_order, names(marijuana)))]
  
  marijuana <- left_join(unique(coordinate.HIDTA[,c(7, 12, 14, 15)]), marijuana, by=c("state_name", "county")) %>%
    select(state_name, county, GEOID, HIDTA, Jan_2018:Dec_2021) %>% arrange(GEOID)
  names(marijuana)[1] <- "state"
}

marijuana %>% filter(is.na(county))
grep("Jan_2018", names(marijuana)) # 5
grep("Dec_2021", names(marijuana)) # 52
marijuana <- marijuana %>% filter(GEOID %in% non_missing_GEOID | !is.na(HIDTA))
marijuana[,5:52][is.na(marijuana[,5:52])] <- 0
marijuana %>% filter(is.na(GEOID))
marijuana # 1518 counties in total.
# write.csv(marijuana, "marijuana count HIDTA (02-21-2023).csv", row.names=F)

# marijuana weight
{
  marijuana <- seizures %>% filter(grepl("Marijuana", Drug), Unit == "Kg")
  marijuana <- marijuana %>% mutate(Month2=month(Month, label=T))
  marijuana <- marijuana %>% group_by(state_name, county, Year, Month2) %>%
    summarise(weight=sum(Quantity)) %>%
    complete(state_name, county, Year, Month2) %>% 
    pivot_wider(names_from=c(Month2, Year), names_sep="_", values_from=weight)
  marijuana <- marijuana[, c(1:2, match(cal_order, names(marijuana)))]
  
  marijuana <- left_join(unique(coordinate.HIDTA[,c(7, 12, 14, 15)]), marijuana, by=c("state_name", "county")) %>%
    select(state_name, county, GEOID, HIDTA, Jan_2018:Dec_2021) %>% arrange(GEOID)
  names(marijuana)[1] <- "state"
}

marijuana %>% filter(is.na(county))
grep("Jan_2018", names(marijuana)) # 5
grep("Dec_2021", names(marijuana)) # 52
marijuana <- marijuana %>% filter(GEOID %in% non_missing_GEOID | !is.na(HIDTA))
marijuana[,5:52][is.na(marijuana[,5:52])] <- 0
marijuana %>% filter(is.na(GEOID))
marijuana # 1518 counties in total.
# write.csv(marijuana, "marijuana weight HIDTA (02-21-2023).csv", row.names=F)