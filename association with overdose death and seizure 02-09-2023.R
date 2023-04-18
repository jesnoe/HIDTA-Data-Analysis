library(fpp2)
library(spdep)
library(readxl)
library(tidyverse)
library(gridExtra)
library(lubridate)
{
  adjacency <- read.csv("county_adjacency2010.csv")
  adjacency$county <- substring(adjacency$countyname, 1, str_locate(adjacency$countyname, ",")[,1]-1)
  adjacency$state1 <- substring(adjacency$countyname, str_locate(adjacency$countyname, ",")[,1]+2)
  adjacency$neighbor_county <- substring(adjacency$neighborname, 1, str_locate(adjacency$neighborname, ",")[,1]-1)
  adjacency$neighbor_state1 <- substring(adjacency$neighborname, str_locate(adjacency$neighborname, ",")[,1]+2)
  
  state_names <- read.csv("state names.csv")
  names(state_names)[1] <- "state"
  names(state_names)[3] <- "state1"
  adjacency <- merge(adjacency, state_names[,c(1,3)], all.x=T)
  
  names(state_names)[1] <- "neighbor_state"
  names(state_names)[3] <- "neighbor_state1"
  adjacency <- merge(adjacency, state_names[,c(1,3)], all.x=T)
  adjacency <- adjacency[, c(4, 9, 7, 6, 8, 10)]
  head(adjacency)
  names(adjacency)[1] <- "county_fips"
  names(adjacency)[4] <- "neighbor_county_fips"
  adjacency <- adjacency[adjacency$county_fips != adjacency$neighbor_county_fips,]
  
  county.fips <- unique(adjacency[,1:3])
  populations <- read.csv("co-est2019-alldata.csv") %>% filter(COUNTY != 0) %>%
    select(SUMLEV, REGION, STATE, COUNTY, STNAME, CTYNAME, POPESTIMATE2018, POPESTIMATE2019)
  populations <- merge(populations,
                       read.csv("co-est2021-alldata.csv") %>% filter(COUNTY != 0) %>% select(STATE, COUNTY, POPESTIMATE2020, POPESTIMATE2021),
                       by=c("STATE", "COUNTY"))
  
  names(populations)[5] <- "state"
  names(populations)[6] <- "county"
  county.fips <- merge(county.fips, populations[,5:10], by=c("state", "county"), all.x=T) %>% as_tibble
  names(county.fips)[3] <- "GEOID"
  
  populations <- populations %>% 
    rename(`2018`=POPESTIMATE2018,
           `2019`=POPESTIMATE2019,
           `2020`=POPESTIMATE2020,
           `2021`=POPESTIMATE2021) %>% 
    pivot_longer(`2018`:`2021`, names_to="Year", values_to="population") %>%
    mutate(Year=as.integer(Year)) %>% 
    left_join(county.fips[,1:3], by=c("state", "county")) %>% 
    select(state:GEOID)
}

overdose <- read.csv("VSRR_Provisional_County-Level_Drug_Overdose_Death_Counts.csv") %>% as_tibble
names(overdose)[1] <- "state_abbv"
overdose <- overdose %>% filter(Year %in% 2018:2021 & !(state_abbv %in% c("AK", "HI", "US")))
overdose$Provisional.Drug.Overdose.Deaths <- as.numeric(overdose$Provisional.Drug.Overdose.Deaths)
overdose$Month <- as.numeric(overdose$Month)
overdose <- overdose %>% complete(Year, Month)
overdose <- overdose %>% 
  mutate(month_year=str_c(Month, Year, sep="-") %>% 
           parse_date("%m-%Y") %>% format("%b_%Y"))
overdose <- overdose %>% rename(state=STATE_NAME, county=COUNTYNAME, GEOID=FIPS)
overdose <- overdose %>% 
  select(GEOID, month_year, Provisional.Drug.Overdose.Deaths) %>% rename(Month_Year=month_year) %>% 
  rename(overdose=Provisional.Drug.Overdose.Deaths)

crack.rel <- read.csv("crack count KNN5 R codes 999 upper tail (relative).csv") %>% as_tibble
crack.LISA_C.rel <- crack.rel[,c(1:3, grep("LISA_C", names(crack.rel)))]
crack.LISA_C.rel <- crack.LISA_C.rel[,-(4:27)]
names(crack.LISA_C.rel)[-(1:3)] <- unique(overdose$Month_Year)
crack.LISA_C.rel <- crack.LISA_C.rel %>% 
  pivot_longer(-c(state, county, GEOID), names_to="Month_Year", values_to="crack_HH") %>% 
  select(GEOID:crack_HH)

cocaine.rel <- read.csv("cocaine other counts KNN5 R codes 999 upper tail (relative).csv") %>% as_tibble
cocaine.LISA_C.rel <- cocaine.rel[,c(1:3, grep("LISA_C", names(cocaine.rel)))]
cocaine.LISA_C.rel <- cocaine.LISA_C.rel[,-(4:27)]
names(cocaine.LISA_C.rel)[-(1:3)] <- unique(overdose$Month_Year)
cocaine.LISA_C.rel <- cocaine.LISA_C.rel %>%
  pivot_longer(-c(state, county, GEOID), names_to="Month_Year", values_to="cocaine_HH") %>% 
  select(GEOID:cocaine_HH)

crack.rel <- crack.rel %>% select(GEOID, Jan_2020:Dec_2021) %>% 
  pivot_longer(-GEOID, names_to="Month_Year", values_to="crack_count")
cocaine.rel <- cocaine.rel %>% select(GEOID, Jan_2020:Dec_2021) %>% 
  pivot_longer(-GEOID, names_to="Month_Year", values_to="cocine_count")

all.drugs.count <- read.csv("all_drugs count HIDTA (02-09-2023).csv") %>% as_tibble %>% 
  select(GEOID, Jan_2020:Dec_2021) %>% 
  pivot_longer(-GEOID, names_to="Month_Year", values_to="all_count")
all.drugs.weight <- read.csv("all_drugs weight HIDTA (02-09-2023).csv") %>% as_tibble %>% 
  select(GEOID, Jan_2020:Dec_2021) %>% 
  pivot_longer(-GEOID, names_to="Month_Year", values_to="all_weight")

populations <- populations %>% filter(Year == 2020) %>% select(-Year)
overdose_integrated <- left_join(overdose, crack.LISA_C.rel, by=c("GEOID", "Month_Year")) %>% 
  left_join(populations, by=c("GEOID")) %>% 
  mutate(overdose.per.pop=overdose/population) %>% 
  left_join(cocaine.LISA_C.rel, by=c("GEOID", "Month_Year")) %>% 
  left_join(crack.rel, by=c("GEOID", "Month_Year")) %>% 
  left_join(cocaine.rel, by=c("GEOID", "Month_Year")) %>% 
  left_join(all.drugs.count, by=c("GEOID", "Month_Year")) %>% 
  left_join(all.drugs.weight, by=c("GEOID", "Month_Year"))


cor(overdose_integrated %>% select(overdose, crack_count:all_weight) %>% na.omit)
cor(cbind(overdose_integrated$overdose.per.pop, 
          overdose_integrated$crack_count/overdose_integrated$population) %>% na.omit)
cor(cbind(overdose_integrated$overdose.per.pop, 
          overdose_integrated$cocine_count/overdose_integrated$population) %>% na.omit)
cor(cbind(overdose_integrated$overdose.per.pop, 
          overdose_integrated$all_count/overdose_integrated$population) %>% na.omit)
cor(cbind(overdose_integrated$overdose.per.pop, 
          overdose_integrated$all_weight/overdose_integrated$population) %>% na.omit)

overdose_integrated %>% group_by(crack_HH) %>%
  summarise(avg.overdose=mean(overdose, na.rm=T),
            med.overdose=median(overdose, na.rm=T))

overdose_integrated %>% group_by(crack_HH) %>%
  summarise(avg.overdose.per.pop=mean(overdose.per.pop, na.rm=T),
            med.overdose.per.pop=median(overdose.per.pop, na.rm=T),
            var.overdose.per.pop=var(overdose.per.pop, na.rm=T))
var(overdose_integrated$overdose.per.pop, na.rm=T)

t.test(overdose_integrated %>% filter(crack_HH==0) %>% pull(overdose.per.pop),
       overdose_integrated %>% filter(crack_HH==1) %>% pull(overdose.per.pop))
