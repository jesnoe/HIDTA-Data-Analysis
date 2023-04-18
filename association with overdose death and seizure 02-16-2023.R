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
populations <- populations %>% filter(Year == 2020) %>% select(-Year)

add_drug_info <- function(drug_data, drug_name, data_to_add) {
  drug.LISA_C.rel <- drug_data[,c(1:3, grep("LISA_C", names(drug_data)))]
  drug.LISA_C.rel <- drug.LISA_C.rel[,-(4:27)]
  names(drug.LISA_C.rel)[-(1:3)] <- unique(overdose$Month_Year)
  drug.LISA_C.rel <- drug.LISA_C.rel %>% 
    pivot_longer(-c(state, county, GEOID), names_to="Month_Year", values_to=paste(drug_name, "HH", sep="_")) %>% 
    select(3:5)
  drug.count <- drug_data %>% select(GEOID, Jan_2020:Dec_2021) %>% 
    pivot_longer(-GEOID, names_to="Month_Year", values_to=paste(drug_name, "count", sep="_"))
  result <- left_join(data_to_add, drug.LISA_C.rel, by=c("GEOID", "Month_Year")) %>%
    left_join(drug.count, by=c("GEOID", "Month_Year"))
  return(result)
}

crack <- read.csv("crack count KNN5 R codes 999 upper tail (relative).csv") %>% as_tibble
other_cocaine <- read.csv("cocaine other counts KNN5 R codes 999 upper tail (relative).csv") %>% as_tibble
meth <- read.csv("meth counts KNN5 R codes 999 upper tail (relative).csv") %>% as_tibble
marijuana <- read.csv("marijuana count KNN5 R codes 999 upper tail (relative).csv") %>% as_tibble

overdose_integrated <- left_join(overdose, populations, by=c("GEOID")) %>% 
  mutate(overdose.per.pop=overdose/population)
  
overdose_integrated <- add_drug_info(crack, "crack", overdose_integrated)
overdose_integrated <- add_drug_info(other_cocaine, "cocaine", overdose_integrated)
overdose_integrated <- add_drug_info(meth, "meth", overdose_integrated)
overdose_integrated <- add_drug_info(marijuana, "marijuana", overdose_integrated)

all.drugs.count <- read.csv("all_drugs count HIDTA (02-09-2023).csv") %>% as_tibble %>% 
  select(GEOID, Jan_2020:Dec_2021) %>% 
  pivot_longer(-GEOID, names_to="Month_Year", values_to="all_count")
all.drugs.weight <- read.csv("all_drugs weight HIDTA (02-09-2023).csv") %>% as_tibble %>% 
  select(GEOID, Jan_2020:Dec_2021) %>% 
  pivot_longer(-GEOID, names_to="Month_Year", values_to="all_weight")


cor(overdose_integrated[,c(3, 7:15)] %>% na.omit) %>% write.csv("correlation overdose vs seizures.csv", row.names=F)

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

overdose_integrated %>% group_by(meth_HH) %>%
  summarise(avg.overdose=mean(overdose, na.rm=T),
            med.overdose=median(overdose, na.rm=T))

overdose_integrated %>% group_by(meth_HH) %>%
  summarise(avg.overdose.per.pop=mean(overdose.per.pop, na.rm=T),
            med.overdose.per.pop=median(overdose.per.pop, na.rm=T),
            var.overdose.per.pop=var(overdose.per.pop, na.rm=T))

overdose_integrated %>% group_by(marijuana_HH) %>%
  summarise(avg.overdose=mean(overdose, na.rm=T),
            med.overdose=median(overdose, na.rm=T))

overdose_integrated %>% group_by(marijuana_HH) %>%
  summarise(avg.overdose.per.pop=mean(overdose.per.pop, na.rm=T),
            med.overdose.per.pop=median(overdose.per.pop, na.rm=T),
            var.overdose.per.pop=var(overdose.per.pop, na.rm=T))

t.test(overdose_integrated %>% filter(crack_HH==0) %>% pull(overdose.per.pop),
       overdose_integrated %>% filter(crack_HH==1) %>% pull(overdose.per.pop))
