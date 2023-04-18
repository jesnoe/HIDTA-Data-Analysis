library(fpp2)
library(spdep)
library(readxl)
library(urbnmapr)
library(tidyverse)
library(gridExtra)

# urbnmapr data
unique(countydata$county_fips) # some ips have 5 digits. need to be converted into integer.
HIDTA.dist <- read.csv("HIDTA Regions.csv", header=T) %>% as_tibble
coordinate.HIDTA <- left_join(counties, HIDTA.dist, by=c("state_name", "county_name")) %>% filter(!(state_name %in% c("Alaska", "Hawaii", "Puerto Rico")))
names(coordinate.HIDTA)[7] <- "GEOID"

HIDTA.regions <- coordinate.HIDTA %>% filter(!is.na(HIDTA)) %>% select(GEOID, HIDTA) %>% unique

## Too small observations for non-all_drugs cocaine
seizures <- read_xlsx("Drug Seizures All HIDTAs All Drugs 2018-2021 Combined.xlsx")
seizures$County <- substring(seizures$County, 1, str_locate(seizures$County, ",")[,1]-1)
seizures$Year <- substring(seizures$SeizureDate, 1, 4)
seizures$Month <- substring(seizures$SeizureDate, 6, 7)
seizures$Year <- as.numeric(seizures$Year)
seizures$Month <- as.numeric(seizures$Month)

# all drugs count
{
all_drugs <- seizures %>% filter(!(State %in% c("Alaska", "Hawaii", "U.S. Virgin Islands", "France")) & !is.na(County) )
calender <- data.frame(Month=1:12, Month2=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
all_drugs <- inner_join(all_drugs, calender, by="Month")

all_drugs <- all_drugs %>% group_by(State, County, Year, Month2) %>% summarise(count=sum(Quantity > 0)) %>%
  pivot_wider(names_from=c(Month2, Year), names_sep="_", values_from=count, values_fill=0)
cal_order <- c(str_c(calender$Month2, 2018, sep="_"),
               str_c(calender$Month2, 2019, sep="_"),
               str_c(calender$Month2, 2020, sep="_"),
               str_c(calender$Month2, 2021, sep="_"))
all_drugs <- all_drugs[, c(1:2, match(cal_order, names(all_drugs)))]
names(all_drugs)[1:2] <- c("state_name", "county")
coordinate.HIDTA$county <- str_split(coordinate.HIDTA$county_name, " ") %>% lapply(function(x) str_c(c(x[(1:(length(x)-1))]), collapse=" ") ) %>% unlist
all_drugs <- left_join(all_drugs, unique(coordinate.HIDTA[,c(7, 12, 14, 15)]), by=c("state_name", "county")) %>%
  select(state_name, county, GEOID, HIDTA, Jan_2018:Dec_2021)
all_drugs <- rbind(all_drugs, left_join(unique(coordinate.HIDTA[!is.na(coordinate.HIDTA$HIDTA),c(7, 12, 14, 15)]), all_drugs,
                                by=c("state_name", "county", "GEOID")) %>%
                       select(state_name, county, GEOID, HIDTA.x, Jan_2018:Dec_2021) %>% rename(HIDTA=HIDTA.x) %>%
                       filter(is.na(Jan_2018))) %>% arrange(GEOID)
names(all_drugs)[1] <- "state"
all_drugs[,5:52][is.na(all_drugs[,5:52])] <- 0
}

grep("Jan_2018", names(all_drugs)) # 5
grep("Dec_2021", names(all_drugs)) # 52
all_drugs <- all_drugs[(all_drugs[,5:52] %>% apply(1, sum) > 0) | (!is.na(all_drugs$HIDTA)),]
all_drugs[which(all_drugs$state == "District of Columbia"), 3] <- 11001
all_drugs <- all_drugs %>% filter(!is.na(GEOID))
all_drugs # 1515 counties in total.
# write.csv(all_drugs, "all_drugs count HIDTA (02-09-2023).csv", row.names=F)

# all drugs weight
{
all_drugs <- seizures %>% filter(!(State %in% c("Alaska", "Hawaii", "U.S. Virgin Islands", "France")) & !is.na(County) )
calender <- data.frame(Month=1:12, Month2=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
all_drugs <- inner_join(all_drugs, calender, by="Month")

all_drugs <- all_drugs %>% group_by(State, County, Year, Month2) %>% summarise(weight=sum(Quantity)) %>%
  pivot_wider(names_from=c(Month2, Year), names_sep="_", values_from=weight, values_fill=0)
cal_order <- c(str_c(calender$Month2, 2018, sep="_"),
               str_c(calender$Month2, 2019, sep="_"),
               str_c(calender$Month2, 2020, sep="_"),
               str_c(calender$Month2, 2021, sep="_"))
all_drugs <- all_drugs[, c(1:2, match(cal_order, names(all_drugs)))]
names(all_drugs)[1:2] <- c("state_name", "county")
coordinate.HIDTA$county <- str_split(coordinate.HIDTA$county_name, " ") %>% lapply(function(x) str_c(c(x[(1:(length(x)-1))]), collapse=" ") ) %>% unlist
all_drugs <- left_join(all_drugs, unique(coordinate.HIDTA[,c(7, 12, 14, 15)]), by=c("state_name", "county")) %>%
  select(state_name, county, GEOID, HIDTA, Jan_2018:Dec_2021)
all_drugs <- rbind(all_drugs, left_join(unique(coordinate.HIDTA[!is.na(coordinate.HIDTA$HIDTA),c(7, 12, 14, 15)]), all_drugs,
                                by=c("state_name", "county", "GEOID")) %>%
                       select(state_name, county, GEOID, HIDTA.x, Jan_2018:Dec_2021) %>% rename(HIDTA=HIDTA.x) %>%
                       filter(is.na(Jan_2018))) %>% arrange(GEOID)
names(all_drugs)[1] <- "state"
all_drugs[,5:52][is.na(all_drugs[,5:52])] <- 0
}

grep("Jan_2018", names(all_drugs)) # 5
grep("Dec_2021", names(all_drugs)) # 52
all_drugs <- all_drugs[(all_drugs[,5:52] %>% apply(1, sum) > 0) | (!is.na(all_drugs$HIDTA)),]
all_drugs[which(all_drugs$state == "District of Columbia"), 3] <- 11001
all_drugs <- all_drugs %>% filter(!is.na(GEOID)) %>% arrange(GEOID)
all_drugs # 1515 counties in total.
# write.csv(all_drugs, "all_drugs weight HIDTA (02-09-2023).csv", row.names=F)

# All drugs seizure
counties.obs <- counties
names(counties.obs)[7] <- "GEOID"
counties.obs <- counties.obs %>% rename(state=state_name, county=county_name)
counties.obs$GEOID <- as.numeric(counties.obs$GEOID)
counties.obs <- counties.obs %>% filter(!(state %in% c("Alaska", "Hawaii")))
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

# All drugs count maps
LISA.rel <- read.csv("all_drugs count HIDTA (02-09-2023).csv") %>% as_tibble

LISA.rel.map <- left_join(counties.obs[, c(1:2,6:7,10,12)], LISA.rel[, c(3, 5:52)], by = "GEOID")
LISA.rel.map <- LISA.rel.map %>% 
  pivot_longer(c(-long, -lat, -group, -GEOID, -state, -county), names_to="Month_Year", values_to="seizure_counts") %>% 
  mutate(Year=as.integer(substr(Month_Year, 5,8)))
LISA.rel.map <- left_join(LISA.rel.map, populations[, -(1:2)], by=c("GEOID", "Year"))
LISA.rel.map$counts_per_pop <- LISA.rel.map$seizure_counts/LISA.rel.map$population

LISA.rel.map$Month_Year <- parse_date(LISA.rel.map$Month_Year, "%b_%Y")

for (year in 2020:2021) { # monthly maps of seizure counts
  LISA.rel.map %>% filter(year(Month_Year) == year & month(Month_Year) %in% 1:4) %>% 
    ggplot(mapping = aes(long, lat, group = group, fill=seizure_counts)) +
    geom_polygon(color = "#000000", size = .05) +
    facet_wrap(. ~ substr(Month_Year, 1, 7)) +
    scale_fill_viridis_c(na.value="white") +
    labs(fill = "Seizure Counts", title="All Drugs Seizure Counts") + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) -> seizure_counts_map
  ggsave(paste("All drugs seizure_counts_", year, " (4 months).jpg", sep=""), seizure_counts_map, width=20, height=15, units="cm")
  
  LISA.rel.map %>% filter(grepl(year, Month_Year) & month(Month_Year) %in% 1:4) %>% 
    ggplot(mapping = aes(long, lat, group = group, fill=counts_per_pop)) +
    geom_polygon(color = "#000000", size = .05) +
    facet_wrap(. ~ substr(Month_Year, 1, 7)) +
    scale_fill_viridis_c(na.value="white") +
    labs(fill = "Counts Per Pop.", title="All Drugs Seizure Counts/Pop.") + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) -> counts_per_pop.map
  ggsave(paste("All drugs counts_per_pop_", year, " (4 months).jpg", sep=""), counts_per_pop.map, width=20, height=15, units="cm")
}

# All drugs weight maps
LISA.rel <- read.csv("all_drugs weight HIDTA (02-09-2023).csv") %>% as_tibble

LISA.rel.map <- left_join(counties.obs[, c(1:2,6:7,10,12)], LISA.rel[, c(3, 5:52)], by = "GEOID")
LISA.rel.map <- LISA.rel.map %>% 
  pivot_longer(c(-long, -lat, -group, -GEOID, -state, -county), names_to="Month_Year", values_to="seizure_counts") %>% 
  mutate(Year=as.integer(substr(Month_Year, 5,8)))
LISA.rel.map <- left_join(LISA.rel.map, populations[, -(1:2)], by=c("GEOID", "Year"))
LISA.rel.map$counts_per_pop <- LISA.rel.map$seizure_counts/LISA.rel.map$population

LISA.rel.map$Month_Year <- parse_date(LISA.rel.map$Month_Year, "%b_%Y")

for (year in 2020:2021) { # monthly maps of seizure weights
  LISA.rel.map %>% filter(year(Month_Year) == year & month(Month_Year) %in% 1:4) %>% 
    ggplot(mapping = aes(long, lat, group = group, fill=seizure_counts)) +
    geom_polygon(color = "#000000", size = .05) +
    facet_wrap(. ~ substr(Month_Year, 1, 7)) +
    scale_fill_viridis_c(na.value="white") +
    labs(fill = "Seizure Weights", title="All Drugs Seizure Weights") + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) -> seizure_counts_map
  ggsave(paste("All drugs seizure_weights_", year, " (4 months).jpg", sep=""), seizure_counts_map, width=20, height=15, units="cm")
  
  LISA.rel.map %>% filter(grepl(year, Month_Year) & month(Month_Year) %in% 1:4) %>% 
    ggplot(mapping = aes(long, lat, group = group, fill=counts_per_pop)) +
    geom_polygon(color = "#000000", size = .05) +
    facet_wrap(. ~ substr(Month_Year, 1, 7)) +
    scale_fill_viridis_c(na.value="white") +
    labs(fill = "Weights Per Pop.", title="All Drugs Seizure Weights/Pop.") + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) -> counts_per_pop.map
  ggsave(paste("All drugs weights_per_pop_", year, " (4 months).jpg", sep=""), counts_per_pop.map, width=20, height=15, units="cm")
}
