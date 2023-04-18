library(fpp2)
library(spdep)
library(readxl)
library(urbnmapr)
library(tidyverse)
library(gridExtra)
# library(maps)
# library(rnaturalearth)

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
}
overdose <- read.csv("VSRR_Provisional_Drug_Overdose_Death_Counts.csv") %>% as_tibble
cocaine <- read.csv("cocaine HIDTA (12-07-2022).csv") %>% as_tibble
crack <- read.csv("crack cocaine HIDTA (12-07-2022).csv") %>% as_tibble
cocaine_count <- read.csv("cocaine count HIDTA (12-07-2022).csv") %>% as_tibble
crack_count <- read.csv("crack count HIDTA (12-07-2022).csv") %>% as_tibble

# urbnmapr data
unique(countydata$county_fips) # some ips have 5 digits. need to be converted into integer.

HIDTA.dist <- read.csv("HIDTA Regions.csv", header=T) %>% as_tibble
coordinate.HIDTA <- left_join(counties, HIDTA.dist, by=c("state_name", "county_name")) %>% filter(!(state_name %in% c("Alaska", "Hawaii", "Puerto Rico")))
names(coordinate.HIDTA)[7] <- "GEOID"

HIDTA.regions <- coordinate.HIDTA %>% filter(!is.na(HIDTA)) %>% select(GEOID, HIDTA) %>% unique

cocaine <- read.csv("CocaineLISAResults.csv") %>% as_tibble %>% arrange(GEOID)
state.fips <- read.csv("us-state-ansi-fips.csv")
names(state.fips)[1] <- "state"
names(state.fips)[2] <- "STATEFP"
cocaine <- merge(cocaine, state.fips[,1:2], by="STATEFP") %>% as_tibble
cocaine <- cocaine[,c(1, 212, 2:211)]

counties.obs <- counties %>% filter(county_fips %in% unique(cocaine$GEOID))
names(counties.obs)[7] <- "GEOID"
grep("Jan_2018", names(cocaine)) # 21
grep("Dec_2021", names(cocaine)) # 68

cocaine$HIDTA <- left_join(cocaine, HIDTA.regions, by="GEOID")$HIDTA
cocaine %>% filter(!is.na(HIDTA)) %>% group_by(state, NAMELSAD) %>% summarise(total_seizures=sum(Jan_2018:Dec_2021)) %>% arrange(total_seizures)

# Counties with 0 kg seizure for 48 months are unobserved
# cocaine[cocaine[,21:68] %>% apply(1, sum) > 0,] %>% nrow # 971
# cocaine[(cocaine[,21:68] %>% apply(1, sum) > 0) | (!is.na(cocaine$HIDTA)),] %>% nrow #1092
# names(cocaine)[8] <- "county"
# cocaine <- cocaine[(cocaine[,21:68] %>% apply(1, sum) > 0) | (!is.na(cocaine$HIDTA)),] %>%
#   select(state, county, GEOID, HIDTA, Jan_2018:Dec_2021)
# write.csv(cocaine, "cocaine HIDTA (12-07-2022).csv", row.names=F)
#  
# {
# seizures <- read_xlsx("Drug Seizures All HIDTAs All Drugs 2018-2021 Combined.xlsx")
# seizures$County <- substring(seizures$County, 1, str_locate(seizures$County, ",")[,1]-1)
# seizures$Year <- substring(seizures$SeizureDate, 1, 4)
# seizures$Month <- substring(seizures$SeizureDate, 6, 7)
# seizures$Year <- as.numeric(seizures$Year)
# seizures$Month <- as.numeric(seizures$Month)
# 
# crack <- seizures %>% filter(Drug=="Crack" & !(State %in% c("Alaska", "Hawaii", "U.S. Virgin Islands", "France")) & !is.na(County) )
# calander <- data.frame(Month=1:12, Month2=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
# crack <- inner_join(crack, calander, by="Month")
# 
# crack <- crack %>% group_by(State, County, Year, Month2) %>% summarise(weight=sum(Quantity)) %>%
#   pivot_wider(names_from=c(Month2, Year), names_sep="_", values_from=weight, values_fill=0)
# cal_order <- c(str_c(calander$Month2, 2018, sep="_"),
#                str_c(calander$Month2, 2019, sep="_"),
#                str_c(calander$Month2, 2020, sep="_"),
#                str_c(calander$Month2, 2021, sep="_"))
# crack <- crack[, c(1:2, match(cal_order, names(crack)))]
# names(crack)[1:2] <- c("state_name", "county")
# coordinate.HIDTA$county <- str_split(coordinate.HIDTA$county_name, " ") %>% lapply(function(x) str_c(c(x[(1:(length(x)-1))]), collapse=" ") ) %>% unlist
# crack <- left_join(crack, unique(coordinate.HIDTA[,c(7, 12, 14, 15)]), by=c("state_name", "county")) %>%
#   select(state_name, county, GEOID, HIDTA, Jan_2018:Dec_2021)
# crack <- rbind(crack, left_join(unique(coordinate.HIDTA[!is.na(coordinate.HIDTA$HIDTA),c(7, 12, 14, 15)]), crack,
#                                 by=c("state_name", "county", "GEOID")) %>%
#                        select(state_name, county, GEOID, HIDTA.x, Jan_2018:Dec_2021) %>% rename(HIDTA=HIDTA.x) %>%
#                        filter(is.na(Jan_2018))) %>% arrange(GEOID)
# names(crack)[1] <- "state"
# crack[,5:52][is.na(crack[,5:52])] <- 0
# }
# 
# grep("Jan_2018", names(crack)) # 5
# grep("Dec_2021", names(crack)) # 52
# crack <- crack[(crack[,5:52] %>% apply(1, sum) > 0) | (!is.na(crack$HIDTA)),]
# crack[which(crack$state == "District of Columbia"), 3] <- 11001
# crack <- crack %>% filter(!is.na(GEOID))
# crack # 813 counties in total.
# write.csv(crack, "crack cocaine HIDTA (12-07-2022).csv", row.names=F)




# {
# seizures <- read_xlsx("Drug Seizures All HIDTAs All Drugs 2018-2021 Combined.xlsx")
# seizures$County <- substring(seizures$County, 1, str_locate(seizures$County, ",")[,1]-1)
# seizures$Year <- substring(seizures$SeizureDate, 1, 4)
# seizures$Month <- substring(seizures$SeizureDate, 6, 7)
# seizures$Year <- as.numeric(seizures$Year)
# seizures$Month <- as.numeric(seizures$Month)
# 
# cocaine_count <- seizures %>% 
#   filter(Drug %in% c("Cocaine", "Crack", "Cocaine Powder", "Cocaine Base", "Cocaine precursor", "Coca Leaves", "Cocaine metabolite") & 
#            !(State %in% c("Alaska", "Hawaii", "U.S. Virgin Islands", "France")) & !is.na(County) )
# calander <- data.frame(Month=1:12, Month2=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
# cal_order <- c(str_c(calander$Month2, 2018, sep="_"),
#                str_c(calander$Month2, 2019, sep="_"),
#                str_c(calander$Month2, 2020, sep="_"),
#                str_c(calander$Month2, 2021, sep="_"))
# coordinate.HIDTA$county <- str_split(coordinate.HIDTA$county_name, " ") %>% lapply(function(x) str_c(c(x[(1:(length(x)-1))]), collapse=" ") ) %>% unlist
# 
# cocaine_count <- inner_join(cocaine_count, calander, by="Month")
# cocaine_count <- cocaine_count %>% group_by(State, County, Year, Month2) %>% summarise(weight=sum(Quantity > 0)) %>%
#   pivot_wider(names_from=c(Month2, Year), names_sep="_", values_from=weight, values_fill=0)
# 
# cocaine_count <- cocaine_count[, c(1:2, match(cal_order, names(cocaine_count)))]
# names(cocaine_count)[1:2] <- c("state_name", "county")
# cocaine_count <- left_join(cocaine_count, unique(coordinate.HIDTA[,c(7, 10, 12, 14, 15)]), by=c("state_name", "county")) %>%
#   select(state_name, county_name, GEOID, HIDTA, Jan_2018:Dec_2021)
# cocaine_count <- cocaine_count[,-1]
# cocaine_count <- rbind(cocaine_count, left_join(unique(coordinate.HIDTA[!is.na(coordinate.HIDTA$HIDTA),c(7, 10, 12, 14)]), cocaine_count, by=c("GEOID")) %>%
#   select(state_name.x, county_name.x, GEOID, HIDTA.x, Jan_2018:Dec_2021) %>% rename(state_name=state_name.x, county_name=county_name.x, HIDTA=HIDTA.x) %>%
#     filter(is.na(Jan_2018))) %>% arrange(GEOID)
# names(cocaine_count)[c(1:2, 4)] <- c("state", "county", "HIDTA")
# cocaine_count[,5:52][is.na(cocaine_count[,5:52])] <- 0
# cocaine_count <- cocaine_count[(cocaine_count[,5:52] %>% apply(1, sum) > 0) | (!is.na(cocaine_count$HIDTA)),] # 1,196 counties in total.
# write.csv(cocaine_count, "cocaine count HIDTA (12-07-2022).csv", row.names=F)
# 
# 
# crack_count <- seizures %>% filter(Drug=="Crack" & !(State %in% c("Alaska", "Hawaii", "U.S. Virgin Islands", "France")) & !is.na(County) )
# crack_count <- inner_join(crack_count, calander, by="Month")
# 
# crack_count <- crack_count %>% group_by(State, County, Year, Month2) %>% summarise(weight=sum(Quantity > 0)) %>%
#   pivot_wider(names_from=c(Month2, Year), names_sep="_", values_from=weight, values_fill=0)
# 
# crack_count <- crack_count[, c(1:2, match(cal_order, names(crack_count)))]
# names(crack_count)[1:2] <- c("state_name", "county")
# crack_count <- left_join(crack_count, unique(coordinate.HIDTA[,c(7, 10, 12, 14, 15)]), by=c("state_name", "county")) %>%
#   select(state_name, county_name, GEOID, HIDTA, Jan_2018:Dec_2021)
# crack_count <- crack_count[,-1]
# crack_count <- rbind(crack_count, left_join(unique(coordinate.HIDTA[!is.na(coordinate.HIDTA$HIDTA),c(7, 10, 12, 14)]), crack_count, by=c("GEOID")) %>%
#                          select(state_name.x, county_name.x, GEOID, HIDTA.x, Jan_2018:Dec_2021) %>% rename(state_name=state_name.x, county_name=county_name.x, HIDTA=HIDTA.x) %>%
#                          filter(is.na(Jan_2018))) %>% arrange(GEOID)
# names(crack_count)[c(1:2, 4)] <- c("state", "county", "HIDTA")
# crack_count[,5:52][is.na(crack_count[,5:52])] <- 0
# }
# 
# grep("Jan_2018", names(crack_count)) # 5
# grep("Dec_2021", names(crack_count)) # 52
# crack_count <- crack_count[(crack_count[,5:52] %>% apply(1, sum) > 0) | (!is.na(crack_count$HIDTA)),] # 849 counties in total.
# write.csv(crack_count, "crack count HIDTA (12-07-2022).csv", row.names=F)


# ts plots
seizure_ts_plot <- function(seizure, state_name, drug) {
  seizure %>% filter(state==state_name) %>% pivot_longer(cols=Jan_2018:Dec_2021, names_to="month", values_to="weight") %>% 
    ggplot(aes(x=month, y=weight, group=GEOID)) +
    geom_line(aes(color=county)) +
    labs(title=str_c(state_name, drug, sep="-"),
         xlab="Month-Year") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}

names(overdose)[1] <- "state_abbv"
cocaine_overdose <- overdose %>% filter(Indicator == "Cocaine (T40.5)" & Year %in% 2018:2021 & !(state_abbv %in% c("AK", "HI", "US")))
cocaine_overdose$Data.Value <- as.numeric(cocaine_overdose$Data.Value)
cocaine_overdose$Month <- str_sub(cocaine_overdose$Month, 1, 3)
cocaine_overdose$month <- str_c(cocaine_overdose$Month, cocaine_overdose$Year, sep="_")
cocaine_overdose$month <- factor(cocaine_overdose$month, levels=names(cocaine)[-(1:4)])
cocaine_overdose <- cocaine_overdose %>% rename(state=State.Name)
cocaine_overdose %>% complete(Year, Month)



state_name <- "California"
grid.arrange(seizure_ts_plot(cocaine, state_name, "cocaine"),
             seizure_ts_plot(crack, state_name, "crack"), ncol=2)
grid.arrange(seizure_ts_plot(cocaine_count, state_name, "cocaine"),
             seizure_ts_plot(crack_count, state_name, "crack"), ncol=2)
grid.arrange(
  cocaine_overdose %>% filter(state == state_name) %>% 
    ggplot(aes(x=month, y=Data.Value)) +
    geom_point(na.rm=T) +
    labs(title=str_c(state_name, "cocaine overdose", sep="-"),
         x="Month-Year", y="Counts") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)),
  cocaine %>% filter(state==state_name) %>% select(Jan_2018:Dec_2021) %>% 
    apply(2,sum) %>%
    ts(start=c(2018, 1), frequency=12) %>% 
    autoplot() + 
    labs(title=str_c(state_name, "all cocaine", sep="-"),
         x="Year", y="weight"),
  cocaine_count %>% filter(state==state_name) %>% select(Jan_2018:Dec_2021) %>% 
    apply(2,sum) %>%
    ts(start=c(2018, 1), frequency=12) %>% 
    autoplot() + 
    labs(title=str_c(state_name, "all cocaine", sep="-"),
         x="Year", y="count"), ncol=3
)

state_name <- "Florida"
grid.arrange(seizure_ts_plot(cocaine, state_name, "cocaine"),
             seizure_ts_plot(crack, state_name, "crack"), ncol=2)
grid.arrange(
  cocaine_overdose %>% filter(state == state_name) %>% 
    ggplot(aes(x=month, y=Data.Value)) +
    geom_point(na.rm=T) +
    labs(title=str_c(state_name, "cocaine overdose", sep="-"),
         x="Month-Year", y="Counts") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)),
  cocaine %>% filter(state==state_name) %>% select(Jan_2018:Dec_2021) %>% 
    apply(2,sum) %>%
    ts(start=c(2018, 1), frequency=12) %>% 
    autoplot() + 
    labs(title=str_c(state_name, "all cocaine", sep="-"),
         x="Year", y="weight"),
  cocaine_count %>% filter(state==state_name) %>% select(Jan_2018:Dec_2021) %>% 
    apply(2,sum) %>%
    ts(start=c(2018, 1), frequency=12) %>% 
    autoplot() + 
    labs(title=str_c(state_name, "all cocaine", sep="-"),
         x="Year", y="count"), ncol=3
)


state_name <- "New York"
grid.arrange(seizure_ts_plot(cocaine, state_name, "cocaine"),
             seizure_ts_plot(crack, state_name, "crack"), ncol=2)
grid.arrange(
  cocaine_overdose %>% filter(state == state_name) %>% 
    ggplot(aes(x=month, y=Data.Value)) +
    geom_point(na.rm=T) +
    labs(title=str_c(state_name, "cocaine overdose", sep="-"),
         x="Month-Year", y="Counts") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)),
  cocaine %>% filter(state==state_name) %>% select(Jan_2018:Dec_2021) %>% 
    apply(2,sum) %>%
    ts(start=c(2018, 1), frequency=12) %>% 
    autoplot() + 
    labs(title=str_c(state_name, "all cocaine", sep="-"),
         x="Year", y="weight"),
  cocaine_count %>% filter(state==state_name) %>% select(Jan_2018:Dec_2021) %>% 
    apply(2,sum) %>%
    ts(start=c(2018, 1), frequency=12) %>% 
    autoplot() + 
    labs(title=str_c(state_name, "all cocaine", sep="-"),
         x="Year", y="count"), ncol=3
)

