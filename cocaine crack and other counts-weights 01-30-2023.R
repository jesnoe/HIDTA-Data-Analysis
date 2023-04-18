library(fpp2)
library(spdep)
library(readxl)
library(urbnmapr)
library(tidyverse)
library(gridExtra)

overdose <- read.csv("VSRR_Provisional_County-Level_Drug_Overdose_Death_Counts.csv") %>% as_tibble
cocaine <- read.csv("cocaine other weight HIDTA (01-30-2023).csv") %>% as_tibble
crack <- read.csv("cocaine crack weight HIDTA (01-30-2023).csv") %>% as_tibble
cocaine_count <- read.csv("cocaine other count HIDTA (01-30-2023).csv") %>% as_tibble
crack_count <- read.csv("cocaine crack count HIDTA (12-20-2022).csv") %>% as_tibble

# urbnmapr data
unique(countydata$county_fips) # some ips have 5 digits. need to be converted into integer.
HIDTA.dist <- read.csv("HIDTA Regions.csv", header=T) %>% as_tibble
coordinate.HIDTA <- left_join(counties, HIDTA.dist, by=c("state_name", "county_name")) %>% filter(!(state_name %in% c("Alaska", "Hawaii", "Puerto Rico")))
names(coordinate.HIDTA)[7] <- "GEOID"

HIDTA.regions <- coordinate.HIDTA %>% filter(!is.na(HIDTA)) %>% select(GEOID, HIDTA) %>% unique

## Too small observations for non-crack cocaine
seizures <- read_xlsx("Drug Seizures All HIDTAs All Drugs 2018-2021 Combined.xlsx")
seizures$County <- substring(seizures$County, 1, str_locate(seizures$County, ",")[,1]-1)
seizures$Year <- substring(seizures$SeizureDate, 1, 4)
seizures$Month <- substring(seizures$SeizureDate, 6, 7)
seizures$Year <- as.numeric(seizures$Year)
seizures$Month <- as.numeric(seizures$Month)

seizures %>%
  filter(Drug=="Cocaine" & !(State %in% c("Alaska", "Hawaii", "U.S. Virgin Islands", "France")) & !is.na(County) ) %>%
  nrow # 52626

seizures %>%
  filter(Drug=="Crack" & !(State %in% c("Alaska", "Hawaii", "U.S. Virgin Islands", "France")) & !is.na(County) ) %>%
  nrow # 22654

seizures %>%
  filter(Drug=="Cocaine, Powder" & !(State %in% c("Alaska", "Hawaii", "U.S. Virgin Islands", "France")) & !is.na(County) ) %>%
  nrow # 54

seizures %>%
  filter(Drug=="Cocaine Base" & !(State %in% c("Alaska", "Hawaii", "U.S. Virgin Islands", "France")) & !is.na(County) ) %>%
  nrow # 33

seizures %>%
  filter(Drug=="Cocaine precursor" & !(State %in% c("Alaska", "Hawaii", "U.S. Virgin Islands", "France")) & !is.na(County) ) %>%
  nrow # 5

seizures %>%
  filter(Drug=="Coca Leaves" & !(State %in% c("Alaska", "Hawaii", "U.S. Virgin Islands", "France")) & !is.na(County) ) %>%
  nrow # 3

seizures %>%
  filter(Drug=="Cocaine metabolite" & !(State %in% c("Alaska", "Hawaii", "U.S. Virgin Islands", "France")) & !is.na(County) ) %>%
  nrow # 2
# 
# crack count
{
crack <- seizures %>% filter(Drug=="Crack" & !(State %in% c("Alaska", "Hawaii", "U.S. Virgin Islands", "France")) & !is.na(County) )
calender <- data.frame(Month=1:12, Month2=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
crack <- inner_join(crack, calender, by="Month")

crack <- crack %>% group_by(State, County, Year, Month2) %>% summarise(count=sum(Quantity > 0)) %>%
  pivot_wider(names_from=c(Month2, Year), names_sep="_", values_from=count, values_fill=0)
cal_order <- c(str_c(calender$Month2, 2018, sep="_"),
               str_c(calender$Month2, 2019, sep="_"),
               str_c(calender$Month2, 2020, sep="_"),
               str_c(calender$Month2, 2021, sep="_"))
crack <- crack[, c(1:2, match(cal_order, names(crack)))]
names(crack)[1:2] <- c("state_name", "county")
coordinate.HIDTA$county <- str_split(coordinate.HIDTA$county_name, " ") %>% lapply(function(x) str_c(c(x[(1:(length(x)-1))]), collapse=" ") ) %>% unlist
crack <- left_join(crack, unique(coordinate.HIDTA[,c(7, 12, 14, 15)]), by=c("state_name", "county")) %>%
  select(state_name, county, GEOID, HIDTA, Jan_2018:Dec_2021)
crack <- rbind(crack, left_join(unique(coordinate.HIDTA[!is.na(coordinate.HIDTA$HIDTA),c(7, 12, 14, 15)]), crack,
                                by=c("state_name", "county", "GEOID")) %>%
                       select(state_name, county, GEOID, HIDTA.x, Jan_2018:Dec_2021) %>% rename(HIDTA=HIDTA.x) %>%
                       filter(is.na(Jan_2018))) %>% arrange(GEOID)
names(crack)[1] <- "state"
crack[,5:52][is.na(crack[,5:52])] <- 0
}

grep("Jan_2018", names(crack)) # 5
grep("Dec_2021", names(crack)) # 52
crack <- crack[(crack[,5:52] %>% apply(1, sum) > 0) | (!is.na(crack$HIDTA)),]
crack[which(crack$state == "District of Columbia"), 3] <- 11001
crack <- crack %>% filter(!is.na(GEOID))
crack # 1518 counties in total.
# write.csv(crack, "cocaine crack count HIDTA (12-20-2022).csv", row.names=F)

# # other cocaine count
# {
#   cocaine <- seizures %>% filter(Drug %in% c("Cocaine", "Cocaine, Powder", "Cocaine Base",
#                                              "Cocaine precursor", "Coca Leaves", "Cocaine metabolite") &
#                                    !(State %in% c("Alaska", "Hawaii", "U.S. Virgin Islands", "France")) &
#                                    !is.na(County) )
#   cocaine <- inner_join(cocaine, calender, by="Month")
# 
#   cocaine <- cocaine %>% group_by(State, County, Year, Month2) %>% summarise(count=sum(Quantity > 0)) %>%
#     pivot_wider(names_from=c(Month2, Year), names_sep="_", values_from=count, values_fill=0)
#   cocaine <- cocaine[, c(1:2, match(cal_order, names(cocaine)))]
#   names(cocaine)[1:2] <- c("state_name", "county")
#   cocaine <- left_join(cocaine, unique(coordinate.HIDTA[,c(7, 10, 12, 14, 15)]), by=c("state_name", "county")) %>%
#     select(state_name, county, GEOID, HIDTA, Jan_2018:Dec_2021)
#   cocaine <- rbind(cocaine, left_join(unique(coordinate.HIDTA[!is.na(coordinate.HIDTA$HIDTA),c(7, 12, 14, 15)]), cocaine,
#                                       by=c("state_name", "county", "GEOID")) %>%
#                      select(state_name, county, GEOID, HIDTA.x, Jan_2018:Dec_2021) %>% rename(HIDTA=HIDTA.x) %>%
#                      filter(is.na(Jan_2018))) %>% arrange(GEOID)
#   names(cocaine)[1] <- "state"
#   cocaine[,5:52][is.na(cocaine[,5:52])] <- 0
# }
# 
# grep("Jan_2018", names(cocaine)) # 5
# grep("Dec_2021", names(cocaine)) # 52
# cocaine <- cocaine[(cocaine[,5:52] %>% apply(1, sum) > 0) | (!is.na(cocaine$HIDTA)),]
# cocaine[which(cocaine$state == "District of Columbia"), 3] <- 11001
# cocaine <- cocaine %>% filter(!is.na(GEOID) & !is.na(county)) %>% arrange(GEOID)
# cocaine # 1070 counties in total.
# write.csv(cocaine, "cocaine other count HIDTA (01-30-2023).csv", row.names=F)


# # crack weight
# {
# crack <- seizures %>% filter(Drug=="Crack" & !(State %in% c("Alaska", "Hawaii", "U.S. Virgin Islands", "France")) & !is.na(County) )
# calender <- data.frame(Month=1:12, Month2=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
# crack <- inner_join(crack, calender, by="Month")
# 
# crack <- crack %>% group_by(State, County, Year, Month2) %>% summarise(weight=sum(Quantity)) %>%
#   pivot_wider(names_from=c(Month2, Year), names_sep="_", values_from=weight, values_fill=0)
# cal_order <- c(str_c(calender$Month2, 2018, sep="_"),
#                str_c(calender$Month2, 2019, sep="_"),
#                str_c(calender$Month2, 2020, sep="_"),
#                str_c(calender$Month2, 2021, sep="_"))
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
# crack <- crack %>% filter(!is.na(GEOID)) %>% arrange(GEOID)
# crack # 813 counties in total.
# write.csv(crack, "cocaine crack weight HIDTA (01-30-2023).csv", row.names=F)

# # other cocaine weight
# {
#   cocaine <- seizures %>% filter(Drug %in% c("Cocaine", "Cocaine, Powder", "Cocaine Base",
#                                              "Cocaine precursor", "Coca Leaves", "Cocaine metabolite") &
#                                    !(State %in% c("Alaska", "Hawaii", "U.S. Virgin Islands", "France")) &
#                                    !is.na(County) )
#   cocaine <- inner_join(cocaine, calender, by="Month")
# 
#   cocaine <- cocaine %>% group_by(State, County, Year, Month2) %>% summarise(weight=sum(Quantity)) %>%
#     pivot_wider(names_from=c(Month2, Year), names_sep="_", values_from=weight, values_fill=0)
#   cocaine <- cocaine[, c(1:2, match(cal_order, names(cocaine)))]
#   names(cocaine)[1:2] <- c("state_name", "county")
#   cocaine <- left_join(cocaine, unique(coordinate.HIDTA[,c(7, 10, 12, 14, 15)]), by=c("state_name", "county")) %>%
#     select(state_name, county, GEOID, HIDTA, Jan_2018:Dec_2021)
#   cocaine <- rbind(cocaine, left_join(unique(coordinate.HIDTA[!is.na(coordinate.HIDTA$HIDTA),c(7, 12, 14, 15)]), cocaine,
#                                       by=c("state_name", "county", "GEOID")) %>%
#                      select(state_name, county, GEOID, HIDTA.x, Jan_2018:Dec_2021) %>% rename(HIDTA=HIDTA.x) %>%
#                      filter(is.na(Jan_2018))) %>% arrange(GEOID)
#   names(cocaine)[1] <- "state"
#   cocaine[,5:52][is.na(cocaine[,5:52])] <- 0
# }
# 
# grep("Jan_2018", names(cocaine)) # 5
# grep("Dec_2021", names(cocaine)) # 52
# cocaine <- cocaine[(cocaine[,5:52] %>% apply(1, sum) > 0) | (!is.na(cocaine$HIDTA)),]
# cocaine[which(cocaine$state == "District of Columbia"), 3] <- 11001
# cocaine <- cocaine %>% filter(!is.na(GEOID) & !is.na(county)) %>% arrange(GEOID)
# cocaine # 1070 counties in total.
# write.csv(cocaine, "cocaine other weight HIDTA (01-30-2023).csv", row.names=F)

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

