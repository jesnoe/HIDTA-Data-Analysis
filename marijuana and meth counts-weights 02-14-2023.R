library(fpp2)
library(spdep)
library(readxl)
library(urbnmapr)
library(tidyverse)
library(gridExtra)
# library(maps)
# library(rnaturalearth)

# {
#   adjacency <- read.csv("county_adjacency2010.csv")
#   adjacency$county <- substring(adjacency$countyname, 1, str_locate(adjacency$countyname, ",")[,1]-1)
#   adjacency$state1 <- substring(adjacency$countyname, str_locate(adjacency$countyname, ",")[,1]+2)
#   adjacency$neighbor_county <- substring(adjacency$neighborname, 1, str_locate(adjacency$neighborname, ",")[,1]-1)
#   adjacency$neighbor_state1 <- substring(adjacency$neighborname, str_locate(adjacency$neighborname, ",")[,1]+2)
#   
#   state_names <- read.csv("state names.csv")
#   names(state_names)[1] <- "state"
#   names(state_names)[3] <- "state1"
#   adjacency <- merge(adjacency, state_names[,c(1,3)], all.x=T)
#   
#   names(state_names)[1] <- "neighbor_state"
#   names(state_names)[3] <- "neighbor_state1"
#   adjacency <- merge(adjacency, state_names[,c(1,3)], all.x=T)
#   adjacency <- adjacency[, c(4, 9, 7, 6, 8, 10)]
#   head(adjacency)
#   names(adjacency)[1] <- "county_fips"
#   names(adjacency)[4] <- "neighbor_county_fips"
#   adjacency <- adjacency[adjacency$county_fips != adjacency$neighbor_county_fips,]
#   
#   county.fips <- unique(adjacency[,1:3])
#   populations <- read.csv("co-est2019-alldata.csv") %>% filter(COUNTY != 0) %>%
#     select(SUMLEV, REGION, STATE, COUNTY, STNAME, CTYNAME, POPESTIMATE2018, POPESTIMATE2019)
#   populations <- merge(populations,
#                        read.csv("co-est2021-alldata.csv") %>% filter(COUNTY != 0) %>% select(STATE, COUNTY, POPESTIMATE2020, POPESTIMATE2021),
#                        by=c("STATE", "COUNTY"))
#   
#   names(populations)[5] <- "state"
#   names(populations)[6] <- "county"
#   county.fips <- merge(county.fips, populations[,5:10], by=c("state", "county"), all.x=T) %>% as_tibble
#   names(county.fips)[3] <- "GEOID"
# }

overdose <- read.csv("VSRR_Provisional_County-Level_Drug_Overdose_Death_Counts.csv") %>% as_tibble

# urbnmapr data
unique(countydata$county_fips) # some ips have 5 digits. need to be converted into integer.
HIDTA.dist <- read.csv("HIDTA Regions.csv", header=T) %>% as_tibble
coordinate.HIDTA <- left_join(counties, HIDTA.dist, by=c("state_name", "county_name")) %>% filter(!(state_name %in% c("Alaska", "Hawaii", "Puerto Rico")))
names(coordinate.HIDTA)[7] <- "GEOID"

HIDTA.regions <- coordinate.HIDTA %>% filter(!is.na(HIDTA)) %>% select(GEOID, HIDTA) %>% unique

## Too small observations for non-meth cocaine
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

# meth count
{
meth <- seizures %>% filter(Drug %in% c("Methamphetamine", "Ice", "Amphetamine", "Methamphetamine in solution",
                                        "Meth precursor: Ephedrine hydrochloride", "Methamphetamine Oil", "Meth Precursor Chemicals",
                                        "Methamphetamine Tablet", "Ephedrine", "Pseudoephedrine") & !(State %in% c("Alaska", "Hawaii", "U.S. Virgin Islands", "France")) & !is.na(County) )
calender <- data.frame(Month=1:12, Month2=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
meth <- inner_join(meth, calender, by="Month")

meth <- meth %>% group_by(State, County, Year, Month2) %>% summarise(count=sum(Quantity > 0)) %>%
  pivot_wider(names_from=c(Month2, Year), names_sep="_", values_from=count, values_fill=0)
cal_order <- c(str_c(calender$Month2, 2018, sep="_"),
               str_c(calender$Month2, 2019, sep="_"),
               str_c(calender$Month2, 2020, sep="_"),
               str_c(calender$Month2, 2021, sep="_"))
meth <- meth[, c(1:2, match(cal_order, names(meth)))]
names(meth)[1:2] <- c("state_name", "county")
coordinate.HIDTA$county <- str_split(coordinate.HIDTA$county_name, " ") %>% lapply(function(x) str_c(c(x[(1:(length(x)-1))]), collapse=" ") ) %>% unlist
meth <- left_join(meth, unique(coordinate.HIDTA[,c(7, 12, 14, 15)]), by=c("state_name", "county")) %>%
  select(state_name, county, GEOID, HIDTA, Jan_2018:Dec_2021)
meth <- rbind(meth, left_join(unique(coordinate.HIDTA[!is.na(coordinate.HIDTA$HIDTA),c(7, 12, 14, 15)]), meth,
                                by=c("state_name", "county", "GEOID")) %>%
                       select(state_name, county, GEOID, HIDTA.x, Jan_2018:Dec_2021) %>% rename(HIDTA=HIDTA.x) %>%
                       filter(is.na(Jan_2018))) %>% arrange(GEOID)
names(meth)[1] <- "state"
meth[,5:52][is.na(meth[,5:52])] <- 0
}

grep("Jan_2018", names(meth)) # 5
grep("Dec_2021", names(meth)) # 52
meth <- meth[(meth[,5:52] %>% apply(1, sum) > 0) | (!is.na(meth$HIDTA)),]
meth[which(meth$state == "District of Columbia"), 3] <- 11001
meth <- meth %>% filter(!is.na(GEOID))
meth # 1318 counties in total.
# write.csv(meth, "meth count HIDTA (02-14-2023).csv", row.names=F)


# meth weight
{
meth <- seizures %>% filter(Drug %in% c("Methamphetamine", "Ice", "Amphetamine", "Methamphetamine in solution",
                                          "Meth precursor: Ephedrine hydrochloride", "Methamphetamine Oil", "Meth Precursor Chemicals",
                                          "Methamphetamine Tablet", "Ephedrine", "Pseudoephedrine") & !(State %in% c("Alaska", "Hawaii", "U.S. Virgin Islands", "France")) & !is.na(County) )
calender <- data.frame(Month=1:12, Month2=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
meth <- inner_join(meth, calender, by="Month")

meth <- meth %>% group_by(State, County, Year, Month2) %>% summarise(weight=sum(Quantity)) %>%
  pivot_wider(names_from=c(Month2, Year), names_sep="_", values_from=weight, values_fill=0)
cal_order <- c(str_c(calender$Month2, 2018, sep="_"),
               str_c(calender$Month2, 2019, sep="_"),
               str_c(calender$Month2, 2020, sep="_"),
               str_c(calender$Month2, 2021, sep="_"))
meth <- meth[, c(1:2, match(cal_order, names(meth)))]
names(meth)[1:2] <- c("state_name", "county")
coordinate.HIDTA$county <- str_split(coordinate.HIDTA$county_name, " ") %>% lapply(function(x) str_c(c(x[(1:(length(x)-1))]), collapse=" ") ) %>% unlist
meth <- left_join(meth, unique(coordinate.HIDTA[,c(7, 12, 14, 15)]), by=c("state_name", "county")) %>%
  select(state_name, county, GEOID, HIDTA, Jan_2018:Dec_2021)
meth <- rbind(meth, left_join(unique(coordinate.HIDTA[!is.na(coordinate.HIDTA$HIDTA),c(7, 12, 14, 15)]), meth,
                                by=c("state_name", "county", "GEOID")) %>%
                       select(state_name, county, GEOID, HIDTA.x, Jan_2018:Dec_2021) %>% rename(HIDTA=HIDTA.x) %>%
                       filter(is.na(Jan_2018))) %>% arrange(GEOID)
names(meth)[1] <- "state"
meth[,5:52][is.na(meth[,5:52])] <- 0
}

grep("Jan_2018", names(meth)) # 5
grep("Dec_2021", names(meth)) # 52
meth <- meth[(meth[,5:52] %>% apply(1, sum) > 0) | (!is.na(meth$HIDTA)),]
meth[which(meth$state == "District of Columbia"), 3] <- 11001
meth <- meth %>% filter(!is.na(GEOID)) %>% arrange(GEOID)
meth # 1318 counties in total.
# write.csv(meth, "meth weight HIDTA (02-14-2023).csv", row.names=F)

# marijuana count
{
marijuana <- seizures %>% filter(grepl("Marijuana", Drug) & !(State %in% c("Alaska", "Hawaii", "U.S. Virgin Islands", "France")) & !is.na(County) )
calender <- data.frame(Month=1:12, Month2=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
marijuana <- inner_join(marijuana, calender, by="Month")

marijuana <- marijuana %>% group_by(State, County, Year, Month2) %>% summarise(count=sum(Quantity > 0)) %>%
  pivot_wider(names_from=c(Month2, Year), names_sep="_", values_from=count, values_fill=0)
cal_order <- c(str_c(calender$Month2, 2018, sep="_"),
               str_c(calender$Month2, 2019, sep="_"),
               str_c(calender$Month2, 2020, sep="_"),
               str_c(calender$Month2, 2021, sep="_"))
marijuana <- marijuana[, c(1:2, match(cal_order, names(marijuana)))]
names(marijuana)[1:2] <- c("state_name", "county")
coordinate.HIDTA$county <- str_split(coordinate.HIDTA$county_name, " ") %>% lapply(function(x) str_c(c(x[(1:(length(x)-1))]), collapse=" ") ) %>% unlist
marijuana <- left_join(marijuana, unique(coordinate.HIDTA[,c(7, 12, 14, 15)]), by=c("state_name", "county")) %>%
  select(state_name, county, GEOID, HIDTA, Jan_2018:Dec_2021)
marijuana <- rbind(marijuana, left_join(unique(coordinate.HIDTA[!is.na(coordinate.HIDTA$HIDTA),c(7, 12, 14, 15)]), marijuana,
                                        by=c("state_name", "county", "GEOID")) %>%
                     select(state_name, county, GEOID, HIDTA.x, Jan_2018:Dec_2021) %>% rename(HIDTA=HIDTA.x) %>%
                     filter(is.na(Jan_2018))) %>% arrange(GEOID)
names(marijuana)[1] <- "state"
marijuana[,5:52][is.na(marijuana[,5:52])] <- 0
}

grep("Jan_2018", names(marijuana)) # 5
grep("Dec_2021", names(marijuana)) # 52
marijuana <- marijuana[(marijuana[,5:52] %>% apply(1, sum) > 0) | (!is.na(marijuana$HIDTA)),]
marijuana[which(marijuana$state == "District of Columbia"), 3] <- 11001
marijuana <- marijuana %>% filter(!is.na(GEOID))
marijuana # 1218 counties in total.
# write.csv(marijuana, "marijuana count HIDTA (02-14-2023).csv", row.names=F)


# marijuana weight
{
marijuana <- seizures %>% filter(grepl("Marijuana", Drug) & !(State %in% c("Alaska", "Hawaii", "U.S. Virgin Islands", "France")) & !is.na(County) )
calender <- data.frame(Month=1:12, Month2=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
marijuana <- inner_join(marijuana, calender, by="Month")

marijuana <- marijuana %>% group_by(State, County, Year, Month2) %>% summarise(weight=sum(Quantity)) %>%
  pivot_wider(names_from=c(Month2, Year), names_sep="_", values_from=weight, values_fill=0)
cal_order <- c(str_c(calender$Month2, 2018, sep="_"),
               str_c(calender$Month2, 2019, sep="_"),
               str_c(calender$Month2, 2020, sep="_"),
               str_c(calender$Month2, 2021, sep="_"))
marijuana <- marijuana[, c(1:2, match(cal_order, names(marijuana)))]
names(marijuana)[1:2] <- c("state_name", "county")
coordinate.HIDTA$county <- str_split(coordinate.HIDTA$county_name, " ") %>% lapply(function(x) str_c(c(x[(1:(length(x)-1))]), collapse=" ") ) %>% unlist
marijuana <- left_join(marijuana, unique(coordinate.HIDTA[,c(7, 12, 14, 15)]), by=c("state_name", "county")) %>%
  select(state_name, county, GEOID, HIDTA, Jan_2018:Dec_2021)
marijuana <- rbind(marijuana, left_join(unique(coordinate.HIDTA[!is.na(coordinate.HIDTA$HIDTA),c(7, 12, 14, 15)]), marijuana,
                                        by=c("state_name", "county", "GEOID")) %>%
                     select(state_name, county, GEOID, HIDTA.x, Jan_2018:Dec_2021) %>% rename(HIDTA=HIDTA.x) %>%
                     filter(is.na(Jan_2018))) %>% arrange(GEOID)
names(marijuana)[1] <- "state"
marijuana[,5:52][is.na(marijuana[,5:52])] <- 0
}

grep("Jan_2018", names(marijuana)) # 5
grep("Dec_2021", names(marijuana)) # 52
marijuana <- marijuana[(marijuana[,5:52] %>% apply(1, sum) > 0) | (!is.na(marijuana$HIDTA)),]
marijuana[which(marijuana$state == "District of Columbia"), 3] <- 11001
marijuana <- marijuana %>% filter(!is.na(GEOID)) %>% arrange(GEOID)
marijuana # 1308 counties in total.
# write.csv(marijuana, "marijuana weight HIDTA (02-14-2023).csv", row.names=F)
