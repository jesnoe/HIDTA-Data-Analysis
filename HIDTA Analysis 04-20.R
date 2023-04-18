library(readxl)
library(tidyverse)
library(urbnmapr)
library(gridExtra)

adjacency <- read.csv("county_adjacency2010.csv")
adjacency$county <- substring(adjacency$countyname, 1, str_locate(adjacency$countyname, " ")[,1]-1)
adjacency$state1 <- substring(adjacency$countyname, str_locate(adjacency$countyname, ",")[,1]+2)
adjacency$neighbor_county <- substring(adjacency$neighborname, 1, str_locate(adjacency$neighborname, " ")[,1]-1)
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
populations <- read.csv("co-est2019-alldata.csv") %>% filter(COUNTY != 0) %>% select(SUMLEV, REGION, STATE, COUNTY, STNAME, CTYNAME, POPESTIMATE2019)
populations$county <- substring(populations$CTYNAME, 1, str_locate(populations$CTYNAME, " ")[,1]-1)
names(populations)[5] <- "state"
names(populations)[7] <- "population"
county.fips <- merge(county.fips, populations[,c(5,7,8)], by=c("state", "county"), all.x=T)


# urbnmapr data
unique(countydata$county_fips) # some ips have 5 digits. need to be converted into integer.

# HIDTA Data
seizures <- read_xlsx("Drug Seizures All HIDTAs All Drugs 2018-2021 Combined.xlsx")
seizures$County <- substring(seizures$County, 1, str_locate(seizures$County, ",")[,1]-1)
seizures$Year <- substring(seizures$SeizureDate, 1, 4)
seizures$Month <- substring(seizures$SeizureDate, 6, 7)
seizures$Year <- as.factor(seizures$Year)
seizures$Month <- as.factor(seizures$Month)
cocaine <- seizures %>% filter(Drug=="Cocaine")
cocaine %>% group_by(State) %>% summarise(log.obs=n()) %>% arrange(desc(log.obs)) %>% as.data.frame


cocaine.FL <- cocaine %>% filter(State=="Florida")
cocaine.AL <- cocaine %>% filter(State=="Alabama")
cocaine.GA <- cocaine %>% filter(State=="Georgia")
adjacency.FL <- adjacency %>% filter(state=="Florida") %>% arrange(county)
unique(cocaine.FL$County)
unique(cocaine.FL$County) %>% length

unique(adjacency.FL$county)
unique(adjacency.FL$county)[!(unique(adjacency.FL$county) %in% unique(cocaine.FL$County))]
adjacency.FL %>% filter(county %in% c("Escambia", "Flagler", "Franklin", "Gadsden")) %>% select(state, county, neighbor_county, neighbor_state)

# Florida map
(cocaine.FL.county <- cocaine.FL %>% group_by(County) %>% summarise(log.obs=log(n()), avg.seizure=mean(Quantity), total.seizure=sum(Quantity), n.obs=n()) %>% arrange(desc(log.obs)))
cocaine.FL.county %>% as.data.frame
names(cocaine.FL.county)[1] <- "county"
FL.data <- left_join(unique(adjacency.FL[,c(1,3)]), cocaine.FL.county)
FL.population <- (county.fips %>% filter(state=="Florida"))[,3:4]
FL.data <- left_join(FL.data, FL.population)

counties$county_fips <- as.integer(counties$county_fips)
counties %>% 
  filter(state_name =="Florida") %>% 
  left_join(FL.data, by = "county_fips") %>% 
  ggplot(mapping = aes(long, lat, group = group, fill = log.obs)) +
  geom_polygon(color = "#ffffff", size = .25) +
  labs(fill = "Log of Obs.")

counties %>% 
  filter(state_name =="Florida") %>% 
  left_join(FL.data, by = "county_fips") %>% 
  ggplot(mapping = aes(long, lat, group = group, fill = log.obs/population)) +
  geom_polygon(color = "#ffffff", size = .25) +
  labs(fill = "Log of Obs.\n/Population")

counties %>% 
  filter(state_name =="Florida") %>% 
  left_join(FL.data, by = "county_fips") %>% 
  ggplot(mapping = aes(long, lat, group = group, fill = log.obs/log(population))) +
  geom_polygon(color = "#ffffff", size = .25) +
  labs(fill = "Log of Obs.\n/Log(Population)")


counties %>% 
  filter(state_name =="Florida") %>% 
  left_join(FL.data, by = "county_fips") %>% 
  ggplot(mapping = aes(long, lat, group = group, fill = n.obs/population)) +
  geom_polygon(color = "#ffffff", size = .25) +
  labs(fill = "# of Obs.\n/Population")

counties %>% 
  filter(state_name =="Florida") %>% 
  left_join(FL.data, by = "county_fips") %>% 
  ggplot(mapping = aes(long, lat, group = group, fill = log.obs!=0)) +
  geom_polygon(color = "#ffffff", size = .25) +
  labs(fill = "Log of Obs.")

log.obs.FL <- counties[,c(7, 12)] %>% 
  filter(state_name =="Florida") %>% 
  left_join(FL.data, by = "county_fips") %>% unique
log.obs.FL[is.na(log.obs.FL$log.obs),] %>% arrange(state_name) %>%  as.data.frame


# Alabama and Georgia map

(cocaine.AL.county <- cocaine.AL %>% group_by(County) %>% summarise(log.obs=log(n()), avg.seizure=mean(Quantity), total.seizure=sum(Quantity)) %>% arrange(desc(log.obs)))
cocaine.AL.county %>% as.data.frame
names(cocaine.AL.county)[1] <- "county"
adjacency.AL <- adjacency %>% filter(state=="Alabama") %>% arrange(county)
AL.data <- left_join(unique(adjacency.AL[,c(1,3)]), cocaine.AL.county)
counties %>% 
  filter(state_name =="Alabama") %>% 
  left_join(AL.data, by = "county_fips") %>% 
  ggplot(mapping = aes(long, lat, group = group, fill = log.obs)) +
  geom_polygon(color = "#ffffff", size = .25) +
  labs(fill = "Log of Obs.")

(cocaine.GA.county <- cocaine.GA %>% group_by(County) %>% summarise(log.obs=log(n()), avg.seizure=mean(Quantity), total.seizure=sum(Quantity)) %>% arrange(desc(log.obs)))
cocaine.GA.county %>% as.data.frame
names(cocaine.GA.county)[1] <- "county"
adjacency.GA <- adjacency %>% filter(state=="Georgia") %>% arrange(county)
GA.data <- left_join(unique(adjacency.GA[,c(1,3)]), cocaine.GA.county)
counties %>% 
  filter(state_name =="Georgia") %>% 
  left_join(GA.data, by = "county_fips") %>% 
  ggplot(mapping = aes(long, lat, group = group, fill = log.obs)) +
  geom_polygon(color = "#ffffff", size = .25) +
  labs(fill = "Log of Obs.")


# LISA Plot
cocaine.FL.monthly <- cocaine.FL %>% group_by(State, County, Year, Month) %>% summarise(log.obs=log(n()), avg.seizure=mean(Quantity), total.seizure=sum(Quantity))
cocaine.FL.monthly$Month <- as.numeric(cocaine.FL.monthly$Month)

FL.counties <- unique(cocaine.FL.monthly$County)
names(cocaine.FL.monthly)[1] <- "state"
names(cocaine.FL.monthly)[2] <- "county"
temp <- cocaine.FL.monthly
cocaine.FL.monthly <- temp[1,]
for (county in FL.counties) {
  for (year in c("2018", "2019", "2020", "2021")) {
    dummy <- data.frame(state=rep("Florida", 12), county=rep(county, 12), Year=rep(year, 12), Month=1:12)
    county.year <- temp %>% filter(county==county, Year==year)
    cocaine.FL.monthly <- rbind(cocaine.FL.monthly, merge(dummy, county.year, all.x=T))
  }
}
cocaine.FL.monthly <- cocaine.FL.monthly[-1,]
cocaine.FL.monthly <- left_join(cocaine.FL.monthly, county.fips, all.x=T)


cocaine.AL.monthly <- cocaine.AL %>% group_by(State, County, Year, Month) %>% summarise(log.obs=log(n()), avg.seizure=mean(Quantity), total.seizure=sum(Quantity))
cocaine.AL.monthly$Month <- as.numeric(cocaine.AL.monthly$Month)

AL.counties <- unique(cocaine.AL.monthly$County)
names(cocaine.AL.monthly)[1] <- "state"
names(cocaine.AL.monthly)[2] <- "county"
temp <- cocaine.AL.monthly
cocaine.AL.monthly <- temp[1,]
for (county in AL.counties) {
  for (year in c("2018", "2019", "2020", "2021")) {
    dummy <- data.frame(state=rep("Alabama", 12), county=rep(county, 12), Year=rep(year, 12), Month=1:12)
    county.year <- temp %>% filter(county==county, Year==year)
    cocaine.AL.monthly <- rbind(cocaine.AL.monthly, merge(dummy, county.year, all.x=T))
  }
}

cocaine.AL.monthly <- cocaine.AL.monthly[-1,]
cocaine.AL.monthly <- left_join(cocaine.AL.monthly, county.fips, all.x=T)
AL.county <- adjacency.AL %>% filter(neighbor_state == "Florida") %>% select(county) %>% unique
cocaine.AL.monthly <- cocaine.AL.monthly %>% filter(county %in% AL.county$county)


cocaine.GA.monthly <- cocaine.GA %>% group_by(State, County, Year, Month) %>% summarise(log.obs=log(n()), avg.seizure=mean(Quantity), total.seizure=sum(Quantity))
cocaine.GA.monthly$Month <- as.numeric(cocaine.GA.monthly$Month)

GA.counties <- unique(cocaine.GA.monthly$County)
names(cocaine.GA.monthly)[1] <- "state"
names(cocaine.GA.monthly)[2] <- "county"
temp <- cocaine.GA.monthly
cocaine.GA.monthly <- temp[1,]
for (county in GA.counties) {
  for (year in c("2018", "2019", "2020", "2021")) {
    dummy <- data.frame(state=rep("Georgia", 12), county=rep(county, 12), Year=rep(year, 12), Month=1:12)
    county.year <- temp %>% filter(county==county, Year==year)
    cocaine.GA.monthly <- rbind(cocaine.GA.monthly, merge(dummy, county.year, all.x=T))
  }
}

cocaine.GA.monthly <- cocaine.GA.monthly[-1,]
cocaine.GA.monthly <- left_join(cocaine.GA.monthly, county.fips, all.x=T)
GA.county <- adjacency.GA %>% filter(neighbor_state == "Florida") %>% select(county) %>% unique
cocaine.GA.monthly <- cocaine.GA.monthly %>% filter(county %in% GA.county$county)

n.FL <- nrow(cocaine.FL.monthly) # 1824
n.AL <- nrow(cocaine.AL.monthly) # 1344
n.GA <- nrow(cocaine.GA.monthly) # 2112
cocaine.FL.monthly[is.na(cocaine.FL.monthly)] <- 0
cocaine.AL.monthly[is.na(cocaine.AL.monthly)] <- 0
cocaine.GA.monthly[is.na(cocaine.GA.monthly)] <- 0
log.obs.neigbors.std <- scale(c(cocaine.FL.monthly$log.obs, cocaine.AL.monthly$log.obs, cocaine.GA.monthly$log.obs)) %>% as.vector

cocaine.FL.monthly$log.obs.std <- scale(cocaine.FL.monthly$log.obs) %>% as.vector
cocaine.FL.monthly$log.obs.neigbors.std <- log.obs.neigbors.std[1:n.FL]
cocaine.AL.monthly$log.obs.neigbors.std <- log.obs.neigbors.std[(n.FL+1):(n.FL+n.AL)]
cocaine.GA.monthly$log.obs.neigbors.std <- log.obs.neigbors.std[(n.FL+n.AL+1):(n.FL+n.AL+n.GA)]

cocaine.FL.monthly %>% select(state, county, Year, Month, log.obs, log.obs.std, log.obs.neigbors.std)
rbind(head(cocaine.FL.monthly),
      head(cocaine.AL.monthly),
      head(cocaine.GA.monthly)) %>% select(state, county, Year, Month, log.obs)

FL.neighbors <- list()
for (county in unique(adjacency.FL$county)) {
  neighbor_county_fips <- adjacency.FL$neighbor_county_fips[adjacency.FL$county == county]

  FL.neighbors[[county]] <- rbind(cocaine.FL.monthly %>% filter(county_fips %in% neighbor_county_fips),
                                  cocaine.AL.monthly %>% filter(county_fips %in% neighbor_county_fips),
                                  cocaine.GA.monthly %>% filter(county_fips %in% neighbor_county_fips))
}
head(adjacency.FL)
head(FL.neighbors$Alachua, 12)
head(FL.neighbors$Washington, 12)




lisa <- function(row) {
  county <- row[[2]]
  year <- row[[3]]
  month <- row[[4]]
  if ( is.null(FL.neighbors[[county]]) ) {return(NA)}
  neighbors.df <- FL.neighbors[[county]] %>% filter(Year == year, Month == month)
  log.obs <- neighbors.df$log.obs.neigbors.std
  result <- mean(log.obs, na.rm=T)
  if (is.nan(result)) {return(NA)}
  return(result)
}

LISA.neighbors <- c()
for (i in 1:nrow(cocaine.FL.monthly)) {
  LISA.neighbors <- c(LISA.neighbors, lisa(cocaine.FL.monthly[i,]))
}
LISA.neighbors

cocaine.FL.monthly$LISA.neighbors <- LISA.neighbors
cocaine.FL.monthly

min(cocaine.FL.monthly$log.obs.std, na.rm=T); max(cocaine.FL.monthly$log.obs.std, na.rm=T) # c(-1.546014, 3.482996)
min(cocaine.FL.monthly$LISA.neighbors, na.rm=T); max(cocaine.FL.monthly$LISA.neighbors, na.rm=T) # c(-1.261443, 2.685786)

par(mfrow=c(1,2))
hist(cocaine.FL.monthly$log.obs.std, xlab="Local Log of Obs.", main=NULL)
hist(cocaine.FL.monthly$LISA.neighbors, xlab="Neighbor Log of Obs.", main=NULL)
par(mfrow=c(1,1))

year <- "2018"; month <- 2
month.data <- cocaine.FL.monthly %>% filter(Year == year, Month == month)
plot(month.data$log.obs.std, month.data$LISA.neighbors,
     main=paste("Local and Neighbor Std Log of Obs. in ", year),
     xlab="Local Std Log of Obs.",
     ylab="Neighbor Std Log of Obs.",
     xlim=c(-1, 4), ylim=c(-1, 4), pch=19)

{
par(mfrow=c(2,2))
for (month in 1:4) {
  month.data <- cocaine.FL.monthly %>% filter(Year == year, Month == month)
  plot(month.data$log.obs.std, month.data$LISA.neighbors,
       main=paste("Local and Neighbor Std Log of Obs. in ", year),
       xlab="Local Std Log of Obs.",
       ylab="Neighbor Std Log of Obs.",
       xlim=c(-1, 4), ylim=c(-1, 4), pch=19)
}
par(mfrow=c(1,1))
}

cocaine.FL.monthly.total <- cocaine.FL %>% group_by(Year, Month) %>% summarise(log.obs=log(n()))
cocaine.FL.monthly.total$Year.Month <- as.numeric(as.character(cocaine.FL.monthly.total$Year)) + 
                                       as.numeric(as.character(cocaine.FL.monthly.total$Month))/12
plot(cocaine.FL.monthly.total$Year.Month,
     cocaine.FL.monthly.total$log.obs,
     xlab="Year", ylab="Log of Obs.", main="Total Log of Obs.", pch=19, type="b")
year2 <- "2020"
windows(width=9,  height=8)
{
  par(mfrow=c(2,2))
  for (month in 1:4) {
    month.data <- cocaine.FL.monthly %>% filter(Year == year, Month == month)
    month.data2 <- cocaine.FL.monthly %>% filter(Year == year2, Month == month)
    year1.nobs <- nrow(month.data)
    plot(month.data$log.obs.std, month.data$LISA.neighbors,
         main=paste(month, "-th Month in", year, "and", year2),
         xlab="Local Std Log of Obs.",
         ylab="Neighbor Std Log of Obs.",
         xlim=c(-2, 4), ylim=c(-2, 4), pch=19)
    lines(rep(0,year1.nobs), seq(-2, 4, len=year1.nobs))
    lines(seq(-2, 4, len=year1.nobs), rep(0,year1.nobs))
    points(month.data2$log.obs.std, month.data2$LISA.neighbors, col=2, pch=19)
  }
  legend("topright", c("2018", "2020"), col=1:2, pch=19, cex=1, text.width=0.5)
  par(mfrow=c(1,1))
}
