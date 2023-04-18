library(fpp2)
library(readxl)
library(urbnmapr)
library(tidyverse)

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

avg.seizure <- read.csv("All Drugs Average Seizures.csv")
names(avg.seizure)[2] <- "avg.seizure.all"
avg.seizure <- left_join(avg.seizure, read.csv("Cocaine Average Seizures.csv"), by="state")
names(avg.seizure)[3] <- "avg.seizure.cocaine"

state.centroid <- counties %>% group_by(state_name) %>% summarise(long=mean(long), lat=mean(lat))
names(state.centroid)[1] <- "state"
avg.seizure <- left_join(avg.seizure, state.centroid, by="state")
n_state <- nrow(avg.seizure)

avg.seizure$dist_from_FL <- numeric(n_state)
avg.seizure$dist_from_AZ <- numeric(n_state)
avg.seizure$dist_from_CA <- numeric(n_state)
avg.seizure$dist_from_NM <- numeric(n_state)
avg.seizure$dist_from_TX <- numeric(n_state)
for (i in 1:n_state) {
  avg.seizure$dist_from_FL[i] <- sqrt(sum((avg.seizure[i,] %>% select(long, lat) - avg.seizure[7,] %>% select(long, lat))^2))
  avg.seizure$dist_from_AZ[i] <- sqrt(sum((avg.seizure[i,] %>% select(long, lat) - avg.seizure[6,] %>% select(long, lat))^2))
  avg.seizure$dist_from_CA[i] <- sqrt(sum((avg.seizure[i,] %>% select(long, lat) - avg.seizure[1,] %>% select(long, lat))^2))
  avg.seizure$dist_from_NM[i] <- sqrt(sum((avg.seizure[i,] %>% select(long, lat) - avg.seizure[18,] %>% select(long, lat))^2))
  avg.seizure$dist_from_TX[i] <- sqrt(sum((avg.seizure[i,] %>% select(long, lat) - avg.seizure[4,] %>% select(long, lat))^2))
}

  # Southern border is the minimum distances from Arizona, California, New Mexico, and Texas
avg.seizure$dist_from_SB <- numeric(n_state)
Southern.Border.centroid <- avg.seizure %>% filter(state %in% c("Arizona", "California", "New Mexico", "Texas")) %>% select(long, lat)
names(Southern.Border.centroid) <- NULL
for (i in 1:n_state) {
  if (avg.seizure$state[i] %in% c("Arizona", "California", "New Mexico", "Texas")) next
  coordinate.i <- avg.seizure[i,] %>% select(long, lat) %>% as.numeric
  avg.seizure$dist_from_SB[i] <- sqrt(
    apply((matrix(coordinate.i, 4,2, byrow=T) - as.matrix(Southern.Border.centroid, 4, 2))^2, 1, sum)
    ) %>% min
}

avg.seizure
state.populations <- populations %>% group_by(state) %>% summarise(population=sum(POPESTIMATE2018))
avg.seizure <- left_join(avg.seizure, state.populations, by="state")
# interstate <- read.csv("interstate.csv") %>% as_tibble
# interstate$interstate <- interstate$RURAL.INTERSTATE + interstate$URBAN.INTERSTATE
# names(interstate)[1] <- "state"
# avg.seizure <- left_join(avg.seizure, interstate[,c(1,17)], by="state")


## 1. Average seizures vs. populations, locations, ...
pairs(avg.seizure[,-c(1,3:4)])
lm(avg.seizure.all~population+dist_from_SB+dist_from_FL, data=avg.seizure) %>% summary
lm(avg.seizure.cocaine~population+dist_from_SB+dist_from_FL, data=avg.seizure) %>% summary
  # distance from Florida is significant for all drugs, not for cocaine.

lm(avg.seizure.all~population+dist_from_FL+dist_from_AZ+dist_from_CA+dist_from_NM+dist_from_TX,
   data=avg.seizure) %>% summary
lm(avg.seizure.cocaine~population+dist_from_FL+dist_from_AZ+dist_from_CA+dist_from_NM+dist_from_TX,
   data=avg.seizure) %>% summary
  # But Florida becomes insignificant when regressed on individual southern border states. Instead California becomes significant

lm(avg.seizure.all~population+dist_from_FL+dist_from_AZ+dist_from_CA+dist_from_NM+dist_from_TX,
   data=avg.seizure) %>% summary
lm(avg.seizure.cocaine~population+dist_from_FL+dist_from_AZ+dist_from_CA+dist_from_NM+dist_from_TX,
   data=avg.seizure) %>% summary

reg.pop <- lm(avg.seizure.all~population, data=avg.seizure[-1,])
reg.pop.res <- data.frame(state=avg.seizure[-1,1], residuals=reg.pop$residuals)

reg.pop.coc <- lm(avg.seizure.cocaine~population, data=avg.seizure[-1,])
reg.pop.res.coc <- data.frame(state=avg.seizure[-1,1], residuals=reg.pop.coc$residuals)

MainStates <- map_data("state")
names(MainStates)[5] <- "state"
reg.pop.res$state <- str_to_lower(reg.pop.res$state)
reg.pop.res <- left_join(MainStates, reg.pop.res, by="state")
ggplot(mapping = aes(long, lat, group = state, fill=residuals), data=reg.pop.res) +
  geom_polygon(color = "#000000", size = .05)



# Regression without southern border states including Florida
lm(avg.seizure.all~population+dist_from_SB+dist_from_FL, data=avg.seizure[-c(1,4,6,7,19),]) %>% summary
lm(avg.seizure.cocaine~population+dist_from_SB+dist_from_FL, data=avg.seizure[-c(1,4,6,7,19),]) %>% summary

lm(avg.seizure.all~population+dist_from_FL+dist_from_AZ+dist_from_CA+dist_from_NM+dist_from_TX,
   data=avg.seizure[-c(1,4,6,7,19),]) %>% summary
lm(avg.seizure.cocaine~population+dist_from_FL+dist_from_AZ+dist_from_CA+dist_from_NM+dist_from_TX,
   data=avg.seizure[-c(1,4,6,7,19),]) %>% summary


# Regression on a single state
lm(avg.seizure.all~population+dist_from_FL, data=avg.seizure) %>% summary
lm(avg.seizure.cocaine~population+dist_from_FL, data=avg.seizure) %>% summary
  # Florida is significant again for all drugs. For cocaine, p value of Florida is close to 0.05

lm(avg.seizure.all~population+dist_from_AZ, data=avg.seizure) %>% summary
lm(avg.seizure.cocaine~population+dist_from_AZ, data=avg.seizure) %>% summary
lm(avg.seizure.all~population+dist_from_CA, data=avg.seizure) %>% summary
lm(avg.seizure.cocaine~population+dist_from_CA, data=avg.seizure) %>% summary
lm(avg.seizure.all~population+dist_from_NM, data=avg.seizure) %>% summary
lm(avg.seizure.cocaine~population+dist_from_NM, data=avg.seizure) %>% summary
lm(avg.seizure.all~population+dist_from_TX, data=avg.seizure) %>% summary
lm(avg.seizure.cocaine~population+dist_from_TX, data=avg.seizure) %>% summary
lm(avg.seizure.all~population+dist_from_SB, data=avg.seizure) %>% summary
lm(avg.seizure.cocaine~population+dist_from_SB, data=avg.seizure) %>% summary
  # Only California is significant

# Regression on coordinates
lm(avg.seizure.all~population+long+lat+long*lat, data=avg.seizure) %>% summary
lm(avg.seizure.cocaine~population+long+lat+long*lat, data=avg.seizure) %>% summary
  # For all drugs, average seizure weights increases as a state is close to west and north

# Average seizures for each year
LISA.all <- read.csv("AllKGLISAResults.csv") %>% as_tibble
state.fips <- read.csv("us-state-ansi-fips.csv")
names(state.fips)[1] <- "state"
names(state.fips)[2] <- "STATEFP"
LISA.all <- merge(LISA.all, state.fips[,1:2], by="STATEFP") %>% as_tibble
LISA.all <- LISA.all[,c(1, 212, 2:211)]
LISA.all <- LISA.all[LISA.all[,21:68] %>% apply(1, sum) > 0,]

LISA.cocaine <- read.csv("CountyKNN3.csv") %>% as_tibble
LISA.cocaine <- merge(LISA.cocaine, state.fips[,1:2], by="STATEFP") %>% as_tibble
LISA.cocaine <- LISA.cocaine[,c(1, 215, 2:214)]

avg.seizure.per.state <- function(seizures.state, year, drug) {
  monthly.seizures <- data.frame()
  for (state.name in unique(seizures.state$state)) {
    seizures.i <- seizures.state %>% filter(state == state.name)
    monthly.seizures <- rbind(monthly.seizures, c(state.name, apply(seizures.i[,-(1:3)], 2, sum)))
  }
  names(monthly.seizures) <- names(seizures.state)[-(2:3)]
  monthly.seizures[,-1] <- apply(monthly.seizures[,-1], 2, as.numeric)
  monthly.seizures[[paste("avg.seizure", drug, year, sep=".")]] <- apply(monthly.seizures[,-1], 1, mean)
  return(monthly.seizures[,c(1,ncol(monthly.seizures))])
}


LISA.all.2018 <- LISA.all[,c(2,6,8,21:32)]
LISA.all.2019 <- LISA.all[,c(2,6,8,33:44)]
LISA.all.2020 <- LISA.all[,c(2,6,8,45:56)]
LISA.all.2021 <- LISA.all[,c(2,6,8,57:68)]

LISA.cocaine.2018 <- LISA.cocaine[,c(2,6,8,21:32)]
LISA.cocaine.2019 <- LISA.cocaine[,c(2,6,8,33:44)]
LISA.cocaine.2020 <- LISA.cocaine[,c(2,6,8,45:56)]
LISA.cocaine.2021 <- LISA.cocaine[,c(2,6,8,57:68)]

avg.seizure <- left_join(avg.seizure, avg.seizure.per.state(LISA.all.2018, 2018, "all"), by="state") %>% 
  left_join(avg.seizure.per.state(LISA.all.2019, 2019, "all"), by="state") %>% 
  left_join(avg.seizure.per.state(LISA.all.2020, 2020, "all"), by="state") %>% 
  left_join(avg.seizure.per.state(LISA.all.2021, 2021, "all"), by="state") %>% 
  left_join(avg.seizure.per.state(LISA.cocaine.2018, 2018, "cocaine"), by="state") %>% 
  left_join(avg.seizure.per.state(LISA.cocaine.2019, 2019, "cocaine"), by="state") %>% 
  left_join(avg.seizure.per.state(LISA.cocaine.2020, 2020, "cocaine"), by="state") %>% 
  left_join(avg.seizure.per.state(LISA.cocaine.2021, 2021, "cocaine"), by="state")
avg.seizure

lm(avg.seizure.all.2018~population+dist_from_SB+dist_from_FL, data=avg.seizure) %>% summary
lm(avg.seizure.all.2019~population+dist_from_SB+dist_from_FL, data=avg.seizure) %>% summary
lm(avg.seizure.all.2020~population+dist_from_SB+dist_from_FL, data=avg.seizure) %>% summary
lm(avg.seizure.all.2021~population+dist_from_SB+dist_from_FL, data=avg.seizure) %>% summary
lm(avg.seizure.cocaine.2018~population+dist_from_SB+dist_from_FL, data=avg.seizure) %>% summary
lm(avg.seizure.cocaine.2019~population+dist_from_SB+dist_from_FL, data=avg.seizure) %>% summary
lm(avg.seizure.cocaine.2020~population+dist_from_SB+dist_from_FL, data=avg.seizure) %>% summary
lm(avg.seizure.cocaine.2021~population+dist_from_SB+dist_from_FL, data=avg.seizure) %>% summary
  # similar to 4-year data

lm(avg.seizure.all.2018~population+dist_from_FL, data=avg.seizure) %>% summary
lm(avg.seizure.all.2019~population+dist_from_FL, data=avg.seizure) %>% summary
lm(avg.seizure.all.2020~population+dist_from_FL, data=avg.seizure) %>% summary
lm(avg.seizure.all.2021~population+dist_from_FL, data=avg.seizure) %>% summary
lm(avg.seizure.cocaine.2018~population+dist_from_FL, data=avg.seizure) %>% summary
lm(avg.seizure.cocaine.2019~population+dist_from_FL, data=avg.seizure) %>% summary
lm(avg.seizure.cocaine.2020~population+dist_from_FL, data=avg.seizure) %>% summary
lm(avg.seizure.cocaine.2021~population+dist_from_FL, data=avg.seizure) %>% summary
  # For cocaine, Florida is significant only in 2020

lm(avg.seizure.all.2018~population+dist_from_SB, data=avg.seizure) %>% summary
lm(avg.seizure.all.2019~population+dist_from_SB, data=avg.seizure) %>% summary
lm(avg.seizure.all.2020~population+dist_from_SB, data=avg.seizure) %>% summary
lm(avg.seizure.all.2021~population+dist_from_SB, data=avg.seizure) %>% summary
lm(avg.seizure.cocaine.2018~population+dist_from_SB, data=avg.seizure) %>% summary
lm(avg.seizure.cocaine.2019~population+dist_from_SB, data=avg.seizure) %>% summary
lm(avg.seizure.cocaine.2020~population+dist_from_SB, data=avg.seizure) %>% summary
lm(avg.seizure.cocaine.2021~population+dist_from_SB, data=avg.seizure) %>% summary


lm(avg.seizure.all.2018~population+long+lat+long*lat, data=avg.seizure) %>% summary
lm(avg.seizure.all.2019~population+long+lat+long*lat, data=avg.seizure) %>% summary
lm(avg.seizure.all.2020~population+long+lat+long*lat, data=avg.seizure) %>% summary
lm(avg.seizure.all.2021~population+long+lat+long*lat, data=avg.seizure) %>% summary
lm(avg.seizure.cocaine.2018~population+long+lat+long*lat, data=avg.seizure) %>% summary
lm(avg.seizure.cocaine.2019~population+long+lat+long*lat, data=avg.seizure) %>% summary
lm(avg.seizure.cocaine.2020~population+long+lat+long*lat, data=avg.seizure) %>% summary
lm(avg.seizure.cocaine.2021~population+long+lat+long*lat, data=avg.seizure) %>% summary


## 2. Number of zero-seizure counties per month
n_zeros <- LISA.all[,21:68] %>% apply(2, function(x) sum(x==0))

plot(ts(n_zeros, start=c(2018, 1), frequency=12),
     xlab="Year", ylab="# of zero seizure counties",
     main="1422 counties in total")
 # because of end of fiscal year (Sep) ?


## 3. Seizure patterns comparison
state_name <- "Kentucky"
plot(ts((LISA.all %>% filter(state==state_name))[, 21:68] %>% apply(2, sum), start=c(2018, 1), frequency=12),
     xlab="Year",
     ylab="Seizure Weights (kg)",
     main=paste("Monthly Seizure Weights in", state_name, "(All Drugs)"),
     pch=19, type="b")
plot(ts((LISA.cocaine %>% filter(state==state_name))[, 21:68] %>% apply(2, sum), start=c(2018, 1), frequency=12),
     xlab="Year",
     ylab="Seizure Weights (kg)",
     main=paste("Monthly Seizure Weights in", state_name, "(Cocaine)"),
     pch=19, type="b")

plot(ts((LISA.all %>% filter(state==state_name, NAMELSAD=="Jefferson County"))[, 21:68] %>% apply(2, sum), start=c(2018, 1), frequency=12),
     xlab="Year",
     ylab=paste("Seizure Weights in", state_name),
     pch=19, type="b")
plot(ts((LISA.all %>% filter(state==state_name, NAMELSAD=="Pike County"))[, 21:68] %>% apply(2, sum), start=c(2018, 1), frequency=12),
     xlab="Year",
     ylab=paste("Seizure Weights in", state_name),
     pch=19, type="b")


state_name <- "California"
plot(ts((LISA.all %>% filter(state==state_name))[, 21:68] %>% apply(2, sum), start=c(2018, 1), frequency=12),
     xlab="Year",
     ylab="Seizure Weights (kg)",
     main=paste("Monthly Seizure Weights in", state_name),
     pch=19, type="b")

state_name <- "Texas"
plot(ts((LISA.all %>% filter(state==state_name))[, 21:68] %>% apply(2, sum), start=c(2018, 1), frequency=12),
     xlab="Year",
     ylab="Seizure Weights (kg)",
     main=paste("Monthly Seizure Weights in", state_name),
     pch=19, type="b")
# patterns are very different by states

KT.seizures <- (LISA.all %>% filter(state=="Kentucky"))[,c(2,8,21:68)]
cbind(KT.seizures[,2], apply(KT.seizures[,-(1:2)], 1, sum))
par(mfrow=c(5,1))
for (i in 1:4) {
  county <- KT.seizures[i,]
  plot(ts(t(county[1, 3:50]), start=c(2018, 1), frequency=12),
       xlab="Year",
       ylab="Seizure Weights (kg)",
       main=county$NAMELSAD,
       pch=19, type="b")
}
state_name <- "Kentucky"
plot(ts((LISA.all %>% filter(state==state_name))[, 21:68] %>% apply(2, sum), start=c(2018, 1), frequency=12),
     xlab="Year",
     ylab=paste("Seizure Weights in", state_name),
     main="Kentucky",
     pch=19, type="b")
par(mfrow=c(1,1))

HIDTA.dist <- read.csv("HIDTA Regions.csv", header=T) %>% as_tibble
coordinate.HIDTA <- left_join(counties, HIDTA.dist, by=c("state_name", "county_name")) %>% filter(!(state_name %in% c("Alaska", "Hawaii", "Puerto Rico")))
names(coordinate.HIDTA)[7] <- "GEOID"


centroids.HIDTA <- coordinate.HIDTA %>% filter(!(state_name %in% c("Alaska", "Hawaii")), !is.na(HIDTA)) %>%
  group_by(HIDTA) %>% summarise(long=mean(long), lat=mean(lat))

HIDTA.regions <- coordinate.HIDTA %>% filter(!is.na(HIDTA)) %>% select(GEOID, HIDTA) %>% unique
LISA.all$HIDTA <- left_join(LISA.all, HIDTA.regions, by="GEOID")$HIDTA

seizures.HIDTA <- (LISA.all %>% group_by(HIDTA))[,c(2,8,213,21:68)] %>% na.omit
seizures.Appalachia <- seizures.HIDTA %>% filter(HIDTA=="Appalachia")

j <- 2
{
par(mfrow=c(5,1))
for (i in (j:(j+4))) {
  plot(ts(t(seizures.Appalachia[i, 4:51]), start=c(2018, 1), frequency=12),
       xlab="Year",
       ylab="Seizure Weights (kg)",
       main=paste(seizures.Appalachia[i, 1], seizures.Appalachia[i, 2]),
       pch=19, type="b")
}
par(mfrow=c(1,1))
}

seizures.non_App <- seizures.HIDTA %>% filter(HIDTA!="Appalachia")
n_plots <- 4
{
  county.index <- sample(1:441, n_plots)
  par(mfrow=c(n_plots,1))
  for (i in county.index) {
    plot(ts(t(seizures.non_App[i, 4:51]), start=c(2018, 1), frequency=12),
         xlab="Year",
         ylab="Seizure Weights (kg)",
         main=paste(seizures.non_App[i, 1], seizures.non_App[i, 2]),
         pch=19, type="b")
  }
  par(mfrow=c(1,1))
}
