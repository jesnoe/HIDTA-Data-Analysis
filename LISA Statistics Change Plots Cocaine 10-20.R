library(fpp2)
library(readxl)
library(maps)
library(urbnmapr)
library(tidyverse)
library(gridExtra)
library(rnaturalearth)

# 0: Insignificant
# 1: High/High
# 2: Low/Low
# 3: Low/High
# 4: High/Low

LISA.rel <- read.csv("CountyKNN3 R codes 999 (relative).csv") %>% as_tibble
LISA.abs.1 <- read.csv("CountyKNN3 R codes 999 (absolute base 1).csv") %>% as_tibble
LISA.abs.t_1 <- read.csv("CountyKNN3 R codes 999 (absolute base t_1).csv") %>% as_tibble

LISA.rel[,grep("LISA_I", names(LISA.rel))]; LISA.abs.1[,grep("LISA_I", names(LISA.rel))]; LISA.abs.t_1[,grep("LISA_I", names(LISA.rel))]

LISA_C.rel <- LISA.rel[,grep("LISA_C", names(LISA.rel))]
LISA_C.abs.1 <- LISA.abs.1[,grep("LISA_C", names(LISA.rel))]
LISA_C.abs.t_1 <- LISA.abs.t_1[,grep("LISA_C", names(LISA.rel))]
county.names <- LISA.rel[, c(2,6,8)]

# LISA_C.annual.rel <- county.names
# LISA_C.annual.abs.1 <- county.names
# LISA_C.annual.abs.t_1 <- county.names
# for (i in 1:4) {
#   LISA_C.annual.rel[[paste("Insig", 2017+i, sep="_")]] <- apply(LISA_C.rel[(12*i-11):(12*i)], 1, function(x) sum(x == 0))
#   LISA_C.annual.rel[[paste("LL", 2017+i, sep="_")]] <- apply(LISA_C.rel[(12*i-11):(12*i)], 1, function(x) sum(x == 2))
#   LISA_C.annual.rel[[paste("LH", 2017+i, sep="_")]] <- apply(LISA_C.rel[(12*i-11):(12*i)], 1, function(x) sum(x == 3))
#   LISA_C.annual.rel[[paste("HL", 2017+i, sep="_")]] <- apply(LISA_C.rel[(12*i-11):(12*i)], 1, function(x) sum(x == 4))
#   LISA_C.annual.rel[[paste("HH", 2017+i, sep="_")]] <- apply(LISA_C.rel[(12*i-11):(12*i)], 1, function(x) sum(x == 1))
#   
#   LISA_C.annual.abs.1[[paste("Insig", 2017+i, sep="_")]] <- apply(LISA_C.abs.1[(12*i-11):(12*i)], 1, function(x) sum(x == 0))
#   LISA_C.annual.abs.1[[paste("LL", 2017+i, sep="_")]] <- apply(LISA_C.abs.1[(12*i-11):(12*i)], 1, function(x) sum(x == 2))
#   LISA_C.annual.abs.1[[paste("LH", 2017+i, sep="_")]] <- apply(LISA_C.abs.1[(12*i-11):(12*i)], 1, function(x) sum(x == 3))
#   LISA_C.annual.abs.1[[paste("HL", 2017+i, sep="_")]] <- apply(LISA_C.abs.1[(12*i-11):(12*i)], 1, function(x) sum(x == 4))
#   LISA_C.annual.abs.1[[paste("HH", 2017+i, sep="_")]] <- apply(LISA_C.abs.1[(12*i-11):(12*i)], 1, function(x) sum(x == 1))
#   
#   LISA_C.annual.abs.t_1[[paste("Insig", 2017+i, sep="_")]] <- apply(LISA_C.abs.t_1[(12*i-11):(12*i)], 1, function(x) sum(x == 0))
#   LISA_C.annual.abs.t_1[[paste("LL", 2017+i, sep="_")]] <- apply(LISA_C.abs.t_1[(12*i-11):(12*i)], 1, function(x) sum(x == 2))
#   LISA_C.annual.abs.t_1[[paste("LH", 2017+i, sep="_")]] <- apply(LISA_C.abs.t_1[(12*i-11):(12*i)], 1, function(x) sum(x == 3))
#   LISA_C.annual.abs.t_1[[paste("HL", 2017+i, sep="_")]] <- apply(LISA_C.abs.t_1[(12*i-11):(12*i)], 1, function(x) sum(x == 4))
#   LISA_C.annual.abs.t_1[[paste("HH", 2017+i, sep="_")]] <- apply(LISA_C.abs.t_1[(12*i-11):(12*i)], 1, function(x) sum(x == 1))
#   
# }
# LISA_C.annual.rel %>% write.csv("KNN3_relative_annual_LISA_classes_999.csv", row.names=F)
# LISA_C.annual.abs.1 %>% write.csv("KNN3_absolute_b1_annual_LISA_classes_999.csv", row.names=F)
# LISA_C.annual.abs.t_1 %>% write.csv("KNN3_absolute_t1_annual_LISA_classes_999.csv", row.names=F)

LISA_C.rel <- cbind(county.names, LISA_C.rel)
LISA_C.abs.1 <- cbind(county.names, LISA_C.abs.1)
LISA_C.abs.t_1 <- cbind(county.names, LISA_C.abs.t_1)


## Monthly or Quarterly Heatmap

grep("Jan_2018", names(LISA.rel)) # 21
nrow(LISA.rel) # 792
seizure.acf <- apply(LISA.rel[, 21:68], 1, function(x) sum(acf(x, plot=F)$acf[-1] >= 0.1))
seizure.pacf <- apply(LISA.rel[, 21:68], 1, function(x) sum(pacf(x, plot=F)$acf[-1] >= 0.1))

par(mar=c(5,5,0,0), mfrow=c(1,2))
plot(as.factor(seizure.acf), xlab="# of lags with acf >= 0.1", ylab="# of counties")
plot(as.factor(seizure.pacf), xlab="# of lags with pacf >= 0.1", ylab="# of counties")


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

counties.original <- counties #<- counties.original
names(counties)[7] <- "GEOID"
counties$GEOID <- as.numeric(counties$GEOID)
LISA.rel <- left_join(LISA.rel, county.fips[,-(1:2)], by="GEOID")
# 216:219
LISA.rel.map <- left_join(counties[, c(1:2,7)], LISA.rel[, c(2, 6, 8, 21:68)], by = "GEOID")
LISA.rel.map <- LISA.rel.map %>% pivot_longer(c(-long, -lat, -GEOID, -state, -NAMELSAD), names_to="Month_Year", values_to="Seizure_weights")
LISA.rel.map <- left_join(LISA.rel.map, county.fips[,-(1:2)])
LISA.rel.map$Month_Year <- factor(LISA.rel.map$Month_Year, levels=unique(LISA.rel.map$Month_Year))
LISA.rel.map$weights_per_pop. <- 
  ifelse(substr(LISA.rel.map$Month_Year, 5, 8) == "2018", LISA.rel.map$Seizure_weights/LISA.rel.map$POPESTIMATE2018,
         ifelse(substr(LISA.rel.map$Month_Year, 5, 8) == "2019", LISA.rel.map$Seizure_weights/LISA.rel.map$POPESTIMATE2019,
                ifelse(substr(LISA.rel.map$Month_Year, 5, 8) == "2010", LISA.rel.map$Seizure_weights/LISA.rel.map$POPESTIMATE2020,
                       LISA.rel.map$Seizure_weights/LISA.rel.map$POPESTIMATE2021)))

for (year in 2018:2021) { # monthly maps of seziure weights
  LISA.rel.map %>% filter(grepl(year, Month_Year)) %>% 
    ggplot(mapping = aes(long, lat, group = GEOID, fill=log(Seizure_weights))) +
    geom_polygon(color = "#000000", size = .05) +
    facet_wrap(. ~ Month_Year) +
    coord_equal() +
    scale_fill_viridis_c() +
    labs(fill = "log Seizure Weights") -> seizure_weights_map
  ggsave(paste("seizure_weights_", year, ".jpg", sep=""), seizure_weights_map, width=55, height=30, units="cm")
  
  LISA.rel.map %>% filter(grepl("2018", Month_Year)) %>% 
    ggplot(mapping = aes(long, lat, group = GEOID, fill=log(weights_per_pop.))) +
    geom_polygon(color = "#000000", size = .05) +
    facet_wrap(. ~ Month_Year) +
    coord_equal() +
    scale_fill_viridis_c() +
    labs(fill = "log Weights Per Pop.") -> weights_per_pop.map
  ggsave(paste("log weights_per_pop_", year, ".jpg", sep=""), weights_per_pop.map, width=55, height=30, units="cm")
}

LISA_C.rel.map <- left_join(counties[, c(1:2,7)], LISA_C.rel, by = "GEOID")
LISA_C.rel.map <- LISA_C.rel.map %>% pivot_longer(c(-long, -lat, -GEOID, -state, -NAMELSAD), names_to="Month_Year", values_to="LISA_C")
ref <- tibble(LISA_C=0:4, label=c("Insig", "HH", "LL", "LH", "HL"))
LISA_C.rel.map$LISA_C <- left_join(LISA_C.rel.map[, 6:7], ref, by="LISA_C")$label
LISA_C.rel.map$LISA_C <- factor(LISA_C.rel.map$LISA_C)
LISA_C.rel.map$Month_Year <- factor(LISA_C.rel.map$Month_Year, levels=unique(LISA_C.rel.map$Month_Year))

for (year in 2018:2021) { # monthly maps of LISA_C.rel
  LISA_C.rel.map %>% filter(substr(Month_Year,10,10) == substr(year,4,4)) %>% 
    ggplot(mapping = aes(long, lat, group = GEOID, fill=LISA_C)) +
    geom_polygon(color = "#000000", size = .05) +
    facet_wrap(. ~ Month_Year) +
    scale_fill_manual(values = c("Insig"="grey90",
                                  "LL"="blue",
                                "LH"="steelblue",
                                "HL"="orange",
                                "HH"="red")) +
    labs(fill = "LISA_C_rel") -> LISA_C_map
  ggsave(paste("LISA_C_rel_", year, ".jpg", sep=""), LISA_C_map, width=55, height=30, units="cm")
}

LISA_C.abs1.map <- left_join(counties[, c(1:2,7)], LISA_C.abs.1, by = "GEOID")
LISA_C.abs1.map <- LISA_C.abs1.map %>% pivot_longer(c(-long, -lat, -GEOID, -state, -NAMELSAD), names_to="Month_Year", values_to="LISA_C")
LISA_C.abs1.map$LISA_C <- left_join(LISA_C.abs1.map[, 6:7], ref, by="LISA_C")$label
LISA_C.abs1.map$LISA_C <- factor(LISA_C.abs1.map$LISA_C)
LISA_C.abs1.map$Month_Year <- factor(LISA_C.abs1.map$Month_Year, levels=unique(LISA_C.abs1.map$Month_Year))

for (year in 2018:2021) { # monthly maps of LISA_C.abs1
  LISA_C.abs1.map %>% filter(substr(Month_Year,10,10) == substr(year,4,4)) %>% 
    ggplot(mapping = aes(long, lat, group = GEOID, fill=LISA_C)) +
    geom_polygon(color = "#000000", size = .05) +
    facet_wrap(. ~ Month_Year) +
    scale_fill_manual(values = c("Insig"="grey90",
                                 "LL"="blue",
                                 "LH"="steelblue",
                                 "HL"="orange",
                                 "HH"="red")) +
    labs(fill = "LISA_C_abs1") -> LISA_C_map
  ggsave(paste("LISA_C_abs1_", year, ".jpg", sep=""), LISA_C_map, width=55, height=30, units="cm")
}

LISA_C.abs_t1.map <- left_join(counties[, c(1:2,7)], LISA_C.abs.t_1, by = "GEOID")
LISA_C.abs_t1.map <- LISA_C.abs_t1.map %>% pivot_longer(c(-long, -lat, -GEOID, -state, -NAMELSAD), names_to="Month_Year", values_to="LISA_C")
LISA_C.abs_t1.map$LISA_C <- left_join(LISA_C.abs_t1.map[, 6:7], ref, by="LISA_C")$label
LISA_C.abs_t1.map$LISA_C <- factor(LISA_C.abs_t1.map$LISA_C)
LISA_C.abs_t1.map$Month_Year <- factor(LISA_C.abs_t1.map$Month_Year, levels=unique(LISA_C.abs_t1.map$Month_Year))

for (year in 2018:2021) { # monthly maps of LISA_C.abs_t1
  LISA_C.abs_t1.map %>% filter(substr(Month_Year,10,10) == substr(year,4,4)) %>% 
    ggplot(mapping = aes(long, lat, group = GEOID, fill=LISA_C)) +
    geom_polygon(color = "#000000", size = .05) +
    facet_wrap(. ~ Month_Year) +
    scale_fill_manual(values = c("Insig"="grey90",
                                 "LL"="blue",
                                 "LH"="steelblue",
                                 "HL"="orange",
                                 "HH"="red")) +
    labs(fill = "LISA_C_abs_t1") -> LISA_C_map
  ggsave(paste("LISA_C_abs_t1_", year, ".jpg", sep=""), LISA_C_map, width=55, height=30, units="cm")
}

LISA_C.abs_t1.map %>% 
  ggplot(mapping = aes(long, lat, group = GEOID, fill=!is.na(LISA_C))) +
  geom_polygon(color = "#000000", size = .05) +
  scale_fill_manual(values = c("FALSE"="grey90",
                               "TRUE"="steelblue")) +
  labs(fill = "Observed?") -> LISA_C_map
ggsave("Observed Counties.jpg", LISA_C_map, width=15, height=10, units="cm")
