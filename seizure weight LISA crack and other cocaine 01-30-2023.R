library(fpp2)
library(spdep)
library(readxl)
library(urbnmapr)
library(tidyverse)
library(gridExtra)

cocaine <- read.csv("cocaine HIDTA (12-07-2022).csv") %>% as_tibble
crack <- read.csv("crack cocaine HIDTA (12-07-2022).csv") %>% as_tibble
LISA3 <- read.csv("CountyKNN3.csv") %>% as_tibble %>% arrange(GEOID)

overdose <- read.csv("VSRR_Provisional_Drug_Overdose_Death_Counts.csv") %>% as_tibble
names(overdose)[1] <- "state_abbv"
cocaine_overdose <- overdose %>% filter(Indicator == "Cocaine (T40.5)" & Year %in% 2018:2021 & !(state_abbv %in% c("AK", "HI", "US")))
cocaine_overdose$Data.Value <- as.numeric(cocaine_overdose$Data.Value)
cocaine_overdose$Month <- str_sub(cocaine_overdose$Month, 1, 3)
cocaine_overdose$month <- str_c(cocaine_overdose$Month, cocaine_overdose$Year, sep="_")
cocaine_overdose$month <- factor(cocaine_overdose$month, levels=names(cocaine)[-(1:4)])
cocaine_overdose <- cocaine_overdose %>% rename(state=State.Name)
cocaine_overdose %>% complete(Year, Month)

counties.obs <- counties
names(counties.obs)[7] <- "GEOID"
coords.crack <- counties.obs %>% filter(GEOID %in% crack$GEOID) %>%
  group_by(GEOID) %>% summarise(x=mean(long), y=mean(lat))
GEOIDS.crack <- coords.crack$GEOID
coords.crack <- coords.crack[,-1]
nb_crack <- knn2nb(knearneigh(coords.crack, k=5), row.names=GEOIDS.crack)

alpha <- 0.05
nperm <- 999

## For crack
crack <- cbind(crack, matrix(0, nrow(crack), 48*3)) %>% as_tibble
names(crack)[(which(names(crack) == "1"):which(names(crack) == "144"))] <- names(LISA3)[71:214]
Jan_2018_index <- grep("Jan_2018", names(crack))[1]
LISA_I.index <- grep("LISA_I", names(crack))[1]
nb.obj.crack <- nb2listw(nb_crack, style="B")
nrow(crack) # 813
seizures.crack <- crack[, 5:52]

grep("LISA_I", names(crack))[1] # 53

relative.MoransI <- crack
set.seed(100)
for (i in 1:48) {
  seizure.crack <- t(seizures.crack)[i,]
  localM.month <- localmoran_abs(seizure.crack, nb.obj.crack, nsim=nperm, zero.policy=T, xx=NULL)
  localM.month$LISA_C <- ifelse(localM.month$quadr=="High-High", 1,
                                ifelse(localM.month$quadr=="Low-Low", 2,
                                       ifelse(localM.month$quadr=="Low-High", 3, 4)))
  localM.month$LISA_C <- ifelse(localM.month$`Pr(folded) Sim` <= alpha, localM.month$LISA_C, 0)
  relative.MoransI[,(LISA_I.index+3*i-3):(LISA_I.index+3*i-1)] <- localM.month[,c(1,13,7)]
}

# write.csv(relative.MoransI, "crack KNN5 R codes 999 (relative).csv", row.names=F)

absolute.b1.MoransI <- crack
x.bar_p <- mean(absolute.b1.MoransI$Jan_2018)
set.seed(300)
for (i in 1:48) {
  seizure.crack <- t(seizures.crack)[i,]
  localM.month <- localmoran_abs(seizure.crack, nb.obj.crack, nsim=nperm, zero.policy=T, xx=x.bar_p)
  localM.month$LISA_C <- ifelse(localM.month$quadr=="High-High", 1,
                                ifelse(localM.month$quadr=="Low-Low", 2,
                                       ifelse(localM.month$quadr=="Low-High", 3, 4)))
  localM.month$LISA_C <- ifelse(localM.month$`Pr(folded) Sim` <= alpha, localM.month$LISA_C, 0)
  absolute.b1.MoransI[,(LISA_I.index+3*i-3):(LISA_I.index+3*i-1)] <- localM.month[,c(1,13,7)]
}

# write.csv(absolute.b1.MoransI, "crack KNN5 R codes 999 (absolute base 1).csv", row.names=F)

absolute.t1.MoransI <- crack
x.bar_p <- mean(absolute.t1.MoransI$Jan_2018)
set.seed(500)
for (i in 1:48) {
  seizure.crack <- t(seizures.crack)[i,]
  localM.month <- localmoran_abs(seizure.crack, nb.obj.crack, nsim=nperm, zero.policy=T, xx=x.bar_p)
  localM.month$LISA_C <- ifelse(localM.month$quadr=="High-High", 1,
                                ifelse(localM.month$quadr=="Low-Low", 2,
                                       ifelse(localM.month$quadr=="Low-High", 3, 4)))
  localM.month$LISA_C <- ifelse(localM.month$`Pr(folded) Sim` <= alpha, localM.month$LISA_C, 0)
  absolute.t1.MoransI[,(LISA_I.index+3*i-3):(LISA_I.index+3*i-1)] <- localM.month[,c(1,13,7)]
  x.bar_p <- mean(as.data.frame(crack)[,Jan_2018_index+i-1])
}

# write.csv(absolute.t1.MoransI, "crack KNN5 R codes 999 (absolute base t_1).csv", row.names=F)

# crack Plots
LISA.rel <- read.csv("crack KNN5 R codes 999 (relative).csv") %>% as_tibble
LISA.abs.1 <- read.csv("crack KNN5 R codes 999 (absolute base 1).csv") %>% as_tibble
LISA.abs.t_1 <- read.csv("crack KNN5 R codes 999 (absolute base t_1).csv") %>% as_tibble

LISA_C.rel <- LISA.rel[,c(1:3, grep("LISA_C", names(LISA.rel)))]
LISA_C.abs.1 <- LISA.abs.1[,c(1:3, grep("LISA_C", names(LISA.rel)))]
LISA_C.abs.t_1 <- LISA.abs.t_1[,c(1:3, grep("LISA_C", names(LISA.rel)))]
county.names <- LISA.rel[, 1:3]

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

names(counties.obs)[7] <- "GEOID"
counties.obs <- counties.obs %>% rename(state=state_name, county=county_name)
counties.obs$GEOID <- as.numeric(counties.obs$GEOID)
counties.obs <- counties.obs %>% filter(!(state %in% c("Alaska", "Hawaii")))

LISA.rel.map <- left_join(counties.obs[, c(1:2,7,10,12)], LISA.rel[, c(3, 5:52)], by = "GEOID")
LISA.rel.map <- LISA.rel.map %>% pivot_longer(c(-long, -lat, -GEOID, -state, -county), names_to="Month_Year", values_to="Seizure_weights")
LISA.rel.map <- left_join(LISA.rel.map, county.fips[,-(1:2)])
LISA.rel.map$Month_Year <- factor(LISA.rel.map$Month_Year, levels=unique(LISA.rel.map$Month_Year))
LISA.rel.map$weights_per_pop. <- 
  ifelse(substr(LISA.rel.map$Month_Year, 5, 8) == "2018", LISA.rel.map$Seizure_weights/LISA.rel.map$POPESTIMATE2018,
         ifelse(substr(LISA.rel.map$Month_Year, 5, 8) == "2019", LISA.rel.map$Seizure_weights/LISA.rel.map$POPESTIMATE2019,
                ifelse(substr(LISA.rel.map$Month_Year, 5, 8) == "2010", LISA.rel.map$Seizure_weights/LISA.rel.map$POPESTIMATE2020,
                       LISA.rel.map$Seizure_weights/LISA.rel.map$POPESTIMATE2021)))

overdose_map <- cocaine_overdose %>% select(state, month, Data.Value) %>% rename(Month_Year=month)
LISA.rel.map <- left_join(LISA.rel.map, overdose_map, by=c("state", "Month_Year"))
LISA.rel.map <- LISA.rel.map %>% rename(overdose=Data.Value)
LISA.rel.map$overdose_per_pop. <- 
  ifelse(substr(LISA.rel.map$Month_Year, 5, 8) == "2018", LISA.rel.map$overdose/LISA.rel.map$POPESTIMATE2018,
         ifelse(substr(LISA.rel.map$Month_Year, 5, 8) == "2019", LISA.rel.map$overdose/LISA.rel.map$POPESTIMATE2019,
                ifelse(substr(LISA.rel.map$Month_Year, 5, 8) == "2010", LISA.rel.map$overdose/LISA.rel.map$POPESTIMATE2020,
                       LISA.rel.map$overdose/LISA.rel.map$POPESTIMATE2021)))

for (year in 2018:2021) { # monthly maps of overdose count
  LISA.rel.map %>% filter(grepl(year, Month_Year)) %>% 
    ggplot(mapping = aes(long, lat, group = GEOID, fill=overdose)) +
    geom_polygon(color = "#FFFFFF", size = .05) +
    facet_wrap(. ~ Month_Year) +
    scale_fill_viridis_c() +
    labs(fill = "Overdose Count") -> seizure_weights_map
  ggsave(paste("cocaine overdose_", year, ".jpg", sep=""), seizure_weights_map, width=50, height=30, units="cm")
  
  LISA.rel.map %>% filter(grepl(year, Month_Year)) %>% 
    ggplot(mapping = aes(long, lat, group = GEOID, fill=overdose_per_pop.)) +
    geom_polygon(color = "#FFFFFF", size = .05) +
    facet_wrap(. ~ Month_Year) +
    labs(fill = "Overdose Count Per Pop.") -> weights_per_pop.map
  ggsave(paste("cocaine overdose_per_pop_", year, ".jpg", sep=""), weights_per_pop.map, width=50, height=30, units="cm")
}


for (year in 2018:2021) { # monthly maps of seziure weights
  LISA.rel.map %>% filter(grepl(year, Month_Year)) %>% 
    ggplot(mapping = aes(long, lat, group = GEOID, fill=log(Seizure_weights))) +
    geom_polygon(color = "#FFFFFF", size = .05) +
    facet_wrap(. ~ Month_Year) +
    # coord_equal() +
    scale_fill_viridis_c() +
    labs(fill = "log Seizure Weights") -> seizure_weights_map
  ggsave(paste("crack seizure_weights_", year, ".jpg", sep=""), seizure_weights_map, width=50, height=30, units="cm")
  
  LISA.rel.map %>% filter(grepl(year, Month_Year)) %>% 
    ggplot(mapping = aes(long, lat, group = GEOID, fill=log(weights_per_pop.))) +
    geom_polygon(color = "#FFFFFF", size = .05) +
    facet_wrap(. ~ Month_Year) +
    # coord_equal() +
    scale_fill_viridis_c() +
    labs(fill = "log Weights Per Pop.") -> weights_per_pop.map
  ggsave(paste("crack log weights_per_pop_", year, ".jpg", sep=""), weights_per_pop.map, width=50, height=30, units="cm")
}


LISA_C.rel.map <- left_join(counties.obs[, c(1:2,7,10,12)], LISA_C.rel[,-(1:2)], by = "GEOID")
LISA_C.rel.map <- LISA_C.rel.map %>% pivot_longer(c(-long, -lat, -GEOID, -state, -county), names_to="Month_Year", values_to="LISA_C")
ref <- tibble(LISA_C=0:4, label=c("Insig", "HH", "LL", "LH", "HL"))
LISA_C.rel.map$LISA_C <- left_join(LISA_C.rel.map[, 6:7], ref, by="LISA_C")$label
LISA_C.rel.map$LISA_C <- factor(LISA_C.rel.map$LISA_C)
LISA_C.rel.map$Month_Year <- factor(LISA_C.rel.map$Month_Year, levels=unique(LISA_C.rel.map$Month_Year))

for (year in 2018:2021) { # monthly maps of LISA_C.rel
  LISA_C.rel.map %>% filter(substr(Month_Year,10,10) == substr(year,4,4)) %>% 
    ggplot(mapping = aes(long, lat, group = GEOID, fill=LISA_C)) +
    geom_polygon(color = "#FFFFFF", size = .05) +
    facet_wrap(. ~ Month_Year) +
    scale_fill_manual(values = c("Insig"="grey90",
                                 "LL"="blue",
                                 "LH"="steelblue",
                                 "HL"="orange",
                                 "HH"="red")) +
    labs(fill = "LISA_C_rel") -> LISA_C_map
  ggsave(paste("Crack_LISA_C_rel_", year, ".jpg", sep=""), LISA_C_map, width=50, height=30, units="cm")
}

LISA_C.abs1.map <- left_join(counties.obs[, c(1:2,7,10,12)], LISA_C.abs.1[,-(1:2)], by = "GEOID")
LISA_C.abs1.map <- LISA_C.abs1.map %>% pivot_longer(c(-long, -lat, -GEOID, -state, -county), names_to="Month_Year", values_to="LISA_C")
LISA_C.abs1.map$LISA_C <- left_join(LISA_C.abs1.map[, 6:7], ref, by="LISA_C")$label
LISA_C.abs1.map$LISA_C <- factor(LISA_C.abs1.map$LISA_C)
LISA_C.abs1.map$Month_Year <- factor(LISA_C.abs1.map$Month_Year, levels=unique(LISA_C.abs1.map$Month_Year))

for (year in 2018:2021) { # monthly maps of LISA_C.abs1
  LISA_C.abs1.map %>% filter(substr(Month_Year,10,10) == substr(year,4,4)) %>% 
    ggplot(mapping = aes(long, lat, group = GEOID, fill=LISA_C)) +
    geom_polygon(color = "#FFFFFF", size = .05) +
    facet_wrap(. ~ Month_Year) +
    scale_fill_manual(values = c("Insig"="grey90",
                                 "LL"="blue",
                                 "LH"="steelblue",
                                 "HL"="orange",
                                 "HH"="red")) +
    labs(fill = "LISA_C_abs1") -> LISA_C_map
  ggsave(paste("Crack_LISA_C_abs1_", year, ".jpg", sep=""), LISA_C_map, width=50, height=30, units="cm")
}

LISA_C.abs_t1.map <- left_join(counties.obs[, c(1:2,7,10,12)], LISA_C.abs.t_1[,-(1:2)], by = "GEOID")
LISA_C.abs_t1.map <- LISA_C.abs_t1.map %>% pivot_longer(c(-long, -lat, -GEOID, -state, -county), names_to="Month_Year", values_to="LISA_C")
LISA_C.abs_t1.map$LISA_C <- left_join(LISA_C.abs_t1.map[, 6:7], ref, by="LISA_C")$label
LISA_C.abs_t1.map$LISA_C <- factor(LISA_C.abs_t1.map$LISA_C)
LISA_C.abs_t1.map$Month_Year <- factor(LISA_C.abs_t1.map$Month_Year, levels=unique(LISA_C.abs_t1.map$Month_Year))

for (year in 2018:2021) { # monthly maps of LISA_C.abs_t1
  LISA_C.abs_t1.map %>% filter(substr(Month_Year,10,10) == substr(year,4,4)) %>% 
    ggplot(mapping = aes(long, lat, group = GEOID, fill=LISA_C)) +
    geom_polygon(color = "#FFFFFF", size = .05) +
    facet_wrap(. ~ Month_Year) +
    scale_fill_manual(values = c("Insig"="grey90",
                                 "LL"="blue",
                                 "LH"="steelblue",
                                 "HL"="orange",
                                 "HH"="red")) +
    labs(fill = "LISA_C_abs_t1") -> LISA_C_map
  ggsave(paste("Crack_LISA_C_abs_t1_", year, ".jpg", sep=""), LISA_C_map, width=50, height=30, units="cm")
}



## For all cocaine
coords.cocaine <- counties.obs %>% filter(GEOID %in% cocaine$GEOID) %>%
  group_by(GEOID) %>% summarise(x=mean(long), y=mean(lat))
GEOIDS.cocaine <- coords.cocaine$GEOID
coords.cocaine <- coords.cocaine[,-1]
nb_cocaine <- knn2nb(knearneigh(coords.cocaine, k=5), row.names=GEOIDS.cocaine)

cocaine <- cbind(cocaine, matrix(0, nrow(cocaine), 48*3)) %>% as_tibble
names(cocaine)[(which(names(cocaine) == "1"):which(names(cocaine) == "144"))] <- names(LISA3)[71:214]
Jan_2018_index <- grep("Jan_2018", names(cocaine))[1]
LISA_I.index <- grep("LISA_I", names(cocaine))[1]
nb.obj.cocaine <- nb2listw(nb_cocaine, style="B")
nrow(cocaine) # 1092
seizures.cocaine <- cocaine[, 5:52]

grep("LISA_I", names(cocaine))[1] # 53

relative.MoransI <- cocaine
set.seed(100)
for (i in 1:48) {
  seizure.cocaine <- t(seizures.cocaine)[i,]
  localM.month <- localmoran_abs(seizure.cocaine, nb.obj.cocaine, nsim=nperm, zero.policy=T, xx=NULL)
  localM.month$LISA_C <- ifelse(localM.month$quadr=="High-High", 1,
                                ifelse(localM.month$quadr=="Low-Low", 2,
                                       ifelse(localM.month$quadr=="Low-High", 3, 4)))
  localM.month$LISA_C <- ifelse(localM.month$`Pr(folded) Sim` <= alpha, localM.month$LISA_C, 0)
  relative.MoransI[,(LISA_I.index+3*i-3):(LISA_I.index+3*i-1)] <- localM.month[,c(1,13,7)]
}

# write.csv(relative.MoransI, "cocaine KNN5 R codes 999 (relative).csv", row.names=F)

absolute.b1.MoransI <- cocaine
x.bar_p <- mean(absolute.b1.MoransI$Jan_2018)
set.seed(300)
for (i in 1:48) {
  seizure.cocaine <- t(seizures.cocaine)[i,]
  localM.month <- localmoran_abs(seizure.cocaine, nb.obj.cocaine, nsim=nperm, zero.policy=T, xx=x.bar_p)
  localM.month$LISA_C <- ifelse(localM.month$quadr=="High-High", 1,
                                ifelse(localM.month$quadr=="Low-Low", 2,
                                       ifelse(localM.month$quadr=="Low-High", 3, 4)))
  localM.month$LISA_C <- ifelse(localM.month$`Pr(folded) Sim` <= alpha, localM.month$LISA_C, 0)
  absolute.b1.MoransI[,(LISA_I.index+3*i-3):(LISA_I.index+3*i-1)] <- localM.month[,c(1,13,7)]
}

# write.csv(absolute.b1.MoransI, "cocaine KNN5 R codes 999 (absolute base 1).csv", row.names=F)

absolute.t1.MoransI <- cocaine
x.bar_p <- mean(absolute.t1.MoransI$Jan_2018)
set.seed(500)
for (i in 1:48) {
  seizure.cocaine <- t(seizures.cocaine)[i,]
  localM.month <- localmoran_abs(seizure.cocaine, nb.obj.cocaine, nsim=nperm, zero.policy=T, xx=x.bar_p)
  localM.month$LISA_C <- ifelse(localM.month$quadr=="High-High", 1,
                                ifelse(localM.month$quadr=="Low-Low", 2,
                                       ifelse(localM.month$quadr=="Low-High", 3, 4)))
  localM.month$LISA_C <- ifelse(localM.month$`Pr(folded) Sim` <= alpha, localM.month$LISA_C, 0)
  absolute.t1.MoransI[,(LISA_I.index+3*i-3):(LISA_I.index+3*i-1)] <- localM.month[,c(1,13,7)]
  x.bar_p <- mean(as.data.frame(cocaine)[,Jan_2018_index+i-1])
}

# write.csv(absolute.t1.MoransI, "cocaine KNN5 R codes 999 (absolute base t_1).csv", row.names=F)

# Cocaine Plots
LISA.rel <- read.csv("cocaine KNN5 R codes 999 (relative).csv") %>% as_tibble
LISA.abs.1 <- read.csv("cocaine KNN5 R codes 999 (absolute base 1).csv") %>% as_tibble
LISA.abs.t_1 <- read.csv("cocaine KNN5 R codes 999 (absolute base t_1).csv") %>% as_tibble

LISA_C.rel <- LISA.rel[,c(1:3, grep("LISA_C", names(LISA.rel)))]
LISA_C.abs.1 <- LISA.abs.1[,c(1:3, grep("LISA_C", names(LISA.rel)))]
LISA_C.abs.t_1 <- LISA.abs.t_1[,c(1:3, grep("LISA_C", names(LISA.rel)))]
county.names <- LISA.rel[, 1:3]


LISA.rel.map <- left_join(counties.obs[, c(1:2,7,10,12)], LISA.rel[, c(3, 5:52)], by = "GEOID")
LISA.rel.map <- LISA.rel.map %>% pivot_longer(c(-long, -lat, -GEOID, -state, -county), names_to="Month_Year", values_to="Seizure_weights")
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
    geom_polygon(color = "#FFFFFF", size = .05) +
    facet_wrap(. ~ Month_Year) +
    # coord_equal() +
    scale_fill_viridis_c() +
    labs(fill = "log Seizure Weights") -> seizure_weights_map
  ggsave(paste("cocaine seizure_weights_", year, ".jpg", sep=""), seizure_weights_map, width=50, height=30, units="cm")
  
  LISA.rel.map %>% filter(grepl(year, Month_Year)) %>% 
    ggplot(mapping = aes(long, lat, group = GEOID, fill=log(weights_per_pop.))) +
    geom_polygon(color = "#FFFFFF", size = .05) +
    facet_wrap(. ~ Month_Year) +
    # coord_equal() +
    scale_fill_viridis_c() +
    labs(fill = "log Weights Per Pop.") -> weights_per_pop.map
  ggsave(paste("cocaine log weights_per_pop_", year, ".jpg", sep=""), weights_per_pop.map, width=50, height=30, units="cm")
}


LISA_C.rel.map <- left_join(counties.obs[, c(1:2,7,10,12)], LISA_C.rel[,-(1:2)], by = "GEOID")
LISA_C.rel.map <- LISA_C.rel.map %>% pivot_longer(c(-long, -lat, -GEOID, -state, -county), names_to="Month_Year", values_to="LISA_C")
ref <- tibble(LISA_C=0:4, label=c("Insig", "HH", "LL", "LH", "HL"))
LISA_C.rel.map$LISA_C <- left_join(LISA_C.rel.map[, 6:7], ref, by="LISA_C")$label
LISA_C.rel.map$LISA_C <- factor(LISA_C.rel.map$LISA_C)
LISA_C.rel.map$Month_Year <- factor(LISA_C.rel.map$Month_Year, levels=unique(LISA_C.rel.map$Month_Year))

for (year in 2018:2021) { # monthly maps of LISA_C.rel
  LISA_C.rel.map %>% filter(substr(Month_Year,10,10) == substr(year,4,4)) %>% 
    ggplot(mapping = aes(long, lat, group = GEOID, fill=LISA_C)) +
    geom_polygon(color = "#FFFFFF", size = .05) +
    facet_wrap(. ~ Month_Year) +
    scale_fill_manual(values = c("Insig"="grey90",
                                 "LL"="blue",
                                 "LH"="steelblue",
                                 "HL"="orange",
                                 "HH"="red")) +
    labs(fill = "LISA_C_rel") -> LISA_C_map
  ggsave(paste("cocaine_LISA_C_rel_", year, ".jpg", sep=""), LISA_C_map, width=50, height=30, units="cm")
}

LISA_C.abs1.map <- left_join(counties.obs[, c(1:2,7,10,12)], LISA_C.abs.1[,-(1:2)], by = "GEOID")
LISA_C.abs1.map <- LISA_C.abs1.map %>% pivot_longer(c(-long, -lat, -GEOID, -state, -county), names_to="Month_Year", values_to="LISA_C")
LISA_C.abs1.map$LISA_C <- left_join(LISA_C.abs1.map[, 6:7], ref, by="LISA_C")$label
LISA_C.abs1.map$LISA_C <- factor(LISA_C.abs1.map$LISA_C)
LISA_C.abs1.map$Month_Year <- factor(LISA_C.abs1.map$Month_Year, levels=unique(LISA_C.abs1.map$Month_Year))

for (year in 2018:2021) { # monthly maps of LISA_C.abs1
  LISA_C.abs1.map %>% filter(substr(Month_Year,10,10) == substr(year,4,4)) %>% 
    ggplot(mapping = aes(long, lat, group = GEOID, fill=LISA_C)) +
    geom_polygon(color = "#FFFFFF", size = .05) +
    facet_wrap(. ~ Month_Year) +
    scale_fill_manual(values = c("Insig"="grey90",
                                 "LL"="blue",
                                 "LH"="steelblue",
                                 "HL"="orange",
                                 "HH"="red")) +
    labs(fill = "LISA_C_abs1") -> LISA_C_map
  ggsave(paste("cocaine_LISA_C_abs1_", year, ".jpg", sep=""), LISA_C_map, width=50, height=30, units="cm")
}

LISA_C.abs_t1.map <- left_join(counties.obs[, c(1:2,7,10,12)], LISA_C.abs.t_1[,-(1:2)], by = "GEOID")
LISA_C.abs_t1.map <- LISA_C.abs_t1.map %>% pivot_longer(c(-long, -lat, -GEOID, -state, -county), names_to="Month_Year", values_to="LISA_C")
LISA_C.abs_t1.map$LISA_C <- left_join(LISA_C.abs_t1.map[, 6:7], ref, by="LISA_C")$label
LISA_C.abs_t1.map$LISA_C <- factor(LISA_C.abs_t1.map$LISA_C)
LISA_C.abs_t1.map$Month_Year <- factor(LISA_C.abs_t1.map$Month_Year, levels=unique(LISA_C.abs_t1.map$Month_Year))

for (year in 2018:2021) { # monthly maps of LISA_C.abs_t1
  LISA_C.abs_t1.map %>% filter(substr(Month_Year,10,10) == substr(year,4,4)) %>% 
    ggplot(mapping = aes(long, lat, group = GEOID, fill=LISA_C)) +
    geom_polygon(color = "#FFFFFF", size = .05) +
    facet_wrap(. ~ Month_Year) +
    scale_fill_manual(values = c("Insig"="grey90",
                                 "LL"="blue",
                                 "LH"="steelblue",
                                 "HL"="orange",
                                 "HH"="red")) +
    labs(fill = "LISA_C_abs_t1") -> LISA_C_map
  ggsave(paste("cocaine_LISA_C_abs_t1_", year, ".jpg", sep=""), LISA_C_map, width=50, height=30, units="cm")
}