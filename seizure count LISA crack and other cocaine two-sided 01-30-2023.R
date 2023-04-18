library(fpp2)
library(spdep)
library(readxl)
library(urbnmapr)
library(tidyverse)
library(gridExtra)
library(lubridate)

cocaine <- read.csv("cocaine other count HIDTA (01-30-2023).csv") %>% as_tibble
crack <- read.csv("cocaine crack count HIDTA (12-20-2022).csv") %>% as_tibble
LISA3 <- read.csv("CountyKNN3.csv") %>% as_tibble %>% arrange(GEOID)
overdose <- read.csv("VSRR_Provisional_County-Level_Drug_Overdose_Death_Counts.csv") %>% as_tibble
names(overdose)[1] <- "state_abbv"
cocaine_overdose <- overdose %>% filter(Year %in% 2018:2021 & !(state_abbv %in% c("AK", "HI", "US")))
cocaine_overdose$Provisional.Drug.Overdose.Deaths <- as.numeric(cocaine_overdose$Provisional.Drug.Overdose.Deaths)
cocaine_overdose$Month <- as.numeric(cocaine_overdose$Month)
cocaine_overdose <- cocaine_overdose %>% complete(Year, Month)
cocaine_overdose <- cocaine_overdose %>% 
  mutate(month_year=str_c(Month, Year, sep="-") %>% 
           parse_date("%m-%Y") %>% format("%b_%Y"))
# cocaine_overdose$month_year <- factor(cocaine_overdose$month_year, levels=names(cocaine)[-(1:4)])
cocaine_overdose <- cocaine_overdose %>% rename(state=STATE_NAME, county=COUNTYNAME, GEOID=FIPS)


counties.obs <- counties
names(counties.obs)[7] <- "GEOID"
counties.obs <- counties.obs %>% rename(state=state_name, county=county_name)
counties.obs$GEOID <- as.numeric(counties.obs$GEOID)
counties.obs <- counties.obs %>% filter(!(state %in% c("Alaska", "Hawaii")))
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
  localM.month <- localmoran_abs(seizure.crack, nb.obj.crack, nsim=nperm, zero.policy=T, xx=NULL, alternative="two.sided")
  localM.month$LISA_C <- ifelse(localM.month$quadr=="High-High", 1,
                                ifelse(localM.month$quadr=="Low-Low", 2,
                                       ifelse(localM.month$quadr=="Low-High", 3, 4)))
  localM.month$LISA_C <- ifelse(localM.month$`Pr(folded) Sim` <= alpha, localM.month$LISA_C, 0)
  relative.MoransI[,(LISA_I.index+3*i-3):(LISA_I.index+3*i-1)] <- localM.month[,c(1,13,7)]
}

# write.csv(relative.MoransI, "crack counts KNN5 R codes 999 two-sided (relative).csv", row.names=F)

absolute.b1.MoransI <- crack
x.bar_p <- mean(absolute.b1.MoransI$Jan_2018)
set.seed(300)
for (i in 1:48) {
  seizure.crack <- t(seizures.crack)[i,]
  localM.month <- localmoran_abs(seizure.crack, nb.obj.crack, nsim=nperm, zero.policy=T, xx=x.bar_p, alternative="two.sided")
  localM.month$LISA_C <- ifelse(localM.month$quadr=="High-High", 1,
                                ifelse(localM.month$quadr=="Low-Low", 2,
                                       ifelse(localM.month$quadr=="Low-High", 3, 4)))
  localM.month$LISA_C <- ifelse(localM.month$`Pr(folded) Sim` <= alpha, localM.month$LISA_C, 0)
  absolute.b1.MoransI[,(LISA_I.index+3*i-3):(LISA_I.index+3*i-1)] <- localM.month[,c(1,13,7)]
}

# write.csv(absolute.b1.MoransI, "crack counts KNN5 R codes 999 two-sided (absolute base 1).csv", row.names=F)

absolute.t1.MoransI <- crack
x.bar_p <- mean(absolute.t1.MoransI$Jan_2018)
set.seed(500)
for (i in 1:48) {
  seizure.crack <- t(seizures.crack)[i,]
  localM.month <- localmoran_abs(seizure.crack, nb.obj.crack, nsim=nperm, zero.policy=T, xx=x.bar_p, alternative="two.sided")
  localM.month$LISA_C <- ifelse(localM.month$quadr=="High-High", 1,
                                ifelse(localM.month$quadr=="Low-Low", 2,
                                       ifelse(localM.month$quadr=="Low-High", 3, 4)))
  localM.month$LISA_C <- ifelse(localM.month$`Pr(folded) Sim` <= alpha, localM.month$LISA_C, 0)
  absolute.t1.MoransI[,(LISA_I.index+3*i-3):(LISA_I.index+3*i-1)] <- localM.month[,c(1,13,7)]
  x.bar_p <- mean(as.data.frame(crack)[,Jan_2018_index+i-1])
}

# write.csv(absolute.t1.MoransI, "crack counts KNN5 R codes 999 two-sided (absolute base t_1).csv", row.names=F)

# crack Plots
LISA.rel <- read.csv("crack counts KNN5 R codes 999 two-sided (relative).csv") %>% as_tibble
LISA.abs.1 <- read.csv("crack counts KNN5 R codes 999 two-sided (absolute base 1).csv") %>% as_tibble
LISA.abs.t_1 <- read.csv("crack counts KNN5 R codes 999 two-sided (absolute base t_1).csv") %>% as_tibble

LISA_C.rel <- LISA.rel[,c(1:3, grep("LISA_C", names(LISA.rel)))]
LISA_C.abs.1 <- LISA.abs.1[,c(1:3, grep("LISA_C", names(LISA.rel)))]
LISA_C.abs.t_1 <- LISA.abs.t_1[,c(1:3, grep("LISA_C", names(LISA.rel)))]
county.names <- LISA.rel[, 1:3]

LISA_C_to_Month <- function(col_names) {
  col_names <- gsub("0", "2020", col_names)
  col_names <- gsub("1", "2021", col_names)
  col_names <- gsub("8", "2018", col_names)
  col_names <- gsub("9", "2019", col_names)
  col_names <- parse_date(col_names, "LISA_C%b%Y") %>% substr(1,7)
  return(col_names)
}

names(LISA_C.rel)[4:51] <- LISA_C_to_Month(names(LISA_C.rel)[4:51])
names(LISA_C.abs.1)[4:51] <- LISA_C_to_Month(names(LISA_C.abs.1)[4:51])
names(LISA_C.abs.t_1)[4:51] <- LISA_C_to_Month(names(LISA_C.abs.t_1)[4:51])


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


# Number of HL maps
n_HL_map <- function(row) {
  result <- ifelse(sum(is.na(row) == length(row)),
                   NA, sum(row==4, na.rm=T))
  return(result)
}
LISA_n_HL <- LISA_C.rel
for (year in 2018:2021) {
  LISA_C_index <- grep(year, names(LISA_n_HL))
  LISA_n_HL[[paste("n_HL", year, sep="_")]] <- apply(LISA_n_HL[,LISA_C_index], 1, n_HL_map)
}
LISA_n_HL <- LISA_n_HL %>% select(state:GEOID, n_HL_2018:n_HL_2021)
LISA_n_HL$n_HL_All <- apply(LISA_n_HL[, 4:7], 1, sum)
LISA_n_HL.map <- left_join(counties.obs[, c(1:2,6:7,10,12)], LISA_n_HL[,c(3:8)], by = "GEOID")
LISA_n_HL.map %>% pivot_longer(-c(long, lat, group, GEOID, state, county, n_HL_All), names_to="year", values_to="n_HL_year") %>% 
  ggplot(mapping = aes(long, lat, group = group, fill=n_HL_year)) +
  geom_polygon(color = "#000000", size = .05) +
  facet_wrap(. ~ year) +
  scale_fill_viridis_c(na.value="white") +
  labs(fill = "# of HL") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) -> n_HL_map
# ggsave(paste("crack n_HL annual.jpg", sep=""), n_HL_map, width=20, height=15, units="cm")

LISA_n_HL.map %>% 
  ggplot(mapping = aes(long, lat, group = group, fill=n_HL_All)) +
  geom_polygon(color = "#000000", size = .05) +
  scale_fill_viridis_c(na.value="white") +
  labs(fill = "# of HL", title="# of HL 2018-2021 (Crack)") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) -> n_HL_map_all
# ggsave(paste("crack n_HL 2018-2021.jpg", sep=""), n_HL_map_all, width=20, height=15, units="cm")


LISA_C.rel.map <- left_join(counties.obs[, c(1:2,7,10,12)], LISA_C.rel[,-(1:2)], by = "GEOID")
LISA_C.rel.map <- LISA_C.rel.map %>% pivot_longer(c(-long, -lat, -GEOID, -state, -county), names_to="Month_Year", values_to="LISA_C")
ref <- tibble(LISA_C=0:4, label=c("Insig", "HH", "LL", "LH", "HL"))
LISA_C.rel.map$LISA_C <- left_join(LISA_C.rel.map[, 6:7], ref, by="LISA_C")$label
LISA_C.rel.map$LISA_C <- factor(LISA_C.rel.map$LISA_C)
LISA_C.rel.map$Month_Year <- parse_date(LISA_C.rel.map$Month_Year, "%Y-%m")

for (year in 2018:2021) { # monthly maps of LISA_C.rel
  LISA_C.rel.map %>% filter(year(Month_Year) == year & month(Month_Year) %in% 1:4) %>% 
    ggplot(mapping = aes(long, lat, group = group, fill=LISA_C)) +
    geom_polygon(color = "#000000", size = .05) +
    facet_wrap(. ~ substr(Month_Year, 1, 7)) +
    scale_fill_manual(values = c("Insig"="grey60",
                                 "LL"="blue",
                                 "LH"="steelblue",
                                 "HL"="orange",
                                 "HH"="red"),
                      na.value = "white") +
    labs(fill = "LISA_C_rel") + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) -> LISA_C_map
  ggsave(paste("Crack_Count_LISA_C_rel_", year, " two-sided (4 months).jpg", sep=""), LISA_C_map, width=20, height=15, units="cm")
}

LISA_C.abs1.map <- left_join(counties.obs[, c(1:2,7,10,12)], LISA_C.abs.1[,-(1:2)], by = "GEOID")
LISA_C.abs1.map <- LISA_C.abs1.map %>% pivot_longer(c(-long, -lat, -GEOID, -state, -county), names_to="Month_Year", values_to="LISA_C")
LISA_C.abs1.map$LISA_C <- left_join(LISA_C.abs1.map[, 6:7], ref, by="LISA_C")$label
LISA_C.abs1.map$LISA_C <- factor(LISA_C.abs1.map$LISA_C)
LISA_C.abs1.map$Month_Year <- parse_date(LISA_C.abs1.map$Month_Year, "%Y-%m")

for (year in 2018:2021) { # monthly maps of LISA_C.abs1
  LISA_C.abs1.map %>% filter(year(Month_Year) == year & month(Month_Year) %in% 1:4) %>% 
    ggplot(mapping = aes(long, lat, group = group, fill=LISA_C)) +
    geom_polygon(color = "#000000", size = .05) +
    facet_wrap(. ~ substr(Month_Year, 1, 7)) +
    scale_fill_manual(values = c("Insig"="grey60",
                                 "LL"="blue",
                                 "LH"="steelblue",
                                 "HL"="orange",
                                 "HH"="red"),
                      na.value = "white") +
    labs(fill = "LISA_C_abs1") + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) -> LISA_C_map
  ggsave(paste("Crack_Count_LISA_C_abs1_", year, " two-sided (4 months).jpg", sep=""), LISA_C_map, width=20, height=15, units="cm")
}

LISA_C.abs_t1.map <- left_join(counties.obs[, c(1:2,7,10,12)], LISA_C.abs.t_1[,-(1:2)], by = "GEOID")
LISA_C.abs_t1.map <- LISA_C.abs_t1.map %>% pivot_longer(c(-long, -lat, -GEOID, -state, -county), names_to="Month_Year", values_to="LISA_C")
LISA_C.abs_t1.map$LISA_C <- left_join(LISA_C.abs_t1.map[, 6:7], ref, by="LISA_C")$label
LISA_C.abs_t1.map$LISA_C <- factor(LISA_C.abs_t1.map$LISA_C)
LISA_C.abs_t1.map$Month_Year <- parse_date(LISA_C.abs_t1.map$Month_Year, "%Y-%m")

for (year in 2018:2021) { # monthly maps of LISA_C.abs_t1
  LISA_C.abs_t1.map %>% filter(year(Month_Year) == year & month(Month_Year) %in% 1:4) %>% 
    ggplot(mapping = aes(long, lat, group = group, fill=LISA_C)) +
    geom_polygon(color = "#000000", size = .05) +
    facet_wrap(. ~ substr(Month_Year, 1, 7)) +
    scale_fill_manual(values = c("Insig"="grey60",
                                 "LL"="blue",
                                 "LH"="steelblue",
                                 "HL"="orange",
                                 "HH"="red"),
                      na.value = "white") +
    labs(fill = "LISA_C_abs_t1") + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) -> LISA_C_map
  ggsave(paste("Crack_Count_LISA_C_abs_t1_", year, " two-sided (4 months).jpg", sep=""), LISA_C_map, width=20, height=15, units="cm")
}



### For other cocaine
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
nrow(cocaine) # 1070
seizures.cocaine <- cocaine[, 5:52]

relative.MoransI <- cocaine
set.seed(100)
for (i in 1:48) {
  seizure.cocaine <- t(seizures.cocaine)[i,]
  localM.month <- localmoran_abs(seizure.cocaine, nb.obj.cocaine, nsim=nperm, zero.policy=T, xx=NULL, alternative="two.sided")
  localM.month$LISA_C <- ifelse(localM.month$quadr=="High-High", 1,
                                ifelse(localM.month$quadr=="Low-Low", 2,
                                       ifelse(localM.month$quadr=="Low-High", 3, 4)))
  localM.month$LISA_C <- ifelse(localM.month$`Pr(folded) Sim` <= alpha, localM.month$LISA_C, 0)
  relative.MoransI[,(LISA_I.index+3*i-3):(LISA_I.index+3*i-1)] <- localM.month[,c(1,13,7)]
}

# write.csv(relative.MoransI, "cocaine other counts KNN5 R codes 999 two-sided (relative).csv", row.names=F)

absolute.b1.MoransI <- cocaine
x.bar_p <- mean(absolute.b1.MoransI$Jan_2018)
set.seed(300)
for (i in 1:48) {
  seizure.cocaine <- t(seizures.cocaine)[i,]
  localM.month <- localmoran_abs(seizure.cocaine, nb.obj.cocaine, nsim=nperm, zero.policy=T, xx=x.bar_p, alternative="two.sided")
  localM.month$LISA_C <- ifelse(localM.month$quadr=="High-High", 1,
                                ifelse(localM.month$quadr=="Low-Low", 2,
                                       ifelse(localM.month$quadr=="Low-High", 3, 4)))
  localM.month$LISA_C <- ifelse(localM.month$`Pr(folded) Sim` <= alpha, localM.month$LISA_C, 0)
  absolute.b1.MoransI[,(LISA_I.index+3*i-3):(LISA_I.index+3*i-1)] <- localM.month[,c(1,13,7)]
}

# write.csv(absolute.b1.MoransI, "cocaine other counts KNN5 R codes 999 two-sided (absolute base 1).csv", row.names=F)

absolute.t1.MoransI <- cocaine
x.bar_p <- mean(absolute.t1.MoransI$Jan_2018)
set.seed(500)
for (i in 1:48) {
  seizure.cocaine <- t(seizures.cocaine)[i,]
  localM.month <- localmoran_abs(seizure.cocaine, nb.obj.cocaine, nsim=nperm, zero.policy=T, xx=x.bar_p, alternative="two.sided")
  localM.month$LISA_C <- ifelse(localM.month$quadr=="High-High", 1,
                                ifelse(localM.month$quadr=="Low-Low", 2,
                                       ifelse(localM.month$quadr=="Low-High", 3, 4)))
  localM.month$LISA_C <- ifelse(localM.month$`Pr(folded) Sim` <= alpha, localM.month$LISA_C, 0)
  absolute.t1.MoransI[,(LISA_I.index+3*i-3):(LISA_I.index+3*i-1)] <- localM.month[,c(1,13,7)]
  x.bar_p <- mean(as.data.frame(cocaine)[,Jan_2018_index+i-1])
}

# write.csv(absolute.t1.MoransI, "cocaine other counts KNN5 R codes 999 two-sided (absolute base t_1).csv", row.names=F)

# Other Cocaine Plots
LISA.rel <- read.csv("cocaine other counts KNN5 R codes 999 two-sided (relative).csv") %>% as_tibble
LISA.abs.1 <- read.csv("cocaine other counts KNN5 R codes 999 two-sided (absolute base 1).csv") %>% as_tibble
LISA.abs.t_1 <- read.csv("cocaine other counts KNN5 R codes 999 two-sided (absolute base t_1).csv") %>% as_tibble

LISA_C.rel <- LISA.rel[,c(1:3, grep("LISA_C", names(LISA.rel)))]
LISA_C.abs.1 <- LISA.abs.1[,c(1:3, grep("LISA_C", names(LISA.rel)))]
LISA_C.abs.t_1 <- LISA.abs.t_1[,c(1:3, grep("LISA_C", names(LISA.rel)))]
county.names <- LISA.rel[, 1:3]

names(LISA_C.rel)[4:51] <- LISA_C_to_Month(names(LISA_C.rel)[4:51])
names(LISA_C.abs.1)[4:51] <- LISA_C_to_Month(names(LISA_C.abs.1)[4:51])
names(LISA_C.abs.t_1)[4:51] <- LISA_C_to_Month(names(LISA_C.abs.t_1)[4:51])


# Number of HL maps
n_HL_map <- function(row) {
  result <- ifelse(sum(is.na(row) == length(row)),
                   NA, sum(row==4, na.rm=T))
  return(result)
}
LISA_n_HL <- LISA_C.rel
for (year in 2018:2021) {
  LISA_C_index <- grep(year, names(LISA_n_HL))
  LISA_n_HL[[paste("n_HL", year, sep="_")]] <- apply(LISA_n_HL[,LISA_C_index], 1, n_HL_map)
}
LISA_n_HL <- LISA_n_HL %>% select(state:GEOID, n_HL_2018:n_HL_2021)
LISA_n_HL$n_HL_All <- apply(LISA_n_HL[, 4:7], 1, sum)
LISA_n_HL.map <- left_join(counties.obs[, c(1:2,6:7,10,12)], LISA_n_HL[,c(3:8)], by = "GEOID")
LISA_n_HL.map %>% pivot_longer(-c(long, lat, group, GEOID, state, county, n_HL_All), names_to="year", values_to="n_HL_year") %>% 
  ggplot(mapping = aes(long, lat, group = group, fill=n_HL_year)) +
  geom_polygon(color = "#000000", size = .05) +
  facet_wrap(. ~ year) +
  scale_fill_viridis_c(na.value="white") +
  labs(fill = "# of HL") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) -> n_HL_map
# ggsave(paste("other cocaine n_HL annual.jpg", sep=""), n_HL_map, width=20, height=15, units="cm")

LISA_n_HL.map %>% 
  ggplot(mapping = aes(long, lat, group = group, fill=n_HL_All)) +
  geom_polygon(color = "#000000", size = .05) +
  scale_fill_viridis_c(na.value="white") +
  labs(fill = "# of HL", title="# of HL 2018-2021 (other cocaine )") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) -> n_HL_map_all
# ggsave(paste("other cocaine n_HL 2018-2021.jpg", sep=""), n_HL_map_all, width=20, height=15, units="cm")

LISA_C.rel.map <- left_join(counties.obs[, c(1:2,7,10,12)], LISA_C.rel[,-(1:2)], by = "GEOID")
LISA_C.rel.map <- LISA_C.rel.map %>% pivot_longer(c(-long, -lat, -GEOID, -state, -county), names_to="Month_Year", values_to="LISA_C")
LISA_C.rel.map$LISA_C <- left_join(LISA_C.rel.map[, 6:7], ref, by="LISA_C")$label
LISA_C.rel.map$LISA_C <- factor(LISA_C.rel.map$LISA_C)
LISA_C.rel.map$Month_Year <- parse_date(LISA_C.rel.map$Month_Year, "%Y-%m")

for (year in 2018:2021) { # monthly maps of LISA_C.rel
  LISA_C.rel.map %>% filter(year(Month_Year) == year & month(Month_Year) %in% 1:4) %>% 
    ggplot(mapping = aes(long, lat, group = group, fill=LISA_C)) +
    geom_polygon(color = "#000000", size = .05) +
    facet_wrap(. ~ substr(Month_Year, 1, 7)) +
    scale_fill_manual(values = c("Insig"="grey60",
                                 "LL"="blue",
                                 "LH"="steelblue",
                                 "HL"="orange",
                                 "HH"="red"),
                      na.value = "white") +
    labs(fill = "LISA_C_rel") + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) -> LISA_C_map
  ggsave(paste("cocaine_other_counts_LISA_C_rel_", year, " two-sided (4 months).jpg", sep=""), LISA_C_map, width=20, height=15, units="cm")
}

LISA_C.abs1.map <- left_join(counties.obs[, c(1:2,7,10,12)], LISA_C.abs.1[,-(1:2)], by = "GEOID")
LISA_C.abs1.map <- LISA_C.abs1.map %>% pivot_longer(c(-long, -lat, -GEOID, -state, -county), names_to="Month_Year", values_to="LISA_C")
LISA_C.abs1.map$LISA_C <- left_join(LISA_C.abs1.map[, 6:7], ref, by="LISA_C")$label
LISA_C.abs1.map$LISA_C <- factor(LISA_C.abs1.map$LISA_C)
LISA_C.abs1.map$Month_Year <- parse_date(LISA_C.abs1.map$Month_Year, "%Y-%m")

for (year in 2018:2021) { # monthly maps of LISA_C.abs1
  LISA_C.abs1.map %>% filter(year(Month_Year) == year & month(Month_Year) %in% 1:4) %>% 
    ggplot(mapping = aes(long, lat, group = group, fill=LISA_C)) +
    geom_polygon(color = "#000000", size = .05) +
    facet_wrap(. ~ substr(Month_Year, 1, 7)) +
    scale_fill_manual(values = c("Insig"="grey60",
                                 "LL"="blue",
                                 "LH"="steelblue",
                                 "HL"="orange",
                                 "HH"="red"),
                      na.value = "white") +
    labs(fill = "LISA_C_abs1") + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) -> LISA_C_map
  ggsave(paste("cocaine_other_counts_LISA_C_abs1_", year, " two-sided (4 months).jpg", sep=""), LISA_C_map, width=20, height=15, units="cm")
}

LISA_C.abs_t1.map <- left_join(counties.obs[, c(1:2,7,10,12)], LISA_C.abs.t_1[,-(1:2)], by = "GEOID")
LISA_C.abs_t1.map <- LISA_C.abs_t1.map %>% pivot_longer(c(-long, -lat, -GEOID, -state, -county), names_to="Month_Year", values_to="LISA_C")
LISA_C.abs_t1.map$LISA_C <- left_join(LISA_C.abs_t1.map[, 6:7], ref, by="LISA_C")$label
LISA_C.abs_t1.map$LISA_C <- factor(LISA_C.abs_t1.map$LISA_C)
LISA_C.abs_t1.map$Month_Year <- parse_date(LISA_C.abs_t1.map$Month_Year, "%Y-%m")

for (year in 2018:2021) { # monthly maps of LISA_C.abs_t1
  LISA_C.abs_t1.map %>% filter(year(Month_Year) == year & month(Month_Year) %in% 1:4) %>% 
    ggplot(mapping = aes(long, lat, group = group, fill=LISA_C)) +
    geom_polygon(color = "#000000", size = .05) +
    facet_wrap(. ~ substr(Month_Year, 1, 7)) +
    scale_fill_manual(values = c("Insig"="grey60",
                                 "LL"="blue",
                                 "LH"="steelblue",
                                 "HL"="orange",
                                 "HH"="red"),
                      na.value = "white") +
    labs(fill = "LISA_C_abs_t1") + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) -> LISA_C_map
  ggsave(paste("cocaine_other_counts_LISA_C_abs_t1_", year, " two-sided (4 months).jpg", sep=""), LISA_C_map, width=20, height=15, units="cm")
}