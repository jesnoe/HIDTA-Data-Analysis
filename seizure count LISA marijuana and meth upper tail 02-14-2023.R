library(fpp2)
library(spdep)
library(readxl)
library(urbnmapr)
library(tidyverse)
library(gridExtra)
library(lubridate)

marijuana <- read.csv("marijuana count HIDTA (02-14-2023).csv") %>% as_tibble
LISA3 <- read.csv("CountyKNN3.csv") %>% as_tibble %>% arrange(GEOID)

counties.obs <- counties
names(counties.obs)[7] <- "GEOID"
counties.obs <- counties.obs %>% rename(state=state_name, county=county_name)
counties.obs$GEOID <- as.numeric(counties.obs$GEOID)
counties.obs <- counties.obs %>% filter(!(state %in% c("Alaska", "Hawaii")))
coords.marijuana <- counties.obs %>% filter(GEOID %in% marijuana$GEOID) %>%
  group_by(GEOID) %>% summarise(x=mean(long), y=mean(lat))
GEOIDS.marijuana <- coords.marijuana$GEOID
coords.marijuana <- coords.marijuana[,-1]
nb_marijuana <- knn2nb(knearneigh(coords.marijuana, k=5), row.names=GEOIDS.marijuana)

alpha <- 0.05
nperm <- 999

## For marijuana
marijuana <- cbind(marijuana, matrix(0, nrow(marijuana), 48*3)) %>% as_tibble

names(marijuana)[(which(names(marijuana) == "1"):which(names(marijuana) == "144"))] <- names(LISA3)[71:214]
Jan_2018_index <- grep("Jan_2018", names(marijuana))[1]
LISA_I.index <- grep("LISA_I", names(marijuana))[1]
nb.obj.marijuana <- nb2listw(nb_marijuana, style="B")
nrow(marijuana) # 1218
seizures.marijuana <- marijuana[, 5:52]

grep("LISA_I", names(marijuana))[1] # 53

relative.MoransI <- marijuana
set.seed(100)
for (i in 1:48) {
  seizure.marijuana <- t(seizures.marijuana)[i,]
  localM.month <- localmoran_abs(seizure.marijuana, nb.obj.marijuana, nsim=nperm, zero.policy=T, xx=NULL, alternative="greater")
  localM.month$LISA_C <- ifelse(localM.month$quadr=="High-High", 1,
                                ifelse(localM.month$quadr=="Low-Low", 2,
                                       ifelse(localM.month$quadr=="Low-High", 3, 4)))
  localM.month$LISA_C <- ifelse(localM.month$`Pr(folded) Sim` <= alpha, localM.month$LISA_C, 0)
  relative.MoransI[,(LISA_I.index+3*i-3):(LISA_I.index+3*i-1)] <- localM.month[,c(1,13,7)]
}

# write.csv(relative.MoransI, "marijuana count KNN5 R codes 999 upper tail (relative).csv", row.names=F)

absolute.b1.MoransI <- marijuana
x.bar_p <- mean(absolute.b1.MoransI$Jan_2018)
set.seed(300)
for (i in 1:48) {
  seizure.marijuana <- t(seizures.marijuana)[i,]
  localM.month <- localmoran_abs(seizure.marijuana, nb.obj.marijuana, nsim=nperm, zero.policy=T, xx=x.bar_p, alternative="greater")
  localM.month$LISA_C <- ifelse(localM.month$quadr=="High-High", 1,
                                ifelse(localM.month$quadr=="Low-Low", 2,
                                       ifelse(localM.month$quadr=="Low-High", 3, 4)))
  localM.month$LISA_C <- ifelse(localM.month$`Pr(folded) Sim` <= alpha, localM.month$LISA_C, 0)
  absolute.b1.MoransI[,(LISA_I.index+3*i-3):(LISA_I.index+3*i-1)] <- localM.month[,c(1,13,7)]
}

# write.csv(absolute.b1.MoransI, "marijuana count KNN5 R codes 999 upper tail (absolute base 1).csv", row.names=F)

absolute.t1.MoransI <- marijuana
x.bar_p <- mean(absolute.t1.MoransI$Jan_2018)
set.seed(500)
for (i in 1:48) {
  seizure.marijuana <- t(seizures.marijuana)[i,]
  localM.month <- localmoran_abs(seizure.marijuana, nb.obj.marijuana, nsim=nperm, zero.policy=T, xx=x.bar_p, alternative="greater")
  localM.month$LISA_C <- ifelse(localM.month$quadr=="High-High", 1,
                                ifelse(localM.month$quadr=="Low-Low", 2,
                                       ifelse(localM.month$quadr=="Low-High", 3, 4)))
  localM.month$LISA_C <- ifelse(localM.month$`Pr(folded) Sim` <= alpha, localM.month$LISA_C, 0)
  absolute.t1.MoransI[,(LISA_I.index+3*i-3):(LISA_I.index+3*i-1)] <- localM.month[,c(1,13,7)]
  x.bar_p <- mean(as.data.frame(marijuana)[,Jan_2018_index+i-1])
}

# write.csv(absolute.t1.MoransI, "marijuana count KNN5 R codes 999 upper tail (absolute base t_1).csv", row.names=F)

# marijuana Plots
LISA.rel <- read.csv("marijuana count KNN5 R codes 999 upper tail (relative).csv") %>% as_tibble
LISA.abs.1 <- read.csv("marijuana count KNN5 R codes 999 upper tail (absolute base 1).csv") %>% as_tibble
LISA.abs.t_1 <- read.csv("marijuana count KNN5 R codes 999 upper tail (absolute base t_1).csv") %>% as_tibble

LISA_C.rel <- LISA.rel[,c(1:3, grep("LISA_C", names(LISA.rel)))]
LISA_C.abs.1 <- LISA.abs.1[,c(1:3, grep("LISA_C", names(LISA.rel)))]
LISA_C.abs.t_1 <- LISA.abs.t_1[,c(1:3, grep("LISA_C", names(LISA.rel)))]
county.names <- LISA.rel[, 1:3]

# cbind(LISA_C.rel[,1:3], n_HH=apply(LISA_C.rel[, -(1:3)], 1, sum)) %>% write.csv("marijuana count rel HH summary (02-02-2023).csv", row.names=F)

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


LISA.rel.map <- left_join(counties.obs[, c(1:2,6:7,10,12)], LISA.rel[, c(3, 5:52)], by = "GEOID")
LISA.rel.map <- LISA.rel.map %>% 
  pivot_longer(c(-long, -lat, -group, -GEOID, -state, -county), names_to="Month_Year", values_to="seizure_counts") %>% 
  mutate(Year=as.integer(substr(Month_Year, 5,8)))
LISA.rel.map <- left_join(LISA.rel.map, populations[, -(1:2)], by=c("GEOID", "Year"))
LISA.rel.map$counts_per_pop <- LISA.rel.map$seizure_counts/LISA.rel.map$population
LISA.rel.map$Month_Year <- parse_date(LISA.rel.map$Month_Year, "%b_%Y")

for (year in 2018:2021) { # monthly maps of seizure counts
  LISA.rel.map %>% filter(year(Month_Year) == year & month(Month_Year) %in% 1:4) %>% 
    ggplot(mapping = aes(long, lat, group = group, fill=seizure_counts)) +
    geom_polygon(color = "#000000", size = .05) +
    facet_wrap(. ~ substr(Month_Year, 1, 7)) +
    scale_fill_viridis_c(na.value="white") +
    labs(fill = "Seizure Counts") + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) -> seizure_counts_map
  ggsave(paste("marijuana seizure_counts_", year, " (4 months).jpg", sep=""), seizure_counts_map, width=20, height=15, units="cm")
  
  LISA.rel.map %>% filter(grepl(year, Month_Year) & month(Month_Year) %in% 1:4) %>% 
    ggplot(mapping = aes(long, lat, group = group, fill=counts_per_pop)) +
    geom_polygon(color = "#000000", size = .05) +
    facet_wrap(. ~ substr(Month_Year, 1, 7)) +
    scale_fill_viridis_c(na.value="white") +
    labs(fill = "Counts Per Pop.") + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) -> counts_per_pop.map
  ggsave(paste("marijuana counts_per_pop_", year, " (4 months).jpg", sep=""), counts_per_pop.map, width=20, height=15, units="cm")
}

# Number of HH Maps
n_HH_map <- function(row) {
  result <- ifelse(sum(is.na(row) == length(row)),
                   NA, sum(row, na.rm=T))
  return(result)
}
LISA_n_HH <- LISA_C.rel
for (year in 2018:2021) {
  LISA_C_index <- grep(year, names(LISA_n_HH))
  LISA_n_HH[[paste("n_HH", year, sep="_")]] <- apply(LISA_n_HH[,LISA_C_index], 1, n_HH_map)
}
LISA_n_HH <- LISA_n_HH %>% select(state:GEOID, n_HH_2018:n_HH_2021)
LISA_n_HH$n_HH_All <- apply(LISA_n_HH[, 4:7], 1, sum)
LISA_n_HH.map <- left_join(counties.obs[, c(1:2,6:7,10,12)], LISA_n_HH[,c(3:8)], by = "GEOID")
LISA_n_HH.map %>% pivot_longer(-c(long, lat, group, GEOID, state, county, n_HH_All), names_to="year", values_to="n_HH_year") %>% 
  ggplot(mapping = aes(long, lat, group = group, fill=n_HH_year)) +
  geom_polygon(color = "#000000", size = .05) +
  facet_wrap(. ~ year) +
  scale_fill_viridis_c(na.value="white") +
  labs(fill = "# of HH") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) -> n_HH_map
# ggsave(paste("marijuana n_HH annual.jpg", sep=""), n_HH_map, width=20, height=15, units="cm")

LISA_n_HH.map %>% 
  ggplot(mapping = aes(long, lat, group = group, fill=n_HH_All)) +
  geom_polygon(color = "#000000", size = .05) +
  scale_fill_viridis_c(na.value="white") +
  labs(fill = "# of HH", title="# of HH 2018-2021 (marijuana)") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) -> n_HH_map_all
# ggsave(paste("marijuana n_HH 2018-2021.jpg", sep=""), n_HH_map_all, width=20, height=15, units="cm")


LISA_C.rel.map <- left_join(counties.obs[, c(1:2,6:7,10,12)], LISA_C.rel[,-(1:2)], by = "GEOID")
LISA_C.rel.map <- LISA_C.rel.map %>% pivot_longer(c(-long, -lat, -group, -GEOID, -state, -county), names_to="Month_Year", values_to="LISA_C")
ref <- tibble(LISA_C=0:4, label=c("Insig", "HH", "LL", "LH", "HL"))
LISA_C.rel.map$LISA_C <- left_join(LISA_C.rel.map[, 7:8], ref, by="LISA_C")$label
LISA_C.rel.map$LISA_C <- factor(LISA_C.rel.map$LISA_C)
LISA_C.rel.map$Month_Year <- parse_date(LISA_C.rel.map$Month_Year, "%Y-%m")

for (year in 2018:2021) { # monthly maps of LISA_C.rel
  LISA_C.rel.map %>% filter(year(Month_Year) == year & month(Month_Year) %in% 1:4) %>% 
    ggplot(mapping = aes(long, lat, group = group, fill=LISA_C)) +
    geom_polygon(color = "#000000", size = .05) +
    facet_wrap(. ~ substr(Month_Year, 1, 7)) +
    scale_fill_manual(values = c("Insig"="grey60",
                                 "HH"="red"),
                      na.value = "white") +
    labs(fill = "LISA_C_rel") + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) -> LISA_C_map
  ggsave(paste("marijuana_Count_LISA_C_rel_upper_tail_", year, " (4 months).jpg", sep=""), LISA_C_map, width=20, height=15, units="cm")
}

LISA_C.abs1.map <- left_join(counties.obs[, c(1:2,6:7,10,12)], LISA_C.abs.1[,-(1:2)], by = "GEOID")
LISA_C.abs1.map <- LISA_C.abs1.map %>% pivot_longer(c(-long, -lat, -group, -GEOID, -state, -county), names_to="Month_Year", values_to="LISA_C")
LISA_C.abs1.map$LISA_C <- left_join(LISA_C.abs1.map[, 7:8], ref, by="LISA_C")$label
LISA_C.abs1.map$LISA_C <- factor(LISA_C.abs1.map$LISA_C)
LISA_C.abs1.map$Month_Year <- parse_date(LISA_C.abs1.map$Month_Year, "%Y-%m")

for (year in 2018:2021) { # monthly maps of LISA_C.abs1
  LISA_C.abs1.map %>% filter(year(Month_Year) == year & month(Month_Year) %in% 1:4) %>% 
    ggplot(mapping = aes(long, lat, group = group, fill=LISA_C)) +
    geom_polygon(color = "#000000", size = .05) +
    facet_wrap(. ~ substr(Month_Year, 1, 7)) +
    scale_fill_manual(values = c("Insig"="grey60",
                                 "HH"="red"),
                      na.value = "white") +
    labs(fill = "LISA_C_abs1") + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) -> LISA_C_map
  ggsave(paste("marijuana_Count_LISA_C_abs1_upper_tail_", year, " (4 months).jpg", sep=""), LISA_C_map, width=20, height=15, units="cm")
}

LISA_C.abs_t1.map <- left_join(counties.obs[, c(1:2,6:7,10,12)], LISA_C.abs.t_1[,-(1:2)], by = "GEOID")
LISA_C.abs_t1.map <- LISA_C.abs_t1.map %>% pivot_longer(c(-long, -lat, -group, -GEOID, -state, -county), names_to="Month_Year", values_to="LISA_C")
LISA_C.abs_t1.map$LISA_C <- left_join(LISA_C.abs_t1.map[, 7:8], ref, by="LISA_C")$label
LISA_C.abs_t1.map$LISA_C <- factor(LISA_C.abs_t1.map$LISA_C)
LISA_C.abs_t1.map$Month_Year <- parse_date(LISA_C.abs_t1.map$Month_Year, "%Y-%m")

for (year in 2018:2021) { # monthly maps of LISA_C.abs_t1
  LISA_C.abs_t1.map %>% filter(year(Month_Year) == year & month(Month_Year) %in% 1:4) %>% 
    ggplot(mapping = aes(long, lat, group = group, fill=LISA_C)) +
    geom_polygon(color = "#000000", size = .05) +
    facet_wrap(. ~ substr(Month_Year, 1, 7)) +
    scale_fill_manual(values = c("Insig"="grey60",
                                 "HH"="red"),
                      na.value = "white") +
    labs(fill = "LISA_C_abs_t1") + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) -> LISA_C_map
  ggsave(paste("marijuana_Count_LISA_C_abs_t1_upper_tail_", year, " (4 months).jpg", sep=""), LISA_C_map, width=20, height=15, units="cm")
}



### For meth
meth <- read.csv("meth count HIDTA (02-14-2023).csv") %>% as_tibble
coords.meth <- counties.obs %>% filter(GEOID %in% meth$GEOID) %>%
  group_by(GEOID) %>% summarise(x=mean(long), y=mean(lat))
GEOIDS.meth <- coords.meth$GEOID
coords.meth <- coords.meth[,-1]
nb_meth <- knn2nb(knearneigh(coords.meth, k=5), row.names=GEOIDS.meth)

alpha <- 0.05
nperm <- 999


meth <- cbind(meth, matrix(0, nrow(meth), 48*3)) %>% as_tibble
names(meth)[(which(names(meth) == "1"):which(names(meth) == "144"))] <- names(LISA3)[71:214]
Jan_2018_index <- grep("Jan_2018", names(meth))[1]
LISA_I.index <- grep("LISA_I", names(meth))[1]
nb.obj.meth <- nb2listw(nb_meth, style="B")
nrow(meth) # 1318
seizures.meth <- meth[, 5:52]

relative.MoransI <- meth
set.seed(100)
for (i in 1:48) {
  seizure.meth <- t(seizures.meth)[i,]
  localM.month <- localmoran_abs(seizure.meth, nb.obj.meth, nsim=nperm, zero.policy=T, xx=NULL, alternative="greater")
  localM.month$LISA_C <- ifelse(localM.month$quadr=="High-High", 1,
                                ifelse(localM.month$quadr=="Low-Low", 2,
                                       ifelse(localM.month$quadr=="Low-High", 3, 4)))
  localM.month$LISA_C <- ifelse(localM.month$`Pr(folded) Sim` <= alpha, localM.month$LISA_C, 0)
  relative.MoransI[,(LISA_I.index+3*i-3):(LISA_I.index+3*i-1)] <- localM.month[,c(1,13,7)]
}

# write.csv(relative.MoransI, "meth counts KNN5 R codes 999 upper tail (relative).csv", row.names=F)

absolute.b1.MoransI <- meth
x.bar_p <- mean(absolute.b1.MoransI$Jan_2018)
set.seed(300)
for (i in 1:48) {
  seizure.meth <- t(seizures.meth)[i,]
  localM.month <- localmoran_abs(seizure.meth, nb.obj.meth, nsim=nperm, zero.policy=T, xx=x.bar_p, alternative="greater")
  localM.month$LISA_C <- ifelse(localM.month$quadr=="High-High", 1,
                                ifelse(localM.month$quadr=="Low-Low", 2,
                                       ifelse(localM.month$quadr=="Low-High", 3, 4)))
  localM.month$LISA_C <- ifelse(localM.month$`Pr(folded) Sim` <= alpha, localM.month$LISA_C, 0)
  absolute.b1.MoransI[,(LISA_I.index+3*i-3):(LISA_I.index+3*i-1)] <- localM.month[,c(1,13,7)]
}

# write.csv(absolute.b1.MoransI, "meth counts KNN5 R codes 999 upper tail (absolute base 1).csv", row.names=F)

absolute.t1.MoransI <- meth
x.bar_p <- mean(absolute.t1.MoransI$Jan_2018)
set.seed(500)
for (i in 1:48) {
  seizure.meth <- t(seizures.meth)[i,]
  localM.month <- localmoran_abs(seizure.meth, nb.obj.meth, nsim=nperm, zero.policy=T, xx=x.bar_p, alternative="greater")
  localM.month$LISA_C <- ifelse(localM.month$quadr=="High-High", 1,
                                ifelse(localM.month$quadr=="Low-Low", 2,
                                       ifelse(localM.month$quadr=="Low-High", 3, 4)))
  localM.month$LISA_C <- ifelse(localM.month$`Pr(folded) Sim` <= alpha, localM.month$LISA_C, 0)
  absolute.t1.MoransI[,(LISA_I.index+3*i-3):(LISA_I.index+3*i-1)] <- localM.month[,c(1,13,7)]
  x.bar_p <- mean(as.data.frame(meth)[,Jan_2018_index+i-1])
}

# write.csv(absolute.t1.MoransI, "meth counts KNN5 R codes 999 upper tail (absolute base t_1).csv", row.names=F)

# meth Plots
LISA.rel <- read.csv("meth counts KNN5 R codes 999 upper tail (relative).csv") %>% as_tibble
LISA.abs.1 <- read.csv("meth counts KNN5 R codes 999 upper tail (absolute base 1).csv") %>% as_tibble
LISA.abs.t_1 <- read.csv("meth counts KNN5 R codes 999 upper tail (absolute base t_1).csv") %>% as_tibble

LISA_C.rel <- LISA.rel[,c(1:3, grep("LISA_C", names(LISA.rel)))]
LISA_C.abs.1 <- LISA.abs.1[,c(1:3, grep("LISA_C", names(LISA.rel)))]
LISA_C.abs.t_1 <- LISA.abs.t_1[,c(1:3, grep("LISA_C", names(LISA.rel)))]
county.names <- LISA.rel[, 1:3]

names(LISA_C.rel)[4:51] <- LISA_C_to_Month(names(LISA_C.rel)[4:51])
names(LISA_C.abs.1)[4:51] <- LISA_C_to_Month(names(LISA_C.abs.1)[4:51])
names(LISA_C.abs.t_1)[4:51] <- LISA_C_to_Month(names(LISA_C.abs.t_1)[4:51])

LISA.rel.map <- left_join(counties.obs[, c(1:2,6:7,10,12)], LISA.rel[, c(3, 5:52)], by = "GEOID")
LISA.rel.map <- LISA.rel.map %>% 
  pivot_longer(c(-long, -lat, -group, -GEOID, -state, -county), names_to="Month_Year", values_to="seizure_counts") %>% 
  mutate(Year=as.integer(substr(Month_Year, 5,8)))
LISA.rel.map <- left_join(LISA.rel.map, populations[, -(1:2)], by=c("GEOID", "Year"))
LISA.rel.map$Month_Year <- parse_date(LISA.rel.map$Month_Year, "%b_%Y")
LISA.rel.map$counts_per_pop <- LISA.rel.map$seizure_counts/LISA.rel.map$population

for (year in 2018:2021) { # monthly maps of seizure counts
  LISA.rel.map %>% filter(year(Month_Year) == year & month(Month_Year) %in% 1:4) %>% 
    ggplot(mapping = aes(long, lat, group = group, fill=seizure_counts)) +
    geom_polygon(color = "#000000", size = .05) +
    facet_wrap(. ~ substr(Month_Year, 1, 7)) +
    # coord_equal() +
    scale_fill_viridis_c(na.value="white") +
    labs(fill = "Seizure Counts") + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) -> seizure_counts_map
  ggsave(paste("meth seizure_counts_", year, " (4 months).jpg", sep=""), seizure_counts_map, width=20, height=15, units="cm")
  
  LISA.rel.map %>% filter(grepl(year, Month_Year) & month(Month_Year) %in% 1:4) %>% 
    ggplot(mapping = aes(long, lat, group = group, fill=counts_per_pop)) +
    geom_polygon(color = "#000000", size = .05) +
    facet_wrap(. ~ substr(Month_Year, 1, 7)) +
    # coord_equal() +
    scale_fill_viridis_c(na.value="white") +
    labs(fill = "Counts Per Pop.") + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) -> counts_per_pop.map
  ggsave(paste("meth counts_per_pop_", year, " (4 months).jpg", sep=""), counts_per_pop.map, width=20, height=15, units="cm")
}

# Number of HH Maps
n_HH_map <- function(row) {
  result <- ifelse(sum(is.na(row) == length(row)),
                   NA, sum(row, na.rm=T))
  return(result)
}
LISA_n_HH <- LISA_C.rel
for (year in 2018:2021) {
  LISA_C_index <- grep(year, names(LISA_n_HH))
  LISA_n_HH[[paste("n_HH", year, sep="_")]] <- apply(LISA_n_HH[,LISA_C_index], 1, n_HH_map)
}
LISA_n_HH <- LISA_n_HH %>% select(state:GEOID, n_HH_2018:n_HH_2021)
LISA_n_HH$n_HH_All <- apply(LISA_n_HH[, 4:7], 1, sum)
LISA_n_HH.map <- left_join(counties.obs[, c(1:2,6:7,10,12)], LISA_n_HH[,c(3:8)], by = "GEOID")
LISA_n_HH.map %>% pivot_longer(-c(long, lat, group, GEOID, state, county, n_HH_All), names_to="year", values_to="n_HH_year") %>% 
  ggplot(mapping = aes(long, lat, group = group, fill=n_HH_year)) +
  geom_polygon(color = "#000000", size = .05) +
  facet_wrap(. ~ year) +
  scale_fill_viridis_c(na.value="white") +
  labs(fill = "# of HH") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) -> n_HH_map
# ggsave(paste(" meth n_HH annual.jpg", sep=""), n_HH_map, width=20, height=15, units="cm")

LISA_n_HH.map %>% 
  ggplot(mapping = aes(long, lat, group = group, fill=n_HH_All)) +
  geom_polygon(color = "#000000", size = .05) +
  scale_fill_viridis_c(na.value="white") +
  labs(fill = "# of HH", title="# of HH 2018-2021 (meth)") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) -> n_HH_map_all
# ggsave(paste(" meth n_HH 2018-2021.jpg", sep=""), n_HH_map_all, width=20, height=15, units="cm")

LISA_C.rel.map <- left_join(counties.obs[, c(1:2,6:7,10,12)], LISA_C.rel[,-(1:2)], by = "GEOID")
LISA_C.rel.map <- LISA_C.rel.map %>% pivot_longer(c(-long, -lat, -group, -GEOID, -state, -county), names_to="Month_Year", values_to="LISA_C")
LISA_C.rel.map$LISA_C <- left_join(LISA_C.rel.map[, 7:8], ref, by="LISA_C")$label
LISA_C.rel.map$LISA_C <- factor(LISA_C.rel.map$LISA_C)
LISA_C.rel.map$Month_Year <- parse_date(LISA_C.rel.map$Month_Year, "%Y-%m")

for (year in 2018:2021) { # monthly maps of LISA_C.rel
  LISA_C.rel.map %>% filter(year(Month_Year) == year & month(Month_Year) %in% 1:4) %>% 
    ggplot(mapping = aes(long, lat, group = group, fill=LISA_C)) +
    geom_polygon(color = "#000000", size = .05) +
    facet_wrap(. ~ substr(Month_Year, 1, 7)) +
    scale_fill_manual(values = c("Insig"="grey60",
                                 "HH"="red"),
                      na.value = "white") +
    labs(fill = "LISA_C_rel") + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) -> LISA_C_map
  ggsave(paste("meth_counts_LISA_C_rel_", year, " upper tail (4 months).jpg", sep=""), LISA_C_map, width=20, height=15, units="cm")
}

LISA_C.abs1.map <- left_join(counties.obs[, c(1:2,6:7,10,12)], LISA_C.abs.1[,-(1:2)], by = "GEOID")
LISA_C.abs1.map <- LISA_C.abs1.map %>% pivot_longer(c(-long, -lat, -group, -GEOID, -state, -county), names_to="Month_Year", values_to="LISA_C")
LISA_C.abs1.map$LISA_C <- left_join(LISA_C.abs1.map[, 7:8], ref, by="LISA_C")$label
LISA_C.abs1.map$LISA_C <- factor(LISA_C.abs1.map$LISA_C)
LISA_C.abs1.map$Month_Year <- parse_date(LISA_C.abs1.map$Month_Year, "%Y-%m")

for (year in 2018:2021) { # monthly maps of LISA_C.abs1
  LISA_C.abs1.map %>% filter(year(Month_Year) == year & month(Month_Year) %in% 1:4) %>% 
    ggplot(mapping = aes(long, lat, group = group, fill=LISA_C)) +
    geom_polygon(color = "#000000", size = .05) +
    facet_wrap(. ~ substr(Month_Year, 1, 7)) +
    scale_fill_manual(values = c("Insig"="grey60",
                                 "HH"="red"),
                      na.value = "white") +
    labs(fill = "LISA_C_abs1") + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) -> LISA_C_map
  ggsave(paste("meth_counts_LISA_C_abs1_", year, " upper tail (4 months).jpg", sep=""), LISA_C_map, width=20, height=15, units="cm")
}

LISA_C.abs_t1.map <- left_join(counties.obs[, c(1:2,6:7,10,12)], LISA_C.abs.t_1[,-(1:2)], by = "GEOID")
LISA_C.abs_t1.map <- LISA_C.abs_t1.map %>% pivot_longer(c(-long, -lat, -group, -GEOID, -state, -county), names_to="Month_Year", values_to="LISA_C")
LISA_C.abs_t1.map$LISA_C <- left_join(LISA_C.abs_t1.map[, 7:8], ref, by="LISA_C")$label
LISA_C.abs_t1.map$LISA_C <- factor(LISA_C.abs_t1.map$LISA_C)
LISA_C.abs_t1.map$Month_Year <- parse_date(LISA_C.abs_t1.map$Month_Year, "%Y-%m")

for (year in 2018:2021) { # monthly maps of LISA_C.abs_t1
  LISA_C.abs_t1.map %>% filter(year(Month_Year) == year & month(Month_Year) %in% 1:4) %>% 
    ggplot(mapping = aes(long, lat, group = group, fill=LISA_C)) +
    geom_polygon(color = "#000000", size = .05) +
    facet_wrap(. ~ substr(Month_Year, 1, 7)) +
    scale_fill_manual(values = c("Insig"="grey60",
                                 "HH"="red"),
                      na.value = "white") +
    labs(fill = "LISA_C_abs_t1") + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) -> LISA_C_map
  ggsave(paste("meth_counts_LISA_C_abs_t1_", year, " upper tail (4 months).jpg", sep=""), LISA_C_map, width=20, height=15, units="cm")
}