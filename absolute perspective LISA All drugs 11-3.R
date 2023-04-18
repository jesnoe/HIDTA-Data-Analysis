library(fpp2)
library(spdep)
library(readxl)
library(maps)
library(urbnmapr)
library(tidyverse)
library(rnaturalearth)

counties$county_fips <- as.numeric(counties$county_fips)

LISA3 <- read.csv("AllKGLISAResults.csv") %>% as_tibble %>% arrange(GEOID)
state.fips <- read.csv("us-state-ansi-fips.csv")
names(state.fips)[1] <- "state"
names(state.fips)[2] <- "STATEFP"
LISA3 <- merge(LISA3, state.fips[,1:2], by="STATEFP") %>% as_tibble
LISA3 <- LISA3[,c(1, 212, 2:211)]

counties.obs <- counties %>% filter(county_fips %in% unique(LISA3$GEOID))
names(counties.obs)[7] <- "GEOID"

grep("Jan_2018", names(LISA3)) # 21
grep("Dec_2021", names(LISA3)) # 68

plot(ts((LISA3 %>% filter(state=="Kentucky"))[, 21:68] %>% apply(2, sum), start=c(2018, 1), frequency=12),
     xlab="Year",
     ylab="Seizure Weights in Kentucky",
     pch=19, type="b")

# Counties with 0 kg seizure for 48 months are unobserved
LISA3 <- LISA3[LISA3[,21:68] %>% apply(1, sum) > 0,]
nrow(LISA3) # 1422, 630 more counties are observed compared to cocaine only data (792 counties)

seizures.state <- (LISA3 %>% group_by(state))[,c(2, 6, 8, 21:68)]
monthly.seizures <- data.frame()
for (state.name in unique(seizures.state$state)) {
  seizures.i <- seizures.state %>% filter(state == state.name)
  monthly.seizures <- rbind(monthly.seizures, c(state.name, apply(seizures.i[,-(1:3)], 2, sum)))
}
names(monthly.seizures) <- names(seizures.state)[-(2:3)]
monthly.seizures[,-1] <- apply(monthly.seizures[,-1], 2, as.numeric)
monthly.seizures$avg.seizure <- apply(monthly.seizures[,-1], 1, mean)
monthly.seizures %>% select(state, avg.seizure) %>% arrange(desc(avg.seizure))# %>% write.csv("All Drugs Average Seizures.csv", row.names=F)


coords.US <- counties.obs %>% group_by(GEOID) %>% summarise(x=mean(long), y=mean(lat))
GEOIDS <- coords.US$GEOID
coords.US <- coords.US[,-1]
LISA3_nb <- knn2nb(knearneigh(coords.US, k=3), row.names=GEOIDS)
LISA5_nb <- knn2nb(knearneigh(coords.US, k=5), row.names=GEOIDS)
LISA10_nb <- knn2nb(knearneigh(coords.US, k=10), row.names=GEOIDS)



# oopar <- par(mfrow=c(1,3), mar=c(1,1,1,1)+0.1)
# plot(spdf_US, border="grey60")
# plot(LISA3_nb, coords.US, add=TRUE, pch=".")
# text(bbox(Syracuse)[1,1], bbox(Syracuse)[2,2], labels="k=3", cex=0.8)
# plot(spdf_US, border="grey60")
# plot(Sy5_nb, coords, add=TRUE, pch=".")
# text(bbox(Syracuse)[1,1], bbox(Syracuse)[2,2], labels="k=5", cex=0.8)
# plot(Syracuse, border="grey60")
# plot(Sy10_nb, coords, add=TRUE, pch=".")
# text(bbox(Syracuse)[1,1], bbox(Syracuse)[2,2], labels="k=10", cex=0.8)


## No California
LISA3.XCal <- LISA3 %>% filter(state != "California")
coords.XCal <- counties.obs %>% filter(GEOID %in% LISA3.XCal$GEOID) %>%
  group_by(GEOID) %>% summarise(x=mean(long), y=mean(lat))
GEOIDS.XCal <- coords.XCal$GEOID
coords.XCal <- coords.XCal[,-1]
LISA3_nb_XCal <- knn2nb(knearneigh(coords.XCal, k=3), row.names=GEOIDS.XCal)

alpha <- 0.05
nperm <- 999
LISA_I.index <- grep("LISA_I", names(LISA3.XCal))[1]-1
nb.obj.XCal <- nb2listw(LISA3_nb_XCal, style="B")
nrow(LISA3.XCal) # 1386, 1422-1386=36 counties in Cal removed
seizures.XCal <- LISA3.XCal[, 21:68]

grep("LISA_I", names(LISA3.XCal))[1] # 69

relative.MoransI <- LISA3.XCal
set.seed(100)
for (i in 1:48) {
  seizure.XCal <- t(seizures.XCal)[i,]
  localM.month <- localmoran_abs(seizure.XCal, nb.obj.XCal, nsim=nperm, zero.policy=T, xx=NULL)
  localM.month$LISA_C <- ifelse(localM.month$mean=="High-High", 1,
                                ifelse(localM.month$mean=="Low-Low", 2,
                                       ifelse(localM.month$mean=="Low-High", 3, 4)))
  localM.month$LISA_C <- ifelse(localM.month$`Pr(folded) Sim` <= alpha, localM.month$LISA_C, 0)
  relative.MoransI[,(69+3*i-3):(69+3*i-1)] <- localM.month[,c(1,13,7)]
}

# write.csv(relative.MoransI, "CountyKNN3 XCal R codes 999 (relative).csv", row.names=F)

absolute.b1.MoransI <- LISA3.XCal
x.bar_p <- mean(absolute.b1.MoransI$Jan_2018)
set.seed(300)
for (i in 1:48) {
  seizure.XCal <- t(seizures.XCal)[i,]
  localM.month <- localmoran_abs(seizure.XCal, nb.obj.XCal, nsim=nperm, zero.policy=T, xx=NULL)
  localM.month$LISA_C <- ifelse(localM.month$mean=="High-High", 1,
                                ifelse(localM.month$mean=="Low-Low", 2,
                                       ifelse(localM.month$mean=="Low-High", 3, 4)))
  localM.month$LISA_C <- ifelse(localM.month$`Pr(folded) Sim` <= alpha, localM.month$LISA_C, 0)
  absolute.b1.MoransI[,(69+3*i-3):(69+3*i-1)] <- localM.month[,c(1,13,7)]
}

# write.csv(absolute.b1.MoransI, "CountyKNN3 XCal R codes 999 (absolute base 1).csv", row.names=F)

absolute.t1.MoransI <- LISA3.XCal
x.bar_p <- mean(absolute.t1.MoransI$Jan_2018)
set.seed(500)
for (i in 1:48) {
  seizure.XCal <- t(seizures.XCal)[i,]
  localM.month <- localmoran_abs(seizure.XCal, nb.obj.XCal, nsim=nperm, zero.policy=T, xx=NULL)
  localM.month$LISA_C <- ifelse(localM.month$mean=="High-High", 1,
                                ifelse(localM.month$mean=="Low-Low", 2,
                                       ifelse(localM.month$mean=="Low-High", 3, 4)))
  localM.month$LISA_C <- ifelse(localM.month$`Pr(folded) Sim` <= alpha, localM.month$LISA_C, 0)
  absolute.t1.MoransI[,(69+3*i-3):(69+3*i-1)] <- localM.month[,c(1,13,7)]
  x.bar_p <- mean(as.data.frame(LISA3.XCal)[,21+i-1])
}

# write.csv(absolute.t1.MoransI, "CountyKNN3 XCal R codes 999 (absolute base t_1).csv", row.names=F)


# XCal Plots
LISA.rel <- read.csv("CountyKNN3 XCal R codes 999 (relative).csv") %>% as_tibble
LISA.abs.1 <- read.csv("CountyKNN3 XCal R codes 999 (absolute base 1).csv") %>% as_tibble
LISA.abs.t_1 <- read.csv("CountyKNN3 XCal R codes 999 (absolute base t_1).csv") %>% as_tibble

LISA.rel[,grep("LISA_I", names(LISA.rel))]; LISA.abs.1[,grep("LISA_I", names(LISA.rel))]; LISA.abs.t_1[,grep("LISA_I", names(LISA.rel))]

LISA.cocaine <- read.csv("CountyKNN3.csv")
LISA_C_names <- names(LISA.cocaine)[grep("LISA_C", names(LISA.cocaine))]
names(LISA.rel)[grep("LISA_C", names(LISA.rel))] <- LISA_C_names
names(LISA.abs.1)[grep("LISA_C", names(LISA.rel))] <- LISA_C_names
names(LISA.abs.t_1)[grep("LISA_C", names(LISA.rel))] <- LISA_C_names
LISA_C.rel <- LISA.rel[,grep("LISA_C", names(LISA.rel))]
LISA_C.abs.1 <- LISA.abs.1[,grep("LISA_C", names(LISA.rel))]
LISA_C.abs.t_1 <- LISA.abs.t_1[,grep("LISA_C", names(LISA.rel))]
county.names <- LISA.rel[, c(2,6,8)]

LISA_C.rel <- cbind(county.names, LISA_C.rel)
LISA_C.abs.1 <- cbind(county.names, LISA_C.abs.1)
LISA_C.abs.t_1 <- cbind(county.names, LISA_C.abs.t_1)

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
counties <- counties %>% filter(!(state_name %in% c("Alaska", "Hawaii")))

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
  ggsave(paste("all drugs seizure_weights_", year, ".jpg", sep=""), seizure_weights_map, width=55, height=30, units="cm")
  
  LISA.rel.map %>% filter(grepl(year, Month_Year)) %>% 
    ggplot(mapping = aes(long, lat, group = GEOID, fill=log(weights_per_pop.))) +
    geom_polygon(color = "#000000", size = .05) +
    facet_wrap(. ~ Month_Year) +
    coord_equal() +
    scale_fill_viridis_c() +
    labs(fill = "log Weights Per Pop.") -> weights_per_pop.map
  ggsave(paste("all drugs log weights_per_pop_", year, ".jpg", sep=""), weights_per_pop.map, width=55, height=30, units="cm")
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
  ggsave(paste("All_Drugs_LISA_C_rel_", year, "_XCal.jpg", sep=""), LISA_C_map, width=40, height=30, units="cm")
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
  ggsave(paste("All_Drugs_LISA_C_abs1_", year, "_XCal.jpg", sep=""), LISA_C_map, width=40, height=30, units="cm")
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
  ggsave(paste("All_Drugs_LISA_C_abs_t1_", year, "_XCal.jpg", sep=""), LISA_C_map, width=40, height=30, units="cm")
}

#######
## Exclude Top 4 States (California, Oregon, Kentucky, Texas)
LISA3.XCal <- LISA3 %>% filter(state %in% c("Tennessee", "Louisiana", "Pennsylvania", "Alabama","New York", "Illinois", "New Hampshire",
                                            "Ohio", "Georgia", "North Carolina", "Michigan", "New Jersey", "Wisconsin", "Virginia",
                                            "South Carolina", "Mississippi", "Massachusetts", "Rhode Island", "Indiana",
                                            "Connecticut", "Delaware", "Arkansas", "West Virginia", "Florida"))
coords.XCal <- counties.obs %>% filter(GEOID %in% LISA3.XCal$GEOID) %>%
  group_by(GEOID) %>% summarise(x=mean(long), y=mean(lat))
GEOIDS.XCal <- coords.XCal$GEOID
coords.XCal <- coords.XCal[,-1]
LISA3_nb_XCal <- knn2nb(knearneigh(coords.XCal, k=3), row.names=GEOIDS.XCal)

alpha <- 0.05
nperm <- 999
LISA_I.index <- grep("LISA_I", names(LISA3.XCal))[1]-1
nb.obj.XCal <- nb2listw(LISA3_nb_XCal, style="B")
nrow(LISA3.XCal) # 1386, 1422-1386=36 counties in Cal removed
seizures.XCal <- LISA3.XCal[, 21:68]

grep("LISA_I", names(LISA3.XCal))[1] # 69

relative.MoransI <- LISA3.XCal
set.seed(100)
for (i in 1:48) {
  seizure.XCal <- t(seizures.XCal)[i,]
  localM.month <- localmoran_abs(seizure.XCal, nb.obj.XCal, nsim=nperm, zero.policy=T, xx=NULL)
  localM.month$LISA_C <- ifelse(localM.month$mean=="High-High", 1,
                                ifelse(localM.month$mean=="Low-Low", 2,
                                       ifelse(localM.month$mean=="Low-High", 3, 4)))
  localM.month$LISA_C <- ifelse(localM.month$`Pr(folded) Sim` <= alpha, localM.month$LISA_C, 0)
  relative.MoransI[,(69+3*i-3):(69+3*i-1)] <- localM.month[,c(1,13,7)]
}

# write.csv(relative.MoransI, "CountyKNN3 East_NoKentucky R codes 999 (relative).csv", row.names=F)

absolute.b1.MoransI <- LISA3.XCal
x.bar_p <- mean(absolute.b1.MoransI$Jan_2018)
set.seed(300)
for (i in 1:48) {
  seizure.XCal <- t(seizures.XCal)[i,]
  localM.month <- localmoran_abs(seizure.XCal, nb.obj.XCal, nsim=nperm, zero.policy=T, xx=NULL)
  localM.month$LISA_C <- ifelse(localM.month$mean=="High-High", 1,
                                ifelse(localM.month$mean=="Low-Low", 2,
                                       ifelse(localM.month$mean=="Low-High", 3, 4)))
  localM.month$LISA_C <- ifelse(localM.month$`Pr(folded) Sim` <= alpha, localM.month$LISA_C, 0)
  absolute.b1.MoransI[,(69+3*i-3):(69+3*i-1)] <- localM.month[,c(1,13,7)]
}

# write.csv(absolute.b1.MoransI, "CountyKNN3 East_NoKentucky R codes 999 (absolute base 1).csv", row.names=F)

absolute.t1.MoransI <- LISA3.XCal
x.bar_p <- mean(absolute.t1.MoransI$Jan_2018)
set.seed(500)
for (i in 1:48) {
  seizure.XCal <- t(seizures.XCal)[i,]
  localM.month <- localmoran_abs(seizure.XCal, nb.obj.XCal, nsim=nperm, zero.policy=T, xx=NULL)
  localM.month$LISA_C <- ifelse(localM.month$mean=="High-High", 1,
                                ifelse(localM.month$mean=="Low-Low", 2,
                                       ifelse(localM.month$mean=="Low-High", 3, 4)))
  localM.month$LISA_C <- ifelse(localM.month$`Pr(folded) Sim` <= alpha, localM.month$LISA_C, 0)
  absolute.t1.MoransI[,(69+3*i-3):(69+3*i-1)] <- localM.month[,c(1,13,7)]
  x.bar_p <- mean(as.data.frame(LISA3.XCal)[,21+i-1])
}

# write.csv(absolute.t1.MoransI, "CountyKNN3 East_NoKentucky R codes 999 (absolute base t_1).csv", row.names=F)

##########
## Monthly Plots
LISA.rel <- read.csv("CountyKNN3 East_NoKentucky R codes 999 (relative).csv") %>% as_tibble
LISA.abs.1 <- read.csv("CountyKNN3 East_NoKentucky R codes 999 (absolute base 1).csv") %>% as_tibble
LISA.abs.t_1 <- read.csv("CountyKNN3 East_NoKentucky R codes 999 (absolute base t_1).csv") %>% as_tibble

LISA.cocaine <- read.csv("CountyKNN3.csv")
LISA_C_names <- names(LISA.cocaine)[grep("LISA_C", names(LISA.cocaine))]
names(LISA.rel)[grep("LISA_C", names(LISA.rel))] <- LISA_C_names
names(LISA.abs.1)[grep("LISA_C", names(LISA.rel))] <- LISA_C_names
names(LISA.abs.t_1)[grep("LISA_C", names(LISA.rel))] <- LISA_C_names
LISA_C.rel <- LISA.rel[,grep("LISA_C", names(LISA.rel))]
LISA_C.abs.1 <- LISA.abs.1[,grep("LISA_C", names(LISA.rel))]
LISA_C.abs.t_1 <- LISA.abs.t_1[,grep("LISA_C", names(LISA.rel))]
county.names <- LISA.rel[, c(2,6,8)]

LISA_C.rel <- cbind(county.names, LISA_C.rel)
LISA_C.abs.1 <- cbind(county.names, LISA_C.abs.1)
LISA_C.abs.t_1 <- cbind(county.names, LISA_C.abs.t_1)

counties.original <- counties #<- counties.original
names(counties)[7] <- "GEOID"
counties$GEOID <- as.numeric(counties$GEOID)
counties <- counties %>% filter(!(state_name %in% c("Alaska", "Hawaii")))


LISA_C.rel.map <- left_join(counties[, c(1:2,7)], LISA_C.rel, by = "GEOID")
LISA_C.rel.map <- LISA_C.rel.map %>% pivot_longer(c(-long, -lat, -GEOID, -state, -NAMELSAD), names_to="Month_Year", values_to="LISA_C")
ref <- tibble(LISA_C=0:4, label=c("Insig", "HH", "LL", "LH", "HL"))
LISA_C.rel.map$LISA_C <- left_join(LISA_C.rel.map[, 6:7], ref, by="LISA_C")$label
LISA_C.rel.map$LISA_C <- factor(LISA_C.rel.map$LISA_C)
LISA_C.rel.map$Month_Year <- factor(LISA_C.rel.map$Month_Year, levels=unique(LISA_C.rel.map$Month_Year))

for (year in 2018:2021) { # monthly maps of LISA_C.rel
  LISA_C.rel.map %>% filter(substr(Month_Year,10,10) == substr(year,4,4)) %>% 
    ggplot(mapping = aes(long, lat, group = GEOID, fill=LISA_C)) +
    xlim(c(-97,-72)) + ylim(c(22, 50)) +
    geom_polygon(color = "#000000", size = .05) +
    facet_wrap(. ~ Month_Year) +
    scale_fill_manual(values = c("Insig"="grey90",
                                 "LL"="blue",
                                 "LH"="steelblue",
                                 "HL"="orange",
                                 "HH"="red")) +
    labs(fill = "LISA_C_rel") -> LISA_C_map
  ggsave(paste("All_Drugs_LISA_C_rel_", year, "_East.jpg", sep=""), LISA_C_map, width=40, height=30, units="cm")
}

LISA_C.abs1.map <- left_join(counties[, c(1:2,7)], LISA_C.abs.1, by = "GEOID")
LISA_C.abs1.map <- LISA_C.abs1.map %>% pivot_longer(c(-long, -lat, -GEOID, -state, -NAMELSAD), names_to="Month_Year", values_to="LISA_C")
LISA_C.abs1.map$LISA_C <- left_join(LISA_C.abs1.map[, 6:7], ref, by="LISA_C")$label
LISA_C.abs1.map$LISA_C <- factor(LISA_C.abs1.map$LISA_C)
LISA_C.abs1.map$Month_Year <- factor(LISA_C.abs1.map$Month_Year, levels=unique(LISA_C.abs1.map$Month_Year))

for (year in 2018:2021) { # monthly maps of LISA_C.abs1
  LISA_C.abs1.map %>% filter(substr(Month_Year,10,10) == substr(year,4,4)) %>% 
    ggplot(mapping = aes(long, lat, group = GEOID, fill=LISA_C)) +
    xlim(c(-97,-72)) + ylim(c(22, 50)) +
    geom_polygon(color = "#000000", size = .05) +
    facet_wrap(. ~ Month_Year) +
    scale_fill_manual(values = c("Insig"="grey90",
                                 "LL"="blue",
                                 "LH"="steelblue",
                                 "HL"="orange",
                                 "HH"="red")) +
    labs(fill = "LISA_C_abs1") -> LISA_C_map
  ggsave(paste("All_Drugs_LISA_C_abs1_", year, "_East.jpg", sep=""), LISA_C_map, width=40, height=30, units="cm")
}

LISA_C.abs_t1.map <- left_join(counties[, c(1:2,7)], LISA_C.abs.t_1, by = "GEOID")
LISA_C.abs_t1.map <- LISA_C.abs_t1.map %>% pivot_longer(c(-long, -lat, -GEOID, -state, -NAMELSAD), names_to="Month_Year", values_to="LISA_C")
LISA_C.abs_t1.map$LISA_C <- left_join(LISA_C.abs_t1.map[, 6:7], ref, by="LISA_C")$label
LISA_C.abs_t1.map$LISA_C <- factor(LISA_C.abs_t1.map$LISA_C)
LISA_C.abs_t1.map$Month_Year <- factor(LISA_C.abs_t1.map$Month_Year, levels=unique(LISA_C.abs_t1.map$Month_Year))

for (year in 2018:2021) { # monthly maps of LISA_C.abs_t1
  LISA_C.abs_t1.map %>% filter(substr(Month_Year,10,10) == substr(year,4,4)) %>% 
    ggplot(mapping = aes(long, lat, group = GEOID, fill=LISA_C)) +
    xlim(c(-97,-72)) + ylim(c(22, 50)) +
    geom_polygon(color = "#000000", size = .05) +
    facet_wrap(. ~ Month_Year) +
    scale_fill_manual(values = c("Insig"="grey90",
                                 "LL"="blue",
                                 "LH"="steelblue",
                                 "HL"="orange",
                                 "HH"="red")) +
    labs(fill = "LISA_C_abs_t1") -> LISA_C_map
  ggsave(paste("All_Drugs_LISA_C_abs_t1_", year, "_East.jpg", sep=""), LISA_C_map, width=40, height=30, units="cm")
}


## Explanation for average seizure weights per state
avg.seizure <- read.csv("All Drugs Average Seizures.csv")
state.centroid <- counties %>% group_by(state_name) %>% summarise(long=mean(long), lat=mean(lat))
names(state.centroid)[1] <- "state"
avg.seizure <- left_join(avg.seizure, state.centroid, by="state")
n_state <- nrow(avg.seizure)
avg.seizure$dist_from_FL <- numeric(n_state)

for (i in 1:n_state) {
  if (avg.seizure$state[i] == "Florida") next
  avg.seizure$dist_from_FL[i] <- sqrt(sum((avg.seizure[i,3:4] - avg.seizure[7,3:4])^2))
}

avg.seizure$dist_from_SB <- numeric(n_state)
Southern.Border.centroid <- avg.seizure[c(1,4,6,18),3:4] %>% apply(2, mean)
for (i in 1:n_state) {
  if (avg.seizure$state[i] %in% c("Arizona", "California", "New Mexico", "Texas")) next
  avg.seizure$dist_from_SB[i] <- sqrt(sum((avg.seizure[i,3:4] - Southern.Border.centroid)^2))
}

avg.seizure
state.populations <- populations %>% group_by(state) %>% summarise(population=sum(POPESTIMATE2018))
interstate <- read.csv("interstate.csv") %>% as_tibble
interstate$interstate <- interstate$RURAL.INTERSTATE + interstate$URBAN.INTERSTATE
names(interstate)[1] <- "state"
avg.seizure <- left_join(avg.seizure, state.populations, by="state")
avg.seizure <- left_join(avg.seizure, interstate[,c(1,17)], by="state")

pairs(avg.seizure[,-c(1,3:4)])
lm(avg.seizure~population+interstate+dist_from_SB+dist_from_FL, data=avg.seizure) %>% summary
lm(avg.seizure~population+dist_from_SB+dist_from_FL, data=avg.seizure) %>% summary
lm(avg.seizure~interstate+dist_from_SB+dist_from_FL, data=avg.seizure) %>% summary
