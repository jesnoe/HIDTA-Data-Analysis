library(fpp2)
library(spdep)
library(readxl)
library(maps)
library(urbnmapr)
library(tidyverse)
library(gridExtra)
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


# 1  District of Columbia Washington/Baltimore District of Columbia County
# 2             Louisiana           Gulf Coast              Bossier County -> Parish
# 3             Louisiana           Gulf Coast            Calcasieu County
# 4             Louisiana           Gulf Coast     East Baton Rouge County
# 5              Maryland Washington/Baltimore       Baltimore City County -> just city
# 6             Minnesota        North Central             Olmstead County -> Olmsted
# 7             Minnesota        North Central          Saint Louis County -> St. Louis
# 8              Missouri              Midwest       St. Louis City County -> just city
# 9              Virginia           Appalachia                Galax County -> city
# 10             Virginia Washington/Baltimore           Alexandria County -> city
# 11             Virginia Washington/Baltimore           Chesapeake County -> city
# 12             Virginia Washington/Baltimore         Fairfax City County
# 13             Virginia Washington/Baltimore         Falls Church County -> city
# 14             Virginia Washington/Baltimore             Hopewell County -> city
# 15             Virginia Washington/Baltimore             Manassas County -> city
# 16             Virginia Washington/Baltimore        Manassas Park County -> city
# 17             Virginia Washington/Baltimore         Newport News County -> city
# 18             Virginia Washington/Baltimore           Petersburg County -> city
# 19             Virginia Washington/Baltimore           Portsmouth County -> city
# 20             Virginia Washington/Baltimore        Richmond City County
# 21             Virginia Washington/Baltimore       Virginia Beach County -> city

HIDTA.dist <- read.csv("HIDTA Regions.csv", header=T) %>% as_tibble
# HIDTA.dist$county_name <- ifelse(HIDTA.dist$state_name == "Louisiana", paste(HIDTA.dist$county_name, "Parish"),
#                                  ifelse(grepl("city", HIDTA.dist$county_name), HIDTA.dist$county_name, paste(HIDTA.dist$county_name, "County")))
# HIDTA.dist %>% filter(!(state_name %in% c("Alaska", "Hawaii", "Puerto Rico")), !(county_name %in% counties$county_name)) %>% as.data.frame
# write.csv(HIDTA.dist, "HIDTA Regions.csv")
coordinate.HIDTA <- left_join(counties, HIDTA.dist, by=c("state_name", "county_name")) %>% filter(!(state_name %in% c("Alaska", "Hawaii", "Puerto Rico")))
names(coordinate.HIDTA)[7] <- "GEOID"


centroids.HIDTA <- coordinate.HIDTA %>% filter(!(state_name %in% c("Alaska", "Hawaii")), !is.na(HIDTA)) %>%
  group_by(HIDTA) %>% summarise(long=mean(long), lat=mean(lat))

coordinate.HIDTA %>% ggplot() +
  geom_polygon(aes(long, lat, group=GEOID, fill=HIDTA),
               color = "#000000",
               size = .05)  +
  geom_point(aes(x=long, y=lat, fill=HIDTA),
             size=3, pch=21,
             data=centroids.HIDTA) + 
  geom_text(aes(long, lat, label=HIDTA), 
            data=centroids.HIDTA,
            vjust=-0.5,
            fontface="bold") +
  xlim(-132, -68) -> HIDTA.map.w.NA

coordinate.HIDTA %>% filter(coordinate.HIDTA$GEOID %in% unique(LISA3$GEOID)) %>% ggplot() +
  geom_polygon(aes(long, lat, group=GEOID, fill=HIDTA),
               color = "#000000",
               size = .05) +
  geom_point(aes(x=long, y=lat, fill=HIDTA),
             size=3, pch=21,
             data=centroids.HIDTA) + 
  geom_text(aes(long, lat, label=HIDTA), 
            data=centroids.HIDTA,
            vjust=-0.5,
            fontface="bold") +
  xlim(-132, -68) -> HIDTA.map.wo.NA

HIDTA.regions <- grid.arrange(HIDTA.map.w.NA, HIDTA.map.wo.NA, nrow=2)
ggsave("HIDTA Regions.jpg", HIDTA.regions, width=28, height=25, units="cm")



HIDTA.regions <- coordinate.HIDTA %>% filter(!is.na(HIDTA)) %>% select(GEOID, HIDTA) %>% unique
LISA3$HIDTA <- left_join(LISA3, HIDTA.regions, by="GEOID")$HIDTA

seizures.HIDTA <- (LISA3 %>% group_by(HIDTA))[,c(213,21:68)]
seizures.HIDTA$HIDTA[is.na(seizures.HIDTA$HIDTA)] <- "NA"
monthly.seizures <- data.frame()
for (HIDTA.region in unique(seizures.HIDTA$HIDTA)) {
  seizures.i <- seizures.HIDTA %>% filter(HIDTA == HIDTA.region)
  monthly.seizures <- rbind(monthly.seizures, c(HIDTA.region, apply(seizures.i[,-1], 2, sum)))
}
names(monthly.seizures) <- names(seizures.HIDTA)
monthly.seizures[,-1] <- apply(monthly.seizures[,-1], 2, as.numeric)
monthly.seizures$avg.seizure <- apply(monthly.seizures[,-1], 1, mean)
avg.seizures <- monthly.seizures %>% select(HIDTA, avg.seizure) %>% arrange(desc(avg.seizure))
#write.csv(avg.seizures, "All Drugs Average Seizures per HIDTA.csv", row.names=F)



coords.US <- coordinate.HIDTA %>% filter(!is.na(HIDTA)) %>% group_by(HIDTA) %>% summarise(x=mean(long), y=mean(lat))
HIDTA.regions <- coords.US$HIDTA
coords.US <- coords.US[,-1]
LISA3_nb <- knn2nb(knearneigh(coords.US, k=3), row.names=HIDTA.regions)
LISA5_nb <- knn2nb(knearneigh(coords.US, k=5), row.names=HIDTA.regions)


## All HIDTA
LISA3.All_HIDTA <- monthly.seizures[,-50] %>% filter(HIDTA != "NA")
nb.obj.All_HIDTA <- nb2listw(LISA3_nb, style="B")
alpha <- 0.05
nperm <- 999
seizures.All_HIDTA <- LISA3.All_HIDTA[, -1]

relative.MoransI <- LISA3.All_HIDTA
set.seed(100)
for (i in 1:48) {
  seizure.All_HIDTA <- t(seizures.All_HIDTA)[i,]
  localM.month <- localmoran_abs(seizure.All_HIDTA, nb.obj.All_HIDTA, nsim=nperm, zero.policy=T, xx=NULL)
  localM.month$LISA_C <- ifelse(localM.month$mean=="High-High", 1,
                                ifelse(localM.month$mean=="Low-Low", 2,
                                       ifelse(localM.month$mean=="Low-High", 3, 4)))
  localM.month$LISA_C <- ifelse(localM.month$`Pr(folded) Sim` <= alpha, localM.month$LISA_C, 0)
  relative.MoransI <- cbind(relative.MoransI, localM.month[,c(1,13,7)])
}
LISA.cocaine <- read.csv("CountyKNN3.csv")
LISA_colnames <- names(LISA.cocaine)[grep("LISA", names(LISA.cocaine))]
names(relative.MoransI)[50:193] <- LISA_colnames
# write.csv(relative.MoransI, "CountyKNN3 All_HIDTA R codes 999 (relative).csv", row.names=F)

absolute.b1.MoransI <- LISA3.All_HIDTA
x.bar_p <- mean(absolute.b1.MoransI$Jan_2018)
set.seed(300)
for (i in 1:48) {
  seizure.All_HIDTA <- t(seizures.All_HIDTA)[i,]
  localM.month <- localmoran_abs(seizure.All_HIDTA, nb.obj.All_HIDTA, nsim=nperm, zero.policy=T, xx=x.bar_p)
  localM.month$LISA_C <- ifelse(localM.month$mean=="High-High", 1,
                                ifelse(localM.month$mean=="Low-Low", 2,
                                       ifelse(localM.month$mean=="Low-High", 3, 4)))
  localM.month$LISA_C <- ifelse(localM.month$`Pr(folded) Sim` <= alpha, localM.month$LISA_C, 0)
  absolute.b1.MoransI <- cbind(absolute.b1.MoransI, localM.month[,c(1,13,7)])
}
names(absolute.b1.MoransI)[50:193] <- LISA_colnames
# write.csv(absolute.b1.MoransI, "CountyKNN3 All_HIDTA R codes 999 (absolute base 1).csv", row.names=F)

absolute.t1.MoransI <- LISA3.All_HIDTA
x.bar_p <- mean(absolute.t1.MoransI$Jan_2018)
set.seed(500)
for (i in 1:48) {
  seizure.All_HIDTA <- t(seizures.All_HIDTA)[i,]
  localM.month <- localmoran_abs(seizure.All_HIDTA, nb.obj.All_HIDTA, nsim=nperm, zero.policy=T, xx=x.bar_p)
  localM.month$LISA_C <- ifelse(localM.month$mean=="High-High", 1,
                                ifelse(localM.month$mean=="Low-Low", 2,
                                       ifelse(localM.month$mean=="Low-High", 3, 4)))
  localM.month$LISA_C <- ifelse(localM.month$`Pr(folded) Sim` <= alpha, localM.month$LISA_C, 0)
  absolute.t1.MoransI <- cbind(absolute.t1.MoransI, localM.month[,c(1,13,7)])
  x.bar_p <- mean(as.data.frame(LISA3.All_HIDTA)[,2+i-1])
}
names(absolute.t1.MoransI)[50:193] <- LISA_colnames
# write.csv(absolute.t1.MoransI, "CountyKNN3 All_HIDTA R codes 999 (absolute base t_1).csv", row.names=F)


# All_HIDTA Plots
LISA.rel <- read.csv("CountyKNN3 All_HIDTA R codes 999 (relative).csv") %>% as_tibble
LISA.abs.1 <- read.csv("CountyKNN3 All_HIDTA R codes 999 (absolute base 1).csv") %>% as_tibble
LISA.abs.t_1 <- read.csv("CountyKNN3 All_HIDTA R codes 999 (absolute base t_1).csv") %>% as_tibble

LISA.rel[,grep("LISA_I", names(LISA.rel))]; LISA.abs.1[,grep("LISA_I", names(LISA.rel))]; LISA.abs.t_1[,grep("LISA_I", names(LISA.rel))]

LISA_C.rel <- LISA.rel[,grep("LISA_C", names(LISA.rel))]
LISA_C.abs.1 <- LISA.abs.1[,grep("LISA_C", names(LISA.rel))]
LISA_C.abs.t_1 <- LISA.abs.t_1[,grep("LISA_C", names(LISA.rel))]
HIDTA.names <- LISA.rel[, 1]

LISA_C.rel <- cbind(HIDTA.names, LISA_C.rel)
LISA_C.abs.1 <- cbind(HIDTA.names, LISA_C.abs.1)
LISA_C.abs.t_1 <- cbind(HIDTA.names, LISA_C.abs.t_1)

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
  HIDTA.pop <- left_join(county.fips, coordinate.HIDTA %>% filter(!is.na(HIDTA)) %>% select(HIDTA, GEOID) %>% unique) %>% 
    group_by(HIDTA) %>% summarise(POPESTIMATE2018=sum(POPESTIMATE2018),
                                  POPESTIMATE2019=sum(POPESTIMATE2019),
                                  POPESTIMATE2020=sum(POPESTIMATE2020),
                                  POPESTIMATE2021=sum(POPESTIMATE2021))
}

HIDTAs <- coordinate.HIDTA

LISA_C.rel.map <- left_join(HIDTAs[, c(1:2,7,13)], LISA_C.rel, by = "HIDTA")
LISA_C.rel.map <- LISA_C.rel.map %>% pivot_longer(c(-long, -lat, -GEOID, -HIDTA), names_to="Month_Year", values_to="LISA_C")
ref <- tibble(LISA_C=0:4, label=c("Insig", "HH", "LL", "LH", "HL"))
LISA_C.rel.map$LISA_C <- left_join(LISA_C.rel.map[, 5:6], ref, by="LISA_C")$label
LISA_C.rel.map$LISA_C <- factor(LISA_C.rel.map$LISA_C)
LISA_C.rel.map$Month_Year <- factor(LISA_C.rel.map$Month_Year, levels=unique(LISA_C.rel.map$Month_Year))

for (year in 2018:2021) { # monthly maps of LISA_C.rel
  LISA_C.rel.map %>% filter(substr(Month_Year,10,10) == substr(year,4,4)) %>% 
    ggplot() + geom_polygon(mapping = aes(long, lat, group = GEOID, fill=LISA_C),
                            color = "#000000", size = .01)   +
    geom_point(aes(x=long, y=lat),
               size=2, pch=19, color="green",
               data=centroids.HIDTA) +
    xlim(-132, -68) +
    facet_wrap(. ~ Month_Year) +
    scale_fill_manual(values = c("Insig"="grey90",
                                 "LL"="blue",
                                 "LH"="steelblue",
                                 "HL"="orange",
                                 "HH"="red")) +
    labs(fill = "LISA_C_rel") -> LISA_C_map
  ggsave(paste("All_Drugs_LISA_C_rel_", year, "_All_HIDTA.jpg", sep=""), LISA_C_map, width=55, height=30, units="cm")
}

LISA_C.abs1.map <- left_join(HIDTAs[, c(1:2,7,13)], LISA_C.abs.1, by = "HIDTA")
LISA_C.abs1.map <- LISA_C.abs1.map %>% pivot_longer(c(-long, -lat, -GEOID, -HIDTA), names_to="Month_Year", values_to="LISA_C")
LISA_C.abs1.map$LISA_C <- left_join(LISA_C.abs1.map[, 5:6], ref, by="LISA_C")$label
LISA_C.abs1.map$LISA_C <- factor(LISA_C.abs1.map$LISA_C)
LISA_C.abs1.map$Month_Year <- factor(LISA_C.abs1.map$Month_Year, levels=unique(LISA_C.abs1.map$Month_Year))

for (year in 2018:2021) { # monthly maps of LISA_C.abs1
  LISA_C.abs1.map %>% filter(substr(Month_Year,10,10) == substr(year,4,4)) %>% 
    ggplot() + geom_polygon(mapping = aes(long, lat, group = GEOID, fill=LISA_C),
                            color = "#000000", size = .01)   +
    geom_point(aes(x=long, y=lat),
               size=2, pch=19, color="green",
               data=centroids.HIDTA) +
    xlim(-132, -68) +
    facet_wrap(. ~ Month_Year) +
    scale_fill_manual(values = c("Insig"="grey90",
                                 "LL"="blue",
                                 "LH"="steelblue",
                                 "HL"="orange",
                                 "HH"="red")) +
    labs(fill = "LISA_C_abs1") -> LISA_C_map
  ggsave(paste("All_Drugs_LISA_C_abs1_", year, "_All_HIDTA.jpg", sep=""), LISA_C_map, width=55, height=30, units="cm")
}

LISA_C.abs_t1.map <- left_join(HIDTAs[, c(1:2,7,13)], LISA_C.abs.t_1, by = "HIDTA")
LISA_C.abs_t1.map <- LISA_C.abs_t1.map %>% pivot_longer(c(-long, -lat, -GEOID, -HIDTA), names_to="Month_Year", values_to="LISA_C")
LISA_C.abs_t1.map$LISA_C <- left_join(LISA_C.abs_t1.map[, 5:6], ref, by="LISA_C")$label
LISA_C.abs_t1.map$LISA_C <- factor(LISA_C.abs_t1.map$LISA_C)
LISA_C.abs_t1.map$Month_Year <- factor(LISA_C.abs_t1.map$Month_Year, levels=unique(LISA_C.abs_t1.map$Month_Year))

for (year in 2018:2021) { # monthly maps of LISA_C.abs_t1
  LISA_C.abs_t1.map %>% filter(substr(Month_Year,10,10) == substr(year,4,4)) %>% 
    ggplot() + geom_polygon(mapping = aes(long, lat, group = GEOID, fill=LISA_C),
                            color = "#000000", size = .01)   +
    geom_point(aes(x=long, y=lat),
               size=2, pch=19, color="green",
               data=centroids.HIDTA) +
    xlim(-132, -68) +
    facet_wrap(. ~ Month_Year) +
    scale_fill_manual(values = c("Insig"="grey90",
                                 "LL"="blue",
                                 "LH"="steelblue",
                                 "HL"="orange",
                                 "HH"="red")) +
    labs(fill = "LISA_C_abs_t1") -> LISA_C_map
  ggsave(paste("All_Drugs_LISA_C_abs_t1_", year, "_All_HIDTA.jpg", sep=""), LISA_C_map, width=55, height=30, units="cm")
}


## No California
LISA3.XCal <- monthly.seizures[,-50] %>% filter(!(HIDTA %in% c("NA", "Central Valley California", "Los Angeles")))
coords.XCal <- coordinate.HIDTA %>% filter(!is.na(HIDTA),!(HIDTA %in% c("Central Valley California", "Los Angeles"))) %>%
  group_by(HIDTA) %>% summarise(x=mean(long), y=mean(lat))
HIDTA.XCal <- coords.XCal$HIDTA
coords.XCal <- coords.XCal[,-1]
LISA3_nb_XCal <- knn2nb(knearneigh(coords.XCal, k=3), row.names=HIDTA.XCal)
nb.obj.XCal <- nb2listw(LISA3_nb_XCal, style="B")
alpha <- 0.05
nperm <- 999
seizures.XCal <- LISA3.XCal[, -1]

relative.MoransI <- LISA3.XCal
set.seed(100)
for (i in 1:48) {
  seizure.XCal <- t(seizures.XCal)[i,]
  localM.month <- localmoran_abs(seizure.XCal, nb.obj.XCal, nsim=nperm, zero.policy=T, xx=NULL)
  localM.month$LISA_C <- ifelse(localM.month$mean=="High-High", 1,
                                ifelse(localM.month$mean=="Low-Low", 2,
                                       ifelse(localM.month$mean=="Low-High", 3, 4)))
  localM.month$LISA_C <- ifelse(localM.month$`Pr(folded) Sim` <= alpha, localM.month$LISA_C, 0)
  relative.MoransI <- cbind(relative.MoransI, localM.month[,c(1,13,7)])
}
LISA_colnames <- names(LISA.cocaine)[grep("LISA", names(LISA.cocaine))]
names(relative.MoransI)[50:193] <- LISA_colnames
# write.csv(relative.MoransI, "CountyKNN3 XCal_HIDTA R codes 999 (relative).csv", row.names=F)

absolute.b1.MoransI <- LISA3.XCal
x.bar_p <- mean(absolute.b1.MoransI$Jan_2018)
set.seed(300)
for (i in 1:48) {
  seizure.XCal <- t(seizures.XCal)[i,]
  localM.month <- localmoran_abs(seizure.XCal, nb.obj.XCal, nsim=nperm, zero.policy=T, xx=x.bar_p)
  localM.month$LISA_C <- ifelse(localM.month$mean=="High-High", 1,
                                ifelse(localM.month$mean=="Low-Low", 2,
                                       ifelse(localM.month$mean=="Low-High", 3, 4)))
  localM.month$LISA_C <- ifelse(localM.month$`Pr(folded) Sim` <= alpha, localM.month$LISA_C, 0)
  absolute.b1.MoransI <- cbind(absolute.b1.MoransI, localM.month[,c(1,13,7)])
}
names(absolute.b1.MoransI)[50:193] <- LISA_colnames
# write.csv(absolute.b1.MoransI, "CountyKNN3 XCal_HIDTA R codes 999 (absolute base 1).csv", row.names=F)

absolute.t1.MoransI <- LISA3.XCal
x.bar_p <- mean(absolute.t1.MoransI$Jan_2018)
set.seed(500)
for (i in 1:48) {
  seizure.XCal <- t(seizures.XCal)[i,]
  localM.month <- localmoran_abs(seizure.XCal, nb.obj.XCal, nsim=nperm, zero.policy=T, xx=x.bar_p)
  localM.month$LISA_C <- ifelse(localM.month$mean=="High-High", 1,
                                ifelse(localM.month$mean=="Low-Low", 2,
                                       ifelse(localM.month$mean=="Low-High", 3, 4)))
  localM.month$LISA_C <- ifelse(localM.month$`Pr(folded) Sim` <= alpha, localM.month$LISA_C, 0)
  absolute.t1.MoransI <- cbind(absolute.t1.MoransI, localM.month[,c(1,13,7)])
  x.bar_p <- mean(as.data.frame(LISA3.XCal)[,2+i-1])
}
names(absolute.t1.MoransI)[50:193] <- LISA_colnames
# write.csv(absolute.t1.MoransI, "CountyKNN3 XCal_HIDTA R codes 999 (absolute base t_1).csv", row.names=F)


# XCal Plots
LISA.rel <- read.csv("CountyKNN3 XCal_HIDTA R codes 999 (relative).csv") %>% as_tibble
LISA.abs.1 <- read.csv("CountyKNN3 XCal_HIDTA R codes 999 (absolute base 1).csv") %>% as_tibble
LISA.abs.t_1 <- read.csv("CountyKNN3 XCal_HIDTA R codes 999 (absolute base t_1).csv") %>% as_tibble

LISA.rel[,grep("LISA_I", names(LISA.rel))]; LISA.abs.1[,grep("LISA_I", names(LISA.rel))]; LISA.abs.t_1[,grep("LISA_I", names(LISA.rel))]

LISA_C.rel <- LISA.rel[,grep("LISA_C", names(LISA.rel))]
LISA_C.abs.1 <- LISA.abs.1[,grep("LISA_C", names(LISA.rel))]
LISA_C.abs.t_1 <- LISA.abs.t_1[,grep("LISA_C", names(LISA.rel))]
HIDTA.names <- LISA.rel[, 1]

LISA_C.rel <- cbind(HIDTA.names, LISA_C.rel)
LISA_C.abs.1 <- cbind(HIDTA.names, LISA_C.abs.1)
LISA_C.abs.t_1 <- cbind(HIDTA.names, LISA_C.abs.t_1)

HIDTAs <- coordinate.HIDTA

LISA_C.rel.map <- left_join(HIDTAs[, c(1:2,7,13)], LISA_C.rel, by = "HIDTA")
LISA_C.rel.map <- LISA_C.rel.map %>% pivot_longer(c(-long, -lat, -GEOID, -HIDTA), names_to="Month_Year", values_to="LISA_C")
ref <- tibble(LISA_C=0:4, label=c("Insig", "HH", "LL", "LH", "HL"))
LISA_C.rel.map$LISA_C <- left_join(LISA_C.rel.map[, 5:6], ref, by="LISA_C")$label
LISA_C.rel.map$LISA_C <- factor(LISA_C.rel.map$LISA_C)
LISA_C.rel.map$Month_Year <- factor(LISA_C.rel.map$Month_Year, levels=unique(LISA_C.rel.map$Month_Year))

for (year in 2018:2021) { # monthly maps of LISA_C.rel
  LISA_C.rel.map %>% filter(substr(Month_Year,10,10) == substr(year,4,4)) %>% 
    ggplot() + geom_polygon(mapping = aes(long, lat, group = GEOID, fill=LISA_C),
                            color = "#000000", size = .01)   +
    geom_point(aes(x=long, y=lat),
               size=2, pch=19, color="green",
               data=centroids.HIDTA) +
    xlim(-132, -68) +
    facet_wrap(. ~ Month_Year) +
    scale_fill_manual(values = c("Insig"="grey90",
                                 "LL"="blue",
                                 "LH"="steelblue",
                                 "HL"="orange",
                                 "HH"="red")) +
    labs(fill = "LISA_C_rel") -> LISA_C_map
  ggsave(paste("All_Drugs_LISA_C_rel_", year, "_XCal_HIDTA.jpg", sep=""), LISA_C_map, width=55, height=30, units="cm")
}

LISA_C.abs1.map <- left_join(HIDTAs[, c(1:2,7,13)], LISA_C.abs.1, by = "HIDTA")
LISA_C.abs1.map <- LISA_C.abs1.map %>% pivot_longer(c(-long, -lat, -GEOID, -HIDTA), names_to="Month_Year", values_to="LISA_C")
LISA_C.abs1.map$LISA_C <- left_join(LISA_C.abs1.map[, 5:6], ref, by="LISA_C")$label
LISA_C.abs1.map$LISA_C <- factor(LISA_C.abs1.map$LISA_C)
LISA_C.abs1.map$Month_Year <- factor(LISA_C.abs1.map$Month_Year, levels=unique(LISA_C.abs1.map$Month_Year))

for (year in 2018:2021) { # monthly maps of LISA_C.abs1
  LISA_C.abs1.map %>% filter(substr(Month_Year,10,10) == substr(year,4,4)) %>% 
    ggplot() + geom_polygon(mapping = aes(long, lat, group = GEOID, fill=LISA_C),
                            color = "#000000", size = .01)   +
    geom_point(aes(x=long, y=lat),
               size=2, pch=19, color="green",
               data=centroids.HIDTA) +
    xlim(-132, -68) +
    facet_wrap(. ~ Month_Year) +
    scale_fill_manual(values = c("Insig"="grey90",
                                 "LL"="blue",
                                 "LH"="steelblue",
                                 "HL"="orange",
                                 "HH"="red")) +
    labs(fill = "LISA_C_abs1") -> LISA_C_map
  ggsave(paste("All_Drugs_LISA_C_abs1_", year, "_XCal_HIDTA.jpg", sep=""), LISA_C_map, width=55, height=30, units="cm")
}

LISA_C.abs_t1.map <- left_join(HIDTAs[, c(1:2,7,13)], LISA_C.abs.t_1, by = "HIDTA")
LISA_C.abs_t1.map <- LISA_C.abs_t1.map %>% pivot_longer(c(-long, -lat, -GEOID, -HIDTA), names_to="Month_Year", values_to="LISA_C")
LISA_C.abs_t1.map$LISA_C <- left_join(LISA_C.abs_t1.map[, 5:6], ref, by="LISA_C")$label
LISA_C.abs_t1.map$LISA_C <- factor(LISA_C.abs_t1.map$LISA_C)
LISA_C.abs_t1.map$Month_Year <- factor(LISA_C.abs_t1.map$Month_Year, levels=unique(LISA_C.abs_t1.map$Month_Year))

for (year in 2018:2021) { # monthly maps of LISA_C.abs_t1
  LISA_C.abs_t1.map %>% filter(substr(Month_Year,10,10) == substr(year,4,4)) %>% 
    ggplot() + geom_polygon(mapping = aes(long, lat, group = GEOID, fill=LISA_C),
                            color = "#000000", size = .01)   +
    geom_point(aes(x=long, y=lat),
               size=2, pch=19, color="green",
               data=centroids.HIDTA) +
    xlim(-132, -68) +
    facet_wrap(. ~ Month_Year) +
    scale_fill_manual(values = c("Insig"="grey90",
                                 "LL"="blue",
                                 "LH"="steelblue",
                                 "HL"="orange",
                                 "HH"="red")) +
    labs(fill = "LISA_C_abs_t1") -> LISA_C_map
  ggsave(paste("All_Drugs_LISA_C_abs_t1_", year, "_XCal_HIDTA.jpg", sep=""), LISA_C_map, width=55, height=30, units="cm")
}

