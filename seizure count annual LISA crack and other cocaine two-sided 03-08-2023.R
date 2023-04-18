library(fpp2)
library(spdep)
library(readxl)
library(urbnmapr)
library(tidyverse)
library(gridExtra)
library(lubridate)

cocaine <- read.csv("cocaine other count HIDTA (02-21-2023).csv") %>% as_tibble
crack <- read.csv("cocaine crack count HIDTA (02-21-2023).csv") %>% as_tibble
LISA3 <- read.csv("CountyKNN3.csv") %>% as_tibble %>% arrange(GEOID)
overdose <- read.csv("VSRR_Provisional_County-Level_Drug_Overdose_Death_Counts.csv") %>% as_tibble
names(overdose)[1] <- "state_abbv"
overdose_filtered <- overdose %>% filter(Year %in% 2018:2021 & !(state_abbv %in% c("AK", "HI", "US")))
overdose_filtered$Provisional.Drug.Overdose.Deaths <- as.numeric(overdose_filtered$Provisional.Drug.Overdose.Deaths)
overdose_filtered$Month <- as.numeric(overdose_filtered$Month)
overdose_filtered <- overdose_filtered %>% complete(Year, Month)
overdose_filtered <- overdose_filtered %>% 
  mutate(month_year=str_c(Month, Year, sep="-") %>% 
           parse_date("%m-%Y") %>% format("%b_%Y"))
# overdose_filtered$month_year <- factor(overdose_filtered$month_year, levels=names(cocaine)[-(1:4)])
overdose_filtered <- overdose_filtered %>% rename(state=STATE_NAME, county=COUNTYNAME, GEOID=FIPS)


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
nrow(crack) # 1518
seizures.crack <- crack[, 5:52]

grep("LISA_I", names(crack))[1] # 53

relative.MoransI <- crack
set.seed(100)

localM_2020 <- seizures.crack %>%
  select(Jan_2020:Dec_2020) %>%
  apply(1, sum) %>% 
  localmoran_abs(nb.obj.crack, nsim=nperm, zero.policy=T, xx=NULL, alternative="two.sided")
localM_2020$LISA_C <- ifelse(localM_2020$quadr=="High-High", 1,
                              ifelse(localM_2020$quadr=="Low-Low", 2,
                                     ifelse(localM_2020$quadr=="Low-High", 3, 4)))
localM_2020$LISA_C <- ifelse(localM_2020$`Pr(folded) Sim` <= alpha, localM_2020$LISA_C, 0)
relative.MoransI_2020 <- cbind(relative.MoransI[,1:4], localM_2020[,c(1,13,7)]) %>%
  as_tibble %>% 
  rename(LISA_I=Ii, LISA_C=LISA_C, LISA_P=`Pr(folded) Sim`)

# for other years
# for (year in 2018:2021) {
#   localM_annual <- seizures.crack %>%
#     select(paste("Jan", year):Dec_2020) %>%
#     apply(1, sum) %>% 
#     localmoran_abs(nb.obj.crack, nsim=nperm, zero.policy=T, xx=NULL, alternative="two.sided")
#   localM_annual$LISA_C <- ifelse(localM_annual$quadr=="High-High", 1,
#                                  ifelse(localM_annual$quadr=="Low-Low", 2,
#                                         ifelse(localM_annual$quadr=="Low-High", 3, 4)))
#   localM_annual$LISA_C <- ifelse(localM_annual$`Pr(folded) Sim` <= alpha, localM_annual$LISA_C, 0)
#   relative.MoransI_2020 <- cbind(relative.MoransI[,1:4], localM_annual[,c(1,13,7)]) %>%
#     as_tibble %>% 
#     rename(LISA_I=Ii, LISA_C=LISA_C, LISA_P=`Pr(folded) Sim`)
# }

# crack Plots
LISA.rel <- relative.MoransI_2020
LISA_C.rel <- LISA.rel[,c(1:3, grep("LISA_C", names(LISA.rel)))]
county.names <- LISA.rel[, 1:3]

LISA_C.rel.map <- left_join(counties.obs[, c(1:2,6:7,10,12)], LISA_C.rel[,-(1:2)], by = "GEOID")
ref <- tibble(LISA_C=0:4, label=c("Insig", "HH", "LL", "LH", "HL"))
LISA_C.rel.map$LISA_C <- left_join(LISA_C.rel.map[, 6:7], ref, by="LISA_C")$label
LISA_C.rel.map$LISA_C <- as.factor(LISA_C.rel.map$LISA_C)
LISA_C.rel.map %>% 
    ggplot(mapping = aes(long, lat, group = group, fill=LISA_C)) +
    geom_polygon(color = "#000000", size = .05) +
    scale_fill_manual(values = c("Insig"="grey60",
                                 "LL"="blue",
                                 "LH"="steelblue",
                                 "HL"="orange",
                                 "HH"="red"),
                      na.value = "white") +
    labs(fill = "LISA_C_rel", title="Crack LISA Classes Aggregated in 2020") + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) -> LISA_C_map
# ggsave(paste("Crack_Count_LISA_C_rel_", year, " two-sided (annual).jpg", sep=""), LISA_C_map, width=20, height=15, units="cm")

LISA_C.abs1.map <- left_join(counties.obs[, c(1:2,6:7,10,12)], LISA_C.abs.1[,-(1:2)], by = "GEOID")
LISA_C.abs1.map <- LISA_C.abs1.map %>% pivot_longer(c(-long, -lat, -GEOID, -group, -state, -county), names_to="Month_Year", values_to="LISA_C")
LISA_C.abs1.map$LISA_C <- left_join(LISA_C.abs1.map[, 7:8], ref, by="LISA_C")$label
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

LISA_C.abs_t1.map <- left_join(counties.obs[, c(1:2,6:7,10,12)], LISA_C.abs.t_1[,-(1:2)], by = "GEOID")
LISA_C.abs_t1.map <- LISA_C.abs_t1.map %>% pivot_longer(c(-long, -lat, -GEOID, -group, -state, -county), names_to="Month_Year", values_to="LISA_C")
LISA_C.abs_t1.map$LISA_C <- left_join(LISA_C.abs_t1.map[, 7:8], ref, by="LISA_C")$label
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

# write.csv(relative.MoransI, "cocaine other counts KNN5 R codes 999 two-sided relative (02-21-2023).csv", row.names=F)

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

# write.csv(absolute.b1.MoransI, "cocaine other counts KNN5 R codes 999 two-sided absolute base_1 (02-21-2023).csv", row.names=F)

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

# write.csv(absolute.t1.MoransI, "cocaine other counts KNN5 R codes 999 two-sided absolute base_t1 (02-21-2023).csv", row.names=F)

# Other Cocaine Plots
LISA.rel <- read.csv("cocaine other counts KNN5 R codes 999 two-sided relative (02-21-2023).csv") %>% as_tibble
LISA.abs.1 <- read.csv("cocaine other counts KNN5 R codes 999 two-sided absolute base_1 (02-21-2023).csv") %>% as_tibble
LISA.abs.t_1 <- read.csv("cocaine other counts KNN5 R codes 999 two-sided absolute base_t1 (02-21-2023).csv") %>% as_tibble

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
# ggsave("other cocaine n_HL annual (rel).jpg", n_HL_map, width=20, height=15, units="cm")

LISA_n_HL.map %>% 
  ggplot(mapping = aes(long, lat, group = group, fill=n_HL_All)) +
  geom_polygon(color = "#000000", size = .05) +
  scale_fill_viridis_c(na.value="white") +
  labs(fill = "# of HL", title="# of HL 2018-2021 (other cocaine )") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) -> n_HL_map_all
# ggsave("other cocaine n_HL 2018-2021 (rel).jpg", n_HL_map_all, width=20, height=15, units="cm")

LISA_C.rel.map <- left_join(counties.obs[, c(1:2,6:7,10,12)], LISA_C.rel[,-(1:2)], by = "GEOID")
LISA_C.rel.map <- LISA_C.rel.map %>% pivot_longer(c(-long, -lat, -GEOID, -group, -state, -county), names_to="Month_Year", values_to="LISA_C")
LISA_C.rel.map$LISA_C <- left_join(LISA_C.rel.map[, 7:8], ref, by="LISA_C")$label
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

LISA_C.abs1.map <- left_join(counties.obs[, c(1:2,6:7,10,12)], LISA_C.abs.1[,-(1:2)], by = "GEOID")
LISA_C.abs1.map <- LISA_C.abs1.map %>% pivot_longer(c(-long, -lat, -GEOID, -group, -state, -county), names_to="Month_Year", values_to="LISA_C")
LISA_C.abs1.map$LISA_C <- left_join(LISA_C.abs1.map[, 7:8], ref, by="LISA_C")$label
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

LISA_C.abs_t1.map <- left_join(counties.obs[, c(1:2,6:7,10,12)], LISA_C.abs.t_1[,-(1:2)], by = "GEOID")
LISA_C.abs_t1.map <- LISA_C.abs_t1.map %>% pivot_longer(c(-long, -lat, -GEOID, -group, -state, -county), names_to="Month_Year", values_to="LISA_C")
LISA_C.abs_t1.map$LISA_C <- left_join(LISA_C.abs_t1.map[, 7:8], ref, by="LISA_C")$label
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

# check seizure counts in SW area
for (drug in c("crack", "cocaine other", "meth", "marijuana")) {
  LISA.rel <- read.csv(paste(drug, "counts KNN5 R codes 999 two-sided relative (02-21-2023).csv")) %>% as_tibble
  LISA.rel.map <- left_join(counties.obs[, c(1:2,6:7,10,12)], LISA.rel[, c(3, 5:52)], by = "GEOID")
  LISA.rel.map <- LISA.rel.map %>% 
    pivot_longer(c(-long, -lat, -group, -GEOID, -state, -county), names_to="Month_Year", values_to="seizure_counts") %>% 
    mutate(Year=as.integer(substr(Month_Year, 5,8)))
  LISA.rel.map <- left_join(LISA.rel.map, populations[, -(1:2)], by=c("GEOID", "Year"))
  LISA.rel.map$counts_per_pop <- LISA.rel.map$seizure_counts/LISA.rel.map$population
  LISA.rel.map$Month_Year <- parse_date(LISA.rel.map$Month_Year, "%b_%Y")
  
  seizure_counts <- aggregate(cbind(long, lat) ~ GEOID+seizure_counts,
                              data=LISA.rel.map %>% filter(state %in% c("Arizona", "California", "Nevada") & Month_Year=="2020-01-01"), 
                              FUN=function(x) mean(range(x)))
  LISA.rel.map %>% filter(state %in% c("Arizona", "California", "Nevada")) %>% 
    ggplot(mapping = aes(x=long, y=lat)) +
    geom_polygon(aes(group = group), fill=NA, color = "#000000", size = .05) +
    geom_text(data=seizure_counts, aes(x=long, y=lat, label=seizure_counts), size=2) +
    labs(title=paste(drug, "seizure counts in Jan 2020"))+
    theme_bw() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) -> seizure_map
  ggsave(paste(drug, "seizure counts of AZ-CA-NV Jan_2020.jpg"), seizure_map, width=8, height=10, units="cm")
}

for (drug in c("crack", "cocaine other", "meth", "marijuana")) {
  LISA.rel <- read.csv(paste(drug, "counts KNN5 R codes 999 two-sided (relative).csv")) %>% as_tibble
  LISA.rel.map <- left_join(counties.obs[, c(1:2,6:7,10,12)], LISA.rel[, c(3, 5:52)], by = "GEOID")
  LISA.rel.map <- LISA.rel.map %>% 
    pivot_longer(c(-long, -lat, -group, -GEOID, -state, -county), names_to="Month_Year", values_to="seizure_counts") %>% 
    mutate(Year=as.integer(substr(Month_Year, 5,8)))
  LISA.rel.map <- left_join(LISA.rel.map, populations[, -(1:2)], by=c("GEOID", "Year"))
  LISA.rel.map$counts_per_pop <- LISA.rel.map$seizure_counts/LISA.rel.map$population
  LISA.rel.map$Month_Year <- parse_date(LISA.rel.map$Month_Year, "%b_%Y")
  
  seizure_counts <- aggregate(cbind(long, lat) ~ GEOID+seizure_counts,
                              data=LISA.rel.map %>% filter(state %in% c("Arizona", "California", "Nevada") & Month_Year=="2020-01-01"), 
                              FUN=function(x) mean(range(x)))
  LISA.rel.map %>% filter(state %in% c("Arizona", "California", "Nevada")) %>% 
    ggplot(mapping = aes(x=long, y=lat)) +
    geom_polygon(aes(group = group), fill=NA, color = "#000000", size = .05) +
    geom_text(data=seizure_counts, aes(x=long, y=lat, label=seizure_counts), size=2) +
    labs(title=paste(drug, "seizure counts in Jan 2020 (less missing)"))+
    theme_bw() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) -> seizure_map
  ggsave(paste(drug, "seizure counts of AZ-CA-NV Jan_2020 (less missing).jpg"), seizure_map, width=8, height=10, units="cm")
}

crack <- read.csv(paste("crack", "counts KNN5 R codes 999 two-sided relative (02-21-2023).csv")) %>% as_tibble
crack_less <- read.csv(paste("crack", "counts KNN5 R codes 999 two-sided (relative).csv")) %>% as_tibble
crack_less %>% filter(!(GEOID %in% crack$GEOID))
