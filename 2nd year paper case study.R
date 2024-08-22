# setwd("/Users/euseongjang/Documents/R")
# setwd("C:/Users/gkfrj/Documents/R/Improved LISA")
library(fpp2)
library(spdep)
library(readxl)
library(urbnmapr)
library(tidyverse)
library(gridExtra)
library(lubridate)
library(xtable)

crack <- read.csv("cocaine crack count HIDTA (06-28-2023).csv") %>% as_tibble
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
nb.obj.crack <- nb2listw(nb_crack, style="B")

alpha <- 0.05
nperm <- 999

## For crack
crack <- cbind(crack, matrix(0, nrow(crack), 48*3)) %>% as_tibble

names(crack)[(which(names(crack) == "1"):which(names(crack) == "144"))] <- names(LISA3)[71:214]
Jan_2018_index <- grep("Jan_2018", names(crack))[1]
LISA_I.index <- grep("LISA_I", names(crack))[1]
nrow(crack) # 1518
seizures.crack <- crack[, Jan_2018_index:(LISA_I.index-1)]

grep("LISA_I", names(crack))[1] # 53

original.MoransI <- crack
set.seed(100)
t1 <- Sys.time()
# for (i in 1:48) {
  seizure.crack <- t(seizures.crack)[i,]
  localM.month <- localmoran_abs(seizure.crack, nb.obj.crack, nsim=nperm, zero.policy=T, xx=NULL, alternative="two.sided")
  localM.month$LISA_C <- as.character(localM.month$quadr_ps)
  localM.month$LISA_C <- ifelse(localM.month$`Pr(folded) Sim` <= alpha, localM.month$LISA_C, "Insig")
  original.MoransI[,(LISA_I.index+3*i-3):(LISA_I.index+3*i-1)] <- localM.month[,c(1,13,7)]
# }
t2 <- Sys.time()
# write.csv(original.MoransI, "crack counts KNN5 R codes 999 two-sided original (08-01-2023).csv", row.names=F)

moderate.MoransI <- crack
x.bar_p <- mean(moderate.MoransI$Jan_2018)
set.seed(100)
t3 <- Sys.time()
for (i in 1:48) {
  seizure.crack <- t(seizures.crack)[i,]
  localM.month <- localmoran_abs(seizure.crack, nb.obj.crack, nsim=nperm, zero.policy=T, xx=NULL, alternative="two.sided", moderate=T)
  localM.month$LISA_C <- as.character(localM.month$quadr_ps)
  localM.month$LISA_C <- ifelse(localM.month$`Pr(folded) Sim` <= alpha, localM.month$LISA_C, "Insig")
  moderate.MoransI[,(LISA_I.index+3*i-3):(LISA_I.index+3*i-1)] <- localM.month[,c(1,13,7)]
}
t4 <- Sys.time()
# write.csv(moderate.MoransI, "crack counts KNN5 R codes 999 two-sided moderate (08-01-2023).csv", row.names=F)

perm.i.MoransI <- crack
x.bar_p <- mean(perm.i.MoransI$Jan_2018)
set.seed(100)
t5 <- Sys.time()
# for (i in 1:48) {
  seizure.crack <- t(seizures.crack)[i,]
  localM.month <- localmoran_abs(seizure.crack, nb.obj.crack, nsim=nperm, zero.policy=T, xx=NULL, alternative="two.sided", perm.i=T)
  localM.month$LISA_C <- as.character(localM.month$quadr_ps)
  localM.month$LISA_C <- ifelse(localM.month$`Pr(folded) Sim` <= alpha, localM.month$LISA_C, "Insig")
  perm.i.MoransI[,(LISA_I.index+3*i-3):(LISA_I.index+3*i-1)] <- localM.month[,c(1,13,7)]
# }
t6 <- Sys.time()
# write.csv(perm.i.MoransI, "crack counts KNN5 R codes 999 two-sided permute i (08-01-2023).csv", row.names=F)

both.MoransI <- crack
x.bar_p <- mean(both.MoransI$Jan_2018)
set.seed(100)
for (i in 1:48) {
  seizure.crack <- t(seizures.crack)[i,]
  localM.month <- localmoran_abs(seizure.crack, nb.obj.crack, nsim=nperm, zero.policy=T, xx=NULL, alternative="two.sided", moderate=T, perm.i=T)
  localM.month$LISA_C <- as.character(localM.month$quadr_ps)
  localM.month$LISA_C <- ifelse(localM.month$`Pr(folded) Sim` <= alpha, localM.month$LISA_C, "Insig")
  both.MoransI[,(LISA_I.index+3*i-3):(LISA_I.index+3*i-1)] <- localM.month[,c(1,13,7)]
}

# write.csv(both.MoransI, "crack counts KNN5 R codes 999 two-sided moderate permute i (08-01-2023).csv", row.names=F)

# crack Plots
LISA.org <- read.csv("crack counts KNN5 R codes 999 two-sided original (08-01-2023).csv") %>% as_tibble
LISA.mod <- read.csv("crack counts KNN5 R codes 999 two-sided moderate (08-01-2023).csv") %>% as_tibble
LISA.perm.i <- read.csv("crack counts KNN5 R codes 999 two-sided permute i (08-01-2023).csv") %>% as_tibble
LISA.both <- read.csv("crack counts KNN5 R codes 999 two-sided moderate permute i (08-01-2023).csv") %>% as_tibble

LISA_C.org <- LISA.org[,c(1:3, grep("LISA_C", names(LISA.org)))]
LISA_C.mod <- LISA.mod[,c(1:3, grep("LISA_C", names(LISA.mod)))]
LISA_C.perm.i <- LISA.perm.i[,c(1:3, grep("LISA_C", names(LISA.perm.i)))]
LISA_C.both <- LISA.both[,c(1:3, grep("LISA_C", names(LISA.both)))]
county.names <- LISA_C.org[, 1:3]

LISA_C_to_Month <- function(col_names) {
  col_names <- gsub("0", "2020", col_names)
  col_names <- gsub("1", "2021", col_names)
  col_names <- gsub("8", "2018", col_names)
  col_names <- gsub("9", "2019", col_names)
  col_names <- parse_date(col_names, "LISA_C%b%Y") %>% substr(1,7)
  return(col_names)
}

names(LISA_C.org)[4:51] <- LISA_C_to_Month(names(LISA_C.org)[4:51])
names(LISA_C.mod)[4:51] <- LISA_C_to_Month(names(LISA_C.mod)[4:51])
names(LISA_C.perm.i)[4:51] <- LISA_C_to_Month(names(LISA_C.perm.i)[4:51])
names(LISA_C.both)[4:51] <- LISA_C_to_Month(names(LISA_C.both)[4:51])


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

{
LISA_C.org.map <- left_join(counties.obs[, c(1:2,6:7,10,12)], LISA_C.org[,-(1:2)], by = "GEOID")
LISA_C.org.map <- LISA_C.org.map %>% pivot_longer(c(-long, -lat, -GEOID, -group, -state, -county), names_to="Month_Year", values_to="LISA_C")
LISA_C.org.map$LISA_C <- factor(LISA_C.org.map$LISA_C, levels=c("HH", "HL", "LH", "LL", "Insig"))
LISA_C.org.map$Month_Year <- parse_date(LISA_C.org.map$Month_Year, "%Y-%m")

LISA_C.mod.map <- left_join(counties.obs[, c(1:2,6:7,10,12)], LISA_C.mod[,-(1:2)], by = "GEOID")
LISA_C.mod.map <- LISA_C.mod.map %>% pivot_longer(c(-long, -lat, -GEOID, -group, -state, -county), names_to="Month_Year", values_to="LISA_C")
LISA_C.mod.map$LISA_C <- factor(LISA_C.mod.map$LISA_C, levels=c("HH", "HL", "MH", "ML", "LH", "LL", "Insig"))
LISA_C.mod.map$Month_Year <- parse_date(LISA_C.mod.map$Month_Year, "%Y-%m")

LISA_C.perm.i.map <- left_join(counties.obs[, c(1:2,6:7,10,12)], LISA_C.perm.i[,-(1:2)], by = "GEOID")
LISA_C.perm.i.map <- LISA_C.perm.i.map %>% pivot_longer(c(-long, -lat, -GEOID, -group, -state, -county), names_to="Month_Year", values_to="LISA_C")
LISA_C.perm.i.map$LISA_C <- factor(LISA_C.perm.i.map$LISA_C, levels=c("HH", "HL", "LH", "LL", "Insig"))
LISA_C.perm.i.map$Month_Year <- parse_date(LISA_C.perm.i.map$Month_Year, "%Y-%m")

LISA_C.both.map <- left_join(counties.obs[, c(1:2,6:7,10,12)], LISA_C.both[,-(1:2)], by = "GEOID")
LISA_C.both.map <- LISA_C.both.map %>% pivot_longer(c(-long, -lat, -GEOID, -group, -state, -county), names_to="Month_Year", values_to="LISA_C")
LISA_C.both.map$LISA_C <- factor(LISA_C.both.map$LISA_C, levels=c("HH", "HL", "MH", "ML", "LH", "LL", "Insig"))
LISA_C.both.map$Month_Year <- parse_date(LISA_C.both.map$Month_Year, "%Y-%m")
}


# comparison of results for a month, Jan 2020
LISA_C.org.map %>% filter(Month_Year == "2020-01-01") %>% 
  ggplot(mapping = aes(long, lat, group = group, fill=LISA_C)) +
  geom_polygon(color = "#000000", linewidth = .05) +
  scale_fill_manual(values = c("Insig"="grey60",
                               "LL"="blue",
                               "LH"="steelblue",
                               "HL"="orange",
                               "HH"="red"),
                    na.value = "white") +
  labs(fill = "", title="orignial", x="", y="") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) -> LISA_C_org_map

LISA_C.mod.map %>% filter(Month_Year == "2020-01-01") %>% 
  ggplot(mapping = aes(long, lat, group = group, fill=LISA_C)) +
  geom_polygon(color = "#000000", linewidth = .05) +
  scale_fill_manual(values = c("Insig"="grey60",
                               "LL"="blue",
                               "LH"="steelblue",
                               "ML"="#abd9e9",
                               "MH"="#fee090",
                               "HL"="orange",
                               "HH"="red"),
                    na.value = "white") +
  labs(fill = "", title="moderate", x="", y="") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) -> LISA_C_mod_map

LISA_C.perm.i.map %>% filter(Month_Year == "2020-01-01") %>% 
  ggplot(mapping = aes(long, lat, group = group, fill=LISA_C)) +
  geom_polygon(color = "#000000", linewidth = .05) +
  scale_fill_manual(values = c("Insig"="grey60",
                               "LL"="blue",
                               "LH"="steelblue",
                               "HL"="orange",
                               "HH"="red"),
                    na.value = "white") +
  labs(fill = "", title="permute i", x="", y="") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) -> LISA_C_perm.i_map

LISA_C.both.map %>% filter(Month_Year == "2020-01-01") %>% 
  ggplot(mapping = aes(long, lat, group = group, fill=LISA_C)) +
  geom_polygon(color = "#000000", linewidth = .05) +
  scale_fill_manual(values = c("Insig"="grey60",
                               "LL"="blue",
                               "LH"="steelblue",
                               "ML"="#abd9e9",
                               "MH"="#fee090",
                               "HL"="orange",
                               "HH"="red"),
                    na.value = "white") +
  labs(fill = "", title="both", x="", y="") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) -> LISA_C_both_map

grid.arrange(LISA_C_org_map, LISA_C_mod_map, LISA_C_perm.i_map, LISA_C_both_map, ncol=2)
# ggsave("Crack_Count_LISA_C original two-sided (Jan 2020).pdf", LISA_C_org_map, width=15, height=10, units="cm")
# ggsave("Crack_Count_LISA_C moderate two-sided (Jan 2020).pdf", LISA_C_mod_map, width=15, height=10, units="cm")
# ggsave("Crack_Count_LISA_C permute i two-sided (Jan 2020).pdf", LISA_C_perm.i_map, width=15, height=10, units="cm")
# ggsave("Crack_Count_LISA_C moderate permute i two-sided (Jan 2020).pdf", LISA_C_both_map, width=15, height=10, units="cm")

# NE states
NE_states <- c("Michigan", "Illinois", "Ohio", "Indiana", "Pennsylvania", "New York", "New Hampshire", "New Jersey", "Rhode Island", "Massachusetts",
               "Maryland", "Maine", "District of Columbia", "Delaware", "Vermont", "Virginia", "West Virginia", "Wisconsin", "Connecticut", "Kentucky")

LISA_C.org.map %>% filter(state=="Massachusetts") %>%  select(county) %>% pull %>% unique()

Baltimore_centroid <- LISA_C.org.map %>% filter(state == "Maryland" & county == "Baltimore city") %>% select(long, lat) %>% apply(2, mean)
Boston_centroid <- LISA_C.org.map %>% filter(state == "Massachusetts" & county == "Suffolk County") %>% select(long, lat) %>% apply(2, mean)
Chicago_centroid <- LISA_C.org.map %>% filter(state == "Illinois" & county == "Cook County") %>% select(long, lat) %>% apply(2, mean)
Detroit_centroid <- LISA_C.org.map %>% filter(state == "Michigan" & county == "Wayne County") %>% select(long, lat) %>% apply(2, mean)
Pittsburgh_centroid <- LISA_C.org.map %>% filter(state == "Pennsylvania" & county == "Allegheny County") %>% select(long, lat) %>% apply(2, mean)
VA_WV_entroid <- LISA_C.org.map %>% filter(state == "Virginia" & GEOID == 51161| # Roanoke County
                                           state == "West Virginia" & GEOID %in% c(54055, 54081) ) %>% select(long, lat) %>% apply(2, mean) # Mercer, Raleigh County

LISA_C.perm.i %>% filter(state %in% c("Virginia", "West Virginia") & `2020-01` != "Insig") %>% select(state,county,GEOID,`2020-01`)

NE_cities_centroid <- rbind(Baltimore_centroid, Boston_centroid, Chicago_centroid, Detroit_centroid, Pittsburgh_centroid, VA_WV_entroid) %>% as.data.frame
NE_cities_centroid$LISA_C <- rep("1",6)

NE_cities_centroid_1 <- NE_cities_centroid[1,-3]
NE_cities_centroid_2 <- NE_cities_centroid[2,-3]
NE_cities_centroid_3 <- NE_cities_centroid[3,-3]
NE_cities_centroid_4 <- NE_cities_centroid[4,-3]
NE_cities_centroid_5 <- NE_cities_centroid[5,-3]
NE_cities_centroid_6 <- NE_cities_centroid[6,-3]

NE_cities_centroid_1 <- rbind(NE_cities_centroid_1, NE_cities_centroid_1 + 0.9*c(5,0))
NE_cities_centroid_2 <- rbind(NE_cities_centroid_2, NE_cities_centroid_2 + 0.8*c(2,1))
NE_cities_centroid_3 <- rbind(NE_cities_centroid_3, NE_cities_centroid_3 + 0.9*c(-5,0))
NE_cities_centroid_4 <- rbind(NE_cities_centroid_4, NE_cities_centroid_4 + 0.8*c(2,1))
NE_cities_centroid_5 <- rbind(NE_cities_centroid_5, NE_cities_centroid_5 + 0.9*c(0,4))

LISA_C.org.map %>% filter(Month_Year == "2020-01-01" & state %in% NE_states) %>% 
  ggplot(mapping = aes(long, lat)) +
  geom_polygon(mapping = aes(group = group, fill=LISA_C), color = "#000000", linewidth = .05) +
  scale_fill_manual(values = c("Insig"="grey60",
                               "LL"="blue",
                               "LH"="steelblue",
                               "HL"="orange",
                               "HH"="red"),
                    na.value = "white") +
  labs(fill = "", title="", x="", y="") + 
  geom_path(data=NE_cities_centroid_1,
            aes(x=long, y=lat),
            linewidth=0.3) +
  geom_path(data=NE_cities_centroid_2,
            aes(x=long, y=lat),
            linewidth=0.3) +
  geom_path(data=NE_cities_centroid_3,
            aes(x=long, y=lat),
            linewidth=0.3) +
  geom_path(data=NE_cities_centroid_4,
            aes(x=long, y=lat),
            linewidth=0.3) +
  geom_path(data=NE_cities_centroid_5,
            aes(x=long, y=lat),
            linewidth=0.3) +
  geom_text(data=NE_cities_centroid[c(2,4),],
            aes(x=long, y=lat, group=LISA_C, label=as.character(c(2,4))),
            nudge_x = 2, nudge_y = 1,
            size=5) +
  geom_text(data=NE_cities_centroid[1,],
            aes(x=long, y=lat, group=LISA_C, label="1"),
            nudge_x = 5,
            size=5) +
  geom_text(data=NE_cities_centroid[3,],
            aes(x=long, y=lat, group=LISA_C, label="3"),
            nudge_x = -5,
            size=5) +
  geom_text(data=NE_cities_centroid[5,],
            aes(x=long, y=lat, group=LISA_C, label="5"),
            nudge_y = 4,
            size=5) +
  geom_text(data=NE_cities_centroid[6,],
            aes(x=long, y=lat, group=LISA_C, label="6"),
            nudge_x = 0.2,
            size=5) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) -> NE_org_map

LISA_C.mod.map %>% filter(Month_Year == "2020-01-01" & state %in% NE_states) %>% 
  ggplot(mapping = aes(long, lat)) +
  geom_polygon(mapping = aes(group = group, fill=LISA_C), color = "#000000", linewidth = .05) +
  scale_fill_manual(values = c("Insig"="grey60",
                               "LL"="blue",
                               "LH"="steelblue",
                               "ML"="#abd9e9",
                               "MH"="#fee090",
                               "HL"="orange",
                               "HH"="red"),
                    na.value = "white") +
  labs(fill = "", title="", x="", y="") + 
  geom_path(data=NE_cities_centroid_1,
            aes(x=long, y=lat),
            linewidth=0.3) +
  geom_path(data=NE_cities_centroid_2,
            aes(x=long, y=lat),
            linewidth=0.3) +
  geom_path(data=NE_cities_centroid_3,
            aes(x=long, y=lat),
            linewidth=0.3) +
  geom_path(data=NE_cities_centroid_4,
            aes(x=long, y=lat),
            linewidth=0.3) +
  geom_path(data=NE_cities_centroid_5,
            aes(x=long, y=lat),
            linewidth=0.3) +
  geom_text(data=NE_cities_centroid[c(2,4),],
            aes(x=long, y=lat, group=LISA_C, label=as.character(c(2,4))),
            nudge_x = 2, nudge_y = 1,
            size=5) +
  geom_text(data=NE_cities_centroid[1,],
            aes(x=long, y=lat, group=LISA_C, label="1"),
            nudge_x = 5,
            size=5) +
  geom_text(data=NE_cities_centroid[3,],
            aes(x=long, y=lat, group=LISA_C, label="3"),
            nudge_x = -5,
            size=5) +
  geom_text(data=NE_cities_centroid[5,],
            aes(x=long, y=lat, group=LISA_C, label="5"),
            nudge_y = 4,
            size=5) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) -> NE_mod_map

LISA_C.perm.i.map %>% filter(Month_Year == "2020-01-01" & state %in% NE_states) %>% 
  ggplot(mapping = aes(long, lat)) +
  geom_polygon(mapping = aes(group = group, fill=LISA_C), color = "#000000", linewidth = .05) +
  scale_fill_manual(values = c("Insig"="grey60",
                               "LL"="blue",
                               "LH"="steelblue",
                               "HL"="orange",
                               "HH"="red"),
                    na.value = "white") +
  labs(fill = "", title="", x="", y="") + 
  geom_path(data=NE_cities_centroid_1,
            aes(x=long, y=lat),
            linewidth=0.3) +
  geom_path(data=NE_cities_centroid_2,
            aes(x=long, y=lat),
            linewidth=0.3) +
  geom_path(data=NE_cities_centroid_3,
            aes(x=long, y=lat),
            linewidth=0.3) +
  geom_path(data=NE_cities_centroid_4,
            aes(x=long, y=lat),
            linewidth=0.3) +
  geom_path(data=NE_cities_centroid_5,
            aes(x=long, y=lat),
            linewidth=0.3) +
  geom_text(data=NE_cities_centroid[c(2,4),],
            aes(x=long, y=lat, group=LISA_C, label=as.character(c(2,4))),
            nudge_x = 2, nudge_y = 1,
            size=5) +
  geom_text(data=NE_cities_centroid[1,],
            aes(x=long, y=lat, group=LISA_C, label="1"),
            nudge_x = 5,
            size=5) +
  geom_text(data=NE_cities_centroid[3,],
            aes(x=long, y=lat, group=LISA_C, label="3"),
            nudge_x = -5,
            size=5) +
  geom_text(data=NE_cities_centroid[5,],
            aes(x=long, y=lat, group=LISA_C, label="5"),
            nudge_y = 4,
            size=5) +
  geom_text(data=NE_cities_centroid[6,],
            aes(x=long, y=lat, group=LISA_C, label="6"),
            nudge_x = 0.2,
            size=5) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) -> NE_perm.i_map

LISA_C.both.map %>% filter(Month_Year == "2020-01-01" & state %in% NE_states) %>% 
  ggplot(mapping = aes(long, lat)) +
  geom_polygon(mapping = aes(group = group, fill=LISA_C), color = "#000000", linewidth = .05) +
  scale_fill_manual(values = c("Insig"="grey60",
                               "LL"="blue",
                               "LH"="steelblue",
                               "ML"="#abd9e9",
                               "MH"="#fee090",
                               "HL"="orange",
                               "HH"="red"),
                    na.value = "white") +
  labs(fill = "", title="", x="", y="") + 
  geom_path(data=NE_cities_centroid_1,
            aes(x=long, y=lat),
            linewidth=0.3) +
  geom_path(data=NE_cities_centroid_2,
            aes(x=long, y=lat),
            linewidth=0.3) +
  geom_path(data=NE_cities_centroid_3,
            aes(x=long, y=lat),
            linewidth=0.3) +
  geom_path(data=NE_cities_centroid_4,
            aes(x=long, y=lat),
            linewidth=0.3) +
  geom_path(data=NE_cities_centroid_5,
            aes(x=long, y=lat),
            linewidth=0.3) +
  geom_text(data=NE_cities_centroid[c(2,4),],
            aes(x=long, y=lat, group=LISA_C, label=as.character(c(2,4))),
            nudge_x = 2, nudge_y = 1,
            size=5) +
  geom_text(data=NE_cities_centroid[1,],
            aes(x=long, y=lat, group=LISA_C, label="1"),
            nudge_x = 5,
            size=5) +
  geom_text(data=NE_cities_centroid[3,],
            aes(x=long, y=lat, group=LISA_C, label="3"),
            nudge_x = -5,
            size=5) +
  geom_text(data=NE_cities_centroid[5,],
            aes(x=long, y=lat, group=LISA_C, label="5"),
            nudge_y = 4,
            size=5) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) -> NE_both_map

# ggsave("Crack_Count_NE original two-sided (Jan 2020).pdf", NE_org_map, width=15, height=10, units="cm")
# ggsave("Crack_Count_NE moderate two-sided (Jan 2020).pdf", NE_mod_map, width=15, height=10, units="cm")
# ggsave("Crack_Count_NE permute i two-sided (Jan 2020).pdf", NE_perm.i_map, width=15, height=10, units="cm")
# ggsave("Crack_Count_NE combined two-sided (Jan 2020).pdf", NE_both_map, width=15, height=10, units="cm")


LISA_C.org.map %>% filter(Month_Year == "2020-01-01" & state == "Florida") %>% 
  ggplot(mapping = aes(long, lat, group = group, fill=LISA_C)) +
  geom_polygon(color = "#000000", linewidth = .05) +
  scale_fill_manual(values = c("Insig"="grey60",
                               "LL"="blue",
                               "LH"="steelblue",
                               "HL"="orange",
                               "HH"="red"),
                    na.value = "white") +
  labs(fill = "", title="", x="", y="") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) -> FL_org_map

LISA_C.mod.map %>% filter(Month_Year == "2020-01-01" & state == "Florida") %>% 
  ggplot(mapping = aes(long, lat, group = group, fill=LISA_C)) +
  geom_polygon(color = "#000000", linewidth = .05) +
  scale_fill_manual(values = c("Insig"="grey60",
                               "LL"="blue",
                               "LH"="steelblue",
                               "ML"="#abd9e9",
                               "MH"="#fee090",
                               "HL"="orange",
                               "HH"="red"),
                    na.value = "white") +
  labs(fill = "", title="", x="", y="") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) -> FL_mod_map

LISA_C.perm.i.map %>% filter(Month_Year == "2020-01-01" & state == "Florida") %>% 
  ggplot(mapping = aes(long, lat, group = group, fill=LISA_C)) +
  geom_polygon(color = "#000000", linewidth = .05) +
  scale_fill_manual(values = c("Insig"="grey60",
                               "LL"="blue",
                               "LH"="steelblue",
                               "HL"="orange",
                               "HH"="red"),
                    na.value = "white") +
  labs(fill = "", title="", x="", y="") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) -> FL_perm.i_map

LISA_C.both.map %>% filter(Month_Year == "2020-01-01" & state == "Florida") %>% 
  ggplot(mapping = aes(long, lat, group = group, fill=LISA_C)) +
  geom_polygon(color = "#000000", linewidth = .05) +
  scale_fill_manual(values = c("Insig"="grey60",
                               "LL"="blue",
                               "LH"="steelblue",
                               "ML"="#abd9e9",
                               "MH"="#fee090",
                               "HL"="orange",
                               "HH"="red"),
                    na.value = "white") +
  labs(fill = "", title="", x="", y="") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) -> FL_both_map

# ggsave("Crack_Count_FL original two-sided (Jan 2020).pdf", FL_org_map, width=10, height=8, units="cm")
# ggsave("Crack_Count_FL moderate two-sided (Jan 2020).pdf", FL_mod_map, width=10, height=8, units="cm")
# ggsave("Crack_Count_FL permute i two-sided (Jan 2020).pdf", FL_perm.i_map, width=10, height=8, units="cm")
# ggsave("Crack_Count_FL combined two-sided (Jan 2020).pdf", FL_both_map, width=10, height=8, units="cm")


# comparison of summaries for 2018-2021
most_frequent_label <- LISA_C.org[,1:3]
most_frequent_label$original <- LISA_C.org[,-(1:3)] %>% apply(1, table) %>% sapply(function(x) return(names(x)[order(x, decreasing=T)][1]))
most_frequent_label$moderate <- LISA_C.mod[,-(1:3)] %>% apply(1, table) %>% sapply(function(x) return(names(x)[order(x, decreasing=T)][1]))
most_frequent_label$perm.i <- LISA_C.perm.i[,-(1:3)] %>% apply(1, table) %>% sapply(function(x) return(names(x)[order(x, decreasing=T)][1]))
most_frequent_label$both <- LISA_C.both[,-(1:3)] %>% apply(1, table) %>% sapply(function(x) return(names(x)[order(x, decreasing=T)][1]))
most_frequent_label[,-(1:3)] %>% apply(2, table)

# perm.i has the highest number of HH
# moderate and both have slightly higher number of Insignificant due to the introduction of a new label "M",
# which can reduce the number of each significant label and thus can make Insignificant the most frequent label for some counties
changed_counties_moderate <- which(most_frequent_label$original != most_frequent_label$moderate)
changed_counties_perm.i <- which(most_frequent_label$original != most_frequent_label$perm.i)
changed_counties_both <- which(most_frequent_label$original != most_frequent_label$both)

most_frequent_label[changed_counties_moderate, -c(6, 7)]
most_frequent_label[changed_counties_perm.i, -c(5, 7)]
most_frequent_label[changed_counties_both, -c(5, 6)]

LISA_C.org[,-(1:3)] %>% flatten %>% unlist %>% table
LISA_C.mod[,-(1:3)] %>% flatten %>% unlist %>% table
LISA_C.perm.i[,-(1:3)] %>% flatten %>% unlist %>% table
LISA_C.both[,-(1:3)] %>% flatten %>% unlist %>% table

perm.i_HH_index <- which(LISA_C.perm.i == "HH", arr.ind=T) %>% as.data.frame
label_check <- data.frame(org_label=perm.i_HH_index %>% apply(1, function(x) return(LISA_C.org[x[1], x[2]])) %>% unlist)
label_check$perm.i_label <- perm.i_HH_index %>% apply(1, function(x) return(LISA_C.perm.i[x[1], x[2]])) %>% unlist
new_HH_index <- perm.i_HH_index[with(label_check, which(org_label == "Insig" & perm.i_label == "HH")),]
new_HH_seizure_counts <- data.frame(xi=new_HH_index %>% apply(1, function(x) return(LISA.perm.i[x[1], x[2]+1])) %>% unlist, row.names=NULL)
new_HH_seizure_counts$sum_of_neighbors <- new_HH_index %>% apply(1, function(x) return(
  lag.listw(nb.obj.crack, LISA.perm.i[,x[2]+1] %>% pull, zero.policy=T)[x[1]]
  ))
new_HH_seizure_counts$month <- new_HH_index %>% apply(1, function(x) return(names(LISA.perm.i)[x[2]+1])) %>% unlist

head(new_HH_index)
head(LISA.org[head(new_HH_index$row), 5] %>% pull)
lag.listw(nb.obj.crack, LISA.perm.i[,5] %>% pull, zero.policy=T)[head(new_HH_index$row)]
head(new_HH_seizure_counts)


2*sd(LISA.perm.i$Jan_2020) # two standard error = 8.724866

summary(lag.listw(nb.obj.crack, LISA.perm.i$Jan_2020, zero.policy=T))
new_HH_seizure_counts %>% 
  filter(month == "Jan_2020") %>%
  ggplot(aes(x=sum_of_neighbors, y=xi)) +
  geom_point() + xlim(0, 170)

summary(lag.listw(nb.obj.crack, LISA.perm.i$Jan_2021, zero.policy=T))
new_HH_seizure_counts %>% 
  filter(month == "Jan_2021") %>%
  ggplot(aes(x=sum_of_neighbors, y=xi)) +
  geom_point() + xlim(0, 46)

label_check_Jan_2020 <- tibble(org_label=LISA_C.org$`2020-01`,
                               perm.i_label=LISA_C.perm.i$`2020-01`,
                               moderate_label=LISA_C.mod$`2020-01`,
                               both_label=LISA_C.both$`2020-01`)

label_check_Jan_2020$state <- LISA.org$state
label_check_Jan_2020$county <- LISA.org$county
label_check_Jan_2020$seizure_count <- LISA.org$Jan_2020
label_check_Jan_2020$sum_of_neighbors <- lag.listw(nb.obj.crack, LISA.org$Jan_2020, zero.policy=T)

with(label_check_Jan_2020, org_label[which(org_label != perm.i_label)]) %>% table
with(label_check_Jan_2020, perm.i_label[which(org_label != perm.i_label)]) %>% table

label_check_Jan_2020 %>%
  filter(org_label != perm.i_label) %>% 
  select(state, county, org_label, perm.i_label, seizure_count, sum_of_neighbors) %>% 
  arrange(org_label, perm.i_label, seizure_count, sum_of_neighbors) %>% 
  as.data.frame -> Jan_2020_table

label_check_Jan_2020$org_label %>% table
label_check_Jan_2020$moderate_label %>% table
label_check_Jan_2020$perm.i_label %>% table
label_check_Jan_2020$both_label %>% table

names(Jan_2020_table) <- c("State", "County", "Original", "Permute i", "Seizure Count ($x_i$)", "Sum of neighbors")
Jan_2020_table
xtable(Jan_2020_table %>% filter(Original == "HL"))
xtable(Jan_2020_table %>% filter(Original == "LH"))
xtable(Jan_2020_table %>% filter(Original == "Insig"))

label_check_Jan_2020 %>%
  filter(org_label == "HL" | perm.i_label == "HL") %>% 
  select(state, county, org_label, perm.i_label, seizure_count, sum_of_neighbors) %>% 
  arrange(org_label, perm.i_label, seizure_count, sum_of_neighbors) %>% 
  as.data.frame

label_check_Jan_2020 %>% filter(seizure_count > 8 & sum_of_neighbors > 3) %>% nrow # 31
label_check_Jan_2020 %>% filter(org_label != "Insig" & seizure_count > 8 & sum_of_neighbors > 3) %>% nrow # 17
label_check_Jan_2020 %>% filter(perm.i_label != "Insig" & seizure_count > 8 & sum_of_neighbors > 3) %>% nrow # 30



GEOID_LH_to_Insig <- LISA.org$GEOID[which(LISA.org$LISA_CJan0 == "LH" & LISA.perm.i$LISA_CJan0 == "Insig")]
LISA.org %>% filter(GEOID %in% GEOID_LH_to_Insig) %>% select(state, county)

LISA_C.org.map %>% filter(Month_Year == "2020-01-01") %>% 
  mutate(LISA_C = as.character(LISA_C)) %>% 
  mutate(LISA_C = as.factor(ifelse(GEOID %in% GEOID_LH_to_Insig, "LH to Insig", LISA_C))) %>% 
  ggplot(mapping = aes(long, lat, group = group, fill=LISA_C)) +
  geom_polygon(color = "#000000", linewidth = .05) +
  scale_fill_manual(values = c("Insig"="grey60",
                               "LL"="blue",
                               "LH"="steelblue",
                               "HL"="orange",
                               "HH"="red",
                               "LH to Insig"="green"),
                    na.value = "white") +
  labs(fill = "", title="orignial", x="", y="") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) -> LISA_C_org_map_test

LISA_C.perm.i.map %>% filter(Month_Year == "2020-01-01") %>% 
  mutate(LISA_C = as.character(LISA_C)) %>% 
  mutate(LISA_C = as.factor(ifelse(GEOID %in% GEOID_LH_to_Insig, "LH to Insig", LISA_C))) %>% 
  ggplot(mapping = aes(long, lat, group = group, fill=LISA_C)) +
  geom_polygon(color = "#000000", linewidth = .05) +
  scale_fill_manual(values = c("Insig"="grey60",
                               "LL"="blue",
                               "LH"="steelblue",
                               "HL"="orange",
                               "HH"="red",
                               "LH to Insig"="green"),
                    na.value = "white") +
  labs(fill = "", title="permute i", x="", y="") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) -> LISA_C_perm.i_map_test

LH_to_Insig_check <- grid.arrange(LISA_C_org_map_test, LISA_C_perm.i_map_test, ncol=2)
ggsave("LH to Insig check.pdf", plot=LH_to_Insig_check, width=30, height=10, units="cm")
