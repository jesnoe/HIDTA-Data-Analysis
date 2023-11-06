# setwd("C:/Users/gkfrj/Documents/R")
library(fpp2)
library(spdep)
library(readxl)
library(urbnmapr)
library(tidyverse)
library(gridExtra)
library(lubridate)
library(ggbreak)

# urbnmapr data
unique(countydata$county_fips) # some ips have 5 digits. need to be converted into integer.
HIDTA.dist <- read.csv("HIDTA Regions.csv", header=T) %>% as_tibble
coordinate.HIDTA <- left_join(counties, HIDTA.dist, by=c("state_name", "county_name")) %>% filter(!(state_name %in% c("Alaska", "Hawaii", "Puerto Rico")))
names(coordinate.HIDTA)[7] <- "GEOID"
coordinate.HIDTA$county <- str_split(coordinate.HIDTA$county_name, " ") %>%
  lapply(function(x) str_c(c(x[(1:(length(x)-1))]), collapse=" ") ) %>%
  unlist
HIDTA.regions <- coordinate.HIDTA %>% filter(!is.na(HIDTA)) %>% select(GEOID, HIDTA) %>% unique

seizures <- read_xlsx("Drug Seizures All HIDTAs All Drugs 2018-2021 Combined.xlsx") %>% filter(State %in% unique(coordinate.HIDTA$state_name))
seizures <- seizures %>% filter(County != "Michigan")
seizures$County <- substring(seizures$County, 1, str_locate(seizures$County, ",")[,1]-1)
# seizures$County <- sapply(seizures$County, function(x) paste(x[1], str_to_lower(substring(x, 2, length(x)-1)), sep=""))
seizures$Year <- substring(seizures$SeizureDate, 1, 4)
seizures$Month <- substring(seizures$SeizureDate, 6, 7)
seizures$Year <- as.numeric(seizures$Year)
seizures$Month <- as.numeric(seizures$Month)
seizures$County[seizures$State == "District of Columbia"] <- "District of"

names(seizures)[1:2] <- c("state_name", "county")

all_counts <- seizures %>% left_join(unique(coordinate.HIDTA[,c(7, 12, 14, 15)]), by=c("state_name", "county")) %>% 
  group_by(state_name, county, GEOID, Year, Month) %>%
  summarise(count=n()) %>% 
  ungroup %>%
  select(GEOID:count) %>%
  complete(GEOID, Year, Month)

all_counts <- all_counts %>% mutate(Month2=month(Month, label=T))
all_counts_missingness <- all_counts %>% select(-Month) %>% 
  pivot_wider(names_from=c(Month2, Year), values_from=count) %>% 
  arrange(GEOID)
all_counts_missingness$missing <- all_counts_missingness[,-1] %>%
  apply(1, function(x) ifelse(sum(x, na.rm=T) == 0, 1, 0)) # 1 if missing (no observation for 48 months)
non_missing_GEOID <- all_counts_missingness$GEOID[all_counts_missingness$missing == 0]

# All drugs except Cannabinoids
{
drugs <- seizures$Drug %>% unique %>% sort
Cannabinoids <- c("Cannabis", "Hashish", grep("Marijuana", drugs, value=T), grep("THC", drugs, value=T))
all_drugs <- seizures %>% filter(!(Drug %in% Cannabinoids))
all_drugs <- all_drugs %>% mutate(Month2=month(Month, label=T))
calender <- data.frame(Month=1:12, Month2=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
all_drugs <- all_drugs %>% group_by(state_name, county, Year, Month2) %>%
  summarise(count=sum(Quantity > 0), .groups="drop") %>%
  complete(state_name, county, Year, Month2) %>% 
  pivot_wider(names_from=c(Month2, Year), names_sep="_", values_from=count)
cal_order <- c(str_c(calender$Month2, 2018, sep="_"),
               str_c(calender$Month2, 2019, sep="_"),
               str_c(calender$Month2, 2020, sep="_"),
               str_c(calender$Month2, 2021, sep="_"))
all_drugs <- all_drugs[, c(1:2, match(cal_order, names(all_drugs)))]

all_drugs <- left_join(unique(coordinate.HIDTA[,c(7, 12, 14, 15)]), all_drugs, by=c("state_name", "county")) %>%
  select(state_name, county, GEOID, HIDTA, Jan_2018:Dec_2021) %>% arrange(GEOID)
names(all_drugs)[1] <- "state"
}

all_drugs %>% filter(is.na(county))
grep("Jan_2018", names(all_drugs)) # 5
grep("Dec_2021", names(all_drugs)) # 52
all_drugs <- all_drugs %>% filter(GEOID %in% non_missing_GEOID | !is.na(HIDTA))
all_drugs[,5:52][is.na(all_drugs[,5:52])] <- 0
all_drugs %>% filter(is.na(GEOID))
all_drugs # 1518 counties in total.
# write.csv(all_drugs, "all_drugs (X Cannabinoids) count HIDTA (11-05-2023).csv", row.names=F)

# all_drugs LISA
counties.obs <- counties
names(counties.obs)[7] <- "GEOID"
counties.obs <- counties.obs %>% rename(state=state_name, county=county_name)
counties.obs$GEOID <- as.numeric(counties.obs$GEOID)
counties.obs <- counties.obs %>% filter(!(state %in% c("Alaska", "Hawaii")))
coords.all_drugs <- counties.obs %>% filter(GEOID %in% as.numeric(all_drugs$GEOID)) %>%
  group_by(GEOID) %>% summarise(x=mean(long), y=mean(lat))
GEOIDS.all_drugs <- coords.all_drugs$GEOID
coords.all_drugs <- coords.all_drugs[,-1]
nb_all_drugs <- knn2nb(knearneigh(coords.all_drugs, k=5), row.names=GEOIDS.all_drugs)
nb.obj.all_drugs <- nb2listw(nb_all_drugs, style="B")

alpha <- 0.05
nperm <- 9999

all_drugs <- cbind(all_drugs %>% select(-(Jan_2018:Dec_2019)), matrix(0, nrow(all_drugs), 24*3)) %>% as_tibble
LISA.crack <- read.csv("crack counts KNN5 R codes 999 upper tail relative (02-21-2023).csv") %>% as_tibble
which(names(LISA.crack) == "LISA_IJan0")

names(all_drugs)[(which(names(all_drugs) == "1"):which(names(all_drugs) == "72"))] <- names(LISA.crack)[125:196]
Jan_2020_index <- grep("Jan_2020", names(all_drugs))[1]
LISA_I.index <- grep("LISA_I", names(all_drugs))[1]
LISA_I.index # 29
nrow(all_drugs) # 1518
seizures.all_drugs <- all_drugs[, Jan_2020_index:(LISA_I.index-1)]

original.MoransI <- all_drugs
set.seed(100)
for (i in 1:24) {
  seizure.all_drugs <- t(seizures.all_drugs)[i,]
  localM.month <- localmoran_abs(seizure.all_drugs, nb.obj.all_drugs, nsim=nperm, zero.policy=T, xx=NULL, alternative="two.sided")
  localM.month$LISA_C <- as.character(localM.month$quadr_ps)
  localM.month$LISA_C <- ifelse(localM.month$`Pr(folded) Sim` <= alpha, localM.month$LISA_C, "Insig")
  original.MoransI[,(LISA_I.index+3*i-3):(LISA_I.index+3*i-1)] <- localM.month[,c(1,13,7)]
}

# write.csv(original.MoransI, "all_drugs counts KNN5 R codes 999 two-sided LISA original (11-05-2023).csv", row.names=F)

perm.i.MoransI <- all_drugs
x.bar_p <- mean(perm.i.MoransI$Jan_2020)
set.seed(100)
for (i in 1:24) {
  seizure.all_drugs <- t(seizures.all_drugs)[i,]
  localM.month <- localmoran_abs(seizure.all_drugs, nb.obj.all_drugs, nsim=nperm, zero.policy=T, xx=NULL, alternative="two.sided", perm.i=T)
  localM.month$LISA_C <- as.character(localM.month$quadr_ps)
  localM.month$LISA_C <- ifelse(localM.month$`Pr(folded) Sim` <= alpha, localM.month$LISA_C, "Insig")
  perm.i.MoransI[,(LISA_I.index+3*i-3):(LISA_I.index+3*i-1)] <- localM.month[,c(1,13,7)]
}

# write.csv(perm.i.MoransI, "all_drugs counts KNN5 R codes 999 two-sided LISA permute i (11-05-2023).csv", row.names=F)


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

overdose <- read.csv("VSRR_Provisional_County-Level_Drug_Overdose_Death_Counts.csv") %>% as_tibble
names(overdose)[4] <- "state_abbv"
overdose <- overdose %>% filter(!(state_abbv %in% c("AK", "HI", "US")))
overdose$Provisional.Drug.Overdose.Deaths <- as.numeric(overdose$Provisional.Drug.Overdose.Deaths)
overdose$Month <- as.numeric(overdose$Month)
overdose <- overdose %>% complete(Year, Month)
overdose <- overdose %>% 
  mutate(month_year=str_c(Month, Year, sep="-") %>% 
           parse_date("%m-%Y") %>% format("%b_%Y"))
overdose <- overdose %>% rename(state=STATE_NAME, county=COUNTYNAME, GEOID=FIPS)
overdose <- overdose %>% filter(!(county %in% c("Bedford (city)", "Clifton Forge (city)")) & !is.na(GEOID))
overdose <- overdose %>% 
  select(GEOID, month_year, Provisional.Drug.Overdose.Deaths) %>% rename(Month_Year=month_year) %>% 
  rename(overdose=Provisional.Drug.Overdose.Deaths)
overdose$overdose[is.na(overdose$overdose)] <- 5
populations <- populations %>% filter(Year == 2020) %>% select(-Year)

add_drug_info <- function(drug_data, drug_name, data_to_add) {
  drug.LISA_C.rel <- drug_data[,c(1:3, grep("LISA_C", names(drug_data)))]
  drug.LISA_C.rel <- drug.LISA_C.rel[,-(4:27)]
  names(drug.LISA_C.rel)[-(1:3)] <- unique(overdose$Month_Year)
  drug.LISA_C.rel <- drug.LISA_C.rel %>% 
    pivot_longer(-c(state, county, GEOID), names_to="Month_Year", values_to=paste(drug_name, "HH", sep="_")) %>% 
    select(3:5)
  drug.count <- drug_data %>% select(GEOID, Jan_2020:Dec_2021) %>% 
    pivot_longer(-GEOID, names_to="Month_Year", values_to=paste(drug_name, "count", sep="_"))
  result <- left_join(data_to_add, drug.LISA_C.rel, by=c("GEOID", "Month_Year")) %>%
    left_join(drug.count, by=c("GEOID", "Month_Year"))
  return(result)
}

overdose_integrated <- left_join(overdose, populations, by=c("GEOID")) %>% 
  mutate(overdose.per.pop=overdose/population)
  
overdose_integrated <- add_drug_info(all_drugs, "all_drugs", overdose_integrated)

# overdose LISA
overdose_LISA <- overdose_integrated %>%
  filter(substr(Month_Year, 5, 8) != "2022") %>% 
  pivot_wider(c(state, county, GEOID), names_from = Month_Year, values_from = overdose)
coords.overdose <- counties.obs %>% filter(GEOID %in% as.numeric(overdose_LISA$GEOID)) %>%
  group_by(GEOID) %>% summarise(x=mean(long), y=mean(lat))
GEOIDS.overdose <- coords.overdose$GEOID
coords.overdose <- coords.overdose[,-1]
nb_overdose <- knn2nb(knearneigh(coords.overdose, k=5), row.names=GEOIDS.overdose)
nb.obj.overdose <- nb2listw(nb_overdose, style="B")

alpha <- 0.05
nperm <- 9999

overdose <- cbind(overdose_LISA, matrix(0, nrow(overdose_LISA), 24*3)) %>% as_tibble
LISA.crack <- read.csv("crack counts KNN5 R codes 999 upper tail relative (02-21-2023).csv") %>% as_tibble
which(names(LISA.crack) == "LISA_IJan0")

names(overdose)[(which(names(overdose) == "1"):which(names(overdose) == "72"))] <- names(LISA.crack)[125:196]
Jan_2020_index <- grep("Jan_2020", names(overdose))[1]
LISA_I.index <- grep("LISA_I", names(overdose))[1]
LISA_I.index # 29
nrow(overdose) # 3108
seizures.overdose <- overdose[, Jan_2020_index:(LISA_I.index-1)]

original.MoransI <- overdose
set.seed(100)
for (i in 1:24) {
  seizure.overdose <- t(seizures.overdose)[i,]
  localM.month <- localmoran_abs(seizure.overdose, nb.obj.overdose, nsim=nperm, zero.policy=T, xx=NULL, alternative="two.sided")
  localM.month$LISA_C <- as.character(localM.month$quadr_ps)
  localM.month$LISA_C <- ifelse(localM.month$`Pr(folded) Sim` <= alpha, localM.month$LISA_C, "Insig")
  original.MoransI[,(LISA_I.index+3*i-3):(LISA_I.index+3*i-1)] <- localM.month[,c(1,13,7)]
}

# write.csv(original.MoransI, "overdose counts KNN5 R codes 999 two-sided LISA original (11-05-2023).csv", row.names=F)

perm.i.MoransI <- overdose
x.bar_p <- mean(perm.i.MoransI$Jan_2020)
set.seed(100)
for (i in 1:24) {
  seizure.overdose <- t(seizures.overdose)[i,]
  localM.month <- localmoran_abs(seizure.overdose, nb.obj.overdose, nsim=nperm, zero.policy=T, xx=NULL, alternative="two.sided", perm.i=T)
  localM.month$LISA_C <- as.character(localM.month$quadr_ps)
  localM.month$LISA_C <- ifelse(localM.month$`Pr(folded) Sim` <= alpha, localM.month$LISA_C, "Insig")
  perm.i.MoransI[,(LISA_I.index+3*i-3):(LISA_I.index+3*i-1)] <- localM.month[,c(1,13,7)]
}

# write.csv(perm.i.MoransI, "overdose counts KNN5 R codes 999 two-sided LISA permute i (11-05-2023).csv", row.names=F)

# all_drugs seizure count map
coordinate_map <- coordinate.HIDTA
coordinate_map$GEOID <- as.numeric(coordinate_map$GEOID)
names(coordinate_map)[10] <- "county"
names(coordinate_map)[12] <- "state"

all_drugs.org <- read.csv("all_drugs counts KNN5 R codes 999 two-sided LISA original (11-05-2023).csv") %>% as_tibble
all_drugs.perm.i <- read.csv("all_drugs counts KNN5 R codes 999 two-sided LISA permute i (11-05-2023).csv") %>% as_tibble
all_drugs.map <- left_join(coordinate_map[, c(1:2,6:7,10,12)], all_drugs.org[, c(3, 5:52)], by = "GEOID")

all_drugs.map <- all_drugs.map %>% 
  select(long:state, Jan_2020) %>% 
  pivot_longer(c(-long, -lat, -group, -GEOID, -state, -county), names_to="Month_Year", values_to="seizure_counts")

all_drugs.map %>%
  arrange(seizure_counts) %>% 
  ggplot(mapping = aes(long, lat, group = group, fill=seizure_counts)) +
  geom_polygon(color = "#000000", linewidth = .05) +
  scale_fill_viridis_c(alpha=0.9, na.value="white") +
  labs(fill = "Seizure Counts", x="", y="") + 
  theme_bw() + 
  theme(legend.position="bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank()) -> seizure_counts_map
unique(ggplot_build(seizure_counts_map)$data[[1]]$fill)


all_drugs.map %>%
  mutate(seizure_counts=as.factor(seizure_counts)) %>% 
  ggplot(mapping = aes(long, lat, group = group, fill=seizure_counts)) +
  geom_polygon(color = "#000000", linewidth = .05) +
  scale_fill_manual(values=c("gray70", "#440154E6", "#440556E6", "#440757E6", "#450958E6", "#450B59E6", "#450D5AE6", "#450F5BE6", "#45115CE6", "#45135DE6",
                             "#45145EE6", "#45165FE6", "#451760E6", "#451961E6", "#461A62E6", "#461C63E6", "#461D64E6", "#461E65E6", "#462066E6", "#462167E6",
                             "#462268E6", "#462469E6", "#46256AE6", "#46266BE6", "#46276CE6", "#46296EE6", "#462A6FE6", "#452C71E6", "#452F73E6", "#453074E6",
                             "#453175E6", "#453276E6", "#453377E6", "#453478E6", "#443579E6", "#433B7EE6", "#433E82E6", "#424083E6", "#414487E6", "#404988E6",
                             "#404B88E6", "#404C88E6", "#404D88E6", "#3F5089E6", "#3F5289E6", "#3E5389E6", "#3E5489E6", "#3E5589E6", "#3D598AE6", "#3D5A8AE6",
                             "#3B5E8BE6", "#3B5F8BE6", "#37668CE6", "#36688CE6", "#336E8DE6", "#31708DE6", "#2E748DE6", "#2B778EE6", "#2B818CE6", "#4BB575E6",
                             "#AAD946E6", "#FDE725E6"),
                    na.value="white") +
  labs(fill = "Seizure Counts", x="", y="") + 
  theme_bw() + 
  theme(legend.position="bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank()) -> seizure_counts_map
# ggsave("all_drugs seizure counts map Jan_2020 (gray 0).pdf", seizure_counts_map, width=18, height=23, units="cm")


all.drugs.count <- read.csv("all_drugs count HIDTA (02-09-2023).csv") %>% as_tibble %>% 
  select(GEOID, Jan_2020:Dec_2021) %>% 
  pivot_longer(-GEOID, names_to="Month_Year", values_to="all_count")
all.drugs.weight <- read.csv("all_drugs weight HIDTA (02-09-2023).csv") %>% as_tibble %>% 
  select(GEOID, Jan_2020:Dec_2021) %>% 
  pivot_longer(-GEOID, names_to="Month_Year", values_to="all_weight")

# overdose death count map
overdose.org <- read.csv("overdose counts KNN5 R codes 999 two-sided LISA original (11-05-2023).csv") %>% as_tibble
overdose.perm.i <- read.csv("overdose counts KNN5 R codes 999 two-sided LISA permute i (11-05-2023).csv") %>% as_tibble
overdose.map <- left_join(coordinate_map[, c(1:2,6:7,10,12)], overdose.org[, c(3, 4:52)], by = "GEOID")

overdose.map <- overdose.map %>% 
  select(long:state, Jan_2020) %>% 
  pivot_longer(c(-long, -lat, -group, -GEOID, -state, -county), names_to="Month_Year", values_to="death_counts")

overdose.map %>%
  arrange(death_counts) %>% 
  mutate(death_counts=as.factor(death_counts)) %>% 
  ggplot(mapping = aes(long, lat, group = group, fill=death_counts)) +
  geom_polygon(color = "#000000", linewidth = .05) +
  scale_fill_viridis_d(alpha=0.9, na.value="white") +
  labs(fill = "Death Counts", x="", y="") + 
  theme_bw() + 
  theme(legend.position="bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank()) -> seizure_counts_map
unique(ggplot_build(seizure_counts_map)$data[[1]]$fill)


overdose.map %>%
  mutate(death_counts=as.factor(death_counts)) %>% 
  ggplot(mapping = aes(long, lat, group = group, fill=death_counts)) +
  geom_polygon(color = "#000000", linewidth = .05) +
  scale_fill_manual(values=c("gray70", "#440154E6", "#440256E6", "#450458E6", "#46065AE6", "#46085CE6", "#460A5DE6", "#460B5EE6", "#470D60E6", "#470F62E6", "#471164E6", "#471365E6", "#481467E6", "#481668E6",
                             "#481769E6", "#48196BE6", "#481B6DE6", "#481C6EE6", "#481D6FE6", "#482070E6", "#482172E6", "#482374E6", "#482475E6", "#482576E6", "#482777E6", "#482979E6", "#472A7AE6",
                             "#472C7AE6", "#472D7BE6", "#472E7CE6", "#462F7EE6", "#46327EE6", "#46337FE6", "#463480E6", "#453681E6", "#453881E6", "#443983E6", "#443A83E6", "#443B84E6", "#433D84E6",
                             "#433E85E6", "#424086E6", "#424186E6", "#414287E6", "#414487E6", "#404588E6", "#3F4788E6", "#3F4889E6", "#3E4989E6", "#3E4A89E6", "#3E4C8AE6", "#3D4E8AE6", "#3C4F8AE6",
                             "#3C508BE6", "#3B518BE6", "#3B528BE6", "#3A548CE6", "#39558CE6", "#39568CE6", "#38588CE6", "#38598CE6", "#375B8CE6", "#365C8DE6", "#365D8DE6", "#355E8DE6", "#355F8DE6",
                             "#34608DE6", "#33628DE6", "#33638DE6", "#32648EE6", "#32658EE6", "#31668EE6", "#31688EE6", "#30698EE6", "#306A8EE6", "#2F6B8EE6", "#2F6C8EE6", "#2E6E8EE6", "#2E6F8EE6",
                             "#2D708EE6", "#2D718EE6", "#2C718EE6", "#2C738EE6", "#2B748EE6", "#2B758EE6", "#2A768EE6", "#2A778EE6", "#2A788EE6", "#297A8EE6", "#297B8EE6", "#287C8EE6", "#287D8EE6",
                             "#277E8EE6", "#27808EE6", "#26818EE6", "#26828EE6", "#25838EE6", "#25858EE6", "#24868EE6", "#24878EE6", "#23888EE6", "#23898EE6", "#228B8DE6", "#228C8DE6", "#228D8DE6",
                             "#218E8DE6", "#218F8DE6", "#21908CE6", "#20928CE6", "#20938CE6", "#1F948CE6", "#1F958BE6", "#1F978BE6", "#1F988BE6", "#1F998AE6", "#1F9A8AE6", "#1E9B8AE6", "#1E9D89E6",
                             "#1F9E89E6", "#1F9F88E6", "#1FA088E6", "#1FA188E6", "#1FA187E6", "#20A386E6", "#20A486E6", "#21A585E6", "#21A685E6", "#22A885E6", "#23A983E6", "#24AA83E6", "#25AB82E6",
                             "#25AC82E6", "#27AD81E6", "#28AE80E6", "#29AF7FE6", "#2AB07FE6", "#2CB17EE6", "#2DB27CE6", "#2FB47CE6", "#31B57BE6", "#32B67AE6", "#34B679E6", "#36B779E6", "#38B977E6",
                             "#3ABA76E6", "#3BBB75E6", "#3DBC74E6", "#3FBC73E6", "#41BE71E6", "#43BF70E6", "#46C06FE6", "#48C16EE6", "#4BC16DE6", "#4DC36BE6", "#4FC46AE6", "#52C569E6", "#54C568E6",
                             "#57C666E6", "#59C864E6", "#5BC863E6", "#5EC962E6", "#60CA60E6", "#64CB5FE6", "#66CB5DE6", "#68CD5BE6", "#6BCD5AE6", "#6ECE58E6", "#71CF57E6", "#74D055E6", "#76D153E6",
                             "#79D151E6", "#7CD250E6", "#7FD34EE6", "#82D34CE6", "#85D54AE6", "#88D548E6", "#8BD646E6", "#8ED645E6", "#91D742E6", "#94D741E6", "#97D83FE6", "#9BD93CE6", "#9DD93BE6",
                             "#A1DA38E6", "#A3DA37E6", "#A7DB35E6", "#AADC32E6", "#ADDC30E6", "#B0DD2FE6", "#B3DD2CE6", "#B7DE2AE6", "#BADE28E6", "#BDDF26E6", "#C0DF25E6", "#C3DF22E6", "#C7E020E6",
                             "#C9E11FE6", "#CDE11DE6", "#D0E11CE6", "#D3E21BE6", "#D7E219E6", "#D9E319E6", "#DDE318E6", "#DFE318E6", "#E3E418E6", "#E6E419E6", "#E9E51AE6", "#ECE51BE6", "#EFE51CE6",
                             "#F2E51DE6", "#F5E61FE6", "#F7E621E6", "#FAE723E6", "#FDE725E6", "red"),
                    na.value="white") +
  labs(fill = "Death Counts", x="", y="") +
  guides(fill=guide_legend(ncol=11)) +
  theme_bw() + 
  theme(legend.position="bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank()) -> seizure_counts_map
# ggsave("overdose seizure counts map Jan_2020 (gray 0).pdf", seizure_counts_map, width=20, height=29, units="cm")


# LISA label maps
all_drugs_LISA_C.org <- all_drugs.org[,c(1:3, grep("LISA_C", names(all_drugs.org)))]
all_drugs_LISA_C.perm.i <- all_drugs.perm.i[,c(1:3, grep("LISA_C", names(all_drugs.perm.i)))]
overdose_LISA_C.org <- overdose.org[,c(1:3, grep("LISA_C", names(overdose.org)))]
overdose_LISA_C.perm.i <- overdose.perm.i[,c(1:3, grep("LISA_C", names(overdose.perm.i)))]
county.names <- all_drugs_LISA_C.org[, 1:3]

LISA_C_to_Month <- function(col_names) {
  col_names <- gsub("0", "2020", col_names)
  col_names <- gsub("1", "2021", col_names)
  col_names <- gsub("8", "2018", col_names)
  col_names <- gsub("9", "2019", col_names)
  col_names <- parse_date(col_names, "LISA_C%b%Y") %>% substr(1,7)
  return(col_names)
}

names(all_drugs_LISA_C.org)[4:51] <- LISA_C_to_Month(names(all_drugs_LISA_C.org)[4:51])
names(all_drugs_LISA_C.perm.i)[4:51] <- LISA_C_to_Month(names(all_drugs_LISA_C.perm.i)[4:51])
names(overdose_LISA_C.org)[4:51] <- LISA_C_to_Month(names(overdose_LISA_C.org)[4:51])
names(overdose_LISA_C.perm.i)[4:51] <- LISA_C_to_Month(names(overdose_LISA_C.perm.i)[4:51])

all_drugs_C.org.map <- left_join(counties.obs[, c(1:2,6:7,10,12)], all_drugs_LISA_C.org[,-(1:2)], by = "GEOID")
all_drugs_C.org.map <- all_drugs_C.org.map %>% pivot_longer(c(-long, -lat, -group, -GEOID, -state, -county), names_to="Month_Year", values_to="LISA_C")
all_drugs_C.org.map$LISA_C <- factor(all_drugs_C.org.map$LISA_C)

all_drugs_C.perm.i.map <- left_join(counties.obs[, c(1:2,6:7,10,12)], all_drugs_LISA_C.perm.i[,-(1:2)], by = "GEOID")
all_drugs_C.perm.i.map <- all_drugs_C.perm.i.map %>% pivot_longer(c(-long, -lat, -group, -GEOID, -state, -county), names_to="Month_Year", values_to="LISA_C")
all_drugs_C.perm.i.map$LISA_C <- factor(all_drugs_C.perm.i.map$LISA_C)

overdose_C.org.map <- left_join(counties.obs[, c(1:2,6:7,10,12)], overdose_LISA_C.org[,-(1:2)], by = "GEOID")
overdose_C.org.map <- overdose_C.org.map %>% pivot_longer(c(-long, -lat, -group, -GEOID, -state, -county), names_to="Month_Year", values_to="LISA_C")
overdose_C.org.map$LISA_C <- factor(overdose_C.org.map$LISA_C)

overdose_C.perm.i.map <- left_join(counties.obs[, c(1:2,6:7,10,12)], overdose_LISA_C.perm.i[,-(1:2)], by = "GEOID")
overdose_C.perm.i.map <- overdose_C.perm.i.map %>% pivot_longer(c(-long, -lat, -group, -GEOID, -state, -county), names_to="Month_Year", values_to="LISA_C")
overdose_C.perm.i.map$LISA_C <- factor(overdose_C.perm.i.map$LISA_C)

all_drugs_C.org.map %>% filter(Month_Year == "2020-01") %>% 
  ggplot(mapping = aes(long, lat, group = group, fill=LISA_C)) +
  geom_polygon(color = "#000000", linewidth = .05) +
  scale_fill_manual(values = c("Insig"="grey60",
                               "LL"="blue",
                               "LH"="steelblue",
                               "HL"="orange",
                               "HH"="red"),
                    na.value = "white") +
  labs(fill = "LISA_C_rel") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) -> LISA_C_map1
# ggsave(paste("all drugs count LISA label original Jan 2020.png", sep=""), LISA_C_map1, width=12, height=8, units="cm")

all_drugs_C.perm.i.map %>% filter(Month_Year == "2020-01") %>% 
  ggplot(mapping = aes(long, lat, group = group, fill=LISA_C)) +
  geom_polygon(color = "#000000", linewidth = .05) +
  scale_fill_manual(values = c("Insig"="grey60",
                               "LL"="blue",
                               "LH"="steelblue",
                               "HL"="orange",
                               "HH"="red"),
                    na.value = "white") +
  labs(fill = "LISA_C_rel") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) -> LISA_C_map2
# ggsave(paste("all drugs count LISA label permute i Jan 2020.png", sep=""), LISA_C_map2, width=12, height=8, units="cm")

overdose_C.org.map %>% filter(Month_Year == "2020-01") %>% 
  ggplot(mapping = aes(long, lat, group = group, fill=LISA_C)) +
  geom_polygon(color = "#000000", linewidth = .05) +
  scale_fill_manual(values = c("Insig"="grey60",
                               "LL"="blue",
                               "LH"="steelblue",
                               "HL"="orange",
                               "HH"="red"),
                    na.value = "white") +
  labs(fill = "LISA_C_rel") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) -> LISA_C_map3
# ggsave(paste("overdose count LISA label original Jan 2020.png", sep=""), LISA_C_map3, width=12, height=8, units="cm")

overdose_C.perm.i.map %>% filter(Month_Year == "2020-01") %>% 
  ggplot(mapping = aes(long, lat, group = group, fill=LISA_C)) +
  geom_polygon(color = "#000000", linewidth = .05) +
  scale_fill_manual(values = c("Insig"="grey60",
                               "LL"="blue",
                               "LH"="steelblue",
                               "HL"="orange",
                               "HH"="red"),
                    na.value = "white") +
  labs(fill = "LISA_C_rel") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) -> LISA_C_map4
# ggsave(paste("overdose count LISA label permute i Jan 2020.png", sep=""), LISA_C_map4, width=12, height=8, units="cm")