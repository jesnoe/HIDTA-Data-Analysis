library(fpp2)
library(spdep)
library(readxl)
library(scales)
library(stringr)
library(urbnmapr)
library(tidyverse)
library(gridExtra)
library(lubridate)

counties.obs <- counties
names(counties.obs)[7] <- "GEOID"
counties.obs <- counties.obs %>% rename(state=state_name, county=county_name)
counties.obs$GEOID <- as.numeric(counties.obs$GEOID)
counties.obs <- counties.obs %>% filter(!(state %in% c("Alaska", "Hawaii")))

crack <- read.csv("crack counts KNN5 R codes 999 two-sided relative (02-21-2023).csv") %>% as_tibble
cocaine <- read.csv("cocaine other counts KNN5 R codes 999 two-sided relative (02-21-2023).csv") %>% as_tibble
marijuana <- read.csv("marijuana counts KNN5 R codes 999 two-sided relative (02-21-2023).csv") %>% as_tibble
meth <- read.csv("meth counts KNN5 R codes 999 two-sided relative (02-21-2023).csv") %>% as_tibble

crack_LISA_C <- crack[,c(1:3, grep("LISA_C", names(crack)))]
cocaine_LISA_C <- cocaine[,c(1:3, grep("LISA_C", names(cocaine)))]
marijuana_LISA_C <- marijuana[,c(1:3, grep("LISA_C", names(marijuana)))]
meth_LISA_C <- meth[,c(1:3, grep("LISA_C", names(meth)))]
county.names <- crack[, 1:3]

LISA_C_to_Month <- function(col_names) {
  col_names <- gsub("0", "2020", col_names)
  col_names <- gsub("1", "2021", col_names)
  col_names <- gsub("8", "2018", col_names)
  col_names <- gsub("9", "2019", col_names)
  col_names <- parse_date(col_names, "LISA_C%b%Y") %>% substr(1,7)
  return(col_names)
}

names(crack_LISA_C)[4:51] <- LISA_C_to_Month(names(crack_LISA_C)[4:51])
names(cocaine_LISA_C)[4:51] <- LISA_C_to_Month(names(cocaine_LISA_C)[4:51])
names(marijuana_LISA_C)[4:51] <- LISA_C_to_Month(names(marijuana_LISA_C)[4:51])
names(meth_LISA_C)[4:51] <- LISA_C_to_Month(names(meth_LISA_C)[4:51])

crack_counts <- crack %>% select(Jan_2018:Dec_2021)
cocaine_counts <- cocaine %>% select(Jan_2018:Dec_2021)
marijuana_counts <- marijuana %>% select(Jan_2018:Dec_2021)
meth_counts <- meth %>% select(Jan_2018:Dec_2021)

seizure_count_monthly <- crack %>% select(state:Dec_2021)
seizure_count_monthly[,-(1:4)] <- crack_counts + cocaine_counts + marijuana_counts + meth_counts
seizure_count_monthly <- seizure_count_monthly %>%
  pivot_longer(c(-state, -county, -GEOID, -HIDTA), names_to="Month_Year", values_to="seizure_count")

seizure_count_monthly_map <- left_join(counties.obs[, c(1:2,6:7,10,12)], seizure_count_monthly[,-(1:2)], by = "GEOID")
seizure_count_monthly_map$Month_Year <- parse_date(seizure_count_monthly_map$Month_Year, "%b_%Y")

for (year in 2018:2021) { # Monthly Seizure Count Out of 4 Drugs
  seizure_count_monthly_map %>% filter(year(Month_Year) == year & month(Month_Year) %in% 1:4) %>% 
    ggplot(mapping = aes(long, lat, group = group, fill=seizure_count)) +
    geom_polygon(color = "#000000", size = .05) +
    facet_wrap(. ~ substr(Month_Year, 1, 7)) +
    labs(fill = "Seizure Count", title="Monthly Seizure Count Out of 4 Drugs") + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) -> LISA_C_map
  ggsave(paste("seizure_count_monthly_map_", year, " (4 months).jpg", sep=""), LISA_C_map, width=20, height=15, units="cm")
}

# Number of LISA_C maps (Monthly)
n_LISA_C_monthly <- function(drug_LISA_C, LISA_C) {
  drug_n_LISA_C_monthly <- drug_LISA_C[,-(1:3)] == LISA_C
  return(drug_n_LISA_C_monthly)
}

LISA_class <- 1
n_HH_monthly <- cbind(county.names,
              n_LISA_C_monthly(crack_LISA_C, LISA_class) +
                n_LISA_C_monthly(cocaine_LISA_C, LISA_class) +
                n_LISA_C_monthly(marijuana_LISA_C, LISA_class) +
                n_LISA_C_monthly(meth_LISA_C, LISA_class)) %>% 
  as_tibble

LISA_class <- 2
n_LL_monthly <- cbind(county.names,
                      n_LISA_C_monthly(crack_LISA_C, LISA_class) +
                        n_LISA_C_monthly(cocaine_LISA_C, LISA_class) +
                        n_LISA_C_monthly(marijuana_LISA_C, LISA_class) +
                        n_LISA_C_monthly(meth_LISA_C, LISA_class)) %>% 
  as_tibble

LISA_class <- 3
n_LH_monthly <- cbind(county.names,
                      n_LISA_C_monthly(crack_LISA_C, LISA_class) +
                        n_LISA_C_monthly(cocaine_LISA_C, LISA_class) +
                        n_LISA_C_monthly(marijuana_LISA_C, LISA_class) +
                        n_LISA_C_monthly(meth_LISA_C, LISA_class)) %>% 
  as_tibble

LISA_class <- 4
n_HL_monthly <- cbind(county.names,
                      n_LISA_C_monthly(crack_LISA_C, LISA_class) +
                        n_LISA_C_monthly(cocaine_LISA_C, LISA_class) +
                        n_LISA_C_monthly(marijuana_LISA_C, LISA_class) +
                        n_LISA_C_monthly(meth_LISA_C, LISA_class)) %>% 
  as_tibble


n_HH_monthly_map <- left_join(counties.obs[, c(1:2,6:7,10,12)], n_HH_monthly[,-(1:2)], by = "GEOID")
n_HH_monthly_map <- n_HH_monthly_map %>% pivot_longer(c(-long, -lat, -GEOID, -group, -state, -county), names_to="Month_Year", values_to="n_HH")
n_HH_monthly_map$Month_Year <- parse_date(n_HH_monthly_map$Month_Year, "%Y-%m")

for (year in 2018:2021) { # monthly Number of HH Out of 4 Drugs
  n_HH_monthly_map %>% filter(year(Month_Year) == year & month(Month_Year) %in% 1:4) %>% 
    ggplot(mapping = aes(long, lat, group = group, fill=n_HH)) +
    geom_polygon(color = "#000000", size = .05) +
    facet_wrap(. ~ substr(Month_Year, 1, 7)) +
    labs(fill = "# of HH", title="Monthly Number of HH Out of 4 Drugs") + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) -> LISA_C_map
  ggsave(paste("n_HH_monthly_map_", year, " (4 months).jpg", sep=""), LISA_C_map, width=20, height=15, units="cm")
}

n_HH_annual <- n_HH_monthly %>% pivot_longer(c(-GEOID, -state, -county), names_to="Month_Year", values_to="n_HH")
n_HH_annual$Year <- year(parse_date(n_HH_annual$Month_Year, format="%Y-%m"))
n_HH_annual <- n_HH_annual %>% group_by(GEOID, Year) %>% summarise(n_HH=sum(n_HH))
n_HH_annual_map <- left_join(counties.obs[, c(1:2,6:7,10,12)], n_HH_annual, by = "GEOID") %>% filter(!is.na(Year))
n_HH_annual_map
n_HH_annual_map %>% 
  ggplot(mapping = aes(long, lat, group = group, fill=n_HH)) +
  geom_polygon(color = "#000000", size = .05) +
  facet_wrap(. ~ substr(Year, 1, 4)) +
  labs(fill = "# of HH", title="Annual Number of HH Out of 4 Drugs") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) -> LISA_C_map
# ggsave(paste("n_HH_annual_map.jpg", sep=""), LISA_C_map, width=20, height=15, units="cm")

n_HL_monthly_map <- left_join(counties.obs[, c(1:2,6:7,10,12)], n_HL_monthly[,-(1:2)], by = "GEOID")
n_HL_monthly_map <- n_HL_monthly_map %>% pivot_longer(c(-long, -lat, -GEOID, -group, -state, -county), names_to="Month_Year", values_to="n_HH")
n_HL_monthly_map$Month_Year <- parse_date(n_HL_monthly_map$Month_Year, "%Y-%m")

for (year in 2018:2021) { # monthly Number of HL Out of 4 Drugs
  n_HL_monthly_map %>% filter(year(Month_Year) == year & month(Month_Year) %in% 1:4) %>% 
    ggplot(mapping = aes(long, lat, group = group, fill=n_HH)) +
    geom_polygon(color = "#000000", size = .05) +
    facet_wrap(. ~ substr(Month_Year, 1, 7)) +
    labs(fill = "# of HL", title="Monthly Number of HL Out of 4 Drugs") + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) -> LISA_C_map
  ggsave(paste("n_HL_monthly_map_", year, " (4 months).jpg", sep=""), LISA_C_map, width=20, height=15, units="cm")
}

n_LH_monthly_map <- left_join(counties.obs[, c(1:2,6:7,10,12)], n_LH_monthly[,-(1:2)], by = "GEOID")
n_LH_monthly_map <- n_LH_monthly_map %>% pivot_longer(c(-long, -lat, -GEOID, -group, -state, -county), names_to="Month_Year", values_to="n_LH")
n_LH_monthly_map$Month_Year <- parse_date(n_LH_monthly_map$Month_Year, "%Y-%m")

for (year in 2018:2021) { # monthly Number of LH Out of 4 Drugs
  n_LH_monthly_map %>% filter(year(Month_Year) == year & month(Month_Year) %in% 1:4) %>% 
    ggplot(mapping = aes(long, lat, group = group, fill=n_LH)) +
    geom_polygon(color = "#000000", size = .05) +
    facet_wrap(. ~ substr(Month_Year, 1, 7)) +
    labs(fill = "# of LH", title="Monthly Number of LH Out of 4 Drugs") + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) -> LISA_C_map
  ggsave(paste("n_LH_monthly_map_", year, " (4 months).jpg", sep=""), LISA_C_map, width=20, height=15, units="cm")
}

n_HH_monthly$total_n_HH <- n_HH_monthly %>% select(`2018-01`:`2021-12`) %>% apply(1, sum)
n_HH_total_map <- left_join(counties.obs[, c(1:2,6:7,10,12)], n_HH_monthly %>% select(GEOID, total_n_HH), by = "GEOID")
n_HH_total_map %>% ggplot(mapping = aes(long, lat, group = group, fill=total_n_HH)) +
  geom_polygon(color = "#000000", size = .05) +
  labs(fill = "# of HH", title="Total Number of HH Out of 4 Drugs") + 
  scale_fill_viridis_c(na.value="white") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) -> n_HH_total_ggplot
ggsave(paste("n_HH_total_map_", year, " (4 months).jpg", sep=""), n_HH_total_ggplot, width=20, height=15, units="cm")

n_HL_monthly$total_n_HL <- n_HL_monthly %>% select(`2018-01`:`2021-12`) %>% apply(1, sum)
n_HL_total_map <- left_join(counties.obs[, c(1:2,6:7,10,12)], n_HL_monthly %>% select(GEOID, total_n_HL), by = "GEOID")
n_HL_total_map %>% ggplot(mapping = aes(long, lat, group = group, fill=total_n_HL)) +
  geom_polygon(color = "#000000", size = .05) +
  labs(fill = "# of HL", title="Total Number of HL Out of 4 Drugs") + 
  scale_fill_viridis_c(na.value="white") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) -> n_HL_total_ggplot
ggsave(paste("n_HL_total_map_", year, " (4 months).jpg", sep=""), n_HL_total_ggplot, width=20, height=15, units="cm")

n_LH_monthly$total_n_LH <- n_LH_monthly %>% select(`2018-01`:`2021-12`) %>% apply(1, sum)
n_LH_total_map <- left_join(counties.obs[, c(1:2,6:7,10,12)], n_LH_monthly %>% select(GEOID, total_n_LH), by = "GEOID")
n_LH_total_map %>% ggplot(mapping = aes(long, lat, group = group, fill=total_n_LH)) +
  geom_polygon(color = "#000000", size = .05) +
  labs(fill = "# of LH", title="Total Number of LH Out of 4 Drugs") + 
  scale_fill_viridis_c(na.value="white") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) -> n_LH_total_ggplot
ggsave(paste("n_LH_total_map_", year, " (4 months).jpg", sep=""), n_LH_total_ggplot, width=20, height=15, units="cm")

# Moran's I plot
crack_LISA_I <- crack[,c(1:3, grep("LISA_I", names(crack)))]
LISA_I_to_Month <- function(col_names) {
  col_names <- gsub("0", "2020", col_names)
  col_names <- gsub("1", "2021", col_names)
  col_names <- gsub("8", "2018", col_names)
  col_names <- gsub("9", "2019", col_names)
  col_names <- parse_date(col_names, "LISA_I%b%Y") %>% substr(1,7)
  return(col_names)
}
names(crack_LISA_I)[4:51] <- LISA_I_to_Month(names(crack_LISA_I)[4:51])
crack_LISA_I_map <- left_join(counties.obs[, c(1:2,6:7,10,12)], crack_LISA_I[,-(1:2)], by = "GEOID")
crack_LISA_I_map <- crack_LISA_I_map %>% pivot_longer(c(-long, -lat, -GEOID, -group, -state, -county), names_to="Month_Year", values_to="LISA_I")
crack_LISA_I_map$Month_Year <- parse_date(crack_LISA_I_map$Month_Year, "%Y-%m")

for (year in 2018:2021) { # monthly Moran's I of Crack
  crack_LISA_I_map %>% filter(year(Month_Year) == year & month(Month_Year) %in% 1:4) %>% 
    ggplot(mapping = aes(long, lat, group = group, fill=LISA_I)) +
    geom_polygon(color = "#000000", size = .05) +
    facet_wrap(. ~ substr(Month_Year, 1, 7)) +
    labs(fill = "LISA_I", title="Monthly Moran's I of Crack") + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) -> LISA_C_map
  ggsave(paste("crack_LISA_I_map", year, " (4 months).png", sep=""), LISA_C_map, width=20, height=15, units="cm")
}
