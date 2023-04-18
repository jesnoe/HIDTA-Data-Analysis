library(fpp2)
library(spdep)
library(readxl)
library(urbnmapr)
library(tidyverse)
library(gridExtra)
library(lubridate)

# urbnmapr data
unique(countydata$county_fips) # some ips have 5 digits. need to be converted into integer.
HIDTA.dist <- read.csv("HIDTA Regions.csv", header=T) %>% as_tibble
coordinate.HIDTA <- left_join(counties, HIDTA.dist, by=c("state_name", "county_name")) %>% filter(!(state_name %in% c("Alaska", "Hawaii", "Puerto Rico")))
names(coordinate.HIDTA)[7] <- "GEOID"
coordinate.HIDTA$county <- str_split(coordinate.HIDTA$county_name, " ") %>%
  lapply(function(x) str_c(c(x[(1:(length(x)-1))]), collapse=" ") ) %>%
  unlist
HIDTA.regions <- coordinate.HIDTA %>% filter(!is.na(HIDTA)) %>% select(GEOID, HIDTA) %>% unique

seizures <- read_xlsx("Drug Seizures All HIDTAs All Drugs 2018-2021 Combined.xlsx") %>% filter(State %in% coordinate.HIDTA$state_name)
seizures <- seizures %>% filter(County != "Michigan")
seizures$County <- substring(seizures$County, 1, str_locate(seizures$County, ",")[,1]-1)
# seizures$County <- sapply(seizures$County, function(x) paste(x[1], str_to_lower(substring(x, 2, length(x)-1)), sep=""))
seizures$Year <- substring(seizures$SeizureDate, 1, 4)
seizures$Month <- substring(seizures$SeizureDate, 6, 7)
seizures$Year <- as.numeric(seizures$Year)
seizures$Month <- as.numeric(seizures$Month)
names(seizures)[1:2] <- c("state_name", "county")
seizure_counts <- seizures %>% 
  filter(Unit=="Kg") %>% 
  group_by(state_name, county, Drug) %>%
  summarise(count=sum(Quantity > 0))
# write.csv(seizure_counts, "overall seizure counts.csv", row.names=F)
  

most_seized_drug <- seizure_counts %>% 
  group_by(state_name, county) %>% 
  summarise(most_seized_drug=Drug[which.max(count)], count=count[which.max(count)])
# write.csv(most_seized_drug, "most seized drug for each county during 2018-2021.csv", row.names=F)
most_seized_drug_map <- left_join(unique(coordinate.HIDTA[,c(1:2, 6:7, 12, 14, 15)]), most_seized_drug, by=c("state_name", "county"))

most_seized_drug_map %>% ggplot(mapping = aes(long, lat, group = group, fill=most_seized_drug)) +
  geom_polygon(color = "#000000", size = .05) +
  labs(fill = "Drug", title="Most Seized Drugs") + 
  scale_fill_hue(na.value="white") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) -> most_seized_drug_ggplot
# ggsave(paste("most seized drug for each county during 2018-2021.png", sep=""), most_seized_drug_ggplot, width=30, height=15, units="cm")


seizure_counts_reduced <- seizure_counts
seizure_counts_reduced$Drug <- ifelse(seizure_counts_reduced$Drug %in% c("Cocaine", "Cocaine, Powder", "Cocaine Base",
                                                                         "Cocaine precursor", "Coca Leaves", "Cocaine metabolite"),
                                      "Cocaine", seizure_counts_reduced$Drug)

seizure_counts_reduced$Drug <- ifelse(seizure_counts_reduced$Drug %in% c("Methamphetamine", "Ice", "Amphetamine", "Methamphetamine in solution",
                                                                         "Meth precursor: Ephedrine hydrochloride", "Methamphetamine Oil", "Meth Precursor Chemicals",
                                                                         "Methamphetamine Tablet", "Ephedrine", "Pseudoephedrine"),
                                      "Methamphetamine", seizure_counts_reduced$Drug)
seizure_counts_reduced$Drug <- ifelse(grepl("Marijuana", seizure_counts_reduced$Drug),
                                      "Marijuana", seizure_counts_reduced$Drug)
seizure_counts_reduced$Drug <- ifelse(grepl("Heroin", seizure_counts_reduced$Drug),
                                      "Heroin", seizure_counts_reduced$Drug)
seizure_counts_reduced$Drug <- ifelse(grepl("THC", seizure_counts_reduced$Drug),
                                      "THC (Tetrahydrocannabinols)", seizure_counts_reduced$Drug)
most_seized_drug_reduced <- seizure_counts_reduced %>% 
  group_by(state_name, county) %>% 
  summarise(most_seized_drug=Drug[which.max(count)], count=count[which.max(count)])
# write.csv(most_seized_drug_reduced, "most seized drug (reduced) for each county during 2018-2021.csv", row.names=F)
most_seized_drug_reduced_map <- left_join(unique(coordinate.HIDTA[,c(1:2, 6:7, 12, 14, 15)]), most_seized_drug_reduced, by=c("state_name", "county"))

most_seized_drug_reduced_map %>% ggplot(mapping = aes(long, lat, group = group, fill=most_seized_drug)) +
  geom_polygon(color = "#000000", size = .05) +
  labs(fill = "Drug", title="Most Seized Drugs") + 
  scale_fill_hue(na.value="white") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) -> most_seized_drug_ggplot_reduced
# ggsave(paste("most seized drug (reduced) for each county during 2018-2021.png", sep=""), most_seized_drug_ggplot_reduced, width=28, height=15, units="cm")

most_seized_drug_reduced2 <- most_seized_drug_reduced %>% 
  filter(most_seized_drug %in% c("Cocaine", "Crack", "Fentanyl", "Heroin", "Marijuana", "Methamphetamine"))
most_seized_drug_reduced_map2 <- left_join(unique(coordinate.HIDTA[,c(1:2, 6:7, 12, 14, 15)]), most_seized_drug_reduced2, by=c("state_name", "county"))
most_seized_drug_reduced_map2 %>%
  ggplot(mapping = aes(long, lat, group = group, fill=most_seized_drug)) +
  geom_polygon(color = "#000000", size = .05) +
  labs(fill = "Drug", title="Most Seized Drugs") + 
  scale_fill_hue(na.value="white") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) -> most_seized_drug_ggplot_reduced2
# ggsave(paste("most seized drug (reduced2) for each county during 2018-2021.png", sep=""), most_seized_drug_ggplot_reduced2, width=23, height=15, units="cm")
