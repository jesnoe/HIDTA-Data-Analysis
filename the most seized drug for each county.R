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
  group_by(state_name, county, Drug) %>%
  summarise(count=sum(Quantity > 0))
# write.csv(seizure_counts, "overall seizure counts.csv", row.names=F)
  

seizure_count_per_drug <- seizure_counts %>% 
  group_by(state_name, county, Drug) %>% 
  summarise(count=sum(count))
# write.csv(seizure_count_per_drug, "seizure count per drug for each county during 2018-2021.csv", row.names=F)

most_seized_drug <- seizure_counts %>% 
  group_by(state_name, county) %>% 
  summarise(most_seized_drug=Drug[which.max(count)], count=count[which.max(count)])
# write.csv(most_seized_drug, "most seized drug for each county during 2018-2021.csv", row.names=F)
most_seized_drug_map <- left_join(unique(coordinate.HIDTA[,c(1:2, 6:7, 12, 14, 15)]), most_seized_drug, by=c("state_name", "county"))

most_seized_drug_map %>% ggplot(mapping = aes(long, lat, group = group, fill=most_seized_drug)) +
  geom_polygon(color = "#000000", linewidth = .05) +
  labs(fill = "Drug", title="Most Seized Drugs") + 
  scale_fill_hue(na.value="white") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) -> most_seized_drug_ggplot
# ggsave(paste("most seized drug for each county during 2018-2021.png", sep=""), most_seized_drug_ggplot, width=30, height=15, units="cm")


seizures_reduced <- seizures
seizures_reduced$Drug <- ifelse(seizures_reduced$Drug %in% c("Cocaine", "Cocaine, Powder", "Cocaine Base",
                                                                         "Cocaine precursor", "Coca Leaves", "Cocaine metabolite"),
                                      "Cocaine", seizures_reduced$Drug)

seizures_reduced$Drug <- ifelse(seizures_reduced$Drug %in% c("Methamphetamine", "Ice", "Amphetamine", "Methamphetamine in solution",
                                                                         "Meth precursor: Ephedrine hydrochloride", "Methamphetamine Oil", "Meth Precursor Chemicals",
                                                                         "Methamphetamine Tablet", "Ephedrine", "Pseudoephedrine"),
                                      "Methamphetamine", seizures_reduced$Drug)
seizures_reduced$Drug <- ifelse(grepl("Marijuana", seizures_reduced$Drug),
                                      "Marijuana", seizures_reduced$Drug)
seizures_reduced$Drug <- ifelse(grepl("THC", seizures_reduced$Drug),
                                      "Marijuana", seizures_reduced$Drug)
seizures_reduced$Drug <- ifelse(grepl("Heroin", seizures_reduced$Drug),
                                      "Heroin", seizures_reduced$Drug)
seizure_counts_reduced <- seizures_reduced %>% 
  group_by(state_name, county, Drug) %>%
  summarise(count=sum(Quantity > 0))

seizure_count_per_drug_reduced <- seizure_counts_reduced %>% 
  group_by(state_name, county, Drug) %>% 
  summarise(count=sum(count))
# write.csv(seizure_count_per_drug_reduced, "seizure count per drug (reduced) for each county during 2018-2021.csv", row.names=F)

most_seized_drug_reduced <- seizure_count_per_drug_reduced %>% 
  group_by(state_name, county) %>% 
  summarise(most_seized_drug=Drug[which.max(count)], count=count[which.max(count)])
# write.csv(most_seized_drug_reduced, "most seized drug (reduced) for each county during 2018-2021.csv", row.names=F)

most_seized_drug_reduced_map <- left_join(unique(coordinate.HIDTA[,c(1:2, 6:7, 12, 14, 15)]), most_seized_drug_reduced, by=c("state_name", "county"))
most_seized_drug_reduced_map %>% ggplot(mapping = aes(long, lat, group = group, fill=most_seized_drug)) +
  geom_polygon(color = "#000000", linewidth = .05) +
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
  geom_polygon(color = "#000000", linewidth = .05) +
  labs(fill = "Drug", title="Most Seized Drugs") + 
  scale_fill_hue(na.value="white") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) -> most_seized_drug_ggplot_reduced2
# ggsave(paste("most seized drug (reduced2) for each county during 2018-2021.png", sep=""), most_seized_drug_ggplot_reduced2, width=23, height=15, units="cm")

## Regroup by drug classifications from https://www.addictioncenter.com/drugs/drug-classifications/
drugs <- seizures$Drug %>% unique %>% sort
Barbiturates <- c("Barbiturates not specifically listed", "Phenobarbital")
Benzodiazepines <-  c("Ativan", "Benzodiazepine", "Halcion", "Klonopin", "Librium", "Valium", "Xanax", grep("gamma", drugs, value=T))
Cannabinoids <- c("Cannabis", "Hashish", grep("Marijuana", drugs, value=T), grep("THC", drugs, value=T))
Cocaine <- c("Crack", grep("Cocaine", drugs, value=T)) # Should include Coca leaves?
Meth <- c("Methamphetamine", "Ice", "Amphetamine", "Methamphetamine in solution",
          "Meth precursor: Ephedrine hydrochloride", "Methamphetamine Oil", "Meth Precursor Chemicals",
          "Methamphetamine Tablet", "Ephedrine", "Pseudoephedrine", "MDMA")
Opioids <- c("Fentanyl", grep("Heroin", drugs, value=T), "Opiates", grep("Opium", drugs, value=T), "Oxycodone")

  # by effect
Depressants <- c(Barbiturates, "GHB", "Ketamine", "Morphine", Opioids)
Stimulants <- c(Cocaine, "Dexedrine", "Ephedrine", Meth)
Hallucinogens <- c("LSD", "PCP", "Psilocybin")

seizures_class <- seizures
seizures_class$Drug <- ifelse(seizures_class$Drug %in% Barbiturates, "Barbiturates", seizures_class$Drug)
seizures_class$Drug <- ifelse(seizures_class$Drug %in% Benzodiazepines, "Benzodiazepines", seizures_class$Drug)
seizures_class$Drug <- ifelse(seizures_class$Drug %in% Cannabinoids, "Cannabinoids", seizures_class$Drug)
seizures_class$Drug <- ifelse(seizures_class$Drug %in% Depressants, "Depressants", seizures_class$Drug)
seizures_class$Drug <- ifelse(seizures_class$Drug %in% Stimulants, "Stimulants", seizures_class$Drug)
seizures_class$Drug <- ifelse(seizures_class$Drug %in% Hallucinogens, "Hallucinogens", seizures_class$Drug)

seizure_counts_class <- seizures_class %>% 
  group_by(state_name, county, Drug) %>%
  summarise(count=sum(Quantity > 0))

most_seized_drug_class <- seizure_counts_class %>% 
  group_by(state_name, county) %>% 
  summarise(Stimulants_prop=ifelse("Stimulants" %in% Drug, count[Drug == "Stimulants"]/sum(count), 0.001), 
            most_seized_drug=Drug[which.max(count)],
            count=max(count))
most_seized_drug_class$most_seized_drug %>% unique %>% sort

most_seized_drug_class %>% select(-Stimulants_prop) %>% 
  filter(!(most_seized_drug %in% c("Cannabinoids", "Depressants", "Stimulants", "Hallucinogens"))) %>% arrange(desc(count))

most_seized_drug_class2 <- most_seized_drug_class
most_seized_drug_class2$most_seized_drug <- ifelse(most_seized_drug_class2$most_seized_drug %in% c("Cannabinoids", "Depressants", "Stimulants", "Hallucinogens"),
                                                   most_seized_drug_class2$most_seized_drug, "Other")

most_seized_drug_class_map <- left_join(unique(coordinate.HIDTA[,c(1:2, 6:7, 12, 14, 15)]), most_seized_drug_class2, by=c("state_name", "county"))
most_seized_drug_class_map %>% ggplot(mapping = aes(long, lat, group = group, fill=most_seized_drug)) +
  geom_polygon(color = "#000000", linewidth = .05) +
  labs(fill = "Drug", title="Most Seized Drugs") + 
  scale_fill_manual(values= c("Cannabinoids"="green",
                              "Depressants"="blue",
                              "Stimulants"="red",
                              "Hallucinogens"="orange",
                              "Other"="purple"),
                    na.value="white") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) -> most_seized_drug_ggplot_class
# ggsave(paste("most seized drug (classified) for each county during 2018-2021.png", sep=""), most_seized_drug_ggplot_class, width=22, height=15, units="cm")

most_seized_drug_class_map %>% ggplot(mapping = aes(long, lat, group = group, fill=most_seized_drug, alpha=Stimulants_prop)) +
  geom_polygon(color = "#000000", linewidth = .05) +
  labs(fill = "Drug", title="Most Seized Drugs") + 
  scale_fill_manual(values= c("Cannabinoids"="green",
                              "Depressants"="blue",
                              "Stimulants"="red",
                              "Hallucinogens"="orange",
                              "Other"="purple"),
                    na.value="white") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) -> most_seized_drug_ggplot_class_alpha
# ggsave(paste("most seized drug (classified alpha) for each county during 2018-2021.png", sep=""), most_seized_drug_ggplot_class_alpha, width=22, height=15, units="cm")

second_seized_drug_class <- seizure_counts_class %>% 
  group_by(state_name, county) %>% 
  arrange(count, .by_group=T) %>% 
  summarise(second_seized_drug=Drug[ifelse(n()-1 > 0, n()-1, NA)],
            count=ifelse(n()-1 > 0, count[n()-1], NA))
second_seized_drug_class$second_seized_drug <- ifelse(second_seized_drug_class$second_seized_drug %in% c("Cannabinoids", "Depressants", "Stimulants", "Hallucinogens"),
                                                       second_seized_drug_class$second_seized_drug, "Other")

second_seized_drug_class_map <- left_join(unique(coordinate.HIDTA[,c(1:2, 6:7, 12, 14, 15)]), second_seized_drug_class, by=c("state_name", "county"))
second_seized_drug_class_map %>% ggplot(mapping = aes(long, lat, group = group, fill=second_seized_drug)) +
  geom_polygon(color = "#000000", linewidth = .05) +
  labs(fill = "Drug", title="Second-Most Seized Drugs") + 
  scale_fill_manual(values= c("Cannabinoids"="green",
                              "Depressants"="blue",
                              "Stimulants"="red",
                              "Hallucinogens"="orange",
                              "Other"="purple"),
                    na.value="white") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) -> second_seized_drug_ggplot_class
# ggsave(paste("second seized drug (classified) for each county during 2018-2021.png", sep=""), second_seized_drug_ggplot_class, width=22, height=15, units="cm")

 # Without Cannabinoids
most_seized_drug_class3 <- seizure_counts_class %>% 
  filter(Drug != "Cannabinoids") %>% 
  group_by(state_name, county) %>% 
  summarise(Stimulants_prop=ifelse("Stimulants" %in% Drug, count[Drug == "Stimulants"]/sum(count), 0.001), 
            most_seized_drug=Drug[which.max(count)],
            count=max(count))

most_seized_drug_class3$most_seized_drug <- ifelse(most_seized_drug_class3$most_seized_drug %in% c("Depressants", "Stimulants", "Hallucinogens"),
                                                   most_seized_drug_class3$most_seized_drug, "Other")

most_seized_drug_class_map3 <- left_join(unique(coordinate.HIDTA[,c(1:2, 6:7, 12, 14, 15)]), most_seized_drug_class3, by=c("state_name", "county"))
most_seized_drug_class_map3 %>% ggplot(mapping = aes(long, lat, group = group, fill=most_seized_drug)) +
  geom_polygon(color = "#000000", linewidth = .05) +
  labs(fill = "Drug", title="Most Seized Drugs") + 
  scale_fill_manual(values= c("Depressants"="blue",
                              "Stimulants"="red",
                              "Hallucinogens"="orange",
                              "Other"="purple"),
                    na.value="white") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) -> most_seized_drug_ggplot_class
# ggsave(paste("most seized drug (classified, without Cannabinoids) for each county during 2018-2021.png", sep=""), most_seized_drug_ggplot_class, width=22, height=15, units="cm")

most_seized_drug_class_map3 %>% ggplot(mapping = aes(long, lat, group = group, fill=most_seized_drug, alpha=Stimulants_prop)) +
  geom_polygon(color = "#000000", linewidth = .05) +
  labs(fill = "Drug", title="Most Seized Drugs") + 
  scale_fill_manual(values= c("Depressants"="blue",
                    "Stimulants"="red",
                    "Hallucinogens"="orange",
                    "Other"="purple"),
                    na.value="white") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) -> most_seized_drug_ggplot_class_alpha
# ggsave(paste("most seized drug (classified alpha, without Cannabinoids) for each county during 2018-2021.png", sep=""), most_seized_drug_ggplot_class_alpha, width=22, height=15, units="cm")

second_seized_drug_class2 <- seizure_counts_class %>% 
  filter(Drug != "Cannabinoids") %>% 
  group_by(state_name, county) %>% 
  arrange(count, .by_group=T) %>% 
  summarise(second_seized_drug=Drug[ifelse(n()-1 > 0, n()-1, NA)],
            count=ifelse(n()-1 > 0, count[n()-1], NA))
second_seized_drug_class2$second_seized_drug <- ifelse(second_seized_drug_class2$second_seized_drug %in% c("Depressants", "Stimulants", "Hallucinogens"),
                                                      second_seized_drug_class2$second_seized_drug, "Other")

second_seized_drug_class_map2 <- left_join(unique(coordinate.HIDTA[,c(1:2, 6:7, 12, 14, 15)]), second_seized_drug_class2, by=c("state_name", "county"))
second_seized_drug_class_map2 %>% ggplot(mapping = aes(long, lat, group = group, fill=second_seized_drug)) +
  geom_polygon(color = "#000000", linewidth = .05) +
  labs(fill = "Drug", title="Second-Most Seized Drugs") + 
  scale_fill_manual(values= c("Depressants"="blue",
                              "Stimulants"="red",
                              "Hallucinogens"="orange",
                              "Other"="purple"),
                    na.value="white") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) -> second_seized_drug_ggplot_class
# ggsave(paste("second seized drug (classified, without Cannabinoids) for each county during 2018-2021.png", sep=""), second_seized_drug_ggplot_class, width=22, height=15, units="cm")

## Correlation Cannabinoids vs. Crack / Cocaine / Meth
seizure_counts_for_cor <- seizure_counts
seizure_counts_for_cor$Drug <- ifelse(seizure_counts_for_cor$Drug %in% Cannabinoids, "Cannabinoids", seizure_counts_for_cor$Drug)
seizure_counts_for_cor$Drug <- ifelse(seizure_counts_for_cor$Drug %in% grep("Cocaine", drugs, value=T), "Other_Cocaine", seizure_counts_for_cor$Drug)
seizure_counts_for_cor$Drug <- ifelse(seizure_counts_for_cor$Drug %in% Meth, "Meth", seizure_counts_for_cor$Drug)
seizure_counts_for_cor <- seizure_counts_for_cor %>% 
  filter(Drug %in% c("Cannabinoids", "Crack", "Other_Cocaine", "Meth")) %>% 
  group_by(state_name, county, Drug) %>%
  summarise(count=sum(count)) %>% 
  pivot_wider(names_from="Drug", values_from="count") %>% 
  replace_na(list(Cannabinoids=0, Crack=0, Other_Cocaine=0, Meth=0))

cor(seizure_counts_for_cor[, 3:6])

## Meth vs. Ice
seizure_counts_meth <- seizure_counts
seizure_counts_meth$Drug <- ifelse(seizure_counts_meth$Drug %in% Meth[Meth != "Ice"], "Meth", seizure_counts_meth$Drug)
seizure_counts_meth <- seizure_counts_meth %>% 
  filter(Drug %in% c("Meth", "Ice")) %>% 
  group_by(state_name, county, Drug) %>%
  summarise(count=sum(count)) %>% 
  pivot_wider(names_from="Drug", values_from="count") %>% 
  replace_na(list(Meth=0, Ice=0))
seizure_counts_meth

cor(seizure_counts_meth[,3:4])

seizure_counts_meth[,3:4] %>%
  group_by(Ice, Meth) %>% 
  summarise(num_of_cases=n()) %>%
  ggplot(aes(x=Meth, y=Ice, size=num_of_cases)) +
  geom_point() +
  labs(xlab="Meth Seizure Count", ylab="Ice Seizure Count")

seizure_counts_meth %>% arrange(desc(Ice))
seizure_counts_meth %>% arrange(desc(Meth))
seizure_counts_meth_map <- left_join(unique(coordinate.HIDTA[,c(1:2, 6:7, 12, 14, 15)]), seizure_counts_meth, by=c("state_name", "county"))
seizure_counts_meth_map %>% ggplot(mapping = aes(long, lat, group = group, fill=Meth)) +
  geom_polygon(color = "#000000", linewidth = .05) +
  labs(fill = "Drug", title="Meth Seizure Counts") + 
  scale_fill_viridis_c(na.value="white") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

seizure_counts_meth_map %>% ggplot(mapping = aes(long, lat, group = group, fill=Ice)) +
  geom_polygon(color = "#000000", linewidth = .05) +
  labs(fill = "Drug", title="Ice Seizure Counts") + 
  scale_fill_viridis_c(na.value="white") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

seizure_counts_meth_map2 <- left_join(unique(coordinate.HIDTA[,c(1:2, 6:7, 12, 14, 15)]),
                                     seizure_counts_meth %>% filter(Ice < 3000),
                                     by=c("state_name", "county"))
seizure_counts_meth_map2 %>% ggplot(mapping = aes(long, lat, group = group, fill=Ice)) +
  geom_polygon(color = "#000000", linewidth = .05) +
  labs(fill = "Drug", title="Ice Seizure Counts") + 
  scale_fill_viridis_c(na.value="white") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

seizure_counts_meth %>% arrange(desc(Meth))
seizure_counts_meth %>% arrange(desc(Ice))
