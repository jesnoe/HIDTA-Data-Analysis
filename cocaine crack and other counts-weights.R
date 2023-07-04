# setwd("/Users/euseongjang/Documents/R")
# setwd("C:/Users/gkfrj/Documents/R")
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

# all cocaine count
{
all_cocaine <- seizures %>% filter(Drug %in% c("Cocaine", "Cocaine, Powder", "Cocaine Base", "Crack",
                                               "Cocaine precursor", "Coca Leaves", "Cocaine metabolite"))
all_cocaine <- all_cocaine %>% mutate(Month2=month(Month, label=T))
calender <- data.frame(Month=1:12, Month2=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
all_cocaine <- all_cocaine %>% group_by(state_name, county, Year, Month2) %>%
  summarise(count=sum(Quantity > 0), .groups="drop") %>%
  complete(state_name, county, Year, Month2) %>% 
  pivot_wider(names_from=c(Month2, Year), names_sep="_", values_from=count)
cal_order <- c(str_c(calender$Month2, 2018, sep="_"),
               str_c(calender$Month2, 2019, sep="_"),
               str_c(calender$Month2, 2020, sep="_"),
               str_c(calender$Month2, 2021, sep="_"))
all_cocaine <- all_cocaine[, c(1:2, match(cal_order, names(all_cocaine)))]

all_cocaine <- left_join(unique(coordinate.HIDTA[,c(7, 12, 14, 15)]), all_cocaine, by=c("state_name", "county")) %>%
  select(state_name, county, GEOID, HIDTA, Jan_2018:Dec_2021) %>% arrange(GEOID)
names(all_cocaine)[1] <- "state"


all_cocaine %>% filter(is.na(county))
grep("Jan_2018", names(all_cocaine)) # 5
grep("Dec_2021", names(all_cocaine)) # 52
all_cocaine <- all_cocaine %>% filter(GEOID %in% non_missing_GEOID | !is.na(HIDTA))
all_cocaine[,5:52][is.na(all_cocaine[,5:52])] <- 0
all_cocaine %>% filter(is.na(GEOID))
}
all_cocaine # 1518 counties in total.
# write.csv(all_cocaine, "all_cocaine count HIDTA (06-21-2023).csv", row.names=F)

# crack count
{
crack <- seizures %>% filter(Drug=="Crack")
crack <- crack %>% mutate(Month2=month(Month, label=T))
calender <- data.frame(Month=1:12, Month2=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
crack <- crack %>% group_by(state_name, county, Year, Month2) %>%
  summarise(count=sum(Quantity > 0), .groups="drop") %>%
  complete(state_name, county, Year, Month2) %>% 
  pivot_wider(names_from=c(Month2, Year), names_sep="_", values_from=count)
cal_order <- c(str_c(calender$Month2, 2018, sep="_"),
               str_c(calender$Month2, 2019, sep="_"),
               str_c(calender$Month2, 2020, sep="_"),
               str_c(calender$Month2, 2021, sep="_"))
crack <- crack[, c(1:2, match(cal_order, names(crack)))]

crack <- left_join(unique(coordinate.HIDTA[,c(7, 12, 14, 15)]), crack, by=c("state_name", "county")) %>%
  select(state_name, county, GEOID, HIDTA, Jan_2018:Dec_2021) %>% arrange(GEOID)
names(crack)[1] <- "state"
}

crack %>% filter(is.na(county))
grep("Jan_2018", names(crack)) # 5
grep("Dec_2021", names(crack)) # 52
crack <- crack %>% filter(GEOID %in% non_missing_GEOID | !is.na(HIDTA))
crack[,5:52][is.na(crack[,5:52])] <- 0
crack %>% filter(is.na(GEOID))
crack # 1518 counties in total.
# write.csv(crack, "cocaine crack count HIDTA (06-28-2023).csv", row.names=F)

# other cocaine count
{
  cocaine <- seizures %>% filter(Drug %in% c("Cocaine", "Cocaine, Powder", "Cocaine Base",
                                             "Cocaine precursor", "Coca Leaves", "Cocaine metabolite"), Unit == "Kg")
  cocaine <- cocaine %>% mutate(Month2=month(Month, label=T))
  cocaine <- cocaine %>% group_by(state_name, county, Year, Month2) %>%
    summarise(count=sum(Quantity > 0), .groups="drop") %>%
    complete(state_name, county, Year, Month2) %>% 
    pivot_wider(names_from=c(Month2, Year), names_sep="_", values_from=count)
  cocaine <- cocaine[, c(1:2, match(cal_order, names(cocaine)))]
  
  cocaine <- left_join(unique(coordinate.HIDTA[,c(7, 12, 14, 15)]), cocaine, by=c("state_name", "county")) %>%
    select(state_name, county, GEOID, HIDTA, Jan_2018:Dec_2021) %>% arrange(GEOID)
  names(cocaine)[1] <- "state"
}

cocaine %>% filter(is.na(county))
grep("Jan_2018", names(cocaine)) # 5
grep("Dec_2021", names(cocaine)) # 52
cocaine <- cocaine %>% filter(GEOID %in% non_missing_GEOID | !is.na(HIDTA))
cocaine[,5:52][is.na(cocaine[,5:52])] <- 0
cocaine %>% filter(is.na(GEOID))
cocaine # 1518 counties in total.
# write.csv(cocaine, "cocaine other count HIDTA (06-28-2023).csv", row.names=F)


# crack weight
{
  crack <- seizures %>% filter(Drug=="Crack", Unit == "Kg")
  crack <- crack %>% mutate(Month2=month(Month, label=T))
  crack <- crack %>% group_by(state_name, county, Year, Month2) %>%
    summarise(weight=sum(Quantity)) %>%
    complete(state_name, county, Year, Month2) %>% 
    pivot_wider(names_from=c(Month2, Year), names_sep="_", values_from=weight)
  crack <- crack[, c(1:2, match(cal_order, names(crack)))]
  
  crack <- left_join(unique(coordinate.HIDTA[,c(7, 12, 14, 15)]), crack, by=c("state_name", "county")) %>%
    select(state_name, county, GEOID, HIDTA, Jan_2018:Dec_2021) %>% arrange(GEOID)
  names(crack)[1] <- "state"
}

crack %>% filter(is.na(county))
grep("Jan_2018", names(crack)) # 5
grep("Dec_2021", names(crack)) # 52
crack <- crack %>% filter(GEOID %in% non_missing_GEOID | !is.na(HIDTA))
crack[,5:52][is.na(crack[,5:52])] <- 0
crack %>% filter(is.na(GEOID))
crack # 1459 counties in total.
# write.csv(crack, "cocaine crack weight HIDTA (02-21-2023).csv", row.names=F)

# other cocaine weight
{
  cocaine <- seizures %>% filter(Drug %in% c("Cocaine", "Cocaine, Powder", "Cocaine Base",
                                             "Cocaine precursor", "Coca Leaves", "Cocaine metabolite"))
  cocaine <- cocaine %>% mutate(Month2=month(Month, label=T))
  cocaine <- cocaine %>% group_by(state_name, county, Year, Month2) %>%
    summarise(weight=sum(Quantity)) %>%
    complete(state_name, county, Year, Month2) %>% 
    pivot_wider(names_from=c(Month2, Year), names_sep="_", values_from=weight)
  cocaine <- cocaine[, c(1:2, match(cal_order, names(cocaine)))]
  
  cocaine <- left_join(unique(coordinate.HIDTA[,c(7, 12, 14, 15)]), cocaine, by=c("state_name", "county")) %>%
    select(state_name, county, GEOID, HIDTA, Jan_2018:Dec_2021) %>% arrange(GEOID)
  names(cocaine)[1] <- "state"
}

cocaine %>% filter(is.na(county))
grep("Jan_2018", names(cocaine)) # 5
grep("Dec_2021", names(cocaine)) # 52
cocaine <- cocaine %>% filter(GEOID %in% non_missing_GEOID | !is.na(HIDTA))
cocaine[,5:52][is.na(cocaine[,5:52])] <- 0
cocaine %>% filter(is.na(GEOID))
cocaine # 1518 counties in total.
# write.csv(cocaine, "cocaine other weight HIDTA (02-21-2023).csv", row.names=F)

# Annual figures of crack
coordinate_map <- coordinate.HIDTA
names(coordinate_map)[10] <- "county"
names(coordinate_map)[12] <- "state"
crack$`2018` <- apply(crack[,grep("2018", names(crack))], 1, sum)
crack$`2019` <- apply(crack[,grep("2019", names(crack))], 1, sum)
crack$`2020` <- apply(crack[,grep("2020", names(crack))], 1, sum)
crack$`2021` <- apply(crack[,grep("2021", names(crack))], 1, sum)

crack.map <- left_join(coordinate_map[, c(1:2,6:7,10,12)], crack[, c(3, 5:56)], by = "GEOID")

crack.map <- crack.map %>% 
  select(long:state, `2018`:`2021`) %>% 
  pivot_longer(c(-long, -lat, -group, -GEOID, -state, -county), names_to="Year", values_to="seizure_counts")
# crack.map$counts_per_pop <- crack.map$seizure_counts/all_cocaine.map$population

crack.map %>%
  ggplot(mapping = aes(long, lat, group = group, fill=seizure_counts)) +
  geom_polygon(color = "#000000", linewidth = .05) +
  facet_wrap(. ~ Year) +
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
# ggsave("crack annual seizure counts map 2020.pdf", seizure_counts_map, width=18, height=16, units="cm")


crack.map %>% 
  select(GEOID, `2018`:`2021`) %>% 
  pivot_longer(-GEOID, names_to="Year", values_to="seizure_counts") %>% 
  ggplot(mapping = aes(seizure_counts)) +
  geom_histogram(bins=100) +
  facet_wrap(. ~ Year) -> seizure_counts_hist
# ggsave("crack annual seizure counts histogram.pdf", seizure_counts_hist, width=20, height=15, units="cm")

  # 2020 plots
crack.map %>%
  filter(Year=="2020") %>% 
  ggplot(mapping = aes(long, lat, group = group, fill=seizure_counts)) +
  geom_polygon(color = "#000000", linewidth = .05) +
  scale_fill_viridis_c(alpha=0.9, na.value="white") +
  labs(fill = "Seizure Counts", x="", y="") + 
  theme_bw() + 
  theme(legend.position="bottom",
        legend.key.size = unit(0.3, 'cm'),
        legend.title = element_text(size=5),
        legend.text = element_text(size=5),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank()) -> seizure_counts_map

crack %>% 
  select(GEOID, `2020`) %>% 
  ggplot(mapping = aes(`2020`)) +
  geom_histogram(bins=100) +
  xlab("Seizure Count") -> seizure_counts_hist
# ggsave("crack annual seizure counts histogram 2020.pdf", seizure_counts_hist, width=10, height=7, units="cm")


# ts plots
seizure_ts_plot <- function(seizure, state_name, drug) {
  seizure %>% filter(state==state_name) %>% pivot_longer(cols=Jan_2018:Dec_2021, names_to="month", values_to="weight") %>% 
    ggplot(aes(x=month, y=weight, group=GEOID)) +
    geom_line(aes(color=county)) +
    labs(title=str_c(state_name, drug, sep="-"),
         xlab="Month-Year") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}

names(overdose)[1] <- "state_abbv"
cocaine_overdose <- overdose %>% filter(Indicator == "Cocaine (T40.5)" & Year %in% 2018:2021 & !(state_abbv %in% c("AK", "HI", "US")))
cocaine_overdose$Data.Value <- as.numeric(cocaine_overdose$Data.Value)
cocaine_overdose$Month <- str_sub(cocaine_overdose$Month, 1, 3)
cocaine_overdose$month <- str_c(cocaine_overdose$Month, cocaine_overdose$Year, sep="_")
cocaine_overdose$month <- factor(cocaine_overdose$month, levels=names(cocaine)[-(1:4)])
cocaine_overdose <- cocaine_overdose %>% rename(state=State.Name)
cocaine_overdose %>% complete(Year, Month)



state_name <- "California"
grid.arrange(seizure_ts_plot(cocaine, state_name, "cocaine"),
             seizure_ts_plot(crack, state_name, "crack"), ncol=2)
grid.arrange(seizure_ts_plot(cocaine_count, state_name, "cocaine"),
             seizure_ts_plot(crack_count, state_name, "crack"), ncol=2)
grid.arrange(
  cocaine_overdose %>% filter(state == state_name) %>% 
    ggplot(aes(x=month, y=Data.Value)) +
    geom_point(na.rm=T) +
    labs(title=str_c(state_name, "cocaine overdose", sep="-"),
         x="Month-Year", y="Counts") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)),
  cocaine %>% filter(state==state_name) %>% select(Jan_2018:Dec_2021) %>% 
    apply(2,sum) %>%
    ts(start=c(2018, 1), frequency=12) %>% 
    autoplot() + 
    labs(title=str_c(state_name, "all cocaine", sep="-"),
         x="Year", y="weight"),
  cocaine_count %>% filter(state==state_name) %>% select(Jan_2018:Dec_2021) %>% 
    apply(2,sum) %>%
    ts(start=c(2018, 1), frequency=12) %>% 
    autoplot() + 
    labs(title=str_c(state_name, "all cocaine", sep="-"),
         x="Year", y="count"), ncol=3
)

state_name <- "Florida"
grid.arrange(seizure_ts_plot(cocaine, state_name, "cocaine"),
             seizure_ts_plot(crack, state_name, "crack"), ncol=2)
grid.arrange(
  cocaine_overdose %>% filter(state == state_name) %>% 
    ggplot(aes(x=month, y=Data.Value)) +
    geom_point(na.rm=T) +
    labs(title=str_c(state_name, "cocaine overdose", sep="-"),
         x="Month-Year", y="Counts") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)),
  cocaine %>% filter(state==state_name) %>% select(Jan_2018:Dec_2021) %>% 
    apply(2,sum) %>%
    ts(start=c(2018, 1), frequency=12) %>% 
    autoplot() + 
    labs(title=str_c(state_name, "all cocaine", sep="-"),
         x="Year", y="weight"),
  cocaine_count %>% filter(state==state_name) %>% select(Jan_2018:Dec_2021) %>% 
    apply(2,sum) %>%
    ts(start=c(2018, 1), frequency=12) %>% 
    autoplot() + 
    labs(title=str_c(state_name, "all cocaine", sep="-"),
         x="Year", y="count"), ncol=3
)


state_name <- "New York"
grid.arrange(seizure_ts_plot(cocaine, state_name, "cocaine"),
             seizure_ts_plot(crack, state_name, "crack"), ncol=2)
grid.arrange(
  cocaine_overdose %>% filter(state == state_name) %>% 
    ggplot(aes(x=month, y=Data.Value)) +
    geom_point(na.rm=T) +
    labs(title=str_c(state_name, "cocaine overdose", sep="-"),
         x="Month-Year", y="Counts") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)),
  cocaine %>% filter(state==state_name) %>% select(Jan_2018:Dec_2021) %>% 
    apply(2,sum) %>%
    ts(start=c(2018, 1), frequency=12) %>% 
    autoplot() + 
    labs(title=str_c(state_name, "all cocaine", sep="-"),
         x="Year", y="weight"),
  cocaine_count %>% filter(state==state_name) %>% select(Jan_2018:Dec_2021) %>% 
    apply(2,sum) %>%
    ts(start=c(2018, 1), frequency=12) %>% 
    autoplot() + 
    labs(title=str_c(state_name, "all cocaine", sep="-"),
         x="Year", y="count"), ncol=3
)

