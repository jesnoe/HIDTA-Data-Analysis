library(tidyverse)
library(fpp2)

knn3.4x4 <- read.csv("county KNN3 diffusion summary 4x4 (count).csv")
knn5.4x4 <- read.csv("county KNN5 diffusion summary 4x4 (count).csv")
knn10.4x4 <- read.csv("county KNN10 diffusion summary 4x4 (count).csv")


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
populations <- read.csv("co-est2019-alldata.csv") %>% filter(COUNTY != 0) %>% select(SUMLEV, REGION, STATE, COUNTY, STNAME, CTYNAME, POPESTIMATE2019)
names(populations)[5] <- "state"
names(populations)[6] <- "county"
names(populations)[7] <- "population"
county.fips <- merge(county.fips, populations[,5:7], by=c("state", "county"), all.x=T)
names(county.fips)[3] <- "geiod"

knn3.4x4 <- merge(knn3.4x4, county.fips[,3:4], by="geiod")
knn5.4x4 <- merge(knn5.4x4, county.fips[,3:4], by="geiod")
knn10.4x4 <- merge(knn10.4x4, county.fips[,3:4], by="geiod")
pairs(knn3.4x4[,c(9:12, 13, 16)], pch=19)


plot(knn3.4x4$population/1000,
     knn3.4x4$local_High,
     xlab="Population (Thousands)",
     ylab="Counts of High",
     main="KNN3 LISA Results",
     pch=19)
plot(knn3.4x4$population/1000,
     knn3.4x4$local_High,
     xlab="Population (Thousands)",
     ylab="Counts of High",
     main="KNN3 LISA Results",
     pch=19,
     xlim=c(0, 2000))

(knn3.4x4 %>% filter(local_High > 10, population < 500000))[,c(2:3, 9:12, 13:16)]
(knn5.4x4 %>% filter(local_High > 10, population < 500000))[,c(2:3, 9:12, 13:16)]
(knn10.4x4 %>% filter(local_High > 10, population < 500000))[,c(2:3, 9:12, 13:16)]
knn3.4x4 %>% arrange(desc(population))

plot(knn3.4x4$population/1000,
     knn3.4x4$adj_inc,
     xlab="Population (Thousands)",
     ylab="Counts of adj_inc",
     main="KNN3 LISA Results",
     pch=19)
plot(knn3.4x4$population/1000,
     knn3.4x4$adj_inc,
     xlab="Population (Thousands)",
     ylab="Counts of adj_inc",
     main="KNN3 LISA Results",
     pch=19,
     xlim=c(0, 2000))
(knn3.4x4 %>% filter(adj_inc > 1, population < 500000))[,c(2:3, 9:12, 13:16)]
(knn5.4x4 %>% filter(adj_inc > 1, population < 500000))[,c(2:3, 9:12, 13:16)]
(knn10.4x4 %>% filter(adj_inc > 1, population < 500000))[,c(2:3, 9:12, 13:16)]

knn3.LISA <- read.csv("CountyKNN3.csv")
state.fips <- read.csv("us-state-ansi-fips.csv")
names(state.fips)[1] <- "state"
names(state.fips)[2] <- "STATEFP"
knn3.LISA <- merge(knn3.LISA, state.fips[,1:2], by="STATEFP") %>% as.tibble
knn3.LISA <- knn3.LISA[,c(1, 215, 2:214)] # Bring state name to column 2
knn3.LISA[knn3.LISA$NAMELSAD=="Kenedy County", 21:68] %>% t
autoplot(ts(t(knn3.LISA[knn3.LISA$NAMELSAD=="Kenedy County", 21:68]),
            start=c(2018, 1), frequency=12),
         xlab="Year",
         ylab="Seizure Weights (kg)",
         main="Seizure Wieghts in Kenedy County, TX")
knn3.LISA[knn3.LISA$NAMELSAD=="Kenedy County", seq(73, by=3, 215)] %>% t
autoplot(ts(t(knn3.LISA[knn3.LISA$NAMELSAD=="Kenedy County", seq(73, by=3, 215)]),
             start=c(2018, 2), frequency=12),
         xlab="Year",
         ylab="LISA Class")

## OH, SW Results
OH.4x4 <- merge(OH.4x4, county.fips[,3:4], by="geiod")
pairs(OH.4x4[,c(9:12, 13, 16)], pch=19)
plot(OH.4x4$population/1000,
     OH.4x4$local_High,
     xlab="Population (Thousands)",
     ylab="Counts of High",
     main="Ohio 4x4",
     pch=19)
plot(OH.4x4$population/1000,
     OH.4x4$local_High,
     xlab="Population (Thousands)",
     ylab="Counts of High",
     main="Ohio 4x4",
     pch=19,
     xlim=c(0, 2000))

(OH.4x4 %>% filter(local_High > 10, population < 1000000))[,c(2:3, 9:12, 13:16)]
OH.4x4 %>% arrange(desc(population))

SW.4x4 <- merge(SW.4x4, county.fips[,3:4], by="geiod")
pairs(SW.4x4[,c(9:12, 13, 16)], pch=19)
plot(SW.4x4$population/1000,
     SW.4x4$local_High,
     xlab="Population (Thousands)",
     ylab="Counts of High",
     main="Southwest 4x4",
     pch=19)
plot(SW.4x4$population/1000,
     SW.4x4$local_High,
     xlab="Population (Thousands)",
     ylab="Counts of High",
     main="Southwest 4x4",
     pch=19,
     xlim=c(0, 2000))

(SW.4x4 %>% filter(local_High > 4, population < 1000000))[,c(2:3, 9:12, 13:16)]
SW.4x4 %>% arrange(desc(population))

SW.5x5 <- merge(SW.5x5, county.fips[,3:4], by="geiod")
plot(SW.5x5$population/1000,
     SW.5x5$adj_inc,
     xlab="Population (Thousands)",
     ylab="Counts of Adjacent Increases",
     main="Southwest 5x5",
     pch=19)
plot(SW.5x5$population/1000,
     SW.5x5$adj_inc,
     xlab="Population (Thousands)",
     ylab="Counts of Adjacent Increases",
     main="Southwest 5x5",
     pch=19,
     xlim=c(0, 2000))
(SW.5x5 %>% filter(adj_inc > 4, population < 1000000))[,c(2:3, 9:12, 13:16)]

plot(SW.5x5$population/1000,
     SW.5x5$non_adj_inc,
     xlab="Population (Thousands)",
     ylab="Counts of Non-adjacent Increases",
     main="Southwest 5x5",
     pch=19)
plot(SW.5x5$population/1000,
     SW.5x5$non_adj_inc,
     xlab="Population (Thousands)",
     ylab="Counts of Non-adjacent Increases",
     main="Southwest 5x5",
     pch=19,
     xlim=c(0, 2000))
(SW.5x5 %>% filter(non_adj_inc > 2, population < 1000000))[,c(2:3, 9:12, 13:16)]
(SW.5x5 %>% filter(non_adj_inc > 4))[,c(2:3, 9:12, 13:16)]
SW.ts$`New Mexico`$`Lincoln County`
SW[SW$NAMELSAD == "Lincoln County", 21:66] %>% t %>% as.vector
