library(fpp2)
library(spdep)
library(readxl)
library(maps)
library(urbnmapr)
library(tidyverse)
library(rnaturalearth)


# counties.map <- st_as_sf(map("county", plot = FALSE, fill = TRUE))
# counties.map <- subset(counties.map, !grepl("alaska", counties.map$ID))
# 
# spdf_US <- ne_countries(country="United States of America")
# sfdf_US <- get_urbn_map(map="counties", sf=TRUE)
# spdf_US <- as(sfdf_US, "Spatial") %>% st_transform

counties$county_fips <- as.numeric(counties$county_fips)

LISA3 <- read.csv("CountyKNN3.csv") %>% as_tibble %>% arrange(GEOID)
state.fips <- read.csv("us-state-ansi-fips.csv")
names(state.fips)[1] <- "state"
names(state.fips)[2] <- "STATEFP"
LISA3 <- merge(LISA3, state.fips[,1:2], by="STATEFP") %>% as_tibble
LISA3 <- LISA3[,c(1, 215, 2:214)]

counties.obs <- counties %>% filter(county_fips %in% unique(LISA3$GEOID))
names(counties.obs)[7] <- "GEOID"

# conties.sp <- Polygons(counties.obs[,1:2], ID=counties.obs[,7])
# conties.sp <- SpatialPolygonsDataFrame(counties.obs[,1:2], counties.obs[,7])

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

grep("Jan_2018", names(LISA3)) # 21
grep("Dec_2021", names(LISA3)) # 68
grep("LISA_I", names(LISA3))[1] # 72
names(LISA3)[71+3*48] # 215

alpha <- 0.05
nperm <- 999
LISA_I.index <- grep("LISA_I", names(LISA3))[1]-1
nb.obj <- nb2listw(LISA3_nb, style="B")
seizures <- data.frame(Q1_2018=apply(LISA3[, 21:24], 1, sum))
for (i in 1:15) {
  seizures <- cbind(seizures, apply(LISA3[,(3*i + 21):(3*i + 23)], 1, sum))
}
quarter.names <- c("Q1_2018", "Q2_2018", "Q3_2018", "Q4_2018",
                   "Q1_2019", "Q2_2019", "Q3_2019", "Q4_2019",
                   "Q1_2020", "Q2_2020", "Q3_2020", "Q4_2020",
                   "Q1_2021", "Q2_2021", "Q3_2021", "Q4_2021")
names(seizures) <- quarter.names

absolute.MoransI <- LISA3[,1:71]
x.bar_p <- mean(seizures$Q1_2018)
set.seed(1000)
for (i in 1:16) {
  seizure <- t(seizures)[i,]
  localM.month <- localmoran_abs(seizure, nb.obj, nsim=nperm, zero.policy=T, xx=x.bar_p)
  localM.month$LISA_C <- ifelse(localM.month$mean=="High-High", 1,
                                ifelse(localM.month$mean=="Low-Low", 2,
                                       ifelse(localM.month$mean=="Low-High", 3, 4)))
  localM.month$LISA_C <- ifelse(localM.month$`Pr(folded) Sim` <= alpha, localM.month$LISA_C, 0)
  absolute.MoransI[,(69+3*i):(69+3*i+2)] <- localM.month[,c(1,13,7)]
  # x.bar_p <- mean(seizures[,i])
}

quarter.names2 <- rep(quarter.names, each=3)
for (j in 72:119) {
  prefix <- ifelse(j%%3 == 0, "LISA_I",
                   ifelse(j%%3 == 1, "LISA_C", "LISA_P"))
  names(absolute.MoransI)[j] <- paste(prefix, quarter.names2[j-71], sep="_")
}

LISA3[,72:82]; absolute.MoransI[,72:82]

# write.csv(absolute.MoransI, "CountyKNN3 Quarterly R codes 999 (relative).csv", row.names=F)
# write.csv(absolute.MoransI, "CountyKNN3 Quarterly R codes 9999 (relative).csv", row.names=F)
# write.csv(absolute.MoransI, "CountyKNN3 Quarterly R codes 99999 (relative).csv", row.names=F)
# write.csv(absolute.MoransI, "CountyKNN3 Quarterly R codes 999 (absolute base 1).csv", row.names=F)
# write.csv(absolute.MoransI, "CountyKNN3 Quarterly R codes 9999 (absolute).csv", row.names=F)
# write.csv(absolute.MoransI, "CountyKNN3 Quarterly R codes 99999 (absolute).csv", row.names=F)


## East Coast
coords.east <- counties.obs %>% filter(state_name %in% c("Tennessee", "Louisiana", "Pennsylvania", "Alabama","New York", "Illinois", "New Hampshire",
                                                         "Ohio", "Georgia", "North Carolina", "Michigan", "New Jersey", "Wisconsin", "Virginia",
                                                         "South Carolina", "Mississippi", "Massachusetts", "Rhode Island", "Indiana",
                                                         "Connecticut", "Delaware", "Arkansas", "West Virginia", "Kentucky")) %>%
  group_by(GEOID) %>% summarise(x=mean(long), y=mean(lat))
GEOIDS.east <- coords.east$GEOID
coords.east <- coords.east[,-1]
LISA3_nb_east <- knn2nb(knearneigh(coords.east, k=3), row.names=GEOIDS.east)

alpha <- 0.05
nperm <- 999
LISA_I.index <- grep("LISA_I", names(LISA3))[1]-1
nb.obj.east <- nb2listw(LISA3_nb_east, style="B")
LISA3.east <- LISA3 %>% filter(GEOID %in% GEOIDS.east)
seizures.east <- LISA3.east[, 21:68]

absolute.MoransI <- LISA3.east
x.bar_p <- mean(absolute.MoransI$Jan_2018)
set.seed(1000)
for (i in 1:48) {
  seizure.east <- t(seizures.east)[i,]
  localM.month <- localmoran_abs(seizure.east, nb.obj.east, nsim=nperm, zero.policy=T, xx=NULL)
  localM.month$LISA_C <- ifelse(localM.month$mean=="High-High", 1,
                                ifelse(localM.month$mean=="Low-Low", 2,
                                       ifelse(localM.month$mean=="Low-High", 3, 4)))
  localM.month$LISA_C <- ifelse(localM.month$`Pr(folded) Sim` <= alpha, localM.month$LISA_C, 0)
  absolute.MoransI[,(69+3*i):(69+3*i+2)] <- localM.month[,c(1,13,7)]
  x.bar_p <- mean(as.data.frame(LISA3.east)[,21+i-1])
}
LISA3.east[,72:82]; absolute.MoransI[,72:82]

# write.csv(absolute.MoransI, "CountyKNN3 East R codes 999 (relative).csv", row.names=F)
# write.csv(absolute.MoransI, "CountyKNN3 East R codes 999 (absolute base 1).csv", row.names=F)
# write.csv(absolute.MoransI, "CountyKNN3 East R codes 999 (absolute base t_1).csv", row.names=F)

## Diffusion Pathways Counts
LISA.list <- function(LISA, lag.t=1) {
  LISA.ts <- list()
  diffusion <- list()
  diffusion_summary <- list()
  state.summary <- list()
  for (state_name in unique(LISA$state)) {
    LISA_state <- LISA %>% filter(state == state_name)
    counties <- c()
    for (county in sort(unique(LISA_state$NAMELSAD))) { # County vector
      counties <- c(counties, LISA_state$NAMELSAD[grepl(county, LISA$NAMELSAD[LISA$state==state_name])])
    }
    
    LISA.ts[[state_name]] <- list()
    for (county in counties) {
      LISA.i <- LISA_state %>% filter(NAMELSAD == county)
      LISA.df <- data.frame(seizure.weights=t(LISA.i[,grep("Jan_2018", names(LISA.i)):grep("Dec_2021", names(LISA.i))]),
                            seizure.weights.lag=lag(t(LISA.i[,grep("Jan_2018", names(LISA.i)):grep("Dec_2021", names(LISA.i))])),
                            LISA_I=t(LISA.i[,grep("LISA_I", names(LISA.i))]),
                            LISA_I_lag=lag(t(LISA.i[,grep("LISA_I", names(LISA.i))]),lag.t),
                            LISA_C=t(LISA.i[,grep("LISA_C", names(LISA.i))]) %>% as.factor,
                            LISA_C_lag=lag(t(LISA.i[,grep("LISA_C", names(LISA.i))]),lag.t) %>% as.factor)
      row.names(LISA.df) <- NULL
      LISA.ts[[state_name]][[county]] <- LISA.df
    }
    
    
    diffusion[[state_name]] <- list()
    for (county in counties) {
      LISA_C_i <- LISA.ts[[state_name]][[county]]$LISA_C
      LISA_C_lag_i <- LISA.ts[[state_name]][[county]]$LISA_C_lag
      diffusion_pathways <- matrix(0, 5, 5) # Row: LISA_C at t, Col: LISA_C at t+1
      row.names(diffusion_pathways) <- c("Mid", "HH", "LL", "LH", "HL")
      colnames(diffusion_pathways) <- c("Mid", "HH", "LL", "LH", "HL")
      diffusion_i <- table(LISA_C_lag_i, LISA_C_i)
      rows_i <- row.names(diffusion_i) %>% as.numeric + 1
      cols_i <- colnames(diffusion_i) %>% as.numeric + 1
      for (j in 1:length(rows_i)) {
        for (k in 1:length(cols_i)){
          diffusion_pathways[rows_i[j], cols_i[k]] <- diffusion_i[j,k]
        }
      }
      diffusion_pathways <- diffusion_pathways[c(3,5,4,2,1), c(3,5,4,2,1)]
      N <- diffusion_pathways
      N.tot <- sum(N[-5,-5])
      n_stationary <- sum(c(diag(N)[-5], N[1,3], N[2,4], N[3,1], N[4,2]))
      
      diffusion_count <- data.frame(
        state=state_name,
        county=county,
        geiod=LISA_state$GEOID[LISA_state$NAMELSAD == county],
        n_total=N.tot,
        stationary=n_stationary
      )
      diffusion_count$change <- N.tot-diffusion_count$stationary
      n_increase <- sum(N[c(1,3), c(2,4)])
      n_decrease <- sum(N[c(2,4), c(1,3)])
      diffusion_count$increase <- n_increase
      diffusion_count$decrease <- n_decrease
      diffusion_count$non_adj_inc <- sum(N[1,c(2,4)])
      diffusion_count$adj_inc <- sum(N[3,c(2,4)])
      diffusion_count$non_adj_dec <- sum(N[4,c(1,3)])
      diffusion_count$adj_dec <- sum(N[2,c(1,3)])
      diffusion_count$local_High <- sum(LISA_C_i %in% c(1,4))
      diffusion_count$local_Low <- sum(LISA_C_i %in% c(2,3))
      diffusion_count$local_Ins <- sum(LISA_C_i==0)
      
      
      result <- list(diffusion_matrix=diffusion_pathways, count=diffusion_count)
      diffusion[[state_name]][[county]] <- result
    }
    
    diffusion_count_state <- data.frame(diffusion[[state_name]][[counties[1]]]$count)
    for (county in counties[-1]) {
      county_count <- data.frame(diffusion[[state_name]][[county]]$count)
      diffusion_count_state <- rbind(diffusion_count_state, county_count)
    }
    diffusion_summary[[state_name]][["count"]] <- diffusion_count_state
  }
  return(list(LISA.ts=LISA.ts, diffusion_summary=diffusion_summary))
}
LISA3.summary <- LISA.list(absolute.MoransI)

diffusion_summary <- LISA3.summary$diffusion_summary
diffusion_summary_count <- diffusion_summary$Alabama$count
for (state in unique(LISA3$state)[-1]) {
  diffusion_summary_count <- rbind(diffusion_summary_count, diffusion_summary[[state]]$count)
}
# write.csv(diffusion_summary_count, "KNN3 relative diffusion summary 4x4 R codes 999 (count).csv", row.names=F)
# write.csv(diffusion_summary_count, "KNN3 relative diffusion summary 4x4 R codes 9999 (count).csv", row.names=F)
# write.csv(diffusion_summary_count, "KNN3 relative diffusion summary 4x4 R codes 99999 (count).csv", row.names=F)
# write.csv(diffusion_summary_count, "KNN3 absolute base t_1 diffusion summary 4x4 R codes 999 (count).csv", row.names=F)
# write.csv(diffusion_summary_count, "KNN3 absolute diffusion summary 4x4 R codes 9999 (count).csv", row.names=F)
# write.csv(diffusion_summary_count, "KNN3 absolute diffusion summary 4x4 R codes 99999 (count).csv", row.names=F)

relative.MoransI <- read.csv("CountyKNN3 R codes 999 (relative).csv") %>% as_tibble
absolute.MoransI <- read.csv("CountyKNN3 R codes 999 (absolute).csv") %>% as_tibble


