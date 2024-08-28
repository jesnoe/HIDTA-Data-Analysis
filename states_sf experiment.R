# setwd("C:/Users/gkfrj/Documents/R")
library(spdep)
library(readxl)
library(scales)
library(stringr)
library(urbnmapr)
library(tidyverse)
library(gridExtra)
library(lubridate)
library(sp)
library(rgdal)

states_sf <- get_urbn_map(map = "states", sf = TRUE)

states_sf %>%
  left_join(statedata %>% filter(year == 2015) %>% select(state_fips, -state_name, horate), by="state_fips") %>% 
  filter(!(state_name %in% c("Alaska", "Hawaii"))) %>% 
  ggplot() +
  geom_sf(aes(fill=horate), 
          color = "#ffffff", size = 0.25)

counties_for_sf <- counties %>% mutate(GEOID=as.numeric(county_fips))
counties_sf <- st_as_sf(counties_for_sf %>% select(group, GEOID, long, lat, state_name, county_name), coords=c("long", "lat"))
counties_GEOID <- counties_for_sf$GEOID %>% unique %>% sort
counties_group <- counties_for_sf$group %>% unique %>% sort
counties_poly <- list()
for (i in 1:length(counties_group)) {
  # GEOID_i <- counties_GEOID[i]
  group_i <- counties_group[i]
  counties_poly[[as.character(group_i)]] <- st_polygon(counties_for_sf %>% filter(group == group_i) %>% select(long, lat) %>% as.matrix %>% list)
}
counties_sfc <- st_sfc(counties_poly)
counties_sf <- st_as_sf(counties_sfc, id=names(counties_poly)) %>% 
  mutate(GEOID = substr(id, 1, 5) %>% as.numeric)
counties_sf %>% filter(grepl("01003", id))
st_multipolygon(counties_sf$x[1:2])

counties_multipoly_sf <- counties_sf[1,] %>% select(-id)
for (j in 2:length(counties_GEOID)) {
  GEOID_j <- counties_GEOID[j]
  counties_sf_j <- counties_sf %>% filter(GEOID == GEOID_j) %>% select(-id)
  counties_multipoly_j <- st_multipolygon(counties_sf_j$x)
  counties_sf_j$x[1] <- counties_multipoly_j
  counties_multipoly_sf <- rbind(counties_multipoly_sf, counties_sf_j[1,])
}
counties_multipoly_sf

counties_neighbors <- poly2nb(counties_multipoly_sf)
counties_listw <- nb2listw(counties_neighbors, style="W", zero.policy=T)
counties_multipoly_sf$x_sample <- sample(seizures.crack$Jan_2020, nrow(counties_multipoly_sf), replace=T)

nperm <- 9999
t3 <- Sys.time()
localmoran_abs(counties_multipoly_sf$x_sample, counties_listw, nsim=nperm, zero.policy=T, xx=NULL, alternative="two.sided", perm.i=T)
t4 <- Sys.time()

t7 <- Sys.time()
localmoran_abs(counties_multipoly_sf$x_sample, counties_listw, nsim=nperm, zero.policy=T, xx=NULL, alternative="two.sided", perm.i=F)
t8 <- Sys.time()

sapply(counties_listw$weights, length) %>% table

counties_WGS84 <- st_transform(counties_sf, crs=4326)
coordinates(counties) <- c("long", "lat")
proj4string(counties) <- CRS("+proj=county +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84")



