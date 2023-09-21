# setwd("/Users/euseongjang/Documents/R")
# setwd("C:/Users/gkfrj/Documents/R")
library(fpp2)
library(spdep)
library(readxl)
library(urbnmapr)
library(tidyverse)
library(gridExtra)
library(lubridate)
library(xtable)
library(scales)
library(stringr)

crack <- read.csv("cocaine crack count HIDTA (06-28-2023).csv") %>% as_tibble
LISA3 <- read.csv("CountyKNN3.csv") %>% as_tibble %>% arrange(GEOID)

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
nperm <- 9999

## For crack
crack <- cbind(crack, matrix(0, nrow(crack), 48*3)) %>% as_tibble

names(crack)[(which(names(crack) == "1"):which(names(crack) == "144"))] <- names(LISA3)[71:214]
Jan_2018_index <- grep("Jan_2018", names(crack))[1]
LISA_I.index <- grep("LISA_I", names(crack))[1]
nrow(crack) # 1518
seizure.crack_4years <- crack[, Jan_2018_index:(LISA_I.index-1)] %>% apply(1, sum)
seizure.crack_2019 <- crack[, (Jan_2018_index+12):(Jan_2018_index+23)] %>% apply(1, sum)
seizure.crack_2020 <- crack[, (Jan_2018_index+24):(Jan_2018_index+35)] %>% apply(1, sum)

grep("LISA_I", names(crack))[1] # 53
LISA.org_agg <- crack[,1:3]
LISA.org_agg <- cbind(LISA.org_agg, seizure.crack_4years, seizure.crack_2019, seizure.crack_2020) %>% as_tibble

set.seed(100)
LISA.org_4years <- localmoran_abs(seizure.crack_4years, nb.obj.crack, nsim=nperm, zero.policy=T, xx=NULL, alternative="two.sided")
LISA.org_2019 <- localmoran_abs(seizure.crack_2019, nb.obj.crack, nsim=nperm, zero.policy=T, xx=NULL, alternative="two.sided")
LISA.org_2020 <- localmoran_abs(seizure.crack_2020, nb.obj.crack, nsim=nperm, zero.policy=T, xx=NULL, alternative="two.sided")

LISA.org_4years$LISA_C <- as.character(LISA.org_4years$quadr_ps)
LISA.org_4years$LISA_C <- ifelse(LISA.org_4years$`Pr(folded) Sim` <= alpha, LISA.org_4years$LISA_C, "Insig")
LISA.org_2019$LISA_C <- as.character(LISA.org_2019$quadr_ps)
LISA.org_2019$LISA_C <- ifelse(LISA.org_2019$`Pr(folded) Sim` <= alpha, LISA.org_2019$LISA_C, "Insig")
LISA.org_2020$LISA_C <- as.character(LISA.org_2020$quadr_ps)
LISA.org_2020$LISA_C <- ifelse(LISA.org_2020$`Pr(folded) Sim` <= alpha, LISA.org_2020$LISA_C, "Insig")

LISA.org_agg$LISA_C_4years <- LISA.org_4years$LISA_C
LISA.org_agg$LISA_C_2019 <- LISA.org_2019$LISA_C
LISA.org_agg$LISA_C_2020 <- LISA.org_2020$LISA_C


LISA.mod_agg <- crack[,1:3]
LISA.mod_agg <- cbind(LISA.mod_agg, seizure.crack_4years, seizure.crack_2019, seizure.crack_2020) %>% as_tibble

set.seed(100)
LISA.mod_4years <- localmoran_abs(seizure.crack_4years, nb.obj.crack, nsim=nperm, zero.policy=T, xx=NULL, alternative="two.sided", moderate=T)
LISA.mod_2019 <- localmoran_abs(seizure.crack_2019, nb.obj.crack, nsim=nperm, zero.policy=T, xx=NULL, alternative="two.sided", moderate=T)
LISA.mod_2020 <- localmoran_abs(seizure.crack_2020, nb.obj.crack, nsim=nperm, zero.policy=T, xx=NULL, alternative="two.sided", moderate=T)

LISA.mod_4years$LISA_C <- as.character(LISA.mod_4years$quadr_ps)
LISA.mod_4years$LISA_C <- ifelse(LISA.mod_4years$`Pr(folded) Sim` <= alpha, LISA.mod_4years$LISA_C, "Insig")
LISA.mod_2019$LISA_C <- as.character(LISA.mod_2019$quadr_ps)
LISA.mod_2019$LISA_C <- ifelse(LISA.mod_2019$`Pr(folded) Sim` <= alpha, LISA.mod_2019$LISA_C, "Insig")
LISA.mod_2020$LISA_C <- as.character(LISA.mod_2020$quadr_ps)
LISA.mod_2020$LISA_C <- ifelse(LISA.mod_2020$`Pr(folded) Sim` <= alpha, LISA.mod_2020$LISA_C, "Insig")

LISA.mod_agg$LISA_C_4years <- LISA.mod_4years$LISA_C
LISA.mod_agg$LISA_C_2019 <- LISA.mod_2019$LISA_C
LISA.mod_agg$LISA_C_2020 <- LISA.mod_2020$LISA_C


LISA.perm.i_agg <- crack[,1:3]
LISA.perm.i_agg <- cbind(LISA.perm.i_agg, seizure.crack_4years, seizure.crack_2019, seizure.crack_2020) %>% as_tibble

set.seed(100)
LISA.perm.i_4years <- localmoran_abs(seizure.crack_4years, nb.obj.crack, nsim=nperm, zero.policy=T, xx=NULL, alternative="two.sided", perm.i=T)
LISA.perm.i_2019 <- localmoran_abs(seizure.crack_2019, nb.obj.crack, nsim=nperm, zero.policy=T, xx=NULL, alternative="two.sided", perm.i=T)
LISA.perm.i_2020 <- localmoran_abs(seizure.crack_2020, nb.obj.crack, nsim=nperm, zero.policy=T, xx=NULL, alternative="two.sided", perm.i=T)

LISA.perm.i_4years$LISA_C <- as.character(LISA.perm.i_4years$quadr_ps)
LISA.perm.i_4years$LISA_C <- ifelse(LISA.perm.i_4years$`Pr(folded) Sim` <= alpha, LISA.perm.i_4years$LISA_C, "Insig")
LISA.perm.i_2019$LISA_C <- as.character(LISA.perm.i_2019$quadr_ps)
LISA.perm.i_2019$LISA_C <- ifelse(LISA.perm.i_2019$`Pr(folded) Sim` <= alpha, LISA.perm.i_2019$LISA_C, "Insig")
LISA.perm.i_2020$LISA_C <- as.character(LISA.perm.i_2020$quadr_ps)
LISA.perm.i_2020$LISA_C <- ifelse(LISA.perm.i_2020$`Pr(folded) Sim` <= alpha, LISA.perm.i_2020$LISA_C, "Insig")

LISA.perm.i_agg$LISA_C_4years <- LISA.perm.i_4years$LISA_C
LISA.perm.i_agg$LISA_C_2019 <- LISA.perm.i_2019$LISA_C
LISA.perm.i_agg$LISA_C_2020 <- LISA.perm.i_2020$LISA_C


LISA.both_agg <- crack[,1:3]
LISA.both_agg <- cbind(LISA.both_agg, seizure.crack_4years, seizure.crack_2019, seizure.crack_2020) %>% as_tibble

set.seed(100)
LISA.both_4years <- localmoran_abs(seizure.crack_4years, nb.obj.crack, nsim=nperm, zero.policy=T, xx=NULL, alternative="two.sided", moderate=T, perm.i=T)
LISA.both_2019 <- localmoran_abs(seizure.crack_2019, nb.obj.crack, nsim=nperm, zero.policy=T, xx=NULL, alternative="two.sided", moderate=T, perm.i=T)
LISA.both_2020 <- localmoran_abs(seizure.crack_2020, nb.obj.crack, nsim=nperm, zero.policy=T, xx=NULL, alternative="two.sided", moderate=T, perm.i=T)

LISA.both_4years$LISA_C <- as.character(LISA.both_4years$quadr_ps)
LISA.both_4years$LISA_C <- ifelse(LISA.both_4years$`Pr(folded) Sim` <= alpha, LISA.both_4years$LISA_C, "Insig")
LISA.both_2019$LISA_C <- as.character(LISA.both_2019$quadr_ps)
LISA.both_2019$LISA_C <- ifelse(LISA.both_2019$`Pr(folded) Sim` <= alpha, LISA.both_2019$LISA_C, "Insig")
LISA.both_2020$LISA_C <- as.character(LISA.both_2020$quadr_ps)
LISA.both_2020$LISA_C <- ifelse(LISA.both_2020$`Pr(folded) Sim` <= alpha, LISA.both_2020$LISA_C, "Insig")

LISA.both_agg$LISA_C_4years <- LISA.both_4years$LISA_C
LISA.both_agg$LISA_C_2019 <- LISA.both_2019$LISA_C
LISA.both_agg$LISA_C_2020 <- LISA.both_2020$LISA_C

# write.csv(LISA.org_agg, "crack counts KNN5 R codes 999 two-sided original (aggregated).csv", row.names=F)
# write.csv(LISA.mod_agg, "crack counts KNN5 R codes 999 two-sided moderate (aggregated).csv", row.names=F)
# write.csv(LISA.perm.i_agg, "crack counts KNN5 R codes 999 two-sided permute i (aggregated).csv", row.names=F)
# write.csv(LISA.both_agg, "crack counts KNN5 R codes 999 two-sided moderate permute i (aggregated).csv", row.names=F)

# crack Plots
LISA.org <- read.csv("crack counts KNN5 R codes 999 two-sided original (aggregated).csv") %>% as_tibble
LISA.mod <- read.csv("crack counts KNN5 R codes 999 two-sided moderate (aggregated).csv") %>% as_tibble
LISA.perm.i <- read.csv("crack counts KNN5 R codes 999 two-sided permute i (aggregated).csv") %>% as_tibble
LISA.both <- read.csv("crack counts KNN5 R codes 999 two-sided moderate permute i (aggregated).csv") %>% as_tibble

LISA_C.org <- LISA.org[,c(1:3, 7:9)]
LISA_C.mod <- LISA.mod[,c(1:3, 7:9)]
LISA_C.perm.i <- LISA.perm.i[,c(1:3, 7:9)]
LISA_C.both <- LISA.both[,c(1:3, 7:9)]
county.names <- LISA_C.org[, 1:3]


{
  LISA_C.org.map <- left_join(counties.obs[, c(1:2,6:7,10,12)], LISA_C.org[,-(1:2)], by = "GEOID")
  LISA_C.mod.map <- left_join(counties.obs[, c(1:2,6:7,10,12)], LISA_C.mod[,-(1:2)], by = "GEOID")
  LISA_C.perm.i.map <- left_join(counties.obs[, c(1:2,6:7,10,12)], LISA_C.perm.i[,-(1:2)], by = "GEOID")
  LISA_C.both.map <- left_join(counties.obs[, c(1:2,6:7,10,12)], LISA_C.both[,-(1:2)], by = "GEOID")
}


# comparison of results for 2018-2021
LISA_C.org.map %>% 
  ggplot(mapping = aes(long, lat, group = group, fill=LISA_C_4years)) +
  geom_polygon(color = "#000000", linewidth = .05) +
  scale_fill_manual(values = c("Insig"="grey60",
                               "LL"="blue",
                               "LH"="steelblue",
                               "HL"="orange",
                               "HH"="red"),
                    na.value = "white") +
  labs(fill = "LISA Labels", title="orignial 2018-2021", x="", y="") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) -> LISA_C_org_map_4years

LISA_C.org.map %>% 
  ggplot(mapping = aes(long, lat, group = group, fill=LISA_C_2019)) +
  geom_polygon(color = "#000000", linewidth = .05) +
  scale_fill_manual(values = c("Insig"="grey60",
                               "LL"="blue",
                               "LH"="steelblue",
                               "HL"="orange",
                               "HH"="red"),
                    na.value = "white") +
  labs(fill = "LISA Labels", title="orignial 2019", x="", y="") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) -> LISA_C_org_map_2019

LISA_C.org.map %>% 
  ggplot(mapping = aes(long, lat, group = group, fill=LISA_C_2020)) +
  geom_polygon(color = "#000000", linewidth = .05) +
  scale_fill_manual(values = c("Insig"="grey60",
                               "LL"="blue",
                               "LH"="steelblue",
                               "HL"="orange",
                               "HH"="red"),
                    na.value = "white") +
  labs(fill = "LISA Labels", title="orignial 2020", x="", y="") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) -> LISA_C_org_map_2020

LISA_C.mod.map %>%
  ggplot(mapping = aes(long, lat, group = group, fill=LISA_C_4years)) +
  geom_polygon(color = "#000000", linewidth = .05) +
  scale_fill_manual(values = c("Insig"="grey60",
                               "LL"="blue",
                               "LH"="steelblue",
                               "ML"="#abd9e9",
                               "MH"="#fee090",
                               "HL"="orange",
                               "HH"="red"),
                    na.value = "white") +
  labs(fill = "LISA Labels", title="moderate 2018-2021", x="", y="") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) -> LISA_C_mod_map_4years

LISA_C.mod.map %>%
  ggplot(mapping = aes(long, lat, group = group, fill=LISA_C_2019)) +
  geom_polygon(color = "#000000", linewidth = .05) +
  scale_fill_manual(values = c("Insig"="grey60",
                               "LL"="blue",
                               "LH"="steelblue",
                               "ML"="#abd9e9",
                               "MH"="#fee090",
                               "HL"="orange",
                               "HH"="red"),
                    na.value = "white") +
  labs(fill = "LISA Labels", title="moderate 2019", x="", y="") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) -> LISA_C_mod_map_2019

LISA_C.mod.map %>%
  ggplot(mapping = aes(long, lat, group = group, fill=LISA_C_2020)) +
  geom_polygon(color = "#000000", linewidth = .05) +
  scale_fill_manual(values = c("Insig"="grey60",
                               "LL"="blue",
                               "LH"="steelblue",
                               "ML"="#abd9e9",
                               "MH"="#fee090",
                               "HL"="orange",
                               "HH"="red"),
                    na.value = "white") +
  labs(fill = "LISA Labels", title="moderate 2020", x="", y="") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) -> LISA_C_mod_map_2020

LISA_C.perm.i.map %>%
  ggplot(mapping = aes(long, lat, group = group, fill=LISA_C_4years)) +
  geom_polygon(color = "#000000", linewidth = .05) +
  scale_fill_manual(values = c("Insig"="grey60",
                               "LL"="blue",
                               "LH"="steelblue",
                               "HL"="orange",
                               "HH"="red"),
                    na.value = "white") +
  labs(fill = "LISA Labels", title="permute i 2018-2021", x="", y="") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) -> LISA_C_perm.i_map_4years

LISA_C.perm.i.map %>%
  ggplot(mapping = aes(long, lat, group = group, fill=LISA_C_2019)) +
  geom_polygon(color = "#000000", linewidth = .05) +
  scale_fill_manual(values = c("Insig"="grey60",
                               "LL"="blue",
                               "LH"="steelblue",
                               "HL"="orange",
                               "HH"="red"),
                    na.value = "white") +
  labs(fill = "LISA Labels", title="permute i 2019", x="", y="") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) -> LISA_C_perm.i_map_2019

LISA_C.perm.i.map %>%
  ggplot(mapping = aes(long, lat, group = group, fill=LISA_C_2020)) +
  geom_polygon(color = "#000000", linewidth = .05) +
  scale_fill_manual(values = c("Insig"="grey60",
                               "LL"="blue",
                               "LH"="steelblue",
                               "HL"="orange",
                               "HH"="red"),
                    na.value = "white") +
  labs(fill = "LISA Labels", title="permute i 2020", x="", y="") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) -> LISA_C_perm.i_map_2020

LISA_C.both.map %>%
  ggplot(mapping = aes(long, lat, group = group, fill=LISA_C_4years)) +
  geom_polygon(color = "#000000", linewidth = .05) +
  scale_fill_manual(values = c("Insig"="grey60",
                               "LL"="blue",
                               "LH"="steelblue",
                               "ML"="#abd9e9",
                               "MH"="#fee090",
                               "HL"="orange",
                               "HH"="red"),
                    na.value = "white") +
  labs(fill = "LISA Labels", title="both 2018-2021", x="", y="") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) -> LISA_C_both_map_4years

LISA_C.both.map %>%
  ggplot(mapping = aes(long, lat, group = group, fill=LISA_C_2019)) +
  geom_polygon(color = "#000000", linewidth = .05) +
  scale_fill_manual(values = c("Insig"="grey60",
                               "LL"="blue",
                               "LH"="steelblue",
                               "ML"="#abd9e9",
                               "MH"="#fee090",
                               "HL"="orange",
                               "HH"="red"),
                    na.value = "white") +
  labs(fill = "LISA Labels", title="both 2019", x="", y="") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) -> LISA_C_both_map_2019

LISA_C.both.map %>%
  ggplot(mapping = aes(long, lat, group = group, fill=LISA_C_2020)) +
  geom_polygon(color = "#000000", linewidth = .05) +
  scale_fill_manual(values = c("Insig"="grey60",
                               "LL"="blue",
                               "LH"="steelblue",
                               "ML"="#abd9e9",
                               "MH"="#fee090",
                               "HL"="orange",
                               "HH"="red"),
                    na.value = "white") +
  labs(fill = "LISA Labels", title="both 2020", x="", y="") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) -> LISA_C_both_map_2020

# ggsave("Aggregated LISA Results/Crack_Count_LISA_C original two-sided (2018-2021).pdf", LISA_C_org_map_4years, width=15, height=10, units="cm")
# ggsave("Aggregated LISA Results/Crack_Count_LISA_C moderate two-sided (2018-2021).pdf", LISA_C_mod_map_4years, width=15, height=10, units="cm")
# ggsave("Aggregated LISA Results/Crack_Count_LISA_C permute i two-sided (2018-2021).pdf", LISA_C_perm.i_map_4years, width=15, height=10, units="cm")
# ggsave("Aggregated LISA Results/Crack_Count_LISA_C moderate permute i two-sided (2018-2021).pdf", LISA_C_both_map_4years, width=15, height=10, units="cm")
# ggsave("Aggregated LISA Results/Crack_Count_LISA_C original two-sided (2019).pdf", LISA_C_org_map_2019, width=15, height=10, units="cm")
# ggsave("Aggregated LISA Results/Crack_Count_LISA_C moderate two-sided (2019).pdf", LISA_C_mod_map_2019, width=15, height=10, units="cm")
# ggsave("Aggregated LISA Results/Crack_Count_LISA_C permute i two-sided (2019).pdf", LISA_C_perm.i_map_2019, width=15, height=10, units="cm")
# ggsave("Aggregated LISA Results/Crack_Count_LISA_C moderate permute i two-sided (2019).pdf", LISA_C_both_map_2019, width=15, height=10, units="cm")
# ggsave("Aggregated LISA Results/Crack_Count_LISA_C original two-sided (2020).pdf", LISA_C_org_map_2020, width=15, height=10, units="cm")
# ggsave("Aggregated LISA Results/Crack_Count_LISA_C moderate two-sided (2020).pdf", LISA_C_mod_map_2020, width=15, height=10, units="cm")
# ggsave("Aggregated LISA Results/Crack_Count_LISA_C permute i two-sided (2020).pdf", LISA_C_perm.i_map_2020, width=15, height=10, units="cm")
# ggsave("Aggregated LISA Results/Crack_Count_LISA_C moderate permute i two-sided (2020).pdf", LISA_C_both_map_2020, width=15, height=10, units="cm")

# Significance Region
permI_dist <- function(zi, z_i, crdi, wtsi, nsim, Ii, replacement) {
  if (replacement==T) {
    sz_i_w_rep <- matrix(sample(z_i, size = crdi * nsim, replace = T), 
                         ncol = crdi, nrow = nsim)
    lz_i_w_rep <- sz_i_w_rep %*% wtsi
    I_perm <- zi*lz_i_w_rep/s2
  } else {
    sz_i_wo_rep <- matrix(rep(0,crdi), ncol=crdi)
    for (i in 1:nsim){
      sz_i_wo_rep <- rbind(sz_i_wo_rep, matrix(sample(z_i, size=crdi, replace=F),
                                               ncol = crdi))
    }
    sz_i_wo_rep <- sz_i_wo_rep[-1,]
    lz_i_w_rep <- sz_i_wo_rep %*% wtsi
    I_perm <- zi*lz_i_w_rep/s2
  }
  I_perm <- c(I_perm, Ii)
  result <- data.frame(I_perm=I_perm,
                       observation=c(rep(0,nsim), 1))
  return(result)
}

permI_dist_perm_z <- function(zi, z_i, crdi, wtsi, nsim, Ii, replacement) {
  if (replacement==T) {
    zi <- sample(c(zi, z_i), size = nsim, replace = T)
    sz_i_w_rep <- matrix(sample(z_i, size = crdi * nsim, replace = T), 
                         ncol = crdi, nrow = nsim)
    lz_i_w_rep <- sz_i_w_rep %*% wtsi
    I_perm <- (zi/s2) * lz_i_w_rep
  } else {
    zi <- sample(c(zi, z_i), size = nsim, replace = F)
    sz_i_wo_rep <- matrix(rep(0,crdi), ncol=crdi)
    for (i in 1:nsim){
      sz_i_wo_rep <- rbind(sz_i_wo_rep, matrix(sample(z_i, size=crdi, replace=F),
                                               ncol = crdi))
    }
    sz_i_wo_rep <- sz_i_wo_rep[-1,]
    lz_i_w_rep <- sz_i_wo_rep %*% wtsi
    I_perm <- zi*lz_i_w_rep/s2
  }
  I_perm <- c(I_perm, Ii)
  result <- data.frame(I_perm=I_perm,
                       observation=c(rep(0,nsim), 1))
  return(result)
}

sig_region <- function(x, listw, data_period, nsim, alpha_sim, perm.i=F) {
  n <- length(x)
  xx <- mean(x)
  z <- x - xx
  lz <- lag.listw(listw, z)
  lx <- lag.listw(listw, x)
  lxx <- mean(lx)
  s2 <- sum(z^2)/n
  
  lbs_sim <- c("L", "H")
  max_z <- max(z)
  min_z <- min(z)
  
  lz_simul <- lag.listw(listw, z)
  max_sum_of_z <- max(lz_simul)
  min_sum_of_z <- min(lz_simul)
  simulated_z <- seq(min_z, max_z, by=(max_z-min_z)/100)
  simulated_sum_of_z <- seq(min_sum_of_z, max_sum_of_z, by=(max_sum_of_z-min_sum_of_z)/100)
  simulated_z_pairs <- merge(simulated_z, simulated_sum_of_z) %>%
    mutate(z=x, sum_of_z_neigh=y) %>% 
    select(z, sum_of_z_neigh)
  simulated_z_pairs$z_label <- cut(simulated_z_pairs$z, c(-Inf, 0, Inf), labels = lbs_sim)
  simulated_z_pairs$sum_of_z_neigh_label <- cut(simulated_z_pairs$sum_of_z_neigh, c(-Inf, 0, Inf), labels = lbs_sim)
  
  
  simulated_z_pairs$LISA_C <- apply(simulated_z_pairs, 1,
                                    function(x) paste(x[3], x[4], sep=""))
  
  simulated_z_pairs$pseudo_p <- numeric(nrow(simulated_z_pairs))
  crd_sim <- length(listw$weights[[1]])
  wts_sim <- listw$weights[[1]]
  
  simulated_z_pairs_tested <- simulated_z_pairs
  set.seed(100)
  
  if (perm.i) {
    for (i in 1:nrow(simulated_z_pairs_tested)) {
      zi <- simulated_z_pairs_tested$z[i]
      sum_of_znj <- simulated_z_pairs_tested$sum_of_z_neigh[i]
      Ii <- zi*sum_of_znj/s2
      I_perm_w_rep_sim <- permI_dist_perm_z(zi, z, crd_sim, wts_sim, nsim, Ii, replacement=T)
      R_plus <- sum(I_perm_w_rep_sim$I_perm[-(nsim+1)] >= Ii)
      pseudo_p <- min(R_plus, nsim-R_plus)/(nsim+1)
      current_LISA_C <- simulated_z_pairs_tested$LISA_C[i]
      simulated_z_pairs_tested$pseudo_p[i] <- pseudo_p
      simulated_z_pairs_tested$LISA_C[i] <- ifelse(pseudo_p > alpha_sim, "Insig", current_LISA_C)
    }
  }else{
    for (i in 1:nrow(simulated_z_pairs_tested)) {
      zi <- simulated_z_pairs_tested$z[i]
      sum_of_znj <- simulated_z_pairs_tested$sum_of_z_neigh[i]
      Ii <- zi*sum_of_znj/s2
      I_perm_w_rep_sim <- permI_dist(zi, z, crd_sim, wts_sim, nsim, Ii, replacement=T)
      R_plus <- sum(I_perm_w_rep_sim$I_perm[-(nsim+1)] >= Ii)
      pseudo_p <- (min(R_plus, nsim-R_plus)+1)/(nsim+1)
      current_LISA_C <- simulated_z_pairs_tested$LISA_C[i]
      simulated_z_pairs_tested$pseudo_p[i] <- pseudo_p
      simulated_z_pairs_tested$LISA_C[i] <- ifelse(pseudo_p > alpha_sim, "Insig", current_LISA_C)
    }
  }
  
  simulated_z_pairs_tested$LISA_C <- as.factor(simulated_z_pairs_tested$LISA_C)
  
  observed_z_sum <- data.frame(z=z, lz=lz)
  simulated_z_pairs_tested %>% 
    ggplot(aes(x=sum_of_z_neigh, y=z, color=LISA_C)) +
    geom_point(size=0.9) +
    labs(
      title=paste0(data_period, " (k=5, M=", nsim, ")"),
      x=expression(sum(paste(w[ij],z[j]), "j=1", N)),
      y=expression(z[i])
    ) +
    scale_color_manual(values = c("Insig"="grey60",
                                  "LL"="blue",
                                  "LH"="steelblue",
                                  "HL"="orange",
                                  "HH"="red",
                                  "Obs."="black")) +
    geom_point(data=observed_z_sum, aes(x=lz, y=z, color="Obs.")) -> z_sig_region
  return(z_sig_region)
}

org_sig_region_4years <- sig_region(x=seizure.crack_4years, listw=nb.obj.crack, data_period="2018-2021", nsim=9999, alpha_sim=0.05, perm.i=F)
org_sig_region_2019 <- sig_region(x=seizure.crack_2019, listw=nb.obj.crack, data_period="2019", nsim=9999, alpha_sim=0.05, perm.i=F)
org_sig_region_2020 <- sig_region(x=seizure.crack_2020, listw=nb.obj.crack, data_period="2020", nsim=9999, alpha_sim=0.05, perm.i=F)

perm.i_sig_region_4years <- sig_region(x=seizure.crack_4years, listw=nb.obj.crack, data_period="2018-2021", nsim=9999, alpha_sim=0.05, perm.i=T)
perm.i_sig_region_2019 <- sig_region(x=seizure.crack_2019, listw=nb.obj.crack, data_period="2019", nsim=9999, alpha_sim=0.05, perm.i=T)
perm.i_sig_region_2020 <- sig_region(x=seizure.crack_2020, listw=nb.obj.crack, data_period="2020", nsim=9999, alpha_sim=0.05, perm.i=T)
# ggsave("Aggregated LISA Results/Crack Significance Region Original in 2018-2021.png", org_sig_region_4years, width=15, height=10, units="cm")
# ggsave("Aggregated LISA Results/Crack Significance Region Original in 2019.png", org_sig_region_2019, width=15, height=10, units="cm")
# ggsave("Aggregated LISA Results/Crack Significance Region Original in 2020.png", org_sig_region_2020, width=15, height=10, units="cm")
# ggsave("Aggregated LISA Results/Crack Significance Region Permute i in 2018-2021.png", perm.i_sig_region_4years, width=15, height=10, units="cm")
# ggsave("Aggregated LISA Results/Crack Significance Region Permute i in 2019.png", perm.i_sig_region_2019, width=15, height=10, units="cm")
# ggsave("Aggregated LISA Results/Crack Significance Region Permute i in 2020.png", perm.i_sig_region_2020, width=15, height=10, units="cm")




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



