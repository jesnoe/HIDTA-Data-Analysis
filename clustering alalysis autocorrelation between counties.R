library(fpp2)
library(spdep)
library(readxl)
library(scales)
library(stringr)
library(urbnmapr)
library(tidyverse)
library(gridExtra)
library(lubridate)

crack <- read.csv("cocaine crack count HIDTA (02-21-2023).csv") %>% as_tibble
cocaine <- read.csv("cocaine other count HIDTA (02-21-2023).csv") %>% as_tibble
marijuana <- read.csv("marijuana count HIDTA (02-21-2023).csv") %>% as_tibble
meth <- read.csv("meth count HIDTA (02-21-2023).csv") %>% as_tibble

counties.obs <- counties
names(counties.obs)[7] <- "GEOID"
counties.obs <- counties.obs %>% rename(state=state_name, county=county_name)
counties.obs$GEOID <- as.numeric(counties.obs$GEOID)
counties.obs <- counties.obs %>% filter(!(state %in% c("Alaska", "Hawaii")))



add_ccf_cols <- function(seizure_data, ref_states, lag_max) {
  
  n_counties <- nrow(seizure_data)
  n_ref <- length(ref_states)
  ref_seizure_list <- list()
  for (j in 1:n_ref) {
    ref_state <- ref_states[j]
    ref_seizure_list[[ref_state]] <- seizure_data %>% filter(state == ref_state) %>%
      select(Jan_2018:Dec_2021) %>% 
      apply(2, sum) %>% t %>% as.vector
  }
  
  cross_cor <- matrix(0, n_counties, n_ref*(lag_max+1))
  for (i in 1:n_counties) {
    ccf_i <- c()
    county_seizure_counts <- seizure_data[i,] %>% select(Jan_2018:Dec_2021) %>% t %>% as.vector
    if (var(county_seizure_counts) == 0) {
      ccf_i <- c(ccf_i, rep(0, n_ref*(lag_max+1)))
    } else {
      for (j in 1:n_ref) {
        ccf_ij <- ccf(x=ref_seizure_list[[j]],
                      y=county_seizure_counts,
                      lag.max=lag_max,
                      type="correlation",
                      plot=F)
        ccf_i <- c(ccf_i, ccf_ij$acf[(lag_max+1):(2*lag_max+1)])
      }
    }
  cross_cor[i,] <- ccf_i
  }
  cross_cor <- cross_cor %>% as.data.frame %>% as_tibble
  col_names <- c()
  for (j in 1:n_ref) {
    col_names <- c(col_names, paste(ref_states[j], "ccf", 0:lag_max, sep="_"))
  }
  names(cross_cor) <- col_names
  seizure_data <- cbind(seizure_data, cross_cor) %>% as_tibble
  return(seizure_data)
}

lag_max_input <- 6
crack_ccf <- add_ccf_cols(crack, ref_states=c("California", "Texas"), lag_max=lag_max_input)
cocaine_ccf <- add_ccf_cols(cocaine, ref_states=c("California", "Texas"), lag_max=lag_max_input)
marijuana_ccf <- add_ccf_cols(marijuana, ref_states=c("California", "Texas"), lag_max=lag_max_input)
meth_ccf <- add_ccf_cols(meth, ref_states=c("California", "Texas"), lag_max=lag_max_input)

crack_ccf_pca <- prcomp(crack_ccf %>% select(California_ccf_0:Texas_ccf_6), scale.=T, center=T)
cocaine_ccf_pca <- prcomp(cocaine_ccf %>% select(California_ccf_0:Texas_ccf_6), scale.=T, center=T)
marijuana_ccf_pca <- prcomp(marijuana_ccf %>% select(California_ccf_0:Texas_ccf_6), scale.=T, center=T)
meth_ccf_pca <- prcomp(meth_ccf %>% select(California_ccf_0:Texas_ccf_6), scale.=T, center=T)

crack_ccf_pca$sdev/sum(crack_ccf_pca$sdev); cocaine_ccf_pca$sdev/sum(cocaine_ccf_pca$sdev);
marijuana_ccf_pca$sdev/sum(marijuana_ccf_pca$sdev); meth_ccf_pca$sdev/sum(meth_ccf_pca$sdev)

crack_ccf_pca$rotation; cocaine_ccf_pca$rotation; marijuana_ccf_pca$rotation; meth_ccf_pca$rotation

crack_ccf_pca$x %>% as.data.frame %>% 
  ggplot(aes(PC1, PC2)) +
  geom_point() +
  labs(title="Crack PC1 vs. PC2") +
  ylim(-20, 20) -> p1

cocaine_ccf_pca$x %>% as.data.frame %>% 
  ggplot(aes(PC1, PC2)) +
  geom_point() +
  labs(title="Cocaine") +
  ylim(-20, 20) -> p2

marijuana_ccf_pca$x %>% as.data.frame %>% 
  ggplot(aes(PC1, PC2)) +
  geom_point() +
  labs(title="Marijuana") +
  ylim(-20, 20) -> p3

meth_ccf_pca$x %>% as.data.frame %>% 
  ggplot(aes(PC1, PC2)) +
  geom_point() +
  labs(title="Meth") +
  ylim(-20, 20) -> p4
grid.arrange(p1, p2, p3, p4, ncol=2)

crack_ccf %>% select(California_ccf_0:Texas_ccf_6) %>% summary
crack_ccf$California_ccf_1 %>% table

{
crack_ccf_map <- left_join(counties.obs[, c(1:2,6:7)], crack_ccf, by = "GEOID")
crack_ccf_map %>% ggplot(aes(x=long, y=lat, group=group, fill=California_ccf_0)) +
  geom_polygon(color = "#000000", size = .05) +
  scale_fill_gradientn(limits = c(-1,1),
                       breaks=c(-1, 0, 1),
                       colors=c("blue", "white", "red")) + 
  labs(fill = "CC_0 with CA", title="Crack") -> cor_p1

cocaine_ccf_map <- left_join(counties.obs[, c(1:2,6:7)], cocaine_ccf, by = "GEOID")
cocaine_ccf_map %>% ggplot(aes(x=long, y=lat, group=group, fill=California_ccf_0)) +
  geom_polygon(color = "#000000", size = .05) +
  scale_fill_gradientn(limits = c(-1,1),
                       breaks=c(-1, 0, 1),
                       colors=c("blue", "white", "red")) + 
  labs(fill = "CC_0 with CA", title="Other Cocaine") -> cor_p2

marijuana_ccf_map <- left_join(counties.obs[, c(1:2,6:7)], marijuana_ccf, by = "GEOID")
marijuana_ccf_map %>% ggplot(aes(x=long, y=lat, group=group, fill=California_ccf_0)) +
  geom_polygon(color = "#000000", size = .05) +
  scale_fill_gradientn(limits = c(-1,1),
                       breaks=c(-1, 0, 1),
                       colors=c("blue", "white", "red")) + 
  labs(fill = "CC_0 with CA", title="Marijuana") -> cor_p3

meth_ccf_map <- left_join(counties.obs[, c(1:2,6:7)], meth_ccf, by = "GEOID")
meth_ccf_map %>% ggplot(aes(x=long, y=lat, group=group, fill=California_ccf_0)) +
  geom_polygon(color = "#000000", size = .05) +
  scale_fill_gradientn(limits = c(-1,1),
                       breaks=c(-1, 0, 1),
                       colors=c("blue", "white", "red")) + 
  labs(fill="CC_0 with CA", title="Meth") -> cor_p4
}
grid.arrange(cor_p1, cor_p2, cor_p3, cor_p4, ncol=2)

crack[which(crack_ccf_pca$x[,2] < -5),] %>% select(state, county) %>% as.data.frame
cocaine[which(cocaine_ccf_pca$x[,2] < -5),] %>% select(state, county) %>% as.data.frame
meth[which(meth_ccf_pca$x[,2] < -5),] %>% select(state, county) %>% as.data.frame

# CCF for all pairs
seizure_data <- crack
n_counties <- nrow(seizure_data)
ccf_list <- list()
start <- Sys.time()
for (i in 1:n_rows) {
  seizure_data_i <- seizure_data[i,] %>% select(Jan_2018:Dec_2021) %>% t %>% as.vector
  
  for (j in (i+1):n_rows) {
    seizure_data_j <- seizure_data[j,] %>% select(Jan_2018:Dec_2021) %>% t %>% as.vector
    ccf_list[[paste(i, j, sep=",")]] <- ccf(seizure_data_i, 
                                            seizure_data_j,
                                            lag.max=lag_max,
                                            type="covariance",
                                            plot=F)$acf
  }
}
Sys.time() - start
