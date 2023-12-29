# setwd("/Users/euseongjang/Documents/R")
# setwd("C:/Users/gkfrj/Documents/R")
library(fpp2)
library(spdep)
library(readxl)
library(scales)
library(stringr)
library(urbnmapr)
library(tidyverse)
library(gridExtra)
library(lubridate)
{
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
  alpha <- 0.05
  nperm <- 999
  crack <- cbind(crack, matrix(0, nrow(crack), 48*3)) %>% as_tibble
  names(crack)[(which(names(crack) == "1"):which(names(crack) == "144"))] <- names(LISA3)[71:214]
  Jan_2018_index <- grep("Jan_2018", names(crack))[1]
  LISA_I.index <- grep("LISA_I", names(crack))[1]
  nb.obj.crack <- nb2listw(nb_crack, style="B")
  nrow(crack)
  seizures.crack <- crack[, 5:52]
  
  crack.rel <- read.csv("crack counts KNN5 R codes 999 two-sided relative (02-21-2023).csv") %>% as_tibble
  crack.rel.sig <- crack.rel %>%
    select(state:GEOID, Jan_2020, LISA_IJan0, LISA_PJan0, LISA_CJan0) %>%
    filter(LISA_CJan0 != 0) %>%
    arrange(desc(LISA_IJan0)) %>% 
    as.data.frame # Among HH regions, Highest LISA_I (Miami-Dade, Florida, 12086), Lowest LISA_I (Worcester, Massachusetts, 25027)
  crack.rel.insig <- crack.rel %>%
    select(state:GEOID, Jan_2020, LISA_IJan0, LISA_PJan0, LISA_CJan0) %>%
    filter(LISA_CJan0 == 0) %>%
    arrange(desc(LISA_IJan0)) %>% 
    as.data.frame
  {
    x = crack.rel$Jan_2020; listw=nb.obj.crack; nsim = 999; zero.policy = NULL; na.action = na.fail; 
    alternative = "two.sided"; mlvar = TRUE; spChk = NULL; xx=NULL;
    adjust.x = FALSE; sample_Ei = TRUE; iseed = NULL
  }
  alternative <- match.arg(alternative, c("greater", 
                                          "less", "two.sided"))
  NAOK <- deparse(substitute(na.action)) == "na.pass"
  x <- na.action(x)
  na.act <- attr(x, "na.action")
  rn <- attr(listw, "region.id")
  if (!is.null(na.act)) {
    subset <- !(1:length(listw$neighbours) %in% na.act)
    listw <- subset(listw, subset, zero.policy = zero.policy)
    excl <- class(na.act) == "exclude"
  }
  n <- length(listw$neighbours)
  res <- matrix(nrow = n, ncol = 9) #### res def
  gr <- punif((1:(nsim + 1))/(nsim + 1), 0, 1)
  ls <- rev(gr)
  ts <- (ifelse(gr > ls, ls, gr)) * 2
  if (alternative == "two.sided") {
    probs <- ts
    Prname <- "Pr(z != E(Ii))"
  }else if (alternative == "greater") {
    Prname <- "Pr(z > E(Ii))"
    probs <- gr
  }else {
    Prname <- "Pr(z < E(Ii))"
    probs <- ls
  }
  Prname_rank <- paste0(Prname, " Sim")
  Prname_sim <- "Pr(folded) Sim"
  colnames(res) <- c("Ii", "E.Ii", "Var.Ii", 
                     "Z.Ii", Prname, Prname_rank, Prname_sim, "Skewness", 
                     "Kurtosis")
  if (is.null(xx)) {
    if (adjust.x) {
      nc <- card(listw$neighbours) > 0L
      xx <- mean(x[nc], na.rm = NAOK)
    }
    else {
      xx <- mean(x, na.rm = NAOK)
    }
  }
  
  if (xx == "binary") {
    z <- x
    EIc <- -n*(z^2 * sapply(listw$weights, sum))
  } else {
    z <- x - xx
    EIc <- -(z^2 * sapply(listw$weights, sum))/((n - 1) * (sum(z * 
                                                                 z)/n))
  }
  
  lz <- lag.listw(listw, z, zero.policy = zero.policy, NAOK = NAOK)
  lbs <- c("Low", "High")
  quadr_ps <- interaction(cut(z, c(-Inf, 0, Inf), labels = lbs), 
                          cut(lz, c(-Inf, 0, Inf), labels = lbs), sep = "-")
  lx <- lag.listw(listw, x, zero.policy = zero.policy, NAOK = NAOK)
  lxx <- mean(lx, na.rm = NAOK)
  quadr <- interaction(cut(x, c(-Inf, xx, Inf), labels = lbs), 
                       cut(lx, c(-Inf, lxx, Inf), labels = lbs), sep = "-")
  xmed <- median(x, na.rm = NAOK)
  lxmed <- median(lx, na.rm = NAOK)
  quadr_med <- interaction(cut(x, c(-Inf, xmed, Inf), labels = lbs), 
                           cut(lx, c(-Inf, lxmed, Inf), labels = lbs), sep = "-")
  if (mlvar) {
    if (adjust.x) {
      s2 <- sum(z[nc]^2, na.rm = NAOK)/sum(nc)
    }
    else {
      s2 <- sum(z^2, na.rm = NAOK)/n
    }
  }else {
    if (adjust.x) {
      s2 <- sum(z[nc]^2, na.rm = NAOK)/(sum(nc) - 1)
    }
    else {
      s2 <- sum(z^2, na.rm = NAOK)/(n - 1)
    }
  }
  
  if (xx == "binary") {
    res[, 1] <- z * lz
  } else {
    res[, 1] <- (z/s2) * lz
  }
  crd <- card(listw$neighbours)
  permI_int <- function(i, zi, z_i, crdi, wtsi, nsim, Ii, binary) {
    res_i <- rep(as.numeric(NA), 8)
    if (crdi > 0) {
      sz_i <- matrix(sample(z_i, size = crdi * nsim, replace = TRUE), 
                     ncol = crdi, nrow = nsim)
      lz_i <- sz_i %*% wtsi
      
      if (binary == "binary") {
        res_p <- zi * lz_i
      } else {
        res_p <- (zi/s2) * lz_i
      }
      res_i[1] <- mean(res_p)
      res_i[2] <- var(res_p)
      xrank <- rank(c(res_p, Ii))[(nsim + 1L)]
      res_i[5] <- xrank
      rnk0 <- as.integer(sum(res_p >= Ii))
      drnk0 <- nsim - rnk0
      rnk <- ifelse(drnk0 < rnk0, drnk0, rnk0)
      res_i[6] <- rnk0
      res_i[7] <- e1071::skewness(res_p)
      res_i[8] <- e1071::kurtosis(res_p)
    }
    res_i
  }
  
  lww <- listw$weights
  Iis <- res[, 1]
  
  perm_index <- function(z_, crd_, j) {
    n_counties <- length(z_)
    return(sample((1:n_counties)[-j], size = crd_, replace = F))
  }
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
}



permI_dist_perm_z <- function(zi, z_i, crdi, wtsi, nsim, Ii, replacement) {
  if (replacement==T) {
    zi <- sample(c(zi, z_i), size = nsim, replace = T)
    sz_i_w_rep <- matrix(sample(c(zi, z_i), size = crdi * nsim, replace = T), 
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

lbs_sim <- c("L", "H")
max_z <- max(z)
min_z <- min(z)
nb_crack_k <- knn2nb(knearneigh(coords.crack, k=5), row.names=GEOIDS.crack)
listw_k <- nb2listw(nb_crack_k, style="B")
lz_simul <- lag.listw(listw_k, z, zero.policy = zero.policy, NAOK = NAOK)
upper_bound <- 17.18
min_sum_of_z <- min(lz_simul)
simulated_z_Jan2020 <- seq(min_z, max_z, by=0.7)
simulated_sum_of_z_Jan2020 <- seq(min_sum_of_z, upper_bound, by=0.25)
simulated_z_pairs <- merge(simulated_z_Jan2020, simulated_sum_of_z_Jan2020) %>%
  mutate(z=x, sum_of_z_neigh=y) %>% 
  select(z, sum_of_z_neigh)
simulated_z_pairs$z_label <- cut(simulated_z_pairs$z, c(-Inf, 0, Inf), labels = lbs_sim)
simulated_z_pairs$sum_of_z_neigh_label <- cut(simulated_z_pairs$sum_of_z_neigh, c(-Inf, 0, Inf), labels = lbs_sim)
simulated_z_pairs$LISA_C <- apply(simulated_z_pairs, 1,
                                  function(x) paste(x[3], x[4], sep=""))
simulated_z_pairs$pseudo_p <- numeric(nrow(simulated_z_pairs))
alpha_sim <- 0.05
crd_sim <- length(listw_k$weights[[1]])
wts_sim <- listw_k$weights[[1]]

nsim <- 9999
simulated_z_pairs_tested <- simulated_z_pairs
set.seed(100)
for (i in 1:nrow(simulated_z_pairs_tested)) {
  zi <- simulated_z_pairs_tested$z[i]
  sum_of_znj <- simulated_z_pairs_tested$sum_of_z_neigh[i]
  Ii <- zi*sum_of_znj/s2
  I_perm_w_rep_sim <- permI_dist(zi, z, crd_sim, wts_sim, nsim, Ii, replacement=T)
  R_plus <- sum(I_perm_w_rep_sim$I_perm[-(nsim+1)] >= Ii)
  pseudo_p <- (min(R_plus, nsim-R_plus)+1)/(nsim+1) # p-value for either lower or upper tail due to min
  current_LISA_C <- simulated_z_pairs_tested$LISA_C[i]
  simulated_z_pairs_tested$pseudo_p[i] <- pseudo_p
  simulated_z_pairs_tested$LISA_C[i] <- ifelse(pseudo_p > alpha_sim, "Insig", current_LISA_C)
}
simulated_z_pairs_tested$LISA_C <- as.factor(simulated_z_pairs_tested$LISA_C)


simulated_z_pairs_tested_perm_i <- simulated_z_pairs
set.seed(100)
for (i in 1:nrow(simulated_z_pairs_tested_perm_i)) {
  zi <- simulated_z_pairs_tested_perm_i$z[i]
  sum_of_znj <- simulated_z_pairs_tested_perm_i$sum_of_z_neigh[i]
  Ii <- zi*sum_of_znj/s2
  I_perm_w_rep_sim <- permI_dist_perm_z(zi, z, crd_sim, wts_sim, nsim, Ii, replacement=T)
  R_plus <- sum(I_perm_w_rep_sim$I_perm[-(nsim+1)] >= Ii)
  pseudo_p <- min(R_plus, nsim-R_plus)/(nsim+1)
  current_LISA_C <- simulated_z_pairs_tested_perm_i$LISA_C[i]
  simulated_z_pairs_tested_perm_i$pseudo_p[i] <- pseudo_p
  simulated_z_pairs_tested_perm_i$LISA_C[i] <- ifelse(pseudo_p > alpha_sim, "Insig", current_LISA_C)
}
simulated_z_pairs_tested_perm_i$LISA_C <- as.factor(simulated_z_pairs_tested_perm_i$LISA_C)

df0 <- data.frame(x=min_sum_of_z, y=1)
df1 <- data.frame(x=c(upper_bound-3, upper_bound+0.5), y=c(2,2))
df2 <- data.frame(x=c(upper_bound+0.5, upper_bound+0.5), y=c(2,-3))
df3 <- data.frame(x=c(upper_bound-3, upper_bound-3), y=c(2,-3))
df4 <- data.frame(x=c(upper_bound-3, upper_bound+0.5), y=c(-3,-3))
simulated_z_pairs_tested %>% 
  ggplot() +
  geom_point(aes(x=sum_of_z_neigh, y=z, color=LISA_C)) +
  geom_point(data=df0, aes(x=x, y=y),
             shape=1, size=6) +
  geom_line(data=df1, aes(x=x, y=y)) +
  geom_line(data=df2, aes(x=x, y=y)) +
  geom_line(data=df3, aes(x=x, y=y)) +
  geom_line(data=df4, aes(x=x, y=y)) +
  labs(color="",
    # title="Centered Seizure Counts vs. Sum of Neighbors' in Jan 2020 (k=5)",
    x=expression(sum(paste(w[ij],z[j]), "j=1", N)),
    y=expression(z[i]),
  ) +
  # geom_text(aes(x=-2.5, y=50, label="Minor\n island"), color="black") +
  # geom_text(aes(x=10, y=50, label="Minor cluster"), color="black") +
  scale_color_manual(values = c("Insig"="grey60",
                                "LL"="blue",
                                "LH"="steelblue",
                                "HL"="orange",
                                "HH"="red")) -> crack_Jan2020_zi_sig_region
# ggsave("Crack z_i Significance Region in Jan 2020 (zoom).png", crack_Jan2020_zi_sig_region, width=12, height=8, units="cm")

simulated_z_pairs_tested_perm_i %>% 
  ggplot() +
  geom_point(aes(x=sum_of_z_neigh, y=z, color=LISA_C)) +
  labs(
    # title="Centered Seizure Counts vs. Sum of Neighbors' in Jan 2020 (k=5)",
    x=expression(sum(paste(w[ij],z[j]), "j=1", N)),
    y=expression(z[i]),
  ) +
  xlim(-5, 17.19) +
  scale_color_manual(values = c("Insig"="grey60",
                                "LL"="blue",
                                "LH"="steelblue",
                                "HL"="orange",
                                "HH"="red")) -> crack_Jan2020_zi_perm_sig_region
# ggsave("Crack z_i Permuted Significance Region in Jan 2020 (zoom).png", crack_Jan2020_zi_perm_sig_region, width=12, height=8, units="cm")