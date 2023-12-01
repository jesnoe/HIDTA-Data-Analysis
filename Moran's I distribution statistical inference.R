# setwd("/Users/euseongjang/Documents/R")
# setwd("C:/Users/gkfrj/Documents/R")
{
  library(fpp2)
  library(spdep)
  library(readxl)
  library(scales)
  library(stringr)
  library(urbnmapr)
  library(fitdistrplus)
  library(tidyverse)
  library(gridExtra)
  library(lubridate)
  library(lamW)
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

## plot of simulated z vs. sum of neighbors' z
lbs_sim <- c("L", "H")
max_z <- max(z)
min_z <- min(z)

nb_crack_k <- knn2nb(knearneigh(coords.crack, k=5), row.names=GEOIDS.crack)
listw_k <- nb2listw(nb_crack_k, style="B")
lz_simul <- lag.listw(listw_k, z, zero.policy = zero.policy, NAOK = NAOK)
max_sum_of_z <- max(lz_simul)
min_sum_of_z <- min(lz_simul)
simulated_z_Jan2018 <- seq(min_z, max_z, by=1)
simulated_sum_of_z_Jan2018 <- seq(min_sum_of_z, max_sum_of_z, by=1)
simulated_z_pairs <- merge(simulated_z_Jan2018, simulated_sum_of_z_Jan2018) %>%
  mutate(z=x, sum_of_z_neigh=y) %>% 
  select(z, sum_of_z_neigh)
# simulated_z_pairs$z_label <- cut(simulated_z_pairs$z, c(-Inf, 0, Inf), labels = lbs_sim)
# simulated_z_pairs$sum_of_z_neigh_label <- cut(simulated_z_pairs$sum_of_z_neigh, c(-Inf, 0, Inf), labels = lbs_sim)

# label Low-Med-High
lbs3_sim <- c("L", "M", "H")
x_qunatile <- data.frame(table(x))
x_qunatile$quantile <- cumsum(x_qunatile$Freq) / sum(x_qunatile$Freq)
x_qunatile #%>% write.csv("crack_Jan2020 quantile.csv", row.names=F)
lx_qunatile <- data.frame(table(lx))
lx_qunatile$quantile <- cumsum(lx_qunatile$Freq) / sum(lx_qunatile$Freq)
lx_qunatile #%>% write.csv("crack_Jan2020 sum of neighbors quantile.csv", row.names=F)

simulated_z_pairs$z_label <- cut(simulated_z_pairs$z, c(-Inf, 0, 2*sd(x), Inf), labels = lbs3_sim)
simulated_z_pairs$sum_of_z_neigh_label <- cut(simulated_z_pairs$sum_of_z_neigh, c(-Inf, 0, 36-lxx, Inf), labels = lbs3_sim)

# Other threshold for High/Low
# simulated_z_pairs$z_label <- cut(simulated_z_pairs$z, c(-Inf, 4-xx, Inf), labels = lbs_sim)
# simulated_z_pairs$sum_of_z_neigh_label <- cut(simulated_z_pairs$sum_of_z_neigh, c(-Inf, 30-lxx, Inf), labels = lbs_sim)

simulated_z_pairs$LISA_C <- apply(simulated_z_pairs, 1,
                                  function(x) paste(x[3], x[4], sep=""))

simulated_z_pairs$pseudo_p <- numeric(nrow(simulated_z_pairs))
alpha_sim <- 0.05
crd_sim <- length(listw_k$weights[[1]])
wts_sim <- listw_k$weights[[1]]

nsim <- 999
simulated_z_pairs_tested <- simulated_z_pairs
set.seed(100)
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
simulated_z_pairs_tested$LISA_C <- as.factor(simulated_z_pairs_tested$LISA_C)

observed_z_sum <- data.frame(z=z, lz=lz)
simulated_z_pairs_tested %>% 
  ggplot(aes(x=sum_of_z_neigh, y=z, color=LISA_C)) +
  geom_point(size=0.9) +
  labs(
    title=paste0("Centered Seizure Counts vs. Sum of Neighbors' in Jan 2020 (k=5, M=", nsim, ")"),
    x=expression(sum(paste(w[ij],z[j]), "j=1", N)),
    y=expression(z[i])
    ) +
  scale_color_manual(values = c("Insig"="grey60",
                                "LL"="blue",
                                "LH"="steelblue",
                                "HL"="orange",
                                "HH"="red",
                                "Obs."="black")) +
  geom_point(data=observed_z_sum, aes(x=lz, y=z, color="Obs."))# -> crack_Jan2020_sig_region
# ggsave("Crack Significance Region in Jan 2020.png", crack_Jan2020_sig_region, width=15, height=10, units="cm")

simulated_z_pairs_tested %>% 
  ggplot(aes(x=sum_of_z_neigh, y=z, color=LISA_C)) +
  geom_point(size=0.9) +
  labs(color="",
    # title=paste0("Centered Seizure Counts vs. Sum of Neighbors' in Jan 2020 (k=5, M=", nsim, ")"),
    x=expression(sum(paste(w[ij],z[j]), "j=1", N)),
    y=expression(z[i])
  ) +
  scale_color_manual(breaks = c("Insig", "LL", "LM", "LH", "ML", "MM", "MH", "HL", "HM", "HH", "Obs."),
                     values = c("Insig"="grey60",
                                "LL"="blue",
                                # "LM"="#74add1",
                                "LH"="steelblue",
                                "ML"="#e0f3f8",
                                # "MM"="#ffffbf",
                                "MH"="#fee090",
                                "HL"="orange",
                                # "HM"="#f46d43",
                                "HH"="red",
                                "Obs."="black")) +
  geom_point(data=observed_z_sum, aes(x=lz, y=z, color="Obs."), size=0.7) -> crack_Jan2020_sig_region_extended
# ggsave("Crack Extended Significance Region in Jan 2020.png", crack_Jan2020_sig_region_extended, width=15, height=10, units="cm")
# ggsave("Crack z_i Extended Permuted Significance Region in Jan 2020.png", crack_Jan2020_sig_region_extended, width=15, height=10, units="cm")

simulated_z_pairs_tested %>% 
  ggplot(aes(x=sum_of_z_neigh, y=z, color=LISA_C)) +
  geom_point() +
  xlim(-10, 170) +
  labs(title="Centered Seizure Counts vs. Sum of Neighbors' in Jan 2020 (k=10)",
       x=expression(sum(paste(w[ij],z[j]), "j=1", N)),
       y=expression(z[i])
       ) +
  scale_color_manual(values = c("Insig"="grey60",
                                "LL"="#4575b4",
                                "LH"="#abd9e9",
                                "HL"="#fdae61",
                                "HH"="#d73027"))

simulated_z_pairs_tested %>% filter(pseudo_p > 0.05 & sum_of_z_neigh < -3.5 & z > 50)
which(simulated_z_pairs_tested$pseudo_p > 0.05 & simulated_z_pairs_tested$sum_of_z_neigh < -3.5 & simulated_z_pairs_tested$z > 50)
sum(x>0)/length(x)
1-(1-0.08)^5 # 0.3409185

## plot of simulated normal z vs. sum of neighbors' z
summary(x) # min=0, max=92
n_counties <- length(x)

simulated_z_norm <- function(x_mean, x_sd, n_points, alpha_sim, listw, nsim_, perm_z=F) {
  x_norm <- rnorm(n_points, x_mean, x_sd)
  xx_norm <- mean(x_norm)
  z_norm <- x_norm - xx_norm
  lz_norm <- lag.listw(listw, z_norm, zero.policy = zero.policy, NAOK = NAOK)
  
  max_z_norm <- max(z_norm)
  min_z_norm <- min(z_norm)
  max_sum_of_z_norm <- max(lz_norm)
  min_sum_of_z_norm <- min(lz_norm)
  simulated_z_norm <- seq(min_z_norm, max_z_norm, by=0.5)
  simulated_sum_of_z_norm <- seq(min_sum_of_z_norm, max_sum_of_z_norm, by=0.5)
  simulated_z_norm_pairs <- merge(simulated_z_norm, simulated_sum_of_z_norm) %>%
    mutate(z=x, sum_of_z_norm_neigh=y) %>% 
    select(z, sum_of_z_norm_neigh)
  simulated_z_norm_pairs$z_label <- cut(simulated_z_norm_pairs$z, c(-Inf, 0, Inf), labels = lbs_sim)
  simulated_z_norm_pairs$sum_of_z_norm_neigh_label <- cut(simulated_z_norm_pairs$sum_of_z_norm_neigh, c(-Inf, 0, Inf), labels = lbs_sim)
  simulated_z_norm_pairs$LISA_C <- apply(simulated_z_norm_pairs, 1,
                                         function(x) paste(x[3], x[4], sep=""))
  simulated_z_norm_pairs$pseudo_p <- numeric(nrow(simulated_z_norm_pairs))
  
  crd_sim <- 5
  wts_sim <- lww[[1]]
  
  simulated_z_norm_pairs_tested <- simulated_z_norm_pairs
  for (i in 1:nrow(simulated_z_norm_pairs_tested)) {
    zi <- simulated_z_norm_pairs_tested$z[i]
    sum_of_z_normnj <- simulated_z_norm_pairs_tested$sum_of_z_norm_neigh[i]
    Ii <- zi*sum_of_z_normnj/s2
    if (perm_z) {
      I_perm_w_rep_sim <- permI_dist_perm_z(zi, z_norm, crd_sim, wts_sim, nsim_, Ii, replacement=T)
    }else{
      I_perm_w_rep_sim <- permI_dist(zi, z_norm, crd_sim, wts_sim, nsim_, Ii, replacement=T)
    }
    R_plus <- sum(I_perm_w_rep_sim$I_perm[-(nsim_+1)] >= Ii)
    pseudo_p <- min(R_plus, nsim_-R_plus)/(nsim_+1)
    current_LISA_C <- simulated_z_norm_pairs_tested$LISA_C[i]
    simulated_z_norm_pairs_tested$pseudo_p[i] <- pseudo_p
    simulated_z_norm_pairs_tested$LISA_C[i] <- ifelse(pseudo_p > alpha_sim, "Insig", current_LISA_C)
  }
  simulated_z_norm_pairs_tested$LISA_C <- as.factor(simulated_z_norm_pairs_tested$LISA_C)
  return(simulated_z_norm_pairs_tested)
}


# coords.all <- counties.obs %>% group_by(GEOID) %>% summarise(x=mean(long), y=mean(lat))
# GEOIDS.all <- coords.all$GEOID
# coords.all <- coords.all[,-1]
nb_crack_k <- knn2nb(knearneigh(coords.crack, k=5), row.names=GEOIDS.crack)
listw_k <- nb2listw(nb_crack_k, style="B")
set.seed(5839)
simulated_z_norm_10_5_0.05 <- simulated_z_norm(x_mean=10, x_sd=5, n_points=n_counties, alpha_sim=0.05, listw=listw_k, nsim_=9999)
set.seed(5839)
simulated_z_norm_10_5_0.10 <- simulated_z_norm(x_mean=10, x_sd=5, n_points=n_counties, alpha_sim=0.10, listw=listw_k, nsim_=999)

simulated_z_norm_10_5_0.05 %>% 
  ggplot(aes(x=sum_of_z_norm_neigh, y=z, color=LISA_C)) +
  geom_point() +
  labs(
    # title=expression(paste("Simulated Normal Z vs. Sum of Neighbors' ", alpha, "=0.05, ", mu, "=10, ", sigma, "=5 (k=5)")),
    x=expression(sum(paste(w[ij],z[j]), "j=1", N)),
    y=expression(z[i])
    ) +
  xlim(-40, 40) +
  scale_color_manual(values = c("Insig"="grey60",
                                "LL"="blue",
                                "LH"="steelblue",
                                "HL"="orange",
                                "HH"="red"))

simulated_z_norm_10_5_0.10 %>% 
  ggplot(aes(x=sum_of_z_norm_neigh, y=z, color=LISA_C)) +
  geom_point() +
  xlim(-40, 40) +
  labs(title=expression(paste("Simulated Normal Z vs. Sum of Neighbors' ",
                              alpha,
                              "=0.10, ",
                              mu,
                              "=10, ",
                              sigma,
                              "=5 (k=5)")),
       ) +
  scale_color_manual(values = c("Insig"="grey60",
                                "LL"="blue",
                                "LH"="steelblue",
                                "HL"="orange",
                                "HH"="red"))
set.seed(5839)
x_norm <- rnorm(n_counties, 10, 5)
xx_norm <- mean(x_norm)
z_norm <- x_norm - xx_norm
lz_norm <- lag.listw(listw_k, z_norm, zero.policy = zero.policy, NAOK = NAOK)
simulated_z_sum <- data.frame(z=z_norm, lz=lz_norm)
sum_of_z_neigh_LB <- 5*sqrt(5)*qnorm(0.05, 0, 1, lower.tail=T)
sum_of_z_neigh_UB <- 5*sqrt(5)*qnorm(0.05, 0, 1, lower.tail=F)

simulated_z_norm_10_5_0.05 %>% 
  ggplot(aes(x=sum_of_z_norm_neigh, y=z, color=LISA_C)) +
  geom_point() +
  labs(
    # title="Centered Seizure Counts vs. Sum of Neighbors' in Jan 2020",
    x=expression(sum(paste(w[ij],z[j]), "j=1", N)),
    y=expression(z[i]),
    color="LISA_C"
  ) +
  scale_color_manual(values = c("Insig"="grey60",
                                "LL"="blue",
                                "LH"="steelblue",
                                "HL"="orange",
                                "HH"="red",
                                "Obs."="black")) +
  # scale_color_manual(values = c("Insig"="grey60",
  #                               "LL"="#4575b4",
  #                               "LH"="#abd9e9",
  #                               "HL"="#fdae61",
  #                               "HH"="#d73027",
  #                               "Obs."="black")) +
  geom_vline(xintercept=c(sum_of_z_neigh_LB, sum_of_z_neigh_UB)) +
  geom_point(data=simulated_z_sum, aes(x=lz, y=z, color="Obs.")) -> z_norm_0.05_z_plot
# ggsave("z norm 0.05 Exact Significance Region in Jan 2020.png", z_norm_0.05_z_plot, width=15, height=10, units="cm")

simulated_z_norm_10_5_0.05 %>% filter(sum_of_z_norm_neigh > sum_of_z_neigh_LB & sum_of_z_norm_neigh < sum_of_z_neigh_UB) %>% 
  pull(pseudo_p)

set.seed(5839)
x_norm <- rnorm(n_counties, 10, 5)
xx_norm <- mean(x_norm)
z_norm <- x_norm - xx_norm
lz_norm <- lag.listw(listw_k, z_norm, zero.policy = zero.policy, NAOK = NAOK)
simulated_z_sum <- data.frame(z=z_norm, lz=lz_norm)
simulated_z_norm_10_5_0.05 %>% 
  ggplot(aes(x=sum_of_z_norm_neigh, y=z, color=LISA_C)) +
  geom_point() +
  labs(title=expression(paste("Simulated Normal Z vs. Sum of Neighbors' ",
                              alpha,
                              "=0.05, ",
                              mu,
                              "=10, ",
                              sigma,
                              "=5 (k=5)")),
  ) +
  xlim(-40, 40) +
  scale_color_manual(values = c("Insig"="grey60",
                                "LL"="blue",
                                "LH"="steelblue",
                                "HL"="orange",
                                "HH"="red",
                                "Obs."="black")) +
  geom_point(data=simulated_z_sum, aes(x=lz, y=z, color="Obs."))

## plot of simulated z (and normal z) vs. sum of neighbors' z with permuting z
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
nb_crack_k <- knn2nb(knearneigh(coords.crack, k=5), row.names=GEOIDS.crack)
listw_k <- nb2listw(nb_crack_k, style="B")
lz_simul <- lag.listw(listw_k, z, zero.policy = zero.policy, NAOK = NAOK)
max_sum_of_z <- max(lz_simul)
min_sum_of_z <- min(lz_simul)
simulated_z_Jan2018 <- seq(min_z, max_z, by=1)
simulated_sum_of_z_Jan2018 <- seq(min_sum_of_z, max_sum_of_z, by=1)
simulated_z_pairs <- merge(simulated_z_Jan2018, simulated_sum_of_z_Jan2018) %>%
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

simulated_z_pairs_tested <- simulated_z_pairs
set.seed(100)
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
simulated_z_pairs_tested$LISA_C <- as.factor(simulated_z_pairs_tested$LISA_C)

simulated_z_pairs_tested %>% 
  ggplot(aes(x=sum_of_z_neigh, y=z, color=LISA_C)) +
  geom_point() +
  labs(
    # title="Centered Seizure Counts vs. Sum of Neighbors' in Jan 2020 (k=5)",
    x=expression(sum(paste(w[ij],z[j]), "j=1", N)),
    y=expression(z[i]),
    ) +
  scale_color_manual(values = c("Insig"="grey60",
                                "LL"="blue",
                                "LH"="steelblue",
                                "HL"="orange",
                                "HH"="red")) -> crack_Jan2020_zi_perm_sig_region
# ggsave("Crack z_i Permuted Significance Region in Jan 2020.png", crack_Jan2020_zi_perm_sig_region, width=15, height=10, units="cm")


nb_crack_k <- knn2nb(knearneigh(coords.crack, k=5), row.names=GEOIDS.crack)
listw_k <- nb2listw(nb_crack_k, style="B")
set.seed(5839)
simulated_z_norm_10_5_0.05 <- simulated_z_norm(x_mean=10, x_sd=5, n_points=n_counties, alpha_sim=0.05, listw=listw_k, perm_z=T)

simulated_z_norm_10_5_0.05 %>% 
  ggplot(aes(x=sum_of_z_norm_neigh, y=z, color=LISA_C)) +
  geom_point() +
  labs(
    # title=expression(paste("Simulated Normal Z vs. Sum of Neighbors' ", alpha, "=0.05, ", mu, "=10, ", sigma, "=5 (k=5)")),
    x=expression(sum(paste(w[ij],z[j]), "j=1", N)),
    y=expression(z[i]),
  ) +
  scale_color_manual(values = c("Insig"="grey60",
                                "LL"="blue",
                                "LH"="steelblue",
                                "HL"="orange",
                                "HH"="red")) -> z_norm_0.05_zi_perm_sig_region
# ggsave("z norm 0.05 z_i Permuted Significance Region in Jan 2020.png", z_norm_0.05_zi_perm_sig_region, width=15, height=10, units="cm")

# Significance regions of Chi-square random samples
x_var <- c(seq(0.1, 0.9, by=0.1), 1:90)
chi_sq_1 <- dchisq(x_var, 1)
chi_sq_5 <- dchisq(x_var, 5)
chi_sq_10 <- dchisq(x_var, 10)
chi_sq_20 <- dchisq(x_var, 20)
df <- data.frame(x_var, chi_sq_1, chi_sq_5, chi_sq_10, chi_sq_20)
df %>% ggplot(aes(x_var, chi_sq_1)) +
  geom_line() -> p1
df %>% ggplot(aes(x_var, chi_sq_5)) +
  geom_line() -> p2
df %>% ggplot(aes(x_var, chi_sq_10)) +
  geom_line() -> p3
df %>% ggplot(aes(x_var, chi_sq_20)) +
  geom_line() -> p4
grid.arrange(p1, p2, p3, p4, ncol=2)

simulated_z_chisq <- function(df, n_points, alpha_sim, listw, perm_z=F, nsim_) {
  x_chisq <- rchisq(n_points, df)
  xx_chisq <- mean(x_chisq)
  z_chisq <- x_chisq - xx_chisq
  lz_chisq <- lag.listw(listw, z_chisq, zero.policy = zero.policy, NAOK = NAOK)
  
  max_z_chisq <- max(z_chisq)
  min_z_chisq <- min(z_chisq)
  max_sum_of_z_chisq <- max(lz_chisq)
  min_sum_of_z_chisq <- min(lz_chisq)
  simulated_z_chisq <- seq(min_z_chisq, max_z_chisq, by=0.2)
  simulated_sum_of_z_chisq <- seq(min_sum_of_z_chisq, max_sum_of_z_chisq, by=0.2)
  simulated_z_chisq_pairs <- merge(simulated_z_chisq, simulated_sum_of_z_chisq) %>%
    mutate(z=x, sum_of_z_chisq_neigh=y) %>% 
    select(z, sum_of_z_chisq_neigh)
  simulated_z_chisq_pairs$z_label <- cut(simulated_z_chisq_pairs$z, c(-Inf, 0, Inf), labels = lbs_sim)
  simulated_z_chisq_pairs$sum_of_z_chisq_neigh_label <- cut(simulated_z_chisq_pairs$sum_of_z_chisq_neigh, c(-Inf, 0, Inf), labels = lbs_sim)
  simulated_z_chisq_pairs$LISA_C <- apply(simulated_z_chisq_pairs, 1,
                                          function(x) paste(x[3], x[4], sep=""))
  simulated_z_chisq_pairs$pseudo_p <- numeric(nrow(simulated_z_chisq_pairs))
  
  crd_sim <- 5
  wts_sim <- lww[[1]]
  
  simulated_z_chisq_pairs_tested <- simulated_z_chisq_pairs
  for (i in 1:nrow(simulated_z_chisq_pairs_tested)) {
    zi <- simulated_z_chisq_pairs_tested$z[i]
    sum_of_z_chisqnj <- simulated_z_chisq_pairs_tested$sum_of_z_chisq_neigh[i]
    Ii <- zi*sum_of_z_chisqnj/s2
    if (perm_z) {
      I_perm_w_rep_sim <- permI_dist_perm_z(zi, z_chisq, crd_sim, wts_sim, nsim_, Ii, replacement=T)
    }else{
      I_perm_w_rep_sim <- permI_dist(zi, z_chisq, crd_sim, wts_sim, nsim_, Ii, replacement=T)
    }
    R_plus <- sum(I_perm_w_rep_sim$I_perm[-(nsim_+1)] >= Ii)
    pseudo_p <- min(R_plus, nsim_-R_plus)/(nsim_+1)
    current_LISA_C <- simulated_z_chisq_pairs_tested$LISA_C[i]
    simulated_z_chisq_pairs_tested$pseudo_p[i] <- pseudo_p
    simulated_z_chisq_pairs_tested$LISA_C[i] <- ifelse(pseudo_p > alpha_sim, "Insig", current_LISA_C)
  }
  simulated_z_chisq_pairs_tested$LISA_C <- as.factor(simulated_z_chisq_pairs_tested$LISA_C)
  return(simulated_z_chisq_pairs_tested)
}

nb_crack_k <- knn2nb(knearneigh(coords.crack, k=5), row.names=GEOIDS.crack)
listw_k <- nb2listw(nb_crack_k, style="B")
nsim <- 9999
set.seed(5839)
simulated_z_chisq_1_0.05 <- simulated_z_chisq(df=1, n_points=n_counties, alpha_sim=0.05, listw=listw_k, perm_z=F, nsim_=nsim)
simulated_z_chvisq_5_0.05 <- simulated_z_chisq(df=5, n_points=n_counties, alpha_sim=0.05, listw=listw_k, perm_z=T, nsim_=nsim)
simulated_z_chisq_10_0.05 <- simulated_z_chisq(df=10, n_points=n_counties, alpha_sim=0.05, listw=listw_k, perm_z=T, nsim_=nsim)
simulated_z_chisq_20_0.05 <- simulated_z_chisq(df=20, n_points=n_counties, alpha_sim=0.05, listw=listw_k, perm_z=T, nsim_=nsim)

chisq_sum_of_z_neigh_LB <- qchisq(0.05, 5, lower.tail=T) #  1.145476
chisq_sum_of_z_neigh_UB <- qchisq(0.05, 5, lower.tail=F) # 11.0705
set.seed(5839)
x_chisq <- rchisq(n_counties, 1)
xx_chisq <- mean(x_chisq)
lx_chisq <- lag.listw(listw, x_chisq, zero.policy = zero.policy, NAOK = NAOK)
lxx_chisq <- mean(lx_chisq)
simulated_z_chisq_1_0.05 %>% 
  ggplot(aes(x=sum_of_z_chisq_neigh, y=z, color=LISA_C)) +
  geom_point() +
  labs(
    title=expression(paste("Simulated Chisq Z vs. Sum of Neighbors' ", alpha, "=0.05, ", nu, "=1, (k=5)")),
    x=expression(sum(paste(w[ij],z[j]), "j=1", N)),
    y=expression(z[i]),
  ) +
  geom_vline(xintercept=c(chisq_sum_of_z_neigh_LB-lxx_chisq, chisq_sum_of_z_neigh_UB-lxx_chisq)) +
  scale_color_manual(values = c("Insig"="grey60",
                                "LL"="#4575b4",
                                "LH"="#abd9e9",
                                "HL"="#fdae61",
                                "HH"="#d73027")) -> chisq_p1

simulated_z_chisq_5_0.05 %>% 
  ggplot(aes(x=sum_of_z_chisq_neigh, y=z, color=LISA_C)) +
  geom_point() +
  labs(
    title=expression(paste("Simulated Chisq Z vs. Sum of Neighbors' ", alpha, "=0.05, ", nu, "=5, (k=5)")),
    x=expression(sum(paste(w[ij],z[j]), "j=1", N)),
    y=expression(z[i]),
  ) +
  scale_color_manual(values = c("Insig"="grey60",
                                "LL"="blue",
                                "LH"="steelblue",
                                "HL"="orange",
                                "HH"="red")) -> chisq_p2
simulated_z_chisq_10_0.05 %>% 
  ggplot(aes(x=sum_of_z_chisq_neigh, y=z, color=LISA_C)) +
  geom_point() +
  labs(
    title=expression(paste("Simulated Chisq Z vs. Sum of Neighbors' ", alpha, "=0.05, ", nu, "=10, (k=5)")),
    x=expression(sum(paste(w[ij],z[j]), "j=1", N)),
    y=expression(z[i]),
  ) +
  scale_color_manual(values = c("Insig"="grey60",
                                "LL"="blue",
                                "LH"="steelblue",
                                "HL"="orange",
                                "HH"="red")) -> chisq_p3
simulated_z_chisq_20_0.05 %>% 
  ggplot(aes(x=sum_of_z_chisq_neigh, y=z, color=LISA_C)) +
  geom_point() +
  labs(
    title=expression(paste("Simulated Chisq Z vs. Sum of Neighbors' ", alpha, "=0.05, ", nu, "=20, (k=5)")),
    x=expression(sum(paste(w[ij],z[j]), "j=1", N)),
    y=expression(z[i]),
  ) +
  scale_color_manual(values = c("Insig"="grey60",
                                "LL"="blue",
                                "LH"="steelblue",
                                "HL"="orange",
                                "HH"="red")) -> chisq_p4

chisq_sim_plots <- grid.arrange(chisq_p1, chisq_p2, chisq_p3, chisq_p4, ncol=2)
# ggsave("Simulated Chisq Z vs. Sum of Neighbors.png", chisq_sim_plots, width=15, height=10, units="cm")


# Significance regions of gamma random samples
crack.rel %>% ggplot(aes(Jan_2020)) +
  geom_density() +
  labs(x="Seizure Count", title="Dist. of Seizure Counts I in Jan 2020")

x_var <- c(seq(0.1, 0.9, by=0.1), 1:90)
gamma_1 <- dgamma(x_var, shape=.5, rate=1)
gamma_2 <- dgamma(x_var, shape=.5, rate=0.5)
gamma_3 <- dgamma(x_var, shape=.25, rate=1)
gamma_4 <- dgamma(x_var, shape=.25, rate=.25)
df <- data.frame(x_var, gamma_1, gamma_2, gamma_3, gamma_4)
df %>% ggplot(aes(x_var, gamma_1)) +
  geom_line()
df %>% ggplot(aes(x_var, gamma_2)) +
  geom_line()
df %>% ggplot(aes(x_var, gamma_3)) +
  geom_line()
df %>% ggplot(aes(x_var, gamma_4)) +
  geom_line()

crack.rel %>%
  ggplot(aes(Jan_2020)) +
  geom_histogram(binwidth=1) +
  # xlim(1, 90) +
  labs(x="Seizure Count", title="Dist. of Seizure Counts I in Jan 2020")

nonzero_x_var <- c(1:90)
nonzero_gamma_1 <- pgamma(nonzero_x_var, shape=1, rate=.5)
nonzero_gamma_2 <- pgamma(nonzero_x_var, shape=2, rate=.5)
nonzero_gamma_3 <- pgamma(nonzero_x_var, shape=3, rate=.5)
nonzero_gamma_4 <- pgamma(nonzero_x_var, shape=3, rate=.5)
nonzero_df <- data.frame(nonzero_x_var, nonzero_gamma_1, nonzero_gamma_2, nonzero_gamma_3, nonzero_gamma_4)
nonzero_df %>% ggplot(aes(nonzero_x_var, nonzero_gamma_1)) +
  geom_line()
nonzero_df %>% ggplot(aes(nonzero_x_var, nonzero_gamma_2)) +
  geom_line()
nonzero_df %>% ggplot(aes(nonzero_x_var, nonzero_gamma_3)) +
  geom_line()
nonzero_df %>% ggplot(aes(nonzero_x_var, nonzero_gamma_4)) +
  geom_line()

Jan_2020_all_tb <- table(crack.rel %>% pull(Jan_2020))
Jan_2020_all_df <- data.frame(seizure_count=as.numeric(rownames(Jan_2020_all_tb)),
                          density=cumsum(Jan_2020_all_tb)/sum(Jan_2020_all_tb))
Jan_2020_all_df %>% ggplot(aes(seizure_count, density)) +
  geom_line()

Jan_2020_tb <- table(crack.rel %>% filter(Jan_2020 > 0) %>% pull(Jan_2020))
Jan_2020_df <- data.frame(seizure_count=as.numeric(rownames(Jan_2020_tb)),
                          density=cumsum(Jan_2020_tb)/sum(Jan_2020_tb))

Jan_2020_df
pgamma(Jan_2020_df$seizure_count, shape=5, rate=1)

Jan_2020_df %>% ggplot(aes(seizure_count, density)) +
  geom_line() +
  geom_line(data=nonzero_df, 
            mapping=aes(nonzero_x_var, nonzero_gamma_1),
            color=2) +
  geom_line(data=nonzero_df, 
            mapping=aes(nonzero_x_var, nonzero_gamma_2),
            color=3) +
  geom_line(data=nonzero_df, 
            mapping=aes(nonzero_x_var, nonzero_gamma_3),
            color=4) +
  geom_line(data=nonzero_df, 
            mapping=aes(nonzero_x_var, nonzero_gamma_4),
            color=5)

x_seizure <- Jan_2020_df$seizure_count
nonzero_gamma_1 <- pgamma(x_seizure, shape=.8, rate=.1)
nonzero_gamma_2 <- pgamma(x_seizure, shape=.8, rate=.12)
nonzero_gamma_3 <- pgamma(x_seizure, shape=.8, rate=.13)
nonzero_gamma_4 <- pgamma(x_seizure, shape=.8, rate=.14)
nonzero_df <- data.frame(x_seizure, nonzero_gamma_1, nonzero_gamma_2, nonzero_gamma_3, nonzero_gamma_4)

Jan_2020_df %>% ggplot(aes(seizure_count, density)) +
  geom_line() +
  geom_line(data=nonzero_df, 
            mapping=aes(x_seizure, nonzero_gamma_1),
            color=2) +
  geom_line(data=nonzero_df, 
            mapping=aes(x_seizure, nonzero_gamma_2),
            color=3) +
  geom_line(data=nonzero_df, 
            mapping=aes(x_seizure, nonzero_gamma_3),
            color=4) +
  geom_line(data=nonzero_df, 
            mapping=aes(x_seizure, nonzero_gamma_4),
            color=5)

# Find alpha and beta of gamma dist. with MLE
nonzero_Jan_2020 <- crack.rel %>% filter(Jan_2020 > 0) %>% pull(Jan_2020)
gamma_mle <- fitdist(nonzero_Jan_2020, distr = "gamma", method = "mle")
gamma_mle # shape=1.1415833, rate=0.1198246
shape_mle <- gamma_mle$estimate[1]
rate_mle <- gamma_mle$estimate[2]

zero_inflated_gamma <- function(xvar, shape, rate, k, n_nonzero, pi) {
  result <- (prod(k:(k-n_nonzero+1))/prod(1:n_nonzero)) * (1-pi)^n_nonzero * pi^(k-n_nonzero) * pgamma(xvar, n_nonzero*shape, rate, lower.tail=F)
  return(result)
}

sum_of_neigh_gamma <- function(xvar, shape, rate, k, pi) {
  result <- c()
  for (i in 1:k) {
    result <- sum(result, zero_inflated_gamma(xvar, shape, rate, k, i, pi))
  }
  return(result)
}
sum_of_neigh_gamma <- Vectorize(sum_of_neigh_gamma, vectorize.args="xvar")

gamma_upper_bound_df <- data.frame(z=5:30, 
                             upper_tail_probability=sum_of_neigh_gamma(5:30, shape=shape_mle, rate=rate_mle, k=5, pi=0.92))
gamma_upper_bound_df %>% ggplot(aes(z, upper_tail_probability)) +
  geom_line()

ub_binary <- function(ub, step, alpha) {
  upper_tail_probability <- sum_of_neigh_gamma(ub, shape=shape_mle, rate=rate_mle, k=5, pi=0.92)
  if (abs(upper_tail_probability - alpha) < 0.00001) {return(ub)}
  if (upper_tail_probability > alpha) {
    ub_binary(ub+step, step/2, alpha)
  }else{
    ub_binary(ub-step, step/2, alpha)
  }
}

ub_binary(20, 5, 0.05)

data.frame(z=seq(20, 22, by=0.01), 
           upper_tail_probability=sum_of_neigh_gamma(seq(20, 22, by=0.01), shape=shape_mle, rate=rate_mle, k=5, pi=0.92))[99:105,]

gamma_upper_tail <- 21.01


lbs_sim <- c("L", "H")
max_z <- max(z)
min_z <- min(z)

nb_crack_k <- knn2nb(knearneigh(coords.crack, k=5), row.names=GEOIDS.crack)
listw_k <- nb2listw(nb_crack_k, style="B")
lz_simul <- lag.listw(listw_k, z, zero.policy = zero.policy, NAOK = NAOK)
max_sum_of_z <- max(lz_simul)
min_sum_of_z <- min(lz_simul)
simulated_z_Jan2018 <- seq(min_z, max_z, by=1)
simulated_sum_of_z_Jan2018 <- seq(min_sum_of_z, max_sum_of_z, by=1)
simulated_z_pairs <- merge(simulated_z_Jan2018, simulated_sum_of_z_Jan2018) %>%
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
  pseudo_p <- (min(R_plus, nsim-R_plus)+1)/(nsim+1)
  current_LISA_C <- simulated_z_pairs_tested$LISA_C[i]
  simulated_z_pairs_tested$pseudo_p[i] <- pseudo_p
  simulated_z_pairs_tested$LISA_C[i] <- ifelse(pseudo_p > alpha_sim, "Insig", current_LISA_C)
}
simulated_z_pairs_tested$LISA_C <- as.factor(simulated_z_pairs_tested$LISA_C)
observed_z_sum <- data.frame(z=z, lz=lz)

simulated_z_pairs_tested %>% 
  ggplot(aes(x=sum_of_z_neigh, y=z, color=LISA_C)) +
  geom_point(size=0.9) +
  labs(
    # title=paste0("Centered Seizure Counts vs. Sum of Neighbors' in Jan 2020 (k=5, M=", nsim, ")"),
    x=expression(sum(paste(w[ij],z[j]), "j=1", N)),
    y=expression(z[i])
  ) +
  geom_vline(xintercept=gamma_upper_tail-5*xx) +
  scale_color_manual(values = c("Insig"="grey60",
                                "LL"="blue",
                                "LH"="steelblue",
                                "HL"="orange",
                                "HH"="red",
                                "Obs."="black")) +
  geom_point(data=observed_z_sum, aes(x=lz, y=z, color="Obs.")) -> crack_Jan2020_sig_region
# ggsave("Crack Exact Significance Region in Jan 2020.png", crack_Jan2020_sig_region, width=15, height=10, units="cm")


# gamma_mme <- fitdist(nonzero_Jan_2020, distr = "gamma", method = "mme")
# gamma_mge <- fitdist(nonzero_Jan_2020, distr = "gamma", method = "mge")
# gamma_mse <- fitdist(nonzero_Jan_2020, distr = "gamma", method = "mse")
# 
x_seizure <- Jan_2020_df$seizure_count
nonzero_gamma_1 <- pgamma(x_seizure, shape=.8, rate=.1)
nonzero_gamma_mle <- pgamma(x_seizure, shape=shape_mle, rate=rate_mle)
# nonzero_gamma_mge <- pgamma(x_seizure, shape=gamma_mge$estimate[1], rate=gamma_mge$estimate[2])
# nonzero_gamma_mse <- pgamma(x_seizure, shape=gamma_mse$estimate[1], rate=gamma_mse$estimate[2])
# nonzero_df <- data.frame(x_seizure, nonzero_gamma_mge, nonzero_gamma_mse)
# 
Jan_2020_df %>% ggplot(aes(seizure_count, density)) +
  geom_line() +
  geom_line(data=nonzero_df,
            mapping=aes(x_seizure, nonzero_gamma_1),
            color=2) +
  geom_line(data=nonzero_df,
            mapping=aes(x_seizure, nonzero_gamma_mle),
            color=3)
#   geom_line(data=nonzero_df, 
#             mapping=aes(x_seizure, nonzero_gamma_mge),
#             color=4) +
#   geom_line(data=nonzero_df, 
#             mapping=aes(x_seizure, nonzero_gamma_mse),
#             color=5)



# Find lambda and pi of zero-inflated Poisson dist. with MME (method of moment estimator )
var_x <- var(x)
lambda_mme <- (var_x+xx^2)/xx - 1
pi_mme <- (var_x-xx)/(var_x+xx^2-xx)

ZIP_mme_density <- c(pi_mme, rep(0, 10)) + (1-pi_mme)*dpois(0:10, lambda_mme)
plot(0:10, ZIP_mme_density)

# Find lambda and pi of zero-inflated Poisson dist. with MLE
s <- xx / (1- sum(x==0)/length(x))
lambda_mle <- lambertW0(-s*exp(-s)) + s
pi_mle <- 1 - xx/lambda_mle
ZIP_mle_density <- c(pi_mle, rep(0, 10)) + (1-pi_mle)*dpois(0:10, lambda_mle)
plot(0:10, ZIP_mle_density)
dgamma(1:10, shape_mle, rate=rate_mle)

zero_inflated_ZIP <- function(xvar, lambda, k, n_nonzero, pi) {
  result <- (prod(k:(k-n_nonzero+1))/prod(1:n_nonzero)) * (pi+(1-pi)*exp(-(k-n_nonzero)*lambda))^(k-n_nonzero) *
            (1-pi)^n_nonzero * ppois(xvar, n_nonzero*lambda, lower.tail=F)
  return(result)
}

sum_of_neigh_ZIP <- function(xvar, lambda, k, pi) {
  result <- c()
  for (i in 1:k) {
    result <- sum(result, zero_inflated_ZIP(xvar, lambda, k, i, pi))
  }
  return(result)
}

sum_of_neigh_ZIP <- Vectorize(sum_of_neigh_ZIP, vectorize.args="xvar")

ZIP_upper_bound_df <- data.frame(z=5:30, 
                                   upper_tail_probability=sum_of_neigh_ZIP(5:30, lambda=lambda_mle, k=5, pi=pi_mle))
ZIP_upper_bound_df %>% ggplot(aes(z, upper_tail_probability)) +
  geom_line()

ZIP_upper_tail <- 16

## plot of simulated x vs. sum of neighbors' x
lbs_sim <- c("L", "H")
max_x <- max(x)
min_x <- min(x)

nb_crack_k <- knn2nb(knearneigh(coords.crack, k=5), row.names=GEOIDS.crack)
listw_k <- nb2listw(nb_crack_k, style="B")
lx_simul <- lag.listw(listw_k, x, xero.policy = xero.policy, NAOK = NAOK)
max_sum_of_x <- max(lx_simul)
min_sum_of_x <- min(lx_simul)
simulated_x_Jan2018 <- seq(min_x, max_x, by=1)
simulated_sum_of_x_Jan2018 <- seq(min_sum_of_x, max_sum_of_x, by=1)
simulated_x_pairs <- merge(simulated_x_Jan2018, simulated_sum_of_x_Jan2018) %>%
  rename(sum_of_x_neigh=y)
simulated_x_pairs$x_label <- cut(simulated_x_pairs$x, c(-Inf, xx, Inf), labels = lbs_sim)
simulated_x_pairs$sum_of_x_neigh_label <- cut(simulated_x_pairs$sum_of_x_neigh, c(-Inf, xx, Inf), labels = lbs_sim)

simulated_x_pairs$LISA_C <- apply(simulated_x_pairs, 1,
                                  function(x) paste(x[3], x[4], sep=""))

simulated_x_pairs$pseudo_p <- numeric(nrow(simulated_x_pairs))
alpha_sim <- 0.05
crd_sim <- length(listw_k$weights[[1]])
wts_sim <- listw_k$weights[[1]]

nsim <- 9999
simulated_x_pairs_tested <- simulated_x_pairs
set.seed(100)
for (i in 1:nrow(simulated_x_pairs_tested)) {
  xi <- simulated_x_pairs_tested$x[i]
  sum_of_xnj <- simulated_x_pairs_tested$sum_of_x_neigh[i]
  Ii <- xi*sum_of_xnj/s2
  I_perm_w_rep_sim <- permI_dist(xi, x, crd_sim, wts_sim, nsim, Ii, replacement=T)
  R_plus <- sum(I_perm_w_rep_sim$I_perm[-(nsim+1)] >= Ii)
  pseudo_p <- (min(R_plus, nsim-R_plus)+1)/(nsim+1)
  current_LISA_C <- simulated_x_pairs_tested$LISA_C[i]
  simulated_x_pairs_tested$pseudo_p[i] <- pseudo_p
  simulated_x_pairs_tested$LISA_C[i] <- ifelse(pseudo_p > alpha_sim, "Insig", current_LISA_C)
}
simulated_x_pairs_tested$LISA_C <- as.factor(simulated_x_pairs_tested$LISA_C)

simulated_x_pairs_tested %>% 
  ggplot(aes(x=sum_of_x_neigh, y=x, color=LISA_C)) +
  geom_point(size=1.5) +
  labs(
    title=paste0("Centered Seizure Counts vs. Sum of Neighbors' in Jan 2020 (k=5, M=", nsim, ")"),
    x=expression(sum(paste(w[ij],x[j]), "j=1", N)),
    y=expression(x[i])
  ) +
  # geom_vline(xintercept=gamma_upper_tail) +
  geom_vline(xintercept=ZIP_upper_tail) +
  scale_color_manual(values = c("Insig"="grey60",
                                "LL"="blue",
                                "LH"="steelblue",
                                "HL"="orange",
                                "HH"="red"))


ref <- tibble(LISA_CJan0=0:4, label=c("Insig", "HH", "LL", "LH", "HL"))
crack.rel %>% 
  select(state:GEOID, Jan_2020, LISA_IJan0, LISA_PJan0, LISA_CJan0) %>% 
  mutate(x=x, sum_of_x_neigh=lx) %>% 
  left_join( ref, by="LISA_CJan0") %>%
  ggplot(aes(x=sum_of_x_neigh, y=x, color=label)) +
  geom_point() +
  labs(
    # title="Centered Seizure Counts vs. Sum of Neighbors' in Jan 2020",
    x=expression(sum(paste(w[ij],x[j]), "j=1", N)),
    y=expression(x[i]),
    color="LISA_C"
  ) +
  geom_vline(xintercept=gamma_upper_tail) +
  # geom_vline(xintercept=ZIP_upper_tail) +
  scale_color_manual(values = c("Insig"="grey60",
                                "LL"="blue",
                                "LH"="steelblue",
                                "HL"="orange",
                                "HH"="red"))

Jan_2020_result_9999 <- localmoran_abs(x, listw_k, nsim=9999, zero.policy=T, xx=NULL, alternative="two.sided")
Jan_2020_result_9999$LISA_C <- ifelse(Jan_2020_result_9999$quadr=="High-High", "HH",
                                      ifelse(Jan_2020_result_9999$quadr=="Low-Low", "LL",
                                             ifelse(Jan_2020_result_9999$quadr=="Low-High", "LH", "HL")))
Jan_2020_result_9999$LISA_C <- ifelse(Jan_2020_result_9999$`Pr(folded) Sim` <= 0.05, Jan_2020_result_9999$LISA_C, "Insig")

crack_Jan_2020 <- crack.rel
crack_Jan_2020$LISA_CJan0_9999 <- Jan_2020_result_9999$LISA_C
crack_Jan_2020 %>% 
  select(state:GEOID, Jan_2020, LISA_CJan0_9999) %>% 
  mutate(x=x, sum_of_x_neigh=lx, label=LISA_CJan0_9999) %>% 
  ggplot(aes(x=sum_of_x_neigh, y=x, color=label)) +
  geom_point() +
  labs(
    # title="Centered Seizure Counts vs. Sum of Neighbors' in Jan 2020",
    x=expression(sum(paste(w[ij],x[j]), "j=1", N)),
    y=expression(x[i]),
    color="LISA_C"
  ) +
  # geom_vline(xintercept=gamma_upper_tail) +
  geom_vline(xintercept=ZIP_upper_tail) +
  scale_color_manual(values = c("Insig"="grey60",
                                "LL"="blue",
                                "LH"="steelblue",
                                "HL"="orange",
                                "HH"="red"))
