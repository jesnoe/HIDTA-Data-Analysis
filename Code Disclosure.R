library(spdep)
library(readxl)
library(scales)
library(stringr)
library(urbnmapr)
library(fitdistrplus)
library(tidyverse)
library(gridExtra)
library(lubridate)

# Setting the same neighborhood as HIDTA data we used in our paper (crack cocaine seizures in January 2020)
{
  crack.counties <- read.csv("Data/HIDTA crack cocaine counties.csv") %>% as_tibble
  counties.obs <- counties
  names(counties.obs)[7] <- "GEOID"
  counties.obs <- counties.obs %>% rename(state=state_name, county=county_name)
  counties.obs$GEOID <- as.numeric(counties.obs$GEOID)
  counties.obs <- counties.obs %>% filter(!(state %in% c("Alaska", "Hawaii")))
  coords.crack <- counties.obs %>% filter(GEOID %in% crack.counties$GEOID) %>%
    group_by(GEOID) %>% summarise(x=mean(long), y=mean(lat))
  GEOIDS.crack <- coords.crack$GEOID
  coords.crack <- coords.crack[,-1]
  nb_crack <- knn2nb(knearneigh(coords.crack, k=5), row.names=GEOIDS.crack)
  
  {
    listw=nb2listw(nb_crack, style="B"); zero.policy = NULL; na.action = na.fail; 
    alternative = "two.sided"; mlvar = TRUE; spChk = NULL; xx=NULL;
    adjust.x = FALSE; sample_Ei = TRUE; iseed = NULL
  }
  alternative <- match.arg(alternative, c("greater", "less", "two.sided"))
  lww <- listw$weights
  NAOK <- deparse(substitute(na.action)) == "na.pass"
}

# Results for a normal simulation data
set.seed(5839)
n_counties <- nrow(crack.counties)
x_norm <- rnorm(n_counties, 10, 5)
xx_norm <- mean(x_norm)
z_norm <- x_norm - xx_norm
nb_crack_k <- knn2nb(knearneigh(coords.crack, k=5), row.names=GEOIDS.crack)
listw_k <- nb2listw(nb_crack_k, style="B")
lz_norm <- lag.listw(listw_k, z_norm, zero.policy = zero.policy, NAOK = NAOK)
simulated_z_sum <- data.frame(z=z_norm, lz=lz_norm)

simulated_z_sum %>% 
  ggplot(aes(x=lz, y=z)) +
  geom_point() +
  labs(
    x=expression(sum(paste(w[ij],z[j]), "j=1", N)),
    y=expression(z[i])
  ) -> plot_z_norm_simulation_data
ggsave("Result/scatter plot of z norm simulation data.pdf", plot_z_norm_simulation_data, width=15, height=10, units="cm")


lbs_sim <- c("L", "H")

permI_dist <- function(zi, z_i, crdi, wtsi, nsim, Ii, replacement) {
  z <- c(zi, z_i)
  s2 <- sum(z^2)
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
  z <- c(zi, z_i)
  s2 <- sum(z^2)
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

simulated_z_norm <- function(z_data, lz_norm, alpha_sim, listw, nsim_, perm_z=F) {
  z_norm <- z_data$z
  lz_norm <- z_data$lz
  s2 <- sum(z_norm^2)
  
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
  lww <- listw$weights
  wts_sim <- lww[[1]]
  
  simulated_z_norm_pairs_tested <- simulated_z_norm_pairs
  for (i in 1:nrow(simulated_z_norm_pairs_tested)) {
    zi <- simulated_z_norm_pairs_tested$z[i]
    sum_of_z_norm_neigh_i <- simulated_z_norm_pairs_tested$sum_of_z_norm_neigh[i]
    Ii <- zi*sum_of_z_norm_neigh_i/s2
    if (perm_z) {
      I_perm_w_rep_sim <- permI_dist_perm_z(zi, z_norm[-i], crd_sim, wts_sim, nsim_, Ii, replacement=T)
    }else{
      I_perm_w_rep_sim <- permI_dist(zi, z_norm[-i], crd_sim, wts_sim, nsim_, Ii, replacement=T)
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

set.seed(100)
simulated_z_norm_10_5_0.05 <- simulated_z_norm(z_data=simulated_z_sum, alpha_sim=0.05, listw=listw_k, nsim_=9999)

# Upper/lower bounds of the two-tail test for the significance of local Moran's I described in Section 3.4.1 of our paper
sum_of_z_neigh_LB <- 5*sqrt(5)*qnorm(0.05, 0, 1, lower.tail=T)
sum_of_z_neigh_UB <- 5*sqrt(5)*qnorm(0.05, 0, 1, lower.tail=F)

simulated_z_norm_10_5_0.05 %>% 
  ggplot(aes(x=sum_of_z_norm_neigh, y=z, color=LISA_C)) +
  geom_point() +
  labs(
    x=expression(sum(paste(w[ij],z[j]), "j=1", N)),
    y=expression(z[i]),
    color=""
  ) +
  scale_color_manual(values = c("Insig"="grey60",
                                "LL"="blue",
                                "LH"="steelblue",
                                "HL"="orange",
                                "HH"="red",
                                "Obs."="black")) +
  geom_vline(xintercept=c(sum_of_z_neigh_LB, sum_of_z_neigh_UB)) +
  geom_point(data=simulated_z_sum, aes(x=lz, y=z, color="Obs.")) -> z_norm_0.05_z_plot
ggsave("Result/z norm 0.05 Exact Significance Region (Standard).pdf", z_norm_0.05_z_plot, width=15, height=10, units="cm")

simulated_z_norm_10_5_0.05_perm.i <- simulated_z_norm(z_data=simulated_z_sum, alpha_sim=0.05, listw=listw_k, nsim_=9999, perm_z=T)
simulated_z_norm_10_5_0.05_perm.i %>% 
  ggplot(aes(x=sum_of_z_norm_neigh, y=z, color=LISA_C)) +
  geom_point() +
  labs(
    x=expression(sum(paste(w[ij],z[j]), "j=1", N)),
    y=expression(z[i]),
    color=""
  ) +
  scale_color_manual(values = c("Insig"="grey60",
                                "LL"="blue",
                                "LH"="steelblue",
                                "HL"="orange",
                                "HH"="red",
                                "Obs."="black")) +
  geom_point(data=simulated_z_sum, aes(x=lz, y=z, color="Obs.")) -> z_norm_0.05_z_plot_perm.i
ggsave("Result/z norm 0.05 Exact Significance Region (Permute i).pdf", z_norm_0.05_z_plot_perm.i, width=15, height=10, units="cm")


# Results for a gamma simulation data
zero_rate <- 0.9202899; shape_hat <- 1.1415833; rate_hat<-0.1198246 # parameters of zero-inflated gamma distribution estimated from HIDTA crack cocaine seizure data in January 2020
set.seed(2024)
zero_samples <- sample(c(0,1), n_counties, replace=T, prob=c(zero_rate, 1-zero_rate))
gamma_samples <- rgamma(sum(zero_samples), shape=shape_hat, rate_hat)
zero_gamma_samples <- zero_samples
zero_gamma_samples[which(zero_samples==1)] <- gamma_samples

  # Upper/lower bounds of the two-tail test for the significance of local Moran's I described in Section 3.4.2 of our paper
gamma_mle <- fitdist(zero_gamma_samples[which(zero_gamma_samples>0)], distr = "gamma", method = "mle")
(shape_mle <- gamma_mle$estimate[1]) # shape_mle=1.2426780  / true shape: 1.1415833
(rate_mle <- gamma_mle$estimate[2]) # rate_mle=0.1249428  / true rate: 0.1198246

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

ub_binary <- function(ub, step, alpha) {
  upper_tail_probability <- sum_of_neigh_gamma(ub, shape=shape_mle, rate=rate_mle, k=5, pi=0.92)
  if (abs(upper_tail_probability - alpha) < 0.00001) {return(ub)}
  if (upper_tail_probability > alpha) {
    ub_binary(ub+step, step/2, alpha)
  }else{
    ub_binary(ub-step, step/2, alpha)
  }
}

gamma_upper_bound_df %>% ggplot(aes(z, upper_tail_probability)) + geom_line()
# picked ub as an integer value of "z" with "upper_tail_probability" close to 0.05 in the plot generated by line 503 (just above)
gamma_upper_tail <- ub_binary(ub=20, 5, 0.05)

  # Significance region for zero-inflated gamma simulation data
xx_gamma <- mean(zero_gamma_samples)
z_gamma <- zero_gamma_samples - xx_gamma
lz_gamma <- lag.listw(listw_k, z_gamma, zero.policy = zero.policy, NAOK = NAOK)
simulated_z_gamma_sum <- data.frame(z=z_gamma, lz=lz_gamma)
simulated_z_gamma_sum %>% 
  ggplot(aes(x=lz, y=z)) +
  geom_point() +
  labs(
    x=expression(sum(paste(w[ij],z[j]), "j=1", N)),
    y=expression(z[i])
  ) -> plot_z_gamma_simulation_data
ggsave("Result/scatter plot of z zero-inflated gamma simulation data.pdf", plot_z_gamma_simulation_data, width=15, height=10, units="cm")

zero_gamma_significance_region <- function(z_data, lz_gamma, alpha_sim, listw, nsim_, perm_z=F, moderate=F) {
  z_gamma <- z_data$z
  lz_gamma <- z_data$lz
  interval <- 2*sd(z_gamma)
  s2 <- sum(z_gamma^2)
  
  max_z_gamma <- max(z_gamma)
  min_z_gamma <- min(z_gamma)
  max_sum_of_z_gamma <- max(lz_gamma)
  min_sum_of_z_gamma <- min(lz_gamma)
  simulated_z_gamma <- seq(min_z_gamma, max_z_gamma, by=0.5)
  simulated_sum_of_z_gamma <- seq(min_sum_of_z_gamma, max_sum_of_z_gamma, by=0.5)
  simulated_z_gamma_pairs <- merge(simulated_z_gamma, simulated_sum_of_z_gamma) %>%
    mutate(z=x, sum_of_z_gamma_neigh=y) %>% 
    select(z, sum_of_z_gamma_neigh)
  
  if (moderate) {
    simulated_z_gamma_pairs$z_label <- cut(simulated_z_gamma_pairs$z, c(-Inf, 0, interval, Inf), labels = c("L", "M", "H"))
  } else {
    simulated_z_gamma_pairs$z_label <- cut(simulated_z_gamma_pairs$z, c(-Inf, 0, Inf), labels = c("L", "H"))
  }
  
  simulated_z_gamma_pairs$sum_of_z_gamma_neigh_label <- cut(simulated_z_gamma_pairs$sum_of_z_gamma_neigh, c(-Inf, 0, Inf), labels = lbs_sim)
  simulated_z_gamma_pairs$LISA_C <- apply(simulated_z_gamma_pairs, 1,
                                          function(x) paste(x[3], x[4], sep=""))
  simulated_z_gamma_pairs$pseudo_p <- numeric(nrow(simulated_z_gamma_pairs))
  
  crd_sim <- 5
  lww <- listw$weights
  wts_sim <- lww[[1]]
  
  simulated_z_gamma_pairs_tested <- simulated_z_gamma_pairs
  for (i in 1:nrow(simulated_z_gamma_pairs_tested)) {
    zi <- simulated_z_gamma_pairs_tested$z[i]
    sum_of_z_gamma_neigh_i <- simulated_z_gamma_pairs_tested$sum_of_z_gamma_neigh[i]
    Ii <- zi*sum_of_z_gamma_neigh_i/s2
    if (perm_z) {
      I_perm_w_rep_sim <- permI_dist_perm_z(zi, z_gamma[-i], crd_sim, wts_sim, nsim_, Ii, replacement=T)
    }else{
      I_perm_w_rep_sim <- permI_dist(zi, z_gamma[-i], crd_sim, wts_sim, nsim_, Ii, replacement=T)
    }
    R_plus <- sum(I_perm_w_rep_sim$I_perm[-(nsim_+1)] >= Ii)
    pseudo_p <- min(R_plus, nsim_-R_plus)/(nsim_+1)
    current_LISA_C <- simulated_z_gamma_pairs_tested$LISA_C[i]
    simulated_z_gamma_pairs_tested$pseudo_p[i] <- pseudo_p
    simulated_z_gamma_pairs_tested$LISA_C[i] <- ifelse(pseudo_p > alpha_sim, "Insig", current_LISA_C)
  }
  simulated_z_gamma_pairs_tested$LISA_C <- as.factor(simulated_z_gamma_pairs_tested$LISA_C)
  return(simulated_z_gamma_pairs_tested)
}

set.seed(100)
simulated_z_gamma_standard <- zero_gamma_significance_region(z_data=simulated_z_gamma_sum, alpha_sim=0.05, listw=listw_k, nsim_=9999)

simulated_z_gamma_standard %>% 
  ggplot(aes(x=sum_of_z_gamma_neigh, y=z, color=LISA_C)) +
  geom_point() +
  labs(
    x=expression(sum(paste(w[ij],z[j]), "j=1", N)),
    y=expression(z[i]),
    color=""
  ) +
  scale_color_manual(values = c("Insig"="grey60",
                                "LL"="blue",
                                "LH"="steelblue",
                                "HL"="orange",
                                "HH"="red",
                                "Obs."="black")) +
  geom_vline(xintercept=gamma_upper_tail-5*xx_gamma) +
  geom_point(data=simulated_z_gamma_sum, aes(x=lz, y=z, color="Obs.")) -> z_gamma_standard_plot
ggsave("Result/z gamma Exact Significance Region (Standard).pdf", z_gamma_standard_plot, width=15, height=10, units="cm")

simulated_z_gamma_perm.i <- zero_gamma_significance_region(z_data=simulated_z_gamma_sum, alpha_sim=0.05, listw=listw_k, nsim_=9999, perm_z=T)
simulated_z_gamma_perm.i %>% 
  ggplot(aes(x=sum_of_z_gamma_neigh, y=z, color=LISA_C)) +
  geom_point() +
  labs(
    x=expression(sum(paste(w[ij],z[j]), "j=1", N)),
    y=expression(z[i]),
    color=""
  ) +
  scale_color_manual(values = c("Insig"="grey60",
                                "LL"="blue",
                                "LH"="steelblue",
                                "HL"="orange",
                                "HH"="red",
                                "Obs."="black")) +
  geom_point(data=simulated_z_gamma_sum, aes(x=lz, y=z, color="Obs.")) -> z_gamma_perm.i_plot
ggsave("Result/z gamma Exact Significance Region (Permute i).pdf", z_gamma_perm.i_plot, width=15, height=10, units="cm")

simulated_z_gamma_moderate <- zero_gamma_significance_region(z_data=simulated_z_gamma_sum, alpha_sim=0.05, listw=listw_k, nsim_=9999, perm_z=F, moderate=T)
simulated_z_gamma_moderate %>% 
  ggplot(aes(x=sum_of_z_gamma_neigh, y=z, color=LISA_C)) +
  geom_point() +
  labs(
    x=expression(sum(paste(w[ij],z[j]), "j=1", N)),
    y=expression(z[i]),
    color=""
  ) +
  scale_color_manual(values = c("Insig"="grey60",
                                "LL"="blue",
                                "LH"="steelblue",
                                "ML"="#abd9e9",
                                "MH"="#fee090",
                                "HL"="orange",
                                "HH"="red",
                                "Obs."="black")) +
  geom_vline(xintercept=gamma_upper_tail-5*xx_gamma) +
  geom_point(data=simulated_z_gamma_sum, aes(x=lz, y=z, color="Obs.")) -> z_gamma_moderate_plot
ggsave("Result/z gamma Exact Significance Region (Moderate Label).pdf", z_gamma_moderate_plot, width=15, height=10, units="cm")

simulated_z_gamma_combined <- zero_gamma_significance_region(z_data=simulated_z_gamma_sum, alpha_sim=0.05, listw=listw_k, nsim_=9999, perm_z=T, moderate=T)
simulated_z_gamma_combined %>% 
  ggplot(aes(x=sum_of_z_gamma_neigh, y=z, color=LISA_C)) +
  geom_point() +
  labs(
    x=expression(sum(paste(w[ij],z[j]), "j=1", N)),
    y=expression(z[i]),
    color=""
  ) +
  scale_color_manual(values = c("Insig"="grey60",
                                "LL"="blue",
                                "LH"="steelblue",
                                "ML"="#abd9e9",
                                "MH"="#fee090",
                                "HL"="orange",
                                "HH"="red",
                                "Obs."="black")) +
  geom_point(data=simulated_z_gamma_sum, aes(x=lz, y=z, color="Obs.")) -> z_gamma_combined_plot
ggsave("Result/z gamma Exact Significance Region (Combined).pdf", z_gamma_combined_plot, width=15, height=10, units="cm")