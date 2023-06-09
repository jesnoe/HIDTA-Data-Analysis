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

crack.rel.sig

set.seed(100)
which(crack.rel$GEOID == 12086) # Among HH regions, Highest LISA_I (Miami-Dade, Florida, 12086)
i <- 191; zi <- z[i]; z_i <- z[-i]; crdi <- crd[i]; wtsi <- lww[[i]]; Ii <- Iis[i]
I_perm_w_rep <- permI_dist(zi, z_i, crdi, wtsi, nsim, Ii, replacement=T)
I_perm_wo_rep <- permI_dist(zi, z_i, crdi, wtsi, nsim, Ii, replacement=F)

R_plus <- sum(I_perm_w_rep$I_perm[-(nsim+1)] >= Ii)
min(R_plus, nsim-R_plus)/(nsim+1)
I_perm_w_rep %>% ggplot() +
  geom_histogram(aes(x=I_perm), bins=100) +
  geom_vline(xintercept=Ii, color="red") +
  geom_text(aes(x=Ii, y=100, label=round(Ii, 2))) +
  geom_vline(xintercept=quantile(I_perm_w_rep$I_perm, probs=c(0.05, 0.95)),
             color="blue") +
  labs(x="Permuted Moran's I", y="counts", title="With Replacement Miami-Dade, Florida in Jan 2020")

R_plus <- sum(I_perm_wo_rep$I_perm[-(nsim+1)] >= Ii)
min(R_plus, nsim-R_plus)/(nsim+1)
I_perm_wo_rep %>% ggplot() +
  geom_histogram(aes(x=I_perm), bins=100) +
  geom_vline(xintercept=Ii, color="red") +
  geom_text(aes(x=Ii, y=100, label=round(Ii, 2))) +
  geom_vline(xintercept=quantile(I_perm_wo_rep$I_perm, probs=c(0.05, 0.95)),
             color="blue") +
  labs(x="Permuted Moran's I", y="counts", title="Without Replacement Miami-Dade, Florida in Jan 2020")

which(crack.rel$GEOID == 25027) # Among HH regions, Lowest LISA_I (Worcester, Massachusetts, 25027)
i <- 591; zi <- z[i]; z_i <- z[-i]; crdi <- crd[i]; wtsi <- lww[[i]]; Ii <- Iis[i]
I_perm_w_rep <- permI_dist(zi, z_i, crdi, wtsi, nsim, Ii, replacement=T)
I_perm_wo_rep <- permI_dist(zi, z_i, crdi, wtsi, nsim, Ii, replacement=F)

R_plus <- sum(I_perm_w_rep$I_perm[-(nsim+1)] >= Ii)
min(R_plus, nsim-R_plus)/(nsim+1)
I_perm_w_rep %>% ggplot() +
  geom_histogram(aes(x=I_perm), bins=100) +
  geom_vline(xintercept=Ii, color="red") +
  geom_text(aes(x=Ii, y=100, label=round(Ii, 2))) +
  geom_vline(xintercept=quantile(I_perm_w_rep$I_perm, probs=c(0.05, 0.95)),
             color="blue") +
  labs(x="Permuted Moran's I", y="counts", title="With Replacement Worcester, MA in Jan 2020")

R_plus <- sum(I_perm_wo_rep$I_perm[-(nsim+1)] >= Ii)
min(R_plus, nsim-R_plus)/(nsim+1)
I_perm_wo_rep %>% ggplot() +
  geom_histogram(aes(x=I_perm), bins=100) +
  geom_vline(xintercept=Ii, color="red") +
  geom_text(aes(x=Ii, y=100, label=round(Ii, 2))) +
  geom_vline(xintercept=quantile(I_perm_wo_rep$I_perm, probs=c(0.05, 0.95)),
             color="blue") +
  labs(x="Permuted Moran's I", y="counts", title="Without Replacement Worcester, MA in Jan 2020")

# plot for HL
which(crack.rel$GEOID == 26121) # Among HL regions, Lowest LISA_I (Muskegon, Michigan, 26121)
set.seed(3345)
i <- 622; zi <- z[i]; z_i <- z[-i]; crdi <- crd[i]; wtsi <- lww[[i]]; Ii <- Iis[i]
I_perm_w_rep <- permI_dist(zi, z_i, crdi, wtsi, nsim, Ii, replacement=T)
I_perm_wo_rep <- permI_dist(zi, z_i, crdi, wtsi, nsim, Ii, replacement=F)

R_plus <- sum(I_perm_w_rep$I_perm[-(nsim+1)] >= Ii)
min(R_plus, nsim-R_plus)/(nsim+1)
I_perm_w_rep %>% ggplot() +
  geom_histogram(aes(x=I_perm), bins=100) +
  geom_vline(xintercept=Ii, color="red") +
  geom_text(aes(x=Ii, y=100, label=round(Ii, 2))) +
  geom_vline(xintercept=quantile(I_perm_w_rep$I_perm, probs=c(0.05, 0.95)),
             color="blue") +
  labs(x="Permuted Moran's I", y="counts", title="With Replacement Muskegon, Michigan in Jan 2020")

R_plus <- sum(I_perm_wo_rep$I_perm[-(nsim+1)] >= Ii)
min(R_plus, nsim-R_plus)/(nsim+1)
I_perm_wo_rep %>% ggplot() +
  geom_histogram(aes(x=I_perm), bins=100) +
  geom_vline(xintercept=Ii, color="red") +
  geom_text(aes(x=Ii, y=100, label=round(Ii, 2))) +
  geom_vline(xintercept=quantile(I_perm_wo_rep$I_perm, probs=c(0.05, 0.95)),
             color="blue") +
  labs(x="Permuted Moran's I", y="counts", title="Without Replacement Muskegon, Michigan in Jan 2020")


which(crack.rel$GEOID == 1073) # Among HL regions, Highest LISA_I (Jefferson, Alabama, 1073)
set.seed(3345)
i <- 29; zi <- z[i]; z_i <- z[-i]; crdi <- crd[i]; wtsi <- lww[[i]]; Ii <- Iis[i]
I_perm_w_rep <- permI_dist(zi, z_i, crdi, wtsi, nsim, Ii, replacement=T)
I_perm_wo_rep <- permI_dist(zi, z_i, crdi, wtsi, nsim, Ii, replacement=F)

R_plus <- sum(I_perm_w_rep$I_perm[-(nsim+1)] >= Ii)
min(R_plus, nsim-R_plus)/(nsim+1)
I_perm_w_rep %>% ggplot() +
  geom_histogram(aes(x=I_perm), bins=100) +
  geom_vline(xintercept=Ii, color="red") +
  geom_text(aes(x=Ii, y=100, label=round(Ii, 2))) +
  geom_vline(xintercept=quantile(I_perm_w_rep$I_perm, probs=c(0.05, 0.95)),
             color="blue") +
  labs(x="Permuted Moran's I", y="counts", title="With ReplacementJefferson, Alabama in Jan 2020")

R_plus <- sum(I_perm_wo_rep$I_perm[-(nsim+1)] >= Ii)
min(R_plus, nsim-R_plus)/(nsim+1)
I_perm_wo_rep %>% ggplot() +
  geom_histogram(aes(x=I_perm), bins=100) +
  geom_vline(xintercept=Ii, color="red") +
  geom_text(aes(x=Ii, y=100, label=round(Ii, 2))) +
  geom_vline(xintercept=quantile(I_perm_wo_rep$I_perm, probs=c(0.05, 0.95)),
             color="blue") +
  labs(x="Permuted Moran's I", y="counts", title="Without Replacement Jefferson, Alabama in Jan 2020")

# plot for Insignificant
crack.rel.insig
which(crack.rel$GEOID == 12103) # Among Insig. regions, Highest LISA_I (Pinellas, Florida, 26121)
set.seed(989)
i <- 199; zi <- z[i]; z_i <- z[-i]; crdi <- crd[i]; wtsi <- lww[[i]]; Ii <- Iis[i]
I_perm_w_rep <- permI_dist(zi, z_i, crdi, wtsi, nsim, Ii, replacement=T)
I_perm_wo_rep <- permI_dist(zi, z_i, crdi, wtsi, nsim, Ii, replacement=F)

R_plus <- sum(I_perm_w_rep$I_perm[-(nsim+1)] >= Ii)
min(R_plus, nsim-R_plus)/(nsim+1)
I_perm_w_rep %>% ggplot() +
  geom_histogram(aes(x=I_perm), bins=100) +
  geom_vline(xintercept=Ii, color="red") +
  geom_text(aes(x=Ii, y=100, label=round(Ii, 2))) +
  geom_vline(xintercept=quantile(I_perm_w_rep$I_perm, probs=c(0.05, 0.95)),
             color="blue") +
  xlim(-20, 450) +
  labs(x="Permuted Moran's I", y="counts", title="With Replacement Pinellas, Florida in Jan 2020")

R_plus <- sum(I_perm_wo_rep$I_perm[-(nsim+1)] >= Ii)
min(R_plus, nsim-R_plus)/(nsim+1)
I_perm_wo_rep %>% ggplot() +
  geom_histogram(aes(x=I_perm), bins=100) +
  geom_vline(xintercept=Ii, color="red") +
  geom_text(aes(x=Ii, y=100, label=round(Ii, 2))) +
  geom_vline(xintercept=quantile(I_perm_wo_rep$I_perm, probs=c(0.05, 0.95)),
             color="blue") +
  xlim(-20, 450) +
  labs(x="Permuted Moran's I", y="counts", title="Without Replacement Pinellas, Florida in Jan 2020")


x_Jun <- crack.rel$Jun_2020
xx_Jun <- mean(x_Jun, na.rm = NAOK)
z_Jun <- x_Jun - xx_Jun
lz_Jun <- lag.listw(listw, z_Jun, zero.policy = zero.policy, NAOK = NAOK)
ref <- tibble(LISA_CJun0=0:4, label=c("Insig", "HH", "LL", "LH", "HL"))
crack.rel %>% 
  select(state:GEOID, Jun_2020, LISA_IJun0, LISA_PJun0, LISA_CJun0) %>% 
  mutate(z=z_Jun, sum_of_z_neigh=lz_Jun) %>% 
  left_join( ref, by="LISA_CJun0") %>%
  ggplot(aes(x=sum_of_z_neigh, y=z, color=label)) +
  geom_point() +
  labs(title="Centered Seizure Counts vs. Sum of Neighbors' in Jun 2020",
       color="LISA_C")

# LISA_C in Jan 2020 via permutation without replacement
LISA_CJan0_wo_rep <- c()
for (i in 1:length(z)) {
  I_perm_wo_rep <- permI_dist(z[i], z[-i], crd[i], lww[[i]], nsim, Iis[i], replacement=F)
  R_plus <- sum(I_perm_wo_rep$I_perm[-(nsim+1)] >= Iis[i])
  pseudo_p <- min(R_plus, nsim-R_plus)/(nsim+1)
  LISA_C_i <- ifelse(pseudo_p <= 0.05, as.character(quadr_ps), "Insig")
  LISA_CJan0_wo_rep <- c(LISA_CJan0_wo_rep, LISA_C_i)
}
table(LISA_CJan0_wo_rep)

crack.rel %>% 
  select(state:GEOID, Jan_2020, LISA_IJan0, LISA_PJan0, LISA_CJan0) %>% 
  mutate(z=z, sum_of_z_neigh=lz, LISA_CJan0_wo_rep=LISA_CJan0_wo_rep) %>% 
  ggplot(aes(x=sum_of_z_neigh, y=z, color=LISA_CJan0_wo_rep)) +
  geom_point() +
  labs(title="Centered Seizure Counts vs. Sum of Neighbors' in Jan 2020",
       color="LISA_C")

set.seed(989)
which(crack.rel$GEOID == 6109) # Among Insig. regions, one of zero-seizure counties (Tuolumne, California, 6109)
i <- 136; zi <- z[i]; z_i <- z[-i]; crdi <- crd[i]; wtsi <- lww[[i]]; Ii <- Iis[i]
I_perm_w_rep <- permI_dist(zi, z_i, crdi, wtsi, nsim, Ii, replacement=T)
I_perm_wo_rep <- permI_dist(zi, z_i, crdi, wtsi, nsim, Ii, replacement=F)

R_plus <- sum(I_perm_w_rep$I_perm[-(nsim+1)] >= Ii)
min(R_plus, nsim-R_plus)/(nsim+1)
I_perm_w_rep %>% ggplot() +
  geom_histogram(aes(x=I_perm), bins=100) +
  geom_vline(xintercept=Ii, color="red") +
  geom_text(aes(x=Ii, y=100, label=round(Ii, 2))) +
  geom_vline(xintercept=quantile(I_perm_w_rep$I_perm, probs=c(0.05, 0.95)),
             color="blue") +
  labs(x="Permuted Moran's I", y="counts", title="With Replacement Tuolumne, California in Jan 2020")

R_plus <- sum(I_perm_wo_rep$I_perm[-(nsim+1)] >= Ii)
min(R_plus, nsim-R_plus)/(nsim+1)
I_perm_wo_rep %>% ggplot() +
  geom_histogram(aes(x=I_perm), bins=100) +
  geom_vline(xintercept=Ii, color="red") +
  geom_text(aes(x=Ii, y=100, label=round(Ii, 2))) +
  geom_vline(xintercept=quantile(I_perm_wo_rep$I_perm, probs=c(0.05, 0.95)),
             color="blue") +
  labs(x="Permuted Moran's I", y="counts", title="Without Replacement Tuolumne, California in Jan 2020")

## plot of z vs. sum of neighbors' z
ref <- tibble(LISA_CJan0=0:4, label=c("Insig", "HH", "LL", "LH", "HL"))
crack.rel %>% 
  select(state:GEOID, Jan_2020, LISA_IJan0, LISA_PJan0, LISA_CJan0) %>% 
  mutate(z=z, sum_of_z_neigh=lz) %>% 
  left_join( ref, by="LISA_CJan0") %>%
  ggplot(aes(x=sum_of_z_neigh, y=z, color=label)) +
  geom_point() +
  labs(title="Centered Seizure Counts vs. Sum of Neighbors' in Jan 2020") +
  scale_color_manual(values = c("Insig"="grey60",
                               "LL"="blue",
                               "LH"="steelblue",
                               "HL"="orange",
                               "HH"="red"))

freq_z_space <- crack.rel %>% 
  select(state:GEOID, Jan_2020, LISA_IJan0, LISA_PJan0, LISA_CJan0) %>% 
  mutate(z=z, sum_of_z_neigh=lz) %>% 
  left_join( ref, by="LISA_CJan0") %>%
  group_by(z, sum_of_z_neigh) %>% 
  summarise(num_of_cases=n()) %>%
  arrange(num_of_cases)
freq_z_space$cum_dist <- round(cumsum(freq_z_space$num_of_cases)/1518, 3)
freq_z_space %>% as.data.frame

alpha <- 0.01
freq_z_space$significance <- ifelse(freq_z_space$cum_dist <= alpha, "sig", "insig")
freq_z_space %>% 
  ggplot(aes(x=sum_of_z_neigh, y=z, color=significance)) +
  geom_point() +
  labs(title="Centered Seizure Counts vs. Sum of Neighbors' in Jan 2020")

alpha <- 0.03
freq_z_space$significance <- ifelse(freq_z_space$cum_dist <= alpha, "sig", "insig")
freq_z_space %>% 
  ggplot(aes(x=sum_of_z_neigh, y=z, color=significance)) +
  geom_point() +
  labs(title="Centered Seizure Counts vs. Sum of Neighbors' in Jan 2020")

alpha <- 0.05
freq_z_space$significance <- ifelse(freq_z_space$cum_dist <= alpha, "sig", "insig")
freq_z_space %>% 
  ggplot(aes(x=sum_of_z_neigh, y=z, color=significance)) +
  geom_point() +
  labs(title="Centered Seizure Counts vs. Sum of Neighbors' in Jan 2020")

freq_z_space %>% 
  ggplot(aes(x=sum_of_z_neigh, y=z, size=num_of_cases)) +
  geom_point() +
  labs(title="Centered Seizure Counts vs. Sum of Neighbors' in Jan 2020")


crack.rel %>% 
  select(state:GEOID, Jan_2020, LISA_IJan0, LISA_PJan0, LISA_CJan0) %>% 
  mutate(z=z, sum_of_z_neigh=lz) %>% 
  left_join( ref, by="LISA_CJan0") %>% 
  filter(z < -0.758 & sum_of_z_neigh %in% c(14.2, 16.2))

crack.rel %>% 
  select(state:GEOID, Jan_2020, LISA_IJan0, LISA_PJan0, LISA_CJan0) %>% 
  mutate(z=z, sum_of_z_neigh=lz) %>% 
  left_join( ref, by="LISA_CJan0") %>%
  ggplot(aes(x=sum_of_z_neigh, y=z, color=label)) +
  geom_point(position=position_jitter(w=1, h=0)) +
  labs(title="Centered Seizure Counts vs. Sum of Neighbors' in Jan 2020")

crack.rel %>% 
  select(state:GEOID, Jan_2020, LISA_IJan0, LISA_PJan0, LISA_CJan0) %>% 
  mutate(z=z, sum_of_z_neigh=lz) %>% 
  left_join( ref, by="LISA_CJan0") %>% 
  filter(label == "HL") %>% 
  as.data.frame

crack.rel %>% 
  select(state:GEOID, Jan_2020, LISA_IJan0, LISA_PJan0, LISA_CJan0) %>% 
  mutate(z=z, sum_of_z_neigh=lz) %>% 
  left_join( ref, by="LISA_CJan0") %>% 
  filter(label == "LH") %>% 
  as.data.frame

ref <- tibble(LISA_CJan0=0:4, label=c("Insig", "HH", "LL", "LH", "HL"))
crack.rel %>% 
  select(state:GEOID, Jan_2020, LISA_IJan0, LISA_PJan0, LISA_CJan0) %>% 
  mutate(z=z, sum_of_z_neigh=lz) %>% 
  left_join( ref, by="LISA_CJan0") %>%
  ggplot(aes(x=sum_of_z_neigh, y=z, color=label)) +
  geom_point() +
  labs(title="Centered Seizure Counts vs. Sum of Neighbors' in Jan 2020")

freq_z_space <- crack.rel %>% 
  select(state:GEOID, Jan_2020, LISA_IJan0, LISA_PJan0, LISA_CJan0) %>% 
  mutate(z=z, sum_of_z_neigh=lz) %>% 
  left_join( ref, by="LISA_CJan0") %>%
  group_by(z, sum_of_z_neigh) %>% 
  summarise(num_of_cases=n()) %>%
  arrange(num_of_cases)
freq_z_space$cum_dist <- round(cumsum(freq_z_space$num_of_cases)/1518, 3)
freq_z_space %>% as.data.frame

alpha <- 0.01
freq_z_space$significance <- ifelse(freq_z_space$cum_dist <= alpha, "sig", "insig")
freq_z_space %>% 
  ggplot(aes(x=sum_of_z_neigh, y=z, color=significance)) +
  geom_point() +
  labs(title="Centered Seizure Counts vs. Sum of Neighbors' in Jan 2020")

alpha <- 0.03
freq_z_space$significance <- ifelse(freq_z_space$cum_dist <= alpha, "sig", "insig")
freq_z_space %>% 
  ggplot(aes(x=sum_of_z_neigh, y=z, color=significance)) +
  geom_point() +
  labs(title="Centered Seizure Counts vs. Sum of Neighbors' in Jan 2020")

alpha <- 0.05
freq_z_space$significance <- ifelse(freq_z_space$cum_dist <= alpha, "sig", "insig")
freq_z_space %>% 
  ggplot(aes(x=sum_of_z_neigh, y=z, color=significance)) +
  geom_point() +
  labs(title="Centered Seizure Counts vs. Sum of Neighbors' in Jan 2020")

freq_z_space %>% 
  ggplot(aes(x=sum_of_z_neigh, y=z, size=num_of_cases)) +
  geom_point() +
  labs(title="Centered Seizure Counts vs. Sum of Neighbors' in Jan 2020")


crack.rel %>% 
  select(state:GEOID, Jan_2020, LISA_IJan0, LISA_PJan0, LISA_CJan0) %>% 
  mutate(z=z, sum_of_z_neigh=lz) %>% 
  left_join( ref, by="LISA_CJan0") %>% 
  filter(z < -0.758 & sum_of_z_neigh %in% c(14.2, 16.2))

crack.rel %>% 
  select(state:GEOID, Jan_2020, LISA_IJan0, LISA_PJan0, LISA_CJan0) %>% 
  mutate(z=z, sum_of_z_neigh=lz) %>% 
  left_join( ref, by="LISA_CJan0") %>%
  ggplot(aes(x=sum_of_z_neigh, y=z, color=label)) +
  geom_point(position=position_jitter(w=1, h=0)) +
  labs(title="Centered Seizure Counts vs. Sum of Neighbors' in Jan 2020")

crack.rel %>% 
  select(state:GEOID, Jan_2020, LISA_IJan0, LISA_PJan0, LISA_CJan0) %>% 
  mutate(z=z, sum_of_z_neigh=lz) %>% 
  left_join( ref, by="LISA_CJan0") %>% 
  filter(label == "HL") %>% 
  as.data.frame

## plot of simulated z vs. sum of neighbors' z
lbs_sim <- c("L", "H")
max_z <- max(z)
min_z <- min(z)

nb_crack_k <- knn2nb(knearneigh(coords.crack, k=10), row.names=GEOIDS.crack)
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
  I_perm_w_rep_sim <- permI_dist(zi, z, crd_sim, wts_sim, nsim, Ii, replacement=T)
  R_plus <- sum(I_perm_w_rep_sim$I_perm[-(nsim+1)] >= Ii)
  pseudo_p <- min(R_plus, nsim-R_plus)/(nsim+1)
  current_LISA_C <- simulated_z_pairs_tested$LISA_C[i]
  simulated_z_pairs_tested$pseudo_p[i] <- pseudo_p
  simulated_z_pairs_tested$LISA_C[i] <- ifelse(pseudo_p > alpha_sim, "Insig", current_LISA_C)
}
simulated_z_pairs_tested$LISA_C <- as.factor(simulated_z_pairs_tested$LISA_C)

observed_z_sum <- data.frame(z=z, lz=lz)
simulated_z_pairs_tested %>% 
  ggplot(aes(x=sum_of_z_neigh, y=z, color=LISA_C)) +
  geom_point() +
  labs(title="Centered Seizure Counts vs. Sum of Neighbors' in Jan 2020 (k=5)") +
  scale_color_manual(values = c("Insig"="grey60",
                                "LL"="blue",
                                "LH"="steelblue",
                                "HL"="orange",
                                "HH"="red",
                                "Obs."="black")) +
  geom_point(data=observed_z_sum, aes(x=lz, y=z, color="Obs."))

simulated_z_pairs_tested %>% 
  ggplot(aes(x=sum_of_z_neigh, y=z, color=LISA_C)) +
  geom_point() +
  xlim(-10, 170) +
  labs(title="Centered Seizure Counts vs. Sum of Neighbors' in Jan 2020 (k=10)") +
  scale_color_manual(values = c("Insig"="grey60",
                                "LL"="blue",
                                "LH"="steelblue",
                                "HL"="orange",
                                "HH"="red"))

simulated_z_pairs_tested %>% filter(pseudo_p > 0.05 & sum_of_z_neigh < -3.5 & z > 50)
which(simulated_z_pairs_tested$pseudo_p > 0.05 & simulated_z_pairs_tested$sum_of_z_neigh < -3.5 & simulated_z_pairs_tested$z > 50)
sum(x>0)/length(x)
1-(1-0.08)^5 # 0.3409185

## plot of simulated normal z vs. sum of neighbors' z
summary(x) # min=0, max=92
n_counties <- length(x)

simulated_z_norm <- function(x_mean, x_sd, n_points, alpha_sim, listw, perm_z=F) {
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
      I_perm_w_rep_sim <- permI_dist_perm_z(zi, z_norm, crd_sim, wts_sim, nsim, Ii, replacement=T)
    }else{
      I_perm_w_rep_sim <- permI_dist(zi, z_norm, crd_sim, wts_sim, nsim, Ii, replacement=T)
    }
    R_plus <- sum(I_perm_w_rep_sim$I_perm[-(nsim+1)] >= Ii)
    pseudo_p <- min(R_plus, nsim-R_plus)/(nsim+1)
    current_LISA_C <- simulated_z_norm_pairs_tested$LISA_C[i]
    simulated_z_norm_pairs_tested$pseudo_p[i] <- pseudo_p
    simulated_z_norm_pairs_tested$LISA_C[i] <- ifelse(pseudo_p > alpha_sim, "Insig", current_LISA_C)
  }
  simulated_z_norm_pairs_tested$LISA_C <- as.factor(simulated_z_norm_pairs_tested$LISA_C)
  return(simulated_z_norm_pairs_tested)
}

nb_crack_k <- knn2nb(knearneigh(coords.crack, k=5), row.names=GEOIDS.crack)
listw_k <- nb2listw(nb_crack_k, style="B")
set.seed(5839)
simulated_z_norm_10_5_0.05 <- simulated_z_norm(x_mean=10, x_sd=5, n_points=n_counties, alpha_sim=0.05, listw=listw_k)
set.seed(5839)
simulated_z_norm_10_5_0.10 <- simulated_z_norm(x_mean=10, x_sd=5, n_points=n_counties, alpha_sim=0.10, listw=listw_k)

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
  labs(title="Centered Seizure Counts vs. Sum of Neighbors' in Jan 2020 (k=5)") +
  scale_color_manual(values = c("Insig"="grey60",
                                "LL"="blue",
                                "LH"="steelblue",
                                "HL"="orange",
                                "HH"="red"))


nb_crack_k <- knn2nb(knearneigh(coords.crack, k=5), row.names=GEOIDS.crack)
listw_k <- nb2listw(nb_crack_k, style="B")
set.seed(5839)
simulated_z_norm_10_5_0.05 <- simulated_z_norm(x_mean=10, x_sd=5, n_points=n_counties, alpha_sim=0.05, listw=listw_k, perm_z=T)

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
  scale_color_manual(values = c("Insig"="grey60",
                                "LL"="blue",
                                "LH"="steelblue",
                                "HL"="orange",
                                "HH"="red"))

