library(readxl)
library(tidyverse)
library(fpp2)

# Diffusion Pathways
LISA <- read.csv("CokeWeight7313.csv") %>% as.tibble
names(LISA)
names(LISA)[5] <- "state"

grep("y1973", names(LISA)) # 14
grep("LISA_I", names(LISA))[1] # 55
LISA[,14:54]

lag.t <- 1
LISA.ts <- list()
diffusion <- list()
diffusion_summary <- list()
state.summary <- list()
for (state_name in unique(LISA$state)) {
  LISA_state <- LISA %>% filter(state == state_name)
  
  LISA.df <- data.frame(LISA_I=t(LISA_state[,grep("LISA_I", names(LISA_state))]),
                        LISA_C=t(LISA_state[,grep("LISA_C", names(LISA_state))]),
                        LISA_P=t(LISA_state[,grep("LISA_P", names(LISA_state))]),
                        LISA_I_lag=lag(t(LISA_state[,grep("LISA_I", names(LISA_state))]),lag.t),
                        LISA_C_lag=lag(t(LISA_state[,grep("LISA_C", names(LISA_state))]),lag.t),
                        LISA_P_lag=lag(t(LISA_state[,grep("LISA_P", names(LISA_state))]),lag.t))
  row.names(LISA.df) <- NULL
  LISA.ts[[state_name]] <- LISA.df
  
  
  diffusion[[state_name]] <- list()
  LISA_C_i <- LISA.ts[[state_name]]$LISA_C
  LISA_C_lag_i <- LISA.ts[[state_name]]$LISA_C_lag
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
  # diffusion_ratio <- data.frame(
  #   state=state_name,
  #   county=county,
  #   geiod=LISA_state$GEOID[LISA_state$NAMELSAD == county],
  #   n_total=N.tot,
  #   stationary=n_stationary/N.tot
  # )
  # diffusion_ratio$change <- 1-diffusion_ratio$stationary
  # n_increase <- sum(N[c(1,3), c(2,4)])
  # n_decrease <- sum(N[c(2,4), c(1,3)])
  # diffusion_ratio$increase <- n_increase/(N.tot-n_stationary)
  # diffusion_ratio$decrease <- n_decrease/(N.tot-n_stationary)
  # diffusion_ratio$non_adj_inc <- sum(N[1,c(2,4)])/n_increase
  # diffusion_ratio$adj_inc <- sum(N[3,c(2,4)])/n_increase
  # diffusion_ratio$non_adj_dec <- sum(N[4,c(1,3)])/n_decrease
  # diffusion_ratio$adj_dec <- sum(N[2,c(1,3)])/n_decrease
  # diffusion_ratio$local_High <- sum(LISA_C_i %in% c(1,4))/48
  # diffusion_ratio$local_Low <- sum(LISA_C_i %in% c(2,3))/48
  # diffusion_ratio$local_Ins <- sum(LISA_C_i==0)/48
  
  diffusion_count <- data.frame(
    state=state_name,
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
  diffusion[[state_name]] <- result
  
  # diffusion_ratio_state <- data.frame(diffusion[[state_name]][[counties[1]]]$ratio)
  # for (county in counties[-1]) {
  #   county_ratio <- data.frame(diffusion[[state_name]][[county]]$ratio)
  #   diffusion_ratio_state <- rbind(diffusion_ratio_state, county_ratio)
  # }
  # diffusion_summary[[state_name]] <- list()
  # diffusion_summary[[state_name]][["ratio"]] <- diffusion_ratio_state

}
diffusion_summary <- diffusion$Alabama$count
for (state_name in unique(LISA$state)[-1]) {
  diffusion_summary <- rbind(diffusion_summary, diffusion[[state_name]][["count"]])
}
diffusion_summary
diffusion_summary_cocaine <- diffusion_summary


diffusion_summary <- diffusion_summary_cocaine
diffusion_summary$California$count %>% select(-state, -geiod) %>% arrange(desc(change))
diffusion_summary$Florida$count %>% select(-state, -geiod) %>% arrange(desc(change))
diffusion_summary$Texas$count %>% select(-state, -geiod) %>%arrange(desc(change))

# write.csv(diffusion_summary, "STRIDE cocaine diffusion summary 4x4 (count).csv", row.names=F)

diffusion_summary_count %>% arrange(desc(non_adj_inc))
diffusion_summary_count %>% arrange(desc(local_High))



# Monthly changes in national level diffusion factors
# 0=insignificant, 1=HH, 2=LL, 3=LH, 4=HL
start <- grep("LISA_C", names(LISA))[1] # 73
end <- tail(grep("LISA_C", names(LISA)), 1) # 214
LISA_C_index <- seq(start, end, by=3)
adj_inc_monthly <- c()
non_adj_inc_monthly <- c()
adj_dec_monthly <- c()
non_adj_dec_monthly <- c()
for (i in 1:47) {
  LISA_C_t0 <- LISA[, LISA_C_index[i]] %>% as.matrix %>% as.vector
  LISA_C_t1 <- LISA[, LISA_C_index[i]+3] %>% as.matrix %>% as.vector
  adj_inc_m <- ifelse(LISA_C_t0 == 3 & LISA_C_t1 %in% c(1, 4), 1, 0) %>% sum
  non_adj_inc_m <- ifelse(LISA_C_t0 == 2 & LISA_C_t1 %in% c(1, 4), 1, 0) %>% sum
  adj_dec_m <- ifelse(LISA_C_t0 == 4 & LISA_C_t1 %in% c(2, 3), 1, 0) %>% sum
  non_adj_dec_m <- ifelse(LISA_C_t0 == 1 & LISA_C_t1 %in% c(2, 3), 1, 0) %>% sum
  adj_inc_monthly <- c(adj_inc_monthly, adj_inc_m)
  non_adj_inc_monthly <- c(non_adj_inc_monthly, non_adj_inc_m)
  adj_dec_monthly <- c(adj_dec_monthly, adj_dec_m)
  non_adj_dec_monthly <- c(non_adj_dec_monthly, non_adj_dec_m)
}

adj_inc_monthly
non_adj_inc_monthly
adj_dec_monthly
non_adj_dec_monthly

diffusion_monthly <- data.frame(adj_inc=adj_inc_monthly,
                                non_adj_inc=non_adj_inc_monthly,
                                adj_dec=adj_dec_monthly,
                                non_adj_dec=non_adj_dec_monthly)
diffusion_monthly_cocaine <- diffusion_monthly

##
diffusion_monthly <- diffusion_monthly_cocaine
diffusion_monthly
autoplot(ts(diffusion_monthly$adj_inc, start=c(2018, 2), frequency=12),
         main="County Data Adj & Non-adj Increases (Cocaine)",
         xlab="Year",
         ylab="Counts",
         series="adj_inc") + 
  autolayer(ts(diffusion_monthly$non_adj_inc, start=c(2018, 2), frequency=12), series="non_adj_inc")

autoplot(ts(diffusion_monthly$adj_dec, start=c(2018, 2), frequency=12),
         main="County Data Adj & Non-adj Decreases (Cocaine)",
         xlab="Year",
         ylab="Counts",
         series="adj_dec") +
  autolayer(ts(diffusion_monthly$non_adj_dec, start=c(2018, 2), frequency=12), series="non_adj_dec")

autoplot(ts(diffusion_monthly$adj_inc + diffusion_monthly$non_adj_inc,
            start=c(2018, 2), frequency=12),
         main="County Data Increases & Decreases (Cocaine)",
         xlab="Year",
         ylab="Counts",
         series="increase") +
  autolayer(ts(diffusion_monthly$adj_dec + diffusion_monthly$non_adj_dec, start=c(2018, 2), frequency=12), series="decrease")

autoplot(ts(apply(diffusion_monthly, 1, sum),
            start=c(2018, 2), frequency=12),
         main="County Data Significant Changes (Cocaine)",
         xlab="Year",
         ylab="Counts")

plot(diffusion_monthly$adj_inc, diffusion_monthly$non_adj_inc, xlab="adj_inc", ylab="non_adj_inc", main="cocaine", pch=19)
plot(diffusion_monthly$adj_dec, diffusion_monthly$non_adj_dec, xlab="adj_dec", ylab="non_adj_dec", main="cocaine", pch=19)
plot(diffusion_monthly$adj_inc+diffusion_monthly$non_adj_inc, diffusion_monthly$adj_dec+diffusion_monthly$non_adj_dec,
     main="cocaine", xlab="increases", ylab="decreases", pch=19)


{
par(mfrow=c(2,2))
plot(diffusion_monthly$adj_inc[1:11], diffusion_monthly$non_adj_inc[1:11], xlab="adj_inc", ylab="non_adj_inc", main="cocaine 2018", pch=19)
plot(diffusion_monthly$adj_inc[12:23], diffusion_monthly$non_adj_inc[12:23], xlab="adj_inc", ylab="non_adj_inc", main="cocaine 2019", pch=19)
plot(diffusion_monthly$adj_inc[24:35], diffusion_monthly$non_adj_inc[24:35], xlab="adj_inc", ylab="non_adj_inc", main="cocaine 2020", pch=19)
plot(diffusion_monthly$adj_inc[36:47], diffusion_monthly$non_adj_inc[36:47], xlab="adj_inc", ylab="non_adj_inc", main="cocaine 2021", pch=19)
par(mfrow=c(1,1))
}

{
par(mfrow=c(2,2))
plot(diffusion_monthly$adj_dec[1:11], diffusion_monthly$non_adj_dec[1:11], xlab="adj_dec", ylab="non_adj_dec", main="cocaine 2018", pch=19)
plot(diffusion_monthly$adj_dec[12:23], diffusion_monthly$non_adj_dec[12:23], xlab="adj_dec", ylab="non_adj_dec", main="cocaine 2019", pch=19)
plot(diffusion_monthly$adj_dec[24:35], diffusion_monthly$non_adj_dec[24:35], xlab="adj_dec", ylab="non_adj_dec", main="cocaine 2020", pch=19)
plot(diffusion_monthly$adj_dec[36:47], diffusion_monthly$non_adj_dec[36:47], xlab="adj_dec", ylab="non_adj_dec", main="cocaine 2021", pch=19)
par(mfrow=c(1,1))
}

{
par(mfrow=c(2,2))
plot(diffusion_monthly$adj_inc[1:11]+diffusion_monthly$non_adj_inc[1:11], diffusion_monthly$adj_dec[1:11]+diffusion_monthly$non_adj_dec[1:11],
     xlab="increases", ylab="decreases", main="cocaine 2018", pch=19)
plot(diffusion_monthly$adj_inc[12:23]+diffusion_monthly$non_adj_inc[12:23], diffusion_monthly$adj_dec[12:23]+diffusion_monthly$non_adj_dec[12:23],
     xlab="increases", ylab="decreases", main="cocaine 2019", pch=19)
plot(diffusion_monthly$adj_inc[24:35]+diffusion_monthly$non_adj_inc[24:35], diffusion_monthly$adj_dec[24:35]+diffusion_monthly$non_adj_dec[24:35],
     xlab="increases", ylab="decreases", main="cocaine 2020", pch=19)
plot(diffusion_monthly$adj_inc[36:47]+diffusion_monthly$non_adj_inc[36:47], diffusion_monthly$adj_dec[36:47]+diffusion_monthly$non_adj_dec[36:47],
     xlab="increases", ylab="decreases", main="cocaine 2021", pch=19)
par(mfrow=c(1,1))
}


plot(diffusion_monthly$non_adj_inc, diffusion_monthly$adj_dec, xlab="non_adj_inc", ylab="adj_dec", main="cocaine", pch=19)
plot(diffusion_monthly$non_adj_dec, diffusion_monthly$adj_inc, xlab="non_adj_dec", ylab="adj_inc", main="cocaine", pch=19)

pairs(diffusion_summary_count[,9:12], pch=19)
# numbers of counties are different for diffusion_summary_count and LISA

start <- grep("Jan_2018", names(LISA))[1] # 21
end <- grep("Dec_2021", names(LISA))[1] # 68
seizure.weights.cocaine <- ts(apply(LISA[,start:end], 2, sum), start=c(2018,1), frequency=12)
seizure.weights <- seizure.weights.cocaine
autoplot(seizure.weights, xlab="Year", ylab="Seizure Weights (kg)", main="National Seizure Weights (Cocaine)")

{
par(mfrow=c(2,2))
plot(as.vector(seizure.weights)[-48], diffusion_monthly$adj_inc,
     main="cocaine", ylab="adj_inc at t+1", xlab="Seizure Weights at t (kg)", pch=19)
plot(as.vector(seizure.weights)[-48], diffusion_monthly$non_adj_inc,
     ylab="non_adj_inc at t+1", xlab="", pch=19)

plot(as.vector(seizure.weights)[-48], diffusion_monthly$adj_dec,
     ylab="adj_dec at t+1", xlab="Seizure Weights at t (kg)", pch=19)
plot(as.vector(seizure.weights)[-48], diffusion_monthly$non_adj_dec,
     ylab="non_adj_dec at t+1", xlab="", pch=19)
par(mfrow=c(1,1))
}

{
par(mfrow=c(2,2))
plot(log(as.vector(seizure.weights)[-48]), diffusion_monthly$adj_inc,
     main="cocaine", ylab="adj_inc at t+1", xlab="Seizure Weights at t (log kg)", pch=19)
plot(log(as.vector(seizure.weights)[-48]), diffusion_monthly$non_adj_inc,
     ylab="non_adj_inc at t+1", xlab="", pch=19)

plot(log(as.vector(seizure.weights)[-48]), diffusion_monthly$adj_dec,
     ylab="adj_dec at t+1", xlab="Seizure Weights at t (log kg)", pch=19)
plot(log(as.vector(seizure.weights)[-48]), diffusion_monthly$non_adj_dec,
     ylab="non_adj_dec at t+1", xlab="", pch=19)
par(mfrow=c(1,1))
}

# The impact of seizure weights in a state
# LISA <- read.csv("CocaineLISAResults.csv") %>% as.tibble # GEOID is county fips
# state.fips <- read.csv("us-state-ansi-fips.csv")
# names(state.fips)[1] <- "state"
# names(state.fips)[2] <- "STATEFP"
# LISA <- merge(LISA, state.fips[,1:2], by="STATEFP") %>% as.tibble
# LISA <- LISA[,c(1, 212, 2:211)]

state <- "California"
start <- grep("LISA_C", names(LISA))[1] # 70
end <- tail(grep("LISA_C", names(LISA)), 1) # 211
LISA_C_index <- seq(start, end, by=3)
adj_inc_monthly <- c()
non_adj_inc_monthly <- c()
adj_dec_monthly <- c()
non_adj_dec_monthly <- c()
for (i in 1:47) {
  LISA_C_t0 <- LISA[LISA$state == state, LISA_C_index[i]] %>% as.matrix %>% as.vector
  LISA_C_t1 <- LISA[LISA$state == state, LISA_C_index[i]+3] %>% as.matrix %>% as.vector
  adj_inc_m <- ifelse(LISA_C_t0 == 3 & LISA_C_t1 %in% c(1, 4), 1, 0) %>% sum
  non_adj_inc_m <- ifelse(LISA_C_t0 == 2 & LISA_C_t1 %in% c(1, 4), 1, 0) %>% sum
  adj_dec_m <- ifelse(LISA_C_t0 == 4 & LISA_C_t1 %in% c(2, 3), 1, 0) %>% sum
  non_adj_dec_m <- ifelse(LISA_C_t0 == 1 & LISA_C_t1 %in% c(2, 3), 1, 0) %>% sum
  adj_inc_monthly <- c(adj_inc_monthly, adj_inc_m)
  non_adj_inc_monthly <- c(non_adj_inc_monthly, non_adj_inc_m)
  adj_dec_monthly <- c(adj_dec_monthly, adj_dec_m)
  non_adj_dec_monthly <- c(non_adj_dec_monthly, non_adj_dec_m)
}

diffusion_state_monthly <- data.frame(adj_inc=adj_inc_monthly,
                                      non_adj_inc=non_adj_inc_monthly,
                                      adj_dec=adj_dec_monthly,
                                      non_adj_dec=non_adj_dec_monthly)
diffusion_state_monthly_cocaine <- diffusion_state_monthly
seizure.start <- grep("Jan_2018", names(LISA))[1] # 21
seizure.end <- grep("Dec_2021", names(LISA))[1] # 68
seizure.weights.state.cocaine <- apply(LISA[LISA$state == state, seizure.start:seizure.end], 2, sum)[-48]

##
diffusion_state_monthly <- diffusion_state_monthly_cocaine
autoplot(ts(diffusion_state_monthly$adj_inc, start=c(2018, 2), frequency=12),
         main="County Data Adj & Non-adj Increases (Cocaine)",
         xlab="Year",
         ylab="Counts",
         series="adj_inc") + 
  autolayer(ts(diffusion_state_monthly$non_adj_inc, start=c(2018, 2), frequency=12), series="non_adj_inc")

autoplot(ts(diffusion_state_monthly$adj_dec, start=c(2018, 2), frequency=12),
         main="County Data Adj & Non-adj Decreases (Cocaine)",
         xlab="Year",
         ylab="Counts",
         series="adj_dec") +
  autolayer(ts(diffusion_state_monthly$non_adj_dec, start=c(2018, 2), frequency=12), series="non_adj_dec")

autoplot(ts(diffusion_state_monthly$adj_inc + diffusion_state_monthly$non_adj_inc,
            start=c(2018, 2), frequency=12),
         main="County Data Increases & Decreases (Cocaine)",
         xlab="Year",
         ylab="Counts",
         series="increase") +
  autolayer(ts(diffusion_state_monthly$adj_dec + diffusion_state_monthly$non_adj_dec, start=c(2018, 2), frequency=12), series="decrease")

autoplot(ts(apply(diffusion_state_monthly, 1, sum),
            start=c(2018, 2), frequency=12),
         main="County Data Significant Changes (Cocaine)",
         xlab="Year",
         ylab="Counts")




{
par(mfrow=c(2,2))
plot(as.vector(seizure.weights.state.cocaine)[-48], diffusion_state_monthly$adj_inc, 
     main="cocaine", ylab="adj_inc at t+1", xlab="Seizure Weights at t (kg)", pch=19)
plot(as.vector(seizure.weights.state.cocaine)[-48], diffusion_state_monthly$non_adj_inc, 
     ylab="non_adj_inc at t+1", xlab="", pch=19)

plot(as.vector(seizure.weights.state.cocaine)[-48], diffusion_state_monthly$adj_dec, 
     ylab="adj_dec at t+1", xlab="Seizure Weights at t (kg)", pch=19)
plot(as.vector(seizure.weights.state.cocaine)[-48], diffusion_state_monthly$non_adj_dec, 
     ylab="non_adj_dec at t+1", xlab="", pch=19)
par(mfrow=c(1,1))
}

{
par(mfrow=c(2,2))
plot(log(as.vector(seizure.weights.state.cocaine)[-48]), diffusion_state_monthly$adj_inc,
     main="cocaine", ylab="adj_inc at t+1", xlab="Seizure Weights at t (log kg)", pch=19)
plot(log(as.vector(seizure.weights.state.cocaine)[-48]), diffusion_state_monthly$non_adj_inc,
     ylab="non_adj_inc at t+1", xlab="", pch=19)

plot(log(as.vector(seizure.weights.state.cocaine)[-48]), diffusion_state_monthly$adj_dec,
     ylab="adj_dec at t+1", xlab="Seizure Weights at t (log kg)", pch=19)
plot(log(as.vector(seizure.weights.state.cocaine)[-48]), diffusion_state_monthly$non_adj_dec,
     ylab="non_adj_dec at t+1", xlab="", pch=19)
par(mfrow=c(1,1))
}


# Linear regression coefficients
start <- grep("LISA_C", names(LISA))[1] # 70
end <- tail(grep("LISA_C", names(LISA)), 1) # 211
LISA_C_index <- seq(start, end, by=3)
seizure.start <- grep("Jan_2018", names(LISA))[1] # 21
seizure.end <- grep("Dec_2021", names(LISA))[1] # 68

state.names <- unique(LISA$state)
adj_inc_reg_cocaine <- data.frame(state=NA, coef=NA, p_value=NA, sig.=NA, R_sq=NA)
nadj_inc_reg_cocaine <- data.frame(state=NA, coef=NA, p_value=NA, sig.=NA, R_sq=NA)
adj_dec_reg_cocaine <- data.frame(state=NA, coef=NA, p_value=NA, sig.=NA, R_sq=NA)
nadj_dec_reg_cocaine <- data.frame(state=NA, coef=NA, p_value=NA, sig.=NA, R_sq=NA)
for (i in 1:length(state.names)) {
  state <- state.names[i]
  adj_inc_monthly <- c()
  non_adj_inc_monthly <- c()
  adj_dec_monthly <- c()
  non_adj_dec_monthly <- c()
  for (j in 1:47) {
    LISA_C_t0 <- LISA[LISA$state == state, LISA_C_index[j]] %>% as.matrix %>% as.vector
    LISA_C_t1 <- LISA[LISA$state == state, LISA_C_index[j]+3] %>% as.matrix %>% as.vector
    n_stationary_insig <- sum(LISA_C_t0 + LISA_C_t1 == 0)
    adj_inc_m <- ifelse(LISA_C_t0 == 3 & LISA_C_t1 %in% c(1, 4), 1, 0) %>% sum
    non_adj_inc_m <- ifelse(LISA_C_t0 == 2 & LISA_C_t1 %in% c(1, 4), 1, 0) %>% sum
    adj_dec_m <- ifelse(LISA_C_t0 == 4 & LISA_C_t1 %in% c(2, 3), 1, 0) %>% sum
    non_adj_dec_m <- ifelse(LISA_C_t0 == 1 & LISA_C_t1 %in% c(2, 3), 1, 0) %>% sum
    adj_inc_monthly <- c(adj_inc_monthly, adj_inc_m)
    non_adj_inc_monthly <- c(non_adj_inc_monthly, non_adj_inc_m)
    adj_dec_monthly <- c(adj_dec_monthly, adj_dec_m)
    non_adj_dec_monthly <- c(non_adj_dec_monthly, non_adj_dec_m)
  }
  
  seizure.weights.state.cocaine <- apply(LISA[LISA$state == state, seizure.start:seizure.end], 2, sum)[-48]
  diffusion_state_monthly <- data.frame(adj_inc=adj_inc_monthly,
                                        non_adj_inc=non_adj_inc_monthly,
                                        adj_dec=adj_dec_monthly,
                                        non_adj_dec=non_adj_dec_monthly,
                                        seizure.weight=seizure.weights.state.cocaine)
  if (sum(seizure.weights.state.cocaine) == 0) {
    result_i <- data.frame(state=state, coef=NA, p_value=NA, sig.="No Seizures", R_sq=NA)
    adj_inc_reg_cocaine <- rbind(adj_inc_reg_cocaine, result_i)
    nadj_inc_reg_cocaine <- rbind(nadj_inc_reg_cocaine, result_i)
    adj_dec_reg_cocaine <- rbind(adj_dec_reg_cocaine, result_i)
    nadj_dec_reg_cocaine <- rbind(nadj_dec_reg_cocaine, result_i)
    next 
  }
  
  if (sum(diffusion_state_monthly$adj_inc) == 0) result_i <- data.frame(state=state, coef=NA, p_value=NA, sig.="No adj_inc", R_sq=NA) else{
    adj_inc_reg <- lm(adj_inc~seizure.weight, data=diffusion_state_monthly) %>% summary
    c2 <- adj_inc_reg$coefficients[2,1] %>% round(digits=5)
    c3 <- adj_inc_reg$coefficients[2,4]
    c4 <- ifelse(c3 > 0.05, " ", "*")
    c5 <- adj_inc_reg$r.squared %>% round(digits=3)
    result_i <- data.frame(state=state, coef=c2, p_value=c3, sig.=c4, R_sq=c5)}
  
  adj_inc_reg_cocaine <- rbind(adj_inc_reg_cocaine, result_i)
  
  if (sum(diffusion_state_monthly$non_adj_inc) == 0) result_i <- data.frame(state=state, coef=NA, p_value=NA, sig.="No non_adj_inc", R_sq=NA) else{
    nadj_inc_reg <- lm(non_adj_inc~seizure.weight, data=diffusion_state_monthly) %>% summary
    c2 <- nadj_inc_reg$coefficients[2,1] %>% round(digits=5)
    c3 <- nadj_inc_reg$coefficients[2,4]
    c4 <- ifelse(c3 > 0.05, " ", "*")
    c5 <- nadj_inc_reg$r.squared %>% round(digits=3)
    result_i <- data.frame(state=state, coef=c2, p_value=c3, sig.=c4, R_sq=c5)}
  
  nadj_inc_reg_cocaine <- rbind(nadj_inc_reg_cocaine, result_i)
  
  if (sum(diffusion_state_monthly$adj_dec) == 0) result_i <- data.frame(state=state, coef=NA, p_value=NA, sig.="No adj_dec", R_sq=NA) else{
    adj_dec_reg <- lm(adj_dec~seizure.weight, data=diffusion_state_monthly) %>% summary
    c2 <- adj_dec_reg$coefficients[2,1] %>% round(digits=5)
    c3 <- adj_dec_reg$coefficients[2,4]
    c4 <- ifelse(c3 > 0.05, " ", "*")
    c5 <- adj_dec_reg$r.squared %>% round(digits=3)
    result_i <- data.frame(state=state, coef=c2, p_value=c3, sig.=c4, R_sq=c5)}
  
  adj_dec_reg_cocaine <- rbind(adj_dec_reg_cocaine, result_i)
  
  if (sum(diffusion_state_monthly$non_adj_dec) == 0) result_i <- data.frame(state=state, coef=NA, p_value=NA, sig.="No non_adj_dec", R_sq=NA) else{
    nadj_dec_reg <- lm(non_adj_dec~seizure.weight, data=diffusion_state_monthly) %>% summary
    c2 <- nadj_dec_reg$coefficients[2,1] %>% round(digits=5)
    c3 <- nadj_dec_reg$coefficients[2,4]
    c4 <- ifelse(c3 > 0.05, " ", "*")
    c5 <- nadj_dec_reg$r.squared %>% round(digits=3)
    result_i <- data.frame(state=state, coef=c2, p_value=c3, sig.=c4, R_sq=c5)}
  
  nadj_dec_reg_cocaine <- rbind(nadj_dec_reg_cocaine, result_i)
}
adj_inc_reg_cocaine <- adj_inc_reg_cocaine[-1,]
nadj_inc_reg_cocaine <- nadj_inc_reg_cocaine[-1,]
adj_dec_reg_cocaine <- adj_dec_reg_cocaine[-1,]
nadj_dec_reg_cocaine <- nadj_dec_reg_cocaine[-1,]

adj_inc_reg_cocaine
adj_inc_reg_cocaine %>% arrange(desc(R_sq))
nadj_inc_reg_cocaine
nadj_inc_reg_cocaine %>% arrange(desc(R_sq))
adj_dec_reg_cocaine
adj_dec_reg_cocaine %>% arrange(desc(R_sq))
nadj_dec_reg_cocaine
nadj_dec_reg_cocaine %>% arrange(desc(R_sq))

high_rsq_states <- (adj_inc_reg_cocaine %>% arrange(desc(R_sq)))$state[1:5]
diffusion_summary_count[diffusion_summary_count$state %in% high_rsq_states,]

# write.csv(adj_inc_reg_cocaine, "adj_inc_reg_cocaine (4x4).csv", row.names=F)
# write.csv(nadj_inc_reg_cocaine, "nadj_inc_reg_cocaine (4x4).csv", row.names=F)

LISA_all <- read.csv("AllKGLISAResults.csv") %>% as.tibble
seizure_all <- LISA_all[,20:67]
seizure_all_monthly <- apply(seizure_all, 2, sum)

{
par(mfrow=c(2,2))
plot(seizure_all_monthly[-48], diffusion_monthly$adj_inc,
     main="cocaine", ylab="adj_inc at t+1", xlab="Seizure Weights at t (kg)", pch=19)
plot(seizure_all_monthly[-48], diffusion_monthly$non_adj_inc,
     ylab="non_adj_inc at t+1", xlab="", pch=19)

plot(seizure_all_monthly[-48], diffusion_monthly$adj_dec,
     ylab="adj_dec at t+1", xlab="Seizure Weights at t (kg)", pch=19)
plot(seizure_all_monthly[-48], diffusion_monthly$non_adj_dec,
     ylab="non_adj_dec at t+1", xlab="", pch=19)
par(mfrow=c(1,1))
}

{
par(mfrow=c(2,2))
plot(log(seizure_all_monthly[-48]), diffusion_monthly$adj_inc,
     main="cocaine", ylab="adj_inc at t+1", xlab="Seizure Weights at t (log kg)", pch=19)
plot(log(seizure_all_monthly[-48]), diffusion_monthly$non_adj_inc,
     ylab="non_adj_inc at t+1", xlab="", pch=19)

plot(log(seizure_all_monthly[-48]), diffusion_monthly$adj_dec,
     ylab="adj_dec at t+1", xlab="Seizure Weights at t (log kg)", pch=19)
plot(log(seizure_all_monthly[-48]), diffusion_monthly$non_adj_dec,
     ylab="non_adj_dec at t+1", xlab="", pch=19)
par(mfrow=c(1,1))
}
