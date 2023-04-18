# devtools::install_github("synth-inference/synthdid")
library(synthdid)
library(tidyverse)
library(readxl)
regional.covariates <- read_xlsx("regional_covariates.xlsx")
regional.covariates %>% filter(name=="Jena, Stadt") # 16053
regional.covariates %>% 
  filter((name %in% c("Baden-Baden, Stadt", "Berlin, Stadt", "Brandenburg an der Havel, Stadt", "Bremen, Stadt", "Hamburg, Stadt",
                      "Main-Kinzig-Kreis", "Rottweil", "Wolfsburg, Stadt", "Braunschweig, Stadt", "Nordhausen"))) %>% select(idlandkreis) -> mandatory.states
regional.covariates %>% 
  filter(!(name %in% c("Baden-Baden, Stadt", "Berlin, Stadt", "Brandenburg an der Havel, Stadt", "Bremen, Stadt", "Hamburg, Stadt",
                       "Main-Kinzig-Kreis", "Rottweil", "Wolfsburg, Stadt", "Braunschweig, Stadt", "Nordhausen", "Jena, Stadt"))) %>%
  select(idlandkreis) -> non.mandatory.states
covid <- read_xlsx("01_Covid19_cases_by_regions_overall.xlsx")
covid <- covid %>% filter(!(idlandkreis %in% mandatory.states$idlandkreis), datum > "2020-03-05")
covid$day <- covid$day - 38
covid$treatment <- ifelse(covid$idlandkreis == 16053 & covid$datum > "2020-04-05", 1, 0)
covid$idlandkreis <- as.factor(covid$idlandkreis)
covid %>% filter(idlandkreis == 16053) %>% as.data.frame
covid <- as.data.frame(covid)

covid2 <- covid %>% filter(id %in% c(sample(non.mandatory.states$idlandkreis, 100), 16053)) %>% as.data.frame
input <- panel.matrices(covid, unit=1, time=5, outcome=6, treatment=9)
tau.sdid <- synthdid_estimate(input$Y, input$N0, input$T0)
summary(tau.sdid)
plot(tau.sdid)
plot(tau.sdid, overlay=-5)

tau.sc <- sc_estimate(input$Y, input$N0, input$T0)
tau.did <- did_estimate(input$Y, input$N0, input$T0)
estimates <- list(tau.did, tau.sc, tau.sdid)
unlist(estimates)
names(estimates) <- c('Diff-in-Diff', 'Synthetic Control', 'Synthetic Diff-in-Diff')
synthdid_plot(estimates) + ylim(c(-100, 400))
synthdid_plot(estimates, overlay=-5)
synthdid_units_plot(estimates)

tau.did <- c()
tau.sc <- c()
tau.sdid <- c()
set.seed(1000)
for (i in 1:100) {
  covid2 <- covid %>% filter(id %in% c(sample(non.mandatory.states$idlandkreis, 100), 16053)) %>% as.data.frame
  input <- panel.matrices(covid2, unit=1, time=5, outcome=6, treatment=9)
  tau.sdid <- c(tau.sdid, synthdid_estimate(input$Y, input$N0, input$T0))
  tau.sc <- c(tau.sc, sc_estimate(input$Y, input$N0, input$T0))
  tau.did <- c(tau.did, did_estimate(input$Y, input$N0, input$T0))
}
means <- list(mean(tau.did), mean(tau.sc),  mean(tau.sdid))
names(means) <- c('Diff-in-Diff', 'Synthetic Control', 'Synthetic Diff-in-Diff')
sds <- list(sd(tau.did), sd(tau.sc),  sd(tau.sdid))
names(sds) <- c('Diff-in-Diff', 'Synthetic Control', 'Synthetic Diff-in-Diff')
unlist(means)
unlist(sds)
