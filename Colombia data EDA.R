# setwd("/Users/euseongjang/Documents/R")
# setwd("C:/Users/gkfrj/Documents/R")
library(fpp2)
library(readxl)
library(urbnmapr)
library(tidyverse)
library(gridExtra)
library(lubridate)

cultivation <- read_xlsx("Colombia Data/Colombia Coca Cultivation 1999-2016 (Ha).xlsx")
eradication_aerial <- read_xlsx("Colombia Data/Colombia-Coca Erradication-1994-2021 (Aerial).xlsx")
eradication_manual <- read_xlsx("Colombia Data/Colombia-Coca Erradication-1994-2021 (Manual).xlsx")

cultivation$n_obs <- cultivation %>% apply(1, function(x) sum(!is.na(x[5:22])))
eradication_aerial$n_obs <- eradication_aerial %>% apply(1, function(x) sum(!is.na(x[5:26])))
eradication_manual$n_obs <- eradication_manual %>% apply(1, function(x) sum(!is.na(x[5:29])))

cultivation %>%
  arrange(desc(n_obs), DEPARTAMENTO, MUNICIPIO) %>% 
  select(DEPARTAMENTO, MUNICIPIO, n_obs)
eradication_aerial %>%
  arrange(Departamento, MUNICIPIO) %>% 
  select(Departamento, MUNICIPIO, n_obs) %>% 
  as.data.frame
eradication_manual %>%
  arrange(desc(n_obs), Departamento, MUNICIPIO) %>% 
  select(Departamento, MUNICIPIO, n_obs)

cultivation %>% filter(DEPARTAMENTO == "ANTIOQUIA" & MUNICIPIO %in% (eradication_aerial %>% filter(Departamento == "ANTIOQUIA") %>% pull(MUNICIPIO) %>% unique))

cultivation %>% filter(DEPARTAMENTO == "ANTIOQUIA" & MUNICIPIO == "ANORÍ") -> cultivation_ANORÍ
eradication_aerial %>% filter(Departamento == "ANTIOQUIA" & MUNICIPIO == "ANORÍ") -> eradication_aerial_ANORÍ
cultivation_ANORÍ_ts <- ts(t(cultivation_ANORÍ[9:21]), start=2003)
eradication_aerial_ANORÍ_ts <- ts(t(eradication_aerial_ANORÍ[14:26]), start=2003)
cor(cultivation_ANORÍ_ts, eradication_aerial_ANORÍ_ts)
