# setwd("/Users/R")
# setwd("C:/Users/gkfrj/Documents/R")
library(fpp2)
library(readxl)
library(urbnmapr)
library(tidyverse)
library(gridExtra)
library(lubridate)
library(colmaps)

## check misspelling in regional names and correct
cultivation <- read_xlsx("Colombia Data/Colombia Coca Cultivation 1999-2016 (Ha).xlsx")
labs_HCl <- read_xlsx("Colombia Data/Colombia-Laboratories-1997-2022 (COCAINE HYDROCHLORIDE).xlsx")
labs_PPI <- read_xlsx("Colombia Data/Colombia-Laboratories-1997-2022 (PRIMARY PRODUCTION INFRASTRUCTURE).xlsx")
armed_groups <- list()
armed_groups$y2008 <- read_xlsx("Colombia Data/Colombia-Armed groups-Paramilitar 2008.xlsx") %>% filter(Pais == "Colombia") %>% select(-Pais)
armed_groups$y2010 <- read_xlsx("Colombia Data/Colombia-Armed groups-Paramilitar 2010.xlsx") %>% filter(Pais == "Colombia") %>% select(-Pais)
armed_groups$y2011 <- read_xlsx("Colombia Data/Colombia-Armed groups-Paramilitar 2011.xlsx") %>% filter(Pais == "Colombia") %>% select(-Pais)
armed_groups$y2012 <- read_xlsx("Colombia Data/Colombia-Armed groups-Paramilitar 2012.xlsx") %>% filter(Pais == "Colombia") %>% select(-Pais)
armed_groups$y2013 <- read_xlsx("Colombia Data/Colombia-Armed groups-Paramilitar 2013.xlsx") %>% filter(Pais == "Colombia") %>% select(-Pais)
armed_groups$y2014 <- read_xlsx("Colombia Data/Colombia-Armed groups-Paramilitar 2014.xlsx") %>% filter(Pais == "Colombia") %>% select(-Pais)
armed_groups$y2016 <- read_xlsx("Colombia Data/Colombia-Armed groups-Paramilitar 2016.xlsx") %>% filter(Pais == "Colombia") %>% select(-Pais)

armed_groups$y2008 %>% apply(1, function(x) sum(!is.na(x[-(1:4)])))

cultivation$DEPARTAMENTO <- str_to_title(cultivation$DEPARTAMENTO, locale="sp")
cultivation$DEPARTAMENTO <- gsub("Del ", "del ", cultivation$DEPARTAMENTO)
cultivation$DEPARTAMENTO <- gsub("De ", "de ", cultivation$DEPARTAMENTO)
cultivation$MUNICIPIO <- str_to_title(cultivation$MUNICIPIO, locale="sp")
cultivation$MUNICIPIO <- gsub("Del ", "del ", cultivation$MUNICIPIO)
cultivation$MUNICIPIO <- gsub("De ", "de ", cultivation$MUNICIPIO)
labs_HCl$DEPARTAMENTO <- str_to_title(labs_HCl$DEPARTAMENTO, locale="sp")
labs_HCl$DEPARTAMENTO <- gsub("Del ", "del ", labs_HCl$DEPARTAMENTO)
labs_HCl$DEPARTAMENTO[which(labs_HCl$DEPARTAMENTO == "Bogotá D.c.")] <- "Bogotá D.C."
labs_HCl$DEPARTAMENTO <- gsub("De ", "de ", labs_HCl$DEPARTAMENTO)
labs_HCl$MUNICIPIO <- str_to_title(labs_HCl$MUNICIPIO, locale="sp")
labs_HCl$MUNICIPIO <- gsub("Del ", "del ", labs_HCl$MUNICIPIO)
labs_HCl$MUNICIPIO <- gsub("De ", "de ", labs_HCl$MUNICIPIO)
labs_HCl$MUNICIPIO[which(labs_HCl$MUNICIPIO == "Bogotá D.c.")] <- "Bogotá D.C."
labs_PPI$DEPARTAMENTO <- str_to_title(labs_PPI$DEPARTAMENTO, locale="sp")
labs_PPI$DEPARTAMENTO <- gsub("Del ", "del ", labs_PPI$DEPARTAMENTO)
labs_PPI$DEPARTAMENTO <- gsub("De ", "de ", labs_PPI$DEPARTAMENTO)
labs_PPI$DEPARTAMENTO[which(labs_PPI$DEPARTAMENTO == "Bogotá D.c.")] <- "Bogotá D.C."
labs_PPI$MUNICIPIO <- str_to_title(labs_PPI$MUNICIPIO, locale="sp")
labs_PPI$MUNICIPIO <- gsub("Del ", "del ", labs_PPI$MUNICIPIO)
labs_PPI$MUNICIPIO <- gsub("De ", "de ", labs_PPI$MUNICIPIO)
labs_PPI$MUNICIPIO[which(labs_PPI$MUNICIPIO == "Bogotá D.c.")] <- "Bogotá D.C."

cultivation$MUNICIPIO[grep("\\(", cultivation$MUNICIPIO)] <- substr(cultivation$MUNICIPIO[grep("\\(", cultivation$MUNICIPIO)],
                                start=1,
                                stop=regexpr("\\(", cultivation$MUNICIPIO[grep("\\(", cultivation$MUNICIPIO)])-2)

labs_HCl$MUNICIPIO[grep("\\(", labs_HCl$MUNICIPIO)] <- substr(labs_HCl$MUNICIPIO[grep("\\(", labs_HCl$MUNICIPIO)],
                                                              start=1,
                                                              stop=regexpr("\\(", labs_HCl$MUNICIPIO[grep("\\(", labs_HCl$MUNICIPIO)])-2)
labs_PPI$MUNICIPIO[grep("\\(", labs_PPI$MUNICIPIO)] <- substr(labs_PPI$MUNICIPIO[grep("\\(", labs_PPI$MUNICIPIO)],
                                                              start=1,
                                                              stop=regexpr("\\(", labs_PPI$MUNICIPIO[grep("\\(", labs_PPI$MUNICIPIO)])-2)

cultivation_longer <- cultivation %>% pivot_longer(-(CODDEPTO:MUNICIPIO), names_to="year", values_to="cultivation")
cultivation_longer$cultivation %>% summary

armed_groups_regions <- armed_groups$y2008 %>% select(DEPARTAMENTO, MUNICIPIO) %>% unique %>% arrange(DEPARTAMENTO, MUNICIPIO)
for (i in 2:7) {
  armed_groups_regions <- rbind(armed_groups_regions, armed_groups[[i]] %>% select(DEPARTAMENTO, MUNICIPIO) %>% unique %>% arrange(DEPARTAMENTO, MUNICIPIO)) %>% 
    unique %>% arrange(DEPARTAMENTO, MUNICIPIO)
}
# armed_groups_regions %>% write.csv("Colombia Data/armed groups Dep-Muni names.csv", row.names=F, fileEncoding="UTF-8")

cultivation %>% 
  filter(!(MUNICIPIO %in% armed_groups_regions$MUNICIPIO)) %>% 
  select(DEPARTAMENTO, MUNICIPIO) %>%
  arrange(DEPARTAMENTO, MUNICIPIO) -> cultivation_regions
labs_HCl %>% 
  filter(!(MUNICIPIO %in% armed_groups_regions$MUNICIPIO)) %>% 
  select(DEPARTAMENTO, MUNICIPIO) %>%
  arrange(DEPARTAMENTO, MUNICIPIO) -> labs_HCl_regions
labs_PPI %>% 
  filter(!(MUNICIPIO %in% armed_groups_regions$MUNICIPIO)) %>% 
  select(DEPARTAMENTO, MUNICIPIO) %>%
  arrange(DEPARTAMENTO, MUNICIPIO) -> labs_PPI_regions
# finished by Norosi

cultivation %>% filter(!(DEPARTAMENTO %in% armed_groups_regions$DEPARTAMENTO)) %>% pull(DEPARTAMENTO)
labs_HCl %>% filter(!(DEPARTAMENTO %in% armed_groups_regions$DEPARTAMENTO)) %>% pull(DEPARTAMENTO)
labs_PPI %>% filter(!(DEPARTAMENTO %in% armed_groups_regions$DEPARTAMENTO)) %>% pull(DEPARTAMENTO)

years <- c(2008, 2010:2014, 2016)
