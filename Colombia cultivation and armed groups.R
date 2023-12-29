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
{
# cultivation <- read_xlsx("Colombia Data/Colombia Coca Cultivation 1999-2016 (Ha).xlsx")
# labs_HCl <- read_xlsx("Colombia Data/Colombia-Laboratories-1997-2022 (COCAINE HYDROCHLORIDE).xlsx")
# labs_PPI <- read_xlsx("Colombia Data/Colombia-Laboratories-1997-2022 (PRIMARY PRODUCTION INFRASTRUCTURE).xlsx")
armed_groups <- list()
armed_groups$y2008 <- read_xlsx("Colombia Data/Colombia-Armed groups-Paramilitar 2008.xlsx") %>% filter(Pais == "Colombia") %>% select(-Pais)
armed_groups$y2010 <- read_xlsx("Colombia Data/Colombia-Armed groups-Paramilitar 2010.xlsx") %>% filter(Pais == "Colombia") %>% select(-Pais)
armed_groups$y2011 <- read_xlsx("Colombia Data/Colombia-Armed groups-Paramilitar 2011.xlsx") %>% filter(Pais == "Colombia") %>% select(-Pais)
armed_groups$y2012 <- read_xlsx("Colombia Data/Colombia-Armed groups-Paramilitar 2012.xlsx") %>% filter(Pais == "Colombia") %>% select(-Pais)
armed_groups$y2013 <- read_xlsx("Colombia Data/Colombia-Armed groups-Paramilitar 2013.xlsx") %>% filter(Pais == "Colombia") %>% select(-Pais)
armed_groups$y2014 <- read_xlsx("Colombia Data/Colombia-Armed groups-Paramilitar 2014.xlsx") %>% filter(Pais == "Colombia") %>% select(-Pais)
armed_groups$y2016 <- read_xlsx("Colombia Data/Colombia-Armed groups-Paramilitar 2016.xlsx") %>% filter(Pais == "Colombia") %>% select(-Pais)

# 
# cultivation$DEPARTAMENTO <- str_to_title(cultivation$DEPARTAMENTO, locale="sp")
# cultivation$DEPARTAMENTO <- gsub("Del ", "del ", cultivation$DEPARTAMENTO)
# cultivation$DEPARTAMENTO <- gsub("De ", "de ", cultivation$DEPARTAMENTO)
# cultivation$MUNICIPIO <- str_to_title(cultivation$MUNICIPIO, locale="sp")
# cultivation$MUNICIPIO <- gsub("Del ", "del ", cultivation$MUNICIPIO)
# cultivation$MUNICIPIO <- gsub("De ", "de ", cultivation$MUNICIPIO)
# labs_HCl$DEPARTAMENTO <- str_to_title(labs_HCl$DEPARTAMENTO, locale="sp")
# labs_HCl$DEPARTAMENTO <- gsub("Del ", "del ", labs_HCl$DEPARTAMENTO)
# labs_HCl$DEPARTAMENTO[which(labs_HCl$DEPARTAMENTO == "Bogotá D.c.")] <- "Bogotá D.C."
# labs_HCl$DEPARTAMENTO <- gsub("De ", "de ", labs_HCl$DEPARTAMENTO)
# labs_HCl$MUNICIPIO <- str_to_title(labs_HCl$MUNICIPIO, locale="sp")
# labs_HCl$MUNICIPIO <- gsub("Del ", "del ", labs_HCl$MUNICIPIO)
# labs_HCl$MUNICIPIO <- gsub("De ", "de ", labs_HCl$MUNICIPIO)
# labs_HCl$MUNICIPIO[which(labs_HCl$MUNICIPIO == "Bogotá D.c.")] <- "Bogotá D.C."
# labs_PPI$DEPARTAMENTO <- str_to_title(labs_PPI$DEPARTAMENTO, locale="sp")
# labs_PPI$DEPARTAMENTO <- gsub("Del ", "del ", labs_PPI$DEPARTAMENTO)
# labs_PPI$DEPARTAMENTO <- gsub("De ", "de ", labs_PPI$DEPARTAMENTO)
# labs_PPI$DEPARTAMENTO[which(labs_PPI$DEPARTAMENTO == "Bogotá D.c.")] <- "Bogotá D.C."
# labs_PPI$MUNICIPIO <- str_to_title(labs_PPI$MUNICIPIO, locale="sp")
# labs_PPI$MUNICIPIO <- gsub("Del ", "del ", labs_PPI$MUNICIPIO)
# labs_PPI$MUNICIPIO <- gsub("De ", "de ", labs_PPI$MUNICIPIO)
# labs_PPI$MUNICIPIO[which(labs_PPI$MUNICIPIO == "Bogotá D.c.")] <- "Bogotá D.C."
# 
# cultivation$MUNICIPIO[grep("\\(", cultivation$MUNICIPIO)] <- substr(cultivation$MUNICIPIO[grep("\\(", cultivation$MUNICIPIO)],
#                                 start=1,
#                                 stop=regexpr("\\(", cultivation$MUNICIPIO[grep("\\(", cultivation$MUNICIPIO)])-2)
# 
# labs_HCl$MUNICIPIO[grep("\\(", labs_HCl$MUNICIPIO)] <- substr(labs_HCl$MUNICIPIO[grep("\\(", labs_HCl$MUNICIPIO)],
#                                                               start=1,
#                                                               stop=regexpr("\\(", labs_HCl$MUNICIPIO[grep("\\(", labs_HCl$MUNICIPIO)])-2)
# labs_PPI$MUNICIPIO[grep("\\(", labs_PPI$MUNICIPIO)] <- substr(labs_PPI$MUNICIPIO[grep("\\(", labs_PPI$MUNICIPIO)],
#                                                               start=1,
#                                                               stop=regexpr("\\(", labs_PPI$MUNICIPIO[grep("\\(", labs_PPI$MUNICIPIO)])-2)
}
# cultivation_longer <- cultivation %>% pivot_longer(-(CODDEPTO:MUNICIPIO), names_to="year", values_to="cultivation")
# cultivation_longer$cultivation %>% summary
# 
# armed_groups_regions <- armed_groups$y2008 %>% select(DEPARTAMENTO, MUNICIPIO) %>% unique %>% arrange(DEPARTAMENTO, MUNICIPIO)
# for (i in 2:7) {
#   armed_groups_regions <- rbind(armed_groups_regions, armed_groups[[i]] %>% select(DEPARTAMENTO, MUNICIPIO) %>% unique %>% arrange(DEPARTAMENTO, MUNICIPIO)) %>% 
#     unique %>% arrange(DEPARTAMENTO, MUNICIPIO)
# }
# cultivation %>% 
#   filter(!(MUNICIPIO %in% armed_groups_regions$MUNICIPIO)) %>% 
#   select(DEPARTAMENTO, MUNICIPIO) %>%
#   arrange(DEPARTAMENTO, MUNICIPIO)
# labs_HCl %>% 
#   filter(!(MUNICIPIO %in% armed_groups_regions$MUNICIPIO)) %>% 
#   select(DEPARTAMENTO, MUNICIPIO) %>%
#   arrange(DEPARTAMENTO, MUNICIPIO)
# labs_PPI %>% 
#   filter(!(MUNICIPIO %in% armed_groups_regions$MUNICIPIO)) %>% 
#   select(DEPARTAMENTO, MUNICIPIO) %>%
#   arrange(DEPARTAMENTO, MUNICIPIO)


# armed_groups_regions %>% write.csv("Colombia Data/armed groups Dep-Muni names.csv", row.names=F, fileEncoding="UTF-8")
# cultivation %>% write.csv("Colombia Data/Colombia Coca Cultivation 1999-2016 renamed (Ha).csv", row.names=F, fileEncoding="UTF-8")
# labs_HCl %>% write.csv("Colombia Data/Colombia-Laboratories-1997-2022 renamed (COCAINE HYDROCHLORIDE).csv", row.names=F, fileEncoding="UTF-8")
# labs_PPI %>% write.csv("Colombia Data/Colombia-Laboratories-1997-2022 renamed (PRIMARY PRODUCTION INFRASTRUCTURE).csv", row.names=F, fileEncoding="UTF-8")
cultivation <- read.csv("Colombia Data/Colombia Coca Cultivation 1999-2016 renamed (Ha).csv") %>% as_tibble
labs_HCl <- read.csv("Colombia Data/Colombia-Laboratories-1997-2022 renamed (COCAINE HYDROCHLORIDE).csv") %>% as_tibble
labs_PPI <- read.csv("Colombia Data/Colombia-Laboratories-1997-2022 renamed (PRIMARY PRODUCTION INFRASTRUCTURE).csv") %>% as_tibble

armed_groups$y2008 %>% apply(1, function(x) sum(!is.na(x[-(1:4)])))
armed_groups$y2008[,-(1:3)] %>% as.vector %>% unlist %>% table

armed_groups_table <- list()
armed_groups_combined <- list()
years <- c(2008, 2010:2014, 2016)
for (year in years) {
  year_name <- paste0("y", year)
  
  armed_groups_vec <- armed_groups[[year_name]][,-(1:3)] %>% apply(2, function(x) sum(x, na.rm=T))
  armed_groups_tb_year <- as_tibble(as.list(c(armed_groups_vec)))
  armed_groups_tb_year <- armed_groups_tb_year %>% pivot_longer(cols=everything(),
                                                      names_to="armed_group",
                                                      values_to="total_num_of_Muni.")
  armed_groups_table[[year_name]] <- armed_groups_tb_year
  
  armed_groups_year <- armed_groups[[year_name]]
  armed_groups_combined[[year_name]] <- full_join(armed_groups_year, cultivation[,c(2, 4, year-1999+5)], by=c("DEPARTAMENTO", "MUNICIPIO"))
  armed_groups_combined[[year_name]] <- full_join(armed_groups_combined[[year_name]], labs_HCl[,c(2, 4, year-1997+5)], by=c("DEPARTAMENTO", "MUNICIPIO"))
  armed_groups_combined[[year_name]] <- full_join(armed_groups_combined[[year_name]], labs_PPI[,c(2, 4, year-1997+5)], by=c("DEPARTAMENTO", "MUNICIPIO"))
  num_col <- ncol(armed_groups_combined[[year_name]])
  names(armed_groups_combined[[year_name]])[(num_col-2):num_col] <- c(paste0("cultivation_", year), paste0("labs_HCl_", year), paste0("labs_PPI_", year))
  armed_groups_combined[[year_name]] <- armed_groups_combined[[year_name]][,c(1:3, (num_col-2):num_col, 4:(num_col-3))]
}

armed_groups_table$y2008
armed_groups_combined$y2008

count_concurring_cultivation <- function(x) {
  n_columns <- length(x)
  x <- ifelse(is.na(x), 0, 1)
  result <- x[[4]] * x[8:n_columns]
  return(result)
}

count_concurring_labs_HCl <- function(x) {
  n_columns <- length(x)
  x <- ifelse(is.na(x), 0, 1)
  result <- x[[5]] * x[8:n_columns]
  return(result)
}

count_concurring_labs_PPI <- function(x) {
  n_columns <- length(x)
  x <- ifelse(is.na(x), 0, 1)
  result <- x[[6]] * x[8:n_columns]
  return(result)
}

concurrence_table <- tibble(year=0,
                            n_municipios=0,
                            n_cultivation=0,
                            n_labs_HCl=0,
                            n_labs_PPI=0,
                            cultivation_labs_HCl=0,
                            cultivation_labs_PPI=0,
                            labs_HCl_PPI=0)
for (year in years) {
  year_name <- paste0("y", year)
  armed_groups_combined_year <- armed_groups_combined[[year_name]]
  n_municipios <- nrow(armed_groups_combined_year)
  cultivation_year <- armed_groups_combined_year[,4] %>% as.vector %>% unlist
  labs_HCl_year <- armed_groups_combined_year[,5] %>% as.vector %>% unlist
  labs_PPI_year <- armed_groups_combined_year[,6] %>% as.vector %>% unlist
  concurrence_table_i <- c(year,
                           n_municipios,
                           sum(as.logical(cultivation_year), na.rm=T),
                           sum(as.logical(labs_HCl_year), na.rm=T),
                           sum(as.logical(labs_PPI_year), na.rm=T),
                           sum(as.logical(cultivation_year) & as.logical(labs_HCl_year), na.rm=T),
                           sum(as.logical(cultivation_year) & as.logical(labs_PPI_year), na.rm=T),
                           sum(as.logical(labs_HCl_year) & as.logical(labs_PPI_year), na.rm=T)
                           )
  concurrence_table <- rbind(concurrence_table, concurrence_table_i)
  
  n_columns <- ncol(armed_groups_combined_year)
  
  cultivation_concurrence <- armed_groups_combined_year %>% 
    apply(1, count_concurring_cultivation) %>% 
    apply(1, sum) %>% as.list %>% as_tibble %>% 
    pivot_longer(cols=everything(),
                 names_to="armed_group",
                 values_to="cultivation_concurrence")
  labs_HCl_concurrence <- armed_groups_combined_year %>% 
    apply(1, count_concurring_labs_HCl) %>% 
    apply(1, sum) %>% as.list %>% as_tibble%>% 
    pivot_longer(cols=everything(),
                 names_to="armed_group",
                 values_to="labs_HCl_concurrence")
  labs_PPI_concurrence <- armed_groups_combined_year %>% 
    apply(1, count_concurring_labs_PPI) %>% 
    apply(1, sum) %>% as.list %>% as_tibble%>% 
    pivot_longer(cols=everything(),
                 names_to="armed_group",
                 values_to="labs_PPI_concurrence")
  concurrence_year <- inner_join(cultivation_concurrence, labs_HCl_concurrence, by="armed_group") %>% 
    inner_join(labs_PPI_concurrence, by="armed_group")
  armed_groups_table[[year_name]] <- inner_join(armed_groups_table[[year_name]], concurrence_year, by="armed_group")
}

(concurrence_table <- concurrence_table[-1,])
armed_groups_table$y2008 %>% as.data.frame
armed_groups_table$y2016 %>% as.data.frame


armed_groups_combined$y2016 %>% 
  select(REGION:labs_PPI_2016, `Clan del Golfo (Formerly Los Urabeños)`) %>% 
  filter(!is.na(`Clan del Golfo (Formerly Los Urabeños)`) &
           !is.na(cultivation_2016))# %>% write.csv("Colombia Data/Clan del Golfo region with cultivation 2016.csv", row.names=F)

armed_groups_combined$y2016 %>% 
  select(REGION:labs_PPI_2016, `Clan del Golfo (Formerly Los Urabeños)`) %>% 
  filter(!is.na(`Clan del Golfo (Formerly Los Urabeños)`) &
           is.na(cultivation_2016) &
           is.na(labs_HCl_2016) &
           is.na(labs_PPI_2016))# %>% write.csv("Colombia Data/Clan del Golfo region without activities 2016.csv", row.names=F)
