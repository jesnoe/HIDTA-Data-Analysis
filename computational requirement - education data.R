library(spdep)
library(tidyverse)
library(gridExtra)
library(lubridate)
library(educationdata)

# In the 2022 Census of Governments, the United States Census Bureau reported 12,546 school district governments in the United States
# locations are point coordinates
edu_data <- get_education_data(level = "school-districts",
                               source = "ccd",
                               topic = "directory",
                               filters = list(year = 2022)) %>% as_tibble %>% 
  filter(!(state_mailing %in% c("AK", "HI")))
edu_data
edu_data$leaid %>% unique %>% length
names(edu_data)
edu_data_dist <- dist(edu_data %>% select(longitude, latitude))
edu_data_dist_tri <- lower.tri(edu_data_dist)
edu_data_dist_vec <- as.vector(edu_data_dist[edu_data_dist_tri])
summary(edu_data_dist_vec)

nb_edu <- knn2nb(knearneigh(edu_data %>% select(longitude, latitude), k=5), row.names=edu_data$leaid)
nb_obj_edu <- nb2listw(nb_edu, style="B")
alpha <- 0.05
nperm <- 9999

# teachers_elementary_fte: Number of full-time equivalent elementary school teachers
input_data <- edu_data$teachers_elementary_fte
input_data[is.na(input_data)] <- 0
set.seed(100)
t1 <- Sys.time()
localM.month <- localmoran_abs(input_data, nb_obj_edu, nsim=nperm, zero.policy=T, xx=NULL, alternative="two.sided")
t2 <- Sys.time()

set.seed(100)
t3 <- Sys.time()
localM.month <- localmoran_abs(input_data, nb_obj_edu, nsim=nperm, zero.policy=T, xx=NULL, alternative="two.sided", perm.i=T)
t4 <- Sys.time()

t2 - t1
t4 - t3
