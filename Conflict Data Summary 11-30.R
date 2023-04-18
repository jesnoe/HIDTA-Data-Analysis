library(readxl)
library(tidyverse)
library(data.table)
library(ggplot2)
library(colmaps) # For loading departamentos

conflict <- read_xlsx("UCDP.xlsx")
setDT(conflict)
colnames(conflict)
summary(conflict %>% select(year, active_year, type_of_violence, deaths_a, deaths_b, deaths_civilians, deaths_unknown, best, high, low, gwnoa, gwnob))


department <- gsub(" department", "", unique(conflict$adm_1))
departamentos@data$depto[11] <- "Norte de Santander"
departamentos@data$depto[17] <- "Valle del Cauca"
departamentos@data$depto[29] <- "Bogotá"
sum(department %in% departamentos@data$depto)
department[which(!(department %in% departamentos@data$depto))] # "Archipiélago De San Andrés, Providencia Y Santa Catalina" or missing value?
table(conflict$adm_1)
conflict$adm_1 <- as.factor(conflict$adm_1)
summary(conflict %>% select(adm_1))

sum(is.na(conflict$adm_1))
na.location <- data.frame(x=conflict$longitude[is.na(conflict$adm_1)], y=conflict$latitude[is.na(conflict$adm_1)]) %>% group_by(x, y) %>% summarise(n.points=n())
# Check where these are on the Colombian map


departamentos@data$id <- as.factor(departamentos@data$id)
positions <- data.frame()
centers <- data.frame()
for (i in 1:33) {
  positions <- rbind(positions, data.frame(departamentos@polygons[[i]]@Polygons[[1]]@coords,
                                           id=rep(departamentos@polygons[[i]]@ID, nrow(departamentos@polygons[[i]]@Polygons[[1]]@coords))))
  centers <- rbind(centers, data.frame(x=(max(positions$x[positions$id==departamentos@polygons[[i]]@ID]) + min(positions$x[positions$id==departamentos@polygons[[i]]@ID]))/2,
                                       y=(max(positions$y[positions$id==departamentos@polygons[[i]]@ID]) + min(positions$y[positions$id==departamentos@polygons[[i]]@ID]))/2))
}
departamentos@data <- cbind(departamentos@data, x=centers$x, y=centers$y)
positions$id <- as.factor(positions$id)
dep.id <- factor(departamentos@data$id)
dep.cols <- data.frame(id=dep.id, col=(1:length(dep.id)))
ggplot(dep.cols, fill=1) + geom_polygon(data=positions, aes(x=x, y=y, group=id), fill=NA, colour="black") + ggtitle("Locations of Points With Department Missing") +
  geom_point(data=na.location, aes(x=x,y=y), colour="red") + geom_text(data=na.location, aes(x=x,, y=y, label=n.points), nudge_x=-0.2, nudge_y=0.2)+
  geom_text(aes(x=x, y=y, label=adm_1), data=departamentos@data, size=3, angle=20, fontface = "bold")
#  geom_map(aes(map_id=dep.id), map=positions) + expand_limits(positions) +

which(colnames(conflict) == "deaths_a")
conflict <- conflict[-is.na(conflict$adm_1),]
conflict$adm_1 <- gsub(" department", "", conflict$adm_1)
conflict$n.victims <- apply(conflict[, 41:47], 1, sum)
colnames(departamentos@data)[2] <- "adm_1"
conflict <- merge(conflict, departamentos@data[,1:2], by="adm_1")
colnames(conflict)[length(conflict)] <- "dep.id"
conflict.summary <- conflict %>% group_by(year, dep.id, adm_1) %>% summarise(n.conflicts=n(), n.victims=sum(n.victims))
plot(conflict.summary$n.conflicts, conflict.summary$n.victims, main="Correlation Between the Number of Conflicts and Victims", 
     xlab="Number of Conflicts", ylab="Number of Victims", pch=19, cex=0.5)

# Plot n.conflicts
n.conflicts.i <- (conflict.summary %>% filter(year==2000) %>% select(dep.id, n.conflicts))[,-1]
colnames(n.conflicts.i)[1] <- "id"
missing.index <- which(!(departamentos@data$id %in% n.conflicts.i$id))
n.conflicts.i <- rbind(n.conflicts.i, data.frame(id=departamentos@data$id[missing.index], n.conflicts=rep(0, length(missing.index))))
ggplot(data=n.conflicts.i, aes(fill=n.conflicts)) +  geom_map(aes(map_id=id), map=positions) + expand_limits(positions) +
  ggtitle("Number of Conflicts In 2000")

max(conflict.summary$n.conflicts) # 355
min(conflict.summary$n.conflicts) # 1
colnames(conflict.summary)[2] <- "id"
ggplot(data=conflict.summary %>% filter(year <= 1994), aes(fill=n.conflicts)) +geom_map(aes(map_id=id), map=positions) + expand_limits(positions) +
  ggtitle("Number of Conflicts") + facet_wrap(vars(year)) + geom_polygon(data=positions, aes(x=x, y=y, group=id), fill=NA, colour="black") +
  scale_fill_gradient(low="green", high="red", limits=c(0, 355), oob=scales::squish)

ggplot(data=conflict.summary %>% filter(year > 1994, year <= 2003), aes(fill=n.conflicts)) +geom_map(aes(map_id=id), map=positions) + expand_limits(positions) +
  ggtitle("Number of Conflicts") + facet_wrap(vars(year)) + geom_polygon(data=positions, aes(x=x, y=y, group=id), fill=NA, colour="black") +
  scale_fill_gradient(low="green", high="red", limits=c(0, 355), oob=scales::squish)

ggplot(data=conflict.summary %>% filter(year > 2003, year <= 2012), aes(fill=n.conflicts)) +geom_map(aes(map_id=id), map=positions) + expand_limits(positions) +
  ggtitle("Number of Conflicts") + facet_wrap(vars(year)) + geom_polygon(data=positions, aes(x=x, y=y, group=id), fill=NA, colour="black") +
  scale_fill_gradient(low="green", high="red", limits=c(0, 355), oob=scales::squish)

ggplot(data=conflict.summary %>% filter(year > 2012), aes(fill=n.conflicts)) +geom_map(aes(map_id=id), map=positions) + expand_limits(positions) +
  ggtitle("Number of Conflicts") + facet_wrap(vars(year)) + geom_polygon(data=positions, aes(x=x, y=y, group=id), fill=NA, colour="black") +
  scale_fill_gradient(low="green", high="red", limits=c(0, 355), oob=scales::squish)

  # Lowering maximum
ggplot(data=conflict.summary %>% filter(year > 2003, year <= 2012), aes(fill=n.conflicts)) +geom_map(aes(map_id=id), map=positions) + expand_limits(positions) +
  ggtitle("Number of Conflicts") + facet_wrap(vars(year)) + geom_polygon(data=positions, aes(x=x, y=y, group=id), fill=NA, colour="black") +
  scale_fill_gradient(low="green", high="red", limits=c(0, 150), oob=scales::squish)

ggplot(data=conflict.summary %>% filter(year > 2012), aes(fill=n.conflicts)) +geom_map(aes(map_id=id), map=positions) + expand_limits(positions) +
  ggtitle("Number of Conflicts") + facet_wrap(vars(year)) + geom_polygon(data=positions, aes(x=x, y=y, group=id), fill=NA, colour="black") +
  scale_fill_gradient(low="green", high="red", limits=c(0, 150), oob=scales::squish)

# Plot n.victims
max(conflict.summary$n.victims) # 4790
min(conflict.summary$n.victims) # 1
ggplot(data=conflict.summary %>% filter(year <= 1994), aes(fill=n.victims)) +geom_map(aes(map_id=id), map=positions) + expand_limits(positions) +
  ggtitle("Number of Victims") + facet_wrap(vars(year)) + geom_polygon(data=positions, aes(x=x, y=y, group=id), fill=NA, colour="black") +
  scale_fill_gradient(low="green", high="red", limits=c(0, 4800), oob=scales::squish)

ggplot(data=conflict.summary %>% filter(year > 1994, year <= 2003), aes(fill=n.victims)) +geom_map(aes(map_id=id), map=positions) + expand_limits(positions) +
  ggtitle("Number of Victims") + facet_wrap(vars(year)) + geom_polygon(data=positions, aes(x=x, y=y, group=id), fill=NA, colour="black") +
  scale_fill_gradient(low="green", high="red", limits=c(0, 4800), oob=scales::squish)

ggplot(data=conflict.summary %>% filter(year > 2003, year <= 2012), aes(fill=n.victims)) +geom_map(aes(map_id=id), map=positions) + expand_limits(positions) +
  ggtitle("Number of Victims") + facet_wrap(vars(year)) + geom_polygon(data=positions, aes(x=x, y=y, group=id), fill=NA, colour="black") +
  scale_fill_gradient(low="green", high="red", limits=c(0, 4800), oob=scales::squish)

ggplot(data=conflict.summary %>% filter(year > 2012), aes(fill=n.victims)) +geom_map(aes(map_id=id), map=positions) + expand_limits(positions) +
  ggtitle("Number of Victims") + facet_wrap(vars(year)) + geom_polygon(data=positions, aes(x=x, y=y, group=id), fill=NA, colour="black") +
  scale_fill_gradient(low="green", high="red", limits=c(0, 4800), oob=scales::squish)

  # Lowering maximum
ggplot(data=conflict.summary %>% filter(year > 2003, year <= 2012), aes(fill=n.victims)) +geom_map(aes(map_id=id), map=positions) + expand_limits(positions) +
  ggtitle("Number of Victims") + facet_wrap(vars(year)) + geom_polygon(data=positions, aes(x=x, y=y, group=id), fill=NA, colour="black") +
  scale_fill_gradient(low="green", high="red", limits=c(0, 1000), oob=scales::squish)

ggplot(data=conflict.summary %>% filter(year > 2012), aes(fill=n.victims)) +geom_map(aes(map_id=id), map=positions) + expand_limits(positions) +
  ggtitle("Number of Victims") + facet_wrap(vars(year)) + geom_polygon(data=positions, aes(x=x, y=y, group=id), fill=NA, colour="black") +
  scale_fill_gradient(low="green", high="red", limits=c(0, 250), oob=scales::squish)


conflict.summary.national <- conflict.summary %>% group_by(year) %>% summarise(n.conflicts=sum(n.conflicts), n.victims=sum(n.victims))
conflict.summary.FARC <- conflict %>% filter(side_a == "FARC" | side_b == "FARC") %>% group_by(year) %>% summarise(n.conflicts=n(), n.victims=sum(n.victims))
par(mfrow=c(2,2))
{
a <- barplot(conflict.summary.national$n.conflicts, main="Number of Conflicts In Colombia")
axis(1, at=a, labels=conflict.summary.national$year)
c <- barplot(conflict.summary.FARC$n.conflicts, main="Number of Conflicts With FARC In Colombia")
axis(1, at=c, labels=conflict.summary.FARC$year)

b <- barplot(conflict.summary.national$n.victims, main="Number of Victims In Colombia")
axis(1, at=b, labels=conflict.summary.national$year)
d <- barplot(conflict.summary.FARC$n.victims, main="Number of Victims With FARC In Colombia")
axis(1, at=d, labels=conflict.summary.FARC$year)
}
par(mfrow=c(1,1))