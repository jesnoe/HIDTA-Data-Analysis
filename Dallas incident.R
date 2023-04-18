library(zoo)
library(readxl)
library(tidyverse)
dallas.incidents <- read.csv("Police_Incidents.csv") %>% as.tibble
names(dallas.incidents)
dallas.incidents$Type.of.Incident %>% unique
grep("DRUG|TRAFFICKING", unique(dallas.incidents$Type.of.Incident), value=T)# %>% write.csv("Dallas incident types.csv", row.names=F)
types.of.interest <- grep("MARIJUANA", unique(dallas.incidents$Type.of.Incident), value=T)
dallas.incidents.drug <- dallas.incidents %>% filter(Type.of.Incident %in% types.of.interest, City=="DALLAS")
# write.csv(dallas.incidents.drug, "Dallas Incidents-DRUG.csv", row.names=F)
dallas.incidents.drug <- read.csv("Dallas Incidents-DRUG.csv")
dallas.incidents.drug$Month1.of.Occurence <- match(dallas.incidents.drug$Month1.of.Occurence, month.name)
dallas.drug.occurence <- dallas.incidents.drug %>% filter(Year1.of.Occurrence > 2017, Year1.of.Occurrence < 2022) %>% 
  group_by(Year1.of.Occurrence, Month1.of.Occurence) %>% summarise(n.occurence=n())
year.month <- data.frame(Year1.of.Occurrence=rep(2018:2021, each=12),
                         Month1.of.Occurence=rep(1:12, 4))
dallas.drug.occurence <- merge(dallas.drug.occurence, year.month, all.y=T)

par(mar = c(5, 4, 4, 4) + 0.3)
plot(ts(dallas.drug.occurence$n.occurence, start=c(2018, 1), frequency=12),
     type="l",
     xlab="Year",
     ylab="# of Drug Incidents",
     pch=19)
par(new=T)
lines(ts(mortality.Dallas$Deaths, start=c(2018, 1), frequency=12),
     type="l",
     xlab="",
     ylab="",
     pch=19,
     col=2,
     axes=F)
axis(side = 4)
mtext("# of Deaths by Cocaine in Dallas", side = 4, line = 3)
legend("topleft", c("Drug Incidents", "# of Deaths"), cex=0.7, col=1:2, lty=1)


plot(seq(2018, 2022, length.out=48), dallas.drug.occurence$n.occurence,
     type="l",
     xlab="Year",
     ylab="# of Drug Incidents",
     pch=19)
par(new=T)
lines(seq(2018, 2022, length.out=48), mortality.Dallas$Deaths,
      type="l",
      xlab="",
      ylab="",
      pch=19,
      col=2)
axis(side = 4)
mtext("# of Deaths by Cocaine in Dallas", side = 4, line = 3)
legend("topleft", c("Drug Incidents", "# of Deaths"), cex=0.7, col=1:2, lty=1)
