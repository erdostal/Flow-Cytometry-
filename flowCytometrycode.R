library(ggplot2)
library(lattice)
library(ggmap)
library(dplyr)
library(broom)

#Import data set
flow <- read.csv("Hirsutum.csv", header = TRUE)


#Plot data unmodified 
plot (GenomeSize.mgb ~ Variety, data = flow,
     xlab = "Variety: Domesticated or Wild",
     ylab = "Corrected Genome size (mgb)",
     main = "Genome Size by Variety Type"
)

#Plot data unmodified 
plot(GenomeSize.pg ~ Variety, data = flow,
     xlab = "Variety: Domesticate or Wild",
     ylab = "Corrected Genome size (pg)",
     main = "Genome Size by Variety Type"
)

#Fit linear model to data
flow.modmgb <- lm(GenomeSize.mgb ~ Variety, data = flow)
summary(flow.modmgb)
flow.modpg <- lm(GenomeSize.pg ~ Variety, data = flow)
summary(flow.modpg)
#Plot linear model 
plot(resid(flow.modmgb) ~ fitted(flow.modmgb),
     xlab = "Fitted Values",
     ylab = "Residuals",
     main = "Residual Diagnostic Plot",
)
plot(resid(flow.modpg) ~ fitted(flow.modpg),
     xlab = "Fitted Values",
     ylab = "Residuals",
     main = "Residual Diagnostic Plot",
)

#Plot residuals to determine normality 
qqmath( ~ resid(flow.modmgb),
        xlab = "Theoretical Quantiles",
        ylab = "Residuals"
)
qqmath( ~ resid(flow.modpg),
        xlab = "Theoretical Quantiles",
        ylab = "Residuals"
)

#Test significance between variety and genome size
flowANOVA <- aov(GenomeSize.mgb ~ Variety, data = flow)
summary(flowANOVA)
flowANOVA <- aov(GenomeSize.pg ~ Variety, data = flow)
summary(flowANOVA)

#Chi squared, check this accuracy 
chisq.test(table(flow$GenomeSize.mgb, flow$Variety))
chisq.test(table(flow$GenomeSize.pg, flow$Variety))

#Test significance with two levels
t.test(GenomeSize.mgb ~ Variety, data = flow)

##### Create map showing location of accessions #####
map <- get_map(location = c(lon = -85.6024, lat = 12.7690), zoom = 3)

flow$Longitude <- flow$Longitude * -1

mapPoints <- ggmap(map) + 
  geom_point(aes(x = Longitude, y = Latitude, color = Variety.1), data = flow, alpha = .5)

### Find averages of each accessions
flow.averages <- aggregate(flow[, 8], list(flow$Accession), mean)
flow.SD <- aggregate(flow[, 8], list(flow$Accession), sd)
group <- c("Acala", "Coker_315", "DP90", "FM966", "Hopi", "TM1", "TX_6", "TX_12", "TX_44", "TX_109", "TX_180", "TX_279")
mean <- flow.averages$x
mean <-as.integer(mean)
SD <- flow.SD$x
SD <- as.integer(SD)
Variety <- c(1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
stats <- cbind(mean, SD, Variety)
SummaryStatistics <- as.matrix(SummaryStatistics)

Accessions <- read.csv("Accessions.csv")
TX_6SingleT <- t.test(Accessions$TX_0006, mu=2396)
TX_12SingleT <- t.test(Accessions$TX_0012, mu=2396)


###Finding how many reads in 1% of the genome ###
flow.averages$basepairs <- flow.averages$x * 1000000
flow.averages$NumberReads <- flow.averages$basepairs / 95
flow.averages$OnePercent <- flow.averages$NumberReads * .01
write.csv(flow.averages, "OnePercentReads.csv")

##################### REPEAT EXPLORER ##########################
Hcluster <- read.csv("HirsutumCluster.csv", header = TRUE)

library(ggplot2)
library(scales)
library(ggrepel)
library(factoextra)
library(reshape2)
library(gridExtra)
library(geomorph)
library(testit)

library(ape)
library(phytools)
library(phylotools)
library(geiger)

