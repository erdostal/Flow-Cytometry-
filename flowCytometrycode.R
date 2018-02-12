library(ggplot2)
library(lattice)
library(ggmap)
library(dplyr)
library(broom)


#Import data set
flow <- read.csv("Hirsutum.csv", header = TRUE)


#Plot data unmodified 
plot (GenomeSize.mgb ~ Variety, data = flow,
     xlab = "Variety",
     ylab = "Genome size (mgb)",
     main = "Genome Size by Variety Type"
)
plot (GenomeSize.mgb ~ Accession, data = flow,
      xlab = "Accession",
      ylab = "Genome size (mgb)",
      main = "Genome Size by Accession"
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
stats <- cbind(flow.averages$x,flow.SD$x)
stats <- as.data.frame(stats)
stats$CV <- stats$SD / stats$Average
Accessions <- unique(as.character(flow$Accession))
Variety <- c("Domesticated", "Domesticated", "Domesticated", "Domesticated", "Domesticated","Domesticated", 
             "Landrace", "Landrace", "Landrace", "Landrace", "Domesticated", "Landrace", "Landrace", "Landrace", 
             "Landrace", "Landrace", "Landrace", "Landrace", "Landrace", "Landrace", "Wild", "Landrace", "Wild", 
             "Landrace", "Landrace", "Landrace", "Landrace", "Landrace", "Landrace", "Landrace", "Landrace", 
             "Landrace", "Landrace", "Landrace", "Domesticated", "Landrace", "Landrace", "Wild", "Wild", 
             "Landrace", "Landrace")
stats$Accession <- Accessions
stats$Variety <- Variety
## Write new table with Standard Deviation and CV into a file
write.csv(stats, "HirsutumCV.csv")

# Plot Accession genome size averages with error bars
library(Hmisc)
yplus <- stats$Average + stats$SD
yminus <- stats$Average - stats$SD

errbar(stats$Accession, stats$Average, yplus, yminus)
plot(flow$Accession, flow$GenomeSize.mgb, )

## Find average/sd for wild/landrace/domesticated
Variety.averages <- aggregate(flow[, 8], list(flow$Variety), mean)
Variety.SD <- aggregate(flow[, 8], list(flow$Variety), sd)

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

