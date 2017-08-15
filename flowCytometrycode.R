library(ggplot2)
library(lattice)
library(ggmap)

#Import data set
flow <- read.csv("FlowCytometery.csv", header = TRUE)

#Plot data unmodified 
plot(Corrected.Genome ~ Variety, data = flow,
     xlab = "Variety: Cultivar, Obsolete Cultivar, Wild",
     ylab = "Genome Size by Variety Type",
     main = "Corrected Genome size (picograms)"
)

#Plot data unmodified 
plot(Corrected.Genome ~ Variety.1, data = flow,
     xlab = "Variety: Domesticate or Wild",
     ylab = "Corrected Genome size (picograms)",
     main = "Genome Size by Variety Type"
)

#Fit linear model to data
flow.mod <- lm(Corrected.Genome ~ Variety, data = flow)
summary(flow.mod)
#Plot linear model 
plot(resid(flow.mod) ~ fitted(flow.mod),
     xlab = "Fitted Values",
     ylab = "Residuals",
     main = "Residual Diagnostic Plot",
)

#Plot residuals to determine normality 
qqmath( ~ resid(flow.mod),
        xlab = "Theoretical Quantiles",
        ylab = "Residuals"
)

#Test significance between variety and genome size
flowANOVA <- aov(Corrected.Genome ~ Variety, data = flow)
summary(flowANOVA)

#Chi squared, check this accuracy 
chisq.test(table(flow$Corrected.Genome, flow$Variety))

#Test significance with two levels
t.test(Corrected.Genome ~ Variety.1, data = flow)

map <- get_map(location = c(lon = -85.6024, lat = 12.7690), zoom = 3)

flow$Longitude <- flow$Longitude * -1

mapPoints <- ggmap(map) + 
  geom_point(aes(x = Longitude, y = Latitude, color = Variety.1), data = flow, alpha = .5)



