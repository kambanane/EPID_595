# Assessing spatial correlation for areal data
# Assessing spatial correlation for areal data in R
# This code shows how to calculate Moran’s I and create a correlogram in R.
# The R package that is needed for the analysis of areal data is the spdep package. 
# In addition to the spdep package, the packages maps and maptools are also useful 
# in that they contains shapefiles for counties and states in the United States and 
# functions that are useful to transform the information in shapefiles into quantities 
# that are needed for the analysis of areal data (e.g. notion of adjacency). 
# All these packages can be downloaded from a CRAN repository.
# We also use the packages RColorBrewer and classInt to make chloropleth maps in R.
# knitr::opts_chunk$set(fig.width=8, fig.height=8) 


# Installing the necessary packages
install.packages(c("maps","maptools","spdep","RColorBrewer","classInt"),repos="https://cloud.r-project.org")

library(maps)
library(maptools)
library(spdep)
library(RColorBrewer)
library(classInt)



# Reading in the SAT dataset
sat.data <- read.table("state-sat.csv",sep=",",header=T)
# To look at the names of the variables in the SAT dataset, we can simply do the following:
names(sat.data)
## [1] "state"  "verbal" "math"   "pct"
# Here we are taking the individual columns in the SAT dataset, and create variables within the R environment.
sat.state <- sat.data$state
sat.math <- sat.data$math


# The following command takes the shapefile for the states within the United States and creates a spatial polygon.
# Here we load the shapefile with the states in the US.
us.states <-map("state",fill=T,plot=F)
# This creates an ID for each US state
us.IDs <- sapply(strsplit(us.states$names,":"),function(x) x[1])
# This creates a spatial polygon from a map object.
us.poly <- map2SpatialPolygons(us.states,IDs=us.IDs,proj4string=CRS("+proj=longlat +datum=WGS84"))

# Creating a chloropleth map of the Average SAT Math score for all the 
# United states except Alaska and Hawaii.
# Here we define the variable that we want to create a plot for:
sat.continental <- sat.math[-c(which(as.character(sat.state)=="alaska"),which(as.character(sat.state)=="hawaii"))]
plotvar <- sat.continental
# Here we define how many colors we want to use for the map.
nclr <- 5
# This command defines the 5 different colors to use for the map within 
# the blue palette.
plotclr <- brewer.pal(nclr,"Blues")
# This command assigns each observation to one of 5 intervals 
# represented by one of 5 different shades of blue. The 5 intervals that 
# we are using to assign each observation to a color are: [475, 500), 
# [500, 525), [525, 550), [550, 575), and [575,605].
class <- classIntervals(plotvar,nclr,style="fixed",fixedBreaks=c(475,500,525,550,575,605))
# This command assigns one the 5 colors to each of the observations based
# on the interval in which the observation falls.
colcode <- findColours(class,plotclr)

# This creates a map of the US with borders for each state, but the 
# states are not filled with colors.
plot(us.poly,border="black",axes=TRUE)
title(xlab="Longitude",ylab="Latitude",main="Average SAT Math score")
# This command fills the map with colors. For this command to work, it is 
# important that the spatial polygon lists the US states in the same 
# order as in the dataset with the SAT data.
plot(us.poly,col=colcode,add=T)

# This command adds a legend to the plot.
leg.txt<-c("[475; 500)","[500; 525)","[525; 550)","[550; 575)","[575; 605]")
legend(-122,27,legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")

# Here we look at the distribution of the Average SAT Math score.
hist(sat.continental,breaks=20,col="coral",xlab="Average SAT Math score",freq=F,main="Histogram of Avg SAT Math score")
lines(density(sat.continental),lwd=2,col="black")

# Here we calculate Moran's I. The function that derives Moran's I is the function
# moran.test in the spdep package.
# To calculate Moran's I, it is necessary to have defined the adjacency weights
# w_{ij}. These weights need to be passed in a list form.

# The function poly2nb takes a spatial polygon object and creates an nb object.
# The default notion of adjacency used by the poly2nb function is the Queen's
# notion of adjacency. If "queen" is set equal to FALSE, then the notion 
# of adjacency used is Rook's.
us.nb <- poly2nb(us.poly)
# The nb does not provide much information except for overall summary
# statistics. 
us.nb
## Neighbour list object:
## Number of regions: 49 
## Number of nonzero links: 218 
## Percentage nonzero weights: 9.07955 
## Average number of links: 4.44898
# The nb2listw transforms the nb object into a list with 3 components with
# the weights. If style is "B", the adjacency weights w_{ij} are binary 
# weights, while if style is "W", the binary adjacency weights w_{ij}
# for each unit i are divided by the number of neighbors areal unit i has.
# The option zero.policy=TRUE is to allow the calculation of the adjacency
# weights even in case some areal units have no neighbors.
us.listw <- nb2listw(us.nb,style="B",zero.policy=TRUE)

# The following command calculated Moran's I. The option rank is to 
# distinguish situations where the variable for which we want to calculate
# Moran's I is continuous (in this case, rank=FALSE) vs discrete (in this 
# case, rank=TRUE). The option zero.policy=TRUE is to allow calculations 
# of Moran's I when there are areal units with no neighbors.
# The function moran.test also reports a p-value for testing whether 
# Moran's I is equal to what it would be under the assumption of independence.
moran.i.sat <- moran.test(sat.continental, us.listw, rank=FALSE,zero.policy=TRUE)
moran.i.sat
## 
##  Moran I test under randomisation
## 
## data:  sat.continental  
## weights: us.listw    
## 
## Moran I statistic standard deviate = 7.0926, p-value = 6.582e-13
## alternative hypothesis: greater
## sample estimates:
## Moran I statistic       Expectation          Variance 
##       0.623652403      -0.020833333       0.008256943
# This code shows how to calculate a correlogram. 
# The function that calculates the different Moran's I at various order of neighborhood is 
# the function sp.correlogram in the package spdep.
# The function requires as input
# (i) an object of the class nb that lists which areal unit is adjacent to which
# (ii) "var": the name of the variable for which to calculate the correlogram
# (iii) "order": the order up to which calculate the spatial statistics
# (iv) "method": the type of spatial statistics to calculate for the correlogram: "I" is for Moran's I
# (v) "style": the type of adjacency weights to use. "B" is for binary adjacency weights, while 
# "W" is for adjacency weights obtained by taking the binary adjacency weights divided by the number
# of neighbors each areal unit has.
# (vi) "randomisation": this option determines whether the standard errors on Moran's I (if "method" equal "I") 
# should be calculated using a bootstrap approach (if "randomisation" is equal to TRUE), or using 
# the asymptotic distribution of Moran's I.
# (vii) "zero.policy": this option instruct R how to proceed in case some areal units do not have a neighbor.
# If the zero.policy is TRUE, the command will still compute Moran's I even in the case of areal units without 
# neighbors.
correl.sat <- sp.correlogram(us.nb,sat.continental, order = 5, method = "I",
                             style = "B", randomisation = FALSE, zero.policy = TRUE)
## This command is to create the correlogram plot
plot(correl.sat, main="Correlogram for Avg SAT MATh score", xlab="Order",ylab="Moran's I")
