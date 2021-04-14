# Calculating local Moran’s I
# Calculating local Moran’s I in R
# This code shows how to calculate local Moran’s I in R.
# The R package that contains the function that calculates Moran’s I is the spdep package. 
# In addition to the spdep package, other packages that can be useful are maps and maptools. 
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

# Reading in the SAT dataset.
sat.data <- read.table("state-sat.csv",sep=",",header=T)

# To look at the names of the variables in the SAT dataset, we can simply do the following:
names(sat.data)
## [1] "state"  "verbal" "math"   "pct"

# Here we are taking the individual columns in the SAT dataset, and create variables within the R environment
sat.state <- sat.data$state
sat.math <- sat.data$math
sat.continental <- sat.math[-c(which(as.character(sat.state)=="alaska"),which(as.character(sat.state)=="hawaii"))]


# The following command takes the shapefile for the states within the United States and creates a spatial polygon.
# Here we load the shapefile with the states in the US.
us.states <-map("state",fill=T,plot=F)
# This creates an ID for each US state.
us.IDs <- sapply(strsplit(us.states$names,":"),function(x) x[1])
# This creates a spatial polygon from a map object.
us.poly <- map2SpatialPolygons(us.states,IDs=us.IDs,proj4string=CRS("+proj=longlat +datum=WGS84"))
# This creates an nb object from a spatial polygon.
us.nb <- poly2nb(us.poly)
# This command creates a list with binary (since style="B") adjacency weights.
us.listw <- nb2listw(us.nb,style="B",zero.policy=TRUE)


# The command "localmoran" calculates local Moran's I.
# The function returns a matrix with multiple columns: the first column reports the values of local Moran's I, 
# the second
# column reports the expected value of local Moran's I would the data be independent, the third column reports 
# the variance of local Moran's I would the data be independent.
# The fourth column reports the Z-value associated with each areal unit and the fifth column 
# reports the p-value.
localmoran.i <- localmoran(sat.continental, us.listw,zero.policy=FALSE)
I.i <- localmoran.i[,1]
pval.Ii <- localmoran.i[,5]
## Plotting the local Moran's I for each state.
plotvar <- as.numeric(I.i)
nclr <- 5
plotclr <- brewer.pal(nclr,"YlGnBu")

summary(I.i)
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  -2.935   0.272   1.598   2.775   4.021  16.174
class <- classIntervals(plotvar,nclr,style="fixed",fixedBreaks=c(-3.0,0,1,2,3,16.5))
colcode <- findColours(class,plotclr)
plot(us.poly,border="black",axes=TRUE)
title(xlab="Longitude",ylab="Latitude",main="Local Moran's I \n Average SAT Math score")
plot(us.poly,col=colcode,add=T)
leg.txt<-c("[-0.75; 0)","[0; 1)","[1; 2)","[2; 3)","[3; 3.5]")
legend(-122,27,legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")

## Plot of the p-values of local Moran I in each US state.
plotvar <- as.numeric(pval.Ii)
nclr <- 5
plotclr <- brewer.pal(nclr,"RdGy")
class <- classIntervals(plotvar,nclr,style="fixed",fixedBreaks=c(0,0.01,0.05,0.1,0.5,1.0))
colcode <- findColours(class,plotclr)
plot(us.poly,border="black",axes=TRUE)
title(xlab="Longitude",ylab="Latitude",main="Local Moran's I \n P-values")
plot(us.poly,col=colcode,add=T)
leg.txt<-c("[0; 0.01)","[0.01; 0.05)","[0.05; 0.1)","[0.1; 0.5)","[0.5; 1.0]")
legend(-122,27,legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")
