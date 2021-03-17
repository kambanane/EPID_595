#### EPID 595 R code set 1: Variograms ####


## First install all the packages that we need
install.packages("gstat",repos="https://cloud.r-project.org")

## 
## The downloaded binary packages are in
##  /var/folders/9n/gn30q9z56gv9qbdhfqjyfqn80000gn/T//RtmphyO1pV/downloaded_packages

## Get the gstat library 
library(gstat)

# Reading in the daily average PM2.5 concentration dataset presented in Unit 1
pm25.data <- read.table("Data/Avg_PM25_Feb2019_California_and_other_vars.csv",sep=",",header=T)
# To look at the names of the variables in the PM2.5 dataset, we can simply do the following:

ames(pm25.data)
##  [1] "Site.ID"          "Longitude"        "Latitude"         "UTM_Easting"     
##  [5] "UTM_Northing"     "Location_setting" "Avg_PM25"         "Avg_temperature" 
##  [9] "Average_ppt_amt"  "Density_pop"      "Elevation"
# Here we are taking the individual columns in the PM2.5 dataset, and creating variables within the R environment
pm25 <- pm25.data$Avg_PM25
temperature <- pm25.data$Avg_temperature
ppt <- pm25.data$Average_ppt_amt
x.coord <- pm25.data$UTM_Easting
y.coord <- pm25.data$UTM_Northing
lon.pm25 <- pm25.data$Longitude
lat.pm25 <- pm25.data$Latitude
pop.dens <- pm25.data$Density_pop
elev <- pm25.data$Elevation

# The following command obtains summary statistics of PM2.5 concentration and creates a histogram with overlaid the estimated density
summary(pm25)
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   1.337   3.379   4.778   5.234   6.687  16.056
# Here we calculate the histogram with probabilities shown on the vertical axis
hist(pm25,prob=T,col="coral",xlab="Average PM2.5 concentration",main="Average PM2.5 concentration \n in California, February 2019",breaks=30)
# Here we are estimating the density
dens.pm25 <- density(pm25)
# This displays the density for PM2.5 concentration as an overlaid line on the histogram
lines(dens.pm25,col="black",lwd=2,lty=1)

# The following command creates a spatial plot of the daily average PM2.5 concentration. Similar commands can be used to generate spatial plots of the covariates - i.e. daily average temperature, rainfall amount, population density, and elevation - used to model daily average PM2.5 concentration.
# To create spatial maps, we use the function image.plot that is contained in the R package fields.
# Since we want to add the border of California to the map, we also need to install and load two additional packages: maps and maptools.
# Hence:
install.packages("fields",repos="https://cloud.r-project.org")
## 
## The downloaded binary packages are in
##  /var/folders/9n/gn30q9z56gv9qbdhfqjyfqn80000gn/T//RtmphyO1pV/downloaded_packages
install.packages("maps",repos="https://cloud.r-project.org")
## 
## The downloaded binary packages are in
##  /var/folders/9n/gn30q9z56gv9qbdhfqjyfqn80000gn/T//RtmphyO1pV/downloaded_packages
install.packages("maptools",repos="https://cloud.r-project.org")
## 
## The downloaded binary packages are in
##  /var/folders/9n/gn30q9z56gv9qbdhfqjyfqn80000gn/T//RtmphyO1pV/downloaded_packages
library(fields)
## Loading required package: spam
## Loading required package: dotCall64
## Loading required package: grid
## Spam version 2.6-0 (2020-12-14) is loaded.
## Type 'help( Spam)' or 'demo( spam)' for a short introduction 
## and overview of this package.
## Help for individual functions is also obtained by adding the
## suffix '.spam' to the function name, e.g. 'help( chol.spam)'.
## 
## Attaching package: 'spam'
## The following objects are masked from 'package:base':
## 
##     backsolve, forwardsolve
## See https://github.com/NCAR/Fields for
##  an extensive vignette, other supplements and source code
library(maps)
library(maptools)
## Loading required package: sp
## Checking rgeos availability: TRUE
# This determines the interval used to display the color scale
zlim <- round(range(pm25,na.rm=T),1)
# To create a spatial map, we use the command image.plot. The function image.plot asks that the data be provided in the following format: an increasing vector of coordinates for the horizontal axis (e.g. longitudes or Eastings) separated by regular intervals, an increasing vector of coordinates for the vertical axis (e.g. latitudes or Northings) separated by regular intervals, and a matrix of dimensions, the length of the vector of coordinates on the horizontal axis times the length of the vector of coordinates on the vertical axis.
# If the location of the observations are not arranged on a regular grid, we fill this matrix with values outside the zlim values. This will force R to create just a box for the spatial domain, filled with blank values.
# Here we create a vector of length 200 of regularly spaced longitudes and latitudes spanning the range of longitudes and latitudes of the PM2.5 monitors.
lon.rg <- range(lon.pm25)
lon.rg
## [1] -124.1795 -115.4831
lon.grid <- seq(lon.rg[1], lon.rg[2], len = 200)

lat.rg <- range(lat.pm25)
lat.rg
## [1] 32.63124 41.72689
lat.grid <- seq(lat.rg[1], lat.rg[2], len = 200)

# Here we create the matrix of dimension 200 by 200 filled with values of -1 (since we know PM2.5 concentration has to be positive)
z.mat <- matrix(-1,200,200)
# This command assigns to each observed PM2.5 concentration a color
col.pm25 <- color.scale(pm25)
# This command actually generates the map filled with blank values
image.plot(lon.grid, lat.grid,z.mat,zlim=zlim,xlab="Longitude",ylab="Latitude",main="Daily average PM2.5 concentration \n February 2019",xlim=c(-125,-114),ylim=c(32,43))
# This adds the California boundaries
map("state", add = T, col = "gray")
# This adds to the map the actual point locations colored according to the observed PM2.5 concentration at the location. Note that here to generate the plot, we use the latitude and longitude coordinates.
points(lon.pm25,lat.pm25,col=col.pm25,pch=19)

## This command constructs an empirical semi-variogram: the syntax requests the response variable on the left hand side of the ~ sign, while on the right hand side, we have to list all the predictors we want to include in the linear regression model. The residuals of the linear regression model is what we use to construct the empirical semi-variogram.
# In order for R to be able to calculate distances between observation sites and create the distance bins, R needs to be provided with the geographical coordinates of the locations of the observations. These coordinates are specified after the option "locations=~".
# The variogram function has a default option as for how to define the distance bins in an empirical semi-variogram. A user can change the definition of the bins by specifying a vector of increasing distance bins upper boundaries using the option "boundaries" in the function variogram.
# To run the variogram function, the user need to create first a dataframe containing all the needed information used when running the function variogram.
pm25.df <- data.frame(cbind(x.coord,y.coord,pm25,temperature,ppt,pop.dens,elev))
emp.variog.pm25 <- variogram(pm25~temperature+ppt+pop.dens+elev,locations=~x.coord+y.coord,data=pm25.df)
# The function variogram reports, the number of pairs that fall within each bin (np), the value of the semi-variance for each bin (gamma), and each bin's representative, or the bins' distance midpoints (dist).
emp.variog.pm25
##     np      dist    gamma dir.hor dir.ver   id
## 1  114  17829.22 2.170668       0       0 var1
## 2  177  42708.42 3.752572       0       0 var1
## 3  254  71327.45 4.364023       0       0 var1
## 4  298 100013.65 3.740616       0       0 var1
## 5  361 128391.32 5.367459       0       0 var1
## 6  383 156508.50 5.987330       0       0 var1
## 7  320 183865.75 5.589148       0       0 var1
## 8  295 212837.82 5.125567       0       0 var1
## 9  317 241715.64 6.308610       0       0 var1
## 10 321 269801.75 7.717730       0       0 var1
## 11 287 297441.66 5.743690       0       0 var1
## 12 243 326845.24 6.559499       0       0 var1
## 13 257 355195.50 5.385841       0       0 var1
## 14 219 382693.72 7.209512       0       0 var1
## 15 227 412235.45 4.390867       0       0 var1
# To plot an empirical semi-variogram, we simply use the command plot on the output of the function variogram.
plot(emp.variog.pm25)