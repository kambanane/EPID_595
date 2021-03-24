# Performing spatial interpolation in R
# Performing spatial interpolation of point referenced data in R
# This code shows how to estimate a spatial variable (e.g. an environmental exposure) at a location without observations.
# Specifically the code will implement the methods of:(i) nearest monitor
# (ii) inverse distance weighting (IDW)
# and (iii) kriging:Ordinary and Universal.
# The R package that are needed are:fields (for calculating distances and mapping), gstat (for the empirical semi -
#                                                                                            variogram and fitting a parametric semi - variogram to it), geoR (for kriging).
# We will use as example the daily average PM2.5 concentration in California in February 2019.


#install.packages(c("fields","gstat","geoR"),repos="https://cloud.r-project.org")
## 
##   There is a binary version available but the source version is later:
##       binary source needs_compilation
## gstat  2.0-6  2.0-7              TRUE
## 
## 
## The downloaded binary packages are in
##  /var/folders/9n/gn30q9z56gv9qbdhfqjyfqn80000gn/T//RtmpSmibXl/downloaded_packages
## installing the source package 'gstat'
## Warning in install.packages(c("fields", "gstat", "geoR"), repos = "https://
## cloud.r-project.org"): installation of package 'gstat' had non-zero exit status

library(fields)
library(gstat)
library(geoR)

## --------------------------------------------------------------
##  Analysis of Geostatistical Data
##  For an Introduction to geoR go to http://www.leg.ufpr.br/geoR
##  geoR version 1.8-1 (built on 2020-02-08) is now loaded
## --------------------------------------------------------------

# Reading in the daily average PM2.5 concentration dataset 
pm25.data <- read.table("Data/Avg_PM25_Feb2019_California_and_other_vars.csv",sep=",",header=T)

pm25 <- pm25.data$Avg_PM25
temperature <- pm25.data$Avg_temperature
ppt <- pm25.data$Average_ppt_amt
x.coord <- pm25.data$UTM_Easting
y.coord <- pm25.data$UTM_Northing
lon.pm25 <- pm25.data$Longitude
lat.pm25 <- pm25.data$Latitude
pop.dens <- pm25.data$Density_pop
elev <- pm25.data$Elevation


# This is the fictitious dataset with the residential addresses of the babies for which we have birthweight data information. Additionally, the dataset provides the values of the explanatory variables (daily average temperature, daily avg rainfall amount, population density and elevation) at the birthweight locations.
bw.data <- read.csv("Data/Birthwt_locs_Cali.csv",sep=",",header=T)
names(bw.data)
## [1] "longitude"   "latitude"    "temperature" "ppt"         "pop.density"
## [6] "elevation"
bw.lon <- bw.data$longitude
bw.lat <- bw.data$latitude
bw.temperature <- bw.data$temperature
bw.ppt <- bw.data$ppt
bw.pop.dens <- bw.data$pop.density
bw.elev <- bw.data$elevation


# Taking a random location in the birtweight dataset for which we are going to estimate the environmental exposure, e.g. the daily average PM2.5 concentration.
index.loc <- 36


### Method of nearest monitor
# For this method, we need to first calculate the distance between the selected location and the environmental exposure locations.
Dist.pm25.bw <- rdist.earth(matrix(cbind(lon.pm25,lat.pm25),nrow=length(lon.pm25),ncol=2),matrix(c(bw.lon[index.loc],bw.lat[index.loc]),nrow=1,ncol=2),miles=F)
dim(Dist.pm25.bw)
## [1] 114   1
# This will return which entry in the vector with the coordinates of the environmental exposure locations has the shortest distance from the selected location.
which.min(Dist.pm25.bw)
## [1] 76
# The estimated environmental exposure at the selected location according to the method of nearest monitor is:
pm25[which.min(Dist.pm25.bw)]
## [1] 5.045757
### Inverse distance weighting 
# To calculated the estimated environmental exposure at the selected location according to the method of IDW, we need to take a weighted average of the environmental exposure data where the weights are inversely proportional to the distance the selected sites has to each environmental exposure location. The weights are then normalized so that they add up to 1. Once the weights have been normalized, the estimated environmental exposure at the selected location is a weighted average of the environmental exposure data with normalixed weights.

# Calculating the weights as the inverse distance b/w the selected location and the environmental exposure location
idw.weights.loc <- 1/as.numeric(as.numeric(Dist.pm25.bw))
# This takes care of the normalization
idw.weights.norm.loc <- idw.weights.loc/sum(idw.weights.loc)
# The estimated exposure is then given by
idw.pm25.sel.loc <- sum(idw.weights.norm.loc*pm25)
idw.pm25.sel.loc
## [1] 5.853388
# Here we perform kriging.
# To this goal, we first  construct the empirical semi-variogram for the residuals of the multiple linear regression model regressing daily average PM2.5 concentration on the covariates, i.e. daily average temperature, rainfall amount, population density and elevation.
# Then we fit to it an parametric semi-variogram (here, an exponential) to the empirical semi-variogram.
# And we finally we run a spatial linear regression model to it.
# We start by transforming Easting and Northing from  meters to kilometers.
x.coord.km <- x.coord/1000
y.coord.km <- y.coord/1000
pm25.df <- data.frame(cbind(x.coord.km,y.coord.km,pm25,temperature,ppt,pop.dens,elev))
emp.variog.pm25 <- variogram(pm25~temperature+ppt+pop.dens+elev,locations=~x.coord.km+y.coord.km,data=pm25.df)
emp.variog.pm25
##     np      dist    gamma dir.hor dir.ver   id
## 1  114  17.82922 2.170668       0       0 var1
## 2  177  42.70842 3.752572       0       0 var1
## 3  254  71.32745 4.364023       0       0 var1
## 4  298 100.01365 3.740616       0       0 var1
## 5  361 128.39132 5.367459       0       0 var1
## 6  383 156.50850 5.987330       0       0 var1
## 7  320 183.86575 5.589148       0       0 var1
## 8  295 212.83782 5.125567       0       0 var1
## 9  317 241.71564 6.308610       0       0 var1
## 10 321 269.80175 7.717730       0       0 var1
## 11 287 297.44166 5.743690       0       0 var1
## 12 243 326.84524 6.559499       0       0 var1
## 13 257 355.19550 5.385841       0       0 var1
## 14 219 382.69372 7.209512       0       0 var1
## 15 227 412.23545 4.390867       0       0 var1
plot(emp.variog.pm25)

# Then, we fit an exponential semi-variogram to it via WLS
variog.pm25 <- fit.variogram(emp.variog.pm25,vgm(psill=4.0,"Exp",300,2.0))
est.nugget <- variog.pm25$psill[1]
est.nugget
## [1] 1.048932
est.psill <- variog.pm25$psill[2]
est.psill
## [1] 4.944584
est.range <- variog.pm25$range[2]
est.range
## [1] 66.25919
# Finally, we fit a spatial linear regression model to the daily average PM2.5 concentration.
# The mean model includes the explanatory variables: daily average temperature, daily average rainfall amount, population density and elevation.
# The variance model is given by the exponential covariance function.
# The function we use to fit a spatial linear regression model via ML or REML is the function likfit.
pm25.df <- data.frame(cbind(x.coord.km,y.coord.km,pm25,temperature,ppt,pop.dens,elev),nrow=length(x.coord.km),ncol=7)
pm25.geo <- as.geodata(pm25.df,coords.col=c(1,2),data.col=3,covar.cols=c(4,7))
pm25.reml <- likfit(pm25.geo, trend = ~temperature+ppt+pop.dens+elev,cov.model="exponential",ini=c(est.psill,est.range),nugget=est.nugget, fix.nug = FALSE, lik.met="REML")
## kappa not used for the exponential correlation function
## ---------------------------------------------------------------
## likfit: likelihood maximisation using the function optim.
## likfit: Use control() to pass additional
##          arguments for the maximisation function.
##         For further details see documentation for optim.
## likfit: It is highly advisable to run this function several
##         times with different initial values for the parameters.
## likfit: WARNING: This step can be time demanding!
## ---------------------------------------------------------------
## likfit: end of numerical maximisation.
# Having fitted the spatial linear regression model to the data, we can now perform kriging and estimate PM2.5 concentration at the selected location.
# In order to do this, we first have to transform the coordinate of the selected location into Easting and Northing in km.
# To transform latitude and longitude into Easting and Northing, we use the packages sp and rgdal.

## Make sure you have these packages
#install.packages(c("sp","rgdal"),repos="https://cloud.r-project.org")

## 
## The downloaded binary packages are in
##  /var/folders/9n/gn30q9z56gv9qbdhfqjyfqn80000gn/T//RtmpSmibXl/downloaded_packages
library(sp)
library(rgdal)
## rgdal: version: 1.5-23, (SVN revision 1121)
## Geospatial Data Abstraction Library extensions to R successfully loaded
## Loaded GDAL runtime: GDAL 3.1.4, released 2020/10/20
## Path to GDAL shared files: /Library/Frameworks/R.framework/Versions/4.0/Resources/library/rgdal/gdal
## GDAL binary built with GEOS: TRUE 
## Loaded PROJ runtime: Rel. 6.3.1, February 10th, 2020, [PJ_VERSION: 631]
## Path to PROJ shared files: /Library/Frameworks/R.framework/Versions/4.0/Resources/library/rgdal/proj
## Linking to sp version:1.4-5
## To mute warnings of possible GDAL/OSR exportToProj4() degradation,
## use options("rgdal_show_exportToProj4_warnings"="none") before loading rgdal.
SP.longlat.bw <- SpatialPoints(coords=cbind(as.numeric(bw.lon),as.numeric(bw.lat)),proj4string=CRS("+proj=longlat +datum=WGS84"))
# This transforms longitude and latitude in UTM coordinates
SP.utm.bw <- spTransform(SP.longlat.bw,CRS("+proj=utm +zone=10 +datum=WGS84"))
bw.x <- coordinates(SP.utm.bw)[,1]
bw.y <- coordinates(SP.utm.bw)[,2]
bw.x.km <- bw.x/1000
bw.y.km <- bw.y/1000
### Ordinary Kriging
# The function that actually performs kriging is the function krige.conv in geoR.
# The function that specifies the type of kriging and provides all the necessary specification and information to implement kriging is the function krige.control.

# In the function krige.control, we need to specify:
# (i) type of kriging (type.krige): for both ordinary kriging and universal kriging, here we should type "ok";
# (ii) the mean model at the locations where we have the environmental exposure data (trend.d): for ordinary kriging, the mean model is a equal to an unknown, but constant value beta0. For ordinary kriging, we write "cte" in trend.d;
# (iii) the mean model at the location where we want to make predictions (trend.l): for ordinary kriging, this is again "cte";
# (iv) the variance model (cov.model): here we express the name of the semi-variogram model or covariance function model to fit to the data. Choices include "exponential" or "gaussian" (among others).
# (v) the parameters of the covariance, or semi-variogram, function (cov.pars). This is a 2-dimensional vector that should provide in order the value of the partial sill (o sigma^2) and of the range parameter (or phi).
# (vi) that value of the nugget effect (nugget). This should be either a variable or a number that provides the estimated nugget effect.
kc.ok.control <- krige.control(type.krige="ok",trend.d="cte",trend.l="cte",cov.model="exponential",cov.pars=c(pm25.reml$cov.pars[1],pm25.reml$cov.pars[2]),nugget=pm25.reml$nugget)

# Having provided all the input for the type of kriging, this is the command that actually performs kriging.
# First we place the locations where we want to generate estimates of the environmental exposure into an n by 2 matrix.
loc.kriging <- matrix(c(bw.x[index.loc],bw.y[index.loc]),nrow=1,ncol=2)
kc.ok.s0 <- krige.conv(pm25.geo,locations=loc.kriging,krige=kc.ok.control)
## krige.conv: model with constant mean
## krige.conv: Kriging performed using global neighbourhood
names(kc.ok.s0)
## [1] "predict"      "krige.var"    "beta.est"     "distribution" "message"     
## [6] "call"
# The predicted environmental exposure is stored into the part called "predict", while the kriging variance is what is stored into the part called "krige.var".
kc.ok.s0$predict
##     data 
## 5.021546
kc.ok.s0$krige.var
## [1] 6.58312
# We can create a 95% CI for the estimated environmental exposure by adding and subtracting 1.96*square root of the kriging variance to then estimated environmental expsure.
ci.95.kc.ok <-  c(kc.ok.s0$predict-1.96*sqrt(kc.ok.s0$krige.var),kc.ok.s0$predict+1.96*sqrt(kc.ok.s0$krige.var))
ci.95.kc.ok    
##         data         data 
## -0.007341619 10.050434480
# If the lower bound of the confidence interval for the selected location is less than 0, we transform the lower bound of the interval to 0, the lowest possibile value for the PM2.5 concentration.
### Universal kriging
# To perform universal kriging, the krige.control function needs to change from the form it had under ordinary kriging.
# Specifically, we need to change the form of the mean model at the environmental exposure data locations and at the locations where we want to estimate the environmental exposure.
# Hence, now:
# trend.d needs to provide the form of the mean model at the locations where we have observed the data: this is specified by "~" followed by the name of the explanatory variables at the environmental exposure locations. 
# trend.l needs to provide the form of the mean model where we want to estimate environmental exposure. The format is still "~" followd by the name of the explanatory variables at the location(s) where we want to predict the environmental exposure.
# The remainder of the specifications for the krige.control are exactly the same as under ordinary kriging.
kc.uk.control <- krige.control(type.krige="ok",trend.d=~temperature+ppt+pop.dens+elev,trend.l=~bw.temperature[index.loc]+bw.ppt[index.loc]+bw.pop.dens[index.loc]+bw.elev[index.loc],cov.model="exponential",cov.pars=c(pm25.reml$cov.pars[1],pm25.reml$cov.pars[2]),nugget=pm25.reml$nugget)

loc.kriging <- matrix(c(bw.x[index.loc],bw.y[index.loc]),nrow=1,ncol=2)
kc.uk.s0 <- krige.conv(pm25.geo,locations=loc.kriging,krige=kc.uk.control)
## krige.conv: model with mean defined by covariates provided by the user
## krige.conv: Kriging performed using global neighbourhood
names(kc.uk.s0)
## [1] "predict"      "krige.var"    "beta.est"     "distribution" "message"     
## [6] "call"
# The predicted environmental exposure is stored into the part called "predict", while the kriging variance is what is stored into the part called "krige.var".
kc.uk.s0$predict
##     data 
## 5.469714
kc.uk.s0$krige.var
## [1] 6.899578
# A 95% CI for the estimated environmental exposure by adding and subtracting 1.96*square root of the kriging variance to then estimated environmental expsure.
ci.95.kc.uk <-  c(kc.uk.s0$predict-1.96*sqrt(kc.uk.s0$krige.var),kc.uk.s0$predict+1.96*sqrt(kc.uk.s0$krige.var))
ci.95.kc.uk    
##       data       data 
##  0.3213726 10.6180557