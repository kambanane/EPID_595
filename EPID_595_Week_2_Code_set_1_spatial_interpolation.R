
# Performing spatial interpolation in R
# Performing spatial interpolation of point referenced data in R
# This code shows how to estimate a spatial variable (e.g. an environmental exposure) at a location without observations.
# Specifically the code will implement the methods of:
# (i) nearest monitor
# (ii) inverse distance weighting (IDW)
# and (iii) kriging:Ordinary and Universal.

# The R package that are needed are:fields (for calculating distances and mapping), gstat (for the empirical semi -
# variogram and fitting a parametric semi - variogram to it), geoR (for kriging).
# We will use as example the daily average PM2.5 concentration in California in February 2019.


#install.packages(c("fields","gstat","geoR"),repos="https://cloud.r-project.org")

library(fields)
library(gstat)
library(geoR)

# Reading in the daily average PM2.5 concentration dataset
pm25.data <-
  read.table(
    "Data/Avg_PM25_Feb2019_California_and_other_vars.csv",
    sep = ",",
    header = T
  )

pm25 <- pm25.data$Avg_PM25
temperature <- pm25.data$Avg_temperature
ppt <- pm25.data$Average_ppt_amt
x.coord <- pm25.data$UTM_Easting
y.coord <- pm25.data$UTM_Northing
lon.pm25 <- pm25.data$Longitude
lat.pm25 <- pm25.data$Latitude
pop.dens <- pm25.data$Density_pop
elev <- pm25.data$Elevation

# This is the fictitious dataset with the residential addresses of the babies for which we have birthweight data information. 
## Additionally, the dataset provides the values of the explanatory variables (daily average temperature, daily avg rainfall amount, population density and elevation) at the birthweight locations.
bw.data <- read.csv("Data/Birthwt_locs_Cali.csv",
                    sep = ",",
                    header = T)

names(bw.data)

### Make a map 
bw.data.sp <- bw.data
library(sp)
coordinates(bw.data.sp) <- ~longitude+latitude

pm25.data.sp <- pm25.data
coordinates(pm25.data.sp) <- ~Longitude+Latitude

plot(bw.data.sp, pch=16)
plot(pm25.data.sp, pch=16, col="blue", add=TRUE)

bw.lon <- bw.data$longitude
bw.lat <- bw.data$latitude
bw.temperature <- bw.data$temperature
bw.ppt <- bw.data$ppt
bw.pop.dens <- bw.data$pop.density
bw.elev <- bw.data$elevation


# Taking a random location in the birthweight dataset for which we are going to estimate the environmental exposure, 
# e.g. the daily average PM2.5 concentration.
index.loc <- 36

#### Method of nearest monitor ####

# For this method, we need to first calculate the distance between the selected location and the 
# environmental exposure locations.
Dist.pm25.bw <-
  rdist.earth(matrix(
    cbind(lon.pm25, lat.pm25),
    nrow = length(lon.pm25),
    ncol = 2
  ),
  matrix(c(bw.lon[index.loc], bw.lat[index.loc]), nrow = 1, ncol = 2),
  miles = F)
dim(Dist.pm25.bw)

## [1] 114   1
# This will return which entry in the vector with the coordinates of the environmental 
#exposure locations has the shortest distance from the selected location.

which.min(Dist.pm25.bw)
## [1] 76

# The estimated environmental exposure at the selected location according to the 
# method of nearest monitor is:
pm25[which.min(Dist.pm25.bw)]
## [1] 5.045757

#### Inverse distance weighting ####

# To calculate the estimated environmental exposure at the selected location according to the method of IDW, 
# we need to take a weighted average of the environmental exposure data where the weights are 
# inversely proportional to the distance the selected sites has to each environmental exposure location. 
# The weights are then normalized so that they add up to 1. Once the weights have been normalized, 
# the estimated environmental exposure at the selected location is a weighted average of the 
# environmental exposure data with normalized weights.

# Calculating the weights as the inverse distance b/w the selected location and the environmental exposure location
idw.weights.loc <- 1 / as.numeric(as.numeric(Dist.pm25.bw))

# This takes care of the normalization
idw.weights.norm.loc <- idw.weights.loc / sum(idw.weights.loc)

# The estimated exposure is then given by
idw.pm25.sel.loc <- sum(idw.weights.norm.loc * pm25)
idw.pm25.sel.loc

## [1] 5.853388

#### KRIGING ####

# Here we perform kriging.
# To this goal, we first  construct the empirical semi-variogram for the residuals 
# of the multiple linear regression model regressing daily average PM2.5 concentration on the covariates, 
# i.e. daily average temperature, rainfall amount, population density and elevation.
# Then we fit to it a parametric semi-variogram (here, an exponential) to the empirical semi-variogram.
# And we finally we run a spatial linear regression model to it.

# We start by transforming Easting and Northing from  meters to kilometers.

x.coord.km <- x.coord / 1000
y.coord.km <- y.coord / 1000

pm25.df <-
  data.frame(cbind(x.coord.km, y.coord.km, pm25, temperature, ppt, pop.dens, elev))

emp.variog.pm25 <-
  variogram(
    pm25 ~ temperature + ppt + pop.dens + elev,
    locations =  ~ x.coord.km + y.coord.km,
    data = pm25.df
  )

emp.variog.pm25

plot(emp.variog.pm25)

# Then, we fit an exponential semi-variogram to it via WLS
variog.pm25 <-
  fit.variogram(emp.variog.pm25, vgm(psill = 4.0, "Exp", 300, 2.0))

# Get the nugget
est.nugget <- variog.pm25$psill[1]
est.nugget

# Get the sill
est.psill <- variog.pm25$psill[2]
est.psill

# Get the range
est.range <- variog.pm25$range[2]
est.range

# Finally, we fit a spatial linear regression model to the daily average PM2.5 concentration.
# The mean model includes the explanatory variables: daily average temperature, 
# daily average rainfall amount, population density and elevation.
# The variance model is given by the exponential covariance function.
# The function we use to fit a spatial linear regression model via ML or REML is the function likfit. 

pm25.df <-
  data.frame(
    cbind(x.coord.km, y.coord.km, pm25, temperature, ppt, pop.dens, elev),
    nrow = length(x.coord.km),
    ncol = 7
  )

pm25.geo <-
  as.geodata(
    pm25.df,
    coords.col = c(1, 2),
    data.col = 3,
    covar.cols = c(4, 7)
  )

pm25.reml <-
  likfit(
    pm25.geo,
    trend = ~ temperature + ppt + pop.dens + elev,
    cov.model = "exponential",
    ini = c(est.psill, est.range),
    nugget = est.nugget,
    fix.nug = FALSE,
    lik.met = "REML"
  )


# Having fitted the spatial linear regression model to the data, 
# we can now perform kriging and estimate PM2.5 concentration at the selected location.
# In order to do this, we first have to transform the coordinate of the selected location into Easting and Northing in km.
# To transform latitude and longitude into Easting and Northing, we use the packages sp and rgdal.

## Make sure you have these packages
#install.packages(c("sp","rgdal"),repos="https://cloud.r-project.org")

library(sp)
library(rgdal)


SP.longlat.bw <-
  SpatialPoints(
    coords = cbind(as.numeric(bw.lon), as.numeric(bw.lat)),
    proj4string = CRS("+proj=longlat +datum=WGS84")
  )
# This transforms longitude and latitude in UTM coordinates
SP.utm.bw <-
  spTransform(SP.longlat.bw, CRS("+proj=utm +zone=10 +datum=WGS84"))

bw.x <- coordinates(SP.utm.bw)[, 1]
bw.y <- coordinates(SP.utm.bw)[, 2]

bw.x.km <- bw.x / 1000
bw.y.km <- bw.y / 1000

#### Ordinary Kriging ####
# The function that actually performs kriging is the function krige.conv in geoR.
# The function that specifies the type of kriging and provides all the necessary specification and information to implement kriging is the function krige.control.

# In the function krige.control, we need to specify:
# (i) type of kriging (type.krige): for both ordinary kriging and universal kriging, here we should type "ok";
# (ii) the mean model at the locations where we have the environmental exposure data (trend.d): for ordinary kriging, the mean model is a equal to an unknown, but constant value beta0. For ordinary kriging, we write "cte" in trend.d;
# (iii) the mean model at the location where we want to make predictions (trend.l): for ordinary kriging, this is again "cte";
# (iv) the variance model (cov.model): here we express the name of the semi-variogram model or covariance function model to fit to the data. Choices include "exponential" or "gaussian" (among others).
# (v) the parameters of the covariance, or semi-variogram, function (cov.pars). This is a 2-dimensional vector that should provide in order the value of the partial sill (o sigma^2) and of the range parameter (or phi).
# (vi) that value of the nugget effect (nugget). This should be either a variable or a number that provides the estimated nugget effect.

kc.ok.control <-
  krige.control(
    type.krige = "ok",
    trend.d = "cte",
    trend.l = "cte",
    cov.model = "exponential",
    cov.pars = c(pm25.reml$cov.pars[1], pm25.reml$cov.pars[2]),
    nugget = pm25.reml$nugget
  )

# Having provided all the input for the type of kriging, this is the command that actually performs kriging.
# First we place the locations where we want to generate estimates of the environmental exposure into an n by 2 matrix.

loc.kriging <-
  matrix(c(bw.x[index.loc], bw.y[index.loc]), nrow = 1, ncol = 2)

kc.ok.s0 <-
  krige.conv(pm25.geo, locations = loc.kriging, krige = kc.ok.control)

## krige.conv: model with constant mean
## krige.conv: Kriging performed using global neighbourhood
names(kc.ok.s0)

## THIS IS THE PREDICTION ##
# The predicted environmental exposure is stored into the part called "predict", 
# while the kriging variance is what is stored into the part called "krige.var".
kc.ok.s0$predict

kc.ok.s0$krige.var

# We can create a 95% CI for the estimated environmental exposure by adding and subtracting
# 1.96*square root of the kriging variance to then estimated environmental expsure.
ci.95.kc.ok <-
  c(
    kc.ok.s0$predict - 1.96 * sqrt(kc.ok.s0$krige.var),
    kc.ok.s0$predict + 1.96 * sqrt(kc.ok.s0$krige.var)
  )
ci.95.kc.ok

# If the lower bound of the confidence interval for the selected location is less than 0, we transform the lower bound of the interval to 0, the lowest possibile value for the PM2.5 concentration.
### Universal kriging
# To perform universal kriging, the krige.control function needs to change from the form it had under ordinary kriging.
# Specifically, we need to change the form of the mean model at the environmental exposure data locations and at the locations where we want to estimate the environmental exposure.
# Hence, now:
# trend.d needs to provide the form of the mean model at the locations where we have observed the data: this is specified by "~" followed by the name of the explanatory variables at the environmental exposure locations.
# trend.l needs to provide the form of the mean model where we want to estimate environmental exposure. The format is still "~" followd by the name of the explanatory variables at the location(s) where we want to predict the environmental exposure.
# The remainder of the specifications for the krige.control are exactly the same as under ordinary kriging.

#### UNIVERSAL KRIGING ####

kc.uk.control <-
  krige.control(
    type.krige = "ok",
    trend.d =  ~ temperature + ppt + pop.dens + elev,
    trend.l =  ~ bw.temperature[index.loc] + bw.ppt[index.loc] + bw.pop.dens[index.loc] +
      bw.elev[index.loc],
    cov.model = "exponential",
    cov.pars = c(pm25.reml$cov.pars[1], pm25.reml$cov.pars[2]),
    nugget = pm25.reml$nugget
  )

loc.kriging <-
  matrix(c(bw.x[index.loc], bw.y[index.loc]), nrow = 1, ncol = 2)
kc.uk.s0 <-
  krige.conv(pm25.geo, locations = loc.kriging, krige = kc.uk.control)
## krige.conv: model with mean defined by covariates provided by the user
## krige.conv: Kriging performed using global neighbourhood
names(kc.uk.s0)

# The predicted environmental exposure is stored into the part called "predict", while the kriging variance is what is stored into the part called "krige.var".
kc.uk.s0$predict
##     data
## 5.469714
kc.uk.s0$krige.var
## [1] 6.899578
# A 95% CI for the estimated environmental exposure by adding and subtracting 1.96*square root of the kriging variance to then estimated environmental expsure.
ci.95.kc.uk <-
  c(
    kc.uk.s0$predict - 1.96 * sqrt(kc.uk.s0$krige.var),
    kc.uk.s0$predict + 1.96 * sqrt(kc.uk.s0$krige.var)
  )
ci.95.kc.uk
##       data       data
##  0.3213726 10.6180557