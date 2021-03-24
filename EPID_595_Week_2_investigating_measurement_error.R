# Investigating measurement error
# Investigating measurement error in environmental epidemiology
# This code investigate the issue of measurement error that arises in environmental epidemiology when there is spatial misalignment between the environmental exposure data and the health outcome data.
# For this examination, we will use the daily average PM2.5 concentration in California in February 2019 and a simulated dataset reporting the birthweight (in kilograms) of children born in California in March 2019.
# The birthweight dataset reports each baby’s birthweight and its geographical coordinates - latitude and logitude but also Easting and Northing in kms.
# Our goal in this activity are:to estimate daily average PM2.5 concentration at the babies’ locations using your preferred spatial interpolation method.
# estimate the effect of exposure to daily average PM2.5 concentration on birthweight via linear regression.
# The birthweight data has been simulated under the assumption of a reduction in birthweight of 12grams (e.g. gamma =
#                                                                                                          -0.012) per micrograms per cubic meter of PM2.5 concentration (the unit of measure of PM2.5 concentration). Effects of these magnitude have been reported in part of the scientific literature.


#install.packages(c("fields","gstat","geoR"),repos="https://cloud.r-project.org")
## 
## The downloaded binary packages are in
##  /var/folders/9n/gn30q9z56gv9qbdhfqjyfqn80000gn/T//RtmpF0Ne3Z/downloaded_packages
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


# This is the fictitious dataset with the geographical coordinates of the babies for which we have birthweight data information. The dataset provides information on the babies' birthweight in kilogram, the values of the explanatory variables (daily average temperature, daily avg rainfall amount, population density and elevation) at the birthweight locations, and the geographical coordinates in two formats: latitude and longitude and Easting and Northing in kms.
bw.data <- read.csv("Data/Simulated_birthwt_Cali.csv",sep=",",header=T)
names(bw.data)
## [1] "Longitude"   "Latitude"    "x_UTM_km"    "y_UTM_km"    "temperature"
## [6] "ppt"         "pop.density" "elevation"   "sim.bw"
bw.lon <- bw.data$longitude
bw.lat <- bw.data$latitude
bw.x <- bw.data$x_UTM_km
bw.y <- bw.data$y_UTM_km
bw.temperature <- bw.data$temperature
bw.ppt <- bw.data$ppt
bw.pop.dens <- bw.data$pop.density
bw.elev <- bw.data$elevation
bw <- bw.data$sim.bw
# Here we address the issue of spatial misalignment between the two datasets by performing Universal Kriging and kriging daily average PM2.5 at the babies' birthweight locations.
# To perform universal kriging we first need to estimate the spatial covariance function parameters.
# For this goal, we follow previous code's examples, and we: (i) work with the Easting and Northing coordintes in kms, (ii) construct an empirical semi-variogram and (iii) fit an exponential or Gaussian variogram to it.
x.coord.km <- x.coord/1000
y.coord.km <- y.coord/1000
pm25.df <- data.frame(cbind(x.coord.km,y.coord.km,pm25,temperature,ppt,pop.dens,elev))
emp.variog.pm25 <- variogram(pm25~temperature+ppt+pop.dens+elev,locations=~x.coord.km+y.coord.km,data=pm25.df)
emp.variog.pm25

plot(emp.variog.pm25)


# Fitting an exponential semi-variogram to it via WLS and saving the estimates
variog.pm25 <- fit.variogram(emp.variog.pm25,vgm(psill=4.0,"Exp",300,2.0))
est.nugget <- variog.pm25$psill[1]
est.psill <- variog.pm25$psill[2]
est.range <- variog.pm25$range[2]

# For the univeral kriging of PM2.5, we proceed as follows using an exponential covariance 
# function with parameter equal to the WLS estimates of the covariance function parameters 
# obtained above.
# First we create the geodata object with all the data.
pm25.df <- data.frame(cbind(x.coord.km,y.coord.km,pm25,temperature,ppt,pop.dens,elev),nrow=length(x.coord.km),ncol=7)
pm25.geo <- as.geodata(pm25.df,coords.col=c(1,2),data.col=3,covar.cols=c(4,7))
# Then we specify all the options regarding the type of kriging, the mean model, 
# where we want to predict the environmental exposure, the covariance function to use (e.g. the variance 
# model), and the estimates of parameters of the covariance function.
kc.uk.control <- krige.control(type.krige="ok",trend.d=~temperature+ppt+pop.dens+elev,trend.l=~bw.temperature+bw.ppt+bw.pop.dens+bw.elev,cov.model="exponential",cov.pars=c(est.psill,est.range),nugget=est.nugget)

loc.kriging <- matrix(cbind(bw.x,bw.y),nrow=length(bw.x),ncol=2)
kc.uk.bw.locs <- krige.conv(pm25.geo,locations=loc.kriging,krige=kc.uk.control)
## krige.conv: model with mean defined by covariates provided by the user
## krige.conv: Kriging performed using global neighbourhood
names(kc.uk.bw.locs)
## [1] "predict"      "krige.var"    "beta.est"     "distribution" "message"     
## [6] "call"
# The predicted environmental exposure is stored into the part called "predict".
pm25.bw.locs <- kc.uk.bw.locs$predict


# Now that we have the estimated environmental exposure, e.g. the daily average PM2.5 concentration, 
# we can estimate the health effect of PM2.5 daily average concentration on birthweight by 
# doing the following:

health.effect <- lm(bw~pm25.bw.locs)
summary(health.effect)
## 
## Call:
## lm(formula = bw ~ pm25.bw.locs)
## 
## Residuals:
##       Min        1Q    Median        3Q       Max 
## -0.288633 -0.066588  0.003844  0.075108  0.230157 
## 
## Coefficients:
##               Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   3.091759   0.023805 129.881   <2e-16 ***
## pm25.bw.locs -0.001824   0.004235  -0.431    0.667    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.09514 on 248 degrees of freedom
## Multiple R-squared:  0.0007479,  Adjusted R-squared:  -0.003281 
## F-statistic: 0.1856 on 1 and 248 DF,  p-value: 0.667
# What did you obtain as estimate of the health effects?
# How much do the results change in terms of spatial interpolation if we use the REML estimates of 
# the covariance parameters instead of the WLS estimates as done here? Look at the interpolated values 
# of daily average PM2.5 concentration at the first 3 birthweight locations that you obtain in the 
# two cases.
# How do the results change in terms of health effects if we change the form of the spatial linear 
# regresison model used, specifically the mean model (the explanatory variables) in Universal Kriging?
# How do the results change in terms of health effects if we change the form of the spatial linear 
# regression model, particularly the variance model (e.g. change the form of the covariance function) 
# in Universal Kriging?