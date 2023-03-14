#########################################################
#### EPID 595 Cod set 2 #################################
#### Fitting semi-variograms

# Fitting semi - variogram models to empirical semi - variograms via Weighted Least Squares in R
# This code shows how to fit a semi - variogram model to an empirical semi -
#   variogram constructed from point referenced spatial data using the method of WLS in R.
# The R package that is needed to fit a semi - variogram model to an empirical semi -
#   variogram is contained in the gstat package. You can download the package from a CRAN repository.


#install.packages("gstat",repos="https://cloud.r-project.org")
## 
## The downloaded binary packages are in
##  /var/folders/9n/gn30q9z56gv9qbdhfqjyfqn80000gn/T//RtmpKaLeoC/downloaded_packages
library(gstat)

# Reading in the daily average PM2.5 concentration dataset presented in Unit 1
pm25.data <- read.table("Avg_PM25_Feb2019_California_and_other_vars.csv",sep=",",header=T)
# To look at the names of the variables in the PM2.5 dataset, we can simply do the following:
names(pm25.data)
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


# Here we construct the empirical semi-variogram for the residuals of the multiple linear regression model regressing daily average PM2.5 concentration on the covariates, i.e. daily average temperature, rainfall amount, population density and elevation.
pm25.df <- data.frame(cbind(x.coord,y.coord,pm25,temperature,ppt,pop.dens,elev))
emp.variog.pm25 <- variogram(pm25~temperature+ppt+pop.dens+elev,locations=~x.coord+y.coord,data=pm25.df)
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
# The emp.variog.pm25 object contains, among others, a vector (dist) reporting the midpoints of the distance bins used to construct the empirical semi-variogram, and a vector (gamma) reporting the values of the semi-variances at the bin midpoints according to the empirical semi-variogram 


# To fit a semi-variogram model to an empirical semi-variogram via the method of Weighted Least Squares, we use the function fit.variogram.
# Note that since Easting and Northing are given in terms of meters (and not kms), if we want to use as initial value for the range parameter, phi, 300 km, we have to specify as initial value for phi 300,000 meters. We call this Option 1.
# Alternatively (Option 2), we can divide the Easting and Northing coordinates by 1000, thus obtaining the Easting and Northing coordinates in km, recalculate the empirical semi-variogram using Easting and Northing in kms rather than in meters and then fit a semi-variogram model to the empirical semi-variogram providing the initial value for phi on the km scale.

# Fitting a semi-variogram using Option 1: We specify as initial values for sigma^2, tau^2 and phi the following values (as discussed in Unit 1): 6 for the partial sill (psill), 2 for the nugget effect (nugget) and 300km or 300,000km for the range parameter (range). We choose an exponential semi-variogram model. We provide all these choices in the option vgm within the fit.variogram function, specifying in order, first, the value of the partial sill, then the semi-variogram model chosen (either "Exp" for exponential or "Gau" for Gaussian), the range parameter on a km scale, and finally the nugget effect.
variog.pm25 <- fit.variogram(emp.variog.pm25,vgm(psill=4.0,"Exp",300000,2.0))
variog.pm25
##   model    psill    range
## 1   Nug 1.048898     0.00
## 2   Exp 4.944573 66257.04
# The estimates of the semi-variogram parameters are reported in a 2x2 matrix, where the top left number provides the estimated nugget effect, the bottom left number reports the estimated partial sill, while the bottom right number is the estimated range parameter in meters.
# To plot the fitted semi-variogram on top of the empirical semi-variogram, we use the function variogramLine.
# The function variogramLine reports the value of a semi-variogram model at a given distance when the semi-variogram parameters are equal to provided values.
# Setting the semi-variogram parameters equal to their WLS values we have:
est.nugget <- variog.pm25$psill[1]
est.nugget
## [1] 1.048898
est.psill <- variog.pm25$psill[2]
est.psill
## [1] 4.944573
est.range <- variog.pm25$range[2]
est.range
## [1] 66257.04
# The fitted exponential semi-variogram is obtained as follows.
# The syntax of the variogramLine, asks first to write "vgm()" and then provide, in order, the value of the estimated partial sill, the name of the semi-variogram model (here "Exp" for exponential), the estimated range parameters, and finally the estimated nugget effect. After this, the variogramLine function requests the user to provide a vector, called  dist_vector, with the distances where R should evaluate the semi-variogram model. Here we use the midpoints of the bins used to construct the empirical semi-variogram, which are stored in emp.variog.pm25$dist.
fitted.variog <- variogramLine(vgm(est.psill,"Exp",est.range,est.nugget),dist_vector=emp.variog.pm25$dist)

# Finally to plot both the empirical semi-variogram and the fitted exponential semi-variogram we proceed as follows
plot(emp.variog.pm25$dist,emp.variog.pm25$gamma,type="p",pch=19,xlab="Distance (in meters)",ylab="Semi-variance")
lines(emp.variog.pm25$dist,fitted.variog$gamma,col="red")

# In Option 2, we repeat the same steps as above, the only difference being that we divide the Easting and Northing coordinates by 1000, so that they are now provided in km's rather than in meters.
x.coord.km <- x.coord/1000
y.coord.km <- y.coord/1000

# Then we repeat all the steps above, making sure to replace x.coord and y.coord with x.coord.km and y.coord.km whenever we used them in the empirical semi-variogram construction.
pm25.df.km <- data.frame(cbind(x.coord.km,y.coord.km,pm25,temperature,ppt,pop.dens,elev))
emp.variog.pm25.km <- variogram(pm25~temperature+ppt+pop.dens+elev,locations=~x.coord.km+y.coord.km,data=pm25.df.km)
emp.variog.pm25.km
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
variog.pm25.km <- fit.variogram(emp.variog.pm25.km,vgm(psill=4.0,"Exp",300,2.0))
variog.pm25.km
##   model    psill    range
## 1   Nug 1.048932  0.00000
## 2   Exp 4.944584 66.25919
fitted.variog.km <- variogramLine(vgm(est.psill,"Exp",est.range,est.nugget),dist_vector=emp.variog.pm25.km$dist)

# Finally to plot both the empirical semi-variogram and the fitted exponential semi-variogram we proceed as follows
plot(emp.variog.pm25.km$dist,emp.variog.pm25.km$gamma,type="p",pch=19,xlab="Distance (in meters)",ylab="Semi-variance")
lines(emp.variog.pm25.km$dist,fitted.variog$gamma,col="red")