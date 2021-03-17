############################################
#### EPID 595 Code set 3 ###################
#### Fitting spatial linear regressions in R

# Fitting a spatial linear regression to point referenced data in R
# This code shows how to fit a spatial linear regression model to point referenced data in R 
# using the method of Maximum Likelihood Estimation or the method of 
# Restricted Maximum Likelihood Estimation (REML).
# The R package that is needed to obtain REML or MLE estimates of parameters in a spatial 
# linear regression is the package geoR. You can download the package from a CRAN repository.
# If your computer is a Mac, you will need to install one or two potential other software:XQuartz and Xcode. 
# The latter takes a little bit of time to install.
# 


install.packages("gstat", repos = "https://cloud.r-project.org")
## 
## The downloaded binary packages are in
##  /var/folders/9n/gn30q9z56gv9qbdhfqjyfqn80000gn/T//RtmpHTfyXs/downloaded_packages
library(gstat)
install.packages("geoR",repos="https://cloud.r-project.org")
## 
## The downloaded binary packages are in
##  /var/folders/9n/gn30q9z56gv9qbdhfqjyfqn80000gn/T//RtmpHTfyXs/downloaded_packages
library(geoR)
## --------------------------------------------------------------
##  Analysis of Geostatistical Data
##  For an Introduction to geoR go to http://www.leg.ufpr.br/geoR
##  geoR version 1.8-1 (built on 2020-02-08) is now loaded
## --------------------------------------------------------------
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
# Before doing so, we transform Easting and Northing from  meters to kilometers.
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

# Here we fit a semi-variogram model to an empirical semi-variogram via the method of Weighted Least Squares, we use the function fit.variogram. 
# We specify as initial values for sigma^2, tau^2 and phi the following values (as discussed in Unit 1): 6 for the partial sill (psill), 2 for the nugget effect (nugget) and 300km for the range parameter (range). We choose an exponential semi-variogram model. 
variog.pm25 <- fit.variogram(emp.variog.pm25,vgm(psill=4.0,"Exp",300,2.0))
variog.pm25
##   model    psill    range
## 1   Nug 1.048932  0.00000
## 2   Exp 4.944584 66.25919
# Setting the semi-variogram parameters equal to their WLS values we have:
est.nugget <- variog.pm25$psill[1]
est.nugget
## [1] 1.048932
est.psill <- variog.pm25$psill[2]
est.psill
## [1] 4.944584
est.range <- variog.pm25$range[2]
est.range
## [1] 66.25919
# Here plot the empirical semi-variogram and the fitted exponential semi-variogram
fitted.variog <- variogramLine(vgm(est.psill,"Exp",est.range,est.nugget),dist_vector=emp.variog.pm25$dist)
plot(emp.variog.pm25$dist,emp.variog.pm25$gamma,type="p",pch=19,xlab="Distance (in meters)",ylab="Semi-variance")
lines(emp.variog.pm25$dist,fitted.variog$gamma,col="red")

# Here we fit a spatial linear regression model to the daily average PM2.5 concentration.
# The mean model includes the explanatory variables: daily average temperature, daily average rainfall amount, population density and elevation.
# The variance model is given by the exponential covariance function.
# The function we use to fit a spatial linear regression model via ML or REML is the function likfit.
# In order to use the function likfit, the data has to be put in a geodata format. A geodata format is a dataframe where it is clearly specified what columns in the dataframe refer to the geographical coordinates, what column contains the response variable and what columns contain the covariates.
# Here we create the dataframe with all the data.
pm25.df <- data.frame(cbind(x.coord.km,y.coord.km,pm25,temperature,ppt,pop.dens,elev),nrow=length(x.coord.km),ncol=7)
pm25.geo <- as.geodata(pm25.df,coords.col=c(1,2),data.col=3,covar.cols=c(4,7))

# Having created the geodata, to fit a spatial linear regression model, we use the function likfit.
# The function requires the following input:
# (i) the name of the geodata object containing the data
# (ii) trend: this specifies the mean model. The format for trend is: "trend=~", to the RHS of the "~" the user needs to list the covariates that are included in the mean model.
# (iii) cov.model: this specifies the variance model. Specifically, the user has to indicate what parametric semi-variogram (equivalently, covariance function) intends to use. Choices are: "exponential", "gaussian".
#(iv) ini: this specifies the initial values for the parameters of the parametric semi-variogram. The order in which these initial values must be provided is: initial value for the partial sill (or sigma^2) and inital value for the range parameter (or phi).
# (v) nugget: this specifies the initial value for the nugget effect.
# (vi) fix.nugget: this can be either TRUE or FALSE. It asks whether the nugget should be estimated or kept fixed to the initial value provided by the user.
# (vii) lik.met: this specifies what type of estimation R should perform. Choices are "ML" for Maximum Likelihood or "REML" for Restricted Maximum Likelihood.
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
# This provides just the summary output of the likfit command which provides the REML estimates of the regression coefficients and the covariance function parameters
pm25.reml
## likfit: estimated model parameters:
##     beta0     beta1     beta2     beta3     beta4     tausq   sigmasq       phi 
## " 6.4499" "-0.0387" "-0.0107" " 0.0001" " 0.0002" " 0.4250" " 5.9767" "36.8937" 
## Practical Range with cor=0.05 for asymptotic range: 110.5237
## 
## likfit: maximised log-likelihood = -233.5
# Using the function summary on the output of the likfit function provides a longer view of the output generated by the likfit function.
summary(pm25.reml)
## Summary of the parameter estimation
## -----------------------------------
## Estimation method: restricted maximum likelihood 
## 
## Parameters of the mean component (trend):
##   beta0   beta1   beta2   beta3   beta4 
##  6.4499 -0.0387 -0.0107  0.0001  0.0002 
## 
## Parameters of the spatial component:
##    correlation function: exponential
##       (estimated) variance parameter sigmasq (partial sill) =  5.977
##       (estimated) cor. fct. parameter phi (range parameter)  =  36.89
##    anisotropy parameters:
##       (fixed) anisotropy angle = 0  ( 0 degrees )
##       (fixed) anisotropy ratio = 1
## 
## Parameter of the error component:
##       (estimated) nugget =  0.425
## 
## Transformation parameter:
##       (fixed) Box-Cox parameter = 1 (no transformation)
## 
## Practical Range with cor=0.05 for asymptotic range: 110.5237
## 
## Maximised Likelihood:
##    log.L n.params      AIC      BIC 
## "-233.5"      "8"    "483"  "504.9" 
## 
## non spatial model:
##    log.L n.params      AIC      BIC 
## "-250.8"      "6"  "513.7"  "530.1" 
## 
## Call:
## likfit(geodata = pm25.geo, trend = ~temperature + ppt + pop.dens + 
##     elev, ini.cov.pars = c(est.psill, est.range), fix.nugget = FALSE, 
##     nugget = est.nugget, cov.model = "exponential", lik.method = "REML")
# To see the full list of information that the likfit function returns, we can apply the command names onto the oject, output of the likfit function.
names(pm25.reml)
##  [1] "cov.model"                  "nugget"                    
##  [3] "cov.pars"                   "sigmasq"                   
##  [5] "phi"                        "kappa"                     
##  [7] "beta"                       "beta.var"                  
##  [9] "lambda"                     "aniso.pars"                
## [11] "tausq"                      "practicalRange"            
## [13] "method.lik"                 "trend"                     
## [15] "loglik"                     "npars"                     
## [17] "AIC"                        "BIC"                       
## [19] "parameters.summary"         "info.minimisation.function"
## [21] "max.dist"                   "trend"                     
## [23] "trend.matrix"               "transform.info"            
## [25] "nospatial"                  "model.components"          
## [27] "call"
# For example, one could also obtain the AIC value corresponding to a spatial linear regresison model.
pm25.reml$AIC
## [1] 482.9675
# We could also create 95% confidence interval for the regression coefficients.
# For this, we use two parts of the output of the likfit function: the beta part which contains the estimates of the regression coefficients; and the betea.var part which provides a matrix with the estimated variance of the regression coefficients in the diagonal part of the matrix, or the squared standard errors of the regression coefficients.
# We extract the diagonal part of the matrix by doing the following:
squared.se.beta <- diag(pm25.reml$beta.var)
# We simply take the standard errors of the beta regression coefficients (by taking the square root of the squared se's),multiply them by the appropriate multiplier (e.g. 1.65 for a 90% CI and 1.96 for a 95% CI) and add and subtract them to the estimated regression coefficients.
# Example: 95% CI for the regression coefficients:
ci.beta <- cbind(pm25.reml$beta-1.96*sqrt(squared.se.beta),pm25.reml$beta+1.96*sqrt(squared.se.beta))
ci.beta
##                    [,1]         [,2]
## intercept -8.358363e+00 2.125825e+01
## covar1    -3.287034e-01 2.512392e-01
## covar2    -1.379820e-01 1.165964e-01
## covar3    -4.686407e-05 3.418452e-04
## covar4    -8.521013e-04 1.228524e-03