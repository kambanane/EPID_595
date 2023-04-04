# Fitting a spatial regression to areal data using a CAR model
# Fitting a spatial regression to areal data using a CAR model in R
# This code shows how to fit a spatial regression model to areal data using spatial random
# effects modeled using a Conditionally AutoRegressive (CAR) model in R.
# The R package that is needed for the analysis of areal data is the CARBayes package. 
# In addition to the CARBayes package, the same packages that are needed to creating adjacency weights and to plot spatial areal data are needed (e.g. maps, maptools, spdep, RColorBrewer, and classInt). All these packages can be downloaded from a CRAN repository.
# knitr::opts_chunk$set(fig.width=12, fig.height=8) 
# Installing the necessary packages

sf::sf_use_s2(FALSE)
# 
# install.packages(c(
#   "maps",
#   "maptools",
#   "spdep",
#   "RColorBrewer",
#   "classInt",
#   "CARBayes"
# ),
# repos = "https://cloud.r-project.org")
## 
## The downloaded binary packages are in
##  /var/folders/9n/gn30q9z56gv9qbdhfqjyfqn80000gn/T//RtmpeDTZeO/downloaded_packages
library(maps)
library(maptools)
library(spdep)
library(RColorBrewer)
library(classInt)
library(CARBayes)

# Reading in the data on average violent crime arrest rate per 100,000 people between the ages of 10-17 years in Michigan 
# counties during years 2011-2012, and exposure to high blood lead levels in children during years 2005-2006.
mich.data <- read.table("Data/mich_crime_lead.csv",
                        sep = ",",
                        header = T)
# To look at the names of the variables in the Michigan crime and lead dataset, we can simply do the following:
names(mich.data)
##  [1] "county"                 "ebll"                   "crime_rate"            
##  [4] "p_foodstamp"            "p_male_juve"            "p_whiteundr18"         
##  [7] "p_blackundr18"          "p_aianundr18"           "p_asianundr18"         
## [10] "p_hispanicundr18"       "pct_insur"              "belowpovertylastyr_pct"
## [13] "medianage"              "unemploy_rate"
ebll <- mich.data$ebll
unempt <- mich.data$unemploy_rate
poverty <- mich.data$belowpovertylastyr_pct

# The outcome variable is the average violent crime arrest rate in Michigan counties.
# Here we create a histogram of the response variable to see if it is distributed according to an approximate normal distribution.
hist(
  mich.data$crime_rate,
  breaks = 30,
  xlab = "Average violent crime arrest rate per 100,000 people",
  col = "dodgerblue",
  main = "Histogram of avg. violent crime arrest rate"
)

# Given the long tail to the right, we decide to transform the variable and work instead with 
# log(average rate). This transformation will address: the right skewness in the distribution and the fact that the crime rate
# can only be positive; working on the log, we can also have negative log violent crime arrest rate.
# Since some counties have no violent crime arrests, and log of 0 is -infinity, to not have this issue, we replace 
# any observation of 0 with a 1. This is not a major problem since, when we take the logarithm, log of 1 is 0.
y.rev <- mich.data$crime_rate
y.rev[which(mich.data$crime_rate == 0)] <-
  rep(1, length(which(mich.data$crime_rate == 0)))

# Here we take the logarithm of the variable where we replaced the 0's with 1's
l.yrev <- log(y.rev)
# Here we create the histogram of the new variable. As the histogram shows, the new variable has a distribution that 
# is more closely Gaussian.
hist(l.yrev,breaks=40,xlab="Average log violent crime arrest rate per 100,000 people",col="orchid3",main="Histogram of log avg. violent crime arrest rate")

# This part of the code create chloropleth maps of the response variables and of the explanatory variables
mich.county <- map("county", "michigan", fill = T, plot = F)
mich.IDs <- sapply(strsplit(mich.county$names, ":"), function(x)
  x[1])
mich.poly <-
  map2SpatialPolygons(mich.county,
                      IDs = mich.IDs,
                      proj4string = CRS("+proj=longlat +datum=WGS84"))

# Plotting the response variable: log average violent crime rates per 100,000 people between 2011-2012
plotvar <- l.yrev
nclr <- 5
plotclr <- brewer.pal(nclr,"Purples")
class <- classIntervals(plotvar,nclr,style="jenks")
colcode <- findColours(class,plotclr)

plot(mich.poly,
     border = "black",
     axes = TRUE,
     xlim = c(-92, -81))
title(xlab = "Longitude", ylab = "Latitude", main = "Log average violent crime arrest rate, years 2011-2012")
plot(mich.poly, col = colcode, add = T)

leg.txt<-c("0",paste("(",round(class$brks[2],3),"; ",round(class$brks[3],3),")",sep=""),paste("[",round(class$brks[3],3),"; ",round(class$brks[4],3),")",sep=""),paste("[",round(class$brks[4],3),"; ",round(class$brks[5],3),")",sep=""),paste("[",round(class$brks[5],3),"; ",round(class$brks[6],3),")",sep=""))
legend(-92,44,legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")

# Plotting the exposure: proportion of children that in 2005-2006 reported elevated blood lead level
plotvar <- ebll
nclr <- 5
plotclr <- brewer.pal(nclr,"Greens")
class <- classIntervals(plotvar,nclr,style="jenks")
colcode <- findColours(class,plotclr)

plot(mich.poly,border="black",axes=TRUE,xlim=c(-92,-81))
title(xlab="Longitude",ylab="Latitude",main="Proportion of children with elevated blood lead level in 2005-2006")
plot(mich.poly,col=colcode,add=T)

leg.txt<-c(paste("[",round(class$brks[1],3),"; ",round(class$brks[2],3),")",sep=""),paste("[",round(class$brks[2],3),"; ",round(class$brks[3],3),")",sep=""),paste("[",round(class$brks[3],3),"; ",round(class$brks[4],3),")",sep=""),paste("[",round(class$brks[4],3),"; ",round(class$brks[5],3),")",sep=""),paste("[",round(class$brks[5],3),"; ",round(class$brks[6],3),")",sep=""))
legend(-92,44,legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")

# Plotting unemployment rate
plotvar <- unempt
nclr <- 5
plotclr <- brewer.pal(nclr,"Reds")
class <- classIntervals(plotvar,nclr,style="jenks")
colcode <- findColours(class,plotclr)

plot(mich.poly,border="black",axes=TRUE,xlim=c(-92,-81))
title(xlab="Longitude",ylab="Latitude",main="Unemployment rate in 2011")
plot(mich.poly,col=colcode,add=T)

leg.txt<-c(paste("[",round(class$brks[1],3),"; ",round(class$brks[2],3),")",sep=""),paste("[",round(class$brks[2],3),"; ",round(class$brks[3],3),")",sep=""),paste("[",round(class$brks[3],3),"; ",round(class$brks[4],3),")",sep=""),paste("[",round(class$brks[4],3),"; ",round(class$brks[5],3),")",sep=""),paste("[",round(class$brks[5],3),"; ",round(class$brks[6],3),")",sep=""))
legend(-92,44,legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")

# Plotting poverty
plotvar <- poverty
nclr <- 5
plotclr <- brewer.pal(nclr,"Greys")
class <- classIntervals(plotvar,nclr,style="jenks")
colcode <- findColours(class,plotclr)

plot(mich.poly,border="black",axes=TRUE,xlim=c(-92,-81))
title(xlab="Longitude",ylab="Latitude",main="Proportion of household living under poverty line")
plot(mich.poly,col=colcode,add=T)

leg.txt<-c(paste("[",round(class$brks[1],3),"; ",round(class$brks[2],3),")",sep=""),paste("[",round(class$brks[2],3),"; ",round(class$brks[3],3),")",sep=""),paste("[",round(class$brks[3],3),"; ",round(class$brks[4],3),")",sep=""),paste("[",round(class$brks[4],3),"; ",round(class$brks[5],3),")",sep=""),paste("[",round(class$brks[5],3),"; ",round(class$brks[6],3),")",sep=""))
legend(-92,44,legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")

# Here we perform some exploratory analysis: computing Moran's I on the response variable directly to see if the response variable
# shows spatial correlation, and computing Moran's I on the residuals of a linear regression regressing log average violent crime arrest 
# rate on the explanatory variables.
mi.nb <- poly2nb(mich.poly)
mi.nb
## Neighbour list object:
## Number of regions: 83 
## Number of nonzero links: 420 
## Percentage nonzero weights: 6.096676 
## Average number of links: 5.060241
mi.listw <- nb2listw(mi.nb,style="B",zero.policy=TRUE)
moran.i.crime <- moran.test(l.yrev, mi.listw, rank=FALSE,zero.policy=TRUE)
moran.i.crime
## 
##  Moran I test under randomisation
## 
## data:  l.yrev  
## weights: mi.listw    
## 
## Moran I statistic standard deviate = -0.71477, p-value = 0.7626
## alternative hypothesis: greater
## sample estimates:
## Moran I statistic       Expectation          Variance 
##      -0.059301822      -0.012195122       0.004343463
# Based on Moran's I, the log average violent crime arrest rate does not present significant amount of spatial dependence.
# We consider the residuals of the linear regression of the log average violent crime arrest rate on unemployment rate, proportion of 
# households living under the poverty line and the proportion of children that in 2005-2006 had elevated blood lead levels.
lm.crime <- lm(l.yrev~ebll+unempt+poverty)
res.lm.crime <- as.numeric(lm.crime$residuals)
moran.i.res.crime <- moran.test(res.lm.crime, mi.listw, rank=FALSE,zero.policy=TRUE)
moran.i.res.crime
## 
##  Moran I test under randomisation
## 
## data:  res.lm.crime  
## weights: mi.listw    
## 
## Moran I statistic standard deviate = -0.54129, p-value = 0.7058
## alternative hypothesis: greater
## sample estimates:
## Moran I statistic       Expectation          Variance 
##      -0.047831587      -0.012195122       0.004334415
# Also the residuals of the linear regression do not present significant spatial correlation, based on the value of Moran's I.

# Despite the non-significant Moran's I on the residuals of the linear regression, we fit a spatial linear regression model
# on the log average violent crime arrest rate with spatial random effects that are modeled using a CAR model.
# To fit such a model we use the function S.CARleroux. This function requires that the user provides R a matrix W with as many rows and as many 
# columns as the number of areal unit, where the (i,j)-th entry is the binary adjacency weight for areal unit i and areal unit j.
# Here we create the adjacency matrix W of dimension N by N where N is the number of areal units, or number of observations.
N <- length(l.yrev)
# To create the adjaceny matrix, we first need to know which areal unit is adjacent to which and how many neighbors an areal unit has.
mi.weights <- nb2WB(mi.nb)
# This is the vector that lists the IDs of all the areal units neighbor to areal unit 1, 2,.... and N.
adj.mi <- mi.weights$adj
# This is the vector that indicates how many neighbors an areal unit has.
num.mi <- mi.weights$num

# The following command creates a vector where the unit 1 is repeated as many times as its number of neighbors, and so forth 
# for unit 2, 3, .., N.
rep.mi <- rep(1:N,num.mi)
W <- matrix(0,N,N)
for(i in 1:N){
  # This guarantees that there is a 1 on the i-th row at all the columns j, that are neighor of areal unit i, for i=1,2,...,N
  W[i,adj.mi[rep.mi==i]] <- rep(1,num.mi[i])
}

# Having created the adjacency matrix W, we can pass onto running the function S.CARleroux to fit the spatial linear regression model 
# to the data. This first requires the following information:
# 1. "formula": this requires the formula of the model that explains how the mean of the response variable depends on the explanatory 
# variables.
# 2. "W": this requires the adjacency matrix W
# 3. "family": this requires the name of the distribution of the data. Possibilities include: "binomial",
# "gaussian","poisson", "zip".
# 4. "burnin": this requires the user to specify the number of iterations to run for burn-in.
# 5: "n.sample": this asks for the number of MCMC iteration to run in total.
# 6: "thin": this asks for a number that indicates whether the MCMC output should be used every iteration (if "thin"=1) 
# or used every certain number k of iteration (if "thin"=k).
# 7: "prior.mean.beta": this option asks for the prior mean of the regression coefficients. Specifying "NULL" tells
# R to use the default, which is a vector of all zeroes.
# 8: "prior.var.beta": this option asks for a vector with the variances for all the regresison coefficients. Specifying
# "NULL" tells R to use the default, which is vector where all the variances are equal to 100,000.
# 9: "prior.nu2": this option asls for a vector with the shape and scale for the Inverse Gamma distribution used as 
# prior for the non-spatial variance. Specifying "NULL" tells R to use the default, which is: shape equal to 1
# and scale equal to 0.01.
# 10" "prior.tau2": as in 9, this asks for the vector with the shape and scale for the Inverse Gamma distribution 
# used as prior for the spatial variance. Specifying "NULL" tells R to use the default, which is: shape equal to 1
# and scale equal to 0.01.
# 11: "rho": this should always be set equal to 1, if the spatial random effects are modeled using
# the Conditionally AutoRegressive (CAR) model.
# 12: "verbose": this options asks the user whether they want to be updated on the progress. Setting equal to "TRUE"
# provides the user with updates on the MCMC algorithm.
# 
# Here, we use all the defaults on the prior distributions. We ran the algorithm for 100,000 iterations, of which 
# half are thrown away for burn-in. We use every iteration (thin=1). We chose 100,000 iterations for the MCMC algorithm
# based on the convergence diagnostics results, which did not indicate convergence until more iterations were considered.
formula <- l.yrev~ebll+unempt+poverty
model.car <- S.CARleroux(
  formula = formula,
  W = W,
  family = "gaussian",
  burnin = 50000,
  n.sample = 100000,
  thin = 1,
  prior.mean.beta = NULL,
  prior.var.beta = NULL,
  prior.nu2 = NULL,
  prior.tau2 = NULL,
  rho = 1,
  verbose = TRUE
)

## Summarising results.
## Finished in  28.3 seconds.
# By typing the name of the S.CARleroux followed by $summary.results, we obtain the summary statistics on the 
# MCMC sample post burn-in. The summary statistics indicates the median of the samples post-burnin, the 
# lower and upper boundary of a 95% credible intervals. How many samples were used post-burnin, the percentage 
# of samples accepted (here all, since the MCMC algorithm does not require the Metropolis-Hasting algorithm).
# The last column reports the value of an MCMC convergence diagnostics, Geweked diagnostics. This test uses a 
# z-test, thus the values reported are z-scores. If any of the values reported in this column is above 2.0,
# it means that the parameter has not converged, thus the number of iterations has to be increased for both
# n.sample and burnin.
model.car$summary.results
##              Median     2.5%   97.5% n.sample % accept n.effective Geweke.diag
## (Intercept)  3.8419   2.0007  5.7027    50000      100      1685.7         2.1
## ebll        22.3677 -28.5755 72.8003    50000      100      4984.5        -2.1
## unempt      -0.1020  -0.2921  0.0876    50000      100      1475.9        -1.9
## poverty      7.8488  -2.7316 18.4969    50000      100      5101.9         1.9
## nu2          2.7613   2.0311  3.8355    50000      100     13192.8        -2.3
## tau2         0.0143   0.0027  0.3540    50000      100       516.4         1.5
## rho          1.0000   1.0000  1.0000       NA       NA          NA          NA
# If we want to have access to the MCMC samples for spatial random effects, we look at the object samples and at 
# the element called "phi" within that object.
# The samples are organized in a matrix with as many rows as MCMC samples post-burnin and as many columns as areal
# units.
samples.eta <- model.car$samples$phi
# This comand calculates the median of the MCMC samples for the spatial random effects by taking the median 
# of each column.
post.median.eta <- as.numeric(apply(samples.eta,2,median))

# Here we create a chloropleth map of the posterior median of the spatial random effects
plotvar <- post.median.eta
nclr <- 5
# We use a palette with 5 shades that goes from red to blue where we invert the colors so that the 2 shades of 
# blue correspond to lower (hopefully) negative values ad the 2 shades of red correspond to larger (hopefully) 
# positive values.
plotclr <- brewer.pal(nclr,"RdBu")[5:1]
class <- classIntervals(plotvar,nclr,style="jenks")
colcode <- findColours(class,plotclr)

plot(mich.poly,border="black",axes=TRUE,xlim=c(-90,-82))
title(xlab="Longitude",ylab="Latitude",main="Spatial random effects \n CAR model")
plot(mich.poly,col=colcode,add=T)

leg.txt<-c("[-1.35;-1.0)","[-1.0;-0.5)","[-0.5;0.0)","[0.0;0.01)","[0.01;0.025]")
legend(-90.3,44,legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")

# Here we calculate a 95% credible interval for the spatial random effects for each areal unit.
# We start by calculating the lower bound of the 95% credible interval. For this, for each areal unit we compute 
# the 2.5-th percentile of the MCMC samples of the spatial random effects for the given areal unit.
# Hence, we calculate the 2.5-th percentile of each column.
low.bd95ci.eta <- as.numeric(apply(samples.eta,2,quantile,0.025))
summary(low.bd95ci.eta)
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## -2.6146 -0.3380 -0.2594 -0.6359 -0.2304 -0.1827
# Here we create a chloropleth map of the lower bound of the 95% credible interval for the spatial random 
# effects.
plotvar <- low.bd95ci.eta
nclr <- 5

plotclr <- brewer.pal(nclr,"Blues")[5:1]
class <- classIntervals(plotvar,nclr,style="jenks")
colcode <- findColours(class,plotclr)

plot(mich.poly,border="black",axes=TRUE,xlim=c(-90,-82))
title(xlab="Longitude",ylab="Latitude",main="Lower bound of 95% CI for spatial random effects \n CAR model")
plot(mich.poly,col=colcode,add=T)

leg.txt<-c("0",paste("(",round(class$brks[2],3),"; ",round(class$brks[3],3),")",sep=""),paste("[",round(class$brks[3],3),"; ",round(class$brks[4],3),")",sep=""),paste("[",round(class$brks[4],3),"; ",round(class$brks[5],3),")",sep=""),paste("[",round(class$brks[5],3),"; ",round(class$brks[6],3),")",sep=""))
legend(-90.3,44,legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")

# Similarly, here we calculate the upper bound of the 95% credible interval. For this, for each areal unit we compute 
# the 97.5-th percentile of the MCMC samples of the spatial random effects for the given areal unit.
# Hence, we calculate the 97.5-th percentile of each column.
upp.bd95ci.eta <- as.numeric(apply(samples.eta,2,quantile,0.975))
summary(upp.bd95ci.eta)
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## -0.1048  0.1676  0.2323  0.2015  0.2854  0.4664
# Here we create a chloropleth map of the upper bound of the 95% credible interval for the spatial random 
# effects of each areal unit.
plotvar <- upp.bd95ci.eta
nclr <- 5


plotclr <- brewer.pal(nclr,"Reds")
class <- classIntervals(plotvar,nclr,style="jenks")
colcode <- findColours(class,plotclr)

plot(mich.poly,border="black",axes=TRUE,xlim=c(-90,-82))
title(xlab="Longitude",ylab="Latitude",main="Upper bound of 95% CI for spatial random effects \n CAR model")
plot(mich.poly,col=colcode,add=T)

leg.txt<-c("0",paste("(",round(class$brks[2],3),"; ",round(class$brks[3],3),")",sep=""),paste("[",round(class$brks[3],3),"; ",round(class$brks[4],3),")",sep=""),paste("[",round(class$brks[4],3),"; ",round(class$brks[5],3),")",sep=""),paste("[",round(class$brks[5],3),"; ",round(class$brks[6],3),")",sep=""))
legend(-90.3,44,legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")

# If we want to look at the smoothed map for the response variable, that is obtained by taking the estimate 
# of the regression coefficients multiplied by the explanatory variables, to which we add the estimate
# of the spatial random effects, we can look at "model.car$fitted.values". This is a vector with a smoothed value 
# of the response variable for each areal unit. We can create a chloropleht map of the smoothed values as 
# done above.
plotvar <- model.car$fitted.values
nclr <- 5
plotclr <- brewer.pal(nclr,"Purples")
class <- classIntervals(plotvar,nclr,style="jenks")
colcode <- findColours(class,plotclr)

plot(mich.poly,border="black",axes=TRUE,xlim=c(-92,-81))
title(xlab="Longitude",ylab="Latitude",main="Fitted values for \n Log average violent crime arrest rate, years 2011-2012")
plot(mich.poly,col=colcode,add=T)

leg.txt<-c("0",paste("(",round(class$brks[2],3),"; ",round(class$brks[3],3),")",sep=""),paste("[",round(class$brks[3],3),"; ",round(class$brks[4],3),")",sep=""),paste("[",round(class$brks[4],3),"; ",round(class$brks[5],3),")",sep=""),paste("[",round(class$brks[5],3),"; ",round(class$brks[6],3),")",sep=""))
legend(-92,44,legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")
