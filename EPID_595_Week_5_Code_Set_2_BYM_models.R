# Fitting a BYM model to disease incidence data
# Fitting a BYM model to disease incidence data in R

# This code shows how to fit a disease mapping model that involves a Poisson likelihood and
# includes explanatory variables as well as spatial random effects in the model for the log relative risks.
# The code also shows how to fit a Besag-York-Mollie’ or BYM disease mapping model which includes
# both spatial and non-spatial random effects in the model for the log relative risks.

# For the analysis, we will need the R packages: maps, maptools, spdep, RColorBrewer, classInt,
# CARBayes, and SpatialEpi. All these packages can be downloaded from a CRAN repository.
# The dataset that we will is the famous lip cancer dataset already loaded in R in the SpatialEpi package.

# In our analysis, we will remove the counties in Scotland that don’t have any neighbor based on Queen’s
# definition of adjacency.

# knitr::opts_chunk$set(fig.width=6, fig.height=8)
# Installing the necessary packages
# install.packages(
#   c(
#     "maps",
#     "maptools",
#     "spdep",
#     "RColorBrewer",
#     "classInt",
#     "CARBayes",
#     "SpatialEpi"
#   ),
#   repos = "https://cloud.r-project.org"
# )
##
## The downloaded binary packages are in
##  /var/folders/9n/gn30q9z56gv9qbdhfqjyfqn80000gn/T//RtmpkX5Dzt/downloaded_packages
library(maps)
library(maptools)
## Loading required package: sp
## Checking rgeos availability: TRUE
library(spdep)
## Loading required package: spData
## To access larger datasets in this package, install the spDataLarge
## package with: `install.packages('spDataLarge',
## repos='https://nowosad.github.io/drat/', type='source')`
## Loading required package: sf
## Linking to GEOS 3.8.1, GDAL 3.1.4, PROJ 6.3.1
library(RColorBrewer)
library(classInt)
library(CARBayes)
## Loading required package: MASS
## Loading required package: Rcpp
## Registered S3 method overwritten by 'GGally':
##   method from
##   +.gg   ggplot2
library(SpatialEpi)


#### Reading in the data. ####

scotland.spatial.polygon<- st_read("Data/Lip_cancer/scotlip/scotlip.shp")


# 
# 
# data(scotland)
# # Looking at the names of the variables in the dataset.
# names(scotland)
# ## [1] "geo"             "data"            "spatial.polygon" "polygon"
# # The lip cancer data is stored in the object scotland$data.
# scotland.data <- scotland$data
# # Here we give R names to different variables in the lip cancer dataset.
# Y <- scotland.data$cases
# E <- scotland.data$expected
# # The covariate is the proportion of the population that is involved in agriculture, fisheries
# # and forestries.
# X <- scotland.data$AFF
# 
# N <- length(Y)
# 
# # The geographical information for lip cancer data is contained in the object called
# # polygon. This object contains a list "polygon" that includes all the points
# # that are needed to draw the boundaries of each county. The vector "nrepeats"
# # lists for each county, how many subpolygons are needed to draw the boundary of the county.
# 
# scotland.polygon <- scotland$polygon$polygon
# scotland.nrepeats <- scotland$polygon$nrepeats
# scotland.names <- scotland$data$county.names
# # The command polygon2spatial_polygons transforms the "spatial_polygon" object that is
# # a data structure defined within the package SpatialEpi into a polygon object on
# # which we can then run functions to obtain the list of adjacency, the number of neighbors
# # of each county and so forth.
# scotland.spatial.polygon <-
#   polygon2spatial_polygon(scotland.polygon,
#                           coordinate.system = "+proj=utm",
#                           scotland.names,
#                           scotland.nrepeats)

# The function poly2nb derives the list of neighbors of each areal unit. This list of neighbors
# is contained into an object of type nb.
scotland.nb <- poly2nb(scotland.spatial.polygon)

# The function nb2WB transforms the nb object into a list composed of 3 items: an "adj" vector
# that reports the list of neighboring areal units to each areal unit; a "weights" vector
# that reports the binary adjacency vector for the neighboring areal units of an areal unit;
# and a "num" vector that provides the number of neighboring areal units to an areal unit.
scotland.weights <- nb2WB(scotland.nb)
adj <- scotland.weights$adj
weights <- scotland.weights$weights
num <- scotland.weights$num
# Here we identify which counties do not have neighbors.
index.islands <- which(num == 0)
index.islands
## [1]  6  8 11
# There are 3 islands that don't have neighbors. We remove them from the dataset.
#Y.noisland <- Y[-index.islands]
Y.noisland <- scotland.spatial.polygon$CANCER[-index.islands]
E.noisland <- scotland.spatial.polygon$CEXP[-index.islands]
X.noisland <- scotland.spatial.polygon$AFF[-index.islands]
# E.noisland <- E[-index.islands]
# X.noisland <- X[-index.islands]
N.noisland <- length(Y.noisland)

# We also remove these 3 counties from the spatial polygon. We need this new polygon
# for plotting purposes, among others.
scotland.noisland.spatial.polygon <-
  scotland.spatial.polygon[-index.islands,]
# We create the nb object with the information on which areal unit is neighbor of which other
# areal unit, the number of neighbors to each areal unit and a vector with the binary adjacency
# weights.
scotland.noisland.nb <- poly2nb(scotland.noisland.spatial.polygon)
scotland.noisland.weights <- nb2WB(scotland.noisland.nb)
adj.noisland <- scotland.noisland.weights$adj
weights.noisland <- scotland.noisland.weights$weights
num.noisland <- scotland.noisland.weights$num
N.noisland <- length(num.noisland)


# Making maps of the observed number, the expected number of lip cancer cases in Scotland, and
# the explanatory variable (AFF).


tmap_arrange(
  tm_shape(scotland.noisland.spatial.polygon) +
    tm_polygons("CANCER", style = "quantile") +
    tm_layout(main.title = "Cancer cases")
  ,
  tm_shape(scotland.noisland.spatial.polygon) +
    tm_polygons("CEXP", style = "quantile") +
    tm_layout(main.title = "Expected number of cases")
  ,
  tm_shape(scotland.noisland.spatial.polygon) +
    tm_polygons("AFF", style = "quantile") +
    tm_layout(main.title = "Agriculture, Forest and Fisheries (percent)")
)


# Here we fit the disease mapping model where the observed number of cases in an areal unit
# are modeled as following a Poisson distribution with mean equal to the expected number of
# cases in the areal unit times the relative risk in the areal unit.
# In turn, the relative risks in the areal units are modeled to depend on the
# explanatory variables through a linear model that expresses the log relative risk in the
# i-th areal unit as a linear combination of the values of the explanatory variables
# in the i-th areal unit plus a spatial random effect for areal unit i.
# The spatial random effects are modeled to follow a Conditionally AutoRegressive (CAR) model.
# The function that we can use to fit this model is the function S.CARleroux
# which we use also to model spatial areal data for a continuous variable. The function
# is contained in the CARBayes package.

# First we create the adjacency matrix W that has as many rows as many columns as the number
# of counties in Scotland from which the counties without neighbors have been removed.

W <- matrix(0, N.noisland, N.noisland)
rep.scotland.noisland <- rep(1:N.noisland, num.noisland)
for (i in 1:N.noisland) {
  W[i, adj.noisland[rep.scotland.noisland == i]] <-
    rep(1, num.noisland[i])
}

# Here change the format of the explanatory variable on the percent of the population
# that is involved in agriculture, fishery and forestry by multiplying it by 100.
X.noisland2 <- X.noisland * 100

# Having defined the adjacency matrix W, the next step is to define the equation of the model
# that states how the log relative risks of the disease are associated with the explanatory
# variable(s). Because in the disease mapping model, the observed number of cases
# is assumed to follow a Poisson distribution, to avoid incorrect conclusions regarding the
# relative risks, we also need to include in the model an offest which accounts for the size
# of the population  at risk.

formula <- Y.noisland ~ X.noisland2 + offset(log(E.noisland))

### Finally, the command that fits the Poisson model to the data is S.CARleroux.
### The function takes as input, the formula for the log relative risk that we just specified
### above, the adjacency matrix, the number of iterations for the MCMC
### algorithm (n.sample) and the number of iterations discarded for
### burnin (burnin).
model.scotland <-
  S.CARleroux(
    formula = formula,
    family = "poisson",
    W = W,
    burnin = 20000,
    n.sample = 50000,
    rho = 1
  )
## Setting up the model.
## Generating 30000 post burnin and thinned (if requested) samples.
##

## Summarising results.
## Finished in  7.1 seconds.
# A summary table reporting the posterior median, and the boundary of the 95% credible interval
# for all model parameters, along with Geweke diagnostic can be obtained by writing
# the command below.
model.scotland$summary.results
##              Median    2.5%   97.5% n.sample % accept n.effective Geweke.diag
## (Intercept) -0.3204 -0.5594 -0.0745    30000     31.8       912.4        -0.6
## X.noisland2  0.0431  0.0158  0.0690    30000     31.8       696.1         0.0
## tau2         0.4079  0.1999  0.8223    30000    100.0      2121.4         0.4
## rho          1.0000  1.0000  1.0000       NA       NA          NA          NA
# To obtain the samples of the spatial random effects, derive the posterior median of the
# spatial random effects and create a chloropleth map of the estimated spatial random effects,
# we use the approach below.
# The samples for the spatial random effects are contained in the object "sample" under the
# name "phi".
spatial.raneff <- model.scotland$samples$phi
# The matrix with the samples of the spatial random effects has as many rows as the number of MCMC
# iterations post burn-in and as many columns as the number of areal units.
scotland.noisland.spatial.polygon$spatial_raneff <- apply(spatial.raneff, 2, median)

#Here we create a map of the estimated spatial random effects.
tmap_arrange(
  tm_shape(scotland.noisland.spatial.polygon) +
    tm_polygons("spatial_raneff", style = "quantile") +
    tm_layout(main.title = "Spatial random effects")
)



# Here we make a map of the estimated relative risks. To calculate the relative risks, we need
# to take samples of the spatial random effects, and add to them samples of the regression
# coefficients multiplied by the explanatory variables. This sum creates an object with
# as many columns as areal units and as many rows as the number of MCMC iterations post burnin.
# This matrix provides the value of the log relative risk at each areal unit for each
# MCMC iteration after burn-in.
# To obtain the relative risks of the disease, we exponentiate the matrix with the samples of the
# log relative risks, and then we calculate the median across each column.
samples.covariate <-
  model.scotland$samples$beta %*% t(matrix(cbind(rep(1, N.noisland), X.noisland2), N.noisland, 2))
samples.RR <- exp(samples.covariate + spatial.raneff)
scotland.noisland.spatial.polygon$RR<- apply(samples.RR, 2, median)

#Here we create a map of the estimated spatial random effects.
tmap_arrange(
  tm_shape(scotland.noisland.spatial.polygon) +
    tm_polygons("RR", style = "quantile") +
    tm_layout(main.title = "Estimated relative risks")
)


# Here we see how to fit a BYM model to the data. The function that we will use for this task is
# the function S.CARbym in the CARBayes package.
# The syntax is the same as that of the S.CARleroux function.
# Thus, we need to specify the model for the disease incidence data, which remains the same as
# for the S.CARleroux function.
formula <- Y.noisland ~ X.noisland2 + offset(log(E.noisland))

# The BYM disease mapping model is fit through the function S.CARbym. The syntax
# of the function S.CARbym is the same as that of the S.CARleroux function.
model.scotland.bym <-
  S.CARbym(
    formula = formula,
    family = "poisson",
    W = W,
    burnin = 5000,
    n.sample = 10000
  )
## Setting up the model.
## Generating 5000 post burnin and thinned (if requested) samples.
##

## Summarising results.
## Finished in  1.7 seconds.
# To see a summary of the model results, we can access them in the object called
# summary results. This object, as in the previous disease mapping model, contains the
# posterior median, the 95% credible interval and the Geweke's convergence diagnostic.
model.scotland.bym$summary.results
##              Median    2.5%   97.5% n.sample % accept n.effective Geweke.diag
## (Intercept) -0.3442 -0.5407 -0.0880     5000     31.7       164.9         0.2
## X.noisland2  0.0462  0.0164  0.0672     5000     31.7       139.3        -0.5
## tau2         0.3593  0.1474  0.7402     5000    100.0       270.9        -0.5
## sigma2       0.0078  0.0020  0.0499     5000    100.0        43.5         1.4
# To see what are the samples that the function S.CARbym reports as part of its output,
# we can run the function "names" on the "model.scotland$samples" object.
names(model.scotland.bym$samples)
## [1] "beta"   "psi"    "tau2"   "sigma2" "fitted" "Y"
# Here we make a map of the estimated relative risks. The code to calcuate an estimate of the
# relative risk for each areal unit is very similar to the code we used to calculate the relative
# risks in the case of a disease mapping model that only includes spatial random effects.
# The difference between the previous code (for a disease mapping model with only spatial
# random effect) is that when fitting a BYM model, the function S.CARbym does not produce samples
# for the spatial random effects (phi_i) nor for the non-spatial random effects (theta_i),
# but it provides samples for the sum of the two random effects. This sum is
# called psi. Hence, to crease samples for the relative risks, we need
# to take samples of the summed random effects (e.g. samples of the psi_i's), and
# add to them samples of the regression coefficients multiplied by the explanatory variables.
# This sum creates an object with as many columns as areal units and as many rows as the number of
# MCMC iterations post burnin. This matrix provides the value of the log relative risk at each
# areal unit for each MCMC iteration after burn-in.
# To obtain the relative risks of the disease, we exponentiate the matrix with the samples of the
# log relative risks, and then we calculate the median across each column.
samples.covariate.bym <-
  model.scotland.bym$samples$beta %*% t(matrix(cbind(rep(1, N.noisland), X.noisland2), N.noisland, 2))
raneff.bym <- model.scotland.bym$samples$psi
samples.RR.bym <- exp(samples.covariate.bym + raneff.bym)
scotland.noisland.spatial.polygon$RR.bym <- apply(samples.RR.bym, 2, median)

## Now let's look at the BYM derived relative risks alongside our other maps.
tmap_arrange(
  tm_shape(scotland.noisland.spatial.polygon) +
    tm_polygons("CANCER", style = "quantile") +
    tm_layout(main.title = "Cancer cases")
  ,
  tm_shape(scotland.noisland.spatial.polygon) +
    tm_polygons("CEXP", style = "quantile") +
    tm_layout(main.title = "Expected number of cases")
  ,
  tm_shape(scotland.noisland.spatial.polygon) +
    tm_polygons("AFF", style = "quantile") +
    tm_layout(main.title = "Agriculture, Forest and Fisheries (percent)"),
  tm_shape(scotland.noisland.spatial.polygon) +
    tm_polygons("RR.bym", style = "quantile") +
    tm_layout(main.title = "BYM derived relative risks"),
  ncol = 2
)


