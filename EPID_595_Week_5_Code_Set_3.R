# Finding clusters and hotspots in disease incidence data
# Finding clusters and hotspots in disease incidence data in R

# This code shows how to identify clusters and hoptspots in disease incidence data using disease mapping models.
# Specifically, we focus on the two disease mapping models that model the observed number of cases of a disease
# as following a Poisson distribution with mean equal to the expected number of cases times a log relative risk.

# The relative risk, in one case, is modeled as a linear combination of explanatory variables and spatial random
# effects. This is the disease mapping model with CAR spatial random effects.

# The second disease mapping model is Besag-York-Mollie’ or the BYM disease mapping model which includes both spatial and non-spatial random effects in the model for the log relative risks.
# For the analysis, we will need the R packages: maps, maptools, spdep, RColorBrewer, classInt, and CARBayes. All these packages can be downloaded from a CRAN repository.
# The dataset that we will use is a dataset reporting the number of cases of lung cancer in Pennsylvania in 2002 for different strata of the population. Here we focus on 60-69 years olds.
# knitr::opts_chunk$set(fig.width=6, fig.height=8)
# Installing the necessary packages
# install.packages(c(
#   "maps",
#   "maptools",
#   "spdep",
#   "RColorBrewer",
#   "classInt",
#   "CARBayes"
# ),
# repos = "https://cloud.r-project.org")
# ##
## The downloaded binary packages are in
##  /var/folders/9n/gn30q9z56gv9qbdhfqjyfqn80000gn/T//RtmpzLNsBP/downloaded_packages
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
# Reading in the data on lung cancer in counties in Pennsylvania.
penn.lc <-
  read.table("Data/Penn_lung_cancer_60-69yo.csv",
             sep = ",",
             header = T)
# Looking at the names of the variables in the dataset.
names(penn.lc)
## [1] "county"  "cases"   "E"       "race"    "gender"  "smoking"
# Here we give R names to different variables in the lung cancer dataset.
Y.all <- penn.lc$cases
E.all <- penn.lc$E
race <- penn.lc$race
gender <- penn.lc$gender
smoking <- penn.lc$smoking

# Here we see what are the demographic strata that the dataset contains.
table(race)
## race
##   o   w
## 134 134
table(gender)
## gender
##   f   m
## 134 134
# Since the data refers to different race and gender groups, we are going to
# create different variables that report the number of observed cases and the
# number of expected cases of lung cancer for the different subsets of the population
Y.o.f <-
  Y.all[which(as.character(race) == "o" & as.character(gender) == "f")]
Y.o.m <-
  Y.all[which(as.character(race) == "o" & as.character(gender) == "m")]
Y.w.f <-
  Y.all[which(as.character(race) == "w" & as.character(gender) == "f")]
Y.w.m <-
  Y.all[which(as.character(race) == "w" & as.character(gender) == "m")]
E.o.f <-
  E.all[which(as.character(race) == "o" & as.character(gender) == "f")]
E.o.m <-
  E.all[which(as.character(race) == "o" & as.character(gender) == "m")]
E.w.f <-
  E.all[which(as.character(race) == "w" & as.character(gender) == "f")]
E.w.m <-
  E.all[which(as.character(race) == "w" & as.character(gender) == "m")]

# For this analysis, we focus on lung cancer in white female.
# Calculating the standardized mortality ratio of lung cancer for white female.
SMR.w.f <- Y.w.f / E.w.f

# Here we extract the geographical information for the state of Pennsylvania and its counties.
# We use the maps package to get the polygon that correspond to counties in Pennsylvania.
penn.county <- map("county", "pennsylvania", fill = T, plot = F)
penn.IDs <- sapply(strsplit(penn.county$names, ":"), function(x)
  x[1])
penn.poly <-
  map2SpatialPolygons(penn.county,
                      IDs = penn.IDs,
                      proj4string = CRS("+proj=longlat +datum=WGS84"))


# Making a map of the SMR for white females in Pennsylvania, and for the explanatory variable
# representing the exposure of interest: the proportion of the population that smokes.

# SMR of lung cancer for white, females of 60-69 years of age.
plotvar <- SMR.w.f
nclr <- 5
plotclr <- brewer.pal(nclr, "Purples")
class <- classIntervals(plotvar, nclr, style = "jenks")
colcode <- findColours(class, plotclr)

plot(
  penn.poly,
  border = "black",
  axes = TRUE,
  xlim = c(-82, -74),
  ylim = c(39, 43.5)
)
title(xlab = "Longitude", ylab = "Latitude", main = "SMR white, female")
plot(penn.poly, col = colcode, add = T)

leg.txt <-
  c(
    "0",
    paste("(", round(class$brks[2], 3), "; ", round(class$brks[3], 3), ")", sep =
            ""),
    paste("[", round(class$brks[3], 3), "; ", round(class$brks[4], 3), ")", sep =
            ""),
    paste("[", round(class$brks[4], 3), "; ", round(class$brks[5], 3), ")", sep =
            ""),
    paste("[", round(class$brks[5], 3), "; ", round(class$brks[6], 3), ")", sep =
            "")
  )
legend(
  -81,
  43.5,
  legend = leg.txt,
  fill = plotclr,
  cex = 1,
  ncol = 2,
  bty = "n"
)

# Proportion of population that self-identifies as smokers in Pennsylvania counties in year 2002.
# Since smoking is a variable that is defined for various strata and we just want to plot it
# for the strata corresponding to white, females, we proceed as follows.
smoking.w.f <-
  smoking[which(as.character(race) == "w" &
                  as.character(gender) == "f")]
plotvar <- smoking.w.f
nclr <- 5
plotclr <- brewer.pal(nclr, "Reds")
class <- classIntervals(plotvar, nclr, style = "jenks")
colcode <- findColours(class, plotclr)

plot(
  penn.poly,
  border = "black",
  axes = TRUE,
  xlim = c(-82, -74),
  ylim = c(39, 43.5)
)
title(xlab = "Longitude", ylab = "Latitude", main = "Smoking")
plot(penn.poly, col = colcode, add = T)

leg.txt <-
  c(
    "0",
    paste("(", round(class$brks[2], 3), "; ", round(class$brks[3], 3), ")", sep =
            ""),
    paste("[", round(class$brks[3], 3), "; ", round(class$brks[4], 3), ")", sep =
            ""),
    paste("[", round(class$brks[4], 3), "; ", round(class$brks[5], 3), ")", sep =
            ""),
    paste("[", round(class$brks[5], 3), "; ", round(class$brks[6], 3), ")", sep =
            "")
  )
legend(
  -81,
  43.5,
  legend = leg.txt,
  fill = plotclr,
  cex = 1,
  ncol = 2,
  bty = "n"
)

# Here we identify clusters and hotspots using the results of a disease mapping model with
# spatial random effects modeled to follow a Conditionally AutoRegressive (CAR) model.
# To fit this model we use the function S.CARleroux in the CARBayes package.

# We start by creating he adjacency matrix W that has as many rows as many columns as the number
# of counties in Pennsylvania.
# First we need to create the object nb that contains the information of which county
# is neighbor to which county.
penn.nb <- poly2nb(penn.poly)
penn.nb
## Neighbour list object:
## Number of regions: 67
## Number of nonzero links: 346
## Percentage nonzero weights: 7.70773
## Average number of links: 5.164179
# Having derived the nb onject we can get the vectors with the list of counties that
# are adjacent to each county in Pennsylvania using the command nb2WB.
penn.weights <- nb2WB(penn.nb)
adj <- penn.weights$adj
weights <- penn.weights$weights
num <- penn.weights$num
N <- length(num)

# Here we create the adjacency matrix W.
W <- matrix(0, N, N)
rep.penn <- rep(1:N, num)
for (i in 1:N) {
  W[i, adj[rep.penn == i]] <- rep(1, num[i])
}

# Fitting the disease mapping model with spatial random effects modeled as having a CAR model.
formula <- Y.w.f ~ smoking.w.f + offset(log(E.w.f))
model.penn.w.f <-
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
## Finished in  7.9 seconds.
# Here are the results' summary. Based on these results there is an association that is mostly
# positive between smoking and risk of lung cancer in white females 60-69 years old.
model.penn.w.f$summary.results
##              Median    2.5%   97.5% n.sample % accept n.effective Geweke.diag
## (Intercept) -0.7129 -1.4005 -0.0227    30000     36.3      1702.6         0.4
## smoking.w.f  2.8203 -0.0348  5.7148    30000     36.3      1879.7         0.0
## tau2         0.0190  0.0036  0.0738    30000    100.0       303.5        -1.9
## rho          1.0000  1.0000  1.0000       NA       NA          NA          NA
# Making a map of the estimated spatial random effects. To calculate the posterior median
# of the spatial random effects, we first take the posterior sample for the spatoal random effects,
# and then we calculate the posterior median of the samples.
spat.raneff <- model.penn.w.f$samples$phi

post.median.phi <- apply(spat.raneff, 2, median)
plotvar <- post.median.phi
nclr <- 5

range(plotvar)
## [1] -0.06472656  0.16937937
plotclr <- brewer.pal(nclr, "RdBu")[5:1]
class <-
  classIntervals(
    plotvar,
    nclr,
    style = "fixed",
    fixedBreaks = c(-0.07, -0.03, 0, 0.05, 0.10, 0.20)
  )
class
## style: fixed
##   one of 720,720 possible partitions of this variable into 5 classes
## [-0.07,-0.03)     [-0.03,0)      [0,0.05)    [0.05,0.1)     [0.1,0.2]
##            17            28            11             6             5
colcode <- findColours(class, plotclr)

plot(
  penn.poly,
  border = "black",
  axes = TRUE,
  xlim = c(-82, -74),
  ylim = c(39, 43.5)
)
title(xlab = "Longitude", ylab = "Latitude", main = "Estimated spatial random effects")
plot(penn.poly, col = colcode, add = T)

leg.txt <-
  c("[-0.07; -0.03)",
    "[-0.03; 0)",
    "[0; 0.05)",
    "[0.05; 0.10)",
    "[0.10; 0.20)")
legend(
  -81,
  43.5,
  legend = leg.txt,
  fill = plotclr,
  cex = 1,
  ncol = 2,
  bty = "n"
)

# To identify clusters and hotspots based on the spatial random effects, we calculate 95%
# credible intervals for the spatial random effects for each county in Pennsylvania and
# we determine which county has a spatial random effect whose credible interval is
# either completely lower than 0 or completely greater than 0.
# The credible intervals for the spatial random effects are calculated by taking the 2.75-th
# and the 97.5-th percentiles of the samples of the spatial random effects.

post.lb95.phi <- apply(spat.raneff, 2, quantile, 0.025)
plotvar <- post.lb95.phi
nclr <- 5

range(plotvar)
## [1] -0.29855930  0.01111612
plotclr <- brewer.pal(nclr, "RdBu")[5:1]
class <-
  classIntervals(
    plotvar,
    nclr,
    style = "fixed",
    fixedBreaks = c(-0.31, -0.15, 0, 0.005, 0.01, 0.025)
  )
class
## style: fixed
##   one of 720,720 possible partitions of this variable into 5 classes
## [-0.31,-0.15)     [-0.15,0)     [0,0.005)  [0.005,0.01)  [0.01,0.025]
##            49            16             1             0             1
colcode <- findColours(class, plotclr)

plot(
  penn.poly,
  border = "black",
  axes = TRUE,
  xlim = c(-82, -74),
  ylim = c(39, 43.5)
)
title(xlab = "Longitude", ylab = "Latitude", main = "Lower bound 95% CI spatial random effects")
plot(penn.poly, col = colcode, add = T)

leg.txt <-
  c("[-0.31; -0.15)",
    "[-0.15; 0)",
    "[0; 0.005)",
    "[0.005; 0.01)",
    "[0.01; 0.025)")
legend(
  -81,
  43.5,
  legend = leg.txt,
  fill = plotclr,
  cex = 1,
  ncol = 2,
  bty = "n"
)

#
post.ub95.phi <- apply(spat.raneff, 2, quantile, 0.975)
plotvar <- post.ub95.phi
nclr <- 5

range(plotvar)
## [1] 0.05744437 0.39083668
plotclr <- brewer.pal(nclr, "Reds")
class <-
  classIntervals(
    plotvar,
    nclr,
    style = "fixed",
    fixedBreaks = c(0.05, 0.10, 0.15, 0.20, 0.30, 0.40)
  )
class
## style: fixed
##   one of 720,720 possible partitions of this variable into 5 classes
## [0.05,0.1) [0.1,0.15) [0.15,0.2)  [0.2,0.3)  [0.3,0.4]
##         10         22         20          9          6
colcode <- findColours(class, plotclr)

plot(
  penn.poly,
  border = "black",
  axes = TRUE,
  xlim = c(-82, -74),
  ylim = c(39, 43.5)
)
title(xlab = "Longitude", ylab = "Latitude", main = "Upper bound 95% CI spatial random effects")
plot(penn.poly, col = colcode, add = T)

leg.txt <-
  c("[0.05; 0.10)",
    "[0.10; 0.15)",
    "[0.15; 0.20)",
    "[0.20; 0.30)",
    "[0.30; 0.40)")
legend(
  -81,
  43.5,
  legend = leg.txt,
  fill = plotclr,
  cex = 1,
  ncol = 2,
  bty = "n"
)

# This code identifies counties for which the spatial random effect is estimated to be
# lower (resp. greater) than 0 with 95% probability given the data.
hot.spat.raneff <- which(post.lb95.phi > 0 & post.ub95.phi > 0)
cold.spat.raneff <- which(post.lb95.phi < 0 & post.ub95.phi < 0)


# Plotting counties for which the spatial random effect is estimated to be positive, that is,
# counties for which the relative risk of lung cancer for white females, aged 60-69,
# have a larger risk of lung cancer that expected given the proportion of population
# in the county that smokes.
plot(
  penn.poly,
  border = "black",
  axes = TRUE,
  xlim = c(-82, -74),
  ylim = c(39, 43.5)
)
title(xlab = "Longitude", ylab = "Latitude", main = "Counties with positive spatial random effects")
plot(penn.poly[hot.spat.raneff], col = "red", add = T)

# Here we identify clusters and hotspots based on the estimated relative risk of lung cancer,
# whether the relative risks (RRs) is significantly greater than 1.
# To estimate the RRs we need to create posterior samples for the relative risks.
# We can generate them easily by simply taking the posterior samples for the spatial random effects
# and posterior samples obtained by multiplying the sampled regression coefficients tumes
# the matrix with the explanatory variables.
RR.random <- model.penn.w.f$samples$phi
RR.covariate <-
  model.penn.w.f$samples$beta %*% t(matrix(cbind(rep(1, N), smoking.w.f), N, 2))
RR.all <- exp(RR.covariate + RR.random)

# Estimating the relative risks via the posterior median.
post.median.RR <- apply(RR.all, 2, median)
plotvar <- post.median.RR
nclr <- 5
range(plotvar)
## [1] 0.8450948 1.2674452
plotclr <- brewer.pal(nclr, "YlOrRd")
class <-
  classIntervals(plotvar,
                 nclr,
                 style = "fixed",
                 fixedBreaks = c(0.8, 0.9, 1.0, 1.1, 1.2, 1.3))
class
## style: fixed
##   one of 720,720 possible partitions of this variable into 5 classes
## [0.8,0.9)   [0.9,1)   [1,1.1) [1.1,1.2) [1.2,1.3]
##        13        31        22         0         1
colcode <- findColours(class, plotclr)

plot(
  penn.poly,
  border = "black",
  axes = TRUE,
  xlim = c(-82, -74),
  ylim = c(39, 43.5)
)
title(xlab = "Longitude", ylab = "Latitude", main = "Estimated relative risk \n Lung cancer, white female 60-69 y.o.")
plot(penn.poly, col = colcode, add = T)

leg.txt <-
  c("[0.8; 0.9)",
    "[0.9; 1.0)",
    "[1.0; 1.1)",
    "[1.1; 1.2)",
    "[1.2; 1.3)")
legend(
  -81,
  43.5,
  legend = leg.txt,
  fill = plotclr,
  cex = 1,
  ncol = 2,
  bty = "n"
)

# Here we determine the clusters and hotspots by looking at counties for which the 95%
# credible interval for the relative risks is all above 1.
# Determining the lower and upper bound of the 95% credible interval for the relative risk
# of lung cancer for white females, aged 60-69 years.
# Deriving the lower bound.
post.lb95.RR <- apply(RR.all, 2, quantile, 0.025)
plotvar <- post.lb95.RR
nclr <- 5

range(plotvar)
## [1] 0.6702196 1.0608808
plotclr <- brewer.pal(nclr, "RdBu")[5:1]
class <-
  classIntervals(plotvar,
                 nclr,
                 style = "fixed",
                 fixedBreaks = c(0.6, 0.7, 0.8, 0.9, 1.0, 1.1))
class
## style: fixed
##   one of 720,720 possible partitions of this variable into 5 classes
## [0.6,0.7) [0.7,0.8) [0.8,0.9)   [0.9,1)   [1,1.1]
##         7        32        25         2         1
colcode <- findColours(class, plotclr)

plot(
  penn.poly,
  border = "black",
  axes = TRUE,
  xlim = c(-82, -74),
  ylim = c(39, 43.5)
)
title(xlab = "Longitude", ylab = "Latitude", main = "Lower bound 95% CI relative risks")
plot(penn.poly, col = colcode, add = T)

leg.txt <-
  c("[0.6; 0.7)",
    "[0.7; 0.8)",
    "[0.8; 0.9)",
    "[0.9; 1.0)",
    "[1.0; 1.1)")
legend(
  -81,
  43.5,
  legend = leg.txt,
  fill = plotclr,
  cex = 1,
  ncol = 2,
  bty = "n"
)

# Deriving the upper bound
post.ub95.RR <- apply(RR.all, 2, quantile, 0.975)
plotvar <- post.ub95.RR
nclr <- 5

range(plotvar)
## [1] 1.001992 1.534370
plotclr <- brewer.pal(nclr, "Reds")
class <-
  classIntervals(
    plotvar,
    nclr,
    style = "fixed",
    fixedBreaks = c(0.99, 1.1, 1.2, 1.3, 1.4, 1.55)
  )
class
## style: fixed
##   one of 720,720 possible partitions of this variable into 5 classes
## [0.99,1.1)  [1.1,1.2)  [1.2,1.3)  [1.3,1.4) [1.4,1.55]
##         27         16         17          6          1
colcode <- findColours(class, plotclr)

plot(
  penn.poly,
  border = "black",
  axes = TRUE,
  xlim = c(-82, -74),
  ylim = c(39, 43.5)
)
title(xlab = "Longitude", ylab = "Latitude", main = "Upper bound 95% CI relative risks")
plot(penn.poly, col = colcode, add = T)

leg.txt <-
  c("[0.99; 1.1)",
    "[1.1; 1.2)",
    "[1.2; 1.3)",
    "[1.3; 1.4)",
    "[1.4; 1.55)")
legend(
  -81,
  43.5,
  legend = leg.txt,
  fill = plotclr,
  cex = 1,
  ncol = 2,
  bty = "n"
)

# Identifying coldspots and hotspots, e.g. counties where the 95% credible interval
# for the relative risk is, respectively, below 1 (coldspots), and counties where the 95%
# credible interval for the relative risk is all above 1.
coldspot.rr <- which(post.lb95.RR < 1 & post.ub95.RR < 1)
hotspot.rr <- which(post.lb95.RR > 1 & post.ub95.RR > 1)


plot(
  penn.poly,
  border = "black",
  axes = TRUE,
  xlim = c(-82, -74),
  ylim = c(39, 43.5)
)
title(xlab = "Longitude", ylab = "Latitude", main = "Counties with relative risk \n of lung cancer significantly different from 1")
plot(penn.poly[coldspot.rr], col = "blue", add = T)
plot(penn.poly[hotspot.rr], col = "red", add = T)

leg.txt <- c("RR significantly < 1", "RR significantly > 1")
legend(
  -81,
  43.5,
  legend = leg.txt,
  fill = c("blue", "red"),
  cex = 1,
  ncol = 1,
  bty = "n"
)

# Now we see how to identify hotspots and coldspots by fitting a BYM model to the data.
# In this case the hotspots and coldspots are identified only based on the relative risks,
# since we have discussed how in the BYM model the spatial and non-spatial random effects
# cannot be estimated separately.

# Fitting a BYM model to the data.
formula <- Y.w.f ~ smoking.w.f + offset(log(E.w.f))
model.bym.penn.w.f <-
  S.CARbym(
    formula = formula,
    family = "poisson",
    W = W,
    burnin = 80000,
    n.sample = 100000
  )
## Setting up the model.
## Generating 20000 post burnin and thinned (if requested) samples.
##

## Summarising results.
## Finished in  16.6 seconds.
model.bym.penn.w.f$summary.results
##              Median    2.5%  97.5% n.sample % accept n.effective Geweke.diag
## (Intercept) -0.7079 -1.4647 0.0332    20000     34.9      1026.4         0.8
## smoking.w.f  2.7967 -0.3158 5.9282    20000     34.9      1068.3        -1.0
## tau2         0.0158  0.0033 0.0710    20000    100.0       231.2         0.1
## sigma2       0.0055  0.0019 0.0225    20000    100.0       241.0         1.3
# Estimated relative risks
RR.random <- model.bym.penn.w.f$samples$psi
RR.covariate <-
  model.bym.penn.w.f$samples$beta %*% t(matrix(cbind(rep(1, N), smoking.w.f), N, 2))
RR.all <- exp(RR.covariate + RR.random)

post.median.RR <- apply(RR.all, 2, median)
plotvar <- post.median.RR
nclr <- 5
range(plotvar)
## [1] 0.8224604 1.2765376
plotclr <- brewer.pal(nclr, "YlOrRd")
class <-
  classIntervals(plotvar,
                 nclr,
                 style = "fixed",
                 fixedBreaks = c(0.8, 0.9, 1.0, 1.1, 1.2, 1.3))
class
## style: fixed
##   one of 720,720 possible partitions of this variable into 5 classes
## [0.8,0.9)   [0.9,1)   [1,1.1) [1.1,1.2) [1.2,1.3]
##        19        26        20         1         1
colcode <- findColours(class, plotclr)

plot(
  penn.poly,
  border = "black",
  axes = TRUE,
  xlim = c(-82, -74),
  ylim = c(39, 43.5)
)
title(xlab = "Longitude", ylab = "Latitude", main = "Estimated relative risk \n Lung cancer, white female 60-69 y.o. \n BYM model")
plot(penn.poly, col = colcode, add = T)

leg.txt <-
  c("[0.8; 0.9)",
    "[0.9; 1.0)",
    "[1.0; 1.1)",
    "[1.1; 1.2)",
    "[1.2; 1.3)")
legend(
  -81,
  43.5,
  legend = leg.txt,
  fill = plotclr,
  cex = 1,
  ncol = 2,
  bty = "n"
)

# 95% credible interval for the estimated relative risks..
post.lb95.RR <- apply(RR.all, 2, quantile, 0.025)
plotvar <- post.lb95.RR
nclr <- 5

range(plotvar)
## [1] 0.6176619 1.0576242
plotclr <- brewer.pal(nclr, "RdBu")[5:1]
class <-
  classIntervals(plotvar,
                 nclr,
                 style = "fixed",
                 fixedBreaks = c(0.6, 0.7, 0.8, 0.9, 1.0, 1.1))
class
## style: fixed
##   one of 720,720 possible partitions of this variable into 5 classes
## [0.6,0.7) [0.7,0.8) [0.8,0.9)   [0.9,1)   [1,1.1]
##        24        27        14         1         1
colcode <- findColours(class, plotclr)

plot(
  penn.poly,
  border = "black",
  axes = TRUE,
  xlim = c(-82, -74),
  ylim = c(39, 43.5)
)
title(xlab = "Longitude", ylab = "Latitude", main = "Lower bound 95% CI relative risks \n BYM model")
plot(penn.poly, col = colcode, add = T)

leg.txt <-
  c("[0.6; 0.7)",
    "[0.7; 0.8)",
    "[0.8; 0.9)",
    "[0.9; 1.0)",
    "[1.0; 1.1)")
legend(
  -81,
  43.5,
  legend = leg.txt,
  fill = plotclr,
  cex = 1,
  ncol = 2,
  bty = "n"
)

# Upper bound 95% creidble intēerval for the relative risk
post.ub95.RR <- apply(RR.all, 2, quantile, 0.975)
plotvar <- post.ub95.RR
nclr <- 5

range(plotvar)
## [1] 1.025240 1.541444
plotclr <- brewer.pal(nclr, "Reds")
class <-
  classIntervals(plotvar,
                 nclr,
                 style = "fixed",
                 fixedBreaks = c(1.0, 1.1, 1.2, 1.3, 1.4, 1.55))
class
## style: fixed
##   one of 720,720 possible partitions of this variable into 5 classes
##    [1,1.1)  [1.1,1.2)  [1.2,1.3)  [1.3,1.4) [1.4,1.55]
##          9         29         12         15          2
colcode <- findColours(class, plotclr)

plot(
  penn.poly,
  border = "black",
  axes = TRUE,
  xlim = c(-82, -74),
  ylim = c(39, 43.5)
)
title(xlab = "Longitude", ylab = "Latitude", main = "Upper bound 95% CI relative risks \n BYM model")
plot(penn.poly, col = colcode, add = T)

leg.txt <-
  c("[1.0; 1.1)",
    "[1.1; 1.2)",
    "[1.2; 1.3)",
    "[1.3; 1.4)",
    "[1.4; 1.55)")
legend(
  -81,
  43.5,
  legend = leg.txt,
  fill = plotclr,
  cex = 1,
  ncol = 2,
  bty = "n"
)

cold.RR.bym <- which(post.lb95.RR < 1 & post.ub95.RR < 1)
hot.RR.bym <- which(post.lb95.RR > 1 & post.ub95.RR > 1)


plot(
  penn.poly,
  border = "black",
  axes = TRUE,
  xlim = c(-82, -74),
  ylim = c(39, 43.5)
)
title(xlab = "Longitude", ylab = "Latitude", main = "Counties with relative risk \n of lung cancer significantly greater than 1 \n BYM")
plot(penn.poly[cold.RR.bym], col = "red", add = T)
plot(penn.poly[hot.RR.bym], col = "red", add = T)
