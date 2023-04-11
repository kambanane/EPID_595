# Fitting a Poisson-Gamma model to disease incidence data
# Fitting a Poisson-Gamma model to disease incidence data in R
# This code shows how to fit the disease mapping model that involves a 
# Poisson likelihood with a Gamma distribution for the relative risks in R.
# For the analysis, we will need the R packages: maps, maptools, spdep, RColorBrewer, classInt, and SpatialEpi. 
# All these packages can be downloaded from a CRAN repository.
# The dataset that we will is the famous lip cancer dataset already loaded in R in the SpatialEpi package.
# In our analysis, we will remove the counties in Scotland that don’t have any neighbor based on Queen’s 
# definition of adjacency.

# knitr::opts_chunk$set(fig.width=6, fig.height=8)

# Installing the necessary packages
# install.packages(c(
#   "maps",
#   "maptools",
#   "spdep",
#   "RColorBrewer",
#   "classInt",
#   "SpatialEpi"
# ),
# repos = "https://cloud.r-project.org")
# 
# ##
## The downloaded binary packages are in
##  /var/folders/9n/gn30q9z56gv9qbdhfqjyfqn80000gn/T//Rtmp5barlF/downloaded_packages
library(maps)
library(maptools)
library(spdep)
library(RColorBrewer)
library(classInt)
library(SpatialEpi)


########################### PART 1 ############################################

# Reading in the data.
data(scotland)
# Looking at the names of the variables in the dataset.
names(scotland)
## [1] "geo"             "data"            "spatial.polygon" "polygon"
# The lip cancer data is stored in the object scotland$data.
scotland.data <- scotland$data
# Here we give R names to different variables in the lip cancer dataset.
Y <- scotland.data$cases
E <- scotland.data$expected
X <- scotland.data$AFF

N <- length(Y)

# The geographical information for lip cancer data is contained in the object called
# polygon. This object contains a list "polygon" that includes all the points
# that are needed to draw the boundaries of each county. The vector "nrepeats"
# lists for each county, how many subpolygons are needed to draw the boundary of the county.

scotland.polygon <- scotland$polygon$polygon
scotland.nrepeats <- scotland$polygon$nrepeats
scotland.names <- scotland$data$county.names
# The command polygon2spatial_polygons transforms the "spatial_polygon" object that is
# a data structure defined within the package SpatialEpi into a polygon object on
# which we can then run functions to obtain the list of adjacency, the number of neighbors
# of each county and so forth.
scotland.spatial.polygon <-
  polygon2spatial_polygon(scotland.polygon,
                          coordinate.system = "+proj=utm",
                          scotland.names,
                          scotland.nrepeats)

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
Y.noisland <- Y[-index.islands]
E.noisland <- E[-index.islands]
# We also remove these 3 counties from the spatial polygon. We need this new polygon
# for plotting purposes, among others.
scotland.noisland.spatial.polygon <-
  scotland.spatial.polygon[seq(1:N)[-index.islands]]
# Here we derive (they will be of use later) the NB object with information on which areal unit
# is neighbor of which areal unit in the new spatial polygon made of the counties of Scotland
# from which we have removed the counties with no neighbors.
scotland.noisland.nb <- poly2nb(scotland.noisland.spatial.polygon)
# We then use the nb2WB function and we get
# the three vectors with the list of which areal unit are adjacent to each areal unit,
# the number of neighbors each areal unit has, and vector (quite uselees) of binary adjacency
# weights.
scotland.noisland.weights <- nb2WB(scotland.noisland.nb)
adj.noisland <- scotland.noisland.weights$adj
weights.noisland <- scotland.noisland.weights$weights
num.noisland <- scotland.noisland.weights$num
N.noisland <- length(num.noisland)


# Making maps of the observed number and of the expected number of lip cancer cases in Scotland.

# Observed number of cases:
plotvar <- Y.noisland
nclr <- 5

plotclr <- brewer.pal(nclr, "Purples")
class <-
  classIntervals(plotvar,
                 nclr,
                 style = "fixed",
                 fixedBreaks = quantile(plotvar, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))
class
## style: fixed
##   one of 4,845 possible partitions of this variable into 5 classes
##   [0,3)   [3,7)   [7,9)  [9,15) [15,39]
##       9      11       9      12      12
colcode <- findColours(class, plotclr)
plot(
  scotland.noisland.spatial.polygon,
  border = "black",
  axes = TRUE,
  xlim = c(-150, 550)
)
title(xlab = "Eastings (km)", ylab = "Northings (km)", main = "Observed number of cases of \n lip cancer in Scotland")
plot(scotland.noisland.spatial.polygon,
     col = colcode,
     add = T)

leg.txt <- c("0-2", "3-6", "7-8", "9-14", "15-39")
legend(
  -150,
  1000,
  legend = leg.txt,
  fill = plotclr,
  cex = 1,
  ncol = 1,
  bty = "n"
)

# Expected number of cases:
plotvar <- E.noisland
nclr <- 5

plotclr <- brewer.pal(nclr, "OrRd")
class <-
  classIntervals(plotvar,
                 nclr,
                 style = "fixed",
                 fixedBreaks = quantile(plotvar, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))
class
## style: fixed
##   one of 163,185 possible partitions of this variable into 5 classes
##   [1.1,3.48)  [3.48,5.58)  [5.58,8.26) [8.26,12.42) [12.42,88.7]
##           11           10           11           10           11
colcode <- findColours(class, plotclr)
plot(
  scotland.noisland.spatial.polygon,
  border = "black",
  axes = TRUE,
  xlim = c(-150, 550)
)
title(xlab = "Eastings (km)", ylab = "Northings (km)", main = "Expected number of cases of \n lip cancer in Scotland")
plot(scotland.noisland.spatial.polygon,
     col = colcode,
     add = T)

leg.txt <-
  c("[1.0; 3.48)",
    "[3.48; 5.58)",
    "[5.58; 8.26)",
    "[8.26; 12.42)",
    "[12.42; 88.7]")
legend(
  -150,
  1000,
  legend = leg.txt,
  fill = plotclr,
  cex = 1,
  ncol = 1,
  bty = "n"
)

# Calculating the standardized mortality ratios (SMRs)
SMR.noisland <- Y.noisland / E.noisland
# Deriving the lower and upper bound of the 95% CI
lb.95ci.SMR <-
  SMR.noisland + qnorm(0.025, 0, 1) * (sqrt(Y.noisland) / E.noisland)
ub.95ci.SMR <-
  SMR.noisland + qnorm(0.975, 0, 1) * (sqrt(Y.noisland) / E.noisland)


# Making a plot of the SMRs
plotvar <- SMR.noisland
nclr <- 5

plotclr <- brewer.pal(nclr, "Greens")
class <-
  classIntervals(plotvar,
                 nclr,
                 style = "fixed",
                 fixedBreaks = quantile(plotvar, c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))
class
## style: fixed
##   one of 211,876 possible partitions of this variable into 5 classes
##        [0,0.3694136) [0.3694136,0.893665)  [0.893665,1.230645)
##                   11                   10                   11
##   [1.230645,2.32634)   [2.32634,6.428571]
##                   10                   11
colcode <- findColours(class, plotclr)

### Make of plot of the SMRs
plot(
  scotland.noisland.spatial.polygon,
  border = "black",
  axes = TRUE,
  xlim = c(-150, 550)
)
title(xlab = "Eastings (km)", ylab = "Northings (km)", main = "Standardized mortality ratio")
plot(scotland.noisland.spatial.polygon,
     col = colcode,
     add = T)

leg.txt <- c("[0; 0.37)",
             "[0.37; 0.89)",
             "[0.89; 1.23)",
             "[1.23; 2.33)",
             "[2.33; 6.43]")
legend(
  -150,
  1000,
  legend = leg.txt,
  fill = plotclr,
  cex = 1,
  ncol = 1,
  bty = "n"
)

# Here we identify which counties in Scotland have SMRs significantly greater than 1 and which
# counties have SMRs significantly lower than 1. After having identified them, we plot them
# in different color: we use red for the counties where the SMRs are significantly greater than 1,
# while we color in blue the counties where the SMRs are significantly lower than 1.
SMR.lt1 <- which(ub.95ci.SMR < 1)
SMR.gt1 <- which(lb.95ci.SMR > 1)

plot(
  scotland.noisland.spatial.polygon,
  border = "black",
  axes = TRUE,
  xlim = c(-150, 550)
)
title(xlab = "Eastings (km)", ylab = "Northings (km)", main = "Counties with SMR significantly different from 1")
plot(
  scotland.noisland.spatial.polygon[SMR.lt1],
  col = rep("dodgerblue", length(SMR.lt1)),
  add = T
)
plot(
  scotland.noisland.spatial.polygon[SMR.gt1],
  col = rep("red1", length(SMR.gt1)),
  add = T
)
leg.txt <- c("SMR significantly > 1", "SMR significantly < 1")
legend(
  -150,
  950,
  legend = leg.txt,
  fill = c("red1", "dodgerblue"),
  cex = 1,
  ncol = 1,
  bty = "n"
)


############################### POISSON-GAMMA MODEL ################################################


# Here we show how to fit a Poisson-Gamma model, where the observed data in an areal unit
# are modeled to follow a Poisson distribution with mean equal to the product of the expected
# number of cases in the areal unti and the relative risk of the disease in the areal unit.
# The relative risks of the disease are then modeled to each follow a Gamma distribution.
# The command that fit this model in the package SpatialEpi is the function eBayes.
# The function eBayes needs as input:
# (i) the observed number of cases; and
# (ii) the expected number of cases.
# Since we are not including any covariate in the model, we write "Xmat=NULL".
Poisson.Gamma.model <- eBayes(Y.noisland, E.noisland, Xmat = NULL)

# The estimated relative risks are contained in the object called RRmed.
# Here we take them and plot them.
plotvar <- Poisson.Gamma.model$RRmed
plotclr <- brewer.pal(nclr, "Greys")
class <-
  classIntervals(
    plotvar,
    nclr,
    style = "fixed",
    fixedBreaks = c(0, 0.37, 0.89, 1.23, 2.33, 4.01)
  )
class
## style: fixed
##   one of 270,725 possible partitions of this variable into 5 classes
##    [0,0.37) [0.37,0.89) [0.89,1.23) [1.23,2.33) [2.33,4.01]
##           5          14          15          12           7
colcode <- findColours(class, plotclr)


plot(
  scotland.noisland.spatial.polygon,
  border = "black",
  axes = TRUE,
  xlim = c(-150, 550)
)
title(xlab = "Eastings (km)", ylab = "Northings (km)", main = "Estimated relative risk of lip cancer \n Poisson-Gamma model")
plot(scotland.noisland.spatial.polygon,
     col = colcode,
     add = T)

leg.txt <- c("[0; 0.37)",
             "[0.37; 0.89)",
             "[0.89; 1.23)",
             "[1.23; 2.33)",
             "[2.33; 4.01]")
legend(
  -150,
  1000,
  legend = leg.txt,
  fill = plotclr,
  cex = 1,
  ncol = 1,
  bty = "n"
)
