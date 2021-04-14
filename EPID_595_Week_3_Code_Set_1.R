# Fitting a Bayesian spatial linear regressions to house sale price data

# Fitting a Bayesian spatial linear regression to point referenced data on house sale price

# This code guides and ask to fit a Bayesian spatial linear regression model to point referenced 
# data in R on house sale prices in Baton Rouge.
# As noted in other instances, you will need the packages spBayes, gstat, fields, and maps to 
# perform exploratory analysis and make maps.

# Additionally we download and install the package akima to obtain an algorithmic method for spatial interpolation.

# The dataset that we will use is a dataset providing information on the sale price (on the log scale) of 
# 70 houses in Baton Rouge, Louisiana.
# Besides information on the log selling price, the dataset also provides information on: 
# the Easting and Northing of the house, the latitude and longitude (for plotting purposes), 
# the size of the living area within the house, the size of other areas in the house, age, number of bedrooms, 
# number of bathrooms, number of half baths.


# Install these packages

#install.packages(c("fields","gstat","maps","spBayes","akima"),repos="https://cloud.r-project.org")

## Call libraries
library(fields)
library(gstat)
library(maps)
library(spBayes)
library(akima)


# Reading in the house sale price-data
baton.data <-
  read.table("Data/Baton_Rouge_houses.csv", sep = ",", header = TRUE)

dim(baton.data)

## [1] 70 11
#Giving a name to all the variables in the dataset that we are going to
# to use in our spatial analysis.

x.coord <- baton.data$Easting
y.coord <- baton.data$Northing
lat <- baton.data$Latitude
lon <- baton.data$Longitude
living <- baton.data$LivingArea
age <- baton.data$Age
other <- baton.data$OtherArea
baths <- baton.data$Bathrooms
logprice <- baton.data$logSellingPr

### Plotting log selling price
# Creating a fictitious grid just for plotting purposes
lon.grid <- seq(min(lon), max(lon), len = 200)
lat.grid <- seq(min(lat), max(lat), len = 200)
z.mat <- matrix(NA, 200, 200)
zlim <- range(logprice, na.rm = T)

col <- color.scale(logprice, col = heat.colors(100)[100:1])
image.plot(
  lon.grid,
  lat.grid,
  z.mat,
  zlim = zlim,
  xlab = "Longitude",
  ylab = "Latitude",
  main = "Log selling prices ",
  xlim = c(range(lon)[1] - 0.1, range(lon)[2] + 0.1),
  ylim = c(range(lat)[1] - 0.01, range(lat)[2] + 0.01)
)
## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning Inf
## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning -
## Inf
points(lon, lat, col = col, pch = 19)

# For this analysis, we are are going to use all the data, except for the last 2 observations.
# We are going to hold these two out, fit the model to the first 68 observations and then use them
# to predict out-of-sample the log selling price of the last two houses.
# Constructing an empirical semi-variogram of the dataframe containing only the first 68
# observations.

# In order to be able to reference to these variables within the variogram command it is safer to give
# a name to each of the vectors with only the first 68 observations.
x.coord.68 <- x.coord[1:68]
y.coord.68 <- y.coord[1:68]
logprice.68 <- logprice[1:68]
living.68 <- living[1:68]
age.68 <- age[1:68]
other.68 <- other[1:68]
baths.68 <- baths[1:68]

price68.df <-
  data.frame(cbind(
    x.coord.68,
    y.coord.68,
    logprice.68,
    living.68,
    age.68,
    other.68,
    baths.68
  ))

price68.df[1:3, ]
##   x.coord.68 y.coord.68 logprice.68 living.68 age.68 other.68 baths.68
## 1   690396.2    3363182    10.87616      1606     11      794        2
## 2   684576.0    3374180    10.54534       979      2      452        1
## 3   681534.1    3376912    10.59663      1144      0      584        2

# The empirical semi-variogram is problematic if it is not truncated. Here we use only distances up to
# 5000m and we use specific bins for distances.

emp.variog.price68 <-
  variogram(
    logprice.68 ~ living.68 + age.68 + other.68 + baths.68,
    locations =  ~ x.coord.68 + y.coord.68,
    data = price68.df,
    cutoff = 5000,
    boundaries = c(5, 10, 40, 70, 500, seq(1500, 5000, by = 500))
  )
plot(emp.variog.price68)

# To convince you that without these choices the empirical semi-variogram looks problematic, construct
# the empirical semi-variogram using the default choices and compare the two empirical semi-variograms.

# Fitting an exponential semi-variogram to the empirical one using initial values informed by the plot

exp.variog.price68 <-
  fit.variogram(emp.variog.price68,
                vgm(
                  psill = 0.02,
                  "Exp",
                  range = 100,
                  nugget = 0.0
                ),
                fit.method = 2)

exp.variog.price68
##   model      psill    range
## 1   Nug 0.00000000   0.0000
## 2   Exp 0.01938659 198.7277

## Partial Sill value that you will use for the scale parameter for the Gamma prior !! #
est.psill.68 <- exp.variog.price68$psill[2]
est.psill.68
## [1] 0.01938659

est.range.68 <- exp.variog.price68$range[2]
est.range.68
## [1] 198.7277

est.nugget.68 <- exp.variog.price68$psill[1]
est.nugget.68
## [1] 0
# Plotting the empirical and the fitted semi-variogram

fitted.variog68 <-
  variogramLine(vgm(est.psill.68, "Exp", est.range.68, est.nugget.68),
                dist_vector = emp.variog.price68$dist)

plot(
  emp.variog.price68$dist,
  emp.variog.price68$gamma,
  xlab = "Distance in meters",
  ylab = "Semi-variance",
  type = "p",
  pch = 19,
  col = "coral",
  ylim = c(0, 0.03)
)
lines(emp.variog.price68$dist, fitted.variog68$gamma, col = "dodgerblue")

# We use the results above to specify the prior distributions of the parameters
# of the spatial regression parameters that we need for our Bayesian analysis.
# We use the package spBayes.
# First we define a matrix where we store all the coordinates of the data, of the first 68
# observations.
coords.68 <- as.matrix(cbind(x.coord.68, y.coord.68),
                       nrow = 68,
                       ncol = 2)

# Then we pass onto specifying priors for the parameters of the spatial linear regression model.

# sigma^2: We specify an Inverse Gamma prior for the partial sill.
# We set the shape to be equal to 2 and the scale parameter to be equal to the WLS estimate of the
# partial sill.
shape.sigma2 <- 2
scale.sigma2 <- est.psill.68 ## Partial sill of the exponential variogram from above

# tau^2: We specify an Inverse Gamma prior for the nugget effect.
# We choose the shape to be equal to 2 and the scale parameter to be equal to its WLS estimate.
# However, since the estimate of nugget effect that we obtained is 0, and the scale parameter HAS to be
# positive, we add a small number of it.
shape.tau2 <- 2
scale.tau2 <- est.nugget.68 + 0.01

# phi: We specify a Uniform prior on the range parameter. Hence, we need to identify the interval
# where this Uniform prior is defined. Remembering that in the package spBayes, phi denotes
# 1/range parameter, rather than the range parameter, we proceed as follows.
# First we determine a range of possible values for the range parameter, and then we invert the
# boundaries of the interval.
# We choose the interval for the range parameter so that in one scenario, the range parameter is
# equal to 5 meters (which means that the correlation vanishes at 5*3=15 meters).
# In the other scenario, the range parameter is equal to the maximum distance between the houses in
# the (full) dataset.
max.dist <-
  max(rdist(matrix(
    cbind(x.coord, y.coord), nrow = 70, ncol = 2
  )))
max.dist
## [1] 28103.33

# Hence we take as prior on phi (=1/range parameter) a Uniform distribution defined on the interval
# (1/max.dist; 1/5)
int.prior.phi <- c(1 / max.dist, 1 / 5)
int.prior.phi
## [1] 3.558298e-05 2.000000e-01
# We use a uniform prior on the entire real line for the regression coefficients. This is also called
# flat prior.

# Finally we need to specify initial values for the spatial covariance function parameters.
# For this, we use the WLS estimates of spatial covariance function parameters, remembering that the
# nugget effect was estimated to be equal to 0 via WLS and that is not an allowed initial value.
sigma2.ini <- est.psill.68
tau2.ini <- est.nugget.68 + 0.01
phi.ini <- 1 / est.range.68

# Now we fit the Bayesian spatial regression model to the log selling price of the first 68 houses.
# We run the algorithm for 50,000 iterations.
# To determine the tuning or size of the neighborhood to use in the Metropolis-Hastings algorithm, it i
# is good to get an idea of the magnitude of each parameter.
# sigma^2 is estimated to be 0.02, which means that maybe a good tuning parameter for sigma^2 is 0.001.
# Similarly phi ranges (a priori) between and , which means that a good tuning parameter
# for phi is . Same applied to tau^2.
price68.bayes <-
  spLM(
    logprice.68 ~ living.68 + age.68 + other.68 + baths.68,
    coords = coords.68,
    starting = list(
      "phi" = phi.ini,
      "sigma.sq" = sigma2.ini,
      "tau.sq" = tau2.ini
    ),
    tuning = list(
      "phi" = 0.1,
      "sigma.sq" = 0.01,
      "tau.sq" = 0.01
    ),
    priors = list(
      "phi.Unif" = int.prior.phi,
      "sigma.sq.IG" = c(2, scale.sigma2),
      "tau.sq.IG" = c(2, scale.tau2),
      "beta.Flat"
    ),
    cov.model = "exponential",
    n.samples = 50000,
    verbose = TRUE,
    n.report = 500
  )


names(price68.bayes)
##  [1] "p.theta.samples" "acceptance"      "Y"               "X"
##  [5] "coords"          "is.pp"           "cov.model"       "nugget"
##  [9] "beta.prior"      "beta.Norm"       "x.names"         "run.time"

# Trace plots of spatial covariance parameters
par(mai = rep(0.4, 4))
plot(price68.bayes$p.theta.samples[, 1:3])

# We not get samples of the regression coefficients.
# Based on the trace plots, we decide to throw away 80% of the iterations that we have obtained
# for burnin.
n.samples <- 50000
burn.in <- 0.80 * n.samples
price68.bayes.other.pars <-
  spRecover(price68.bayes, start = burn.in, verbose = FALSE)

names(price68.bayes.other.pars)

par(mai = rep(0.4, 4), mfrow = c(3, 2))
plot(price68.bayes.other.pars$p.beta.recover.samples[, 1:5])

# Posterior ingerence.
# Here we obtain estimates and 95% credible intervals for the spatial covariance function parameters.
round(summary(price68.bayes.other.pars$p.theta.samples)$quantiles[, c(1, 3, 5)],
      4)

round(summary(price68.bayes.other.pars$p.beta.recover.samples)$quantiles[, c(1, 3, 5)],
      4)

w.summary <-
  apply(t(price68.bayes.other.pars$p.w.recover.samples),
        2,
        quantile,
        c(0.025, 0.5, 0.975))
w.summary

# Creating the interpolated surface of spatial random effects estimates using an
# algorithmic method.
surf.w <- interp(x.coord.68, y.coord.68, as.numeric(w.summary[2, ]))

# Here we create a plot of the interpolated estimates of spatial random effects.
# Black dots denote locations where we have data.
par(mfrow = c(1, 1), mai = c(1, 1, 1, 0.5))
## Run this first
image.plot(
  surf.w$x,
  surf.w$y,
  surf.w$z,
  col = heat.colors(100)[100:1],
  zlim = range(surf.w$z, na.rm = T),
  xlab = "Easting in km",
  ylab = "Northing in km",
  main = "Interpolated posterior median \n of spatial random effects"
)

## Run this second
contour(
  surf.w$x,
  surf.w$y,
  surf.w$z,
  nlevels = 10,
  add = T,
  col = "black"
)
points(x.coord.68, y.coord.68, pch = 19)

## Making predictions.
## Now we predict the log selling price at the last 2 houses.

## This is the matrix with the explanatory variables (living size, age, size of other areas in the
## house, number of baths) for the last 2 hourses in the datasets. The explanatory variables
## should appear in the same order as used in the spLM function.
## Additionally, since the model includes an intercept term, the first 2 columns should be equal to 1
predcov <-
  matrix(cbind(rep(1, 2), living[69:70], age[69:70], other[69:70], baths[69:70]),
         nrow = 2,
         ncol = 5)

# These are the coordinates of the sites where we want to make predictions
pred.coords <- cbind(x.coord[69:70], y.coord[69:70])

# Here we make predictions using all the iterations after burnin.
pred.last2 <-
  spPredict(
    price68.bayes,
    pred.coords = pred.coords,
    pred.covars = predcov,
    start = burn.in,
    thin = 1
  )

# To calculate an estimated predicted value at each of the two locations, we take all the many
# predictions made above and we compute the mean.
post.pred.mean.last2 <-
  apply(pred.last2$p.y.predictive.samples, 1, mean)
post.pred.mean.last2
## [1] 11.15024 11.02179
# We can compare them to the actual observed values at the last 2 houses.
price.last2 <- logprice[69:70]
price.last2
## [1] 11.39392 10.66663
exp(price.last2)

# How did you think the predictions fared to the actual true selling price?
##  Now we are going to perform a similar analysis, except that we are changing the mean model.
## Specifically: our goal is now to fit a Bayesian spatial linear regression model with only
## age as explanatory variable (besides the intercept term).
## We will proceed following the same steps as above, using again only the first 68 observations.

# 1. construct the empirical semi-variogram for the residuals of the linear regression regressing
#    log-selling price for the first 68 houses on the house's age.
#    Start by using the default in R. If the semi-variogram does not display the expected trend that
#    we are used to, change the definition of the bins, and cutoff the semi-variogram after a distance
#    of 5000 meters.

# 2. To the empirical semi-variogram derived in 1, fit an exponential semi-variogram. What are the
#    WLS estimates of the nugget effect, the range parameter and the partial sill?

# 3. Use the WLS estimates to define the priors of the parameters. Note that if the nugget effect
#    is not estimated to be 0, then some of the code we used above will have to be modified (e.g. no
#    need to add 0.01 to the estimate of the nugget effect to define the scale parameter of the Inverse
#    Gamma distribution for tau^2).

#    Define the priors for the spatial covariance function parameters using the same procedure
#    outlined above. Run the MCMC algorithm for 50,000 iterations. To test the influence of the
#    tuning values, use as tuning: 1 for the tuning of phi, 0.1 for sigma^2 and 0.5 for tau^2.
#    What happens to Overall Metropolis Acceptance rate?

# 4. Now use the same tuning values as in the code above. What is the acceptance rate?

# 5. Once you believe you have achieved convergence, derive the estimates of the regression coefficients
#    (the median of the samples post-burnin) and the estimates of the covariance parameters. Use the
#    output for the MCMC with the tuning values that you think works best.
#    Remember the difference in interpretation between phi in the package spBayes (and function spLM)
#    and the range parameter, as wel learned it (phi=1/range parameter). Report the estimates
#    and 95% credible interval for the regression coefficients and the spatial covariance parameter.
#    How would you go about calculating a 90% credible interval?

# 6. Using again the MCMC output with the tuning values that you think works best. Predict the
#    log selling price at the last two homes. Use the median of the posterior samples
#    (obtained after burn-in) as estimate for the predictions. Compare the new predictions you
#    obtained before and the recent ones, and compare it to the true observed values. Which
#    predictions are better?
#
#
#
# 