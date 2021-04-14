# 
# ---
#   title: "Spatial analysis Michigan burglaries"
# output: rmarkdown::html_vignette
# fig_width: 8
# fig_height: 10
# ---
  
#   ## Performing a spatial analysis of burglaries occurring in Michigan in 2006
#   
#   * This code guides and asks you to fit a spatial regression model with spatial random effects modeled using a Conditionally AutoRegressive (CAR) model to a dataset reporting the rate of burglaries in Michigan counties in year 2006.
# * For the analysis, we will need the R packages:CARBayes, maps, maptools, spdep, RColorBrewer, and classInt. All these packages can be downloaded from a CRAN repository.
# * The dataset that we will use for this analysis is a dataset providing information on the rates of burglaries per 1,000 people in Michigan during year 2006.
# * Besides information on the rate of burglaries per 1,000 people in Michigan during year 2006, the dataset also reports information on unemployment rate and on the the median household income in 1,000 dollars. 

# Installing the necessary packages
install.packages(c("maps","maptools","spdep","RColorBrewer","classInt","CARBayes"),repos="https://cloud.r-project.org")
library(maps)
library(maptools)
library(spdep)
library(RColorBrewer)
library(classInt)
library(CARBayes)


# Reading in the data.
mich.crime <- read.table("Michigan_crime_2006.csv",sep=",",header=T)
# To look at the names of the variables in the Michigan crime dataset, we can simply do the following:
names(mich.crime)

y <- mich.crime$burglary_rate
unempt <- mich.crime$unempt_rate
income <- mich.crime$median_income_1000s


# Q1: Create a histogram of the outcome variable, the burglary rate per 1,000 people in Michigan in 2006.
# Does it seem like a reasonable assumption to use a Gaussian distribution to model the response variable?
hist(y,breaks=30,xlab="",col="coral",main="")

#Q2: If you decide to transform the data to get a new response variable that is more approximately 
# normally distributed, what transformation would you use?


# Here you are asked to create chloropleth maps of the response variable (the burglary rate per 1,000 people) and of the explanatory variables (unemployment rate and median household income in 1,000 dollars).
mich.county <-map("county","michigan",fill=T,plot=F)
mich.IDs <- sapply(strsplit(mich.county$names,":"),function(x) x[1])
mich.poly <- map2SpatialPolygons(mich.county,IDs=mich.IDs,proj4string=CRS("+proj=longlat +datum=WGS84"))

# Plotting the response variable: (work on the original scale). Choose the number of different color shades,
# the color palette that you prefer.
plotvar <-
  nclr <- 
  plotclr <- brewer.pal(nclr,"Purples")
class <- classIntervals(plotvar,nclr,style="jenks")
colcode <- findColours(class,plotclr)

# Fix the name of labels and the main title accordingly
plot(mich.poly,border="black",axes=TRUE,xlim=c(-92,-81))
title(xlab="",ylab="",main="")
plot(mich.poly,col=colcode,add=T)

# Fix the legend to reflect the number of colors that you have chosen to use
leg.txt<-c(paste("(",round(class$brks[1],3),"; ",round(class$brks[2],3),")",sep=""),paste("(",round(class$brks[2],3),"; ",round(class$brks[3],3),")",sep=""),paste("[",round(class$brks[3],3),"; ",round(class$brks[4],3),")",sep=""),...)
legend(-92,44,legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")


# Plotting the explanatory variables: unemployment rate and median household income in 1,000 dollars.
# Use as many number of levels as you retain best.
plotvar <- 
  nclr <- 
  plotclr <- brewer.pal(nclr,"Greens")
class <- classIntervals(plotvar,nclr,style="jenks")
colcode <- findColours(class,plotclr)

plot(mich.poly,border="black",axes=TRUE,xlim=c(-92,-81))
title(xlab="",ylab="",main="")
plot(mich.poly,col=colcode,add=T)

leg.txt<-c(paste("[",round(class$brks[1],3),"; ",round(class$brks[2],3),")",sep=""),paste("[",round(class$brks[2],3),"; ",round(class$brks[3],3),")",sep=""),....)
legend(-92,44,legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")


# 
plotvar <- 
  nclr <- 
  plotclr <- brewer.pal(nclr,"Reds")
class <- classIntervals(plotvar,nclr,style="jenks")
colcode <- findColours(class,plotclr)

plot(mich.poly,border="black",axes=TRUE,xlim=c(-92,-81))
title(xlab="Longitude",ylab="Latitude",main="Unemployment rate in 2011")
plot(mich.poly,col=colcode,add=T)

leg.txt<-c(paste("[",round(class$brks[1],3),"; ",round(class$brks[2],3),")",sep=""),paste("[",round(class$brks[2],3),"; ",round(class$brks[3],3),")",sep=""),...)
legend(-92,44,legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")

# Here we compute Moran's I on the response variable directly and on the residuals of a linear regression regressing 
# the rate of burglaries per 1,000 people on unemployment rate and income.
mi.nb <- poly2nb(mich.poly)
mi.nb

# Q3: Calculate Moran's I for the rate of burglaries per 1,000 people. Is the value of Moran's I significantly
# different from the value we would expect if the data were independent?

# Q4: Now regress the rate of burglaries per 1,000 people on unemployment rate and the median household income 
# in 1,000 dollars. Derive the residuals and calculate Moran's I on the residuals.
# Is the value of Moran's I on the residuals significantly different from the value we would expect if the 
# data were independent?
lm.crime <- lm(y~unempt+income)
res.lm.crime <- as.numeric(lm.crime$residuals)
```

```{r}
# Now we fit a spatial linear regression model to the rate of burglaries per 1,000 people in Michigan in 2016.
# The model includes the explanatory variables (unemployment rate and income in 1,000 dollars) and spatial random 
# effects that are modeled using a CAR model.
# Use the function S.CARleroux to fit the model. Run enough MCMC iterations to ensure convergence.
N <- length(y)
# Creating the adjaceny matrix W.
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

# Here, we use all the defaults on the prior distributions.
formula <- 
  model.car <- S.CARleroux(formula=formula, W=W, family="gaussian",
                           burnin=, n.sample=,thin=1, prior.mean.beta=NULL, prior.var.beta=NULL,prior.nu2=NULL,prior.tau2=NULL,rho=1, verbose=TRUE)

# Q5. Provide the estimates of the regression coefficients and of the spatial and non-spatial variance.
# What is the association between unemployment rate and burglary rate? What about the association between median
# household income and burglary rate? Are they in the direction you expected them?
model.car$summary.results



# Here we calculate the posterior median of the spatial random effects and we make a plot
# of the estimated spatial random effects.
samples.eta <- model.car$samples$phi
post.median.eta <- as.numeric(apply(samples.eta,2,median))

# Q6. Create a chloropleth map of the posterior median of the spatial random effects. Use as many 
# color you retain appropriate. Use a "divergent" color palette to more easily distinguish negative from 
# positive spatial random effects. What does a positive - respectively, a negative - spatial random effect
# means here? Where are the positive spatial random effects located? And what about the negative spatial 
# random effects?
plotvar <- post.median.eta
nclr <- 
  plotclr <- brewer.pal(nclr,)
class <- classIntervals(plotvar,nclr,style="jenks")
colcode <- findColours(class,plotclr)

plot(mich.poly,border="black",axes=TRUE,xlim=c(-90,-82))
title(xlab="",ylab="",main="")
plot(mich.poly,col=colcode,add=T)

leg.txt<-c(..)
legend(-90.3,44,legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")





# Q7. Create a smoothed map for the burglary rate per 1,000 people in Michigan in 2006. Use as 
# many color classes as you think it's appropriate.
# Compare it with the map for the raw data. What do you observe are the diferences between the two maps?
plotvar <- model.car$fitted.values
nclr <- 
  plotclr <- brewer.pal(nclr,"Purples")
class <- classIntervals(plotvar,nclr,style="jenks")
colcode <- findColours(class,plotclr)

plot(mich.poly,border="black",axes=TRUE,xlim=c(-92,-81))
title(xlab="",ylab="",main="")
plot(mich.poly,col=colcode,add=T)

leg.txt<-c(paste("(",round(class$brks[1],3),"; ",round(class$brks[2],3),")",sep=""),paste("[",round(class$brks[2],3),"; ",round(class$brks[3],3),")",sep=""),...)
legend(-92,44,legend=leg.txt,fill=plotclr,cex=1,ncol=1,bty="n")

```


