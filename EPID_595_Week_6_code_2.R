# HEADER --------------------------------------------
#
# Project Name: EPID 595 Week 6
#
# Script Name:
#
# Script Description:
#
# Author: Peter S. Larson
# Copyright (c) Peter S. Larson, 2023
# Email:  anfangen@umich.edu
#
# Date: 2023-04-18
#
#
# Notes:
#
#



if (!require("rspatial")) remotes::install_github('rspatial/rspatial')

library(rspatial)
city <- sp_data('city')
crime <- sp_data('crime')

plot(city, col='light blue')
points(crime, col='red', cex=.5, pch='+')

tb <- sort(table(crime$CATEGORY))[-1]

xy <- coordinates(crime)
dim(xy)
## [1] 2661    2
xy <- unique(xy)
dim(xy)
## [1] 1208    2
head(xy)


# mean center
mc <- apply(xy, 2, mean)
# standard distance
sd <- sqrt(sum((xy[,1] - mc[1])^2 + (xy[,2] - mc[2])^2) / nrow(xy))


plot(city, col='light blue')
points(crime, cex=.5)
points(cbind(mc[1], mc[2]), pch='*', col='red', cex=5)
# make a circle
bearing <- 1:360 * pi/180
cx <- mc[1] + sd * cos(bearing)
cy <- mc[2] + sd * sin(bearing)
circle <- cbind(cx, cy)
lines(circle, col='red', lwd=2)


##### Finding a simple density ###########################

CityArea <- raster::area(city)
dens <- nrow(xy) / CityArea


## To compute quadrat counts I first create quadrats (a RasterLayer). I get the extent for the raster from the city polygon, 
## and then assign an an arbitrary resolution of 1000. (In real life one should always try a range of resolutions, I think).

r <- raster(city)
res(r) <- 1000
r

## To find the cells that are in the city, and for easy display, I create polygons from the RasterLayer.

r <- rasterize(city, r)
plot(r)
quads <- as(r, 'SpatialPolygons')
plot(quads, add=TRUE)
points(crime, col='red', cex=.5)


## The number of events in each quadrat can be counted using the ‘rasterize’ function. 
## That function can be used to summarize the number of points within each cell, but also to compute statistics based on the ‘marks’ (attributes). 
## For example we could compute the number of different crime types) by changing the ‘fun’ argument to another function (see ?rasterize).

nc <- rasterize(coordinates(crime), r, fun='count', background=0)
plot(nc)
plot(city, add=TRUE)


ncrimes <- mask(nc, r)
plot(ncrimes)
plot(city, add=TRUE)

f <- freq(ncrimes, useNA='no')
head(f)
##      value count
## [1,]     0    48
## [2,]     1    29
## [3,]     2    24
## [4,]     3    22
## [5,]     4    19
## [6,]     5    16
plot(f, pch=20)

## Now compute average number of cases per quadrat.

# number of quadrats
quadrats <- sum(f[,2])
# number of cases
cases <- sum(f[,1] * f[,2])
mu <- cases / quadrats
mu
## [1] 9.261484

#### Distance based measures

## As we are using a planar coordinate system we can use the dist function to compute the distances between pairs of points. 
## If we were using longitude/latitude we could compute distance via spherical trigonometry functions. 
## These are available in the sp, raster, and notably the geosphere package (among others). For example, see raster::pointDistance.


d <- dist(xy)
class(d)
## [1] "dist"

dm <- as.matrix(d)
dm[1:5, 1:5]

## To get, for each point, the minimum distance to another event, we can use the ‘apply’ function. 
## Think of the rows as each point, and the columns of all other points (vice versa could also work).

dmin <- apply(dm, 1, min, na.rm=TRUE)
head(dmin)
##         1         2         3         4         5         6
## 266.07892 293.58874  47.90260 140.80688  40.06865 510.41231

mdmin <- mean(dmin)

## Do you want to know, for each point, Which point is its nearest neighbour? 
## Use the ‘which.min’ function (but note that this ignores the possibility of multiple points at the same minimum distance).

wdmin <- apply(dm, 1, which.min)

### And what are the most isolated cases? That is the furtest away from their nearest neigbor. I plot the top 25. 

plot(city)
points(crime, cex=.1)
ord <- rev(order(dmin))
far25 <- ord[1:25]
neighbors <- wdmin[far25]
points(xy[far25, ], col='blue', pch=20)
points(xy[neighbors, ], col='red')
# drawing the lines, easiest via a loop
for (i in far25) {
  lines(rbind(xy[i, ], xy[wdmin[i], ]), col='red')
}

#### The spatstat package ####

library(spatstat)

##We start with making make a Kernel Density raster. I first create a ‘ppp’ (point pattern) object, as defined in the spatstat package.

## A ppp object has the coordinates of the points and the analysis ‘window’ (study region). To assign the points locations we need to extract the coordinates from our SpatialPoints object. To set the window, we first need to to coerce our SpatialPolygons into an ‘owin’ object. We need a function from the maptools package for this coercion.

## Coerce from SpatialPolygons to an object of class “owin” (observation window)

library(maptools)
cityOwin <- as.owin(city)
class(cityOwin)
## [1] "owin"
cityOwin
## window: polygonal boundary

### Extract coordinates from SpatialPointsDataFrame:

pts <- coordinates(crime)
head(pts)
##      coords.x1 coords.x2
## [1,]   6628868   1963718
## [2,]   6632796   1964362
## [3,]   6636855   1964873
## [4,]   6626493   1964343
## [5,]   6639506   1966094
## [6,]   6640478   1961983
## enclosing rectangle: [6620591, 6654380] x [1956729.8, 1971518.9] units


## Now we can create a ‘ppp’ (point pattern) object
p <- ppp(pts[,1], pts[,2], window=cityOwin)
## Warning: 20 points were rejected as lying outside the specified window
## Warning: data contain duplicated points
class(p)
## [1] "ppp"
p
## Planar point pattern: 2641 points
## window: polygonal boundary
## enclosing rectangle: [6620591, 6654380] x [1956729.8, 1971518.9] units
## *** 20 illegal points stored in attr(,"rejects") ***
plot(p)
## Warning in plot.ppp(p): 20 illegal points also plotted

#### Having all the data well organized, it is now easy to compute Kernel Density
ds <- density(p)
class(ds)
## [1] "im"
plot(ds, main='crime density')


### Maps with the city limits and the incidence of ‘auto-theft’, ‘drunk in public’, ‘DUI’, and ‘Arson’.

par(mfrow=c(2,2), mai=c(0.25, 0.25, 0.25, 0.25))
for (offense in c("Auto Theft", "Drunk in Public", "DUI", "Arson")) {
  plot(city, col='grey')
  acrime <- crime[crime$CATEGORY == offense, ]
  points(acrime, col = "red")
  title(offense)
}

## Create a marked point pattern object (ppp) for all crimes. It is important to coerce the marks to a factor variable.
crime$fcat <- as.factor(crime$CATEGORY)
w <- as.owin(city)
xy <- coordinates(crime)
mpp <- ppp(xy[,1], xy[,2], window = w, marks=crime$fcat)
## Warning: 20 points were rejected as lying outside the specified window
## Warning: data contain duplicated points


##We can split the mpp object by category (crime)
spp <- split(mpp)
plot(spp[1:4], main='')


### Make a map of the densities of each type of crime. 
plot(density(spp[1:4]), main='')