# HEADER --------------------------------------------
#
# Project Name: EPID 595 Week 6 code
#
# Script Name: Point pattern analysis
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

### Load in some libraries, if you don't have them, please install them.

library(tidyverse)
#library(GISTools)
library(sp)
library(rgeos)
library(tmap)
library(tmaptools)
# load the spatial libraries
library(raster)
library(adehabitatHR)

## Do some stuff to make our lives easier
tmap_options(check.and.fix = TRUE)


## Kernel density estimation using R

## load in the newhaven database
load("Data/newhaven.RData")

# look at it
# select 'view' mode
tmap_mode('view')

# Create the map of blocks and incidents
tm_shape(blocks) + tm_borders() +
  tm_shape(breach) +
  tm_dots(col = 'navyblue')

# Function to choose bandwidth according Bowman and Azzalini / Scott's rule
# for use with <smooth_map> in <tmaptools>

choose_bw <- function(spdf) {
  X <- coordinates(spdf)
  sigma <-
    c(sd(X[, 1]), sd(X[, 2]))  * (2 / (3 * nrow(X))) ^ (1 / 6)
  return(sigma / 1000)
  
}

# library(tmaptools)
# tmap_mode('view')
# breach_dens <-
#   oldtmaptools::smooth_map(breach, cover = blocks, bandwidth = choose_bw(breach))
#

# runs the kernel density estimation,look up the function parameters for more options
kde.output <- kernelUD(breach, h = "href", grid = 1000)

# converts to raster
kde <- raster(kde.output)
# sets projection to British National Grid
projection(kde) <- CRS("+init=EPSG:27700")

# maps the raster in tmap, "ud" is the density variable
tm_shape(kde) + tm_raster("ud")


### Hexagonal Binning

install.packages("fMultivar", depend = TRUE)

### Here we have a custom function to bin our observations into hexagons.
hexbin_map <- function(spdf, ...) {
  hbins <- fMultivar::hexBinning(coordinates(spdf), ...)
  
  # Hex binning code block
  # Set up the hexagons to plot,  as polygons
  u <- c(1, 0,-1,-1, 0, 1)
  u <- u * min(diff(unique(sort(hbins$x))))
  v <- c(1, 2, 1, -1, -2, -1)
  v <- v * min(diff(unique(sort(hbins$y)))) / 3
  
  # Construct each polygon in the sp model
  hexes_list <- vector(length(hbins$x), mode = 'list')
  for (i in 1:length(hbins$x)) {
    pol <- Polygon(cbind(u + hbins$x[i], v + hbins$y[i]), hole = FALSE)
    hexes_list[[i]] <- Polygons(list(pol), i)
  }
  
  # Build the spatial polygons data frame
  hex_cover_sp <-
    SpatialPolygons(hexes_list, proj4string = CRS(proj4string(spdf)))
  hex_cover <- SpatialPolygonsDataFrame(hex_cover_sp,
                                        data.frame(z = hbins$z), match.ID =
                                          FALSE)
  # Return the result
  return(hex_cover)
}

## Let's create the hexagon binned map.
tmap_mode('view')
breach_hex <- hexbin_map(breach, bins = 20)
tm_shape(breach_hex) +
  tm_fill(col = 'z', title = 'Count', alpha = 0.7)

## Making a hexagon binned map of burglaries.

## Yet another custom function to help us create our map
hexprop_map <- function(spdf, ...) {
  hbins <- fMultivar::hexBinning(coordinates(spdf), ...)
  
  # Hex binning code block
  # Set up the hexagons to plot,  as polygons
  u <- c(1, 0,-1,-1, 0, 1)
  u <- u * min(diff(unique(sort(hbins$x))))
  v <- c(1, 2, 1, -1, -2, -1)
  v <- v * min(diff(unique(sort(hbins$y)))) / 3
  
  scaler <- sqrt(hbins$z / max(hbins$z))
  # Construct each polygon in the sp model
  hexes_list <- vector(length(hbins$x), mode = 'list')
  for (i in 1:length(hbins$x)) {
    pol <-
      Polygon(cbind(u * scaler[i] + hbins$x[i], v * scaler[i] + hbins$y[i]), hole =
                FALSE)
    hexes_list[[i]] <- Polygons(list(pol), i)
  }
  
  # Build the spatial polygons data frame
  hex_cover_sp <-
    SpatialPolygons(hexes_list, proj4string = CRS(proj4string(spdf)))
  hex_cover <- SpatialPolygonsDataFrame(hex_cover_sp,
                                        data.frame(z = hbins$z), match.ID =
                                          FALSE)
  # Return the result
  return(hex_cover)
}

### Make the map.
tmap_mode('plot')

breach_prop <- hexprop_map(breach, bins = 20)

tm_shape(blocks) + tm_borders(col = 'grey') +
  tm_shape(breach_prop) +
  tm_fill(col = 'indianred', alpha = 0.7) +
  tm_layout("Breach of Peace Incidents", title.position = c('right', 'top'))


