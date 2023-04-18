# HEADER --------------------------------------------
#
# Project Name: EPID 595 
#
# Script Name: Week 6 Code set 3
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


## The K-function is a method used in spatial Point Pattern Analysis (PPA) to inspect the spatial distribution of a set of points. 
## It allows the user to assess if the set of points is more or less clustered that what we could expect from a given distribution. 
## Most of the time, the set of point is compared with a random distribution.

library(spNetwork)
library(tmap)

data(main_network_mtl)
data(mtl_libraries)
data(mtl_theatres)

tm_shape(main_network_mtl) + 
  tm_lines("black") + 
  tm_shape(mtl_libraries) + 
  tm_dots(col = "red", size = 0.2) +
  tm_shape(mtl_theatres) + 
  tm_dots(col = "blue", size = 0.2)

## As one can see, the theatres seems to be more clustered than the libraries.

kfun_theatre <- kfunctions(main_network_mtl, mtl_theatres,
                           start = 0, end = 5000, step = 50, 
                           width = 1000, nsim = 50, resolution = 50,
                           verbose = FALSE, conf_int = 0.05)
kfun_theatre$plotk

## The blue line is the empirical network K-function of the theatres in Montreal. 
## The gray area represents the results of the 50 simulations in the interval 2.5% - 97.5%. 
## Because the blue line is way above the gray area, we can conclude that the theatres are more clustered 
## than what we can expect from a random distribution. 
## (Note: usually, more simulations are required for inference).

kfun_theatre$plotg

## The G-function is also indicating a clustering situation, which is maximum between two and three kilometers. 
## This is consistent with the fact that we have a high concentration of theatres in the central neighbourhoods and then more dispersed points. 
## We can perform the same analysis for libraries.

kfun_biblio <- kfunctions(main_network_mtl, mtl_libraries,
                          start = 0, end = 5000, step = 50,
                          width = 1000, nsim = 50, verbose = FALSE)

kfun_biblio$plotk

## For distances bellow two kilometres, the libraries are a bit dispersed. 
## Above two kilometres, the libraries tend to be randomly located. 
## The chart of the G-function confirms this observation, but with a lower limit (1.3 km) :

kfun_biblio$plotg

## The network cross-K-function in spNetwork

## The cross-K-function is used to determine if two set of points A and B tend to be close or far away one from each other.

## Note that the cross-K-function A to B is not necessarily the same results as the cross-K-function B to A. 
## Again, in spNetwork, the inference is based on Monte Carlo Simulations. The locations of the reference set 
## of points are randomized to estimate if the current situation is more or less clustered that what we could expect at random.

cross_biblio_theatre <- cross_kfunctions(main_network_mtl, mtl_libraries,
                                         mtl_theatres, start = 0, end = 5000, step = 50,
                                         width = 1000, nsim = 50, verbose = FALSE)

cross_biblio_theatre$plotk

## One can conclude from the chart that the libraries are randomly 
## located around the theatres and display nor clustering nor dispersion around theatres.

cross_theatre_biblio <- cross_kfunctions(main_network_mtl, mtl_theatres,
                                         mtl_libraries, start = 0, end = 5000,
                                         step = 50, width = 1000, nsim = 50, verbose = FALSE)

cross_theatre_biblio$plotk

## However, this second chart shows that the theatres tend to be clustered around libraries. 
## This is coherent with the map above. A random library is often located far from the theatres. 
## But the theatres are concentrated in the city centre and close to some specific libraries.