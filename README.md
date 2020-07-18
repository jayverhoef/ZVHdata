[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)

[![minimal R version](https://img.shields.io/badge/R%3E%3D-3.6.3-6666ff.svg)](https://cran.r-project.org/) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/kotzeb0912)](https://cran.r-project.org/package=kotzeb0912) [![packageversion](https://img.shields.io/badge/Package%20version-1.0.0-orange.svg?style=flat-square)](commits/master)

[![Last-changedate](https://img.shields.io/badge/last%20change-2020--7--13-yellowgreen.svg)](/commits/master)

# ZVHdata 
## An R data package in support of the book, "Spatial Linear Models for Environmental Data." 

#### Dale. L. Zimmerman<sup>a</sup> and Jay M. Ver Hoef<sup>b</sup> 

#### <sup>a</sup>Department of Statistics and Actuarial Science, University of Iowa
#### <sup>b</sup>Alaska Fisheries Science Center, NOAA Fisheries (NMFS) 

As a scientific work, and in keeping with common scientific practicies, we kindly request that you cite our book and applicable publications if you use our work(s) or data in your publications or presentations. Additionally, we strongly encourage and welcome collaboration to promote use of these data in the proper context and scope.  The book is currently in development:

#### Zimmerman, Dale L. and Ver Hoef, Jay. M. Spatial Linear Models for Environmental Data.  *Chapman & Hall/CRC Press*.


Executive Summary
-----------------

Data sets that are used in the book.

Installation
------------

Installation of this R data package is done through the `devtools::install_github()` function or by downloading the [source package from the latest release](https://github.com/jayverhoef/ZVHdata).

```
library("devtools")
install_github("jayverhoef/ZVHdata")
```

Examine the Example Data
------------------------

There are 4 main data sets that are used throughout the book. Here is how to find out more about each of them.

### sealPolys

This is a set of areal data with seal trends in polygons, and there some polygons having missing data.  The seals are in 5 different stocks, which serves as a covariate. The help file for this example data set can found by typing

```
library(ZVHdata)
help(sealPolys)
```
in R. Use the examples in the help file to access the data. 

### caribouDF

This is a set of areal data with plant biomass in a grid of plots.  There is a 2 x 3 factorial design of treatments applied. The help file for this example data set can found by typing

```
library(ZVHdata)
help(caribouDF)
```
in R. Use the examples in the help file to access the data. 

### SO4obs and SO4preds

This is a set of point-referenced data with SO4 air quality measurements throughout the United States.  The data are in two files, one with observed data, and the other with locations to be used for predictions. The help files for these data sets can found by typing

```
library(ZVHdata)
help(SO4obs)
help(SO4preds)
```
in R. Use the examples in the help files to access the data. 

### MOSSobs and MOSSpreds

This is a set point-referenced data with concentrations of 3 different heavy metals in mosses from Cape Krusenstern National Monument, where a road from the Reddog mine traverses the Monument.  Data were collected in two years, 2001 and 2006, and mitigation for dust from road improved between the years. Covariates include year, distance from road, and side of road.   The help files for these data sets can found by typing

```
library(ZVHdata)
help(MOSSobs)
help(MOSSpreds)
```
in R. Use the examples in the help files to access the data. 


-------------
##### Disclaimer

<sub>This repository is a scientific product and is not official communication of the Alaska Fisheries Science Center, the National Oceanic and Atmospheric Administration, or the United States Department of Commerce. All AFSC Marine Mammal Laboratory (AFSC-MML) GitHub project code is provided on an ‘as is’ basis and the user assumes responsibility for its use. AFSC-MML has relinquished control of the information and no longer has responsibility to protect the integrity, confidentiality, or availability of the information. Any claims against the Department of Commerce or Department of Commerce bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by DOC or the United States Government.</sub>
