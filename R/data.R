#' A boundary file for Alaska, USA
#'
#' A boundary file for Alaska, USA, of object type \code{sf} from \code{sf} package.
#'
#' @docType data
#' @usage data(AKboundary)
#' @format A \code{sf} object. The projection is in Alaska Albers, EPSG:3338.
#' \describe{}
#' @source Shapefile downloaded from
#' \cr \cr
#' \href{https://www.arcgis.com/home/item.html?id=b07a9393ecbd430795a6f6218443dccc}{link to: Shapefile of 50 US States hosted by ArcGIS}. 
# \cr \cr
#' Imported to \code{R} with the \code{rgdal} package as a \code{SpatialPolygonsDataFrame} object from package \code{sp}, subsetted to just the Alaska polygon, and then simplified to a \code{SpatialPolygons} object.  The code that creates the data file from the raw data can be found at \code{system.file("raw data/polyMaps/create_AKboundary.R",package = "ZVHdata")}. It was subsequently converted to \code{sf} using \code{st_as_sf} function.
#' @examples
#' library(ZVHdata)
#' library(sf)
#' data(AKboundary)
#' summary(AKboundary)
#' plot(AKboundary)
"AKboundary"

#' A boundary file for Cape Krusenstern National Moment, Alaska, USA
#'
#' A boundary file for Cape Krusenstern National Moment, Alaska, USA, with object type \code{sf} from \code{sf} package.
#'
#' @docType data
#' @usage data(CAKRboundary)
#' @format A \code{sf} object. The projection is in Alaska Albers, EPSG:3338.
#' \describe{}
#' @source Data downloaded from
#' \cr \cr
#' \href{https://catalog.data.gov/dataset/national-parks}{data.gov}. 
#' \cr \cr
#' The National Parks dataset September 30, 2019 is part of the U.S. Department of Transportation (USDOT)/Bureau of Transportation Statistics's (BTS's) National Transportation Atlas Database (NTAD). These are the National Park Service unit boundaries. These park boundaries signify legislative boundary definitions and local park names have been consolidated according to the legislation.  The Cape Krusenstern National Monument was subsetted from the all NPS boundaries by importing with the \code{rgdal} package as a \code{SpatialPolygonDataFrame} object from package \code{sp}, and then simplified to a \code{SpatialPolygons} object. The R code that creates the data file from raw data can be found at \code{system.file("raw data/polyMaps/create_CAKRboundary.R",package = "ZVHdata")}. It was subsequently converted to \code{sf} using \code{st_as_sf} function.
#' @examples
#' library(ZVHdata)
#' library(sp)
#' data(CAKRboundary)
#' summary(CAKRboundary)
#' plot(CAKRboundary)
"CAKRboundary"

#' Data for a caribou forage experiment
#'
#' Data for a caribou forage experiment in Alaska, where 30 plots were placed in the field in a 5 x 6 grid and a two factor experiment, one with water added or not, and another with 3 different tarp treatments, were applied.
#'
#' @docType data
#' @usage data(caribouDF)
#' @format A \code{data.frame} with 30 records and 7 variables:
#' \describe{
#'   \item{i}{The ith row if envisioning the grid of plots as a matrix}
#'   \item{j}{The jth column if envisioning the grid of plots as a matrix}
#'   \item{water}{A factor with level "N" if no water was added, and a "Y" if water was added.}
#'   \item{tarp}{A factor with 3 levels: "clear" for a clear tarp over the plot, "none" as a control, and "shade" for a shade tarp over the plot}
#'   \item{trt}{A factor with 6 levels, treating each unique combination of water and tarp factors as a treatment level}
#'   \item{z}{The response variable, which was the percent nitrogen in prostrate willows in the 3rd clipping of the second year of the study}
#'   \item{y}{The y-coordinates of the plots in 2-D space}
#'   \item{x}{The x-coordinates of the plots in 2-D space}
#' }
#' @source These data were provided by Elizabeth Lenart of the Alaska Department of Fish and Game.  The data were used in the publication below.
#' \cr \cr
#' The R code that creates the data file from raw data can be found at \code{system.file("raw data/caribouForage/create_caribouDF.R",package = "ZVHdata")}
#' @references Lenart, E.A., Bowyer, R.T., Ver Hoef, J.M. and Ruess, R.W. 2002. Climate Change and Caribou: Effects of Summer Weather on Forage. Canadian Journal of Zoology 80: 664-678.
#' (\href{https://www.nrcresearchpress.com/doi/abs/10.1139/z02-034#.XvweB3VKjmF}{link to: Climate Change and Caribou: Effects of Summer Weather on Forage})

#' @examples
#' library(ZVHdata)
#' data(caribouDF)
#' names(caribouDF)
#' summary(caribouDF)
#' summary(lm(z ~ water + tarp + water:tarp, data = caribouDF))
"caribouDF"

#' A line file for the  haul road through Cape Krusenstern National Monument, Alaska, USA
#'
#' A line file for the  haul road through Cape Krusenstern National Monument, Alaska, USA, with object type \code{sf}.
#'
#' @docType data
#' @usage data(haulRoad)
#' @format A \code{sf} object. The projection is in Alaska Albers, EPSG:3338.
#' \describe{}
#' @source Shapefile provided by Peter Neitlich of the National Park Service. The file was developed for the publication below. 
#' \cr \cr
#' Imported to \code{R} with the \code{rgdal} package as a \code{SpatialLines} object from package \code{sp}. The code that creates the data file from the raw data can be found at \code{system.file("raw data/mossHeavyMetals/create_haulRoad.R",package = "ZVHdata")}. It was subsequently converted to \code{sf} using \code{st_as_sf} function. 
#' @references 
#' Neitlich, P.N., Ver Hoef, J.M., Berryman, S. D., Mines, A., Geiser, L.H., Hasselbach, L.M., and Shiel, A. E. 2017. Trends in Spatial Patterns of Heavy Metal Deposition on National Park Service Lands Along the Red Dog Mine Haul Road, Alaska, 2001-2006. PLOS ONE 12(5):e0177936 DOI:10.1371/journal.pone.0177936
#' (\href{https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0177936}{link to: Trends in Spatial Patterns of Heavy Metal Deposition on National Park Service Lands Along the Red Dog Mine Haul Road, Alaska, 2001-2006.})
#' @examples
#' data(CAKRboundary)
#' library(sf)
#' data(haulRoad)
#' plot(CAKRboundary)
#' plot(haulRoad, add = TRUE, lwd = 2, col = 'red')
"haulRoad"

#' Data on heavy metals in mosses near a mining road in Alaska, USA
#'
#' Data on heavy metal concentrations in mosses near a mining road in Alaska, USA with object type \code{SpatialPointsDataFrame} from \code{sp} package. Data on three metals: 1) Cd (cadmium), 2) Pb (lead), and 3) Zn (zinc) are included, along with covariates that have relationships to heavy metal concentrations.
#'
#' @docType data
#' @usage data(MOSSobs)
#' @format A \code{sf} object. The projection is in Alaska Albers, EPSG:3338.  The dataframe has 365 observations and 11 variables.
#' \describe{
#'   \item{sample}{A factor variable with location identifier.  Some samples were duplicated in the field, and replicated in the lab, so there are 318 unique spatial locations}
#'   \item{lon}{Longitude, in decimal degrees}
#'   \item{lat}{Latitude, in decimal degrees}
#'   \item{dist2road}{Distance to haul road, in meters}
#'   \item{sideroad}{A factor variable with two levels: "N" for north of the haul road, and "S" for south of the haul road}
#'   \item{field_dup}{Integer for any field duplicate samples. All locations start with 1, and any additional samples continue with integers 2, ...}
#'   \item{lab_rep}{Integer for any replicated samples for laboratory analysis, where the sample was split and a separate laboratory analysis was performed. All analyses start with 1, and any additional analyses continue with integers 2, ...}
#'   \item{Cd}{Cadmium concentration in moss tissue, in mg/kg}
#'   \item{Pb}{Lead concentration in moss tissue, in mg/kg}
#'   \item{Zn}{Zinc concentration in moss tissue, in mg/kg}
#'   \item{year}{year of data collection, either 2001 or 2006}
#' }
#' @source Data obtained from Peter Neitlich and Linda Hasselbach of the National Park Service.  Data were used in the publications listed below. 
#' \cr \cr 
#' The code that creates the data file from the raw data can be found at \code{system.file("raw data/MossHeavyMetals/create_MOSSobs.R",package = "ZVHdata")}. It was subsequently converted to \code{sf} using \code{st_as_sf} function. 
#' @references 
#' Neitlich, P.N., Ver Hoef, J.M., Berryman, S. D., Mines, A., Geiser, L.H., Hasselbach, L.M., and Shiel, A. E. 2017. Trends in Spatial Patterns of Heavy Metal Deposition on National Park Service Lands Along the Red Dog Mine Haul Road, Alaska, 2001-2006. PLOS ONE 12(5):e0177936 DOI:10.1371/journal.pone.0177936
#' (\href{https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0177936}{link to: Trends in Spatial Patterns of Heavy Metal Deposition on National Park Service Lands Along the Red Dog Mine Haul Road, Alaska, 2001-2006.}) 
#' \cr \cr 
#' \cr \cr
#' Hasselbach, L., Ver Hoef, J.M., Ford, J., Neitlich, P., Berryman, S., Wolk B. and Bohle, T. 2005. Spatial Patterns of Cadmium, Lead and Zinc Deposition on National Park Service Lands in the Vicinity of Red Dog Mine, Alaska. Science of the Total Environment 348: 211-230.
#' (\href{https://www.sciencedirect.com/science/article/pii/S0048969705000082}{link to: Spatial Patterns of Cadmium, Lead and Zinc Deposition on National Park Service Lands in the Vicinity of Red Dog Mine, Alaska.})
#' @examples
#' library(ZVHdata)
#' library(sf)
#' data(MOSSobs)
#' names(MOSSobs)
#' summary(MOSSobs)
#' data(CAKRboundary)
#' plot(CAKRboundary)
#' MOSSobs[,'lnPb'] = log(MOSSobs$Pb)
#' plot(MOSSobs[,'lnPb'], add = TRUE, pch = 19)
"MOSSobs"

#' Prediction sites for heavy metals in mosses near a mining road in Alaska, USA
#'
#' Prediction sites for heavy metal concentrations in mosses near a mining road in Alaska, USA with object type \code{SpatialPointsDataFrame} from \code{sp} package. Data on covariates used to model observed data are included with site location. 
#'
#' @docType data
#' @usage data(MOSSpreds)
#' @format A \code{sf} object. The projection is in Alaska Albers, EPSG:3338.  The dataframe has 2357 observations and 5 variables.
#' \describe{
#'   \item{lon}{Longitude, in decimal degrees}
#'   \item{lat}{Latitude, in decimal degrees}
#'   \item{dist2road}{Distance to haul road, in meters}
#'   \item{sideroad}{A factor variable with two levels: "N" for north of the haul road, and "S" for south of the haul road}
#'   \item{strat}{A stratification variable for prediction sites, in distance classes from the haul road, where 1 is the closest, and 5 is the farthest, from the road}
#' }
#' @source This prediction grid was provided by Peter Neitlich and used in the publication given below. 
#' \cr \cr
#' Imported to \code{R} with the \code{rgdal} package as a \code{SpatialPoints} object from package \code{sp}. The code that creates the data file from the raw data can be found at \code{system.file("raw data/MOSSHeavyMetals/create_MOSSpreds.R",package = "ZVHdata")}. It was subsequently converted to \code{sf} using \code{st_as_sf} function. 
#' @references
#' Neitlich, P.N., Ver Hoef, J.M., Berryman, S. D., Mines, A., Geiser, L.H., Hasselbach, L.M., and Shiel, A. E. 2017. Trends in Spatial Patterns of Heavy Metal Deposition on National Park Service Lands Along the Red Dog Mine Haul Road, Alaska, 2001-2006. PLOS ONE 12(5):e0177936 DOI:10.1371/journal.pone.0177936
#' (\href{https://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0177936&type=printable}{link to: Trends in Spatial Patterns of Heavy Metal Deposition on National Park Service Lands Along the Red Dog Mine Haul Road, Alaska, 2001-2006.}) 
#' @examples
#' library(ZVHdata)
#' library(sf)
#' data(MOSSpreds)
#' names(MOSSpreds)
#' summary(MOSSpreds)
#' data(CAKRboundary)
#' plot(CAKRboundary)
#' plot(MOSSpreds[MOSSpreds$strat==5,], add = TRUE, pch = 19, cex = .5, col = 'green')
#' plot(MOSSpreds[MOSSpreds$strat==2,], add = TRUE, pch = 19, cex = .2, col = 'blue')
#' plot(MOSSpreds[MOSSpreds$strat==1,], add = TRUE, pch = 19, cex = .2, col = 'red')
"MOSSpreds"

#' A polygon file with data on harbor seal trends in southeast Alaska, USA.
#'
#' A polygon file for southeast Alaska, USA, of object type \code{sf}. The data.frame contains 306 polygons with estimated trends of harbor seals for each polygon, and 157 polygons with missing data.
#'
#' @docType data
#' @usage data(sealPolys)
#' @format A \code{sf} object. The projection is in Alaska Alber, EPSG:3338.  The data.frame has 463 records and 7 variables:
#' \describe{
#'   \item{polyid}{A factor variable with 1307 levels for all harbor seal polygons in Alaska. Only 463 of these occur in southeast Alaska.  Each polygon has a unique polyid.}
#'   \item{stockid}{An integer from 8 to 12. There are 5 distinct genetic stocks for harbor seals in southeast Alaska, and each polygon is associated with one of the 5 stocks. Numbers 1 to 7, and 13, are stocks elsewhere in Alaska}
#'   \item{stockname}{A factor variablewith 13 levels, which are the names of the stocks.}
#'   \item{Estimate}{The estimated trend of harbor seal abundance in each polygon, which was a linear estimate on the log scale using Poisson regression. We treat this as the response variable. 157 polygons have missing data on trends.}
#'   \item{StdErr}{The estimated standard error of the estimated trend (on the log scale).  We ignore this for demonstration purposes.}
#'   \item{x}{x-coordinate of the polygon centroid.}
#'   \item{y}{y-coordinate of the polygon centroid.}
#' }
#' @source These data were collected by the Polar Ecosystem Program of the Marine Mammal Laboratory of the Alaska Fisheries Science Center of NOAA Fisheries. The data were used in the publication given below.
#' \cr \cr
#' The data and code for the analyses were permanently archived as an \code{R} package at \href{https://zenodo.org/record/1035987#.XxJRzXVKikA}{link to: Zenodo.}  Here, we simplified the data into a single \code{SpatialPolygonsDataFrame}. It was subsequently converted to \code{sf} using \code{st_as_sf} function. 
#' @references
#' Ver Hoef, J.M., Peterson, E. E., Hooten, M. B., Hanks, E. M., and Fortin, M.-J. 2018. Spatial Autoregressive Models for Statistical Inference from Ecological Data. Ecological Monographs, 88: 36-59. DOI: 10.1002/ecm.1283 \cr
#' Link to:\href{https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecm.1283}{link to: Spatial Autoregressive Models for Statistical Inference from Ecological Data.}.
#' @examples
#' library(ZVHdata)
#' library(sf)
#' data(sealPolys)
#' names(sealPolys)
#' summary(sealPolys)
#' plot(sealPolys)
#' boxplot(Estimate ~ as.factor(stockid), data = sealPolys)
#' summary(lm(Estimate ~ as.factor(stockid), data = sealPolys))
"sealPolys"

#' Data on SO4 atmospheric deposition throughout the conterminous USA
#'
#' Data on SO4 atmospheric deposition throughout the conterminous USA with object type \code{sf} from \code{sf} package.
#'
#' @docType data
#' @usage data(SO4obs)
#' @format A \code{sf} object. The projection is in Conus (Continental US) Albers, EPSG:5070.  The dataframe has 197 observations of a single variable.
#' \describe{
#'   \item{SO4}{Total wet deposition amounts of sulfate (SO4), in kilograms per hectare in 1987, at sites of the National Atmospheric Deposition Program/National Trends Network (NADP/NTN).}
#' }
#' @source 
#' These data were used in the publication listed below. Data were downloaded from the website \cr \cr
#' (\href{http://nadp.slh.wisc.edu/NTN/}{link to: National Atmospheric Deposition Program}). It was subsequently converted to \code{sf} using \code{st_as_sf} function. 
#' \cr \cr
#' The code that creates the data file from the raw data can be found at \code{system.file("raw data/wetSulfateDep/create_SO4obs.R",package = "ZVHdata")}
#' @references
#' Zimmerman, D.L. (1994). Statistical analysis of spatial data. Pages 375-402 in \emph{ Statistical Methods for Physical Science}, J. Stanford and S. Vardeman (eds.), Academic Press: New York. \cr
#' Link to:\href{https://www.elsevier.com/books/statistical-methods-for-physical-science/stanford/978-0-12-475973-2}{link to: Statistical analysis of spatial data}.
#' @examples
#' library(ZVHdata)
#' library(sf)
#' data(SO4obs)
#' names(SO4obs)
#' summary(SO4obs)
#' data(USboundary)
#' plot(st_geometry(USboundary))
#' plot(st_geometry(SO4obs), add = TRUE, pch = 19, cex = .5)
"SO4obs"

#' Prediction sites for SO4 atmospheric deposition throughout the conterminous USA
#'
#' Prediction sites for SO4 atmospheric deposition throughout the conterminous USA with object type \code{sf} from \code{sf} package.
#'
#' @docType data
#' @usage data(SO4preds)
#' @format A \code{sf} object. The projection is in Conus (Continental US) Albers, EPSG:5070. There are 9558 locations for prediction.
#' @source The code that creates the data file from a systematic grid within the conterminous US can be found at \code{system.file("raw data/wetSulfateDep/create_SO4preds.R",package = "ZVHdata")}. It was subsequently converted to \code{sf} using \code{st_as_sf} function. 
#' @examples
#' library(ZVHdata)
#' library(sf)
#' data(SO4preds)
#' str(SO4preds)
#' summary(SO4preds)
#' data(USboundary)
#' plot(st_geometry(USboundary))
#' plot(st_geometry(SO4preds), add = TRUE, pch = 19, cex = .5)
"SO4preds"

#' A boundary file for conterminous United States
#'
#' A boundary file for conterminous United States, of object type \code{SpatialPolygonsDataFrame} from \code{sp} package.
#'
#' @docType data
#' @usage data(USboundary)
#' @format A \code{sf} object. The projection is in Conus (Continental US) Albers, EPSG:5070.  There are 49 polygons with a single variable for the state name.
#' \describe{
#'   \item{STATE_NAME}{Names of the 48 contiguous states, plus District of Columbia.}
#' }
#' @source Shapefile downloaded from
#' \cr \cr 
#' \href{https://www.arcgis.com/home/item.html?id=b07a9393ecbd430795a6f6218443dccc}{link to: Shapefile of 50 US States hosted by ArcGIS}). 
#' \cr \cr
#' Imported to \code{R} with the \code{rgdal} package as a \code{SpatialPolygonsDataFrame} object from package \code{sp}, and then the contiguous US states were subsetted from that object that included other US jurisdictions. The code that creates the data file from the raw data can be found at \code{system.file("raw data/polyMaps/create_USboundary.R",package = "ZVHdata")}. It was subsequently converted to \code{sf} using \code{st_as_sf} function. 
#' @examples
#' library(ZVHdata)
#' library(sf)
#' data(USboundary)
#' summary(USboundary)
#' plot(USboundary)
"USboundary"
