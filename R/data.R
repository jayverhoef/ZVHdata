#' A boundary file for Alaska, USA
#'
#' A boundary file for Alaska, USA, of object type \code{SpatialPolygons} from \code{sp} package.
#'
#' @docType data
#' @usage data(Alaska)
#' @format A \code{SpatialPolygons} object with slots @polygons, @plotOrder, @bbox, and @proj4string. The projection is in Alaska Alber, EPSG:3338.
#' @examples
#' data(Alaska)
#' summary(Alaska)
#' library(sp)
#' plot(Alaska)
"Alaska"

#' A boundary file for Cape Krusenstern National Moment, Alaska, USA
#'
#' A boundary file for Cape Krusenstern National Moment, Alaska, USA, with object type \code{SpatialPolygonDataFrame} from \code{sp} package.
#'
#' @docType data
#' @usage data(CAKR)
#' @format A \code{SpatialPolygonDataFrame} object with slots @polygons, @plotOrder, @bbox, and @proj4string. The projection is in Alaska Albers, EPSG:3338.
#' @source Obtained from \href{https://catalog.data.gov/dataset/national-parks}{data.gov}. The National Parks dataset September 30, 2019 is part of the U.S. Department of Transportation (USDOT)/Bureau of Transportation Statistics's (BTS's) National Transportation Atlas Database (NTAD). These are the National Park Service unit boundaries. These park boundaries signify legislative boundary definitions and local park names have been consolidated according to the legislation.  The Cape Krusenstern National Monument was subsetted from the all NPS boundaries by importing with the \code{rgdal} package as a \code{SpatialPolygonDataFrame} object from package \code{sp}.
#' @examples
#' data(CAKR)
#' summary(CAKR)
#' library(sp)
#' plot(CAKR)
"CAKR"

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
#' @references Lenart, E.A., Bowyer, R.T., Ver Hoef, J.M. and Ruess, R.W. 2002. Climate Change and Caribou: Effects of Summer Weather on Forage. Canadian Journal of Zoology 80: 664-678.
#' (\href{https://www.nrcresearchpress.com/doi/abs/10.1139/z02-034#.XvweB3VKjmF}{link to: Climate Change and Caribou: Effects of Summer Weather on Forage})
#' @examples
#' data(caribouDF)
#' names(caribouDF)
#' summary(caribouDF)
#' summary(lm(z ~ water + tarp + water:tarp, data = caribouDF))
"caribouDF"

#' A line file for the  haul road to Red Dog Mine, Alaska
#'
#' A line file for the  haul road to Red Dog Mine, Alaska, USA, with object type \code{SpatialLinesDataFrame} from \code{sp} package.
#'
#' @docType data
#' @usage data(haulRoad)
#' @format A \code{SpatialLines} object with slots @lines, @bbox, and @proj4string. The projection is in Alaska Albers, EPSG:3338.
#' @examples
#' data(CAKR)
#' data(haulRoad)
#' library(sp)
#' plot(CAKR)
#' plot(haulRoad, add = TRUE, lwd = 2, col = 'red')
"haulRoad"

#' Data on heavy metals in lichens near a mining road in Alaska, USA
#'
#' Data on heavy metal concentrations in mosses near a mining road in Alaska, USA with object type \code{SpatialPointsDataFrame} from \code{sp} package. Data on three metals: 1) Cd (cadmium), 2) Pb (lead), and 3) Zn (zinc) are included, along with covariates that have relationships to heavy metal concentrations.
#'
#' @docType data
#' @usage data(rdobs)
#' @format A \code{SpatialPointsDataFrame} object with slots @data, @coords.nrs, @coords, @bbox, and @proj4string. The projection is in Alaska Albers, EPSG:3338.  The @data is a dataframe with 365 observations and 15 variables.
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
#' }
#' @references Neitlich, P.N., Ver Hoef, J.M., Berryman, S. D., Mines, A., Geiser, L.H., Hasselbach, L.M., and Shiel, A. E. 2017. Trends in Spatial Patterns of Heavy Metal Deposition on National Park Service Lands Along the Red Dog Mine Haul Road, Alaska, 2001–2006. PLOS ONE 12(5):e0177936 DOI:10.1371/journal.pone.0177936
#' (\href{https://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0177936&type=printable}{link to: Trends in Spatial Patterns of Heavy Metal Deposition on National Park Service Lands Along the Red Dog Mine Haul Road, Alaska, 2001–2006.}) \cr Hasselbach, L., Ver Hoef, J.M., Ford, J., Neitlich, P., Berryman, S., Wolk B. and Bohle, T. 2005. Spatial Patterns of Cadmium, Lead and Zinc Deposition on National Park Service Lands in the Vicinity of Red Dog Mine, Alaska. Science of the Total Environment 348: 211–230.
#' (\href{https://www.sciencedirect.com/science/article/pii/S0048969705000082}{link to: Spatial Patterns of Cadmium, Lead and Zinc Deposition on National Park Service Lands in the Vicinity of Red Dog Mine, Alaska.})
#' @examples
#' data(rdobs)
#' names(rdobs@data)
#' summary(rdobs)
#' library(sp)
#' data(CAKR)
#' plot(CAKR)
#' plot(rdobs, add = TRUE, pch = 19)
"rdobs"

#' Prediction sites for heavy metals in mosses near a mining road in Alaska, USA
#'
#' Prediction sites for heavy metal concentrations in mosses near a mining road in Alaska, USA with object type \code{SpatialPointsDataFrame} from \code{sp} package. Data on covariates used to model observed data are included with site location. 
#'
#' @docType data
#' @usage data(rdpreds)
#' @format A \code{SpatialPointsDataFrame} object with slots @data, @coords.nrs, @coords, @bbox, and @proj4string. The projection is in Alaska Albers, EPSG:3338.  The @data is a dataframe with 2357 observations and 4 variables.
#' \describe{
#'   \item{lon}{Longitude, in decimal degrees}
#'   \item{lat}{Latitude, in decimal degrees}
#'   \item{dist2road}{Distance to haul road, in meters}
#'   \item{sideroad}{A factor variable with two levels: "N" for north of the haul road, and "S" for south of the haul road}
#'   \item{strat}{A stratification variable for prediction sites, in distance classes from the haul road, where 1 is the closest, and 5 is the farthest, from the road}
#' }
#' @references Neitlich, P.N., Ver Hoef, J.M., Berryman, S. D., Mines, A., Geiser, L.H., Hasselbach, L.M., and Shiel, A. E. 2017. Trends in Spatial Patterns of Heavy Metal Deposition on National Park Service Lands Along the Red Dog Mine Haul Road, Alaska, 2001–2006. PLOS ONE 12(5):e0177936 DOI:10.1371/journal.pone.0177936
#' (\href{https://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0177936&type=printable}{link to: Trends in Spatial Patterns of Heavy Metal Deposition on National Park Service Lands Along the Red Dog Mine Haul Road, Alaska, 2001–2006.}) \cr Hasselbach, L., Ver Hoef, J.M., Ford, J., Neitlich, P., Berryman, S., Wolk B. and Bohle, T. 2005. Spatial Patterns of Cadmium, Lead and Zinc Deposition on National Park Service Lands in the Vicinity of Red Dog Mine, Alaska. Science of the Total Environment 348: 211–230.
#' (\href{https://www.sciencedirect.com/science/article/pii/S0048969705000082}{link to: Spatial Patterns of Cadmium, Lead and Zinc Deposition on National Park Service Lands in the Vicinity of Red Dog Mine, Alaska.})
#' @examples
#' data(rdpreds)
#' names(rdpreds@data)
#' summary(rdpreds)
#' library(sp)
#' data(CAKR)
#' plot(CAKR)
#' plot(rdpreds[rdpreds$strat==5,], add = TRUE, pch = 19, cex = .5, col = 'green')
#' plot(rdpreds[rdpreds$strat==2,], add = TRUE, pch = 19, cex = .2, col = 'blue')
#' plot(rdpreds[rdpreds$strat==1,], add = TRUE, pch = 19, cex = .2, col = 'red')
"rdpreds"

#' A polygon file with data on harbor seal trends in southeast Alaska, USA.
#'
#' A polygon file for southeast Alaska, USA, of object type \code{SpatialPolygonsDataFrame} from \code{sp} package. The data.frame contains 306 polygons with estimated trends of harbor seals for each polygon, and 157 polygons with missing data.
#'
#' @docType data
#' @usage data(sealPolys)
#' @format A \code{SpatialPolygonsDataFrame} object with slots @data, @polygons, @plotOrder, @bbox, and @proj4string. The projection is in Alaska Alber, EPSG:3338.  The data.frame has 463 records and 7 variables:
#' \describe{
#'   \item{polyid}{A factor variable with 1307 levels for all harbor seal polygons in Alaska. Only 463 of these occur in southeast Alaska.  Each polygon has a unique polyid.}
#'   \item{stockid}{An integer from 8 to 12. There are 5 distinct genetic stocks for harbor seals in southeast Alaska, and each polygon is associated with one of the 5 stocks. Numbers 1 to 7, and 13, are stocks elsewhere in Alaska}
#'   \item{stockname}{A factor variablewith 13 levels, which are the names of the stocks.}
#'   \item{Estimate}{The estimated trend of harbor seal abundance in each polygon, which was a linear estimate on the log scale using Poisson regression. We treat this as the response variable. 157 polygons have missing data on trends.}
#'   \item{StdErr}{The estimated standard error of the estimated trend (on the log scale).  We ignore this for demonstration purposes.}
#'   \item{x}{x-coordinate of the polygon centroid.}
#'   \item{y}{y-coordinate of the polygon centroid.}
#' }
#' @references Ver Hoef, J.M., Peterson, E. E., Hooten, M. B., Hanks, E. M., and Fortin, M.-J. 2018. Spatial Autoregressive Models for Statistical Inference from Ecological Data. Ecological Monographs, 88: 36–59. DOI: 10.1002/ecm.1283
#' (\href{https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecm.1283}{link to: Spatial Autoregressive Models for Statistical Inference from Ecological Data.})
#' @examples
#' data(sealPolys)
#' names(sealPolys)
#' summary(sealPolys)
#' library(sp)
#' plot(sealPolys)
#' boxplot(Estimate ~ as.factor(stockid), data = sealPolys@data)
#' summary(lm(Estimate ~ as.factor(stockid), data = sealPolys@data))
"sealPolys"
