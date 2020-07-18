library(sp)
library(rgdal)

# download site
# https://www.arcgis.com/home/item.html?id=b07a9393ecbd430795a6f6218443dccc
# downloaded 10 Nov 2019

# path to directory where package is being developed
path1 = '/media/jay/data/desktop_data/2019_packages/ZVHdata_package/ZVHdata/'
path2 = '/inst/raw data/polyMaps/states_21basic/'
# boundary shapefile
ShapeFile = 'states'
# path to shapefile
shape.path.filename <- paste0(path1, path2, ShapeFile,'.shp')
# read data, automatically creates SpatialPolygonsDataFrame
AKboundary <- rgdal::readOGR(shape.path.filename)
# get rid of jurisdictions not in continental US
AKboundary = AKboundary[AKboundary$STATE_NAME == 'Alaska',]
# reproject into lat/long decimal degrees
AKboundary = sp::spTransform(AKboundary, sp::CRS("+proj=longlat +datum=WGS84"))
# reproject into Alaska Albers
AKboundary = sp::spTransform(AKboundary, sp::CRS("+init=epsg:3338"))
# simplify to SpatialPolygons, no need for data.frame information
Polys = AKboundary@polygons
AKboundary = SpatialPolygons(Polys, proj4string=AKboundary@proj4string)

# write object to disk as a data object for the package
path3 = 'data/'
save(AKboundary, file = paste0(path1, path3, 'AKboundary.rda'))

# test after recompiling package
library(ZVHdata)
library(sp)
data(AKboundary)
plot(AKboundary)
 
