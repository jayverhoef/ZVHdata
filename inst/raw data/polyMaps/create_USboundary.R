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
USboundary <- rgdal::readOGR(shape.path.filename)
# get rid of jurisdictions not in continental US
USboundary = USboundary[!USboundary$STATE_NAME %in% c('Alaska','Hawaii'),]
# reproject into lat/long decimal degrees
USboundary = sp::spTransform(USboundary, sp::CRS("+proj=longlat +datum=WGS84"))
# reproject into Albers equal area conic
USboundary = sp::spTransform(USboundary, sp::CRS("+init=epsg:5070"))
# reproject the observations in CONUS Albers EPSG:5070
# simplify @data to only include state name
DF = USboundary@data
# get rid of factor levels that were eliminated
DF[,'STATE_NAME'] = as.factor(as.character(DF[,'STATE_NAME']))
USboundary@data = DF
USboundary = USboundary[,'STATE_NAME']
path3 = '/data/'
save(USboundary, file = paste0(path1, path3, 'USboundary.rda'))

# test after recompiling package
library(ZVHdata)
library(sp)
data(USboundary)
plot(USboundary)
 
