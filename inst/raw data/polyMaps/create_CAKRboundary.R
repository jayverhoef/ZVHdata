library(sp)
library(rgdal)

# download site
# https://catalog.data.gov/dataset/national-parks
# downloaded 10 Nov 2019

# path to directory where package is being developed
path1 = '/media/jay/data/desktop_data/2019_packages/ZVHdata_package/ZVHdata/'
path2 = 'inst/raw data/polyMaps/NPS_boundaries/'
# boundary shapefile
ShapeFile = 'National_Park_Service__Park_Unit_Boundaries'
# path to shapefile
shape.path.filename <- paste0(path1, path2, ShapeFile,'.shp')
# read data, automatically creates SpatialPolygonsDataFrame
CAKRboundary <- rgdal::readOGR(shape.path.filename)
# Get just Cape Krusenstern National Monument
CAKRboundary = CAKRboundary[CAKRboundary@data$UNIT_NAME == 
  'Cape Krusenstern National Monument',]
# reproject into lat/long decimal degrees
CAKRboundary = sp::spTransform(CAKRboundary, sp::CRS("+proj=longlat +datum=WGS84"))
# reproject into Alaska Albers equal area conic
CAKRboundary = sp::spTransform(CAKRboundary, sp::CRS("+init=epsg:3338"))
# simplify to SpatialPolygons, no need for data.frame information
Polys = CAKRboundary@polygons
CAKRboundary = SpatialPolygons(Polys, proj4string=CAKRboundary@proj4string)

# write object to disk as a data object for the package
path3 = 'data/'
save(CAKRboundary, file = paste0(path1, path3, 'CAKRboundary.rda'))
