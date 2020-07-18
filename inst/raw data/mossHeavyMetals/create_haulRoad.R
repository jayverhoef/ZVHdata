library(rgdal)
library(sp)

# path to directory where package is being developed
path1 = '/media/jay/data/desktop_data/2019_packages/ZVHdata_package/ZVHdata/'
path2 = 'inst/raw data/mossHeavyMetals/'
# Road shapefile
ShapeFile = 'DMTS_CL_Only'
shape.path.filename <- paste0(path1, path2, ShapeFile,'.shp')
# import shapefile
haulRoad <- rgdal::readOGR(shape.path.filename)
# reproject into lat/long decimal degrees
haulRoad = sp::spTransform(haulRoad, sp::CRS("+proj=longlat +datum=WGS84"))
# reproject into Alaska Albers equal area conic
haulRoad = sp::spTransform(haulRoad, sp::CRS("+init=epsg:3338"))
# simplify to SpatialLines, no need for data.frame information
Lines = haulRoad@lines
haulRoad = SpatialLines(Lines, proj4string=haulRoad@proj4string)

# write object to disk as a data object for the package
path3 = 'data/'
save(haulRoad, file = paste0(path1, path3, 'haulRoad.rda'))
