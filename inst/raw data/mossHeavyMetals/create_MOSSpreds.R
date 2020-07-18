library(rgdal)
library(sp)

# path to directory where package is being developed
path1 = '/media/jay/data/desktop_data/2019_packages/ZVHdata_package/ZVHdata/'
path2 = 'inst/raw data/mossHeavyMetals/'

ShapeFile = 'Stratum1_fin'
shape.path.filename <- paste0(path1, path2, ShapeFile,'.shp')
stratum1 <- rgdal::readOGR(shape.path.filename)
stratum1$strat = 1

ShapeFile = 'Stratum2_fin'
shape.path.filename <- paste0(path1, path2, ShapeFile,'.shp')
stratum2 <- rgdal::readOGR(shape.path.filename)
stratum2$strat = 2

ShapeFile = 'Stratum3_fin'
shape.path.filename <- paste0(path1, path2, ShapeFile,'.shp')
stratum3 <- rgdal::readOGR(shape.path.filename)
stratum3$strat = 3

ShapeFile = 'Stratum4_fin'
shape.path.filename <- paste0(path1, path2, ShapeFile,'.shp')
stratum4 <- rgdal::readOGR(shape.path.filename)
stratum4$strat = 4

ShapeFile = 'Stratum5_fin'
shape.path.filename <- paste0(path1, path2, ShapeFile,'.shp')
stratum5 <- rgdal::readOGR(shape.path.filename)
stratum5$strat = 5

MOSSpreds = rbind(stratum1, stratum2, stratum3, stratum4, stratum5)

# reproject into Alaska Albers
MOSSpreds = sp::spTransform(MOSSpreds, sp::CRS("+init=epsg:3338"))

# rename covariates to match the observed data
MOSSpreds$dist2road = MOSSpreds$NEAR_DIST
MOSSpreds$sideroad = MOSSpreds$Side_Rd
MOSSpreds$lat = MOSSpreds$Lat
MOSSpreds$lon = MOSSpreds$Lon

MOSSpreds = MOSSpreds[,c('lon','lat','dist2road','sideroad','strat')]

# write object to disk as a data object for the package
path3 = 'data/'
save(MOSSpreds, file = paste0(path1, path3, 'MOSSpreds.rda'))
