library(rgdal)
library(sp)

# path to directory where package is being developed
path1 = '/media/jay/data/desktop_data/2019_packages/ZVHdata_package/ZVHdata/'
path2 = 'inst/raw data/mossHeavyMetals/'
# Moss heavy metal data
dataFile = 'MOSSobs.csv'
# read in the CSV file
MOSSobs = read.table(paste0(path1,path2,dataFile), sep = ',', header = TRUE)
# keep only relevant columns
MOSSobs = MOSSobs[,c(1:10,14)]
# get the coordinates to create SpatialPointsDataFrame
coords = MOSSobs[,c('lon','lat')]
# create proj4 projection as lat/lon
crs = sp::CRS("+proj=longlat +datum=WGS84")
# create SpatialPointsDataFrame
MOSSobs = sp::SpatialPointsDataFrame(coords = coords,
	data = MOSSobs, proj4string = crs)
# reproject into Alaska Albers
MOSSobs = sp::spTransform(MOSSobs, sp::CRS("+init=epsg:3338"))

# write object to disk as a data object for the package
path3 = 'data/'
save(MOSSobs, file = paste0(path1, path3, 'MOSSobs.rda'))
