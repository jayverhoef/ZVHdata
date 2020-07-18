library(sp)

#download site: http://nadp.slh.wisc.edu/data/NTN/ntnAllsites.aspx
#downloaded 17 July 2020

# path to directory where package is being developed
path1 = '/media/jay/data/desktop_data/2019_packages/ZVHdata_package/ZVHdata/'
path2 = 'inst/raw data/wetSulfateDep/'
# atmospheric deposition data
data1File = 'NTN-All-cydep.csv'
# site-level data
data2File = 'NTNsites.csv'

d1 = read.csv(paste0(path1, path2, data1File))
d2 = read.csv(paste0(path1, path2, data2File))
d1 = d1[d1$yr == 1987,c('siteID','SO4')]
colnames(d1) = c('siteid','SO4')
d1$siteid = as.factor(as.character(d1$siteid))
# add site-level data to the SO4 data
d1 = merge(d1, d2, by = 'siteid') 
# copy lat and lon to variables named x and y
d1$x = d1$longitude
d1$y = d1$latitude
#use locator() to elimate sites that are not in contental USA
d1 = d1[d1$x > -132.7, ]
d1 = d1[d1$y > 22.6, ]
plot(d1$x, d1$y)
# remove missing values, which are indicated by a -9 value
d1 = d1[d1$SO4 > 0, ]

coords = d1[,c('x','y')]
data = data.frame(SO4 = d1[,'SO4'])
# create proj4 projection as lat/lon
crs = sp::CRS("+proj=longlat +datum=WGS84")
# create a SpatialPointsDataFrame
d1 = sp::SpatialPointsDataFrame(coords = coords,
	data = data, proj4string = crs)
# transform the projection to CONUS Albers EPSG:5070
d1 = sp::spTransform(d1, sp::CRS("+init=epsg:5070"))
# rename object
SO4obs = d1

# write object to disk as a data object for the package
path3 = '/data/'
save(SO4obs, file = paste0(path1, path3, 'SO4obs.rda'))

# test after recompiling package
library(ZVHdata)
library(sp)
data(USboundary)
data(SO4obs)
plot(USboundary)
plot(SO4obs, add = TRUE, pch = 19, col = 'red')
