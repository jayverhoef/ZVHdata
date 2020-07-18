library(ZVHdata)
library(sp)
data(USboundary)

# path to directory where package is being developed
path1 = '/media/jay/data/desktop_data/2019_packages/ZVHdata_package/ZVHdata/'
path2 = 'inst/raw data/wetSulfateDep/'

# create a systematic grid of prediction points
# get ratio of bounding box for x- and y-coordinates
xyratio = (USboundary@bbox[1,2] - USboundary@bbox[1,1])/
	(USboundary@bbox[2,2] - USboundary@bbox[2,1])
# create systematic grid
preds = pointSimSyst(nrow = 100, ncol = round(100*xyratio), 
	lower.x.lim =USboundary@bbox[1,1], upper.x.lim = USboundary@bbox[1,2],
	lower.y.lim = USboundary@bbox[2,1], USboundary@bbox[2,2])
# turn them into a SpatialPoints object
coords = preds
crs = sp::CRS("+init=epsg:5070")
preds = sp::SpatialPoints(coords = coords,
	proj4string = crs)
# clip the grid to just those points falling within US boundaries
SO4preds = preds[USboundary,]

# write object to disk as a data object for the package
path3 = 'data/'
save(SO4preds, file = paste0(path1, path3, 'SO4preds.rda'))

# test after recompiling package
library(ZVHdata)
library(sp)
data(USboundary)
data(SO4obs)
data(SO4preds)
plot(USboundary)
plot(SO4obs, add = TRUE, pch = 19, col = 'red')
plot(SO4preds, add = TRUE, pch = 1, cex = .5)
