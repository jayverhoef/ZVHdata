library(slm)

path1 = '/media/jay/Hitachi2GB/00NMML/ActiveRPack/slm_package'
#path1 = '/media/jay/ExtraDrive1/transfer/slm_package'
path2 = '/slm/inst/raw data/wetSulfateDep/so4.txt'
d1 = read.table(paste0(path1,path2), sep = '',
	col.names = c('nlat.deg', 'nlat.min', 'nlat.sec', 'wlong.deg', 
	'wlong.min', 'wlong.sec', 'SO4'))
d1$y = d1$nlat.deg + (60*d1$nlat.min + d1$nlat.sec)/60^2
d1$x = -d1$wlong.deg - (60*d1$wlong.min + d1$wlong.sec)/60^2
coords = d1[,c('x','y')]
data = data.frame(SO4 = d1[,'SO4'])
crs = sp::CRS("+proj=longlat +datum=WGS84")
d2 = sp::SpatialPointsDataFrame(coords = coords,
	data = data, proj4string = crs)

# USA Map
path1 = '/media/jay/data/desktop_data/2019_packages/ZVHdata_package/ZVHdata'
path2 = '/inst/raw data/wetSulfateDep/states_21basic/'
ShapeFile = 'states'
shape.path.filename <- paste0(path1, path2, ShapeFile,'.shp')
states <- rgdal::readOGR(shape.path.filename)
# get rid of Hawaii and Alaska
states = states[!states$STATE_NAME %in% c('Alaska','Hawaii'),]
# reproject into lat/long decimal degrees
states = sp::spTransform(states, sp::CRS("+proj=longlat +datum=WGS84"))
# reproject into Albers equal area conic
states = sp::spTransform(states, sp::CRS("+init=epsg:5070"))
# reproject the observations
d2 = sp::spTransform(d2, sp::CRS("+init=epsg:5070"))

# create a systematic grid of prediction points
# get ratio of bounding box for x- and y-coordinates
xyratio = (max(sp::coordinates(d2)[,1]) - min(sp::coordinates(d2)[,1]))/
	(max(sp::coordinates(d2)[,2]) - min(sp::coordinates(d2)[,2]))
# create systematic grid
preds = pointSimSyst(nrow = 50, ncol = round(50*xyratio), 
	lower.x.lim = d2@bbox[1,1], upper.x.lim = d2@bbox[1,2],
	lower.y.lim = d2@bbox[2,1], upper.y.lim = d2@bbox[2,2])
# turn them into a SpatialPoints object
coords = preds
crs = sp::CRS("+init=epsg:5070")
preds = sp::SpatialPoints(coords = coords,
	proj4string = crs)
# clip the grid to just those points falling within US boundaries
preds_sub = preds[states,]

# sp::plot(states)
# sp::plot(preds_sub, add = TRUE, pch = 19, cex = .5)
# sp::plot(d2, add = TRUE, pch = 19, col = 'red')

USboundary = states
SO4obs = d2
SO4pred = preds_sub
path3 = '/data/'
save(USboundary, file = paste0(path1, path3, 'USboundary.rda'))
save(SO4obs, file = paste0(path1, path3, 'SO4obs.rda'))
save(SO4pred, file = paste0(path1, path3, 'SO4pred.rda'))
