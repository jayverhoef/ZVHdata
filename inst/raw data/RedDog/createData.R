library(slm)

path1 = '/media/jay/Hitachi2GB/00NMML/ActiveRPack/slm_package'
path1 = '/media/jay/data/desktop_data/2019_packages/slm_package'
path2 = '/slm/inst/raw data/wetSulfateDep/so4.txt'
path2 = '/slm/inst/raw data/RedDog/reddog.csv'
d1 = read.table(paste0(path1,path2), sep = ',', header = TRUE)
coords = d1[,c('lon','lat')]
crs = sp::CRS("+proj=longlat +datum=WGS84")
d2 = sp::SpatialPointsDataFrame(coords = coords,
	data = d1, proj4string = crs)

# USA Map
path2 = '/slm/inst/raw data/wetSulfateDep/states_21basic/'
ShapeFile = 'states'
shape.path.filename <- paste0(path1, path2, ShapeFile,'.shp')
states <- rgdal::readOGR(shape.path.filename)
# get rid of Hawaii and Alaska
Alaska = states[states$STATE_NAME == 'Alaska',]
# reproject into lat/long decimal degrees
Alaska = sp::spTransform(Alaska, sp::CRS("+proj=longlat +datum=WGS84"))
# reproject into Albers equal area conic
Alaska = sp::spTransform(Alaska, sp::CRS("+init=epsg:3338"))

# CAKR Map
path2 = '/slm/inst/raw data/RedDog/'
ShapeFile = 'National_Park_Service__Park_Unit_Boundaries'
shape.path.filename <- paste0(path1, path2, ShapeFile,'.shp')
NPSbounds <- rgdal::readOGR(shape.path.filename)
# get rid of Hawaii and Alaska
CAKR = NPSbounds[NPSbounds@data$UNIT_NAME == 'Cape Krusenstern National Monument',]
# reproject into lat/long decimal degrees
CAKR = sp::spTransform(CAKR, sp::CRS("+proj=longlat +datum=WGS84"))
# reproject into Alaska Albers equal area conic
CAKR = sp::spTransform(CAKR, sp::CRS("+init=epsg:3338"))
# reproject the observations
rdobs = sp::spTransform(d2, sp::CRS("+init=epsg:3338"))

# Road Map
path2 = '/slm/inst/raw data/RedDog/'
ShapeFile = 'DMTS_CL_Only'
shape.path.filename <- paste0(path1, path2, ShapeFile,'.shp')
haulRoad <- rgdal::readOGR(shape.path.filename)
# reproject into lat/long decimal degrees
haulRoad = sp::spTransform(haulRoad, sp::CRS("+proj=longlat +datum=WGS84"))
# reproject into Alaska Albers equal area conic
haulRoad = sp::spTransform(haulRoad, sp::CRS("+init=epsg:3338"))

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

rdpreds = rbind(stratum1, stratum2, stratum3, stratum4, stratum5)

path3 = '/slm/data/'
save(Alaska, file = paste0(path1, path3, 'Alaska.rda'))
save(CAKR, file = paste0(path1, path3, 'CAKR.rda'))
save(haulRoad, file = paste0(path1, path3, 'haulRoad.rda'))
save(rdobs, file = paste0(path1, path3, 'rdobs.rda'))
save(rdpreds, file = paste0(path1, path3, 'rdpreds.rda'))


layout(matrix(c(1,1,2,3), nrow = 2, ncol = 2, byrow = TRUE))

par(mar = c(2,0,0,0))
plot(Alaska, lwd = 2)
plot(CAKR, add = TRUE, col = 'red')

lines(c(-450000, -450000, -370000, -370000, -450000),
  c(1900000, 2030000, 2030000, 1900000, 1900000), type = 'l',
  lwd = 5, col = 'red')
text(-1520000, 2200000, 'A', cex = 4)

par(mar = c(0,0,2,0))
plot(CAKR, lwd = 2)
plot(rdobs[rdobs$year == 2001,], add = TRUE, pch = 19, cex = .8, col = 'darkorchid2')
plot(rdobs[rdobs$year == 2006,], add = TRUE, pch = 19, cex = .8, col = 'darkolivegreen3')
plot(haulRoad, add = TRUE, lwd = 2, col = 'red2')

text(-450000, 2010000, 'B', cex = 4)
legend(-390000, 1990000, legend = c('2001', '2006'),
  pch = c(19, 19), col = c('darkorchid2','darkolivegreen3'),
  cex = 2)

pal = viridis(7)
par(mar = c(0,0,2,0))
plot(CAKR, lwd = 2)
plot(haulRoad, add = TRUE, lwd = 2, col = 'red2')
plot(rdpreds[rdpreds$strat==1,], add = TRUE, pch = 19, cex = .2, col=pal[1])
plot(rdpreds[rdpreds$strat==2,], add = TRUE, pch = 19, cex = .2, col=pal[2])
plot(rdpreds[rdpreds$strat==3,], add = TRUE, pch = 19, cex = .2, col=pal[3])
plot(rdpreds[rdpreds$strat==4,], add = TRUE, pch = 19, cex = .5, col=pal[6])
plot(rdpreds[rdpreds$strat==5,], add = TRUE, pch = 19, cex = 1, col=pal[7])
text(-450000, 2010000, 'C', cex = 4)
legend(-390000, 1990000, 
  legend = c('strat 1', 'strat 2', 'strat 3', 'strat 4', 'strat 5'),
  pch = c(19, 19), col = c(pal[1:3],pal[6:7]),
  cex = 2)
