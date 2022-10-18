setwd(paste0('/mnt/ExtraDrive1/Work/desktop_data/ActiveBooks',
	'/GeostatDale/ZVHdata/ZVHdata/data/'))
library(ZVHdata)
library(sf)

load('AKboundary.rda')
AKboundary = st_as_sf(AKboundary)
save(AKboundary, file = 'AKboundary.rda')

load('CAKRboundary.rda')
CAKRboundary = st_as_sf(CAKRboundary)
save(CAKRboundary, file = 'CAKRboundary.rda')

load('haulRoad.rda')
class(haulRoad)
haulRoad = st_as_sf(haulRoad)
save(haulRoad, file = 'haulRoad.rda')

load('MOSSobs.rda')
class(MOSSobs)
MOSSobs = st_as_sf(MOSSobs)
save(MOSSobs, file = 'MOSSobs.rda')

load('MOSSpreds.rda')
class(MOSSpreds)
MOSSpreds = st_as_sf(MOSSpreds)
save(MOSSpreds, file = 'MOSSpreds.rda')

load('sealPolys.rda')
class(sealPolys)
sealPolys = st_as_sf(sealPolys)
save(sealPolys, file = 'sealPolys.rda')

load('SO4obs.rda')
class(SO4obs)
SO4obs = st_as_sf(SO4obs)
save(SO4obs, file = 'SO4obs.rda')

load('SO4preds.rda')
class(SO4preds)
SO4preds = st_as_sf(SO4preds)
save(SO4preds, file = 'SO4preds.rda')

load('USboundary.rda')
class(USboundary)
USboundary = st_as_sf(USboundary)
save(USboundary, file = 'USboundary.rda')
