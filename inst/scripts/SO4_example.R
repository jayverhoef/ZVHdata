library(slm)
#library(sp)
library(classInt)
library(viridis)
library(colorspace)

# load data from package
data(SO4obs)
data(SO4pred)
data(USboundary)

# Estimate covariance parameters for constant mean model using REML
#source('/media/jay/Hitachi2GB/00NMML/ActiveRPack/slm_package/slm/R/cope.R')
#source('/media/jay/Hitachi2GB/00NMML/ActiveRPack/slm_package/slm/R/utils.R')
#undebug(cope)
#undebug(m2LLg)

cope_out = cope(SO4 ~ 1, SO4obs, spatial_model = spherical)
cope_out

cope_out_ani = cope(SO4 ~ 1, SO4obs, spatial_model = spherical,
	use_anisotropy = TRUE)
cope_out_ani


cope_outa = cope(SO4 ~ 1, SO4obs, spatial_model = besselk,
	use_grid = TRUE)
cope_outa


cope_outb = cope(SO4 ~ 1, SO4obs, spatial_model = besselk,
	fixed_parms = cope_outa$theta, use_anisotropy = TRUE)
cope_outb


cope_outc = cope(SO4 ~ 1, SO4obs, spatial_model = besselk,
	use_anisotropy = TRUE, 
	thetaini = cope_outb$theta,
	use_grid = TRUE)
cope_outc

cope(SO4 ~ 1, SO4obs, spatial_model = besselk,
	use_anisotropy = TRUE, use_grid = TRUE)

fp = c(643000, 70)
attr(fp, 'label') = c('range', 'rotate')
ti = c(.3, .1, .7)
attr(ti, 'label') = c('partial sill', 'nugget', 'minorp', 'rotate')
source('/media/jay/Hitachi2GB/00NMML/ActiveRPack/slm_package/slm/R/cope.R')
#undebug(cope)
cope_outb = cope(SO4 ~ 1, SO4obs, spatial_model = exponential,
	use_anisotropy = TRUE, use_grid = TRUE, thetaini = ti,
	fixed_parms = fp, profile_sigma = TRUE)
cope_outb

	
out = cope_outb
-2*log_lik(out$theta, z = out$z, X = out$X, 
	xcoords = out$coord[,1], ycoords = out$coords[,2], 
	spatial_model = out$spatial_model,
	estMeth = out$cope_est_meth, 
	attr_theta = attributes(out$theta), 
	grpindx = out$gi, Zs = out$Zs)

AIC(cope_outb)

# Ordinary Kriging
source('/media/jay/Hitachi2GB/00NMML/ActiveRPack/slm_package/slm/R/pulo.R')
#undebug(pulo)
pulo_out = pulo(cope_out_ani, SO4pred)

# Put the predictions into the prediction data set
SO4p = SpatialPointsDataFrame(SO4pred, as.data.frame(pulo_out))


# Make traditional plots of predictions and standard errors
cip = classIntervals(pulo_out[,1], n = 10, style = 'fisher')
palp = viridis(20)[11:20]
cip_colors = findColours(cip, palp)

cise = classIntervals(pulo_out[,2], n = 10, style = 'fisher')
palse = viridis(20)[1:10]
cise_colors = findColours(cise, palse)

layout(matrix(1:4, nrow = 2, byrow = TRUE), widths = c(3,1,3,1))
par(mar = c(0,0,5,0))
plot(SO4p, col = cip_colors, pch = 15)
plot(USboundary, add = TRUE, border = 'black')
par(mar = c(0,0,0,0))
text(-400000, 3100000, pos = 4, 'Predictions', cex = 4)
plot(c(0,1),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
	xlab = '', ylab = '', bty = 'n')
addBreakColorLegend(xleft = 0, ybottom = .1, xright = .2, ytop = .8,
	breaks = cip$brks, colors = palp, cex = 1.5)
par(mar = c(0,0,5,0))
plot(SO4p, col = cise_colors, pch = 15)
plot(USboundary, add = TRUE, border = 'white')
text(-700000, 3100000, pos = 4, 'Standard Errors', cex = 4)
plot(SO4obs, add = TRUE, col = 'white')
par(mar = c(0,0,0,0))
plot(c(0,1),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
	xlab = '', ylab = '', bty = 'n')
addBreakColorLegend(xleft = 0, ybottom = .1, xright = .2, ytop = .8,
	breaks = cise$brks, colors = palse, cex = 1.5)


# Make plots of predictions with transparency proportional to std. error
pal = viridis(10)
ci = classIntervals(pulo_out[,1], n = 10, style = 'fisher')
ci_colors = findColours(cip, pal)
cise = classIntervals(pulo_out[,2], n = 10, style = 'fisher')
ci_alpha = 1 - as.numeric(findCols(cise))/10 + .1

alp = ( (pulo_out[,2] - min(pulo_out[,2])) / 
	(max(pulo_out[,2]) - min(pulo_out[,2])) )*.8 + .2
layout(matrix(1:2, nrow = 1, byrow = TRUE), widths = c(3,1,3,1))
par(mar = c(0,0,5,0))
plot(SO4p, col = rgb(red = hex2RGB(ci_colors)@coords[,1], 
	green = hex2RGB(ci_colors)@coords[,2], 
	blue = hex2RGB(ci_colors)@coords[,3],
	alpha = ci_alpha), pch = 15)
plot(SO4obs, add = TRUE, col = 'black', pch = 19)
plot(USboundary, add = TRUE, border = 'black')
par(mar = c(0,0,0,0))
plot(c(0,1),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
	xlab = '', ylab = '', bty = 'n')
addBreakColorLegend(xleft = 0, ybottom = .3, xright = .2, ytop = .7,
	breaks = cise$brks, colors = pal, cex = 1.5)

# ----------------------------------------------------------------------
#                 Universal Kriging
# ----------------------------------------------------------------------

library(slm)
#library(sp)
library(classInt)
library(viridis)
library(colorspace)

# load data from package
data(SO4obs)
data(SO4pred)
data(USboundary)

SO4uk = SO4obs
DF = SO4uk@data
DF[,'x'] = coordinates(SO4obs)[,1]/10e+6
DF[,'y'] = coordinates(SO4obs)[,2]/10e+6
SO4uk@data = DF
SO4ukp = SpatialPointsDataFrame(SO4pred, 
	data.frame(x = coordinates(SO4pred)[,1]/10e+6, 
	y = coordinates(SO4pred)[,2]/10e+6))


#source('/media/jay/Hitachi2GB/00NMML/ActiveRPack/slm_package/slm/R/utils.R')
#undebug(m2LLa)
cope_out = cope(SO4 ~ x + y + I(x^2) + I(y^2) + I(x*y), data = SO4uk, 
	spatial_model = spherical, estMeth = 'ML')
#source('/media/jay/Hitachi2GB/00NMML/ActiveRPack/slm_package/slm/R/summary.slm_cope.R')
#undebug(summary.slm_cope)
cope_out
cope_out$opt_time

#source('/media/jay/Hitachi2GB/00NMML/ActiveRPack/slm_package/slm/R/AIC.R')
#undebug(AIC.slm_cope)
AIC(cope_out)
logLik(cope_out)

source('/media/jay/Hitachi2GB/00NMML/ActiveRPack/slm_package/slm/R/GR2.R')
#undebug(GR2)
GR2(cope_out)

#source('/media/jay/Hitachi2GB/00NMML/ActiveRPack/slm_package/slm/R/fefe.R')
#undebug(fefe)
fefe_out = fefe(cope_out)
#source('/media/jay/Hitachi2GB/00NMML/ActiveRPack/slm_package/slm/R/summary.slm_fefe.R')
#undebug(print.summary.slm_fefe)
fefe_out

source('/media/jay/Hitachi2GB/00NMML/ActiveRPack/slm_package/slm/R/pulo.R')
#undebug(pulo)
pulo_out = pulo(cope_out, SO4ukp)


# Make traditional plots of predictions and standard errors
cip = classIntervals(pulo_out[,1], n = 10, style = 'fisher')
palp = viridis(20)[11:20]
cip_colors = findColours(cip, palp)

cise = classIntervals(pulo_out[,2], n = 10, style = 'fisher')
palse = viridis(20)[1:10]
cise_colors = findColours(cise, palse)

layout(matrix(1:4, nrow = 2, byrow = TRUE), widths = c(3,1,3,1))
par(mar = c(0,0,5,0))
plot(SO4ukp, col = cip_colors, pch = 15)
plot(USboundary, add = TRUE, border = 'black')
par(mar = c(0,0,0,0))
text(-400000, 3100000, pos = 4, 'Predictions', cex = 4)
plot(c(0,1),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
	xlab = '', ylab = '', bty = 'n')
addBreakColorLegend(xleft = 0, ybottom = .1, xright = .2, ytop = .8,
	breaks = cip$brks, colors = palp, cex = 1.5)
par(mar = c(0,0,5,0))
plot(SO4ukp, col = cise_colors, pch = 15)
plot(USboundary, add = TRUE, border = 'white')
text(-700000, 3100000, pos = 4, 'Standard Errors', cex = 4)
plot(SO4obs, add = TRUE, col = 'white')
par(mar = c(0,0,0,0))
plot(c(0,1),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
	xlab = '', ylab = '', bty = 'n')
addBreakColorLegend(xleft = 0, ybottom = .1, xright = .2, ytop = .8,
	breaks = cise$brks, colors = palse, cex = 1.5)
	
