library(slm)
library(viridis)

# load data for graphics and analysis
# these are all sp objects
# outline of Alaska
data(Alaska)
# outline of Cape Krusenstern National Preserve
data(CAKR)
data(rdobs)
data(rdpreds)

# make maps showing study area, data, and prediction locations
layout(matrix(c(1,1,2,3), nrow = 2, ncol = 2, byrow = TRUE))

# first plot Alaska, with Cape Krusenstern as red polygon, with big red
# box around it
par(mar = c(2,0,0,0))
plot(Alaska, lwd = 2)
plot(CAKR, add = TRUE, col = 'red')
lines(c(-450000, -450000, -370000, -370000, -450000),
  c(1900000, 2030000, 2030000, 1900000, 1900000), type = 'l',
  lwd = 5, col = 'red')
text(-1520000, 2200000, 'A', cex = 4)

# show the sampling locations in 
par(mar = c(0,0,2,0))
plot(CAKR, lwd = 2)
plot(rdobs[rdobs$year == 2001,], add = TRUE, pch = 19, cex = .8, col = 'darkorchid2')
plot(rdobs[rdobs$year == 2006,], add = TRUE, pch = 19, cex = .8, col = 'darkolivegreen3')
text(-450000, 2010000, 'B', cex = 4)
legend(-390000, 1990000, legend = c('2001', '2006'),
  pch = c(19, 19), col = c('darkorchid2','darkolivegreen3'),
  cex = 2)

pal = viridis(7)
par(mar = c(0,0,2,0))
plot(CAKR, lwd = 2)
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

names(rdobs)
rdobs$logPb = log(rdobs$Pb)
rdobs$grp = as.integer(as.factor(rdobs$year))
source('/media/jay/ExtraDrive1/transfer/slm_package/slm/R/cope.R')
source('/media/jay/ExtraDrive1/transfer/slm_package/slm/R/utils.R')
#undebug(cope)
#undebug(m2LLa)
cope_out_prof = cope(logPb ~ I(log10(dist2road)) + sideroad + 
  I(log10(dist2road)):sideroad + year + year:I(log10(dist2road)) + 
  year:sideroad + year:I(log10(dist2road)):sideroad,
  rdobs, spatial_model = spherical, subSampCol = 'grp',
  random_formula = list(site = ~ -1 + site, 
    dup_within_site = ~ -1 + site:field_dup), thetaini = cope_out_nopr$theta,
  profile_sigma = TRUE)
cope_out_prof
cope_out_prof$m2LL

out = cope_out_prof
-2*log_lik(out$theta, z = out$z, X = out$X, 
	xcoords = out$coord[,1], ycoords = out$coords[,2], 
	spatial_model = out$spatial_model,
	estMeth = out$cope_est_meth, 
	attr_theta = attributes(out$theta), 
	grpindx = out$gi, Zs = out$Zs)

#undebug(cope)
cope_out_nopr = cope(logPb ~ I(log10(dist2road)) + sideroad + 
  I(log10(dist2road)):sideroad + year + year:I(log10(dist2road)) + 
  year:sideroad + year:I(log10(dist2road)):sideroad,
  rdobs, spatial_model = spherical, subSampCol = 'grp',
  random_formula = list(site = ~ -1 + site, 
    dup_within_site = ~ -1 + site:field_dup),
  profile_sigma = FALSE, thetaini = cope_out_prof$theta)
cope_out_nopr
cope_out_nopr$m2LL
