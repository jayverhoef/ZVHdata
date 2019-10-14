library(slm)
library(classInt)
library(viridis)
library(colorspace)
	
# ----------------------------------------------------------------------
#                 Simulated Data
# ----------------------------------------------------------------------

	
xy_obs = pointSimCSR(npair = 100)
xy_preds = pointSimSyst(nrow = 30, ncol = 30)
plot(xy_obs, pch = 19)
points(xy_preds)

d1 = rbind(xy_obs, xy_preds)

d1$re1 = as.factor(trunc(10*runif(dim(d1)[1]))+1)
d1$re2 = as.factor(trunc(4*runif(dim(d1)[1]))+1)
d1$x1 = rnorm(dim(d1)[1])
d1$x2 = rnorm(dim(d1)[1])
d1$f1 = as.factor(trunc(2*runif(dim(d1)[1]))+1)
d1$f2 = as.factor(trunc(4*runif(dim(d1)[1]))+1)


Z1 = model.matrix( ~ -1 + re1, d1)
Z2 = model.matrix( ~ -1 + re2, d1)


#source('/media/jay/ExtraDrive1/transfer/slm_package/slm/R/simulate_geostat.R')
#undebug(simulate_geostat)
geo_error_sim = simulate_geostat(loc.data = d1, 
	parsil = 1, range = 1, nugget = 10,
	minorp = .3, rotate = 10, extrap = 1,
	spatial_model = besselk,
	random_effects_Z = list(Z1, Z2),
	random_effects_varcomps = c(5, 10))$z
	
d1$z = 2*d1$x1 + 0*d1$x2 + (as.numeric(d1$f1)-1) + 
	0.333*(as.numeric(d1$f2)-1) + geo_error_sim
 
cid1 = classIntervals(d1[,'z'], n = 10, style = 'fisher')
palcid1 = viridis(10)
cid1_colors = findColours(cid1, palcid1)

plot(d1[,c('x','y')], col = cid1_colors, pch = 19)

sim_obs = d1[1:dim(xy_obs)[1],]

summary(lm(z ~ x1 + x2 + f1 + f2, data = sim_obs))

fp = c(.01, 10, .5)
attr(fp, 'label') = c('nugget','rotate', 'extrap')
source('/media/jay/ExtraDrive1/transfer/slm_package/slm/R/cope.R')
source('/media/jay/ExtraDrive1/transfer/slm_package/slm/R/utils.R')
#source('/media/jay/Hitachi2GB/00NMML/ActiveRPack/slm_package/slm/R/summary.slm_cope.R')
#undebug(cope)
cope_out1 = cope(z ~ x1 + x2 + f1 + f2, data = sim_obs, 
		x_column = 'x', y_column = 'y', spatial_model = besselk, 
		estMeth = 'REML', max_range = 10, use_anisotropy = TRUE,
		random_formula = list(raneff1 = ~ -1 + re1, raneff2 = ~ -1 + re2),
		profile_sigma = FALSE, use_grid = TRUE)
cope_out1
cope_out1$m2LL
cope_out1$opt_time

out = cope_out1
-2*log_lik(out$theta, z = out$z, X = out$X, 
	xcoords = out$coord[,1], ycoords = out$coords[,2], 
	spatial_model = out$spatial_model,
	estMeth = out$cope_est_meth, 
	attr_theta = attributes(out$theta), 
	grpindx = out$gi, Zs = out$Zs)

cope_out1a = cope(z ~ x1 + x2 + f1 + f2, data = sim_obs, 
		x_column = 'x', y_column = 'y', spatial_model = besselk, 
		estMeth = 'REML', max_range = 10, use_anisotropy = TRUE,
		random_formula = list(raneff1 = ~ -1 + re1, raneff2 = ~ -1 + re2),
		profile_sigma = TRUE) #, use_grid = TRUE)
cope_out1a
cope_out1a$m2LL
cope_out1a$opt_time

out = cope_out1a
-2*log_lik(out$theta, z = out$z, X = out$X, 
	xcoords = out$coord[,1], ycoords = out$coords[,2], 
	spatial_model = out$spatial_model,
	estMeth = out$cope_est_meth, 
	attr_theta = attributes(out$theta), 
	grpindx = out$gi, Zs = out$Zs)

fp = cope_out1$theta[attr(cope_out1$theta,'label') == 'extrap']
attr(fp, 'label') = c('extrap')
cope_out2 = cope(z ~ x1 + x2 + f1 + f2, data = sim_obs, 
		x_column = 'x', y_column = 'y', spatial_model = besselk, 
		estMeth = 'REML', max_range = 10, use_anisotropy = TRUE,
		random_formula = list(raneff1 = ~ -1 + re1, raneff2 = ~ -1 + re2),
		profile_sigma = TRUE, use_grid = TRUE, thetaini = cope_out1$theta,
    fixed_parms = fp)
cope_out2
cope_out2$m2LL
cope_out2$opt_time


# ----------------------------------------------------------------------
#                 Simulated Data with Covariates
# ----------------------------------------------------------------------

xy_obs = pointSimCSR(npair = 200)
xy_preds = pointSimSyst(nrow = 30, ncol = 30)
plot(xy_obs, pch = 19)
points(xy_preds)

d1 = rbind(xy_obs, xy_preds)

d1$re1 = as.factor(trunc(10*runif(dim(d1)[1]))+1)
d1$re2 = as.factor(trunc(4*runif(dim(d1)[1]))+1)
d1$x1 = rnorm(dim(d1)[1])
d1$x2 = rnorm(dim(d1)[1])
d1$f1 = as.factor(trunc(2*runif(dim(d1)[1]))+1)
d1$f2 = as.factor(trunc(4*runif(dim(d1)[1]))+1)


Z1 = model.matrix( ~ -1 + re1, d1)
Z2 = model.matrix( ~ -1 + re2, d1)


source('/media/jay/Hitachi2GB/00NMML/ActiveRPack/slm_package/slm/R/simulate_geostat.R')
#undebug(simulate_geostat)
geo_error_sim = simulate_geostat(loc.data = d1, 
	parsil = 5, range = 1, nugget = 5,
	minorp = 1, rotate = 10, extrap = NULL,
	spatial_model = spherical
)$z
	
d1$z = 
	2*d1$x1 + 0*d1$x2 + (as.numeric(d1$f1)-1) + 
		0.333*(as.numeric(d1$f2)-1) + 
		geo_error_sim
 
cid1 = classIntervals(d1_sim[,'z'], n = 10, style = 'fisher')
palcid1 = viridis(10)
cid1_colors = findColours(cid1, palcid1)

plot(d1[,c('x','y')], col = cid1_colors, pch = 19)

sim_obs = d1[1:dim(xy_obs)[1],]

summary(lm(z ~ x1 + x2 + f1 + f2, data = sim_obs))

fp = .01
attr(fp, 'label') = 'nugget'
source('/media/jay/Hitachi2GB/00NMML/ActiveRPack/slm_package/slm/R/cope.R')
source('/media/jay/Hitachi2GB/00NMML/ActiveRPack/slm_package/slm/R/utils.R')
#source('/media/jay/Hitachi2GB/00NMML/ActiveRPack/slm_package/slm/R/summary.slm_cope.R')
#undebug(cope)
cope_out1 = cope(z ~ x1 + x2 + f1 + f2, data = sim_obs, 
		x_column = 'x', y_column = 'y', spatial_model = spherical, 
		estMeth = 'REML', max_range = 10, use_anisotropy = FALSE,
		profile_sigma = FALSE, use_grid = TRUE)
cope_out1
cope_out1$opt_time

cope_out2 = cope(z ~ x1 + x2 + f1 + f2, data = sim_obs, 
		x_column = 'x', y_column = 'y', spatial_model = spherical, 
		estMeth = 'REML', max_range = 10, use_anisotropy = FALSE,
		profile_sigma = TRUE, use_grid = TRUE)
cope_out2
cope_out2$opt_time

library(geoR)

cope_out1 = cope(z ~ x1 + x2 + f1 + f2, data = sim_obs, 
		x_column = 'x', y_column = 'y', spatial_model = exponential, 
		estMeth = 'ML', max_range = 10, use_anisotropy = FALSE,
		profile_sigma = FALSE)
cope_out1
cope_out1$opt_time

cope_out2 = cope(z ~ x1 + x2 + f1 + f2, data = sim_obs, 
		x_column = 'x', y_column = 'y', spatial_model = cauchy, 
		estMeth = 'REML', max_range = 10, use_anisotropy = TRUE,
		profile_sigma = TRUE)
cope_out2
cope_out2$opt_time
-cope_out2$m2LL/2

log_lik(cope_out2$theta, z = cope_out1$z, X = cope_out1$X, 
	xcoords = cope_out1$coord[,1], ycoords = cope_out1$coords[,2], 
	spatial_model = cope_out1$spatial_model,
	estMeth = cope_out1$cope_est_meth, 
	attr_theta = attributes(cope_out1$theta), 
	grpindx = cope_out1$gi, Zs = cope_out1$Zs)

geo_data = as.geodata(sim_obs, coords.col = 1:2, data.col = 9, covar.col = 5:8)
geoR_out <- likfit(geo_data, cov.model = 'cauchy', kappa = 1, fix.kappa = FALSE,
	ini.cov.pars=c(5, 1), fix.nug = FALSE, lik.method = 'REML',
	psiA = .5, fix.psiA = FALSE, psiR = 3, fix.psiR = FALSE,
	trend = trend.spatial(~ x1 + x2 + f1 + f2, geo_data))
geoR_out$parameters.summary
geoR_out$loglik
#geoR_out$parameters.summary[row.names(geoR_out$parameters.summary) == 'psiA','values']*57.2958
#1/geoR_out$parameters.summary[row.names(geoR_out$parameters.summary) == 'psiR','values']
#geoR_out$AIC

geoR_out_theta = geoR_out$parameters.summary[
	rownames(geoR_out$parameters.summary) %in% 
	c('sigmasq', 'phi', 'tausq'), 'values' ]
geoRnames = rownames(geoR_out$parameters.summary)[
	rownames(geoR_out$parameters.summary) %in% 
	c('sigmasq', 'phi', 'tausq') ]
geoRnames[geoRnames == 'sigmasq'] = 'partial sill'
geoRnames[geoRnames == 'tausq'] = 'nugget'
geoRnames[geoRnames == 'phi'] = 'range'
attr(geoR_out_theta,'label') = geoRnames

geoR_out_theta
log_lik(geoR_out_theta, z = cope_out1$z, X = cope_out1$X, 
	xcoords = cope_out1$coord[,1], ycoords = cope_out1$coords[,2], 
	spatial_model = cope_out1$spatial_model,
	estMeth = cope_out1$cope_est_meth, 
	attr_theta = attributes(cope_out1$theta), 
	grpindx = cope_out1$gi, Zs = cope_out1$Zs)


# ----------------------------------------------------------------------
#          Simulated Data with Grouping and Repeated Measures
# ----------------------------------------------------------------------

library(slm)
library(classInt)
library(viridis)
library(colorspace)

xy_obs1 = pointSimCSR(npair = 100)
xy_obs1$locID = 1:dim(xy_obs1)[1]
xy_obs1$rep = 1
xy_obs2 = pointSimCSR(npair = 100)
xy_obs2$locID = 101:(100 + dim(xy_obs2)[1])
xy_obs2$rep = 2
#repeated measurments per location
xy_obs3 = xy_obs1
xy_obs3$rep = 1
xy_obs4 = xy_obs2
xy_obs4$rep = 2
xy_obs = rbind(xy_obs1, xy_obs2, xy_obs3, xy_obs4)
xy_obs$rep_F = as.factor(xy_obs$rep)
xy_obs$locID_F = as.factor(xy_obs$locID)

xy_preds = pointSimSyst(nrow = 30, ncol = 30)
plot(xy_obs1[,1:2], pch = 19, col = 'red', cex = 2)
points(xy_obs2[,1:2], pch = 19, col = 'green', cex = 2)
points(xy_obs3[,1:2], pch = 19, col = 'orange')
points(xy_obs4[,1:2], pch = 19, col = 'blue')
points(xy_preds)

x1_1 = rnorm(dim(xy_obs1)[1])
x1_2 = rnorm(dim(xy_obs2)[1])
x1 = c(x1_1, x1_2, x1_1, x1_2)
x2_1 = rnorm(dim(xy_obs1)[1])
x2_2 = rnorm(dim(xy_obs1)[1])
x2 = c(x2_1, x2_2, x2_1, x2_2)

Z1 = model.matrix( ~ -1 + locID_F, xy_obs)


#source('/media/jay/Hitachi2GB/00NMML/ActiveRPack/slm_package/slm/R/simulate_geostat.R')
#undebug(simulate_geostat)
geo_error_sim1 = simulate_geostat(loc.data = xy_obs1[,1:2], 
	parsil = 5, range = 1, nugget = 5,
	minorp = 1, rotate = 10, extrap = NULL,
	spatial_model = spherical
)$z
geo_error_sim2 = simulate_geostat(loc.data = xy_obs2[,1:2], 
	parsil = 5, range = 1, nugget = 5,
	minorp = 1, rotate = 10, extrap = NULL,
	spatial_model = spherical
)$z
	
geo_error_sim = c(geo_error_sim1, geo_error_sim2, 
  geo_error_sim1, geo_error_sim2)

ran_eff_sim = Z1 %*% rnorm(dim(Z1)[2],0,sqrt(2))

z = 2*x1 + 0*x2 + ran_eff_sim + geo_error_sim

d1 = cbind(xy_obs, x1, x2, z)
 
summary(lm(z ~ x1 + x2, data = d1))

fp = .01
attr(fp, 'label') = 'nugget'
source('/media/jay/ExtraDrive1/transfer/slm_package/slm/R/cope.R')
source('/media/jay/ExtraDrive1/transfer/slm_package/slm/R/utils.R')
#source('/media/jay/Hitachi2GB/00NMML/ActiveRPack/slm_package/slm/R/summary.slm_cope.R')
#undebug(cope)
cope_out1 = cope(z ~ x1 + x2, data = d1, 
		x_column = 'x', y_column = 'y', spatial_model = spherical, 
		estMeth = 'REML', max_range = 10, use_anisotropy = FALSE,
		profile_sigma = FALSE, random_formula = list(locID = ~ -1 + locID_F),
    subSampCol = 'rep')
cope_out1
cope_out1$opt_time

cope_out2 = cope(z ~ x1 + x2 + f1 + f2, data = sim_obs, 
		x_column = 'x', y_column = 'y', spatial_model = spherical, 
		estMeth = 'REML', max_range = 10, use_anisotropy = FALSE,
		profile_sigma = TRUE, use_grid = TRUE)
cope_out2
cope_out2$opt_time

library(geoR)

cope_out1 = cope(z ~ x1 + x2 + f1 + f2, data = sim_obs, 
		x_column = 'x', y_column = 'y', spatial_model = exponential, 
		estMeth = 'ML', max_range = 10, use_anisotropy = FALSE,
		profile_sigma = FALSE)
cope_out1
cope_out1$opt_time

cope_out2 = cope(z ~ x1 + x2 + f1 + f2, data = sim_obs, 
		x_column = 'x', y_column = 'y', spatial_model = cauchy, 
		estMeth = 'REML', max_range = 10, use_anisotropy = TRUE,
		profile_sigma = TRUE)
cope_out2
cope_out2$opt_time
-cope_out2$m2LL/2

log_lik(cope_out2$theta, z = cope_out1$z, X = cope_out1$X, 
	xcoords = cope_out1$coord[,1], ycoords = cope_out1$coords[,2], 
	spatial_model = cope_out1$spatial_model,
	estMeth = cope_out1$cope_est_meth, 
	attr_theta = attributes(cope_out1$theta), 
	grpindx = cope_out1$gi, Zs = cope_out1$Zs)

geo_data = as.geodata(sim_obs, coords.col = 1:2, data.col = 9, covar.col = 5:8)
geoR_out <- likfit(geo_data, cov.model = 'cauchy', kappa = 1, fix.kappa = FALSE,
	ini.cov.pars=c(5, 1), fix.nug = FALSE, lik.method = 'REML',
	psiA = .5, fix.psiA = FALSE, psiR = 3, fix.psiR = FALSE,
	trend = trend.spatial(~ x1 + x2 + f1 + f2, geo_data))
geoR_out$parameters.summary
geoR_out$loglik
#geoR_out$parameters.summary[row.names(geoR_out$parameters.summary) == 'psiA','values']*57.2958
#1/geoR_out$parameters.summary[row.names(geoR_out$parameters.summary) == 'psiR','values']
#geoR_out$AIC

geoR_out_theta = geoR_out$parameters.summary[
	rownames(geoR_out$parameters.summary) %in% 
	c('sigmasq', 'phi', 'tausq'), 'values' ]
geoRnames = rownames(geoR_out$parameters.summary)[
	rownames(geoR_out$parameters.summary) %in% 
	c('sigmasq', 'phi', 'tausq') ]
geoRnames[geoRnames == 'sigmasq'] = 'partial sill'
geoRnames[geoRnames == 'tausq'] = 'nugget'
geoRnames[geoRnames == 'phi'] = 'range'
attr(geoR_out_theta,'label') = geoRnames

geoR_out_theta
log_lik(geoR_out_theta, z = cope_out1$z, X = cope_out1$X, 
	xcoords = cope_out1$coord[,1], ycoords = cope_out1$coords[,2], 
	spatial_model = cope_out1$spatial_model,
	estMeth = cope_out1$cope_est_meth, 
	attr_theta = attributes(cope_out1$theta), 
	grpindx = cope_out1$gi, Zs = cope_out1$Zs)
