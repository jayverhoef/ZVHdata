library(slm)
library(classInt)
library(viridis)
library(colorspace)

# ----------------------------------------------------------------------
#                 Simulated Data
# ----------------------------------------------------------------------


start_time = Sys.time()
nsim = 500
store_rot_est = matrix(nrow = nsim, ncol = 2)
for(kk in 1:nsim) {
	cat('\nSimulation number:', kk)
	xy_obs = pointSimCSR(npair = 100)
	xy_preds = pointSimSyst(nrow = 30, ncol = 30)
	plot(xy_obs, pch = 19)
	points(xy_preds)

	d1 = rbind(xy_obs, xy_preds)

	d1_sim = simulate_geostat(loc.data = d1, 
		parsil = 1, range = 1, nugget = 0.01,
		minorp = .3, rotate = 175, extrap = NULL,
		spatial_model = spherical)

	cid1 = classIntervals(d1_sim[,'z'], n = 10, style = 'fisher')
	palcid1 = viridis(10)
	cid1_colors = findColours(cid1, palcid1)

	plot(d1, col = cid1_colors, pch = 19)

	sim_obs = d1_sim[1:dim(xy_obs)[1],]

	source('/media/jay/Hitachi2GB/00NMML/ActiveRPack/slm_package/slm/R/cope.R')
	fp = c(0, 1, .90, 1)
	attr(fp,'label') = c('partial sill','range', 'rotate', 'minorp')
	#undebug(cope)
	cope_out = cope(z ~ 1, data = sim_obs, 
		x_column = 'x', y_column = 'y', spatial_model = exponential, 
		estMeth = 'REML', max_range = 10, use_anisotropy = TRUE,
		fixed_parms = fp)
  summary(cope_out, digits = 8)

  summary(lm(z ~ 1, data = sim_obs))
  
  store_rot_est[kk,] = c(
		cope_out$theta[attr(cope_out$theta,'label') == 'rotate'],
		cope_out$theta[attr(cope_out$theta,'label') == 'minorp'])
}
end_time = Sys.time()
end_time - start_time

store_rot_est[,1] = round(store_rot_est[,1], digits = 1)
store_rot_est[,2] = round(store_rot_est[,2], digits = 3)
store_rot_est = data.frame(rotate = store_rot_est[,1],
	minorp = store_rot_est[,2])
#store_rot_est

#hist(store_rot_est[store_rot_est$minorp < .9,1])
store_rot_est[,1][store_rot_est[,1] <= 180 & store_rot_est[,1] > 90] = 
	store_rot_est[,1][store_rot_est[,1] <= 180 & store_rot_est[,1] > 90] - 
	180
hist(store_rot_est[,1], breaks = -45:45)
mean(store_rot_est[,1])


