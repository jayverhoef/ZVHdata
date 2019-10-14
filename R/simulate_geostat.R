#-------------------------------------------------------------------------------
#
#           geostatSim
#
#-------------------------------------------------------------------------------

#' Simulate geostatistical data on set of given locations
#'
#' simulates geostatistical data on set of given locations
#'
#' @param loc.data data.frame with x- and y-coordinates of locations for simulated data 
#' @param xcol name of the column in loc.data with x-coordinates, default is "x"
#' @param ycol name of the column loc.data with y-coordinates, default is "y"
#' @param parsil partial sill of autocorrelation model, default = 1
#' @param range range of autocorrelation model, default = 1
#' @param nugget range of autocorrelation model, default = 0
#' @param minorp proportion of range in x direction to that of y direction for unrotated anisotropic model, default = 1
#' @param rotate rotation of anisotropic axes, default = 90
#' @param CorModel autocorrelation model, default = "Exponential".  Other possibilities are "Spherical".
#'
#' @return data.frame of three columns, the original loc.data appended with a 3rd column of simulated geostatistical data
#'
#' @author Jay Ver Hoef
#' @export
simulate_geostat <- function(loc.data, xcol = "x", ycol = "y", 
	parsil = 1, range = 1, nugget = 0,
	minorp = 1, rotate = 90, extrap = NULL,
	spatial_model = exponential,
	random_effects_Z = NULL,
	random_effects_varcomps = NULL)
{
	xcoord <- loc.data[,xcol]
	ycoord <- loc.data[,ycol]
	n <- length(xcoord)
	dismat <- distGeoAni(xcoord, ycoord, xcoord, ycoord, rotate, range, minorp)
	# compute correlation matrix for scaled distance matrix
	CovMat <- spatial_model(dismat, extrap)
	CovMat <- parsil*CovMat + diag(nugget, nrow = n, ncol = n)
	if(!is.null(random_effects_Z))
		for(i in 1:length(random_effects_Z))
			CovMat = CovMat + random_effects_varcomps[i] * 
				random_effects_Z[[i]] %*% t(random_effects_Z[[i]])
	data.frame(loc.data, z = t(chol(CovMat))%*%rnorm(n))
}

