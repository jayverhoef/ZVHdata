#-------------------------------------------------------------------------------
#
#           fefe
#
#-------------------------------------------------------------------------------

#' Fits geostatistical linear model
#'
#' fits a geostatistical linear model
#'
#' @param formula an R linear model formula
#' @param spdata an sp SpatialPointsDataFrame
#' @param estMeth estimation method.  Default is "REML" for restricted maximum likelihood.  Other options are "ML" for maximum likelihood
#' @param varComps a list of variance components, including spatial autocorrelation models 
#' for errors and traditional random effects.  The list of spatial autocorrelation 
#' models is "exponential","expRadon2","expRadon4","gaussian","stable",
#' "rationalQuad","cauchyGrav","cauchyMag","cauchy","circular","spherical",
#' "cubic","penta","cardinalSine","besselK","besselJ"  Default is "exponential".
#' Any names in the list not given above will be searched among the columns in the data set and
#' used as a factor variable for levels of a traditional random effect.
#' @param useAnisotropy include anistropy in parameter estimation?  Default is "FALSE"
#'
#' @return a list of class "splmm".  The functions "summary" and "print" are used to obtain and print a summary. "anova" returns just the analysis of variance table...
#'
#' @author Jay Ver Hoef
#' @export
fefe <- function(cope_out, data = NULL, spatial_model = NULL, 
	theta = NULL, z = NULL, X = NULL, xcoords = NULL, ycoords = NULL, 
	extrap = NULL, par = FALSE, subsampindx = NULL) 
{
	grpindx = subsampindx
	Vilist = NULL
	bhatList = NULL
	Sxx =  NULL
	Sxy =  NULL
	Wxx = NULL
	ViX = NULL
	Viz = NULL

	cope_est_meth = 'User'
	if(!is.null(cope_out)) {
		theta = cope_out$theta
		extrap = cope_out$extrap
		xcoords = cope_out$coords[,1]
		ycoords = cope_out$coords[,2]
		z = cope_out$z
		X = cope_out$X
		spatial_model = cope_out$spatial_model
		spatial_name = cope_out$spatial_name
		cope_est_meth = cope_out$cope_est_meth
  }
	if(is.null(subsampindx)) {
		dismat <- distGeoAni(xcoords, ycoords, xcoords, ycoords, 
			rotate = 90, range = theta[2], minorp = 1)	
		if(length(theta) == 4) extrap = theta[4]
		covMat <- theta[1]*spatial_model(dismat, extrap) +
				theta[3]*diag(dim(dismat)[1])
		QR = qr(covMat)
		ViX <- solve(covMat, X)
		Viz = solve(covMat,z)
		XViX <- crossprod(X, ViX)
		covb <- solve(XViX)
		bhat <- covb %*% crossprod(ViX, z)
	}
	if(!is.null(subsampindx)) {
		Wxxloop = function(k, ijpair, theta, X, Vilist, 
			xcoords, ycoords, grpindx)
		{ 
			t(X[grpindx == ijpair[k,1],]) %*% Vilist[[ijpair[k,1]]] %*% (theta[1]*
				corModelExponential(
					distGeoAni(
						xcoords[grpindx == ijpair[k,1]], ycoords[grpindx == ijpair[k,1]],
						xcoords[grpindx == ijpair[k,2]], ycoords[grpindx == ijpair[k,2]], 
						rotate = 90, range = theta[2], minorp = 1
					)
				)) %*% Vilist[[ijpair[k,2]]]%*% X[grpindx == ijpair[k,2],]
		}
		if(par == TRUE) {
			Vilist = vector("list", max(grpindx)) 
			Sxx = matrix(0, nrow = p, ncol = p)
			Sxy = matrix(0, nrow = p, ncol = 1)
			for(i in 1:max(grpindx)) {
				dismat <- distGeoAni(xcoords[grpindx == i],
					ycoords[grpindx == i], xcoords[grpindx == i], 
					ycoords[grpindx == i], rotate = 90, 
					range = theta[2], minorp = 1)	
				covMatg <- theta[1]*corModelExponential(dismat) +
					theta[3]*diag(dim(dismat)[1])
				Vilist[[i]] = solve(covMatg)
				Sxx = Sxx + t(X[grpindx == i,]) %*% Vilist[[i]] %*% 
					X[grpindx == i,]
				Sxy = Sxy + t(X[grpindx == i,]) %*% Vilist[[i]] %*% 
					z[grpindx == i]
			}
			bhat = solve(Sxx) %*% Sxy

			ijpair = NULL
			for(i in 1:(max(grpindx)-1)) {
				ijpair = rbind(ijpair, 
					cbind(rep(i, times = max(grpindx) - i), (i+1):max(grpindx))
				)
			}
			Wxxlist = foreach(k=1:length(ijpair[,1])) %dopar% {
				Wxxloop(k, ijpair, theta, X, Vilist, xcoords, ycoords, grpindx)
			}
			Wxx = Reduce('+', Wxxlist)
			covb = solve(Sxx) + 2*solve(Sxx) %*% Wxx %*% solve(Sxx)
		}
		if(par == FALSE) {
			Vilist = vector("list", max(grpindx))
			bhatList = vector("list", max(grpindx))
			Sxx = matrix(0, nrow = p, ncol = p)
			Sxy = matrix(0, nrow = p, ncol = 1)
			for(i in 1:max(grpindx)) {
				dismat <- distGeoAni(xcoords[grpindx == i],
					ycoords[grpindx == i], xcoords[grpindx == i], 
					ycoords[grpindx == i], rotate = 90, 
					range = theta[2], minorp = 1)	
				covMatg <- theta[1]*corModelExponential(dismat) +
					theta[3]*diag(dim(dismat)[1])
				Vilist[[i]] = solve(covMatg)
				XViX = t(X[grpindx == i,]) %*% Vilist[[i]] %*% 
					X[grpindx == i,]
				XViz = t(X[grpindx == i,]) %*% Vilist[[i]] %*% 
					z[grpindx == i]
				bhatList[[i]] = cbind(solve(XViX) %*% XViz, diag(solve(XViX)))
				Sxx = Sxx + XViX
				Sxy = Sxy + XViz
			}
			bhat = solve(Sxx) %*% Sxy

			ijpair = NULL
			for(i in 1:(max(grpindx)-1)) 
				ijpair = rbind(ijpair, 
					cbind(rep(i, times = max(grpindx) - i), (i+1):max(grpindx))
				)
			Wxx = 0
			for(k in 1:length(ijpair[,1])) {
				Wxx = Wxx + Wxxloop(k, ijpair, theta, X, Vilist, 
					xcoords, ycoords, grpindx)
			}
			covb = solve(Sxx) + 2*solve(Sxx) %*% Wxx %*% solve(Sxx)
		}
	}

	outpt <- list(
	  call =  match.call(expand.dots = FALSE),
	  cope_est_meth = cope_est_meth,
	  spatial_name = spatial_name,
	  theta = theta,
	  extrap = extrap,
		bhat = bhat,
		covb = covb,
		p = dim(X)[2],
		n = length(z),
		X = X,
		z = z,
		theta = theta,
		Vilist = Vilist,
		bhatList = bhatList,
    Sxx = Sxx,
    Sxy = Sxy,
    Wxx = Wxx,
    ViX = ViX,
    Viz = Viz
	)

	class(outpt) <- "slm_fefe"
	outpt
}

