#-------------------------------------------------------------------------------
#
#          m2LLg
#
#-------------------------------------------------------------------------------

# minus 2 times loglikelihood for big data using subsampling
#
# fits a geostatistical linear model
#
# @param theta vector of estimated covariance parameters
# @param z vector of data
# @param X design matrix for fixed effects
# @param Z list of design matrices for each random effect 
# @param xcoords vector with the x-coordinates
# @param ycoords vector with the y-coordinates
# @param estMeth estimation method.  Default is "REML" for restricted maximum likelihood.  Other options are "ML" for maximum likelihood
#
# @return minus 2 times the loglikelihood
#
# @author Jay Ver Hoef
#'@export
#'@export
m2LLprof <- function(theta_star, z, X, xcoords, ycoords, spatial_model,
	estMeth, grpindx, max_range, attr_theta, fixed_parms, Zs)
{
  if(any(abs(theta_star) > 20)) return(1e+30)
  attributes(theta_star) = attr_theta
  theta = unlist(japply(theta_star, attr(theta_star, "transfType")))*
			attr(theta_star,'scaleFactors')
  attributes(theta) = attr_theta
	if(!is.null(fixed_parms)) {
		theta = c(theta, fixed_parms)
		attr(theta,'label') = c(attr(theta_star,'label'), attr(fixed_parms,'label'))
  }
  if(!any('rotate' %in% attr(theta,'label'))) {
		rotate = 90
		minorp = 1
	} else {
		rotate = theta[attr(theta,'label') == 'rotate']
		minorp = theta[attr(theta,'label') == 'minorp']
	}
	if(theta[attr(theta,'label') =='range'] > max_range) return (1e+30)
	extrap = NULL
	if(any('extrap' %in% attr(theta,'label'))) 
		extrap = theta[attr(theta,'label') == 'extrap']
	parsil = theta[attr(theta,'label') == 'partial sill']
	if(is.null(Zs)) {
		nugget = 1 - parsil
	} else {
		namesZ = names(Zs)
		ps = c(1 - parsil, (1 - parsil)*
			cumprod(1 - theta[attr(theta,'label') %in% namesZ]))*
			c(theta[attr(theta,'label') %in% namesZ],1)
		vcs = ps[1:length(namesZ)]
		nugget = ps[length(ps)]
	}
  p = length(X[1,])
  n = length(X[,1])
  qrlist = vector("list", max(grpindx))
  Sxx = matrix(0, nrow = p, ncol = p)
  Sxy = matrix(0, nrow = p, ncol = 1)
  logDetV = 0
  logDetXViX = 0
  for(i in 1:max(grpindx)) {
	  dismat <- distGeoAni(xcoords[grpindx == i], ycoords[grpindx == i], 
      xcoords[grpindx == i], ycoords[grpindx == i], 
      rotate = rotate, 
      range = theta[attr(theta,'label') == 'range'], 
      minorp = minorp)
	  covMat <- parsil*spatial_model(dismat, extrap) +
      (nugget + 1e-6)*diag(dim(dismat)[1])
    if(!is.null(Zs)) {
			for(j in 1:length(Zs))
			covMat = covMat + vcs[j]*Zs[[j]][grpindx == i,] %*% 
        t(Zs[[j]][grpindx == i,])
		}
    qrlist[[i]] = qr(covMat, LAPACK = TRUE)
    ViX = solve(qrlist[[i]], X[grpindx == i,])
    XViX = crossprod(X[grpindx == i,],ViX)
    Sxx = Sxx + XViX
    Sxy = Sxy + t(crossprod(z[grpindx == i],ViX)) 
    logDetV = logDetV + sum(log(abs(diag(qr.R(qrlist[[i]])))))
    logDetXViX = logDetXViX + as.numeric(determinant(XViX, 
      logarithm = TRUE)$modulus)
  }
  betaHat = solve(Sxx) %*% Sxy
  rVir = 0
  for(i in 1:max(grpindx)) 
    rVir = rVir + t(z[grpindx == i] - 
      X[grpindx == i,] %*% betaHat) %*%
      solve(qrlist[[i]], (z[grpindx == i] - 
      X[grpindx == i,] %*% betaHat))
	minus2LL1 <- logDetV + n*log(rVir) + n + n*log(2*pi/n)
	if(estMeth == "REML") minus2LL1 <- minus2LL1 + logDetXViX - 
		p*log(rVir) - p - n*log(2*pi/n) + (n - p)*log(2*pi/(n-p))
	m2LL = minus2LL1
	
	if(any('rotate' %in% attr(theta_star,'label'))) {
		Sxx = matrix(0, nrow = p, ncol = p)
		Sxy = matrix(0, nrow = p, ncol = 1)
		logDetV = 0
		logDetXViX = 0
		for(i in 1:max(grpindx)) {
			dismat <- distGeoAni(xcoords[grpindx == i], ycoords[grpindx == i], 
				xcoords[grpindx == i], ycoords[grpindx == i], 
				rotate = 180 - rotate, 
				range = theta[attr(theta,'label') == 'range'], 
				minorp = minorp)
			covMat <- parsil*spatial_model(dismat, extrap) +
				(nugget + 1e-6)*diag(dim(dismat)[1])
			if(!is.null(Zs)) {
				for(j in 1:length(Zs))
				covMat = covMat + vcs[j]*Zs[[j]][grpindx == i,] %*% 
          t(Zs[[j]][grpindx == i,])
			}
			qrlist[[i]] = qr(covMat, LAPACK = TRUE)
			ViX = solve(qrlist[[i]], X[grpindx == i,])
			XViX = crossprod(X[grpindx == i,],ViX)
			Sxx = Sxx + XViX
			Sxy = Sxy + t(crossprod(z[grpindx == i],ViX)) 
			logDetV = logDetV + sum(log(abs(diag(qr.R(qrlist[[i]])))))
			logDetXViX = logDetXViX + as.numeric(determinant(XViX, 
				logarithm = TRUE)$modulus)
		}
		betaHat = solve(Sxx) %*% Sxy
		rVir = 0
		for(i in 1:max(grpindx)) 
			rVir = rVir + t(z[grpindx == i] - 
				X[grpindx == i,] %*% betaHat) %*%
				solve(qrlist[[i]], (z[grpindx == i] - 
				X[grpindx == i,] %*% betaHat))
    minus2LL2 <- logDetV + n*log(rVir) + n + n*log(2*pi/n)		
    if(estMeth == "REML") minus2LL2 <- minus2LL2 + logDetXViX - 
			p*log(rVir) - p - n*log(2*pi/n) + (n - p)*log(2*pi/(n-p))
		
		m2LL =  as.numeric(min(minus2LL1,minus2LL2))
		attr(m2LL, 'which') = which(c(minus2LL1,minus2LL2) == m2LL)[1]
	}
	return(m2LL)
}

#'@export
m2LLa <- function(theta_star, z, X, xcoords, ycoords, spatial_model,
	estMeth, grpindx, max_range, attr_theta, fixed_parms, Zs)
{
  if(any(abs(theta_star) > 30)) return(1e+30)
  attributes(theta_star) = attr_theta
  theta = unlist(japply(theta_star, attr(theta_star, "transfType")))*
			attr(theta_star,'scaleFactors')
  attributes(theta) = attr_theta
	if(!is.null(fixed_parms)) {
		theta = c(theta, fixed_parms)
		attr(theta,'label') = c(attr(theta_star,'label'), attr(fixed_parms,'label'))
  }
  if(!any('rotate' %in% attr(theta,'label'))) {
		rotate = 90
		minorp = 1
	} else {
		rotate = theta[attr(theta,'label') == 'rotate']
		minorp = theta[attr(theta,'label') == 'minorp']
	}
	if(theta[attr(theta,'label') =='range'] > max_range) return (1e+30)
	extrap = NULL
	if(any('extrap' %in% attr(theta,'label'))) 
		extrap = theta[attr(theta,'label') == 'extrap']
  p = length(X[1,])
  n = length(X[,1])
  qrlist = vector("list", max(grpindx))
  Sxx = matrix(0, nrow = p, ncol = p)
  Sxy = matrix(0, nrow = p, ncol = 1)
  logDetV = 0
  logDetXViX = 0
  for(i in 1:max(grpindx)) {
	  dismat <- distGeoAni(xcoords[grpindx == i], ycoords[grpindx == i], 
      xcoords[grpindx == i], ycoords[grpindx == i], 
      rotate = rotate, 
      range = theta[attr(theta,'label') == 'range'], 
      minorp = minorp)
	  covMat <- theta[attr(theta,'label') == 'partial sill']*
			spatial_model(dismat, extrap) +
      (theta[attr(theta,'label') == 'nugget'] + 1e-6)*
      diag(dim(dismat)[1])
    if(!is.null(Zs)) {
			namesZ = names(Zs)
			for(j in 1:length(Zs))
			covMat = covMat + theta[attr(theta,'label') == namesZ[j]]*
				Zs[[j]][grpindx == i,] %*% t(Zs[[j]][grpindx == i,])
		}
    qrlist[[i]] = qr(covMat, LAPACK = TRUE)
    ViX = solve(qrlist[[i]], X[grpindx == i,])
    XViX = crossprod(X[grpindx == i,],ViX)
    Sxx = Sxx + XViX
    Sxy = Sxy + t(crossprod(z[grpindx == i],ViX)) 
    logDetV = logDetV + sum(log(abs(diag(qr.R(qrlist[[i]])))))
    logDetXViX = logDetXViX + as.numeric(determinant(XViX, 
      logarithm = TRUE)$modulus)
  }
  betaHat = solve(Sxx) %*% Sxy
  rVir = 0
  for(i in 1:max(grpindx)) 
    rVir = rVir + t(z[grpindx == i] - 
      X[grpindx == i,] %*% betaHat) %*%
      solve(qrlist[[i]], (z[grpindx == i] - 
      X[grpindx == i,] %*% betaHat))
	minus2LL1 <- logDetV + rVir + n*log(2*pi)
	if(estMeth == "REML") minus2LL1 <- minus2LL1 + logDetXViX - p*log(2*pi)
	m2LL = minus2LL1
	
	if(any('rotate' %in% attr(theta_star,'label'))) {
		Sxx = matrix(0, nrow = p, ncol = p)
		Sxy = matrix(0, nrow = p, ncol = 1)
		logDetV = 0
		logDetXViX = 0
		for(i in 1:max(grpindx)) {
			dismat <- distGeoAni(xcoords[grpindx == i], ycoords[grpindx == i], 
				xcoords[grpindx == i], ycoords[grpindx == i], 
				rotate = 180 - theta[attr(theta,'label') == 'rotate'], 
				range = theta[attr(theta,'label') == 'range'], 
				minorp = theta[attr(theta,'label') == 'minorp'])
			covMat <- theta[attr(theta,'label') == 'partial sill']*
				spatial_model(dismat, extrap) +
				(theta[attr(theta,'label') == 'nugget'] + 1e-6)*
				diag(dim(dismat)[1])
			if(!is.null(Zs)) {
				namesZ = names(Zs)
				for(j in 1:length(Zs))
				covMat = covMat + theta[attr(theta,'label') == namesZ[j]]*
					Zs[[j]][grpindx == i,] %*% t(Zs[[j]][grpindx == i,])
			}			
			qrlist[[i]] = qr(covMat, LAPACK = TRUE)
			ViX = solve(qrlist[[i]], X[grpindx == i,])
			XViX = crossprod(X[grpindx == i,],ViX)
			Sxx = Sxx + XViX
			Sxy = Sxy + t(crossprod(z[grpindx == i],ViX)) 
			logDetV = logDetV + sum(log(abs(diag(qr.R(qrlist[[i]])))))
			logDetXViX = logDetXViX + as.numeric(determinant(XViX, 
				logarithm = TRUE)$modulus)
		}
		betaHat = solve(Sxx) %*% Sxy
		rVir = 0
		for(i in 1:max(grpindx)) 
			rVir = rVir + t(z[grpindx == i] - 
				X[grpindx == i,] %*% betaHat) %*%
				solve(qrlist[[i]], (z[grpindx == i] - 
				X[grpindx == i,] %*% betaHat))
		minus2LL2 <- logDetV + rVir + n*log(2*pi)
		if(estMeth == "REML") minus2LL2 <- minus2LL2 + logDetXViX - p*log(2*pi)
		
		m2LL =  as.numeric(min(minus2LL1,minus2LL2))
		attr(m2LL, 'which') = which(c(minus2LL1,minus2LL2) == m2LL)[1]
	}
	return(m2LL)
}

#-------------------------------------------------------------------------------
#
#           distGeoAni
#
#-------------------------------------------------------------------------------

# Compute anistropy corrected distance between two sets of data
#
# computes anistropy corrected distance between two sets of data
#
# @param xrow vector with x-coordinates that will form rows of distance matrix 
# @param yrow vector with y-coordinates that will form rows of distance matrix, must be of same length as xrow
# @param xcol vector with x-coordinates that will form columns of distance matrix 
# @param ycol vector with y-coordinates that will form columns of distance matrix, must be of same length as xcol
# @param rotate rotation of anisotropic axes, default = 0
# @param range range of autocorrelation model, default = 1
# @param minorp proportion of range in x direction to that of y direction for unrotated anisotropic model, default = 1
#
# @return matrix of distances
#
# @author Jay Ver Hoef
#'@export
distGeoAni <- function(xrow, yrow, xcol, ycol, rotate = 0, range = 1, minorp = 1)
{
	# total number of observations for each set of coordinates  
		n.rows <- length(xrow)
		n.cols <- length(xcol)
	# expand all x-coordinates
		sxr <- matrix(xrow, ncol = 1) %*% 
			matrix(rep(1,times = n.cols), nrow = 1)
		sxc <- matrix(rep(1,times = n.rows), ncol = 1) %*% 
			matrix(xcol, nrow = 1)
		syr <- matrix(yrow,ncol = 1) %*% 
			matrix(rep(1,times = n.cols), nrow = 1)
		syc <- matrix(rep(1,times = n.rows), ncol = 1) %*% 
			matrix(ycol, nrow = 1)
	# find difference in coordinates between all pairwise locations
		sxdif <- sxr - sxc
		sydif <- syr - syc
	# rotate coordinates
		newx <- cos(rotate*2*pi/360)*sxdif - sin(rotate*2*pi/360)*sydif
		newy <- sin(rotate*2*pi/360)*sxdif + cos(rotate*2*pi/360)*sydif
	# scale coordinates by minor and major axes */
		newx <- newx/(range*minorp)
		newy <- newy/range
	# compute distance for the scaled and rotated coordinates */
		sqrt(newx^2 + newy^2)
}



# --------------- A BUNCH OF SPATIAL CORRELATION MODELS

#'@export
exponential <- function(distance.matrix, extrap)
{
	exp(-distance.matrix) 
}

#'@export
expRadon2 <- function(distance.matrix, extrap)
{
	(1 + distance.matrix)*exp(-distance.matrix) 
}

#'@export
expRadon4 <- function(distance.matrix, extrap)
{
	(1 + distance.matrix + distance.matrix^2/3)*exp(-distance.matrix) 
}

#'@export
Gaussian <- function(distance.matrix, extrap)
{
	exp(-distance.matrix^2) 
}

#'@export
stable <- function(distance.matrix, extrap)
{
	exp(-distance.matrix^(2*extrap/(1 + extrap))) 
}

#'@export
rationalQuad <- function(distance.matrix, extrap)
{
	1/(1+distance.matrix^2)
}

#'@export
cauchyGrav <- function(distance.matrix, extrap)
{
	1/sqrt(1+distance.matrix^2)
}

#'@export
cauchyMag <- function(distance.matrix, extrap)
{
	1/(sqrt(1+distance.matrix^2))^3
}

#'@export
cauchy <- function(distance.matrix, extrap)
{
	1/(1+distance.matrix^2)^extrap
}

#'@export
circular <- function(distance.matrix, extrap)
{
	d <- distance.matrix
	d[distance.matrix > 1] <- 0
	CovMat <- 2*(acos(d) - d*sqrt(1 - d^2))/pi
	CovMat[distance.matrix >= 1] <- 0
	CovMat
}

#'@export
spherical <- function(distance.matrix, extrap)
{
	CovMat <- (1 - 1.5*distance.matrix + 0.5*distance.matrix^3)
	CovMat[distance.matrix > 1] <- 0
	CovMat
}

#'@export
cubic <- function(distance.matrix, extrap)
{
	CovMat <- (1 - 7*distance.matrix^2 + 35*distance.matrix^3/4 - 7*distance.matrix^5/2
		+ 3*distance.matrix^7/4)
	CovMat[distance.matrix > 1] <- 0
	CovMat
}

#'@export
penta <- function(distance.matrix, extrap)
{
	CovMat <- (1 - 22*distance.matrix^2/3 + 33*distance.matrix^4 - 77*distance.matrix^5/2
		+ 33*distance.matrix^7/2 - 11*distance.matrix^9/2 + 5*distance.matrix^11/6)
	CovMat[distance.matrix > 1] <- 0
	CovMat
}

#'@export
cardinalSine <- function(distance.matrix, extrap)
{
	d <- distance.matrix
	d[distance.matrix == 0] <- 1
	CorMat <- sin(d)/d
	CorMat[distance.matrix == 0] <- 1
	CorMat
}

#'@export
besselj <- function(distance.matrix, extrap)
{
	d <- distance.matrix
	d[distance.matrix == 0] <- 1
	extrap <- extrap + 2.0
	CorMat <- d^(1-extrap/2)*besselJ(d, extrap/2 - 1)*2^(extrap/2 - 1)*gamma(extrap/2)
	CorMat[distance.matrix == 0] <- 1
	CorMat
}

#'@export
besselk <- function(distance.matrix, extrap)
{
	d <- distance.matrix
	d[distance.matrix == 0] <- 1
	CorMat <- d^extrap*besselK(d, extrap)/(2^(extrap - 1)*gamma(extrap))
	CorMat[distance.matrix == 0] <- 1
	CorMat
}

#'@export
expit = function(x)
{
	exp(x)/(1 + exp(x))
}

#'@export
logit = function(x)
{
	log(x/(1-x))
}

#'@export
japply = function(listvals, funs)
{
	out = vector('list',length(listvals))
	for(i in 1:length(listvals))
	  out[[i]] = do.call(funs[[i]],as.list(listvals[[i]]))
	out
}

#'@export
log_lik <- function(theta, z, X, xcoords, ycoords, spatial_model,
	estMeth, attr_theta, grpindx, Zs)
{

#  attributes(theta) = attr_theta
  if(!any('rotate' %in% attr(theta,'label'))) {
		rotate = 90
		minorp = 1
	} else {
		rotate = theta[attr(theta,'label') == 'rotate']
		minorp = theta[attr(theta,'label') == 'minorp']
	}
	extrap = NULL
	if(any('extrap' %in% attr(theta,'label'))) 
		extrap = theta[attr(theta,'label') == 'extrap']
	parsil = 0
	if(any('partial sill' %in% attr(theta,'label'))) 
		parsil = theta[attr(theta,'label') == 'partial sill']
	nugget = 0
	if(any('nugget' %in% attr(theta,'label'))) 
		nugget = theta[attr(theta,'label') == 'nugget']
	rang = 0
	if(any('range' %in% attr(theta,'label'))) 
		rang = theta[attr(theta,'label') == 'range']

  p = length(X[1,])
  n = length(X[,1])
  qrlist = vector("list", max(grpindx))
  Sxx = matrix(0, nrow = p, ncol = p)
  Sxy = matrix(0, nrow = p, ncol = 1)
  logDetV = 0
  logDetXViX = 0
  for(i in 1:max(grpindx)) {
	  dismat <- distGeoAni(xcoords[grpindx == i], ycoords[grpindx == i], 
      xcoords[grpindx == i], ycoords[grpindx == i], 
      rotate = rotate, 
      range = theta[attr(theta,'label') == 'range'], 
      minorp = minorp)
	  covMat <- parsil*spatial_model(dismat, extrap) +
      (nugget + 1e-6)*diag(dim(dismat)[1])
    if(!is.null(Zs)) {
			namesZ = names(Zs)
			for(j in 1:length(Zs))
        covMat = covMat + theta[attr(theta,'label') == namesZ[j]]*
          Zs[[j]][grpindx == i,] %*% t(Zs[[j]][grpindx == i,])
		}
    qrlist[[i]] = qr(covMat, LAPACK = TRUE)
    ViX = solve(qrlist[[i]], X[grpindx == i,])
    XViX = crossprod(X[grpindx == i,],ViX)
    Sxx = Sxx + XViX
    Sxy = Sxy + t(crossprod(z[grpindx == i],ViX)) 
    logDetV = logDetV + sum(log(abs(diag(qr.R(qrlist[[i]])))))
    logDetXViX = logDetXViX + as.numeric(determinant(XViX, 
      logarithm = TRUE)$modulus)
  }
  betaHat = solve(Sxx) %*% Sxy
  rVir = 0
  for(i in 1:max(grpindx)) 
    rVir = rVir + t(z[grpindx == i] - 
      X[grpindx == i,] %*% betaHat) %*%
      solve(qrlist[[i]], (z[grpindx == i] - 
      X[grpindx == i,] %*% betaHat))
	m2LL <- logDetV + rVir + n*log(2*pi)
	if(estMeth == "REML") m2LL <- m2LL + logDetXViX - p*log(2*pi)
	loglik = -m2LL/2
	
	return(loglik)
}
