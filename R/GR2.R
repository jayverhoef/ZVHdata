GR2 <-
function(x)
{
	theta = x$theta
	extrap = x$extrap
	xcoords = x$coords[,1]
	ycoords = x$coords[,2]
	z = x$z
	X = x$X
  spatial_model = x$spatial_model
  
	dismat <- distGeoAni(xcoords, ycoords, xcoords, ycoords, 
		rotate = 90, range = theta[2], minorp = 1)	
	if(length(theta) == 4) extrap = theta[4]
	covMat <- theta[1]*spatial_model(dismat, extrap) +
			theta[3]*diag(dim(dismat)[1])

	QR = qr(covMat)
	ViX <- solve(covMat, X)
	Vi1 <- solve(covMat, rep(1, times = length(z)))
	Viz = solve(covMat,z)
	XViX <- crossprod(X, ViX)
	covb <- solve(XViX)
	bhat <- covb %*% crossprod(ViX, z)
  muhat <- sum(Viz)/sum(Vi1)
  1 - t(z - X %*% bhat) %*% (Viz - ViX %*% bhat)/
    t(z - rep(muhat, times = length(z))) %*% (Viz - Vi1 * muhat)
}

