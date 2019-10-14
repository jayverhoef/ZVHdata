#-------------------------------------------------------------------------------
#
#               pulo
#
#-------------------------------------------------------------------------------

#' Predict Method for Spatial Linear Mixed Model Fits
#'
#' Makes predictions for missing data or using external prediction data set. 
#'
#' @param x A fitted splm object
#' @param predData A data.frame object to be used for prediction.  It must have the same columns and names for the response variable, any covariates, and x-coordinates and y-coordinates.
#'
#' @details These are the universal kriging equations found in Cressie (1993, pg. 154-155).
#'
#' @return \code{predict.splm} produces a data.frame with 2 columns for predictions, and prediction standard errors for missing response variable, one with the response variable name appended by \code{.pred} for predictions and by \code{.predSE} for the prediction standard errors.
#'
#' @references \cite{Cressie, N.A.C. (1993) Statistics for Spatial Data. Wiley.}
#'
#' @author Jay Ver Hoef
#' @export
pulo <- function(cope_out, pred_data, theta = cope_out$theta, 
	extrap = cope_out$extrap, spatial_model = cope_out$spatial_model, 
	z = NULL, X = NULL, Xp = NULL, Zs = NULL, Zsp = NULL,
	xcoords = NULL, ycoords = NULL, xcoordsp = NULL, ycoordsp = NULL,
  predmeth = 'alldata', nNN = 50, par = FALSE)
{
	if(!is.null(cope_out)) {
		if(class(pred_data) == 'SpatialPointsDataFrame') {
			coords = sp::coordinates(pred_data)
			xcoordsp = coords[,1]
			ycoordsp = coords[,2]
			pDF = pred_data@data
			pDF[,cope_out$respCol] = rep(0, times = length(xcoordsp))
		}
		if(class(pred_data) == 'SpatialPoints') {
			coords = sp::coordinates(pred_data)
			xcoordsp = coords[,1]
			ycoordsp = coords[,2]
			pDF = data.frame(resp =  rep(0, times = length(xcoordsp)))
			colnames(pDF) = cope_out$respCol
		}
		trms <- terms(cope_out$formula, data = pDF)
		covList <- attr(trms,"term.labels")
		#design matrix
		Xp <- model.matrix(cope_out$formula, pDF)

		np <- length(xcoordsp)
		
		theta = cope_out$theta
		if(length(theta) == 4) extrap = theta[4]
		xcoords = cope_out$coords[,1]
		ycoords = cope_out$coords[,2]
		z = cope_out$z
		X = cope_out$X
	}

	parsil = theta[attr(theta, 'label') == 'partial sill']
	nugget = theta[attr(theta, 'label') == 'nugget']
	rang = theta[attr(theta, 'label') == 'range']
	rotate = 90
	if(any('rotate' %in% attr(theta, 'label')))
		rotate = theta[attr(theta, 'label') == 'rotate']
	minorp = 1
	if(any('minorp' %in% attr(theta, 'label')))
		minorp = theta[attr(theta, 'label') == 'minorp']
	extrap = NULL
	if(any('extrap' %in% attr(theta, 'label')))
		extrap = theta[attr(theta, 'label') == 'extrap']
		
	if(predmeth == 'alldata') {
		dismat <- distGeoAni(xcoords, ycoords, xcoords, ycoords, 
			rotate = rotate, range = rang, minorp = minorp)	
		covMat <- parsil*spatial_model(dismat, extrap) +
				nugget*diag(dim(dismat)[1])
		dismat <- distGeoAni(xcoords, ycoords, xcoordsp, 
				ycoordsp, rotate = rotate, 
				range = rang, minorp = minorp)
		Vpred <- parsil*spatial_model(dismat, extrap)
    qrlist = qr(covMat, LAPACK = TRUE)
    ViX = solve(qrlist, X)
    Viz = solve(qrlist, z)
    ViVpred = solve(qrlist, Vpred)
		XViX <- crossprod(X, ViX)
		covb <- solve(XViX)
		bhat <- covb %*% crossprod(ViX, z)
		sill <- parsil + nugget
		preds <- matrix(NA, nrow = length(xcoordsp), ncol = 2)
		preds[,1] <- apply(as.vector(Viz) * Vpred, 2, sum) +
			Xp %*% bhat - t(Vpred) %*% (ViX %*% bhat)	
		preds[,2] <- sqrt(rep(sill, times = nrow(Xp)) - 
			apply(ViVpred * Vpred, 2, sum) +
			apply((covb %*% t(Xp)) * t(Xp), 2, sum) -
			2*apply((covb %*% t(Xp)) * (t(X) %*% ViVpred), 2, sum) +
			apply(((covb %*% t(ViX)) %*% Vpred) * (t(X) %*% ViVpred), 2, sum))
		colnames(preds) <- c(paste0(cope_out$respCol, ".pred"), 
			paste(cope_out$respCol, ".predSE", sep = ""))

	}
	if(predmeth == 'nearnei') {
		dxy = as.matrix(cbind(xcoords,ycoords))
		pxy = as.matrix(cbind(xcoordsp,ycoordsp))
		nearxy = knn(data = dxy, query = pxy, k = nNN)
		if(par == FALSE) {
		  preds <- matrix(NA, nrow = length(xcoordsp), ncol = 2)
			for(i in 1:length(xcoordsp)) {
				dismat <- distGeoAni(xcoords[nearxy$nn.idx[i,]], 
					ycoords[nearxy$nn.idx[i,]], xcoordsp[i], 
					ycoordsp[i], rotate = rotate, 
					range = rang, minorp = minorp)
				Vpred <- parsil*spatial_model(dismat, extrap)
				dismat <- distGeoAni(xcoords[nearxy$nn.idx[i,]], 
					ycoords[nearxy$nn.idx[i,]], xcoords[nearxy$nn.idx[i,]], 
					ycoords[nearxy$nn.idx[i,]], rotate = rotate, 
					range = rang, minorp = minorp)
				V <- parsil*spatial_model(dismat, extrap)
				diag(V) = diag(V) + nugget					
				Xi <- X[nearxy$nn.idx[i,],]
				zi <- z[nearxy$nn.idx[i,]]
				qrlist = qr(V, LAPACK = TRUE)
				ViX = solve(qrlist, Xi)
				Viz = solve(qrlist, zi)
				ViVpred = solve(qrlist, Vpred)
				XViX <- crossprod(Xi, ViX)
				covb <- solve(XViX)
				bhat <- covb %*% crossprod(ViX, zi)
				sill <- parsil + nugget
				preds[i,1] = sum(as.vector((Viz)) * Vpred) +
					Xp[i,] %*% bhat - t(ViVpred) %*% (Xi %*% bhat)
				preds[i,2] = sqrt(sill -
					sum(ViVpred * Vpred) +
					sum((covb %*% Xp[i,]) * t(Xp[i,, drop = F])) -
					2*sum((covb %*% Xp[i,]) * (t(Xi) %*% ViVpred)) +
					sum((covb %*% t(Xi) %*% ViVpred) * 
					(t(Xi) %*% ViVpred)))
			}
		}
		if(par == TRUE) {
			NNloop = function(i, theta, X, Xp, xcoords, ycoords,
				xcoordsp, ycoordsp, nearxy) 
			{
		  predi <- matrix(NA, nrow = 1, ncol = 2)
				dismat <- distGeoAni(xcoords[nearxy$nn.idx[i,]], 
					ycoords[nearxy$nn.idx[i,]], xcoordsp[i], 
					ycoordsp[i], rotate = rotate, 
					range = rang, minorp = minorp)
				Vpred <- parsil*spatial_model(dismat, extrap)
				dismat <- distGeoAni(xcoords[nearxy$nn.idx[i,]], 
					ycoords[nearxy$nn.idx[i,]], xcoords[nearxy$nn.idx[i,]], 
					ycoords[nearxy$nn.idx[i,]], rotate = rotate, 
					range = rang, minorp = minorp)
				V <- parsil*spatial_model(dismat, extrap)
				diag(V) = diag(V) + nugget					
				Xi <- X[nearxy$nn.idx[i,],]
				zi <- z[nearxy$nn.idx[i,]]
				qrlist = qr(V, LAPACK = TRUE)
				ViX = solve(qrlist, Xi)
				Viz = solve(qrlist, zi)
				ViVpred = solve(qrlist, Vpred)
				XViX <- crossprod(Xi, ViX)
				covb <- solve(XViX)
				bhat <- covb %*% crossprod(ViX, zi)
				sill <- parsil + nugget
				predi[1,1] = sum(as.vector(Viz) * Vpred) +
					Xp[i,] %*% bhat - t(ViVpred) %*% (Xi %*% bhat)
				predi[1,2] = sqrt(sill -
					sum(ViVpred * Vpred) +
					sum((covb %*% Xp[i,]) * t(Xp[i,, drop = F])) -
					2*sum((covb %*% Xp[i,]) * (t(Xi) %*% ViVpred)) +
					sum((covb %*% t(Xi) %*% ViVpred) * 
					(t(Xi) %*% ViVpred)))
				predi
			}
			NNout = foreach(j=1:length(xcoordsp)) %dopar% {
				NNloop(j, theta, X, Xp, xcoords, ycoords, xcoordsp, 
					ycoordsp, nearxy)
			}
			preds = t(matrix(unlist(NNout), nrow = 2))
		}
	}
	preds
}


