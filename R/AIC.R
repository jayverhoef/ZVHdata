#'@export logLik.slm_cope
logLik.slm_cope <- function(object,...) {
    out = -object$m2LL/2
    attr(out,'Est_Meth') = object$call$estMeth
    out
}

#'@export AIC.slm_cope
AIC.slm_cope <- function(object,...,k = 2) {
    nparmstheta <- length(object$theta_ini)
    if(object$profile_sigma == TRUE) nparmstheta = nparmstheta + 1
    rankX <- sum(svd(object$X)$d>1e-10)
    out = object$m2LL + k*(nparmstheta + rankX)
    attr(out,'Est_Meth') = object$call$estMeth
    out
}

