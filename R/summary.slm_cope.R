#' @export summary.slm_cope
summary.slm_cope <-
function(object, ...)
{
	theta = object$theta
#	attributes(theta) = attr_theta
	covp_est = data.frame(type = attr(theta,'label'),
		est = theta, fixed = attr(theta,'fixed'))
	names(covp_est) = c('Covariance Parm', 'Estimate', 'Fixed?')
  outpt = list(catCall = object$call,
    cov_parm_estimates = covp_est,
    m2LL = object$m2LL,
    cope_est_meth = object$cope_est_meth,
    AIC = AIC(object))
  class(outpt) <- "summary.slm_cope"
  outpt
}

#' @export print.slm_cope
print.slm_cope <- function(x,...) {
    print(summary(x,...),...)
}

#' @export print.summary.slm_cope
print.summary.slm_cope <- function(x, 
  digits = max(3L, getOption("digits")), ...)
{
  cat("\nCall:\n", paste(deparse(x$catCall), sep = "\n", collapse = "\n"), 
        "\n", sep = "")
  cat("\nSpatial Model:\n", as.character(x$catCall$spatial_model),"\n")
  cat("\nCovariance Parameter Estimation Method:\n", x$cope_est_meth, "\n") 
  cat("\nCovariance Parameter Estimates:\n")
  cpe = x$cov_parm_estimates
  cpe[,2] = format(round(cpe[,2], digits = digits), 
		scientific  = FALSE)
	print(cpe)
  cat("\nAIC:\n", x$AIC)
  cat("\n")
  invisible(x)
}
