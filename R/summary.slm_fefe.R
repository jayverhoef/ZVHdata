#' @export summary.slm_fefe
summary.slm_fefe <-
function(object, ...)
{

	b.hat <- object$bhat
	bhat.se <- sqrt(diag(object$covb))
	n <- object$sampinfo$obs.sample.size
	p <- object$p
	tvec = b.hat/bhat.se
	ptvec <- round(100000*(1 - pt(abs(b.hat/bhat.se), 
		df = object$n - object$p))*2)/100000
	fixed.effects.estimates <- data.frame(FactorLevel = row.names(b.hat), Estimate = b.hat,
			std.err = bhat.se, t.value = tvec, prob.t = ptvec)
	res = object$z - object$X %*% b.hat

  outpt = list(catCall = object$call,
    fixed.effects.estimates = fixed.effects.estimates,
    covariance.parameter.estimates = object$theta,
    res = res,
    spatial_name = object$spatial_name,
    cope_est_meth = object$cope_est_meth,
    theta = object$theta,
    extrap = object$extrap)
#    Rsquared = GR2(object))
  class(outpt) <- "summary.slm_fefe"
  outpt
}

#' @export print.slm_fefe
print.slm_fefe <- function(x,...) {
    print(summary(x,...))
}

#' @export print.summary.slm_fefe
print.summary.slm_fefe <- function(x, 
  digits = max(3L, getOption("digits") - 3L),
  signif.stars = getOption("show.signif.stars"), ...)
{
  cat("\nCall:\n", paste(deparse(x$catCall), sep = "\n", collapse = "\n"), 
        "\n", sep = "")
  cat("\nResiduals:\n")
  resQ = c(min(x$res), quantile(x$res, p = c(.25, .5, .75), 
    na.rm = TRUE), max(x$res))
  names(resQ) <- c("Min", "1Q", "Median", "3Q", "Max")
  print(resQ, digits = digits)
  cat("\nFixed Effects Estimation Method:\n GLS\n")
  cat("\nCoefficients:\n")
  coef1 = x$fixed.effects.estimates
  coefs = coef1[,2:5]
  row.names(coefs) = coef1[,1]
  colnames(coefs) = c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  printCoefmat(coefs, digits = digits, signif.stars = signif.stars, 
            na.print = "NA", ...)
  cat("\nSpatial Model:\n", as.character(x$spatial_name),"\n")
  cat("\nCovariance Parameters:\n")
	cpe = data.frame(type = c('parsil','range','nugget'),
		est = x$theta[1:3])
	if(!is.null(x$extrap) | length(x$theta) == 4) {
		if(!is.null(x$extrap)) {extrap = x$extrap } else
			{extrap = x$theta[4] }
		cpe = rbind(cpe, data.frame(type = 'extrap',
			est = extrap))
	}
	names(cpe) = c('Covariance Parm', 'Estimate')
	cpe[,2] = format(cpe[,2], scientific  = FALSE)
	print(cpe)
	rse = sum(x$theta[c(1,3)])
#  rse = sqrt(sum(cpe[cpe[,"Parameter"] == "parsill","Estimate"]))
  cat("\nResidual standard error:", rse)
#  cat("\nGeneralized R-squared:", x$Rsquared)

  cat("\n")
  invisible(x)
}
