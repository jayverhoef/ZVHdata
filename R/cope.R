#-------------------------------------------------------------------------------
#
#           cope
#
#-------------------------------------------------------------------------------

#' Covariance parameter estimation for a geostatistical linear model
#'
#' Covariance parameter estimation for a geostatistical linear model
#'
#' @param formula an R linear model formula
#' @param data an data object with spatial coordinates.  If a plain data.frame, then x_column and y_column must be specified.  If an sp object, x_column and y_column should be NULL (the default).
#' @param x_column name, in quotes, of the column containing the x-coordinate (if it was not possible to obtain the coordinates from the data class). Default is NULL. 
#' @param y_column name, in quotes, of the column containing the x-coordinate (if it was not possible to obtain the coordinates from the data class). Default is NULL. 
#' @param estMeth estimation method.  Default is "REML" for restricted maximum likelihood.  Other options are "ML" for maximum likelihood
#' @param spatial_model spatial autocorrelation models 
#' for random errors.  The list of spatial autocorrelation 
#' models is "exponential","expRadon2","expRadon4","gaussian","stable",
#' "rationalQuad","cauchyGrav","cauchyMag","cauchy","circular","spherical",
#' "cubic","penta","cardinalSine","besselK","besselJ"  Default is "exponential".
#' @param random_formula a list of variance components, 
#' Any names in the list not given above will be searched among the columns in the data set and
#' used as a factor variable for levels of a traditional random effect.
#' @param useAnisotropy include anistropy in parameter estimation?  Default is "FALSE"
#' @param subsample_col A column of factors indicating grouping for use in subsampling for large data sets.  Default is "NULL," which creates a single grouping.

#'
#' @return a list of class "slm_cope".  The functions "summary" and "print" are used to obtain and print a summary. "anova" returns just the analysis of variance table...
#'
#' @author Jay Ver Hoef
#' @export
cope <- function(formula, data, x_column = NULL, y_column = NULL,
	spatial_model = exponential, extrap = NULL, max_range = NULL, 
	use_anisotropy = FALSE, random_formula = NULL,
	fixed_parms = NULL, thetaini = NULL, estMeth = "REML", 
  profile_sigma = TRUE, use_grid = FALSE, subSampCol = NULL, 
  par = FALSE, ...) 
{

	# ----------------------------------------------------------------------------
	# prepare data
	# ----------------------------------------------------------------------------

	# handle data from sp package
	if(class(data) == 'SpatialPointsDataFrame') {
		coords = sp::coordinates(data)
		xcoords = coords[,1]
		ycoords = coords[,2]
		DF = data@data
	# handle data for data.frame
	} else {
		DF = as.data.frame(data)
		xcoords = DF[,x_column]
		ycoords = DF[,y_column]
		coords = cbind(xcoords, ycoords)
	}
	# get the column name for the response variable
	trms <- terms(formula, data = DF)
	respCol <- as.character(as.list(attr(trms,"variables")[-1]))[1]
	# design matrix
	X <- model.matrix(formula, DF)
	# column of observations of response variable
	z <- as.matrix(DF[, respCol])
	# if formulas for random effects, create Z matrices
	Zs = NULL
	if(!is.null(random_formula)) {
		Zs = vector('list', length(random_formula))
		for(i in 1:length(random_formula))
			Zs[[i]] = model.matrix( random_formula[[i]], data)
			names(Zs) = names(random_formula)
	}
	# sample size (number of non-missing records
	n <- length(z)
	# dimension of the design matrix
	p <- sum(svd(X)$d>1e-10)
	# creating grouping index for big data methods of all 1's if missing
	if(is.null(subSampCol)) {grpindx = rep(1, times = n)
	} else { grpindx = DF[,subSampCol] }
	
	# keep track of optimization time
	starttime = Sys.time()
	#do not profile overall variance parameter if any variance components
	# are fixed
	Test2 = FALSE
	if(!is.null(fixed_parms))  
		Test2 = any(attr(fixed_parms,'label') %in% 
		c('partial sill', 'nugget', names(random_formula)))
	if(!profile_sigma | Test2) {
		if(profile_sigma & Test2)
			cat('Profiling disabled when holding any variance components fixed.\r
			... working ...\n')
		profile_sigma = FALSE
		# create initial parameter estimate
		if(is.null(thetaini)) {
			# start with 3 usual spatial covariance parameters
			tlength = 3
			label = c('partial sill', 'range', 'nugget')
			# scaling factors. Will expand parameters by these factors inside
			# optimization routine
			varz = summary(lm(formula, data))$sigma^2
			scaleFactors <- c(.4*varz, 
					sqrt((var(xcoords) + var(ycoords))/2), .4*varz)
			# transformatin functions.  Will transform parameters by these functions
			# inside optimization routine
			transfType = c(exp, exp, exp)
			# if using anisotropy, add two parameters and scaling factors and
			# transformation functions
			if(use_anisotropy == TRUE) {
				tlength = tlength + 2
				label = c(label, 'minorp', 'rotate')
				scaleFactors = c(scaleFactors, 1, 90)
				transfType = c(transfType, expit, expit)
			}
			# check for models with extra parameters, and add it to parameter
			# initialization, with scaling factors and transformation functions
			if(any(as.character(match.call()$spatial_model) %in% c('cauchy', 'stable',
				'besselk', 'besselj'))) {
					tlength = tlength + 1
					label = c(label, 'extrap')
					scaleFactors = c(scaleFactors, 1)
					transfType = c(transfType, exp)
			}
			# check for random effects
			if(!is.null(random_formula)) {
        if(is.null(names(random_formula)))
          stop('random_formula argument must be a NAMED list')
				tlength = tlength + length(random_formula)
				label = c(label, names(random_formula))
				scaleFactors = c(scaleFactors, 
					rep(.2*varz, times = length(random_formula)))
				for(i in 1:length(random_formula))
					transfType = c(transfType, exp)
			}		
			# start with all parameters at 0 (essentially unscaled with an 
			# inverse transformation
			thetaini <- matrix(rep(0, times = tlength), ncol  = 1)
			rownames(thetaini) <- 1:length(thetaini)
			colnames(thetaini) <- "parms"
			attr(thetaini,'label') = label
			attr(thetaini,'scaleFactors') = scaleFactors
			attr(thetaini,'transfType') = transfType 
		} else {
		# do inverse tranform to get scaled parameters
			scaleFactors = NULL
			transfType = NULL
			label = NULL
      tini_new = NULL
			varz = summary(lm(formula, data))$sigma^2
			if(any(attr(thetaini,'label') == 'partial sill')) {
				tini_new = c(tini_new,
					log(thetaini[attr(thetaini,'label') == 'partial sill']/(.4*varz)))
        label = c(label, 'partial sill')
				scaleFactors = c(scaleFactors, .4*varz)
				transfType = c(transfType, exp)
			}
			if(any(attr(thetaini,'label') == 'range')) {
				tini_new = c(tini_new, 
					log(thetaini[attr(thetaini,'label') == 'range']/
					(sqrt((var(xcoords) + var(ycoords))/2))))
        label = c(label, 'range')
				scaleFactors = c(scaleFactors, sqrt((var(xcoords) + var(ycoords))/2))
				transfType = c(transfType, exp)
			}
			if(any(attr(thetaini,'label') == 'nugget')) {
				tini_new = c(tini_new,
					log(thetaini[attr(thetaini,'label') == 'nugget']/(.4*varz)))
        label = c(label, 'nugget')
				scaleFactors = c(scaleFactors, .4*varz)
				transfType = c(transfType, exp)
			}
			if(any(attr(thetaini,'label') == 'minorp')) {
				tini_new = c(tini_new, 
					logit(thetaini[attr(thetaini,'label') == 'minorp']))
        label = c(label, 'minorp')
				scaleFactors = c(scaleFactors, 1)
				transfType = c(transfType, expit)
			}
			if(any(attr(thetaini,'label') == 'rotate')) {
				tini_new = c(tini_new, 
					logit(min(thetaini[attr(thetaini,'label') == 'rotate'],
					(180 - thetaini[attr(thetaini,'label') == 'rotate']))/90))
        label = c(label, 'rotate')
				scaleFactors = c(scaleFactors, 90)
				transfType = c(transfType, expit)
			}
			if(any(attr(thetaini,'label') == 'extrap')) {
				tini_new = c(tini_new, 
					log(thetaini[attr(thetaini,'label') == 'extrap']))
        label = c(label, 'extrap')
				scaleFactors = c(scaleFactors, 1)
				transfType = c(transfType, exp)
			}
			# check for random effects
			if(!is.null(random_formula)) {
				namesZ = names(random_formula)
				for(j in 1:length(namesZ)) {
					tini_new = c(tini_new,
					log(thetaini[attr(thetaini,'label') == namesZ[j]]/(.2*varz)))
          label = c(label, namesZ[j])
					scaleFactors = c(scaleFactors, .2*varz)
					transfType = c(transfType, exp)
				}
			}
      thetaini = tini_new
			attr(thetaini,'label') = label
			attr(thetaini,'scaleFactors') = scaleFactors
			attr(thetaini,'transfType') = transfType 
		}
		# set maximum range value if none given
		if(is.null(max_range)) max_range = 10*
					sqrt((var(xcoords) + var(ycoords))/2)
		# if some parameters are fixed, set those here and subset the
		# initial parameter estimate
		if(!is.null(fixed_parms)) {
			ind4opt = !attr(thetaini,'label') %in% attr(fixed_parms, 'label')
			label = attr(thetaini,'label')
			thetaini = thetaini[ind4opt]
			label = label[ind4opt]
			scaleFactors = scaleFactors[ind4opt]
			transfType = transfType[ind4opt] 
      # store the attributes to pass to optimization function
      attr(thetaini,'label') = label
      attr(thetaini,'scaleFactors') = scaleFactors
      attr(thetaini,'transfType') = transfType 
		}
    attrs = attributes(thetaini)
		# if no parallel processing...
		if(par == FALSE) {
			# if only 1 optimization parameter, use optimize()
			if(length(thetaini) == 1) {
				# estimate parameters
				parmest = optimize(m2LLa, lower = -30, upper = 30,
					z = z, X = X, xcoords = xcoords, ycoords = ycoords, 
					spatial_model = spatial_model, estMeth = estMeth, 
					grpindx = grpindx, max_range = max_range, 
					attr_theta = attrs, fixed_parms = fixed_parms,
					Zs = Zs)
				# get the minimized objective function value, equal to 
				# minus 2 times the log-likelihood
				m2LL <- parmest$objective
				# parameter estimate at minimized log-likelihood
				theta_parms = parmest$minimum
			# if multiple parameters, use optim()
			} else {
				if(use_grid == TRUE & length(thetaini) > 1 & length(thetaini) < 9) {
					# do a search grid to get good starting values
					# set the dimensions of grid search depending on number of parameters
					if(length(thetaini) == 2) gbase = 7
					if(length(thetaini) == 3) gbase = 5
					if(length(thetaini) %in% 4:5) gbase = 3
					if(length(thetaini) > 5) gbase = 2
					# get current loglik value for initial theta values
					min_m2LL = m2LLa(thetaini, z = z, X = X, 
						xcoords = xcoords, ycoords = ycoords, 
						spatial_model = spatial_model, estMeth = estMeth, 
						grpindx = grpindx, max_range = max_range, 
						attr_theta = attrs, fixed_parms = fixed_parms,
						Zs = Zs)
					# create a matrix of parameters values to test
					theta_test = NULL	
					for(j in 0:(length(thetaini)-1))	
            theta_test = cbind(theta_test,
              rep(1, times = 2^j) %x% 0:(gbase-1)  %x% 
                rep(1, times = 2^(length(thetaini)-j-1)))
          theta_test = 2*(theta_test/(gbase-1) - 0.5)
					# theta_test is centered on zero.  Now, center it on thetaini
					theta_test = as.vector(thetaini) + t(theta_test)
					# loop through grids values to find smallest loglik
					# keep theta with smallest loglike
					for(j in 1:dim(theta_test)[2]) {
						theta_j = theta_test[,j]
						attributes(theta_j) = attributes(thetaini)
						m2LL_test = m2LLa(theta_j, z = z, X = X, 
							xcoords = xcoords, ycoords = ycoords, 
							spatial_model = spatial_model, estMeth = estMeth, 
							grpindx = grpindx, max_range = max_range, 
							attr_theta = attrs, fixed_parms = fixed_parms,
							Zs = Zs)
						if (m2LL_test < min_m2LL) {
								thetaini = theta_j
								min_m2LL = m2LL_test
						}
					}
				}
				# estimate parameters with optim
				parmest <-optim(thetaini, m2LLa,  
					z = z, X = X, xcoords = xcoords, ycoords = ycoords, 
					spatial_model = spatial_model, estMeth = estMeth, 
					grpindx = grpindx, max_range = max_range, 
					attr_theta = attrs, fixed_parms = fixed_parms,
					Zs = Zs, ...)
				# get the minimized objective function value, equal to 
				# minus 2 times the log-likelihood
				m2LL <- parmest$value
				# parameter estimate at minimized log-likelihood
				theta_parms = parmest$par
			}
			# now apply the transformation functions, and scale, to get back to
			# natural scaling for parameters
			theta_opt = unlist(japply(theta_parms, attr(thetaini, "transfType")))*
							attr(thetaini,'scaleFactors')
			# attach attributes to parameter estimates
			attributes(theta_opt) = attrs
			# if optimizing for rotate parameter, find which quadrant it is in
			if('rotate' %in% attrs$label) {
				if(attr(m2LLa(theta = theta_parms, z = z,
					X = X, xcoords = xcoords, ycoords = ycoords, 
					spatial_model = spatial_model, estMeth = estMeth, 
					grpindx = grpindx, max_range = max_range, attr_theta = attrs,
					fixed_parms = fixed_parms, Zs = Zs), 
					'which') == 2)
						theta_opt[attr(theta_opt,'label') == 'rotate'] = 180 - 
							theta_opt[attr(theta_opt,'label') == 'rotate']
			}
			# if there are fixed parameters, concatenate with optimized ones
			# to create all parameter values, with attributes
			if(!is.null(fixed_parms)) {
				ind4opt = !attr(thetaini,'label') %in% attr(fixed_parms, 'label')
				theta = c(theta_opt,fixed_parms)
				attr(theta,'label') = c(attrs$label, attr(fixed_parms, 'label'))
				attr(theta,'scaleFactors') = c(attrs$scaleFactors, 
					rep(NA, times = length(fixed_parms)))
				attr(theta,'transfType') = c(attrs$transfType, 
					rep(NA, times = length(fixed_parms)))
				attr(theta,'fixed') = c(rep('No', times = length(theta_opt)), 
					rep('Yes', times = length(fixed_parms)))
			} else {
				theta = theta_opt
				attributes(theta) = attrs
				attr(theta,'fixed') = rep('No', times = length(theta))
			}
		}
	} else { # profile overall variance parameter if possible
		# create initial parameter estimate
		if(is.null(thetaini)) {
			# start with 2 spatial covariance parameters for correlation matrix
			# nugget will be leftover proportion of partial sill
			# (0 <= partial sill <= 1)
			tlength = 2
			label = c('partial sill', 'range')
			# scaling factors. Will expand parameters by these factors inside
			# optimization routine
			scaleFactors <- c(1, 2*
					sqrt((var(xcoords) + var(ycoords))/2))
			# transformatin functions.  Will transform parameters by these functions
			# inside optimization routine
			transfType = c(expit, exp)
			# if using anisotropy, add two parameters and scaling factors and
			# transformation functions
			if(use_anisotropy == TRUE) {
				tlength = tlength + 2
				label = c(label, 'minorp', 'rotate')
				scaleFactors = c(scaleFactors, 1, 90)
				transfType = c(transfType, expit, expit)
			}
			# check for models with extra parameters, and add it to parameter
			# initialization, with scaling factors and transformation functions
			if(any(as.character(match.call()$spatial_model) %in% c('cauchy', 'stable',
				'besselk', 'besselj'))) {
					tlength = tlength + 1
					label = c(label, 'extrap')
					scaleFactors = c(scaleFactors, 1)
					transfType = c(transfType, exp)
			}
			# check for random effects
			if(!is.null(random_formula)) {
				tlength = tlength + length(random_formula)
				label = c(label, names(random_formula))
				scaleFactors = c(scaleFactors, 
					rep(1, times = length(random_formula)))
				for(i in 1:length(random_formula))
					transfType = c(transfType, expit)
			}		
			# start with all parameters at 0 (essentially unscaled with an 
			# inverse transformation
			thetaini <- matrix(rep(0, times = tlength), ncol  = 1)
			rownames(thetaini) <- 1:length(thetaini)
			colnames(thetaini) <- "parms"
			attr(thetaini,'label') = label
			attr(thetaini,'scaleFactors') = scaleFactors
			attr(thetaini,'transfType') = transfType 
		} else {
		# if initial parameters are provided, inverse scale and transform
			sill = 0
			pm = 1
			vec = NULL
			label = NULL
			scaleFactors = NULL
			transfType = NULL
			if(any(attr(thetaini,'label') == 'partial sill')) sill = sill +
				thetaini[attr(thetaini,'label') == 'partial sill']
			if(any(attr(thetaini,'label') == 'nugget')) sill = sill +
				thetaini[attr(thetaini,'label') == 'nugget']
			if(!is.null(Zs)) {
				namesZ = names(Zs)
				for(j in 1:length(namesZ)) 
					sill = sill + thetaini[attr(thetaini,'label') == namesZ[j]]
			}
			if(any(attr(thetaini,'label') == 'partial sill'))
			{				
				val = logit(thetaini[attr(thetaini,'label') == 
					'partial sill']/(pm*sill))
				vec = c(vec, val)
				label = c(label, 'partial sill')
				scaleFactors = c(scaleFactors, 1)
				transfType = c(transfType, expit)
				pm = pm*(1 - expit(val))
			}	
			if(!is.null(Zs)) {
				namesZ = names(Zs)
				for(j in 1:length(namesZ)) {
					val = logit(thetaini[attr(thetaini,'label') == 
						namesZ[j]]/(sill*pm))
					vec = c(vec, val)
					pm = pm*(1 - expit(val))
					label = c(label, namesZ[j])
					scaleFactors = c(scaleFactors, 1)
					transfType = c(transfType, expit)
				}
			}
			if(any(attr(thetaini,'label') == 'range')) {
				vec = c(vec, 
				log(thetaini[attr(thetaini,'label') == 'range']/
					(sqrt((var(xcoords) + var(ycoords))/2))))
				label = c(label, 'range')
				scaleFactors = c(scaleFactors, 
					sqrt((var(xcoords) + var(ycoords))/2))
				transfType = c(transfType, exp)
			}
			if(any(attr(thetaini,'label') == 'minorp')) {
				vec = c(vec, logit(thetaini[attr(thetaini,'label') == 'minorp']))
				label = c(label, 'minorp')
				scaleFactors = c(scaleFactors, 1)
				transfType = c(transfType, expit)
			}
			if(any(attr(thetaini,'label') == 'rotate')) {
				vec = c(vec, 
					logit(min(thetaini[attr(thetaini,'label') == 'rotate'],
					(180 - thetaini[attr(thetaini,'label') == 'rotate']))/90))
				label = c(label, 'rotate')
				scaleFactors = c(scaleFactors, 90)
				transfType = c(transfType, expit)
			}
			if(any(attr(thetaini,'label') == 'extrap')) {
				vec = c(vec, log(thetaini[attr(thetaini,'label') == 'extrap']))
				label = c(label, 'extrap')
				scaleFactors = c(scaleFactors, 1)
				transfType = c(transfType, exp)
			}
			thetaini = vec
      attr(thetaini,'label') = label
      attr(thetaini, 'scaleFactors') = scaleFactors
      attr(thetaini, 'transfType') = transfType
		}
		# set maximum range value if none given
		if(is.null(max_range)) max_range = 10*
					sqrt((var(xcoords) + var(ycoords))/2)
		# if some parameters are fixed, set those here and subset the
		# initial parameter estimate
		if(!is.null(fixed_parms)) {
			ind4opt = !attr(thetaini,'label') %in% attr(fixed_parms, 'label')
			label = attr(thetaini, 'label')
			scaleFactors = attr(thetaini, 'scaleFactors')
			transfType = attr(thetaini, 'transfType')
			thetaini = thetaini[ind4opt]
			attr(thetaini,'label') = label[ind4opt]
			attr(thetaini,'scaleFactors') = scaleFactors[ind4opt]
			attr(thetaini,'transfType') = transfType[ind4opt] 
		}
		# store the attributes to pass to optimization function
		attrs = attributes(thetaini)
    # estimate parameters
    # if only 1 optimization parameter, use optimize()
		if(length(thetaini) == 1) {
			# estimate parameters
			parmest = optimize(m2LLprof, lower = -30, upper = 30,
				z = z, X = X, xcoords = xcoords, ycoords = ycoords, 
				spatial_model = spatial_model, estMeth = estMeth, 
				grpindx = grpindx, max_range = max_range, 
				attr_theta = attrs, fixed_parms = fixed_parms,
				Zs = Zs)
			# get the minimized objective function value, equal to 
			# minus 2 times the log-likelihood
			m2LL <- parmest$objective
			# parameter estimate at minimized log-likelihood
			theta_parms = parmest$minimum
			attributes(theta_parms) = attrs
			# if multiple parameters, use optim()
			} else {
				# do a search grid to get good starting values
				if(use_grid == TRUE & length(thetaini) > 1 & length(thetaini) < 9) {
					if(length(thetaini) == 2) gbase = 7
					if(length(thetaini) == 3) gbase = 5
					if(length(thetaini) %in% 4:5) gbase = 3
					if(length(thetaini) > 5) gbase = 2
					min_m2LL = m2LLprof(thetaini, z = z, X = X, 
							xcoords = xcoords, ycoords = ycoords, 
							spatial_model = spatial_model, estMeth = estMeth, 
							grpindx = grpindx, max_range = max_range, 
							attr_theta = attrs, fixed_parms = fixed_parms,
							Zs = Zs)
					theta_test = NULL	
					for(j in 0:(length(thetaini)-1))	
            theta_test = cbind(theta_test,
              rep(1, times = 2^j) %x% 0:(gbase-1)  %x% 
                rep(1, times = 2^(length(thetaini)-j-1)))
          theta_test = 2*(theta_test/(gbase-1) - 0.5)
					# theta_test is centered on zero.  Now, center it on thetaini
					theta_test = as.vector(thetaini) + t(theta_test)
					for(j in 1:dim(theta_test)[2]) {
						theta_j = theta_test[,j]
						attributes(theta_j) = attributes(thetaini)
						m2LL_test = m2LLprof(theta_j, z = z, X = X, 
								xcoords = xcoords, ycoords = ycoords, 
								spatial_model = spatial_model, estMeth = estMeth, 
								grpindx = grpindx, max_range = max_range, 
								attr_theta = attrs, fixed_parms = fixed_parms,
								Zs = Zs)
						if (m2LL_test < min_m2LL) {
							thetaini = theta_j
							min_m2LL = m2LL_test
						}
					}
				}
				parmest <-optim(thetaini, m2LLprof,  
					z = z, X = X, xcoords = xcoords, ycoords = ycoords, 
					spatial_model = spatial_model, estMeth = estMeth, 
					grpindx = grpindx, max_range = max_range, 
					attr_theta = attrs, fixed_parms = fixed_parms,
					Zs = Zs, ...)
				# get the minimized objective function value, equal to 
				# minus 2 times the log-likelihood
				m2LL <- parmest$value
				# parameter estimate at minimized log-likelihood
				theta_parms = parmest$par
		}
		theta_opt = unlist(japply(theta_parms, attr(theta_parms, "transfType")))*
			attr(theta_parms,'scaleFactors')
		attributes(theta_opt) = attrs
		# if optimizing for rotate parameter, find which quadrant it is in
		if('rotate' %in% attrs$label) {
			if(attr(m2LLprof(theta = theta_parms, z = z,
				X = X, xcoords = xcoords, ycoords = ycoords, 
				spatial_model = spatial_model, estMeth = estMeth, 
				grpindx = grpindx, max_range = max_range, attr_theta = attrs,
				fixed_parms = fixed_parms, Zs = Zs), 
				'which') == 2)
					theta_opt[attr(theta_opt,'label') == 'rotate'] = 180 - 
						theta_opt[attr(theta_opt,'label') == 'rotate']
		}
		if(!is.null(fixed_parms)) {
			theta_opt = c(theta_opt, fixed_parms)
			attr(theta_opt,'label') = c(attrs$label, attr(fixed_parms,'label'))
			attr(theta_opt,'fixed') = c(rep('No', 
				times = length(theta_opt) - length(fixed_parms)), 
				rep('Yes', times = length(fixed_parms)))
			attrs = attributes(theta_opt)
		} else {
      attr(theta_opt,'fixed') = c(rep('No', 
				times = length(theta_opt)))
    }
		if(!any('rotate' %in% attr(theta_opt,'label'))) {
			rotate = 90
			minorp = 1
		} else {
			rotate = theta_opt[attr(theta_opt,'label') == 'rotate']
			minorp = theta_opt[attr(theta_opt,'label') == 'minorp']
		}
		extrap = NULL
		if(any('extrap' %in% attr(theta_opt,'label'))) 
			extrap = theta_opt[attr(theta_opt,'label') == 'extrap']
		parsil = theta_opt[attr(theta_opt,'label') == 'partial sill']					
		if(is.null(Zs)) {
			nugget = 1 - parsil
		} else {
			namesZ = names(Zs)
			ps = c(1 - parsil, (1 - parsil)*
				cumprod(1 - theta_opt[attr(theta_opt,'label') %in% namesZ]))*
				c(theta_opt[attr(theta_opt,'label') %in% namesZ],1)
			vcs = ps[1:length(namesZ)]
			nugget = ps[length(ps)]
			theta_opt[attr(theta_opt,'label') %in% namesZ] = vcs
		}
		theta = c(theta_opt, nugget)
		attr(theta,'label') = c(attr(theta_opt, 'label'), 'nugget')
		attr(theta,'fixed') = c(attr(theta_opt, 'fixed'), 'no')
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
				range = theta[attr(theta_opt,'label') == 'range'], 
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
		if(estMeth == 'ML') sig2 = as.numeric(rVir/n)
		if(estMeth == 'REML') sig2 = as.numeric(rVir/(n - p))
		theta[attr(theta,'label') == 'partial sill'] = 
			theta[attr(theta,'label') == 'partial sill']*sig2
		if(!is.null(Zs)) {
			theta[attr(theta,'label') %in% namesZ] = 
				theta[attr(theta,'label') %in% namesZ]*sig2
		}
	  theta[attr(theta,'label') == 'nugget'] = 
			theta[attr(theta,'label') == 'nugget']*sig2
	}
	# if the user chooses parallel processing...
  if(par == TRUE) {
	parmest <-optim(log(thetaini), m2LLgPAR, z = z, X = X, 
			xcoords = xcoords, ycoords = ycoords, 
			estMeth = 'REML', grpindx = grpindx)
	}
	stoptime = Sys.time()

	# create output object
	out = list(call =  match.call(expand.dots = FALSE),
		theta = theta, m2LL = m2LL, z = z, X = X, Zs = Zs,
		gi = grpindx, coords = coords, formula = formula,
		respCol = respCol, spatial_model = spatial_model,
		spatial_name = as.character(match.call()$spatial_model),
		cope_est_meth = estMeth,
		fixed_parms = fixed_parms,
		theta_ini = thetaini, max_range = max_range,
		parmest = parmest, profile_sigma = profile_sigma,
		opt_time = stoptime - starttime
	)
	# give it a class
	class(out) = 'slm_cope'
	# return
	return(out)
}
