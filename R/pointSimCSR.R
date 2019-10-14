#-------------------------------------------------------------------------------
#
#           pointSimCSR
#
#-------------------------------------------------------------------------------

#' Simulate clustered spatial point pattern
#'
#' simulates a clustered spatial point pattern
#'
#' @param npair number of points to add that are completely spatially random (CSR), default = 100 
#' @param lower.x.lim left limit of boundary, default = -1
#' @param upper.x.lim right limit of boundary, default = 1
#' @param lower.y.lim lower limit of boundary, default = -1
#' @param upper.y.lim upper limit of boundary, default = 1
#'
#' @return data.frame of two columns, x-coordinate in the first, and y-coordinate in the second
#'
#' @author Jay Ver Hoef
#' @export
pointSimCSR <-
function(npair = 100, lower.x.lim = -1, upper.x.lim = 1,
	lower.y.lim = -1, upper.y.lim = 1)
{
	x.range <- upper.x.lim - lower.x.lim
	y.range <- upper.y.lim - lower.y.lim
	data.frame(x=lower.x.lim + runif(npair)*x.range,
		y=lower.y.lim + runif(npair)*y.range)
}

