#' Creates a systematic grid of points
#'
#' Creates a systematic grid of points 
#'
#' @param lower.x.lim the lower limit for x-coordinate
#' @param upper.x.lim the upper limit for x-coordinate
#' @param lower.y.lim the lower limit for y-coordinate
#' @param upper.y.lim the upper limit for y-coordinate
#' @param ncol the number of cols in the systematic grid
#' @param nrow the number of rows in the systematic grid

#' @return A data.frame with x- and y-coordinates of simulated locations
#' @author Jay Ver Hoef \email{jay.verhoef@@noaa.gov}
#' @export

pointSimSyst <- function (nrow = 10, ncol = 10, lower.x.lim = -1, upper.x.lim = 1, 
    lower.y.lim = -1, upper.y.lim = 1) 
{
    x.range <- upper.x.lim - lower.x.lim
    y.range <- upper.y.lim - lower.y.lim
    y.mat <- lower.y.lim + y.range * (nrow - matrix(rep(1:nrow, 
        times = ncol), nrow = nrow, ncol = ncol))/(nrow) + y.range/(2 * 
        nrow)
    x.mat <- lower.x.lim + x.range * (t(matrix(rep(1:ncol, times = nrow), 
        nrow = ncol, ncol = nrow)) - 1)/(ncol) + x.range/(2 * 
        ncol)
    data.frame(x = matrix(x.mat, ncol = 1), y = matrix(y.mat, 
        ncol = 1))
}

