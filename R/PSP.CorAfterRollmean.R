#' @title Numerical correlation of random timeseries with different filtering (running mean)
#' @param N Number of points of the timeseries
#' @param betaSignal powerlaw slope of the signal
#' @param betaNoise powerlaw slope of the noise
#' @param R expected correlation of the whole timeseries
#' @return correlation of unfiltered, 10point mean and 50 point mean
#' @export
#' @examples
#' temp <- replicate(1000,PSP.CorAfterRollmean(1000,1,0,0.5))
#' rowMeans(temp)
#' PSP.CorUntilF(c(0.5,0.5/10,0.5/50),1,0,1000,0.5)
#' @author Thomas Laepple
#'
PSP.CorAfterRollmean <- function(N, betaSignal, betaNoise, R)
{
    x <- SimulatePowerlawSignalPair(N, betaSignal, betaNoise, R)
    c.1 <- cor(x$y1, x$y2)
    c.10 <- cor(rollmean(x$y1, 10), rollmean(x$y2, 10))
    c.50 <- cor(rollmean(x$y1, 50), rollmean(x$y2, 50))
    return(c(c.1, c.10, c.50))
}
