#' @title lowpass filtered expected correlation of powerlaw signal pair
#' @description Correlation of two timeseries with powerlaw signal and powerlaw
#' noise evaluated until f
#' this equals the correlation of linearly coupled lowpass filtered powerlaw
#' process
#' @param f frequency until which to integrate
#' @param betaSignal powerlaw slope of the signal
#' @param betaNoise powerlaw slope of the noise
#' @param N Number of points per timeseries
#' @param r expected correlation of the unfiltered
#' @return expected correlation of the timeseries
#' @export
#' @examples
#' temp <- replicate(1000,PSP.CorAfterRollmean(1000,1,0,0.5))
#' rowMeans(temp)
#' PSP.CorUntilF(f=c(0.5,0.5/10,0.5/50),betaSignal=1,betaNoise=0,N=1000,r=0.5)
#' @author Thomas Laepple
PSP.CorUntilF <- function(f, betaSignal, betaNoise, N,r) {

    scale.var.noise <- ((1-r^2)/(r^2))
    var.signal <- PS.VarUntilF(f, betaSignal, N)/ PS.VarUntilF(0.5,
        betaSignal, N)
    var.noise <- PS.VarUntilF(f, betaNoise, N)/ PS.VarUntilF(0.5,
        betaNoise, N)
    sqrt(var.signal/(var.signal + var.noise * scale.var.noise))
}
