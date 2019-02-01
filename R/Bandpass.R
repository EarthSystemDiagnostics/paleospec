
#' @title calculate  weights for a bandpass filter
#' @description Derive the (smoothed) least square bandpass based on Bloomfield 1976
#' @param omega.upper upper cutoff frequency
#' @param omega.lower lower cutoff frequency
#' @param n length of the filter, has to be odd
#' @param sample sampling rate of the timeseries on which the filter will be applied (1/deltat)
#' @param convergence TRUE: smoothed least square lowpass; FALSE = unsmoothed
#' @param omega.c cutoff frequency
#' @return vector of filter weights
#' @author Thomas Laepple
#' @export
Bandpass<-function(omega.upper,omega.lower,n,sample=1,convergence=T)
{
if ((n %% 2) == 0) stop("N must be odd, this function calculates only symetrical = phase preserving filters")


	return(lowpass(omega.upper,n,sample,convergence)-lowpass(omega.lower,n,sample,convergence))
}
