#' @title calculate  weights for a bandpass filter
#' @description Derive the (smoothed) least square highpass based on Bloomfield 1976
#' @param omega.c cutoff frequency
#' @param n length of the filter, has to be odd
#' @param sample sampling rate of the timeseries on which the filter will be applied (1/deltat)
#' @param convergence TRUE: smoothed least square lowpass; FALSE = unsmoothed
#' @param omega.c cutoff frequency
#' @return vector of filter weights
#' @author Thomas Laepple
#' @export

Highpass<-function(omega.c,n=9,sample=1,convergence=T)
{
    if ((n %% 2) == 0) stop("N must be odd, this function calculates only symetrical = phase preserving filters")

	directtransfer<-rep(0,n)
	directtransfer[(n+1)/2]<-1
	return(directtransfer-Lowpass(omega.c,n=n,sample=sample,convergence=convergence))
}

