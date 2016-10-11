##' MTM spectral estimator
##' calls spec.mtm from library multitaper
##'  see spec.mtm from library mulitaper ?spec.mtm
##' @title MTM spectral estimator
##' @param timeSeries A time series of equally spaced data, this can be created
##' by the ts() function where deltat is specified.
##' @param k a positive integer, the number of tapers, often 2*nw.
##' @param nw  a positive double precision number, the time-bandwidth
##'          parameter.
##' @param nFFT  This function pads the data before computing the fft. nFFT
##'          indicates the total length of the data after padding.
##' @param centre 
##' @param dpssIN 
##' @param returnZeroFreq 
##' @param Ftest 
##' @param jackknife 
##' @param jkCIProb 
##' @param maxAdaptiveIterations 
##' @param plot 
##' @param na.action 
##' @param returnInternals 
##' @param detrend 
##' @param bPad 
##' @param ... 
##' @return spectra object list(freq,spec,dof)
##' examples
##' x<-ts(arima.sim(list(ar = 0.9),1000))
##' spec<-SpecMTM(x)
##' LPlot(spec,col="grey")
##' LLines(LogSmooth(spec),lwd=2)
##' @author Thomas Laepple
##' @export
SpecMTM<-function (timeSeries, k=3, nw=2, nFFT = "default", centre = c("Slepian"),
    dpssIN = NULL, returnZeroFreq = FALSE, Ftest = FALSE, jackknife = FALSE,
    jkCIProb = 0.95, maxAdaptiveIterations = 100, plot = FALSE,
    na.action = na.fail, returnInternals = FALSE,detrend=TRUE,bPad=FALSE, ...)
{

    if (sum(is.na(timeSeries))>0) stop("missing data")
    if (!bPad) nFFT=length(timeSeries)
    if (detrend) timeSeries[]<-lm(timeSeries~seq(timeSeries))$residuals
    result<-multitaper::spec.mtm(timeSeries=timeSeries, k=k, nw=nw, nFFT=nFFT, centre=centre,
    dpssIN=dpssIN, returnZeroFreq=returnZeroFreq, Ftest=Ftest, jackknife=jackknife,
    jkCIProb=jkCIProb,  maxAdaptiveIterations=maxAdaptiveIterations, plot=plot,
   na.action= na.action, returnInternals=returnInternals, ...)
    result$dof<-result$mtm$dofs
    return(result)

}
