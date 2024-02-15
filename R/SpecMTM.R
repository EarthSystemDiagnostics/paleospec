#' @title MTM spectral estimator
#' @description calls \code{\link[multitaper]{spec.mtm}} from library multitaper
#' @param timeSeries A time series of equally spaced data, this can be created
#' by the ts() function where deltat is specified.
#' @param k a positive integer, the number of tapers, often 2*nw.
#' @param nw  a positive double precision number, the time-bandwidth
#'          parameter.
#' @param nFFT  This function pads the data before computing the fft. nFFT
#'          indicates the total length of the data after padding.
#' @param detrend logical, detrend timeseries before estimating the spectrum
#' @param bPad if FALSE (the default) nFFT is set to the length of the timeseries
#' @inheritParams multitaper::spec.mtm
#' @param ... additional arguments to multitaper::spec.mtm
#' @return spectra object list(freq, spec, dof)
#' @family functions to estimate power spectra
#' @importFrom multitaper spec.mtm
#' @examples
#' x <- ts(arima.sim(list(ar = 0.9), 1000))
#' spec <- SpecMTM(x)
#' LPlot(spec, col='grey')
#' LLines(LogSmooth(spec), lwd=2)
#' @author Thomas Laepple
#' @export

SpecMTM <- function(timeSeries, k = 3, nw = 2, nFFT = "default",
  centre = c("Slepian"), dpssIN = NULL, returnZeroFreq = FALSE,
  Ftest = FALSE, jackknife = FALSE, jkCIProb = 0.95, maxAdaptiveIterations = 100,
  plot = FALSE, na.action = na.fail, returnInternals = FALSE,
  detrend = TRUE, bPad = FALSE, ...) {

  if (bPad != FALSE) warning("bPad doesn't do anything unless FALSE")

  if (sum(is.na(timeSeries)) > 0)
    stop("missing data")
  if (!bPad)
    nFFT = length(timeSeries)
  if (detrend)
    timeSeries[] <- lm(timeSeries ~ seq(timeSeries))$residuals

  result <- multitaper::spec.mtm(timeSeries = timeSeries, k = k,
    nw = nw, nFFT = nFFT, centre = centre, dpssIN = dpssIN,
    returnZeroFreq = returnZeroFreq, Ftest = Ftest, jackknife = jackknife,
    jkCIProb = jkCIProb, maxAdaptiveIterations = maxAdaptiveIterations,
    plot = plot, na.action = na.action, returnInternals = returnInternals,
    ...)
  result$dof <- result$mtm$dofs

  # if all dof are the same (e.g. when k = 1) then spec.mtm returns a vector of
  # length 1, this can cause issues with later functions, so replicate to length
  # of freq vector
  if (length(result$dof) == 1){
    result$dof <- rep(result$dof, length(result$freq))
  }

  return(result)
}



