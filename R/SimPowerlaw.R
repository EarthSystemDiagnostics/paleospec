#' @title A PSD(freq) for a powerlaw with variance 1
#' @param beta slope of the powerlaw
#' @param freq frequency vector
#' @return vector containing the PSD
#' @author Thomas Laepple
#' @export
AnPowerlaw<-function(beta,freq,return.scaling=FALSE)
{
  #power law PSD. variance 1, slope beta, on the frequencies freq
       power=1/(freq^beta);
       powerScale<-sum(power)*mean(diff(freq))*2
       return(power/powerScale)
   }



#' Simulate a random timeseries with a powerlaw spectrum
#'
#' Method: FFT white noise, rescale, FFT back, the result is scaled to variance 1
#' @title Simulate a random timeseries with a powerlaw spectrum
#' @param beta slope
#' @param N length of timeseries to be generated
#' @return vector containing the timeseries
#' @author Thomas Laepple
#' @export
#' @family SimPowerlaw SimPLS SimFromEmpiricalSpectrum
SimPowerlaw <- function(beta, N)
{
  N2 <- (3^ceiling(log(N, base = 3)))
  df  <- 1 / N2
  f <- seq(from = df, to = 1/2, by = df)
  Filter <- sqrt(1/(f^beta))
  Filter <- c(max(Filter), Filter, rev(Filter))
  x   <- rnorm(N2, 1)
  fx  <- fft(x)
  ffx <- fx * Filter
  result <- Re(fft(ffx, inverse = TRUE))[1:N]
  return(scale(result)[1:N])
}


#' Simulate a random timeseries with a powerlaw spectrum
#'
#' @param beta Slope of the powerlaw. beta = 1 produces timeseries with -1 slope
#'   when plotted on log-log power ~ frequency axes
#' @param N length of timeseries to be generated
#' @param alpha the constant. If alpha > 0 this is the parameter alpha * f^(-beta).
#'   If alpha < 0, the variance of the returned timeseries is scaled so that its
#'   expected value is abs(alpha)
#' @author Torben Kunz
#' @return a vector containing the timeseries
#' @description This function creates a power-law series. It has the problem
#'   that it effectively produces (fractional) Brownian bridges, that is, the
#'   end is close to the beginning (cyclic signal), rather than true fBm or fGn
#'   series.
#'
#'   If alpha>0, then the EXPECTED PSD is equal to alpha*f^(-beta).
#'
#'   If alpha<0, then the timeseries is normalized such that it has EXPECTED
#'   variance abs(alpha), and the EXPECTED PSD is proportional to f^(-beta).
#' @family SimPLS SimPowerlaw SimFromEmpiricalSpectrum
#' @export
#'
#' @examples
#' # With a beta = 1 and alpha = 0.1
#' set.seed(202010312)
#' ts1 <- ts(SimPLS(N = 1000, beta = 1, alpha = 0.1))
#' plot(ts1)
#' sp1 <- SpecMTM(ts1)
#' LPlot(sp1)
#' abline(log10(0.1), -1, col = "Red")
#'
#' # beta = 0.5, alpha = 0.4
#' ts2 <- ts(SimPLS(1000, beta = 0.5, alpha = 0.4))
#' plot(ts2)
#' sp2 <- SpecMTM(ts2)
#' LPlot(sp2)
#' abline(log10(0.4), -0.5, col = "Red")
#'
#' # beta = 1, alpha = -2
#' ts3 <- ts(SimPLS(1000, 1, alpha = -2))
#' plot(ts3)
#' var(ts3)
#'
#' # the EXPECTED variance is -2, for a given random timeseries the actual value will differ
#' rep.var <- replicate(100, {
#'  var(SimPLS(1000, 1, -2))
#' })

#' hist(rep.var)
#' abline(v  = 2, col = "Red")
#' mean(rep.var)
SimPLS <- function(N, beta, alpha = -1){

  frequency.axis <- function(N){
    fax <- 0:(N-1) / N
    fax[fax>0.5] <- fax[fax>0.5] - 1
    fax
  }

  # Pad the length of the timeseries so that it is highly composite - this speeds
  # up the FFT operations.
  N2 <- (3^ceiling(log(N, base = 3)))

  x2 <- rnorm(N2)
  xfft <- fft(x2)

  fax <- frequency.axis(N2)

  P <- c(0, abs(alpha) * abs(fax[2:N2])^(-beta))

  if (alpha<0) P <- P * abs(alpha) * (N2-1) / sum(P[1:N2])

  xfft <- xfft * sqrt(P)
  x2 <- fft(xfft, inverse=TRUE) / N2
  x2 <- Re(x2)

  # trim the timeseries back to the requested length
  x <- x2[1:N]

  # scale the variance of timeseries at requested length
  if (alpha<0) {
    sdx2 <- sd(x2)
    x <- sdx2 * x / sd(x)
  }

  return(x)
}



#' @title Simulate a random timeseries consistent with an arbitrary numerical power spectrum
#' @param spec Numerical power spectrum consisting of a list with components $freq and $spec
#' @param N length of timeseries to be generated
#' @description Adapted from SimPowerlaw
#' @return vector containing the timeseries
#' @author Thomas Laepple and Andrew Dolman
#' @family SimPowerlaw SimPLS SimFromEmpiricalSpectrum
#' @export
#' @examples
#' # Create a piecewise spectrum
#'
#' ## helper function to generate continuous piecewise spectrum
#'
#' PiecewiseLinear <- function(x, val.at.min.x, breaks, slopes){
#'
#' breaks <- c(-Inf, breaks, Inf)
#' slp.vec <- slopes[findInterval(x, breaks)]
#' d.x <- diff(x)
#' d.y <- c(d.x * tail(slp.vec, -1))
#'
#' y <- cumsum(c(val.at.min.x, d.y))
#'
#' data.frame(x, y)
#'
#' }
#' slps <- c(-1, -0.5, -1)
#' brks <- c(1e-03, 1e-02)
#' emp.spec <- PiecewiseLinear(log(seq(1/1e05, 1/2, 1/1e05)), 0, log(brks), slps)
#' emp.spec <- exp(emp.spec)
#' names(emp.spec) <- c("freq", "spec")
#'
#' plot(emp.spec, type = "l", log = "xy")
#'
#' # Sample consistent with spectrum
#' ts1 <- ts(SimFromEmpiricalSpec(emp.spec, 50000))
#'
#' # re-estimate power spectrum
#' spec1 <- SpecMTM(ts1)
#' LPlot(spec1)
#' lines(emp.spec, col = "Red")
#' abline(v = brks, col = "Green")
SimFromEmpiricalSpec <- function(spec, N)
{
  N2 <- (3^ceiling(log(N, base = 3)))
  df  <- 1 / N2
  f <- seq(from = df, to = 1/2, by = df)

  if (max(f) > max(spec$freq))
    stop("Need higher maximum frequency in input spectrum")

  if (min(f) < min(spec$freq))
    stop("Need lower minimum frequency in input spectrum")

  Filter <- sqrt(approx(spec$freq, spec$spec, f)$y)

  Filter <- c(max(Filter), Filter, rev(Filter))
  x   <- rnorm(N2, 1)
  fx  <- fft(x)
  ffx <- fx * Filter
  result <- Re(fft(ffx, inverse = TRUE))[1:N]

  # scale variance
  # integrate spectrum to 1/N
  t.var <- sqrt(2*sum(spec$spec[spec$freq >= 1/(2*N)] * abs(diff(spec$freq)[1])))

  out <- scale(result) * t.var
  out <- as.vector(out)

  return(out)
}




