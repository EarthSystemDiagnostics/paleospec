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

#' @title Simulate a random timeseries consistent with an arbitrary numerical power spectrum
#' @param spec Numerical power spectrum consisting of a list with components $freq and $spec
#' @param N length of timeseries to be generated
#' @description Adapted from SimPowerlaw
#' @return vector containing the timeseries
#' @author Thomas Laepple and Andrew Dolman
#' @export
#' @examples
#' # Create a piecewise spectrum
#'
#' ## helper function to generate continuous piecewise spectrum
#' new.scl <- function(f, prev.scl, new.slp, prev.slp){
#'   prev.scl*f^(prev.slp-new.slp)
#' }
#'
#' emp.spec <- data.frame(freq = seq(1/1e05, 1/2, 1/1e05))
#' emp.spec$spec <- with(emp.spec, {
#'   0.1 * freq ^ -1 * (freq > 1e-02) +
#'     new.scl(1e-02, 0.1,-2,-1) * freq ^ -2 * (freq <= 1e-02 & freq > 1e-03) +
#'     new.scl(1e-03, new.scl(1e-02, 0.1,-2,-1),-1,-2) * freq ^ -1 * (freq <= 1e-03)
#'
#' })
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
#' abline(v = c(1e-02, 1e-03), col = "Green")
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
  t.var <- sqrt(2*sum(spec$spec[spec$freq >= 1/N] * abs(diff(spec$freq)[1])))

  return(scale(result) * t.var)
}




