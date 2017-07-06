##' @title A PSD(freq) for a powerlaw with variance 1
##' @param beta slope of the powerlaw
##' @param freq frequency vector
##' @return vector containing the PSD
##' @author Thomas Laepple
##' @export
AnPowerlaw<-function(beta,freq,return.scaling=FALSE)
{
  #power law PSD. variance 1, slope beta, on the frequencies freq
       power=1/(freq^beta);
       powerScale<-sum(power)*mean(diff(freq))*2
       return(power/powerScale)
   }


##' Simulate a random timeseries with a powerlaw spectrum
##'
##' Method: FFT white noise, rescale, FFT back, the result is scaled to variance 1
##' @title Simulate a random timeseries with a powerlaw spectrum
##' @param beta slope
##' @param N length of timeseries to be generated
##' @return vector containing the timeseries
##' @author Thomas Laepple
##' @export
SimPowerlaw <- function(beta, N)
{
  N2 <- (3^ceiling(log(N, base = 3)))
  df  <- 1 / N2
  f <- seq(from = df, to = 1/2, by = df)
  Filter <- sqrt(1/(f^beta))
  Filter <- c(max(Filter), Filter, rev(Filter))
  x   <- scale(rnorm(N2, 1))
  fx  <- fft(x)
  ffx <- fx * Filter
  result <- Re(fft(ffx, inverse = TRUE))[1:N]
  return(scale(result))
}

