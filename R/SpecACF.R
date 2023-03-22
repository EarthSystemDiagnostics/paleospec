#' Estimate Power Spectra via the Autocovariance Function
#'
#' @param x a vector or matrix of binned values, possibly with gaps
#' @param bin.width the width of the bins, effectively delta_t
#' @param pos.f.only return only positive frequencies, defaults to TRUE If TRUE,
#'   freq == 0, and frequencies higher than 1/(2*bin.width) which correspond to
#'   the negative frequencies are removed
#' @param demean remove the mean from each record (column) in x, defaults to
#'   TRUE. If detrend is TRUE, mean will be removed during detrending regardless
#'   of the value of demean
#' @param detrend remove the mean and any linear trend from each record (column)
#'   in x, defaults to FALSE
#' @description Estimates the power spectrum from a single time series, or the
#'   mean spectrum of a set of timeseries stored as the columns of a matrix.
#'   Timeseries can contain (some) gaps coded as NA values. Gaps results in
#'   additional estimation error so that the power estimates are no longer
#'   chi-square distributed and can contain additional additive error, to the
#'   extent that power at some frequencies can be negative. We do not have a
#'   full understanding of this estimation uncertainty, but simulation testing
#'   indicates that the estimates are unbiased such that smoothing across
#'   frequencies to remove negative estimates results in an unbiased power
#'   spectrum.
#' @author Torben Kunz and Andrew Dolman <andrew.dolman@awi.de>
#' @return a spec object (list)
#' @export
#'
#' @examples
#' set.seed(20230312)
#' x <- cumsum(rgamma(200, shape = 1.5, rate = 1.5/10))
#' y <- SimProxySeries(a = 0.1, b = 1, t.smpl = x, nt = 2000,
#'  smth.lab = list(type = "rect", tau = 1))
#' y_binned <- BinTimeseries(x, y, bin.width = 15)
#' sp1 <- SpecACF(y_binned$mean.value, bin.width = 15)
#' sp2 <- LogSmooth(sp1)
#' LPlot(sp1)
#' LLines(sp2, col = "red")
SpecACF <- function(x, bin.width,
                    demean = TRUE, detrend = TRUE,
                    TrimNA = TRUE,
                    pos.f.only = TRUE,
                    return.working = FALSE){

   # Convert ts to vector

  if (is.ts(x)){
    d <- dim(x)
    x <- as.vector(x)

    if (is.null(d) == FALSE){
      dim(x) <- d
    }

  }


  # Convert vector to matrix
  if (is.vector(x)){
    x <- matrix(x, ncol = 1)
  }



  if (is.data.frame(x) == TRUE){
    x <- as.matrix(x)
  }

  if (TrimNA){
    x <- TrimNA(x)
  }

  if (detrend){
    i <- seq_along(x[,1])
    x <- apply(x, 2, function(y) {
      stats::residuals(stats::lm(y ~ i, na.action = "na.exclude"))
    })
  }

  # remove mean from each record
  if (demean){
    x <- x - colMeans(x, na.rm = TRUE)
  }

  lag = (0:(nrow(x)-1))

  if (return.working == TRUE){
    working <- mean.ACF(x, return.working = TRUE)
    acf <- working[["acf"]]
  } else{
    acf <- mean.ACF(x)
  }

  freq = (lag / (nrow(x))) / bin.width
  spec = Re(stats::fft(acf)) * bin.width

  rfreq <- (utils::tail(freq, 1) + freq[2])/2

  if (pos.f.only){
    spec <- spec[freq > 0 & freq <= 1/(2*bin.width)]
    freq <- freq[freq > 0 & freq <= 1/(2*bin.width)]
  }

  out <- list(bin.width = bin.width,
              rfreq = rfreq,
              nrec = ncol(x),
              lag = lag, acf = acf,
              freq = freq, spec = spec, f.length = 1,
              dof = rep(2*ncol(x), length(freq)))

  class(out) <-  c("SpecACF", "spec")

  if (return.working){
    out <- list(working = working, spec = out)
  }

  return(out)
}


#' Bin a Timeseries Preserving Empty Bins
#'
#' @param time.var time variable, a vector
#' @param value.var value variable, a vector
#' @param bin.width with of time bins in output
#' @param strt.time first time value in output, defaults to min(time.var) - bin.width/2
#' @param end.time last time values in output, defaults to max(time.var) + bin.width/2
#' @description Convert a continuous time series to a discrete (regular)
#'   timeseries by binning, preserving any empty bins. Additionally returns the
#'   number of data points in each bin
#' @return a dataframe
#' @export
#' @author Andrew Dolman <andrew.dolman@awi.de>
#' @examples
#' set.seed(20230312)
#' x <- cumsum(rgamma(20, shape = 1.5, rate = 1.5/10))
#' y <- SimProxySeries(a = 0.1, b = 1, t.smpl = x)
#' BinTimeseries(x, y, bin.width = 15)
BinTimeseries <- function(time.var, value.var, bin.width,
                          strt.time = NULL, end.time = NULL){

  if (is.null(strt.time)){
    strt.time <-  round((min(time.var) - bin.width/2))
  }

  if (is.null(end.time)){
    end.time <-   round((max(time.var) + bin.width/2))
  }

  breaks <- seq(strt.time, end.time + bin.width/2, by = bin.width)

  time <- seq(strt.time + bin.width/2, end.time, by = bin.width)

  f <- cut(time.var, breaks = breaks)

  mean.value <- tapply(value.var, f, mean)
  mean.time <- tapply(time.var, f, mean)
  n.bin <- tapply(value.var, f, length)

  data.frame(time, bin = levels(f), mean.time, mean.value, n.bin, bin.width)
}


#' Title
#'
#' @param x a matrix
#' @param return.working return intermediate steps, defaults to FALSE
#'
#' @return a matrix
#' @keywords internal
mean.ACF <- function(x, return.working = FALSE){
  # Arguments:
  # x is a matrix, each column of which contains a timeseries (with missing values indicated by NA)

  # Value:
  # the circular ACF, with length(ACF)==dim(x)[1], being the length of the timeseries

  ix <- !is.na(x)
  x[is.na(x)] <- 0

  mvacfx <- mvacf.by.fft(x)
  mvacfix <- mvacf.by.fft(ix)

  sum.x.mvacf <- rowSums(mvacfx)
  sum.ix.mvacf <- rowSums(mvacfix)

  if (min(round(sum.ix.mvacf * dim(x)[1]))==0){
    warning("Some lags cannot be computed!")
    sum.ix.mvacf[round(sum.ix.mvacf * dim(x)[1])==0] <- NA
  }
  # This means the ACF could not be computed for at least one lag.
  # One could set ix.mvacf[round(ix.mvacf)==0] <- NA
  # but from an ACF with NAs one cannot compute the PSD via the FFT
  # unless one first interpolates the ACF across the missing lags.

  if (return.working == TRUE){
    out <- list(x = x, ix = ix, mvacfx = mvacfx, mvacfix = mvacfix,
                sum.x.mvacf = sum.x.mvacf, sum.ix.mvacf = sum.ix.mvacf,
                acf = sum.x.mvacf / sum.ix.mvacf)
  } else {
    out <- sum.x.mvacf / sum.ix.mvacf
  }

  return(out)

}


#' Title
#'
#' @param x a matrix
#'
#' @return a matrix
#' @keywords internal
mvacf.by.fft <- function(x){
  x.fft <- stats::mvfft(x)
  (Re(stats::mvfft(x.fft * Conj(x.fft), inverse=TRUE))) / dim(x)[1]^2
}

#' Remove leading and trailing rows of all NA
#'
#' @param m a numeric matrix
#'
#' @return a numeric matrix
#' @keywords internal
#' @examples
#' m <- matrix(c(NA, NA, NA, 1:9, NA,NA,NA, 10:12, NA,NA,NA), ncol = 3, byrow = TRUE)
#' m
#' PaleoSpec:::TrimNA(m)
TrimNA <- function(m){

  empty.row <- is.nan(rowMeans(m, na.rm = TRUE))
  rank.good <- (empty.row == FALSE) * 1:length(empty.row)

  first.good <- which.min(empty.row * 1:length(empty.row))
  last.good <- which.max(rank.good)

  m[first.good:last.good, , drop = FALSE]

}


