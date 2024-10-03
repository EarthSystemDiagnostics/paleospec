#' Estimate Power Spectra via the Autocovariance Function With Optional Slepian
#' Tapers
#'
#' @param x a vector or matrix of binned values, possibly with gaps
#' @param deltat,bin.width the time-step of the timeseries, equivalently the
#' width of the bins in a binned timeseries, set only one
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
#' @inheritParams SpecMTM
#' @importFrom multitaper spec.mtm
#' @author Torben Kunz and Andrew Dolman <andrew.dolman@awi.de>
#' @return a spec object (list)
#' @family functions to estimate power spectra
#' @export
#'
#' @examples
#' set.seed(20230312)
#'
#' # Comparison with SpecMTM
#'
#' tsM <- replicate(2, SimPLS(1e03, 1, 0.1))
#' spMk3 <- SpecACF(tsM, bin.width = 1, k = 3, nw = 2)
#' spMk1 <- SpecACF(tsM, bin.width = 1, k = 1, nw = 0)
#'
#' spMTMa <- SpecMTM(tsM[,1], deltat = 1)
#' spMTMb <- SpecMTM(tsM[,2], deltat = 1)
#' spMTM <- spMTMa
#' spMTM$spec <- (spMTMa$spec + spMTMb$spec)/2
#'
#' gg_spec(list(
#'   `ACF k=1` = spMk1,
#'   `ACF k=3` = spMk3,
#'   `MTM k=3` = spMTM
#' ), alpha.line = 0.75) +
#'   ggplot2::facet_wrap(~spec_id)
#'
#' ## No gaps
#'
#' ts1 <- SimPLS(1000, 1, 0.1)
#'
#' sp_ACF1 <- SpecACF(ts1, 1, k = 1)
#' sp_MTM7 <- SpecMTM(ts1, nw = 4, k = 7, deltat = 1)
#' sp_ACF7 <- SpecACF(ts1, 1, k = 7, nw = 4)
#'
#' gg_spec(list(
#'   `ACF k=1` = sp_ACF1, `ACF k=7` = sp_ACF7, `MTM k=7` = sp_MTM7
#' ))
#'
#' # With Gaps
#'
#' gaps <- (arima.sim(list(ar = 0.5), n = length(ts1))) > 1
#' table(gaps)
#' ts1_g <- ts1
#' ts1_g[gaps] <- NA
#'
#' sp_ACF1_g <- SpecACF(ts1_g, 1)
#' sp_ACFMTM1_g <- SpecACF(ts1_g, bin.width = 1, nw = 4, k = 7)
#'
#' gg_spec(list(
#'   ACF_g = sp_ACF1_g,
#'   ACF_g_smoothed = FilterSpecLog(sp_ACF1_g),
#'   ACF_g_tapered = sp_ACFMTM1_g
#' ), conf = FALSE) +
#'   ggplot2::geom_abline(intercept = log10(0.1), slope = -1, lty = 2)
#'
#'
#'
#' ## AR4
#' arc_spring <- c(2.7607, -3.8106, 2.6535, -0.9238)
#'
#' tsAR4 <- arima.sim(list(ar = arc_spring), n = 1e03) + rnorm(1e03, 0, 10)
#' plot(tsAR4)
#' spAR4_ACF <- SpecACF(tsAR4, 1)
#' spAR4_MTACF <- SpecACF(as.numeric(tsAR4), 1, k = 15, nw = 8)
#'
#' gg_spec(list(#'
#'   `ACF k=1` = spAR4_ACF,
#'   `ACF k=15` = spAR4_MTACF)
#' )
#'
#' ## Add gaps to timeseries
#'
#' gaps <- (arima.sim(list(ar = 0.5), n = length(tsAR4))) > 2
#' table(gaps)
#' tsAR4_g <- tsAR4
#' tsAR4_g[gaps] <- NA
#'
#' plot(tsAR4, col = "green")
#' lines(tsAR4_g, col = "blue")
#'
#' table(tsAR4_g > 0, useNA = "always")
#'
#' spAR4_ACF_g <- SpecACF(as.numeric(tsAR4_g), 1)
#' spAR4_MTACF_g <- SpecACF(as.numeric(tsAR4_g), 1, nw = 8, k = 15)
#'
#' table(spAR4_ACF_g$spec < 0)
#' table(spAR4_MTACF_g$spec < 0)
#'
#' gg_spec(list(
#'   `ACF gaps k=1` = spAR4_ACF_g,
#'   `ACF gaps k = 15` = spAR4_MTACF_g,
#'   `ACF full k = 15` = spAR4_MTACF
#' )
#' )
SpecACF <- function(x,
                      deltat = NULL, bin.width = NULL,
                      k = 1, nw = 0,
                      demean = TRUE, detrend = TRUE,
                      TrimNA = TRUE,
                      pos.f.only = TRUE,
                      return.working = FALSE) {
  if (is.null(deltat) & is.null(bin.width) & is.ts(x) == FALSE) {
    stop("One of deltat or bin.width must be set")
  }

  # Convert ts to vector
  if (is.ts(x)) {
    d <- dim(x)
    dt_ts <- deltat(x)

    x <- as.vector(x)

    if (is.null(deltat) == FALSE) {
      if (dt_ts != deltat) stop("timeseries deltat does not match argument deltat")
    }

    if (is.null(bin.width) == FALSE) {
      if (dt_ts != bin.width) stop("timeseries deltat does not match argument bin.width")
    }

    if (is.null(deltat)) deltat <- dt_ts
    if (is.null(bin.width)) bin.width <- dt_ts

    if (is.null(d) == FALSE) {
      dim(x) <- d
    }
  }

  if (is.null(bin.width) & is.null(deltat) == FALSE) {
    bin.width <- deltat
  }

  if (is.null(bin.width) == FALSE & is.null(deltat)) {
    deltat <- bin.width
  }

  # Convert vector to matrix
  if (is.vector(x)) {
    x <- matrix(x, ncol = 1)
  }

  if (is.data.frame(x) == TRUE) {
    x <- as.matrix(x)
  }

  if (TrimNA) {
    x <- TrimNA(x)
  }

  if (detrend) {
    i <- seq_along(x[, 1])
    x <- apply(x, 2, function(y) {
      stats::residuals(stats::lm(y ~ i, na.action = "na.exclude"))
    })
  }

  # remove mean from each record
  if (demean) {
    x <- x - colMeans(x, na.rm = TRUE)
  }

  lag <- (0:(nrow(x) - 1))

  ncolx <- ncol(x)

  ## Tapering
  if (k > 1) {
    n <- nrow(x)

    dpssIN <- multitaper:::dpss(n,
                                k = k, nw = nw,
                                returnEigenvalues = TRUE
    )
    dw <- dpssIN$v #* sqrt(bin.width)
    ev <- dpssIN$eigen

    x <- matrix(unlist(lapply(1:ncol(x), function(i) {
      lapply(1:ncol(dw), function(j) {
        x[, i] * dw[, j]
      })
    })), nrow = nrow(x), byrow = FALSE)
  }

  if (return.working == TRUE) {
    working <- PaleoSpec:::mean.ACF(x, return.working = TRUE)
    acf <- working[["acf"]]
  } else {
    acf <- PaleoSpec:::mean.ACF(x)
  }

  freq <- (lag / (nrow(x))) / bin.width

  if (k > 1){
    spec <- Re(stats::fft(acf)) * bin.width * nrow(x)
  } else {
    spec <- Re(stats::fft(acf)) * bin.width
  }




  rfreq <- max(freq)

  if (pos.f.only) {
    spec <- spec[freq > 0 & freq <= 1 / (2 * bin.width)]
    freq <- freq[freq > 0 & freq <= 1 / (2 * bin.width)]
  }

  out <- list(
    bin.width = bin.width,
    rfreq = rfreq,
    nrec = ncolx,
    lag = lag, acf = acf,
    freq = freq, spec = spec, f.length = 1,
    dof = rep(2 * ncolx * k, length(freq))
  )

  class(out) <- c("SpecACF", "spec")

  if (return.working) {
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
#' @param m a numeric matrix, data.frame or vector
#' @param trim trim leading and trailing rows of "all" NA or containing "any" NA values
#' @return a numeric matrix, data.frame or vector
#' @export
#' @examples
#' m <- matrix(c(NA, NA, NA, 1, NA, NA, NA, 1, 1, NA, NA, NA, 1:9, NA,NA,NA, 10:12, NA, 1, NA, NA,NA,NA), ncol = 3, byrow = TRUE)
#' m
#' TrimNA(m)
#' TrimNA(m, trim = "any")
TrimNA <- function(m, trim = c("all", "any")) {
  trim <- match.arg(trim)

  # make it work on vectors
  class_m <- class(m)
  if (class_m[1] == "numeric") {
    m <- cbind(m)
  }

  if (trim == "all") {
    empty.row <- is.nan(rowMeans(m, na.rm = TRUE))
    rank.good <- (empty.row == FALSE) * 1:length(empty.row)
  } else if (trim == "any") {
    empty.row <- is.na(rowMeans(m))
    rank.good <- (empty.row == FALSE) * 1:length(empty.row)
  }

  first.good <- which.min(empty.row * 1:length(empty.row))
  last.good <- which.max(rank.good)

  m <- m[first.good:last.good, , drop = FALSE]

  # return to a vector
  if (class_m[1] == "numeric") {
    m <- as.numeric(m)
  }

  return(m)
}


