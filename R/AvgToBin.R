#' Bin averaging
#'
#' Average a vector into bins.
#'
#' This function averages the vector \code{y} into bins according to the positon
#' of \code{x} within the breaks. You can either specify a desired number N of
#' breaks which are used to calculate the actual breaks via \code{pretty(x, N)},
#' or directly specify the N + 1 break positions. For \code{right = TRUE} (the
#' default) the averaging bins are defined via \code{x > breaks[i]} and \code{x
#' <= breaks[i + 1]}, else they are defined via \code{x >= breaks[i]} and
#' \code{x < breaks[i + 1]}. If \code{bFill = TRUE}, empty bins are filled using
#' linear interpolation from the neighbours to the center of the bin.
#'
#' Probably the binning could be considerably speeded up by using \code{?cut}.
#'
#' @param x vector of values on which the data in \code{y} is tabulated;
#'   e.g. depth or time points.
#' @param y vector of observation values to be averaged into bins. Must have the
#'   same length as \code{x}.
#' @param N desired number of breaks (ignored if \code{breaks} are supplied
#'   directly).
#' @param breaks vector of break point positions to define the averagig bins; if
#'   omitted, break point positions are calculated from the range of \code{x}
#'   and the desired number of breaks given by \code{N}.
#' @param right logical; indicate whether the bin intervals should be closed on
#'   the right and open on the left (\code{TRUE}, the default), or vice versa
#'   (\code{FALSE}).
#' @param bFill logical; if \code{TRUE}, fill empty bins using linear
#'   interpolation from the neighbours to the center of the bin.
#'
#' @return a list with four elements:
#' \describe{
#' \item{\code{breaks}:}{numeric vector of the used break point positions.}
#' \item{\code{centers}:}{numeric vector with the positions of the bin centers.}
#' \item{\code{avg}:}{numeric vector with the bin-averaged values.}
#' \item{\code{nobs}:}{numeric vector with the number of observations
#'   contributing to each bin average.}
#' }
#'
#' @author Thomas Laepple
#' @export
#' @examples
#' N <- 100
#' x1 <- seq(1, N, by = 1)
#' y1 <- SimPLS(N, a = 0.1, b = 1)
#' plot(x1, y1, type = "l")
#'
#' y2 <- AvgToBin(x1, y1, 13)
#' lines(y2$centers, y2$avg, col = "green")
#'
#' #Add some NA values to the timeseries
#' y1[(N/2):(N/2 +5)] <- NA
#' plot(x1, y1, type = "l")
#'
#' y2 <- AvgToBin(x1, y1, 25)
#' lines(y2$centers, y2$avg, col = "green")
#'
#' # Large enough bins will average across the gap
#' y3 <- AvgToBin(x1, y1, 10)
#' lines(y3$centers, y3$avg, col = "red")


#' # Or interpolate to nearest neighbour
#' y3 <- AvgToBin(x1, y1, 25, bFill = TRUE)
#' lines(y3$centers, y3$avg, col = "blue")

AvgToBin <- function(x, y, N = 2, breaks = pretty(x, N),
                     right = TRUE, bFill = FALSE) {

  if (length(x) != length(y)) {
    stop("'x' and 'y' must have the same length.", call. = FALSE)
  }

  nBins <- length(breaks) - 1

  centers <- (breaks[1 : nBins] + breaks[2 : (nBins + 1)]) / 2
  nObs <- avg <- rep(NA, nBins)

  for (i in 1 : nBins) {

    if (right) {
      selection <- y[which((x > breaks[i]) & (x <= breaks[i + 1]))]
    } else {
      selection <- y[which((x >= breaks[i]) & (x < breaks[i + 1]))]
    }

    avg[i]  <- mean(na.omit(selection))
    nObs[i] <- sum(!is.na(selection))

  }

  if ((sum(missing <- is.na(avg)) > 0) & (bFill)) {

    avg[missing] <- (approx(x, y, centers)$y)[missing]

  }

  list(breaks = breaks, centers = centers, avg = avg, nobs = nObs)

}
