#' Filter time series
#'
#' Apply a given filter to a time series using different endpoint constraints.
#'
#' Note that when passing objects of class \code{ts}, the time step provided is
#' not used; thus, for time series with a time step different from 1, the filter
#' has to be adapted accordingly.
#'
#' Leading and trailing NA values are automatically stripped from the input
#' vector so that they do not spread into the filtered data when applying the
#' endpoint constraints, but added in again after filtering so that the output
#' vector has the same length as the input. This does not apply to any internal
#' NA values, which instead are handled by \code{na.rm}.
#'
#' The function applies endpoint constrains following Mann et al., GRL, 2004;
#' available methods are:
#' \itemize{
#'   \item method = 0: no constraint (loss at both ends);
#'   \item method = 1: minimum norm constraint;
#'   \item method = 2: minimum slope constraint;
#'   \item method = 3: minimum roughness constraint;
#'   \item method = 4: circular filtering.
#' }
#'
#' @param data numeric vector with the input timeseries (standard or ts object).
#' @param filter numeric vector of filter weights.
#' @param method single integer for choosing an endpoint constraint method;
#'   available choices are integers 0-4, see details.
#' @param na.rm logical; control the handling of internal NA values in
#'   \code{data}. If set to \code{TRUE}, any internal NA values are removed by
#'   linear interpolation from the neighbouring values; defaults to
#'   \code{FALSE}.
#' @return a ts object with the filtered timeseries.
#' @author Thomas Laepple
#' @source The endpoint constraint methods are based on the study:\cr
#'   Michael E. Mann, On smoothing potentially non‐stationary climate time
#'   series, Geophys. Res. Lett., 31, L07214, doi:10.1029/2004GL019569, 2004.
#' @examples
#' # Simple running mean filter across three bins
#'
#' x <- 1 : 10
#' filter <- rep(1 / 3, 3)
#'
#' # no endpoint constraints lead to loss at both ends
#' ApplyFilter(x, filter, method = 0)
#'
#' # circular filtering avoids end losses, so as the other methods
#' ApplyFilter(x, filter, method = 4)
#'
#' # leading and trailing NA's are ignored but added in again afterwards
#' x <- c(NA, 1 : 10, NA)
#' ApplyFilter(x, filter, method = 4)
#'
#' # ... but not internal NA's
#' x <- c(1 : 5, NA, 7 : 10)
#' ApplyFilter(x, filter, method = 4)
#'
#' # if not explicitly removed by linear interpolation
#' ApplyFilter(x, filter, method = 4, na.rm = TRUE)
#'
#'
#' # Visual comparison of methods ----------
#' set.seed(20220302)
#'
#' x <- PaleoSpec::SimPowerlaw(1, 1e02)
#' x <- x + 0.1 * (1:length(x))
#'
#' filt <- rep(1/30, 30)
#'
#' plot(x, type = "l")
#'
#' x0 <- ApplyFilter(x, filt, method = 0)
#' lines(x0, col = "Green")
#'
#' x1 <- ApplyFilter(x, filt, method = 1)
#' lines(x1, col = "blue")
#'
#' x2 <- ApplyFilter(x, filt, method = 2)
#' lines(x2, col = "red")
#'
#' x3 <- ApplyFilter(x, filt, method = 3)
#' lines(x3, col = "orange")
#'
#' x4 <- ApplyFilter(x, filt, method = 4)
#' lines(x4, col = "Purple", lty = 2)
#'
#' lines(x0, col = "Green")
#'
#' legend(x = "topleft",
#'        legend = c("0: no constraint (lose ends)",
#'                   "1: min norm (pad with mean)",
#'                   "2: min slope (reflect ends)",
#'                   "3: min roughness (reflect and invert ends)",
#'                   "4: circular"),
#'        col = c("Green", "Blue", "Red", "Orange", "Purple"),
#'        lwd = 2, lty = c(1,1,1,1,2))
#'
#'
#' # Repeat with linear trend, no noise
#'
#' x <- 1:100
#' filt <- rep(1/30, 30)
#'
#' plot(x, type = "l")
#'
#' x0 <- ApplyFilter(x, filt, method = 0)
#' lines(x0, col = "Green")
#'
#' x1 <- ApplyFilter(x, filt, method = 1)
#' lines(x1, col = "blue")
#'
#' x2 <- ApplyFilter(x, filt, method = 2)
#' lines(x2, col = "red")
#'
#' x3 <- ApplyFilter(x, filt, method = 3)
#' lines(x3, col = "orange")
#'
#' x4 <- ApplyFilter(x, filt, method = 4)
#' lines(x4, col = "Purple", lty = 2)
#'
#' lines(x0, col = "Green")
#'
#' legend(x = "topleft",
#'        legend = c("0: no constraint (lose ends)",
#'                   "1: min norm (pad with mean)",
#'                   "2: min slope (reflect ends)",
#'                   "3: min roughness (reflect and invert ends)",
#'                   "4: circular"),
#'        col = c("Green", "Blue", "Red", "Orange", "Purple"),
#'        lwd = 2, lty = c(1,1,1,1,2))
#' @export
ApplyFilter <- function(data, filter, method = 0, na.rm = FALSE) {

  if (!method %in% (0 : 4))
    stop("Unknown method; only 0 : 4 available.")

  result <- rep(NA, length(data))

  # remove leading and trailing NA's
  x <- c(zoo::na.trim(data))
  n <- length(x)

  # linearly interpolate internal NA's if requested
  if (na.rm) {x <- stats::approx(1 : n, x, 1 : n)$y}

  circular = FALSE

  if (method == 0 | method == 4) {

    if (method == 4) {circular = TRUE}

    xf <- stats::filter(x, filter, circular = circular)

  } else {

    N <- floor(length(filter) / 2)

    if (method == 1) {

      before <- rep(mean(x), N)
      after  <- rep(mean(x), N)

    } else if (method == 2 | method == 3) {

      before <- x[N : 1]
      after  <- x[n : (n - N + 1)]

      if (method == 3) {

        before <- x[1] - (before - mean(before))
        after  <- x[n] - (after - mean(after))

      }
    }

    xf <- stats::filter(c(before, x, after), filter, circular = circular)
    xf <- xf[(N + 1) : (N + n)]

  }

  i <- seq(match(x[1], data), by = 1, length.out = n)
  result[i] <- xf

  return(ts(result, frequency = frequency(data)))

}
