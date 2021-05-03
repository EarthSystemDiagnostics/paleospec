#' Filter time series
#'
#' Apply a given filter to a time series using different endpoint constraints as
#' described below.
#'
#' Note that when passing objects of class \code{ts}, the time step provided is
#' not used; thus, for time series with a time step different from 1, the filter
#' has to be adapted accordingly.
#'
#' The function applies endpoint constrains following Mann et al., GRL, 2003;
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
#' @return filtered timeseries as a ts object.
#' @author Thomas Laepple
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
#' @export
ApplyFilter <- function(data, filter, method = 0) {

  if (!method %in% (0 : 4))
    stop("Unknown method; only 0 : 4 available.")

  circular = FALSE

  if (method == 0 | method == 4) {

    if (method == 4) {circular = TRUE}

    result <- stats::filter(c(data), filter, circular = circular)

  } else {

    N <- floor(length(filter) / 2)

    if (method == 1) {

      before <- rep(mean(data), N)
      after  <- rep(mean(data), N)

    } else if (method == 2 | method == 3) {

      before <- c(data)[N : 1]
      after  <- c(data)[length(data) : (length(data) - N + 1)]

      if (method == 3) {

        before <- c(data)[1] - (before - mean(before))
        after  <- c(data)[length(data)] - (after - mean(after))

      }
    }

    result <- stats::filter(c(before, data, after), filter, circular = circular)
    result <- result[(N + 1) : (N + length(data))]

  }

  return(ts(result, frequency = frequency(data)))

}
