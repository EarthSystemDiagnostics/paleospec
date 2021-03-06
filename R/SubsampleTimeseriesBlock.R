#' @title Subsample (downsample) timeseries using block averaging#'
#' @param ts  ts object or vector containing the equidistant timeseries
#' @param timepoints  vector with the points in time
#' @description  Resample a equidistant timeseries (e.g. model result) at the
#'   'timepoints' using block averaging. The blocks are divided at 1/2 time
#'   between the requested output points. For the first (and last) timepoint,
#'   the interval starting mean(diff(timepoints)) before (ending after) are
#'   used. Example usage is to downsample a model timeseries to mimick an
#'   integrating proxy (e.g. water isotopes that are measured by melting pieces
#'   of ice).
#' @return values at timepoints
#' @examples
#' input <- ts(SimPowerlaw(0.5, 1000))
#' timepoints <- seq(from = 50, to = 950, by = 50)
#' result <- SubsampleTimeseriesBlock(input, timepoints)
#' plot(input, main = "Comparison of block avg. vs. simple interpolation",
#'   ylab = "unitless")
#' points(timepoints, result, pch = 19, col = "red", lwd = 3)
#' points(approx(time(input), c(input), timepoints), col = "green",
#'   pch = 10, lwd = 3)
#' legend("bottom", col = c("black", "red", "green"), lwd = 2, bty = "n",
#'   c("High-resolution timeseries (input)", " Block Avg", "interpolated values"))
#' @author Thomas Laepple
#' @export
SubsampleTimeseriesBlock <- function(ts, timepoints) {
  result <- list()
  dt <- mean(diff(timepoints))
  timepoints.bound <- c(head(timepoints, 1) - dt/2, timepoints[-length(timepoints)] +
    diff(timepoints)/2, tail(timepoints + dt/2, 1))
  breaks <- .bincode(c(time(ts)), breaks = timepoints.bound,
    TRUE, TRUE)
  temp <- tapply(c(ts), breaks, mean, na.rm = TRUE)
  return(temp)
}


#' @title Bin daily values to monthly values
#' @description Assumes months of equal length and a 365 day long year
#' @param ts.daily vector of 365 values
#' @return vector of 12 values
#' @author Thomas Laepple
#' @export
MonthlyFromDaily <- function(ts.daily) {
  mLength <- (365/12)
  monthBreaks <- mLength/2 + (0:11) * mLength
  monthly <- SubsampleTimeseriesBlock(ts.daily, timepoints = monthBreaks)
  return(monthly)
}

