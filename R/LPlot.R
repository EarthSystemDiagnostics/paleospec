#' Log-log spectral plot.
#'
#' @description This function plots a spectrum on a double-logarithmic scale and
#'  optionally adds a transparent confidence interval.
#' @param x a spectral object resulting from a call to \code{\link{SpecMTM}}.
#' @param conf if \code{TRUE} (the default) add a transparent confidence
#' interval (suppressed if \code{x} contains no error limits).
#' @param bPeriod if \code{TRUE} the x-axis is displayed in units of period
#'     (inverse frequency), increasing to the left. Defaults to \code{FALSE}.
#' @param bNoPlot if \code{TRUE} only produce the plot frame (\code{type = "n"}
#' behaviour of function \code{\link{plot}}). Defaults to \code{FALSE}.
#' @param axes if \code{FALSE} the plotting of the x and y axes is
#' suppressed. Defaults to \code{TRUE}.
#' @param col color for the line plot and the confidence interval.
#' @param alpha transparency level (between 0 and 1) for the confidence
#' interval. Defaults to \code{0.3}.
#' @param removeFirst omit \code{removeFirst} values on the low frequency side.
#' @param removeLast omit \code{removeLast} values on the high frequency side.
#' @param xlab character string for labelling the x-axis.
#' @param ylab character string for labelling the y-axis.
#' @param xlim range of x-axis values; if \code{NULL} (the default) it is
#'     calculated internally and automatically reversed for
#'     \code{bPeriod = TRUE}.
#' @param ylim range of y-axis values; if \code{NULL} (the default) it is
#'     calculated internally.
#' @inheritParams graphics::plot.default
#' @param ... further graphical parameters passed to \code{plot}.
#' @family functions to plot power spectra
#' @examples
#' x <- ts(arima.sim(list(ar = 0.9), 1000))
#' spec <- SpecMTM(x)
#' LPlot(spec, col = "grey")
#' LLines(LogSmooth(spec), lwd = 2)
#' @author Thomas Laepple
#' @export
LPlot <- function(x, conf = TRUE, bPeriod = FALSE, bNoPlot = FALSE, axes = TRUE,
                  col = "black", alpha = 0.3, removeFirst = 0, removeLast = 0,
                  xlab = "f", ylab = "PSD", xlim = NULL, ylim = NULL, log = "xy", ...) {

    is.spectrum(x, dof = FALSE)

    if (bPeriod) {
        x$freq <- 1/x$freq
        if (is.null(xlim)) xlim <- rev(range(x$freq))
        if (xlab == "f") xlab <- "period"
    }

    x <- remove.lowestFreq(x, iRemove = removeFirst)
    x <- remove.highestFreq(x, iRemove = removeLast)

    plot(x$freq, x$spec, type = "n", log = log, xlab = xlab, ylab = ylab,
         xlim = xlim, ylim = ylim, axes = axes, ...)

    if (conf & has.limits(x) & !bNoPlot) {
        polygon(c(x$freq, rev(x$freq)), c(x$lim.1, rev(x$lim.2)),
                col = ColTransparent(col, alpha), border = NA)
    }

    if (!bNoPlot) lines(x$freq, x$spec, col = col, ...)

}


#' Add a spectrum to an existing log-log spectral plot.
#'
#' @description This function adds a spectrum to an existing double-logarithmic plot and
#' optionally adds a transparent confidence interval.
#' @param ... further graphical parameters passed to \code{lines}.
#' @inheritParams LPlot
#' @family functions to plot power spectra
#' @examples
#' x <- ts(arima.sim(list(ar = 0.9), 1000))
#' spec <- SpecMTM(x)
#' LPlot(spec, col = "grey")
#' LLines(LogSmooth(spec), lwd = 2)
#' @author Thomas Laepple
#' @export
LLines<-function(x, conf = TRUE, bPeriod = FALSE, col = "black", alpha = 0.3,
                 removeFirst = 0, removeLast = 0, ...) {

  is.spectrum(x, dof = FALSE)

  if (bPeriod) x$freq <- 1/x$freq

  x <- remove.lowestFreq(x, iRemove = removeFirst)
  x <- remove.highestFreq(x, iRemove = removeLast)

  if (conf & has.limits(x)) {
    polygon(c(x$freq, rev(x$freq)), c(x$lim.1, rev(x$lim.2)),
            col = ColTransparent(col, alpha), border = NA)
  }

  lines(x$freq, x$spec, col = col, ...)

}

