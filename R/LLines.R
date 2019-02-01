#' Add a spectrum to an existing log-log spectral plot.
#'
#' @description This function adds a spectrum to an existing double-logarithmic plot and
#' optionally adds a transparent confidence interval.
#' @param x  a spectral object resulting from a call to \code{\link{SpecMTM}}.
#' @param conf if \code{TRUE} (the default) add a transparent confidence
#' interval (suppressed if \code{x} contains no error limits).
#' @param bPeriod if \code{TRUE} treat the x-axis values in units of period
#'     (inverse frequency). Defaults to \code{FALSE}.
#' @param col color for the line plot and the confidence interval.
#' @param alpha transparency level (between 0 and 1) for the confidence
#' interval. Defaults to \code{0.3}.
#' @param removeFirst omit \code{removeFirst} values on the low frequency side.
#' @param removeLast omit \code{removeLast} values on the high frequency side.
#' @param ... further graphical parameters passed to \code{lines}.
#' @examples
#' x <- ts(arima.sim(list(ar = 0.9), 1000))
#' spec <- SpecMTM(x)
#' LPlot(spec, col = "grey")
#' LLines(LogSmooth(spec), lwd = 2)
#' @author Thomas Laepple
#' @export
LLines<-function(x, conf = TRUE, bPeriod = FALSE, col = "black", alpha = 0.3,
                 removeFirst = 0, removeLast = 0, ...) {

    if (bPeriod) x$freq <- 1/x$freq

    index <- (removeFirst + 1) : (length(x$freq) - removeLast)

    x$freq  <- x$freq[index]
    x$spec  <- x$spec[index]
    x$lim.1 <- x$lim.1[index]
    x$lim.2 <- x$lim.2[index]

    lim <- !is.null(x$lim.1) & !is.null(x$lim.2)
    if (conf & lim) {
        polygon(c(x$freq, rev(x$freq)), c(x$lim.1, rev(x$lim.2)),
                col = ColTransparent(col, alpha), border = NA)
    }

    lines(x$freq, x$spec, col = col, ...)

}
