##' Log-log spectral plot.
##'
##' This function plots a spectrum on a double-logarithmic scale and optionally
##' adds a transparent confidence interval.
##' @param x a spectral object resulting from a call to \code{\link{SpecMTM}}.
##' @param conf if \code{TRUE} (the default) add a transparent confidence
##' interval (suppressed if \code{x} contains no error limits).
##' @param bPeriod if \code{TRUE} the x-axis is displayed in units of period
##'     (inverse frequency), increasing to the left. Defaults to \code{FALSE}.
##' @param col color for the line plot and the confidence interval.
##' @param alpha transparency level (between 0 and 1) for the confidence
##' interval. Defaults to \code{0.3}.
##' @param removeFirst omit \code{removeFirst} values on the low frequency side.
##' @param removeLast omit \code{removeLast} values on the high frequency side.
##' @param xlab character string for labelling the x-axis.
##' @param ylab character string for labelling the y-axis.
##' @param xlim range of x-axis values; if \code{NULL} (the default) it is
##'     calculated internally and automatically reversed for
##'     \code{bPeriod = TRUE}.
##' @param ylim range of y-axis values; if \code{NULL} (the default) it is
##'     calculated internally.
##' @param ... further graphical parameters passed to \code{plot}.
##' @examples
##' x <- ts(arima.sim(list(ar = 0.9), 1000))
##' spec <- SpecMTM(x)
##' LPlot(spec, col = "grey")
##' LLines(LogSmooth(spec), lwd = 2)
##' @author Thomas Laepple
##' @export
LPlot <- function(x, conf = TRUE, bPeriod = FALSE, col = "black", alpha = 0.3,
                  removeFirst = 0, removeLast = 0, xlab = "f", ylab = "PSD",
                  xlim = NULL, ylim = NULL, ...) {
    
    if (bPeriod) {
        x$freq <- 1/x$freq
        if (is.null(xlim)) xlim <- rev(range(x$freq))
        if (xlab == "f") xlab <- "period"
    }

    index <- (removeFirst + 1) : (length(x$freq) - removeLast)
    
    x$freq  <- x$freq[index]
    x$spec  <- x$spec[index]
    x$lim.1 <- x$lim.1[index]
    x$lim.2 <- x$lim.2[index]

    plot(x$freq, x$spec, type = "n", log = "xy", xlab = xlab, ylab = ylab,
         xlim = xlim, ylim = ylim, ...)

    lim <- !is.null(x$lim.1) & !is.null(x$lim.2)
    if (conf & lim) {
        polygon(c(x$freq, rev(x$freq)), c(x$lim.1, rev(x$lim.2)),
                col = ColTransparent(col, alpha), border = NA)
    }

    lines(x$freq, x$spec, col = col, ...)
    
}
