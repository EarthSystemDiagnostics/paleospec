#' @title Derives and plots the transfer function (given a filter)
#' @description Get the transfer function of a symetric filter, page 122 in Bloomfield 1976,
#' page 135 in Bloomfield 2000
#' @param g.u a filter, a numeric vector, values should sum to 1
#' @param resolution the number of frequencies at which to evaluate the
#'  transfer function
#' @param freq the specific frequency(s) at which to evaluate the transfer function
#' if NULL, transfer function is evaluated at 1:resolution frequencies
#' @param bPlot logical, plot the transfer function, defaults to FALSE if frequencies
#' are specified, TRUE otherwise
#' @param add logical, add to a previous plot
#' @param ... other arguments to pass to the plotting function
#' @return list(freq, y) containing the transfer function
#' @author Thomas Laepple
#' @examples
#' l <- 11
#' tf <- GetTransferFunction(rep(1/l, l), resolution = 1000)
#' tf <- GetTransferFunction(rep(1/l, l), freq = 1/c(100, 10, 2))
#' @export
GetTransferFunction<-function(g.u, resolution=100, freq = NULL,
                              bPlot=is.null(freq)==TRUE, add=FALSE, ...)
{
  n <- length(g.u)
  n.side <- (n - 1) / 2

  if (is.null(freq)){
    omega = (1:resolution) * pi / resolution}else{
      if (max(freq)>=1) warning("Maximum frequency must be < 1.
                                Returning frequencies < 1")
      freq[freq>=1] <- NA
      omega = 2 * pi * freq
      }

  yt <- rep(0, length(omega))

  for (u in (-1 * n.side):n.side)
    yt <- yt + (g.u[u + n.side + 1] * exp(-1 * 1i * omega * u))

  if (bPlot) {
    if (!add) {
      plot(omega / 2 / pi, abs(yt)^2, type = "l",
           xlab = "frequency", main = "Transfer function", ...)
    }
    else
      lines(omega / 2 / pi, abs(yt)^2, col = "blue", ...)
  }
  return(list(freq = omega / 2 / pi, y = abs(yt)^2))
}
