#' Interpolate spectrum
#'
#' Interpolate a spectrum to a given frequency resolution. The interpolation
#' includes the spectral estimates themselves as well as their degrees of
#' freedom.
#'
#' @param freqRef numeric vector with the target frequency resolution.
#' @param spec a spectral object of class \code{"spec"} or a list with the
#'   minimum components \code{freq}, \code{spec} and \code{dof} which are
#'   vectors of the same length giving the original frequency resolution and the
#'   corresponding spectral estimates and degrees of freedom.
#' @return a spectral object of class \code{"spec"} with the spectrum
#'   interpolated to the target resolution.
#' @author Thomas Laepple
#' @examples
#'
#' freqRef <- seq(0.1, 0.5, 0.1)
#' spec <- list(freq = c(0.1, 0.2, 0.4, 0.5), spec = c(1, 2, 4, 5),
#'              dof = rep(1, 4))
#'
#' SpecInterpolate(freqRef, spec)
#' @export
SpecInterpolate <- function(freqRef, spec) {

  result <- list(
    freq = freqRef,
    spec = approx(spec$freq, spec$spec, freqRef)$y,
    dof  = approx(spec$freq, spec$dof, freqRef)$y
  )

  class(result) <- "spec"

  return(result)

}
