#' Weighted mean spectrum
#'
#' Calculate the mean spectrum from a list of individual spectra including
#' weighting of the individual spectra and, if needed, interpolation to
#' the highest resolution frequency grid. The default weighting produces the
#' simple arithmetic mean. Interpolation is only performed when the frequency
#' axes of the individual spectra have different lengths and/or differ in
#' frequency discretization.
#'
#' @param specList list of spectra, i.e. objects of class \code{"spec"} (see
#'   \code{\link{SpecMTM}} for details), where each spectrum has to be a list of
#'   the vectors \code{freq}, \code{spec} and \code{dof} of a common length.
#' @param iRemoveLowest integer; number of lowest frequencies to remove from
#'   each individual spectral estimate (e.g. to remove detrending bias) prior to
#'   the interpolation and averaging.
#' @param weights numeric vector of weights; its length must match the number of
#'   elements in \code{specList}. The default setting of identical weights
#'   yields the simple arithmetic average of the input spectra; change this
#'   parameter accordingly to apply non-arithmetic averaging (see example).
#' @return object of class \code{"spec"} with the weighted mean spectrum,
#'   amended by the element \code{nRecord} which gives the number of records
#'   contributing to each mean spectral estimate.
#' @author Thomas Laepple and Thomas Münch
#' @examples
#'
#' # Simple arithmetic average
#' f1 <- 1 : 5
#' f2 <- f1
#' s1 <- rep(1, length(f1))
#' s2 <- rep(3, length(f2))
#' dof1 <- rep(1, length(f1))
#' dof2 <- rep(1, length(f2))
#'
#' spectra <- list(list(freq = f1, spec = s1, dof = dof1),
#'                 list(freq = f2, spec = s2, dof = dof2))
#'
#' MeanSpectrum(spectra, iRemoveLowest = 0)
#'
#' # Weighted mean with interpolation
#' f1 <- 1 : 5
#' f2 <- 3 : 8
#' s1 <- rep(1, length(f1))
#' s2 <- rep(3, length(f2))
#' dof1 <- rep(1, length(f1))
#' dof2 <- rep(1, length(f2))
#'
#' spectra <- list(list(freq = f1, spec = s1, dof = dof1),
#'                 list(freq = f2, spec = s2, dof = dof2))
#'
#' MeanSpectrum(spectra, iRemoveLowest = 0, weights = c(1, 2))
#' # with some detrending bias removal
#' MeanSpectrum(spectra, iRemoveLowest = 1, weights = c(1, 2))
#' @export
#'
MeanSpectrum <- function(specList, iRemoveLowest = 1,
                         weights = rep(1, length(specList))) {

  if (!is.list(specList)) {
    stop("'specList' input needs to be a list.")
  }

  ns <- length(specList)

  if (!all(sapply(specList, is.spectrum, check.only = TRUE))) {
    stop("Non-spectral objects in input list.")
  }

  if (length(weights) != ns) {
    stop("Need as many weights as input spectra.")
  }

  # remove lowest/biased frequencies from spectral estimates
  specList <- lapply(specList, remove.lowestFreq, iRemove = iRemoveLowest)

  # interpolate if spectra have different lengths
  interpolate <- stats::var(sapply(specList, get.length)) > 0

  if (!interpolate) {
    # equal lengths but check if frequency axes are different
    freqs <- sapply(specList, get.freq)
    interpolate <- !all(apply(freqs, 1, stats::var) == 0.)
    # later, one could introduce some tolerance here for the deviation between
    # the individual frequency axes...
  }

  if (interpolate) {

    # use the longest run for the reference spectrum
    freqRef <- seq(from = min(sapply(specList, get.fstart.existing)),
                   to = max(sapply(specList, get.fend.existing)),
                   by = min(sapply(specList, get.df)))

    specList <- lapply(specList, SpecInterpolate,
                       freqRef = freqRef, check = FALSE)

  } else {

    freqRef <- get.freq(specList[[1]])
  }

  # Build matrices from input data

  nf <- length(freqRef)

  specMatrix <- sapply(specList, get.spec)
  dofMatrix  <- sapply(specList, get.dofs)

  # weights at each frequency
  weightMatrix <- matrix(weights, nrow = nf, ncol = ns, byrow = TRUE)
  # missing spectral estimates
  missingObs <- is.na(specMatrix)
  # contributing number of records
  nRecord <- rowSums(!missingObs)

  # put weights to NA at missing estimates
  weightMatrix[missingObs] <- NA

  # normalise weights across spectra
  weightMatrix <- t(apply(weightMatrix, 1, function(x) {
    x / sum(x, na.rm = TRUE)}))

  # mean weighted spectra
  result <- list()
  result$freq <- freqRef
  result$spec <- rowSums(specMatrix * weightMatrix, na.rm = TRUE)
  result$dof <- rowSums(dofMatrix * weightMatrix, na.rm = TRUE) * nRecord
  result$nRecord <- nRecord
  class(result) <- "spec"

  return(AddConfInterval(result))
}
