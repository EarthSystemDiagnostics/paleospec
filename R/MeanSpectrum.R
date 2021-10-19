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
#'   elements in \code{specList}.
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

    specList <- lapply(specList, SpecInterpolate, freqRef = freqRef)

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


# Helper functions

#' Remove low frequencies
#'
#' Remove the lowest frequencies from a spectra, since they are biased by
#' detrending and the MTM estimator.
#'
#' @param spec a list object of class \code{"spec"}.
#' @param iRemove integer; number of lowest frequencies to remove.
#' @return the input data with the data at the \code{iRemove} lowest
#'   frequencies removed.
#' @author Thomas Laepple
remove.lowestFreq <- function(spec, iRemove) {

  if (iRemove == 0) {

    return(spec)

  } else {

    index <- 1 : iRemove
    spec$spec <- spec$spec[-index]
    spec$freq <- spec$freq[-index]
    spec$dof  <- spec$dof[-index]

    if (has.limits(spec)) {

      spec$lim.1 <- spec$lim.1[-index]
      spec$lim.2 <- spec$lim.2[-index]
    }

    return(spec)

  }
}

#' Remove high frequencies
#'
#' Remove the highest frequencies from a spectra, since they may be biased by
#' the MTM estimator.
#'
#' @param spec a list object of class \code{"spec"}.
#' @param iRemove integer; number of highest frequencies to remove.
#' @return the input data with the data at the \code{iRemove} highest
#'   frequencies removed.
#' @author Thomas Münch
remove.highestFreq <- function(spec, iRemove) {

  if (iRemove == 0) {

    return(spec)

  } else {

    index <- 1 : (length(spec$freq) - iRemove)
    spec$spec <- spec$spec[index]
    spec$freq <- spec$freq[index]
    spec$dof  <- spec$dof[index]

    if (has.limits(spec)) {

      spec$lim.1 <- spec$lim.1[index]
      spec$lim.2 <- spec$lim.2[index]
    }

    return(spec)

  }
}

# test whether spectral object has confidence interval included
has.limits <- function(spec) !is.null(spec$lim.1) & !is.null(spec$lim.2)

# Helper functions to access the list elements to determine
# the target frequency discretisation
get.df <- function(x) return(mean(diff(x$freq)))
get.fend.existing <- function(x) return(max(x$freq[!is.na(x$spec)]))
get.fstart.existing <- function(x) return(min(x$freq[!is.na(x$spec)]))

# helper functions to obtain certain elements from spectral object
get.length <- function(x) return(length(x$freq))
get.freq <- function(x) return(x$freq)
get.spec <- function(x) return(x$spec)
get.dofs <- function(x) return(x$dof)

#' Check for spectral object
#'
#' In \code{PaleoSpec} a "spectral object" is defined to be a list consisting as
#' the minimum requirement of the three elements \code{freq}, \code{spec}, and
#' \code{dof}, which are numeric vectors of the same length giving the frequency
#' axis, the power spectral density, and the degrees of freedom of the spectral
#' estimate. Most functions in this package require such an object as their main
#' input. This function checks an object for the minimum requirements.
#'
#' @param x an object to be checked.
#' @param check.only logical; the default means to issue an error when a
#'   requirement is not met. If set to \code{TRUE}, all requirements are
#'   checked, and the function returns \code{TRUE} if they all pass, or
#'   \code{FALSE} if one or more checks fail.
#' @param dof logical; per default, \code{x} is also checked for containing a
#'   \code{dof} vector, which can be turned off here when unnecessary.
#' @return a logical value for \code{check.only = TRUE}.
#' @author Thomas Münch
is.spectrum <- function(x, check.only = FALSE, dof = TRUE) {

  isList <- is.list(x)

  if (isList) {

    hasFreq <- !is.null(get.freq(x))
    hasSpec <- !is.null(get.spec(x))
    hasDofs <- ifelse(dof, !is.null(get.dofs(x)), TRUE)

    hasAll <- all(hasFreq, hasSpec, hasDofs)

    if (hasAll) {
      lengths <- c(get.length(x), length(get.spec(x)))
      if (dof) lengths <- c(lengths, length(get.dofs(x)))
      hasEqualLength <- stats::var(lengths) == 0

    } else {
      hasEqualLength <- FALSE
    }
  }

  if (check.only) {

    if (isList) {
      hasAll & hasEqualLength
    } else {
      isList
    }

  } else {

    if (!isList) {
      stop("Passed argument is not a spectral list object.", call. = FALSE)
    }
    if (!hasFreq) {
      stop("Passed object has no frequency vector.", call. = FALSE)
    }
    if (!hasSpec) {
      stop("Passed object has no spectral density vector.", call. = FALSE)
    }
    if (!hasDofs) {
      stop("Passed object has no dof vector.", call. = FALSE)
    }
    if (!hasEqualLength) {
      stop("Frequency, PSD and DOF vectors have different lengths.",
           call. = FALSE)
    }
  }

}

# f (bNormalize) Normalize to the mean of the common interval
# fmin<-min(unlist(lapply(temp,get.fend.existing)))
# fmax<-max(unlist(lapply(temp,get.fstart.existing)))
# i.min<-ClosestElement(freqRef,fmin)
# i.max<-ClosestElement(freqRef,fmax) var.band<-vector() for
# (i in 1:length(temp)) {
# var.band[i]<-mean(temp[[i]]$spec[i.min:i.max]) }
# rescale<-mean(var.band)/var.band if (!is.null(weights))
# rescale=weights for (i in 1:length(temp))
# temp[[i]]$spec<-temp[[i]]$spec*rescale[i]
