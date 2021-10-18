#' @title average spectra with weighting
#' @description Calculate the weighted mean spectrum of all spectra by interpolating them to
#' the highest resolution frequency grid and averaging them.
#'
#' Spectra can have different resolution and span a different freq range.
#' @param speclist list of spectra
#' @param iRemoveLowest number of lowest frequencies to remove (e.g. to remove detrending bias)
#' @param weights vector of weights (same length as elements in speclist)
#' @return list(spec,nRecords) spec=average spectrum, nRecords = number of records contributing to each spectral estimate
#' @author Thomas Laepple
#' @export
#'
MeanSpectrum <- function(specList, iRemoveLowest = 1, weights = rep(1,
  length(specList))) {

  if (length(weights) != length(specList))
    stop("specList and weights have a different number of elements")

  # Remove the lowest/biased frequencies from the spectral
  # estimates
  specList <- lapply(specList, remove.lowestFreq, iRemove = iRemoveLowest)

  # Check if interpolation is needed

  # interpolate if spectra have different lengths
  interpolate <- stats::var(sapply(specList, get.length)) > 0

  if (!interpolate) {
    # equal lengths but check if frequency axes are different
    freqs <- sapply(specList, function(x) {x$freq})
    interpolate <- !all(apply(freqs, 1, stats::var) == 0.)
    # later, one could introduce some tolerance here for the deviation between
    # the individual frequency axes...
  }

  specList.interpolated <- specList

  if (interpolate) {

    # use the longest run for the reference spectrum
    freqRef <- seq(from = min(unlist(lapply(specList, get.fstart.existing))),
                   to = max(unlist(lapply(specList, get.fend.existing))),
                   by = min(unlist(lapply(specList, get.df))))

    for (i in 1:length(specList)) specList.interpolated[[i]] <- SpecInterpolate(freqRef,
    specList[[i]])

  }

  # Build matrices from input data

  ns <- length(specList.interpolated)
  nf <- get.length(specList.interpolated[[1]])

  specMatrix <- sapply(specList.interpolated, function(x) {x$spec})
  dofMatrix <- sapply(specList.interpolated, function(x) {x$dof})

  # weights at each frequency
  weightMatrix <- matrix(weights, nrow = nf, ncol = ns, byrow = TRUE)
  # missing spectral estimates
  missingObs <- is.na(specMatrix)

  # put weights to NA at missing estimates
  weightMatrix[missingObs] <- NA

  # normalise weights across spectra
  weightMatrix <- t(apply(weightMatrix, 1, function(x) {
    x / sum(x, na.rm = TRUE)}))

  # mean weighted spectra
  result <- list()
  result$freq <- specList.interpolated[[1]]$freq
  result$spec <- rowSums(specMatrix * weightMatrix, na.rm = TRUE)
  result$dof <- rowSums(dofMatrix, na.rm = TRUE)
  result$nRecord <- rowSums(!missingObs)
  class(result) <- "spec"

  return(AddConfInterval(result))
}


# Helper functions


remove.lowestFreq <- function(spec, iRemove) # Remove lowest frequencies from a spectra (as they are
# biased by detrending and the MTM estimator)
{
  if (iRemove == 0)
    index = seq(spec$spec) else index <- (-(1:iRemove))
  spec$spec <- spec$spec[index]
  spec$freq <- spec$freq[index]
  spec$dof <- spec$dof[index]
  return(spec)
}

# Helper functions to access the list elements to determine
# the target frequency discretisation
get.df <- function(x) return(mean(diff(x$freq)))
get.fend.existing <- function(x) return(max(x$freq[!is.na(x$spec)]))
get.fstart.existing <- function(x) return(min(x$freq[!is.na(x$spec)]))

# helper function to obtain length of frequency axis
get.length <- function(x) return(length(x$freq))

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
