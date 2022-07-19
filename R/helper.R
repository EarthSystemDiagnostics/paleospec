# -------------------------------------------------------------
#
# various unexported helper functions for the PaleoSpec package
#
# -------------------------------------------------------------

# -------------------------------------------------------------
# Documented functions

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

# -------------------------------------------------------------
# Undocumented functions

# Check for confidence interval
#
# @param spec spectral object.
# @return logical value.
# @author Thomas Münch
has.limits <- function(spec) {
  !is.null(spec$lim.1) & !is.null(spec$lim.2)
}
# Average frequency discretization
#
# @param x spectral object
# @return numerical value.
# @author Thomas Münch
get.df <- function(x) {
  mean(diff(x$freq))
}
# Frequency of highest existing spectral estimate
#
# @param x spectral object
# @return numerical value.
# @author Thomas Münch
get.fend.existing <- function(x) {
  max(x$freq[!is.na(x$spec)])
}
# Frequency of lowest existing spectral estimate
#
# @param x spectral object
# @return numerical value.
# @author Thomas Münch
get.fstart.existing <- function(x) {
  min(x$freq[!is.na(x$spec)])
}
# Number of spectral observations
#
# @param x spectral object
# @return numerical value.
# @author Thomas Münch
get.length <- function(x) {
  length(x$freq)
}
# Extract frequency vector
#
# @param x spectral object
# @return numerical vector of frequency axis.
# @author Thomas Münch
get.freq <- function(x) {
  x$freq
}
# Extract PSD vector
#
# @param x spectral object
# @return numerical vector of power spectral density values.
# @author Thomas Münch
get.spec <- function(x) {
  x$spec
}
# Extract DOF vector
#
# @param x spectral object
# @return numerical vector of degrees of freedom.
# @author Thomas Münch
get.dofs <- function(x) {
  x$dof
}
