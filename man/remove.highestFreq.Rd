% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper.R
\name{remove.highestFreq}
\alias{remove.highestFreq}
\title{Remove high frequencies}
\usage{
remove.highestFreq(spec, iRemove)
}
\arguments{
\item{spec}{a list object of class \code{"spec"}.}

\item{iRemove}{integer; number of highest frequencies to remove.}
}
\value{
the input data with the data at the \code{iRemove} highest
  frequencies removed.
}
\description{
Remove the highest frequencies from a spectra, since they may be biased by
the MTM estimator.
}
\author{
Thomas Münch
}
