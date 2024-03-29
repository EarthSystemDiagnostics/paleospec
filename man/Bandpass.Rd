% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Bandpass.R
\name{Bandpass}
\alias{Bandpass}
\title{Calculate  Weights for a Bandpass Filter}
\usage{
Bandpass(omega.upper, omega.lower, n, sample = 1, convergence = T)
}
\arguments{
\item{omega.upper}{upper cutoff frequency}

\item{omega.lower}{lower cutoff frequency}

\item{n}{length of the filter, has to be odd}

\item{sample}{sampling rate of the timeseries on which the filter will be applied (1/deltat)}

\item{convergence}{TRUE: smoothed least square lowpass; FALSE = unsmoothed}
}
\value{
vector of filter weights
}
\description{
Derive the (smoothed) least square bandpass based on Bloomfield 1976
}
\author{
Thomas Laepple
}
