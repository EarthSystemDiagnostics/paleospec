% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SimulatePowerlawSignalPair.R
\name{SimulatePowerlawSignalPair}
\alias{SimulatePowerlawSignalPair}
\title{Create a pair of random signals with powerlaw signal and powerlaw noise}
\usage{
SimulatePowerlawSignalPair(n, beta.signal, beta.noise, r)
}
\arguments{
\item{r}{expected correlation between both vectors}

\item{N}{Number of points per timeseries}

\item{betaSignal}{powerlaw slope of the signal}

\item{betaNoise}{powerlaw slope of the noise}
}
\value{
list containing both vectors y1 and y2
}
\description{
The timeseries have an expected variance of 1 and an expected
correlation of r
}
\examples{
mean(replicate(1000,{test <- SimulatePowerlawSignalPair(200,1,1,0.5);cor(test$y1,test$y2)}))
}
\author{
Thomas Laepple
}
