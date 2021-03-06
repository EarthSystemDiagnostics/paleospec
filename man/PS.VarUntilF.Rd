% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PS.VarUntilF.R
\name{PS.VarUntilF}
\alias{PS.VarUntilF}
\title{Variance of a powerlaw process if integrated from until frequency f}
\usage{
PS.VarUntilF(f, beta, N)
}
\arguments{
\item{f}{frequency until which to integrate}

\item{beta}{powerlaw slope}

\item{N}{length of the timeseries}
}
\value{
non-normalized variance
}
\description{
Integral of PSD=f^(-beta) from f1=1/N to f2=f
this equals the variance of a lowpass filtered powerlaw process
WARNING: The result is not normalized
}
\examples{
beta <- 1
signal <- ts(SimPowerlaw(beta,100000))
spec <- SpecMTM(signal)
v1 <- GetVarFromSpectra(spec,f=c(1/length(signal),0.5))
v2 <- GetVarFromSpectra(spec,f=c(1/length(signal),0.01))
PS.VarUntilF(0.01,beta,length(signal))/PS.VarUntilF(0.5,beta,length(signal))
v2$var/v1$var
}
\author{
Thomas Laepple
}
