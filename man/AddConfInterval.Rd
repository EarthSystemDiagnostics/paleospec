% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AddConfInterval.R
\name{AddConfInterval}
\alias{AddConfInterval}
\title{Add confidence intervals to a spectrum}
\usage{
AddConfInterval(spec, MINVALUE = 1e-10, pval = 0.05)
}
\arguments{
\item{spec}{spectrum list(spec,freq,dof)}

\item{MINVALUE}{Minimum value to which the confidence interval is limited}

\item{pval}{Interval from (pval/2 to 1-pval/2) is constructed}
}
\value{
spectrum as the input but including lim.1 and lim.2 as new list elements
}
\description{
Add confidence intervals to a spectrum
}
\examples{
N.R=1000
N.T=100
save.spec<-matrix(NA,N.T/2,N.R)
for (i.R in 1:N.R) {
save.spec[,i.R]<-SpecMTM(ts(SimPowerlaw(1, N.T)))$spec
}

q.empirical<-apply(save.spec,1,quantile,c(0.025,0.975))
testspec<-SpecMTM(ts(SimPowerlaw(1, N.T)))
LPlot(AddConfInterval(testspec),ylim=c(0.05,10))
lines(testspec$freq,q.empirical[1,],col="red")
lines(testspec$freq,q.empirical[2,],col="red")
legend("bottomleft",lwd=2,col=c("black","red"),
c("one realization with chisq conf intervals","MC confidence intervals"))    

}
\author{
Thomas Laepple
}
