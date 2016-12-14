##'  Add confidence intervals to a spectrum
##'  
##' 
##' @title  Add confidence intervals to a spectrum
##' @param spec spectrum list(spec,freq,dof)
##' @param MINVALUE Minimum value to which the confidence interval is limited
##' @param pval Interval from (pval/2 to 1-pval/2) is constructed
##' @return spectrum as the input but including lim.1 and lim.2 as new list elements
##' @author Thomas Laepple
##' @examples
##' N.R=1000
##' N.T=100
##' save.spec<-matrix(NA,N.T/2,N.R)
##' for (i.R in 1:N.R) {
##' save.spec[,i.R]<-SpecMTM(ts(SimPowerlaw(1, N.T)))$spec
##' }
##'
##' q.empirical<-apply(save.spec,1,quantile,c(0.025,0.975))
##' testspec<-SpecMTM(ts(SimPowerlaw(1, N.T)))
##' LPlot(AddConfInterval(testspec),ylim=c(0.05,10))
##' lines(testspec$freq,q.empirical[1,],col="red")
##' lines(testspec$freq,q.empirical[2,],col="red")
##' legend("bottomleft",lwd=2,col=c("black","red"),
##' c("one realization with chisq conf intervals","MC confidence intervals"))    
##' 
##' @export
AddConfInterval <- function(spec,MINVALUE = 1e-10,pval = 0.05)
{
    spec$lim.1 <- spec$spec*qchisq(c(1-pval/2),spec$dof)/(spec$dof)
    spec$lim.2 <- spec$spec*qchisq(c(pval/2),spec$dof)/(spec$dof)
    spec$lim.1[spec$lim.1 < MINVALUE] <- MINVALUE
    spec$lim.2[spec$lim.2 < MINVALUE] <- MINVALUE
    return(spec)
}



