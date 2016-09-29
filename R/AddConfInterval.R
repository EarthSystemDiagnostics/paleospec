##'  Add confidence intervals to a spectrum
##'
##' 
##' @title  Add confidence intervals to a spectrum
##' @param spec spectrum list(spec,freq,dof)
##' @param MINVALUE Minimum value to which the confidence interval is limited
##' @param pval Interval from (pval/2 to 1-pval/2) is constructed
##' @return spectrum as the input but including lim.1 and lim.2 as new list elements
##' @author Thomas Laepple 
AddConfInterval<-function(spec,MINVALUE=1e-10,pval=0.05)
{
    spec$lim.1<-spec$spec*qchisq(c(1-pval/2),spec$dof)/(spec$dof)
    spec$lim.2<-spec$spec*qchisq(c(pval/2),spec$dof)/(spec$dof)
    spec$lim.1[spec$lim.1<MINVALUE]<-MINVALUE
    spec$lim.2[spec$lim.2<MINVALUE]<-MINVALUE
    return(spec)
}
