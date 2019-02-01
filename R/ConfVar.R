#' @title Provide ChiSquared confidence intervals for ratios
#' @param varlist list(var,dof)
#' @param pval requested p-value
#' @return Output: confidence intervals
#' @author Thomas Laepple
ConfVar<-function(varlist,pval=0.05)
{

 lim.1<-varlist$var*qchisq(c(1-pval/2),varlist$dof)/(varlist$dof)
 lim.2<-varlist$var*qchisq(c(pval/2),varlist$dof)/(varlist$dof)
 return(list(lim.1=lim.1,lim.2=lim.2))
}
