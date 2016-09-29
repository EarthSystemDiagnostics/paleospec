##' Confidence Interval of ratios
##' based on a ChiSquare Distribution
##'
##' 
##' @title Confidence Interval of ratios
##' @param varratio 
##' @param df.1 degree of freedom of denominator
##' @param df.2 degree of freedom of numerator
##' @param pval 
##' @return lower and upper confidence intervals
##' @author Thomas Laepple 
ConfRatio<-function(varratio,df.1,df.2,pval=0.1)
{
     return(varratio*qf(c(pval/2,(1-pval/2)),df1=df.1,df2=df.2))
}
