##' @title 
##' @param x 
##' @return 
##' @author Thomas Laepple 
NaFillTs<-function(x) {
    if (sum(is.na(x))==0) return(x)
    nonMissing<-!is.na(x)
    temp<-approx(c(time(x))[nonMissing],c(x)[nonMissing],c(time(x)),rule=2)$y
     return(ts(temp,start=start(x),deltat=deltat(x)))
               }
