##' add Logplot + transparent confidence interval for the spectral plotting
##'
##' 
##' @title add Logplot + transparent confidence interval for the spectral plotting
##' @param x  spectra object
##' @param conf TRUE: Plot confidence interval
##' @param bPeriod if \code{TRUE} treat the x-axis values in units of period
##'     (inverse frequency). Defaults to \code{FALSE}.
##' @param col color
##' @param alpha  transparency
##' @param removeFirst omit removeFirst values on the low frequency side
##' @param removeLast omit removeFirst values on the high frequency side
##' @param ... other parameters to be passed to the line function
##' @return none
##' @examples
##' x<-ts(arima.sim(list(ar = 0.9),1000))
##' spec<-SpecMTM(x)
##' LPlot(spec,col="grey")
##' LLines(LogSmooth(spec),lwd=3)
##' @author Thomas Laepple
##' @export
LLines<-function(x,conf=TRUE,bPeriod=FALSE,col="black",alpha=0.3,removeFirst=0,removeLast=0,...)
{

    if (bPeriod) x$freq=1/x$freq
    
    index<-(removeFirst+1):(length(x$freq)-removeLast)  
    x$freq<-x$freq[index]
    x$spec<-x$spec[index]
    x$lim.1<-x$lim.1[index]
    x$lim.2<-x$lim.2[index]

   if (conf) polygon(c(x$freq,rev(x$freq)),c(x$lim.1,
        rev(x$lim.2)),col=ColTransparent(col,alpha),border=NA)
    lines(x$freq,x$spec,col=col,...)
}


