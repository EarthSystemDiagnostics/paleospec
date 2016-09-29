##' Derives and plots the transfer function (given a filter)
##'
##' Get the transfer function of a symetric filter, page 122 in Bloomfield 1976, 
##' @title Derives and plots the transfer function
##' @param g.u 
##' @param resolution 
##' @param bPlot 
##' @param add 
##' @param ... 
##' @return list(omega,y) containing the transfer function
##' @author Thomas Laepple 
GetTransferFunction<-function(g.u,resolution=100,bPlot=TRUE,add=FALSE, ...)
{
    n<-length(g.u)
    n.side<-(n-1)/2
    omega=(1:resolution)*pi/resolution
    yt<-rep(0,length(omega))

    for (u in (-1*n.side):n.side)
	yt<-yt+(g.u[u+n.side+1]*exp(-1*1i*omega*u))

    if (bPlot)
        {
            if (!add) plot(omega/2/pi,abs(yt)^2,type="l",xlab="omega",main="Transfer function", ...)
            else lines(omega/2/pi,abs(yt)^2,col="blue", h...)
        }

    return(list(omega=omega/2/pi,y=abs(yt)^2))
}
