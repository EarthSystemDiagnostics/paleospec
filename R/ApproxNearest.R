
##' extends approx which always takes the right or left neighbour 
##' or the weighted mean between both if f>0<1
##' 
##' @title approximate a timeseries using the nearest neighbour
##' @param x numeric vector giving the coordinates of the points to be interpolate
##' @param y corresponding y values
##' @param xout set of numeric values specifying where interpolation is to take place.
##' @param rule an integer (of length 1 or 2) describing how interpolation is
##' to take place outside the interval see ?approx
##' @return  a list with components 'x' and 'y', containing
##' length(xouth) coordinates which interpolate the given data points 
##' @author Thomas Laepple
##' @examples
##' x<-1:10
##' y<-1:10
##' xout<-seq(from=0,to=11,by=0.01)
##' plot(x,y,type="b",pch=19,xlim=range(xout))
##' result<-ApproxNearest(x,y,xout)
##' lines(result,col="red")
##' @export
ApproxNearest <- function(x, y, xout,rule = 1)
{
    result <- list()
    result$x = xout
    result$y = approx(c(x[1], x + c(diff(x)/2, 0)), c(y[1], y), 
        xout, method = "constant", f = 1,rule=rule)$y
    return(result)
}

