
##' The timeseries have an expected variance of 1 and an expected correlation of r
##' 
##' @title Create a pair of random signals with powerlaw signal and powerlaw noise
##' @param N Number of points per timeseries
##' @param betaSignal powerlaw slope of the signal
##' @param betaNoise powerlaw slope of the noise
##' @param r expected correlation between both vectors
##' @return list containing both vectors y1 and y2
##' @examples
##' mean(replicate(1000,{test <- SimulatePowerlawSignalPair(200,1,1,0.5);cor(test$y1,test$y2)}))
##' @author Thomas Laepple
##' @export
SimulatePowerlawSignalPair<-function(n,beta.signal,beta.noise,r)
    {

        signal<-SimPowerlaw(beta.signal,n)
        noise<-(SimPowerlaw(beta.noise,n)*sqrt((1-r^2)/(r^2))) 
        
        y1<-signal
        y2<-(signal+noise)*r  #Normalize by r to get variance 1
            
        return(list(y1=y1,y2=y2))
    }
