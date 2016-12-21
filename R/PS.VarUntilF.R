

##' Integral of PSD=f^(-beta) from f1=1/N to f2=f
##' this equals the variance of a lowpass filtered powerlaw process
##' WARNING: The result is not normalized
##' @title Variance of a powerlaw process if integrated from until frequency f
##' @param f frequency until which to integrate
##' @param beta powerlaw slope 
##' @param N length of the timeseries
##' @return non-normalized variance
##' examples
##' beta <- 1
##' signal <- ts(SimPowerlaw(beta,100000))
##' spec <- SpecMTM(signal)
##' v1 <- GetVarFromSpectra(spec,f=c(1/length(signal),0.5))
##' v2 <- GetVarFromSpectra(spec,f=c(1/length(signal),0.01))
##' PS.VarUntilF(0.01,beta,length(signal))/PS.VarUntilF(0.5,beta,length(signal))
##' v2$var/v1$var
##' @author Thomas Laepple
##' @export
PS.VarUntilF<-function(f,beta,N)
    {
        if (beta==1) return(log(f)-log(1/N))
        else return((f^(1-beta)-(1/N)^(1-beta))/(1-beta))
    }
