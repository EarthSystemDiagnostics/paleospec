##' Variance estimate by integrating a part of the spectrum
##'
##' 
##' @title Variance estimate by integrating a part of the spectrum
##' @param spec spectrum (list of spec,freq,dof) to be analysed
##' @param f f[1],f[2]: frequency interval to be analysed
##' @param dfreq frequency discretisation used in the temporary  interpolation
##' @param df.log if > 0, smooth the spectra prior to integrating
##' @param bw the bandwidth assumed for the confinterval calculation (from the multitaper spectral estimate)
##' @return list(var,dof)  variance and corresponding dof
##' @author Thomas Laepple
##' @examples
##' x<-ts(rnorm(100))
##' spec<-SpecMTM(x)
##' var(x) #Sample variance of the timeseries
##' GetVarFromSpectra(spec,c(1/100,0.5))
##' GetVarFromSpectra(spec,c(0.25,0.5))
##' @export
GetVarFromSpectra <- function(spec,f,dfreq=NULL,df.log = 0,bw = 3)
{
    
    if (f[1] >= f[2]) stop("f1 must be < f2")
    freqVector <- spec$freq
  
     ## Test it both frequencies are included 
    if (f[1] < FirstElement(freqVector)) {
        warning("f[1] is smaller than the lowest frequency in the spectrum, set to the lowest frequency")
        f[1] <- FirstElement(freqVector)
        }
     if (f[2] > LastElement(freqVector))  {
        warning("f[2] is larger than the highest frequency in the spectrum, set to the highest frequency")
        f[2] <- LastElement(freqVector)
     }

    if (is.null(dfreq)) dfreq <- min(diff(spec$freq)[1]/5, (f[2]-f[1])/100)
    newFreq <- seq(from = f[1],to = f[2],by = dfreq) # temporary frequency vector

    ## For spectra fromn the periodogram, the dof are supplied as a scalar named df; for MTN as a vector called DOF
    if (is.null(spec$dof)) spec$dof <- rep(spec$df,length(spec$freq))
        
    
     dof.original <- mean(SpecInterpolate(newFreq,spec)$dof)
     ## DOF before smoothing
    
    if (df.log > 0) spec <- LogSmooth(spec,removeLast = 0,df.log = df.log)
    
     vars <- mean(SpecInterpolate(newFreq,spec)$spec)*diff(f)*2
     dof <- mean(SpecInterpolate(newFreq,spec)$dof)

     ## Estimate the DOF by calculating how many independent
     ## spectral estimates contribute to the calculated mean value
     dfreq <- mean(diff(spec$freq))
     nSpecEstimate <- (f[2]-f[1])/dfreq
     dof <- dof*nSpecEstimate/bw
    
     if (dof < dof.original) dof = dof.original

     return(list(var = vars,dof = dof))
}

