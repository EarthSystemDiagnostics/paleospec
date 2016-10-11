##' @title  Interpolates the spectrum spec to the specRef frequency resolution
##' @param freqRef  frequency vector of the target resolution
##' @param spec list(spec,freq,dof)
##' @return one spectum as list(spec,freq,dof) (spec on the specRef resolution)
##' @author Thomas Laepple
##' @export
SpecInterpolate<-function(freqRef,spec)
{
    
 result<-list()
 result$freq<-freqRef


 result$spec<-approx(spec$freq,spec$spec,freqRef)$y
 result$dof<-approx(spec$freq,spec$dof,freqRef)$y

 class(result)<-"spec"
 return(result)
}


