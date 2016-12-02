
##' average spectra with weighting, spectra can have different resolution and span a different freq range, in this case, they get linearly interpolate to the same (maximum range and highes resolution) frequency resolution first.
##' @title average spectra, optionally with weighting and interpolation
##' @param speclist list of spectra
##' @param iRemoveLowest number of lowest frequencies to remove (e.g. to remove detrending bias)
##' @param weights vector of weights (same length as elements in speclist)
##' @return list(spec,nRecords) spec=average spectrum, nRecords = number of records contributing to each spectral estimate
##' @author Thomas Laepple
##' @example
##' spec1<-SpecMTM(ts(rnorm(100),deltat=1))
##' spec2<-SpecMTM(ts(rnorm(100),deltat=1))
##' spec3<-SpecMTM(ts(rnorm(100),deltat=2))
##' spec.avg.12<-MeanSpectrum(list(spec1,spec2))
##' spec.avg.123<--MeanSpectrum(list(spec1,spec2,spec3)) 
##' @export
##' 
MeanSpectrum<-function(specList,iRemoveLowest=0,weights=rep(1,length(specList)))
{

#### Average the  spectra together
    print("Average spectra")
    if (length(weights) != length(specList)) stop("specList and weights have a different number of elements")
    weights<-weights/sum(weights) #
    #Remove the lowest/biased frequencies from the spectral estimates
    specList<-lapply(specList,remove.lowestFreq,iRemove=iRemoveLowest)

                                        #Check if all the spectra have the same frequency vector
    df<-sapply(specList,get.df)
    fStart<-sapply(specList,get.fstart.existing)
    fEnd<-sapply(specList,get.fend.existing)

    if (all(c(df==mean(df),fStart==mean(fStart),fEnd==mean(fEnd)))) {
        specList.interpolated=specList
    } else {       #Use the longest run for the reference spectrum and interpolate
        freqRef<-seq(from=min(fStart),to=max(fEnd),by=min(df))
        specList.interpolated<-list()
        for (i in 1:length(specList)) specList.interpolated[[i]]<-SpecInterpolate(freqRef,specList[[i]])
        for (i in 1:length(specList)) specList.interpolated[[i]]$spec<-specList.interpolated[[i]]$spec*weights[i]
    }


    NSpectra<-length(specList.interpolated)

    result<-list(freq=specList.interpolated[[1]]$freq,spec=rep(0,length(specList.interpolated[[1]]$spec)))
    specMatrix<-matrix(NA,NSpectra,length(specList.interpolated[[1]]$spec))
    dofMatrix<-matrix(NA,NSpectra,length(specList.interpolated[[1]]$spec))

    for (i in 1:length(specList.interpolated))
        {
            if (sum((result$freq-specList.interpolated[[i]]$freq)^2) > 0.1) stop("Different spectra length or resolutions")
            specMatrix[i,]<-specList.interpolated[[i]]$spec
            dofMatrix[i,]<-specList.interpolated[[i]]$dof    
        }

    result$spec<-colSums(specMatrix,na.rm=TRUE)
    nRecord=colSums(!is.na(specMatrix))
    result$dof<-colSums(dofMatrix,na.rm=TRUE)
    class(result)<-"spec"
    
    return(list(spec=AddConfInterval(result),nRecord=nRecord)) 
}

spec1<-SpecMTM(ts(rnorm(100),deltat=1))
spec2<-SpecMTM(ts(rnorm(100),deltat=1))
spec3<-SpecMTM(ts(rnorm(100),deltat=2))
spec.avg.12<-MeanSpectrum(list(spec1,spec2))
spec.avg.123<-MeanSpectrum(list(spec1,spec2,spec3))

LPlot(spec.avg.123$spec)
plot(spec.avg.123$spec$freq,spec.avg.123$nRecord)

                                        #Helper functions


remove.lowestFreq<-function(spec,iRemove)
                                        #Remove lowest frequencies from a spectra (as they are biased by detrending and the MTM estimator)
{
    if (iRemove==0) index=seq(spec$spec) else index<-(-(1:iRemove))
    spec$spec<-spec$spec[index]
    spec$freq<-spec$freq[index]    
    spec$dof<-spec$dof[index]
    return(spec)
}

                                        #Helper functions to access the list elements to determine the target frequency discretisation
get.df<-function(x) return(mean(diff(x$freq)))
get.fend.existing<-function(x) return(max(x$freq[!is.na(x$spec)]))
get.fstart.existing<-function(x) return(min(x$freq[!is.na(x$spec)]))


                                        #f (bNormalize) Normalize to the mean of the common interval
                                        #
                                        #   fmin<-min(unlist(lapply(temp,get.fend.existing)))
                                        #   fmax<-max(unlist(lapply(temp,get.fstart.existing)))
                                        #   i.min<-closest.element(freqRef,fmin)
                                        #   i.max<-closest.element(freqRef,fmax)
                                        #
                                        #   var.band<-vector()
                                        #   for (i in 1:length(temp))
                                        #   {
                                        #       var.band[i]<-mean(temp[[i]]$spec[i.min:i.max])
                                        #   }
                                        #
                                        #   rescale<-mean(var.band)/var.band
                                        #   if (!is.null(weights)) rescale=weights
                                        #   for (i in 1:length(temp)) temp[[i]]$spec<-temp[[i]]$spec*rescale[i]
                                        #
