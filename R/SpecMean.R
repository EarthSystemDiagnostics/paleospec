SpecMean<-function(speclist,weight=FALSE,dof=TRUE)
{
# Returns the mean of the spectra
#
# Inputs:
# speclist[[]], each containing
# a list(spec,freq,dof)
#                     spec[specIndex]: spectra density vector
#                     freq[specIndex]: frequency vector
#                     dof[specIndex]: DOF ... or a single value df
# weight:             weight by the uncertainty
# Return:
# the mean spectra

#Weighting by 1/sigma^2;  sigma^2 = variance  is proportional to 1/DOF; This
#is true if we assume a constant spectral density... if not, higher spectral densities #have higher uncertainty

if (weight) warning("Uncertainty weighting uses DOF; has to be checked")

 index<-which(unlist(lapply(speclist,length))>0) ##Ignore empty spectra
 NSpectra<-length(speclist)

 result<-list(freq=speclist[[index[1]]]$freq,spec=rep(0,length(speclist[[index[1]]]$spec)))

 specMatrix<-matrix(NA,NSpectra,length(speclist[[index[1]]]$spec))
 if (dof) dofMatrix<-matrix(NA,NSpectra,length(speclist[[index[1]]]$spec))

 weightMatrix<-matrix(NA,NSpectra,length(speclist[[index[1]]]$spec))

 for (i in index)
 {
    if (sum((result$freq-speclist[[i]]$freq)^2) > 0.1) stop("Different spectra length or resolutions")

     if (weight) w<-(speclist[[i]]$dof) else w<-rep(1,length(speclist[[i]]$spec)) #weight by dof

     w[is.na(speclist[[i]]$spec)]<-NA #weights are NA in areas with missing spectral estimates
     specMatrix[i,]<-speclist[[i]]$spec*w

     if (dof) dofMatrix[i,]<-speclist[[i]]$dof*w
     weightMatrix[i,]<-w  #Keep track of the weights for the normalizing at the end

 }

 if (weight) wMean<-colMeans(weightMatrix,na.rm=TRUE)  else wMean<-1 #Mean of the weights

 result$spec<-colMeans(specMatrix,na.rm=TRUE)/wMean
 if (dof) result$dof<-colSums(dofMatrix,na.rm=TRUE)/wMean
 class(result)<-"spec"
 return(result)
}





