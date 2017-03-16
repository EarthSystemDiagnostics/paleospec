## Aim: Test filter, inverse filter and power spectra
## Status/belongs to: research process of spectral uncertainty paper


library(PaleoSpec)
library(ecustools)

BETA <- 1
ts.input <- SimPowerlaw(BETA,100000)
filter.weights <- ImpulseResponse.BergerHeath(seq(from=-200,to=200),d=40)

n<-length(filter.weights)
directtransfer<-rep(0,n)
directtransfer[(n+1)/2]<-1

filter.weights.inv <- directtransfer-filter.weights

ts.filtered <- ApplyFilter(ts.input,filter.weights)
ts.remain <-   na.omit(ts.input - ts.filtered)
ts.filtered <- na.omit(ts.filtered)


transfer     <- GetTransferFunction(filter.weights,resolution=100000)
transfer.inv <- GetTransferFunction(filter.weights.inv,resolution=100000)



freq <- transfer$omega
index <- (freq>=(1/length(ts.input))) #only look at frequencies beyond the Raleigh frequency 
#if not the normalization changes the result
freq <- freq[index]
transfer$y <- transfer$y[index]
transfer.inv$y <- transfer.inv$y[index]
psd <- AnPowerlaw(BETA,freq)



LPlot(SpecMTM(ts.filtered,k=40,nw=20),ylim=c(1e-3,5000))
LLines(SpecMTM(ts.remain,k=40,nw=20),col="blue")
LLines(SpecMTM(ts.input,k=40,nw=20),col="grey")

lines(freq,psd*transfer$y,col="red",lwd=3)
lines(freq,psd*(1-transfer$y),col="green",lwd=3)
lines(freq,psd*(transfer.inv$y),col="cyan",lwd=3)
lines(freq,psd,col="grey",lwd=3)

legend("bottomleft",col=c("grey","blue","black"),lwd=1,c("Input","Input-Lowpass(Input)","Lowpass(Input)"),bty="n")

legend("bottom",col=c("grey","green","cyan","red"),lwd=3,c("Input theoretical","PSD(Input)*(1-T_Lowpass)","PSD(Input)*(T_Highpass)","PSD(Input)*T_Lowpass"),bty="n")

