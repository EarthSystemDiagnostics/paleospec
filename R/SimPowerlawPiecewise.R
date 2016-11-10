##' @title #Simulate a timeseries with length N which has a spectra consisting of two powerlaws
##' @param beta1 slope for frequencies lower than breakpoint
##' @param beta2 slope for frequencies higher than breakpoint
##' @param N Number of points to simulate
##' @param deltat  timestep of the timeseries
##' @param breakpoint frequency of the breakpoint
##' @return N random numbers drawn according to the piecewise powerlaw PSD
##' @export
##' @author Thomas Laepple 
SimPowerlawPiecewise<-function(beta1,beta2,N,deltat=1,breakpoint=1/50)
{
#Simulate a timeseries with length N which has a spectra consisting of two powerlaws
#beta1 for frequencies lower than breakpoint
#beta 1 for frequencies higher than the breakpoint
#freq = 1/dt
    Norg<-N
    N<-ceiling(N/2)*2
    df  = 1/(N);
    f=seq(from=df,to=1/(2),by=df)
    index.breakpoint<-closest.element(f*(1/deltat),breakpoint)

    #Intercept to cross at breakpoint
    a<- (-1*beta1*log(f[index.breakpoint]))-(-1*beta2*log(f[index.breakpoint]))
    FilterLow <- exp((-1)*log(f[1:index.breakpoint])*beta1-a)
    FilterHigh <- exp((-1)*log(f[-1*(1:index.breakpoint)])*beta2)

    Filter<-c(sqrt(FilterLow),sqrt(FilterHigh))
    Filter = c(max(Filter), Filter,rev(Filter))

    x   = scale(rnorm(N+1,1))
    fx  =fft(x)
    ffx =fx*Filter;
    result<-scale(Re(fft(ffx,inverse=TRUE)))
    return(result[1:Norg])
}


