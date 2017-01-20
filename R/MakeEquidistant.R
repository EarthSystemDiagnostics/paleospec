##' Average an irregular timeseries to a regular timeseries
##'
##' Make an irregular timeseries equidistant by
##' interpolating to high resolution, lowpass filtering to the Nyquist
##' frequency, and subsampling; e.g. as used in Huybers and Laepple, EPSL 2014
##' @title Average an irregular timeseries to a regular timeseries
##' @param t.x vector of timepoints
##' @param t.y vector of corresponding values
##' @param dt target timestep; can be omitted if time.target is supplied
##' @param time.target time vector to which timeseries should be averaged/interpolated to
##' by default the same range as t.x with a timestep dt
##' @param dt.hres timestep of the intermediate high-resolution interpolation.
##' Should be smaller than the smallest timestep
##' @param bFilter (TRUE) low passs filter the data to avoid aliasing, (FALSE) just interpolate
##' @param k scaling factor for the Length of the filter (increasing creates
                                        #a sharper filter, thus less aliasing)
##' @param kf  scaling factor for the lowpass frequency; 1 = Nyquist, 1.2 =
##' 1.2xNyquist is a tradeoff between reducing variance loss and keeping
##' aliasing small
##' @param method.interpolation
##' "linear" or "constant", see approx
##' @param method.filter
##' To avoid loosing data at the ends of the dataset, endpoint constrains are used (see ApplyFilter)
##' no constraint (loss at both ends) (method=0), only works if t.x covers more time than time.target
##' minimum norm constraint (method=1)
##' minimum slope constraint (method=2)
##' minimum roughness constraint (method=3)
##' circular filtering (method=4) 
##' @return ts object with the equidistant timeseries
##' @author Thomas Laepple
##' @export
MakeEquidistant<-function(t.x,t.y,dt=NULL,time.target=seq(from=t.x[1],to=t.x[length(t.x)],by=dt),dt.hres=NULL,bFilter=TRUE,k=5,kf=1.2,method.interpolation="linear",method.filter=2)
{
    index<-!is.na(t.x)
    t.x<-t.x[index]
    t.y<-t.y[index]

    if (is.null(dt)) dt <- mean(diff(time.target))

    
    if (is.null(dt.hres)) #Choose dt.hres if not supplied
        {
            dt.hres<-dt/10

            #Get the minimum timestep not considering a zero timestep 
            dt.x<-diff(t.x)
            minTimeStep<-min(dt.x[dt.x>0])
            if (dt.hres>minTimeStep) dt.hres=minTimeStep
        }
    

    if (dt.hres > min(diff(t.x),na.rm=TRUE)) warning("dt.hres is lower
than the minimum timestep")

    index<-(!is.na(t.y))
    

    time.hres<-seq(from=first(t.x),to=last(t.x),by=dt.hres)
    data.hres<-approx(t.x[index],t.y[index],time.hres,method=method.interpolation)

    index<-!is.na(data.hres$y)
    data.hres$x<-data.hres$x[index]
    data.hres$y<-data.hres$y[index]

    filterLength=round(k*(dt/(2*dt.hres)))
    if ((filterLength %% 2)==0) filterLength=filterLength+1

    if (bFilter)
        {
            f.lowpass<-Lowpass(1/(2*dt)*kf,filterLength,sample=1/dt.hres)
            meanvalue<-mean(data.hres$y)
            data.hres.filtered<-ApplyFilter(data.hres$y-meanvalue,f.lowpass,method=method.filter)+meanvalue
        } else
            {
                data.hres.filtered<-data.hres$y
            }

    
    index<-!is.na(data.hres.filtered)
    data.target<-approx(data.hres$x[index],data.hres.filtered[index],time.target,method = method.interpolation)

    return(ts(data.target$y,start=first(data.target$x),deltat=diff(data.target$x)[1]))
}
