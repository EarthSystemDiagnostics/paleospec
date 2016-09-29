##' Apply a filter to a timeseries
##' 
##' Using endpoint constrains as describen in  Mann et al., GRL 2003
##' minimum norm constraint (method=1)
##' minimum slope constraint (method=2)
##' minimum roughness constraint (method=3)
##' @title  Apply a filter to a timeseries
##' @param data Input timeseries (ts object)
##' @param filter vector of filter weights
##' @param method constraint method choice 1-3
##' @return filtered timeseries (ts object)
##' @author Thomas Laepple 
ApplyFilter <- function(data,filter,method=1)
{
N<-floor(length(filter)/2)
	if (method == 1)  #Minimum Norm
	{
		before<-rep(mean(data),N)
		after<-rep(mean(data),N)
	}
	if (method == 2)
	{
		before<-c(data)[N:1]
		after<- c(data)[length(data):(length(data)-N+1)]
	}
	if (method == 3)
	{
		before<-c(data)[N:1]
		after<- c(data)[length(data):(length(data)-N+1)]

		before<-c(data)[1]-(before - mean(before))

                after<-c(data)[length(data)]-(after - mean(after))



	}

	result<-filter(c(before,data,after),filter,circular=F)[(N+1):(N+length(data))]
	return(ts(result,frequency=frequency(data)))
}
