
#' @title Construct the inverse filter in the time domain
#' @param filter.weights 
#' @return filter weights for the inverse filter
#' @author Thomas Laepple
#' @export
InverseFilter <- function(filter.weights) {
     if ((length(filter.weights)%%2) == 0) 
        stop("N must be odd")
    n <- length(filter.weights)
    directtransfer <- rep(0, n)
    directtransfer[(n + 1)/2] <- 1
    filter.weights.inv <- directtransfer - filter.weights
    return(filter.weights.inv)
}

