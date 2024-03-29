#' Get closest element of a vector
#'
#' @param xvector a vector of values
#' @param x the value to find the closest match to
#' @param type methods N, M or L
#'
#' @return a vector of length 1
#' @export
#'
#' @examples
#' a <- 1:10
#' ClosestElement(a, 3.4)
ClosestElement <- function(xvector, x, type = "N")
{
  if (type == "N")
    return(which.min(abs(x - xvector)))

  if (min(diff(xvector)) < 0)
    stop("Vector must be monotonically increasing for the the methods M and L")
  if (type == "M")
    return(first(which(x <= xvector)))
  if (type == "L")
    return(last(which(x >= xvector)))
}
