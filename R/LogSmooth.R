#' @title Smooths the spectrum using a log smoother
#' @param spectra spectra: list(spec,freq) spec[specIndex]: spectra density
#'   vector freq[specIndex]: frequency vector
#' @param df.log width of the smoother in log units
#' @param removeFirst elements to remove on the slow side (one element
#'   recommended because of the detrending)
#' @param removeLast elements to remove on the fast side
#' @param bLog TRUE: average in the log space of the power, FALSE: arithmetic
#'   average
#' @return smoothed spectrum
#' @examples
#' x<-ts(arima.sim(list(ar = 0.9),1000))
#' spec<-SpecMTM(x)
#' LPlot(spec,col='grey')
#' LLines(LogSmooth(spec,df.log=0.01),lwd=2,col='green')
#' LLines(LogSmooth(spec,df.log=0.05),lwd=2,col='blue')
#' LLines(LogSmooth(spec,df.log=0.1),lwd=2,col='red')
#' legend('bottomleft', col=c('grey','green','blue','red'),
#' lwd=2,c('raw','smoothed 0.01',
#'  'smoothed 0.05', 'smoothed 0.1'), bty='n')
#'
#' # Removal of lower (first) and higher (last) frequencies
#' LPlot(spec,col='grey')
#' LLines(LogSmooth(spec,df.log=0.01, removeFirst = 1),lwd=2,col='green')
#' LLines(LogSmooth(spec,df.log=0.05, removeLast = 20),lwd=2,col='blue')
#' LLines(LogSmooth(spec,df.log=0.1, removeFirst = 3, removeLast = 20),lwd=2,col='red')
#' @author Thomas Laepple
#' @export
LogSmooth <- function(spectra, df.log = 0.05,
                      removeFirst = 0, removeLast = 0,
                      bLog = FALSE) {

  result <- list()
  removeFirst <- removeFirst+1
  nfreq <- length(spectra$spec)
  result$freq <- spectra$freq[removeFirst:nfreq]

  if (bLog) {

    temp <- smoothlog.cutEnd(log(spectra$spec[removeFirst:nfreq]),
                             result$freq, df.log,
                             dof = spectra$dof[removeFirst:nfreq])

    temp$spec <- exp(temp$spec)

  } else {

    temp <- smoothlog.cutEnd(spectra$spec[removeFirst:nfreq],
                             result$freq, df.log,
                             dof = spectra$dof[removeFirst:nfreq])
  }

  result$spec <- temp$spec[1:(length(temp$spec) - removeLast)]
  result$freq <- result$freq[1:(length(result$freq) - removeLast)]
  temp$dof <- temp$dof[1:(length(temp$dof) - removeLast)]

  ### Calculate the confidence intervals
  result$dof <- temp$dof
  result <- AddConfInterval(result)

  #######

  class(result) <- "spec"
  return(result)
}



# Helper functions to smooth power spectra

#' weights
#'
#' @return
#' @keywords internal
weights <- function(x, sigma) {
  1/sqrt(2 * pi * sigma^2) * exp(-x^2/(2 * sigma^2))
  }

#' @title weights
#' @return weight vector
#' @author Thomas Laepple
#' @keywords internal
fweights <- function(ftarget, f, df.log) {
  sigma <- ftarget * (exp(df.log) - exp(-df.log))
  return(weights(f - ftarget, sigma))
}

#' @title fweights.lin
#' @return  weight vector
#' @author Thomas Laepple
#' @keywords internal
fweights.lin <- function(ftarget, f, df.log) {
  sigma <- df.log
  return(weights(f - ftarget, sigma))
}


#' @title smoothlog
#' @return  smoothed x
#' @author Thomas Laepple
#' @keywords internal
smoothlog <- function(x, f, df.log) {
  x.smooth <- vector()
  for (i in 1:length(f)) {
    w <- fweights(f[i], f, df.log)
    x.smooth[i] <- sum(x * (w/sum(w)))
  }
  return(x.smooth)
}



#' @title smoothlog.cutEnd
#' @return  smoothed x
#' @author Thomas Laepple
#' @keywords internal
smoothlog.cutEnd <- function(x, f, df.log, dof = 1) {
  x.smooth <- vector()
  dof.smooth <- vector()
  for (i in 1:length(f)) {
    w <- fweights(f[i], f, df.log)
    ## Cut the weights down
    DistanceSlowEnd <- i - 1
    DistanceFastEnd <- length(f) - i

    if ((i + DistanceSlowEnd + 1) <= length(f))
      w[(i + DistanceSlowEnd + 1):length(f)] <- 0
    if ((i - DistanceFastEnd - 1) >= 1)
      w[1:(i - DistanceFastEnd - 1)] <- 0
    w <- w/f
    w <- w/sum(w)  #normalize to 1

    x.smooth[i] <- sum(x * w)
    dof.smooth[i] <- sum(w * dof)/sum(w^2)
  }
  return(list(spec = x.smooth, dof = dof.smooth))
}



#' @title smoothlin.cutEnd
#' @return  smoothed x
#' @author Thomas Laepple
#' @keywords internal
smoothlin.cutEnd <- function(x, f, df.log, dof = 1) {
  x.smooth <- vector()
  dof.smooth <- vector()
  for (i in 1:length(f)) {
    w <- fweights.lin(f[i], f, df.log)
    ## Cut the weights down
    DistanceSlowEnd <- i - 1
    DistanceFastEnd <- length(f) - i

    if ((i + DistanceSlowEnd + 1) <= length(f))
      w[(i + DistanceSlowEnd + 1):length(f)] <- 0
    if ((i - DistanceFastEnd - 1) >= 1)
      w[1:(i - DistanceFastEnd - 1)] <- 0
    w <- w/f
    w <- w/sum(w)  #normalize to 1
    x.smooth[i] <- sum(x * w)
    dof.smooth[i] <- sum(w * dof)/sum(w^2)
  }
  return(list(spec = x.smooth, dof = dof.smooth))
}

