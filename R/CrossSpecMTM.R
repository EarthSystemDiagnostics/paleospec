#' @title MTM cross spectral estimator
#' @description estimates cross spectral densities using slepian tapers with adaptive weighting. Function is an adaption of Peter Huybers Matlab coherence function cmtm.m
#' @param x,y input data vectors (or array if y=NULL)
#' @param dt time stepping (default = 1)
#' @param NW number of slepian taper (default = 8)
#' @param k order (default = NULL)
#' @param dpssIN allows to use an predefined dpss object to reduce computation time when applied multiple times (default=NULL). If specified NW and k are obsolete.
#' @param detrend (default=T)
#'
#' @return spectra object
#' @export
CrossSpecMTM <- function(x,y=NULL,dt=1,NW=8,k=NULL,dpssIN=NULL,detrend=T) {

# NW      - number of windows (default = 8)
# how does NW in CrossSpecMTM fit to k,nw in SpecMTM?!


  # Checking Input
  if ( !is.null(dim(x)[2]) & is.null(y) ) {
    y <- x[,2]
    x <- x[,1]
  }
  # substract mean
  x <- x - mean(x)
  y <- y - mean(y)
  # Define parameters
  if (detrend) {
    x <- lm(x ~ seq(x))$residuals
    y <- lm(x ~ seq(x))$residuals
  }
  N <- length(x)

  if ( is.null(dpssIN) ) {

    if ( is.null(k) ) {
      k <- min(c(round(2*NW),N))
      k <- max(c(k-1,1))
    }
    s <- seq(0,1/dt-1/N,1/(N*dt))
    v <- 2*NW -1
    dpssIN <- dpss(n = N,k = k,nw = NW)
  } else {
    k = length(dpssIN$eigen)
    v=NULL
  }

  pls <- seq(2,(N+1)/2+1,1)
  if ( (N %% 2) == 1 ) pls <- pls[-length(pls)]


  E <- dpssIN$v
  V <- dpssIN$eigen
  fkx <- mvfft(E[,1:k]*matrix(rep(x,k),ncol=k))
  fky <- mvfft(E[,1:k]*matrix(rep(y,k),ncol=k))

  Pkx <- abs(fkx)^2
  Pky <- abs(fky)^2

  for ( i1 in 1:2 ) {
    if ( i1 == 1) {vari <- t(x)%*%x/N; Pk=Pkx}
    if ( i1 == 2) {vari <-  t(y)%*%y/N; Pk=Pky}
    P <- (Pk[,1]+Pk[,2])/2 # initial spectrum estimate
    Ptemp <-array(data=0,dim=c(N))
    P1 <- Ptemp
    tol  <- .0005*vari[1,1]/N          # usually within 'tol'erance in about three iterations, see equations from [2] (P&W pp 368-370).
    a <-vari[1,1]*abs(1-V)
    while ( sum(abs(P-P1)/N) > tol ) {
      b=(P%*%t(array(1,dim=c(k)))/(P%*%t(V)+array(1,dim=c(N))%*%t(a))) # weights
      wk=(b^2)*(array(1,dim=N)%*%t(V))           # new spectral estimate
      P1=colSums(t(wk)*t(Pk))/colSums(t(wk))
      Ptemp=P1
      P1=P
      P=Ptemp              # swap P and P1
    }
    if ( i1 == 1 ) {
    fkx=sqrt(k)*sqrt(wk)*fkx/matrix(rep(colSums(sqrt(t(wk))),k),ncol=k);
    Fx=P;  #Power density spectral estimate of x
    }
    if ( i1 == 2 ) {
    fky=sqrt(k)*sqrt(wk)*fky/matrix(rep(colSums(sqrt(t(wk))),k),ncol=k);
    Fy=P;  #Power density spectral estimate of y
    }
  }

  # Cross Spectrum
  Cxy <- rowSums(fkx*Conj(fky))
  resultFreqs = c(0,seq(1/N,1/2,1/N),seq(1/2-1/N,1/N,-1/N))
  # Individual Spectra
  # Fx
  # Fy
  # Phase spectrum
  ph  <- atan(Im(Cxy)/Re(Cxy))*180/pi*-1
  # Coherence
  c <- abs(Cxy)/sqrt(rowSums(apply(fkx,2,abs)^2)*rowSums(apply(fky,2,abs)^2))

  auxiliary <- list(dpss=dpssIN,
                    eigenCoefs=NULL,
                    eigenCoefWt=NULL,
                    taper="dpss",
                    v=v)

  spec.out <- list(crossspec = Cxy[pls],
                   spec = cbind(Fx[pls],Fy[pls]),
                   coh = c,
                   ph = ph,
                   freq=resultFreqs[pls],
                   series=cbind(x,y),
                   mtm=auxiliary)

  class(spec.out) = c("spec","mtm")
  return(spec.out)

} # EOF

