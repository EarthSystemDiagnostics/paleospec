#' @title MTM coherence estimator
#' @description estimates coherence using slepian tapers with adaptive
#'   weighting. Function is an adaption of Peter Huybers Matlab coherence
#'   function cmtm.m
#' @param x,y input data vectors (or array if y=NULL)
#' @param dt time stepping (default = 1)
#' @param nw number of Slepian taper (default = 8)
#' @param k order (default = NULL)
#' @param dpssIN allows to use an predefined dpss object to reduce computation
#'   time when applied multiple times (default=NULL). If specified nw and k are
#'   obsolete.
#' @param detrend (default=T)
#' @param qbias bias correction fo coherence estimate (default = 0/no)
#' @param confn number of interactions to use in MC for uncertainty estimation
#'   of the phase (default = 0)
#'
#' @return spectra object
#' @export
#' @examples
Coherence <- function(x, y = NULL, dt = 1, nw = 8, k = NULL, detrend = T, qbias = 0, confn = 0) {
  if (!is.null(dim(x)[2]) & is.null(y)) {
    y <- x[, 2]
    x <- x[, 1]
  }

  coh <- PaleoSpec::CrossSpecMTM(x, y, dt = dt, nw = nw, k = k, detrend = detrend)

  N <- length(x)
  c <- coh$c
  ph <- coh$ph
  v <- coh$mtm$v
  s <- seq(0, 1 / dt - 1 / N, 1 / (N * dt))

  pls <- seq(2, (N + 1) / 2 + 1, 1)

  if (qbias > 0) {
    c <- CoherenceBias(dof = v, cb = c)
  }

  if (confn > 0) {
    cb <- t(CoherenceBias(dof = v, cb = c)) # doubled bias correction when qbias and confn is active?!....

    iter.steps <- confn
    phi <- array(data = NaN, dim = c(iter.steps, length(c)))
    ciph <- phi
    for (iter in 1:iter.steps) {
      fx <- fft(rnorm(length(x)) + 1)
      fx <- fx / sum(abs(fx))
      fy <- fft(rnorm(length(y)) + 1)
      fy <- fy / sum(abs(fy))
      ys <- Re(fft(fy * sqrt(1 - t(cb)^2), inverse = T))
      ys <- ys + Re(fft(fx * t(cb), inverse = T))
      xs <- Re(fft(fx, inverse = T))

      tmp <- PaleoSpec::CrossSpecMTM(x = xs, y = ys, dt = dt, nw = nw, k = k, detrend = detrend)
      phi[iter, ] <- tmp$ph
      ciph[iter, ] <- tmp$c
    }

    phi <- apply(phi, 2, function(x) {
      quantile(x, probs = c(0.025, 0.975))
    })
    phi <- (phi[2, ] + phi[1, ] * -1) / 2
    phi <- stats::filter(phi, filter = c(1, 1, 1) / 3, method = "convolution")
  }
  c <- c[pls]
  s <- s[pls]
  ph <- ph[pls]

  output <- list(
    coh = c, # return ,spec' instead of crossspec to allow further handling e.g. with LogSmooth...
    specs = coh$specs, # individual specs....
    ph = ph,
    freq = s[pls],
    series = cbind(x, y),
    mtm = coh$mtm
  )

  if (confn > 0) {
    output <- c(output, list(ph.conf = phi))
  }

  # Coherence confidence level
  if (all(confn != 0)) {
    conf.lev <- c(0.95)
    c.conf <-
      sapply(conf.lev, function(y) {
        sapply(c, function(x) {
          CoherenceConf(v, lev = y, unbias = 1, c = 0)
        })
      })
    output <- c(output, list(c.conf = c(c.conf)))
  }

  return(output)
}

CoherenceBias <- function(cb, dof) {
  if (dof < 2) stop("dof has to be >=2")
  if ((sum(cb < 0) + sum(cb > 1)) > 0) stop("Coherence must be between 0 and 1")
  if (dof > 50) {
    dof <- 50
    warning("Using 50 degrees of freedom")
  }

  ### The Estimation of the look-up-table seems to work. The results are equal
  ### to those of Huybers Matlab script. Nevertheless the provided lot looks
  ### different to the one estimated here...

  n <- seq(2, 50)
  c <- seq(0.1, 1, 0.1)
  expect <- array(data = NaN, dim = c(length(n), length(c)))
  z <- seq(0, 1, .1)
  for (i3 in 1:length(n)) {
    for (i2 in 1:(length(c) - 1)) {
      A <- M <- f <- vector()
      for (i1 in 1:length(z)) {
        A[1] <- 1
        # Calculated according to: Amos and Koopmans, "Tables of the distribution of the
        # coefficient of coherence for stationary bivariate Gaussian processes", Sandia
        # Corporation, 1963
        #
        # Also see the manuscript of Wunsch, C. "Time-Series Analysis.  A Heuristic Primer".
        for (k in 1:(n[i3] - 1)) {
          A[k + 1] <- A[k] * (n[i3] - k) * (2 * k - 1) / ((2 * n[i3] - (2 * k + 1)) * k) * ((1 - c[i2] * z[i1]) / (1 + c[i2] * z[i1]))^2
        }

        f[i1] <- 2 * (n[i3] - 1) * (1 - c[i2]^2)^n[i3] * z[i1] * (1 - z[i1]^2)^(n[i3] - 2) / ((1 + c[i2] * z[i1]) * (1 - c[i2] * z[i1])^(2 * n[i3] - 1)) *
          gamma(n[i3] - .5) / (sqrt(pi) * gamma(n[i3])) * sum(A)
      }

      for (i1 in 2:(length(f) / 2)) {
        M[i1] <- (f[2 * (i1 - 1) + 1] + 4 * f[2 * i1] + f[2 * i1 + 1]) * z[2 * i1]
      }

      expect[i3, i2] <- sum(c(M, 1), na.rm = T) / (6 * (length(M)))
    }
    expect[i3, i2 + 1] <- 1
  }


  ec <- vector()
  ### Interpolieren zum richtigen Freiheitsgrad
  for (i in 1:length(c)) ec[i] <- approx(n, expect[, i], dof)$y
  # Interpolieren zu cb
  result <- approx(ec, c, cb)$y
  result[is.na(result)] <- 0 ### NA coherencies = below the table values are set to zero
  return(result)
}


CoherenceConf <- function(v, lev = 0.95, unbias = 1, c = 0) {
  z <- seq(0, 1, .0005)

  f <- vector()
  for (i1 in 1:length(z)) {
    A <- vector()
    A[1] <- 1
    for (k in 1:(v - 1)) {
      A[k + 1] <- A[k] * (v - k) * (2 * k - 1) / ((2 * v - (2 * k + 1)) * k) * ((1 - c * z[i1]) / (1 + c * z[i1]))^2
    }
    f[i1] <- 2 * (v - 1) * (1 - c^2)^v * z[i1] * (1 - z[i1]^2)^(v - 2) / ((1 + c * z[i1]) * (1 - c * z[i1])^(2 * v - 1)) *
      gamma(v - .5) / (sqrt(pi) * gamma(v)) * sum(A)
  }

  F <- vector()
  for (i1 in seq(2, length(f) / 2)) F[i1] <- f[2 * (i1 - 1) + 1] + 4 * f[2 * i1] + f[2 * i1 + 1]

  F <- F / (6 * length(F))
  F <- 1 - cumsum(F[length(F):1] / sum(F, na.rm = T))
  F <- c(F[length(F):1], 1)
  F[is.na(F)] <- 0
  Fz <- c(z[seq(1, length(z) - 2, 2)], 1)
  pl <- which(diff(F) > 0)
  pl <- c(1, pl + 1)
  cl <- approx(F[pl], Fz[pl], lev)$y

  if (unbias) cl <- CoherenceBias(cb = cl, dof = v)

  return(cl)
}
