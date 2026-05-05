#' Simulate a Proxy Time Series Assuming a Power-Law Power Spectrum of the
#' Climate
#'
#' @description This function creates a time series that is designed such that
#'   it can be interpreted as a typical paleo proxy record. It corresponds to
#'   what one may expect when an underlying signal, exhibiting power-law
#'   frequency scaling over a finite range of frequencies (from f_low to
#'   f_high), is subject to archive and laboratory smoothing and then subsampled
#'   at potentially unequally spaced sampling times, with a single noisy
#'   measurement made on a finite number of signal carriers retrieved from the
#'   same physical sample.
#'
#' @param a If a >= 0: proportionality constant of the underlying expected power
#'   spectral density P = a * f^-b.
#'
#'   If a < 0: expected mean square (variance) of the final time series (the
#'   value of the function). Thus, the underlying power spectral density is
#'   given by P ~ f^-b.
#' @param b Scaling exponent of the underlying power spectral density P = a *
#'   f^-b.
#' @param nt Determines the lower cut-off frequency: f_low = f.scl / nt.
#'
#'   If t.smpl == NULL: number of sampling times (length of the final time
#'   series).
#' @param f.scl If t.smpl != NULL: Scales the cut-off frequencies (lower and
#'   upper): f_low = f.scl / nt, and f_high = f.scl * 0.5.
#'
#'   If t.smpl == NULL: f.scl has no effect.
#' @param smth.arch Sets type and time scale of archive smoothing. Possible
#'   values of type are "n": no archive smoothing; "bioturbation": 'Berger and
#'   Heath' impulse response function, with a Lorentzian-shaped squared spectral
#'   transfer function 1 / (1 + (2*pi*f*tau)^2); "diffusion": Gaussian-shaped
#'   squared spectral transfer function exp(-(2*pi*f*tau)^2).
#'
#'   Possible values of tau are positive real numbers. Setting tau == 0 is
#'   equivalent to setting type="n".
#' @param smth.lab Sets type and time scale of laboratory smoothing. Possible
#'   values of type are "n": no archive smoothing "rect": a rectangular window
#'   in the time domain (archive slice) -> Squared-sinc-shaped squared spectral
#'   transfer function sinc^2(pi*f*tau)
#'
#'   Possible values of tau are positive real numbers. Setting tau == 0 is
#'   equivalent to setting type="n".
#' @param N N > 0: Finite number of signal carriers used for a single
#'   measurement, implying an additional sample white noise scaled by 1/N. N ==
#'   0: Infinite number of signal carriers, implying no additional sample white
#'   noise.
#' @param t.smpl The desired (arbitrary) sampling times (may contain duplicated
#'   entries and does not need to be ordered in time). If t.smpl == NULL:
#'   sampling times are 0:(nt-1), using the FFT algorithm for fast computation.
#'   If t.smpl != NULL: specifies the desired arbitrary (and potentially
#'   unequally spaced) sampling times. When set to 0:(nt-1), the result is
#'   identical to that obtained with t.smpl == NULL, but computation is
#'   (potentially much) slower because FFT algorithm is not used.
#' @param var.noise Variance of additional sample white noise (which may
#'   represent, for example, measurement noise).
#' @param val If val == 0: the value of the function is the time series.
#'
#'   If val > 0: the value of the function is a list that holds information
#'   regarding the true PSD. It has the components
#'
#'   $fax: the frequency axis $psd: the power spectral density $var: the
#'   variance of the total sample white noise (component exists only if the
#'   sample white noise is not included in the component $psd, see below)
#'
#'   Possible positive values of val are
#'
#'   1: $psd holds the full PSD, including the sample white noise; vectors $psd
#'   and $fax are limited to those elements where $fax > 0
#'
#'   2: $psd holds the PSD excluding the sample white noise, the variance of
#'   which is provided in the component $var; vectors $psd and $fax are limited
#'   to those elements where $fax > 0
#'
#'   3: same as '1', but the vectors $psd and $fax extend over all frequencies
#'   (including zero and the negative ones)
#'
#'   4: same as '2', but the vectors $psd and $fax extend over all frequencies
#'   (including zero and the negative ones)
#'
#'   Note: If t.smpl != NULL and val == 1 or val == 3, then the function
#'   internally sets val <- val + 1, because t.smpl might specify unequally
#'   spaced sampling times such that the PSD of the sample white noise might not
#'   be defined.
#' @param PSD.format When val > 0; if PSD.format = "torb", the returned PSD is
#' in the format list(fax, psd); if PSD.format = "spec", a PSD in the standard
#' PaleoSpec format is returned.
#' @param rseeds vector of length 3, integer or NA. Random seeds to fix the
#' realisation of specific random processes. There are three random processes:
#'
#'   process-1 (controlled by rseeds[1]) represents the climate
#'   process-2 (controlled by rseeds[2]) represents the signal carrier mixing
#'   paths
#'   process-3 (controlled by rseeds[3]) represents the additional sample white
#'   noise (e.g., measurement noise)
#'
#'   If rseeds[i] == NA: a random sequence is drawn with rnorm() without setting
#'   a specific random seed before.
#'
#'   If rseeds[i] == integer: a random seed is set with set.seed(rseeds[i]) and
#'   then a random sequence is drawn with rnorm. After this draw the random
#'   number generator is reset to its previous state in order to avoid any
#'   interference with other calls to the random number generator inside or
#'   outside of sim.proxy.series.
#' @param nser integer (>=1) If val == 0: number of time series generated (all
#'   with the same parameters). If val > 0: nser has no effect.
#' @return A vector or list.
#' @family functions to generate timeseries with powerlaw like spectra
#' @author Torben Kunz
#' @export
#'
#' @examples
#'
#' ## Generate a power-law timeseries with no proxy archiving or measurement effects
#'
#' n = 1e04
#'
#' set.seed(1)
#' ts1 <- SimProxySeries(a = 0.1, b = 1, nt = n)
#' ts1 <- ts(ts1, deltat = 1)
#' plot(ts1)
#'
#' sp1 <- SpecMTM(ts1)
#' LPlot(sp1, col = "grey")
#' abline(log10(0.1), -1, col = "black")
#'
#' ## Generate the same climate timeseries (by setting seed the same) but with smoothing
#' ## from bioturbation
#'
#' tau_b <- 50
#'
#' set.seed(1)
#' ts2 <- SimProxySeries(a = 0.1, b = 1, nt = n,
#'  smth.arch = list(type = "bioturbation", tau = tau_b))
#' ts2 <- ts(ts2, deltat = 1)
#' plot(ts1)
#' lines(ts2, col = "blue")
#'
#' sp2 <- SpecMTM(ts2)
#' LPlot(sp1, ylim = range(c(sp1$spec, sp2$spec)), col = "grey")
#' abline(log10(0.1), -1, col = "black")
#' LLines(sp2, col = "blue")
#'
#' ## Then add measurement error
#'
#' set.seed(1)
#' ts3 <- SimProxySeries(a = 0.1, b = 1, nt = n,
#' smth.arch = list(type = "bioturbation", tau = tau_b),
#'  var.noise = 0.1^2)
#' ts3 <- ts(ts3, deltat = 1)
#'
#' sp3 <- SpecMTM(ts3)
#' LLines(sp3, col = "pink")
#' abline(h = 1/2*(0.1^2 / diff(range(sp3$freq))), col = "red")
#'
#' ## sample time series at arbitrary times
#'
#' #' set.seed(1)
#' time <- runif(100, 0, n)
#' ts4 <- SimProxySeries(a = 0.1, b = 1, nt = n, t.smpl = time,
#' smth.arch = list(type = "bioturbation", tau = tau_b),
#'  var.noise = 0.1^2)
#'
#' dat <- data.frame(time = time, ts4 = ts4)
#' dat <- dat[order(time), ]
#' plot(ts4~time, data = dat, type = "b")
#'
#' # Note: even though seed was set, these points will not be from the same
#' # timeseries as ts3.
#'
#'
#' ## Alternatively, generate the whole timeseries and subset at the desired
#' ## timepoints if integer timepoints are acceptable. This can often be faster
#' ## as it allows the fast FFT algorithm to be used.
#'
#' time <- round(sort(time))
#' ts3_sub <- ts3[time]
#'
#' plot(ts1, col = "grey")
#' lines(ts2, col = "blue")
#' lines(time, ts3_sub, col = "Orange", pch = 16)
#' points(time, ts3_sub, col = "Orange", pch = 16)
#'
#' # Return a PaleoSpec standard spec object or Torben Kunz format
#'
#' expected_spec_torb <- SimProxySeries(a = 0.1, b = 1, nt = n, t.smpl = time,
#'                                     smth.arch = list(type = "bioturbation", tau = tau_b),
#'                                     var.noise = 0.1^2,
#'                                     val = 1)
#'
#' expected_spec_std <- SimProxySeries(a = 0.1, b = 1, nt = n, t.smpl = time,
#'                                      smth.arch = list(type = "bioturbation", tau = tau_b),
#'                                     var.noise = 0.1^2,
#'                                     val = 1,
#'                                     PSD.format = "spec")
#'
#' gg_spec(expected_spec_std)
#' plot(expected_spec_torb$fax, expected_spec_torb$psd, log = "xy", type = "l")
#'
SimProxySeries <- function(a = -1, b = 1,
                           nt = 100, f.scl = 1,
                           smth.arch = list(type = "n", tau = 0),
                           smth.lab = list(type = "n", tau = 0),
                           N = 0,
                           t.smpl = NULL,
                           var.noise = 0,
                           val = 0,
                           PSD.format = c("torb", "spec"),
                           rseeds = c(NA, NA, NA),
                           nser = 1) {


  # validate inputs
  stopifnot(length(val) == 1)
  stopifnot(val %in% 0:4)

  smth.arch$type <- match.arg(smth.arch$type, c("n", "bioturbation", "diffusion"))
  smth.lab$type <- match.arg(smth.lab$type, c("n", "rect"))

  PSD.format <- match.arg(PSD.format)

  if (is.null(t.smpl)) {
    f.scl <- 1
  } else {
    if (val == 1 | val == 3) val <- val + 1
  }

  fax <- 0:(nt - 1) / nt
  fax[fax > 0.5] <- fax[fax > 0.5] - 1
  fax <- fax * f.scl

  P <- c(0, abs(a) * abs(fax[2:nt])^(-b))

  var.noise.N <- 0

  if (smth.arch$type != "n" | smth.lab$type != "n") {
    tf2 <- rep(1, length(fax))
    if (smth.arch$type == "bioturbation") tf2[-1] <- tf2[-1] / (1 + (2 * pi * fax[-1] * smth.arch$tau)^2)
    if (smth.arch$type == "diffusion") tf2[-1] <- tf2[-1] * exp(-(2 * pi * fax[-1] * smth.arch$tau)^2)
    if (smth.lab$type == "rect" & smth.lab$tau > 0) {
      y <- pi * fax[-1] * smth.lab$tau
      tf2[-1] <- tf2[-1] * (sin(y) / y)^2
    }
    if (N > 0) var.noise.N <- mean(P * (1 - tf2)) * f.scl / N
    P <- P * tf2
  }

  if (a < 0) {
    var.tot <- var.noise + var.noise.N + mean(P) * f.scl
    P <- abs(a) * P / var.tot
    var.noise <- abs(a) * var.noise / var.tot
    var.noise.N <- abs(a) * var.noise.N / var.tot
  }

  if (val == 0) {
    xfft <- mvfft(gen.ran.seq(nt, nser, rseeds[1])) * sqrt(P)
    if (is.null(t.smpl)) {
      x <- mvfft(xfft, inverse = TRUE) / nt
    } else {
      x <- matrix(0, length(t.smpl), nser)
      for (iif in 1:nt) {
        x <- x + exp(1i * 2 * pi * fax[iif] * t.smpl) %*% t(xfft[iif, ])
      }
      x <- x / nt
    }
    x <- Re(x) * sqrt(f.scl)
    if (var.noise > 0) x <- x + sqrt(var.noise) * gen.ran.seq(dim(x)[1], nser, rseeds[3])
    if (var.noise.N > 0) x <- x + sqrt(var.noise.N) * gen.ran.seq(dim(x)[1], nser, rseeds[2])
    if (nser == 1) x <- as.vector(x) else x <- x
  } else if (val == 1) {
    x <- list(fax = fax[fax > 0], psd = P[fax > 0] + var.noise + var.noise.N)
  } else if (val == 2) {
    x <- list(fax = fax[fax > 0], psd = P[fax > 0], var = var.noise + var.noise.N)
  } else if (val == 3) {
    x <- list(fax = fax, psd = P + var.noise + var.noise.N)
  } else if (val == 4) {
    x <- list(fax = fax, psd = P, var = var.noise + var.noise.N)
  }

  if (val > 0) {
    if(PSD.format == "spec") {
      x <- Torb2Spec(x)
    }
  }
    return(x)
  }

#' @rdname SimProxySeries
#' @export
sim.proxy.series <- SimProxySeries


#' Generate Random Numbers Without Affecting the State of .Random.seed
#'
#' @param n number of observations
#' @param nser number of replicate series
#' @param s NA or integer, passed to set.seed() if integer
#' @return a matrix
#' @keywords internal
gen.ran.seq <- function(n, nser, s) {
  if (is.na(s)) {
    x <- matrix(rnorm(n * nser), n, nser)
  } else {
    if (exists(".Random.seed")) rs <- .Random.seed
    set.seed(s)
    x <- matrix(rep(rnorm(n), nser), n, nser)
    if (exists("rs")){.Random.seed <<- rs}    else {
      set.seed(NULL)
    }
  }
  x
}


#' Convert Torben Kunz' spec output to standard PaleoSpec spec object
#'
#' @param x
#'
#' @return spec object
#' @keywords internal
Torb2Spec <- function(x){
  as.spec(list(freq = x$fax, spec = x$psd))
}


