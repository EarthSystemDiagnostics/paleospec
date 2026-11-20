# Cache the Julia Multitaper module import and the R-facing Julia wrapper across
# calls within a session
.julia_cache <- new.env(parent = emptyenv())

.CheckJuliaConnectoR <- function() {
  if (!requireNamespace("JuliaConnectoR", quietly = TRUE)) {
    stop(
      "The 'JuliaConnectoR' package is needed but is not installed.\n",
      "Install it with:\n",
      "  install.packages(\"JuliaConnectoR\")\n\n",
      "Julia itself must also be installed:\n",
      "  https://julialang.org/downloads/\n\n",
      "Then install the Julia Multitaper package (once per Julia installation):\n",
      "  JuliaConnectoR::juliaEval('using Pkg; Pkg.add(\"Multitaper\")')",
      call. = FALSE
    )
  }
}

# Define (once per session) a Julia function that returns the output of
# Multitaper.mdmultispec() as a NamedTuple, which converts cleanly to a named R
# list via juliaGet(); a Dict does not preserve key order/pairing across the
# R/Julia boundary.
.get_md_spec_for_r <- function() {
  if (!exists("md_spec_for_r", envir = .julia_cache)) {
    .CheckJuliaConnectoR()
    JuliaConnectoR::juliaEval("using Multitaper")
    JuliaConnectoR::juliaEval('
function md_spec_for_r(tt, x; bw=5/length(tt), k=7, dt=1.0, nz=0)
    pkg, nu1 = Multitaper.mdmultispec(tt, x; bw=bw, k=k, dt=dt, nz=nz, dof=true, Ftest=false)
    return (freq = collect(pkg.f), spec = pkg.S, dof = nu1)
end
')
    .julia_cache$md_spec_for_r <- JuliaConnectoR::juliaFun("md_spec_for_r")
  }
  .julia_cache$md_spec_for_r
}


#' Multitaper Spectrum for Irregularly Sampled or Gappy Series (Chave 2019)
#'
#' Multitaper spectral estimate for irregularly sampled series, or regularly
#' sampled series with gaps, after Chave (2019), using a non-uniform FFT and
#' Slepian tapers computed for the actual sample times. The method is
#' implemented in the Julia package \code{Multitaper.jl} (Haley and Geoga 2020),
#' whose \code{mdmultispec()} function is called here via \code{JuliaConnectoR}.
#' Pass either a \code{ts} (optionally with \code{NA}s, for a gappy-regular
#' series) or a raw value vector plus explicit \code{times} (for an irregular
#' series).
#'
#' @param x A \code{ts} object (possibly with \code{NA}s) or, if \code{times}
#'   is given, a numeric vector of values.
#' @param times Optional numeric vector of sample times, the same length as
#'   \code{x}, for an irregular series. If \code{NULL}, the times and missing
#'   values are taken from the \code{ts}.
#' @param deltat The nominal sample spacing, used to scale the time axis
#'   (Nyquist is \code{1 / (2 * deltat)}). For an irregular series use e.g. the
#'   mean time step. Ignored for the \code{ts} method in the sense that the
#'   time axis is taken from the \code{ts} itself; see Details.
#' @param nw A positive double, the time-bandwidth product (default 4). The
#'   \code{mdmultispec} bandwidth is \code{bw = nw / N} with \code{N} the number
#'   of observations, matching Thomson's \eqn{NW = bandwidth \times N}
#'   convention.
#' @param k A positive integer, the number of Slepian tapers (default 7).
#' @param nz Zero-padding factor passed to \code{mdmultispec()}.
#' @param detrend logical; remove a linear trend (against \code{times}) before
#'   estimating the spectrum. Default \code{FALSE}: detrend the series
#'   beforehand if required.
#'
#' @return A \code{spec} object with \code{freq}, \code{spec} and \code{dof}.
#'   The zero frequency is dropped. Degrees of freedom are estimated adaptively
#'   by Julia and capped at \code{2 * k}.
#'
#' @details
#' \code{mdmultispec()} evaluates its non-uniform FFT over \eqn{(-\pi, \pi)}
#' with the times used as given, and its taper bandwidth is in cycles per time
#' unit, so the times must be in units of the nominal sample spacing (Nyquist =
#' 0.5 cycles per unit). The times are therefore rescaled to units of
#' \code{deltat}, the estimate is made with \code{dt = 1} and the frequencies
#' and spectrum are converted back to the original time units. For an irregular
#' series pass \code{deltat} = the nominal spacing (e.g. the mean time step).
#'
#' The degrees of freedom are capped because \code{k} tapers can give at most
#' \code{2 * k}; \code{mdmultispec}'s adaptive-weight estimate can greatly
#' exceed that where the spectrum is large relative to the series variance (up
#' to about 1e6 has been seen), which makes \code{\link{FilterSpec}} build
#' enormous prior filters and run very slowly.
#'
#' Requires Julia, the \code{JuliaConnectoR} R package and the Julia
#' \code{Multitaper} package:
#'
#' \enumerate{
#'   \item Install Julia from \url{https://julialang.org/downloads/}
#'   \item Install the \code{JuliaConnectoR} R package:
#'     \code{install.packages("JuliaConnectoR")}
#'   \item Install the Julia \code{Multitaper} package (once per Julia
#'     installation):
#'     \preformatted{JuliaConnectoR::juliaEval(
#'       'using Pkg; Pkg.add("Multitaper")'
#'     )}
#' }
#'
#' The Julia module is imported once per R session and cached.
#'
#' @section Credit:
#' The spectral estimator is the work of the authors of the Julia package
#' \code{Multitaper.jl}, Charlotte L. Haley and Christopher J. Geoga
#' (\url{https://github.com/lootie/Multitaper.jl}, GPL-2); this function is only
#' a wrapper that calls it. If you use it, please cite Haley and Geoga (2020)
#' as requested by the package, and Chave (2019) for the method.
#'
#' @references
#' Chave, A. D. (2019). Multitaper spectral analysis for irregularly sampled
#' and gappy data. \emph{Geophysical Journal International}, 218(3),
#' 2165-2178. \doi{10.1093/gji/ggz280}
#'
#' Haley, C. L. and Geoga, C. J. (2020). Multitaper.jl: A Julia package for
#' frequency domain analysis of time series. \emph{Journal of Open Source
#' Software}, 5(55), 2463. \doi{10.21105/joss.02463}
#'
#' @family functions to estimate power spectra
#' @seealso \code{\link{SpecMTM}}, \code{\link{SpecACF}},
#'   \code{\link{SpecMTMChaveMat}}
#' @author Andrew Dolman
#' @examples
#' \dontrun{
#' set.seed(42)
#' N <- 1000
#' x_full <- ts(SimPLS(N, beta = 1, alpha = 1))
#'
#' # Spectrum from the complete series
#' sp_full <- SpecMTMChave(x_full, nw = 2, k = 3)
#'
#' # Insert 20% random gaps
#' set.seed(7)
#' x_gap <- x_full
#' x_gap[sample(N, size = 0.2 * N)] <- NA
#' sp_gap <- SpecMTMChave(x_gap, nw = 2, k = 3)
#'
#' # An irregularly sampled series with explicit times
#' t_irr <- sort(runif(N, 0, N))
#' sp_irr <- SpecMTMChave(as.numeric(x_full), times = t_irr,
#'                        deltat = mean(diff(t_irr)), nw = 2, k = 3)
#'
#' gg_spec(list(`full series` = sp_full, `20% gaps` = sp_gap))
#' }
#' @export
SpecMTMChave <- function(x, times = NULL, deltat = 1, nw = 4, k = 7, nz = 0,
                         detrend = FALSE) {

  md_spec_for_r <- .get_md_spec_for_r()

  if (is.null(times)) {
    tt_full <- as.numeric(stats::time(x))
    xx_full <- as.numeric(x)
    keep <- !is.na(xx_full)
    tt <- tt_full[keep]
    xx <- xx_full[keep]
  } else {
    tt <- as.numeric(times)
    xx <- as.numeric(x)
  }

  if (detrend) xx <- as.numeric(stats::residuals(stats::lm(xx ~ tt)))

  bw <- nw / length(tt)

  # see Details: times in units of the nominal sample spacing
  tt <- (tt - tt[1]) / deltat

  out <- JuliaConnectoR::juliaGet(
    md_spec_for_r(tt, xx, bw = bw, k = as.integer(k), dt = 1, nz = nz))
  out$freq <- out$freq / deltat
  out$spec <- out$spec * deltat

  # Drop the zero (DC) frequency bin: FilterSpecLog() takes log(freq), and
  # SpecMTM/SpecACF outputs likewise start at the fundamental frequency, not 0.
  keep <- out$freq > 0
  # k tapers can give at most 2k degrees of freedom (see Details)
  dof <- pmin(out$dof[keep], 2 * k)
  as.spec(list(
    freq = out$freq[keep],
    spec = out$spec[keep],
    dof  = dof
  ))
}


#' Mean Multitaper Spectrum of Several Binned Records (Chave 2019)
#'
#' Multi-record counterpart to \code{\link{SpecACF}(m, ...)}: \code{m} is a
#' matrix of records already binned onto a \emph{shared} regular grid (one row
#' per bin, one column per record, \code{NA} = gap), the same representation
#' \code{SpecACF}/\code{BinTimeseries} produce. The spectrum of each record is
#' estimated with \code{\link{SpecMTMChave}} (i.e. \code{Multitaper.jl}'s
#' \code{mdmultispec()}, see there for credit and references) and the mean over
#' records is returned.
#'
#' @param m A matrix of binned records, one row per bin and one column per
#'   record, with \code{NA} for empty bins.
#' @param bin.width The bin width, the sampling interval of the rows of \code{m}.
#' @param nw,k Time-bandwidth product and number of tapers, as for
#'   \code{\link{SpecMTMChave}}.
#' @param min_pts Records (columns) with fewer populated bins than this are
#'   dropped.
#' @param detrend logical; remove a linear trend (against bin index) from each
#'   record first, as \code{SpecACF} does.
#'
#' @details
#' \code{mdmultispec()}'s frequency grid is sized from the number of
#' non-missing points it is actually given, so two records with different gap
#' fractions would otherwise come back on different frequency axes even though
#' they share a bin width and time origin. This is fixed here by choosing each
#' column's zero-padding factor \code{nz} so every column's FFT length pads out
#' to the same target regardless of how many of its bins are populated, after
#' which every record's spectrum lands on the same frequency axis and the mean
#' spectrum across records is a plain elementwise mean (no
#' \code{\link{SpecInterpolate}} needed). \code{nz} is solved exactly
#' (accounting for \code{mdmultispec}'s internal rounding up of the point count
#' to even) rather than approximately, since an approximate \code{nz} can round
#' to an off-by-one FFT length for some column lengths and not others.
#'
#' Real (non-simulated) record sets can have some records with very few
#' populated bins after binning, not enough for \code{k} well-concentrated
#' tapers at bandwidth \code{nw}. \code{k} is capped per column (both by the
#' standard \code{2 * nw} taper limit and by the column's own point count);
#' columns with fewer than \code{min_pts} points are dropped rather than
#' crashing.
#'
#' \code{rfreq} is the Nyquist frequency \code{1 / (2 * bin.width)}, matching
#' the \code{rfreq} convention that the Stan model
#' (\code{s_n = 0.5 * sigma_noise_adj^2 / rfreq}) was calibrated against: the
#' local \code{SpecACF()} this pipeline originally used (before it was replaced
#' by \code{PaleoSpec::SpecACF()} in March 2023) computed \code{rfreq} as
#' exactly the Nyquist frequency. \code{PaleoSpec::SpecACF()}'s current
#' \code{rfreq} is instead the full sampling frequency (about
#' \code{1 / bin.width}, computed as \code{max(freq)} before truncating to
#' positive frequencies), nearly twice the value this model expects, which
#' inflates fitted \code{sigma_noise} by about \code{sqrt(2)} to compensate.
#' \code{mdmultispec}'s own frequency axis already stops at Nyquist, so
#' \code{max(freq)} here is already the right quantity with no correction.
#'
#' @inheritSection SpecMTMChave Credit
#' @inherit SpecMTMChave references
#'
#' @return A \code{spec} object with \code{freq}, \code{spec} (mean over
#'   records), \code{dof} (summed over records), \code{bin.width}, \code{rfreq}
#'   and \code{nrec}.
#' @family functions to estimate power spectra
#' @seealso \code{\link{SpecMTMChave}}, \code{\link{SpecACF}}
#' @author Andrew Dolman
#' @export
SpecMTMChaveMat <- function(m, bin.width = 1, nw = 4, k = 7, min_pts = 4,
                            detrend = TRUE) {

  md_spec_for_r <- .get_md_spec_for_r()

  m <- as.matrix(m)
  n_bins <- nrow(m)
  target_nfft <- n_bins + (n_bins %% 2)

  # Matches PaleoSpec::SpecACF()'s own detrending: a linear fit against bin
  # index per column, residuals taken with NAs preserved in place.
  if (detrend) {
    i <- seq_len(n_bins)
    m <- apply(m, 2, function(y) stats::residuals(stats::lm(y ~ i, na.action = "na.exclude")))
  }

  specs <- lapply(seq_len(ncol(m)), function(j) {
    col <- m[, j]
    keep <- !is.na(col)
    tt <- as.numeric(which(keep))
    xx <- as.numeric(col[keep])
    if (length(xx) < min_pts) return(NULL)

    k_col <- max(2, min(k, floor(2 * nw) - 1, length(xx) %/% 2))

    base <- length(xx) + (length(xx) %% 2)
    nz <- target_nfft / base - 1
    bw <- nw / length(xx)

    JuliaConnectoR::juliaGet(md_spec_for_r(tt, xx, bw = bw, k = as.integer(k_col),
                                           dt = bin.width, nz = nz))
  })
  specs <- Filter(Negate(is.null), specs)

  freq_lengths <- vapply(specs, function(s) length(s$freq), integer(1))
  stopifnot(length(unique(freq_lengths)) == 1)

  spec_mat <- vapply(specs, function(s) s$spec, numeric(freq_lengths[1]))
  dof_mat  <- vapply(specs, function(s) s$dof,  numeric(freq_lengths[1]))

  freq <- specs[[1]]$freq
  keep <- freq > 0
  as.spec(list(
    freq = freq[keep],
    spec = rowMeans(spec_mat)[keep],
    dof  = rowSums(dof_mat)[keep],
    bin.width = bin.width,
    rfreq = max(freq),
    nrec = length(specs)
  ))
}
