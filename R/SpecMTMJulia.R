# Cache the Julia Multitaper module import across calls within a session
.julia_cache <- new.env(parent = emptyenv())

.get_julia_multitaper <- function() {
  if (!exists("mtj", envir = .julia_cache)) {
    .julia_cache$mtj <- JuliaConnectoR::juliaImport("Multitaper")
  }
  .julia_cache$mtj
}

#' MTM Spectral Estimator via Julia
#'
#' Calls \code{mdmultispec} from the Julia \code{Multitaper} package via
#' \code{JuliaConnectoR}. Returns the same structure as \code{\link{SpecMTM}}.
#' Unlike \code{\link{SpecMTM}}, this function handles time series with missing
#' values (gaps), which are passed to Julia as a non-uniform time axis so that
#' the non-uniform FFT can be used for spectral estimation.
#'
#' @inheritParams SpecMTM
#' @param timeSeries A time series of equally spaced observations, possibly
#'   containing \code{NA} values representing gaps. Can be a \code{ts} object
#'   (in which case \code{deltat} is extracted automatically) or a plain numeric
#'   vector (in which case a sampling interval of 1 is assumed).
#' @param nw a positive double, the time-bandwidth product (default 2).
#'   Converted to Julia's \code{bw} parameter as \code{nw / length(timeSeries)}.
#' @param k a positive integer, the number of Slepian tapers (default 3).
#' @param detrend logical; remove the mean and linear trend before estimating
#'   the spectrum (default \code{TRUE}).
#'
#' @return A \code{spec} object: a list with at minimum \code{freq}, \code{spec},
#'   and \code{dof} vectors of equal length, plus \code{dt} and \code{n}.
#'   Degrees of freedom are estimated adaptively by Julia and will be lower at
#'   frequencies where gaps reduce the effective number of observations.
#'
#' @details
#' Requires Julia to be installed and the Julia \code{Multitaper} package to be
#' available. On first use, follow these steps:
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
#' The Julia module is imported once per R session and cached, so subsequent
#' calls within the same session have minimal overhead.
#'
#' @family functions to estimate power spectra
#' @seealso \code{\link{SpecMTM}}, \code{\link{SpecACF}}
#' @author Andrew Dolman
#' @examples
#' \dontrun{
#' set.seed(42)
#' N <- 1000
#' x_full <- ts(SimPLS(N, beta = 1, alpha = 1))
#'
#' # Spectrum from the complete series
#' sp_full <- SpecMTMJulia(x_full)
#'
#' # Insert 20% random gaps — SpecMTMJulia handles these via non-uniform FFT
#' set.seed(7)
#' x_gap <- x_full
#' x_gap[sample(N, size = 0.2 * N)] <- NA
#' sp_gap <- SpecMTMJulia(x_gap)
#'
#' gg_spec(list(`full series` = sp_full, `20% gaps` = sp_gap))
#' }
#' @export
SpecMTMJulia <- function(timeSeries, nw = 2, k = 3, detrend = TRUE) {

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

  dt <- if (is.ts(timeSeries)) deltat(timeSeries) else 1
  x  <- as.numeric(timeSeries)
  N  <- length(x)

  if (detrend)
    x <- residuals(lm(x ~ seq_along(x), na.action = na.exclude))

  bw     <- nw / N
  t_axis <- as.numeric(seq_len(N))
  obs    <- !is.na(x)

  mtj <- .get_julia_multitaper()
  # dof=TRUE makes mdmultispec return a 2-tuple: (MTSpectrum, Vector{Float64})
  raw <- JuliaConnectoR::juliaGet(
    mtj$mdmultispec(t_axis[obs], x[obs], dt = dt, bw = bw, k = as.integer(k), dof = TRUE)
  )

  freq <- JuliaConnectoR::juliaCall("collect", raw[[1]]$f)
  spec <- as.numeric(raw[[1]]$S)
  dof  <- as.numeric(raw[[2]])

  keep <- freq > 0
  freq <- freq[keep]
  spec <- spec[keep]
  dof  <- dof[keep]

  result <- list(
    freq = freq,
    spec = spec,
    dof  = dof,
    dt   = dt,
    n    = N
  )

  class(result) <- c("SpecMTMJulia", "spec")
  return(result)
}
