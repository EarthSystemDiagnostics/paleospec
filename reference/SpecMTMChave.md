# Multitaper Spectrum for Irregularly Sampled or Gappy Series (Chave 2019)

Multitaper spectral estimate for irregularly sampled series, or
regularly sampled series with gaps, after Chave (2019), using a
non-uniform FFT and Slepian tapers computed for the actual sample times.
The method is implemented in the Julia package `Multitaper.jl` (Haley
and Geoga 2020), whose `mdmultispec()` function is called here via
`JuliaConnectoR`. Pass either a `ts` (optionally with `NA`s, for a
gappy-regular series) or a raw value vector plus explicit `times` (for
an irregular series).

## Usage

``` r
SpecMTMChave(
  x,
  times = NULL,
  deltat = 1,
  nw = 4,
  k = 7,
  nz = 0,
  detrend = FALSE
)
```

## Arguments

- x:

  A `ts` object (possibly with `NA`s) or, if `times` is given, a numeric
  vector of values.

- times:

  Optional numeric vector of sample times, the same length as `x`, for
  an irregular series. If `NULL`, the times and missing values are taken
  from the `ts`.

- deltat:

  The nominal sample spacing, used to scale the time axis (Nyquist is
  `1 / (2 * deltat)`). For an irregular series use e.g. the mean time
  step. Ignored for the `ts` method in the sense that the time axis is
  taken from the `ts` itself; see Details.

- nw:

  A positive double, the time-bandwidth product (default 4). The
  `mdmultispec` bandwidth is `bw = nw / N` with `N` the number of
  observations, matching Thomson's \\NW = bandwidth \times N\\
  convention.

- k:

  A positive integer, the number of Slepian tapers (default 7).

- nz:

  Zero-padding factor passed to `mdmultispec()`.

- detrend:

  logical; remove a linear trend (against `times`) before estimating the
  spectrum. Default `FALSE`: detrend the series beforehand if required.

## Value

A `spec` object with `freq`, `spec` and `dof`. The zero frequency is
dropped. Degrees of freedom are estimated adaptively by Julia and capped
at `2 * k`.

## Details

`mdmultispec()` evaluates its non-uniform FFT over \\(-\pi, \pi)\\ with
the times used as given, and its taper bandwidth is in cycles per time
unit, so the times must be in units of the nominal sample spacing
(Nyquist = 0.5 cycles per unit). The times are therefore rescaled to
units of `deltat`, the estimate is made with `dt = 1` and the
frequencies and spectrum are converted back to the original time units.
For an irregular series pass `deltat` = the nominal spacing (e.g. the
mean time step).

The degrees of freedom are capped because `k` tapers can give at most
`2 * k`; `mdmultispec`'s adaptive-weight estimate can greatly exceed
that where the spectrum is large relative to the series variance (up to
about 1e6 has been seen), which makes
[`FilterSpec`](https://earthsystemdiagnostics.github.io/paleospec/reference/FilterSpec.md)
build enormous prior filters and run very slowly.

Requires Julia, the `JuliaConnectoR` R package and the Julia
`Multitaper` package:

1.  Install Julia from <https://julialang.org/downloads/>

2.  Install the `JuliaConnectoR` R package:
    `install.packages("JuliaConnectoR")`

3.  Install the Julia `Multitaper` package (once per Julia
    installation):

        JuliaConnectoR::juliaEval(
              'using Pkg; Pkg.add("Multitaper")'
            )

The Julia module is imported once per R session and cached.

## Credit

The spectral estimator is the work of the authors of the Julia package
`Multitaper.jl`, Charlotte L. Haley and Christopher J. Geoga
(<https://github.com/lootie/Multitaper.jl>, GPL-2); this function is
only a wrapper that calls it. If you use it, please cite Haley and Geoga
(2020) as requested by the package, and Chave (2019) for the method.

## References

Chave, A. D. (2019). Multitaper spectral analysis for irregularly
sampled and gappy data. *Geophysical Journal International*, 218(3),
2165-2178. [doi:10.1093/gji/ggz280](https://doi.org/10.1093/gji/ggz280)

Haley, C. L. and Geoga, C. J. (2020). Multitaper.jl: A Julia package for
frequency domain analysis of time series. *Journal of Open Source
Software*, 5(55), 2463.
[doi:10.21105/joss.02463](https://doi.org/10.21105/joss.02463)

## See also

[`SpecMTM`](https://earthsystemdiagnostics.github.io/paleospec/reference/SpecMTM.md),
[`SpecACF`](https://earthsystemdiagnostics.github.io/paleospec/reference/SpecACF.md),
[`SpecMTMChaveMat`](https://earthsystemdiagnostics.github.io/paleospec/reference/SpecMTMChaveMat.md)

Other functions to estimate power spectra:
[`SpecACF()`](https://earthsystemdiagnostics.github.io/paleospec/reference/SpecACF.md),
[`SpecMTM()`](https://earthsystemdiagnostics.github.io/paleospec/reference/SpecMTM.md),
[`SpecMTMChaveMat()`](https://earthsystemdiagnostics.github.io/paleospec/reference/SpecMTMChaveMat.md)

## Author

Andrew Dolman

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(42)
N <- 1000
x_full <- ts(SimPLS(N, beta = 1, alpha = 1))

# Spectrum from the complete series
sp_full <- SpecMTMChave(x_full, nw = 2, k = 3)

# Insert 20% random gaps
set.seed(7)
x_gap <- x_full
x_gap[sample(N, size = 0.2 * N)] <- NA
sp_gap <- SpecMTMChave(x_gap, nw = 2, k = 3)

# An irregularly sampled series with explicit times
t_irr <- sort(runif(N, 0, N))
sp_irr <- SpecMTMChave(as.numeric(x_full), times = t_irr,
                       deltat = mean(diff(t_irr)), nw = 2, k = 3)

gg_spec(list(`full series` = sp_full, `20% gaps` = sp_gap))
} # }
```
