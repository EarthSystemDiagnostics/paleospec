# Mean Multitaper Spectrum of Several Binned Records (Chave 2019)

Multi-record counterpart to
[`SpecACF`](https://earthsystemdiagnostics.github.io/paleospec/reference/SpecACF.md)`(m, ...)`:
`m` is a matrix of records already binned onto a *shared* regular grid
(one row per bin, one column per record, `NA` = gap), the same
representation `SpecACF`/`BinTimeseries` produce. The spectrum of each
record is estimated with
[`SpecMTMChave`](https://earthsystemdiagnostics.github.io/paleospec/reference/SpecMTMChave.md)
(i.e. `Multitaper.jl`'s `mdmultispec()`, see there for credit and
references) and the mean over records is returned.

## Usage

``` r
SpecMTMChaveMat(m, bin.width = 1, nw = 4, k = 7, min_pts = 4, detrend = TRUE)
```

## Arguments

- m:

  A matrix of binned records, one row per bin and one column per record,
  with `NA` for empty bins.

- bin.width:

  The bin width, the sampling interval of the rows of `m`.

- nw, k:

  Time-bandwidth product and number of tapers, as for
  [`SpecMTMChave`](https://earthsystemdiagnostics.github.io/paleospec/reference/SpecMTMChave.md).

- min_pts:

  Records (columns) with fewer populated bins than this are dropped.

- detrend:

  logical; remove a linear trend (against bin index) from each record
  first, as `SpecACF` does.

## Value

A `spec` object with `freq`, `spec` (mean over records), `dof` (summed
over records), `bin.width`, `rfreq` and `nrec`.

## Details

`mdmultispec()`'s frequency grid is sized from the number of non-missing
points it is actually given, so two records with different gap fractions
would otherwise come back on different frequency axes even though they
share a bin width and time origin. This is fixed here by choosing each
column's zero-padding factor `nz` so every column's FFT length pads out
to the same target regardless of how many of its bins are populated,
after which every record's spectrum lands on the same frequency axis and
the mean spectrum across records is a plain elementwise mean (no
[`SpecInterpolate`](https://earthsystemdiagnostics.github.io/paleospec/reference/SpecInterpolate.md)
needed). `nz` is solved exactly (accounting for `mdmultispec`'s internal
rounding up of the point count to even) rather than approximately, since
an approximate `nz` can round to an off-by-one FFT length for some
column lengths and not others.

Real (non-simulated) record sets can have some records with very few
populated bins after binning, not enough for `k` well-concentrated
tapers at bandwidth `nw`. `k` is capped per column (both by the standard
`2 * nw` taper limit and by the column's own point count); columns with
fewer than `min_pts` points are dropped rather than crashing.

`rfreq` is the Nyquist frequency `1 / (2 * bin.width)`, matching the
`rfreq` convention that the Stan model
(`s_n = 0.5 * sigma_noise_adj^2 / rfreq`) was calibrated against: the
local
[`SpecACF()`](https://earthsystemdiagnostics.github.io/paleospec/reference/SpecACF.md)
this pipeline originally used (before it was replaced by
[`PaleoSpec::SpecACF()`](https://earthsystemdiagnostics.github.io/paleospec/reference/SpecACF.md)
in March 2023) computed `rfreq` as exactly the Nyquist frequency.
[`PaleoSpec::SpecACF()`](https://earthsystemdiagnostics.github.io/paleospec/reference/SpecACF.md)'s
current `rfreq` is instead the full sampling frequency (about
`1 / bin.width`, computed as `max(freq)` before truncating to positive
frequencies), nearly twice the value this model expects, which inflates
fitted `sigma_noise` by about `sqrt(2)` to compensate. `mdmultispec`'s
own frequency axis already stops at Nyquist, so `max(freq)` here is
already the right quantity with no correction.

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

[`SpecMTMChave`](https://earthsystemdiagnostics.github.io/paleospec/reference/SpecMTMChave.md),
[`SpecACF`](https://earthsystemdiagnostics.github.io/paleospec/reference/SpecACF.md)

Other functions to estimate power spectra:
[`SpecACF()`](https://earthsystemdiagnostics.github.io/paleospec/reference/SpecACF.md),
[`SpecMTM()`](https://earthsystemdiagnostics.github.io/paleospec/reference/SpecMTM.md),
[`SpecMTMChave()`](https://earthsystemdiagnostics.github.io/paleospec/reference/SpecMTMChave.md)

## Author

Andrew Dolman
