# Weighted mean spectrum

Calculate the mean spectrum from a list of individual spectra including
weighting of the individual spectra and, if needed, interpolation to the
highest resolution frequency grid. The default weighting produces the
simple arithmetic mean. Interpolation is only performed when the
frequency axes of the individual spectra have different lengths and/or
differ in frequency discretization.

## Usage

``` r
MeanSpectrum(specList, iRemoveLowest = 1, weights = rep(1, length(specList)))
```

## Arguments

- specList:

  list of spectra, i.e. objects of class `"spec"` (see
  [`SpecMTM`](https://earthsystemdiagnostics.github.io/paleospec/reference/SpecMTM.md)
  for details), where each spectrum has to be a list of the vectors
  `freq`, `spec` and `dof` of a common length.

- iRemoveLowest:

  integer; number of lowest frequencies to remove from each individual
  spectral estimate (e.g. to remove detrending bias) prior to the
  interpolation and averaging.

- weights:

  numeric vector of weights; its length must match the number of
  elements in `specList`. The default setting of identical weights
  yields the simple arithmetic average of the input spectra; change this
  parameter accordingly to apply non-arithmetic averaging (see example).

## Value

object of class `"spec"` with the weighted mean spectrum, amended by the
element `nRecord` which gives the number of records contributing to each
mean spectral estimate.

## Author

Thomas Laepple and Thomas Münch

## Examples

``` r

# Simple arithmetic average
f1 <- 1 : 5
f2 <- f1
s1 <- rep(1, length(f1))
s2 <- rep(3, length(f2))
dof1 <- rep(1, length(f1))
dof2 <- rep(1, length(f2))

spectra <- list(list(freq = f1, spec = s1, dof = dof1),
                list(freq = f2, spec = s2, dof = dof2))

MeanSpectrum(spectra, iRemoveLowest = 0)
#> $freq
#> [1] 1 2 3 4 5
#> 
#> $spec
#> [1] 2 2 2 2 2
#> 
#> $dof
#> [1] 2 2 2 2 2
#> 
#> $nRecord
#> [1] 2 2 2 2 2
#> 
#> $lim.1
#> [1] 78.99578 78.99578 78.99578 78.99578 78.99578
#> 
#> $lim.2
#> [1] 0.5421701 0.5421701 0.5421701 0.5421701 0.5421701
#> 
#> $pval
#> [1] 0.05
#> 
#> attr(,"class")
#> [1] "spec"

# Weighted mean with interpolation
f1 <- 1 : 5
f2 <- 3 : 8
s1 <- rep(1, length(f1))
s2 <- rep(3, length(f2))
dof1 <- rep(1, length(f1))
dof2 <- rep(1, length(f2))

spectra <- list(list(freq = f1, spec = s1, dof = dof1),
                list(freq = f2, spec = s2, dof = dof2))

MeanSpectrum(spectra, iRemoveLowest = 0, weights = c(1, 2))
#> $freq
#> [1] 1 2 3 4 5 6 7 8
#> 
#> $spec
#> [1] 1.000000 1.000000 2.333333 2.333333 2.333333 3.000000 3.000000 3.000000
#> 
#> $dof
#> [1] 1 1 2 2 2 1 1 1
#> 
#> $nRecord
#> [1] 1 1 2 2 2 1 1 1
#> 
#> $lim.1
#> [1] 1018.25827 1018.25827   92.16174   92.16174   92.16174 3054.77481 3054.77481
#> [8] 3054.77481
#> 
#> $lim.2
#> [1] 0.1990491 0.1990491 0.6325317 0.6325317 0.6325317 0.5971473 0.5971473
#> [8] 0.5971473
#> 
#> $pval
#> [1] 0.05
#> 
#> attr(,"class")
#> [1] "spec"
# with some detrending bias removal
MeanSpectrum(spectra, iRemoveLowest = 1, weights = c(1, 2))
#> $freq
#> [1] 2 3 4 5 6 7 8
#> 
#> $spec
#> [1] 1.000000 1.000000 2.333333 2.333333 3.000000 3.000000 3.000000
#> 
#> $dof
#> [1] 1 1 2 2 1 1 1
#> 
#> $nRecord
#> [1] 1 1 2 2 1 1 1
#> 
#> $lim.1
#> [1] 1018.25827 1018.25827   92.16174   92.16174 3054.77481 3054.77481 3054.77481
#> 
#> $lim.2
#> [1] 0.1990491 0.1990491 0.6325317 0.6325317 0.5971473 0.5971473 0.5971473
#> 
#> $pval
#> [1] 0.05
#> 
#> attr(,"class")
#> [1] "spec"
```
