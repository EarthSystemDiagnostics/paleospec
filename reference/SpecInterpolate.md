# Interpolate spectrum

Interpolate a spectrum to a given frequency resolution. The
interpolation includes the spectral estimates themselves as well as
their degrees of freedom.

## Usage

``` r
SpecInterpolate(spec, freqRef, check = TRUE)
```

## Arguments

- spec:

  a spectral object of class `"spec"` or a list with the minimum
  components `freq`, `spec` and `dof` which are vectors of the same
  length giving the original frequency resolution and the corresponding
  spectral estimates and degrees of freedom.

- freqRef:

  numeric vector with the target frequency resolution.

- check:

  logical; if `TRUE` (the default) `spec` is checked for being a proper
  spectral object. This can be turned off in programmatic use when the
  respective check has already been performed.

## Value

a spectral object of class `"spec"` with the spectrum interpolated to
the target resolution.

## Author

Thomas Laepple

## Examples

``` r

freqRef <- seq(0.1, 0.5, 0.1)
spec <- list(freq = c(0.1, 0.2, 0.4, 0.5), spec = c(1, 2, 4, 5),
             dof = rep(1, 4))

SpecInterpolate(spec, freqRef)
#> $freq
#> [1] 0.1 0.2 0.3 0.4 0.5
#> 
#> $spec
#> [1] 1 2 3 4 5
#> 
#> $dof
#> [1] 1 1 1 1 1
#> 
#> attr(,"class")
#> [1] "spec"
```
