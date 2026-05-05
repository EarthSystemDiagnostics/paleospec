# Variance estimate by integrating a part of the spectrum

Variance estimate by integrating a part of the spectrum

## Usage

``` r
GetVarFromSpectra(spec, f, dfreq = NULL, df.log = 0, bw = 3)
```

## Arguments

- spec:

  spectrum (list of spec,freq,dof) to be analysed

- f:

  f\[1\],f\[2\]: frequency interval to be analysed

- dfreq:

  frequency discretisation used in the temporary interpolation

- df.log:

  if \> 0, smooth the spectra prior to integrating

- bw:

  the bandwidth assumed for the confinterval calculation (from the
  multitaper spectral estimate)

## Value

list(var,dof) variance and corresponding dof

## Author

Thomas Laepple

## Examples

``` r
x<-ts(rnorm(100))
spec<-SpecMTM(x)
var(x) #Sample variance of the timeseries
#> [1] 1.004234
GetVarFromSpectra(spec,c(1/100,0.5))
#> $var
#> [1] 0.9748779
#> 
#> $dof
#> [1] 94.67882
#> 
GetVarFromSpectra(spec,c(0.25,0.5))
#> $var
#> [1] 0.4913773
#> 
#> $dof
#> [1] 48.16226
#> 
```
