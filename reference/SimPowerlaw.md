# Simulate a random timeseries with a powerlaw spectrum

Simulate a random timeseries with a powerlaw spectrum

## Usage

``` r
SimPowerlaw(beta, N)
```

## Arguments

- beta:

  slope

- N:

  length of timeseries to be generated

## Value

vector containing the timeseries

## Details

Method: FFT white noise, rescale, FFT back, the result is scaled to
variance 1

## See also

Other functions to generate timeseries with powerlaw like spectra:
[`SimFromEmpiricalSpec()`](https://earthsystemdiagnostics.github.io/paleospec/reference/SimFromEmpiricalSpec.md),
[`SimPLS()`](https://earthsystemdiagnostics.github.io/paleospec/reference/SimPLS.md),
[`SimProxySeries()`](https://earthsystemdiagnostics.github.io/paleospec/reference/SimProxySeries.md)

## Author

Thomas Laepple
