# Numerical correlation of random timeseries with different filtering (running mean)

Numerical correlation of random timeseries with different filtering
(running mean)

## Usage

``` r
PSP.CorAfterRollmean(N, betaSignal, betaNoise, R)
```

## Arguments

- N:

  Number of points of the timeseries

- betaSignal:

  powerlaw slope of the signal

- betaNoise:

  powerlaw slope of the noise

- R:

  expected correlation of the whole timeseries

## Value

correlation of unfiltered, 10point mean and 50 point mean

## Author

Thomas Laepple

## Examples

``` r
temp <- replicate(1000,PSP.CorAfterRollmean(1000,1,0,0.5))
rowMeans(temp)
#> [1] 0.4991412 0.8204650 0.9296877
PSP.CorUntilF(c(0.5,0.5/10,0.5/50),1,0,1000,0.5)
#> [1] 0.5000000 0.8253528 0.9341159
```
