# lowpass filtered expected correlation of powerlaw signal pair

Correlation of two timeseries with powerlaw signal and powerlaw noise
evaluated until f this equals the correlation of linearly coupled
lowpass filtered powerlaw process

## Usage

``` r
PSP.CorUntilF(f, betaSignal, betaNoise, N, r)
```

## Arguments

- f:

  frequency until which to integrate

- betaSignal:

  powerlaw slope of the signal

- betaNoise:

  powerlaw slope of the noise

- N:

  Number of points per timeseries

- r:

  expected correlation of the unfiltered

## Value

expected correlation of the timeseries

## Author

Thomas Laepple

## Examples

``` r
temp <- replicate(1000,PSP.CorAfterRollmean(1000,1,0,0.5))
rowMeans(temp)
#> [1] 0.4998948 0.8218859 0.9308612
PSP.CorUntilF(f=c(0.5,0.5/10,0.5/50),betaSignal=1,betaNoise=0,N=1000,r=0.5)
#> [1] 0.5000000 0.8253528 0.9341159
```
