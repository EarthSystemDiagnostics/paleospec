# Create a pair of random signals with powerlaw signal and powerlaw noise

The timeseries have an expected variance of 1 and an expected
correlation of r

## Usage

``` r
SimulatePowerlawSignalPair(n, beta.signal, beta.noise, r)
```

## Arguments

- r:

  expected correlation between both vectors

- N:

  Number of points per timeseries

- betaSignal:

  powerlaw slope of the signal

- betaNoise:

  powerlaw slope of the noise

## Value

list containing both vectors y1 and y2

## Author

Thomas Laepple

## Examples

``` r
mean(replicate(1000,{test <- SimulatePowerlawSignalPair(200,1,1,0.5);cor(test$y1,test$y2)}))
#> [1] 0.4931165
```
