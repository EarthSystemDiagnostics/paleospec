# calculate weights for a bandpass filter

Derive the (smoothed) least square highpass based on Bloomfield 1976

## Usage

``` r
Highpass(omega.c, n = 9, sample = 1, convergence = T)
```

## Arguments

- omega.c:

  cutoff frequency

- n:

  length of the filter, has to be odd

- sample:

  sampling rate of the timeseries on which the filter will be applied
  (1/deltat)

- convergence:

  TRUE: smoothed least square lowpass; FALSE = unsmoothed

## Value

vector of filter weights

## Author

Thomas Laepple
