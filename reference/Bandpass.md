# Calculate Weights for a Bandpass Filter

Derive the (smoothed) least square bandpass based on Bloomfield 1976

## Usage

``` r
Bandpass(omega.upper, omega.lower, n, sample = 1, convergence = T)
```

## Arguments

- omega.upper:

  upper cutoff frequency

- omega.lower:

  lower cutoff frequency

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
