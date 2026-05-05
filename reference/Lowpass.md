# Calculate weights for lowpass filter

Derive the (smoothed) least square lowpass, given the cutoff frequency
omega.c and the length of the filter n

based on Bloomfield 1976

## Usage

``` r
Lowpass(omega.c, n = 9, sample = 1, convergence = TRUE)
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
