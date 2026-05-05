# Simulate a timeseries with length N which has a spectra consisting of two powerlaws

Simulate a timeseries with length N which has a spectra consisting of
two powerlaws

## Usage

``` r
SimPowerlawPiecewise(beta1, beta2, N, deltat = 1, breakpoint = 1/50)
```

## Arguments

- beta1:

  slope for frequencies lower than breakpoint

- beta2:

  slope for frequencies higher than breakpoint

- N:

  Number of points to simulate

- deltat:

  timestep of the timeseries

- breakpoint:

  frequency of the breakpoint

## Value

N random numbers drawn according to the piecewise powerlaw PSD

## Author

Thomas Laepple
