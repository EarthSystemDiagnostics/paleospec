# Variance of a powerlaw process if integrated from until frequency f

Integral of PSD=f^(-beta) from f1=1/N to f2=f this equals the variance
of a lowpass filtered powerlaw process WARNING: The result is not
normalized

## Usage

``` r
PS.VarUntilF(f, beta, N)
```

## Arguments

- f:

  frequency until which to integrate

- beta:

  powerlaw slope

- N:

  length of the timeseries

## Value

non-normalized variance

## Author

Thomas Laepple

## Examples

``` r
beta <- 1
signal <- ts(SimPowerlaw(beta,100000))
spec <- SpecMTM(signal)
v1 <- GetVarFromSpectra(spec,f=c(1/length(signal),0.5))
v2 <- GetVarFromSpectra(spec,f=c(1/length(signal),0.01))
PS.VarUntilF(0.01,beta,length(signal))/PS.VarUntilF(0.5,beta,length(signal))
#> [1] 0.6384378
v2$var/v1$var
#> [1] 0.6110537
```
