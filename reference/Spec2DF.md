# Transform Spec Object(s) Into a Dataframe

Transform Spec Object(s) Into a Dataframe

## Usage

``` r
Spec2DF(x)
```

## Arguments

- x:

  A spec object or (optionally named) list of spec objects.

## Value

A spec_df object

## See also

Other TidySpec:
[`DF2Spec()`](https://earthsystemdiagnostics.github.io/paleospec/reference/DF2Spec.md),
[`as.data.frame.spec()`](https://earthsystemdiagnostics.github.io/paleospec/reference/as.data.frame.spec.md),
[`as.spec()`](https://earthsystemdiagnostics.github.io/paleospec/reference/as.spec.md),
[`as_spec_df()`](https://earthsystemdiagnostics.github.io/paleospec/reference/as_spec_df.md)

## Author

Andrew Dolman \<andrew.dolman@awi.de\>

## Examples

``` r
library(PaleoSpec)
ts1 <- ts(SimPLS(1000, beta = 1))
sp1 <- SpecMTM(ts1)
sp1 <- AddConfInterval(sp1)
ts2 <- ts(rnorm(1000))
sp2 <- SpecMTM(ts2)
sp_lst <- list(sp1 = sp1, sp2 = sp2)
sp_df <- Spec2DF(sp_lst)
sp_df
#> # A tibble: 1,000 × 6
#>    spec_id  freq  spec   dof  lim.1  lim.2
#>    <chr>   <dbl> <dbl> <dbl>  <dbl>  <dbl>
#>  1 sp1     0.001 16.5   5.98  80.3   6.84 
#>  2 sp1     0.002 31.7   5.96 155.   13.1  
#>  3 sp1     0.003 21.2   5.99 103.    8.82 
#>  4 sp1     0.004 31.7   5.99 154.   13.2  
#>  5 sp1     0.005 28.9   5.98 140.   12.0  
#>  6 sp1     0.006 11.6   5.98  56.6   4.83 
#>  7 sp1     0.007  6.74  5.95  33.0   2.79 
#>  8 sp1     0.008  6.76  5.95  33.0   2.80 
#>  9 sp1     0.009  1.72  5.90   8.47  0.709
#> 10 sp1     0.01   2.98  5.96  14.6   1.23 
#> # ℹ 990 more rows
```
