# Transform a spec_df Object into a spec Object or List of spec Objects

Transform a spec_df Object into a spec Object or List of spec Objects

## Usage

``` r
DF2Spec(spec_df)
```

## Arguments

- spec_df:

  A spec_df object

## Value

A spec object or a list of spec objects

## See also

Other TidySpec:
[`Spec2DF()`](https://earthsystemdiagnostics.github.io/paleospec/reference/Spec2DF.md),
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
ts2 <- ts(rnorm(1000))
sp2 <- SpecMTM(ts2)
sp_lst <- list(sp1 = sp1, sp2 = sp2)
sp_df <- Spec2DF(sp_lst)
sp_df
#> # A tibble: 1,000 × 4
#>    spec_id  freq  spec   dof
#>    <chr>   <dbl> <dbl> <dbl>
#>  1 sp1     0.001 13.9   5.99
#>  2 sp1     0.002 20.3   5.99
#>  3 sp1     0.003 10.2   5.99
#>  4 sp1     0.004 11.5   5.99
#>  5 sp1     0.005  8.47  5.97
#>  6 sp1     0.006  2.67  5.95
#>  7 sp1     0.007  2.56  5.94
#>  8 sp1     0.008  6.25  5.96
#>  9 sp1     0.009  9.67  5.99
#> 10 sp1     0.01   8.85  5.98
#> # ℹ 990 more rows
str(DF2Spec(sp_df))
#> List of 2
#>  $ sp1:List of 4
#>   ..$ spec_id: chr [1:500] "sp1" "sp1" "sp1" "sp1" ...
#>   ..$ freq   : num [1:500] 0.001 0.002 0.003 0.004 0.005 0.006 0.007 0.008 0.009 0.01 ...
#>   ..$ spec   : num [1:500] 13.85 20.33 10.16 11.49 8.47 ...
#>   ..$ dof    : num [1:500] 5.99 5.99 5.99 5.99 5.97 ...
#>   ..- attr(*, "class")= chr [1:2] "spec" "list"
#>  $ sp2:List of 4
#>   ..$ spec_id: chr [1:500] "sp2" "sp2" "sp2" "sp2" ...
#>   ..$ freq   : num [1:500] 0.001 0.002 0.003 0.004 0.005 0.006 0.007 0.008 0.009 0.01 ...
#>   ..$ spec   : num [1:500] 0.694 1.142 1.687 0.872 1.309 ...
#>   ..$ dof    : num [1:500] 5.75 5.84 5.9 5.79 5.86 ...
#>   ..- attr(*, "class")= chr [1:2] "spec" "list"
#>  - attr(*, "class")= chr [1:2] "list" "spec_list"
```
