# Remove leading and trailing rows of all NA

Remove leading and trailing rows of all NA

## Usage

``` r
TrimNA(m, trim = c("all", "any"))
```

## Arguments

- m:

  a numeric matrix, data.frame or vector

- trim:

  trim leading and trailing rows of "all" NA or containing "any" NA
  values

## Value

a numeric matrix, data.frame or vector

## Examples

``` r
m <- matrix(c(NA, NA, NA, 1, NA, NA, NA, 1, 1, NA, NA, NA, 1:9, NA,NA,NA, 10:12, NA, 1, NA, NA,NA,NA), ncol = 3, byrow = TRUE)
m
#>       [,1] [,2] [,3]
#>  [1,]   NA   NA   NA
#>  [2,]    1   NA   NA
#>  [3,]   NA    1    1
#>  [4,]   NA   NA   NA
#>  [5,]    1    2    3
#>  [6,]    4    5    6
#>  [7,]    7    8    9
#>  [8,]   NA   NA   NA
#>  [9,]   10   11   12
#> [10,]   NA    1   NA
#> [11,]   NA   NA   NA
TrimNA(m)
#>       [,1] [,2] [,3]
#>  [1,]    1   NA   NA
#>  [2,]   NA    1    1
#>  [3,]   NA   NA   NA
#>  [4,]    1    2    3
#>  [5,]    4    5    6
#>  [6,]    7    8    9
#>  [7,]   NA   NA   NA
#>  [8,]   10   11   12
#>  [9,]   NA    1   NA
TrimNA(m, trim = "any")
#>      [,1] [,2] [,3]
#> [1,]    1    2    3
#> [2,]    4    5    6
#> [3,]    7    8    9
#> [4,]   NA   NA   NA
#> [5,]   10   11   12
```
