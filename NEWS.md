# PaleoSpec 0.2.6

* Return vector of DOF from LogSmooth even if they are all the same
* Add tests for freq trimming in LogSmooth
* Minor documentation fixes

# PaleoSpec 0.2.5

* New function SimProxySeries, sim.proxy.series

# PaleoSpec 0.2.4

* Bug fix to `Bandpass()` was calling lowpass instead of Lowpass

# PaleoSpec 0.2.3

* Overhaul of spectral averaging and interpolation functions

# PaleoSpec 0.2.2

* Bug fix to `ApplyFilter()` function following up on #3:
  
  A bug occurred upon matching the filtered values back on the index positions
  of the original input time series vector after removing any leading/trailing
  NA values, which happened only in the case of duplicate values in the input
  series. This is fixed now by an update to the index matching.

# PaleoSpec 0.2.1

* Filter function `ApplyFilter()` now handles NA values in the input time
  series: 
  - leading and trailing NA values are automatically stripped from the
    input so that they do not migrate into the filtered result upon applying any
    of the endpoint constraint methods, but they are added in again after
    applying the filter, so that the output result has the same length as the
    input;
 - internal NA values are still ignored per default, but can be removed by
   linear interpolation from the neighbouring values by toggling the respective
   new function parameter.
* Source code for the function was tidied, the function documentation updated,
  and unit code testing was added.

# PaleoSpec 0.2.0

* New function `SimPLS()` to simulate powerlaw timeseries with a prescribed
  alpha, i.e. alpha*f^(-beta).

# PaleoSpec 0.1.0

* Proper package version.
  
# PaleoSpec 0.0.0.9000

* Rough migration from bundle of functions based on `paleolibrary` to
  `PaleoSpec` package.
