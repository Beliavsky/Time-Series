# robustarima

[Back to the implemented-package index](../../README.md#implemented-packages).

## Overview

`robustarima.f90` translates the
[robustarima](https://cran.r-project.org/web/packages/robustarima/index.html) R
package.

## Algorithms and Procedures

It provides filtered tau-estimation for regression with ARIMA errors and the
package's bounded redescending rho and psi functions, MAD-initialized M- and
tau-scales, clipped pseudo-observations for robust correlation diagnostics, and
a bounded-innovation ARMA filter.

Joint fitting uses partial-autocorrelation transforms for stationary AR and
invertible MA parameters, supports ordinary and seasonal differencing,
regression terms, and one seasonal MA factor, and can select an AR order by
robust AIC.

Forecasts, parameter covariance estimates, and iterative IO/AO/level-shift
scoring and cleaning are included.

Regression inference uses the package's filtered-design tau sandwich covariance
and reports its inverse-efficiency correction; the finite-difference Hessian
remains the fallback for a singular estimating equation and supplies inference
for ARMA parameters.

The filter propagates the innovations-form state and prediction covariance,
reports time-varying prediction scales, and implements the package's
rewind-and-accept rule for runs of extreme innovations that mark level shifts.

Outlier effect regressions use those prediction variances, and forecasts start
from the reconstructed robustly cleaned path.

## Licensing

The translation is licensed under BSD-3-Clause; see `LICENSE-ROBUSTARIMA`.
