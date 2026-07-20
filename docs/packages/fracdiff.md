# fracdiff

[Back to the implemented-package index](../../README.md#implemented-packages).

## Overview

`fracdiff.f90` translates the
[fracdiff](https://cran.r-project.org/web/packages/fracdiff/index.html) R
package.

## Algorithms and Procedures

It provides Jensen-Nielsen FFT fractional differencing, GPH and tapered Sperio
semiparametric memory estimators, the Haslett-Raftery truncated likelihood
filter and stable BFGS fit with Hessian inference, and fracdiff-compatible
simulation.

Simulation reuses `arfima_mod`'s Durbin-Levinson fractional-noise generator,
while the FFT is centralized in `fourier_mod`.

Generic ARFIMA covariance, forecasting, and diagnostics are not duplicated.

## Licensing

The translation is licensed under GPL-2.0-or-later.
