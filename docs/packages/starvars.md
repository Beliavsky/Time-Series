# starvars

[Back to the implemented-package index](../../README.md#implemented-packages).

## Overview

`starvars.f90` translates the
[starvars](https://cran.r-project.org/web/packages/starvars/index.html) R
package.

## Algorithms and Procedures

It estimates vector logistic smooth-transition autoregressions by profile
nonlinear least squares or concentrated Gaussian likelihood, including
grid-selected transition starts, conditional coefficient inference, and AIC and
BIC.

Forecasting supports recursive point predictions, Gaussian Monte Carlo paths,
and residual bootstrap paths with central intervals.

The module also implements the package's joint multivariate linearity test,
Bartlett long-run variance estimator, recursive common covariance-break CUMSUM
test, and grouped realized covariance construction from prices or returns.

Aggregation periods can be supplied directly as integer group identifiers or
derived from `date_t` values at daily, monthly, quarterly, or yearly frequency.

Principal modeling interfaces are available through `multivariate_mod`,
forecasts through `forecasting_mod`, tests through `diagnostics_mod`, and
realized covariance through `volatility_mod`.

## Licensing

The translation is licensed under GPL-3.0-or-later; see `LICENSE-STARVARS`.
