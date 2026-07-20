# EXPARMA

[Back to the implemented-package index](../../README.md#implemented-packages).

## Overview

`exparma.f90` translates the
[EXPARMA](https://cran.r-project.org/web/packages/EXPARMA/index.html) R
package.

## Algorithms and Procedures

It evaluates exponential autoregressive moving-average recursions in which AR
coefficients depend on the previous observation amplitude and MA coefficients
depend on the previous residual amplitude.

Conditional RSS estimation uses the shared finite-difference BFGS optimizer,
starting values reuse Hannan-Rissanen ARMA estimates, and a two-dimensional AIC
search selects AR and MA orders.

The principal interfaces are also available through `arma_mod`.

## Licensing

The translation is licensed under GPL-3.0-only; see `LICENSE-EXPARMA`.
