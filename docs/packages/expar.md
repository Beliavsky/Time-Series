# EXPAR

[Back to the implemented-package index](../../README.md#implemented-packages).

## Overview

`expar.f90` translates the
[EXPAR](https://cran.r-project.org/web/packages/EXPAR/index.html) R package.

## Algorithms and Procedures

It evaluates amplitude-dependent exponential autoregressions, estimates their
phi, pi, and exponential-scale parameters by conditional RSS minimization with
finite-difference BFGS, constructs AR-based starting values, selects the order
by AIC, corrected AIC, or BIC, and computes recursive point forecasts.

The principal estimation interfaces are also available through `arma_mod`, and
forecasting is available through `forecasting_mod`.

## Licensing

The translation is licensed under GPL-3.0-only; see `LICENSE-EXPAR`.
