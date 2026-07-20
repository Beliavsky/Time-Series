# garma

[Back to the implemented-package index](../../README.md#implemented-packages).

## Overview

`garma.f90` translates the
[garma](https://cran.r-project.org/web/packages/garma/index.html) R package.

## Algorithms and Procedures

It implements the Gegenbauer ARMA numerical layer: Gegenbauer expansions, raw
periodograms and pole selection, Arteche GSP and LPR semiparametric estimates,
CSS, Whittle, and WLL objectives and fits, spectral inference, long-memory
removal, fitted values, Godet-style forecasts, and the Bartlett periodogram
diagnostic.

Discontinuous pole estimation uses the shared derivative-free Nelder-Mead
optimizer in `optimization_mod`.

The high-level regression workflow adds intercepts, drift, external regressors,
integer differencing and reintegration, future-regressor forecasts, CSS Hessian
inference, WLL exponent errors, and standard forecast-accuracy measures.

## Licensing

The translation is licensed under GPL-3.0-only.
