# echos

[Back to the implemented-package index](../../README.md#implemented-packages).

## Overview

`echos.f90` translates the
[echos](https://cran.r-project.org/web/packages/echos/index.html) R package.

## Algorithms and Procedures

It constructs exact-density random recurrent reservoirs, rescales them to a
requested general-matrix spectral radius, and generates leaky tanh echo states
from lagged inputs.

Only the linear readout is estimated: intercept-unpenalized ridge candidates
are selected by AIC, AICc, BIC, or HQC using effective ridge degrees of
freedom.

Preprocessing provides interval scaling and optional ordinary differencing
selected by the package's KPSS procedure.

Recursive forecasts retain the final reservoir state, propagate predictions or
shocks through future lagged inputs, and optionally form simulation paths,
standard deviations, and central intervals from a moving-block bootstrap of
centered readout residuals.

Expanding-window tuning evaluates grids of leakage rates, spectral radii, and
reservoir-size fractions with horizon-specific MSE and MAE.

Principal modeling interfaces are available through `forecasting_mod`, while
KPSS diagnostics are also available through `diagnostics_mod`.

## Licensing

The translation is licensed under GPL-3.0-only; see `LICENSE-ECHOS`.
