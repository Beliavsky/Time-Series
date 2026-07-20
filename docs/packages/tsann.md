# TSANN

[Back to the implemented-package index](../../README.md#implemented-packages).

## Overview

`tsann.f90` translates the
[TSANN](https://cran.r-project.org/web/packages/TSANN/index.html) R package.

## Algorithms and Procedures

It derives a maximum autoregressive order from the sample ACF and searches the
Cartesian grid of lag orders and single-hidden-layer widths using the shared
feed-forward neural-network ensemble.

The search reports training, validation, and test RMSE for every candidate,
refits each candidate on the combined training and validation prefix, and
retains aligned fitted values and recursive test forecasts.

Chronological blocks are disjoint and validation RMSE is the default selection
criterion, correcting the upstream routine's overlapping validation sample and
test-set selection; explicit test-RMSE selection remains available for
compatibility.

Principal interfaces are available through `forecasting_mod`.

## Licensing

The translation is licensed under GPL-3.0-only; see `LICENSE-TSANN`.
