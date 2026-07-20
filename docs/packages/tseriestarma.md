# tseriesTARMA

[Back to the implemented-package index](../../README.md#implemented-packages).

## Overview

`tseriestarma.f90` translates the
[tseriesTARMA](https://cran.r-project.org/web/packages/tseriesTARMA/index.html)
R package.

## Algorithms and Procedures

The foundational layer represents two-regime TARMA models with regime-specific
sparse autoregressive and moving-average lag sets, intercepts, innovation
scales, a delayed threshold variable, and a threshold value. It provides
conditional fitted-value and innovation recursions, regime paths,
supplied-standard-innovation and shared-RNG simulation with burn-in, regimewise
autoregressive starting estimates, conditional least-squares fitting at a fixed
threshold, and sample-trimmed threshold profiling with AIC and BIC.

Analytic conditional-innovation derivatives now drive ordinary least-squares
and frozen-weight optimization.

Conditional covariance estimates accompany ordinary fits.

Gaussian and fixed-degrees-of-freedom Student density-power weights support
trimmed-start IRLS robust refinement, regime-aware residual reconstruction,
robust scale correction, convergence diagnostics, and sandwich covariance
estimates.

The nonlinear-diagnostic layer provides classic and heteroskedasticity-robust
supremum LM tests of an autoregressive null against a threshold-autoregressive
alternative, including the complete threshold profile and maximizing threshold
estimates. It also supports exact ARMA null fits for ARMA-versus-TARMA tests,
with either fixed common MA parameters or joint AR and MA threshold-effect
scores and optional explicit threshold grids.

Bootstrap inference for AR-versus-TAR tests includes IID residual, Rademacher
and Gaussian wild, and Hansen score-perturbation schemes, with supplied-random
and shared-RNG interfaces.

Supremum LM unit-root testing covers an integrated MA(1) null against a
stationary TARMA(1,1) alternative, including the package's complete tabulated
critical-value array and IID, Rademacher-wild, and Gaussian-wild bootstrap
inference.

The GARCH-aware layer jointly estimates Gaussian ARMA-GARCH or ARMA-ARCH null
models, propagates threshold effects through both the innovation and
conditional-variance scores, and supplies the package's `ACValues` critical
values for ARMA-versus-TARMA supLM tests.

Principal forecasting interfaces provide zero-future-innovation point paths and
simulation-based means, standard deviations, arbitrary marginal quantiles,
complete paths, and regime paths from supplied or shared-RNG standard-normal
draws.

The `TARMA.fit2` layer builds common AR, regime-specific threshold AR, optional
common or regime intercept, and external-regressor designs, then reuses
`astsa_mod` for exact common regular and seasonal MA likelihood fitting,
conditional covariance inference, and warm-started threshold profiling. It
returns an ordinary TARMA forecast adapter when an internal threshold is used.

Interfaces are available through `arma_mod`, `forecasting_mod`, and
`regression_time_series_mod`.

## Licensing

The translation is licensed under GPL-3.0-or-later; see `LICENSE-TSERIESTARMA`.
