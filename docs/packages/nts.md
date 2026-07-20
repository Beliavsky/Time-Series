# NTS

[Back to the implemented-package index](../../README.md#implemented-packages).

## Overview

`nts.f90` translates the
[NTS](https://cran.r-project.org/web/packages/NTS/index.html) R package.

## Algorithms and Procedures

Its threshold autoregressive layer provides univariate and multivariate
self-exciting TAR simulation for two or more regimes, conditional least-squares
estimation with fixed thresholds, exhaustive trimmed threshold search, and
simulation-based forecasts with pointwise prediction intervals.

Multivariate models support regime-specific VAR orders and innovation
covariances, delayed component or external threshold variables, and AIC or
determinant threshold selection.

Supplied-innovation and supplied-forecast-draw cores support deterministic
testing, while public wrappers use the shared random stream.

The module also provides Markov-switching autoregressive simulation with
arbitrary regime transition matrices, the arranged recursive least-squares
threshold test, Tsay's quadratic nonlinearity test, rolling-origin SETAR
evaluation, rank-based and determinant portmanteau tests, and a selected-lag
quadratic F test.

Time-varying AR models use shared state-space filtering, smoothing, and
likelihood optimization, while random-coefficient AR models provide likelihood
estimates, Hessian errors, and standardized residuals.

Multivariate TAR fits can be restricted by coefficient t ratios.

Autoregressive conditional-mean count models support exogenous predictors and
Poisson, negative-binomial, and double-Poisson likelihoods.

Continuous functional autoregression supports shared-RNG and supplied-draw
simulation, spline-convolution estimation with an OU spatial covariance,
sequential order tests, recursive and partial-curve forecasting, and
registration-based estimation and tests for irregular curves.

Generic sequential Monte Carlo supports model callbacks, supplied proposal and
resampling draws, scheduled stratified resampling, effective sample sizes,
delayed state estimates, retained ancestry, genealogical smoothing, and
shared-RNG execution.

A companion Rao-Blackwellized filter propagates nonlinear particles together
with conditional Gaussian means and covariances.

Relevant entry points are also available through `arma_mod`,
`multivariate_mod`, `markov_switching_mod`, `state_space_mod`,
`diagnostics_mod`, `count_time_series_mod`, and `functional_time_series_mod`.

## Licensing

The translation is licensed under GPL-2.0-or-later; see `LICENSE-NTS`.
