# itsmr

[Back to the implemented-package index](../../README.md#implemented-packages).

## Overview

`itsmr.f90` translates the
[itsmr](https://cran.r-project.org/web/packages/itsmr/index.html) R package.

## Algorithms and Procedures

It provides theoretical ARMA autocovariances from an innovations-form state
representation and the innovations algorithm for MA estimation. It also
provides Hannan-Rissanen ARMA estimation with ITSMR's preliminary-order rule,
regression standard errors, innovation variance, and corrected AIC.

Exact innovations maximum likelihood uses stable partial-autocorrelation
coordinates, BFGS optimization, Hessian standard errors, and corrected-AIC
order selection.

Yule-Walker estimation is centralized in `time_series_stats.f90` for reuse by
ASTSA and ITSMR; existing sample correlations, simulation, smoothing, and
periodogram code are not duplicated.

The graphics-independent ARMA forecast returns recovered innovations,
MA-infinity weights, recursive point forecasts, forecast standard errors, and
configurable normal prediction intervals.

The ARAR implementation performs ITSMR's iterative memory shortening,
exhaustive four-lag sparse AR selection through lag 26, composed-filter
forecasting, and impulse-response prediction intervals. Its result retains the
selected lags, coefficients, memory polynomial, and final filter.

Burg autoregression is implemented in the shared statistics module using
forward-backward errors and reflection coefficients.

The ITSMR adapter adds asymptotic coefficient standard errors and
innovations-based variance and AICc.

ARMA models also expose AR-infinity polynomial coefficients and exact residual
output containing fitted values, innovations, time-specific prediction
variances, and standardized innovations for downstream diagnostics.

Shared harmonic regression fits optional polynomial trend terms together with
ITSMR-phase cosine and sine pairs. It returns OLS inference, fitted signal,
residuals, and direct out-of-sample extrapolation.

Typed transformed forecasting composes regular or seasonal differencing,
harmonic regression, polynomial trend, classical seasonal adjustment, and an
optional outer log transform.

Difference polynomials are included in forecast uncertainty before every fitted
transform is reversed.

The shared Fourier layer provides normalized direct transforms for real series.

ITSMR spectral rank filtering retains the strongest positive-frequency bins
with conjugate symmetry and returns their indices, frequencies, amplitudes,
coefficients, and real reconstruction.

ITSMR residual diagnostics provide Ljung-Box, McLeod-Li, turning-point,
difference-sign, and rank tests with numerical p-values and no plotting layer.

Compatibility smoothers preserve ITSMR's replicated-endpoint moving average,
recursive exponential initialization, symmetric-index Fourier low-pass filter,
and classical seasonal-component behavior.

The distinct numerical exports translated in `itsmr.f90` are `aacvf`, `ia`,
`hannan`, `arma`, `autofit`, `forecast`, `burg`, `Resid`, `test`, `arar`,
`smooth.rank`, `smooth.fft`, `smooth.ma`, `smooth.exp`, `season`, `hr`,
`ma.inf`, and `ar.inf`. Their Fortran APIs use descriptive names and derived
results rather than R lists and string vectors.

Sample autocovariance, Yule-Walker estimation, ARMA stability checks,
simulation, periodograms, polynomial trends, and Fourier construction reuse
`forecast_mod`, `time_series_stats_mod`, `astsa_mod`, or the shared Fourier
module.

`plotc`, `plota`, and `plots` are graphics-only and intentionally omitted.

`specify` is represented directly by `itsmr_arma_model_t`; package datasets and
`selftest` are not library algorithms.

## Licensing

The translation is licensed under BSD-2-Clause.
