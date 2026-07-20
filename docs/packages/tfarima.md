# tfarima

[Back to the implemented-package index](../../README.md#implemented-packages).

## Overview

`tfarima.f90` translates the
[tfarima](https://cran.r-project.org/web/packages/tfarima/index.html) R
package.

## Algorithms and Procedures

It provides the transfer-function ARIMA layer with increasing-lag polynomial
multiplication, division, GCD, powers, rational expansions, and derivative
evaluation. It also provides delayed rational filtering, impulse and step
responses, lag differencing, pulse/step/ramp interventions, seasonal dummies,
harmonic regressors, ARIMA forecast and backcast recursions, theoretical
autocovariances and partial autocorrelations, and Newton autocovariance-to-MA
factorization.

ARMA covariance calculations reuse `itsmr_mod`.

Restricted sparse and powered lag polynomials use explicit offset and loading
matrices for linear parameter constraints.

Unrestricted and restricted rational transfer functions can be estimated by a
profiled conditional Gaussian likelihood using a Nelder-Mead warm start, shared
BFGS refinement, Hessian inference, and a caller-supplied ARMA noise filter.

The outlier layer constructs innovation, additive-outlier, level-shift, and
temporary-change responses in the ARIMA residual domain. It supports automatic
or caller-selected timing, repeated residual cleanup, configurable screening
and retention thresholds, and joint least-squares effect and t-ratio refitting.

Monthly calendar support provides Sunday-through-Saturday counts, all six
`CalendarVar` coding forms, reference-day and working-day contrasts,
month-length and leap-year adjustments, Gregorian Easter and Easter-Monday
windows, and direct forecast-horizon extension through the requested number of
monthly observations.

Calendar arithmetic reuses `calendar_mod`.

UCARIMA support represents normalized independent ARIMA components, constructs
their common denominator and lifted finite numerators, rejects non-identifiable
shared denominator factors, and obtains the aggregate innovation model through
Cramer-Wold autocovariance factorization.

Wiener-Kolmogorov filters expose the exact rational symmetric numerator and
aggregate MA denominator, expand finite symmetric weights, and decompose
complete samples using ARIMA forecast and backcast endpoint extensions.

Extended Euclidean polynomial GCDs return Bezout coefficients, and a dense
coprime partial-fraction solver supports Hillmer-Tiao conversion of an
aggregate ARIMA spectrum into admissible or canonical UCARIMA components.

The conversion records its Wold-domain fractions and verifies spectral
reconstruction after Cramer-Wold factorization.

UCARIMA components can also be assembled into a block innovations state-space
model.

Stationary blocks use iterated Lyapunov covariances, integrated blocks use KFAS
diffuse initialization, and the shared filter and smoother return component
estimates, conditional variances, forecasts, and normal intervals.

Exact UCARIMA estimation uses the same ordinary or diffuse KFAS likelihood.

Component variances use log-scale coordinates, optional masks select AR and MA
coefficients, optional regressors are estimated jointly, and finite-difference
Hessian inference is reported on the natural parameter scale.

Exact transfer-function estimation supports multiple delayed rational inputs,
joint deterministic regressors, free transfer and ARIMA-noise masks, and
ordinary or diffuse KFAS likelihoods.

Known future inputs and regressors feed forecast means, while the fitted ARIMA
noise supplies forecast uncertainty and normal prediction intervals.

Prewhitened cross-correlations expose signed lags, significance flags, and
scaled impulse weights.

Transfer identification can select the first significant nonnegative delay and
recover rational numerator and denominator starting values, while fitted
multi-input models provide an input-by-input residual cross-correlation check.

Backward transfer-model selection uses the exact-fit covariance matrix for
two-sided normal tests, removes the least significant eligible regressor or
complete dynamic input, and refits after every step.

Caller keep masks protect required terms, and the result retains
original-variable masks, p-values, and the ordered removal history.

Transfer diagnostics combine standardized-residual ACF and PACF values,
classical and weighted Ljung-Box tests, Jarque-Bera normality, a cumulative
periodogram white-noise check, and prewhitened residual CCFs for every dynamic
input.

Portmanteau calculations are shared with `tsissm` through
`time_series_diagnostics_mod`.

Transfer simulation evaluates complete multi-input and regression signals and
combines them with recursively generated ARIMA noise.

Exact fitted models can use fixed inputs or independently simulated ARIMA
inputs, multiple paths, burn-in, reproducible shared random seeds,
caller-supplied innovations, and presample series and innovation histories.

Results retain output, noise, innovations, inputs, and each deterministic
signal component.

Robust Cramer-Wold factorization supports Newton moment solving, minimum-phase
palindromic roots, Bauer initialization, Laurie AS 175, Wilson iteration, and
automatic residual-based selection.

UCARIMA construction and ARIMA-to-UCARIMA conversion use the automatic fallback
and report the chosen method, iterations, convergence, and reconstruction norm.

Public conversions map between palindromic and Wold polynomial coordinates.

ARIMA model algebra adds or subtracts independent model spectra by forming
least-common AR and differencing denominators, lifting the MA numerators,
factoring the resulting covariance polynomial, and cancelling common factors.

Finite PSI and PI expansions optionally include differencing, report cumulative
weights, and propagate an AR-then-MA coefficient covariance matrix to ordinary
and cumulative standard errors by the delta method.

Time-invariant univariate state-space models can be reduced in the opposite
direction through the Leverrier-Faddeev characteristic polynomial and adjoint.

The reduction separates stationary and unit roots, constructs observation and
state-disturbance numerator polynomials with explicit KFAS timing, aggregates
their joint covariance by Cramer-Wold factorization, and cancels common ARIMA
factors.

Covariance-aware state-space forms can also switch between contemporaneous and
one-lag state-disturbance conventions.

The alternative root decomposition groups inverse AR and differencing roots as
trend, seasonal, exponential, or cyclical effects and constructs the TFARIMA
deterministic root basis.

Exact Gaussian presample ARMA residuals feed forecast, backcast, and mixed
effect recursions.

Results include all four components, the irregular effect, reconstruction
diagnostics, and additive or log-scale seasonal adjustment.

A validated Aberth root fallback handles sparse seasonal polynomials that are
unstable under plain Durand-Kerner iteration.

Direct ARIMA-to-structural-state-space conversion uses the eventual-forecast
root basis and PSI weights to construct either single-source or multiple-source
disturbance forms.

Single-source models retain their rank-one joint disturbance covariance.

Multiple-source models match the finite MA autocovariances with optional
variance grouping, report unrestricted admissibility, and use a nonnegative
least-squares fallback when variance matching would otherwise produce negative
disturbances. Both contemporaneous and one-lag forms are available and
round-trip through the state-space-to-ARIMA reduction.

Joint-disturbance structural models can be initialized by filtered GLS and run
through a lag-timed Kalman innovations filter that retains observation/state
disturbance correlation in its gain.

The associated information smoother returns observation-aligned conditional
states and covariance matrices.

Filtering reports exact Gaussian likelihoods and standardized innovations;
multi-step forecasts provide latent moments, optional deterministic regression,
and lognormal original-scale moments.

Known regressors and log transformations are supported consistently by
initialization, filtering, and smoothing.

The alternative exact reduced-form likelihood applies the aggregate stationary
AR polynomial to the differenced observations and factors the resulting banded
covariance directly. It supports multiple correlated finite disturbance
numerators, profiles the common innovation scale, and returns raw whitened and
determinant-normalized residuals.

Public wrappers accept built UCARIMA models or structural state-space forms,
while reusable banded Cholesky and forward-solve routines expose the underlying
numerical operations.

## Licensing

The translation is licensed under GPL-2.0-or-later.
