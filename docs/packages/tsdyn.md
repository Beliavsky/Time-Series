# tsDyn

[Back to the implemented-package index](../../README.md#implemented-packages).

## Overview

`tsdyn.f90` translates the
[tsDyn](https://cran.r-project.org/web/packages/tsDyn/index.html) R package.

## Algorithms and Procedures

The smooth-transition layer provides the logistic transition function,
supplied-innovation and shared-RNG LSTAR simulation, variable-projection least
squares estimation with a trimmed threshold and slope grid followed by local
refinement, fitted transition weights, information criteria, and recursive
point forecasts.

Internal lagged or externally supplied transition variables are supported
during estimation.

Numerical-Hessian covariance estimates, standard errors, confidence intervals,
parameter gradients, auxiliary-regression tests for additional
smooth-transition regimes, and exhaustive order and delay selection complete
the two-regime layer.

General additive STAR models support two or more logistic components, profiled
estimation, gradients, additional-regime tests, recursive forecasts, and
refitting after regime expansion.

Local-linear autoregression provides delay embedding, epsilon-grid diagnostics,
Theiler-window exclusion, normalized leave-one-out errors, automatic radius
selection, neighbor counts, fitted values, and recursive forecasts with
optional neighborhood enlargement for isolated states.

Additive autoregression uses centered penalized B-spline lag components,
componentwise GCV smoothing selection and effective degrees of freedom, fitted
component evaluation, and recursive forecasts.

Threshold VECM estimation supports one or two thresholds, grid searches over
cointegrating coefficients and thresholds for two-variable systems, and fixed
normalized cointegrating vectors for systems with three or more variables.

Full-regime or error-correction-only switching, trimming, conditional least
squares, supplied-innovation simulation, level forecasts, residual-bootstrap
refitting, regime classification, covariance inference, and generalized impulse
responses all support the multivariate form.

Exact thresholds and bounded threshold searches are available, along with
supplied deterministic or exogenous contributions to the cointegrating relation
and deterministic terms shared across regimes.

SETAR and threshold VAR routines provide trimmed one- or two-threshold grid
searches, level or momentum transition variables, regime-specific lag orders,
transition-delay and TVAR transition-variable selection, supplied-draw
simulation, and probabilistic forecasts.

Hansen-Seo linearity against threshold cointegration and Seo no-cointegration
tests include draw-driven wild and residual bootstrap distributions.

Paired nonlinear generalized impulse responses are available for threshold VAR
and threshold VECM fits, with normalized forecast-error variance
decompositions.

Fixed-threshold SETAR fits support shared deterministic terms, shared lag
effects, combined sharing, and outer-regime symmetry; TVAR fits support a
jointly estimated common intercept.

Hansen SETAR statistics and multivariate TVAR likelihood-ratio statistics
compare linear, two-regime, and three-regime models and evaluate supplied
null-bootstrap series.

The nonlinear unit-root layer provides the Bec-Ben Salem-Carrasco supremum LR,
Wald, and LM tests over symmetric three-regime SETAR alternatives, together
with the Kapetanios-Shin supremum, average, and exponential-average Wald
statistics over asymmetric thresholds. Both accept supplied null-bootstrap
series for reproducible p-values and are re-exported by `diagnostics_mod`.

Neural-network autoregression uses delay embeddings, sigmoid hidden units, a
linear output, analytic-gradient BFGS training, recursive one-step forecasts,
and AIC/BIC hidden-size selection. Its reusable dense-regression kernel is in
`neural_network.f90`, while the time-series interface is re-exported by
`arma_mod`.

Draw-driven probabilistic forecasts are available for SETAR, LSTAR, general
STAR, local-linear AR, additive AR, and neural-network autoregressions. They
return complete simulated paths, path means, sample standard errors, and
pointwise type-7 empirical intervals, with an option to suppress the first
future innovation.

Gaussian Monte Carlo, indexed residual bootstrap, and block bootstrap forecasts
are obtained by constructing the supplied innovation matrix with `random_mod`
and `resampling_mod`.

Rolling-origin evaluation for the same six nonlinear model families supports
multiple horizons, sequential updates with realized observations, and optional
periodic expanding-window refits.

The shared `rolling_forecast_result_t` aligns targets, origins, forecasts,
actuals, signed errors, validity, and refit events; `rolling_forecast_accuracy`
summarizes each horizon with the library's existing accuracy measures. These
interfaces are re-exported by `forecasting_mod`.

Manzan conditional-independence diagnostics provide correlation-integral delta
statistics with delay embeddings, maximum-norm neighborhoods, and
Theiler-window exclusion.

Grid tests use caller-supplied permutations, while the linearity test compares
the nonparametric statistic with its covariance-eigenvalue benchmark and
simulates an AIC-selected autoregressive null from supplied standard-normal
draws. These diagnostics are re-exported by `diagnostics_mod`.

Residual-bootstrap paths for SETAR, TVAR, and TVECM fits retain the observed
initial lag block and accept caller-constructed innovation sequences. Each path
can optionally re-estimate its thresholds and coefficients.

Indexed, circular moving-block, additive Gaussian or Rademacher, and
conventional multiplicative wild resampling are provided by the reusable
`resampling_mod`; package-specific bootstrap entry points are re-exported by
`arma_mod` and `multivariate_mod`.

Linear AR, SETAR, LSTAR, and additive STAR regime diagnostics report long-run
means, characteristic roots and moduli, effective AR orders, unit-root flags,
and stability outside the unit circle.

Trend models and unit-root regimes mark their equilibrium means as undefined.

The shared ascending-power polynomial root solver is provided by
`polynomial_mod`, and the diagnostics are re-exported by `arma_mod` and
`diagnostics_mod`.

Johansen rank testing supports the five Doornik deterministic-term cases, trace
and maximum-eigenvalue statistics, gamma-approximation p-values,
finite-sample-adjusted trace p-values, specific null ranks, and automatic rank
selection.

Joint rank and lag selection evaluates AIC, BIC, and Hannan-Quinn grids using
likelihood or residual-covariance scoring, with common or lag-specific samples.
These interfaces are re-exported by `multivariate_mod`, and the rank test is
also available through `diagnostics_mod`.

SETAR residual inference provides pooled and regime-specific innovation
variances with OLS or maximum-likelihood denominators.

Conditional SETAR coefficient covariance blocks may use pooled or regime
variances.

TVAR inference combines pooled or regime residual covariance with inverse
regime-design cross-products through Kronecker covariance blocks. Both return
mapped standard errors, active-parameter masks and counts, and accept external
threshold series.

TVECM inference supplies pooled Kronecker or regime-covariance sandwich
estimates for fully switching and shared-short-run specifications, together
with coefficient standard errors, t statistics, and two-sided Student p-values.

The SETAR interfaces are available through `arma_mod` and `diagnostics_mod`;
TVAR and TVECM inference are available through `multivariate_mod`.

Unified regime paths are available for SETAR, LSTAR, additive STAR, TVAR, and
TVECM models. Each path retains full-series alignment, leading-position
validity, discrete regime labels, and continuous logistic transition weights
for smooth-transition models.

Replacement observations and external threshold drivers are supported where
applicable.

Univariate paths are re-exported by `arma_mod`, while multivariate paths are
re-exported by `multivariate_mod`.

Principal entry points are also available through `arma_mod` and
`multivariate_mod`.

## Licensing

The translation is licensed under GPL-2.0-or-later; see `LICENSE-TSDYN`.
