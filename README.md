# Fortran Time-Series Library

Fortran 2008 translations of numerical algorithms from R time-series packages.
The library currently contains package modules including `forecast_mod`,
`kfas_mod`, `urca_mod`, `astsa_mod`, `itsmr_mod`, `arfima_mod`, and
`fracdiff_mod`, `fracdist_mod`, `nsarfima_mod`, `garma_mod`, `esemifar_mod`, `tfarima_mod`,
`robustarima_mod`, `bsts_mod`, `bssm_mod`, `mar_mod`, `fcvar_mod`, and
`bigtime_mod`, `bigvar_mod`, `varshrink_mod`, `vars_mod`, `var_etp_mod`,
`bvartools_mod`, `nts_mod`, `mswm_mod`, `tsdyn_mod`, `expar_mod`,
`exparma_mod`, `bentcablear_mod`, `baystar_mod`, `mixar_mod`, `count_time_series_mod`, and
`functional_time_series_mod`.

## Shared Infrastructure

- `kind.f90`: shared `kind_mod` module defining the double-precision kind `dp`.
- `utils.f90`: IEEE NaN and option-string helpers.
- `calendar.f90`: Gregorian `date_t`, parsing, arithmetic,
  comparisons, ISO weekdays, ordinal days, leap years, and Easter dates.
- `linalg.f90`: inversion, log determinants, Cholesky and symmetric
  eigen decompositions, Householder-reduced general eigenvalues, rank
  estimation, and matrix helpers.
- `stats.f90`: sorting, quantiles, descriptive statistics, ordinary least
  squares, and regression RSS.
- `spline.f90`: quantile-knot B-spline bases, finite-difference penalties,
  penalized least squares, covariance estimates, effective degrees of freedom,
  GCV smoothing selection, and prediction.
- `random.f90`: seeding, uniform and standard-normal generation,
  normal matrices, and multivariate-normal transformations.
- `fourier.f90`, `polynomial.f90`, and `special_functions.f90`: reusable
  transforms, polynomial operations, and distribution-related functions.
- `optimization.f90`: finite-difference BFGS with Armijo search.
- `time_series_stats.f90`: ACF, PACF, CCF, AR estimation, and harmonic
  regression.
- `time_series_diagnostics.f90`: weighted univariate portmanteau tests and
  FCVAR-style univariate and multivariate white-noise Q and robust LM tests.

Package implementations reuse this layer and earlier package modules. In
particular, the astsa state-space adapter is implemented over `kfas_mod`.

## Topic Facades

Package modules preserve translation provenance and expose complete
package-specific APIs. The following implementation-free facade modules provide
curated, topic-based entry points for applications and examples:

- `arma_mod`: AR, threshold AR, MA, ARMA, ARIMA, and SARIMA models.
- `long_memory_mod`: fractional differencing and long-memory models.
- `state_space_mod`: filtering, smoothing, structural models, and particle methods.
- `multivariate_mod`: VAR, VARMA, VECM, cointegration, and factor models.
- `volatility_mod`: ARCH, multivariate GARCH, and stochastic volatility.
- `spectral_mod`: Fourier transforms, periodograms, and spectral estimation.
- `diagnostics_mod`: residual tests, portmanteau tests, and forecast diagnostics.
- `forecasting_mod`: univariate, multivariate, and structural forecasting.
- `regression_time_series_mod`: dynamic regression, transfer functions,
  interventions, and calendar regressors.
- `bayesian_time_series_mod`: Bayesian estimation, MCMC, and predictive methods.
- `markov_switching_mod`: Markov-switching filtering, smoothing, and estimation.
- `functional_time_series_mod`: functional autoregression, estimation, tests,
  and forecasting.

Facade modules contain no numerical implementations. New package translations
should remain in package modules and add their principal user-facing procedures
to the relevant facades. Because some algorithms belong to multiple topics,
applications importing several facades should normally use explicit `only:`
lists.

Run `python tools/fortran_style_audit.py` to audit all Fortran sources for the
project's mechanical style rules and CMake coverage. Add `--suggest-purity` for
non-failing, heuristic `pure` and `elemental` candidates, or `--json` for a
machine-readable report. The command exits with status 1 when definite style
errors are found.

## Source Licensing

This is a multi-license source tree. Each Fortran source starts with an
`SPDX-License-Identifier` naming its applicable license and an
`SPDX-FileComment` recording its origin. Package translations retain the
upstream package license. Original shared infrastructure is MIT licensed, and
tests use the license of the translated module they exercise.

New Fortran sources must put these two comments on lines 1 and 2. The style
audit reports `F011` or `F012` when either comment is missing. The
Package-specific notices and license texts are retained under `licenses/`.

## Implemented Packages

`tsdyn.f90` translates the
[tsDyn](https://cran.r-project.org/web/packages/tsDyn/index.html) R package.
The smooth-transition layer provides the logistic transition function,
supplied-innovation and shared-RNG LSTAR simulation, variable-projection least
squares estimation with a trimmed threshold and slope grid followed by local
refinement, fitted transition weights, information criteria, and recursive
point forecasts. Internal lagged or externally supplied transition variables
are supported during estimation. Numerical-Hessian covariance estimates,
standard errors, confidence intervals, parameter gradients, auxiliary-regression
tests for additional smooth-transition regimes, and exhaustive order and delay
selection complete the two-regime layer. General additive STAR models support
two or more logistic components, profiled estimation, gradients, additional-
regime tests, recursive forecasts, and refitting after regime expansion.
Local-linear autoregression provides delay embedding, epsilon-grid diagnostics,
Theiler-window exclusion, normalized leave-one-out errors, automatic radius
selection, neighbor counts, fitted values, and recursive forecasts with optional
neighborhood enlargement for isolated states.
Additive autoregression uses centered penalized B-spline lag components,
componentwise GCV smoothing selection and effective degrees of freedom, fitted
component evaluation, and recursive forecasts. Threshold VECM estimation
supports one or two thresholds, grid searches over cointegrating coefficients
and thresholds for two-variable systems, and fixed normalized cointegrating
vectors for systems with three or more variables. Full-regime or error-
correction-only switching, trimming, conditional least squares, supplied-
innovation simulation, level forecasts, residual-bootstrap refitting, regime
classification, covariance inference, and generalized impulse responses all
support the multivariate form. Exact thresholds and bounded threshold searches
are available, along with supplied deterministic or exogenous contributions to
the cointegrating relation and deterministic terms shared across regimes.
SETAR and threshold VAR routines provide trimmed one- or two-threshold grid
searches, level or momentum transition variables, regime-specific lag orders,
transition-delay and TVAR transition-variable selection, supplied-draw
simulation, and probabilistic forecasts. Hansen-Seo linearity
against threshold cointegration and Seo no-cointegration tests include
draw-driven wild and residual bootstrap distributions. Paired nonlinear
generalized impulse responses are available for threshold VAR and threshold
VECM fits, with normalized forecast-error variance decompositions.
Fixed-threshold SETAR fits support shared deterministic terms, shared lag
effects, combined sharing, and outer-regime symmetry; TVAR fits support a
jointly estimated common intercept. Hansen SETAR statistics and multivariate
TVAR likelihood-ratio statistics compare linear, two-regime, and three-regime
models and evaluate supplied null-bootstrap series.
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
future innovation. Gaussian Monte Carlo, indexed residual bootstrap, and block
bootstrap forecasts are obtained by constructing the supplied innovation matrix
with `random_mod` and `resampling_mod`.
Rolling-origin evaluation for the same six nonlinear model families supports
multiple horizons, sequential updates with realized observations, and optional
periodic expanding-window refits. The shared `rolling_forecast_result_t` aligns
targets, origins, forecasts, actuals, signed errors, validity, and refit events;
`rolling_forecast_accuracy` summarizes each horizon with the library's existing
accuracy measures. These interfaces are re-exported by `forecasting_mod`.
Manzan conditional-independence diagnostics provide correlation-integral delta
statistics with delay embeddings, maximum-norm neighborhoods, and Theiler-window
exclusion. Grid tests use caller-supplied permutations, while the linearity test
compares the nonparametric statistic with its covariance-eigenvalue benchmark
and simulates an AIC-selected autoregressive null from supplied standard-normal
draws. These diagnostics are re-exported by `diagnostics_mod`.
Residual-bootstrap paths for SETAR, TVAR, and TVECM fits retain the observed
initial lag block and accept caller-constructed innovation sequences. Each path
can optionally re-estimate its thresholds and coefficients. Indexed, circular
moving-block, additive Gaussian or Rademacher, and conventional multiplicative
wild resampling are provided by the reusable `resampling_mod`; package-specific
bootstrap entry points are re-exported by `arma_mod` and `multivariate_mod`.
Linear AR, SETAR, LSTAR, and additive STAR regime diagnostics report long-run
means, characteristic roots and moduli, effective AR orders, unit-root flags,
and stability outside the unit circle. Trend models and unit-root regimes mark
their equilibrium means as undefined. The shared ascending-power polynomial
root solver is provided by `polynomial_mod`, and the diagnostics are re-exported
by `arma_mod` and `diagnostics_mod`.
Johansen rank testing supports the five Doornik deterministic-term cases,
trace and maximum-eigenvalue statistics, gamma-approximation p-values,
finite-sample-adjusted trace p-values, specific null ranks, and automatic rank
selection. Joint rank and lag selection evaluates AIC, BIC, and Hannan-Quinn
grids using likelihood or residual-covariance scoring, with common or
lag-specific samples. These interfaces are re-exported by `multivariate_mod`,
and the rank test is also available through `diagnostics_mod`.
SETAR residual inference provides pooled and regime-specific innovation
variances with OLS or maximum-likelihood denominators. Conditional SETAR
coefficient covariance blocks may use pooled or regime variances. TVAR
inference combines pooled or regime residual covariance with inverse
regime-design cross-products through Kronecker covariance blocks. Both return
mapped standard errors, active-parameter masks and counts, and accept external
threshold series. TVECM inference supplies pooled Kronecker or regime-covariance
sandwich estimates for fully switching and shared-short-run specifications,
together with coefficient standard errors, t statistics, and two-sided Student
p-values. The SETAR interfaces are available through `arma_mod` and
`diagnostics_mod`; TVAR and TVECM inference are available through
`multivariate_mod`.
Unified regime paths are available for SETAR, LSTAR, additive STAR, TVAR, and
TVECM models. Each path retains full-series alignment, leading-position
validity, discrete regime labels, and continuous logistic transition weights
for smooth-transition models. Replacement observations and external threshold
drivers are supported where applicable. Univariate paths are re-exported by
`arma_mod`, while multivariate paths are re-exported by `multivariate_mod`.
Principal entry points are also available
through `arma_mod` and `multivariate_mod`. The translation is licensed under
GPL-2.0-or-later; see
`LICENSE-TSDYN`.

`setartree.f90` translates the
[setartree](https://cran.r-project.org/web/packages/setartree/index.html) R
package. It fits global SETAR trees to ordinary regression matrices or lag
windows pooled across multiple series. Every node searches equally spaced
thresholds over all features using binned cumulative `X'X`, `X'y`, `y'y`, and
sample-count sufficient statistics, and terminal
pooled regressions are selected by the package's nested F test, relative SSE
improvement, or both criteria with depth-dependent significance. Flat binary
node storage supports direct traversal, point prediction, leaf residual
scales, normal prediction intervals, and recursive multi-series forecasts.
SETAR forests use sampling without replacement, optional randomized stopping
hyperparameters with independent controls for significance, its per-depth
divider, and the error-improvement threshold, averaged point predictions, and
the package's pooled
within-leaf variance intervals. Optional observation weights propagate through
node regression, cumulative split statistics, effective-size stopping tests,
leaf residual scales, forest subsampling, and pooled forest variances.
Integer-coded categorical predictors use first-seen level ordering and omit the
last level as the reference. Trees and forests retain the training levels,
reproduce their indicator columns during prediction, and reject unseen levels.
Per-series mean normalization and per-window
normalization are supported by the pooled-series interfaces.
Rectangular observation-by-series matrices and reusable ragged
`real_vector_t` collections support equal- or unequal-length input histories.
Principal
interfaces are available through `arma_mod`, `forecasting_mod`, and
`regression_time_series_mod`. The translation is licensed under MIT; see
`LICENSE-SETARTREE`.

`tseriestarma.f90` translates the
[tseriesTARMA](https://cran.r-project.org/web/packages/tseriesTARMA/index.html)
R package. The foundational layer represents two-regime TARMA models with
regime-specific sparse autoregressive and moving-average lag sets, intercepts,
innovation scales, a delayed threshold variable, and a threshold value. It
provides conditional fitted-value and innovation recursions, regime paths,
supplied-standard-innovation and shared-RNG simulation with burn-in, regimewise
autoregressive starting estimates, conditional least-squares fitting at a fixed
threshold, and sample-trimmed threshold profiling with AIC and BIC. Analytic
conditional-innovation derivatives now drive ordinary least-squares and
frozen-weight optimization. Conditional covariance estimates accompany
ordinary fits. Gaussian and fixed-degrees-of-freedom Student density-power
weights support trimmed-start IRLS robust refinement, regime-aware residual
reconstruction, robust scale correction, convergence diagnostics, and sandwich
covariance estimates. The nonlinear-diagnostic layer provides classic and
heteroskedasticity-robust supremum LM tests of an autoregressive null against a
threshold-autoregressive alternative, including the complete threshold profile
and maximizing threshold estimates. It also supports exact ARMA null fits for
ARMA-versus-TARMA tests, with either fixed common MA parameters or joint AR and
MA threshold-effect scores and optional explicit threshold grids. Bootstrap
inference for AR-versus-TAR tests includes IID residual, Rademacher and Gaussian
wild, and Hansen score-perturbation schemes, with supplied-random and shared-RNG
interfaces. Supremum LM unit-root testing covers an integrated MA(1) null
against a stationary TARMA(1,1) alternative, including the package's complete
tabulated critical-value array and IID, Rademacher-wild, and Gaussian-wild
bootstrap inference. The GARCH-aware layer jointly estimates Gaussian
ARMA-GARCH or ARMA-ARCH null models, propagates threshold effects through both
the innovation and conditional-variance scores, and supplies the package's
`ACValues` critical values for ARMA-versus-TARMA supLM tests. Principal
forecasting interfaces provide zero-future-innovation point paths and
simulation-based means, standard deviations, arbitrary marginal quantiles,
complete paths, and regime paths from supplied or shared-RNG standard-normal
draws. The `TARMA.fit2` layer builds common AR, regime-specific threshold AR,
optional common or regime intercept, and external-regressor designs, then
reuses `astsa_mod` for exact common regular and seasonal MA likelihood fitting,
conditional covariance inference, and warm-started threshold profiling. It
returns an ordinary TARMA forecast adapter when an internal threshold is used.
Interfaces are available through `arma_mod`, `forecasting_mod`, and
`regression_time_series_mod`. The
translation is licensed under GPL-3.0-or-later; see `LICENSE-TSERIESTARMA`.

`gmdh.f90` translates the
[GMDH](https://cran.r-project.org/web/packages/GMDH/index.html) R package. It
fits classical group-method-of-data-handling networks by generating pairwise
quadratic neurons, selecting polynomial, sigmoid, radial-basis, or tangent
response transformations, choosing unpenalized-intercept ridge coefficients
by a sequential validation split, and retaining the best neurons at each
self-organized layer. The revised GMDH algorithm adds increasing-prefix linear
feedback neurons to every candidate layer. Separate fitting and recursive
prediction interfaces accompany the package-compatible short-term forecast,
which estimates horizon-specific standard errors from rolling forecast errors
and constructs normal prediction intervals for as many as five steps.
Fitted-value, residual, coefficient, and display interfaces are included, and
the principal routines are also available through `forecasting_mod`. The
translation is licensed under GPL-2.0-or-later; see `LICENSE-GMDH`.

`nnfor.f90` translates the
[nnfor](https://cran.r-project.org/web/packages/nnfor/index.html) R package.
The foundational layer supplies reversible min-max scaling and the package's
symmetric fast sigmoid, then fits fast extreme-learning-machine ensembles from
supplied or shared-RNG hidden weights. Output weights support least squares,
validation-selected ridge and lasso regression, or backward-AIC stepwise
selection, with optional direct linear connections and mean, median, or
Gaussian-KDE mode ensemble combination. Lagged univariate adapters construct
autoregressive training matrices and recursively forecast each ELM member.
One-hidden-layer MLP ensembles reuse `neural_network_mod`, accept supplied or
shared-RNG initial parameters, preserve memberwise fitted and forecast paths,
and provide recursive lagged forecasts. The preprocessing layer applies
sequential ordinary or seasonal differences, builds binary or trigonometric
seasonal inputs, aligns independently lagged exogenous regressors, and retains
every differencing level for exact one-step fitted-value and recursive forecast
reconstruction. Backward-AIC response-lag selection supports forced lags while
retaining deterministic and exogenous inputs. Both ELM and MLP ensembles use
this preprocessing path, and terminal-sample validation selects an MLP hidden
size across repeated deterministic starts. A pure supplied-fold interface and
shared-RNG random holdout or five-fold interfaces provide the package's other
hidden-size selection paths. The ELM-specific automatic heuristic builds an
oversized random hidden layer, applies backward-AIC output selection, counts
Student-t-significant hidden coefficients in every repetition, and selects
their rounded median. Centered-moving-average detrending and a tie-corrected
Friedman diagnostic detect deterministic seasonality. The `tsutils`
magnitude-level correlation test automatically chooses additive or
multiplicative detrending from the strongest seasonal position. Joint
trigonometric Canova-Hansen testing then distinguishes stable seasonal inputs
from a required seasonal difference, using Bartlett Newey-West covariance
estimation and the `uroot` response-surface critical values. A dummy-season
formulation with the published raw critical table is also available. A
level-versus-trend exponential-smoothing AIC comparison selects ordinary
differencing.
High-level automatic ELM and MLP fitters combine supplied or selected
differences, default or supplied candidate lags, forced-lag-aware backward-AIC
selection, deterministic seasonal inputs, exogenous regressors, hidden-size
selection, final ensemble fitting, and recursive forecasting. Their result
objects retain every selection diagnostic and support `display`. Principal
fitters also accept multiple seasonal periods with an independent binary,
trigonometric, or disabled input specification for each period. Automatic
multi-season fitting assesses each undifferenced period, uses trigonometric
inputs by default, and permits deterministic or exogenous inputs to remain when
lag selection removes every response lag. All seasonal phases are propagated
through recursive ELM and MLP forecasts. The shared neural-network engine and
the package-facing MLP fitters support vectors of hidden-layer sizes, analytic
multilayer backpropagation gradients, and reproducible supplied starts.
Multilayer ELM fitters propagate fixed random or supplied weights through every
sigmoid hidden layer and estimate the final output with least squares, ridge,
lasso, or stepwise selection. Direct linear connections, ensemble combination,
preprocessing, and recursive forecasts are available for both multilayer model
families. Shared-RNG ELM constructors optionally orthogonalize every feasible
bias-augmented hidden-weight matrix using twice-reorthogonalized modified
Gram-Schmidt, while leaving supplied weights unchanged. Fixed-weight refitting
recomputes scaling, fitted values, residuals, and MSE on compatible new data;
specification-preserving retraining carries lags, differencing, seasonal and
exogenous inputs, hidden architecture, ensemble settings, and regularization
into a newly estimated ELM or MLP ensemble. ELM and MLP THieF callback adapters
fit and forecast individual temporal-aggregation levels and align fitted values
and residuals to the complete input sample for subsequent reconciliation by a
temporal-hierarchy implementation. Principal time-series interfaces are also
available through
`forecasting_mod`. The translation is licensed under GPL-3.0-only; see
`LICENSE-NNFOR`.

`narfima.f90` translates the
[narfima](https://cran.r-project.org/web/packages/narfima/index.html) R package.
Its common neural autoregressive model combines consecutive response lags,
ordinary and seasonal baseline-innovation lags, and optional contemporaneous
regressors. Response, innovation, and regressor transformations retain their
scales for forecasting, and transformations may be disabled. Repeated
one-hidden-layer networks are combined by their mean, with optional direct
linear connections optimized jointly with the nonlinear weights. Automatic
lag selection uses seasonally adjusted Yule-Walker fits. Package-facing
constructors obtain regression-aware baseline innovations from the library's
ARFIMA, ARIMA, jointly backfitted Bayesian structural time-series, or naive
implementations without duplicating the neural fitting layer. Missing baseline
inputs are interpolated, while the common neural model retains complete-case
alignment for lagged predictors.
Recursive point forecasts use zero future innovations; supplied innovations,
Gaussian simulation, and residual bootstrap produce complete forecast paths,
sample standard deviations, and type-8 empirical prediction intervals.
Principal interfaces are available through `long_memory_mod` and
`forecasting_mod`. The translation is licensed under GPL-3.0-only; see
`LICENSE-NARFIMA`.

`nlints.f90` translates the
[NlinTS](https://cran.r-project.org/web/packages/NlinTS/index.html) R package.
It provides discrete Shannon and joint entropy, bivariate and multivariate
mutual information, and lagged transfer entropy with selectable logarithm
bases and package-compatible normalization. Continuous information measures
include Kozachenko maximum-norm entropy, both NlinTS KSG mutual-information
variants, and KSG transfer entropy. The multivariate neural autoregression
retains variable-major lag construction and columnwise min-max scaling, supports
arbitrary hidden-layer sizes, optional biases, linear, sigmoid, ReLU, and tanh
activations, and mini-batch SGD or bias-corrected Adam training. Models retain
aligned fitted values, residuals, normalized residual sums of squares, retained
scales, incremental training, prediction, next-step forecasts, and package-shaped
rolling forecast tables whose final row is out of sample. Versioned text
persistence retains the complete architecture, scaling state, weights, and SGD
or Adam optimizer state. Classical and neural Granger tests report the causality
index, F statistic, degrees of freedom, and beta-ratio p-value. Classical tests
may infer separate integration orders by repeated trend ADF decisions, difference
each series accordingly, and align them to the common maximum-order sample; the
package ADF entry point reuses `urca_mod`.
Principal interfaces are available through `diagnostics_mod`, `forecasting_mod`,
and `multivariate_mod`. The translation is licensed under GPL-2.0-or-later; see
`LICENSE-NLINTS`.

`tslstmplus.f90` translates the
[TSLSTMplus](https://cran.r-project.org/web/packages/TSLSTMplus/index.html) R
package without requiring Keras or TensorFlow. The package adapter prepares
response and exogenous lags as recurrent timesteps or flattened features,
supports independent exogenous lag orders, and retains standard or min-max
input and output scaling. One or more LSTM layers use configurable candidate,
cell, and recurrent-gate activations, optional input dropout, and stateful or
independent sequence execution. Optional hidden dense layers precede the final
linear forecast output. Analytic truncated backpropagation through time trains
the complete network with mini-batch SGD, Adam, or RMSprop and MSE or MAE loss;
terminal validation samples, minimum improvement, and patience provide early
stopping. Rolling fitted values and residuals retain their original alignment,
while recursive forecasts accept known future exogenous observations. The
reusable stand-alone recurrent engine is implemented in
`recurrent_network.f90`. Principal interfaces are available through
`forecasting_mod` and `regression_time_series_mod`. The translation is licensed
under GPL-3.0-only; see `LICENSE-TSLSTMPLUS`.

`tsann.f90` translates the
[TSANN](https://cran.r-project.org/web/packages/TSANN/index.html) R package.
It derives a maximum autoregressive order from the sample ACF and searches the
Cartesian grid of lag orders and single-hidden-layer widths using the shared
feed-forward neural-network ensemble. The search reports training, validation,
and test RMSE for every candidate, refits each candidate on the combined
training and validation prefix, and retains aligned fitted values and recursive
test forecasts. Chronological blocks are disjoint and validation RMSE is the
default selection criterion, correcting the upstream routine's overlapping
validation sample and test-set selection; explicit test-RMSE selection remains
available for compatibility. Principal interfaces are available through
`forecasting_mod`. The translation is licensed under GPL-3.0-only; see
`LICENSE-TSANN`.

`echos.f90` translates the
[echos](https://cran.r-project.org/web/packages/echos/index.html) R package.
It constructs exact-density random recurrent reservoirs, rescales them to a
requested general-matrix spectral radius, and generates leaky tanh echo states
from lagged inputs. Only the linear readout is estimated: intercept-unpenalized
ridge candidates are selected by AIC, AICc, BIC, or HQC using effective ridge
degrees of freedom. Preprocessing provides interval scaling and optional
ordinary differencing selected by the package's KPSS procedure. Recursive
forecasts retain the final reservoir state, propagate predictions or shocks
through future lagged inputs, and optionally form simulation paths, standard
deviations, and central intervals from a moving-block bootstrap of centered
readout residuals. Expanding-window tuning evaluates grids of leakage rates,
spectral radii, and reservoir-size fractions with horizon-specific MSE and MAE.
Principal modeling interfaces are available through `forecasting_mod`, while
KPSS diagnostics are also available through `diagnostics_mod`. The translation
is licensed under GPL-3.0-only; see `LICENSE-ECHOS`.

`starvars.f90` translates the
[starvars](https://cran.r-project.org/web/packages/starvars/index.html) R
package. It estimates vector logistic smooth-transition autoregressions by
profile nonlinear least squares or concentrated Gaussian likelihood, including
grid-selected transition starts, conditional coefficient inference, and AIC
and BIC. Forecasting supports recursive point predictions, Gaussian Monte
Carlo paths, and residual bootstrap paths with central intervals. The module
also implements the package's joint multivariate linearity test, Bartlett
long-run variance estimator, recursive common covariance-break CUMSUM test,
and grouped realized covariance construction from prices or returns. Aggregation
periods can be supplied directly as integer group identifiers or derived from
`date_t` values at daily, monthly, quarterly, or yearly frequency. Principal
modeling interfaces are available through `multivariate_mod`, forecasts through
`forecasting_mod`, tests through `diagnostics_mod`, and realized covariance
through `volatility_mod`. The translation is licensed under GPL-3.0-or-later;
see `LICENSE-STARVARS`.

`rugarch.f90` translates the numerical core of the
[rugarch](https://cran.r-project.org/web/packages/rugarch/index.html) R package.
The initial layer provides `ugarchspec`, filtering, likelihood estimation,
Hessian inference, analytic forecasting, and simulation for sGARCH, integrated
GARCH, exponential GARCH, GJR-GARCH, asymmetric power ARCH, truncated
fiGARCH(1,d,1), component GARCH, log-linear realized GARCH, and the Hentschel
fGARCH omnibus family. The
nonfractional variance families support general ARCH and GARCH orders where the
model permits them. Realized GARCH jointly evaluates returns and a positive
realized-variance series through its leverage-aware measurement equation.
Conditional means support ARMA and fractional ARFIMA dynamics, optional means,
ARCH-in-mean effects, and external regressors; external variance regressors are
also supported. Gaussian, variance-standardized
Student-t, generalized-error, Fernandez-Steel skew-normal, skew-Student,
skew-GED, Johnson SU, standardized normal-inverse-Gaussian, generalized
hyperbolic, and generalized-hyperbolic skew-Student innovations are available,
together with persistence, unconditional-variance, half-life, and
information-criterion calculations. Standalone diagnostics translate the
Berkowitz density-calibration test, Pesaran-Timmermann and Anatolyev-Gerko
directional-accuracy tests, Kupiec and Christoffersen VaR tests, and the
conditional expected-shortfall test. Core variance formulas were adapted from
the MIT-licensed GARCH-BFGS project and integrated with this library's shared
optimizer, random-number generator, statistics, and linear algebra. Modeling
interfaces are available through `volatility_mod`, and forecasts through
`forecasting_mod`; calibration and risk tests are available through
`diagnostics_mod`. The translation is licensed under GPL-3.0-only; see
`LICENSE-RUGARCH`.

`rugarch_extensions.f90` adds parametric bootstrap intervals, expanding and
moving-window rolling refits, multi-specification fitting, asymptotic parameter
simulation, multi-model filtering and forecasting, loss-based model confidence
sets, news-impact curves, and Engle-Ng sign and size-bias diagnostics. Bootstrap
paths recursively update the fitted conditional mean and variance after each
simulated innovation. These graphics-independent workflows are also
exported through `volatility_mod`. The translation is licensed under
GPL-3.0-only; see `LICENSE-RUGARCH`.

`rugarch_diagnostics.f90` translates the Hansen-Nyblom parameter-stability,
grouped probability-transform goodness-of-fit, censored-Weibull VaR-duration,
GMM orthogonality, and quartic-kernel Hong-Li specification tests. The tests are
exported through `diagnostics_mod`; fitted rugarch objects retain numerical
per-observation likelihood scores for direct Nyblom testing. The translation is licensed under
GPL-3.0-only; see `LICENSE-RUGARCH`.

`distribution.f90` centralizes the standardized innovation densities,
parameter transformations, and random generation used by `rugarch.f90`.
Modified-Bessel approximations and the generalized inverse-Gaussian
ratio-of-uniforms sampler reuse MIT-licensed numerical work from GARCH-BFGS;
the combined distribution translation is licensed under GPL-3.0-only.

`distribution_fit.f90` provides constrained maximum-likelihood location-scale
fitting and covariance inference for the standardized innovation families in
`distribution.f90`. The translation is licensed under GPL-3.0-only; see
`LICENSE-RUGARCH`.

`expar.f90` translates the
[EXPAR](https://cran.r-project.org/web/packages/EXPAR/index.html) R package.
It evaluates amplitude-dependent exponential autoregressions, estimates their
phi, pi, and exponential-scale parameters by conditional RSS minimization with
finite-difference BFGS, constructs AR-based starting values, selects the order
by AIC, corrected AIC, or BIC, and computes recursive point forecasts. The
principal estimation interfaces are also available through `arma_mod`, and
forecasting is available through `forecasting_mod`. The translation is
licensed under GPL-3.0-only; see `LICENSE-EXPAR`.

`exparma.f90` translates the
[EXPARMA](https://cran.r-project.org/web/packages/EXPARMA/index.html) R package.
It evaluates exponential autoregressive moving-average recursions in which AR
coefficients depend on the previous observation amplitude and MA coefficients
depend on the previous residual amplitude. Conditional RSS estimation uses the
shared finite-difference BFGS optimizer, starting values reuse Hannan-Rissanen
ARMA estimates, and a two-dimensional AIC search selects AR and MA orders. The
principal interfaces are also available through `arma_mod`. The translation
is licensed under GPL-3.0-only; see `LICENSE-EXPARMA`.

`bentcablear.f90` translates the
[bentcableAR](https://cran.r-project.org/web/packages/bentcableAR/index.html)
R package. It provides linear-quadratic-linear bent-cable and broken-stick
bases, regression design matrices, conditional residual and AR-innovation
recursions, profile-deviance surfaces over transition centers and widths,
joint conditional-RSS estimation with independent or autoregressive errors,
AR stationarity checks, and an alternating Yule-Walker fitting path. Analytic
conditional Fisher information supports delta-method confidence intervals for
the critical transition point. The principal interfaces are also available
through `regression_time_series_mod`. The translation is licensed under
GPL-3.0-or-later; see `LICENSE-BENTCABLEAR`.

`baystar.f90` translates the
[BAYSTAR](https://cran.r-project.org/web/packages/BAYSTAR/index.html) R
package. It provides sparse-lag two-regime TAR likelihoods, conjugate Gaussian
coefficient and inverse-gamma variance posteriors, decreasing-prior delay-lag
probabilities, bounded random-walk threshold Metropolis updates, shared-RNG
Gibbs sampling, supplied-innovation and shared-RNG simulation, posterior
summaries, regime means, residual reconstruction, threshold acceptance rates,
modal delay selection, and DIC. Internal or externally supplied threshold
variables are supported. Principal interfaces are available through both
`arma_mod` and `bayesian_time_series_mod`. The translation is licensed under
GPL-2.0-or-later; see `LICENSE-BAYSTAR`.

`mixar.f90` translates the
[mixAR](https://cran.r-project.org/web/packages/mixAR/index.html) R package.
The foundational Gaussian mixture-autoregression layer supports ragged
component orders, conditional component locations, mixture densities and
distribution functions, posterior component probabilities, fitted values,
conditional variances, residuals, log likelihoods, supplied-draw and
shared-RNG simulation, and the package's probability-weighted Kronecker
stability test. Fixed-point EM estimation uses weighted componentwise AR
regressions and reports convergence, AIC, and BIC. Components may independently
use standard-normal or unit-variance Student-t innovations, with mixed-family
filtering, simulation, conditional distributions, and a likelihood-damped
generalized EM estimator with selectively fixed or estimated degrees of
freedom. Exact multi-step Gaussian forecasts enumerate regime paths and
propagate shifts, AR states, innovation scales, and path probabilities;
simulation forecasts support every innovation family and provide empirical
means, variances, quantiles, and central intervals. Principal interfaces are
available through `arma_mod` and `forecasting_mod`. The regression layer fits
shared or component-specific covariate effects with Gaussian or Student-t MAR
errors, and supplies filtering, likelihood evaluation, simulation, exact
Gaussian forecasts, and simulation forecasts for known future regressors. Its
interfaces are also available through `regression_time_series_mod`. Additive
seasonal MAR models retain separate
ordinary and seasonal ragged coefficients while reusing a sparse expanded-lag
representation for filtering, stability, simulation, and forecasting. Gaussian
and Student-t generalized EM fitting, exact Gaussian forecasts, and simulated
forecast distributions are supported. Observed-Hessian inference covers
ordinary, seasonal, and regression models; residual diagnostics include PIT
and normal-quantile residuals, ACF and Ljung-Box checks, posterior
classification confidence and entropy, and conditional-BIC model comparison.
The multivariate layer represents ragged Gaussian mixture VAR models, evaluates
posterior component probabilities with a stabilized conditional likelihood,
fits component intercepts, VAR matrices, probabilities, and innovation
covariances by EM, and applies the probability-weighted Kronecker stability
criterion. It also provides supplied-draw and shared-RNG simulation,
multi-step path forecasts with point, covariance, and interval summaries, and
multivariate standardized-residual, white-noise, and classification
diagnostics. These interfaces are available through `multivariate_mod`,
`forecasting_mod`, and `diagnostics_mod`. A Bayesian Gaussian MAR layer adds
Dirichlet allocation-weight updates, hierarchical Gamma precision sampling,
Gaussian component-mean updates, stability-constrained random-walk Metropolis
AR updates, reproducible supplied-random and shared-RNG interfaces, posterior
summaries, equal-order label correction, and a posterior-ordinate KDE
approximation to the marginal likelihood. Reversible-jump MCMC selects the
component-specific AR orders under flat, ratio, or Poisson birth/death
proposals, retaining the complete order and parameter trace, posterior order
frequencies, the modal order vector, and a reconstructed modal model. Both
supplied-random and shared-RNG order-selection interfaces are provided.
Data-driven initialization fits component regressions on reproducible
without-replacement subsamples, converts their residuals into soft component
memberships, performs a weighted Gaussian M-step, and contracts coefficients
when needed to obtain a stable starting model. Supplied-index and shared-RNG
initializers support multistart EM fitting, retaining every initialization and
fit while selecting the largest-likelihood converged solution.
Analytical moment routines cover Gaussian and standardized Student-t absolute
and ordinary moments, innovation-mixture moments, arbitrary-order one-step raw
and central conditional moments, kurtosis and excess kurtosis, and explicit
detection of nonexistent heavy-tail moments. Stable MAR stationary means and
covariances are obtained from companion-state first- and second-moment
equations.
Bayesian interfaces are also
available through `arma_mod` and `bayesian_time_series_mod`. The translation is licensed under
GPL-2.0-or-later; see `LICENSE-MIXAR`.

`mswm.f90` translates the
[MSwM](https://cran.r-project.org/web/packages/MSwM/index.html) R package. The
Gaussian linear-model layer provides scaled Hamilton filtering, Kim smoothing,
expected transition probabilities, and EM estimation for arbitrary regime
counts. Regression coefficients may be shared or regime-specific, innovation
variance may be common or switching, and callers may provide initial
coefficients, scales, and transition probabilities. Generalized switching
regression supports Poisson with a log link, Bernoulli-binomial with a logit
link, and Gamma models with log or inverse links, using weighted IRLS within
each EM M-step. Results include conditional means, residuals, filtered and
smoothed state probabilities, decoded states, transition estimates, log
likelihood, AIC, and BIC. Gaussian and generalized fits provide finite-difference
Hessian covariance estimates, natural-scale standard errors, and approximate
confidence intervals. Diagnostics include raw Gaussian residuals and
family-standardized generalized residuals, returned by regime or combined using
smoothed state probabilities. Gaussian fitting also supports selection over
supplied or reproducible random EM starts, with equivalent multistart fitting for
generalized models. Autoregressive wrappers construct response lags alongside
optional exogenous predictors and retain the lag order and terminal response
history for subsequent forecasting. Gaussian and generalized results share
overloaded residual, state-decoding, and concise display interfaces; observation
arrays are omitted from displays unless explicitly requested. The principal
entry points are also available through `markov_switching_mod`. The translation
is licensed under GPL-2.0-or-later; see `LICENSE-MSWM`.

`nts.f90` translates the
[NTS](https://cran.r-project.org/web/packages/NTS/index.html) R package. Its
threshold autoregressive layer provides univariate and multivariate
self-exciting TAR simulation for two or more regimes, conditional least-squares
estimation with fixed thresholds, exhaustive trimmed threshold search, and
simulation-based forecasts with pointwise prediction intervals. Multivariate
models support regime-specific VAR orders and innovation covariances, delayed
component or external threshold variables, and AIC or determinant threshold
selection. Supplied-innovation and supplied-forecast-draw cores support
deterministic testing, while public wrappers use the shared random stream. The
module also provides Markov-switching autoregressive simulation with arbitrary
regime transition matrices, the arranged recursive least-squares threshold
test, Tsay's quadratic nonlinearity test, rolling-origin SETAR evaluation,
rank-based and determinant portmanteau tests, and a selected-lag quadratic F
test. Time-varying AR models use shared state-space filtering, smoothing, and
likelihood optimization, while random-coefficient AR models provide likelihood
estimates, Hessian errors, and standardized residuals. Multivariate TAR fits can
be restricted by coefficient t ratios. Autoregressive conditional-mean count
models support exogenous predictors and Poisson, negative-binomial, and
double-Poisson likelihoods. Continuous functional autoregression supports
shared-RNG and supplied-draw simulation, spline-convolution estimation with an
OU spatial covariance, sequential order tests, recursive and partial-curve
forecasting, and registration-based estimation and tests for irregular curves.
Generic sequential Monte Carlo supports model callbacks, supplied proposal and
resampling draws, scheduled stratified resampling, effective sample sizes,
delayed state estimates, retained ancestry, genealogical smoothing, and
shared-RNG execution. A companion Rao-Blackwellized filter propagates nonlinear
particles together with conditional Gaussian means and covariances.
Relevant entry points are also available through
`arma_mod`, `multivariate_mod`, `markov_switching_mod`, `state_space_mod`,
`diagnostics_mod`, `count_time_series_mod`, and `functional_time_series_mod`.
The translation is licensed under
GPL-2.0-or-later; see `LICENSE-NTS`.

`bvar.f90` translates the [BVAR](https://cran.r-project.org/web/packages/BVAR/index.html)
R package. It
provides the package's diagonal Minnesota row variances, sum-of-coefficients
and single-unit-root dummy observations, conjugate posterior sufficient
statistics, and closed-form log marginal likelihood. Hierarchical lambda,
alpha, innovation-scale, SOC, and SUR tightness estimation uses bounded
random-walk Metropolis sampling with Gamma and inverse-Gamma hyperpriors,
optional burn-in acceptance tuning, a supplied-randomness reproducible core,
and a shared-RNG wrapper. The translation is licensed under GPL-3.0-or-later;
see `LICENSE-BVAR`.

`bvartools.f90` translates the
[bvartools](https://cran.r-project.org/web/packages/bvartools/index.html) R
package. The numerical layer provides Minnesota and coefficient-
inclusion priors, semiautomatic-compatible SSVS updates, Korobilis BVS updates,
multivariate and observation-varying SUR Gaussian posterior moments and draws,
gamma precision updates for measurement and random-walk state innovations,
Primiceri lower-triangular covariance regressions, constant and time-varying
covariance-coefficient posteriors, covariance-vector reconstruction, and
constant-to-TVP SUR expansion. Fixed and George-Sun-Ni semiautomatic SSVS
prior scales use multivariate OLS standard errors, with fixed fallback scales
for contemporaneous covariance coefficients. Constant BVAR and BVEC prior
bundles dimension regular or SSVS coefficient blocks, Wishart or gamma
innovation priors, cointegration priors, and OLS initial values directly from
the raw-series constructors. End-to-end constant BVAR and BVEC fit wrappers
dispatch those bundles to the appropriate covariance and selection samplers;
rank-zero BVEC models are routed through the unrestricted differenced-BVAR
kernel. TVP-BVAR and positive-rank TVP-BVEC bundles and fit wrappers add
random-walk coefficient priors, persistent cointegration states, trajectory
BVS, lower-triangular covariance states with optional BVS, and KSC or OCSN
stochastic volatility. It reuses the shared linear algebra and random
modules. Constant-parameter BVAR and reduced-rank BVEC Gibbs drivers, the pure
KLS cointegration posterior draw for constant or observation-specific
innovation precision, the random-walk coefficient core of the
TVP-BVAR Gibbs sampler, and a reusable time-varying lower-triangular covariance
Gibbs block are also available. Both TVP blocks support trajectory-level BVS.
The joint TVP-BVAR driver alternates coefficient paths, lower-triangular
covariance coefficients, and orthogonal innovation variances in one Gibbs loop.
Covariance coefficients may follow random walks with BVS or remain constant
with SSVS, and their period-specific covariance matrices feed back into every
coefficient-path update.
An end-to-end structural TVP-BVAR wrapper augments the SUR design with the
identified lower-triangular contemporaneous A-model block, provides separate
BVS controls for reduced-form and structural trajectories, and returns those
state blocks separately. Its covariance block supports TVP BVS or constant
SSVS specifications.
Constant structural BVAR and BVEC samplers use the same identified A-model
design. They update contemporaneous coefficients conditionally within the Gibbs
loop, support independent structural SSVS, and retain structural draws and
inverse contemporaneous impact matrices ready for prediction, IRF, and FEVD
routines.
The TVP-BVEC core also accepts the identified contemporaneous SUR block. Its
structural adapter expands the loading and unrestricted-coefficient state
prior, removes structural effects before each cointegration-path update, and
retains structural trajectories, period-specific impact matrices, and BVS
indicators separately. It supports inverse-Wishart errors, diagonal gamma
precisions, gamma errors with time-varying covariance states, and stochastic
volatility.
The TVP-BVAR driver supports selectable KSC or OCSN random-walk stochastic
volatility.
A TVP-BVEC core supports time-varying loading, cointegration, unrestricted
coefficient paths, and lower-triangular time-varying covariance states, with
trajectory-level BVS for unrestricted coefficients and covariance states. Its
TVP-BVEC stochastic volatility also supports both mixture approximations.
A Bayesian dynamic-factor Gibbs driver provides jointly sampled factor paths,
identified factor loadings, factor VAR dynamics, and diagonal measurement and
factor innovation variances. Posterior predictive simulation supports constant
and terminal-state TVP or stochastic-volatility BVAR draws, future exogenous
and deterministic regressors, optional structural contemporaneous transforms,
equal-tail credible intervals, a pure supplied-normal core, and a shared-RNG
wrapper. Draw-wise posterior impulse responses cover forecast-error,
orthogonalized, generalized, structural, and structural-generalized definitions
for constant or selected-period TVP draws, with configurable shock scaling,
cumulative responses, and credible bands. The moving-average recursion reuses
`mts_var_psi`. Posterior orthogonalized, generalized, structural, and
structural-generalized FEVDs reuse those draw-wise responses, provide credible
bands for constant or selected-period TVP models, and optionally normalize
generalized shares within every posterior draw. BVEC posterior coefficients can
be converted to level-VAR form for constant or fully time-varying draws,
including endogenous lag identities, restricted and unrestricted deterministic
terms, exogenous error-correction and differenced-exogenous blocks, rank-zero
and one-lag cases, and reconstruction of levels from first differences.
Raw-series BVAR and BVEC constructors provide aligned endogenous, exogenous,
deterministic, seasonal, holdout, SUR, TVP-SUR, and off-diagonal structural
design matrices while reusing the shared VAR lag builders. The raw-series DFM
constructor reproduces `gen_dfm` sample standardization and model-grid
enumeration while retaining the location and scale needed to transform results
back to the observed units. A dimension-aware DFM prior factory reproduces the
`add_priors.dfmodel` defaults, and a grid Gibbs driver fits every prepared
factor-count and lag-order combination. Posterior model comparison reproduces
the `summary.bvarlist` observation-wise draw averaging and AIC, BIC, and HQ
criteria for constant or time-varying covariance draws. A pure multivariate
Gaussian log-likelihood interface accepts constant or vertically stacked
time-varying covariances. The `kalman_dk` compatibility adapter maps stacked
bvartools arrays to the shared bssm Durbin-Koopman simulation smoother and
returns the full state path through the post-observation state. The translation
is licensed under GPL-2.0-or-later; see
`LICENSE-BVARTOOLS`.

`var_etp.f90` translates the
[VAR.etp](https://cran.r-project.org/web/packages/VAR.etp/index.html) R package.
It provides
Nicholls-Pope asymptotic and residual-bootstrap VAR bias
corrections with Kilian stationarity adjustment; analytic forecast MSE including
coefficient uncertainty; forward/backward bootstrap and bootstrap-after-bootstrap
prediction intervals; constrained system estimates; Wald and nested likelihood-
ratio tests with iid or Mammen wild bootstrap inference; and Kim's improved
augmented predictive regression with Shaman-Stine predictor correction, joint
tests, order selection, and dynamic forecasts. Ordinary VAR estimation, roots,
MA coefficients, and impulse responses reuse `vars_mod`. The translation is
licensed under GPL-2.0-only; see `LICENSE-VAR-ETP`.

`bvars.f90` translates the
[bvars](https://cran.r-project.org/web/packages/bvars/index.html) R package. It
provides
lag-major BVAR data preparation, the package's default random-walk or stationary
matrix-normal/inverse-Wishart prior, matrix-normal simulation, and the complete
homoskedastic Gaussian posterior sampler. Posterior residuals, fitted-density
draws, recursive forecasts with contemporaneous exogenous inputs and common
variance multipliers, and orthogonalized FEVDs are included. FEVD computation
reuses `bvartools_mod`. The homoskedastic multivariate Student-t extension
samples observation-specific inverse-gamma variance multipliers and degrees of
freedom with an adaptive transformed Metropolis step; matching future scale
simulation feeds directly into the forecast routine. Centred and non-centred common stochastic volatility
implement the package's mixture-indicator, persistence, innovation-variance,
hierarchical variance, and latent-state updates. The non-centred sampler uses
ASIS interweaving and generalized-inverse-Gaussian slice draws. Both accept the
dimension-specific auxiliary mixture table explicitly, retain every state and
hyperparameter draw, and simulate future common-variance paths for forecasts.
The Student-t scale-mixture and degrees-of-freedom updates can be combined with
either stochastic-volatility parameterization in a single Gibbs sampler.
Forecasts can condition individual variables and horizons on supplied values;
use NaN for entries left free.
`bvars_sv_auxiliary_mixture` constructs the required auxiliary mixture table
for any number of variables with a reusable deterministic Gaussian-mixture EM
fit; its simulation count and convergence settings are configurable. The
translation is licensed under GPL-3.0-or-later; see
`LICENSE-BVARS`.

`bvarsv.f90` translates the
[bvarsv](https://cran.r-project.org/web/packages/bvarsv/index.html) R package.
It provides the Primiceri time-varying structural VAR
algorithms, including the multivariate
Carter-Kohn random-walk state sampler, packed contemporaneous-matrix and
structural-covariance reconstruction, and the seven-component KSC update for
correlated random-walk log-volatility states. Equation-block Carter-Kohn updates
draw the packed time-varying contemporaneous coefficients. Reusable
inverse-Wishart updates draw the coefficient and volatility innovation
covariances and the equation-block contemporaneous innovation covariances. The
package-specific equation-stacked VAR design and time-varying coefficient-state
update are also provided. The training-sample OLS prior reproduces the package's
GLS coefficient moments, recursive covariance decomposition, initial log
variances, and inverse-Wishart Monte Carlo covariance for contemporaneous states.
`bvarsv_gibbs` assembles these kernels into the complete Primiceri sampler with
configurable burn-in, thinning, and innovation-prior scales, retaining all state
paths and the Q, S, and W covariance draws.
`bvarsv_forecast` produces recursive draw-wise predictive means, covariance
matrices, and realizations with either random-walk parameter drift or terminal
parameters held fixed.
`bvarsv_irf` computes draw-wise responses at any retained state time with the
package's identity, Cholesky, or Primiceri average-volatility impact scenarios.
`bvarsv_predictive_density` evaluates the draw-averaged Gaussian predictive
PDF or CDF at user-supplied points.
`bvarsv_simulate_var1` generates the package's TVP-VAR(1) process with
random-walk coefficients, contemporaneous relations, and log variances.
`bvarsv_predictive_draws` and `bvarsv_parameter_draws` extract forecast draws
and selected time-varying parameter paths without the R package's plotting
layer. The translation is licensed under GPL-2.0-or-later.

`gmvarkit.f90` translates the
[gmvarkit](https://cran.r-project.org/web/packages/gmvarkit/index.html) R
package. It
provides stationary regime means and lag-vector covariances, endogenous
Gaussian and Student-t regime weights, conditional regime and mixture moments,
ARCH scaling for Student-t regimes, and conditional or exact mixture log
likelihoods for arbitrary VAR order.
Recursive simulation retains sampled regimes and future mixing weights for each
path. Monte Carlo forecasting provides means, medians, empirical type-7
quantiles, and expected future regime weights.
`gmvarkit_estimate` performs local finite-difference BFGS likelihood refinement
with Cholesky covariance, softmax weight, and transformed degrees-of-freedom
parameters. Invalid and nonstationary trial models are penalized.
`gmvarkit_genetic_estimate` supplies the package's global-search layer with
elitist selection, crossover, decreasing transformed-space mutation, a recorded
best-objective path, and final BFGS refinement.
`gmvarkit_estimate_constrained` supports general linear maps for stacked AR
coefficients, shared unconditional means across regime groups, and fixed
mixture weights while optimizing only the remaining free coordinates.
Hessian inversion supplies transformed-scale covariance matrices and standard
errors. General linear Wald tests and nested-model likelihood-ratio tests use
chi-square reference distributions.
The inference vector stores each regime's intercept, column-major AR entries,
and rowwise lower-Cholesky entries, followed by weight logits and log-transformed
Student-t degrees of freedom; diagonal Cholesky entries are stored on log scale.
Observationwise finite-difference scores provide OPG information and Rao score
tests. Sequential multivariate quantile residuals support Gaussian and Student-t
regimes; moment diagnostics test normality, serial correlation, and dependence
in squared residuals. Supplying the fitted model and data applies the
Kalliovirta-Saikkonen parameter-estimation correction using residual-moment
derivatives, OPG Fisher information, and moment-score cross covariances;
otherwise the tests use empirical moment covariance matrices.
`gmvarkit_girf` estimates nonlinear generalized impulse responses by paired
simulation with common regimes and innovations, lower-Cholesky identification,
and structural-shock replacement. It also returns responses of the endogenous
regime probabilities. `gmvarkit_gfevd` forms horizon-specific decompositions
from normalized cumulative squared mean GIRFs.
`gmvarkit_linear_irf` computes deterministic regime-specific VAR responses
under lower-Cholesky recursive identification. Selected variables can be
accumulated over horizons, and selected shocks can be normalized to a requested
impact response of a specified variable.
`gmvarkit_unconditional_moments` combines regime means and autocovariances into
the mixture mean, lag-zero-through-p autocovariances, and autocorrelations.
`gmvarkit_pearson_residuals` returns raw residuals or applies the inverse
symmetric square root of each conditional covariance. The elemental
`gmvarkit_information_criteria` evaluates AIC, HQIC, and BIC from a likelihood,
free-parameter count, and effective sample size.
Structural covariance support represents each regime as
`Omega_m = W Lambda_m W'`. `gmvarkit_identify_structural` recovers a common
impact matrix and relative variances through covariance heteroskedasticity and
checks simultaneous diagonalization. Companion procedures construct implied
covariances, reorder or sign-reverse shocks, and change the unit-variance
reference regime. Linear and generalized impulse responses accept this
identification as an optional argument.
`gmvarkit_linear_irf_bootstrap` provides fixed-design Rademacher wild-bootstrap
confidence intervals when regime dynamics are linear. Starting observations
remain fixed, fitted residual vectors receive common rowwise signs, and mixture
weights are held fixed during re-estimation. Structural bootstrap solutions are
matched to the baseline lambda profiles and impact-column signs before their
responses are summarized.
`gmvarkit_girf_inference` adds the outer initial-state distribution used by the
package's full GIRF interface. It supports every observed length-p history,
fixed histories, or stationary Gaussian and Student histories drawn from
selected regimes; it returns individual responses, point estimates, and
equal-tail bounds after optional accumulation and instantaneous or peak
scaling. `gmvarkit_gfevd_inference` calculates a decomposition for every
history and averages them, including decompositions of mixing-weight responses.
`gmvarkit_estimate_structural` optimizes the common impact matrix and positive
relative regime variances directly in the likelihood. Structural restrictions
support exact fixed and zero entries of `W`, positive or negative sign
constraints on free entries, fully fixed lambda matrices, and nonnegative
linear lambda mappings from positive free coordinates. The result includes the
compatible reduced-form model and optional transformed-coordinate Hessian
covariance, gradient, and standard errors.
`gmvarkit_profile_likelihood` evaluates fixed-nuisance likelihood slices for
selected estimation coordinates. Reduced-form and structural overloads use the
same transformed parameterizations as their estimators and return grid values,
log likelihoods, and a validity mask without imposing a graphics dependency.
`gmvarkit_convert_student_regimes` replaces Student regimes above a selected
degrees-of-freedom threshold by Gaussian regimes, orders each regime family by
mixing weight, and returns both directions of the resulting permutation.
Structural variance columns and their reference regime follow the same mapping,
and the converted reduced-form or structural model can optionally be re-fitted.
`gmvarkit_companion_eigenvalues` returns the complex companion roots, moduli,
spectral radii, stationarity flags, and tolerance-based near-unit-root flags for
each regime. `gmvarkit_covariance_eigenvalues` reports ordinary covariance
eigenvalues and near-singularity flags, together with Cholesky-whitened
generalized eigenvalues for every covariance pair. Relative eigenvalue
separations flag weak heteroskedastic identification without printing warnings.
`gmvarkit_multistart_estimate` optimizes an array of starting models, retains
all local fits and convergence information, and ranks successful solutions by
log likelihood. Regimes are ordered canonically within Gaussian and Student
families before parameter comparisons, allowing label-switched and numerically
equivalent fits to be linked through `duplicate_of` rather than counted as
distinct local maxima.
`gmvarkit_structural_multistart_estimate` applies common impact and relative
variance restrictions to paired reduced-form and structural starts. It ranks
fits by structural likelihood, compares their canonical implied reduced forms
to remove regime, shock-order, and shock-sign aliases, and calculates numerical
Hessian inference only for the selected best solution.
`gmvarkit_location_parameters` returns both intercepts and their equivalent
stationary regime means. `gmvarkit_model_from_regime_means` performs the inverse
mapping and returns an ordinary intercept-parameterized model, providing the
package's reversible mean/intercept workflow without introducing ambiguous
state into the likelihood model type. The translation is licensed under
GPL-3.0-only.

`ugmar.f90` translates the
[uGMAR](https://cran.r-project.org/web/packages/uGMAR/index.html) R package.
It specializes the shared Gaussian and Student-t mixture autoregressive engine
to a univariate interface for GMAR, StMAR, and G-StMAR models. The module
provides validated model construction, stationarity tests, endogenous mixing
weights, conditional regime and mixture moments, conditional or exact log
likelihoods, stationary moments, BFGS, genetic, multistart, and restricted
estimation, stationary random and local smart-mutation initializers, recursive
simulation, and Monte Carlo forecasts. Reversible
intercept and stationary-mean parameterizations, numerical Hessian inference,
Wald and likelihood-ratio tests, transformed-coordinate likelihood profiles,
AR characteristic-root reporting, and high-degree Student-to-Gaussian regime
conversion are included. Gaussian
quantile residuals support residual tests with optional parameter-estimation
correction. Principal interfaces are available through `arma_mod`,
`forecasting_mod`, and `diagnostics_mod`. The translation is licensed under
GPL-3.0-only; see `LICENSE-UGMAR`.

`vars.f90` translates the
[vars](https://cran.r-project.org/web/packages/vars/index.html) R package. It
provides
reduced-form VAR estimation with constant, trend,
centered seasonal, and exogenous terms; common-sample AIC, HQ, SC, and FPE lag
selection; manual and sequential t-ratio restrictions; MA coefficients,
orthogonalized responses, and companion roots; asymptotic and adjusted
portmanteau, Breusch-Godfrey, Edgerton-Shukur, multivariate Jarque-Bera, and
multivariate and equationwise ARCH tests; instantaneous-causality statistics;
and OLS-CUSUM fluctuation paths. Structural routines include restricted A, B,
and AB estimation by direct Gaussian likelihood or the package's Fisher-scoring
recursion, parameter standard errors, Blanchard-Quah identification, and SVEC
short- and long-run zero restrictions built from `urca_mod` Johansen results.
The Johansen-to-level-VAR conversion, reduced-form IRF intervals, structural
impact intervals, SVEC bootstrap standard errors, and bootstrap Granger test
are also included. SVAR and SVEC structural impulse responses and forecast-error
variance decompositions share one level-VAR response kernel. Principal entry
points are re-exported by `multivariate_mod` and `diagnostics_mod`. The
translation is licensed under GPL-2.0-or-later; see
`LICENSE-VARS`.

`varshrink.f90` translates the
[VARshrink](https://cran.r-project.org/web/packages/VARshrink/index.html) R
package. It provides
multivariate ridge regression over supplied shrinkage paths, effective degrees
of freedom, generalized cross-validation selection, optional predictor
standardization, and direct VAR-form ridge fitting. The semiparametric Bayesian
regression core implements fixed-lambda conjugate and non-conjugate covariance
priors, Gaussian or multivariate Student-t innovations, iterative latent
weights, covariance estimation, fitted values, residuals, and convergence
diagnostics. The nonparametric estimator includes the Opgen-Rhein/Strimmer
James-Stein covariance calculation used by VARshrink, with separately supplied
or estimated correlation and variance shrinkage intensities. Shrinkage-specific
Gaussian and Student-t log likelihoods and the package's Stein-type
marginal-variance intensity are also included.
Simulation, companion stability, MA representations, impulse responses, and
serial-correlation tests reuse `mts_mod`, `bigtime_mod`, and the shared
diagnostics instead of duplicating those algorithms. The translation is
licensed under GPL-3.0-or-later; see
`LICENSE-VARSHRINK`.

`bigvar.f90` translates the
[BigVAR](https://cran.r-project.org/web/packages/BigVAR/index.html) R package.
It provides the package's
distinct convex Lag and Own/Other structured penalties for VAR models. The Lag
penalty groups complete coefficient matrices by lag. Own/Other separates each
lag into diagonal autoregressive and cross-series groups with BigVAR's
cardinality weights. SparseLag and SparseOO add within-group elementwise
sparsity using BigVAR's `alpha` mixture and its default `1/(k+1)`.
Basic supplies BigVAR's elementwise lasso. Basic, BasicEN, componentwise,
Own/Other, and elementwise HLAG, MCP, and SCAD support separate penalty
parameters, paths, and rolling selections for each response equation.
Relaxed VAR and VARX refitting re-estimates each selected support by least
squares using the shared pseudoinverse and blends it with the penalized fit
through a `refit_fraction` from zero to one.
Direct multi-step VAR fits align each lag vector at origin `t` with the
horizon-specific response at `t+h`. Scalar and response-specific fits, paths,
lambda grids, relaxed refits, rolling validation, and evaluation retain the
selected direct horizon; recursive forecasts remain the default.
An optional per-series Minnesota target places the supplied values on the
own-series first-lag diagonal and estimates penalized deviations from that
matrix. It is supported by VAR, response-specific VAR, VARX, contemporaneous
VARX, direct fits, paths, lambda grids, relaxed refits, rolling validation, and
leave-one-out selection. Targeted fits use a zero intercept and report both the
restored coefficients and target metadata. This shrinkage target is distinct
from the BGR dummy-observation prior and cannot be combined with BGR validation.
Forecast uncertainty includes the fitted innovation covariance, full
horizon-by-horizon VAR forecast-error covariance propagation, standard errors,
and configurable marginal normal intervals. VARX intervals are conditional on
the supplied exogenous path, while direct fits use their horizon-specific
residual covariance.
Componentwise HLAG permits a different maximum lag for each response equation.
Own/Other HLAG additionally prioritizes each response's own lag over
cross-series lags through two nested suffix groups at every lag. Elementwise
HLAG permits a separate maximum lag for every response-predictor pair.
BasicEN combines elementwise lasso and squared-Frobenius shrinkage using the
upstream closed-form proximal update. Tapered implements BigVAR's lag-weighted
lasso directly with weights `lag**alpha`, encouraging stronger shrinkage at
more distant lags. MCP and SCAD use BigVAR's piecewise coordinate updates,
lambda-sized ridge stabilization, centered residual updates, and warm starts
along the regularization path to manage their non-convex objectives.
Structured VARX estimation supports Lag, Own/Other, SparseLag, SparseOO,
BasicEN, MCP, and SCAD with unequal endogenous and exogenous lag orders.
Exogenous series-lag columns form response-vector groups with BigVAR's
`sqrt(k)` weight, while endogenous coefficients retain their corresponding
VAR structure. Joint lambda grids, warm-started paths, fitted diagnostics, and
recursive forecasts are included; ordinary forecasting reuses `bigtime_mod`.
Transfer-function models set the endogenous lag order to zero and estimate
responses solely from lagged or contemporaneous exogenous predictors. Basic,
BasicEN, MCP, and SCAD transfer functions support fixed fits, lambda grids,
warm-started and alpha-specific paths, relaxed refits, conditional forecasts
and intervals, rolling and leave-one-out validation, and rolling reselection.
Fit and path results identify these models through `transfer_function`.
The optional contemporaneous VARX form places `X_t` first in the exogenous
coefficient matrix, followed by `X_(t-1)` through `X_(t-s)`. It supports
current-only predictors with `s=0`, relaxed refits, lambda grids, paths,
rolling validation, forecasts, and conditional forecast intervals.
EFX adds BigVAR's endogenous-first nested VARX penalty. At each paired lag it
shrinks exogenous coefficients alone and then the combined endogenous and
exogenous block, requiring the exogenous order not to exceed the endogenous
order. EFX does not support contemporaneous predictors. Rolling multi-step
validation is available for VAR and VARX fits with
expanding or fixed windows, L1, L2, and Huber losses, joint lambda-alpha
searches, BigVAR's one-standard-error rule, and retained candidate forecasts.
Dual-grid selection constructs a separate descending lambda path for every
alpha value used by SparseLag, SparseOO, BasicEN, or Tapered penalties.
Lambda matrices use rows for lambda positions and columns for alpha values.
Tapered remains VAR-only. Dual VAR and VARX path types retain
four-dimensional coefficient arrays, while
rolling and leave-one-out validation return both the compatible flattened
candidate vectors and explicit loss, standard-error, and lambda surfaces with
selected row and column indices.
VAR, response-specific VAR, VARX, and dual alpha-specific grid constructors
accept `linear=.true.` to space penalties evenly between the unchanged
zero-model upper bound and `lambda_max/grid_ratio`. Geometric spacing remains
the default.
Rolling out-of-sample reselection updates the penalty choice before every
evaluation forecast using candidate losses whose targets are already observed
at that origin. Ordinary VAR, response-specific VAR, dual-grid VAR, ordinary
VARX, and dual-grid VARX evaluators retain forecasts, component or aggregate
losses, sparsity, convergence status, and complete lambda, alpha, and candidate
index histories. `selection_window` limits the prior validation origins used
for reselection, while `window_size` independently limits the raw estimation
history used for each candidate fit.
EFX fixed fits, paths, lambda bounds, and diagnostics share the VARX interfaces.
Leave-one-out validation removes each selected raw timestamp, reconstructs the
lagged estimation sample, and evaluates the omitted value from its original
history. VAR, VARX, contemporaneous VARX, relaxed refits, joint lambda-alpha
searches, and response-specific penalties are supported. Each fold uses
warm-started paths, and optional first and last omitted observations permit
bounded runs when full LOO is too expensive.
BigVAR's unpenalized benchmark layer fits VAR and VARX models by least
squares, including optional intercept suppression, intercept-only candidates, and direct
multi-step response alignment. Joint endogenous and exogenous lag searches
retain the complete AIC or BIC surface and the selected fit. Rolling benchmark
evaluation reselects the lag orders at every origin and supports direct or
iterated forecasts with L1, L2, and Huber losses.
Unconditional-mean and no-drift random-walk benchmarks return every forecast
and loss together with the mean loss and its standard error. Mean forecasts
use expanding history by default or an optional fixed trailing window; both
benchmarks support arbitrary horizons and L1, L2, or Huber loss.
Accelerated proximal estimation includes centered intercept recovery, fitted
values, residuals, objective and sparsity
diagnostics, automatic zero-model lambda bounds, and warm-started paths.
Structured VAR, response-specific VAR, and VARX fits, paths, validation,
reselection, and lambda grids accept `include_intercept=.false.`. Minnesota
targets already imply a zero intercept, while BGR retains its intercept dummy.
Forecasting and companion-root diagnostics reuse `bigtime_mod` through typed
adapters. The existing `bigtime_mod` implementation remains the shared source
for ordinary lasso and elementwise hierarchical-lag behavior.
BGR implements BigVAR's Banbura-Giannone-Reichlin Bayesian VAR through dummy
observations. It estimates equation-specific prior scales from univariate AR
fits and includes random-walk, sum-of-coefficients, covariance, and intercept
dummies. Fixed tightness fits, optional per-series random-walk indicators, and
the package-default 161-point grid from `sqrt(k*p)` through `5*sqrt(k*p)` are
available through the standard fit, path, and forecast types.
`bigvar_var_to_companion` provides BigVAR's `VarptoVar1MC` conversion while
reusing the shared `bigtime_mod` companion constructor. VAR simulation is
available from supplied innovations or correlated Gaussian innovations.
Random simulation checks stationarity, accepts a full innovation covariance,
and uses BigVAR's 500-observation burn-in by default. The translation is
licensed under GPL-2.0-or-later; see `LICENSE-BIGVAR`.

`bigtime.f90` translates the
[bigtime](https://cran.r-project.org/web/packages/bigtime/index.html) R package.
It provides sparse VAR
estimation under elementwise lasso and hierarchical lag penalties. It includes
the package's centered-data intercept recovery, row-wise accelerated proximal
algorithm, warm-started regularization paths, and automatic geometric lambda
grids whose upper endpoint produces an all-zero coefficient estimate. Sparse
coefficients retain bigtime's lag-major `k x (k*p)` layout.
Sparse VARX estimation jointly updates endogenous and exogenous lag blocks
using the spectral step size of their combined design. Both blocks support L1
or hierarchical lag penalties, separate penalty values and grids, paired
warm-start paths, unequal lag orders, and bigtime's optional squared-Frobenius
shrinkage.
Sparse VARMA estimation implements bigtime's two-stage procedure by fitting a
long sparse VAR, extending and centering its residuals into an innovation
proxy, and passing that proxy to the shared sparse VARX solver as the moving
average regressors. User-supplied innovation proxies remain unchanged. Fixed
fits and Phase II regularization paths retain both stage results for
diagnostics and reuse one Phase I estimate across each path.
Expanding-window one-step cross-validation supports VAR, paired-penalty VARX,
and VARMA Phase II paths. Results include every fold error, mean MSFE, standard
errors, mean sparsity, the minimum-MSFE index, and bigtime's one-standard-error
sparsity choice. VARX and VARMA require both one-standard-error penalties to be
at least as strong as the minimum-MSFE pair. Path information criteria use the
residual covariance determinant and nonzero coefficients for AIC, BIC, and HQ;
selection helpers materialize complete fitted VAR or VARX results from any
chosen path slice.
Post-estimation helpers recursively forecast selected VAR, VARX, and VARMA
models and complete VAR or VARX paths. VARX forecasts accept exogenous values
aligned with the endogenous history and extending through the rows needed by
the forecast horizon. VARMA forecasts retain observed innovation proxies
through the sample end and set future innovations to zero. Companion-matrix
eigenvalues provide stability flags and maximum root moduli, while active-lag
matrices report the largest nonzero lag for every response-predictor pair in
AR, exogenous, or moving-average coefficient blocks.
VAR simulation completes the numerical bigtime translation. Deterministic
coefficient construction accepts supplied Gaussian draws plus arbitrary L1
zero indices or hierarchical trailing-zero counts, applies lag decay, and
contracts the AR block by the upstream `0.99` factor until the requested
companion-root bound is met. Deterministic simulation accepts innovations,
intercepts, a companion-form initial state, and burn-in. Shared-RNG wrappers
generate Gaussian coefficients, random sparse patterns, and Gaussian
innovations through `random_mod`, so `set_random_seed` reproduces
complete coefficient and simulation results. The translation is licensed under
GPL-2.0-or-later; see `LICENSE-BIGTIME`.

`fcvar.f90` translates the
[FCVAR](https://cran.r-project.org/web/packages/FCVAR/index.html) R package. It
provides its level-preserving fast
fractional difference, powers of the fractional lag operator, regression-array
transformation, and concentrated reduced-rank Gaussian estimator at fixed
fractional orders. It supports ranks from zero through full rank and restricted
or unrestricted constants while reusing the shared FFT and linear algebra.
Bounded one- or two-dimensional fractional-order estimation evaluates the
concentrated likelihood on a grid and refines its maximum with the shared
finite-difference BFGS optimizer, retaining boundary grid optima when needed.
Grid results report the global maximum and every strict eight-neighbor local
maximum. Setting `prefer_high_b_local_max` selects the local maximum with the
largest `b`, following FCVAR's identification strategy, and keeps subsequent
refinement in that local basin. Level-parameter grids reuse the preceding
profile estimate as the next simplex start and retain the complete level grid.
The package-default equality restriction `d = b` has a dedicated bounded
one-dimensional likelihood estimator.
General fractional-order estimation accepts affine equalities
`R*[d,b]' = r` and inequalities `C*[d,b]' >= c` together with separate box
bounds and the optional `b <= d` constraint. Two independent equalities are
evaluated directly, one free dimension uses an exact feasible interval and
golden search, and two free dimensions combine feasible-grid initialization,
constrained simplex refinement, and explicit searches of every polygon edge
so boundary optima are retained.
An opt-in level parameter subtracts a fitted vector `mu` before transformation,
profiles its concentrated likelihood with the shared Nelder-Mead optimizer,
and restores the level automatically in forecasts, simulations, and bootstrap
paths. Fixed levels can also be supplied directly. Level coefficients are
included in parameter packing, standard errors, free-parameter counts, and
exact linear restrictions.
Recursive forecasting, supplied-innovation and standard-normal simulation,
and centered-residual Rademacher wild bootstrap reproduce the upstream FCVAR
recursion while using the shared random-number stream for stochastic wrappers.
Ordinary nested-model restrictions have chi-square likelihood-ratio tests using
the shared incomplete-gamma implementation. Rank testing estimates every rank,
reports rank-versus-full likelihood-ratio statistics, and selects ranks by AIC
and BIC. Eligible fractional rank tests use the embedded `fracdist` response
surfaces for p-values through rank difference 12; unsupported specifications
remain explicitly unavailable and can use bootstrap rank inference.
Lag-order selection estimates orders `0:kmax`, applies sequential chi-square
LR tests to added short-run coefficient matrices, and selects an order by AIC
and BIC. It also reports per-series Ljung-Box and heteroskedasticity-robust LM
p-values and a multivariate Ljung-Box p-value for every candidate order.
Bootstrap rank inference estimates nested null and alternative ranks, simulates
under the null with centered-residual Rademacher disturbances, and re-estimates
both ranks for every replication. A supplied-sign pure core supports exact
reproduction, while the stochastic wrapper uses the shared random stream. The
reported bootstrap p-value follows FCVAR's strict-exceedance count divided by
the number of replications.
Postestimation inference packs the complete identified mean-parameter vector,
evaluates a centered finite-difference Hessian of the full concentrated
likelihood, and returns the observed-information covariance and mapped standard
errors. This includes free normalized-`beta` entries in addition to the
parameters differentiated by upstream FCVAR. Characteristic roots are inverse
eigenvalues of the FCVAR block companion matrix and include modulus and
unit-circle diagnostics. The general eigenvalue implementation is shared with
`mar_mod`.
General hypothesis testing supports exact linear restrictions on `(d,b)`,
the level vector, column-major `alpha`, and column-major normalized `beta`
augmented by a restricted constant when present. Restrictions can also be supplied directly
against `fcvar_pack_parameters` output. A pseudoinverse projection and
null-space parameterization enforce the restrictions exactly, including
redundant rows, and LR degrees of freedom use the effective restriction rank.
At fixed fractional orders, coefficient-only `alpha` and `beta` restrictions
use the Boswijk-Doornik switching algorithm: alternating constrained GLS
updates, covariance re-estimation, and Doornik's extrapolation line search.
The high-level restricted estimator routes applicable models through this
solver while mixed packed, order, or level restrictions retain the general
null-space optimizer.
Ordinary tests use the chi-square tail; supplied-sign and shared-RNG wild
bootstrap wrappers simulate under the refreshed restricted fit and report
FCVAR's strict-exceedance bootstrap p-value. The translation is licensed under
GPL-3.0-only; see `LICENSE-FCVAR`.

`mar.f90` translates the
[mAr](https://cran.r-project.org/web/packages/mAr/index.html) R package. It
provides its
augmented and column-scaled QR estimator for multivariate autoregressions,
phase-normalized companion eigenmodes with period and damping diagnostics,
PCA-subspace fitting, and stationary-mean Gaussian VAR simulation. It reuses
the shared linear algebra and random-number modules. The translation is
licensed under GPL-2.0-or-later.

`forecast.f90` translates the
[forecast](https://cran.r-project.org/web/packages/forecast/index.html) R
package. It contains benchmark forecasts, exponential smoothing, Croston and
Theta methods, transforms, Fourier regressors, correlations, accuracy measures,
and the Diebold-Mariano statistic. The translation is licensed under
GPL-3.0-only.

`kfas.f90` translates the
[KFAS](https://cran.r-project.org/web/packages/KFAS/index.html) R package. It
contains typed Gaussian state-space models, Kalman filtering,
smoothing, prediction, exact diffuse filtering for diagonal observation
covariance, standardized innovations, and disturbance means. Covariance
smoothing returns lag-one and forward conditional state covariances using
singular-safe pseudoinverses, while a fast smoother can omit covariance output.
The translation is licensed under GPL-2.0-or-later.

`urca.f90` translates the
[urca](https://cran.r-project.org/web/packages/urca/index.html) R package. It
contains ADF, KPSS, Phillips-Perron, ERS/DF-GLS, Zivot-Andrews, and
Johansen cointegration estimation. Johansen supports deterministic terms,
long-run/transitory specifications, seasonal dummies, and external regressors.
The translation is licensed under GPL-2.0-or-later.

`astsa.f90` translates the
[astsa](https://cran.r-project.org/web/packages/astsa/index.html) R package. It
contains astsa-compatible version-1 `Kfilter` and `Ksmooth` adapters,
polynomial multiplication, ARMA-to-infinite-AR conversion, symmetric matrix
powers, and FDR cutoff behavior. Version-2 filtering and smoothing support
contemporaneously correlated state and observation disturbances using the
shared linear-algebra layer. The translation is licensed under GPL-3.0-only.

The no-input Gaussian state-space EM estimator reuses `astsa_ksmooth`, exposes
RTS gains, computes lag-one smoothed covariances, and estimates the transition,
process covariance, diagonal observation covariance, and initial state moments.

Forward-filter backward-sampling supports caller-provided standard-normal draws
for deterministic or application-controlled RNG streams, plus an intrinsic-RNG
convenience interface using Box-Muller normal generation.

The graphics-independent spectral layer provides theoretical ARMA spectra and
univariate or multivariate periodograms with demeaning or linear detrending,
split tapering, padding, modified-Daniell smoothing, complex spectral matrices,
coherence, phase, bandwidth, and degrees-of-freedom metadata.

SARIMA simulation supports ordinary and seasonal AR, MA, and differencing
terms, astsa-compatible burn-in defaults, caller-supplied innovations for pure
deterministic use, shared-library normal generation, and polynomial-root
causality and invertibility checks.

ARMA diagnostics expose ordinary and seasonal polynomial roots, causality and
invertibility flags, and approximate ordinary or seasonal common-factor checks
without depending on graphics.

Prewhitening supports repeated differencing, Yule-Walker AR order selection by
AIC, aligned filtering of two series, and cross-correlation values with their
lag indices.

The first SARIMA estimation layer evaluates pure conditional Gaussian
likelihoods for supplied ordinary and seasonal ARMA parameters, ordinary and
seasonal differencing, intercept or drift effects, and regression terms. It
returns transformed data, fitted values, innovations, variance, log likelihood,
and AIC, AICc, and BIC.

Conditional SARIMA parameter estimation uses the shared finite-difference BFGS
optimizer with caller-provided starting values, causal and invertible parameter
validation, convergence status, and the final likelihood result. The optimizer
is adapted from the MIT-licensed GARCH-BFGS implementation; its notice is in
`licenses/LICENSE-BFGS`.

The estimator also supports fixed-parameter masks, intercept and drift terms,
regression coefficients, a finite-difference Hessian, parameter covariance and
standard errors, and coefficient-to-standard-error statistics.

Unconstrained optimizer coordinates are mapped through partial
autocorrelations and the Levinson recursion so estimated AR and seasonal AR
blocks remain causal and MA blocks remain invertible. Covariances are returned
on the reported coefficient scale through the transformation Jacobian. The
transform is disabled when dynamic coefficients are fixed individually.

An optional exact Gaussian likelihood converts ordinary or seasonal ARMA
polynomials to an innovations-form companion state-space model and initializes
its covariance from the stationary Lyapunov equation. Ordinary and seasonal
differencing add explicit integration states with exact diffuse initialization
through `kfas_mod`. The innovation scale is profiled analytically after the
diffuse phase.

SARIMA forecasting propagates ordinary and seasonal ARMA recursions, reverses
ordinary and seasonal differencing, continues intercept, drift, and future
regression effects, and computes Gaussian forecast intervals from the full
integrated impulse response.

Graphics-independent SARIMA diagnostics provide standardized effective
residuals, residual autocorrelations, Ljung-Box statistics with fitted-order
degrees-of-freedom corrections, chi-squared tail probabilities, and normal
Q-Q coordinates. Conditional and exact diffuse fits share the same interface.

The LagReg numerical layer estimates a smoothed frequency-domain transfer
function, inverts it to two-sided lag coefficients, selects forward or inverse
lags by threshold, aligns the usable observations, and returns fitted values,
residuals, an intercept, and MSE without plotting.

Signal extraction implements astsa's cosine-tapered finite approximation to an
ideal low-frequency band-pass response. It returns symmetric two-sided filter
coefficients, the filtered series with explicit valid endpoints, and desired
and attained frequency-response arrays.

Stochastic regression implements astsa's multivariate spectral full-versus-
reduced regression analysis. It returns residual power spectra, two-sided
frequency-domain regression coefficients, partial F statistics, squared
coherence, degrees of freedom, and configurable critical values without plots.

Stochastic-volatility maximum likelihood implements astsa's two-component
Gaussian-mixture filter, optional return feedback and leverage correlation,
constrained BFGS estimation, and Hessian covariance and standard errors. Its
filter structure reuses ideas from the MIT-licensed GARCH-BFGS `sv.f90`; the
notice is preserved in `licenses/LICENSE-GARCH-SV`.

The stochastic-volatility MCMC foundation provides a pure conditional particle
filter with ancestor sampling, caller-supplied normal and uniform draws,
normalized particle weights, likelihood estimates, terminal selection, and
genealogy tracing. A convenience wrapper obtains draws from the shared random
module.

The complete stochastic-volatility particle-Gibbs wrapper adds correlated
random-walk Metropolis updates for persistence and state volatility,
inverse-gamma observation-scale updates using centralized gamma generation,
burn-in removal, retained latent trajectories, acceptance rates, and effective
sample-size estimates.

The scalar `ssm` translation estimates state persistence, intercept, process
noise, and observation noise for astsa's linear Gaussian model. It supports
fixed persistence and returns predicted, filtered, and smoothed states and
variances together with Hessian covariance and standard errors.

AR residual bootstrapping reuses the shared Yule-Walker implementation, accepts
caller-supplied residual indices for pure reproducible operation, and provides
a shared-RNG wrapper. Results include simulated series, coefficient draws,
means, and configurable R type-7 quantiles.

Bayesian AR sampling implements astsa's conjugate intercept-plus-lag regression
with Gaussian coefficient updates and inverse-gamma innovation-variance
updates. Pure supplied-draw and shared-RNG interfaces return retained posterior
draws, means, standard deviations, quantiles, and effective sample sizes.

AR spectral information-criterion selection evaluates all Yule-Walker orders
through a configurable maximum, returns relative AIC and BIC tables, and builds
the theoretical spectrum for the selected model with optional linear
detrending.

The frequency-domain linearity test computes block DFTs, normalized bispectral
coefficients, astsa's nonlinearity statistic matrix, an estimated
noncentrality, and noncentral chi-squared p-values without contour plotting.

`itsmr.f90` translates the
[itsmr](https://cran.r-project.org/web/packages/itsmr/index.html) R package. It
provides theoretical ARMA
autocovariances from an innovations-form state representation and the
innovations algorithm for MA estimation. It also provides Hannan-Rissanen ARMA
estimation with ITSMR's preliminary-order rule, regression standard errors,
innovation variance, and corrected AIC. Exact innovations maximum likelihood
uses stable partial-autocorrelation coordinates, BFGS optimization, Hessian
standard errors, and corrected-AIC order selection. Yule-Walker estimation is
centralized in `time_series_stats.f90` for reuse by ASTSA and ITSMR; existing
sample correlations, simulation, smoothing, and periodogram code are not
duplicated. The graphics-independent ARMA forecast returns recovered
innovations, MA-infinity weights, recursive point forecasts, forecast standard
errors, and configurable normal prediction intervals. The translation is
licensed under BSD-2-Clause.

The ARAR implementation performs ITSMR's iterative memory shortening,
exhaustive four-lag sparse AR selection through lag 26, composed-filter
forecasting, and impulse-response prediction intervals. Its result retains the
selected lags, coefficients, memory polynomial, and final filter.

Burg autoregression is implemented in the shared statistics module using
forward-backward errors and reflection coefficients. The ITSMR adapter adds
asymptotic coefficient standard errors and innovations-based variance and AICc.

ARMA models also expose AR-infinity polynomial coefficients and exact residual
output containing fitted values, innovations, time-specific prediction
variances, and standardized innovations for downstream diagnostics.

Shared harmonic regression fits optional polynomial trend terms together with
ITSMR-phase cosine and sine pairs. It returns OLS inference, fitted signal,
residuals, and direct out-of-sample extrapolation.

Typed transformed forecasting composes regular or seasonal differencing,
harmonic regression, polynomial trend, classical seasonal adjustment, and an
optional outer log transform. Difference polynomials are included in forecast
uncertainty before every fitted transform is reversed.

The shared Fourier layer provides normalized direct transforms for real
series. ITSMR spectral rank filtering retains the strongest positive-frequency
bins with conjugate symmetry and returns their indices, frequencies,
amplitudes, coefficients, and real reconstruction.

ITSMR residual diagnostics provide Ljung-Box, McLeod-Li, turning-point,
difference-sign, and rank tests with numerical p-values and no plotting layer.

Compatibility smoothers preserve ITSMR's replicated-endpoint moving average,
recursive exponential initialization, symmetric-index Fourier low-pass filter,
and classical seasonal-component behavior.

### ITSMR coverage

The distinct numerical exports translated in `itsmr.f90` are `aacvf`, `ia`,
`hannan`, `arma`, `autofit`, `forecast`, `burg`, `Resid`, `test`, `arar`,
`smooth.rank`, `smooth.fft`, `smooth.ma`, `smooth.exp`, `season`, `hr`,
`ma.inf`, and `ar.inf`. Their Fortran APIs use descriptive names and derived
results rather than R lists and string vectors.

Sample autocovariance, Yule-Walker estimation, ARMA stability checks,
simulation, periodograms, polynomial trends, and Fourier construction reuse
`forecast_mod`, `time_series_stats_mod`, `astsa_mod`, or the shared Fourier
module. `plotc`, `plota`, and `plots` are graphics-only and intentionally
omitted. `specify` is represented directly by `itsmr_arma_model_t`; package
datasets and `selftest` are not library algorithms.

`arima2.f90` translates the
[arima2](https://cran.r-project.org/web/packages/arima2/index.html) R package. It
builds likelihood inference on the shared ASTSA SARIMA
engine and provides AR/MA polynomial roots, inverse-root coefficient recovery,
Durbin-Levinson coefficient generation, explicit multi-start exact-likelihood
selection with root constraints, AIC/AICc order tables, and fixed-parameter
likelihood profiles. Its random coefficient sampler uses bounded
Durbin-Levinson partial autocorrelations, optional seasonal blocks, root
separation, and intercept draws. The user-facing fitter combines a conditional
SARIMA baseline with these starts and selects the best constrained exact
likelihood. It supports CSS, ML, and CSS-ML selection, regression covariates,
user initial values, and fixed parameter masks. Both Durbin-Levinson and
uniform inverse-root coefficient sampling are available, together with MA root
inversion. The translation is licensed under GPL-3.0-or-later; see
`LICENSE-ARIMA2`.

`fracdiff.f90` translates the
[fracdiff](https://cran.r-project.org/web/packages/fracdiff/index.html) R
package. It provides Jensen-Nielsen
FFT fractional differencing, GPH and tapered Sperio semiparametric memory
estimators, the Haslett-Raftery truncated likelihood filter and stable BFGS
fit with Hessian inference, and fracdiff-compatible simulation. Simulation
reuses `arfima_mod`'s Durbin-Levinson fractional-noise generator, while the FFT
is centralized in `fourier_mod`. Generic ARFIMA covariance,
forecasting, and diagnostics are not duplicated. The translation is licensed
under GPL-2.0-or-later.

`fracdist.f90` translates the
[fracdist](https://cran.r-project.org/web/packages/fracdist/index.html) R
package. It provides response-surface procedures for
fractional unit-root and cointegration-rank distributions. It includes local
quadratic interpolation across the fractional order, statistic-to-chi-square
and chi-square-to-critical-value response surfaces, low-order chi-square
fallback, and pure p-value and critical-value APIs for rank differences 1
through 12. The 24 upstream simulation tables are embedded by the generated
`fracdist_tables.f90`, making the library self-contained. The translation is
licensed under GPL-3.0-only; see
`LICENSE-FRACDIST`.

`nsarfima.f90` translates the
[nsarfima](https://cran.r-project.org/web/packages/nsarfima/index.html) R
package. It adds nonstationary ARFIMA filtering and estimation: zero-padded
causal FFT convolution, Mayoral residual autocorrelations and minimum-distance
estimation, Beran residual pseudo-likelihood, and the package's two integer
integration simulation conventions. Spectral covariance calculations reuse
`arfima_mod`, and FFT, optimization, linear algebra, and random generation use
the shared numerical modules. The translation is licensed under
GPL-3.0-or-later.

`garma.f90` translates the
[garma](https://cran.r-project.org/web/packages/garma/index.html) R package. It
implements the Gegenbauer ARMA numerical layer: Gegenbauer
expansions, raw periodograms and pole selection, Arteche GSP and LPR
semiparametric estimates, CSS, Whittle, and WLL objectives and fits, spectral
inference, long-memory removal, fitted values, Godet-style forecasts, and the
Bartlett periodogram diagnostic. Discontinuous pole estimation uses the shared
derivative-free Nelder-Mead optimizer in `optimization_mod`. The
high-level regression workflow adds intercepts, drift, external regressors,
integer differencing and reintegration, future-regressor forecasts, CSS Hessian
inference, WLL exponent errors, and standard forecast-accuracy measures. The
translation is licensed under GPL-3.0-only.

`esemifar.f90` translates the
[esemifar](https://cran.r-project.org/web/packages/esemifar/index.html) R
package. It implements boundary-aware local-polynomial trend and derivative
smoothing, ESEMIFAR iterative plug-in bandwidth selection, FARIMA order grids,
finite and infinite filter conversions, and analytic, residual-bootstrap, and
FARIMA-refitted predictive-root forecasts. Parametric fitting and simulation
reuse `fracdiff_mod`; random sampling and numerical utilities use the shared
library modules. The translation is licensed under GPL-3.0-only.

`tfarima.f90` translates the
[tfarima](https://cran.r-project.org/web/packages/tfarima/index.html) R package.
It provides the transfer-function ARIMA layer with increasing-lag
polynomial multiplication, division, GCD, powers, rational expansions, and
derivative evaluation. It also provides delayed rational filtering, impulse and
step responses, lag differencing, pulse/step/ramp interventions, seasonal
dummies, harmonic regressors, ARIMA forecast and backcast recursions, theoretical
autocovariances and partial autocorrelations, and Newton autocovariance-to-MA
factorization. ARMA covariance calculations reuse `itsmr_mod`.
Restricted sparse and powered lag polynomials use explicit offset and loading
matrices for linear parameter constraints. Unrestricted and restricted rational
transfer functions can be estimated by a profiled conditional Gaussian
likelihood using a Nelder-Mead warm start, shared BFGS refinement, Hessian
inference, and a caller-supplied ARMA noise filter.
The outlier layer constructs innovation, additive-outlier, level-shift, and
temporary-change responses in the ARIMA residual domain. It supports automatic
or caller-selected timing, repeated residual cleanup, configurable screening
and retention thresholds, and joint least-squares effect and t-ratio refitting.
Monthly calendar support provides Sunday-through-Saturday counts, all six
`CalendarVar` coding forms, reference-day and working-day contrasts,
month-length and leap-year adjustments, Gregorian Easter and Easter-Monday
windows, and direct forecast-horizon extension through the requested number of
monthly observations. Calendar arithmetic reuses `calendar_mod`.
UCARIMA support represents normalized independent ARIMA components, constructs
their common denominator and lifted finite numerators, rejects non-identifiable
shared denominator factors, and obtains the aggregate innovation model through
Cramer-Wold autocovariance factorization. Wiener-Kolmogorov filters expose the
exact rational symmetric numerator and aggregate MA denominator, expand finite
symmetric weights, and decompose complete samples using ARIMA forecast and
backcast endpoint extensions.
Extended Euclidean polynomial GCDs return Bezout coefficients, and a dense
coprime partial-fraction solver supports Hillmer-Tiao conversion of an
aggregate ARIMA spectrum into admissible or canonical UCARIMA components. The
conversion records its Wold-domain fractions and verifies spectral
reconstruction after Cramer-Wold factorization.
UCARIMA components can also be assembled into a block innovations state-space
model. Stationary blocks use iterated Lyapunov covariances, integrated blocks
use KFAS diffuse initialization, and the shared filter and smoother return
component estimates, conditional variances, forecasts, and normal intervals.
Exact UCARIMA estimation uses the same ordinary or diffuse KFAS likelihood.
Component variances use log-scale coordinates, optional masks select AR and MA
coefficients, optional regressors are estimated jointly, and finite-difference
Hessian inference is reported on the natural parameter scale.
Exact transfer-function estimation supports multiple delayed rational inputs,
joint deterministic regressors, free transfer and ARIMA-noise masks, and
ordinary or diffuse KFAS likelihoods. Known future inputs and regressors feed
forecast means, while the fitted ARIMA noise supplies forecast uncertainty and
normal prediction intervals.
Prewhitened cross-correlations expose signed lags, significance flags, and
scaled impulse weights. Transfer identification can select the first
significant nonnegative delay and recover rational numerator and denominator
starting values, while fitted multi-input models provide an input-by-input
residual cross-correlation check.
Backward transfer-model selection uses the exact-fit covariance matrix for
two-sided normal tests, removes the least significant eligible regressor or
complete dynamic input, and refits after every step. Caller keep masks protect
required terms, and the result retains original-variable masks, p-values, and
the ordered removal history.
Transfer diagnostics combine standardized-residual ACF and PACF values,
classical and weighted Ljung-Box tests, Jarque-Bera normality, a cumulative
periodogram white-noise check, and prewhitened residual CCFs for every dynamic
input. Portmanteau calculations are shared with `tsissm` through
`time_series_diagnostics_mod`.
Transfer simulation evaluates complete multi-input and regression signals and
combines them with recursively generated ARIMA noise. Exact fitted models can
use fixed inputs or independently simulated ARIMA inputs, multiple paths,
burn-in, reproducible shared random seeds, caller-supplied innovations, and
presample series and innovation histories. Results retain output, noise,
innovations, inputs, and each deterministic signal component.
Robust Cramer-Wold factorization supports Newton moment solving,
minimum-phase palindromic roots, Bauer initialization, Laurie AS 175, Wilson
iteration, and automatic residual-based selection. UCARIMA construction and
ARIMA-to-UCARIMA conversion use the automatic fallback and report the chosen
method, iterations, convergence, and reconstruction norm. Public conversions
map between palindromic and Wold polynomial coordinates.
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
factors. Covariance-aware state-space forms can also switch between
contemporaneous and one-lag state-disturbance conventions.
The alternative root decomposition groups inverse AR and differencing roots
as trend, seasonal, exponential, or cyclical effects and constructs the
TFARIMA deterministic root basis. Exact Gaussian presample ARMA residuals feed
forecast, backcast, and mixed effect recursions. Results include all four
components, the irregular effect, reconstruction diagnostics, and additive or
log-scale seasonal adjustment. A validated Aberth root fallback handles sparse
seasonal polynomials that are unstable under plain Durand-Kerner iteration.
Direct ARIMA-to-structural-state-space conversion uses the eventual-forecast
root basis and PSI weights to construct either single-source or multiple-source
disturbance forms. Single-source models retain their rank-one joint disturbance
covariance. Multiple-source models match the finite MA autocovariances with
optional variance grouping, report unrestricted admissibility, and use a
nonnegative least-squares fallback when variance matching would otherwise
produce negative disturbances. Both contemporaneous and one-lag forms are
available and round-trip through the state-space-to-ARIMA reduction.
Joint-disturbance structural models can be initialized by filtered GLS and run
through a lag-timed Kalman innovations filter that retains observation/state
disturbance correlation in its gain. The associated information smoother
returns observation-aligned conditional states and covariance matrices.
Filtering reports exact Gaussian likelihoods and standardized innovations;
multi-step forecasts provide latent moments, optional deterministic regression,
and lognormal original-scale moments. Known regressors and log transformations
are supported consistently by initialization, filtering, and smoothing.
The alternative exact reduced-form likelihood applies the aggregate stationary
AR polynomial to the differenced observations and factors the resulting banded
covariance directly. It supports multiple correlated finite disturbance
numerators, profiles the common innovation scale, and returns raw whitened and
determinant-normalized residuals. Public wrappers accept built UCARIMA models or
structural state-space forms, while reusable banded Cholesky and forward-solve
routines expose the underlying numerical operations. The translation is
licensed under GPL-2.0-or-later.

`robustarima.f90` translates the
[robustarima](https://cran.r-project.org/web/packages/robustarima/index.html) R
package. It provides filtered
tau-estimation for regression with ARIMA
errors and the package's bounded redescending rho and psi functions,
MAD-initialized M- and tau-scales, clipped pseudo-observations for robust
correlation diagnostics, and a bounded-innovation ARMA filter. Joint fitting
uses partial-autocorrelation transforms for stationary AR and invertible MA
parameters, supports ordinary and seasonal differencing, regression terms, and
one seasonal MA factor, and can select an AR order by robust AIC. Forecasts,
parameter covariance estimates, and iterative IO/AO/level-shift scoring and
cleaning are included. Regression inference uses the package's filtered-design
tau sandwich covariance and reports its inverse-efficiency correction; the
finite-difference Hessian remains the fallback for a singular estimating
equation and supplies inference for ARMA parameters. The filter propagates the
innovations-form state and
prediction covariance, reports time-varying prediction scales, and implements
the package's rewind-and-accept rule for runs of extreme innovations that mark
level shifts. Outlier effect regressions use those prediction variances, and
forecasts start from the reconstructed robustly cleaned path. The translation
is licensed under BSD-3-Clause; see
`LICENSE-ROBUSTARIMA`.

`bsts.f90` translates the
[bsts](https://cran.r-project.org/web/packages/bsts/index.html) R package. It reuses
the shared KFAS filter and random-number modules to provide reproducible
draw-driven and random-stream Gibbs samplers for local-level and local-linear-
trend models. Posterior outputs retain complete state, observation-variance,
and component-variance draws. Posterior predictive simulation returns draws,
means, standard deviations, and pointwise 95 percent intervals. The package's
static-intercept state has a Gaussian prior, remains constant over time, and
supports posterior observation-variance updates, iteration-specific offsets,
missing responses, and offset-aware forecasts without state innovations. Its
duration-aware sum-to-zero seasonal component and harmonic trigonometric
component are supported with their lower-rank and shared-variance disturbance
structures. Forecast simulation advances the appropriate transition and
innovation-loading schedule. Semilocal linear trends add a random-walk level
and a nonzero-mean AR(1) slope, with posterior draws for the long-run slope,
constrained AR coefficient, and both disturbance variances. The package's
static Gaussian spike-and-slab regression uses collapsed Bernoulli indicator
updates, a full conjugate Gaussian slab, inverse-gamma residual variance,
forced inclusion or exclusion, model-size and model-change limits, structural
offset draws, posterior inclusion probabilities, and regression prediction.
Dynamic-intercept regression supports multiple Gaussian observations at each
ordered time point. It jointly samples a shared local-level intercept, static
spike-and-slab regression, observation and level innovation variances, missing
responses, fitted values and residuals. Grouped forecasts evolve the intercept
once per future time point and share that draw across its observation rows.
Mixed-frequency local-level regression samples latent fine-scale observations
conditioned exactly on observed coarse flow totals through Harvey aggregation
constraints. It supports boundary membership fractions, missing coarse totals,
sparse static regression, retained cumulator paths, and posterior forecasts at
both fine and coarse frequencies. A composite variant adds a local-linear
trend and duration-aware sum-to-zero seasonal block, with separate level,
slope, and seasonal innovation variances and phase-aware structural forecasts.
Gaussian structural and static-regression fits provide posterior one-step
prediction errors, forecast variances, optional standardized innovations,
posterior mean errors, RMSE, and MAE. A local-level holdout routine refits on a
specified prefix before filtering the complete series, matching the package's
out-of-sample cutpoint semantics.
Compatible prediction-error results can be compared over a common full or
restricted interval. The comparison retains cumulative absolute posterior-mean
errors, RMSE and MAE scores, stable tied ranks, and best-model indices while
rejecting inconsistent time dimensions or standardization.
Random-walk dynamic regression uses time-varying coefficient paths sampled by
KFAS forward filtering and backward sampling. It supports independent
inverse-gamma coefficient innovation variances or a predictor-scale-adjusted
hierarchy with a shared gamma rate, structural offset draws, missing responses,
and posterior forecasts that propagate both coefficient and observation noise.
AR(p) dynamic regression adds independent stationary companion processes for
the coefficients, truncated-Gaussian posterior updates for their AR vectors,
predictor-second-moment innovation scaling, and recursive posterior forecasts.
Ordinary Bayesian AR(p) state components use the same stationary companion
machinery for a latent process observed with Gaussian noise, retaining complete
state, AR-parameter, observation-variance, and innovation-variance draws.
Automatic AR components add spike-and-slab lag indicators, geometrically
decreasing default inclusion probabilities and slab scales, optional limits on
indicator changes, stationary coefficient truncation, and posterior lag
inclusion probabilities.
The robust Student-t local-linear trend uses separate normal-gamma mixtures for
level and slope disturbances, conjugate latent-precision and scale updates,
bounded Metropolis degrees-of-freedom draws, optional saved weights, and
Student-t posterior forecasts.
The monthly annual cycle provides an 11-state sum-to-zero seasonal model for
consecutive daily data. Calendar-derived transitions and innovations occur
only on entry to a new month, including leap-year boundaries, and its forecast
routine continues the schedule from the fitted series' final date.
Holiday support defines fixed-date, nth-weekday, last-weekday, irregular
date-range, and standard named US holiday calendars with configurable influence
windows. Random-walk holiday states persist one effect per relative window day,
receive innovations only when that day recurs, and forecast through recurring
or explicitly provided future windows. Fixed holiday regression concatenates
multiple calendars, including unequal-width windows, and samples one constant
coefficient per relative holiday day under a shared normal prior. It supports
structural offset draws, missing responses, inverse-gamma residual variance,
and calendar-aware posterior forecasts. Hierarchical holiday regression pools
three or more equal-width holiday patterns through a learned multivariate mean
and covariance, using multivariate-normal and inverse-Wishart hyperpriors while
retaining holiday-specific coefficient and hyperparameter draws.
Shared local-level models provide multivariate random-walk factors with
independent innovation variances and a learned rectangular loading matrix.
Lower-triangular zero restrictions and a unit diagonal identify the factors;
the sampler supports series-specific observation variances, Gaussian loading
priors, missing observations, structural offset draws, and multivariate
posterior forecasts. Optional spike-and-slab updates select only the free
lower-triangular loadings, with forced inclusion or exclusion, limits on model
changes, retained indicator draws, and posterior inclusion probabilities.
The `mbsts` composition layer augments shared factors with jointly sampled
static regression coefficients for each observed series. Its default
series-specific spike-and-slab update supports forced inclusion or exclusion,
limits on model size and indicator changes, retained indicator draws, and
posterior inclusion probabilities; dense Gaussian regression remains
available. It accepts series-specific predictor arrays, structural offset
draws and missing responses, retains shared and regression contributions,
fitted values and residuals, and produces joint forecasts from future
predictors. Optional series-specific local levels add independent random-walk
states and innovation variances to each response alongside the shared factors;
their posterior states, contributions, variances, and forecast paths are
retained separately. Series-specific local-linear trends add level and slope
states with separate innovation variances. Duration-aware sum-to-zero seasonal
blocks support a common number of seasons and season duration while retaining
independent states and innovation variances for each response.
Non-Gaussian local-level models cover binomial-logit successes and trials and
Poisson counts with observation-specific exposures. Their self-contained
Metropolis-within-Gibbs samplers retain latent log-rate or log-odds states,
fitted means, and innovation variances; supplied-draw and random forecasts
return count-valued posterior predictive paths. Static logit and Poisson
regressions add birth/death spike-and-slab selection, forced predictor states,
model-size and change limits, active-coefficient Metropolis updates, retained
indicators and posterior inclusion probabilities, and predictor-aware count
forecasts.
Logit and Poisson structural models also combine a local-linear trend with a
duration-aware sum-to-zero seasonal state. Exact observation likelihoods are
used in prior-path Metropolis updates, followed by separate conjugate updates
for level, slope, and seasonal innovation variances and structural count
forecasts.
Geometric-sequence and Harvey mixed-frequency cumulator utilities are also
included. Numeric and Gregorian timestamps support duplicate and gap checks,
regular-grid expansion, and observation-to-grid mappings. Integer-labelled
multivariate series can be reshaped between long and wide layouts. General
Harvey aggregation supports vectors and time-by-series matrices, with calendar
helpers for month and quarter boundaries and fractionally apportioned weekly-
to-monthly aggregation. The translation is licensed under MIT; see
`LICENSE-BSTS`.

`bssm.f90` translates the
[bssm](https://cran.r-project.org/web/packages/bssm/index.html) R package. It
provides stable observation
log densities for stochastic-volatility, Poisson, binomial, negative-binomial,
Gamma, and Gaussian models. Its bootstrap particle filter supports fixed or
time-varying observation loadings, transition matrices, transition offsets,
and Gaussian disturbance loadings. A callback interface supports nonlinear
observation densities and state-dependent nonlinear Gaussian transitions.
Filters handle missing observations, stratified resampling, stable likelihood
accumulation, particle ancestry, and predicted and filtered state moments.
Both supplied-random-draw and shared-stream interfaces are available.
Genealogical particle smoothing returns traced state paths and weighted
smoothed moments. An extended Kalman particle filter uses analytic observation
Jacobians to construct particle-specific Gaussian proposals, with exact
observation and transition-to-proposal importance corrections. Durbin-Koopman
simulation smoothing supplies conditional linear-Gaussian state paths,
posterior summaries, missing-observation draws, approximate diffuse
initialization, and optional antithetic shared-stream sampling. A
pseudoinverse RTS correction supports rank-deficient state dynamics. Scalar
and mixed-family non-Gaussian simulation smoothers draw from their Laplace
Gaussian approximations through supplied-normal or antithetic shared-stream
interfaces. Posterior prediction propagates batches of terminal state draws
and returns state, linear-signal, inverse-link mean, and response samples.
Supplied-draw and shared-stream variants cover Gaussian, stochastic-volatility,
Poisson, binomial, negative-binomial, and Gamma observations, time-varying
state-space matrices, exposures, offsets, correlated Gaussian noise, and
pairs each posterior state path with parameter-dependent observation schedules
for in-sample fitted means and response replication. A pure model-update
callback supports mixed families and correlated Gaussian errors. Nonlinear
Gaussian prediction
pairs batches of posterior parameters and terminal states, supports vector
observation and state-dependent noise-loading callbacks, and provides both
future and in-sample supplied-normal or shared-stream response draws. Ordinary
and iterated extended Kalman filters use analytic transition Jacobians,
Joseph-form
covariance updates, convergence diagnostics, and approximate innovation
likelihoods. Multivariate nonlinear Gaussian filtering adds vector observation
Jacobians, general noise-loading covariances, partial component missingness, and
full or mean-only extended Kalman smoothing. Multivariate nonlinear bootstrap
and EKF-proposal particle filters add reproducible supplied-draw and
shared-stream interfaces, ancestry, weighted state summaries, and exact
Gaussian proposal corrections. EKPF proposals can optionally use the iterated
observation update. Multivariate nonlinear bootstrap likelihood estimates also
feed supplied-draw and shared-stream adaptive PMMH parameter samplers. Matching
multivariate EKF-proposal PMMH samplers support ordinary or iterated proposal
updates. Multivariate delayed-acceptance PMMH screens parameter proposals with
the deterministic IEKF likelihood before evaluating a bootstrap particle
likelihood. A matching delayed-acceptance variant uses the multivariate EKPF
for the corrected second-stage likelihood.
The scalar and multivariate unscented Kalman filters provide scaled symmetric
sigma points with configurable `alpha`, `beta`, and `kappa`, nonlinear
observation and transition propagation, partial component missingness, general
Gaussian observation covariances, and innovation likelihoods. Laplace
Gaussian approximations provide distribution-specific pseudo-observations,
conditional modes, scaling corrections, and Gaussian proposal moments while
reusing `kfas_mod` for filtering and smoothing. The corresponding psi auxiliary
particle filter has reproducible supplied-draw and shared-stream interfaces and
returns likelihood estimates, ancestry, and weighted state summaries. SPDK
non-sequential importance sampling draws complete Gaussian-approximation state
trajectories, supports antithetic shared-stream draws, and returns normalized
weights, effective sample size, corrected likelihood, and weighted state
moments. A global iterated extended Kalman smoother linearizes nonlinear
Gaussian observation and transition equations through analytic Jacobian
callbacks. Its nonlinear psi filter corrects both observation and transition
density approximations and provides supplied-draw and shared-stream interfaces.
The vector-observation global approximation retains full observation
covariances, partial component missingness, conditional simulation factors,
and exact-versus-linearized observation and transition corrections. Its psi
filter provides supplied-draw and shared-stream likelihood estimates with
vector observation and nonlinear transition corrections. Adaptive PMMH uses
that estimator for parameter inference with reproducible supplied draws or the
shared random stream. A delayed-acceptance variant screens parameter proposals
with the global Gaussian likelihood before evaluating the psi correction.
Model-aware scalar and vector nonlinear post-correction drivers evaluate psi
likelihoods at approximate-chain states and apply IS1, IS2, or IS3 weights,
with supplied-draw and shared-stream interfaces.
Terminal psi genealogy smoothing retains all `n+1` states. Chain-weighted
helpers combine these smoothers into corrected state moments or resampled
latent trajectories for scalar and vector models.
Model-aware `suggest_N` drivers run replicated scalar or vector nonlinear psi
filters over candidate particle counts and select the first count meeting a
configurable log-likelihood standard-deviation target.
Nonlinear fitted prediction summaries combine MCMC frequency counts and
post-correction weights to return signal and observation means, standard
deviations, and configurable weighted empirical quantiles.
Mixed-family multivariate models support series-specific Poisson, binomial,
negative-binomial, Gamma, and Gaussian observations, time-varying loadings,
partial missingness, and shared linear Gaussian states. Their multivariate
Laplace approximation, bootstrap filter, and psi filter provide supplied-draw
and shared-stream interfaces. Multivariate SPDK importance sampling adds
complete conditional state trajectories, antithetic draws, normalized joint
weights, effective sample size, corrected likelihood, and weighted state
moments. Scalar continuous-discrete SDE models provide Euler-Maruyama and
Milstein substeps, optional positivity reflection, arbitrary observation
log-density callbacks, and supplied-draw and shared-stream bootstrap particle
filters. IS2 state sampling runs a fine filter for each posterior parameter
draw, traces an `n+1` state trajectory through final resampling ancestry, and
returns stable fine-versus-approximate correction weights and ESS diagnostics.
Particle marginal Metropolis-Hastings estimates SDE parameters with
fixed or robust adaptive Gaussian proposals and retains likelihood and
acceptance diagnostics. Delayed-acceptance PMMH screens proposals with a
coarse SDE discretization before applying a fine-discretization correction.
The generic two-estimator delayed-acceptance kernel also supports arbitrary
model-specific likelihood estimators. Direct nonlinear interfaces combine a
Gaussian approximation with bootstrap correction or an EKF likelihood with
EKF-proposal particle correction. Ordinary and iterated extended Kalman
smoothers return nonlinear smoothed state means and covariances, with a
mean-only variant that avoids covariance output.
The reusable supplied-draw PMMH kernel also supports user-defined likelihood
estimators, with direct nonlinear bootstrap, psi, and EKF-proposal particle
filter interfaces. Approximate-likelihood MCMC has direct nonlinear Gaussian
approximation and extended-Kalman adapters. IS1, IS2, and IS3 post-correction
produce normalized weights, effective sample sizes, and corrected parameter
moments. State post-correction combines conditional means and covariances,
transforms them to fitted or predictive linear signals, and resamples weighted
particle trajectories. Replicated log-likelihood diagnostics select a particle
count with standard deviation below a requested threshold.
Sokal adaptive-window and Geyer initial-monotone autocorrelation diagnostics
provide IACT, weighted asymptotic variance, Monte Carlo standard error, and
autocorrelation-adjusted ESS for ordinary, jump-chain, and IS-corrected output.
Particle-filter output is compatible with the common genealogical particle
smoother. Gaussian-approximation proposals reuse the KFAS conditional covariance
smoother, including singular state-covariance support. The translation is
licensed under GPL-2.0-or-later; see
`LICENSE-BSSM`.

`mts.f90` translates the
[MTS](https://cran.r-project.org/web/packages/MTS/index.html) R package. It
provides consecutive and
selected-lag VAR estimation, equation-specific fixed masks, MTS information
criteria, MA-representation matrices, recursive forecasts, innovation
uncertainty, common-sample order selection, orthogonalized and generalized
impulse responses, cumulative responses, and forecast-error variance
decomposition. VARX support includes distributed exogenous lags, fixed masks,
order grids, forecasts from future inputs, and endogenous and exogenous dynamic
responses. The VARMA/VMA analytical layer provides PSI and PI weights, impulse
responses, truncated theoretical covariance, supplied or Gaussian innovation
simulation, and forecasts with error covariance. Conditional Gaussian VMA
estimation supports consecutive or selected lags, fixed masks, Hessian standard
errors, AIC/BIC order selection, and an invertibility-radius diagnostic. Full
conditional Gaussian VARMA estimation adds joint AR/MA optimization,
high-order VAR initialization, fixed parameter masks, selected AR and MA lags,
Hessian inference, and stationarity and invertibility diagnostics.
Multiplicative seasonal VARMA support expands regular and
seasonal matrix polynomials with either MTS multiplication order and estimates
their constrained component matrices by conditional Gaussian likelihood. The
seasonal estimator expects an already stationary input series; differencing is
left to preprocessing. Ordinary and seasonal VARMA refinement iteratively
removes the weakest coefficient below an MTS-style t-ratio threshold, supports
protected coefficients, and records the active-count and removal history.
The VECM layer fits error-correction regressions from supplied cointegration
vectors or vectors selected from `urca_mod`'s Johansen estimator. It supports
equation-specific fixed masks, short-run difference dynamics, inference and
criteria, exact conversion to a level VAR, and level and difference forecasts.
The factor layer extracts standardized or covariance-scale principal components,
reports explained variance and three Bai-Ng factor-count criteria, reconstructs
common and idiosyncratic components, and forecasts the original panel through a
factor VAR with combined factor and residual uncertainty.
The MTS `hfactor` translation estimates weighted least-squares factors under a
caller-supplied loading constraint matrix, returning normalized scores,
coefficient matrix Omega, constrained loadings, residual covariance diagnostics,
explained variation, and forecasts through the shared factor-VAR path.
The conjugate BVAR implementation accepts matrix-normal coefficient and
inverse-Wishart covariance priors, returns posterior coefficient and covariance
moments, and converts posterior means to the shared VAR representation. A
Minnesota-style constructor supplies own-first-lag means, lag-decay shrinkage,
diffuse intercepts, and a data-scaled covariance prior.
The MTS `comVol` translation optionally prewhitens returns by VAR, whitens their
contemporaneous covariance, aggregates lagged dependence among quadratic cross
products, extracts normalized common-volatility directions, and reports exact
component-wise ARCH F tests at caller-selected lag orders.
The MTS `MCHdiag` translation checks an existing conditional covariance path
using radial-residual Ljung-Box and rank tests plus ordinary and robust
multivariate portmanteau tests on squared standardized residuals.
The MTS `BEKK11` translation estimates Gaussian two- or three-dimensional
BEKK(1,1) models with optional means and fixed parameter masks, Hessian
inference, a second-moment persistence guard, conditional covariance and
standardized-residual paths, covariance forecasts, and direct `MCHdiag`
compatibility.
Gaussian Engle DCC(1,1) estimation accepts standardized marginal residuals and
optional marginal variance paths, uses a constrained multinomial-logit
parameterization, returns Q, correlation, residual, and covariance paths with
delta-method inference, forecasts correlations for supplied marginal variance
forecasts, and interoperates directly with `MCHdiag`.
Gaussian and standardized Student-t ADCC extend the shared DCC recursion with
negative-shock dependence, constrained arch, garch, and asymmetry coefficients,
an estimated Student-t degrees of freedom, delta-method inference, covariance
forecasts, and `MCHdiag`-compatible covariance paths.
The MTS-default Tse-Tsui DCC is available with Gaussian and standardized
Student-t likelihoods, configurable rolling correlation windows, constrained
previous and local-correlation weights, degrees-of-freedom estimation,
delta-method inference, held-target covariance forecasts, and MCH diagnostics.
`EWMAvol` support provides fixed or likelihood-estimated decay, Gaussian
covariance paths, inference, and forecasts. `MCholV` combines optional VAR
prewhitening, expanding-window recursive triangular regressions, EWMA-smoothed
Cholesky coefficients, constrained component GARCH(1,1) fits, positive-definite
covariance reconstruction, forecasts, and MCH diagnostics.
`SCCor` returns unconstrained and exact group-averaged correlation matrices over
a selected window. Public `archTest` and `MarchTest` equivalents provide the
MTS McLeod-Li, rank-based, radial, multivariate, and robust multivariate ARCH
statistics with degrees of freedom and p-values.
The `mtCopula` layer converts correlations to and from hyperspherical angles,
builds positive-definite group-equicorrelation matrices, estimates dynamic
baseline/local angle recursions with a standardized Student-t copula likelihood,
supports fixed or estimated baseline angles, and returns Hessian inference and
complete correlation paths.
The MTS `Vmiss` and `Vpmiss` translations estimate fully or partially missing
vector observations by conditional GLS over the current and future PI-weight
innovation equations, returning conditional covariance and a completed series.
Granger causality is available as a joint VAR coefficient Wald test for selected
target equations. `mq` and `MTSdiag` equivalents provide cumulative
multivariate Ljung-Box statistics, residual cross-correlation paths, and
significance probabilities. Rolling-origin VAR backtesting records forecasts
and errors with horizon-specific RMSE and mean absolute error.
Scalar-component translation has begun with first-stage `SCMid`: canonical
correlations over the AR/MA order grid, serial-correlation variance correction,
sequential zero-correlation tests, and diagonal-difference identification tables.
The translation is licensed under Artistic-2.0; see
`LICENSE-MTS`.

## Build

```sh
cmake -S . -B build -G Ninja -DCMAKE_Fortran_COMPILER=gfortran
cmake --build build
ctest --test-dir build --output-on-failure
```

Run the VAR(2) simulation and fitting example with:

```sh
./build/example_var2_fit
```

The translated source packages use GPL licenses. Redistribution of derived code
must preserve the applicable license terms.

The shared date implementation is adapted from the MIT-licensed DataFrame
project. Its notice is preserved in `licenses/LICENSE-DATE`.
