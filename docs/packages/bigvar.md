# BigVAR

[Back to the implemented-package index](../../README.md#implemented-packages).

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
