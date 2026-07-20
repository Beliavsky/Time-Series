# nnfor

[Back to the implemented-package index](../../README.md#implemented-packages).

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
