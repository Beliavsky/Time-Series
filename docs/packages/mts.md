# MTS

[Back to the implemented-package index](../../README.md#implemented-packages).

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
