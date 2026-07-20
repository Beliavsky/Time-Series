# bvartools

[Back to the implemented-package index](../../README.md#implemented-packages).

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
