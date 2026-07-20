# bigtime

[Back to the implemented-package index](../../README.md#implemented-packages).

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
