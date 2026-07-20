# FCVAR

[Back to the implemented-package index](../../README.md#implemented-packages).

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
