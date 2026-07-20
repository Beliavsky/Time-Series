# gmvarkit

[Back to the implemented-package index](../../README.md#implemented-packages).

## Overview

`gmvarkit.f90` translates the
[gmvarkit](https://cran.r-project.org/web/packages/gmvarkit/index.html) R
package.

## Algorithms and Procedures

It provides stationary regime means and lag-vector covariances, endogenous
Gaussian and Student-t regime weights, conditional regime and mixture moments,
ARCH scaling for Student-t regimes, and conditional or exact mixture log
likelihoods for arbitrary VAR order.

Recursive simulation retains sampled regimes and future mixing weights for each
path.

Monte Carlo forecasting provides means, medians, empirical type-7 quantiles,
and expected future regime weights.

`gmvarkit_estimate` performs local finite-difference BFGS likelihood refinement
with Cholesky covariance, softmax weight, and transformed degrees-of-freedom
parameters.

Invalid and nonstationary trial models are penalized.

`gmvarkit_genetic_estimate` supplies the package's global-search layer with
elitist selection, crossover, decreasing transformed-space mutation, a recorded
best-objective path, and final BFGS refinement.

`gmvarkit_estimate_constrained` supports general linear maps for stacked AR
coefficients, shared unconditional means across regime groups, and fixed
mixture weights while optimizing only the remaining free coordinates.

Hessian inversion supplies transformed-scale covariance matrices and standard
errors.

General linear Wald tests and nested-model likelihood-ratio tests use
chi-square reference distributions.

The inference vector stores each regime's intercept, column-major AR entries,
and rowwise lower-Cholesky entries, followed by weight logits and
log-transformed Student-t degrees of freedom; diagonal Cholesky entries are
stored on log scale.

Observationwise finite-difference scores provide OPG information and Rao score
tests.

Sequential multivariate quantile residuals support Gaussian and Student-t
regimes; moment diagnostics test normality, serial correlation, and dependence
in squared residuals.

Supplying the fitted model and data applies the Kalliovirta-Saikkonen
parameter-estimation correction using residual-moment derivatives, OPG Fisher
information, and moment-score cross covariances; otherwise the tests use
empirical moment covariance matrices.

`gmvarkit_girf` estimates nonlinear generalized impulse responses by paired
simulation with common regimes and innovations, lower-Cholesky identification,
and structural-shock replacement. It also returns responses of the endogenous
regime probabilities.

`gmvarkit_gfevd` forms horizon-specific decompositions from normalized
cumulative squared mean GIRFs.

`gmvarkit_linear_irf` computes deterministic regime-specific VAR responses
under lower-Cholesky recursive identification.

Selected variables can be accumulated over horizons, and selected shocks can be
normalized to a requested impact response of a specified variable.

`gmvarkit_unconditional_moments` combines regime means and autocovariances into
the mixture mean, lag-zero-through-p autocovariances, and autocorrelations.

`gmvarkit_pearson_residuals` returns raw residuals or applies the inverse
symmetric square root of each conditional covariance.

The elemental `gmvarkit_information_criteria` evaluates AIC, HQIC, and BIC from
a likelihood, free-parameter count, and effective sample size.

Structural covariance support represents each regime as `Omega_m = W Lambda_m
W'`.

`gmvarkit_identify_structural` recovers a common impact matrix and relative
variances through covariance heteroskedasticity and checks simultaneous
diagonalization.

Companion procedures construct implied covariances, reorder or sign-reverse
shocks, and change the unit-variance reference regime.

Linear and generalized impulse responses accept this identification as an
optional argument.

`gmvarkit_linear_irf_bootstrap` provides fixed-design Rademacher wild-bootstrap
confidence intervals when regime dynamics are linear.

Starting observations remain fixed, fitted residual vectors receive common
rowwise signs, and mixture weights are held fixed during re-estimation.

Structural bootstrap solutions are matched to the baseline lambda profiles and
impact-column signs before their responses are summarized.

`gmvarkit_girf_inference` adds the outer initial-state distribution used by the
package's full GIRF interface. It supports every observed length-p history,
fixed histories, or stationary Gaussian and Student histories drawn from
selected regimes; it returns individual responses, point estimates, and
equal-tail bounds after optional accumulation and instantaneous or peak
scaling.

`gmvarkit_gfevd_inference` calculates a decomposition for every history and
averages them, including decompositions of mixing-weight responses.

`gmvarkit_estimate_structural` optimizes the common impact matrix and positive
relative regime variances directly in the likelihood.

Structural restrictions support exact fixed and zero entries of `W`, positive
or negative sign constraints on free entries, fully fixed lambda matrices, and
nonnegative linear lambda mappings from positive free coordinates. The result
includes the compatible reduced-form model and optional transformed-coordinate
Hessian covariance, gradient, and standard errors.

`gmvarkit_profile_likelihood` evaluates fixed-nuisance likelihood slices for
selected estimation coordinates.

Reduced-form and structural overloads use the same transformed
parameterizations as their estimators and return grid values, log likelihoods,
and a validity mask without imposing a graphics dependency.

`gmvarkit_convert_student_regimes` replaces Student regimes above a selected
degrees-of-freedom threshold by Gaussian regimes, orders each regime family by
mixing weight, and returns both directions of the resulting permutation.

Structural variance columns and their reference regime follow the same mapping,
and the converted reduced-form or structural model can optionally be re-fitted.

`gmvarkit_companion_eigenvalues` returns the complex companion roots, moduli,
spectral radii, stationarity flags, and tolerance-based near-unit-root flags
for each regime.

`gmvarkit_covariance_eigenvalues` reports ordinary covariance eigenvalues and
near-singularity flags, together with Cholesky-whitened generalized eigenvalues
for every covariance pair.

Relative eigenvalue separations flag weak heteroskedastic identification
without printing warnings.

`gmvarkit_multistart_estimate` optimizes an array of starting models, retains
all local fits and convergence information, and ranks successful solutions by
log likelihood.

Regimes are ordered canonically within Gaussian and Student families before
parameter comparisons, allowing label-switched and numerically equivalent fits
to be linked through `duplicate_of` rather than counted as distinct local
maxima.

`gmvarkit_structural_multistart_estimate` applies common impact and relative
variance restrictions to paired reduced-form and structural starts. It ranks
fits by structural likelihood, compares their canonical implied reduced forms
to remove regime, shock-order, and shock-sign aliases, and calculates numerical
Hessian inference only for the selected best solution.

`gmvarkit_location_parameters` returns both intercepts and their equivalent
stationary regime means.

`gmvarkit_model_from_regime_means` performs the inverse mapping and returns an
ordinary intercept-parameterized model, providing the package's reversible
mean/intercept workflow without introducing ambiguous state into the likelihood
model type.

## Licensing

The translation is licensed under GPL-3.0-only.
