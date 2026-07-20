# mixAR

[Back to the implemented-package index](../../README.md#implemented-packages).

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
