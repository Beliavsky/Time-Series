# astsa

[Back to the implemented-package index](../../README.md#implemented-packages).

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
