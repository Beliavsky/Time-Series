# bssm

[Back to the implemented-package index](../../README.md#implemented-packages).

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
