# bsts

[Back to the implemented-package index](../../README.md#implemented-packages).

## Overview

`bsts.f90` translates the
[bsts](https://cran.r-project.org/web/packages/bsts/index.html) R package.

## Algorithms and Procedures

It reuses the shared KFAS filter and random-number modules to provide
reproducible draw-driven and random-stream Gibbs samplers for local-level and
local-linear-trend models.

Posterior outputs retain complete state, observation-variance, and
component-variance draws.

Posterior predictive simulation returns draws, means, standard deviations, and
pointwise 95 percent intervals.

The package's static-intercept state has a Gaussian prior, remains constant
over time, and supports posterior observation-variance updates,
iteration-specific offsets, missing responses, and offset-aware forecasts
without state innovations. Its duration-aware sum-to-zero seasonal component
and harmonic trigonometric component are supported with their lower-rank and
shared-variance disturbance structures.

Forecast simulation advances the appropriate transition and innovation-loading
schedule.

Semilocal linear trends add a random-walk level and a nonzero-mean AR(1) slope,
with posterior draws for the long-run slope, constrained AR coefficient, and
both disturbance variances.

The package's static Gaussian spike-and-slab regression uses collapsed
Bernoulli indicator updates, a full conjugate Gaussian slab, inverse-gamma
residual variance, forced inclusion or exclusion, model-size and model-change
limits, structural offset draws, posterior inclusion probabilities, and
regression prediction.

Dynamic-intercept regression supports multiple Gaussian observations at each
ordered time point. It jointly samples a shared local-level intercept, static
spike-and-slab regression, observation and level innovation variances, missing
responses, fitted values and residuals.

Grouped forecasts evolve the intercept once per future time point and share
that draw across its observation rows.

Mixed-frequency local-level regression samples latent fine-scale observations
conditioned exactly on observed coarse flow totals through Harvey aggregation
constraints. It supports boundary membership fractions, missing coarse totals,
sparse static regression, retained cumulator paths, and posterior forecasts at
both fine and coarse frequencies.

A composite variant adds a local-linear trend and duration-aware sum-to-zero
seasonal block, with separate level, slope, and seasonal innovation variances
and phase-aware structural forecasts.

Gaussian structural and static-regression fits provide posterior one-step
prediction errors, forecast variances, optional standardized innovations,
posterior mean errors, RMSE, and MAE.

A local-level holdout routine refits on a specified prefix before filtering the
complete series, matching the package's out-of-sample cutpoint semantics.

Compatible prediction-error results can be compared over a common full or
restricted interval.

The comparison retains cumulative absolute posterior-mean errors, RMSE and MAE
scores, stable tied ranks, and best-model indices while rejecting inconsistent
time dimensions or standardization.

Random-walk dynamic regression uses time-varying coefficient paths sampled by
KFAS forward filtering and backward sampling. It supports independent
inverse-gamma coefficient innovation variances or a predictor-scale-adjusted
hierarchy with a shared gamma rate, structural offset draws, missing responses,
and posterior forecasts that propagate both coefficient and observation noise.

AR(p) dynamic regression adds independent stationary companion processes for
the coefficients, truncated-Gaussian posterior updates for their AR vectors,
predictor-second-moment innovation scaling, and recursive posterior forecasts.

Ordinary Bayesian AR(p) state components use the same stationary companion
machinery for a latent process observed with Gaussian noise, retaining complete
state, AR-parameter, observation-variance, and innovation-variance draws.

Automatic AR components add spike-and-slab lag indicators, geometrically
decreasing default inclusion probabilities and slab scales, optional limits on
indicator changes, stationary coefficient truncation, and posterior lag
inclusion probabilities.

The robust Student-t local-linear trend uses separate normal-gamma mixtures for
level and slope disturbances, conjugate latent-precision and scale updates,
bounded Metropolis degrees-of-freedom draws, optional saved weights, and
Student-t posterior forecasts.

The monthly annual cycle provides an 11-state sum-to-zero seasonal model for
consecutive daily data.

Calendar-derived transitions and innovations occur only on entry to a new
month, including leap-year boundaries, and its forecast routine continues the
schedule from the fitted series' final date.

Holiday support defines fixed-date, nth-weekday, last-weekday, irregular
date-range, and standard named US holiday calendars with configurable influence
windows.

Random-walk holiday states persist one effect per relative window day, receive
innovations only when that day recurs, and forecast through recurring or
explicitly provided future windows.

Fixed holiday regression concatenates multiple calendars, including
unequal-width windows, and samples one constant coefficient per relative
holiday day under a shared normal prior. It supports structural offset draws,
missing responses, inverse-gamma residual variance, and calendar-aware
posterior forecasts.

Hierarchical holiday regression pools three or more equal-width holiday
patterns through a learned multivariate mean and covariance, using
multivariate-normal and inverse-Wishart hyperpriors while retaining
holiday-specific coefficient and hyperparameter draws.

Shared local-level models provide multivariate random-walk factors with
independent innovation variances and a learned rectangular loading matrix.

Lower-triangular zero restrictions and a unit diagonal identify the factors;
the sampler supports series-specific observation variances, Gaussian loading
priors, missing observations, structural offset draws, and multivariate
posterior forecasts.

Optional spike-and-slab updates select only the free lower-triangular loadings,
with forced inclusion or exclusion, limits on model changes, retained indicator
draws, and posterior inclusion probabilities.

The `mbsts` composition layer augments shared factors with jointly sampled
static regression coefficients for each observed series. Its default
series-specific spike-and-slab update supports forced inclusion or exclusion,
limits on model size and indicator changes, retained indicator draws, and
posterior inclusion probabilities; dense Gaussian regression remains available.
It accepts series-specific predictor arrays, structural offset draws and
missing responses, retains shared and regression contributions, fitted values
and residuals, and produces joint forecasts from future predictors.

Optional series-specific local levels add independent random-walk states and
innovation variances to each response alongside the shared factors; their
posterior states, contributions, variances, and forecast paths are retained
separately.

Series-specific local-linear trends add level and slope states with separate
innovation variances.

Duration-aware sum-to-zero seasonal blocks support a common number of seasons
and season duration while retaining independent states and innovation variances
for each response.

Non-Gaussian local-level models cover binomial-logit successes and trials and
Poisson counts with observation-specific exposures. Their self-contained
Metropolis-within-Gibbs samplers retain latent log-rate or log-odds states,
fitted means, and innovation variances; supplied-draw and random forecasts
return count-valued posterior predictive paths.

Static logit and Poisson regressions add birth/death spike-and-slab selection,
forced predictor states, model-size and change limits, active-coefficient
Metropolis updates, retained indicators and posterior inclusion probabilities,
and predictor-aware count forecasts.

Logit and Poisson structural models also combine a local-linear trend with a
duration-aware sum-to-zero seasonal state.

Exact observation likelihoods are used in prior-path Metropolis updates,
followed by separate conjugate updates for level, slope, and seasonal
innovation variances and structural count forecasts.

Geometric-sequence and Harvey mixed-frequency cumulator utilities are also
included.

Numeric and Gregorian timestamps support duplicate and gap checks, regular-grid
expansion, and observation-to-grid mappings.

Integer-labelled multivariate series can be reshaped between long and wide
layouts.

General Harvey aggregation supports vectors and time-by-series matrices, with
calendar helpers for month and quarter boundaries and fractionally apportioned
weekly-to-monthly aggregation.

## Licensing

The translation is licensed under MIT; see `LICENSE-BSTS`.
