# bvars

[Back to the implemented-package index](../../README.md#implemented-packages).

## Overview

`bvars.f90` translates the
[bvars](https://cran.r-project.org/web/packages/bvars/index.html) R package.

## Algorithms and Procedures

It provides lag-major BVAR data preparation, the package's default random-walk
or stationary matrix-normal/inverse-Wishart prior, matrix-normal simulation,
and the complete homoskedastic Gaussian posterior sampler.

Posterior residuals, fitted-density draws, recursive forecasts with
contemporaneous exogenous inputs and common variance multipliers, and
orthogonalized FEVDs are included.

FEVD computation reuses `bvartools_mod`.

The homoskedastic multivariate Student-t extension samples observation-specific
inverse-gamma variance multipliers and degrees of freedom with an adaptive
transformed Metropolis step; matching future scale simulation feeds directly
into the forecast routine.

Centred and non-centred common stochastic volatility implement the package's
mixture-indicator, persistence, innovation-variance, hierarchical variance, and
latent-state updates.

The non-centred sampler uses ASIS interweaving and generalized-inverse-Gaussian
slice draws. Both accept the dimension-specific auxiliary mixture table
explicitly, retain every state and hyperparameter draw, and simulate future
common-variance paths for forecasts.

The Student-t scale-mixture and degrees-of-freedom updates can be combined with
either stochastic-volatility parameterization in a single Gibbs sampler.

Forecasts can condition individual variables and horizons on supplied values;
use NaN for entries left free.

`bvars_sv_auxiliary_mixture` constructs the required auxiliary mixture table
for any number of variables with a reusable deterministic Gaussian-mixture EM
fit; its simulation count and convergence settings are configurable.

## Licensing

The translation is licensed under GPL-3.0-or-later; see `LICENSE-BVARS`.
