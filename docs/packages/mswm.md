# MSwM

[Back to the implemented-package index](../../README.md#implemented-packages).

## Overview

`mswm.f90` translates the
[MSwM](https://cran.r-project.org/web/packages/MSwM/index.html) R package.

## Algorithms and Procedures

The Gaussian linear-model layer provides scaled Hamilton filtering, Kim
smoothing, expected transition probabilities, and EM estimation for arbitrary
regime counts.

Regression coefficients may be shared or regime-specific, innovation variance
may be common or switching, and callers may provide initial coefficients,
scales, and transition probabilities.

Generalized switching regression supports Poisson with a log link,
Bernoulli-binomial with a logit link, and Gamma models with log or inverse
links, using weighted IRLS within each EM M-step.

Results include conditional means, residuals, filtered and smoothed state
probabilities, decoded states, transition estimates, log likelihood, AIC, and
BIC.

Gaussian and generalized fits provide finite-difference Hessian covariance
estimates, natural-scale standard errors, and approximate confidence intervals.

Diagnostics include raw Gaussian residuals and family-standardized generalized
residuals, returned by regime or combined using smoothed state probabilities.

Gaussian fitting also supports selection over supplied or reproducible random
EM starts, with equivalent multistart fitting for generalized models.

Autoregressive wrappers construct response lags alongside optional exogenous
predictors and retain the lag order and terminal response history for
subsequent forecasting.

Gaussian and generalized results share overloaded residual, state-decoding, and
concise display interfaces; observation arrays are omitted from displays unless
explicitly requested.

The principal entry points are also available through `markov_switching_mod`.

## Licensing

The translation is licensed under GPL-2.0-or-later; see `LICENSE-MSWM`.
