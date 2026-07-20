# BVAR

[Back to the implemented-package index](../../README.md#implemented-packages).

## Overview

`bvar.f90` translates the
[BVAR](https://cran.r-project.org/web/packages/BVAR/index.html) R package.

## Algorithms and Procedures

It provides the package's diagonal Minnesota row variances, sum-of-coefficients
and single-unit-root dummy observations, conjugate posterior sufficient
statistics, and closed-form log marginal likelihood.

Hierarchical lambda, alpha, innovation-scale, SOC, and SUR tightness estimation
uses bounded random-walk Metropolis sampling with Gamma and inverse-Gamma
hyperpriors, optional burn-in acceptance tuning, a supplied-randomness
reproducible core, and a shared-RNG wrapper.

## Licensing

The translation is licensed under GPL-3.0-or-later; see `LICENSE-BVAR`.
