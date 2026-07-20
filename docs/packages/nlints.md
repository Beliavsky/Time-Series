# NlinTS

[Back to the implemented-package index](../../README.md#implemented-packages).

## Overview

`nlints.f90` translates the
[NlinTS](https://cran.r-project.org/web/packages/NlinTS/index.html) R package.

## Algorithms and Procedures

It provides discrete Shannon and joint entropy, bivariate and multivariate
mutual information, and lagged transfer entropy with selectable logarithm bases
and package-compatible normalization.

Continuous information measures include Kozachenko maximum-norm entropy, both
NlinTS KSG mutual-information variants, and KSG transfer entropy.

The multivariate neural autoregression retains variable-major lag construction
and columnwise min-max scaling, supports arbitrary hidden-layer sizes, optional
biases, linear, sigmoid, ReLU, and tanh activations, and mini-batch SGD or
bias-corrected Adam training.

Models retain aligned fitted values, residuals, normalized residual sums of
squares, retained scales, incremental training, prediction, next-step
forecasts, and package-shaped rolling forecast tables whose final row is out of
sample.

Versioned text persistence retains the complete architecture, scaling state,
weights, and SGD or Adam optimizer state.

Classical and neural Granger tests report the causality index, F statistic,
degrees of freedom, and beta-ratio p-value.

Classical tests may infer separate integration orders by repeated trend ADF
decisions, difference each series accordingly, and align them to the common
maximum-order sample; the package ADF entry point reuses `urca_mod`.

Principal interfaces are available through `diagnostics_mod`,
`forecasting_mod`, and `multivariate_mod`.

## Licensing

The translation is licensed under GPL-2.0-or-later; see `LICENSE-NLINTS`.
