# nsarfima

[Back to the implemented-package index](../../README.md#implemented-packages).

## Overview

`nsarfima.f90` translates the
[nsarfima](https://cran.r-project.org/web/packages/nsarfima/index.html) R
package.

## Algorithms and Procedures

It adds nonstationary ARFIMA filtering and estimation: zero-padded causal FFT
convolution, Mayoral residual autocorrelations and minimum-distance estimation,
Beran residual pseudo-likelihood, and the package's two integer integration
simulation conventions.

Spectral covariance calculations reuse `arfima_mod`, and FFT, optimization,
linear algebra, and random generation use the shared numerical modules.

## Licensing

The translation is licensed under GPL-3.0-or-later.
