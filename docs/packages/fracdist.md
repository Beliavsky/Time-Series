# fracdist

[Back to the implemented-package index](../../README.md#implemented-packages).

## Overview

`fracdist.f90` translates the
[fracdist](https://cran.r-project.org/web/packages/fracdist/index.html) R
package.

## Algorithms and Procedures

It provides response-surface procedures for fractional unit-root and
cointegration-rank distributions. It includes local quadratic interpolation
across the fractional order, statistic-to-chi-square and
chi-square-to-critical-value response surfaces, low-order chi-square fallback,
and pure p-value and critical-value APIs for rank differences 1 through 12.

The 24 upstream simulation tables are embedded by the generated
`fracdist_tables.f90`, making the library self-contained.

## Licensing

The translation is licensed under GPL-3.0-only; see `LICENSE-FRACDIST`.
