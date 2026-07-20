# GMDH

[Back to the implemented-package index](../../README.md#implemented-packages).

`gmdh.f90` translates the
[GMDH](https://cran.r-project.org/web/packages/GMDH/index.html) R package. It
fits classical group-method-of-data-handling networks by generating pairwise
quadratic neurons, selecting polynomial, sigmoid, radial-basis, or tangent
response transformations, choosing unpenalized-intercept ridge coefficients
by a sequential validation split, and retaining the best neurons at each
self-organized layer. The revised GMDH algorithm adds increasing-prefix linear
feedback neurons to every candidate layer. Separate fitting and recursive
prediction interfaces accompany the package-compatible short-term forecast,
which estimates horizon-specific standard errors from rolling forecast errors
and constructs normal prediction intervals for as many as five steps.
Fitted-value, residual, coefficient, and display interfaces are included, and
the principal routines are also available through `forecasting_mod`. The
translation is licensed under GPL-2.0-or-later; see `LICENSE-GMDH`.
