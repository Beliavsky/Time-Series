# vars

[Back to the implemented-package index](../../README.md#implemented-packages).

`vars.f90` translates the
[vars](https://cran.r-project.org/web/packages/vars/index.html) R package. It
provides
reduced-form VAR estimation with constant, trend,
centered seasonal, and exogenous terms; common-sample AIC, HQ, SC, and FPE lag
selection; manual and sequential t-ratio restrictions; MA coefficients,
orthogonalized responses, and companion roots; asymptotic and adjusted
portmanteau, Breusch-Godfrey, Edgerton-Shukur, multivariate Jarque-Bera, and
multivariate and equationwise ARCH tests; instantaneous-causality statistics;
and OLS-CUSUM fluctuation paths. Structural routines include restricted A, B,
and AB estimation by direct Gaussian likelihood or the package's Fisher-scoring
recursion, parameter standard errors, Blanchard-Quah identification, and SVEC
short- and long-run zero restrictions built from `urca_mod` Johansen results.
The Johansen-to-level-VAR conversion, reduced-form IRF intervals, structural
impact intervals, SVEC bootstrap standard errors, and bootstrap Granger test
are also included. SVAR and SVEC structural impulse responses and forecast-error
variance decompositions share one level-VAR response kernel. Principal entry
points are re-exported by `multivariate_mod` and `diagnostics_mod`. The
translation is licensed under GPL-2.0-or-later; see
`LICENSE-VARS`.
