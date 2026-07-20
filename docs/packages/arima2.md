# arima2

[Back to the implemented-package index](../../README.md#implemented-packages).

`arima2.f90` translates the
[arima2](https://cran.r-project.org/web/packages/arima2/index.html) R package. It
builds likelihood inference on the shared ASTSA SARIMA
engine and provides AR/MA polynomial roots, inverse-root coefficient recovery,
Durbin-Levinson coefficient generation, explicit multi-start exact-likelihood
selection with root constraints, AIC/AICc order tables, and fixed-parameter
likelihood profiles. Its random coefficient sampler uses bounded
Durbin-Levinson partial autocorrelations, optional seasonal blocks, root
separation, and intercept draws. The user-facing fitter combines a conditional
SARIMA baseline with these starts and selects the best constrained exact
likelihood. It supports CSS, ML, and CSS-ML selection, regression covariates,
user initial values, and fixed parameter masks. Both Durbin-Levinson and
uniform inverse-root coefficient sampling are available, together with MA root
inversion. The translation is licensed under GPL-3.0-or-later; see
`LICENSE-ARIMA2`.
