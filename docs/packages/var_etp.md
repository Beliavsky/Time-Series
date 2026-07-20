# VAR.etp

[Back to the implemented-package index](../../README.md#implemented-packages).

`var_etp.f90` translates the
[VAR.etp](https://cran.r-project.org/web/packages/VAR.etp/index.html) R package.
It provides
Nicholls-Pope asymptotic and residual-bootstrap VAR bias
corrections with Kilian stationarity adjustment; analytic forecast MSE including
coefficient uncertainty; forward/backward bootstrap and bootstrap-after-bootstrap
prediction intervals; constrained system estimates; Wald and nested likelihood-
ratio tests with iid or Mammen wild bootstrap inference; and Kim's improved
augmented predictive regression with Shaman-Stine predictor correction, joint
tests, order selection, and dynamic forecasts. Ordinary VAR estimation, roots,
MA coefficients, and impulse responses reuse `vars_mod`. The translation is
licensed under GPL-2.0-only; see `LICENSE-VAR-ETP`.
