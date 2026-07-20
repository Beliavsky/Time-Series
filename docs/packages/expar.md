# EXPAR

[Back to the implemented-package index](../../README.md#implemented-packages).

`expar.f90` translates the
[EXPAR](https://cran.r-project.org/web/packages/EXPAR/index.html) R package.
It evaluates amplitude-dependent exponential autoregressions, estimates their
phi, pi, and exponential-scale parameters by conditional RSS minimization with
finite-difference BFGS, constructs AR-based starting values, selects the order
by AIC, corrected AIC, or BIC, and computes recursive point forecasts. The
principal estimation interfaces are also available through `arma_mod`, and
forecasting is available through `forecasting_mod`. The translation is
licensed under GPL-3.0-only; see `LICENSE-EXPAR`.
