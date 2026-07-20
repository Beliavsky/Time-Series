# narfima

[Back to the implemented-package index](../../README.md#implemented-packages).

`narfima.f90` translates the
[narfima](https://cran.r-project.org/web/packages/narfima/index.html) R package.
Its common neural autoregressive model combines consecutive response lags,
ordinary and seasonal baseline-innovation lags, and optional contemporaneous
regressors. Response, innovation, and regressor transformations retain their
scales for forecasting, and transformations may be disabled. Repeated
one-hidden-layer networks are combined by their mean, with optional direct
linear connections optimized jointly with the nonlinear weights. Automatic
lag selection uses seasonally adjusted Yule-Walker fits. Package-facing
constructors obtain regression-aware baseline innovations from the library's
ARFIMA, ARIMA, jointly backfitted Bayesian structural time-series, or naive
implementations without duplicating the neural fitting layer. Missing baseline
inputs are interpolated, while the common neural model retains complete-case
alignment for lagged predictors.
Recursive point forecasts use zero future innovations; supplied innovations,
Gaussian simulation, and residual bootstrap produce complete forecast paths,
sample standard deviations, and type-8 empirical prediction intervals.
Principal interfaces are available through `long_memory_mod` and
`forecasting_mod`. The translation is licensed under GPL-3.0-only; see
`LICENSE-NARFIMA`.
