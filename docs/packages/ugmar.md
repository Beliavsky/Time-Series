# uGMAR

[Back to the implemented-package index](../../README.md#implemented-packages).

`ugmar.f90` translates the
[uGMAR](https://cran.r-project.org/web/packages/uGMAR/index.html) R package.
It specializes the shared Gaussian and Student-t mixture autoregressive engine
to a univariate interface for GMAR, StMAR, and G-StMAR models. The module
provides validated model construction, stationarity tests, endogenous mixing
weights, conditional regime and mixture moments, conditional or exact log
likelihoods, stationary moments, BFGS, genetic, multistart, and restricted
estimation, stationary random and local smart-mutation initializers, recursive
simulation, and Monte Carlo forecasts. Reversible
intercept and stationary-mean parameterizations, numerical Hessian inference,
Wald and likelihood-ratio tests, transformed-coordinate likelihood profiles,
AR characteristic-root reporting, and high-degree Student-to-Gaussian regime
conversion are included. Gaussian
quantile residuals support residual tests with optional parameter-estimation
correction. Principal interfaces are available through `arma_mod`,
`forecasting_mod`, and `diagnostics_mod`. The translation is licensed under
GPL-3.0-only; see `LICENSE-UGMAR`.
