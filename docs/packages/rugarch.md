# rugarch

[Back to the implemented-package index](../../README.md#implemented-packages).

`rugarch.f90` translates the numerical core of the
[rugarch](https://cran.r-project.org/web/packages/rugarch/index.html) R package.
The initial layer provides `ugarchspec`, filtering, likelihood estimation,
Hessian inference, analytic forecasting, and simulation for sGARCH, integrated
GARCH, exponential GARCH, GJR-GARCH, asymmetric power ARCH, truncated
fiGARCH(1,d,1), component GARCH, log-linear realized GARCH, and the Hentschel
fGARCH omnibus family. The
nonfractional variance families support general ARCH and GARCH orders where the
model permits them. Realized GARCH jointly evaluates returns and a positive
realized-variance series through its leverage-aware measurement equation.
Conditional means support ARMA and fractional ARFIMA dynamics, optional means,
ARCH-in-mean effects, and external regressors; external variance regressors are
also supported. Gaussian, variance-standardized
Student-t, generalized-error, Fernandez-Steel skew-normal, skew-Student,
skew-GED, Johnson SU, standardized normal-inverse-Gaussian, generalized
hyperbolic, and generalized-hyperbolic skew-Student innovations are available,
together with persistence, unconditional-variance, half-life, and
information-criterion calculations. Standalone diagnostics translate the
Berkowitz density-calibration test, Pesaran-Timmermann and Anatolyev-Gerko
directional-accuracy tests, Kupiec and Christoffersen VaR tests, and the
conditional expected-shortfall test. Core variance formulas were adapted from
the MIT-licensed GARCH-BFGS project and integrated with this library's shared
optimizer, random-number generator, statistics, and linear algebra. Modeling
interfaces are available through `volatility_mod`, and forecasts through
`forecasting_mod`; calibration and risk tests are available through
`diagnostics_mod`. The translation is licensed under GPL-3.0-only; see
`LICENSE-RUGARCH`.

`rugarch_extensions.f90` adds parametric bootstrap intervals, expanding and
moving-window rolling refits, multi-specification fitting, asymptotic parameter
simulation, multi-model filtering and forecasting, loss-based model confidence
sets, news-impact curves, and Engle-Ng sign and size-bias diagnostics. Bootstrap
paths recursively update the fitted conditional mean and variance after each
simulated innovation. These graphics-independent workflows are also
exported through `volatility_mod`. The translation is licensed under
GPL-3.0-only; see `LICENSE-RUGARCH`.

`rugarch_diagnostics.f90` translates the Hansen-Nyblom parameter-stability,
grouped probability-transform goodness-of-fit, censored-Weibull VaR-duration,
GMM orthogonality, and quartic-kernel Hong-Li specification tests. The tests are
exported through `diagnostics_mod`; fitted rugarch objects retain numerical
per-observation likelihood scores for direct Nyblom testing. The translation is licensed under
GPL-3.0-only; see `LICENSE-RUGARCH`.

`distribution.f90` centralizes the standardized innovation densities,
parameter transformations, and random generation used by `rugarch.f90`.
Modified-Bessel approximations and the generalized inverse-Gaussian
ratio-of-uniforms sampler reuse MIT-licensed numerical work from GARCH-BFGS;
the combined distribution translation is licensed under GPL-3.0-only.

`distribution_fit.f90` provides constrained maximum-likelihood location-scale
fitting and covariance inference for the standardized innovation families in
`distribution.f90`. The translation is licensed under GPL-3.0-only; see
`LICENSE-RUGARCH`.
