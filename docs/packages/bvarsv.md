# bvarsv

[Back to the implemented-package index](../../README.md#implemented-packages).

`bvarsv.f90` translates the
[bvarsv](https://cran.r-project.org/web/packages/bvarsv/index.html) R package.
It provides the Primiceri time-varying structural VAR
algorithms, including the multivariate
Carter-Kohn random-walk state sampler, packed contemporaneous-matrix and
structural-covariance reconstruction, and the seven-component KSC update for
correlated random-walk log-volatility states. Equation-block Carter-Kohn updates
draw the packed time-varying contemporaneous coefficients. Reusable
inverse-Wishart updates draw the coefficient and volatility innovation
covariances and the equation-block contemporaneous innovation covariances. The
package-specific equation-stacked VAR design and time-varying coefficient-state
update are also provided. The training-sample OLS prior reproduces the package's
GLS coefficient moments, recursive covariance decomposition, initial log
variances, and inverse-Wishart Monte Carlo covariance for contemporaneous states.
`bvarsv_gibbs` assembles these kernels into the complete Primiceri sampler with
configurable burn-in, thinning, and innovation-prior scales, retaining all state
paths and the Q, S, and W covariance draws.
`bvarsv_forecast` produces recursive draw-wise predictive means, covariance
matrices, and realizations with either random-walk parameter drift or terminal
parameters held fixed.
`bvarsv_irf` computes draw-wise responses at any retained state time with the
package's identity, Cholesky, or Primiceri average-volatility impact scenarios.
`bvarsv_predictive_density` evaluates the draw-averaged Gaussian predictive
PDF or CDF at user-supplied points.
`bvarsv_simulate_var1` generates the package's TVP-VAR(1) process with
random-walk coefficients, contemporaneous relations, and log variances.
`bvarsv_predictive_draws` and `bvarsv_parameter_draws` extract forecast draws
and selected time-varying parameter paths without the R package's plotting
layer. The translation is licensed under GPL-2.0-or-later.
