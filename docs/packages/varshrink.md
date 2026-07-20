# VARshrink

[Back to the implemented-package index](../../README.md#implemented-packages).

`varshrink.f90` translates the
[VARshrink](https://cran.r-project.org/web/packages/VARshrink/index.html) R
package. It provides
multivariate ridge regression over supplied shrinkage paths, effective degrees
of freedom, generalized cross-validation selection, optional predictor
standardization, and direct VAR-form ridge fitting. The semiparametric Bayesian
regression core implements fixed-lambda conjugate and non-conjugate covariance
priors, Gaussian or multivariate Student-t innovations, iterative latent
weights, covariance estimation, fitted values, residuals, and convergence
diagnostics. The nonparametric estimator includes the Opgen-Rhein/Strimmer
James-Stein covariance calculation used by VARshrink, with separately supplied
or estimated correlation and variance shrinkage intensities. Shrinkage-specific
Gaussian and Student-t log likelihoods and the package's Stein-type
marginal-variance intensity are also included.
Simulation, companion stability, MA representations, impulse responses, and
serial-correlation tests reuse `mts_mod`, `bigtime_mod`, and the shared
diagnostics instead of duplicating those algorithms. The translation is
licensed under GPL-3.0-or-later; see
`LICENSE-VARSHRINK`.
