# BAYSTAR

[Back to the implemented-package index](../../README.md#implemented-packages).

`baystar.f90` translates the
[BAYSTAR](https://cran.r-project.org/web/packages/BAYSTAR/index.html) R
package. It provides sparse-lag two-regime TAR likelihoods, conjugate Gaussian
coefficient and inverse-gamma variance posteriors, decreasing-prior delay-lag
probabilities, bounded random-walk threshold Metropolis updates, shared-RNG
Gibbs sampling, supplied-innovation and shared-RNG simulation, posterior
summaries, regime means, residual reconstruction, threshold acceptance rates,
modal delay selection, and DIC. Internal or externally supplied threshold
variables are supported. Principal interfaces are available through both
`arma_mod` and `bayesian_time_series_mod`. The translation is licensed under
GPL-2.0-or-later; see `LICENSE-BAYSTAR`.
