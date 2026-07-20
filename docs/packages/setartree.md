# setartree

[Back to the implemented-package index](../../README.md#implemented-packages).

`setartree.f90` translates the
[setartree](https://cran.r-project.org/web/packages/setartree/index.html) R
package. It fits global SETAR trees to ordinary regression matrices or lag
windows pooled across multiple series. Every node searches equally spaced
thresholds over all features using binned cumulative `X'X`, `X'y`, `y'y`, and
sample-count sufficient statistics, and terminal
pooled regressions are selected by the package's nested F test, relative SSE
improvement, or both criteria with depth-dependent significance. Flat binary
node storage supports direct traversal, point prediction, leaf residual
scales, normal prediction intervals, and recursive multi-series forecasts.
SETAR forests use sampling without replacement, optional randomized stopping
hyperparameters with independent controls for significance, its per-depth
divider, and the error-improvement threshold, averaged point predictions, and
the package's pooled
within-leaf variance intervals. Optional observation weights propagate through
node regression, cumulative split statistics, effective-size stopping tests,
leaf residual scales, forest subsampling, and pooled forest variances.
Integer-coded categorical predictors use first-seen level ordering and omit the
last level as the reference. Trees and forests retain the training levels,
reproduce their indicator columns during prediction, and reject unseen levels.
Per-series mean normalization and per-window
normalization are supported by the pooled-series interfaces.
Rectangular observation-by-series matrices and reusable ragged
`real_vector_t` collections support equal- or unequal-length input histories.
Principal
interfaces are available through `arma_mod`, `forecasting_mod`, and
`regression_time_series_mod`. The translation is licensed under MIT; see
`LICENSE-SETARTREE`.
