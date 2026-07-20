# Fortran Time-Series Library

Fortran 2008 translations of numerical algorithms from R time-series packages.
The implemented-package index links each translation to its detailed coverage,
provenance, and licensing notes.

## Shared Infrastructure

- `kind.f90`: shared `kind_mod` module defining the double-precision kind `dp`.
- `utils.f90`: IEEE NaN and option-string helpers.
- `calendar.f90`: Gregorian `date_t`, parsing, arithmetic,
  comparisons, ISO weekdays, ordinal days, leap years, and Easter dates.
- `linalg.f90`: inversion, log determinants, Cholesky and symmetric
  eigen decompositions, Householder-reduced general eigenvalues, rank
  estimation, and matrix helpers.
- `stats.f90`: sorting, quantiles, descriptive statistics, ordinary least
  squares, and regression RSS.
- `spline.f90`: quantile-knot B-spline bases, finite-difference penalties,
  penalized least squares, covariance estimates, effective degrees of freedom,
  GCV smoothing selection, and prediction.
- `random.f90`: seeding, uniform and standard-normal generation,
  normal matrices, and multivariate-normal transformations.
- `fourier.f90`, `polynomial.f90`, and `special_functions.f90`: reusable
  transforms, polynomial operations, and distribution-related functions.
- `optimization.f90`: finite-difference BFGS with Armijo search.
- `time_series_stats.f90`: ACF, PACF, CCF, AR estimation, and harmonic
  regression.
- `time_series_diagnostics.f90`: weighted univariate portmanteau tests and
  FCVAR-style univariate and multivariate white-noise Q and robust LM tests.

Package implementations reuse this layer and earlier package modules. In
particular, the astsa state-space adapter is implemented over `kfas_mod`.

## Topic Facades

Package modules preserve translation provenance and expose complete
package-specific APIs. The following implementation-free facade modules provide
curated, topic-based entry points for applications and examples:

- `arma_mod`: AR, threshold AR, MA, ARMA, ARIMA, and SARIMA models.
- `long_memory_mod`: fractional differencing and long-memory models.
- `state_space_mod`: filtering, smoothing, structural models, and particle methods.
- `multivariate_mod`: VAR, VARMA, VECM, cointegration, and factor models.
- `volatility_mod`: ARCH, multivariate GARCH, and stochastic volatility.
- `spectral_mod`: Fourier transforms, periodograms, and spectral estimation.
- `diagnostics_mod`: residual tests, portmanteau tests, and forecast diagnostics.
- `forecasting_mod`: univariate, multivariate, and structural forecasting.
- `regression_time_series_mod`: dynamic regression, transfer functions,
  interventions, and calendar regressors.
- `bayesian_time_series_mod`: Bayesian estimation, MCMC, and predictive methods.
- `markov_switching_mod`: Markov-switching filtering, smoothing, and estimation.
- `functional_time_series_mod`: functional autoregression, estimation, tests,
  and forecasting.

Facade modules contain no numerical implementations. New package translations
should remain in package modules and add their principal user-facing procedures
to the relevant facades. Because some algorithms belong to multiple topics,
applications importing several facades should normally use explicit `only:`
lists.

Run `python tools/fortran_style_audit.py` to audit all Fortran sources for the
project's mechanical style rules and CMake coverage. Add `--suggest-purity` for
non-failing, heuristic `pure` and `elemental` candidates, or `--json` for a
machine-readable report. The command exits with status 1 when definite style
errors are found.

## Source Licensing

This is a multi-license source tree. Each Fortran source starts with an
`SPDX-License-Identifier` naming its applicable license and an
`SPDX-FileComment` recording its origin. Package translations retain the
upstream package license. Original shared infrastructure is MIT licensed, and
tests use the license of the translated module they exercise.

New Fortran sources must put these two comments on lines 1 and 2. The style
audit reports `F011` or `F012` when either comment is missing. The
Package-specific notices and license texts are retained under `licenses/`.

## Implemented Packages

[`tsdyn.f90`](docs/packages/tsdyn.md) translates the
[tsDyn](https://cran.r-project.org/package=tsDyn) R package
(*Nonlinear Time Series Models with Regime Switching*).

[`setartree.f90`](docs/packages/setartree.md) translates the
[setartree](https://cran.r-project.org/package=setartree) R package
(*SETAR-Tree - A Novel and Accurate Tree Algorithm for Global Time Series Forecasting*).

[`tseriestarma.f90`](docs/packages/tseriestarma.md) translates the
[tseriesTARMA](https://cran.r-project.org/package=tseriesTARMA) R package
(*Analysis of Nonlinear Time Series Through Threshold Autoregressive Moving Average Models (TARMA) Models*).

[`gmdh.f90`](docs/packages/gmdh.md) translates the
[GMDH](https://cran.r-project.org/package=GMDH) R package
(*Short Term Forecasting via GMDH-Type Neural Network Algorithms*).

[`nnfor.f90`](docs/packages/nnfor.md) translates the
[nnfor](https://cran.r-project.org/package=nnfor) R package
(*Time Series Forecasting with Neural Networks*).

[`narfima.f90`](docs/packages/narfima.md) translates the
[narfima](https://cran.r-project.org/package=narfima) R package
(*Neural AutoRegressive Fractionally Integrated Moving Average Model*).

[`nlints.f90`](docs/packages/nlints.md) translates the
[NlinTS](https://cran.r-project.org/package=NlinTS) R package
(*Models for Non Linear Causality Detection in Time Series*).

[`tslstmplus.f90`](docs/packages/tslstmplus.md) translates the
[TSLSTMplus](https://cran.r-project.org/package=TSLSTMplus) R package
(*Long-Short Term Memory for Time-Series Forecasting, Enhanced*).

[`tsann.f90`](docs/packages/tsann.md) translates the
[TSANN](https://cran.r-project.org/package=TSANN) R package
(*Time Series Artificial Neural Network*).

[`echos.f90`](docs/packages/echos.md) translates the
[echos](https://cran.r-project.org/package=echos) R package
(*Echo State Networks for Time Series Modeling and Forecasting*).

[`starvars.f90`](docs/packages/starvars.md) translates the
[starvars](https://cran.r-project.org/package=starvars) R package
(*Vector Logistic Smooth Transition Models Estimation and Prediction*).

[`rugarch.f90`](docs/packages/rugarch.md) translates the
[rugarch](https://cran.r-project.org/package=rugarch) R package
(*Univariate GARCH Models*).

[`expar.f90`](docs/packages/expar.md) translates the
[EXPAR](https://cran.r-project.org/package=EXPAR) R package
(*Fitting of Exponential Autoregressive (EXPAR) Model*).

[`exparma.f90`](docs/packages/exparma.md) translates the
[EXPARMA](https://cran.r-project.org/package=EXPARMA) R package
(*Fitting of Exponential Autoregressive Moving Average (EXPARMA) Model*).

[`bentcablear.f90`](docs/packages/bentcablear.md) translates the
[bentcableAR](https://cran.r-project.org/package=bentcableAR) R package
(*Bent-Cable Regression for Independent Data or Autoregressive Time Series*).

[`baystar.f90`](docs/packages/baystar.md) translates the
[BAYSTAR](https://cran.r-project.org/package=BAYSTAR) R package
(*On Bayesian Analysis of Threshold Autoregressive Models*).

[`mixar.f90`](docs/packages/mixar.md) translates the
[mixAR](https://cran.r-project.org/package=mixAR) R package
(*Mixture Autoregressive Models*).

[`mswm.f90`](docs/packages/mswm.md) translates the
[MSwM](https://cran.r-project.org/package=MSwM) R package
(*Fitting Markov Switching Models*).

[`nts.f90`](docs/packages/nts.md) translates the
[NTS](https://cran.r-project.org/package=NTS) R package
(*Nonlinear Time Series Analysis*).

[`bvar.f90`](docs/packages/bvar.md) translates the
[BVAR](https://cran.r-project.org/package=BVAR) R package
(*Hierarchical Bayesian Vector Autoregression*).

[`bvartools.f90`](docs/packages/bvartools.md) translates the
[bvartools](https://cran.r-project.org/package=bvartools) R package
(*Bayesian Inference of Vector Autoregressive and Error Correction Models*).

[`var_etp.f90`](docs/packages/var_etp.md) translates the
[VAR.etp](https://cran.r-project.org/package=VAR.etp) R package
(*VAR Modelling: Estimation, Testing, and Prediction*).

[`bvars.f90`](docs/packages/bvars.md) translates the
[bvars](https://cran.r-project.org/package=bvars) R package
(*Bayesian Forecasting with Large Vector Autoregressions*).

[`bvarsv.f90`](docs/packages/bvarsv.md) translates the
[bvarsv](https://cran.r-project.org/package=bvarsv) R package
(*Bayesian Analysis of a Vector Autoregressive Model with Stochastic Volatility and Time-Varying Parameters*).

[`gmvarkit.f90`](docs/packages/gmvarkit.md) translates the
[gmvarkit](https://cran.r-project.org/package=gmvarkit) R package
(*Estimate Gaussian and Student's t Mixture Vector Autoregressive Models*).

[`ugmar.f90`](docs/packages/ugmar.md) translates the
[uGMAR](https://cran.r-project.org/package=uGMAR) R package
(*Estimate Univariate Gaussian and Student's t Mixture Autoregressive Models*).

[`vars.f90`](docs/packages/vars.md) translates the
[vars](https://cran.r-project.org/package=vars) R package
(*VAR Modelling*).

[`varshrink.f90`](docs/packages/varshrink.md) translates the
[VARshrink](https://cran.r-project.org/package=VARshrink) R package
(*Shrinkage Estimation Methods for Vector Autoregressive Models*).

[`bigvar.f90`](docs/packages/bigvar.md) translates the
[BigVAR](https://cran.r-project.org/package=BigVAR) R package
(*Dimension Reduction Methods for Multivariate Time Series*).

[`bigtime.f90`](docs/packages/bigtime.md) translates the
[bigtime](https://cran.r-project.org/package=bigtime) R package
(*Sparse Estimation of Large Time Series Models*).

[`fcvar.f90`](docs/packages/fcvar.md) translates the
[FCVAR](https://cran.r-project.org/package=FCVAR) R package
(*Estimation and Inference for the Fractionally Cointegrated VAR*).

[`mar.f90`](docs/packages/mar.md) translates the
[mAr](https://cran.r-project.org/package=mAr) R package
(*Multivariate AutoRegressive Analysis*).

[`forecast.f90`](docs/packages/forecast.md) translates the
[forecast](https://cran.r-project.org/package=forecast) R package
(*Forecasting Functions for Time Series and Linear Models*).

[`kfas.f90`](docs/packages/kfas.md) translates the
[KFAS](https://cran.r-project.org/package=KFAS) R package
(*Kalman Filter and Smoother for Exponential Family State Space Models*).

[`urca.f90`](docs/packages/urca.md) translates the
[urca](https://cran.r-project.org/package=urca) R package
(*Unit Root and Cointegration Tests for Time Series Data*).

[`astsa.f90`](docs/packages/astsa.md) translates the
[astsa](https://cran.r-project.org/package=astsa) R package
(*Applied Statistical Time Series Analysis*).

[`itsmr.f90`](docs/packages/itsmr.md) translates the
[itsmr](https://cran.r-project.org/package=itsmr) R package
(*Time Series Analysis Using the Innovations Algorithm*).

[`arima2.f90`](docs/packages/arima2.md) translates the
[arima2](https://cran.r-project.org/package=arima2) R package
(*Likelihood Based Inference for ARIMA Modeling*).

[`fracdiff.f90`](docs/packages/fracdiff.md) translates the
[fracdiff](https://cran.r-project.org/package=fracdiff) R package
(*Fractionally Differenced ARIMA aka ARFIMA(P,d,q) Models*).

[`fracdist.f90`](docs/packages/fracdist.md) translates the
[fracdist](https://cran.r-project.org/package=fracdist) R package
(*Numerical CDFs for Fractional Unit Root and Cointegration Tests*).

[`nsarfima.f90`](docs/packages/nsarfima.md) translates the
[nsarfima](https://cran.r-project.org/package=nsarfima) R package
(*Methods for Fitting and Simulating Non-Stationary ARFIMA Models*).

[`garma.f90`](docs/packages/garma.md) translates the
[garma](https://cran.r-project.org/package=garma) R package
(*Fitting and Forecasting Gegenbauer ARMA Time Series Models*).

[`esemifar.f90`](docs/packages/esemifar.md) translates the
[esemifar](https://cran.r-project.org/package=esemifar) R package
(*Smoothing Long-Memory Time Series*).

[`tfarima.f90`](docs/packages/tfarima.md) translates the
[tfarima](https://cran.r-project.org/package=tfarima) R package
(*Transfer Function and ARIMA Models*).

[`robustarima.f90`](docs/packages/robustarima.md) translates the
[robustarima](https://cran.r-project.org/package=robustarima) R package
(*Robust ARIMA Modeling*).

[`bsts.f90`](docs/packages/bsts.md) translates the
[bsts](https://cran.r-project.org/package=bsts) R package
(*Bayesian Structural Time Series*).

[`bssm.f90`](docs/packages/bssm.md) translates the
[bssm](https://cran.r-project.org/package=bssm) R package
(*Bayesian Inference of Non-Linear and Non-Gaussian State Space Models*).

[`mts.f90`](docs/packages/mts.md) translates the
[MTS](https://cran.r-project.org/package=MTS) R package
(*All-Purpose Toolkit for Analyzing Multivariate Time Series (MTS) and Estimating Multivariate Volatility Models*).

## Build

```sh
cmake -S . -B build -G Ninja -DCMAKE_Fortran_COMPILER=gfortran
cmake --build build
ctest --test-dir build --output-on-failure
```

Run the VAR(2) simulation and fitting example with:

```sh
./build/example_var2_fit
```

The translated source packages use GPL licenses. Redistribution of derived code
must preserve the applicable license terms.

The shared date implementation is adapted from the MIT-licensed DataFrame
project. Its notice is preserved in `licenses/LICENSE-DATE`.
