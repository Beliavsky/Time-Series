! SPDX-License-Identifier: Artistic-2.0
! SPDX-FileComment: Algorithms translated from the R MTS package.
! Numerical translations from the Artistic-2.0 MTS package.
module mts_mod
   use kind_mod, only: dp
   use linalg_mod, only: invert_matrix, inverse_logdet, cholesky_lower, &
      symmetric_eigen, identity_matrix
   use random_mod, only: random_multivariate_normal
   use stats_mod, only: normal_quantile, sort
   use stats_mod, only: sample_covariance_shared => covariance
   use stats_mod, only: data_correlation_matrix => correlation_matrix
   use optimization_mod, only: optimization_result_t, bfgs_minimize_fd, &
      finite_difference_hessian
   use urca_mod, only: johansen_result_t, johansen_test
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   private

   type, public :: mts_var_fit_t
      ! Fitted VAR coefficients, inference, residuals, covariance, and criteria.
      real(dp), allocatable :: coefficients(:, :), standard_errors(:, :)
      real(dp), allocatable :: ar(:, :, :), intercept(:), residuals(:, :), sigma(:, :)
      integer, allocatable :: lags(:)
      real(dp) :: aic = 0.0_dp
      real(dp) :: bic = 0.0_dp
      real(dp) :: hq = 0.0_dp
      integer :: max_lag = 0
      integer :: info = 0
      logical :: includes_mean = .true.
   end type

   type, public :: mts_var_forecast_t
      ! VAR point forecasts and innovation forecast-error uncertainty.
      real(dp), allocatable :: mean(:, :), standard_error(:, :), covariance(:, :, :)
      integer :: info = 0
   end type

   type, public :: mts_var_order_t
      ! MTS VAR order criteria and sequential likelihood-ratio statistics.
      real(dp), allocatable :: aic(:), bic(:), hq(:), statistic(:), p_value(:)
      integer :: aic_order = 0
      integer :: bic_order = 0
      integer :: hq_order = 0
      integer :: info = 0
   end type

   type, public :: mts_var_irf_t
      ! Raw, orthogonalized, generalized, and cumulative VAR impulse responses.
      real(dp), allocatable :: psi(:, :, :), orthogonal(:, :, :), generalized(:, :, :)
      real(dp), allocatable :: cumulative_orthogonal(:, :, :), cumulative_generalized(:, :, :)
      real(dp), allocatable :: shock_factor(:, :)
      integer :: info = 0
   end type

   type, public :: mts_fevd_t
      ! Forecast-error variance contributions and total variances by horizon.
      real(dp), allocatable :: contribution(:, :, :), variance(:, :), standard_error(:, :)
      logical :: generalized = .false.
      integer :: info = 0
   end type

   type, public :: mts_varx_fit_t
      ! Fitted VARX coefficients, residual inference, and information criteria.
      real(dp), allocatable :: ar(:, :, :), exogenous(:, :, :), intercept(:)
      real(dp), allocatable :: coefficients(:, :), standard_errors(:, :)
      real(dp), allocatable :: residuals(:, :), sigma(:, :)
      integer :: ar_order = 0
      integer :: exogenous_order = 0
      integer :: info = 0
      logical :: includes_mean = .true.
      real(dp) :: aic = 0.0_dp
      real(dp) :: bic = 0.0_dp
      real(dp) :: hq = 0.0_dp
   end type

   type, public :: mts_varx_irf_t
      ! Endogenous and exogenous dynamic responses and cumulative responses.
      real(dp), allocatable :: endogenous(:, :, :), exogenous(:, :, :)
      real(dp), allocatable :: cumulative_endogenous(:, :, :), cumulative_exogenous(:, :, :)
      integer :: info = 0
   end type

   type, public :: mts_varx_order_t
      ! VARX information-criterion grids and selected lag pairs.
      real(dp), allocatable :: aic(:, :), bic(:, :), hq(:, :)
      integer :: aic_order(2) = 0
      integer :: bic_order(2) = 0
      integer :: hq_order(2) = 0
      integer :: info = 0
   end type

   type, public :: mts_varma_model_t
      ! VARMA coefficient matrices under MTS's subtractive MA convention.
      real(dp), allocatable :: ar(:, :, :), ma(:, :, :), intercept(:), sigma(:, :)
      integer :: info = 0
   end type

   type, public :: mts_varma_covariance_t
      ! Theoretical VARMA autocovariance and cross-correlation matrices.
      real(dp), allocatable :: covariance(:, :, :), correlation(:, :, :)
      integer :: truncation = 0
      integer :: info = 0
   end type

   type, public :: mts_varma_simulation_t
      ! Simulated VARMA observations, retained innovations, and status.
      real(dp), allocatable :: series(:, :), innovations(:, :)
      integer :: info = 0
   end type

   type, public :: mts_vma_fit_t
      ! Conditional Gaussian VMA estimates, inference, residuals, and diagnostics.
      type(mts_varma_model_t) :: model
      real(dp), allocatable :: coefficients(:), standard_errors(:), covariance(:, :)
      real(dp), allocatable :: residuals(:, :)
      logical, allocatable :: estimated(:)
      integer, allocatable :: lags(:)
      real(dp) :: log_likelihood = 0.0_dp
      real(dp) :: aic = 0.0_dp
      real(dp) :: bic = 0.0_dp
      real(dp) :: invertibility_radius = 0.0_dp
      integer :: iterations = 0
      integer :: info = 0
      logical :: converged = .false.
   end type

   type, public :: mts_varma_fit_t
      ! Conditional Gaussian VARMA estimates, inference, residuals, and diagnostics.
      type(mts_varma_model_t) :: model
      real(dp), allocatable :: coefficients(:), standard_errors(:), covariance(:, :)
      real(dp), allocatable :: residuals(:, :)
      logical, allocatable :: estimated(:)
      integer, allocatable :: ar_lags(:), ma_lags(:)
      real(dp) :: log_likelihood = 0.0_dp
      real(dp) :: aic = 0.0_dp
      real(dp) :: bic = 0.0_dp
      real(dp) :: stationarity_radius = 0.0_dp
      real(dp) :: invertibility_radius = 0.0_dp
      integer :: iterations = 0
      integer :: info = 0
      logical :: converged = .false.
   end type

   type, public :: mts_svarma_model_t
      ! Multiplicative regular and seasonal VARMA components and expanded model.
      real(dp), allocatable :: regular_ar(:, :, :), seasonal_ar(:, :, :)
      real(dp), allocatable :: regular_ma(:, :, :), seasonal_ma(:, :, :)
      real(dp), allocatable :: intercept(:), sigma(:, :)
      type(mts_varma_model_t) :: expanded
      integer :: period = 0
      logical :: switched = .false.
      integer :: info = 0
   end type

   type, public :: mts_svarma_fit_t
      ! Conditional Gaussian seasonal VARMA estimates and diagnostics.
      type(mts_svarma_model_t) :: model
      real(dp), allocatable :: coefficients(:), standard_errors(:), covariance(:, :)
      real(dp), allocatable :: residuals(:, :)
      logical, allocatable :: estimated(:)
      real(dp) :: log_likelihood = 0.0_dp
      real(dp) :: aic = 0.0_dp
      real(dp) :: bic = 0.0_dp
      real(dp) :: stationarity_radius = 0.0_dp
      real(dp) :: invertibility_radius = 0.0_dp
      integer :: iterations = 0
      integer :: info = 0
      logical :: converged = .false.
   end type

   type, public :: mts_varma_refinement_t
      ! Refined VARMA fit and backward-elimination history.
      type(mts_varma_fit_t) :: fit
      integer, allocatable :: active_count(:), removed_index(:)
      integer :: steps = 0
      integer :: info = 0
   end type

   type, public :: mts_svarma_refinement_t
      ! Refined seasonal VARMA fit and backward-elimination history.
      type(mts_svarma_fit_t) :: fit
      integer, allocatable :: active_count(:), removed_index(:)
      integer :: steps = 0
      integer :: info = 0
   end type

   type, public :: mts_vecm_fit_t
      ! Fitted VECM coefficients, inference, residuals, and level-VAR form.
      real(dp), allocatable :: cointegration(:, :), loading(:, :), gamma(:, :, :)
      real(dp), allocatable :: intercept(:), coefficients(:, :), standard_errors(:, :)
      real(dp), allocatable :: residuals(:, :), sigma(:, :)
      logical, allocatable :: estimated(:, :)
      type(mts_var_fit_t) :: level_var
      real(dp) :: aic = 0.0_dp
      real(dp) :: bic = 0.0_dp
      integer :: rank = 0
      integer :: level_order = 0
      integer :: info = 0
      logical :: includes_constant = .false.
   end type

   type, public :: mts_vecm_forecast_t
      ! VECM level and differenced forecasts with level uncertainty.
      type(mts_var_forecast_t) :: level
      real(dp), allocatable :: difference(:, :)
      integer :: info = 0
   end type

   type, public :: mts_factor_model_t
      ! Principal-component factor model and factor-count criteria.
      real(dp), allocatable :: mean(:), scale(:), loadings(:, :), scores(:, :)
      real(dp), allocatable :: common(:, :), residuals(:, :), eigenvalues(:)
      real(dp), allocatable :: explained(:), ic1(:), ic2(:), ic3(:)
      integer :: factors = 0
      integer :: ic1_factors = 0
      integer :: ic2_factors = 0
      integer :: ic3_factors = 0
      integer :: info = 0
      logical :: standardized = .true.
   end type

   type, public :: mts_factor_forecast_t
      ! Factor-VAR forecasts reconstructed in the original variable scale.
      type(mts_var_fit_t) :: factor_var
      type(mts_var_forecast_t) :: factor_forecast
      real(dp), allocatable :: mean(:, :), standard_error(:, :), covariance(:, :, :)
      integer :: info = 0
   end type

   type, public :: mts_constrained_factor_t
      ! Weighted least-squares factor model under caller-supplied constraints.
      type(mts_factor_model_t) :: factor_model
      real(dp), allocatable :: constraint(:, :), omega(:, :), psi(:, :)
      real(dp), allocatable :: constrained_eigenvalues(:)
      real(dp) :: explained = 0.0_dp
      integer :: info = 0
   end type

   type, public :: mts_bvar_prior_t
      ! Matrix-normal coefficient and inverse-Wishart covariance prior.
      real(dp), allocatable :: mean(:, :), precision(:, :), scale(:, :)
      real(dp) :: degrees_of_freedom = 0.0_dp
      integer :: info = 0
   end type

   type, public :: mts_bvar_fit_t
      ! Bayesian VAR posterior moments, residuals, and shared VAR representation.
      type(mts_var_fit_t) :: model
      type(mts_bvar_prior_t) :: prior
      real(dp), allocatable :: coefficient_mean(:, :), coefficient_covariance(:, :)
      real(dp), allocatable :: standard_errors(:, :), posterior_precision(:, :)
      real(dp), allocatable :: residuals(:, :), sigma(:, :)
      real(dp) :: posterior_degrees_of_freedom = 0.0_dp
      integer :: info = 0
   end type

   type, public :: mts_common_volatility_t
      ! Common-volatility eigensystem, transformations, and ARCH diagnostics.
      real(dp), allocatable :: residuals(:, :), whitened(:, :), aggregate(:, :)
      real(dp), allocatable :: eigenvalues(:), proportions(:), eigenvectors(:, :)
      real(dp), allocatable :: directions(:, :), components(:, :)
      real(dp), allocatable :: arch_statistic(:, :), arch_p_value(:, :)
      integer, allocatable :: arch_lags(:)
      integer :: prewhiten_order = 0
      integer :: max_lag = 0
      integer :: info = 0
      logical :: standardized = .false.
   end type

   type, public :: mts_mch_diagnostic_t
      ! Diagnostics for multivariate conditional covariance model residuals.
      real(dp), allocatable :: standardized_residuals(:, :), radial_residual(:)
      real(dp) :: radial_q = 0.0_dp
      real(dp) :: rank_q = 0.0_dp
      real(dp) :: multivariate_q = 0.0_dp
      real(dp) :: robust_q = 0.0_dp
      real(dp) :: p_value(4) = 1.0_dp
      integer :: degrees_of_freedom(4) = 0
      integer :: max_lag = 0
      integer :: robust_observations = 0
      integer :: info = 0
   end type

   type, public :: mts_bekk_fit_t
      ! Gaussian BEKK(1,1) estimates, inference, and conditional covariance path.
      real(dp), allocatable :: mean(:), constant(:, :), arch(:, :), garch(:, :)
      real(dp), allocatable :: coefficients(:), standard_errors(:), parameter_covariance(:, :)
      real(dp), allocatable :: residuals(:, :), standardized_residuals(:, :)
      real(dp), allocatable :: covariance(:, :, :)
      logical, allocatable :: estimated(:)
      real(dp) :: log_likelihood = 0.0_dp
      real(dp) :: aic = 0.0_dp
      real(dp) :: bic = 0.0_dp
      real(dp) :: persistence = 0.0_dp
      integer :: iterations = 0
      integer :: info = 0
      logical :: converged = .false.
      logical :: includes_mean = .true.
   end type

   type, public :: mts_dcc_fit_t
      ! Gaussian DCC(1,1) correlation estimates and covariance paths.
      real(dp), allocatable :: standardized_residuals(:, :), residuals(:, :)
      real(dp), allocatable :: marginal_variance(:, :), unconditional(:, :)
      real(dp), allocatable :: q(:, :, :), correlation(:, :, :), covariance(:, :, :)
      real(dp), allocatable :: parameter_covariance(:, :), standard_errors(:)
      real(dp) :: arch = 0.0_dp
      real(dp) :: garch = 0.0_dp
      real(dp) :: log_likelihood = 0.0_dp
      real(dp) :: aic = 0.0_dp
      real(dp) :: bic = 0.0_dp
      integer :: iterations = 0
      integer :: info = 0
      logical :: converged = .false.
   end type

   type, public :: mts_adcc_fit_t
      ! Gaussian or Student-t asymmetric DCC estimates and covariance paths.
      real(dp), allocatable :: standardized_residuals(:, :), residuals(:, :)
      real(dp), allocatable :: marginal_variance(:, :), unconditional(:, :), negative_unconditional(:, :)
      real(dp), allocatable :: q(:, :, :), correlation(:, :, :), covariance(:, :, :)
      real(dp), allocatable :: parameter_covariance(:, :), standard_errors(:)
      real(dp) :: arch = 0.0_dp
      real(dp) :: garch = 0.0_dp
      real(dp) :: asymmetry = 0.0_dp
      real(dp) :: degrees_of_freedom = 0.0_dp
      real(dp) :: log_likelihood = 0.0_dp
      real(dp) :: aic = 0.0_dp
      real(dp) :: bic = 0.0_dp
      integer :: iterations = 0
      integer :: info = 0
      logical :: converged = .false.
      logical :: student_t = .false.
   end type

   type, public :: mts_tse_tsui_fit_t
      ! Gaussian or Student-t Tse-Tsui rolling-correlation DCC estimates.
      real(dp), allocatable :: standardized_residuals(:, :), residuals(:, :)
      real(dp), allocatable :: marginal_variance(:, :), unconditional(:, :)
      real(dp), allocatable :: local_correlation(:, :, :), correlation(:, :, :), covariance(:, :, :)
      real(dp), allocatable :: parameter_covariance(:, :), standard_errors(:)
      real(dp) :: previous_weight = 0.0_dp
      real(dp) :: rolling_weight = 0.0_dp
      real(dp) :: degrees_of_freedom = 0.0_dp
      real(dp) :: log_likelihood = 0.0_dp
      real(dp) :: aic = 0.0_dp
      real(dp) :: bic = 0.0_dp
      integer :: window = 0
      integer :: iterations = 0
      integer :: info = 0
      logical :: converged = .false.
      logical :: student_t = .false.
   end type

   type, public :: mts_ewma_fit_t
      ! Multivariate EWMA covariance path and optional decay estimate.
      real(dp), allocatable :: mean(:), residuals(:, :), covariance(:, :, :)
      real(dp) :: decay = 0.96_dp
      real(dp) :: standard_error = 0.0_dp
      real(dp) :: log_likelihood = 0.0_dp
      integer :: iterations = 0
      integer :: info = 0
      logical :: estimated = .false.
      logical :: converged = .true.
   end type

   type, public :: mts_mchol_fit_t
      ! Moving-Cholesky coefficients, orthogonal GARCH components, and covariances.
      real(dp), allocatable :: mean(:), prewhitened(:, :), smoothed_coefficients(:, :)
      real(dp), allocatable :: orthogonal_residuals(:, :), component_variance(:, :)
      real(dp), allocatable :: garch_parameters(:, :), covariance(:, :, :)
      integer :: window = 0
      integer :: var_order = 0
      real(dp) :: decay = 0.96_dp
      integer :: info = 0
   end type

   type, public :: mts_sccor_t
      ! Unconstrained and group-constrained sample correlation matrices.
      real(dp), allocatable :: unconstrained(:, :), constrained(:, :)
      integer :: start = 0
      integer :: end = 0
      integer :: span = 0
      integer :: info = 0
   end type

   type, public :: mts_arch_test_t
      ! McLeod-Li and rank-based univariate ARCH tests.
      real(dp) :: statistic(2) = 0.0_dp
      real(dp) :: p_value(2) = 1.0_dp
      integer :: lag = 0
      integer :: info = 0
   end type

   type, public :: mts_march_test_t
      ! Radial, rank, multivariate, and robust multivariate ARCH tests.
      real(dp) :: statistic(4) = 0.0_dp
      real(dp) :: p_value(4) = 1.0_dp
      integer :: degrees_of_freedom(4) = 0
      integer :: lag = 0
      integer :: robust_observations = 0
      integer :: info = 0
   end type

   type, public :: mts_copula_fit_t
      ! Dynamic grouped Student-t copula angles and correlation paths.
      real(dp), allocatable :: baseline_angles(:), local_angles(:, :), angles(:, :)
      real(dp), allocatable :: correlation(:, :, :), coefficients(:), standard_errors(:)
      real(dp), allocatable :: parameter_covariance(:, :)
      integer, allocatable :: groups(:)
      real(dp) :: degrees_of_freedom = 0.0_dp
      real(dp) :: previous_weight = 0.0_dp
      real(dp) :: local_weight = 0.0_dp
      real(dp) :: log_likelihood = 0.0_dp
      integer :: window = 0
      integer :: iterations = 0
      integer :: info = 0
      logical :: converged = .false.
      logical :: estimates_baseline = .true.
   end type

   type, public :: mts_missing_result_t
      ! Conditional GLS estimates and uncertainty for one missing observation.
      real(dp), allocatable :: estimate(:), covariance(:, :), completed(:, :)
      logical, allocatable :: missing(:)
      integer :: time_index = 0
      integer :: info = 0
   end type

   type, public :: mts_granger_test_t
      ! Wald test of excluded predictor blocks in designated VAR equations.
      type(mts_var_fit_t) :: unrestricted
      real(dp), allocatable :: restrictions(:), restriction_covariance(:, :)
      integer, allocatable :: targets(:), predictors(:)
      real(dp) :: statistic = 0.0_dp
      real(dp) :: p_value = 1.0_dp
      integer :: degrees_of_freedom = 0
      integer :: info = 0
   end type

   type, public :: mts_mq_t
      ! Multivariate Ljung-Box statistics through each requested lag.
      real(dp), allocatable :: statistic(:), p_value(:)
      integer, allocatable :: degrees_of_freedom(:)
      integer :: adjustment = 0
      integer :: info = 0
   end type

   type, public :: mts_diagnostic_t
      ! Residual cross-correlations, significance p-values, and MQ statistics.
      real(dp), allocatable :: cross_correlation(:, :, :), p_value(:)
      type(mts_mq_t) :: mq
      integer :: info = 0
   end type

   type, public :: mts_var_backtest_t
      ! Rolling-origin VAR forecasts, errors, and horizon loss summaries.
      real(dp), allocatable :: forecast(:, :, :), error(:, :, :)
      real(dp), allocatable :: rmse(:, :), mean_absolute_error(:, :)
      integer :: origin = 0
      integer :: reestimate = 1
      integer :: info = 0
   end type

   type, public :: mts_scm_identification_t
      ! First-stage scalar-component zero-correlation and diagonal-difference tables.
      integer, allocatable :: zero_count(:, :), diagonal_difference(:, :)
      integer :: max_ar = 0
      integer :: max_ma = 0
      integer :: extra_lags = 0
      real(dp) :: significance = 0.05_dp
      integer :: info = 0
   end type

   type, public :: mts_scm_structure_t
      ! Detailed scalar-component orders and transformation matrix.
      integer, allocatable :: order(:, :)
      real(dp), allocatable :: transformation(:, :)
      integer :: components = 0
      integer :: info = 0
   end type

   type, public :: mts_scm_spec_t
      ! SCM transformation, AR, and MA parameter indicators.
      integer, allocatable :: transformation(:, :), ar(:, :, :), ma(:, :, :)
      integer, allocatable :: order(:, :), pivot(:)
      integer :: ar_order = 0
      integer :: ma_order = 0
      integer :: info = 0
   end type

   type, public :: mts_kronecker_identification_t
      ! Kronecker indices and sequential canonical-correlation test results.
      integer, allocatable :: index(:)
      real(dp), allocatable :: statistic(:, :), p_value(:, :)
      logical, allocatable :: tested(:, :)
      integer :: past_lag = 0
      real(dp) :: significance = 0.05_dp
      integer :: info = 0
   end type

   type, public :: mts_kronecker_spec_t
      ! Fixed-zero, fixed-one, and estimated structural coefficient indicators.
      integer, allocatable :: ar(:, :, :), ma(:, :, :)
      integer, allocatable :: index(:)
      integer :: order = 0
      integer :: info = 0
   end type

   type, public :: mts_kronecker_fit_t
      ! Constrained structural estimates and their reduced-form VARMA model.
      type(mts_varma_model_t) :: model
      type(mts_kronecker_spec_t) :: specification
      real(dp), allocatable :: structural_ar(:, :, :), structural_ma(:, :, :)
      real(dp), allocatable :: coefficients(:), standard_errors(:), covariance(:, :)
      real(dp), allocatable :: residuals(:, :)
      logical, allocatable :: estimated(:)
      real(dp) :: log_likelihood = 0.0_dp
      real(dp) :: aic = 0.0_dp
      real(dp) :: bic = 0.0_dp
      integer :: iterations = 0
      integer :: info = 0
      logical :: converged = .false.
      logical :: includes_mean = .true.
   end type

   type, public :: mts_kronecker_refinement_t
      ! Refined Kronecker fit and backward-elimination history.
      type(mts_kronecker_fit_t) :: fit
      integer, allocatable :: active_count(:), removed_index(:)
      integer :: steps = 0
      integer :: info = 0
   end type

   type, public :: mts_scm_fit_t
      ! Fitted SCM restrictions and shared structural VARMA estimates.
      type(mts_scm_spec_t) :: specification
      type(mts_kronecker_fit_t) :: fit
      integer :: info = 0
   end type

   type, public :: mts_scm_refinement_t
      ! Refined SCM fit and backward-elimination history.
      type(mts_scm_fit_t) :: fit
      integer, allocatable :: active_count(:), removed_index(:)
      integer :: steps = 0
      integer :: info = 0
   end type

   type, public :: mts_transfer_fit_t
      ! Dynamic-regression transfer coefficients, ARMA noise, and inference.
      real(dp), allocatable :: numerator(:), denominator(:), ar(:), ma(:)
      real(dp), allocatable :: coefficients(:), standard_errors(:), covariance(:, :)
      real(dp), allocatable :: filtered_input(:), noise(:), residuals(:)
      real(dp) :: intercept = 0.0_dp
      real(dp) :: sigma2 = 0.0_dp
      real(dp) :: log_likelihood = 0.0_dp
      real(dp) :: aic = 0.0_dp
      real(dp) :: bic = 0.0_dp
      integer :: delay = 0
      integer :: differences = 0
      integer :: iterations = 0
      integer :: info = 0
      logical :: converged = .false.
   end type

   type, public :: mts_transfer2_fit_t
      ! Two-input transfer-function estimates with seasonal ARMA disturbances.
      real(dp), allocatable :: numerator1(:), denominator1(:)
      real(dp), allocatable :: numerator2(:), denominator2(:)
      real(dp), allocatable :: ar(:), ma(:), seasonal_ar(:), seasonal_ma(:)
      real(dp), allocatable :: coefficients(:), standard_errors(:), covariance(:, :)
      real(dp), allocatable :: filtered_input1(:), filtered_input2(:), noise(:), residuals(:)
      real(dp) :: intercept = 0.0_dp
      real(dp) :: deterministic_coefficient = 0.0_dp
      real(dp) :: equilibrium_coefficient = 0.0_dp
      real(dp) :: sigma2 = 0.0_dp
      real(dp) :: log_likelihood = 0.0_dp
      real(dp) :: aic = 0.0_dp
      real(dp) :: bic = 0.0_dp
      integer :: delay1 = 0
      integer :: delay2 = 0
      integer :: differences = 0
      integer :: seasonal_differences = 0
      integer :: season = 1
      integer :: iterations = 0
      integer :: info = 0
      logical :: has_second_input = .false.
      logical :: has_deterministic = .false.
      logical :: has_equilibrium = .false.
      logical :: converged = .false.
   end type

   type, public :: mts_transfer_backtest_t
      ! Rolling one-step transfer-function forecast errors and loss summaries.
      real(dp), allocatable :: error(:)
      real(dp) :: bias = 0.0_dp
      real(dp) :: mse = 0.0_dp
      real(dp) :: rmse = 0.0_dp
      real(dp) :: mean_absolute_error = 0.0_dp
      integer :: origin = 0
      integer :: info = 0
   end type

   type, public :: mts_transfer_forecast_t
      ! Transfer-function point forecasts under zero future innovations.
      real(dp), allocatable :: mean(:)
      integer :: info = 0
   end type

   type, public :: mts_regts_fit_t
      ! Multivariate regression estimates with VAR error dynamics.
      real(dp), allocatable :: beta(:, :), ar(:, :, :), coefficients(:)
      real(dp), allocatable :: standard_errors(:), covariance(:, :), residuals(:, :), sigma(:, :)
      logical, allocatable :: estimated(:, :)
      real(dp) :: log_likelihood = 0.0_dp
      real(dp) :: aic = 0.0_dp
      real(dp) :: bic = 0.0_dp
      integer :: order = 0
      integer :: iterations = 0
      integer :: info = 0
      logical :: includes_mean = .true.
      logical :: converged = .false.
   end type

   type, public :: mts_regts_refinement_t
      ! Refined regression-with-VAR-errors fit and elimination history.
      type(mts_regts_fit_t) :: fit
      integer, allocatable :: active_count(:), removed_index(:)
      integer :: steps = 0
      integer :: info = 0
   end type

   type, public :: mts_ecm_known_fit_t
      ! Error-correction estimates for caller-supplied cointegrating processes.
      real(dp), allocatable :: loading(:, :), gamma(:, :, :), intercept(:)
      real(dp), allocatable :: coefficients(:, :), standard_errors(:, :)
      real(dp), allocatable :: residuals(:, :), sigma(:, :)
      logical, allocatable :: estimated(:, :)
      real(dp) :: aic = 0.0_dp
      real(dp) :: bic = 0.0_dp
      integer :: level_order = 0
      integer :: info = 0
      logical :: includes_constant = .false.
   end type

   type, public :: mts_ecm_normalized_fit_t
      ! Joint normalized cointegration, loading, and short-run ECM estimates.
      real(dp), allocatable :: loading(:, :), cointegration(:, :), gamma(:, :, :), intercept(:)
      real(dp), allocatable :: coefficients(:), standard_errors(:), covariance(:, :)
      real(dp), allocatable :: short_run(:, :), short_run_standard_errors(:, :)
      real(dp), allocatable :: residuals(:, :), sigma(:, :)
      logical, allocatable :: short_run_estimated(:, :)
      real(dp) :: log_likelihood = 0.0_dp
      real(dp) :: aic = 0.0_dp
      real(dp) :: bic = 0.0_dp
      integer :: rank = 0
      integer :: level_order = 0
      integer :: iterations = 0
      integer :: info = 0
      logical :: includes_constant = .false.
      logical :: converged = .false.
   end type

   type, public :: mts_reverse_mq_t
      ! Lag cross-correlations and reversed multivariate Ljung-Box tests.
      real(dp), allocatable :: cross_correlation(:, :, :)
      real(dp), allocatable :: cumulative(:), statistic(:), p_value(:)
      integer, allocatable :: degrees_of_freedom(:)
      integer :: max_lag = 0
      integer :: info = 0
   end type

   type, public :: mts_eccm_t
      ! Extended cross-correlation matrices and AR/MA order p-value table.
      real(dp), allocatable :: cross_correlation(:, :, :, :), p_value(:, :)
      integer :: max_ar = 0
      integer :: max_ma = 0
      integer :: info = 0
      logical :: reversed = .true.
   end type

   type, public :: mts_corner_t
      ! Transfer-function corner determinants and significance indicators.
      real(dp), allocatable :: value(:, :)
      logical, allocatable :: significant(:, :)
      real(dp), allocatable :: normalized_cross_correlation(:)
      real(dp) :: threshold = 0.0_dp
      integer :: info = 0
   end type

   type, public :: mts_vma_exact_fit_t
      ! Exact finite-sample Gaussian VMA estimates and inference.
      type(mts_varma_model_t) :: model
      real(dp), allocatable :: coefficients(:), standard_errors(:), covariance(:, :)
      real(dp), allocatable :: residuals(:, :)
      logical, allocatable :: estimated(:)
      real(dp) :: log_likelihood = 0.0_dp
      real(dp) :: aic = 0.0_dp
      real(dp) :: bic = 0.0_dp
      integer :: iterations = 0
      integer :: info = 0
      logical :: includes_mean = .true.
      logical :: converged = .false.
   end type

   type, public :: mts_vma_exact_refinement_t
      ! Refined exact VMA fit and backward-elimination history.
      type(mts_vma_exact_fit_t) :: fit
      integer, allocatable :: active_count(:), removed_index(:)
      integer :: steps = 0
      integer :: info = 0
   end type

   type, public :: mts_apca_t
      ! Asymptotic principal-component factors, loadings, and deviations.
      real(dp), allocatable :: factors(:, :), loadings(:, :), standard_deviation(:)
      integer :: components = 0
      integer :: info = 0
      logical :: transposed = .false.
   end type

   type, public :: mts_diffusion_forecast_t
      ! Stock-Watson diffusion-index regression and out-of-sample forecasts.
      real(dp), allocatable :: coefficients(:), forecast(:), loadings(:, :), index(:, :)
      real(dp), allocatable :: predictor_mean(:), predictor_scale(:)
      real(dp) :: mse = 0.0_dp
      integer :: origin = 0
      integer :: factors = 0
      integer :: info = 0
   end type

   type, public :: mts_multivariate_regression_t
      ! Multivariate least-squares coefficients, inference, and residual covariance.
      real(dp), allocatable :: coefficients(:, :), standard_errors(:, :)
      real(dp), allocatable :: covariance(:, :), residuals(:, :), sigma(:, :)
      integer :: residual_degrees_of_freedom = 0
      integer :: info = 0
      logical :: includes_constant = .true.
   end type

   type, public :: mts_var_chi_t
      ! Joint chi-square test of individually weak VAR coefficients.
      type(mts_var_fit_t) :: fit
      integer, allocatable :: coefficient_index(:)
      real(dp), allocatable :: values(:), covariance(:, :)
      real(dp) :: statistic = 0.0_dp
      real(dp) :: p_value = 1.0_dp
      real(dp) :: threshold = 1.645_dp
      integer :: degrees_of_freedom = 0
      integer :: info = 0
   end type

   type, public :: mts_vma_order_t
      ! VMA order criteria and selected AIC and BIC orders.
      real(dp), allocatable :: aic(:), bic(:)
      integer :: aic_order = 0
      integer :: bic_order = 0
      integer :: info = 0
   end type

   public :: mts_var, mts_vars, mts_var_psi, mts_var_forecast, mts_var_order
   public :: mts_var_irf, mts_fevd
   public :: mts_varx, mts_varx_forecast, mts_varx_irf, mts_varx_order
   public :: mts_varma_psi, mts_varma_pi, mts_varma_irf, mts_varma_covariance
   public :: mts_varma_simulate_from_innovations, mts_varma_simulate, mts_varma_forecast
   public :: mts_varma_fit, mts_varmas_fit
   public :: mts_svarma_expand, mts_svarma_fit
   public :: mts_refine_varma, mts_refine_svarma
   public :: mts_vecm_fit, mts_vecm_forecast
   public :: mts_factor_fit, mts_factor_forecast
   public :: mts_constrained_factor_fit, mts_constrained_factor_forecast
   public :: mts_bvar_fit, mts_minnesota_prior
   public :: mts_common_volatility
   public :: mts_mch_diagnostic
   public :: mts_bekk_fit, mts_bekk_forecast
   public :: mts_dcc_fit, mts_dcc_forecast
   public :: mts_adcc_fit, mts_adcc_t_fit, mts_adcc_forecast
   public :: mts_tse_tsui_fit, mts_tse_tsui_t_fit, mts_tse_tsui_forecast
   public :: mts_ewma_fit, mts_ewma_forecast, mts_mchol_fit, mts_mchol_forecast
   public :: mts_sccor, mts_arch_test, mts_march_test
   public :: mts_correlation_to_angles, mts_angles_to_correlation, mts_copula_fit
   public :: regularized_beta_mts
   public :: mts_var_missing, mts_var_partial_missing
   public :: mts_granger_test, mts_mq, mts_diagnostic, mts_var_backtest
   public :: mts_scm_identify
   public :: mts_scm_identify_details, mts_scm_specification
   public :: mts_scm_fit, mts_scm_refine, mts_scm_forecast
   public :: mts_transfer_fit
   public :: mts_transfer2_fit, mts_transfer2_forecast, mts_transfer2_backtest
   public :: mts_regts_fit, mts_regts_refine, mts_regts_forecast
   public :: mts_matrix_filter, mts_ecm_known_fit, mts_ecm_known_refine
   public :: mts_ecm_normalized_fit, mts_ecm_normalized_refine
   public :: mts_reverse_mq, mts_eccm, mts_corner
   public :: mts_vma_exact_fit, mts_vma_exact_refine
   public :: mts_apca, mts_diffusion_forecast
   public :: mts_multivariate_regression, mts_var_chi, mts_var_fore
   public :: mts_kronecker_identify, mts_kronecker_specification
   public :: mts_kronecker_fit, mts_kronecker_refine, mts_kronecker_forecast
   public :: mts_vma_fit, mts_vmas_fit, mts_vma_order

contains

   pure function diagonal_values(matrix) result(values)
      !! Return the main diagonal of a rectangular matrix.
      real(dp), intent(in) :: matrix(:, :) !! Input matrix.
      real(dp) :: values(min(size(matrix, 1), size(matrix, 2)))
      integer :: i

      do i = 1, size(values)
         values(i) = matrix(i, i)
      end do
   end function diagonal_values

   pure function mts_matrix_filter(data, weights, initial) result(filtered)
      !! Apply a recursive multivariate matrix-polynomial filter.
      real(dp), intent(in) :: data(:, :) !! Data.
      real(dp), intent(in) :: weights(:, :, :) !! Observation or objective weights.
      real(dp), intent(in), optional :: initial(:, :) !! Initial value.
      real(dp), allocatable :: filtered(:, :)
      real(dp), allocatable :: work(:, :)
      integer :: n, dimension, order, initial_count, t, lag, source

      n = size(data, 1)
      dimension = size(data, 2)
      order = size(weights, 3)
      initial_count = 0
      if (present(initial)) initial_count = size(initial, 1)
      allocate(work(initial_count + n, dimension), filtered(n, dimension))
      if (initial_count > 0) work(:initial_count, :) = initial
      do t = 1, n
         source = initial_count + t
         work(source, :) = data(t, :)
         do lag = 1, min(order, source - 1)
            work(source, :) = work(source, :) + &
               matmul(weights(:, :, lag), work(source - lag, :))
         end do
         filtered(t, :) = work(source, :)
      end do
   end function mts_matrix_filter

   pure function mts_ecm_known_fit(series, level_order, cointegrating_process, &
      include_constant, estimated) result(out)
      !! Fit an ECM by least squares for supplied cointegrating processes.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      real(dp), intent(in) :: cointegrating_process(:, :) !! Cointegrating process.
      integer, intent(in) :: level_order !! Level order.
      logical, intent(in), optional :: include_constant !! Whether to include the constant.
      logical, intent(in), optional :: estimated(:, :) !! Flag controlling estimated.
      type(mts_ecm_known_fit_t) :: out
      real(dp), allocatable :: differences(:, :), adjusted_process(:, :), design(:, :), response(:, :)
      real(dp), allocatable :: selected_design(:, :), xtx(:, :), inverse(:, :), beta(:)
      integer, allocatable :: selected(:)
      logical :: use_constant
      integer :: n, dimension, rank, regressors, rows, start, column, lag, equation, count_selected
      integer :: status, i, parameter_count
      real(dp) :: logdet
      real(dp), allocatable :: covariance_inverse(:, :)

      n = size(series, 1)
      dimension = size(series, 2)
      rank = size(cointegrating_process, 2)
      use_constant = .false.
      if (present(include_constant)) use_constant = include_constant
      regressors = rank + merge(1, 0, use_constant) + (level_order - 1)*dimension
      if (size(cointegrating_process, 1) /= n .or. level_order < 1 .or. rank < 1 .or. &
         n <= level_order + regressors + 1) then
         out%info = 1
         return
      end if
      allocate(differences(n, dimension))
      differences = 0.0_dp
      differences(2:, :) = series(2:, :) - series(:n - 1, :)
      adjusted_process = cointegrating_process
      if (.not. use_constant) adjusted_process = adjusted_process - &
         spread(sum(adjusted_process, 1)/real(n, dp), 1, n)
      start = max(1, level_order)
      rows = n - start + 1
      allocate(design(rows, regressors), response(rows, dimension))
      response = differences(start:, :)
      design(:, :rank) = adjusted_process(start - 1:n - 1, :)
      column = rank
      if (use_constant) then
         column = column + 1
         design(:, column) = 1.0_dp
      end if
      do lag = 1, level_order - 1
         design(:, column + 1:column + dimension) = differences(start - lag:n - lag, :)
         column = column + dimension
      end do
      allocate(out%coefficients(regressors, dimension), out%standard_errors(regressors, dimension))
      allocate(out%estimated(regressors, dimension), out%residuals(rows, dimension))
      out%coefficients = 0.0_dp
      out%standard_errors = 0.0_dp
      out%estimated = .true.
      if (present(estimated)) then
         if (any(shape(estimated) /= [regressors, dimension])) then
            out%info = 2
            return
         end if
         out%estimated = estimated
      end if
      do equation = 1, dimension
         count_selected = count(out%estimated(:, equation))
         out%residuals(:, equation) = response(:, equation)
         if (count_selected == 0) cycle
         selected = pack([(i, i=1, regressors)], out%estimated(:, equation))
         selected_design = design(:, selected)
         xtx = matmul(transpose(selected_design), selected_design)
         call invert_matrix(xtx, inverse, status)
         if (status /= 0) then
            out%info = 10 + equation
            return
         end if
         beta = matmul(inverse, matmul(transpose(selected_design), response(:, equation)))
         out%coefficients(selected, equation) = beta
         out%residuals(:, equation) = response(:, equation) - matmul(selected_design, beta)
         out%standard_errors(selected, equation) = sqrt(max(0.0_dp, diagonal_values(inverse)* &
            sum(out%residuals(:, equation)**2)/real(rows, dp)))
      end do
      out%sigma = matmul(transpose(out%residuals), out%residuals)/real(rows, dp)
      allocate(covariance_inverse(dimension, dimension))
      call inverse_logdet(out%sigma, covariance_inverse, logdet, status, 100.0_dp*epsilon(1.0_dp))
      if (status /= 0) then
         out%info = 20 + status
         return
      end if
      parameter_count = count(out%estimated)
      out%aic = logdet + 2.0_dp*real(parameter_count, dp)/real(n, dp)
      out%bic = logdet + log(real(n, dp))*real(parameter_count, dp)/real(n, dp)
      out%loading = transpose(out%coefficients(:rank, :))
      allocate(out%intercept(dimension), out%gamma(dimension, dimension, level_order - 1))
      out%intercept = 0.0_dp
      column = rank
      if (use_constant) then
         column = column + 1
         out%intercept = out%coefficients(column, :)
      end if
      do lag = 1, level_order - 1
         out%gamma(:, :, lag) = transpose(out%coefficients(column + 1:column + dimension, :))
         column = column + dimension
      end do
      out%level_order = level_order
      out%includes_constant = use_constant
   end function mts_ecm_known_fit

   pure function mts_ecm_known_refine(series, cointegrating_process, fit, threshold) result(out)
      !! Refit a known-process ECM after thresholding coefficient t ratios.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      real(dp), intent(in) :: cointegrating_process(:, :) !! Cointegrating process.
      type(mts_ecm_known_fit_t), intent(in) :: fit !! Previously fitted model.
      real(dp), intent(in), optional :: threshold !! Decision or truncation threshold.
      type(mts_ecm_known_fit_t) :: out
      logical, allocatable :: mask(:, :)
      real(dp) :: cutoff

      cutoff = 1.0_dp
      if (present(threshold)) cutoff = threshold
      if (fit%info /= 0 .or. cutoff < 0.0_dp) then
         out%info = 1
         return
      end if
      mask = fit%estimated .and. fit%standard_errors > 0.0_dp .and. &
         abs(fit%coefficients) > cutoff*fit%standard_errors
      out = mts_ecm_known_fit(series, fit%level_order, cointegrating_process, &
         fit%includes_constant, mask)
   end function mts_ecm_known_refine

   pure function mts_ecm_normalized_fit(series, level_order, initial_cointegration, &
      include_constant, short_run_estimated, initial, max_iterations, tolerance) result(out)
      !! Jointly fit an identity-normalized reduced-rank ECM.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      real(dp), intent(in) :: initial_cointegration(:, :) !! Initial cointegration.
      integer, intent(in) :: level_order !! Level order.
      logical, intent(in), optional :: include_constant !! Whether to include the constant.
      logical, intent(in), optional :: short_run_estimated(:, :) !! Flag controlling short run estimated.
      real(dp), intent(in), optional :: initial(:) !! Initial value.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      type(mts_ecm_normalized_fit_t) :: out
      type(mts_ecm_known_fit_t) :: preliminary
      type(optimization_result_t) :: optimization
      real(dp), allocatable :: process(:, :), parameters(:), hessian(:, :), inverse(:, :)
      real(dp), allocatable :: free_values(:)
      integer, allocatable :: free_index(:)
      logical, allocatable :: mask(:), short_mask(:, :)
      logical :: use_constant
      integer :: n, dimension, rank, short_rows, alpha_count, beta_count, total_count
      integer :: free_count, offset, equation, i, status, limit
      real(dp) :: gradient_tolerance

      n = size(series, 1)
      dimension = size(series, 2)
      rank = size(initial_cointegration, 2)
      use_constant = .false.
      if (present(include_constant)) use_constant = include_constant
      short_rows = merge(1, 0, use_constant) + (level_order - 1)*dimension
      alpha_count = dimension*rank
      beta_count = (dimension - rank)*rank
      if (size(initial_cointegration, 1) /= dimension .or. rank < 1 .or. rank >= dimension .or. &
         level_order < 1 .or. n <= level_order + dimension + 2) then
         out%info = 1
         return
      end if
      allocate(short_mask(short_rows, dimension))
      short_mask = .true.
      if (present(short_run_estimated)) then
         if (any(shape(short_run_estimated) /= [short_rows, dimension])) then
            out%info = 2
            return
         end if
         short_mask = short_run_estimated
      end if
      process = matmul(series, initial_cointegration)
      preliminary = mts_ecm_known_fit(series, level_order, process, use_constant)
      if (preliminary%info /= 0) then
         out%info = 10 + preliminary%info
         return
      end if
      total_count = alpha_count + beta_count + short_rows*dimension
      allocate(parameters(total_count), mask(total_count))
      parameters = 0.0_dp
      parameters(:alpha_count) = reshape(preliminary%loading, [alpha_count])
      parameters(alpha_count + 1:alpha_count + beta_count) = &
         reshape(initial_cointegration(rank + 1:, :), [beta_count])
      if (short_rows > 0) parameters(alpha_count + beta_count + 1:) = &
         reshape(preliminary%coefficients(rank + 1:, :), [short_rows*dimension])
      if (present(initial)) then
         if (size(initial) /= total_count) then
            out%info = 2
            return
         end if
         parameters = initial
      end if
      mask = .true.
      if (short_rows > 0) mask(alpha_count + beta_count + 1:) = &
         reshape(short_mask, [short_rows*dimension])
      free_count = count(mask)
      allocate(free_index(free_count), free_values(free_count))
      free_index = pack([(i, i=1, total_count)], mask)
      free_values = pack(parameters, mask)
      limit = 200
      gradient_tolerance = 1.0e-6_dp
      if (present(max_iterations)) limit = max_iterations
      if (present(tolerance)) gradient_tolerance = tolerance
      optimization = bfgs_minimize_fd(objective, free_values, limit, gradient_tolerance)
      parameters(free_index) = optimization%parameters
      out%iterations = optimization%iterations
      out%converged = optimization%converged
      if (optimization%info /= 0) out%info = 100 + optimization%info
      call evaluate(parameters, out%loading, out%cointegration, out%short_run, &
         out%residuals, out%sigma, out%log_likelihood, status)
      if (status /= 0) then
         out%info = 20 + status
         return
      end if
      out%coefficients = parameters
      out%rank = rank
      out%level_order = level_order
      out%includes_constant = use_constant
      out%short_run_estimated = short_mask
      out%aic = -2.0_dp*out%log_likelihood + 2.0_dp*real(free_count + 1, dp)
      out%bic = -2.0_dp*out%log_likelihood + log(real(n - level_order + 1, dp))* &
         real(free_count + 1, dp)
      allocate(out%covariance(total_count, total_count), out%standard_errors(total_count))
      out%covariance = 0.0_dp
      out%standard_errors = 0.0_dp
      hessian = finite_difference_hessian(objective, pack(parameters, mask))
      call invert_matrix(hessian, inverse, status)
      if (status == 0) then
         do i = 1, free_count
            out%covariance(free_index, free_index(i)) = inverse(:, i)
            out%standard_errors(free_index(i)) = sqrt(max(0.0_dp, inverse(i, i)))
         end do
      else if (out%info == 0) then
         out%info = 200 + status
      end if
      allocate(out%short_run_standard_errors(short_rows, dimension))
      out%short_run_standard_errors = 0.0_dp
      if (short_rows > 0) out%short_run_standard_errors = reshape(&
         out%standard_errors(alpha_count + beta_count + 1:), [short_rows, dimension])
      allocate(out%intercept(dimension), out%gamma(dimension, dimension, level_order - 1))
      out%intercept = 0.0_dp
      offset = 0
      if (use_constant) then
         offset = 1
         out%intercept = out%short_run(1, :)
      end if
      do i = 1, level_order - 1
         out%gamma(:, :, i) = transpose(out%short_run(offset + 1:offset + dimension, :))
         offset = offset + dimension
      end do

   contains

      pure function objective(free_parameters) result(value)
         !! Return the concentrated negative likelihood for normalized ECM parameters.
         real(dp), intent(in) :: free_parameters(:) !! Free parameters.
         real(dp) :: value, likelihood
         real(dp) :: full_parameters(total_count)
         real(dp), allocatable :: loading(:, :), cointegration(:, :), short_run(:, :)
         real(dp), allocatable :: residual_values(:, :), sigma(:, :)
         integer :: evaluation_status

         full_parameters = parameters
         full_parameters(free_index) = free_parameters
         call evaluate(full_parameters, loading, cointegration, short_run, residual_values, &
            sigma, likelihood, evaluation_status)
         if (evaluation_status == 0 .and. ieee_is_finite(likelihood)) then
            value = -likelihood
         else
            value = 1.0e30_dp + dot_product(free_parameters, free_parameters)
         end if
      end function objective

      pure subroutine evaluate(values, loading, cointegration, short_run, residual_values, &
         sigma, likelihood, evaluation_status)
         !! Unpack normalized ECM parameters and evaluate innovations.
         real(dp), intent(in) :: values(:) !! Input values.
         real(dp), allocatable, intent(out) :: loading(:, :) !! Loading.
         real(dp), allocatable, intent(out) :: cointegration(:, :) !! Cointegration.
         real(dp), allocatable, intent(out) :: short_run(:, :) !! Short run.
         real(dp), allocatable, intent(out) :: residual_values(:, :) !! Residual values.
         real(dp), allocatable, intent(out) :: sigma(:, :) !! Scale parameter or standard deviation.
         real(dp), intent(out) :: likelihood !! Likelihood.
         integer, intent(out) :: evaluation_status !! Evaluation status.
         real(dp), allocatable :: differences(:, :), covariance_inverse(:, :)
         real(dp) :: logdet
         integer :: rows, start, t, lag, column

         loading = reshape(values(:alpha_count), [dimension, rank])
         allocate(cointegration(dimension, rank))
         cointegration = 0.0_dp
         do t = 1, rank
            cointegration(t, t) = 1.0_dp
         end do
         cointegration(rank + 1:, :) = reshape(&
            values(alpha_count + 1:alpha_count + beta_count), [dimension - rank, rank])
         allocate(short_run(short_rows, dimension))
         if (short_rows > 0) short_run = reshape(values(alpha_count + beta_count + 1:), &
            [short_rows, dimension])
         allocate(differences(n, dimension))
         differences = 0.0_dp
         differences(2:, :) = series(2:, :) - series(:n - 1, :)
         start = max(2, level_order)
         rows = n - start + 1
         allocate(residual_values(rows, dimension))
         do t = start, n
            residual_values(t - start + 1, :) = differences(t, :) - &
               matmul(loading, matmul(transpose(cointegration), series(t - 1, :)))
            column = 0
            if (use_constant) then
               column = 1
               residual_values(t - start + 1, :) = residual_values(t - start + 1, :) - &
                  short_run(1, :)
            end if
            do lag = 1, level_order - 1
               residual_values(t - start + 1, :) = residual_values(t - start + 1, :) - &
                  matmul(transpose(short_run(column + 1:column + dimension, :)), &
                  differences(t - lag, :))
               column = column + dimension
            end do
         end do
         sigma = matmul(transpose(residual_values), residual_values)/real(rows, dp)
         allocate(covariance_inverse(dimension, dimension))
         call inverse_logdet(sigma, covariance_inverse, logdet, evaluation_status, &
            100.0_dp*epsilon(1.0_dp))
         if (evaluation_status /= 0) then
            likelihood = -huge(1.0_dp)
            return
         end if
         likelihood = -0.5_dp*real(rows, dp)*(real(dimension, dp)* &
            (log(2.0_dp*acos(-1.0_dp)) + 1.0_dp) + logdet)
      end subroutine evaluate
   end function mts_ecm_normalized_fit

   pure function mts_ecm_normalized_refine(series, fit, threshold, &
      max_iterations, tolerance) result(out)
      !! Refit a normalized ECM after thresholding short-run parameters.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      type(mts_ecm_normalized_fit_t), intent(in) :: fit !! Previously fitted model.
      real(dp), intent(in), optional :: threshold !! Decision or truncation threshold.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      type(mts_ecm_normalized_fit_t) :: out
      logical, allocatable :: mask(:, :)
      real(dp) :: cutoff

      cutoff = 1.0_dp
      if (present(threshold)) cutoff = threshold
      if (fit%info /= 0 .or. cutoff < 0.0_dp) then
         out%info = 1
         return
      end if
      mask = fit%short_run_estimated .and. fit%short_run_standard_errors > 0.0_dp .and. &
         abs(fit%short_run) > cutoff*fit%short_run_standard_errors
      out = mts_ecm_normalized_fit(series, fit%level_order, fit%cointegration, &
         fit%includes_constant, mask, fit%coefficients, max_iterations, tolerance)
   end function mts_ecm_normalized_refine

   pure function mts_reverse_mq(series, max_lag) result(out)
      !! Compute reversed multivariate Ljung-Box tests through a maximum lag.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      integer, intent(in) :: max_lag !! Maximum lag to consider.
      type(mts_reverse_mq_t) :: out
      real(dp), allocatable :: centered(:, :), covariance(:, :), inverse(:, :), lag_covariance(:, :)
      real(dp), allocatable :: scale(:), product(:, :)
      integer :: n, dimension, lag, status

      n = size(series, 1)
      dimension = size(series, 2)
      if (dimension < 1 .or. max_lag < 1 .or. n <= max_lag + dimension) then
         out%info = 1
         return
      end if
      centered = series - spread(sum(series, 1)/real(n, dp), 1, n)
      covariance = matmul(transpose(centered), centered)/real(n - 1, dp)
      call invert_matrix(covariance, inverse, status)
      if (status /= 0 .or. any(diagonal_values(covariance) <= 0.0_dp)) then
         out%info = 2
         return
      end if
      scale = sqrt(diagonal_values(covariance))
      allocate(out%cross_correlation(dimension, dimension, max_lag))
      allocate(out%cumulative(max_lag), out%statistic(max_lag), out%p_value(max_lag))
      allocate(out%degrees_of_freedom(max_lag))
      out%cumulative = 0.0_dp
      do lag = 1, max_lag
         lag_covariance = matmul(transpose(centered(lag + 1:, :)), centered(:n - lag, :))/ &
            real(n - 1, dp)
         out%cross_correlation(:, :, lag) = lag_covariance/ &
            spread(scale, 2, dimension)/spread(scale, 1, dimension)
         product = matmul(transpose(lag_covariance), matmul(inverse, &
            matmul(lag_covariance, inverse)))
         if (lag > 1) out%cumulative(lag) = out%cumulative(lag - 1)
         out%cumulative(lag) = out%cumulative(lag) + &
            real(n*n, dp)*sum(diagonal_values(product))/real(n - lag, dp)
      end do
      do lag = 1, max_lag
         if (lag == 1) then
            out%statistic(lag) = out%cumulative(max_lag)
         else
            out%statistic(lag) = out%cumulative(max_lag) - out%cumulative(lag - 1)
         end if
         out%degrees_of_freedom(lag) = dimension*dimension*(max_lag - lag + 1)
         out%p_value(lag) = chi_square_survival(out%statistic(lag), out%degrees_of_freedom(lag))
      end do
      out%max_lag = max_lag
   end function mts_reverse_mq

   pure function mts_eccm(series, max_ar, max_ma, include_mean, reversed) result(out)
      !! Compute iterated-regression extended cross-correlation order tables.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      integer, intent(in) :: max_ar !! Maximum autoregressive.
      integer, intent(in) :: max_ma !! Maximum moving-average.
      logical, intent(in), optional :: include_mean !! Whether to include a mean term.
      logical, intent(in), optional :: reversed !! Flag controlling reversed.
      type(mts_eccm_t) :: out
      type(mts_var_fit_t) :: var_fit
      type(mts_reverse_mq_t) :: reverse_test
      type(mts_mq_t) :: forward_test
      real(dp), allocatable :: working(:, :), previous_residual(:, :), innovations(:, :)
      real(dp), allocatable :: response(:, :), design(:, :), xtx(:, :), inverse(:, :), beta(:, :), wt(:, :)
      logical :: center, use_reverse
      integer :: n, dimension, p, q, start, rows, columns, lag, status

      n = size(series, 1)
      dimension = size(series, 2)
      center = .false.
      use_reverse = .true.
      if (present(include_mean)) center = include_mean
      if (present(reversed)) use_reverse = reversed
      if (dimension < 1 .or. max_ar < 0 .or. max_ma < 0 .or. &
         n <= max_ar + max_ma + dimension + 3) then
         out%info = 1
         return
      end if
      working = series
      if (center) working = working - spread(sum(working, 1)/real(n, dp), 1, n)
      allocate(out%cross_correlation(dimension, dimension, 0:max_ma, 0:max_ar))
      allocate(out%p_value(0:max_ar, 0:max_ma))
      out%cross_correlation = 0.0_dp
      out%p_value = 1.0_dp
      do p = 0, max_ar
         allocate(previous_residual(n, dimension))
         previous_residual = 0.0_dp
         if (p == 0) then
            previous_residual = working
         else
            var_fit = mts_var(working, p, .false.)
            if (var_fit%info /= 0) then
               out%info = 10 + p
               return
            end if
            previous_residual(p + 1:, :) = var_fit%residuals
         end if
         reverse_test = mts_reverse_mq(previous_residual(p + 1:, :), max_ma + 1)
         if (reverse_test%info /= 0) then
            out%info = 20 + p
            return
         end if
         out%cross_correlation(:, :, 0, p) = reverse_test%cross_correlation(:, :, 1)
         if (use_reverse) then
            out%p_value(p, 0) = reverse_test%p_value(1)
         else
            forward_test = mts_mq(previous_residual(p + 1:, :), 1)
            out%p_value(p, 0) = forward_test%p_value(1)
         end if
         do q = 1, max_ma
            start = p + q + 1
            rows = n - start + 1
            columns = (p + q)*dimension
            allocate(response(rows, dimension), design(rows, columns))
            response = working(start:, :)
            do lag = 1, p
               design(:, (lag - 1)*dimension + 1:lag*dimension) = &
                  working(start - lag:n - lag, :)
            end do
            do lag = 1, q
               design(:, (p + lag - 1)*dimension + 1:(p + lag)*dimension) = &
                  previous_residual(start - lag:n - lag, :)
            end do
            xtx = matmul(transpose(design), design)
            call invert_matrix(xtx, inverse, status)
            if (status /= 0) then
               out%info = 30 + p*(max_ma + 1) + q
               return
            end if
            beta = matmul(inverse, matmul(transpose(design), response))
            if (p > 0) then
               wt = response - matmul(design(:, :p*dimension), beta(:p*dimension, :))
            else
               wt = response
            end if
            allocate(innovations(n, dimension))
            innovations = 0.0_dp
            innovations(start:, :) = response - matmul(design, beta)
            previous_residual = innovations
            reverse_test = mts_reverse_mq(wt, max_ma + 1)
            if (reverse_test%info /= 0) then
               out%info = 40 + p*(max_ma + 1) + q
               return
            end if
            out%cross_correlation(:, :, q, p) = reverse_test%cross_correlation(:, :, q + 1)
            if (use_reverse) then
               out%p_value(p, q) = reverse_test%p_value(q + 1)
            else
               forward_test = mts_mq(wt, q + 1)
               out%p_value(p, q) = forward_test%p_value(q + 1)
            end if
            deallocate(response, design, xtx, inverse, beta, wt, innovations)
         end do
         deallocate(previous_residual)
      end do
      out%max_ar = max_ar
      out%max_ma = max_ma
      out%reversed = use_reverse
   end function mts_eccm

   pure function mts_corner(output_series, input_series, row_count, column_count) result(out)
      !! Compute transfer-function corner determinants from prewhitened series.
      real(dp), intent(in) :: output_series(:) !! Output series.
      real(dp), intent(in) :: input_series(:) !! Input series.
      integer, intent(in) :: row_count !! Number of row.
      integer, intent(in) :: column_count !! Number of column.
      type(mts_corner_t) :: out
      real(dp), allocatable :: y(:), x(:), matrix(:, :)
      real(dp) :: sy, sx, maximum
      integer :: n, lag_count, observations, lag, row, column, i, j, index

      n = min(size(output_series), size(input_series))
      lag_count = row_count + column_count + 1
      if (row_count < 1 .or. column_count < 1 .or. n <= lag_count + 2) then
         out%info = 1
         return
      end if
      y = output_series(:n) - sum(output_series(:n))/real(n, dp)
      x = input_series(:n) - sum(input_series(:n))/real(n, dp)
      sy = sqrt(sum(y**2)/real(n - 1, dp))
      sx = sqrt(sum(x**2)/real(n - 1, dp))
      if (sx <= tiny(1.0_dp) .or. sy <= tiny(1.0_dp)) then
         out%info = 2
         return
      end if
      observations = n - lag_count + 1
      allocate(out%normalized_cross_correlation(lag_count))
      do lag = 0, lag_count - 1
         out%normalized_cross_correlation(lag + 1) = &
            dot_product(y(lag_count:n), x(lag_count - lag:n - lag))/ &
            real(observations - 1, dp)/sx**2
      end do
      maximum = maxval(abs(out%normalized_cross_correlation))
      if (maximum <= tiny(1.0_dp)) then
         out%info = 2
         return
      end if
      out%normalized_cross_correlation = out%normalized_cross_correlation/maximum
      allocate(out%value(row_count, column_count), out%significant(row_count, column_count))
      out%value(:, 1) = out%normalized_cross_correlation(:row_count)
      do column = 2, column_count
         allocate(matrix(column, column))
         do row = 1, row_count
            matrix = 0.0_dp
            do i = 1, column
               matrix(i, i) = out%normalized_cross_correlation(row)
            end do
            do i = 2, column
               do j = 1, i - 1
                  matrix(i, j) = out%normalized_cross_correlation(row + j)
               end do
            end do
            do j = 2, column
               do i = 1, j - 1
                  index = row - j + 1
                  if (index > 0) matrix(i, j) = out%normalized_cross_correlation(index)
               end do
            end do
            out%value(row, column) = matrix_determinant(matrix)
         end do
         deallocate(matrix)
      end do
      out%threshold = 2.0_dp/sqrt(real(n, dp))
      out%significant = abs(out%value) > out%threshold
   end function mts_corner

   pure function matrix_determinant(matrix) result(value)
      !! Compute a square-matrix determinant by pivoted elimination.
      real(dp), intent(in) :: matrix(:, :) !! Input matrix.
      real(dp) :: value
      real(dp), allocatable :: work(:, :), row_values(:)
      real(dp) :: pivot
      integer :: n, column, pivot_row, row

      n = size(matrix, 1)
      work = matrix
      value = 1.0_dp
      do column = 1, n
         pivot_row = column - 1 + maxloc(abs(work(column:, column)), dim=1)
         pivot = work(pivot_row, column)
         if (abs(pivot) <= tiny(1.0_dp)) then
            value = 0.0_dp
            return
         end if
         if (pivot_row /= column) then
            row_values = work(column, :)
            work(column, :) = work(pivot_row, :)
            work(pivot_row, :) = row_values
            value = -value
         end if
         value = value*work(column, column)
         do row = column + 1, n
            work(row, column + 1:) = work(row, column + 1:) - &
               work(row, column)/work(column, column)*work(column, column + 1:)
         end do
      end do
   end function matrix_determinant

   pure function mts_vma_exact_fit(series, order, include_mean, initial, estimated, &
      max_iterations, tolerance) result(out)
      !! Fit a VMA model by its exact finite-sample Gaussian likelihood.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      integer, intent(in) :: order !! Model or polynomial order.
      logical, intent(in), optional :: include_mean !! Whether to include a mean term.
      logical, intent(in), optional :: estimated(:) !! Flag controlling estimated.
      real(dp), intent(in), optional :: initial(:) !! Initial value.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      type(mts_vma_exact_fit_t) :: out
      type(optimization_result_t) :: optimization
      real(dp), allocatable :: parameters(:), free_values(:), hessian(:, :), inverse(:, :)
      integer, allocatable :: free_index(:)
      logical :: use_mean
      integer :: n, dimension, block_size, parameter_count, free_count, i, status, limit
      real(dp) :: gradient_tolerance

      n = size(series, 1)
      dimension = size(series, 2)
      use_mean = .true.
      if (present(include_mean)) use_mean = include_mean
      block_size = dimension*dimension
      parameter_count = order*block_size + merge(dimension, 0, use_mean)
      if (dimension < 1 .or. order < 1 .or. n <= order + dimension + 1) then
         out%info = 1
         return
      end if
      allocate(parameters(parameter_count), out%estimated(parameter_count))
      parameters = 0.0_dp
      if (use_mean) parameters(:dimension) = sum(series, 1)/real(n, dp)
      if (present(initial)) then
         if (size(initial) /= parameter_count) then
            out%info = 2
            return
         end if
         parameters = initial
      end if
      out%estimated = .true.
      if (present(estimated)) then
         if (size(estimated) /= parameter_count) then
            out%info = 2
            return
         end if
         out%estimated = estimated
      end if
      free_count = count(out%estimated)
      allocate(free_index(free_count), free_values(free_count))
      free_index = pack([(i, i=1, parameter_count)], out%estimated)
      free_values = pack(parameters, out%estimated)
      limit = 120
      gradient_tolerance = 1.0e-6_dp
      if (present(max_iterations)) limit = max_iterations
      if (present(tolerance)) gradient_tolerance = tolerance
      if (free_count > 0) then
         optimization = bfgs_minimize_fd(objective, free_values, limit, gradient_tolerance)
         parameters(free_index) = optimization%parameters
         out%iterations = optimization%iterations
         out%converged = optimization%converged
         if (optimization%info /= 0) out%info = 100 + optimization%info
      else
         out%converged = .true.
      end if
      call evaluate(parameters, out%model, out%residuals, out%log_likelihood, status)
      if (status /= 0) then
         out%info = 10 + status
         return
      end if
      out%coefficients = parameters
      out%includes_mean = use_mean
      out%aic = -2.0_dp*out%log_likelihood + 2.0_dp*real(free_count + 1, dp)
      out%bic = -2.0_dp*out%log_likelihood + log(real(n, dp))*real(free_count + 1, dp)
      allocate(out%covariance(parameter_count, parameter_count), out%standard_errors(parameter_count))
      out%covariance = 0.0_dp
      out%standard_errors = 0.0_dp
      if (free_count > 0) then
         hessian = finite_difference_hessian(objective, pack(parameters, out%estimated))
         call invert_matrix(hessian, inverse, status)
         if (status == 0) then
            do i = 1, free_count
               out%covariance(free_index, free_index(i)) = inverse(:, i)
               out%standard_errors(free_index(i)) = sqrt(max(0.0_dp, inverse(i, i)))
            end do
         else if (out%info == 0) then
            out%info = 200 + status
         end if
      end if

   contains

      pure function objective(free_parameters) result(value)
         !! Return the exact negative Gaussian log likelihood.
         real(dp), intent(in) :: free_parameters(:) !! Free parameters.
         real(dp) :: value, likelihood
         real(dp) :: full_parameters(parameter_count)
         real(dp), allocatable :: residual_values(:, :)
         type(mts_varma_model_t) :: model
         integer :: evaluation_status

         full_parameters = parameters
         full_parameters(free_index) = free_parameters
         call evaluate(full_parameters, model, residual_values, likelihood, evaluation_status)
         if (evaluation_status == 0 .and. ieee_is_finite(likelihood)) then
            value = -likelihood
         else
            value = 1.0e30_dp + dot_product(free_parameters, free_parameters)
         end if
      end function objective

      pure subroutine evaluate(values, model, residual_values, likelihood, evaluation_status)
         !! Build the block covariance and evaluate the finite-sample likelihood.
         real(dp), intent(in) :: values(:) !! Input values.
         type(mts_varma_model_t), intent(out) :: model !! Model specification.
         real(dp), allocatable, intent(out) :: residual_values(:, :) !! Residual values.
         real(dp), intent(out) :: likelihood !! Likelihood.
         integer, intent(out) :: evaluation_status !! Evaluation status.
         real(dp), allocatable :: centered(:, :), coefficient(:, :, :), gamma(:, :, :)
         real(dp), allocatable :: full_covariance(:, :), full_inverse(:, :), vector(:)
         real(dp) :: logdet, quadratic
         integer :: offset, lag, t, j, row_start, column_start

         allocate(model%ar(dimension, dimension, 0), model%ma(dimension, dimension, order))
         allocate(model%intercept(dimension), model%sigma(dimension, dimension))
         model%intercept = 0.0_dp
         offset = 0
         if (use_mean) then
            model%intercept = values(:dimension)
            offset = dimension
         end if
         do lag = 1, order
            model%ma(:, :, lag) = reshape(values(offset + 1:offset + block_size), &
               [dimension, dimension])
            offset = offset + block_size
         end do
         if (companion_radius(model%ma) >= 1.0_dp) then
            evaluation_status = 2
            likelihood = -huge(1.0_dp)
            return
         end if
         centered = series - spread(model%intercept, 1, n)
         allocate(residual_values(n, dimension))
         do t = 1, n
            residual_values(t, :) = centered(t, :)
            do lag = 1, min(order, t - 1)
               residual_values(t, :) = residual_values(t, :) + &
                  matmul(model%ma(:, :, lag), residual_values(t - lag, :))
            end do
         end do
         model%sigma = matmul(transpose(residual_values), residual_values)/real(n, dp)
         allocate(coefficient(dimension, dimension, 0:order), gamma(dimension, dimension, 0:order))
         coefficient = 0.0_dp
         coefficient(:, :, 0) = identity_matrix(dimension)
         coefficient(:, :, 1:) = -model%ma
         do lag = 0, order
            gamma(:, :, lag) = 0.0_dp
            do j = 0, order - lag
               gamma(:, :, lag) = gamma(:, :, lag) + matmul(coefficient(:, :, j + lag), &
                  matmul(model%sigma, transpose(coefficient(:, :, j))))
            end do
         end do
         allocate(full_covariance(n*dimension, n*dimension), full_inverse(n*dimension, n*dimension))
         full_covariance = 0.0_dp
         do t = 1, n
            row_start = (t - 1)*dimension + 1
            do j = 1, n
               column_start = (j - 1)*dimension + 1
               lag = abs(t - j)
               if (lag > order) cycle
               if (t >= j) then
                  full_covariance(row_start:row_start + dimension - 1, &
                     column_start:column_start + dimension - 1) = gamma(:, :, lag)
               else
                  full_covariance(row_start:row_start + dimension - 1, &
                     column_start:column_start + dimension - 1) = transpose(gamma(:, :, lag))
               end if
            end do
         end do
         call inverse_logdet(full_covariance, full_inverse, logdet, evaluation_status, &
            1000.0_dp*epsilon(1.0_dp))
         if (evaluation_status /= 0) then
            likelihood = -huge(1.0_dp)
            return
         end if
         vector = reshape(transpose(centered), [n*dimension])
         quadratic = dot_product(vector, matmul(full_inverse, vector))
         likelihood = -0.5_dp*(real(n*dimension, dp)*log(2.0_dp*acos(-1.0_dp)) + &
            logdet + quadratic)
      end subroutine evaluate
   end function mts_vma_exact_fit

   pure function mts_vma_exact_refine(series, fit, threshold, max_steps, &
      max_iterations, tolerance) result(out)
      !! Refine exact VMA estimates by backward t-ratio elimination.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      type(mts_vma_exact_fit_t), intent(in) :: fit !! Previously fitted model.
      real(dp), intent(in), optional :: threshold !! Decision or truncation threshold.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      integer, intent(in), optional :: max_steps !! Maximum steps.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      type(mts_vma_exact_refinement_t) :: out
      logical, allocatable :: mask(:)
      real(dp) :: cutoff, ratio, smallest
      integer :: limit, iterations_limit, step, weakest, i

      cutoff = 1.0_dp
      limit = size(fit%coefficients)
      iterations_limit = 120
      if (present(threshold)) cutoff = threshold
      if (present(max_steps)) limit = min(limit, max_steps)
      if (present(max_iterations)) iterations_limit = max_iterations
      if (fit%info /= 0 .or. cutoff < 0.0_dp) then
         out%info = 1
         return
      end if
      mask = fit%estimated
      out%fit = fit
      allocate(out%active_count(0:limit), out%removed_index(limit))
      out%active_count(0) = count(mask)
      do step = 1, limit
         weakest = 0
         smallest = huge(1.0_dp)
         do i = 1, size(mask)
            if (.not. mask(i) .or. out%fit%standard_errors(i) <= 0.0_dp) cycle
            ratio = abs(out%fit%coefficients(i)/out%fit%standard_errors(i))
            if (ratio < cutoff .and. ratio < smallest) then
               weakest = i
               smallest = ratio
            end if
         end do
         if (weakest == 0) exit
         mask(weakest) = .false.
         out%fit = mts_vma_exact_fit(series, size(fit%model%ma, 3), fit%includes_mean, &
            out%fit%coefficients, mask, iterations_limit, tolerance)
         if (out%fit%info /= 0 .and. out%fit%info < 100) then
            out%info = 10 + step
            return
         end if
         out%removed_index(step) = weakest
         out%active_count(step) = count(mask)
         out%steps = step
      end do
      out%active_count = out%active_count(:out%steps)
      out%removed_index = out%removed_index(:out%steps)
   end function mts_vma_exact_refine

   pure function mts_apca(data, component_count) result(out)
      !! Perform asymptotic PCA by interchanging observations and variables.
      real(dp), intent(in) :: data(:, :) !! Data.
      integer, intent(in) :: component_count !! Number of component.
      type(mts_apca_t) :: out
      real(dp), allocatable :: work(:, :), centered(:, :), covariance(:, :)
      real(dp), allocatable :: eigenvalues(:), eigenvectors(:, :)
      integer :: observations, variables, components, status

      if (size(data, 1) < 2 .or. size(data, 2) < 2 .or. component_count < 1) then
         out%info = 1
         return
      end if
      work = data
      if (size(data, 2) <= size(data, 1)) then
         work = transpose(data)
         out%transposed = .true.
      end if
      observations = size(work, 1)
      variables = size(work, 2)
      components = min(component_count, observations)
      centered = work - spread(sum(work, 2)/real(variables, dp), 2, variables)
      covariance = matmul(centered, transpose(centered))/real(variables - 1, dp)
      call symmetric_eigen(covariance, eigenvalues, eigenvectors, status)
      if (status /= 0) then
         out%info = 10 + status
         return
      end if
      out%standard_deviation = sqrt(max(0.0_dp, eigenvalues))
      out%factors = eigenvectors(:, :components)
      out%loadings = matmul(transpose(centered), out%factors)
      out%components = components
   end function mts_apca

   pure function mts_diffusion_forecast(response, predictors, origin, factor_count) result(out)
      !! Produce Stock-Watson diffusion-index out-of-sample forecasts.
      real(dp), intent(in) :: response(:) !! Response observations.
      real(dp), intent(in) :: predictors(:, :) !! Predictor matrix.
      integer, intent(in) :: origin !! Origin.
      integer, intent(in) :: factor_count !! Number of factor.
      type(mts_diffusion_forecast_t) :: out
      real(dp), allocatable :: standardized(:, :), covariance(:, :), eigenvalues(:), eigenvectors(:, :)
      real(dp), allocatable :: design(:, :), xtx(:, :), inverse(:, :), fitted(:)
      integer :: n, variables, factors, status

      n = size(predictors, 1)
      variables = size(predictors, 2)
      factors = min(factor_count, variables)
      if (size(response) /= n .or. variables < 1 .or. origin < 2 .or. origin > n .or. &
         factor_count < 1 .or. origin <= factors + 1) then
         out%info = 1
         return
      end if
      out%predictor_mean = sum(predictors(:origin, :), 1)/real(origin, dp)
      out%predictor_scale = sqrt(sum((predictors(:origin, :) - &
         spread(out%predictor_mean, 1, origin))**2, 1)/real(origin - 1, dp))
      if (any(out%predictor_scale <= tiny(1.0_dp))) then
         out%info = 2
         return
      end if
      standardized = (predictors - spread(out%predictor_mean, 1, n))/ &
         spread(out%predictor_scale, 1, n)
      covariance = matmul(transpose(standardized(:origin, :)), standardized(:origin, :))/ &
         real(origin - 1, dp)
      call symmetric_eigen(covariance, eigenvalues, eigenvectors, status)
      if (status /= 0) then
         out%info = 10 + status
         return
      end if
      out%loadings = eigenvectors(:, :factors)
      out%index = matmul(standardized, out%loadings)
      allocate(design(origin, factors + 1))
      design(:, 1) = 1.0_dp
      design(:, 2:) = out%index(:origin, :)
      xtx = matmul(transpose(design), design)
      call invert_matrix(xtx, inverse, status)
      if (status /= 0) then
         out%info = 20 + status
         return
      end if
      out%coefficients = matmul(inverse, matmul(transpose(design), response(:origin)))
      allocate(out%forecast(n - origin))
      if (origin < n) then
         fitted = out%coefficients(1) + matmul(out%index(origin + 1:, :), out%coefficients(2:))
         out%forecast = fitted
         out%mse = sum((response(origin + 1:) - fitted)**2)/real(n - origin, dp)
      end if
      out%origin = origin
      out%factors = factors
   end function mts_diffusion_forecast

   pure function mts_multivariate_regression(response, predictors, include_constant) result(out)
      !! Fit a multivariate multiple linear regression by least squares.
      real(dp), intent(in) :: response(:, :) !! Response observations.
      real(dp), intent(in) :: predictors(:, :) !! Predictor matrix.
      logical, intent(in), optional :: include_constant !! Whether to include the constant.
      type(mts_multivariate_regression_t) :: out
      real(dp), allocatable :: design(:, :), xtx(:, :), inverse(:, :)
      logical :: use_constant
      integer :: n, responses, regressors, status, first, second

      n = size(response, 1)
      responses = size(response, 2)
      use_constant = .true.
      if (present(include_constant)) use_constant = include_constant
      regressors = size(predictors, 2) + merge(1, 0, use_constant)
      if (size(predictors, 1) /= n .or. responses < 1 .or. n <= regressors) then
         out%info = 1
         return
      end if
      allocate(design(n, regressors))
      if (use_constant) then
         design(:, 1) = 1.0_dp
         design(:, 2:) = predictors
      else
         design = predictors
      end if
      xtx = matmul(transpose(design), design)
      call invert_matrix(xtx, inverse, status)
      if (status /= 0) then
         out%info = 2
         return
      end if
      out%coefficients = matmul(inverse, matmul(transpose(design), response))
      out%residuals = response - matmul(design, out%coefficients)
      out%residual_degrees_of_freedom = n - regressors
      out%sigma = matmul(transpose(out%residuals), out%residuals)/ &
         real(out%residual_degrees_of_freedom, dp)
      allocate(out%covariance(regressors*responses, regressors*responses))
      do second = 1, responses
         do first = 1, responses
            out%covariance((first - 1)*regressors + 1:first*regressors, &
               (second - 1)*regressors + 1:second*regressors) = out%sigma(first, second)*inverse
         end do
      end do
      out%standard_errors = reshape(sqrt(max(0.0_dp, diagonal_values(out%covariance))), &
         [regressors, responses])
      out%includes_constant = use_constant
   end function mts_multivariate_regression

   pure function mts_var_chi(series, order, include_mean, threshold) result(out)
      !! Test the joint null that individually weak VAR coefficients are zero.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      integer, intent(in) :: order !! Model or polynomial order.
      logical, intent(in), optional :: include_mean !! Whether to include a mean term.
      real(dp), intent(in), optional :: threshold !! Decision or truncation threshold.
      type(mts_var_chi_t) :: out
      real(dp), allocatable :: design(:, :), response(:, :), xtx(:, :), inverse(:, :)
      real(dp), allocatable :: full_covariance(:, :), flattened(:), standard_errors(:), target_inverse(:, :)
      real(dp), allocatable :: regression_residual(:, :), regression_sigma(:, :)
      logical :: use_mean
      integer :: n, dimension, rows, regressors, lag, first, second, i, status

      n = size(series, 1)
      dimension = size(series, 2)
      use_mean = .true.
      if (present(include_mean)) use_mean = include_mean
      if (present(threshold)) out%threshold = threshold
      regressors = order*dimension + merge(1, 0, use_mean)
      rows = n - order
      if (order < 1 .or. dimension < 1 .or. rows <= regressors .or. out%threshold < 0.0_dp) then
         out%info = 1
         return
      end if
      out%fit = mts_var(series, order, use_mean)
      if (out%fit%info /= 0) then
         out%info = 10 + out%fit%info
         return
      end if
      allocate(design(rows, regressors), response(rows, dimension))
      response = series(order + 1:, :)
      i = 0
      if (use_mean) then
         i = 1
         design(:, 1) = 1.0_dp
      end if
      do lag = 1, order
         design(:, i + (lag - 1)*dimension + 1:i + lag*dimension) = &
            series(order + 1 - lag:n - lag, :)
      end do
      xtx = matmul(transpose(design), design)
      call invert_matrix(xtx, inverse, status)
      if (status /= 0) then
         out%info = 20 + status
         return
      end if
      regression_residual = response - matmul(design, out%fit%coefficients)
      regression_sigma = matmul(transpose(regression_residual), regression_residual)/ &
         real(rows - regressors, dp)
      allocate(full_covariance(regressors*dimension, regressors*dimension))
      do second = 1, dimension
         do first = 1, dimension
            full_covariance((first - 1)*regressors + 1:first*regressors, &
               (second - 1)*regressors + 1:second*regressors) = &
               regression_sigma(first, second)*inverse
         end do
      end do
      flattened = reshape(out%fit%coefficients, [regressors*dimension])
      standard_errors = sqrt(max(0.0_dp, diagonal_values(full_covariance)))
      out%coefficient_index = pack([(i, i=1, size(flattened))], &
         abs(flattened) < out%threshold*standard_errors)
      out%degrees_of_freedom = size(out%coefficient_index)
      allocate(out%values(out%degrees_of_freedom), &
         out%covariance(out%degrees_of_freedom, out%degrees_of_freedom))
      if (out%degrees_of_freedom == 0) return
      out%values = flattened(out%coefficient_index)
      out%covariance = full_covariance(out%coefficient_index, out%coefficient_index)
      call invert_matrix(out%covariance, target_inverse, status)
      if (status /= 0) then
         out%info = 30 + status
         return
      end if
      out%statistic = dot_product(out%values, matmul(target_inverse, out%values))
      out%p_value = chi_square_survival(out%statistic, out%degrees_of_freedom)
   end function mts_var_chi

   pure function mts_var_fore(fit, history, horizon) result(out)
      !! Provide the MTS VARfore interface through the shared VAR forecast engine.
      type(mts_var_fit_t), intent(in) :: fit !! Previously fitted model.
      real(dp), intent(in) :: history(:, :) !! History.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      type(mts_var_forecast_t) :: out

      out = mts_var_forecast(fit, history, horizon)
   end function mts_var_fore

   pure function mts_transfer_fit(response, input, noise_order, transfer_order, &
      initial, estimated, max_iterations, tolerance) result(out)
      !! Fit a one-input rational transfer function with ARMA disturbances.
      real(dp), intent(in) :: response(:) !! Response observations.
      real(dp), intent(in) :: input(:) !! Input.
      integer, intent(in) :: noise_order(3) !! Noise order.
      integer, intent(in) :: transfer_order(3) !! Transfer order.
      real(dp), intent(in), optional :: initial(:) !! Initial value.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      logical, intent(in), optional :: estimated(:) !! Flag controlling estimated.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      type(mts_transfer_fit_t) :: out
      type(optimization_result_t) :: optimization
      real(dp), allocatable :: y(:), x(:), parameters(:), free_values(:), hessian(:, :), inverse(:, :)
      integer, allocatable :: free_index(:)
      logical, allocatable :: mask(:)
      integer :: n, p, d, q, r, s, delay, count_parameters, free_count, i, status, limit
      real(dp) :: gradient_tolerance

      p = noise_order(1)
      d = noise_order(2)
      q = noise_order(3)
      r = transfer_order(1)
      s = transfer_order(2)
      delay = transfer_order(3)
      if (size(response) /= size(input) .or. min(p, d, q, r, s, delay) < 0 .or. &
         size(response) <= d + delay + s + max(p, q) + 3) then
         out%info = 1
         return
      end if
      y = response
      x = input
      do i = 1, d
         y = y(2:) - y(:size(y) - 1)
         x = x(2:) - x(:size(x) - 1)
      end do
      n = size(y)
      count_parameters = 1 + s + 1 + r + p + q
      allocate(parameters(count_parameters), mask(count_parameters))
      parameters = 0.0_dp
      parameters(1) = sum(y)/real(n, dp)
      if (present(initial)) then
         if (size(initial) /= count_parameters) then
            out%info = 2
            return
         end if
         parameters = initial
      end if
      mask = .true.
      if (present(estimated)) then
         if (size(estimated) /= count_parameters) then
            out%info = 2
            return
         end if
         mask = estimated
      end if
      free_count = count(mask)
      allocate(free_index(free_count), free_values(free_count))
      free_index = pack([(i, i=1, count_parameters)], mask)
      free_values = pack(parameters, mask)
      limit = 200
      gradient_tolerance = 1.0e-6_dp
      if (present(max_iterations)) limit = max_iterations
      if (present(tolerance)) gradient_tolerance = tolerance
      if (free_count > 0) then
         optimization = bfgs_minimize_fd(objective, free_values, limit, gradient_tolerance)
         parameters(free_index) = optimization%parameters
         out%iterations = optimization%iterations
         out%converged = optimization%converged
         if (optimization%info /= 0) out%info = 100 + optimization%info
      else
         out%converged = .true.
      end if
      call evaluate(parameters, out%filtered_input, out%noise, out%residuals, &
         out%sigma2, out%log_likelihood, status)
      if (status /= 0) then
         out%info = 10 + status
         return
      end if
      out%coefficients = parameters
      out%intercept = parameters(1)
      out%numerator = parameters(2:s + 2)
      out%denominator = parameters(s + 3:s + 2 + r)
      out%ar = parameters(s + 3 + r:s + 2 + r + p)
      out%ma = parameters(s + 3 + r + p:count_parameters)
      out%delay = delay
      out%differences = d
      out%aic = -2.0_dp*out%log_likelihood + 2.0_dp*real(free_count + 1, dp)
      out%bic = -2.0_dp*out%log_likelihood + log(real(size(out%residuals), dp))* &
         real(free_count + 1, dp)
      allocate(out%covariance(count_parameters, count_parameters), &
         out%standard_errors(count_parameters))
      out%covariance = 0.0_dp
      out%standard_errors = 0.0_dp
      if (free_count > 0) then
         free_values = pack(parameters, mask)
         hessian = finite_difference_hessian(objective, free_values)
         call invert_matrix(hessian, inverse, status)
         if (status == 0) then
            do i = 1, free_count
               out%covariance(free_index, free_index(i)) = inverse(:, i)
               out%standard_errors(free_index(i)) = sqrt(max(0.0_dp, inverse(i, i)))
            end do
         else if (out%info == 0) then
            out%info = 200 + status
         end if
      end if

   contains

      pure function objective(free_parameters) result(value)
         !! Return the concentrated negative Gaussian log likelihood.
         real(dp), intent(in) :: free_parameters(:) !! Free parameters.
         real(dp) :: value, variance, likelihood
         real(dp) :: full_parameters(count_parameters)
         real(dp), allocatable :: filtered(:), disturbance(:), innovations(:)
         integer :: evaluation_status

         full_parameters = parameters
         full_parameters(free_index) = free_parameters
         call evaluate(full_parameters, filtered, disturbance, innovations, variance, &
            likelihood, evaluation_status)
         if (evaluation_status == 0 .and. ieee_is_finite(likelihood)) then
            value = -likelihood
         else
            value = 1.0e30_dp + dot_product(free_parameters, free_parameters)
         end if
      end function objective

      pure subroutine evaluate(values, filtered, disturbance, innovations, variance, &
         likelihood, evaluation_status)
         !! Filter the input and recursively compute ARMA innovations.
         real(dp), intent(in) :: values(:) !! Input values.
         real(dp), allocatable, intent(out) :: filtered(:) !! Filtered.
         real(dp), allocatable, intent(out) :: disturbance(:) !! Disturbance.
         real(dp), allocatable, intent(out) :: innovations(:) !! Model innovations.
         real(dp), intent(out) :: variance !! Variance value or matrix.
         real(dp), intent(out) :: likelihood !! Likelihood.
         integer, intent(out) :: evaluation_status !! Evaluation status.
         integer :: t, lag, start, effective, offset

         allocate(filtered(n), disturbance(n), innovations(n))
         filtered = x
         do t = 1, n
            do lag = 1, min(r, t - 1)
               filtered(t) = filtered(t) + values(s + 2 + lag)*filtered(t - lag)
            end do
         end do
         disturbance = y - values(1)
         do t = 1, n
            do lag = 0, s
               if (t - delay - lag >= 1) disturbance(t) = disturbance(t) - &
                  values(2 + lag)*filtered(t - delay - lag)
            end do
         end do
         innovations = 0.0_dp
         offset = s + 2 + r
         do t = 1, n
            innovations(t) = disturbance(t)
            do lag = 1, min(p, t - 1)
               innovations(t) = innovations(t) - values(offset + lag)*disturbance(t - lag)
            end do
            do lag = 1, min(q, t - 1)
               innovations(t) = innovations(t) + values(offset + p + lag)*innovations(t - lag)
            end do
         end do
         start = max(delay + s + 1, max(p, q) + 1)
         effective = n - start + 1
         variance = sum(innovations(start:)**2)/real(effective, dp)
         if (variance <= tiny(1.0_dp) .or. .not. ieee_is_finite(variance)) then
            evaluation_status = 1
            likelihood = -huge(1.0_dp)
            return
         end if
         innovations = innovations(start:)
         likelihood = -0.5_dp*real(effective, dp)*(log(2.0_dp*acos(-1.0_dp)*variance) + 1.0_dp)
         evaluation_status = 0
      end subroutine evaluate
   end function mts_transfer_fit

   pure function mts_transfer2_fit(response, input1, noise_order, seasonal_order, season, &
      transfer_order1, input2, transfer_order2, deterministic, equilibrium, initial, &
      estimated, max_iterations, tolerance) result(out)
      !! Fit a two-input transfer model with multiplicative seasonal ARMA noise.
      real(dp), intent(in) :: response(:) !! Response observations.
      real(dp), intent(in) :: input1(:) !! Input1.
      integer, intent(in) :: noise_order(3) !! Noise order.
      integer, intent(in) :: seasonal_order(3) !! Seasonal order.
      integer, intent(in) :: season !! Season.
      integer, intent(in) :: transfer_order1(3) !! Transfer order1.
      real(dp), intent(in), optional :: input2(:) !! Input2.
      real(dp), intent(in), optional :: deterministic(:) !! Deterministic.
      real(dp), intent(in), optional :: equilibrium(:) !! Equilibrium.
      integer, intent(in), optional :: transfer_order2(3) !! Transfer order2.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      real(dp), intent(in), optional :: initial(:) !! Initial value.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      logical, intent(in), optional :: estimated(:) !! Flag controlling estimated.
      type(mts_transfer2_fit_t) :: out
      type(optimization_result_t) :: optimization
      real(dp), allocatable :: data(:, :), parameters(:), free_values(:), hessian(:, :), inverse(:, :)
      integer, allocatable :: free_index(:)
      logical, allocatable :: mask(:)
      logical :: second, trend, equilibrating
      integer :: n, p, d, q, sp, sd, sq, r1, s1, b1, r2, s2, b2
      integer :: parameter_count, free_count, offset, i, status, limit
      real(dp) :: gradient_tolerance

      n = size(response)
      second = present(input2)
      trend = present(deterministic)
      equilibrating = present(equilibrium)
      if (size(input1) /= n .or. (second .neqv. present(transfer_order2))) then
         out%info = 1
         return
      end if
      if (second) then
         if (size(input2) /= n) then
            out%info = 1
            return
         end if
      end if
      if (trend) then
         if (size(deterministic) /= n) then
            out%info = 1
            return
         end if
      end if
      if (equilibrating) then
         if (size(equilibrium) /= n) then
            out%info = 1
            return
         end if
      end if
      p = noise_order(1)
      d = noise_order(2)
      q = noise_order(3)
      sp = seasonal_order(1)
      sd = seasonal_order(2)
      sq = seasonal_order(3)
      r1 = transfer_order1(1)
      s1 = transfer_order1(2)
      b1 = transfer_order1(3)
      r2 = 0
      s2 = 0
      b2 = 0
      if (second) then
         r2 = transfer_order2(1)
         s2 = transfer_order2(2)
         b2 = transfer_order2(3)
      end if
      if (min(p, d, q, sp, sd, sq, r1, s1, b1, r2, s2, b2) < 0 .or. &
         ((sp + sd + sq > 0) .and. season < 2)) then
         out%info = 1
         return
      end if
      allocate(data(n, 5))
      data = 0.0_dp
      data(:, 1) = response
      data(:, 2) = input1
      if (second) data(:, 3) = input2
      if (trend) data(:, 4) = deterministic
      if (equilibrating) data(:, 5) = equilibrium
      do i = 1, d
         data = data(2:, :) - data(:size(data, 1) - 1, :)
      end do
      do i = 1, sd
         if (size(data, 1) <= season) then
            out%info = 1
            return
         end if
         data = data(season + 1:, :) - data(:size(data, 1) - season, :)
      end do
      n = size(data, 1)
      parameter_count = 1 + s1 + 1 + r1 + merge(1, 0, trend) + merge(1, 0, equilibrating) + &
         merge(s2 + 1 + r2, 0, second) + p + q + sp + sq
      allocate(parameters(parameter_count), mask(parameter_count))
      parameters = 0.0_dp
      parameters(1) = sum(data(:, 1))/real(n, dp)
      if (present(initial)) then
         if (size(initial) /= parameter_count) then
            out%info = 2
            return
         end if
         parameters = initial
      end if
      mask = .true.
      if (present(estimated)) then
         if (size(estimated) /= parameter_count) then
            out%info = 2
            return
         end if
         mask = estimated
      end if
      free_count = count(mask)
      allocate(free_index(free_count), free_values(free_count))
      free_index = pack([(i, i=1, parameter_count)], mask)
      free_values = pack(parameters, mask)
      limit = 200
      gradient_tolerance = 1.0e-6_dp
      if (present(max_iterations)) limit = max_iterations
      if (present(tolerance)) gradient_tolerance = tolerance
      if (free_count > 0) then
         optimization = bfgs_minimize_fd(objective, free_values, limit, gradient_tolerance)
         parameters(free_index) = optimization%parameters
         out%iterations = optimization%iterations
         out%converged = optimization%converged
         if (optimization%info /= 0) out%info = 100 + optimization%info
      else
         out%converged = .true.
      end if
      call evaluate(parameters, out%filtered_input1, out%filtered_input2, out%noise, &
         out%residuals, out%sigma2, out%log_likelihood, status)
      if (status /= 0) then
         out%info = 10 + status
         return
      end if
      out%coefficients = parameters
      out%intercept = parameters(1)
      offset = 1
      out%numerator1 = parameters(offset + 1:offset + s1 + 1)
      offset = offset + s1 + 1
      out%denominator1 = parameters(offset + 1:offset + r1)
      offset = offset + r1
      if (trend) then
         out%deterministic_coefficient = parameters(offset + 1)
         offset = offset + 1
      end if
      if (equilibrating) then
         out%equilibrium_coefficient = parameters(offset + 1)
         offset = offset + 1
      end if
      if (second) then
         out%numerator2 = parameters(offset + 1:offset + s2 + 1)
         offset = offset + s2 + 1
         out%denominator2 = parameters(offset + 1:offset + r2)
         offset = offset + r2
      else
         allocate(out%numerator2(0), out%denominator2(0))
      end if
      out%ar = parameters(offset + 1:offset + p)
      offset = offset + p
      out%ma = parameters(offset + 1:offset + q)
      offset = offset + q
      out%seasonal_ar = parameters(offset + 1:offset + sp)
      offset = offset + sp
      out%seasonal_ma = parameters(offset + 1:offset + sq)
      out%delay1 = b1
      out%delay2 = b2
      out%differences = d
      out%seasonal_differences = sd
      out%season = season
      out%has_second_input = second
      out%has_deterministic = trend
      out%has_equilibrium = equilibrating
      out%aic = -2.0_dp*out%log_likelihood + 2.0_dp*real(free_count + 1, dp)
      out%bic = -2.0_dp*out%log_likelihood + log(real(size(out%residuals), dp))* &
         real(free_count + 1, dp)
      allocate(out%covariance(parameter_count, parameter_count), out%standard_errors(parameter_count))
      out%covariance = 0.0_dp
      out%standard_errors = 0.0_dp
      if (free_count > 0) then
         free_values = pack(parameters, mask)
         hessian = finite_difference_hessian(objective, free_values)
         call invert_matrix(hessian, inverse, status)
         if (status == 0) then
            do i = 1, free_count
               out%covariance(free_index, free_index(i)) = inverse(:, i)
               out%standard_errors(free_index(i)) = sqrt(max(0.0_dp, inverse(i, i)))
            end do
         else if (out%info == 0) then
            out%info = 200 + status
         end if
      end if

   contains

      pure function objective(free_parameters) result(value)
         !! Return the concentrated negative likelihood for the two-input model.
         real(dp), intent(in) :: free_parameters(:) !! Free parameters.
         real(dp) :: value, variance, likelihood
         real(dp) :: full_parameters(parameter_count)
         real(dp), allocatable :: filtered1(:), filtered2(:), disturbance(:), innovations(:)
         integer :: evaluation_status

         full_parameters = parameters
         full_parameters(free_index) = free_parameters
         call evaluate(full_parameters, filtered1, filtered2, disturbance, innovations, &
            variance, likelihood, evaluation_status)
         if (evaluation_status == 0 .and. ieee_is_finite(likelihood)) then
            value = -likelihood
         else
            value = 1.0e30_dp + dot_product(free_parameters, free_parameters)
         end if
      end function objective

      pure subroutine evaluate(values, filtered1, filtered2, disturbance, innovations, &
         variance, likelihood, evaluation_status)
         !! Evaluate filtered inputs and multiplicative seasonal ARMA innovations.
         real(dp), intent(in) :: values(:) !! Input values.
         real(dp), allocatable, intent(out) :: filtered1(:) !! Filtered1.
         real(dp), allocatable, intent(out) :: filtered2(:) !! Filtered2.
         real(dp), allocatable, intent(out) :: disturbance(:) !! Disturbance.
         real(dp), allocatable, intent(out) :: innovations(:) !! Model innovations.
         real(dp), intent(out) :: variance !! Variance value or matrix.
         real(dp), intent(out) :: likelihood !! Likelihood.
         integer, intent(out) :: evaluation_status !! Evaluation status.
         real(dp), allocatable :: regular(:)
         integer :: t, lag, position, ar_position, ma_position, sar_position, sma_position
         integer :: start, effective

         allocate(filtered1(n), filtered2(n), disturbance(n), innovations(n), regular(n))
         filtered1 = data(:, 2)
         filtered2 = data(:, 3)
         position = 1 + s1 + 1
         do t = 1, n
            do lag = 1, min(r1, t - 1)
               filtered1(t) = filtered1(t) + values(position + lag)*filtered1(t - lag)
            end do
         end do
         position = position + r1 + merge(1, 0, trend) + merge(1, 0, equilibrating)
         if (second) then
            position = position + s2 + 1
            do t = 1, n
               do lag = 1, min(r2, t - 1)
                  filtered2(t) = filtered2(t) + values(position + lag)*filtered2(t - lag)
               end do
            end do
            position = position + r2
         end if
         disturbance = data(:, 1) - values(1)
         do t = 1, n
            do lag = 0, s1
               if (t - b1 - lag >= 1) disturbance(t) = disturbance(t) - &
                  values(2 + lag)*filtered1(t - b1 - lag)
            end do
         end do
         position = 1 + s1 + 1 + r1
         if (trend) then
            disturbance = disturbance - values(position + 1)*data(:, 4)
            position = position + 1
         end if
         if (equilibrating) then
            disturbance = disturbance - values(position + 1)*data(:, 5)
            position = position + 1
         end if
         if (second) then
            do t = 1, n
               do lag = 0, s2
                  if (t - b2 - lag >= 1) disturbance(t) = disturbance(t) - &
                     values(position + 1 + lag)*filtered2(t - b2 - lag)
               end do
            end do
            position = position + s2 + 1 + r2
         end if
         ar_position = position
         ma_position = ar_position + p
         sar_position = ma_position + q
         sma_position = sar_position + sp
         regular = 0.0_dp
         do t = 1, n
            regular(t) = disturbance(t)
            do lag = 1, min(p, t - 1)
               regular(t) = regular(t) - values(ar_position + lag)*disturbance(t - lag)
            end do
            do lag = 1, min(q, t - 1)
               regular(t) = regular(t) + values(ma_position + lag)*regular(t - lag)
            end do
         end do
         innovations = 0.0_dp
         do t = 1, n
            innovations(t) = regular(t)
            do lag = 1, sp
               if (t - lag*season >= 1) innovations(t) = innovations(t) - &
                  values(sar_position + lag)*regular(t - lag*season)
            end do
            do lag = 1, sq
               if (t - lag*season >= 1) innovations(t) = innovations(t) + &
                  values(sma_position + lag)*innovations(t - lag*season)
            end do
         end do
         start = max(max(b1 + s1, b2 + s2) + 1, &
            max(max(p, q), max(sp, sq)*season) + 1)
         effective = n - start + 1
         variance = sum(innovations(start:)**2)/real(effective, dp)
         if (effective < 2 .or. variance <= tiny(1.0_dp) .or. .not. ieee_is_finite(variance)) then
            evaluation_status = 1
            likelihood = -huge(1.0_dp)
            return
         end if
         innovations = innovations(start:)
         likelihood = -0.5_dp*real(effective, dp)*(log(2.0_dp*acos(-1.0_dp)*variance) + 1.0_dp)
         evaluation_status = 0
      end subroutine evaluate
   end function mts_transfer2_fit

   pure function mts_transfer2_forecast(fit, response, input1, future_input1, &
      noise_order, seasonal_order, transfer_order1, input2, future_input2, &
      transfer_order2, deterministic, future_deterministic, equilibrium, &
      future_equilibrium) result(out)
      !! Forecast a fitted two-input transfer model with zero future innovations.
      type(mts_transfer2_fit_t), intent(in) :: fit !! Previously fitted model.
      real(dp), intent(in) :: response(:) !! Response observations.
      real(dp), intent(in) :: input1(:) !! Input1.
      real(dp), intent(in) :: future_input1(:) !! Future input1.
      integer, intent(in) :: noise_order(3) !! Noise order.
      integer, intent(in) :: seasonal_order(3) !! Seasonal order.
      integer, intent(in) :: transfer_order1(3) !! Transfer order1.
      real(dp), intent(in), optional :: input2(:) !! Input2.
      real(dp), intent(in), optional :: future_input2(:) !! Future input2.
      integer, intent(in), optional :: transfer_order2(3) !! Transfer order2.
      real(dp), intent(in), optional :: deterministic(:) !! Deterministic.
      real(dp), intent(in), optional :: future_deterministic(:) !! Future deterministic.
      real(dp), intent(in), optional :: equilibrium(:) !! Equilibrium.
      real(dp), intent(in), optional :: future_equilibrium(:) !! Future equilibrium.
      type(mts_transfer_forecast_t) :: out
      type(mts_transfer2_fit_t) :: candidate
      real(dp), allocatable :: y(:), x1(:), x2(:), det(:), eq(:)
      logical, allocatable :: fixed(:)
      integer :: n, horizon, step

      n = size(response)
      horizon = size(future_input1)
      if (fit%info /= 0 .or. size(input1) /= n .or. horizon < 1 .or. &
         (fit%has_second_input .neqv. present(input2)) .or. &
         (fit%has_second_input .neqv. present(future_input2)) .or. &
         (fit%has_deterministic .neqv. present(deterministic)) .or. &
         (fit%has_deterministic .neqv. present(future_deterministic)) .or. &
         (fit%has_equilibrium .neqv. present(equilibrium)) .or. &
         (fit%has_equilibrium .neqv. present(future_equilibrium))) then
         out%info = 1
         return
      end if
      if (fit%has_second_input) then
         if (.not. present(transfer_order2) .or. size(input2) /= n .or. &
            size(future_input2) /= horizon) then
            out%info = 1
            return
         end if
      end if
      allocate(y(n + horizon), x1(n + horizon), out%mean(horizon))
      y(:n) = response
      x1(:n) = input1
      x1(n + 1:) = future_input1
      if (fit%has_second_input) then
         allocate(x2(n + horizon))
         x2(:n) = input2
         x2(n + 1:) = future_input2
      end if
      if (fit%has_deterministic) then
         allocate(det(n + horizon))
         det(:n) = deterministic
         det(n + 1:) = future_deterministic
      end if
      if (fit%has_equilibrium) then
         allocate(eq(n + horizon))
         eq(:n) = equilibrium
         eq(n + 1:) = future_equilibrium
      end if
      allocate(fixed(size(fit%coefficients)))
      fixed = .false.
      do step = 1, horizon
         y(n + step) = 0.0_dp
         if (fit%has_second_input .and. fit%has_deterministic .and. fit%has_equilibrium) then
            candidate = mts_transfer2_fit(y(:n + step), x1(:n + step), noise_order, &
               seasonal_order, fit%season, transfer_order1, x2(:n + step), transfer_order2, &
               det(:n + step), eq(:n + step), fit%coefficients, fixed)
         else if (fit%has_second_input .and. fit%has_deterministic) then
            candidate = mts_transfer2_fit(y(:n + step), x1(:n + step), noise_order, &
               seasonal_order, fit%season, transfer_order1, x2(:n + step), transfer_order2, &
               deterministic=det(:n + step), initial=fit%coefficients, estimated=fixed)
         else if (fit%has_second_input .and. fit%has_equilibrium) then
            candidate = mts_transfer2_fit(y(:n + step), x1(:n + step), noise_order, &
               seasonal_order, fit%season, transfer_order1, x2(:n + step), transfer_order2, &
               equilibrium=eq(:n + step), initial=fit%coefficients, estimated=fixed)
         else if (fit%has_second_input) then
            candidate = mts_transfer2_fit(y(:n + step), x1(:n + step), noise_order, &
               seasonal_order, fit%season, transfer_order1, x2(:n + step), transfer_order2, &
               initial=fit%coefficients, estimated=fixed)
         else if (fit%has_deterministic .and. fit%has_equilibrium) then
            candidate = mts_transfer2_fit(y(:n + step), x1(:n + step), noise_order, &
               seasonal_order, fit%season, transfer_order1, deterministic=det(:n + step), &
               equilibrium=eq(:n + step), initial=fit%coefficients, estimated=fixed)
         else if (fit%has_deterministic) then
            candidate = mts_transfer2_fit(y(:n + step), x1(:n + step), noise_order, &
               seasonal_order, fit%season, transfer_order1, deterministic=det(:n + step), &
               initial=fit%coefficients, estimated=fixed)
         else if (fit%has_equilibrium) then
            candidate = mts_transfer2_fit(y(:n + step), x1(:n + step), noise_order, &
               seasonal_order, fit%season, transfer_order1, equilibrium=eq(:n + step), &
               initial=fit%coefficients, estimated=fixed)
         else
            candidate = mts_transfer2_fit(y(:n + step), x1(:n + step), noise_order, &
               seasonal_order, fit%season, transfer_order1, initial=fit%coefficients, &
               estimated=fixed)
         end if
         if (candidate%info /= 0) then
            out%info = 10 + step
            return
         end if
         y(n + step) = -candidate%residuals(size(candidate%residuals))
         out%mean(step) = y(n + step)
      end do
   end function mts_transfer2_forecast

   pure function mts_transfer2_backtest(response, input1, origin, noise_order, &
      seasonal_order, season, transfer_order1, initial, estimated, input2, &
      transfer_order2, deterministic, equilibrium, max_iterations, tolerance) result(out)
      !! Run rolling one-step re-estimation and forecasts for a transfer model.
      real(dp), intent(in) :: response(:) !! Response observations.
      real(dp), intent(in) :: input1(:) !! Input1.
      real(dp), intent(in) :: initial(:) !! Initial value.
      integer, intent(in) :: origin !! Origin.
      integer, intent(in) :: noise_order(3) !! Noise order.
      integer, intent(in) :: seasonal_order(3) !! Seasonal order.
      integer, intent(in) :: season !! Season.
      integer, intent(in) :: transfer_order1(3) !! Transfer order1.
      logical, intent(in) :: estimated(:) !! Flag controlling estimated.
      real(dp), intent(in), optional :: input2(:) !! Input2.
      real(dp), intent(in), optional :: deterministic(:) !! Deterministic.
      real(dp), intent(in), optional :: equilibrium(:) !! Equilibrium.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      integer, intent(in), optional :: transfer_order2(3) !! Transfer order2.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      type(mts_transfer_backtest_t) :: out
      type(mts_transfer2_fit_t) :: fitted, candidate
      real(dp), allocatable :: y(:), start(:)
      logical, allocatable :: fixed(:)
      integer :: n, current, index

      n = size(response)
      if (size(input1) /= n .or. origin < 2 .or. origin >= n .or. &
         size(initial) /= size(estimated)) then
         out%info = 1
         return
      end if
      allocate(out%error(n - origin), y(n), start(size(initial)), fixed(size(initial)))
      y = response
      start = initial
      fixed = .false.
      do current = origin, n - 1
         fitted = fit_sample(response(:current), current, start, estimated)
         if (fitted%info /= 0 .and. fitted%info < 100) then
            out%info = 10 + current - origin
            return
         end if
         y(:current) = response(:current)
         y(current + 1) = 0.0_dp
         candidate = fit_sample(y(:current + 1), current + 1, fitted%coefficients, fixed)
         if (candidate%info /= 0) then
            out%info = 20 + current - origin
            return
         end if
         index = current - origin + 1
         out%error(index) = response(current + 1) + &
            candidate%residuals(size(candidate%residuals))
         start = fitted%coefficients
      end do
      out%bias = sum(out%error)/real(size(out%error), dp)
      out%mse = sum(out%error**2)/real(size(out%error), dp)
      out%rmse = sqrt(out%mse)
      out%mean_absolute_error = sum(abs(out%error))/real(size(out%error), dp)
      out%origin = origin

   contains

      pure function fit_sample(y_values, observations, starting, mask_values) result(model)
         !! Fit one rolling sample while preserving the optional regressor layout.
         real(dp), intent(in) :: y_values(:) !! Y values.
         real(dp), intent(in) :: starting(:) !! Starting.
         integer, intent(in) :: observations !! Observed time-series values.
         logical, intent(in) :: mask_values(:) !! Flag controlling mask values.
         type(mts_transfer2_fit_t) :: model

         if (present(input2) .and. present(deterministic) .and. present(equilibrium)) then
            model = mts_transfer2_fit(y_values, input1(:observations), noise_order, &
               seasonal_order, season, transfer_order1, input2(:observations), &
               transfer_order2, deterministic(:observations), equilibrium(:observations), &
               starting, mask_values, max_iterations, tolerance)
         else if (present(input2) .and. present(deterministic)) then
            model = mts_transfer2_fit(y_values, input1(:observations), noise_order, &
               seasonal_order, season, transfer_order1, input2(:observations), &
               transfer_order2, deterministic=deterministic(:observations), initial=starting, &
               estimated=mask_values, max_iterations=max_iterations, tolerance=tolerance)
         else if (present(input2) .and. present(equilibrium)) then
            model = mts_transfer2_fit(y_values, input1(:observations), noise_order, &
               seasonal_order, season, transfer_order1, input2(:observations), &
               transfer_order2, equilibrium=equilibrium(:observations), initial=starting, &
               estimated=mask_values, max_iterations=max_iterations, tolerance=tolerance)
         else if (present(input2)) then
            model = mts_transfer2_fit(y_values, input1(:observations), noise_order, &
               seasonal_order, season, transfer_order1, input2(:observations), &
               transfer_order2, initial=starting, estimated=mask_values, &
               max_iterations=max_iterations, tolerance=tolerance)
         else if (present(deterministic) .and. present(equilibrium)) then
            model = mts_transfer2_fit(y_values, input1(:observations), noise_order, &
               seasonal_order, season, transfer_order1, deterministic=deterministic(:observations), &
               equilibrium=equilibrium(:observations), initial=starting, estimated=mask_values, &
               max_iterations=max_iterations, tolerance=tolerance)
         else if (present(deterministic)) then
            model = mts_transfer2_fit(y_values, input1(:observations), noise_order, &
               seasonal_order, season, transfer_order1, deterministic=deterministic(:observations), &
               initial=starting, estimated=mask_values, max_iterations=max_iterations, &
               tolerance=tolerance)
         else if (present(equilibrium)) then
            model = mts_transfer2_fit(y_values, input1(:observations), noise_order, &
               seasonal_order, season, transfer_order1, equilibrium=equilibrium(:observations), &
               initial=starting, estimated=mask_values, max_iterations=max_iterations, &
               tolerance=tolerance)
         else
            model = mts_transfer2_fit(y_values, input1(:observations), noise_order, &
               seasonal_order, season, transfer_order1, initial=starting, estimated=mask_values, &
               max_iterations=max_iterations, tolerance=tolerance)
         end if
      end function fit_sample
   end function mts_transfer2_backtest

   pure function mts_regts_fit(response, regressors, order, include_mean, initial, &
      estimated, max_iterations, tolerance) result(out)
      !! Fit multivariate regression with jointly estimated VAR disturbances.
      real(dp), intent(in) :: response(:, :) !! Response observations.
      real(dp), intent(in) :: regressors(:, :) !! Regression design matrix.
      integer, intent(in) :: order !! Model or polynomial order.
      logical, intent(in), optional :: include_mean !! Whether to include a mean term.
      logical, intent(in), optional :: estimated(:, :) !! Flag controlling estimated.
      real(dp), intent(in), optional :: initial(:) !! Initial value.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      type(mts_regts_fit_t) :: out
      type(optimization_result_t) :: optimization
      type(mts_var_fit_t) :: var_start
      real(dp), allocatable :: design(:, :), beta_start(:, :), regression_residual(:, :)
      real(dp), allocatable :: xtx(:, :), xtx_inverse(:, :), parameters(:), free_values(:)
      real(dp), allocatable :: hessian(:, :), inverse(:, :)
      integer, allocatable :: free_index(:)
      logical, allocatable :: flat_mask(:)
      logical :: use_mean
      integer :: n, dimension, regressor_count, rows_per_equation, parameter_count
      integer :: free_count, equation, lag, offset, i, status, limit
      real(dp) :: gradient_tolerance

      n = size(response, 1)
      dimension = size(response, 2)
      use_mean = .true.
      if (present(include_mean)) use_mean = include_mean
      regressor_count = size(regressors, 2) + merge(1, 0, use_mean)
      rows_per_equation = regressor_count + order*dimension
      parameter_count = rows_per_equation*dimension
      if (size(regressors, 1) /= n .or. dimension < 1 .or. order < 0 .or. &
         n <= order + rows_per_equation + 1) then
         out%info = 1
         return
      end if
      allocate(design(n, regressor_count))
      if (use_mean) then
         design(:, 1) = 1.0_dp
         design(:, 2:) = regressors
      else
         design = regressors
      end if
      xtx = matmul(transpose(design), design)
      call invert_matrix(xtx, xtx_inverse, status)
      if (status /= 0) then
         out%info = 2
         return
      end if
      beta_start = matmul(xtx_inverse, matmul(transpose(design), response))
      regression_residual = response - matmul(design, beta_start)
      allocate(parameters(parameter_count), out%estimated(rows_per_equation, dimension))
      parameters = 0.0_dp
      do equation = 1, dimension
         offset = (equation - 1)*rows_per_equation
         parameters(offset + 1:offset + regressor_count) = beta_start(:, equation)
      end do
      if (order > 0) then
         var_start = mts_var(regression_residual, order, .false.)
         if (var_start%info == 0) then
            do equation = 1, dimension
               offset = (equation - 1)*rows_per_equation + regressor_count
               do lag = 1, order
                  parameters(offset + (lag - 1)*dimension + 1:offset + lag*dimension) = &
                     var_start%ar(equation, :, lag)
               end do
            end do
         end if
      end if
      if (present(initial)) then
         if (size(initial) /= parameter_count) then
            out%info = 3
            return
         end if
         parameters = initial
      end if
      out%estimated = .true.
      if (present(estimated)) then
         if (any(shape(estimated) /= [rows_per_equation, dimension])) then
            out%info = 3
            return
         end if
         out%estimated = estimated
      end if
      flat_mask = reshape(out%estimated, [parameter_count])
      free_count = count(flat_mask)
      allocate(free_index(free_count), free_values(free_count))
      free_index = pack([(i, i=1, parameter_count)], flat_mask)
      free_values = pack(parameters, flat_mask)
      limit = 200
      gradient_tolerance = 1.0e-6_dp
      if (present(max_iterations)) limit = max_iterations
      if (present(tolerance)) gradient_tolerance = tolerance
      if (free_count > 0) then
         optimization = bfgs_minimize_fd(objective, free_values, limit, gradient_tolerance)
         parameters(free_index) = optimization%parameters
         out%iterations = optimization%iterations
         out%converged = optimization%converged
         if (optimization%info /= 0) out%info = 100 + optimization%info
      else
         out%converged = .true.
      end if
      call evaluate(parameters, out%beta, out%ar, out%residuals, out%sigma, &
         out%log_likelihood, status)
      if (status /= 0) then
         out%info = 10 + status
         return
      end if
      out%coefficients = parameters
      out%order = order
      out%includes_mean = use_mean
      out%aic = -2.0_dp*out%log_likelihood + 2.0_dp*real(free_count + 1, dp)
      out%bic = -2.0_dp*out%log_likelihood + log(real(n - order, dp))*real(free_count + 1, dp)
      allocate(out%covariance(parameter_count, parameter_count), out%standard_errors(parameter_count))
      out%covariance = 0.0_dp
      out%standard_errors = 0.0_dp
      if (free_count > 0) then
         free_values = pack(parameters, flat_mask)
         hessian = finite_difference_hessian(objective, free_values)
         call invert_matrix(hessian, inverse, status)
         if (status == 0) then
            do i = 1, free_count
               out%covariance(free_index, free_index(i)) = inverse(:, i)
               out%standard_errors(free_index(i)) = sqrt(max(0.0_dp, inverse(i, i)))
            end do
         else if (out%info == 0) then
            out%info = 200 + status
         end if
      end if

   contains

      pure function objective(free_parameters) result(value)
         !! Return the concentrated negative likelihood for free parameters.
         real(dp), intent(in) :: free_parameters(:) !! Free parameters.
         real(dp) :: value, likelihood
         real(dp) :: full_parameters(parameter_count)
         real(dp), allocatable :: beta(:, :), ar(:, :, :), residual_values(:, :), sigma(:, :)
         integer :: evaluation_status

         full_parameters = parameters
         full_parameters(free_index) = free_parameters
         call evaluate(full_parameters, beta, ar, residual_values, sigma, likelihood, evaluation_status)
         if (evaluation_status == 0 .and. ieee_is_finite(likelihood)) then
            value = -likelihood
         else
            value = 1.0e30_dp + dot_product(free_parameters, free_parameters)
         end if
      end function objective

      pure subroutine evaluate(values, beta, ar, residual_values, sigma, likelihood, &
         evaluation_status)
         !! Unpack coefficients and evaluate regression VAR innovations.
         real(dp), intent(in) :: values(:) !! Input values.
         real(dp), allocatable, intent(out) :: beta(:, :) !! Regression or model coefficients.
         real(dp), allocatable, intent(out) :: ar(:, :, :) !! Autoregressive coefficients.
         real(dp), allocatable, intent(out) :: residual_values(:, :) !! Residual values.
         real(dp), allocatable, intent(out) :: sigma(:, :) !! Scale parameter or standard deviation.
         real(dp), intent(out) :: likelihood !! Likelihood.
         integer, intent(out) :: evaluation_status !! Evaluation status.
         real(dp), allocatable :: raw_residual(:, :), covariance_inverse(:, :)
         real(dp) :: logdet
         integer :: current_equation, current_lag, current_offset, t

         allocate(beta(regressor_count, dimension), ar(dimension, dimension, order))
         do current_equation = 1, dimension
            current_offset = (current_equation - 1)*rows_per_equation
            beta(:, current_equation) = values(current_offset + 1:current_offset + regressor_count)
            current_offset = current_offset + regressor_count
            do current_lag = 1, order
               ar(current_equation, :, current_lag) = &
                  values(current_offset + (current_lag - 1)*dimension + 1: &
                  current_offset + current_lag*dimension)
            end do
         end do
         raw_residual = response - matmul(design, beta)
         allocate(residual_values(n - order, dimension))
         do t = order + 1, n
            residual_values(t - order, :) = raw_residual(t, :)
            do current_lag = 1, order
               residual_values(t - order, :) = residual_values(t - order, :) - &
                  matmul(ar(:, :, current_lag), raw_residual(t - current_lag, :))
            end do
         end do
         sigma = matmul(transpose(residual_values), residual_values)/real(n - order, dp)
         allocate(covariance_inverse(dimension, dimension))
         call inverse_logdet(sigma, covariance_inverse, logdet, evaluation_status, &
            100.0_dp*epsilon(1.0_dp))
         if (evaluation_status /= 0) then
            likelihood = -huge(1.0_dp)
            return
         end if
         likelihood = -0.5_dp*real(n - order, dp)*(real(dimension, dp)* &
            (log(2.0_dp*acos(-1.0_dp)) + 1.0_dp) + logdet)
      end subroutine evaluate
   end function mts_regts_fit

   pure function mts_regts_refine(response, regressors, fit, threshold, max_steps, &
      max_iterations, tolerance) result(out)
      !! Remove insignificant regression or VAR-error parameters and refit.
      real(dp), intent(in) :: response(:, :) !! Response observations.
      real(dp), intent(in) :: regressors(:, :) !! Regression design matrix.
      type(mts_regts_fit_t), intent(in) :: fit !! Previously fitted model.
      real(dp), intent(in), optional :: threshold !! Decision or truncation threshold.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      integer, intent(in), optional :: max_steps !! Maximum steps.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      type(mts_regts_refinement_t) :: out
      logical, allocatable :: mask(:, :), flat_mask(:)
      real(dp) :: cutoff, ratio, smallest
      integer :: limit, iteration_limit, step, weakest, i

      cutoff = 1.0_dp
      limit = size(fit%coefficients)
      iteration_limit = 200
      if (present(threshold)) cutoff = threshold
      if (present(max_steps)) limit = min(limit, max_steps)
      if (present(max_iterations)) iteration_limit = max_iterations
      if (fit%info /= 0 .or. cutoff < 0.0_dp) then
         out%info = 1
         return
      end if
      mask = fit%estimated
      out%fit = fit
      allocate(out%active_count(0:limit), out%removed_index(limit))
      out%active_count(0) = count(mask)
      do step = 1, limit
         flat_mask = reshape(mask, [size(mask)])
         weakest = 0
         smallest = huge(1.0_dp)
         do i = 1, size(flat_mask)
            if (.not. flat_mask(i) .or. out%fit%standard_errors(i) <= 0.0_dp) cycle
            ratio = abs(out%fit%coefficients(i)/out%fit%standard_errors(i))
            if (ratio < cutoff .and. ratio < smallest) then
               weakest = i
               smallest = ratio
            end if
         end do
         if (weakest == 0) exit
         flat_mask(weakest) = .false.
         mask = reshape(flat_mask, shape(mask))
         out%fit = mts_regts_fit(response, regressors, fit%order, fit%includes_mean, &
            out%fit%coefficients, mask, iteration_limit, tolerance)
         if (out%fit%info /= 0 .and. out%fit%info < 100) then
            out%info = 10 + step
            return
         end if
         out%removed_index(step) = weakest
         out%active_count(step) = count(mask)
         out%steps = step
      end do
      out%active_count = out%active_count(:out%steps)
      out%removed_index = out%removed_index(:out%steps)
   end function mts_regts_refine

   pure function mts_regts_forecast(fit, response, regressors, future_regressors) result(out)
      !! Forecast regression responses and VAR-error uncertainty.
      type(mts_regts_fit_t), intent(in) :: fit !! Previously fitted model.
      real(dp), intent(in) :: response(:, :) !! Response observations.
      real(dp), intent(in) :: regressors(:, :) !! Regression design matrix.
      real(dp), intent(in) :: future_regressors(:, :) !! Future regressors.
      type(mts_var_forecast_t) :: out
      type(mts_var_fit_t) :: error_model
      type(mts_var_forecast_t) :: error_forecast
      real(dp), allocatable :: design(:, :), future_design(:, :), errors(:, :)
      integer :: n, dimension, horizon

      n = size(response, 1)
      dimension = size(response, 2)
      horizon = size(future_regressors, 1)
      if (fit%info /= 0 .or. size(regressors, 1) /= n .or. &
         size(regressors, 2) /= size(future_regressors, 2) .or. horizon < 1) then
         out%info = 1
         return
      end if
      allocate(design(n, size(fit%beta, 1)), future_design(horizon, size(fit%beta, 1)))
      if (fit%includes_mean) then
         design(:, 1) = 1.0_dp
         design(:, 2:) = regressors
         future_design(:, 1) = 1.0_dp
         future_design(:, 2:) = future_regressors
      else
         design = regressors
         future_design = future_regressors
      end if
      errors = response - matmul(design, fit%beta)
      error_model%ar = fit%ar
      allocate(error_model%intercept(dimension), error_model%sigma(dimension, dimension))
      error_model%intercept = 0.0_dp
      error_model%sigma = fit%sigma
      error_model%max_lag = fit%order
      error_model%includes_mean = .false.
      error_forecast = mts_var_forecast(error_model, errors, horizon)
      if (error_forecast%info /= 0) then
         out%info = 10 + error_forecast%info
         return
      end if
      out = error_forecast
      out%mean = out%mean + matmul(future_design, fit%beta)
   end function mts_regts_forecast

   pure function mts_var(series, order, include_mean, estimated) result(out)
      !! Fit a consecutive-lag VAR model by equation-wise least squares.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      integer, intent(in) :: order !! Model or polynomial order.
      logical, intent(in), optional :: include_mean !! Whether to include a mean term.
      logical, intent(in), optional :: estimated(:, :) !! Flag controlling estimated.
      type(mts_var_fit_t) :: out
      integer, allocatable :: lags(:)
      integer :: i

      if (order < 1) then
         out%info = 1
         return
      end if
      lags = [(i, i=1, order)]
      out = mts_vars(series, lags, include_mean, estimated)
   end function mts_var

   pure function mts_vars(series, lags, include_mean, estimated) result(out)
      !! Fit a selected-lag VAR model with optional equation-specific masks.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      integer, intent(in) :: lags(:) !! Lags.
      logical, intent(in), optional :: include_mean !! Whether to include a mean term.
      logical, intent(in), optional :: estimated(:, :) !! Flag controlling estimated.
      type(mts_var_fit_t) :: out
      real(dp), allocatable :: design(:, :), response(:, :), xtx(:, :), inverse(:, :), beta(:)
      real(dp), allocatable :: equation_residual(:)
      logical, allocatable :: mask(:, :)
      integer, allocatable :: index(:)
      real(dp) :: rss, variance, logdet
      integer :: n, dimension, observations, columns, equation, lag_index, row, status, count

      n = size(series, 1)
      dimension = size(series, 2)
      out%includes_mean = .true.
      if (present(include_mean)) out%includes_mean = include_mean
      if (dimension < 1 .or. size(lags) < 1 .or. any(lags < 1)) then
         out%info = 1
         return
      end if
      out%max_lag = maxval(lags)
      observations = n - out%max_lag
      columns = dimension*size(lags) + merge(1, 0, out%includes_mean)
      if (observations <= columns) then
         out%info = 2
         return
      end if
      out%lags = lags
      allocate(design(observations, columns), response(observations, dimension))
      design = 0.0_dp
      if (out%includes_mean) design(:, 1) = 1.0_dp
      do row = 1, observations
         response(row, :) = series(out%max_lag + row, :)
         do lag_index = 1, size(lags)
            count = merge(1, 0, out%includes_mean) + (lag_index - 1)*dimension
            design(row, count + 1:count + dimension) = &
               series(out%max_lag + row - lags(lag_index), :)
         end do
      end do
      allocate(mask(columns, dimension))
      mask = .true.
      if (present(estimated)) then
         if (any(shape(estimated) /= [columns, dimension])) then
            out%info = 3
            return
         end if
         mask = estimated
      end if
      allocate(out%coefficients(columns, dimension), out%standard_errors(columns, dimension))
      allocate(out%residuals(observations, dimension))
      out%coefficients = 0.0_dp
      out%standard_errors = 0.0_dp
      do equation = 1, dimension
         index = pack([(row, row=1, columns)], mask(:, equation))
         out%residuals(:, equation) = response(:, equation)
         if (size(index) == 0) cycle
         xtx = matmul(transpose(design(:, index)), design(:, index))
         call invert_matrix(xtx, inverse, status)
         if (status /= 0) then
            out%info = 10 + equation
            return
         end if
         beta = matmul(inverse, matmul(transpose(design(:, index)), response(:, equation)))
         out%coefficients(index, equation) = beta
         equation_residual = response(:, equation) - matmul(design(:, index), beta)
         out%residuals(:, equation) = equation_residual
         rss = sum(equation_residual**2)
         variance = rss/real(observations - size(index), dp)
         do row = 1, size(index)
            out%standard_errors(index(row), equation) = sqrt(max(0.0_dp, variance*inverse(row, row)))
         end do
      end do
      out%sigma = matmul(transpose(out%residuals), out%residuals)/real(observations, dp)
      if (allocated(inverse)) deallocate(inverse)
      allocate(inverse(dimension, dimension))
      call inverse_logdet(out%sigma, inverse, logdet, status, 100.0_dp*epsilon(1.0_dp))
      if (status /= 0) then
         out%info = 20 + status
         return
      end if
      count = sum(merge(1, 0, mask))
      if (out%includes_mean) count = count - &
         sum(merge(1, 0, abs(out%coefficients(1, :)) > 1.0e-8_dp))
      out%aic = logdet + 2.0_dp*real(count, dp)/real(n, dp)
      out%bic = logdet + log(real(n, dp))*real(count, dp)/real(n, dp)
      out%hq = logdet + 2.0_dp*log(log(real(n, dp)))*real(count, dp)/real(n, dp)
      allocate(out%intercept(dimension), out%ar(dimension, dimension, out%max_lag))
      out%intercept = 0.0_dp
      out%ar = 0.0_dp
      count = 0
      if (out%includes_mean) then
         out%intercept = out%coefficients(1, :)
         count = 1
      end if
      do lag_index = 1, size(lags)
         out%ar(:, :, lags(lag_index)) = transpose( &
            out%coefficients(count + 1:count + dimension, :))
         count = count + dimension
      end do
   end function mts_vars

   pure function mts_var_psi(ar, max_lag) result(psi)
      !! Compute VAR MA-representation matrices from lag zero through max_lag.
      real(dp), intent(in) :: ar(:, :, :) !! Autoregressive coefficients.
      integer, intent(in) :: max_lag !! Maximum lag to consider.
      real(dp), allocatable :: psi(:, :, :)
      integer :: dimension, p, lag, j, i

      dimension = size(ar, 1)
      p = size(ar, 3)
      if (dimension < 1 .or. size(ar, 2) /= dimension .or. max_lag < 0) then
         allocate(psi(0, 0, 0))
         return
      end if
      allocate(psi(dimension, dimension, max_lag + 1))
      psi = 0.0_dp
      do i = 1, dimension
         psi(i, i, 1) = 1.0_dp
      end do
      do lag = 1, max_lag
         if (lag <= p) psi(:, :, lag + 1) = ar(:, :, lag)
         do j = 1, min(lag - 1, p)
            psi(:, :, lag + 1) = psi(:, :, lag + 1) + &
               matmul(ar(:, :, j), psi(:, :, lag - j + 1))
         end do
      end do
   end function mts_var_psi

   pure function mts_var_forecast(model, history, horizon) result(out)
      !! Forecast a fitted VAR and accumulate innovation forecast covariance.
      type(mts_var_fit_t), intent(in) :: model !! Model specification.
      real(dp), intent(in) :: history(:, :) !! History.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      type(mts_var_forecast_t) :: out
      real(dp), allocatable :: extended(:, :), psi(:, :, :)
      integer :: n, dimension, step, lag, i

      n = size(history, 1)
      dimension = size(history, 2)
      if (model%info /= 0 .or. horizon < 1 .or. n < model%max_lag .or. &
         dimension /= size(model%ar, 1)) then
         out%info = 1
         return
      end if
      allocate(extended(n + horizon, dimension), out%mean(horizon, dimension))
      allocate(out%standard_error(horizon, dimension), out%covariance(dimension, dimension, horizon))
      extended(:n, :) = history
      do step = 1, horizon
         extended(n + step, :) = model%intercept
         do lag = 1, model%max_lag
            extended(n + step, :) = extended(n + step, :) + &
               matmul(model%ar(:, :, lag), extended(n + step - lag, :))
         end do
         out%mean(step, :) = extended(n + step, :)
      end do
      psi = mts_var_psi(model%ar, horizon - 1)
      out%covariance = 0.0_dp
      do step = 1, horizon
         if (step > 1) out%covariance(:, :, step) = out%covariance(:, :, step - 1)
         out%covariance(:, :, step) = out%covariance(:, :, step) + &
            matmul(matmul(psi(:, :, step), model%sigma), transpose(psi(:, :, step)))
         do i = 1, dimension
            out%standard_error(step, i) = sqrt(max(0.0_dp, out%covariance(i, i, step)))
         end do
      end do
   end function mts_var_forecast

   pure function mts_var_irf(model, max_lag, shock_factor) result(out)
      !! Compute orthogonalized and generalized VAR impulse responses.
      type(mts_var_fit_t), intent(in) :: model !! Model specification.
      integer, intent(in) :: max_lag !! Maximum lag to consider.
      real(dp), intent(in), optional :: shock_factor(:, :) !! Shock factor.
      type(mts_var_irf_t) :: out
      real(dp), allocatable :: factor(:, :)
      integer :: dimension, lag, shock, status

      if (model%info /= 0 .or. max_lag < 0) then
         out%info = 1
         return
      end if
      dimension = size(model%sigma, 1)
      if (present(shock_factor)) then
         if (any(shape(shock_factor) /= [dimension, dimension])) then
            out%info = 2
            return
         end if
         factor = shock_factor
      else
         call cholesky_lower(model%sigma, factor, status)
         if (status /= 0) then
            out%info = 10 + status
            return
         end if
      end if
      out%psi = mts_var_psi(model%ar, max_lag)
      allocate(out%orthogonal(dimension, dimension, max_lag + 1))
      allocate(out%generalized(dimension, dimension, max_lag + 1))
      allocate(out%cumulative_orthogonal(dimension, dimension, max_lag + 1))
      allocate(out%cumulative_generalized(dimension, dimension, max_lag + 1))
      do lag = 1, max_lag + 1
         out%orthogonal(:, :, lag) = matmul(out%psi(:, :, lag), factor)
         do shock = 1, dimension
            if (model%sigma(shock, shock) <= tiny(1.0_dp)) then
               out%info = 3
               return
            end if
            out%generalized(:, shock, lag) = &
               matmul(out%psi(:, :, lag), model%sigma(:, shock))/sqrt(model%sigma(shock, shock))
         end do
         if (lag == 1) then
            out%cumulative_orthogonal(:, :, lag) = out%orthogonal(:, :, lag)
            out%cumulative_generalized(:, :, lag) = out%generalized(:, :, lag)
         else
            out%cumulative_orthogonal(:, :, lag) = &
               out%cumulative_orthogonal(:, :, lag - 1) + out%orthogonal(:, :, lag)
            out%cumulative_generalized(:, :, lag) = &
               out%cumulative_generalized(:, :, lag - 1) + out%generalized(:, :, lag)
         end if
      end do
      out%shock_factor = factor
   end function mts_var_irf

   pure function mts_fevd(model, max_horizon, generalized) result(out)
      !! Decompose VAR forecast-error variance into normalized shock contributions.
      type(mts_var_fit_t), intent(in) :: model !! Model specification.
      integer, intent(in) :: max_horizon !! Maximum horizon.
      logical, intent(in), optional :: generalized !! Flag controlling generalized.
      type(mts_fevd_t) :: out
      type(mts_var_irf_t) :: responses
      real(dp), allocatable :: numerator(:, :), total(:), covariance(:, :)
      integer :: dimension, horizon, lag, response, shock

      out%generalized = .false.
      if (present(generalized)) out%generalized = generalized
      if (model%info /= 0 .or. max_horizon < 1) then
         out%info = 1
         return
      end if
      dimension = size(model%sigma, 1)
      responses = mts_var_irf(model, max_horizon - 1)
      if (responses%info /= 0) then
         out%info = 10 + responses%info
         return
      end if
      allocate(out%contribution(dimension, dimension, max_horizon))
      allocate(out%variance(dimension, max_horizon), out%standard_error(dimension, max_horizon))
      allocate(numerator(dimension, dimension), total(dimension), covariance(dimension, dimension))
      numerator = 0.0_dp
      covariance = 0.0_dp
      do horizon = 1, max_horizon
         lag = horizon
         covariance = covariance + matmul(matmul(responses%psi(:, :, lag), model%sigma), &
            transpose(responses%psi(:, :, lag)))
         do response = 1, dimension
            out%variance(response, horizon) = covariance(response, response)
            out%standard_error(response, horizon) = sqrt(max(0.0_dp, covariance(response, response)))
            do shock = 1, dimension
               if (out%generalized) then
                  numerator(response, shock) = numerator(response, shock) + &
                     responses%generalized(response, shock, lag)**2
               else
                  numerator(response, shock) = numerator(response, shock) + &
                     responses%orthogonal(response, shock, lag)**2
               end if
            end do
         end do
         total = sum(numerator, 2)
         do response = 1, dimension
            if (total(response) <= tiny(1.0_dp)) then
               out%info = 2
               return
            end if
            out%contribution(response, :, horizon) = numerator(response, :)/total(response)
         end do
      end do
   end function mts_fevd

   pure function mts_varx(endogenous, exogenous, ar_order, exogenous_order, &
      include_mean, estimated) result(out)
      !! Fit a VARX model with contemporaneous and distributed exogenous lags.
      real(dp), intent(in) :: endogenous(:, :) !! Endogenous time-series observations.
      real(dp), intent(in) :: exogenous(:, :) !! Exogenous predictor observations.
      integer, intent(in) :: ar_order !! Autoregressive order.
      integer, intent(in) :: exogenous_order !! Exogenous order.
      logical, intent(in), optional :: include_mean !! Whether to include a mean term.
      logical, intent(in), optional :: estimated(:, :) !! Flag controlling estimated.
      type(mts_varx_fit_t) :: out
      real(dp), allocatable :: design(:, :), response(:, :), xtx(:, :), inverse(:, :), beta(:), residual(:)
      logical, allocatable :: mask(:, :)
      integer, allocatable :: index(:)
      real(dp) :: variance, logdet
      integer :: n, dimension, exogenous_dimension, first, observations, columns
      integer :: row, lag, equation, offset, i, status, parameter_count

      n = min(size(endogenous, 1), size(exogenous, 1))
      dimension = size(endogenous, 2)
      exogenous_dimension = size(exogenous, 2)
      out%includes_mean = .true.
      if (present(include_mean)) out%includes_mean = include_mean
      if (dimension < 1 .or. exogenous_dimension < 1 .or. ar_order < 0 .or. exogenous_order < 0) then
         out%info = 1
         return
      end if
      first = max(ar_order, exogenous_order) + 1
      observations = n - first + 1
      columns = merge(1, 0, out%includes_mean) + &
         (exogenous_order + 1)*exogenous_dimension + ar_order*dimension
      if (observations <= columns) then
         out%info = 2
         return
      end if
      allocate(design(observations, columns), response(observations, dimension))
      do row = 1, observations
         response(row, :) = endogenous(first + row - 1, :)
         offset = 0
         if (out%includes_mean) then
            design(row, 1) = 1.0_dp
            offset = 1
         end if
         do lag = 0, exogenous_order
            design(row, offset + 1:offset + exogenous_dimension) = &
               exogenous(first + row - 1 - lag, :)
            offset = offset + exogenous_dimension
         end do
         do lag = 1, ar_order
            design(row, offset + 1:offset + dimension) = endogenous(first + row - 1 - lag, :)
            offset = offset + dimension
         end do
      end do
      allocate(mask(columns, dimension))
      mask = .true.
      if (present(estimated)) then
         if (any(shape(estimated) /= [columns, dimension])) then
            out%info = 3
            return
         end if
         mask = estimated
      end if
      allocate(out%coefficients(columns, dimension), out%standard_errors(columns, dimension))
      allocate(out%residuals(observations, dimension))
      out%coefficients = 0.0_dp
      out%standard_errors = 0.0_dp
      do equation = 1, dimension
         index = pack([(i, i=1, columns)], mask(:, equation))
         out%residuals(:, equation) = response(:, equation)
         if (size(index) == 0) cycle
         xtx = matmul(transpose(design(:, index)), design(:, index))
         call invert_matrix(xtx, inverse, status)
         if (status /= 0) then
            out%info = 10 + equation
            return
         end if
         beta = matmul(inverse, matmul(transpose(design(:, index)), response(:, equation)))
         out%coefficients(index, equation) = beta
         residual = response(:, equation) - matmul(design(:, index), beta)
         out%residuals(:, equation) = residual
         variance = sum(residual**2)/real(observations - size(index), dp)
         do i = 1, size(index)
            out%standard_errors(index(i), equation) = sqrt(max(0.0_dp, variance*inverse(i, i)))
         end do
      end do
      out%sigma = matmul(transpose(out%residuals), out%residuals)/real(observations, dp)
      if (allocated(inverse)) deallocate(inverse)
      allocate(inverse(dimension, dimension))
      call inverse_logdet(out%sigma, inverse, logdet, status, 100.0_dp*epsilon(1.0_dp))
      if (status /= 0) then
         out%info = 20 + status
         return
      end if
      parameter_count = sum(merge(1, 0, mask))
      out%aic = logdet + 2.0_dp*real(parameter_count, dp)/real(observations, dp)
      out%bic = logdet + log(real(observations, dp))*real(parameter_count, dp)/real(observations, dp)
      out%hq = logdet + 2.0_dp*log(log(real(observations, dp)))* &
         real(parameter_count, dp)/real(observations, dp)
      out%ar_order = ar_order
      out%exogenous_order = exogenous_order
      allocate(out%intercept(dimension), out%ar(dimension, dimension, ar_order))
      allocate(out%exogenous(dimension, exogenous_dimension, 0:exogenous_order))
      out%intercept = 0.0_dp
      offset = 0
      if (out%includes_mean) then
         out%intercept = out%coefficients(1, :)
         offset = 1
      end if
      do lag = 0, exogenous_order
         out%exogenous(:, :, lag) = transpose( &
            out%coefficients(offset + 1:offset + exogenous_dimension, :))
         offset = offset + exogenous_dimension
      end do
      do lag = 1, ar_order
         out%ar(:, :, lag) = transpose(out%coefficients(offset + 1:offset + dimension, :))
         offset = offset + dimension
      end do
   end function mts_varx

   pure function mts_varx_forecast(model, endogenous_history, exogenous_history, &
      future_exogenous, horizon) result(out)
      !! Forecast VARX using supplied future exogenous observations.
      type(mts_varx_fit_t), intent(in) :: model !! Model specification.
      real(dp), intent(in) :: endogenous_history(:, :) !! Endogenous history.
      real(dp), intent(in) :: exogenous_history(:, :) !! Exogenous history.
      real(dp), intent(in) :: future_exogenous(:, :) !! Future exogenous.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      type(mts_var_forecast_t) :: out
      real(dp), allocatable :: extended_y(:, :), extended_x(:, :), psi(:, :, :)
      integer :: n, dimension, exogenous_dimension, step, lag, i

      n = size(endogenous_history, 1)
      dimension = size(endogenous_history, 2)
      exogenous_dimension = size(exogenous_history, 2)
      if (model%info /= 0 .or. horizon < 1 .or. size(future_exogenous, 1) < horizon .or. &
         size(future_exogenous, 2) /= exogenous_dimension .or. n < model%ar_order .or. &
         n < model%exogenous_order .or. size(exogenous_history, 1) /= n .or. &
         dimension /= size(model%ar, 1) .or. &
         exogenous_dimension /= size(model%exogenous, 2)) then
         out%info = 1
         return
      end if
      allocate(extended_y(n + horizon, dimension), extended_x(n + horizon, exogenous_dimension))
      allocate(out%mean(horizon, dimension), out%standard_error(horizon, dimension))
      allocate(out%covariance(dimension, dimension, horizon))
      extended_y(:n, :) = endogenous_history
      extended_x(:n, :) = exogenous_history
      extended_x(n + 1:, :) = future_exogenous(:horizon, :)
      do step = 1, horizon
         extended_y(n + step, :) = model%intercept
         do lag = 1, model%ar_order
            extended_y(n + step, :) = extended_y(n + step, :) + &
               matmul(model%ar(:, :, lag), extended_y(n + step - lag, :))
         end do
         do lag = 0, model%exogenous_order
            extended_y(n + step, :) = extended_y(n + step, :) + &
               matmul(model%exogenous(:, :, lag), extended_x(n + step - lag, :))
         end do
         out%mean(step, :) = extended_y(n + step, :)
      end do
      psi = mts_var_psi(model%ar, horizon - 1)
      out%covariance = 0.0_dp
      do step = 1, horizon
         if (step > 1) out%covariance(:, :, step) = out%covariance(:, :, step - 1)
         out%covariance(:, :, step) = out%covariance(:, :, step) + &
            matmul(matmul(psi(:, :, step), model%sigma), transpose(psi(:, :, step)))
         do i = 1, dimension
            out%standard_error(step, i) = sqrt(max(0.0_dp, out%covariance(i, i, step)))
         end do
      end do
   end function mts_varx_forecast

   pure function mts_varx_irf(model, max_lag) result(out)
      !! Compute endogenous and exogenous VARX dynamic responses.
      type(mts_varx_fit_t), intent(in) :: model !! Model specification.
      integer, intent(in) :: max_lag !! Maximum lag to consider.
      type(mts_varx_irf_t) :: out
      real(dp), allocatable :: psi(:, :, :)
      integer :: dimension, exogenous_dimension, lag, j

      if (model%info /= 0 .or. max_lag < 0) then
         out%info = 1
         return
      end if
      dimension = size(model%ar, 1)
      exogenous_dimension = size(model%exogenous, 2)
      psi = mts_var_psi(model%ar, max_lag)
      out%endogenous = psi
      allocate(out%exogenous(dimension, exogenous_dimension, max_lag + 1))
      allocate(out%cumulative_endogenous(dimension, dimension, max_lag + 1))
      allocate(out%cumulative_exogenous(dimension, exogenous_dimension, max_lag + 1))
      out%exogenous = 0.0_dp
      do lag = 0, max_lag
         if (lag <= model%exogenous_order) out%exogenous(:, :, lag + 1) = &
            model%exogenous(:, :, lag)
         do j = 1, min(lag, model%ar_order)
            out%exogenous(:, :, lag + 1) = out%exogenous(:, :, lag + 1) + &
               matmul(model%ar(:, :, j), out%exogenous(:, :, lag - j + 1))
         end do
         if (lag == 0) then
            out%cumulative_endogenous(:, :, 1) = out%endogenous(:, :, 1)
            out%cumulative_exogenous(:, :, 1) = out%exogenous(:, :, 1)
         else
            out%cumulative_endogenous(:, :, lag + 1) = &
               out%cumulative_endogenous(:, :, lag) + out%endogenous(:, :, lag + 1)
            out%cumulative_exogenous(:, :, lag + 1) = &
               out%cumulative_exogenous(:, :, lag) + out%exogenous(:, :, lag + 1)
         end if
      end do
   end function mts_varx_irf

   pure function mts_varx_order(endogenous, exogenous, max_ar_order, max_exogenous_order) result(out)
      !! Select endogenous and exogenous VARX lag orders by information criteria.
      real(dp), intent(in) :: endogenous(:, :) !! Endogenous time-series observations.
      real(dp), intent(in) :: exogenous(:, :) !! Exogenous predictor observations.
      integer, intent(in) :: max_ar_order !! Maximum autoregressive order.
      integer, intent(in) :: max_exogenous_order !! Maximum exogenous order.
      type(mts_varx_order_t) :: out
      type(mts_varx_fit_t) :: fitted
      integer :: p, m, location(2)

      if (max_ar_order < 0 .or. max_exogenous_order < 0) then
         out%info = 1
         return
      end if
      allocate(out%aic(0:max_ar_order, 0:max_exogenous_order))
      allocate(out%bic(0:max_ar_order, 0:max_exogenous_order))
      allocate(out%hq(0:max_ar_order, 0:max_exogenous_order))
      do p = 0, max_ar_order
         do m = 0, max_exogenous_order
            fitted = mts_varx(endogenous, exogenous, p, m)
            if (fitted%info /= 0) then
               out%info = 10 + p*(max_exogenous_order + 1) + m
               return
            end if
            out%aic(p, m) = fitted%aic
            out%bic(p, m) = fitted%bic
            out%hq(p, m) = fitted%hq
         end do
      end do
      location = minloc(out%aic) - 1
      out%aic_order = location
      location = minloc(out%bic) - 1
      out%bic_order = location
      location = minloc(out%hq) - 1
      out%hq_order = location
   end function mts_varx_order

   pure function mts_varma_psi(model, max_lag) result(psi)
      !! Compute MTS VARMA PSI weights from lag zero through max_lag.
      type(mts_varma_model_t), intent(in) :: model !! Model specification.
      integer, intent(in) :: max_lag !! Maximum lag to consider.
      real(dp), allocatable :: psi(:, :, :)
      integer :: dimension, p, q, lag, j, i

      if (.not. allocated(model%ar) .or. .not. allocated(model%ma) .or. max_lag < 0) then
         allocate(psi(0, 0, 0))
         return
      end if
      dimension = size(model%ar, 1)
      p = size(model%ar, 3)
      q = size(model%ma, 3)
      if (dimension < 1 .or. size(model%ar, 2) /= dimension .or. &
         size(model%ma, 1) /= dimension .or. size(model%ma, 2) /= dimension) then
         allocate(psi(0, 0, 0))
         return
      end if
      allocate(psi(dimension, dimension, max_lag + 1))
      psi = 0.0_dp
      do i = 1, dimension
         psi(i, i, 1) = 1.0_dp
      end do
      do lag = 1, max_lag
         if (lag <= q) psi(:, :, lag + 1) = -model%ma(:, :, lag)
         do j = 1, min(lag, p)
            psi(:, :, lag + 1) = psi(:, :, lag + 1) + &
               matmul(model%ar(:, :, j), psi(:, :, lag - j + 1))
         end do
      end do
   end function mts_varma_psi

   pure function mts_varma_pi(model, max_lag) result(pi_weights)
      !! Compute MTS VARMA infinite autoregressive PI weights.
      type(mts_varma_model_t), intent(in) :: model !! Model specification.
      integer, intent(in) :: max_lag !! Maximum lag to consider.
      real(dp), allocatable :: pi_weights(:, :, :)
      integer :: dimension, p, q, lag, j, i

      if (.not. allocated(model%ar) .or. .not. allocated(model%ma) .or. max_lag < 0) then
         allocate(pi_weights(0, 0, 0))
         return
      end if
      dimension = size(model%ar, 1)
      p = size(model%ar, 3)
      q = size(model%ma, 3)
      allocate(pi_weights(dimension, dimension, max_lag + 1))
      pi_weights = 0.0_dp
      do i = 1, dimension
         pi_weights(i, i, 1) = 1.0_dp
      end do
      do lag = 1, max_lag
         if (lag <= p) pi_weights(:, :, lag + 1) = model%ar(:, :, lag)
         do j = 1, min(lag, q)
            pi_weights(:, :, lag + 1) = pi_weights(:, :, lag + 1) + &
               matmul(model%ma(:, :, j), pi_weights(:, :, lag - j + 1))
         end do
      end do
   end function mts_varma_pi

   pure function mts_varma_irf(model, max_lag, shock_factor) result(out)
      !! Compute raw, orthogonalized, generalized, and cumulative VARMA responses.
      type(mts_varma_model_t), intent(in) :: model !! Model specification.
      integer, intent(in) :: max_lag !! Maximum lag to consider.
      real(dp), intent(in), optional :: shock_factor(:, :) !! Shock factor.
      type(mts_var_irf_t) :: out
      real(dp), allocatable :: factor(:, :)
      integer :: dimension, lag, shock, status

      if (model%info /= 0 .or. .not. allocated(model%sigma) .or. max_lag < 0) then
         out%info = 1
         return
      end if
      dimension = size(model%sigma, 1)
      if (present(shock_factor)) then
         if (any(shape(shock_factor) /= [dimension, dimension])) then
            out%info = 2
            return
         end if
         factor = shock_factor
      else
         call cholesky_lower(model%sigma, factor, status)
         if (status /= 0) then
            out%info = 10 + status
            return
         end if
      end if
      out%psi = mts_varma_psi(model, max_lag)
      if (size(out%psi, 1) /= dimension) then
         out%info = 3
         return
      end if
      allocate(out%orthogonal(dimension, dimension, max_lag + 1))
      allocate(out%generalized(dimension, dimension, max_lag + 1))
      allocate(out%cumulative_orthogonal(dimension, dimension, max_lag + 1))
      allocate(out%cumulative_generalized(dimension, dimension, max_lag + 1))
      do lag = 1, max_lag + 1
         out%orthogonal(:, :, lag) = matmul(out%psi(:, :, lag), factor)
         do shock = 1, dimension
            out%generalized(:, shock, lag) = matmul(out%psi(:, :, lag), &
               model%sigma(:, shock))/sqrt(model%sigma(shock, shock))
         end do
         if (lag == 1) then
            out%cumulative_orthogonal(:, :, lag) = out%orthogonal(:, :, lag)
            out%cumulative_generalized(:, :, lag) = out%generalized(:, :, lag)
         else
            out%cumulative_orthogonal(:, :, lag) = out%cumulative_orthogonal(:, :, lag - 1) + &
               out%orthogonal(:, :, lag)
            out%cumulative_generalized(:, :, lag) = out%cumulative_generalized(:, :, lag - 1) + &
               out%generalized(:, :, lag)
         end if
      end do
      out%shock_factor = factor
   end function mts_varma_irf

   pure function mts_varma_covariance(model, max_lag, truncation) result(out)
      !! Approximate theoretical VARMA covariance and correlation by PSI truncation.
      type(mts_varma_model_t), intent(in) :: model !! Model specification.
      integer, intent(in) :: max_lag !! Maximum lag to consider.
      integer, intent(in), optional :: truncation !! Truncation.
      type(mts_varma_covariance_t) :: out
      real(dp), allocatable :: psi(:, :, :)
      real(dp) :: scale_i, scale_j
      integer :: limit, dimension, lag, j, i, k

      limit = 120
      if (present(truncation)) limit = truncation
      if (model%info /= 0 .or. max_lag < 0 .or. limit < max_lag) then
         out%info = 1
         return
      end if
      dimension = size(model%sigma, 1)
      psi = mts_varma_psi(model, limit)
      allocate(out%covariance(dimension, dimension, max_lag + 1))
      allocate(out%correlation(dimension, dimension, max_lag + 1))
      out%covariance = 0.0_dp
      do lag = 0, max_lag
         do j = 0, limit - lag
            out%covariance(:, :, lag + 1) = out%covariance(:, :, lag + 1) + &
               matmul(matmul(psi(:, :, j + lag + 1), model%sigma), transpose(psi(:, :, j + 1)))
         end do
      end do
      do lag = 0, max_lag
         do i = 1, dimension
            scale_i = sqrt(out%covariance(i, i, 1))
            do k = 1, dimension
               scale_j = sqrt(out%covariance(k, k, 1))
               out%correlation(i, k, lag + 1) = &
                  out%covariance(i, k, lag + 1)/(scale_i*scale_j)
            end do
         end do
      end do
      out%truncation = limit
   end function mts_varma_covariance

   pure function mts_varma_simulate_from_innovations(model, innovations, burnin) result(out)
      !! Simulate VARMA observations from caller-supplied innovations.
      type(mts_varma_model_t), intent(in) :: model !! Model specification.
      real(dp), intent(in) :: innovations(:, :) !! Model innovations.
      integer, intent(in), optional :: burnin !! Number of initial simulation draws to discard.
      type(mts_varma_simulation_t) :: out
      real(dp), allocatable :: work(:, :)
      integer :: discard, n, dimension, p, q, t, lag

      discard = 0
      if (present(burnin)) discard = burnin
      n = size(innovations, 1)
      dimension = size(innovations, 2)
      p = size(model%ar, 3)
      q = size(model%ma, 3)
      if (model%info /= 0 .or. discard < 0 .or. discard >= n .or. &
         dimension /= size(model%ar, 1)) then
         out%info = 1
         return
      end if
      allocate(work(n, dimension))
      work = 0.0_dp
      do t = 1, n
         work(t, :) = model%intercept + innovations(t, :)
         do lag = 1, min(p, t - 1)
            work(t, :) = work(t, :) + matmul(model%ar(:, :, lag), work(t - lag, :))
         end do
         do lag = 1, min(q, t - 1)
            work(t, :) = work(t, :) - matmul(model%ma(:, :, lag), innovations(t - lag, :))
         end do
      end do
      out%series = work(discard + 1:, :)
      out%innovations = innovations(discard + 1:, :)
   end function mts_varma_simulate_from_innovations

   function mts_varma_simulate(model, observation_count, burnin) result(out)
      !! Simulate Gaussian VARMA observations using the shared random stream.
      type(mts_varma_model_t), intent(in) :: model !! Model specification.
      integer, intent(in) :: observation_count !! Observation count.
      integer, intent(in), optional :: burnin !! Number of initial simulation draws to discard.
      type(mts_varma_simulation_t) :: out
      real(dp), allocatable :: innovations(:, :), zero(:)
      integer :: discard, total, dimension, t, status

      discard = 200
      if (present(burnin)) discard = burnin
      dimension = size(model%sigma, 1)
      total = observation_count + discard
      if (observation_count < 1 .or. discard < 0) then
         out%info = 1
         return
      end if
      allocate(innovations(total, dimension), zero(dimension))
      zero = 0.0_dp
      do t = 1, total
         call random_multivariate_normal(zero, model%sigma, innovations(t, :), status)
         if (status /= 0) then
            out%info = 10 + status
            return
         end if
      end do
      out = mts_varma_simulate_from_innovations(model, innovations, discard)
   end function mts_varma_simulate

   pure function mts_varma_forecast(model, history, residuals, horizon) result(out)
      !! Forecast VARMA observations using historical innovations and zero future shocks.
      type(mts_varma_model_t), intent(in) :: model !! Model specification.
      real(dp), intent(in) :: history(:, :) !! History.
      real(dp), intent(in) :: residuals(:, :) !! Model residuals.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      type(mts_var_forecast_t) :: out
      real(dp), allocatable :: extended(:, :), psi(:, :, :)
      integer :: n, dimension, p, q, step, lag, i

      n = size(history, 1)
      dimension = size(history, 2)
      p = size(model%ar, 3)
      q = size(model%ma, 3)
      if (model%info /= 0 .or. horizon < 1 .or. size(residuals, 1) /= n .or. &
         size(residuals, 2) /= dimension .or. n < max(p, q)) then
         out%info = 1
         return
      end if
      allocate(extended(n + horizon, dimension), out%mean(horizon, dimension))
      allocate(out%standard_error(horizon, dimension), out%covariance(dimension, dimension, horizon))
      extended(:n, :) = history
      do step = 1, horizon
         extended(n + step, :) = model%intercept
         do lag = 1, p
            extended(n + step, :) = extended(n + step, :) + &
               matmul(model%ar(:, :, lag), extended(n + step - lag, :))
         end do
         do lag = 1, q
            if (n + step - lag <= n) extended(n + step, :) = extended(n + step, :) - &
               matmul(model%ma(:, :, lag), residuals(n + step - lag, :))
         end do
         out%mean(step, :) = extended(n + step, :)
      end do
      psi = mts_varma_psi(model, horizon - 1)
      out%covariance = 0.0_dp
      do step = 1, horizon
         if (step > 1) out%covariance(:, :, step) = out%covariance(:, :, step - 1)
         out%covariance(:, :, step) = out%covariance(:, :, step) + &
            matmul(matmul(psi(:, :, step), model%sigma), transpose(psi(:, :, step)))
         do i = 1, dimension
            out%standard_error(step, i) = sqrt(max(0.0_dp, out%covariance(i, i, step)))
         end do
      end do
   end function mts_varma_forecast

   pure function mts_vecm_fit(series, level_order, rank, include_constant, &
      cointegration, estimated) result(out)
      !! Fit a VECM with supplied or Johansen-estimated cointegration vectors.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      integer, intent(in) :: level_order !! Level order.
      integer, intent(in) :: rank !! Matrix or cointegration rank.
      logical, intent(in), optional :: include_constant !! Whether to include the constant.
      logical, intent(in), optional :: estimated(:, :) !! Flag controlling estimated.
      real(dp), intent(in), optional :: cointegration(:, :) !! Cointegration.
      type(mts_vecm_fit_t) :: out
      type(johansen_result_t) :: johansen
      real(dp), allocatable :: differences(:, :), design(:, :), response(:, :)
      real(dp), allocatable :: xtx(:, :), inverse(:, :), xty(:), beta(:), covariance_inverse(:, :)
      integer, allocatable :: selected(:)
      real(dp) :: logdet
      integer :: n, dimension, rows, regressors, t, lag, column, equation, count_selected, status, i
      integer :: parameter_count
      logical :: use_constant

      n = size(series, 1)
      dimension = size(series, 2)
      use_constant = .false.
      if (present(include_constant)) use_constant = include_constant
      if (dimension < 2 .or. rank < 1 .or. rank >= dimension .or. level_order < 1 .or. &
         n <= level_order + dimension + 1) then
         out%info = 1
         return
      end if
      if (present(cointegration)) then
         if (size(cointegration, 1) /= dimension .or. size(cointegration, 2) < rank) then
            out%info = 2
            return
         end if
         out%cointegration = cointegration(:, :rank)
      else
         if (level_order < 2) then
            out%info = 3
            return
         end if
         johansen = johansen_test(series, 'trace', 'none', level_order, 'longrun')
         if (johansen%info /= 0) then
            out%info = 10 + johansen%info
            return
         end if
         out%cointegration = johansen%cointegration(:dimension, :rank)
      end if
      rows = n - level_order
      regressors = rank + merge(1, 0, use_constant) + dimension*(level_order - 1)
      allocate(differences(n - 1, dimension), design(rows, regressors), response(rows, dimension))
      differences = series(2:, :) - series(:n - 1, :)
      do i = 1, rows
         t = level_order + i
         design(i, :rank) = matmul(series(t - 1, :), out%cointegration)
         column = rank
         if (use_constant) then
            column = column + 1
            design(i, column) = 1.0_dp
         end if
         do lag = 1, level_order - 1
            design(i, column + 1:column + dimension) = differences(t - 1 - lag, :)
            column = column + dimension
         end do
         response(i, :) = differences(t - 1, :)
      end do
      allocate(out%coefficients(regressors, dimension), out%standard_errors(regressors, dimension))
      allocate(out%estimated(regressors, dimension), out%residuals(rows, dimension))
      out%coefficients = 0.0_dp
      out%standard_errors = 0.0_dp
      out%estimated = .true.
      if (present(estimated)) then
         if (any(shape(estimated) /= [regressors, dimension])) then
            out%info = 4
            return
         end if
         out%estimated = estimated
      end if
      parameter_count = 0
      do equation = 1, dimension
         count_selected = count(out%estimated(:, equation))
         parameter_count = parameter_count + count_selected
         if (count_selected == 0) then
            out%residuals(:, equation) = response(:, equation)
            cycle
         end if
         allocate(selected(count_selected))
         selected = pack([(i, i=1, regressors)], out%estimated(:, equation))
         xtx = matmul(transpose(design(:, selected)), design(:, selected))
         call invert_matrix(xtx, inverse, status)
         if (status /= 0) then
            out%info = 20 + equation
            return
         end if
         xty = matmul(transpose(design(:, selected)), response(:, equation))
         beta = matmul(inverse, xty)
         out%coefficients(selected, equation) = beta
         out%residuals(:, equation) = response(:, equation) - matmul(design(:, selected), beta)
         out%standard_errors(selected, equation) = sqrt(max(0.0_dp, &
            sum(out%residuals(:, equation)**2)/real(rows, dp))* &
            [(inverse(i, i), i=1, count_selected)])
         deallocate(selected)
      end do
      allocate(out%loading(dimension, rank), out%intercept(dimension))
      allocate(out%gamma(dimension, dimension, level_order - 1))
      out%loading = transpose(out%coefficients(:rank, :))
      out%intercept = 0.0_dp
      column = rank
      if (use_constant) then
         column = column + 1
         out%intercept = out%coefficients(column, :)
      end if
      do lag = 1, level_order - 1
         out%gamma(:, :, lag) = transpose(out%coefficients(&
            column + 1:column + dimension, :))
         column = column + dimension
      end do
      out%sigma = matmul(transpose(out%residuals), out%residuals)/real(rows, dp)
      allocate(covariance_inverse(dimension, dimension))
      call inverse_logdet(out%sigma, covariance_inverse, logdet, status, &
         100.0_dp*epsilon(1.0_dp))
      if (status /= 0) then
         out%info = 30 + status
         return
      end if
      out%aic = logdet + 2.0_dp*real(parameter_count, dp)/real(n, dp)
      out%bic = logdet + log(real(n, dp))*real(parameter_count, dp)/real(n, dp)
      call build_level_var(out)
      out%rank = rank
      out%level_order = level_order
      out%includes_constant = use_constant

   contains

      pure subroutine build_level_var(fitted)
         !! Convert VECM coefficients to an equivalent level VAR.
         type(mts_vecm_fit_t), intent(inout) :: fitted !! Fitted, updated in place.
         real(dp) :: identity(dimension, dimension), pi_matrix(dimension, dimension)
         integer :: current_lag

         identity = 0.0_dp
         do current_lag = 1, dimension
            identity(current_lag, current_lag) = 1.0_dp
         end do
         pi_matrix = matmul(fitted%loading, transpose(fitted%cointegration))
         allocate(fitted%level_var%ar(dimension, dimension, level_order))
         allocate(fitted%level_var%intercept(dimension), fitted%level_var%sigma(dimension, dimension))
         fitted%level_var%ar = 0.0_dp
         fitted%level_var%ar(:, :, 1) = identity + pi_matrix
         if (level_order > 1) then
            fitted%level_var%ar(:, :, 1) = fitted%level_var%ar(:, :, 1) + fitted%gamma(:, :, 1)
            do current_lag = 2, level_order - 1
               fitted%level_var%ar(:, :, current_lag) = fitted%gamma(:, :, current_lag) - &
                  fitted%gamma(:, :, current_lag - 1)
            end do
            fitted%level_var%ar(:, :, level_order) = -fitted%gamma(:, :, level_order - 1)
         end if
         fitted%level_var%intercept = fitted%intercept
         fitted%level_var%sigma = fitted%sigma
         fitted%level_var%residuals = fitted%residuals
         fitted%level_var%lags = [(current_lag, current_lag=1, level_order)]
         fitted%level_var%max_lag = level_order
         fitted%level_var%includes_mean = use_constant
      end subroutine build_level_var
   end function mts_vecm_fit

   pure function mts_vecm_forecast(fitted, history, horizon) result(out)
      !! Forecast VECM levels and first differences through its level-VAR form.
      type(mts_vecm_fit_t), intent(in) :: fitted !! Fitted.
      real(dp), intent(in) :: history(:, :) !! History.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      type(mts_vecm_forecast_t) :: out
      integer :: step

      if (fitted%info /= 0 .or. horizon < 1) then
         out%info = 1
         return
      end if
      out%level = mts_var_forecast(fitted%level_var, history, horizon)
      if (out%level%info /= 0) then
         out%info = 10 + out%level%info
         return
      end if
      allocate(out%difference(horizon, size(history, 2)))
      out%difference(1, :) = out%level%mean(1, :) - history(size(history, 1), :)
      do step = 2, horizon
         out%difference(step, :) = out%level%mean(step, :) - out%level%mean(step - 1, :)
      end do
   end function mts_vecm_forecast

   pure function mts_factor_fit(series, factors, max_factors, standardize) result(out)
      !! Extract principal-component factors and Bai-Ng factor-count criteria.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      integer, intent(in), optional :: factors !! Factors.
      integer, intent(in), optional :: max_factors !! Maximum factors.
      logical, intent(in), optional :: standardize !! Flag controlling standardize.
      type(mts_factor_model_t) :: out
      real(dp), allocatable :: normalized(:, :), covariance(:, :), vectors(:, :)
      real(dp), allocatable :: trial_scores(:, :), trial_common(:, :)
      real(dp) :: variance, penalty1, penalty2, penalty3, total_variance
      integer :: n, variables, maximum, selected, count_factors, status, i
      logical :: use_scaling

      n = size(series, 1)
      variables = size(series, 2)
      use_scaling = .true.
      if (present(standardize)) use_scaling = standardize
      if (n < 3 .or. variables < 2) then
         out%info = 1
         return
      end if
      maximum = min(variables - 1, max(1, min(n - 1, 10)))
      if (present(max_factors)) maximum = min(variables - 1, min(n - 1, max_factors))
      if (maximum < 1) then
         out%info = 2
         return
      end if
      allocate(out%mean(variables), out%scale(variables), normalized(n, variables))
      out%mean = sum(series, 1)/real(n, dp)
      normalized = series - spread(out%mean, 1, n)
      out%scale = 1.0_dp
      if (use_scaling) then
         do i = 1, variables
            out%scale(i) = sqrt(sum(normalized(:, i)**2)/real(n - 1, dp))
            if (out%scale(i) <= sqrt(epsilon(1.0_dp))) then
               out%info = 3
               return
            end if
            normalized(:, i) = normalized(:, i)/out%scale(i)
         end do
      end if
      covariance = matmul(transpose(normalized), normalized)/real(n, dp)
      call symmetric_eigen(covariance, out%eigenvalues, vectors, status)
      if (status /= 0) then
         out%info = 10 + status
         return
      end if
      out%eigenvalues = max(0.0_dp, out%eigenvalues)
      total_variance = sum(out%eigenvalues)
      allocate(out%explained(variables), out%ic1(0:maximum), out%ic2(0:maximum), &
         out%ic3(0:maximum))
      if (total_variance > tiny(1.0_dp)) then
         out%explained = out%eigenvalues/total_variance
      else
         out%info = 4
         return
      end if
      penalty1 = real(n + variables, dp)/real(n*variables, dp)* &
         log(real(n*variables, dp)/real(n + variables, dp))
      penalty2 = real(n + variables, dp)/real(n*variables, dp)* &
         log(real(min(n, variables), dp))
      penalty3 = log(real(min(n, variables), dp))/real(min(n, variables), dp)
      do count_factors = 0, maximum
         if (count_factors == 0) then
            variance = sum(normalized**2)/real(n*variables, dp)
         else
            trial_scores = matmul(normalized, vectors(:, :count_factors))
            trial_common = matmul(trial_scores, transpose(vectors(:, :count_factors)))
            variance = sum((normalized - trial_common)**2)/real(n*variables, dp)
         end if
         variance = max(variance, tiny(1.0_dp))
         out%ic1(count_factors) = log(variance) + real(count_factors, dp)*penalty1
         out%ic2(count_factors) = log(variance) + real(count_factors, dp)*penalty2
         out%ic3(count_factors) = log(variance) + real(count_factors, dp)*penalty3
      end do
      out%ic1_factors = minloc(out%ic1, dim=1) - 1
      out%ic2_factors = minloc(out%ic2, dim=1) - 1
      out%ic3_factors = minloc(out%ic3, dim=1) - 1
      selected = max(1, out%ic2_factors)
      if (present(factors)) selected = factors
      if (selected < 1 .or. selected > maximum) then
         out%info = 5
         return
      end if
      out%loadings = vectors(:, :selected)
      out%scores = matmul(normalized, out%loadings)
      out%common = matmul(out%scores, transpose(out%loadings))
      out%common = spread(out%mean, 1, n) + out%common*spread(out%scale, 1, n)
      out%residuals = series - out%common
      out%factors = selected
      out%standardized = use_scaling
   end function mts_factor_fit

   pure function mts_factor_forecast(model, var_order, horizon) result(out)
      !! Fit factor dynamics by VAR and reconstruct forecasts in original units.
      type(mts_factor_model_t), intent(in) :: model !! Model specification.
      integer, intent(in) :: var_order !! Var order.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      type(mts_factor_forecast_t) :: out
      real(dp), allocatable :: idiosyncratic_covariance(:, :), scaled_loadings(:, :)
      integer :: variables, factors, step, i

      variables = size(model%loadings, 1)
      factors = size(model%loadings, 2)
      if (model%info /= 0 .or. var_order < 1 .or. horizon < 1 .or. &
         size(model%scores, 1) <= var_order + 1) then
         out%info = 1
         return
      end if
      out%factor_var = mts_var(model%scores, var_order)
      if (out%factor_var%info /= 0) then
         out%info = 10 + out%factor_var%info
         return
      end if
      out%factor_forecast = mts_var_forecast(out%factor_var, model%scores, horizon)
      if (out%factor_forecast%info /= 0) then
         out%info = 20 + out%factor_forecast%info
         return
      end if
      allocate(scaled_loadings(variables, factors))
      scaled_loadings = model%loadings*spread(model%scale, 2, factors)
      allocate(idiosyncratic_covariance(variables, variables))
      idiosyncratic_covariance = matmul(transpose(model%residuals), model%residuals)/ &
         real(size(model%residuals, 1), dp)
      allocate(out%mean(horizon, variables), out%standard_error(horizon, variables))
      allocate(out%covariance(variables, variables, horizon))
      do step = 1, horizon
         out%mean(step, :) = model%mean + matmul(out%factor_forecast%mean(step, :), &
            transpose(scaled_loadings))
         out%covariance(:, :, step) = matmul(scaled_loadings, &
            matmul(out%factor_forecast%covariance(:, :, step), transpose(scaled_loadings))) + &
            idiosyncratic_covariance
         do i = 1, variables
            out%standard_error(step, i) = sqrt(max(0.0_dp, out%covariance(i, i, step)))
         end do
      end do
   end function mts_factor_forecast

   pure function mts_constrained_factor_fit(series, constraint, factors) result(out)
      !! Fit the MTS hfactor weighted least-squares constrained factor model.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      real(dp), intent(in) :: constraint(:, :) !! Constraint.
      integer, intent(in) :: factors !! Factors.
      type(mts_constrained_factor_t) :: out
      real(dp), allocatable :: normalized(:, :), gram(:, :), eigenvalues(:), eigenvectors(:, :)
      real(dp), allocatable :: inverse_sqrt(:, :), projected(:, :), covariance(:, :)
      real(dp), allocatable :: directions(:, :), scores(:, :), constrained_loadings(:, :)
      real(dp) :: score_scale, tolerance
      integer :: n, variables, constraints, status, i

      n = size(series, 1)
      variables = size(series, 2)
      constraints = size(constraint, 2)
      if (n < 3 .or. variables < 2 .or. size(constraint, 1) /= variables .or. &
         factors < 1 .or. factors > constraints) then
         out%info = 1
         return
      end if
      allocate(out%factor_model%mean(variables), out%factor_model%scale(variables))
      allocate(normalized(n, variables))
      out%factor_model%mean = sum(series, 1)/real(n, dp)
      normalized = series - spread(out%factor_model%mean, 1, n)
      do i = 1, variables
         out%factor_model%scale(i) = sqrt(sum(normalized(:, i)**2)/real(n - 1, dp))
         if (out%factor_model%scale(i) <= sqrt(epsilon(1.0_dp))) then
            out%info = 2
            return
         end if
         normalized(:, i) = normalized(:, i)/out%factor_model%scale(i)
      end do
      gram = matmul(transpose(constraint), constraint)
      call symmetric_eigen(gram, eigenvalues, eigenvectors, status)
      if (status /= 0) then
         out%info = 10 + status
         return
      end if
      tolerance = 100.0_dp*epsilon(1.0_dp)*max(1.0_dp, maxval(eigenvalues))
      if (minval(eigenvalues) <= tolerance) then
         out%info = 3
         return
      end if
      allocate(inverse_sqrt(constraints, constraints))
      inverse_sqrt = matmul(eigenvectors, matmul(&
         diag_values(1.0_dp/sqrt(eigenvalues)), transpose(eigenvectors)))
      projected = matmul(matmul(normalized, constraint), inverse_sqrt)
      covariance = matmul(transpose(projected), projected)/real(n, dp)
      call symmetric_eigen(covariance, out%constrained_eigenvalues, directions, status)
      if (status /= 0) then
         out%info = 20 + status
         return
      end if
      scores = matmul(projected, directions(:, :factors))
      do i = 1, factors
         score_scale = sqrt(sum(scores(:, i)**2)/real(n - 1, dp))
         if (score_scale <= sqrt(epsilon(1.0_dp))) then
            out%info = 4
            return
         end if
         scores(:, i) = scores(:, i)/score_scale
      end do
      out%omega = matmul(matmul(inverse_sqrt, inverse_sqrt), &
         matmul(transpose(matmul(normalized, constraint)), scores))/real(n, dp)
      constrained_loadings = matmul(constraint, out%omega)
      out%constraint = constraint
      out%factor_model%loadings = constrained_loadings
      out%factor_model%scores = scores
      out%factor_model%common = matmul(scores, transpose(constrained_loadings))
      out%factor_model%common = spread(out%factor_model%mean, 1, n) + &
         out%factor_model%common*spread(out%factor_model%scale, 1, n)
      out%factor_model%residuals = series - out%factor_model%common
      out%factor_model%factors = factors
      out%factor_model%standardized = .true.
      out%psi = matmul(transpose(normalized), normalized)/real(n, dp) - &
         matmul(constrained_loadings, transpose(constrained_loadings))
      out%explained = sum(constrained_loadings**2)/real(variables, dp)

   contains

      pure function diag_values(values) result(matrix)
         !! Form a diagonal matrix from a vector.
         real(dp), intent(in) :: values(:) !! Input values.
         real(dp) :: matrix(size(values), size(values))
         integer :: diagonal

         matrix = 0.0_dp
         do diagonal = 1, size(values)
            matrix(diagonal, diagonal) = values(diagonal)
         end do
      end function diag_values
   end function mts_constrained_factor_fit

   pure function mts_constrained_factor_forecast(model, var_order, horizon) result(out)
      !! Forecast a constrained factor model through its normalized factor scores.
      type(mts_constrained_factor_t), intent(in) :: model !! Model specification.
      integer, intent(in) :: var_order !! Var order.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      type(mts_factor_forecast_t) :: out

      if (model%info /= 0) then
         out%info = 1
         return
      end if
      out = mts_factor_forecast(model%factor_model, var_order, horizon)
   end function mts_constrained_factor_forecast

   pure function mts_bvar_fit(series, order, prior, include_mean) result(out)
      !! Estimate the conjugate matrix-normal/inverse-Wishart MTS BVAR.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      integer, intent(in) :: order !! Model or polynomial order.
      type(mts_bvar_prior_t), intent(in) :: prior !! Prior-distribution specification.
      logical, intent(in), optional :: include_mean !! Whether to include a mean term.
      type(mts_bvar_fit_t) :: out
      real(dp), allocatable :: design(:, :), response(:, :), cross_x(:, :), cross_y(:, :)
      real(dp), allocatable :: posterior_inverse(:, :), weighted_mean(:, :), difference(:, :)
      real(dp) :: denominator
      integer :: n, variables, rows, regressors, offset, lag, t, status, equation, other, i
      logical :: use_mean

      n = size(series, 1)
      variables = size(series, 2)
      use_mean = .true.
      if (present(include_mean)) use_mean = include_mean
      regressors = variables*order + merge(1, 0, use_mean)
      rows = n - order
      if (prior%info /= 0 .or. .not. allocated(prior%mean) .or. &
         .not. allocated(prior%precision) .or. .not. allocated(prior%scale)) then
         out%info = 1
         return
      end if
      if (order < 1 .or. variables < 1 .or. rows <= variables + 1 .or. &
         any(shape(prior%mean) /= [regressors, variables]) .or. &
         any(shape(prior%precision) /= [regressors, regressors]) .or. &
         any(shape(prior%scale) /= [variables, variables])) then
         out%info = 1
         return
      end if
      denominator = prior%degrees_of_freedom + real(rows - variables - 1, dp)
      if (denominator <= 0.0_dp) then
         out%info = 2
         return
      end if
      allocate(design(rows, regressors), response(rows, variables))
      design = 0.0_dp
      offset = 0
      if (use_mean) then
         design(:, 1) = 1.0_dp
         offset = 1
      end if
      do lag = 1, order
         design(:, offset + (lag - 1)*variables + 1:offset + lag*variables) = &
            series(order + 1 - lag:n - lag, :)
      end do
      response = series(order + 1:, :)
      cross_x = matmul(transpose(design), design)
      cross_y = matmul(transpose(design), response)
      out%posterior_precision = cross_x + prior%precision
      call invert_matrix(out%posterior_precision, posterior_inverse, status)
      if (status /= 0) then
         out%info = 10 + status
         return
      end if
      weighted_mean = cross_y + matmul(prior%precision, prior%mean)
      out%coefficient_mean = matmul(posterior_inverse, weighted_mean)
      out%residuals = response - matmul(design, out%coefficient_mean)
      difference = out%coefficient_mean - prior%mean
      out%sigma = (prior%scale + matmul(transpose(out%residuals), out%residuals) + &
         matmul(transpose(difference), matmul(prior%precision, difference)))/denominator
      allocate(out%standard_errors(regressors, variables))
      allocate(out%coefficient_covariance(regressors*variables, regressors*variables))
      do equation = 1, variables
         do i = 1, regressors
            out%standard_errors(i, equation) = &
               sqrt(max(0.0_dp, out%sigma(equation, equation)*posterior_inverse(i, i)))
         end do
         do other = 1, variables
            out%coefficient_covariance((equation - 1)*regressors + 1:equation*regressors, &
               (other - 1)*regressors + 1:other*regressors) = &
               out%sigma(equation, other)*posterior_inverse
         end do
      end do
      out%prior = prior
      out%posterior_degrees_of_freedom = prior%degrees_of_freedom + real(rows, dp)
      call build_bvar_model(out)

   contains

      pure subroutine build_bvar_model(fitted)
         !! Convert posterior coefficient means to the shared VAR representation.
         type(mts_bvar_fit_t), intent(inout) :: fitted !! Fitted, updated in place.
         integer :: current_lag, start

         allocate(fitted%model%ar(variables, variables, order))
         allocate(fitted%model%intercept(variables), fitted%model%sigma(variables, variables))
         fitted%model%intercept = 0.0_dp
         start = 0
         if (use_mean) then
            fitted%model%intercept = fitted%coefficient_mean(1, :)
            start = 1
         end if
         do current_lag = 1, order
            fitted%model%ar(:, :, current_lag) = transpose(fitted%coefficient_mean(&
               start + (current_lag - 1)*variables + 1:start + current_lag*variables, :))
         end do
         fitted%model%sigma = fitted%sigma
         fitted%model%residuals = fitted%residuals
         fitted%model%lags = [(current_lag, current_lag=1, order)]
         fitted%model%max_lag = order
         fitted%model%includes_mean = use_mean
      end subroutine build_bvar_model
   end function mts_bvar_fit

   pure function mts_minnesota_prior(series, order, tightness, lag_decay, &
      intercept_variance, random_walk, degrees_of_freedom, include_mean) result(out)
      !! Construct a common-precision Minnesota-style prior for MTS BVAR.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      integer, intent(in) :: order !! Model or polynomial order.
      real(dp), intent(in), optional :: tightness !! Tightness.
      real(dp), intent(in), optional :: lag_decay !! Lag decay.
      real(dp), intent(in), optional :: intercept_variance !! Intercept variance.
      real(dp), intent(in), optional :: degrees_of_freedom !! Degrees of freedom.
      logical, intent(in), optional :: random_walk !! Flag controlling random walk.
      logical, intent(in), optional :: include_mean !! Whether to include a mean term.
      type(mts_bvar_prior_t) :: out
      real(dp), allocatable :: centered(:, :)
      real(dp) :: lambda, decay, constant_variance
      integer :: n, variables, regressors, offset, lag, variable
      logical :: own_unit_root, use_mean

      n = size(series, 1)
      variables = size(series, 2)
      use_mean = .true.
      own_unit_root = .true.
      if (present(include_mean)) use_mean = include_mean
      if (present(random_walk)) own_unit_root = random_walk
      lambda = 0.2_dp
      decay = 1.0_dp
      constant_variance = 1.0e6_dp
      if (present(tightness)) lambda = tightness
      if (present(lag_decay)) decay = lag_decay
      if (present(intercept_variance)) constant_variance = intercept_variance
      if (order < 1 .or. n < 3 .or. variables < 1 .or. lambda <= 0.0_dp .or. &
         decay < 0.0_dp .or. constant_variance <= 0.0_dp) then
         out%info = 1
         return
      end if
      regressors = variables*order + merge(1, 0, use_mean)
      allocate(out%mean(regressors, variables), out%precision(regressors, regressors))
      allocate(out%scale(variables, variables), centered(n - 1, variables))
      out%mean = 0.0_dp
      out%precision = 0.0_dp
      offset = 0
      if (use_mean) then
         out%precision(1, 1) = 1.0_dp/constant_variance
         offset = 1
      end if
      do lag = 1, order
         do variable = 1, variables
            out%precision(offset + (lag - 1)*variables + variable, &
               offset + (lag - 1)*variables + variable) = real(lag, dp)**(2.0_dp*decay)/lambda**2
         end do
      end do
      if (own_unit_root) then
         do variable = 1, variables
            out%mean(offset + variable, variable) = 1.0_dp
         end do
      end if
      centered = series(2:, :) - series(:n - 1, :)
      out%scale = matmul(transpose(centered), centered)/real(n - 1, dp)
      out%degrees_of_freedom = real(variables + 2, dp)
      if (present(degrees_of_freedom)) out%degrees_of_freedom = degrees_of_freedom
      if (out%degrees_of_freedom <= real(variables - 1, dp)) out%info = 2
   end function mts_minnesota_prior

   pure function mts_common_volatility(returns, max_lag, var_order, standardized, &
      arch_lags) result(out)
      !! Extract MTS common-volatility directions from quadratic lag dependence.
      real(dp), intent(in) :: returns(:, :) !! Returns.
      integer, intent(in) :: max_lag !! Maximum lag to consider.
      integer, intent(in), optional :: var_order !! Var order.
      integer, intent(in), optional :: arch_lags(:) !! Arch lags.
      logical, intent(in), optional :: standardized !! Flag controlling standardized.
      type(mts_common_volatility_t) :: out
      type(mts_var_fit_t) :: prewhitened
      real(dp), allocatable :: covariance(:, :), values(:), vectors(:, :), inverse_sqrt(:, :)
      real(dp), allocatable :: cross_covariance(:, :), first_product(:), second_product(:)
      real(dp), allocatable :: diagonal_scale(:)
      real(dp) :: total, norm
      integer :: n, variables, order, lag, first, second, row, column, status, i, j
      logical :: use_standardization

      n = size(returns, 1)
      variables = size(returns, 2)
      order = 1
      use_standardization = .false.
      if (present(var_order)) order = var_order
      if (present(standardized)) use_standardization = standardized
      if (variables < 2 .or. n < 5 .or. max_lag < 1 .or. max_lag >= n .or. order < 0) then
         out%info = 1
         return
      end if
      if (order > 0) then
         prewhitened = mts_var(returns, order)
         if (prewhitened%info /= 0) then
            out%info = 10 + prewhitened%info
            return
         end if
         out%residuals = prewhitened%residuals
      else
         out%residuals = returns - spread(sum(returns, 1)/real(n, dp), 1, n)
      end if
      n = size(out%residuals, 1)
      if (max_lag >= n) then
         out%info = 2
         return
      end if
      covariance = matmul(transpose(out%residuals), out%residuals)/real(n - 1, dp)
      call symmetric_eigen(covariance, values, vectors, status)
      if (status /= 0 .or. minval(values) <= 100.0_dp*epsilon(1.0_dp)*maxval(values)) then
         out%info = 20 + status
         return
      end if
      inverse_sqrt = matmul(vectors, matmul(diagonal_matrix(1.0_dp/sqrt(values)), transpose(vectors)))
      out%whitened = matmul(out%residuals, inverse_sqrt)
      allocate(out%aggregate(variables, variables), cross_covariance(variables, variables))
      out%aggregate = 0.0_dp
      do lag = 1, max_lag
         do first = 1, variables
            do second = first, variables
               second_product = out%whitened(:n - lag, first)*out%whitened(:n - lag, second)
               do row = 1, variables
                  do column = row, variables
                     first_product = out%whitened(lag + 1:, row)*out%whitened(lag + 1:, column)
                     cross_covariance(row, column) = sample_covariance_shared(first_product, second_product)* &
                        (real(n - lag, dp)/real(n, dp))**2
                     cross_covariance(column, row) = cross_covariance(row, column)
                  end do
               end do
               out%aggregate = out%aggregate + matmul(cross_covariance, cross_covariance)
            end do
         end do
      end do
      if (use_standardization) then
         allocate(diagonal_scale(variables))
         do i = 1, variables
            if (out%aggregate(i, i) <= tiny(1.0_dp)) then
               out%info = 3
               return
            end if
            diagonal_scale(i) = 1.0_dp/sqrt(out%aggregate(i, i))
         end do
         out%aggregate = out%aggregate*spread(diagonal_scale, 1, variables)* &
            spread(diagonal_scale, 2, variables)
      else
         out%aggregate = out%aggregate/real(variables*(variables + 1)/2, dp)
      end if
      call symmetric_eigen(out%aggregate, out%eigenvalues, out%eigenvectors, status)
      if (status /= 0) then
         out%info = 30 + status
         return
      end if
      out%eigenvalues = max(0.0_dp, out%eigenvalues)
      total = sum(out%eigenvalues)
      allocate(out%proportions(variables))
      if (total > tiny(1.0_dp)) then
         out%proportions = out%eigenvalues/total
      else
         out%proportions = 0.0_dp
      end if
      out%directions = matmul(inverse_sqrt, out%eigenvectors)
      do j = 1, variables
         norm = norm2(out%directions(:, j))
         if (norm > tiny(1.0_dp)) out%directions(:, j) = out%directions(:, j)/norm
      end do
      out%components = matmul(out%residuals, out%directions)
      if (present(arch_lags)) then
         if (size(arch_lags) < 1 .or. any(arch_lags < 1) .or. any(arch_lags >= n - 2)) then
            out%info = 4
            return
         end if
         out%arch_lags = arch_lags
      else
         out%arch_lags = [10, 20, 30]
         if (any(out%arch_lags >= n - 2)) out%arch_lags = [min(10, max(1, n/4))]
      end if
      allocate(out%arch_statistic(variables, size(out%arch_lags)))
      allocate(out%arch_p_value(variables, size(out%arch_lags)))
      do j = 1, variables
         do i = 1, size(out%arch_lags)
            call arch_f_test(out%components(:, j), out%arch_lags(i), &
               out%arch_statistic(j, i), out%arch_p_value(j, i), status)
            if (status /= 0) then
               out%info = 40 + j
               return
            end if
         end do
      end do
      out%prewhiten_order = order
      out%max_lag = max_lag
      out%standardized = use_standardization

   contains

      pure function diagonal_matrix(diagonal) result(matrix)
         !! Form a diagonal matrix from supplied values.
         real(dp), intent(in) :: diagonal(:) !! Diagonal.
         real(dp) :: matrix(size(diagonal), size(diagonal))
         integer :: index

         matrix = 0.0_dp
         do index = 1, size(diagonal)
            matrix(index, index) = diagonal(index)
         end do
      end function diagonal_matrix

   end function mts_common_volatility

   pure function mts_mch_diagnostic(residuals, covariance, max_lag, robust_probability) result(out)
      !! Diagnose a fitted multivariate conditional covariance path as in MCHdiag.
      real(dp), intent(in) :: residuals(:, :) !! Model residuals.
      real(dp), intent(in) :: covariance(:, :, :) !! Covariance matrix.
      integer, intent(in) :: max_lag !! Maximum lag to consider.
      real(dp), intent(in), optional :: robust_probability !! Robust probability.
      type(mts_mch_diagnostic_t) :: out
      real(dp), allocatable :: values(:), vectors(:, :), inverse_sqrt(:, :)
      real(dp), allocatable :: autocorrelation(:), ranks(:), squared(:, :), robust_squared(:, :)
      real(dp), allocatable :: sorted(:)
      real(dp) :: probability, cutoff, mean_value, denominator, mu, variance_rank
      integer :: n, variables, lag, t, status, retained, index

      n = size(residuals, 1)
      variables = size(residuals, 2)
      probability = 0.95_dp
      if (present(robust_probability)) probability = robust_probability
      if (n < 5 .or. variables < 1 .or. any(shape(covariance) /= [variables, variables, n]) .or. &
         max_lag < 1 .or. max_lag >= n - variables - 1 .or. &
         probability <= 0.0_dp .or. probability >= 1.0_dp) then
         out%info = 1
         return
      end if
      allocate(out%standardized_residuals(n, variables), out%radial_residual(n))
      do t = 1, n
         call symmetric_eigen(covariance(:, :, t), values, vectors, status)
         if (status /= 0 .or. minval(values) <= &
            100.0_dp*epsilon(1.0_dp)*max(1.0_dp, maxval(values))) then
            out%info = 10 + t
            return
         end if
         inverse_sqrt = matmul(vectors, matmul(diagonal_path(1.0_dp/sqrt(values)), transpose(vectors)))
         out%standardized_residuals(t, :) = matmul(residuals(t, :), inverse_sqrt)
         out%radial_residual(t) = sum(out%standardized_residuals(t, :)**2) - real(variables, dp)
      end do
      allocate(autocorrelation(max_lag))
      mean_value = sum(out%radial_residual)/real(n, dp)
      denominator = sum((out%radial_residual - mean_value)**2)
      if (denominator <= tiny(1.0_dp)) then
         out%info = 2
         return
      end if
      do lag = 1, max_lag
         autocorrelation(lag) = sum((out%radial_residual(lag + 1:) - mean_value)* &
            (out%radial_residual(:n - lag) - mean_value))/denominator
         out%radial_q = out%radial_q + autocorrelation(lag)**2/real(n - lag, dp)
      end do
      out%radial_q = out%radial_q*real(n, dp)*real(n + 2, dp)
      out%degrees_of_freedom(1) = max_lag
      out%p_value(1) = chi_square_survival(out%radial_q, max_lag)

      allocate(ranks(n))
      do t = 1, n
         ranks(t) = 1.0_dp + real(count(out%radial_residual < out%radial_residual(t)), dp) + &
            0.5_dp*real(count(out%radial_residual == out%radial_residual(t)) - 1, dp)
      end do
      mean_value = sum(ranks)/real(n, dp)
      denominator = sum((ranks - mean_value)**2)
      do lag = 1, max_lag
         autocorrelation(lag) = sum((ranks(lag + 1:) - mean_value)*(ranks(:n - lag) - mean_value))/denominator
         mu = -real(n - lag, dp)/(real(n, dp)*real(n - 1, dp))
         variance_rank = (5.0_dp*real(n, dp)**4 - real(5*lag + 9, dp)*real(n, dp)**3 + &
            9.0_dp*real(lag - 2, dp)*real(n, dp)**2 + &
            2.0_dp*real(lag*(5*lag + 8), dp)*real(n, dp) + 16.0_dp*real(lag*lag, dp))/ &
            (5.0_dp*real(n - 1, dp)**2*real(n, dp)**2*real(n + 1, dp))
         out%rank_q = out%rank_q + (autocorrelation(lag) - mu)**2/variance_rank
      end do
      out%degrees_of_freedom(2) = max_lag
      out%p_value(2) = chi_square_survival(out%rank_q, max_lag)

      squared = out%standardized_residuals**2
      call multivariate_portmanteau(squared, max_lag, out%multivariate_q, status)
      if (status /= 0) then
         out%info = 20 + status
         return
      end if
      out%degrees_of_freedom(3) = variables*variables*max_lag
      out%p_value(3) = chi_square_survival(out%multivariate_q, out%degrees_of_freedom(3))

      sorted = out%radial_residual
      call sort(sorted)
      index = max(1, min(n, ceiling(probability*real(n, dp))))
      cutoff = sorted(index)
      retained = count(out%radial_residual <= cutoff)
      allocate(robust_squared(retained, variables))
      index = 0
      do t = 1, n
         if (out%radial_residual(t) <= cutoff) then
            index = index + 1
            robust_squared(index, :) = squared(t, :)
         end if
      end do
      if (retained <= max_lag + variables + 1) then
         out%info = 3
         return
      end if
      call multivariate_portmanteau(robust_squared, max_lag, out%robust_q, status)
      if (status /= 0) then
         out%info = 30 + status
         return
      end if
      out%degrees_of_freedom(4) = variables*variables*max_lag
      out%p_value(4) = chi_square_survival(out%robust_q, out%degrees_of_freedom(4))
      out%max_lag = max_lag
      out%robust_observations = retained

   contains

      pure function diagonal_path(diagonal) result(matrix)
         !! Form a diagonal matrix for covariance standardization.
         real(dp), intent(in) :: diagonal(:) !! Diagonal.
         real(dp) :: matrix(size(diagonal), size(diagonal))
         integer :: i

         matrix = 0.0_dp
         do i = 1, size(diagonal)
            matrix(i, i) = diagonal(i)
         end do
      end function diagonal_path

      pure subroutine multivariate_portmanteau(observations, lag_count, statistic, evaluation_status)
         !! Compute the MTS squared-standardized-residual portmanteau statistic.
         real(dp), intent(in) :: observations(:, :) !! Observed time-series values.
         integer, intent(in) :: lag_count !! Number of lag.
         real(dp), intent(out) :: statistic !! Statistic.
         integer, intent(out) :: evaluation_status !! Evaluation status.
         real(dp), allocatable :: covariance_zero(:, :), covariance_inverse(:, :), covariance_lag(:, :)
         real(dp), allocatable :: centered(:, :), product(:, :)
         integer :: observation_count, current_lag, diagonal

         observation_count = size(observations, 1)
         centered = observations - spread(sum(observations, 1)/real(observation_count, dp), &
            1, observation_count)
         covariance_zero = matmul(transpose(centered), centered)/real(observation_count - 1, dp)
         call invert_matrix(covariance_zero, covariance_inverse, evaluation_status)
         if (evaluation_status /= 0) return
         statistic = 0.0_dp
         do current_lag = 1, lag_count
            covariance_lag = matmul(transpose(centered(current_lag + 1:, :)), &
               centered(:observation_count - current_lag, :))/real(observation_count - 1, dp)
            product = matmul(transpose(covariance_lag), matmul(covariance_inverse, &
               matmul(covariance_lag, covariance_inverse)))
            statistic = statistic + real(observation_count, dp)**2* &
               sum([(product(diagonal, diagonal), diagonal=1, size(observations, 2))])/ &
               real(observation_count - current_lag, dp)
         end do
      end subroutine multivariate_portmanteau

   end function mts_mch_diagnostic

   pure function mts_bekk_fit(returns, include_mean, initial, estimated, &
      max_iterations, tolerance) result(out)
      !! Estimate the Gaussian bivariate or trivariate MTS BEKK(1,1) model.
      real(dp), intent(in) :: returns(:, :) !! Returns.
      logical, intent(in), optional :: include_mean !! Whether to include a mean term.
      logical, intent(in), optional :: estimated(:) !! Flag controlling estimated.
      real(dp), intent(in), optional :: initial(:) !! Initial value.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      type(mts_bekk_fit_t) :: out
      type(optimization_result_t) :: optimization
      real(dp), allocatable :: parameters(:), free_parameters(:), covariance(:, :)
      real(dp), allocatable :: lower(:, :), hessian(:, :), inverse(:, :)
      integer, allocatable :: free_index(:)
      integer :: n, dimension, triangular_count, parameter_count, free_count
      integer :: offset, row, column, i, status, limit
      real(dp) :: gradient_tolerance
      logical :: use_mean

      n = size(returns, 1)
      dimension = size(returns, 2)
      use_mean = .true.
      if (present(include_mean)) use_mean = include_mean
      triangular_count = dimension*(dimension + 1)/2
      parameter_count = merge(dimension, 0, use_mean) + triangular_count + 2*dimension*dimension
      if ((dimension /= 2 .and. dimension /= 3) .or. n < 20) then
         out%info = 1
         return
      end if
      allocate(parameters(parameter_count), out%estimated(parameter_count))
      parameters = 0.0_dp
      offset = 0
      if (use_mean) then
         parameters(:dimension) = sum(returns, 1)/real(n, dp)
         offset = dimension
      end if
      covariance = matmul(transpose(returns - spread(sum(returns, 1)/real(n, dp), 1, n)), &
         returns - spread(sum(returns, 1)/real(n, dp), 1, n))/real(n - 1, dp)
      call cholesky_lower(covariance, lower, status)
      if (status /= 0) then
         out%info = 10 + status
         return
      end if
      do column = 1, dimension
         do row = column, dimension
            offset = offset + 1
            parameters(offset) = lower(row, column)
         end do
      end do
      do column = 1, dimension
         do row = 1, dimension
            offset = offset + 1
            parameters(offset) = merge(0.1_dp, 0.02_dp, row == column)
         end do
      end do
      do column = 1, dimension
         do row = 1, dimension
            offset = offset + 1
            parameters(offset) = merge(merge(0.9_dp, 0.8_dp, dimension == 2), &
               0.02_dp, row == column)
         end do
      end do
      if (present(initial)) then
         if (size(initial) /= parameter_count) then
            out%info = 2
            return
         end if
         parameters = initial
      end if
      out%estimated = .true.
      if (present(estimated)) then
         if (size(estimated) /= parameter_count) then
            out%info = 2
            return
         end if
         out%estimated = estimated
      end if
      free_count = count(out%estimated)
      allocate(free_index(free_count), free_parameters(free_count))
      free_index = pack([(i, i=1, parameter_count)], out%estimated)
      free_parameters = pack(parameters, out%estimated)
      limit = 250
      gradient_tolerance = 1.0e-5_dp
      if (present(max_iterations)) limit = max_iterations
      if (present(tolerance)) gradient_tolerance = tolerance
      if (free_count > 0) then
         optimization = bfgs_minimize_fd(objective, free_parameters, limit, gradient_tolerance)
         parameters(free_index) = optimization%parameters
         out%iterations = optimization%iterations
         out%converged = optimization%converged
         if (optimization%info /= 0) out%info = 100 + optimization%info
      else
         out%converged = .true.
      end if
      call evaluate(parameters, out%mean, out%constant, out%arch, out%garch, &
         out%residuals, out%covariance, out%standardized_residuals, out%log_likelihood, status)
      if (status /= 0) then
         out%info = 20 + status
         return
      end if
      out%coefficients = parameters
      out%persistence = bekk_persistence(out%arch, out%garch)
      out%aic = -2.0_dp*out%log_likelihood + 2.0_dp*real(free_count, dp)
      out%bic = -2.0_dp*out%log_likelihood + log(real(n, dp))*real(free_count, dp)
      out%includes_mean = use_mean
      allocate(out%parameter_covariance(parameter_count, parameter_count))
      allocate(out%standard_errors(parameter_count))
      out%parameter_covariance = 0.0_dp
      out%standard_errors = 0.0_dp
      if (free_count > 0) then
         free_parameters = pack(parameters, out%estimated)
         hessian = finite_difference_hessian(objective, free_parameters)
         call invert_matrix(hessian, inverse, status)
         if (status == 0) then
            do i = 1, free_count
               out%parameter_covariance(free_index, free_index(i)) = inverse(:, i)
               out%standard_errors(free_index(i)) = sqrt(max(0.0_dp, inverse(i, i)))
            end do
         else if (out%info == 0) then
            out%info = 200 + status
         end if
      end if

   contains

      pure function objective(free_values) result(value)
         !! Return the BEKK Gaussian negative log likelihood.
         real(dp), intent(in) :: free_values(:) !! Free values.
         real(dp) :: value, log_likelihood
         real(dp) :: full_values(parameter_count)
         real(dp), allocatable :: mean(:), constant(:, :), arch(:, :), garch(:, :)
         real(dp), allocatable :: residuals(:, :), covariance(:, :, :), standardized(:, :)
         integer :: evaluation_status

         full_values = parameters
         full_values(free_index) = free_values
         call evaluate(full_values, mean, constant, arch, garch, residuals, covariance, &
            standardized, log_likelihood, evaluation_status)
         if (evaluation_status == 0 .and. ieee_is_finite(log_likelihood)) then
            value = -log_likelihood
         else
            value = 1.0e30_dp + dot_product(free_values, free_values)
         end if
      end function objective

      pure subroutine evaluate(values, mean, constant, arch, garch, residuals, covariance, &
         standardized, log_likelihood, evaluation_status)
         !! Unpack BEKK parameters and evaluate its covariance recursion.
         real(dp), intent(in) :: values(:) !! Input values.
         real(dp), allocatable, intent(out) :: mean(:) !! Mean value or vector.
         real(dp), allocatable, intent(out) :: constant(:, :) !! Constant.
         real(dp), allocatable, intent(out) :: arch(:, :) !! Arch.
         real(dp), allocatable, intent(out) :: garch(:, :) !! Garch.
         real(dp), allocatable, intent(out) :: residuals(:, :) !! Model residuals.
         real(dp), allocatable, intent(out) :: covariance(:, :, :) !! Covariance matrix.
         real(dp), allocatable, intent(out) :: standardized(:, :) !! Standardized.
         real(dp), intent(out) :: log_likelihood !! Log-likelihood value.
         integer, intent(out) :: evaluation_status !! Evaluation status.
         real(dp), allocatable :: covariance_inverse(:, :), cholesky(:, :), shock(:)
         real(dp) :: logdet, quadratic
         integer :: value_offset, current_row, current_column, t

         allocate(mean(dimension), constant(dimension, dimension))
         allocate(arch(dimension, dimension), garch(dimension, dimension))
         mean = 0.0_dp
         constant = 0.0_dp
         value_offset = 0
         if (use_mean) then
            mean = values(:dimension)
            value_offset = dimension
         end if
         do current_column = 1, dimension
            do current_row = current_column, dimension
               value_offset = value_offset + 1
               constant(current_row, current_column) = values(value_offset)
            end do
         end do
         arch = reshape(values(value_offset + 1:value_offset + dimension*dimension), &
            [dimension, dimension])
         value_offset = value_offset + dimension*dimension
         garch = reshape(values(value_offset + 1:value_offset + dimension*dimension), &
            [dimension, dimension])
         do current_row = 1, dimension
            if (constant(current_row, current_row) <= sqrt(epsilon(1.0_dp))) then
               evaluation_status = 1
               log_likelihood = -huge(1.0_dp)
               return
            end if
         end do
         if (bekk_persistence(arch, garch) >= 0.999_dp) then
            evaluation_status = 2
            log_likelihood = -huge(1.0_dp)
            return
         end if
         residuals = returns - spread(mean, 1, n)
         allocate(covariance(dimension, dimension, n), standardized(n, dimension))
         covariance(:, :, 1) = matmul(transpose(residuals), residuals)/real(n - 1, dp)
         standardized = 0.0_dp
         log_likelihood = 0.0_dp
         allocate(covariance_inverse(dimension, dimension), shock(dimension))
         do t = 2, n
            shock = matmul(arch, residuals(t - 1, :))
            covariance(:, :, t) = matmul(constant, transpose(constant)) + &
               spread(shock, 2, dimension)*spread(shock, 1, dimension) + &
               matmul(garch, matmul(covariance(:, :, t - 1), transpose(garch)))
            covariance(:, :, t) = 0.5_dp*(covariance(:, :, t) + transpose(covariance(:, :, t)))
            call inverse_logdet(covariance(:, :, t), covariance_inverse, logdet, &
               evaluation_status, 100.0_dp*epsilon(1.0_dp))
            if (evaluation_status /= 0) then
               log_likelihood = -huge(1.0_dp)
               return
            end if
            quadratic = dot_product(residuals(t, :), matmul(covariance_inverse, residuals(t, :)))
            log_likelihood = log_likelihood - 0.5_dp*(real(dimension, dp)*log(2.0_dp*acos(-1.0_dp)) + &
               logdet + quadratic)
         end do
         do t = 1, n
            call cholesky_lower(covariance(:, :, t), cholesky, evaluation_status)
            if (evaluation_status /= 0) return
            call invert_matrix(cholesky, covariance_inverse, evaluation_status)
            if (evaluation_status /= 0) return
            standardized(t, :) = matmul(covariance_inverse, residuals(t, :))
         end do
      end subroutine evaluate
   end function mts_bekk_fit

   pure function mts_bekk_forecast(model, horizon) result(covariance_forecast)
      !! Forecast BEKK conditional covariance matrices with zero future shocks.
      type(mts_bekk_fit_t), intent(in) :: model !! Model specification.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      real(dp), allocatable :: covariance_forecast(:, :, :)
      real(dp) :: shock(size(model%mean))
      integer :: step, dimension

      if (.not. allocated(model%mean) .or. .not. allocated(model%covariance) .or. &
         .not. allocated(model%residuals) .or. horizon < 1) then
         allocate(covariance_forecast(0, 0, 0))
         return
      end if
      dimension = size(model%mean)
      allocate(covariance_forecast(dimension, dimension, horizon))
      shock = matmul(model%arch, model%residuals(size(model%residuals, 1), :))
      covariance_forecast(:, :, 1) = matmul(model%constant, transpose(model%constant)) + &
         spread(shock, 2, dimension)*spread(shock, 1, dimension) + matmul(model%garch, &
         matmul(model%covariance(:, :, size(model%covariance, 3)), transpose(model%garch)))
      do step = 2, horizon
         covariance_forecast(:, :, step) = matmul(model%constant, transpose(model%constant)) + &
            matmul(model%arch, matmul(covariance_forecast(:, :, step - 1), transpose(model%arch))) + &
            matmul(model%garch, matmul(covariance_forecast(:, :, step - 1), transpose(model%garch)))
      end do
   end function mts_bekk_forecast

   pure real(dp) function bekk_persistence(arch, garch) result(radius)
      !! Approximate the dominant modulus of the BEKK second-moment operator.
      real(dp), intent(in) :: arch(:, :) !! Arch.
      real(dp), intent(in) :: garch(:, :) !! Garch.
      real(dp), allocatable :: operator(:, :), vector(:), next_vector(:)
      real(dp) :: scale
      integer :: dimension, state_size, i, j, row, column, iteration

      dimension = size(arch, 1)
      state_size = dimension*dimension
      allocate(operator(state_size, state_size), vector(state_size), next_vector(state_size))
      do j = 1, dimension
         do i = 1, dimension
            row = i + (j - 1)*dimension
            do column = 1, dimension
               operator(row, (column - 1)*dimension + 1:column*dimension) = &
                  arch(i, :)*arch(j, column) + garch(i, :)*garch(j, column)
            end do
         end do
      end do
      vector = 1.0_dp/sqrt(real(state_size, dp))
      radius = 0.0_dp
      do iteration = 1, 200
         next_vector = matmul(operator, vector)
         scale = norm2(next_vector)
         if (scale <= tiny(1.0_dp)) then
            radius = 0.0_dp
            return
         end if
         vector = next_vector/scale
         radius = scale
      end do
   end function bekk_persistence

   pure function mts_dcc_fit(standardized_residuals, marginal_variance, initial, &
      max_iterations, tolerance) result(out)
      !! Estimate Gaussian Engle DCC(1,1) correlations by stage-two likelihood.
      real(dp), intent(in) :: standardized_residuals(:, :) !! Standardized residuals.
      real(dp), intent(in), optional :: marginal_variance(:, :) !! Marginal variance.
      real(dp), intent(in), optional :: initial(:) !! Initial value.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      type(mts_dcc_fit_t) :: out
      type(optimization_result_t) :: optimization
      real(dp) :: raw_initial(2), physical_initial(2), gradient_tolerance
      real(dp), allocatable :: raw_covariance(:, :)
      real(dp) :: jacobian(2, 2), hessian(2, 2)
      integer :: n, dimension, limit, status, i

      n = size(standardized_residuals, 1)
      dimension = size(standardized_residuals, 2)
      if (n < 20 .or. dimension < 2) then
         out%info = 1
         return
      end if
      allocate(out%marginal_variance(n, dimension))
      out%marginal_variance = 1.0_dp
      if (present(marginal_variance)) then
         if (any(shape(marginal_variance) /= [n, dimension]) .or. &
            any(marginal_variance <= 0.0_dp)) then
            out%info = 2
            return
         end if
         out%marginal_variance = marginal_variance
      end if
      physical_initial = [0.02_dp, 0.95_dp]
      if (present(initial)) then
         if (size(initial) /= 2 .or. any(initial <= 0.0_dp) .or. sum(initial) >= 1.0_dp) then
            out%info = 3
            return
         end if
         physical_initial = initial
      end if
      raw_initial = dcc_inverse_transform(physical_initial)
      limit = 200
      gradient_tolerance = 1.0e-6_dp
      if (present(max_iterations)) limit = max_iterations
      if (present(tolerance)) gradient_tolerance = tolerance
      optimization = bfgs_minimize_fd(objective, raw_initial, limit, gradient_tolerance)
      call dcc_transform(optimization%parameters, out%arch, out%garch)
      out%iterations = optimization%iterations
      out%converged = optimization%converged
      if (optimization%info /= 0) out%info = 100 + optimization%info
      call evaluate(optimization%parameters, out%unconditional, out%q, out%correlation, &
         out%covariance, out%residuals, out%log_likelihood, status)
      if (status /= 0) then
         out%info = 10 + status
         return
      end if
      out%standardized_residuals = standardized_residuals
      out%aic = -2.0_dp*out%log_likelihood + 4.0_dp
      out%bic = -2.0_dp*out%log_likelihood + 2.0_dp*log(real(n, dp))
      hessian = finite_difference_hessian(objective, optimization%parameters)
      call invert_matrix(hessian, raw_covariance, status)
      allocate(out%parameter_covariance(2, 2), out%standard_errors(2))
      out%parameter_covariance = 0.0_dp
      out%standard_errors = 0.0_dp
      if (status == 0) then
         jacobian = reshape([out%arch*(1.0_dp - out%arch), -out%arch*out%garch, &
            -out%arch*out%garch, out%garch*(1.0_dp - out%garch)], [2, 2])
         out%parameter_covariance = matmul(jacobian, matmul(raw_covariance, transpose(jacobian)))
         do i = 1, 2
            out%standard_errors(i) = sqrt(max(0.0_dp, out%parameter_covariance(i, i)))
         end do
      else if (out%info == 0) then
         out%info = 200 + status
      end if

   contains

      pure function objective(raw_parameters) result(value)
         !! Return the Gaussian DCC stage-two negative log likelihood.
         real(dp), intent(in) :: raw_parameters(:) !! Raw parameters.
         real(dp) :: value, log_likelihood
         real(dp), allocatable :: unconditional(:, :), q(:, :, :), correlation(:, :, :)
         real(dp), allocatable :: covariance(:, :, :), residuals(:, :)
         integer :: evaluation_status

         call evaluate(raw_parameters, unconditional, q, correlation, covariance, residuals, &
            log_likelihood, evaluation_status)
         if (evaluation_status == 0 .and. ieee_is_finite(log_likelihood)) then
            value = -log_likelihood
         else
            value = 1.0e30_dp + dot_product(raw_parameters, raw_parameters)
         end if
      end function objective

      pure subroutine evaluate(raw_parameters, unconditional, q, correlation, covariance, &
         residuals, log_likelihood, evaluation_status)
         !! Evaluate DCC Q, correlation, covariance, and likelihood paths.
         real(dp), intent(in) :: raw_parameters(:) !! Raw parameters.
         real(dp), allocatable, intent(out) :: unconditional(:, :) !! Unconditional.
         real(dp), allocatable, intent(out) :: q(:, :, :) !! Model order, dimension, or parameter.
         real(dp), allocatable, intent(out) :: correlation(:, :, :) !! Correlation.
         real(dp), allocatable, intent(out) :: covariance(:, :, :) !! Covariance matrix.
         real(dp), allocatable, intent(out) :: residuals(:, :) !! Model residuals.
         real(dp), intent(out) :: log_likelihood !! Log-likelihood value.
         integer, intent(out) :: evaluation_status !! Evaluation status.
         real(dp), allocatable :: correlation_inverse(:, :)
         real(dp) :: arch_value, garch_value, logdet, quadratic, independent
         real(dp) :: scale(dimension), shock(dimension)
         integer :: t, row, column

         call dcc_transform(raw_parameters, arch_value, garch_value)
         unconditional = matmul(transpose(standardized_residuals), standardized_residuals)/real(n, dp)
         allocate(q(dimension, dimension, n), correlation(dimension, dimension, n))
         allocate(covariance(dimension, dimension, n), residuals(n, dimension))
         q(:, :, 1) = unconditional
         residuals = standardized_residuals*sqrt(out%marginal_variance)
         log_likelihood = 0.0_dp
         allocate(correlation_inverse(dimension, dimension))
         do t = 1, n
            do row = 1, dimension
               scale(row) = sqrt(max(q(row, row, t), tiny(1.0_dp)))
            end do
            do column = 1, dimension
               do row = 1, dimension
                  correlation(row, column, t) = q(row, column, t)/(scale(row)*scale(column))
               end do
            end do
            do row = 1, dimension
               correlation(row, row, t) = 1.0_dp
            end do
            call inverse_logdet(correlation(:, :, t), correlation_inverse, logdet, &
               evaluation_status, 100.0_dp*epsilon(1.0_dp))
            if (evaluation_status /= 0) then
               log_likelihood = -huge(1.0_dp)
               return
            end if
            quadratic = dot_product(standardized_residuals(t, :), &
               matmul(correlation_inverse, standardized_residuals(t, :)))
            independent = dot_product(standardized_residuals(t, :), standardized_residuals(t, :))
            log_likelihood = log_likelihood - 0.5_dp*(logdet + quadratic - independent)
            scale = sqrt(out%marginal_variance(t, :))
            covariance(:, :, t) = correlation(:, :, t)* &
               spread(scale, 1, dimension)*spread(scale, 2, dimension)
            if (t < n) then
               shock = standardized_residuals(t, :)
               q(:, :, t + 1) = (1.0_dp - arch_value - garch_value)*unconditional + &
                  spread(shock, 2, dimension)*spread(shock, 1, dimension)*arch_value + &
                  garch_value*q(:, :, t)
            end if
         end do
      end subroutine evaluate
   end function mts_dcc_fit

   pure function mts_dcc_forecast(model, future_variance) result(covariance_forecast)
      !! Forecast DCC covariance matrices for supplied marginal variance forecasts.
      type(mts_dcc_fit_t), intent(in) :: model !! Model specification.
      real(dp), intent(in) :: future_variance(:, :) !! Future variance.
      real(dp), allocatable :: covariance_forecast(:, :, :)
      real(dp), allocatable :: q_forecast(:, :, :)
      real(dp) :: scale(size(future_variance, 2)), shock(size(future_variance, 2))
      integer :: horizon, dimension, step, row, column

      horizon = size(future_variance, 1)
      dimension = size(future_variance, 2)
      if (.not. allocated(model%q) .or. horizon < 1 .or. &
         dimension /= size(model%unconditional, 1) .or. any(future_variance <= 0.0_dp)) then
         allocate(covariance_forecast(0, 0, 0))
         return
      end if
      allocate(covariance_forecast(dimension, dimension, horizon))
      allocate(q_forecast(dimension, dimension, horizon))
      shock = model%standardized_residuals(size(model%standardized_residuals, 1), :)
      q_forecast(:, :, 1) = (1.0_dp - model%arch - model%garch)*model%unconditional + &
         model%arch*spread(shock, 2, dimension)*spread(shock, 1, dimension) + &
         model%garch*model%q(:, :, size(model%q, 3))
      do step = 1, horizon
         do row = 1, dimension
            scale(row) = sqrt(q_forecast(row, row, step))
         end do
         do column = 1, dimension
            do row = 1, dimension
               covariance_forecast(row, column, step) = q_forecast(row, column, step)/ &
                  (scale(row)*scale(column))*sqrt(future_variance(step, row)* &
                  future_variance(step, column))
            end do
         end do
         if (step < horizon) q_forecast(:, :, step + 1) = &
            (1.0_dp - model%arch - model%garch)*model%unconditional + &
            (model%arch + model%garch)*q_forecast(:, :, step)
      end do
   end function mts_dcc_forecast

   pure subroutine dcc_transform(raw_parameters, arch, garch)
      !! Map unconstrained DCC parameters to positive coefficients summing below one.
      real(dp), intent(in) :: raw_parameters(:) !! Raw parameters.
      real(dp), intent(out) :: arch !! Arch.
      real(dp), intent(out) :: garch !! Garch.
      real(dp) :: maximum, weights(3), denominator

      maximum = max(0.0_dp, maxval(raw_parameters))
      weights = [exp(-maximum), exp(raw_parameters(1) - maximum), &
         exp(raw_parameters(2) - maximum)]
      denominator = sum(weights)
      arch = weights(2)/denominator
      garch = weights(3)/denominator
   end subroutine dcc_transform

   pure function dcc_inverse_transform(parameters) result(raw_parameters)
      !! Map physical DCC coefficients to unconstrained optimizer parameters.
      real(dp), intent(in) :: parameters(:) !! Model parameter values.
      real(dp) :: raw_parameters(2), remainder

      remainder = max(1.0_dp - sum(parameters), tiny(1.0_dp))
      raw_parameters = log(max(parameters, tiny(1.0_dp))/remainder)
   end function dcc_inverse_transform

   pure function mts_adcc_fit(standardized_residuals, marginal_variance, initial, &
      max_iterations, tolerance) result(out)
      !! Estimate Gaussian asymmetric DCC correlations.
      real(dp), intent(in) :: standardized_residuals(:, :) !! Standardized residuals.
      real(dp), intent(in), optional :: marginal_variance(:, :) !! Marginal variance.
      real(dp), intent(in), optional :: initial(:) !! Initial value.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      type(mts_adcc_fit_t) :: out

      out = fit_adcc_core(standardized_residuals, .false., marginal_variance, initial, &
         max_iterations, tolerance)
   end function mts_adcc_fit

   pure function mts_adcc_t_fit(standardized_residuals, marginal_variance, initial, &
      max_iterations, tolerance) result(out)
      !! Estimate standardized Student-t asymmetric DCC correlations.
      real(dp), intent(in) :: standardized_residuals(:, :) !! Standardized residuals.
      real(dp), intent(in), optional :: marginal_variance(:, :) !! Marginal variance.
      real(dp), intent(in), optional :: initial(:) !! Initial value.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      type(mts_adcc_fit_t) :: out

      out = fit_adcc_core(standardized_residuals, .true., marginal_variance, initial, &
         max_iterations, tolerance)
   end function mts_adcc_t_fit

   pure function fit_adcc_core(standardized_residuals, use_student_t, marginal_variance, &
      initial, max_iterations, tolerance) result(out)
      !! Fit Gaussian or Student-t ADCC through a shared constrained recursion.
      real(dp), intent(in) :: standardized_residuals(:, :) !! Standardized residuals.
      logical, intent(in) :: use_student_t !! Whether to use the student t.
      real(dp), intent(in), optional :: marginal_variance(:, :) !! Marginal variance.
      real(dp), intent(in), optional :: initial(:) !! Initial value.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      type(mts_adcc_fit_t) :: out
      type(optimization_result_t) :: optimization
      real(dp), allocatable :: raw_initial(:), raw_covariance(:, :), hessian(:, :), jacobian(:, :)
      real(dp), allocatable :: physical_initial(:)
      real(dp) :: gradient_tolerance
      integer :: n, dimension, parameter_count, limit, status, i

      n = size(standardized_residuals, 1)
      dimension = size(standardized_residuals, 2)
      parameter_count = merge(4, 3, use_student_t)
      if (n < 20 .or. dimension < 2) then
         out%info = 1
         return
      end if
      allocate(out%marginal_variance(n, dimension))
      out%marginal_variance = 1.0_dp
      if (present(marginal_variance)) then
         if (any(shape(marginal_variance) /= [n, dimension]) .or. &
            any(marginal_variance <= 0.0_dp)) then
            out%info = 2
            return
         end if
         out%marginal_variance = marginal_variance
      end if
      allocate(physical_initial(parameter_count), raw_initial(parameter_count))
      physical_initial(:3) = [0.02_dp, 0.9_dp, 0.02_dp]
      if (use_student_t) physical_initial(4) = 8.0_dp
      if (present(initial)) then
         if (size(initial) /= parameter_count .or. any(initial(:3) <= 0.0_dp) .or. &
            sum(initial(:3)) >= 1.0_dp) then
            out%info = 3
            return
         end if
         if (use_student_t) then
            if (initial(4) <= 2.0_dp) then
               out%info = 3
               return
            end if
         end if
         physical_initial = initial
      end if
      raw_initial(:3) = adcc_inverse_transform(physical_initial(:3))
      if (use_student_t) raw_initial(4) = log(physical_initial(4) - 2.0_dp)
      limit = 200
      gradient_tolerance = 1.0e-6_dp
      if (present(max_iterations)) limit = max_iterations
      if (present(tolerance)) gradient_tolerance = tolerance
      optimization = bfgs_minimize_fd(objective, raw_initial, limit, gradient_tolerance)
      call adcc_transform(optimization%parameters, out%arch, out%garch, out%asymmetry)
      if (use_student_t) out%degrees_of_freedom = 2.0_dp + exp(min(50.0_dp, optimization%parameters(4)))
      out%iterations = optimization%iterations
      out%converged = optimization%converged
      if (optimization%info /= 0) out%info = 100 + optimization%info
      call evaluate(optimization%parameters, out%unconditional, out%negative_unconditional, &
         out%q, out%correlation, out%covariance, out%residuals, out%log_likelihood, status)
      if (status /= 0) then
         out%info = 10 + status
         return
      end if
      out%standardized_residuals = standardized_residuals
      out%student_t = use_student_t
      out%aic = -2.0_dp*out%log_likelihood + 2.0_dp*real(parameter_count, dp)
      out%bic = -2.0_dp*out%log_likelihood + log(real(n, dp))*real(parameter_count, dp)
      hessian = finite_difference_hessian(objective, optimization%parameters)
      call invert_matrix(hessian, raw_covariance, status)
      allocate(out%parameter_covariance(parameter_count, parameter_count))
      allocate(out%standard_errors(parameter_count), jacobian(parameter_count, parameter_count))
      out%parameter_covariance = 0.0_dp
      out%standard_errors = 0.0_dp
      if (status == 0) then
         jacobian = 0.0_dp
         jacobian(1, 1) = out%arch*(1.0_dp - out%arch)
         jacobian(1, 2) = -out%arch*out%garch
         jacobian(1, 3) = -out%arch*out%asymmetry
         jacobian(2, 1) = -out%garch*out%arch
         jacobian(2, 2) = out%garch*(1.0_dp - out%garch)
         jacobian(2, 3) = -out%garch*out%asymmetry
         jacobian(3, 1) = -out%asymmetry*out%arch
         jacobian(3, 2) = -out%asymmetry*out%garch
         jacobian(3, 3) = out%asymmetry*(1.0_dp - out%asymmetry)
         if (use_student_t) jacobian(4, 4) = out%degrees_of_freedom - 2.0_dp
         out%parameter_covariance = matmul(jacobian, matmul(raw_covariance, transpose(jacobian)))
         do i = 1, parameter_count
            out%standard_errors(i) = sqrt(max(0.0_dp, out%parameter_covariance(i, i)))
         end do
      else if (out%info == 0) then
         out%info = 200 + status
      end if

   contains

      pure function objective(raw_parameters) result(value)
         !! Return the selected ADCC negative log likelihood.
         real(dp), intent(in) :: raw_parameters(:) !! Raw parameters.
         real(dp) :: value, log_likelihood
         real(dp), allocatable :: unconditional(:, :), negative_unconditional(:, :)
         real(dp), allocatable :: q(:, :, :), correlation(:, :, :), covariance(:, :, :), residuals(:, :)
         integer :: evaluation_status

         call evaluate(raw_parameters, unconditional, negative_unconditional, q, correlation, &
            covariance, residuals, log_likelihood, evaluation_status)
         if (evaluation_status == 0 .and. ieee_is_finite(log_likelihood)) then
            value = -log_likelihood
         else
            value = 1.0e30_dp + dot_product(raw_parameters, raw_parameters)
         end if
      end function objective

      pure subroutine evaluate(raw_parameters, unconditional, negative_unconditional, q, &
         correlation, covariance, residuals, log_likelihood, evaluation_status)
         !! Evaluate Gaussian or Student-t ADCC paths and likelihood.
         real(dp), intent(in) :: raw_parameters(:) !! Raw parameters.
         real(dp), allocatable, intent(out) :: unconditional(:, :) !! Unconditional.
         real(dp), allocatable, intent(out) :: negative_unconditional(:, :) !! Negative unconditional.
         real(dp), allocatable, intent(out) :: q(:, :, :) !! Model order, dimension, or parameter.
         real(dp), allocatable, intent(out) :: correlation(:, :, :) !! Correlation.
         real(dp), allocatable, intent(out) :: covariance(:, :, :) !! Covariance matrix.
         real(dp), allocatable, intent(out) :: residuals(:, :) !! Model residuals.
         real(dp), intent(out) :: log_likelihood !! Log-likelihood value.
         integer, intent(out) :: evaluation_status !! Evaluation status.
         real(dp), allocatable :: correlation_inverse(:, :), negative(:, :)
         real(dp) :: arch_value, garch_value, asymmetry_value, degrees, constant
         real(dp) :: logdet, quadratic, independent, scale(dimension), shock(dimension), negative_shock(dimension)
         integer :: t, row, column

         call adcc_transform(raw_parameters, arch_value, garch_value, asymmetry_value)
         degrees = 0.0_dp
         if (use_student_t) degrees = 2.0_dp + exp(min(50.0_dp, raw_parameters(4)))
         unconditional = matmul(transpose(standardized_residuals), standardized_residuals)/real(n, dp)
         negative = min(standardized_residuals, 0.0_dp)
         negative_unconditional = matmul(transpose(negative), negative)/real(n, dp)
         allocate(q(dimension, dimension, n), correlation(dimension, dimension, n))
         allocate(covariance(dimension, dimension, n), residuals(n, dimension))
         q(:, :, 1) = unconditional
         residuals = standardized_residuals*sqrt(out%marginal_variance)
         log_likelihood = 0.0_dp
         if (use_student_t) constant = log_gamma(0.5_dp*degrees) - &
            log_gamma(0.5_dp*(degrees + real(dimension, dp))) + &
            0.5_dp*real(dimension, dp)*log(acos(-1.0_dp)*(degrees - 2.0_dp))
         allocate(correlation_inverse(dimension, dimension))
         do t = 1, n
            do row = 1, dimension
               scale(row) = sqrt(max(q(row, row, t), tiny(1.0_dp)))
            end do
            do column = 1, dimension
               do row = 1, dimension
                  correlation(row, column, t) = q(row, column, t)/(scale(row)*scale(column))
               end do
            end do
            do row = 1, dimension
               correlation(row, row, t) = 1.0_dp
            end do
            call inverse_logdet(correlation(:, :, t), correlation_inverse, logdet, &
               evaluation_status, 100.0_dp*epsilon(1.0_dp))
            if (evaluation_status /= 0) then
               log_likelihood = -huge(1.0_dp)
               return
            end if
            quadratic = dot_product(standardized_residuals(t, :), &
               matmul(correlation_inverse, standardized_residuals(t, :)))
            if (use_student_t) then
               log_likelihood = log_likelihood - constant - 0.5_dp*logdet - &
                  0.5_dp*(degrees + real(dimension, dp))*log(1.0_dp + quadratic/(degrees - 2.0_dp))
            else
               independent = dot_product(standardized_residuals(t, :), standardized_residuals(t, :))
               log_likelihood = log_likelihood - 0.5_dp*(logdet + quadratic - independent)
            end if
            scale = sqrt(out%marginal_variance(t, :))
            covariance(:, :, t) = correlation(:, :, t)* &
               spread(scale, 1, dimension)*spread(scale, 2, dimension)
            if (t < n) then
               shock = standardized_residuals(t, :)
               negative_shock = min(shock, 0.0_dp)
               q(:, :, t + 1) = (1.0_dp - arch_value - garch_value)*unconditional - &
                  asymmetry_value*negative_unconditional + &
                  arch_value*spread(shock, 2, dimension)*spread(shock, 1, dimension) + &
                  garch_value*q(:, :, t) + asymmetry_value* &
                  spread(negative_shock, 2, dimension)*spread(negative_shock, 1, dimension)
            end if
         end do
      end subroutine evaluate
   end function fit_adcc_core

   pure function mts_adcc_forecast(model, future_variance) result(covariance_forecast)
      !! Forecast Gaussian or Student-t ADCC covariance matrices.
      type(mts_adcc_fit_t), intent(in) :: model !! Model specification.
      real(dp), intent(in) :: future_variance(:, :) !! Future variance.
      real(dp), allocatable :: covariance_forecast(:, :, :), q_forecast(:, :, :)
      real(dp) :: scale(size(future_variance, 2)), shock(size(future_variance, 2))
      real(dp) :: negative_shock(size(future_variance, 2))
      integer :: horizon, dimension, step, row, column

      horizon = size(future_variance, 1)
      dimension = size(future_variance, 2)
      if (.not. allocated(model%q) .or. horizon < 1 .or. &
         dimension /= size(model%unconditional, 1) .or. any(future_variance <= 0.0_dp)) then
         allocate(covariance_forecast(0, 0, 0))
         return
      end if
      allocate(covariance_forecast(dimension, dimension, horizon))
      allocate(q_forecast(dimension, dimension, horizon))
      shock = model%standardized_residuals(size(model%standardized_residuals, 1), :)
      negative_shock = min(shock, 0.0_dp)
      q_forecast(:, :, 1) = (1.0_dp - model%arch - model%garch)*model%unconditional - &
         model%asymmetry*model%negative_unconditional + model%arch* &
         spread(shock, 2, dimension)*spread(shock, 1, dimension) + &
         model%garch*model%q(:, :, size(model%q, 3)) + model%asymmetry* &
         spread(negative_shock, 2, dimension)*spread(negative_shock, 1, dimension)
      do step = 1, horizon
         do row = 1, dimension
            scale(row) = sqrt(q_forecast(row, row, step))
         end do
         do column = 1, dimension
            do row = 1, dimension
               covariance_forecast(row, column, step) = q_forecast(row, column, step)/ &
                  (scale(row)*scale(column))*sqrt(future_variance(step, row)* &
                  future_variance(step, column))
            end do
         end do
         if (step < horizon) q_forecast(:, :, step + 1) = &
            (1.0_dp - model%arch - model%garch)*model%unconditional + &
            (model%arch + model%garch)*q_forecast(:, :, step)
      end do
   end function mts_adcc_forecast

   pure subroutine adcc_transform(raw_parameters, arch, garch, asymmetry)
      !! Map unconstrained ADCC parameters to positive coefficients summing below one.
      real(dp), intent(in) :: raw_parameters(:) !! Raw parameters.
      real(dp), intent(out) :: arch !! Arch.
      real(dp), intent(out) :: garch !! Garch.
      real(dp), intent(out) :: asymmetry !! Asymmetry.
      real(dp) :: maximum, weights(4), denominator

      maximum = max(0.0_dp, maxval(raw_parameters(:3)))
      weights = [exp(-maximum), exp(raw_parameters(1) - maximum), &
         exp(raw_parameters(2) - maximum), exp(raw_parameters(3) - maximum)]
      denominator = sum(weights)
      arch = weights(2)/denominator
      garch = weights(3)/denominator
      asymmetry = weights(4)/denominator
   end subroutine adcc_transform

   pure function adcc_inverse_transform(parameters) result(raw_parameters)
      !! Map physical ADCC coefficients to unconstrained optimizer parameters.
      real(dp), intent(in) :: parameters(:) !! Model parameter values.
      real(dp) :: raw_parameters(3), remainder

      remainder = max(1.0_dp - sum(parameters), tiny(1.0_dp))
      raw_parameters = log(max(parameters, tiny(1.0_dp))/remainder)
   end function adcc_inverse_transform

   pure function mts_tse_tsui_fit(standardized_residuals, window, marginal_variance, &
      initial, max_iterations, tolerance) result(out)
      !! Estimate Gaussian Tse-Tsui rolling-correlation DCC.
      real(dp), intent(in) :: standardized_residuals(:, :) !! Standardized residuals.
      integer, intent(in), optional :: window !! Window.
      real(dp), intent(in), optional :: marginal_variance(:, :) !! Marginal variance.
      real(dp), intent(in), optional :: initial(:) !! Initial value.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      type(mts_tse_tsui_fit_t) :: out

      out = fit_tse_tsui_core(standardized_residuals, .false., window, marginal_variance, &
         initial, max_iterations, tolerance)
   end function mts_tse_tsui_fit

   pure function mts_tse_tsui_t_fit(standardized_residuals, window, marginal_variance, &
      initial, max_iterations, tolerance) result(out)
      !! Estimate standardized Student-t Tse-Tsui rolling-correlation DCC.
      real(dp), intent(in) :: standardized_residuals(:, :) !! Standardized residuals.
      integer, intent(in), optional :: window !! Window.
      real(dp), intent(in), optional :: marginal_variance(:, :) !! Marginal variance.
      real(dp), intent(in), optional :: initial(:) !! Initial value.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      type(mts_tse_tsui_fit_t) :: out

      out = fit_tse_tsui_core(standardized_residuals, .true., window, marginal_variance, &
         initial, max_iterations, tolerance)
   end function mts_tse_tsui_t_fit

   pure function fit_tse_tsui_core(standardized_residuals, use_student_t, window, &
      marginal_variance, initial, max_iterations, tolerance) result(out)
      !! Fit Gaussian or Student-t Tse-Tsui DCC through one rolling recursion.
      real(dp), intent(in) :: standardized_residuals(:, :) !! Standardized residuals.
      logical, intent(in) :: use_student_t !! Whether to use the student t.
      integer, intent(in), optional :: window !! Window.
      real(dp), intent(in), optional :: marginal_variance(:, :) !! Marginal variance.
      real(dp), intent(in), optional :: initial(:) !! Initial value.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      type(mts_tse_tsui_fit_t) :: out
      type(optimization_result_t) :: optimization
      real(dp), allocatable :: physical_initial(:), raw_initial(:), hessian(:, :)
      real(dp), allocatable :: raw_covariance(:, :), jacobian(:, :)
      real(dp) :: gradient_tolerance
      integer :: n, dimension, parameter_count, rolling_window, limit, status, i

      n = size(standardized_residuals, 1)
      dimension = size(standardized_residuals, 2)
      parameter_count = merge(3, 2, use_student_t)
      rolling_window = dimension + 1
      if (present(window)) rolling_window = window
      if (n < 20 .or. dimension < 2 .or. rolling_window < 2 .or. rolling_window >= n - 2) then
         out%info = 1
         return
      end if
      allocate(out%marginal_variance(n, dimension))
      out%marginal_variance = 1.0_dp
      if (present(marginal_variance)) then
         if (any(shape(marginal_variance) /= [n, dimension]) .or. &
            any(marginal_variance <= 0.0_dp)) then
            out%info = 2
            return
         end if
         out%marginal_variance = marginal_variance
      end if
      allocate(physical_initial(parameter_count), raw_initial(parameter_count))
      physical_initial(:2) = [0.9_dp, 0.02_dp]
      if (use_student_t) physical_initial(3) = 8.0_dp
      if (present(initial)) then
         if (size(initial) /= parameter_count .or. any(initial(:2) <= 0.0_dp) .or. &
            sum(initial(:2)) >= 1.0_dp) then
            out%info = 3
            return
         end if
         if (use_student_t) then
            if (initial(3) <= 2.0_dp) then
               out%info = 3
               return
            end if
         end if
         physical_initial = initial
      end if
      raw_initial(:2) = dcc_inverse_transform(physical_initial(:2))
      if (use_student_t) raw_initial(3) = log(physical_initial(3) - 2.0_dp)
      limit = 200
      gradient_tolerance = 1.0e-6_dp
      if (present(max_iterations)) limit = max_iterations
      if (present(tolerance)) gradient_tolerance = tolerance
      optimization = bfgs_minimize_fd(objective, raw_initial, limit, gradient_tolerance)
      call dcc_transform(optimization%parameters, out%previous_weight, out%rolling_weight)
      if (use_student_t) out%degrees_of_freedom = 2.0_dp + exp(min(50.0_dp, optimization%parameters(3)))
      out%iterations = optimization%iterations
      out%converged = optimization%converged
      if (optimization%info /= 0) out%info = 100 + optimization%info
      call evaluate(optimization%parameters, out%unconditional, out%local_correlation, &
         out%correlation, out%covariance, out%residuals, out%log_likelihood, status)
      if (status /= 0) then
         out%info = 10 + status
         return
      end if
      out%standardized_residuals = standardized_residuals
      out%student_t = use_student_t
      out%window = rolling_window
      out%aic = -2.0_dp*out%log_likelihood + 2.0_dp*real(parameter_count, dp)
      out%bic = -2.0_dp*out%log_likelihood + log(real(n, dp))*real(parameter_count, dp)
      hessian = finite_difference_hessian(objective, optimization%parameters)
      call invert_matrix(hessian, raw_covariance, status)
      allocate(out%parameter_covariance(parameter_count, parameter_count))
      allocate(out%standard_errors(parameter_count), jacobian(parameter_count, parameter_count))
      out%parameter_covariance = 0.0_dp
      out%standard_errors = 0.0_dp
      if (status == 0) then
         jacobian = 0.0_dp
         jacobian(1, 1) = out%previous_weight*(1.0_dp - out%previous_weight)
         jacobian(1, 2) = -out%previous_weight*out%rolling_weight
         jacobian(2, 1) = -out%rolling_weight*out%previous_weight
         jacobian(2, 2) = out%rolling_weight*(1.0_dp - out%rolling_weight)
         if (use_student_t) jacobian(3, 3) = out%degrees_of_freedom - 2.0_dp
         out%parameter_covariance = matmul(jacobian, matmul(raw_covariance, transpose(jacobian)))
         do i = 1, parameter_count
            out%standard_errors(i) = sqrt(max(0.0_dp, out%parameter_covariance(i, i)))
         end do
      else if (out%info == 0) then
         out%info = 200 + status
      end if

   contains

      pure function objective(raw_parameters) result(value)
         !! Return the selected Tse-Tsui negative log likelihood.
         real(dp), intent(in) :: raw_parameters(:) !! Raw parameters.
         real(dp) :: value, log_likelihood
         real(dp), allocatable :: unconditional(:, :), local(:, :, :), correlation(:, :, :)
         real(dp), allocatable :: covariance(:, :, :), residuals(:, :)
         integer :: evaluation_status

         call evaluate(raw_parameters, unconditional, local, correlation, covariance, &
            residuals, log_likelihood, evaluation_status)
         if (evaluation_status == 0 .and. ieee_is_finite(log_likelihood)) then
            value = -log_likelihood
         else
            value = 1.0e30_dp + dot_product(raw_parameters, raw_parameters)
         end if
      end function objective

      pure subroutine evaluate(raw_parameters, unconditional, local, correlation, covariance, &
         residuals, log_likelihood, evaluation_status)
         !! Evaluate rolling correlations, covariance paths, and likelihood.
         real(dp), intent(in) :: raw_parameters(:) !! Raw parameters.
         real(dp), allocatable, intent(out) :: unconditional(:, :) !! Unconditional.
         real(dp), allocatable, intent(out) :: local(:, :, :) !! Local.
         real(dp), allocatable, intent(out) :: correlation(:, :, :) !! Correlation.
         real(dp), allocatable, intent(out) :: covariance(:, :, :) !! Covariance matrix.
         real(dp), allocatable, intent(out) :: residuals(:, :) !! Model residuals.
         real(dp), intent(out) :: log_likelihood !! Log-likelihood value.
         integer, intent(out) :: evaluation_status !! Evaluation status.
         real(dp), allocatable :: inverse(:, :)
         real(dp) :: previous, rolling, degrees, constant, logdet, quadratic, independent
         real(dp) :: scale(dimension)
         integer :: t, row, column

         call dcc_transform(raw_parameters, previous, rolling)
         degrees = 0.0_dp
         if (use_student_t) degrees = 2.0_dp + exp(min(50.0_dp, raw_parameters(3)))
         unconditional = data_correlation_matrix(standardized_residuals)
         allocate(local(dimension, dimension, n), correlation(dimension, dimension, n))
         allocate(covariance(dimension, dimension, n), residuals(n, dimension))
         local = spread(unconditional, 3, n)
         correlation = spread(unconditional, 3, n)
         do t = rolling_window + 1, n
            local(:, :, t) = data_correlation_matrix(&
               standardized_residuals(t - rolling_window:t - 1, :))
            correlation(:, :, t) = (1.0_dp - previous - rolling)*unconditional + &
               previous*correlation(:, :, t - 1) + rolling*local(:, :, t)
         end do
         residuals = standardized_residuals*sqrt(out%marginal_variance)
         log_likelihood = 0.0_dp
         if (use_student_t) constant = log_gamma(0.5_dp*degrees) - &
            log_gamma(0.5_dp*(degrees + real(dimension, dp))) + &
            0.5_dp*real(dimension, dp)*log(acos(-1.0_dp)*(degrees - 2.0_dp))
         allocate(inverse(dimension, dimension))
         do t = 1, n
            scale = sqrt(out%marginal_variance(t, :))
            covariance(:, :, t) = correlation(:, :, t)* &
               spread(scale, 1, dimension)*spread(scale, 2, dimension)
            if (t <= rolling_window) cycle
            call inverse_logdet(correlation(:, :, t), inverse, logdet, evaluation_status, &
               100.0_dp*epsilon(1.0_dp))
            if (evaluation_status /= 0) then
               log_likelihood = -huge(1.0_dp)
               return
            end if
            quadratic = dot_product(standardized_residuals(t, :), &
               matmul(inverse, standardized_residuals(t, :)))
            if (use_student_t) then
               log_likelihood = log_likelihood - constant - 0.5_dp*logdet - &
                  0.5_dp*(degrees + real(dimension, dp))*log(1.0_dp + quadratic/(degrees - 2.0_dp))
            else
               independent = dot_product(standardized_residuals(t, :), standardized_residuals(t, :))
               log_likelihood = log_likelihood - 0.5_dp*(logdet + quadratic - independent)
            end if
         end do
         evaluation_status = 0
      end subroutine evaluate

   end function fit_tse_tsui_core

   pure function mts_tse_tsui_forecast(model, future_variance) result(covariance_forecast)
      !! Forecast Tse-Tsui covariance while holding the latest rolling target fixed.
      type(mts_tse_tsui_fit_t), intent(in) :: model !! Model specification.
      real(dp), intent(in) :: future_variance(:, :) !! Future variance.
      real(dp), allocatable :: covariance_forecast(:, :, :), correlation(:, :, :)
      real(dp) :: scale(size(future_variance, 2))
      integer :: horizon, dimension, step, row, column

      horizon = size(future_variance, 1)
      dimension = size(future_variance, 2)
      if (.not. allocated(model%correlation) .or. horizon < 1 .or. &
         dimension /= size(model%unconditional, 1) .or. any(future_variance <= 0.0_dp)) then
         allocate(covariance_forecast(0, 0, 0))
         return
      end if
      allocate(covariance_forecast(dimension, dimension, horizon))
      allocate(correlation(dimension, dimension, horizon))
      correlation(:, :, 1) = (1.0_dp - model%previous_weight - model%rolling_weight)* &
         model%unconditional + model%previous_weight* &
         model%correlation(:, :, size(model%correlation, 3)) + &
         model%rolling_weight*model%local_correlation(:, :, size(model%local_correlation, 3))
      do step = 1, horizon
         scale = sqrt(future_variance(step, :))
         covariance_forecast(:, :, step) = correlation(:, :, step)* &
            spread(scale, 1, dimension)*spread(scale, 2, dimension)
         if (step < horizon) correlation(:, :, step + 1) = &
            (1.0_dp - model%previous_weight - model%rolling_weight)*model%unconditional + &
            model%previous_weight*correlation(:, :, step) + model%rolling_weight* &
            model%local_correlation(:, :, size(model%local_correlation, 3))
      end do
   end function mts_tse_tsui_forecast

   pure function mts_ewma_fit(returns, decay, max_iterations, tolerance) result(out)
      !! Compute or estimate the MTS multivariate EWMA covariance model.
      real(dp), intent(in) :: returns(:, :) !! Returns.
      real(dp), intent(in), optional :: decay !! Decay.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      type(mts_ewma_fit_t) :: out
      type(optimization_result_t) :: optimization
      real(dp) :: requested_decay, raw_initial(1), hessian(1, 1), variance_raw
      real(dp) :: gradient_tolerance
      integer :: n, variables, limit, status

      n = size(returns, 1)
      variables = size(returns, 2)
      requested_decay = 0.96_dp
      if (present(decay)) requested_decay = decay
      if (n < 5 .or. variables < 1 .or. requested_decay >= 1.0_dp) then
         out%info = 1
         return
      end if
      out%mean = sum(returns, 1)/real(n, dp)
      out%residuals = returns - spread(out%mean, 1, n)
      if (requested_decay > 0.0_dp) then
         out%decay = requested_decay
      else
         raw_initial(1) = log(0.96_dp/0.04_dp)
         limit = 150
         gradient_tolerance = 1.0e-7_dp
         if (present(max_iterations)) limit = max_iterations
         if (present(tolerance)) gradient_tolerance = tolerance
         optimization = bfgs_minimize_fd(objective, raw_initial, limit, gradient_tolerance)
         out%decay = logistic(optimization%parameters(1))
         out%iterations = optimization%iterations
         out%converged = optimization%converged
         out%estimated = .true.
         if (optimization%info /= 0) out%info = 100 + optimization%info
         hessian = finite_difference_hessian(objective, optimization%parameters)
         if (hessian(1, 1) > tiny(1.0_dp)) then
            variance_raw = 1.0_dp/hessian(1, 1)
            out%standard_error = out%decay*(1.0_dp - out%decay)*sqrt(variance_raw)
         else if (out%info == 0) then
            out%info = 200
         end if
      end if
      call ewma_path(out%residuals, out%decay, out%covariance, out%log_likelihood, status)
      if (status /= 0) out%info = 10 + status

   contains

      pure function objective(raw_parameter) result(value)
         !! Return the EWMA Gaussian negative log likelihood.
         real(dp), intent(in) :: raw_parameter(:) !! Raw parameter.
         real(dp) :: value, log_likelihood
         real(dp), allocatable :: covariance(:, :, :)
         integer :: evaluation_status

         call ewma_path(out%residuals, logistic(raw_parameter(1)), covariance, &
            log_likelihood, evaluation_status)
         if (evaluation_status == 0) then
            value = -log_likelihood
         else
            value = 1.0e30_dp + dot_product(raw_parameter, raw_parameter)
         end if
      end function objective
   end function mts_ewma_fit

   pure function mts_ewma_forecast(model, horizon) result(covariance_forecast)
      !! Forecast EWMA covariance matrices after the final observed shock.
      type(mts_ewma_fit_t), intent(in) :: model !! Model specification.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      real(dp), allocatable :: covariance_forecast(:, :, :)
      real(dp), allocatable :: shock(:)
      integer :: variables, step

      if (.not. allocated(model%covariance) .or. horizon < 1) then
         allocate(covariance_forecast(0, 0, 0))
         return
      end if
      variables = size(model%mean)
      allocate(covariance_forecast(variables, variables, horizon), shock(variables))
      shock = model%residuals(size(model%residuals, 1), :)
      covariance_forecast(:, :, 1) = (1.0_dp - model%decay)* &
         spread(shock, 2, variables)*spread(shock, 1, variables) + &
         model%decay*model%covariance(:, :, size(model%covariance, 3))
      do step = 2, horizon
         covariance_forecast(:, :, step) = covariance_forecast(:, :, step - 1)
      end do
   end function mts_ewma_forecast

   pure subroutine ewma_path(residuals, decay, covariance, log_likelihood, info)
      !! Evaluate an EWMA covariance recursion and Gaussian likelihood.
      real(dp), intent(in) :: residuals(:, :) !! Model residuals.
      real(dp), intent(in) :: decay !! Decay.
      real(dp), allocatable, intent(out) :: covariance(:, :, :) !! Covariance matrix.
      real(dp), intent(out) :: log_likelihood !! Log-likelihood value.
      integer, intent(out) :: info !! Status code; zero indicates success.
      real(dp), allocatable :: inverse(:, :)
      real(dp) :: shock(size(residuals, 2)), logdet, quadratic
      integer :: n, variables, t

      n = size(residuals, 1)
      variables = size(residuals, 2)
      allocate(covariance(variables, variables, n), inverse(variables, variables))
      covariance(:, :, 1) = matmul(transpose(residuals), residuals)/real(n - 1, dp)
      log_likelihood = 0.0_dp
      do t = 2, n
         shock = residuals(t - 1, :)
         covariance(:, :, t) = (1.0_dp - decay)*spread(shock, 2, variables)* &
            spread(shock, 1, variables) + decay*covariance(:, :, t - 1)
         call inverse_logdet(covariance(:, :, t), inverse, logdet, info, &
            100.0_dp*epsilon(1.0_dp))
         if (info /= 0) return
         quadratic = dot_product(residuals(t, :), matmul(inverse, residuals(t, :)))
         log_likelihood = log_likelihood - 0.5_dp*(real(variables, dp)*log(2.0_dp*acos(-1.0_dp)) + &
            logdet + quadratic)
      end do
      info = 0
   end subroutine ewma_path

   pure elemental real(dp) function logistic(value) result(probability)
      !! Return a numerically stable logistic transform.
      real(dp), intent(in) :: value !! Input value.
      if (value >= 0.0_dp) then
         probability = 1.0_dp/(1.0_dp + exp(-value))
      else
         probability = exp(value)/(1.0_dp + exp(value))
      end if
   end function logistic

   pure function mts_mchol_fit(returns, window, decay, var_order, max_iterations, tolerance) result(out)
      !! Fit the MTS moving-Cholesky model with component GARCH variances.
      real(dp), intent(in) :: returns(:, :) !! Returns.
      integer, intent(in), optional :: window !! Window.
      integer, intent(in), optional :: var_order !! Var order.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      real(dp), intent(in), optional :: decay !! Decay.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      type(mts_mchol_fit_t) :: out
      type(mts_var_fit_t) :: prefit
      real(dp), allocatable :: raw_coefficients(:, :), smoothed(:, :), component(:), variance(:)
      real(dp), allocatable :: triangular(:, :), triangular_inverse(:, :), diagonal(:, :)
      real(dp) :: component_parameters(3), smoothing, gradient_tolerance
      integer :: n, variables, initial_window, order, rows, coefficient_count
      integer :: variable, start, column, t, status, limit

      n = size(returns, 1)
      variables = size(returns, 2)
      initial_window = 36
      order = 0
      smoothing = 0.96_dp
      limit = 150
      gradient_tolerance = 1.0e-6_dp
      if (present(window)) initial_window = window
      if (present(var_order)) order = var_order
      if (present(decay)) smoothing = decay
      if (present(max_iterations)) limit = max_iterations
      if (present(tolerance)) gradient_tolerance = tolerance
      if (variables < 2 .or. order < 0 .or. smoothing <= 0.0_dp .or. smoothing >= 1.0_dp) then
         out%info = 1
         return
      end if
      out%mean = sum(returns, 1)/real(n, dp)
      if (order > 0) then
         prefit = mts_var(returns, order)
         if (prefit%info /= 0) then
            out%info = 10 + prefit%info
            return
         end if
         out%prewhitened = prefit%residuals
      else
         out%prewhitened = returns - spread(out%mean, 1, n)
      end if
      n = size(out%prewhitened, 1)
      if (initial_window <= variables .or. initial_window >= n - 5) then
         out%info = 2
         return
      end if
      rows = n - initial_window
      coefficient_count = variables*(variables - 1)/2
      allocate(out%orthogonal_residuals(rows, variables), out%component_variance(rows, variables))
      allocate(out%garch_parameters(variables, 3), out%smoothed_coefficients(rows, coefficient_count))
      out%orthogonal_residuals(:, 1) = out%prewhitened(initial_window + 1:, 1)
      call fit_component_garch(out%orthogonal_residuals(:, 1), limit, gradient_tolerance, &
         component_parameters, variance, status)
      if (status /= 0) then
         out%info = 20 + status
         return
      end if
      out%component_variance(:, 1) = variance
      out%garch_parameters(1, :) = component_parameters
      start = 0
      do variable = 2, variables
         call recursive_coefficients(out%prewhitened(:, variable), &
            out%prewhitened(:, :variable - 1), initial_window, raw_coefficients, status)
         if (status /= 0) then
            out%info = 30 + variable
            return
         end if
         smoothed = raw_coefficients
         do column = 1, variable - 1
            smoothed(1, column) = raw_coefficients(1, column)
            do t = 2, rows
               smoothed(t, column) = smoothing*smoothed(t - 1, column) + &
                  (1.0_dp - smoothing)*raw_coefficients(t, column)
            end do
         end do
         out%smoothed_coefficients(:, start + 1:start + variable - 1) = -smoothed
         component = out%prewhitened(initial_window + 1:, variable)
         do column = 1, variable - 1
            component = component - smoothed(:, column)* &
               out%prewhitened(initial_window + 1:, column)
         end do
         out%orthogonal_residuals(:, variable) = component
         call fit_component_garch(component, limit, gradient_tolerance, component_parameters, variance, status)
         if (status /= 0) then
            out%info = 40 + variable
            return
         end if
         out%component_variance(:, variable) = variance
         out%garch_parameters(variable, :) = component_parameters
         start = start + variable - 1
      end do
      allocate(out%covariance(variables, variables, rows), triangular(variables, variables))
      allocate(diagonal(variables, variables))
      do t = 1, rows
         triangular = 0.0_dp
         do variable = 1, variables
            triangular(variable, variable) = 1.0_dp
         end do
         start = 0
         do variable = 2, variables
            triangular(variable, :variable - 1) = &
               out%smoothed_coefficients(t, start + 1:start + variable - 1)
            start = start + variable - 1
         end do
         call invert_matrix(triangular, triangular_inverse, status)
         if (status /= 0) then
            out%info = 50 + t
            return
         end if
         diagonal = 0.0_dp
         do variable = 1, variables
            diagonal(variable, variable) = out%component_variance(t, variable)
         end do
         out%covariance(:, :, t) = matmul(triangular_inverse, &
            matmul(diagonal, transpose(triangular_inverse)))
      end do
      out%window = initial_window
      out%var_order = order
      out%decay = smoothing
   end function mts_mchol_fit

   pure function mts_mchol_forecast(model, horizon) result(covariance_forecast)
      !! Forecast moving-Cholesky covariances with final smoothed coefficients.
      type(mts_mchol_fit_t), intent(in) :: model !! Model specification.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      real(dp), allocatable :: covariance_forecast(:, :, :)
      real(dp), allocatable :: triangular(:, :), inverse(:, :), diagonal(:, :), variance(:, :)
      integer :: variables, variable, start, step, status

      if (.not. allocated(model%covariance) .or. horizon < 1) then
         allocate(covariance_forecast(0, 0, 0))
         return
      end if
      variables = size(model%garch_parameters, 1)
      allocate(covariance_forecast(variables, variables, horizon), variance(horizon, variables))
      do variable = 1, variables
         variance(1, variable) = model%garch_parameters(variable, 1) + &
            model%garch_parameters(variable, 2)* &
            model%orthogonal_residuals(size(model%orthogonal_residuals, 1), variable)**2 + &
            model%garch_parameters(variable, 3)* &
            model%component_variance(size(model%component_variance, 1), variable)
         do step = 2, horizon
            variance(step, variable) = model%garch_parameters(variable, 1) + &
               sum(model%garch_parameters(variable, 2:3))*variance(step - 1, variable)
         end do
      end do
      allocate(triangular(variables, variables), diagonal(variables, variables))
      triangular = 0.0_dp
      do variable = 1, variables
         triangular(variable, variable) = 1.0_dp
      end do
      start = 0
      do variable = 2, variables
         triangular(variable, :variable - 1) = model%smoothed_coefficients(&
            size(model%smoothed_coefficients, 1), start + 1:start + variable - 1)
         start = start + variable - 1
      end do
      call invert_matrix(triangular, inverse, status)
      if (status /= 0) then
         deallocate(covariance_forecast)
         allocate(covariance_forecast(0, 0, 0))
         return
      end if
      do step = 1, horizon
         diagonal = 0.0_dp
         do variable = 1, variables
            diagonal(variable, variable) = variance(step, variable)
         end do
         covariance_forecast(:, :, step) = matmul(inverse, matmul(diagonal, transpose(inverse)))
      end do
   end function mts_mchol_forecast

   pure subroutine recursive_coefficients(response, predictors, window, coefficients, info)
      !! Compute expanding-window recursive least-squares coefficients.
      real(dp), intent(in) :: response(:) !! Response observations.
      real(dp), intent(in) :: predictors(:, :) !! Predictor matrix.
      integer, intent(in) :: window !! Window.
      real(dp), allocatable, intent(out) :: coefficients(:, :) !! Model coefficients.
      integer, intent(out) :: info !! Status code; zero indicates success.
      real(dp), allocatable :: cross(:, :), inverse(:, :), beta(:), gain(:), predictor(:)
      real(dp) :: denominator, error
      integer :: n, count_predictors, t, row, status

      n = size(response)
      count_predictors = size(predictors, 2)
      cross = matmul(transpose(predictors(:window, :)), predictors(:window, :))
      call invert_matrix(cross, inverse, status)
      if (status /= 0) then
         info = status
         return
      end if
      beta = matmul(inverse, matmul(transpose(predictors(:window, :)), response(:window)))
      allocate(coefficients(n - window, count_predictors), predictor(count_predictors))
      do t = window + 1, n
         row = t - window
         predictor = predictors(t, :)
         gain = matmul(inverse, predictor)
         denominator = 1.0_dp + dot_product(predictor, gain)
         gain = gain/denominator
         error = response(t) - dot_product(predictor, beta)
         beta = beta + gain*error
         inverse = inverse - spread(gain, 2, count_predictors)* &
            spread(matmul(predictor, inverse), 1, count_predictors)
         coefficients(row, :) = beta
      end do
      info = 0
   end subroutine recursive_coefficients

   pure subroutine fit_component_garch(series, max_iterations, tolerance, parameters, variance, info)
      !! Fit a Gaussian GARCH(1,1) variance path for one orthogonal component.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      integer, intent(in) :: max_iterations !! Maximum number of algorithm iterations.
      real(dp), intent(in) :: tolerance !! Numerical convergence tolerance.
      real(dp), intent(out) :: parameters(3) !! Model parameter values.
      real(dp), allocatable, intent(out) :: variance(:) !! Variance value or matrix.
      integer, intent(out) :: info !! Status code; zero indicates success.
      type(optimization_result_t) :: optimization
      real(dp) :: initial(3), sample_variance

      sample_variance = sum(series**2)/real(size(series), dp)
      initial(1) = log(max(0.05_dp*sample_variance, tiny(1.0_dp)))
      initial(2:3) = dcc_inverse_transform([0.08_dp, 0.88_dp])
      optimization = bfgs_minimize_fd(objective, initial, max_iterations, tolerance)
      call transform_garch(optimization%parameters, parameters)
      call variance_path(parameters, variance, info)

   contains

      pure function objective(raw_parameters) result(value)
         !! Return a univariate Gaussian GARCH negative log likelihood.
         real(dp), intent(in) :: raw_parameters(:) !! Raw parameters.
         real(dp) :: value, physical(3)
         real(dp), allocatable :: path(:)
         integer :: evaluation_status

         call transform_garch(raw_parameters, physical)
         call variance_path(physical, path, evaluation_status)
         if (evaluation_status == 0) then
            value = 0.5_dp*sum(log(2.0_dp*acos(-1.0_dp)*path) + series**2/path)
         else
            value = 1.0e30_dp + dot_product(raw_parameters, raw_parameters)
         end if
      end function objective

      pure subroutine transform_garch(raw_parameters, physical)
         !! Map unconstrained values to positive stationary GARCH parameters.
         real(dp), intent(in) :: raw_parameters(:) !! Raw parameters.
         real(dp), intent(out) :: physical(3) !! Physical.
         real(dp) :: arch_value, garch_value

         physical(1) = exp(min(50.0_dp, raw_parameters(1)))
         call dcc_transform(raw_parameters(2:3), arch_value, garch_value)
         physical(2:3) = [arch_value, garch_value]
      end subroutine transform_garch

      pure subroutine variance_path(physical, path, evaluation_status)
         !! Evaluate the GARCH variance recursion.
         real(dp), intent(in) :: physical(3) !! Physical.
         real(dp), allocatable, intent(out) :: path(:) !! Path.
         integer, intent(out) :: evaluation_status !! Evaluation status.
         integer :: t

         allocate(path(size(series)))
         path(1) = max(sample_variance, physical(1)/(1.0_dp - physical(2) - physical(3)))
         do t = 2, size(series)
            path(t) = physical(1) + physical(2)*series(t - 1)**2 + physical(3)*path(t - 1)
            if (.not. ieee_is_finite(path(t)) .or. path(t) <= tiny(1.0_dp)) then
               evaluation_status = 1
               return
            end if
         end do
         evaluation_status = 0
      end subroutine variance_path
   end subroutine fit_component_garch

   pure function mts_sccor(series, end_index, span, groups) result(out)
      !! Estimate MTS group-constrained correlations over a selected sample window.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      integer, intent(in) :: end_index !! Index of end.
      integer, intent(in) :: span !! Span.
      integer, intent(in) :: groups(:) !! Groups.
      type(mts_sccor_t) :: out
      integer :: n, variables, ending, window, required, group_count
      integer :: first_group, second_group, first_start, second_start
      integer :: i, j, count_pairs
      real(dp) :: average

      n = size(series, 1)
      variables = size(series, 2)
      if (size(groups) < 1 .or. any(groups < 1) .or. sum(groups) /= variables) then
         out%info = 1
         return
      end if
      ending = end_index
      window = span
      if (ending < window .or. ending > n) ending = n
      group_count = size(groups)
      required = group_count*(group_count - 1)/2
      required = required + count(groups > 1)
      if (maxval(groups) > 1) window = max(window, required + 1)
      if (window < 2 .or. window > n) then
         out%info = 2
         return
      end if
      if (ending < window) ending = n
      out%start = ending - window + 1
      out%end = ending
      out%span = window
      out%unconstrained = data_correlation_matrix(series(out%start:ending, :))
      out%constrained = out%unconstrained
      first_start = 1
      do first_group = 1, group_count
         if (groups(first_group) > 1) then
            average = 0.0_dp
            count_pairs = 0
            do i = first_start, first_start + groups(first_group) - 2
               do j = i + 1, first_start + groups(first_group) - 1
                  average = average + out%unconstrained(i, j)
                  count_pairs = count_pairs + 1
               end do
            end do
            average = average/real(count_pairs, dp)
            do i = first_start, first_start + groups(first_group) - 2
               do j = i + 1, first_start + groups(first_group) - 1
                  out%constrained(i, j) = average
                  out%constrained(j, i) = average
               end do
            end do
         end if
         second_start = first_start + groups(first_group)
         do second_group = first_group + 1, group_count
            average = sum(out%unconstrained(first_start:first_start + groups(first_group) - 1, &
               second_start:second_start + groups(second_group) - 1))/ &
               real(groups(first_group)*groups(second_group), dp)
            out%constrained(first_start:first_start + groups(first_group) - 1, &
               second_start:second_start + groups(second_group) - 1) = average
            out%constrained(second_start:second_start + groups(second_group) - 1, &
               first_start:first_start + groups(first_group) - 1) = average
            second_start = second_start + groups(second_group)
         end do
         first_start = first_start + groups(first_group)
      end do
   end function mts_sccor

   pure function mts_arch_test(series, lag) result(out)
      !! Perform the MTS McLeod-Li and rank-based ARCH tests.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      integer, intent(in) :: lag !! Lag index or number of lags.
      type(mts_arch_test_t) :: out

      call serial_dependence_tests(series**2, lag, out%statistic, out%p_value, out%info)
      if (out%info == 0) out%lag = lag
   end function mts_arch_test

   pure function mts_march_test(series, lag, robust_probability) result(out)
      !! Perform the four MTS multivariate ARCH tests.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      integer, intent(in) :: lag !! Lag index or number of lags.
      real(dp), intent(in), optional :: robust_probability !! Robust probability.
      type(mts_march_test_t) :: out
      type(mts_mch_diagnostic_t) :: diagnostic
      real(dp), allocatable :: centered(:, :), covariance(:, :), covariance_path(:, :, :)
      integer :: n

      n = size(series, 1)
      if (n < 5 .or. size(series, 2) < 2) then
         out%info = 1
         return
      end if
      centered = series - spread(sum(series, 1)/real(n, dp), 1, n)
      covariance = matmul(transpose(centered), centered)/real(n - 1, dp)
      covariance_path = spread(covariance, 3, n)
      diagnostic = mts_mch_diagnostic(centered, covariance_path, lag, robust_probability)
      if (diagnostic%info /= 0) then
         out%info = diagnostic%info
         return
      end if
      out%statistic = [diagnostic%radial_q, diagnostic%rank_q, &
         diagnostic%multivariate_q, diagnostic%robust_q]
      out%p_value = diagnostic%p_value
      out%degrees_of_freedom = diagnostic%degrees_of_freedom
      out%lag = lag
      out%robust_observations = diagnostic%robust_observations
   end function mts_march_test

   pure subroutine serial_dependence_tests(values, lag, statistic, p_value, info)
      !! Compute Ljung-Box and rank-based serial-dependence tests.
      real(dp), intent(in) :: values(:) !! Input values.
      integer, intent(in) :: lag !! Lag index or number of lags.
      real(dp), intent(out) :: statistic(2) !! Statistic.
      real(dp), intent(out) :: p_value(2) !! P value.
      integer, intent(out) :: info !! Status code; zero indicates success.
      real(dp), allocatable :: autocorrelation(:), ranks(:)
      real(dp) :: mean_value, denominator, mu, rank_variance
      integer :: n, current_lag, t

      n = size(values)
      statistic = 0.0_dp
      p_value = 1.0_dp
      if (lag < 1 .or. lag >= n - 1) then
         info = 1
         return
      end if
      allocate(autocorrelation(lag), ranks(n))
      mean_value = sum(values)/real(n, dp)
      denominator = sum((values - mean_value)**2)
      if (denominator <= tiny(1.0_dp)) then
         info = 2
         return
      end if
      do current_lag = 1, lag
         autocorrelation(current_lag) = sum((values(current_lag + 1:) - mean_value)* &
            (values(:n - current_lag) - mean_value))/denominator
         statistic(1) = statistic(1) + autocorrelation(current_lag)**2/real(n - current_lag, dp)
      end do
      statistic(1) = statistic(1)*real(n, dp)*real(n + 2, dp)
      do t = 1, n
         ranks(t) = 1.0_dp + real(count(values < values(t)), dp) + &
            0.5_dp*real(count(values == values(t)) - 1, dp)
      end do
      mean_value = sum(ranks)/real(n, dp)
      denominator = sum((ranks - mean_value)**2)
      do current_lag = 1, lag
         autocorrelation(current_lag) = sum((ranks(current_lag + 1:) - mean_value)* &
            (ranks(:n - current_lag) - mean_value))/denominator
         mu = -real(n - current_lag, dp)/(real(n, dp)*real(n - 1, dp))
         rank_variance = (5.0_dp*real(n, dp)**4 - real(5*current_lag + 9, dp)*real(n, dp)**3 + &
            9.0_dp*real(current_lag - 2, dp)*real(n, dp)**2 + &
            2.0_dp*real(current_lag*(5*current_lag + 8), dp)*real(n, dp) + &
            16.0_dp*real(current_lag*current_lag, dp))/ &
            (5.0_dp*real(n - 1, dp)**2*real(n, dp)**2*real(n + 1, dp))
         statistic(2) = statistic(2) + (autocorrelation(current_lag) - mu)**2/rank_variance
      end do
      p_value(1) = chi_square_survival(statistic(1), lag)
      p_value(2) = chi_square_survival(statistic(2), lag)
      info = 0
   end subroutine serial_dependence_tests

   pure function mts_correlation_to_angles(correlation) result(angles)
      !! Convert a positive-definite correlation matrix to hyperspherical angles.
      real(dp), intent(in) :: correlation(:, :) !! Correlation.
      real(dp), allocatable :: angles(:)
      real(dp), allocatable :: lower(:, :)
      real(dp) :: product
      integer :: dimension, row, column, index, status

      dimension = size(correlation, 1)
      allocate(angles(dimension*(dimension - 1)/2))
      call cholesky_lower(correlation, lower, status)
      if (status /= 0) then
         angles = 0.0_dp
         return
      end if
      index = 0
      do row = 2, dimension
         product = 1.0_dp
         do column = 1, row - 1
            index = index + 1
            angles(index) = acos(max(-1.0_dp, min(1.0_dp, lower(row, column)/product)))
            product = product*sin(angles(index))
         end do
      end do
   end function mts_correlation_to_angles

   pure function mts_angles_to_correlation(angles, dimension) result(correlation)
      !! Convert full hyperspherical angles to a correlation matrix.
      real(dp), intent(in) :: angles(:) !! Angles.
      integer, intent(in) :: dimension !! Dimension.
      real(dp), allocatable :: correlation(:, :)
      real(dp), allocatable :: lower(:, :)
      real(dp) :: product
      integer :: row, column, index

      if (dimension < 1 .or. size(angles) /= dimension*(dimension - 1)/2) then
         allocate(correlation(0, 0))
         return
      end if
      allocate(lower(dimension, dimension))
      lower = 0.0_dp
      lower(1, 1) = 1.0_dp
      index = 0
      do row = 2, dimension
         product = 1.0_dp
         do column = 1, row - 1
            index = index + 1
            lower(row, column) = product*cos(angles(index))
            product = product*sin(angles(index))
         end do
         lower(row, row) = product
      end do
      correlation = matmul(lower, transpose(lower))
   end function mts_angles_to_correlation

   pure function mts_copula_fit(standardized_residuals, groups, window, baseline_angles, &
      estimate_baseline, initial, max_iterations, tolerance) result(out)
      !! Estimate the grouped dynamic Student-t copula from MTS mtCopula.
      real(dp), intent(in) :: standardized_residuals(:, :) !! Standardized residuals.
      integer, intent(in) :: groups(:) !! Groups.
      integer, intent(in), optional :: window !! Window.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      real(dp), intent(in), optional :: baseline_angles(:) !! Baseline angles.
      real(dp), intent(in), optional :: initial(:) !! Initial value.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      logical, intent(in), optional :: estimate_baseline !! Whether to estimate the baseline.
      type(mts_copula_fit_t) :: out
      type(optimization_result_t) :: optimization
      type(mts_sccor_t) :: grouped
      real(dp), allocatable :: raw_initial(:), hessian(:, :), inverse(:, :)
      real(dp), allocatable :: full_angles(:), local_correlation(:, :)
      real(dp) :: gradient_tolerance
      integer :: n, dimension, angle_count, rolling_window, parameter_count
      integer :: t, limit, status, i
      logical :: fit_baseline

      n = size(standardized_residuals, 1)
      dimension = size(standardized_residuals, 2)
      angle_count = size(groups)*(size(groups) - 1)/2 + count(groups > 1)
      rolling_window = angle_count + 1
      fit_baseline = .true.
      if (present(window)) rolling_window = window
      if (present(estimate_baseline)) fit_baseline = estimate_baseline
      if (n < 20 .or. size(groups) < 1 .or. any(groups < 1) .or. &
         sum(groups) /= dimension .or. rolling_window < 2 .or. rolling_window >= n - 2) then
         out%info = 1
         return
      end if
      if (present(baseline_angles)) then
         if (size(baseline_angles) /= angle_count) then
            out%info = 2
            return
         end if
         out%baseline_angles = baseline_angles
      else
         allocate(out%baseline_angles(angle_count))
         grouped = mts_sccor(standardized_residuals, n, n, groups)
         full_angles = mts_correlation_to_angles(grouped%constrained)
         call extract_group_angles(full_angles, dimension, groups, out%baseline_angles)
      end if
      allocate(out%local_angles(n, angle_count))
      out%local_angles = spread(out%baseline_angles, 1, n)
      do t = rolling_window + 1, n
         grouped = mts_sccor(standardized_residuals, t - 1, rolling_window, groups)
         full_angles = mts_correlation_to_angles(grouped%constrained)
         call extract_group_angles(full_angles, dimension, groups, out%local_angles(t, :))
      end do
      parameter_count = 3 + merge(angle_count, 0, fit_baseline)
      allocate(raw_initial(parameter_count))
      raw_initial(1) = log(7.0_dp - 2.0_dp)
      raw_initial(2:3) = dcc_inverse_transform([0.9_dp, 0.02_dp])
      if (fit_baseline) then
         out%baseline_angles = max(1.0e-6_dp, min(acos(-1.0_dp) - 1.0e-6_dp, out%baseline_angles))
         raw_initial(4:) = log(out%baseline_angles/(acos(-1.0_dp) - out%baseline_angles))
      end if
      if (present(initial)) then
         if (size(initial) /= parameter_count) then
            out%info = 3
            return
         end if
         raw_initial = initial
      end if
      limit = 250
      gradient_tolerance = 1.0e-6_dp
      if (present(max_iterations)) limit = max_iterations
      if (present(tolerance)) gradient_tolerance = tolerance
      optimization = bfgs_minimize_fd(objective, raw_initial, limit, gradient_tolerance)
      out%iterations = optimization%iterations
      out%converged = optimization%converged
      if (optimization%info /= 0) out%info = 100 + optimization%info
      call evaluate(optimization%parameters, out%angles, out%correlation, out%log_likelihood, status)
      if (status /= 0) then
         out%info = 10 + status
         return
      end if
      out%degrees_of_freedom = 2.0_dp + exp(min(50.0_dp, optimization%parameters(1)))
      call dcc_transform(optimization%parameters(2:3), out%previous_weight, out%local_weight)
      if (fit_baseline) out%baseline_angles = acos(-1.0_dp)*logistic(optimization%parameters(4:))
      if (fit_baseline) then
         out%coefficients = [out%degrees_of_freedom, out%previous_weight, &
            out%local_weight, out%baseline_angles]
      else
         out%coefficients = [out%degrees_of_freedom, out%previous_weight, out%local_weight]
      end if
      out%groups = groups
      out%window = rolling_window
      out%estimates_baseline = fit_baseline
      hessian = finite_difference_hessian(objective, optimization%parameters)
      call invert_matrix(hessian, inverse, status)
      allocate(out%parameter_covariance(parameter_count, parameter_count), &
         out%standard_errors(parameter_count))
      out%parameter_covariance = 0.0_dp
      out%standard_errors = 0.0_dp
      if (status == 0) then
         out%parameter_covariance = inverse
         out%standard_errors = sqrt(max(0.0_dp, [(inverse(i, i), i=1, parameter_count)]))
      else if (out%info == 0) then
         out%info = 200 + status
      end if

   contains

      pure function objective(raw_parameters) result(value)
         !! Return the dynamic Student-t copula negative log likelihood.
         real(dp), intent(in) :: raw_parameters(:) !! Raw parameters.
         real(dp) :: value, log_likelihood
         real(dp), allocatable :: angles(:, :), correlation(:, :, :)
         integer :: evaluation_status

         call evaluate(raw_parameters, angles, correlation, log_likelihood, evaluation_status)
         if (evaluation_status == 0 .and. ieee_is_finite(log_likelihood)) then
            value = -log_likelihood
         else
            value = 1.0e30_dp + dot_product(raw_parameters, raw_parameters)
         end if
      end function objective

      pure subroutine evaluate(raw_parameters, angle_path, correlation_path, log_likelihood, evaluation_status)
         !! Evaluate dynamic grouped angles and the standardized t-copula likelihood.
         real(dp), intent(in) :: raw_parameters(:) !! Raw parameters.
         real(dp), allocatable, intent(out) :: angle_path(:, :) !! Angle path.
         real(dp), allocatable, intent(out) :: correlation_path(:, :, :) !! Correlation path.
         real(dp), intent(out) :: log_likelihood !! Log-likelihood value.
         integer, intent(out) :: evaluation_status !! Evaluation status.
         real(dp), allocatable :: baseline(:), inverse(:, :)
         real(dp) :: degrees, previous, local, logdet, quadratic, constant, marginal
         integer :: observation

         degrees = 2.0_dp + exp(min(50.0_dp, raw_parameters(1)))
         call dcc_transform(raw_parameters(2:3), previous, local)
         if (fit_baseline) then
            baseline = acos(-1.0_dp)*logistic(raw_parameters(4:))
         else
            baseline = out%baseline_angles
         end if
         allocate(angle_path(n, angle_count), correlation_path(dimension, dimension, n))
         angle_path = spread(baseline, 1, n)
         call grouped_angles_to_correlation(baseline, groups, correlation_path(:, :, 1), evaluation_status)
         if (evaluation_status /= 0) return
         do observation = 2, rolling_window
            correlation_path(:, :, observation) = correlation_path(:, :, 1)
         end do
         do observation = rolling_window + 1, n
            angle_path(observation, :) = (1.0_dp - previous - local)*baseline + &
               previous*angle_path(observation - 1, :) + local*out%local_angles(observation, :)
            call grouped_angles_to_correlation(angle_path(observation, :), groups, &
               correlation_path(:, :, observation), evaluation_status)
            if (evaluation_status /= 0) return
         end do
         allocate(inverse(dimension, dimension))
         constant = log_gamma(0.5_dp*(degrees + real(dimension, dp))) + &
            real(dimension - 1, dp)*log_gamma(0.5_dp*degrees) - &
            real(dimension, dp)*log_gamma(0.5_dp*(degrees + 1.0_dp))
         log_likelihood = 0.0_dp
         do observation = rolling_window + 1, n
            call inverse_logdet(correlation_path(:, :, observation), inverse, logdet, &
               evaluation_status, 100.0_dp*epsilon(1.0_dp))
            if (evaluation_status /= 0) return
            quadratic = dot_product(standardized_residuals(observation, :), &
               matmul(inverse, standardized_residuals(observation, :)))
            marginal = sum(log(1.0_dp + standardized_residuals(observation, :)**2/(degrees - 2.0_dp)))
            log_likelihood = log_likelihood + constant - 0.5_dp*logdet - &
               0.5_dp*(degrees + real(dimension, dp))*log(1.0_dp + quadratic/(degrees - 2.0_dp)) + &
               0.5_dp*(degrees + 1.0_dp)*marginal
         end do
         evaluation_status = 0
      end subroutine evaluate
   end function mts_copula_fit

   pure subroutine extract_group_angles(full_angles, dimension, groups, group_angles)
      !! Extract within- and between-group representative correlation angles.
      real(dp), intent(in) :: full_angles(:) !! Full angles.
      integer, intent(in) :: dimension !! Dimension.
      integer, intent(in) :: groups(:) !! Groups.
      real(dp), intent(out) :: group_angles(:) !! Group angles.
      real(dp), allocatable :: correlation(:, :)
      real(dp) :: average
      integer :: angle_count, first_group, second_group, first_start, second_start, index

      angle_count = size(groups)*(size(groups) - 1)/2 + count(groups > 1)
      correlation = mts_angles_to_correlation(full_angles, dimension)
      index = 0
      first_start = 1
      do first_group = 1, size(groups)
         if (groups(first_group) > 1) then
            index = index + 1
            average = sum(correlation(first_start:first_start + groups(first_group) - 1, &
               first_start:first_start + groups(first_group) - 1)) - real(groups(first_group), dp)
            average = average/real(groups(first_group)*(groups(first_group) - 1), dp)
            group_angles(index) = acos(max(-1.0_dp, min(1.0_dp, average)))
         end if
         second_start = first_start + groups(first_group)
         do second_group = first_group + 1, size(groups)
            index = index + 1
            average = sum(correlation(first_start:first_start + groups(first_group) - 1, &
               second_start:second_start + groups(second_group) - 1))/ &
               real(groups(first_group)*groups(second_group), dp)
            group_angles(index) = acos(max(-1.0_dp, min(1.0_dp, average)))
            second_start = second_start + groups(second_group)
         end do
         first_start = first_start + groups(first_group)
      end do
   end subroutine extract_group_angles

   pure subroutine grouped_angles_to_correlation(angles, groups, correlation, info)
      !! Construct a positive-definite block-equicorrelation matrix from group angles.
      real(dp), intent(in) :: angles(:) !! Angles.
      integer, intent(in) :: groups(:) !! Groups.
      real(dp), intent(out) :: correlation(:, :) !! Correlation.
      integer, intent(out) :: info !! Status code; zero indicates success.
      real(dp), allocatable :: lower(:, :)
      real(dp) :: value
      integer :: first_group, second_group, first_start, second_start, index, i

      correlation = 0.0_dp
      do i = 1, size(correlation, 1)
         correlation(i, i) = 1.0_dp
      end do
      index = 0
      first_start = 1
      do first_group = 1, size(groups)
         if (groups(first_group) > 1) then
            index = index + 1
            value = cos(angles(index))
            correlation(first_start:first_start + groups(first_group) - 1, &
               first_start:first_start + groups(first_group) - 1) = value
            do i = first_start, first_start + groups(first_group) - 1
               correlation(i, i) = 1.0_dp
            end do
         end if
         second_start = first_start + groups(first_group)
         do second_group = first_group + 1, size(groups)
            index = index + 1
            value = cos(angles(index))
            correlation(first_start:first_start + groups(first_group) - 1, &
               second_start:second_start + groups(second_group) - 1) = value
            correlation(second_start:second_start + groups(second_group) - 1, &
               first_start:first_start + groups(first_group) - 1) = value
            second_start = second_start + groups(second_group)
         end do
         first_start = first_start + groups(first_group)
      end do
      call cholesky_lower(correlation, lower, info)
   end subroutine grouped_angles_to_correlation

   pure function mts_var_missing(series, pi_weights, sigma, time_index, intercept) result(out)
      !! Estimate a completely missing vector observation by conditional GLS.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      real(dp), intent(in) :: pi_weights(:, :, :) !! Pi weights.
      real(dp), intent(in) :: sigma(:, :) !! Scale parameter or standard deviation.
      integer, intent(in) :: time_index !! Index of time.
      real(dp), intent(in), optional :: intercept(:) !! Model intercept.
      type(mts_missing_result_t) :: out
      logical :: missing(size(series, 2))

      missing = .true.
      out = mts_var_partial_missing(series, pi_weights, sigma, time_index, missing, intercept)
   end function mts_var_missing

   pure function mts_var_partial_missing(series, pi_weights, sigma, time_index, &
      missing, intercept) result(out)
      !! Estimate selected missing components of one vector observation by GLS.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      real(dp), intent(in) :: pi_weights(:, :, :) !! Pi weights.
      real(dp), intent(in) :: sigma(:, :) !! Scale parameter or standard deviation.
      integer, intent(in) :: time_index !! Index of time.
      logical, intent(in) :: missing(:) !! Flag controlling missing.
      real(dp), intent(in), optional :: intercept(:) !! Model intercept.
      type(mts_missing_result_t) :: out
      real(dp), allocatable :: sigma_inverse(:, :), normal(:, :), right_hand(:)
      real(dp), allocatable :: design(:, :), base(:), constant(:), work_series(:, :)
      integer, allocatable :: missing_index(:)
      integer :: n, variables, lag_count, missing_count, lag, affected_time, status, i

      n = size(series, 1)
      variables = size(series, 2)
      lag_count = size(pi_weights, 3)
      if (time_index < 1 .or. time_index > n .or. size(missing) /= variables .or. &
         count(missing) < 1 .or. any(shape(sigma) /= [variables, variables]) .or. &
         size(pi_weights, 1) /= variables .or. size(pi_weights, 2) /= variables) then
         out%info = 1
         return
      end if
      allocate(constant(variables))
      constant = 0.0_dp
      if (present(intercept)) then
         if (size(intercept) /= variables) then
            out%info = 2
            return
         end if
         constant = intercept
      end if
      call invert_matrix(sigma, sigma_inverse, status)
      if (status /= 0) then
         out%info = 10 + status
         return
      end if
      missing_count = count(missing)
      allocate(missing_index(missing_count))
      missing_index = pack([(i, i=1, variables)], missing)
      allocate(normal(missing_count, missing_count), right_hand(missing_count))
      allocate(design(variables, missing_count), base(variables), work_series(n, variables))
      normal = 0.0_dp
      right_hand = 0.0_dp
      work_series = series
      work_series(time_index, missing_index) = 0.0_dp

      design = 0.0_dp
      do i = 1, missing_count
         design(missing_index(i), i) = 1.0_dp
      end do
      base = work_series(time_index, :) - constant
      do lag = 1, min(lag_count, time_index - 1)
         base = base - matmul(pi_weights(:, :, lag), work_series(time_index - lag, :))
      end do
      normal = normal + matmul(transpose(design), matmul(sigma_inverse, design))
      right_hand = right_hand - matmul(transpose(design), matmul(sigma_inverse, base))

      do affected_time = time_index + 1, min(n, time_index + lag_count)
         lag = affected_time - time_index
         design = -pi_weights(:, missing_index, lag)
         base = work_series(affected_time, :) - constant
         do i = 1, min(lag_count, affected_time - 1)
            base = base - matmul(pi_weights(:, :, i), work_series(affected_time - i, :))
         end do
         normal = normal + matmul(transpose(design), matmul(sigma_inverse, design))
         right_hand = right_hand - matmul(transpose(design), matmul(sigma_inverse, base))
      end do
      call invert_matrix(normal, out%covariance, status)
      if (status /= 0) then
         out%info = 20 + status
         return
      end if
      out%estimate = matmul(out%covariance, right_hand)
      out%completed = series
      out%completed(time_index, missing_index) = out%estimate
      out%missing = missing
      out%time_index = time_index
   end function mts_var_partial_missing

   pure function mts_granger_test(series, order, targets, include_mean) result(out)
      !! Wald-test whether complementary variables Granger-cause target equations.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      integer, intent(in) :: order !! Model or polynomial order.
      integer, intent(in) :: targets(:) !! Targets.
      logical, intent(in), optional :: include_mean !! Whether to include a mean term.
      type(mts_granger_test_t) :: out
      real(dp), allocatable :: design(:, :), cross(:, :), inverse(:, :), covariance_inverse(:, :)
      integer :: n, variables, rows, offset, lag, i, j, restriction, status
      integer :: target, predictor, first_regressor, restriction_i, restriction_j
      logical :: use_mean

      n = size(series, 1)
      variables = size(series, 2)
      use_mean = .true.
      if (present(include_mean)) use_mean = include_mean
      if (order < 1 .or. size(targets) < 1 .or. any(targets < 1) .or. &
         any(targets > variables) .or. size(unique_indices(targets)) /= size(targets) .or. &
         size(targets) >= variables .or. n <= order + variables + 1) then
         out%info = 1
         return
      end if
      out%targets = targets
      out%predictors = pack([(i, i=1, variables)], &
         [(all(targets /= i), i=1, variables)])
      out%unrestricted = mts_var(series, order, use_mean)
      if (out%unrestricted%info /= 0) then
         out%info = 10 + out%unrestricted%info
         return
      end if
      rows = n - order
      offset = merge(1, 0, use_mean)
      allocate(design(rows, offset + order*variables))
      if (use_mean) design(:, 1) = 1.0_dp
      do lag = 1, order
         design(:, offset + (lag - 1)*variables + 1:offset + lag*variables) = &
            series(order + 1 - lag:n - lag, :)
      end do
      cross = matmul(transpose(design), design)
      call invert_matrix(cross, inverse, status)
      if (status /= 0) then
         out%info = 20 + status
         return
      end if
      out%degrees_of_freedom = size(targets)*size(out%predictors)*order
      allocate(out%restrictions(out%degrees_of_freedom))
      allocate(out%restriction_covariance(out%degrees_of_freedom, out%degrees_of_freedom))
      restriction = 0
      do i = 1, size(targets)
         target = targets(i)
         do lag = 1, order
            do j = 1, size(out%predictors)
               predictor = out%predictors(j)
               restriction = restriction + 1
               out%restrictions(restriction) = out%unrestricted%ar(target, predictor, lag)
            end do
         end do
      end do
      restriction_i = 0
      do i = 1, size(targets)
         do lag = 1, order
            do j = 1, size(out%predictors)
               restriction_i = restriction_i + 1
               first_regressor = offset + (lag - 1)*variables + out%predictors(j)
               restriction_j = 0
               do target = 1, size(targets)
                  do predictor = 1, order
                     do restriction = 1, size(out%predictors)
                        restriction_j = restriction_j + 1
                        out%restriction_covariance(restriction_i, restriction_j) = &
                           out%unrestricted%sigma(targets(i), targets(target))* &
                           inverse(first_regressor, offset + (predictor - 1)*variables + &
                           out%predictors(restriction))
                     end do
                  end do
               end do
            end do
         end do
      end do
      call invert_matrix(out%restriction_covariance, covariance_inverse, status)
      if (status /= 0) then
         out%info = 30 + status
         return
      end if
      out%statistic = dot_product(out%restrictions, matmul(covariance_inverse, out%restrictions))
      out%p_value = chi_square_survival(out%statistic, out%degrees_of_freedom)
   end function mts_granger_test

   pure function unique_indices(values) result(unique)
      !! Return the first occurrence of each integer value.
      integer, intent(in) :: values(:) !! Input values.
      integer, allocatable :: unique(:)
      integer :: work(size(values)), count_unique, i

      count_unique = 0
      do i = 1, size(values)
         if (count_unique == 0 .or. all(work(:count_unique) /= values(i))) then
            count_unique = count_unique + 1
            work(count_unique) = values(i)
         end if
      end do
      unique = work(:count_unique)
   end function unique_indices

   pure function mts_mq(series, max_lag, adjustment) result(out)
      !! Compute multivariate Ljung-Box statistics for all lags through max_lag.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      integer, intent(in) :: max_lag !! Maximum lag to consider.
      integer, intent(in), optional :: adjustment !! Adjustment.
      type(mts_mq_t) :: out
      real(dp), allocatable :: centered(:, :), covariance(:, :), inverse(:, :), lag_covariance(:, :), product(:, :)
      integer :: n, variables, lag, status, adjust, i

      n = size(series, 1)
      variables = size(series, 2)
      adjust = 0
      if (present(adjustment)) adjust = adjustment
      if (max_lag < 1 .or. max_lag >= n - 1 .or. adjust < 0) then
         out%info = 1
         return
      end if
      centered = series - spread(sum(series, 1)/real(n, dp), 1, n)
      covariance = matmul(transpose(centered), centered)/real(n - 1, dp)
      call invert_matrix(covariance, inverse, status)
      if (status /= 0) then
         out%info = 10 + status
         return
      end if
      allocate(out%statistic(max_lag), out%p_value(max_lag), out%degrees_of_freedom(max_lag))
      out%statistic = 0.0_dp
      do lag = 1, max_lag
         lag_covariance = matmul(transpose(centered(lag + 1:, :)), centered(:n - lag, :))/real(n - 1, dp)
         product = matmul(transpose(lag_covariance), matmul(inverse, matmul(lag_covariance, inverse)))
         if (lag > 1) out%statistic(lag) = out%statistic(lag - 1)
         out%statistic(lag) = out%statistic(lag) + real(n, dp)**2* &
            sum([(product(i, i), i=1, variables)])/real(n - lag, dp)
         out%degrees_of_freedom(lag) = lag*variables*variables - adjust
         out%p_value(lag) = 1.0_dp
         if (out%degrees_of_freedom(lag) > variables*variables - 1) then
            out%p_value(lag) = chi_square_survival(out%statistic(lag), out%degrees_of_freedom(lag))
         end if
      end do
      out%adjustment = adjust
   end function mts_mq

   pure function mts_diagnostic(residuals, max_lag, adjustment) result(out)
      !! Compute residual cross-correlations and multivariate Ljung-Box diagnostics.
      real(dp), intent(in) :: residuals(:, :) !! Model residuals.
      integer, intent(in) :: max_lag !! Maximum lag to consider.
      integer, intent(in), optional :: adjustment !! Adjustment.
      type(mts_diagnostic_t) :: out
      real(dp), allocatable :: centered(:, :), covariance(:, :), inverse(:, :), standardized(:, :)
      real(dp), allocatable :: vector(:), kron_inverse(:, :)
      integer :: n, variables, lag, status, i, j, row, column

      n = size(residuals, 1)
      variables = size(residuals, 2)
      if (max_lag < 1 .or. max_lag >= n - 1) then
         out%info = 1
         return
      end if
      centered = residuals - spread(sum(residuals, 1)/real(n, dp), 1, n)
      covariance = matmul(transpose(centered), centered)/real(n - 1, dp)
      call invert_matrix(covariance, inverse, status)
      if (status /= 0) then
         out%info = 10 + status
         return
      end if
      standardized = centered/spread(sqrt([(covariance(i, i), i=1, variables)]), 1, n)
      allocate(out%cross_correlation(variables, variables, 0:max_lag), out%p_value(max_lag))
      out%cross_correlation(:, :, 0) = matmul(transpose(standardized), standardized)/real(n - 1, dp)
      allocate(kron_inverse(variables*variables, variables*variables), vector(variables*variables))
      do column = 1, variables
         do row = 1, variables
            do j = 1, variables
               do i = 1, variables
                  kron_inverse(i + (j - 1)*variables, row + (column - 1)*variables) = &
                     inverse(i, row)*inverse(j, column)
               end do
            end do
         end do
      end do
      do lag = 1, max_lag
         out%cross_correlation(:, :, lag) = matmul(transpose(standardized(lag + 1:, :)), &
            standardized(:n - lag, :))/real(n, dp)
         vector = reshape(out%cross_correlation(:, :, lag), [variables*variables])
         out%p_value(lag) = chi_square_survival(real(n, dp)**2* &
            dot_product(vector, matmul(kron_inverse, vector))/real(n - lag, dp), variables*variables)
      end do
      out%mq = mts_mq(residuals, max_lag, adjustment)
      if (out%mq%info /= 0) out%info = 20 + out%mq%info
   end function mts_diagnostic

   pure function mts_var_backtest(series, order, origin, horizon, reestimate, include_mean) result(out)
      !! Run rolling-origin VAR forecasts with periodic model re-estimation.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      integer, intent(in) :: order !! Model or polynomial order.
      integer, intent(in) :: origin !! Origin.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      integer, intent(in), optional :: reestimate !! Reestimate.
      logical, intent(in), optional :: include_mean !! Whether to include a mean term.
      type(mts_var_backtest_t) :: out
      type(mts_var_fit_t) :: fitted
      type(mts_var_forecast_t) :: prediction
      integer :: n, variables, origins, cadence, origin_index, current, available, step, variable
      integer :: valid_count
      logical :: use_mean

      n = size(series, 1)
      variables = size(series, 2)
      cadence = 1
      use_mean = .true.
      if (present(reestimate)) cadence = reestimate
      if (present(include_mean)) use_mean = include_mean
      if (order < 1 .or. origin <= order + variables .or. origin >= n .or. &
         horizon < 1 .or. cadence < 1) then
         out%info = 1
         return
      end if
      origins = n - origin
      allocate(out%forecast(origins, horizon, variables), out%error(origins, horizon, variables))
      allocate(out%rmse(horizon, variables), out%mean_absolute_error(horizon, variables))
      out%forecast = 0.0_dp
      out%error = 0.0_dp
      do origin_index = 1, origins
         current = origin + origin_index - 1
         if (origin_index == 1 .or. mod(origin_index - 1, cadence) == 0) then
            fitted = mts_var(series(:current, :), order, use_mean)
            if (fitted%info /= 0) then
               out%info = 10 + origin_index
               return
            end if
         end if
         prediction = mts_var_forecast(fitted, series(:current, :), horizon)
         if (prediction%info /= 0) then
            out%info = 20 + origin_index
            return
         end if
         available = min(horizon, n - current)
         out%forecast(origin_index, :available, :) = prediction%mean(:available, :)
         out%error(origin_index, :available, :) = series(current + 1:current + available, :) - &
            prediction%mean(:available, :)
      end do
      do step = 1, horizon
         valid_count = origins - step + 1
         do variable = 1, variables
            out%rmse(step, variable) = sqrt(sum(out%error(:valid_count, step, variable)**2)/ &
               real(valid_count, dp))
            out%mean_absolute_error(step, variable) = &
               sum(abs(out%error(:valid_count, step, variable)))/real(valid_count, dp)
         end do
      end do
      out%origin = origin
      out%reestimate = cadence
   end function mts_var_backtest

   pure function mts_scm_identify(series, max_ar, max_ma, extra_lags, significance) result(out)
      !! Identify first-stage scalar components using corrected canonical correlations.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      integer, intent(in) :: max_ar !! Maximum autoregressive.
      integer, intent(in) :: max_ma !! Maximum moving-average.
      integer, intent(in), optional :: extra_lags !! Extra lags.
      real(dp), intent(in), optional :: significance !! Significance.
      type(mts_scm_identification_t) :: out
      real(dp), allocatable :: current_block(:, :), past_block(:, :)
      real(dp), allocatable :: squared_correlation(:), x_coefficients(:, :), y_coefficients(:, :)
      real(dp), allocatable :: x_component(:), y_component(:), variance_correction(:)
      real(dp) :: level, normal_cutoff, threshold, adjusted, statistic, p_value
      integer :: n, variables, h, ar_order, ma_order, start, rows, x_columns, y_columns
      integer :: block, canonical_count, tested_count, index, status, lag, count_zero

      n = size(series, 1)
      variables = size(series, 2)
      h = 0
      level = 0.05_dp
      if (present(extra_lags)) h = extra_lags
      if (present(significance)) level = significance
      if (max_ar < 0 .or. max_ma < 0 .or. h < 0 .or. variables < 1 .or. &
         level <= 0.0_dp .or. level >= 1.0_dp .or. &
         n <= max_ar + max_ma + h + 3) then
         out%info = 1
         return
      end if
      allocate(out%zero_count(0:max_ar, 0:max_ma), &
         out%diagonal_difference(0:max_ar, 0:max_ma))
      out%zero_count = 0
      normal_cutoff = normal_quantile(1.0_dp - 0.5_dp*level)
      threshold = normal_cutoff**2/real(n, dp)
      do ar_order = 0, max_ar
         do ma_order = 0, max_ma
            start = ar_order + ma_order + h + 2
            rows = n - start + 1
            x_columns = (ar_order + 1)*variables
            y_columns = (ar_order + h + 1)*variables
            allocate(current_block(rows, x_columns), past_block(rows, y_columns))
            do block = 0, ar_order
               current_block(:, block*variables + 1:(block + 1)*variables) = &
                  series(start - block:n - block, :)
            end do
            do block = 0, ar_order + h
               past_block(:, block*variables + 1:(block + 1)*variables) = &
                  series(start - ma_order - 1 - block:n - ma_order - 1 - block, :)
            end do
            call canonical_correlations(current_block, past_block, squared_correlation, &
               x_coefficients, y_coefficients, status)
            if (status /= 0) then
               out%info = 10 + ar_order*(max_ma + 1) + ma_order
               return
            end if
            canonical_count = size(squared_correlation)
            allocate(variance_correction(canonical_count))
            variance_correction = 1.0_dp
            if (ma_order > 0) then
               do index = 1, canonical_count
                  x_component = matmul(current_block, &
                     x_coefficients(:, index)/max(norm2(x_coefficients(:, index)), tiny(1.0_dp)))
                  y_component = matmul(past_block, &
                     y_coefficients(:, index)/max(norm2(y_coefficients(:, index)), tiny(1.0_dp)))
                  do lag = 1, ma_order
                     variance_correction(index) = variance_correction(index) + 2.0_dp* &
                        sample_autocorrelation(x_component, lag)*sample_autocorrelation(y_component, lag)
                  end do
               end do
            end if
            where (squared_correlation > threshold) variance_correction = 1.0_dp
            tested_count = min((ma_order + 1)*variables, canonical_count)
            statistic = 0.0_dp
            count_zero = 0
            do index = 1, tested_count
               adjusted = squared_correlation(canonical_count - index + 1)/ &
                  max(variance_correction(canonical_count - index + 1), tiny(1.0_dp))
               adjusted = min(0.999_dp, max(0.0_dp, adjusted))
               statistic = statistic - real(n - ar_order - ma_order, dp)*log(1.0_dp - adjusted)
               p_value = chi_square_survival(statistic, index*(index + h*variables))
               if (p_value >= level) count_zero = count_zero + 1
            end do
            out%zero_count(ar_order, ma_order) = count_zero
            deallocate(current_block, past_block, squared_correlation, x_coefficients, &
               y_coefficients, variance_correction)
         end do
      end do
      out%diagonal_difference = out%zero_count
      do ar_order = 1, max_ar
         do ma_order = 1, max_ma
            out%diagonal_difference(ar_order, ma_order) = min(variables, &
               out%zero_count(ar_order, ma_order) - out%zero_count(ar_order - 1, ma_order - 1))
         end do
      end do
      out%max_ar = max_ar
      out%max_ma = max_ma
      out%extra_lags = h
      out%significance = level
   end function mts_scm_identify

   pure function mts_scm_identify_details(series, max_ar, max_ma, extra_lags, &
      significance) result(out)
      !! Find independent scalar components and their individual ARMA orders.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      integer, intent(in) :: max_ar !! Maximum autoregressive.
      integer, intent(in) :: max_ma !! Maximum moving-average.
      integer, intent(in), optional :: extra_lags !! Extra lags.
      real(dp), intent(in), optional :: significance !! Significance.
      type(mts_scm_structure_t) :: out
      real(dp), allocatable :: current(:, :), past(:, :), correlations(:)
      real(dp), allocatable :: x_coefficients(:, :), y_coefficients(:, :), direction(:), basis(:, :)
      real(dp) :: level, rho, statistic, projection
      integer :: n, dimension, h, total, ar_order, ma_order, start, rows, lag, index, status, found

      n = size(series, 1)
      dimension = size(series, 2)
      h = 0
      level = 0.05_dp
      if (present(extra_lags)) h = extra_lags
      if (present(significance)) level = significance
      if (dimension < 1 .or. max_ar < 0 .or. max_ma < 0 .or. h < 0 .or. &
         n <= max_ar + max_ma + h + dimension + 2) then
         out%info = 1
         return
      end if
      allocate(out%order(dimension, 2), out%transformation(dimension, dimension))
      allocate(basis(dimension, dimension), direction(dimension))
      out%order = 0
      out%transformation = 0.0_dp
      basis = 0.0_dp
      found = 0
      do total = 0, max_ar + max_ma
         do ma_order = 0, min(max_ma, total)
            ar_order = total - ma_order
            if (ar_order > max_ar .or. found == dimension) cycle
            start = ar_order + ma_order + h + 2
            rows = n - start + 1
            allocate(current(rows, (ar_order + 1)*dimension))
            allocate(past(rows, (ar_order + h + 1)*dimension))
            do lag = 0, ar_order
               current(:, lag*dimension + 1:(lag + 1)*dimension) = &
                  series(start - lag:n - lag, :)
            end do
            do lag = 0, ar_order + h
               past(:, lag*dimension + 1:(lag + 1)*dimension) = &
                  series(start - ma_order - 1 - lag:n - ma_order - 1 - lag, :)
            end do
            call canonical_correlations(current, past, correlations, x_coefficients, &
               y_coefficients, status)
            if (status /= 0) then
               deallocate(current, past)
               cycle
            end if
            do index = size(correlations), 1, -1
               rho = min(1.0_dp - epsilon(1.0_dp), max(0.0_dp, correlations(index)))
               statistic = -real(n - ar_order - ma_order, dp)*log(1.0_dp - rho)
               if (chi_square_survival(statistic, max(1, size(correlations) - index + 1)) < level) cycle
               direction = x_coefficients(:dimension, index)
               do lag = 1, found
                  projection = dot_product(direction, basis(:, lag))
                  direction = direction - projection*basis(:, lag)
               end do
               if (norm2(direction) <= 100.0_dp*epsilon(1.0_dp)) cycle
               found = found + 1
               basis(:, found) = direction/norm2(direction)
               out%transformation(found, :) = basis(:, found)
               out%order(found, :) = [ar_order, ma_order]
               if (found == dimension) exit
            end do
            deallocate(current, past, correlations, x_coefficients, y_coefficients)
         end do
      end do
      out%components = found
      if (found < dimension) out%info = 2
   end function mts_scm_identify_details

   pure function mts_scm_specification(order, pivot) result(out)
      !! Build SCM parameter indicators for supplied component orders and pivots.
      integer, intent(in) :: order(:, :) !! Model or polynomial order.
      integer, intent(in) :: pivot(:) !! Pivot.
      type(mts_scm_spec_t) :: out
      integer :: dimension, row, column, lag, redundant

      dimension = size(order, 1)
      if (dimension < 1 .or. size(order, 2) /= 2 .or. size(pivot) /= dimension .or. &
         any(order < 0) .or. any(pivot < 1) .or. any(pivot > dimension)) then
         out%info = 1
         return
      end if
      if (any([(count(pivot == row), row=1, dimension)] /= 1)) then
         out%info = 2
         return
      end if
      out%ar_order = maxval(order(:, 1))
      out%ma_order = maxval(order(:, 2))
      allocate(out%order(dimension, 2), out%pivot(dimension))
      allocate(out%transformation(dimension, dimension))
      allocate(out%ar(dimension, dimension, out%ar_order))
      allocate(out%ma(dimension, dimension, out%ma_order))
      out%order = order
      out%pivot = pivot
      out%transformation = 2
      out%ar = 0
      out%ma = 0
      do row = 1, dimension
         out%transformation(row, pivot(row)) = 1
         do lag = 1, order(row, 1)
            out%ar(row, :, lag) = 2
         end do
         do lag = 1, order(row, 2)
            out%ma(row, :, lag) = 2
         end do
         do column = 1, row
            redundant = min(order(row, 1) - order(column, 1), &
               order(row, 2) - order(column, 2))
            do lag = 1, max(0, redundant)
               out%ma(row, column, lag) = 0
            end do
         end do
      end do
      do row = 2, dimension
         do column = 1, row - 1
            redundant = min(order(row, 1) - order(column, 1), &
               order(row, 2) - order(column, 2))
            if (redundant >= 0) out%transformation(row, pivot(column)) = 0
         end do
      end do
   end function mts_scm_specification

   pure function mts_scm_fit(series, order, pivot, include_mean, initial, estimated, &
      max_iterations, tolerance) result(out)
      !! Fit a scalar-component VARMA model with joint transformation estimates.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      integer, intent(in) :: order(:, :) !! Model or polynomial order.
      integer, intent(in) :: pivot(:) !! Pivot.
      logical, intent(in), optional :: include_mean !! Whether to include a mean term.
      logical, intent(in), optional :: estimated(:) !! Flag controlling estimated.
      real(dp), intent(in), optional :: initial(:) !! Initial value.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      type(mts_scm_fit_t) :: out
      type(mts_kronecker_spec_t) :: structural_specification
      integer :: dimension, overall_order

      dimension = size(series, 2)
      out%specification = mts_scm_specification(order, pivot)
      if (out%specification%info /= 0 .or. size(order, 1) /= dimension) then
         out%info = 1
         return
      end if
      overall_order = max(out%specification%ar_order, out%specification%ma_order)
      allocate(structural_specification%index(dimension))
      allocate(structural_specification%ar(dimension, dimension, 0:overall_order))
      allocate(structural_specification%ma(dimension, dimension, 0:overall_order))
      structural_specification%index = maxval(order, 2)
      structural_specification%ar = 0
      structural_specification%ma = 0
      structural_specification%ar(:, :, 0) = out%specification%transformation
      structural_specification%ma(:, :, 0) = out%specification%transformation
      if (out%specification%ar_order > 0) structural_specification%ar(:, :, &
         1:out%specification%ar_order) = out%specification%ar
      if (out%specification%ma_order > 0) structural_specification%ma(:, :, &
         1:out%specification%ma_order) = out%specification%ma
      structural_specification%order = overall_order
      out%fit = mts_kronecker_fit(series, structural_specification%index, include_mean, &
         initial, estimated, max_iterations, tolerance, structural_specification)
      out%info = out%fit%info
   end function mts_scm_fit

   pure function mts_scm_refine(series, fit, threshold, max_steps, &
      max_iterations, tolerance) result(out)
      !! Refine an SCM fit by backward elimination of insignificant parameters.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      type(mts_scm_fit_t), intent(in) :: fit !! Previously fitted model.
      real(dp), intent(in), optional :: threshold !! Decision or truncation threshold.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      integer, intent(in), optional :: max_steps !! Maximum steps.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      type(mts_scm_refinement_t) :: out
      type(mts_kronecker_refinement_t) :: refined

      if (fit%info /= 0) then
         out%info = 1
         return
      end if
      refined = mts_kronecker_refine(series, fit%fit, threshold, max_steps, &
         max_iterations, tolerance)
      out%fit%specification = fit%specification
      out%fit%fit = refined%fit
      out%fit%info = refined%info
      out%active_count = refined%active_count
      out%removed_index = refined%removed_index
      out%steps = refined%steps
      out%info = refined%info
   end function mts_scm_refine

   pure function mts_scm_forecast(fit, history, horizon) result(out)
      !! Forecast an SCM fit using its shared reduced-form VARMA representation.
      type(mts_scm_fit_t), intent(in) :: fit !! Previously fitted model.
      real(dp), intent(in) :: history(:, :) !! History.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      type(mts_var_forecast_t) :: out

      if (fit%info /= 0) then
         out%info = 1
         return
      end if
      out = mts_kronecker_forecast(fit%fit, history, horizon)
   end function mts_scm_forecast

   pure function mts_kronecker_specification(indices) result(out)
      !! Construct the echelon-form restrictions implied by Kronecker indices.
      integer, intent(in) :: indices(:) !! Indices.
      type(mts_kronecker_spec_t) :: out
      integer :: dimension, order, row, column, lag

      dimension = size(indices)
      if (dimension < 1 .or. any(indices < 0)) then
         out%info = 1
         return
      end if
      order = maxval(indices)
      allocate(out%index(dimension), out%ar(dimension, dimension, 0:order))
      allocate(out%ma(dimension, dimension, 0:order))
      out%index = indices
      out%ar = 2
      out%ma = 2
      do row = 1, dimension
         out%ar(row, row, 0) = 1
         out%ma(row, row, 0) = 1
         do lag = indices(row) + 1, order
            out%ar(row, :, lag) = 0
            out%ma(row, :, lag) = 0
         end do
         if (row < dimension) out%ar(row, row + 1:, 0) = 0
      end do
      do row = 2, dimension
         do column = 1, row - 1
            if (indices(column) <= indices(row)) out%ar(row, column, 0) = 0
         end do
      end do
      out%ma(:, :, 0) = out%ar(:, :, 0)
      do row = 1, dimension
         do column = 1, dimension
            do lag = 1, max(0, indices(row) - indices(column))
               out%ar(row, column, lag) = 0
            end do
         end do
      end do
      out%order = order
   end function mts_kronecker_specification

   pure function mts_kronecker_identify(series, past_lag, significance) result(out)
      !! Estimate Kronecker indices by sequential canonical-correlation tests.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      integer, intent(in), optional :: past_lag !! Past lag.
      real(dp), intent(in), optional :: significance !! Significance.
      type(mts_kronecker_identification_t) :: out
      real(dp), allocatable :: past(:, :), future(:, :), correlations(:)
      real(dp), allocatable :: x_coefficients(:, :), y_coefficients(:, :)
      real(dp) :: level, rho, statistic
      integer :: n, dimension, lag_count, height, variable, row, lag, rows, status, df

      n = size(series, 1)
      dimension = size(series, 2)
      lag_count = 5
      level = 0.05_dp
      if (present(past_lag)) lag_count = past_lag
      if (present(significance)) level = significance
      if (dimension < 1 .or. lag_count < 1 .or. n <= 2*lag_count + dimension + 2 .or. &
         level <= 0.0_dp .or. level >= 1.0_dp) then
         out%info = 1
         return
      end if
      allocate(out%index(dimension), out%statistic(dimension, 0:lag_count))
      allocate(out%p_value(dimension, 0:lag_count), out%tested(dimension, 0:lag_count))
      out%index = lag_count
      out%statistic = 0.0_dp
      out%p_value = 1.0_dp
      out%tested = .false.
      do height = 0, lag_count
         rows = n - lag_count - height
         allocate(past(rows, lag_count*dimension))
         do lag = 1, lag_count
            past(:, (lag - 1)*dimension + 1:lag*dimension) = &
               series(lag_count + 1 - lag:n - height - lag, :)
         end do
         do variable = 1, dimension
            if (height > 0 .and. out%index(variable) < height) cycle
            allocate(future(rows, height + 1))
            do lag = 0, height
               do row = 1, rows
                  future(row, lag + 1) = series(lag_count + row + lag, variable)
               end do
            end do
            call canonical_correlations(past, future, correlations, x_coefficients, &
               y_coefficients, status)
            if (status /= 0) then
               out%info = 10 + variable
               return
            end if
            rho = min(1.0_dp - epsilon(1.0_dp), max(0.0_dp, correlations(size(correlations))))
            df = max(1, size(past, 2) - size(future, 2) + 1)
            statistic = -real(rows - 1, dp)*log(1.0_dp - rho)
            out%statistic(variable, height) = statistic
            out%p_value(variable, height) = chi_square_survival(statistic, df)
            out%tested(variable, height) = .true.
            if (out%p_value(variable, height) > level) out%index(variable) = height
            deallocate(future, correlations, x_coefficients, y_coefficients)
         end do
         deallocate(past)
      end do
      out%past_lag = lag_count
      out%significance = level
   end function mts_kronecker_identify

   pure function mts_kronecker_fit(series, indices, include_mean, initial, estimated, &
      max_iterations, tolerance, specification) result(out)
      !! Fit an echelon-form VARMA model by conditional Gaussian likelihood.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      integer, intent(in) :: indices(:) !! Indices.
      logical, intent(in), optional :: include_mean !! Whether to include a mean term.
      logical, intent(in), optional :: estimated(:) !! Flag controlling estimated.
      real(dp), intent(in), optional :: initial(:) !! Initial value.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      type(mts_kronecker_spec_t), intent(in), optional :: specification !! Specification.
      type(mts_kronecker_fit_t) :: out
      type(optimization_result_t) :: optimization
      real(dp), allocatable :: parameters(:), free_values(:), hessian(:, :), inverse(:, :)
      integer, allocatable :: free_index(:)
      logical :: use_mean
      integer :: n, dimension, order, parameter_count, free_count, i, status, limit
      real(dp) :: gradient_tolerance

      n = size(series, 1)
      dimension = size(series, 2)
      use_mean = .true.
      if (present(include_mean)) use_mean = include_mean
      out%includes_mean = use_mean
      if (present(specification)) then
         out%specification = specification
      else
         out%specification = mts_kronecker_specification(indices)
      end if
      if (out%specification%info /= 0 .or. size(indices) /= dimension) then
         out%info = 1
         return
      end if
      order = out%specification%order
      if (n <= order + dimension + 2) then
         out%info = 1
         return
      end if
      parameter_count = merge(dimension, 0, use_mean) + count(out%specification%ar == 2) + &
         count(out%specification%ma(:, :, 1:) == 2)
      allocate(parameters(parameter_count), out%estimated(parameter_count))
      parameters = 0.0_dp
      if (use_mean) parameters(:dimension) = sum(series, 1)/real(n, dp)
      if (present(initial)) then
         if (size(initial) /= parameter_count) then
            out%info = 2
            return
         end if
         parameters = initial
      end if
      out%estimated = .true.
      if (present(estimated)) then
         if (size(estimated) /= parameter_count) then
            out%info = 2
            return
         end if
         out%estimated = estimated
      end if
      free_count = count(out%estimated)
      allocate(free_index(free_count), free_values(free_count))
      free_index = pack([(i, i=1, parameter_count)], out%estimated)
      free_values = pack(parameters, out%estimated)
      limit = 200
      gradient_tolerance = 1.0e-6_dp
      if (present(max_iterations)) limit = max_iterations
      if (present(tolerance)) gradient_tolerance = tolerance
      if (free_count > 0) then
         optimization = bfgs_minimize_fd(objective, free_values, limit, gradient_tolerance)
         parameters(free_index) = optimization%parameters
         out%iterations = optimization%iterations
         out%converged = optimization%converged
         if (optimization%info /= 0) out%info = 100 + optimization%info
      else
         out%converged = .true.
      end if
      call evaluate(parameters, out%structural_ar, out%structural_ma, out%model, &
         out%residuals, out%log_likelihood, status)
      if (status /= 0) then
         out%info = 10 + status
         return
      end if
      out%coefficients = parameters
      out%aic = -2.0_dp*out%log_likelihood + 2.0_dp*real(free_count + 1, dp)
      out%bic = -2.0_dp*out%log_likelihood + log(real(n - order, dp))*real(free_count + 1, dp)
      allocate(out%covariance(parameter_count, parameter_count), out%standard_errors(parameter_count))
      out%covariance = 0.0_dp
      out%standard_errors = 0.0_dp
      if (free_count > 0) then
         free_values = pack(parameters, out%estimated)
         hessian = finite_difference_hessian(objective, free_values)
         call invert_matrix(hessian, inverse, status)
         if (status == 0) then
            do i = 1, free_count
               out%covariance(free_index, free_index(i)) = inverse(:, i)
               out%standard_errors(free_index(i)) = sqrt(max(0.0_dp, inverse(i, i)))
            end do
         else if (out%info == 0) then
            out%info = 200 + status
         end if
      end if

   contains

      pure function objective(free_parameters) result(value)
         !! Evaluate the negative conditional log likelihood for free parameters.
         real(dp), intent(in) :: free_parameters(:) !! Free parameters.
         real(dp) :: value, likelihood
         real(dp) :: full_parameters(parameter_count)
         real(dp), allocatable :: structural_ar(:, :, :), structural_ma(:, :, :), residual_values(:, :)
         type(mts_varma_model_t) :: candidate
         integer :: evaluation_status

         full_parameters = parameters
         full_parameters(free_index) = free_parameters
         call evaluate(full_parameters, structural_ar, structural_ma, candidate, residual_values, &
            likelihood, evaluation_status)
         if (evaluation_status == 0 .and. ieee_is_finite(likelihood)) then
            value = -likelihood
         else
            value = 1.0e30_dp + dot_product(free_parameters, free_parameters)
         end if
      end function objective

      pure subroutine evaluate(values, structural_ar, structural_ma, model, residual_values, &
         likelihood, evaluation_status)
         !! Map structural parameters to reduced form and evaluate residuals.
         real(dp), intent(in) :: values(:) !! Input values.
         real(dp), allocatable, intent(out) :: structural_ar(:, :, :) !! Structural autoregressive.
         real(dp), allocatable, intent(out) :: structural_ma(:, :, :) !! Structural moving-average.
         type(mts_varma_model_t), intent(out) :: model !! Model specification.
         real(dp), allocatable, intent(out) :: residual_values(:, :) !! Residual values.
         real(dp), intent(out) :: likelihood !! Likelihood.
         integer, intent(out) :: evaluation_status !! Evaluation status.
         real(dp), allocatable :: lag_zero_inverse(:, :), covariance_inverse(:, :)
         real(dp) :: logdet
         integer :: offset, row, column, lag, t, effective

         allocate(structural_ar(dimension, dimension, 0:order))
         allocate(structural_ma(dimension, dimension, 0:order))
         structural_ar = 0.0_dp
         structural_ma = 0.0_dp
         offset = merge(dimension, 0, use_mean)
         do lag = 0, order
            do column = 1, dimension
               do row = 1, dimension
                  if (out%specification%ar(row, column, lag) == 1) then
                     structural_ar(row, column, lag) = 1.0_dp
                  else if (out%specification%ar(row, column, lag) == 2) then
                     offset = offset + 1
                     structural_ar(row, column, lag) = values(offset)
                  end if
               end do
            end do
         end do
         structural_ma(:, :, 0) = structural_ar(:, :, 0)
         do lag = 1, order
            do column = 1, dimension
               do row = 1, dimension
                  if (out%specification%ma(row, column, lag) == 2) then
                     offset = offset + 1
                     structural_ma(row, column, lag) = values(offset)
                  end if
               end do
            end do
         end do
         call invert_matrix(structural_ar(:, :, 0), lag_zero_inverse, evaluation_status)
         if (evaluation_status /= 0) then
            likelihood = -huge(1.0_dp)
            return
         end if
         allocate(model%ar(dimension, dimension, order), model%ma(dimension, dimension, order))
         allocate(model%intercept(dimension), model%sigma(dimension, dimension))
         model%intercept = 0.0_dp
         if (use_mean) model%intercept = matmul(lag_zero_inverse, values(:dimension))
         do lag = 1, order
            model%ar(:, :, lag) = matmul(lag_zero_inverse, structural_ar(:, :, lag))
            model%ma(:, :, lag) = -matmul(lag_zero_inverse, structural_ma(:, :, lag))
         end do
         allocate(residual_values(n, dimension))
         residual_values = 0.0_dp
         do t = 1, n
            residual_values(t, :) = series(t, :) - model%intercept
            do lag = 1, min(order, t - 1)
               residual_values(t, :) = residual_values(t, :) - &
                  matmul(model%ar(:, :, lag), series(t - lag, :)) + &
                  matmul(model%ma(:, :, lag), residual_values(t - lag, :))
            end do
         end do
         effective = n - order
         model%sigma = matmul(transpose(residual_values(order + 1:, :)), &
            residual_values(order + 1:, :))/real(effective, dp)
         allocate(covariance_inverse(dimension, dimension))
         call inverse_logdet(model%sigma, covariance_inverse, logdet, evaluation_status, &
            100.0_dp*epsilon(1.0_dp))
         if (evaluation_status /= 0) then
            likelihood = -huge(1.0_dp)
            return
         end if
         likelihood = -0.5_dp*real(effective, dp)*(real(dimension, dp)* &
            (log(2.0_dp*acos(-1.0_dp)) + 1.0_dp) + logdet)
      end subroutine evaluate
   end function mts_kronecker_fit

   pure function mts_kronecker_refine(series, fit, threshold, max_steps, &
      max_iterations, tolerance) result(out)
      !! Remove individually insignificant Kronecker parameters and refit.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      type(mts_kronecker_fit_t), intent(in) :: fit !! Previously fitted model.
      real(dp), intent(in), optional :: threshold !! Decision or truncation threshold.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      integer, intent(in), optional :: max_steps !! Maximum steps.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      type(mts_kronecker_refinement_t) :: out
      logical, allocatable :: mask(:)
      real(dp) :: cutoff
      integer :: limit, iteration_limit, step, weakest

      cutoff = 1.96_dp
      limit = size(fit%coefficients)
      iteration_limit = 200
      if (present(threshold)) cutoff = threshold
      if (present(max_steps)) limit = min(limit, max_steps)
      if (present(max_iterations)) iteration_limit = max_iterations
      if (fit%info /= 0 .or. cutoff < 0.0_dp) then
         out%info = 1
         return
      end if
      mask = fit%estimated
      out%fit = fit
      allocate(out%active_count(0:limit), out%removed_index(limit))
      out%active_count(0) = count(mask)
      do step = 1, limit
         weakest = weakest_parameter(out%fit%coefficients, out%fit%standard_errors, mask, cutoff)
         if (weakest == 0) exit
         mask(weakest) = .false.
         out%fit = mts_kronecker_fit(series, fit%specification%index, include_mean=fit%includes_mean, &
            initial=out%fit%coefficients, estimated=mask, max_iterations=iteration_limit, &
            tolerance=tolerance, specification=fit%specification)
         if (out%fit%info /= 0 .and. out%fit%info < 100) then
            out%info = 10 + step
            return
         end if
         out%removed_index(step) = weakest
         out%active_count(step) = count(mask)
         out%steps = step
      end do
      out%active_count = out%active_count(:out%steps)
      out%removed_index = out%removed_index(:out%steps)
   end function mts_kronecker_refine

   pure function mts_kronecker_forecast(fit, history, horizon) result(out)
      !! Forecast a fitted Kronecker-index model through its reduced VARMA form.
      type(mts_kronecker_fit_t), intent(in) :: fit !! Previously fitted model.
      real(dp), intent(in) :: history(:, :) !! History.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      type(mts_var_forecast_t) :: out

      if (fit%info /= 0 .or. size(fit%residuals, 1) /= size(history, 1)) then
         out%info = 1
         return
      end if
      out = mts_varma_forecast(fit%model, history, fit%residuals, horizon)
   end function mts_kronecker_forecast

   pure function weakest_parameter(coefficients, standard_errors, mask, threshold) result(index)
      !! Return the active parameter with the smallest sub-threshold t ratio.
      real(dp), intent(in) :: coefficients(:) !! Model coefficients.
      real(dp), intent(in) :: standard_errors(:) !! Standard errors.
      real(dp), intent(in) :: threshold !! Decision or truncation threshold.
      logical, intent(in) :: mask(:) !! Flag controlling mask.
      integer :: index
      real(dp) :: ratio, smallest
      integer :: i

      index = 0
      smallest = huge(1.0_dp)
      do i = 1, size(coefficients)
         if (.not. mask(i) .or. standard_errors(i) <= 0.0_dp) cycle
         ratio = abs(coefficients(i)/standard_errors(i))
         if (ratio < threshold .and. ratio < smallest) then
            smallest = ratio
            index = i
         end if
      end do
   end function weakest_parameter

   pure subroutine canonical_correlations(x, y, squared_correlation, x_coefficients, &
      y_coefficients, info)
      !! Compute squared canonical correlations and coefficient vectors.
      real(dp), intent(in) :: x(:, :) !! Input data or predictor values.
      real(dp), intent(in) :: y(:, :) !! Response or time-series observations.
      real(dp), allocatable, intent(out) :: squared_correlation(:) !! Squared correlation.
      real(dp), allocatable, intent(out) :: x_coefficients(:, :) !! X coefficients.
      real(dp), allocatable, intent(out) :: y_coefficients(:, :) !! Y coefficients.
      integer, intent(out) :: info !! Status code; zero indicates success.
      real(dp), allocatable :: centered_x(:, :), centered_y(:, :)
      real(dp), allocatable :: sxx(:, :), syy(:, :), sxy(:, :), wx(:, :), wy(:, :)
      real(dp), allocatable :: kernel(:, :), eigenvalues(:), eigenvectors(:, :), syy_inverse(:, :)
      integer :: rows, count_canonical, index, status

      rows = size(x, 1)
      centered_x = x - spread(sum(x, 1)/real(rows, dp), 1, rows)
      centered_y = y - spread(sum(y, 1)/real(rows, dp), 1, rows)
      sxx = matmul(transpose(centered_x), centered_x)/real(rows - 1, dp)
      syy = matmul(transpose(centered_y), centered_y)/real(rows - 1, dp)
      sxy = matmul(transpose(centered_x), centered_y)/real(rows - 1, dp)
      call symmetric_inverse_sqrt(sxx, wx, status)
      if (status /= 0) then
         info = 1
         return
      end if
      call symmetric_inverse_sqrt(syy, wy, status)
      if (status /= 0) then
         info = 2
         return
      end if
      kernel = matmul(wx, matmul(sxy, matmul(matmul(wy, wy), matmul(transpose(sxy), wx))))
      kernel = 0.5_dp*(kernel + transpose(kernel))
      call symmetric_eigen(kernel, eigenvalues, eigenvectors, status)
      if (status /= 0) then
         info = 3
         return
      end if
      count_canonical = min(size(x, 2), size(y, 2))
      squared_correlation = max(0.0_dp, min(1.0_dp, eigenvalues(:count_canonical)))
      x_coefficients = matmul(wx, eigenvectors(:, :count_canonical))
      syy_inverse = matmul(wy, wy)
      allocate(y_coefficients(size(y, 2), count_canonical))
      do index = 1, count_canonical
         if (squared_correlation(index) > sqrt(epsilon(1.0_dp))) then
            y_coefficients(:, index) = matmul(syy_inverse, &
               matmul(transpose(sxy), x_coefficients(:, index)))/sqrt(squared_correlation(index))
         else
            y_coefficients(:, index) = 0.0_dp
         end if
      end do
      info = 0
   end subroutine canonical_correlations

   pure subroutine symmetric_inverse_sqrt(matrix, inverse_sqrt, info)
      !! Compute a symmetric positive-definite inverse square root.
      real(dp), intent(in) :: matrix(:, :) !! Input matrix.
      real(dp), allocatable, intent(out) :: inverse_sqrt(:, :) !! Inverse sqrt.
      integer, intent(out) :: info !! Status code; zero indicates success.
      real(dp), allocatable :: values(:), vectors(:, :), diagonal(:, :)
      real(dp) :: tolerance
      integer :: dimension, index

      dimension = size(matrix, 1)
      call symmetric_eigen(matrix, values, vectors, info)
      if (info /= 0) return
      tolerance = 100.0_dp*epsilon(1.0_dp)*max(1.0_dp, maxval(values))
      if (minval(values) <= tolerance) then
         info = 1
         return
      end if
      allocate(diagonal(dimension, dimension))
      diagonal = 0.0_dp
      do index = 1, dimension
         diagonal(index, index) = 1.0_dp/sqrt(values(index))
      end do
      inverse_sqrt = matmul(vectors, matmul(diagonal, transpose(vectors)))
   end subroutine symmetric_inverse_sqrt

   pure real(dp) function sample_autocorrelation(values, lag) result(correlation)
      !! Return the conventional sample autocorrelation at one lag.
      real(dp), intent(in) :: values(:) !! Input values.
      integer, intent(in) :: lag !! Lag index or number of lags.
      real(dp) :: mean_value, denominator

      mean_value = sum(values)/real(size(values), dp)
      denominator = sum((values - mean_value)**2)
      if (denominator <= tiny(1.0_dp)) then
         correlation = 0.0_dp
      else
         correlation = sum((values(lag + 1:) - mean_value)* &
            (values(:size(values) - lag) - mean_value))/denominator
      end if
   end function sample_autocorrelation

   pure subroutine arch_f_test(series, lag_count, statistic, p_value, info)
      !! Test squared observations on an intercept and their own lags.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      integer, intent(in) :: lag_count !! Number of lag.
      real(dp), intent(out) :: statistic !! Statistic.
      real(dp), intent(out) :: p_value !! P value.
      integer, intent(out) :: info !! Status code; zero indicates success.
      real(dp), allocatable :: design(:, :), response(:), xtx(:, :), inverse(:, :), beta(:), residual(:)
      real(dp) :: centered_sum, residual_sum, numerator_df, denominator_df, beta_argument
      integer :: n, rows, lag, status

      n = size(series)
      rows = n - lag_count
      if (lag_count < 1 .or. rows <= lag_count + 1) then
         info = 1
         return
      end if
      allocate(design(rows, lag_count + 1), response(rows))
      design(:, 1) = 1.0_dp
      response = series(lag_count + 1:)**2
      do lag = 1, lag_count
         design(:, lag + 1) = series(lag_count + 1 - lag:n - lag)**2
      end do
      xtx = matmul(transpose(design), design)
      call invert_matrix(xtx, inverse, status)
      if (status /= 0) then
         info = 10 + status
         return
      end if
      beta = matmul(inverse, matmul(transpose(design), response))
      residual = response - matmul(design, beta)
      residual_sum = sum(residual**2)
      centered_sum = sum((response - sum(response)/real(rows, dp))**2)
      numerator_df = real(lag_count, dp)
      denominator_df = real(rows - lag_count - 1, dp)
      if (residual_sum <= tiny(1.0_dp) .or. denominator_df <= 0.0_dp) then
         info = 2
         return
      end if
      statistic = max(0.0_dp, ((centered_sum - residual_sum)/numerator_df)/ &
         (residual_sum/denominator_df))
      beta_argument = numerator_df*statistic/(numerator_df*statistic + denominator_df)
      p_value = 1.0_dp - regularized_beta_mts(beta_argument, &
         0.5_dp*numerator_df, 0.5_dp*denominator_df)
      p_value = max(0.0_dp, min(1.0_dp, p_value))
      info = 0
   end subroutine arch_f_test

   pure real(dp) function regularized_beta_mts(value, first_shape, second_shape) result(probability)
      !! Return the regularized incomplete beta function.
      real(dp), intent(in) :: value !! Input value.
      real(dp), intent(in) :: first_shape !! First shape.
      real(dp), intent(in) :: second_shape !! Second shape.
      real(dp) :: factor

      if (value <= 0.0_dp) then
         probability = 0.0_dp
         return
      end if
      if (value >= 1.0_dp) then
         probability = 1.0_dp
         return
      end if
      factor = exp(log_gamma(first_shape + second_shape) - log_gamma(first_shape) - &
         log_gamma(second_shape) + first_shape*log(value) + second_shape*log(1.0_dp - value))
      if (value < (first_shape + 1.0_dp)/(first_shape + second_shape + 2.0_dp)) then
         probability = factor*beta_fraction_mts(value, first_shape, second_shape)/first_shape
      else
         probability = 1.0_dp - factor*beta_fraction_mts(1.0_dp - value, &
            second_shape, first_shape)/second_shape
      end if
   end function regularized_beta_mts

   pure real(dp) function beta_fraction_mts(value, first_shape, second_shape) result(fraction)
      !! Evaluate the incomplete-beta continued fraction.
      real(dp), intent(in) :: value !! Input value.
      real(dp), intent(in) :: first_shape !! First shape.
      real(dp), intent(in) :: second_shape !! Second shape.
      real(dp) :: qab, qap, qam, c, d, h, aa, delta
      integer :: iteration, twice

      qab = first_shape + second_shape
      qap = first_shape + 1.0_dp
      qam = first_shape - 1.0_dp
      c = 1.0_dp
      d = 1.0_dp - qab*value/qap
      if (abs(d) < tiny(1.0_dp)) d = tiny(1.0_dp)
      d = 1.0_dp/d
      h = d
      do iteration = 1, 10000
         twice = 2*iteration
         aa = real(iteration, dp)*(second_shape - real(iteration, dp))*value/ &
            ((qam + real(twice, dp))*(first_shape + real(twice, dp)))
         d = 1.0_dp + aa*d
         if (abs(d) < tiny(1.0_dp)) d = tiny(1.0_dp)
         c = 1.0_dp + aa/c
         if (abs(c) < tiny(1.0_dp)) c = tiny(1.0_dp)
         d = 1.0_dp/d
         h = h*d*c
         aa = -(first_shape + real(iteration, dp))*(qab + real(iteration, dp))*value/ &
            ((first_shape + real(twice, dp))*(qap + real(twice, dp)))
         d = 1.0_dp + aa*d
         if (abs(d) < tiny(1.0_dp)) d = tiny(1.0_dp)
         c = 1.0_dp + aa/c
         if (abs(c) < tiny(1.0_dp)) c = tiny(1.0_dp)
         d = 1.0_dp/d
         delta = d*c
         h = h*delta
         if (abs(delta - 1.0_dp) <= 10.0_dp*epsilon(1.0_dp)) exit
      end do
      fraction = h
   end function beta_fraction_mts

   pure function mts_svarma_expand(regular_ar, seasonal_ar, regular_ma, seasonal_ma, &
      intercept, sigma, period, switched) result(out)
      !! Expand multiplicative seasonal matrix polynomials into ordinary VARMA lags.
      real(dp), intent(in) :: regular_ar(:, :, :) !! Regular autoregressive.
      real(dp), intent(in) :: seasonal_ar(:, :, :) !! Seasonal autoregressive.
      real(dp), intent(in) :: regular_ma(:, :, :) !! Regular moving-average.
      real(dp), intent(in) :: seasonal_ma(:, :, :) !! Seasonal moving-average.
      real(dp), intent(in) :: intercept(:) !! Model intercept.
      real(dp), intent(in) :: sigma(:, :) !! Scale parameter or standard deviation.
      integer, intent(in) :: period !! Seasonal period.
      logical, intent(in), optional :: switched !! Flag controlling switched.
      type(mts_svarma_model_t) :: out
      integer :: dimension, p, q, seasonal_p, seasonal_q, ar_order, ma_order, i, j, lag
      logical :: reverse_products

      dimension = size(intercept)
      p = size(regular_ar, 3)
      seasonal_p = size(seasonal_ar, 3)
      q = size(regular_ma, 3)
      seasonal_q = size(seasonal_ma, 3)
      reverse_products = .false.
      if (present(switched)) reverse_products = switched
      if (period < 1 .or. size(sigma, 1) /= dimension .or. size(sigma, 2) /= dimension .or. &
         size(regular_ar, 1) /= dimension .or. size(regular_ar, 2) /= dimension .or. &
         size(seasonal_ar, 1) /= dimension .or. size(seasonal_ar, 2) /= dimension .or. &
         size(regular_ma, 1) /= dimension .or. size(regular_ma, 2) /= dimension .or. &
         size(seasonal_ma, 1) /= dimension .or. size(seasonal_ma, 2) /= dimension) then
         out%info = 1
         return
      end if
      out%regular_ar = regular_ar
      out%seasonal_ar = seasonal_ar
      out%regular_ma = regular_ma
      out%seasonal_ma = seasonal_ma
      out%intercept = intercept
      out%sigma = sigma
      out%period = period
      out%switched = reverse_products
      ar_order = max(p, seasonal_p*period)
      if (p > 0 .and. seasonal_p > 0) ar_order = p + seasonal_p*period
      ma_order = max(q, seasonal_q*period)
      if (q > 0 .and. seasonal_q > 0) ma_order = q + seasonal_q*period
      allocate(out%expanded%ar(dimension, dimension, ar_order))
      allocate(out%expanded%ma(dimension, dimension, ma_order))
      allocate(out%expanded%intercept(dimension), out%expanded%sigma(dimension, dimension))
      out%expanded%ar = 0.0_dp
      out%expanded%ma = 0.0_dp
      out%expanded%intercept = intercept
      out%expanded%sigma = sigma
      do i = 1, p
         out%expanded%ar(:, :, i) = out%expanded%ar(:, :, i) + regular_ar(:, :, i)
      end do
      do j = 1, seasonal_p
         lag = j*period
         out%expanded%ar(:, :, lag) = out%expanded%ar(:, :, lag) + seasonal_ar(:, :, j)
         do i = 1, p
            if (reverse_products) then
               out%expanded%ar(:, :, lag + i) = out%expanded%ar(:, :, lag + i) - &
                  matmul(seasonal_ar(:, :, j), regular_ar(:, :, i))
            else
               out%expanded%ar(:, :, lag + i) = out%expanded%ar(:, :, lag + i) - &
                  matmul(regular_ar(:, :, i), seasonal_ar(:, :, j))
            end if
         end do
      end do
      do i = 1, q
         out%expanded%ma(:, :, i) = out%expanded%ma(:, :, i) + regular_ma(:, :, i)
      end do
      do j = 1, seasonal_q
         lag = j*period
         out%expanded%ma(:, :, lag) = out%expanded%ma(:, :, lag) + seasonal_ma(:, :, j)
         do i = 1, q
            if (reverse_products) then
               out%expanded%ma(:, :, lag + i) = out%expanded%ma(:, :, lag + i) - &
                  matmul(seasonal_ma(:, :, j), regular_ma(:, :, i))
            else
               out%expanded%ma(:, :, lag + i) = out%expanded%ma(:, :, lag + i) - &
                  matmul(regular_ma(:, :, i), seasonal_ma(:, :, j))
            end if
         end do
      end do
   end function mts_svarma_expand

   pure function mts_svarma_fit(series, regular_ar_order, regular_ma_order, &
      seasonal_ar_order, seasonal_ma_order, period, include_mean, switched, &
      initial, estimated, max_iterations, tolerance) result(out)
      !! Estimate a multiplicative seasonal VARMA model by conditional likelihood.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      integer, intent(in) :: regular_ar_order !! Regular autoregressive order.
      integer, intent(in) :: regular_ma_order !! Regular moving-average order.
      integer, intent(in) :: seasonal_ar_order !! Seasonal autoregressive order.
      integer, intent(in) :: seasonal_ma_order !! Seasonal moving-average order.
      integer, intent(in) :: period !! Seasonal period.
      logical, intent(in), optional :: include_mean !! Whether to include a mean term.
      logical, intent(in), optional :: switched !! Flag controlling switched.
      logical, intent(in), optional :: estimated(:) !! Flag controlling estimated.
      real(dp), intent(in), optional :: initial(:) !! Initial value.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      type(mts_svarma_fit_t) :: out
      type(optimization_result_t) :: optimization
      type(mts_var_fit_t) :: var_start
      real(dp), allocatable :: parameters(:), free_parameters(:), hessian(:, :), inverse(:, :)
      integer, allocatable :: free_index(:)
      integer :: n, dimension, block_size, component_count, parameter_count
      integer :: free_count, offset, lag, i, status, limit
      real(dp) :: gradient_tolerance
      logical :: use_mean, reverse_products

      n = size(series, 1)
      dimension = size(series, 2)
      use_mean = .true.
      reverse_products = .false.
      if (present(include_mean)) use_mean = include_mean
      if (present(switched)) reverse_products = switched
      component_count = regular_ar_order + seasonal_ar_order + &
         regular_ma_order + seasonal_ma_order
      block_size = dimension*dimension
      parameter_count = component_count*block_size + merge(dimension, 0, use_mean)
      if (min(regular_ar_order, regular_ma_order, seasonal_ar_order, seasonal_ma_order) < 0 .or. &
         component_count < 1 .or. period < 2 .or. dimension < 1 .or. &
         n <= max(regular_ar_order + seasonal_ar_order*period, &
         regular_ma_order + seasonal_ma_order*period) + 2) then
         out%info = 1
         return
      end if
      allocate(parameters(parameter_count), out%estimated(parameter_count))
      parameters = 0.0_dp
      offset = 0
      if (use_mean) then
         parameters(:dimension) = sum(series, 1)/real(n, dp)
         offset = dimension
      end if
      if (regular_ar_order > 0) then
         var_start = mts_var(series, regular_ar_order, use_mean)
         if (var_start%info == 0) then
            if (use_mean) parameters(:dimension) = var_start%intercept
            do lag = 1, regular_ar_order
               parameters(offset + 1:offset + block_size) = &
                  reshape(var_start%ar(:, :, lag), [block_size])
               offset = offset + block_size
            end do
         end if
      end if
      if (present(initial)) then
         if (size(initial) /= parameter_count) then
            out%info = 2
            return
         end if
         parameters = initial
      end if
      out%estimated = .true.
      if (present(estimated)) then
         if (size(estimated) /= parameter_count) then
            out%info = 2
            return
         end if
         out%estimated = estimated
      end if
      free_count = sum(merge(1, 0, out%estimated))
      allocate(free_index(free_count), free_parameters(free_count))
      free_index = pack([(i, i=1, parameter_count)], out%estimated)
      free_parameters = pack(parameters, out%estimated)
      limit = 200
      gradient_tolerance = 1.0e-6_dp
      if (present(max_iterations)) limit = max_iterations
      if (present(tolerance)) gradient_tolerance = tolerance
      if (free_count > 0) then
         optimization = bfgs_minimize_fd(objective, free_parameters, limit, gradient_tolerance)
         parameters(free_index) = optimization%parameters
         out%iterations = optimization%iterations
         out%converged = optimization%converged
         if (optimization%info /= 0) out%info = 100 + optimization%info
      else
         out%converged = .true.
      end if
      call evaluate(parameters, out%model, out%residuals, out%log_likelihood, status)
      if (status /= 0) then
         out%info = 10 + status
         return
      end if
      out%coefficients = parameters
      out%aic = -2.0_dp*out%log_likelihood + 2.0_dp*real(free_count + 1, dp)
      out%bic = -2.0_dp*out%log_likelihood + log(real(n, dp))*real(free_count + 1, dp)
      out%stationarity_radius = companion_radius(out%model%expanded%ar)
      out%invertibility_radius = companion_radius(out%model%expanded%ma)
      allocate(out%covariance(parameter_count, parameter_count), out%standard_errors(parameter_count))
      out%covariance = 0.0_dp
      out%standard_errors = 0.0_dp
      if (free_count > 0) then
         free_parameters = pack(parameters, out%estimated)
         hessian = finite_difference_hessian(objective, free_parameters)
         call invert_matrix(hessian, inverse, status)
         if (status == 0) then
            do i = 1, free_count
               out%covariance(free_index, free_index(i)) = inverse(:, i)
               out%standard_errors(free_index(i)) = sqrt(max(0.0_dp, inverse(i, i)))
            end do
         else if (out%info == 0) then
            out%info = 200 + status
         end if
      end if

   contains

      pure function objective(free_values) result(value)
         !! Return the seasonal VARMA negative log likelihood.
         real(dp), intent(in) :: free_values(:) !! Free values.
         real(dp) :: value, log_likelihood
         real(dp) :: full_values(parameter_count)
         real(dp), allocatable :: residual_values(:, :)
         type(mts_svarma_model_t) :: candidate
         integer :: evaluation_status

         full_values = parameters
         full_values(free_index) = free_values
         call evaluate(full_values, candidate, residual_values, log_likelihood, evaluation_status)
         if (evaluation_status == 0 .and. ieee_is_finite(log_likelihood)) then
            value = -log_likelihood
         else
            value = 1.0e30_dp + dot_product(free_values, free_values)
         end if
      end function objective

      pure subroutine evaluate(values, model, residual_values, log_likelihood, evaluation_status)
         !! Expand component parameters and evaluate conditional residuals.
         real(dp), intent(in) :: values(:) !! Input values.
         type(mts_svarma_model_t), intent(out) :: model !! Model specification.
         real(dp), allocatable, intent(out) :: residual_values(:, :) !! Residual values.
         real(dp), intent(out) :: log_likelihood !! Log-likelihood value.
         integer, intent(out) :: evaluation_status !! Evaluation status.
         real(dp), allocatable :: regular_ar(:, :, :), seasonal_ar(:, :, :)
         real(dp), allocatable :: regular_ma(:, :, :), seasonal_ma(:, :, :)
         real(dp), allocatable :: intercept(:), sigma(:, :), covariance_inverse(:, :)
         real(dp) :: logdet
         integer :: value_offset, current_lag, t, p_expanded, q_expanded

         allocate(regular_ar(dimension, dimension, regular_ar_order))
         allocate(seasonal_ar(dimension, dimension, seasonal_ar_order))
         allocate(regular_ma(dimension, dimension, regular_ma_order))
         allocate(seasonal_ma(dimension, dimension, seasonal_ma_order))
         allocate(intercept(dimension), sigma(dimension, dimension))
         intercept = 0.0_dp
         value_offset = 0
         if (use_mean) then
            intercept = values(:dimension)
            value_offset = dimension
         end if
         call unpack_component(values, value_offset, regular_ar)
         call unpack_component(values, value_offset, seasonal_ar)
         call unpack_component(values, value_offset, regular_ma)
         call unpack_component(values, value_offset, seasonal_ma)
         sigma = 0.0_dp
         model = mts_svarma_expand(regular_ar, seasonal_ar, regular_ma, seasonal_ma, &
            intercept, sigma, period, reverse_products)
         p_expanded = size(model%expanded%ar, 3)
         q_expanded = size(model%expanded%ma, 3)
         allocate(residual_values(n, dimension))
         do t = 1, n
            residual_values(t, :) = series(t, :) - intercept
            do current_lag = 1, min(p_expanded, t - 1)
               residual_values(t, :) = residual_values(t, :) - &
                  matmul(model%expanded%ar(:, :, current_lag), series(t - current_lag, :))
            end do
            do current_lag = 1, min(q_expanded, t - 1)
               residual_values(t, :) = residual_values(t, :) + &
                  matmul(model%expanded%ma(:, :, current_lag), residual_values(t - current_lag, :))
            end do
         end do
         model%sigma = matmul(transpose(residual_values), residual_values)/real(n, dp)
         model%expanded%sigma = model%sigma
         allocate(covariance_inverse(dimension, dimension))
         call inverse_logdet(model%sigma, covariance_inverse, logdet, evaluation_status, &
            100.0_dp*epsilon(1.0_dp))
         if (evaluation_status /= 0) then
            log_likelihood = -huge(1.0_dp)
            return
         end if
         log_likelihood = -0.5_dp*real(n, dp)*(real(dimension, dp)*(log(2.0_dp*acos(-1.0_dp)) + &
            1.0_dp) + logdet)
      end subroutine evaluate

      pure subroutine unpack_component(values, value_offset, component)
         !! Read one lag-major component block from the parameter vector.
         real(dp), intent(in) :: values(:) !! Input values.
         integer, intent(inout) :: value_offset !! Value offset, updated in place.
         real(dp), intent(out) :: component(:, :, :) !! Component.
         integer :: current_lag

         do current_lag = 1, size(component, 3)
            component(:, :, current_lag) = reshape(&
               values(value_offset + 1:value_offset + block_size), [dimension, dimension])
            value_offset = value_offset + block_size
         end do
      end subroutine unpack_component
   end function mts_svarma_fit

   pure function mts_varma_fit(series, ar_order, ma_order, include_mean, initial, estimated, &
      max_iterations, tolerance) result(out)
      !! Estimate a VARMA model by conditional Gaussian likelihood.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      integer, intent(in) :: ar_order !! Autoregressive order.
      integer, intent(in) :: ma_order !! Moving-average order.
      logical, intent(in), optional :: include_mean !! Whether to include a mean term.
      logical, intent(in), optional :: estimated(:) !! Flag controlling estimated.
      real(dp), intent(in), optional :: initial(:) !! Initial value.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      type(mts_varma_fit_t) :: out
      type(optimization_result_t) :: optimization
      type(mts_var_fit_t) :: var_start
      real(dp), allocatable :: parameters(:), free_initial(:), hessian(:, :), inverse(:, :)
      integer, allocatable :: free_index(:)
      logical :: use_mean
      integer :: n, dimension, block_size, count, free_count, limit, offset, lag, i, status
      real(dp) :: gradient_tolerance

      n = size(series, 1)
      dimension = size(series, 2)
      use_mean = .true.
      if (present(include_mean)) use_mean = include_mean
      block_size = dimension*dimension
      count = (ar_order + ma_order)*block_size + merge(dimension, 0, use_mean)
      if (ar_order < 0 .or. ma_order < 0 .or. ar_order + ma_order < 1 .or. &
         dimension < 1 .or. n <= max(ar_order, ma_order) + 2) then
         out%info = 1
         return
      end if
      allocate(parameters(count), out%estimated(count))
      parameters = 0.0_dp
      offset = 0
      if (use_mean) then
         parameters(:dimension) = sum(series, 1)/real(n, dp)
         offset = dimension
      end if
      if (ar_order > 0) then
         var_start = mts_var(series, ar_order, use_mean)
         if (var_start%info == 0) then
            if (use_mean) parameters(:dimension) = var_start%intercept
            do lag = 1, ar_order
               parameters(offset + 1:offset + block_size) = reshape(var_start%ar(:, :, lag), [block_size])
               offset = offset + block_size
            end do
         end if
      end if
      offset = merge(dimension, 0, use_mean) + ar_order*block_size
      if (present(initial)) then
         if (size(initial) /= count) then
            out%info = 2
            return
         end if
         parameters = initial
      end if
      out%estimated = .true.
      if (present(estimated)) then
         if (size(estimated) /= count) then
            out%info = 2
            return
         end if
         out%estimated = estimated
         if (.not. present(initial)) then
            where (.not. out%estimated)
               parameters = 0.0_dp
            end where
         end if
      end if
      free_count = sum(merge(1, 0, out%estimated))
      allocate(free_index(free_count), free_initial(free_count))
      free_index = pack([(i, i=1, count)], out%estimated)
      free_initial = pack(parameters, out%estimated)
      limit = 200
      if (present(max_iterations)) limit = max_iterations
      gradient_tolerance = 1.0e-6_dp
      if (present(tolerance)) gradient_tolerance = tolerance
      if (free_count > 0) then
         optimization = bfgs_minimize_fd(objective, free_initial, limit, gradient_tolerance)
         parameters(free_index) = optimization%parameters
         out%iterations = optimization%iterations
         out%converged = optimization%converged
         if (optimization%info /= 0) out%info = 100 + optimization%info
      else
         out%converged = .true.
      end if
      call evaluate(parameters, out%model, out%residuals, out%log_likelihood, status)
      if (status /= 0) then
         out%info = 10 + status
         return
      end if
      out%coefficients = parameters
      out%aic = -2.0_dp*out%log_likelihood + 2.0_dp*real(free_count + 1, dp)
      out%bic = -2.0_dp*out%log_likelihood + log(real(n, dp))*real(free_count + 1, dp)
      out%stationarity_radius = companion_radius(out%model%ar)
      out%invertibility_radius = companion_radius(out%model%ma)
      allocate(out%covariance(count, count), out%standard_errors(count))
      out%covariance = 0.0_dp
      out%standard_errors = 0.0_dp
      if (free_count > 0) then
         free_initial = pack(parameters, out%estimated)
         hessian = finite_difference_hessian(objective, free_initial)
         call invert_matrix(hessian, inverse, status)
         if (status == 0) then
            do i = 1, free_count
               out%covariance(free_index, free_index(i)) = inverse(:, i)
               out%standard_errors(free_index(i)) = sqrt(max(0.0_dp, inverse(i, i)))
            end do
         else if (out%info == 0) then
            out%info = 200 + status
         end if
      end if
      out%ar_lags = [(i, i=1, ar_order)]
      out%ma_lags = [(i, i=1, ma_order)]

   contains

      pure function objective(free_parameters) result(value)
         !! Return the conditional VARMA negative log likelihood.
         real(dp), intent(in) :: free_parameters(:) !! Free parameters.
         real(dp) :: value, log_likelihood
         real(dp) :: full_parameters(count)
         real(dp), allocatable :: residual_values(:, :)
         type(mts_varma_model_t) :: candidate
         integer :: evaluation_status

         full_parameters = parameters
         full_parameters(free_index) = free_parameters
         call evaluate(full_parameters, candidate, residual_values, log_likelihood, evaluation_status)
         if (evaluation_status == 0 .and. ieee_is_finite(log_likelihood)) then
            value = -log_likelihood
         else
            value = 1.0e30_dp + dot_product(free_parameters, free_parameters)
         end if
      end function objective

      pure subroutine evaluate(values, model, residual_values, log_likelihood, evaluation_status)
         !! Unpack parameters and evaluate conditional VARMA residual likelihood.
         real(dp), intent(in) :: values(:) !! Input values.
         type(mts_varma_model_t), intent(out) :: model !! Model specification.
         real(dp), allocatable, intent(out) :: residual_values(:, :) !! Residual values.
         real(dp), intent(out) :: log_likelihood !! Log-likelihood value.
         integer, intent(out) :: evaluation_status !! Evaluation status.
         real(dp), allocatable :: covariance_inverse(:, :)
         real(dp) :: logdet
         integer :: parameter_offset, current_lag, t

         allocate(model%ar(dimension, dimension, ar_order), model%ma(dimension, dimension, ma_order))
         allocate(model%intercept(dimension), model%sigma(dimension, dimension))
         model%intercept = 0.0_dp
         parameter_offset = 0
         if (use_mean) then
            model%intercept = values(:dimension)
            parameter_offset = dimension
         end if
         do current_lag = 1, ar_order
            model%ar(:, :, current_lag) = reshape(&
               values(parameter_offset + 1:parameter_offset + block_size), [dimension, dimension])
            parameter_offset = parameter_offset + block_size
         end do
         do current_lag = 1, ma_order
            model%ma(:, :, current_lag) = reshape(&
               values(parameter_offset + 1:parameter_offset + block_size), [dimension, dimension])
            parameter_offset = parameter_offset + block_size
         end do
         allocate(residual_values(n, dimension))
         do t = 1, n
            residual_values(t, :) = series(t, :) - model%intercept
            do current_lag = 1, min(ar_order, t - 1)
               residual_values(t, :) = residual_values(t, :) - &
                  matmul(model%ar(:, :, current_lag), series(t - current_lag, :))
            end do
            do current_lag = 1, min(ma_order, t - 1)
               residual_values(t, :) = residual_values(t, :) + &
                  matmul(model%ma(:, :, current_lag), residual_values(t - current_lag, :))
            end do
         end do
         model%sigma = matmul(transpose(residual_values), residual_values)/real(n, dp)
         allocate(covariance_inverse(dimension, dimension))
         call inverse_logdet(model%sigma, covariance_inverse, logdet, evaluation_status, &
            100.0_dp*epsilon(1.0_dp))
         if (evaluation_status /= 0) then
            log_likelihood = -huge(1.0_dp)
            return
         end if
         log_likelihood = -0.5_dp*real(n, dp)*(real(dimension, dp)*(log(2.0_dp*acos(-1.0_dp)) + &
            1.0_dp) + logdet)
      end subroutine evaluate
   end function mts_varma_fit

   pure function mts_varmas_fit(series, ar_lags, ma_lags, include_mean, &
      max_iterations, tolerance) result(out)
      !! Estimate a VARMA model with selected nonzero AR and MA lags.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      integer, intent(in) :: ar_lags(:) !! Autoregressive lags.
      integer, intent(in) :: ma_lags(:) !! Moving-average lags.
      logical, intent(in), optional :: include_mean !! Whether to include a mean term.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      type(mts_varma_fit_t) :: out
      logical, allocatable :: mask(:)
      integer :: dimension, ar_order, ma_order, block_size, offset, i
      logical :: use_mean

      if ((size(ar_lags) > 0 .and. any(ar_lags < 1)) .or. &
         (size(ma_lags) > 0 .and. any(ma_lags < 1)) .or. &
         size(ar_lags) + size(ma_lags) < 1) then
         out%info = 1
         return
      end if
      dimension = size(series, 2)
      ar_order = 0
      ma_order = 0
      if (size(ar_lags) > 0) ar_order = maxval(ar_lags)
      if (size(ma_lags) > 0) ma_order = maxval(ma_lags)
      block_size = dimension*dimension
      use_mean = .true.
      if (present(include_mean)) use_mean = include_mean
      allocate(mask((ar_order + ma_order)*block_size + merge(dimension, 0, use_mean)))
      mask = .false.
      offset = 0
      if (use_mean) then
         mask(:dimension) = .true.
         offset = dimension
      end if
      do i = 1, size(ar_lags)
         mask(offset + (ar_lags(i) - 1)*block_size + 1:offset + ar_lags(i)*block_size) = .true.
      end do
      offset = offset + ar_order*block_size
      do i = 1, size(ma_lags)
         mask(offset + (ma_lags(i) - 1)*block_size + 1:offset + ma_lags(i)*block_size) = .true.
      end do
      out = mts_varma_fit(series, ar_order, ma_order, use_mean, estimated=mask, &
         max_iterations=max_iterations, tolerance=tolerance)
      if (out%model%info == 0) then
         out%ar_lags = ar_lags
         out%ma_lags = ma_lags
      end if
   end function mts_varmas_fit

   pure function mts_refine_varma(series, ar_order, ma_order, threshold, include_mean, &
      protected, max_steps, max_iterations, tolerance) result(out)
      !! Refine a VARMA model by iterative t-ratio backward elimination.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      integer, intent(in) :: ar_order !! Autoregressive order.
      integer, intent(in) :: ma_order !! Moving-average order.
      real(dp), intent(in), optional :: threshold !! Decision or truncation threshold.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      logical, intent(in), optional :: include_mean !! Whether to include a mean term.
      logical, intent(in), optional :: protected(:) !! Flag controlling protected.
      integer, intent(in), optional :: max_steps !! Maximum steps.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      type(mts_varma_refinement_t) :: out
      logical, allocatable :: active(:), keep(:)
      real(dp), allocatable :: starts(:)
      real(dp) :: cutoff, ratio, weakest_ratio, parameter_cutoff
      integer :: count_parameters, dimension, limit, step, i, weakest
      logical :: use_mean

      dimension = size(series, 2)
      use_mean = .true.
      if (present(include_mean)) use_mean = include_mean
      count_parameters = (ar_order + ma_order)*dimension*dimension + merge(dimension, 0, use_mean)
      cutoff = 1.5_dp
      if (present(threshold)) cutoff = threshold
      if (count_parameters < 1 .or. cutoff < 0.0_dp) then
         out%info = 1
         return
      end if
      allocate(active(count_parameters), keep(count_parameters))
      active = .true.
      keep = .false.
      if (present(protected)) then
         if (size(protected) /= count_parameters) then
            out%info = 2
            return
         end if
         keep = protected
      end if
      limit = count_parameters
      if (present(max_steps)) limit = min(limit, max(0, max_steps))
      allocate(out%active_count(0:limit), out%removed_index(limit))
      out%fit = mts_varma_fit(series, ar_order, ma_order, use_mean, &
         max_iterations=max_iterations, tolerance=tolerance)
      if (out%fit%model%info /= 0) then
         out%info = 10 + out%fit%model%info
         return
      end if
      out%active_count(0) = count(active)
      do step = 1, limit
         weakest = 0
         weakest_ratio = huge(1.0_dp)
         do i = 1, count_parameters
            if (.not. active(i) .or. keep(i)) cycle
            if (out%fit%standard_errors(i) <= tiny(1.0_dp)) cycle
            ratio = abs(out%fit%coefficients(i)/out%fit%standard_errors(i))
            parameter_cutoff = cutoff
            if (use_mean .and. i <= dimension) parameter_cutoff = min(cutoff, 1.0_dp)
            if (ratio < parameter_cutoff .and. ratio < weakest_ratio) then
               weakest = i
               weakest_ratio = ratio
            end if
         end do
         if (weakest == 0) exit
         active(weakest) = .false.
         starts = out%fit%coefficients
         starts(weakest) = 0.0_dp
         out%fit = mts_varma_fit(series, ar_order, ma_order, use_mean, starts, active, &
            max_iterations=max_iterations, tolerance=tolerance)
         out%steps = step
         out%removed_index(step) = weakest
         out%active_count(step) = count(active)
         if (out%fit%model%info /= 0) then
            out%info = 20 + out%fit%model%info
            exit
         end if
      end do
      out%active_count = out%active_count(:out%steps)
      out%removed_index = out%removed_index(:out%steps)
   end function mts_refine_varma

   pure function mts_refine_svarma(series, regular_ar_order, regular_ma_order, &
      seasonal_ar_order, seasonal_ma_order, period, threshold, include_mean, &
      switched, protected, max_steps, max_iterations, tolerance) result(out)
      !! Refine a seasonal VARMA model by iterative t-ratio elimination.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      integer, intent(in) :: regular_ar_order !! Regular autoregressive order.
      integer, intent(in) :: regular_ma_order !! Regular moving-average order.
      integer, intent(in) :: seasonal_ar_order !! Seasonal autoregressive order.
      integer, intent(in) :: seasonal_ma_order !! Seasonal moving-average order.
      integer, intent(in) :: period !! Seasonal period.
      real(dp), intent(in), optional :: threshold !! Decision or truncation threshold.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      logical, intent(in), optional :: include_mean !! Whether to include a mean term.
      logical, intent(in), optional :: switched !! Flag controlling switched.
      logical, intent(in), optional :: protected(:) !! Flag controlling protected.
      integer, intent(in), optional :: max_steps !! Maximum steps.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      type(mts_svarma_refinement_t) :: out
      logical, allocatable :: active(:), keep(:)
      real(dp), allocatable :: starts(:)
      real(dp) :: cutoff, ratio, weakest_ratio, parameter_cutoff
      integer :: count_parameters, component_count, dimension, limit, step, i, weakest
      logical :: use_mean, reverse_products

      dimension = size(series, 2)
      use_mean = .true.
      reverse_products = .false.
      if (present(include_mean)) use_mean = include_mean
      if (present(switched)) reverse_products = switched
      component_count = regular_ar_order + regular_ma_order + &
         seasonal_ar_order + seasonal_ma_order
      count_parameters = component_count*dimension*dimension + merge(dimension, 0, use_mean)
      cutoff = 0.8_dp
      if (present(threshold)) cutoff = threshold
      if (count_parameters < 1 .or. cutoff < 0.0_dp) then
         out%info = 1
         return
      end if
      allocate(active(count_parameters), keep(count_parameters))
      active = .true.
      keep = .false.
      if (present(protected)) then
         if (size(protected) /= count_parameters) then
            out%info = 2
            return
         end if
         keep = protected
      end if
      limit = count_parameters
      if (present(max_steps)) limit = min(limit, max(0, max_steps))
      allocate(out%active_count(0:limit), out%removed_index(limit))
      out%fit = mts_svarma_fit(series, regular_ar_order, regular_ma_order, &
         seasonal_ar_order, seasonal_ma_order, period, use_mean, reverse_products, &
         max_iterations=max_iterations, tolerance=tolerance)
      if (out%fit%model%info /= 0) then
         out%info = 10 + out%fit%model%info
         return
      end if
      out%active_count(0) = count(active)
      do step = 1, limit
         weakest = 0
         weakest_ratio = huge(1.0_dp)
         do i = 1, count_parameters
            if (.not. active(i) .or. keep(i)) cycle
            if (out%fit%standard_errors(i) <= tiny(1.0_dp)) cycle
            ratio = abs(out%fit%coefficients(i)/out%fit%standard_errors(i))
            parameter_cutoff = cutoff
            if (use_mean .and. i <= dimension) parameter_cutoff = min(cutoff, 1.0_dp)
            if (ratio < parameter_cutoff .and. ratio < weakest_ratio) then
               weakest = i
               weakest_ratio = ratio
            end if
         end do
         if (weakest == 0) exit
         active(weakest) = .false.
         starts = out%fit%coefficients
         starts(weakest) = 0.0_dp
         out%fit = mts_svarma_fit(series, regular_ar_order, regular_ma_order, &
            seasonal_ar_order, seasonal_ma_order, period, use_mean, reverse_products, &
            starts, active, max_iterations, tolerance)
         out%steps = step
         out%removed_index(step) = weakest
         out%active_count(step) = count(active)
         if (out%fit%model%info /= 0) then
            out%info = 20 + out%fit%model%info
            exit
         end if
      end do
      out%active_count = out%active_count(:out%steps)
      out%removed_index = out%removed_index(:out%steps)
   end function mts_refine_svarma

   pure function mts_vma_fit(series, order, include_mean, initial, estimated, &
      max_iterations, tolerance) result(out)
      !! Estimate a consecutive-lag VMA model by conditional Gaussian likelihood.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      integer, intent(in) :: order !! Model or polynomial order.
      logical, intent(in), optional :: include_mean !! Whether to include a mean term.
      logical, intent(in), optional :: estimated(:) !! Flag controlling estimated.
      real(dp), intent(in), optional :: initial(:) !! Initial value.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      type(mts_vma_fit_t) :: out
      type(optimization_result_t) :: optimization
      real(dp), allocatable :: parameters(:), free_initial(:), hessian(:, :), inverse(:, :)
      integer, allocatable :: free_index(:)
      logical :: use_mean
      integer :: n, dimension, count, free_count, limit, i, status
      real(dp) :: gradient_tolerance

      n = size(series, 1)
      dimension = size(series, 2)
      use_mean = .true.
      if (present(include_mean)) use_mean = include_mean
      count = order*dimension*dimension + merge(dimension, 0, use_mean)
      if (order < 1 .or. dimension < 1 .or. n <= order + 2) then
         out%info = 1
         return
      end if
      allocate(parameters(count), out%estimated(count))
      parameters = 0.0_dp
      if (use_mean) parameters(:dimension) = sum(series, 1)/real(n, dp)
      if (present(initial)) then
         if (size(initial) /= count) then
            out%info = 2
            return
         end if
         parameters = initial
      end if
      out%estimated = .true.
      if (present(estimated)) then
         if (size(estimated) /= count) then
            out%info = 2
            return
         end if
         out%estimated = estimated
      end if
      free_count = sum(merge(1, 0, out%estimated))
      allocate(free_index(free_count), free_initial(free_count))
      free_index = pack([(i, i=1, count)], out%estimated)
      free_initial = pack(parameters, out%estimated)
      limit = 200
      if (present(max_iterations)) limit = max_iterations
      gradient_tolerance = 1.0e-6_dp
      if (present(tolerance)) gradient_tolerance = tolerance
      if (free_count > 0) then
         optimization = bfgs_minimize_fd(objective, free_initial, limit, gradient_tolerance)
         parameters(free_index) = optimization%parameters
         out%iterations = optimization%iterations
         out%converged = optimization%converged
         if (optimization%info /= 0) out%info = 100 + optimization%info
      else
         out%converged = .true.
      end if
      call evaluate(parameters, out%model, out%residuals, out%log_likelihood, status)
      if (status /= 0) then
         out%info = 10 + status
         return
      end if
      out%coefficients = parameters
      out%aic = -2.0_dp*out%log_likelihood + 2.0_dp*real(free_count + 1, dp)
      out%bic = -2.0_dp*out%log_likelihood + log(real(n, dp))*real(free_count + 1, dp)
      out%invertibility_radius = companion_radius(out%model%ma)
      allocate(out%covariance(count, count), out%standard_errors(count))
      out%covariance = 0.0_dp
      out%standard_errors = 0.0_dp
      if (free_count > 0) then
         free_initial = pack(parameters, out%estimated)
         hessian = finite_difference_hessian(objective, free_initial)
         call invert_matrix(hessian, inverse, status)
         if (status == 0) then
            do i = 1, free_count
               out%covariance(free_index, free_index(i)) = inverse(:, i)
               out%standard_errors(free_index(i)) = sqrt(max(0.0_dp, inverse(i, i)))
            end do
         else if (out%info == 0) then
            out%info = 200 + status
         end if
      end if
      out%lags = [(i, i=1, order)]

   contains

      pure function objective(free_parameters) result(value)
         !! Return the conditional VMA negative log likelihood.
         real(dp), intent(in) :: free_parameters(:) !! Free parameters.
         real(dp) :: value, log_likelihood
         real(dp) :: full_parameters(count)
         real(dp), allocatable :: residual_values(:, :)
         type(mts_varma_model_t) :: candidate
         integer :: evaluation_status

         full_parameters = parameters
         full_parameters(free_index) = free_parameters
         call evaluate(full_parameters, candidate, residual_values, log_likelihood, evaluation_status)
         if (evaluation_status == 0 .and. ieee_is_finite(log_likelihood)) then
            value = -log_likelihood
         else
            value = 1.0e30_dp + dot_product(free_parameters, free_parameters)
         end if
      end function objective

      pure subroutine evaluate(values, model, residual_values, log_likelihood, evaluation_status)
         !! Unpack parameters and evaluate conditional VMA residual likelihood.
         real(dp), intent(in) :: values(:) !! Input values.
         type(mts_varma_model_t), intent(out) :: model !! Model specification.
         real(dp), allocatable, intent(out) :: residual_values(:, :) !! Residual values.
         real(dp), intent(out) :: log_likelihood !! Log-likelihood value.
         integer, intent(out) :: evaluation_status !! Evaluation status.
         real(dp), allocatable :: centered(:, :), covariance_inverse(:, :)
         real(dp) :: logdet
         integer :: offset, lag, t

         allocate(model%ar(dimension, dimension, 0), model%ma(dimension, dimension, order))
         allocate(model%intercept(dimension), model%sigma(dimension, dimension))
         model%intercept = 0.0_dp
         offset = 0
         if (use_mean) then
            model%intercept = values(:dimension)
            offset = dimension
         end if
         do lag = 1, order
            model%ma(:, :, lag) = reshape(values(offset + 1:offset + dimension*dimension), &
               [dimension, dimension])
            offset = offset + dimension*dimension
         end do
         centered = series - spread(model%intercept, 1, n)
         allocate(residual_values(n, dimension))
         do t = 1, n
            residual_values(t, :) = centered(t, :)
            do lag = 1, min(order, t - 1)
               residual_values(t, :) = residual_values(t, :) + &
                  matmul(model%ma(:, :, lag), residual_values(t - lag, :))
            end do
         end do
         model%sigma = matmul(transpose(residual_values), residual_values)/real(n, dp)
         allocate(covariance_inverse(dimension, dimension))
         call inverse_logdet(model%sigma, covariance_inverse, logdet, evaluation_status, &
            100.0_dp*epsilon(1.0_dp))
         if (evaluation_status /= 0) then
            log_likelihood = -huge(1.0_dp)
            return
         end if
         log_likelihood = -0.5_dp*real(n, dp)*(real(dimension, dp)*(log(2.0_dp*acos(-1.0_dp)) + &
            1.0_dp) + logdet)
      end subroutine evaluate
   end function mts_vma_fit

   pure function mts_vmas_fit(series, lags, include_mean, max_iterations, tolerance) result(out)
      !! Estimate a VMA model with selected nonzero MA lags.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      integer, intent(in) :: lags(:) !! Lags.
      logical, intent(in), optional :: include_mean !! Whether to include a mean term.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      type(mts_vma_fit_t) :: out
      logical, allocatable :: mask(:)
      integer :: dimension, maximum_lag, offset, lag, i
      logical :: use_mean

      if (size(lags) < 1 .or. any(lags < 1)) then
         out%info = 1
         return
      end if
      dimension = size(series, 2)
      maximum_lag = maxval(lags)
      use_mean = .true.
      if (present(include_mean)) use_mean = include_mean
      allocate(mask(maximum_lag*dimension*dimension + merge(dimension, 0, use_mean)))
      mask = .false.
      offset = 0
      if (use_mean) then
         mask(:dimension) = .true.
         offset = dimension
      end if
      do i = 1, size(lags)
         lag = lags(i)
         mask(offset + (lag - 1)*dimension*dimension + 1: &
            offset + lag*dimension*dimension) = .true.
      end do
      out = mts_vma_fit(series, maximum_lag, use_mean, estimated=mask, &
         max_iterations=max_iterations, tolerance=tolerance)
      if (out%info == 0 .or. out%info >= 100) out%lags = lags
   end function mts_vmas_fit

   pure function mts_vma_order(series, max_order, include_mean, max_iterations, tolerance) result(out)
      !! Select consecutive VMA order by conditional-likelihood AIC and BIC.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      integer, intent(in) :: max_order !! Maximum order.
      logical, intent(in), optional :: include_mean !! Whether to include a mean term.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      type(mts_vma_order_t) :: out
      type(mts_vma_fit_t) :: fitted
      integer :: order

      if (max_order < 1) then
         out%info = 1
         return
      end if
      allocate(out%aic(1:max_order), out%bic(1:max_order))
      do order = 1, max_order
         fitted = mts_vma_fit(series, order, include_mean, max_iterations=max_iterations, tolerance=tolerance)
         if (fitted%model%info /= 0 .or. .not. ieee_is_finite(fitted%log_likelihood)) then
            out%info = 10 + order
            return
         end if
         out%aic(order) = fitted%aic
         out%bic(order) = fitted%bic
      end do
      out%aic_order = minloc(out%aic, dim=1)
      out%bic_order = minloc(out%bic, dim=1)
   end function mts_vma_order

   pure real(dp) function companion_radius(ma) result(radius)
      !! Estimate the dominant modulus of the VMA inverse-recursion companion.
      real(dp), intent(in) :: ma(:, :, :) !! Moving-average coefficients.
      real(dp), allocatable :: companion(:, :), vector(:), next_vector(:)
      real(dp) :: scale
      integer :: dimension, order, state_size, lag, i, iteration

      dimension = size(ma, 1)
      order = size(ma, 3)
      state_size = dimension*order
      if (state_size == 0) then
         radius = 0.0_dp
         return
      end if
      allocate(companion(state_size, state_size), vector(state_size), next_vector(state_size))
      companion = 0.0_dp
      do lag = 1, order
         companion(:dimension, (lag - 1)*dimension + 1:lag*dimension) = ma(:, :, lag)
      end do
      do i = 1, state_size - dimension
         companion(dimension + i, i) = 1.0_dp
      end do
      vector = 1.0_dp/sqrt(real(state_size, dp))
      radius = 0.0_dp
      do iteration = 1, 200
         next_vector = matmul(companion, vector)
         scale = norm2(next_vector)
         if (scale <= tiny(1.0_dp)) then
            radius = 0.0_dp
            return
         end if
         vector = next_vector/scale
         radius = scale
      end do
   end function companion_radius

   pure function mts_var_order(series, max_order, common_sample) result(out)
      !! Select VAR order using MTS AIC, BIC, and Hannan-Quinn criteria.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      integer, intent(in) :: max_order !! Maximum order.
      logical, intent(in), optional :: common_sample !! Flag controlling common sample.
      type(mts_var_order_t) :: out
      real(dp), allocatable :: comparison(:, :), covariance(:, :), inverse(:, :), log_determinant(:)
      type(mts_var_fit_t) :: fitted
      real(dp) :: logdet
      integer :: n, dimension, order, status
      logical :: use_common

      n = size(series, 1)
      dimension = size(series, 2)
      use_common = .true.
      if (present(common_sample)) use_common = common_sample
      if (max_order < 1 .or. n - max_order <= dimension*max_order + 1) then
         out%info = 1
         return
      end if
      allocate(out%aic(0:max_order), out%bic(0:max_order), out%hq(0:max_order))
      allocate(out%statistic(max_order), out%p_value(max_order))
      allocate(log_determinant(0:max_order))
      comparison = series
      if (use_common) comparison = series(max_order + 1:, :)
      covariance = matmul(transpose(comparison - spread(sum(comparison, 1)/ &
         real(size(comparison, 1), dp), 1, size(comparison, 1))), &
         comparison - spread(sum(comparison, 1)/real(size(comparison, 1), dp), 1, size(comparison, 1)))/ &
         real(size(comparison, 1), dp)
      allocate(inverse(dimension, dimension))
      call inverse_logdet(covariance, inverse, logdet, status, 100.0_dp*epsilon(1.0_dp))
      if (status /= 0) then
         out%info = 2
         return
      end if
      out%aic(0) = logdet
      out%bic(0) = logdet
      out%hq(0) = logdet
      log_determinant(0) = logdet
      do order = 1, max_order
         if (use_common) then
            fitted = fit_var_common_sample(series, order, max_order)
         else
            fitted = mts_var(series, order)
         end if
         if (fitted%info /= 0) then
            out%info = 10 + order
            return
         end if
         call inverse_logdet(fitted%sigma, inverse, logdet, status, 100.0_dp*epsilon(1.0_dp))
         log_determinant(order) = logdet
         out%aic(order) = logdet + 2.0_dp*real(order*dimension*dimension, dp)/real(n, dp)
         out%bic(order) = logdet + log(real(n, dp))*real(order*dimension*dimension, dp)/real(n, dp)
         out%hq(order) = logdet + 2.0_dp*log(log(real(n, dp)))* &
            real(order*dimension*dimension, dp)/real(n, dp)
      end do
      out%aic_order = minloc(out%aic, dim=1) - 1
      out%bic_order = minloc(out%bic, dim=1) - 1
      out%hq_order = minloc(out%hq, dim=1) - 1
      do order = 1, max_order
         out%statistic(order) = (real(n - max_order - dimension*order, dp) - 1.5_dp)* &
            (log_determinant(order - 1) - log_determinant(order))
         out%p_value(order) = chi_square_survival(out%statistic(order), dimension*dimension)
      end do
   end function mts_var_order

   pure function fit_var_common_sample(series, order, maximum_order) result(out)
      !! Fit a VAR order using the response sample aligned to maximum_order.
      real(dp), intent(in) :: series(:, :) !! Time-series observations.
      integer, intent(in) :: order !! Model or polynomial order.
      integer, intent(in) :: maximum_order !! Maximum order.
      type(mts_var_fit_t) :: out
      real(dp), allocatable :: truncated(:, :)
      truncated = series(maximum_order - order + 1:, :)
      out = mts_var(truncated, order)
   end function fit_var_common_sample

   pure real(dp) function chi_square_survival(value, degrees) result(probability)
      !! Approximate a chi-square upper tail by Wilson-Hilferty transformation.
      real(dp), intent(in) :: value !! Input value.
      integer, intent(in) :: degrees !! Degrees.
      real(dp) :: z, d
      d = real(degrees, dp)
      z = ((max(value, 0.0_dp)/d)**(1.0_dp/3.0_dp) - (1.0_dp - 2.0_dp/(9.0_dp*d)))/ &
         sqrt(2.0_dp/(9.0_dp*d))
      probability = 0.5_dp*erfc(z/sqrt(2.0_dp))
   end function chi_square_survival
end module mts_mod
