! SPDX-License-Identifier: GPL-2.0-only
! SPDX-FileComment: Algorithms translated from the R tsissm package.
! Distinct innovations state-space algorithms translated from the tsissm package.
module tsissm_mod
   use kind_mod, only: dp
   use forecast_mod, only: box_cox, accuracy_result_t, &
      forecast_accuracy
   use itsmr_mod, only: itsmr_randomness_tests_t, residual_randomness_tests
   use special_functions_mod, only: regularized_gamma_q
   use mts_mod, only: mts_arch_test_t, mts_arch_test, regularized_beta_mts
   use optimization_mod, only: optimization_result_t, bfgs_minimize_fd, &
      finite_difference_hessian
   use linalg_mod, only: invert_matrix, symmetric_eigen, outer_product
   use stats_mod, only: sort, quantile
   use time_series_stats_mod, only: acf_values, pacf_values
   use random_mod, only: set_random_seed, random_uniform, random_standard_normal, &
      random_standard_student, random_standard_johnson_su
   use time_series_diagnostics_mod, only: &
      tsissm_box_test_box_pierce => box_test_box_pierce, &
      tsissm_box_test_ljung_box => box_test_ljung_box, &
      tsissm_box_test_monti => box_test_monti, &
      tsissm_residual_raw => residual_raw, &
      tsissm_residual_squared => residual_squared, &
      tsissm_residual_log_squared => residual_log_squared, &
      tsissm_residual_absolute => residual_absolute, &
      tsissm_weighted_box_test_t => weighted_box_test_t, &
      shared_weighted_box_test => weighted_box_test
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   private

   integer, parameter, public :: tsissm_distribution_gaussian = 1
   integer, parameter, public :: tsissm_distribution_student = 2
   integer, parameter, public :: tsissm_distribution_johnson_su = 7
   integer, parameter, public :: tsissm_selection_aic = 1
   integer, parameter, public :: tsissm_selection_aicc = 2
   integer, parameter, public :: tsissm_selection_bic = 3
   integer, parameter, public :: tsissm_selection_external = 4
   integer, parameter, public :: tsissm_covariance_hessian = 1
   integer, parameter, public :: tsissm_covariance_opg = 2
   integer, parameter, public :: tsissm_covariance_qmle = 3
   integer, parameter, public :: tsissm_covariance_newey_west = 4
   integer, parameter, public :: tsissm_variance_initial_full = 1
   integer, parameter, public :: tsissm_variance_initial_sample = 2
   integer, parameter, public :: tsissm_hac_bartlett = 1
   integer, parameter, public :: tsissm_hac_parzen = 2
   integer, parameter, public :: tsissm_hac_quadratic_spectral = 3
   integer, parameter, public :: tsissm_hac_truncated = 4
   integer, parameter, public :: tsissm_hac_tukey_hanning = 5

   type, public :: tsissm_filter_t
      ! Innovations state-space filter output and system matrices.
      real(dp), allocatable :: state(:, :), fitted(:), transformed(:), innovation(:)
      real(dp), allocatable :: transition(:, :), discount_transition(:, :)
      real(dp), allocatable :: observation_loading(:), persistence(:)
      real(dp), allocatable :: conditional_sd(:)
      integer :: observations = 0
      integer :: info = 0
   end type

   type, public :: tsissm_prediction_t
      ! Simulated future observations and innovations state paths.
      real(dp), allocatable :: observation(:, :), state(:, :, :)
      real(dp), allocatable :: innovation(:, :), conditional_sd(:, :)
      integer :: info = 0
   end type

   type, public :: tsissm_garch_t
      ! Conditional GARCH variances and standard deviations.
      real(dp), allocatable :: variance(:), standard_deviation(:)
      integer :: info = 0
   end type

   type, public :: tsissm_variance_initialization_t
      ! Initial GARCH variance with the number of residuals used.
      real(dp) :: variance = 0.0_dp
      integer :: observations = 0
      integer :: method = tsissm_variance_initial_full
      integer :: info = 0
   end type

   type, public :: tsissm_model_t
      ! Structural ISSM matrices and component boundaries.
      real(dp), allocatable :: transition_base(:, :), transition_scale(:, :)
      real(dp), allocatable :: transition_parameter(:, :), observation_loading(:)
      real(dp), allocatable :: persistence(:), initial_state(:)
      integer, allocatable :: seasonal_start(:), seasonal_end(:)
      integer :: level_index = 1
      integer :: slope_index = 0
      integer :: ar_start = 0
      integer :: ar_end = 0
      integer :: ma_start = 0
      integer :: ma_end = 0
      integer :: info = 0
   end type

   type, public :: tsissm_likelihood_t
      ! Innovation likelihood contributions, scales, residuals, and criteria.
      real(dp), allocatable :: contribution(:), standardized_residual(:), scale(:)
      real(dp) :: objective = 0.0_dp
      real(dp) :: log_likelihood = 0.0_dp
      real(dp) :: package_log_likelihood = 0.0_dp
      real(dp) :: sigma = 0.0_dp
      real(dp) :: aic = 0.0_dp
      real(dp) :: bic = 0.0_dp
      integer :: observations = 0
      integer :: info = 0
   end type

   type, public :: tsissm_constraint_t
      ! Constraint residuals using nonnegative values for feasible parameters.
      real(dp), allocatable :: residual(:)
      logical :: feasible = .false.
      integer :: info = 0
   end type

   type, public :: tsissm_initialization_t
      ! Least-squares seasonal/ARMA seed states and status.
      real(dp), allocatable :: state(:)
      integer :: info = 0
   end type

   type, public :: tsissm_fit_t
      ! Constant-variance ISSM estimates, filter output, and likelihood.
      real(dp), allocatable :: persistence(:), coefficients(:), arch(:), garch(:)
      type(tsissm_filter_t) :: filter
      type(tsissm_likelihood_t) :: likelihood
      real(dp) :: objective = huge(1.0_dp)
      real(dp) :: variance_intercept = 0.0_dp
      real(dp) :: initial_variance = 0.0_dp
      integer :: iterations = 0
      integer :: info = 0
      logical :: converged = .false.
   end type

   type, public :: tsissm_structural_fit_t
      ! Joint structural estimates and their final model, filter, and likelihood.
      real(dp), allocatable :: parameters(:), arch(:), garch(:)
      type(tsissm_model_t) :: model
      type(tsissm_filter_t) :: filter
      type(tsissm_likelihood_t) :: likelihood
      real(dp) :: lambda = 1.0_dp
      real(dp) :: skew = 0.0_dp
      real(dp) :: shape = 0.0_dp
      real(dp) :: variance_intercept = 0.0_dp
      real(dp) :: initial_variance = 0.0_dp
      real(dp) :: objective = huge(1.0_dp)
      integer :: iterations = 0
      integer :: info = 0
      logical :: converged = .false.
      logical :: dynamic_variance = .false.
   end type

   type, public :: tsissm_moments_t
      ! Analytical forecast means and variances on transformed and original scales.
      real(dp), allocatable :: transformed_mean(:), transformed_variance(:)
      real(dp), allocatable :: mean(:), variance(:), garch_variance(:)
      integer :: info = 0
   end type

   type, public :: tsissm_hresiduals_t
      ! Multi-horizon in-sample forecasts, errors, and cross-horizon summaries.
      real(dp), allocatable :: forecast(:, :), actual(:, :), residual(:, :)
      real(dp), allocatable :: covariance(:, :), mean_residual(:), rmse(:)
      integer, allocatable :: pair_count(:, :)
      logical, allocatable :: valid(:, :)
      integer :: first_origin = 0
      integer :: horizon = 0
      integer :: info = 0
      logical :: transformed = .false.
   end type

   type, public :: tsissm_covariance_t
      ! Parameter covariance components, standard errors, and correlations.
      real(dp), allocatable :: covariance(:, :), bread(:, :), meat(:, :)
      real(dp), allocatable :: standard_error(:), correlation(:, :)
      real(dp), allocatable :: scores(:, :), hessian(:, :)
      real(dp) :: bandwidth = 0.0_dp
      real(dp) :: bias_correction = 1.0_dp
      real(dp) :: effective_degrees = 0.0_dp
      integer :: covariance_type = tsissm_covariance_hessian
      integer :: observations = 0
      integer :: parameters = 0
      integer :: lag = 0
      integer :: kernel = tsissm_hac_bartlett
      integer :: prewhite_order = 0
      integer :: info = 0
      logical :: adjusted = .false.
   end type

   type, public :: tsissm_decomposition_t
      ! Filtered state contributions and reconstructed observations.
      real(dp), allocatable :: level(:), slope(:), seasonal(:, :)
      real(dp), allocatable :: autoregressive(:), moving_average(:), regression(:)
      real(dp), allocatable :: irregular(:), sigma(:), trend(:), seasonal_total(:), arma(:)
      real(dp), allocatable :: fitted_transformed(:), fitted_original(:)
      real(dp), allocatable :: reconstructed_transformed(:), reconstructed_original(:)
      real(dp), allocatable :: residual_original(:)
      integer :: info = 0
   end type

   type, public :: tsissm_prediction_decomposition_t
      ! Simulated state contributions with simulation in the first dimension.
      real(dp), allocatable :: level(:, :), slope(:, :), seasonal(:, :, :)
      real(dp), allocatable :: autoregressive(:, :), moving_average(:, :), regression(:, :)
      real(dp), allocatable :: irregular(:, :), trend(:, :), seasonal_total(:, :), arma(:, :)
      real(dp), allocatable :: reconstructed_transformed(:, :), reconstructed_original(:, :)
      integer :: info = 0
   end type

   type, public :: tsissm_forecast_t
      ! Simulated forecast distribution, summaries, states, and analytic moments.
      type(tsissm_prediction_t) :: prediction
      type(tsissm_moments_t) :: moments
      real(dp), allocatable :: transformed(:, :), distribution(:, :)
      real(dp), allocatable :: mean(:), probabilities(:), quantile(:, :)
      integer :: info = 0
   end type

   type, public :: tsissm_diagnostics_t
      ! Residual views, serial-dependence tests, normality, and ARCH diagnostics.
      real(dp), allocatable :: raw_residual(:), transformed_residual(:), standardized_residual(:)
      real(dp), allocatable :: residual_acf(:), squared_residual_acf(:)
      type(itsmr_randomness_tests_t) :: randomness
      type(mts_arch_test_t) :: arch
      real(dp) :: weighted_ljung_box = 0.0_dp
      real(dp) :: weighted_ljung_box_p_value = 1.0_dp
      real(dp) :: normality_statistic = 0.0_dp
      real(dp) :: normality_p_value = 1.0_dp
      integer :: lag = 0
      integer :: info = 0
   end type

   type, public :: tsissm_structural_diagnostics_t
      ! Residual diagnostics with weighted tests, roots, stability, and outliers.
      type(tsissm_diagnostics_t) :: residual
      real(dp), allocatable :: weighted_ljung_box(:), weighted_ljung_box_p_value(:)
      real(dp), allocatable :: outlier_statistic(:), outlier_critical(:)
      real(dp), allocatable :: stability_modulus(:)
      complex(dp), allocatable :: ar_inverse_root(:), ma_inverse_root(:)
      integer, allocatable :: weighted_lag(:), outlier_candidate_index(:), outlier_index(:)
      integer :: outlier_count = 0
      integer :: info = 0
   end type

   type, public :: tsissm_accuracy_t
      ! Point and distributional forecast accuracy measures.
      type(accuracy_result_t) :: point
      real(dp) :: mse = 0.0_dp
      real(dp) :: smape = 0.0_dp
      real(dp) :: mslre = 0.0_dp
      real(dp) :: bias = 0.0_dp
      real(dp) :: log_score = 0.0_dp
      real(dp) :: interval_coverage = 0.0_dp
      real(dp) :: interval_score = 0.0_dp
      real(dp) :: crps = 0.0_dp
      integer :: observations = 0
      integer :: info = 0
   end type

   type, public :: tsissm_backtest_t
      ! Rolling-origin forecasts, errors, scores, parameters, and horizon summaries.
      real(dp), allocatable :: forecast(:, :), actual(:, :), error(:, :)
      real(dp), allocatable :: lower(:, :), upper(:, :), crps(:, :), log_score(:, :)
      real(dp), allocatable :: persistence(:, :), coefficients(:, :)
      real(dp), allocatable :: parameters(:, :)
      real(dp), allocatable :: rmse(:), mean_absolute_error(:), bias(:)
      real(dp), allocatable :: interval_coverage(:), interval_score(:)
      real(dp), allocatable :: mean_crps(:), mean_log_score(:)
      logical, allocatable :: valid(:, :), reestimated(:)
      integer, allocatable :: training_start(:), training_end(:)
      integer :: first_origin = 0
      integer :: horizon = 0
      integer :: window_length = 0
      integer :: reestimate_every = 0
      integer :: info = 0
   end type

   type, public :: tsissm_profile_t
      ! Simulation-refit parameter and holdout forecast distributions.
      real(dp), allocatable :: parameter(:, :), true_parameter(:)
      real(dp), allocatable :: actual(:, :), predicted(:, :)
      real(dp), allocatable :: mape(:, :), percent_bias(:, :), mase(:, :), crps(:, :)
      logical, allocatable :: valid(:)
      integer :: simulations = 0
      integer :: successful = 0
      integer :: horizon = 0
      integer :: info = 0
   end type

   type, public :: tsissm_candidate_t
      ! One fitted automatic-selection candidate and its structural specification.
      type(tsissm_model_t) :: model
      type(tsissm_filter_t) :: filter
      type(tsissm_likelihood_t) :: likelihood
      real(dp), allocatable :: parameters(:), persistence(:), coefficients(:), arch(:), garch(:)
      integer, allocatable :: seasonal_harmonics(:)
      real(dp) :: lambda = 1.0_dp
      real(dp) :: skew = 0.0_dp
      real(dp) :: shape = 0.0_dp
      real(dp) :: variance_intercept = 0.0_dp
      real(dp) :: aic = huge(1.0_dp)
      real(dp) :: aicc = huge(1.0_dp)
      real(dp) :: bic = huge(1.0_dp)
      real(dp) :: score = huge(1.0_dp)
      real(dp) :: weight = 0.0_dp
      integer :: ar_order = 0
      integer :: ma_order = 0
      integer :: distribution = tsissm_distribution_gaussian
      integer :: parameter_count = 0
      integer :: info = 0
      logical :: slope = .false.
      logical :: damped_slope = .false.
      logical :: seasonal = .false.
      logical :: dynamic_variance = .false.
      logical :: successful = .false.
   end type

   type, public :: tsissm_selection_t
      ! Full candidate grid with ranked top-model indices and criterion weights.
      type(tsissm_candidate_t), allocatable :: candidates(:)
      integer, allocatable :: top_index(:)
      real(dp), allocatable :: weight(:)
      integer :: criterion = tsissm_selection_aic
      integer :: evaluated = 0
      integer :: successful = 0
      integer :: info = 0
   end type

   type, public :: tsissm_ensemble_t
      ! Selected-model forecasts and their pathwise weighted ensemble.
      type(tsissm_forecast_t), allocatable :: model_forecast(:)
      real(dp), allocatable :: distribution(:, :), mean(:), analytic_mean(:)
      real(dp), allocatable :: probabilities(:), quantile(:, :), weight(:), correlation(:, :)
      real(dp), allocatable :: trend(:, :), seasonal(:, :), arma(:, :)
      real(dp), allocatable :: regression(:, :), irregular(:, :)
      logical :: has_decomposition = .false.
      integer :: info = 0
   end type

   public :: tsissm_filter_constant, tsissm_predict_constant
   public :: tsissm_garch_recursion, tsissm_filter_dynamic, tsissm_predict_dynamic
   public :: tsissm_initialize_variance
   public :: tsissm_model
   public :: tsissm_log_density, tsissm_likelihood
   public :: tsissm_ar_constraint, tsissm_ma_constraint, tsissm_stability_constraint
   public :: tsissm_parameter_constraints
   public :: tsissm_initialize_states, tsissm_fit_constant
   public :: tsissm_fit_dynamic
   public :: tsissm_fit_structural_constant
   public :: tsissm_fit_structural_dynamic
   public :: tsissm_moments_constant, tsissm_moments_dynamic
   public :: tsissm_hresiduals_constant, tsissm_hresiduals_dynamic
   public :: tsissm_fd_scores, tsissm_parameter_covariance
   public :: tsissm_covariance_constant, tsissm_covariance_dynamic
   public :: tsissm_covariance_structural
   public :: tsissm_covariance_structural_dynamic
   public :: tsissm_decompose_filter, tsissm_decompose_prediction
   public :: tsissm_simulate_constant, tsissm_simulate_dynamic
   public :: tsissm_profile_constant, tsissm_profile_dynamic
   public :: tsissm_profile_structural_constant, tsissm_profile_structural_dynamic
   public :: tsissm_diagnose, tsissm_forecast_accuracy
   public :: tsissm_diagnose_structural
   public :: tsissm_weighted_box_test
   public :: tsissm_weighted_box_test_t
   public :: tsissm_box_test_box_pierce, tsissm_box_test_ljung_box
   public :: tsissm_box_test_monti, tsissm_residual_raw
   public :: tsissm_residual_squared, tsissm_residual_log_squared
   public :: tsissm_residual_absolute
   public :: tsissm_backtest_constant, tsissm_backtest_dynamic
   public :: tsissm_backtest_structural_constant, tsissm_backtest_structural_dynamic
   public :: tsissm_select_models
   public :: tsissm_rerank_models
   public :: tsissm_selection_correlation, tsissm_ensemble_forecast

contains

   pure function tsissm_initialize_variance(innovations, observed, method, sample_size) result(out)
      !! Initialize GARCH variance from all or the first selected valid innovations.
      real(dp), intent(in) :: innovations(:) !! Model innovations.
      logical, intent(in), optional :: observed(:) !! Flag controlling observed.
      integer, intent(in), optional :: method !! Algorithm or estimation method.
      integer, intent(in), optional :: sample_size !! Sample size.
      type(tsissm_variance_initialization_t) :: out
      logical, allocatable :: use_observation(:)
      real(dp), allocatable :: selected(:)
      integer :: requested

      out%method = tsissm_variance_initial_full
      if (present(method)) out%method = method
      if (size(innovations) < 1 .or. &
         (out%method /= tsissm_variance_initial_full .and. &
         out%method /= tsissm_variance_initial_sample)) then
         out%info = 1
         return
      end if
      allocate(use_observation(size(innovations)))
      use_observation = .true.
      if (present(observed)) then
         if (size(observed) /= size(innovations)) then
            out%info = 1
            return
         end if
         use_observation = observed
      end if
      selected = pack(innovations, use_observation)
      if (size(selected) < 1) then
         out%info = 1
         return
      end if
      requested = size(selected)
      if (out%method == tsissm_variance_initial_sample) then
         if (.not. present(sample_size)) then
            out%info = 1
            return
         end if
         requested = sample_size
         if (requested < 1 .or. requested > size(selected)) then
            out%info = 1
            return
         end if
      end if
      out%observations = requested
      out%variance = sum(selected(:requested)**2)/real(requested, dp)
      if (.not. ieee_is_finite(out%variance) .or. out%variance < 0.0_dp) out%info = 2
   end function tsissm_initialize_variance

   pure function tsissm_fd_scores(contribution_minus, contribution_plus, step) result(scores)
      !! Differentiate per-observation negative log-likelihood contributions.
      real(dp), intent(in) :: contribution_minus(:, :) !! Contribution minus.
      real(dp), intent(in) :: contribution_plus(:, :) !! Contribution plus.
      real(dp), intent(in) :: step(:) !! Step.
      real(dp), allocatable :: scores(:, :)
      integer :: i

      if (any(shape(contribution_minus) /= shape(contribution_plus)) .or. &
         size(contribution_minus, 2) /= size(step) .or. any(step <= 0.0_dp)) then
         allocate(scores(0, 0))
         return
      end if
      allocate(scores(size(contribution_minus, 1), size(contribution_minus, 2)))
      do i = 1, size(step)
         scores(:, i) = (contribution_plus(:, i) - contribution_minus(:, i))/(2.0_dp*step(i))
      end do
   end function tsissm_fd_scores

   pure function tsissm_parameter_covariance(scores, hessian, covariance_type, adjust, &
      newey_west_lag, hac_kernel, prewhite_order, hac_bandwidth) result(out)
      !! Estimate Hessian, OPG, QMLE, or Bartlett Newey-West parameter covariance.
      real(dp), intent(in) :: scores(:, :) !! Scores.
      real(dp), intent(in) :: hessian(:, :) !! Hessian.
      integer, intent(in) :: covariance_type !! Covariance type.
      logical, intent(in), optional :: adjust !! Flag controlling adjust.
      integer, intent(in), optional :: newey_west_lag !! Newey west lag.
      integer, intent(in), optional :: hac_kernel !! Hac kernel.
      integer, intent(in), optional :: prewhite_order !! Prewhite order.
      real(dp), intent(in), optional :: hac_bandwidth !! Hac bandwidth.
      type(tsissm_covariance_t) :: out
      real(dp), allocatable :: inverse_meat(:, :)
      real(dp) :: correction
      integer :: i, inversion_info, n, parameters

      n = size(scores, 1)
      parameters = size(scores, 2)
      out%covariance_type = covariance_type
      out%observations = n
      out%parameters = parameters
      if (present(adjust)) out%adjusted = adjust
      if (n < 1 .or. parameters < 1 .or. any(shape(hessian) /= [parameters, parameters]) .or. &
         covariance_type < tsissm_covariance_hessian .or. &
         covariance_type > tsissm_covariance_newey_west) then
         out%info = 1
         return
      end if
      if (out%adjusted .and. n <= parameters) then
         out%info = 1
         return
      end if
      call invert_matrix(hessian, out%bread, inversion_info)
      if (inversion_info /= 0) then
         out%info = 2
         return
      end if
      out%bread = 0.5_dp*(out%bread + transpose(out%bread))
      correction = 1.0_dp
      if (out%adjusted) correction = real(n, dp)/real(n - parameters, dp)
      select case (covariance_type)
      case (tsissm_covariance_hessian)
         allocate(out%meat(parameters, parameters))
         out%meat = 0.0_dp
         out%covariance = out%bread
      case (tsissm_covariance_opg)
         out%meat = matmul(transpose(scores), scores)
         call invert_matrix(out%meat, inverse_meat, inversion_info)
         if (inversion_info /= 0) then
            out%info = 3
            return
         end if
         out%covariance = correction*inverse_meat
      case (tsissm_covariance_qmle)
         out%meat = correction*matmul(transpose(scores), scores)
         out%covariance = matmul(matmul(out%bread, out%meat), out%bread)
      case (tsissm_covariance_newey_west)
         if (present(hac_kernel)) out%kernel = hac_kernel
         if (present(prewhite_order)) out%prewhite_order = prewhite_order
         call hac_score_meat(scores, out%kernel, out%prewhite_order, newey_west_lag, &
            hac_bandwidth, out%meat, out%lag, out%bandwidth, out%bias_correction, &
            out%effective_degrees, inversion_info)
         if (inversion_info /= 0) then
            out%info = 4 + inversion_info
            return
         end if
         out%meat = correction*out%meat
         out%covariance = matmul(matmul(out%bread, out%meat), out%bread)
      end select
      out%covariance = 0.5_dp*(out%covariance + transpose(out%covariance))
      allocate(out%standard_error(parameters), out%correlation(parameters, parameters))
      out%standard_error = sqrt(max([(out%covariance(i, i), i=1, parameters)], 0.0_dp))
      out%correlation = 0.0_dp
      do i = 1, parameters
         if (out%standard_error(i) > 0.0_dp) then
            out%correlation(i, :) = out%covariance(i, :)/ &
               (out%standard_error(i)*max(out%standard_error, tiny(1.0_dp)))
         end if
      end do
      out%correlation = max(-1.0_dp, min(1.0_dp, out%correlation))
   end function tsissm_parameter_covariance

   pure subroutine hac_score_meat(scores, kernel, prewhite_order, requested_lag, &
      requested_bandwidth, meat, lag, bandwidth, bias_correction, effective_degrees, info)
      !! Form a kernel HAC score covariance with optional VAR prewhitening.
      real(dp), intent(in) :: scores(:, :) !! Scores.
      integer, intent(in) :: kernel !! Kernel.
      integer, intent(in) :: prewhite_order !! Prewhite order.
      integer, intent(in), optional :: requested_lag !! Requested lag.
      real(dp), intent(in), optional :: requested_bandwidth !! Requested bandwidth.
      real(dp), allocatable, intent(out) :: meat(:, :) !! Meat.
      integer, intent(out) :: lag !! Lag index or number of lags.
      integer, intent(out) :: info !! Status code; zero indicates success.
      real(dp), intent(out) :: bandwidth !! Smoothing or spectral bandwidth.
      real(dp), intent(out) :: bias_correction !! Bias correction.
      real(dp), intent(out) :: effective_degrees !! Effective degrees.
      real(dp), allocatable :: coefficients(:, :), design(:, :), inverse(:, :)
      real(dp), allocatable :: recolor(:, :), response(:, :), working(:, :), xtx(:, :)
      real(dp), allocatable :: ar_sum(:, :), identity(:, :)
      real(dp) :: weight, weight_sum, weight_square_sum
      integer :: block, i, k, m, n, status

      n = size(scores, 1)
      k = size(scores, 2)
      info = 0
      if (kernel < tsissm_hac_bartlett .or. kernel > tsissm_hac_tukey_hanning .or. &
         prewhite_order < 0 .or. prewhite_order >= n .or. n < 1 .or. k < 1) then
         info = 1
         return
      end if
      working = scores
      allocate(identity(k, k))
      identity = 0.0_dp
      do i = 1, k
         identity(i, i) = 1.0_dp
      end do
      recolor = identity
      if (prewhite_order > 0) then
         m = n - prewhite_order
         if (m <= k*prewhite_order) then
            info = 1
            return
         end if
         allocate(design(m, k*prewhite_order), response(m, k))
         response = scores(prewhite_order + 1:, :)
         do block = 1, prewhite_order
            design(:, (block - 1)*k + 1:block*k) = &
               scores(prewhite_order - block + 1:n - block, :)
         end do
         xtx = matmul(transpose(design), design)
         call invert_matrix(xtx, inverse, status)
         if (status /= 0) then
            info = 2
            return
         end if
         coefficients = matmul(matmul(inverse, transpose(design)), response)
         working = response - matmul(design, coefficients)
         allocate(ar_sum(k, k))
         ar_sum = 0.0_dp
         do block = 1, prewhite_order
            ar_sum = ar_sum + transpose(coefficients((block - 1)*k + 1:block*k, :))
         end do
         call invert_matrix(identity - ar_sum, recolor, status)
         if (status /= 0) then
            info = 3
            return
         end if
      end if
      m = size(working, 1)
      if (present(requested_lag)) then
         lag = requested_lag
         if (lag < 0 .or. lag >= m) then
            info = 1
            return
         end if
         bandwidth = real(lag + 1, dp)
      else
         bandwidth = automatic_hac_bandwidth(working, kernel)
         if (present(requested_bandwidth)) bandwidth = requested_bandwidth
         if (bandwidth <= 0.0_dp) then
            info = 1
            return
         end if
         if (kernel == tsissm_hac_quadratic_spectral) then
            lag = m - 1
         else
            lag = min(m - 1, max(0, ceiling(bandwidth) - 1))
         end if
      end if
      meat = matmul(transpose(working), working)
      weight_sum = real(m, dp)
      weight_square_sum = real(m, dp)
      do i = 1, lag
         weight = hac_kernel_weight(real(i, dp)/bandwidth, kernel)
         meat = meat + weight*(matmul(transpose(working(1:m - i, :)), &
            working(i + 1:m, :)) + matmul(transpose(working(i + 1:m, :)), &
            working(1:m - i, :)))
         weight_sum = weight_sum + 2.0_dp*real(m - i, dp)*weight
         weight_square_sum = weight_square_sum + 2.0_dp*real(m - i, dp)*weight**2
      end do
      if (prewhite_order > 0) meat = matmul(matmul(recolor, meat), transpose(recolor))
      bias_correction = real(m, dp)**2/max(real(m, dp)**2 - weight_sum, tiny(1.0_dp))
      effective_degrees = real(m, dp)**2/max(weight_square_sum, tiny(1.0_dp))
   end subroutine hac_score_meat

   pure real(dp) function automatic_hac_bandwidth(scores, kernel) result(bandwidth)
      !! Estimate a plug-in HAC bandwidth from aggregate score autocovariances.
      real(dp), intent(in) :: scores(:, :) !! Scores.
      integer, intent(in) :: kernel !! Kernel.
      real(dp), allocatable :: series(:)
      real(dp) :: gamma, long_run, moment, ratio
      integer :: i, pilot

      series = sum(scores, dim=2)
      series = series - sum(series)/real(size(series), dp)
      pilot = min(size(series) - 1, max(1, floor(sqrt(real(size(series), dp)))))
      long_run = dot_product(series, series)/real(size(series), dp)
      moment = 0.0_dp
      do i = 1, pilot
         gamma = dot_product(series(1:size(series) - i), series(i + 1:))/ &
            real(size(series), dp)
         long_run = long_run + 2.0_dp*gamma
         if (kernel == tsissm_hac_parzen .or. &
            kernel == tsissm_hac_quadratic_spectral) then
            moment = moment + 2.0_dp*real(i*i, dp)*gamma
         else
            moment = moment + 2.0_dp*real(i, dp)*gamma
         end if
      end do
      ratio = (moment/max(abs(long_run), sqrt(epsilon(1.0_dp))))**2
      select case (kernel)
      case (tsissm_hac_parzen)
         bandwidth = 2.6614_dp*(ratio*real(size(series), dp))**0.2_dp
      case (tsissm_hac_quadratic_spectral)
         bandwidth = 1.3221_dp*(ratio*real(size(series), dp))**0.2_dp
      case default
         bandwidth = 1.1447_dp*(ratio*real(size(series), dp))**(1.0_dp/3.0_dp)
      end select
      bandwidth = max(1.0_dp, bandwidth)
   end function automatic_hac_bandwidth

   pure elemental real(dp) function hac_kernel_weight(argument, kernel) result(weight)
      !! Evaluate a supported HAC kernel at a nonnegative scaled lag.
      real(dp), intent(in) :: argument !! Argument.
      integer, intent(in) :: kernel !! Kernel.
      real(dp) :: scaled

      select case (kernel)
      case (tsissm_hac_bartlett)
         weight = max(0.0_dp, 1.0_dp - argument)
      case (tsissm_hac_parzen)
         if (argument <= 0.5_dp) then
            weight = 1.0_dp - 6.0_dp*argument**2 + 6.0_dp*argument**3
         else if (argument <= 1.0_dp) then
            weight = 2.0_dp*(1.0_dp - argument)**3
         else
            weight = 0.0_dp
         end if
      case (tsissm_hac_quadratic_spectral)
         if (argument <= sqrt(epsilon(1.0_dp))) then
            weight = 1.0_dp
         else
            scaled = 6.0_dp*acos(-1.0_dp)*argument/5.0_dp
            weight = 25.0_dp/(12.0_dp*acos(-1.0_dp)**2*argument**2)* &
               (sin(scaled)/scaled - cos(scaled))
         end if
      case (tsissm_hac_truncated)
         weight = merge(1.0_dp, 0.0_dp, argument <= 1.0_dp)
      case (tsissm_hac_tukey_hanning)
         weight = merge(0.5_dp*(1.0_dp + cos(acos(-1.0_dp)*argument)), &
            0.0_dp, argument <= 1.0_dp)
      case default
         weight = 0.0_dp
      end select
   end function hac_kernel_weight

   pure function tsissm_covariance_constant(observations, filtered, regressors, persistence, &
      coefficients, lambda, distribution, skew, shape, covariance_type, observed, adjust, &
      newey_west_lag, difference_step, hac_kernel, prewhite_order, hac_bandwidth) result(out)
      !! Estimate parameter covariance for a fitted constant-variance ISSM.
      real(dp), intent(in) :: observations(:) !! Observed time-series values.
      real(dp), intent(in) :: regressors(:, :) !! Regression design matrix.
      real(dp), intent(in) :: persistence(:) !! Persistence.
      real(dp), intent(in) :: coefficients(:) !! Model coefficients.
      type(tsissm_filter_t), intent(in) :: filtered !! Filtered.
      real(dp), intent(in) :: lambda !! Penalty or shrinkage parameter.
      real(dp), intent(in) :: skew !! Skew.
      real(dp), intent(in) :: shape !! Shape.
      integer, intent(in) :: distribution !! Probability-distribution specification.
      integer, intent(in) :: covariance_type !! Covariance type.
      logical, intent(in), optional :: observed(:) !! Flag controlling observed.
      logical, intent(in), optional :: adjust !! Flag controlling adjust.
      integer, intent(in), optional :: newey_west_lag !! Newey west lag.
      real(dp), intent(in), optional :: difference_step !! Difference step.
      integer, intent(in), optional :: hac_kernel !! Hac kernel.
      integer, intent(in), optional :: prewhite_order !! Prewhite order.
      real(dp), intent(in), optional :: hac_bandwidth !! Hac bandwidth.
      type(tsissm_covariance_t) :: out
      logical, allocatable :: use_observation(:)
      real(dp), allocatable :: all_scores(:, :), minus(:, :), parameters(:), plus(:, :), steps(:)
      real(dp), allocatable :: hessian(:, :), identity(:, :)
      real(dp) :: step
      integer :: i, parameter_count

      parameter_count = size(persistence) + size(coefficients)
      if (.not. covariance_fit_inputs_valid(observations, filtered, regressors, coefficients) .or. &
         size(persistence) /= size(filtered%persistence) .or. parameter_count < 1) then
         out%info = 1
         return
      end if
      allocate(use_observation(size(observations)))
      use_observation = .true.
      if (present(observed)) then
         if (size(observed) /= size(observations)) then
            out%info = 1
            return
         end if
         use_observation = observed
      end if
      if (count(use_observation) <= parameter_count) then
         out%info = 1
         return
      end if
      allocate(identity(size(persistence), size(persistence)))
      identity = 1.0_dp
      parameters = [persistence, coefficients]
      step = epsilon(1.0_dp)**0.25_dp
      if (present(difference_step)) step = difference_step
      if (step <= 0.0_dp) then
         out%info = 1
         return
      end if
      allocate(minus(size(observations), parameter_count), plus(size(observations), parameter_count))
      allocate(steps(parameter_count))
      do i = 1, parameter_count
         steps(i) = step*max(1.0_dp, abs(parameters(i)))
         parameters(i) = parameters(i) - steps(i)
         minus(:, i) = contributions(parameters)
         parameters(i) = parameters(i) + 2.0_dp*steps(i)
         plus(:, i) = contributions(parameters)
         parameters(i) = parameters(i) - steps(i)
      end do
      if (.not. score_perturbations_valid(minus, plus)) then
         out%info = 4
         return
      end if
      all_scores = tsissm_fd_scores(minus, plus, steps)
      hessian = finite_difference_hessian(objective, parameters, step)
      out = tsissm_parameter_covariance(pack_rows(all_scores, use_observation), hessian, &
         covariance_type, adjust, newey_west_lag, hac_kernel, prewhite_order, hac_bandwidth)
      out%scores = pack_rows(all_scores, use_observation)
      out%hessian = hessian

   contains

      pure function contributions(trial_parameters) result(value)
         !! Return negative log-likelihood contributions at trial parameters.
         real(dp), intent(in) :: trial_parameters(:) !! Trial parameters.
         real(dp) :: value(size(observations))
         type(tsissm_filter_t) :: trial_filter
         type(tsissm_likelihood_t) :: trial_likelihood

         trial_filter = tsissm_filter_constant(observations, filtered%transition, identity, identity, &
            filtered%observation_loading, trial_parameters(:size(persistence)), regressors, &
            trial_parameters(size(persistence) + 1:), filtered%state(:, 0), use_observation, lambda)
         if (trial_filter%info /= 0) then
            value = huge(1.0_dp)
            return
         end if
         trial_likelihood = tsissm_likelihood(trial_filter%innovation, observations, lambda, &
            distribution, skew, shape, parameter_count, use_observation)
         if (trial_likelihood%info /= 0) then
            value = huge(1.0_dp)
         else
            value = trial_likelihood%contribution
         end if
      end function contributions

      pure real(dp) function objective(trial_parameters) result(value)
         !! Sum negative log-likelihood contributions at trial parameters.
         real(dp), intent(in) :: trial_parameters(:) !! Trial parameters.

         value = sum(contributions(trial_parameters))
      end function objective
   end function tsissm_covariance_constant

   pure function tsissm_covariance_dynamic(observations, filtered, regressors, persistence, &
      coefficients, arch, garch, lambda, distribution, skew, shape, covariance_type, &
      observed, adjust, newey_west_lag, difference_step, variance_initialization, &
      variance_sample_size, hac_kernel, prewhite_order, hac_bandwidth) result(out)
      !! Estimate parameter covariance for a fitted variance-targeted GARCH ISSM.
      real(dp), intent(in) :: observations(:) !! Observed time-series values.
      real(dp), intent(in) :: regressors(:, :) !! Regression design matrix.
      real(dp), intent(in) :: persistence(:) !! Persistence.
      real(dp), intent(in) :: coefficients(:) !! Model coefficients.
      real(dp), intent(in) :: arch(:) !! Arch.
      real(dp), intent(in) :: garch(:) !! Garch.
      real(dp), intent(in) :: lambda !! Penalty or shrinkage parameter.
      real(dp), intent(in) :: skew !! Skew.
      real(dp), intent(in) :: shape !! Shape.
      type(tsissm_filter_t), intent(in) :: filtered !! Filtered.
      integer, intent(in) :: distribution !! Probability-distribution specification.
      integer, intent(in) :: covariance_type !! Covariance type.
      logical, intent(in), optional :: observed(:) !! Flag controlling observed.
      logical, intent(in), optional :: adjust !! Flag controlling adjust.
      integer, intent(in), optional :: newey_west_lag !! Newey west lag.
      real(dp), intent(in), optional :: difference_step !! Difference step.
      integer, intent(in), optional :: variance_initialization !! Variance initialization.
      integer, intent(in), optional :: variance_sample_size !! Variance sample size.
      integer, intent(in), optional :: hac_kernel !! Hac kernel.
      integer, intent(in), optional :: prewhite_order !! Prewhite order.
      real(dp), intent(in), optional :: hac_bandwidth !! Hac bandwidth.
      type(tsissm_covariance_t) :: out
      logical, allocatable :: use_observation(:)
      real(dp), allocatable :: all_scores(:, :), minus(:, :), parameters(:), plus(:, :), steps(:)
      real(dp), allocatable :: hessian(:, :), identity(:, :)
      real(dp) :: step
      integer :: i, parameter_count

      parameter_count = size(persistence) + size(coefficients) + size(arch) + size(garch)
      if (.not. covariance_fit_inputs_valid(observations, filtered, regressors, coefficients) .or. &
         size(persistence) /= size(filtered%persistence) .or. max(size(arch), size(garch)) < 1 .or. &
         any(arch < 0.0_dp) .or. any(garch < 0.0_dp)) then
         out%info = 1
         return
      end if
      allocate(use_observation(size(observations)))
      use_observation = .true.
      if (present(observed)) then
         if (size(observed) /= size(observations)) then
            out%info = 1
            return
         end if
         use_observation = observed
      end if
      if (count(use_observation) <= parameter_count) then
         out%info = 1
         return
      end if
      allocate(identity(size(persistence), size(persistence)))
      identity = 1.0_dp
      parameters = [persistence, coefficients, arch, garch]
      step = epsilon(1.0_dp)**0.25_dp
      if (present(difference_step)) step = difference_step
      if (step <= 0.0_dp) then
         out%info = 1
         return
      end if
      allocate(minus(size(observations), parameter_count), plus(size(observations), parameter_count))
      allocate(steps(parameter_count))
      do i = 1, parameter_count
         steps(i) = step*max(1.0_dp, abs(parameters(i)))
         parameters(i) = parameters(i) - steps(i)
         minus(:, i) = contributions(parameters)
         parameters(i) = parameters(i) + 2.0_dp*steps(i)
         plus(:, i) = contributions(parameters)
         parameters(i) = parameters(i) - steps(i)
      end do
      if (.not. score_perturbations_valid(minus, plus)) then
         out%info = 4
         return
      end if
      all_scores = tsissm_fd_scores(minus, plus, steps)
      hessian = finite_difference_hessian(objective, parameters, step)
      out = tsissm_parameter_covariance(pack_rows(all_scores, use_observation), hessian, &
         covariance_type, adjust, newey_west_lag, hac_kernel, prewhite_order, hac_bandwidth)
      out%scores = pack_rows(all_scores, use_observation)
      out%hessian = hessian

   contains

      pure function contributions(trial_parameters) result(value)
         !! Return dynamic negative log-likelihood contributions at trial parameters.
         real(dp), intent(in) :: trial_parameters(:) !! Trial parameters.
         real(dp) :: value(size(observations))
         type(tsissm_filter_t) :: constant_filter, trial_filter
         type(tsissm_likelihood_t) :: trial_likelihood
         type(tsissm_variance_initialization_t) :: initialized_variance
         real(dp), allocatable :: initial_arch(:), initial_variance(:)
         real(dp) :: intercept, variance
         integer :: arch_start, garch_start, history

         arch_start = size(persistence) + size(coefficients) + 1
         garch_start = arch_start + size(arch)
         history = max(size(arch), size(garch))
         constant_filter = tsissm_filter_constant(observations, filtered%transition, identity, &
            identity, filtered%observation_loading, trial_parameters(:size(persistence)), &
            regressors, trial_parameters(size(persistence) + 1:arch_start - 1), &
            filtered%state(:, 0), use_observation, lambda)
         if (constant_filter%info /= 0) then
            value = huge(1.0_dp)
            return
         end if
         initialized_variance = tsissm_initialize_variance(constant_filter%innovation, &
            use_observation, variance_initialization, variance_sample_size)
         if (initialized_variance%info /= 0) then
            value = huge(1.0_dp)
            return
         end if
         variance = initialized_variance%variance
         intercept = max(1.0e-12_dp, variance*(1.0_dp - sum(trial_parameters(arch_start:))))
         allocate(initial_arch(history), initial_variance(history))
         initial_arch = variance
         initial_variance = variance
         trial_filter = tsissm_filter_dynamic(observations, filtered%transition, identity, identity, &
            filtered%observation_loading, trial_parameters(:size(persistence)), regressors, &
            trial_parameters(size(persistence) + 1:arch_start - 1), filtered%state(:, 0), &
            use_observation, lambda, trial_parameters(arch_start:garch_start - 1), &
            trial_parameters(garch_start:), initial_arch, initial_variance, intercept)
         if (trial_filter%info /= 0) then
            value = huge(1.0_dp)
            return
         end if
         trial_likelihood = tsissm_likelihood(trial_filter%innovation, observations, lambda, &
            distribution, skew, shape, parameter_count, use_observation, trial_filter%conditional_sd)
         if (trial_likelihood%info /= 0) then
            value = huge(1.0_dp)
         else
            value = trial_likelihood%contribution
         end if
      end function contributions

      pure real(dp) function objective(trial_parameters) result(value)
         !! Sum dynamic negative log-likelihood contributions at trial parameters.
         real(dp), intent(in) :: trial_parameters(:) !! Trial parameters.

         value = sum(contributions(trial_parameters))
      end function objective
   end function tsissm_covariance_dynamic

   pure function tsissm_covariance_structural(observations, slope, damped_slope, &
      seasonal_frequency, seasonal_harmonics, regular_seasonal, ar_order, ma_order, &
      regressors, initial_state, parameters, distribution, covariance_type, observed, &
      adjust, newey_west_lag, difference_step, hac_kernel, prewhite_order, &
      hac_bandwidth) result(out)
      !! Estimate joint structural ISSM parameter covariance from likelihood scores.
      real(dp), intent(in) :: observations(:) !! Observed time-series values.
      real(dp), intent(in) :: seasonal_frequency(:) !! Seasonal frequency.
      real(dp), intent(in) :: regressors(:, :) !! Regression design matrix.
      real(dp), intent(in) :: initial_state(:) !! Initial state vector.
      real(dp), intent(in) :: parameters(:) !! Model parameter values.
      logical, intent(in) :: slope !! Flag controlling slope.
      logical, intent(in) :: damped_slope !! Flag controlling damped slope.
      logical, intent(in) :: regular_seasonal !! Flag controlling regular seasonal.
      integer, intent(in) :: seasonal_harmonics(:) !! Seasonal harmonics.
      integer, intent(in) :: ar_order !! Autoregressive order.
      integer, intent(in) :: ma_order !! Moving-average order.
      integer, intent(in) :: distribution !! Probability-distribution specification.
      integer, intent(in) :: covariance_type !! Covariance type.
      logical, intent(in), optional :: observed(:) !! Flag controlling observed.
      logical, intent(in), optional :: adjust !! Flag controlling adjust.
      integer, intent(in), optional :: newey_west_lag !! Newey west lag.
      real(dp), intent(in), optional :: difference_step !! Difference step.
      integer, intent(in), optional :: hac_kernel !! Hac kernel.
      integer, intent(in), optional :: prewhite_order !! Prewhite order.
      real(dp), intent(in), optional :: hac_bandwidth !! Hac bandwidth.
      type(tsissm_covariance_t) :: out
      logical, allocatable :: use_observation(:)
      real(dp), allocatable :: all_scores(:, :), minus(:, :), plus(:, :), steps(:)
      real(dp), allocatable :: hessian(:, :), trial_parameters(:)
      real(dp) :: step
      integer :: distribution_count, expected, i

      distribution_count = 0
      if (distribution == tsissm_distribution_student) distribution_count = 1
      if (distribution == tsissm_distribution_johnson_su) distribution_count = 2
      expected = 1 + merge(1, 0, slope) + merge(1, 0, damped_slope) + &
         size(seasonal_frequency) + ar_order + ma_order + size(regressors, 2) + &
         1 + distribution_count
      if (distribution /= tsissm_distribution_gaussian .and. &
         distribution /= tsissm_distribution_student .and. &
         distribution /= tsissm_distribution_johnson_su) then
         out%info = 1
         return
      end if
      if (size(parameters) /= expected .or. size(observations) < 1 .or. &
         size(regressors, 1) /= size(observations) .or. any(observations <= 0.0_dp) .or. &
         ar_order < 0 .or. ma_order < 0 .or. size(seasonal_harmonics) /= &
         size(seasonal_frequency)) then
         out%info = 1
         return
      end if
      if (.not. valid_distribution_parameters(distribution, &
         structural_shape(parameters, distribution))) then
         out%info = 1
         return
      end if
      allocate(use_observation(size(observations)))
      use_observation = .true.
      if (present(observed)) then
         if (size(observed) /= size(observations)) then
            out%info = 1
            return
         end if
         use_observation = observed
      end if
      if (count(use_observation) <= expected) then
         out%info = 1
         return
      end if
      step = epsilon(1.0_dp)**0.25_dp
      if (present(difference_step)) step = difference_step
      if (step <= 0.0_dp) then
         out%info = 1
         return
      end if
      allocate(minus(size(observations), expected), plus(size(observations), expected))
      allocate(steps(expected))
      trial_parameters = parameters
      do i = 1, expected
         steps(i) = step*max(1.0_dp, abs(parameters(i)))
         trial_parameters(i) = parameters(i) - steps(i)
         minus(:, i) = contributions(trial_parameters)
         trial_parameters(i) = parameters(i) + steps(i)
         plus(:, i) = contributions(trial_parameters)
         trial_parameters(i) = parameters(i)
      end do
      if (.not. score_perturbations_valid(minus, plus)) then
         out%info = 4
         return
      end if
      all_scores = tsissm_fd_scores(minus, plus, steps)
      hessian = finite_difference_hessian(objective, parameters, step)
      out = tsissm_parameter_covariance(pack_rows(all_scores, use_observation), hessian, &
         covariance_type, adjust, newey_west_lag, hac_kernel, prewhite_order, hac_bandwidth)
      out%scores = pack_rows(all_scores, use_observation)
      out%hessian = hessian

   contains

      pure function contributions(trial) result(value)
         !! Decode structural parameters and return likelihood contributions.
         real(dp), intent(in) :: trial(:) !! Trial.
         real(dp) :: value(size(observations))
         type(tsissm_filter_t) :: filtered
         type(tsissm_likelihood_t) :: likelihood
         type(tsissm_model_t) :: model
         real(dp), allocatable :: ar(:), coefficients(:), ma(:), seasonal_persistence(:)
         real(dp) :: alpha, beta, damping, lambda, shape, skew
         integer :: index

         index = 1
         alpha = trial(index)
         index = index + 1
         beta = 0.0_dp
         if (slope) then
            beta = trial(index)
            index = index + 1
         end if
         damping = 1.0_dp
         if (damped_slope) then
            damping = trial(index)
            index = index + 1
         end if
         seasonal_persistence = trial(index:index + size(seasonal_frequency) - 1)
         index = index + size(seasonal_frequency)
         ar = trial(index:index + ar_order - 1)
         index = index + ar_order
         ma = trial(index:index + ma_order - 1)
         index = index + ma_order
         coefficients = trial(index:index + size(regressors, 2) - 1)
         index = index + size(regressors, 2)
         lambda = trial(index)
         index = index + 1
         skew = 0.0_dp
         shape = 0.0_dp
         if (distribution == tsissm_distribution_student) then
            shape = trial(index)
         else if (distribution == tsissm_distribution_johnson_su) then
            skew = trial(index)
            shape = trial(index + 1)
         end if
         model = tsissm_model(slope, damped_slope, seasonal_frequency, seasonal_harmonics, &
            regular_seasonal, ar, ma, alpha, beta, seasonal_persistence, damping)
         if (model%info /= 0) then
            value = huge(1.0_dp)
            return
         end if
         if (size(model%initial_state) /= size(initial_state)) then
            value = huge(1.0_dp)
            return
         end if
         filtered = tsissm_filter_constant(observations, model%transition_base, &
            model%transition_scale, model%transition_parameter, model%observation_loading, &
            model%persistence, regressors, coefficients, initial_state, use_observation, lambda)
         if (filtered%info /= 0) then
            value = huge(1.0_dp)
            return
         end if
         likelihood = tsissm_likelihood(filtered%innovation, observations, lambda, distribution, &
            skew, shape, size(trial), use_observation)
         if (likelihood%info /= 0) then
            value = huge(1.0_dp)
         else
            value = likelihood%contribution
         end if
      end function contributions

      pure real(dp) function objective(trial) result(value)
         !! Sum joint structural negative log-likelihood contributions.
         real(dp), intent(in) :: trial(:) !! Trial.

         value = sum(contributions(trial))
      end function objective
   end function tsissm_covariance_structural

   pure function tsissm_covariance_structural_dynamic(observations, slope, damped_slope, &
      seasonal_frequency, seasonal_harmonics, regular_seasonal, ar_order, ma_order, &
      regressors, initial_state, parameters, arch_order, garch_order, distribution, &
      covariance_type, observed, adjust, newey_west_lag, difference_step, &
      variance_initialization, variance_sample_size, hac_kernel, prewhite_order, &
      hac_bandwidth) result(out)
      !! Estimate full joint structural and GARCH parameter covariance.
      real(dp), intent(in) :: observations(:) !! Observed time-series values.
      real(dp), intent(in) :: seasonal_frequency(:) !! Seasonal frequency.
      real(dp), intent(in) :: regressors(:, :) !! Regression design matrix.
      real(dp), intent(in) :: initial_state(:) !! Initial state vector.
      real(dp), intent(in) :: parameters(:) !! Model parameter values.
      logical, intent(in) :: slope !! Flag controlling slope.
      logical, intent(in) :: damped_slope !! Flag controlling damped slope.
      logical, intent(in) :: regular_seasonal !! Flag controlling regular seasonal.
      integer, intent(in) :: seasonal_harmonics(:) !! Seasonal harmonics.
      integer, intent(in) :: ar_order !! Autoregressive order.
      integer, intent(in) :: ma_order !! Moving-average order.
      integer, intent(in) :: arch_order !! Arch order.
      integer, intent(in) :: garch_order !! Garch order.
      integer, intent(in) :: distribution !! Probability-distribution specification.
      integer, intent(in) :: covariance_type !! Covariance type.
      logical, intent(in), optional :: observed(:) !! Flag controlling observed.
      logical, intent(in), optional :: adjust !! Flag controlling adjust.
      integer, intent(in), optional :: newey_west_lag !! Newey west lag.
      real(dp), intent(in), optional :: difference_step !! Difference step.
      integer, intent(in), optional :: variance_initialization !! Variance initialization.
      integer, intent(in), optional :: variance_sample_size !! Variance sample size.
      integer, intent(in), optional :: hac_kernel !! Hac kernel.
      integer, intent(in), optional :: prewhite_order !! Prewhite order.
      real(dp), intent(in), optional :: hac_bandwidth !! Hac bandwidth.
      type(tsissm_covariance_t) :: out
      logical, allocatable :: use_observation(:)
      real(dp), allocatable :: all_scores(:, :), minus(:, :), plus(:, :), steps(:)
      real(dp), allocatable :: hessian(:, :), trial(:)
      real(dp) :: step
      integer :: distribution_count, i, structural_count, total_count

      distribution_count = 0
      if (distribution == tsissm_distribution_student) distribution_count = 1
      if (distribution == tsissm_distribution_johnson_su) distribution_count = 2
      structural_count = 1 + merge(1, 0, slope) + merge(1, 0, damped_slope) + &
         size(seasonal_frequency) + ar_order + ma_order + size(regressors, 2) + &
         1 + distribution_count
      total_count = structural_count + arch_order + garch_order
      if (distribution /= tsissm_distribution_gaussian .and. &
         distribution /= tsissm_distribution_student .and. &
         distribution /= tsissm_distribution_johnson_su) then
         out%info = 1
         return
      end if
      if (size(parameters) /= total_count .or. max(arch_order, garch_order) < 1 .or. &
         size(observations) < 1 .or. size(regressors, 1) /= size(observations) .or. &
         any(observations <= 0.0_dp) .or. ar_order < 0 .or. ma_order < 0 .or. &
         arch_order < 0 .or. garch_order < 0 .or. &
         size(seasonal_harmonics) /= size(seasonal_frequency)) then
         out%info = 1
         return
      end if
      if (.not. valid_distribution_parameters(distribution, &
         structural_shape(parameters(:structural_count), distribution))) then
         out%info = 1
         return
      end if
      allocate(use_observation(size(observations)))
      use_observation = .true.
      if (present(observed)) then
         if (size(observed) /= size(observations)) then
            out%info = 1
            return
         end if
         use_observation = observed
      end if
      if (count(use_observation) <= total_count) then
         out%info = 1
         return
      end if
      step = epsilon(1.0_dp)**0.25_dp
      if (present(difference_step)) step = difference_step
      if (step <= 0.0_dp) then
         out%info = 1
         return
      end if
      allocate(minus(size(observations), total_count), plus(size(observations), total_count))
      allocate(steps(total_count))
      trial = parameters
      do i = 1, total_count
         steps(i) = step*max(1.0_dp, abs(parameters(i)))
         trial(i) = parameters(i) - steps(i)
         minus(:, i) = contributions(trial)
         trial(i) = parameters(i) + steps(i)
         plus(:, i) = contributions(trial)
         trial(i) = parameters(i)
      end do
      if (.not. score_perturbations_valid(minus, plus)) then
         out%info = 4
         return
      end if
      all_scores = tsissm_fd_scores(minus, plus, steps)
      hessian = finite_difference_hessian(objective, parameters, step)
      out = tsissm_parameter_covariance(pack_rows(all_scores, use_observation), hessian, &
         covariance_type, adjust, newey_west_lag, hac_kernel, prewhite_order, hac_bandwidth)
      out%scores = pack_rows(all_scores, use_observation)
      out%hessian = hessian

   contains

      pure function contributions(trial_parameters) result(value)
         !! Return joint dynamic negative log-likelihood contributions.
         real(dp), intent(in) :: trial_parameters(:) !! Trial parameters.
         real(dp) :: value(size(observations))
         type(tsissm_filter_t) :: constant_filter, filtered
         type(tsissm_likelihood_t) :: likelihood
         type(tsissm_model_t) :: model
         type(tsissm_variance_initialization_t) :: initialized_variance
         real(dp), allocatable :: ar(:), arch(:), coefficients(:), garch(:), ma(:)
         real(dp), allocatable :: initial_arch(:), initial_variance(:), seasonal_persistence(:)
         real(dp) :: alpha, beta, damping, intercept, lambda, shape, skew, variance
         integer :: history, index

         index = 1
         alpha = trial_parameters(index)
         index = index + 1
         beta = 0.0_dp
         if (slope) then
            beta = trial_parameters(index)
            index = index + 1
         end if
         damping = 1.0_dp
         if (damped_slope) then
            damping = trial_parameters(index)
            index = index + 1
         end if
         seasonal_persistence = trial_parameters(index:index + size(seasonal_frequency) - 1)
         index = index + size(seasonal_frequency)
         ar = trial_parameters(index:index + ar_order - 1)
         index = index + ar_order
         ma = trial_parameters(index:index + ma_order - 1)
         index = index + ma_order
         coefficients = trial_parameters(index:index + size(regressors, 2) - 1)
         index = index + size(regressors, 2)
         lambda = trial_parameters(index)
         index = index + 1
         skew = 0.0_dp
         shape = 0.0_dp
         if (distribution == tsissm_distribution_student) then
            shape = trial_parameters(index)
            index = index + 1
         else if (distribution == tsissm_distribution_johnson_su) then
            skew = trial_parameters(index)
            shape = trial_parameters(index + 1)
            index = index + 2
         end if
         arch = trial_parameters(index:index + arch_order - 1)
         index = index + arch_order
         garch = trial_parameters(index:)
         model = tsissm_model(slope, damped_slope, seasonal_frequency, seasonal_harmonics, &
            regular_seasonal, ar, ma, alpha, beta, seasonal_persistence, damping)
         if (model%info /= 0) then
            value = huge(1.0_dp)
            return
         end if
         if (size(model%initial_state) /= size(initial_state)) then
            value = huge(1.0_dp)
            return
         end if
         constant_filter = tsissm_filter_constant(observations, model%transition_base, &
            model%transition_scale, model%transition_parameter, model%observation_loading, &
            model%persistence, regressors, coefficients, initial_state, use_observation, lambda)
         if (constant_filter%info /= 0) then
            value = huge(1.0_dp)
            return
         end if
         initialized_variance = tsissm_initialize_variance(constant_filter%innovation, &
            use_observation, variance_initialization, variance_sample_size)
         if (initialized_variance%info /= 0) then
            value = huge(1.0_dp)
            return
         end if
         variance = initialized_variance%variance
         intercept = max(1.0e-12_dp, variance*(1.0_dp - sum(arch) - sum(garch)))
         history = max(arch_order, garch_order)
         allocate(initial_arch(history), initial_variance(history))
         initial_arch = variance
         initial_variance = variance
         filtered = tsissm_filter_dynamic(observations, model%transition_base, &
            model%transition_scale, model%transition_parameter, model%observation_loading, &
            model%persistence, regressors, coefficients, initial_state, use_observation, lambda, &
            arch, garch, initial_arch, initial_variance, intercept)
         if (filtered%info /= 0) then
            value = huge(1.0_dp)
            return
         end if
         likelihood = tsissm_likelihood(filtered%innovation, observations, lambda, distribution, &
            skew, shape, total_count, use_observation, filtered%conditional_sd)
         if (likelihood%info /= 0) then
            value = huge(1.0_dp)
         else
            value = likelihood%contribution
         end if
      end function contributions

      pure real(dp) function objective(trial_parameters) result(value)
         !! Sum joint dynamic negative log-likelihood contributions.
         real(dp), intent(in) :: trial_parameters(:) !! Trial parameters.

         value = sum(contributions(trial_parameters))
      end function objective
   end function tsissm_covariance_structural_dynamic

   pure real(dp) function structural_shape(parameters, distribution) result(value)
      !! Extract the structural parameter vector's distribution shape.
      real(dp), intent(in) :: parameters(:) !! Model parameter values.
      integer, intent(in) :: distribution !! Probability-distribution specification.

      value = 0.0_dp
      if (distribution == tsissm_distribution_student) then
         value = parameters(size(parameters))
      else if (distribution == tsissm_distribution_johnson_su) then
         value = parameters(size(parameters))
      end if
   end function structural_shape

   pure logical function score_perturbations_valid(minus, plus) result(valid)
      !! Check finite-difference likelihood contributions before forming curvature.
      real(dp), intent(in) :: minus(:, :) !! Minus.
      real(dp), intent(in) :: plus(:, :) !! Plus.

      valid = all(ieee_is_finite(minus)) .and. all(ieee_is_finite(plus)) .and. &
         all(abs(minus) < sqrt(huge(1.0_dp))) .and. &
         all(abs(plus) < sqrt(huge(1.0_dp)))
   end function score_perturbations_valid

   pure logical function covariance_fit_inputs_valid(observations, filtered, regressors, &
      coefficients) result(valid)
      !! Validate data and stored filter dimensions used by covariance wrappers.
      real(dp), intent(in) :: observations(:) !! Observed time-series values.
      real(dp), intent(in) :: regressors(:, :) !! Regression design matrix.
      real(dp), intent(in) :: coefficients(:) !! Model coefficients.
      type(tsissm_filter_t), intent(in) :: filtered !! Filtered.

      valid = .false.
      if (.not. allocated(filtered%state) .or. .not. allocated(filtered%transition) .or. &
         .not. allocated(filtered%observation_loading) .or. .not. allocated(filtered%persistence)) return
      valid = size(observations) > 0 .and. all(observations > 0.0_dp) .and. &
         size(regressors, 1) == size(observations) .and. size(regressors, 2) == size(coefficients) .and. &
         lbound(filtered%state, 2) == 0 .and. ubound(filtered%state, 2) >= size(observations) .and. &
         size(filtered%transition, 1) == size(filtered%state, 1) .and. &
         size(filtered%transition, 2) == size(filtered%state, 1) .and. &
         size(filtered%observation_loading) == size(filtered%state, 1)
   end function covariance_fit_inputs_valid

   pure function tsissm_select_models(observations, regressors, slope_options, damped_options, &
      include_nonseasonal, seasonal_frequency, seasonal_harmonic_grid, regular_seasonal, &
      ar_orders, ma_orders, include_constant, include_dynamic, garch_order, distribution, &
      initial_lambda, initial_skew, initial_shape, top_n, criterion, max_iterations, &
      tolerance, variance_initialization, variance_sample_size) result(out)
      !! Enumerate, fit, screen, and rank a configurable tsissm structural grid.
      real(dp), intent(in) :: observations(:) !! Observed time-series values.
      real(dp), intent(in) :: regressors(:, :) !! Regression design matrix.
      real(dp), intent(in) :: seasonal_frequency(:) !! Seasonal frequency.
      logical, intent(in) :: slope_options(:) !! Flag controlling slope options.
      logical, intent(in) :: damped_options(:) !! Flag controlling damped options.
      logical, intent(in) :: include_nonseasonal !! Whether to include the nonseasonal.
      integer, intent(in) :: seasonal_harmonic_grid(:, :) !! Seasonal harmonic grid.
      integer, intent(in) :: ar_orders(:) !! Autoregressive orders.
      integer, intent(in) :: ma_orders(:) !! Moving-average orders.
      logical, intent(in) :: regular_seasonal !! Flag controlling regular seasonal.
      logical, intent(in) :: include_constant !! Whether to include the constant.
      logical, intent(in) :: include_dynamic !! Whether to include the dynamic.
      integer, intent(in) :: garch_order(2) !! Garch order.
      integer, intent(in) :: distribution !! Probability-distribution specification.
      integer, intent(in) :: top_n !! Top n.
      integer, intent(in) :: criterion !! Criterion.
      integer, intent(in) :: max_iterations !! Maximum number of algorithm iterations.
      real(dp), intent(in) :: initial_lambda !! Initial lambda.
      real(dp), intent(in) :: initial_skew !! Initial skew.
      real(dp), intent(in) :: initial_shape !! Initial shape.
      real(dp), intent(in) :: tolerance !! Numerical convergence tolerance.
      integer, intent(in), optional :: variance_initialization !! Variance initialization.
      integer, intent(in), optional :: variance_sample_size !! Variance sample size.
      type(tsissm_selection_t) :: out
      type(tsissm_model_t) :: initial_model
      type(tsissm_structural_fit_t) :: dynamic_structural_fit, structural_fit
      real(dp), allocatable :: initial_parameters(:), initial_state(:)
      integer, allocatable :: harmonics(:), no_harmonics(:)
      real(dp), allocatable :: no_frequency(:)
      integer :: ar_index, candidate_index, configuration, damped_index, i, ma_index
      integer :: seasonal_cases, slope_index, trend_cases, variance_cases
      logical :: damped, seasonal, slope

      if (size(observations) < 3 .or. any(observations <= 0.0_dp) .or. &
         size(regressors, 1) /= size(observations) .or. size(slope_options) < 1 .or. &
         size(damped_options) < 1 .or. size(ar_orders) < 1 .or. size(ma_orders) < 1 .or. &
         any(ar_orders < 0) .or. any(ma_orders < 0) .or. &
         size(seasonal_harmonic_grid, 2) /= size(seasonal_frequency) .or. &
         (.not. include_nonseasonal .and. size(seasonal_harmonic_grid, 1) < 1) .or. &
         (.not. include_constant .and. .not. include_dynamic) .or. &
         (include_dynamic .and. maxval(garch_order) < 1) .or. any(garch_order < 0) .or. &
         top_n < 1 .or. criterion < tsissm_selection_aic .or. &
         criterion > tsissm_selection_bic .or. max_iterations < 1 .or. tolerance <= 0.0_dp .or. &
         .not. valid_distribution_parameters(distribution, initial_shape)) then
         out%info = 1
         return
      end if
      if (size(seasonal_harmonic_grid, 1) > 0) then
         if (any(seasonal_frequency <= 1.0_dp) .or. any(seasonal_harmonic_grid < 1) .or. &
            (.not. regular_seasonal .and. &
            any(2*seasonal_harmonic_grid > spread(nint(seasonal_frequency), 1, &
            size(seasonal_harmonic_grid, 1))))) then
            out%info = 1
            return
         end if
      end if
      trend_cases = 0
      do slope_index = 1, size(slope_options)
         do damped_index = 1, size(damped_options)
            if (.not. slope_options(slope_index) .and. damped_options(damped_index)) cycle
            trend_cases = trend_cases + 1
         end do
      end do
      seasonal_cases = size(seasonal_harmonic_grid, 1) + merge(1, 0, include_nonseasonal)
      variance_cases = merge(1, 0, include_constant) + merge(1, 0, include_dynamic)
      allocate(out%candidates(trend_cases*seasonal_cases*size(ar_orders)*size(ma_orders)* &
         variance_cases))
      allocate(no_frequency(0), no_harmonics(0))
      candidate_index = 0
      do slope_index = 1, size(slope_options)
         slope = slope_options(slope_index)
         do damped_index = 1, size(damped_options)
            damped = damped_options(damped_index)
            if (.not. slope .and. damped) cycle
            do configuration = 0, seasonal_cases - 1
               seasonal = configuration >= merge(1, 0, include_nonseasonal)
               if (seasonal) then
                  harmonics = seasonal_harmonic_grid(configuration + &
                     merge(0, 1, include_nonseasonal), :)
               else
                  harmonics = no_harmonics
               end if
               do ar_index = 1, size(ar_orders)
                  do ma_index = 1, size(ma_orders)
                     if (seasonal) then
                        initial_model = tsissm_model(slope, damped, seasonal_frequency, harmonics, &
                           regular_seasonal, [(0.0_dp, i=1, ar_orders(ar_index))], &
                           [(0.0_dp, i=1, ma_orders(ma_index))], 0.2_dp, 0.1_dp, &
                           [(0.1_dp, i=1, size(seasonal_frequency))], 0.95_dp)
                     else
                        initial_model = tsissm_model(slope, damped, no_frequency, no_harmonics, &
                           regular_seasonal, [(0.0_dp, i=1, ar_orders(ar_index))], &
                           [(0.0_dp, i=1, ma_orders(ma_index))], 0.2_dp, 0.1_dp, &
                           [real(dp) ::], 0.95_dp)
                     end if
                     if (initial_model%info /= 0) then
                        out%info = 2
                        return
                     end if
                     allocate(initial_state(size(initial_model%initial_state)))
                     initial_state = 0.0_dp
                     initial_parameters = selection_initial_parameters(slope, damped, &
                        merge(size(seasonal_frequency), 0, seasonal), ar_orders(ar_index), &
                        ma_orders(ma_index), size(regressors, 2), initial_lambda, distribution, &
                        initial_skew, initial_shape)
                     if (seasonal) then
                        structural_fit = tsissm_fit_structural_constant(observations, slope, &
                           damped, seasonal_frequency, harmonics, regular_seasonal, &
                           ar_orders(ar_index), ma_orders(ma_index), regressors, initial_state, &
                           initial_parameters, distribution=distribution, &
                           max_iterations=max_iterations, tolerance=tolerance)
                     else
                        structural_fit = tsissm_fit_structural_constant(observations, slope, &
                           damped, no_frequency, no_harmonics, regular_seasonal, &
                           ar_orders(ar_index), ma_orders(ma_index), regressors, initial_state, &
                           initial_parameters, distribution=distribution, &
                           max_iterations=max_iterations, tolerance=tolerance)
                     end if
                     if (include_constant) then
                        candidate_index = candidate_index + 1
                        call store_selection_candidate(out%candidates(candidate_index), &
                           structural_fit, slope, damped, seasonal, harmonics, ar_orders(ar_index), &
                           ma_orders(ma_index), .false., distribution, garch_order, regressors, &
                           observations, max_iterations, tolerance)
                     end if
                     if (include_dynamic) then
                        if (seasonal) then
                           dynamic_structural_fit = tsissm_fit_structural_dynamic(observations, &
                              slope, damped, seasonal_frequency, harmonics, regular_seasonal, &
                              ar_orders(ar_index), ma_orders(ma_index), regressors, initial_state, &
                              initial_parameters, selection_variance_initial(garch_order(1), 0.05_dp), &
                              selection_variance_initial(garch_order(2), 0.8_dp), &
                              distribution=distribution, max_iterations=max_iterations, &
                              tolerance=tolerance, variance_initialization=variance_initialization, &
                              variance_sample_size=variance_sample_size)
                        else
                           dynamic_structural_fit = tsissm_fit_structural_dynamic(observations, &
                              slope, damped, no_frequency, no_harmonics, regular_seasonal, &
                              ar_orders(ar_index), ma_orders(ma_index), regressors, initial_state, &
                              initial_parameters, selection_variance_initial(garch_order(1), 0.05_dp), &
                              selection_variance_initial(garch_order(2), 0.8_dp), &
                              distribution=distribution, max_iterations=max_iterations, &
                              tolerance=tolerance, variance_initialization=variance_initialization, &
                              variance_sample_size=variance_sample_size)
                        end if
                        candidate_index = candidate_index + 1
                        call store_selection_candidate(out%candidates(candidate_index), &
                           dynamic_structural_fit, slope, damped, seasonal, harmonics, &
                           ar_orders(ar_index), &
                           ma_orders(ma_index), .true., distribution, garch_order, regressors, &
                           observations, max_iterations, tolerance)
                     end if
                     deallocate(initial_state)
                  end do
               end do
            end do
         end do
      end do
      out%evaluated = candidate_index
      out%successful = count(out%candidates%successful)
      out%criterion = criterion
      if (out%successful < 1) then
         out%info = 3
         return
      end if
      call rank_selection(out, top_n, criterion)
   end function tsissm_select_models

   pure function tsissm_rerank_models(selection, scores, top_n) result(out)
      !! Re-rank fitted candidates using supplied backtest or validation losses.
      type(tsissm_selection_t), intent(in) :: selection !! Selection.
      real(dp), intent(in) :: scores(:) !! Scores.
      integer, intent(in) :: top_n !! Top n.
      type(tsissm_selection_t) :: out
      integer :: i

      if (size(scores) /= size(selection%candidates) .or. top_n < 1) then
         out%info = 1
         return
      end if
      out = selection
      if (allocated(out%top_index)) deallocate(out%top_index)
      if (allocated(out%weight)) deallocate(out%weight)
      do i = 1, size(out%candidates)
         out%candidates(i)%weight = 0.0_dp
         if (.not. out%candidates(i)%successful) cycle
         if (.not. ieee_is_finite(scores(i))) then
            out%candidates(i)%successful = .false.
         else
            out%candidates(i)%score = scores(i)
         end if
      end do
      out%successful = count(out%candidates%successful)
      out%criterion = tsissm_selection_external
      if (out%successful < 1) then
         out%info = 2
         return
      end if
      out%info = 0
      call rank_selection(out, top_n, tsissm_selection_external)
   end function tsissm_rerank_models

   pure function tsissm_selection_correlation(selection) result(correlation)
      !! Estimate the selected models' Gaussian-copula correlation from Kendall tau.
      type(tsissm_selection_t), intent(in) :: selection !! Selection.
      real(dp), allocatable :: correlation(:, :)
      real(dp), allocatable :: residual(:, :)
      real(dp) :: concordance, denominator, difference_first, difference_second, tau
      integer :: first, i, j, model_count, n, second

      if (.not. allocated(selection%top_index)) then
         allocate(correlation(0, 0))
         return
      end if
      model_count = size(selection%top_index)
      if (model_count < 1) then
         allocate(correlation(0, 0))
         return
      end if
      n = size(selection%candidates(selection%top_index(1))%filter%innovation)
      if (n < 2) then
         allocate(correlation(0, 0))
         return
      end if
      allocate(correlation(model_count, model_count), residual(n, model_count))
      correlation = 0.0_dp
      do i = 1, model_count
         if (allocated(selection%candidates(selection%top_index(i))%filter%conditional_sd)) then
            residual(:, i) = selection%candidates(selection%top_index(i))%filter%innovation/ &
               selection%candidates(selection%top_index(i))%filter%conditional_sd
         else
            residual(:, i) = selection%candidates(selection%top_index(i))%filter%innovation/ &
               selection%candidates(selection%top_index(i))%likelihood%sigma
         end if
         correlation(i, i) = 1.0_dp
      end do
      denominator = 0.5_dp*real(n, dp)*real(n - 1, dp)
      do first = 1, model_count - 1
         do second = first + 1, model_count
            concordance = 0.0_dp
            do i = 1, n - 1
               do j = i + 1, n
                  difference_first = residual(j, first) - residual(i, first)
                  difference_second = residual(j, second) - residual(i, second)
                  concordance = concordance + sign(1.0_dp, difference_first*difference_second)
                  if (difference_first == 0.0_dp .or. difference_second == 0.0_dp) &
                     concordance = concordance - sign(1.0_dp, difference_first*difference_second)
               end do
            end do
            tau = concordance/denominator
            correlation(first, second) = sin(0.5_dp*acos(-1.0_dp)*tau)
            correlation(second, first) = correlation(first, second)
         end do
      end do
   end function tsissm_selection_correlation

   function tsissm_ensemble_forecast(selection, regressors, simulations, probabilities, seed, &
      weights, correlation_matrix) result(out)
      !! Forecast selected models with Gaussian-copula innovations and combine paths.
      type(tsissm_selection_t), intent(in) :: selection !! Selection.
      real(dp), intent(in) :: regressors(:, :) !! Regression design matrix.
      integer, intent(in) :: simulations !! Number of simulation draws.
      real(dp), intent(in), optional :: probabilities(:) !! Probability values.
      real(dp), intent(in), optional :: weights(:) !! Observation or objective weights.
      real(dp), intent(in), optional :: correlation_matrix(:, :) !! Correlation matrix.
      integer, intent(in), optional :: seed !! Random-number seed.
      type(tsissm_ensemble_t) :: out
      type(tsissm_prediction_decomposition_t) :: decomposition
      real(dp), allocatable :: copula_normal(:), eigenvalues(:), eigenvectors(:, :)
      real(dp), allocatable :: innovations(:, :, :), selected_probabilities(:), standard(:)
      real(dp) :: lambda_maximum, lambda_minimum, total_weight
      integer :: h, i, info, model, model_count, simulation

      if (.not. allocated(selection%top_index) .or. simulations < 2 .or. &
         size(regressors, 1) < 1) then
         out%info = 1
         return
      end if
      model_count = size(selection%top_index)
      if (model_count < 1) then
         out%info = 1
         return
      end if
      if (size(regressors, 2) /= &
         size(selection%candidates(selection%top_index(1))%coefficients)) then
         out%info = 1
         return
      end if
      call prepare_probabilities(probabilities, selected_probabilities, info)
      if (info /= 0) then
         out%info = 2
         return
      end if
      if (present(weights)) then
         if (size(weights) /= model_count .or. any(weights < 0.0_dp) .or. sum(weights) <= 0.0_dp) then
            out%info = 3
            return
         end if
         out%weight = weights/sum(weights)
      else
         if (size(selection%weight) /= model_count) then
            out%info = 3
            return
         end if
         out%weight = selection%weight
      end if
      if (present(correlation_matrix)) then
         if (any(shape(correlation_matrix) /= [model_count, model_count])) then
            out%info = 4
            return
         end if
         out%correlation = correlation_matrix
      else
         out%correlation = tsissm_selection_correlation(selection)
      end if
      call copula_eigensystem(out%correlation, eigenvalues, eigenvectors, info)
      if (info /= 0) then
         out%info = 5
         return
      end if
      if (present(seed)) call set_random_seed(seed)
      h = size(regressors, 1)
      allocate(innovations(simulations, h, model_count), standard(model_count))
      allocate(copula_normal(model_count))
      do i = 1, h
         do simulation = 1, simulations
            do model = 1, model_count
               standard(model) = random_standard_normal()
            end do
            copula_normal = matmul(eigenvectors, sqrt(max(eigenvalues, 0.0_dp))*standard)
            do model = 1, model_count
               innovations(simulation, i, model) = standardized_innovation_from_normal( &
                  copula_normal(model), &
                  selection%candidates(selection%top_index(model))%distribution, &
                  selection%candidates(selection%top_index(model))%skew, &
                  selection%candidates(selection%top_index(model))%shape)
            end do
         end do
      end do
      allocate(out%model_forecast(model_count))
      allocate(out%distribution(simulations, h), out%analytic_mean(h))
      out%distribution = 0.0_dp
      out%analytic_mean = 0.0_dp
      do model = 1, model_count
         out%model_forecast(model) = forecast_selection_candidate( &
            selection%candidates(selection%top_index(model)), regressors, &
            innovations(:, :, model), selected_probabilities)
         if (out%model_forecast(model)%info /= 0) then
            out%info = 10 + model
            return
         end if
         out%distribution = out%distribution + out%weight(model)* &
            out%model_forecast(model)%distribution
         out%analytic_mean = out%analytic_mean + out%weight(model)* &
            out%model_forecast(model)%moments%mean
      end do
      call summarize_ensemble(out, selected_probabilities)
      lambda_minimum = selection%candidates(selection%top_index(1))%lambda
      lambda_maximum = lambda_minimum
      do model = 2, model_count
         lambda_minimum = min(lambda_minimum, selection%candidates(selection%top_index(model))%lambda)
         lambda_maximum = max(lambda_maximum, selection%candidates(selection%top_index(model))%lambda)
      end do
      if (lambda_maximum - lambda_minimum <= sqrt(epsilon(1.0_dp))) then
         allocate(out%trend(simulations, h), out%seasonal(simulations, h), out%arma(simulations, h))
         allocate(out%regression(simulations, h), out%irregular(simulations, h))
         out%trend = 0.0_dp
         out%seasonal = 0.0_dp
         out%arma = 0.0_dp
         out%regression = 0.0_dp
         out%irregular = 0.0_dp
         do model = 1, model_count
            decomposition = tsissm_decompose_prediction( &
               selection%candidates(selection%top_index(model))%model, &
               out%model_forecast(model)%prediction, regressors, &
               selection%candidates(selection%top_index(model))%coefficients, lambda_minimum)
            if (decomposition%info /= 0) then
               out%info = 20 + model
               return
            end if
            out%trend = out%trend + out%weight(model)*decomposition%trend
            out%seasonal = out%seasonal + out%weight(model)*decomposition%seasonal_total
            out%arma = out%arma + out%weight(model)*decomposition%arma
            out%regression = out%regression + out%weight(model)*decomposition%regression
            out%irregular = out%irregular + out%weight(model)*decomposition%irregular
         end do
         out%has_decomposition = .true.
      end if
      total_weight = sum(out%weight)
      if (abs(total_weight - 1.0_dp) > 10.0_dp*epsilon(1.0_dp)) out%info = 6
   end function tsissm_ensemble_forecast

   pure function forecast_selection_candidate(candidate, regressors, standardized, &
      probabilities) result(out)
      !! Forecast one selected candidate from supplied standardized innovations.
      type(tsissm_candidate_t), intent(in) :: candidate !! Candidate.
      real(dp), intent(in) :: regressors(:, :) !! Regression design matrix.
      real(dp), intent(in) :: standardized(:, :) !! Standardized.
      real(dp), intent(in) :: probabilities(:) !! Probability values.
      type(tsissm_forecast_t) :: out
      real(dp), allocatable :: arch_history(:), innovations(:, :), variance_history(:)
      integer :: history, last_state

      last_state = ubound(candidate%filter%state, 2)
      if (candidate%dynamic_variance) then
         history = max(size(candidate%arch), size(candidate%garch))
         if (history < 1 .or. size(candidate%filter%innovation) < history .or. &
            .not. allocated(candidate%filter%conditional_sd)) then
            out%info = 1
            return
         end if
         arch_history = candidate%filter%innovation( &
            size(candidate%filter%innovation) - history + 1:)**2
         variance_history = candidate%filter%conditional_sd( &
            size(candidate%filter%conditional_sd) - history + 1:)**2
         out%prediction = tsissm_predict_dynamic(candidate%model%transition_base, &
            candidate%model%transition_scale, candidate%model%transition_parameter, &
            candidate%model%observation_loading, candidate%persistence, regressors, &
            candidate%coefficients, candidate%filter%state(:, last_state), standardized, &
            candidate%arch, candidate%garch, arch_history, variance_history, &
            candidate%variance_intercept)
         out%moments = tsissm_moments_dynamic(candidate%model%transition_base* &
            candidate%model%transition_scale*candidate%model%transition_parameter, &
            candidate%model%observation_loading, candidate%persistence, &
            candidate%filter%state(:, last_state), regressors, candidate%coefficients, &
            candidate%arch, candidate%garch, candidate%variance_intercept, arch_history, &
            variance_history, candidate%lambda, .true.)
      else
         innovations = candidate%likelihood%sigma*standardized
         out%prediction = tsissm_predict_constant(candidate%model%transition_base, &
            candidate%model%transition_scale, candidate%model%transition_parameter, &
            candidate%model%observation_loading, candidate%persistence, regressors, &
            candidate%coefficients, candidate%filter%state(:, last_state), innovations)
         out%moments = tsissm_moments_constant(candidate%model%transition_base* &
            candidate%model%transition_scale*candidate%model%transition_parameter, &
            candidate%model%observation_loading, candidate%persistence, &
            candidate%filter%state(:, last_state), regressors, candidate%coefficients, &
            candidate%likelihood%sigma**2, candidate%lambda, .true.)
      end if
      if (out%prediction%info /= 0 .or. out%moments%info /= 0) then
         out%info = 2
         return
      end if
      call summarize_forecast(out, probabilities, candidate%lambda)
   end function forecast_selection_candidate

   pure subroutine summarize_ensemble(out, probabilities)
      !! Calculate ensemble path means and type-7 quantiles.
      type(tsissm_ensemble_t), intent(inout) :: out !! Procedure result, updated in place.
      real(dp), intent(in) :: probabilities(:) !! Probability values.
      real(dp), allocatable :: ordered(:)
      integer :: horizon, probability_index

      out%probabilities = probabilities
      allocate(out%mean(size(out%distribution, 2)))
      allocate(out%quantile(size(probabilities), size(out%distribution, 2)))
      out%mean = sum(out%distribution, dim=1)/real(size(out%distribution, 1), dp)
      allocate(ordered(size(out%distribution, 1)))
      do horizon = 1, size(out%distribution, 2)
         ordered = out%distribution(:, horizon)
         call sort(ordered)
         do probability_index = 1, size(probabilities)
            out%quantile(probability_index, horizon) = &
               quantile(ordered, probabilities(probability_index))
         end do
      end do
   end subroutine summarize_ensemble

   pure subroutine copula_eigensystem(correlation, eigenvalues, eigenvectors, info)
      !! Project a correlation matrix to positive semidefinite form and factor it.
      real(dp), intent(inout) :: correlation(:, :) !! Correlation, updated in place.
      real(dp), allocatable, intent(out) :: eigenvalues(:) !! Eigenvalues.
      real(dp), allocatable, intent(out) :: eigenvectors(:, :) !! Eigenvectors.
      integer, intent(out) :: info !! Status code; zero indicates success.
      real(dp), allocatable :: projected(:, :), scale(:), values(:), vectors(:, :)
      integer :: i, j, n

      n = size(correlation, 1)
      if (n < 1 .or. size(correlation, 2) /= n .or. any(.not. ieee_is_finite(correlation))) then
         info = 1
         return
      end if
      correlation = 0.5_dp*(correlation + transpose(correlation))
      do i = 1, n
         if (correlation(i, i) <= 0.0_dp) then
            info = 1
            return
         end if
      end do
      scale = sqrt([(correlation(i, i), i=1,n)])
      do i = 1, n
         do j = 1, n
            correlation(i, j) = correlation(i, j)/(scale(i)*scale(j))
         end do
      end do
      call symmetric_eigen(correlation, values, vectors, info)
      if (info /= 0) return
      values = max(values, 100.0_dp*epsilon(1.0_dp))
      projected = matmul(vectors*spread(values, 1, n), transpose(vectors))
      scale = sqrt([(projected(i, i), i=1,n)])
      do i = 1, n
         do j = 1, n
            correlation(i, j) = projected(i, j)/(scale(i)*scale(j))
         end do
         correlation(i, i) = 1.0_dp
      end do
      call symmetric_eigen(correlation, eigenvalues, eigenvectors, info)
      if (info == 0) eigenvalues = max(eigenvalues, 0.0_dp)
   end subroutine copula_eigensystem

   pure real(dp) function standardized_innovation_from_normal(normal, distribution, skew, &
      shape) result(innovation)
      !! Map a Gaussian-copula variate to one standardized innovation marginal.
      real(dp), intent(in) :: normal !! Normal.
      real(dp), intent(in) :: skew !! Skew.
      real(dp), intent(in) :: shape !! Shape.
      integer, intent(in) :: distribution !! Probability-distribution specification.
      real(dp) :: c, omega, probability, reciprocal_shape, w

      select case (distribution)
      case (tsissm_distribution_gaussian)
         innovation = normal
      case (tsissm_distribution_student)
         probability = max(1.0e-12_dp, min(1.0_dp - 1.0e-12_dp, &
            0.5_dp*erfc(-normal/sqrt(2.0_dp))))
         innovation = student_quantile(probability, shape)*sqrt((shape - 2.0_dp)/shape)
      case (tsissm_distribution_johnson_su)
         reciprocal_shape = 1.0_dp/shape
         w = exp(reciprocal_shape**2)
         omega = -skew*reciprocal_shape
         c = sqrt(1.0_dp/(0.5_dp*(w - 1.0_dp)*(w*cosh(2.0_dp*omega) + 1.0_dp)))
         innovation = c*sinh(reciprocal_shape*(normal + skew)) + c*sqrt(w)*sinh(omega)
      case default
         innovation = 0.0_dp
      end select
   end function standardized_innovation_from_normal

   pure real(dp) function student_quantile(probability, degrees) result(value)
      !! Invert the Student-t CDF by symmetric monotone bisection.
      real(dp), intent(in) :: probability !! Probability value.
      real(dp), intent(in) :: degrees !! Degrees.
      real(dp) :: lower, middle, target, upper
      integer :: iteration

      if (probability == 0.5_dp) then
         value = 0.0_dp
         return
      end if
      target = max(probability, 1.0_dp - probability)
      lower = 0.0_dp
      upper = 1.0_dp
      do while (student_cdf(upper, degrees) < target .and. upper < 1.0e6_dp)
         upper = 2.0_dp*upper
      end do
      do iteration = 1, 100
         middle = 0.5_dp*(lower + upper)
         if (student_cdf(middle, degrees) < target) then
            lower = middle
         else
            upper = middle
         end if
      end do
      value = 0.5_dp*(lower + upper)
      if (probability < 0.5_dp) value = -value
   end function student_quantile

   pure real(dp) function student_cdf(value, degrees) result(probability)
      !! Evaluate a positive Student-t CDF through the regularized beta function.
      real(dp), intent(in) :: value !! Input value.
      real(dp), intent(in) :: degrees !! Degrees.
      real(dp) :: beta_argument

      beta_argument = degrees/(degrees + value**2)
      probability = 1.0_dp - 0.5_dp*regularized_beta_mts(beta_argument, &
         0.5_dp*degrees, 0.5_dp)
   end function student_cdf

   pure function selection_initial_parameters(slope, damped, seasonal_count, ar_order, ma_order, &
      regression_count, lambda, distribution, skew, shape) result(parameters)
      !! Build the structural estimator's flat initial-parameter vector.
      logical, intent(in) :: slope !! Flag controlling slope.
      logical, intent(in) :: damped !! Flag controlling damped.
      integer, intent(in) :: seasonal_count !! Number of seasonal.
      integer, intent(in) :: ar_order !! Autoregressive order.
      integer, intent(in) :: ma_order !! Moving-average order.
      integer, intent(in) :: regression_count !! Number of regression.
      integer, intent(in) :: distribution !! Probability-distribution specification.
      real(dp), intent(in) :: lambda !! Penalty or shrinkage parameter.
      real(dp), intent(in) :: skew !! Skew.
      real(dp), intent(in) :: shape !! Shape.
      real(dp), allocatable :: parameters(:)
      integer :: count, i, index

      count = 1 + merge(1, 0, slope) + merge(1, 0, damped) + seasonal_count + &
         ar_order + ma_order + regression_count + 1
      if (distribution == tsissm_distribution_student) count = count + 1
      if (distribution == tsissm_distribution_johnson_su) count = count + 2
      allocate(parameters(count))
      parameters = 0.0_dp
      index = 1
      parameters(index) = 0.2_dp
      index = index + 1
      if (slope) then
         parameters(index) = 0.1_dp
         index = index + 1
      end if
      if (damped) then
         parameters(index) = 0.95_dp
         index = index + 1
      end if
      do i = 1, seasonal_count
         parameters(index) = 0.1_dp
         index = index + 1
      end do
      index = index + ar_order + ma_order + regression_count
      parameters(index) = lambda
      index = index + 1
      if (distribution == tsissm_distribution_student) then
         parameters(index) = shape
      else if (distribution == tsissm_distribution_johnson_su) then
         parameters(index) = skew
         parameters(index + 1) = shape
      end if
   end function selection_initial_parameters

   pure function selection_variance_initial(order, total) result(parameters)
      !! Split an initial ARCH or GARCH total equally across its lags.
      integer, intent(in) :: order !! Model or polynomial order.
      real(dp), intent(in) :: total !! Total.
      real(dp), allocatable :: parameters(:)

      allocate(parameters(order))
      if (order > 0) parameters = total/real(order, dp)
   end function selection_variance_initial

   pure subroutine store_selection_candidate(candidate, structural, slope, damped, seasonal, &
      harmonics, ar_order, ma_order, dynamic, distribution, garch_order, regressors, &
      observations, max_iterations, tolerance)
      !! Convert a structural fit into a screened constant or dynamic candidate.
      type(tsissm_candidate_t), intent(out) :: candidate !! Candidate.
      type(tsissm_structural_fit_t), intent(in) :: structural !! Structural.
      logical, intent(in) :: slope !! Flag controlling slope.
      logical, intent(in) :: damped !! Flag controlling damped.
      logical, intent(in) :: seasonal !! Flag controlling seasonal.
      logical, intent(in) :: dynamic !! Flag controlling dynamic.
      integer, intent(in) :: harmonics(:) !! Harmonics.
      integer, intent(in) :: ar_order !! Autoregressive order.
      integer, intent(in) :: ma_order !! Moving-average order.
      integer, intent(in) :: distribution !! Probability-distribution specification.
      integer, intent(in) :: garch_order(2) !! Garch order.
      real(dp), intent(in) :: regressors(:, :) !! Regression design matrix.
      real(dp), intent(in) :: observations(:) !! Observed time-series values.
      real(dp), intent(in) :: tolerance !! Numerical convergence tolerance.
      integer, intent(in) :: max_iterations !! Maximum number of algorithm iterations.
      integer :: coefficient_start, parameter_count

      candidate%slope = slope
      candidate%damped_slope = damped
      candidate%seasonal = seasonal
      candidate%dynamic_variance = dynamic
      candidate%ar_order = ar_order
      candidate%ma_order = ma_order
      candidate%distribution = distribution
      candidate%seasonal_harmonics = harmonics
      if (structural%info < 0 .or. structural%info > 4 .or. &
         .not. allocated(structural%parameters)) then
         candidate%info = 10 + max(0, structural%info)
         return
      end if
      candidate%model = structural%model
      candidate%filter = structural%filter
      candidate%likelihood = structural%likelihood
      candidate%parameters = structural%parameters
      candidate%persistence = structural%model%persistence
      candidate%lambda = structural%lambda
      candidate%skew = structural%skew
      candidate%shape = structural%shape
      coefficient_start = 2 + merge(1, 0, slope) + merge(1, 0, damped) + &
         merge(size(harmonics), 0, seasonal) + ar_order + ma_order
      candidate%coefficients = structural%parameters( &
         coefficient_start:coefficient_start + size(regressors, 2) - 1)
      if (dynamic) then
         if (.not. structural%dynamic_variance .or. .not. allocated(structural%arch) .or. &
            .not. allocated(structural%garch)) then
            candidate%info = 20
            return
         end if
         candidate%arch = structural%arch
         candidate%garch = structural%garch
         candidate%variance_intercept = structural%variance_intercept
      else
         allocate(candidate%arch(0), candidate%garch(0))
      end if
      parameter_count = size(structural%parameters)
      candidate%parameter_count = parameter_count
      candidate%aic = -2.0_dp*candidate%likelihood%log_likelihood + &
         2.0_dp*real(parameter_count, dp)
      candidate%aicc = candidate%aic
      if (size(observations) > parameter_count + 1) then
         candidate%aicc = candidate%aic + 2.0_dp*real(parameter_count*(parameter_count + 1), dp)/ &
            real(size(observations) - parameter_count - 1, dp)
      else
         candidate%aicc = huge(1.0_dp)
      end if
      candidate%bic = -2.0_dp*candidate%likelihood%log_likelihood + &
         log(real(size(observations), dp))*real(parameter_count, dp)
      candidate%successful = selection_candidate_feasible(candidate, structural%parameters)
      if (.not. candidate%successful) candidate%info = 30
   end subroutine store_selection_candidate

   pure logical function selection_candidate_feasible(candidate, structural_parameters) result(feasible)
      !! Check likelihood finiteness and ISSM, AR, MA, and GARCH stability.
      type(tsissm_candidate_t), intent(in) :: candidate !! Candidate.
      real(dp), intent(in) :: structural_parameters(:) !! Structural parameters.
      integer :: ar_start, seasonal_count

      feasible = candidate%likelihood%info == 0 .and. &
         ieee_is_finite(candidate%likelihood%log_likelihood)
      if (.not. feasible) return
      feasible = all(candidate%persistence >= 0.0_dp) .and. &
         all(candidate%persistence <= 1.0_dp) .and. candidate%lambda >= 0.0_dp .and. &
         candidate%lambda <= 1.5_dp .and. &
         valid_distribution_parameters(candidate%distribution, candidate%shape)
      if (.not. feasible) return
      feasible = tsissm_stability_constraint(candidate%model%transition_base, &
         candidate%model%transition_scale, candidate%model%transition_parameter, &
         candidate%model%observation_loading, candidate%persistence, 0.001_dp) >= 0.0_dp
      if (.not. feasible) return
      seasonal_count = merge(size(candidate%seasonal_harmonics), 0, candidate%seasonal)
      ar_start = 2 + merge(1, 0, candidate%slope) + merge(1, 0, candidate%damped_slope) + &
         seasonal_count
      if (candidate%ar_order > 0) feasible = tsissm_ar_constraint( &
         structural_parameters(ar_start:ar_start + candidate%ar_order - 1), 0.001_dp) >= 0.0_dp
      if (.not. feasible) return
      if (candidate%ma_order > 0) feasible = tsissm_ma_constraint( &
         structural_parameters(ar_start + candidate%ar_order: &
         ar_start + candidate%ar_order + candidate%ma_order - 1), 0.001_dp) >= 0.0_dp
      if (candidate%dynamic_variance) feasible = feasible .and. &
         all(candidate%arch >= 0.0_dp) .and. all(candidate%garch >= 0.0_dp) .and. &
         sum(candidate%arch) + sum(candidate%garch) < 1.0_dp
   end function selection_candidate_feasible

   pure subroutine rank_selection(out, top_n, criterion)
      !! Sort successful candidates and normalize information-criterion weights.
      type(tsissm_selection_t), intent(inout) :: out !! Procedure result, updated in place.
      integer, intent(in) :: top_n !! Top n.
      integer, intent(in) :: criterion !! Criterion.
      integer, allocatable :: order(:)
      integer :: held, i, j, retained
      real(dp) :: minimum_score, total

      allocate(order(out%successful))
      j = 0
      do i = 1, size(out%candidates)
         if (.not. out%candidates(i)%successful) cycle
         j = j + 1
         order(j) = i
         select case (criterion)
         case (tsissm_selection_aic)
            out%candidates(i)%score = out%candidates(i)%aic
         case (tsissm_selection_aicc)
            out%candidates(i)%score = out%candidates(i)%aicc
         case (tsissm_selection_bic)
            out%candidates(i)%score = out%candidates(i)%bic
         case (tsissm_selection_external)
            continue
         end select
      end do
      do i = 2, size(order)
         held = order(i)
         j = i - 1
         do while (j >= 1)
            if (out%candidates(order(j))%score <= out%candidates(held)%score) exit
            order(j + 1) = order(j)
            j = j - 1
         end do
         order(j + 1) = held
      end do
      retained = min(top_n, out%successful)
      allocate(out%top_index(retained), out%weight(retained))
      out%top_index = order(:retained)
      minimum_score = out%candidates(out%top_index(1))%score
      do i = 1, retained
         out%weight(i) = exp(-0.5_dp*(out%candidates(out%top_index(i))%score - minimum_score))
      end do
      total = sum(out%weight)
      out%weight = out%weight/total
      do i = 1, retained
         out%candidates(out%top_index(i))%weight = out%weight(i)
      end do
   end subroutine rank_selection

   pure function tsissm_diagnose(filtered, original_observations, lambda, lag, &
      fitted_parameter_count, constant_sigma) result(out)
      !! Diagnose transformed and standardized one-step innovation residuals.
      type(tsissm_filter_t), intent(in) :: filtered !! Filtered.
      real(dp), intent(in) :: original_observations(:) !! Original observations.
      real(dp), intent(in) :: lambda !! Penalty or shrinkage parameter.
      integer, intent(in) :: lag !! Lag index or number of lags.
      integer, intent(in) :: fitted_parameter_count !! Number of fitted parameter.
      real(dp), intent(in), optional :: constant_sigma !! Constant sigma.
      type(tsissm_diagnostics_t) :: out
      real(dp), allocatable :: scale(:)
      real(dp) :: centered, kurtosis, mean_value, second_moment, skewness
      real(dp) :: gamma_scale, gamma_shape, sigma
      integer :: i, n

      n = size(original_observations)
      if (filtered%info /= 0 .or. n < 4 .or. size(filtered%innovation) /= n .or. &
         size(filtered%fitted) /= n .or. lag < 1 .or. lag >= n - 1 .or. &
         fitted_parameter_count < 0 .or. lag < max(1, 3*fitted_parameter_count - 1) .or. &
         any(original_observations <= 0.0_dp)) then
         out%info = 1
         return
      end if
      allocate(out%raw_residual(n), out%transformed_residual(n), out%standardized_residual(n))
      allocate(scale(n))
      if (allocated(filtered%conditional_sd)) then
         if (size(filtered%conditional_sd) /= n .or. any(filtered%conditional_sd <= 0.0_dp)) then
            out%info = 2
            return
         end if
         scale = filtered%conditional_sd
      else
         sigma = sqrt(sum(filtered%innovation**2)/real(n, dp))
         if (present(constant_sigma)) sigma = constant_sigma
         if (sigma <= 0.0_dp) then
            out%info = 2
            return
         end if
         scale = sigma
      end if
      do i = 1, n
         if (abs(lambda) > sqrt(epsilon(1.0_dp)) .and. &
            1.0_dp + lambda*filtered%fitted(i) <= 0.0_dp) then
            out%info = 3
            return
         end if
         out%raw_residual(i) = original_observations(i) - &
            inverse_box_cox_value(filtered%fitted(i), lambda)
      end do
      out%transformed_residual = filtered%innovation
      out%standardized_residual = filtered%innovation/scale
      out%residual_acf = acf_values(out%standardized_residual, lag)
      out%squared_residual_acf = acf_values(out%standardized_residual**2, lag)
      out%randomness = residual_randomness_tests(out%standardized_residual, lag)
      out%arch = mts_arch_test(out%standardized_residual, lag)
      call weighted_ljung_box(out%standardized_residual, lag, fitted_parameter_count, &
         out%weighted_ljung_box, gamma_shape, gamma_scale)
      out%weighted_ljung_box_p_value = &
         regularized_gamma_q(gamma_shape, out%weighted_ljung_box/gamma_scale)
      mean_value = sum(out%standardized_residual)/real(n, dp)
      second_moment = sum((out%standardized_residual - mean_value)**2)/real(n, dp)
      if (second_moment <= tiny(1.0_dp)) then
         out%info = 4
         return
      end if
      skewness = sum((out%standardized_residual - mean_value)**3)/ &
         (real(n, dp)*second_moment**1.5_dp)
      kurtosis = sum((out%standardized_residual - mean_value)**4)/ &
         (real(n, dp)*second_moment**2)
      centered = kurtosis - 3.0_dp
      out%normality_statistic = real(n, dp)*(skewness**2 + 0.25_dp*centered**2)/6.0_dp
      out%normality_p_value = exp(-0.5_dp*out%normality_statistic)
      out%lag = lag
      if (out%randomness%info /= 0 .or. out%arch%info /= 0) out%info = 5
   end function tsissm_diagnose

   pure function tsissm_diagnose_structural(filtered, original_observations, lambda, ar, ma, &
      max_outliers, outlier_alpha, constant_sigma) result(out)
      !! Diagnose fitted structural residuals, roots, stability, and generalized ESD outliers.
      type(tsissm_filter_t), intent(in) :: filtered !! Filtered.
      real(dp), intent(in) :: original_observations(:) !! Original observations.
      real(dp), intent(in) :: lambda !! Penalty or shrinkage parameter.
      real(dp), intent(in) :: ar(:) !! Autoregressive coefficients.
      real(dp), intent(in) :: ma(:) !! Moving-average coefficients.
      integer, intent(in) :: max_outliers !! Maximum outliers.
      real(dp), intent(in) :: outlier_alpha !! Outlier alpha.
      real(dp), intent(in), optional :: constant_sigma !! Constant sigma.
      type(tsissm_structural_diagnostics_t) :: out
      real(dp) :: gamma_scale, gamma_shape
      integer :: degrees, diagnostic_lag, i, j, n

      n = size(original_observations)
      degrees = size(ar) + size(ma)
      j = merge(0, 1, degrees > 0)
      allocate(out%weighted_lag(4), out%weighted_ljung_box(4))
      allocate(out%weighted_ljung_box_p_value(4))
      out%weighted_lag(1) = 1
      out%weighted_lag(2) = max(3*degrees - 1, 1 + degrees + j)
      out%weighted_lag(3) = max(3*degrees - 1, 2 + degrees + j)
      out%weighted_lag(4) = max(3*degrees - 1, 3 + degrees + j)
      diagnostic_lag = maxval(out%weighted_lag)
      if (filtered%info /= 0 .or. n < 4 .or. diagnostic_lag >= n - 1 .or. &
         max_outliers < 0 .or. max_outliers > (n - 1)/2 .or. &
         outlier_alpha <= 0.0_dp .or. outlier_alpha >= 1.0_dp) then
         out%info = 1
         return
      end if
      out%residual = tsissm_diagnose(filtered, original_observations, lambda, diagnostic_lag, &
         degrees, constant_sigma)
      if (out%residual%info /= 0) then
         out%info = 10 + out%residual%info
         return
      end if
      do i = 1, 4
         call weighted_ljung_box(out%residual%standardized_residual, out%weighted_lag(i), &
            merge(0, degrees, i == 1), out%weighted_ljung_box(i), gamma_shape, gamma_scale)
         out%weighted_ljung_box_p_value(i) = &
            regularized_gamma_q(gamma_shape, out%weighted_ljung_box(i)/gamma_scale)
      end do
      call generalized_esd(out%residual%transformed_residual, max_outliers, outlier_alpha, &
         out%outlier_candidate_index, out%outlier_statistic, out%outlier_critical, &
         out%outlier_count)
      out%outlier_index = out%outlier_candidate_index(:out%outlier_count)
      out%ar_inverse_root = inverse_arma_roots(ar, -1.0_dp)
      out%ma_inverse_root = inverse_arma_roots(ma, 1.0_dp)
      out%stability_modulus = abs(matrix_polynomial_roots(filtered%discount_transition))
   end function tsissm_diagnose_structural

   pure subroutine generalized_esd(values, max_outliers, alpha, candidate_index, statistic, &
      critical, outlier_count)
      !! Apply Rosner's generalized ESD test and retain original candidate indices.
      real(dp), intent(in) :: values(:) !! Input values.
      real(dp), intent(in) :: alpha !! Significance, smoothing, or model coefficient.
      integer, intent(in) :: max_outliers !! Maximum outliers.
      integer, allocatable, intent(out) :: candidate_index(:) !! Index of candidate.
      real(dp), allocatable, intent(out) :: statistic(:) !! Statistic.
      real(dp), allocatable, intent(out) :: critical(:) !! Critical.
      integer, intent(out) :: outlier_count !! Number of outlier.
      logical, allocatable :: retained(:)
      real(dp) :: average, probability, scale, t_value
      integer :: candidate, current_count, i

      allocate(candidate_index(max_outliers), statistic(max_outliers), critical(max_outliers))
      allocate(retained(size(values)))
      retained = .true.
      candidate_index = 0
      statistic = 0.0_dp
      critical = 0.0_dp
      outlier_count = 0
      do i = 1, max_outliers
         current_count = count(retained)
         average = sum(pack(values, retained))/real(current_count, dp)
         scale = sqrt(sum(pack((values - average)**2, retained))/real(current_count - 1, dp))
         if (scale <= tiny(1.0_dp)) exit
         candidate = maxloc(abs(values - average), dim=1, mask=retained)
         candidate_index(i) = candidate
         statistic(i) = abs(values(candidate) - average)/scale
         probability = 1.0_dp - alpha/(2.0_dp*real(current_count, dp))
         t_value = student_quantile(probability, real(current_count - 2, dp))
         critical(i) = real(current_count - 1, dp)*t_value/ &
            sqrt(real(current_count, dp)*(real(current_count - 2, dp) + t_value**2))
         retained(candidate) = .false.
         if (statistic(i) > critical(i)) outlier_count = i
      end do
   end subroutine generalized_esd

   pure function inverse_arma_roots(parameters, sign_value) result(inverse_roots)
      !! Return inverse roots of an AR or MA lag polynomial.
      real(dp), intent(in) :: parameters(:) !! Model parameter values.
      real(dp), intent(in) :: sign_value !! Sign value.
      complex(dp), allocatable :: inverse_roots(:)
      real(dp), allocatable :: polynomial(:)
      complex(dp), allocatable :: roots(:)

      allocate(polynomial(size(parameters) + 1))
      polynomial(1) = 1.0_dp
      polynomial(2:) = sign_value*parameters
      roots = polynomial_roots(polynomial)
      allocate(inverse_roots(size(roots)))
      if (size(roots) > 0) inverse_roots = 1.0_dp/roots
   end function inverse_arma_roots

   pure function tsissm_forecast_accuracy(actual, forecast, training, period, alpha) result(out)
      !! Calculate point, interval, CRPS, and kernel log scores for forecast paths.
      real(dp), intent(in) :: actual(:) !! Observed values used for evaluation.
      real(dp), intent(in) :: training(:) !! Training observations.
      real(dp), intent(in) :: alpha !! Significance, smoothing, or model coefficient.
      type(tsissm_forecast_t), intent(in) :: forecast !! Forecast values.
      integer, intent(in) :: period !! Seasonal period.
      type(tsissm_accuracy_t) :: out
      real(dp), allocatable :: lower(:), ordered(:), predicted(:), upper(:)
      real(dp) :: bandwidth, density, difference, scale, standard_deviation
      integer :: h, i, simulation, simulations, valid_percentage

      h = size(actual)
      if (forecast%info /= 0 .or. h < 1 .or. size(forecast%distribution, 2) /= h .or. &
         size(forecast%distribution, 1) < 2 .or. size(training) <= period .or. period < 1 .or. &
         alpha <= 0.0_dp .or. alpha >= 1.0_dp) then
         out%info = 1
         return
      end if
      simulations = size(forecast%distribution, 1)
      predicted = sum(forecast%distribution, dim=1)/real(simulations, dp)
      out%point = forecast_accuracy(actual, predicted, training, period)
      out%mse = sum((actual - predicted)**2)/real(h, dp)
      out%bias = sum(predicted - actual)/real(h, dp)
      scale = sum(abs(training(period + 1:) - training(:size(training) - period)))/ &
         real(size(training) - period, dp)
      if (scale <= tiny(1.0_dp)) out%point%mase = huge(1.0_dp)
      valid_percentage = count(abs(actual) > tiny(1.0_dp))
      if (valid_percentage > 0) then
         out%point%mape = 100.0_dp*sum(pack(abs((actual - predicted)/actual), &
            abs(actual) > tiny(1.0_dp)))/real(valid_percentage, dp)
         out%point%mpe = 100.0_dp*sum(pack((actual - predicted)/actual, &
            abs(actual) > tiny(1.0_dp)))/real(valid_percentage, dp)
      else
         out%point%mape = huge(1.0_dp)
         out%point%mpe = huge(1.0_dp)
      end if
      out%smape = 200.0_dp*sum(abs(actual - predicted)/ &
         max(abs(actual) + abs(predicted), tiny(1.0_dp)))/real(h, dp)
      if (any(actual <= -1.0_dp) .or. any(predicted <= -1.0_dp)) then
         out%mslre = huge(1.0_dp)
      else
         out%mslre = sum((log(1.0_dp + actual) - log(1.0_dp + predicted))**2)/real(h, dp)
      end if
      allocate(lower(h), upper(h), ordered(simulations))
      do i = 1, h
         ordered = forecast%distribution(:, i)
         call sort(ordered)
         lower(i) = quantile(ordered, 0.5_dp*alpha)
         upper(i) = quantile(ordered, 1.0_dp - 0.5_dp*alpha)
         if (actual(i) >= lower(i) .and. actual(i) <= upper(i)) &
            out%interval_coverage = out%interval_coverage + 1.0_dp
         out%interval_score = out%interval_score + upper(i) - lower(i)
         if (actual(i) < lower(i)) then
            out%interval_score = out%interval_score + 2.0_dp*(lower(i) - actual(i))/alpha
         else if (actual(i) > upper(i)) then
            out%interval_score = out%interval_score + 2.0_dp*(actual(i) - upper(i))/alpha
         end if
         out%crps = out%crps + sum(abs(ordered - actual(i)))/real(simulations, dp)
         do simulation = 1, simulations
            out%crps = out%crps - real(2*simulation - simulations - 1, dp)* &
               ordered(simulation)/(real(simulations, dp)*real(simulations, dp))
         end do
         standard_deviation = sqrt(sum((ordered - sum(ordered)/real(simulations, dp))**2)/ &
            real(simulations - 1, dp))
         bandwidth = max(1.06_dp*standard_deviation*real(simulations, dp)**(-0.2_dp), &
            sqrt(epsilon(1.0_dp))*max(1.0_dp, abs(actual(i))))
         density = 0.0_dp
         do simulation = 1, simulations
            difference = (actual(i) - ordered(simulation))/bandwidth
            density = density + exp(-0.5_dp*difference**2)
         end do
         density = density/(real(simulations, dp)*bandwidth*sqrt(2.0_dp*acos(-1.0_dp)))
         out%log_score = out%log_score - log(max(density, tiny(1.0_dp)))
      end do
      out%interval_coverage = out%interval_coverage/real(h, dp)
      out%interval_score = out%interval_score/real(h, dp)
      out%crps = out%crps/real(h, dp)
      out%log_score = out%log_score/real(h, dp)
      out%observations = h
   end function tsissm_forecast_accuracy

   function tsissm_backtest_constant(observations, regressors, transition_base, transition_scale, &
      transition_parameter, observation_loading, initial_persistence, initial_coefficients, &
      initial_state, seed_index, first_origin, horizon, simulations, window_length, &
      reestimate_every, lambda, distribution, skew, shape, alpha, seed, max_iterations, &
      tolerance) result(out)
      !! Run expanding- or fixed-window constant-variance rolling forecasts.
      real(dp), intent(in) :: observations(:) !! Observed time-series values.
      real(dp), intent(in) :: regressors(:, :) !! Regression design matrix.
      real(dp), intent(in) :: transition_base(:, :) !! Transition base.
      real(dp), intent(in) :: transition_scale(:, :) !! Transition scale.
      real(dp), intent(in) :: transition_parameter(:, :) !! Transition parameter.
      real(dp), intent(in) :: observation_loading(:) !! Observation loading matrix.
      real(dp), intent(in) :: initial_persistence(:) !! Initial persistence.
      real(dp), intent(in) :: initial_coefficients(:) !! Initial coefficients.
      real(dp), intent(in) :: initial_state(:) !! Initial state vector.
      integer, intent(in) :: seed_index(:) !! Index of seed.
      integer, intent(in) :: first_origin !! First origin.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      integer, intent(in) :: simulations !! Number of simulation draws.
      integer, intent(in) :: window_length !! Window length.
      integer, intent(in) :: reestimate_every !! Reestimate every.
      integer, intent(in) :: distribution !! Probability-distribution specification.
      integer, intent(in) :: max_iterations !! Maximum number of algorithm iterations.
      real(dp), intent(in) :: lambda !! Penalty or shrinkage parameter.
      real(dp), intent(in) :: skew !! Skew.
      real(dp), intent(in) :: shape !! Shape.
      real(dp), intent(in) :: alpha !! Significance, smoothing, or model coefficient.
      real(dp), intent(in) :: tolerance !! Numerical convergence tolerance.
      integer, intent(in), optional :: seed !! Random-number seed.
      type(tsissm_backtest_t) :: out
      type(tsissm_filter_t) :: filtered
      type(tsissm_fit_t) :: fitted
      type(tsissm_forecast_t) :: forecast
      real(dp), allocatable :: coefficients(:), persistence(:), training_regressors(:, :)
      real(dp) :: sigma
      integer :: available, current, effective_horizon, forecast_seed, n, origin_index, origins, start
      logical :: reestimate_now

      n = size(observations)
      if (.not. valid_backtest_inputs(observations, regressors, observation_loading, &
         initial_persistence, initial_coefficients, initial_state, first_origin, horizon, &
         simulations, window_length, reestimate_every, alpha, max_iterations, tolerance)) then
         out%info = 1
         return
      end if
      origins = n - first_origin
      effective_horizon = min(horizon, origins)
      call allocate_backtest(out, origins, effective_horizon, size(initial_persistence), &
         size(initial_coefficients))
      persistence = initial_persistence
      coefficients = initial_coefficients
      do origin_index = 1, origins
         current = first_origin + origin_index - 1
         start = 1
         if (window_length > 0) start = max(1, current - window_length + 1)
         training_regressors = regressors(start:current, :)
         out%training_start(origin_index) = start
         out%training_end(origin_index) = current
         reestimate_now = .false.
         if (reestimate_every > 0) then
            reestimate_now = origin_index == 1 .or. mod(origin_index - 1, reestimate_every) == 0
         end if
         if (reestimate_now) then
            fitted = tsissm_fit_constant(observations(start:current), transition_base, &
               transition_scale, transition_parameter, observation_loading, persistence, &
               training_regressors, coefficients, initial_state, seed_index, lambda=lambda, &
               distribution=distribution, skew=skew, shape=shape, max_iterations=max_iterations, &
               tolerance=tolerance)
            if (fitted%info < 0 .or. fitted%info > 4) then
               out%info = 10 + origin_index
               return
            end if
            persistence = fitted%persistence
            coefficients = fitted%coefficients
            filtered = fitted%filter
            sigma = fitted%likelihood%sigma
            out%reestimated(origin_index) = .true.
         else
            filtered = tsissm_filter_constant(observations(start:current), transition_base, &
               transition_scale, transition_parameter, observation_loading, persistence, &
               training_regressors, coefficients, initial_state, lambda=lambda)
            if (filtered%info /= 0) then
               out%info = 20 + origin_index
               return
            end if
            sigma = sqrt(sum(filtered%innovation**2)/real(size(filtered%innovation), dp))
         end if
         if (sigma <= tiny(1.0_dp)) sigma = sqrt(epsilon(1.0_dp))
         out%persistence(origin_index, :) = persistence
         out%coefficients(origin_index, :) = coefficients
         available = min(effective_horizon, n - current)
         if (present(seed)) then
            forecast_seed = seed + origin_index - 1
            forecast = tsissm_simulate_constant(transition_base, transition_scale, &
               transition_parameter, observation_loading, persistence, &
               regressors(current + 1:current + available, :), coefficients, &
               filtered%state(:, size(filtered%state, 2) - 1), sigma, lambda, distribution, skew, &
               shape, simulations, [0.5_dp*alpha, 1.0_dp - 0.5_dp*alpha], forecast_seed)
         else
            forecast = tsissm_simulate_constant(transition_base, transition_scale, &
               transition_parameter, observation_loading, persistence, &
               regressors(current + 1:current + available, :), coefficients, &
               filtered%state(:, size(filtered%state, 2) - 1), sigma, lambda, distribution, skew, &
               shape, simulations, [0.5_dp*alpha, 1.0_dp - 0.5_dp*alpha])
         end if
         if (forecast%info /= 0) then
            out%info = 30 + origin_index
            return
         end if
         call store_backtest_origin(out, origin_index, observations(current + 1:current + available), &
            forecast)
      end do
      call summarize_backtest(out, alpha)
      out%first_origin = first_origin
      out%horizon = effective_horizon
      out%window_length = window_length
      out%reestimate_every = reestimate_every
   end function tsissm_backtest_constant

   function tsissm_backtest_dynamic(observations, regressors, transition_base, transition_scale, &
      transition_parameter, observation_loading, initial_persistence, initial_coefficients, &
      initial_arch, initial_garch, initial_state, seed_index, initial_arch_history, &
      initial_variance_history, variance_intercept, first_origin, horizon, simulations, &
      window_length, reestimate_every, lambda, distribution, skew, shape, alpha, seed, &
      max_iterations, tolerance, variance_initialization, variance_sample_size) result(out)
      !! Run expanding- or fixed-window GARCH rolling forecasts.
      real(dp), intent(in) :: observations(:) !! Observed time-series values.
      real(dp), intent(in) :: regressors(:, :) !! Regression design matrix.
      real(dp), intent(in) :: transition_base(:, :) !! Transition base.
      real(dp), intent(in) :: transition_scale(:, :) !! Transition scale.
      real(dp), intent(in) :: transition_parameter(:, :) !! Transition parameter.
      real(dp), intent(in) :: observation_loading(:) !! Observation loading matrix.
      real(dp), intent(in) :: initial_persistence(:) !! Initial persistence.
      real(dp), intent(in) :: initial_coefficients(:) !! Initial coefficients.
      real(dp), intent(in) :: initial_arch(:) !! Initial arch.
      real(dp), intent(in) :: initial_garch(:) !! Initial GARCH.
      real(dp), intent(in) :: initial_state(:) !! Initial state vector.
      real(dp), intent(in) :: initial_arch_history(:) !! Initial arch history.
      real(dp), intent(in) :: initial_variance_history(:) !! Initial variance history.
      real(dp), intent(in) :: variance_intercept !! Variance intercept.
      real(dp), intent(in) :: lambda !! Penalty or shrinkage parameter.
      real(dp), intent(in) :: skew !! Skew.
      real(dp), intent(in) :: shape !! Shape.
      real(dp), intent(in) :: alpha !! Significance, smoothing, or model coefficient.
      real(dp), intent(in) :: tolerance !! Numerical convergence tolerance.
      integer, intent(in) :: seed_index(:) !! Index of seed.
      integer, intent(in) :: first_origin !! First origin.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      integer, intent(in) :: simulations !! Number of simulation draws.
      integer, intent(in) :: window_length !! Window length.
      integer, intent(in) :: reestimate_every !! Reestimate every.
      integer, intent(in) :: distribution !! Probability-distribution specification.
      integer, intent(in) :: max_iterations !! Maximum number of algorithm iterations.
      integer, intent(in), optional :: seed !! Random-number seed.
      integer, intent(in), optional :: variance_initialization !! Variance initialization.
      integer, intent(in), optional :: variance_sample_size !! Variance sample size.
      type(tsissm_backtest_t) :: out
      type(tsissm_filter_t) :: filtered
      type(tsissm_fit_t) :: fitted
      type(tsissm_forecast_t) :: forecast
      real(dp), allocatable :: arch(:), arch_history(:), coefficients(:), garch(:)
      real(dp), allocatable :: persistence(:), training_regressors(:, :), variance_history(:)
      real(dp) :: intercept
      integer :: available, current, effective_horizon, forecast_seed, history
      integer :: n, origin_index, origins, start
      logical :: reestimate_now

      n = size(observations)
      history = max(size(initial_arch), size(initial_garch))
      if (.not. valid_backtest_inputs(observations, regressors, observation_loading, &
         initial_persistence, initial_coefficients, initial_state, first_origin, horizon, &
         simulations, window_length, reestimate_every, alpha, max_iterations, tolerance) .or. &
         history < 1 .or. first_origin < history .or. &
         (window_length > 0 .and. window_length < history) .or. &
         size(initial_arch_history) < history .or. &
         size(initial_variance_history) < history .or. variance_intercept < 0.0_dp) then
         out%info = 1
         return
      end if
      origins = n - first_origin
      effective_horizon = min(horizon, origins)
      call allocate_backtest(out, origins, effective_horizon, size(initial_persistence), &
         size(initial_coefficients))
      persistence = initial_persistence
      coefficients = initial_coefficients
      arch = initial_arch
      garch = initial_garch
      intercept = variance_intercept
      do origin_index = 1, origins
         current = first_origin + origin_index - 1
         start = 1
         if (window_length > 0) start = max(1, current - window_length + 1)
         training_regressors = regressors(start:current, :)
         out%training_start(origin_index) = start
         out%training_end(origin_index) = current
         reestimate_now = .false.
         if (reestimate_every > 0) then
            reestimate_now = origin_index == 1 .or. mod(origin_index - 1, reestimate_every) == 0
         end if
         if (reestimate_now) then
            fitted = tsissm_fit_dynamic(observations(start:current), transition_base, &
               transition_scale, transition_parameter, observation_loading, persistence, &
               training_regressors, coefficients, arch, garch, initial_state, seed_index, &
               lambda=lambda, distribution=distribution, skew=skew, shape=shape, &
               max_iterations=max_iterations, tolerance=tolerance, &
               variance_initialization=variance_initialization, &
               variance_sample_size=variance_sample_size)
            if (fitted%info < 0 .or. fitted%info > 4) then
               out%info = 10 + origin_index
               return
            end if
            persistence = fitted%persistence
            coefficients = fitted%coefficients
            arch = fitted%arch
            garch = fitted%garch
            intercept = fitted%variance_intercept
            filtered = fitted%filter
            out%reestimated(origin_index) = .true.
         else
            filtered = tsissm_filter_dynamic(observations(start:current), transition_base, &
               transition_scale, transition_parameter, observation_loading, persistence, &
               training_regressors, coefficients, initial_state, lambda=lambda, arch=arch, &
               garch=garch, initial_arch=initial_arch_history, &
               initial_variance=initial_variance_history, variance_intercept=intercept)
            if (filtered%info /= 0) then
               out%info = 20 + origin_index
               return
            end if
         end if
         out%persistence(origin_index, :) = persistence
         out%coefficients(origin_index, :) = coefficients
         arch_history = filtered%innovation(size(filtered%innovation) - history + 1:)**2
         variance_history = filtered%conditional_sd(size(filtered%conditional_sd) - history + 1:)**2
         available = min(effective_horizon, n - current)
         if (present(seed)) then
            forecast_seed = seed + origin_index - 1
            forecast = tsissm_simulate_dynamic(transition_base, transition_scale, &
               transition_parameter, observation_loading, persistence, &
               regressors(current + 1:current + available, :), coefficients, &
               filtered%state(:, size(filtered%state, 2) - 1), arch, garch, intercept, &
               arch_history, variance_history, lambda, distribution, skew, shape, simulations, &
               [0.5_dp*alpha, 1.0_dp - 0.5_dp*alpha], forecast_seed)
         else
            forecast = tsissm_simulate_dynamic(transition_base, transition_scale, &
               transition_parameter, observation_loading, persistence, &
               regressors(current + 1:current + available, :), coefficients, &
               filtered%state(:, size(filtered%state, 2) - 1), arch, garch, intercept, &
               arch_history, variance_history, lambda, distribution, skew, shape, simulations, &
               [0.5_dp*alpha, 1.0_dp - 0.5_dp*alpha])
         end if
         if (forecast%info /= 0) then
            out%info = 30 + origin_index
            return
         end if
         call store_backtest_origin(out, origin_index, observations(current + 1:current + available), &
            forecast)
      end do
      call summarize_backtest(out, alpha)
      out%first_origin = first_origin
      out%horizon = effective_horizon
      out%window_length = window_length
      out%reestimate_every = reestimate_every
   end function tsissm_backtest_dynamic

   function tsissm_backtest_structural_constant(observations, regressors, slope, &
      damped_slope, seasonal_frequency, seasonal_harmonics, regular_seasonal, ar_order, &
      ma_order, initial_state, initial_parameters, first_origin, horizon, simulations, &
      window_length, reestimate_every, distribution, alpha, max_iterations, tolerance, &
      seed) result(out)
      !! Run rolling forecasts with joint structural constant-variance refits.
      real(dp), intent(in) :: observations(:) !! Observed time-series values.
      real(dp), intent(in) :: regressors(:, :) !! Regression design matrix.
      real(dp), intent(in) :: seasonal_frequency(:) !! Seasonal frequency.
      real(dp), intent(in) :: initial_state(:) !! Initial state vector.
      real(dp), intent(in) :: initial_parameters(:) !! Initial parameter values.
      real(dp), intent(in) :: alpha !! Significance, smoothing, or model coefficient.
      real(dp), intent(in) :: tolerance !! Numerical convergence tolerance.
      logical, intent(in) :: slope !! Flag controlling slope.
      logical, intent(in) :: damped_slope !! Flag controlling damped slope.
      logical, intent(in) :: regular_seasonal !! Flag controlling regular seasonal.
      integer, intent(in) :: seasonal_harmonics(:) !! Seasonal harmonics.
      integer, intent(in) :: ar_order !! Autoregressive order.
      integer, intent(in) :: ma_order !! Moving-average order.
      integer, intent(in) :: first_origin !! First origin.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      integer, intent(in) :: simulations !! Number of simulation draws.
      integer, intent(in) :: window_length !! Window length.
      integer, intent(in) :: reestimate_every !! Reestimate every.
      integer, intent(in) :: distribution !! Probability-distribution specification.
      integer, intent(in) :: max_iterations !! Maximum number of algorithm iterations.
      integer, intent(in), optional :: seed !! Random-number seed.
      type(tsissm_backtest_t) :: out
      type(tsissm_filter_t) :: filtered
      type(tsissm_forecast_t) :: forecast
      type(tsissm_model_t) :: model
      type(tsissm_structural_fit_t) :: fitted
      real(dp), allocatable :: coefficients(:), parameters(:), training_regressors(:, :)
      real(dp) :: lambda, shape, sigma, skew
      integer :: available, current, effective_horizon, forecast_seed, n
      integer :: origin_index, origins, start, status
      logical :: reestimate_now

      n = size(observations)
      call decode_structural_parameters(slope, damped_slope, seasonal_frequency, &
         seasonal_harmonics, regular_seasonal, ar_order, ma_order, size(regressors, 2), &
         initial_parameters, distribution, model, coefficients, lambda, skew, shape, status)
      if (.not. valid_structural_backtest_inputs(observations, regressors, initial_state, &
         initial_parameters, model, status, first_origin, horizon, simulations, window_length, &
         reestimate_every, alpha, max_iterations, tolerance)) then
         out%info = 1
         return
      end if
      origins = n - first_origin
      effective_horizon = min(horizon, origins)
      call allocate_backtest(out, origins, effective_horizon, size(model%persistence), &
         size(coefficients), size(initial_parameters))
      parameters = initial_parameters
      do origin_index = 1, origins
         current = first_origin + origin_index - 1
         start = 1
         if (window_length > 0) start = max(1, current - window_length + 1)
         training_regressors = regressors(start:current, :)
         out%training_start(origin_index) = start
         out%training_end(origin_index) = current
         reestimate_now = reestimate_every > 0 .and. &
            (origin_index == 1 .or. mod(origin_index - 1, reestimate_every) == 0)
         if (reestimate_now) then
            fitted = tsissm_fit_structural_constant(observations(start:current), slope, &
               damped_slope, seasonal_frequency, seasonal_harmonics, regular_seasonal, &
               ar_order, ma_order, training_regressors, initial_state, parameters, &
               distribution=distribution, max_iterations=max_iterations, tolerance=tolerance)
            if (fitted%info < 0 .or. fitted%info > 4) then
               out%info = 10 + origin_index
               return
            end if
            parameters = fitted%parameters
            model = fitted%model
            filtered = fitted%filter
            lambda = fitted%lambda
            skew = fitted%skew
            shape = fitted%shape
            call structural_coefficients(parameters, slope, damped_slope, &
               size(seasonal_frequency), ar_order, ma_order, size(regressors, 2), coefficients)
            sigma = fitted%likelihood%sigma
            out%reestimated(origin_index) = .true.
         else
            call decode_structural_parameters(slope, damped_slope, seasonal_frequency, &
               seasonal_harmonics, regular_seasonal, ar_order, ma_order, size(regressors, 2), &
               parameters, distribution, model, coefficients, lambda, skew, shape, status)
            if (status /= 0) then
               out%info = 20 + origin_index
               return
            end if
            filtered = tsissm_filter_constant(observations(start:current), &
               model%transition_base, model%transition_scale, model%transition_parameter, &
               model%observation_loading, model%persistence, training_regressors, coefficients, &
               initial_state, lambda=lambda)
            if (filtered%info /= 0) then
               out%info = 20 + origin_index
               return
            end if
            sigma = sqrt(sum(filtered%innovation**2)/real(size(filtered%innovation), dp))
         end if
         sigma = max(sigma, sqrt(epsilon(1.0_dp)))
         out%parameters(origin_index, :) = parameters
         out%persistence(origin_index, :) = model%persistence
         out%coefficients(origin_index, :) = coefficients
         available = min(effective_horizon, n - current)
         if (present(seed)) then
            forecast_seed = seed + origin_index - 1
            forecast = tsissm_simulate_constant(model%transition_base, model%transition_scale, &
               model%transition_parameter, model%observation_loading, model%persistence, &
               regressors(current + 1:current + available, :), coefficients, &
               filtered%state(:, ubound(filtered%state, 2)), sigma, lambda, distribution, &
               skew, shape, simulations, [0.5_dp*alpha, 1.0_dp - 0.5_dp*alpha], forecast_seed)
         else
            forecast = tsissm_simulate_constant(model%transition_base, model%transition_scale, &
               model%transition_parameter, model%observation_loading, model%persistence, &
               regressors(current + 1:current + available, :), coefficients, &
               filtered%state(:, ubound(filtered%state, 2)), sigma, lambda, distribution, &
               skew, shape, simulations, [0.5_dp*alpha, 1.0_dp - 0.5_dp*alpha])
         end if
         if (forecast%info /= 0) then
            out%info = 30 + origin_index
            return
         end if
         call store_backtest_origin(out, origin_index, &
            observations(current + 1:current + available), forecast)
      end do
      call summarize_backtest(out, alpha)
      out%first_origin = first_origin
      out%horizon = effective_horizon
      out%window_length = window_length
      out%reestimate_every = reestimate_every
   end function tsissm_backtest_structural_constant

   function tsissm_backtest_structural_dynamic(observations, regressors, slope, &
      damped_slope, seasonal_frequency, seasonal_harmonics, regular_seasonal, ar_order, &
      ma_order, initial_state, initial_parameters, initial_arch, initial_garch, &
      first_origin, horizon, simulations, window_length, reestimate_every, distribution, &
      alpha, max_iterations, tolerance, seed, variance_initialization, &
      variance_sample_size) result(out)
      !! Run rolling forecasts with joint structural and GARCH refits.
      real(dp), intent(in) :: observations(:) !! Observed time-series values.
      real(dp), intent(in) :: regressors(:, :) !! Regression design matrix.
      real(dp), intent(in) :: seasonal_frequency(:) !! Seasonal frequency.
      real(dp), intent(in) :: initial_state(:) !! Initial state vector.
      real(dp), intent(in) :: initial_parameters(:) !! Initial parameter values.
      real(dp), intent(in) :: initial_arch(:) !! Initial arch.
      real(dp), intent(in) :: initial_garch(:) !! Initial GARCH.
      real(dp), intent(in) :: alpha !! Significance, smoothing, or model coefficient.
      real(dp), intent(in) :: tolerance !! Numerical convergence tolerance.
      logical, intent(in) :: slope !! Flag controlling slope.
      logical, intent(in) :: damped_slope !! Flag controlling damped slope.
      logical, intent(in) :: regular_seasonal !! Flag controlling regular seasonal.
      integer, intent(in) :: seasonal_harmonics(:) !! Seasonal harmonics.
      integer, intent(in) :: ar_order !! Autoregressive order.
      integer, intent(in) :: ma_order !! Moving-average order.
      integer, intent(in) :: first_origin !! First origin.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      integer, intent(in) :: simulations !! Number of simulation draws.
      integer, intent(in) :: window_length !! Window length.
      integer, intent(in) :: reestimate_every !! Reestimate every.
      integer, intent(in) :: distribution !! Probability-distribution specification.
      integer, intent(in) :: max_iterations !! Maximum number of algorithm iterations.
      integer, intent(in), optional :: seed !! Random-number seed.
      integer, intent(in), optional :: variance_initialization !! Variance initialization.
      integer, intent(in), optional :: variance_sample_size !! Variance sample size.
      type(tsissm_backtest_t) :: out
      type(tsissm_filter_t) :: constant_filter, filtered
      type(tsissm_forecast_t) :: forecast
      type(tsissm_model_t) :: model
      type(tsissm_structural_fit_t) :: fitted
      type(tsissm_variance_initialization_t) :: initialized_variance
      real(dp), allocatable :: arch(:), arch_history(:), coefficients(:), garch(:)
      real(dp), allocatable :: initial_arch_history(:), initial_variance_history(:)
      real(dp), allocatable :: parameters(:), training_regressors(:, :), variance_history(:)
      real(dp) :: intercept, lambda, shape, skew, variance
      integer :: available, current, effective_horizon, forecast_seed, history, n
      integer :: origin_index, origins, start, status, structural_count
      logical :: reestimate_now

      n = size(observations)
      history = max(size(initial_arch), size(initial_garch))
      structural_count = size(initial_parameters)
      call decode_structural_parameters(slope, damped_slope, seasonal_frequency, &
         seasonal_harmonics, regular_seasonal, ar_order, ma_order, size(regressors, 2), &
         initial_parameters, distribution, model, coefficients, lambda, skew, shape, status)
      if (.not. valid_structural_backtest_inputs(observations, regressors, initial_state, &
         initial_parameters, model, status, first_origin, horizon, simulations, window_length, &
         reestimate_every, alpha, max_iterations, tolerance) .or. history < 1 .or. &
         first_origin < history .or. (window_length > 0 .and. window_length < history) .or. &
         any(initial_arch < 0.0_dp) .or. any(initial_garch < 0.0_dp)) then
         out%info = 1
         return
      end if
      origins = n - first_origin
      effective_horizon = min(horizon, origins)
      call allocate_backtest(out, origins, effective_horizon, size(model%persistence), &
         size(coefficients), structural_count + size(initial_arch) + size(initial_garch))
      parameters = [initial_parameters, initial_arch, initial_garch]
      allocate(initial_arch_history(history), initial_variance_history(history))
      do origin_index = 1, origins
         current = first_origin + origin_index - 1
         start = 1
         if (window_length > 0) start = max(1, current - window_length + 1)
         training_regressors = regressors(start:current, :)
         out%training_start(origin_index) = start
         out%training_end(origin_index) = current
         reestimate_now = reestimate_every > 0 .and. &
            (origin_index == 1 .or. mod(origin_index - 1, reestimate_every) == 0)
         if (reestimate_now) then
            fitted = tsissm_fit_structural_dynamic(observations(start:current), slope, &
               damped_slope, seasonal_frequency, seasonal_harmonics, regular_seasonal, &
               ar_order, ma_order, training_regressors, initial_state, &
               parameters(:structural_count), &
               parameters(structural_count + 1:structural_count + size(initial_arch)), &
               parameters(structural_count + size(initial_arch) + 1:), &
               distribution=distribution, max_iterations=max_iterations, tolerance=tolerance, &
               variance_initialization=variance_initialization, &
               variance_sample_size=variance_sample_size)
            if (fitted%info < 0 .or. fitted%info > 4) then
               out%info = 10 + origin_index
               return
            end if
            parameters = fitted%parameters
            model = fitted%model
            filtered = fitted%filter
            arch = fitted%arch
            garch = fitted%garch
            intercept = fitted%variance_intercept
            lambda = fitted%lambda
            skew = fitted%skew
            shape = fitted%shape
            call structural_coefficients(parameters(:structural_count), slope, damped_slope, &
               size(seasonal_frequency), ar_order, ma_order, size(regressors, 2), coefficients)
            out%reestimated(origin_index) = .true.
         else
            call decode_structural_parameters(slope, damped_slope, seasonal_frequency, &
               seasonal_harmonics, regular_seasonal, ar_order, ma_order, size(regressors, 2), &
               parameters(:structural_count), distribution, model, coefficients, lambda, &
               skew, shape, status)
            if (status /= 0) then
               out%info = 20 + origin_index
               return
            end if
            arch = parameters(structural_count + 1: &
               structural_count + size(initial_arch))
            garch = parameters(structural_count + size(initial_arch) + 1:)
            constant_filter = tsissm_filter_constant(observations(start:current), &
               model%transition_base, model%transition_scale, model%transition_parameter, &
               model%observation_loading, model%persistence, training_regressors, coefficients, &
               initial_state, lambda=lambda)
            if (constant_filter%info /= 0) then
               out%info = 20 + origin_index
               return
            end if
            initialized_variance = tsissm_initialize_variance(constant_filter%innovation, &
               method=variance_initialization, sample_size=variance_sample_size)
            if (initialized_variance%info /= 0) then
               out%info = 20 + origin_index
               return
            end if
            variance = initialized_variance%variance
            intercept = max(1.0e-12_dp, variance*(1.0_dp - sum(arch) - sum(garch)))
            initial_arch_history = variance
            initial_variance_history = variance
            filtered = tsissm_filter_dynamic(observations(start:current), &
               model%transition_base, model%transition_scale, model%transition_parameter, &
               model%observation_loading, model%persistence, training_regressors, coefficients, &
               initial_state, lambda=lambda, arch=arch, garch=garch, &
               initial_arch=initial_arch_history, initial_variance=initial_variance_history, &
               variance_intercept=intercept)
            if (filtered%info /= 0) then
               out%info = 20 + origin_index
               return
            end if
         end if
         out%parameters(origin_index, :) = parameters
         out%persistence(origin_index, :) = model%persistence
         out%coefficients(origin_index, :) = coefficients
         arch_history = filtered%innovation(size(filtered%innovation) - history + 1:)**2
         variance_history = filtered%conditional_sd( &
            size(filtered%conditional_sd) - history + 1:)**2
         available = min(effective_horizon, n - current)
         if (present(seed)) then
            forecast_seed = seed + origin_index - 1
            forecast = tsissm_simulate_dynamic(model%transition_base, model%transition_scale, &
               model%transition_parameter, model%observation_loading, model%persistence, &
               regressors(current + 1:current + available, :), coefficients, &
               filtered%state(:, ubound(filtered%state, 2)), arch, garch, intercept, &
               arch_history, variance_history, lambda, distribution, skew, shape, simulations, &
               [0.5_dp*alpha, 1.0_dp - 0.5_dp*alpha], forecast_seed)
         else
            forecast = tsissm_simulate_dynamic(model%transition_base, model%transition_scale, &
               model%transition_parameter, model%observation_loading, model%persistence, &
               regressors(current + 1:current + available, :), coefficients, &
               filtered%state(:, ubound(filtered%state, 2)), arch, garch, intercept, &
               arch_history, variance_history, lambda, distribution, skew, shape, simulations, &
               [0.5_dp*alpha, 1.0_dp - 0.5_dp*alpha])
         end if
         if (forecast%info /= 0) then
            out%info = 30 + origin_index
            return
         end if
         call store_backtest_origin(out, origin_index, &
            observations(current + 1:current + available), forecast)
      end do
      call summarize_backtest(out, alpha)
      out%first_origin = first_origin
      out%horizon = effective_horizon
      out%window_length = window_length
      out%reestimate_every = reestimate_every
   end function tsissm_backtest_structural_dynamic

   pure subroutine decode_structural_parameters(slope, damped_slope, seasonal_frequency, &
      seasonal_harmonics, regular_seasonal, ar_order, ma_order, regression_count, parameters, &
      distribution, model, coefficients, lambda, skew, shape, status)
      !! Decode the common flat structural parameter representation.
      logical, intent(in) :: slope !! Flag controlling slope.
      logical, intent(in) :: damped_slope !! Flag controlling damped slope.
      logical, intent(in) :: regular_seasonal !! Flag controlling regular seasonal.
      real(dp), intent(in) :: seasonal_frequency(:) !! Seasonal frequency.
      real(dp), intent(in) :: parameters(:) !! Model parameter values.
      integer, intent(in) :: seasonal_harmonics(:) !! Seasonal harmonics.
      integer, intent(in) :: ar_order !! Autoregressive order.
      integer, intent(in) :: ma_order !! Moving-average order.
      integer, intent(in) :: regression_count !! Number of regression.
      integer, intent(in) :: distribution !! Probability-distribution specification.
      type(tsissm_model_t), intent(out) :: model !! Model specification.
      real(dp), allocatable, intent(out) :: coefficients(:) !! Model coefficients.
      real(dp), intent(out) :: lambda !! Penalty or shrinkage parameter.
      real(dp), intent(out) :: skew !! Skew.
      real(dp), intent(out) :: shape !! Shape.
      integer, intent(out) :: status !! Status.
      real(dp), allocatable :: ar(:), ma(:), seasonal_persistence(:)
      real(dp) :: alpha, beta, damping
      integer :: expected, index

      expected = structural_parameter_count(slope, damped_slope, size(seasonal_frequency), &
         ar_order, ma_order, regression_count, distribution)
      status = 1
      if (size(parameters) /= expected .or. ar_order < 0 .or. ma_order < 0 .or. &
         regression_count < 0 .or. size(seasonal_harmonics) /= size(seasonal_frequency)) return
      index = 1
      alpha = parameters(index)
      index = index + 1
      beta = 0.0_dp
      if (slope) then
         beta = parameters(index)
         index = index + 1
      end if
      damping = 1.0_dp
      if (damped_slope) then
         damping = parameters(index)
         index = index + 1
      end if
      seasonal_persistence = parameters(index:index + size(seasonal_frequency) - 1)
      index = index + size(seasonal_frequency)
      ar = parameters(index:index + ar_order - 1)
      index = index + ar_order
      ma = parameters(index:index + ma_order - 1)
      index = index + ma_order
      coefficients = parameters(index:index + regression_count - 1)
      index = index + regression_count
      lambda = parameters(index)
      index = index + 1
      skew = 0.0_dp
      shape = 0.0_dp
      if (distribution == tsissm_distribution_student) then
         shape = parameters(index)
      else if (distribution == tsissm_distribution_johnson_su) then
         skew = parameters(index)
         shape = parameters(index + 1)
      end if
      if (.not. valid_distribution_parameters(distribution, shape)) return
      model = tsissm_model(slope, damped_slope, seasonal_frequency, seasonal_harmonics, &
         regular_seasonal, ar, ma, alpha, beta, seasonal_persistence, damping)
      status = model%info
   end subroutine decode_structural_parameters

   pure subroutine structural_coefficients(parameters, slope, damped_slope, seasonal_count, &
      ar_order, ma_order, regression_count, coefficients)
      !! Extract regression coefficients from a structural parameter vector.
      real(dp), intent(in) :: parameters(:) !! Model parameter values.
      logical, intent(in) :: slope !! Flag controlling slope.
      logical, intent(in) :: damped_slope !! Flag controlling damped slope.
      integer, intent(in) :: seasonal_count !! Number of seasonal.
      integer, intent(in) :: ar_order !! Autoregressive order.
      integer, intent(in) :: ma_order !! Moving-average order.
      integer, intent(in) :: regression_count !! Number of regression.
      real(dp), allocatable, intent(out) :: coefficients(:) !! Model coefficients.
      integer :: first

      first = 2 + merge(1, 0, slope) + merge(1, 0, damped_slope) + &
         seasonal_count + ar_order + ma_order
      coefficients = parameters(first:first + regression_count - 1)
   end subroutine structural_coefficients

   pure integer function structural_parameter_count(slope, damped_slope, seasonal_count, &
      ar_order, ma_order, regression_count, distribution) result(count_parameters)
      !! Count parameters in the common flat structural representation.
      logical, intent(in) :: slope !! Flag controlling slope.
      logical, intent(in) :: damped_slope !! Flag controlling damped slope.
      integer, intent(in) :: seasonal_count !! Number of seasonal.
      integer, intent(in) :: ar_order !! Autoregressive order.
      integer, intent(in) :: ma_order !! Moving-average order.
      integer, intent(in) :: regression_count !! Number of regression.
      integer, intent(in) :: distribution !! Probability-distribution specification.

      count_parameters = 1 + merge(1, 0, slope) + merge(1, 0, damped_slope) + &
         seasonal_count + ar_order + ma_order + regression_count + 1
      if (distribution == tsissm_distribution_student) count_parameters = count_parameters + 1
      if (distribution == tsissm_distribution_johnson_su) count_parameters = count_parameters + 2
   end function structural_parameter_count

   pure logical function valid_structural_backtest_inputs(observations, regressors, initial_state, &
      parameters, model, model_status, first_origin, horizon, simulations, window_length, &
      reestimate_every, alpha, max_iterations, tolerance) result(valid)
      !! Check dimensions and controls for joint structural rolling backtests.
      real(dp), intent(in) :: observations(:) !! Observed time-series values.
      real(dp), intent(in) :: regressors(:, :) !! Regression design matrix.
      real(dp), intent(in) :: initial_state(:) !! Initial state vector.
      real(dp), intent(in) :: parameters(:) !! Model parameter values.
      type(tsissm_model_t), intent(in) :: model !! Model specification.
      integer, intent(in) :: model_status !! Model status.
      integer, intent(in) :: first_origin !! First origin.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      integer, intent(in) :: simulations !! Number of simulation draws.
      integer, intent(in) :: window_length !! Window length.
      integer, intent(in) :: reestimate_every !! Reestimate every.
      integer, intent(in) :: max_iterations !! Maximum number of algorithm iterations.
      real(dp), intent(in) :: alpha !! Significance, smoothing, or model coefficient.
      real(dp), intent(in) :: tolerance !! Numerical convergence tolerance.

      valid = .false.
      if (model_status /= 0) return
      valid = size(observations) >= 3 .and. all(observations > 0.0_dp) .and. &
         size(regressors, 1) == size(observations) .and. size(parameters) > 0 .and. &
         size(initial_state) == size(model%initial_state) .and. first_origin >= 2 .and. &
         first_origin < size(observations) .and. horizon >= 1 .and. simulations >= 2 .and. &
         window_length >= 0 .and. (window_length == 0 .or. window_length >= 2) .and. &
         reestimate_every >= 0 .and. alpha > 0.0_dp .and. alpha < 1.0_dp .and. &
         max_iterations >= 1 .and. tolerance > 0.0_dp
   end function valid_structural_backtest_inputs

   pure logical function valid_backtest_inputs(observations, regressors, observation_loading, &
      persistence, coefficients, initial_state, first_origin, horizon, simulations, window_length, &
      reestimate_every, alpha, max_iterations, tolerance) result(valid)
      !! Check dimensions and common rolling-backtest controls.
      real(dp), intent(in) :: observations(:) !! Observed time-series values.
      real(dp), intent(in) :: regressors(:, :) !! Regression design matrix.
      real(dp), intent(in) :: observation_loading(:) !! Observation loading matrix.
      real(dp), intent(in) :: persistence(:) !! Persistence.
      real(dp), intent(in) :: coefficients(:) !! Model coefficients.
      real(dp), intent(in) :: initial_state(:) !! Initial state vector.
      real(dp), intent(in) :: alpha !! Significance, smoothing, or model coefficient.
      real(dp), intent(in) :: tolerance !! Numerical convergence tolerance.
      integer, intent(in) :: first_origin !! First origin.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      integer, intent(in) :: simulations !! Number of simulation draws.
      integer, intent(in) :: window_length !! Window length.
      integer, intent(in) :: reestimate_every !! Reestimate every.
      integer, intent(in) :: max_iterations !! Maximum number of algorithm iterations.

      valid = size(observations) >= 3 .and. all(observations > 0.0_dp) .and. &
         size(regressors, 1) == size(observations) .and. &
         size(regressors, 2) == size(coefficients) .and. &
         size(observation_loading) == size(initial_state) .and. &
         size(persistence) == size(initial_state) .and. first_origin >= 2 .and. &
         first_origin < size(observations) .and. horizon >= 1 .and. simulations >= 2 .and. &
         window_length >= 0 .and. (window_length == 0 .or. window_length >= 2) .and. &
         reestimate_every >= 0 .and. alpha > 0.0_dp .and. alpha < 1.0_dp .and. &
         max_iterations >= 1 .and. tolerance > 0.0_dp
   end function valid_backtest_inputs

   pure subroutine allocate_backtest(out, origins, horizon, state_count, coefficient_count, &
      parameter_count)
      !! Allocate and clear a rolling-backtest result.
      type(tsissm_backtest_t), intent(out) :: out !! Procedure result.
      integer, intent(in) :: origins !! Origins.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      integer, intent(in) :: state_count !! State count.
      integer, intent(in) :: coefficient_count !! Number of coefficient.
      integer, intent(in), optional :: parameter_count !! Number of parameter.
      integer :: stored_parameters

      allocate(out%forecast(origins, horizon), out%actual(origins, horizon))
      allocate(out%error(origins, horizon), out%lower(origins, horizon), out%upper(origins, horizon))
      allocate(out%crps(origins, horizon), out%log_score(origins, horizon))
      allocate(out%persistence(origins, state_count), out%coefficients(origins, coefficient_count))
      stored_parameters = 0
      if (present(parameter_count)) stored_parameters = parameter_count
      allocate(out%parameters(origins, stored_parameters))
      allocate(out%valid(origins, horizon), out%reestimated(origins))
      allocate(out%training_start(origins), out%training_end(origins))
      allocate(out%rmse(horizon), out%mean_absolute_error(horizon), out%bias(horizon))
      allocate(out%interval_coverage(horizon), out%interval_score(horizon))
      allocate(out%mean_crps(horizon), out%mean_log_score(horizon))
      out%forecast = 0.0_dp
      out%actual = 0.0_dp
      out%error = 0.0_dp
      out%lower = 0.0_dp
      out%upper = 0.0_dp
      out%crps = 0.0_dp
      out%log_score = 0.0_dp
      out%parameters = 0.0_dp
      out%valid = .false.
      out%reestimated = .false.
   end subroutine allocate_backtest

   pure subroutine store_backtest_origin(out, origin_index, actual, forecast)
      !! Store one origin and calculate its horizon-specific distribution scores.
      type(tsissm_backtest_t), intent(inout) :: out !! Procedure result, updated in place.
      integer, intent(in) :: origin_index !! Index of origin.
      real(dp), intent(in) :: actual(:) !! Observed values used for evaluation.
      type(tsissm_forecast_t), intent(in) :: forecast !! Forecast values.
      integer :: step

      do step = 1, size(actual)
         out%forecast(origin_index, step) = forecast%moments%mean(step)
         out%actual(origin_index, step) = actual(step)
         out%error(origin_index, step) = actual(step) - out%forecast(origin_index, step)
         out%lower(origin_index, step) = forecast%quantile(1, step)
         out%upper(origin_index, step) = forecast%quantile(2, step)
         call distribution_scores(actual(step), forecast%distribution(:, step), &
            out%crps(origin_index, step), out%log_score(origin_index, step))
         out%valid(origin_index, step) = .true.
      end do
   end subroutine store_backtest_origin

   pure subroutine distribution_scores(actual, draws, crps, log_score)
      !! Calculate one observation's ensemble CRPS and kernel log score.
      real(dp), intent(in) :: actual !! Observed values used for evaluation.
      real(dp), intent(in) :: draws(:) !! Draws.
      real(dp), intent(out) :: crps !! Crps.
      real(dp), intent(out) :: log_score !! Log score.
      real(dp), allocatable :: ordered(:)
      real(dp) :: bandwidth, density, difference, mean_value, standard_deviation
      integer :: i, simulations

      simulations = size(draws)
      ordered = draws
      call sort(ordered)
      crps = sum(abs(ordered - actual))/real(simulations, dp)
      do i = 1, simulations
         crps = crps - real(2*i - simulations - 1, dp)*ordered(i)/ &
            (real(simulations, dp)*real(simulations, dp))
      end do
      mean_value = sum(ordered)/real(simulations, dp)
      standard_deviation = sqrt(sum((ordered - mean_value)**2)/real(simulations - 1, dp))
      bandwidth = max(1.06_dp*standard_deviation*real(simulations, dp)**(-0.2_dp), &
         sqrt(epsilon(1.0_dp))*max(1.0_dp, abs(actual)))
      density = 0.0_dp
      do i = 1, simulations
         difference = (actual - ordered(i))/bandwidth
         density = density + exp(-0.5_dp*difference**2)
      end do
      density = density/(real(simulations, dp)*bandwidth*sqrt(2.0_dp*acos(-1.0_dp)))
      log_score = -log(max(density, tiny(1.0_dp)))
   end subroutine distribution_scores

   pure subroutine summarize_backtest(out, alpha)
      !! Aggregate rolling errors and distribution scores by forecast horizon.
      type(tsissm_backtest_t), intent(inout) :: out !! Procedure result, updated in place.
      real(dp), intent(in) :: alpha !! Significance, smoothing, or model coefficient.
      real(dp) :: score
      integer :: count_valid, origin, step

      do step = 1, size(out%forecast, 2)
         count_valid = count(out%valid(:, step))
         out%rmse(step) = sqrt(sum(pack(out%error(:, step)**2, out%valid(:, step)))/ &
            real(count_valid, dp))
         out%mean_absolute_error(step) = sum(pack(abs(out%error(:, step)), &
            out%valid(:, step)))/real(count_valid, dp)
         out%bias(step) = sum(pack(out%forecast(:, step) - out%actual(:, step), &
            out%valid(:, step)))/real(count_valid, dp)
         out%interval_coverage(step) = 0.0_dp
         out%interval_score(step) = 0.0_dp
         do origin = 1, size(out%forecast, 1)
            if (.not. out%valid(origin, step)) cycle
            if (out%actual(origin, step) >= out%lower(origin, step) .and. &
               out%actual(origin, step) <= out%upper(origin, step)) then
               out%interval_coverage(step) = out%interval_coverage(step) + 1.0_dp
            end if
            score = out%upper(origin, step) - out%lower(origin, step)
            if (out%actual(origin, step) < out%lower(origin, step)) then
               score = score + 2.0_dp*(out%lower(origin, step) - out%actual(origin, step))/alpha
            else if (out%actual(origin, step) > out%upper(origin, step)) then
               score = score + 2.0_dp*(out%actual(origin, step) - out%upper(origin, step))/alpha
            end if
            out%interval_score(step) = out%interval_score(step) + score
         end do
         out%interval_coverage(step) = out%interval_coverage(step)/real(count_valid, dp)
         out%interval_score(step) = out%interval_score(step)/real(count_valid, dp)
         out%mean_crps(step) = sum(pack(out%crps(:, step), out%valid(:, step)))/real(count_valid, dp)
         out%mean_log_score(step) = sum(pack(out%log_score(:, step), &
            out%valid(:, step)))/real(count_valid, dp)
      end do
   end subroutine summarize_backtest

   pure function tsissm_weighted_box_test(residuals, lag, test_type, &
      fitted_parameter_count, weighted, transform) result(out)
      !! Compute weighted or classical Box-Pierce, Ljung-Box, and Monti tests.
      real(dp), intent(in) :: residuals(:) !! Model residuals.
      integer, intent(in) :: lag !! Lag index or number of lags.
      integer, intent(in), optional :: test_type !! Test type.
      integer, intent(in), optional :: fitted_parameter_count !! Number of fitted parameter.
      integer, intent(in), optional :: transform !! Transform.
      logical, intent(in), optional :: weighted !! Flag controlling weighted.
      type(tsissm_weighted_box_test_t) :: out

      out = shared_weighted_box_test(residuals, lag, test_type, &
         fitted_parameter_count, weighted, transform)
   end function tsissm_weighted_box_test

   pure subroutine weighted_ljung_box(residuals, lag, fitted_parameter_count, statistic, &
      gamma_shape, gamma_scale)
      !! Compute Fisher-Gallagher's weighted Ljung-Box gamma approximation.
      real(dp), intent(in) :: residuals(:) !! Model residuals.
      integer, intent(in) :: lag !! Lag index or number of lags.
      integer, intent(in) :: fitted_parameter_count !! Number of fitted parameter.
      real(dp), intent(out) :: statistic !! Statistic.
      real(dp), intent(out) :: gamma_shape !! Gamma shape.
      real(dp), intent(out) :: gamma_scale !! Gamma scale.
      type(tsissm_weighted_box_test_t) :: test

      test = tsissm_weighted_box_test(residuals, lag, tsissm_box_test_ljung_box, &
         fitted_parameter_count, .true., tsissm_residual_raw)
      statistic = test%statistic
      gamma_shape = test%gamma_shape
      gamma_scale = test%gamma_scale
   end subroutine weighted_ljung_box

   function tsissm_simulate_constant(transition_base, transition_scale, transition_parameter, &
      observation_loading, persistence, regressors, coefficients, initial_state, sigma, lambda, &
      distribution, skew, shape, simulations, probabilities, seed, bootstrap_residuals) result(out)
      !! Simulate and summarize constant-variance parametric or bootstrap forecasts.
      real(dp), intent(in) :: transition_base(:, :) !! Transition base.
      real(dp), intent(in) :: transition_scale(:, :) !! Transition scale.
      real(dp), intent(in) :: transition_parameter(:, :) !! Transition parameter.
      real(dp), intent(in) :: observation_loading(:) !! Observation loading matrix.
      real(dp), intent(in) :: persistence(:) !! Persistence.
      real(dp), intent(in) :: regressors(:, :) !! Regression design matrix.
      real(dp), intent(in) :: coefficients(:) !! Model coefficients.
      real(dp), intent(in) :: initial_state(:) !! Initial state vector.
      real(dp), intent(in) :: sigma !! Scale parameter or standard deviation.
      real(dp), intent(in) :: lambda !! Penalty or shrinkage parameter.
      real(dp), intent(in) :: skew !! Skew.
      real(dp), intent(in) :: shape !! Shape.
      integer, intent(in) :: distribution !! Probability-distribution specification.
      integer, intent(in) :: simulations !! Number of simulation draws.
      real(dp), intent(in), optional :: probabilities(:) !! Probability values.
      real(dp), intent(in), optional :: bootstrap_residuals(:) !! Bootstrap residuals.
      integer, intent(in), optional :: seed !! Random-number seed.
      type(tsissm_forecast_t) :: out
      real(dp), allocatable :: innovations(:, :), selected_probabilities(:)

      if (simulations < 1 .or. size(regressors, 1) < 1 .or. sigma <= 0.0_dp .or. &
         .not. valid_distribution_parameters(distribution, shape)) then
         out%info = 1
         return
      end if
      call prepare_probabilities(probabilities, selected_probabilities, out%info)
      if (out%info /= 0) return
      if (present(seed)) call set_random_seed(seed)
      allocate(innovations(simulations, size(regressors, 1)))
      if (present(bootstrap_residuals)) then
         if (size(bootstrap_residuals) < 1) then
            out%info = 1
            return
         end if
         call resample_innovations(bootstrap_residuals, innovations)
      else
         call random_standardized_innovations(innovations, distribution, skew, shape)
         innovations = sigma*innovations
      end if
      out%prediction = tsissm_predict_constant(transition_base, transition_scale, &
         transition_parameter, observation_loading, persistence, regressors, coefficients, &
         initial_state, innovations)
      if (out%prediction%info /= 0) then
         out%info = 10 + out%prediction%info
         return
      end if
      out%moments = tsissm_moments_constant(transition_base*transition_scale*transition_parameter, &
         observation_loading, persistence, initial_state, regressors, coefficients, sigma**2, &
         lambda, .true.)
      call summarize_forecast(out, selected_probabilities, lambda)
   end function tsissm_simulate_constant

   function tsissm_simulate_dynamic(transition_base, transition_scale, transition_parameter, &
      observation_loading, persistence, regressors, coefficients, initial_state, arch, garch, &
      variance_intercept, initial_arch, initial_variance, lambda, distribution, skew, shape, &
      simulations, probabilities, seed, bootstrap_standardized_residuals) result(out)
      !! Simulate and summarize GARCH-scaled parametric or bootstrap forecasts.
      real(dp), intent(in) :: transition_base(:, :) !! Transition base.
      real(dp), intent(in) :: transition_scale(:, :) !! Transition scale.
      real(dp), intent(in) :: transition_parameter(:, :) !! Transition parameter.
      real(dp), intent(in) :: observation_loading(:) !! Observation loading matrix.
      real(dp), intent(in) :: persistence(:) !! Persistence.
      real(dp), intent(in) :: regressors(:, :) !! Regression design matrix.
      real(dp), intent(in) :: coefficients(:) !! Model coefficients.
      real(dp), intent(in) :: initial_state(:) !! Initial state vector.
      real(dp), intent(in) :: arch(:) !! Arch.
      real(dp), intent(in) :: garch(:) !! Garch.
      real(dp), intent(in) :: variance_intercept !! Variance intercept.
      real(dp), intent(in) :: initial_arch(:) !! Initial arch.
      real(dp), intent(in) :: initial_variance(:) !! Initial variance.
      real(dp), intent(in) :: lambda !! Penalty or shrinkage parameter.
      real(dp), intent(in) :: skew !! Skew.
      real(dp), intent(in) :: shape !! Shape.
      integer, intent(in) :: distribution !! Probability-distribution specification.
      integer, intent(in) :: simulations !! Number of simulation draws.
      real(dp), intent(in), optional :: probabilities(:) !! Probability values.
      real(dp), intent(in), optional :: bootstrap_standardized_residuals(:) !! Bootstrap standardized residuals.
      integer, intent(in), optional :: seed !! Random-number seed.
      type(tsissm_forecast_t) :: out
      real(dp), allocatable :: selected_probabilities(:), standardized(:, :)

      if (simulations < 1 .or. size(regressors, 1) < 1 .or. variance_intercept < 0.0_dp .or. &
         .not. valid_distribution_parameters(distribution, shape)) then
         out%info = 1
         return
      end if
      call prepare_probabilities(probabilities, selected_probabilities, out%info)
      if (out%info /= 0) return
      if (present(seed)) call set_random_seed(seed)
      allocate(standardized(simulations, size(regressors, 1)))
      if (present(bootstrap_standardized_residuals)) then
         if (size(bootstrap_standardized_residuals) < 1) then
            out%info = 1
            return
         end if
         call resample_innovations(bootstrap_standardized_residuals, standardized)
      else
         call random_standardized_innovations(standardized, distribution, skew, shape)
      end if
      out%prediction = tsissm_predict_dynamic(transition_base, transition_scale, &
         transition_parameter, observation_loading, persistence, regressors, coefficients, &
         initial_state, standardized, arch, garch, initial_arch, initial_variance, variance_intercept)
      if (out%prediction%info /= 0) then
         out%info = 10 + out%prediction%info
         return
      end if
      out%moments = tsissm_moments_dynamic(transition_base*transition_scale*transition_parameter, &
         observation_loading, persistence, initial_state, regressors, coefficients, arch, garch, &
         variance_intercept, initial_arch, initial_variance, lambda, .true.)
      call summarize_forecast(out, selected_probabilities, lambda)
   end function tsissm_simulate_dynamic

   function tsissm_profile_constant(transition_base, transition_scale, transition_parameter, &
      observation_loading, persistence, regressors, coefficients, initial_state, sigma, lambda, &
      distribution, skew, shape, horizon, profile_simulations, forecast_simulations, period, &
      seed_index, max_iterations, tolerance, seed) result(out)
      !! Simulate, refit, and evaluate constant-variance ISSM paths.
      real(dp), intent(in) :: transition_base(:, :) !! Transition base.
      real(dp), intent(in) :: transition_scale(:, :) !! Transition scale.
      real(dp), intent(in) :: transition_parameter(:, :) !! Transition parameter.
      real(dp), intent(in) :: observation_loading(:) !! Observation loading matrix.
      real(dp), intent(in) :: persistence(:) !! Persistence.
      real(dp), intent(in) :: regressors(:, :) !! Regression design matrix.
      real(dp), intent(in) :: coefficients(:) !! Model coefficients.
      real(dp), intent(in) :: initial_state(:) !! Initial state vector.
      real(dp), intent(in) :: sigma !! Scale parameter or standard deviation.
      real(dp), intent(in) :: lambda !! Penalty or shrinkage parameter.
      real(dp), intent(in) :: skew !! Skew.
      real(dp), intent(in) :: shape !! Shape.
      real(dp), intent(in) :: tolerance !! Numerical convergence tolerance.
      integer, intent(in) :: distribution !! Probability-distribution specification.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      integer, intent(in) :: profile_simulations !! Profile simulations.
      integer, intent(in) :: forecast_simulations !! Forecast simulations.
      integer, intent(in) :: period !! Seasonal period.
      integer, intent(in) :: seed_index(:) !! Index of seed.
      integer, intent(in) :: max_iterations !! Maximum number of algorithm iterations.
      integer, intent(in), optional :: seed !! Random-number seed.
      type(tsissm_profile_t) :: out
      type(tsissm_fit_t) :: fitted
      type(tsissm_forecast_t) :: generated, forecast
      real(dp), allocatable :: training(:)
      integer :: n, simulation

      n = size(regressors, 1) - horizon
      if (n <= period .or. horizon < 1 .or. profile_simulations < 1 .or. &
         forecast_simulations < 2 .or. period < 1 .or. sigma <= 0.0_dp .or. &
         max_iterations < 1 .or. tolerance <= 0.0_dp .or. &
         size(regressors, 2) /= size(coefficients) .or. &
         .not. valid_distribution_parameters(distribution, shape)) then
         out%info = 1
         return
      end if
      call allocate_profile(out, profile_simulations, horizon, &
         size(persistence) + size(coefficients))
      out%true_parameter = [persistence, coefficients]
      if (present(seed)) call set_random_seed(seed)
      do simulation = 1, profile_simulations
         generated = tsissm_simulate_constant(transition_base, transition_scale, &
            transition_parameter, observation_loading, persistence, regressors, coefficients, &
            initial_state, sigma, lambda, distribution, skew, shape, 1)
         if (generated%info /= 0) cycle
         training = generated%distribution(1, :n)
         if (any(training <= 0.0_dp)) cycle
         fitted = tsissm_fit_constant(training, transition_base, transition_scale, &
            transition_parameter, observation_loading, persistence, regressors(:n, :), &
            coefficients, initial_state, seed_index, lambda=lambda, distribution=distribution, &
            skew=skew, shape=shape, max_iterations=max_iterations, tolerance=tolerance)
         if (fitted%info < 0 .or. fitted%info > 4) cycle
         forecast = tsissm_simulate_constant(transition_base, transition_scale, &
            transition_parameter, observation_loading, fitted%persistence, regressors(n + 1:, :), &
            fitted%coefficients, fitted%filter%state(:, n), fitted%likelihood%sigma, lambda, &
            distribution, skew, shape, forecast_simulations)
         if (forecast%info /= 0) cycle
         out%parameter(simulation, :) = [fitted%persistence, fitted%coefficients]
         call store_profile_metrics(out, simulation, generated%distribution(1, n + 1:), &
            forecast, training, period)
      end do
      if (out%successful == 0) out%info = 2
   end function tsissm_profile_constant

   function tsissm_profile_dynamic(transition_base, transition_scale, transition_parameter, &
      observation_loading, persistence, regressors, coefficients, initial_state, arch, garch, &
      variance_intercept, initial_arch, initial_variance, lambda, distribution, skew, shape, &
      horizon, profile_simulations, forecast_simulations, period, seed_index, max_iterations, &
      tolerance, seed, variance_initialization, variance_sample_size) result(out)
      !! Simulate, refit, and evaluate variance-targeted GARCH ISSM paths.
      real(dp), intent(in) :: transition_base(:, :) !! Transition base.
      real(dp), intent(in) :: transition_scale(:, :) !! Transition scale.
      real(dp), intent(in) :: transition_parameter(:, :) !! Transition parameter.
      real(dp), intent(in) :: observation_loading(:) !! Observation loading matrix.
      real(dp), intent(in) :: persistence(:) !! Persistence.
      real(dp), intent(in) :: regressors(:, :) !! Regression design matrix.
      real(dp), intent(in) :: coefficients(:) !! Model coefficients.
      real(dp), intent(in) :: initial_state(:) !! Initial state vector.
      real(dp), intent(in) :: arch(:) !! Arch.
      real(dp), intent(in) :: garch(:) !! Garch.
      real(dp), intent(in) :: variance_intercept !! Variance intercept.
      real(dp), intent(in) :: initial_arch(:) !! Initial arch.
      real(dp), intent(in) :: initial_variance(:) !! Initial variance.
      real(dp), intent(in) :: lambda !! Penalty or shrinkage parameter.
      real(dp), intent(in) :: skew !! Skew.
      real(dp), intent(in) :: shape !! Shape.
      real(dp), intent(in) :: tolerance !! Numerical convergence tolerance.
      integer, intent(in) :: distribution !! Probability-distribution specification.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      integer, intent(in) :: profile_simulations !! Profile simulations.
      integer, intent(in) :: forecast_simulations !! Forecast simulations.
      integer, intent(in) :: period !! Seasonal period.
      integer, intent(in) :: seed_index(:) !! Index of seed.
      integer, intent(in) :: max_iterations !! Maximum number of algorithm iterations.
      integer, intent(in), optional :: seed !! Random-number seed.
      integer, intent(in), optional :: variance_initialization !! Variance initialization.
      integer, intent(in), optional :: variance_sample_size !! Variance sample size.
      type(tsissm_profile_t) :: out
      type(tsissm_fit_t) :: fitted
      type(tsissm_forecast_t) :: generated, forecast
      real(dp), allocatable :: forecast_arch(:), forecast_variance(:), training(:)
      integer :: history, n, simulation

      n = size(regressors, 1) - horizon
      history = max(size(arch), size(garch))
      if (n <= max(period, history) .or. horizon < 1 .or. profile_simulations < 1 .or. &
         forecast_simulations < 2 .or. period < 1 .or. history < 1 .or. &
         size(initial_arch) < history .or. size(initial_variance) < history .or. &
         variance_intercept < 0.0_dp .or. max_iterations < 1 .or. tolerance <= 0.0_dp .or. &
         size(regressors, 2) /= size(coefficients) .or. any(arch < 0.0_dp) .or. &
         any(garch < 0.0_dp) .or. .not. valid_distribution_parameters(distribution, shape)) then
         out%info = 1
         return
      end if
      call allocate_profile(out, profile_simulations, horizon, &
         size(persistence) + size(coefficients) + size(arch) + size(garch))
      out%true_parameter = [persistence, coefficients, arch, garch]
      allocate(forecast_arch(history), forecast_variance(history))
      if (present(seed)) call set_random_seed(seed)
      do simulation = 1, profile_simulations
         generated = tsissm_simulate_dynamic(transition_base, transition_scale, &
            transition_parameter, observation_loading, persistence, regressors, coefficients, &
            initial_state, arch, garch, variance_intercept, initial_arch, initial_variance, &
            lambda, distribution, skew, shape, 1)
         if (generated%info /= 0) cycle
         training = generated%distribution(1, :n)
         if (any(training <= 0.0_dp)) cycle
         fitted = tsissm_fit_dynamic(training, transition_base, transition_scale, &
            transition_parameter, observation_loading, persistence, regressors(:n, :), &
            coefficients, arch, garch, initial_state, seed_index, lambda=lambda, &
            distribution=distribution, skew=skew, shape=shape, max_iterations=max_iterations, &
            tolerance=tolerance, variance_initialization=variance_initialization, &
            variance_sample_size=variance_sample_size)
         if (fitted%info < 0 .or. fitted%info > 4) cycle
         forecast_arch = fitted%filter%innovation(n - history + 1:n)**2
         forecast_variance = fitted%filter%conditional_sd(n - history + 1:n)**2
         forecast = tsissm_simulate_dynamic(transition_base, transition_scale, &
            transition_parameter, observation_loading, fitted%persistence, regressors(n + 1:, :), &
            fitted%coefficients, fitted%filter%state(:, n), fitted%arch, fitted%garch, &
            fitted%variance_intercept, forecast_arch, forecast_variance, lambda, distribution, &
            skew, shape, forecast_simulations)
         if (forecast%info /= 0) cycle
         out%parameter(simulation, :) = [fitted%persistence, fitted%coefficients, &
            fitted%arch, fitted%garch]
         call store_profile_metrics(out, simulation, generated%distribution(1, n + 1:), &
            forecast, training, period)
      end do
      if (out%successful == 0) out%info = 2
   end function tsissm_profile_dynamic

   function tsissm_profile_structural_constant(regressors, slope, damped_slope, &
      seasonal_frequency, seasonal_harmonics, regular_seasonal, ar_order, ma_order, &
      initial_state, parameters, sigma, horizon, profile_simulations, forecast_simulations, &
      period, distribution, max_iterations, tolerance, seed) result(out)
      !! Profile complete constant-variance structural fits by simulation and refitting.
      real(dp), intent(in) :: regressors(:, :) !! Regression design matrix.
      real(dp), intent(in) :: seasonal_frequency(:) !! Seasonal frequency.
      real(dp), intent(in) :: initial_state(:) !! Initial state vector.
      real(dp), intent(in) :: parameters(:) !! Model parameter values.
      real(dp), intent(in) :: sigma !! Scale parameter or standard deviation.
      real(dp), intent(in) :: tolerance !! Numerical convergence tolerance.
      logical, intent(in) :: slope !! Flag controlling slope.
      logical, intent(in) :: damped_slope !! Flag controlling damped slope.
      logical, intent(in) :: regular_seasonal !! Flag controlling regular seasonal.
      integer, intent(in) :: seasonal_harmonics(:) !! Seasonal harmonics.
      integer, intent(in) :: ar_order !! Autoregressive order.
      integer, intent(in) :: ma_order !! Moving-average order.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      integer, intent(in) :: profile_simulations !! Profile simulations.
      integer, intent(in) :: forecast_simulations !! Forecast simulations.
      integer, intent(in) :: period !! Seasonal period.
      integer, intent(in) :: distribution !! Probability-distribution specification.
      integer, intent(in) :: max_iterations !! Maximum number of algorithm iterations.
      integer, intent(in), optional :: seed !! Random-number seed.
      type(tsissm_profile_t) :: out
      type(tsissm_forecast_t) :: forecast, generated
      type(tsissm_model_t) :: model
      type(tsissm_structural_fit_t) :: fitted
      real(dp), allocatable :: coefficients(:), training(:)
      real(dp) :: lambda, shape, skew
      integer :: n, simulation, status

      n = size(regressors, 1) - horizon
      call decode_structural_parameters(slope, damped_slope, seasonal_frequency, &
         seasonal_harmonics, regular_seasonal, ar_order, ma_order, size(regressors, 2), &
         parameters, distribution, model, coefficients, lambda, skew, shape, status)
      if (status /= 0) then
         out%info = 1
         return
      end if
      if (n <= period .or. horizon < 1 .or. profile_simulations < 1 .or. &
         forecast_simulations < 2 .or. period < 1 .or. sigma <= 0.0_dp .or. &
         size(initial_state) /= size(model%initial_state) .or. max_iterations < 1 .or. &
         tolerance <= 0.0_dp) then
         out%info = 1
         return
      end if
      call allocate_profile(out, profile_simulations, horizon, size(parameters))
      out%true_parameter = parameters
      if (present(seed)) call set_random_seed(seed)
      do simulation = 1, profile_simulations
         generated = tsissm_simulate_constant(model%transition_base, model%transition_scale, &
            model%transition_parameter, model%observation_loading, model%persistence, regressors, &
            coefficients, initial_state, sigma, lambda, distribution, skew, shape, 1)
         if (generated%info /= 0) cycle
         training = generated%distribution(1, :n)
         if (any(training <= 0.0_dp)) cycle
         fitted = tsissm_fit_structural_constant(training, slope, damped_slope, &
            seasonal_frequency, seasonal_harmonics, regular_seasonal, ar_order, ma_order, &
            regressors(:n, :), initial_state, parameters, distribution=distribution, &
            max_iterations=max_iterations, tolerance=tolerance)
         if (fitted%info < 0 .or. fitted%info > 4) cycle
         call structural_coefficients(fitted%parameters, slope, damped_slope, &
            size(seasonal_frequency), ar_order, ma_order, size(regressors, 2), coefficients)
         forecast = tsissm_simulate_constant(fitted%model%transition_base, &
            fitted%model%transition_scale, fitted%model%transition_parameter, &
            fitted%model%observation_loading, fitted%model%persistence, regressors(n + 1:, :), &
            coefficients, fitted%filter%state(:, n), fitted%likelihood%sigma, fitted%lambda, &
            distribution, fitted%skew, fitted%shape, forecast_simulations)
         if (forecast%info /= 0) cycle
         out%parameter(simulation, :) = fitted%parameters
         call store_profile_metrics(out, simulation, generated%distribution(1, n + 1:), &
            forecast, training, period)
      end do
      if (out%successful == 0) out%info = 2
   end function tsissm_profile_structural_constant

   function tsissm_profile_structural_dynamic(regressors, slope, damped_slope, &
      seasonal_frequency, seasonal_harmonics, regular_seasonal, ar_order, ma_order, &
      initial_state, parameters, arch, garch, variance_intercept, initial_arch, &
      initial_variance, horizon, profile_simulations, forecast_simulations, period, &
      distribution, max_iterations, tolerance, seed, variance_initialization, &
      variance_sample_size) result(out)
      !! Profile complete structural and GARCH fits by simulation and refitting.
      real(dp), intent(in) :: regressors(:, :) !! Regression design matrix.
      real(dp), intent(in) :: seasonal_frequency(:) !! Seasonal frequency.
      real(dp), intent(in) :: initial_state(:) !! Initial state vector.
      real(dp), intent(in) :: parameters(:) !! Model parameter values.
      real(dp), intent(in) :: arch(:) !! Arch.
      real(dp), intent(in) :: garch(:) !! Garch.
      real(dp), intent(in) :: variance_intercept !! Variance intercept.
      real(dp), intent(in) :: initial_arch(:) !! Initial arch.
      real(dp), intent(in) :: initial_variance(:) !! Initial variance.
      real(dp), intent(in) :: tolerance !! Numerical convergence tolerance.
      logical, intent(in) :: slope !! Flag controlling slope.
      logical, intent(in) :: damped_slope !! Flag controlling damped slope.
      logical, intent(in) :: regular_seasonal !! Flag controlling regular seasonal.
      integer, intent(in) :: seasonal_harmonics(:) !! Seasonal harmonics.
      integer, intent(in) :: ar_order !! Autoregressive order.
      integer, intent(in) :: ma_order !! Moving-average order.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      integer, intent(in) :: profile_simulations !! Profile simulations.
      integer, intent(in) :: forecast_simulations !! Forecast simulations.
      integer, intent(in) :: period !! Seasonal period.
      integer, intent(in) :: distribution !! Probability-distribution specification.
      integer, intent(in) :: max_iterations !! Maximum number of algorithm iterations.
      integer, intent(in), optional :: seed !! Random-number seed.
      integer, intent(in), optional :: variance_initialization !! Variance initialization.
      integer, intent(in), optional :: variance_sample_size !! Variance sample size.
      type(tsissm_profile_t) :: out
      type(tsissm_forecast_t) :: forecast, generated
      type(tsissm_model_t) :: model
      type(tsissm_structural_fit_t) :: fitted
      real(dp), allocatable :: coefficients(:), forecast_arch(:), forecast_variance(:), training(:)
      real(dp) :: lambda, shape, skew
      integer :: history, n, simulation, status

      n = size(regressors, 1) - horizon
      history = max(size(arch), size(garch))
      call decode_structural_parameters(slope, damped_slope, seasonal_frequency, &
         seasonal_harmonics, regular_seasonal, ar_order, ma_order, size(regressors, 2), &
         parameters, distribution, model, coefficients, lambda, skew, shape, status)
      if (status /= 0) then
         out%info = 1
         return
      end if
      if (n <= max(period, history) .or. horizon < 1 .or. &
         profile_simulations < 1 .or. forecast_simulations < 2 .or. period < 1 .or. &
         history < 1 .or. size(initial_arch) < history .or. size(initial_variance) < history .or. &
         variance_intercept < 0.0_dp .or. any(arch < 0.0_dp) .or. any(garch < 0.0_dp) .or. &
         size(initial_state) /= size(model%initial_state) .or. max_iterations < 1 .or. &
         tolerance <= 0.0_dp) then
         out%info = 1
         return
      end if
      call allocate_profile(out, profile_simulations, horizon, &
         size(parameters) + size(arch) + size(garch))
      out%true_parameter = [parameters, arch, garch]
      allocate(forecast_arch(history), forecast_variance(history))
      if (present(seed)) call set_random_seed(seed)
      do simulation = 1, profile_simulations
         generated = tsissm_simulate_dynamic(model%transition_base, model%transition_scale, &
            model%transition_parameter, model%observation_loading, model%persistence, regressors, &
            coefficients, initial_state, arch, garch, variance_intercept, initial_arch, &
            initial_variance, lambda, distribution, skew, shape, 1)
         if (generated%info /= 0) cycle
         training = generated%distribution(1, :n)
         if (any(training <= 0.0_dp)) cycle
         fitted = tsissm_fit_structural_dynamic(training, slope, damped_slope, seasonal_frequency, &
            seasonal_harmonics, regular_seasonal, ar_order, ma_order, regressors(:n, :), &
            initial_state, parameters, arch, garch, distribution=distribution, &
            max_iterations=max_iterations, tolerance=tolerance, &
            variance_initialization=variance_initialization, &
            variance_sample_size=variance_sample_size)
         if (fitted%info < 0 .or. fitted%info > 4) cycle
         call structural_coefficients(fitted%parameters(:size(parameters)), slope, damped_slope, &
            size(seasonal_frequency), ar_order, ma_order, size(regressors, 2), coefficients)
         forecast_arch = fitted%filter%innovation(n - history + 1:n)**2
         forecast_variance = fitted%filter%conditional_sd(n - history + 1:n)**2
         forecast = tsissm_simulate_dynamic(fitted%model%transition_base, &
            fitted%model%transition_scale, fitted%model%transition_parameter, &
            fitted%model%observation_loading, fitted%model%persistence, regressors(n + 1:, :), &
            coefficients, fitted%filter%state(:, n), fitted%arch, fitted%garch, &
            fitted%variance_intercept, forecast_arch, forecast_variance, fitted%lambda, &
            distribution, fitted%skew, fitted%shape, forecast_simulations)
         if (forecast%info /= 0) cycle
         out%parameter(simulation, :) = fitted%parameters
         call store_profile_metrics(out, simulation, generated%distribution(1, n + 1:), &
            forecast, training, period)
      end do
      if (out%successful == 0) out%info = 2
   end function tsissm_profile_structural_dynamic

   pure subroutine allocate_profile(out, simulations, horizon, parameters)
      !! Allocate and initialize a simulation profile result.
      type(tsissm_profile_t), intent(out) :: out !! Procedure result.
      integer, intent(in) :: simulations !! Number of simulation draws.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      integer, intent(in) :: parameters !! Model parameter values.

      allocate(out%parameter(simulations, parameters), out%true_parameter(parameters))
      allocate(out%actual(simulations, horizon), out%predicted(simulations, horizon))
      allocate(out%mape(simulations, horizon), out%percent_bias(simulations, horizon))
      allocate(out%mase(simulations, horizon), out%crps(simulations, horizon))
      allocate(out%valid(simulations))
      out%parameter = 0.0_dp
      out%true_parameter = 0.0_dp
      out%actual = 0.0_dp
      out%predicted = 0.0_dp
      out%mape = 0.0_dp
      out%percent_bias = 0.0_dp
      out%mase = 0.0_dp
      out%crps = 0.0_dp
      out%valid = .false.
      out%simulations = simulations
      out%horizon = horizon
   end subroutine allocate_profile

   pure subroutine store_profile_metrics(out, simulation, actual, forecast, training, period)
      !! Store point metrics and cumulative CRPS for one successful profile path.
      type(tsissm_profile_t), intent(inout) :: out !! Procedure result, updated in place.
      integer, intent(in) :: simulation !! Simulation.
      integer, intent(in) :: period !! Seasonal period.
      real(dp), intent(in) :: actual(:) !! Observed values used for evaluation.
      real(dp), intent(in) :: training(:) !! Training observations.
      type(tsissm_forecast_t), intent(in) :: forecast !! Forecast values.
      real(dp) :: crps_value, log_score, scale
      integer :: i

      scale = sum(abs(training(period + 1:) - training(:size(training) - period)))/ &
         real(size(training) - period, dp)
      out%actual(simulation, :) = actual
      out%predicted(simulation, :) = forecast%mean
      out%mape(simulation, :) = 100.0_dp*abs((forecast%mean - actual)/actual)
      out%percent_bias(simulation, :) = 100.0_dp*(forecast%mean - actual)/actual
      if (scale > tiny(1.0_dp)) then
         out%mase(simulation, :) = abs(forecast%mean - actual)/scale
      else
         out%mase(simulation, :) = huge(1.0_dp)
      end if
      do i = 1, size(actual)
         call distribution_scores(actual(i), forecast%distribution(:, i), crps_value, log_score)
         out%crps(simulation, i) = crps_value
         if (i > 1) then
            out%crps(simulation, i) = (real(i - 1, dp)*out%crps(simulation, i - 1) + &
               crps_value)/real(i, dp)
         end if
      end do
      out%valid(simulation) = .true.
      out%successful = out%successful + 1
   end subroutine store_profile_metrics

   subroutine random_standardized_innovations(innovations, distribution, skew, shape)
      !! Fill a matrix from a standardized tsissm innovation distribution.
      real(dp), intent(out) :: innovations(:, :) !! Model innovations.
      integer, intent(in) :: distribution !! Probability-distribution specification.
      real(dp), intent(in) :: skew !! Skew.
      real(dp), intent(in) :: shape !! Shape.
      integer :: i, j

      do j = 1, size(innovations, 2)
         do i = 1, size(innovations, 1)
            select case (distribution)
            case (tsissm_distribution_gaussian)
               innovations(i, j) = random_standard_normal()
            case (tsissm_distribution_student)
               innovations(i, j) = random_standard_student(shape)
            case (tsissm_distribution_johnson_su)
               innovations(i, j) = random_standard_johnson_su(skew, shape)
            end select
         end do
      end do
   end subroutine random_standardized_innovations

   subroutine resample_innovations(residuals, innovations)
      !! Sample empirical residuals independently with replacement.
      real(dp), intent(in) :: residuals(:) !! Model residuals.
      real(dp), intent(out) :: innovations(:, :) !! Model innovations.
      integer :: i, index, j

      do j = 1, size(innovations, 2)
         do i = 1, size(innovations, 1)
            index = min(size(residuals), 1 + int(random_uniform()*real(size(residuals), dp)))
            innovations(i, j) = residuals(index)
         end do
      end do
   end subroutine resample_innovations

   pure subroutine prepare_probabilities(probabilities, selected, info)
      !! Validate requested probabilities or supply the usual 95 percent interval.
      real(dp), intent(in), optional :: probabilities(:) !! Probability values.
      real(dp), allocatable, intent(out) :: selected(:) !! Selected.
      integer, intent(out) :: info !! Status code; zero indicates success.

      info = 0
      if (present(probabilities)) then
         if (size(probabilities) < 1 .or. any(probabilities < 0.0_dp) .or. &
            any(probabilities > 1.0_dp)) then
            info = 2
            return
         end if
         selected = probabilities
      else
         selected = [0.025_dp, 0.5_dp, 0.975_dp]
      end if
   end subroutine prepare_probabilities

   pure subroutine summarize_forecast(out, probabilities, lambda)
      !! Back-transform paths and calculate means and type-7 sample quantiles.
      type(tsissm_forecast_t), intent(inout) :: out !! Procedure result, updated in place.
      real(dp), intent(in) :: probabilities(:) !! Probability values.
      real(dp), intent(in) :: lambda !! Penalty or shrinkage parameter.
      real(dp), allocatable :: ordered(:)
      real(dp) :: base
      integer :: horizon, i, probability_index, simulation

      out%transformed = out%prediction%observation
      allocate(out%distribution(size(out%transformed, 1), size(out%transformed, 2)))
      do simulation = 1, size(out%transformed, 1)
         do horizon = 1, size(out%transformed, 2)
            base = 1.0_dp + lambda*out%transformed(simulation, horizon)
            if (abs(lambda) > sqrt(epsilon(1.0_dp)) .and. base <= 0.0_dp) then
               out%info = 3
               return
            end if
            out%distribution(simulation, horizon) = &
               inverse_box_cox_value(out%transformed(simulation, horizon), lambda)
         end do
      end do
      out%probabilities = probabilities
      allocate(out%mean(size(out%distribution, 2)))
      allocate(out%quantile(size(probabilities), size(out%distribution, 2)))
      out%mean = sum(out%distribution, dim=1)/real(size(out%distribution, 1), dp)
      allocate(ordered(size(out%distribution, 1)))
      do horizon = 1, size(out%distribution, 2)
         ordered = out%distribution(:, horizon)
         call sort(ordered)
         do probability_index = 1, size(probabilities)
            out%quantile(probability_index, horizon) = &
               quantile(ordered, probabilities(probability_index))
         end do
      end do
   end subroutine summarize_forecast

   pure logical function valid_distribution_parameters(distribution, shape) result(valid)
      !! Check the distribution code and its required shape parameter.
      integer, intent(in) :: distribution !! Probability-distribution specification.
      real(dp), intent(in) :: shape !! Shape.

      valid = distribution == tsissm_distribution_gaussian .or. &
         (distribution == tsissm_distribution_student .and. shape > 2.0_dp) .or. &
         (distribution == tsissm_distribution_johnson_su .and. shape > 0.0_dp)
   end function valid_distribution_parameters

   pure function tsissm_decompose_filter(model, filtered, regressors, coefficients, &
      original_observations, lambda, updated_states) result(out)
      !! Decompose filtered values into structural, regression, and irregular terms.
      type(tsissm_model_t), intent(in) :: model !! Model specification.
      type(tsissm_filter_t), intent(in) :: filtered !! Filtered.
      real(dp), intent(in) :: regressors(:, :) !! Regression design matrix.
      real(dp), intent(in) :: coefficients(:) !! Model coefficients.
      real(dp), intent(in) :: original_observations(:) !! Original observations.
      real(dp), intent(in) :: lambda !! Penalty or shrinkage parameter.
      logical, intent(in), optional :: updated_states !! Flag controlling updated states.
      type(tsissm_decomposition_t) :: out
      logical :: use_updated
      integer :: i, n, seasonal_count, state_column
      real(dp) :: base

      n = size(original_observations)
      seasonal_count = size(model%seasonal_start)
      use_updated = .false.
      if (present(updated_states)) use_updated = updated_states
      if (model%info /= 0 .or. filtered%info /= 0 .or. n < 1 .or. &
         size(filtered%state, 1) /= size(model%observation_loading) .or. &
         size(filtered%state, 2) /= n + 1 .or. size(filtered%innovation) /= n .or. &
         size(regressors, 1) /= n .or. size(regressors, 2) /= size(coefficients) .or. &
         size(model%seasonal_end) /= seasonal_count .or. any(original_observations <= 0.0_dp)) then
         out%info = 1
         return
      end if
      allocate(out%level(n), out%slope(n), out%seasonal(n, seasonal_count))
      allocate(out%autoregressive(n), out%moving_average(n), out%regression(n))
      allocate(out%irregular(n), out%sigma(n), out%trend(n), out%seasonal_total(n), out%arma(n))
      allocate(out%fitted_transformed(n), out%fitted_original(n))
      allocate(out%reconstructed_transformed(n), out%reconstructed_original(n))
      allocate(out%residual_original(n))
      out%slope = 0.0_dp
      out%seasonal = 0.0_dp
      out%autoregressive = 0.0_dp
      out%moving_average = 0.0_dp
      out%sigma = 0.0_dp
      if (allocated(filtered%conditional_sd)) out%sigma = filtered%conditional_sd
      do i = 1, n
         state_column = i - 1
         if (use_updated) state_column = i
         out%level(i) = model%observation_loading(model%level_index)* &
            filtered%state(model%level_index, state_column)
         if (model%slope_index > 0) then
            out%slope(i) = model%observation_loading(model%slope_index)* &
               filtered%state(model%slope_index, state_column)
         end if
         call decompose_seasonal_state(model, filtered%state(:, state_column), out%seasonal(i, :))
         if (model%ar_start > 0) then
            out%autoregressive(i) = dot_product(model%observation_loading(model%ar_start:model%ar_end), &
               filtered%state(model%ar_start:model%ar_end, state_column))
         end if
         if (model%ma_start > 0) then
            out%moving_average(i) = dot_product(model%observation_loading(model%ma_start:model%ma_end), &
               filtered%state(model%ma_start:model%ma_end, state_column))
         end if
         out%regression(i) = dot_product(regressors(i, :), coefficients)
      end do
      out%irregular = filtered%innovation
      out%trend = out%level + out%slope
      out%seasonal_total = sum(out%seasonal, dim=2)
      out%arma = out%autoregressive + out%moving_average
      out%fitted_transformed = out%trend + out%seasonal_total + out%arma + out%regression
      out%reconstructed_transformed = out%fitted_transformed + out%irregular
      do i = 1, n
         base = 1.0_dp + lambda*out%fitted_transformed(i)
         if (abs(lambda) > sqrt(epsilon(1.0_dp)) .and. base <= 0.0_dp) then
            out%info = 2
            return
         end if
         out%fitted_original(i) = inverse_box_cox_value(out%fitted_transformed(i), lambda)
         base = 1.0_dp + lambda*out%reconstructed_transformed(i)
         if (abs(lambda) > sqrt(epsilon(1.0_dp)) .and. base <= 0.0_dp) then
            out%info = 2
            return
         end if
         out%reconstructed_original(i) = inverse_box_cox_value(out%reconstructed_transformed(i), lambda)
      end do
      out%residual_original = original_observations - out%fitted_original
   end function tsissm_decompose_filter

   pure function tsissm_decompose_prediction(model, prediction, regressors, coefficients, &
      lambda, updated_states) result(out)
      !! Decompose simulated predictions using pre-update or updated state paths.
      type(tsissm_model_t), intent(in) :: model !! Model specification.
      type(tsissm_prediction_t), intent(in) :: prediction !! Prediction.
      real(dp), intent(in) :: regressors(:, :) !! Regression design matrix.
      real(dp), intent(in) :: coefficients(:) !! Model coefficients.
      real(dp), intent(in) :: lambda !! Penalty or shrinkage parameter.
      logical, intent(in), optional :: updated_states !! Flag controlling updated states.
      type(tsissm_prediction_decomposition_t) :: out
      logical :: use_updated
      integer :: h, i, seasonal_count, simulation, simulations, state_column
      real(dp) :: base

      simulations = size(prediction%observation, 1)
      h = size(prediction%observation, 2)
      seasonal_count = size(model%seasonal_start)
      use_updated = .false.
      if (present(updated_states)) use_updated = updated_states
      if (model%info /= 0 .or. prediction%info /= 0 .or. simulations < 1 .or. h < 1 .or. &
         size(prediction%state, 1) /= size(model%observation_loading) .or. &
         size(prediction%state, 2) /= h + 1 .or. size(prediction%state, 3) /= simulations .or. &
         any(shape(prediction%innovation) /= [simulations, h]) .or. &
         size(regressors, 1) /= h .or. size(regressors, 2) /= size(coefficients) .or. &
         size(model%seasonal_end) /= seasonal_count) then
         out%info = 1
         return
      end if
      allocate(out%level(simulations, h), out%slope(simulations, h))
      allocate(out%seasonal(simulations, h, seasonal_count))
      allocate(out%autoregressive(simulations, h), out%moving_average(simulations, h))
      allocate(out%regression(simulations, h), out%irregular(simulations, h))
      allocate(out%trend(simulations, h), out%seasonal_total(simulations, h), out%arma(simulations, h))
      allocate(out%reconstructed_transformed(simulations, h), out%reconstructed_original(simulations, h))
      out%slope = 0.0_dp
      out%seasonal = 0.0_dp
      out%autoregressive = 0.0_dp
      out%moving_average = 0.0_dp
      do simulation = 1, simulations
         do i = 1, h
            state_column = i - 1
            if (use_updated) state_column = i
            out%level(simulation, i) = model%observation_loading(model%level_index)* &
               prediction%state(model%level_index, state_column, simulation)
            if (model%slope_index > 0) then
               out%slope(simulation, i) = model%observation_loading(model%slope_index)* &
                  prediction%state(model%slope_index, state_column, simulation)
            end if
            call decompose_seasonal_state(model, prediction%state(:, state_column, simulation), &
               out%seasonal(simulation, i, :))
            if (model%ar_start > 0) then
               out%autoregressive(simulation, i) = dot_product( &
                  model%observation_loading(model%ar_start:model%ar_end), &
                  prediction%state(model%ar_start:model%ar_end, state_column, simulation))
            end if
            if (model%ma_start > 0) then
               out%moving_average(simulation, i) = dot_product( &
                  model%observation_loading(model%ma_start:model%ma_end), &
                  prediction%state(model%ma_start:model%ma_end, state_column, simulation))
            end if
            out%regression(simulation, i) = dot_product(regressors(i, :), coefficients)
         end do
      end do
      out%irregular = prediction%innovation
      out%trend = out%level + out%slope
      out%seasonal_total = sum(out%seasonal, dim=3)
      out%arma = out%autoregressive + out%moving_average
      out%reconstructed_transformed = out%trend + out%seasonal_total + out%arma + &
         out%regression + out%irregular
      do simulation = 1, simulations
         do i = 1, h
            base = 1.0_dp + lambda*out%reconstructed_transformed(simulation, i)
            if (abs(lambda) > sqrt(epsilon(1.0_dp)) .and. base <= 0.0_dp) then
               out%info = 2
               return
            end if
            out%reconstructed_original(simulation, i) = &
               inverse_box_cox_value(out%reconstructed_transformed(simulation, i), lambda)
         end do
      end do
   end function tsissm_decompose_prediction

   pure subroutine decompose_seasonal_state(model, state, contribution)
      !! Reduce each seasonal state block through its observation loading.
      type(tsissm_model_t), intent(in) :: model !! Model specification.
      real(dp), intent(in) :: state(:) !! State vector or state sequence.
      real(dp), intent(out) :: contribution(:) !! Contribution.
      integer :: i

      do i = 1, size(contribution)
         contribution(i) = dot_product( &
            model%observation_loading(model%seasonal_start(i):model%seasonal_end(i)), &
            state(model%seasonal_start(i):model%seasonal_end(i)))
      end do
   end subroutine decompose_seasonal_state

   pure elemental real(dp) function inverse_box_cox_value(value, lambda) result(inverse)
      !! Invert the Box-Cox transform for a valid transformed value.
      real(dp), intent(in) :: value !! Input value.
      real(dp), intent(in) :: lambda !! Penalty or shrinkage parameter.

      if (abs(lambda) <= sqrt(epsilon(1.0_dp))) then
         inverse = exp(value)
      else
         inverse = (1.0_dp + lambda*value)**(1.0_dp/lambda)
      end if
   end function inverse_box_cox_value

   pure real(dp) function tsissm_log_density(value, distribution, skew, shape) result(log_density)
      !! Evaluate a standardized Gaussian, Student-t, or Johnson SU log density.
      real(dp), intent(in) :: value !! Input value.
      real(dp), intent(in) :: skew !! Skew.
      real(dp), intent(in) :: shape !! Shape.
      integer, intent(in) :: distribution !! Probability-distribution specification.
      real(dp) :: c, omega, r, reciprocal_shape, scaled, w, z

      select case (distribution)
      case (tsissm_distribution_gaussian)
         log_density = -0.5_dp*(log(2.0_dp*acos(-1.0_dp)) + value**2)
      case (tsissm_distribution_student)
         if (shape <= 2.0_dp) then
            log_density = -huge(1.0_dp)
            return
         end if
         scaled = sqrt(shape/(shape - 2.0_dp))
         log_density = log_gamma(0.5_dp*(shape + 1.0_dp)) - log_gamma(0.5_dp*shape) - &
            0.5_dp*log(shape*acos(-1.0_dp)) - 0.5_dp*(shape + 1.0_dp)* &
            log(1.0_dp + (value*scaled)**2/shape) + log(scaled)
      case (tsissm_distribution_johnson_su)
         if (shape <= 0.0_dp) then
            log_density = -huge(1.0_dp)
            return
         end if
         reciprocal_shape = 1.0_dp/shape
         w = exp(reciprocal_shape**2)
         omega = -skew*reciprocal_shape
         c = sqrt(1.0_dp/(0.5_dp*(w - 1.0_dp)*(w*cosh(2.0_dp*omega) + 1.0_dp)))
         z = (value - c*sqrt(w)*sinh(omega))/c
         r = -skew + asinh(z)/reciprocal_shape
         log_density = -log(c) - log(reciprocal_shape) - 0.5_dp*log(z**2 + 1.0_dp) - &
            0.5_dp*log(2.0_dp*acos(-1.0_dp)) - 0.5_dp*r**2
      case default
         log_density = -huge(1.0_dp)
      end select
   end function tsissm_log_density

   pure function tsissm_likelihood(innovations, original_observations, lambda, distribution, &
      skew, shape, parameter_count, observed, conditional_sd) result(out)
      !! Evaluate constant- or dynamic-scale tsissm innovation likelihoods.
      real(dp), intent(in) :: innovations(:) !! Model innovations.
      real(dp), intent(in) :: original_observations(:) !! Original observations.
      real(dp), intent(in) :: lambda !! Penalty or shrinkage parameter.
      real(dp), intent(in) :: skew !! Skew.
      real(dp), intent(in) :: shape !! Shape.
      integer, intent(in) :: distribution !! Probability-distribution specification.
      integer, intent(in) :: parameter_count !! Number of parameter.
      logical, intent(in), optional :: observed(:) !! Flag controlling observed.
      real(dp), intent(in), optional :: conditional_sd(:) !! Conditional standard deviation.
      type(tsissm_likelihood_t) :: out
      logical, allocatable :: use_observation(:)
      real(dp) :: density
      integer :: i, n

      n = size(innovations)
      if (n < 1 .or. size(original_observations) /= n .or. any(original_observations <= 0.0_dp) .or. &
         parameter_count < 0 .or. &
         (distribution /= tsissm_distribution_gaussian .and. &
         distribution /= tsissm_distribution_student .and. &
         distribution /= tsissm_distribution_johnson_su) .or. &
         (distribution == tsissm_distribution_student .and. shape <= 2.0_dp) .or. &
         (distribution == tsissm_distribution_johnson_su .and. shape <= 0.0_dp)) then
         out%info = 1
         return
      end if
      if (present(observed)) then
         if (size(observed) /= n) then
            out%info = 1
            return
         end if
      end if
      if (present(conditional_sd)) then
         if (size(conditional_sd) /= n .or. any(conditional_sd <= 0.0_dp)) then
            out%info = 1
            return
         end if
      end if
      allocate(use_observation(n), out%contribution(n), out%standardized_residual(n), out%scale(n))
      use_observation = .true.
      if (present(observed)) use_observation = observed
      out%observations = count(use_observation)
      if (out%observations < 1) then
         out%info = 1
         return
      end if
      if (present(conditional_sd)) then
         out%scale = conditional_sd
         out%sigma = sqrt(sum(conditional_sd**2)/real(n, dp))
      else
         out%sigma = sqrt(sum(pack(innovations**2, use_observation))/real(n, dp))
         if (out%sigma <= tiny(1.0_dp)) then
            out%info = 2
            return
         end if
         out%scale = out%sigma
      end if
      out%contribution = 0.0_dp
      out%standardized_residual = 0.0_dp
      do i = 1, n
         if (.not. use_observation(i)) cycle
         out%standardized_residual(i) = innovations(i)/out%scale(i)
         density = tsissm_log_density(out%standardized_residual(i), distribution, skew, shape)
         out%contribution(i) = -density + log(out%scale(i)) - &
            (lambda - 1.0_dp)*log(original_observations(i))
      end do
      out%objective = real(out%observations, dp)/real(n, dp)*sum(out%contribution)
      out%package_log_likelihood = -out%objective
      out%log_likelihood = -sum(out%contribution)
      out%aic = -2.0_dp*out%log_likelihood + 2.0_dp*real(parameter_count, dp)
      out%bic = -2.0_dp*out%log_likelihood + log(real(out%observations, dp))* &
         real(parameter_count, dp)
   end function tsissm_likelihood

   pure real(dp) function tsissm_ar_constraint(coefficients, margin) result(residual)
      !! Return the AR companion stability margin.
      real(dp), intent(in) :: coefficients(:) !! Model coefficients.
      real(dp), intent(in) :: margin !! Margin.
      real(dp), allocatable :: companion(:, :)
      integer :: i

      if (size(coefficients) < 1 .or. margin < 0.0_dp .or. margin >= 1.0_dp) then
         residual = -huge(1.0_dp)
         return
      end if
      allocate(companion(size(coefficients), size(coefficients)))
      companion = 0.0_dp
      companion(1, :) = coefficients
      do i = 2, size(coefficients)
         companion(i, i - 1) = 1.0_dp
      end do
      residual = 1.0_dp - margin - spectral_radius(companion)
   end function tsissm_ar_constraint

   pure real(dp) function tsissm_ma_constraint(coefficients, margin) result(residual)
      !! Return the MA inverse-companion stability margin.
      real(dp), intent(in) :: coefficients(:) !! Model coefficients.
      real(dp), intent(in) :: margin !! Margin.
      real(dp), allocatable :: companion(:, :)
      integer :: i

      if (size(coefficients) < 1 .or. margin < 0.0_dp .or. margin >= 1.0_dp) then
         residual = -huge(1.0_dp)
         return
      end if
      allocate(companion(size(coefficients), size(coefficients)))
      companion = 0.0_dp
      companion(1, :) = -coefficients
      do i = 2, size(coefficients)
         companion(i, i - 1) = 1.0_dp
      end do
      residual = 1.0_dp - margin - spectral_radius(companion)
   end function tsissm_ma_constraint

   pure real(dp) function tsissm_stability_constraint(transition_base, transition_scale, &
      transition_parameter, observation_loading, persistence, margin) result(residual)
      !! Return the discount-transition stability margin.
      real(dp), intent(in) :: transition_base(:, :) !! Transition base.
      real(dp), intent(in) :: transition_scale(:, :) !! Transition scale.
      real(dp), intent(in) :: transition_parameter(:, :) !! Transition parameter.
      real(dp), intent(in) :: observation_loading(:) !! Observation loading matrix.
      real(dp), intent(in) :: persistence(:) !! Persistence.
      real(dp), intent(in) :: margin !! Margin.
      real(dp), allocatable :: discount(:, :)
      integer :: dimension

      dimension = size(observation_loading)
      if (dimension < 1 .or. any(shape(transition_base) /= [dimension, dimension]) .or. &
         any(shape(transition_scale) /= [dimension, dimension]) .or. &
         any(shape(transition_parameter) /= [dimension, dimension]) .or. &
         size(persistence) /= dimension .or. margin < 0.0_dp .or. margin >= 1.0_dp) then
         residual = -huge(1.0_dp)
         return
      end if
      discount = transition_base*transition_scale*transition_parameter - &
         outer_product(persistence, observation_loading)
      residual = 1.0_dp - margin - spectral_radius(discount)
   end function tsissm_stability_constraint

   pure function tsissm_parameter_constraints(transition_base, transition_scale, &
      transition_parameter, observation_loading, persistence, ar_coefficients, ma_coefficients, &
      level_persistence, slope_persistence, arch, garch, distribution, skew, shape, lambda, &
      lambda_lower, lambda_upper, margin) result(out)
      !! Collect ISSM, ARMA, slope, GARCH, distribution, and bound constraints.
      real(dp), intent(in) :: transition_base(:, :) !! Transition base.
      real(dp), intent(in) :: transition_scale(:, :) !! Transition scale.
      real(dp), intent(in) :: transition_parameter(:, :) !! Transition parameter.
      real(dp), intent(in) :: observation_loading(:) !! Observation loading matrix.
      real(dp), intent(in) :: persistence(:) !! Persistence.
      real(dp), intent(in) :: ar_coefficients(:) !! Autoregressive coefficients.
      real(dp), intent(in) :: ma_coefficients(:) !! Moving-average coefficients.
      real(dp), intent(in) :: level_persistence !! Level persistence.
      real(dp), intent(in) :: slope_persistence !! Slope persistence.
      real(dp), intent(in) :: arch(:) !! Arch.
      real(dp), intent(in) :: garch(:) !! Garch.
      real(dp), intent(in) :: skew !! Skew.
      real(dp), intent(in) :: shape !! Shape.
      real(dp), intent(in) :: lambda !! Penalty or shrinkage parameter.
      real(dp), intent(in) :: lambda_lower !! Lambda lower.
      real(dp), intent(in) :: lambda_upper !! Lambda upper.
      real(dp), intent(in) :: margin !! Margin.
      integer, intent(in) :: distribution !! Probability-distribution specification.
      type(tsissm_constraint_t) :: out
      integer :: i, position

      if (lambda_lower > lambda_upper .or. margin < 0.0_dp .or. margin >= 1.0_dp) then
         out%info = 1
         return
      end if
      allocate(out%residual(11 + 2*size(persistence)))
      out%residual = huge(1.0_dp)
      if (size(ar_coefficients) > 0) out%residual(1) = tsissm_ar_constraint(ar_coefficients, margin)
      if (size(ma_coefficients) > 0) out%residual(2) = tsissm_ma_constraint(ma_coefficients, margin)
      out%residual(3) = tsissm_stability_constraint(transition_base, transition_scale, &
         transition_parameter, observation_loading, persistence, margin)
      out%residual(4) = slope_persistence - level_persistence - 0.01_dp
      out%residual(5) = 0.99_dp - sum(arch) - sum(garch)
      if (size(arch) > 0) out%residual(6) = min(out%residual(6), minval(arch))
      if (size(garch) > 0) out%residual(6) = min(out%residual(6), minval(garch))
      select case (distribution)
      case (tsissm_distribution_gaussian)
         out%residual(7) = 1.0_dp
      case (tsissm_distribution_student)
         out%residual(7) = shape - 2.0_dp
      case (tsissm_distribution_johnson_su)
         out%residual(7) = shape
      case default
         out%residual(7) = -1.0_dp
      end select
      out%residual(8) = lambda - lambda_lower
      out%residual(9) = lambda_upper - lambda
      if (distribution == tsissm_distribution_johnson_su) then
         out%residual(10) = skew + 20.0_dp
         out%residual(11) = 20.0_dp - skew
      end if
      position = 12
      do i = 1, size(persistence)
         out%residual(position) = persistence(i)
         out%residual(position + 1) = 1.0_dp - persistence(i)
         position = position + 2
      end do
      out%feasible = all(out%residual >= 0.0_dp)
   end function tsissm_parameter_constraints

   pure function tsissm_initialize_states(observations, transition_base, transition_scale, &
      transition_parameter, observation_loading, persistence, regressors, coefficients, &
      initial_state, seed_index, observed, lambda) result(out)
      !! Estimate selected initial states by the tsissm propagated-loading least squares.
      real(dp), intent(in) :: observations(:) !! Observed time-series values.
      real(dp), intent(in) :: transition_base(:, :) !! Transition base.
      real(dp), intent(in) :: transition_scale(:, :) !! Transition scale.
      real(dp), intent(in) :: transition_parameter(:, :) !! Transition parameter.
      real(dp), intent(in) :: observation_loading(:) !! Observation loading matrix.
      real(dp), intent(in) :: persistence(:) !! Persistence.
      real(dp), intent(in) :: regressors(:, :) !! Regression design matrix.
      real(dp), intent(in) :: coefficients(:) !! Model coefficients.
      real(dp), intent(in) :: initial_state(:) !! Initial state vector.
      real(dp), intent(in) :: lambda !! Penalty or shrinkage parameter.
      integer, intent(in) :: seed_index(:) !! Index of seed.
      logical, intent(in), optional :: observed(:) !! Flag controlling observed.
      type(tsissm_initialization_t) :: out
      logical, allocatable :: use_observation(:)
      real(dp), allocatable :: design(:, :), inverse(:, :), state(:), target(:), transformed(:)
      real(dp), allocatable :: transition(:, :), discount(:, :), loading(:), normal(:, :), rhs(:)
      real(dp) :: error, previous_error, fitted
      integer :: i, info, n, seeds

      n = size(observations)
      seeds = size(seed_index)
      if (n < 1 .or. size(initial_state) < 1 .or. any(observations <= 0.0_dp) .or. &
         any(seed_index < 1) .or. any(seed_index > size(initial_state)) .or. &
         size(observation_loading) /= size(initial_state) .or. &
         size(persistence) /= size(initial_state) .or. size(regressors, 1) /= n .or. &
         size(regressors, 2) /= size(coefficients) .or. &
         any(shape(transition_base) /= [size(initial_state), size(initial_state)]) .or. &
         any(shape(transition_scale) /= [size(initial_state), size(initial_state)]) .or. &
         any(shape(transition_parameter) /= [size(initial_state), size(initial_state)])) then
         out%info = 1
         return
      end if
      if (present(observed)) then
         if (size(observed) /= n) then
            out%info = 1
            return
         end if
      end if
      allocate(out%state(size(initial_state)))
      out%state = initial_state
      if (seeds == 0) return
      allocate(use_observation(n), design(n, seeds), target(n), transformed(n))
      allocate(state(size(initial_state)), transition(size(initial_state), size(initial_state)))
      allocate(discount(size(initial_state), size(initial_state)), loading(size(initial_state)))
      use_observation = .true.
      if (present(observed)) use_observation = observed
      transition = transition_base*transition_scale*transition_parameter
      discount = transition - outer_product(persistence, observation_loading)
      transformed = box_cox(observations, lambda)
      state = initial_state
      loading = observation_loading
      previous_error = 0.0_dp
      do i = 1, n
         fitted = dot_product(observation_loading, state) + dot_product(regressors(i, :), coefficients)
         if (use_observation(i)) then
            error = transformed(i) - fitted
         else
            error = previous_error
         end if
         design(i, :) = loading(seed_index)
         target(i) = error
         state = matmul(transition, state) + persistence*error
         loading = matmul(loading, discount)
         previous_error = error
      end do
      normal = matmul(transpose(pack_rows(design, use_observation)), pack_rows(design, use_observation))
      rhs = matmul(transpose(pack_rows(design, use_observation)), pack(target, use_observation))
      call invert_matrix(normal, inverse, info)
      if (info /= 0) then
         out%info = 2
         return
      end if
      out%state(seed_index) = matmul(inverse, rhs)
   end function tsissm_initialize_states

   pure function tsissm_fit_constant(observations, transition_base, transition_scale, &
      transition_parameter, observation_loading, initial_persistence, regressors, &
      initial_coefficients, initial_state, seed_index, observed, lambda, distribution, &
      skew, shape, max_iterations, tolerance) result(out)
      !! Estimate ISSM persistence and regression coefficients by penalized BFGS.
      real(dp), intent(in) :: observations(:) !! Observed time-series values.
      real(dp), intent(in) :: transition_base(:, :) !! Transition base.
      real(dp), intent(in) :: transition_scale(:, :) !! Transition scale.
      real(dp), intent(in) :: transition_parameter(:, :) !! Transition parameter.
      real(dp), intent(in) :: observation_loading(:) !! Observation loading matrix.
      real(dp), intent(in) :: initial_persistence(:) !! Initial persistence.
      real(dp), intent(in) :: regressors(:, :) !! Regression design matrix.
      real(dp), intent(in) :: initial_coefficients(:) !! Initial coefficients.
      real(dp), intent(in) :: initial_state(:) !! Initial state vector.
      real(dp), intent(in) :: lambda !! Penalty or shrinkage parameter.
      real(dp), intent(in) :: skew !! Skew.
      real(dp), intent(in) :: shape !! Shape.
      real(dp), intent(in) :: tolerance !! Numerical convergence tolerance.
      integer, intent(in) :: seed_index(:) !! Index of seed.
      integer, intent(in) :: distribution !! Probability-distribution specification.
      integer, intent(in) :: max_iterations !! Maximum number of algorithm iterations.
      logical, intent(in), optional :: observed(:) !! Flag controlling observed.
      type(tsissm_fit_t) :: out
      type(optimization_result_t) :: optimized
      type(tsissm_initialization_t) :: initialized
      logical, allocatable :: use_observation(:)
      real(dp), allocatable :: initial(:)
      integer :: dimension

      dimension = size(initial_persistence)
      if (dimension < 1 .or. size(initial_coefficients) /= size(regressors, 2) .or. &
         max_iterations < 1 .or. tolerance <= 0.0_dp) then
         out%info = 1
         return
      end if
      allocate(use_observation(size(observations)))
      use_observation = .true.
      if (present(observed)) then
         if (size(observed) /= size(observations)) then
            out%info = 1
            return
         end if
         use_observation = observed
      end if
      if (count(use_observation) < 1) then
         out%info = 1
         return
      end if
      initialized = tsissm_initialize_states(observations, transition_base, transition_scale, &
         transition_parameter, observation_loading, initial_persistence, regressors, &
         initial_coefficients, initial_state, seed_index, use_observation, lambda)
      if (initialized%info /= 0) then
         out%info = 10 + initialized%info
         return
      end if
      initial = [initial_persistence, initial_coefficients]
      optimized = bfgs_minimize_fd(objective, initial, max_iterations, tolerance)
      out%iterations = optimized%iterations
      out%converged = optimized%converged
      out%objective = optimized%objective
      if (.not. allocated(optimized%parameters)) then
         out%info = 20 + optimized%info
         return
      end if
      out%persistence = optimized%parameters(:dimension)
      out%coefficients = optimized%parameters(dimension + 1:)
      out%filter = tsissm_filter_constant(observations, transition_base, transition_scale, &
         transition_parameter, observation_loading, out%persistence, regressors, out%coefficients, &
         initialized%state, use_observation, lambda)
      if (out%filter%info /= 0) then
         out%info = 30 + out%filter%info
         return
      end if
      out%likelihood = tsissm_likelihood(out%filter%innovation, observations, lambda, distribution, &
         skew, shape, size(optimized%parameters), use_observation)
      if (out%likelihood%info /= 0) then
         out%info = 40 + out%likelihood%info
         return
      end if
      out%info = optimized%info

   contains

      pure real(dp) function objective(parameters) result(value)
         !! Evaluate penalized likelihood for the current persistence and regressors.
         real(dp), intent(in) :: parameters(:) !! Model parameter values.
         type(tsissm_filter_t) :: trial_filter
         type(tsissm_likelihood_t) :: trial_likelihood
         real(dp) :: stability, penalty

         penalty = 1.0e6_dp*sum(max(0.0_dp, -parameters(:dimension))**2 + &
            max(0.0_dp, parameters(:dimension) - 1.0_dp)**2)
         stability = tsissm_stability_constraint(transition_base, transition_scale, &
            transition_parameter, observation_loading, parameters(:dimension), 0.001_dp)
         penalty = penalty + 1.0e6_dp*max(0.0_dp, -stability)**2
         trial_filter = tsissm_filter_constant(observations, transition_base, transition_scale, &
            transition_parameter, observation_loading, parameters(:dimension), regressors, &
            parameters(dimension + 1:), initialized%state, use_observation, lambda)
         if (trial_filter%info /= 0) then
            value = huge(1.0_dp)
            return
         end if
         trial_likelihood = tsissm_likelihood(trial_filter%innovation, observations, lambda, &
            distribution, skew, shape, size(parameters), use_observation)
         if (trial_likelihood%info /= 0) then
            value = huge(1.0_dp)
         else
            value = trial_likelihood%objective + penalty
         end if
      end function objective
   end function tsissm_fit_constant

   pure function tsissm_fit_dynamic(observations, transition_base, transition_scale, &
      transition_parameter, observation_loading, initial_persistence, regressors, &
      initial_coefficients, initial_arch, initial_garch, initial_state, seed_index, observed, &
      lambda, distribution, skew, shape, max_iterations, tolerance, variance_initialization, &
      variance_sample_size) result(out)
      !! Estimate ISSM, regression, and GARCH parameters by penalized BFGS.
      real(dp), intent(in) :: observations(:) !! Observed time-series values.
      real(dp), intent(in) :: transition_base(:, :) !! Transition base.
      real(dp), intent(in) :: transition_scale(:, :) !! Transition scale.
      real(dp), intent(in) :: transition_parameter(:, :) !! Transition parameter.
      real(dp), intent(in) :: observation_loading(:) !! Observation loading matrix.
      real(dp), intent(in) :: initial_persistence(:) !! Initial persistence.
      real(dp), intent(in) :: regressors(:, :) !! Regression design matrix.
      real(dp), intent(in) :: initial_coefficients(:) !! Initial coefficients.
      real(dp), intent(in) :: initial_arch(:) !! Initial arch.
      real(dp), intent(in) :: initial_garch(:) !! Initial GARCH.
      real(dp), intent(in) :: initial_state(:) !! Initial state vector.
      real(dp), intent(in) :: lambda !! Penalty or shrinkage parameter.
      real(dp), intent(in) :: skew !! Skew.
      real(dp), intent(in) :: shape !! Shape.
      real(dp), intent(in) :: tolerance !! Numerical convergence tolerance.
      integer, intent(in) :: seed_index(:) !! Index of seed.
      integer, intent(in) :: distribution !! Probability-distribution specification.
      integer, intent(in) :: max_iterations !! Maximum number of algorithm iterations.
      integer, intent(in), optional :: variance_initialization !! Variance initialization.
      integer, intent(in), optional :: variance_sample_size !! Variance sample size.
      logical, intent(in), optional :: observed(:) !! Flag controlling observed.
      type(tsissm_fit_t) :: out
      type(optimization_result_t) :: optimized
      type(tsissm_initialization_t) :: initialized
      logical, allocatable :: use_observation(:)
      real(dp), allocatable :: initial(:)
      integer :: arch_count, coefficient_count, dimension, garch_count, history

      dimension = size(initial_persistence)
      coefficient_count = size(initial_coefficients)
      arch_count = size(initial_arch)
      garch_count = size(initial_garch)
      history = max(arch_count, garch_count)
      if (dimension < 1 .or. coefficient_count /= size(regressors, 2) .or. history < 1 .or. &
         max_iterations < 1 .or. tolerance <= 0.0_dp .or. any(initial_arch < 0.0_dp) .or. &
         any(initial_garch < 0.0_dp)) then
         out%info = 1
         return
      end if
      allocate(use_observation(size(observations)))
      use_observation = .true.
      if (present(observed)) then
         if (size(observed) /= size(observations)) then
            out%info = 1
            return
         end if
         use_observation = observed
      end if
      if (count(use_observation) < 1) then
         out%info = 1
         return
      end if
      initialized = tsissm_initialize_states(observations, transition_base, transition_scale, &
         transition_parameter, observation_loading, initial_persistence, regressors, &
         initial_coefficients, initial_state, seed_index, use_observation, lambda)
      if (initialized%info /= 0) then
         out%info = 10 + initialized%info
         return
      end if
      initial = [initial_persistence, initial_coefficients, initial_arch, initial_garch]
      optimized = bfgs_minimize_fd(objective, initial, max_iterations, tolerance)
      out%iterations = optimized%iterations
      out%converged = optimized%converged
      out%objective = optimized%objective
      if (.not. allocated(optimized%parameters)) then
         out%info = 20 + optimized%info
         return
      end if
      out%persistence = optimized%parameters(:dimension)
      out%coefficients = optimized%parameters(dimension + 1:dimension + coefficient_count)
      out%arch = optimized%parameters(dimension + coefficient_count + 1: &
         dimension + coefficient_count + arch_count)
      out%garch = optimized%parameters(dimension + coefficient_count + arch_count + 1:)
      call variance_target(out%persistence, out%coefficients, out%arch, out%garch, &
         out%initial_variance, out%variance_intercept, out%filter, out%likelihood)
      if (out%filter%info /= 0) then
         out%info = 30 + out%filter%info
         return
      end if
      if (out%likelihood%info /= 0) then
         out%info = 40 + out%likelihood%info
         return
      end if
      out%info = optimized%info

   contains

      pure real(dp) function objective(parameters) result(value)
         !! Evaluate penalized dynamic-scale likelihood.
         real(dp), intent(in) :: parameters(:) !! Model parameter values.
         type(tsissm_filter_t) :: trial_filter
         type(tsissm_likelihood_t) :: trial_likelihood
         real(dp) :: initial_variance_value, intercept, penalty, stability
         integer :: arch_start, garch_start

         arch_start = dimension + coefficient_count + 1
         garch_start = arch_start + arch_count
         penalty = 1.0e6_dp*sum(max(0.0_dp, -parameters(:dimension))**2 + &
            max(0.0_dp, parameters(:dimension) - 1.0_dp)**2)
         penalty = penalty + 1.0e6_dp*sum(max(0.0_dp, -parameters(arch_start:))**2)
         penalty = penalty + 1.0e6_dp*max(0.0_dp, &
            sum(parameters(arch_start:garch_start - 1)) + sum(parameters(garch_start:)) - 0.99_dp)**2
         stability = tsissm_stability_constraint(transition_base, transition_scale, &
            transition_parameter, observation_loading, parameters(:dimension), 0.001_dp)
         penalty = penalty + 1.0e6_dp*max(0.0_dp, -stability)**2
         call variance_target(parameters(:dimension), &
            parameters(dimension + 1:dimension + coefficient_count), &
            parameters(arch_start:garch_start - 1), parameters(garch_start:), &
            initial_variance_value, intercept, trial_filter, trial_likelihood)
         if (trial_filter%info /= 0 .or. trial_likelihood%info /= 0) then
            value = huge(1.0_dp)
         else
            value = trial_likelihood%objective + penalty
         end if
      end function objective

      pure subroutine variance_target(persistence, coefficients, arch, garch, &
         initial_variance_value, intercept, filtered, likelihood)
         !! Filter once, target GARCH variance, then evaluate dynamic likelihood.
         real(dp), intent(in) :: persistence(:) !! Persistence.
         real(dp), intent(in) :: coefficients(:) !! Model coefficients.
         real(dp), intent(in) :: arch(:) !! Arch.
         real(dp), intent(in) :: garch(:) !! Garch.
         real(dp), intent(out) :: initial_variance_value !! Initial variance value.
         real(dp), intent(out) :: intercept !! Model intercept.
         type(tsissm_filter_t), intent(out) :: filtered !! Filtered.
         type(tsissm_likelihood_t), intent(out) :: likelihood !! Likelihood.
         type(tsissm_filter_t) :: constant_filter
         type(tsissm_variance_initialization_t) :: initialized_variance
         real(dp), allocatable :: initial_arch_history(:), initial_variance_history(:)

         constant_filter = tsissm_filter_constant(observations, transition_base, transition_scale, &
            transition_parameter, observation_loading, persistence, regressors, coefficients, &
            initialized%state, use_observation, lambda)
         if (constant_filter%info /= 0) then
            filtered%info = constant_filter%info
            likelihood%info = 1
            return
         end if
         initialized_variance = tsissm_initialize_variance(constant_filter%innovation, &
            use_observation, variance_initialization, variance_sample_size)
         if (initialized_variance%info /= 0) then
            filtered%info = 10 + initialized_variance%info
            likelihood%info = 1
            return
         end if
         initial_variance_value = initialized_variance%variance
         intercept = max(1.0e-12_dp, initial_variance_value*(1.0_dp - sum(arch) - sum(garch)))
         allocate(initial_arch_history(history), initial_variance_history(history))
         initial_arch_history = initial_variance_value
         initial_variance_history = initial_variance_value
         filtered = tsissm_filter_dynamic(observations, transition_base, transition_scale, &
            transition_parameter, observation_loading, persistence, regressors, coefficients, &
            initialized%state, use_observation, lambda, arch, garch, initial_arch_history, &
            initial_variance_history, intercept)
         if (filtered%info /= 0) then
            likelihood%info = 1
            return
         end if
         likelihood = tsissm_likelihood(filtered%innovation, observations, lambda, distribution, &
            skew, shape, dimension + coefficient_count + arch_count + garch_count, &
            use_observation, filtered%conditional_sd)
      end subroutine variance_target
   end function tsissm_fit_dynamic

   pure function tsissm_fit_structural_constant(observations, slope, damped_slope, &
      seasonal_frequency, seasonal_harmonics, regular_seasonal, ar_order, ma_order, &
      regressors, initial_state, initial_parameters, observed, distribution, &
      max_iterations, tolerance) result(out)
      !! Jointly estimate a constant-variance structural ISSM from a flat parameter vector.
      real(dp), intent(in) :: observations(:) !! Observed time-series values.
      real(dp), intent(in) :: seasonal_frequency(:) !! Seasonal frequency.
      real(dp), intent(in) :: regressors(:, :) !! Regression design matrix.
      integer, intent(in) :: seasonal_harmonics(:) !! Seasonal harmonics.
      integer, intent(in) :: ar_order !! Autoregressive order.
      integer, intent(in) :: ma_order !! Moving-average order.
      integer, intent(in) :: distribution !! Probability-distribution specification.
      integer, intent(in) :: max_iterations !! Maximum number of algorithm iterations.
      logical, intent(in) :: slope !! Flag controlling slope.
      logical, intent(in) :: damped_slope !! Flag controlling damped slope.
      logical, intent(in) :: regular_seasonal !! Flag controlling regular seasonal.
      real(dp), intent(in) :: initial_state(:) !! Initial state vector.
      real(dp), intent(in) :: initial_parameters(:) !! Initial parameter values.
      real(dp), intent(in) :: tolerance !! Numerical convergence tolerance.
      logical, intent(in), optional :: observed(:) !! Flag controlling observed.
      type(tsissm_structural_fit_t) :: out
      type(optimization_result_t) :: optimized
      logical, allocatable :: use_observation(:)
      integer :: distribution_count, expected

      distribution_count = 0
      if (distribution == tsissm_distribution_student) distribution_count = 1
      if (distribution == tsissm_distribution_johnson_su) distribution_count = 2
      expected = 1 + merge(1, 0, slope) + merge(1, 0, damped_slope) + &
         size(seasonal_frequency) + ar_order + ma_order + size(regressors, 2) + 1 + distribution_count
      if (size(initial_parameters) /= expected .or. size(regressors, 1) /= size(observations) .or. &
         ar_order < 0 .or. ma_order < 0 .or. max_iterations < 1 .or. tolerance <= 0.0_dp .or. &
         (distribution /= tsissm_distribution_gaussian .and. &
         distribution /= tsissm_distribution_student .and. &
         distribution /= tsissm_distribution_johnson_su)) then
         out%info = 1
         return
      end if
      allocate(use_observation(size(observations)))
      use_observation = .true.
      if (present(observed)) then
         if (size(observed) /= size(observations)) then
            out%info = 1
            return
         end if
         use_observation = observed
      end if
      optimized = bfgs_minimize_fd(objective, initial_parameters, max_iterations, tolerance)
      out%iterations = optimized%iterations
      out%converged = optimized%converged
      out%objective = optimized%objective
      if (.not. allocated(optimized%parameters)) then
         out%info = 20 + optimized%info
         return
      end if
      out%parameters = optimized%parameters
      call evaluate(out%parameters, out%model, out%filter, out%likelihood, out%lambda, &
         out%skew, out%shape, expected)
      if (expected /= 0) then
         out%info = 30 + expected
         return
      end if
      out%info = optimized%info

   contains

      pure real(dp) function objective(parameters) result(value)
         !! Evaluate joint structural likelihood with parameter-domain penalties.
         real(dp), intent(in) :: parameters(:) !! Model parameter values.
         type(tsissm_model_t) :: model
         type(tsissm_filter_t) :: filtered
         type(tsissm_likelihood_t) :: likelihood
         real(dp) :: lambda, skew, shape
         integer :: status

         call evaluate(parameters, model, filtered, likelihood, lambda, skew, shape, status)
         if (status /= 0) then
            value = huge(1.0_dp)
         else
            value = likelihood%objective + structural_penalty(parameters, model, lambda, skew, shape)
         end if
      end function objective

      pure subroutine evaluate(parameters, model, filtered, likelihood, lambda, skew, shape, &
         status)
         !! Decode the structural parameter order and run filter plus likelihood.
         real(dp), intent(in) :: parameters(:) !! Model parameter values.
         type(tsissm_model_t), intent(out) :: model !! Model specification.
         type(tsissm_filter_t), intent(out) :: filtered !! Filtered.
         type(tsissm_likelihood_t), intent(out) :: likelihood !! Likelihood.
         real(dp), intent(out) :: lambda !! Penalty or shrinkage parameter.
         real(dp), intent(out) :: skew !! Skew.
         real(dp), intent(out) :: shape !! Shape.
         integer, intent(out) :: status !! Status.
         real(dp), allocatable :: ar(:), coefficients(:), ma(:), seasonal_persistence(:)
         real(dp) :: alpha, beta, damping
         integer :: index

         index = 1
         alpha = parameters(index)
         index = index + 1
         beta = 0.0_dp
         if (slope) then
            beta = parameters(index)
            index = index + 1
         end if
         damping = 1.0_dp
         if (damped_slope) then
            damping = parameters(index)
            index = index + 1
         end if
         seasonal_persistence = parameters(index:index + size(seasonal_frequency) - 1)
         index = index + size(seasonal_frequency)
         ar = parameters(index:index + ar_order - 1)
         index = index + ar_order
         ma = parameters(index:index + ma_order - 1)
         index = index + ma_order
         coefficients = parameters(index:index + size(regressors, 2) - 1)
         index = index + size(regressors, 2)
         lambda = parameters(index)
         index = index + 1
         skew = 0.0_dp
         shape = 0.0_dp
         if (distribution == tsissm_distribution_student) then
            shape = parameters(index)
         else if (distribution == tsissm_distribution_johnson_su) then
            skew = parameters(index)
            shape = parameters(index + 1)
         end if
         model = tsissm_model(slope, damped_slope, seasonal_frequency, seasonal_harmonics, &
            regular_seasonal, ar, ma, alpha, beta, seasonal_persistence, damping)
         status = model%info
         if (status /= 0) return
         if (size(model%initial_state) /= size(initial_state)) then
            status = 1
            return
         end if
         filtered = tsissm_filter_constant(observations, model%transition_base, model%transition_scale, &
            model%transition_parameter, model%observation_loading, model%persistence, regressors, &
            coefficients, initial_state, use_observation, lambda)
         status = filtered%info
         if (status /= 0) return
         likelihood = tsissm_likelihood(filtered%innovation, observations, lambda, distribution, &
            skew, shape, size(parameters), use_observation)
         status = likelihood%info
      end subroutine evaluate

      pure real(dp) function structural_penalty(parameters, model, lambda, skew, shape) result(penalty)
         !! Penalize structural, transformation, and distribution constraint violations.
         real(dp), intent(in) :: parameters(:) !! Model parameter values.
         real(dp), intent(in) :: lambda !! Penalty or shrinkage parameter.
         real(dp), intent(in) :: skew !! Skew.
         real(dp), intent(in) :: shape !! Shape.
         type(tsissm_model_t), intent(in) :: model !! Model specification.
         real(dp) :: stability
         integer :: ar_start, ma_start

         penalty = 1.0e6_dp*(max(0.0_dp, -lambda)**2 + max(0.0_dp, lambda - 1.5_dp)**2)
         penalty = penalty + 1.0e6_dp*sum(max(0.0_dp, -model%persistence)**2 + &
            max(0.0_dp, model%persistence - 1.0_dp)**2)
         stability = tsissm_stability_constraint(model%transition_base, model%transition_scale, &
            model%transition_parameter, model%observation_loading, model%persistence, 0.001_dp)
         penalty = penalty + 1.0e6_dp*max(0.0_dp, -stability)**2
         ar_start = 2 + merge(1, 0, slope) + merge(1, 0, damped_slope) + size(seasonal_frequency)
         ma_start = ar_start + ar_order
         if (ar_order > 0) penalty = penalty + 1.0e6_dp*max(0.0_dp, &
            -tsissm_ar_constraint(parameters(ar_start:ar_start + ar_order - 1), 0.001_dp))**2
         if (ma_order > 0) penalty = penalty + 1.0e6_dp*max(0.0_dp, &
            -tsissm_ma_constraint(parameters(ma_start:ma_start + ma_order - 1), 0.001_dp))**2
         if (distribution == tsissm_distribution_student) then
            penalty = penalty + 1.0e6_dp*max(0.0_dp, 2.01_dp - shape)**2
         else if (distribution == tsissm_distribution_johnson_su) then
            penalty = penalty + 1.0e6_dp*(max(0.0_dp, -shape)**2 + max(0.0_dp, abs(skew) - 20.0_dp)**2)
         end if
      end function structural_penalty
   end function tsissm_fit_structural_constant

   pure function tsissm_fit_structural_dynamic(observations, slope, damped_slope, &
      seasonal_frequency, seasonal_harmonics, regular_seasonal, ar_order, ma_order, &
      regressors, initial_state, initial_parameters, initial_arch, initial_garch, &
      observed, distribution, max_iterations, tolerance, variance_initialization, &
      variance_sample_size) result(out)
      !! Jointly estimate structural, transformation, distribution, and GARCH parameters.
      real(dp), intent(in) :: observations(:) !! Observed time-series values.
      real(dp), intent(in) :: seasonal_frequency(:) !! Seasonal frequency.
      real(dp), intent(in) :: regressors(:, :) !! Regression design matrix.
      real(dp), intent(in) :: initial_state(:) !! Initial state vector.
      real(dp), intent(in) :: initial_parameters(:) !! Initial parameter values.
      real(dp), intent(in) :: initial_arch(:) !! Initial arch.
      real(dp), intent(in) :: initial_garch(:) !! Initial GARCH.
      real(dp), intent(in) :: tolerance !! Numerical convergence tolerance.
      logical, intent(in) :: slope !! Flag controlling slope.
      logical, intent(in) :: damped_slope !! Flag controlling damped slope.
      logical, intent(in) :: regular_seasonal !! Flag controlling regular seasonal.
      integer, intent(in) :: seasonal_harmonics(:) !! Seasonal harmonics.
      integer, intent(in) :: ar_order !! Autoregressive order.
      integer, intent(in) :: ma_order !! Moving-average order.
      integer, intent(in) :: distribution !! Probability-distribution specification.
      integer, intent(in) :: max_iterations !! Maximum number of algorithm iterations.
      logical, intent(in), optional :: observed(:) !! Flag controlling observed.
      integer, intent(in), optional :: variance_initialization !! Variance initialization.
      integer, intent(in), optional :: variance_sample_size !! Variance sample size.
      type(tsissm_structural_fit_t) :: out
      type(optimization_result_t) :: optimized
      logical, allocatable :: use_observation(:)
      real(dp), allocatable :: initial(:)
      integer :: distribution_count, status, structural_count

      distribution_count = 0
      if (distribution == tsissm_distribution_student) distribution_count = 1
      if (distribution == tsissm_distribution_johnson_su) distribution_count = 2
      structural_count = 1 + merge(1, 0, slope) + merge(1, 0, damped_slope) + &
         size(seasonal_frequency) + ar_order + ma_order + size(regressors, 2) + &
         1 + distribution_count
      if (size(initial_parameters) /= structural_count .or. &
         max(size(initial_arch), size(initial_garch)) < 1 .or. &
         size(regressors, 1) /= size(observations) .or. ar_order < 0 .or. ma_order < 0 .or. &
         max_iterations < 1 .or. tolerance <= 0.0_dp .or. any(initial_arch < 0.0_dp) .or. &
         any(initial_garch < 0.0_dp) .or. &
         (distribution /= tsissm_distribution_gaussian .and. &
         distribution /= tsissm_distribution_student .and. &
         distribution /= tsissm_distribution_johnson_su)) then
         out%info = 1
         return
      end if
      allocate(use_observation(size(observations)))
      use_observation = .true.
      if (present(observed)) then
         if (size(observed) /= size(observations)) then
            out%info = 1
            return
         end if
         use_observation = observed
      end if
      if (count(use_observation) < 1) then
         out%info = 1
         return
      end if
      initial = [initial_parameters, initial_arch, initial_garch]
      optimized = bfgs_minimize_fd(objective, initial, max_iterations, tolerance)
      out%iterations = optimized%iterations
      out%converged = optimized%converged
      out%objective = optimized%objective
      out%dynamic_variance = .true.
      if (.not. allocated(optimized%parameters)) then
         out%info = 20 + optimized%info
         return
      end if
      out%parameters = optimized%parameters
      out%arch = optimized%parameters(structural_count + 1: &
         structural_count + size(initial_arch))
      out%garch = optimized%parameters(structural_count + size(initial_arch) + 1:)
      call evaluate(out%parameters, out%model, out%filter, out%likelihood, out%lambda, &
         out%skew, out%shape, out%initial_variance, out%variance_intercept, status)
      if (status /= 0) then
         out%info = 30 + status
         return
      end if
      out%info = optimized%info

   contains

      pure real(dp) function objective(parameters) result(value)
         !! Evaluate penalized joint dynamic structural likelihood.
         real(dp), intent(in) :: parameters(:) !! Model parameter values.
         type(tsissm_model_t) :: model
         type(tsissm_filter_t) :: filtered
         type(tsissm_likelihood_t) :: likelihood
         real(dp) :: initial_variance, intercept, lambda, penalty, shape, skew
         integer :: evaluation_status

         call evaluate(parameters, model, filtered, likelihood, lambda, skew, shape, &
            initial_variance, intercept, evaluation_status)
         if (evaluation_status /= 0) then
            value = huge(1.0_dp)
            return
         end if
         penalty = joint_penalty(parameters, model, lambda, skew, shape)
         value = likelihood%objective + penalty
      end function objective

      pure subroutine evaluate(parameters, model, filtered, likelihood, lambda, skew, shape, &
         initial_variance, intercept, evaluation_status)
         !! Decode all parameters, target GARCH variance, and evaluate likelihood.
         real(dp), intent(in) :: parameters(:) !! Model parameter values.
         type(tsissm_model_t), intent(out) :: model !! Model specification.
         type(tsissm_filter_t), intent(out) :: filtered !! Filtered.
         type(tsissm_likelihood_t), intent(out) :: likelihood !! Likelihood.
         real(dp), intent(out) :: lambda !! Penalty or shrinkage parameter.
         real(dp), intent(out) :: skew !! Skew.
         real(dp), intent(out) :: shape !! Shape.
         real(dp), intent(out) :: initial_variance !! Initial variance.
         real(dp), intent(out) :: intercept !! Model intercept.
         integer, intent(out) :: evaluation_status !! Evaluation status.
         type(tsissm_filter_t) :: constant_filter
         type(tsissm_variance_initialization_t) :: initialized_variance
         real(dp), allocatable :: ar(:), arch(:), coefficients(:), garch(:), ma(:)
         real(dp), allocatable :: initial_arch_history(:), initial_variance_history(:)
         real(dp), allocatable :: seasonal_persistence(:)
         real(dp) :: alpha, beta, damping
         integer :: history, index

         index = 1
         alpha = parameters(index)
         index = index + 1
         beta = 0.0_dp
         if (slope) then
            beta = parameters(index)
            index = index + 1
         end if
         damping = 1.0_dp
         if (damped_slope) then
            damping = parameters(index)
            index = index + 1
         end if
         seasonal_persistence = parameters(index:index + size(seasonal_frequency) - 1)
         index = index + size(seasonal_frequency)
         ar = parameters(index:index + ar_order - 1)
         index = index + ar_order
         ma = parameters(index:index + ma_order - 1)
         index = index + ma_order
         coefficients = parameters(index:index + size(regressors, 2) - 1)
         index = index + size(regressors, 2)
         lambda = parameters(index)
         index = index + 1
         skew = 0.0_dp
         shape = 0.0_dp
         if (distribution == tsissm_distribution_student) then
            shape = parameters(index)
            index = index + 1
         else if (distribution == tsissm_distribution_johnson_su) then
            skew = parameters(index)
            shape = parameters(index + 1)
            index = index + 2
         end if
         arch = parameters(index:index + size(initial_arch) - 1)
         index = index + size(initial_arch)
         garch = parameters(index:)
         model = tsissm_model(slope, damped_slope, seasonal_frequency, seasonal_harmonics, &
            regular_seasonal, ar, ma, alpha, beta, seasonal_persistence, damping)
         evaluation_status = model%info
         if (evaluation_status /= 0) return
         if (size(model%initial_state) /= size(initial_state)) then
            evaluation_status = 1
            return
         end if
         constant_filter = tsissm_filter_constant(observations, model%transition_base, &
            model%transition_scale, model%transition_parameter, model%observation_loading, &
            model%persistence, regressors, coefficients, initial_state, use_observation, lambda)
         evaluation_status = constant_filter%info
         if (evaluation_status /= 0) return
         initialized_variance = tsissm_initialize_variance(constant_filter%innovation, &
            use_observation, variance_initialization, variance_sample_size)
         if (initialized_variance%info /= 0) then
            evaluation_status = 10 + initialized_variance%info
            return
         end if
         initial_variance = initialized_variance%variance
         intercept = max(1.0e-12_dp, initial_variance*(1.0_dp - sum(arch) - sum(garch)))
         history = max(size(arch), size(garch))
         allocate(initial_arch_history(history), initial_variance_history(history))
         initial_arch_history = initial_variance
         initial_variance_history = initial_variance
         filtered = tsissm_filter_dynamic(observations, model%transition_base, &
            model%transition_scale, model%transition_parameter, model%observation_loading, &
            model%persistence, regressors, coefficients, initial_state, use_observation, lambda, &
            arch, garch, initial_arch_history, initial_variance_history, intercept)
         evaluation_status = filtered%info
         if (evaluation_status /= 0) return
         likelihood = tsissm_likelihood(filtered%innovation, observations, lambda, distribution, &
            skew, shape, size(parameters), use_observation, filtered%conditional_sd)
         evaluation_status = likelihood%info
      end subroutine evaluate

      pure real(dp) function joint_penalty(parameters, model, lambda, skew, shape) result(penalty)
         !! Penalize joint structural, distribution, and GARCH constraint violations.
         real(dp), intent(in) :: parameters(:) !! Model parameter values.
         real(dp), intent(in) :: lambda !! Penalty or shrinkage parameter.
         real(dp), intent(in) :: skew !! Skew.
         real(dp), intent(in) :: shape !! Shape.
         type(tsissm_model_t), intent(in) :: model !! Model specification.
         real(dp) :: stability
         integer :: ar_start, arch_start, ma_start

         penalty = 1.0e6_dp*(max(0.0_dp, -lambda)**2 + max(0.0_dp, lambda - 1.5_dp)**2)
         penalty = penalty + 1.0e6_dp*sum(max(0.0_dp, -model%persistence)**2 + &
            max(0.0_dp, model%persistence - 1.0_dp)**2)
         stability = tsissm_stability_constraint(model%transition_base, model%transition_scale, &
            model%transition_parameter, model%observation_loading, model%persistence, 0.001_dp)
         penalty = penalty + 1.0e6_dp*max(0.0_dp, -stability)**2
         ar_start = 2 + merge(1, 0, slope) + merge(1, 0, damped_slope) + &
            size(seasonal_frequency)
         ma_start = ar_start + ar_order
         if (ar_order > 0) penalty = penalty + 1.0e6_dp*max(0.0_dp, &
            -tsissm_ar_constraint(parameters(ar_start:ar_start + ar_order - 1), 0.001_dp))**2
         if (ma_order > 0) penalty = penalty + 1.0e6_dp*max(0.0_dp, &
            -tsissm_ma_constraint(parameters(ma_start:ma_start + ma_order - 1), 0.001_dp))**2
         arch_start = structural_count + 1
         penalty = penalty + 1.0e6_dp*sum(max(0.0_dp, -parameters(arch_start:))**2)
         penalty = penalty + 1.0e6_dp*max(0.0_dp, sum(parameters(arch_start:)) - 0.99_dp)**2
         if (distribution == tsissm_distribution_student) then
            penalty = penalty + 1.0e6_dp*max(0.0_dp, 2.01_dp - shape)**2
         else if (distribution == tsissm_distribution_johnson_su) then
            penalty = penalty + 1.0e6_dp*(max(0.0_dp, -shape)**2 + &
               max(0.0_dp, abs(skew) - 20.0_dp)**2)
         end if
      end function joint_penalty
   end function tsissm_fit_structural_dynamic

   pure function tsissm_moments_constant(transition, observation_loading, persistence, &
      initial_state, regressors, coefficients, innovation_variance, lambda, transform) result(out)
      !! Compute constant-scale ISSM analytical forecast moments.
      real(dp), intent(in) :: transition(:, :) !! State transition matrix.
      real(dp), intent(in) :: observation_loading(:) !! Observation loading matrix.
      real(dp), intent(in) :: persistence(:) !! Persistence.
      real(dp), intent(in) :: initial_state(:) !! Initial state vector.
      real(dp), intent(in) :: regressors(:, :) !! Regression design matrix.
      real(dp), intent(in) :: coefficients(:) !! Model coefficients.
      real(dp), intent(in) :: innovation_variance !! Innovation variance.
      real(dp), intent(in) :: lambda !! Penalty or shrinkage parameter.
      logical, intent(in) :: transform !! Flag controlling transform.
      type(tsissm_moments_t) :: out
      real(dp) :: accumulated, impulse(size(initial_state)), state(size(initial_state))
      integer :: horizon, i

      horizon = size(regressors, 1)
      if (horizon < 1 .or. innovation_variance < 0.0_dp .or. &
         any(shape(transition) /= [size(initial_state), size(initial_state)]) .or. &
         size(observation_loading) /= size(initial_state) .or. &
         size(persistence) /= size(initial_state) .or. size(regressors, 2) /= size(coefficients)) then
         out%info = 1
         return
      end if
      allocate(out%transformed_mean(horizon), out%transformed_variance(horizon))
      allocate(out%mean(horizon), out%variance(horizon))
      state = initial_state
      impulse = persistence
      accumulated = 0.0_dp
      do i = 1, horizon
         out%transformed_mean(i) = dot_product(observation_loading, state) + &
            dot_product(regressors(i, :), coefficients)
         if (i > 1) then
            accumulated = accumulated + dot_product(observation_loading, impulse)**2
            impulse = matmul(transition, impulse)
         end if
         out%transformed_variance(i) = innovation_variance*(1.0_dp + accumulated)
         state = matmul(transition, state)
      end do
      out%mean = out%transformed_mean
      out%variance = out%transformed_variance
      if (transform) call box_cox_moments(out%transformed_mean, out%transformed_variance, &
         lambda, out%mean, out%variance, out%info)
   end function tsissm_moments_constant

   pure function tsissm_moments_dynamic(transition, observation_loading, persistence, &
      initial_state, regressors, coefficients, arch, garch, variance_intercept, &
      initial_arch, initial_variance, lambda, transform) result(out)
      !! Compute GARCH-scale ISSM analytical forecast moments.
      real(dp), intent(in) :: transition(:, :) !! State transition matrix.
      real(dp), intent(in) :: observation_loading(:) !! Observation loading matrix.
      real(dp), intent(in) :: persistence(:) !! Persistence.
      real(dp), intent(in) :: initial_state(:) !! Initial state vector.
      real(dp), intent(in) :: regressors(:, :) !! Regression design matrix.
      real(dp), intent(in) :: coefficients(:) !! Model coefficients.
      real(dp), intent(in) :: arch(:) !! Arch.
      real(dp), intent(in) :: garch(:) !! Garch.
      real(dp), intent(in) :: variance_intercept !! Variance intercept.
      real(dp), intent(in) :: initial_arch(:) !! Initial arch.
      real(dp), intent(in) :: initial_variance(:) !! Initial variance.
      real(dp), intent(in) :: lambda !! Penalty or shrinkage parameter.
      logical, intent(in) :: transform !! Flag controlling transform.
      type(tsissm_moments_t) :: out
      real(dp), allocatable :: expected_square(:), variance_history(:)
      real(dp) :: accumulated, impulse(size(initial_state)), state(size(initial_state))
      real(dp) :: unconditional
      integer :: history, horizon, i, j, lag

      horizon = size(regressors, 1)
      history = max(size(arch), size(garch))
      if (horizon < 1 .or. history < 1 .or. size(initial_arch) < history .or. &
         size(initial_variance) < history .or. variance_intercept < 0.0_dp .or. &
         any(shape(transition) /= [size(initial_state), size(initial_state)]) .or. &
         size(observation_loading) /= size(initial_state) .or. &
         size(persistence) /= size(initial_state) .or. size(regressors, 2) /= size(coefficients)) then
         out%info = 1
         return
      end if
      allocate(expected_square(history + horizon), variance_history(history + horizon))
      expected_square = 0.0_dp
      variance_history = 0.0_dp
      expected_square(1:history) = initial_arch(1:history)
      variance_history(1:history) = initial_variance(1:history)
      do i = 1, horizon
         variance_history(history + i) = variance_intercept
         do j = 1, size(arch)
            lag = history + i - j
            variance_history(history + i) = variance_history(history + i) + arch(j)*expected_square(lag)
         end do
         do j = 1, size(garch)
            lag = history + i - j
            variance_history(history + i) = variance_history(history + i) + garch(j)*variance_history(lag)
         end do
         expected_square(history + i) = variance_history(history + i)
      end do
      allocate(out%transformed_mean(horizon), out%transformed_variance(horizon))
      allocate(out%mean(horizon), out%variance(horizon), out%garch_variance(horizon))
      out%garch_variance = variance_history(history + 1:)
      state = initial_state
      impulse = persistence
      accumulated = 0.0_dp
      do i = 1, horizon
         out%transformed_mean(i) = dot_product(observation_loading, state) + &
            dot_product(regressors(i, :), coefficients)
         if (i > 1) then
            accumulated = accumulated + dot_product(observation_loading, impulse)**2
            impulse = matmul(transition, impulse)
         end if
         out%transformed_variance(i) = out%garch_variance(i)*(1.0_dp + accumulated)
         state = matmul(transition, state)
      end do
      out%mean = out%transformed_mean
      out%variance = out%transformed_variance
      if (transform) then
         unconditional = variance_intercept/max(1.0_dp - sum(arch) - sum(garch), tiny(1.0_dp))
         call box_cox_moments(out%transformed_mean, out%transformed_variance, lambda, &
            out%mean, out%variance, out%info, out%garch_variance, unconditional)
      end if
   end function tsissm_moments_dynamic

   pure function tsissm_hresiduals_constant(observations, filtered, regressors, coefficients, &
      innovation_variance, lambda, horizon, first_origin, transformed) result(out)
      !! Compute expanding-origin residuals for every requested forecast horizon.
      real(dp), intent(in) :: observations(:) !! Observed time-series values.
      real(dp), intent(in) :: regressors(:, :) !! Regression design matrix.
      real(dp), intent(in) :: coefficients(:) !! Model coefficients.
      type(tsissm_filter_t), intent(in) :: filtered !! Filtered.
      real(dp), intent(in) :: innovation_variance !! Innovation variance.
      real(dp), intent(in) :: lambda !! Penalty or shrinkage parameter.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      integer, intent(in) :: first_origin !! First origin.
      logical, intent(in) :: transformed !! Flag controlling transformed.
      type(tsissm_hresiduals_t) :: out
      type(tsissm_moments_t) :: moments
      integer :: available, current, i, origins

      if (.not. hresidual_inputs_valid(observations, filtered, regressors, coefficients, &
         horizon, first_origin) .or. innovation_variance < 0.0_dp) then
         out%info = 1
         return
      end if
      origins = size(observations) - first_origin
      call allocate_hresiduals(out, origins, min(horizon, origins), first_origin, transformed)
      do i = 1, origins
         current = first_origin + i - 1
         available = min(out%horizon, size(observations) - current)
         moments = tsissm_moments_constant(filtered%transition, filtered%observation_loading, &
            filtered%persistence, filtered%state(:, current), &
            regressors(current + 1:current + available, :), coefficients, &
            innovation_variance, lambda, .not. transformed)
         if (moments%info /= 0) then
            out%info = 2
            return
         end if
         out%forecast(i, 1:available) = moments%mean
         if (transformed) then
            out%actual(i, 1:available) = box_cox(observations(current + 1:current + available), lambda)
         else
            out%actual(i, 1:available) = observations(current + 1:current + available)
         end if
         out%residual(i, 1:available) = out%actual(i, 1:available) - out%forecast(i, 1:available)
         out%valid(i, 1:available) = .true.
      end do
      call summarize_hresiduals(out)
   end function tsissm_hresiduals_constant

   pure function tsissm_hresiduals_dynamic(observations, filtered, regressors, coefficients, &
      arch, garch, variance_intercept, initial_arch, initial_variance, lambda, horizon, &
      first_origin, transformed) result(out)
      !! Compute expanding-origin residuals with origin-specific GARCH histories.
      real(dp), intent(in) :: observations(:) !! Observed time-series values.
      real(dp), intent(in) :: regressors(:, :) !! Regression design matrix.
      real(dp), intent(in) :: coefficients(:) !! Model coefficients.
      type(tsissm_filter_t), intent(in) :: filtered !! Filtered.
      real(dp), intent(in) :: arch(:) !! Arch.
      real(dp), intent(in) :: garch(:) !! Garch.
      real(dp), intent(in) :: variance_intercept !! Variance intercept.
      real(dp), intent(in) :: initial_arch(:) !! Initial arch.
      real(dp), intent(in) :: initial_variance(:) !! Initial variance.
      real(dp), intent(in) :: lambda !! Penalty or shrinkage parameter.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      integer, intent(in) :: first_origin !! First origin.
      logical, intent(in) :: transformed !! Flag controlling transformed.
      type(tsissm_hresiduals_t) :: out
      type(tsissm_moments_t) :: moments
      real(dp), allocatable :: arch_history(:), variance_history(:)
      integer :: available, current, history, i, origins

      history = max(size(arch), size(garch))
      if (.not. hresidual_inputs_valid(observations, filtered, regressors, coefficients, &
         horizon, first_origin) .or. history < 1 .or. .not. allocated(filtered%conditional_sd) .or. &
         size(filtered%conditional_sd) /= size(observations) .or. size(initial_arch) < history .or. &
         size(initial_variance) < history .or. variance_intercept < 0.0_dp .or. &
         any(arch < 0.0_dp) .or. any(garch < 0.0_dp)) then
         out%info = 1
         return
      end if
      origins = size(observations) - first_origin
      allocate(arch_history(history), variance_history(history))
      call allocate_hresiduals(out, origins, min(horizon, origins), first_origin, transformed)
      do i = 1, origins
         current = first_origin + i - 1
         available = min(out%horizon, size(observations) - current)
         call garch_history_at_origin(filtered, initial_arch, initial_variance, current, &
            arch_history, variance_history)
         moments = tsissm_moments_dynamic(filtered%transition, filtered%observation_loading, &
            filtered%persistence, filtered%state(:, current), &
            regressors(current + 1:current + available, :), coefficients, arch, garch, &
            variance_intercept, arch_history, variance_history, lambda, .not. transformed)
         if (moments%info /= 0) then
            out%info = 2
            return
         end if
         out%forecast(i, 1:available) = moments%mean
         if (transformed) then
            out%actual(i, 1:available) = box_cox(observations(current + 1:current + available), lambda)
         else
            out%actual(i, 1:available) = observations(current + 1:current + available)
         end if
         out%residual(i, 1:available) = out%actual(i, 1:available) - out%forecast(i, 1:available)
         out%valid(i, 1:available) = .true.
      end do
      call summarize_hresiduals(out)
   end function tsissm_hresiduals_dynamic

   pure logical function hresidual_inputs_valid(observations, filtered, regressors, &
      coefficients, horizon, first_origin) result(valid)
      !! Validate dimensions shared by constant and dynamic horizon residuals.
      real(dp), intent(in) :: observations(:) !! Observed time-series values.
      real(dp), intent(in) :: regressors(:, :) !! Regression design matrix.
      real(dp), intent(in) :: coefficients(:) !! Model coefficients.
      type(tsissm_filter_t), intent(in) :: filtered !! Filtered.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      integer, intent(in) :: first_origin !! First origin.

      valid = .false.
      if (.not. allocated(filtered%state) .or. .not. allocated(filtered%transition) .or. &
         .not. allocated(filtered%observation_loading) .or. .not. allocated(filtered%persistence)) return
      valid = size(observations) > 0 .and. all(observations > 0.0_dp) .and. horizon > 0 .and. &
         first_origin >= 0 .and. first_origin < size(observations) .and. &
         size(regressors, 1) == size(observations) .and. &
         size(regressors, 2) == size(coefficients) .and. &
         lbound(filtered%state, 2) == 0 .and. ubound(filtered%state, 2) >= size(observations) .and. &
         size(filtered%state, 1) == size(filtered%observation_loading) .and. &
         size(filtered%persistence) == size(filtered%observation_loading) .and. &
         all(shape(filtered%transition) == [size(filtered%state, 1), size(filtered%state, 1)])
   end function hresidual_inputs_valid

   pure subroutine allocate_hresiduals(out, origins, horizon, first_origin, transformed)
      !! Allocate and initialize a multi-horizon residual result.
      type(tsissm_hresiduals_t), intent(out) :: out !! Procedure result.
      integer, intent(in) :: origins !! Origins.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      integer, intent(in) :: first_origin !! First origin.
      logical, intent(in) :: transformed !! Flag controlling transformed.

      allocate(out%forecast(origins, horizon), out%actual(origins, horizon))
      allocate(out%residual(origins, horizon), out%valid(origins, horizon))
      allocate(out%covariance(horizon, horizon), out%pair_count(horizon, horizon))
      allocate(out%mean_residual(horizon), out%rmse(horizon))
      out%forecast = 0.0_dp
      out%actual = 0.0_dp
      out%residual = 0.0_dp
      out%valid = .false.
      out%covariance = 0.0_dp
      out%pair_count = 0
      out%mean_residual = 0.0_dp
      out%rmse = 0.0_dp
      out%first_origin = first_origin
      out%horizon = horizon
      out%transformed = transformed
   end subroutine allocate_hresiduals

   pure subroutine garch_history_at_origin(filtered, initial_arch, initial_variance, &
      origin, arch_history, variance_history)
      !! Form the latest squared-innovation and variance histories at an origin.
      type(tsissm_filter_t), intent(in) :: filtered !! Filtered.
      real(dp), intent(in) :: initial_arch(:) !! Initial arch.
      real(dp), intent(in) :: initial_variance(:) !! Initial variance.
      integer, intent(in) :: origin !! Origin.
      real(dp), intent(out) :: arch_history(:) !! Arch history.
      real(dp), intent(out) :: variance_history(:) !! Variance history.
      integer :: history, retained

      history = size(arch_history)
      if (origin >= history) then
         arch_history = filtered%innovation(origin - history + 1:origin)**2
         variance_history = filtered%conditional_sd(origin - history + 1:origin)**2
      else
         retained = history - origin
         if (retained > 0) then
            arch_history(1:retained) = initial_arch(size(initial_arch) - retained + 1:)
            variance_history(1:retained) = initial_variance(size(initial_variance) - retained + 1:)
         end if
         if (origin > 0) then
            arch_history(retained + 1:) = filtered%innovation(1:origin)**2
            variance_history(retained + 1:) = filtered%conditional_sd(1:origin)**2
         end if
      end if
   end subroutine garch_history_at_origin

   pure subroutine summarize_hresiduals(out)
      !! Summarize horizon errors using pairwise-complete forecast origins.
      type(tsissm_hresiduals_t), intent(inout) :: out !! Procedure result, updated in place.
      real(dp) :: mean_i, mean_j
      integer :: i, j, pair_total

      do i = 1, out%horizon
         pair_total = count(out%valid(:, i))
         if (pair_total > 0) then
            out%mean_residual(i) = sum(out%residual(:, i), mask=out%valid(:, i))/real(pair_total, dp)
            out%rmse(i) = sqrt(sum(out%residual(:, i)**2, &
               mask=out%valid(:, i))/real(pair_total, dp))
         end if
         do j = i, out%horizon
            pair_total = count(out%valid(:, i) .and. out%valid(:, j))
            out%pair_count(i, j) = pair_total
            out%pair_count(j, i) = pair_total
            if (pair_total > 1) then
               mean_i = sum(out%residual(:, i), mask=out%valid(:, i) .and. &
                  out%valid(:, j))/real(pair_total, dp)
               mean_j = sum(out%residual(:, j), mask=out%valid(:, i) .and. &
                  out%valid(:, j))/real(pair_total, dp)
               out%covariance(i, j) = sum((out%residual(:, i) - mean_i)* &
                  (out%residual(:, j) - mean_j), mask=out%valid(:, i) .and. &
                  out%valid(:, j))/real(pair_total - 1, dp)
               out%covariance(j, i) = out%covariance(i, j)
            end if
         end do
      end do
   end subroutine summarize_hresiduals

   pure subroutine box_cox_moments(mean, variance, lambda, corrected_mean, corrected_variance, &
      info, garch_variance, unconditional_variance)
      !! Apply tsissm Box-Cox mean and higher-order variance corrections.
      real(dp), intent(in) :: mean(:) !! Mean value or vector.
      real(dp), intent(in) :: variance(:) !! Variance value or matrix.
      real(dp), intent(in) :: lambda !! Penalty or shrinkage parameter.
      real(dp), intent(out) :: corrected_mean(:) !! Corrected mean.
      real(dp), intent(out) :: corrected_variance(:) !! Corrected variance.
      integer, intent(out) :: info !! Status code; zero indicates success.
      real(dp), intent(in), optional :: garch_variance(:) !! Garch variance.
      real(dp), intent(in), optional :: unconditional_variance !! Unconditional variance.
      real(dp) :: base, scale
      integer :: i

      info = 0
      do i = 1, size(mean)
         if (abs(lambda) < epsilon(1.0_dp)) then
            corrected_mean(i) = exp(mean(i))*(1.0_dp + 0.5_dp*variance(i))
            corrected_variance(i) = exp(2.0_dp*mean(i))*(exp(variance(i)) - 1.0_dp)
         else
            base = lambda*mean(i) + 1.0_dp
            if (base <= 0.0_dp) then
               info = 2
               return
            end if
            corrected_mean(i) = base**(1.0_dp/lambda)*(1.0_dp + &
               variance(i)*(1.0_dp - lambda)/(2.0_dp*base**2))
            corrected_variance(i) = variance(i)*base**(2.0_dp*(1.0_dp - lambda)/lambda) + &
               variance(i)**2*(1.0_dp - lambda)**2/(2.0_dp*base**4)* &
               base**(4.0_dp*(1.0_dp - lambda)/lambda - 4.0_dp)
            if (present(garch_variance) .and. present(unconditional_variance)) then
               scale = garch_variance(i)/max(unconditional_variance, tiny(1.0_dp))
               corrected_variance(i) = corrected_variance(i)*scale
            end if
         end if
      end do
   end subroutine box_cox_moments

   pure function tsissm_model(slope, damped_slope, seasonal_frequency, seasonal_harmonics, &
      regular_seasonal, ar_coefficients, ma_coefficients, level_persistence, slope_persistence, &
      seasonal_persistence, damping) result(out)
      !! Construct level, slope, multiple-seasonal, AR, and MA innovations matrices.
      logical, intent(in) :: slope !! Flag controlling slope.
      logical, intent(in) :: damped_slope !! Flag controlling damped slope.
      logical, intent(in) :: regular_seasonal !! Flag controlling regular seasonal.
      real(dp), intent(in) :: seasonal_frequency(:) !! Seasonal frequency.
      real(dp), intent(in) :: ar_coefficients(:) !! Autoregressive coefficients.
      real(dp), intent(in) :: ma_coefficients(:) !! Moving-average coefficients.
      integer, intent(in) :: seasonal_harmonics(:) !! Seasonal harmonics.
      real(dp), intent(in) :: level_persistence !! Level persistence.
      real(dp), intent(in) :: slope_persistence !! Slope persistence.
      real(dp), intent(in) :: seasonal_persistence(:) !! Seasonal persistence.
      real(dp), intent(in) :: damping !! Damping.
      type(tsissm_model_t) :: out
      real(dp) :: angle, phi
      integer :: ar_order, dimension, frequency_count, i, j, k, ma_order, offset, seasonal_dimension

      frequency_count = size(seasonal_frequency)
      ar_order = size(ar_coefficients)
      ma_order = size(ma_coefficients)
      if (size(seasonal_harmonics) /= frequency_count .or. &
         size(seasonal_persistence) /= frequency_count .or. &
         any(seasonal_frequency <= 1.0_dp) .or. any(seasonal_harmonics < 1) .or. &
         level_persistence < 0.0_dp .or. slope_persistence < 0.0_dp .or. &
         damping <= 0.0_dp .or. damping > 1.0_dp) then
         out%info = 1
         return
      end if
      if (regular_seasonal) then
         seasonal_dimension = sum(nint(seasonal_frequency))
         if (any(abs(seasonal_frequency - real(nint(seasonal_frequency), dp)) > &
            sqrt(epsilon(1.0_dp)))) then
            out%info = 1
            return
         end if
      else
         seasonal_dimension = 2*sum(seasonal_harmonics)
         if (any(2*seasonal_harmonics > nint(seasonal_frequency))) then
            out%info = 1
            return
         end if
      end if
      dimension = 1 + merge(1, 0, slope) + seasonal_dimension + ar_order + ma_order
      allocate(out%transition_base(dimension, dimension), out%transition_scale(dimension, dimension))
      allocate(out%transition_parameter(dimension, dimension), out%observation_loading(dimension))
      allocate(out%persistence(dimension), out%initial_state(dimension))
      allocate(out%seasonal_start(frequency_count), out%seasonal_end(frequency_count))
      out%transition_base = 0.0_dp
      out%transition_scale = 1.0_dp
      out%transition_parameter = 1.0_dp
      out%observation_loading = 0.0_dp
      out%persistence = 0.0_dp
      out%initial_state = 0.0_dp
      phi = 1.0_dp
      if (damped_slope) phi = damping
      out%transition_base(1, 1) = 1.0_dp
      out%observation_loading(1) = 1.0_dp
      out%persistence(1) = level_persistence
      offset = 1
      if (slope) then
         out%slope_index = 2
         out%transition_base(1, 2) = phi
         out%transition_base(2, 2) = phi
         out%observation_loading(2) = phi
         out%persistence(2) = slope_persistence
         offset = 2
      end if
      do i = 1, frequency_count
         out%seasonal_start(i) = offset + 1
         if (regular_seasonal) then
            k = nint(seasonal_frequency(i))
            out%transition_base(offset + 1, offset + k) = 1.0_dp
            do j = 2, k
               out%transition_base(offset + j, offset + j - 1) = 1.0_dp
            end do
            out%observation_loading(offset + k) = 1.0_dp
            out%persistence(offset + 1) = seasonal_persistence(i)
            offset = offset + k
         else
            k = seasonal_harmonics(i)
            do j = 1, k
               angle = 2.0_dp*acos(-1.0_dp)*real(j, dp)/seasonal_frequency(i)
               out%transition_base(offset + j, offset + j) = cos(angle)
               out%transition_base(offset + j, offset + k + j) = sin(angle)
               out%transition_base(offset + k + j, offset + j) = -sin(angle)
               out%transition_base(offset + k + j, offset + k + j) = cos(angle)
               out%observation_loading(offset + j) = 1.0_dp
               out%persistence(offset + j) = seasonal_persistence(i)
               out%persistence(offset + k + j) = seasonal_persistence(i)
            end do
            offset = offset + 2*k
         end if
         out%seasonal_end(i) = offset
      end do
      if (ar_order > 0) then
         out%ar_start = offset + 1
         out%ar_end = offset + ar_order
         do j = 1, ar_order
            out%transition_base(1:offset, offset + j) = out%persistence(1:offset)*ar_coefficients(j)
            out%transition_base(offset + 1, offset + j) = ar_coefficients(j)
            out%observation_loading(offset + j) = ar_coefficients(j)
         end do
         out%persistence(offset + 1) = 1.0_dp
         do j = 2, ar_order
            out%transition_base(offset + j, offset + j - 1) = 1.0_dp
         end do
         offset = offset + ar_order
      end if
      if (ma_order > 0) then
         out%ma_start = offset + 1
         out%ma_end = offset + ma_order
         do j = 1, ma_order
            out%transition_base(1:offset, offset + j) = out%persistence(1:offset)*ma_coefficients(j)
            out%observation_loading(offset + j) = ma_coefficients(j)
         end do
         out%persistence(offset + 1) = 1.0_dp
         do j = 2, ma_order
            out%transition_base(offset + j, offset + j - 1) = 1.0_dp
         end do
      end if
   end function tsissm_model

   pure function tsissm_filter_constant(observations, transition_base, transition_scale, &
      transition_parameter, observation_loading, persistence, regressors, coefficients, &
      initial_state, observed, lambda) result(out)
      !! Filter a constant-variance single-source-of-error innovations model.
      real(dp), intent(in) :: observations(:) !! Observed time-series values.
      real(dp), intent(in) :: transition_base(:, :) !! Transition base.
      real(dp), intent(in) :: transition_scale(:, :) !! Transition scale.
      real(dp), intent(in) :: transition_parameter(:, :) !! Transition parameter.
      real(dp), intent(in) :: observation_loading(:) !! Observation loading matrix.
      real(dp), intent(in) :: persistence(:) !! Persistence.
      real(dp), intent(in) :: regressors(:, :) !! Regression design matrix.
      real(dp), intent(in) :: coefficients(:) !! Model coefficients.
      real(dp), intent(in) :: initial_state(:) !! Initial state vector.
      real(dp), intent(in) :: lambda !! Penalty or shrinkage parameter.
      logical, intent(in), optional :: observed(:) !! Flag controlling observed.
      type(tsissm_filter_t) :: out
      logical, allocatable :: use_observation(:)
      integer :: dimension, i, n

      n = size(observations)
      dimension = size(initial_state)
      if (n < 1 .or. dimension < 1 .or. size(transition_base, 1) /= dimension .or. &
         size(transition_base, 2) /= dimension .or. &
         any(shape(transition_scale) /= [dimension, dimension]) .or. &
         any(shape(transition_parameter) /= [dimension, dimension]) .or. &
         size(observation_loading) /= dimension .or. size(persistence) /= dimension .or. &
         size(regressors, 1) /= n .or. size(regressors, 2) /= size(coefficients) .or. &
         any(observations <= 0.0_dp)) then
         out%info = 1
         return
      end if
      if (present(observed)) then
         if (size(observed) /= n) then
            out%info = 1
            return
         end if
      end if
      allocate(use_observation(n))
      use_observation = .true.
      if (present(observed)) use_observation = observed
      allocate(out%state(dimension, 0:n), out%fitted(n), out%transformed(n), out%innovation(n))
      allocate(out%transition(dimension, dimension), out%discount_transition(dimension, dimension))
      allocate(out%observation_loading(dimension), out%persistence(dimension))
      out%transition = transition_base*transition_scale*transition_parameter
      out%discount_transition = out%transition - outer_product(persistence, observation_loading)
      out%observation_loading = observation_loading
      out%persistence = persistence
      out%state(:, 0) = initial_state
      out%transformed = box_cox(observations, lambda)
      do i = 1, n
         out%fitted(i) = dot_product(observation_loading, out%state(:, i - 1)) + &
            dot_product(regressors(i, :), coefficients)
         if (use_observation(i)) then
            out%innovation(i) = out%transformed(i) - out%fitted(i)
            out%observations = out%observations + 1
         else
            out%innovation(i) = 0.0_dp
         end if
         out%state(:, i) = matmul(out%transition, out%state(:, i - 1)) + &
            persistence*out%innovation(i)
      end do
   end function tsissm_filter_constant

   pure function tsissm_predict_constant(transition_base, transition_scale, transition_parameter, &
      observation_loading, persistence, regressors, coefficients, initial_state, &
      innovations) result(out)
      !! Simulate constant-variance ISSM forecasts from supplied innovations.
      real(dp), intent(in) :: transition_base(:, :) !! Transition base.
      real(dp), intent(in) :: transition_scale(:, :) !! Transition scale.
      real(dp), intent(in) :: transition_parameter(:, :) !! Transition parameter.
      real(dp), intent(in) :: observation_loading(:) !! Observation loading matrix.
      real(dp), intent(in) :: persistence(:) !! Persistence.
      real(dp), intent(in) :: regressors(:, :) !! Regression design matrix.
      real(dp), intent(in) :: coefficients(:) !! Model coefficients.
      real(dp), intent(in) :: initial_state(:) !! Initial state vector.
      real(dp), intent(in) :: innovations(:, :) !! Model innovations.
      type(tsissm_prediction_t) :: out
      real(dp), allocatable :: transition(:, :)
      integer :: dimension, horizon, i, simulation

      dimension = size(initial_state)
      horizon = size(innovations, 2)
      if (dimension < 1 .or. horizon < 1 .or. size(innovations, 1) < 1 .or. &
         any(shape(transition_base) /= [dimension, dimension]) .or. &
         any(shape(transition_scale) /= [dimension, dimension]) .or. &
         any(shape(transition_parameter) /= [dimension, dimension]) .or. &
         size(observation_loading) /= dimension .or. size(persistence) /= dimension .or. &
         size(regressors, 1) /= horizon .or. size(regressors, 2) /= size(coefficients)) then
         out%info = 1
         return
      end if
      transition = transition_base*transition_scale*transition_parameter
      allocate(out%observation(size(innovations, 1), horizon))
      allocate(out%state(dimension, 0:horizon, size(innovations, 1)))
      allocate(out%innovation(size(innovations, 1), horizon))
      out%innovation = innovations
      do simulation = 1, size(innovations, 1)
         out%state(:, 0, simulation) = initial_state
         do i = 1, horizon
            out%observation(simulation, i) = dot_product(observation_loading, &
               out%state(:, i - 1, simulation)) + dot_product(regressors(i, :), coefficients) + &
               innovations(simulation, i)
            out%state(:, i, simulation) = matmul(transition, out%state(:, i - 1, simulation)) + &
               persistence*innovations(simulation, i)
         end do
      end do
   end function tsissm_predict_constant

   pure function tsissm_garch_recursion(innovations, arch, garch, initial_arch, &
      initial_variance, variance_intercept) result(out)
      !! Compute GARCH variance using chronological pre-sample histories.
      real(dp), intent(in) :: innovations(:) !! Model innovations.
      real(dp), intent(in) :: arch(:) !! Arch.
      real(dp), intent(in) :: garch(:) !! Garch.
      real(dp), intent(in) :: initial_arch(:) !! Initial arch.
      real(dp), intent(in) :: initial_variance(:) !! Initial variance.
      real(dp), intent(in) :: variance_intercept !! Variance intercept.
      type(tsissm_garch_t) :: out
      integer :: history, i, j, lag

      history = max(size(arch), size(garch))
      if (size(innovations) < 1 .or. history < 1 .or. size(initial_arch) < history .or. &
         size(initial_variance) < history .or. variance_intercept < 0.0_dp .or. &
         any(arch < 0.0_dp) .or. any(garch < 0.0_dp)) then
         out%info = 1
         return
      end if
      allocate(out%variance(size(innovations)), out%standard_deviation(size(innovations)))
      do i = 1, size(innovations)
         out%variance(i) = variance_intercept
         do j = 1, size(arch)
            lag = i - j
            if (lag >= 1) then
               out%variance(i) = out%variance(i) + arch(j)*innovations(lag)**2
            else
               out%variance(i) = out%variance(i) + arch(j)*initial_arch(history + lag)
            end if
         end do
         do j = 1, size(garch)
            lag = i - j
            if (lag >= 1) then
               out%variance(i) = out%variance(i) + garch(j)*out%variance(lag)
            else
               out%variance(i) = out%variance(i) + garch(j)*initial_variance(history + lag)
            end if
         end do
         if (out%variance(i) < 0.0_dp) then
            out%info = 2
            return
         end if
         out%standard_deviation(i) = sqrt(out%variance(i))
      end do
   end function tsissm_garch_recursion

   pure function tsissm_filter_dynamic(observations, transition_base, transition_scale, &
      transition_parameter, observation_loading, persistence, regressors, coefficients, &
      initial_state, observed, lambda, arch, garch, initial_arch, initial_variance, &
      variance_intercept) result(out)
      !! Filter an ISSM and attach its observation-driven GARCH scale path.
      real(dp), intent(in) :: observations(:) !! Observed time-series values.
      real(dp), intent(in) :: transition_base(:, :) !! Transition base.
      real(dp), intent(in) :: transition_scale(:, :) !! Transition scale.
      real(dp), intent(in) :: transition_parameter(:, :) !! Transition parameter.
      real(dp), intent(in) :: observation_loading(:) !! Observation loading matrix.
      real(dp), intent(in) :: persistence(:) !! Persistence.
      real(dp), intent(in) :: regressors(:, :) !! Regression design matrix.
      real(dp), intent(in) :: coefficients(:) !! Model coefficients.
      real(dp), intent(in) :: initial_state(:) !! Initial state vector.
      real(dp), intent(in) :: lambda !! Penalty or shrinkage parameter.
      logical, intent(in), optional :: observed(:) !! Flag controlling observed.
      real(dp), intent(in) :: arch(:) !! Arch.
      real(dp), intent(in) :: garch(:) !! Garch.
      real(dp), intent(in) :: initial_arch(:) !! Initial arch.
      real(dp), intent(in) :: initial_variance(:) !! Initial variance.
      real(dp), intent(in) :: variance_intercept !! Variance intercept.
      type(tsissm_filter_t) :: out
      type(tsissm_garch_t) :: variance_path

      out = tsissm_filter_constant(observations, transition_base, transition_scale, &
         transition_parameter, observation_loading, persistence, regressors, coefficients, &
         initial_state, observed, lambda)
      if (out%info /= 0) return
      variance_path = tsissm_garch_recursion(out%innovation, arch, garch, initial_arch, &
         initial_variance, variance_intercept)
      if (variance_path%info /= 0) then
         out%info = 10 + variance_path%info
         return
      end if
      out%conditional_sd = variance_path%standard_deviation
   end function tsissm_filter_dynamic

   pure function tsissm_predict_dynamic(transition_base, transition_scale, transition_parameter, &
      observation_loading, persistence, regressors, coefficients, initial_state, &
      standardized_innovations, arch, garch, initial_arch, initial_variance, &
      variance_intercept) result(out)
      !! Simulate ISSM forecasts with path-dependent GARCH innovations.
      real(dp), intent(in) :: transition_base(:, :) !! Transition base.
      real(dp), intent(in) :: transition_scale(:, :) !! Transition scale.
      real(dp), intent(in) :: transition_parameter(:, :) !! Transition parameter.
      real(dp), intent(in) :: observation_loading(:) !! Observation loading matrix.
      real(dp), intent(in) :: persistence(:) !! Persistence.
      real(dp), intent(in) :: regressors(:, :) !! Regression design matrix.
      real(dp), intent(in) :: coefficients(:) !! Model coefficients.
      real(dp), intent(in) :: initial_state(:) !! Initial state vector.
      real(dp), intent(in) :: standardized_innovations(:, :) !! Standardized innovations.
      real(dp), intent(in) :: arch(:) !! Arch.
      real(dp), intent(in) :: garch(:) !! Garch.
      real(dp), intent(in) :: initial_arch(:) !! Initial arch.
      real(dp), intent(in) :: initial_variance(:) !! Initial variance.
      real(dp), intent(in) :: variance_intercept !! Variance intercept.
      type(tsissm_prediction_t) :: out
      real(dp), allocatable :: innovation_history(:), variance_history(:), transition(:, :)
      integer :: dimension, history, horizon, i, j, lag, simulation

      dimension = size(initial_state)
      horizon = size(standardized_innovations, 2)
      history = max(size(arch), size(garch))
      if (dimension < 1 .or. horizon < 1 .or. size(standardized_innovations, 1) < 1 .or. &
         history < 1 .or. size(initial_arch) < history .or. size(initial_variance) < history .or. &
         any(shape(transition_base) /= [dimension, dimension]) .or. &
         any(shape(transition_scale) /= [dimension, dimension]) .or. &
         any(shape(transition_parameter) /= [dimension, dimension]) .or. &
         size(observation_loading) /= dimension .or. size(persistence) /= dimension .or. &
         size(regressors, 1) /= horizon .or. size(regressors, 2) /= size(coefficients) .or. &
         variance_intercept < 0.0_dp .or. any(arch < 0.0_dp) .or. any(garch < 0.0_dp)) then
         out%info = 1
         return
      end if
      transition = transition_base*transition_scale*transition_parameter
      allocate(out%observation(size(standardized_innovations, 1), horizon))
      allocate(out%innovation(size(standardized_innovations, 1), horizon))
      allocate(out%conditional_sd(size(standardized_innovations, 1), horizon))
      allocate(out%state(dimension, 0:horizon, size(standardized_innovations, 1)))
      allocate(innovation_history(history + horizon), variance_history(history + horizon))
      do simulation = 1, size(standardized_innovations, 1)
         innovation_history = 0.0_dp
         variance_history = 0.0_dp
         innovation_history(1:history) = sqrt(max(initial_arch(1:history), 0.0_dp))
         variance_history(1:history) = initial_variance(1:history)
         out%state(:, 0, simulation) = initial_state
         do i = 1, horizon
            variance_history(history + i) = variance_intercept
            do j = 1, size(arch)
               lag = history + i - j
               variance_history(history + i) = variance_history(history + i) + &
                  arch(j)*innovation_history(lag)**2
            end do
            do j = 1, size(garch)
               lag = history + i - j
               variance_history(history + i) = variance_history(history + i) + &
                  garch(j)*variance_history(lag)
            end do
            if (variance_history(history + i) < 0.0_dp) then
               out%info = 2
               return
            end if
            out%conditional_sd(simulation, i) = sqrt(variance_history(history + i))
            out%innovation(simulation, i) = standardized_innovations(simulation, i)* &
               out%conditional_sd(simulation, i)
            innovation_history(history + i) = out%innovation(simulation, i)
            out%observation(simulation, i) = dot_product(observation_loading, &
               out%state(:, i - 1, simulation)) + dot_product(regressors(i, :), coefficients) + &
               out%innovation(simulation, i)
            out%state(:, i, simulation) = matmul(transition, out%state(:, i - 1, simulation)) + &
               persistence*out%innovation(simulation, i)
         end do
      end do
   end function tsissm_predict_dynamic

   pure function pack_rows(matrix, mask) result(packed)
      !! Pack selected matrix rows while retaining every column.
      real(dp), intent(in) :: matrix(:, :) !! Input matrix.
      logical, intent(in) :: mask(:) !! Flag controlling mask.
      real(dp), allocatable :: packed(:, :)
      integer :: i, row

      allocate(packed(count(mask), size(matrix, 2)))
      row = 0
      do i = 1, size(mask)
         if (mask(i)) then
            row = row + 1
            packed(row, :) = matrix(i, :)
         end if
      end do
   end function pack_rows

   pure function matrix_polynomial_roots(matrix) result(roots)
      !! Compute matrix eigenvalues from Faddeev-LeVerrier coefficients.
      real(dp), intent(in) :: matrix(:, :) !! Input matrix.
      complex(dp), allocatable :: roots(:)
      real(dp) :: coefficient(size(matrix, 1)), work(size(matrix, 1), size(matrix, 1))
      real(dp) :: identity(size(matrix, 1), size(matrix, 1)), polynomial(size(matrix, 1) + 1)
      integer :: i, j, n

      n = size(matrix, 1)
      identity = 0.0_dp
      do i = 1, n
         identity(i, i) = 1.0_dp
      end do
      work = identity
      do i = 1, n
         work = matmul(matrix, work)
         coefficient(i) = -sum([(work(j, j), j=1, n)])/real(i, dp)
         work = work + coefficient(i)*identity
      end do
      polynomial(n + 1) = 1.0_dp
      do i = 1, n
         polynomial(i) = coefficient(n - i + 1)
      end do
      roots = polynomial_roots(polynomial)
   end function matrix_polynomial_roots

   pure function polynomial_roots(polynomial) result(roots)
      !! Find all roots of a real ascending-power polynomial by Durand-Kerner iteration.
      real(dp), intent(in) :: polynomial(:) !! Polynomial.
      complex(dp), allocatable :: roots(:)
      complex(dp), allocatable :: next_roots(:)
      complex(dp) :: denominator, value
      real(dp) :: root_bound
      integer :: degree, i, iteration, j

      degree = size(polynomial) - 1
      allocate(roots(max(0, degree)), next_roots(max(0, degree)))
      if (degree < 1) return
      if (abs(polynomial(degree + 1)) <= tiny(1.0_dp)) then
         roots = cmplx(huge(1.0_dp), 0.0_dp, dp)
         return
      end if
      root_bound = 1.0_dp + maxval(abs(polynomial(:degree)/polynomial(degree + 1)))
      do i = 1, degree
         roots(i) = root_bound*exp(cmplx(0.0_dp, 2.0_dp*acos(-1.0_dp)* &
            real(i - 1, dp)/real(degree, dp), dp))
      end do
      do iteration = 1, 500
         do i = 1, degree
            value = cmplx(polynomial(degree + 1), 0.0_dp, dp)
            do j = degree, 1, -1
               value = value*roots(i) + polynomial(j)
            end do
            denominator = cmplx(polynomial(degree + 1), 0.0_dp, dp)
            do j = 1, degree
               if (j /= i) denominator = denominator*(roots(i) - roots(j))
            end do
            if (abs(denominator) <= tiny(1.0_dp)) denominator = cmplx(tiny(1.0_dp), 0.0_dp, dp)
            next_roots(i) = roots(i) - value/denominator
         end do
         if (maxval(abs(next_roots - roots)) <= 100.0_dp*epsilon(1.0_dp)) exit
         roots = next_roots
      end do
      roots = next_roots
   end function polynomial_roots

   pure real(dp) function spectral_radius(matrix) result(radius)
      !! Compute the largest matrix eigenvalue modulus.
      real(dp), intent(in) :: matrix(:, :) !! Input matrix.
      complex(dp), allocatable :: roots(:)

      roots = matrix_polynomial_roots(matrix)
      radius = maxval(abs(roots))
   end function spectral_radius
end module tsissm_mod
