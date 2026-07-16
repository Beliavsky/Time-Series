! SPDX-License-Identifier: GPL-2.0-only
! SPDX-FileComment: Regression tests for the tsissm translation.
program test_tsissm
   use kind_mod, only: dp
   use forecast_mod, only: acf_values, pacf_values
   use tsissm_mod
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   type(tsissm_filter_t) :: filtered
   type(tsissm_prediction_t) :: predicted
   type(tsissm_filter_t) :: dynamic_filtered
   type(tsissm_prediction_t) :: dynamic_predicted
   type(tsissm_garch_t) :: variance_path
   type(tsissm_variance_initialization_t) :: full_variance, sample_variance
   type(tsissm_model_t) :: structural, regular
   type(tsissm_likelihood_t) :: gaussian_likelihood, student_likelihood, jsu_likelihood
   type(tsissm_constraint_t) :: constraints, bad_constraints
   type(tsissm_initialization_t) :: initialized
   type(tsissm_fit_t) :: estimated
   type(tsissm_fit_t) :: dynamic_estimated
   type(tsissm_structural_fit_t) :: structural_estimated
   type(tsissm_structural_fit_t) :: dynamic_structural_estimated
   type(tsissm_moments_t) :: constant_moments, dynamic_moments, log_moments
   type(tsissm_hresiduals_t) :: constant_hresiduals, dynamic_hresiduals
   type(tsissm_covariance_t) :: hessian_covariance, opg_covariance
   type(tsissm_covariance_t) :: qmle_covariance, newey_west_covariance
   type(tsissm_covariance_t) :: parzen_covariance, qs_covariance
   type(tsissm_covariance_t) :: truncated_covariance, tukey_covariance, prewhite_covariance
   type(tsissm_covariance_t) :: constant_fit_covariance, dynamic_fit_covariance
   type(tsissm_covariance_t) :: fitted_hac_covariance
   type(tsissm_covariance_t) :: structural_fit_covariance
   type(tsissm_covariance_t) :: dynamic_structural_covariance
   type(tsissm_filter_t) :: covariance_dynamic_filter
   type(tsissm_profile_t) :: constant_profile, dynamic_profile
   type(tsissm_profile_t) :: structural_constant_profile, structural_dynamic_profile
   type(tsissm_filter_t) :: complete_filter
   type(tsissm_decomposition_t) :: decomposition
   type(tsissm_prediction_decomposition_t) :: prediction_decomposition
   type(tsissm_filter_t) :: structural_filter
   type(tsissm_prediction_t) :: structural_prediction
   type(tsissm_forecast_t) :: constant_forecast, repeated_forecast, bootstrap_forecast
   type(tsissm_forecast_t) :: student_forecast, jsu_forecast, dynamic_forecast
   type(tsissm_diagnostics_t) :: diagnostics
   type(tsissm_structural_diagnostics_t) :: structural_diagnostics
   type(tsissm_structural_diagnostics_t) :: outlier_diagnostics
   type(tsissm_weighted_box_test_t) :: box_tests(6), transformed_test, direct_test
   type(tsissm_filter_t) :: outlier_filter
   type(tsissm_accuracy_t) :: forecast_scores
   type(tsissm_backtest_t) :: constant_backtest, refitted_backtest, dynamic_backtest
   type(tsissm_backtest_t) :: structural_constant_backtest, structural_dynamic_backtest
   type(tsissm_selection_t) :: selected_models, variance_selection
   type(tsissm_selection_t) :: reranked_models
   type(tsissm_selection_t) :: compatible_selection
   type(tsissm_ensemble_t) :: ensemble_forecast
   real(dp) :: transition(1, 1), loading(1), persistence(1), initial(1)
   real(dp) :: regressors(3, 0), coefficients(0), future_regressors(2, 0)
   real(dp) :: expected
   real(dp) :: fit_regressors(8, 0)
   real(dp) :: moment_regressors(3, 0)
   real(dp) :: complete_regressors(4, 0)
   real(dp) :: covariance_regressors(8, 0)
   real(dp) :: profile_regressors(10, 0)
   real(dp) :: structural_backtest_regressors(12, 0)
   real(dp) :: decomposition_regressors(3, 1), prediction_regressors(2, 1)
   integer :: i
   logical :: observed(3)
   real(dp), allocatable :: ensemble_components(:, :), ensemble_expected(:, :)
   real(dp), allocatable :: ensemble_distribution_expected(:, :)
   real(dp) :: score_minus(2, 2), score_plus(2, 2), score_step(2)
   real(dp) :: score_matrix(4, 2), score_hessian(2, 2)
   real(dp), allocatable :: finite_scores(:, :)
   real(dp), allocatable :: box_acf(:), box_pacf(:)
   real(dp) :: box_residuals(8)

   transition = 1.0_dp
   loading = 1.0_dp
   persistence = 0.5_dp
   initial = 0.0_dp
   observed = [.true., .false., .true.]
   filtered = tsissm_filter_constant([1.0_dp, 2.0_dp, 3.0_dp], transition, transition, &
      transition, loading, persistence, regressors, coefficients, initial, observed, 1.0_dp)
   call check(filtered%info == 0 .and. filtered%observations == 2 .and. &
      all(abs(filtered%transformed - [0.0_dp, 1.0_dp, 2.0_dp]) < 1.0e-12_dp) .and. &
      all(abs(filtered%innovation - [0.0_dp, 0.0_dp, 2.0_dp]) < 1.0e-12_dp) .and. &
      abs(filtered%state(1, 3) - 1.0_dp) < 1.0e-12_dp .and. &
      abs(filtered%discount_transition(1, 1) - 0.5_dp) < 1.0e-12_dp, &
      'constant innovations filter')

   predicted = tsissm_predict_constant(transition, transition, transition, loading, persistence, &
      future_regressors, coefficients, initial, reshape([1.0_dp, 2.0_dp], [1, 2]))
   call check(predicted%info == 0 .and. all(shape(predicted%observation) == [1, 2]) .and. &
      abs(predicted%observation(1, 1) - 1.0_dp) < 1.0e-12_dp .and. &
      abs(predicted%observation(1, 2) - 2.5_dp) < 1.0e-12_dp .and. &
      abs(predicted%state(1, 2, 1) - 1.5_dp) < 1.0e-12_dp, &
      'constant innovations prediction')

   variance_path = tsissm_garch_recursion([2.0_dp, 3.0_dp], [0.2_dp], [0.5_dp], &
      [4.0_dp], [9.0_dp], 1.0_dp)
   call check(variance_path%info == 0 .and. &
      maxval(abs(variance_path%variance - [6.3_dp, 4.95_dp])) < 1.0e-12_dp, &
      'GARCH variance recursion')

   full_variance = tsissm_initialize_variance([1.0_dp, 2.0_dp, 3.0_dp])
   sample_variance = tsissm_initialize_variance([1.0_dp, 2.0_dp, 3.0_dp], &
      method=tsissm_variance_initial_sample, sample_size=2)
   call check(full_variance%info == 0 .and. full_variance%observations == 3 .and. &
      abs(full_variance%variance - 14.0_dp/3.0_dp) < 1.0e-12_dp .and. &
      sample_variance%info == 0 .and. sample_variance%observations == 2 .and. &
      abs(sample_variance%variance - 2.5_dp) < 1.0e-12_dp, &
      'configurable GARCH variance initialization')
   sample_variance = tsissm_initialize_variance([1.0_dp, 2.0_dp, 3.0_dp], &
      [.false., .true., .true.], tsissm_variance_initial_sample, 1)
   call check(sample_variance%info == 0 .and. sample_variance%observations == 1 .and. &
      abs(sample_variance%variance - 4.0_dp) < 1.0e-12_dp, &
      'valid-observation GARCH initialization window')

   dynamic_filtered = tsissm_filter_dynamic([1.0_dp, 2.0_dp, 3.0_dp], transition, transition, &
      transition, loading, persistence, regressors, coefficients, initial, lambda=1.0_dp, &
      arch=[0.2_dp], garch=[0.5_dp], initial_arch=[4.0_dp], initial_variance=[9.0_dp], &
      variance_intercept=1.0_dp)
   call check(dynamic_filtered%info == 0 .and. size(dynamic_filtered%conditional_sd) == 3 .and. &
      abs(dynamic_filtered%conditional_sd(1)**2 - 6.3_dp) < 1.0e-12_dp, &
      'dynamic innovations filter')

   complete_filter = tsissm_filter_constant([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], &
      transition, transition, transition, loading, persistence, complete_regressors, &
      coefficients, initial, lambda=1.0_dp)
   constant_hresiduals = tsissm_hresiduals_constant([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], &
      complete_filter, complete_regressors, coefficients, 1.0_dp, 1.0_dp, 3, 0, .true.)
   call check(constant_hresiduals%info == 0 .and. &
      all(shape(constant_hresiduals%residual) == [4, 3]) .and. &
      all(constant_hresiduals%pair_count == reshape([4, 3, 2, 3, 3, 2, 2, 2, 2], [3, 3])) .and. &
      maxval(abs(constant_hresiduals%residual(:, 1) - complete_filter%innovation)) < 1.0e-12_dp .and. &
      maxval(abs(constant_hresiduals%covariance - &
      transpose(constant_hresiduals%covariance))) < 1.0e-12_dp, &
      'constant multi-horizon residuals')

   dynamic_hresiduals = tsissm_hresiduals_dynamic([1.0_dp, 2.0_dp, 3.0_dp], &
      dynamic_filtered, regressors, coefficients, [0.2_dp], [0.5_dp], 1.0_dp, &
      [4.0_dp], [9.0_dp], 1.0_dp, 2, 0, .true.)
   call check(dynamic_hresiduals%info == 0 .and. &
      all(shape(dynamic_hresiduals%residual) == [3, 2]) .and. &
      all(dynamic_hresiduals%pair_count == reshape([3, 2, 2, 2], [2, 2])) .and. &
      maxval(abs(dynamic_hresiduals%residual(:, 1) - dynamic_filtered%innovation)) < 1.0e-12_dp .and. &
      all(ieee_is_finite(dynamic_hresiduals%covariance)), &
      'dynamic multi-horizon residuals')

   score_minus = reshape([0.0_dp, 1.0_dp, 2.0_dp, 4.0_dp], [2, 2])
   score_plus = reshape([2.0_dp, 5.0_dp, 8.0_dp, 12.0_dp], [2, 2])
   score_step = [1.0_dp, 2.0_dp]
   finite_scores = tsissm_fd_scores(score_minus, score_plus, score_step)
   call check(all(shape(finite_scores) == [2, 2]) .and. &
      maxval(abs(finite_scores - reshape([1.0_dp, 2.0_dp, 1.5_dp, 2.0_dp], [2, 2]))) < 1.0e-12_dp, &
      'finite-difference observation scores')

   score_matrix = reshape([1.0_dp, -1.0_dp, 1.0_dp, -1.0_dp, &
      0.0_dp, 1.0_dp, -1.0_dp, 0.0_dp], [4, 2])
   score_hessian = reshape([4.0_dp, 0.0_dp, 0.0_dp, 8.0_dp], [2, 2])
   hessian_covariance = tsissm_parameter_covariance(score_matrix, score_hessian, &
      tsissm_covariance_hessian)
   opg_covariance = tsissm_parameter_covariance(score_matrix, score_hessian, &
      tsissm_covariance_opg)
   qmle_covariance = tsissm_parameter_covariance(score_matrix, score_hessian, &
      tsissm_covariance_qmle, .true.)
   newey_west_covariance = tsissm_parameter_covariance(score_matrix, score_hessian, &
      tsissm_covariance_newey_west, newey_west_lag=1)
   parzen_covariance = tsissm_parameter_covariance(score_matrix, score_hessian, &
      tsissm_covariance_newey_west, newey_west_lag=1, hac_kernel=tsissm_hac_parzen)
   qs_covariance = tsissm_parameter_covariance(score_matrix, score_hessian, &
      tsissm_covariance_newey_west, hac_kernel=tsissm_hac_quadratic_spectral)
   truncated_covariance = tsissm_parameter_covariance(score_matrix, score_hessian, &
      tsissm_covariance_newey_west, newey_west_lag=1, hac_kernel=tsissm_hac_truncated)
   tukey_covariance = tsissm_parameter_covariance(score_matrix, score_hessian, &
      tsissm_covariance_newey_west, hac_kernel=tsissm_hac_tukey_hanning, hac_bandwidth=3.0_dp)
   prewhite_covariance = tsissm_parameter_covariance(score_matrix, score_hessian, &
      tsissm_covariance_newey_west, newey_west_lag=0, &
      hac_kernel=tsissm_hac_bartlett, prewhite_order=1)
   call check(hessian_covariance%info == 0 .and. &
      maxval(abs(hessian_covariance%covariance - &
      reshape([0.25_dp, 0.0_dp, 0.0_dp, 0.125_dp], [2, 2]))) < 1.0e-12_dp .and. &
      maxval(abs(matmul(opg_covariance%meat, opg_covariance%covariance) - &
      reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2]))) < 1.0e-12_dp .and. &
      maxval(abs(qmle_covariance%covariance - matmul(matmul(qmle_covariance%bread, &
      qmle_covariance%meat), qmle_covariance%bread))) < 1.0e-12_dp .and. &
      newey_west_covariance%info == 0 .and. newey_west_covariance%lag == 1 .and. &
      maxval(abs(newey_west_covariance%covariance - &
      transpose(newey_west_covariance%covariance))) < 1.0e-12_dp, &
      'parameter covariance estimators')
   call check(parzen_covariance%info == 0 .and. &
      parzen_covariance%kernel == tsissm_hac_parzen .and. &
      qs_covariance%info == 0 .and. qs_covariance%lag == 3 .and. &
      qs_covariance%bandwidth >= 1.0_dp .and. &
      truncated_covariance%info == 0 .and. &
      truncated_covariance%kernel == tsissm_hac_truncated .and. &
      tukey_covariance%info == 0 .and. tukey_covariance%lag == 2 .and. &
      prewhite_covariance%info == 0 .and. prewhite_covariance%prewhite_order == 1 .and. &
      prewhite_covariance%bias_correction > 0.0_dp .and. &
      prewhite_covariance%effective_degrees > 0.0_dp .and. &
      all(ieee_is_finite(prewhite_covariance%covariance)), &
      'HAC kernels bandwidth and prewhitening')

   constant_fit_covariance = tsissm_covariance_constant([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], &
      complete_filter, complete_regressors, persistence, coefficients, 1.0_dp, &
      tsissm_distribution_gaussian, 0.0_dp, 0.0_dp, tsissm_covariance_hessian)
   call check(constant_fit_covariance%info == 0 .and. &
      all(shape(constant_fit_covariance%scores) == [4, 1]) .and. &
      all(shape(constant_fit_covariance%hessian) == [1, 1]) .and. &
      all(ieee_is_finite(constant_fit_covariance%covariance)), &
      'constant fitted-model covariance')
   fitted_hac_covariance = tsissm_covariance_constant([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], &
      complete_filter, complete_regressors, persistence, coefficients, 1.0_dp, &
      tsissm_distribution_gaussian, 0.0_dp, 0.0_dp, tsissm_covariance_newey_west, &
      newey_west_lag=1, hac_kernel=tsissm_hac_parzen)
   call check(fitted_hac_covariance%info == 0 .and. &
      fitted_hac_covariance%kernel == tsissm_hac_parzen .and. &
      fitted_hac_covariance%lag == 1 .and. fitted_hac_covariance%bandwidth == 2.0_dp .and. &
      all(ieee_is_finite(fitted_hac_covariance%covariance)), &
      'fitted-model HAC covariance controls')

   covariance_dynamic_filter = tsissm_filter_dynamic([1.0_dp, 1.8_dp, 1.4_dp, 2.1_dp, &
      1.7_dp, 2.4_dp, 2.0_dp, 2.6_dp], transition, transition, transition, loading, &
      [0.3_dp], covariance_regressors, coefficients, initial, lambda=1.0_dp, &
      arch=[0.1_dp], garch=[0.5_dp], initial_arch=[1.0_dp], initial_variance=[1.0_dp], &
      variance_intercept=0.4_dp)
   dynamic_fit_covariance = tsissm_covariance_dynamic([1.0_dp, 1.8_dp, 1.4_dp, 2.1_dp, &
      1.7_dp, 2.4_dp, 2.0_dp, 2.6_dp], covariance_dynamic_filter, covariance_regressors, &
      [0.3_dp], coefficients, [0.1_dp], [0.5_dp], 1.0_dp, tsissm_distribution_gaussian, &
      0.0_dp, 0.0_dp, tsissm_covariance_qmle)
   call check(dynamic_fit_covariance%info == 0 .and. &
      all(shape(dynamic_fit_covariance%scores) == [8, 3]) .and. &
      all(shape(dynamic_fit_covariance%hessian) == [3, 3]) .and. &
      all(ieee_is_finite(dynamic_fit_covariance%covariance)), &
      'dynamic fitted-model covariance')

   constant_profile = tsissm_profile_constant(transition, transition, transition, loading, &
      [0.2_dp], profile_regressors, coefficients, [1.0_dp], 0.05_dp, 1.0_dp, &
      tsissm_distribution_gaussian, 0.0_dp, 0.0_dp, 2, 2, 20, 1, [integer ::], &
      40, 1.0e-5_dp, 713)
   call check(constant_profile%info == 0 .and. constant_profile%successful == 2 .and. &
      all(constant_profile%valid) .and. all(shape(constant_profile%parameter) == [2, 1]) .and. &
      all(shape(constant_profile%actual) == [2, 2]) .and. &
      all(ieee_is_finite(constant_profile%mape)) .and. &
      all(ieee_is_finite(constant_profile%crps)), &
      'constant simulation-refit profile')

   dynamic_profile = tsissm_profile_dynamic(transition, transition, transition, loading, &
      [0.2_dp], profile_regressors, coefficients, [1.0_dp], [0.1_dp], [0.5_dp], &
      0.01_dp, [0.025_dp], [0.025_dp], 1.0_dp, tsissm_distribution_gaussian, &
      0.0_dp, 0.0_dp, 2, 1, 20, 1, [integer ::], 50, 1.0e-5_dp, 917)
   call check(dynamic_profile%info == 0 .and. dynamic_profile%successful == 1 .and. &
      dynamic_profile%valid(1) .and. all(shape(dynamic_profile%parameter) == [1, 3]) .and. &
      all(ieee_is_finite(dynamic_profile%predicted)) .and. &
      all(ieee_is_finite(dynamic_profile%mase)) .and. &
      all(ieee_is_finite(dynamic_profile%crps)), &
      'dynamic simulation-refit profile')

   structural_constant_profile = tsissm_profile_structural_constant(profile_regressors, &
      .false., .false., [real(dp) ::], [integer ::], .false., 0, 0, [1.0_dp], &
      [0.2_dp, 1.0_dp], 0.05_dp, 2, 1, 20, 1, tsissm_distribution_gaussian, &
      50, 1.0e-5_dp, 1013)
   call check(structural_constant_profile%info == 0 .and. &
      structural_constant_profile%successful == 1 .and. &
      all(shape(structural_constant_profile%parameter) == [1, 2]) .and. &
      structural_constant_profile%valid(1) .and. &
      all(ieee_is_finite(structural_constant_profile%mape)) .and. &
      all(ieee_is_finite(structural_constant_profile%crps)), &
      'joint structural constant simulation profile')

   structural_dynamic_profile = tsissm_profile_structural_dynamic(profile_regressors, &
      .false., .false., [real(dp) ::], [integer ::], .false., 0, 0, [1.0_dp], &
      [0.2_dp, 1.0_dp], [0.1_dp], [0.5_dp], 0.01_dp, [0.025_dp], [0.025_dp], &
      2, 1, 20, 1, tsissm_distribution_gaussian, 60, 1.0e-5_dp, 1021, &
      variance_initialization=tsissm_variance_initial_sample, variance_sample_size=4)
   call check(structural_dynamic_profile%info == 0 .and. &
      structural_dynamic_profile%successful == 1 .and. &
      all(shape(structural_dynamic_profile%parameter) == [1, 4]) .and. &
      structural_dynamic_profile%valid(1) .and. &
      all(ieee_is_finite(structural_dynamic_profile%mase)) .and. &
      all(ieee_is_finite(structural_dynamic_profile%crps)), &
      'joint structural dynamic simulation profile')

   dynamic_predicted = tsissm_predict_dynamic(transition, transition, transition, loading, &
      persistence, future_regressors, coefficients, initial, reshape([0.0_dp, 0.0_dp], [1, 2]), &
      [0.2_dp], [0.5_dp], [4.0_dp], [9.0_dp], 1.0_dp)
   call check(dynamic_predicted%info == 0 .and. &
      maxval(abs(dynamic_predicted%conditional_sd(1, :)**2 - [6.3_dp, 4.15_dp])) < 1.0e-12_dp .and. &
      maxval(abs(dynamic_predicted%observation)) == 0.0_dp, &
      'dynamic innovations prediction')

   structural = tsissm_model(.true., .true., [12.0_dp], [2], .false., [0.4_dp], &
      [0.3_dp], 0.2_dp, 0.1_dp, [0.05_dp], 0.9_dp)
   call check(structural%info == 0 .and. all(shape(structural%transition_base) == [8, 8]) .and. &
      structural%slope_index == 2 .and. structural%seasonal_start(1) == 3 .and. &
      structural%seasonal_end(1) == 6 .and. structural%ar_start == 7 .and. &
      structural%ma_start == 8 .and. &
      maxval(abs(structural%observation_loading - &
         [1.0_dp, 0.9_dp, 1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.4_dp, 0.3_dp])) < 1.0e-12_dp .and. &
      abs(structural%transition_base(1, 7) - 0.08_dp) < 1.0e-12_dp .and. &
      abs(structural%transition_base(7, 8) - 0.3_dp) < 1.0e-12_dp, &
      'trigonometric structural matrix construction')

   regular = tsissm_model(.false., .false., [4.0_dp], [1], .true., &
      [real(dp) ::], [real(dp) ::], 0.2_dp, 0.0_dp, [0.1_dp], 1.0_dp)
   call check(regular%info == 0 .and. all(shape(regular%transition_base) == [5, 5]) .and. &
      regular%transition_base(2, 5) == 1.0_dp .and. regular%transition_base(3, 2) == 1.0_dp .and. &
      regular%observation_loading(5) == 1.0_dp, 'regular seasonal matrix construction')

   gaussian_likelihood = tsissm_likelihood([1.0_dp, -1.0_dp], [1.0_dp, 1.0_dp], 1.0_dp, &
      tsissm_distribution_gaussian, 0.0_dp, 0.0_dp, 2)
   expected = log(2.0_dp*acos(-1.0_dp)) + 1.0_dp
   call check(gaussian_likelihood%info == 0 .and. &
      abs(gaussian_likelihood%objective - expected) < 1.0e-12_dp .and. &
      abs(gaussian_likelihood%log_likelihood + expected) < 1.0e-12_dp .and. &
      abs(gaussian_likelihood%sigma - 1.0_dp) < 1.0e-12_dp, &
      'Gaussian innovation likelihood')

   student_likelihood = tsissm_likelihood([0.0_dp], [2.0_dp], 1.0_dp, &
      tsissm_distribution_student, 0.0_dp, 5.0_dp, 1, conditional_sd=[1.0_dp])
   expected = log_gamma(3.0_dp) - log_gamma(2.5_dp) - 0.5_dp*log(5.0_dp*acos(-1.0_dp)) + &
      0.5_dp*log(5.0_dp/3.0_dp)
   call check(student_likelihood%info == 0 .and. &
      abs(student_likelihood%log_likelihood - expected) < 1.0e-12_dp, &
      'standardized Student innovation likelihood')

   jsu_likelihood = tsissm_likelihood([0.0_dp, 0.5_dp], [1.0_dp, 1.0_dp], 1.0_dp, &
      tsissm_distribution_johnson_su, 0.4_dp, 2.0_dp, 2, conditional_sd=[1.0_dp, 1.0_dp])
   call check(jsu_likelihood%info == 0 .and. ieee_is_finite(jsu_likelihood%log_likelihood) .and. &
      all(ieee_is_finite(jsu_likelihood%contribution)), 'Johnson SU innovation likelihood')

   call check(abs(tsissm_ar_constraint([0.5_dp], 0.01_dp) - 0.49_dp) < 1.0e-12_dp .and. &
      tsissm_ar_constraint([1.1_dp], 0.01_dp) < 0.0_dp .and. &
      abs(tsissm_ma_constraint([0.4_dp], 0.01_dp) - 0.59_dp) < 1.0e-12_dp, &
      'AR and MA companion constraints')
   call check(abs(tsissm_stability_constraint(reshape([0.9_dp], [1, 1]), &
      reshape([1.0_dp], [1, 1]), reshape([1.0_dp], [1, 1]), [1.0_dp], [0.1_dp], &
      0.01_dp) - 0.19_dp) < 1.0e-12_dp, 'discount transition constraint')

   constraints = tsissm_parameter_constraints(reshape([0.8_dp], [1, 1]), &
      reshape([1.0_dp], [1, 1]), reshape([1.0_dp], [1, 1]), [1.0_dp], [0.1_dp], &
      [real(dp) ::], [real(dp) ::], 0.1_dp, 0.2_dp, [0.1_dp], [0.8_dp], &
      tsissm_distribution_student, 0.0_dp, 5.0_dp, 0.5_dp, 0.0_dp, 1.0_dp, 0.01_dp)
   call check(constraints%info == 0 .and. constraints%feasible .and. &
      all(constraints%residual >= 0.0_dp), 'aggregate feasible constraints')
   bad_constraints = tsissm_parameter_constraints(reshape([0.8_dp], [1, 1]), &
      reshape([1.0_dp], [1, 1]), reshape([1.0_dp], [1, 1]), [1.0_dp], [0.1_dp], &
      [real(dp) ::], [real(dp) ::], 0.1_dp, 0.2_dp, [0.6_dp], [0.5_dp], &
      tsissm_distribution_student, 0.0_dp, 5.0_dp, 0.5_dp, 0.0_dp, 1.0_dp, 0.01_dp)
   call check(.not. bad_constraints%feasible .and. bad_constraints%residual(5) < 0.0_dp, &
      'aggregate infeasible constraints')

   initialized = tsissm_initialize_states([3.0_dp, 3.0_dp, 3.0_dp], transition, transition, &
      transition, loading, [0.0_dp], regressors, coefficients, [0.0_dp], [1], lambda=1.0_dp)
   call check(initialized%info == 0 .and. abs(initialized%state(1) - 2.0_dp) < 1.0e-12_dp, &
      'least-squares state initialization')

   estimated = tsissm_fit_constant([1.0_dp, 2.0_dp, 1.5_dp, 2.2_dp, 1.8_dp, 2.4_dp, &
      2.0_dp, 2.5_dp], transition, transition, transition, loading, [0.3_dp], fit_regressors, &
      [real(dp) ::], [0.0_dp], [integer ::], lambda=1.0_dp, &
      distribution=tsissm_distribution_gaussian, skew=0.0_dp, shape=0.0_dp, &
      max_iterations=80, tolerance=1.0e-6_dp)
   call check((estimated%info == 0 .or. estimated%info == 4) .and. &
      allocated(estimated%persistence) .and. estimated%persistence(1) >= 0.0_dp .and. &
      estimated%persistence(1) <= 1.0_dp .and. estimated%likelihood%info == 0 .and. &
      ieee_is_finite(estimated%likelihood%log_likelihood), 'constant ISSM parameter estimation')

   dynamic_estimated = tsissm_fit_dynamic([1.0_dp, 2.0_dp, 1.5_dp, 2.2_dp, 1.8_dp, &
      2.4_dp, 2.0_dp, 2.5_dp], transition, transition, transition, loading, [0.3_dp], &
      fit_regressors, [real(dp) ::], [0.1_dp], [0.7_dp], [0.0_dp], [integer ::], &
      lambda=1.0_dp, distribution=tsissm_distribution_gaussian, skew=0.0_dp, shape=0.0_dp, &
      max_iterations=80, tolerance=1.0e-6_dp)
   call check(dynamic_estimated%info >= 0 .and. dynamic_estimated%info <= 4 .and. &
      allocated(dynamic_estimated%arch) .and. allocated(dynamic_estimated%garch) .and. &
      dynamic_estimated%arch(1) >= 0.0_dp .and. dynamic_estimated%garch(1) >= 0.0_dp .and. &
      sum(dynamic_estimated%arch) + sum(dynamic_estimated%garch) < 1.0_dp .and. &
      dynamic_estimated%variance_intercept > 0.0_dp .and. &
      all(dynamic_estimated%filter%conditional_sd > 0.0_dp) .and. &
      ieee_is_finite(dynamic_estimated%likelihood%log_likelihood), &
      'dynamic ISSM parameter estimation')

   structural_estimated = tsissm_fit_structural_constant([1.0_dp, 2.0_dp, 1.5_dp, &
      2.2_dp, 1.8_dp, 2.4_dp, 2.0_dp, 2.5_dp], .false., .false., [real(dp) ::], &
      [integer ::], .false., 0, 0, fit_regressors, [0.0_dp], [0.3_dp, 1.0_dp], &
      distribution=tsissm_distribution_gaussian, max_iterations=80, tolerance=1.0e-6_dp)
   call check(structural_estimated%info >= 0 .and. structural_estimated%info <= 4 .and. &
      allocated(structural_estimated%parameters) .and. &
      all(shape(structural_estimated%model%transition_base) == [1, 1]) .and. &
      structural_estimated%model%persistence(1) >= 0.0_dp .and. &
      structural_estimated%model%persistence(1) <= 1.0_dp .and. &
      structural_estimated%lambda >= 0.0_dp .and. structural_estimated%lambda <= 1.5_dp .and. &
      ieee_is_finite(structural_estimated%likelihood%log_likelihood), &
      'joint structural parameter estimation')

   structural_fit_covariance = tsissm_covariance_structural([1.0_dp, 2.0_dp, 1.5_dp, &
      2.2_dp, 1.8_dp, 2.4_dp, 2.0_dp, 2.5_dp], .false., .false., [real(dp) ::], &
      [integer ::], .false., 0, 0, fit_regressors, [0.0_dp], &
      structural_estimated%parameters, tsissm_distribution_gaussian, &
      tsissm_covariance_qmle)
   call check(structural_fit_covariance%info == 0 .and. &
      all(shape(structural_fit_covariance%scores) == [8, 2]) .and. &
      all(shape(structural_fit_covariance%hessian) == [2, 2]) .and. &
      all(shape(structural_fit_covariance%covariance) == [2, 2]) .and. &
      all(ieee_is_finite(structural_fit_covariance%standard_error)), &
      'joint structural parameter covariance')

   dynamic_structural_estimated = tsissm_fit_structural_dynamic([1.0_dp, 2.0_dp, &
      1.5_dp, 2.2_dp, 1.8_dp, 2.4_dp, 2.0_dp, 2.5_dp], .false., .false., &
      [real(dp) ::], [integer ::], .false., 0, 0, fit_regressors, [0.0_dp], &
      [0.3_dp, 1.0_dp], [0.1_dp], [0.5_dp], distribution=tsissm_distribution_gaussian, &
      max_iterations=100, tolerance=1.0e-6_dp, &
      variance_initialization=tsissm_variance_initial_sample, variance_sample_size=3)
   call check(dynamic_structural_estimated%info >= 0 .and. &
      dynamic_structural_estimated%info <= 4 .and. &
      dynamic_structural_estimated%dynamic_variance .and. &
      all(shape(dynamic_structural_estimated%parameters) == [4]) .and. &
      size(dynamic_structural_estimated%arch) == 1 .and. &
      size(dynamic_structural_estimated%garch) == 1 .and. &
      dynamic_structural_estimated%arch(1) >= 0.0_dp .and. &
      dynamic_structural_estimated%garch(1) >= 0.0_dp .and. &
      sum(dynamic_structural_estimated%arch) + sum(dynamic_structural_estimated%garch) < 1.0_dp .and. &
      dynamic_structural_estimated%variance_intercept > 0.0_dp .and. &
      abs(dynamic_structural_estimated%initial_variance - &
      sum(dynamic_structural_estimated%filter%innovation(:3)**2)/3.0_dp) < 1.0e-10_dp .and. &
      all(dynamic_structural_estimated%filter%conditional_sd > 0.0_dp) .and. &
      ieee_is_finite(dynamic_structural_estimated%likelihood%log_likelihood), &
      'joint dynamic structural estimation')

   dynamic_structural_covariance = tsissm_covariance_structural_dynamic([1.0_dp, &
      2.0_dp, 1.5_dp, 2.2_dp, 1.8_dp, 2.4_dp, 2.0_dp, 2.5_dp], .false., .false., &
      [real(dp) ::], [integer ::], .false., 0, 0, fit_regressors, [0.0_dp], &
      [0.3_dp, 1.0_dp, 0.1_dp, 0.5_dp], 1, 1, tsissm_distribution_gaussian, &
      tsissm_covariance_qmle, variance_initialization=tsissm_variance_initial_sample, &
      variance_sample_size=3)
   call check(dynamic_structural_covariance%info == 0 .and. &
      all(shape(dynamic_structural_covariance%scores) == [8, 4]) .and. &
      all(shape(dynamic_structural_covariance%hessian) == [4, 4]) .and. &
      all(shape(dynamic_structural_covariance%covariance) == [4, 4]) .and. &
      all(ieee_is_finite(dynamic_structural_covariance%standard_error)), &
      'joint dynamic structural covariance')

   constant_moments = tsissm_moments_constant(transition, loading, persistence, [2.0_dp], &
      moment_regressors, coefficients, 4.0_dp, 1.0_dp, .true.)
   call check(constant_moments%info == 0 .and. &
      maxval(abs(constant_moments%transformed_mean - 2.0_dp)) < 1.0e-12_dp .and. &
      maxval(abs(constant_moments%transformed_variance - [4.0_dp, 5.0_dp, 6.0_dp])) < 1.0e-12_dp .and. &
      maxval(abs(constant_moments%mean - 3.0_dp)) < 1.0e-12_dp .and. &
      maxval(abs(constant_moments%variance - constant_moments%transformed_variance)) < 1.0e-12_dp, &
      'constant analytical forecast moments')

   dynamic_moments = tsissm_moments_dynamic(transition, loading, persistence, [2.0_dp], &
      moment_regressors, coefficients, [0.2_dp], [0.5_dp], 1.0_dp, [4.0_dp], [9.0_dp], &
      1.0_dp, .false.)
   call check(dynamic_moments%info == 0 .and. &
      abs(dynamic_moments%garch_variance(1) - 6.3_dp) < 1.0e-12_dp .and. &
      abs(dynamic_moments%garch_variance(2) - 5.41_dp) < 1.0e-12_dp .and. &
      abs(dynamic_moments%transformed_variance(2) - 6.7625_dp) < 1.0e-12_dp, &
      'dynamic analytical forecast moments')

   log_moments = tsissm_moments_constant(transition, loading, persistence, [0.0_dp], &
      moment_regressors, coefficients, 0.2_dp, 0.0_dp, .true.)
   call check(log_moments%info == 0 .and. &
      abs(log_moments%mean(1) - 1.1_dp) < 1.0e-12_dp .and. &
      abs(log_moments%variance(1) - (exp(0.2_dp) - 1.0_dp)) < 1.0e-12_dp, &
      'Box-Cox forecast bias correction')

   decomposition_regressors(:, 1) = [1.0_dp, 2.0_dp, 3.0_dp]
   structural_filter = tsissm_filter_constant([20.0_dp, 22.0_dp, 24.0_dp], &
      structural%transition_base, structural%transition_scale, structural%transition_parameter, &
      structural%observation_loading, structural%persistence, decomposition_regressors, [0.5_dp], &
      [10.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp, 7.0_dp, 8.0_dp], lambda=1.0_dp)
   decomposition = tsissm_decompose_filter(structural, structural_filter, &
      decomposition_regressors, [0.5_dp], [20.0_dp, 22.0_dp, 24.0_dp], 1.0_dp)
   call check(decomposition%info == 0 .and. size(decomposition%seasonal, 2) == 1 .and. &
      abs(decomposition%level(1) - 10.0_dp) < 1.0e-12_dp .and. &
      abs(decomposition%slope(1) - 1.8_dp) < 1.0e-12_dp .and. &
      abs(decomposition%seasonal(1, 1) - 7.0_dp) < 1.0e-12_dp .and. &
      abs(decomposition%autoregressive(1) - 2.8_dp) < 1.0e-12_dp .and. &
      abs(decomposition%moving_average(1) - 2.4_dp) < 1.0e-12_dp .and. &
      maxval(abs(decomposition%reconstructed_original - [20.0_dp, 22.0_dp, 24.0_dp])) < 1.0e-11_dp, &
      'filtered structural decomposition')

   prediction_regressors(:, 1) = [1.0_dp, 2.0_dp]
   structural_prediction = tsissm_predict_constant(structural%transition_base, &
      structural%transition_scale, structural%transition_parameter, structural%observation_loading, &
      structural%persistence, prediction_regressors, [0.5_dp], &
      [10.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp, 7.0_dp, 8.0_dp], &
      reshape([0.2_dp, -0.1_dp, 0.3_dp, 0.4_dp], [2, 2]))
   prediction_decomposition = tsissm_decompose_prediction(structural, structural_prediction, &
      prediction_regressors, [0.5_dp], 1.0_dp)
   call check(prediction_decomposition%info == 0 .and. &
      all(shape(prediction_decomposition%seasonal) == [2, 2, 1]) .and. &
      maxval(abs(prediction_decomposition%reconstructed_transformed - &
      structural_prediction%observation)) < 1.0e-12_dp, 'simulated structural decomposition')

   constant_forecast = tsissm_simulate_constant(transition, transition, transition, loading, &
      persistence, future_regressors, coefficients, initial, 0.2_dp, 0.0_dp, &
      tsissm_distribution_gaussian, 0.0_dp, 0.0_dp, 200, [0.1_dp, 0.5_dp, 0.9_dp], 731)
   repeated_forecast = tsissm_simulate_constant(transition, transition, transition, loading, &
      persistence, future_regressors, coefficients, initial, 0.2_dp, 0.0_dp, &
      tsissm_distribution_gaussian, 0.0_dp, 0.0_dp, 200, [0.1_dp, 0.5_dp, 0.9_dp], 731)
   call check(constant_forecast%info == 0 .and. &
      all(shape(constant_forecast%distribution) == [200, 2]) .and. &
      all(shape(constant_forecast%quantile) == [3, 2]) .and. &
      maxval(abs(constant_forecast%distribution - repeated_forecast%distribution)) == 0.0_dp .and. &
      all(constant_forecast%quantile(1, :) <= constant_forecast%quantile(2, :)) .and. &
      all(constant_forecast%quantile(2, :) <= constant_forecast%quantile(3, :)) .and. &
      all(constant_forecast%distribution > 0.0_dp), 'constant parametric forecast simulation')

   bootstrap_forecast = tsissm_simulate_constant(transition, transition, transition, loading, &
      persistence, future_regressors, coefficients, initial, 1.0_dp, 0.0_dp, &
      tsissm_distribution_gaussian, 0.0_dp, 0.0_dp, 10, seed=19, bootstrap_residuals=[0.1_dp])
   call check(bootstrap_forecast%info == 0 .and. &
      maxval(abs(bootstrap_forecast%transformed(:, 1) - 0.1_dp)) < 1.0e-12_dp .and. &
      maxval(abs(bootstrap_forecast%transformed(:, 2) - 0.15_dp)) < 1.0e-12_dp, &
      'constant residual bootstrap simulation')

   student_forecast = tsissm_simulate_constant(transition, transition, transition, loading, &
      persistence, future_regressors, coefficients, initial, 0.1_dp, 0.0_dp, &
      tsissm_distribution_student, 0.0_dp, 6.0_dp, 100, seed=41)
   jsu_forecast = tsissm_simulate_constant(transition, transition, transition, loading, &
      persistence, future_regressors, coefficients, initial, 0.1_dp, 0.0_dp, &
      tsissm_distribution_johnson_su, 0.4_dp, 2.0_dp, 100, seed=43)
   call check(student_forecast%info == 0 .and. jsu_forecast%info == 0 .and. &
      all(ieee_is_finite(student_forecast%distribution)) .and. &
      all(ieee_is_finite(jsu_forecast%distribution)), 'Student and Johnson SU forecast simulation')

   dynamic_forecast = tsissm_simulate_dynamic(transition, transition, transition, loading, &
      persistence, future_regressors, coefficients, initial, [0.2_dp], [0.5_dp], 1.0_dp, &
      [4.0_dp], [9.0_dp], 0.0_dp, tsissm_distribution_gaussian, 0.0_dp, 0.0_dp, 10, &
      seed=29, bootstrap_standardized_residuals=[0.0_dp])
   call check(dynamic_forecast%info == 0 .and. &
      maxval(abs(dynamic_forecast%distribution - 1.0_dp)) < 1.0e-12_dp .and. &
      abs(dynamic_forecast%prediction%conditional_sd(1, 1)**2 - 6.3_dp) < 1.0e-12_dp .and. &
      abs(dynamic_forecast%moments%garch_variance(1) - 6.3_dp) < 1.0e-12_dp, &
      'dynamic standardized-residual bootstrap simulation')

   diagnostics = tsissm_diagnose(estimated%filter, [1.0_dp, 2.0_dp, 1.5_dp, 2.2_dp, &
      1.8_dp, 2.4_dp, 2.0_dp, 2.5_dp], 1.0_dp, 2, 0, estimated%likelihood%sigma)
   call check(diagnostics%info == 0 .and. diagnostics%lag == 2 .and. &
      all(shape(diagnostics%residual_acf) == [3]) .and. &
      all(shape(diagnostics%squared_residual_acf) == [3]) .and. &
      maxval(abs(diagnostics%standardized_residual - &
      estimated%filter%innovation/estimated%likelihood%sigma)) < 1.0e-12_dp .and. &
      diagnostics%weighted_ljung_box_p_value >= 0.0_dp .and. &
      diagnostics%weighted_ljung_box_p_value <= 1.0_dp .and. &
      diagnostics%normality_p_value >= 0.0_dp .and. diagnostics%normality_p_value <= 1.0_dp .and. &
      diagnostics%arch%info == 0, 'residual and conditional-variance diagnostics')

   box_residuals = [1.0_dp, -0.5_dp, 0.25_dp, -0.75_dp, 0.6_dp, -0.2_dp, &
      0.4_dp, -0.1_dp]
   box_acf = acf_values(box_residuals, 2)
   box_pacf = pacf_values(box_residuals, 2)
   box_tests(1) = tsissm_weighted_box_test(box_residuals, 2, &
      tsissm_box_test_box_pierce)
   box_tests(2) = tsissm_weighted_box_test(box_residuals, 2, &
      tsissm_box_test_ljung_box)
   box_tests(3) = tsissm_weighted_box_test(box_residuals, 2, tsissm_box_test_monti)
   box_tests(4) = tsissm_weighted_box_test(box_residuals, 2, &
      tsissm_box_test_box_pierce, weighted=.false.)
   box_tests(5) = tsissm_weighted_box_test(box_residuals, 2, &
      tsissm_box_test_ljung_box, weighted=.false.)
   box_tests(6) = tsissm_weighted_box_test(box_residuals, 2, &
      tsissm_box_test_monti, weighted=.false.)
   call check(all(box_tests%info == 0) .and. &
      abs(box_tests(1)%statistic - 8.0_dp*(box_acf(2)**2 + &
      0.5_dp*box_acf(3)**2)) < 1.0e-12_dp .and. &
      abs(box_tests(2)%statistic - 80.0_dp*(box_acf(2)**2/7.0_dp + &
      0.5_dp*box_acf(3)**2/6.0_dp)) < 1.0e-12_dp .and. &
      abs(box_tests(3)%statistic - 80.0_dp*(box_pacf(1)**2/7.0_dp + &
      0.5_dp*box_pacf(2)**2/6.0_dp)) < 1.0e-12_dp .and. &
      abs(box_tests(4)%statistic - 8.0_dp*sum(box_acf(2:3)**2)) < 1.0e-12_dp .and. &
      abs(box_tests(5)%statistic - 80.0_dp*(box_acf(2)**2/7.0_dp + &
      box_acf(3)**2/6.0_dp)) < 1.0e-12_dp .and. &
      abs(box_tests(6)%statistic - 80.0_dp*(box_pacf(1)**2/7.0_dp + &
      box_pacf(2)**2/6.0_dp)) < 1.0e-12_dp .and. &
      all(box_tests%p_value >= 0.0_dp) .and. all(box_tests%p_value <= 1.0_dp) .and. &
      all(box_tests%degrees_of_freedom == 2), &
      'weighted and classical Box-Pierce Ljung-Box and Monti variants')

   transformed_test = tsissm_weighted_box_test(box_residuals, 2, &
      tsissm_box_test_ljung_box, 1, transform=tsissm_residual_squared)
   direct_test = tsissm_weighted_box_test(box_residuals**2, 2, &
      tsissm_box_test_ljung_box)
   call check(transformed_test%info == 0 .and. transformed_test%fitted_parameters == 0 .and. &
      abs(transformed_test%statistic - direct_test%statistic) < 1.0e-12_dp, &
      'weighted box test squared-residual transform')
   transformed_test = tsissm_weighted_box_test(box_residuals, 2, &
      tsissm_box_test_ljung_box, transform=tsissm_residual_log_squared)
   direct_test = tsissm_weighted_box_test(log(box_residuals**2), 2, &
      tsissm_box_test_ljung_box)
   call check(transformed_test%info == 0 .and. &
      abs(transformed_test%statistic - direct_test%statistic) < 1.0e-12_dp, &
      'weighted box test log-squared-residual transform')
   transformed_test = tsissm_weighted_box_test(box_residuals, 2, &
      tsissm_box_test_ljung_box, transform=tsissm_residual_absolute)
   direct_test = tsissm_weighted_box_test(abs(box_residuals), 2, &
      tsissm_box_test_ljung_box)
   call check(transformed_test%info == 0 .and. &
      abs(transformed_test%statistic - direct_test%statistic) < 1.0e-12_dp, &
      'weighted box test absolute-residual transform')

   structural_diagnostics = tsissm_diagnose_structural(estimated%filter, [1.0_dp, &
      2.0_dp, 1.5_dp, 2.2_dp, 1.8_dp, 2.4_dp, 2.0_dp, 2.5_dp], 1.0_dp, &
      [0.5_dp], [0.25_dp], 2, 0.05_dp, estimated%likelihood%sigma)
   call check(structural_diagnostics%info == 0 .and. &
      all(structural_diagnostics%weighted_lag == [1, 5, 5, 5]) .and. &
      all(structural_diagnostics%weighted_ljung_box_p_value >= 0.0_dp) .and. &
      all(structural_diagnostics%weighted_ljung_box_p_value <= 1.0_dp) .and. &
      size(structural_diagnostics%outlier_statistic) == 2 .and. &
      size(structural_diagnostics%outlier_critical) == 2 .and. &
      structural_diagnostics%outlier_count >= 0 .and. &
      structural_diagnostics%outlier_count <= 2 .and. &
      abs(structural_diagnostics%ar_inverse_root(1) - 0.5_dp) < 1.0e-10_dp .and. &
      abs(structural_diagnostics%ma_inverse_root(1) + 0.25_dp) < 1.0e-10_dp .and. &
      abs(structural_diagnostics%stability_modulus(1) - &
      abs(estimated%filter%discount_transition(1, 1))) < 1.0e-10_dp, &
      'structural roots outliers and weighted diagnostics')

   outlier_filter = tsissm_filter_constant([1.0_dp, 1.1_dp, 0.9_dp, 1.05_dp, &
      0.95_dp, 10.0_dp, 1.02_dp, 0.98_dp, 1.08_dp, 0.92_dp, 1.03_dp, 0.97_dp], &
      transition, transition, transition, loading, [0.0_dp], structural_backtest_regressors, &
      coefficients, initial, lambda=1.0_dp)
   outlier_diagnostics = tsissm_diagnose_structural(outlier_filter, [1.0_dp, &
      1.1_dp, 0.9_dp, 1.05_dp, 0.95_dp, 10.0_dp, 1.02_dp, 0.98_dp, 1.08_dp, &
      0.92_dp, 1.03_dp, 0.97_dp], 1.0_dp, [real(dp) ::], [real(dp) ::], 2, 0.05_dp)
   call check(outlier_diagnostics%info == 0 .and. &
      outlier_diagnostics%outlier_count >= 1 .and. &
      outlier_diagnostics%outlier_index(1) == 6, &
      'generalized ESD outlier identification')

   forecast_scores = tsissm_forecast_accuracy(exp([0.1_dp, 0.15_dp]), bootstrap_forecast, &
      [1.0_dp, 1.1_dp, 1.2_dp, 1.3_dp], 1, 0.1_dp)
   call check(forecast_scores%info == 0 .and. forecast_scores%observations == 2 .and. &
      forecast_scores%mse < 1.0e-24_dp .and. forecast_scores%point%mae < 1.0e-12_dp .and. &
      forecast_scores%smape < 1.0e-12_dp .and. forecast_scores%interval_coverage == 1.0_dp .and. &
      forecast_scores%interval_score < 1.0e-12_dp .and. &
      abs(forecast_scores%crps) < 1.0e-12_dp .and. &
      ieee_is_finite(forecast_scores%log_score), 'point and distributional forecast accuracy')

   constant_backtest = tsissm_backtest_constant([1.0_dp, 2.0_dp, 1.5_dp, 2.2_dp, &
      1.8_dp, 2.4_dp, 2.0_dp, 2.5_dp], fit_regressors, transition, transition, transition, &
      loading, [0.3_dp], [real(dp) ::], [0.0_dp], [integer ::], 5, 3, 30, 4, 0, &
      0.0_dp, tsissm_distribution_gaussian, 0.0_dp, 0.0_dp, 0.2_dp, 301, 20, 1.0e-6_dp)
   call check(constant_backtest%info == 0 .and. &
      all(shape(constant_backtest%forecast) == [3, 3]) .and. &
      all(count(constant_backtest%valid, dim=1) == [3, 2, 1]) .and. &
      all(constant_backtest%training_start == [2, 3, 4]) .and. &
      .not. any(constant_backtest%reestimated) .and. &
      abs(constant_backtest%rmse(2) - &
      sqrt(sum(constant_backtest%error(:2, 2)**2)/2.0_dp)) < 1.0e-12_dp .and. &
      all(constant_backtest%interval_coverage >= 0.0_dp) .and. &
      all(constant_backtest%interval_coverage <= 1.0_dp), &
      'fixed-window constant rolling backtest')

   refitted_backtest = tsissm_backtest_constant([1.0_dp, 2.0_dp, 1.5_dp, 2.2_dp, &
      1.8_dp, 2.4_dp, 2.0_dp, 2.5_dp], fit_regressors, transition, transition, transition, &
      loading, [0.3_dp], [real(dp) ::], [0.0_dp], [integer ::], 5, 1, 10, 0, 2, &
      0.0_dp, tsissm_distribution_gaussian, 0.0_dp, 0.0_dp, 0.2_dp, 401, 30, 1.0e-5_dp)
   call check(refitted_backtest%info == 0 .and. &
      all(refitted_backtest%reestimated .eqv. [.true., .false., .true.]) .and. &
      all(refitted_backtest%training_start == 1), 'expanding-window periodic re-estimation')

   dynamic_backtest = tsissm_backtest_dynamic([1.0_dp, 2.0_dp, 1.5_dp, 2.2_dp, &
      1.8_dp, 2.4_dp, 2.0_dp, 2.5_dp], fit_regressors, transition, transition, transition, &
      loading, [0.3_dp], [real(dp) ::], [0.2_dp], [0.5_dp], [0.0_dp], [integer ::], &
      [4.0_dp], [9.0_dp], 1.0_dp, 5, 2, 20, 0, 0, 0.0_dp, &
      tsissm_distribution_gaussian, 0.0_dp, 0.0_dp, 0.2_dp, 501, 20, 1.0e-6_dp)
   call check(dynamic_backtest%info == 0 .and. &
      all(shape(dynamic_backtest%forecast) == [3, 2]) .and. &
      all(count(dynamic_backtest%valid, dim=1) == [3, 2]) .and. &
      all(ieee_is_finite(dynamic_backtest%mean_crps)) .and. &
      all(ieee_is_finite(dynamic_backtest%mean_log_score)), &
      'dynamic rolling backtest distribution scores')

   structural_constant_backtest = tsissm_backtest_structural_constant([1.0_dp, &
      1.3_dp, 1.1_dp, 1.5_dp, 1.2_dp, 1.6_dp, 1.4_dp, 1.8_dp, 1.5_dp, 1.9_dp, &
      1.7_dp, 2.0_dp], structural_backtest_regressors, .false., .false., &
      [real(dp) ::], [integer ::], .false., 0, 0, [0.0_dp], [0.3_dp, 0.0_dp], &
      8, 2, 20, 0, 2, tsissm_distribution_gaussian, 0.1_dp, 60, 1.0e-5_dp, 1201)
   call check(structural_constant_backtest%info == 0 .and. &
      all(shape(structural_constant_backtest%parameters) == [4, 2]) .and. &
      all(structural_constant_backtest%reestimated .eqv. [.true., .false., .true., .false.]) .and. &
      all(structural_constant_backtest%valid(:, 1)) .and. &
      count(structural_constant_backtest%valid(:, 2)) == 3 .and. &
      all(ieee_is_finite(structural_constant_backtest%mean_crps)), &
      'joint structural constant rolling backtest')

   structural_dynamic_backtest = tsissm_backtest_structural_dynamic([1.0_dp, &
      1.3_dp, 1.1_dp, 1.5_dp, 1.2_dp, 1.6_dp, 1.4_dp, 1.8_dp, 1.5_dp, 1.9_dp, &
      1.7_dp, 2.0_dp], structural_backtest_regressors, .false., .false., &
      [real(dp) ::], [integer ::], .false., 0, 0, [0.0_dp], [0.3_dp, 0.0_dp], &
      [0.1_dp], [0.5_dp], 8, 2, 20, 0, 2, tsissm_distribution_gaussian, &
      0.1_dp, 60, 1.0e-5_dp, 1301, &
      variance_initialization=tsissm_variance_initial_sample, variance_sample_size=3)
   call check(structural_dynamic_backtest%info == 0 .and. &
      all(shape(structural_dynamic_backtest%parameters) == [4, 4]) .and. &
      all(structural_dynamic_backtest%reestimated .eqv. [.true., .false., .true., .false.]) .and. &
      all(structural_dynamic_backtest%valid(:, 1)) .and. &
      count(structural_dynamic_backtest%valid(:, 2)) == 3 .and. &
      all(ieee_is_finite(structural_dynamic_backtest%mean_log_score)), &
      'joint structural dynamic rolling backtest')

   selected_models = tsissm_select_models([1.0_dp, 2.0_dp, 1.5_dp, 2.2_dp, &
      1.8_dp, 2.4_dp, 2.0_dp, 2.5_dp], fit_regressors, [.false., .true.], &
      [.false., .true.], .true., [4.0_dp], reshape([1], [1, 1]), .false., [0], [0], &
      .true., .false., [1, 1], tsissm_distribution_gaussian, 1.0_dp, 0.0_dp, 0.0_dp, &
      2, tsissm_selection_bic, 40, 1.0e-5_dp)
   call check(selected_models%info == 0 .and. selected_models%evaluated == 6 .and. &
      selected_models%successful >= 1 .and. size(selected_models%top_index) <= 2 .and. &
      abs(sum(selected_models%weight) - 1.0_dp) < 1.0e-12_dp .and. &
      count(selected_models%candidates%seasonal) == 3 .and. &
      selected_models%candidates(selected_models%top_index(1))%score <= &
      selected_models%candidates(selected_models%top_index(size(selected_models%top_index)))%score .and. &
      all([(selected_models%candidates(selected_models%top_index(i))%successful, &
      i=1,size(selected_models%top_index))]), 'automatic structural model selection')

   reranked_models = tsissm_rerank_models(selected_models, &
      [(real(size(selected_models%candidates) - i, dp), i=1,size(selected_models%candidates))], 2)
   call check(reranked_models%info == 0 .and. &
      reranked_models%criterion == tsissm_selection_external .and. &
      abs(sum(reranked_models%weight) - 1.0_dp) < 1.0e-12_dp .and. &
      reranked_models%candidates(reranked_models%top_index(1))%score <= &
      reranked_models%candidates(reranked_models%top_index(size(reranked_models%top_index)))%score, &
      'external backtest-score model ranking')

   variance_selection = tsissm_select_models([1.0_dp, 2.0_dp, 1.5_dp, 2.2_dp, &
      1.8_dp, 2.4_dp, 2.0_dp, 2.5_dp], fit_regressors, [.false.], [.false.], .true., &
      [real(dp) ::], reshape([integer ::], [0, 0]), .false., [0], [0], .true., .true., &
      [1, 1], tsissm_distribution_gaussian, 1.0_dp, 0.0_dp, 0.0_dp, 1, &
      tsissm_selection_aic, 40, 1.0e-5_dp, &
      variance_initialization=tsissm_variance_initial_sample, variance_sample_size=3)
   call check(variance_selection%info == 0 .and. variance_selection%evaluated == 2 .and. &
      .not. variance_selection%candidates(1)%dynamic_variance .and. &
      variance_selection%candidates(2)%dynamic_variance .and. &
      variance_selection%candidates(2)%parameter_count == 4 .and. &
      size(variance_selection%candidates(2)%parameters) == 4 .and. &
      size(variance_selection%candidates(2)%arch) == 1 .and. &
      size(variance_selection%candidates(2)%garch) == 1 .and. &
      allocated(variance_selection%candidates(2)%filter%conditional_sd), &
      'constant and joint dynamic variance selection')

   compatible_selection = selected_models
   do i = 1, size(compatible_selection%top_index)
      compatible_selection%candidates(compatible_selection%top_index(i))%lambda = 0.0_dp
   end do
   compatible_selection%candidates(compatible_selection%top_index(1))%distribution = &
      tsissm_distribution_student
   compatible_selection%candidates(compatible_selection%top_index(1))%shape = 6.0_dp
   if (size(compatible_selection%top_index) > 1) then
      compatible_selection%candidates(compatible_selection%top_index(2))%distribution = &
         tsissm_distribution_johnson_su
      compatible_selection%candidates(compatible_selection%top_index(2))%skew = 0.4_dp
      compatible_selection%candidates(compatible_selection%top_index(2))%shape = 2.0_dp
   end if
   ensemble_forecast = tsissm_ensemble_forecast(compatible_selection, future_regressors, &
      100, [0.1_dp, 0.9_dp], seed=811)
   ensemble_components = ensemble_forecast%trend + ensemble_forecast%seasonal + &
      ensemble_forecast%arma + ensemble_forecast%regression + ensemble_forecast%irregular
   allocate(ensemble_expected(100, 2))
   allocate(ensemble_distribution_expected(100, 2))
   ensemble_expected = 0.0_dp
   ensemble_distribution_expected = 0.0_dp
   do i = 1, size(compatible_selection%top_index)
      ensemble_expected = ensemble_expected + ensemble_forecast%weight(i)* &
         ensemble_forecast%model_forecast(i)%prediction%observation
      ensemble_distribution_expected = ensemble_distribution_expected + ensemble_forecast%weight(i)* &
         ensemble_forecast%model_forecast(i)%distribution
   end do
   call check(ensemble_forecast%info == 0 .and. ensemble_forecast%has_decomposition .and. &
      all(shape(ensemble_forecast%distribution) == [100, 2]) .and. &
      all(shape(ensemble_forecast%correlation) == &
      [size(compatible_selection%top_index), size(compatible_selection%top_index)]) .and. &
      all([(abs(ensemble_forecast%correlation(i, i) - 1.0_dp) < 1.0e-12_dp, &
      i=1,size(compatible_selection%top_index))]) .and. &
      all(ensemble_forecast%quantile(1, :) <= ensemble_forecast%quantile(2, :)) .and. &
      all(ieee_is_finite(ensemble_forecast%distribution)) .and. &
      maxval(abs(ensemble_forecast%distribution - ensemble_distribution_expected)) < 1.0e-12_dp .and. &
      maxval(abs(ensemble_components - ensemble_expected)) < 1.0e-10_dp, &
      'Gaussian-copula weighted forecast ensemble')
   print '(a)', 'All tsissm_mod tests passed.'

contains

   subroutine check(ok, name)
      ! Stop the test program when a named assertion fails.
      logical, intent(in) :: ok
      character(len=*), intent(in) :: name

      if (.not. ok) then
         print '(a)', 'FAILED: '//name
         error stop 1
      end if
   end subroutine check
end program test_tsissm
