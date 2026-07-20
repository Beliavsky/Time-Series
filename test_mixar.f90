! SPDX-License-Identifier: GPL-2.0-or-later
! SPDX-FileComment: Regression tests for the R mixAR translation.
program test_mixar
   use kind_mod, only: dp
   use mixar_mod
   use random_mod, only: set_random_seed
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   type(mixar_model_t) :: model, initial
   type(mixar_simulation_t) :: simulation
   type(mixar_filter_t) :: filtered
   type(mixar_fit_t) :: fit
   type(mixar_fit_t) :: general_fit
   type(mixar_model_t) :: general_model, general_initial
   type(mixar_simulation_t) :: general_simulation
   type(mixar_model_t) :: single_model
   type(mixar_exact_forecast_t) :: exact_forecast
   type(mixar_forecast_sample_t) :: sampled_forecast
   type(mixar_model_t) :: regression_mar, regression_initial_mar
   type(mixar_regression_model_t) :: regression_model, regression_initial
   type(mixar_regression_model_t) :: component_regression_model
   type(mixar_regression_fit_t) :: regression_fit, component_regression_fit
   type(mixar_regression_simulation_t) :: regression_simulation
   type(mixar_regression_simulation_t) :: component_regression_simulation
   type(mixar_seasonal_model_t) :: seasonal_model, seasonal_initial
   type(mixar_seasonal_model_t) :: seasonal_student_model
   type(mixar_seasonal_fit_t) :: seasonal_fit, seasonal_student_fit
   type(mixar_seasonal_simulation_t) :: seasonal_simulation
   type(mixar_seasonal_simulation_t) :: seasonal_student_simulation
   type(mixar_inference_t) :: inference, seasonal_inference, regression_inference
   type(mixar_diagnostics_t) :: diagnostics
   type(mixar_bic_selection_t) :: bic_selection
   type(mixar_var_model_t) :: var_model, var_initial
   type(mixar_var_simulation_t) :: var_simulation
   type(mixar_var_filter_t) :: var_filtered
   type(mixar_var_fit_t) :: var_fit
   type(mixar_var_forecast_t) :: var_forecast
   type(mixar_var_diagnostics_t) :: var_diagnostics
   type(mixar_bayesian_random_t) :: bayesian_random
   type(mixar_bayesian_draws_t) :: bayesian_draws, bayesian_repeat
   type(mixar_bayesian_draws_t) :: relabeled_draws
   type(mixar_marginal_likelihood_t) :: marginal_likelihood
   type(mixar_order_random_t) :: order_random
   type(mixar_order_selection_t) :: order_selection, order_repeat
   type(mixar_initialization_t) :: random_initialization
   type(mixar_multistart_t) :: multistart
   type(mixar_moment_t) :: normal_fourth, student_fourth, undefined_fourth
   type(mixar_moment_t) :: innovation_second, conditional_third
   type(mixar_moment_t) :: conditional_variance, conditional_kurtosis
   type(mixar_stationary_moments_t) :: stationary_moments
   type(mixar_model_t) :: bayesian_model
   type(mixar_simulation_t) :: bayesian_simulation
   type(mixar_model_t) :: candidate_models(2)
   real(dp) :: ar(2, 2), initial_ar(2, 2), density, cdf
   real(dp), allocatable :: uniforms(:), normals(:), locations(:, :)
   real(dp), allocatable :: forecast_uniforms(:, :), forecast_innovations(:, :)
   real(dp) :: single_ar(1, 1), expected_mean, expected_variance, forecast_quantile
   real(dp) :: regression_ar(1, 1), component_ar(1, 2)
   real(dp) :: seasonal_ar(1, 1), seasonal_initial_ar(1, 1)
   real(dp), allocatable :: regression_x(:, :), regression_history_x(:, :)
   real(dp), allocatable :: component_x(:, :), component_uniforms(:)
   real(dp), allocatable :: component_innovations(:)
   real(dp) :: var_ar(2, 2, 2, 2), var_initial_ar(2, 2, 2, 2)
   real(dp) :: var_covariance(2, 2, 2), var_intercept(2, 2)
   real(dp), allocatable :: var_uniforms(:), var_normals(:, :)
   real(dp), allocatable :: var_forecast_uniforms(:, :)
   real(dp), allocatable :: var_forecast_normals(:, :, :)
   integer, allocatable :: subsample_indices(:, :, :)
   integer, allocatable :: repeated_indices(:, :, :)
   integer :: i, j, order_status
   real(dp) :: birth_probability, death_probability

   ar(:, 1) = [0.55_dp, 0.0_dp]
   ar(:, 2) = [-0.15_dp, 0.25_dp]
   model = mixar_model([0.65_dp, 0.35_dp], [1, 2], &
      [0.30_dp, -0.45_dp], [0.22_dp, 0.32_dp], ar)
   call check(model%info == 0 .and. mixar_is_stable(model), &
      'valid stable ragged Gaussian MAR model')

   allocate(uniforms(1300), normals(1300))
   do i = 1, size(uniforms)
      uniforms(i) = modulo(0.6180339887498949_dp*real(i, dp), 1.0_dp)
      normals(i) = sin(0.731_dp*real(i, dp)) + 0.55_dp*cos(1.117_dp*real(i, dp))
   end do
   simulation = mixar_simulate_from_draws(model, 1100, [0.0_dp, 0.0_dp], &
      uniforms, normals, 200)
   call check(simulation%info == 0 .and. size(simulation%series) == 1100 .and. &
      all(simulation%regime >= 1) .and. all(simulation%regime <= 2), &
      'supplied-draw Gaussian MAR simulation')

   locations = mixar_component_locations(model, simulation%series)
   call check(all(shape(locations) == [1098, 2]) .and. &
      abs(locations(1, 1) - (model%shift(1) + &
      model%ar(1, 1)*simulation%series(2))) < 1.0e-13_dp, &
      'ragged component conditional locations')
   density = mixar_conditional_density(model, simulation%series(3), &
      simulation%series(1:2))
   cdf = mixar_conditional_cdf(model, simulation%series(3), simulation%series(1:2))
   call check(density > 0.0_dp .and. cdf > 0.0_dp .and. cdf < 1.0_dp, &
      'conditional Gaussian mixture density and cdf')

   filtered = mixar_filter(model, simulation%series)
   call check(filtered%info == 0 .and. ieee_is_finite(filtered%log_likelihood) .and. &
      maxval(abs(sum(filtered%responsibility, dim=2) - 1.0_dp)) < 1.0e-13_dp .and. &
      all(filtered%variance > 0.0_dp), 'filter and posterior component probabilities')

   initial_ar(:, 1) = [0.40_dp, 0.0_dp]
   initial_ar(:, 2) = [-0.05_dp, 0.10_dp]
   initial = mixar_model([0.55_dp, 0.45_dp], [1, 2], &
      [0.10_dp, -0.20_dp], [0.35_dp, 0.40_dp], initial_ar)
   fit = mixar_fit(simulation%series, initial, max_iterations=250, tolerance=1.0e-9_dp)
   call check(fit%info == 0 .and. fit%converged .and. &
      ieee_is_finite(fit%log_likelihood) .and. &
      fit%log_likelihood > mixar_log_likelihood(initial, simulation%series) .and. &
      maxval(abs(fit%model%probability - model%probability)) < 0.15_dp .and. &
      maxval(abs(fit%model%shift - model%shift)) < 0.20_dp, &
      'Gaussian MAR fixed-point EM estimation')

   call check(abs(mixar_standard_density(0.0_dp, 0.0_dp) - &
      1.0_dp/sqrt(2.0_dp*acos(-1.0_dp))) < 1.0e-14_dp .and. &
      abs(mixar_standard_cdf(-0.7_dp, 5.0_dp) + &
      mixar_standard_cdf(0.7_dp, 5.0_dp) - 1.0_dp) < 1.0e-13_dp, &
      'standard normal and unit-variance Student-t distributions')

   general_model = mixar_model([0.60_dp, 0.40_dp], [1, 2], &
      [0.25_dp, -0.35_dp], [0.20_dp, 0.28_dp], ar, [0.0_dp, 6.0_dp])
   call set_random_seed(9271)
   general_simulation = mixar_simulate(general_model, 1400, [0.0_dp, 0.0_dp], 250)
   call check(general_simulation%info == 0 .and. &
      all(ieee_is_finite(general_simulation%series)), &
      'mixed Gaussian and Student-t shared-RNG simulation')
   filtered = mixar_filter(general_model, general_simulation%series)
   call check(filtered%info == 0 .and. &
      maxval(abs(sum(filtered%responsibility, dim=2) - 1.0_dp)) < 1.0e-13_dp, &
      'mixed-family filtering and responsibilities')

   general_initial = mixar_model([0.52_dp, 0.48_dp], [1, 2], &
      [0.10_dp, -0.15_dp], [0.32_dp, 0.38_dp], initial_ar, [0.0_dp, 10.0_dp])
   general_fit = mixar_general_fit(general_simulation%series, general_initial, &
      max_iterations=250, tolerance=1.0e-8_dp, &
      estimate_degrees=[.false., .true.])
   call check(general_fit%info == 0 .and. general_fit%converged .and. &
      general_fit%log_likelihood >= &
      mixar_log_likelihood(general_initial, general_simulation%series) .and. &
      general_fit%model%degrees_of_freedom(1) == 0.0_dp .and. &
      general_fit%model%degrees_of_freedom(2) > 2.05_dp .and. &
      general_fit%model%degrees_of_freedom(2) < 30.0_dp .and. &
      maxval(abs(general_fit%model%probability - &
      general_model%probability)) < 0.16_dp, &
      'generalized EM with selective Student-t degree estimation')

   single_ar(1, 1) = 0.6_dp
   single_model = mixar_model([1.0_dp], [1], [0.4_dp], [0.3_dp], single_ar)
   exact_forecast = mixar_exact_forecast(single_model, [1.25_dp], 3, 0.90_dp)
   expected_mean = 0.4_dp*(1.0_dp + 0.6_dp) + 0.6_dp**2*1.25_dp
   expected_variance = 0.3_dp**2*(1.0_dp + 0.6_dp**2)
   call check(exact_forecast%info == 0 .and. &
      abs(exact_forecast%mean(2) - expected_mean) < 1.0e-13_dp .and. &
      abs(exact_forecast%variance(2) - expected_variance) < 1.0e-13_dp .and. &
      exact_forecast%lower(2) < exact_forecast%mean(2) .and. &
      exact_forecast%upper(2) > exact_forecast%mean(2), &
      'exact AR(1) shift and innovation propagation')

   exact_forecast = mixar_exact_forecast(model, simulation%series(1099:1100), 3)
   forecast_quantile = mixar_predictive_quantile(exact_forecast%distribution(3), 0.8_dp)
   call check(exact_forecast%info == 0 .and. &
      size(exact_forecast%distribution(1)%probability) == 2 .and. &
      size(exact_forecast%distribution(3)%probability) == 8 .and. &
      abs(sum(exact_forecast%distribution(3)%probability) - 1.0_dp) < 1.0e-13_dp .and. &
      mixar_predictive_density(exact_forecast%distribution(3), forecast_quantile) > 0.0_dp .and. &
      abs(mixar_predictive_cdf(exact_forecast%distribution(3), &
      forecast_quantile) - 0.8_dp) < 1.0e-11_dp, &
      'exact multi-component regime-path forecast distribution')

   allocate(forecast_uniforms(4, 600), forecast_innovations(4, 600))
   do j = 1, size(forecast_uniforms, 2)
      do i = 1, size(forecast_uniforms, 1)
         forecast_uniforms(i, j) = modulo(0.4142135623730950_dp* &
            real(i + 4*(j - 1), dp), 1.0_dp)
         forecast_innovations(i, j) = sin(0.319_dp*real(i + 4*(j - 1), dp)) + &
            0.45_dp*cos(0.811_dp*real(i + 4*(j - 1), dp))
      end do
   end do
   sampled_forecast = mixar_forecast_from_draws(general_model, &
      general_simulation%series(1399:1400), forecast_uniforms, &
      forecast_innovations, 0.90_dp)
   call check(sampled_forecast%info == 0 .and. &
      all(shape(sampled_forecast%paths) == [4, 600]) .and. &
      all(ieee_is_finite(sampled_forecast%mean)) .and. &
      all(sampled_forecast%variance > 0.0_dp) .and. &
      all(sampled_forecast%lower < sampled_forecast%upper), &
      'supplied-draw mixed-family forecast paths and summaries')

   regression_ar(1, 1) = 0.45_dp
   regression_mar = mixar_model([1.0_dp], [1], [0.0_dp], [0.18_dp], &
      regression_ar)
   regression_model = mixar_regression_model(regression_mar, &
      reshape([1.2_dp, 0.7_dp], [2, 1]))
   allocate(regression_history_x(1, 2), regression_x(900, 2))
   regression_history_x(1, :) = [1.0_dp, 0.0_dp]
   do i = 1, size(regression_x, 1)
      regression_x(i, :) = [1.0_dp, sin(0.031_dp*real(i, dp))]
   end do
   call set_random_seed(6517)
   regression_simulation = mixar_regression_simulate(regression_model, &
      [1.2_dp], regression_history_x, regression_x)
   call check(regression_simulation%info == 0 .and. &
      all(ieee_is_finite(regression_simulation%series)), &
      'shared-coefficient MAR regression simulation')
   regression_ar(1, 1) = 0.20_dp
   regression_initial_mar = mixar_model([1.0_dp], [1], [0.0_dp], &
      [0.30_dp], regression_ar)
   regression_initial = mixar_regression_model(regression_initial_mar, &
      reshape([0.8_dp, 0.3_dp], [2, 1]))
   regression_fit = mixar_regression_fit(regression_simulation%series, &
      regression_x, regression_initial, max_iterations=200, tolerance=1.0e-8_dp)
   call check(regression_fit%info == 0 .and. regression_fit%converged .and. &
      regression_fit%log_likelihood > mixar_regression_log_likelihood( &
      regression_initial, regression_simulation%series, regression_x) .and. &
      maxval(abs(regression_fit%model%coefficient(:, 1) - &
      [1.2_dp, 0.7_dp])) < 0.10_dp .and. &
      abs(regression_fit%model%mar%ar(1, 1) - 0.45_dp) < 0.10_dp, &
      'shared-coefficient MAR regression estimation')

   exact_forecast = mixar_regression_exact_forecast(regression_model, &
      regression_simulation%series(900:900), regression_x(900:900, :), &
      reshape([1.0_dp, sin(0.031_dp*901.0_dp), 1.0_dp, &
      sin(0.031_dp*902.0_dp)], [2, 2], order=[2, 1]))
   call check(exact_forecast%info == 0 .and. &
      size(exact_forecast%distribution(2)%probability) == 1 .and. &
      exact_forecast%variance(2) > exact_forecast%variance(1), &
      'exact MAR regression forecast with future regressors')

   component_ar = 0.0_dp
   regression_mar = mixar_model([0.58_dp, 0.42_dp], [1, 1], &
      [0.0_dp, 0.0_dp], [0.16_dp, 0.22_dp], component_ar, [0.0_dp, 7.0_dp])
   component_regression_model = mixar_regression_model(regression_mar, &
      reshape([0.4_dp, 0.8_dp, 1.6_dp, -0.5_dp], [2, 2]), .true.)
   allocate(component_x(1200, 2), component_uniforms(1200))
   allocate(component_innovations(1200))
   do i = 1, size(component_x, 1)
      component_x(i, :) = [1.0_dp, sin(0.023_dp*real(i, dp))]
      component_uniforms(i) = modulo(0.7548776662466927_dp*real(i, dp), 1.0_dp)
      component_innovations(i) = sin(0.619_dp*real(i, dp)) + &
         0.4_dp*cos(1.037_dp*real(i, dp))
   end do
   component_regression_simulation = mixar_regression_simulate_from_draws( &
      component_regression_model, [0.4_dp], reshape([1.0_dp, 0.0_dp], [1, 2]), &
      component_x, component_uniforms, component_innovations)
   filtered = mixar_regression_filter(component_regression_model, &
      component_regression_simulation%series, component_x)
   call check(component_regression_simulation%info == 0 .and. filtered%info == 0 .and. &
      maxval(abs(sum(filtered%responsibility, dim=2) - 1.0_dp)) < 1.0e-13_dp, &
      'component-specific mixed-family MAR regression filtering')

   regression_initial_mar = mixar_model([0.52_dp, 0.48_dp], [1, 1], &
      [0.0_dp, 0.0_dp], [0.25_dp, 0.28_dp], component_ar, [0.0_dp, 7.0_dp])
   regression_initial = mixar_regression_model(regression_initial_mar, &
      reshape([0.2_dp, 0.5_dp, 1.3_dp, -0.2_dp], [2, 2]), .true.)
   component_regression_fit = mixar_regression_fit( &
      component_regression_simulation%series, component_x, regression_initial, &
      max_iterations=250, tolerance=1.0e-8_dp)
   call check(component_regression_fit%info == 0 .and. &
      component_regression_fit%converged .and. &
      component_regression_fit%log_likelihood >= mixar_regression_log_likelihood( &
      regression_initial, component_regression_simulation%series, component_x) .and. &
      maxval(abs(component_regression_fit%model%coefficient - &
      component_regression_model%coefficient)) < 0.25_dp, &
      'component-specific mixed-family MAR regression estimation')

   sampled_forecast = mixar_regression_forecast_from_draws( &
      component_regression_model, component_regression_simulation%series(1200:1200), &
      component_x(1200:1200, :), component_x(1:4, :), forecast_uniforms, &
      forecast_innovations, 0.90_dp)
   call check(sampled_forecast%info == 0 .and. &
      all(shape(sampled_forecast%paths) == [4, 600]) .and. &
      all(sampled_forecast%lower < sampled_forecast%upper), &
      'mixed-family MAR regression forecast with future regressors')

   regression_ar(1, 1) = 0.30_dp
   seasonal_ar(1, 1) = 0.42_dp
   regression_mar = mixar_model([1.0_dp], [1], [0.15_dp], [0.20_dp], &
      regression_ar)
   seasonal_model = mixar_seasonal_model(regression_mar, 4, [1], seasonal_ar)
   regression_initial_mar = mixar_seasonal_expanded_model(seasonal_model)
   call check(seasonal_model%info == 0 .and. mixar_seasonal_is_stable(seasonal_model) .and. &
      abs(regression_initial_mar%ar(4, 1) - 0.42_dp) < 1.0e-14_dp, &
      'seasonal MAR construction, lag expansion, and stability')
   call set_random_seed(8129)
   seasonal_simulation = mixar_seasonal_simulate(seasonal_model, 1300, &
      [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], 250)
   filtered = mixar_seasonal_filter(seasonal_model, seasonal_simulation%series)
   call check(seasonal_simulation%info == 0 .and. filtered%info == 0 .and. &
      ieee_is_finite(mixar_seasonal_log_likelihood(seasonal_model, &
      seasonal_simulation%series)), 'seasonal MAR simulation and filtering')

   regression_ar(1, 1) = 0.12_dp
   seasonal_initial_ar(1, 1) = 0.20_dp
   regression_initial_mar = mixar_model([1.0_dp], [1], [0.15_dp], [0.32_dp], &
      regression_ar)
   seasonal_initial = mixar_seasonal_model(regression_initial_mar, 4, [1], &
      seasonal_initial_ar)
   seasonal_fit = mixar_seasonal_fit(seasonal_simulation%series, seasonal_initial, &
      max_iterations=250, tolerance=1.0e-8_dp)
   call check(seasonal_fit%info == 0 .and. seasonal_fit%converged .and. &
      seasonal_fit%log_likelihood > mixar_seasonal_log_likelihood( &
      seasonal_initial, seasonal_simulation%series) .and. &
      abs(seasonal_fit%model%mar%ar(1, 1) - 0.30_dp) < 0.10_dp .and. &
      abs(seasonal_fit%model%seasonal_ar(1, 1) - 0.42_dp) < 0.10_dp, &
      'Gaussian seasonal MAR generalized EM estimation')

   exact_forecast = mixar_seasonal_exact_forecast(seasonal_model, &
      seasonal_simulation%series(1297:1300), 5, 0.90_dp)
   sampled_forecast = mixar_seasonal_forecast_from_draws(seasonal_model, &
      seasonal_simulation%series(1297:1300), forecast_uniforms, &
      forecast_innovations, 0.90_dp)
   call check(exact_forecast%info == 0 .and. sampled_forecast%info == 0 .and. &
      size(exact_forecast%mean) == 5 .and. all(shape(sampled_forecast%paths) == [4, 600]), &
      'exact and supplied-draw seasonal MAR forecasts')

   regression_mar = mixar_model([1.0_dp], [1], [0.0_dp], [0.22_dp], &
      reshape([0.25_dp], [1, 1]), [6.0_dp])
   seasonal_student_model = mixar_seasonal_model(regression_mar, 4, [1], &
      reshape([0.35_dp], [1, 1]))
   call set_random_seed(9137)
   seasonal_student_simulation = mixar_seasonal_simulate(seasonal_student_model, &
      1000, [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], 200)
   regression_initial_mar = mixar_model([1.0_dp], [1], [0.0_dp], [0.30_dp], &
      reshape([0.10_dp], [1, 1]), [10.0_dp])
   seasonal_initial = mixar_seasonal_model(regression_initial_mar, 4, [1], &
      reshape([0.15_dp], [1, 1]))
   seasonal_student_fit = mixar_seasonal_fit(seasonal_student_simulation%series, &
      seasonal_initial, max_iterations=250, tolerance=1.0e-8_dp, &
      estimate_degrees=[.true.])
   call check(seasonal_student_fit%info == 0 .and. seasonal_student_fit%converged .and. &
      seasonal_student_fit%model%mar%degrees_of_freedom(1) > 2.05_dp .and. &
      seasonal_student_fit%model%mar%degrees_of_freedom(1) < 30.0_dp, &
      'Student-t seasonal MAR fitting with estimated degrees')

   inference = mixar_observed_inference(fit%model, simulation%series)
   call check(inference%info == 0 .and. size(inference%estimate) > 0 .and. &
      size(inference%probability_standard_error) == 2 .and. &
      all(ieee_is_finite(inference%standard_error)) .and. &
      all(inference%standard_error >= 0.0_dp), &
      'ordinary MAR observed-information inference')

   seasonal_inference = mixar_seasonal_observed_inference(seasonal_fit%model, &
      seasonal_simulation%series, estimate_shift=.false.)
   call check(seasonal_inference%info == 0 .and. &
      any(index(seasonal_inference%parameter, 'seasonal_ar_') == 1) .and. &
      all(ieee_is_finite(seasonal_inference%standard_error)), &
      'seasonal MAR observed-information inference')

   regression_inference = mixar_regression_observed_inference(regression_fit%model, &
      regression_simulation%series, regression_x, estimate_shift=.false.)
   call check(regression_inference%info == 0 .and. &
      count(index(regression_inference%parameter, 'regression_') == 1) == 2 .and. &
      all(ieee_is_finite(regression_inference%standard_error)), &
      'MAR regression observed-information inference')

   diagnostics = mixar_diagnose(fit%model, simulation%series, 12)
   call check(diagnostics%info == 0 .and. &
      size(diagnostics%residual_acf) == 13 .and. &
      all(diagnostics%uniform_residual > 0.0_dp) .and. &
      all(diagnostics%uniform_residual < 1.0_dp) .and. &
      all(ieee_is_finite(diagnostics%quantile_residual)) .and. &
      all(diagnostics%classification >= 1) .and. &
      all(diagnostics%classification <= 2) .and. &
      diagnostics%residual_test%info == 0 .and. &
      diagnostics%squared_residual_test%info == 0, &
      'residual, quantile, portmanteau, and classification diagnostics')

   candidate_models(1) = fit%model
   candidate_models(2) = initial
   bic_selection = mixar_select_bic(candidate_models, simulation%series)
   call check(bic_selection%info == 0 .and. bic_selection%selected == 1 .and. &
      bic_selection%bic(1) < bic_selection%bic(2), &
      'conditional BIC model comparison')

   var_ar = 0.0_dp
   var_ar(:, :, 1, 1) = reshape([0.42_dp, 0.08_dp, -0.12_dp, 0.30_dp], [2, 2])
   var_ar(:, :, 1, 2) = reshape([-0.18_dp, 0.05_dp, 0.10_dp, 0.22_dp], [2, 2])
   var_ar(:, :, 2, 2) = reshape([0.08_dp, 0.02_dp, -0.03_dp, 0.06_dp], [2, 2])
   var_intercept = reshape([0.20_dp, -0.10_dp, -0.25_dp, 0.18_dp], [2, 2])
   var_covariance(:, :, 1) = reshape([0.18_dp, 0.04_dp, 0.04_dp, 0.14_dp], [2, 2])
   var_covariance(:, :, 2) = reshape([0.25_dp, -0.03_dp, -0.03_dp, 0.20_dp], [2, 2])
   var_model = mixar_var_model([0.62_dp, 0.38_dp], [1, 2], var_intercept, &
      var_ar, var_covariance)
   call check(var_model%info == 0 .and. mixar_var_is_stable(var_model), &
      'valid stable ragged Gaussian mixture VAR model')

   allocate(var_uniforms(1400), var_normals(1400, 2))
   do i = 1, size(var_uniforms)
      var_uniforms(i) = modulo(0.414213562373095_dp*real(i, dp), 1.0_dp)
      var_normals(i, 1) = sin(0.517_dp*real(i, dp)) + 0.3_dp*cos(1.113_dp*real(i, dp))
      var_normals(i, 2) = cos(0.739_dp*real(i, dp)) - 0.2_dp*sin(1.271_dp*real(i, dp))
   end do
   var_simulation = mixar_var_simulate_from_draws(var_model, 1200, &
      reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], [2, 2]), &
      var_uniforms, var_normals, 200)
   var_filtered = mixar_var_filter(var_model, var_simulation%series)
   call check(var_simulation%info == 0 .and. var_filtered%info == 0 .and. &
      all(shape(var_simulation%series) == [1200, 2]) .and. &
      all(shape(var_filtered%location) == [1198, 2, 2]) .and. &
      maxval(abs(sum(var_filtered%responsibility, dim=2) - 1.0_dp)) < 1.0e-12_dp, &
      'mixture VAR simulation, likelihood, and posterior probabilities')

   var_initial_ar = 0.7_dp*var_ar
   var_initial = mixar_var_model([0.55_dp, 0.45_dp], [1, 2], &
      0.6_dp*var_intercept, var_initial_ar, 1.3_dp*var_covariance)
   var_fit = mixar_var_fit(var_simulation%series, var_initial, &
      max_iterations=250, tolerance=1.0e-8_dp)
   call check(var_fit%info == 0 .and. var_fit%converged .and. &
      var_fit%log_likelihood > mixar_var_log_likelihood( &
      var_initial, var_simulation%series) .and. &
      abs(sum(var_fit%model%probability) - 1.0_dp) < 1.0e-12_dp, &
      'mixture VAR expectation-maximization estimation')

   allocate(var_forecast_uniforms(5, 300), var_forecast_normals(5, 2, 300))
   do j = 1, 300
      do i = 1, 5
         var_forecast_uniforms(i, j) = modulo(0.2718281828_dp*real(i + 5*j, dp), 1.0_dp)
         var_forecast_normals(i, 1, j) = sin(0.37_dp*real(i + 5*j, dp))
         var_forecast_normals(i, 2, j) = cos(0.43_dp*real(i + 5*j, dp))
      end do
   end do
   var_forecast = mixar_var_forecast_from_draws(var_model, &
      var_simulation%series(1199:1200, :), var_forecast_uniforms, &
      var_forecast_normals, 0.90_dp)
   var_diagnostics = mixar_var_diagnose(var_model, var_simulation%series, 10)
   call check(var_forecast%info == 0 .and. &
      all(shape(var_forecast%paths) == [5, 2, 300]) .and. &
      all(var_forecast%upper > var_forecast%lower) .and. &
      all(ieee_is_finite(var_forecast%covariance)) .and. &
      var_diagnostics%info == 0 .and. &
      var_diagnostics%white_noise_test%info == 0 .and. &
      all(var_diagnostics%classification >= 1) .and. &
      all(var_diagnostics%classification <= 2), &
      'mixture VAR path forecasts and multivariate residual diagnostics')

   bayesian_model = mixar_model([0.58_dp, 0.42_dp], [1, 1], &
      [-0.35_dp, 0.45_dp], [0.28_dp, 0.34_dp], &
      reshape([0.45_dp, -0.20_dp], [1, 2]))
   call set_random_seed(19731)
   bayesian_simulation = mixar_simulate(bayesian_model, 500, [0.0_dp], 150)
   call set_random_seed(48291)
   bayesian_random = mixar_bayesian_random(499, 2, 1, 180)
   bayesian_draws = mixar_bayesian_sample_from_random( &
      bayesian_simulation%series, bayesian_model, 120, 60, bayesian_random, &
      proposal_scale=[0.06_dp, 0.06_dp])
   bayesian_repeat = mixar_bayesian_sample_from_random( &
      bayesian_simulation%series, bayesian_model, 120, 60, bayesian_random, &
      proposal_scale=[0.06_dp, 0.06_dp])
   call check(bayesian_draws%info == 0 .and. &
      all(shape(bayesian_draws%probability) == [120, 2]) .and. &
      maxval(abs(sum(bayesian_draws%probability, dim=2) - 1.0_dp)) < 1.0e-12_dp .and. &
      all(bayesian_draws%scale > 0.0_dp) .and. &
      all(bayesian_draws%acceptance_rate >= 0.0_dp) .and. &
      all(bayesian_draws%acceptance_rate <= 1.0_dp) .and. &
      maxval(abs(bayesian_draws%ar - bayesian_repeat%ar)) < 1.0e-15_dp .and. &
      bayesian_draws%posterior_mean_model%info == 0, &
      'supplied-random Bayesian MAR Gibbs and Metropolis sampling')

   relabeled_draws = mixar_bayesian_relabel(bayesian_draws)
   marginal_likelihood = mixar_marginal_likelihood( &
      bayesian_simulation%series, relabeled_draws)
   call check(relabeled_draws%info == 0 .and. &
      all(relabeled_draws%component_mean(:, 1) <= &
      relabeled_draws%component_mean(:, 2)) .and. &
      marginal_likelihood%info == 0 .and. &
      ieee_is_finite(marginal_likelihood%log_marginal_likelihood), &
      'Bayesian MAR label correction and marginal likelihood')

   call mixar_order_birth_death_probability(mixar_order_prior_poisson, &
      1.5_dp, 2, birth_probability, death_probability, order_status)
   call check(order_status == 0 .and. &
      abs(birth_probability - 0.5_dp) < 1.0e-14_dp .and. &
      abs(death_probability - 0.5_dp) < 1.0e-14_dp, &
      'Poisson component-order birth and death probabilities')
   call set_random_seed(73591)
   order_random = mixar_order_random(500, 2, 2, 80)
   order_selection = mixar_order_select_from_random( &
      bayesian_simulation%series, bayesian_model, 2, 80, order_random, &
      [0.06_dp, 0.06_dp], mixar_order_prior_flat)
   order_repeat = mixar_order_select_from_random( &
      bayesian_simulation%series, bayesian_model, 2, 80, order_random, &
      [0.06_dp, 0.06_dp], mixar_order_prior_flat)
   call check(order_selection%info == 0 .and. &
      all(shape(order_selection%order_trace) == [80, 2]) .and. &
      all(order_selection%order_trace >= 0) .and. &
      all(order_selection%order_trace <= 2) .and. &
      abs(sum(order_selection%posterior_probability) - 1.0_dp) < 1.0e-14_dp .and. &
      all(order_selection%modal_order >= 0) .and. &
      all(order_selection%modal_order <= 2) .and. &
      order_selection%modal_model%info == 0 .and. &
      all(order_selection%order_trace == order_repeat%order_trace) .and. &
      maxval(abs(order_selection%ar_trace - order_repeat%ar_trace)) < 1.0e-15_dp, &
      'reproducible reversible-jump component-order selection')

   call set_random_seed(91427)
   subsample_indices = mixar_random_subsample_indices(500, 1, 2, 5, 20)
   call set_random_seed(91427)
   repeated_indices = mixar_random_subsample_indices(500, 1, 2, 5, 20)
   random_initialization = mixar_initialize_from_indices( &
      bayesian_simulation%series, [1, 1], subsample_indices(:, :, 1))
   multistart = mixar_multistart_fit(bayesian_simulation%series, [1, 1], &
      subsample_indices, max_iterations=250, tolerance=1.0e-8_dp)
   call check(all(subsample_indices == repeated_indices) .and. &
      random_initialization%info == 0 .and. &
      random_initialization%model%info == 0 .and. &
      mixar_is_stable(random_initialization%model) .and. &
      all(random_initialization%model%scale > 0.0_dp) .and. &
      abs(sum(random_initialization%model%probability) - 1.0_dp) < 1.0e-13_dp .and. &
      multistart%info == 0 .and. multistart%successful >= 1 .and. &
      multistart%best >= 1 .and. multistart%best <= 5 .and. &
      multistart%fit(multistart%best)%converged .and. &
      multistart%fit(multistart%best)%log_likelihood > &
      multistart%initial_log_likelihood(multistart%best), &
      'stable data-driven initialization and Gaussian MAR multistart fitting')

   normal_fourth = mixar_standard_moment(0.0_dp, 4)
   student_fourth = mixar_standard_moment(6.0_dp, 4)
   undefined_fourth = mixar_standard_moment(4.0_dp, 4)
   innovation_second = mixar_innovation_moment(bayesian_model, 2)
   conditional_third = mixar_conditional_moment(single_model, [0.5_dp], 3)
   conditional_variance = mixar_conditional_central_moment( &
      single_model, [0.5_dp], 2)
   conditional_kurtosis = mixar_conditional_kurtosis(single_model, [0.5_dp])
   stationary_moments = mixar_stationary_moments(single_model)
   call check(normal_fourth%exists .and. &
      abs(normal_fourth%value - 3.0_dp) < 1.0e-13_dp .and. &
      student_fourth%exists .and. &
      abs(student_fourth%value - 6.0_dp) < 1.0e-12_dp .and. &
      .not. undefined_fourth%exists .and. &
      innovation_second%exists .and. &
      abs(innovation_second%value - sum(bayesian_model%probability* &
      bayesian_model%scale**2)) < 1.0e-13_dp .and. &
      conditional_third%exists .and. &
      abs(conditional_third%value - (0.7_dp**3 + &
      3.0_dp*0.7_dp*0.3_dp**2)) < 1.0e-13_dp .and. &
      abs(conditional_variance%value - 0.3_dp**2) < 1.0e-13_dp .and. &
      abs(conditional_kurtosis%value - 3.0_dp) < 1.0e-12_dp .and. &
      stationary_moments%info == 0 .and. &
      abs(stationary_moments%mean - 1.0_dp) < 1.0e-12_dp .and. &
      abs(stationary_moments%variance - 0.140625_dp) < 1.0e-12_dp, &
      'innovation, conditional, kurtosis, and stationary MAR moments')

   print '(a)', 'mixAR tests passed'

contains

   subroutine check(condition, message)
      !! Stop the test program when a condition fails.
      logical, intent(in) :: condition !! Test condition.
      character(len=*), intent(in) :: message !! Failure message.

      if (.not. condition) then
         print '(a)', 'FAILED: '//trim(message)
         error stop 1
      end if
   end subroutine check

end program test_mixar
