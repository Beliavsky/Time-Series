! SPDX-License-Identifier: GPL-3.0-only
! SPDX-FileComment: Numerical tests for algorithms translated from R gmvarkit.
program test_gmvarkit
   use kind_mod, only: dp
   use gmvarkit_mod, only: gmvarkit_model_t, gmvarkit_regime_moments_t, &
      gmvarkit_evaluation_t, gmvarkit_regime_moments, gmvarkit_evaluate
   use gmvarkit_mod, only: gmvarkit_location_parameters_t, &
      gmvarkit_location_parameters, gmvarkit_model_from_regime_means
   use gmvarkit_mod, only: gmvarkit_simulation_t, gmvarkit_forecast_t, &
      gmvarkit_simulate, gmvarkit_forecast
   use gmvarkit_mod, only: gmvarkit_girf_t, gmvarkit_gfevd_t, &
      gmvarkit_girf, gmvarkit_gfevd
   use gmvarkit_mod, only: gmvarkit_girf_inference_t, &
      gmvarkit_gfevd_inference_t, gmvarkit_girf_inference, &
      gmvarkit_gfevd_inference
   use gmvarkit_mod, only: gmvarkit_linear_irf_t, gmvarkit_linear_irf
   use gmvarkit_mod, only: gmvarkit_unconditional_moments_t, &
      gmvarkit_pearson_residuals_t, gmvarkit_information_criteria_t, &
      gmvarkit_unconditional_moments, gmvarkit_pearson_residuals, &
      gmvarkit_information_criteria
   use gmvarkit_mod, only: gmvarkit_structural_t, &
      gmvarkit_structural_covariances, gmvarkit_identify_structural, &
      gmvarkit_reorder_structural, gmvarkit_swap_structural_signs, &
      gmvarkit_rebase_structural
   use gmvarkit_mod, only: gmvarkit_linear_irf_bootstrap_t, &
      gmvarkit_linear_irf_bootstrap
   use gmvarkit_mod, only: gmvarkit_structural_constraints_t, &
      gmvarkit_structural_fit_t, gmvarkit_estimate_structural
   use gmvarkit_mod, only: gmvarkit_structural_multistart_fit_t, &
      gmvarkit_structural_multistart_estimate
   use gmvarkit_mod, only: gmvarkit_likelihood_profile_t, &
      gmvarkit_profile_likelihood
   use gmvarkit_mod, only: gmvarkit_regime_conversion_t, &
      gmvarkit_convert_student_regimes
   use gmvarkit_mod, only: gmvarkit_companion_eigen_t, &
      gmvarkit_covariance_eigen_t, gmvarkit_companion_eigenvalues, &
      gmvarkit_covariance_eigenvalues
   use gmvarkit_mod, only: gmvarkit_fit_t, gmvarkit_estimate
   use gmvarkit_mod, only: gmvarkit_genetic_fit_t, gmvarkit_genetic_estimate
   use gmvarkit_mod, only: gmvarkit_multistart_fit_t, &
      gmvarkit_multistart_estimate
   use gmvarkit_mod, only: gmvarkit_constraints_t, &
      gmvarkit_estimate_constrained
   use gmvarkit_mod, only: gmvarkit_inference_t, &
      gmvarkit_hypothesis_test_t, gmvarkit_inference, gmvarkit_wald_test, &
      gmvarkit_likelihood_ratio
   use gmvarkit_mod, only: gmvarkit_score_t, gmvarkit_score_matrix, &
      gmvarkit_rao_test, gmvarkit_quantile_residuals_t, &
      gmvarkit_residual_tests_t, gmvarkit_quantile_residuals, &
      gmvarkit_quantile_residual_tests
   use random_mod, only: set_random_seed
   implicit none

   type(gmvarkit_model_t) :: model
   type(gmvarkit_regime_moments_t) :: moments
   type(gmvarkit_location_parameters_t) :: location, changed_location
   type(gmvarkit_model_t) :: mean_model
   type(gmvarkit_evaluation_t) :: conditional, exact
   type(gmvarkit_simulation_t) :: simulation
   type(gmvarkit_forecast_t) :: forecast
   type(gmvarkit_girf_t) :: girf
   type(gmvarkit_gfevd_t) :: gfevd
   type(gmvarkit_girf_inference_t) :: girf_inference
   type(gmvarkit_gfevd_inference_t) :: gfevd_inference
   type(gmvarkit_linear_irf_t) :: linear_irf
   type(gmvarkit_unconditional_moments_t) :: unconditional_moments
   type(gmvarkit_pearson_residuals_t) :: pearson_residuals
   type(gmvarkit_information_criteria_t) :: criteria
   type(gmvarkit_structural_t) :: identified, reconstructed
   type(gmvarkit_structural_t) :: reordered, sign_swapped, rebased
   type(gmvarkit_structural_t) :: bootstrap_structural
   type(gmvarkit_linear_irf_bootstrap_t) :: irf_bootstrap
   type(gmvarkit_structural_constraints_t) :: structural_constraints
   type(gmvarkit_structural_fit_t) :: structural_fit
   type(gmvarkit_structural_multistart_fit_t) :: structural_multistart
   type(gmvarkit_model_t) :: structural_start_model(2)
   type(gmvarkit_structural_t) :: structural_start(2)
   type(gmvarkit_likelihood_profile_t) :: profile, structural_profile
   type(gmvarkit_regime_conversion_t) :: conversion, reestimated_conversion
   type(gmvarkit_regime_conversion_t) :: unchanged_conversion
   type(gmvarkit_model_t) :: student_model
   type(gmvarkit_structural_t) :: student_structural
   type(gmvarkit_companion_eigen_t) :: companion_eigen, companion_alert
   type(gmvarkit_covariance_eigen_t) :: covariance_eigen
   type(gmvarkit_model_t) :: initial_model
   type(gmvarkit_fit_t) :: fit
   type(gmvarkit_genetic_fit_t) :: genetic_fit
   type(gmvarkit_multistart_fit_t) :: multistart_fit
   type(gmvarkit_model_t) :: starting_model(3)
   type(gmvarkit_model_t) :: constrained_model
   type(gmvarkit_constraints_t) :: constraints
   type(gmvarkit_inference_t) :: inference
   type(gmvarkit_hypothesis_test_t) :: hypothesis
   type(gmvarkit_fit_t) :: unrestricted_fit, restricted_fit
   type(gmvarkit_score_t) :: score
   type(gmvarkit_quantile_residuals_t) :: quantile_residuals
   type(gmvarkit_residual_tests_t) :: residual_tests
   type(gmvarkit_evaluation_t) :: initial_evaluation, fitted_evaluation
   real(dp) :: series(5, 2)
   real(dp) :: estimation_series(40, 1)
   integer :: time

   allocate(model%intercept(2, 2), model%ar(2, 2, 1, 2))
   allocate(model%covariance(2, 2, 2), model%weight(2))
   allocate(model%degrees_of_freedom(2))
   model%intercept = reshape([0.0_dp, 0.0_dp, 1.0_dp, -0.5_dp], [2, 2])
   model%ar = 0.0_dp
   model%ar(1, 1, 1, 1) = 0.5_dp
   model%ar(2, 2, 1, 1) = 0.5_dp
   model%ar(1, 1, 1, 2) = 0.2_dp
   model%ar(2, 2, 1, 2) = 0.2_dp
   model%covariance = 0.0_dp
   model%covariance(1, 1, 1) = 1.0_dp
   model%covariance(2, 2, 1) = 1.0_dp
   model%covariance(1, 1, 2) = 2.0_dp
   model%covariance(2, 2, 2) = 0.5_dp
   model%weight = [0.6_dp, 0.4_dp]
   model%degrees_of_freedom = [0.0_dp, 8.0_dp]
   model%gaussian_regimes = 1
   location = gmvarkit_location_parameters(model)
   mean_model = gmvarkit_model_from_regime_means(model, location%mean)
   changed_location = gmvarkit_location_parameters( &
      gmvarkit_model_from_regime_means(model, location%mean + 1.0_dp))
   call check(location%info == 0 .and. mean_model%info == 0 .and. &
      maxval(abs(location%intercept - model%intercept)) < 1.0e-14_dp .and. &
      maxval(abs(location%mean - reshape([0.0_dp, 0.0_dp, 1.25_dp, &
      -0.625_dp], [2, 2]))) < 1.0e-14_dp .and. &
      maxval(abs(mean_model%intercept - model%intercept)) < 1.0e-14_dp .and. &
      maxval(abs(changed_location%mean - (location%mean + 1.0_dp))) < &
      1.0e-14_dp, 'intercept and stationary-mean parameterization round trip')
   companion_eigen = gmvarkit_companion_eigenvalues(model)
   companion_alert = gmvarkit_companion_eigenvalues(model, tolerance=0.6_dp)
   covariance_eigen = gmvarkit_covariance_eigenvalues(model, &
      positive_definite_tolerance=0.75_dp, identification_tolerance=0.8_dp)
   call check(companion_eigen%info == 0 .and. &
      all(shape(companion_eigen%eigenvalue) == [2, 2]) .and. &
      maxval(abs(companion_eigen%spectral_radius - [0.5_dp, 0.2_dp])) < &
      1.0e-12_dp .and. all(companion_eigen%stationary) .and. &
      .not. any(companion_eigen%near_unit_root) .and. &
      all(companion_alert%near_unit_root .eqv. [.true., .false.]) .and. &
      covariance_eigen%info == 0 .and. &
      all(shape(covariance_eigen%eigenvalue) == [2, 2]) .and. &
      all(covariance_eigen%pair(:, 1) == [1, 2]) .and. &
      maxval(abs(covariance_eigen%ratio_eigenvalue(:, 1) - &
      [2.0_dp, 0.5_dp])) < 1.0e-12_dp .and. &
      abs(covariance_eigen%minimum_separation(1) - 0.75_dp) < 1.0e-12_dp .and. &
      all(covariance_eigen%near_singular .eqv. [.false., .true.]) .and. &
      covariance_eigen%weakly_identified(1), &
      'companion stability and covariance eigenvalue diagnostics')
   moments = gmvarkit_regime_moments(model)
   call check(moments%info == 0 .and. &
      all(abs(moments%mean(:, 1)) < 1.0e-14_dp) .and. &
      maxval(abs(moments%mean(:, 2) - [1.25_dp, -0.625_dp])) < 1.0e-14_dp .and. &
      abs(moments%lag_covariance(1, 1, 1) - 4.0_dp/3.0_dp) < 1.0e-12_dp, &
      'stationary regime moments')

   series = reshape([0.1_dp, -0.2_dp, 0.3_dp, 0.1_dp, -0.1_dp, 0.4_dp, &
      0.2_dp, 0.3_dp, 0.5_dp, -0.1_dp], [5, 2])
   conditional = gmvarkit_evaluate(series, model)
   exact = gmvarkit_evaluate(series, model, conditional=.false.)
   call check(conditional%info == 0 .and. &
      all(shape(conditional%mixing_weight) == [4, 2]) .and. &
      maxval(abs(sum(conditional%mixing_weight, dim=2) - 1.0_dp)) < &
      1.0e-12_dp .and. all(conditional%mixing_weight > 0.0_dp), &
      'endogenous regime mixing weights')
   call check(all(abs(conditional%arch_scalar(:, 1) - 1.0_dp) < 1.0e-14_dp) .and. &
      all(conditional%arch_scalar(:, 2) > 0.0_dp) .and. &
      all(conditional%conditional_covariance(1, 1, :) > 0.0_dp) .and. &
      all(conditional%conditional_covariance(2, 2, :) > 0.0_dp), &
      'Gaussian and Student-t conditional covariances')
   call check(abs(conditional%log_likelihood) < huge(1.0_dp) .and. &
      exact%info == 0 .and. exact%log_likelihood < conditional%log_likelihood, &
      'conditional and exact mixture log likelihood')
   unconditional_moments = gmvarkit_unconditional_moments(model)
   call check(unconditional_moments%info == 0 .and. &
      maxval(abs(unconditional_moments%mean - [0.5_dp, -0.25_dp])) < &
      1.0e-14_dp .and. &
      all(shape(unconditional_moments%autocovariance) == [2, 2, 2]) .and. &
      all(shape(unconditional_moments%regime_autocovariance) == [2, 2, 2, 2]) .and. &
      abs(unconditional_moments%autocovariance(1, 1, 1) - &
      241.0_dp/120.0_dp) < 1.0e-12_dp .and. &
      maxval(abs([(unconditional_moments%autocorrelation(time, time, 1), &
      time=1, 2)] - 1.0_dp)) < 1.0e-14_dp, &
      'unconditional mixture autocovariances and autocorrelations')
   pearson_residuals = gmvarkit_pearson_residuals(series, model, &
      standardize=.false.)
   call check(pearson_residuals%info == 0 .and. &
      .not. pearson_residuals%standardized .and. &
      all(shape(pearson_residuals%residual) == [4, 2]) .and. &
      maxval(abs(pearson_residuals%residual - &
      (series(2:, :) - transpose(conditional%conditional_mean)))) < &
      1.0e-14_dp, 'raw Pearson residuals')
   pearson_residuals = gmvarkit_pearson_residuals(series, model)
   call check(pearson_residuals%info == 0 .and. &
      pearson_residuals%standardized .and. &
      all(abs(pearson_residuals%residual) < huge(1.0_dp)), &
      'symmetric-covariance-standardized Pearson residuals')
   criteria = gmvarkit_information_criteria(-10.0_dp, 5, 100)
   call check(criteria%info == 0 .and. abs(criteria%aic - 30.0_dp) < &
      1.0e-14_dp .and. abs(criteria%hqic - &
      (20.0_dp + 10.0_dp*log(log(100.0_dp)))) < 1.0e-14_dp .and. &
      abs(criteria%bic - (20.0_dp + 5.0_dp*log(100.0_dp))) < 1.0e-14_dp, &
      'AIC, HQIC, and BIC')
   identified = gmvarkit_identify_structural(model)
   reconstructed = gmvarkit_structural_covariances(identified%impact, &
      identified%relative_variance)
   call check(identified%info == 0 .and. reconstructed%info == 0 .and. &
      identified%reconstruction_error < 1.0e-12_dp .and. &
      maxval(abs(identified%relative_variance(:, 1) - 1.0_dp)) < &
      1.0e-12_dp .and. &
      maxval(abs(reconstructed%covariance - model%covariance)) < 1.0e-12_dp, &
      'heteroskedastic structural covariance identification')
   reordered = gmvarkit_reorder_structural(identified, [2, 1])
   sign_swapped = gmvarkit_swap_structural_signs(identified, [1])
   rebased = gmvarkit_rebase_structural(identified, 2)
   call check(reordered%info == 0 .and. sign_swapped%info == 0 .and. &
      rebased%info == 0 .and. &
      maxval(abs(reordered%impact(:, 1) - identified%impact(:, 2))) < &
      1.0e-14_dp .and. &
      maxval(abs(sign_swapped%impact(:, 1) + identified%impact(:, 1))) < &
      1.0e-14_dp .and. &
      maxval(abs(rebased%relative_variance(:, 2) - 1.0_dp)) < 1.0e-14_dp .and. &
      maxval(abs(rebased%covariance - model%covariance)) < 1.0e-14_dp, &
      'structural shock ordering, signs, and reference regime')
   call set_random_seed(1207)
   girf_inference = gmvarkit_girf_inference(model, series, 2, &
      inner_simulations=20, outer_replications=4, confidence_level=0.8_dp, &
      initial_value_mode='random', cumulative_variables=[1], &
      scale_shocks=[1], scale_variables=[1], scale_values=[1.0_dp], &
      scale_type='instant', structural=identified)
   gfevd_inference = gmvarkit_gfevd_inference(girf_inference)
   call check(girf_inference%info == 0 .and. &
      girf_inference%outer_replications == 4 .and. &
      all(shape(girf_inference%response_draw) == [3, 2, 2, 4]) .and. &
      all(girf_inference%lower_response <= girf_inference%upper_response) .and. &
      maxval(abs(girf_inference%response_draw(1, 1, 1, :) - 1.0_dp)) < &
      1.0e-12_dp, 'random-history GIRF inference, accumulation, and scaling')
   call check(gfevd_inference%info == 0 .and. &
      all(shape(gfevd_inference%decomposition) == [3, 2, 2]) .and. &
      all(shape(gfevd_inference%mixing_weight_decomposition) == [3, 2, 2]) .and. &
      maxval(abs(sum(gfevd_inference%decomposition, dim=3) - 1.0_dp)) < &
      1.0e-12_dp .and. minval(gfevd_inference%mixing_weight_decomposition) >= &
      0.0_dp .and. maxval(gfevd_inference%mixing_weight_decomposition) <= &
      1.0_dp, 'history-averaged variable and mixing-weight GFEVDs')
   call set_random_seed(1208)
   girf_inference = gmvarkit_girf_inference(model, series, 1, &
      inner_simulations=5, initial_value_mode='data')
   call check(girf_inference%info == 0 .and. &
      girf_inference%outer_replications == 5, &
      'GIRF inference over every observed initial history')
   call set_random_seed(1209)
   girf_inference = gmvarkit_girf_inference(model, series, 1, &
      inner_simulations=5, initial_value_mode='fixed', &
      fixed_initial_values=series)
   call check(girf_inference%info == 0 .and. &
      girf_inference%outer_replications == 1 .and. &
      maxval(abs(girf_inference%lower_response - &
      girf_inference%upper_response)) < 1.0e-14_dp, &
      'fixed-history GIRF inference')
   quantile_residuals = gmvarkit_quantile_residuals(series, model)
   call check(quantile_residuals%info == 0 .and. &
      all(shape(quantile_residuals%residual) == [4, 2]) .and. &
      all(quantile_residuals%cdf > 0.0_dp) .and. &
      all(quantile_residuals%cdf < 1.0_dp) .and. &
      all(abs(quantile_residuals%residual) < huge(1.0_dp)), &
      'sequential mixed-regime quantile residuals')
   call set_random_seed(1201)
   simulation = gmvarkit_simulate(model, series, 6, paths=3)
   call check(simulation%info == 0 .and. &
      all(shape(simulation%series) == [6, 2, 3]) .and. &
      all(shape(simulation%mixing_weight) == [6, 2, 3]) .and. &
      all(simulation%regime >= 1) .and. all(simulation%regime <= 2) .and. &
      maxval(abs(sum(simulation%mixing_weight, dim=2) - 1.0_dp)) < 1.0e-12_dp, &
      'recursive mixed-regime simulation')
   call set_random_seed(1202)
   forecast = gmvarkit_forecast(model, series, 3, simulations=100, &
      probabilities=[0.1_dp, 0.9_dp])
   call check(forecast%info == 0 .and. forecast%simulations == 100 .and. &
      all(shape(forecast%mean) == [3, 2]) .and. &
      all(shape(forecast%quantile) == [3, 2, 2]) .and. &
      all(forecast%quantile(:, 1, :) <= forecast%quantile(:, 2, :)) .and. &
      maxval(abs(sum(forecast%mixing_weight_mean, dim=2) - 1.0_dp)) < &
      1.0e-12_dp, 'simulation-based forecasts and quantiles')
   call set_random_seed(1204)
   girf = gmvarkit_girf(model, series, 3, simulations=200, shock_size=1.0_dp)
   gfevd = gmvarkit_gfevd(girf)
   call check(girf%info == 0 .and. girf%simulations == 200 .and. &
      all(shape(girf%response) == [4, 2, 2]) .and. &
      all(shape(girf%mixing_weight_response) == [4, 2, 2]) .and. &
      abs(girf%response(1, 1, 2)) < 1.0e-12_dp .and. &
      maxval(abs(girf%mixing_weight_response(1, :, :))) < 1.0e-14_dp, &
      'paired nonlinear generalized impulse responses')
   call check(gfevd%info == 0 .and. &
      all(shape(gfevd%decomposition) == [4, 2, 2]) .and. &
      minval(gfevd%decomposition) >= 0.0_dp .and. &
      maxval(gfevd%decomposition) <= 1.0_dp .and. &
      maxval(abs(sum(gfevd%decomposition, dim=3) - 1.0_dp)) < 1.0e-12_dp, &
      'cumulative generalized forecast-error variance decomposition')
   linear_irf = gmvarkit_linear_irf(model, 3, regime=1, &
      cumulative_variables=[1], scale_shocks=[2], scale_variables=[2], &
      scale_values=[2.0_dp])
   call check(linear_irf%info == 0 .and. linear_irf%regime == 1 .and. &
      all(shape(linear_irf%response) == [4, 2, 2]) .and. &
      maxval(abs(linear_irf%response(:, 1, 1) - &
      [1.0_dp, 1.5_dp, 1.75_dp, 1.875_dp])) < 1.0e-14_dp .and. &
      maxval(abs(linear_irf%response(:, 2, 2) - &
      [2.0_dp, 1.0_dp, 0.5_dp, 0.25_dp])) < 1.0e-14_dp .and. &
      abs(linear_irf%response(1, 1, 2)) < 1.0e-14_dp, &
      'regime linear impulse responses with accumulation and scaling')
   linear_irf = gmvarkit_linear_irf(model, 2, regime=2)
   call check(linear_irf%info == 0 .and. &
      abs(linear_irf%response(1, 1, 1) - sqrt(2.0_dp)) < 1.0e-14_dp .and. &
      abs(linear_irf%response(2, 1, 1) - 0.2_dp*sqrt(2.0_dp)) < &
      1.0e-14_dp, 'selected-regime Cholesky impulse responses')
   linear_irf = gmvarkit_linear_irf(model, 2, regime=2, structural=reordered)
   call check(linear_irf%info == 0 .and. &
      maxval(abs(linear_irf%impact - reordered%impact*spread( &
      sqrt(reordered%relative_variance(:, 2)), 1, 2))) < 1.0e-14_dp .and. &
      abs(linear_irf%response(1, 1, 1)) < 1.0e-14_dp, &
      'heteroskedastically identified linear impulse responses')
   linear_irf = gmvarkit_linear_irf(model, 2, scale_shocks=[2], &
      scale_variables=[1], scale_values=[1.0_dp])
   call check(linear_irf%info == 9, 'reject scaling a zero impact response')
   call set_random_seed(1205)
   girf = gmvarkit_girf(model, series, 2, simulations=100, &
      structural=reordered)
   call check(girf%info == 0 .and. &
      abs(girf%response(1, 1, 1)) < 1.0e-12_dp .and. &
      maxval(abs(girf%mixing_weight_response(1, :, :))) < 1.0e-14_dp, &
      'heteroskedastically identified generalized impulse responses')
   allocate(initial_model%intercept(1, 1), initial_model%ar(1, 1, 1, 1))
   allocate(initial_model%covariance(1, 1, 1), initial_model%weight(1))
   allocate(initial_model%degrees_of_freedom(1))
   initial_model%intercept = 0.0_dp
   initial_model%ar = 0.1_dp
   initial_model%covariance = 0.5_dp
   initial_model%weight = 1.0_dp
   initial_model%degrees_of_freedom = 0.0_dp
   initial_model%gaussian_regimes = 1
   estimation_series(1, 1) = 0.0_dp
   do time = 2, size(estimation_series, 1)
      estimation_series(time, 1) = 0.2_dp + &
         0.55_dp*estimation_series(time - 1, 1) + &
         0.12_dp*sin(1.7_dp*real(time, dp))
   end do
   allocate(student_model%intercept(1, 3), student_model%ar(1, 1, 1, 3))
   allocate(student_model%covariance(1, 1, 3), student_model%weight(3))
   allocate(student_model%degrees_of_freedom(3))
   student_model%intercept = reshape([0.1_dp, 0.2_dp, 0.3_dp], [1, 3])
   student_model%ar = reshape([0.1_dp, 0.2_dp, 0.3_dp], [1, 1, 1, 3])
   student_model%covariance = reshape([1.0_dp, 4.0_dp, 9.0_dp], [1, 1, 3])
   student_model%weight = [0.2_dp, 0.5_dp, 0.3_dp]
   student_model%degrees_of_freedom = [0.0_dp, 150.0_dp, 8.0_dp]
   student_model%gaussian_regimes = 1
   student_structural = gmvarkit_structural_covariances( &
      reshape([1.0_dp], [1, 1]), reshape([1.0_dp, 4.0_dp, 9.0_dp], &
      [1, 3]), 1)
   conversion = gmvarkit_convert_student_regimes(student_model, 100.0_dp, &
      structural=student_structural)
   unchanged_conversion = gmvarkit_convert_student_regimes(student_model, &
      200.0_dp)
   reestimated_conversion = gmvarkit_convert_student_regimes(student_model, &
      100.0_dp, series=estimation_series, estimate=.true., max_iterations=1, &
      tolerance=1.0e-4_dp)
   call check(conversion%info == 0 .and. conversion%has_structural .and. &
      .not. conversion%reestimated .and. &
      all(conversion%converted .eqv. [.false., .true., .false.]) .and. &
      all(conversion%new_to_old == [2, 1, 3]) .and. &
      all(conversion%old_to_new == [2, 1, 3]) .and. &
      conversion%fit%model%gaussian_regimes == 2 .and. &
      all(conversion%fit%model%degrees_of_freedom == &
      [0.0_dp, 0.0_dp, 8.0_dp]) .and. &
      maxval(abs(conversion%fit%model%intercept - &
      reshape([0.2_dp, 0.1_dp, 0.3_dp], [1, 3]))) < 1.0e-14_dp .and. &
      conversion%structural%reference_regime == 2 .and. &
      maxval(abs(conversion%structural%relative_variance - &
      reshape([4.0_dp, 1.0_dp, 9.0_dp], [1, 3]))) < 1.0e-14_dp .and. &
      .not. any(unchanged_conversion%converted) .and. &
      all(unchanged_conversion%new_to_old == [1, 2, 3]) .and. &
      unchanged_conversion%fit%model%gaussian_regimes == 1 .and. &
      reestimated_conversion%reestimated .and. &
      (reestimated_conversion%info == 0 .or. &
      reestimated_conversion%info == 4) .and. &
      reestimated_conversion%fit%model%gaussian_regimes == 2, &
      'large-df Student conversion, regime mapping, and optional refit')
   initial_evaluation = gmvarkit_evaluate(estimation_series, initial_model)
   fit = gmvarkit_estimate(estimation_series, initial_model, max_iterations=80, &
      tolerance=1.0e-5_dp)
   fitted_evaluation = gmvarkit_evaluate(estimation_series, fit%model)
   call check(fitted_evaluation%info == 0 .and. &
      fit%log_likelihood >= initial_evaluation%log_likelihood - 1.0e-10_dp .and. &
      fit%model%covariance(1, 1, 1) > 0.0_dp .and. &
      abs(fit%model%weight(1) - 1.0_dp) < 1.0e-14_dp .and. &
      fit%iterations > 0, 'local transformed-parameter likelihood estimation')
   bootstrap_structural = gmvarkit_identify_structural(fit%model)
   allocate(structural_constraints%impact_fixed(1, 1), &
      structural_constraints%impact_value(1, 1), &
      structural_constraints%impact_sign(1, 1))
   structural_constraints%impact_fixed = .false.
   structural_constraints%impact_value = 0.0_dp
   structural_constraints%impact_sign = 1
   structural_fit = gmvarkit_estimate_structural(estimation_series, fit%model, &
      bootstrap_structural, structural_constraints, max_iterations=30, &
      tolerance=1.0e-5_dp)
   call check((structural_fit%info == 0 .or. structural_fit%info == 4) .and. &
      structural_fit%structural%impact(1, 1) > 0.0_dp .and. &
      abs(structural_fit%fit%model%covariance(1, 1, 1) - &
      structural_fit%structural%impact(1, 1)**2) < 1.0e-12_dp .and. &
      structural_fit%fit%log_likelihood >= fit%log_likelihood - 1.0e-8_dp .and. &
      structural_fit%inference%info == 0 .and. &
      all(structural_fit%inference%standard_error >= 0.0_dp), &
      'direct sign-constrained structural likelihood and Hessian inference')
   profile = gmvarkit_profile_likelihood(estimation_series, fit%model, &
      parameter=[2], scale=0.05_dp, points=5)
   structural_profile = gmvarkit_profile_likelihood(estimation_series, &
      structural_fit%fit%model, structural_fit%structural, &
      constraints=structural_constraints, parameter=[3], scale=0.05_dp, &
      points=5)
   call check(profile%info == 0 .and. .not. profile%structural .and. &
      all(shape(profile%value) == [5, 1]) .and. all(profile%valid) .and. &
      profile%parameter(1) == 2 .and. &
      abs(profile%value(3, 1) - profile%center(1)) < 1.0e-14_dp .and. &
      abs(profile%log_likelihood(3, 1) - fit%log_likelihood) < 1.0e-10_dp .and. &
      structural_profile%info == 0 .and. structural_profile%structural .and. &
      all(structural_profile%valid) .and. &
      abs(structural_profile%log_likelihood(3, 1) - &
      structural_fit%fit%log_likelihood) < 1.0e-10_dp, &
      'reduced-form and structural transformed-coordinate likelihood profiles')
   structural_start_model = structural_fit%fit%model
   structural_start = structural_fit%structural
   structural_multistart = gmvarkit_structural_multistart_estimate( &
      estimation_series, structural_start_model, structural_start, &
      constraints=structural_constraints, max_iterations=20, &
      tolerance=1.0e-5_dp)
   call check(structural_multistart%info == 0 .and. &
      structural_multistart%successful_count == 2 .and. &
      structural_multistart%best_index == 1 .and. &
      structural_multistart%distinct(1) .and. &
      .not. structural_multistart%distinct(2) .and. &
      structural_multistart%duplicate_of(2) == 1 .and. &
      allocated(structural_multistart%fit(1)%inference%hessian) .and. &
      .not. allocated(structural_multistart%fit(2)%inference%hessian) .and. &
      structural_multistart%fit(1)%inference%info == 0, &
      'structural multistart ranking, duplicates, and best-only inference')
   call set_random_seed(1206)
   irf_bootstrap = gmvarkit_linear_irf_bootstrap(estimation_series, fit, 3, &
      replications=8, confidence_level=0.8_dp, &
      structural=bootstrap_structural, max_iterations=30, tolerance=1.0e-4_dp)
   call check(irf_bootstrap%info == 0 .and. &
      irf_bootstrap%requested_replications == 8 .and. &
      irf_bootstrap%successful_replications >= 2 .and. &
      all(shape(irf_bootstrap%point) == [4, 1, 1]) .and. &
      size(irf_bootstrap%draw, 4) == irf_bootstrap%successful_replications .and. &
      all(irf_bootstrap%lower <= irf_bootstrap%upper), &
      'fixed-design wild-bootstrap linear impulse-response intervals')
   deallocate(structural_constraints%impact_fixed, &
      structural_constraints%impact_value, structural_constraints%impact_sign)
   allocate(structural_constraints%impact_fixed(2, 2), &
      structural_constraints%impact_value(2, 2), &
      structural_constraints%impact_sign(2, 2), &
      structural_constraints%lambda_mapping(2, 2))
   structural_constraints%impact_fixed = .false.
   structural_constraints%impact_fixed(1, 2) = .true.
   structural_constraints%impact_fixed(2, 1) = .true.
   structural_constraints%impact_value = 0.0_dp
   structural_constraints%impact_sign = 0
   structural_constraints%impact_sign(1, 1) = 1
   structural_constraints%impact_sign(2, 2) = 1
   structural_constraints%lambda_mapping = 0.0_dp
   structural_constraints%lambda_mapping(1, 1) = 1.0_dp
   structural_constraints%lambda_mapping(2, 2) = 1.0_dp
   structural_fit = gmvarkit_estimate_structural(series, model, identified, &
      structural_constraints, max_iterations=2, tolerance=1.0e-4_dp, &
      calculate_inference=.false.)
   call check((structural_fit%info == 0 .or. structural_fit%info == 4) .and. &
      structural_fit%structural%impact(1, 2) == 0.0_dp .and. &
      structural_fit%structural%impact(2, 1) == 0.0_dp .and. &
      all(structural_fit%structural%relative_variance > 0.0_dp) .and. &
      maxval(abs(structural_fit%fit%model%covariance - &
      structural_fit%structural%covariance)) < 1.0e-14_dp, &
      'zero-restricted W and positive linear lambda mapping')
   inference = gmvarkit_inference(estimation_series, fit, &
      difference_step=1.0e-4_dp)
   call check(inference%info == 0 .and. size(inference%parameter) == 3 .and. &
      all(shape(inference%covariance) == [3, 3]) .and. &
      all(inference%standard_error >= 0.0_dp), &
      'Hessian covariance and transformed-scale standard errors')
   hypothesis = gmvarkit_wald_test(inference, &
      reshape([0.0_dp, 1.0_dp, 0.0_dp], [1, 3]), [0.0_dp])
   call check(hypothesis%info == 0 .and. hypothesis%statistic >= 0.0_dp .and. &
      hypothesis%degrees_of_freedom == 1 .and. hypothesis%p_value >= 0.0_dp .and. &
      hypothesis%p_value <= 1.0_dp, 'linear Wald hypothesis test')
   score = gmvarkit_score_matrix(estimation_series, fit%model, &
      difference_step=1.0e-4_dp)
   hypothesis = gmvarkit_rao_test(score, 1)
   call check(score%info == 0 .and. &
      all(shape(score%observation) == [39, 3]) .and. &
      all(shape(score%opg) == [3, 3]) .and. hypothesis%info == 0 .and. &
      hypothesis%statistic >= 0.0_dp .and. hypothesis%p_value >= 0.0_dp .and. &
      hypothesis%p_value <= 1.0_dp, 'observationwise scores and Rao test')
   quantile_residuals = gmvarkit_quantile_residuals(estimation_series, fit%model)
   residual_tests = gmvarkit_quantile_residual_tests( &
      quantile_residuals%residual, [1, 2], [1, 2], &
      series=estimation_series, model=fit%model, difference_step=1.0e-4_dp)
   call check(residual_tests%info == 0 .and. &
      residual_tests%parameter_corrected .and. &
      residual_tests%normality%info == 0 .and. &
      all([(residual_tests%autocorrelation(time)%info == 0, time=1, 2)]) .and. &
      all([(residual_tests%heteroskedasticity(time)%info == 0, time=1, 2)]) .and. &
      residual_tests%autocorrelation(2)%degrees_of_freedom == 2 .and. &
      residual_tests%heteroskedasticity(2)%degrees_of_freedom == 2, &
      'quantile-residual moment diagnostics')
   unrestricted_fit%log_likelihood = -10.0_dp
   unrestricted_fit%parameter_count = 5
   restricted_fit%log_likelihood = -12.0_dp
   restricted_fit%parameter_count = 3
   hypothesis = gmvarkit_likelihood_ratio(unrestricted_fit, restricted_fit)
   call check(hypothesis%info == 0 .and. &
      abs(hypothesis%statistic - 4.0_dp) < 1.0e-14_dp .and. &
      hypothesis%degrees_of_freedom == 2 .and. &
      abs(hypothesis%p_value - exp(-2.0_dp)) < 1.0e-12_dp, &
      'nested-model likelihood-ratio test')
   starting_model(1) = initial_model
   starting_model(2) = initial_model
   starting_model(3) = initial_model
   starting_model(3)%intercept = 0.7_dp
   multistart_fit = gmvarkit_multistart_estimate(estimation_series, &
      starting_model, max_iterations=30, tolerance=1.0e-5_dp)
   call check(multistart_fit%info == 0 .and. &
      size(multistart_fit%fit) == 3 .and. &
      multistart_fit%successful_count == 3 .and. &
      multistart_fit%best_index == multistart_fit%order(1) .and. &
      multistart_fit%fit(multistart_fit%order(1))%log_likelihood >= &
      multistart_fit%fit(multistart_fit%order(2))%log_likelihood .and. &
      multistart_fit%distinct(1) .and. .not. multistart_fit%distinct(2) .and. &
      multistart_fit%duplicate_of(2) == 1 .and. &
      multistart_fit%distinct_count >= 1, &
      'ranked multistart estimation and canonical duplicate detection')
   call set_random_seed(1203)
   genetic_fit = gmvarkit_genetic_estimate(estimation_series, initial_model, &
      population_size=8, generations=4, mutation_scale=0.15_dp, &
      local_iterations=30)
   fitted_evaluation = gmvarkit_evaluate(estimation_series, genetic_fit%fit%model)
   call check(genetic_fit%info == 0 .and. fitted_evaluation%info == 0 .and. &
      genetic_fit%evaluations == 40 .and. &
      all(genetic_fit%best_objective(2:) <= &
      genetic_fit%best_objective(:3) + 1.0e-12_dp) .and. &
      genetic_fit%fit%log_likelihood >= initial_evaluation%log_likelihood - &
      1.0e-10_dp, 'elitist genetic search and local refinement')
   allocate(constrained_model%intercept(1, 2), &
      constrained_model%ar(1, 1, 1, 2))
   allocate(constrained_model%covariance(1, 1, 2), &
      constrained_model%weight(2), constrained_model%degrees_of_freedom(2))
   constrained_model%intercept = reshape([0.1_dp, 0.15_dp], [1, 2])
   constrained_model%ar = reshape([0.3_dp, 0.4_dp], [1, 1, 1, 2])
   constrained_model%covariance = reshape([0.3_dp, 0.6_dp], [1, 1, 2])
   constrained_model%weight = [0.6_dp, 0.4_dp]
   constrained_model%degrees_of_freedom = 0.0_dp
   constrained_model%gaussian_regimes = 2
   allocate(constraints%ar_mapping(2, 1), constraints%mean_group(2), &
      constraints%fixed_weight(2))
   constraints%ar_mapping = 1.0_dp
   constraints%mean_group = [1, 1]
   constraints%fixed_weight = [0.6_dp, 0.4_dp]
   fit = gmvarkit_estimate_constrained(estimation_series, constrained_model, &
      constraints, max_iterations=40, tolerance=1.0e-5_dp)
   moments = gmvarkit_regime_moments(fit%model)
   fitted_evaluation = gmvarkit_evaluate(estimation_series, fit%model)
   call check(fitted_evaluation%info == 0 .and. moments%info == 0 .and. &
      abs(fit%model%ar(1, 1, 1, 1) - fit%model%ar(1, 1, 1, 2)) < &
      1.0e-14_dp .and. abs(moments%mean(1, 1) - moments%mean(1, 2)) < &
      1.0e-12_dp .and. &
      maxval(abs(fit%model%weight - [0.6_dp, 0.4_dp])) < 1.0e-14_dp, &
      'linear AR, shared-mean, and fixed-weight constraints')
   print '(a)', 'gmvarkit tests passed'

contains

   subroutine check(condition, message)
      !! Stop the test program when a logical assertion fails.
      logical, intent(in) :: condition !! Assertion condition.
      character(len=*), intent(in) :: message !! Failure description.

      if (.not. condition) error stop 'FAIL: '//trim(message)
   end subroutine check

end program test_gmvarkit
