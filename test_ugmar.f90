! SPDX-License-Identifier: GPL-3.0-only
! SPDX-FileComment: Numerical tests for algorithms translated from R uGMAR.
program test_ugmar
   use kind_mod, only: dp
   use random_mod, only: set_random_seed
   use ugmar_mod, only: ugmar_model_t, ugmar_evaluation_t, &
      ugmar_stationary_moments_t, ugmar_simulation_t, ugmar_forecast_t, &
      ugmar_quantile_residuals_t, ugmar_residual_tests_t, &
      ugmar_fit_t, ugmar_genetic_fit_t, ugmar_multistart_fit_t, &
      ugmar_constraints_t, ugmar_inference_t, ugmar_hypothesis_test_t, &
      ugmar_likelihood_profile_t, ugmar_regime_conversion_t, &
      ugmar_ar_roots_t, &
      ugmar_model, ugmar_is_stationary, &
      ugmar_random_model, ugmar_smart_model, ugmar_ar_roots, &
      ugmar_regime_means, ugmar_model_from_regime_means, &
      ugmar_evaluate, ugmar_log_likelihood, ugmar_mixing_weights, &
      ugmar_stationary_moments, ugmar_simulate, ugmar_forecast, &
      ugmar_quantile_residuals, ugmar_quantile_residual_tests, &
      ugmar_estimate_constrained, ugmar_genetic_estimate, &
      ugmar_multistart_estimate, ugmar_inference, ugmar_wald_test, &
      ugmar_likelihood_ratio, ugmar_profile_likelihood, &
      ugmar_convert_student_regimes
   implicit none

   type(ugmar_model_t) :: model, unstable, mean_model, high_df_model
   type(ugmar_model_t) :: random_model, smart_model
   type(ugmar_model_t) :: starts(2), gaussian_model
   type(ugmar_evaluation_t) :: evaluation
   type(ugmar_stationary_moments_t) :: moments
   type(ugmar_simulation_t) :: simulation
   type(ugmar_forecast_t) :: forecast
   type(ugmar_quantile_residuals_t) :: residuals
   type(ugmar_residual_tests_t) :: residual_tests
   type(ugmar_fit_t) :: fit, restricted_fit, constrained_fit
   type(ugmar_genetic_fit_t) :: genetic_fit
   type(ugmar_multistart_fit_t) :: multistart_fit
   type(ugmar_constraints_t) :: constraints
   type(ugmar_inference_t) :: inference
   type(ugmar_hypothesis_test_t) :: hypothesis
   type(ugmar_likelihood_profile_t) :: profile
   type(ugmar_regime_conversion_t) :: conversion
   type(ugmar_ar_roots_t) :: roots
   real(dp), allocatable :: weights(:, :)
   real(dp), allocatable :: regime_means(:)
   real(dp) :: series(12), long_series(40)
   integer :: failures, observation

   failures = 0
   series = [0.10_dp, -0.20_dp, 0.05_dp, 0.30_dp, -0.10_dp, 0.15_dp, &
      0.25_dp, -0.05_dp, 0.08_dp, 0.18_dp, -0.12_dp, 0.04_dp]
   model = ugmar_model([0.1_dp, -0.2_dp], reshape([0.4_dp, -0.3_dp], &
      [1, 2]), [0.5_dp, 0.8_dp], [0.6_dp, 0.4_dp], &
      gaussian_regimes=1, degrees_of_freedom=[0.0_dp, 8.0_dp])
   call check(model%info == 0, "construct mixed model", failures)
   call check(ugmar_is_stationary(model), "stationarity", failures)
   roots = ugmar_ar_roots(model)
   call check(roots%info == 0 .and. minval(roots%modulus) > 1.0_dp, &
      "characteristic roots", failures)
   regime_means = ugmar_regime_means(model)
   mean_model = ugmar_model_from_regime_means(model, regime_means)
   call check(mean_model%info == 0 .and. &
      maxval(abs(mean_model%intercept - model%intercept)) < 1.0e-12_dp, &
      "location parametrization round trip", failures)

   moments = ugmar_stationary_moments(model)
   call check(moments%info == 0, "stationary moments status", failures)
   call check(abs(moments%regime_mean(1) - 1.0_dp/6.0_dp) < 1.0e-10_dp, &
      "first regime mean", failures)
   call check(abs(moments%regime_variance(1) - 0.5_dp/0.84_dp) < 1.0e-10_dp, &
      "first regime variance", failures)
   call check(abs(moments%variance - moments%autocovariance(1)) < 1.0e-12_dp, &
      "mixture variance", failures)

   evaluation = ugmar_evaluate(series, model)
   call check(evaluation%info == 0, "likelihood evaluation status", failures)
   call check(size(evaluation%conditional_mean) == 11, &
      "conditional mean extent", failures)
   call check(maxval(abs(sum(evaluation%mixing_weight, dim=2) - 1.0_dp)) < &
      1.0e-12_dp, "mixing weights sum", failures)
   call check(abs(ugmar_log_likelihood(series, model) - &
      evaluation%log_likelihood) < 1.0e-12_dp, "log likelihood adapter", failures)
   weights = ugmar_mixing_weights(series, model)
   call check(all(shape(weights) == [11, 2]), "mixing weight shape", failures)

   residuals = ugmar_quantile_residuals(series, model)
   call check(residuals%info == 0, "quantile residual status", failures)
   call check(size(residuals%residual) == 11, "quantile residual extent", failures)
   call check(all(residuals%cdf > 0.0_dp .and. residuals%cdf < 1.0_dp), &
      "conditional probabilities", failures)
   residual_tests = ugmar_quantile_residual_tests(residuals%residual, &
      [1], [1])
   call check(residual_tests%info == 0, "residual tests status", failures)
   call check(residual_tests%normality%degrees_of_freedom == 3, &
      "residual normality degrees of freedom", failures)

   call set_random_seed(7301)
   simulation = ugmar_simulate(model, series(:1), 8, paths=3)
   call check(simulation%info == 0, "simulation status", failures)
   call check(all(shape(simulation%series) == [8, 3]), &
      "simulation shape", failures)
   call check(all(simulation%regime >= 1 .and. simulation%regime <= 2), &
      "simulated regimes", failures)
   call set_random_seed(7302)
   forecast = ugmar_forecast(model, series(:1), 3, simulations=40, &
      probabilities=[0.1_dp, 0.9_dp])
   call check(forecast%info == 0, "forecast status", failures)
   call check(all(shape(forecast%quantile) == [3, 2]), &
      "forecast quantile shape", failures)
   call check(all(forecast%quantile(:, 1) <= forecast%quantile(:, 2)), &
      "ordered forecast quantiles", failures)

   profile = ugmar_profile_likelihood(series, model, parameter=[1], &
      scale=0.01_dp, points=5)
   call check(profile%info == 0, "profile likelihood status", failures)
   call check(all(shape(profile%log_likelihood) == [5, 1]), &
      "profile likelihood shape", failures)

   starts = [model, model]
   multistart_fit = ugmar_multistart_estimate(series, starts, &
      max_iterations=1)
   call check(allocated(multistart_fit%fit) .and. &
      size(multistart_fit%fit) == 2, "multistart fits", failures)
   call check(multistart_fit%successful_count >= 1, &
      "multistart successful fit", failures)
   call set_random_seed(7303)
   genetic_fit = ugmar_genetic_estimate(series, model, population_size=4, &
      generations=1, mutation_scale=0.01_dp, local_iterations=1)
   call check(genetic_fit%info == 0, "genetic estimate status", failures)
   call check(size(genetic_fit%best_objective) == 1 .and. &
      genetic_fit%evaluations == 8, "genetic search accounting", failures)

   constraints%fixed_weight = model%weight
   constrained_fit = ugmar_estimate_constrained(series, model, constraints, &
      max_iterations=1)
   call check(allocated(constrained_fit%model%weight), &
      "constrained fit model", failures)
   call check(maxval(abs(constrained_fit%model%weight - model%weight)) < &
      1.0e-12_dp, "fixed mixing weights", failures)

   long_series = [(0.15_dp*sin(0.4_dp*real(observation, dp)) + &
      0.01_dp*real(observation, dp), observation=1, 40)]
   gaussian_model = ugmar_model([0.05_dp], reshape([0.3_dp], [1, 1]), &
      [0.4_dp], [1.0_dp])
   fit%model = gaussian_model
   fit%log_likelihood = ugmar_log_likelihood(long_series, gaussian_model)
   fit%parameter_count = 3
   inference = ugmar_inference(long_series, fit)
   call check(inference%info == 0, "numerical inference status", failures)
   if (inference%info == 0) then
      hypothesis = ugmar_wald_test(inference, reshape([1.0_dp, 0.0_dp, &
         0.0_dp], [1, 3]), [inference%parameter(1)])
      call check(hypothesis%info == 0 .and. hypothesis%statistic < 1.0e-10_dp, &
         "Wald test", failures)
   end if
   restricted_fit = fit
   restricted_fit%parameter_count = 2
   restricted_fit%log_likelihood = fit%log_likelihood - 1.0_dp
   hypothesis = ugmar_likelihood_ratio(fit, restricted_fit)
   call check(hypothesis%info == 0 .and. &
      abs(hypothesis%statistic - 2.0_dp) < 1.0e-12_dp, &
      "likelihood-ratio test", failures)

   high_df_model = ugmar_model([0.0_dp, 0.1_dp], &
      reshape([0.2_dp, 0.4_dp], [1, 2]), [0.5_dp, 0.7_dp], &
      [0.45_dp, 0.55_dp], gaussian_regimes=0, &
      degrees_of_freedom=[12.0_dp, 120.0_dp])
   conversion = ugmar_convert_student_regimes(high_df_model, &
      maximum_degrees_of_freedom=100.0_dp, series=series)
   call check(conversion%info == 0, "Student regime conversion status", failures)
   call check(count(conversion%converted) == 1 .and. &
      conversion%fit%model%gaussian_regimes == 1, &
      "Student regime conversion", failures)

   call set_random_seed(7304)
   random_model = ugmar_random_model(2, 3, gaussian_regimes=1, &
      location_mean=0.5_dp, location_standard_deviation=0.2_dp, &
      variance_scale=0.4_dp)
   call check(random_model%info == 0 .and. ugmar_is_stationary(random_model), &
      "random stationary model", failures)
   smart_model = ugmar_smart_model(random_model, accuracy=20.0_dp, &
      random_regimes=[2], location_mean=0.5_dp, &
      location_standard_deviation=0.2_dp, variance_scale=0.4_dp)
   call check(smart_model%info == 0 .and. ugmar_is_stationary(smart_model), &
      "smart model mutation", failures)

   unstable = ugmar_model([0.0_dp], reshape([1.05_dp], [1, 1]), &
      [1.0_dp], [1.0_dp])
   call check(unstable%info == 4, "reject unstable model", failures)
   call check(.not. ugmar_is_stationary(unstable), &
      "unstable stationarity result", failures)

   if (failures /= 0) error stop 1
   print *, "uGMAR tests passed"

contains

   subroutine check(condition, label, failure_count)
      !! Record a failed numerical or structural assertion.
      logical, intent(in) :: condition !! Assertion outcome.
      character(len=*), intent(in) :: label !! Assertion description.
      integer, intent(inout) :: failure_count !! Accumulated failure count.

      if (.not. condition) then
         print *, "FAILED: ", trim(label)
         failure_count = failure_count + 1
      end if
   end subroutine check

end program test_ugmar
