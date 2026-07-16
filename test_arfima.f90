! SPDX-License-Identifier: MIT
! SPDX-FileComment: Regression tests for the arfima translation.
program test_arfima
   use kind_mod, only: dp
   use arfima_mod
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   type(arfima_acvf_t) :: fdwn, fgn, hd, farma, short_memory
   type(arfima_acvf_t) :: invalid_fdwn, invalid_fgn, invalid_hd
   type(arfima_model_t) :: seasonal_model, white_model, integrated_model
   type(arfima_model_t) :: short_model, ar_model, fd_model, fd_initial, fixed_model
   type(arfima_model_t) :: coupled_fixed_model
   type(arfima_simulation_t) :: simulated, random_first, random_second
   type(arfima_likelihood_t) :: exact_value, css_value, fixed_mean_value
   type(arfima_fit_t) :: white_fit, ar_fit, fd_fit, fixed_fit, all_fixed_fit
   type(arfima_fit_t) :: coupled_fixed_fit
   type(arfima_multifit_t) :: mode_fit, random_mode_first, random_mode_second
   type(arfima_multifit_t) :: weeded_modes, best_modes, removed_mode, family_modes
   type(arfima_multifit_t) :: wall_modes, wall_weeded
   type(arfima_multifit_t) :: fixed_mode_fit
   type(arfima_forecast_t) :: white_forecast, ar_forecast, integrated_forecast
   type(arfima_forecast_t) :: fitted_forecast, path_first, path_second
   type(arfima_transfer_t) :: simple_transfer, numerator_transfer, recursive_transfer
   type(arfima_transfer_t) :: delayed_transfer, true_transfer, initial_transfer
   type(arfima_transfer_t) :: empty_transfers(0), transfer_array(1)
   type(arfima_regression_fit_t) :: regression_fit, transfer_fit, fixed_regression_fit
   type(arfima_likelihood_t) :: regression_likelihood
   type(arfima_forecast_t) :: regression_forecast, transfer_forecast
   type(arfima_simulation_t) :: regression_simulation, base_simulation
   type(arfima_identifiability_t) :: identifiable, common_factor, unstable_result
   type(arfima_information_t) :: exact_information, approximate_information
   type(arfima_information_t) :: seasonal_information, fitted_information
   type(arfima_information_t) :: singular_information
   type(arfima_information_t) :: mode_information
   type(arfima_diagnostics_t) :: fit_diagnostics, mode_diagnostics
   type(arfima_diagnostics_t) :: regression_diagnostics, transfer_diagnostics
   type(arfima_covariance_t) :: fit_covariance, fixed_covariance
   type(arfima_transfer_t) :: unstable_transfer_array(1)
   real(dp), allocatable :: weights(:), psi(:)
   real(dp), allocatable :: integrated(:)
   real(dp) :: recovered(3)
   real(dp) :: likelihood_series(20)
   real(dp) :: regression_series(20), regression_x(20, 1), transfer_x(20, 1)
   real(dp) :: empty_regressors(20, 0), empty_transfer_regressors(20, 0)
   real(dp) :: empty_future_regressors(3, 0)
   real(dp) :: future_regressors(3, 1), combined_transfer_regressors(23, 1)
   real(dp) :: empty_future_transfer_regressors(23, 0)
   real(dp), allocatable :: transfer_response(:), expected_response(:)
   real(dp), allocatable :: mode_distances(:, :)
   type(arfima_forecast_t) :: average_forecast
   real(dp) :: mode_starts(1, 3)
   real(dp) :: fixed_mode_starts(2, 2)
   integer :: i
   real(dp), parameter :: pi = acos(-1.0_dp)

   weights = arfima_fractional_weights(0.5_dp, 4)
   call check(size(weights) == 5 .and. &
      maxval(abs(weights - [1.0_dp, -0.5_dp, -0.125_dp, -0.0625_dp, &
      -0.0390625_dp])) < 1.0e-14_dp, 'fractional differencing weights')
   weights = arfima_fractional_weights(0.5_dp, 6, 2)
   call check(maxval(abs(weights - [1.0_dp, 0.0_dp, -0.5_dp, 0.0_dp, &
      -0.125_dp, 0.0_dp, -0.0625_dp])) < 1.0e-14_dp, &
      'seasonal fractional differencing weights')

   psi = arfima_psi_weights([0.5_dp], [0.2_dp], 0.0_dp, 4)
   call check(maxval(abs(psi - [1.0_dp, 0.3_dp, 0.15_dp, 0.075_dp, &
      0.0375_dp])) < 1.0e-14_dp, 'FARMA impulse weights')
   psi = arfima_seasonal_psi_weights([real(dp) ::], [real(dp) ::], 0.0_dp, &
      [0.5_dp], [0.2_dp], 0.0_dp, 4, 8)
   call check(maxval(abs(psi - [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.3_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 0.15_dp])) < 1.0e-14_dp, &
      'seasonal FARMA impulse weights')
   ar_model = arfima_model([0.5_dp], [real(dp) ::], 0.0_dp, &
      long_memory_type=arfima_long_memory_none)
   weights = arfima_pi_weights(ar_model, 3)
   call check(maxval(abs(weights - [1.0_dp, -0.5_dp, 0.0_dp, 0.0_dp])) < &
      1.0e-14_dp, 'AR inverse-filter pi weights')
   ar_model = arfima_model([real(dp) ::], [0.2_dp], 0.0_dp, &
      long_memory_type=arfima_long_memory_none)
   weights = arfima_pi_weights(ar_model, 3)
   call check(maxval(abs(weights - [1.0_dp, 0.2_dp, 0.04_dp, 0.008_dp])) < &
      1.0e-14_dp, 'MA inverse-filter pi weights')
   fd_model = arfima_model([real(dp) ::], [real(dp) ::], 0.25_dp)
   weights = arfima_pi_weights(fd_model, 2)
   call check(maxval(abs(weights - [1.0_dp, -0.25_dp, -0.09375_dp])) < &
      1.0e-14_dp, 'fractional inverse-filter pi weights')

   fdwn = arfima_fdwn_acvf(0.25_dp, 4)
   call check(fdwn%info == 0 .and. &
      abs(fdwn%covariance(0) - exp(log_gamma(0.5_dp) - &
      2.0_dp*log_gamma(0.75_dp))) < 1.0e-14_dp .and. &
      abs(fdwn%covariance(1)/fdwn%covariance(0) - 1.0_dp/3.0_dp) < 1.0e-14_dp, &
      'fractional white-noise autocovariance')

   fgn = arfima_fgn_acvf(0.5_dp, 4)
   call check(fgn%info == 0 .and. abs(fgn%covariance(0) - 1.0_dp) < 1.0e-14_dp .and. &
      maxval(abs(fgn%covariance(1:))) < 1.0e-14_dp, &
      'fractional Gaussian noise at H one half')

   hd = arfima_hd_acvf(2.0_dp, 4)
   call check(hd%info == 0 .and. abs(hd%covariance(0) - 1.0_dp) < 1.0e-14_dp .and. &
      abs(hd%covariance(1) + 3.0_dp/pi**2) < 1.0e-11_dp .and. &
      abs(hd%covariance(2) - hd%covariance(1)/4.0_dp) < 1.0e-13_dp, &
      'hyperbolic-decay autocovariance')

   short_memory = arfima_farma_acvf([0.5_dp], [real(dp) ::], 0.0_dp, 4, 2.0_dp, 32)
   call check(short_memory%info == 0 .and. &
      maxval(abs(short_memory%covariance - [2.0_dp/0.75_dp, 1.0_dp/0.75_dp, &
      0.5_dp/0.75_dp, 0.25_dp/0.75_dp, 0.125_dp/0.75_dp])) < 1.0e-11_dp, &
      'FARMA covariance reuses ARMA covariance')
   farma = arfima_farma_acvf([real(dp) ::], [real(dp) ::], 0.25_dp, 4, 1.5_dp, 32)
   call check(farma%info == 0 .and. &
      maxval(abs(farma%covariance - 1.5_dp*fdwn%covariance)) < 1.0e-12_dp, &
      'pure fractional FARMA covariance')
   farma = arfima_farma_acvf([0.5_dp], [real(dp) ::], 0.25_dp, 4, 1.0_dp)
   call check(farma%info == 0 .and. maxval(abs(farma%covariance - &
      [2.423171_dp, 1.848623_dp, 1.411492_dp, 1.118041_dp, 0.9237375_dp])) < &
      2.0e-6_dp, 'mixed FARMA covariance against arfima reference')

   seasonal_model = arfima_model([0.2_dp], [real(dp) ::], 0.2_dp, &
      seasonal_ar=[0.3_dp], seasonal_d=0.1_dp, period=4)
   farma = arfima_model_acvf(seasonal_model, 6)
   call check(seasonal_model%info == 0 .and. farma%info == 0 .and. &
      maxval(abs(farma%covariance - [1.942011_dp, 1.099388_dp, 0.8465211_dp, &
      0.8520393_dp, 1.151747_dp, 0.7728994_dp, 0.6490468_dp])) < 3.0e-6_dp, &
      'seasonal ARFIMA covariance against arfima reference')

   integrated = arfima_integrate([1.0_dp, 2.0_dp, 3.0_dp], [10.0_dp], 1, 0, 0)
   call check(maxval(abs(integrated - [10.0_dp, 11.0_dp, 13.0_dp, 16.0_dp])) < &
      1.0e-14_dp, 'ordinary inverse differencing')
   integrated = arfima_integrate([1.0_dp, 2.0_dp, 3.0_dp], &
      [10.0_dp, 20.0_dp], 0, 1, 2)
   call check(maxval(abs(integrated - [10.0_dp, 20.0_dp, 11.0_dp, 22.0_dp, &
      14.0_dp])) < 1.0e-14_dp, 'seasonal inverse differencing')
   integrated = arfima_integrate([1.0_dp, 2.0_dp, 3.0_dp], &
      [5.0_dp, 7.0_dp, 9.0_dp], 1, 1, 2)
   do i = 4, 6
      recovered(i - 3) = (integrated(i) - integrated(i - 2)) - &
         (integrated(i - 1) - integrated(i - 3))
   end do
   call check(maxval(abs(recovered - [1.0_dp, 2.0_dp, 3.0_dp])) < 1.0e-14_dp, &
      'combined ordinary and seasonal inverse differencing')

   simulated = arfima_durbin_levinson_simulate([1.0_dp, -2.0_dp, 0.5_dp], &
      [2.0_dp, 0.0_dp, 0.0_dp])
   call check(simulated%info == 0 .and. &
      maxval(abs(simulated%series - sqrt(2.0_dp)*[1.0_dp, -2.0_dp, 0.5_dp])) < &
      1.0e-14_dp .and. maxval(abs(simulated%innovation_variance - 2.0_dp)) < &
      1.0e-14_dp, 'Durbin-Levinson covariance simulation')

   white_model = arfima_model([real(dp) ::], [real(dp) ::], 0.0_dp, &
      innovation_variance=4.0_dp, mean=5.0_dp)
   simulated = arfima_simulate_from_innovations(white_model, &
      [1.0_dp, -1.0_dp, 2.0_dp, -2.0_dp])
   call check(simulated%info == 0 .and. &
      maxval(abs(simulated%stationary_series - [7.0_dp, 3.0_dp, 9.0_dp, &
      1.0_dp])) < 1.0e-14_dp .and. &
      maxval(abs(simulated%series - simulated%stationary_series)) < 1.0e-14_dp, &
      'high-level deterministic ARFIMA simulation')

   integrated_model = arfima_model([real(dp) ::], [real(dp) ::], 0.0_dp, &
      difference_order=1)
   simulated = arfima_simulate_from_innovations(integrated_model, &
      [1.0_dp, -1.0_dp], [10.0_dp])
   call check(simulated%info == 0 .and. &
      maxval(abs(simulated%stationary_series - [1.0_dp, -1.0_dp])) < 1.0e-14_dp .and. &
      maxval(abs(simulated%series - [11.0_dp, 10.0_dp])) < 1.0e-14_dp, &
      'integer integration in ARFIMA simulation')

   random_first = arfima_simulate(white_model, 8, seed=791)
   random_second = arfima_simulate(white_model, 8, seed=791)
   call check(random_first%info == 0 .and. random_second%info == 0 .and. &
      maxval(abs(random_first%series - random_second%series)) < 1.0e-14_dp, &
      'seeded random ARFIMA simulation')

   exact_value = arfima_likelihood([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], white_model)
   css_value = arfima_likelihood([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], white_model, &
      arfima_likelihood_css)
   fixed_mean_value = arfima_likelihood([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], &
      white_model, estimate_mean=.false.)
   call check(exact_value%info == 0 .and. css_value%info == 0 .and. &
      abs(exact_value%mean - 2.5_dp) < 1.0e-14_dp .and. &
      abs(exact_value%innovation_variance - 1.25_dp) < 1.0e-14_dp .and. &
      abs(css_value%log_likelihood - exact_value%log_likelihood) < 1.0e-13_dp .and. &
      abs(fixed_mean_value%mean - 5.0_dp) < 1.0e-14_dp .and. &
      abs(fixed_mean_value%innovation_variance - 7.5_dp) < 1.0e-14_dp, &
      'exact and conditional profiled Gaussian likelihoods')

   short_model = arfima_model([real(dp) ::], [real(dp) ::], 0.0_dp, &
      long_memory_type=arfima_long_memory_none)
   white_fit = arfima_fit([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], short_model)
   call check(white_fit%info == 0 .and. white_fit%converged .and. &
      white_fit%parameter_count == 0 .and. abs(white_fit%model%mean - 2.5_dp) < &
      1.0e-14_dp .and. ieee_is_finite(white_fit%aic) .and. &
      ieee_is_finite(white_fit%aicc) .and. ieee_is_finite(white_fit%bic), &
      'closed-form white-noise ARFIMA fit')

   likelihood_series(1) = 0.4_dp
   do i = 2, size(likelihood_series)
      likelihood_series(i) = 0.55_dp*likelihood_series(i - 1) + &
         sin(1.3_dp*real(i, dp))
   end do
   ar_model = arfima_model([0.1_dp], [real(dp) ::], 0.0_dp, &
      long_memory_type=arfima_long_memory_none)
   ar_fit = arfima_fit(likelihood_series, ar_model, max_iterations=100, tolerance=1.0e-5_dp)
   call check(ar_fit%likelihood%info == 0 .and. ar_fit%converged .and. &
      ar_fit%parameter_count == 1 .and. &
      abs(ar_fit%model%ar(1)) < 1.0_dp .and. &
      ieee_is_finite(ar_fit%likelihood%log_likelihood) .and. &
      size(ar_fit%parameter_covariance, 1) == 1 .and. &
      size(ar_fit%standard_error) == 1, 'BFGS ARFIMA fit and inference')
   fixed_model = arfima_model([0.2_dp], [0.1_dp], 0.0_dp, &
      long_memory_type=arfima_long_memory_none)
   fixed_fit = arfima_fit_fixed(likelihood_series, fixed_model, [.true., .false.], &
      max_iterations=100, tolerance=1.0e-5_dp)
   all_fixed_fit = arfima_fit_fixed(likelihood_series, fixed_model, [.true., .true.])
   call check(fixed_fit%likelihood%info == 0 .and. fixed_fit%converged .and. &
      fixed_fit%parameter_count == 2 .and. fixed_fit%estimated_parameter_count == 1 .and. &
      abs(fixed_fit%model%ar(1) - 0.2_dp) < 1.0e-14_dp .and. &
      fixed_fit%standard_error(1) == 0.0_dp .and. &
      maxval(abs(fixed_fit%parameter_covariance(1, :))) == 0.0_dp .and. &
      maxval(abs(fixed_fit%parameter_covariance(:, 1))) == 0.0_dp .and. &
      all_fixed_fit%converged .and. all_fixed_fit%estimated_parameter_count == 0 .and. &
      maxval(abs(all_fixed_fit%model%ar - fixed_model%ar)) < 1.0e-14_dp .and. &
      maxval(abs(all_fixed_fit%model%theta - fixed_model%theta)) < 1.0e-14_dp, &
      'partial and fully fixed ARFIMA likelihood fits')
   fixed_covariance = arfima_fit_covariance(fixed_fit, resolution=1024)
   call check(fixed_covariance%expected_available .and. &
      maxval(abs(fixed_covariance%observed(1, :))) == 0.0_dp .and. &
      maxval(abs(fixed_covariance%expected(1, :))) == 0.0_dp .and. &
      fixed_covariance%observed(2, 2) > 0.0_dp .and. &
      fixed_covariance%expected(2, 2) > 0.0_dp, &
      'zero-padded covariance for fixed ARFIMA parameters')
   coupled_fixed_model = arfima_model([0.2_dp, 0.1_dp], [real(dp) ::], 0.0_dp, &
      long_memory_type=arfima_long_memory_none)
   coupled_fixed_fit = arfima_fit_fixed(likelihood_series, coupled_fixed_model, &
      [.true., .false.], max_iterations=100, tolerance=1.0e-5_dp)
   call check(coupled_fixed_fit%likelihood%info == 0 .and. &
      coupled_fixed_fit%converged .and. &
      abs(coupled_fixed_fit%model%ar(1) - 0.2_dp) < 1.0e-14_dp, &
      'fixed operating coefficient under coupled PACF mapping')
   fixed_mode_starts = reshape([atanh(0.2_dp), atanh(-0.7_dp), &
      atanh(0.8_dp), atanh(0.7_dp)], [2, 2])
   fixed_mode_fit = arfima_fit_modes_fixed(likelihood_series, fixed_model, &
      [.true., .false.], fixed_mode_starts, max_iterations=100, tolerance=1.0e-5_dp)
   call check(fixed_mode_fit%info == 0 .and. size(fixed_mode_fit%modes) == 2 .and. &
      abs(fixed_mode_fit%modes(1)%model%ar(1) - 0.2_dp) < 1.0e-14_dp .and. &
      abs(fixed_mode_fit%modes(2)%model%ar(1) - 0.2_dp) < 1.0e-14_dp, &
      'multimode fitting with fixed ARFIMA parameters')
   fit_diagnostics = arfima_diagnostics(likelihood_series, ar_fit, max_lag=5)
   fit_covariance = arfima_fit_covariance(ar_fit, resolution=1024)
   call check(fit_diagnostics%info == 0 .and. fit_diagnostics%observations == 20 .and. &
      maxval(abs(fit_diagnostics%fitted + fit_diagnostics%residuals - &
      likelihood_series)) < 1.0e-13_dp .and. &
      abs(fit_diagnostics%residual_acf(1) - 1.0_dp) < 1.0e-14_dp .and. &
      all(fit_diagnostics%ljung_box_p_value >= 0.0_dp) .and. &
      all(fit_diagnostics%ljung_box_p_value <= 1.0_dp) .and. &
      all(fit_diagnostics%squared_ljung_box_p_value >= 0.0_dp) .and. &
      all(fit_diagnostics%squared_ljung_box_p_value <= 1.0_dp) .and. &
      all(fit_diagnostics%qq_sample(2:) >= fit_diagnostics%qq_sample(:19)), &
      'ARFIMA fitted values and residual diagnostics')
   call check(fit_covariance%info == 0 .and. fit_covariance%expected_available .and. &
      all(shape(fit_covariance%observed) == [1, 1]) .and. &
      all(shape(fit_covariance%expected) == [1, 1]) .and. &
      abs(fit_covariance%observed_correlation(1, 1) - 1.0_dp) < 1.0e-14_dp .and. &
      abs(fit_covariance%expected_correlation(1, 1) - 1.0_dp) < 1.0e-14_dp, &
      'observed and expected ARFIMA covariance extraction')

   mode_starts = reshape([atanh(0.1_dp), atanh(0.8_dp), atanh(-0.7_dp)], [1, 3])
   mode_fit = arfima_fit_modes(likelihood_series, ar_model, mode_starts, &
      max_iterations=100, tolerance=1.0e-5_dp)
   call check(mode_fit%info == 0 .and. mode_fit%attempted == 3 .and. &
      mode_fit%converged == 3 .and. size(mode_fit%modes) == 3 .and. &
      abs(sum(mode_fit%weight) - 1.0_dp) < 1.0e-14_dp .and. &
      all(mode_fit%starting_index == [1, 2, 3]), &
      'independent ARFIMA fits from supplied starts')
   random_mode_first = arfima_multistart_fit(likelihood_series, ar_model, 3, &
      seed=481, max_iterations=100, tolerance=1.0e-5_dp)
   random_mode_second = arfima_multistart_fit(likelihood_series, ar_model, 3, &
      seed=481, max_iterations=100, tolerance=1.0e-5_dp)
   call check(all(shape(random_mode_first%starting_parameters) == [1, 3]) .and. &
      maxval(abs(random_mode_first%starting_parameters - &
      random_mode_second%starting_parameters)) < 1.0e-14_dp, &
      'seeded stable random ARFIMA starts')
   mode_distances = arfima_mode_distances(mode_fit, transformed=.true.)
   call check(maxval(abs(mode_distances - transpose(mode_distances))) < 1.0e-14_dp .and. &
      maxval(abs([(mode_distances(i, i), i=1, 3)])) < 1.0e-14_dp, &
      'symmetric ARFIMA mode distances in PACF space')
   weeded_modes = arfima_weed_modes(mode_fit, arfima_weed_both, 1.0e-3_dp, &
      adaptive=.false.)
   best_modes = arfima_best_modes(mode_fit, 2)
   removed_mode = arfima_remove_mode(best_modes, 1)
   call check(weeded_modes%weeded .and. size(weeded_modes%modes) == 1 .and. &
      size(best_modes%modes) == 2 .and. &
      best_modes%modes(1)%likelihood%log_likelihood >= &
      best_modes%modes(2)%likelihood%log_likelihood .and. &
      size(removed_mode%modes) == 1 .and. &
      abs(sum(removed_mode%weight) - 1.0_dp) < 1.0e-14_dp, &
      'weed, rank, and remove ARFIMA likelihood modes')
   wall_modes = mode_fit
   wall_modes%modes(1)%model%ar(1) = 0.999_dp
   wall_modes%modes(2)%model%ar(1) = 0.995_dp
   wall_modes%modes(3)%model%ar(1) = -0.5_dp
   wall_weeded = arfima_weed_modes(wall_modes, arfima_weed_none, walls=.true., &
      wall_tolerance=0.01_dp)
   call check(wall_weeded%weeded .and. size(wall_weeded%modes) == 2, &
      'boundary-aware ARFIMA mode weeding')
   fitted_forecast = arfima_mode_forecast(likelihood_series, best_modes, 1, 2)
   average_forecast = arfima_average_forecast(likelihood_series, best_modes, 2)
   mode_information = arfima_mode_information(best_modes, 1, resolution=1024)
   identifiable = arfima_mode_identifiability(best_modes, 1)
   mode_diagnostics = arfima_mode_diagnostics(likelihood_series, best_modes, 1, 5)
   call check(fitted_forecast%info == 0 .and. average_forecast%info == 0 .and. &
      mode_information%info == 0 .and. mode_diagnostics%info == 0 .and. &
      identifiable%info == 0 .and. identifiable%identifiable, &
      'mode-specific ARFIMA forecast and diagnostics')
   call check(maxval(abs(average_forecast%covariance - &
      transpose(average_forecast%covariance))) < 1.0e-13_dp .and. &
      all(average_forecast%standard_error > 0.0_dp) .and. &
      all(average_forecast%lower < average_forecast%mean) .and. &
      all(average_forecast%upper > average_forecast%mean), &
      'likelihood-weighted multimode ARFIMA forecast')
   family_modes = arfima_select_long_memory(likelihood_series, ar_model, 1, &
      seed=712, max_iterations=100, tolerance=1.0e-5_dp)
   call check(family_modes%attempted == 4 .and. size(family_modes%modes) >= 1 .and. &
      family_modes%weeded .and. abs(sum(family_modes%weight) - 1.0_dp) < 1.0e-14_dp .and. &
      all(family_modes%long_memory_type >= arfima_long_memory_none) .and. &
      all(family_modes%long_memory_type <= arfima_long_memory_hd), &
      'selection across ARFIMA long-memory families')

   do i = 1, size(likelihood_series)
      likelihood_series(i) = sin(1.7_dp*real(i, dp)) + &
         0.3_dp*cos(0.41_dp*real(i, dp))
   end do
   fd_model = arfima_model([real(dp) ::], [real(dp) ::], 0.2_dp)
   simulated = arfima_simulate_from_innovations(fd_model, likelihood_series)
   fd_initial = arfima_model([real(dp) ::], [real(dp) ::], 0.05_dp)
   fd_fit = arfima_fit(simulated%series, fd_initial, max_iterations=100, &
      tolerance=1.0e-5_dp)
   call check(fd_fit%likelihood%info == 0 .and. fd_fit%converged .and. &
      fd_fit%parameter_count == 1 .and. &
      fd_fit%model%long_memory_parameter > -1.0_dp .and. &
      fd_fit%model%long_memory_parameter < 0.5_dp, &
      'bounded fractional-difference likelihood fit')

   white_forecast = arfima_forecast([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], &
      white_model, 3, 0.9_dp)
   call check(white_forecast%info == 0 .and. white_forecast%horizon == 3 .and. &
      maxval(abs(white_forecast%mean - 5.0_dp)) < 1.0e-13_dp .and. &
      maxval(abs(white_forecast%standard_error - 2.0_dp)) < 1.0e-13_dp .and. &
      maxval(abs(white_forecast%covariance - &
      reshape([4.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 4.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 4.0_dp], [3, 3]))) < 1.0e-12_dp .and. &
      all(white_forecast%lower < white_forecast%mean) .and. &
      all(white_forecast%upper > white_forecast%mean), &
      'exact white-noise forecasts and intervals')

   ar_model = arfima_model([0.5_dp], [real(dp) ::], 0.0_dp, &
      innovation_variance=2.0_dp, long_memory_type=arfima_long_memory_none)
   ar_forecast = arfima_forecast([1.0_dp, 0.5_dp, 0.25_dp], ar_model, 2)
   call check(ar_forecast%info == 0 .and. &
      maxval(abs(ar_forecast%mean - [0.125_dp, 0.0625_dp])) < 1.0e-12_dp .and. &
      maxval(abs(ar_forecast%covariance - &
      reshape([2.0_dp, 1.0_dp, 1.0_dp, 2.5_dp], [2, 2]))) < 1.0e-11_dp, &
      'exact AR forecast means and cross-horizon covariance')

   integrated_model = arfima_model([real(dp) ::], [real(dp) ::], 0.0_dp, &
      difference_order=1, long_memory_type=arfima_long_memory_none)
   integrated_forecast = arfima_forecast([10.0_dp, 12.0_dp, 15.0_dp], &
      integrated_model, 3)
   call check(integrated_forecast%info == 0 .and. &
      maxval(abs(integrated_forecast%mean - 15.0_dp)) < 1.0e-12_dp .and. &
      maxval(abs(integrated_forecast%covariance - &
      reshape([1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 2.0_dp, 2.0_dp, &
      1.0_dp, 2.0_dp, 3.0_dp], [3, 3]))) < 1.0e-12_dp, &
      'integrated forecasts and accumulated uncertainty')

   fitted_forecast = arfima_forecast([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], &
      white_fit, 2)
   call check(fitted_forecast%info == 0 .and. &
      maxval(abs(fitted_forecast%mean - 2.5_dp)) < 1.0e-13_dp, &
      'forecast directly from ARFIMA fit')
   path_first = arfima_forecast_paths([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], &
      white_fit, 3, 10, seed=917)
   path_second = arfima_forecast_paths([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], &
      white_fit, 3, 10, seed=917)
   call check(path_first%info == 0 .and. path_second%info == 0 .and. &
      all(shape(path_first%paths) == [3, 10]) .and. &
      maxval(abs(path_first%paths - path_second%paths)) < 1.0e-14_dp, &
      'seeded conditional ARFIMA forecast paths')

   simple_transfer = arfima_transfer([real(dp) ::], [2.0_dp])
   numerator_transfer = arfima_transfer([real(dp) ::], [2.0_dp, 0.5_dp])
   recursive_transfer = arfima_transfer([0.5_dp], [2.0_dp])
   delayed_transfer = arfima_transfer([real(dp) ::], [2.0_dp], delay=1)
   transfer_response = arfima_transfer_response([1.0_dp, 2.0_dp, 3.0_dp], &
      simple_transfer)
   call check(maxval(abs(transfer_response - [2.0_dp, 4.0_dp, 6.0_dp])) < &
      1.0e-14_dp, 'static transfer-function response')
   transfer_response = arfima_transfer_response([1.0_dp, 2.0_dp, 3.0_dp], &
      numerator_transfer)
   call check(maxval(abs(transfer_response - [0.0_dp, 3.5_dp, 5.0_dp])) < &
      1.0e-14_dp, 'transfer-function numerator signs')
   transfer_response = arfima_transfer_response([1.0_dp, 2.0_dp, 3.0_dp], &
      recursive_transfer)
   call check(maxval(abs(transfer_response - [0.0_dp, 4.0_dp, 8.0_dp])) < &
      1.0e-14_dp, 'recursive transfer-function denominator')
   transfer_response = arfima_transfer_response([1.0_dp, 2.0_dp, 3.0_dp], &
      delayed_transfer)
   call check(maxval(abs(transfer_response - [0.0_dp, 2.0_dp, 4.0_dp])) < &
      1.0e-14_dp, 'delayed transfer-function response')

   do i = 1, 20
      regression_x(i, 1) = sin(0.31_dp*real(i, dp))
      transfer_x(i, 1) = regression_x(i, 1)
      regression_series(i) = 1.5_dp*regression_x(i, 1) + &
         0.2_dp*cos(1.7_dp*real(i, dp))
   end do
   regression_likelihood = arfima_regression_likelihood(regression_series, short_model, &
      regression_x, [1.5_dp], empty_transfers, empty_transfer_regressors)
   call check(regression_likelihood%info == 0 .and. &
      regression_likelihood%innovation_variance < 0.1_dp, &
      'known static-regression ARFIMA likelihood')
   regression_fit = arfima_regression_fit(regression_series, short_model, regression_x, &
      [0.0_dp], empty_transfers, empty_transfer_regressors, max_iterations=100, &
      tolerance=1.0e-6_dp)
   fixed_model = arfima_model([0.1_dp], [real(dp) ::], 0.0_dp, &
      long_memory_type=arfima_long_memory_none)
   fixed_regression_fit = arfima_regression_fit_fixed(regression_series, fixed_model, &
      regression_x, [1.5_dp], empty_transfers, empty_transfer_regressors, &
      [.false., .true.], max_iterations=100, tolerance=1.0e-6_dp)
   call check(regression_fit%likelihood%info == 0 .and. regression_fit%converged .and. &
      abs(regression_fit%coefficients(1) - 1.5_dp) < 0.1_dp .and. &
      size(regression_fit%parameter_covariance, 1) == 1, &
      'joint static-regression ARFIMA fit')
   call check(fixed_regression_fit%likelihood%info == 0 .and. &
      fixed_regression_fit%converged .and. &
      fixed_regression_fit%parameter_count == 2 .and. &
      fixed_regression_fit%estimated_parameter_count == 1 .and. &
      abs(fixed_regression_fit%coefficients(1) - 1.5_dp) < 1.0e-14_dp .and. &
      fixed_regression_fit%standard_error(2) == 0.0_dp .and. &
      maxval(abs(fixed_regression_fit%parameter_covariance(2, :))) == 0.0_dp .and. &
      fixed_regression_fit%standard_error(1) > 0.0_dp, &
      'fixed regression ARFIMA coefficient')
   regression_diagnostics = arfima_diagnostics(regression_series, regression_fit, &
      regression_x, empty_transfer_regressors, max_lag=5)
   call check(regression_diagnostics%info == 0 .and. &
      maxval(abs(regression_diagnostics%fitted + &
      regression_diagnostics%residuals - regression_series)) < 1.0e-13_dp .and. &
      maxval(abs(regression_diagnostics%regression_residuals - regression_series + &
      regression_fit%coefficients(1)*regression_x(:, 1))) < 1.0e-13_dp, &
      'regression ARFIMA fitted values and residual diagnostics')

   true_transfer = arfima_transfer([real(dp) ::], [1.7_dp])
   initial_transfer = arfima_transfer([real(dp) ::], [0.5_dp])
   transfer_response = arfima_transfer_response(transfer_x(:, 1), true_transfer)
   regression_series = transfer_response + &
      [(0.15_dp*cos(1.9_dp*real(i, dp)), i=1, 20)]
   transfer_array(1) = initial_transfer
   transfer_fit = arfima_regression_fit(regression_series, short_model, empty_regressors, &
      [real(dp) ::], transfer_array, transfer_x, max_iterations=100, &
      tolerance=1.0e-6_dp)
   call check(transfer_fit%likelihood%info == 0 .and. transfer_fit%converged .and. &
      abs(transfer_fit%transfers(1)%numerator(1) - 1.7_dp) < 0.1_dp, &
      'joint dynamic-regression transfer-function fit')
   transfer_diagnostics = arfima_diagnostics(regression_series, transfer_fit, &
      empty_regressors, transfer_x, max_lag=5)
   expected_response = arfima_transfer_response(transfer_x(:, 1), &
      transfer_fit%transfers(1))
   call check(transfer_diagnostics%info == 0 .and. &
      maxval(abs(transfer_diagnostics%regression_residuals - regression_series + &
      expected_response)) < 1.0e-13_dp, &
      'transfer-function ARFIMA residual diagnostics')

   future_regressors(:, 1) = [0.2_dp, -0.4_dp, 0.7_dp]
   regression_forecast = arfima_regression_forecast(regression_series, regression_fit, &
      regression_x, future_regressors, empty_future_transfer_regressors, 3)
   call check(regression_forecast%info == 0 .and. &
      maxval(abs(regression_forecast%mean - regression_fit%model%mean - &
      regression_fit%coefficients(1)*future_regressors(:, 1))) < 1.0e-12_dp, &
      'regression ARFIMA forecast with future regressors')

   combined_transfer_regressors(:20, 1) = transfer_x(:, 1)
   combined_transfer_regressors(21:, 1) = [0.2_dp, -0.4_dp, 0.7_dp]
   transfer_forecast = arfima_regression_forecast(regression_series, transfer_fit, &
      empty_regressors, empty_future_regressors, &
      combined_transfer_regressors, 3)
   expected_response = arfima_transfer_response(combined_transfer_regressors(:, 1), &
      transfer_fit%transfers(1))
   call check(transfer_forecast%info == 0 .and. &
      maxval(abs(transfer_forecast%regression_mean - expected_response(21:23))) < &
      1.0e-12_dp, 'transfer-function ARFIMA forecast')

   regression_simulation = arfima_regression_simulate_from_innovations(regression_fit, &
      likelihood_series, regression_x, empty_transfer_regressors)
   base_simulation = arfima_simulate_from_innovations(regression_fit%model, &
      likelihood_series)
   call check(regression_simulation%info == 0 .and. base_simulation%info == 0 .and. &
      maxval(abs(regression_simulation%series - base_simulation%series - &
      regression_fit%coefficients(1)*regression_x(:, 1))) < 1.0e-12_dp, &
      'simulation from fitted regression ARFIMA model')

   ar_model = arfima_model([0.3_dp], [0.2_dp], 0.0_dp, &
      long_memory_type=arfima_long_memory_none)
   identifiable = arfima_identifiability(ar_model)
   ar_model = arfima_model([0.3_dp], [0.3_dp], 0.0_dp, &
      long_memory_type=arfima_long_memory_none)
   common_factor = arfima_identifiability(ar_model)
   unstable_transfer_array(1) = arfima_transfer([1.2_dp], [1.0_dp])
   unstable_result = arfima_identifiability(short_model, unstable_transfer_array)
   call check(identifiable%info == 0 .and. identifiable%identifiable .and. &
      identifiable%minimum_root_distance > 1.0_dp .and. &
      .not. common_factor%identifiable .and. &
      common_factor%minimum_root_distance < 1.0e-10_dp .and. &
      .not. unstable_result%transfer_stable .and. .not. unstable_result%identifiable, &
      'ARFIMA stability and common-factor identifiability')

   fd_model = arfima_model([real(dp) ::], [real(dp) ::], 0.2_dp)
   exact_information = arfima_fisher_information(fd_model, resolution=8192)
   approximate_information = arfima_fisher_information(fd_model, exact=.false., &
      resolution=10000)
   call check(exact_information%info == 0 .and. exact_information%positive_definite .and. &
      abs(exact_information%information(1, 1) - pi**2/6.0_dp) < 2.0e-3_dp .and. &
      abs(approximate_information%information(1, 1) - pi**2/6.0_dp) < &
      1.1e-4_dp, 'exact and truncated fractional Fisher information')

   ar_model = arfima_model([0.5_dp], [real(dp) ::], 0.0_dp, &
      long_memory_type=arfima_long_memory_none)
   exact_information = arfima_fisher_information(ar_model, resolution=4096)
   approximate_information = arfima_fisher_information(ar_model, exact=.false., &
      resolution=4096)
   call check(abs(exact_information%information(1, 1) - 4.0_dp/3.0_dp) < &
      1.0e-10_dp .and. abs(approximate_information%information(1, 1) - &
      4.0_dp/3.0_dp) < 1.0e-10_dp, 'AR(1) Fisher information')

   seasonal_model = arfima_model([real(dp) ::], [real(dp) ::], 0.1_dp, &
      seasonal_d=0.2_dp, period=4)
   seasonal_information = arfima_fisher_information(seasonal_model, resolution=16384)
   call check(seasonal_information%info == 0 .and. &
      abs(seasonal_information%information(1, 1) - pi**2/6.0_dp) < 2.0e-3_dp .and. &
      abs(seasonal_information%information(2, 2) - pi**2/6.0_dp) < 8.0e-3_dp .and. &
      abs(seasonal_information%information(1, 2) - pi**2/24.0_dp) < 3.0e-3_dp, &
      'ordinary and seasonal fractional Fisher information')

   fitted_information = arfima_fit_information(fd_fit, exact=.false., resolution=4096)
   call check(fitted_information%info == 0 .and. fitted_information%observations == 20 .and. &
      fitted_information%positive_definite .and. &
      fitted_information%covariance(1, 1) > 0.0_dp, &
      'information covariance from fitted ARFIMA model')
   ar_model = arfima_model([0.3_dp], [0.3_dp], 0.0_dp, &
      long_memory_type=arfima_long_memory_none)
   singular_information = arfima_fisher_information(ar_model, resolution=4096)
   call check(.not. singular_information%positive_definite, &
      'common ARMA factor gives singular Fisher information')

   invalid_fdwn = arfima_fdwn_acvf(0.5_dp, 3)
   invalid_fgn = arfima_fgn_acvf(1.0_dp, 3)
   invalid_hd = arfima_hd_acvf(0.0_dp, 3)
   call check(invalid_fdwn%info /= 0 .and. invalid_fgn%info /= 0 .and. &
      invalid_hd%info /= 0, 'long-memory parameter validation')

contains

   subroutine check(condition, label)
      ! Stop the test program when one assertion fails.
      logical, intent(in) :: condition
      character(*), intent(in) :: label

      if (.not. condition) then
         write (*, '(a)') 'FAILED: '//label
         error stop 1
      end if
   end subroutine check

end program test_arfima
