! SPDX-License-Identifier: GPL-3.0-only
! SPDX-FileComment: Tests for algorithms translated from the R echos package.
program test_echos
   use kind_mod, only: dp
   use echos_mod, only: echos_ridge_t, echos_model_t, echos_forecast_t, &
      echos_tuning_t, echos_kpss_t, echos_run_reservoir, echos_fit_ridge, echos_kpss, &
      echos_estimate_differences, echos_fit, echos_forecast, echos_tune, display
   use linalg_mod, only: general_eigenvalues
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none

   real(dp) :: inputs(4, 1), input_weights(2, 1), reservoir_weights(2, 2)
   real(dp) :: design(5, 2), response(5), series(60)
   real(dp), allocatable :: states(:, :)
   complex(dp), allocatable :: eigenvalues(:)
   type(echos_ridge_t) :: ridge
   type(echos_model_t) :: model
   type(echos_forecast_t) :: forecast
   type(echos_tuning_t) :: tuning
   type(echos_kpss_t) :: kpss
   integer :: observation, failures, eigen_status

   failures = 0
   inputs(:, 1) = [0.2_dp, -0.1_dp, 0.4_dp, 0.3_dp]
   input_weights(:, 1) = [0.5_dp, -0.25_dp]
   reservoir_weights = reshape([0.2_dp, -0.1_dp, 0.3_dp, 0.15_dp], [2, 2])
   states = echos_run_reservoir(inputs, input_weights, reservoir_weights, 0.8_dp)
   call check(size(states, 1) == 4 .and. size(states, 2) == 2, &
      'reservoir state shape')
   call check(all(abs(states(1, :)) < epsilon(1.0_dp)) .and. &
      all(ieee_is_finite(states)), 'package-compatible reservoir initialization')

   design(:, 1) = 1.0_dp
   design(:, 2) = [-2.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp]
   response = 2.0_dp + 3.0_dp*design(:, 2)
   ridge = echos_fit_ridge(design, response, 0.0_dp)
   call check(ridge%info == 0 .and. &
      maxval(abs(ridge%coefficients - [2.0_dp, 3.0_dp])) < 1.0e-10_dp, &
      'ridge readout coefficients')

   series(1:2) = [0.2_dp, -0.1_dp]
   do observation = 3, size(series)
      series(observation) = 0.72_dp*series(observation - 1) - &
         0.18_dp*series(observation - 2) + &
         0.12_dp*sin(0.43_dp*real(observation, dp))
   end do
   kpss = echos_kpss(series)
   call check(kpss%info == 0, 'KPSS calculation')
   call check(echos_estimate_differences(series) >= 0, &
      'automatic ordinary differencing')

   model = echos_fit(series, lags=[1, 2], difference_order=0, &
      model_count=5, state_count=8, initial_count=2, seed=31415, &
      alpha=0.8_dp, rho=0.75_dp, density=0.4_dp)
   call check(model%info == 0, 'ESN fit status')
   call check(size(model%states, 2) == 8 .and. &
      size(model%output_weights) == 9, 'ESN retained weights and states')
   call general_eigenvalues(model%reservoir_weights, eigenvalues, eigen_status)
   call check(eigen_status == 0 .and. &
      abs(maxval(abs(eigenvalues)) - 0.75_dp) < 1.0e-8_dp .and. &
      count(abs(model%reservoir_weights) > 0.0_dp) == 26, &
      'exact reservoir density and spectral radius')
   call check(count(ieee_is_finite(model%fitted)) == &
      size(series) - 4, 'aligned ESN fitted values')
   call check(model%lambda > 0.0_dp .and. &
      model%degrees_of_freedom > 0.0_dp .and. &
      size(model%candidate_bic) == 5, 'selected ridge readout and metrics')

   forecast = echos_forecast(model, 6, levels=[80.0_dp, 95.0_dp], &
      simulations=30, seed=2718)
   call check(forecast%info == 0 .and. size(forecast%point) == 6, &
      'recursive ESN point forecast')
   call check(size(forecast%actual) == size(series) .and. &
      size(forecast%fitted) == size(series), 'forecast training data retention')
   call check(size(forecast%simulation, 1) == 6 .and. &
      size(forecast%simulation, 2) == 30 .and. &
      size(forecast%interval, 2) == 4, 'bootstrap forecast dimensions')
   call check(all(ieee_is_finite(forecast%point)) .and. &
      all(ieee_is_finite(forecast%interval)) .and. &
      all(forecast%interval(:, 1) <= forecast%interval(:, 3)), &
      'finite ordered forecast intervals')

   tuning = echos_tune(series, horizon=3, split_count=2, &
      alphas=[0.6_dp, 0.9_dp], rhos=[0.7_dp], taus=[0.2_dp], &
      minimum_train=30, lags=[1, 2], difference_order=0, model_count=3, &
      state_count=6, initial_count=1, seed=1618, density=0.5_dp)
   call check(tuning%info == 0 .and. size(tuning%mse) == 4, &
      'expanding-window tuning grid')
   call check(tuning%best >= 1 .and. tuning%best <= 2 .and. &
      all(ieee_is_finite(tuning%mse)), 'tuning selection and errors')
   call display(model)
   call display(forecast)
   call display(tuning)

   if (failures > 0) error stop 'echos tests failed'
   print '(a)', 'echos tests passed'

contains

   subroutine check(condition, label)
      !! Record a failed logical test.
      logical, intent(in) :: condition !! Test condition.
      character(len=*), intent(in) :: label !! Test label.

      if (condition) return
      failures = failures + 1
      print '(a)', 'FAILED: '//trim(label)
   end subroutine check

end program test_echos
