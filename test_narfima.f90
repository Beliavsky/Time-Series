! SPDX-License-Identifier: GPL-3.0-only
! SPDX-FileComment: Regression tests for the R narfima package translation.
program test_narfima
   use kind_mod, only: dp
   use narfima_mod
   use neural_network_mod, only: neural_network_from_parameters, &
      neural_network_predict
   use random_mod, only: set_random_seed
   use utils_mod, only: quiet_nan
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan
   implicit none

   type(narfima_model_t) :: model
   type(narfima_model_t) :: naive_model, arfima_model, arima_model, bsts_model
   type(narfima_model_t) :: incomplete_model
   type(narfima_model_t) :: reference_model
   type(narfima_forecast_t) :: point, distribution
   type(narfima_forecast_t) :: reference_forecast
   real(dp) :: series(96), errors(96), regressors(96, 1)
   real(dp) :: future_regressors(4, 1), innovations(4, 3)
   real(dp) :: incomplete_series(96), incomplete_regressors(96, 1)
   real(dp) :: reference_series(5), reference_errors(5)
   real(dp) :: reference_predictors(4, 2), reference_parameters(7)
   real(dp) :: reference_fitted(4), reference_future(3)
   real(dp), allocatable :: network_fitted(:, :)
   integer :: time, scratch_unit, display_size

   do time = 1, size(series)
      regressors(time, 1) = sin(0.11_dp*real(time, dp))
      errors(time) = 0.08_dp*sin(0.73_dp*real(time, dp))
      series(time) = 20.0_dp + 0.04_dp*real(time, dp) + &
         0.9_dp*sin(2.0_dp*acos(-1.0_dp)*real(time, dp)/12.0_dp) + &
         0.5_dp*regressors(time, 1) + errors(time)
   end do
   call set_random_seed(38171)
   model = narfima_fit(series, errors, p=2, q=2, period=12, &
      seasonal_order=2, hidden_count=2, repetitions=2, direct=.true., &
      response_lambda=0.5_dp, error_lambda=0.5_dp, regressors=regressors, &
      max_iterations=300, tolerance=1.0e-6_dp)
   call check(model%info == 0 .and. all(model%response_lags == [1, 2]) .and. &
      all(model%error_lags == [1, 2, 12, 24]) .and. &
      model%maximum_lag == 24 .and. ieee_is_nan(model%fitted(1)) .and. &
      all(ieee_is_finite(model%fitted(25:))) .and. &
      ieee_is_finite(model%mse) .and. model%members(1)%network%direct, &
      "NARFIMA transformed ensemble fit with joint direct weights")

   call set_random_seed(27183)
   naive_model = narfima_auto_nnaive(series, p=2, q=1, period=1, &
      seasonal_order=0, hidden_count=1, repetitions=1, &
      max_iterations=150, tolerance=1.0e-6_dp)
   call set_random_seed(27183)
   arfima_model = narfima_auto_narfima(series, p=2, q=1, period=1, &
      seasonal_order=0, hidden_count=1, repetitions=1, &
      regressors=regressors, max_iterations=150, tolerance=1.0e-6_dp)
   call set_random_seed(27183)
   arima_model = narfima_auto_narima(series, ar_order=1, &
      difference_order=0, ma_order=0, p=2, q=1, period=1, &
      seasonal_order=0, hidden_count=1, repetitions=1, &
      regressors=regressors, max_iterations=150, tolerance=1.0e-6_dp)
   call set_random_seed(27183)
   bsts_model = narfima_auto_nbsts(series, iterations=8, burn=2, p=2, q=1, &
      period=4, seasonal_order=1, hidden_count=1, repetitions=1, &
      regressors=regressors, max_iterations=150, tolerance=1.0e-6_dp)
   call check(naive_model%info == 0 .and. arfima_model%info == 0 .and. &
      arima_model%info == 0 .and. bsts_model%info == 0, &
      "naive, ARFIMA, ARIMA, and BSTS residual constructors")

   incomplete_series = series
   incomplete_regressors = regressors
   incomplete_series(30) = quiet_nan()
   incomplete_regressors(40, 1) = quiet_nan()
   call set_random_seed(38171)
   incomplete_model = narfima_fit(incomplete_series, errors, p=2, q=2, &
      period=12, seasonal_order=2, hidden_count=1, repetitions=1, &
      regressors=incomplete_regressors, max_iterations=200, &
      tolerance=1.0e-6_dp)
   call check(incomplete_model%info == 0 .and. &
      ieee_is_nan(incomplete_model%fitted(30)) .and. &
      ieee_is_nan(incomplete_model%fitted(31)) .and. &
      ieee_is_nan(incomplete_model%fitted(40)) .and. &
      size(incomplete_model%fitted_indices) < size(series) - &
      incomplete_model%maximum_lag, "complete-case lag fitting and alignment")

   reference_series = [1.0_dp, 1.2_dp, 0.8_dp, 1.5_dp, 1.1_dp]
   reference_errors = [0.1_dp, -0.2_dp, 0.05_dp, 0.3_dp, -0.1_dp]
   reference_parameters = [0.2_dp, 0.4_dp, -0.3_dp, 0.1_dp, 1.2_dp, &
      0.5_dp, -0.25_dp]
   reference_predictors(:, 1) = reference_series(:4)
   reference_predictors(:, 2) = reference_errors(:4)
   reference_fitted = [1.341515810178610_dp, 1.562395027486228_dp, &
      1.235860353957871_dp, 1.579481391770642_dp]
   reference_future = [1.468803791043543_dp, 1.659160580928185_dp, &
      1.773691068740423_dp]
   allocate(reference_model%members(1))
   reference_model%members(1)%network = neural_network_from_parameters( &
      2, 1, 1, reference_parameters, direct=.true.)
   allocate(reference_model%members(1)%direct_coefficients(0))
   network_fitted = neural_network_predict( &
      reference_model%members(1)%network, reference_predictors)
   call check(maxval(abs(network_fitted(:, 1) - reference_fitted)) < &
      1.0e-12_dp, "fixed nnet weights reproduce R fitted values")
   reference_model%series = reference_series
   reference_model%errors = reference_errors
   allocate(reference_model%regressors(5, 0))
   allocate(reference_model%regressor_scales(0))
   reference_model%response_lags = [1]
   reference_model%error_lags = [1]
   reference_model%maximum_lag = 1
   reference_model%hidden_count = 1
   reference_model%repetitions = 1
   reference_model%direct = .true.
   reference_model%transform_response = .false.
   reference_model%transform_errors = .false.
   reference_model%scale_inputs = .false.
   reference_forecast = narfima_forecast_from_innovations(reference_model, &
      reshape([0.0_dp, 0.0_dp, 0.0_dp], [3, 1]))
   call check(reference_forecast%info == 0 .and. &
      maxval(abs(reference_forecast%mean - reference_future)) < 1.0e-12_dp, &
      "fixed nnet weights reproduce R recursive forecasts")

   do time = 1, 4
      future_regressors(time, 1) = sin(0.11_dp*real(96 + time, dp))
   end do
   point = narfima_forecast(model, 4, future_regressors=future_regressors)
   innovations(:, 1) = [-0.10_dp, 0.00_dp, 0.08_dp, -0.03_dp]
   innovations(:, 2) = 0.0_dp
   innovations(:, 3) = -innovations(:, 1)
   distribution = narfima_forecast_from_innovations(model, innovations, &
      future_regressors, [80.0_dp, 95.0_dp])
   call check(point%info == 0 .and. size(point%mean) == 4 .and. &
      all(ieee_is_finite(point%mean)) .and. distribution%info == 0 .and. &
      all(shape(distribution%paths) == [4, 3]) .and. &
      all(shape(distribution%lower) == [4, 2]) .and. &
      all(distribution%lower <= distribution%upper), &
      "recursive point and supplied-innovation forecasts")

   open(newunit=scratch_unit, status="scratch", action="readwrite")
   call display(model, scratch_unit)
   call display(distribution, scratch_unit)
   inquire(unit=scratch_unit, size=display_size)
   close(scratch_unit)
   call check(display_size > 0, "NARFIMA display methods")

   print '(a)', "narfima tests passed"

contains

   subroutine check(condition, message)
      !! Stop the test program when a condition fails.
      logical, intent(in) :: condition !! Test condition.
      character(len=*), intent(in) :: message !! Failure message.

      if (.not. condition) then
         print '(a)', "FAILED: "//message
         error stop 1
      end if
   end subroutine check

end program test_narfima
