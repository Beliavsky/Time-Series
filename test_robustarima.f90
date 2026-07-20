! SPDX-License-Identifier: BSD-3-Clause
! SPDX-FileComment: Regression tests for the robustarima translation.
program test_robustarima
   use kind_mod, only: dp
   use robustarima_mod
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none

   type(robustarima_filter_t) :: filtered
   type(robustarima_fit_t) :: fit
   type(robustarima_forecast_t) :: forecast
   type(robustarima_outliers_t) :: outliers
   type(robustarima_order_selection_t) :: selection
   type(robustarima_tau_inference_t) :: inference
   real(dp), allocatable :: series(:), regressors(:, :), pseudo(:)
   real(dp) :: innovation
   integer :: i

   call check(abs(robustarima_rho(1.5_dp) - 1.125_dp) < 1.0e-14_dp, &
      'rho quadratic branch')
   call check(abs(robustarima_rho(4.0_dp) - 3.25_dp) < 1.0e-14_dp, &
      'rho bounded branch')
   call check(abs(robustarima_psi(1.5_dp) - 1.5_dp) < 1.0e-14_dp .and. &
      robustarima_psi(4.0_dp) == 0.0_dp, 'psi branches')
   call check(abs(robustarima_psi_derivative(1.5_dp) - 1.0_dp) < &
      1.0e-14_dp, 'psi derivative')
   call check(robustarima_m_scale([-2.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, &
      2.0_dp]) > 0.0_dp, 'M-scale positivity')
   call check(robustarima_tau_scale([-2.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, &
      2.0_dp]) > 0.0_dp, 'tau-scale positivity')
   pseudo = robustarima_correlation_series([0.0_dp, 1.0_dp, 100.0_dp, &
      -100.0_dp])
   call check(size(pseudo) == 4 .and. maxval(abs(pseudo)) <= 2.5_dp, &
      'robust correlation pseudo-series')

   filtered = robustarima_bounded_filter([0.0_dp, 1.0_dp, 20.0_dp, &
      0.5_dp], [0.5_dp], [real(dp) ::], tuning_constant=1.0_dp, &
      scale=1.0_dp)
   call check(filtered%info == 0, 'bounded filter status')
   call check(filtered%bounded_innovations(3) == 0.0_dp .and. &
      abs(filtered%cleaned(3) - 0.5_dp) < 1.0e-14_dp, &
      'bounded filter rejects gross innovation')
   call check(all(shape(filtered%state) == [1, 4]) .and. &
      all(filtered%prediction_scale > 0.0_dp) .and. &
      all(ieee_is_finite(filtered%state_covariance)), &
      'bounded state filter moments')
   filtered = robustarima_bounded_filter([0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 6.0_dp, 6.0_dp, 6.0_dp, 6.0_dp], [0.9_dp], &
      [real(dp) ::], tuning_constant=1.0_dp, scale=1.0_dp)
   call check(filtered%info == 0 .and. filtered%level_shift_candidate(6) .and. &
      filtered%weights(6) == 1.0_dp .and. &
      abs(filtered%cleaned(6) - 6.0_dp) < 1.0e-12_dp, &
      'bounded filter level-shift rewind')

   allocate(series(120), regressors(120, 1))
   regressors = 1.0_dp
   series(1) = 2.0_dp
   do i = 2, size(series)
      innovation = 0.15_dp*sin(0.73_dp*real(i, dp)) + &
         0.08_dp*cos(1.37_dp*real(i, dp))
      series(i) = 2.0_dp + 0.6_dp*(series(i - 1) - 2.0_dp) + innovation
   end do
   series(65) = series(65) + 12.0_dp
   fit = robustarima_fit(series, 1, 0, regressors=regressors, &
      tuning_constant=1.0_dp, max_iterations=300)
   call check(fit%info == 0 .and. fit%converged, 'robust AR fit status')
   call check(abs(fit%ar(1) - 0.6_dp) < 0.2_dp .and. &
      abs(fit%regression_coefficients(1) - 2.0_dp) < 0.2_dp, &
      'robust AR and regression estimates')
   call check(fit%filter%weights(65) < 0.1_dp .and. &
      fit%filter%innovation_scale > 0.0_dp, 'robust AR outlier weight')
   call check(all(ieee_is_finite(fit%standard_error)), &
      'robust parameter inference')
   inference = robustarima_tau_inference(fit, regressors)
   call check(inference%info == 0 .and. &
      all(shape(inference%filtered_design) == [120, 1]) .and. &
      size(inference%implicit_weight) > 0, 'tau filtered regression design')
   call check(inference%inverse_efficiency > 0.0_dp .and. &
      inference%covariance(1, 1) > 0.0_dp .and. &
      all(ieee_is_finite(inference%covariance)), &
      'tau sandwich covariance')
   call check(maxval(abs(inference%covariance - &
      transpose(inference%covariance))) < 1.0e-12_dp .and. &
      abs(fit%regression_covariance(1, 1) - &
      inference%covariance(1, 1)) < 1.0e-12_dp .and. &
      fit%regression_standard_error(1) > 0.0_dp .and. &
      fit%tau_inverse_efficiency > 0.0_dp, &
      'tau inference stored in fit')
   call check(abs(fit%parameter_covariance(2, 2) - &
      inference%covariance(1, 1)) < 1.0e-12_dp .and. &
      abs(fit%standard_error(2) - inference%standard_error(1)) < &
      1.0e-12_dp, 'tau covariance replaces regression Hessian block')

   outliers = robustarima_detect_outliers(fit, series, 3.0_dp)
   call check(outliers%info == 0 .and. any(outliers%position == 65), &
      'robust additive outlier detection')
   call check(abs(outliers%cleaned_series(65) - &
      (series(65) - 12.0_dp)) < 1.0_dp, 'robust outlier cleaning')

   forecast = robustarima_forecast(fit, series, 3, regressors, &
      reshape([1.0_dp, 1.0_dp, 1.0_dp], [3, 1]), [0.95_dp])
   call check(forecast%info == 0 .and. size(forecast%mean) == 3 .and. &
      all(forecast%standard_error > 0.0_dp), 'robust forecast status')
   call check(all(forecast%lower(:, 1) < forecast%mean) .and. &
      all(forecast%upper(:, 1) > forecast%mean), &
      'robust forecast intervals')

   selection = robustarima_select_ar_order(series, 2, &
      regressors=regressors, tuning_constant=1.0_dp, max_iterations=300)
   call check(selection%info == 0 .and. selection%selected_order >= 1 .and. &
      selection%criterion(1) < selection%criterion(0), &
      'robust automatic AR selection')

   deallocate(series, regressors)
   allocate(series(40))
   series(:4) = [0.2_dp, -0.1_dp, 0.05_dp, 0.15_dp]
   do i = 5, size(series)
      series(i) = series(i - 4) + 0.03_dp*sin(0.91_dp*real(i, dp))
   end do
   fit = robustarima_fit(series, 0, 0, seasonal_period=4, &
      seasonal_difference_order=1)
   call check(fit%info == 0 .and. size(fit%difference_polynomial) == 5 .and. &
      fit%seasonal_difference_order == 1, 'seasonal differencing fit')

contains

   subroutine check(condition, label)
      !! Stop the test when an assertion fails.
      logical, intent(in) :: condition !! Flag controlling condition.
      character(*), intent(in) :: label !! Label.

      if (.not. condition) then
         write (*, '(a)') 'FAILED: '//label
         error stop 1
      end if
   end subroutine check

end program test_robustarima
