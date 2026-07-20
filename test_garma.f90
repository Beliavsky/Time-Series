! SPDX-License-Identifier: GPL-3.0-only
! SPDX-FileComment: Regression tests for the garma translation.
program test_garma
   use kind_mod, only: dp
   use garma_mod
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   real(dp), parameter :: series(16) = [0.4_dp, -0.2_dp, 0.7_dp, 1.1_dp, &
      -0.5_dp, 0.3_dp, 0.9_dp, -0.8_dp, 0.2_dp, 0.6_dp, 0.1_dp, -0.4_dp, &
      0.8_dp, 0.5_dp, -0.1_dp, 0.2_dp]
   real(dp), parameter :: coefficient_reference(8) = [1.0_dp, 0.2_dp, &
      -0.08_dp, -0.152_dp, -0.0736_dp, 0.041536_dp, 0.0899712_dp, &
      0.04764672_dp]
   real(dp), parameter :: spectrum_reference(8) = [0.028480242930789608_dp, &
      0.072560167798748423_dp, 0.075850283571967020_dp, &
      0.103049999999999975_dp, 1.573563517740589335_dp, &
      0.149139832201251543_dp, 0.089905955756653999_dp, &
      0.078399999999999984_dp]
   type(garma_periodogram_t) :: periodogram
   type(garma_factors_t) :: gsp, lpr, peaks
   type(garma_fit_t) :: fit, css_fit, wll_fit, estimated_fit
   type(garma_forecast_t) :: prediction
   type(garma_gof_t) :: diagnostic
   type(garma_regression_fit_t) :: regression_fit, integrated_fit
   type(garma_accuracy_t) :: accuracy
   real(dp), allocatable :: coefficients(:), residuals(:), generated(:)
   real(dp), allocatable :: regressors(:, :), future_regressors(:, :), regression_series(:)
   real(dp) :: whittle, full_whittle, wll
   integer :: i

   coefficients = garma_gegenbauer_coefficients(8, 0.2_dp, 0.5_dp)
   call check(maxval(abs(coefficients - coefficient_reference)) < 2.0e-14_dp, &
      'Gegenbauer coefficients against R')
   periodogram = garma_periodogram(series)
   call check(periodogram%info == 0 .and. maxval(abs(periodogram%spectrum - &
      spectrum_reference)) < 3.0e-14_dp, 'raw detrended periodogram against R')

   peaks = garma_semiparametric(series, factor_count=2, peak_bandwidth=1)
   call check(peaks%info == 0 .and. maxval(abs(peaks%frequency - &
      [0.3125_dp, 0.4375_dp])) < 1.0e-14_dp, 'Yajima peak selection against R')
   gsp = garma_semiparametric(series, periods=[4.0_dp], alpha=0.8_dp, &
      method=garma_semiparametric_gsp)
   lpr = garma_semiparametric(series, periods=[4.0_dp], alpha=0.8_dp, &
      method=garma_semiparametric_lpr)
   call check(abs(gsp%d(1) - 0.9005713940505506_dp) < 3.0e-5_dp, &
      'Gaussian semiparametric estimate against R')
   call check(abs(lpr%d(1) - 0.5618209552150502_dp) < 2.0e-12_dp, &
      'log-periodogram regression estimate against R')

   whittle = garma_whittle_objective(series, [-0.5_dp], [0.2_dp], &
      [0.3_dp], [-0.1_dp])
   full_whittle = garma_whittle_objective(series, [-0.5_dp], [0.2_dp], &
      [0.3_dp], [-0.1_dp], include_log_term=.true.)
   wll = garma_wll_objective(series, [-0.5_dp], [0.2_dp], [0.3_dp], &
      [-0.1_dp], 0.7_dp)
   call check(abs(whittle - 0.04029314806077091_dp) < 2.0e-14_dp .and. &
      abs(full_whittle + 0.7794167515029308_dp) < 2.0e-13_dp .and. &
      abs(wll - 9.984005658158706_dp) < 2.0e-12_dp, &
      'Whittle and WLL objectives against R')

   residuals = garma_css_residuals(series, [real(dp) ::], [real(dp) ::], &
      [real(dp) ::], [real(dp) ::])
   call check(maxval(abs(residuals - series)) < 1.0e-14_dp, &
      'identity CSS filter')
   residuals = garma_css_residuals(series, [real(dp) ::], [real(dp) ::], &
      [0.3_dp], [real(dp) ::])
   call check(abs(residuals(2) + 0.32_dp) < 1.0e-14_dp .and. &
      abs(residuals(3) - 0.76_dp) < 1.0e-14_dp, 'short-memory CSS filter')

   allocate(generated(96))
   do i = 1, 96
      generated(i) = sin(2.0_dp*acos(-1.0_dp)*real(i, dp)/12.0_dp) + &
         0.35_dp*sin(1.17_dp*real(i, dp))
   end do
   fit = garma_fit(generated, [cos(2.0_dp*acos(-1.0_dp)/12.0_dp)], &
      [0.2_dp], [real(dp) ::], [real(dp) ::], method=garma_method_whittle, &
      max_iterations=1000, tolerance=1.0e-7_dp)
   call check(fit%info == 0 .and. fit%d(1) >= 0.0_dp .and. fit%d(1) <= 0.5_dp .and. &
      fit%innovation_variance > 0.0_dp .and. size(fit%residuals) == 96 .and. &
      all(ieee_is_finite(fit%standard_error)), 'fixed-frequency Whittle GARMA fit')
   css_fit = garma_fit(generated, [cos(2.0_dp*acos(-1.0_dp)/12.0_dp)], &
      [0.2_dp], [real(dp) ::], [real(dp) ::], method=garma_method_css, &
      max_iterations=500)
   wll_fit = garma_fit(generated, [cos(2.0_dp*acos(-1.0_dp)/12.0_dp)], &
      [0.2_dp], [real(dp) ::], [real(dp) ::], method=garma_method_wll, &
      max_iterations=700)
   estimated_fit = garma_fit(generated, [0.85_dp], [0.2_dp], [real(dp) ::], &
      [real(dp) ::], method=garma_method_whittle, estimate_frequencies=.true., &
      max_iterations=800)
   call check(css_fit%info == 0 .and. css_fit%innovation_variance > 0.0_dp .and. &
      all(shape(css_fit%covariance) == [1, 1]) .and. &
      css_fit%standard_error(1) > 0.0_dp .and. wll_fit%info == 0 .and. &
      wll_fit%innovation_variance > 0.0_dp .and. wll_fit%standard_error(1) > 0.0_dp, &
      'CSS and WLL GARMA fits')
   call check(estimated_fit%info == 0 .and. estimated_fit%u(1) > 0.0_dp .and. &
      estimated_fit%u(1) < 1.0_dp .and. estimated_fit%d(1) >= 0.0_dp .and. &
      estimated_fit%d(1) <= 0.5_dp, 'joint Gegenbauer pole estimation')
   prediction = garma_forecast(generated, fit%u, fit%d, fit%ar, fit%ma, 6)
   call check(prediction%info == 0 .and. size(prediction%mean) == 6 .and. &
      all(ieee_is_finite(prediction%mean)), 'Godet-style GARMA forecasts')
   diagnostic = garma_goodness_of_fit(fit%residuals)
   call check(diagnostic%info == 0 .and. size(diagnostic%p_value) == 47 .and. &
      all(diagnostic%p_value >= 0.0_dp) .and. all(diagnostic%p_value <= 1.0_dp), &
      'Bartlett periodogram goodness-of-fit diagnostic')

   allocate(regressors(96, 1), regression_series(96), future_regressors(4, 1))
   do i = 1, 96
      regressors(i, 1) = cos(0.37_dp*real(i, dp))
      regression_series(i) = 2.0_dp + 0.03_dp*real(i, dp) + &
         1.5_dp*regressors(i, 1) + 0.2_dp*generated(i)
   end do
   do i = 1, 4
      future_regressors(i, 1) = cos(0.37_dp*real(96 + i, dp))
   end do
   regression_fit = garma_regression_fit(regression_series, &
      [cos(2.0_dp*acos(-1.0_dp)/12.0_dp)], [0.2_dp], [real(dp) ::], &
      [real(dp) ::], regressors, include_mean=.true., include_drift=.true., &
      method=garma_method_whittle, max_iterations=700)
   prediction = garma_regression_forecast(regression_fit, future_regressors, 4)
   call check(regression_fit%info == 0 .and. &
      size(regression_fit%regression_coefficients) == 3 .and. &
      size(regression_fit%fitted) == 96 .and. size(regression_fit%residuals) == 96 .and. &
      prediction%info == 0 .and. all(ieee_is_finite(prediction%mean)), &
      'mean, drift, and external-regressor GARMA workflow')

   integrated_fit = garma_regression_fit(regression_series, [real(dp) ::], &
      [real(dp) ::], [0.1_dp], [real(dp) ::], difference_order=1, &
      method=garma_method_css, max_iterations=500)
   prediction = garma_regression_forecast(integrated_fit, horizon=4)
   call check(integrated_fit%info == 0 .and. &
      size(integrated_fit%differenced_series) == 95 .and. &
      size(integrated_fit%fitted) == 96 .and. prediction%info == 0 .and. &
      all(ieee_is_finite(prediction%mean)), 'integer differencing and reintegration')

   accuracy = garma_accuracy([2.0_dp, 4.0_dp, 5.0_dp], &
      [1.0_dp, 5.0_dp, 4.0_dp], [1.0_dp, 2.0_dp, 4.0_dp, 7.0_dp])
   call check(accuracy%info == 0 .and. abs(accuracy%mean_error - &
      1.0_dp/3.0_dp) < 1.0e-14_dp .and. abs(accuracy%root_mean_squared_error - &
      1.0_dp) < 1.0e-14_dp .and. abs(accuracy%mean_absolute_error - 1.0_dp) < &
      1.0e-14_dp .and. abs(accuracy%mean_absolute_scaled_error - 0.5_dp) < &
      1.0e-14_dp, 'forecast accuracy measures')

contains

   subroutine check(condition, label)
      !! Stop the test program when an assertion fails.
      logical, intent(in) :: condition !! Flag controlling condition.
      character(*), intent(in) :: label !! Label.

      if (.not. condition) then
         write (*, '(a)') 'FAILED: '//label
         error stop 1
      end if
   end subroutine check

end program test_garma
