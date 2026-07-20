! SPDX-License-Identifier: GPL-3.0-or-later
! SPDX-FileComment: Regression tests for the arima2 translation.
program test_arima2
   use kind_mod, only: dp
   use arima2_mod
   use random_mod, only: set_random_seed
   implicit none
   type(arima2_roots_t) :: roots
   type(arima2_fit_t) :: restart_fit
   type(arima2_aic_table_t) :: aic_table
   type(arima2_profile_t) :: profile
   type(arima2_coefficient_samples_t) :: samples
   complex(dp) :: inverse_roots(2)
   real(dp), allocatable :: coefficients(:)
   real(dp) :: series(60), starts(3, 3), matching_error, regressor(60, 1)
   integer :: i

   inverse_roots = [cmplx(0.6_dp, 0.2_dp, dp), cmplx(0.6_dp, -0.2_dp, dp)]
   coefficients = inverse_roots_to_coefficients(inverse_roots)
   roots = arma_polynomial_roots(coefficients)
   matching_error = min(max(abs(1.0_dp/roots%roots(1) - inverse_roots(1)), &
      abs(1.0_dp/roots%roots(2) - inverse_roots(2))), &
      max(abs(1.0_dp/roots%roots(1) - inverse_roots(2)), &
      abs(1.0_dp/roots%roots(2) - inverse_roots(1))))
   call check(roots%info == 0 .and. matching_error < 1.0e-10_dp, &
      'inverse-root coefficient round trip')
   coefficients = durbin_levinson_coefficients([0.5_dp, -0.2_dp])
   call check(maxval(abs(coefficients - [0.6_dp, -0.2_dp])) < 1.0e-12_dp, &
      'Durbin-Levinson AR coefficients')
   coefficients = durbin_levinson_coefficients([0.5_dp, -0.2_dp], .true.)
   call check(maxval(abs(coefficients - [-0.6_dp, 0.2_dp])) < 1.0e-12_dp, &
      'Durbin-Levinson MA coefficients')

   series(1) = 0.2_dp
   do i = 2, 60
      series(i) = 0.65_dp*series(i - 1) + sin(1.3_dp*real(i, dp))
   end do
   starts(:, 1) = [0.0_dp, 0.0_dp, sum(series)/60.0_dp]
   starts(:, 2) = [0.4_dp, 0.1_dp, sum(series)/60.0_dp]
   starts(:, 3) = [-0.4_dp, -0.1_dp, sum(series)/60.0_dp]
   restart_fit = arima2_fit_from_starts(series, 1, 0, 1, 0, 0, 0, 1, starts, &
      max_repeats=3, max_iterations=80, tolerance=1.0e-5_dp)
   call check(restart_fit%info == 0 .and. restart_fit%starts_attempted >= 1 .and. &
      restart_fit%fit%likelihood%info == 0, 'multi-start exact ARIMA fit')
   call check(restart_fit%fit%likelihood%log_likelihood >= restart_fit%log_likelihoods(1) - 1.0e-8_dp, &
      'restart likelihood selection')
   aic_table = arima2_aic_table(series, 1, 0, 1, corrected=.true., &
      max_iterations=60, tolerance=1.0e-5_dp)
   call check(aic_table%info == 0 .and. all(aic_table%values < huge(1.0_dp)), &
      'ARIMA AICc table')
   profile = arima2_profile(series, restart_fit%fit, 1, &
      restart_fit%fit%coefficients(1) + [-0.1_dp, 0.0_dp, 0.1_dp], 60, 1.0e-5_dp)
   call check(profile%info == 0 .and. all(profile%log_likelihood > -huge(1.0_dp)) .and. &
      maxval(abs(profile%coefficients(1, :) - profile%parameter_values)) < 1.0e-12_dp, &
      'fixed-parameter likelihood profile')
   call set_random_seed(24680)
   samples = sample_arma_coefficients(2, 2, 1, 1, 8, 0.05_dp, .true., &
      intercept_mean=2.0_dp, max_attempts=5000)
   call check(samples%info == 0 .and. all(shape(samples%values) == [7, 8]), &
      'Durbin-Levinson coefficient sampling')
   do i = 1, 8
      roots = arma_polynomial_roots(samples%values(1:2, i))
      call check(roots%info == 0 .and. all(abs(roots%roots) > 1.0_dp), &
         'sampled AR stationarity')
      roots = arma_polynomial_roots(samples%values(3:4, i), .true.)
      call check(roots%info == 0 .and. all(abs(roots%roots) > 1.0_dp), &
         'sampled MA invertibility')
   end do
   coefficients = invert_ma_coefficients([2.0_dp])
   call check(size(coefficients) == 1 .and. abs(coefficients(1) - 0.5_dp) < 1.0e-12_dp, &
      'MA root inversion')
   call set_random_seed(86420)
   samples = sample_arma_coefficients(2, 2, 0, 0, 5, sampling_method= &
      arima2_sampling_uniform_roots, modulus_bounds=[0.1_dp, 0.9_dp])
   call check(samples%info == 0 .and. all(shape(samples%values) == [4, 5]), &
      'uniform inverse-root coefficient sampling')
   do i = 1, 5
      roots = arma_polynomial_roots(samples%values(1:2, i))
      call check(minval(1.0_dp/abs(roots%roots)) >= 0.1_dp .and. &
         maxval(1.0_dp/abs(roots%roots)) <= 0.9_dp, 'uniform AR inverse-root bounds')
   end do
   call set_random_seed(13579)
   restart_fit = arima2_fit(series, 1, 0, 1, 0, 0, 0, 1, max_starts=4, &
      max_repeats=3, max_iterations=60, tolerance=1.0e-5_dp)
   call check(restart_fit%info == 0 .and. restart_fit%starts_attempted >= 1 .and. &
      restart_fit%fit%uses_exact_likelihood, 'random-restart ARIMA2 fit')
   restart_fit = arima2_fit(series, 1, 0, 0, 0, 0, 0, 1, max_starts=1, &
      method=arima2_method_css, initial=[0.0_dp, 1.25_dp], estimated=[.true., .false.], &
      max_iterations=60, tolerance=1.0e-5_dp)
   call check(restart_fit%info == 0 .and. .not. restart_fit%fit%uses_exact_likelihood .and. &
      abs(restart_fit%fit%coefficients(2) - 1.25_dp) < 1.0e-12_dp, &
      'CSS method and fixed mean')
   do i = 1, 60
      regressor(i, 1) = cos(0.2_dp*real(i, dp))
   end do
   call set_random_seed(97531)
   restart_fit = arima2_fit(series, 1, 0, 1, 0, 0, 0, 1, max_starts=3, &
      method=arima2_method_ml, initial=[0.0_dp, 0.0_dp, 0.5_dp, 0.25_dp], &
      estimated=[.true., .true., .true., .false.], regressors=regressor, &
      max_repeats=3, max_iterations=60, tolerance=1.0e-5_dp)
   call check(restart_fit%info == 0 .and. restart_fit%fit%regression_count == 1 .and. &
      restart_fit%fit%uses_exact_likelihood .and. &
      abs(restart_fit%fit%coefficients(4) - 0.25_dp) < 1.0e-12_dp, &
      'ML regressors and fixed coefficient mask')
   print '(a)', 'All arima2_mod tests passed.'

contains

   subroutine check(ok, name)
      !! Stop the test program when a named assertion fails.
      logical, intent(in) :: ok !! Flag controlling ok.
      character(len=*), intent(in) :: name !! Name.
      if (.not. ok) then
         print '(a)', 'FAILED: '//name
         error stop 1
      end if
   end subroutine check
end program test_arima2
