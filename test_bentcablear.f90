! SPDX-License-Identifier: GPL-3.0-or-later
! SPDX-FileComment: Regression tests for the R bentcableAR translation.
program test_bentcablear
   use kind_mod, only: dp
   use bentcablear_mod
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   type(bentcable_residuals_t) :: residual
   type(bentcable_profile_t) :: profile
   type(bentcable_fit_t) :: fit
   type(bentcable_change_t) :: change
   type(bentcable_ar_covariance_t) :: dependence
   type(bentcable_fit_t) :: iterative_fit
   real(dp), parameter :: beta(3) = [2.0_dp, 0.15_dp, -0.35_dp]
   real(dp), parameter :: tau = 40.0_dp, gamma = 5.0_dp
   real(dp) :: data(100), time(100), errors(100), innovations(100)
   real(dp) :: basis_values(5), expected_change
   real(dp), allocatable :: fisher(:, :)
   integer :: i

   basis_values = bentcable_basis([-2.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, &
      2.0_dp], 1.0_dp)
   call check(maxval(abs(basis_values - &
      [0.0_dp, 0.0_dp, 0.25_dp, 1.0_dp, 2.0_dp])) < 1.0e-14_dp, &
      'linear-quadratic-linear cable basis')
   call check(abs(bentcable_basis(2.0_dp, 0.0_dp) - 2.0_dp) < 1.0e-14_dp .and. &
      bentcable_basis(-2.0_dp, 0.0_dp) == 0.0_dp, 'broken-stick basis')

   time = [(real(i - 1, dp), i=1, size(time))]
   errors = 0.0_dp
   innovations = 0.0_dp
   do i = 1, size(data)
      innovations(i) = 0.06_dp*sin(0.71_dp*real(i, dp)) + &
         0.025_dp*cos(0.29_dp*real(i, dp))
      errors(i) = innovations(i)
      if (i > 1) errors(i) = errors(i) + 0.4_dp*errors(i - 1)
      data(i) = bentcable_value(time(i), beta(1), beta(2), beta(3), &
         tau, gamma) + errors(i)
   end do
   residual = bentcable_residuals(data, time, beta, tau, gamma, [0.4_dp])
   call check(residual%info == 0 .and. &
      maxval(abs(residual%innovations - innovations(2:))) < 1.0e-13_dp .and. &
      abs(residual%rss - sum(innovations(2:)**2)) < 1.0e-13_dp, &
      'conditional AR cable residuals')
   call check(bentcable_stationary([0.4_dp]) .and. &
      .not. bentcable_stationary([1.1_dp]), 'AR stationarity check')
   dependence = bentcable_ar_covariance(residual%residuals, 1)
   call check(dependence%info == 0 .and. &
      all(shape(dependence%covariance) == [100, 100]) .and. &
      maxval(abs(dependence%covariance - transpose(dependence%covariance))) < &
      1.0e-12_dp .and. abs(dependence%covariance(1, 1) - 1.0_dp) < &
      1.0e-14_dp, 'Yule-Walker AR correlation matrix')

   profile = bentcable_profile(data, [35.0_dp, 40.0_dp, 45.0_dp], &
      [3.0_dp, 5.0_dp, 7.0_dp], 0, time)
   call check(profile%info == 0 .and. all(ieee_is_finite(profile%deviance)) .and. &
      abs(maxval(profile%deviance)) < 1.0e-12_dp .and. &
      size(profile%initial_cable) == 5, 'bent-cable profile deviance surface')

   fit = bentcable_fit(data, 1, [beta, tau, gamma], [0.4_dp], time, &
      max_iterations=500, tolerance=1.0e-7_dp)
   call check(fit%info == 0 .and. fit%rss <= residual%rss + 1.0e-10_dp .and. &
      fit%stationary .and. all(ieee_is_finite(fit%beta)) .and. &
      all(ieee_is_finite(fit%innovations)), 'joint bent-cable AR CSS fit')
   iterative_fit = bentcable_fit_iterative_yw(data, 1, [beta, tau, gamma], &
      [0.4_dp], time, max_iterations=20, tolerance=1.0e-4_dp)
   call check(iterative_fit%info == 0 .and. iterative_fit%stationary .and. &
      all(ieee_is_finite(iterative_fit%beta)) .and. &
      ieee_is_finite(iterative_fit%rss), 'iterative Yule-Walker cable fit')

   fisher = bentcable_fisher_information(fit)
   change = bentcable_change_interval(fit, 0.95_dp)
   expected_change = fit%tau - fit%gamma - &
      2.0_dp*fit%beta(2)*fit%gamma/fit%beta(3)
   call check(all(shape(fisher) == [5, 5]) .and. &
      maxval(abs(fisher - transpose(fisher))) < 1.0e-10_dp .and. &
      change%info == 0 .and. abs(change%estimate - expected_change) < 1.0e-12_dp .and. &
      change%lower < change%estimate .and. change%upper > change%estimate, &
      'Fisher information and critical-transition-point interval')

   fit = bentcable_fit(data, 0, [beta, tau], time=time, stick=.true., &
      max_iterations=400)
   call check(fit%info == 0 .and. fit%stick .and. fit%gamma == 0.0_dp .and. &
      all(ieee_is_finite(fit%fitted)), 'broken-stick regression fit')

   print '(a)', 'bentcableAR tests passed'

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

end program test_bentcablear
