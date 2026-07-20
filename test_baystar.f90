! SPDX-License-Identifier: GPL-2.0-or-later
! SPDX-FileComment: Regression tests for the R BAYSTAR translation.
program test_baystar
   use kind_mod, only: dp
   use baystar_mod
   use random_mod, only: set_random_seed
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   type(baystar_simulation_t) :: simulation
   type(baystar_prior_t) :: prior
   type(baystar_coefficient_posterior_t) :: coefficient_posterior
   type(baystar_variance_posterior_t) :: variance_posterior
   type(baystar_summary_t) :: summary
   type(baystar_fit_t) :: fit
   real(dp) :: innovations(700), probability_sum, likelihood
   real(dp), allocatable :: probability(:), draws(:, :)
   integer :: i

   do i = 1, size(innovations)
      innovations(i) = sin(0.73_dp*real(i, dp)) + &
         0.35_dp*cos(0.29_dp*real(i, dp))
   end do
   simulation = baystar_simulate_from_innovations(500, [0.25_dp, 0.40_dp], &
      [-0.25_dp, 0.35_dp], 0.04_dp, 0.06_dp, 1, 0.0_dp, [1], [1], &
      innovations, 200)
   call check(simulation%info == 0 .and. size(simulation%series) == 500 .and. &
      all(ieee_is_finite(simulation%series)), 'supplied-innovation TAR simulation')

   likelihood = baystar_log_likelihood(simulation%series, &
      [0.25_dp, 0.40_dp], [-0.25_dp, 0.35_dp], 0.04_dp, 0.06_dp, 1, &
      0.0_dp, [1], [1])
   call check(ieee_is_finite(likelihood), 'two-regime Gaussian likelihood')
   prior = baystar_default_prior(simulation%series, [1], [1])
   call check(prior%info == 0 .and. prior%variance_scale > 0.0_dp .and. &
      all(shape(prior%precision1) == [2, 2]), 'default BAYSTAR prior')
   coefficient_posterior = baystar_coefficient_posterior(1, simulation%series, &
      0.04_dp, 1, 0.0_dp, [1], [1], prior%mean1, prior%precision1)
   call check(coefficient_posterior%info == 0 .and. &
      all(shape(coefficient_posterior%covariance) == [2, 2]) .and. &
      all(ieee_is_finite(coefficient_posterior%mean)), &
      'conditional Gaussian coefficient posterior')
   variance_posterior = baystar_variance_posterior(2, simulation%series, &
      [-0.25_dp, 0.35_dp], 1, 0.0_dp, [1], [1], &
      prior%variance_degrees, prior%variance_scale)
   call check(variance_posterior%info == 0 .and. &
      variance_posterior%shape > 0.0_dp .and. variance_posterior%rate > 0.0_dp, &
      'conditional inverse-gamma variance posterior')
   probability = baystar_delay_probabilities(simulation%series, &
      [0.25_dp, 0.40_dp], [-0.25_dp, 0.35_dp], 0.04_dp, 0.06_dp, &
      0.0_dp, [1], [1], 2)
   probability_sum = sum(probability)
   call check(size(probability) == 2 .and. &
      abs(probability_sum - 1.0_dp) < 1.0e-13_dp .and. &
      all(probability >= 0.0_dp), 'delay-lag posterior probabilities')

   allocate(draws(5, 2))
   draws(:, 1) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
   draws(:, 2) = [2.0_dp, 4.0_dp, 6.0_dp, 8.0_dp, 10.0_dp]
   summary = baystar_summary(draws)
   call check(summary%info == 0 .and. &
      maxval(abs(summary%mean - [3.0_dp, 6.0_dp])) < 1.0e-14_dp .and. &
      maxval(abs(summary%median - [3.0_dp, 6.0_dp])) < 1.0e-14_dp, &
      'posterior summary statistics')

   call set_random_seed(7319)
   fit = baystar_fit(simulation%series, [1], [1], 180, 60, 0.03_dp, 2)
   call check(fit%info == 0 .and. fit%selected_delay >= 1 .and. &
      fit%selected_delay <= 2 .and. size(fit%residuals) > 0 .and. &
      ieee_is_finite(fit%dic) .and. fit%threshold_acceptance >= 0.0_dp .and. &
      fit%threshold_acceptance <= 1.0_dp, 'complete BAYSTAR MCMC fit')

   print '(a)', 'BAYSTAR tests passed'

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

end program test_baystar
