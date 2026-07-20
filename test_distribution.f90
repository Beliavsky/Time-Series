! SPDX-License-Identifier: MIT
! SPDX-FileComment: Tests for shared standardized probability distributions.
program test_distribution
   use kind_mod, only: dp
   use distribution_mod, only: distribution_normal, distribution_student, &
      distribution_ged, distribution_skew_normal, distribution_skew_student, &
      distribution_skew_ged, distribution_johnson_su, distribution_nig, &
      distribution_ghyp, distribution_gh_skew_student, &
      standardized_log_density, standardized_cdf, standardized_quantile, &
      standardized_skewness, standardized_excess_kurtosis, &
      distribution_cdf, distribution_quantile, random_standardized, &
      distribution_name
   use distribution_fit_mod, only: distribution_fit_t, fit_distribution
   use random_mod, only: set_random_seed
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   integer, parameter :: distribution_count = 10
   integer, parameter :: draws = 2000
   integer :: codes(distribution_count), distribution, draw
   real(dp) :: shape, skew, lambda, value, total, squares
   real(dp) :: sample(300), probability, quantile
   type(distribution_fit_t) :: fitted

   codes = [distribution_normal, distribution_student, distribution_ged, &
      distribution_skew_normal, distribution_skew_student, &
      distribution_skew_ged, distribution_johnson_su, distribution_nig, &
      distribution_ghyp, distribution_gh_skew_student]
   call set_random_seed(271828)
   do distribution = 1, distribution_count
      call test_parameters(codes(distribution), shape, skew, lambda)
      if (.not. ieee_is_finite(standardized_log_density(0.2_dp, &
         codes(distribution), shape, skew, lambda))) &
         error stop 'nonfinite standardized density'
      if (len_trim(distribution_name(codes(distribution))) == 0) &
         error stop 'missing distribution name'
      total = 0.0_dp
      squares = 0.0_dp
      do draw = 1, draws
         value = random_standardized(codes(distribution), shape, skew, lambda)
         total = total + value
         squares = squares + value**2
      end do
      if (abs(total/real(draws, dp)) > 0.15_dp) &
         error stop 'standardized distribution is not centered'
      if (abs(squares/real(draws, dp) - 1.0_dp) > 0.35_dp) &
         error stop 'standardized distribution does not have unit variance'
   end do

   probability = standardized_cdf(0.0_dp, distribution_normal, 0.0_dp, &
      1.0_dp, -0.5_dp)
   if (abs(probability - 0.5_dp) > 2.0e-4_dp) &
      error stop 'standardized normal CDF failed'
   if (abs(distribution_cdf(0.7_dp, distribution_normal, 0.7_dp, 1.8_dp, &
      0.0_dp, 1.0_dp, -0.5_dp) - 0.5_dp) > 1.0e-12_dp) &
      error stop 'location-scale CDF failed'
   if (abs(distribution_quantile(0.5_dp, distribution_normal, 0.7_dp, &
      1.8_dp, 0.0_dp, 1.0_dp, -0.5_dp) - 0.7_dp) > 1.0e-10_dp) &
      error stop 'location-scale quantile failed'
   quantile = standardized_quantile(0.95_dp, distribution_student, 8.0_dp, &
      1.0_dp, -0.5_dp)
   if (abs(standardized_cdf(quantile, distribution_student, 8.0_dp, 1.0_dp, &
      -0.5_dp) - 0.95_dp) > 2.0e-4_dp) &
      error stop 'standardized quantile inversion failed'
   if (abs(standardized_skewness(distribution_normal, 0.0_dp, 1.0_dp, &
      -0.5_dp)) > 1.0e-8_dp) error stop 'normal skewness failed'
   if (abs(standardized_excess_kurtosis(distribution_normal, 0.0_dp, 1.0_dp, &
      -0.5_dp)) > 2.0e-3_dp) error stop 'normal kurtosis failed'
   call set_random_seed(314159)
   do draw = 1, size(sample)
      sample(draw) = 0.7_dp + 1.8_dp*random_standardized( &
         distribution_normal, 0.0_dp, 1.0_dp, -0.5_dp)
   end do
   fitted = fit_distribution(sample, distribution_normal, max_iterations=100)
   if (fitted%info /= 0 .or. abs(fitted%location - 0.7_dp) > 0.25_dp .or. &
      abs(fitted%scale - 1.8_dp) > 0.25_dp) &
      error stop 'distribution fitting failed'

   print '(a)', 'distribution tests passed'

contains

   pure subroutine test_parameters(distribution, shape, skew, lambda)
      !! Return stable test parameters for one distribution.
      integer, intent(in) :: distribution !! Distribution code.
      real(dp), intent(out) :: shape !! Shape parameter.
      real(dp), intent(out) :: skew !! Skew parameter.
      real(dp), intent(out) :: lambda !! GH lambda parameter.

      shape = 1.5_dp
      skew = 1.25_dp
      lambda = -0.5_dp
      if (distribution == distribution_student .or. &
         distribution == distribution_skew_student) shape = 8.0_dp
      if (distribution == distribution_johnson_su) then
         shape = 1.5_dp
         skew = -0.25_dp
      else if (distribution == distribution_nig) then
         shape = 3.0_dp
         skew = -0.20_dp
      else if (distribution == distribution_ghyp) then
         shape = 3.0_dp
         skew = 0.20_dp
         lambda = 1.0_dp
      else if (distribution == distribution_gh_skew_student) then
         shape = 8.0_dp
         skew = 0.30_dp
      end if
   end subroutine test_parameters

end program test_distribution
