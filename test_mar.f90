! SPDX-License-Identifier: GPL-2.0-or-later
! SPDX-FileComment: Regression tests for the mAr translation.
program test_mar
   use kind_mod, only: dp
   use mar_mod
   use time_series_random_mod, only: set_random_seed
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none

   type(mar_fit_t) :: fit
   type(mar_modes_t) :: modes
   type(mar_pca_fit_t) :: pca_fit
   type(mar_simulation_t) :: simulation, random_simulation
   real(dp) :: x(12, 2), expected_intercept(2), expected_ar(2, 2)
   real(dp) :: expected_covariance(2, 2), coefficients(2, 2)
   real(dp) :: intercept(2), ar(2, 2, 1), innovations(5, 2)
   real(dp) :: covariance(2, 2), pca_data(12, 3)
   complex(dp) :: expected_root
   integer :: row

   x = reshape([ &
      1.0_dp, 1.3_dp, 1.1_dp, 0.8_dp, 0.6_dp, 0.5_dp, &
      0.2_dp, 0.1_dp, -0.1_dp, -0.2_dp, -0.15_dp, -0.3_dp, &
      2.0_dp, 1.7_dp, 1.4_dp, 1.2_dp, 1.0_dp, 0.7_dp, &
      0.6_dp, 0.3_dp, 0.2_dp, 0.0_dp, -0.1_dp, -0.05_dp], [12, 2])
   expected_intercept = [-0.28462474909104962_dp, -0.06610626295833047_dp]
   expected_ar = reshape([ &
      0.0084954404048515064_dp, 0.18620209898798246_dp, &
      0.77079119132406504_dp, 0.74646978697261479_dp], [2, 2])
   expected_covariance = reshape([ &
      0.0058567664764433407_dp, 0.0020020065425835447_dp, &
      0.0020020065425835447_dp, 0.0067896047444019041_dp], [2, 2])
   fit = mar_estimate(x, 1)
   call check(fit%info == 0, 'stabilized estimator status')
   call check(maxval(abs(fit%intercept - expected_intercept)) < 2.0e-11_dp, &
      'stabilized estimator intercept')
   call check(maxval(abs(fit%ar(:, :, 1) - expected_ar)) < 2.0e-11_dp, &
      'stabilized estimator coefficients')
   call check(maxval(abs(fit%covariance - expected_covariance)) < 2.0e-11_dp, &
      'stabilized estimator covariance')
   call check(abs(fit%sbc + 4.7838631570196117_dp) < 2.0e-11_dp, &
      'stabilized estimator SBC')
   call check(maxval(abs(fit%residuals(1, :) - &
      [0.03454692603806797_dp, 0.08696459002511836_dp])) < 2.0e-11_dp, &
      'stabilized estimator residuals')

   coefficients = 0.0_dp
   coefficients(1, 1) = 0.8_dp
   coefficients(2, 2) = -0.5_dp
   modes = mar_eigenmodes(coefficients)
   call check(modes%info == 0, 'real eigenmode status')
   call check(modes%stable, 'stable real eigenmodes')
   call check(any(abs(modes%eigenvalues - cmplx(0.8_dp, 0.0_dp, dp)) < &
      1.0e-10_dp), 'positive real root')
   call check(any(abs(modes%eigenvalues - cmplx(-0.5_dp, 0.0_dp, dp)) < &
      1.0e-10_dp), 'negative real root')
   call check(abs(modes%damping_times(1) + 1.0_dp/log(0.8_dp)) < &
      1.0e-9_dp, 'sorted real damping time')

   coefficients = reshape([0.0_dp, 1.0_dp, -0.81_dp, 0.0_dp], [2, 2])
   modes = mar_eigenmodes(coefficients)
   call check(modes%info == 0, 'oscillatory eigenmode status')
   expected_root = cmplx(0.0_dp, 0.9_dp, dp)
   call check(any(abs(modes%eigenvalues - expected_root) < 1.0e-8_dp), &
      'positive imaginary root')
   call check(any(abs(modes%eigenvalues - conjg(expected_root)) < 1.0e-8_dp), &
      'negative imaginary root')
   call check(maxval(abs(modes%periods - 4.0_dp)) < 1.0e-8_dp, &
      'oscillatory periods')

   intercept = [0.2_dp, 0.4_dp]
   ar = 0.0_dp
   ar(1, 1, 1) = 0.5_dp
   ar(2, 2, 1) = 0.25_dp
   innovations = 0.0_dp
   simulation = mar_simulate_from_innovations(intercept, ar, innovations, 2)
   call check(simulation%info == 0, 'deterministic simulation status')
   call check(size(simulation%series, 1) == 3, 'simulation burn-in size')
   call check(maxval(abs(simulation%series(:, 1) - 0.4_dp)) < 1.0e-13_dp, &
      'first stationary mean')
   call check(maxval(abs(simulation%series(:, 2) - 0.4_dp/0.75_dp)) < &
      1.0e-13_dp, 'second stationary mean')

   do row = 1, 12
      pca_data(row, 1) = sin(0.4_dp*real(row, dp)) + 0.03_dp*real(row, dp)
      pca_data(row, 2) = cos(0.3_dp*real(row, dp)) - 0.02_dp*real(row, dp)
      pca_data(row, 3) = 0.4_dp*pca_data(row, 1) - &
         0.7_dp*pca_data(row, 2) + 0.1_dp*sin(real(row, dp))
   end do
   pca_fit = mar_pca(pca_data, 1, 2)
   call check(pca_fit%info == 0, 'PCA mAr status')
   call check(all(pca_fit%fraction_variance > 0.0_dp) .and. &
      all(pca_fit%fraction_variance <= 1.0_dp), 'PCA explained variance range')
   call check(pca_fit%fraction_variance(2) >= pca_fit%fraction_variance(1), &
      'PCA cumulative variance order')
   call check(all(ieee_is_finite(real(pca_fit%eigenvectors, dp))), &
      'mapped PCA eigenvectors')

   covariance = reshape([0.5_dp, 0.1_dp, 0.1_dp, 0.4_dp], [2, 2])
   call set_random_seed(731)
   random_simulation = mar_simulate(intercept, ar, covariance, 8, 10)
   call check(random_simulation%info == 0, 'Gaussian simulation status')
   call check(all(shape(random_simulation%series) == [8, 2]), &
      'Gaussian simulation shape')
   call check(all(ieee_is_finite(random_simulation%series)), &
      'Gaussian simulation finite')

   print '(a)', 'mAr tests passed'

contains

   subroutine check(condition, message)
      ! Stop the test program when a condition fails.
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message

      if (.not. condition) then
         print '(a)', 'FAILED: '//trim(message)
         error stop 1
      end if
   end subroutine check

end program test_mar
