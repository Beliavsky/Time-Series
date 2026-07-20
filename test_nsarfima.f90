! SPDX-License-Identifier: GPL-3.0-or-later
! SPDX-FileComment: Regression tests for the nsarfima translation.
program test_nsarfima
   use kind_mod, only: dp
   use nsarfima_mod
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   real(dp), parameter :: series(10) = [0.4_dp, -0.2_dp, 0.7_dp, 1.1_dp, &
      -0.5_dp, 0.3_dp, 0.9_dp, -0.8_dp, 0.2_dp, 0.6_dp]
   real(dp), parameter :: innovations(8) = [1.0_dp, -0.5_dp, 0.7_dp, 1.2_dp, &
      -0.3_dp, 0.4_dp, -0.8_dp, 0.9_dp]
   real(dp), parameter :: residual_reference(10) = [0.13_dp, -0.535_dp, &
      0.6455_dp, 0.67549_dp, -1.219225_dp, 0.27490386_dp, 0.65827519_dp, &
      -1.372921923_dp, 0.37527715878_dp, 0.47874722213_dp]
   type(nsarfima_filter_t) :: filtered
   type(nsarfima_simulation_t) :: simulation, stationary_simulation, first, second
   type(nsarfima_simulation_t) :: generated_simulation
   type(nsarfima_fit_t) :: mde, pml, fixed_pml
   real(dp), allocatable :: residuals(:), convolved(:), generated(:), noise(:)
   integer :: i

   convolved = nsarfima_convolve([1.0_dp, 2.0_dp, 3.0_dp], &
      [4.0_dp, 5.0_dp, 6.0_dp])
   call check(maxval(abs(convolved - [4.0_dp, 13.0_dp, 28.0_dp])) < 1.0e-13_dp, &
      'causal FFT convolution')

   residuals = nsarfima_residuals(series, 0.6_dp, [0.2_dp], [-0.3_dp])
   call check(maxval(abs(residuals - residual_reference)) < 2.0e-9_dp, &
      'nonstationary ARFIMA residuals against R')
   filtered = nsarfima_residual_acf(series, 0.6_dp, [0.2_dp], [-0.3_dp], 3)
   call check(filtered%info == 0 .and. maxval(abs(filtered%autocorrelation - &
      [1.0_dp, -0.4037757592285856_dp, -0.4544593952661507_dp, &
      0.6358704663999711_dp])) < 2.0e-12_dp, 'Mayoral residual ACF against R')

   simulation = nsarfima_simulate_from_innovations(4, 0.6_dp, [0.2_dp], &
      [-0.3_dp], innovations, mean=0.5_dp, burn_in=2)
   stationary_simulation = nsarfima_simulate_from_innovations(4, 0.6_dp, &
      [0.2_dp], [-0.3_dp], innovations, mean=0.5_dp, burn_in=2, &
      stationary_integration=.true.)
   call check(maxval(abs(simulation%series - [1.224_dp, 1.612528_dp, &
      0.555712_dp, 1.7691456_dp])) < 2.0e-12_dp, &
      'default nsarfima simulation against R')
   call check(maxval(abs(stationary_simulation%series - [3.224_dp, 4.112528_dp, &
      3.555712_dp, 5.2691456_dp])) < 2.0e-12_dp, &
      'stationary integer-integration simulation against R')
   first = nsarfima_simulate(12, 0.3_dp, [0.2_dp], [-0.1_dp], burn_in=4, seed=92)
   second = nsarfima_simulate(12, 0.3_dp, [0.2_dp], [-0.1_dp], burn_in=4, seed=92)
   call check(first%info == 0 .and. maxval(abs(first%series - second%series)) < &
      1.0e-14_dp, 'seeded Gaussian simulation')

   allocate(noise(80), generated(80))
   do i = 1, 80
      noise(i) = sin(1.31_dp*real(i, dp)) + 0.25_dp*cos(0.23_dp*real(i, dp))
   end do
   generated_simulation = nsarfima_simulate_from_innovations(80, 0.32_dp, &
      [0.25_dp], [-0.15_dp], noise, burn_in=0)
   generated = generated_simulation%series
   mde = nsarfima_mde(generated, [0.1_dp], [0.05_dp], 0.2_dp, &
      [0.0_dp, 0.8_dp], lag_max=8, estimate_mean=.false., max_iterations=150)
   pml = nsarfima_pml(generated, [0.1_dp], [0.05_dp], 0.2_dp, &
      [0.0_dp, 0.8_dp], estimate_mean=.false., max_iterations=150, &
      information_resolution=1024)
   fixed_pml = nsarfima_pml(generated, [0.1_dp], [0.05_dp], 0.0_dp, &
      [0.0_dp], estimate_mean=.false., max_iterations=100, &
      information_resolution=512)
   call check(mde%info == 0 .and. all(ieee_is_finite(mde%standard_error)) .and. &
      all(shape(mde%covariance) == [3, 3]) .and. size(mde%residuals) == 80 .and. &
      mde%d >= 0.0_dp .and. mde%d <= 0.8_dp, 'Mayoral minimum-distance fit')
   call check(pml%info == 0 .and. pml%innovation_variance > 0.0_dp .and. &
      all(ieee_is_finite(pml%standard_error)) .and. &
      all(shape(pml%covariance) == [4, 4]) .and. size(pml%residuals) == 80 .and. &
      pml%d >= 0.0_dp .and. pml%d <= 0.8_dp, 'Beran pseudo-likelihood fit')
   call check(fixed_pml%info == 0 .and. all(shape(fixed_pml%covariance) == [3, 3]) .and. &
      all(ieee_is_finite(fixed_pml%standard_error)), &
      'fixed fractional order pseudo-likelihood inference')

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

end program test_nsarfima
