! SPDX-License-Identifier: GPL-2.0-or-later
! SPDX-FileComment: Regression tests for the fracdiff translation.
program test_fracdiff
   use kind_mod, only: dp
   use fracdiff_mod
   use fourier_mod, only: fft_transform
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   type(fracdiff_semiparametric_t) :: gph, sperio
   type(fracdiff_filter_t) :: filtered
   type(fracdiff_fit_t) :: fd_fit, arma_fit
   type(fracdiff_simulation_t) :: deterministic_simulation
   type(fracdiff_simulation_t) :: random_first, random_second
   real(dp), parameter :: reference(10) = [0.4_dp, -0.2_dp, 0.7_dp, 1.1_dp, &
      -0.5_dp, 0.3_dp, 0.9_dp, -0.8_dp, 0.2_dp, 0.6_dp]
   real(dp), parameter :: difference_reference(10) = [0.13_dp, -0.496_dp, &
      0.5136_dp, 0.77536_dp, -0.952208_dp, 0.10943232_dp, 0.640658176_dp, &
      -1.1929074176_dp, 0.09413789184_dp, 0.398921535488_dp]
   real(dp), allocatable :: differenced(:), simulated(:), innovation(:)
   complex(dp), allocatable :: transformed(:), recovered(:)
   integer :: i

   transformed = fft_transform(cmplx([1.0_dp, 2.0_dp, -1.0_dp, 0.5_dp], &
      0.0_dp, dp))
   recovered = fft_transform(transformed, inverse=.true.)
   call check(maxval(abs(recovered - cmplx([1.0_dp, 2.0_dp, -1.0_dp, &
      0.5_dp], 0.0_dp, dp))) < 1.0e-13_dp, 'shared radix-2 FFT round trip')

   differenced = fracdiff_difference(reference, 0.2_dp)
   call check(size(differenced) == 10 .and. &
      maxval(abs(differenced - difference_reference)) < 2.0e-13_dp, &
      'Jensen-Nielsen FFT fractional differencing against R')
   differenced = fracdiff_difference(reference, 0.0_dp)
   call check(maxval(abs(differenced - (reference - 0.27_dp))) < 1.0e-13_dp, &
      'zero-order fractional differencing demeans the series')

   gph = fracdiff_gph(reference, 0.6_dp)
   sperio = fracdiff_sperio(reference, 0.6_dp, 0.8_dp)
   call check(gph%info == 0 .and. gph%bandwidth == 3 .and. &
      abs(gph%d + 0.779275418246826_dp) < 2.0e-13_dp .and. &
      abs(gph%asymptotic_standard_error - 0.9250491659969612_dp) < 2.0e-13_dp .and. &
      abs(gph%regression_standard_error - 1.818443780471036_dp) < 2.0e-13_dp, &
      'GPH semiparametric estimator against R')
   call check(sperio%info == 0 .and. sperio%bandwidth == 3 .and. &
      abs(sperio%d + 0.9548676734623545_dp) < 2.0e-13_dp .and. &
      abs(sperio%asymptotic_standard_error - 0.4102752926967372_dp) < &
      2.0e-13_dp .and. abs(sperio%regression_standard_error - &
      0.01200570055115985_dp) < 2.0e-13_dp, &
      'Sperio tapered estimator against R')

   filtered = fracdiff_hr_filter(reference, 0.0_dp, 5)
   call check(filtered%info == 0 .and. abs(filtered%mean - 0.27_dp) < &
      1.0e-14_dp .and. maxval(abs(filtered%filtered - reference + 0.27_dp)) < &
      1.0e-14_dp .and. maxval(abs(filtered%prediction_variance - 1.0_dp)) < &
      1.0e-14_dp .and. abs(filtered%log_variance) < 1.0e-14_dp, &
      'Haslett-Raftery filter at zero fractional order')
   filtered = fracdiff_hr_filter(reference, 0.2_dp, 5)
   call check(filtered%info == 0 .and. all(ieee_is_finite(filtered%filtered)) .and. &
      all(filtered%prediction_variance > 0.0_dp), &
      'truncated Haslett-Raftery fractional filter')

   deterministic_simulation = fracdiff_simulate_from_innovations(4, &
      [real(dp) ::], [real(dp) ::], 0.0_dp, [1.0_dp, -1.0_dp, 2.0_dp, -2.0_dp], &
      mean=3.0_dp, burn_in=0)
   call check(deterministic_simulation%info == 0 .and. &
      maxval(abs(deterministic_simulation%series - &
      [4.0_dp, 2.0_dp, 5.0_dp, 1.0_dp])) < 1.0e-14_dp, &
      'fracdiff simulation reuses Durbin-Levinson generation')
   deterministic_simulation = fracdiff_simulate_from_innovations(3, &
      [real(dp) ::], [0.2_dp], 0.0_dp, [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], &
      burn_in=0)
   call check(maxval(abs(deterministic_simulation%series - &
      [1.8_dp, 2.6_dp, 3.4_dp])) < 1.0e-14_dp, &
      'fracdiff opposite-sign moving-average simulation')
   deterministic_simulation = fracdiff_simulate_from_innovations(4, &
      [0.3_dp], [0.2_dp], 0.2_dp, [1.0_dp, -0.5_dp, 0.7_dp, 1.2_dp, -0.3_dp], &
      burn_in=0)
   call check(maxval(abs(deterministic_simulation%series - &
      [-0.4550398358646283_dp, 0.6805288390371383_dp, &
      1.473872696554201_dp, 0.2654470888018319_dp])) < 2.0e-14_dp, &
      'fracdiff simulation against R reference')
   random_first = fracdiff_simulate(8, [0.3_dp], [real(dp) ::], 0.2_dp, &
      burn_in=2, seed=419)
   random_second = fracdiff_simulate(8, [0.3_dp], [real(dp) ::], 0.2_dp, &
      burn_in=2, seed=419)
   call check(random_first%info == 0 .and. random_second%info == 0 .and. &
      maxval(abs(random_first%series - random_second%series)) < 1.0e-14_dp, &
      'seeded fracdiff Gaussian simulation')

   allocate(innovation(40), simulated(40))
   do i = 1, 40
      innovation(i) = sin(1.37_dp*real(i, dp)) + 0.2_dp*cos(0.31_dp*real(i, dp))
   end do
   simulated = innovation
   do i = 2, 40
      simulated(i) = 0.45_dp*simulated(i - 1) + innovation(i)
      if (i > 2) simulated(i) = simulated(i) - 0.2_dp*innovation(i - 1)
   end do
   fd_fit = fracdiff_fit(simulated, [real(dp) ::], [real(dp) ::], 0.1_dp, &
      truncation=20, max_iterations=100, tolerance=1.0e-5_dp)
   call check(fd_fit%info == 0 .and. fd_fit%converged .and. &
      fd_fit%model%long_memory_parameter >= 0.0_dp .and. &
      fd_fit%model%long_memory_parameter < 0.5_dp .and. &
      abs(fd_fit%log_likelihood + 44.10035379836742_dp) < 5.0e-4_dp .and. &
      abs(sqrt(fd_fit%innovation_variance) - 0.7287743305677344_dp) < &
      1.0e-5_dp .and. &
      all(shape(fd_fit%covariance) == [1, 1]), &
      'Haslett-Raftery fractional-noise fit')
   arma_fit = fracdiff_fit(simulated, [0.1_dp], [0.1_dp], 0.05_dp, &
      truncation=20, max_iterations=150, tolerance=1.0e-5_dp)
   call check(arma_fit%info == 0 .and. arma_fit%converged .and. &
      abs(arma_fit%model%ar(1) - 0.3094384468691782_dp) < 2.0e-3_dp .and. &
      abs(arma_fit%model%theta(1) + 0.8100580930696726_dp) < 2.0e-3_dp .and. &
      abs(arma_fit%log_likelihood + 32.165639642953998_dp) < 4.0e-2_dp .and. &
      all(shape(arma_fit%covariance) == [3, 3]) .and. &
      size(arma_fit%residuals) == 40 .and. size(arma_fit%fitted) == 40, &
      'joint Haslett-Raftery ARFIMA fit and inference')

contains

   subroutine check(condition, label)
      !! Stop the test program when one assertion fails.
      logical, intent(in) :: condition !! Flag controlling condition.
      character(*), intent(in) :: label !! Label.

      if (.not. condition) then
         write (*, '(a)') 'FAILED: '//label
         error stop 1
      end if
   end subroutine check

end program test_fracdiff
