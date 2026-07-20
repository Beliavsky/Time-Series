! SPDX-License-Identifier: MIT
! SPDX-FileComment: Numerical tests for reusable finite-mixture algorithms.
program test_mixture
   use kind_mod, only: dp
   use mixture_mod, only: gaussian_mixture_t, fit_gaussian_mixture
   implicit none

   real(dp) :: values(12)
   type(gaussian_mixture_t) :: fit

   values = [-2.2_dp, -2.1_dp, -2.0_dp, -1.9_dp, -1.8_dp, -1.7_dp, &
      1.7_dp, 1.8_dp, 1.9_dp, 2.0_dp, 2.1_dp, 2.2_dp]
   fit = fit_gaussian_mixture(values, 2)
   call assert_true(fit%info == 0, 'Gaussian mixture fit status')
   call assert_true(abs(sum(fit%weight) - 1.0_dp) < 1.0e-12_dp, &
      'Gaussian mixture weights')
   call assert_true(fit%mean(1) < -1.5_dp .and. fit%mean(2) > 1.5_dp, &
      'Gaussian mixture separation')
   call assert_true(all(fit%variance > 0.0_dp), &
      'Gaussian mixture positive variances')

   print '(a)', 'mixture tests passed'

contains

   subroutine assert_true(condition, message)
      !! Stop the test program when a logical assertion fails.
      logical, intent(in) :: condition !! Assertion condition.
      character(len=*), intent(in) :: message !! Failure description.

      if (.not. condition) error stop 'FAIL: '//trim(message)
   end subroutine assert_true

end program test_mixture
