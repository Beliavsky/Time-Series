! SPDX-License-Identifier: GPL-3.0-only
! SPDX-FileComment: Regression tests for the R EXPARMA package translation.
program test_exparma
   use kind_mod, only: dp
   use exparma_mod
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   type(exparma_fit_t) :: evaluated, fitted
   type(exparma_selection_t) :: selection
   real(dp), parameter :: parameters(8) = [0.45_dp, -0.10_dp, &
      0.20_dp, 0.05_dp, 0.30_dp, 0.10_dp, 0.70_dp, 0.50_dp]
   real(dp) :: series(180), innovations(180), ar_weight, ma_weight
   real(dp), allocatable :: initial(:)
   integer :: time

   series = 0.0_dp
   innovations = 0.0_dp
   series(1:2) = [0.35_dp, -0.15_dp]
   do time = 3, size(series)
      innovations(time) = 0.035_dp*sin(0.79_dp*real(time, dp)) + &
         0.018_dp*cos(0.27_dp*real(time, dp))
      ar_weight = exp(-parameters(7)*series(time - 1)**2)
      ma_weight = exp(-parameters(8)*innovations(time - 1)**2)
      series(time) = (parameters(1) + parameters(3)*ar_weight)* &
         series(time - 1) + (parameters(2) + parameters(4)*ar_weight)* &
         series(time - 2) + (parameters(5) + parameters(6)*ma_weight)* &
         innovations(time - 1) + innovations(time)
   end do

   evaluated = exparma_evaluate(series, 2, 1, parameters)
   call check(evaluated%info == 0 .and. evaluated%ar_order == 2 .and. &
      evaluated%ma_order == 1 .and. &
      maxval(abs(evaluated%residuals(3:) - innovations(3:))) < 1.0e-13_dp .and. &
      abs(evaluated%rss - sum(innovations**2)) < 1.0e-13_dp, &
      'EXPARMA recursive residual evaluation')

   initial = exparma_initial_parameters(series, 2, 1)
   call check(size(initial) == 8 .and. all(ieee_is_finite(initial)) .and. &
      abs(initial(3) - 0.5_dp) < 1.0e-14_dp .and. &
      abs(initial(6) - 0.5_dp) < 1.0e-14_dp .and. &
      all(abs(initial(7:8) - 0.5_dp) < 1.0e-14_dp), &
      'Hannan-Rissanen EXPARMA initial values')

   fitted = exparma_fit(series, 2, 1, parameters, 250, 1.0e-7_dp)
   call check(fitted%info == 0 .and. fitted%rss <= evaluated%rss + 1.0e-10_dp .and. &
      all(ieee_is_finite(fitted%phi)) .and. all(ieee_is_finite(fitted%pi)) .and. &
      all(ieee_is_finite(fitted%theta)) .and. &
      all(ieee_is_finite(fitted%delta)), 'EXPARMA BFGS estimation')

   selection = exparma_select(series, 2, 2, 180, 1.0e-6_dp)
   call check(selection%info == 0 .and. selection%selected_ar_order >= 1 .and. &
      selection%selected_ar_order <= 2 .and. &
      selection%selected_ma_order >= 1 .and. &
      selection%selected_ma_order <= 2 .and. &
      all(ieee_is_finite(selection%aic)), 'EXPARMA AIC order selection')

   print '(a)', 'EXPARMA tests passed'

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

end program test_exparma
