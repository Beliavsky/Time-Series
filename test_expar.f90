! SPDX-License-Identifier: GPL-3.0-only
! SPDX-FileComment: Regression tests for the R EXPAR package translation.
program test_expar
   use kind_mod, only: dp
   use expar_mod
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   type(expar_fit_t) :: evaluated, fitted
   type(expar_selection_t) :: selection
   type(expar_forecast_t) :: forecast
   real(dp), parameter :: parameters(5) = [0.55_dp, -0.15_dp, &
      0.25_dp, 0.10_dp, 0.70_dp]
   real(dp) :: series(180), innovations(180), expected, transition
   real(dp), allocatable :: initial(:)
   integer :: time

   innovations = 0.0_dp
   series = 0.0_dp
   series(1:2) = [0.4_dp, -0.2_dp]
   do time = 3, size(series)
      innovations(time) = 0.04_dp*sin(0.83_dp*real(time, dp)) + &
         0.02_dp*cos(0.31_dp*real(time, dp))
      transition = exp(-parameters(5)*series(time - 1)**2)
      series(time) = (parameters(1) + parameters(3)*transition)* &
         series(time - 1) + (parameters(2) + parameters(4)*transition)* &
         series(time - 2) + innovations(time)
   end do

   evaluated = expar_evaluate(series, parameters)
   call check(evaluated%info == 0 .and. evaluated%order == 2 .and. &
      maxval(abs(evaluated%residuals(3:) - innovations(3:))) < 1.0e-13_dp .and. &
      abs(evaluated%rss - sum(innovations**2)) < 1.0e-13_dp, &
      'EXPAR conditional fit evaluation')

   initial = expar_initial_parameters(series, 2)
   call check(size(initial) == 5 .and. all(ieee_is_finite(initial)) .and. &
      abs(initial(5) - 0.5_dp) < 1.0e-14_dp, &
      'EXPAR autoregressive initial values')

   fitted = expar_fit(series, 2, parameters, 250, 1.0e-7_dp)
   call check(fitted%info == 0 .and. fitted%rss <= evaluated%rss + 1.0e-10_dp .and. &
      all(ieee_is_finite(fitted%phi)) .and. all(ieee_is_finite(fitted%pi)) .and. &
      ieee_is_finite(fitted%gamma), 'EXPAR finite-difference BFGS estimation')

   forecast = expar_forecast(evaluated, 5)
   transition = exp(-parameters(5)*series(size(series))**2)
   expected = (parameters(1) + parameters(3)*transition)*series(size(series)) + &
      (parameters(2) + parameters(4)*transition)*series(size(series) - 1)
   call check(forecast%info == 0 .and. size(forecast%mean) == 5 .and. &
      abs(forecast%mean(1) - expected) < 1.0e-14_dp .and. &
      all(ieee_is_finite(forecast%mean)), 'recursive EXPAR forecast')

   selection = expar_select(series, 3, 'BIC', 180, 1.0e-6_dp)
   call check(selection%info == 0 .and. selection%selected_order >= 1 .and. &
      selection%selected_order <= 3 .and. trim(selection%criterion) == 'BIC' .and. &
      all(ieee_is_finite(selection%bic)), 'EXPAR information-criterion selection')

   print '(a)', 'EXPAR tests passed'

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

end program test_expar
