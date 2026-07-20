! SPDX-License-Identifier: GPL-3.0-only
! SPDX-FileComment: Regression tests for the forecast translation.
program test_forecast
   use kind_mod, only: dp
   use forecast_mod
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
   implicit none
   real(dp) :: y(6) = [1, 2, 3, 4, 5, 6], demand(8) = [0, 2, 0, 0, 4, 0, 0, 0]
   real(dp), allocatable :: f(:, :), a(:), p(:), ma(:)
   type(forecast_result_t) :: fit
   type(accuracy_result_t) :: acc
   fit = meanf(y, 3)
   call check(maxval(abs(fit%mean - 3.5_dp)) < 1e-12, 'meanf')
   fit = rwf(y, 2, .true.)
   call check(maxval(abs(fit%mean - [7._dp, 8._dp])) < 1e-12, 'rwf drift')
   fit = snaive(y, 3, 4)
   call check(maxval(abs(fit%mean - [4._dp, 5._dp, 6._dp, 4._dp])) < 1e-12, 'snaive')
   fit = ses(y, 2, .5_dp)
   call check(abs(fit%mean(1) - 5.03125_dp) < 1e-12, 'ses')
   fit = croston(demand, 2, .1_dp, 'sba')
   call check(all(fit%mean > 0), 'croston')
   f = fourier(4, 4, 1)
   call check(abs(f(1, 1) - 1._dp) < 1e-12, 'fourier')
   ma = moving_average(y, 3)
   call check(abs(ma(3) - 3._dp) < 1e-12 .and. ieee_is_nan(ma(1)), 'moving average')
   a = acf_values(y, 2)
   call check(abs(a(1) - 1._dp) < 1e-12, 'acf')
   p = pacf_values(y, 2)
   call check(size(p) == 2, 'pacf')
   acc = forecast_accuracy(y, y)
   call check(abs(acc%rmse) < epsilon(1.0_dp), 'accuracy')
   print '(a)', 'All forecast_mod tests passed.'
contains
   subroutine check(ok, name)
      !! Stop the test program when a named assertion fails.
      logical, intent(in) :: ok !! Flag controlling ok.
      character(*), intent(in) :: name !! Name.
      if (.not. ok) then
         print '(a)', 'FAILED: '//name
         error stop 1
      end if
   end subroutine check
end program
