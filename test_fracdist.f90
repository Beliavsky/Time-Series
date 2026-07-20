! SPDX-License-Identifier: GPL-3.0-only
! SPDX-FileComment: Regression tests for the fracdist translation.
program test_fracdist
   use kind_mod, only: dp
   use fracdist_mod
   implicit none

   type(fracdist_probability_t) :: probability
   type(fracdist_critical_values_t) :: critical
   type(fracdist_table_t) :: table
   real(dp) :: local

   table = fracdist_get_table(3, 0)
   call check(table%info == 0 .and. all(shape(table%quantile) == [221, 31]), &
      'embedded fracdist table')
   local = fracdist_blocal(0.75_dp, table%quantile(201, :), table%b)
   call check(abs(local - 19.9825543227042_dp) < 2.0e-10_dp, &
      'local fractional-order response surface')

   probability = fracdist_p_value(1, 0, 0.73_dp, 3.84_dp)
   call check(probability%info == 0 .and. &
      abs(probability%p_value - 0.0461_dp) < 1.0e-14_dp, &
      'rank-one fracdist p-value')
   probability = fracdist_p_value(3, 1, 1.27_dp, 32.84_dp)
   call check(probability%info == 0 .and. &
      abs(probability%p_value - 0.1882_dp) < 1.0e-14_dp, &
      'rank-three fracdist p-value')
   probability = fracdist_p_value(12, 1, 1.27_dp, 412.84_dp)
   call check(probability%info == 0 .and. &
      abs(probability%p_value - 0.0320_dp) < 1.0e-14_dp, &
      'rank-twelve fracdist p-value')

   critical = fracdist_critical_values(1, 0, 0.73_dp, [0.05_dp])
   call check(critical%info == 0 .and. &
      abs(critical%critical_value(1) - 3.7066_dp) < 1.0e-14_dp, &
      'rank-one fracdist critical value')
   critical = fracdist_critical_values(3, 1, 1.27_dp, [0.05_dp])
   call check(critical%info == 0 .and. &
      abs(critical%critical_value(1) - 38.7691_dp) < 1.0e-14_dp, &
      'rank-three fracdist critical value')
   critical = fracdist_critical_values(12, 1, 1.27_dp, [0.05_dp])
   call check(critical%info == 0 .and. &
      abs(critical%critical_value(1) - 407.8485_dp) < 1.0e-14_dp, &
      'rank-twelve fracdist critical value')

   probability = fracdist_p_value(3, 1, 0.27_dp, 32.84_dp)
   call check(probability%info == 0 .and. &
      abs(probability%p_value - 0.00014239371259106_dp) < 2.0e-15_dp, &
      'low-order chi-square probability')
   critical = fracdist_critical_values(1, 0, 0.43_dp, [0.05_dp])
   call check(critical%info == 0 .and. &
      abs(critical%critical_value(1) - 3.84145882069412_dp) < 2.0e-12_dp, &
      'low-order chi-square critical value')

   print '(a)', 'fracdist tests passed'

contains

   subroutine check(condition, message)
      !! Stop the test program when a condition fails.
      logical, intent(in) :: condition !! Flag controlling condition.
      character(len=*), intent(in) :: message !! Message.

      if (.not. condition) then
         print '(a)', 'FAILED: '//trim(message)
         error stop 1
      end if
   end subroutine check

end program test_fracdist
