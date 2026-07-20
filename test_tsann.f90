! SPDX-License-Identifier: GPL-3.0-only
! SPDX-FileComment: Tests for algorithms translated from the R TSANN package.
program test_tsann
   use kind_mod, only: dp
   use tsann_mod, only: tsann_fit_t, tsann_auto_fit, tsann_maximum_lag, &
      tsann_select_validation, display
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none

   real(dp) :: series(60)
   type(tsann_fit_t) :: fit
   integer :: time, failures

   failures = 0
   series(1:2) = [0.25_dp, -0.1_dp]
   do time = 3, size(series)
      series(time) = 0.65_dp*series(time - 1) - &
         0.2_dp*series(time - 2) + 0.08_dp*sin(0.7_dp*real(time, dp))
   end do

   call check(tsann_maximum_lag(series, lag_max=8) >= 1, &
      'ACF-derived maximum lag')
   fit = tsann_auto_fit(series, 1, 2, 0.75_dp, maximum_lag=2, &
      validation_fraction=0.2_dp, repetitions=1, max_iterations=80, &
      tolerance=1.0e-5_dp, decay=1.0e-4_dp, seed=2718)
   call check(fit%info == 0, 'automatic TSANN fit status')
   call check(fit%selection_metric == tsann_select_validation, &
      'validation selection default')
   call check(fit%selected_lag >= 1 .and. fit%selected_lag <= 2 .and. &
      fit%selected_hidden >= 1 .and. fit%selected_hidden <= 2, &
      'selected grid coordinates')
   call check(size(fit%trace%lag_order) == 4 .and. &
      size(fit%trace%hidden_size) == 4, 'complete search trace')
   call check(size(fit%predicted) == 15 .and. size(fit%fitted) == 60, &
      'chronological fitted and test shapes')
   call check(all(ieee_is_finite(fit%predicted)) .and. &
      ieee_is_finite(fit%train_rmse) .and. &
      ieee_is_finite(fit%validation_rmse) .and. &
      ieee_is_finite(fit%test_rmse), 'finite TSANN accuracy results')
   call display(fit)

   if (failures > 0) error stop 'TSANN tests failed'
   print '(a)', 'TSANN tests passed'

contains

   subroutine check(condition, label)
      !! Record a failed logical test.
      logical, intent(in) :: condition !! Test condition.
      character(len=*), intent(in) :: label !! Test label.

      if (condition) return
      failures = failures + 1
      print '(a)', 'FAILED: '//trim(label)
   end subroutine check

end program test_tsann
