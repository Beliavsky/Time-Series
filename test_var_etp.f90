! SPDX-License-Identifier: GPL-2.0-only
! SPDX-FileComment: Numerical tests for algorithms translated from R VAR.etp.
program test_var_etp
   use kind_mod, only: dp
   use random_mod, only: set_random_seed
   use vars_mod, only: vars_fit_t
   use var_etp_mod
   implicit none

   real(dp) :: series(140, 2), innovation(2)
   type(var_etp_bias_result_t) :: pope, bootstrap
   type(var_etp_forecast_result_t) :: normal_forecast, bootstrap_forecast
   type(var_etp_test_result_t) :: wald, lr
   type(var_etp_predictive_result_t) :: predictive
   type(var_etp_order_result_t) :: predictive_order
   type(vars_fit_t) :: restricted
   integer :: restrictions(1, 3), no_restrictions(0, 3)
   real(dp) :: slope_restriction(1, 2), predictive_forecast(3, 3)
   integer :: t

   series = 0.0_dp
   do t = 2, size(series, 1)
      innovation = [0.12_dp*sin(0.37_dp*real(t, dp)), &
         0.09_dp*cos(0.29_dp*real(t, dp))]
      series(t, 1) = 0.25_dp + 0.55_dp*series(t - 1, 1) - &
         0.12_dp*series(t - 1, 2) + innovation(1)
      series(t, 2) = -0.1_dp + 0.18_dp*series(t - 1, 1) + &
         0.42_dp*series(t - 1, 2) + innovation(2)
   end do

   pope = var_etp_pope(series, 1)
   call assert_true(pope%info == 0, 'Pope correction status')
   call assert_true(all(shape(pope%bias) == [2, 3]), 'Pope bias shape')
   call assert_true(pope%adjustment >= 0.0_dp .and. pope%adjustment <= 1.0_dp, &
      'Pope adjustment range')

   call set_random_seed(418)
   bootstrap = var_etp_bootstrap_bias(series, 1, 20)
   call assert_true(bootstrap%info == 0, 'bootstrap correction status')
   call assert_true(bootstrap%successful_runs == 20, 'bootstrap correction runs')

   normal_forecast = var_etp_forecast(series, 1, 5, level=0.9_dp)
   call assert_true(normal_forecast%info == 0, 'normal forecast status')
   call assert_true(all(shape(normal_forecast%mse) == [2, 2, 5]), 'forecast MSE shape')
   call assert_true(all(normal_forecast%lower <= normal_forecast%upper), &
      'normal interval ordering')

   call set_random_seed(917)
   bootstrap_forecast = var_etp_bootstrap_prediction(series, 1, 4, 20, level=0.9_dp)
   call assert_true(bootstrap_forecast%info == 0, 'bootstrap prediction status')
   call assert_true(bootstrap_forecast%successful_runs == 20, 'bootstrap prediction runs')
   call assert_true(all(bootstrap_forecast%lower <= bootstrap_forecast%upper), &
      'bootstrap interval ordering')

   restrictions(1, :) = [1, 1, 2]
   restricted = var_etp_restrict(series, 1, restrictions)
   call assert_true(restricted%info == 0, 'restricted VAR status')
   call assert_true(abs(restricted%ar(1, 2, 1)) < 1.0e-12_dp, &
      'restricted coefficient value')
   call set_random_seed(121)
   wald = var_etp_wald_test(series, 1, restrictions, bootstrap_method=2, runs=20)
   call assert_true(wald%info == 0, 'wild-bootstrap Wald status')
   call assert_probability(wald%p_value, 'Wald p-value')
   call assert_probability(wald%bootstrap_p_value, 'bootstrap Wald p-value')
   call set_random_seed(122)
   lr = var_etp_lr_test(series, 1, restrictions, no_restrictions, &
      bootstrap_method=1, runs=20)
   call assert_true(lr%info == 0, 'iid-bootstrap LR status')
   call assert_probability(lr%p_value, 'LR p-value')
   call assert_probability(lr%bootstrap_p_value, 'bootstrap LR p-value')

   slope_restriction = reshape([1.0_dp, 0.0_dp], [1, 2])
   predictive = var_etp_predictive_regression(series(:, 1:2), &
      0.3_dp*series(:, 1) - 0.1_dp*series(:, 2), 1, slope_restriction)
   call assert_true(predictive%info == 0, 'predictive regression status')
   call assert_true(size(predictive%improved_coefficients) == 5, &
      'improved predictive coefficient count')
   call assert_probability(predictive%improved_test%p_value, &
      'improved predictive p-value')
   predictive_order = var_etp_predictive_order(series(:, 1:2), &
      0.3_dp*series(:, 1) - 0.1_dp*series(:, 2), 3)
   call assert_true(predictive_order%info == 0, 'predictive order status')
   call assert_true(predictive_order%bic_order >= 1 .and. &
      predictive_order%bic_order <= 3, 'predictive selected order')
   predictive_forecast = var_etp_predictive_forecast(series(:, 1:2), &
      0.3_dp*series(:, 1) - 0.1_dp*series(:, 2), predictive, 3)
   call assert_true(all(abs(predictive_forecast) < huge(1.0_dp)), &
      'predictive forecast values')

   print '(a)', 'All VAR.etp tests passed.'

contains

   subroutine assert_true(condition, label)
      !! Stop when a logical test condition is false.
      logical, intent(in) :: condition !! Condition expected to be true.
      character(len=*), intent(in) :: label !! Failure label.

      if (.not. condition) error stop 'FAIL: '//label
   end subroutine assert_true

   subroutine assert_probability(value, label)
      !! Stop when a value is outside the unit interval.
      real(dp), intent(in) :: value !! Candidate probability.
      character(len=*), intent(in) :: label !! Failure label.

      if (value < 0.0_dp .or. value > 1.0_dp) error stop 'FAIL: '//label
   end subroutine assert_probability

end program test_var_etp
