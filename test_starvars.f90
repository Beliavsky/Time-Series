! SPDX-License-Identifier: MIT
! SPDX-FileComment: Tests for algorithms translated from the R starvars package.
program test_starvars
   use kind_mod, only: dp
   use calendar_mod, only: date_t
   use starvars_mod, only: starvars_start_t, starvars_fit_t, &
      starvars_forecast_t, starvars_joint_test_t, &
      starvars_long_run_variance_t, starvars_cumsum_t, &
      starvars_realized_covariance_t, starvars_method_nls, &
      starvars_forecast_naive, starvars_forecast_monte_carlo, &
      starvars_forecast_bootstrap, starvars_logistic, starvars_starting, &
      starvars_frequency_monthly, starvars_fit, &
      starvars_forecast, starvars_joint_linearity_test, &
      starvars_long_run_variance, starvars_multiple_cumsum, &
      starvars_realized_covariance
   implicit none
   integer, parameter :: observations = 90
   real(dp) :: series(observations, 2), transition(observations)
   real(dp) :: candidates(observations, 2), prices(12, 2)
   integer :: groups(12), time
   type(date_t) :: dates(12)
   type(starvars_start_t) :: starting, automatic_start
   type(starvars_fit_t) :: fit
   type(starvars_forecast_t) :: forecast
   type(starvars_joint_test_t) :: joint
   type(starvars_long_run_variance_t) :: long_run
   type(starvars_cumsum_t) :: breaks
   type(starvars_realized_covariance_t) :: realized

   series = 0.0_dp
   transition = [(sin(0.13_dp*real(time, dp)), time=1, observations)]
   do time = 2, observations
      series(time, 1) = 0.2_dp + 0.45_dp*series(time - 1, 1) - &
         0.12_dp*series(time - 1, 2) + &
         0.08_dp*sin(0.71_dp*real(time, dp))
      series(time, 2) = -0.1_dp + 0.18_dp*series(time - 1, 1) + &
         0.35_dp*series(time - 1, 2) + &
         0.07_dp*cos(0.53_dp*real(time, dp))
   end do
   allocate(starting%gamma(2, 1), starting%location(2, 1))
   starting%gamma = 1.0_dp
   starting%location = 0.0_dp
   automatic_start = starvars_starting(series, transition, combinations=3)
   if (automatic_start%info /= 0) error stop 'starting values failed'
   fit = starvars_fit(series, transition, order=1, regimes=2, &
      starting=starting, method=starvars_method_nls, max_iterations=40)
   if (fit%info /= 0) then
      error stop 'starvars_fit failed'
   end if
   if (size(fit%coefficients, 1) /= 3) error stop 'coefficient shape failed'
   if (abs(starvars_logistic(0.0_dp, 1.0_dp, 0.0_dp) - 0.5_dp) > &
      1.0e-12_dp) error stop 'logistic failed'

   forecast = starvars_forecast(fit, 4, method=starvars_forecast_naive)
   if (forecast%info /= 0) error stop 'starvars_forecast failed'
   if (size(forecast%point, 1) /= 4) error stop 'forecast shape failed'
   forecast = starvars_forecast(fit, 3, method=starvars_forecast_monte_carlo, &
      simulations=20, seed=123, keep_paths=.true.)
   if (forecast%info /= 0 .or. .not. allocated(forecast%paths)) &
      error stop 'Monte Carlo forecast failed'
   forecast = starvars_forecast(fit, 3, method=starvars_forecast_bootstrap, &
      simulations=20, seed=321)
   if (forecast%info /= 0 .or. allocated(forecast%paths)) &
      error stop 'bootstrap forecast failed'

   candidates(:, 1) = transition
   candidates(:, 2) = [(cos(0.09_dp*real(time, dp)), time=1, observations)]
   joint = starvars_joint_linearity_test(series, candidates, order=1)
   if (joint%info /= 0) error stop 'joint test failed'
   if (any(joint%p_value < 0.0_dp) .or. any(joint%p_value > 1.0_dp)) &
      error stop 'joint p-value failed'

   long_run = starvars_long_run_variance(series(:, 1))
   if (long_run%info /= 0 .or. long_run%variance <= 0.0_dp) &
      error stop 'long-run variance failed'

   breaks = starvars_multiple_cumsum(series, max_breaks=2)
   if (breaks%info /= 0) error stop 'multiple CUMSUM failed'
   if (size(breaks%lambda) /= 2) error stop 'CUMSUM shape failed'

   do time = 1, 12
      prices(time, 1) = 100.0_dp*exp(0.01_dp*real(time, dp))
      prices(time, 2) = 80.0_dp*exp(0.008_dp*real(time, dp) + &
         0.002_dp*sin(real(time, dp)))
      groups(time) = 1 + (time - 1)/4
      dates(time) = date_t(2026, 1 + (time - 1)/6, 1 + mod(time - 1, 6))
   end do
   realized = starvars_realized_covariance(prices, groups)
   if (realized%info /= 0) error stop 'realized covariance failed'
   if (size(realized%covariance, 1) /= 3 .or. &
      size(realized%covariance, 2) /= 3) &
      error stop 'realized covariance shape failed'
   realized = starvars_realized_covariance(prices, dates, &
      frequency=starvars_frequency_monthly)
   if (realized%info /= 0 .or. size(realized%covariance, 1) /= 2) &
      error stop 'date-indexed realized covariance failed'

   print '(a)', 'starvars tests passed'
end program test_starvars
