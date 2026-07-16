! SPDX-License-Identifier: GPL-2.0-only
! SPDX-FileComment: Regression tests for the bayesforecast translation.
program test_bayesforecast
   use kind_mod, only: dp
   use bayesforecast_mod
   use time_series_random_mod, only: set_random_seed
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   type(bf_filter_t) :: sarima, garch
   type(bf_ets_filter_t) :: ets
   type(bf_sv_filter_t) :: sv
   type(bf_mcmc_t) :: sampler
   type(bf_information_criteria_t) :: criteria
   type(bf_difference_t) :: differenced
   type(bf_predictive_t) :: prediction
   real(dp) :: series(8), empty(0), seasonal_empty(0), log_variance(8)
   real(dp) :: expected, normal_prior(4), uniform_prior(4)
   real(dp) :: simple_draws(5, 2), log_likelihood_draws(3, 2)
   real(dp), allocatable :: interval(:, :), prediction_error(:, :)
   real(dp), allocatable :: restored(:)
   integer :: t

   normal_prior = [0.0_dp, 2.0_dp, 0.0_dp, real(bf_prior_normal, dp)]
   uniform_prior = [-1.0_dp, 3.0_dp, 0.0_dp, real(bf_prior_uniform, dp)]
   call check(abs(bf_log_prior(0.0_dp, normal_prior) + &
      0.5_dp*log(8.0_dp*acos(-1.0_dp))) < 1.0e-12_dp, 'normal prior density')
   call check(abs(bf_log_prior(1.0_dp, uniform_prior) + log(4.0_dp)) < 1.0e-12_dp, &
      'uniform prior density')

   series(1) = 0.2_dp
   do t = 2, 8
      series(t) = 0.1_dp + 0.5_dp*series(t - 1) + 0.2_dp*sin(real(t, dp))
   end do
   sarima = bf_sarima_filter(series, 0.4_dp, 0.1_dp, [0.5_dp], empty, &
      seasonal_empty, seasonal_empty, 1)
   expected = 0.1_dp + 0.5_dp*series(1)
   call check(sarima%info == 0 .and. abs(sarima%mean(2) - expected) < 1.0e-12_dp .and. &
      maxval(abs(sarima%residual - series + sarima%mean)) < 1.0e-12_dp .and. &
      abs(sarima%total_log_likelihood - sum(sarima%log_likelihood)) < 1.0e-12_dp, &
      'Bayesian SARIMA recursion')

   ets = bf_ets_filter(series, 0.5_dp, 0.3_dp, series(1))
   expected = 0.3_dp*series(2) + 0.7_dp*series(1)
   call check(ets%info == 0 .and. abs(ets%level(2) - expected) < 1.0e-12_dp .and. &
      abs(ets%observation%mean(2) - series(1)) < 1.0e-12_dp .and. &
      all(ieee_is_finite(ets%observation%log_likelihood)), 'Bayesian ETS recursion')

   garch = bf_garch_filter(series, 0.2_dp, 0.0_dp, empty, empty, [0.15_dp], &
      [0.7_dp], empty, 0, 0.0_dp, 0.0_dp)
   expected = sqrt(0.2_dp + 0.15_dp*series(1)**2 + 0.7_dp*0.2_dp)
   call check(garch%info == 0 .and. abs(garch%scale(2) - expected) < 1.0e-12_dp .and. &
      all(garch%scale > 0.0_dp) .and. ieee_is_finite(garch%total_log_likelihood), &
      'Bayesian GARCH recursion')

   log_variance = log(0.4_dp)
   sv = bf_stochastic_volatility_filter(series, log_variance, log(0.4_dp), &
      0.8_dp, 0.2_dp, 0.0_dp, empty, empty)
   call check(sv%info == 0 .and. maxval(abs(sv%observation%scale - sqrt(0.4_dp))) < 1.0e-12_dp .and. &
      abs(sv%total_log_likelihood - sv%observation%total_log_likelihood - &
      sum(sv%state_log_likelihood)) < 1.0e-12_dp, 'Bayesian stochastic-volatility density')

   call set_random_seed(24680)
   sampler = bf_metropolis_sample(standard_normal_posterior, [0.0_dp], [0.8_dp], 1000, 200, 2)
   call check(sampler%info == 0 .and. all(shape(sampler%draws) == [1000, 1]) .and. &
      sampler%acceptance_rate > 0.2_dp .and. sampler%acceptance_rate < 0.9_dp .and. &
      abs(sum(sampler%draws(:, 1))/1000.0_dp) < 0.25_dp, 'Metropolis posterior sampling')
   simple_draws(:, 1) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
   simple_draws(:, 2) = 2.0_dp*simple_draws(:, 1)
   interval = bf_posterior_interval(simple_draws, 0.8_dp)
   call check(all(shape(interval) == [2, 2]) .and. interval(1, 1) < interval(1, 2) .and. &
      interval(2, 1) == 2.0_dp*interval(1, 1), 'posterior equal-tail intervals')
   prediction_error = bf_predictive_error(simple_draws(:, :2), [3.0_dp, 6.0_dp])
   call check(prediction_error(3, 1) == 0.0_dp .and. prediction_error(3, 2) == 0.0_dp, &
      'posterior predictive errors')
   log_likelihood_draws = reshape([-1.0_dp, -1.1_dp, -0.9_dp, &
      -2.0_dp, -2.2_dp, -1.8_dp], [3, 2])
   criteria = bf_waic(log_likelihood_draws)
   call check(criteria%info == 0 .and. criteria%effective_parameters > 0.0_dp .and. &
      ieee_is_finite(criteria%waic), 'widely applicable information criterion')
   criteria = bf_information_criteria(-20.0_dp, 3, 50)
   call check(criteria%info == 0 .and. criteria%aic == 46.0_dp .and. &
      criteria%aicc > criteria%aic .and. criteria%bic > criteria%aic, &
      'classical information criteria')
   call check(abs(bf_bayes_factor(-10.0_dp, -11.0_dp) - exp(1.0_dp)) < 1.0e-12_dp .and. &
      bf_bayes_factor(-10.0_dp, -11.0_dp, .true.) == 1.0_dp, 'Bayes factor conversion')
   call check(ieee_is_finite(bf_sarima_log_posterior(series, 0.4_dp, 0.1_dp, &
      [0.5_dp], empty, seasonal_empty, seasonal_empty, 1, [0.4_dp, 0.1_dp], &
      reshape([0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp], [2, 4]))), &
      'SARIMA log posterior')
   differenced = bf_difference([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp], 1, 0, 1)
   restored = bf_inverse_difference([1.0_dp, 1.0_dp], differenced)
   call check(differenced%info == 0 .and. all(restored == [6.0_dp, 7.0_dp]), &
      'ordinary forecast inverse differencing')
   call set_random_seed(13579)
   prediction = bf_sarima_predict(series, 0.1_dp, 0.1_dp, [0.5_dp], empty, &
      seasonal_empty, seasonal_empty, 1, 2, 2000, 0.9_dp)
   expected = 0.1_dp + 0.5_dp*series(8)
   call check(prediction%info == 0 .and. all(shape(prediction%draws) == [2000, 2]) .and. &
      all(shape(prediction%interval) == [2, 2]) .and. abs(prediction%mean(1) - expected) < 0.02_dp, &
      'SARIMA posterior predictive simulation')
   call set_random_seed(97531)
   prediction = bf_ets_predict(series, 0.1_dp, 0.3_dp, series(1), 2, 1000, 0.9_dp)
   call check(prediction%info == 0 .and. all(shape(prediction%draws) == [1000, 2]) .and. &
      all(ieee_is_finite(prediction%mean)), 'ETS posterior predictive simulation')
   prediction = bf_garch_predict(series, 0.2_dp, 0.0_dp, empty, empty, [0.15_dp], &
      [0.7_dp], empty, 0, 0.0_dp, 0.0_dp, 2, 1000, 0.9_dp)
   call check(prediction%info == 0 .and. all(shape(prediction%interval) == [2, 2]) .and. &
      all(ieee_is_finite(prediction%draws)), 'GARCH posterior predictive simulation')
   prediction = bf_stochastic_volatility_predict(series, log_variance, log(0.4_dp), &
      0.8_dp, 0.2_dp, 0.0_dp, empty, empty, 2, 1000, 0.9_dp)
   call check(prediction%info == 0 .and. prediction%horizon == 2 .and. &
      all(ieee_is_finite(prediction%draws)), 'stochastic-volatility predictive simulation')
   print '(a)', 'All bayesforecast_mod tests passed.'

contains

   pure function standard_normal_posterior(parameters) result(value)
      ! Evaluate an unnormalized standard-normal posterior.
      real(dp), intent(in) :: parameters(:)
      real(dp) :: value

      value = -0.5_dp*dot_product(parameters, parameters)
   end function standard_normal_posterior

   subroutine check(ok, name)
      ! Stop the test program when a named assertion fails.
      logical, intent(in) :: ok
      character(len=*), intent(in) :: name

      if (.not. ok) then
         print '(a)', 'FAILED: '//name
         error stop 1
      end if
   end subroutine check
end program test_bayesforecast
