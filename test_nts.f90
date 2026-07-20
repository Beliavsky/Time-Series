! SPDX-License-Identifier: GPL-2.0-or-later
! SPDX-FileComment: Numerical tests for algorithms translated from R NTS.
program test_nts
   use kind_mod, only: dp
   use nts_mod
   use random_mod, only: set_random_seed
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none

   real(dp) :: intercept(2), ar(2, 1), thresholds(1), innovation_sd(2)
   integer :: ar_order(2)
   type(nts_tar_simulation_t) :: simulation, repeated
   type(nts_tar_model_t) :: fit
   type(nts_tar_search_t) :: search
   type(nts_tar_forecast_t) :: forecast
   type(nts_msar_simulation_t) :: msar_simulation, repeated_msar
   type(nts_nonlinearity_test_t) :: threshold_test, tsay_test
   type(nts_nonlinearity_test_t) :: quadratic_test
   type(nts_tar_backtest_t) :: backtest, repeated_backtest
   type(nts_rank_portmanteau_t) :: rank_test
   type(nts_prnd_test_t) :: prnd_test
   type(nts_tvar_fit_t) :: tvar_fit
   type(nts_rcar_fit_t) :: rcar_fit
   type(nts_mtar_refinement_t) :: multivariate_refinement
   type(nts_acmx_fit_t) :: poisson_acmx, negative_binomial_acmx
   type(nts_acmx_fit_t) :: double_poisson_acmx
   type(nts_cfar_simulation_t) :: cfar_simulation, repeated_cfar
   type(nts_cfar_irregular_simulation_t) :: irregular_cfar_simulation
   type(nts_cfar_model_t) :: cfar_fit, irregular_cfar_fit
   type(nts_cfar_order_tests_t) :: cfar_tests, irregular_cfar_tests
   type(nts_cfar_forecast_t) :: cfar_forecast, cfar_partial
   type(nts_smc_filter_t) :: smc_filter, repeated_smc
   type(nts_smc_smoother_t) :: smc_smoother
   type(nts_smc_marginal_smoother_t) :: smc_marginal_smoother
   type(nts_rb_smc_filter_t) :: rb_smc_filter
   integer, parameter :: smc_particles = 80, smc_times = 40
   real(dp) :: smc_observations(smc_times, 1)
   real(dp) :: smc_initial(1, smc_particles)
   real(dp) :: smc_draws(1, smc_particles, smc_times)
   real(dp) :: smc_uniforms(smc_particles, smc_times)
   real(dp) :: rb_initial_mean(1, smc_particles)
   real(dp) :: rb_initial_covariance(1, 1, smc_particles)
   logical :: smc_schedule(smc_times)
   real(dp) :: cfar_kernel(21), coordinate
   real(dp) :: irregular_values(90, 6), irregular_positions(90, 6)
   integer :: irregular_counts(90)
   integer, parameter :: count_observations = 320
   real(dp) :: count_series(count_observations)
   real(dp) :: count_exogenous(count_observations, 1)
   real(dp) :: latent_mean(count_observations), conditional_count_mean
   integer :: time, particle
   real(dp) :: msar_transition(2, 2)
   real(dp) :: multivariate_intercept(2, 2)
   real(dp) :: multivariate_ar(2, 2, 1, 2)
   real(dp) :: multivariate_covariance(2, 2, 2)
   type(nts_mtar_simulation_t) :: multivariate_simulation
   type(nts_mtar_simulation_t) :: repeated_multivariate
   type(nts_mtar_model_t) :: multivariate_fit
   type(nts_mtar_search_t) :: multivariate_search
   type(nts_mtar_forecast_t) :: multivariate_forecast

   intercept = [0.20_dp, -0.20_dp]
   ar(:, 1) = [0.60_dp, -0.50_dp]
   ar_order = [1, 1]
   thresholds = 0.0_dp
   innovation_sd = [0.20_dp, 0.25_dp]

   call set_random_seed(9137)
   simulation = nts_utar_simulate(intercept, ar, ar_order, thresholds, &
      innovation_sd, 1, 1200, 500)
   call check(simulation%info == 0, 'SETAR simulation status')
   call check(size(simulation%series) == 1200, 'SETAR simulation size')
   call check(all(ieee_is_finite(simulation%series)), 'SETAR simulation finite')

   call set_random_seed(9137)
   repeated = nts_utar_simulate(intercept, ar, ar_order, thresholds, &
      innovation_sd, 1, 1200, 500)
   call check(maxval(abs(repeated%series - simulation%series)) < 1.0e-14_dp, &
      'SETAR shared RNG reproducibility')

   fit = nts_utar_estimate(simulation%series, ar_order, thresholds, 1)
   call check(fit%info == 0, 'fixed-threshold SETAR fit status')
   call check(all(fit%regime_observations > 100), 'SETAR regime sample sizes')
   call check(maxval(abs(fit%coefficients(:, 1) - intercept)) < 0.08_dp, &
      'SETAR intercept recovery')
   call check(maxval(abs(fit%coefficients(:, 2) - ar(:, 1))) < 0.10_dp, &
      'SETAR AR recovery')
   call check(maxval(abs(fit%innovation_sd - innovation_sd)) < 0.04_dp, &
      'SETAR scale recovery')

   search = nts_utar_threshold_search(simulation%series, 1, 1, 1, &
      trim=[0.20_dp, 0.80_dp])
   call check(search%info == 0, 'SETAR threshold search status')
   call check(search%selected >= 1, 'SETAR selected threshold index')
   call check(abs(search%model%thresholds(1)) < 0.15_dp, &
      'SETAR threshold recovery')

   call set_random_seed(3181)
   forecast = nts_utar_forecast(fit, size(fit%data), 4, 500, 0.90_dp)
   call check(forecast%info == 0, 'SETAR forecast status')
   call check(all(shape(forecast%simulations) == [4, 500]), &
      'SETAR forecast simulation shape')
   call check(all(forecast%lower <= forecast%mean .and. &
      forecast%mean <= forecast%upper), 'SETAR forecast interval ordering')

   msar_transition = reshape([0.95_dp, 0.10_dp, 0.05_dp, 0.90_dp], [2, 2])
   call set_random_seed(4173)
   msar_simulation = nts_msar_simulate(intercept, ar, ar_order, &
      msar_transition, innovation_sd, 1600, 500)
   call check(msar_simulation%info == 0, 'Markov-switching AR simulation status')
   call check(size(msar_simulation%series) == 1600 .and. &
      size(msar_simulation%state) == 1600, &
      'Markov-switching AR simulation dimensions')
   call check(count(msar_simulation%state == 1) > 200 .and. &
      count(msar_simulation%state == 2) > 200, &
      'Markov-switching AR regime representation')
   call set_random_seed(4173)
   repeated_msar = nts_msar_simulate(intercept, ar, ar_order, &
      msar_transition, innovation_sd, 1600, 500)
   call check(maxval(abs(repeated_msar%series - msar_simulation%series)) < &
      1.0e-14_dp .and. all(repeated_msar%state == msar_simulation%state), &
      'Markov-switching AR shared RNG reproducibility')

   threshold_test = nts_threshold_test(simulation%series, 1, 1, &
      initial_count=40, include_mean=.true.)
   call check(threshold_test%info == 0, 'threshold nonlinearity test status')
   call check(ieee_is_finite(threshold_test%statistic) .and. &
      threshold_test%statistic >= 0.0_dp .and. &
      threshold_test%p_value >= 0.0_dp .and. threshold_test%p_value <= 1.0_dp, &
      'threshold nonlinearity test values')
   call check(threshold_test%p_value < 0.05_dp, &
      'threshold nonlinearity detection')
   tsay_test = nts_tsay_test(simulation%series, 1)
   call check(tsay_test%info == 0, 'Tsay nonlinearity test status')
   call check(ieee_is_finite(tsay_test%statistic) .and. &
      tsay_test%statistic >= 0.0_dp .and. &
      tsay_test%p_value >= 0.0_dp .and. tsay_test%p_value <= 1.0_dp, &
      'Tsay nonlinearity test values')

   call set_random_seed(9231)
   backtest = nts_tar_backtest(fit, 1190, 2, 50)
   call check(backtest%info == 0, 'SETAR rolling backtest status')
   call check(all(shape(backtest%error) == [10, 2]) .and. &
      all(ieee_is_finite(backtest%rmse)), 'SETAR rolling backtest results')
   call set_random_seed(9231)
   repeated_backtest = nts_tar_backtest(fit, 1190, 2, 50)
   call check(maxval(abs(repeated_backtest%error - backtest%error)) < 1.0e-14_dp, &
      'SETAR rolling backtest RNG reproducibility')

   rank_test = nts_rank_portmanteau(simulation%series, 10)
   call check(rank_test%info == 0 .and. size(rank_test%statistic) == 10, &
      'rank portmanteau status')
   call check(all(rank_test%p_value >= 0.0_dp .and. &
      rank_test%p_value <= 1.0_dp), 'rank portmanteau p-values')
   quadratic_test = nts_quadratic_f_test(simulation%series, 2, 0.0_dp)
   call check(quadratic_test%info == 0 .and. &
      ieee_is_finite(quadratic_test%statistic), 'quadratic F test')
   prnd_test = nts_prnd_test(simulation%series, 10, 0, 0)
   call check(prnd_test%info == 0 .and. ieee_is_finite(prnd_test%statistic) .and. &
      prnd_test%p_value >= 0.0_dp .and. prnd_test%p_value <= 1.0_dp, &
      'Pena-Rodriguez determinant test')

   tvar_fit = nts_tvar_fit(simulation%series(:240), [1], .true., 60, 1.0e-4_dp)
   call check(tvar_fit%info == 0, 'time-varying AR fit status')
   call check(all(shape(tvar_fit%smoothed_coefficients) == [239, 2]) .and. &
      tvar_fit%observation_variance > 0.0_dp .and. &
      all(tvar_fit%state_variance > 0.0_dp), 'time-varying AR fit results')
   rcar_fit = nts_rcar_fit(simulation%series(:400), [1], .true., 100, 1.0e-5_dp)
   call check(rcar_fit%info == 0, 'random-coefficient AR fit status')
   call check(size(rcar_fit%residuals) == 399 .and. &
      rcar_fit%innovation_variance > 0.0_dp .and. &
      all(rcar_fit%coefficient_variance > 0.0_dp), &
      'random-coefficient AR fit results')

   multivariate_intercept(:, 1) = [0.15_dp, -0.10_dp]
   multivariate_intercept(:, 2) = [-0.15_dp, 0.10_dp]
   multivariate_ar(:, :, 1, 1) = reshape([ &
      0.45_dp, -0.05_dp, &
      0.10_dp, 0.35_dp], [2, 2])
   multivariate_ar(:, :, 1, 2) = reshape([ &
      -0.35_dp, 0.12_dp, &
      0.08_dp, 0.25_dp], [2, 2])
   multivariate_covariance(:, :, 1) = reshape([ &
      0.12_dp, 0.03_dp, &
      0.03_dp, 0.10_dp], [2, 2])
   multivariate_covariance(:, :, 2) = reshape([ &
      0.18_dp, -0.02_dp, &
      -0.02_dp, 0.14_dp], [2, 2])

   call set_random_seed(5179)
   multivariate_simulation = nts_mtar_simulate(multivariate_intercept, &
      multivariate_ar, ar_order, multivariate_covariance, thresholds, &
      1, 1, 1800, 500)
   call check(multivariate_simulation%info == 0, &
      'multivariate SETAR simulation status')
   call check(all(shape(multivariate_simulation%series) == [1800, 2]), &
      'multivariate SETAR simulation shape')
   call check(all(multivariate_simulation%regime_observations > 200), &
      'multivariate SETAR regime sample sizes')

   call set_random_seed(5179)
   repeated_multivariate = nts_mtar_simulate(multivariate_intercept, &
      multivariate_ar, ar_order, multivariate_covariance, thresholds, &
      1, 1, 1800, 500)
   call check(maxval(abs(repeated_multivariate%series - &
      multivariate_simulation%series)) < 1.0e-14_dp, &
      'multivariate SETAR RNG reproducibility')

   multivariate_fit = nts_mtar_estimate(multivariate_simulation%series, &
      ar_order, thresholds, 1, 1)
   call check(multivariate_fit%info == 0, &
      'fixed-threshold multivariate SETAR fit status')
   call check(maxval(abs(multivariate_fit%intercept - &
      multivariate_intercept)) < 0.07_dp, &
      'multivariate SETAR intercept recovery')
   call check(maxval(abs(multivariate_fit%ar - multivariate_ar)) < 0.08_dp, &
      'multivariate SETAR AR recovery')
   call check(maxval(abs(multivariate_fit%covariance - &
      multivariate_covariance)) < 0.035_dp, &
      'multivariate SETAR covariance recovery')

   multivariate_search = nts_mtar_threshold_search( &
      multivariate_simulation%series, 1, 1, 1, 1, &
      trim=[0.20_dp, 0.80_dp], criterion='aic')
   call check(multivariate_search%info == 0, &
      'multivariate SETAR threshold search status')
   call check(abs(multivariate_search%model%thresholds(1)) < 0.15_dp, &
      'multivariate SETAR threshold recovery')

   call set_random_seed(6103)
   multivariate_forecast = nts_mtar_forecast(multivariate_fit, &
      size(multivariate_fit%data, 1), 3, 300, 0.90_dp)
   call check(multivariate_forecast%info == 0, &
      'multivariate SETAR forecast status')
   call check(all(shape(multivariate_forecast%simulations) == [3, 2, 300]), &
      'multivariate SETAR forecast shape')
   call check(all(multivariate_forecast%lower <= multivariate_forecast%mean .and. &
      multivariate_forecast%mean <= multivariate_forecast%upper), &
      'multivariate SETAR forecast interval ordering')

   multivariate_refinement = nts_mtar_refine(multivariate_fit, 0.5_dp)
   call check(multivariate_refinement%info == 0, &
      'multivariate SETAR coefficient refinement status')
   call check(count(multivariate_refinement%retained) > 0 .and. &
      all(ieee_is_finite(multivariate_refinement%model%covariance)), &
      'multivariate SETAR coefficient refinement results')

   count_series(1) = 2.0_dp
   latent_mean(1) = 2.0_dp
   count_exogenous(1, 1) = 0.0_dp
   do time = 2, count_observations
      count_exogenous(time, 1) = 0.2_dp*sin(0.07_dp*real(time, dp))
      latent_mean(time) = 0.8_dp + 0.18_dp*count_series(time - 1) + &
         0.48_dp*latent_mean(time - 1)
      conditional_count_mean = exp(0.35_dp*count_exogenous(time, 1))* &
         latent_mean(time)
      count_series(time) = real(max(0, nint(conditional_count_mean + &
         0.55_dp*sqrt(conditional_count_mean)*sin(1.31_dp*real(time, dp)))), dp)
   end do
   poisson_acmx = nts_acmx_fit(count_series, 1, 1, 'poisson', &
      count_exogenous, max_iterations=120, tolerance=1.0e-4_dp)
   call check(poisson_acmx%info == 0 .and. &
      all(ieee_is_finite(poisson_acmx%conditional_mean)) .and. &
      all(poisson_acmx%conditional_mean > 0.0_dp), 'Poisson ACMx fit')
   negative_binomial_acmx = nts_acmx_fit(count_series, 1, 1, &
      'negative_binomial', count_exogenous, max_iterations=120, &
      tolerance=1.0e-4_dp)
   call check(negative_binomial_acmx%info == 0 .and. &
      negative_binomial_acmx%dispersion > 0.0_dp .and. &
      ieee_is_finite(negative_binomial_acmx%log_likelihood), &
      'negative-binomial ACMx fit')
   double_poisson_acmx = nts_acmx_fit(count_series, 1, 1, &
      'double_poisson', count_exogenous, max_iterations=120, &
      tolerance=1.0e-4_dp)
   call check(double_poisson_acmx%info == 0 .and. &
      double_poisson_acmx%dispersion > 0.0_dp .and. &
      ieee_is_finite(double_poisson_acmx%log_likelihood), &
      'double-Poisson ACMx fit')

   do time = 1, size(cfar_kernel)
      coordinate = real(time - 11, dp)/10.0_dp
      cfar_kernel(time) = 0.45_dp*exp(-8.0_dp*coordinate*coordinate)
   end do
   call set_random_seed(7301)
   cfar_simulation = nts_cfar1_simulate(cfar_kernel, 2.0_dp, 0.25_dp, 90, 80)
   call check(cfar_simulation%info == 0 .and. &
      all(shape(cfar_simulation%series) == [90, 11]) .and. &
      all(ieee_is_finite(cfar_simulation%series)), 'CFAR simulation')
   call set_random_seed(7301)
   repeated_cfar = nts_cfar1_simulate(cfar_kernel, 2.0_dp, 0.25_dp, 90, 80)
   call check(maxval(abs(repeated_cfar%series - cfar_simulation%series)) < &
      1.0e-14_dp, 'CFAR shared RNG reproducibility')
   call set_random_seed(8107)
   irregular_cfar_simulation = nts_cfar_irregular_simulate( &
      reshape(cfar_kernel, [21, 1]), 2.0_dp, 0.25_dp, 40, 4, 2.0_dp, 50)
   call check(irregular_cfar_simulation%info == 0 .and. &
      all(irregular_cfar_simulation%counts >= 4) .and. &
      all(irregular_cfar_simulation%counts <= 11), &
      'heterogeneous CFAR simulation')
   cfar_fit = nts_cfar_estimate(cfar_simulation%series, 1, 3)
   call check(cfar_fit%info == 0 .and. cfar_fit%rho > 0.0_dp .and. &
      cfar_fit%sigma > 0.0_dp .and. all(shape(cfar_fit%kernel) == [21, 1]) .and. &
      all(ieee_is_finite(cfar_fit%kernel)), 'CFAR spline GLS estimation')
   cfar_tests = nts_cfar_order_tests(cfar_simulation%series, 1, 3)
   call check(cfar_tests%info == 0 .and. &
      cfar_tests%p_value(1) >= 0.0_dp .and. cfar_tests%p_value(1) <= 1.0_dp, &
      'CFAR sequential order test')
   cfar_forecast = nts_cfar_forecast(cfar_fit, cfar_simulation%series, 3)
   call check(cfar_forecast%info == 0 .and. &
      all(shape(cfar_forecast%mean) == [3, 11]) .and. &
      all(ieee_is_finite(cfar_forecast%mean)), 'CFAR recursive forecast')
   cfar_partial = nts_cfar_partial_forecast(cfar_fit, cfar_simulation%series, &
      cfar_simulation%series(90, :3))
   call check(cfar_partial%info == 0 .and. maxval(abs(cfar_partial%mean(1, :3) - &
      cfar_simulation%series(90, :3))) < 1.0e-14_dp, 'CFAR partial forecast')
   do time = 1, 90
      irregular_counts(time) = 6
      irregular_positions(time, :) = [0.0_dp, 0.2_dp, 0.4_dp, 0.6_dp, &
         0.8_dp, 1.0_dp]
      irregular_values(time, :) = cfar_simulation%series(time, 1:11:2)
   end do
   irregular_cfar_fit = nts_cfar_irregular_estimate(irregular_values, &
      irregular_positions, irregular_counts, 1, 3, 10)
   call check(irregular_cfar_fit%info == 0 .and. &
      all(ieee_is_finite(irregular_cfar_fit%kernel)), &
      'irregular CFAR registration and estimation')
   irregular_cfar_tests = nts_cfar_irregular_order_tests(irregular_values, &
      irregular_positions, irregular_counts, 1, 3, 10)
   call check(irregular_cfar_tests%info == 0 .and. &
      irregular_cfar_tests%p_value(1) >= 0.0_dp .and. &
      irregular_cfar_tests%p_value(1) <= 1.0_dp, 'irregular CFAR order test')

   do time = 1, smc_particles
      smc_initial(1, time) = 0.01_dp*real(time - 1 - smc_particles/2, dp)
   end do
   do time = 1, smc_times
      smc_observations(time, 1) = 0.4_dp*sin(0.21_dp*real(time, dp))
      smc_draws(1, :, time) = [(sin(0.37_dp*real(time + particle, dp)), &
         particle=1, smc_particles)]
      smc_uniforms(:, time) = [(modulo(0.61803398875_dp*real(particle + &
         3*time, dp), 1.0_dp), particle=1, smc_particles)]
      smc_schedule(time) = modulo(time, 4) == 0
   end do
   smc_filter = nts_smc_filter_draws(test_smc_step, smc_observations, &
      [0.82_dp, 0.18_dp, 0.30_dp], smc_initial, smc_draws, smc_uniforms, &
      smc_schedule, 2)
   call check(smc_filter%info == 0 .and. &
      all(shape(smc_filter%particles) == [1, smc_particles, smc_times]) .and. &
      all(abs(sum(smc_filter%weights, dim=1) - 1.0_dp) < 1.0e-12_dp) .and. &
      all(smc_filter%effective_sample_size >= 1.0_dp), &
      'generic SMC supplied-draw filtering')
   call check(maxval(abs(smc_filter%delayed_mean(:, :, 0) - &
      smc_filter%filtered_mean)) < 1.0e-14_dp .and. &
      ieee_is_finite(smc_filter%log_likelihood), 'generic SMC delayed estimates')
   smc_smoother = nts_smc_smooth(smc_filter)
   call check(smc_smoother%info == 0 .and. &
      all(shape(smc_smoother%mean) == [1, smc_times]) .and. &
      maxval(abs(smc_smoother%mean(:, smc_times) - &
      smc_filter%filtered_mean(:, smc_times))) < 1.0e-12_dp, &
      'generic SMC genealogical smoothing')
   smc_marginal_smoother = nts_smc_marginal_smooth(smc_filter, &
      [0.82_dp, 0.18_dp], test_smc_transition_density)
   call check(smc_marginal_smoother%info == 0 .and. &
      all(abs(sum(smc_marginal_smoother%weights, dim=1) - 1.0_dp) < &
      1.0e-12_dp) .and. maxval(abs(smc_marginal_smoother%mean(:, smc_times) - &
      smc_filter%filtered_mean(:, smc_times))) < 1.0e-12_dp, &
      'generic SMC marginal backward smoothing')
   call set_random_seed(9021)
   repeated_smc = nts_smc_filter(test_smc_step, smc_observations, &
      [0.82_dp, 0.18_dp, 0.30_dp], smc_initial, 1, smc_schedule, 1)
   call set_random_seed(9021)
   smc_filter = nts_smc_filter(test_smc_step, smc_observations, &
      [0.82_dp, 0.18_dp, 0.30_dp], smc_initial, 1, smc_schedule, 1)
   call check(repeated_smc%info == 0 .and. smc_filter%info == 0 .and. &
      maxval(abs(repeated_smc%particles - smc_filter%particles)) < 1.0e-14_dp, &
      'generic SMC shared RNG reproducibility')
   rb_initial_mean = 0.0_dp
   rb_initial_covariance = 1.0_dp
   rb_smc_filter = nts_rb_smc_filter_draws(test_rb_smc_step, smc_observations, &
      [0.75_dp, 0.12_dp, 0.08_dp, 0.25_dp], smc_initial, rb_initial_mean, &
      rb_initial_covariance, smc_draws, smc_uniforms, smc_schedule)
   call check(rb_smc_filter%info == 0 .and. &
      all(shape(rb_smc_filter%gaussian_mean) == [1, smc_particles, smc_times]) .and. &
      all(rb_smc_filter%gaussian_covariance > 0.0_dp) .and. &
      all(ieee_is_finite(rb_smc_filter%filtered_gaussian_mean)), &
      'Rao-Blackwellized generic SMC filtering')

   print '(a)', 'NTS tests passed'

contains

   pure subroutine test_smc_step(previous, observation, parameters, &
      proposal_draws, proposed, log_weight_increment)
      !! Propagate and weight a scalar Gaussian state for generic SMC tests.
      real(dp), intent(in) :: previous(:, :) !! Previous scalar particles.
      real(dp), intent(in) :: observation(:) !! Current scalar observation.
      real(dp), intent(in) :: parameters(:) !! AR, state-scale, and observation-scale values.
      real(dp), intent(in) :: proposal_draws(:, :) !! Standard-normal proposal draws.
      real(dp), intent(out) :: proposed(:, :) !! Propagated scalar particles.
      real(dp), intent(out) :: log_weight_increment(:) !! Gaussian observation log weights.

      proposed(1, :) = parameters(1)*previous(1, :) + &
         parameters(2)*proposal_draws(1, :)
      log_weight_increment = -0.5_dp*((observation(1) - proposed(1, :))/ &
         parameters(3))**2 - log(parameters(3)*sqrt(2.0_dp*acos(-1.0_dp)))
   end subroutine test_smc_step

   pure subroutine test_rb_smc_step(previous, previous_mean, previous_covariance, &
      observation, parameters, proposal_draws, proposed, proposed_mean, &
      proposed_covariance, log_weight_increment)
      !! Propagate a particle and update a conditional scalar Gaussian state.
      real(dp), intent(in) :: previous(:, :) !! Previous nonlinear particles.
      real(dp), intent(in) :: previous_mean(:, :) !! Previous Gaussian means.
      real(dp), intent(in) :: previous_covariance(:, :, :) !! Previous covariances.
      real(dp), intent(in) :: observation(:) !! Current scalar observation.
      real(dp), intent(in) :: parameters(:) !! AR, particle scale, state variance, and noise variance.
      real(dp), intent(in) :: proposal_draws(:, :) !! Standard-normal proposal draws.
      real(dp), intent(out) :: proposed(:, :) !! Propagated nonlinear particles.
      real(dp), intent(out) :: proposed_mean(:, :) !! Updated conditional means.
      real(dp), intent(out) :: proposed_covariance(:, :, :) !! Updated covariances.
      real(dp), intent(out) :: log_weight_increment(:) !! Predictive observation log weights.
      real(dp) :: predicted_variance, innovation_variance, gain
      integer :: particle

      do particle = 1, size(previous, 2)
         proposed(1, particle) = parameters(1)*previous(1, particle) + &
            parameters(2)*proposal_draws(1, particle)
         predicted_variance = previous_covariance(1, 1, particle) + parameters(3)
         innovation_variance = predicted_variance + parameters(4)
         gain = predicted_variance/innovation_variance
         proposed_mean(1, particle) = previous_mean(1, particle) + gain*( &
            observation(1) - proposed(1, particle) - previous_mean(1, particle))
         proposed_covariance(1, 1, particle) = &
            (1.0_dp - gain)*predicted_variance
         log_weight_increment(particle) = -0.5_dp*(log(2.0_dp*acos(-1.0_dp)* &
            innovation_variance) + (observation(1) - proposed(1, particle) - &
            previous_mean(1, particle))**2/innovation_variance)
      end do
   end subroutine test_rb_smc_step

   pure subroutine test_smc_transition_density(previous, next, parameters, &
      log_density)
      !! Evaluate scalar Gaussian transition densities for SMC smoothing tests.
      real(dp), intent(in) :: previous(:, :) !! Previous scalar particles.
      real(dp), intent(in) :: next(:, :) !! Next scalar particles.
      real(dp), intent(in) :: parameters(:) !! AR coefficient and state scale.
      real(dp), intent(out) :: log_density(:, :) !! Previous-by-next log densities.
      integer :: previous_particle, next_particle

      do next_particle = 1, size(next, 2)
         do previous_particle = 1, size(previous, 2)
            log_density(previous_particle, next_particle) = -0.5_dp*( &
               (next(1, next_particle) - parameters(1)* &
               previous(1, previous_particle))/parameters(2))**2 - &
               log(parameters(2)*sqrt(2.0_dp*acos(-1.0_dp)))
         end do
      end do
   end subroutine test_smc_transition_density

   subroutine check(condition, message)
      !! Stop the test program when a condition fails.
      logical, intent(in) :: condition !! Condition expected to be true.
      character(len=*), intent(in) :: message !! Failure message.

      if (.not. condition) then
         print '(a)', 'FAILED: '//trim(message)
         error stop 1
      end if
   end subroutine check

end program test_nts
