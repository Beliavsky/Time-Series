! SPDX-License-Identifier: GPL-2.0-or-later
! SPDX-FileComment: Regression tests for the bssm translation.
program test_bssm
   use kind_mod, only: dp
   use bssm_mod
   use kfas_mod, only: ssm_model_t
   use time_series_random_mod, only: set_random_seed
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value
   use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan
   implicit none

   type(bssm_particle_filter_t) :: filter, random_filter
   type(bssm_particle_filter_t) :: post_filters(2)
   type(bssm_particle_smoother_t) :: smoother
   type(bssm_simulation_smoother_t) :: simulation_smoother
   type(bssm_simulation_smoother_t) :: random_simulation_smoother
   type(bssm_prediction_t) :: prediction, random_prediction, replicates
   type(bssm_prediction_summary_t) :: fitted_summary, random_fitted_summary
   type(ssm_model_t) :: simulation_model
   type(bssm_ekf_t) :: ekf_fit, iekf_fit
   type(bssm_ekf_t) :: ukf_fit
   type(bssm_multivariate_ekf_t) :: multivariate_ekf_fit
   type(bssm_ekf_smoother_t) :: ekf_smoother, fast_ekf_smoother
   type(bssm_gaussian_approximation_t) :: approximation
   type(bssm_importance_sample_t) :: importance, random_importance
   type(bssm_nonlinear_approximation_t) :: nonlinear_approximation
   type(bssm_multivariate_nonlinear_approximation_t) :: &
      multivariate_nonlinear_approximation
   type(bssm_multivariate_approximation_t) :: multivariate_approximation
   type(bssm_mcmc_t) :: chain, random_chain
   type(bssm_mcmc_t) :: nonlinear_chain, nonlinear_psi_chain
   type(bssm_mcmc_t) :: nonlinear_ekpf_chain
   type(bssm_mcmc_t) :: approximate_chain, nonlinear_approximate_chain
   type(bssm_mcmc_t) :: nonlinear_ekf_chain
   type(bssm_da_mcmc_t) :: da_chain, random_da_chain
   type(bssm_da_mcmc_t) :: nonlinear_da_chain, nonlinear_ekpf_da_chain
   type(bssm_post_correction_t) :: correction_is1, correction_is2
   type(bssm_post_correction_t) :: correction_is3
   type(bssm_sde_state_sample_t) :: sde_states, random_sde_states
   type(bssm_state_posterior_t) :: corrected_states, corrected_signals
   type(bssm_trajectory_sample_t) :: corrected_trajectories
   type(bssm_particle_count_t) :: particle_suggestion
   type(bssm_mcmc_diagnostics_t) :: sokal_diagnostics, geyer_diagnostics
   type(bssm_mcmc_diagnostics_t) :: weighted_diagnostics, chain_diagnostics
   type(bssm_mcmc_diagnostics_t) :: da_chain_diagnostics
   real(dp) :: initial_normals(1, 4), innovations(1, 4, 2)
   real(dp) :: uniforms(4, 2), expected
   real(dp) :: past_normals(2, 2, 2), past_uniforms(2, 2, 2)
   real(dp) :: ekpf_proposals(1, 4, 1), ekpf_uniforms(4, 1)
   real(dp) :: nonlinear_ekpf_proposals(1, 4, 2)
   real(dp) :: multivariate_initial_normals(1, 4)
   real(dp) :: multivariate_innovation_normals(1, 4, 2)
   real(dp) :: multivariate_proposal_normals(1, 4, 2)
   real(dp) :: multivariate_resampling_uniforms(4, 2)
   real(dp) :: observation_schedule(1, 2), transition_schedule(1, 1, 2)
   real(dp) :: noise_schedule(1, 1, 2)
   real(dp) :: state_offset_schedule(1, 2)
   real(dp) :: pseudo_observation, pseudo_variance
   real(dp) :: sde_normals(4, 2, 3)
   real(dp) :: sde_is2_normals(4, 2, 3, 2)
   real(dp) :: sde_is2_uniforms(4, 2, 2), sde_terminal_uniforms(2)
   real(dp) :: pmmh_brownian(4, 1, 2, 3), pmmh_resampling(4, 1, 3)
   real(dp) :: pmmh_proposals(1, 2), pmmh_uniforms(2)
   real(dp) :: da_coarse_brownian(4, 1, 2, 3)
   real(dp) :: da_fine_brownian(4, 2, 2, 3)
   real(dp) :: da_coarse_resampling(4, 1, 3)
   real(dp) :: da_fine_resampling(4, 1, 3)
   real(dp) :: da_proposals(1, 2), da_first_uniforms(2)
   real(dp) :: da_second_uniforms(2)
   real(dp) :: diagnostic_samples(2, 200), diagnostic_weights(200)
   real(dp) :: pmmh_initial_normals(1, 4, 3)
   real(dp) :: pmmh_innovation_normals(1, 4, 2, 3)
   real(dp) :: pmmh_nonlinear_resampling(4, 2, 3)
   real(dp) :: post_initial_normals(1, 4, 2, 3)
   real(dp) :: post_state_normals(1, 4, 2, 2, 3)
   real(dp) :: post_terminal_normals(1, 4, 2, 3)
   real(dp) :: post_resampling_uniforms(4, 2, 2, 3)
   real(dp) :: suggest_initial_normals(1, 4, 2, 2)
   real(dp) :: suggest_state_normals(1, 4, 2, 2, 2)
   real(dp) :: suggest_terminal_normals(1, 4, 2, 2)
   real(dp) :: suggest_resampling_uniforms(4, 2, 2, 2)
   real(dp) :: simulation_initial_normals(1, 2)
   real(dp) :: simulation_observation_normals(1, 1, 2)
   real(dp), allocatable :: simulation_state_normals(:, :, :)
   real(dp), allocatable :: approximation_normals(:, :, :)
   real(dp) :: prediction_state_normals(1, 1, 1)
   real(dp) :: prediction_response_normals(6, 2, 1)
   real(dp) :: prediction_response_uniforms(6, 2, 1)
   real(dp) :: prediction_auxiliary(6, 2)
   real(dp) :: nonlinear_prediction_parameters(3, 2)
   real(dp) :: nonlinear_prediction_state_normals(1, 1, 2)
   real(dp) :: nonlinear_prediction_response_normals(2, 2, 2)
   integer :: approximation_info, diagnostic_index
   integer, allocatable :: ancestors(:)

   prediction_state_normals = 1.0_dp
   prediction_response_normals = 0.0_dp
   prediction_response_normals(1, 1, 1) = 3.0_dp
   prediction_response_normals(6, 1, 1) = 2.0_dp
   prediction_response_uniforms = 0.5_dp
   prediction_response_uniforms(4, :, :) = 0.6_dp
   prediction_auxiliary = 1.0_dp
   prediction_auxiliary(3, :) = 2.0_dp
   prediction = bssm_predictive_draws(reshape([0.0_dp], [1, 1]), &
      reshape([1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], &
      [6, 1, 1]), reshape([1.0_dp], [1, 1, 1]), &
      reshape([1.0_dp], [1, 1, 1]), &
      [bssm_svm, bssm_poisson, bssm_binomial, &
      bssm_negative_binomial, bssm_gamma, bssm_gaussian], &
      [2.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 0.5_dp], &
      prediction_state_normals, prediction_response_normals, &
      prediction_response_uniforms, auxiliary=prediction_auxiliary)
   call check(prediction%info == 0 .and. &
      all(shape(prediction%state) == [1, 2, 1]) .and. &
      maxval(abs(prediction%state(1, :, 1) - [0.0_dp, 1.0_dp])) < &
      1.0e-14_dp .and. &
      maxval(abs(prediction%mean(:, 1, 1) - &
      [0.0_dp, 1.0_dp, 0.5_dp, 1.0_dp, 1.0_dp, 0.0_dp])) < &
      1.0e-14_dp .and. &
      maxval(abs(prediction%response(:, 1, 1) - &
      [6.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, log(2.0_dp), 1.0_dp])) < &
      1.0e-11_dp, 'mixed-family supplied-draw posterior prediction')
   replicates = bssm_predictive_replicates_draws(prediction%state, &
      reshape([1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], &
      [6, 1, 1]), [bssm_svm, bssm_poisson, bssm_binomial, &
      bssm_negative_binomial, bssm_gamma, bssm_gaussian], &
      [2.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 0.5_dp], &
      prediction_response_normals, prediction_response_uniforms, &
      auxiliary=prediction_auxiliary)
   call check(replicates%info == 0 .and. &
      maxval(abs(replicates%state - prediction%state)) < 1.0e-14_dp .and. &
      maxval(abs(replicates%response - prediction%response)) < 1.0e-12_dp, &
      'in-sample posterior predictive replication')
   replicates = bssm_predictive_replicates_draws( &
      reshape([0.0_dp], [1, 1, 1]), reshape([0.0_dp, 0.0_dp], [2, 1, 1]), &
      [bssm_gaussian, bssm_gaussian], [0.0_dp, 0.0_dp], &
      reshape([1.0_dp, 2.0_dp], [2, 1, 1]), &
      reshape([0.5_dp, 0.5_dp], [2, 1, 1]), &
      observation_noise_loading=reshape([1.0_dp, 0.25_dp, 0.5_dp, &
      2.0_dp], [2, 2, 1]))
   call check(replicates%info == 0 .and. &
      maxval(abs(replicates%response(:, 1, 1) - &
      [2.0_dp, 4.25_dp])) < 1.0e-14_dp, &
      'correlated Gaussian posterior predictive response')
   past_normals = 0.0_dp
   past_normals(2, :, 2) = 1.0_dp
   past_uniforms = 0.5_dp
   replicates = bssm_predict_past_draws( &
      reshape([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], [1, 2, 2]), &
      reshape([1.0_dp, 0.5_dp, 0.0_dp, 2.0_dp, -1.0_dp, 1.0_dp], [3, 2]), &
      [bssm_poisson, bssm_gaussian], posterior_observation_model, &
      past_normals, past_uniforms)
   call check(replicates%info == 0 .and. &
      maxval(abs(replicates%signal(:, :, 1) - &
      reshape([1.0_dp, 1.0_dp, 2.0_dp, 2.0_dp], [2, 2]))) < &
      1.0e-14_dp .and. &
      maxval(abs(replicates%signal(:, :, 2) - &
      reshape([7.0_dp, 7.0_dp, 9.0_dp, 9.0_dp], [2, 2]))) < &
      1.0e-14_dp .and. &
      maxval(abs(replicates%mean(1, :, 2) - exp([7.0_dp, 9.0_dp]))) < &
      1.0e-10_dp .and. &
      maxval(abs(replicates%response(2, :, 2) - [8.0_dp, 10.0_dp])) < &
      1.0e-14_dp, 'paired-parameter posterior in-sample prediction')
   call set_random_seed(1867)
   random_prediction = bssm_predict_past( &
      reshape([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], [1, 2, 2]), &
      reshape([1.0_dp, 0.5_dp, 0.0_dp, 2.0_dp, -1.0_dp, 1.0_dp], [3, 2]), &
      [bssm_poisson, bssm_gaussian], posterior_observation_model)
   call check(random_prediction%info == 0 .and. &
      all(shape(random_prediction%response) == [2, 2, 2]) .and. &
      all(ieee_is_finite(random_prediction%response)), &
      'shared-stream posterior in-sample prediction')
   call set_random_seed(1879)
   random_prediction = bssm_predictive_sample( &
      reshape([0.0_dp, 0.5_dp], [1, 2]), &
      reshape([1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], &
      [6, 1, 1]), reshape([1.0_dp], [1, 1, 1]), &
      reshape([0.2_dp], [1, 1, 1]), &
      [bssm_svm, bssm_poisson, bssm_binomial, &
      bssm_negative_binomial, bssm_gamma, bssm_gaussian], &
      [2.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 0.5_dp], 2, &
      auxiliary=prediction_auxiliary)
   call check(random_prediction%info == 0 .and. &
      all(shape(random_prediction%response) == [6, 2, 2]) .and. &
      all(ieee_is_finite(random_prediction%state)) .and. &
      all(ieee_is_finite(random_prediction%response)), &
      'shared-stream mixed-family posterior prediction')
   nonlinear_prediction_parameters = reshape( &
      [2.0_dp, 1.0_dp, 0.5_dp, 1.0_dp, -1.0_dp, 0.25_dp], [3, 2])
   nonlinear_prediction_state_normals = &
      reshape([2.0_dp, -4.0_dp], [1, 1, 2])
   nonlinear_prediction_response_normals = 0.0_dp
   nonlinear_prediction_response_normals(:, 1, 1) = [2.0_dp, 1.0_dp]
   prediction = bssm_nonlinear_predictive_draws( &
      reshape([1.0_dp, 2.0_dp], [1, 2]), &
      nonlinear_prediction_parameters, 2, 1, &
      nonlinear_prediction_observation, nonlinear_prediction_transition, &
      nonlinear_prediction_state_normals, &
      nonlinear_prediction_response_normals)
   call check(prediction%info == 0 .and. &
      all(shape(prediction%state) == [1, 2, 2]) .and. &
      maxval(abs(prediction%state(1, :, 1) - [1.0_dp, 4.0_dp])) < &
      1.0e-14_dp .and. &
      maxval(abs(prediction%state(1, :, 2) - [2.0_dp, 0.0_dp])) < &
      1.0e-14_dp .and. &
      maxval(abs(prediction%mean(:, 1, 1) - [3.0_dp, 3.0_dp])) < &
      1.0e-14_dp .and. &
      maxval(abs(prediction%response(:, 1, 1) - [4.0_dp, 4.5_dp])) < &
      1.0e-14_dp .and. &
      maxval(abs(prediction%mean(:, 2, 1) - [18.0_dp, 9.0_dp])) < &
      1.0e-14_dp, 'supplied-draw nonlinear Gaussian prediction')
   replicates = bssm_nonlinear_predictive_replicates_draws( &
      prediction%state, nonlinear_prediction_parameters, &
      nonlinear_prediction_observation, &
      nonlinear_prediction_response_normals)
   call check(replicates%info == 0 .and. &
      maxval(abs(replicates%response - prediction%response)) < 1.0e-14_dp, &
      'nonlinear Gaussian posterior predictive replication')
   replicates = bssm_nonlinear_predict_past_draws(prediction%state, &
      nonlinear_prediction_parameters, nonlinear_prediction_observation, &
      nonlinear_prediction_response_normals)
   call check(replicates%info == 0 .and. &
      maxval(abs(replicates%response - prediction%response)) < 1.0e-14_dp, &
      'nonlinear posterior in-sample prediction')
   fitted_summary = bssm_summarize_prediction(replicates, &
      [0.25_dp, 0.75_dp], [1.0_dp, 3.0_dp], [2, 1])
   call check(fitted_summary%info == 0 .and. &
      all(shape(fitted_summary%signal_quantile) == [2, 2, 2]) .and. &
      maxval(abs(fitted_summary%signal_mean(:, 1) - &
      [4.2_dp, 3.0_dp])) < 1.0e-14_dp .and. &
      maxval(abs(fitted_summary%observation_mean(:, 1) - &
      [4.6_dp, 3.6_dp])) < 1.0e-14_dp .and. &
      maxval(abs(fitted_summary%signal_quantile(1, 1, :) - &
      [3.0_dp, 5.0_dp])) < 1.0e-14_dp, &
      'weighted nonlinear fitted prediction summary')
   fitted_summary = bssm_nonlinear_fitted_summary_draws(prediction%state, &
      nonlinear_prediction_parameters, nonlinear_prediction_observation, &
      nonlinear_prediction_response_normals, [0.25_dp, 0.75_dp], &
      [1.0_dp, 3.0_dp], [2, 1])
   call check(fitted_summary%info == 0 .and. &
      maxval(abs(fitted_summary%observation_mean(:, 1) - &
      [4.6_dp, 3.6_dp])) < 1.0e-14_dp, &
      'supplied-draw nonlinear fitted summary')
   call set_random_seed(1882)
   random_fitted_summary = bssm_nonlinear_fitted_summary(prediction%state, &
      nonlinear_prediction_parameters, 2, nonlinear_prediction_observation, &
      [0.025_dp, 0.975_dp], [1.0_dp, 1.0_dp])
   call check(random_fitted_summary%info == 0 .and. &
      all(ieee_is_finite(random_fitted_summary%signal_mean)) .and. &
      all(ieee_is_finite(random_fitted_summary%observation_sd)) .and. &
      all(ieee_is_finite(random_fitted_summary%observation_quantile)), &
      'shared-stream nonlinear fitted summary')
   call set_random_seed(1883)
   random_prediction = bssm_nonlinear_predictive_sample( &
      reshape([1.0_dp, 2.0_dp], [1, 2]), &
      nonlinear_prediction_parameters, 2, 1, 3, &
      nonlinear_prediction_observation, nonlinear_prediction_transition)
   call check(random_prediction%info == 0 .and. &
      all(shape(random_prediction%response) == [2, 3, 2]) .and. &
      all(ieee_is_finite(random_prediction%state)) .and. &
      all(ieee_is_finite(random_prediction%response)), &
      'shared-stream nonlinear Gaussian prediction')
   random_prediction = bssm_nonlinear_predictive_replicates( &
      prediction%state, nonlinear_prediction_parameters, 2, &
      nonlinear_prediction_observation)
   call check(random_prediction%info == 0 .and. &
      all(shape(random_prediction%response) == [2, 2, 2]) .and. &
      all(ieee_is_finite(random_prediction%response)), &
      'shared-stream nonlinear predictive replication')

   allocate(simulation_model%y(1, 1), simulation_model%z(1, 1, 1))
   allocate(simulation_model%h(1, 1, 1))
   allocate(simulation_model%transition(1, 1, 1))
   allocate(simulation_model%r(1, 1, 1), simulation_model%q(1, 1, 1))
   allocate(simulation_model%a1(1), simulation_model%p1(1, 1))
   allocate(simulation_model%p1inf(1, 1), simulation_model%missing(1, 1))
   simulation_model%y = 1.0_dp
   simulation_model%z = 1.0_dp
   simulation_model%h = 1.0_dp
   simulation_model%transition = 1.0_dp
   simulation_model%r = 1.0_dp
   simulation_model%q = 1.0_dp
   simulation_model%a1 = 0.0_dp
   simulation_model%p1 = 1.0_dp
   simulation_model%p1inf = 0.0_dp
   simulation_model%missing = .false.
   simulation_initial_normals = reshape([1.0_dp, -1.0_dp], [1, 2])
   simulation_observation_normals = reshape([-1.0_dp, 1.0_dp], [1, 1, 2])
   allocate(simulation_state_normals(1, 0, 2))
   simulation_smoother = bssm_simulation_smoother_draws(simulation_model, &
      simulation_initial_normals, simulation_observation_normals, &
      simulation_state_normals)
   call check(simulation_smoother%info == 0 .and. &
      maxval(abs(simulation_smoother%trajectories(1, 1, :) - &
      [1.5_dp, -0.5_dp])) < 1.0e-13_dp .and. &
      abs(simulation_smoother%mean(1, 1) - 0.5_dp) < 1.0e-13_dp .and. &
      abs(simulation_smoother%covariance(1, 1, 1) - 1.0_dp) < &
      1.0e-13_dp .and. &
      maxval(abs(simulation_smoother%observations - 1.0_dp)) < &
      1.0e-13_dp, 'Durbin-Koopman supplied-draw simulation smoother')
   call set_random_seed(1889)
   random_simulation_smoother = bssm_simulation_smoother(simulation_model, &
      4, antithetic=.true.)
   call check(random_simulation_smoother%info == 0 .and. &
      abs(random_simulation_smoother%mean(1, 1) - 0.5_dp) < 1.0e-13_dp .and. &
      all(ieee_is_finite(random_simulation_smoother%trajectories)), &
      'antithetic Gaussian simulation smoother')
   simulation_model%y = ieee_value(0.0_dp, ieee_quiet_nan)
   simulation_model%missing = .true.
   simulation_smoother = bssm_simulation_smoother_draws(simulation_model, &
      simulation_initial_normals, simulation_observation_normals, &
      simulation_state_normals)
   call check(simulation_smoother%info == 0 .and. &
      maxval(abs(simulation_smoother%trajectories(1, 1, :) - &
      [1.0_dp, -1.0_dp])) < 1.0e-13_dp .and. &
      maxval(abs(simulation_smoother%observations)) < 1.0e-13_dp, &
      'simulation-smoother missing-observation completion')

   expected = 2.0_dp*log(3.0_dp) - 3.0_dp - log(2.0_dp)
   call check(abs(bssm_observation_log_density(2.0_dp, 0.0_dp, &
      bssm_poisson, 1.0_dp, 3.0_dp) - expected) < 1.0e-13_dp, &
      'Poisson observation density')
   expected = log(10.0_dp) - 5.0_dp*log(2.0_dp)
   call check(abs(bssm_observation_log_density(2.0_dp, 0.0_dp, &
      bssm_binomial, 1.0_dp, 5.0_dp) - expected) < 1.0e-13_dp, &
      'binomial observation density')
   expected = log_gamma(6.0_dp) - log_gamma(4.0_dp) - log_gamma(3.0_dp) + &
      4.0_dp*log(4.0_dp) + 2.0_dp*log(3.0_dp) - 6.0_dp*log(7.0_dp)
   call check(abs(bssm_observation_log_density(2.0_dp, 0.0_dp, &
      bssm_negative_binomial, 4.0_dp, 3.0_dp) - expected) < 1.0e-13_dp, &
      'negative-binomial observation density')
   expected = 3.0_dp*log(3.0_dp) - log_gamma(3.0_dp) + &
      2.0_dp*log(2.0_dp) - 3.0_dp*log(2.0_dp) - 3.0_dp
   call check(abs(bssm_observation_log_density(2.0_dp, 0.0_dp, &
      bssm_gamma, 3.0_dp, 2.0_dp) - expected) < 1.0e-13_dp, &
      'Gamma observation density')
   expected = -0.5_dp*(log(8.0_dp*acos(-1.0_dp)) + 0.25_dp)
   call check(abs(bssm_observation_log_density(1.0_dp, 0.0_dp, &
      bssm_gaussian, 2.0_dp, 1.0_dp) - expected) < 1.0e-13_dp, &
      'Gaussian observation density')

   ancestors = bssm_stratified_resample([0.1_dp, 0.2_dp, 0.7_dp], &
      [0.5_dp, 0.5_dp, 0.5_dp])
   call check(all(ancestors == [2, 3, 3]), 'stratified resampling')

   call check(abs(bssm_sde_euler_step(0.1_dp, 1.0_dp, -1.0_dp, &
      [0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp], sde_drift, sde_diffusion, &
      .true.) - 0.9_dp) < 1.0e-14_dp, 'positive Euler-Maruyama step')
   call check(abs(bssm_sde_milstein_step(2.0_dp, 0.25_dp, 1.5_dp, &
      [0.0_dp, 0.0_dp, 0.3_dp, 1.0_dp], sde_drift, &
      sde_multiplicative_diffusion, sde_diffusion_derivative, .false.) - &
      2.478125_dp) < 1.0e-14_dp, 'Milstein derivative correction')
   sde_normals = 0.0_dp
   uniforms = 0.5_dp
   filter = bssm_sde_bootstrap_filter_draws([1.4375_dp, 1.68359375_dp], &
      1.0_dp, [0.5_dp, 2.0_dp, 0.0_dp, 1.0_dp], 2, &
      sde_observation_log_density, sde_drift, sde_diffusion, .false., &
      sde_normals, uniforms, sde_zero_diffusion_derivative)
   call check(filter%info == 0 .and. &
      maxval(abs(filter%predicted_mean(1, :) - [1.4375_dp, &
      1.68359375_dp, 1.822021484375_dp])) < 1.0e-14_dp .and. &
      maxval(abs(filter%weights - 0.25_dp)) < 1.0e-14_dp .and. &
      abs(filter%log_likelihood + log(2.0_dp*acos(-1.0_dp))) < 1.0e-13_dp, &
      'draw-driven SDE bootstrap filter')
   smoother = bssm_particle_smoother(filter)
   call check(smoother%info == 0 .and. &
      all(shape(smoother%trajectories) == [1, 2, 4]) .and. &
      maxval(abs(smoother%mean(1, :) - [1.4375_dp, 1.68359375_dp])) < &
      1.0e-14_dp, 'SDE genealogical particle smoother')
   smoother = bssm_psi_particle_smoother(filter)
   call check(smoother%info == 0 .and. &
      all(shape(smoother%trajectories) == [1, 3, 4]) .and. &
      maxval(abs(smoother%mean(1, :) - &
      [1.4375_dp, 1.68359375_dp, 1.822021484375_dp])) < 1.0e-14_dp, &
      'terminal psi particle smoother')
   post_filters(1) = filter
   post_filters(2) = filter
   corrected_states = bssm_post_corrected_particle_moments(post_filters, &
      [0.25_dp, 0.75_dp])
   call check(corrected_states%info == 0 .and. &
      all(shape(corrected_states%mean) == [1, 3]) .and. &
      maxval(abs(corrected_states%mean(1, :) - &
      [1.4375_dp, 1.68359375_dp, 1.822021484375_dp])) < 1.0e-14_dp, &
      'post-corrected terminal particle moments')
   corrected_trajectories = &
      bssm_post_corrected_particle_trajectories_draws(post_filters, &
      [0.25_dp, 0.75_dp], reshape([0.0_dp, 0.0_dp, 0.999_dp, 0.999_dp], &
      [2, 2]))
   call check(corrected_trajectories%info == 0 .and. &
      all(shape(corrected_trajectories%trajectories) == [1, 3, 2]) .and. &
      all(corrected_trajectories%source_sample == [1, 2]) .and. &
      all(corrected_trajectories%source_particle == [1, 4]), &
      'draw-driven post-corrected terminal trajectories')
   call set_random_seed(1701)
   corrected_trajectories = bssm_post_corrected_particle_trajectories( &
      post_filters, [0.25_dp, 0.75_dp], 3)
   call check(corrected_trajectories%info == 0 .and. &
      all(shape(corrected_trajectories%trajectories) == [1, 3, 3]) .and. &
      all(ieee_is_finite(corrected_trajectories%trajectories)), &
      'shared-stream post-corrected terminal trajectories')
   sde_is2_normals = 0.0_dp
   sde_is2_uniforms = 0.5_dp
   sde_terminal_uniforms = [0.0_dp, 0.99_dp]
   sde_states = bssm_sde_is2_state_sampler_draws( &
      [1.4375_dp, 1.68359375_dp], 1.0_dp, &
      reshape([0.5_dp, 2.0_dp, 0.0_dp, 1.0_dp, &
      0.5_dp, 2.0_dp, 0.0_dp, 1.0_dp], [4, 2]), 2, &
      [-log(2.0_dp*acos(-1.0_dp)), &
      -log(2.0_dp*acos(-1.0_dp)) - log(2.0_dp)], &
      sde_observation_log_density, sde_drift, sde_diffusion, .false., &
      sde_is2_normals, sde_is2_uniforms, sde_terminal_uniforms, &
      sde_zero_diffusion_derivative)
   call check(sde_states%info == 0 .and. &
      all(shape(sde_states%trajectories) == [1, 3, 2]) .and. &
      maxval(abs(sde_states%trajectories(1, :, 1) - &
      [1.4375_dp, 1.68359375_dp, 1.822021484375_dp])) < 1.0e-14_dp .and. &
      maxval(abs(sde_states%fine_log_likelihood + &
      log(2.0_dp*acos(-1.0_dp)))) < 1.0e-13_dp .and. &
      maxval(abs(sde_states%weight - [1.0_dp/3.0_dp, 2.0_dp/3.0_dp])) < &
      1.0e-13_dp .and. &
      abs(sde_states%log_mean_weight - log(1.5_dp)) < 1.0e-13_dp .and. &
      abs(sde_states%effective_sample_size - 1.8_dp) < 1.0e-13_dp, &
      'SDE IS2 supplied-draw state sampler')
   call set_random_seed(1897)
   random_sde_states = bssm_sde_is2_state_sampler( &
      [1.4375_dp, 1.68359375_dp], 1.0_dp, &
      reshape([0.5_dp, 2.0_dp, 0.0_dp, 1.0_dp, &
      0.5_dp, 2.0_dp, 0.0_dp, 1.0_dp], [4, 2]), 2, &
      [0.0_dp, 0.0_dp], sde_observation_log_density, sde_drift, &
      sde_diffusion, .false., 4, sde_zero_diffusion_derivative)
   call check(random_sde_states%info == 0 .and. &
      all(shape(random_sde_states%trajectories) == [1, 3, 2]) .and. &
      all(ieee_is_finite(random_sde_states%trajectories)) .and. &
      abs(sum(random_sde_states%weight) - 1.0_dp) < 1.0e-14_dp, &
      'SDE IS2 shared-stream state sampler')
   call set_random_seed(1901)
   random_filter = bssm_sde_bootstrap_filter([1.4_dp, 1.7_dp], 1.0_dp, &
      [0.5_dp, 2.0_dp, 0.2_dp, 1.0_dp], 4, &
      sde_observation_log_density, sde_drift, sde_diffusion, .false., 30)
   call check(random_filter%info == 0 .and. &
      ieee_is_finite(random_filter%log_likelihood) .and. &
      all(ieee_is_finite(random_filter%particles)), &
      'random-stream SDE bootstrap filter')
   pmmh_brownian = 0.0_dp
   pmmh_resampling = 0.5_dp
   pmmh_proposals = reshape([1.0_dp, 0.2_dp], [1, 2])
   pmmh_uniforms = 0.5_dp
   chain = bssm_sde_pmmh_draws([0.0_dp], 0.0_dp, [0.0_dp], 1, 4, &
      pmmh_observation_log_density, pmmh_drift, pmmh_diffusion, .false., &
      pmmh_prior, reshape([1.0_dp], [1, 1]), pmmh_proposals, &
      pmmh_uniforms, pmmh_brownian, pmmh_resampling)
   expected = -0.5_dp*log(2.0_dp*acos(-1.0_dp))
   call check(chain%info == 0 .and. &
      all(shape(chain%parameters) == [1, 3]) .and. &
      maxval(abs(chain%parameters(1, :) - [0.0_dp, 0.0_dp, 0.2_dp])) < &
      1.0e-14_dp .and. all(chain%accepted .eqv. [.false., .true.]) .and. &
      maxval(abs(chain%acceptance_probability - &
      [exp(-1.0_dp), exp(-0.04_dp)])) < 1.0e-14_dp .and. &
      maxval(abs(chain%log_likelihood - &
      [expected, expected, expected - 0.02_dp])) < 1.0e-13_dp .and. &
      abs(chain%acceptance_rate - 0.5_dp) < 1.0e-14_dp, &
      'draw-driven SDE particle marginal Metropolis-Hastings')
   call set_random_seed(1903)
   random_chain = bssm_sde_pmmh([0.0_dp], 0.0_dp, [0.0_dp], 1, 10, 5, &
      pmmh_observation_log_density, pmmh_drift, pmmh_diffusion, .false., &
      pmmh_prior, reshape([0.1_dp], [1, 1]), &
      target_acceptance=0.234_dp, adaptation_exponent=0.6_dp)
   call check(random_chain%info == 0 .and. &
      all(ieee_is_finite(random_chain%parameters)) .and. &
      all(ieee_is_finite(random_chain%log_likelihood)) .and. &
      all(shape(random_chain%final_proposal_factor) == [1, 1]) .and. &
      all(ieee_is_finite(random_chain%final_proposal_factor)) .and. &
      random_chain%final_proposal_factor(1, 1) > 0.0_dp .and. &
      random_chain%acceptance_rate >= 0.0_dp .and. &
      random_chain%acceptance_rate <= 1.0_dp, &
      'adaptive random-stream SDE particle marginal Metropolis-Hastings')
   da_coarse_brownian = 0.0_dp
   da_fine_brownian = 0.0_dp
   da_coarse_resampling = 0.5_dp
   da_fine_resampling = 0.5_dp
   da_proposals = reshape([-1.0_dp, 2.0_dp], [1, 2])
   da_first_uniforms = 0.5_dp
   da_second_uniforms = 0.5_dp
   da_chain = bssm_sde_da_pmmh_draws([0.0_dp], 1.0_dp, [0.0_dp], &
      1, 2, 4, pmmh_observation_log_density, da_pmmh_drift, &
      pmmh_diffusion, .false., pmmh_prior, reshape([1.0_dp], [1, 1]), &
      da_proposals, da_first_uniforms, da_second_uniforms, &
      da_coarse_brownian, da_coarse_resampling, da_fine_brownian, &
      da_fine_resampling)
   expected = -0.5_dp*log(2.0_dp*acos(-1.0_dp))
   call check(da_chain%info == 0 .and. &
      maxval(abs(da_chain%parameters(1, :) - &
      [0.0_dp, -1.0_dp, -1.0_dp])) < 1.0e-14_dp .and. &
      maxval(abs(da_chain%coarse_log_likelihood - &
      [expected - 0.5_dp, expected, expected])) < 1.0e-13_dp .and. &
      maxval(abs(da_chain%fine_log_likelihood - &
      [expected - 0.5_dp, expected - 0.03125_dp, &
      expected - 0.03125_dp])) < 1.0e-13_dp .and. &
      maxval(abs(da_chain%first_stage_probability - &
      [1.0_dp, exp(-2.0_dp)])) < 1.0e-13_dp .and. &
      maxval(abs(da_chain%second_stage_probability - &
      [exp(-0.03125_dp), 0.0_dp])) < 1.0e-13_dp .and. &
      all(da_chain%first_stage_accepted .eqv. [.true., .false.]) .and. &
      all(da_chain%accepted .eqv. [.true., .false.]) .and. &
      abs(da_chain%first_stage_acceptance_rate - 0.5_dp) < 1.0e-14_dp .and. &
      abs(da_chain%acceptance_rate - 0.5_dp) < 1.0e-14_dp, &
      'draw-driven delayed-acceptance SDE PMMH')
   call set_random_seed(1905)
   random_da_chain = bssm_sde_da_pmmh([0.0_dp], 1.0_dp, [0.0_dp], &
      1, 2, 10, 5, pmmh_observation_log_density, da_pmmh_drift, &
      pmmh_diffusion, .false., pmmh_prior, reshape([0.1_dp], [1, 1]), &
      target_acceptance=0.234_dp, adaptation_exponent=0.6_dp)
   call check(random_da_chain%info == 0 .and. &
      all(ieee_is_finite(random_da_chain%parameters)) .and. &
      all(ieee_is_finite(random_da_chain%coarse_log_likelihood)) .and. &
      all(ieee_is_finite(random_da_chain%fine_log_likelihood)) .and. &
      all(ieee_is_finite(random_da_chain%final_proposal_factor)) .and. &
      random_da_chain%first_stage_acceptance_rate >= &
      random_da_chain%acceptance_rate .and. &
      random_da_chain%acceptance_rate >= 0.0_dp .and. &
      random_da_chain%first_stage_acceptance_rate <= 1.0_dp, &
      'adaptive random-stream delayed-acceptance SDE PMMH')
   call set_random_seed(1906)
   nonlinear_da_chain = bssm_nonlinear_da_pmmh([0.0_dp, 0.0_dp], &
      [0.0_dp], reshape([0.1_dp], [1, 1]), [0.0_dp], 1, &
      pmmh_psi_observation, pmmh_psi_transition, &
      pmmh_observation_log_density, pmmh_ekpf_transition, pmmh_prior, &
      20, 1.0e-10_dp, reshape([0.1_dp], [1, 1]), 10, 5, &
      target_acceptance=0.234_dp, adaptation_exponent=0.6_dp)
   call check(nonlinear_da_chain%info == 0 .and. &
      all(ieee_is_finite(nonlinear_da_chain%coarse_log_likelihood)) .and. &
      all(ieee_is_finite(nonlinear_da_chain%fine_log_likelihood)) .and. &
      nonlinear_da_chain%first_stage_acceptance_rate >= &
      nonlinear_da_chain%acceptance_rate, &
      'Gaussian-screened nonlinear bootstrap PMMH')
   call set_random_seed(1908)
   nonlinear_ekpf_da_chain = bssm_nonlinear_ekpf_da_pmmh( &
      [0.0_dp, 0.0_dp], [0.0_dp], reshape([0.1_dp], [1, 1]), &
      [0.0_dp], 1, pmmh_psi_observation, pmmh_psi_transition, &
      pmmh_ekpf_transition, pmmh_prior, 1, 1.0e-10_dp, &
      reshape([0.1_dp], [1, 1]), 10, 5, target_acceptance=0.234_dp, &
      adaptation_exponent=0.6_dp)
   call check(nonlinear_ekpf_da_chain%info == 0 .and. &
      all(ieee_is_finite(nonlinear_ekpf_da_chain%coarse_log_likelihood)) .and. &
      all(ieee_is_finite(nonlinear_ekpf_da_chain%fine_log_likelihood)) .and. &
      nonlinear_ekpf_da_chain%first_stage_acceptance_rate >= &
      nonlinear_ekpf_da_chain%acceptance_rate, &
      'EKF-screened nonlinear EKF-proposal PMMH')
   pmmh_initial_normals = 0.0_dp
   pmmh_innovation_normals = 0.0_dp
   pmmh_nonlinear_resampling = 0.5_dp
   pmmh_proposals = reshape([1.0_dp, 0.2_dp], [1, 2])
   pmmh_uniforms = 0.5_dp
   nonlinear_chain = bssm_nonlinear_pmmh_draws([0.0_dp, 0.0_dp], &
      [0.0_dp], reshape([0.0_dp], [1, 1]), [0.0_dp], 1, &
      pmmh_observation_log_density, pmmh_nonlinear_transition, pmmh_prior, &
      reshape([1.0_dp], [1, 1]), pmmh_proposals, pmmh_uniforms, &
      pmmh_initial_normals, pmmh_innovation_normals, &
      pmmh_nonlinear_resampling)
   call check(nonlinear_chain%info == 0 .and. &
      maxval(abs(nonlinear_chain%parameters(1, :) - &
      [0.0_dp, 0.0_dp, 0.2_dp])) < 1.0e-14_dp .and. &
      maxval(abs(nonlinear_chain%acceptance_probability - &
      [exp(-1.0_dp), exp(-0.04_dp)])) < 1.0e-13_dp .and. &
      all(nonlinear_chain%accepted .eqv. [.false., .true.]), &
      'draw-driven nonlinear bootstrap-filter PMMH')
   nonlinear_chain = bssm_multivariate_nonlinear_pmmh_draws( &
      reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], [2, 2]), [0.0_dp], &
      reshape([0.0_dp], [1, 1]), [0.0_dp], 2, 1, &
      multivariate_linear_gaussian_observation, pmmh_nonlinear_transition, &
      pmmh_prior, reshape([1.0_dp], [1, 1]), pmmh_proposals, &
      pmmh_uniforms, pmmh_initial_normals, pmmh_innovation_normals, &
      pmmh_nonlinear_resampling)
   call check(nonlinear_chain%info == 0 .and. &
      maxval(abs(nonlinear_chain%parameters(1, :) - &
      [0.0_dp, 0.0_dp, 0.2_dp])) < 1.0e-14_dp .and. &
      maxval(abs(nonlinear_chain%acceptance_probability - &
      [exp(-1.5_dp), exp(-0.06_dp)])) < 1.0e-13_dp .and. &
      all(nonlinear_chain%accepted .eqv. [.false., .true.]), &
      'draw-driven multivariate nonlinear bootstrap PMMH')
   call set_random_seed(1911)
   nonlinear_chain = bssm_multivariate_nonlinear_pmmh( &
      reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], [2, 2]), [0.0_dp], &
      reshape([0.1_dp], [1, 1]), [0.0_dp], 2, 1, &
      multivariate_linear_gaussian_observation, pmmh_nonlinear_transition, &
      pmmh_prior, reshape([0.1_dp], [1, 1]), 10, 5, &
      target_acceptance=0.234_dp, adaptation_exponent=0.6_dp)
   call check(nonlinear_chain%info == 0 .and. &
      all(ieee_is_finite(nonlinear_chain%parameters)) .and. &
      all(ieee_is_finite(nonlinear_chain%log_likelihood)) .and. &
      all(ieee_is_finite(nonlinear_chain%final_proposal_factor)) .and. &
      nonlinear_chain%acceptance_rate >= 0.0_dp .and. &
      nonlinear_chain%acceptance_rate <= 1.0_dp, &
      'adaptive random-stream multivariate nonlinear bootstrap PMMH')
   nonlinear_psi_chain = bssm_multivariate_nonlinear_psi_pmmh_draws( &
      reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], [2, 2]), [0.0_dp], &
      reshape([0.1_dp], [1, 1]), [0.0_dp], 2, 1, &
      multivariate_linear_gaussian_observation, pmmh_psi_transition, &
      pmmh_prior, 20, 1.0e-12_dp, reshape([1.0_dp], [1, 1]), &
      pmmh_proposals, pmmh_uniforms, pmmh_initial_normals, &
      pmmh_innovation_normals, pmmh_initial_normals, &
      pmmh_nonlinear_resampling)
   call check(nonlinear_psi_chain%info == 0 .and. &
      all(shape(nonlinear_psi_chain%parameters) == [1, 3]) .and. &
      all(ieee_is_finite(nonlinear_psi_chain%parameters)) .and. &
      all(ieee_is_finite(nonlinear_psi_chain%log_likelihood)) .and. &
      all(ieee_is_finite(nonlinear_psi_chain%acceptance_probability)), &
      'draw-driven multivariate nonlinear psi-filter PMMH')
   call set_random_seed(1912)
   nonlinear_psi_chain = bssm_multivariate_nonlinear_psi_pmmh( &
      reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], [2, 2]), [0.0_dp], &
      reshape([0.1_dp], [1, 1]), [0.0_dp], 2, 1, &
      multivariate_linear_gaussian_observation, pmmh_psi_transition, &
      pmmh_prior, 20, 1.0e-10_dp, reshape([0.1_dp], [1, 1]), 10, 5, &
      target_acceptance=0.234_dp, adaptation_exponent=0.6_dp)
   call check(nonlinear_psi_chain%info == 0 .and. &
      all(ieee_is_finite(nonlinear_psi_chain%parameters)) .and. &
      all(ieee_is_finite(nonlinear_psi_chain%log_likelihood)) .and. &
      all(ieee_is_finite(nonlinear_psi_chain%final_proposal_factor)) .and. &
      nonlinear_psi_chain%acceptance_rate >= 0.0_dp .and. &
      nonlinear_psi_chain%acceptance_rate <= 1.0_dp, &
      'adaptive random-stream multivariate nonlinear psi-filter PMMH')
   nonlinear_da_chain = bssm_multivariate_nonlinear_psi_da_pmmh_draws( &
      reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], [2, 2]), [0.0_dp], &
      reshape([0.1_dp], [1, 1]), [0.0_dp], 2, 1, &
      multivariate_linear_gaussian_observation, pmmh_psi_transition, &
      pmmh_prior, 20, 1.0e-12_dp, reshape([1.0_dp], [1, 1]), &
      pmmh_proposals, pmmh_uniforms, pmmh_uniforms, pmmh_initial_normals, &
      pmmh_innovation_normals, pmmh_initial_normals, &
      pmmh_nonlinear_resampling)
   call check(nonlinear_da_chain%info == 0 .and. &
      all(ieee_is_finite(nonlinear_da_chain%coarse_log_likelihood)) .and. &
      all(ieee_is_finite(nonlinear_da_chain%fine_log_likelihood)) .and. &
      maxval(abs(nonlinear_da_chain%coarse_log_likelihood - &
      nonlinear_da_chain%fine_log_likelihood)) < 1.0e-12_dp .and. &
      nonlinear_da_chain%first_stage_acceptance_rate >= &
      nonlinear_da_chain%acceptance_rate, &
      'draw-driven Gaussian-screened multivariate psi PMMH')
   call set_random_seed(1914)
   nonlinear_da_chain = bssm_multivariate_nonlinear_psi_da_pmmh( &
      reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], [2, 2]), [0.0_dp], &
      reshape([0.1_dp], [1, 1]), [0.0_dp], 2, 1, &
      multivariate_linear_gaussian_observation, pmmh_psi_transition, &
      pmmh_prior, 20, 1.0e-10_dp, reshape([0.1_dp], [1, 1]), 10, 5, &
      target_acceptance=0.234_dp, adaptation_exponent=0.6_dp)
   call check(nonlinear_da_chain%info == 0 .and. &
      all(ieee_is_finite(nonlinear_da_chain%coarse_log_likelihood)) .and. &
      all(ieee_is_finite(nonlinear_da_chain%fine_log_likelihood)) .and. &
      all(ieee_is_finite(nonlinear_da_chain%final_proposal_factor)) .and. &
      nonlinear_da_chain%first_stage_acceptance_rate >= &
      nonlinear_da_chain%acceptance_rate .and. &
      nonlinear_da_chain%acceptance_rate >= 0.0_dp, &
      'adaptive Gaussian-screened multivariate psi PMMH')
   nonlinear_ekpf_chain = bssm_multivariate_ekpf_pmmh_draws( &
      reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], [2, 2]), [0.0_dp], &
      reshape([0.1_dp], [1, 1]), [0.0_dp], 2, 1, &
      multivariate_linear_gaussian_observation, pmmh_ekpf_transition, &
      pmmh_prior, 1, 1.0e-12_dp, reshape([1.0_dp], [1, 1]), &
      pmmh_proposals, pmmh_uniforms, pmmh_initial_normals, &
      pmmh_innovation_normals, pmmh_nonlinear_resampling)
   expected = 0.5_dp + 1.0_dp/1.02_dp
   call check(nonlinear_ekpf_chain%info == 0 .and. &
      maxval(abs(nonlinear_ekpf_chain%parameters(1, :) - &
      [0.0_dp, 0.0_dp, 0.2_dp])) < 1.0e-14_dp .and. &
      maxval(abs(nonlinear_ekpf_chain%acceptance_probability - &
      [exp(-expected), exp(-0.04_dp*expected)])) < 1.0e-12_dp .and. &
      all(nonlinear_ekpf_chain%accepted .eqv. [.false., .true.]), &
      'draw-driven multivariate nonlinear EKPF-PMMH')
   call set_random_seed(1913)
   nonlinear_ekpf_chain = bssm_multivariate_ekpf_pmmh( &
      reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], [2, 2]), [0.0_dp], &
      reshape([0.1_dp], [1, 1]), [0.0_dp], 2, 1, &
      multivariate_linear_gaussian_observation, pmmh_ekpf_transition, &
      pmmh_prior, 2, 1.0e-10_dp, reshape([0.1_dp], [1, 1]), 10, 5, &
      target_acceptance=0.234_dp, adaptation_exponent=0.6_dp)
   call check(nonlinear_ekpf_chain%info == 0 .and. &
      all(ieee_is_finite(nonlinear_ekpf_chain%parameters)) .and. &
      all(ieee_is_finite(nonlinear_ekpf_chain%log_likelihood)) .and. &
      all(ieee_is_finite(nonlinear_ekpf_chain%final_proposal_factor)) .and. &
      nonlinear_ekpf_chain%acceptance_rate >= 0.0_dp .and. &
      nonlinear_ekpf_chain%acceptance_rate <= 1.0_dp, &
      'adaptive random-stream multivariate nonlinear EKPF-PMMH')
   nonlinear_da_chain = bssm_multivariate_nonlinear_da_pmmh_draws( &
      reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], [2, 2]), [0.0_dp], &
      reshape([0.0_dp], [1, 1]), [0.0_dp], 2, 1, &
      multivariate_linear_gaussian_observation, &
      pmmh_deterministic_transition_jacobian, pmmh_nonlinear_transition, &
      pmmh_prior, 1, 1.0e-12_dp, reshape([1.0_dp], [1, 1]), &
      pmmh_proposals, pmmh_uniforms, pmmh_uniforms, pmmh_initial_normals, &
      pmmh_innovation_normals, pmmh_nonlinear_resampling)
   call check(nonlinear_da_chain%info == 0 .and. &
      maxval(abs(nonlinear_da_chain%parameters(1, :) - &
      [0.0_dp, 0.0_dp, 0.2_dp])) < 1.0e-14_dp .and. &
      maxval(abs(nonlinear_da_chain%first_stage_probability - &
      [exp(-1.5_dp), exp(-0.06_dp)])) < 1.0e-12_dp .and. &
      maxval(abs(nonlinear_da_chain%second_stage_probability - &
      [0.0_dp, 1.0_dp])) < 1.0e-12_dp .and. &
      all(nonlinear_da_chain%accepted .eqv. [.false., .true.]) .and. &
      maxval(abs(nonlinear_da_chain%coarse_log_likelihood - &
      nonlinear_da_chain%fine_log_likelihood)) < 1.0e-12_dp, &
      'draw-driven IEKF-screened multivariate bootstrap PMMH')
   call set_random_seed(1915)
   nonlinear_da_chain = bssm_multivariate_nonlinear_da_pmmh( &
      reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], [2, 2]), [0.0_dp], &
      reshape([0.1_dp], [1, 1]), [0.0_dp], 2, 1, &
      multivariate_linear_gaussian_observation, &
      pmmh_deterministic_transition_jacobian, pmmh_nonlinear_transition, &
      pmmh_prior, 2, 1.0e-10_dp, reshape([0.1_dp], [1, 1]), 10, 5, &
      target_acceptance=0.234_dp, adaptation_exponent=0.6_dp)
   call check(nonlinear_da_chain%info == 0 .and. &
      all(ieee_is_finite(nonlinear_da_chain%coarse_log_likelihood)) .and. &
      all(ieee_is_finite(nonlinear_da_chain%fine_log_likelihood)) .and. &
      all(ieee_is_finite(nonlinear_da_chain%final_proposal_factor)) .and. &
      nonlinear_da_chain%first_stage_acceptance_rate >= &
      nonlinear_da_chain%acceptance_rate .and. &
      nonlinear_da_chain%acceptance_rate >= 0.0_dp, &
      'adaptive IEKF-screened multivariate bootstrap PMMH')
   nonlinear_ekpf_da_chain = bssm_multivariate_ekpf_da_pmmh_draws( &
      reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], [2, 2]), [0.0_dp], &
      reshape([0.1_dp], [1, 1]), [0.0_dp], 2, 1, &
      multivariate_linear_gaussian_observation, pmmh_psi_transition, &
      pmmh_ekpf_transition, pmmh_prior, 1, 1.0e-12_dp, &
      reshape([1.0_dp], [1, 1]), pmmh_proposals, pmmh_uniforms, &
      pmmh_uniforms, pmmh_initial_normals, pmmh_innovation_normals, &
      pmmh_nonlinear_resampling)
   expected = 0.5_dp + 1.0_dp/1.02_dp
   call check(nonlinear_ekpf_da_chain%info == 0 .and. &
      maxval(abs(nonlinear_ekpf_da_chain%parameters(1, :) - &
      [0.0_dp, 0.0_dp, 0.2_dp])) < 1.0e-14_dp .and. &
      maxval(abs(nonlinear_ekpf_da_chain%first_stage_probability - &
      [exp(-expected), exp(-0.04_dp*expected)])) < 1.0e-12_dp .and. &
      maxval(abs(nonlinear_ekpf_da_chain%second_stage_probability - &
      [0.0_dp, 1.0_dp])) < 1.0e-12_dp .and. &
      all(nonlinear_ekpf_da_chain%accepted .eqv. [.false., .true.]) .and. &
      maxval(abs(nonlinear_ekpf_da_chain%coarse_log_likelihood - &
      nonlinear_ekpf_da_chain%fine_log_likelihood)) < 1.0e-12_dp, &
      'draw-driven IEKF-screened multivariate EKPF-PMMH')
   call set_random_seed(1917)
   nonlinear_ekpf_da_chain = bssm_multivariate_ekpf_da_pmmh( &
      reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], [2, 2]), [0.0_dp], &
      reshape([0.1_dp], [1, 1]), [0.0_dp], 2, 1, &
      multivariate_linear_gaussian_observation, pmmh_psi_transition, &
      pmmh_ekpf_transition, pmmh_prior, 2, 1.0e-10_dp, &
      reshape([0.1_dp], [1, 1]), 10, 5, target_acceptance=0.234_dp, &
      adaptation_exponent=0.6_dp)
   call check(nonlinear_ekpf_da_chain%info == 0 .and. &
      all(ieee_is_finite( &
      nonlinear_ekpf_da_chain%coarse_log_likelihood)) .and. &
      all(ieee_is_finite(nonlinear_ekpf_da_chain%fine_log_likelihood)) .and. &
      all(ieee_is_finite(nonlinear_ekpf_da_chain%final_proposal_factor)) .and. &
      nonlinear_ekpf_da_chain%first_stage_acceptance_rate >= &
      nonlinear_ekpf_da_chain%acceptance_rate .and. &
      nonlinear_ekpf_da_chain%acceptance_rate >= 0.0_dp, &
      'adaptive IEKF-screened multivariate EKPF-PMMH')
   call set_random_seed(1907)
   nonlinear_psi_chain = bssm_nonlinear_psi_pmmh([0.0_dp, 0.0_dp], &
      [0.0_dp], reshape([0.1_dp], [1, 1]), [0.0_dp], 1, &
      pmmh_psi_observation, pmmh_psi_transition, pmmh_prior, 20, &
      1.0e-10_dp, reshape([0.1_dp], [1, 1]), 10, 5, &
      target_acceptance=0.234_dp, adaptation_exponent=0.6_dp)
   call check(nonlinear_psi_chain%info == 0 .and. &
      all(ieee_is_finite(nonlinear_psi_chain%parameters)) .and. &
      all(ieee_is_finite(nonlinear_psi_chain%log_likelihood)) .and. &
      all(ieee_is_finite(nonlinear_psi_chain%final_proposal_factor)) .and. &
      nonlinear_psi_chain%acceptance_rate >= 0.0_dp .and. &
      nonlinear_psi_chain%acceptance_rate <= 1.0_dp, &
      'adaptive random-stream nonlinear psi-filter PMMH')
   call set_random_seed(1909)
   nonlinear_ekpf_chain = bssm_nonlinear_ekpf_pmmh([0.0_dp, 0.0_dp], &
      [0.0_dp], reshape([0.1_dp], [1, 1]), [0.0_dp], 1, &
      pmmh_psi_observation, pmmh_ekpf_transition, pmmh_prior, 1, &
      1.0e-10_dp, reshape([0.1_dp], [1, 1]), 10, 5, &
      target_acceptance=0.234_dp, adaptation_exponent=0.6_dp)
   call check(nonlinear_ekpf_chain%info == 0 .and. &
      all(ieee_is_finite(nonlinear_ekpf_chain%parameters)) .and. &
      all(ieee_is_finite(nonlinear_ekpf_chain%log_likelihood)) .and. &
      all(ieee_is_finite(nonlinear_ekpf_chain%final_proposal_factor)) .and. &
      nonlinear_ekpf_chain%acceptance_rate >= 0.0_dp .and. &
      nonlinear_ekpf_chain%acceptance_rate <= 1.0_dp, &
      'adaptive random-stream nonlinear EKF-proposal PMMH')
   approximate_chain = bssm_approximate_mcmc_draws([0.0_dp], pmmh_prior, &
      pmmh_approximate_likelihood, reshape([1.0_dp], [1, 1]), &
      pmmh_proposals, pmmh_uniforms)
   call check(approximate_chain%info == 0 .and. &
      maxval(abs(approximate_chain%parameters(1, :) - &
      [0.0_dp, 0.0_dp, 0.2_dp])) < 1.0e-14_dp .and. &
      maxval(abs(approximate_chain%acceptance_probability - &
      [exp(-1.0_dp), exp(-0.04_dp)])) < 1.0e-13_dp, &
      'draw-driven deterministic approximate MCMC')
   call set_random_seed(1911)
   nonlinear_approximate_chain = bssm_nonlinear_approximate_mcmc( &
      [0.0_dp, 0.0_dp], [0.0_dp], reshape([0.1_dp], [1, 1]), [0.0_dp], &
      1, pmmh_psi_observation, pmmh_psi_transition, pmmh_prior, 20, &
      1.0e-10_dp, reshape([0.1_dp], [1, 1]), 5, &
      target_acceptance=0.234_dp, adaptation_exponent=0.6_dp)
   call check(nonlinear_approximate_chain%info == 0 .and. &
      all(ieee_is_finite(nonlinear_approximate_chain%log_likelihood)) .and. &
      all(ieee_is_finite( &
      nonlinear_approximate_chain%final_proposal_factor)), &
      'nonlinear Gaussian-approximation MCMC')
   call set_random_seed(1913)
   nonlinear_ekf_chain = bssm_nonlinear_ekf_mcmc([0.0_dp, 0.0_dp], &
      [0.0_dp], reshape([0.1_dp], [1, 1]), [0.0_dp], 1, &
      pmmh_psi_observation, pmmh_psi_transition, pmmh_prior, 1, &
      1.0e-10_dp, reshape([0.1_dp], [1, 1]), 5)
   call check(nonlinear_ekf_chain%info == 0 .and. &
      all(ieee_is_finite(nonlinear_ekf_chain%log_likelihood)), &
      'nonlinear extended-Kalman approximate MCMC')
   correction_is1 = bssm_importance_post_correction( &
      reshape([0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp], [1, 4]), &
      [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], &
      reshape([0.0_dp, log(3.0_dp), log(2.0_dp), 0.0_dp, &
      log(3.0_dp), 0.0_dp, log(4.0_dp), 0.0_dp], [2, 4]), &
      [.false., .true., .false.], bssm_is1)
   correction_is2 = bssm_importance_post_correction( &
      reshape([0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp], [1, 4]), &
      [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], &
      reshape([0.0_dp, log(3.0_dp), log(2.0_dp), 0.0_dp, &
      log(3.0_dp), 0.0_dp, log(4.0_dp), 0.0_dp], [2, 4]), &
      [.false., .true., .false.], bssm_is2)
   correction_is3 = bssm_importance_post_correction( &
      reshape([0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp], [1, 4]), &
      [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], &
      reshape([0.0_dp, log(3.0_dp), log(2.0_dp), 0.0_dp, &
      log(3.0_dp), 0.0_dp, log(4.0_dp), 0.0_dp], [2, 4]), &
      [.false., .true., .false.], bssm_is3)
   call check(correction_is1%info == 0 .and. &
      maxval(abs(correction_is1%weight - 0.25_dp)) < 1.0e-14_dp .and. &
      abs(correction_is1%parameter_mean(1) - 0.5_dp) < 1.0e-14_dp .and. &
      abs(correction_is1%effective_sample_size - 4.0_dp) < 1.0e-14_dp, &
      'IS1 block-proportional post-correction')
   call check(correction_is2%info == 0 .and. &
      maxval(abs(correction_is2%weight - &
      [0.125_dp, 0.125_dp, 0.375_dp, 0.375_dp])) < 1.0e-14_dp .and. &
      abs(correction_is2%parameter_mean(1) - 0.75_dp) < 1.0e-14_dp .and. &
      abs(correction_is2%effective_sample_size - 3.2_dp) < 1.0e-13_dp, &
      'IS2 jump-chain post-correction')
   call check(correction_is3%info == 0 .and. &
      maxval(abs(correction_is3%weight - &
      [0.1_dp, 0.2_dp, 0.3_dp, 0.4_dp])) < 1.0e-14_dp .and. &
      abs(correction_is3%parameter_mean(1) - 0.7_dp) < 1.0e-14_dp .and. &
      abs(correction_is3%effective_sample_size - &
      10.0_dp/3.0_dp) < 1.0e-13_dp, 'IS3 independent post-correction')
   post_initial_normals = 0.0_dp
   post_state_normals = 0.0_dp
   post_terminal_normals = 0.0_dp
   post_resampling_uniforms = 0.5_dp
   correction_is1 = bssm_nonlinear_psi_post_correction_draws( &
      [0.0_dp, 0.0_dp], [0.0_dp], reshape([0.1_dp], [1, 1]), &
      reshape([0.0_dp, 0.1_dp, 0.2_dp], [1, 3]), &
      [0.0_dp, 0.0_dp, 0.0_dp], [.true., .true.], bssm_is1, 1, &
      pmmh_psi_observation, pmmh_psi_transition, 20, 1.0e-12_dp, &
      post_initial_normals, post_state_normals, post_terminal_normals, &
      post_resampling_uniforms)
   call check(correction_is1%info == 0 .and. &
      abs(sum(correction_is1%weight) - 1.0_dp) < 1.0e-14_dp .and. &
      all(ieee_is_finite(correction_is1%weight)) .and. &
      ieee_is_finite(correction_is1%effective_sample_size), &
      'draw-driven nonlinear psi IS1 post-correction')
   call set_random_seed(1921)
   correction_is2 = bssm_nonlinear_psi_post_correction( &
      [0.0_dp, 0.0_dp], [0.0_dp], reshape([0.1_dp], [1, 1]), &
      reshape([0.0_dp, 0.1_dp, 0.2_dp], [1, 3]), &
      [0.0_dp, 0.0_dp, 0.0_dp], [.true., .true.], bssm_is2, 1, &
      pmmh_psi_observation, pmmh_psi_transition, 20, 1.0e-10_dp, 8, 1)
   call check(correction_is2%info == 0 .and. &
      abs(sum(correction_is2%weight) - 1.0_dp) < 1.0e-14_dp .and. &
      all(ieee_is_finite(correction_is2%parameter_mean)), &
      'shared-stream nonlinear psi IS2 post-correction')
   correction_is3 = &
      bssm_multivariate_nonlinear_psi_post_correction_draws( &
      reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], [2, 2]), [0.0_dp], &
      reshape([0.1_dp], [1, 1]), &
      reshape([0.0_dp, 0.1_dp, 0.2_dp], [1, 3]), &
      [0.0_dp, 0.0_dp, 0.0_dp], [.true., .true.], bssm_is3, 2, 1, &
      multivariate_linear_gaussian_observation, pmmh_psi_transition, &
      20, 1.0e-12_dp, post_initial_normals, post_state_normals, &
      post_terminal_normals, post_resampling_uniforms)
   call check(correction_is3%info == 0 .and. &
      abs(sum(correction_is3%weight) - 1.0_dp) < 1.0e-14_dp .and. &
      all(ieee_is_finite(correction_is3%weight)), &
      'draw-driven multivariate nonlinear psi IS3 post-correction')
   call set_random_seed(1923)
   correction_is2 = bssm_multivariate_nonlinear_psi_post_correction( &
      reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], [2, 2]), [0.0_dp], &
      reshape([0.1_dp], [1, 1]), &
      reshape([0.0_dp, 0.1_dp, 0.2_dp], [1, 3]), &
      [0.0_dp, 0.0_dp, 0.0_dp], [.true., .true.], bssm_is2, 2, 1, &
      multivariate_linear_gaussian_observation, pmmh_psi_transition, &
      20, 1.0e-10_dp, 8, 1)
   call check(correction_is2%info == 0 .and. &
      abs(sum(correction_is2%weight) - 1.0_dp) < 1.0e-14_dp .and. &
      all(ieee_is_finite(correction_is2%parameter_covariance)), &
      'shared-stream multivariate nonlinear psi IS2 post-correction')
   corrected_states = bssm_corrected_state_moments( &
      reshape([0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp], [1, 2, 2]), &
      reshape([1.0_dp, 1.0_dp, 3.0_dp, 3.0_dp], [1, 1, 2, 2]), &
      [0.25_dp, 0.75_dp])
   call check(corrected_states%info == 0 .and. &
      maxval(abs(corrected_states%mean(1, :) - &
      [1.5_dp, 2.5_dp])) < 1.0e-14_dp .and. &
      maxval(abs(corrected_states%covariance(1, 1, :) - &
      3.25_dp)) < 1.0e-14_dp, 'post-corrected total state covariance')
   corrected_signals = bssm_linear_signal_moments(corrected_states, &
      reshape([2.0_dp, 2.0_dp], [1, 1, 2]), &
      reshape([1.0_dp, 1.0_dp], [1, 2]))
   call check(corrected_signals%info == 0 .and. &
      maxval(abs(corrected_signals%mean(1, :) - &
      [4.0_dp, 6.0_dp])) < 1.0e-14_dp .and. &
      maxval(abs(corrected_signals%covariance(1, 1, :) - &
      13.0_dp)) < 1.0e-14_dp, 'corrected fitted and predictive signals')
   corrected_trajectories = bssm_corrected_trajectory_draws( &
      reshape([10.0_dp, 11.0_dp, 20.0_dp, 21.0_dp, 30.0_dp, 31.0_dp, &
      40.0_dp, 41.0_dp], [1, 2, 2, 2]), &
      reshape([0.8_dp, 0.2_dp, 0.1_dp, 0.9_dp], [2, 2]), &
      [0.25_dp, 0.75_dp], reshape([0.1_dp, 0.9_dp, 0.5_dp, 0.5_dp], [2, 2]))
   call check(corrected_trajectories%info == 0 .and. &
      all(corrected_trajectories%source_sample == [1, 2]) .and. &
      all(corrected_trajectories%source_particle == [2, 2]) .and. &
      maxval(abs(corrected_trajectories%trajectories(1, :, :) - &
      reshape([20.0_dp, 21.0_dp, 40.0_dp, 41.0_dp], [2, 2]))) < &
      1.0e-14_dp, 'post-corrected trajectory resampling')
   particle_suggestion = bssm_suggest_particles([10, 20, 40], &
      reshape([0.0_dp, 2.0_dp, 4.0_dp, 0.0_dp, 0.5_dp, 1.0_dp, &
      0.0_dp, 0.2_dp, 0.4_dp], [3, 3]))
   call check(particle_suggestion%info == 0 .and. &
      maxval(abs(particle_suggestion%standard_deviation - &
      [2.0_dp, 0.5_dp, 0.2_dp])) < 1.0e-14_dp .and. &
      particle_suggestion%selected_index == 2 .and. &
      particle_suggestion%particle_count == 20, &
      'particle-count log-likelihood variability diagnostic')
   suggest_initial_normals = 0.0_dp
   suggest_state_normals = 0.0_dp
   suggest_terminal_normals = 0.0_dp
   suggest_resampling_uniforms = 0.5_dp
   particle_suggestion = bssm_nonlinear_psi_suggest_particles_draws( &
      [0.0_dp, 0.0_dp], [0.0_dp], reshape([0.1_dp], [1, 1]), &
      [0.1_dp], 1, pmmh_psi_observation, pmmh_psi_transition, &
      20, 1.0e-12_dp, [2, 4], suggest_initial_normals, &
      suggest_state_normals, suggest_terminal_normals, &
      suggest_resampling_uniforms)
   call check(particle_suggestion%info == 0 .and. &
      particle_suggestion%selected_index == 1 .and. &
      particle_suggestion%particle_count == 2 .and. &
      maxval(abs(particle_suggestion%standard_deviation)) < 1.0e-14_dp, &
      'draw-driven nonlinear psi particle suggestion')
   call set_random_seed(1925)
   particle_suggestion = bssm_nonlinear_psi_suggest_particles( &
      [0.0_dp, 0.0_dp], [0.0_dp], reshape([0.1_dp], [1, 1]), &
      [0.1_dp], 1, pmmh_psi_observation, pmmh_psi_transition, &
      20, 1.0e-10_dp, [2, 4], 2, 1.0e6_dp)
   call check(particle_suggestion%info == 0 .and. &
      particle_suggestion%particle_count == 2 .and. &
      all(ieee_is_finite(particle_suggestion%standard_deviation)), &
      'shared-stream nonlinear psi particle suggestion')
   particle_suggestion = &
      bssm_multivariate_nonlinear_psi_suggest_particles_draws( &
      reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], [2, 2]), [0.0_dp], &
      reshape([0.1_dp], [1, 1]), [0.1_dp], 2, 1, &
      multivariate_linear_gaussian_observation, pmmh_psi_transition, &
      20, 1.0e-12_dp, [2, 4], suggest_initial_normals, &
      suggest_state_normals, suggest_terminal_normals, &
      suggest_resampling_uniforms)
   call check(particle_suggestion%info == 0 .and. &
      particle_suggestion%selected_index == 1 .and. &
      particle_suggestion%particle_count == 2 .and. &
      maxval(abs(particle_suggestion%standard_deviation)) < 1.0e-14_dp, &
      'draw-driven multivariate nonlinear psi particle suggestion')
   call set_random_seed(1927)
   particle_suggestion = bssm_multivariate_nonlinear_psi_suggest_particles( &
      reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], [2, 2]), [0.0_dp], &
      reshape([0.1_dp], [1, 1]), [0.1_dp], 2, 1, &
      multivariate_linear_gaussian_observation, pmmh_psi_transition, &
      20, 1.0e-10_dp, [2, 4], 2, 1.0e6_dp)
   call check(particle_suggestion%info == 0 .and. &
      particle_suggestion%particle_count == 2 .and. &
      all(ieee_is_finite(particle_suggestion%standard_deviation)), &
      'shared-stream multivariate nonlinear psi particle suggestion')
   call set_random_seed(1915)
   call random_number(diagnostic_samples(1, :))
   diagnostic_samples(1, :) = diagnostic_samples(1, :) - 0.5_dp
   do diagnostic_index = 2, 200
      diagnostic_samples(1, diagnostic_index) = &
         0.8_dp*diagnostic_samples(1, diagnostic_index - 1) + &
         diagnostic_samples(1, diagnostic_index)
   end do
   diagnostic_samples(2, :) = 2.0_dp*diagnostic_samples(1, :) + 1.0_dp
   diagnostic_weights = [(1.0_dp + 0.002_dp*real(diagnostic_index, dp), &
      diagnostic_index=1, 200)]
   sokal_diagnostics = bssm_mcmc_diagnostics(diagnostic_samples, bssm_sokal)
   geyer_diagnostics = bssm_mcmc_diagnostics(diagnostic_samples, bssm_geyer)
   weighted_diagnostics = bssm_mcmc_diagnostics(diagnostic_samples, &
      bssm_sokal, diagnostic_weights)
   chain_diagnostics = bssm_chain_diagnostics(approximate_chain, bssm_sokal)
   da_chain_diagnostics = bssm_da_chain_diagnostics(random_da_chain, &
      bssm_sokal)
   call check(sokal_diagnostics%info == 0 .and. &
      sokal_diagnostics%iact(1) > 1.0_dp .and. &
      abs(sokal_diagnostics%iact(1) - &
      bssm_integrated_autocorrelation_time(diagnostic_samples(1, :), &
      bssm_sokal)) < 1.0e-12_dp .and. &
      abs(sokal_diagnostics%asymptotic_variance(1) - &
      bssm_asymptotic_variance(diagnostic_samples(1, :), &
      method=bssm_sokal)) < 1.0e-14_dp .and. &
      abs(sokal_diagnostics%ess(1) - &
      200.0_dp/sokal_diagnostics%iact(1)) < 1.0e-10_dp, &
      'Sokal integrated autocorrelation and ESS')
   call check(abs(sokal_diagnostics%mean(2) - &
      (2.0_dp*sokal_diagnostics%mean(1) + 1.0_dp)) < 1.0e-14_dp .and. &
      abs(sokal_diagnostics%variance(2) - &
      4.0_dp*sokal_diagnostics%variance(1)) < 1.0e-13_dp .and. &
      abs(sokal_diagnostics%mcse(2) - &
      2.0_dp*sokal_diagnostics%mcse(1)) < 1.0e-13_dp .and. &
      abs(sokal_diagnostics%ess(2) - sokal_diagnostics%ess(1)) < 1.0e-10_dp, &
      'affine-invariant MCMC diagnostics')
   call check(geyer_diagnostics%info == 0 .and. &
      all(ieee_is_finite(geyer_diagnostics%iact)) .and. &
      all(geyer_diagnostics%iact >= 0.0_dp) .and. &
      all(geyer_diagnostics%mcse >= 0.0_dp), &
      'Geyer initial-monotone diagnostics')
   call check(weighted_diagnostics%info == 0 .and. &
      all(ieee_is_finite(weighted_diagnostics%asymptotic_variance)) .and. &
      abs(weighted_diagnostics%mean(2) - &
      (2.0_dp*weighted_diagnostics%mean(1) + 1.0_dp)) < 1.0e-14_dp .and. &
      abs(weighted_diagnostics%ess(2) - weighted_diagnostics%ess(1)) < &
      1.0e-10_dp .and. chain_diagnostics%info == 0 .and. &
      da_chain_diagnostics%info == 0, &
      'weighted IS-MCMC and retained-chain diagnostics')

   initial_normals = 0.0_dp
   innovations = 0.0_dp
   uniforms = 0.5_dp
   filter = bssm_bootstrap_filter_draws([0.0_dp, 1.0_dp], [1.0_dp], &
      reshape([1.0_dp], [1, 1]), reshape([0.0_dp], [1, 1]), [0.0_dp], &
      reshape([0.0_dp], [1, 1]), bssm_poisson, 1.0_dp, initial_normals, &
      innovations, uniforms)
   call check(filter%info == 0 .and. &
      all(shape(filter%particles) == [1, 4, 3]) .and. &
      all(shape(filter%weights) == [4, 2]) .and. &
      maxval(abs(filter%weights - 0.25_dp)) < 1.0e-14_dp .and. &
      abs(filter%log_likelihood + 2.0_dp) < 1.0e-13_dp .and. &
      maxval(abs(filter%filtered_mean)) < 1.0e-14_dp, &
      'draw-driven bootstrap filter')

   initial_normals(1, :) = [-1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp]
   filter = bssm_bootstrap_filter_draws([0.0_dp, 0.0_dp], [1.0_dp], &
      reshape([1.0_dp], [1, 1]), reshape([0.0_dp], [1, 1]), [0.0_dp], &
      reshape([1.0_dp], [1, 1]), bssm_gaussian, 1.0_dp, initial_normals, &
      innovations, uniforms)
   call check(filter%info == 0 .and. &
      maxval(abs(sum(filter%weights, dim=1) - 1.0_dp)) < 1.0e-14_dp .and. &
      all(filter%ancestors >= 1) .and. all(filter%ancestors <= 4) .and. &
      abs(filter%filtered_mean(1, 1)) < &
      abs(filter%predicted_mean(1, 1)), 'weighted particle summaries')
   smoother = bssm_particle_smoother(filter)
   call check(smoother%info == 0 .and. &
      all(shape(smoother%trajectories) == [1, 2, 4]) .and. &
      all(shape(smoother%covariance) == [1, 1, 2]) .and. &
      all(ieee_is_finite(smoother%mean)) .and. &
      all(smoother%covariance >= 0.0_dp), 'genealogical particle smoother')

   call set_random_seed(481)
   random_filter = bssm_bootstrap_filter([0.0_dp, 1.0_dp], [1.0_dp], &
      reshape([1.0_dp], [1, 1]), reshape([0.1_dp], [1, 1]), [0.0_dp], &
      reshape([1.0_dp], [1, 1]), bssm_poisson, 1.0_dp, 20)
   call check(random_filter%info == 0 .and. &
      ieee_is_finite(random_filter%log_likelihood) .and. &
      all(ieee_is_finite(random_filter%filtered_mean)), &
      'random-stream bootstrap filter')

   initial_normals = 0.0_dp
   observation_schedule = 1.0_dp
   transition_schedule(1, 1, :) = [1.0_dp, 2.0_dp]
   noise_schedule = 0.0_dp
   state_offset_schedule(1, :) = [0.0_dp, 0.5_dp]
   filter = bssm_bootstrap_filter_tv_draws([1.0_dp, 1.0_dp], &
      observation_schedule, transition_schedule, noise_schedule, [1.0_dp], &
      reshape([0.0_dp], [1, 1]), bssm_gaussian, 1.0_dp, initial_normals, &
      innovations, uniforms, state_offset=state_offset_schedule)
   expected = -log(2.0_dp*acos(-1.0_dp))
   call check(filter%info == 0 .and. &
      maxval(abs(filter%predicted_mean(1, :) - &
      [1.0_dp, 1.0_dp, 2.5_dp])) < 1.0e-14_dp .and. &
      abs(filter%log_likelihood - expected) < 1.0e-13_dp, &
      'time-varying linear bootstrap filter')

   filter = bssm_nonlinear_bootstrap_filter_draws([4.0_dp, 2.25_dp], &
      [2.0_dp], reshape([0.0_dp], [1, 1]), &
      [0.5_dp, 0.5_dp, 0.0_dp], 1, nonlinear_log_density, &
      nonlinear_transition, initial_normals, innovations, uniforms)
   call check(filter%info == 0 .and. &
      maxval(abs(filter%predicted_mean(1, :) - &
      [2.0_dp, 1.5_dp, 1.25_dp])) < 1.0e-14_dp .and. &
      abs(filter%log_likelihood - expected) < 1.0e-13_dp, &
      'nonlinear bootstrap filter')

   initial_normals(1, :) = [-1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp]
   ekpf_proposals = 0.0_dp
   ekpf_uniforms = 0.5_dp
   filter = bssm_ekpf_draws([1.0_dp], [0.0_dp], &
      reshape([1.0_dp], [1, 1]), [1.0_dp, 0.0_dp, 1.0_dp], 1, &
      linear_gaussian_observation, nonlinear_transition, initial_normals, &
      ekpf_proposals, ekpf_uniforms)
   expected = -0.5_dp*(log(4.0_dp*acos(-1.0_dp)) + 0.5_dp)
   call check(filter%info == 0 .and. &
      maxval(abs(filter%weights(:, 1) - 0.25_dp)) < 1.0e-12_dp .and. &
      abs(filter%log_likelihood - expected) < 1.0e-12_dp .and. &
      all(filter%ancestors(:, 1) == [1, 2, 3, 4]), &
      'extended Kalman particle filter')
   call set_random_seed(912)
   random_filter = bssm_ekpf([1.0_dp], [0.0_dp], &
      reshape([1.0_dp], [1, 1]), [1.0_dp, 0.0_dp, 1.0_dp], 1, &
      linear_gaussian_observation, nonlinear_transition, 20)
   call check(random_filter%info == 0 .and. &
      ieee_is_finite(random_filter%log_likelihood), &
      'random-stream extended Kalman particle filter')
   nonlinear_ekpf_proposals = 0.0_dp
   filter = bssm_ekpf_draws([1.2_dp, 1.1_dp], [0.0_dp], &
      reshape([0.5_dp], [1, 1]), [0.8_dp, 0.1_dp, 0.2_dp], 1, &
      exponential_gaussian_observation, nonlinear_transition, initial_normals, &
      nonlinear_ekpf_proposals, uniforms)
   call check(filter%info == 0 .and. &
      ieee_is_finite(filter%log_likelihood) .and. &
      maxval(abs(sum(filter%weights, dim=1) - 1.0_dp)) < 1.0e-13_dp .and. &
      all(ieee_is_finite(filter%filtered_mean)), &
      'nonlinear extended Kalman particle filter')

   ekf_fit = bssm_iekf([1.0_dp], [0.0_dp], &
      reshape([1.0_dp], [1, 1]), [1.0_dp, 0.0_dp, 1.0_dp], 1, &
      linear_gaussian_observation, nonlinear_transition_jacobian, 5, &
      1.0e-12_dp)
   call check(ekf_fit%info == 0 .and. &
      abs(ekf_fit%filtered_mean(1, 1) - 0.5_dp) < 1.0e-13_dp .and. &
      abs(ekf_fit%filtered_covariance(1, 1, 1) - 0.5_dp) < 1.0e-13_dp .and. &
      abs(ekf_fit%predicted_covariance(1, 1, 2) - 1.5_dp) < 1.0e-13_dp .and. &
      abs(ekf_fit%log_likelihood - expected) < 1.0e-12_dp .and. &
      ekf_fit%iterations(1) == 1, 'iterated EKF linear equivalence')
   ekf_smoother = bssm_ekf_smoother([1.0_dp, 2.0_dp], [0.0_dp], &
      reshape([1.0_dp], [1, 1]), [1.0_dp, 0.0_dp, 1.0_dp], 1, &
      linear_gaussian_observation, nonlinear_transition_jacobian, 0, &
      1.0e-12_dp)
   fast_ekf_smoother = bssm_ekf_fast_smoother([1.0_dp, 2.0_dp], [0.0_dp], &
      reshape([1.0_dp], [1, 1]), [1.0_dp, 0.0_dp, 1.0_dp], 1, &
      linear_gaussian_observation, nonlinear_transition_jacobian, 0, &
      1.0e-12_dp)
   call check(ekf_smoother%info == 0 .and. &
      maxval(abs(ekf_smoother%state(1, :) - [0.8_dp, 1.4_dp])) < &
      1.0e-12_dp .and. &
      maxval(abs(ekf_smoother%covariance(1, 1, :) - [0.4_dp, 0.6_dp])) < &
      1.0e-12_dp, 'extended Kalman covariance smoother')
   call check(fast_ekf_smoother%info == 0 .and. &
      maxval(abs(fast_ekf_smoother%state - ekf_smoother%state)) < &
      1.0e-12_dp .and. .not. allocated(fast_ekf_smoother%covariance) .and. &
      abs(fast_ekf_smoother%log_likelihood - &
      ekf_smoother%log_likelihood) < 1.0e-12_dp, &
      'extended Kalman fast smoother')
   ekf_fit = bssm_iekf([2.0_dp], [0.0_dp], &
      reshape([0.5_dp], [1, 1]), [0.8_dp, 0.1_dp, 0.2_dp], 1, &
      exponential_gaussian_observation, nonlinear_transition_jacobian, 0, &
      1.0e-10_dp)
   iekf_fit = bssm_iekf([2.0_dp], [0.0_dp], &
      reshape([0.5_dp], [1, 1]), [0.8_dp, 0.1_dp, 0.2_dp], 1, &
      exponential_gaussian_observation, nonlinear_transition_jacobian, 20, &
      1.0e-10_dp)
   call check(iekf_fit%info == 0 .and. iekf_fit%iterations(1) > 0 .and. &
      abs(iekf_fit%filtered_mean(1, 1) - log(2.0_dp)) < &
      abs(ekf_fit%filtered_mean(1, 1) - log(2.0_dp)), &
      'iterated EKF nonlinear refinement')
   ekf_smoother = bssm_ekf_smoother([2.0_dp, 1.5_dp], [0.0_dp], &
      reshape([0.5_dp], [1, 1]), [0.8_dp, 0.1_dp, 0.2_dp], 1, &
      exponential_gaussian_observation, nonlinear_transition_jacobian, 20, &
      1.0e-10_dp)
   fast_ekf_smoother = bssm_ekf_fast_smoother([2.0_dp, 1.5_dp], [0.0_dp], &
      reshape([0.5_dp], [1, 1]), [0.8_dp, 0.1_dp, 0.2_dp], 1, &
      exponential_gaussian_observation, nonlinear_transition_jacobian, 20, &
      1.0e-10_dp)
   call check(ekf_smoother%info == 0 .and. fast_ekf_smoother%info == 0 .and. &
      all(ieee_is_finite(ekf_smoother%state)) .and. &
      all(ieee_is_finite(ekf_smoother%covariance)) .and. &
      maxval(abs(fast_ekf_smoother%state - ekf_smoother%state)) < &
      1.0e-12_dp .and. any(ekf_smoother%iterations > 0), &
      'iterated nonlinear extended Kalman smoother')
   multivariate_ekf_fit = bssm_multivariate_iekf(reshape([1.0_dp, 3.0_dp, &
      2.0_dp, ieee_value(0.0_dp, ieee_quiet_nan)], [2, 2]), [0.0_dp], &
      reshape([1.0_dp], [1, 1]), [1.0_dp, 0.0_dp, 1.0_dp], 2, 1, &
      multivariate_linear_gaussian_observation, &
      nonlinear_transition_jacobian, 5, 1.0e-12_dp)
   call check(multivariate_ekf_fit%info == 0 .and. &
      maxval(abs(multivariate_ekf_fit%filtered_mean(1, :) - &
      [4.0_dp/3.0_dp, 12.0_dp/7.0_dp])) < 1.0e-12_dp .and. &
      maxval(abs(multivariate_ekf_fit%filtered_covariance(1, 1, :) - &
      [1.0_dp/3.0_dp, 4.0_dp/7.0_dp])) < 1.0e-12_dp .and. &
      multivariate_ekf_fit%innovation(2, 2) == 0.0_dp .and. &
      all(multivariate_ekf_fit%innovation_covariance(2, :, 2) == 0.0_dp) .and. &
      all(multivariate_ekf_fit%iterations == 1) .and. &
      ieee_is_finite(multivariate_ekf_fit%log_likelihood), &
      'partially observed multivariate iterated EKF')
   ekf_smoother = bssm_multivariate_ekf_smoother( &
      reshape([1.0_dp, 3.0_dp, 2.0_dp, &
      ieee_value(0.0_dp, ieee_quiet_nan)], [2, 2]), [0.0_dp], &
      reshape([1.0_dp], [1, 1]), [1.0_dp, 0.0_dp, 1.0_dp], 2, 1, &
      multivariate_linear_gaussian_observation, &
      nonlinear_transition_jacobian, 5, 1.0e-12_dp)
   fast_ekf_smoother = bssm_multivariate_ekf_fast_smoother( &
      reshape([1.0_dp, 3.0_dp, 2.0_dp, &
      ieee_value(0.0_dp, ieee_quiet_nan)], [2, 2]), [0.0_dp], &
      reshape([1.0_dp], [1, 1]), [1.0_dp, 0.0_dp, 1.0_dp], 2, 1, &
      multivariate_linear_gaussian_observation, &
      nonlinear_transition_jacobian, 5, 1.0e-12_dp)
   call check(ekf_smoother%info == 0 .and. &
      maxval(abs(ekf_smoother%state(1, :) - &
      [10.0_dp/7.0_dp, 12.0_dp/7.0_dp])) < 1.0e-12_dp .and. &
      maxval(abs(ekf_smoother%covariance(1, 1, :) - &
      [2.0_dp/7.0_dp, 4.0_dp/7.0_dp])) < 1.0e-12_dp .and. &
      fast_ekf_smoother%info == 0 .and. &
      maxval(abs(fast_ekf_smoother%state - ekf_smoother%state)) < &
      1.0e-12_dp .and. .not. allocated(fast_ekf_smoother%covariance), &
      'multivariate extended Kalman smoothers')
   multivariate_ekf_fit = bssm_multivariate_iekf(reshape([2.0_dp, 1.0_dp, &
      1.5_dp, ieee_value(0.0_dp, ieee_quiet_nan)], [2, 2]), [0.0_dp], &
      reshape([0.5_dp], [1, 1]), [0.8_dp, 0.1_dp, 0.2_dp], 2, 1, &
      multivariate_nonlinear_gaussian_observation, &
      nonlinear_transition_jacobian, 20, 1.0e-10_dp)
   ekf_smoother = bssm_multivariate_ekf_smoother( &
      reshape([2.0_dp, 1.0_dp, 1.5_dp, &
      ieee_value(0.0_dp, ieee_quiet_nan)], [2, 2]), [0.0_dp], &
      reshape([0.5_dp], [1, 1]), [0.8_dp, 0.1_dp, 0.2_dp], 2, 1, &
      multivariate_nonlinear_gaussian_observation, &
      nonlinear_transition_jacobian, 20, 1.0e-10_dp)
   fast_ekf_smoother = bssm_multivariate_ekf_fast_smoother( &
      reshape([2.0_dp, 1.0_dp, 1.5_dp, &
      ieee_value(0.0_dp, ieee_quiet_nan)], [2, 2]), [0.0_dp], &
      reshape([0.5_dp], [1, 1]), [0.8_dp, 0.1_dp, 0.2_dp], 2, 1, &
      multivariate_nonlinear_gaussian_observation, &
      nonlinear_transition_jacobian, 20, 1.0e-10_dp)
   call check(multivariate_ekf_fit%info == 0 .and. &
      any(multivariate_ekf_fit%iterations > 0) .and. &
      all(ieee_is_finite(multivariate_ekf_fit%filtered_mean)) .and. &
      ekf_smoother%info == 0 .and. fast_ekf_smoother%info == 0 .and. &
      all(ieee_is_finite(ekf_smoother%covariance)) .and. &
      maxval(abs(fast_ekf_smoother%state - ekf_smoother%state)) < &
      1.0e-12_dp, 'iterated multivariate nonlinear EKF smoothing')
   multivariate_initial_normals = 0.0_dp
   multivariate_innovation_normals = 0.0_dp
   multivariate_proposal_normals = 0.0_dp
   multivariate_resampling_uniforms = 0.5_dp
   filter = bssm_multivariate_nonlinear_bootstrap_filter_draws( &
      reshape([1.0_dp, 3.0_dp, 2.0_dp, &
      ieee_value(0.0_dp, ieee_quiet_nan)], [2, 2]), [0.0_dp], &
      reshape([1.0_dp], [1, 1]), [1.0_dp, 0.0_dp, 1.0_dp], 2, 1, &
      multivariate_linear_gaussian_observation, nonlinear_transition, &
      multivariate_initial_normals, multivariate_innovation_normals, &
      multivariate_resampling_uniforms)
   call check(filter%info == 0 .and. &
      maxval(abs(filter%weights - 0.25_dp)) < 1.0e-14_dp .and. &
      all(filter%ancestors >= 1) .and. all(filter%ancestors <= 4) .and. &
      all(ieee_is_finite(filter%filtered_mean)) .and. &
      ieee_is_finite(filter%log_likelihood), &
      'multivariate nonlinear bootstrap particle filter')
   filter = bssm_multivariate_ekpf_draws( &
      reshape([2.0_dp, 1.0_dp, 1.5_dp, &
      ieee_value(0.0_dp, ieee_quiet_nan)], [2, 2]), [0.0_dp], &
      reshape([0.5_dp], [1, 1]), [0.8_dp, 0.1_dp, 0.2_dp], 2, 1, &
      multivariate_nonlinear_gaussian_observation, nonlinear_transition, &
      multivariate_initial_normals, multivariate_proposal_normals, &
      multivariate_resampling_uniforms, 20, 1.0e-10_dp)
   call check(filter%info == 0 .and. &
      maxval(abs(filter%weights - 0.25_dp)) < 1.0e-13_dp .and. &
      all(filter%ancestors >= 1) .and. all(filter%ancestors <= 4) .and. &
      all(ieee_is_finite(filter%particles)) .and. &
      ieee_is_finite(filter%log_likelihood), &
      'multivariate nonlinear EKF-proposal particle filter')
   call set_random_seed(1949)
   random_filter = bssm_multivariate_nonlinear_bootstrap_filter( &
      reshape([1.0_dp, 3.0_dp, 2.0_dp, &
      ieee_value(0.0_dp, ieee_quiet_nan)], [2, 2]), [0.0_dp], &
      reshape([1.0_dp], [1, 1]), [1.0_dp, 0.0_dp, 1.0_dp], 2, 1, &
      multivariate_linear_gaussian_observation, nonlinear_transition, 4)
   call check(random_filter%info == 0 .and. &
      ieee_is_finite(random_filter%log_likelihood), &
      'shared-stream multivariate nonlinear bootstrap filter')
   random_filter = bssm_multivariate_ekpf( &
      reshape([2.0_dp, 1.0_dp, 1.5_dp, &
      ieee_value(0.0_dp, ieee_quiet_nan)], [2, 2]), [0.0_dp], &
      reshape([0.5_dp], [1, 1]), [0.8_dp, 0.1_dp, 0.2_dp], 2, 1, &
      multivariate_nonlinear_gaussian_observation, nonlinear_transition, 4, &
      20, 1.0e-10_dp)
   call check(random_filter%info == 0 .and. &
      ieee_is_finite(random_filter%log_likelihood), &
      'shared-stream multivariate nonlinear EKPF')
   filter = bssm_ekpf_draws([1.2_dp, 1.1_dp], [0.0_dp], &
      reshape([0.5_dp], [1, 1]), [0.8_dp, 0.1_dp, 0.2_dp], 1, &
      exponential_gaussian_observation, nonlinear_transition, initial_normals, &
      nonlinear_ekpf_proposals, uniforms, 5, 1.0e-10_dp)
   call check(filter%info == 0 .and. &
      ieee_is_finite(filter%log_likelihood), 'iterated EKPF proposal')

   ekf_fit = bssm_iekf([1.0_dp, 1.0_dp], [0.0_dp], &
      reshape([1.0_dp], [1, 1]), [1.0_dp, 0.0_dp, 0.5_dp], 1, &
      linear_gaussian_observation, nonlinear_transition_jacobian, 0, &
      1.0e-12_dp)
   nonlinear_approximation = bssm_nonlinear_gaussian_approximation( &
      [1.0_dp, 1.0_dp], [0.0_dp], reshape([1.0_dp], [1, 1]), &
      [1.0_dp, 0.0_dp, 0.5_dp], 1, linear_gaussian_observation, &
      nonlinear_transition_jacobian, 20, 1.0e-12_dp)
   call check(nonlinear_approximation%info == 0 .and. &
      nonlinear_approximation%converged .and. &
      all(shape(nonlinear_approximation%mode_state) == [1, 3]) .and. &
      maxval(abs(nonlinear_approximation%scaling)) < 1.0e-12_dp .and. &
      maxval(abs(nonlinear_approximation%observation_jacobian - &
      1.0_dp)) < 1.0e-14_dp .and. &
      maxval(abs(nonlinear_approximation%transition_jacobian - &
      1.0_dp)) < 1.0e-14_dp .and. &
      abs(nonlinear_approximation%gaussian_log_likelihood - &
      ekf_fit%log_likelihood) < 1.0e-12_dp, &
      'global iterated-EKS linear equivalence')
   multivariate_ekf_fit = bssm_multivariate_iekf( &
      reshape([1.0_dp, 3.0_dp, 1.0_dp, &
      ieee_value(0.0_dp, ieee_quiet_nan)], [2, 2]), [0.0_dp], &
      reshape([1.0_dp], [1, 1]), [1.0_dp, 0.0_dp, 0.5_dp], 2, 1, &
      multivariate_linear_gaussian_observation, &
      nonlinear_transition_jacobian, 0, 1.0e-12_dp)
   multivariate_nonlinear_approximation = &
      bssm_multivariate_nonlinear_gaussian_approximation( &
      reshape([1.0_dp, 3.0_dp, 1.0_dp, &
      ieee_value(0.0_dp, ieee_quiet_nan)], [2, 2]), [0.0_dp], &
      reshape([1.0_dp], [1, 1]), [1.0_dp, 0.0_dp, 0.5_dp], 2, 1, &
      multivariate_linear_gaussian_observation, &
      nonlinear_transition_jacobian, 20, 1.0e-12_dp)
   call check(multivariate_nonlinear_approximation%info == 0 .and. &
      multivariate_nonlinear_approximation%converged .and. &
      all(shape(multivariate_nonlinear_approximation%mode_state) == &
      [1, 3]) .and. &
      maxval(abs(multivariate_nonlinear_approximation%scaling)) < &
      1.0e-12_dp .and. &
      maxval(abs(multivariate_nonlinear_approximation% &
      observation_jacobian - 1.0_dp)) < 1.0e-14_dp .and. &
      abs(multivariate_nonlinear_approximation%gaussian_log_likelihood - &
      multivariate_ekf_fit%log_likelihood) < 1.0e-12_dp, &
      'global multivariate iterated-EKS linear equivalence')
   multivariate_initial_normals = 0.0_dp
   multivariate_proposal_normals = 0.0_dp
   multivariate_resampling_uniforms = 0.5_dp
   filter = bssm_multivariate_nonlinear_psi_filter_draws( &
      reshape([1.0_dp, 3.0_dp, 1.0_dp, &
      ieee_value(0.0_dp, ieee_quiet_nan)], [2, 2]), [0.0_dp], &
      reshape([1.0_dp], [1, 1]), [1.0_dp, 0.0_dp, 0.5_dp], 2, 1, &
      multivariate_linear_gaussian_observation, &
      nonlinear_transition_jacobian, 20, 1.0e-12_dp, &
      multivariate_initial_normals, multivariate_proposal_normals, &
      multivariate_initial_normals, multivariate_resampling_uniforms)
   call check(filter%info == 0 .and. &
      all(shape(filter%particles) == [1, 4, 3]) .and. &
      maxval(abs(filter%weights - 0.25_dp)) < 1.0e-13_dp .and. &
      abs(filter%log_likelihood - &
      multivariate_nonlinear_approximation%corrected_log_likelihood) < &
      1.0e-12_dp .and. &
      all(filter%ancestors(:, 1) == [1, 2, 3, 4]) .and. &
      all(filter%ancestors(:, 2) == [1, 2, 3, 4]), &
      'draw-driven multivariate nonlinear psi linear equivalence')
   multivariate_nonlinear_approximation = &
      bssm_multivariate_nonlinear_gaussian_approximation( &
      reshape([2.0_dp, 1.0_dp, 1.5_dp, &
      ieee_value(0.0_dp, ieee_quiet_nan)], [2, 2]), [0.0_dp], &
      reshape([0.5_dp], [1, 1]), [0.8_dp, 0.1_dp, 0.2_dp], 2, 1, &
      multivariate_nonlinear_gaussian_observation, &
      curved_transition_jacobian, 100, 1.0e-10_dp)
   call check(multivariate_nonlinear_approximation%info == 0 .and. &
      multivariate_nonlinear_approximation%converged .and. &
      ieee_is_finite(multivariate_nonlinear_approximation% &
      corrected_log_likelihood) .and. &
      all(ieee_is_finite(multivariate_nonlinear_approximation%scaling)), &
      'global multivariate nonlinear Gaussian approximation')
   filter = bssm_multivariate_nonlinear_psi_filter_draws( &
      reshape([2.0_dp, 1.0_dp, 1.5_dp, &
      ieee_value(0.0_dp, ieee_quiet_nan)], [2, 2]), [0.0_dp], &
      reshape([0.5_dp], [1, 1]), [0.8_dp, 0.1_dp, 0.2_dp], 2, 1, &
      multivariate_nonlinear_gaussian_observation, &
      curved_transition_jacobian, 100, 1.0e-10_dp, &
      multivariate_initial_normals, multivariate_proposal_normals, &
      multivariate_initial_normals, multivariate_resampling_uniforms)
   call check(filter%info == 0 .and. &
      ieee_is_finite(filter%log_likelihood) .and. &
      maxval(abs(sum(filter%weights, dim=1) - 1.0_dp)) < 1.0e-13_dp .and. &
      all(ieee_is_finite(filter%particles)), &
      'draw-driven multivariate nonlinear psi corrections')
   call set_random_seed(1599)
   random_filter = bssm_multivariate_nonlinear_psi_filter( &
      reshape([2.0_dp, 1.0_dp, 1.5_dp, &
      ieee_value(0.0_dp, ieee_quiet_nan)], [2, 2]), [0.0_dp], &
      reshape([0.5_dp], [1, 1]), [0.8_dp, 0.1_dp, 0.2_dp], 2, 1, &
      multivariate_nonlinear_gaussian_observation, &
      curved_transition_jacobian, 20, 100, 1.0e-10_dp)
   call check(random_filter%info == 0 .and. &
      ieee_is_finite(random_filter%log_likelihood) .and. &
      all(ieee_is_finite(random_filter%filtered_mean)), &
      'shared-stream multivariate nonlinear psi filter')
   initial_normals = 0.0_dp
   nonlinear_ekpf_proposals = 0.0_dp
   uniforms = 0.5_dp
   filter = bssm_nonlinear_psi_filter_draws([1.0_dp, 1.0_dp], [0.0_dp], &
      reshape([1.0_dp], [1, 1]), [1.0_dp, 0.0_dp, 0.5_dp], 1, &
      linear_gaussian_observation, nonlinear_transition_jacobian, 20, &
      1.0e-12_dp, initial_normals, nonlinear_ekpf_proposals, &
      initial_normals, uniforms)
   call check(filter%info == 0 .and. &
      all(shape(filter%particles) == [1, 4, 3]) .and. &
      maxval(abs(filter%weights - 0.25_dp)) < 1.0e-13_dp .and. &
      abs(filter%log_likelihood - &
      nonlinear_approximation%corrected_log_likelihood) < 1.0e-12_dp .and. &
      all(filter%ancestors(:, 1) == [1, 2, 3, 4]) .and. &
      all(filter%ancestors(:, 2) == [1, 2, 3, 4]), &
      'draw-driven nonlinear psi linear equivalence')
   call set_random_seed(1601)
   random_filter = bssm_nonlinear_psi_filter([1.0_dp, 1.0_dp], [0.0_dp], &
      reshape([1.0_dp], [1, 1]), [1.0_dp, 0.0_dp, 0.5_dp], 1, &
      linear_gaussian_observation, nonlinear_transition_jacobian, 20, 20, &
      1.0e-10_dp)
   call check(random_filter%info == 0 .and. &
      ieee_is_finite(random_filter%log_likelihood) .and. &
      all(ieee_is_finite(random_filter%filtered_mean)), &
      'random-stream nonlinear psi filter')
   nonlinear_approximation = bssm_nonlinear_gaussian_approximation( &
      [1.2_dp, 1.1_dp], [0.0_dp], reshape([0.3_dp], [1, 1]), &
      [0.7_dp, 0.05_dp, 0.2_dp], 1, exponential_gaussian_observation, &
      curved_transition_jacobian, 100, 1.0e-10_dp)
   call check(nonlinear_approximation%info == 0 .and. &
      nonlinear_approximation%converged .and. &
      ieee_is_finite(nonlinear_approximation%corrected_log_likelihood) .and. &
      all(ieee_is_finite(nonlinear_approximation%scaling)), &
      'curved global nonlinear approximation')
   call set_random_seed(1603)
   random_filter = bssm_nonlinear_psi_filter([1.2_dp, 1.1_dp], [0.0_dp], &
      reshape([0.3_dp], [1, 1]), [0.7_dp, 0.05_dp, 0.2_dp], 1, &
      exponential_gaussian_observation, curved_transition_jacobian, 40, &
      100, 1.0e-10_dp)
   call check(random_filter%info == 0 .and. &
      ieee_is_finite(random_filter%log_likelihood) .and. &
      maxval(abs(sum(random_filter%weights, dim=1) - 1.0_dp)) < &
      1.0e-13_dp .and. all(ieee_is_finite(random_filter%particles)), &
      'curved nonlinear psi corrections')
   nonlinear_approximation = bssm_nonlinear_gaussian_approximation( &
      [ieee_value(0.0_dp, ieee_quiet_nan), 1.1_dp], [0.0_dp], &
      reshape([0.3_dp], [1, 1]), [0.7_dp, 0.05_dp, 0.2_dp], 1, &
      exponential_gaussian_observation, curved_transition_jacobian, 100, &
      1.0e-10_dp)
   call check(nonlinear_approximation%info == 0 .and. &
      nonlinear_approximation%converged .and. &
      nonlinear_approximation%scaling(1) == 0.0_dp, &
      'nonlinear approximation missing observation')

   ukf_fit = bssm_ukf([1.0_dp], [0.0_dp], &
      reshape([1.0_dp], [1, 1]), [1.0_dp, 0.0_dp, 1.0_dp], 1, &
      linear_gaussian_observation, nonlinear_transition_jacobian, &
      0.5_dp, 2.0_dp, 0.0_dp)
   call check(ukf_fit%info == 0 .and. &
      abs(ukf_fit%filtered_mean(1, 1) - 0.5_dp) < 1.0e-13_dp .and. &
      abs(ukf_fit%filtered_covariance(1, 1, 1) - 0.5_dp) < 1.0e-13_dp .and. &
      abs(ukf_fit%predicted_mean(1, 2) - 0.5_dp) < 1.0e-13_dp .and. &
      abs(ukf_fit%predicted_covariance(1, 1, 2) - 1.5_dp) < 1.0e-13_dp .and. &
      abs(ukf_fit%log_likelihood - expected) < 1.0e-12_dp, &
      'unscented Kalman linear equivalence')
   ukf_fit = bssm_ukf([1.2_dp, 1.1_dp], [0.0_dp], &
      reshape([0.5_dp], [1, 1]), [0.8_dp, 0.1_dp, 0.2_dp], 1, &
      exponential_gaussian_observation, nonlinear_transition_jacobian, &
      0.5_dp, 2.0_dp, 0.0_dp)
   call check(ukf_fit%info == 0 .and. &
      ieee_is_finite(ukf_fit%log_likelihood) .and. &
      all(ieee_is_finite(ukf_fit%filtered_mean)) .and. &
      all(ukf_fit%filtered_covariance >= 0.0_dp), &
      'nonlinear unscented Kalman filter')
   ukf_fit = bssm_ukf([ieee_value(0.0_dp, ieee_quiet_nan)], [0.3_dp], &
      reshape([0.2_dp], [1, 1]), [1.0_dp, 0.0_dp, 0.1_dp], 1, &
      exponential_gaussian_observation, nonlinear_transition_jacobian, &
      0.5_dp, 2.0_dp, 0.0_dp)
   call check(ukf_fit%info == 0 .and. &
      abs(ukf_fit%filtered_mean(1, 1) - 0.3_dp) < 1.0e-14_dp .and. &
      ukf_fit%log_likelihood == 0.0_dp, 'unscented Kalman missing observation')
   ukf_fit = bssm_ukf([1.0_dp], [0.0_dp], &
      reshape([1.0_dp], [1, 1]), [1.0_dp, 0.0_dp, 1.0_dp], 1, &
      linear_gaussian_observation, nonlinear_transition_jacobian, &
      -0.5_dp, 2.0_dp, 0.0_dp)
   call check(ukf_fit%info == 1, 'unscented Kalman tuning validation')
   multivariate_ekf_fit = bssm_multivariate_ukf( &
      reshape([1.0_dp, 3.0_dp, 2.0_dp, &
      ieee_value(0.0_dp, ieee_quiet_nan)], [2, 2]), [0.0_dp], &
      reshape([1.0_dp], [1, 1]), [1.0_dp, 0.0_dp, 1.0_dp], 2, 1, &
      multivariate_linear_gaussian_observation, &
      nonlinear_transition_jacobian, 0.5_dp, 2.0_dp, 0.0_dp)
   call check(multivariate_ekf_fit%info == 0 .and. &
      maxval(abs(multivariate_ekf_fit%filtered_mean(1, :) - &
      [4.0_dp/3.0_dp, 12.0_dp/7.0_dp])) < 1.0e-12_dp .and. &
      maxval(abs(multivariate_ekf_fit%filtered_covariance(1, 1, :) - &
      [1.0_dp/3.0_dp, 4.0_dp/7.0_dp])) < 1.0e-12_dp .and. &
      multivariate_ekf_fit%innovation(2, 2) == 0.0_dp .and. &
      all(multivariate_ekf_fit%innovation_covariance(2, :, 2) == 0.0_dp) .and. &
      all(multivariate_ekf_fit%iterations == 0) .and. &
      ieee_is_finite(multivariate_ekf_fit%log_likelihood), &
      'partially observed multivariate unscented Kalman filter')
   multivariate_ekf_fit = bssm_multivariate_ukf( &
      reshape([2.0_dp, 1.0_dp, 1.5_dp, &
      ieee_value(0.0_dp, ieee_quiet_nan)], [2, 2]), [0.0_dp], &
      reshape([0.5_dp], [1, 1]), [0.8_dp, 0.1_dp, 0.2_dp], 2, 1, &
      multivariate_nonlinear_gaussian_observation, &
      curved_transition_jacobian, 0.5_dp, 2.0_dp, 0.0_dp)
   call check(multivariate_ekf_fit%info == 0 .and. &
      ieee_is_finite(multivariate_ekf_fit%log_likelihood) .and. &
      all(ieee_is_finite(multivariate_ekf_fit%filtered_mean)) .and. &
      all(ieee_is_finite(multivariate_ekf_fit%filtered_covariance)), &
      'nonlinear multivariate unscented Kalman filter')
   multivariate_ekf_fit = bssm_multivariate_ukf( &
      reshape([1.0_dp, 2.0_dp], [2, 1]), [0.0_dp], &
      reshape([1.0_dp], [1, 1]), [1.0_dp, 0.0_dp, 1.0_dp], 2, 1, &
      multivariate_linear_gaussian_observation, &
      nonlinear_transition_jacobian, -0.5_dp, 2.0_dp, 0.0_dp)
   call check(multivariate_ekf_fit%info == 1, &
      'multivariate unscented Kalman tuning validation')

   call bssm_laplace_pseudo_observation(2.0_dp, 0.0_dp, bssm_poisson, &
      1.0_dp, 1.0_dp, pseudo_observation, pseudo_variance, approximation_info)
   call check(approximation_info == 0 .and. &
      abs(pseudo_observation - 1.0_dp) < 1.0e-14_dp .and. &
      abs(pseudo_variance - 1.0_dp) < 1.0e-14_dp, &
      'Poisson Laplace pseudo-observation')
   call bssm_laplace_pseudo_observation(2.0_dp, 0.0_dp, bssm_binomial, &
      1.0_dp, 4.0_dp, pseudo_observation, pseudo_variance, approximation_info)
   call check(approximation_info == 0 .and. &
      abs(pseudo_observation) < 1.0e-14_dp .and. &
      abs(pseudo_variance - 1.0_dp) < 1.0e-14_dp, &
      'binomial Laplace pseudo-observation')
   call bssm_laplace_pseudo_observation(1.0_dp, 0.0_dp, bssm_svm, &
      1.0_dp, 1.0_dp, pseudo_observation, pseudo_variance, approximation_info)
   call check(approximation_info == 0 .and. &
      abs(pseudo_observation) < 1.0e-14_dp .and. &
      abs(pseudo_variance - 2.0_dp) < 1.0e-14_dp, &
      'SVM Laplace pseudo-observation')
   call bssm_laplace_pseudo_observation(2.0_dp, 0.0_dp, &
      bssm_negative_binomial, 4.0_dp, 1.0_dp, pseudo_observation, &
      pseudo_variance, approximation_info)
   call check(approximation_info == 0 .and. &
      abs(pseudo_observation - 5.0_dp/6.0_dp) < 1.0e-14_dp .and. &
      abs(pseudo_variance - 25.0_dp/24.0_dp) < 1.0e-14_dp, &
      'negative-binomial Laplace pseudo-observation')
   call bssm_laplace_pseudo_observation(2.0_dp, 0.0_dp, bssm_gamma, &
      3.0_dp, 2.0_dp, pseudo_observation, pseudo_variance, approximation_info)
   call check(approximation_info == 0 .and. &
      abs(pseudo_observation) < 1.0e-14_dp .and. &
      abs(pseudo_variance - 1.0_dp/3.0_dp) < 1.0e-14_dp, &
      'Gamma Laplace pseudo-observation')

   approximation = bssm_gaussian_approximation([1.0_dp, &
      ieee_value(0.0_dp, ieee_quiet_nan), 4.0_dp, 3.0_dp], &
      reshape([1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], [1, 4]), &
      reshape([1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], [1, 1, 4]), &
      reshape([0.2_dp, 0.2_dp, 0.2_dp, 0.2_dp], [1, 1, 4]), &
      [log(2.0_dp)], reshape([0.5_dp], [1, 1]), bssm_poisson, 1.0_dp, &
      100, 1.0e-10_dp, offset=[0.1_dp, 0.1_dp, 0.1_dp, 0.1_dp], &
      state_offset=reshape([0.05_dp, 0.05_dp, 0.05_dp, 0.05_dp], [1, 4]))
   call check(approximation%info == 0 .and. approximation%converged .and. &
      approximation%iterations <= 100 .and. &
      all(approximation%observation_variance > 0.0_dp) .and. &
      maxval(abs(approximation%mode_signal - 0.1_dp - &
      approximation%mode_state(1, :))) < 1.0e-10_dp .and. &
      ieee_is_finite(approximation%gaussian_log_likelihood) .and. &
      ieee_is_finite(approximation%corrected_log_likelihood) .and. &
      all(ieee_is_finite(approximation%scaling)), &
      'Poisson Gaussian approximation')
   allocate(approximation_normals(1, 4, 4))
   approximation_normals = 0.0_dp
   simulation_smoother = bssm_approximation_simulation_draws(approximation, &
      approximation_normals)
   call check(simulation_smoother%info == 0 .and. &
      all(shape(simulation_smoother%trajectories) == [1, 4, 4]) .and. &
      maxval(abs(simulation_smoother%mean - &
      approximation%mode_state)) < 1.0e-13_dp .and. &
      maxval(abs(simulation_smoother%covariance)) < 1.0e-14_dp, &
      'supplied-draw non-Gaussian approximation smoother')
   call set_random_seed(1797)
   random_simulation_smoother = bssm_nongaussian_simulation_smoother( &
      [1.0_dp, ieee_value(0.0_dp, ieee_quiet_nan), 4.0_dp, 3.0_dp], &
      reshape([1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], [1, 4]), &
      reshape([1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], [1, 1, 4]), &
      reshape([0.2_dp, 0.2_dp, 0.2_dp, 0.2_dp], [1, 1, 4]), &
      [log(2.0_dp)], reshape([0.5_dp], [1, 1]), bssm_poisson, 1.0_dp, &
      4, 100, 1.0e-10_dp, .true., &
      offset=[0.1_dp, 0.1_dp, 0.1_dp, 0.1_dp], &
      state_offset=reshape([0.05_dp, 0.05_dp, 0.05_dp, 0.05_dp], [1, 4]))
   call check(random_simulation_smoother%info == 0 .and. &
      maxval(abs(random_simulation_smoother%mean - &
      approximation%mode_state)) < 1.0e-12_dp .and. &
      all(ieee_is_finite(random_simulation_smoother%trajectories)), &
      'antithetic non-Gaussian approximation smoother')

   multivariate_approximation = bssm_multivariate_gaussian_approximation( &
      reshape([2.0_dp, 5.0_dp, 0.2_dp, &
      ieee_value(0.0_dp, ieee_quiet_nan), 6.0_dp, 0.3_dp, 3.0_dp, &
      ieee_value(0.0_dp, ieee_quiet_nan), 0.1_dp], [3, 3]), &
      reshape([1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, &
      1.0_dp, 1.0_dp, 1.0_dp], [3, 1, 3]), &
      reshape([1.0_dp, 1.0_dp, 1.0_dp], [1, 1, 3]), &
      reshape([0.2_dp, 0.2_dp, 0.2_dp], [1, 1, 3]), [0.0_dp], &
      reshape([0.5_dp], [1, 1]), &
      [bssm_poisson, bssm_binomial, bssm_gaussian], &
      [1.0_dp, 1.0_dp, 0.5_dp], 100, 1.0e-10_dp, &
      auxiliary=reshape([1.0_dp, 10.0_dp, 1.0_dp, 1.0_dp, 10.0_dp, &
      1.0_dp, 1.0_dp, 10.0_dp, 1.0_dp], [3, 3]))
   call check(multivariate_approximation%info == 0 .and. &
      multivariate_approximation%converged .and. &
      all(shape(multivariate_approximation%mode_signal) == [3, 3]) .and. &
      maxval(abs(multivariate_approximation%mode_signal(1, :) - &
      multivariate_approximation%mode_state(1, :))) < 1.0e-10_dp .and. &
      maxval(abs(multivariate_approximation%mode_signal(2, :) - &
      multivariate_approximation%mode_state(1, :))) < 1.0e-10_dp .and. &
      maxval(abs(multivariate_approximation%observation_variance(3, :) - &
      0.25_dp)) < 1.0e-14_dp .and. &
      multivariate_approximation%scaling(1, 2) == 0.0_dp .and. &
      multivariate_approximation%scaling(2, 3) == 0.0_dp .and. &
      ieee_is_finite(multivariate_approximation%corrected_log_likelihood), &
      'mixed-family multivariate Gaussian approximation')
   deallocate(approximation_normals)
   allocate(approximation_normals(1, 4, 3))
   approximation_normals = 0.0_dp
   simulation_smoother = bssm_approximation_simulation_draws( &
      multivariate_approximation, approximation_normals)
   call check(simulation_smoother%info == 0 .and. &
      all(shape(simulation_smoother%trajectories) == [1, 3, 4]) .and. &
      maxval(abs(simulation_smoother%mean - &
      multivariate_approximation%mode_state)) < 1.0e-13_dp, &
      'supplied-draw multivariate approximation smoother')
   simulation_smoother = bssm_multivariate_simulation_smoother_draws( &
      reshape([2.0_dp, 5.0_dp, 0.2_dp, &
      ieee_value(0.0_dp, ieee_quiet_nan), 6.0_dp, 0.3_dp, 3.0_dp, &
      ieee_value(0.0_dp, ieee_quiet_nan), 0.1_dp], [3, 3]), &
      reshape([1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, &
      1.0_dp, 1.0_dp, 1.0_dp], [3, 1, 3]), &
      reshape([1.0_dp, 1.0_dp, 1.0_dp], [1, 1, 3]), &
      reshape([0.2_dp, 0.2_dp, 0.2_dp], [1, 1, 3]), [0.0_dp], &
      reshape([0.5_dp], [1, 1]), &
      [bssm_poisson, bssm_binomial, bssm_gaussian], &
      [1.0_dp, 1.0_dp, 0.5_dp], 100, 1.0e-10_dp, &
      approximation_normals, auxiliary=reshape([1.0_dp, 10.0_dp, &
      1.0_dp, 1.0_dp, 10.0_dp, 1.0_dp, 1.0_dp, 10.0_dp, 1.0_dp], &
      [3, 3]))
   call check(simulation_smoother%info == 0 .and. &
      maxval(abs(simulation_smoother%mean - &
      multivariate_approximation%mode_state)) < 1.0e-13_dp, &
      'mixed-family high-level approximation smoother')
   initial_normals = 0.0_dp
   innovations = 0.0_dp
   uniforms = 0.5_dp
   filter = bssm_multivariate_bootstrap_filter_draws(reshape([1.0_dp, &
      2.0_dp, 0.0_dp, ieee_value(0.0_dp, ieee_quiet_nan), 3.0_dp, &
      0.0_dp], [3, 2]), reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp], [3, 1, 2]), &
      reshape([1.0_dp, 1.0_dp], [1, 1, 2]), &
      reshape([0.0_dp, 0.0_dp], [1, 1, 2]), [0.0_dp], &
      reshape([0.0_dp], [1, 1]), &
      [bssm_poisson, bssm_binomial, bssm_gaussian], &
      [1.0_dp, 1.0_dp, 1.0_dp], initial_normals, innovations, uniforms, &
      auxiliary=reshape([1.0_dp, 4.0_dp, 1.0_dp, 1.0_dp, 4.0_dp, &
      1.0_dp], [3, 2]))
   expected = bssm_observation_log_density(1.0_dp, 0.0_dp, &
      bssm_poisson, 1.0_dp, 1.0_dp) + &
      bssm_observation_log_density(2.0_dp, 0.0_dp, bssm_binomial, &
      1.0_dp, 4.0_dp) + bssm_observation_log_density(0.0_dp, 0.0_dp, &
      bssm_gaussian, 1.0_dp, 1.0_dp) + &
      bssm_observation_log_density(3.0_dp, 0.0_dp, bssm_binomial, &
      1.0_dp, 4.0_dp) + bssm_observation_log_density(0.0_dp, 0.0_dp, &
      bssm_gaussian, 1.0_dp, 1.0_dp)
   call check(filter%info == 0 .and. &
      maxval(abs(filter%weights - 0.25_dp)) < 1.0e-14_dp .and. &
      abs(filter%log_likelihood - expected) < 1.0e-12_dp, &
      'mixed-family multivariate bootstrap filter')
   call set_random_seed(1801)
   random_filter = bssm_multivariate_bootstrap_filter(reshape([1.0_dp, &
      2.0_dp, 0.0_dp, 2.0_dp, 3.0_dp, 0.1_dp], [3, 2]), &
      reshape([1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], &
      [3, 1, 2]), reshape([1.0_dp, 1.0_dp], [1, 1, 2]), &
      reshape([0.2_dp, 0.2_dp], [1, 1, 2]), [0.0_dp], &
      reshape([0.5_dp], [1, 1]), &
      [bssm_poisson, bssm_binomial, bssm_gaussian], &
      [1.0_dp, 1.0_dp, 1.0_dp], 20, &
      auxiliary=reshape([1.0_dp, 4.0_dp, 1.0_dp, 1.0_dp, 4.0_dp, &
      1.0_dp], [3, 2]))
   call check(random_filter%info == 0 .and. &
      ieee_is_finite(random_filter%log_likelihood), &
      'random-stream multivariate bootstrap filter')
   multivariate_approximation = bssm_multivariate_gaussian_approximation( &
      reshape([2.0_dp, 5.0_dp, 0.2_dp, 3.0_dp, 6.0_dp, 0.3_dp], &
      [3, 2]), reshape([1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, &
      1.0_dp], [3, 1, 2]), reshape([1.0_dp, 1.0_dp], [1, 1, 2]), &
      reshape([0.2_dp, 0.2_dp], [1, 1, 2]), [0.0_dp], &
      reshape([0.0_dp], [1, 1]), &
      [bssm_poisson, bssm_binomial, bssm_gaussian], &
      [1.0_dp, 1.0_dp, 0.5_dp], 100, 1.0e-12_dp, &
      auxiliary=reshape([1.0_dp, 10.0_dp, 1.0_dp, 1.0_dp, 10.0_dp, &
      1.0_dp], [3, 2]))
   initial_normals = 0.0_dp
   innovations = 0.0_dp
   uniforms = 0.5_dp
   filter = bssm_multivariate_psi_filter_draws(reshape([2.0_dp, 5.0_dp, &
      0.2_dp, 3.0_dp, 6.0_dp, 0.3_dp], [3, 2]), &
      reshape([1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], &
      [3, 1, 2]), reshape([1.0_dp, 1.0_dp], [1, 1, 2]), &
      reshape([0.2_dp, 0.2_dp], [1, 1, 2]), [0.0_dp], &
      reshape([0.0_dp], [1, 1]), &
      [bssm_poisson, bssm_binomial, bssm_gaussian], &
      [1.0_dp, 1.0_dp, 0.5_dp], 100, 1.0e-12_dp, innovations, &
      initial_normals, uniforms, auxiliary=reshape([1.0_dp, 10.0_dp, &
      1.0_dp, 1.0_dp, 10.0_dp, 1.0_dp], [3, 2]))
   call check(filter%info == 0 .and. &
      maxval(abs(filter%weights - 0.25_dp)) < 1.0e-12_dp .and. &
      abs(filter%log_likelihood - &
      multivariate_approximation%corrected_log_likelihood) < 1.0e-10_dp, &
      'draw-driven multivariate psi filter')
   call set_random_seed(1803)
   random_filter = bssm_multivariate_psi_filter(reshape([2.0_dp, 5.0_dp, &
      0.2_dp, 3.0_dp, 6.0_dp, 0.3_dp], [3, 2]), &
      reshape([1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], &
      [3, 1, 2]), reshape([1.0_dp, 1.0_dp], [1, 1, 2]), &
      reshape([0.2_dp, 0.2_dp], [1, 1, 2]), [0.0_dp], &
      reshape([0.5_dp], [1, 1]), &
      [bssm_poisson, bssm_binomial, bssm_gaussian], &
      [1.0_dp, 1.0_dp, 0.5_dp], 30, 100, 1.0e-10_dp, &
      auxiliary=reshape([1.0_dp, 10.0_dp, 1.0_dp, 1.0_dp, 10.0_dp, &
      1.0_dp], [3, 2]))
   call check(random_filter%info == 0 .and. &
      ieee_is_finite(random_filter%log_likelihood) .and. &
      maxval(abs(sum(random_filter%weights, dim=1) - 1.0_dp)) < &
      1.0e-13_dp, 'random-stream multivariate psi filter')
   importance = bssm_multivariate_spdk_importance_draws(reshape([2.0_dp, &
      5.0_dp, 0.2_dp, 3.0_dp, 6.0_dp, 0.3_dp], [3, 2]), &
      reshape([1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], &
      [3, 1, 2]), reshape([1.0_dp, 1.0_dp], [1, 1, 2]), &
      reshape([0.2_dp, 0.2_dp], [1, 1, 2]), [0.0_dp], &
      reshape([0.0_dp], [1, 1]), &
      [bssm_poisson, bssm_binomial, bssm_gaussian], &
      [1.0_dp, 1.0_dp, 0.5_dp], 100, 1.0e-12_dp, innovations, &
      initial_normals, auxiliary=reshape([1.0_dp, 10.0_dp, 1.0_dp, &
      1.0_dp, 10.0_dp, 1.0_dp], [3, 2]))
   call check(importance%info == 0 .and. &
      all(shape(importance%trajectories) == [1, 3, 4]) .and. &
      maxval(abs(importance%weight - 0.25_dp)) < 1.0e-12_dp .and. &
      abs(importance%log_likelihood - &
      multivariate_approximation%corrected_log_likelihood) < 1.0e-10_dp .and. &
      abs(importance%effective_sample_size - 4.0_dp) < 1.0e-12_dp, &
      'draw-driven multivariate SPDK importance sampling')
   call set_random_seed(1805)
   random_importance = bssm_multivariate_spdk_importance(reshape([2.0_dp, &
      5.0_dp, 0.2_dp, 3.0_dp, 6.0_dp, 0.3_dp], [3, 2]), &
      reshape([1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], &
      [3, 1, 2]), reshape([1.0_dp, 1.0_dp], [1, 1, 2]), &
      reshape([0.2_dp, 0.2_dp], [1, 1, 2]), [0.0_dp], &
      reshape([0.5_dp], [1, 1]), &
      [bssm_poisson, bssm_binomial, bssm_gaussian], &
      [1.0_dp, 1.0_dp, 0.5_dp], 4, 100, 1.0e-10_dp, .true., &
      auxiliary=reshape([1.0_dp, 10.0_dp, 1.0_dp, 1.0_dp, 10.0_dp, &
      1.0_dp], [3, 2]))
   call check(random_importance%info == 0 .and. &
      abs(sum(random_importance%weight) - 1.0_dp) < 1.0e-13_dp .and. &
      ieee_is_finite(random_importance%log_likelihood) .and. &
      random_importance%effective_sample_size >= 1.0_dp .and. &
      random_importance%effective_sample_size <= 4.0_dp .and. &
      maxval(abs(random_importance%trajectories(:, 1:2, 1) + &
      random_importance%trajectories(:, 1:2, 3) - &
      random_importance%trajectories(:, 1:2, 2) - &
      random_importance%trajectories(:, 1:2, 4))) < 1.0e-12_dp, &
      'antithetic multivariate SPDK importance sampling')

   innovations = 0.0_dp
   initial_normals = 0.0_dp
   uniforms = 0.5_dp
   approximation = bssm_gaussian_approximation([1.0_dp, 3.0_dp], &
      reshape([1.0_dp, 1.0_dp], [1, 2]), &
      reshape([1.0_dp, 1.0_dp], [1, 1, 2]), &
      reshape([0.2_dp, 0.2_dp], [1, 1, 2]), [log(2.0_dp)], &
      reshape([0.0_dp], [1, 1]), bssm_poisson, 1.0_dp, 100, 1.0e-12_dp)
   filter = bssm_psi_filter_draws([1.0_dp, 3.0_dp], &
      reshape([1.0_dp, 1.0_dp], [1, 2]), &
      reshape([1.0_dp, 1.0_dp], [1, 1, 2]), &
      reshape([0.2_dp, 0.2_dp], [1, 1, 2]), [log(2.0_dp)], &
      reshape([0.0_dp], [1, 1]), bssm_poisson, 1.0_dp, 100, 1.0e-12_dp, &
      innovations, initial_normals, uniforms)
   call check(filter%info == 0 .and. &
      all(shape(filter%particles) == [1, 4, 3]) .and. &
      maxval(abs(filter%weights - 0.25_dp)) < 1.0e-13_dp .and. &
      abs(filter%log_likelihood - approximation%corrected_log_likelihood) < &
      1.0e-10_dp .and. &
      all(filter%ancestors(:, 1) == [1, 2, 3, 4]) .and. &
      all(filter%ancestors(:, 2) == [1, 2, 3, 4]), &
      'draw-driven psi auxiliary particle filter')
   call set_random_seed(925)
   random_filter = bssm_psi_filter([1.0_dp, 3.0_dp], &
      reshape([1.0_dp, 1.0_dp], [1, 2]), &
      reshape([1.0_dp, 1.0_dp], [1, 1, 2]), &
      reshape([0.2_dp, 0.2_dp], [1, 1, 2]), [log(2.0_dp)], &
      reshape([0.0_dp], [1, 1]), bssm_poisson, 1.0_dp, 20, 100, 1.0e-10_dp)
   call check(random_filter%info == 0 .and. &
      ieee_is_finite(random_filter%log_likelihood) .and. &
      all(ieee_is_finite(random_filter%filtered_mean)), &
      'random-stream psi auxiliary particle filter')

   importance = bssm_spdk_importance_draws([1.0_dp, 3.0_dp], &
      reshape([1.0_dp, 1.0_dp], [1, 2]), &
      reshape([1.0_dp, 1.0_dp], [1, 1, 2]), &
      reshape([0.2_dp, 0.2_dp], [1, 1, 2]), [log(2.0_dp)], &
      reshape([0.0_dp], [1, 1]), bssm_poisson, 1.0_dp, 100, 1.0e-12_dp, &
      innovations, initial_normals)
   call check(importance%info == 0 .and. &
      all(shape(importance%trajectories) == [1, 3, 4]) .and. &
      maxval(abs(importance%weight - 0.25_dp)) < 1.0e-13_dp .and. &
      abs(importance%log_likelihood - &
      approximation%corrected_log_likelihood) < 1.0e-10_dp .and. &
      abs(importance%effective_sample_size - 4.0_dp) < 1.0e-13_dp .and. &
      all(ieee_is_finite(importance%mean)) .and. &
      all(ieee_is_finite(importance%covariance)), &
      'draw-driven SPDK importance sampling')
   call set_random_seed(1147)
   random_importance = bssm_spdk_importance([1.0_dp, 3.0_dp], &
      reshape([1.0_dp, 1.0_dp], [1, 2]), &
      reshape([1.0_dp, 1.0_dp], [1, 1, 2]), &
      reshape([0.2_dp, 0.2_dp], [1, 1, 2]), [log(2.0_dp)], &
      reshape([0.5_dp], [1, 1]), bssm_poisson, 1.0_dp, 4, 100, &
      1.0e-10_dp, .true.)
   call check(random_importance%info == 0 .and. &
      abs(sum(random_importance%weight) - 1.0_dp) < 1.0e-13_dp .and. &
      ieee_is_finite(random_importance%log_likelihood) .and. &
      random_importance%effective_sample_size >= 1.0_dp .and. &
      random_importance%effective_sample_size <= 4.0_dp .and. &
      maxval(abs(random_importance%trajectories(:, 1:2, 1) + &
      random_importance%trajectories(:, 1:2, 3) - &
      random_importance%trajectories(:, 1:2, 2) - &
      random_importance%trajectories(:, 1:2, 4))) < 1.0e-12_dp, &
      'antithetic SPDK importance sampling')

   print '(a)', 'bssm tests passed'

contains

   pure real(dp) function pmmh_drift(state, parameters) result(value)
      ! Evaluate the PMMH test model's constant drift.
      real(dp), intent(in) :: state, parameters(:)

      value = parameters(1) + 0.0_dp*state
   end function pmmh_drift

   pure real(dp) function pmmh_diffusion(state, parameters) result(value)
      ! Evaluate the PMMH test model's zero diffusion.
      real(dp), intent(in) :: state, parameters(:)

      value = 0.0_dp*state + 0.0_dp*sum(parameters)
   end function pmmh_diffusion

   pure real(dp) function pmmh_prior(parameters) result(value)
      ! Evaluate a standard-normal log prior up to its constant.
      real(dp), intent(in) :: parameters(:)

      value = -0.5_dp*sum(parameters**2)
   end function pmmh_prior

   pure real(dp) function pmmh_approximate_likelihood(parameters) result(value)
      ! Evaluate the deterministic approximate likelihood used by MCMC tests.
      real(dp), intent(in) :: parameters(:)

      value = -0.5_dp*sum(parameters**2)
   end function pmmh_approximate_likelihood

   pure function pmmh_observation_log_density(time, observation, state, &
      parameters) result(log_density)
      ! Evaluate the PMMH test model's Gaussian observation density.
      integer, intent(in) :: time
      real(dp), intent(in) :: observation, state(:), parameters(:)
      real(dp) :: log_density

      log_density = bssm_observation_log_density(observation, state(1), &
         bssm_gaussian, 1.0_dp, 1.0_dp) + 0.0_dp*real(time, dp) + &
         0.0_dp*sum(parameters)
   end function pmmh_observation_log_density

   pure real(dp) function da_pmmh_drift(state, parameters) result(value)
      ! Evaluate multiplicative drift for delayed-acceptance tests.
      real(dp), intent(in) :: state, parameters(:)

      value = parameters(1)*state
   end function da_pmmh_drift

   pure subroutine pmmh_nonlinear_transition(time, state, parameters, mean, &
      noise_loading)
      ! Define deterministic parameter-dependent nonlinear PMMH dynamics.
      integer, intent(in) :: time
      real(dp), intent(in) :: state(:), parameters(:)
      real(dp), intent(out) :: mean(:), noise_loading(:, :)

      mean = parameters(1) + 0.0_dp*state + 0.0_dp*real(time, dp)
      noise_loading = 0.0_dp
   end subroutine pmmh_nonlinear_transition

   pure subroutine pmmh_deterministic_transition_jacobian(time, state, &
      parameters, mean, jacobian, noise_loading)
      ! Define deterministic dynamics and their Jacobian for DA-PMMH tests.
      integer, intent(in) :: time
      real(dp), intent(in) :: state(:), parameters(:)
      real(dp), intent(out) :: mean(:), jacobian(:, :)
      real(dp), intent(out) :: noise_loading(:, :)

      mean = parameters(1) + 0.0_dp*state + 0.0_dp*real(time, dp)
      jacobian = 0.0_dp
      noise_loading = 0.0_dp
   end subroutine pmmh_deterministic_transition_jacobian

   pure subroutine pmmh_psi_observation(time, state, parameters, mean, &
      jacobian, standard_deviation)
      ! Define the Gaussian observation used by psi-filter PMMH tests.
      integer, intent(in) :: time
      real(dp), intent(in) :: state(:), parameters(:)
      real(dp), intent(out) :: mean, jacobian(:), standard_deviation

      mean = state(1) + 0.0_dp*real(time, dp) + 0.0_dp*sum(parameters)
      jacobian = 1.0_dp
      standard_deviation = 1.0_dp
   end subroutine pmmh_psi_observation

   pure subroutine pmmh_psi_transition(time, state, parameters, mean, &
      jacobian, noise_loading)
      ! Define parameter-dependent dynamics for psi-filter PMMH tests.
      integer, intent(in) :: time
      real(dp), intent(in) :: state(:), parameters(:)
      real(dp), intent(out) :: mean(:), jacobian(:, :)
      real(dp), intent(out) :: noise_loading(:, :)

      mean = parameters(1) + 0.0_dp*state + 0.0_dp*real(time, dp)
      jacobian = 0.0_dp
      noise_loading = 0.1_dp
   end subroutine pmmh_psi_transition

   pure subroutine pmmh_ekpf_transition(time, state, parameters, mean, &
      noise_loading)
      ! Define parameter-dependent dynamics for EKF-proposal PMMH tests.
      integer, intent(in) :: time
      real(dp), intent(in) :: state(:), parameters(:)
      real(dp), intent(out) :: mean(:), noise_loading(:, :)

      mean = parameters(1) + 0.0_dp*state + 0.0_dp*real(time, dp)
      noise_loading = 0.1_dp
   end subroutine pmmh_ekpf_transition

   pure real(dp) function sde_drift(state, parameters) result(value)
      ! Evaluate the test Ornstein-Uhlenbeck drift.
      real(dp), intent(in) :: state, parameters(:)

      value = parameters(1)*(parameters(2) - state)
   end function sde_drift

   pure real(dp) function sde_diffusion(state, parameters) result(value)
      ! Evaluate the test constant diffusion coefficient.
      real(dp), intent(in) :: state, parameters(:)

      value = parameters(3) + 0.0_dp*state
   end function sde_diffusion

   pure real(dp) function sde_multiplicative_diffusion(state, parameters) &
      result(value)
      ! Evaluate a multiplicative diffusion coefficient.
      real(dp), intent(in) :: state, parameters(:)

      value = parameters(3)*state
   end function sde_multiplicative_diffusion

   pure real(dp) function sde_diffusion_derivative(state, parameters) &
      result(value)
      ! Evaluate the multiplicative diffusion derivative.
      real(dp), intent(in) :: state, parameters(:)

      value = parameters(3) + 0.0_dp*state
   end function sde_diffusion_derivative

   pure real(dp) function sde_zero_diffusion_derivative(state, parameters) &
      result(value)
      ! Evaluate the constant diffusion derivative.
      real(dp), intent(in) :: state, parameters(:)

      value = 0.0_dp*state + 0.0_dp*sum(parameters)
   end function sde_zero_diffusion_derivative

   pure function sde_observation_log_density(time, observation, state, &
      parameters) result(log_density)
      ! Evaluate the test SDE Gaussian observation density.
      integer, intent(in) :: time
      real(dp), intent(in) :: observation, state(:), parameters(:)
      real(dp) :: log_density

      log_density = bssm_observation_log_density(observation, state(1), &
         bssm_gaussian, parameters(4), 1.0_dp) + 0.0_dp*real(time, dp)
   end function sde_observation_log_density

   pure function nonlinear_log_density(time, observation, state, parameters) &
      result(log_density)
      ! Evaluate the test model's nonlinear Gaussian observation density.
      integer, intent(in) :: time
      real(dp), intent(in) :: observation, state(:), parameters(:)
      real(dp) :: log_density

      log_density = bssm_observation_log_density(observation, state(1)**2, &
         bssm_gaussian, 1.0_dp, 1.0_dp) + &
         0.0_dp*real(time, dp) + 0.0_dp*sum(parameters)
   end function nonlinear_log_density

   pure subroutine nonlinear_transition(time, state, parameters, mean, &
      noise_loading)
      ! Advance the test model through a nonlinear callback interface.
      integer, intent(in) :: time
      real(dp), intent(in) :: state(:), parameters(:)
      real(dp), intent(out) :: mean(:), noise_loading(:, :)

      mean(1) = parameters(1)*state(1) + parameters(2)
      noise_loading = parameters(3) + 0.0_dp*real(time, dp)
   end subroutine nonlinear_transition

   pure subroutine nonlinear_transition_jacobian(time, state, parameters, &
      mean, jacobian, noise_loading)
      ! Advance the test model and return its transition Jacobian.
      integer, intent(in) :: time
      real(dp), intent(in) :: state(:), parameters(:)
      real(dp), intent(out) :: mean(:), jacobian(:, :)
      real(dp), intent(out) :: noise_loading(:, :)

      mean(1) = parameters(1)*state(1) + parameters(2)
      jacobian = parameters(1)
      noise_loading = parameters(3) + 0.0_dp*real(time, dp)
   end subroutine nonlinear_transition_jacobian

   pure subroutine curved_transition_jacobian(time, state, parameters, mean, &
      jacobian, noise_loading)
      ! Advance the test model through a quadratic transition.
      integer, intent(in) :: time
      real(dp), intent(in) :: state(:), parameters(:)
      real(dp), intent(out) :: mean(:), jacobian(:, :)
      real(dp), intent(out) :: noise_loading(:, :)

      mean(1) = parameters(1)*state(1) + parameters(2)*state(1)**2
      jacobian(1, 1) = parameters(1) + 2.0_dp*parameters(2)*state(1)
      noise_loading = parameters(3) + 0.0_dp*real(time, dp)
   end subroutine curved_transition_jacobian

   pure subroutine nonlinear_prediction_transition(time, state, parameters, &
      mean, noise_loading)
      ! Advance the nonlinear prediction fixture for one forecast time.
      integer, intent(in) :: time
      real(dp), intent(in) :: state(:), parameters(:)
      real(dp), intent(out) :: mean(:), noise_loading(:, :)

      mean(1) = parameters(1)*state(1) + parameters(2) + &
         0.0_dp*real(time, dp)
      noise_loading = parameters(3)
   end subroutine nonlinear_prediction_transition

   pure subroutine nonlinear_prediction_observation(time, state, parameters, &
      mean, noise_loading)
      ! Evaluate a bivariate nonlinear Gaussian prediction fixture.
      integer, intent(in) :: time
      real(dp), intent(in) :: state(:), parameters(:)
      real(dp), intent(out) :: mean(:), noise_loading(:, :)

      mean(1) = state(1)**2 + parameters(1) + 0.0_dp*real(time, dp)
      mean(2) = 2.0_dp*state(1) + parameters(2)
      noise_loading = 0.0_dp
      noise_loading(1, 1) = 0.5_dp
      noise_loading(2, 1) = 0.25_dp
      noise_loading(2, 2) = 1.0_dp
   end subroutine nonlinear_prediction_observation

   pure subroutine linear_gaussian_observation(time, state, parameters, mean, &
      jacobian, standard_deviation)
      ! Define the test model's linear Gaussian observation equation.
      integer, intent(in) :: time
      real(dp), intent(in) :: state(:), parameters(:)
      real(dp), intent(out) :: mean, jacobian(:), standard_deviation

      mean = state(1) + 0.0_dp*real(time, dp)
      jacobian = 1.0_dp
      standard_deviation = parameters(3)
   end subroutine linear_gaussian_observation

   pure subroutine multivariate_linear_gaussian_observation(time, state, &
      parameters, mean, jacobian, noise_loading)
      ! Define a bivariate linear Gaussian observation equation.
      integer, intent(in) :: time
      real(dp), intent(in) :: state(:), parameters(:)
      real(dp), intent(out) :: mean(:), jacobian(:, :)
      real(dp), intent(out) :: noise_loading(:, :)
      integer :: component

      mean = state(1) + 0.0_dp*real(time, dp) + 0.0_dp*parameters(1)
      jacobian = 1.0_dp
      noise_loading = 0.0_dp
      do component = 1, min(size(noise_loading, 1), &
         size(noise_loading, 2))
         noise_loading(component, component) = 1.0_dp
      end do
   end subroutine multivariate_linear_gaussian_observation

   pure subroutine multivariate_nonlinear_gaussian_observation(time, state, &
      parameters, mean, jacobian, noise_loading)
      ! Define a nonlinear bivariate Gaussian observation equation.
      integer, intent(in) :: time
      real(dp), intent(in) :: state(:), parameters(:)
      real(dp), intent(out) :: mean(:), jacobian(:, :)
      real(dp), intent(out) :: noise_loading(:, :)

      mean = [exp(state(1)), 2.0_dp*state(1)] + &
         0.0_dp*real(time, dp)
      jacobian(:, 1) = [exp(state(1)), 2.0_dp]
      noise_loading = 0.0_dp
      noise_loading(1, 1) = parameters(3)
      noise_loading(2, 2) = parameters(3)
   end subroutine multivariate_nonlinear_gaussian_observation

   pure subroutine exponential_gaussian_observation(time, state, parameters, &
      mean, jacobian, standard_deviation)
      ! Define an exponential-mean Gaussian observation equation.
      integer, intent(in) :: time
      real(dp), intent(in) :: state(:), parameters(:)
      real(dp), intent(out) :: mean, jacobian(:), standard_deviation

      mean = exp(state(1)) + 0.0_dp*real(time, dp)
      jacobian = exp(state(1))
      standard_deviation = parameters(3)
   end subroutine exponential_gaussian_observation

   pure subroutine posterior_observation_model(parameters, &
      observation_loading, phi, offset, auxiliary, noise_loading, &
      correlated_gaussian, info)
      ! Expand one posterior draw into a mixed observation model.
      real(dp), intent(in) :: parameters(:)
      real(dp), intent(out) :: observation_loading(:, :, :)
      real(dp), intent(out) :: phi(:), offset(:, :), auxiliary(:, :)
      real(dp), intent(out) :: noise_loading(:, :, :)
      logical, intent(out) :: correlated_gaussian
      integer, intent(out) :: info

      observation_loading = parameters(1)
      phi = [1.0_dp, abs(parameters(2))]
      offset = parameters(3)
      auxiliary = 1.0_dp
      noise_loading = 0.0_dp
      correlated_gaussian = parameters(2) < 0.0_dp
      noise_loading(1, 1, :) = abs(parameters(2))
      noise_loading(2, 2, :) = abs(parameters(2))
      info = 0
   end subroutine posterior_observation_model

   subroutine check(condition, label)
      ! Stop at the first failed BSSM assertion.
      logical, intent(in) :: condition
      character(*), intent(in) :: label

      if (.not. condition) then
         print '(a,1x,a)', 'FAILED:', label
         error stop 1
      end if
   end subroutine check

end program test_bssm
