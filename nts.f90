! SPDX-License-Identifier: GPL-2.0-or-later
! SPDX-FileComment: Threshold autoregressive algorithms translated from R NTS.
module nts_mod
   !! Nonlinear time-series algorithms translated from the CRAN NTS package.
   use kind_mod, only: dp
   use stats_mod, only: ols_fit, quantile, sorted
   use random_mod, only: random_standard_normal_matrix, random_uniform, &
      random_standard_normal
   use linalg_mod, only: cholesky_lower, symmetric_eigen, invert_matrix, &
      inverse_logdet
   use special_functions_mod, only: regularized_beta, regularized_gamma_q
   use optimization_mod, only: optimization_result_t, bfgs_minimize_fd, &
      finite_difference_hessian
   use kfas_mod, only: ssm_model_t, kfs_filter_t, kfs_smoother_t, &
      kfs_filter, kfs_smooth
   use bssm_mod, only: bssm_stratified_resample
   use utils_mod, only: lowercase
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   private

   abstract interface
      pure subroutine nts_smc_step_t(previous, observation, parameters, &
         proposal_draws, proposed, log_weight_increment)
         !! Propagate particles and evaluate their incremental log weights.
         import dp
         real(dp), intent(in) :: previous(:, :) !! Previous state by particle.
         real(dp), intent(in) :: observation(:) !! Current observation vector.
         real(dp), intent(in) :: parameters(:) !! Model-specific parameter vector.
         real(dp), intent(in) :: proposal_draws(:, :) !! Supplied proposal random draws.
         real(dp), intent(out) :: proposed(:, :) !! Proposed state by particle.
         real(dp), intent(out) :: log_weight_increment(:) !! Incremental log weights.
      end subroutine nts_smc_step_t

      pure subroutine nts_rb_smc_step_t(previous, previous_mean, &
         previous_covariance, observation, parameters, proposal_draws, &
         proposed, proposed_mean, proposed_covariance, log_weight_increment)
         !! Propagate particles and their conditional Gaussian state distributions.
         import dp
         real(dp), intent(in) :: previous(:, :) !! Previous nonlinear state particles.
         real(dp), intent(in) :: previous_mean(:, :) !! Previous conditional Gaussian means.
         real(dp), intent(in) :: previous_covariance(:, :, :) !! Previous Gaussian covariances.
         real(dp), intent(in) :: observation(:) !! Current observation vector.
         real(dp), intent(in) :: parameters(:) !! Model-specific parameter vector.
         real(dp), intent(in) :: proposal_draws(:, :) !! Supplied proposal random draws.
         real(dp), intent(out) :: proposed(:, :) !! Proposed nonlinear state particles.
         real(dp), intent(out) :: proposed_mean(:, :) !! Updated conditional Gaussian means.
         real(dp), intent(out) :: proposed_covariance(:, :, :) !! Updated Gaussian covariances.
         real(dp), intent(out) :: log_weight_increment(:) !! Incremental log weights.
      end subroutine nts_rb_smc_step_t

      pure subroutine nts_smc_transition_log_density_t(previous, next, &
         parameters, log_density)
         !! Evaluate transition log densities between two particle clouds.
         import dp
         real(dp), intent(in) :: previous(:, :) !! Previous state by particle.
         real(dp), intent(in) :: next(:, :) !! Next state by particle.
         real(dp), intent(in) :: parameters(:) !! Model-specific parameter vector.
         real(dp), intent(out) :: log_density(:, :) !! Previous-by-next log densities.
      end subroutine nts_smc_transition_log_density_t
   end interface

   type, public :: nts_smc_filter_t
      !! Generic SMC particles, weights, ancestry, and delayed state estimates.
      real(dp), allocatable :: particles(:, :, :)
      real(dp), allocatable :: weights(:, :)
      integer, allocatable :: ancestors(:, :)
      real(dp), allocatable :: filtered_mean(:, :)
      real(dp), allocatable :: delayed_mean(:, :, :)
      real(dp), allocatable :: effective_sample_size(:)
      logical, allocatable :: resampled(:)
      real(dp) :: log_likelihood = 0.0_dp
      integer :: delay = 0
      integer :: info = 0
   end type nts_smc_filter_t

   type, public :: nts_smc_smoother_t
      !! Genealogical SMC trajectories and their smoothing moments.
      real(dp), allocatable :: trajectories(:, :, :)
      real(dp), allocatable :: terminal_weights(:)
      real(dp), allocatable :: mean(:, :)
      real(dp), allocatable :: covariance(:, :, :)
      integer :: info = 0
   end type nts_smc_smoother_t

   type, public :: nts_smc_marginal_smoother_t
      !! Marginal backward SMC weights and smoothing moments.
      real(dp), allocatable :: weights(:, :)
      real(dp), allocatable :: mean(:, :)
      real(dp), allocatable :: covariance(:, :, :)
      integer :: info = 0
   end type nts_smc_marginal_smoother_t

   type, public :: nts_rb_smc_filter_t
      !! Rao-Blackwellized SMC particles and conditional Gaussian moments.
      real(dp), allocatable :: particles(:, :, :)
      real(dp), allocatable :: gaussian_mean(:, :, :)
      real(dp), allocatable :: gaussian_covariance(:, :, :, :)
      real(dp), allocatable :: weights(:, :)
      integer, allocatable :: ancestors(:, :)
      real(dp), allocatable :: filtered_particle_mean(:, :)
      real(dp), allocatable :: filtered_gaussian_mean(:, :)
      real(dp), allocatable :: effective_sample_size(:)
      logical, allocatable :: resampled(:)
      real(dp) :: log_likelihood = 0.0_dp
      integer :: info = 0
   end type nts_rb_smc_filter_t

   type, public :: nts_tar_simulation_t
      !! Simulated SETAR observations, innovations, and generating parameters.
      real(dp), allocatable :: series(:)
      real(dp), allocatable :: innovations(:)
      real(dp), allocatable :: intercept(:)
      real(dp), allocatable :: trend_coefficient(:)
      real(dp), allocatable :: adf_coefficient(:)
      real(dp), allocatable :: ar(:, :)
      real(dp), allocatable :: thresholds(:)
      real(dp), allocatable :: innovation_sd(:)
      integer, allocatable :: ar_order(:)
      logical, allocatable :: lag_active(:, :)
      logical, allocatable :: include_trend(:)
      character(len=5) :: representation = 'level'
      integer :: forecast_step = 1
      integer :: delay = 1
      integer :: burnin = 0
      integer :: info = 0
   end type nts_tar_simulation_t

   type, public :: nts_tar_model_t
      !! Conditional least-squares estimates for a univariate SETAR model.
      real(dp), allocatable :: data(:)
      real(dp), allocatable :: coefficients(:, :)
      real(dp), allocatable :: trend_coefficient(:)
      real(dp), allocatable :: adf_coefficient(:)
      real(dp), allocatable :: thresholds(:)
      real(dp), allocatable :: innovation_sd(:)
      real(dp), allocatable :: residuals(:)
      real(dp), allocatable :: standardized_residuals(:)
      integer, allocatable :: ar_order(:)
      integer, allocatable :: lag_count(:)
      logical, allocatable :: lag_active(:, :)
      integer, allocatable :: regime_observations(:)
      logical, allocatable :: include_mean(:)
      logical, allocatable :: include_trend(:)
      character(len=5) :: representation = 'level'
      integer :: forecast_step = 1
      real(dp) :: aic = huge(1.0_dp)
      integer :: delay = 1
      integer :: info = 0
   end type nts_tar_model_t

   type, public :: nts_tar_search_t
      !! Threshold candidates, residual sums of squares, and selected SETAR fit.
      type(nts_tar_model_t) :: model
      real(dp), allocatable :: candidates(:)
      real(dp), allocatable :: rss(:)
      integer :: selected = 0
      integer :: info = 0
   end type nts_tar_search_t

   type, public :: nts_tar_forecast_t
      !! Simulation-based SETAR forecasts and pointwise prediction intervals.
      real(dp), allocatable :: mean(:)
      real(dp), allocatable :: lower(:)
      real(dp), allocatable :: upper(:)
      real(dp), allocatable :: simulations(:, :)
      real(dp) :: level = 0.95_dp
      integer :: origin = 0
      integer :: info = 0
   end type nts_tar_forecast_t

   type, public :: nts_mtar_simulation_t
      !! Simulated multivariate SETAR observations and generating parameters.
      real(dp), allocatable :: series(:, :)
      real(dp), allocatable :: innovations(:, :)
      real(dp), allocatable :: intercept(:, :)
      real(dp), allocatable :: trend(:, :)
      real(dp), allocatable :: ar(:, :, :, :)
      real(dp), allocatable :: covariance(:, :, :)
      real(dp), allocatable :: thresholds(:)
      integer, allocatable :: ar_order(:)
      integer, allocatable :: regime_observations(:)
      integer :: threshold_component = 1
      integer :: delay = 1
      integer :: burnin = 0
      integer :: info = 0
   end type nts_mtar_simulation_t

   type, public :: nts_mtar_model_t
      !! Conditional least-squares estimates for a multivariate SETAR model.
      real(dp), allocatable :: data(:, :)
      real(dp), allocatable :: intercept(:, :)
      real(dp), allocatable :: trend(:, :)
      real(dp), allocatable :: ar(:, :, :, :)
      real(dp), allocatable :: covariance(:, :, :)
      real(dp), allocatable :: thresholds(:)
      real(dp), allocatable :: residuals(:, :)
      real(dp), allocatable :: standardized_residuals(:, :)
      integer, allocatable :: ar_order(:)
      integer, allocatable :: regime_observations(:)
      logical, allocatable :: include_mean(:)
      logical, allocatable :: include_trend(:)
      real(dp) :: aic = huge(1.0_dp)
      integer :: threshold_component = 1
      integer :: delay = 1
      integer :: info = 0
   end type nts_mtar_model_t

   type, public :: nts_mtar_search_t
      !! Candidate scores and selected multivariate SETAR fit.
      type(nts_mtar_model_t) :: model
      real(dp), allocatable :: candidates(:)
      real(dp), allocatable :: score(:)
      integer :: selected = 0
      integer :: info = 0
   end type nts_mtar_search_t

   type, public :: nts_mtar_forecast_t
      !! Simulation-based multivariate SETAR forecasts and prediction intervals.
      real(dp), allocatable :: mean(:, :)
      real(dp), allocatable :: lower(:, :)
      real(dp), allocatable :: upper(:, :)
      real(dp), allocatable :: simulations(:, :, :)
      real(dp) :: level = 0.95_dp
      integer :: origin = 0
      integer :: info = 0
   end type nts_mtar_forecast_t

   type, public :: nts_msar_simulation_t
      !! Simulated Markov-switching AR observations, innovations, and states.
      real(dp), allocatable :: series(:)
      real(dp), allocatable :: innovations(:)
      integer, allocatable :: state(:)
      real(dp), allocatable :: intercept(:)
      real(dp), allocatable :: ar(:, :)
      integer, allocatable :: ar_order(:)
      real(dp), allocatable :: transition(:, :)
      real(dp), allocatable :: innovation_sd(:)
      integer :: burnin = 0
      integer :: info = 0
   end type nts_msar_simulation_t

   type, public :: nts_nonlinearity_test_t
      !! F test for threshold or quadratic autoregressive nonlinearity.
      real(dp) :: statistic = 0.0_dp
      real(dp) :: p_value = 1.0_dp
      integer :: numerator_df = 0
      integer :: denominator_df = 0
      integer :: initial_count = 0
      integer :: info = 0
   end type nts_nonlinearity_test_t

   type, public :: nts_tar_backtest_t
      !! Rolling-origin SETAR errors and overall and regime-specific accuracy.
      real(dp), allocatable :: error(:, :)
      integer, allocatable :: state(:)
      real(dp), allocatable :: rmse(:)
      real(dp), allocatable :: mae(:)
      real(dp), allocatable :: bias(:)
      real(dp), allocatable :: regime_rmse(:, :)
      real(dp), allocatable :: regime_mae(:, :)
      real(dp), allocatable :: regime_bias(:, :)
      integer, allocatable :: regime_count(:, :)
      integer :: origin = 0
      integer :: horizon = 0
      integer :: info = 0
   end type nts_tar_backtest_t

   type, public :: nts_rank_portmanteau_t
      !! Rank autocorrelations and cumulative rank-based portmanteau tests.
      real(dp), allocatable :: correlation(:)
      real(dp), allocatable :: expected_correlation(:)
      real(dp), allocatable :: correlation_variance(:)
      real(dp), allocatable :: statistic(:)
      real(dp), allocatable :: p_value(:)
      integer :: info = 0
   end type nts_rank_portmanteau_t

   type, public :: nts_prnd_test_t
      !! Peña-Rodríguez determinant portmanteau normal approximation.
      real(dp) :: statistic = 0.0_dp
      real(dp) :: p_value = 1.0_dp
      real(dp) :: determinant_statistic = 0.0_dp
      integer :: lag = 0
      integer :: info = 0
   end type nts_prnd_test_t

   type, public :: nts_tvar_fit_t
      !! Time-varying AR coefficients estimated by Gaussian state-space likelihood.
      real(dp), allocatable :: filtered_coefficients(:, :)
      real(dp), allocatable :: smoothed_coefficients(:, :)
      real(dp), allocatable :: filtered_covariance(:, :, :)
      real(dp), allocatable :: smoothed_covariance(:, :, :)
      real(dp), allocatable :: log_variance(:)
      real(dp), allocatable :: state_variance(:)
      integer, allocatable :: lags(:)
      real(dp) :: observation_variance = 0.0_dp
      real(dp) :: log_likelihood = -huge(1.0_dp)
      logical :: include_mean = .true.
      logical :: converged = .false.
      integer :: iterations = 0
      integer :: info = 0
   end type nts_tvar_fit_t

   type, public :: nts_rcar_fit_t
      !! Gaussian random-coefficient AR likelihood estimates and residuals.
      real(dp), allocatable :: coefficients(:)
      real(dp), allocatable :: coefficient_variance(:)
      real(dp), allocatable :: standard_error(:)
      real(dp), allocatable :: residuals(:)
      real(dp), allocatable :: standardized_residuals(:)
      integer, allocatable :: lags(:)
      real(dp) :: innovation_variance = 0.0_dp
      real(dp) :: sample_mean = 0.0_dp
      real(dp) :: log_likelihood = -huge(1.0_dp)
      logical :: include_mean = .true.
      logical :: converged = .false.
      integer :: iterations = 0
      integer :: info = 0
   end type nts_rcar_fit_t

   type, public :: nts_mtar_refinement_t
      !! Restricted multivariate TAR refit selected by coefficient t ratios.
      type(nts_mtar_model_t) :: model
      logical, allocatable :: retained(:, :, :)
      real(dp), allocatable :: standard_error(:, :, :)
      real(dp), allocatable :: t_ratio(:, :, :)
      real(dp) :: aic = huge(1.0_dp)
      real(dp) :: bic = huge(1.0_dp)
      real(dp) :: hq = huge(1.0_dp)
      real(dp) :: threshold = 1.0_dp
      integer :: info = 0
   end type nts_mtar_refinement_t

   type, public :: nts_acmx_fit_t
      !! Autoregressive conditional-mean count-model likelihood estimates.
      real(dp), allocatable :: exogenous_coefficients(:)
      real(dp), allocatable :: alpha(:)
      real(dp), allocatable :: gamma(:)
      real(dp), allocatable :: parameter(:)
      real(dp), allocatable :: standard_error(:)
      real(dp), allocatable :: conditional_mean(:)
      real(dp), allocatable :: residuals(:)
      real(dp), allocatable :: standardized_residuals(:)
      character(len=:), allocatable :: family
      real(dp) :: omega = 0.0_dp
      real(dp) :: dispersion = 1.0_dp
      real(dp) :: log_likelihood = -huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      real(dp) :: bic = huge(1.0_dp)
      integer :: ar_order = 0
      integer :: mean_order = 0
      logical :: converged = .false.
      integer :: iterations = 0
      integer :: info = 0
   end type nts_acmx_fit_t

   type, public :: nts_cfar_simulation_t
      !! Simulated continuous functional autoregression and its parameters.
      real(dp), allocatable :: series(:, :)
      real(dp), allocatable :: innovation(:)
      real(dp), allocatable :: kernel(:, :)
      real(dp) :: rho = 0.0_dp
      real(dp) :: sigma = 1.0_dp
      integer :: burnin = 0
      integer :: info = 0
   end type nts_cfar_simulation_t

   type, public :: nts_cfar_model_t
      !! Spline and OU generalized least-squares fit of a CFAR model.
      real(dp), allocatable :: kernel_coefficients(:, :)
      real(dp), allocatable :: kernel(:, :)
      real(dp), allocatable :: fitted(:, :)
      real(dp), allocatable :: residuals(:, :)
      real(dp) :: rho = 0.0_dp
      real(dp) :: sigma = 0.0_dp
      real(dp) :: log_likelihood = -huge(1.0_dp)
      integer :: order = 0
      integer :: spline_df = 0
      integer :: grid_size = 0
      integer :: info = 0
   end type nts_cfar_model_t

   type, public :: nts_cfar_order_tests_t
      !! Sequential F tests for the order of a continuous functional AR model.
      real(dp), allocatable :: statistic(:)
      real(dp), allocatable :: p_value(:)
      integer, allocatable :: numerator_df(:)
      integer, allocatable :: denominator_df(:)
      integer :: info = 0
   end type nts_cfar_order_tests_t

   type, public :: nts_cfar_forecast_t
      !! Point forecasts from a fitted continuous functional AR model.
      real(dp), allocatable :: mean(:, :)
      integer :: info = 0
   end type nts_cfar_forecast_t

   type, public :: nts_cfar_irregular_simulation_t
      !! Irregular observations sampled from latent complete CFAR curves.
      type(nts_cfar_simulation_t) :: complete
      real(dp), allocatable :: values(:, :)
      real(dp), allocatable :: positions(:, :)
      integer, allocatable :: counts(:)
      integer :: info = 0
   end type nts_cfar_irregular_simulation_t

   public :: nts_utar_simulate_from_innovations, nts_utar_simulate
   public :: nts_utar_estimate, nts_utar_threshold_search
   public :: nts_utar_forecast_draws, nts_utar_forecast
   public :: nts_mtar_simulate_from_standard, nts_mtar_simulate
   public :: nts_mtar_estimate, nts_mtar_threshold_search
   public :: nts_mtar_forecast_draws, nts_mtar_forecast
   public :: nts_msar_simulate_from_draws, nts_msar_simulate
   public :: nts_threshold_test, nts_tsay_test
   public :: nts_tar_backtest_draws, nts_tar_backtest
   public :: nts_rank_portmanteau, nts_quadratic_f_test, nts_prnd_test
   public :: nts_tvar_filter_smooth, nts_tvar_fit, nts_rcar_fit
   public :: nts_mtar_refine, nts_acmx_fit
   public :: nts_cfar_simulate_from_standard, nts_cfar_simulate
   public :: nts_cfar1_simulate, nts_cfar2_simulate, nts_cfar_estimate
   public :: nts_cfar_order_tests, nts_cfar_forecast, nts_cfar_partial_forecast
   public :: nts_cfar_irregular_estimate, nts_cfar_irregular_order_tests
   public :: nts_cfar_irregular_simulate, nts_cfar2_irregular_simulate
   public :: nts_smc_filter_draws, nts_smc_filter, nts_smc_smooth
   public :: nts_smc_marginal_smooth
   public :: nts_rb_smc_filter_draws, nts_rb_smc_filter

contains

   pure function nts_smc_filter_draws(step, observations, parameters, &
      initial_particles, proposal_draws, resampling_uniforms, &
      resample_schedule, delay) result(out)
      !! Run generic sequential Monte Carlo from explicitly supplied random draws.
      procedure(nts_smc_step_t) :: step !! Model-specific propagation and weighting callback.
      real(dp), intent(in) :: observations(:, :) !! Time-by-observation data matrix.
      real(dp), intent(in) :: parameters(:) !! Model-specific parameter vector.
      real(dp), intent(in) :: initial_particles(:, :) !! Initial state by particle.
      real(dp), intent(in) :: proposal_draws(:, :, :) !! Draw-by-particle-by-time proposals.
      real(dp), intent(in) :: resampling_uniforms(:, :) !! Particle-by-time uniforms.
      logical, intent(in) :: resample_schedule(:) !! Resampling decision at each time.
      integer, intent(in), optional :: delay !! Maximum delayed-estimate lag.
      type(nts_smc_filter_t) :: out
      real(dp), allocatable :: current(:, :), proposed(:, :), increment(:)
      real(dp), allocatable :: current_weight(:), log_unnormalized(:)
      integer, allocatable :: parent(:), selected(:)
      real(dp) :: log_normalizer
      integer :: states, particles, times, selected_delay
      integer :: time, state, particle, lag, target, index, trace_time, info

      states = size(initial_particles, 1)
      particles = size(initial_particles, 2)
      times = size(observations, 1)
      selected_delay = 0
      if (present(delay)) selected_delay = delay
      if (states < 1 .or. particles < 2 .or. times < 1 .or. &
         size(observations, 2) < 1 .or. size(proposal_draws, 2) /= particles .or. &
         size(proposal_draws, 3) /= times .or. &
         any(shape(resampling_uniforms) /= [particles, times]) .or. &
         size(resample_schedule) /= times .or. selected_delay < 0 .or. &
         .not. all(ieee_is_finite(observations)) .or. &
         .not. all(ieee_is_finite(parameters)) .or. &
         .not. all(ieee_is_finite(initial_particles)) .or. &
         .not. all(ieee_is_finite(proposal_draws)) .or. &
         .not. all(ieee_is_finite(resampling_uniforms)) .or. &
         any(resampling_uniforms < 0.0_dp .or. resampling_uniforms > 1.0_dp)) then
         out%info = 1
         return
      end if
      allocate(out%particles(states, particles, times))
      allocate(out%weights(particles, times), out%ancestors(particles, times))
      allocate(out%filtered_mean(states, times))
      allocate(out%delayed_mean(states, times, 0:selected_delay))
      allocate(out%effective_sample_size(times), out%resampled(times))
      allocate(current(states, particles), proposed(states, particles))
      allocate(increment(particles), current_weight(particles))
      allocate(log_unnormalized(particles), parent(particles))
      current = initial_particles
      current_weight = 1.0_dp/real(particles, dp)
      parent = [(particle, particle=1, particles)]
      out%delayed_mean = 0.0_dp
      out%ancestors(:, 1) = 0
      out%log_likelihood = 0.0_dp
      do time = 1, times
         call step(current, observations(time, :), parameters, &
            proposal_draws(:, :, time), proposed, increment)
         if (.not. all(ieee_is_finite(proposed)) .or. &
            .not. all(ieee_is_finite(increment))) then
            out%info = 2
            return
         end if
         log_unnormalized = log(max(current_weight, tiny(1.0_dp))) + increment
         call normalize_log_weights(log_unnormalized, out%weights(:, time), &
            log_normalizer, info)
         if (info /= 0) then
            out%info = 3
            return
         end if
         out%log_likelihood = out%log_likelihood + log_normalizer
         out%particles(:, :, time) = proposed
         if (time > 1) out%ancestors(:, time) = parent
         do state = 1, states
            out%filtered_mean(state, time) = dot_product( &
               proposed(state, :), out%weights(:, time))
         end do
         out%delayed_mean(:, time, 0) = out%filtered_mean(:, time)
         do lag = 1, min(selected_delay, time - 1)
            target = time - lag
            out%delayed_mean(:, target, lag) = 0.0_dp
            do particle = 1, particles
               index = particle
               do trace_time = time, target + 1, -1
                  index = out%ancestors(index, trace_time)
               end do
               out%delayed_mean(:, target, lag) = &
                  out%delayed_mean(:, target, lag) + out%weights(particle, time)* &
                  out%particles(:, index, target)
            end do
         end do
         out%effective_sample_size(time) = 1.0_dp/ &
            sum(out%weights(:, time)*out%weights(:, time))
         out%resampled(time) = resample_schedule(time)
         if (resample_schedule(time)) then
            selected = bssm_stratified_resample(out%weights(:, time), &
               resampling_uniforms(:, time))
            if (any(selected < 1) .or. any(selected > particles)) then
               out%info = 4
               return
            end if
            current = proposed(:, selected)
            current_weight = 1.0_dp/real(particles, dp)
            parent = selected
         else
            current = proposed
            current_weight = out%weights(:, time)
            parent = [(particle, particle=1, particles)]
         end if
      end do
      out%delay = selected_delay
   end function nts_smc_filter_draws

   function nts_smc_filter(step, observations, parameters, initial_particles, &
      draw_dimension, resample_schedule, delay) result(out)
      !! Run generic sequential Monte Carlo using the shared random stream.
      procedure(nts_smc_step_t) :: step !! Model-specific propagation and weighting callback.
      real(dp), intent(in) :: observations(:, :) !! Time-by-observation data matrix.
      real(dp), intent(in) :: parameters(:) !! Model-specific parameter vector.
      real(dp), intent(in) :: initial_particles(:, :) !! Initial state by particle.
      integer, intent(in) :: draw_dimension !! Proposal draws supplied per particle and time.
      logical, intent(in) :: resample_schedule(:) !! Resampling decision at each time.
      integer, intent(in), optional :: delay !! Maximum delayed-estimate lag.
      type(nts_smc_filter_t) :: out
      real(dp), allocatable :: proposal_draws(:, :, :), resampling_uniforms(:, :)
      integer :: time, particles, times

      particles = size(initial_particles, 2)
      times = size(observations, 1)
      if (draw_dimension < 1 .or. particles < 2 .or. times < 1) then
         out%info = 1
         return
      end if
      allocate(proposal_draws(draw_dimension, particles, times))
      allocate(resampling_uniforms(particles, times))
      do time = 1, times
         call random_standard_normal_matrix(proposal_draws(:, :, time))
      end do
      call random_number(resampling_uniforms)
      out = nts_smc_filter_draws(step, observations, parameters, &
         initial_particles, proposal_draws, resampling_uniforms, &
         resample_schedule, delay)
   end function nts_smc_filter

   pure function nts_smc_smooth(filter) result(out)
      !! Trace all terminal particles backward through a generic SMC genealogy.
      type(nts_smc_filter_t), intent(in) :: filter !! Completed generic SMC filter.
      type(nts_smc_smoother_t) :: out
      real(dp), allocatable :: difference(:)
      integer :: states, particles, times, terminal, time, index

      if (filter%info /= 0 .or. .not. allocated(filter%particles) .or. &
         .not. allocated(filter%weights) .or. .not. allocated(filter%ancestors)) then
         out%info = 1
         return
      end if
      states = size(filter%particles, 1)
      particles = size(filter%particles, 2)
      times = size(filter%particles, 3)
      allocate(out%trajectories(states, times, particles))
      allocate(out%terminal_weights(particles), out%mean(states, times))
      allocate(out%covariance(states, states, times), difference(states))
      out%terminal_weights = filter%weights(:, times)
      do terminal = 1, particles
         index = terminal
         out%trajectories(:, times, terminal) = &
            filter%particles(:, index, times)
         do time = times, 2, -1
            index = filter%ancestors(index, time)
            out%trajectories(:, time - 1, terminal) = &
               filter%particles(:, index, time - 1)
         end do
      end do
      do time = 1, times
         out%mean(:, time) = 0.0_dp
         do terminal = 1, particles
            out%mean(:, time) = out%mean(:, time) + &
               out%terminal_weights(terminal)*out%trajectories(:, time, terminal)
         end do
         out%covariance(:, :, time) = 0.0_dp
         do terminal = 1, particles
            difference = out%trajectories(:, time, terminal) - out%mean(:, time)
            out%covariance(:, :, time) = out%covariance(:, :, time) + &
               out%terminal_weights(terminal)*outer_product(difference, difference)
         end do
      end do
   end function nts_smc_smooth

   pure function nts_smc_marginal_smooth(filter, parameters, &
      transition_log_density) result(out)
      !! Smooth SMC particles with the marginal backward-weight recursion.
      type(nts_smc_filter_t), intent(in) :: filter !! Completed generic SMC filter.
      real(dp), intent(in) :: parameters(:) !! Model-specific transition parameters.
      procedure(nts_smc_transition_log_density_t) :: transition_log_density !! Transition-density callback.
      type(nts_smc_marginal_smoother_t) :: out
      real(dp), allocatable :: log_density(:, :), scaled_density(:, :)
      real(dp), allocatable :: denominator(:), difference(:)
      real(dp) :: maximum_density, total
      integer :: states, particles, times, time, previous, next

      if (filter%info /= 0 .or. .not. allocated(filter%particles) .or. &
         .not. allocated(filter%weights) .or. &
         .not. all(ieee_is_finite(parameters))) then
         out%info = 1
         return
      end if
      states = size(filter%particles, 1)
      particles = size(filter%particles, 2)
      times = size(filter%particles, 3)
      allocate(out%weights(particles, times), out%mean(states, times))
      allocate(out%covariance(states, states, times))
      allocate(log_density(particles, particles))
      allocate(scaled_density(particles, particles), denominator(particles))
      allocate(difference(states))
      out%weights(:, times) = filter%weights(:, times)
      do time = times - 1, 1, -1
         call transition_log_density(filter%particles(:, :, time), &
            filter%particles(:, :, time + 1), parameters, log_density)
         if (.not. all(ieee_is_finite(log_density))) then
            out%info = 2
            return
         end if
         do next = 1, particles
            maximum_density = maxval(log_density(:, next))
            scaled_density(:, next) = exp(log_density(:, next) - maximum_density)
            denominator(next) = dot_product(filter%weights(:, time), &
               scaled_density(:, next))
         end do
         if (any(denominator <= tiny(1.0_dp))) then
            out%info = 3
            return
         end if
         do previous = 1, particles
            out%weights(previous, time) = filter%weights(previous, time)*sum( &
               out%weights(:, time + 1)*scaled_density(previous, :)/denominator)
         end do
         total = sum(out%weights(:, time))
         if (total <= 0.0_dp .or. .not. ieee_is_finite(total)) then
            out%info = 3
            return
         end if
         out%weights(:, time) = out%weights(:, time)/total
      end do
      do time = 1, times
         out%mean(:, time) = 0.0_dp
         do previous = 1, particles
            out%mean(:, time) = out%mean(:, time) + out%weights(previous, time)* &
               filter%particles(:, previous, time)
         end do
         out%covariance(:, :, time) = 0.0_dp
         do previous = 1, particles
            difference = filter%particles(:, previous, time) - out%mean(:, time)
            out%covariance(:, :, time) = out%covariance(:, :, time) + &
               out%weights(previous, time)*outer_product(difference, difference)
         end do
      end do
   end function nts_smc_marginal_smooth

   pure function nts_rb_smc_filter_draws(step, observations, parameters, &
      initial_particles, initial_mean, initial_covariance, proposal_draws, &
      resampling_uniforms, resample_schedule) result(out)
      !! Run Rao-Blackwellized SMC from explicitly supplied random draws.
      procedure(nts_rb_smc_step_t) :: step !! Rao-Blackwellized propagation callback.
      real(dp), intent(in) :: observations(:, :) !! Time-by-observation data matrix.
      real(dp), intent(in) :: parameters(:) !! Model-specific parameter vector.
      real(dp), intent(in) :: initial_particles(:, :) !! Initial nonlinear particles.
      real(dp), intent(in) :: initial_mean(:, :) !! Initial conditional Gaussian means.
      real(dp), intent(in) :: initial_covariance(:, :, :) !! Initial Gaussian covariances.
      real(dp), intent(in) :: proposal_draws(:, :, :) !! Draw-by-particle-by-time proposals.
      real(dp), intent(in) :: resampling_uniforms(:, :) !! Particle-by-time uniforms.
      logical, intent(in) :: resample_schedule(:) !! Resampling decision at each time.
      type(nts_rb_smc_filter_t) :: out
      real(dp), allocatable :: current(:, :), proposed(:, :)
      real(dp), allocatable :: current_mean(:, :), proposed_mean(:, :)
      real(dp), allocatable :: current_covariance(:, :, :)
      real(dp), allocatable :: proposed_covariance(:, :, :)
      real(dp), allocatable :: current_weight(:), increment(:), log_weight(:)
      integer, allocatable :: parent(:), selected(:)
      real(dp) :: log_normalizer
      integer :: states, gaussian_states, particles, times
      integer :: time, state, particle, info

      states = size(initial_particles, 1)
      particles = size(initial_particles, 2)
      gaussian_states = size(initial_mean, 1)
      times = size(observations, 1)
      if (states < 1 .or. gaussian_states < 1 .or. particles < 2 .or. times < 1 .or. &
         size(initial_mean, 2) /= particles .or. &
         any(shape(initial_covariance) /= [gaussian_states, gaussian_states, particles]) .or. &
         size(proposal_draws, 2) /= particles .or. size(proposal_draws, 3) /= times .or. &
         any(shape(resampling_uniforms) /= [particles, times]) .or. &
         size(resample_schedule) /= times .or. &
         .not. all(ieee_is_finite(observations)) .or. &
         .not. all(ieee_is_finite(initial_particles)) .or. &
         .not. all(ieee_is_finite(initial_mean)) .or. &
         .not. all(ieee_is_finite(initial_covariance))) then
         out%info = 1
         return
      end if
      allocate(out%particles(states, particles, times))
      allocate(out%gaussian_mean(gaussian_states, particles, times))
      allocate(out%gaussian_covariance(gaussian_states, gaussian_states, &
         particles, times))
      allocate(out%weights(particles, times), out%ancestors(particles, times))
      allocate(out%filtered_particle_mean(states, times))
      allocate(out%filtered_gaussian_mean(gaussian_states, times))
      allocate(out%effective_sample_size(times), out%resampled(times))
      allocate(current(states, particles), proposed(states, particles))
      allocate(current_mean(gaussian_states, particles))
      allocate(proposed_mean(gaussian_states, particles))
      allocate(current_covariance(gaussian_states, gaussian_states, particles))
      allocate(proposed_covariance(gaussian_states, gaussian_states, particles))
      allocate(current_weight(particles), increment(particles), log_weight(particles))
      allocate(parent(particles))
      current = initial_particles
      current_mean = initial_mean
      current_covariance = initial_covariance
      current_weight = 1.0_dp/real(particles, dp)
      parent = [(particle, particle=1, particles)]
      out%ancestors(:, 1) = 0
      do time = 1, times
         call step(current, current_mean, current_covariance, &
            observations(time, :), parameters, proposal_draws(:, :, time), &
            proposed, proposed_mean, proposed_covariance, increment)
         if (.not. all(ieee_is_finite(proposed)) .or. &
            .not. all(ieee_is_finite(proposed_mean)) .or. &
            .not. all(ieee_is_finite(proposed_covariance)) .or. &
            .not. all(ieee_is_finite(increment))) then
            out%info = 2
            return
         end if
         log_weight = log(max(current_weight, tiny(1.0_dp))) + increment
         call normalize_log_weights(log_weight, out%weights(:, time), &
            log_normalizer, info)
         if (info /= 0) then
            out%info = 3
            return
         end if
         out%log_likelihood = out%log_likelihood + log_normalizer
         out%particles(:, :, time) = proposed
         out%gaussian_mean(:, :, time) = proposed_mean
         out%gaussian_covariance(:, :, :, time) = proposed_covariance
         if (time > 1) out%ancestors(:, time) = parent
         do state = 1, states
            out%filtered_particle_mean(state, time) = dot_product( &
               proposed(state, :), out%weights(:, time))
         end do
         do state = 1, gaussian_states
            out%filtered_gaussian_mean(state, time) = dot_product( &
               proposed_mean(state, :), out%weights(:, time))
         end do
         out%effective_sample_size(time) = 1.0_dp/ &
            sum(out%weights(:, time)*out%weights(:, time))
         out%resampled(time) = resample_schedule(time)
         if (resample_schedule(time)) then
            selected = bssm_stratified_resample(out%weights(:, time), &
               resampling_uniforms(:, time))
            current = proposed(:, selected)
            current_mean = proposed_mean(:, selected)
            current_covariance = proposed_covariance(:, :, selected)
            current_weight = 1.0_dp/real(particles, dp)
            parent = selected
         else
            current = proposed
            current_mean = proposed_mean
            current_covariance = proposed_covariance
            current_weight = out%weights(:, time)
            parent = [(particle, particle=1, particles)]
         end if
      end do
   end function nts_rb_smc_filter_draws

   function nts_rb_smc_filter(step, observations, parameters, initial_particles, &
      initial_mean, initial_covariance, draw_dimension, resample_schedule) result(out)
      !! Run Rao-Blackwellized SMC using the shared random stream.
      procedure(nts_rb_smc_step_t) :: step !! Rao-Blackwellized propagation callback.
      real(dp), intent(in) :: observations(:, :) !! Time-by-observation data matrix.
      real(dp), intent(in) :: parameters(:) !! Model-specific parameter vector.
      real(dp), intent(in) :: initial_particles(:, :) !! Initial nonlinear particles.
      real(dp), intent(in) :: initial_mean(:, :) !! Initial conditional Gaussian means.
      real(dp), intent(in) :: initial_covariance(:, :, :) !! Initial Gaussian covariances.
      integer, intent(in) :: draw_dimension !! Proposal draws supplied per particle and time.
      logical, intent(in) :: resample_schedule(:) !! Resampling decision at each time.
      type(nts_rb_smc_filter_t) :: out
      real(dp), allocatable :: proposal_draws(:, :, :), resampling_uniforms(:, :)
      integer :: time, particles, times

      particles = size(initial_particles, 2)
      times = size(observations, 1)
      if (draw_dimension < 1 .or. particles < 2 .or. times < 1) then
         out%info = 1
         return
      end if
      allocate(proposal_draws(draw_dimension, particles, times))
      allocate(resampling_uniforms(particles, times))
      do time = 1, times
         call random_standard_normal_matrix(proposal_draws(:, :, time))
      end do
      call random_number(resampling_uniforms)
      out = nts_rb_smc_filter_draws(step, observations, parameters, &
         initial_particles, initial_mean, initial_covariance, proposal_draws, &
         resampling_uniforms, resample_schedule)
   end function nts_rb_smc_filter

   pure function nts_cfar_simulate_from_standard(kernel, rho, sigma, &
      standard_draws, burnin) result(out)
      !! Simulate a CFAR process from supplied independent standard-normal draws.
      real(dp), intent(in) :: kernel(:, :) !! Kernel grid by autoregressive lag.
      real(dp), intent(in) :: rho !! Positive spatial OU decay parameter.
      real(dp), intent(in) :: sigma !! Positive innovation scale.
      real(dp), intent(in) :: standard_draws(:, :) !! Time-by-grid normal draws.
      integer, intent(in), optional :: burnin !! Number of leading curves to discard.
      type(nts_cfar_simulation_t) :: out
      real(dp), allocatable :: complete(:, :), innovation(:, :)
      real(dp) :: adjacent, innovation_scale
      integer :: points, order, total, selected_burnin
      integer :: time, location, lag, source, kernel_index

      points = size(standard_draws, 2)
      order = size(kernel, 2)
      total = size(standard_draws, 1)
      selected_burnin = 0
      if (present(burnin)) selected_burnin = burnin
      if (points < 2 .or. order < 1 .or. size(kernel, 1) /= 2*points - 1 .or. &
         total <= order .or. selected_burnin < 0 .or. selected_burnin >= total .or. &
         rho <= 0.0_dp .or. sigma <= 0.0_dp .or. &
         .not. all(ieee_is_finite(kernel)) .or. &
         .not. all(ieee_is_finite(standard_draws))) then
         out%info = 1
         return
      end if
      allocate(complete(total, points), innovation(total, points))
      adjacent = exp(-rho/real(points - 1, dp))
      innovation_scale = sqrt(max(0.0_dp, 1.0_dp - adjacent*adjacent))
      do time = 1, total
         innovation(time, 1) = standard_draws(time, 1)
         do location = 2, points
            innovation(time, location) = adjacent*innovation(time, location - 1) + &
               innovation_scale*standard_draws(time, location)
         end do
      end do
      complete = innovation
      do time = order + 1, total
         complete(time, :) = innovation(time, :)
         do lag = 1, order
            do location = 1, points
               do source = 1, points
                  kernel_index = points + location - source
                  complete(time, location) = complete(time, location) + &
                     kernel(kernel_index, lag)*complete(time - lag, source)/ &
                     real(points, dp)
               end do
            end do
         end do
      end do
      out%series = sigma*complete(selected_burnin + 1:, :)
      out%innovation = sigma*innovation(total, :)
      out%kernel = kernel
      out%rho = rho
      out%sigma = sigma
      out%burnin = selected_burnin
   end function nts_cfar_simulate_from_standard

   function nts_cfar_simulate(kernel, rho, sigma, observations, burnin) result(out)
      !! Simulate a CFAR process using the shared random-number generator.
      real(dp), intent(in) :: kernel(:, :) !! Kernel grid by autoregressive lag.
      real(dp), intent(in) :: rho !! Positive spatial OU decay parameter.
      real(dp), intent(in) :: sigma !! Positive innovation scale.
      integer, intent(in) :: observations !! Number of retained functional observations.
      integer, intent(in), optional :: burnin !! Number of leading curves to discard.
      type(nts_cfar_simulation_t) :: out
      real(dp), allocatable :: draws(:, :)
      integer :: selected_burnin, total, points

      selected_burnin = 100
      if (present(burnin)) selected_burnin = burnin
      points = (size(kernel, 1) + 1)/2
      if (observations < 1 .or. selected_burnin < 0 .or. points < 2) then
         out%info = 1
         return
      end if
      total = observations + selected_burnin
      allocate(draws(total, points))
      call random_standard_normal_matrix(draws)
      out = nts_cfar_simulate_from_standard(kernel, rho, sigma, draws, &
         selected_burnin)
   end function nts_cfar_simulate

   function nts_cfar1_simulate(kernel, rho, sigma, observations, burnin) result(out)
      !! Simulate a first-order CFAR process using one sampled kernel.
      real(dp), intent(in) :: kernel(:) !! Kernel values on an equally spaced grid.
      real(dp), intent(in) :: rho !! Positive spatial OU decay parameter.
      real(dp), intent(in) :: sigma !! Positive innovation scale.
      integer, intent(in) :: observations !! Number of retained functional observations.
      integer, intent(in), optional :: burnin !! Number of leading curves to discard.
      type(nts_cfar_simulation_t) :: out
      real(dp), allocatable :: kernels(:, :)

      allocate(kernels(size(kernel), 1))
      kernels(:, 1) = kernel
      out = nts_cfar_simulate(kernels, rho, sigma, observations, burnin)
   end function nts_cfar1_simulate

   function nts_cfar2_simulate(kernel1, kernel2, rho, sigma, observations, &
      burnin) result(out)
      !! Simulate a second-order CFAR process using two sampled kernels.
      real(dp), intent(in) :: kernel1(:) !! First-lag kernel on a regular grid.
      real(dp), intent(in) :: kernel2(:) !! Second-lag kernel on the same grid.
      real(dp), intent(in) :: rho !! Positive spatial OU decay parameter.
      real(dp), intent(in) :: sigma !! Positive innovation scale.
      integer, intent(in) :: observations !! Number of retained functional observations.
      integer, intent(in), optional :: burnin !! Number of leading curves to discard.
      type(nts_cfar_simulation_t) :: out
      real(dp), allocatable :: kernels(:, :)

      if (size(kernel1) /= size(kernel2)) then
         out%info = 1
         return
      end if
      allocate(kernels(size(kernel1), 2))
      kernels(:, 1) = kernel1
      kernels(:, 2) = kernel2
      out = nts_cfar_simulate(kernels, rho, sigma, observations, burnin)
   end function nts_cfar2_simulate

   function nts_cfar_irregular_simulate(kernel, rho, sigma, observations, &
      minimum_observations, mean_extra_observations, burnin) result(out)
      !! Simulate CFAR curves and sample heterogeneous observation locations.
      real(dp), intent(in) :: kernel(:, :) !! Kernel grid by autoregressive lag.
      real(dp), intent(in) :: rho !! Positive spatial OU decay parameter.
      real(dp), intent(in) :: sigma !! Positive innovation scale.
      integer, intent(in) :: observations !! Number of retained curves.
      integer, intent(in) :: minimum_observations !! Minimum locations per curve.
      real(dp), intent(in) :: mean_extra_observations !! Poisson mean above the minimum.
      integer, intent(in), optional :: burnin !! Number of leading curves to discard.
      type(nts_cfar_irregular_simulation_t) :: out
      real(dp), allocatable :: sampled(:)
      integer :: points, curve, location, count, maximum_count

      points = (size(kernel, 1) + 1)/2
      if (observations < 1 .or. minimum_observations < 2 .or. &
         minimum_observations > points .or. mean_extra_observations < 0.0_dp) then
         out%info = 1
         return
      end if
      out%complete = nts_cfar_simulate(kernel, rho, sigma, observations, burnin)
      if (out%complete%info /= 0) then
         out%info = out%complete%info
         return
      end if
      allocate(out%counts(observations))
      do curve = 1, observations
         out%counts(curve) = min(points, minimum_observations + &
            random_poisson(mean_extra_observations))
      end do
      maximum_count = maxval(out%counts)
      allocate(out%values(observations, maximum_count))
      allocate(out%positions(observations, maximum_count))
      out%values = 0.0_dp
      out%positions = 0.0_dp
      do curve = 1, observations
         count = out%counts(curve)
         allocate(sampled(count))
         do location = 1, count
            sampled(location) = random_uniform()
         end do
         call sort_small(sampled)
         out%positions(curve, :count) = sampled
         do location = 1, count
            out%values(curve, location) = interpolate_regular_curve( &
               out%complete%series(curve, :), sampled(location))
         end do
         deallocate(sampled)
      end do
   end function nts_cfar_irregular_simulate

   function nts_cfar2_irregular_simulate(kernel1, kernel2, rho, sigma, &
      observations, minimum_observations, mean_extra_observations, burnin) &
      result(out)
      !! Simulate an irregularly observed second-order CFAR process.
      real(dp), intent(in) :: kernel1(:) !! First-lag kernel on a regular grid.
      real(dp), intent(in) :: kernel2(:) !! Second-lag kernel on the same grid.
      real(dp), intent(in) :: rho !! Positive spatial OU decay parameter.
      real(dp), intent(in) :: sigma !! Positive innovation scale.
      integer, intent(in) :: observations !! Number of retained curves.
      integer, intent(in) :: minimum_observations !! Minimum locations per curve.
      real(dp), intent(in) :: mean_extra_observations !! Poisson mean above the minimum.
      integer, intent(in), optional :: burnin !! Number of leading curves to discard.
      type(nts_cfar_irregular_simulation_t) :: out
      real(dp), allocatable :: kernels(:, :)

      if (size(kernel1) /= size(kernel2)) then
         out%info = 1
         return
      end if
      allocate(kernels(size(kernel1), 2))
      kernels(:, 1) = kernel1
      kernels(:, 2) = kernel2
      out = nts_cfar_irregular_simulate(kernels, rho, sigma, observations, &
         minimum_observations, mean_extra_observations, burnin)
   end function nts_cfar2_irregular_simulate

   pure function nts_cfar_estimate(series, order, spline_df) result(model)
      !! Estimate a CFAR model by spline convolution and OU-weighted GLS.
      real(dp), intent(in) :: series(:, :) !! Time-by-grid functional observations.
      integer, intent(in) :: order !! Positive functional autoregressive order.
      integer, intent(in), optional :: spline_df !! Spline degrees of freedom excluding intercept.
      type(nts_cfar_model_t) :: model
      real(dp), allocatable :: design(:, :), response(:), coefficient(:)
      real(dp), allocatable :: fitted_vector(:), residual_vector(:)
      real(dp) :: rho, score, best_score, lower, upper, step
      integer :: selected_df, points, observations, iteration, candidate, info

      selected_df = 4
      if (present(spline_df)) selected_df = spline_df
      observations = size(series, 1)
      points = size(series, 2)
      if (order < 1 .or. selected_df < 1 .or. points < 3 .or. &
         observations <= order + 1 .or. .not. all(ieee_is_finite(series))) then
         model%info = 1
         return
      end if
      call cfar_design_matrix(series, order, selected_df, design, response)
      best_score = huge(1.0_dp)
      lower = log(0.02_dp)
      upper = log(20.0_dp)
      do iteration = 1, 4
         step = (upper - lower)/30.0_dp
         do candidate = 0, 30
            rho = exp(lower + real(candidate, dp)*step)
            call cfar_gls(design, response, observations - order, points, rho, &
               coefficient, fitted_vector, residual_vector, score, info)
            if (info == 0 .and. score < best_score) then
               best_score = score
               model%rho = rho
            end if
         end do
         lower = log(model%rho) - step
         upper = log(model%rho) + step
      end do
      if (.not. ieee_is_finite(best_score)) then
         model%info = 2
         return
      end if
      call cfar_gls(design, response, observations - order, points, model%rho, &
         coefficient, fitted_vector, residual_vector, score, info, model%sigma, &
         model%log_likelihood)
      if (info /= 0) then
         model%info = info
         return
      end if
      model%order = order
      model%spline_df = selected_df
      model%grid_size = points - 1
      model%kernel_coefficients = reshape(coefficient, [selected_df + 1, order])
      model%kernel_coefficients = transpose(model%kernel_coefficients)
      call cfar_sample_kernels(model%kernel_coefficients, points, model%kernel)
      model%fitted = reshape(fitted_vector, [points, observations - order])
      model%fitted = transpose(model%fitted)
      model%residuals = reshape(residual_vector, [points, observations - order])
      model%residuals = transpose(model%residuals)
   end function nts_cfar_estimate

   pure function nts_cfar_order_tests(series, maximum_order, spline_df) result(out)
      !! Perform sequential OU-weighted F tests for CFAR orders one through p.
      real(dp), intent(in) :: series(:, :) !! Time-by-grid functional observations.
      integer, intent(in) :: maximum_order !! Largest order to test.
      integer, intent(in), optional :: spline_df !! Spline degrees of freedom excluding intercept.
      type(nts_cfar_order_tests_t) :: out
      type(nts_cfar_model_t) :: full, restricted
      real(dp) :: full_sse, restricted_sse
      integer :: selected_df, tested, points, blocks, parameters, start

      selected_df = 4
      if (present(spline_df)) selected_df = spline_df
      points = size(series, 2)
      if (maximum_order < 1 .or. selected_df < 1 .or. &
         size(series, 1) <= maximum_order + 1 .or. points < 3 .or. &
         .not. all(ieee_is_finite(series))) then
         out%info = 1
         return
      end if
      allocate(out%statistic(maximum_order), out%p_value(maximum_order))
      allocate(out%numerator_df(maximum_order), out%denominator_df(maximum_order))
      do tested = 1, maximum_order
         full = nts_cfar_estimate(series, tested, selected_df)
         if (full%info /= 0) then
            out%info = full%info
            return
         end if
         blocks = size(series, 1) - tested
         full_sse = cfar_weighted_sse(full%residuals, full%rho)
         if (tested == 1) then
            restricted_sse = cfar_weighted_sse(series(2:, :), full%rho)
         else
            restricted = nts_cfar_estimate(series, tested - 1, selected_df)
            if (restricted%info /= 0) then
               out%info = restricted%info
               return
            end if
            start = size(restricted%residuals, 1) - blocks + 1
            restricted_sse = cfar_weighted_sse( &
               restricted%residuals(start:, :), full%rho)
         end if
         parameters = tested*(selected_df + 1)
         out%numerator_df(tested) = selected_df + 1
         out%denominator_df(tested) = blocks*points - parameters
         if (full_sse <= 0.0_dp .or. out%denominator_df(tested) <= 0) then
            out%info = 2
            return
         end if
         out%statistic(tested) = max(0.0_dp, restricted_sse - full_sse)* &
            real(out%denominator_df(tested), dp)/ &
            (real(out%numerator_df(tested), dp)*full_sse)
         out%p_value(tested) = f_upper_probability(out%statistic(tested), &
            out%numerator_df(tested), out%denominator_df(tested))
      end do
   end function nts_cfar_order_tests

   pure function nts_cfar_forecast(model, series, horizon) result(out)
      !! Recursively forecast complete functional observations from a CFAR fit.
      type(nts_cfar_model_t), intent(in) :: model !! Fitted CFAR model.
      real(dp), intent(in) :: series(:, :) !! Historical time-by-grid curves.
      integer, intent(in) :: horizon !! Positive forecast horizon.
      type(nts_cfar_forecast_t) :: out
      real(dp), allocatable :: history(:, :), aligned(:, :)
      integer :: points, time, lag, location, source, kernel_index

      points = model%grid_size + 1
      if (model%info /= 0 .or. model%order < 1 .or. horizon < 1 .or. &
         size(series, 1) < model%order .or. size(series, 2) < 2 .or. &
         .not. allocated(model%kernel) .or. .not. all(ieee_is_finite(series))) then
         out%info = 1
         return
      end if
      call cfar_resample_curves(series, points, aligned)
      allocate(history(size(series, 1) + horizon, points))
      history(:size(series, 1), :) = aligned
      do time = size(series, 1) + 1, size(series, 1) + horizon
         history(time, :) = 0.0_dp
         do lag = 1, model%order
            do location = 1, points
               do source = 1, points
                  kernel_index = points + location - source
                  history(time, location) = history(time, location) + &
                     model%kernel(kernel_index, lag)*history(time - lag, source)/ &
                     real(points, dp)
               end do
            end do
         end do
      end do
      out%mean = history(size(series, 1) + 1:, :)
   end function nts_cfar_forecast

   pure function nts_cfar_partial_forecast(model, series, partial_values) result(out)
      !! Forecast one curve conditionally on its observed leading values.
      type(nts_cfar_model_t), intent(in) :: model !! Fitted CFAR model.
      real(dp), intent(in) :: series(:, :) !! Historical complete curves.
      real(dp), intent(in) :: partial_values(:) !! Leading values of the next curve.
      type(nts_cfar_forecast_t) :: out
      type(nts_cfar_forecast_t) :: unconditional
      real(dp) :: distance, correction
      integer :: points, observed, location

      unconditional = nts_cfar_forecast(model, series, 1)
      points = model%grid_size + 1
      observed = size(partial_values)
      if (unconditional%info /= 0 .or. observed < 1 .or. observed >= points .or. &
         .not. all(ieee_is_finite(partial_values))) then
         out%info = 1
         return
      end if
      out%mean = unconditional%mean
      out%mean(1, :observed) = partial_values
      correction = partial_values(observed) - unconditional%mean(1, observed)
      do location = observed + 1, points
         distance = real(location - observed, dp)/real(points - 1, dp)
         out%mean(1, location) = unconditional%mean(1, location) + &
            correction*exp(-model%rho*distance)
      end do
   end function nts_cfar_partial_forecast

   pure function nts_cfar_irregular_estimate(values, positions, counts, order, &
      spline_df, grid_size) result(model)
      !! Estimate CFAR after linearly registering irregular curves to a common grid.
      real(dp), intent(in) :: values(:, :) !! Time-by-capacity irregular curve values.
      real(dp), intent(in) :: positions(:, :) !! Matching positions on the unit interval.
      integer, intent(in) :: counts(:) !! Number of valid entries in each row.
      integer, intent(in) :: order !! Positive functional autoregressive order.
      integer, intent(in), optional :: spline_df !! Spline degrees of freedom excluding intercept.
      integer, intent(in), optional :: grid_size !! Number of common-grid intervals.
      type(nts_cfar_model_t) :: model
      real(dp), allocatable :: registered(:, :)
      integer :: selected_df, selected_grid, info

      selected_df = 4
      selected_grid = 50
      if (present(spline_df)) selected_df = spline_df
      if (present(grid_size)) selected_grid = grid_size
      call cfar_register_irregular(values, positions, counts, selected_grid + 1, &
         registered, info)
      if (info /= 0) then
         model%info = info
         return
      end if
      model = nts_cfar_estimate(registered, order, selected_df)
   end function nts_cfar_irregular_estimate

   pure function nts_cfar_irregular_order_tests(values, positions, counts, &
      maximum_order, spline_df, grid_size) result(out)
      !! Test CFAR orders after registering irregularly observed curves.
      real(dp), intent(in) :: values(:, :) !! Time-by-capacity irregular curve values.
      real(dp), intent(in) :: positions(:, :) !! Matching positions on the unit interval.
      integer, intent(in) :: counts(:) !! Number of valid entries in each row.
      integer, intent(in) :: maximum_order !! Largest order to test.
      integer, intent(in), optional :: spline_df !! Spline degrees of freedom excluding intercept.
      integer, intent(in), optional :: grid_size !! Number of common-grid intervals.
      type(nts_cfar_order_tests_t) :: out
      real(dp), allocatable :: registered(:, :)
      integer :: selected_df, selected_grid, info

      selected_df = 4
      selected_grid = 50
      if (present(spline_df)) selected_df = spline_df
      if (present(grid_size)) selected_grid = grid_size
      call cfar_register_irregular(values, positions, counts, selected_grid + 1, &
         registered, info)
      if (info /= 0) then
         out%info = info
         return
      end if
      out = nts_cfar_order_tests(registered, maximum_order, selected_df)
   end function nts_cfar_irregular_order_tests

   pure function nts_msar_simulate_from_draws(intercept, ar, ar_order, &
      transition, innovation_sd, normal_draws, uniform_draws, burnin) result(out)
      !! Simulate a Markov-switching AR process from supplied random draws.
      real(dp), intent(in) :: intercept(:) !! Regime intercepts.
      real(dp), intent(in) :: ar(:, :) !! Regime-by-lag AR coefficients.
      integer, intent(in) :: ar_order(:) !! Regime-specific AR orders.
      real(dp), intent(in) :: transition(:, :) !! Row-stochastic transition matrix.
      real(dp), intent(in) :: innovation_sd(:) !! Regime innovation standard deviations.
      real(dp), intent(in) :: normal_draws(:) !! Supplied independent standard-normal draws.
      real(dp), intent(in) :: uniform_draws(:) !! Supplied independent uniform draws.
      integer, intent(in), optional :: burnin !! Number of leading observations to discard.
      type(nts_msar_simulation_t) :: out
      real(dp), allocatable :: complete_series(:), complete_innovations(:)
      integer, allocatable :: complete_state(:)
      real(dp) :: conditional_mean
      integer :: regimes, total, selected_burnin, memory, time, lag, regime

      regimes = size(intercept)
      total = size(normal_draws)
      selected_burnin = 0
      if (present(burnin)) selected_burnin = burnin
      if (regimes < 1 .or. total < 1 .or. selected_burnin < 0 .or. &
         selected_burnin >= total .or. size(ar, 1) /= regimes .or. &
         size(ar_order) /= regimes .or. size(innovation_sd) /= regimes .or. &
         any(shape(transition) /= [regimes, regimes]) .or. &
         size(uniform_draws) /= total .or. any(ar_order < 0) .or. &
         maxval(ar_order) > size(ar, 2) .or. any(innovation_sd <= 0.0_dp) .or. &
         any(transition < 0.0_dp) .or. &
         maxval(abs(sum(transition, dim=2) - 1.0_dp)) > 1.0e-10_dp .or. &
         any(uniform_draws < 0.0_dp .or. uniform_draws > 1.0_dp) .or. &
         .not. all(ieee_is_finite(intercept)) .or. &
         .not. all(ieee_is_finite(ar)) .or. &
         .not. all(ieee_is_finite(transition)) .or. &
         .not. all(ieee_is_finite(innovation_sd)) .or. &
         .not. all(ieee_is_finite(normal_draws))) then
         out%info = 1
         return
      end if
      memory = max(1, maxval(ar_order))
      if (total < memory) then
         out%info = 1
         return
      end if
      allocate(complete_series(total), complete_innovations(total))
      allocate(complete_state(total))
      do time = 1, memory
         complete_series(time) = normal_draws(time)
         complete_innovations(time) = normal_draws(time)
         complete_state(time) = min(regimes, 1 + &
            int(uniform_draws(time)*real(regimes, dp)))
      end do
      do time = memory + 1, total
         complete_state(time) = categorical_state( &
            transition(complete_state(time - 1), :), uniform_draws(time))
         regime = complete_state(time)
         conditional_mean = intercept(regime)
         do lag = 1, ar_order(regime)
            conditional_mean = conditional_mean + &
               ar(regime, lag)*complete_series(time - lag)
         end do
         complete_innovations(time) = innovation_sd(regime)*normal_draws(time)
         complete_series(time) = conditional_mean + complete_innovations(time)
      end do
      out%series = complete_series(selected_burnin + 1:total)
      out%innovations = complete_innovations(selected_burnin + 1:total)
      out%state = complete_state(selected_burnin + 1:total)
      out%intercept = intercept
      out%ar = ar
      out%ar_order = ar_order
      out%transition = transition
      out%innovation_sd = innovation_sd
      out%burnin = selected_burnin
   end function nts_msar_simulate_from_draws

   function nts_msar_simulate(intercept, ar, ar_order, transition, innovation_sd, &
      observations, burnin) result(out)
      !! Simulate a Markov-switching AR process using the shared RNG.
      real(dp), intent(in) :: intercept(:) !! Regime intercepts.
      real(dp), intent(in) :: ar(:, :) !! Regime-by-lag AR coefficients.
      integer, intent(in) :: ar_order(:) !! Regime-specific AR orders.
      real(dp), intent(in) :: transition(:, :) !! Row-stochastic transition matrix.
      real(dp), intent(in) :: innovation_sd(:) !! Regime innovation standard deviations.
      integer, intent(in) :: observations !! Number of retained observations.
      integer, intent(in), optional :: burnin !! Number of leading observations to discard.
      type(nts_msar_simulation_t) :: out
      real(dp), allocatable :: normal_draws(:), uniform_draws(:)
      integer :: selected_burnin, total, time

      selected_burnin = 500
      if (present(burnin)) selected_burnin = burnin
      if (observations < 1 .or. selected_burnin < 0) then
         out%info = 1
         return
      end if
      total = observations + selected_burnin
      allocate(normal_draws(total), uniform_draws(total))
      do time = 1, total
         normal_draws(time) = random_standard_normal()
         uniform_draws(time) = random_uniform()
      end do
      out = nts_msar_simulate_from_draws(intercept, ar, ar_order, transition, &
         innovation_sd, normal_draws, uniform_draws, selected_burnin)
   end function nts_msar_simulate

   pure function nts_utar_simulate_from_innovations(intercept, ar, ar_order, &
      thresholds, innovation_sd, delay, innovations, burnin) result(out)
      !! Simulate a univariate SETAR model from supplied standard innovations.
      real(dp), intent(in) :: intercept(:) !! Regime intercepts.
      real(dp), intent(in) :: ar(:, :) !! Regime-by-lag AR coefficients.
      integer, intent(in) :: ar_order(:) !! Regime-specific AR orders.
      real(dp), intent(in) :: thresholds(:) !! Ordered regime thresholds.
      real(dp), intent(in) :: innovation_sd(:) !! Regime innovation standard deviations.
      integer, intent(in) :: delay !! Delay of the threshold variable.
      real(dp), intent(in) :: innovations(:) !! Supplied standard-normal innovations.
      integer, intent(in) :: burnin !! Number of initial observations to discard.
      type(nts_tar_simulation_t) :: out
      real(dp), allocatable :: work(:), errors(:)
      integer :: regimes, maximum_order, start, total, time, lag, regime

      regimes = size(ar_order)
      total = size(innovations)
      if (.not. valid_tar_specification(intercept, ar, ar_order, thresholds, &
         innovation_sd, delay) .or. burnin < 0 .or. burnin >= total .or. &
         .not. all(ieee_is_finite(innovations))) then
         out%info = 1
         return
      end if
      maximum_order = max(1, maxval(ar_order))
      start = max(maximum_order, delay) + 1
      if (total < start) then
         out%info = 2
         return
      end if

      allocate(work(total), errors(total))
      work = 0.0_dp
      errors = 0.0_dp
      work(:start - 1) = innovation_sd(1)*innovations(:start - 1)
      errors(:start - 1) = work(:start - 1)
      do time = start, total
         regime = threshold_regime(work(time - delay), thresholds)
         errors(time) = innovation_sd(regime)*innovations(time)
         work(time) = intercept(regime) + errors(time)
         do lag = 1, ar_order(regime)
            work(time) = work(time) + ar(regime, lag)*work(time - lag)
         end do
      end do

      out%series = work(burnin + 1:)
      out%innovations = errors(burnin + 1:)
      out%intercept = intercept
      out%ar = ar
      out%ar_order = ar_order
      out%thresholds = thresholds
      out%innovation_sd = innovation_sd
      out%delay = delay
      out%burnin = burnin
   end function nts_utar_simulate_from_innovations

   function nts_utar_simulate(intercept, ar, ar_order, thresholds, innovation_sd, &
      delay, observations, burnin) result(out)
      !! Simulate a univariate SETAR model using the shared random stream.
      real(dp), intent(in) :: intercept(:) !! Regime intercepts.
      real(dp), intent(in) :: ar(:, :) !! Regime-by-lag AR coefficients.
      integer, intent(in) :: ar_order(:) !! Regime-specific AR orders.
      real(dp), intent(in) :: thresholds(:) !! Ordered regime thresholds.
      real(dp), intent(in) :: innovation_sd(:) !! Regime innovation standard deviations.
      integer, intent(in) :: delay !! Delay of the threshold variable.
      integer, intent(in) :: observations !! Number of retained observations.
      integer, intent(in), optional :: burnin !! Number of initial observations to discard.
      type(nts_tar_simulation_t) :: out
      real(dp), allocatable :: normal_matrix(:, :), innovations(:)
      integer :: discard

      discard = 500
      if (present(burnin)) discard = burnin
      if (observations < 1 .or. discard < 0) then
         out%info = 1
         return
      end if
      allocate(normal_matrix(1, observations + discard))
      call random_standard_normal_matrix(normal_matrix)
      innovations = normal_matrix(1, :)
      out = nts_utar_simulate_from_innovations(intercept, ar, ar_order, &
         thresholds, innovation_sd, delay, innovations, discard)
   end function nts_utar_simulate

   pure function nts_utar_estimate(series, ar_order, thresholds, delay, &
      threshold_variable, include_mean) result(out)
      !! Estimate a SETAR model with fixed thresholds by conditional least squares.
      real(dp), intent(in) :: series(:) !! Observed univariate time series.
      integer, intent(in) :: ar_order(:) !! Regime-specific AR orders.
      real(dp), intent(in) :: thresholds(:) !! Ordered fixed thresholds.
      integer, intent(in) :: delay !! Delay of the threshold variable.
      real(dp), intent(in), optional :: threshold_variable(:) !! External threshold series.
      logical, intent(in), optional :: include_mean(:) !! Regime intercept inclusion flags.
      type(nts_tar_model_t) :: out
      real(dp), allocatable :: lagged(:, :), response(:), threshold_values(:)
      real(dp), allocatable :: x_regime(:, :), y_regime(:)
      real(dp), allocatable :: beta(:), standard_error(:), residual(:)
      logical, allocatable :: means(:)
      integer, allocatable :: indices(:)
      real(dp) :: rss
      integer :: regimes, maximum_order, start, usable, regime, row, time
      integer :: regime_count, parameter_count, column, status

      regimes = size(ar_order)
      if (size(series) < 3 .or. regimes < 2 .or. &
         size(thresholds) /= regimes - 1 .or. delay < 1 .or. &
         any(ar_order < 0) .or. .not. strictly_increasing(thresholds) .or. &
         .not. all(ieee_is_finite(series)) .or. &
         .not. all(ieee_is_finite(thresholds))) then
         out%info = 1
         return
      end if
      if (present(threshold_variable)) then
         if (size(threshold_variable) /= size(series) .or. &
            .not. all(ieee_is_finite(threshold_variable))) then
            out%info = 2
            return
         end if
      end if
      allocate(means(regimes))
      means = .true.
      if (present(include_mean)) then
         if (size(include_mean) /= regimes) then
            out%info = 2
            return
         end if
         means = include_mean
      end if

      maximum_order = max(1, maxval(ar_order))
      start = max(maximum_order, delay) + 1
      usable = size(series) - start + 1
      if (usable < regimes + sum(ar_order) + count(means)) then
         out%info = 3
         return
      end if
      allocate(lagged(usable, maximum_order), response(usable))
      allocate(threshold_values(usable))
      do row = 1, usable
         time = start + row - 1
         response(row) = series(time)
         do column = 1, maximum_order
            lagged(row, column) = series(time - column)
         end do
         if (present(threshold_variable)) then
            threshold_values(row) = threshold_variable(time - delay)
         else
            threshold_values(row) = series(time - delay)
         end if
      end do

      allocate(out%coefficients(regimes, maximum_order + 1))
      allocate(out%innovation_sd(regimes), out%regime_observations(regimes))
      allocate(out%residuals(usable), out%standardized_residuals(usable))
      out%coefficients = 0.0_dp
      out%residuals = 0.0_dp
      out%standardized_residuals = 0.0_dp
      out%aic = 0.0_dp
      do regime = 1, regimes
         regime_count = count_regime(threshold_values, thresholds, regime)
         parameter_count = ar_order(regime) + merge(1, 0, means(regime))
         out%regime_observations(regime) = regime_count
         if (regime_count <= parameter_count .or. parameter_count < 1) then
            out%info = 4
            return
         end if
         allocate(indices(regime_count), x_regime(regime_count, parameter_count))
         allocate(y_regime(regime_count))
         call collect_regime_indices(threshold_values, thresholds, regime, indices)
         column = 0
         if (means(regime)) then
            column = 1
            x_regime(:, 1) = 1.0_dp
         end if
         if (ar_order(regime) > 0) then
            x_regime(:, column + 1:parameter_count) = &
               lagged(indices, :ar_order(regime))
         end if
         y_regime = response(indices)
         call ols_fit(x_regime, y_regime, beta, standard_error, residual, rss, status)
         if (status /= 0 .or. rss <= 0.0_dp) then
            out%info = 10 + status
            return
         end if
         column = 0
         if (means(regime)) then
            out%coefficients(regime, 1) = beta(1)
            column = 1
         end if
         if (ar_order(regime) > 0) then
            out%coefficients(regime, 2:ar_order(regime) + 1) = &
               beta(column + 1:parameter_count)
         end if
         out%innovation_sd(regime) = &
            sqrt(rss/real(regime_count - parameter_count, dp))
         out%residuals(indices) = residual
         out%standardized_residuals(indices) = residual/out%innovation_sd(regime)
         out%aic = out%aic + real(regime_count, dp)* &
            log(rss/real(regime_count, dp)) + &
            2.0_dp*real(parameter_count, dp)
         deallocate(indices, x_regime, y_regime, beta, standard_error, residual)
      end do
      out%data = series
      out%ar_order = ar_order
      out%thresholds = thresholds
      out%include_mean = means
      out%delay = delay
   end function nts_utar_estimate

   pure function nts_utar_threshold_search(series, first_order, second_order, &
      delay, trim, threshold_variable, include_mean) result(out)
      !! Search a two-regime SETAR threshold by exhaustive conditional least squares.
      real(dp), intent(in) :: series(:) !! Observed univariate time series.
      integer, intent(in) :: first_order !! AR order below the threshold.
      integer, intent(in) :: second_order !! AR order above the threshold.
      integer, intent(in) :: delay !! Delay of the threshold variable.
      real(dp), intent(in), optional :: trim(2) !! Lower and upper search quantiles.
      real(dp), intent(in), optional :: threshold_variable(:) !! External threshold series.
      logical, intent(in), optional :: include_mean !! Whether both regimes include intercepts.
      type(nts_tar_search_t) :: out
      type(nts_tar_model_t) :: candidate_model
      real(dp), allocatable :: threshold_values(:), ordered(:), unique(:)
      real(dp) :: selected_trim(2), lower, upper, candidate_rss
      integer :: orders(2), maximum_order, start, usable, row, time, candidates
      integer :: candidate, best
      logical :: means(2), selected_mean

      selected_trim = [0.1_dp, 0.9_dp]
      if (present(trim)) selected_trim = trim
      selected_mean = .true.
      if (present(include_mean)) selected_mean = include_mean
      means = selected_mean
      orders = [first_order, second_order]
      if (size(series) < 5 .or. any(orders < 0) .or. delay < 1 .or. &
         selected_trim(1) < 0.0_dp .or. selected_trim(2) > 1.0_dp .or. &
         selected_trim(1) >= selected_trim(2)) then
         out%info = 1
         return
      end if
      if (present(threshold_variable)) then
         if (size(threshold_variable) /= size(series)) then
            out%info = 2
            return
         end if
      end if

      maximum_order = max(1, maxval(orders))
      start = max(maximum_order, delay) + 1
      usable = size(series) - start + 1
      if (usable < 4) then
         out%info = 3
         return
      end if
      allocate(threshold_values(usable))
      do row = 1, usable
         time = start + row - 1
         if (present(threshold_variable)) then
            threshold_values(row) = threshold_variable(time - delay)
         else
            threshold_values(row) = series(time - delay)
         end if
      end do
      ordered = sorted(threshold_values)
      lower = quantile(ordered, selected_trim(1))
      upper = quantile(ordered, selected_trim(2))
      allocate(unique(usable))
      candidates = 0
      do row = 1, usable
         if (ordered(row) < lower .or. ordered(row) > upper) cycle
         if (candidates == 0) then
            candidates = candidates + 1
            unique(candidates) = ordered(row)
         else if (ordered(row) > unique(candidates)) then
            candidates = candidates + 1
            unique(candidates) = ordered(row)
         end if
      end do
      if (candidates < 1) then
         out%info = 4
         return
      end if
      allocate(out%candidates(candidates), out%rss(candidates))
      out%candidates = unique(:candidates)
      out%rss = huge(1.0_dp)
      best = 0
      do candidate = 1, candidates
         candidate_model = nts_utar_estimate(series, orders, &
            out%candidates(candidate:candidate), delay, threshold_variable, means)
         if (candidate_model%info /= 0) cycle
         candidate_rss = sum(candidate_model%residuals**2)
         out%rss(candidate) = candidate_rss
         if (best == 0) then
            best = candidate
            out%model = candidate_model
         else if (candidate_rss < out%rss(best)) then
            best = candidate
            out%model = candidate_model
         end if
      end do
      if (best == 0) then
         out%info = 5
         return
      end if
      out%selected = best
   end function nts_utar_threshold_search

   pure function nts_utar_forecast_draws(model, origin, normal_draws, level) &
      result(out)
      !! Forecast a fitted self-exciting SETAR model from supplied normal draws.
      type(nts_tar_model_t), intent(in) :: model !! Fitted SETAR model.
      integer, intent(in) :: origin !! Forecast origin in the fitted data.
      real(dp), intent(in) :: normal_draws(:, :) !! Horizon-by-simulation normal draws.
      real(dp), intent(in), optional :: level !! Pointwise interval coverage.
      type(nts_tar_forecast_t) :: out
      real(dp), allocatable :: history(:), ordered(:)
      real(dp) :: selected_level, tail, value
      integer :: horizon, iterations, simulation, step, time, lag, regime

      selected_level = 0.95_dp
      if (present(level)) selected_level = level
      horizon = size(normal_draws, 1)
      iterations = size(normal_draws, 2)
      if (model%info /= 0 .or. .not. allocated(model%data) .or. &
         origin < 1 .or. origin > size(model%data) .or. horizon < 1 .or. &
         iterations < 1 .or. selected_level <= 0.0_dp .or. &
         selected_level >= 1.0_dp .or. .not. all(ieee_is_finite(normal_draws))) then
         out%info = 1
         return
      end if
      if (origin < max(maxval(model%ar_order), model%delay)) then
         out%info = 2
         return
      end if
      allocate(out%simulations(horizon, iterations))
      do simulation = 1, iterations
         allocate(history(origin + horizon))
         history(:origin) = model%data(:origin)
         do step = 1, horizon
            time = origin + step
            regime = threshold_regime(history(time - model%delay), &
               model%thresholds)
            value = model%coefficients(regime, 1)
            do lag = 1, model%ar_order(regime)
               value = value + model%coefficients(regime, lag + 1)* &
                  history(time - lag)
            end do
            history(time) = value + model%innovation_sd(regime)* &
               normal_draws(step, simulation)
            out%simulations(step, simulation) = history(time)
         end do
         deallocate(history)
      end do
      allocate(out%mean(horizon), out%lower(horizon), out%upper(horizon))
      tail = 0.5_dp*(1.0_dp - selected_level)
      do step = 1, horizon
         ordered = sorted(out%simulations(step, :))
         out%mean(step) = sum(ordered)/real(iterations, dp)
         out%lower(step) = quantile(ordered, tail)
         out%upper(step) = quantile(ordered, 1.0_dp - tail)
      end do
      out%level = selected_level
      out%origin = origin
   end function nts_utar_forecast_draws

   function nts_utar_forecast(model, origin, horizon, iterations, level) result(out)
      !! Forecast a fitted SETAR model using the shared random stream.
      type(nts_tar_model_t), intent(in) :: model !! Fitted SETAR model.
      integer, intent(in) :: origin !! Forecast origin in the fitted data.
      integer, intent(in) :: horizon !! Number of forecast steps.
      integer, intent(in) :: iterations !! Number of simulated forecast paths.
      real(dp), intent(in), optional :: level !! Pointwise interval coverage.
      type(nts_tar_forecast_t) :: out
      real(dp), allocatable :: normal_draws(:, :)

      if (horizon < 1 .or. iterations < 1) then
         out%info = 1
         return
      end if
      allocate(normal_draws(horizon, iterations))
      call random_standard_normal_matrix(normal_draws)
      out = nts_utar_forecast_draws(model, origin, normal_draws, level)
   end function nts_utar_forecast

   pure function nts_mtar_simulate_from_standard(intercept, ar, ar_order, &
      covariance, thresholds, threshold_component, delay, standard_draws, &
      burnin) result(out)
      !! Simulate a multivariate SETAR model from supplied standard-normal draws.
      real(dp), intent(in) :: intercept(:, :) !! Variable-by-regime intercepts.
      real(dp), intent(in) :: ar(:, :, :, :) !! VAR coefficients by lag and regime.
      integer, intent(in) :: ar_order(:) !! Regime-specific VAR orders.
      real(dp), intent(in) :: covariance(:, :, :) !! Regime innovation covariances.
      real(dp), intent(in) :: thresholds(:) !! Ordered regime thresholds.
      integer, intent(in) :: threshold_component !! Component driving regime selection.
      integer, intent(in) :: delay !! Delay of the threshold component.
      real(dp), intent(in) :: standard_draws(:, :) !! Time-by-variable standard normals.
      integer, intent(in) :: burnin !! Number of initial observations to discard.
      type(nts_mtar_simulation_t) :: out
      real(dp), allocatable :: work(:, :), errors(:, :), lower(:, :, :)
      real(dp), allocatable :: lower_regime(:, :)
      integer :: variables, regimes, total, maximum_order, start
      integer :: regime, time, lag, status

      variables = size(intercept, 1)
      regimes = size(intercept, 2)
      total = size(standard_draws, 1)
      if (variables < 1 .or. regimes < 2 .or. size(ar_order) /= regimes .or. &
         any(ar_order < 0) .or. size(ar, 1) /= variables .or. &
         size(ar, 2) /= variables .or. size(ar, 4) /= regimes .or. &
         any(ar_order > size(ar, 3)) .or. &
         any(shape(covariance) /= [variables, variables, regimes]) .or. &
         size(thresholds) /= regimes - 1 .or. &
         .not. strictly_increasing(thresholds) .or. threshold_component < 1 .or. &
         threshold_component > variables .or. delay < 1 .or. burnin < 0 .or. &
         burnin >= total .or. size(standard_draws, 2) /= variables .or. &
         .not. all(ieee_is_finite(intercept)) .or. &
         .not. all(ieee_is_finite(ar)) .or. &
         .not. all(ieee_is_finite(covariance)) .or. &
         .not. all(ieee_is_finite(thresholds)) .or. &
         .not. all(ieee_is_finite(standard_draws))) then
         out%info = 1
         return
      end if
      maximum_order = max(1, maxval(ar_order))
      start = max(maximum_order, delay) + 1
      if (total < start) then
         out%info = 2
         return
      end if
      allocate(lower(variables, variables, regimes))
      do regime = 1, regimes
         call cholesky_lower(covariance(:, :, regime), &
            lower_regime, status)
         if (status /= 0) then
            out%info = 10 + regime
            return
         end if
         lower(:, :, regime) = lower_regime
         deallocate(lower_regime)
      end do

      allocate(work(total, variables), errors(total, variables))
      work = 0.0_dp
      errors = 0.0_dp
      do time = 1, start - 1
         errors(time, :) = matmul(lower(:, :, 1), standard_draws(time, :))
         work(time, :) = intercept(:, 1) + errors(time, :)
      end do
      do time = start, total
         regime = threshold_regime(work(time - delay, threshold_component), &
            thresholds)
         errors(time, :) = matmul(lower(:, :, regime), standard_draws(time, :))
         work(time, :) = intercept(:, regime) + errors(time, :)
         do lag = 1, ar_order(regime)
            work(time, :) = work(time, :) + &
               matmul(ar(:, :, lag, regime), work(time - lag, :))
         end do
      end do

      out%series = work(burnin + 1:, :)
      out%innovations = errors(burnin + 1:, :)
      out%intercept = intercept
      out%ar = ar
      out%ar_order = ar_order
      out%covariance = covariance
      out%thresholds = thresholds
      allocate(out%regime_observations(regimes))
      out%regime_observations = 0
      do time = burnin + 1, total
         regime = 1
         if (time > delay) regime = threshold_regime( &
            work(time - delay, threshold_component), thresholds)
         out%regime_observations(regime) = out%regime_observations(regime) + 1
      end do
      out%threshold_component = threshold_component
      out%delay = delay
      out%burnin = burnin
   end function nts_mtar_simulate_from_standard

   function nts_mtar_simulate(intercept, ar, ar_order, covariance, thresholds, &
      threshold_component, delay, observations, burnin) result(out)
      !! Simulate a multivariate SETAR model using the shared random stream.
      real(dp), intent(in) :: intercept(:, :) !! Variable-by-regime intercepts.
      real(dp), intent(in) :: ar(:, :, :, :) !! VAR coefficients by lag and regime.
      integer, intent(in) :: ar_order(:) !! Regime-specific VAR orders.
      real(dp), intent(in) :: covariance(:, :, :) !! Regime innovation covariances.
      real(dp), intent(in) :: thresholds(:) !! Ordered regime thresholds.
      integer, intent(in) :: threshold_component !! Component driving regime selection.
      integer, intent(in) :: delay !! Delay of the threshold component.
      integer, intent(in) :: observations !! Number of retained observations.
      integer, intent(in), optional :: burnin !! Number of initial observations to discard.
      type(nts_mtar_simulation_t) :: out
      real(dp), allocatable :: draws_transposed(:, :), standard_draws(:, :)
      integer :: discard, variables

      discard = 500
      if (present(burnin)) discard = burnin
      variables = size(intercept, 1)
      if (observations < 1 .or. discard < 0 .or. variables < 1) then
         out%info = 1
         return
      end if
      allocate(draws_transposed(variables, observations + discard))
      call random_standard_normal_matrix(draws_transposed)
      standard_draws = transpose(draws_transposed)
      out = nts_mtar_simulate_from_standard(intercept, ar, ar_order, covariance, &
         thresholds, threshold_component, delay, standard_draws, discard)
   end function nts_mtar_simulate

   pure function nts_mtar_estimate(series, ar_order, thresholds, &
      threshold_component, delay, threshold_variable, include_mean) result(out)
      !! Estimate a multivariate SETAR model with fixed thresholds.
      real(dp), intent(in) :: series(:, :) !! Time-by-variable observations.
      integer, intent(in) :: ar_order(:) !! Regime-specific VAR orders.
      real(dp), intent(in) :: thresholds(:) !! Ordered fixed thresholds.
      integer, intent(in) :: threshold_component !! Component driving regime selection.
      integer, intent(in) :: delay !! Delay of the threshold variable.
      real(dp), intent(in), optional :: threshold_variable(:) !! External threshold series.
      logical, intent(in), optional :: include_mean(:) !! Regime intercept inclusion flags.
      type(nts_mtar_model_t) :: out
      real(dp), allocatable :: lagged(:, :), response(:, :), threshold_values(:)
      real(dp), allocatable :: x_regime(:, :), y_regime(:), residual(:, :)
      real(dp), allocatable :: beta(:), standard_error(:), equation_residual(:)
      real(dp), allocatable :: inverse_root(:, :)
      logical, allocatable :: means(:)
      integer, allocatable :: indices(:)
      real(dp) :: rss, log_determinant
      integer :: variables, regimes, maximum_order, start, usable
      integer :: row, time, lag, regime, regime_count, parameter_count
      integer :: equation, offset, status

      variables = size(series, 2)
      regimes = size(ar_order)
      if (size(series, 1) < 3 .or. variables < 1 .or. regimes < 2 .or. &
         size(thresholds) /= regimes - 1 .or. any(ar_order < 0) .or. &
         .not. strictly_increasing(thresholds) .or. threshold_component < 1 .or. &
         threshold_component > variables .or. delay < 1 .or. &
         .not. all(ieee_is_finite(series)) .or. &
         .not. all(ieee_is_finite(thresholds))) then
         out%info = 1
         return
      end if
      if (present(threshold_variable)) then
         if (size(threshold_variable) /= size(series, 1) .or. &
            .not. all(ieee_is_finite(threshold_variable))) then
            out%info = 2
            return
         end if
      end if
      allocate(means(regimes))
      means = .true.
      if (present(include_mean)) then
         if (size(include_mean) /= regimes) then
            out%info = 2
            return
         end if
         means = include_mean
      end if

      maximum_order = max(1, maxval(ar_order))
      start = max(maximum_order, delay) + 1
      usable = size(series, 1) - start + 1
      if (usable < regimes + variables*sum(ar_order) + count(means)) then
         out%info = 3
         return
      end if
      allocate(lagged(usable, variables*maximum_order))
      allocate(response(usable, variables), threshold_values(usable))
      do row = 1, usable
         time = start + row - 1
         response(row, :) = series(time, :)
         do lag = 1, maximum_order
            offset = (lag - 1)*variables
            lagged(row, offset + 1:offset + variables) = series(time - lag, :)
         end do
         if (present(threshold_variable)) then
            threshold_values(row) = threshold_variable(time - delay)
         else
            threshold_values(row) = series(time - delay, threshold_component)
         end if
      end do

      allocate(out%intercept(variables, regimes))
      allocate(out%ar(variables, variables, maximum_order, regimes))
      allocate(out%covariance(variables, variables, regimes))
      allocate(out%residuals(usable, variables))
      allocate(out%standardized_residuals(usable, variables))
      allocate(out%regime_observations(regimes))
      out%intercept = 0.0_dp
      out%ar = 0.0_dp
      out%covariance = 0.0_dp
      out%residuals = 0.0_dp
      out%standardized_residuals = 0.0_dp
      out%aic = 0.0_dp
      do regime = 1, regimes
         regime_count = count_regime(threshold_values, thresholds, regime)
         parameter_count = variables*ar_order(regime) + &
            merge(1, 0, means(regime))
         out%regime_observations(regime) = regime_count
         if (regime_count <= parameter_count .or. parameter_count < 1) then
            out%info = 4
            return
         end if
         allocate(indices(regime_count), x_regime(regime_count, parameter_count))
         allocate(residual(regime_count, variables), y_regime(regime_count))
         call collect_regime_indices(threshold_values, thresholds, regime, indices)
         offset = 0
         if (means(regime)) then
            x_regime(:, 1) = 1.0_dp
            offset = 1
         end if
         if (ar_order(regime) > 0) then
            x_regime(:, offset + 1:parameter_count) = &
               lagged(indices, :variables*ar_order(regime))
         end if
         do equation = 1, variables
            y_regime = response(indices, equation)
            call ols_fit(x_regime, y_regime, beta, standard_error, &
               equation_residual, rss, status)
            if (status /= 0) then
               out%info = 10 + status
               return
            end if
            residual(:, equation) = equation_residual
            offset = 0
            if (means(regime)) then
               out%intercept(equation, regime) = beta(1)
               offset = 1
            end if
            do lag = 1, ar_order(regime)
               out%ar(equation, :, lag, regime) = beta( &
                  offset + (lag - 1)*variables + 1:offset + lag*variables)
            end do
         end do
         out%covariance(:, :, regime) = &
            matmul(transpose(residual), residual)/real(regime_count, dp)
         call symmetric_inverse_root(out%covariance(:, :, regime), inverse_root, &
            log_determinant, status)
         if (status /= 0) then
            out%info = 20 + regime
            return
         end if
         out%residuals(indices, :) = residual
         out%standardized_residuals(indices, :) = matmul(residual, inverse_root)
         out%aic = out%aic + real(regime_count, dp)*log_determinant + &
            2.0_dp*real(variables*parameter_count, dp)
         deallocate(indices, x_regime, residual, y_regime)
         deallocate(beta, standard_error, equation_residual, inverse_root)
      end do
      out%data = series
      out%ar_order = ar_order
      out%thresholds = thresholds
      out%include_mean = means
      out%threshold_component = threshold_component
      out%delay = delay
   end function nts_mtar_estimate

   pure function nts_mtar_threshold_search(series, first_order, second_order, &
      threshold_component, delay, trim, threshold_variable, include_mean, &
      criterion) result(out)
      !! Search a two-regime multivariate SETAR threshold.
      real(dp), intent(in) :: series(:, :) !! Time-by-variable observations.
      integer, intent(in) :: first_order !! VAR order below the threshold.
      integer, intent(in) :: second_order !! VAR order above the threshold.
      integer, intent(in) :: threshold_component !! Component driving regime selection.
      integer, intent(in) :: delay !! Delay of the threshold variable.
      real(dp), intent(in), optional :: trim(2) !! Lower and upper search quantiles.
      real(dp), intent(in), optional :: threshold_variable(:) !! External threshold series.
      logical, intent(in), optional :: include_mean !! Whether regimes include intercepts.
      character(len=*), intent(in), optional :: criterion !! `aic` or `determinant` score.
      type(nts_mtar_search_t) :: out
      type(nts_mtar_model_t) :: candidate_model
      real(dp), allocatable :: threshold_values(:), ordered(:), unique(:)
      real(dp), allocatable :: inverse_root(:, :)
      real(dp) :: selected_trim(2), lower, upper, candidate_score, log_determinant
      integer :: orders(2), maximum_order, start, usable, row, time, candidates
      integer :: candidate, best, regime, status
      logical :: means(2), selected_mean, determinant_score
      character(len=:), allocatable :: selected_criterion

      selected_trim = [0.1_dp, 0.9_dp]
      if (present(trim)) selected_trim = trim
      selected_mean = .true.
      if (present(include_mean)) selected_mean = include_mean
      means = selected_mean
      selected_criterion = 'aic'
      if (present(criterion)) selected_criterion = adjustl(lowercase(criterion))
      determinant_score = selected_criterion == 'determinant'
      orders = [first_order, second_order]
      if (size(series, 1) < 5 .or. size(series, 2) < 1 .or. &
         any(orders < 0) .or. threshold_component < 1 .or. &
         threshold_component > size(series, 2) .or. delay < 1 .or. &
         selected_trim(1) < 0.0_dp .or. selected_trim(2) > 1.0_dp .or. &
         selected_trim(1) >= selected_trim(2) .or. &
         (selected_criterion /= 'aic' .and. .not. determinant_score)) then
         out%info = 1
         return
      end if
      if (present(threshold_variable)) then
         if (size(threshold_variable) /= size(series, 1)) then
            out%info = 2
            return
         end if
      end if

      maximum_order = max(1, maxval(orders))
      start = max(maximum_order, delay) + 1
      usable = size(series, 1) - start + 1
      allocate(threshold_values(usable))
      do row = 1, usable
         time = start + row - 1
         if (present(threshold_variable)) then
            threshold_values(row) = threshold_variable(time - delay)
         else
            threshold_values(row) = series(time - delay, threshold_component)
         end if
      end do
      ordered = sorted(threshold_values)
      lower = quantile(ordered, selected_trim(1))
      upper = quantile(ordered, selected_trim(2))
      allocate(unique(usable))
      candidates = 0
      do row = 1, usable
         if (ordered(row) < lower .or. ordered(row) > upper) cycle
         if (candidates == 0) then
            candidates = candidates + 1
            unique(candidates) = ordered(row)
         else if (ordered(row) > unique(candidates)) then
            candidates = candidates + 1
            unique(candidates) = ordered(row)
         end if
      end do
      if (candidates < 1) then
         out%info = 3
         return
      end if
      allocate(out%candidates(candidates), out%score(candidates))
      out%candidates = unique(:candidates)
      out%score = huge(1.0_dp)
      best = 0
      do candidate = 1, candidates
         candidate_model = nts_mtar_estimate(series, orders, &
            out%candidates(candidate:candidate), threshold_component, delay, &
            threshold_variable, means)
         if (candidate_model%info /= 0) cycle
         if (determinant_score) then
            candidate_score = 0.0_dp
            do regime = 1, 2
               call symmetric_inverse_root(candidate_model%covariance(:, :, regime), &
                  inverse_root, log_determinant, status)
               if (status /= 0) then
                  candidate_score = huge(1.0_dp)
                  exit
               end if
               candidate_score = candidate_score + exp(log_determinant + &
                  real(size(series, 2), dp)* &
                  log(real(candidate_model%regime_observations(regime), dp)))
               deallocate(inverse_root)
            end do
         else
            candidate_score = candidate_model%aic
         end if
         out%score(candidate) = candidate_score
         if (best == 0) then
            best = candidate
            out%model = candidate_model
         else if (candidate_score < out%score(best)) then
            best = candidate
            out%model = candidate_model
         end if
      end do
      if (best == 0) then
         out%info = 4
         return
      end if
      out%selected = best
   end function nts_mtar_threshold_search

   pure function nts_mtar_forecast_draws(model, origin, standard_draws, level) &
      result(out)
      !! Forecast a fitted multivariate SETAR model from supplied normal draws.
      type(nts_mtar_model_t), intent(in) :: model !! Fitted multivariate SETAR model.
      integer, intent(in) :: origin !! Forecast origin in the fitted data.
      real(dp), intent(in) :: standard_draws(:, :, :) !! Variable-by-horizon-by-path draws.
      real(dp), intent(in), optional :: level !! Pointwise interval coverage.
      type(nts_mtar_forecast_t) :: out
      real(dp), allocatable :: history(:, :), lower_factor(:, :, :), ordered(:)
      real(dp), allocatable :: lower_regime(:, :)
      real(dp) :: selected_level, tail
      integer :: variables, horizon, iterations, regimes
      integer :: regime, simulation, step, time, lag, variable, status

      selected_level = 0.95_dp
      if (present(level)) selected_level = level
      variables = size(standard_draws, 1)
      horizon = size(standard_draws, 2)
      iterations = size(standard_draws, 3)
      if (model%info /= 0 .or. .not. allocated(model%data) .or. &
         variables /= size(model%data, 2) .or. origin < 1 .or. &
         origin > size(model%data, 1) .or. horizon < 1 .or. iterations < 1 .or. &
         selected_level <= 0.0_dp .or. selected_level >= 1.0_dp .or. &
         .not. all(ieee_is_finite(standard_draws))) then
         out%info = 1
         return
      end if
      if (origin < max(maxval(model%ar_order), model%delay)) then
         out%info = 2
         return
      end if
      regimes = size(model%ar_order)
      allocate(lower_factor(variables, variables, regimes))
      do regime = 1, regimes
         call cholesky_lower(model%covariance(:, :, regime), &
            lower_regime, status)
         if (status /= 0) then
            out%info = 10 + regime
            return
         end if
         lower_factor(:, :, regime) = lower_regime
         deallocate(lower_regime)
      end do
      allocate(out%simulations(horizon, variables, iterations))
      do simulation = 1, iterations
         allocate(history(origin + horizon, variables))
         history(:origin, :) = model%data(:origin, :)
         do step = 1, horizon
            time = origin + step
            regime = threshold_regime(history(time - model%delay, &
               model%threshold_component), model%thresholds)
            history(time, :) = model%intercept(:, regime) + &
               matmul(lower_factor(:, :, regime), standard_draws(:, step, simulation))
            do lag = 1, model%ar_order(regime)
               history(time, :) = history(time, :) + &
                  matmul(model%ar(:, :, lag, regime), history(time - lag, :))
            end do
            out%simulations(step, :, simulation) = history(time, :)
         end do
         deallocate(history)
      end do
      allocate(out%mean(horizon, variables), out%lower(horizon, variables))
      allocate(out%upper(horizon, variables))
      tail = 0.5_dp*(1.0_dp - selected_level)
      do step = 1, horizon
         do variable = 1, variables
            ordered = sorted(out%simulations(step, variable, :))
            out%mean(step, variable) = sum(ordered)/real(iterations, dp)
            out%lower(step, variable) = quantile(ordered, tail)
            out%upper(step, variable) = quantile(ordered, 1.0_dp - tail)
         end do
      end do
      out%level = selected_level
      out%origin = origin
   end function nts_mtar_forecast_draws

   function nts_mtar_forecast(model, origin, horizon, iterations, level) result(out)
      !! Forecast a fitted multivariate SETAR model using the shared random stream.
      type(nts_mtar_model_t), intent(in) :: model !! Fitted multivariate SETAR model.
      integer, intent(in) :: origin !! Forecast origin in the fitted data.
      integer, intent(in) :: horizon !! Number of forecast steps.
      integer, intent(in) :: iterations !! Number of simulated forecast paths.
      real(dp), intent(in), optional :: level !! Pointwise interval coverage.
      type(nts_mtar_forecast_t) :: out
      real(dp), allocatable :: draw_matrix(:, :), standard_draws(:, :, :)
      integer :: variables

      if (.not. allocated(model%data) .or. horizon < 1 .or. iterations < 1) then
         out%info = 1
         return
      end if
      variables = size(model%data, 2)
      allocate(draw_matrix(variables, horizon*iterations))
      call random_standard_normal_matrix(draw_matrix)
      standard_draws = reshape(draw_matrix, [variables, horizon, iterations])
      out = nts_mtar_forecast_draws(model, origin, standard_draws, level)
   end function nts_mtar_forecast

   pure function nts_threshold_test(series, ar_order, delay, threshold_variable, &
      initial_count, include_mean) result(out)
      !! Apply NTS's arranged recursive least-squares threshold test.
      real(dp), intent(in) :: series(:) !! Univariate time series.
      integer, intent(in), optional :: ar_order !! Null-model AR order.
      integer, intent(in), optional :: delay !! Threshold-variable delay.
      real(dp), intent(in), optional :: threshold_variable(:) !! Aligned threshold values.
      integer, intent(in), optional :: initial_count !! Initial arranged-regression size.
      logical, intent(in), optional :: include_mean !! Whether the null AR includes an intercept.
      type(nts_nonlinearity_test_t) :: out
      real(dp), allocatable :: target(:), design(:, :), threshold(:)
      real(dp), allocatable :: ordered_target(:), ordered_design(:, :)
      real(dp), allocatable :: cross_product(:, :), inverse(:, :), beta(:)
      real(dp), allocatable :: work(:), standardized(:), auxiliary(:, :)
      real(dp), allocatable :: standard_error(:), residual(:)
      integer, allocatable :: index(:)
      real(dp) :: denominator, error, rss, numerator
      integer :: order, selected_delay, selected_initial, start, effective
      integer :: predictors, row, lag, status, h
      logical :: selected_mean

      order = 1
      if (present(ar_order)) order = ar_order
      selected_delay = 1
      if (present(delay)) selected_delay = delay
      selected_initial = 40
      if (present(initial_count)) selected_initial = initial_count
      selected_mean = .true.
      if (present(include_mean)) selected_mean = include_mean
      if (order < 1 .or. selected_delay < 1 .or. selected_initial < 1 .or. &
         .not. all(ieee_is_finite(series))) then
         out%info = 1
         return
      end if
      start = max(order, selected_delay) + 1
      effective = size(series) - start + 1
      predictors = order + merge(1, 0, selected_mean)
      do while (predictors > selected_initial)
         selected_initial = selected_initial + 10
      end do
      if (effective <= selected_initial + predictors) then
         out%info = 1
         return
      end if
      allocate(target(effective), design(effective, predictors), threshold(effective))
      target = series(start:)
      row = 0
      if (selected_mean) then
         row = 1
         design(:, row) = 1.0_dp
      end if
      do lag = 1, order
         row = row + 1
         design(:, row) = series(start - lag:size(series) - lag)
      end do
      if (present(threshold_variable)) then
         if (size(threshold_variable) >= effective) then
            threshold = threshold_variable(:effective)
         else
            threshold = series(start - selected_delay:size(series) - selected_delay)
         end if
      else
         threshold = series(start - selected_delay:size(series) - selected_delay)
      end if
      if (.not. all(ieee_is_finite(threshold))) then
         out%info = 1
         return
      end if
      call ascending_indices(threshold, index)
      allocate(ordered_target(effective), ordered_design(effective, predictors))
      do row = 1, effective
         ordered_target(row) = target(index(row))
         ordered_design(row, :) = design(index(row), :)
      end do
      cross_product = matmul(transpose(ordered_design(:selected_initial, :)), &
         ordered_design(:selected_initial, :))
      call invert_matrix(cross_product, inverse, status)
      if (status /= 0) then
         out%info = 2
         return
      end if
      beta = matmul(inverse, matmul(transpose( &
         ordered_design(:selected_initial, :)), &
         ordered_target(:selected_initial)))
      allocate(standardized(effective - selected_initial), work(predictors))
      do row = selected_initial + 1, effective
         work = matmul(inverse, ordered_design(row, :))
         error = dot_product(ordered_design(row, :), beta) - ordered_target(row)
         denominator = 1.0_dp + dot_product(ordered_design(row, :), work)
         if (denominator <= tiny(1.0_dp)) then
            out%info = 3
            return
         end if
         beta = beta - work*error/denominator
         inverse = inverse - matmul(reshape(work, [predictors, 1]), &
            reshape(work, [1, predictors]))/denominator
         standardized(row - selected_initial) = error/sqrt(denominator)
      end do
      if (selected_mean) then
         auxiliary = ordered_design(selected_initial + 1:, :)
      else
         allocate(auxiliary(size(standardized), predictors + 1))
         auxiliary(:, 1) = 1.0_dp
         auxiliary(:, 2:) = ordered_design(selected_initial + 1:, :)
      end if
      call ols_fit(auxiliary, standardized, beta, standard_error, residual, &
         rss, status)
      h = max(1, order + 1 - selected_delay)
      out%numerator_df = order + 1
      out%denominator_df = size(series) - selected_delay - selected_initial - &
         order - h
      out%initial_count = selected_initial
      if (status /= 0 .or. out%denominator_df <= 0 .or. rss <= 0.0_dp) then
         out%info = 4
         return
      end if
      numerator = max(0.0_dp, sum(standardized**2) - rss)
      out%statistic = (numerator/real(out%numerator_df, dp))/ &
         (rss/real(out%denominator_df, dp))
      out%p_value = f_upper_probability(out%statistic, out%numerator_df, &
         out%denominator_df)
   end function nts_threshold_test

   pure function nts_mtar_refine(model, t_threshold, threshold_variable) result(out)
      !! Refit a multivariate TAR after selecting coefficients by absolute t ratio.
      type(nts_mtar_model_t), intent(in) :: model !! Fitted multivariate TAR model.
      real(dp), intent(in), optional :: t_threshold !! Minimum retained absolute t ratio.
      real(dp), intent(in), optional :: threshold_variable(:) !! Aligned external threshold values.
      type(nts_mtar_refinement_t) :: out
      real(dp), allocatable :: target(:, :), design(:, :), threshold(:)
      real(dp), allocatable :: selected_design(:, :), beta(:), standard_error(:)
      real(dp), allocatable :: residual(:), covariance(:, :), inverse_root(:, :)
      real(dp), allocatable :: pooled(:, :), inverse(:, :)
      integer, allocatable :: indices(:), predictors(:)
      logical, allocatable :: candidate(:), keep(:)
      real(dp) :: rss, log_determinant
      integer :: variables, regimes, maximum_order, start, effective, maximum_predictors
      integer :: regime, equation, row, lag, source, predictor, count_value
      integer :: retained_count, status, total_parameters, total_observations

      out%threshold = 1.0_dp
      if (present(t_threshold)) out%threshold = t_threshold
      if (model%info /= 0 .or. .not. allocated(model%data) .or. &
         out%threshold < 0.0_dp) then
         out%info = 1
         return
      end if
      variables = size(model%data, 2)
      regimes = size(model%ar_order)
      maximum_order = maxval(model%ar_order)
      start = max(maximum_order, model%delay) + 1
      effective = size(model%data, 1) - start + 1
      maximum_predictors = 1 + variables*maximum_order
      if (variables < 1 .or. regimes < 2 .or. effective < maximum_predictors + 1) then
         out%info = 1
         return
      end if
      allocate(target(effective, variables), design(effective, maximum_predictors))
      allocate(threshold(effective))
      target = model%data(start:, :)
      design(:, 1) = 1.0_dp
      predictor = 1
      do lag = 1, maximum_order
         do source = 1, variables
            predictor = predictor + 1
            design(:, predictor) = model%data(start - lag:size(model%data, 1) - &
               lag, source)
         end do
      end do
      if (present(threshold_variable)) then
         if (size(threshold_variable) >= effective) then
            threshold = threshold_variable(:effective)
         else
            threshold = model%data(start - model%delay: &
               size(model%data, 1) - model%delay, model%threshold_component)
         end if
      else
         threshold = model%data(start - model%delay: &
            size(model%data, 1) - model%delay, model%threshold_component)
      end if
      out%model = model
      out%model%intercept = 0.0_dp
      out%model%ar = 0.0_dp
      out%model%residuals = 0.0_dp
      out%model%standardized_residuals = 0.0_dp
      out%model%covariance = 0.0_dp
      allocate(out%retained(maximum_predictors, variables, regimes))
      allocate(out%standard_error(maximum_predictors, variables, regimes))
      allocate(out%t_ratio(maximum_predictors, variables, regimes))
      out%retained = .false.
      out%standard_error = 0.0_dp
      out%t_ratio = 0.0_dp
      total_parameters = 0
      total_observations = 0
      allocate(candidate(maximum_predictors), keep(maximum_predictors))
      do regime = 1, regimes
         count_value = count([(threshold_regime(threshold(row), model%thresholds) == &
            regime, row=1, effective)])
         if (count_value <= variables + 1) then
            out%info = 2
            return
         end if
         allocate(indices(count_value))
         call collect_regime_indices(threshold, model%thresholds, regime, indices)
         candidate = .false.
         candidate(1) = model%include_mean(regime)
         do lag = 1, model%ar_order(regime)
            candidate(2 + (lag - 1)*variables:1 + lag*variables) = .true.
         end do
         allocate(predictors(count(candidate)))
         predictors = pack([(predictor, predictor=1, maximum_predictors)], candidate)
         do equation = 1, variables
            selected_design = design(indices, predictors)
            call ols_fit(selected_design, target(indices, equation), beta, &
               standard_error, residual, rss, status)
            if (status /= 0) then
               out%info = 3
               return
            end if
            keep = .false.
            do predictor = 1, size(predictors)
               out%standard_error(predictors(predictor), equation, regime) = &
                  standard_error(predictor)
               if (standard_error(predictor) > 0.0_dp) then
                  out%t_ratio(predictors(predictor), equation, regime) = &
                     beta(predictor)/standard_error(predictor)
               end if
               keep(predictors(predictor)) = abs(out%t_ratio( &
                  predictors(predictor), equation, regime)) >= out%threshold
            end do
            retained_count = count(keep)
            out%retained(:, equation, regime) = keep
            total_parameters = total_parameters + retained_count
            if (retained_count > 0) then
               deallocate(predictors)
               allocate(predictors(retained_count))
               predictors = pack([(predictor, predictor=1, maximum_predictors)], keep)
               selected_design = design(indices, predictors)
               call ols_fit(selected_design, target(indices, equation), beta, &
                  standard_error, residual, rss, status)
               if (status /= 0) then
                  out%info = 4
                  return
               end if
               do row = 1, retained_count
                  predictor = predictors(row)
                  if (predictor == 1) then
                     out%model%intercept(equation, regime) = beta(row)
                  else
                     lag = (predictor - 2)/variables + 1
                     source = modulo(predictor - 2, variables) + 1
                     out%model%ar(equation, source, lag, regime) = beta(row)
                  end if
               end do
            else
               residual = target(indices, equation)
            end if
            out%model%residuals(indices, equation) = residual
            if (allocated(predictors)) deallocate(predictors)
            allocate(predictors(count(candidate)))
            predictors = pack([(predictor, predictor=1, maximum_predictors)], candidate)
         end do
         covariance = matmul(transpose(out%model%residuals(indices, :)), &
            out%model%residuals(indices, :))/real(count_value, dp)
         out%model%covariance(:, :, regime) = covariance
         call symmetric_inverse_root(covariance, inverse_root, log_determinant, status)
         if (status /= 0) then
            out%info = 5
            return
         end if
         out%model%standardized_residuals(indices, :) = &
            matmul(out%model%residuals(indices, :), inverse_root)
         out%model%regime_observations(regime) = count_value
         total_observations = total_observations + count_value
         deallocate(indices, predictors)
      end do
      allocate(pooled(variables, variables), inverse(variables, variables))
      pooled = 0.0_dp
      do regime = 1, regimes
         pooled = pooled + real(out%model%regime_observations(regime), dp)* &
            out%model%covariance(:, :, regime)
      end do
      pooled = pooled/real(total_observations, dp)
      call inverse_logdet(pooled, inverse, log_determinant, status, 1.0e-12_dp)
      if (status /= 0) then
         out%info = 6
         return
      end if
      out%aic = real(total_observations, dp)*log_determinant + &
         2.0_dp*real(total_parameters, dp)
      out%bic = real(total_observations, dp)*log_determinant + &
         log(real(total_observations, dp))*real(total_parameters, dp)
      out%hq = real(total_observations, dp)*log_determinant + &
         2.0_dp*log(log(real(total_observations, dp)))*real(total_parameters, dp)
      out%model%aic = out%aic
   end function nts_mtar_refine

   pure function nts_tsay_test(series, ar_order) result(out)
      !! Apply Tsay's quadratic arranged-regression test for AR nonlinearity.
      real(dp), intent(in) :: series(:) !! Univariate time series.
      integer, intent(in), optional :: ar_order !! Null-model AR order.
      type(nts_nonlinearity_test_t) :: out
      real(dp), allocatable :: centered(:), target(:), design(:, :)
      real(dp), allocatable :: quadratic(:, :), residualized(:, :)
      real(dp), allocatable :: beta(:), standard_error(:), residual(:)
      real(dp), allocatable :: auxiliary_residual(:)
      real(dp) :: rss_null, rss_alternative, series_mean
      integer :: order, effective, terms, row, first, second, column, status

      order = 1
      if (present(ar_order)) order = ar_order
      if (order < 1 .or. size(series) <= 2*order + 2 .or. &
         .not. all(ieee_is_finite(series))) then
         out%info = 1
         return
      end if
      effective = size(series) - order
      terms = order*(order + 1)/2
      out%numerator_df = terms
      out%denominator_df = size(series) - order - terms - 1
      if (out%denominator_df <= 0) then
         out%info = 1
         return
      end if
      series_mean = sum(series)/real(size(series), dp)
      centered = series - series_mean
      allocate(target(effective), design(effective, order + 1))
      target = series(order + 1:)
      design(:, 1) = 1.0_dp
      do first = 1, order
         design(:, first + 1) = centered(order + 1 - first:size(series) - first)
      end do
      call ols_fit(design, target, beta, standard_error, residual, rss_null, status)
      if (status /= 0) then
         out%info = 2
         return
      end if
      allocate(quadratic(effective, terms), residualized(effective, terms))
      column = 0
      do first = 1, order
         do second = 1, first
            column = column + 1
            quadratic(:, column) = centered( &
               order + 1 - first:size(series) - first)*centered( &
               order + 1 - second:size(series) - second)
         end do
      end do
      do column = 1, terms
         call ols_fit(design, quadratic(:, column), beta, standard_error, &
            auxiliary_residual, rss_alternative, status)
         if (status /= 0) then
            out%info = 3
            return
         end if
         residualized(:, column) = auxiliary_residual
      end do
      call ols_fit(residualized, residual, beta, standard_error, &
         auxiliary_residual, rss_alternative, status)
      if (status /= 0 .or. rss_alternative <= 0.0_dp) then
         out%info = 4
         return
      end if
      out%statistic = (max(0.0_dp, rss_null - rss_alternative)/ &
         real(out%numerator_df, dp))/(rss_alternative/ &
         real(out%denominator_df, dp))
      out%p_value = f_upper_probability(out%statistic, out%numerator_df, &
         out%denominator_df)
   end function nts_tsay_test

   pure function nts_tar_backtest_draws(model, origin, horizon, normal_draws) &
      result(out)
      !! Backtest a SETAR model from supplied forecast simulation draws.
      type(nts_tar_model_t), intent(in) :: model !! Fitted SETAR specification.
      integer, intent(in) :: origin !! First rolling forecast origin.
      integer, intent(in) :: horizon !! Maximum forecast horizon.
      real(dp), intent(in) :: normal_draws(:, :, :) !! Horizon-by-draw-by-origin normals.
      type(nts_tar_backtest_t) :: out
      type(nts_tar_model_t) :: fitted
      type(nts_tar_forecast_t) :: forecast
      real(dp) :: sum_squared, sum_absolute, sum_error
      integer :: observations, origins, iterations, origin_index, time
      integer :: step, available, regime, count_value

      observations = size(model%data)
      origins = observations - origin
      iterations = size(normal_draws, 2)
      if (model%info /= 0 .or. .not. allocated(model%data) .or. &
         origin < max(maxval(model%ar_order), model%delay) .or. &
         origin >= observations .or. horizon < 1 .or. iterations < 1 .or. &
         any(shape(normal_draws) /= [horizon, iterations, origins]) .or. &
         .not. all(ieee_is_finite(normal_draws))) then
         out%info = 1
         return
      end if
      allocate(out%error(origins, horizon), out%state(origins))
      out%error = 0.0_dp
      do origin_index = 1, origins
         time = origin + origin_index - 1
         out%state(origin_index) = threshold_regime( &
            model%data(time - model%delay), model%thresholds)
         fitted = nts_utar_estimate(model%data(:time), model%ar_order, &
            model%thresholds, model%delay, include_mean=model%include_mean)
         if (fitted%info /= 0) then
            out%info = 10 + fitted%info
            return
         end if
         forecast = nts_utar_forecast_draws(fitted, time, &
            normal_draws(:, :, origin_index))
         if (forecast%info /= 0) then
            out%info = 20 + forecast%info
            return
         end if
         available = min(horizon, observations - time)
         out%error(origin_index, :available) = model%data(time + 1:time + &
            available) - forecast%mean(:available)
      end do
      allocate(out%rmse(horizon), out%mae(horizon), out%bias(horizon))
      do step = 1, horizon
         available = origins - step + 1
         sum_squared = sum(out%error(:available, step)**2)
         sum_absolute = sum(abs(out%error(:available, step)))
         sum_error = sum(out%error(:available, step))
         out%rmse(step) = sqrt(sum_squared/real(available, dp))
         out%mae(step) = sum_absolute/real(available, dp)
         out%bias(step) = sum_error/real(available, dp)
      end do
      allocate(out%regime_rmse(size(model%ar_order), horizon))
      allocate(out%regime_mae(size(model%ar_order), horizon))
      allocate(out%regime_bias(size(model%ar_order), horizon))
      allocate(out%regime_count(size(model%ar_order), horizon))
      out%regime_rmse = 0.0_dp
      out%regime_mae = 0.0_dp
      out%regime_bias = 0.0_dp
      out%regime_count = 0
      do regime = 1, size(model%ar_order)
         do step = 1, horizon
            available = origins - step + 1
            count_value = count(out%state(:available) == regime)
            out%regime_count(regime, step) = count_value
            if (count_value < 1) cycle
            sum_squared = sum(merge(out%error(:available, step)**2, 0.0_dp, &
               out%state(:available) == regime))
            sum_absolute = sum(merge(abs(out%error(:available, step)), 0.0_dp, &
               out%state(:available) == regime))
            sum_error = sum(merge(out%error(:available, step), 0.0_dp, &
               out%state(:available) == regime))
            out%regime_rmse(regime, step) = &
               sqrt(sum_squared/real(count_value, dp))
            out%regime_mae(regime, step) = sum_absolute/real(count_value, dp)
            out%regime_bias(regime, step) = sum_error/real(count_value, dp)
         end do
      end do
      out%origin = origin
      out%horizon = horizon
   end function nts_tar_backtest_draws

   function nts_tar_backtest(model, origin, horizon, iterations) result(out)
      !! Backtest a SETAR model using the shared random stream.
      type(nts_tar_model_t), intent(in) :: model !! Fitted SETAR specification.
      integer, intent(in) :: origin !! First rolling forecast origin.
      integer, intent(in) :: horizon !! Maximum forecast horizon.
      integer, intent(in) :: iterations !! Simulation draws per forecast origin.
      type(nts_tar_backtest_t) :: out
      real(dp), allocatable :: normal_draws(:, :, :)
      integer :: origins, origin_index

      if (.not. allocated(model%data) .or. origin < 1 .or. &
         origin >= size(model%data) .or. horizon < 1 .or. iterations < 1) then
         out%info = 1
         return
      end if
      origins = size(model%data) - origin
      allocate(normal_draws(horizon, iterations, origins))
      do origin_index = 1, origins
         call random_standard_normal_matrix(normal_draws(:, :, origin_index))
      end do
      out = nts_tar_backtest_draws(model, origin, horizon, normal_draws)
   end function nts_tar_backtest

   pure function nts_rank_portmanteau(series, maximum_lag) result(out)
      !! Compute NTS rank-based cumulative portmanteau statistics.
      real(dp), intent(in) :: series(:) !! Univariate observations.
      integer, intent(in), optional :: maximum_lag !! Largest tested lag.
      type(nts_rank_portmanteau_t) :: out
      real(dp), allocatable :: rank_value(:), centered_rank(:)
      real(dp) :: rank_mean, sum_square, variance_denominator, leading, cumulative
      integer :: lag_limit, observations, lag

      observations = size(series)
      lag_limit = 10
      if (present(maximum_lag)) lag_limit = maximum_lag
      if (observations < 3 .or. lag_limit < 1 .or. lag_limit >= observations .or. &
         .not. all(ieee_is_finite(series))) then
         out%info = 1
         return
      end if
      call average_ranks(series, rank_value)
      rank_mean = 0.5_dp*real(observations + 1, dp)
      centered_rank = rank_value - rank_mean
      sum_square = real(observations*(observations**2 - 1), dp)/12.0_dp
      variance_denominator = 5.0_dp*real((observations - 1)**2, dp)* &
         real(observations**2, dp)*real(observations + 1, dp)
      leading = 5.0_dp*real(observations, dp)**4
      allocate(out%correlation(lag_limit), out%expected_correlation(lag_limit))
      allocate(out%correlation_variance(lag_limit), out%statistic(lag_limit))
      allocate(out%p_value(lag_limit))
      cumulative = 0.0_dp
      do lag = 1, lag_limit
         out%correlation(lag) = dot_product(centered_rank(:observations - lag), &
            centered_rank(lag + 1:))/sum_square
         out%expected_correlation(lag) = -real(observations - lag, dp)/ &
            real(observations*(observations - 1), dp)
         out%correlation_variance(lag) = (leading - &
            real(5*lag + 9, dp)*real(observations, dp)**3 + &
            9.0_dp*real(lag - 2, dp)*real(observations, dp)**2 + &
            2.0_dp*real(lag*(5*lag + 8)*observations, dp) + &
            16.0_dp*real(lag**2, dp))/variance_denominator
         cumulative = cumulative + (out%correlation(lag) - &
            out%expected_correlation(lag))**2/out%correlation_variance(lag)
         out%statistic(lag) = cumulative
         out%p_value(lag) = regularized_gamma_q(0.5_dp*real(lag, dp), &
            0.5_dp*cumulative)
      end do
   end function nts_rank_portmanteau

   pure function nts_quadratic_f_test(series, ar_order, t_threshold) result(out)
      !! Test added quadratic terms selected by null-model AR t ratios.
      real(dp), intent(in) :: series(:) !! Univariate observations.
      integer, intent(in) :: ar_order !! Null-model AR order.
      real(dp), intent(in), optional :: t_threshold !! Absolute t-ratio selection threshold.
      type(nts_nonlinearity_test_t) :: out
      real(dp), allocatable :: target(:), null_design(:, :), alternative(:, :)
      real(dp), allocatable :: beta(:), standard_error(:), residual(:)
      integer, allocatable :: selected(:)
      real(dp) :: threshold, null_rss, alternative_rss
      integer :: effective, selected_count, terms, first, second, column, status

      threshold = 0.0_dp
      if (present(t_threshold)) threshold = t_threshold
      if (ar_order < 1 .or. threshold < 0.0_dp .or. &
         size(series) <= 2*ar_order + 3 .or. &
         .not. all(ieee_is_finite(series))) then
         out%info = 1
         return
      end if
      effective = size(series) - ar_order
      allocate(target(effective), null_design(effective, ar_order + 1))
      target = series(ar_order + 1:)
      null_design(:, 1) = 1.0_dp
      do first = 1, ar_order
         null_design(:, first + 1) = &
            series(ar_order + 1 - first:size(series) - first)
      end do
      call ols_fit(null_design, target, beta, standard_error, residual, &
         null_rss, status)
      if (status /= 0) then
         out%info = 2
         return
      end if
      allocate(selected(ar_order))
      selected_count = 0
      do first = 1, ar_order
         if (standard_error(first + 1) > 0.0_dp) then
            if (abs(beta(first + 1)/standard_error(first + 1)) > threshold) then
               selected_count = selected_count + 1
               selected(selected_count) = first
            end if
         end if
      end do
      if (selected_count == 0) then
         selected_count = ar_order
         selected = [(first, first=1, ar_order)]
      end if
      terms = selected_count*(selected_count + 1)/2
      allocate(alternative(effective, ar_order + 1 + terms))
      alternative(:, :ar_order + 1) = null_design
      column = ar_order + 1
      do first = 1, selected_count
         do second = 1, first
            column = column + 1
            alternative(:, column) = null_design(:, selected(first) + 1)* &
               null_design(:, selected(second) + 1)
         end do
      end do
      call ols_fit(alternative, target, beta, standard_error, residual, &
         alternative_rss, status)
      out%numerator_df = terms
      out%denominator_df = effective - size(alternative, 2)
      if (status /= 0 .or. alternative_rss <= 0.0_dp .or. &
         out%denominator_df <= 0) then
         out%info = 3
         return
      end if
      out%statistic = (max(0.0_dp, null_rss - alternative_rss)/ &
         real(terms, dp))/(alternative_rss/real(out%denominator_df, dp))
      out%p_value = f_upper_probability(out%statistic, terms, &
         out%denominator_df)
   end function nts_quadratic_f_test

   pure function nts_prnd_test(series, maximum_lag, ar_order, ma_order) result(out)
      !! Compute the Peña-Rodríguez determinant lack-of-fit statistic.
      real(dp), intent(in) :: series(:) !! Residual or observation series.
      integer, intent(in), optional :: maximum_lag !! Maximum autocorrelation lag.
      integer, intent(in), optional :: ar_order !! Fitted AR order adjustment.
      integer, intent(in), optional :: ma_order !! Fitted MA order adjustment.
      type(nts_prnd_test_t) :: out
      real(dp), allocatable :: centered(:), autocorrelation(:), matrix(:, :)
      real(dp), allocatable :: inverse(:, :)
      real(dp) :: denominator, alpha, beta_value, first_term, denominator_term
      real(dp) :: exponent, lambda, log_determinant, scale
      integer :: lag_limit, p, q, pq, observations, lag, row, column, status

      lag_limit = 10
      if (present(maximum_lag)) lag_limit = maximum_lag
      p = 0
      if (present(ar_order)) p = ar_order
      q = 0
      if (present(ma_order)) q = ma_order
      pq = p + q
      observations = size(series)
      if (lag_limit < 1 .or. observations <= lag_limit .or. p < 0 .or. q < 0 .or. &
         lag_limit <= 2*pq .or. .not. all(ieee_is_finite(series))) then
         out%info = 1
         return
      end if
      centered = series - sum(series)/real(observations, dp)
      denominator = dot_product(centered, centered)
      if (denominator <= tiny(1.0_dp)) then
         out%info = 1
         return
      end if
      allocate(autocorrelation(0:lag_limit), matrix(lag_limit + 1, lag_limit + 1))
      autocorrelation(0) = 1.0_dp
      do lag = 1, lag_limit
         autocorrelation(lag) = real(observations + 2, dp)* &
            dot_product(centered(:observations - lag), centered(lag + 1:))/ &
            (real(observations - lag, dp)*denominator)
      end do
      do row = 1, lag_limit + 1
         do column = 1, lag_limit + 1
            matrix(row, column) = autocorrelation(abs(row - column))
         end do
      end do
      allocate(inverse(lag_limit + 1, lag_limit + 1))
      call inverse_logdet(matrix, inverse, log_determinant, status, 1.0e-12_dp)
      if (status /= 0) then
         out%info = 2
         return
      end if
      alpha = 3.0_dp*real(lag_limit + 1, dp)*real(lag_limit - 2*pq, dp)**2/ &
         (4.0_dp*(real(lag_limit*(2*lag_limit + 1), dp) - &
         6.0_dp*real((lag_limit + 1)*pq, dp)))
      beta_value = 3.0_dp*real(lag_limit + 1, dp)* &
         real(lag_limit - 2*pq, dp)/(2.0_dp*real(lag_limit*(2*lag_limit + 1), dp) - &
         12.0_dp*real((lag_limit + 1)*pq, dp))
      first_term = 2.0_dp*(0.5_dp*real(lag_limit, dp) - real(pq, dp))* &
         (real(lag_limit**2, dp)/(4.0_dp*real(lag_limit + 1, dp)) - real(pq, dp))
      denominator_term = 3.0_dp*(real(lag_limit*(2*lag_limit + 1), dp)/ &
         (6.0_dp*real(lag_limit + 1, dp)) - real(pq, dp))**2
      exponent = 1.0_dp - first_term/denominator_term
      if (alpha <= 0.0_dp .or. beta_value <= 0.0_dp .or. exponent <= 0.0_dp) then
         out%info = 3
         return
      end if
      lambda = 1.0_dp/exponent
      out%determinant_statistic = -real(observations, dp)/ &
         real(lag_limit + 1, dp)*log_determinant
      scale = alpha/beta_value
      out%statistic = scale**(-exponent)*(lambda/sqrt(alpha))* &
         (out%determinant_statistic**exponent - scale**exponent* &
         (1.0_dp - (lambda - 1.0_dp)/(2.0_dp*alpha*lambda**2)))
      out%p_value = erfc(abs(out%statistic)/sqrt(2.0_dp))
      out%lag = lag_limit
   end function nts_prnd_test

   pure function nts_tvar_filter_smooth(series, lags, log_variance, include_mean) &
      result(out)
      !! Filter and smooth a time-varying AR model for supplied log variances.
      real(dp), intent(in) :: series(:) !! Univariate observations.
      integer, intent(in) :: lags(:) !! Included positive AR lags.
      real(dp), intent(in) :: log_variance(:) !! Observation then coefficient log variances.
      logical, intent(in), optional :: include_mean !! Whether the intercept varies over time.
      type(nts_tvar_fit_t) :: out
      type(ssm_model_t) :: state_model
      type(kfs_filter_t) :: filtered
      type(kfs_smoother_t) :: smoothed
      real(dp), allocatable :: response(:), design(:, :), beta(:)
      real(dp), allocatable :: standard_error(:), residual(:)
      real(dp) :: rss
      integer :: status, coefficient, time
      logical :: selected_mean

      selected_mean = .true.
      if (present(include_mean)) selected_mean = include_mean
      call build_lagged_regression(series, lags, selected_mean, response, &
         design, status)
      if (status /= 0 .or. size(log_variance) /= size(design, 2) + 1 .or. &
         .not. all(ieee_is_finite(log_variance))) then
         out%info = 1
         return
      end if
      call ols_fit(design, response, beta, standard_error, residual, rss, status)
      if (status /= 0) then
         out%info = 2
         return
      end if
      state_model = make_tvar_model(response, design, log_variance, beta)
      filtered = kfs_filter(state_model)
      if (filtered%info /= 0) then
         out%info = 10 + filtered%info
         return
      end if
      smoothed = kfs_smooth(state_model, filtered)
      if (smoothed%info /= 0) then
         out%info = 20 + smoothed%info
         return
      end if
      allocate(out%filtered_coefficients(size(response), size(design, 2)))
      allocate(out%smoothed_coefficients(size(response), size(design, 2)))
      allocate(out%filtered_covariance(size(design, 2), size(design, 2), &
         size(response)))
      allocate(out%smoothed_covariance(size(design, 2), size(design, 2), &
         size(response)))
      do time = 1, size(response)
         out%filtered_coefficients(time, :) = filtered%a_filt(:, time)
         out%smoothed_coefficients(time, :) = smoothed%state(:, time)
         out%filtered_covariance(:, :, time) = filtered%p_filt(:, :, time)
         out%smoothed_covariance(:, :, time) = smoothed%covariance(:, :, time)
      end do
      out%log_variance = log_variance
      allocate(out%state_variance(size(design, 2)))
      do coefficient = 1, size(design, 2)
         out%state_variance(coefficient) = exp(log_variance(coefficient + 1))
      end do
      out%observation_variance = exp(log_variance(1))
      out%log_likelihood = filtered%log_likelihood
      out%lags = lags
      out%include_mean = selected_mean
   end function nts_tvar_filter_smooth

   pure function nts_tvar_fit(series, lags, include_mean, max_iterations, &
      tolerance) result(out)
      !! Estimate time-varying AR observation and state variances by likelihood.
      real(dp), intent(in) :: series(:) !! Univariate observations.
      integer, intent(in) :: lags(:) !! Included positive AR lags.
      logical, intent(in), optional :: include_mean !! Whether the intercept varies over time.
      integer, intent(in), optional :: max_iterations !! Maximum BFGS iterations.
      real(dp), intent(in), optional :: tolerance !! Gradient convergence tolerance.
      type(nts_tvar_fit_t) :: out
      type(optimization_result_t) :: optimization
      real(dp), allocatable :: response(:), design(:, :), initial(:), beta(:)
      real(dp), allocatable :: standard_error(:), residual(:)
      real(dp) :: rss, selected_tolerance
      integer :: status, limit
      logical :: selected_mean

      selected_mean = .true.
      if (present(include_mean)) selected_mean = include_mean
      limit = 200
      if (present(max_iterations)) limit = max_iterations
      selected_tolerance = 1.0e-6_dp
      if (present(tolerance)) selected_tolerance = tolerance
      call build_lagged_regression(series, lags, selected_mean, response, &
         design, status)
      if (status /= 0 .or. limit < 1 .or. selected_tolerance <= 0.0_dp) then
         out%info = 1
         return
      end if
      call ols_fit(design, response, beta, standard_error, residual, rss, status)
      if (status /= 0) then
         out%info = 2
         return
      end if
      allocate(initial(size(design, 2) + 1))
      initial = -0.5_dp
      optimization = bfgs_minimize_fd(objective, initial, limit, &
         selected_tolerance)
      if (.not. allocated(optimization%parameters) .or. &
         .not. ieee_is_finite(optimization%objective)) then
         out%info = 3
         return
      end if
      out = nts_tvar_filter_smooth(series, lags, optimization%parameters, &
         selected_mean)
      if (out%info /= 0) return
      out%converged = optimization%converged
      out%iterations = optimization%iterations

   contains

      pure function objective(parameters) result(value)
         !! Evaluate negative time-varying AR Gaussian log likelihood.
         real(dp), intent(in) :: parameters(:) !! Observation and state log variances.
         real(dp) :: value
         type(ssm_model_t) :: candidate_model
         type(kfs_filter_t) :: candidate_filter

         candidate_model = make_tvar_model(response, design, parameters, beta)
         candidate_filter = kfs_filter(candidate_model)
         if (candidate_filter%info /= 0) then
            value = huge(1.0_dp)
         else
            value = -candidate_filter%log_likelihood
         end if
      end function objective

   end function nts_tvar_fit

   pure function nts_rcar_fit(series, lags, include_mean, max_iterations, &
      tolerance) result(out)
      !! Estimate a Gaussian random-coefficient AR model by maximum likelihood.
      real(dp), intent(in) :: series(:) !! Univariate observations.
      integer, intent(in) :: lags(:) !! Included positive AR lags.
      logical, intent(in), optional :: include_mean !! Whether to remove and retain the sample mean.
      integer, intent(in), optional :: max_iterations !! Maximum BFGS iterations.
      real(dp), intent(in), optional :: tolerance !! Gradient convergence tolerance.
      type(nts_rcar_fit_t) :: out
      type(optimization_result_t) :: optimization
      real(dp), allocatable :: centered(:), response(:), design(:, :), initial(:)
      real(dp), allocatable :: beta(:), ols_standard_error(:), ols_residual(:)
      real(dp), allocatable :: hessian(:, :), inverse(:, :), variance(:)
      real(dp) :: rss, selected_tolerance
      integer :: status, limit, coefficients
      logical :: selected_mean

      selected_mean = .true.
      if (present(include_mean)) selected_mean = include_mean
      limit = 200
      if (present(max_iterations)) limit = max_iterations
      selected_tolerance = 1.0e-6_dp
      if (present(tolerance)) selected_tolerance = tolerance
      if (size(lags) < 1 .or. limit < 1 .or. selected_tolerance <= 0.0_dp .or. &
         .not. all(ieee_is_finite(series))) then
         out%info = 1
         return
      end if
      centered = series
      if (selected_mean) then
         out%sample_mean = sum(series)/real(size(series), dp)
         centered = centered - out%sample_mean
      end if
      call build_lagged_regression(centered, lags, .false., response, design, &
         status)
      if (status /= 0) then
         out%info = 1
         return
      end if
      call ols_fit(design, response, beta, ols_standard_error, ols_residual, &
         rss, status)
      if (status /= 0 .or. rss <= 0.0_dp) then
         out%info = 2
         return
      end if
      coefficients = size(design, 2)
      allocate(initial(2*coefficients + 1))
      initial(:coefficients) = beta
      initial(coefficients + 1:2*coefficients) = -3.0_dp
      initial(2*coefficients + 1) = log(rss/real(size(response), dp))
      optimization = bfgs_minimize_fd(objective, initial, limit, &
         selected_tolerance)
      if (.not. allocated(optimization%parameters) .or. &
         .not. ieee_is_finite(optimization%objective)) then
         out%info = 3
         return
      end if
      out%coefficients = optimization%parameters(:coefficients)
      out%coefficient_variance = exp(optimization%parameters( &
         coefficients + 1:2*coefficients))
      out%innovation_variance = exp(optimization%parameters(2*coefficients + 1))
      out%residuals = response - matmul(design, out%coefficients)
      variance = out%innovation_variance + &
         matmul(design**2, out%coefficient_variance)
      out%standardized_residuals = out%residuals/sqrt(variance)
      hessian = finite_difference_hessian(objective, optimization%parameters)
      call invert_matrix(hessian, inverse, status)
      if (status == 0 .and. all(ieee_is_finite(inverse))) then
         out%standard_error = sqrt(max(0.0_dp, diagonal_values(inverse)))
      else
         allocate(out%standard_error(0))
      end if
      out%lags = lags
      out%log_likelihood = -optimization%objective
      out%include_mean = selected_mean
      out%converged = optimization%converged
      out%iterations = optimization%iterations

   contains

      pure function objective(parameters) result(value)
         !! Evaluate negative random-coefficient AR Gaussian log likelihood.
         real(dp), intent(in) :: parameters(:) !! Coefficients and log variances.
         real(dp) :: value
         real(dp), allocatable :: conditional_variance(:), error(:)

         conditional_variance = exp(parameters(2*coefficients + 1)) + &
            matmul(design**2, exp(parameters(coefficients + 1:2*coefficients)))
         error = response - matmul(design, parameters(:coefficients))
         if (any(conditional_variance <= tiny(1.0_dp)) .or. &
            .not. all(ieee_is_finite(conditional_variance))) then
            value = huge(1.0_dp)
         else
            value = 0.5_dp*sum(log(2.0_dp*acos(-1.0_dp)) + &
               log(conditional_variance) + error**2/conditional_variance)
         end if
      end function objective

   end function nts_rcar_fit

   pure function nts_acmx_fit(counts, ar_order, mean_order, family, exogenous, &
      initial_parameters, max_iterations, tolerance) result(out)
      !! Estimate an autoregressive conditional-mean count model with covariates.
      real(dp), intent(in) :: counts(:) !! Nonnegative integer observations.
      integer, intent(in) :: ar_order !! Number of lagged-count terms.
      integer, intent(in) :: mean_order !! Number of lagged conditional-mean terms.
      character(len=*), intent(in) :: family !! Poisson, negative binomial, or double Poisson.
      real(dp), intent(in), optional :: exogenous(:, :) !! Observation-aligned covariates.
      real(dp), intent(in), optional :: initial_parameters(:) !! Natural-scale initial parameters.
      integer, intent(in), optional :: max_iterations !! Maximum BFGS iterations.
      real(dp), intent(in), optional :: tolerance !! Gradient convergence tolerance.
      type(nts_acmx_fit_t) :: out
      type(optimization_result_t) :: optimization
      real(dp), allocatable :: initial(:), natural(:), hessian(:, :), inverse(:, :)
      real(dp), allocatable :: jacobian(:, :), natural_covariance(:, :)
      real(dp), allocatable :: shifted(:), plus(:), minus(:), variance(:)
      real(dp) :: selected_tolerance, step, sample_mean, sample_variance
      integer :: covariates, parameter_count, unconstrained_count, limit
      integer :: status, parameter
      logical :: has_dispersion
      character(len=:), allocatable :: selected_family

      selected_family = trim(adjustl(lowercase(family)))
      select case (selected_family)
      case ('po', 'poisson')
         selected_family = 'poisson'
         has_dispersion = .false.
      case ('nb', 'negative_binomial')
         selected_family = 'negative_binomial'
         has_dispersion = .true.
      case ('dp', 'double_poisson')
         selected_family = 'double_poisson'
         has_dispersion = .true.
      case default
         out%info = 1
         return
      end select
      covariates = 0
      if (present(exogenous)) then
         if (size(exogenous, 1) /= size(counts) .or. &
            .not. all(ieee_is_finite(exogenous))) then
            out%info = 1
            return
         end if
         covariates = size(exogenous, 2)
      end if
      limit = 300
      if (present(max_iterations)) limit = max_iterations
      selected_tolerance = 1.0e-6_dp
      if (present(tolerance)) selected_tolerance = tolerance
      if (size(counts) <= ar_order + mean_order + covariates + 2 .or. &
         ar_order < 0 .or. mean_order < 0 .or. limit < 1 .or. &
         selected_tolerance <= 0.0_dp .or. any(counts < 0.0_dp) .or. &
         any(abs(counts - real(nint(counts), dp)) > 1.0e-10_dp) .or. &
         .not. all(ieee_is_finite(counts))) then
         out%info = 1
         return
      end if
      parameter_count = covariates + 1 + ar_order + mean_order + &
         merge(1, 0, has_dispersion)
      unconstrained_count = parameter_count
      allocate(natural(parameter_count), initial(unconstrained_count))
      sample_mean = sum(counts)/real(size(counts), dp)
      sample_variance = sum((counts - sample_mean)**2)/ &
         real(max(1, size(counts) - 1), dp)
      natural = 0.0_dp
      if (covariates > 0) natural(:covariates) = 0.0_dp
      natural(covariates + 1) = max(0.1_dp, 0.4_dp*sample_mean)
      if (ar_order > 0) natural(covariates + 2:covariates + 1 + ar_order) = &
         0.05_dp/real(ar_order, dp)
      if (mean_order > 0) natural(covariates + 2 + ar_order: &
         covariates + 1 + ar_order + mean_order) = 0.5_dp/real(mean_order, dp)
      if (has_dispersion) then
         if (selected_family == 'negative_binomial') then
            natural(parameter_count) = max(0.5_dp, min(100.0_dp, &
               sample_mean**2/max(0.1_dp, sample_variance - sample_mean)))
         else
            natural(parameter_count) = max(1.0e-3_dp, min(10.0_dp, &
               sample_mean/max(sample_variance, 1.0e-3_dp)))
         end if
      end if
      if (present(initial_parameters)) then
         if (size(initial_parameters) /= parameter_count .or. &
            .not. valid_acmx_natural(initial_parameters, covariates, ar_order, &
            mean_order, has_dispersion)) then
            out%info = 2
            return
         end if
         natural = initial_parameters
      end if
      call acmx_unconstrain(natural, covariates, ar_order, mean_order, &
         has_dispersion, initial)
      optimization = bfgs_minimize_fd(objective, initial, limit, &
         selected_tolerance)
      if (.not. allocated(optimization%parameters) .or. &
         .not. ieee_is_finite(optimization%objective)) then
         out%info = 3
         return
      end if
      call acmx_natural(optimization%parameters, covariates, ar_order, &
         mean_order, has_dispersion, natural)
      out%parameter = natural
      if (covariates > 0) then
         out%exogenous_coefficients = natural(:covariates)
      else
         allocate(out%exogenous_coefficients(0))
      end if
      out%omega = natural(covariates + 1)
      out%alpha = natural(covariates + 2:covariates + 1 + ar_order)
      out%gamma = natural(covariates + 2 + ar_order: &
         covariates + 1 + ar_order + mean_order)
      if (has_dispersion) out%dispersion = natural(parameter_count)
      call acmx_conditional_mean(counts, ar_order, mean_order, exogenous, &
         out%exogenous_coefficients, out%omega, out%alpha, out%gamma, &
         out%conditional_mean)
      out%residuals = counts(ar_order + 1:) - out%conditional_mean
      allocate(variance(size(out%conditional_mean)))
      select case (selected_family)
      case ('poisson')
         variance = out%conditional_mean
      case ('negative_binomial')
         variance = out%conditional_mean + out%conditional_mean**2/out%dispersion
      case ('double_poisson')
         variance = out%conditional_mean/out%dispersion
      end select
      out%standardized_residuals = out%residuals/sqrt(max(variance, tiny(1.0_dp)))
      hessian = finite_difference_hessian(objective, optimization%parameters)
      call invert_matrix(hessian, inverse, status)
      if (status == 0 .and. all(ieee_is_finite(inverse))) then
         allocate(jacobian(parameter_count, unconstrained_count))
         allocate(shifted(unconstrained_count), plus(parameter_count), &
            minus(parameter_count))
         do parameter = 1, unconstrained_count
            shifted = optimization%parameters
            step = 1.0e-5_dp*max(1.0_dp, abs(shifted(parameter)))
            shifted(parameter) = optimization%parameters(parameter) + step
            call acmx_natural(shifted, covariates, ar_order, mean_order, &
               has_dispersion, plus)
            shifted(parameter) = optimization%parameters(parameter) - step
            call acmx_natural(shifted, covariates, ar_order, mean_order, &
               has_dispersion, minus)
            jacobian(:, parameter) = (plus - minus)/(2.0_dp*step)
         end do
         natural_covariance = matmul(matmul(jacobian, inverse), transpose(jacobian))
         out%standard_error = sqrt(max(0.0_dp, diagonal_values(natural_covariance)))
      else
         allocate(out%standard_error(0))
      end if
      out%family = selected_family
      out%ar_order = ar_order
      out%mean_order = mean_order
      out%log_likelihood = -optimization%objective
      out%aic = -2.0_dp*out%log_likelihood + 2.0_dp*real(parameter_count, dp)
      out%bic = -2.0_dp*out%log_likelihood + &
         log(real(size(out%conditional_mean), dp))*real(parameter_count, dp)
      out%converged = optimization%converged
      out%iterations = optimization%iterations

   contains

      pure function objective(parameters) result(value)
         !! Evaluate negative ACMx conditional log likelihood.
         real(dp), intent(in) :: parameters(:) !! Unconstrained ACMx parameters.
         real(dp) :: value
         real(dp), allocatable :: candidate(:), candidate_mean(:)

         call acmx_natural(parameters, covariates, ar_order, mean_order, &
            has_dispersion, candidate)
         call acmx_conditional_mean(counts, ar_order, mean_order, exogenous, &
            candidate(:covariates), candidate(covariates + 1), &
            candidate(covariates + 2:covariates + 1 + ar_order), &
            candidate(covariates + 2 + ar_order: &
            covariates + 1 + ar_order + mean_order), candidate_mean)
         value = acmx_negative_log_likelihood(counts(ar_order + 1:), &
            candidate_mean, selected_family, candidate(parameter_count))
      end function objective

   end function nts_acmx_fit

   pure subroutine symmetric_inverse_root(matrix, inverse_root, log_determinant, &
      info)
      !! Compute a symmetric inverse square root and log determinant.
      real(dp), intent(in) :: matrix(:, :) !! Symmetric positive-definite matrix.
      real(dp), allocatable, intent(out) :: inverse_root(:, :) !! Inverse square root.
      real(dp), intent(out) :: log_determinant !! Natural logarithm of determinant.
      integer, intent(out) :: info !! Status code; zero indicates success.
      real(dp), allocatable :: values(:), vectors(:, :), scales(:, :)
      real(dp) :: tolerance
      integer :: dimension, status

      dimension = size(matrix, 1)
      info = 0
      log_determinant = 0.0_dp
      if (dimension < 1 .or. size(matrix, 2) /= dimension) then
         info = 1
         return
      end if
      call symmetric_eigen(matrix, values, vectors, status)
      if (status /= 0) then
         info = 2
         return
      end if
      tolerance = 100.0_dp*epsilon(1.0_dp)*max(1.0_dp, maxval(abs(values)))
      if (any(values <= tolerance)) then
         info = 3
         return
      end if
      allocate(scales(dimension, dimension))
      scales = spread(1.0_dp/sqrt(values), 1, dimension)
      inverse_root = matmul(vectors*scales, transpose(vectors))
      log_determinant = sum(log(values))
   end subroutine symmetric_inverse_root

   pure logical function valid_acmx_natural(parameters, covariates, ar_order, &
      mean_order, has_dispersion) result(valid)
      !! Check a natural-scale ACMx parameter vector.
      real(dp), intent(in) :: parameters(:) !! Natural-scale parameters.
      integer, intent(in) :: covariates !! Number of exogenous coefficients.
      integer, intent(in) :: ar_order !! Number of lagged-count terms.
      integer, intent(in) :: mean_order !! Number of lagged-mean terms.
      logical, intent(in) :: has_dispersion !! Whether a dispersion parameter is present.
      integer :: dynamic_start, dynamic_end

      valid = all(ieee_is_finite(parameters))
      if (.not. valid) return
      dynamic_start = covariates + 2
      dynamic_end = covariates + 1 + ar_order + mean_order
      valid = parameters(covariates + 1) > 0.0_dp
      if (dynamic_end >= dynamic_start) valid = valid .and. &
         all(parameters(dynamic_start:dynamic_end) > 0.0_dp) .and. &
         sum(parameters(dynamic_start:dynamic_end)) < 1.0_dp
      if (has_dispersion) valid = valid .and. parameters(size(parameters)) > 0.0_dp
   end function valid_acmx_natural

   pure subroutine acmx_unconstrain(natural, covariates, ar_order, mean_order, &
      has_dispersion, parameters)
      !! Map natural ACMx parameters to unconstrained optimization coordinates.
      real(dp), intent(in) :: natural(:) !! Natural-scale parameters.
      integer, intent(in) :: covariates !! Number of exogenous coefficients.
      integer, intent(in) :: ar_order !! Number of lagged-count terms.
      integer, intent(in) :: mean_order !! Number of lagged-mean terms.
      logical, intent(in) :: has_dispersion !! Whether a dispersion parameter is present.
      real(dp), allocatable, intent(out) :: parameters(:) !! Unconstrained parameters.
      real(dp) :: remaining
      integer :: dynamic_count, dynamic_start, dynamic_end

      dynamic_count = ar_order + mean_order
      dynamic_start = covariates + 2
      dynamic_end = covariates + 1 + dynamic_count
      allocate(parameters(size(natural)))
      if (covariates > 0) parameters(:covariates) = natural(:covariates)
      parameters(covariates + 1) = log(natural(covariates + 1))
      if (dynamic_count > 0) then
         remaining = 1.0_dp - sum(natural(dynamic_start:dynamic_end))
         parameters(dynamic_start:dynamic_end) = &
            log(natural(dynamic_start:dynamic_end)/remaining)
      end if
      if (has_dispersion) parameters(size(parameters)) = log(natural(size(natural)))
   end subroutine acmx_unconstrain

   pure subroutine acmx_natural(parameters, covariates, ar_order, mean_order, &
      has_dispersion, natural)
      !! Map unconstrained ACMx coordinates to stable natural parameters.
      real(dp), intent(in) :: parameters(:) !! Unconstrained parameters.
      integer, intent(in) :: covariates !! Number of exogenous coefficients.
      integer, intent(in) :: ar_order !! Number of lagged-count terms.
      integer, intent(in) :: mean_order !! Number of lagged-mean terms.
      logical, intent(in) :: has_dispersion !! Whether a dispersion parameter is present.
      real(dp), allocatable, intent(out) :: natural(:) !! Natural-scale parameters.
      real(dp), allocatable :: exponential(:)
      real(dp) :: maximum_value, denominator
      integer :: dynamic_count, dynamic_start, dynamic_end

      dynamic_count = ar_order + mean_order
      dynamic_start = covariates + 2
      dynamic_end = covariates + 1 + dynamic_count
      allocate(natural(size(parameters)))
      natural = 0.0_dp
      if (covariates > 0) natural(:covariates) = parameters(:covariates)
      natural(covariates + 1) = exp(max(-40.0_dp, min(40.0_dp, &
         parameters(covariates + 1))))
      if (dynamic_count > 0) then
         maximum_value = max(0.0_dp, maxval(parameters(dynamic_start:dynamic_end)))
         exponential = exp(parameters(dynamic_start:dynamic_end) - maximum_value)
         denominator = exp(-maximum_value) + sum(exponential)
         natural(dynamic_start:dynamic_end) = exponential/denominator
      end if
      if (has_dispersion) natural(size(natural)) = exp(max(-40.0_dp, &
         min(40.0_dp, parameters(size(parameters)))))
   end subroutine acmx_natural

   pure subroutine acmx_conditional_mean(counts, ar_order, mean_order, exogenous, &
      exogenous_coefficients, omega, alpha, gamma, conditional_mean)
      !! Evaluate the ACMx conditional-mean recursion.
      real(dp), intent(in) :: counts(:) !! Complete count series.
      integer, intent(in) :: ar_order !! Number of lagged-count terms.
      integer, intent(in) :: mean_order !! Number of lagged-mean terms.
      real(dp), intent(in), optional :: exogenous(:, :) !! Observation-aligned covariates.
      real(dp), intent(in) :: exogenous_coefficients(:) !! Log-exposure coefficients.
      real(dp), intent(in) :: omega !! Positive conditional-mean intercept.
      real(dp), intent(in) :: alpha(:) !! Lagged-count coefficients.
      real(dp), intent(in) :: gamma(:) !! Lagged-mean coefficients.
      real(dp), allocatable, intent(out) :: conditional_mean(:) !! Aligned conditional means.
      real(dp), allocatable :: exposure(:), latent(:)
      real(dp) :: initial_mean, value
      integer :: observations, start, time, lag

      observations = size(counts)
      start = ar_order + 1
      allocate(exposure(observations), latent(observations))
      exposure = 1.0_dp
      if (present(exogenous)) exposure = exp(max(-40.0_dp, min(40.0_dp, &
         matmul(exogenous, exogenous_coefficients))))
      initial_mean = sum(counts/exposure)/real(observations, dp)
      latent = max(initial_mean, 1.0e-8_dp)
      allocate(conditional_mean(observations - ar_order))
      do time = start, observations
         value = omega
         do lag = 1, ar_order
            value = value + alpha(lag)*counts(time - lag)
         end do
         do lag = 1, mean_order
            if (time - lag >= 1) then
               value = value + gamma(lag)*latent(time - lag)
            else
               value = value + gamma(lag)*initial_mean
            end if
         end do
         latent(time) = max(value, 1.0e-12_dp)
         conditional_mean(time - ar_order) = max(1.0e-12_dp, &
            exposure(time)*latent(time))
      end do
   end subroutine acmx_conditional_mean

   pure real(dp) function acmx_negative_log_likelihood(counts, conditional_mean, &
      family, dispersion) result(value)
      !! Evaluate an ACMx negative conditional log likelihood.
      real(dp), intent(in) :: counts(:) !! Aligned count observations.
      real(dp), intent(in) :: conditional_mean(:) !! Positive conditional means.
      character(len=*), intent(in) :: family !! Selected count distribution.
      real(dp), intent(in) :: dispersion !! Positive dispersion when used.
      real(dp) :: log_probability, scaled_mean, correction
      integer :: observation

      value = 0.0_dp
      do observation = 1, size(counts)
         select case (family)
         case ('poisson')
            log_probability = counts(observation)*log( &
               conditional_mean(observation)) - conditional_mean(observation) - &
               log_gamma(counts(observation) + 1.0_dp)
         case ('negative_binomial')
            log_probability = log_gamma(counts(observation) + dispersion) - &
               log_gamma(dispersion) - log_gamma(counts(observation) + 1.0_dp) + &
               dispersion*log(dispersion/(dispersion + &
               conditional_mean(observation))) + counts(observation)*log( &
               conditional_mean(observation)/(dispersion + &
               conditional_mean(observation)))
         case ('double_poisson')
            scaled_mean = conditional_mean(observation)*dispersion
            correction = 1.0_dp + (1.0_dp - dispersion)/ &
               (12.0_dp*scaled_mean)*(1.0_dp + 1.0_dp/scaled_mean)
            correction = max(correction, tiny(1.0_dp))
            log_probability = -log(correction) + 0.5_dp*log(dispersion) - &
               scaled_mean
            if (counts(observation) > 0.0_dp) then
               log_probability = log_probability - &
                  log_gamma(counts(observation) + 1.0_dp) - counts(observation) + &
                  counts(observation)*log(counts(observation)) + &
                  dispersion*counts(observation)*(1.0_dp + &
                  log(conditional_mean(observation)) - log(counts(observation)))
            end if
         end select
         if (.not. ieee_is_finite(log_probability)) then
            value = huge(1.0_dp)
            return
         end if
         value = value - log_probability
      end do
   end function acmx_negative_log_likelihood

   pure subroutine build_lagged_regression(series, lags, include_mean, response, &
      design, info)
      !! Construct an aligned response and selected-lag regression matrix.
      real(dp), intent(in) :: series(:) !! Complete univariate series.
      integer, intent(in) :: lags(:) !! Included positive lags.
      logical, intent(in) :: include_mean !! Whether to prepend an intercept.
      real(dp), allocatable, intent(out) :: response(:) !! Aligned response.
      real(dp), allocatable, intent(out) :: design(:, :) !! Aligned regression matrix.
      integer, intent(out) :: info !! Zero on success.
      integer :: maximum_lag, observations, columns, column

      info = 1
      if (size(lags) < 1 .or. any(lags < 1) .or. &
         .not. all(ieee_is_finite(series))) return
      maximum_lag = maxval(lags)
      observations = size(series) - maximum_lag
      columns = size(lags) + merge(1, 0, include_mean)
      if (observations <= columns) return
      allocate(response(observations), design(observations, columns))
      response = series(maximum_lag + 1:)
      column = 0
      if (include_mean) then
         column = 1
         design(:, column) = 1.0_dp
      end if
      do maximum_lag = 1, size(lags)
         column = column + 1
         design(:, column) = series(maxval(lags) + 1 - lags(maximum_lag): &
            size(series) - lags(maximum_lag))
      end do
      info = 0
   end subroutine build_lagged_regression

   pure function make_tvar_model(response, design, log_variance, initial) &
      result(model)
      !! Construct the dynamic-regression state-space model used by NTS tvAR.
      real(dp), intent(in) :: response(:) !! Aligned response.
      real(dp), intent(in) :: design(:, :) !! Time-varying observation rows.
      real(dp), intent(in) :: log_variance(:) !! Observation then state log variances.
      real(dp), intent(in) :: initial(:) !! Initial coefficient mean.
      type(ssm_model_t) :: model
      integer :: observations, coefficients, time, coefficient

      observations = size(response)
      coefficients = size(design, 2)
      allocate(model%y(observations, 1), model%z(1, coefficients, observations))
      allocate(model%h(1, 1, 1), model%transition(coefficients, coefficients, 1))
      allocate(model%r(coefficients, coefficients, 1))
      allocate(model%q(coefficients, coefficients, 1), model%a1(coefficients))
      allocate(model%p1(coefficients, coefficients), &
         model%p1inf(coefficients, coefficients))
      allocate(model%missing(observations, 1))
      model%y(:, 1) = response
      do time = 1, observations
         model%z(1, :, time) = design(time, :)
      end do
      model%h(1, 1, 1) = exp(log_variance(1))
      model%transition = 0.0_dp
      model%r = 0.0_dp
      model%q = 0.0_dp
      model%p1 = 0.0_dp
      do coefficient = 1, coefficients
         model%transition(coefficient, coefficient, 1) = 1.0_dp
         model%r(coefficient, coefficient, 1) = 1.0_dp
         model%q(coefficient, coefficient, 1) = exp(log_variance(coefficient + 1))
         model%p1(coefficient, coefficient) = 1.0_dp
      end do
      model%a1 = initial
      model%p1inf = 0.0_dp
      model%missing = .false.
   end function make_tvar_model

   pure function diagonal_values(matrix) result(values)
      !! Return the main diagonal of a square matrix.
      real(dp), intent(in) :: matrix(:, :) !! Square input matrix.
      real(dp) :: values(min(size(matrix, 1), size(matrix, 2)))
      integer :: index

      do index = 1, size(values)
         values(index) = matrix(index, index)
      end do
   end function diagonal_values

   pure subroutine normalize_log_weights(log_weight, weight, log_normalizer, info)
      !! Normalize log weights and return their log sum.
      real(dp), intent(in) :: log_weight(:) !! Unnormalized log weights.
      real(dp), intent(out) :: weight(:) !! Normalized nonnegative weights.
      real(dp), intent(out) :: log_normalizer !! Log of the unnormalized weight sum.
      integer, intent(out) :: info !! Zero on success.
      real(dp) :: maximum_weight, total

      info = 1
      weight = 0.0_dp
      log_normalizer = -huge(1.0_dp)
      if (size(weight) /= size(log_weight) .or. size(weight) < 1) return
      maximum_weight = maxval(log_weight)
      if (.not. ieee_is_finite(maximum_weight)) return
      weight = exp(log_weight - maximum_weight)
      total = sum(weight)
      if (total <= 0.0_dp .or. .not. ieee_is_finite(total)) return
      weight = weight/total
      log_normalizer = maximum_weight + log(total)
      info = 0
   end subroutine normalize_log_weights

   pure subroutine cfar_design_matrix(series, order, spline_df, design, response)
      !! Construct spline-convolution regressors for a regular functional series.
      real(dp), intent(in) :: series(:, :) !! Time-by-grid functional observations.
      integer, intent(in) :: order !! Functional autoregressive order.
      integer, intent(in) :: spline_df !! Spline degrees of freedom excluding intercept.
      real(dp), allocatable, intent(out) :: design(:, :) !! Convolution design matrix.
      real(dp), allocatable, intent(out) :: response(:) !! Vectorized aligned responses.
      real(dp), allocatable :: basis(:)
      real(dp) :: difference
      integer :: observations, points, basis_count, blocks
      integer :: block, time, location, source, lag, column, row

      observations = size(series, 1)
      points = size(series, 2)
      basis_count = spline_df + 1
      blocks = observations - order
      allocate(design(blocks*points, order*basis_count))
      allocate(response(blocks*points), basis(basis_count))
      design = 0.0_dp
      do block = 1, blocks
         time = order + block
         do location = 1, points
            row = (block - 1)*points + location
            response(row) = series(time, location)
            do lag = 1, order
               do source = 1, points
                  difference = real(location - source, dp)/real(points - 1, dp)
                  call cfar_spline_basis(difference, basis)
                  do column = 1, basis_count
                     design(row, (lag - 1)*basis_count + column) = &
                        design(row, (lag - 1)*basis_count + column) + &
                        basis(column)*series(time - lag, source)/real(points, dp)
                  end do
               end do
            end do
         end do
      end do
   end subroutine cfar_design_matrix

   pure subroutine cfar_gls(design, response, blocks, points, rho, coefficient, &
      fitted, residual, score, info, sigma, log_likelihood)
      !! Fit one CFAR design for a fixed OU decay parameter.
      real(dp), intent(in) :: design(:, :) !! Vectorized convolution design matrix.
      real(dp), intent(in) :: response(:) !! Vectorized functional responses.
      integer, intent(in) :: blocks !! Number of functional response curves.
      integer, intent(in) :: points !! Number of grid points per curve.
      real(dp), intent(in) :: rho !! Positive OU decay parameter.
      real(dp), allocatable, intent(out) :: coefficient(:) !! GLS coefficients.
      real(dp), allocatable, intent(out) :: fitted(:) !! Vectorized fitted values.
      real(dp), allocatable, intent(out) :: residual(:) !! Vectorized residuals.
      real(dp), intent(out) :: score !! Negative profiled log likelihood.
      integer, intent(out) :: info !! Zero on success.
      real(dp), intent(out), optional :: sigma !! Profile innovation standard deviation.
      real(dp), intent(out), optional :: log_likelihood !! Profile log likelihood.
      real(dp), allocatable :: cross_product(:, :), inverse(:, :), right_hand(:)
      real(dp) :: adjacent, diagonal, sse, actual_sse, log_determinant, likelihood
      integer :: columns, block, location, row, next_row, status

      info = 1
      score = huge(1.0_dp)
      columns = size(design, 2)
      if (rho <= 0.0_dp .or. size(response) /= blocks*points .or. &
         size(design, 1) /= size(response) .or. blocks < 1 .or. points < 2) return
      adjacent = exp(-rho/real(points - 1, dp))
      if (adjacent >= 1.0_dp) return
      allocate(cross_product(columns, columns), right_hand(columns))
      cross_product = 0.0_dp
      right_hand = 0.0_dp
      do block = 1, blocks
         do location = 1, points
            row = (block - 1)*points + location
            diagonal = 1.0_dp
            if (location > 1 .and. location < points) diagonal = &
               1.0_dp + adjacent*adjacent
            cross_product = cross_product + diagonal*outer_product( &
               design(row, :), design(row, :))
            right_hand = right_hand + diagonal*design(row, :)*response(row)
            if (location < points) then
               next_row = row + 1
               cross_product = cross_product - adjacent*(outer_product( &
                  design(row, :), design(next_row, :)) + outer_product( &
                  design(next_row, :), design(row, :)))
               right_hand = right_hand - adjacent*(design(row, :)* &
                  response(next_row) + design(next_row, :)*response(row))
            end if
         end do
      end do
      call invert_matrix(cross_product, inverse, status)
      if (status /= 0) then
         info = 2
         return
      end if
      coefficient = matmul(inverse, right_hand)
      fitted = matmul(design, coefficient)
      residual = response - fitted
      sse = cfar_weighted_sse_vector(residual, blocks, points, adjacent)
      if (sse <= tiny(1.0_dp) .or. .not. ieee_is_finite(sse)) then
         info = 3
         return
      end if
      actual_sse = sse/(1.0_dp - adjacent*adjacent)
      log_determinant = -real(points - 1, dp)*log(1.0_dp - adjacent*adjacent)
      likelihood = 0.5_dp*real(blocks, dp)*log_determinant - &
         0.5_dp*real(blocks*points, dp)*(1.0_dp + log(2.0_dp*acos(-1.0_dp)* &
         actual_sse/real(blocks*points, dp)))
      score = -likelihood
      if (present(sigma)) sigma = sqrt(actual_sse/real(blocks*points, dp))
      if (present(log_likelihood)) log_likelihood = likelihood
      info = 0
   end subroutine cfar_gls

   pure subroutine cfar_spline_basis(value, basis)
      !! Evaluate the cubic regression-spline basis used for CFAR kernels.
      real(dp), intent(in) :: value !! Kernel argument on minus one to one.
      real(dp), intent(out) :: basis(:) !! Evaluated basis functions.
      real(dp) :: knot
      integer :: column, knot_count

      basis = 0.0_dp
      if (size(basis) >= 1) basis(1) = 1.0_dp
      if (size(basis) >= 2) basis(2) = value
      if (size(basis) >= 3) basis(3) = value*value
      if (size(basis) >= 4) basis(4) = value*value*value
      knot_count = size(basis) - 4
      do column = 1, knot_count
         knot = -1.0_dp + 2.0_dp*real(column, dp)/real(knot_count + 1, dp)
         basis(column + 4) = max(0.0_dp, value - knot)**3
      end do
   end subroutine cfar_spline_basis

   pure subroutine cfar_sample_kernels(coefficient, points, kernel)
      !! Evaluate fitted spline kernels on the full difference grid.
      real(dp), intent(in) :: coefficient(:, :) !! Lag-by-basis coefficients.
      integer, intent(in) :: points !! Number of functional grid points.
      real(dp), allocatable, intent(out) :: kernel(:, :) !! Lag-by-difference-grid kernels.
      real(dp), allocatable :: basis(:)
      real(dp) :: value
      integer :: lag, index

      allocate(kernel(2*points - 1, size(coefficient, 1)))
      allocate(basis(size(coefficient, 2)))
      do lag = 1, size(coefficient, 1)
         do index = 1, 2*points - 1
            value = real(index - points, dp)/real(points - 1, dp)
            call cfar_spline_basis(value, basis)
            kernel(index, lag) = dot_product(coefficient(lag, :), basis)
         end do
      end do
   end subroutine cfar_sample_kernels

   pure real(dp) function cfar_weighted_sse(residual, rho) result(sse)
      !! Return OU precision-weighted residual sum of squares for complete curves.
      real(dp), intent(in) :: residual(:, :) !! Curve-by-grid residual matrix.
      real(dp), intent(in) :: rho !! Positive OU decay parameter.
      real(dp) :: adjacent

      adjacent = exp(-rho/real(size(residual, 2) - 1, dp))
      sse = cfar_weighted_sse_vector(reshape(transpose(residual), &
         [size(residual)]), size(residual, 1), size(residual, 2), adjacent)
   end function cfar_weighted_sse

   pure real(dp) function cfar_weighted_sse_vector(residual, blocks, points, &
      adjacent) result(sse)
      !! Return a tri-diagonal OU precision quadratic form.
      real(dp), intent(in) :: residual(:) !! Vectorized curve residuals.
      integer, intent(in) :: blocks !! Number of curves.
      integer, intent(in) :: points !! Grid points in each curve.
      real(dp), intent(in) :: adjacent !! Adjacent-grid OU correlation.
      real(dp) :: diagonal
      integer :: block, location, row

      sse = 0.0_dp
      do block = 1, blocks
         do location = 1, points
            row = (block - 1)*points + location
            diagonal = 1.0_dp
            if (location > 1 .and. location < points) diagonal = &
               1.0_dp + adjacent*adjacent
            sse = sse + diagonal*residual(row)*residual(row)
            if (location < points) sse = sse - &
               2.0_dp*adjacent*residual(row)*residual(row + 1)
         end do
      end do
   end function cfar_weighted_sse_vector

   pure function outer_product(left, right) result(product)
      !! Return the outer product of two real vectors.
      real(dp), intent(in) :: left(:) !! Column vector values.
      real(dp), intent(in) :: right(:) !! Row vector values.
      real(dp) :: product(size(left), size(right))
      integer :: column

      do column = 1, size(right)
         product(:, column) = left*right(column)
      end do
   end function outer_product

   pure subroutine cfar_resample_curves(series, points, resampled)
      !! Linearly resample complete curves onto a requested regular grid.
      real(dp), intent(in) :: series(:, :) !! Time-by-source-grid curves.
      integer, intent(in) :: points !! Number of target grid points.
      real(dp), allocatable, intent(out) :: resampled(:, :) !! Resampled curves.
      real(dp) :: coordinate, source_coordinate, fraction
      integer :: curve, location, left, source_points

      source_points = size(series, 2)
      allocate(resampled(size(series, 1), points))
      do curve = 1, size(series, 1)
         do location = 1, points
            coordinate = real(location - 1, dp)/real(points - 1, dp)
            source_coordinate = coordinate*real(source_points - 1, dp)
            left = min(source_points - 1, int(source_coordinate) + 1)
            fraction = source_coordinate - real(left - 1, dp)
            resampled(curve, location) = (1.0_dp - fraction)*series(curve, left) + &
               fraction*series(curve, left + 1)
         end do
      end do
   end subroutine cfar_resample_curves

   pure subroutine cfar_register_irregular(values, positions, counts, points, &
      registered, info)
      !! Linearly register irregular curves on a common unit-interval grid.
      real(dp), intent(in) :: values(:, :) !! Time-by-capacity irregular values.
      real(dp), intent(in) :: positions(:, :) !! Matching increasing positions.
      integer, intent(in) :: counts(:) !! Valid entries per curve.
      integer, intent(in) :: points !! Number of target grid points.
      real(dp), allocatable, intent(out) :: registered(:, :) !! Registered curves.
      integer, intent(out) :: info !! Zero on success.
      real(dp) :: coordinate, fraction, denominator
      integer :: curve, location, right, count

      info = 1
      if (points < 3 .or. size(values, 1) /= size(positions, 1) .or. &
         size(values, 2) /= size(positions, 2) .or. &
         size(counts) /= size(values, 1) .or. any(counts < 2) .or. &
         any(counts > size(values, 2))) return
      allocate(registered(size(values, 1), points))
      do curve = 1, size(values, 1)
         count = counts(curve)
         if (positions(curve, 1) < 0.0_dp .or. &
            positions(curve, count) > 1.0_dp .or. &
            .not. strictly_increasing(positions(curve, :count)) .or. &
            .not. all(ieee_is_finite(values(curve, :count))) .or. &
            .not. all(ieee_is_finite(positions(curve, :count)))) return
         do location = 1, points
            coordinate = real(location - 1, dp)/real(points - 1, dp)
            if (coordinate <= positions(curve, 1)) then
               registered(curve, location) = values(curve, 1)
            else if (coordinate >= positions(curve, count)) then
               registered(curve, location) = values(curve, count)
            else
               right = 2
               do while (positions(curve, right) < coordinate)
                  right = right + 1
               end do
               denominator = positions(curve, right) - positions(curve, right - 1)
               fraction = (coordinate - positions(curve, right - 1))/denominator
               registered(curve, location) = (1.0_dp - fraction)* &
                  values(curve, right - 1) + fraction*values(curve, right)
            end if
         end do
      end do
      info = 0
   end subroutine cfar_register_irregular

   integer function random_poisson(mean_value) result(draw)
      !! Generate one Poisson variate by product inversion.
      real(dp), intent(in) :: mean_value !! Nonnegative Poisson mean.
      real(dp) :: product, limit

      if (mean_value <= 0.0_dp) then
         draw = 0
         return
      end if
      limit = exp(-mean_value)
      product = 1.0_dp
      draw = -1
      do
         draw = draw + 1
         product = product*max(random_uniform(), tiny(1.0_dp))
         if (product <= limit) exit
      end do
   end function random_poisson

   pure subroutine sort_small(values)
      !! Sort a short real vector in ascending order in place.
      real(dp), intent(inout) :: values(:) !! Values to sort.
      real(dp) :: held
      integer :: position, previous

      do position = 2, size(values)
         held = values(position)
         previous = position - 1
         do while (previous >= 1)
            if (values(previous) <= held) exit
            values(previous + 1) = values(previous)
            previous = previous - 1
         end do
         values(previous + 1) = held
      end do
   end subroutine sort_small

   pure real(dp) function interpolate_regular_curve(curve, position) result(value)
      !! Linearly interpolate a regular unit-interval curve at one position.
      real(dp), intent(in) :: curve(:) !! Curve values on a regular grid.
      real(dp), intent(in) :: position !! Position on the unit interval.
      real(dp) :: coordinate, fraction
      integer :: left

      coordinate = max(0.0_dp, min(1.0_dp, position))*real(size(curve) - 1, dp)
      left = min(size(curve) - 1, int(coordinate) + 1)
      fraction = coordinate - real(left - 1, dp)
      value = (1.0_dp - fraction)*curve(left) + fraction*curve(left + 1)
   end function interpolate_regular_curve

   pure subroutine average_ranks(values, ranks)
      !! Assign average ranks with ties.
      real(dp), intent(in) :: values(:) !! Values to rank.
      real(dp), allocatable, intent(out) :: ranks(:) !! One-based average ranks.
      integer, allocatable :: order(:)
      integer :: first, last, position
      real(dp) :: average

      call ascending_indices(values, order)
      allocate(ranks(size(values)))
      first = 1
      do while (first <= size(values))
         last = first
         do while (last < size(values))
            if (values(order(last + 1)) /= values(order(first))) exit
            last = last + 1
         end do
         average = 0.5_dp*real(first + last, dp)
         do position = first, last
            ranks(order(position)) = average
         end do
         first = last + 1
      end do
   end subroutine average_ranks

   pure integer function categorical_state(probabilities, uniform_draw) result(state)
      !! Select a categorical state from one uniform draw.
      real(dp), intent(in) :: probabilities(:) !! Categorical probabilities.
      real(dp), intent(in) :: uniform_draw !! Uniform draw on the unit interval.
      real(dp) :: cumulative
      integer :: candidate

      state = size(probabilities)
      cumulative = 0.0_dp
      do candidate = 1, size(probabilities)
         cumulative = cumulative + probabilities(candidate)
         if (uniform_draw <= cumulative) then
            state = candidate
            return
         end if
      end do
   end function categorical_state

   pure subroutine ascending_indices(values, index)
      !! Return indices that arrange values in ascending order.
      real(dp), intent(in) :: values(:) !! Values to arrange.
      integer, allocatable, intent(out) :: index(:) !! Ascending-order indices.
      integer :: position, previous, held

      allocate(index(size(values)))
      index = [(position, position=1, size(values))]
      do position = 2, size(values)
         held = index(position)
         previous = position - 1
         do while (previous >= 1)
            if (values(index(previous)) <= values(held)) exit
            index(previous + 1) = index(previous)
            previous = previous - 1
         end do
         index(previous + 1) = held
      end do
   end subroutine ascending_indices

   pure real(dp) function f_upper_probability(statistic, numerator_df, &
      denominator_df) result(probability)
      !! Return the upper-tail probability of an F statistic.
      real(dp), intent(in) :: statistic !! Nonnegative F statistic.
      integer, intent(in) :: numerator_df !! Numerator degrees of freedom.
      integer, intent(in) :: denominator_df !! Denominator degrees of freedom.
      real(dp) :: argument

      if (statistic <= 0.0_dp) then
         probability = 1.0_dp
         return
      end if
      argument = real(denominator_df, dp)/(real(denominator_df, dp) + &
         real(numerator_df, dp)*statistic)
      probability = regularized_beta(argument, 0.5_dp*real(denominator_df, dp), &
         0.5_dp*real(numerator_df, dp))
      probability = max(0.0_dp, min(1.0_dp, probability))
   end function f_upper_probability

   pure logical function valid_tar_specification(intercept, ar, ar_order, &
      thresholds, innovation_sd, delay) result(valid)
      !! Check dimensions and values of a SETAR simulation specification.
      real(dp), intent(in) :: intercept(:) !! Regime intercepts.
      real(dp), intent(in) :: ar(:, :) !! Regime-by-lag AR coefficients.
      integer, intent(in) :: ar_order(:) !! Regime-specific AR orders.
      real(dp), intent(in) :: thresholds(:) !! Ordered regime thresholds.
      real(dp), intent(in) :: innovation_sd(:) !! Regime innovation standard deviations.
      integer, intent(in) :: delay !! Delay of the threshold variable.
      integer :: regimes

      regimes = size(ar_order)
      valid = regimes >= 2 .and. size(intercept) == regimes .and. &
         size(ar, 1) == regimes .and. size(thresholds) == regimes - 1 .and. &
         size(innovation_sd) == regimes .and. delay >= 1
      if (.not. valid) return
      valid = all(ar_order >= 0) .and. all(ar_order <= size(ar, 2)) .and. &
         all(innovation_sd > 0.0_dp) .and. strictly_increasing(thresholds) .and. &
         all(ieee_is_finite(intercept)) .and. all(ieee_is_finite(ar)) .and. &
         all(ieee_is_finite(thresholds)) .and. &
         all(ieee_is_finite(innovation_sd))
   end function valid_tar_specification

   pure logical function strictly_increasing(values) result(increasing)
      !! Return whether a vector is strictly increasing.
      real(dp), intent(in) :: values(:) !! Values to inspect.
      integer :: index

      increasing = .true.
      do index = 2, size(values)
         if (values(index) <= values(index - 1)) then
            increasing = .false.
            return
         end if
      end do
   end function strictly_increasing

   pure integer function threshold_regime(value, thresholds) result(regime)
      !! Map one threshold-variable value to its regime index.
      real(dp), intent(in) :: value !! Threshold-variable value.
      real(dp), intent(in) :: thresholds(:) !! Ordered regime thresholds.
      integer :: index

      regime = 1
      do index = 1, size(thresholds)
         if (value > thresholds(index)) regime = index + 1
      end do
   end function threshold_regime

   pure integer function count_regime(values, thresholds, selected) result(number)
      !! Count observations assigned to one threshold regime.
      real(dp), intent(in) :: values(:) !! Threshold-variable values.
      real(dp), intent(in) :: thresholds(:) !! Ordered regime thresholds.
      integer, intent(in) :: selected !! Regime index to count.
      integer :: row

      number = 0
      do row = 1, size(values)
         if (threshold_regime(values(row), thresholds) == selected) &
            number = number + 1
      end do
   end function count_regime

   pure subroutine collect_regime_indices(values, thresholds, selected, indices)
      !! Collect observation indices assigned to one threshold regime.
      real(dp), intent(in) :: values(:) !! Threshold-variable values.
      real(dp), intent(in) :: thresholds(:) !! Ordered regime thresholds.
      integer, intent(in) :: selected !! Regime index to collect.
      integer, intent(out) :: indices(:) !! Collected observation indices.
      integer :: row, position

      position = 0
      do row = 1, size(values)
         if (threshold_regime(values(row), thresholds) == selected) then
            position = position + 1
            indices(position) = row
         end if
      end do
   end subroutine collect_regime_indices

end module nts_mod
