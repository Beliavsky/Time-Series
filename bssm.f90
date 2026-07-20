! SPDX-License-Identifier: GPL-2.0-or-later
! SPDX-FileComment: Algorithms translated from the R bssm package.
! Particle-inference algorithms translated from the GPL bssm package.
module bssm_mod
   use kind_mod, only: dp
   use special_functions_mod, only: normal_log_density, &
      multivariate_normal_log_density
   use linalg_mod, only: symmetric_eigen, outer_product
   use linalg_mod, only: inverse_logdet, identity_matrix
   use random_mod, only: random_standard_normal_matrix
   use kfas_mod, only: ssm_model_t, kfs_filter_t, kfs_smoother_t
   use kfas_mod, only: kfs_filter, kfs_smooth, kfs_fast_smooth, validate_ssm
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value
   use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan
   implicit none
   private

   integer, parameter, public :: bssm_svm = 0
   integer, parameter, public :: bssm_poisson = 1
   integer, parameter, public :: bssm_binomial = 2
   integer, parameter, public :: bssm_negative_binomial = 3
   integer, parameter, public :: bssm_gamma = 4
   integer, parameter, public :: bssm_gaussian = 5
   integer, parameter, public :: bssm_is1 = 1
   integer, parameter, public :: bssm_is2 = 2
   integer, parameter, public :: bssm_is3 = 3
   integer, parameter, public :: bssm_sokal = 1
   integer, parameter, public :: bssm_geyer = 2

   type, public :: bssm_particle_filter_t
      ! Bootstrap-filter particles, weights, ancestry, and state summaries.
      real(dp), allocatable :: particles(:, :, :), weights(:, :)
      integer, allocatable :: ancestors(:, :)
      real(dp), allocatable :: predicted_mean(:, :), filtered_mean(:, :)
      real(dp), allocatable :: predicted_covariance(:, :, :)
      real(dp), allocatable :: filtered_covariance(:, :, :)
      real(dp) :: log_likelihood = 0.0_dp
      integer :: info = 0
   end type bssm_particle_filter_t

   type, public :: bssm_particle_smoother_t
      ! Genealogical particle paths and weighted smoothing summaries.
      real(dp), allocatable :: trajectories(:, :, :), mean(:, :)
      real(dp), allocatable :: covariance(:, :, :)
      integer :: info = 0
   end type bssm_particle_smoother_t

   type, public :: bssm_simulation_smoother_t
      ! Conditional Gaussian state paths and completed observation paths.
      real(dp), allocatable :: trajectories(:, :, :)
      real(dp), allocatable :: observations(:, :, :)
      real(dp), allocatable :: mean(:, :), covariance(:, :, :)
      integer :: info = 0
   end type bssm_simulation_smoother_t

   type, public :: bssm_prediction_t
      ! Predictive state, signal, inverse-link mean, and response samples.
      real(dp), allocatable :: state(:, :, :), signal(:, :, :)
      real(dp), allocatable :: mean(:, :, :), response(:, :, :)
      integer :: info = 0
   end type bssm_prediction_t

   type, public :: bssm_prediction_summary_t
      ! Weighted fitted signal and observation summaries over posterior draws.
      real(dp), allocatable :: probabilities(:)
      real(dp), allocatable :: signal_mean(:, :), signal_sd(:, :)
      real(dp), allocatable :: signal_quantile(:, :, :)
      real(dp), allocatable :: observation_mean(:, :), observation_sd(:, :)
      real(dp), allocatable :: observation_quantile(:, :, :)
      integer :: info = 0
   end type bssm_prediction_summary_t

   type, public :: bssm_ekf_t
      ! Extended Kalman filter states, covariances, and likelihood diagnostics.
      real(dp), allocatable :: predicted_mean(:, :), filtered_mean(:, :)
      real(dp), allocatable :: predicted_covariance(:, :, :)
      real(dp), allocatable :: filtered_covariance(:, :, :)
      real(dp), allocatable :: innovation(:), innovation_variance(:)
      integer, allocatable :: iterations(:)
      real(dp) :: log_likelihood = 0.0_dp
      integer :: info = 0
   end type bssm_ekf_t

   type, public :: bssm_multivariate_ekf_t
      ! Multivariate EKF states, innovations, covariances, and likelihood.
      real(dp), allocatable :: predicted_mean(:, :), filtered_mean(:, :)
      real(dp), allocatable :: predicted_covariance(:, :, :)
      real(dp), allocatable :: filtered_covariance(:, :, :)
      real(dp), allocatable :: innovation(:, :)
      real(dp), allocatable :: innovation_covariance(:, :, :)
      integer, allocatable :: iterations(:)
      real(dp) :: log_likelihood = 0.0_dp
      integer :: info = 0
   end type bssm_multivariate_ekf_t

   type, public :: bssm_ekf_smoother_t
      ! Extended Kalman smoothed states, covariances, and likelihood.
      real(dp), allocatable :: state(:, :), covariance(:, :, :)
      integer, allocatable :: iterations(:)
      real(dp) :: log_likelihood = 0.0_dp
      integer :: info = 0
   end type bssm_ekf_smoother_t

   type, public :: bssm_gaussian_approximation_t
      ! Laplace Gaussian approximation and its conditional mode diagnostics.
      real(dp), allocatable :: pseudo_observation(:)
      real(dp), allocatable :: observation_variance(:)
      real(dp), allocatable :: mode_state(:, :), mode_signal(:)
      real(dp), allocatable :: scaling(:)
      real(dp), allocatable :: proposal_factor(:, :, :)
      real(dp), allocatable :: conditional_matrix(:, :, :)
      real(dp) :: gaussian_log_likelihood = 0.0_dp
      real(dp) :: corrected_log_likelihood = 0.0_dp
      real(dp) :: difference = 0.0_dp
      integer :: iterations = 0
      integer :: info = 0
      logical :: converged = .false.
   end type bssm_gaussian_approximation_t

   type, public :: bssm_importance_sample_t
      ! Gaussian-proposal state trajectories and importance diagnostics.
      real(dp), allocatable :: trajectories(:, :, :)
      real(dp), allocatable :: log_weight(:), weight(:)
      real(dp), allocatable :: mean(:, :), covariance(:, :, :)
      real(dp) :: log_likelihood = 0.0_dp
      real(dp) :: effective_sample_size = 0.0_dp
      integer :: info = 0
   end type bssm_importance_sample_t

   type, public :: bssm_nonlinear_approximation_t
      ! Global iterated-EKS approximation and conditional proposal moments.
      real(dp), allocatable :: mode_state(:, :), proposal_mean(:, :)
      real(dp), allocatable :: observation_intercept(:)
      real(dp), allocatable :: observation_jacobian(:, :)
      real(dp), allocatable :: observation_standard_deviation(:)
      real(dp), allocatable :: transition_intercept(:, :)
      real(dp), allocatable :: transition_jacobian(:, :, :)
      real(dp), allocatable :: noise_loading(:, :, :)
      real(dp), allocatable :: proposal_factor(:, :, :)
      real(dp), allocatable :: conditional_matrix(:, :, :)
      real(dp), allocatable :: scaling(:)
      real(dp) :: gaussian_log_likelihood = 0.0_dp
      real(dp) :: corrected_log_likelihood = 0.0_dp
      real(dp) :: difference = 0.0_dp
      integer :: iterations = 0
      integer :: info = 0
      logical :: converged = .false.
   end type bssm_nonlinear_approximation_t

   type, public :: bssm_multivariate_nonlinear_approximation_t
      ! Global Gaussian approximation to a vector nonlinear state-space model.
      real(dp), allocatable :: mode_state(:, :), proposal_mean(:, :)
      real(dp), allocatable :: observation_intercept(:, :)
      real(dp), allocatable :: observation_jacobian(:, :, :)
      real(dp), allocatable :: observation_noise_loading(:, :, :)
      real(dp), allocatable :: transition_intercept(:, :)
      real(dp), allocatable :: transition_jacobian(:, :, :)
      real(dp), allocatable :: state_noise_loading(:, :, :)
      real(dp), allocatable :: proposal_factor(:, :, :)
      real(dp), allocatable :: conditional_matrix(:, :, :)
      real(dp), allocatable :: scaling(:)
      real(dp) :: gaussian_log_likelihood = 0.0_dp
      real(dp) :: corrected_log_likelihood = 0.0_dp
      real(dp) :: difference = 0.0_dp
      integer :: iterations = 0
      integer :: info = 0
      logical :: converged = .false.
   end type bssm_multivariate_nonlinear_approximation_t

   type, public :: bssm_multivariate_approximation_t
      ! Mixed-family multivariate Laplace approximation and proposal moments.
      real(dp), allocatable :: pseudo_observation(:, :)
      real(dp), allocatable :: observation_variance(:, :)
      real(dp), allocatable :: mode_state(:, :), mode_signal(:, :)
      real(dp), allocatable :: scaling(:, :)
      real(dp), allocatable :: proposal_factor(:, :, :)
      real(dp), allocatable :: conditional_matrix(:, :, :)
      real(dp) :: gaussian_log_likelihood = 0.0_dp
      real(dp) :: corrected_log_likelihood = 0.0_dp
      real(dp) :: difference = 0.0_dp
      integer :: iterations = 0
      integer :: info = 0
      logical :: converged = .false.
   end type bssm_multivariate_approximation_t

   type, public :: bssm_mcmc_t
      ! Retained parameter chain and Metropolis-Hastings diagnostics.
      real(dp), allocatable :: parameters(:, :)
      real(dp), allocatable :: log_likelihood(:), log_prior(:)
      real(dp), allocatable :: acceptance_probability(:)
      real(dp), allocatable :: final_proposal_factor(:, :)
      logical, allocatable :: accepted(:)
      real(dp) :: acceptance_rate = 0.0_dp
      integer :: info = 0
   end type bssm_mcmc_t

   type, public :: bssm_likelihood_estimate_t
      ! Log-likelihood estimate and estimator status.
      real(dp) :: value = -huge(1.0_dp)
      integer :: info = 0
   end type bssm_likelihood_estimate_t

   type, public :: bssm_da_mcmc_t
      ! Delayed-acceptance parameter chain and two-stage diagnostics.
      real(dp), allocatable :: parameters(:, :)
      real(dp), allocatable :: coarse_log_likelihood(:)
      real(dp), allocatable :: fine_log_likelihood(:), log_prior(:)
      real(dp), allocatable :: first_stage_probability(:)
      real(dp), allocatable :: second_stage_probability(:)
      real(dp), allocatable :: final_proposal_factor(:, :)
      logical, allocatable :: first_stage_accepted(:), accepted(:)
      real(dp) :: first_stage_acceptance_rate = 0.0_dp
      real(dp) :: acceptance_rate = 0.0_dp
      integer :: info = 0
   end type bssm_da_mcmc_t

   type, public :: bssm_post_correction_t
      ! Importance weights and corrected parameter summaries.
      real(dp), allocatable :: log_weight(:), weight(:)
      real(dp), allocatable :: parameter_mean(:), parameter_covariance(:, :)
      real(dp) :: log_mean_weight = -huge(1.0_dp)
      real(dp) :: effective_sample_size = 0.0_dp
      integer :: info = 0
   end type bssm_post_correction_t

   type, public :: bssm_sde_state_sample_t
      ! IS2-corrected SDE trajectories and importance diagnostics.
      real(dp), allocatable :: trajectories(:, :, :)
      real(dp), allocatable :: fine_log_likelihood(:)
      real(dp), allocatable :: log_weight(:), weight(:)
      real(dp) :: log_mean_weight = -huge(1.0_dp)
      real(dp) :: effective_sample_size = 0.0_dp
      integer :: info = 0
   end type bssm_sde_state_sample_t

   type, public :: bssm_state_posterior_t
      ! Corrected state or signal means and total covariances.
      real(dp), allocatable :: mean(:, :), covariance(:, :, :)
      integer :: info = 0
   end type bssm_state_posterior_t

   type, public :: bssm_trajectory_sample_t
      ! Resampled corrected trajectories and their source indices.
      real(dp), allocatable :: trajectories(:, :, :)
      integer, allocatable :: source_sample(:), source_particle(:)
      integer :: info = 0
   end type bssm_trajectory_sample_t

   type, public :: bssm_particle_count_t
      ! Log-likelihood variability and selected particle count.
      real(dp), allocatable :: standard_deviation(:)
      integer :: particle_count = 0
      integer :: selected_index = 0
      integer :: info = 0
   end type bssm_particle_count_t

   type, public :: bssm_mcmc_diagnostics_t
      ! Weighted estimates of autocorrelation, Monte Carlo error, and ESS.
      real(dp), allocatable :: mean(:), variance(:), iact(:)
      real(dp), allocatable :: asymptotic_variance(:), mcse(:), ess(:)
      integer :: info = 0
   end type bssm_mcmc_diagnostics_t

   abstract interface
      pure function bssm_nonlinear_log_density_t(time, observation, state, &
         parameters) result(log_density)
         !! Evaluate a nonlinear observation log density.
         import dp
         integer, intent(in) :: time !! Observation times.
         real(dp), intent(in) :: observation !! Observed value or vector.
         real(dp), intent(in) :: state(:) !! State vector or state sequence.
         real(dp), intent(in) :: parameters(:) !! Model parameter values.
         real(dp) :: log_density
      end function bssm_nonlinear_log_density_t

      pure subroutine bssm_nonlinear_transition_t(time, state, parameters, &
         mean, noise_loading)
         !! Evaluate a nonlinear state transition and noise loading.
         import dp
         integer, intent(in) :: time !! Observation times.
         real(dp), intent(in) :: state(:) !! State vector or state sequence.
         real(dp), intent(in) :: parameters(:) !! Model parameter values.
         real(dp), intent(out) :: mean(:) !! Mean value or vector.
         real(dp), intent(out) :: noise_loading(:, :) !! Noise loading.
      end subroutine bssm_nonlinear_transition_t

      pure subroutine bssm_gaussian_observation_t(time, state, parameters, &
         mean, jacobian, standard_deviation)
         !! Evaluate a scalar Gaussian observation model and Jacobian.
         import dp
         integer, intent(in) :: time !! Observation times.
         real(dp), intent(in) :: state(:) !! State vector or state sequence.
         real(dp), intent(in) :: parameters(:) !! Model parameter values.
         real(dp), intent(out) :: mean !! Mean value or vector.
         real(dp), intent(out) :: jacobian(:) !! Jacobian.
         real(dp), intent(out) :: standard_deviation !! Standard deviation.
      end subroutine bssm_gaussian_observation_t

      pure subroutine bssm_multivariate_gaussian_observation_t(time, state, &
         parameters, mean, jacobian, noise_loading)
         !! Evaluate a multivariate Gaussian observation model.
         import dp
         integer, intent(in) :: time !! Observation times.
         real(dp), intent(in) :: state(:) !! State vector or state sequence.
         real(dp), intent(in) :: parameters(:) !! Model parameter values.
         real(dp), intent(out) :: mean(:) !! Mean value or vector.
         real(dp), intent(out) :: jacobian(:, :) !! Jacobian.
         real(dp), intent(out) :: noise_loading(:, :) !! Noise loading.
      end subroutine bssm_multivariate_gaussian_observation_t

      pure subroutine bssm_nonlinear_prediction_observation_t(time, state, &
         parameters, mean, noise_loading)
         !! Evaluate a nonlinear prediction observation model.
         import dp
         integer, intent(in) :: time !! Observation times.
         real(dp), intent(in) :: state(:) !! State vector or state sequence.
         real(dp), intent(in) :: parameters(:) !! Model parameter values.
         real(dp), intent(out) :: mean(:) !! Mean value or vector.
         real(dp), intent(out) :: noise_loading(:, :) !! Noise loading.
      end subroutine bssm_nonlinear_prediction_observation_t

      pure subroutine bssm_posterior_observation_model_t(parameters, &
         observation_loading, phi, offset, auxiliary, noise_loading, &
         correlated_gaussian, info)
         !! Construct an observation model from posterior parameters.
         import dp
         real(dp), intent(in) :: parameters(:) !! Model parameter values.
         real(dp), intent(out) :: observation_loading(:, :, :) !! Observation loading matrix.
         real(dp), intent(out) :: phi(:) !! Autoregressive or model coefficient.
         real(dp), intent(out) :: offset(:, :) !! Known additive offset.
         real(dp), intent(out) :: auxiliary(:, :) !! Auxiliary.
         real(dp), intent(out) :: noise_loading(:, :, :) !! Noise loading.
         logical, intent(out) :: correlated_gaussian !! Flag controlling correlated gaussian.
         integer, intent(out) :: info !! Status code; zero indicates success.
      end subroutine bssm_posterior_observation_model_t

      pure subroutine bssm_nonlinear_transition_jacobian_t(time, state, &
         parameters, mean, jacobian, noise_loading)
         !! Evaluate a nonlinear transition and its Jacobian.
         import dp
         integer, intent(in) :: time !! Observation times.
         real(dp), intent(in) :: state(:) !! State vector or state sequence.
         real(dp), intent(in) :: parameters(:) !! Model parameter values.
         real(dp), intent(out) :: mean(:) !! Mean value or vector.
         real(dp), intent(out) :: jacobian(:, :) !! Jacobian.
         real(dp), intent(out) :: noise_loading(:, :) !! Noise loading.
      end subroutine bssm_nonlinear_transition_jacobian_t

      pure function bssm_sde_coefficient_t(state, parameters) result(value)
         !! Evaluate one state-dependent SDE coefficient.
         import dp
         real(dp), intent(in) :: state !! State vector or state sequence.
         real(dp), intent(in) :: parameters(:) !! Model parameter values.
         real(dp) :: value
      end function bssm_sde_coefficient_t

      pure function bssm_parameter_log_density_t(parameters) result(value)
         !! Evaluate a parameter log density.
         import dp
         real(dp), intent(in) :: parameters(:) !! Model parameter values.
         real(dp) :: value
      end function bssm_parameter_log_density_t

      pure function bssm_likelihood_estimator_t(parameters, draw_index) &
         result(estimate)
         !! Estimate a likelihood for one parameter vector and draw.
         import dp, bssm_likelihood_estimate_t
         real(dp), intent(in) :: parameters(:) !! Model parameter values.
         integer, intent(in) :: draw_index !! Index of draw.
         type(bssm_likelihood_estimate_t) :: estimate
      end function bssm_likelihood_estimator_t
   end interface

   interface bssm_approximation_simulation_draws
      module procedure bssm_scalar_approximation_simulation_draws
      module procedure bssm_multivariate_approximation_simulation_draws
   end interface bssm_approximation_simulation_draws

   interface bssm_approximation_simulation
      module procedure bssm_scalar_approximation_simulation
      module procedure bssm_multivariate_approximation_simulation
   end interface bssm_approximation_simulation

   public :: bssm_observation_log_density, bssm_stratified_resample
   public :: bssm_bootstrap_filter_draws, bssm_bootstrap_filter
   public :: bssm_bootstrap_filter_tv_draws, bssm_bootstrap_filter_tv
   public :: bssm_nonlinear_bootstrap_filter_draws
   public :: bssm_nonlinear_bootstrap_filter
   public :: bssm_multivariate_nonlinear_bootstrap_filter_draws
   public :: bssm_multivariate_nonlinear_bootstrap_filter
   public :: bssm_ekpf_draws, bssm_ekpf
   public :: bssm_multivariate_ekpf_draws, bssm_multivariate_ekpf
   public :: bssm_iekf, bssm_ekf_smoother, bssm_ekf_fast_smoother, bssm_ukf
   public :: bssm_multivariate_ukf
   public :: bssm_multivariate_iekf, bssm_multivariate_ekf_smoother
   public :: bssm_multivariate_ekf_fast_smoother
   public :: bssm_multivariate_gaussian_observation_t
   public :: bssm_laplace_pseudo_observation
   public :: bssm_gaussian_approximation
   public :: bssm_psi_filter_draws, bssm_psi_filter
   public :: bssm_spdk_importance_draws, bssm_spdk_importance
   public :: bssm_nonlinear_gaussian_approximation
   public :: bssm_multivariate_nonlinear_gaussian_approximation
   public :: bssm_nonlinear_psi_filter_draws, bssm_nonlinear_psi_filter
   public :: bssm_multivariate_nonlinear_psi_filter_draws
   public :: bssm_multivariate_nonlinear_psi_filter
   public :: bssm_multivariate_gaussian_approximation
   public :: bssm_multivariate_bootstrap_filter_draws
   public :: bssm_multivariate_bootstrap_filter
   public :: bssm_multivariate_psi_filter_draws
   public :: bssm_multivariate_psi_filter
   public :: bssm_multivariate_spdk_importance_draws
   public :: bssm_multivariate_spdk_importance
   public :: bssm_sde_euler_step, bssm_sde_milstein_step
   public :: bssm_sde_bootstrap_filter_draws, bssm_sde_bootstrap_filter
   public :: bssm_sde_is2_state_sampler_draws
   public :: bssm_sde_is2_state_sampler
   public :: bssm_sde_pmmh_draws, bssm_sde_pmmh
   public :: bssm_sde_da_pmmh_draws, bssm_sde_da_pmmh
   public :: bssm_pmmh_kernel_draws
   public :: bssm_da_pmmh_kernel_draws
   public :: bssm_approximate_mcmc_draws, bssm_approximate_mcmc
   public :: bssm_nonlinear_approximate_mcmc_draws
   public :: bssm_nonlinear_approximate_mcmc
   public :: bssm_nonlinear_ekf_mcmc_draws, bssm_nonlinear_ekf_mcmc
   public :: bssm_importance_post_correction
   public :: bssm_nonlinear_psi_post_correction_draws
   public :: bssm_nonlinear_psi_post_correction
   public :: bssm_multivariate_nonlinear_psi_post_correction_draws
   public :: bssm_multivariate_nonlinear_psi_post_correction
   public :: bssm_corrected_state_moments, bssm_linear_signal_moments
   public :: bssm_corrected_trajectory_draws, bssm_suggest_particles
   public :: bssm_nonlinear_psi_suggest_particles_draws
   public :: bssm_nonlinear_psi_suggest_particles
   public :: bssm_multivariate_nonlinear_psi_suggest_particles_draws
   public :: bssm_multivariate_nonlinear_psi_suggest_particles
   public :: bssm_integrated_autocorrelation_time
   public :: bssm_asymptotic_variance, bssm_mcmc_diagnostics
   public :: bssm_chain_diagnostics, bssm_da_chain_diagnostics
   public :: bssm_nonlinear_pmmh_draws, bssm_nonlinear_pmmh
   public :: bssm_multivariate_nonlinear_pmmh_draws
   public :: bssm_multivariate_nonlinear_pmmh
   public :: bssm_multivariate_nonlinear_da_pmmh_draws
   public :: bssm_multivariate_nonlinear_da_pmmh
   public :: bssm_nonlinear_psi_pmmh_draws, bssm_nonlinear_psi_pmmh
   public :: bssm_multivariate_nonlinear_psi_pmmh_draws
   public :: bssm_multivariate_nonlinear_psi_pmmh
   public :: bssm_multivariate_nonlinear_psi_da_pmmh_draws
   public :: bssm_multivariate_nonlinear_psi_da_pmmh
   public :: bssm_nonlinear_ekpf_pmmh_draws, bssm_nonlinear_ekpf_pmmh
   public :: bssm_multivariate_ekpf_pmmh_draws
   public :: bssm_multivariate_ekpf_pmmh
   public :: bssm_multivariate_ekpf_da_pmmh_draws
   public :: bssm_multivariate_ekpf_da_pmmh
   public :: bssm_nonlinear_da_pmmh_draws, bssm_nonlinear_da_pmmh
   public :: bssm_nonlinear_ekpf_da_pmmh_draws
   public :: bssm_nonlinear_ekpf_da_pmmh
   public :: bssm_particle_smoother
   public :: bssm_psi_particle_smoother
   public :: bssm_post_corrected_particle_moments
   public :: bssm_post_corrected_particle_trajectories_draws
   public :: bssm_post_corrected_particle_trajectories
   public :: bssm_simulation_smoother_draws, bssm_simulation_smoother
   public :: bssm_approximation_simulation_draws
   public :: bssm_approximation_simulation
   public :: bssm_nongaussian_simulation_smoother_draws
   public :: bssm_nongaussian_simulation_smoother
   public :: bssm_multivariate_simulation_smoother_draws
   public :: bssm_multivariate_simulation_smoother
   public :: bssm_predictive_draws, bssm_predictive_sample
   public :: bssm_predictive_replicates_draws
   public :: bssm_predictive_replicates
   public :: bssm_predict_past_draws, bssm_predict_past
   public :: bssm_posterior_observation_model_t
   public :: bssm_nonlinear_predictive_draws
   public :: bssm_nonlinear_predictive_sample
   public :: bssm_nonlinear_predictive_replicates_draws
   public :: bssm_nonlinear_predictive_replicates
   public :: bssm_nonlinear_predict_past_draws
   public :: bssm_nonlinear_predict_past
   public :: bssm_summarize_prediction
   public :: bssm_nonlinear_fitted_summary_draws
   public :: bssm_nonlinear_fitted_summary

contains

   pure function bssm_simulation_smoother_draws(model, initial_normals, &
      observation_normals, state_normals, diffuse_variance) result(out)
      !! Draw Gaussian state paths by Durbin-Koopman simulation smoothing.
      type(ssm_model_t), intent(in) :: model !! Model specification.
      real(dp), intent(in) :: initial_normals(:, :) !! Initial normals.
      real(dp), intent(in) :: observation_normals(:, :, :) !! Observation normals.
      real(dp), intent(in) :: state_normals(:, :, :) !! State normals.
      real(dp), intent(in), optional :: diffuse_variance !! Diffuse variance.
      type(bssm_simulation_smoother_t) :: out
      type(ssm_model_t) :: working_model, correction_model
      type(kfs_filter_t) :: filtered
      type(kfs_smoother_t) :: smoothed
      real(dp), allocatable :: initial_factor(:, :), observation_factor(:, :)
      real(dp), allocatable :: state_factor(:, :), simulated_state(:, :)
      real(dp), allocatable :: simulated_observation(:, :), difference(:)
      real(dp) :: diffuse_scale
      integer :: draws, state, series, times, noise
      integer :: draw, time, component, info, tz, th, tt, tr, tq

      call validate_ssm(model, info)
      if (info /= 0) then
         out%info = 1
         return
      end if
      state = size(model%a1)
      series = size(model%y, 2)
      times = size(model%y, 1)
      noise = size(model%q, 1)
      draws = size(initial_normals, 2)
      diffuse_scale = 1.0e7_dp
      if (present(diffuse_variance)) diffuse_scale = diffuse_variance
      if (draws < 1 .or. diffuse_scale <= 0.0_dp .or. &
         any(shape(initial_normals) /= [state, draws]) .or. &
         any(shape(observation_normals) /= [series, times, draws]) .or. &
         any(shape(state_normals) /= [noise, max(0, times - 1), draws]) .or. &
         .not. all(ieee_is_finite(initial_normals)) .or. &
         .not. all(ieee_is_finite(observation_normals)) .or. &
         .not. all(ieee_is_finite(state_normals))) then
         out%info = 2
         return
      end if

      working_model = model
      if (allocated(working_model%p1inf)) then
         working_model%p1 = working_model%p1 + &
            diffuse_scale*working_model%p1inf
         working_model%p1inf = 0.0_dp
      end if
      call positive_semidefinite_factor(working_model%p1, initial_factor, info)
      if (info /= 0) then
         out%info = 10 + info
         return
      end if
      allocate(out%trajectories(state, times, draws))
      allocate(out%observations(times, series, draws))
      allocate(out%mean(state, times), out%covariance(state, state, times))
      allocate(simulated_state(state, times))
      allocate(simulated_observation(times, series), difference(state))
      out%mean = 0.0_dp
      out%covariance = 0.0_dp

      do draw = 1, draws
         simulated_state(:, 1) = working_model%a1 + &
            matmul(initial_factor, initial_normals(:, draw))
         do time = 1, times
            tz = bssm_time_index(working_model%z, time)
            th = bssm_time_index(working_model%h, time)
            call positive_semidefinite_factor(working_model%h(:, :, th), &
               observation_factor, info)
            if (info /= 0) then
               out%info = 20 + time
               return
            end if
            simulated_observation(time, :) = &
               matmul(working_model%z(:, :, tz), simulated_state(:, time)) + &
               matmul(observation_factor, observation_normals(:, time, draw))
            if (time < times) then
               tt = bssm_time_index(working_model%transition, time)
               tr = bssm_time_index(working_model%r, time)
               tq = bssm_time_index(working_model%q, time)
               call positive_semidefinite_factor(working_model%q(:, :, tq), &
                  state_factor, info)
               if (info /= 0) then
                  out%info = 30 + time
                  return
               end if
               simulated_state(:, time + 1) = &
                  matmul(working_model%transition(:, :, tt), &
                  simulated_state(:, time)) + &
                  matmul(working_model%r(:, :, tr), &
                  matmul(state_factor, state_normals(:, time, draw)))
            end if
         end do

         correction_model = working_model
         correction_model%a1 = 0.0_dp
         do time = 1, times
            do component = 1, series
               if (bssm_is_observed(working_model, time, component)) then
                  correction_model%y(time, component) = &
                     working_model%y(time, component) - &
                     simulated_observation(time, component)
               end if
            end do
         end do
         filtered = kfs_filter(correction_model)
         if (filtered%info /= 0) then
            out%info = 100 + filtered%info
            return
         end if
         smoothed = kfs_smooth(correction_model, filtered)
         if (smoothed%info /= 0) then
            out%info = 200 + smoothed%info
            return
         end if
         out%trajectories(:, :, draw) = simulated_state + smoothed%state
         do time = 1, times
            tz = bssm_time_index(working_model%z, time)
            do component = 1, series
               if (bssm_is_observed(working_model, time, component)) then
                  out%observations(time, component, draw) = &
                     working_model%y(time, component)
               else
                  out%observations(time, component, draw) = &
                     simulated_observation(time, component) + &
                     dot_product(working_model%z(component, :, tz), &
                     smoothed%state(:, time))
               end if
            end do
         end do
         out%mean = out%mean + out%trajectories(:, :, draw)
      end do
      out%mean = out%mean/real(draws, dp)
      do draw = 1, draws
         do time = 1, times
            difference = out%trajectories(:, time, draw) - out%mean(:, time)
            out%covariance(:, :, time) = out%covariance(:, :, time) + &
               outer_product(difference, difference)
         end do
      end do
      out%covariance = out%covariance/real(draws, dp)
   end function bssm_simulation_smoother_draws

   function bssm_simulation_smoother(model, draws, antithetic, &
      diffuse_variance) result(out)
      !! Draw Gaussian smoothed paths using the shared random-number stream.
      type(ssm_model_t), intent(in) :: model !! Model specification.
      integer, intent(in) :: draws !! Draws.
      logical, intent(in), optional :: antithetic !! Flag controlling antithetic.
      real(dp), intent(in), optional :: diffuse_variance !! Diffuse variance.
      type(bssm_simulation_smoother_t) :: out
      real(dp), allocatable :: initial_normals(:, :)
      real(dp), allocatable :: observation_normals(:, :, :)
      real(dp), allocatable :: state_normals(:, :, :)
      logical :: paired
      integer :: state, series, times, noise, draw, pairs, info

      call validate_ssm(model, info)
      paired = .false.
      if (present(antithetic)) paired = antithetic
      if (info /= 0 .or. draws < 1 .or. (paired .and. mod(draws, 2) /= 0)) then
         out%info = 1
         return
      end if
      state = size(model%a1)
      series = size(model%y, 2)
      times = size(model%y, 1)
      noise = size(model%q, 1)
      allocate(initial_normals(state, draws))
      allocate(observation_normals(series, times, draws))
      allocate(state_normals(noise, max(0, times - 1), draws))
      pairs = draws
      if (paired) pairs = draws/2
      call random_standard_normal_matrix(initial_normals(:, 1:pairs))
      do draw = 1, pairs
         call random_standard_normal_matrix(observation_normals(:, :, draw))
         if (times > 1) then
            call random_standard_normal_matrix(state_normals(:, :, draw))
         end if
      end do
      if (paired) then
         initial_normals(:, pairs + 1:draws) = -initial_normals(:, 1:pairs)
         observation_normals(:, :, pairs + 1:draws) = &
            -observation_normals(:, :, 1:pairs)
         state_normals(:, :, pairs + 1:draws) = &
            -state_normals(:, :, 1:pairs)
      end if
      if (present(diffuse_variance)) then
         out = bssm_simulation_smoother_draws(model, initial_normals, &
            observation_normals, state_normals, diffuse_variance)
      else
         out = bssm_simulation_smoother_draws(model, initial_normals, &
            observation_normals, state_normals)
      end if
   end function bssm_simulation_smoother

   pure function bssm_scalar_approximation_simulation_draws(approximation, &
      normals) result(out)
      !! Draw paths from a scalar-observation Gaussian approximation.
      type(bssm_gaussian_approximation_t), intent(in) :: approximation !! Approximation.
      real(dp), intent(in) :: normals(:, :, :) !! Independent standard-normal draws.
      type(bssm_simulation_smoother_t) :: out

      if (approximation%info /= 0 .or. &
         .not. allocated(approximation%mode_state) .or. &
         .not. allocated(approximation%proposal_factor) .or. &
         .not. allocated(approximation%conditional_matrix)) then
         out%info = 1
         return
      end if
      out = bssm_conditional_simulation_draws(approximation%mode_state, &
         approximation%proposal_factor, approximation%conditional_matrix, &
         normals)
   end function bssm_scalar_approximation_simulation_draws

   pure function bssm_multivariate_approximation_simulation_draws( &
      approximation, normals) result(out)
      !! Draw paths from a mixed-family Gaussian approximation.
      type(bssm_multivariate_approximation_t), intent(in) :: approximation !! Approximation.
      real(dp), intent(in) :: normals(:, :, :) !! Independent standard-normal draws.
      type(bssm_simulation_smoother_t) :: out

      if (approximation%info /= 0 .or. &
         .not. allocated(approximation%mode_state) .or. &
         .not. allocated(approximation%proposal_factor) .or. &
         .not. allocated(approximation%conditional_matrix)) then
         out%info = 1
         return
      end if
      out = bssm_conditional_simulation_draws(approximation%mode_state, &
         approximation%proposal_factor, approximation%conditional_matrix, &
         normals)
   end function bssm_multivariate_approximation_simulation_draws

   function bssm_scalar_approximation_simulation(approximation, samples, &
      use_antithetic) result(out)
      !! Randomly draw paths from a scalar-observation approximation.
      type(bssm_gaussian_approximation_t), intent(in) :: approximation !! Approximation.
      integer, intent(in) :: samples !! Samples.
      logical, intent(in), optional :: use_antithetic !! Whether to use the antithetic.
      type(bssm_simulation_smoother_t) :: out
      real(dp), allocatable :: normals(:, :, :)
      logical :: antithetic

      if (.not. allocated(approximation%mode_state) .or. samples < 1) then
         out%info = 1
         return
      end if
      antithetic = .true.
      if (present(use_antithetic)) antithetic = use_antithetic
      call bssm_random_approximation_normals( &
         size(approximation%mode_state, 1), samples, &
         size(approximation%mode_state, 2), antithetic, normals)
      out = bssm_scalar_approximation_simulation_draws(approximation, normals)
   end function bssm_scalar_approximation_simulation

   function bssm_multivariate_approximation_simulation(approximation, &
      samples, use_antithetic) result(out)
      !! Randomly draw paths from a mixed-family approximation.
      type(bssm_multivariate_approximation_t), intent(in) :: approximation !! Approximation.
      integer, intent(in) :: samples !! Samples.
      logical, intent(in), optional :: use_antithetic !! Whether to use the antithetic.
      type(bssm_simulation_smoother_t) :: out
      real(dp), allocatable :: normals(:, :, :)
      logical :: antithetic

      if (.not. allocated(approximation%mode_state) .or. samples < 1) then
         out%info = 1
         return
      end if
      antithetic = .true.
      if (present(use_antithetic)) antithetic = use_antithetic
      call bssm_random_approximation_normals( &
         size(approximation%mode_state, 1), samples, &
         size(approximation%mode_state, 2), antithetic, normals)
      out = bssm_multivariate_approximation_simulation_draws(approximation, &
         normals)
   end function bssm_multivariate_approximation_simulation

   pure function bssm_nongaussian_simulation_smoother_draws(y, &
      observation_loading, transition, state_noise_loading, initial_mean, &
      initial_covariance, distribution, phi, approximation_iterations, &
      convergence_tolerance, normals, offset, auxiliary, state_offset, &
      initial_mode) result(out)
      !! Approximate non-Gaussian smoothing using supplied normal draws.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: observation_loading(:, :) !! Observation loading matrix.
      real(dp), intent(in) :: transition(:, :, :) !! State transition matrix.
      real(dp), intent(in) :: state_noise_loading(:, :, :) !! State noise loading.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      integer, intent(in) :: distribution !! Probability-distribution specification.
      integer, intent(in) :: approximation_iterations !! Number of approximation iterations.
      real(dp), intent(in) :: phi !! Autoregressive or model coefficient.
      real(dp), intent(in) :: convergence_tolerance !! Convergence tolerance.
      real(dp), intent(in) :: normals(:, :, :) !! Independent standard-normal draws.
      real(dp), intent(in), optional :: offset(:) !! Known additive offset.
      real(dp), intent(in), optional :: auxiliary(:) !! Auxiliary.
      real(dp), intent(in), optional :: state_offset(:, :) !! State offset.
      real(dp), intent(in), optional :: initial_mode(:) !! Initial mode.
      type(bssm_simulation_smoother_t) :: out
      type(bssm_gaussian_approximation_t) :: approximation

      approximation = bssm_gaussian_approximation(y, observation_loading, &
         transition, state_noise_loading, initial_mean, initial_covariance, &
         distribution, phi, approximation_iterations, convergence_tolerance, &
         offset, auxiliary, state_offset, initial_mode)
      if (approximation%info /= 0) then
         out%info = 1000 + approximation%info
         return
      end if
      out = bssm_scalar_approximation_simulation_draws(approximation, normals)
   end function bssm_nongaussian_simulation_smoother_draws

   function bssm_nongaussian_simulation_smoother(y, observation_loading, &
      transition, state_noise_loading, initial_mean, initial_covariance, &
      distribution, phi, samples, approximation_iterations, &
      convergence_tolerance, use_antithetic, offset, auxiliary, state_offset, &
      initial_mode) result(out)
      !! Approximate non-Gaussian smoothing using shared random draws.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: observation_loading(:, :) !! Observation loading matrix.
      real(dp), intent(in) :: transition(:, :, :) !! State transition matrix.
      real(dp), intent(in) :: state_noise_loading(:, :, :) !! State noise loading.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      integer, intent(in) :: distribution !! Probability-distribution specification.
      integer, intent(in) :: samples !! Samples.
      integer, intent(in) :: approximation_iterations !! Number of approximation iterations.
      real(dp), intent(in) :: phi !! Autoregressive or model coefficient.
      real(dp), intent(in) :: convergence_tolerance !! Convergence tolerance.
      logical, intent(in), optional :: use_antithetic !! Whether to use the antithetic.
      real(dp), intent(in), optional :: offset(:) !! Known additive offset.
      real(dp), intent(in), optional :: auxiliary(:) !! Auxiliary.
      real(dp), intent(in), optional :: state_offset(:, :) !! State offset.
      real(dp), intent(in), optional :: initial_mode(:) !! Initial mode.
      type(bssm_simulation_smoother_t) :: out
      type(bssm_gaussian_approximation_t) :: approximation

      approximation = bssm_gaussian_approximation(y, observation_loading, &
         transition, state_noise_loading, initial_mean, initial_covariance, &
         distribution, phi, approximation_iterations, convergence_tolerance, &
         offset, auxiliary, state_offset, initial_mode)
      if (approximation%info /= 0) then
         out%info = 1000 + approximation%info
         return
      end if
      out = bssm_scalar_approximation_simulation(approximation, samples, &
         use_antithetic)
   end function bssm_nongaussian_simulation_smoother

   pure function bssm_multivariate_simulation_smoother_draws(y, &
      observation_loading, transition, state_noise_loading, initial_mean, &
      initial_covariance, distribution, phi, approximation_iterations, &
      convergence_tolerance, normals, offset, auxiliary, state_offset, &
      initial_mode) result(out)
      !! Approximate mixed-family smoothing using supplied normal draws.
      real(dp), intent(in) :: y(:, :) !! Response or time-series observations.
      real(dp), intent(in) :: observation_loading(:, :, :) !! Observation loading matrix.
      real(dp), intent(in) :: transition(:, :, :) !! State transition matrix.
      real(dp), intent(in) :: state_noise_loading(:, :, :) !! State noise loading.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      integer, intent(in) :: distribution(:) !! Probability-distribution specification.
      integer, intent(in) :: approximation_iterations !! Number of approximation iterations.
      real(dp), intent(in) :: phi(:) !! Autoregressive or model coefficient.
      real(dp), intent(in) :: convergence_tolerance !! Convergence tolerance.
      real(dp), intent(in) :: normals(:, :, :) !! Independent standard-normal draws.
      real(dp), intent(in), optional :: offset(:, :) !! Known additive offset.
      real(dp), intent(in), optional :: auxiliary(:, :) !! Auxiliary.
      real(dp), intent(in), optional :: state_offset(:, :) !! State offset.
      real(dp), intent(in), optional :: initial_mode(:, :) !! Initial mode.
      type(bssm_simulation_smoother_t) :: out
      type(bssm_multivariate_approximation_t) :: approximation

      approximation = bssm_multivariate_gaussian_approximation(y, &
         observation_loading, transition, state_noise_loading, initial_mean, &
         initial_covariance, distribution, phi, approximation_iterations, &
         convergence_tolerance, offset, auxiliary, state_offset, initial_mode)
      if (approximation%info /= 0) then
         out%info = 1000 + approximation%info
         return
      end if
      out = bssm_multivariate_approximation_simulation_draws(approximation, &
         normals)
   end function bssm_multivariate_simulation_smoother_draws

   function bssm_multivariate_simulation_smoother(y, observation_loading, &
      transition, state_noise_loading, initial_mean, initial_covariance, &
      distribution, phi, samples, approximation_iterations, &
      convergence_tolerance, use_antithetic, offset, auxiliary, state_offset, &
      initial_mode) result(out)
      !! Approximate mixed-family smoothing using shared random draws.
      real(dp), intent(in) :: y(:, :) !! Response or time-series observations.
      real(dp), intent(in) :: observation_loading(:, :, :) !! Observation loading matrix.
      real(dp), intent(in) :: transition(:, :, :) !! State transition matrix.
      real(dp), intent(in) :: state_noise_loading(:, :, :) !! State noise loading.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      integer, intent(in) :: distribution(:) !! Probability-distribution specification.
      integer, intent(in) :: samples !! Samples.
      integer, intent(in) :: approximation_iterations !! Number of approximation iterations.
      real(dp), intent(in) :: phi(:) !! Autoregressive or model coefficient.
      real(dp), intent(in) :: convergence_tolerance !! Convergence tolerance.
      logical, intent(in), optional :: use_antithetic !! Whether to use the antithetic.
      real(dp), intent(in), optional :: offset(:, :) !! Known additive offset.
      real(dp), intent(in), optional :: auxiliary(:, :) !! Auxiliary.
      real(dp), intent(in), optional :: state_offset(:, :) !! State offset.
      real(dp), intent(in), optional :: initial_mode(:, :) !! Initial mode.
      type(bssm_simulation_smoother_t) :: out
      type(bssm_multivariate_approximation_t) :: approximation

      approximation = bssm_multivariate_gaussian_approximation(y, &
         observation_loading, transition, state_noise_loading, initial_mean, &
         initial_covariance, distribution, phi, approximation_iterations, &
         convergence_tolerance, offset, auxiliary, state_offset, initial_mode)
      if (approximation%info /= 0) then
         out%info = 1000 + approximation%info
         return
      end if
      out = bssm_multivariate_approximation_simulation(approximation, samples, &
         use_antithetic)
   end function bssm_multivariate_simulation_smoother

   pure function bssm_predictive_draws(initial_state, observation_loading, &
      transition, state_noise_loading, distribution, phi, state_normals, &
      response_normals, response_uniforms, offset, auxiliary, state_offset, &
      observation_noise_loading) result(out)
      !! Draw future states and responses from supplied random variates.
      real(dp), intent(in) :: initial_state(:, :) !! Initial state vector.
      real(dp), intent(in) :: observation_loading(:, :, :) !! Observation loading matrix.
      real(dp), intent(in) :: transition(:, :, :) !! State transition matrix.
      real(dp), intent(in) :: state_noise_loading(:, :, :) !! State noise loading.
      integer, intent(in) :: distribution(:) !! Probability-distribution specification.
      real(dp), intent(in) :: phi(:) !! Autoregressive or model coefficient.
      real(dp), intent(in) :: state_normals(:, :, :) !! State normals.
      real(dp), intent(in) :: response_normals(:, :, :) !! Response normals.
      real(dp), intent(in) :: response_uniforms(:, :, :) !! Response uniforms.
      real(dp), intent(in), optional :: offset(:, :) !! Known additive offset.
      real(dp), intent(in), optional :: auxiliary(:, :) !! Auxiliary.
      real(dp), intent(in), optional :: state_offset(:, :) !! State offset.
      real(dp), intent(in), optional :: observation_noise_loading(:, :, :) !! Observation noise loading.
      type(bssm_prediction_t) :: out
      real(dp), allocatable :: offset_work(:, :), auxiliary_work(:, :)
      real(dp), allocatable :: state_offset_work(:, :)
      integer :: state, series, noise, horizon, samples
      integer :: sample, time, tt, tr

      state = size(initial_state, 1)
      samples = size(initial_state, 2)
      series = size(observation_loading, 1)
      horizon = size(response_normals, 2)
      noise = size(state_noise_loading, 2)
      if (state < 1 .or. samples < 1 .or. series < 1 .or. horizon < 1 .or. &
         noise < 1 .or. size(observation_loading, 2) /= state .or. &
         size(transition, 1) /= state .or. size(transition, 2) /= state .or. &
         size(state_noise_loading, 1) /= state .or. &
         size(distribution) /= series .or. size(phi) /= series .or. &
         .not. bssm_valid_prediction_extent(size(observation_loading, 3), &
         horizon) .or. &
         .not. bssm_valid_prediction_extent(size(transition, 3), &
         max(1, horizon - 1)) .or. &
         .not. bssm_valid_prediction_extent( &
         size(state_noise_loading, 3), max(1, horizon - 1)) .or. &
         any(shape(state_normals) /= &
         [noise, max(0, horizon - 1), samples]) .or. &
         any(shape(response_normals) /= [series, horizon, samples]) .or. &
         any(shape(response_uniforms) /= [series, horizon, samples]) .or. &
         any(distribution < bssm_svm) .or. &
         any(distribution > bssm_gaussian) .or. &
         .not. all(ieee_is_finite(initial_state)) .or. &
         .not. all(ieee_is_finite(state_normals)) .or. &
         .not. all(ieee_is_finite(response_normals)) .or. &
         .not. all(ieee_is_finite(response_uniforms)) .or. &
         any(response_uniforms < 0.0_dp) .or. &
         any(response_uniforms > 1.0_dp)) then
         out%info = 1
         return
      end if
      call bssm_prediction_options(series, state, horizon, offset, auxiliary, &
         state_offset, offset_work, auxiliary_work, state_offset_work, &
         out%info)
      if (out%info /= 0) return
      allocate(out%state(state, horizon, samples))
      out%state(:, 1, :) = initial_state
      do sample = 1, samples
         do time = 1, horizon - 1
            tt = min(time, size(transition, 3))
            tr = min(time, size(state_noise_loading, 3))
            out%state(:, time + 1, sample) = state_offset_work(:, time) + &
               matmul(transition(:, :, tt), out%state(:, time, sample)) + &
               matmul(state_noise_loading(:, :, tr), &
               state_normals(:, time, sample))
         end do
      end do
      call bssm_form_predictive_response(out%state, observation_loading, &
         distribution, phi, response_normals, response_uniforms, offset_work, &
         auxiliary_work, out%signal, out%mean, out%response, out%info, &
         observation_noise_loading)
   end function bssm_predictive_draws

   function bssm_predictive_sample(initial_state, observation_loading, &
      transition, state_noise_loading, distribution, phi, horizon, offset, &
      auxiliary, state_offset, observation_noise_loading) result(out)
      !! Draw future states and responses using the shared random stream.
      real(dp), intent(in) :: initial_state(:, :) !! Initial state vector.
      real(dp), intent(in) :: observation_loading(:, :, :) !! Observation loading matrix.
      real(dp), intent(in) :: transition(:, :, :) !! State transition matrix.
      real(dp), intent(in) :: state_noise_loading(:, :, :) !! State noise loading.
      integer, intent(in) :: distribution(:) !! Probability-distribution specification.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      real(dp), intent(in) :: phi(:) !! Autoregressive or model coefficient.
      real(dp), intent(in), optional :: offset(:, :) !! Known additive offset.
      real(dp), intent(in), optional :: auxiliary(:, :) !! Auxiliary.
      real(dp), intent(in), optional :: state_offset(:, :) !! State offset.
      real(dp), intent(in), optional :: observation_noise_loading(:, :, :) !! Observation noise loading.
      type(bssm_prediction_t) :: out
      real(dp), allocatable :: state_normals(:, :, :)
      real(dp), allocatable :: response_normals(:, :, :)
      real(dp), allocatable :: response_uniforms(:, :, :)
      integer :: sample, samples

      samples = size(initial_state, 2)
      if (horizon < 1 .or. samples < 1) then
         out%info = 1
         return
      end if
      allocate(state_normals(size(state_noise_loading, 2), &
         max(0, horizon - 1), samples))
      allocate(response_normals(size(observation_loading, 1), horizon, &
         samples), response_uniforms(size(observation_loading, 1), horizon, &
         samples))
      do sample = 1, samples
         if (horizon > 1) then
            call random_standard_normal_matrix(state_normals(:, :, sample))
         end if
         call random_standard_normal_matrix(response_normals(:, :, sample))
      end do
      call random_number(response_uniforms)
      out = bssm_predictive_draws(initial_state, observation_loading, &
         transition, state_noise_loading, distribution, phi, state_normals, &
         response_normals, response_uniforms, offset, auxiliary, state_offset, &
         observation_noise_loading)
   end function bssm_predictive_sample

   pure function bssm_predictive_replicates_draws(state, &
      observation_loading, distribution, phi, response_normals, &
      response_uniforms, offset, auxiliary, observation_noise_loading) &
      result(out)
      !! Draw in-sample replicated responses conditional on state paths.
      real(dp), intent(in) :: state(:, :, :) !! State vector or state sequence.
      real(dp), intent(in) :: observation_loading(:, :, :) !! Observation loading matrix.
      integer, intent(in) :: distribution(:) !! Probability-distribution specification.
      real(dp), intent(in) :: phi(:) !! Autoregressive or model coefficient.
      real(dp), intent(in) :: response_normals(:, :, :) !! Response normals.
      real(dp), intent(in) :: response_uniforms(:, :, :) !! Response uniforms.
      real(dp), intent(in), optional :: offset(:, :) !! Known additive offset.
      real(dp), intent(in), optional :: auxiliary(:, :) !! Auxiliary.
      real(dp), intent(in), optional :: observation_noise_loading(:, :, :) !! Observation noise loading.
      type(bssm_prediction_t) :: out
      real(dp), allocatable :: offset_work(:, :), auxiliary_work(:, :)
      real(dp), allocatable :: unused_state_offset(:, :)
      integer :: state_count, series, times, samples

      state_count = size(state, 1)
      times = size(state, 2)
      samples = size(state, 3)
      series = size(observation_loading, 1)
      if (state_count < 1 .or. times < 1 .or. samples < 1 .or. &
         series < 1 .or. size(observation_loading, 2) /= state_count .or. &
         .not. bssm_valid_prediction_extent(size(observation_loading, 3), &
         times) .or. size(distribution) /= series .or. &
         size(phi) /= series .or. &
         any(shape(response_normals) /= [series, times, samples]) .or. &
         any(shape(response_uniforms) /= [series, times, samples]) .or. &
         any(distribution < bssm_svm) .or. &
         any(distribution > bssm_gaussian) .or. &
         .not. all(ieee_is_finite(state)) .or. &
         .not. all(ieee_is_finite(response_normals)) .or. &
         .not. all(ieee_is_finite(response_uniforms)) .or. &
         any(response_uniforms < 0.0_dp) .or. &
         any(response_uniforms > 1.0_dp)) then
         out%info = 1
         return
      end if
      call bssm_prediction_options(series, state_count, times, offset, &
         auxiliary, offset_work=offset_work, &
         auxiliary_work=auxiliary_work, &
         state_offset_work=unused_state_offset, info=out%info)
      if (out%info /= 0) return
      out%state = state
      call bssm_form_predictive_response(out%state, observation_loading, &
         distribution, phi, response_normals, response_uniforms, offset_work, &
         auxiliary_work, out%signal, out%mean, out%response, out%info, &
         observation_noise_loading)
   end function bssm_predictive_replicates_draws

   function bssm_predictive_replicates(state, observation_loading, &
      distribution, phi, offset, auxiliary, observation_noise_loading) &
      result(out)
      !! Draw in-sample replicated responses using the shared random stream.
      real(dp), intent(in) :: state(:, :, :) !! State vector or state sequence.
      real(dp), intent(in) :: observation_loading(:, :, :) !! Observation loading matrix.
      integer, intent(in) :: distribution(:) !! Probability-distribution specification.
      real(dp), intent(in) :: phi(:) !! Autoregressive or model coefficient.
      real(dp), intent(in), optional :: offset(:, :) !! Known additive offset.
      real(dp), intent(in), optional :: auxiliary(:, :) !! Auxiliary.
      real(dp), intent(in), optional :: observation_noise_loading(:, :, :) !! Observation noise loading.
      type(bssm_prediction_t) :: out
      real(dp), allocatable :: response_normals(:, :, :)
      real(dp), allocatable :: response_uniforms(:, :, :)
      integer :: sample, samples, times, series

      samples = size(state, 3)
      times = size(state, 2)
      series = size(observation_loading, 1)
      if (samples < 1 .or. times < 1 .or. series < 1) then
         out%info = 1
         return
      end if
      allocate(response_normals(series, times, samples))
      allocate(response_uniforms(series, times, samples))
      do sample = 1, samples
         call random_standard_normal_matrix(response_normals(:, :, sample))
      end do
      call random_number(response_uniforms)
      out = bssm_predictive_replicates_draws(state, observation_loading, &
         distribution, phi, response_normals, response_uniforms, offset, &
         auxiliary, observation_noise_loading)
   end function bssm_predictive_replicates

   pure function bssm_predict_past_draws(state, parameters, distribution, &
      update_model, response_normals, response_uniforms) result(out)
      !! Reconstruct in-sample observations for paired posterior draws.
      real(dp), intent(in) :: state(:, :, :) !! State vector or state sequence.
      real(dp), intent(in) :: parameters(:, :) !! Model parameter values.
      integer, intent(in) :: distribution(:) !! Probability-distribution specification.
      procedure(bssm_posterior_observation_model_t) :: update_model !! Update model callback procedure.
      real(dp), intent(in) :: response_normals(:, :, :) !! Response normals.
      real(dp), intent(in) :: response_uniforms(:, :, :) !! Response uniforms.
      type(bssm_prediction_t) :: out
      real(dp), allocatable :: loading(:, :, :), phi(:), offset(:, :)
      real(dp), allocatable :: auxiliary(:, :), noise_loading(:, :, :)
      real(dp), allocatable :: signal(:, :, :), mean(:, :, :)
      real(dp), allocatable :: response(:, :, :)
      logical :: correlated_gaussian
      integer :: state_count, series, times, samples, sample, model_info
      integer :: response_info

      state_count = size(state, 1)
      times = size(state, 2)
      samples = size(state, 3)
      series = size(distribution)
      if (state_count < 1 .or. times < 1 .or. samples < 1 .or. series < 1 .or. &
         size(parameters, 1) < 1 .or. size(parameters, 2) /= samples .or. &
         any(shape(response_normals) /= [series, times, samples]) .or. &
         any(shape(response_uniforms) /= [series, times, samples]) .or. &
         any(distribution < bssm_svm) .or. &
         any(distribution > bssm_gaussian) .or. &
         .not. all(ieee_is_finite(state)) .or. &
         .not. all(ieee_is_finite(parameters)) .or. &
         .not. all(ieee_is_finite(response_normals)) .or. &
         .not. all(ieee_is_finite(response_uniforms)) .or. &
         any(response_uniforms < 0.0_dp) .or. &
         any(response_uniforms > 1.0_dp)) then
         out%info = 1
         return
      end if
      allocate(out%state(state_count, times, samples))
      allocate(out%signal(series, times, samples))
      allocate(out%mean(series, times, samples))
      allocate(out%response(series, times, samples))
      allocate(loading(series, state_count, times), phi(series))
      allocate(offset(series, times), auxiliary(series, times))
      allocate(noise_loading(series, series, times))
      out%state = state
      do sample = 1, samples
         call update_model(parameters(:, sample), loading, phi, offset, &
            auxiliary, noise_loading, correlated_gaussian, model_info)
         if (model_info /= 0) then
            out%info = 100 + model_info
            return
         end if
         if (.not. all(ieee_is_finite(loading)) .or. &
            .not. all(ieee_is_finite(phi)) .or. &
            .not. all(ieee_is_finite(offset)) .or. &
            .not. all(ieee_is_finite(auxiliary)) .or. &
            any(auxiliary < 0.0_dp) .or. &
            (correlated_gaussian .and. &
            .not. all(ieee_is_finite(noise_loading)))) then
            out%info = 200 + sample
            return
         end if
         if (correlated_gaussian) then
            call bssm_form_predictive_response(state(:, :, sample:sample), &
               loading, distribution, phi, &
               response_normals(:, :, sample:sample), &
               response_uniforms(:, :, sample:sample), offset, auxiliary, &
               signal, mean, response, response_info, noise_loading)
         else
            call bssm_form_predictive_response(state(:, :, sample:sample), &
               loading, distribution, phi, &
               response_normals(:, :, sample:sample), &
               response_uniforms(:, :, sample:sample), offset, auxiliary, &
               signal, mean, response, response_info)
         end if
         if (response_info /= 0) then
            out%info = 300 + response_info
            return
         end if
         out%signal(:, :, sample) = signal(:, :, 1)
         out%mean(:, :, sample) = mean(:, :, 1)
         out%response(:, :, sample) = response(:, :, 1)
      end do
   end function bssm_predict_past_draws

   function bssm_predict_past(state, parameters, distribution, update_model) &
      result(out)
      !! Reconstruct posterior observations using the shared random stream.
      real(dp), intent(in) :: state(:, :, :) !! State vector or state sequence.
      real(dp), intent(in) :: parameters(:, :) !! Model parameter values.
      integer, intent(in) :: distribution(:) !! Probability-distribution specification.
      procedure(bssm_posterior_observation_model_t) :: update_model !! Update model callback procedure.
      type(bssm_prediction_t) :: out
      real(dp), allocatable :: response_normals(:, :, :)
      real(dp), allocatable :: response_uniforms(:, :, :)
      integer :: series, times, samples, sample

      series = size(distribution)
      times = size(state, 2)
      samples = size(state, 3)
      if (series < 1 .or. times < 1 .or. samples < 1) then
         out%info = 1
         return
      end if
      allocate(response_normals(series, times, samples))
      allocate(response_uniforms(series, times, samples))
      do sample = 1, samples
         call random_standard_normal_matrix(response_normals(:, :, sample))
      end do
      call random_number(response_uniforms)
      out = bssm_predict_past_draws(state, parameters, distribution, &
         update_model, response_normals, response_uniforms)
   end function bssm_predict_past

   pure function bssm_nonlinear_predictive_draws(initial_state, parameters, &
      observation_dimension, state_noise_dimension, observation, &
      transition_model, state_normals, response_normals) result(out)
      !! Draw nonlinear Gaussian forecasts from supplied standard normals.
      real(dp), intent(in) :: initial_state(:, :) !! Initial state vector.
      real(dp), intent(in) :: parameters(:, :) !! Model parameter values.
      integer, intent(in) :: observation_dimension !! Observation dimension.
      integer, intent(in) :: state_noise_dimension !! State noise dimension.
      procedure(bssm_nonlinear_prediction_observation_t) :: observation !! Observed value or vector.
      procedure(bssm_nonlinear_transition_t) :: transition_model !! Transition model callback procedure.
      real(dp), intent(in) :: state_normals(:, :, :) !! State normals.
      real(dp), intent(in) :: response_normals(:, :, :) !! Response normals.
      type(bssm_prediction_t) :: out
      real(dp), allocatable :: transition_mean(:), noise_loading(:, :)
      integer :: state, samples, horizon, sample, time

      state = size(initial_state, 1)
      samples = size(initial_state, 2)
      horizon = size(response_normals, 2)
      if (state < 1 .or. samples < 1 .or. horizon < 1 .or. &
         observation_dimension < 1 .or. state_noise_dimension < 1 .or. &
         size(parameters, 1) < 1 .or. size(parameters, 2) /= samples .or. &
         any(shape(state_normals) /= &
         [state_noise_dimension, max(0, horizon - 1), samples]) .or. &
         any(shape(response_normals) /= &
         [observation_dimension, horizon, samples]) .or. &
         .not. all(ieee_is_finite(initial_state)) .or. &
         .not. all(ieee_is_finite(parameters)) .or. &
         .not. all(ieee_is_finite(state_normals)) .or. &
         .not. all(ieee_is_finite(response_normals))) then
         out%info = 1
         return
      end if
      allocate(out%state(state, horizon, samples), transition_mean(state))
      allocate(noise_loading(state, state_noise_dimension))
      out%state(:, 1, :) = initial_state
      do sample = 1, samples
         do time = 1, horizon - 1
            call transition_model(time, out%state(:, time, sample), &
               parameters(:, sample), transition_mean, noise_loading)
            if (.not. all(ieee_is_finite(transition_mean)) .or. &
               .not. all(ieee_is_finite(noise_loading))) then
               out%info = 10 + time
               return
            end if
            out%state(:, time + 1, sample) = transition_mean + &
               matmul(noise_loading, state_normals(:, time, sample))
         end do
      end do
      call bssm_form_nonlinear_prediction(out%state, parameters, observation, &
         response_normals, out%signal, out%mean, out%response, out%info)
   end function bssm_nonlinear_predictive_draws

   function bssm_nonlinear_predictive_sample(initial_state, parameters, &
      observation_dimension, state_noise_dimension, horizon, observation, &
      transition_model) result(out)
      !! Draw nonlinear Gaussian forecasts using the shared random stream.
      real(dp), intent(in) :: initial_state(:, :) !! Initial state vector.
      real(dp), intent(in) :: parameters(:, :) !! Model parameter values.
      integer, intent(in) :: observation_dimension !! Observation dimension.
      integer, intent(in) :: state_noise_dimension !! State noise dimension.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      procedure(bssm_nonlinear_prediction_observation_t) :: observation !! Observed value or vector.
      procedure(bssm_nonlinear_transition_t) :: transition_model !! Transition model callback procedure.
      type(bssm_prediction_t) :: out
      real(dp), allocatable :: state_normals(:, :, :)
      real(dp), allocatable :: response_normals(:, :, :)
      integer :: sample, samples

      samples = size(initial_state, 2)
      if (samples < 1 .or. horizon < 1 .or. observation_dimension < 1 .or. &
         state_noise_dimension < 1) then
         out%info = 1
         return
      end if
      allocate(state_normals(state_noise_dimension, max(0, horizon - 1), &
         samples), response_normals(observation_dimension, horizon, samples))
      do sample = 1, samples
         if (horizon > 1) then
            call random_standard_normal_matrix(state_normals(:, :, sample))
         end if
         call random_standard_normal_matrix(response_normals(:, :, sample))
      end do
      out = bssm_nonlinear_predictive_draws(initial_state, parameters, &
         observation_dimension, state_noise_dimension, observation, &
         transition_model, state_normals, response_normals)
   end function bssm_nonlinear_predictive_sample

   pure function bssm_nonlinear_predictive_replicates_draws(state, &
      parameters, observation, response_normals) result(out)
      !! Draw nonlinear Gaussian replications conditional on state paths.
      real(dp), intent(in) :: state(:, :, :) !! State vector or state sequence.
      real(dp), intent(in) :: parameters(:, :) !! Model parameter values.
      procedure(bssm_nonlinear_prediction_observation_t) :: observation !! Observed value or vector.
      real(dp), intent(in) :: response_normals(:, :, :) !! Response normals.
      type(bssm_prediction_t) :: out
      integer :: samples

      samples = size(state, 3)
      if (size(state, 1) < 1 .or. size(state, 2) < 1 .or. samples < 1 .or. &
         size(parameters, 1) < 1 .or. size(parameters, 2) /= samples .or. &
         size(response_normals, 1) < 1 .or. &
         size(response_normals, 2) /= size(state, 2) .or. &
         size(response_normals, 3) /= samples .or. &
         .not. all(ieee_is_finite(state)) .or. &
         .not. all(ieee_is_finite(parameters)) .or. &
         .not. all(ieee_is_finite(response_normals))) then
         out%info = 1
         return
      end if
      out%state = state
      call bssm_form_nonlinear_prediction(out%state, parameters, observation, &
         response_normals, out%signal, out%mean, out%response, out%info)
   end function bssm_nonlinear_predictive_replicates_draws

   function bssm_nonlinear_predictive_replicates(state, parameters, &
      observation_dimension, observation) result(out)
      !! Draw nonlinear Gaussian replications using the shared random stream.
      real(dp), intent(in) :: state(:, :, :) !! State vector or state sequence.
      real(dp), intent(in) :: parameters(:, :) !! Model parameter values.
      integer, intent(in) :: observation_dimension !! Observation dimension.
      procedure(bssm_nonlinear_prediction_observation_t) :: observation !! Observed value or vector.
      type(bssm_prediction_t) :: out
      real(dp), allocatable :: response_normals(:, :, :)
      integer :: sample, samples, times

      samples = size(state, 3)
      times = size(state, 2)
      if (samples < 1 .or. times < 1 .or. observation_dimension < 1) then
         out%info = 1
         return
      end if
      allocate(response_normals(observation_dimension, times, samples))
      do sample = 1, samples
         call random_standard_normal_matrix(response_normals(:, :, sample))
      end do
      out = bssm_nonlinear_predictive_replicates_draws(state, parameters, &
         observation, response_normals)
   end function bssm_nonlinear_predictive_replicates

   pure function bssm_nonlinear_predict_past_draws(state, parameters, &
      observation, response_normals) result(out)
      !! Reconstruct nonlinear Gaussian observations from supplied draws.
      real(dp), intent(in) :: state(:, :, :) !! State vector or state sequence.
      real(dp), intent(in) :: parameters(:, :) !! Model parameter values.
      procedure(bssm_nonlinear_prediction_observation_t) :: observation !! Observed value or vector.
      real(dp), intent(in) :: response_normals(:, :, :) !! Response normals.
      type(bssm_prediction_t) :: out

      out = bssm_nonlinear_predictive_replicates_draws(state, parameters, &
         observation, response_normals)
   end function bssm_nonlinear_predict_past_draws

   function bssm_nonlinear_predict_past(state, parameters, &
      observation_dimension, observation) result(out)
      !! Reconstruct nonlinear Gaussian observations from the random stream.
      real(dp), intent(in) :: state(:, :, :) !! State vector or state sequence.
      real(dp), intent(in) :: parameters(:, :) !! Model parameter values.
      integer, intent(in) :: observation_dimension !! Observation dimension.
      procedure(bssm_nonlinear_prediction_observation_t) :: observation !! Observed value or vector.
      type(bssm_prediction_t) :: out

      out = bssm_nonlinear_predictive_replicates(state, parameters, &
         observation_dimension, observation)
   end function bssm_nonlinear_predict_past

   pure function bssm_summarize_prediction(prediction, probabilities, &
      sample_weight, counts) result(out)
      !! Summarize fitted prediction draws with chain and frequency weights.
      type(bssm_prediction_t), intent(in) :: prediction !! Prediction.
      real(dp), intent(in) :: probabilities(:) !! Probability values.
      real(dp), intent(in) :: sample_weight(:) !! Sample weight.
      integer, intent(in), optional :: counts(:) !! Counts.
      type(bssm_prediction_summary_t) :: out
      real(dp), allocatable :: weight(:), values(:), ordered_weight(:)
      real(dp) :: total, denominator, difference
      integer :: series, times, samples, probability_count
      integer :: component, time, sample, probability

      if (prediction%info /= 0 .or. &
         .not. allocated(prediction%signal) .or. &
         .not. allocated(prediction%response)) then
         out%info = 1
         return
      end if
      series = size(prediction%signal, 1)
      times = size(prediction%signal, 2)
      samples = size(prediction%signal, 3)
      probability_count = size(probabilities)
      if (series < 1 .or. times < 1 .or. samples < 1 .or. &
         any(shape(prediction%response) /= [series, times, samples]) .or. &
         probability_count < 1 .or. size(sample_weight) /= samples .or. &
         .not. all(ieee_is_finite(prediction%signal)) .or. &
         .not. all(ieee_is_finite(prediction%response)) .or. &
         .not. all(ieee_is_finite(probabilities)) .or. &
         .not. all(ieee_is_finite(sample_weight)) .or. &
         any(probabilities < 0.0_dp) .or. any(probabilities > 1.0_dp) .or. &
         any(sample_weight < 0.0_dp)) then
         out%info = 1
         return
      end if
      if (present(counts)) then
         if (size(counts) /= samples .or. any(counts < 1)) then
            out%info = 1
            return
         end if
      end if
      allocate(weight(samples))
      weight = sample_weight
      if (present(counts)) weight = weight*real(counts, dp)
      total = sum(weight)
      if (.not. ieee_is_finite(total) .or. total <= 0.0_dp) then
         out%info = 1
         return
      end if
      denominator = total - sum(weight**2)/total
      allocate(out%probabilities(probability_count))
      allocate(out%signal_mean(series, times), out%signal_sd(series, times))
      allocate(out%signal_quantile(series, times, probability_count))
      allocate(out%observation_mean(series, times))
      allocate(out%observation_sd(series, times))
      allocate(out%observation_quantile(series, times, probability_count))
      allocate(values(samples), ordered_weight(samples))
      out%probabilities = probabilities
      do time = 1, times
         do component = 1, series
            values = prediction%signal(component, time, :)
            out%signal_mean(component, time) = dot_product(weight, values)/total
            if (denominator > 0.0_dp) then
               out%signal_sd(component, time) = sqrt(dot_product(weight, &
                  (values - out%signal_mean(component, time))**2)/denominator)
            else
               out%signal_sd(component, time) = 0.0_dp
            end if
            ordered_weight = weight
            call sort_weighted_values(values, ordered_weight)
            do probability = 1, probability_count
               out%signal_quantile(component, time, probability) = &
                  weighted_quantile(values, ordered_weight, &
                  probabilities(probability))
            end do
            values = prediction%response(component, time, :)
            out%observation_mean(component, time) = &
               dot_product(weight, values)/total
            if (denominator > 0.0_dp) then
               difference = dot_product(weight, &
                  (values - out%observation_mean(component, time))**2)
               out%observation_sd(component, time) = &
                  sqrt(difference/denominator)
            else
               out%observation_sd(component, time) = 0.0_dp
            end if
            ordered_weight = weight
            call sort_weighted_values(values, ordered_weight)
            do probability = 1, probability_count
               out%observation_quantile(component, time, probability) = &
                  weighted_quantile(values, ordered_weight, &
                  probabilities(probability))
            end do
         end do
      end do

   contains

      pure subroutine sort_weighted_values(value, value_weight)
         !! Sort values and carry their weights through insertion moves.
         real(dp), intent(inout) :: value(:) !! Input value, updated in place.
         real(dp), intent(inout) :: value_weight(:) !! Value weight, updated in place.
         real(dp) :: current, current_weight
         integer :: index, position

         do index = 2, size(value)
            current = value(index)
            current_weight = value_weight(index)
            position = index - 1
            do while (position >= 1)
               if (value(position) <= current) exit
               value(position + 1) = value(position)
               value_weight(position + 1) = value_weight(position)
               position = position - 1
            end do
            value(position + 1) = current
            value_weight(position + 1) = current_weight
         end do
      end subroutine sort_weighted_values

      pure real(dp) function weighted_quantile(value, value_weight, &
         probability_value) result(quantile)
         !! Return the inverse weighted empirical distribution function.
         real(dp), intent(in) :: value(:) !! Input value.
         real(dp), intent(in) :: value_weight(:) !! Value weight.
         real(dp), intent(in) :: probability_value !! Probability value.
         real(dp) :: threshold, cumulative
         integer :: index

         threshold = probability_value*sum(value_weight)
         cumulative = 0.0_dp
         quantile = value(size(value))
         do index = 1, size(value)
            cumulative = cumulative + value_weight(index)
            if (cumulative >= threshold) then
               quantile = value(index)
               return
            end if
         end do
      end function weighted_quantile

   end function bssm_summarize_prediction

   pure function bssm_nonlinear_fitted_summary_draws(state, parameters, &
      observation, response_normals, probabilities, sample_weight, counts) &
      result(out)
      !! Summarize supplied nonlinear in-sample prediction draws.
      real(dp), intent(in) :: state(:, :, :) !! State vector or state sequence.
      real(dp), intent(in) :: parameters(:, :) !! Model parameter values.
      procedure(bssm_nonlinear_prediction_observation_t) :: observation !! Observed value or vector.
      real(dp), intent(in) :: response_normals(:, :, :) !! Response normals.
      real(dp), intent(in) :: probabilities(:) !! Probability values.
      real(dp), intent(in) :: sample_weight(:) !! Sample weight.
      integer, intent(in), optional :: counts(:) !! Counts.
      type(bssm_prediction_summary_t) :: out
      type(bssm_prediction_t) :: prediction

      prediction = bssm_nonlinear_predict_past_draws(state, parameters, &
         observation, response_normals)
      if (prediction%info /= 0) then
         out%info = 100 + prediction%info
         return
      end if
      if (present(counts)) then
         out = bssm_summarize_prediction(prediction, probabilities, &
            sample_weight, counts)
      else
         out = bssm_summarize_prediction(prediction, probabilities, &
            sample_weight)
      end if
   end function bssm_nonlinear_fitted_summary_draws

   function bssm_nonlinear_fitted_summary(state, parameters, &
      observation_dimension, observation, probabilities, sample_weight, &
      counts) result(out)
      !! Summarize nonlinear fitted draws using the shared random stream.
      real(dp), intent(in) :: state(:, :, :) !! State vector or state sequence.
      real(dp), intent(in) :: parameters(:, :) !! Model parameter values.
      integer, intent(in) :: observation_dimension !! Observation dimension.
      procedure(bssm_nonlinear_prediction_observation_t) :: observation !! Observed value or vector.
      real(dp), intent(in) :: probabilities(:) !! Probability values.
      real(dp), intent(in) :: sample_weight(:) !! Sample weight.
      integer, intent(in), optional :: counts(:) !! Counts.
      type(bssm_prediction_summary_t) :: out
      real(dp), allocatable :: response_normals(:, :, :)
      integer :: samples, sample, times

      samples = size(state, 3)
      times = size(state, 2)
      if (samples < 1 .or. times < 1 .or. observation_dimension < 1) then
         out%info = 1
         return
      end if
      allocate(response_normals(observation_dimension, times, samples))
      do sample = 1, samples
         call random_standard_normal_matrix(response_normals(:, :, sample))
      end do
      if (present(counts)) then
         out = bssm_nonlinear_fitted_summary_draws(state, parameters, &
            observation, response_normals, probabilities, sample_weight, &
            counts)
      else
         out = bssm_nonlinear_fitted_summary_draws(state, parameters, &
            observation, response_normals, probabilities, sample_weight)
      end if
   end function bssm_nonlinear_fitted_summary

   pure elemental real(dp) function bssm_observation_log_density(y, signal, &
      distribution, phi, auxiliary) result(log_density)
      !! Evaluate a supported BSSM observation log density.
      real(dp), intent(in) :: y !! Response or time-series observations.
      real(dp), intent(in) :: signal !! Signal.
      real(dp), intent(in) :: phi !! Autoregressive or model coefficient.
      real(dp), intent(in) :: auxiliary !! Auxiliary.
      integer, intent(in) :: distribution !! Probability-distribution specification.
      real(dp) :: log_mean, log_denominator, mean, variance
      real(dp), parameter :: log_two_pi = log(2.0_dp*acos(-1.0_dp))

      log_density = -huge(1.0_dp)
      if (.not. ieee_is_finite(y) .or. .not. ieee_is_finite(signal) .or. &
         .not. ieee_is_finite(phi) .or. .not. ieee_is_finite(auxiliary)) return
      select case (distribution)
      case (bssm_svm)
         if (phi <= 0.0_dp) return
         if (-signal > log(huge(1.0_dp))) return
         log_density = -0.5_dp*(log_two_pi + 2.0_dp*log(phi) + signal + &
            (y/phi)**2*exp(-signal))
      case (bssm_poisson)
         if (y < 0.0_dp .or. auxiliary <= 0.0_dp .or. &
            abs(y - anint(y)) > 1.0e-8_dp) return
         log_mean = log(auxiliary) + signal
         if (log_mean > log(huge(1.0_dp))) return
         log_density = y*log_mean - exp(log_mean) - log_gamma(y + 1.0_dp)
      case (bssm_binomial)
         if (auxiliary < 0.0_dp .or. y < 0.0_dp .or. y > auxiliary .or. &
            abs(y - anint(y)) > 1.0e-8_dp .or. &
            abs(auxiliary - anint(auxiliary)) > 1.0e-8_dp) return
         log_density = log_gamma(auxiliary + 1.0_dp) - &
            log_gamma(y + 1.0_dp) - log_gamma(auxiliary - y + 1.0_dp) + &
            y*signal - auxiliary*softplus(signal)
      case (bssm_negative_binomial)
         if (y < 0.0_dp .or. phi <= 0.0_dp .or. auxiliary <= 0.0_dp .or. &
            abs(y - anint(y)) > 1.0e-8_dp) return
         log_mean = log(auxiliary) + signal
         log_denominator = log_add_exp(log(phi), log_mean)
         log_density = log_gamma(y + phi) - log_gamma(phi) - &
            log_gamma(y + 1.0_dp) + phi*log(phi) + y*log_mean - &
            (y + phi)*log_denominator
      case (bssm_gamma)
         if (y <= 0.0_dp .or. phi <= 0.0_dp .or. auxiliary <= 0.0_dp) return
         log_mean = log(auxiliary) + signal
         if (log(y) - log_mean > log(huge(1.0_dp))) return
         log_density = phi*log(phi) - log_gamma(phi) + &
            (phi - 1.0_dp)*log(y) - phi*log_mean - &
            phi*exp(log(y) - log_mean)
      case (bssm_gaussian)
         if (phi <= 0.0_dp) return
         mean = signal
         variance = phi*phi
         log_density = -0.5_dp*(log_two_pi + log(variance) + &
            (y - mean)**2/variance)
      end select
   end function bssm_observation_log_density

   pure elemental subroutine bssm_laplace_pseudo_observation(y, signal, &
      distribution, phi, auxiliary, pseudo_observation, variance, info)
      !! Construct one Gaussian pseudo-observation from a Laplace linearization.
      real(dp), intent(in) :: y !! Response or time-series observations.
      real(dp), intent(in) :: signal !! Signal.
      real(dp), intent(in) :: phi !! Autoregressive or model coefficient.
      real(dp), intent(in) :: auxiliary !! Auxiliary.
      integer, intent(in) :: distribution !! Probability-distribution specification.
      real(dp), intent(out) :: pseudo_observation !! Pseudo observation.
      real(dp), intent(out) :: variance !! Variance value or matrix.
      integer, intent(out) :: info !! Status code; zero indicates success.
      real(dp) :: mean, probability, response

      pseudo_observation = 0.0_dp
      variance = 1.0_dp
      info = 0
      if (.not. ieee_is_finite(y)) return
      if (.not. ieee_is_finite(signal) .or. .not. ieee_is_finite(phi) .or. &
         .not. ieee_is_finite(auxiliary)) then
         info = 1
         return
      end if
      select case (distribution)
      case (bssm_svm)
         if (phi <= 0.0_dp) then
            info = 1
            return
         end if
         response = max(abs(y), 1.0e-4_dp)
         variance = 2.0_dp*bounded_exp(signal)/(response/phi)**2
         pseudo_observation = signal + 1.0_dp - 0.5_dp*variance
      case (bssm_poisson)
         if (y < 0.0_dp .or. auxiliary <= 0.0_dp .or. &
            abs(y - anint(y)) > 1.0e-8_dp) then
            info = 1
            return
         end if
         mean = auxiliary*bounded_exp(signal)
         variance = 1.0_dp/mean
         pseudo_observation = signal + (y - mean)/mean
      case (bssm_binomial)
         if (auxiliary <= 0.0_dp .or. y < 0.0_dp .or. y > auxiliary .or. &
            abs(y - anint(y)) > 1.0e-8_dp .or. &
            abs(auxiliary - anint(auxiliary)) > 1.0e-8_dp) then
            info = 1
            return
         end if
         probability = logistic(signal)
         variance = 1.0_dp/(auxiliary*probability*(1.0_dp - probability))
         pseudo_observation = signal + &
            (y - auxiliary*probability)/(auxiliary*probability* &
            (1.0_dp - probability))
      case (bssm_negative_binomial)
         if (y < 0.0_dp .or. phi <= 0.0_dp .or. auxiliary <= 0.0_dp .or. &
            abs(y - anint(y)) > 1.0e-8_dp) then
            info = 1
            return
         end if
         mean = auxiliary*bounded_exp(signal)
         variance = (phi/mean + 2.0_dp + mean/phi)/(y + phi)
         pseudo_observation = signal + &
            (1.0_dp + phi/mean)*(y - mean)/(y + phi)
      case (bssm_gamma)
         if (y <= 0.0_dp .or. phi <= 0.0_dp .or. auxiliary <= 0.0_dp) then
            info = 1
            return
         end if
         mean = auxiliary*bounded_exp(signal)
         variance = mean/(y*phi)
         pseudo_observation = signal - mean/y + 1.0_dp
      case (bssm_gaussian)
         if (phi <= 0.0_dp) then
            info = 1
            return
         end if
         pseudo_observation = y
         variance = phi*phi
      case default
         info = 1
      end select
      if (.not. ieee_is_finite(pseudo_observation) .or. &
         .not. ieee_is_finite(variance) .or. variance <= 0.0_dp) info = 2
   end subroutine bssm_laplace_pseudo_observation

   pure function bssm_stratified_resample(probability, uniforms) result(index)
      !! Draw one ancestor per stratum from normalized probabilities.
      real(dp), intent(in) :: probability(:) !! Probability value.
      real(dp), intent(in) :: uniforms(:) !! Uniforms.
      integer, allocatable :: index(:)
      real(dp), allocatable :: cumulative(:)
      real(dp) :: target, total
      integer :: particle, category, n

      n = size(uniforms)
      allocate(index(n))
      index = 0
      if (n < 1 .or. size(probability) < 1 .or. &
         any(probability < 0.0_dp)) return
      total = sum(probability)
      if (total <= 0.0_dp) return
      allocate(cumulative(size(probability)))
      cumulative(1) = probability(1)/total
      do category = 2, size(probability)
         cumulative(category) = cumulative(category - 1) + &
            probability(category)/total
      end do
      cumulative(size(cumulative)) = 1.0_dp
      category = 1
      do particle = 1, n
         target = (real(particle - 1, dp) + &
            min(max(uniforms(particle), 0.0_dp), 1.0_dp))/real(n, dp)
         do while (category < size(cumulative) .and. &
            target > cumulative(category))
            category = category + 1
         end do
         index(particle) = category
      end do
   end function bssm_stratified_resample

   pure real(dp) function bssm_sde_euler_step(state, step_size, normal, &
      parameters, drift, diffusion, positive) result(next_state)
      !! Advance one scalar SDE Euler-Maruyama substep.
      real(dp), intent(in) :: state !! State vector or state sequence.
      real(dp), intent(in) :: step_size !! Step size.
      real(dp), intent(in) :: normal !! Normal.
      real(dp), intent(in) :: parameters(:) !! Model parameter values.
      logical, intent(in) :: positive !! Flag controlling positive.
      procedure(bssm_sde_coefficient_t) :: drift !! Drift callback procedure.
      procedure(bssm_sde_coefficient_t) :: diffusion !! Diffusion callback procedure.

      next_state = state + drift(state, parameters)*step_size + &
         diffusion(state, parameters)*sqrt(step_size)*normal
      if (positive) next_state = abs(next_state)
   end function bssm_sde_euler_step

   pure real(dp) function bssm_sde_milstein_step(state, step_size, normal, &
      parameters, drift, diffusion, diffusion_derivative, positive) &
      result(next_state)
      !! Advance one scalar SDE Milstein substep.
      real(dp), intent(in) :: state !! State vector or state sequence.
      real(dp), intent(in) :: step_size !! Step size.
      real(dp), intent(in) :: normal !! Normal.
      real(dp), intent(in) :: parameters(:) !! Model parameter values.
      logical, intent(in) :: positive !! Flag controlling positive.
      procedure(bssm_sde_coefficient_t) :: drift !! Drift callback procedure.
      procedure(bssm_sde_coefficient_t) :: diffusion !! Diffusion callback procedure.
      procedure(bssm_sde_coefficient_t) :: diffusion_derivative !! Diffusion derivative callback procedure.
      real(dp) :: scale

      scale = diffusion(state, parameters)
      next_state = state + drift(state, parameters)*step_size + &
         scale*sqrt(step_size)*normal + 0.5_dp*scale* &
         diffusion_derivative(state, parameters)*step_size*(normal**2 - 1.0_dp)
      if (positive) next_state = abs(next_state)
   end function bssm_sde_milstein_step

   pure function bssm_sde_bootstrap_filter_draws(y, initial_state, &
      parameters, substeps, log_density, drift, diffusion, positive, &
      brownian_normals, resampling_uniforms, diffusion_derivative) result(out)
      !! Run a draw-driven bootstrap filter for a scalar continuous SDE.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: initial_state !! Initial state vector.
      real(dp), intent(in) :: parameters(:) !! Model parameter values.
      integer, intent(in) :: substeps !! Substeps.
      logical, intent(in) :: positive !! Flag controlling positive.
      real(dp), intent(in) :: brownian_normals(:, :, :) !! Brownian normals.
      real(dp), intent(in) :: resampling_uniforms(:, :) !! Resampling uniforms.
      procedure(bssm_nonlinear_log_density_t) :: log_density !! Log-density value.
      procedure(bssm_sde_coefficient_t) :: drift !! Drift callback procedure.
      procedure(bssm_sde_coefficient_t) :: diffusion !! Diffusion callback procedure.
      procedure(bssm_sde_coefficient_t), optional :: diffusion_derivative !! Diffusion derivative callback procedure.
      type(bssm_particle_filter_t) :: out
      real(dp), allocatable :: log_weight(:), probability(:), parent(:)
      real(dp) :: increment, step_size, value
      integer :: particles, times, particle, time, step, info

      particles = size(brownian_normals, 1)
      times = size(y)
      if (particles < 1 .or. times < 1 .or. substeps < 1 .or. &
         any(shape(brownian_normals) /= [particles, substeps, times + 1]) .or. &
         any(shape(resampling_uniforms) /= [particles, times]) .or. &
         .not. ieee_is_finite(initial_state) .or. &
         .not. all(ieee_is_finite(parameters)) .or. &
         .not. all(ieee_is_finite(brownian_normals)) .or. &
         .not. all(ieee_is_finite(resampling_uniforms))) then
         out%info = 1
         return
      end if
      allocate(out%particles(1, particles, times + 1))
      allocate(out%weights(particles, times), out%ancestors(particles, times))
      allocate(out%predicted_mean(1, times + 1))
      allocate(out%filtered_mean(1, times))
      allocate(out%predicted_covariance(1, 1, times + 1))
      allocate(out%filtered_covariance(1, 1, times))
      allocate(log_weight(particles), probability(particles), parent(particles))
      step_size = 1.0_dp/real(substeps, dp)
      do particle = 1, particles
         value = initial_state
         do step = 1, substeps
            if (present(diffusion_derivative)) then
               value = bssm_sde_milstein_step(value, step_size, &
                  brownian_normals(particle, step, 1), parameters, drift, &
                  diffusion, diffusion_derivative, positive)
            else
               value = bssm_sde_euler_step(value, step_size, &
                  brownian_normals(particle, step, 1), parameters, drift, &
                  diffusion, positive)
            end if
            if (.not. ieee_is_finite(value)) then
               out%info = 2
               return
            end if
         end do
         out%particles(1, particle, 1) = value
      end do
      do time = 1, times
         call particle_summary(out%particles(:, :, time), &
            out%predicted_mean(:, time), &
            out%predicted_covariance(:, :, time))
         if (ieee_is_finite(y(time))) then
            do particle = 1, particles
               log_weight(particle) = log_density(time, y(time), &
                  out%particles(:, particle, time), parameters)
            end do
            call normalize_log_weights(log_weight, probability, increment, &
               info)
            if (info /= 0) then
               out%info = 3
               out%log_likelihood = -huge(1.0_dp)
               return
            end if
            out%log_likelihood = out%log_likelihood + increment
         else
            probability = 1.0_dp/real(particles, dp)
         end if
         out%weights(:, time) = probability
         call weighted_particle_summary(out%particles(:, :, time), &
            probability, out%filtered_mean(:, time), &
            out%filtered_covariance(:, :, time))
         out%ancestors(:, time) = bssm_stratified_resample(probability, &
            resampling_uniforms(:, time))
         do particle = 1, particles
            parent(particle) = &
               out%particles(1, out%ancestors(particle, time), time)
            value = parent(particle)
            do step = 1, substeps
               if (present(diffusion_derivative)) then
                  value = bssm_sde_milstein_step(value, step_size, &
                     brownian_normals(particle, step, time + 1), parameters, &
                     drift, diffusion, diffusion_derivative, positive)
               else
                  value = bssm_sde_euler_step(value, step_size, &
                     brownian_normals(particle, step, time + 1), parameters, &
                     drift, diffusion, positive)
               end if
               if (.not. ieee_is_finite(value)) then
                  out%info = 2
                  return
               end if
            end do
            out%particles(1, particle, time + 1) = value
         end do
      end do
      call particle_summary(out%particles(:, :, times + 1), &
         out%predicted_mean(:, times + 1), &
         out%predicted_covariance(:, :, times + 1))
   end function bssm_sde_bootstrap_filter_draws

   function bssm_sde_bootstrap_filter(y, initial_state, parameters, substeps, &
      log_density, drift, diffusion, positive, particles, &
      diffusion_derivative) result(out)
      !! Run an SDE bootstrap filter using the shared random stream.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: initial_state !! Initial state vector.
      real(dp), intent(in) :: parameters(:) !! Model parameter values.
      integer, intent(in) :: substeps !! Substeps.
      integer, intent(in) :: particles !! Number of particles.
      logical, intent(in) :: positive !! Flag controlling positive.
      procedure(bssm_nonlinear_log_density_t) :: log_density !! Log-density value.
      procedure(bssm_sde_coefficient_t) :: drift !! Drift callback procedure.
      procedure(bssm_sde_coefficient_t) :: diffusion !! Diffusion callback procedure.
      procedure(bssm_sde_coefficient_t), optional :: diffusion_derivative !! Diffusion derivative callback procedure.
      type(bssm_particle_filter_t) :: out
      real(dp), allocatable :: brownian_normals(:, :, :), uniforms(:, :)
      integer :: interval

      if (particles < 1 .or. substeps < 1 .or. size(y) < 1) then
         out%info = 1
         return
      end if
      allocate(brownian_normals(particles, substeps, size(y) + 1))
      allocate(uniforms(particles, size(y)))
      do interval = 1, size(y) + 1
         call random_standard_normal_matrix(brownian_normals(:, :, interval))
      end do
      call random_number(uniforms)
      if (present(diffusion_derivative)) then
         out = bssm_sde_bootstrap_filter_draws(y, initial_state, parameters, &
            substeps, log_density, drift, diffusion, positive, &
            brownian_normals, uniforms, diffusion_derivative)
      else
         out = bssm_sde_bootstrap_filter_draws(y, initial_state, parameters, &
            substeps, log_density, drift, diffusion, positive, &
            brownian_normals, uniforms)
      end if
   end function bssm_sde_bootstrap_filter

   pure function bssm_sde_is2_state_sampler_draws(y, initial_state, &
      parameters, substeps, approximate_log_likelihood, log_density, drift, &
      diffusion, positive, brownian_normals, resampling_uniforms, &
      terminal_uniforms, diffusion_derivative) result(out)
      !! Draw SDE paths and form fine-versus-approximate IS2 weights.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: initial_state !! Initial state vector.
      real(dp), intent(in) :: parameters(:, :) !! Model parameter values.
      integer, intent(in) :: substeps !! Substeps.
      real(dp), intent(in) :: approximate_log_likelihood(:) !! Approximate log likelihood.
      logical, intent(in) :: positive !! Flag controlling positive.
      procedure(bssm_nonlinear_log_density_t) :: log_density !! Log-density value.
      procedure(bssm_sde_coefficient_t) :: drift !! Drift callback procedure.
      procedure(bssm_sde_coefficient_t) :: diffusion !! Diffusion callback procedure.
      real(dp), intent(in) :: brownian_normals(:, :, :, :) !! Brownian normals.
      real(dp), intent(in) :: resampling_uniforms(:, :, :) !! Resampling uniforms.
      real(dp), intent(in) :: terminal_uniforms(:) !! Terminal uniforms.
      procedure(bssm_sde_coefficient_t), optional :: diffusion_derivative !! Diffusion derivative callback procedure.
      type(bssm_sde_state_sample_t) :: out
      type(bssm_particle_filter_t) :: filter
      real(dp), allocatable :: trajectory(:, :)
      integer :: particles, times, samples, sample, info

      particles = size(brownian_normals, 1)
      times = size(y)
      samples = size(parameters, 2)
      if (particles < 1 .or. times < 1 .or. samples < 1 .or. substeps < 1 .or. &
         size(parameters, 1) < 1 .or. &
         size(approximate_log_likelihood) /= samples .or. &
         any(shape(brownian_normals) /= &
         [particles, substeps, times + 1, samples]) .or. &
         any(shape(resampling_uniforms) /= [particles, times, samples]) .or. &
         size(terminal_uniforms) /= samples .or. &
         .not. ieee_is_finite(initial_state) .or. &
         .not. all(ieee_is_finite(parameters)) .or. &
         .not. all(ieee_is_finite(approximate_log_likelihood)) .or. &
         .not. all(ieee_is_finite(brownian_normals)) .or. &
         .not. all(ieee_is_finite(resampling_uniforms)) .or. &
         .not. all(ieee_is_finite(terminal_uniforms)) .or. &
         any(resampling_uniforms < 0.0_dp) .or. &
         any(resampling_uniforms > 1.0_dp) .or. &
         any(terminal_uniforms < 0.0_dp) .or. &
         any(terminal_uniforms > 1.0_dp)) then
         out%info = 1
         return
      end if
      allocate(out%trajectories(1, times + 1, samples))
      allocate(out%fine_log_likelihood(samples), out%log_weight(samples))
      allocate(out%weight(samples))
      do sample = 1, samples
         if (present(diffusion_derivative)) then
            filter = bssm_sde_bootstrap_filter_draws(y, initial_state, &
               parameters(:, sample), substeps, log_density, drift, diffusion, &
               positive, brownian_normals(:, :, :, sample), &
               resampling_uniforms(:, :, sample), diffusion_derivative)
         else
            filter = bssm_sde_bootstrap_filter_draws(y, initial_state, &
               parameters(:, sample), substeps, log_density, drift, diffusion, &
               positive, brownian_normals(:, :, :, sample), &
               resampling_uniforms(:, :, sample))
         end if
         if (filter%info /= 0) then
            out%info = 100 + filter%info
            return
         end if
         call bssm_trace_resampled_sde_path(filter, terminal_uniforms(sample), &
            trajectory, info)
         if (info /= 0) then
            out%info = 200 + info
            return
         end if
         out%trajectories(:, :, sample) = trajectory
         out%fine_log_likelihood(sample) = filter%log_likelihood
         out%log_weight(sample) = filter%log_likelihood - &
            approximate_log_likelihood(sample)
      end do
      call normalize_log_weights(out%log_weight, out%weight, &
         out%log_mean_weight, info)
      if (info /= 0) then
         out%info = 300 + info
         return
      end if
      out%effective_sample_size = 1.0_dp/sum(out%weight**2)
   end function bssm_sde_is2_state_sampler_draws

   function bssm_sde_is2_state_sampler(y, initial_state, parameters, &
      substeps, approximate_log_likelihood, log_density, drift, diffusion, &
      positive, particles, diffusion_derivative) result(out)
      !! Draw IS2-corrected SDE states using the shared random stream.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: initial_state !! Initial state vector.
      real(dp), intent(in) :: parameters(:, :) !! Model parameter values.
      integer, intent(in) :: substeps !! Substeps.
      integer, intent(in) :: particles !! Number of particles.
      real(dp), intent(in) :: approximate_log_likelihood(:) !! Approximate log likelihood.
      logical, intent(in) :: positive !! Flag controlling positive.
      procedure(bssm_nonlinear_log_density_t) :: log_density !! Log-density value.
      procedure(bssm_sde_coefficient_t) :: drift !! Drift callback procedure.
      procedure(bssm_sde_coefficient_t) :: diffusion !! Diffusion callback procedure.
      procedure(bssm_sde_coefficient_t), optional :: diffusion_derivative !! Diffusion derivative callback procedure.
      type(bssm_sde_state_sample_t) :: out
      real(dp), allocatable :: brownian_normals(:, :, :, :)
      real(dp), allocatable :: resampling_uniforms(:, :, :)
      real(dp), allocatable :: terminal_uniforms(:)
      integer :: samples, sample, interval

      samples = size(parameters, 2)
      if (particles < 1 .or. substeps < 1 .or. size(y) < 1 .or. samples < 1) then
         out%info = 1
         return
      end if
      allocate(brownian_normals(particles, substeps, size(y) + 1, samples))
      allocate(resampling_uniforms(particles, size(y), samples))
      allocate(terminal_uniforms(samples))
      do sample = 1, samples
         do interval = 1, size(y) + 1
            call random_standard_normal_matrix( &
               brownian_normals(:, :, interval, sample))
         end do
         call random_number(resampling_uniforms(:, :, sample))
      end do
      call random_number(terminal_uniforms)
      if (present(diffusion_derivative)) then
         out = bssm_sde_is2_state_sampler_draws(y, initial_state, parameters, &
            substeps, approximate_log_likelihood, log_density, drift, &
            diffusion, positive, brownian_normals, resampling_uniforms, &
            terminal_uniforms, diffusion_derivative)
      else
         out = bssm_sde_is2_state_sampler_draws(y, initial_state, parameters, &
            substeps, approximate_log_likelihood, log_density, drift, &
            diffusion, positive, brownian_normals, resampling_uniforms, &
            terminal_uniforms)
      end if
   end function bssm_sde_is2_state_sampler

   pure function bssm_pmmh_kernel_draws(initial_parameters, prior, estimator, &
      proposal_factor, proposal_normals, acceptance_uniforms, &
      target_acceptance, adaptation_exponent) result(out)
      !! Run a supplied-draw pseudo-marginal Metropolis-Hastings kernel.
      real(dp), intent(in) :: initial_parameters(:) !! Initial parameter values.
      real(dp), intent(in) :: proposal_factor(:, :) !! Proposal factor.
      real(dp), intent(in) :: proposal_normals(:, :) !! Standard-normal proposal draws.
      real(dp), intent(in) :: acceptance_uniforms(:) !! Acceptance uniforms.
      procedure(bssm_parameter_log_density_t) :: prior !! Prior-distribution specification.
      procedure(bssm_likelihood_estimator_t) :: estimator !! Estimator callback procedure.
      real(dp), intent(in), optional :: target_acceptance !! Target acceptance.
      real(dp), intent(in), optional :: adaptation_exponent !! Adaptation exponent.
      type(bssm_mcmc_t) :: out
      type(bssm_likelihood_estimate_t) :: estimate
      real(dp), allocatable :: current(:), proposed(:), scale(:, :)
      real(dp) :: current_likelihood, proposed_likelihood
      real(dp) :: current_prior, proposed_prior, log_ratio
      real(dp) :: target, exponent
      integer :: parameter_count, iterations, iteration, adaptation_info
      integer :: accepted_count
      logical :: adapt

      parameter_count = size(initial_parameters)
      iterations = size(proposal_normals, 2)
      adapt = present(target_acceptance)
      target = 0.234_dp
      exponent = 0.6_dp
      if (adapt) target = target_acceptance
      if (present(adaptation_exponent)) exponent = adaptation_exponent
      if (parameter_count < 1 .or. iterations < 1 .or. &
         any(shape(proposal_factor) /= [parameter_count, parameter_count]) .or. &
         size(proposal_normals, 1) /= parameter_count .or. &
         size(acceptance_uniforms) /= iterations .or. &
         .not. all(ieee_is_finite(initial_parameters)) .or. &
         .not. all(ieee_is_finite(proposal_factor)) .or. &
         .not. all(ieee_is_finite(proposal_normals)) .or. &
         .not. all(ieee_is_finite(acceptance_uniforms)) .or. &
         any(acceptance_uniforms < 0.0_dp) .or. &
         any(acceptance_uniforms > 1.0_dp) .or. &
         (present(adaptation_exponent) .and. .not. adapt) .or. &
         (adapt .and. (.not. ieee_is_finite(target) .or. &
         target <= 0.0_dp .or. target >= 1.0_dp .or. &
         .not. ieee_is_finite(exponent) .or. exponent <= 0.5_dp .or. &
         exponent > 1.0_dp))) then
         out%info = 1
         return
      end if
      allocate(out%parameters(parameter_count, iterations + 1))
      allocate(out%log_likelihood(iterations + 1))
      allocate(out%log_prior(iterations + 1))
      allocate(out%acceptance_probability(iterations), out%accepted(iterations))
      allocate(current(parameter_count), proposed(parameter_count))
      allocate(scale(parameter_count, parameter_count))
      scale = proposal_factor
      current = initial_parameters
      current_prior = prior(current)
      if (.not. ieee_is_finite(current_prior) .or. &
         current_prior <= -0.5_dp*huge(1.0_dp)) then
         out%info = 2
         return
      end if
      estimate = estimator(current, 1)
      if (estimate%info /= 0 .or. .not. ieee_is_finite(estimate%value) .or. &
         estimate%value <= -0.5_dp*huge(1.0_dp)) then
         out%info = 3
         return
      end if
      current_likelihood = estimate%value
      out%parameters(:, 1) = current
      out%log_likelihood(1) = current_likelihood
      out%log_prior(1) = current_prior
      out%accepted = .false.
      out%acceptance_probability = 0.0_dp
      accepted_count = 0
      do iteration = 1, iterations
         proposed = current + matmul(scale, proposal_normals(:, iteration))
         proposed_prior = prior(proposed)
         if (ieee_is_finite(proposed_prior) .and. &
            proposed_prior > -0.5_dp*huge(1.0_dp)) then
            estimate = estimator(proposed, iteration + 1)
            if (estimate%info == 0 .and. ieee_is_finite(estimate%value) .and. &
               estimate%value > -0.5_dp*huge(1.0_dp)) then
               proposed_likelihood = estimate%value
               log_ratio = proposed_likelihood - current_likelihood + &
                  proposed_prior - current_prior
               out%acceptance_probability(iteration) = &
                  exp(min(0.0_dp, log_ratio))
               if (log(max(acceptance_uniforms(iteration), &
                  tiny(1.0_dp))) < min(0.0_dp, log_ratio)) then
                  current = proposed
                  current_prior = proposed_prior
                  current_likelihood = proposed_likelihood
                  out%accepted(iteration) = .true.
                  accepted_count = accepted_count + 1
               end if
            end if
         end if
         if (adapt) then
            call ram_proposal_update(scale, proposal_normals(:, iteration), &
               out%acceptance_probability(iteration), target, iteration, &
               exponent, adaptation_info)
            if (adaptation_info /= 0) then
               out%info = 4
               return
            end if
         end if
         out%parameters(:, iteration + 1) = current
         out%log_likelihood(iteration + 1) = current_likelihood
         out%log_prior(iteration + 1) = current_prior
      end do
      out%acceptance_rate = real(accepted_count, dp)/real(iterations, dp)
      out%final_proposal_factor = scale
   end function bssm_pmmh_kernel_draws

   pure function bssm_da_pmmh_kernel_draws(initial_parameters, prior, &
      coarse_estimator, fine_estimator, proposal_factor, proposal_normals, &
      first_stage_uniforms, second_stage_uniforms, target_acceptance, &
      adaptation_exponent) result(out)
      !! Run a supplied-draw two-estimator delayed-acceptance PMMH kernel.
      real(dp), intent(in) :: initial_parameters(:) !! Initial parameter values.
      real(dp), intent(in) :: proposal_factor(:, :) !! Proposal factor.
      real(dp), intent(in) :: proposal_normals(:, :) !! Standard-normal proposal draws.
      real(dp), intent(in) :: first_stage_uniforms(:) !! First stage uniforms.
      real(dp), intent(in) :: second_stage_uniforms(:) !! Second stage uniforms.
      procedure(bssm_parameter_log_density_t) :: prior !! Prior-distribution specification.
      procedure(bssm_likelihood_estimator_t) :: coarse_estimator !! Coarse estimator callback procedure.
      procedure(bssm_likelihood_estimator_t) :: fine_estimator !! Fine estimator callback procedure.
      real(dp), intent(in), optional :: target_acceptance !! Target acceptance.
      real(dp), intent(in), optional :: adaptation_exponent !! Adaptation exponent.
      type(bssm_da_mcmc_t) :: out
      type(bssm_likelihood_estimate_t) :: estimate
      real(dp), allocatable :: current(:), proposed(:), scale(:, :)
      real(dp) :: current_coarse, current_fine, proposed_coarse, proposed_fine
      real(dp) :: current_prior, proposed_prior, log_ratio
      real(dp) :: target, exponent
      integer :: parameter_count, iterations, iteration, adaptation_info
      integer :: first_stage_count, accepted_count
      logical :: adapt

      parameter_count = size(initial_parameters)
      iterations = size(proposal_normals, 2)
      adapt = present(target_acceptance)
      target = 0.234_dp
      exponent = 0.6_dp
      if (adapt) target = target_acceptance
      if (present(adaptation_exponent)) exponent = adaptation_exponent
      if (parameter_count < 1 .or. iterations < 1 .or. &
         any(shape(proposal_factor) /= [parameter_count, parameter_count]) .or. &
         size(proposal_normals, 1) /= parameter_count .or. &
         size(first_stage_uniforms) /= iterations .or. &
         size(second_stage_uniforms) /= iterations .or. &
         .not. all(ieee_is_finite(initial_parameters)) .or. &
         .not. all(ieee_is_finite(proposal_factor)) .or. &
         .not. all(ieee_is_finite(proposal_normals)) .or. &
         .not. all(ieee_is_finite(first_stage_uniforms)) .or. &
         .not. all(ieee_is_finite(second_stage_uniforms)) .or. &
         any(first_stage_uniforms < 0.0_dp) .or. &
         any(first_stage_uniforms > 1.0_dp) .or. &
         any(second_stage_uniforms < 0.0_dp) .or. &
         any(second_stage_uniforms > 1.0_dp) .or. &
         (present(adaptation_exponent) .and. .not. adapt) .or. &
         (adapt .and. (.not. ieee_is_finite(target) .or. &
         target <= 0.0_dp .or. target >= 1.0_dp .or. &
         .not. ieee_is_finite(exponent) .or. exponent <= 0.5_dp .or. &
         exponent > 1.0_dp))) then
         out%info = 1
         return
      end if
      allocate(out%parameters(parameter_count, iterations + 1))
      allocate(out%coarse_log_likelihood(iterations + 1))
      allocate(out%fine_log_likelihood(iterations + 1))
      allocate(out%log_prior(iterations + 1))
      allocate(out%first_stage_probability(iterations))
      allocate(out%second_stage_probability(iterations))
      allocate(out%first_stage_accepted(iterations), out%accepted(iterations))
      allocate(current(parameter_count), proposed(parameter_count))
      allocate(scale(parameter_count, parameter_count))
      scale = proposal_factor
      current = initial_parameters
      current_prior = prior(current)
      if (.not. valid_log_value(current_prior)) then
         out%info = 2
         return
      end if
      estimate = coarse_estimator(current, 1)
      if (estimate%info /= 0 .or. .not. valid_log_value(estimate%value)) then
         out%info = 3
         return
      end if
      current_coarse = estimate%value
      estimate = fine_estimator(current, 1)
      if (estimate%info /= 0 .or. .not. valid_log_value(estimate%value)) then
         out%info = 3
         return
      end if
      current_fine = estimate%value
      out%parameters(:, 1) = current
      out%coarse_log_likelihood(1) = current_coarse
      out%fine_log_likelihood(1) = current_fine
      out%log_prior(1) = current_prior
      out%first_stage_probability = 0.0_dp
      out%second_stage_probability = 0.0_dp
      out%first_stage_accepted = .false.
      out%accepted = .false.
      first_stage_count = 0
      accepted_count = 0
      do iteration = 1, iterations
         proposed = current + matmul(scale, proposal_normals(:, iteration))
         proposed_prior = prior(proposed)
         if (valid_log_value(proposed_prior)) then
            estimate = coarse_estimator(proposed, iteration + 1)
            if (estimate%info == 0 .and. valid_log_value(estimate%value)) then
               proposed_coarse = estimate%value
               log_ratio = proposed_coarse - current_coarse + &
                  proposed_prior - current_prior
               out%first_stage_probability(iteration) = &
                  exp(min(0.0_dp, log_ratio))
               if (log(max(first_stage_uniforms(iteration), tiny(1.0_dp))) < &
                  min(0.0_dp, log_ratio)) then
                  out%first_stage_accepted(iteration) = .true.
                  first_stage_count = first_stage_count + 1
                  estimate = fine_estimator(proposed, iteration + 1)
                  if (estimate%info == 0 .and. &
                     valid_log_value(estimate%value)) then
                     proposed_fine = estimate%value
                     log_ratio = proposed_fine + current_coarse - &
                        current_fine - proposed_coarse
                     out%second_stage_probability(iteration) = &
                        exp(min(0.0_dp, log_ratio))
                     if (log(max(second_stage_uniforms(iteration), &
                        tiny(1.0_dp))) < min(0.0_dp, log_ratio)) then
                        current = proposed
                        current_prior = proposed_prior
                        current_coarse = proposed_coarse
                        current_fine = proposed_fine
                        out%accepted(iteration) = .true.
                        accepted_count = accepted_count + 1
                     end if
                  end if
               end if
            end if
         end if
         if (adapt) then
            call ram_proposal_update(scale, proposal_normals(:, iteration), &
               out%first_stage_probability(iteration), target, iteration, &
               exponent, adaptation_info)
            if (adaptation_info /= 0) then
               out%info = 4
               return
            end if
         end if
         out%parameters(:, iteration + 1) = current
         out%coarse_log_likelihood(iteration + 1) = current_coarse
         out%fine_log_likelihood(iteration + 1) = current_fine
         out%log_prior(iteration + 1) = current_prior
      end do
      out%first_stage_acceptance_rate = &
         real(first_stage_count, dp)/real(iterations, dp)
      out%acceptance_rate = real(accepted_count, dp)/real(iterations, dp)
      out%final_proposal_factor = scale
   end function bssm_da_pmmh_kernel_draws

   pure function bssm_approximate_mcmc_draws(initial_parameters, prior, &
      approximate_log_likelihood, proposal_factor, proposal_normals, &
      acceptance_uniforms, target_acceptance, adaptation_exponent) result(out)
      !! Run supplied-draw MCMC using a deterministic approximate likelihood.
      real(dp), intent(in) :: initial_parameters(:) !! Initial parameter values.
      real(dp), intent(in) :: proposal_factor(:, :) !! Proposal factor.
      real(dp), intent(in) :: proposal_normals(:, :) !! Standard-normal proposal draws.
      real(dp), intent(in) :: acceptance_uniforms(:) !! Acceptance uniforms.
      procedure(bssm_parameter_log_density_t) :: prior !! Prior-distribution specification.
      procedure(bssm_parameter_log_density_t) :: approximate_log_likelihood !! Approximate log likelihood callback procedure.
      real(dp), intent(in), optional :: target_acceptance !! Target acceptance.
      real(dp), intent(in), optional :: adaptation_exponent !! Adaptation exponent.
      type(bssm_mcmc_t) :: out
      real(dp) :: exponent

      if (present(adaptation_exponent) .and. &
         .not. present(target_acceptance)) then
         out%info = 1
         return
      end if
      exponent = 0.6_dp
      if (present(adaptation_exponent)) exponent = adaptation_exponent
      if (present(target_acceptance)) then
         out = bssm_pmmh_kernel_draws(initial_parameters, prior, estimator, &
            proposal_factor, proposal_normals, acceptance_uniforms, &
            target_acceptance, exponent)
      else
         out = bssm_pmmh_kernel_draws(initial_parameters, prior, estimator, &
            proposal_factor, proposal_normals, acceptance_uniforms)
      end if

   contains

      pure function estimator(parameters, draw_index) result(estimate)
         !! Evaluate the deterministic approximation through the PMMH kernel.
         real(dp), intent(in) :: parameters(:) !! Model parameter values.
         integer, intent(in) :: draw_index !! Index of draw.
         type(bssm_likelihood_estimate_t) :: estimate

         estimate%value = approximate_log_likelihood(parameters) + &
            0.0_dp*real(draw_index, dp)
      end function estimator

   end function bssm_approximate_mcmc_draws

   function bssm_approximate_mcmc(initial_parameters, prior, &
      approximate_log_likelihood, proposal_factor, iterations, &
      target_acceptance, adaptation_exponent) result(out)
      !! Run approximate-likelihood MCMC using the shared random stream.
      real(dp), intent(in) :: initial_parameters(:) !! Initial parameter values.
      real(dp), intent(in) :: proposal_factor(:, :) !! Proposal factor.
      integer, intent(in) :: iterations !! Number of algorithm iterations.
      procedure(bssm_parameter_log_density_t) :: prior !! Prior-distribution specification.
      procedure(bssm_parameter_log_density_t) :: approximate_log_likelihood !! Approximate log likelihood callback procedure.
      real(dp), intent(in), optional :: target_acceptance !! Target acceptance.
      real(dp), intent(in), optional :: adaptation_exponent !! Adaptation exponent.
      type(bssm_mcmc_t) :: out
      real(dp), allocatable :: proposal_normals(:, :), acceptance_uniforms(:)
      real(dp) :: exponent

      if (size(initial_parameters) < 1 .or. iterations < 1 .or. &
         (present(adaptation_exponent) .and. &
         .not. present(target_acceptance))) then
         out%info = 1
         return
      end if
      allocate(proposal_normals(size(initial_parameters), iterations))
      allocate(acceptance_uniforms(iterations))
      call random_standard_normal_matrix(proposal_normals)
      call random_number(acceptance_uniforms)
      exponent = 0.6_dp
      if (present(adaptation_exponent)) exponent = adaptation_exponent
      if (present(target_acceptance)) then
         out = bssm_approximate_mcmc_draws(initial_parameters, prior, &
            approximate_log_likelihood, proposal_factor, proposal_normals, &
            acceptance_uniforms, target_acceptance, exponent)
      else
         out = bssm_approximate_mcmc_draws(initial_parameters, prior, &
            approximate_log_likelihood, proposal_factor, proposal_normals, &
            acceptance_uniforms)
      end if
   end function bssm_approximate_mcmc

   pure function bssm_importance_post_correction(parameters, &
      approximate_log_likelihood, exact_log_likelihood, accepted, scheme) &
      result(out)
      !! Correct an approximate MCMC chain by IS1, IS2, or IS3 weighting.
      real(dp), intent(in) :: parameters(:, :) !! Model parameter values.
      real(dp), intent(in) :: approximate_log_likelihood(:) !! Approximate log likelihood.
      real(dp), intent(in) :: exact_log_likelihood(:, :) !! Exact log likelihood.
      logical, intent(in) :: accepted(:) !! Flag controlling accepted.
      integer, intent(in) :: scheme !! Scheme.
      type(bssm_post_correction_t) :: out
      real(dp), allocatable :: block_log_weight(:), difference(:)
      real(dp) :: block_mean
      integer :: parameter_count, samples, replications
      integer :: sample, block_end, block_size, info

      parameter_count = size(parameters, 1)
      samples = size(parameters, 2)
      replications = size(exact_log_likelihood, 1)
      if (parameter_count < 1 .or. samples < 1 .or. replications < 1 .or. &
         size(approximate_log_likelihood) /= samples .or. &
         size(exact_log_likelihood, 2) /= samples .or. &
         size(accepted) /= samples - 1 .or. &
         (scheme /= bssm_is1 .and. scheme /= bssm_is2 .and. &
         scheme /= bssm_is3) .or. &
         .not. all(ieee_is_finite(parameters)) .or. &
         .not. all(ieee_is_finite(approximate_log_likelihood)) .or. &
         .not. all(ieee_is_finite(exact_log_likelihood))) then
         out%info = 1
         return
      end if
      allocate(out%log_weight(samples), out%weight(samples))
      select case (scheme)
      case (bssm_is3)
         out%log_weight = exact_log_likelihood(1, :) - &
            approximate_log_likelihood
      case (bssm_is1, bssm_is2)
         sample = 1
         do while (sample <= samples)
            block_end = sample
            do while (block_end < samples)
               if (accepted(block_end)) exit
               block_end = block_end + 1
            end do
            block_size = block_end - sample + 1
            if (scheme == bssm_is1) then
               if (replications < block_size) then
                  out%info = 1
                  return
               end if
               allocate(block_log_weight(block_size))
               allocate(difference(block_size))
               block_log_weight = exact_log_likelihood(1:block_size, sample) - &
                  approximate_log_likelihood(sample)
               call normalize_log_weights(block_log_weight, difference, &
                  block_mean, info)
               deallocate(block_log_weight, difference)
               if (info /= 0) then
                  out%info = 2
                  return
               end if
            else
               block_mean = exact_log_likelihood(1, sample) - &
                  approximate_log_likelihood(sample)
            end if
            out%log_weight(sample:block_end) = block_mean
            sample = block_end + 1
         end do
      end select
      call normalize_log_weights(out%log_weight, out%weight, &
         out%log_mean_weight, info)
      if (info /= 0) then
         out%info = 2
         return
      end if
      out%effective_sample_size = 1.0_dp/sum(out%weight**2)
      allocate(out%parameter_mean(parameter_count))
      allocate(out%parameter_covariance(parameter_count, parameter_count))
      out%parameter_mean = matmul(parameters, out%weight)
      out%parameter_covariance = 0.0_dp
      do sample = 1, samples
         out%parameter_covariance = out%parameter_covariance + &
            out%weight(sample)*outer_product( &
            parameters(:, sample) - out%parameter_mean, &
            parameters(:, sample) - out%parameter_mean)
      end do
   end function bssm_importance_post_correction

   pure function bssm_post_correction_from_estimator(parameters, &
      approximate_log_likelihood, accepted, scheme, replications, estimator) &
      result(out)
      !! Evaluate replicated likelihoods and apply an IS correction scheme.
      real(dp), intent(in) :: parameters(:, :) !! Model parameter values.
      real(dp), intent(in) :: approximate_log_likelihood(:) !! Approximate log likelihood.
      logical, intent(in) :: accepted(:) !! Flag controlling accepted.
      integer, intent(in) :: scheme !! Scheme.
      integer, intent(in) :: replications !! Replications.
      procedure(bssm_likelihood_estimator_t) :: estimator !! Estimator callback procedure.
      type(bssm_post_correction_t) :: out
      type(bssm_likelihood_estimate_t) :: estimate
      real(dp), allocatable :: exact_log_likelihood(:, :)
      integer :: samples, sample, replication, draw_index
      integer :: block_end, block_size, evaluations

      samples = size(parameters, 2)
      if (size(parameters, 1) < 1 .or. samples < 1 .or. replications < 1 .or. &
         size(approximate_log_likelihood) /= samples .or. &
         size(accepted) /= samples - 1) then
         out%info = 1
         return
      end if
      allocate(exact_log_likelihood(replications, samples))
      do sample = 1, samples
         exact_log_likelihood(:, sample) = approximate_log_likelihood(sample)
      end do
      select case (scheme)
      case (bssm_is3)
         do sample = 1, samples
            draw_index = (sample - 1)*replications + 1
            estimate = estimator(parameters(:, sample), draw_index)
            if (estimate%info /= 0 .or. &
               .not. ieee_is_finite(estimate%value)) then
               out%info = 100 + estimate%info
               return
            end if
            exact_log_likelihood(1, sample) = estimate%value
         end do
      case (bssm_is1, bssm_is2)
         sample = 1
         do while (sample <= samples)
            block_end = sample
            do while (block_end < samples)
               if (accepted(block_end)) exit
               block_end = block_end + 1
            end do
            block_size = block_end - sample + 1
            evaluations = 1
            if (scheme == bssm_is1) evaluations = block_size
            if (replications < evaluations) then
               out%info = 1
               return
            end if
            do replication = 1, evaluations
               draw_index = (sample - 1)*replications + replication
               estimate = estimator(parameters(:, sample), draw_index)
               if (estimate%info /= 0 .or. &
                  .not. ieee_is_finite(estimate%value)) then
                  out%info = 100 + estimate%info
                  return
               end if
               exact_log_likelihood(replication, sample) = estimate%value
            end do
            sample = block_end + 1
         end do
      case default
         out%info = 1
         return
      end select
      out = bssm_importance_post_correction(parameters, &
         approximate_log_likelihood, exact_log_likelihood, accepted, scheme)
   end function bssm_post_correction_from_estimator

   pure function bssm_nonlinear_psi_post_correction_draws(y, initial_mean, &
      initial_covariance, parameters, approximate_log_likelihood, accepted, &
      scheme, noise_dimension, observation, transition_model, &
      approximation_iterations, convergence_tolerance, initial_normals, &
      state_normals, terminal_normals, resampling_uniforms) result(out)
      !! Post-correct a chain with supplied scalar nonlinear psi-filter draws.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: parameters(:, :) !! Model parameter values.
      real(dp), intent(in) :: approximate_log_likelihood(:) !! Approximate log likelihood.
      logical, intent(in) :: accepted(:) !! Flag controlling accepted.
      integer, intent(in) :: scheme !! Scheme.
      integer, intent(in) :: noise_dimension !! Noise dimension.
      integer, intent(in) :: approximation_iterations !! Number of approximation iterations.
      real(dp), intent(in) :: convergence_tolerance !! Convergence tolerance.
      real(dp), intent(in) :: initial_normals(:, :, :, :) !! Initial normals.
      real(dp), intent(in) :: state_normals(:, :, :, :, :) !! State normals.
      real(dp), intent(in) :: terminal_normals(:, :, :, :) !! Terminal normals.
      real(dp), intent(in) :: resampling_uniforms(:, :, :, :) !! Resampling uniforms.
      procedure(bssm_gaussian_observation_t) :: observation !! Observed value or vector.
      procedure(bssm_nonlinear_transition_jacobian_t) :: transition_model !! Transition model callback procedure.
      type(bssm_post_correction_t) :: out
      integer :: state, particles, times, samples, replications

      state = size(initial_mean)
      particles = size(initial_normals, 2)
      replications = size(initial_normals, 3)
      samples = size(parameters, 2)
      times = size(y)
      if (state < 1 .or. particles < 1 .or. replications < 1 .or. &
         samples < 1 .or. times < 1 .or. noise_dimension < 1 .or. &
         approximation_iterations < 1 .or. convergence_tolerance < 0.0_dp .or. &
         any(shape(initial_covariance) /= [state, state]) .or. &
         any(shape(initial_normals) /= &
         [state, particles, replications, samples]) .or. &
         any(shape(state_normals) /= &
         [state, particles, times, replications, samples]) .or. &
         any(shape(terminal_normals) /= &
         [noise_dimension, particles, replications, samples]) .or. &
         any(shape(resampling_uniforms) /= &
         [particles, times, replications, samples]) .or. &
         .not. all(ieee_is_finite(initial_normals)) .or. &
         .not. all(ieee_is_finite(state_normals)) .or. &
         .not. all(ieee_is_finite(terminal_normals)) .or. &
         .not. all(ieee_is_finite(resampling_uniforms))) then
         out%info = 1
         return
      end if
      out = bssm_post_correction_from_estimator(parameters, &
         approximate_log_likelihood, accepted, scheme, replications, estimator)

   contains

      pure function estimator(parameter, draw_index) result(estimate)
         !! Evaluate one supplied scalar nonlinear psi likelihood.
         real(dp), intent(in) :: parameter(:) !! Parameter.
         integer, intent(in) :: draw_index !! Index of draw.
         type(bssm_likelihood_estimate_t) :: estimate
         type(bssm_particle_filter_t) :: fit
         integer :: sample, replication

         sample = (draw_index - 1)/replications + 1
         replication = modulo(draw_index - 1, replications) + 1
         fit = bssm_nonlinear_psi_filter_draws(y, initial_mean, &
            initial_covariance, parameter, noise_dimension, observation, &
            transition_model, approximation_iterations, &
            convergence_tolerance, &
            initial_normals(:, :, replication, sample), &
            state_normals(:, :, :, replication, sample), &
            terminal_normals(:, :, replication, sample), &
            resampling_uniforms(:, :, replication, sample))
         estimate%value = fit%log_likelihood
         estimate%info = fit%info
      end function estimator

   end function bssm_nonlinear_psi_post_correction_draws

   function bssm_nonlinear_psi_post_correction(y, initial_mean, &
      initial_covariance, parameters, approximate_log_likelihood, accepted, &
      scheme, noise_dimension, observation, transition_model, &
      approximation_iterations, convergence_tolerance, particles, &
      replications) result(out)
      !! Post-correct a scalar nonlinear chain using shared randomness.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: parameters(:, :) !! Model parameter values.
      real(dp), intent(in) :: approximate_log_likelihood(:) !! Approximate log likelihood.
      logical, intent(in) :: accepted(:) !! Flag controlling accepted.
      integer, intent(in) :: scheme !! Scheme.
      integer, intent(in) :: noise_dimension !! Noise dimension.
      integer, intent(in) :: approximation_iterations !! Number of approximation iterations.
      integer, intent(in) :: particles !! Number of particles.
      integer, intent(in) :: replications !! Replications.
      real(dp), intent(in) :: convergence_tolerance !! Convergence tolerance.
      procedure(bssm_gaussian_observation_t) :: observation !! Observed value or vector.
      procedure(bssm_nonlinear_transition_jacobian_t) :: transition_model !! Transition model callback procedure.
      type(bssm_post_correction_t) :: out
      real(dp), allocatable :: initial_normals(:, :, :, :)
      real(dp), allocatable :: state_normals(:, :, :, :, :)
      real(dp), allocatable :: terminal_normals(:, :, :, :)
      real(dp), allocatable :: resampling_uniforms(:, :, :, :)
      integer :: samples, sample, replication, time

      samples = size(parameters, 2)
      if (samples < 1 .or. particles < 1 .or. replications < 1 .or. &
         size(y) < 1 .or. size(initial_mean) < 1 .or. noise_dimension < 1 .or. &
         approximation_iterations < 1 .or. convergence_tolerance < 0.0_dp) then
         out%info = 1
         return
      end if
      allocate(initial_normals(size(initial_mean), particles, replications, &
         samples))
      allocate(state_normals(size(initial_mean), particles, size(y), &
         replications, samples))
      allocate(terminal_normals(noise_dimension, particles, replications, &
         samples))
      allocate(resampling_uniforms(particles, size(y), replications, samples))
      call random_number(resampling_uniforms)
      do sample = 1, samples
         do replication = 1, replications
            call random_standard_normal_matrix( &
               initial_normals(:, :, replication, sample))
            call random_standard_normal_matrix( &
               terminal_normals(:, :, replication, sample))
            do time = 1, size(y)
               call random_standard_normal_matrix( &
                  state_normals(:, :, time, replication, sample))
            end do
         end do
      end do
      out = bssm_nonlinear_psi_post_correction_draws(y, initial_mean, &
         initial_covariance, parameters, approximate_log_likelihood, accepted, &
         scheme, noise_dimension, observation, transition_model, &
         approximation_iterations, convergence_tolerance, initial_normals, &
         state_normals, terminal_normals, resampling_uniforms)
   end function bssm_nonlinear_psi_post_correction

   pure function bssm_multivariate_nonlinear_psi_post_correction_draws(y, &
      initial_mean, initial_covariance, parameters, &
      approximate_log_likelihood, accepted, scheme, &
      observation_noise_dimension, state_noise_dimension, observation, &
      transition_model, approximation_iterations, convergence_tolerance, &
      initial_normals, state_normals, terminal_normals, &
      resampling_uniforms) result(out)
      !! Post-correct a chain with supplied vector nonlinear psi-filter draws.
      real(dp), intent(in) :: y(:, :) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: parameters(:, :) !! Model parameter values.
      real(dp), intent(in) :: approximate_log_likelihood(:) !! Approximate log likelihood.
      logical, intent(in) :: accepted(:) !! Flag controlling accepted.
      integer, intent(in) :: scheme !! Scheme.
      integer, intent(in) :: observation_noise_dimension !! Observation noise dimension.
      integer, intent(in) :: state_noise_dimension !! State noise dimension.
      integer, intent(in) :: approximation_iterations !! Number of approximation iterations.
      real(dp), intent(in) :: convergence_tolerance !! Convergence tolerance.
      real(dp), intent(in) :: initial_normals(:, :, :, :) !! Initial normals.
      real(dp), intent(in) :: state_normals(:, :, :, :, :) !! State normals.
      real(dp), intent(in) :: terminal_normals(:, :, :, :) !! Terminal normals.
      real(dp), intent(in) :: resampling_uniforms(:, :, :, :) !! Resampling uniforms.
      procedure(bssm_multivariate_gaussian_observation_t) :: observation !! Observed value or vector.
      procedure(bssm_nonlinear_transition_jacobian_t) :: transition_model !! Transition model callback procedure.
      type(bssm_post_correction_t) :: out
      integer :: state, particles, times, samples, replications

      state = size(initial_mean)
      particles = size(initial_normals, 2)
      replications = size(initial_normals, 3)
      samples = size(parameters, 2)
      times = size(y, 2)
      if (size(y, 1) < 1 .or. state < 1 .or. particles < 1 .or. &
         replications < 1 .or. samples < 1 .or. times < 1 .or. &
         observation_noise_dimension < 1 .or. state_noise_dimension < 1 .or. &
         approximation_iterations < 1 .or. convergence_tolerance < 0.0_dp .or. &
         any(shape(initial_covariance) /= [state, state]) .or. &
         any(shape(initial_normals) /= &
         [state, particles, replications, samples]) .or. &
         any(shape(state_normals) /= &
         [state, particles, times, replications, samples]) .or. &
         any(shape(terminal_normals) /= &
         [state_noise_dimension, particles, replications, samples]) .or. &
         any(shape(resampling_uniforms) /= &
         [particles, times, replications, samples]) .or. &
         .not. all(ieee_is_finite(initial_normals)) .or. &
         .not. all(ieee_is_finite(state_normals)) .or. &
         .not. all(ieee_is_finite(terminal_normals)) .or. &
         .not. all(ieee_is_finite(resampling_uniforms))) then
         out%info = 1
         return
      end if
      out = bssm_post_correction_from_estimator(parameters, &
         approximate_log_likelihood, accepted, scheme, replications, estimator)

   contains

      pure function estimator(parameter, draw_index) result(estimate)
         !! Evaluate one supplied vector nonlinear psi likelihood.
         real(dp), intent(in) :: parameter(:) !! Parameter.
         integer, intent(in) :: draw_index !! Index of draw.
         type(bssm_likelihood_estimate_t) :: estimate
         type(bssm_particle_filter_t) :: fit
         integer :: sample, replication

         sample = (draw_index - 1)/replications + 1
         replication = modulo(draw_index - 1, replications) + 1
         fit = bssm_multivariate_nonlinear_psi_filter_draws(y, initial_mean, &
            initial_covariance, parameter, observation_noise_dimension, &
            state_noise_dimension, observation, transition_model, &
            approximation_iterations, convergence_tolerance, &
            initial_normals(:, :, replication, sample), &
            state_normals(:, :, :, replication, sample), &
            terminal_normals(:, :, replication, sample), &
            resampling_uniforms(:, :, replication, sample))
         estimate%value = fit%log_likelihood
         estimate%info = fit%info
      end function estimator

   end function bssm_multivariate_nonlinear_psi_post_correction_draws

   function bssm_multivariate_nonlinear_psi_post_correction(y, initial_mean, &
      initial_covariance, parameters, approximate_log_likelihood, accepted, &
      scheme, observation_noise_dimension, state_noise_dimension, &
      observation, transition_model, approximation_iterations, &
      convergence_tolerance, particles, replications) result(out)
      !! Post-correct a vector nonlinear chain using shared randomness.
      real(dp), intent(in) :: y(:, :) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: parameters(:, :) !! Model parameter values.
      real(dp), intent(in) :: approximate_log_likelihood(:) !! Approximate log likelihood.
      logical, intent(in) :: accepted(:) !! Flag controlling accepted.
      integer, intent(in) :: scheme !! Scheme.
      integer, intent(in) :: observation_noise_dimension !! Observation noise dimension.
      integer, intent(in) :: state_noise_dimension !! State noise dimension.
      integer, intent(in) :: approximation_iterations !! Number of approximation iterations.
      integer, intent(in) :: particles !! Number of particles.
      integer, intent(in) :: replications !! Replications.
      real(dp), intent(in) :: convergence_tolerance !! Convergence tolerance.
      procedure(bssm_multivariate_gaussian_observation_t) :: observation !! Observed value or vector.
      procedure(bssm_nonlinear_transition_jacobian_t) :: transition_model !! Transition model callback procedure.
      type(bssm_post_correction_t) :: out
      real(dp), allocatable :: initial_normals(:, :, :, :)
      real(dp), allocatable :: state_normals(:, :, :, :, :)
      real(dp), allocatable :: terminal_normals(:, :, :, :)
      real(dp), allocatable :: resampling_uniforms(:, :, :, :)
      integer :: samples, sample, replication, time

      samples = size(parameters, 2)
      if (samples < 1 .or. particles < 1 .or. replications < 1 .or. &
         size(y, 1) < 1 .or. size(y, 2) < 1 .or. &
         size(initial_mean) < 1 .or. observation_noise_dimension < 1 .or. &
         state_noise_dimension < 1 .or. approximation_iterations < 1 .or. &
         convergence_tolerance < 0.0_dp) then
         out%info = 1
         return
      end if
      allocate(initial_normals(size(initial_mean), particles, replications, &
         samples))
      allocate(state_normals(size(initial_mean), particles, size(y, 2), &
         replications, samples))
      allocate(terminal_normals(state_noise_dimension, particles, &
         replications, samples))
      allocate(resampling_uniforms(particles, size(y, 2), replications, &
         samples))
      call random_number(resampling_uniforms)
      do sample = 1, samples
         do replication = 1, replications
            call random_standard_normal_matrix( &
               initial_normals(:, :, replication, sample))
            call random_standard_normal_matrix( &
               terminal_normals(:, :, replication, sample))
            do time = 1, size(y, 2)
               call random_standard_normal_matrix( &
                  state_normals(:, :, time, replication, sample))
            end do
         end do
      end do
      out = bssm_multivariate_nonlinear_psi_post_correction_draws(y, &
         initial_mean, initial_covariance, parameters, &
         approximate_log_likelihood, accepted, scheme, &
         observation_noise_dimension, state_noise_dimension, observation, &
         transition_model, approximation_iterations, convergence_tolerance, &
         initial_normals, state_normals, terminal_normals, &
         resampling_uniforms)
   end function bssm_multivariate_nonlinear_psi_post_correction

   pure function bssm_corrected_state_moments(conditional_mean, &
      conditional_covariance, sample_weight) result(out)
      !! Combine conditional state summaries using total expectation and variance.
      real(dp), intent(in) :: conditional_mean(:, :, :) !! Conditional mean.
      real(dp), intent(in) :: conditional_covariance(:, :, :, :) !! Conditional covariance matrix.
      real(dp), intent(in) :: sample_weight(:) !! Sample weight.
      type(bssm_state_posterior_t) :: out
      real(dp), allocatable :: weight(:), difference(:)
      real(dp) :: total
      integer :: state, times, samples, sample, time

      state = size(conditional_mean, 1)
      times = size(conditional_mean, 2)
      samples = size(conditional_mean, 3)
      if (state < 1 .or. times < 1 .or. samples < 1 .or. &
         any(shape(conditional_covariance) /= &
         [state, state, times, samples]) .or. &
         size(sample_weight) /= samples .or. &
         .not. all(ieee_is_finite(conditional_mean)) .or. &
         .not. all(ieee_is_finite(conditional_covariance)) .or. &
         .not. all(ieee_is_finite(sample_weight)) .or. &
         any(sample_weight < 0.0_dp)) then
         out%info = 1
         return
      end if
      total = sum(sample_weight)
      if (.not. ieee_is_finite(total) .or. total <= 0.0_dp) then
         out%info = 1
         return
      end if
      weight = sample_weight/total
      allocate(out%mean(state, times))
      allocate(out%covariance(state, state, times))
      allocate(difference(state))
      out%mean = 0.0_dp
      do sample = 1, samples
         out%mean = out%mean + &
            weight(sample)*conditional_mean(:, :, sample)
      end do
      out%covariance = 0.0_dp
      do sample = 1, samples
         do time = 1, times
            difference = conditional_mean(:, time, sample) - out%mean(:, time)
            out%covariance(:, :, time) = out%covariance(:, :, time) + &
               weight(sample)*(conditional_covariance(:, :, time, sample) + &
               outer_product(difference, difference))
         end do
      end do
   end function bssm_corrected_state_moments

   pure function bssm_linear_signal_moments(state_posterior, loading, offset) &
      result(out)
      !! Transform corrected state moments to linear fitted or predictive signals.
      type(bssm_state_posterior_t), intent(in) :: state_posterior !! State posterior.
      real(dp), intent(in) :: loading(:, :, :) !! Loading.
      real(dp), intent(in), optional :: offset(:, :) !! Known additive offset.
      type(bssm_state_posterior_t) :: out
      integer :: series, state, times, time

      if (state_posterior%info /= 0 .or. &
         .not. allocated(state_posterior%mean) .or. &
         .not. allocated(state_posterior%covariance)) then
         out%info = 1
         return
      end if
      series = size(loading, 1)
      state = size(loading, 2)
      times = size(loading, 3)
      if (series < 1 .or. state < 1 .or. times < 1 .or. &
         any(shape(state_posterior%mean) /= [state, times]) .or. &
         any(shape(state_posterior%covariance) /= [state, state, times]) .or. &
         .not. all(ieee_is_finite(state_posterior%mean)) .or. &
         .not. all(ieee_is_finite(state_posterior%covariance)) .or. &
         .not. all(ieee_is_finite(loading))) then
         out%info = 1
         return
      end if
      if (present(offset)) then
         if (any(shape(offset) /= [series, times]) .or. &
            .not. all(ieee_is_finite(offset))) then
            out%info = 1
            return
         end if
      end if
      allocate(out%mean(series, times))
      allocate(out%covariance(series, series, times))
      do time = 1, times
         out%mean(:, time) = &
            matmul(loading(:, :, time), state_posterior%mean(:, time))
         if (present(offset)) out%mean(:, time) = out%mean(:, time) + &
            offset(:, time)
         out%covariance(:, :, time) = matmul(matmul(loading(:, :, time), &
            state_posterior%covariance(:, :, time)), &
            transpose(loading(:, :, time)))
      end do
   end function bssm_linear_signal_moments

   pure function bssm_corrected_trajectory_draws(trajectories, &
      particle_weight, sample_weight, uniforms) result(out)
      !! Resample trajectories from corrected chain and conditional weights.
      real(dp), intent(in) :: trajectories(:, :, :, :) !! Trajectories.
      real(dp), intent(in) :: particle_weight(:, :) !! Particle weight.
      real(dp), intent(in) :: sample_weight(:) !! Sample weight.
      real(dp), intent(in) :: uniforms(:, :) !! Uniforms.
      type(bssm_trajectory_sample_t) :: out
      real(dp), allocatable :: chain_probability(:), particle_probability(:)
      real(dp) :: total
      integer :: state, times, particles, samples, draws
      integer :: draw, sample, particle

      state = size(trajectories, 1)
      times = size(trajectories, 2)
      particles = size(trajectories, 3)
      samples = size(trajectories, 4)
      draws = size(uniforms, 2)
      if (state < 1 .or. times < 1 .or. particles < 1 .or. samples < 1 .or. &
         draws < 1 .or. any(shape(particle_weight) /= [particles, samples]) .or. &
         size(sample_weight) /= samples .or. size(uniforms, 1) /= 2 .or. &
         .not. all(ieee_is_finite(trajectories)) .or. &
         .not. all(ieee_is_finite(particle_weight)) .or. &
         .not. all(ieee_is_finite(sample_weight)) .or. &
         .not. all(ieee_is_finite(uniforms)) .or. &
         any(particle_weight < 0.0_dp) .or. any(sample_weight < 0.0_dp) .or. &
         any(uniforms < 0.0_dp) .or. any(uniforms > 1.0_dp)) then
         out%info = 1
         return
      end if
      total = sum(sample_weight)
      if (total <= 0.0_dp) then
         out%info = 1
         return
      end if
      chain_probability = sample_weight/total
      do sample = 1, samples
         if (sum(particle_weight(:, sample)) <= 0.0_dp) then
            out%info = 1
            return
         end if
      end do
      allocate(out%trajectories(state, times, draws))
      allocate(out%source_sample(draws), out%source_particle(draws))
      allocate(particle_probability(particles))
      do draw = 1, draws
         sample = discrete_index(chain_probability, uniforms(1, draw))
         particle_probability = particle_weight(:, sample)/ &
            sum(particle_weight(:, sample))
         particle = discrete_index(particle_probability, uniforms(2, draw))
         out%source_sample(draw) = sample
         out%source_particle(draw) = particle
         out%trajectories(:, :, draw) = &
            trajectories(:, :, particle, sample)
      end do
   end function bssm_corrected_trajectory_draws

   pure function bssm_suggest_particles(particle_counts, &
      log_likelihood_estimates, target_standard_deviation) result(out)
      !! Select the first particle count with sufficiently stable log likelihoods.
      integer, intent(in) :: particle_counts(:) !! Particle counts.
      real(dp), intent(in) :: log_likelihood_estimates(:, :) !! Log likelihood estimates.
      real(dp), intent(in), optional :: target_standard_deviation !! Target standard deviation.
      type(bssm_particle_count_t) :: out
      real(dp) :: target, center
      integer :: replications, candidates, candidate

      replications = size(log_likelihood_estimates, 1)
      candidates = size(log_likelihood_estimates, 2)
      target = 1.0_dp
      if (present(target_standard_deviation)) target = target_standard_deviation
      if (replications < 2 .or. candidates < 1 .or. &
         size(particle_counts) /= candidates .or. any(particle_counts < 1) .or. &
         .not. all(ieee_is_finite(log_likelihood_estimates)) .or. &
         .not. ieee_is_finite(target) .or. target <= 0.0_dp) then
         out%info = 1
         return
      end if
      allocate(out%standard_deviation(candidates))
      do candidate = 1, candidates
         center = sum(log_likelihood_estimates(:, candidate))/ &
            real(replications, dp)
         out%standard_deviation(candidate) = sqrt(sum( &
            (log_likelihood_estimates(:, candidate) - center)**2)/ &
            real(replications - 1, dp))
         if (out%selected_index == 0 .and. &
            out%standard_deviation(candidate) < target) then
            out%selected_index = candidate
            out%particle_count = particle_counts(candidate)
         end if
      end do
   end function bssm_suggest_particles

   pure function bssm_nonlinear_psi_suggest_particles_draws(y, initial_mean, &
      initial_covariance, parameters, noise_dimension, observation, &
      transition_model, approximation_iterations, convergence_tolerance, &
      particle_counts, initial_normals, state_normals, terminal_normals, &
      resampling_uniforms, target_standard_deviation) result(out)
      !! Select a scalar nonlinear psi particle count with supplied draws.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: parameters(:) !! Model parameter values.
      real(dp), intent(in) :: convergence_tolerance !! Convergence tolerance.
      integer, intent(in) :: noise_dimension !! Noise dimension.
      integer, intent(in) :: approximation_iterations !! Number of approximation iterations.
      integer, intent(in) :: particle_counts(:) !! Particle counts.
      real(dp), intent(in) :: initial_normals(:, :, :, :) !! Initial normals.
      real(dp), intent(in) :: state_normals(:, :, :, :, :) !! State normals.
      real(dp), intent(in) :: terminal_normals(:, :, :, :) !! Terminal normals.
      real(dp), intent(in) :: resampling_uniforms(:, :, :, :) !! Resampling uniforms.
      procedure(bssm_gaussian_observation_t) :: observation !! Observed value or vector.
      procedure(bssm_nonlinear_transition_jacobian_t) :: transition_model !! Transition model callback procedure.
      real(dp), intent(in), optional :: target_standard_deviation !! Target standard deviation.
      type(bssm_particle_count_t) :: out
      type(bssm_particle_filter_t) :: fit
      real(dp), allocatable :: log_likelihood(:, :)
      integer :: state, times, candidates, replications, maximum_particles
      integer :: candidate, replication, particles

      state = size(initial_mean)
      times = size(y)
      candidates = size(particle_counts)
      replications = size(initial_normals, 3)
      maximum_particles = size(initial_normals, 2)
      if (state < 1 .or. times < 1 .or. candidates < 1 .or. &
         replications < 2 .or. maximum_particles < 1 .or. &
         noise_dimension < 1 .or. approximation_iterations < 1 .or. &
         convergence_tolerance < 0.0_dp .or. any(particle_counts < 1) .or. &
         maxval(particle_counts) > maximum_particles .or. &
         any(shape(initial_covariance) /= [state, state]) .or. &
         any(shape(initial_normals) /= &
         [state, maximum_particles, replications, candidates]) .or. &
         any(shape(state_normals) /= &
         [state, maximum_particles, times, replications, candidates]) .or. &
         any(shape(terminal_normals) /= &
         [noise_dimension, maximum_particles, replications, candidates]) .or. &
         any(shape(resampling_uniforms) /= &
         [maximum_particles, times, replications, candidates]) .or. &
         .not. all(ieee_is_finite(initial_normals)) .or. &
         .not. all(ieee_is_finite(state_normals)) .or. &
         .not. all(ieee_is_finite(terminal_normals)) .or. &
         .not. all(ieee_is_finite(resampling_uniforms))) then
         out%info = 1
         return
      end if
      allocate(log_likelihood(replications, candidates))
      do candidate = 1, candidates
         particles = particle_counts(candidate)
         do replication = 1, replications
            fit = bssm_nonlinear_psi_filter_draws(y, initial_mean, &
               initial_covariance, parameters, noise_dimension, observation, &
               transition_model, approximation_iterations, &
               convergence_tolerance, &
               initial_normals(:, 1:particles, replication, candidate), &
               state_normals(:, 1:particles, :, replication, candidate), &
               terminal_normals(:, 1:particles, replication, candidate), &
               resampling_uniforms(1:particles, :, replication, candidate))
            if (fit%info /= 0 .or. &
               .not. ieee_is_finite(fit%log_likelihood)) then
               out%info = 100 + fit%info
               return
            end if
            log_likelihood(replication, candidate) = fit%log_likelihood
         end do
      end do
      if (present(target_standard_deviation)) then
         out = bssm_suggest_particles(particle_counts, log_likelihood, &
            target_standard_deviation)
      else
         out = bssm_suggest_particles(particle_counts, log_likelihood)
      end if
   end function bssm_nonlinear_psi_suggest_particles_draws

   function bssm_nonlinear_psi_suggest_particles(y, initial_mean, &
      initial_covariance, parameters, noise_dimension, observation, &
      transition_model, approximation_iterations, convergence_tolerance, &
      particle_counts, replications, target_standard_deviation) result(out)
      !! Select a scalar nonlinear psi particle count with shared randomness.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: parameters(:) !! Model parameter values.
      real(dp), intent(in) :: convergence_tolerance !! Convergence tolerance.
      integer, intent(in) :: noise_dimension !! Noise dimension.
      integer, intent(in) :: approximation_iterations !! Number of approximation iterations.
      integer, intent(in) :: particle_counts(:) !! Particle counts.
      integer, intent(in) :: replications !! Replications.
      procedure(bssm_gaussian_observation_t) :: observation !! Observed value or vector.
      procedure(bssm_nonlinear_transition_jacobian_t) :: transition_model !! Transition model callback procedure.
      real(dp), intent(in), optional :: target_standard_deviation !! Target standard deviation.
      type(bssm_particle_count_t) :: out
      real(dp), allocatable :: initial_normals(:, :, :, :)
      real(dp), allocatable :: state_normals(:, :, :, :, :)
      real(dp), allocatable :: terminal_normals(:, :, :, :)
      real(dp), allocatable :: resampling_uniforms(:, :, :, :)
      integer :: candidates, maximum_particles, candidate, replication, time

      candidates = size(particle_counts)
      if (candidates < 1 .or. replications < 2 .or. &
         any(particle_counts < 1) .or. size(y) < 1 .or. &
         size(initial_mean) < 1 .or. noise_dimension < 1 .or. &
         approximation_iterations < 1 .or. convergence_tolerance < 0.0_dp) then
         out%info = 1
         return
      end if
      maximum_particles = maxval(particle_counts)
      allocate(initial_normals(size(initial_mean), maximum_particles, &
         replications, candidates))
      allocate(state_normals(size(initial_mean), maximum_particles, size(y), &
         replications, candidates))
      allocate(terminal_normals(noise_dimension, maximum_particles, &
         replications, candidates))
      allocate(resampling_uniforms(maximum_particles, size(y), replications, &
         candidates))
      call random_number(resampling_uniforms)
      do candidate = 1, candidates
         do replication = 1, replications
            call random_standard_normal_matrix( &
               initial_normals(:, :, replication, candidate))
            call random_standard_normal_matrix( &
               terminal_normals(:, :, replication, candidate))
            do time = 1, size(y)
               call random_standard_normal_matrix( &
                  state_normals(:, :, time, replication, candidate))
            end do
         end do
      end do
      if (present(target_standard_deviation)) then
         out = bssm_nonlinear_psi_suggest_particles_draws(y, initial_mean, &
            initial_covariance, parameters, noise_dimension, observation, &
            transition_model, approximation_iterations, &
            convergence_tolerance, particle_counts, initial_normals, &
            state_normals, terminal_normals, resampling_uniforms, &
            target_standard_deviation)
      else
         out = bssm_nonlinear_psi_suggest_particles_draws(y, initial_mean, &
            initial_covariance, parameters, noise_dimension, observation, &
            transition_model, approximation_iterations, &
            convergence_tolerance, particle_counts, initial_normals, &
            state_normals, terminal_normals, resampling_uniforms)
      end if
   end function bssm_nonlinear_psi_suggest_particles

   pure function bssm_multivariate_nonlinear_psi_suggest_particles_draws(y, &
      initial_mean, initial_covariance, parameters, &
      observation_noise_dimension, state_noise_dimension, observation, &
      transition_model, approximation_iterations, convergence_tolerance, &
      particle_counts, initial_normals, state_normals, terminal_normals, &
      resampling_uniforms, target_standard_deviation) result(out)
      !! Select a vector nonlinear psi particle count with supplied draws.
      real(dp), intent(in) :: y(:, :) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: parameters(:) !! Model parameter values.
      real(dp), intent(in) :: convergence_tolerance !! Convergence tolerance.
      integer, intent(in) :: observation_noise_dimension !! Observation noise dimension.
      integer, intent(in) :: state_noise_dimension !! State noise dimension.
      integer, intent(in) :: approximation_iterations !! Number of approximation iterations.
      integer, intent(in) :: particle_counts(:) !! Particle counts.
      real(dp), intent(in) :: initial_normals(:, :, :, :) !! Initial normals.
      real(dp), intent(in) :: state_normals(:, :, :, :, :) !! State normals.
      real(dp), intent(in) :: terminal_normals(:, :, :, :) !! Terminal normals.
      real(dp), intent(in) :: resampling_uniforms(:, :, :, :) !! Resampling uniforms.
      procedure(bssm_multivariate_gaussian_observation_t) :: observation !! Observed value or vector.
      procedure(bssm_nonlinear_transition_jacobian_t) :: transition_model !! Transition model callback procedure.
      real(dp), intent(in), optional :: target_standard_deviation !! Target standard deviation.
      type(bssm_particle_count_t) :: out
      type(bssm_particle_filter_t) :: fit
      real(dp), allocatable :: log_likelihood(:, :)
      integer :: state, times, candidates, replications, maximum_particles
      integer :: candidate, replication, particles

      state = size(initial_mean)
      times = size(y, 2)
      candidates = size(particle_counts)
      replications = size(initial_normals, 3)
      maximum_particles = size(initial_normals, 2)
      if (size(y, 1) < 1 .or. state < 1 .or. times < 1 .or. &
         candidates < 1 .or. replications < 2 .or. maximum_particles < 1 .or. &
         observation_noise_dimension < 1 .or. state_noise_dimension < 1 .or. &
         approximation_iterations < 1 .or. convergence_tolerance < 0.0_dp .or. &
         any(particle_counts < 1) .or. &
         maxval(particle_counts) > maximum_particles .or. &
         any(shape(initial_covariance) /= [state, state]) .or. &
         any(shape(initial_normals) /= &
         [state, maximum_particles, replications, candidates]) .or. &
         any(shape(state_normals) /= &
         [state, maximum_particles, times, replications, candidates]) .or. &
         any(shape(terminal_normals) /= &
         [state_noise_dimension, maximum_particles, replications, &
         candidates]) .or. &
         any(shape(resampling_uniforms) /= &
         [maximum_particles, times, replications, candidates]) .or. &
         .not. all(ieee_is_finite(initial_normals)) .or. &
         .not. all(ieee_is_finite(state_normals)) .or. &
         .not. all(ieee_is_finite(terminal_normals)) .or. &
         .not. all(ieee_is_finite(resampling_uniforms))) then
         out%info = 1
         return
      end if
      allocate(log_likelihood(replications, candidates))
      do candidate = 1, candidates
         particles = particle_counts(candidate)
         do replication = 1, replications
            fit = bssm_multivariate_nonlinear_psi_filter_draws(y, &
               initial_mean, initial_covariance, parameters, &
               observation_noise_dimension, state_noise_dimension, &
               observation, transition_model, approximation_iterations, &
               convergence_tolerance, &
               initial_normals(:, 1:particles, replication, candidate), &
               state_normals(:, 1:particles, :, replication, candidate), &
               terminal_normals(:, 1:particles, replication, candidate), &
               resampling_uniforms(1:particles, :, replication, candidate))
            if (fit%info /= 0 .or. &
               .not. ieee_is_finite(fit%log_likelihood)) then
               out%info = 100 + fit%info
               return
            end if
            log_likelihood(replication, candidate) = fit%log_likelihood
         end do
      end do
      if (present(target_standard_deviation)) then
         out = bssm_suggest_particles(particle_counts, log_likelihood, &
            target_standard_deviation)
      else
         out = bssm_suggest_particles(particle_counts, log_likelihood)
      end if
   end function bssm_multivariate_nonlinear_psi_suggest_particles_draws

   function bssm_multivariate_nonlinear_psi_suggest_particles(y, initial_mean, &
      initial_covariance, parameters, observation_noise_dimension, &
      state_noise_dimension, observation, transition_model, &
      approximation_iterations, convergence_tolerance, particle_counts, &
      replications, target_standard_deviation) result(out)
      !! Select a vector nonlinear psi particle count with shared randomness.
      real(dp), intent(in) :: y(:, :) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: parameters(:) !! Model parameter values.
      real(dp), intent(in) :: convergence_tolerance !! Convergence tolerance.
      integer, intent(in) :: observation_noise_dimension !! Observation noise dimension.
      integer, intent(in) :: state_noise_dimension !! State noise dimension.
      integer, intent(in) :: approximation_iterations !! Number of approximation iterations.
      integer, intent(in) :: particle_counts(:) !! Particle counts.
      integer, intent(in) :: replications !! Replications.
      procedure(bssm_multivariate_gaussian_observation_t) :: observation !! Observed value or vector.
      procedure(bssm_nonlinear_transition_jacobian_t) :: transition_model !! Transition model callback procedure.
      real(dp), intent(in), optional :: target_standard_deviation !! Target standard deviation.
      type(bssm_particle_count_t) :: out
      real(dp), allocatable :: initial_normals(:, :, :, :)
      real(dp), allocatable :: state_normals(:, :, :, :, :)
      real(dp), allocatable :: terminal_normals(:, :, :, :)
      real(dp), allocatable :: resampling_uniforms(:, :, :, :)
      integer :: candidates, maximum_particles, candidate, replication, time

      candidates = size(particle_counts)
      if (candidates < 1 .or. replications < 2 .or. &
         any(particle_counts < 1) .or. size(y, 1) < 1 .or. &
         size(y, 2) < 1 .or. size(initial_mean) < 1 .or. &
         observation_noise_dimension < 1 .or. state_noise_dimension < 1 .or. &
         approximation_iterations < 1 .or. convergence_tolerance < 0.0_dp) then
         out%info = 1
         return
      end if
      maximum_particles = maxval(particle_counts)
      allocate(initial_normals(size(initial_mean), maximum_particles, &
         replications, candidates))
      allocate(state_normals(size(initial_mean), maximum_particles, &
         size(y, 2), replications, candidates))
      allocate(terminal_normals(state_noise_dimension, maximum_particles, &
         replications, candidates))
      allocate(resampling_uniforms(maximum_particles, size(y, 2), &
         replications, candidates))
      call random_number(resampling_uniforms)
      do candidate = 1, candidates
         do replication = 1, replications
            call random_standard_normal_matrix( &
               initial_normals(:, :, replication, candidate))
            call random_standard_normal_matrix( &
               terminal_normals(:, :, replication, candidate))
            do time = 1, size(y, 2)
               call random_standard_normal_matrix( &
                  state_normals(:, :, time, replication, candidate))
            end do
         end do
      end do
      if (present(target_standard_deviation)) then
         out = bssm_multivariate_nonlinear_psi_suggest_particles_draws(y, &
            initial_mean, initial_covariance, parameters, &
            observation_noise_dimension, state_noise_dimension, observation, &
            transition_model, approximation_iterations, &
            convergence_tolerance, particle_counts, initial_normals, &
            state_normals, terminal_normals, resampling_uniforms, &
            target_standard_deviation)
      else
         out = bssm_multivariate_nonlinear_psi_suggest_particles_draws(y, &
            initial_mean, initial_covariance, parameters, &
            observation_noise_dimension, state_noise_dimension, observation, &
            transition_model, approximation_iterations, &
            convergence_tolerance, particle_counts, initial_normals, &
            state_normals, terminal_normals, resampling_uniforms)
      end if
   end function bssm_multivariate_nonlinear_psi_suggest_particles

   pure real(dp) function bssm_integrated_autocorrelation_time(x, method) &
      result(iact)
      !! Estimate integrated autocorrelation time by Sokal or Geyer truncation.
      real(dp), intent(in) :: x(:) !! Input data or predictor values.
      integer, intent(in), optional :: method !! Algorithm or estimation method.
      real(dp), allocatable :: centered(:), standardized(:)
      real(dp) :: mean, scale, tau, window, correlation
      real(dp) :: gamma_zero, pair, previous_pair
      integer :: selected_method, n, lag, pair_index

      n = size(x)
      selected_method = bssm_sokal
      if (present(method)) selected_method = method
      if (n < 2 .or. .not. all(ieee_is_finite(x)) .or. &
         (selected_method /= bssm_sokal .and. &
         selected_method /= bssm_geyer)) then
         iact = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      mean = sum(x)/real(n, dp)
      centered = x - mean
      scale = sqrt(sum(centered**2)/real(n - 1, dp))
      if (scale <= tiny(1.0_dp)) then
         iact = 1.0_dp
         return
      end if
      select case (selected_method)
      case (bssm_sokal)
         standardized = centered/scale
         window = max(5.0_dp, log10(real(n, dp)))
         tau = 1.0_dp
         do lag = 1, n - 1
            correlation = dot_product(standardized(1:n-lag), &
               standardized(1+lag:n))/real(n - lag, dp)
            tau = tau + 2.0_dp*correlation
            if (real(lag, dp) > window*tau) exit
         end do
         iact = max(0.0_dp, tau)
      case (bssm_geyer)
         gamma_zero = sum(centered**2)/real(n, dp)
         tau = 0.0_dp
         previous_pair = huge(1.0_dp)
         pair_index = 0
         do
            lag = 2*pair_index
            if (lag > n - 1) exit
            pair = dot_product(centered(1:n-lag), centered(1+lag:n))/ &
               real(n, dp)
            if (lag + 1 <= n - 1) pair = pair + &
               dot_product(centered(1:n-lag-1), centered(2+lag:n))/ &
               real(n, dp)
            pair = pair/gamma_zero
            if (pair <= 0.0_dp) exit
            pair = min(pair, previous_pair)
            tau = tau + pair
            previous_pair = pair
            pair_index = pair_index + 1
         end do
         iact = max(0.0_dp, -1.0_dp + 2.0_dp*tau)
      end select
   end function bssm_integrated_autocorrelation_time

   pure real(dp) function bssm_asymptotic_variance(x, weight, method) &
      result(variance)
      !! Estimate the variance of a weighted Markov-chain mean.
      real(dp), intent(in) :: x(:) !! Input data or predictor values.
      real(dp), intent(in), optional :: weight(:) !! Weight.
      integer, intent(in), optional :: method !! Algorithm or estimation method.
      real(dp), allocatable :: work_weight(:), z(:)
      real(dp) :: weight_mean, estimate_mean, z_mean, z_variance, iact
      integer :: n, selected_method

      n = size(x)
      selected_method = bssm_sokal
      if (present(method)) selected_method = method
      if (n < 2 .or. .not. all(ieee_is_finite(x)) .or. &
         (selected_method /= bssm_sokal .and. &
         selected_method /= bssm_geyer)) then
         variance = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      allocate(work_weight(n))
      work_weight = 1.0_dp
      if (present(weight)) then
         if (size(weight) /= n .or. .not. all(ieee_is_finite(weight)) .or. &
            any(weight < 0.0_dp) .or. .not. any(weight > 0.0_dp)) then
            variance = ieee_value(0.0_dp, ieee_quiet_nan)
            return
         end if
         work_weight = weight
      end if
      weight_mean = sum(work_weight)/real(n, dp)
      estimate_mean = dot_product(work_weight, x)/sum(work_weight)
      z = work_weight*(x - estimate_mean)
      z_mean = sum(z)/real(n, dp)
      z_variance = sum((z - z_mean)**2)/real(n - 1, dp)
      iact = bssm_integrated_autocorrelation_time(z, selected_method)
      variance = z_variance*iact/(weight_mean**2*real(n, dp))
   end function bssm_asymptotic_variance

   pure function bssm_mcmc_diagnostics(samples, method, weight) result(out)
      !! Compute weighted MCMC means, errors, autocorrelation, and ESS by row.
      real(dp), intent(in) :: samples(:, :) !! Samples.
      integer, intent(in), optional :: method !! Algorithm or estimation method.
      real(dp), intent(in), optional :: weight(:) !! Weight.
      type(bssm_mcmc_diagnostics_t) :: out
      real(dp), allocatable :: work_weight(:), z(:)
      real(dp) :: total_weight, denominator
      integer :: variables, draws, variable, selected_method

      variables = size(samples, 1)
      draws = size(samples, 2)
      selected_method = bssm_sokal
      if (present(method)) selected_method = method
      if (variables < 1 .or. draws < 2 .or. &
         .not. all(ieee_is_finite(samples)) .or. &
         (selected_method /= bssm_sokal .and. &
         selected_method /= bssm_geyer)) then
         out%info = 1
         return
      end if
      allocate(work_weight(draws))
      work_weight = 1.0_dp
      if (present(weight)) then
         if (size(weight) /= draws .or. &
            .not. all(ieee_is_finite(weight)) .or. any(weight < 0.0_dp) .or. &
            .not. any(weight > 0.0_dp)) then
            out%info = 1
            return
         end if
         work_weight = weight
      end if
      total_weight = sum(work_weight)
      denominator = total_weight - dot_product(work_weight, work_weight)/ &
         total_weight
      if (denominator <= 0.0_dp) then
         out%info = 1
         return
      end if
      allocate(out%mean(variables), out%variance(variables))
      allocate(out%iact(variables), out%asymptotic_variance(variables))
      allocate(out%mcse(variables), out%ess(variables))
      allocate(z(draws))
      do variable = 1, variables
         out%mean(variable) = dot_product(work_weight, samples(variable, :))/ &
            total_weight
         out%variance(variable) = dot_product(work_weight, &
            (samples(variable, :) - out%mean(variable))**2)/denominator
         z = work_weight*(samples(variable, :) - out%mean(variable))
         out%iact(variable) = bssm_integrated_autocorrelation_time(z, &
            selected_method)
         out%asymptotic_variance(variable) = bssm_asymptotic_variance( &
            samples(variable, :), work_weight, selected_method)
         out%mcse(variable) = sqrt(max(0.0_dp, &
            out%asymptotic_variance(variable)))
         if (out%asymptotic_variance(variable) > 0.0_dp) then
            out%ess(variable) = out%variance(variable)/ &
               out%asymptotic_variance(variable)
         else if (out%variance(variable) > 0.0_dp) then
            out%ess(variable) = huge(1.0_dp)
         else
            out%ess(variable) = 0.0_dp
         end if
      end do
   end function bssm_mcmc_diagnostics

   pure function bssm_chain_diagnostics(chain, method, weight) result(out)
      !! Compute diagnostics directly from a retained BSSM parameter chain.
      type(bssm_mcmc_t), intent(in) :: chain !! Chain.
      integer, intent(in), optional :: method !! Algorithm or estimation method.
      real(dp), intent(in), optional :: weight(:) !! Weight.
      type(bssm_mcmc_diagnostics_t) :: out
      integer :: selected_method

      if (chain%info /= 0 .or. .not. allocated(chain%parameters)) then
         out%info = 1
         return
      end if
      selected_method = bssm_sokal
      if (present(method)) selected_method = method
      if (present(weight)) then
         out = bssm_mcmc_diagnostics(chain%parameters, selected_method, weight)
      else
         out = bssm_mcmc_diagnostics(chain%parameters, selected_method)
      end if
   end function bssm_chain_diagnostics

   pure function bssm_da_chain_diagnostics(chain, method, weight) result(out)
      !! Compute diagnostics from a delayed-acceptance parameter chain.
      type(bssm_da_mcmc_t), intent(in) :: chain !! Chain.
      integer, intent(in), optional :: method !! Algorithm or estimation method.
      real(dp), intent(in), optional :: weight(:) !! Weight.
      type(bssm_mcmc_diagnostics_t) :: out
      integer :: selected_method

      if (chain%info /= 0 .or. .not. allocated(chain%parameters)) then
         out%info = 1
         return
      end if
      selected_method = bssm_sokal
      if (present(method)) selected_method = method
      if (present(weight)) then
         out = bssm_mcmc_diagnostics(chain%parameters, selected_method, weight)
      else
         out = bssm_mcmc_diagnostics(chain%parameters, selected_method)
      end if
   end function bssm_da_chain_diagnostics

   pure function bssm_sde_pmmh_draws(y, initial_state, initial_parameters, &
      substeps, particles, log_density, drift, diffusion, positive, prior, &
      proposal_factor, proposal_normals, acceptance_uniforms, &
      brownian_normals, resampling_uniforms, diffusion_derivative, &
      target_acceptance, adaptation_exponent) result(out)
      !! Run a draw-driven particle marginal Metropolis-Hastings chain.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: initial_state !! Initial state vector.
      real(dp), intent(in) :: initial_parameters(:) !! Initial parameter values.
      integer, intent(in) :: substeps !! Substeps.
      integer, intent(in) :: particles !! Number of particles.
      logical, intent(in) :: positive !! Flag controlling positive.
      real(dp), intent(in) :: proposal_factor(:, :) !! Proposal factor.
      real(dp), intent(in) :: proposal_normals(:, :) !! Standard-normal proposal draws.
      real(dp), intent(in) :: acceptance_uniforms(:) !! Acceptance uniforms.
      real(dp), intent(in) :: brownian_normals(:, :, :, :) !! Brownian normals.
      real(dp), intent(in) :: resampling_uniforms(:, :, :) !! Resampling uniforms.
      procedure(bssm_nonlinear_log_density_t) :: log_density !! Log-density value.
      procedure(bssm_sde_coefficient_t) :: drift !! Drift callback procedure.
      procedure(bssm_sde_coefficient_t) :: diffusion !! Diffusion callback procedure.
      procedure(bssm_parameter_log_density_t) :: prior !! Prior-distribution specification.
      procedure(bssm_sde_coefficient_t), optional :: diffusion_derivative !! Diffusion derivative callback procedure.
      real(dp), intent(in), optional :: target_acceptance !! Target acceptance.
      real(dp), intent(in), optional :: adaptation_exponent !! Adaptation exponent.
      type(bssm_mcmc_t) :: out
      real(dp) :: exponent
      integer :: parameter_count, iterations, times

      parameter_count = size(initial_parameters)
      iterations = size(proposal_normals, 2)
      times = size(y)
      exponent = 0.6_dp
      if (present(adaptation_exponent)) exponent = adaptation_exponent
      if (parameter_count < 1 .or. iterations < 1 .or. times < 1 .or. &
         substeps < 1 .or. particles < 1 .or. &
         any(shape(brownian_normals) /= &
         [particles, substeps, times + 1, iterations + 1]) .or. &
         any(shape(resampling_uniforms) /= &
         [particles, times, iterations + 1]) .or. &
         .not. ieee_is_finite(initial_state) .or. &
         .not. all(ieee_is_finite(brownian_normals)) .or. &
         .not. all(ieee_is_finite(resampling_uniforms)) .or. &
         (present(adaptation_exponent) .and. &
         .not. present(target_acceptance))) then
         out%info = 1
         return
      end if
      if (present(target_acceptance)) then
         out = bssm_pmmh_kernel_draws(initial_parameters, prior, estimator, &
            proposal_factor, proposal_normals, acceptance_uniforms, &
            target_acceptance, exponent)
      else
         out = bssm_pmmh_kernel_draws(initial_parameters, prior, estimator, &
            proposal_factor, proposal_normals, acceptance_uniforms)
      end if

   contains

      pure function estimator(parameters, draw_index) result(estimate)
         !! Evaluate one independent SDE particle likelihood estimate.
         real(dp), intent(in) :: parameters(:) !! Model parameter values.
         integer, intent(in) :: draw_index !! Index of draw.
         type(bssm_likelihood_estimate_t) :: estimate
         type(bssm_particle_filter_t) :: fit

         if (present(diffusion_derivative)) then
            fit = bssm_sde_bootstrap_filter_draws(y, initial_state, parameters, &
               substeps, log_density, drift, diffusion, positive, &
               brownian_normals(:, :, :, draw_index), &
               resampling_uniforms(:, :, draw_index), diffusion_derivative)
         else
            fit = bssm_sde_bootstrap_filter_draws(y, initial_state, parameters, &
               substeps, log_density, drift, diffusion, positive, &
               brownian_normals(:, :, :, draw_index), &
               resampling_uniforms(:, :, draw_index))
         end if
         estimate%value = fit%log_likelihood
         estimate%info = fit%info
      end function estimator

   end function bssm_sde_pmmh_draws

   function bssm_sde_pmmh(y, initial_state, initial_parameters, substeps, &
      particles, iterations, log_density, drift, diffusion, positive, prior, &
      proposal_factor, diffusion_derivative, target_acceptance, &
      adaptation_exponent) result(out)
      !! Run SDE particle marginal Metropolis-Hastings with shared randomness.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: initial_state !! Initial state vector.
      real(dp), intent(in) :: initial_parameters(:) !! Initial parameter values.
      integer, intent(in) :: substeps !! Substeps.
      integer, intent(in) :: particles !! Number of particles.
      integer, intent(in) :: iterations !! Number of algorithm iterations.
      logical, intent(in) :: positive !! Flag controlling positive.
      real(dp), intent(in) :: proposal_factor(:, :) !! Proposal factor.
      procedure(bssm_nonlinear_log_density_t) :: log_density !! Log-density value.
      procedure(bssm_sde_coefficient_t) :: drift !! Drift callback procedure.
      procedure(bssm_sde_coefficient_t) :: diffusion !! Diffusion callback procedure.
      procedure(bssm_parameter_log_density_t) :: prior !! Prior-distribution specification.
      procedure(bssm_sde_coefficient_t), optional :: diffusion_derivative !! Diffusion derivative callback procedure.
      real(dp), intent(in), optional :: target_acceptance !! Target acceptance.
      real(dp), intent(in), optional :: adaptation_exponent !! Adaptation exponent.
      type(bssm_mcmc_t) :: out
      real(dp), allocatable :: proposal_normals(:, :), acceptance_uniforms(:)
      real(dp), allocatable :: brownian_normals(:, :, :, :)
      real(dp), allocatable :: resampling_uniforms(:, :, :)
      real(dp) :: exponent
      integer :: estimate, interval

      if (particles < 1 .or. substeps < 1 .or. iterations < 1 .or. &
         size(y) < 1 .or. size(initial_parameters) < 1 .or. &
         (present(adaptation_exponent) .and. &
         .not. present(target_acceptance))) then
         out%info = 1
         return
      end if
      allocate(proposal_normals(size(initial_parameters), iterations))
      allocate(acceptance_uniforms(iterations))
      allocate(brownian_normals(particles, substeps, size(y) + 1, &
         iterations + 1))
      allocate(resampling_uniforms(particles, size(y), iterations + 1))
      call random_standard_normal_matrix(proposal_normals)
      call random_number(acceptance_uniforms)
      call random_number(resampling_uniforms)
      do estimate = 1, iterations + 1
         do interval = 1, size(y) + 1
            call random_standard_normal_matrix( &
               brownian_normals(:, :, interval, estimate))
         end do
      end do
      exponent = 0.6_dp
      if (present(adaptation_exponent)) exponent = adaptation_exponent
      if (present(diffusion_derivative) .and. present(target_acceptance)) then
         out = bssm_sde_pmmh_draws(y, initial_state, initial_parameters, &
            substeps, particles, log_density, drift, diffusion, positive, &
            prior, proposal_factor, proposal_normals, acceptance_uniforms, &
            brownian_normals, resampling_uniforms, diffusion_derivative, &
            target_acceptance, exponent)
      else if (present(diffusion_derivative)) then
         out = bssm_sde_pmmh_draws(y, initial_state, initial_parameters, &
            substeps, particles, log_density, drift, diffusion, positive, &
            prior, proposal_factor, proposal_normals, acceptance_uniforms, &
            brownian_normals, resampling_uniforms, diffusion_derivative)
      else if (present(target_acceptance)) then
         out = bssm_sde_pmmh_draws(y, initial_state, initial_parameters, &
            substeps, particles, log_density, drift, diffusion, positive, &
            prior, proposal_factor, proposal_normals, acceptance_uniforms, &
            brownian_normals, resampling_uniforms, &
            target_acceptance=target_acceptance, &
            adaptation_exponent=exponent)
      else
         out = bssm_sde_pmmh_draws(y, initial_state, initial_parameters, &
            substeps, particles, log_density, drift, diffusion, positive, &
            prior, proposal_factor, proposal_normals, acceptance_uniforms, &
            brownian_normals, resampling_uniforms)
      end if
   end function bssm_sde_pmmh

   pure function bssm_sde_da_pmmh_draws(y, initial_state, &
      initial_parameters, coarse_substeps, fine_substeps, particles, &
      log_density, drift, diffusion, positive, prior, proposal_factor, &
      proposal_normals, first_stage_uniforms, second_stage_uniforms, &
      coarse_brownian_normals, coarse_resampling_uniforms, &
      fine_brownian_normals, fine_resampling_uniforms, &
      diffusion_derivative, target_acceptance, adaptation_exponent) result(out)
      !! Run draw-driven delayed-acceptance PMMH for an SDE model.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: initial_state !! Initial state vector.
      real(dp), intent(in) :: initial_parameters(:) !! Initial parameter values.
      integer, intent(in) :: coarse_substeps !! Coarse substeps.
      integer, intent(in) :: fine_substeps !! Fine substeps.
      integer, intent(in) :: particles !! Number of particles.
      logical, intent(in) :: positive !! Flag controlling positive.
      real(dp), intent(in) :: proposal_factor(:, :) !! Proposal factor.
      real(dp), intent(in) :: proposal_normals(:, :) !! Standard-normal proposal draws.
      real(dp), intent(in) :: first_stage_uniforms(:) !! First stage uniforms.
      real(dp), intent(in) :: second_stage_uniforms(:) !! Second stage uniforms.
      real(dp), intent(in) :: coarse_brownian_normals(:, :, :, :) !! Coarse brownian normals.
      real(dp), intent(in) :: coarse_resampling_uniforms(:, :, :) !! Coarse resampling uniforms.
      real(dp), intent(in) :: fine_brownian_normals(:, :, :, :) !! Fine brownian normals.
      real(dp), intent(in) :: fine_resampling_uniforms(:, :, :) !! Fine resampling uniforms.
      procedure(bssm_nonlinear_log_density_t) :: log_density !! Log-density value.
      procedure(bssm_sde_coefficient_t) :: drift !! Drift callback procedure.
      procedure(bssm_sde_coefficient_t) :: diffusion !! Diffusion callback procedure.
      procedure(bssm_parameter_log_density_t) :: prior !! Prior-distribution specification.
      procedure(bssm_sde_coefficient_t), optional :: diffusion_derivative !! Diffusion derivative callback procedure.
      real(dp), intent(in), optional :: target_acceptance !! Target acceptance.
      real(dp), intent(in), optional :: adaptation_exponent !! Adaptation exponent.
      type(bssm_da_mcmc_t) :: out
      real(dp) :: exponent
      integer :: parameter_count, iterations, times

      parameter_count = size(initial_parameters)
      iterations = size(proposal_normals, 2)
      times = size(y)
      exponent = 0.6_dp
      if (present(adaptation_exponent)) exponent = adaptation_exponent
      if (parameter_count < 1 .or. iterations < 1 .or. times < 1 .or. &
         coarse_substeps < 1 .or. fine_substeps < 1 .or. particles < 1 .or. &
         fine_substeps <= coarse_substeps .or. &
         any(shape(coarse_brownian_normals) /= &
         [particles, coarse_substeps, times + 1, iterations + 1]) .or. &
         any(shape(fine_brownian_normals) /= &
         [particles, fine_substeps, times + 1, iterations + 1]) .or. &
         any(shape(coarse_resampling_uniforms) /= &
         [particles, times, iterations + 1]) .or. &
         any(shape(fine_resampling_uniforms) /= &
         [particles, times, iterations + 1]) .or. &
         .not. ieee_is_finite(initial_state) .or. &
         .not. all(ieee_is_finite(coarse_brownian_normals)) .or. &
         .not. all(ieee_is_finite(fine_brownian_normals)) .or. &
         .not. all(ieee_is_finite(coarse_resampling_uniforms)) .or. &
         .not. all(ieee_is_finite(fine_resampling_uniforms)) .or. &
         (present(adaptation_exponent) .and. &
         .not. present(target_acceptance))) then
         out%info = 1
         return
      end if
      if (present(target_acceptance)) then
         out = bssm_da_pmmh_kernel_draws(initial_parameters, prior, &
            coarse_estimator, fine_estimator, proposal_factor, &
            proposal_normals, first_stage_uniforms, second_stage_uniforms, &
            target_acceptance, exponent)
      else
         out = bssm_da_pmmh_kernel_draws(initial_parameters, prior, &
            coarse_estimator, fine_estimator, proposal_factor, &
            proposal_normals, first_stage_uniforms, second_stage_uniforms)
      end if

   contains

      pure function coarse_estimator(parameters, draw_index) result(estimate)
         !! Evaluate one coarse SDE particle likelihood estimate.
         real(dp), intent(in) :: parameters(:) !! Model parameter values.
         integer, intent(in) :: draw_index !! Index of draw.
         type(bssm_likelihood_estimate_t) :: estimate
         type(bssm_particle_filter_t) :: fit

         if (present(diffusion_derivative)) then
            fit = bssm_sde_bootstrap_filter_draws(y, initial_state, parameters, &
               coarse_substeps, log_density, drift, diffusion, positive, &
               coarse_brownian_normals(:, :, :, draw_index), &
               coarse_resampling_uniforms(:, :, draw_index), &
               diffusion_derivative)
         else
            fit = bssm_sde_bootstrap_filter_draws(y, initial_state, parameters, &
               coarse_substeps, log_density, drift, diffusion, positive, &
               coarse_brownian_normals(:, :, :, draw_index), &
               coarse_resampling_uniforms(:, :, draw_index))
         end if
         estimate%value = fit%log_likelihood
         estimate%info = fit%info
      end function coarse_estimator

      pure function fine_estimator(parameters, draw_index) result(estimate)
         !! Evaluate one fine SDE particle likelihood estimate.
         real(dp), intent(in) :: parameters(:) !! Model parameter values.
         integer, intent(in) :: draw_index !! Index of draw.
         type(bssm_likelihood_estimate_t) :: estimate
         type(bssm_particle_filter_t) :: fit

         if (present(diffusion_derivative)) then
            fit = bssm_sde_bootstrap_filter_draws(y, initial_state, parameters, &
               fine_substeps, log_density, drift, diffusion, positive, &
               fine_brownian_normals(:, :, :, draw_index), &
               fine_resampling_uniforms(:, :, draw_index), &
               diffusion_derivative)
         else
            fit = bssm_sde_bootstrap_filter_draws(y, initial_state, parameters, &
               fine_substeps, log_density, drift, diffusion, positive, &
               fine_brownian_normals(:, :, :, draw_index), &
               fine_resampling_uniforms(:, :, draw_index))
         end if
         estimate%value = fit%log_likelihood
         estimate%info = fit%info
      end function fine_estimator

   end function bssm_sde_da_pmmh_draws

   function bssm_sde_da_pmmh(y, initial_state, initial_parameters, &
      coarse_substeps, fine_substeps, particles, iterations, log_density, &
      drift, diffusion, positive, prior, proposal_factor, &
      diffusion_derivative, target_acceptance, adaptation_exponent) result(out)
      !! Run delayed-acceptance SDE PMMH using the shared random stream.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: initial_state !! Initial state vector.
      real(dp), intent(in) :: initial_parameters(:) !! Initial parameter values.
      integer, intent(in) :: coarse_substeps !! Coarse substeps.
      integer, intent(in) :: fine_substeps !! Fine substeps.
      integer, intent(in) :: particles !! Number of particles.
      integer, intent(in) :: iterations !! Number of algorithm iterations.
      logical, intent(in) :: positive !! Flag controlling positive.
      real(dp), intent(in) :: proposal_factor(:, :) !! Proposal factor.
      procedure(bssm_nonlinear_log_density_t) :: log_density !! Log-density value.
      procedure(bssm_sde_coefficient_t) :: drift !! Drift callback procedure.
      procedure(bssm_sde_coefficient_t) :: diffusion !! Diffusion callback procedure.
      procedure(bssm_parameter_log_density_t) :: prior !! Prior-distribution specification.
      procedure(bssm_sde_coefficient_t), optional :: diffusion_derivative !! Diffusion derivative callback procedure.
      real(dp), intent(in), optional :: target_acceptance !! Target acceptance.
      real(dp), intent(in), optional :: adaptation_exponent !! Adaptation exponent.
      type(bssm_da_mcmc_t) :: out
      real(dp), allocatable :: proposal_normals(:, :)
      real(dp), allocatable :: first_stage_uniforms(:), second_stage_uniforms(:)
      real(dp), allocatable :: coarse_brownian(:, :, :, :)
      real(dp), allocatable :: fine_brownian(:, :, :, :)
      real(dp), allocatable :: coarse_resampling(:, :, :)
      real(dp), allocatable :: fine_resampling(:, :, :)
      real(dp) :: exponent
      integer :: estimate, interval

      if (particles < 1 .or. coarse_substeps < 1 .or. &
         fine_substeps <= coarse_substeps .or. iterations < 1 .or. &
         size(y) < 1 .or. size(initial_parameters) < 1 .or. &
         (present(adaptation_exponent) .and. &
         .not. present(target_acceptance))) then
         out%info = 1
         return
      end if
      allocate(proposal_normals(size(initial_parameters), iterations))
      allocate(first_stage_uniforms(iterations), second_stage_uniforms(iterations))
      allocate(coarse_brownian(particles, coarse_substeps, size(y) + 1, &
         iterations + 1))
      allocate(fine_brownian(particles, fine_substeps, size(y) + 1, &
         iterations + 1))
      allocate(coarse_resampling(particles, size(y), iterations + 1))
      allocate(fine_resampling(particles, size(y), iterations + 1))
      call random_standard_normal_matrix(proposal_normals)
      call random_number(first_stage_uniforms)
      call random_number(second_stage_uniforms)
      call random_number(coarse_resampling)
      call random_number(fine_resampling)
      do estimate = 1, iterations + 1
         do interval = 1, size(y) + 1
            call random_standard_normal_matrix( &
               coarse_brownian(:, :, interval, estimate))
            call random_standard_normal_matrix( &
               fine_brownian(:, :, interval, estimate))
         end do
      end do
      exponent = 0.6_dp
      if (present(adaptation_exponent)) exponent = adaptation_exponent
      if (present(diffusion_derivative) .and. present(target_acceptance)) then
         out = bssm_sde_da_pmmh_draws(y, initial_state, initial_parameters, &
            coarse_substeps, fine_substeps, particles, log_density, drift, &
            diffusion, positive, prior, proposal_factor, proposal_normals, &
            first_stage_uniforms, second_stage_uniforms, coarse_brownian, &
            coarse_resampling, fine_brownian, fine_resampling, &
            diffusion_derivative, target_acceptance, exponent)
      else if (present(diffusion_derivative)) then
         out = bssm_sde_da_pmmh_draws(y, initial_state, initial_parameters, &
            coarse_substeps, fine_substeps, particles, log_density, drift, &
            diffusion, positive, prior, proposal_factor, proposal_normals, &
            first_stage_uniforms, second_stage_uniforms, coarse_brownian, &
            coarse_resampling, fine_brownian, fine_resampling, &
            diffusion_derivative)
      else if (present(target_acceptance)) then
         out = bssm_sde_da_pmmh_draws(y, initial_state, initial_parameters, &
            coarse_substeps, fine_substeps, particles, log_density, drift, &
            diffusion, positive, prior, proposal_factor, proposal_normals, &
            first_stage_uniforms, second_stage_uniforms, coarse_brownian, &
            coarse_resampling, fine_brownian, fine_resampling, &
            target_acceptance=target_acceptance, &
            adaptation_exponent=exponent)
      else
         out = bssm_sde_da_pmmh_draws(y, initial_state, initial_parameters, &
            coarse_substeps, fine_substeps, particles, log_density, drift, &
            diffusion, positive, prior, proposal_factor, proposal_normals, &
            first_stage_uniforms, second_stage_uniforms, coarse_brownian, &
            coarse_resampling, fine_brownian, fine_resampling)
      end if
   end function bssm_sde_da_pmmh

   pure function bssm_bootstrap_filter_draws(y, observation_loading, &
      transition, state_noise_loading, initial_mean, initial_covariance, &
      distribution, phi, initial_normals, innovation_normals, &
      resampling_uniforms, offset, auxiliary, state_offset) result(out)
      !! Run a draw-driven bootstrap filter with linear Gaussian dynamics.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: observation_loading(:) !! Observation loading matrix.
      real(dp), intent(in) :: transition(:, :) !! State transition matrix.
      real(dp), intent(in) :: state_noise_loading(:, :) !! State noise loading.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      integer, intent(in) :: distribution !! Probability-distribution specification.
      real(dp), intent(in) :: phi !! Autoregressive or model coefficient.
      real(dp), intent(in) :: initial_normals(:, :) !! Initial normals.
      real(dp), intent(in) :: innovation_normals(:, :, :) !! Innovation normals.
      real(dp), intent(in) :: resampling_uniforms(:, :) !! Resampling uniforms.
      real(dp), intent(in), optional :: offset(:) !! Known additive offset.
      real(dp), intent(in), optional :: auxiliary(:) !! Auxiliary.
      real(dp), intent(in), optional :: state_offset(:) !! State offset.
      type(bssm_particle_filter_t) :: out
      real(dp), allocatable :: observation_schedule(:, :)
      real(dp), allocatable :: transition_schedule(:, :, :)
      real(dp), allocatable :: noise_schedule(:, :, :)
      real(dp), allocatable :: offset_work(:), auxiliary_work(:)
      real(dp), allocatable :: state_offset_schedule(:, :)

      observation_schedule = spread(observation_loading, 2, size(y))
      transition_schedule = spread(transition, 3, size(y))
      noise_schedule = spread(state_noise_loading, 3, size(y))
      allocate(offset_work(size(y)), auxiliary_work(size(y)))
      allocate(state_offset_schedule(size(initial_mean), size(y)))
      offset_work = 0.0_dp
      auxiliary_work = 1.0_dp
      state_offset_schedule = 0.0_dp
      if (present(offset)) then
         if (size(offset) /= size(y)) then
            out%info = 1
            return
         end if
         offset_work = offset
      end if
      if (present(auxiliary)) then
         if (size(auxiliary) /= size(y)) then
            out%info = 1
            return
         end if
         auxiliary_work = auxiliary
      end if
      if (present(state_offset)) then
         if (size(state_offset) /= size(initial_mean)) then
            out%info = 1
            return
         end if
         state_offset_schedule = spread(state_offset, 2, size(y))
      end if
      out = bssm_bootstrap_filter_tv_draws(y, observation_schedule, &
         transition_schedule, noise_schedule, initial_mean, &
         initial_covariance, distribution, phi, initial_normals, &
         innovation_normals, resampling_uniforms, offset_work, auxiliary_work, &
         state_offset_schedule)
   end function bssm_bootstrap_filter_draws

   pure function bssm_bootstrap_filter_tv_draws(y, observation_loading, &
      transition, state_noise_loading, initial_mean, initial_covariance, &
      distribution, phi, initial_normals, innovation_normals, &
      resampling_uniforms, offset, auxiliary, state_offset) result(out)
      !! Run a draw-driven bootstrap filter with scheduled linear dynamics.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: observation_loading(:, :) !! Observation loading matrix.
      real(dp), intent(in) :: transition(:, :, :) !! State transition matrix.
      real(dp), intent(in) :: state_noise_loading(:, :, :) !! State noise loading.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      integer, intent(in) :: distribution !! Probability-distribution specification.
      real(dp), intent(in) :: phi !! Autoregressive or model coefficient.
      real(dp), intent(in) :: initial_normals(:, :) !! Initial normals.
      real(dp), intent(in) :: innovation_normals(:, :, :) !! Innovation normals.
      real(dp), intent(in) :: resampling_uniforms(:, :) !! Resampling uniforms.
      real(dp), intent(in), optional :: offset(:) !! Known additive offset.
      real(dp), intent(in), optional :: auxiliary(:) !! Auxiliary.
      real(dp), intent(in), optional :: state_offset(:, :) !! State offset.
      type(bssm_particle_filter_t) :: out
      real(dp), allocatable :: factor(:, :), log_weight(:), probability(:)
      real(dp), allocatable :: offset_work(:), auxiliary_work(:)
      real(dp), allocatable :: state_offset_work(:, :)
      real(dp), allocatable :: parent(:, :)
      real(dp) :: increment, signal
      integer :: state, particles, times, noise, particle, time, info

      state = size(initial_mean)
      particles = size(initial_normals, 2)
      times = size(y)
      noise = size(state_noise_loading, 2)
      if (.not. valid_filter_dimensions()) then
         out%info = 1
         return
      end if
      allocate(offset_work(times), auxiliary_work(times))
      allocate(state_offset_work(state, times))
      offset_work = 0.0_dp
      auxiliary_work = 1.0_dp
      state_offset_work = 0.0_dp
      if (present(offset)) then
         if (size(offset) /= times) then
            out%info = 1
            return
         end if
         offset_work = offset
      end if
      if (present(auxiliary)) then
         if (size(auxiliary) /= times) then
            out%info = 1
            return
         end if
         auxiliary_work = auxiliary
      end if
      if (present(state_offset)) then
         if (any(shape(state_offset) /= [state, times])) then
            out%info = 1
            return
         end if
         state_offset_work = state_offset
      end if
      call positive_semidefinite_factor(initial_covariance, factor, info)
      if (info /= 0) then
         out%info = 2
         return
      end if
      allocate(out%particles(state, particles, times + 1))
      allocate(out%weights(particles, times), out%ancestors(particles, times))
      allocate(out%predicted_mean(state, times + 1))
      allocate(out%filtered_mean(state, times))
      allocate(out%predicted_covariance(state, state, times + 1))
      allocate(out%filtered_covariance(state, state, times))
      allocate(log_weight(particles), probability(particles))
      allocate(parent(state, particles))
      do particle = 1, particles
         out%particles(:, particle, 1) = initial_mean + &
            matmul(factor, initial_normals(:, particle))
      end do
      do time = 1, times
         call particle_summary(out%particles(:, :, time), &
            out%predicted_mean(:, time), &
            out%predicted_covariance(:, :, time))
         if (ieee_is_finite(y(time))) then
            do particle = 1, particles
               signal = offset_work(time) + dot_product(observation_loading(:, time), &
                  out%particles(:, particle, time))
               log_weight(particle) = bssm_observation_log_density(y(time), &
                  signal, distribution, phi, auxiliary_work(time))
            end do
            call normalize_log_weights(log_weight, probability, increment, info)
            if (info /= 0) then
               out%info = 3
               out%log_likelihood = -huge(1.0_dp)
               return
            end if
            out%log_likelihood = out%log_likelihood + increment
         else
            probability = 1.0_dp/real(particles, dp)
         end if
         out%weights(:, time) = probability
         call weighted_particle_summary(out%particles(:, :, time), probability, &
            out%filtered_mean(:, time), &
            out%filtered_covariance(:, :, time))
         out%ancestors(:, time) = bssm_stratified_resample(probability, &
            resampling_uniforms(:, time))
         do particle = 1, particles
            parent(:, particle) = &
               out%particles(:, out%ancestors(particle, time), time)
            out%particles(:, particle, time + 1) = &
               state_offset_work(:, time) + &
               matmul(transition(:, :, time), parent(:, particle)) + &
               matmul(state_noise_loading(:, :, time), &
               innovation_normals(:, particle, time))
         end do
      end do
      call particle_summary(out%particles(:, :, times + 1), &
         out%predicted_mean(:, times + 1), &
         out%predicted_covariance(:, :, times + 1))

   contains

      pure logical function valid_filter_dimensions() result(valid)
         !! Validate all state, particle, and random-array dimensions.
         valid = state > 0 .and. particles > 0 .and. times > 0 .and. noise > 0
         valid = valid .and. all(shape(observation_loading) == [state, times])
         valid = valid .and. all(shape(transition) == [state, state, times])
         valid = valid .and. all(shape(state_noise_loading) == &
            [state, noise, times])
         valid = valid .and. all(shape(initial_covariance) == [state, state])
         valid = valid .and. size(initial_normals, 1) == state
         valid = valid .and. all(shape(innovation_normals) == &
            [noise, particles, times])
         valid = valid .and. all(shape(resampling_uniforms) == [particles, times])
         valid = valid .and. distribution >= bssm_svm .and. &
            distribution <= bssm_gaussian
         valid = valid .and. all(ieee_is_finite(initial_mean))
         valid = valid .and. all(ieee_is_finite(initial_covariance))
         valid = valid .and. all(ieee_is_finite(initial_normals))
         valid = valid .and. all(ieee_is_finite(innovation_normals))
         valid = valid .and. all(ieee_is_finite(resampling_uniforms))
      end function valid_filter_dimensions

   end function bssm_bootstrap_filter_tv_draws

   function bssm_bootstrap_filter(y, observation_loading, transition, &
      state_noise_loading, initial_mean, initial_covariance, distribution, &
      phi, particles, offset, auxiliary, state_offset) result(out)
      !! Run a bootstrap filter using the shared pseudo-random stream.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: observation_loading(:) !! Observation loading matrix.
      real(dp), intent(in) :: transition(:, :) !! State transition matrix.
      real(dp), intent(in) :: state_noise_loading(:, :) !! State noise loading.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      integer, intent(in) :: distribution !! Probability-distribution specification.
      integer, intent(in) :: particles !! Number of particles.
      real(dp), intent(in) :: phi !! Autoregressive or model coefficient.
      real(dp), intent(in), optional :: offset(:) !! Known additive offset.
      real(dp), intent(in), optional :: auxiliary(:) !! Auxiliary.
      real(dp), intent(in), optional :: state_offset(:) !! State offset.
      type(bssm_particle_filter_t) :: out
      real(dp), allocatable :: initial_normals(:, :)
      real(dp), allocatable :: innovation_normals(:, :, :), uniforms(:, :)
      real(dp), allocatable :: offset_work(:), auxiliary_work(:)
      real(dp), allocatable :: state_offset_work(:)
      integer :: time

      if (particles < 1 .or. size(y) < 1) then
         out%info = 1
         return
      end if
      allocate(initial_normals(size(initial_mean), particles))
      allocate(innovation_normals(size(state_noise_loading, 2), particles, &
         size(y)), uniforms(particles, size(y)))
      allocate(offset_work(size(y)), auxiliary_work(size(y)))
      allocate(state_offset_work(size(initial_mean)))
      offset_work = 0.0_dp
      auxiliary_work = 1.0_dp
      state_offset_work = 0.0_dp
      if (present(offset)) then
         if (size(offset) /= size(y)) then
            out%info = 1
            return
         end if
         offset_work = offset
      end if
      if (present(auxiliary)) then
         if (size(auxiliary) /= size(y)) then
            out%info = 1
            return
         end if
         auxiliary_work = auxiliary
      end if
      if (present(state_offset)) then
         if (size(state_offset) /= size(initial_mean)) then
            out%info = 1
            return
         end if
         state_offset_work = state_offset
      end if
      call random_standard_normal_matrix(initial_normals)
      do time = 1, size(y)
         call random_standard_normal_matrix(innovation_normals(:, :, time))
      end do
      call random_number(uniforms)
      out = bssm_bootstrap_filter_draws(y, observation_loading, transition, &
         state_noise_loading, initial_mean, initial_covariance, distribution, &
         phi, initial_normals, innovation_normals, uniforms, offset_work, &
         auxiliary_work, state_offset_work)
   end function bssm_bootstrap_filter

   function bssm_bootstrap_filter_tv(y, observation_loading, transition, &
      state_noise_loading, initial_mean, initial_covariance, distribution, &
      phi, particles, offset, auxiliary, state_offset) result(out)
      !! Run a scheduled linear bootstrap filter using the shared random stream.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: observation_loading(:, :) !! Observation loading matrix.
      real(dp), intent(in) :: transition(:, :, :) !! State transition matrix.
      real(dp), intent(in) :: state_noise_loading(:, :, :) !! State noise loading.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      integer, intent(in) :: distribution !! Probability-distribution specification.
      integer, intent(in) :: particles !! Number of particles.
      real(dp), intent(in) :: phi !! Autoregressive or model coefficient.
      real(dp), intent(in), optional :: offset(:) !! Known additive offset.
      real(dp), intent(in), optional :: auxiliary(:) !! Auxiliary.
      real(dp), intent(in), optional :: state_offset(:, :) !! State offset.
      type(bssm_particle_filter_t) :: out
      real(dp), allocatable :: initial_normals(:, :)
      real(dp), allocatable :: innovation_normals(:, :, :), uniforms(:, :)
      real(dp), allocatable :: offset_work(:), auxiliary_work(:)
      real(dp), allocatable :: state_offset_work(:, :)
      integer :: time

      if (particles < 1 .or. size(y) < 1) then
         out%info = 1
         return
      end if
      allocate(initial_normals(size(initial_mean), particles))
      allocate(innovation_normals(size(state_noise_loading, 2), particles, &
         size(y)), uniforms(particles, size(y)))
      allocate(offset_work(size(y)), auxiliary_work(size(y)))
      allocate(state_offset_work(size(initial_mean), size(y)))
      offset_work = 0.0_dp
      auxiliary_work = 1.0_dp
      state_offset_work = 0.0_dp
      if (present(offset)) then
         if (size(offset) /= size(y)) then
            out%info = 1
            return
         end if
         offset_work = offset
      end if
      if (present(auxiliary)) then
         if (size(auxiliary) /= size(y)) then
            out%info = 1
            return
         end if
         auxiliary_work = auxiliary
      end if
      if (present(state_offset)) then
         if (any(shape(state_offset) /= [size(initial_mean), size(y)])) then
            out%info = 1
            return
         end if
         state_offset_work = state_offset
      end if
      call random_standard_normal_matrix(initial_normals)
      do time = 1, size(y)
         call random_standard_normal_matrix(innovation_normals(:, :, time))
      end do
      call random_number(uniforms)
      out = bssm_bootstrap_filter_tv_draws(y, observation_loading, transition, &
         state_noise_loading, initial_mean, initial_covariance, distribution, &
         phi, initial_normals, innovation_normals, uniforms, offset_work, &
         auxiliary_work, state_offset_work)
   end function bssm_bootstrap_filter_tv

   pure function bssm_nonlinear_bootstrap_filter_draws(y, initial_mean, &
      initial_covariance, parameters, noise_dimension, log_density, transition, &
      initial_normals, innovation_normals, resampling_uniforms) result(out)
      !! Run a draw-driven bootstrap filter for a nonlinear Gaussian state model.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: parameters(:) !! Model parameter values.
      real(dp), intent(in) :: initial_normals(:, :) !! Initial normals.
      integer, intent(in) :: noise_dimension !! Noise dimension.
      procedure(bssm_nonlinear_log_density_t) :: log_density !! Log-density value.
      procedure(bssm_nonlinear_transition_t) :: transition !! State transition matrix.
      real(dp), intent(in) :: innovation_normals(:, :, :) !! Innovation normals.
      real(dp), intent(in) :: resampling_uniforms(:, :) !! Resampling uniforms.
      type(bssm_particle_filter_t) :: out
      real(dp), allocatable :: factor(:, :), log_weight(:), probability(:)
      real(dp), allocatable :: parent(:, :), transition_mean(:)
      real(dp), allocatable :: noise_loading(:, :)
      real(dp) :: increment
      integer :: state, particles, times, particle, time, info

      state = size(initial_mean)
      particles = size(initial_normals, 2)
      times = size(y)
      if (state < 1 .or. particles < 1 .or. times < 1 .or. &
         noise_dimension < 1 .or. size(initial_normals, 1) /= state .or. &
         any(shape(initial_covariance) /= [state, state]) .or. &
         any(shape(innovation_normals) /= &
         [noise_dimension, particles, times]) .or. &
         any(shape(resampling_uniforms) /= [particles, times]) .or. &
         .not. all(ieee_is_finite(initial_mean)) .or. &
         .not. all(ieee_is_finite(initial_covariance)) .or. &
         .not. all(ieee_is_finite(parameters)) .or. &
         .not. all(ieee_is_finite(initial_normals)) .or. &
         .not. all(ieee_is_finite(innovation_normals)) .or. &
         .not. all(ieee_is_finite(resampling_uniforms))) then
         out%info = 1
         return
      end if
      call positive_semidefinite_factor(initial_covariance, factor, info)
      if (info /= 0) then
         out%info = 2
         return
      end if
      allocate(out%particles(state, particles, times + 1))
      allocate(out%weights(particles, times), out%ancestors(particles, times))
      allocate(out%predicted_mean(state, times + 1))
      allocate(out%filtered_mean(state, times))
      allocate(out%predicted_covariance(state, state, times + 1))
      allocate(out%filtered_covariance(state, state, times))
      allocate(log_weight(particles), probability(particles))
      allocate(parent(state, particles), transition_mean(state))
      allocate(noise_loading(state, noise_dimension))
      do particle = 1, particles
         out%particles(:, particle, 1) = initial_mean + &
            matmul(factor, initial_normals(:, particle))
      end do
      do time = 1, times
         call particle_summary(out%particles(:, :, time), &
            out%predicted_mean(:, time), &
            out%predicted_covariance(:, :, time))
         if (ieee_is_finite(y(time))) then
            do particle = 1, particles
               log_weight(particle) = log_density(time, y(time), &
                  out%particles(:, particle, time), parameters)
            end do
            call normalize_log_weights(log_weight, probability, increment, info)
            if (info /= 0) then
               out%info = 3
               out%log_likelihood = -huge(1.0_dp)
               return
            end if
            out%log_likelihood = out%log_likelihood + increment
         else
            probability = 1.0_dp/real(particles, dp)
         end if
         out%weights(:, time) = probability
         call weighted_particle_summary(out%particles(:, :, time), probability, &
            out%filtered_mean(:, time), &
            out%filtered_covariance(:, :, time))
         out%ancestors(:, time) = bssm_stratified_resample(probability, &
            resampling_uniforms(:, time))
         do particle = 1, particles
            parent(:, particle) = &
               out%particles(:, out%ancestors(particle, time), time)
            call transition(time, parent(:, particle), parameters, &
               transition_mean, noise_loading)
            if (.not. all(ieee_is_finite(transition_mean)) .or. &
               .not. all(ieee_is_finite(noise_loading))) then
               out%info = 4
               return
            end if
            out%particles(:, particle, time + 1) = transition_mean + &
               matmul(noise_loading, innovation_normals(:, particle, time))
         end do
      end do
      call particle_summary(out%particles(:, :, times + 1), &
         out%predicted_mean(:, times + 1), &
         out%predicted_covariance(:, :, times + 1))
   end function bssm_nonlinear_bootstrap_filter_draws

   function bssm_nonlinear_bootstrap_filter(y, initial_mean, &
      initial_covariance, parameters, noise_dimension, log_density, transition, &
      particles) result(out)
      !! Run a nonlinear bootstrap filter using the shared random stream.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: parameters(:) !! Model parameter values.
      integer, intent(in) :: noise_dimension !! Noise dimension.
      integer, intent(in) :: particles !! Number of particles.
      procedure(bssm_nonlinear_log_density_t) :: log_density !! Log-density value.
      procedure(bssm_nonlinear_transition_t) :: transition !! State transition matrix.
      type(bssm_particle_filter_t) :: out
      real(dp), allocatable :: initial_normals(:, :)
      real(dp), allocatable :: innovation_normals(:, :, :), uniforms(:, :)
      integer :: time

      if (particles < 1 .or. noise_dimension < 1 .or. size(y) < 1) then
         out%info = 1
         return
      end if
      allocate(initial_normals(size(initial_mean), particles))
      allocate(innovation_normals(noise_dimension, particles, size(y)))
      allocate(uniforms(particles, size(y)))
      call random_standard_normal_matrix(initial_normals)
      do time = 1, size(y)
         call random_standard_normal_matrix(innovation_normals(:, :, time))
      end do
      call random_number(uniforms)
      out = bssm_nonlinear_bootstrap_filter_draws(y, initial_mean, &
         initial_covariance, parameters, noise_dimension, log_density, &
         transition, initial_normals, innovation_normals, uniforms)
   end function bssm_nonlinear_bootstrap_filter

   pure function bssm_multivariate_nonlinear_bootstrap_filter_draws(y, &
      initial_mean, initial_covariance, parameters, &
      observation_noise_dimension, state_noise_dimension, observation, &
      transition, initial_normals, innovation_normals, &
      resampling_uniforms) result(out)
      !! Run a supplied-draw multivariate nonlinear Gaussian bootstrap filter.
      real(dp), intent(in) :: y(:, :) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: parameters(:) !! Model parameter values.
      integer, intent(in) :: observation_noise_dimension !! Observation noise dimension.
      integer, intent(in) :: state_noise_dimension !! State noise dimension.
      procedure(bssm_multivariate_gaussian_observation_t) :: observation !! Observed value or vector.
      procedure(bssm_nonlinear_transition_t) :: transition !! State transition matrix.
      real(dp), intent(in) :: initial_normals(:, :) !! Initial normals.
      real(dp), intent(in) :: innovation_normals(:, :, :) !! Innovation normals.
      real(dp), intent(in) :: resampling_uniforms(:, :) !! Resampling uniforms.
      type(bssm_particle_filter_t) :: out
      real(dp), allocatable :: factor(:, :), log_weight(:), probability(:)
      real(dp), allocatable :: parent(:, :), transition_mean(:)
      real(dp), allocatable :: noise_loading(:, :)
      real(dp) :: increment
      integer :: state, series, particles, times, particle, time, info

      series = size(y, 1)
      times = size(y, 2)
      state = size(initial_mean)
      particles = size(initial_normals, 2)
      if (series < 1 .or. state < 1 .or. particles < 1 .or. times < 1 .or. &
         observation_noise_dimension < 1 .or. state_noise_dimension < 1 .or. &
         size(initial_normals, 1) /= state .or. &
         any(shape(initial_covariance) /= [state, state]) .or. &
         any(shape(innovation_normals) /= &
         [state_noise_dimension, particles, times]) .or. &
         any(shape(resampling_uniforms) /= [particles, times]) .or. &
         .not. all(ieee_is_finite(initial_mean)) .or. &
         .not. all(ieee_is_finite(initial_covariance)) .or. &
         .not. all(ieee_is_finite(parameters)) .or. &
         .not. all(ieee_is_finite(initial_normals)) .or. &
         .not. all(ieee_is_finite(innovation_normals)) .or. &
         .not. all(ieee_is_finite(resampling_uniforms))) then
         out%info = 1
         return
      end if
      call positive_semidefinite_factor(initial_covariance, factor, info)
      if (info /= 0) then
         out%info = 2
         return
      end if
      allocate(out%particles(state, particles, times + 1))
      allocate(out%weights(particles, times), out%ancestors(particles, times))
      allocate(out%predicted_mean(state, times + 1))
      allocate(out%filtered_mean(state, times))
      allocate(out%predicted_covariance(state, state, times + 1))
      allocate(out%filtered_covariance(state, state, times))
      allocate(log_weight(particles), probability(particles))
      allocate(parent(state, particles), transition_mean(state))
      allocate(noise_loading(state, state_noise_dimension))
      do particle = 1, particles
         out%particles(:, particle, 1) = initial_mean + &
            matmul(factor, initial_normals(:, particle))
      end do
      do time = 1, times
         call particle_summary(out%particles(:, :, time), &
            out%predicted_mean(:, time), &
            out%predicted_covariance(:, :, time))
         if (any(ieee_is_finite(y(:, time)))) then
            do particle = 1, particles
               log_weight(particle) = &
                  multivariate_gaussian_observation_log_density(time, &
                  y(:, time), out%particles(:, particle, time), parameters, &
                  observation_noise_dimension, observation)
            end do
            call normalize_log_weights(log_weight, probability, increment, &
               info)
            if (info /= 0) then
               out%info = 3
               out%log_likelihood = -huge(1.0_dp)
               return
            end if
            out%log_likelihood = out%log_likelihood + increment
         else
            probability = 1.0_dp/real(particles, dp)
         end if
         out%weights(:, time) = probability
         call weighted_particle_summary(out%particles(:, :, time), &
            probability, out%filtered_mean(:, time), &
            out%filtered_covariance(:, :, time))
         out%ancestors(:, time) = bssm_stratified_resample(probability, &
            resampling_uniforms(:, time))
         do particle = 1, particles
            parent(:, particle) = &
               out%particles(:, out%ancestors(particle, time), time)
            call transition(time, parent(:, particle), parameters, &
               transition_mean, noise_loading)
            if (.not. all(ieee_is_finite(transition_mean)) .or. &
               .not. all(ieee_is_finite(noise_loading))) then
               out%info = 4
               return
            end if
            out%particles(:, particle, time + 1) = transition_mean + &
               matmul(noise_loading, innovation_normals(:, particle, time))
         end do
      end do
      call particle_summary(out%particles(:, :, times + 1), &
         out%predicted_mean(:, times + 1), &
         out%predicted_covariance(:, :, times + 1))
   end function bssm_multivariate_nonlinear_bootstrap_filter_draws

   function bssm_multivariate_nonlinear_bootstrap_filter(y, initial_mean, &
      initial_covariance, parameters, observation_noise_dimension, &
      state_noise_dimension, observation, transition, particles) result(out)
      !! Run a multivariate nonlinear bootstrap filter from the random stream.
      real(dp), intent(in) :: y(:, :) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: parameters(:) !! Model parameter values.
      integer, intent(in) :: observation_noise_dimension !! Observation noise dimension.
      integer, intent(in) :: state_noise_dimension !! State noise dimension.
      integer, intent(in) :: particles !! Number of particles.
      procedure(bssm_multivariate_gaussian_observation_t) :: observation !! Observed value or vector.
      procedure(bssm_nonlinear_transition_t) :: transition !! State transition matrix.
      type(bssm_particle_filter_t) :: out
      real(dp), allocatable :: initial_normals(:, :)
      real(dp), allocatable :: innovation_normals(:, :, :), uniforms(:, :)
      integer :: time

      if (particles < 1 .or. observation_noise_dimension < 1 .or. &
         state_noise_dimension < 1 .or. size(y, 1) < 1 .or. &
         size(y, 2) < 1) then
         out%info = 1
         return
      end if
      allocate(initial_normals(size(initial_mean), particles))
      allocate(innovation_normals(state_noise_dimension, particles, &
         size(y, 2)))
      allocate(uniforms(particles, size(y, 2)))
      call random_standard_normal_matrix(initial_normals)
      do time = 1, size(y, 2)
         call random_standard_normal_matrix(innovation_normals(:, :, time))
      end do
      call random_number(uniforms)
      out = bssm_multivariate_nonlinear_bootstrap_filter_draws(y, &
         initial_mean, initial_covariance, parameters, &
         observation_noise_dimension, state_noise_dimension, observation, &
         transition, initial_normals, innovation_normals, uniforms)
   end function bssm_multivariate_nonlinear_bootstrap_filter

   pure function bssm_nonlinear_pmmh_draws(y, initial_mean, &
      initial_covariance, initial_parameters, noise_dimension, log_density, &
      transition, prior, proposal_factor, proposal_normals, &
      acceptance_uniforms, initial_normals, innovation_normals, &
      resampling_uniforms, target_acceptance, adaptation_exponent) result(out)
      !! Run draw-driven PMMH with a nonlinear bootstrap particle filter.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: initial_parameters(:) !! Initial parameter values.
      real(dp), intent(in) :: proposal_factor(:, :) !! Proposal factor.
      real(dp), intent(in) :: proposal_normals(:, :) !! Standard-normal proposal draws.
      real(dp), intent(in) :: acceptance_uniforms(:) !! Acceptance uniforms.
      real(dp), intent(in) :: initial_normals(:, :, :) !! Initial normals.
      real(dp), intent(in) :: innovation_normals(:, :, :, :) !! Innovation normals.
      real(dp), intent(in) :: resampling_uniforms(:, :, :) !! Resampling uniforms.
      integer, intent(in) :: noise_dimension !! Noise dimension.
      procedure(bssm_nonlinear_log_density_t) :: log_density !! Log-density value.
      procedure(bssm_nonlinear_transition_t) :: transition !! State transition matrix.
      procedure(bssm_parameter_log_density_t) :: prior !! Prior-distribution specification.
      real(dp), intent(in), optional :: target_acceptance !! Target acceptance.
      real(dp), intent(in), optional :: adaptation_exponent !! Adaptation exponent.
      type(bssm_mcmc_t) :: out
      integer :: state, particles, times, iterations
      real(dp) :: exponent

      state = size(initial_mean)
      particles = size(initial_normals, 2)
      times = size(y)
      iterations = size(proposal_normals, 2)
      if (state < 1 .or. particles < 1 .or. times < 1 .or. &
         noise_dimension < 1 .or. iterations < 1 .or. &
         any(shape(initial_covariance) /= [state, state]) .or. &
         any(shape(initial_normals) /= &
         [state, particles, iterations + 1]) .or. &
         any(shape(innovation_normals) /= &
         [noise_dimension, particles, times, iterations + 1]) .or. &
         any(shape(resampling_uniforms) /= &
         [particles, times, iterations + 1]) .or. &
         .not. all(ieee_is_finite(initial_normals)) .or. &
         .not. all(ieee_is_finite(innovation_normals)) .or. &
         .not. all(ieee_is_finite(resampling_uniforms)) .or. &
         (present(adaptation_exponent) .and. &
         .not. present(target_acceptance))) then
         out%info = 1
         return
      end if
      exponent = 0.6_dp
      if (present(adaptation_exponent)) exponent = adaptation_exponent
      if (present(target_acceptance)) then
         out = bssm_pmmh_kernel_draws(initial_parameters, prior, estimator, &
            proposal_factor, proposal_normals, acceptance_uniforms, &
            target_acceptance, exponent)
      else
         out = bssm_pmmh_kernel_draws(initial_parameters, prior, estimator, &
            proposal_factor, proposal_normals, acceptance_uniforms)
      end if

   contains

      pure function estimator(parameters, draw_index) result(estimate)
         !! Estimate the nonlinear bootstrap-filter log likelihood.
         real(dp), intent(in) :: parameters(:) !! Model parameter values.
         integer, intent(in) :: draw_index !! Index of draw.
         type(bssm_likelihood_estimate_t) :: estimate
         type(bssm_particle_filter_t) :: fit

         fit = bssm_nonlinear_bootstrap_filter_draws(y, initial_mean, &
            initial_covariance, parameters, noise_dimension, log_density, &
            transition, initial_normals(:, :, draw_index), &
            innovation_normals(:, :, :, draw_index), &
            resampling_uniforms(:, :, draw_index))
         estimate%value = fit%log_likelihood
         estimate%info = fit%info
      end function estimator

   end function bssm_nonlinear_pmmh_draws

   function bssm_nonlinear_pmmh(y, initial_mean, initial_covariance, &
      initial_parameters, noise_dimension, log_density, transition, prior, &
      proposal_factor, particles, iterations, target_acceptance, &
      adaptation_exponent) result(out)
      !! Run nonlinear bootstrap-filter PMMH using the shared random stream.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: initial_parameters(:) !! Initial parameter values.
      real(dp), intent(in) :: proposal_factor(:, :) !! Proposal factor.
      integer, intent(in) :: noise_dimension !! Noise dimension.
      integer, intent(in) :: particles !! Number of particles.
      integer, intent(in) :: iterations !! Number of algorithm iterations.
      procedure(bssm_nonlinear_log_density_t) :: log_density !! Log-density value.
      procedure(bssm_nonlinear_transition_t) :: transition !! State transition matrix.
      procedure(bssm_parameter_log_density_t) :: prior !! Prior-distribution specification.
      real(dp), intent(in), optional :: target_acceptance !! Target acceptance.
      real(dp), intent(in), optional :: adaptation_exponent !! Adaptation exponent.
      type(bssm_mcmc_t) :: out
      real(dp), allocatable :: proposal_normals(:, :), acceptance_uniforms(:)
      real(dp), allocatable :: initial_normals(:, :, :)
      real(dp), allocatable :: innovation_normals(:, :, :, :)
      real(dp), allocatable :: resampling_uniforms(:, :, :)
      real(dp) :: exponent
      integer :: draw_index, time

      if (particles < 1 .or. iterations < 1 .or. noise_dimension < 1 .or. &
         size(y) < 1 .or. size(initial_mean) < 1 .or. &
         size(initial_parameters) < 1 .or. &
         (present(adaptation_exponent) .and. &
         .not. present(target_acceptance))) then
         out%info = 1
         return
      end if
      allocate(proposal_normals(size(initial_parameters), iterations))
      allocate(acceptance_uniforms(iterations))
      allocate(initial_normals(size(initial_mean), particles, iterations + 1))
      allocate(innovation_normals(noise_dimension, particles, size(y), &
         iterations + 1))
      allocate(resampling_uniforms(particles, size(y), iterations + 1))
      call random_standard_normal_matrix(proposal_normals)
      call random_number(acceptance_uniforms)
      call random_number(resampling_uniforms)
      do draw_index = 1, iterations + 1
         call random_standard_normal_matrix(initial_normals(:, :, draw_index))
         do time = 1, size(y)
            call random_standard_normal_matrix( &
               innovation_normals(:, :, time, draw_index))
         end do
      end do
      exponent = 0.6_dp
      if (present(adaptation_exponent)) exponent = adaptation_exponent
      if (present(target_acceptance)) then
         out = bssm_nonlinear_pmmh_draws(y, initial_mean, &
            initial_covariance, initial_parameters, noise_dimension, &
            log_density, transition, prior, proposal_factor, proposal_normals, &
            acceptance_uniforms, initial_normals, innovation_normals, &
            resampling_uniforms, target_acceptance, exponent)
      else
         out = bssm_nonlinear_pmmh_draws(y, initial_mean, &
            initial_covariance, initial_parameters, noise_dimension, &
            log_density, transition, prior, proposal_factor, proposal_normals, &
            acceptance_uniforms, initial_normals, innovation_normals, &
            resampling_uniforms)
      end if
   end function bssm_nonlinear_pmmh

   pure function bssm_multivariate_nonlinear_pmmh_draws(y, initial_mean, &
      initial_covariance, initial_parameters, observation_noise_dimension, &
      state_noise_dimension, observation, transition, prior, proposal_factor, &
      proposal_normals, acceptance_uniforms, initial_normals, &
      innovation_normals, resampling_uniforms, target_acceptance, &
      adaptation_exponent) result(out)
      !! Run draw-driven multivariate nonlinear bootstrap-filter PMMH.
      real(dp), intent(in) :: y(:, :) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: initial_parameters(:) !! Initial parameter values.
      real(dp), intent(in) :: proposal_factor(:, :) !! Proposal factor.
      real(dp), intent(in) :: proposal_normals(:, :) !! Standard-normal proposal draws.
      real(dp), intent(in) :: acceptance_uniforms(:) !! Acceptance uniforms.
      real(dp), intent(in) :: initial_normals(:, :, :) !! Initial normals.
      real(dp), intent(in) :: innovation_normals(:, :, :, :) !! Innovation normals.
      real(dp), intent(in) :: resampling_uniforms(:, :, :) !! Resampling uniforms.
      integer, intent(in) :: observation_noise_dimension !! Observation noise dimension.
      integer, intent(in) :: state_noise_dimension !! State noise dimension.
      procedure(bssm_multivariate_gaussian_observation_t) :: observation !! Observed value or vector.
      procedure(bssm_nonlinear_transition_t) :: transition !! State transition matrix.
      procedure(bssm_parameter_log_density_t) :: prior !! Prior-distribution specification.
      real(dp), intent(in), optional :: target_acceptance !! Target acceptance.
      real(dp), intent(in), optional :: adaptation_exponent !! Adaptation exponent.
      type(bssm_mcmc_t) :: out
      real(dp) :: exponent
      integer :: state, particles, times, iterations

      state = size(initial_mean)
      particles = size(initial_normals, 2)
      times = size(y, 2)
      iterations = size(proposal_normals, 2)
      if (size(y, 1) < 1 .or. state < 1 .or. particles < 1 .or. &
         times < 1 .or. iterations < 1 .or. &
         observation_noise_dimension < 1 .or. state_noise_dimension < 1 .or. &
         any(shape(initial_covariance) /= [state, state]) .or. &
         any(shape(initial_normals) /= &
         [state, particles, iterations + 1]) .or. &
         any(shape(innovation_normals) /= &
         [state_noise_dimension, particles, times, iterations + 1]) .or. &
         any(shape(resampling_uniforms) /= &
         [particles, times, iterations + 1]) .or. &
         .not. all(ieee_is_finite(initial_normals)) .or. &
         .not. all(ieee_is_finite(innovation_normals)) .or. &
         .not. all(ieee_is_finite(resampling_uniforms)) .or. &
         (present(adaptation_exponent) .and. &
         .not. present(target_acceptance))) then
         out%info = 1
         return
      end if
      exponent = 0.6_dp
      if (present(adaptation_exponent)) exponent = adaptation_exponent
      if (present(target_acceptance)) then
         out = bssm_pmmh_kernel_draws(initial_parameters, prior, estimator, &
            proposal_factor, proposal_normals, acceptance_uniforms, &
            target_acceptance, exponent)
      else
         out = bssm_pmmh_kernel_draws(initial_parameters, prior, estimator, &
            proposal_factor, proposal_normals, acceptance_uniforms)
      end if

   contains

      pure function estimator(parameters, draw_index) result(estimate)
         !! Estimate the multivariate nonlinear bootstrap log likelihood.
         real(dp), intent(in) :: parameters(:) !! Model parameter values.
         integer, intent(in) :: draw_index !! Index of draw.
         type(bssm_likelihood_estimate_t) :: estimate
         type(bssm_particle_filter_t) :: fit

         fit = bssm_multivariate_nonlinear_bootstrap_filter_draws(y, &
            initial_mean, initial_covariance, parameters, &
            observation_noise_dimension, state_noise_dimension, observation, &
            transition, initial_normals(:, :, draw_index), &
            innovation_normals(:, :, :, draw_index), &
            resampling_uniforms(:, :, draw_index))
         estimate%value = fit%log_likelihood
         estimate%info = fit%info
      end function estimator

   end function bssm_multivariate_nonlinear_pmmh_draws

   function bssm_multivariate_nonlinear_pmmh(y, initial_mean, &
      initial_covariance, initial_parameters, observation_noise_dimension, &
      state_noise_dimension, observation, transition, prior, proposal_factor, &
      particles, iterations, target_acceptance, adaptation_exponent) result(out)
      !! Run multivariate nonlinear bootstrap PMMH with shared randomness.
      real(dp), intent(in) :: y(:, :) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: initial_parameters(:) !! Initial parameter values.
      real(dp), intent(in) :: proposal_factor(:, :) !! Proposal factor.
      integer, intent(in) :: observation_noise_dimension !! Observation noise dimension.
      integer, intent(in) :: state_noise_dimension !! State noise dimension.
      integer, intent(in) :: particles !! Number of particles.
      integer, intent(in) :: iterations !! Number of algorithm iterations.
      procedure(bssm_multivariate_gaussian_observation_t) :: observation !! Observed value or vector.
      procedure(bssm_nonlinear_transition_t) :: transition !! State transition matrix.
      procedure(bssm_parameter_log_density_t) :: prior !! Prior-distribution specification.
      real(dp), intent(in), optional :: target_acceptance !! Target acceptance.
      real(dp), intent(in), optional :: adaptation_exponent !! Adaptation exponent.
      type(bssm_mcmc_t) :: out
      real(dp), allocatable :: proposal_normals(:, :), acceptance_uniforms(:)
      real(dp), allocatable :: initial_normals(:, :, :)
      real(dp), allocatable :: innovation_normals(:, :, :, :)
      real(dp), allocatable :: resampling_uniforms(:, :, :)
      real(dp) :: exponent
      integer :: draw_index, time

      if (size(y, 1) < 1 .or. size(y, 2) < 1 .or. &
         size(initial_mean) < 1 .or. size(initial_parameters) < 1 .or. &
         observation_noise_dimension < 1 .or. state_noise_dimension < 1 .or. &
         particles < 1 .or. iterations < 1 .or. &
         (present(adaptation_exponent) .and. &
         .not. present(target_acceptance))) then
         out%info = 1
         return
      end if
      allocate(proposal_normals(size(initial_parameters), iterations))
      allocate(acceptance_uniforms(iterations))
      allocate(initial_normals(size(initial_mean), particles, iterations + 1))
      allocate(innovation_normals(state_noise_dimension, particles, &
         size(y, 2), iterations + 1))
      allocate(resampling_uniforms(particles, size(y, 2), iterations + 1))
      call random_standard_normal_matrix(proposal_normals)
      call random_number(acceptance_uniforms)
      call random_number(resampling_uniforms)
      do draw_index = 1, iterations + 1
         call random_standard_normal_matrix(initial_normals(:, :, draw_index))
         do time = 1, size(y, 2)
            call random_standard_normal_matrix( &
               innovation_normals(:, :, time, draw_index))
         end do
      end do
      exponent = 0.6_dp
      if (present(adaptation_exponent)) exponent = adaptation_exponent
      if (present(target_acceptance)) then
         out = bssm_multivariate_nonlinear_pmmh_draws(y, initial_mean, &
            initial_covariance, initial_parameters, &
            observation_noise_dimension, state_noise_dimension, observation, &
            transition, prior, proposal_factor, proposal_normals, &
            acceptance_uniforms, initial_normals, innovation_normals, &
            resampling_uniforms, target_acceptance, exponent)
      else
         out = bssm_multivariate_nonlinear_pmmh_draws(y, initial_mean, &
            initial_covariance, initial_parameters, &
            observation_noise_dimension, state_noise_dimension, observation, &
            transition, prior, proposal_factor, proposal_normals, &
            acceptance_uniforms, initial_normals, innovation_normals, &
            resampling_uniforms)
      end if
   end function bssm_multivariate_nonlinear_pmmh

   pure function bssm_multivariate_nonlinear_da_pmmh_draws(y, initial_mean, &
      initial_covariance, initial_parameters, observation_noise_dimension, &
      state_noise_dimension, observation, transition_model, transition, prior, &
      iekf_iterations, convergence_tolerance, proposal_factor, &
      parameter_normals, first_stage_uniforms, second_stage_uniforms, &
      initial_normals, innovation_normals, resampling_uniforms, &
      target_acceptance, adaptation_exponent) result(out)
      !! Run draw-driven IEKF-screened multivariate bootstrap PMMH.
      real(dp), intent(in) :: y(:, :) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: initial_parameters(:) !! Initial parameter values.
      real(dp), intent(in) :: convergence_tolerance !! Convergence tolerance.
      real(dp), intent(in) :: proposal_factor(:, :) !! Proposal factor.
      real(dp), intent(in) :: parameter_normals(:, :) !! Parameter normals.
      real(dp), intent(in) :: first_stage_uniforms(:) !! First stage uniforms.
      real(dp), intent(in) :: second_stage_uniforms(:) !! Second stage uniforms.
      real(dp), intent(in) :: initial_normals(:, :, :) !! Initial normals.
      real(dp), intent(in) :: innovation_normals(:, :, :, :) !! Innovation normals.
      real(dp), intent(in) :: resampling_uniforms(:, :, :) !! Resampling uniforms.
      integer, intent(in) :: observation_noise_dimension !! Observation noise dimension.
      integer, intent(in) :: state_noise_dimension !! State noise dimension.
      integer, intent(in) :: iekf_iterations !! Number of iterated extended Kalman filter iterations.
      procedure(bssm_multivariate_gaussian_observation_t) :: observation !! Observed value or vector.
      procedure(bssm_nonlinear_transition_jacobian_t) :: transition_model !! Transition model callback procedure.
      procedure(bssm_nonlinear_transition_t) :: transition !! State transition matrix.
      procedure(bssm_parameter_log_density_t) :: prior !! Prior-distribution specification.
      real(dp), intent(in), optional :: target_acceptance !! Target acceptance.
      real(dp), intent(in), optional :: adaptation_exponent !! Adaptation exponent.
      type(bssm_da_mcmc_t) :: out
      real(dp) :: exponent
      integer :: state, particles, times, iterations

      state = size(initial_mean)
      particles = size(initial_normals, 2)
      times = size(y, 2)
      iterations = size(parameter_normals, 2)
      if (size(y, 1) < 1 .or. state < 1 .or. particles < 1 .or. &
         times < 1 .or. iterations < 1 .or. &
         observation_noise_dimension < 1 .or. state_noise_dimension < 1 .or. &
         iekf_iterations < 0 .or. convergence_tolerance < 0.0_dp .or. &
         any(shape(initial_covariance) /= [state, state]) .or. &
         any(shape(initial_normals) /= &
         [state, particles, iterations + 1]) .or. &
         any(shape(innovation_normals) /= &
         [state_noise_dimension, particles, times, iterations + 1]) .or. &
         any(shape(resampling_uniforms) /= &
         [particles, times, iterations + 1]) .or. &
         .not. all(ieee_is_finite(initial_normals)) .or. &
         .not. all(ieee_is_finite(innovation_normals)) .or. &
         .not. all(ieee_is_finite(resampling_uniforms)) .or. &
         (present(adaptation_exponent) .and. &
         .not. present(target_acceptance))) then
         out%info = 1
         return
      end if
      exponent = 0.6_dp
      if (present(adaptation_exponent)) exponent = adaptation_exponent
      if (present(target_acceptance)) then
         out = bssm_da_pmmh_kernel_draws(initial_parameters, prior, &
            coarse_estimator, fine_estimator, proposal_factor, &
            parameter_normals, first_stage_uniforms, second_stage_uniforms, &
            target_acceptance, exponent)
      else
         out = bssm_da_pmmh_kernel_draws(initial_parameters, prior, &
            coarse_estimator, fine_estimator, proposal_factor, &
            parameter_normals, first_stage_uniforms, second_stage_uniforms)
      end if

   contains

      pure function coarse_estimator(parameters, draw_index) result(estimate)
         !! Evaluate the deterministic multivariate IEKF likelihood.
         real(dp), intent(in) :: parameters(:) !! Model parameter values.
         integer, intent(in) :: draw_index !! Index of draw.
         type(bssm_likelihood_estimate_t) :: estimate
         type(bssm_multivariate_ekf_t) :: fit

         fit = bssm_multivariate_iekf(y, initial_mean, initial_covariance, &
            parameters, observation_noise_dimension, state_noise_dimension, &
            observation, transition_model, iekf_iterations, &
            convergence_tolerance)
         estimate%value = fit%log_likelihood + 0.0_dp*real(draw_index, dp)
         estimate%info = fit%info
      end function coarse_estimator

      pure function fine_estimator(parameters, draw_index) result(estimate)
         !! Evaluate the multivariate bootstrap particle likelihood.
         real(dp), intent(in) :: parameters(:) !! Model parameter values.
         integer, intent(in) :: draw_index !! Index of draw.
         type(bssm_likelihood_estimate_t) :: estimate
         type(bssm_particle_filter_t) :: fit

         fit = bssm_multivariate_nonlinear_bootstrap_filter_draws(y, &
            initial_mean, initial_covariance, parameters, &
            observation_noise_dimension, state_noise_dimension, observation, &
            transition, initial_normals(:, :, draw_index), &
            innovation_normals(:, :, :, draw_index), &
            resampling_uniforms(:, :, draw_index))
         estimate%value = fit%log_likelihood
         estimate%info = fit%info
      end function fine_estimator

   end function bssm_multivariate_nonlinear_da_pmmh_draws

   function bssm_multivariate_nonlinear_da_pmmh(y, initial_mean, &
      initial_covariance, initial_parameters, observation_noise_dimension, &
      state_noise_dimension, observation, transition_model, transition, prior, &
      iekf_iterations, convergence_tolerance, proposal_factor, particles, &
      iterations, target_acceptance, adaptation_exponent) result(out)
      !! Run IEKF-screened multivariate bootstrap PMMH with shared randomness.
      real(dp), intent(in) :: y(:, :) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: initial_parameters(:) !! Initial parameter values.
      real(dp), intent(in) :: convergence_tolerance !! Convergence tolerance.
      real(dp), intent(in) :: proposal_factor(:, :) !! Proposal factor.
      integer, intent(in) :: observation_noise_dimension !! Observation noise dimension.
      integer, intent(in) :: state_noise_dimension !! State noise dimension.
      integer, intent(in) :: iekf_iterations !! Number of iterated extended Kalman filter iterations.
      integer, intent(in) :: particles !! Number of particles.
      integer, intent(in) :: iterations !! Number of algorithm iterations.
      procedure(bssm_multivariate_gaussian_observation_t) :: observation !! Observed value or vector.
      procedure(bssm_nonlinear_transition_jacobian_t) :: transition_model !! Transition model callback procedure.
      procedure(bssm_nonlinear_transition_t) :: transition !! State transition matrix.
      procedure(bssm_parameter_log_density_t) :: prior !! Prior-distribution specification.
      real(dp), intent(in), optional :: target_acceptance !! Target acceptance.
      real(dp), intent(in), optional :: adaptation_exponent !! Adaptation exponent.
      type(bssm_da_mcmc_t) :: out
      real(dp), allocatable :: parameter_normals(:, :)
      real(dp), allocatable :: first_stage_uniforms(:), second_stage_uniforms(:)
      real(dp), allocatable :: initial_normals(:, :, :)
      real(dp), allocatable :: innovation_normals(:, :, :, :)
      real(dp), allocatable :: resampling_uniforms(:, :, :)
      real(dp) :: exponent
      integer :: draw_index, time

      if (size(y, 1) < 1 .or. size(y, 2) < 1 .or. &
         size(initial_mean) < 1 .or. size(initial_parameters) < 1 .or. &
         observation_noise_dimension < 1 .or. state_noise_dimension < 1 .or. &
         iekf_iterations < 0 .or. convergence_tolerance < 0.0_dp .or. &
         particles < 1 .or. iterations < 1 .or. &
         (present(adaptation_exponent) .and. &
         .not. present(target_acceptance))) then
         out%info = 1
         return
      end if
      allocate(parameter_normals(size(initial_parameters), iterations))
      allocate(first_stage_uniforms(iterations), second_stage_uniforms(iterations))
      allocate(initial_normals(size(initial_mean), particles, iterations + 1))
      allocate(innovation_normals(state_noise_dimension, particles, &
         size(y, 2), iterations + 1))
      allocate(resampling_uniforms(particles, size(y, 2), iterations + 1))
      call random_standard_normal_matrix(parameter_normals)
      call random_number(first_stage_uniforms)
      call random_number(second_stage_uniforms)
      call random_number(resampling_uniforms)
      do draw_index = 1, iterations + 1
         call random_standard_normal_matrix(initial_normals(:, :, draw_index))
         do time = 1, size(y, 2)
            call random_standard_normal_matrix( &
               innovation_normals(:, :, time, draw_index))
         end do
      end do
      exponent = 0.6_dp
      if (present(adaptation_exponent)) exponent = adaptation_exponent
      if (present(target_acceptance)) then
         out = bssm_multivariate_nonlinear_da_pmmh_draws(y, initial_mean, &
            initial_covariance, initial_parameters, &
            observation_noise_dimension, state_noise_dimension, observation, &
            transition_model, transition, prior, iekf_iterations, &
            convergence_tolerance, proposal_factor, parameter_normals, &
            first_stage_uniforms, second_stage_uniforms, initial_normals, &
            innovation_normals, resampling_uniforms, target_acceptance, &
            exponent)
      else
         out = bssm_multivariate_nonlinear_da_pmmh_draws(y, initial_mean, &
            initial_covariance, initial_parameters, &
            observation_noise_dimension, state_noise_dimension, observation, &
            transition_model, transition, prior, iekf_iterations, &
            convergence_tolerance, proposal_factor, parameter_normals, &
            first_stage_uniforms, second_stage_uniforms, initial_normals, &
            innovation_normals, resampling_uniforms)
      end if
   end function bssm_multivariate_nonlinear_da_pmmh

   pure function bssm_nonlinear_da_pmmh_draws(y, initial_mean, &
      initial_covariance, initial_parameters, noise_dimension, observation, &
      transition_model, log_density, transition, prior, &
      approximation_iterations, convergence_tolerance, proposal_factor, &
      parameter_normals, first_stage_uniforms, second_stage_uniforms, &
      initial_normals, innovation_normals, resampling_uniforms, &
      target_acceptance, adaptation_exponent) result(out)
      !! Run nonlinear delayed acceptance with Gaussian and bootstrap likelihoods.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: initial_parameters(:) !! Initial parameter values.
      real(dp), intent(in) :: convergence_tolerance !! Convergence tolerance.
      real(dp), intent(in) :: proposal_factor(:, :) !! Proposal factor.
      real(dp), intent(in) :: parameter_normals(:, :) !! Parameter normals.
      real(dp), intent(in) :: first_stage_uniforms(:) !! First stage uniforms.
      real(dp), intent(in) :: second_stage_uniforms(:) !! Second stage uniforms.
      real(dp), intent(in) :: initial_normals(:, :, :) !! Initial normals.
      real(dp), intent(in) :: innovation_normals(:, :, :, :) !! Innovation normals.
      real(dp), intent(in) :: resampling_uniforms(:, :, :) !! Resampling uniforms.
      integer, intent(in) :: noise_dimension !! Noise dimension.
      integer, intent(in) :: approximation_iterations !! Number of approximation iterations.
      procedure(bssm_gaussian_observation_t) :: observation !! Observed value or vector.
      procedure(bssm_nonlinear_transition_jacobian_t) :: transition_model !! Transition model callback procedure.
      procedure(bssm_nonlinear_log_density_t) :: log_density !! Log-density value.
      procedure(bssm_nonlinear_transition_t) :: transition !! State transition matrix.
      procedure(bssm_parameter_log_density_t) :: prior !! Prior-distribution specification.
      real(dp), intent(in), optional :: target_acceptance !! Target acceptance.
      real(dp), intent(in), optional :: adaptation_exponent !! Adaptation exponent.
      type(bssm_da_mcmc_t) :: out
      real(dp) :: exponent
      integer :: state, particles, times, iterations

      state = size(initial_mean)
      particles = size(initial_normals, 2)
      times = size(y)
      iterations = size(parameter_normals, 2)
      if (state < 1 .or. particles < 1 .or. times < 1 .or. iterations < 1 .or. &
         noise_dimension < 1 .or. approximation_iterations < 1 .or. &
         convergence_tolerance < 0.0_dp .or. &
         any(shape(initial_covariance) /= [state, state]) .or. &
         any(shape(initial_normals) /= &
         [state, particles, iterations + 1]) .or. &
         any(shape(innovation_normals) /= &
         [noise_dimension, particles, times, iterations + 1]) .or. &
         any(shape(resampling_uniforms) /= &
         [particles, times, iterations + 1]) .or. &
         .not. all(ieee_is_finite(initial_normals)) .or. &
         .not. all(ieee_is_finite(innovation_normals)) .or. &
         .not. all(ieee_is_finite(resampling_uniforms)) .or. &
         (present(adaptation_exponent) .and. &
         .not. present(target_acceptance))) then
         out%info = 1
         return
      end if
      exponent = 0.6_dp
      if (present(adaptation_exponent)) exponent = adaptation_exponent
      if (present(target_acceptance)) then
         out = bssm_da_pmmh_kernel_draws(initial_parameters, prior, &
            coarse_estimator, fine_estimator, proposal_factor, &
            parameter_normals, first_stage_uniforms, second_stage_uniforms, &
            target_acceptance, exponent)
      else
         out = bssm_da_pmmh_kernel_draws(initial_parameters, prior, &
            coarse_estimator, fine_estimator, proposal_factor, &
            parameter_normals, first_stage_uniforms, second_stage_uniforms)
      end if

   contains

      pure function coarse_estimator(parameters, draw_index) result(estimate)
         !! Evaluate the deterministic Gaussian-approximation likelihood.
         real(dp), intent(in) :: parameters(:) !! Model parameter values.
         integer, intent(in) :: draw_index !! Index of draw.
         type(bssm_likelihood_estimate_t) :: estimate
         type(bssm_nonlinear_approximation_t) :: approximation

         approximation = bssm_nonlinear_gaussian_approximation(y, &
            initial_mean, initial_covariance, parameters, noise_dimension, &
            observation, transition_model, approximation_iterations, &
            convergence_tolerance)
         estimate%value = approximation%corrected_log_likelihood + &
            0.0_dp*real(draw_index, dp)
         estimate%info = approximation%info
      end function coarse_estimator

      pure function fine_estimator(parameters, draw_index) result(estimate)
         !! Evaluate the nonlinear bootstrap particle likelihood.
         real(dp), intent(in) :: parameters(:) !! Model parameter values.
         integer, intent(in) :: draw_index !! Index of draw.
         type(bssm_likelihood_estimate_t) :: estimate
         type(bssm_particle_filter_t) :: fit

         fit = bssm_nonlinear_bootstrap_filter_draws(y, initial_mean, &
            initial_covariance, parameters, noise_dimension, log_density, &
            transition, initial_normals(:, :, draw_index), &
            innovation_normals(:, :, :, draw_index), &
            resampling_uniforms(:, :, draw_index))
         estimate%value = fit%log_likelihood
         estimate%info = fit%info
      end function fine_estimator

   end function bssm_nonlinear_da_pmmh_draws

   function bssm_nonlinear_da_pmmh(y, initial_mean, initial_covariance, &
      initial_parameters, noise_dimension, observation, transition_model, &
      log_density, transition, prior, approximation_iterations, &
      convergence_tolerance, proposal_factor, particles, iterations, &
      target_acceptance, adaptation_exponent) result(out)
      !! Run nonlinear Gaussian-screened bootstrap PMMH with shared randomness.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: initial_parameters(:) !! Initial parameter values.
      real(dp), intent(in) :: convergence_tolerance !! Convergence tolerance.
      real(dp), intent(in) :: proposal_factor(:, :) !! Proposal factor.
      integer, intent(in) :: noise_dimension !! Noise dimension.
      integer, intent(in) :: approximation_iterations !! Number of approximation iterations.
      integer, intent(in) :: particles !! Number of particles.
      integer, intent(in) :: iterations !! Number of algorithm iterations.
      procedure(bssm_gaussian_observation_t) :: observation !! Observed value or vector.
      procedure(bssm_nonlinear_transition_jacobian_t) :: transition_model !! Transition model callback procedure.
      procedure(bssm_nonlinear_log_density_t) :: log_density !! Log-density value.
      procedure(bssm_nonlinear_transition_t) :: transition !! State transition matrix.
      procedure(bssm_parameter_log_density_t) :: prior !! Prior-distribution specification.
      real(dp), intent(in), optional :: target_acceptance !! Target acceptance.
      real(dp), intent(in), optional :: adaptation_exponent !! Adaptation exponent.
      type(bssm_da_mcmc_t) :: out
      real(dp), allocatable :: parameter_normals(:, :)
      real(dp), allocatable :: first_stage_uniforms(:), second_stage_uniforms(:)
      real(dp), allocatable :: initial_normals(:, :, :)
      real(dp), allocatable :: innovation_normals(:, :, :, :)
      real(dp), allocatable :: resampling_uniforms(:, :, :)
      real(dp) :: exponent
      integer :: draw_index, time

      if (particles < 1 .or. iterations < 1 .or. noise_dimension < 1 .or. &
         size(y) < 1 .or. size(initial_mean) < 1 .or. &
         size(initial_parameters) < 1 .or. &
         (present(adaptation_exponent) .and. &
         .not. present(target_acceptance))) then
         out%info = 1
         return
      end if
      allocate(parameter_normals(size(initial_parameters), iterations))
      allocate(first_stage_uniforms(iterations), second_stage_uniforms(iterations))
      allocate(initial_normals(size(initial_mean), particles, iterations + 1))
      allocate(innovation_normals(noise_dimension, particles, size(y), &
         iterations + 1))
      allocate(resampling_uniforms(particles, size(y), iterations + 1))
      call random_standard_normal_matrix(parameter_normals)
      call random_number(first_stage_uniforms)
      call random_number(second_stage_uniforms)
      call random_number(resampling_uniforms)
      do draw_index = 1, iterations + 1
         call random_standard_normal_matrix(initial_normals(:, :, draw_index))
         do time = 1, size(y)
            call random_standard_normal_matrix( &
               innovation_normals(:, :, time, draw_index))
         end do
      end do
      exponent = 0.6_dp
      if (present(adaptation_exponent)) exponent = adaptation_exponent
      if (present(target_acceptance)) then
         out = bssm_nonlinear_da_pmmh_draws(y, initial_mean, &
            initial_covariance, initial_parameters, noise_dimension, &
            observation, transition_model, log_density, transition, prior, &
            approximation_iterations, convergence_tolerance, proposal_factor, &
            parameter_normals, first_stage_uniforms, second_stage_uniforms, &
            initial_normals, innovation_normals, resampling_uniforms, &
            target_acceptance, exponent)
      else
         out = bssm_nonlinear_da_pmmh_draws(y, initial_mean, &
            initial_covariance, initial_parameters, noise_dimension, &
            observation, transition_model, log_density, transition, prior, &
            approximation_iterations, convergence_tolerance, proposal_factor, &
            parameter_normals, first_stage_uniforms, second_stage_uniforms, &
            initial_normals, innovation_normals, resampling_uniforms)
      end if
   end function bssm_nonlinear_da_pmmh

   pure function bssm_nonlinear_gaussian_approximation(y, initial_mean, &
      initial_covariance, parameters, noise_dimension, observation, &
      transition_model, max_iterations, convergence_tolerance) result(out)
      !! Build a global iterated-EKS approximation to a nonlinear model.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: parameters(:) !! Model parameter values.
      integer, intent(in) :: noise_dimension !! Noise dimension.
      integer, intent(in) :: max_iterations !! Maximum number of algorithm iterations.
      real(dp), intent(in) :: convergence_tolerance !! Convergence tolerance.
      procedure(bssm_gaussian_observation_t) :: observation !! Observed value or vector.
      procedure(bssm_nonlinear_transition_jacobian_t) :: transition_model !! Transition model callback procedure.
      type(bssm_nonlinear_approximation_t) :: out
      type(ssm_model_t) :: model
      type(kfs_filter_t) :: filtered
      type(kfs_smoother_t) :: smoothed
      real(dp), allocatable :: mode(:, :), next_mode(:, :), center(:, :)
      real(dp), allocatable :: covariance(:, :), factor(:, :)
      real(dp), allocatable :: exact_mean(:), exact_loading(:, :)
      real(dp), allocatable :: observation_gradient(:)
      real(dp) :: observation_mean, observation_sd
      real(dp) :: exact_log_density, approximate_log_density
      integer :: state, times, time, component, iteration, info

      state = size(initial_mean)
      times = size(y)
      if (state < 1 .or. times < 1 .or. noise_dimension < 1 .or. &
         max_iterations < 1 .or. convergence_tolerance < 0.0_dp .or. &
         any(shape(initial_covariance) /= [state, state]) .or. &
         .not. all(ieee_is_finite(initial_mean)) .or. &
         .not. all(ieee_is_finite(initial_covariance)) .or. &
         .not. all(ieee_is_finite(parameters))) then
         out%info = 1
         return
      end if
      allocate(mode(state, times + 1), next_mode(state, times + 1))
      allocate(center(state, times), exact_mean(state))
      allocate(exact_loading(state, noise_dimension))
      allocate(observation_gradient(state))
      allocate(covariance(state, state))
      mode(:, 1) = initial_mean
      do time = 1, times
         call transition_model(time, mode(:, time), parameters, exact_mean, &
            covariance, exact_loading)
         if (.not. all(ieee_is_finite(exact_mean)) .or. &
            .not. all(ieee_is_finite(covariance)) .or. &
            .not. all(ieee_is_finite(exact_loading))) then
            out%info = 2
            return
         end if
         mode(:, time + 1) = exact_mean
      end do
      allocate(out%mode_state(state, times + 1))
      allocate(out%proposal_mean(state, times + 1))
      allocate(out%observation_intercept(times))
      allocate(out%observation_jacobian(state, times))
      allocate(out%observation_standard_deviation(times))
      allocate(out%transition_intercept(state, times))
      allocate(out%transition_jacobian(state, state, times))
      allocate(out%noise_loading(state, noise_dimension, times))
      allocate(out%proposal_factor(state, state, times + 1))
      allocate(out%conditional_matrix(state, state, times + 1))
      allocate(out%scaling(times))
      allocate(model%y(times, 1), model%z(1, state, times))
      allocate(model%h(1, 1, times))
      allocate(model%transition(state, state, times))
      allocate(model%r(state, noise_dimension, times))
      allocate(model%q(noise_dimension, noise_dimension, 1))
      allocate(model%a1(state), model%p1(state, state))
      allocate(model%p1inf(state, state), model%missing(times, 1))
      model%q = 0.0_dp
      do component = 1, noise_dimension
         model%q(component, component, 1) = 1.0_dp
      end do
      model%a1 = initial_mean
      model%p1 = initial_covariance
      model%p1inf = 0.0_dp
      model%missing(:, 1) = .not. ieee_is_finite(y)
      do iteration = 1, max_iterations
         out%iterations = iteration
         call linearize_model(mode, model, center, &
            out%observation_intercept, out%observation_jacobian, &
            out%observation_standard_deviation, out%transition_intercept, &
            out%transition_jacobian, out%noise_loading, info)
         if (info /= 0) then
            out%info = 10 + info
            return
         end if
         filtered = kfs_filter(model)
         if (filtered%info /= 0) then
            out%info = 100 + filtered%info
            return
         end if
         smoothed = kfs_smooth(model, filtered)
         if (smoothed%info /= 0) then
            out%info = 200 + smoothed%info
            return
         end if
         next_mode(:, 1:times) = smoothed%state + center
         call transition_model(times, next_mode(:, times), parameters, &
            next_mode(:, times + 1), covariance, exact_loading)
         if (.not. all(ieee_is_finite(next_mode(:, times + 1)))) then
            out%info = 3
            return
         end if
         out%difference = sum((next_mode(:, 1:times) - &
            mode(:, 1:times))**2)/real(state*times, dp)
         mode = next_mode
         if (out%difference <= convergence_tolerance) then
            out%converged = .true.
            exit
         end if
      end do
      call linearize_model(mode, model, center, &
         out%observation_intercept, out%observation_jacobian, &
         out%observation_standard_deviation, out%transition_intercept, &
         out%transition_jacobian, out%noise_loading, info)
      if (info /= 0) then
         out%info = 20 + info
         return
      end if
      filtered = kfs_filter(model)
      if (filtered%info /= 0) then
         out%info = 300 + filtered%info
         return
      end if
      smoothed = kfs_smooth(model, filtered)
      if (smoothed%info /= 0) then
         out%info = 400 + smoothed%info
         return
      end if
      out%mode_state = mode
      out%proposal_mean(:, 1:times) = smoothed%state + center
      out%proposal_mean(:, times + 1) = out%transition_intercept(:, times) + &
         matmul(out%transition_jacobian(:, :, times), &
         out%proposal_mean(:, times))
      out%gaussian_log_likelihood = filtered%log_likelihood
      out%conditional_matrix = 0.0_dp
      out%proposal_factor = 0.0_dp
      out%conditional_matrix(:, :, 1:times) = &
         smoothed%conditional_matrix
      do time = 1, times
         call positive_semidefinite_factor( &
            smoothed%conditional_covariance(:, :, time), factor, info)
         if (info /= 0) then
            out%info = 500 + time
            return
         end if
         out%proposal_factor(:, :, time) = factor
      end do
      out%scaling = 0.0_dp
      do time = 1, times
         call observation(time, mode(:, time), parameters, observation_mean, &
            observation_gradient, observation_sd)
         if (ieee_is_finite(y(time))) then
            exact_log_density = normal_log_density(y(time), observation_mean, &
               observation_sd)
            approximate_log_density = normal_log_density(y(time), &
               out%observation_intercept(time) + &
               dot_product(out%observation_jacobian(:, time), mode(:, time)), &
               out%observation_standard_deviation(time))
            out%scaling(time) = exact_log_density - approximate_log_density
         end if
         if (time > 1) then
            call transition_model(time - 1, mode(:, time - 1), parameters, &
               exact_mean, covariance, exact_loading)
            covariance = matmul(exact_loading, transpose(exact_loading))
            exact_log_density = multivariate_normal_log_density(mode(:, time), &
               exact_mean, covariance)
            covariance = matmul(out%noise_loading(:, :, time - 1), &
               transpose(out%noise_loading(:, :, time - 1)))
            exact_mean = out%transition_intercept(:, time - 1) + &
               matmul(out%transition_jacobian(:, :, time - 1), &
               mode(:, time - 1))
            approximate_log_density = multivariate_normal_log_density( &
               mode(:, time), exact_mean, covariance)
            out%scaling(time) = out%scaling(time) + exact_log_density - &
               approximate_log_density
         end if
      end do
      if (.not. all(ieee_is_finite(out%scaling))) then
         out%info = 6
         return
      end if
      out%corrected_log_likelihood = out%gaussian_log_likelihood + &
         sum(out%scaling)

   contains

      pure subroutine linearize_model(current_mode, model_work, center_work, &
         observation_intercept, observation_jacobian, observation_sigma, &
         transition_intercept, transition_jacobian, noise_loading, status)
         !! Linearize both nonlinear equations around one state trajectory.
         real(dp), intent(in) :: current_mode(:, :) !! Current mode.
         type(ssm_model_t), intent(inout) :: model_work !! Model work, updated in place.
         real(dp), intent(out) :: center_work(:, :) !! Center work.
         real(dp), intent(out) :: observation_intercept(:) !! Observation intercept.
         real(dp), intent(out) :: observation_jacobian(:, :) !! Observation jacobian.
         real(dp), intent(out) :: observation_sigma(:) !! Observation sigma.
         real(dp), intent(out) :: transition_intercept(:, :) !! Transition intercept.
         real(dp), intent(out) :: transition_jacobian(:, :, :) !! Transition jacobian.
         real(dp), intent(out) :: noise_loading(:, :, :) !! Noise loading.
         integer, intent(out) :: status !! Status.
         real(dp) :: mean, sigma
         real(dp) :: gradient(state), transition_mean(state)
         real(dp) :: jacobian(state, state)
         real(dp) :: loading(state, noise_dimension)
         integer :: t

         status = 0
         do t = 1, times
            call observation(t, current_mode(:, t), parameters, mean, &
               gradient, sigma)
            call transition_model(t, current_mode(:, t), parameters, &
               transition_mean, jacobian, loading)
            if (.not. ieee_is_finite(mean) .or. &
               .not. ieee_is_finite(sigma) .or. sigma <= 0.0_dp .or. &
               .not. all(ieee_is_finite(gradient)) .or. &
               .not. all(ieee_is_finite(transition_mean)) .or. &
               .not. all(ieee_is_finite(jacobian)) .or. &
               .not. all(ieee_is_finite(loading))) then
               status = 1
               return
            end if
            observation_jacobian(:, t) = gradient
            observation_intercept(t) = mean - &
               dot_product(gradient, current_mode(:, t))
            observation_sigma(t) = sigma
            transition_jacobian(:, :, t) = jacobian
            transition_intercept(:, t) = transition_mean - &
               matmul(jacobian, current_mode(:, t))
            noise_loading(:, :, t) = loading
         end do
         center_work(:, 1) = 0.0_dp
         do t = 1, times - 1
            center_work(:, t + 1) = transition_intercept(:, t) + &
               matmul(transition_jacobian(:, :, t), center_work(:, t))
         end do
         do t = 1, times
            model_work%z(1, :, t) = observation_jacobian(:, t)
            model_work%h(1, 1, t) = observation_sigma(t)**2
            model_work%transition(:, :, t) = transition_jacobian(:, :, t)
            model_work%r(:, :, t) = noise_loading(:, :, t)
            model_work%y(t, 1) = y(t) - observation_intercept(t) - &
               dot_product(observation_jacobian(:, t), center_work(:, t))
         end do
      end subroutine linearize_model

   end function bssm_nonlinear_gaussian_approximation

   pure function bssm_multivariate_nonlinear_gaussian_approximation(y, &
      initial_mean, initial_covariance, parameters, &
      observation_noise_dimension, state_noise_dimension, observation, &
      transition_model, max_iterations, convergence_tolerance) result(out)
      !! Build a global iterated-EKS approximation for vector observations.
      real(dp), intent(in) :: y(:, :) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: parameters(:) !! Model parameter values.
      integer, intent(in) :: observation_noise_dimension !! Observation noise dimension.
      integer, intent(in) :: state_noise_dimension !! State noise dimension.
      integer, intent(in) :: max_iterations !! Maximum number of algorithm iterations.
      real(dp), intent(in) :: convergence_tolerance !! Convergence tolerance.
      procedure(bssm_multivariate_gaussian_observation_t) :: observation !! Observed value or vector.
      procedure(bssm_nonlinear_transition_jacobian_t) :: transition_model !! Transition model callback procedure.
      type(bssm_multivariate_nonlinear_approximation_t) :: out
      type(ssm_model_t) :: model
      type(kfs_filter_t) :: filtered
      type(kfs_smoother_t) :: smoothed
      real(dp), allocatable :: mode(:, :), next_mode(:, :), center(:, :)
      real(dp), allocatable :: factor(:, :), observation_covariance(:, :)
      real(dp), allocatable :: state_covariance(:, :)
      real(dp), allocatable :: exact_observation_mean(:)
      real(dp), allocatable :: exact_observation_jacobian(:, :)
      real(dp), allocatable :: exact_observation_loading(:, :)
      real(dp), allocatable :: exact_transition_mean(:)
      real(dp), allocatable :: exact_transition_jacobian(:, :)
      real(dp), allocatable :: exact_state_loading(:, :)
      real(dp), allocatable :: approximate_mean(:)
      real(dp), allocatable :: approximate_transition_mean(:)
      real(dp) :: exact_log_density, approximate_log_density
      integer :: series, state, times, time, component, iteration, info

      series = size(y, 1)
      times = size(y, 2)
      state = size(initial_mean)
      if (series < 1 .or. state < 1 .or. times < 1 .or. &
         observation_noise_dimension < 1 .or. state_noise_dimension < 1 .or. &
         max_iterations < 1 .or. convergence_tolerance < 0.0_dp .or. &
         any(shape(initial_covariance) /= [state, state]) .or. &
         .not. all(ieee_is_finite(initial_mean)) .or. &
         .not. all(ieee_is_finite(initial_covariance)) .or. &
         .not. all(ieee_is_finite(parameters))) then
         out%info = 1
         return
      end if
      allocate(mode(state, times + 1), next_mode(state, times + 1))
      allocate(center(state, times), observation_covariance(series, series))
      allocate(state_covariance(state, state))
      allocate(exact_observation_mean(series))
      allocate(exact_observation_jacobian(series, state))
      allocate(exact_observation_loading(series, observation_noise_dimension))
      allocate(exact_transition_mean(state))
      allocate(exact_transition_jacobian(state, state))
      allocate(exact_state_loading(state, state_noise_dimension))
      allocate(approximate_mean(series))
      allocate(approximate_transition_mean(state))
      mode(:, 1) = initial_mean
      do time = 1, times
         call transition_model(time, mode(:, time), parameters, &
            exact_transition_mean, exact_transition_jacobian, &
            exact_state_loading)
         if (.not. all(ieee_is_finite(exact_transition_mean)) .or. &
            .not. all(ieee_is_finite(exact_transition_jacobian)) .or. &
            .not. all(ieee_is_finite(exact_state_loading))) then
            out%info = 2
            return
         end if
         mode(:, time + 1) = exact_transition_mean
      end do
      allocate(out%mode_state(state, times + 1))
      allocate(out%proposal_mean(state, times + 1))
      allocate(out%observation_intercept(series, times))
      allocate(out%observation_jacobian(series, state, times))
      allocate(out%observation_noise_loading(series, &
         observation_noise_dimension, times))
      allocate(out%transition_intercept(state, times))
      allocate(out%transition_jacobian(state, state, times))
      allocate(out%state_noise_loading(state, state_noise_dimension, times))
      allocate(out%proposal_factor(state, state, times + 1))
      allocate(out%conditional_matrix(state, state, times + 1))
      allocate(out%scaling(times))
      allocate(model%y(times, series), model%z(series, state, times))
      allocate(model%h(series, series, times))
      allocate(model%transition(state, state, times))
      allocate(model%r(state, state_noise_dimension, times))
      allocate(model%q(state_noise_dimension, state_noise_dimension, 1))
      allocate(model%a1(state), model%p1(state, state))
      allocate(model%p1inf(state, state), model%missing(times, series))
      model%q = 0.0_dp
      do component = 1, state_noise_dimension
         model%q(component, component, 1) = 1.0_dp
      end do
      model%a1 = initial_mean
      model%p1 = initial_covariance
      model%p1inf = 0.0_dp
      model%missing = transpose(.not. ieee_is_finite(y))
      do iteration = 1, max_iterations
         out%iterations = iteration
         call linearize_model(mode, model, center, &
            out%observation_intercept, out%observation_jacobian, &
            out%observation_noise_loading, out%transition_intercept, &
            out%transition_jacobian, out%state_noise_loading, info)
         if (info /= 0) then
            out%info = 10 + info
            return
         end if
         filtered = kfs_filter(model)
         if (filtered%info /= 0) then
            out%info = 100 + filtered%info
            return
         end if
         smoothed = kfs_smooth(model, filtered)
         if (smoothed%info /= 0) then
            out%info = 200 + smoothed%info
            return
         end if
         next_mode(:, 1:times) = smoothed%state + center
         call transition_model(times, next_mode(:, times), parameters, &
            next_mode(:, times + 1), exact_transition_jacobian, &
            exact_state_loading)
         if (.not. all(ieee_is_finite(next_mode(:, times + 1)))) then
            out%info = 3
            return
         end if
         out%difference = sum((next_mode(:, 1:times) - &
            mode(:, 1:times))**2)/real(state*times, dp)
         mode = next_mode
         if (out%difference <= convergence_tolerance) then
            out%converged = .true.
            exit
         end if
      end do
      call linearize_model(mode, model, center, &
         out%observation_intercept, out%observation_jacobian, &
         out%observation_noise_loading, out%transition_intercept, &
         out%transition_jacobian, out%state_noise_loading, info)
      if (info /= 0) then
         out%info = 20 + info
         return
      end if
      filtered = kfs_filter(model)
      if (filtered%info /= 0) then
         out%info = 300 + filtered%info
         return
      end if
      smoothed = kfs_smooth(model, filtered)
      if (smoothed%info /= 0) then
         out%info = 400 + smoothed%info
         return
      end if
      out%mode_state = mode
      out%proposal_mean(:, 1:times) = smoothed%state + center
      out%proposal_mean(:, times + 1) = out%transition_intercept(:, times) + &
         matmul(out%transition_jacobian(:, :, times), &
         out%proposal_mean(:, times))
      out%gaussian_log_likelihood = filtered%log_likelihood
      out%conditional_matrix = 0.0_dp
      out%proposal_factor = 0.0_dp
      out%conditional_matrix(:, :, 1:times) = smoothed%conditional_matrix
      do time = 1, times
         call positive_semidefinite_factor( &
            smoothed%conditional_covariance(:, :, time), factor, info)
         if (info /= 0) then
            out%info = 500 + time
            return
         end if
         out%proposal_factor(:, :, time) = factor
      end do
      out%scaling = 0.0_dp
      do time = 1, times
         call observation(time, mode(:, time), parameters, &
            exact_observation_mean, exact_observation_jacobian, &
            exact_observation_loading)
         observation_covariance = matmul(exact_observation_loading, &
            transpose(exact_observation_loading))
         exact_log_density = observed_log_density(y(:, time), &
            exact_observation_mean, observation_covariance)
         approximate_mean = out%observation_intercept(:, time) + &
            matmul(out%observation_jacobian(:, :, time), mode(:, time))
         observation_covariance = matmul( &
            out%observation_noise_loading(:, :, time), &
            transpose(out%observation_noise_loading(:, :, time)))
         approximate_log_density = observed_log_density(y(:, time), &
            approximate_mean, observation_covariance)
         out%scaling(time) = exact_log_density - approximate_log_density
         if (time > 1) then
            call transition_model(time - 1, mode(:, time - 1), parameters, &
               exact_transition_mean, exact_transition_jacobian, &
               exact_state_loading)
            state_covariance = matmul(exact_state_loading, &
               transpose(exact_state_loading))
            exact_log_density = multivariate_normal_log_density( &
               mode(:, time), exact_transition_mean, state_covariance)
            approximate_transition_mean = &
               out%transition_intercept(:, time - 1) + &
               matmul(out%transition_jacobian(:, :, time - 1), &
               mode(:, time - 1))
            state_covariance = matmul( &
               out%state_noise_loading(:, :, time - 1), &
               transpose(out%state_noise_loading(:, :, time - 1)))
            approximate_log_density = multivariate_normal_log_density( &
               mode(:, time), approximate_transition_mean, state_covariance)
            out%scaling(time) = out%scaling(time) + exact_log_density - &
               approximate_log_density
         end if
      end do
      if (.not. all(ieee_is_finite(out%scaling))) then
         out%info = 6
         return
      end if
      out%corrected_log_likelihood = out%gaussian_log_likelihood + &
         sum(out%scaling)

   contains

      pure subroutine linearize_model(current_mode, model_work, center_work, &
         observation_intercept, observation_jacobian, observation_loading, &
         transition_intercept, transition_jacobian, state_loading, status)
         !! Linearize vector observation and transition equations.
         real(dp), intent(in) :: current_mode(:, :) !! Current mode.
         type(ssm_model_t), intent(inout) :: model_work !! Model work, updated in place.
         real(dp), intent(out) :: center_work(:, :) !! Center work.
         real(dp), intent(out) :: observation_intercept(:, :) !! Observation intercept.
         real(dp), intent(out) :: observation_jacobian(:, :, :) !! Observation jacobian.
         real(dp), intent(out) :: observation_loading(:, :, :) !! Observation loading matrix.
         real(dp), intent(out) :: transition_intercept(:, :) !! Transition intercept.
         real(dp), intent(out) :: transition_jacobian(:, :, :) !! Transition jacobian.
         real(dp), intent(out) :: state_loading(:, :, :) !! State loading.
         integer, intent(out) :: status !! Status.
         real(dp) :: mean(series), jacobian(series, state)
         real(dp) :: obs_loading(series, observation_noise_dimension)
         real(dp) :: transition_mean(state), transition_derivative(state, state)
         real(dp) :: noise_loading(state, state_noise_dimension)
         integer :: t

         status = 0
         do t = 1, times
            call observation(t, current_mode(:, t), parameters, mean, &
               jacobian, obs_loading)
            call transition_model(t, current_mode(:, t), parameters, &
               transition_mean, transition_derivative, noise_loading)
            if (.not. all(ieee_is_finite(mean)) .or. &
               .not. all(ieee_is_finite(jacobian)) .or. &
               .not. all(ieee_is_finite(obs_loading)) .or. &
               .not. all(ieee_is_finite(transition_mean)) .or. &
               .not. all(ieee_is_finite(transition_derivative)) .or. &
               .not. all(ieee_is_finite(noise_loading))) then
               status = 1
               return
            end if
            observation_jacobian(:, :, t) = jacobian
            observation_intercept(:, t) = mean - &
               matmul(jacobian, current_mode(:, t))
            observation_loading(:, :, t) = obs_loading
            transition_jacobian(:, :, t) = transition_derivative
            transition_intercept(:, t) = transition_mean - &
               matmul(transition_derivative, current_mode(:, t))
            state_loading(:, :, t) = noise_loading
         end do
         center_work(:, 1) = 0.0_dp
         do t = 1, times - 1
            center_work(:, t + 1) = transition_intercept(:, t) + &
               matmul(transition_jacobian(:, :, t), center_work(:, t))
         end do
         do t = 1, times
            model_work%z(:, :, t) = observation_jacobian(:, :, t)
            model_work%h(:, :, t) = matmul(observation_loading(:, :, t), &
               transpose(observation_loading(:, :, t)))
            model_work%transition(:, :, t) = transition_jacobian(:, :, t)
            model_work%r(:, :, t) = state_loading(:, :, t)
            model_work%y(t, :) = y(:, t) - observation_intercept(:, t) - &
               matmul(observation_jacobian(:, :, t), center_work(:, t))
         end do
      end subroutine linearize_model

      pure function observed_log_density(value, mean, covariance) result(value_log)
         !! Evaluate a Gaussian vector density over finite components.
         real(dp), intent(in) :: value(:) !! Input value.
         real(dp), intent(in) :: mean(:) !! Mean value or vector.
         real(dp), intent(in) :: covariance(:, :) !! Covariance matrix.
         real(dp) :: value_log
         real(dp), allocatable :: compact_covariance(:, :)
         integer, allocatable :: observed(:)
         integer :: observed_count, index, row, column

         observed_count = count(ieee_is_finite(value))
         if (observed_count == 0) then
            value_log = 0.0_dp
            return
         end if
         allocate(observed(observed_count))
         observed = pack([(index, index = 1, size(value))], &
            ieee_is_finite(value))
         allocate(compact_covariance(observed_count, observed_count))
         do row = 1, observed_count
            do column = 1, observed_count
               compact_covariance(row, column) = &
                  covariance(observed(row), observed(column))
            end do
         end do
         value_log = multivariate_normal_log_density(value(observed), &
            mean(observed), compact_covariance)
      end function observed_log_density

   end function bssm_multivariate_nonlinear_gaussian_approximation

   pure function bssm_nonlinear_approximate_mcmc_draws(y, initial_mean, &
      initial_covariance, initial_parameters, noise_dimension, observation, &
      transition_model, prior, approximation_iterations, &
      convergence_tolerance, proposal_factor, proposal_normals, &
      acceptance_uniforms, target_acceptance, adaptation_exponent) result(out)
      !! Run MCMC using the nonlinear Gaussian-approximation likelihood.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: initial_parameters(:) !! Initial parameter values.
      real(dp), intent(in) :: convergence_tolerance !! Convergence tolerance.
      real(dp), intent(in) :: proposal_factor(:, :) !! Proposal factor.
      real(dp), intent(in) :: proposal_normals(:, :) !! Standard-normal proposal draws.
      real(dp), intent(in) :: acceptance_uniforms(:) !! Acceptance uniforms.
      integer, intent(in) :: noise_dimension !! Noise dimension.
      integer, intent(in) :: approximation_iterations !! Number of approximation iterations.
      procedure(bssm_gaussian_observation_t) :: observation !! Observed value or vector.
      procedure(bssm_nonlinear_transition_jacobian_t) :: transition_model !! Transition model callback procedure.
      procedure(bssm_parameter_log_density_t) :: prior !! Prior-distribution specification.
      real(dp), intent(in), optional :: target_acceptance !! Target acceptance.
      real(dp), intent(in), optional :: adaptation_exponent !! Adaptation exponent.
      type(bssm_mcmc_t) :: out
      real(dp) :: exponent

      if (size(y) < 1 .or. size(initial_mean) < 1 .or. &
         noise_dimension < 1 .or. approximation_iterations < 1 .or. &
         convergence_tolerance < 0.0_dp .or. &
         any(shape(initial_covariance) /= &
         [size(initial_mean), size(initial_mean)]) .or. &
         (present(adaptation_exponent) .and. &
         .not. present(target_acceptance))) then
         out%info = 1
         return
      end if
      exponent = 0.6_dp
      if (present(adaptation_exponent)) exponent = adaptation_exponent
      if (present(target_acceptance)) then
         out = bssm_approximate_mcmc_draws(initial_parameters, prior, &
            approximate_likelihood, proposal_factor, proposal_normals, &
            acceptance_uniforms, target_acceptance, exponent)
      else
         out = bssm_approximate_mcmc_draws(initial_parameters, prior, &
            approximate_likelihood, proposal_factor, proposal_normals, &
            acceptance_uniforms)
      end if

   contains

      pure function approximate_likelihood(parameters) result(value)
         !! Evaluate the nonlinear Gaussian-approximation likelihood.
         real(dp), intent(in) :: parameters(:) !! Model parameter values.
         real(dp) :: value
         type(bssm_nonlinear_approximation_t) :: approximation

         approximation = bssm_nonlinear_gaussian_approximation(y, &
            initial_mean, initial_covariance, parameters, noise_dimension, &
            observation, transition_model, approximation_iterations, &
            convergence_tolerance)
         if (approximation%info == 0) then
            value = approximation%corrected_log_likelihood
         else
            value = -huge(1.0_dp)
         end if
      end function approximate_likelihood

   end function bssm_nonlinear_approximate_mcmc_draws

   function bssm_nonlinear_approximate_mcmc(y, initial_mean, &
      initial_covariance, initial_parameters, noise_dimension, observation, &
      transition_model, prior, approximation_iterations, &
      convergence_tolerance, proposal_factor, iterations, target_acceptance, &
      adaptation_exponent) result(out)
      !! Run nonlinear Gaussian-approximation MCMC with shared randomness.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: initial_parameters(:) !! Initial parameter values.
      real(dp), intent(in) :: convergence_tolerance !! Convergence tolerance.
      real(dp), intent(in) :: proposal_factor(:, :) !! Proposal factor.
      integer, intent(in) :: noise_dimension !! Noise dimension.
      integer, intent(in) :: approximation_iterations !! Number of approximation iterations.
      integer, intent(in) :: iterations !! Number of algorithm iterations.
      procedure(bssm_gaussian_observation_t) :: observation !! Observed value or vector.
      procedure(bssm_nonlinear_transition_jacobian_t) :: transition_model !! Transition model callback procedure.
      procedure(bssm_parameter_log_density_t) :: prior !! Prior-distribution specification.
      real(dp), intent(in), optional :: target_acceptance !! Target acceptance.
      real(dp), intent(in), optional :: adaptation_exponent !! Adaptation exponent.
      type(bssm_mcmc_t) :: out
      real(dp), allocatable :: proposal_normals(:, :), acceptance_uniforms(:)
      real(dp) :: exponent

      if (size(initial_parameters) < 1 .or. iterations < 1 .or. &
         (present(adaptation_exponent) .and. &
         .not. present(target_acceptance))) then
         out%info = 1
         return
      end if
      allocate(proposal_normals(size(initial_parameters), iterations))
      allocate(acceptance_uniforms(iterations))
      call random_standard_normal_matrix(proposal_normals)
      call random_number(acceptance_uniforms)
      exponent = 0.6_dp
      if (present(adaptation_exponent)) exponent = adaptation_exponent
      if (present(target_acceptance)) then
         out = bssm_nonlinear_approximate_mcmc_draws(y, initial_mean, &
            initial_covariance, initial_parameters, noise_dimension, &
            observation, transition_model, prior, approximation_iterations, &
            convergence_tolerance, proposal_factor, proposal_normals, &
            acceptance_uniforms, target_acceptance, exponent)
      else
         out = bssm_nonlinear_approximate_mcmc_draws(y, initial_mean, &
            initial_covariance, initial_parameters, noise_dimension, &
            observation, transition_model, prior, approximation_iterations, &
            convergence_tolerance, proposal_factor, proposal_normals, &
            acceptance_uniforms)
      end if
   end function bssm_nonlinear_approximate_mcmc

   pure function bssm_nonlinear_ekf_mcmc_draws(y, initial_mean, &
      initial_covariance, initial_parameters, noise_dimension, observation, &
      transition_model, prior, iekf_iterations, convergence_tolerance, &
      proposal_factor, proposal_normals, acceptance_uniforms, &
      target_acceptance, adaptation_exponent) result(out)
      !! Run MCMC using an extended-Kalman approximate likelihood.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: initial_parameters(:) !! Initial parameter values.
      real(dp), intent(in) :: convergence_tolerance !! Convergence tolerance.
      real(dp), intent(in) :: proposal_factor(:, :) !! Proposal factor.
      real(dp), intent(in) :: proposal_normals(:, :) !! Standard-normal proposal draws.
      real(dp), intent(in) :: acceptance_uniforms(:) !! Acceptance uniforms.
      integer, intent(in) :: noise_dimension !! Noise dimension.
      integer, intent(in) :: iekf_iterations !! Number of iterated extended Kalman filter iterations.
      procedure(bssm_gaussian_observation_t) :: observation !! Observed value or vector.
      procedure(bssm_nonlinear_transition_jacobian_t) :: transition_model !! Transition model callback procedure.
      procedure(bssm_parameter_log_density_t) :: prior !! Prior-distribution specification.
      real(dp), intent(in), optional :: target_acceptance !! Target acceptance.
      real(dp), intent(in), optional :: adaptation_exponent !! Adaptation exponent.
      type(bssm_mcmc_t) :: out
      real(dp) :: exponent

      if (size(y) < 1 .or. size(initial_mean) < 1 .or. &
         noise_dimension < 1 .or. iekf_iterations < 0 .or. &
         convergence_tolerance < 0.0_dp .or. &
         any(shape(initial_covariance) /= &
         [size(initial_mean), size(initial_mean)]) .or. &
         (present(adaptation_exponent) .and. &
         .not. present(target_acceptance))) then
         out%info = 1
         return
      end if
      exponent = 0.6_dp
      if (present(adaptation_exponent)) exponent = adaptation_exponent
      if (present(target_acceptance)) then
         out = bssm_approximate_mcmc_draws(initial_parameters, prior, &
            approximate_likelihood, proposal_factor, proposal_normals, &
            acceptance_uniforms, target_acceptance, exponent)
      else
         out = bssm_approximate_mcmc_draws(initial_parameters, prior, &
            approximate_likelihood, proposal_factor, proposal_normals, &
            acceptance_uniforms)
      end if

   contains

      pure function approximate_likelihood(parameters) result(value)
         !! Evaluate the extended-Kalman approximate likelihood.
         real(dp), intent(in) :: parameters(:) !! Model parameter values.
         real(dp) :: value
         type(bssm_ekf_t) :: fit

         fit = bssm_iekf(y, initial_mean, initial_covariance, parameters, &
            noise_dimension, observation, transition_model, iekf_iterations, &
            convergence_tolerance)
         if (fit%info == 0) then
            value = fit%log_likelihood
         else
            value = -huge(1.0_dp)
         end if
      end function approximate_likelihood

   end function bssm_nonlinear_ekf_mcmc_draws

   function bssm_nonlinear_ekf_mcmc(y, initial_mean, initial_covariance, &
      initial_parameters, noise_dimension, observation, transition_model, &
      prior, iekf_iterations, convergence_tolerance, proposal_factor, &
      iterations, target_acceptance, adaptation_exponent) result(out)
      !! Run extended-Kalman approximate MCMC with shared randomness.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: initial_parameters(:) !! Initial parameter values.
      real(dp), intent(in) :: convergence_tolerance !! Convergence tolerance.
      real(dp), intent(in) :: proposal_factor(:, :) !! Proposal factor.
      integer, intent(in) :: noise_dimension !! Noise dimension.
      integer, intent(in) :: iekf_iterations !! Number of iterated extended Kalman filter iterations.
      integer, intent(in) :: iterations !! Number of algorithm iterations.
      procedure(bssm_gaussian_observation_t) :: observation !! Observed value or vector.
      procedure(bssm_nonlinear_transition_jacobian_t) :: transition_model !! Transition model callback procedure.
      procedure(bssm_parameter_log_density_t) :: prior !! Prior-distribution specification.
      real(dp), intent(in), optional :: target_acceptance !! Target acceptance.
      real(dp), intent(in), optional :: adaptation_exponent !! Adaptation exponent.
      type(bssm_mcmc_t) :: out
      real(dp), allocatable :: proposal_normals(:, :), acceptance_uniforms(:)
      real(dp) :: exponent

      if (size(initial_parameters) < 1 .or. iterations < 1 .or. &
         (present(adaptation_exponent) .and. &
         .not. present(target_acceptance))) then
         out%info = 1
         return
      end if
      allocate(proposal_normals(size(initial_parameters), iterations))
      allocate(acceptance_uniforms(iterations))
      call random_standard_normal_matrix(proposal_normals)
      call random_number(acceptance_uniforms)
      exponent = 0.6_dp
      if (present(adaptation_exponent)) exponent = adaptation_exponent
      if (present(target_acceptance)) then
         out = bssm_nonlinear_ekf_mcmc_draws(y, initial_mean, &
            initial_covariance, initial_parameters, noise_dimension, &
            observation, transition_model, prior, iekf_iterations, &
            convergence_tolerance, proposal_factor, proposal_normals, &
            acceptance_uniforms, target_acceptance, exponent)
      else
         out = bssm_nonlinear_ekf_mcmc_draws(y, initial_mean, &
            initial_covariance, initial_parameters, noise_dimension, &
            observation, transition_model, prior, iekf_iterations, &
            convergence_tolerance, proposal_factor, proposal_normals, &
            acceptance_uniforms)
      end if
   end function bssm_nonlinear_ekf_mcmc

   pure function bssm_nonlinear_psi_filter_draws(y, initial_mean, &
      initial_covariance, parameters, noise_dimension, observation, &
      transition_model, approximation_iterations, convergence_tolerance, &
      initial_normals, proposal_normals, terminal_normals, &
      resampling_uniforms) result(out)
      !! Run a draw-driven psi filter for a nonlinear Gaussian model.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: parameters(:) !! Model parameter values.
      integer, intent(in) :: noise_dimension !! Noise dimension.
      integer, intent(in) :: approximation_iterations !! Number of approximation iterations.
      real(dp), intent(in) :: convergence_tolerance !! Convergence tolerance.
      real(dp), intent(in) :: initial_normals(:, :) !! Initial normals.
      real(dp), intent(in) :: proposal_normals(:, :, :) !! Standard-normal proposal draws.
      real(dp), intent(in) :: terminal_normals(:, :) !! Terminal normals.
      real(dp), intent(in) :: resampling_uniforms(:, :) !! Resampling uniforms.
      procedure(bssm_gaussian_observation_t) :: observation !! Observed value or vector.
      procedure(bssm_nonlinear_transition_jacobian_t) :: transition_model !! Transition model callback procedure.
      type(bssm_particle_filter_t) :: out
      type(bssm_nonlinear_approximation_t) :: approximation
      real(dp), allocatable :: parent(:, :), probability(:), log_weight(:)
      real(dp), allocatable :: exact_mean(:), approximate_mean(:)
      real(dp), allocatable :: jacobian(:, :), loading(:, :), covariance(:, :)
      real(dp), allocatable :: observation_gradient(:)
      real(dp) :: observation_mean, observation_sd
      real(dp) :: exact_density, approximate_density, increment
      integer :: state, particles, times, particle, time, info

      state = size(initial_mean)
      particles = size(initial_normals, 2)
      times = size(y)
      if (state < 1 .or. particles < 1 .or. times < 1 .or. &
         noise_dimension < 1 .or. &
         any(shape(initial_covariance) /= [state, state]) .or. &
         size(initial_normals, 1) /= state .or. &
         any(shape(proposal_normals) /= [state, particles, times]) .or. &
         any(shape(terminal_normals) /= [noise_dimension, particles]) .or. &
         any(shape(resampling_uniforms) /= [particles, times]) .or. &
         .not. all(ieee_is_finite(initial_normals)) .or. &
         .not. all(ieee_is_finite(proposal_normals)) .or. &
         .not. all(ieee_is_finite(terminal_normals)) .or. &
         .not. all(ieee_is_finite(resampling_uniforms))) then
         out%info = 1
         return
      end if
      approximation = bssm_nonlinear_gaussian_approximation(y, initial_mean, &
         initial_covariance, parameters, noise_dimension, observation, &
         transition_model, approximation_iterations, convergence_tolerance)
      if (approximation%info /= 0) then
         out%info = 1000 + approximation%info
         return
      end if
      allocate(out%particles(state, particles, times + 1))
      allocate(out%weights(particles, times), out%ancestors(particles, times))
      allocate(out%predicted_mean(state, times + 1))
      allocate(out%filtered_mean(state, times))
      allocate(out%predicted_covariance(state, state, times + 1))
      allocate(out%filtered_covariance(state, state, times))
      allocate(parent(state, particles), probability(particles))
      allocate(log_weight(particles), exact_mean(state))
      allocate(approximate_mean(state), jacobian(state, state))
      allocate(loading(state, noise_dimension), covariance(state, state))
      allocate(observation_gradient(state))
      do particle = 1, particles
         out%particles(:, particle, 1) = approximation%proposal_mean(:, 1) + &
            matmul(approximation%proposal_factor(:, :, 1), &
            initial_normals(:, particle))
      end do
      out%log_likelihood = approximation%corrected_log_likelihood
      parent = 0.0_dp
      do time = 1, times
         call particle_summary(out%particles(:, :, time), &
            out%predicted_mean(:, time), &
            out%predicted_covariance(:, :, time))
         do particle = 1, particles
            log_weight(particle) = 0.0_dp
            call observation(time, out%particles(:, particle, time), &
               parameters, observation_mean, observation_gradient, &
               observation_sd)
            if (.not. ieee_is_finite(observation_mean) .or. &
               .not. ieee_is_finite(observation_sd) .or. &
               observation_sd <= 0.0_dp) then
               out%info = 2
               return
            end if
            if (ieee_is_finite(y(time))) then
               exact_density = normal_log_density(y(time), observation_mean, &
                  observation_sd)
               approximate_density = normal_log_density(y(time), &
                  approximation%observation_intercept(time) + &
                  dot_product(approximation%observation_jacobian(:, time), &
                  out%particles(:, particle, time)), &
                  approximation%observation_standard_deviation(time))
               log_weight(particle) = exact_density - approximate_density
            end if
            if (time > 1) then
               call transition_model(time - 1, parent(:, particle), &
                  parameters, exact_mean, jacobian, loading)
               if (.not. all(ieee_is_finite(exact_mean)) .or. &
                  .not. all(ieee_is_finite(loading))) then
                  out%info = 3
                  return
               end if
               covariance = matmul(loading, transpose(loading))
               exact_density = multivariate_normal_log_density( &
                  out%particles(:, particle, time), exact_mean, covariance)
               approximate_mean = &
                  approximation%transition_intercept(:, time - 1) + &
                  matmul(approximation%transition_jacobian(:, :, time - 1), &
                  parent(:, particle))
               covariance = matmul( &
                  approximation%noise_loading(:, :, time - 1), &
                  transpose(approximation%noise_loading(:, :, time - 1)))
               approximate_density = multivariate_normal_log_density( &
                  out%particles(:, particle, time), approximate_mean, &
                  covariance)
               log_weight(particle) = log_weight(particle) + exact_density - &
                  approximate_density
            end if
            log_weight(particle) = log_weight(particle) - &
               approximation%scaling(time)
         end do
         call normalize_log_weights(log_weight, probability, increment, info)
         if (info /= 0) then
            out%info = 4
            out%log_likelihood = -huge(1.0_dp)
            return
         end if
         out%log_likelihood = out%log_likelihood + increment
         out%weights(:, time) = probability
         call weighted_particle_summary(out%particles(:, :, time), &
            probability, out%filtered_mean(:, time), &
            out%filtered_covariance(:, :, time))
         out%ancestors(:, time) = bssm_stratified_resample(probability, &
            resampling_uniforms(:, time))
         do particle = 1, particles
            parent(:, particle) = out%particles(:, &
               out%ancestors(particle, time), time)
            if (time < times) then
               out%particles(:, particle, time + 1) = &
                  approximation%proposal_mean(:, time + 1) + &
                  matmul(approximation%conditional_matrix(:, :, time + 1), &
                  parent(:, particle) - &
                  approximation%proposal_mean(:, time)) + &
                  matmul(approximation%proposal_factor(:, :, time + 1), &
                  proposal_normals(:, particle, time))
            else
               call transition_model(time, parent(:, particle), parameters, &
                  exact_mean, jacobian, loading)
               if (.not. all(ieee_is_finite(exact_mean)) .or. &
                  .not. all(ieee_is_finite(loading))) then
                  out%info = 3
                  return
               end if
               out%particles(:, particle, time + 1) = exact_mean + &
                  matmul(loading, terminal_normals(:, particle))
            end if
         end do
      end do
      call particle_summary(out%particles(:, :, times + 1), &
         out%predicted_mean(:, times + 1), &
         out%predicted_covariance(:, :, times + 1))
   end function bssm_nonlinear_psi_filter_draws

   function bssm_nonlinear_psi_filter(y, initial_mean, initial_covariance, &
      parameters, noise_dimension, observation, transition_model, particles, &
      approximation_iterations, convergence_tolerance) result(out)
      !! Run a nonlinear Gaussian psi filter using the shared random stream.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: parameters(:) !! Model parameter values.
      integer, intent(in) :: noise_dimension !! Noise dimension.
      integer, intent(in) :: particles !! Number of particles.
      integer, intent(in) :: approximation_iterations !! Number of approximation iterations.
      real(dp), intent(in) :: convergence_tolerance !! Convergence tolerance.
      procedure(bssm_gaussian_observation_t) :: observation !! Observed value or vector.
      procedure(bssm_nonlinear_transition_jacobian_t) :: transition_model !! Transition model callback procedure.
      type(bssm_particle_filter_t) :: out
      real(dp), allocatable :: initial_normals(:, :)
      real(dp), allocatable :: proposal_normals(:, :, :)
      real(dp), allocatable :: terminal_normals(:, :), uniforms(:, :)
      integer :: time

      if (particles < 1 .or. noise_dimension < 1 .or. size(y) < 1) then
         out%info = 1
         return
      end if
      allocate(initial_normals(size(initial_mean), particles))
      allocate(proposal_normals(size(initial_mean), particles, size(y)))
      allocate(terminal_normals(noise_dimension, particles))
      allocate(uniforms(particles, size(y)))
      call random_standard_normal_matrix(initial_normals)
      do time = 1, size(y)
         call random_standard_normal_matrix(proposal_normals(:, :, time))
      end do
      call random_standard_normal_matrix(terminal_normals)
      call random_number(uniforms)
      out = bssm_nonlinear_psi_filter_draws(y, initial_mean, &
         initial_covariance, parameters, noise_dimension, observation, &
         transition_model, approximation_iterations, convergence_tolerance, &
         initial_normals, proposal_normals, terminal_normals, uniforms)
   end function bssm_nonlinear_psi_filter

   pure function bssm_multivariate_nonlinear_psi_filter_draws(y, &
      initial_mean, initial_covariance, parameters, &
      observation_noise_dimension, state_noise_dimension, observation, &
      transition_model, approximation_iterations, convergence_tolerance, &
      initial_normals, proposal_normals, terminal_normals, &
      resampling_uniforms) result(out)
      !! Run a draw-driven psi filter for vector nonlinear Gaussian data.
      real(dp), intent(in) :: y(:, :) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: parameters(:) !! Model parameter values.
      integer, intent(in) :: observation_noise_dimension !! Observation noise dimension.
      integer, intent(in) :: state_noise_dimension !! State noise dimension.
      integer, intent(in) :: approximation_iterations !! Number of approximation iterations.
      real(dp), intent(in) :: convergence_tolerance !! Convergence tolerance.
      real(dp), intent(in) :: initial_normals(:, :) !! Initial normals.
      real(dp), intent(in) :: proposal_normals(:, :, :) !! Standard-normal proposal draws.
      real(dp), intent(in) :: terminal_normals(:, :) !! Terminal normals.
      real(dp), intent(in) :: resampling_uniforms(:, :) !! Resampling uniforms.
      procedure(bssm_multivariate_gaussian_observation_t) :: observation !! Observed value or vector.
      procedure(bssm_nonlinear_transition_jacobian_t) :: transition_model !! Transition model callback procedure.
      type(bssm_particle_filter_t) :: out
      type(bssm_multivariate_nonlinear_approximation_t) :: approximation
      real(dp), allocatable :: parent(:, :), probability(:), log_weight(:)
      real(dp), allocatable :: exact_observation_mean(:)
      real(dp), allocatable :: approximate_observation_mean(:)
      real(dp), allocatable :: observation_jacobian(:, :)
      real(dp), allocatable :: observation_loading(:, :)
      real(dp), allocatable :: observation_covariance(:, :)
      real(dp), allocatable :: exact_transition_mean(:)
      real(dp), allocatable :: approximate_transition_mean(:)
      real(dp), allocatable :: transition_jacobian(:, :)
      real(dp), allocatable :: state_loading(:, :), state_covariance(:, :)
      real(dp) :: exact_density, approximate_density, increment
      integer :: series, state, particles, times, particle, time, info

      series = size(y, 1)
      times = size(y, 2)
      state = size(initial_mean)
      particles = size(initial_normals, 2)
      if (series < 1 .or. state < 1 .or. particles < 1 .or. times < 1 .or. &
         observation_noise_dimension < 1 .or. state_noise_dimension < 1 .or. &
         any(shape(initial_covariance) /= [state, state]) .or. &
         size(initial_normals, 1) /= state .or. &
         any(shape(proposal_normals) /= [state, particles, times]) .or. &
         any(shape(terminal_normals) /= &
         [state_noise_dimension, particles]) .or. &
         any(shape(resampling_uniforms) /= [particles, times]) .or. &
         .not. all(ieee_is_finite(initial_normals)) .or. &
         .not. all(ieee_is_finite(proposal_normals)) .or. &
         .not. all(ieee_is_finite(terminal_normals)) .or. &
         .not. all(ieee_is_finite(resampling_uniforms))) then
         out%info = 1
         return
      end if
      approximation = bssm_multivariate_nonlinear_gaussian_approximation( &
         y, initial_mean, initial_covariance, parameters, &
         observation_noise_dimension, state_noise_dimension, observation, &
         transition_model, approximation_iterations, convergence_tolerance)
      if (approximation%info /= 0) then
         out%info = 1000 + approximation%info
         return
      end if
      allocate(out%particles(state, particles, times + 1))
      allocate(out%weights(particles, times), out%ancestors(particles, times))
      allocate(out%predicted_mean(state, times + 1))
      allocate(out%filtered_mean(state, times))
      allocate(out%predicted_covariance(state, state, times + 1))
      allocate(out%filtered_covariance(state, state, times))
      allocate(parent(state, particles), probability(particles))
      allocate(log_weight(particles), exact_observation_mean(series))
      allocate(approximate_observation_mean(series))
      allocate(observation_jacobian(series, state))
      allocate(observation_loading(series, observation_noise_dimension))
      allocate(observation_covariance(series, series))
      allocate(exact_transition_mean(state))
      allocate(approximate_transition_mean(state))
      allocate(transition_jacobian(state, state))
      allocate(state_loading(state, state_noise_dimension))
      allocate(state_covariance(state, state))
      do particle = 1, particles
         out%particles(:, particle, 1) = approximation%proposal_mean(:, 1) + &
            matmul(approximation%proposal_factor(:, :, 1), &
            initial_normals(:, particle))
      end do
      out%log_likelihood = approximation%corrected_log_likelihood
      parent = 0.0_dp
      do time = 1, times
         call particle_summary(out%particles(:, :, time), &
            out%predicted_mean(:, time), &
            out%predicted_covariance(:, :, time))
         do particle = 1, particles
            call observation(time, out%particles(:, particle, time), &
               parameters, exact_observation_mean, observation_jacobian, &
               observation_loading)
            if (.not. all(ieee_is_finite(exact_observation_mean)) .or. &
               .not. all(ieee_is_finite(observation_jacobian)) .or. &
               .not. all(ieee_is_finite(observation_loading))) then
               out%info = 2
               return
            end if
            observation_covariance = matmul(observation_loading, &
               transpose(observation_loading))
            exact_density = observed_density(y(:, time), &
               exact_observation_mean, observation_covariance)
            approximate_observation_mean = &
               approximation%observation_intercept(:, time) + &
               matmul(approximation%observation_jacobian(:, :, time), &
               out%particles(:, particle, time))
            observation_covariance = matmul( &
               approximation%observation_noise_loading(:, :, time), &
               transpose(approximation% &
               observation_noise_loading(:, :, time)))
            approximate_density = observed_density(y(:, time), &
               approximate_observation_mean, observation_covariance)
            log_weight(particle) = exact_density - approximate_density
            if (time > 1) then
               call transition_model(time - 1, parent(:, particle), &
                  parameters, exact_transition_mean, transition_jacobian, &
                  state_loading)
               if (.not. all(ieee_is_finite(exact_transition_mean)) .or. &
                  .not. all(ieee_is_finite(transition_jacobian)) .or. &
                  .not. all(ieee_is_finite(state_loading))) then
                  out%info = 3
                  return
               end if
               state_covariance = matmul(state_loading, &
                  transpose(state_loading))
               exact_density = multivariate_normal_log_density( &
                  out%particles(:, particle, time), &
                  exact_transition_mean, state_covariance)
               approximate_transition_mean = &
                  approximation%transition_intercept(:, time - 1) + &
                  matmul(approximation% &
                  transition_jacobian(:, :, time - 1), parent(:, particle))
               state_covariance = matmul( &
                  approximation%state_noise_loading(:, :, time - 1), &
                  transpose(approximation% &
                  state_noise_loading(:, :, time - 1)))
               approximate_density = multivariate_normal_log_density( &
                  out%particles(:, particle, time), &
                  approximate_transition_mean, state_covariance)
               log_weight(particle) = log_weight(particle) + &
                  exact_density - approximate_density
            end if
            log_weight(particle) = log_weight(particle) - &
               approximation%scaling(time)
         end do
         call normalize_log_weights(log_weight, probability, increment, info)
         if (info /= 0) then
            out%info = 4
            out%log_likelihood = -huge(1.0_dp)
            return
         end if
         out%log_likelihood = out%log_likelihood + increment
         out%weights(:, time) = probability
         call weighted_particle_summary(out%particles(:, :, time), &
            probability, out%filtered_mean(:, time), &
            out%filtered_covariance(:, :, time))
         out%ancestors(:, time) = bssm_stratified_resample(probability, &
            resampling_uniforms(:, time))
         do particle = 1, particles
            parent(:, particle) = out%particles(:, &
               out%ancestors(particle, time), time)
            if (time < times) then
               out%particles(:, particle, time + 1) = &
                  approximation%proposal_mean(:, time + 1) + &
                  matmul(approximation%conditional_matrix(:, :, time + 1), &
                  parent(:, particle) - &
                  approximation%proposal_mean(:, time)) + &
                  matmul(approximation%proposal_factor(:, :, time + 1), &
                  proposal_normals(:, particle, time))
            else
               call transition_model(time, parent(:, particle), parameters, &
                  exact_transition_mean, transition_jacobian, state_loading)
               if (.not. all(ieee_is_finite(exact_transition_mean)) .or. &
                  .not. all(ieee_is_finite(state_loading))) then
                  out%info = 3
                  return
               end if
               out%particles(:, particle, time + 1) = &
                  exact_transition_mean + &
                  matmul(state_loading, terminal_normals(:, particle))
            end if
         end do
      end do
      call particle_summary(out%particles(:, :, times + 1), &
         out%predicted_mean(:, times + 1), &
         out%predicted_covariance(:, :, times + 1))

   contains

      pure function observed_density(value, mean, covariance) result(value_log)
         !! Evaluate a Gaussian density over the finite vector components.
         real(dp), intent(in) :: value(:) !! Input value.
         real(dp), intent(in) :: mean(:) !! Mean value or vector.
         real(dp), intent(in) :: covariance(:, :) !! Covariance matrix.
         real(dp) :: value_log
         real(dp), allocatable :: compact_covariance(:, :)
         integer, allocatable :: observed(:)
         integer :: observed_count, component, row, column

         observed_count = count(ieee_is_finite(value))
         if (observed_count == 0) then
            value_log = 0.0_dp
            return
         end if
         allocate(observed(observed_count))
         observed = pack([(component, component = 1, size(value))], &
            ieee_is_finite(value))
         allocate(compact_covariance(observed_count, observed_count))
         do row = 1, observed_count
            do column = 1, observed_count
               compact_covariance(row, column) = &
                  covariance(observed(row), observed(column))
            end do
         end do
         value_log = multivariate_normal_log_density(value(observed), &
            mean(observed), compact_covariance)
      end function observed_density

   end function bssm_multivariate_nonlinear_psi_filter_draws

   function bssm_multivariate_nonlinear_psi_filter(y, initial_mean, &
      initial_covariance, parameters, observation_noise_dimension, &
      state_noise_dimension, observation, transition_model, particles, &
      approximation_iterations, convergence_tolerance) result(out)
      !! Run a vector nonlinear Gaussian psi filter with shared randomness.
      real(dp), intent(in) :: y(:, :) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: parameters(:) !! Model parameter values.
      integer, intent(in) :: observation_noise_dimension !! Observation noise dimension.
      integer, intent(in) :: state_noise_dimension !! State noise dimension.
      integer, intent(in) :: particles !! Number of particles.
      integer, intent(in) :: approximation_iterations !! Number of approximation iterations.
      real(dp), intent(in) :: convergence_tolerance !! Convergence tolerance.
      procedure(bssm_multivariate_gaussian_observation_t) :: observation !! Observed value or vector.
      procedure(bssm_nonlinear_transition_jacobian_t) :: transition_model !! Transition model callback procedure.
      type(bssm_particle_filter_t) :: out
      real(dp), allocatable :: initial_normals(:, :)
      real(dp), allocatable :: proposal_normals(:, :, :)
      real(dp), allocatable :: terminal_normals(:, :), uniforms(:, :)
      integer :: time

      if (particles < 1 .or. observation_noise_dimension < 1 .or. &
         state_noise_dimension < 1 .or. size(y, 1) < 1 .or. &
         size(y, 2) < 1) then
         out%info = 1
         return
      end if
      allocate(initial_normals(size(initial_mean), particles))
      allocate(proposal_normals(size(initial_mean), particles, size(y, 2)))
      allocate(terminal_normals(state_noise_dimension, particles))
      allocate(uniforms(particles, size(y, 2)))
      call random_standard_normal_matrix(initial_normals)
      do time = 1, size(y, 2)
         call random_standard_normal_matrix(proposal_normals(:, :, time))
      end do
      call random_standard_normal_matrix(terminal_normals)
      call random_number(uniforms)
      out = bssm_multivariate_nonlinear_psi_filter_draws(y, initial_mean, &
         initial_covariance, parameters, observation_noise_dimension, &
         state_noise_dimension, observation, transition_model, &
         approximation_iterations, convergence_tolerance, initial_normals, &
         proposal_normals, terminal_normals, uniforms)
   end function bssm_multivariate_nonlinear_psi_filter

   pure function bssm_nonlinear_psi_pmmh_draws(y, initial_mean, &
      initial_covariance, initial_parameters, noise_dimension, observation, &
      transition_model, prior, approximation_iterations, &
      convergence_tolerance, proposal_factor, parameter_normals, &
      acceptance_uniforms, initial_normals, state_normals, terminal_normals, &
      resampling_uniforms, target_acceptance, adaptation_exponent) result(out)
      !! Run draw-driven PMMH with a nonlinear psi particle filter.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: initial_parameters(:) !! Initial parameter values.
      real(dp), intent(in) :: convergence_tolerance !! Convergence tolerance.
      real(dp), intent(in) :: proposal_factor(:, :) !! Proposal factor.
      real(dp), intent(in) :: parameter_normals(:, :) !! Parameter normals.
      real(dp), intent(in) :: acceptance_uniforms(:) !! Acceptance uniforms.
      real(dp), intent(in) :: initial_normals(:, :, :) !! Initial normals.
      real(dp), intent(in) :: state_normals(:, :, :, :) !! State normals.
      real(dp), intent(in) :: terminal_normals(:, :, :) !! Terminal normals.
      real(dp), intent(in) :: resampling_uniforms(:, :, :) !! Resampling uniforms.
      integer, intent(in) :: noise_dimension !! Noise dimension.
      integer, intent(in) :: approximation_iterations !! Number of approximation iterations.
      procedure(bssm_gaussian_observation_t) :: observation !! Observed value or vector.
      procedure(bssm_nonlinear_transition_jacobian_t) :: transition_model !! Transition model callback procedure.
      procedure(bssm_parameter_log_density_t) :: prior !! Prior-distribution specification.
      real(dp), intent(in), optional :: target_acceptance !! Target acceptance.
      real(dp), intent(in), optional :: adaptation_exponent !! Adaptation exponent.
      type(bssm_mcmc_t) :: out
      integer :: state, particles, times, iterations
      real(dp) :: exponent

      state = size(initial_mean)
      particles = size(initial_normals, 2)
      times = size(y)
      iterations = size(parameter_normals, 2)
      if (state < 1 .or. particles < 1 .or. times < 1 .or. &
         noise_dimension < 1 .or. approximation_iterations < 1 .or. &
         convergence_tolerance < 0.0_dp .or. iterations < 1 .or. &
         any(shape(initial_covariance) /= [state, state]) .or. &
         any(shape(initial_normals) /= &
         [state, particles, iterations + 1]) .or. &
         any(shape(state_normals) /= &
         [state, particles, times, iterations + 1]) .or. &
         any(shape(terminal_normals) /= &
         [noise_dimension, particles, iterations + 1]) .or. &
         any(shape(resampling_uniforms) /= &
         [particles, times, iterations + 1]) .or. &
         .not. all(ieee_is_finite(initial_normals)) .or. &
         .not. all(ieee_is_finite(state_normals)) .or. &
         .not. all(ieee_is_finite(terminal_normals)) .or. &
         .not. all(ieee_is_finite(resampling_uniforms)) .or. &
         (present(adaptation_exponent) .and. &
         .not. present(target_acceptance))) then
         out%info = 1
         return
      end if
      exponent = 0.6_dp
      if (present(adaptation_exponent)) exponent = adaptation_exponent
      if (present(target_acceptance)) then
         out = bssm_pmmh_kernel_draws(initial_parameters, prior, estimator, &
            proposal_factor, parameter_normals, acceptance_uniforms, &
            target_acceptance, exponent)
      else
         out = bssm_pmmh_kernel_draws(initial_parameters, prior, estimator, &
            proposal_factor, parameter_normals, acceptance_uniforms)
      end if

   contains

      pure function estimator(parameters, draw_index) result(estimate)
         !! Estimate the nonlinear psi-filter log likelihood.
         real(dp), intent(in) :: parameters(:) !! Model parameter values.
         integer, intent(in) :: draw_index !! Index of draw.
         type(bssm_likelihood_estimate_t) :: estimate
         type(bssm_particle_filter_t) :: fit

         fit = bssm_nonlinear_psi_filter_draws(y, initial_mean, &
            initial_covariance, parameters, noise_dimension, observation, &
            transition_model, approximation_iterations, &
            convergence_tolerance, initial_normals(:, :, draw_index), &
            state_normals(:, :, :, draw_index), &
            terminal_normals(:, :, draw_index), &
            resampling_uniforms(:, :, draw_index))
         estimate%value = fit%log_likelihood
         estimate%info = fit%info
      end function estimator

   end function bssm_nonlinear_psi_pmmh_draws

   function bssm_nonlinear_psi_pmmh(y, initial_mean, initial_covariance, &
      initial_parameters, noise_dimension, observation, transition_model, &
      prior, approximation_iterations, convergence_tolerance, &
      proposal_factor, particles, iterations, target_acceptance, &
      adaptation_exponent) result(out)
      !! Run nonlinear psi-filter PMMH using the shared random stream.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: initial_parameters(:) !! Initial parameter values.
      real(dp), intent(in) :: convergence_tolerance !! Convergence tolerance.
      real(dp), intent(in) :: proposal_factor(:, :) !! Proposal factor.
      integer, intent(in) :: noise_dimension !! Noise dimension.
      integer, intent(in) :: approximation_iterations !! Number of approximation iterations.
      integer, intent(in) :: particles !! Number of particles.
      integer, intent(in) :: iterations !! Number of algorithm iterations.
      procedure(bssm_gaussian_observation_t) :: observation !! Observed value or vector.
      procedure(bssm_nonlinear_transition_jacobian_t) :: transition_model !! Transition model callback procedure.
      procedure(bssm_parameter_log_density_t) :: prior !! Prior-distribution specification.
      real(dp), intent(in), optional :: target_acceptance !! Target acceptance.
      real(dp), intent(in), optional :: adaptation_exponent !! Adaptation exponent.
      type(bssm_mcmc_t) :: out
      real(dp), allocatable :: parameter_normals(:, :), acceptance_uniforms(:)
      real(dp), allocatable :: initial_normals(:, :, :)
      real(dp), allocatable :: state_normals(:, :, :, :)
      real(dp), allocatable :: terminal_normals(:, :, :)
      real(dp), allocatable :: resampling_uniforms(:, :, :)
      real(dp) :: exponent
      integer :: draw_index, time

      if (particles < 1 .or. iterations < 1 .or. noise_dimension < 1 .or. &
         approximation_iterations < 1 .or. convergence_tolerance < 0.0_dp .or. &
         size(y) < 1 .or. size(initial_mean) < 1 .or. &
         size(initial_parameters) < 1 .or. &
         (present(adaptation_exponent) .and. &
         .not. present(target_acceptance))) then
         out%info = 1
         return
      end if
      allocate(parameter_normals(size(initial_parameters), iterations))
      allocate(acceptance_uniforms(iterations))
      allocate(initial_normals(size(initial_mean), particles, iterations + 1))
      allocate(state_normals(size(initial_mean), particles, size(y), &
         iterations + 1))
      allocate(terminal_normals(noise_dimension, particles, iterations + 1))
      allocate(resampling_uniforms(particles, size(y), iterations + 1))
      call random_standard_normal_matrix(parameter_normals)
      call random_number(acceptance_uniforms)
      call random_number(resampling_uniforms)
      do draw_index = 1, iterations + 1
         call random_standard_normal_matrix(initial_normals(:, :, draw_index))
         call random_standard_normal_matrix(terminal_normals(:, :, draw_index))
         do time = 1, size(y)
            call random_standard_normal_matrix( &
               state_normals(:, :, time, draw_index))
         end do
      end do
      exponent = 0.6_dp
      if (present(adaptation_exponent)) exponent = adaptation_exponent
      if (present(target_acceptance)) then
         out = bssm_nonlinear_psi_pmmh_draws(y, initial_mean, &
            initial_covariance, initial_parameters, noise_dimension, &
            observation, transition_model, prior, approximation_iterations, &
            convergence_tolerance, proposal_factor, parameter_normals, &
            acceptance_uniforms, initial_normals, state_normals, &
            terminal_normals, resampling_uniforms, target_acceptance, exponent)
      else
         out = bssm_nonlinear_psi_pmmh_draws(y, initial_mean, &
            initial_covariance, initial_parameters, noise_dimension, &
            observation, transition_model, prior, approximation_iterations, &
            convergence_tolerance, proposal_factor, parameter_normals, &
            acceptance_uniforms, initial_normals, state_normals, &
            terminal_normals, resampling_uniforms)
      end if
   end function bssm_nonlinear_psi_pmmh

   pure function bssm_multivariate_nonlinear_psi_pmmh_draws(y, initial_mean, &
      initial_covariance, initial_parameters, observation_noise_dimension, &
      state_noise_dimension, observation, transition_model, prior, &
      approximation_iterations, convergence_tolerance, proposal_factor, &
      parameter_normals, acceptance_uniforms, initial_normals, state_normals, &
      terminal_normals, resampling_uniforms, target_acceptance, &
      adaptation_exponent) result(out)
      !! Run draw-driven PMMH with a vector nonlinear psi particle filter.
      real(dp), intent(in) :: y(:, :) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: initial_parameters(:) !! Initial parameter values.
      real(dp), intent(in) :: convergence_tolerance !! Convergence tolerance.
      real(dp), intent(in) :: proposal_factor(:, :) !! Proposal factor.
      real(dp), intent(in) :: parameter_normals(:, :) !! Parameter normals.
      real(dp), intent(in) :: acceptance_uniforms(:) !! Acceptance uniforms.
      real(dp), intent(in) :: initial_normals(:, :, :) !! Initial normals.
      real(dp), intent(in) :: state_normals(:, :, :, :) !! State normals.
      real(dp), intent(in) :: terminal_normals(:, :, :) !! Terminal normals.
      real(dp), intent(in) :: resampling_uniforms(:, :, :) !! Resampling uniforms.
      integer, intent(in) :: observation_noise_dimension !! Observation noise dimension.
      integer, intent(in) :: state_noise_dimension !! State noise dimension.
      integer, intent(in) :: approximation_iterations !! Number of approximation iterations.
      procedure(bssm_multivariate_gaussian_observation_t) :: observation !! Observed value or vector.
      procedure(bssm_nonlinear_transition_jacobian_t) :: transition_model !! Transition model callback procedure.
      procedure(bssm_parameter_log_density_t) :: prior !! Prior-distribution specification.
      real(dp), intent(in), optional :: target_acceptance !! Target acceptance.
      real(dp), intent(in), optional :: adaptation_exponent !! Adaptation exponent.
      type(bssm_mcmc_t) :: out
      integer :: state, particles, times, iterations
      real(dp) :: exponent

      state = size(initial_mean)
      particles = size(initial_normals, 2)
      times = size(y, 2)
      iterations = size(parameter_normals, 2)
      if (size(y, 1) < 1 .or. state < 1 .or. particles < 1 .or. &
         times < 1 .or. iterations < 1 .or. &
         observation_noise_dimension < 1 .or. state_noise_dimension < 1 .or. &
         approximation_iterations < 1 .or. &
         convergence_tolerance < 0.0_dp .or. &
         any(shape(initial_covariance) /= [state, state]) .or. &
         any(shape(initial_normals) /= &
         [state, particles, iterations + 1]) .or. &
         any(shape(state_normals) /= &
         [state, particles, times, iterations + 1]) .or. &
         any(shape(terminal_normals) /= &
         [state_noise_dimension, particles, iterations + 1]) .or. &
         any(shape(resampling_uniforms) /= &
         [particles, times, iterations + 1]) .or. &
         .not. all(ieee_is_finite(initial_normals)) .or. &
         .not. all(ieee_is_finite(state_normals)) .or. &
         .not. all(ieee_is_finite(terminal_normals)) .or. &
         .not. all(ieee_is_finite(resampling_uniforms)) .or. &
         (present(adaptation_exponent) .and. &
         .not. present(target_acceptance))) then
         out%info = 1
         return
      end if
      exponent = 0.6_dp
      if (present(adaptation_exponent)) exponent = adaptation_exponent
      if (present(target_acceptance)) then
         out = bssm_pmmh_kernel_draws(initial_parameters, prior, estimator, &
            proposal_factor, parameter_normals, acceptance_uniforms, &
            target_acceptance, exponent)
      else
         out = bssm_pmmh_kernel_draws(initial_parameters, prior, estimator, &
            proposal_factor, parameter_normals, acceptance_uniforms)
      end if

   contains

      pure function estimator(parameters, draw_index) result(estimate)
         !! Estimate the vector nonlinear psi-filter log likelihood.
         real(dp), intent(in) :: parameters(:) !! Model parameter values.
         integer, intent(in) :: draw_index !! Index of draw.
         type(bssm_likelihood_estimate_t) :: estimate
         type(bssm_particle_filter_t) :: fit

         fit = bssm_multivariate_nonlinear_psi_filter_draws(y, initial_mean, &
            initial_covariance, parameters, observation_noise_dimension, &
            state_noise_dimension, observation, transition_model, &
            approximation_iterations, convergence_tolerance, &
            initial_normals(:, :, draw_index), &
            state_normals(:, :, :, draw_index), &
            terminal_normals(:, :, draw_index), &
            resampling_uniforms(:, :, draw_index))
         estimate%value = fit%log_likelihood
         estimate%info = fit%info
      end function estimator

   end function bssm_multivariate_nonlinear_psi_pmmh_draws

   function bssm_multivariate_nonlinear_psi_pmmh(y, initial_mean, &
      initial_covariance, initial_parameters, observation_noise_dimension, &
      state_noise_dimension, observation, transition_model, prior, &
      approximation_iterations, convergence_tolerance, proposal_factor, &
      particles, iterations, target_acceptance, adaptation_exponent) &
      result(out)
      !! Run vector nonlinear psi-filter PMMH using shared randomness.
      real(dp), intent(in) :: y(:, :) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: initial_parameters(:) !! Initial parameter values.
      real(dp), intent(in) :: convergence_tolerance !! Convergence tolerance.
      real(dp), intent(in) :: proposal_factor(:, :) !! Proposal factor.
      integer, intent(in) :: observation_noise_dimension !! Observation noise dimension.
      integer, intent(in) :: state_noise_dimension !! State noise dimension.
      integer, intent(in) :: approximation_iterations !! Number of approximation iterations.
      integer, intent(in) :: particles !! Number of particles.
      integer, intent(in) :: iterations !! Number of algorithm iterations.
      procedure(bssm_multivariate_gaussian_observation_t) :: observation !! Observed value or vector.
      procedure(bssm_nonlinear_transition_jacobian_t) :: transition_model !! Transition model callback procedure.
      procedure(bssm_parameter_log_density_t) :: prior !! Prior-distribution specification.
      real(dp), intent(in), optional :: target_acceptance !! Target acceptance.
      real(dp), intent(in), optional :: adaptation_exponent !! Adaptation exponent.
      type(bssm_mcmc_t) :: out
      real(dp), allocatable :: parameter_normals(:, :), acceptance_uniforms(:)
      real(dp), allocatable :: initial_normals(:, :, :)
      real(dp), allocatable :: state_normals(:, :, :, :)
      real(dp), allocatable :: terminal_normals(:, :, :)
      real(dp), allocatable :: resampling_uniforms(:, :, :)
      real(dp) :: exponent
      integer :: draw_index, time

      if (size(y, 1) < 1 .or. size(y, 2) < 1 .or. &
         size(initial_mean) < 1 .or. size(initial_parameters) < 1 .or. &
         observation_noise_dimension < 1 .or. state_noise_dimension < 1 .or. &
         approximation_iterations < 1 .or. &
         convergence_tolerance < 0.0_dp .or. &
         particles < 1 .or. iterations < 1 .or. &
         (present(adaptation_exponent) .and. &
         .not. present(target_acceptance))) then
         out%info = 1
         return
      end if
      allocate(parameter_normals(size(initial_parameters), iterations))
      allocate(acceptance_uniforms(iterations))
      allocate(initial_normals(size(initial_mean), particles, iterations + 1))
      allocate(state_normals(size(initial_mean), particles, size(y, 2), &
         iterations + 1))
      allocate(terminal_normals(state_noise_dimension, particles, &
         iterations + 1))
      allocate(resampling_uniforms(particles, size(y, 2), iterations + 1))
      call random_standard_normal_matrix(parameter_normals)
      call random_number(acceptance_uniforms)
      call random_number(resampling_uniforms)
      do draw_index = 1, iterations + 1
         call random_standard_normal_matrix(initial_normals(:, :, draw_index))
         call random_standard_normal_matrix(terminal_normals(:, :, draw_index))
         do time = 1, size(y, 2)
            call random_standard_normal_matrix( &
               state_normals(:, :, time, draw_index))
         end do
      end do
      exponent = 0.6_dp
      if (present(adaptation_exponent)) exponent = adaptation_exponent
      if (present(target_acceptance)) then
         out = bssm_multivariate_nonlinear_psi_pmmh_draws(y, initial_mean, &
            initial_covariance, initial_parameters, &
            observation_noise_dimension, state_noise_dimension, observation, &
            transition_model, prior, approximation_iterations, &
            convergence_tolerance, proposal_factor, parameter_normals, &
            acceptance_uniforms, initial_normals, state_normals, &
            terminal_normals, resampling_uniforms, target_acceptance, exponent)
      else
         out = bssm_multivariate_nonlinear_psi_pmmh_draws(y, initial_mean, &
            initial_covariance, initial_parameters, &
            observation_noise_dimension, state_noise_dimension, observation, &
            transition_model, prior, approximation_iterations, &
            convergence_tolerance, proposal_factor, parameter_normals, &
            acceptance_uniforms, initial_normals, state_normals, &
            terminal_normals, resampling_uniforms)
      end if
   end function bssm_multivariate_nonlinear_psi_pmmh

   pure function bssm_multivariate_nonlinear_psi_da_pmmh_draws(y, &
      initial_mean, initial_covariance, initial_parameters, &
      observation_noise_dimension, state_noise_dimension, observation, &
      transition_model, prior, approximation_iterations, &
      convergence_tolerance, proposal_factor, parameter_normals, &
      first_stage_uniforms, second_stage_uniforms, initial_normals, &
      state_normals, terminal_normals, resampling_uniforms, &
      target_acceptance, adaptation_exponent) result(out)
      !! Run draw-driven delayed acceptance with a vector nonlinear psi filter.
      real(dp), intent(in) :: y(:, :) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: initial_parameters(:) !! Initial parameter values.
      real(dp), intent(in) :: convergence_tolerance !! Convergence tolerance.
      real(dp), intent(in) :: proposal_factor(:, :) !! Proposal factor.
      real(dp), intent(in) :: parameter_normals(:, :) !! Parameter normals.
      real(dp), intent(in) :: first_stage_uniforms(:) !! First stage uniforms.
      real(dp), intent(in) :: second_stage_uniforms(:) !! Second stage uniforms.
      real(dp), intent(in) :: initial_normals(:, :, :) !! Initial normals.
      real(dp), intent(in) :: state_normals(:, :, :, :) !! State normals.
      real(dp), intent(in) :: terminal_normals(:, :, :) !! Terminal normals.
      real(dp), intent(in) :: resampling_uniforms(:, :, :) !! Resampling uniforms.
      integer, intent(in) :: observation_noise_dimension !! Observation noise dimension.
      integer, intent(in) :: state_noise_dimension !! State noise dimension.
      integer, intent(in) :: approximation_iterations !! Number of approximation iterations.
      procedure(bssm_multivariate_gaussian_observation_t) :: observation !! Observed value or vector.
      procedure(bssm_nonlinear_transition_jacobian_t) :: transition_model !! Transition model callback procedure.
      procedure(bssm_parameter_log_density_t) :: prior !! Prior-distribution specification.
      real(dp), intent(in), optional :: target_acceptance !! Target acceptance.
      real(dp), intent(in), optional :: adaptation_exponent !! Adaptation exponent.
      type(bssm_da_mcmc_t) :: out
      integer :: state, particles, times, iterations
      real(dp) :: exponent

      state = size(initial_mean)
      particles = size(initial_normals, 2)
      times = size(y, 2)
      iterations = size(parameter_normals, 2)
      if (size(y, 1) < 1 .or. state < 1 .or. particles < 1 .or. &
         times < 1 .or. iterations < 1 .or. &
         observation_noise_dimension < 1 .or. state_noise_dimension < 1 .or. &
         approximation_iterations < 1 .or. &
         convergence_tolerance < 0.0_dp .or. &
         any(shape(initial_covariance) /= [state, state]) .or. &
         any(shape(initial_normals) /= &
         [state, particles, iterations + 1]) .or. &
         any(shape(state_normals) /= &
         [state, particles, times, iterations + 1]) .or. &
         any(shape(terminal_normals) /= &
         [state_noise_dimension, particles, iterations + 1]) .or. &
         any(shape(resampling_uniforms) /= &
         [particles, times, iterations + 1]) .or. &
         .not. all(ieee_is_finite(initial_normals)) .or. &
         .not. all(ieee_is_finite(state_normals)) .or. &
         .not. all(ieee_is_finite(terminal_normals)) .or. &
         .not. all(ieee_is_finite(resampling_uniforms)) .or. &
         (present(adaptation_exponent) .and. &
         .not. present(target_acceptance))) then
         out%info = 1
         return
      end if
      exponent = 0.6_dp
      if (present(adaptation_exponent)) exponent = adaptation_exponent
      if (present(target_acceptance)) then
         out = bssm_da_pmmh_kernel_draws(initial_parameters, prior, &
            coarse_estimator, fine_estimator, proposal_factor, &
            parameter_normals, first_stage_uniforms, second_stage_uniforms, &
            target_acceptance, exponent)
      else
         out = bssm_da_pmmh_kernel_draws(initial_parameters, prior, &
            coarse_estimator, fine_estimator, proposal_factor, &
            parameter_normals, first_stage_uniforms, second_stage_uniforms)
      end if

   contains

      pure function coarse_estimator(parameters, draw_index) result(estimate)
         !! Evaluate the global Gaussian-approximation likelihood.
         real(dp), intent(in) :: parameters(:) !! Model parameter values.
         integer, intent(in) :: draw_index !! Index of draw.
         type(bssm_likelihood_estimate_t) :: estimate
         type(bssm_multivariate_nonlinear_approximation_t) :: approximation

         approximation = bssm_multivariate_nonlinear_gaussian_approximation( &
            y, initial_mean, initial_covariance, parameters, &
            observation_noise_dimension, state_noise_dimension, observation, &
            transition_model, approximation_iterations, &
            convergence_tolerance)
         estimate%value = approximation%corrected_log_likelihood + &
            0.0_dp*real(draw_index, dp)
         estimate%info = approximation%info
      end function coarse_estimator

      pure function fine_estimator(parameters, draw_index) result(estimate)
         !! Evaluate the vector nonlinear psi particle likelihood.
         real(dp), intent(in) :: parameters(:) !! Model parameter values.
         integer, intent(in) :: draw_index !! Index of draw.
         type(bssm_likelihood_estimate_t) :: estimate
         type(bssm_particle_filter_t) :: fit

         fit = bssm_multivariate_nonlinear_psi_filter_draws(y, initial_mean, &
            initial_covariance, parameters, observation_noise_dimension, &
            state_noise_dimension, observation, transition_model, &
            approximation_iterations, convergence_tolerance, &
            initial_normals(:, :, draw_index), &
            state_normals(:, :, :, draw_index), &
            terminal_normals(:, :, draw_index), &
            resampling_uniforms(:, :, draw_index))
         estimate%value = fit%log_likelihood
         estimate%info = fit%info
      end function fine_estimator

   end function bssm_multivariate_nonlinear_psi_da_pmmh_draws

   function bssm_multivariate_nonlinear_psi_da_pmmh(y, initial_mean, &
      initial_covariance, initial_parameters, observation_noise_dimension, &
      state_noise_dimension, observation, transition_model, prior, &
      approximation_iterations, convergence_tolerance, proposal_factor, &
      particles, iterations, target_acceptance, adaptation_exponent) &
      result(out)
      !! Run delayed-acceptance vector nonlinear psi PMMH with shared randomness.
      real(dp), intent(in) :: y(:, :) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: initial_parameters(:) !! Initial parameter values.
      real(dp), intent(in) :: convergence_tolerance !! Convergence tolerance.
      real(dp), intent(in) :: proposal_factor(:, :) !! Proposal factor.
      integer, intent(in) :: observation_noise_dimension !! Observation noise dimension.
      integer, intent(in) :: state_noise_dimension !! State noise dimension.
      integer, intent(in) :: approximation_iterations !! Number of approximation iterations.
      integer, intent(in) :: particles !! Number of particles.
      integer, intent(in) :: iterations !! Number of algorithm iterations.
      procedure(bssm_multivariate_gaussian_observation_t) :: observation !! Observed value or vector.
      procedure(bssm_nonlinear_transition_jacobian_t) :: transition_model !! Transition model callback procedure.
      procedure(bssm_parameter_log_density_t) :: prior !! Prior-distribution specification.
      real(dp), intent(in), optional :: target_acceptance !! Target acceptance.
      real(dp), intent(in), optional :: adaptation_exponent !! Adaptation exponent.
      type(bssm_da_mcmc_t) :: out
      real(dp), allocatable :: parameter_normals(:, :)
      real(dp), allocatable :: first_stage_uniforms(:), second_stage_uniforms(:)
      real(dp), allocatable :: initial_normals(:, :, :)
      real(dp), allocatable :: state_normals(:, :, :, :)
      real(dp), allocatable :: terminal_normals(:, :, :)
      real(dp), allocatable :: resampling_uniforms(:, :, :)
      real(dp) :: exponent
      integer :: draw_index, time

      if (size(y, 1) < 1 .or. size(y, 2) < 1 .or. &
         size(initial_mean) < 1 .or. size(initial_parameters) < 1 .or. &
         observation_noise_dimension < 1 .or. state_noise_dimension < 1 .or. &
         approximation_iterations < 1 .or. &
         convergence_tolerance < 0.0_dp .or. &
         particles < 1 .or. iterations < 1 .or. &
         (present(adaptation_exponent) .and. &
         .not. present(target_acceptance))) then
         out%info = 1
         return
      end if
      allocate(parameter_normals(size(initial_parameters), iterations))
      allocate(first_stage_uniforms(iterations), second_stage_uniforms(iterations))
      allocate(initial_normals(size(initial_mean), particles, iterations + 1))
      allocate(state_normals(size(initial_mean), particles, size(y, 2), &
         iterations + 1))
      allocate(terminal_normals(state_noise_dimension, particles, &
         iterations + 1))
      allocate(resampling_uniforms(particles, size(y, 2), iterations + 1))
      call random_standard_normal_matrix(parameter_normals)
      call random_number(first_stage_uniforms)
      call random_number(second_stage_uniforms)
      call random_number(resampling_uniforms)
      do draw_index = 1, iterations + 1
         call random_standard_normal_matrix(initial_normals(:, :, draw_index))
         call random_standard_normal_matrix(terminal_normals(:, :, draw_index))
         do time = 1, size(y, 2)
            call random_standard_normal_matrix( &
               state_normals(:, :, time, draw_index))
         end do
      end do
      exponent = 0.6_dp
      if (present(adaptation_exponent)) exponent = adaptation_exponent
      if (present(target_acceptance)) then
         out = bssm_multivariate_nonlinear_psi_da_pmmh_draws(y, initial_mean, &
            initial_covariance, initial_parameters, &
            observation_noise_dimension, state_noise_dimension, observation, &
            transition_model, prior, approximation_iterations, &
            convergence_tolerance, proposal_factor, parameter_normals, &
            first_stage_uniforms, second_stage_uniforms, initial_normals, &
            state_normals, terminal_normals, resampling_uniforms, &
            target_acceptance, exponent)
      else
         out = bssm_multivariate_nonlinear_psi_da_pmmh_draws(y, initial_mean, &
            initial_covariance, initial_parameters, &
            observation_noise_dimension, state_noise_dimension, observation, &
            transition_model, prior, approximation_iterations, &
            convergence_tolerance, proposal_factor, parameter_normals, &
            first_stage_uniforms, second_stage_uniforms, initial_normals, &
            state_normals, terminal_normals, resampling_uniforms)
      end if
   end function bssm_multivariate_nonlinear_psi_da_pmmh

   pure function bssm_iekf(y, initial_mean, initial_covariance, parameters, &
      noise_dimension, observation, transition, max_iterations, &
      convergence_tolerance) result(out)
      !! Run an ordinary or iterated extended Kalman filter.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: parameters(:) !! Model parameter values.
      integer, intent(in) :: noise_dimension !! Noise dimension.
      integer, intent(in) :: max_iterations !! Maximum number of algorithm iterations.
      real(dp), intent(in) :: convergence_tolerance !! Convergence tolerance.
      procedure(bssm_gaussian_observation_t) :: observation !! Observed value or vector.
      procedure(bssm_nonlinear_transition_jacobian_t) :: transition !! State transition matrix.
      type(bssm_ekf_t) :: out
      real(dp), allocatable :: transition_jacobian(:, :), noise_loading(:, :)
      real(dp) :: innovation, innovation_variance
      integer :: state, times, time, info

      state = size(initial_mean)
      times = size(y)
      if (state < 1 .or. times < 1 .or. noise_dimension < 1 .or. &
         max_iterations < 0 .or. convergence_tolerance < 0.0_dp .or. &
         any(shape(initial_covariance) /= [state, state]) .or. &
         .not. all(ieee_is_finite(initial_mean)) .or. &
         .not. all(ieee_is_finite(initial_covariance)) .or. &
         .not. all(ieee_is_finite(parameters))) then
         out%info = 1
         return
      end if
      allocate(out%predicted_mean(state, times + 1))
      allocate(out%filtered_mean(state, times))
      allocate(out%predicted_covariance(state, state, times + 1))
      allocate(out%filtered_covariance(state, state, times))
      allocate(out%innovation(times), out%innovation_variance(times))
      allocate(out%iterations(times))
      allocate(transition_jacobian(state, state))
      allocate(noise_loading(state, noise_dimension))
      out%predicted_mean(:, 1) = initial_mean
      out%predicted_covariance(:, :, 1) = initial_covariance
      do time = 1, times
         if (ieee_is_finite(y(time))) then
            call ekf_proposal_update(time, y(time), &
               out%predicted_mean(:, time), &
               out%predicted_covariance(:, :, time), parameters, observation, &
               out%filtered_mean(:, time), &
               out%filtered_covariance(:, :, time), info, max_iterations, &
               convergence_tolerance, out%iterations(time), innovation, &
               innovation_variance)
            if (info /= 0) then
               out%info = 2
               out%log_likelihood = -huge(1.0_dp)
               return
            end if
            out%innovation(time) = innovation
            out%innovation_variance(time) = innovation_variance
            out%log_likelihood = out%log_likelihood - 0.5_dp* &
               (log(2.0_dp*acos(-1.0_dp)) + log(innovation_variance) + &
               innovation**2/innovation_variance)
         else
            out%filtered_mean(:, time) = out%predicted_mean(:, time)
            out%filtered_covariance(:, :, time) = &
               out%predicted_covariance(:, :, time)
            out%innovation(time) = 0.0_dp
            out%innovation_variance(time) = 0.0_dp
            out%iterations(time) = 0
         end if
         call transition(time, out%filtered_mean(:, time), parameters, &
            out%predicted_mean(:, time + 1), transition_jacobian, &
            noise_loading)
         if (.not. all(ieee_is_finite(out%predicted_mean(:, time + 1))) .or. &
            .not. all(ieee_is_finite(transition_jacobian)) .or. &
            .not. all(ieee_is_finite(noise_loading))) then
            out%info = 3
            return
         end if
         out%predicted_covariance(:, :, time + 1) = &
            matmul(matmul(transition_jacobian, &
            out%filtered_covariance(:, :, time)), &
            transpose(transition_jacobian)) + &
            matmul(noise_loading, transpose(noise_loading))
         out%predicted_covariance(:, :, time + 1) = &
            0.5_dp*(out%predicted_covariance(:, :, time + 1) + &
            transpose(out%predicted_covariance(:, :, time + 1)))
      end do
   end function bssm_iekf

   pure function bssm_ekf_smoother(y, initial_mean, initial_covariance, &
      parameters, noise_dimension, observation, transition, max_iterations, &
      convergence_tolerance) result(out)
      !! Run an ordinary or iterated extended Kalman covariance smoother.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: parameters(:) !! Model parameter values.
      integer, intent(in) :: noise_dimension !! Noise dimension.
      integer, intent(in) :: max_iterations !! Maximum number of algorithm iterations.
      real(dp), intent(in) :: convergence_tolerance !! Convergence tolerance.
      procedure(bssm_gaussian_observation_t) :: observation !! Observed value or vector.
      procedure(bssm_nonlinear_transition_jacobian_t) :: transition !! State transition matrix.
      type(bssm_ekf_smoother_t) :: out
      type(bssm_ekf_t) :: fit
      type(ssm_model_t) :: model
      type(kfs_filter_t) :: filtered
      type(kfs_smoother_t) :: smoothed
      integer :: info

      fit = bssm_iekf(y, initial_mean, initial_covariance, parameters, &
         noise_dimension, observation, transition, max_iterations, &
         convergence_tolerance)
      out%log_likelihood = fit%log_likelihood
      if (fit%info /= 0) then
         out%info = fit%info
         return
      end if
      out%iterations = fit%iterations
      call bssm_ekf_smoothing_inputs(fit, parameters, noise_dimension, &
         transition, model, filtered, info)
      if (info /= 0) then
         out%info = 100 + info
         return
      end if
      smoothed = kfs_smooth(model, filtered)
      if (smoothed%info /= 0) then
         out%info = 200 + smoothed%info
         return
      end if
      out%state = smoothed%state
      out%covariance = smoothed%covariance
   end function bssm_ekf_smoother

   pure function bssm_ekf_fast_smoother(y, initial_mean, &
      initial_covariance, parameters, noise_dimension, observation, &
      transition, max_iterations, convergence_tolerance) result(out)
      !! Run an extended Kalman smoother without smoothed covariances.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: parameters(:) !! Model parameter values.
      integer, intent(in) :: noise_dimension !! Noise dimension.
      integer, intent(in) :: max_iterations !! Maximum number of algorithm iterations.
      real(dp), intent(in) :: convergence_tolerance !! Convergence tolerance.
      procedure(bssm_gaussian_observation_t) :: observation !! Observed value or vector.
      procedure(bssm_nonlinear_transition_jacobian_t) :: transition !! State transition matrix.
      type(bssm_ekf_smoother_t) :: out
      type(bssm_ekf_t) :: fit
      type(ssm_model_t) :: model
      type(kfs_filter_t) :: filtered
      type(kfs_smoother_t) :: smoothed
      integer :: info

      fit = bssm_iekf(y, initial_mean, initial_covariance, parameters, &
         noise_dimension, observation, transition, max_iterations, &
         convergence_tolerance)
      out%log_likelihood = fit%log_likelihood
      if (fit%info /= 0) then
         out%info = fit%info
         return
      end if
      out%iterations = fit%iterations
      call bssm_ekf_smoothing_inputs(fit, parameters, noise_dimension, &
         transition, model, filtered, info)
      if (info /= 0) then
         out%info = 100 + info
         return
      end if
      smoothed = kfs_fast_smooth(model, filtered)
      if (smoothed%info /= 0) then
         out%info = 200 + smoothed%info
         return
      end if
      out%state = smoothed%state
   end function bssm_ekf_fast_smoother

   pure function bssm_multivariate_iekf(y, initial_mean, initial_covariance, &
      parameters, observation_noise_dimension, state_noise_dimension, &
      observation, transition, max_iterations, convergence_tolerance) &
      result(out)
      !! Run a multivariate ordinary or iterated extended Kalman filter.
      real(dp), intent(in) :: y(:, :) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: parameters(:) !! Model parameter values.
      integer, intent(in) :: observation_noise_dimension !! Observation noise dimension.
      integer, intent(in) :: state_noise_dimension !! State noise dimension.
      integer, intent(in) :: max_iterations !! Maximum number of algorithm iterations.
      real(dp), intent(in) :: convergence_tolerance !! Convergence tolerance.
      procedure(bssm_multivariate_gaussian_observation_t) :: observation !! Observed value or vector.
      procedure(bssm_nonlinear_transition_jacobian_t) :: transition !! State transition matrix.
      type(bssm_multivariate_ekf_t) :: out
      real(dp), allocatable :: transition_jacobian(:, :), noise_loading(:, :)
      real(dp) :: log_density
      integer :: series, state, times, time, info

      series = size(y, 1)
      times = size(y, 2)
      state = size(initial_mean)
      if (series < 1 .or. state < 1 .or. times < 1 .or. &
         observation_noise_dimension < 1 .or. state_noise_dimension < 1 .or. &
         max_iterations < 0 .or. convergence_tolerance < 0.0_dp .or. &
         any(shape(initial_covariance) /= [state, state]) .or. &
         .not. all(ieee_is_finite(initial_mean)) .or. &
         .not. all(ieee_is_finite(initial_covariance)) .or. &
         .not. all(ieee_is_finite(parameters))) then
         out%info = 1
         return
      end if
      allocate(out%predicted_mean(state, times + 1))
      allocate(out%filtered_mean(state, times))
      allocate(out%predicted_covariance(state, state, times + 1))
      allocate(out%filtered_covariance(state, state, times))
      allocate(out%innovation(series, times))
      allocate(out%innovation_covariance(series, series, times))
      allocate(out%iterations(times), transition_jacobian(state, state))
      allocate(noise_loading(state, state_noise_dimension))
      out%predicted_mean(:, 1) = initial_mean
      out%predicted_covariance(:, :, 1) = initial_covariance
      out%innovation = 0.0_dp
      out%innovation_covariance = 0.0_dp
      do time = 1, times
         if (any(ieee_is_finite(y(:, time)))) then
            call multivariate_ekf_proposal_update(time, y(:, time), &
               out%predicted_mean(:, time), &
               out%predicted_covariance(:, :, time), parameters, &
               observation_noise_dimension, observation, &
               out%filtered_mean(:, time), &
               out%filtered_covariance(:, :, time), out%innovation(:, time), &
               out%innovation_covariance(:, :, time), log_density, info, &
               max_iterations, convergence_tolerance, out%iterations(time))
            if (info /= 0) then
               out%info = 10 + info
               out%log_likelihood = -huge(1.0_dp)
               return
            end if
            out%log_likelihood = out%log_likelihood + log_density
         else
            out%filtered_mean(:, time) = out%predicted_mean(:, time)
            out%filtered_covariance(:, :, time) = &
               out%predicted_covariance(:, :, time)
            out%iterations(time) = 0
         end if
         call transition(time, out%filtered_mean(:, time), parameters, &
            out%predicted_mean(:, time + 1), transition_jacobian, &
            noise_loading)
         if (.not. all(ieee_is_finite(out%predicted_mean(:, time + 1))) .or. &
            .not. all(ieee_is_finite(transition_jacobian)) .or. &
            .not. all(ieee_is_finite(noise_loading))) then
            out%info = 3
            return
         end if
         out%predicted_covariance(:, :, time + 1) = &
            matmul(matmul(transition_jacobian, &
            out%filtered_covariance(:, :, time)), &
            transpose(transition_jacobian)) + &
            matmul(noise_loading, transpose(noise_loading))
         out%predicted_covariance(:, :, time + 1) = &
            0.5_dp*(out%predicted_covariance(:, :, time + 1) + &
            transpose(out%predicted_covariance(:, :, time + 1)))
      end do
   end function bssm_multivariate_iekf

   pure function bssm_multivariate_ekf_smoother(y, initial_mean, &
      initial_covariance, parameters, observation_noise_dimension, &
      state_noise_dimension, observation, transition, max_iterations, &
      convergence_tolerance) result(out)
      !! Smooth states from a multivariate ordinary or iterated EKF.
      real(dp), intent(in) :: y(:, :) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: parameters(:) !! Model parameter values.
      integer, intent(in) :: observation_noise_dimension !! Observation noise dimension.
      integer, intent(in) :: state_noise_dimension !! State noise dimension.
      integer, intent(in) :: max_iterations !! Maximum number of algorithm iterations.
      real(dp), intent(in) :: convergence_tolerance !! Convergence tolerance.
      procedure(bssm_multivariate_gaussian_observation_t) :: observation !! Observed value or vector.
      procedure(bssm_nonlinear_transition_jacobian_t) :: transition !! State transition matrix.
      type(bssm_ekf_smoother_t) :: out
      type(bssm_multivariate_ekf_t) :: fit
      type(ssm_model_t) :: model
      type(kfs_filter_t) :: filtered
      type(kfs_smoother_t) :: smoothed
      integer :: info

      fit = bssm_multivariate_iekf(y, initial_mean, initial_covariance, &
         parameters, observation_noise_dimension, state_noise_dimension, &
         observation, transition, max_iterations, convergence_tolerance)
      out%log_likelihood = fit%log_likelihood
      if (fit%info /= 0) then
         out%info = fit%info
         return
      end if
      out%iterations = fit%iterations
      call bssm_multivariate_ekf_smoothing_inputs(fit, parameters, &
         state_noise_dimension, transition, model, filtered, info)
      if (info /= 0) then
         out%info = 100 + info
         return
      end if
      smoothed = kfs_smooth(model, filtered)
      if (smoothed%info /= 0) then
         out%info = 200 + smoothed%info
         return
      end if
      out%state = smoothed%state
      out%covariance = smoothed%covariance
   end function bssm_multivariate_ekf_smoother

   pure function bssm_multivariate_ekf_fast_smoother(y, initial_mean, &
      initial_covariance, parameters, observation_noise_dimension, &
      state_noise_dimension, observation, transition, max_iterations, &
      convergence_tolerance) result(out)
      !! Smooth multivariate EKF state means without covariance output.
      real(dp), intent(in) :: y(:, :) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: parameters(:) !! Model parameter values.
      integer, intent(in) :: observation_noise_dimension !! Observation noise dimension.
      integer, intent(in) :: state_noise_dimension !! State noise dimension.
      integer, intent(in) :: max_iterations !! Maximum number of algorithm iterations.
      real(dp), intent(in) :: convergence_tolerance !! Convergence tolerance.
      procedure(bssm_multivariate_gaussian_observation_t) :: observation !! Observed value or vector.
      procedure(bssm_nonlinear_transition_jacobian_t) :: transition !! State transition matrix.
      type(bssm_ekf_smoother_t) :: out
      type(bssm_multivariate_ekf_t) :: fit
      type(ssm_model_t) :: model
      type(kfs_filter_t) :: filtered
      type(kfs_smoother_t) :: smoothed
      integer :: info

      fit = bssm_multivariate_iekf(y, initial_mean, initial_covariance, &
         parameters, observation_noise_dimension, state_noise_dimension, &
         observation, transition, max_iterations, convergence_tolerance)
      out%log_likelihood = fit%log_likelihood
      if (fit%info /= 0) then
         out%info = fit%info
         return
      end if
      out%iterations = fit%iterations
      call bssm_multivariate_ekf_smoothing_inputs(fit, parameters, &
         state_noise_dimension, transition, model, filtered, info)
      if (info /= 0) then
         out%info = 100 + info
         return
      end if
      smoothed = kfs_fast_smooth(model, filtered)
      if (smoothed%info /= 0) then
         out%info = 200 + smoothed%info
         return
      end if
      out%state = smoothed%state
   end function bssm_multivariate_ekf_fast_smoother

   pure function bssm_ukf(y, initial_mean, initial_covariance, parameters, &
      noise_dimension, observation, transition, alpha, beta, kappa) result(out)
      !! Run an unscented Kalman filter for a nonlinear Gaussian model.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: parameters(:) !! Model parameter values.
      real(dp), intent(in) :: alpha !! Significance, smoothing, or model coefficient.
      real(dp), intent(in) :: beta !! Regression or model coefficients.
      real(dp), intent(in) :: kappa !! Kappa.
      integer, intent(in) :: noise_dimension !! Noise dimension.
      procedure(bssm_gaussian_observation_t) :: observation !! Observed value or vector.
      procedure(bssm_nonlinear_transition_jacobian_t) :: transition !! State transition matrix.
      type(bssm_ekf_t) :: out
      real(dp), allocatable :: mean_weight(:), covariance_weight(:)
      real(dp), allocatable :: sigma(:, :), propagated(:, :), sigma_y(:)
      real(dp), allocatable :: factor(:, :), cross_covariance(:), gain(:)
      real(dp), allocatable :: observation_jacobian(:)
      real(dp), allocatable :: transition_jacobian(:, :), noise_loading(:, :)
      real(dp) :: lambda, sigma_scale, predicted_observation
      real(dp) :: observation_sd, difference, innovation_variance
      integer :: state, times, sigma_count, point, time, info

      state = size(initial_mean)
      times = size(y)
      lambda = alpha*alpha*(real(state, dp) + kappa) - real(state, dp)
      if (state < 1 .or. times < 1 .or. noise_dimension < 1 .or. &
         alpha <= 0.0_dp .or. beta < 0.0_dp .or. kappa < 0.0_dp .or. &
         real(state, dp) + lambda <= 0.0_dp .or. &
         any(shape(initial_covariance) /= [state, state]) .or. &
         .not. all(ieee_is_finite(initial_mean)) .or. &
         .not. all(ieee_is_finite(initial_covariance)) .or. &
         .not. all(ieee_is_finite(parameters))) then
         out%info = 1
         return
      end if
      sigma_count = 2*state + 1
      sigma_scale = sqrt(real(state, dp) + lambda)
      allocate(mean_weight(sigma_count), covariance_weight(sigma_count))
      mean_weight(1) = lambda/(real(state, dp) + lambda)
      mean_weight(2:) = 1.0_dp/(2.0_dp*(real(state, dp) + lambda))
      covariance_weight = mean_weight
      covariance_weight(1) = covariance_weight(1) + &
         1.0_dp - alpha*alpha + beta
      allocate(out%predicted_mean(state, times + 1))
      allocate(out%filtered_mean(state, times))
      allocate(out%predicted_covariance(state, state, times + 1))
      allocate(out%filtered_covariance(state, state, times))
      allocate(out%innovation(times), out%innovation_variance(times))
      allocate(out%iterations(times))
      allocate(sigma(state, sigma_count), propagated(state, sigma_count))
      allocate(sigma_y(sigma_count), cross_covariance(state), gain(state))
      allocate(observation_jacobian(state), transition_jacobian(state, state))
      allocate(noise_loading(state, noise_dimension))
      out%predicted_mean(:, 1) = initial_mean
      out%predicted_covariance(:, :, 1) = initial_covariance
      out%iterations = 0
      do time = 1, times
         call positive_semidefinite_factor( &
            out%predicted_covariance(:, :, time), factor, info)
         if (info /= 0) then
            out%info = 2
            return
         end if
         call form_sigma_points(out%predicted_mean(:, time), factor, &
            sigma_scale, sigma)
         if (ieee_is_finite(y(time))) then
            do point = 1, sigma_count
               call observation(time, sigma(:, point), parameters, &
                  sigma_y(point), observation_jacobian, observation_sd)
               if (.not. ieee_is_finite(sigma_y(point))) then
                  out%info = 3
                  return
               end if
            end do
            call observation(time, out%predicted_mean(:, time), parameters, &
               predicted_observation, observation_jacobian, observation_sd)
            if (.not. ieee_is_finite(observation_sd) .or. &
               observation_sd <= 0.0_dp) then
               out%info = 3
               return
            end if
            predicted_observation = dot_product(mean_weight, sigma_y)
            innovation_variance = observation_sd**2
            cross_covariance = 0.0_dp
            do point = 1, sigma_count
               difference = sigma_y(point) - predicted_observation
               innovation_variance = innovation_variance + &
                  covariance_weight(point)*difference**2
               cross_covariance = cross_covariance + &
                  covariance_weight(point)*(sigma(:, point) - &
                  out%predicted_mean(:, time))*difference
            end do
            if (.not. ieee_is_finite(innovation_variance) .or. &
               innovation_variance <= 0.0_dp) then
               out%info = 4
               return
            end if
            gain = cross_covariance/innovation_variance
            out%innovation(time) = y(time) - predicted_observation
            out%innovation_variance(time) = innovation_variance
            out%filtered_mean(:, time) = out%predicted_mean(:, time) + &
               gain*out%innovation(time)
            out%filtered_covariance(:, :, time) = &
               out%predicted_covariance(:, :, time) - &
               innovation_variance*outer_product(gain, gain)
            out%filtered_covariance(:, :, time) = &
               0.5_dp*(out%filtered_covariance(:, :, time) + &
               transpose(out%filtered_covariance(:, :, time)))
            out%log_likelihood = out%log_likelihood - 0.5_dp* &
               (log(2.0_dp*acos(-1.0_dp)) + log(innovation_variance) + &
               out%innovation(time)**2/innovation_variance)
         else
            out%filtered_mean(:, time) = out%predicted_mean(:, time)
            out%filtered_covariance(:, :, time) = &
               out%predicted_covariance(:, :, time)
            out%innovation(time) = 0.0_dp
            out%innovation_variance(time) = 0.0_dp
         end if
         call positive_semidefinite_factor( &
            out%filtered_covariance(:, :, time), factor, info)
         if (info /= 0) then
            out%info = 5
            return
         end if
         call form_sigma_points(out%filtered_mean(:, time), factor, &
            sigma_scale, sigma)
         do point = 1, sigma_count
            call transition(time, sigma(:, point), parameters, &
               propagated(:, point), transition_jacobian, noise_loading)
            if (.not. all(ieee_is_finite(propagated(:, point)))) then
               out%info = 6
               return
            end if
         end do
         out%predicted_mean(:, time + 1) = &
            matmul(propagated, mean_weight)
         call transition(time, out%filtered_mean(:, time), parameters, &
            sigma(:, 1), transition_jacobian, noise_loading)
         if (.not. all(ieee_is_finite(noise_loading))) then
            out%info = 6
            return
         end if
         out%predicted_covariance(:, :, time + 1) = &
            matmul(noise_loading, transpose(noise_loading))
         do point = 1, sigma_count
            out%predicted_covariance(:, :, time + 1) = &
               out%predicted_covariance(:, :, time + 1) + &
               covariance_weight(point)*outer_product( &
               propagated(:, point) - out%predicted_mean(:, time + 1), &
               propagated(:, point) - out%predicted_mean(:, time + 1))
         end do
         out%predicted_covariance(:, :, time + 1) = &
            0.5_dp*(out%predicted_covariance(:, :, time + 1) + &
            transpose(out%predicted_covariance(:, :, time + 1)))
      end do
   end function bssm_ukf

   pure function bssm_multivariate_ukf(y, initial_mean, initial_covariance, &
      parameters, observation_noise_dimension, state_noise_dimension, &
      observation, transition, alpha, beta, kappa) result(out)
      !! Run an unscented Kalman filter with vector Gaussian observations.
      real(dp), intent(in) :: y(:, :) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: parameters(:) !! Model parameter values.
      integer, intent(in) :: observation_noise_dimension !! Observation noise dimension.
      integer, intent(in) :: state_noise_dimension !! State noise dimension.
      procedure(bssm_multivariate_gaussian_observation_t) :: observation !! Observed value or vector.
      procedure(bssm_nonlinear_transition_jacobian_t) :: transition !! State transition matrix.
      real(dp), intent(in) :: alpha !! Significance, smoothing, or model coefficient.
      real(dp), intent(in) :: beta !! Regression or model coefficients.
      real(dp), intent(in) :: kappa !! Kappa.
      type(bssm_multivariate_ekf_t) :: out
      real(dp), allocatable :: mean_weight(:), covariance_weight(:)
      real(dp), allocatable :: sigma(:, :), propagated(:, :)
      real(dp), allocatable :: sigma_y(:, :), predicted_observation(:)
      real(dp), allocatable :: observation_jacobian(:, :)
      real(dp), allocatable :: observation_noise_loading(:, :)
      real(dp), allocatable :: observation_covariance(:, :)
      real(dp), allocatable :: cross_covariance(:, :)
      real(dp), allocatable :: transition_jacobian(:, :)
      real(dp), allocatable :: state_noise_loading(:, :)
      real(dp), allocatable :: transition_mean(:), factor(:, :)
      real(dp), allocatable :: compact_covariance(:, :), inverse(:, :)
      real(dp), allocatable :: compact_innovation(:), gain(:, :)
      integer, allocatable :: observed(:)
      real(dp) :: lambda, sigma_scale, log_determinant
      real(dp), allocatable :: observation_difference(:)
      integer :: series, state, times, sigma_count, observed_count
      integer :: point, time, component, info

      series = size(y, 1)
      times = size(y, 2)
      state = size(initial_mean)
      lambda = alpha*alpha*(real(state, dp) + kappa) - real(state, dp)
      if (series < 1 .or. state < 1 .or. times < 1 .or. &
         observation_noise_dimension < 1 .or. state_noise_dimension < 1 .or. &
         alpha <= 0.0_dp .or. beta < 0.0_dp .or. kappa < 0.0_dp .or. &
         real(state, dp) + lambda <= 0.0_dp .or. &
         any(shape(initial_covariance) /= [state, state]) .or. &
         .not. all(ieee_is_finite(initial_mean)) .or. &
         .not. all(ieee_is_finite(initial_covariance)) .or. &
         .not. all(ieee_is_finite(parameters))) then
         out%info = 1
         return
      end if
      sigma_count = 2*state + 1
      sigma_scale = sqrt(real(state, dp) + lambda)
      allocate(mean_weight(sigma_count), covariance_weight(sigma_count))
      mean_weight(1) = lambda/(real(state, dp) + lambda)
      mean_weight(2:) = 1.0_dp/(2.0_dp*(real(state, dp) + lambda))
      covariance_weight = mean_weight
      covariance_weight(1) = covariance_weight(1) + &
         1.0_dp - alpha*alpha + beta
      allocate(out%predicted_mean(state, times + 1))
      allocate(out%filtered_mean(state, times))
      allocate(out%predicted_covariance(state, state, times + 1))
      allocate(out%filtered_covariance(state, state, times))
      allocate(out%innovation(series, times))
      allocate(out%innovation_covariance(series, series, times))
      allocate(out%iterations(times))
      allocate(sigma(state, sigma_count), propagated(state, sigma_count))
      allocate(sigma_y(series, sigma_count), predicted_observation(series))
      allocate(observation_jacobian(series, state))
      allocate(observation_noise_loading(series, observation_noise_dimension))
      allocate(observation_covariance(series, series))
      allocate(cross_covariance(state, series))
      allocate(transition_jacobian(state, state))
      allocate(state_noise_loading(state, state_noise_dimension))
      allocate(transition_mean(state), observation_difference(series))
      out%predicted_mean(:, 1) = initial_mean
      out%predicted_covariance(:, :, 1) = initial_covariance
      out%innovation = 0.0_dp
      out%innovation_covariance = 0.0_dp
      out%iterations = 0
      do time = 1, times
         call positive_semidefinite_factor( &
            out%predicted_covariance(:, :, time), factor, info)
         if (info /= 0) then
            out%info = 2
            return
         end if
         call form_sigma_points(out%predicted_mean(:, time), factor, &
            sigma_scale, sigma)
         observed_count = count(ieee_is_finite(y(:, time)))
         if (observed_count > 0) then
            do point = 1, sigma_count
               call observation(time, sigma(:, point), parameters, &
                  sigma_y(:, point), observation_jacobian, &
                  observation_noise_loading)
               if (.not. all(ieee_is_finite(sigma_y(:, point)))) then
                  out%info = 3
                  return
               end if
            end do
            call observation(time, out%predicted_mean(:, time), parameters, &
               predicted_observation, observation_jacobian, &
               observation_noise_loading)
            if (.not. all(ieee_is_finite(observation_noise_loading))) then
               out%info = 3
               return
            end if
            predicted_observation = matmul(sigma_y, mean_weight)
            observation_covariance = matmul(observation_noise_loading, &
               transpose(observation_noise_loading))
            cross_covariance = 0.0_dp
            do point = 1, sigma_count
               observation_difference = sigma_y(:, point) - &
                  predicted_observation
               observation_covariance = observation_covariance + &
                  covariance_weight(point)*outer_product( &
                  observation_difference, observation_difference)
               cross_covariance = cross_covariance + &
                  covariance_weight(point)*outer_product( &
                  sigma(:, point) - out%predicted_mean(:, time), &
                  observation_difference)
            end do
            observation_covariance = 0.5_dp*(observation_covariance + &
               transpose(observation_covariance))
            allocate(observed(observed_count))
            observed = pack([(component, component = 1, series)], &
               ieee_is_finite(y(:, time)))
            allocate(compact_covariance(observed_count, observed_count))
            allocate(inverse(observed_count, observed_count))
            allocate(compact_innovation(observed_count))
            allocate(gain(state, observed_count))
            compact_covariance = observation_covariance(observed, observed)
            call inverse_logdet(compact_covariance, inverse, log_determinant, &
               info, 1.0e-12_dp)
            if (info /= 0) then
               out%info = 4
               return
            end if
            compact_innovation = y(observed, time) - &
               predicted_observation(observed)
            gain = matmul(cross_covariance(:, observed), inverse)
            out%innovation(observed, time) = compact_innovation
            do component = 1, observed_count
               out%innovation_covariance(observed, observed(component), &
                  time) = compact_covariance(:, component)
            end do
            out%filtered_mean(:, time) = out%predicted_mean(:, time) + &
               matmul(gain, compact_innovation)
            out%filtered_covariance(:, :, time) = &
               out%predicted_covariance(:, :, time) - &
               matmul(matmul(gain, compact_covariance), transpose(gain))
            out%filtered_covariance(:, :, time) = &
               0.5_dp*(out%filtered_covariance(:, :, time) + &
               transpose(out%filtered_covariance(:, :, time)))
            out%log_likelihood = out%log_likelihood - 0.5_dp* &
               (real(observed_count, dp)*log(2.0_dp*acos(-1.0_dp)) + &
               log_determinant + dot_product(compact_innovation, &
               matmul(inverse, compact_innovation)))
            deallocate(observed, compact_covariance, inverse)
            deallocate(compact_innovation, gain)
         else
            out%filtered_mean(:, time) = out%predicted_mean(:, time)
            out%filtered_covariance(:, :, time) = &
               out%predicted_covariance(:, :, time)
         end if
         call positive_semidefinite_factor( &
            out%filtered_covariance(:, :, time), factor, info)
         if (info /= 0) then
            out%info = 5
            return
         end if
         call form_sigma_points(out%filtered_mean(:, time), factor, &
            sigma_scale, sigma)
         do point = 1, sigma_count
            call transition(time, sigma(:, point), parameters, &
               propagated(:, point), transition_jacobian, state_noise_loading)
            if (.not. all(ieee_is_finite(propagated(:, point)))) then
               out%info = 6
               return
            end if
         end do
         out%predicted_mean(:, time + 1) = matmul(propagated, mean_weight)
         call transition(time, out%filtered_mean(:, time), parameters, &
            transition_mean, transition_jacobian, state_noise_loading)
         if (.not. all(ieee_is_finite(state_noise_loading))) then
            out%info = 6
            return
         end if
         out%predicted_covariance(:, :, time + 1) = &
            matmul(state_noise_loading, transpose(state_noise_loading))
         do point = 1, sigma_count
            out%predicted_covariance(:, :, time + 1) = &
               out%predicted_covariance(:, :, time + 1) + &
               covariance_weight(point)*outer_product( &
               propagated(:, point) - out%predicted_mean(:, time + 1), &
               propagated(:, point) - out%predicted_mean(:, time + 1))
         end do
         out%predicted_covariance(:, :, time + 1) = &
            0.5_dp*(out%predicted_covariance(:, :, time + 1) + &
            transpose(out%predicted_covariance(:, :, time + 1)))
      end do
   end function bssm_multivariate_ukf

   pure function bssm_gaussian_approximation(y, observation_loading, &
      transition, state_noise_loading, initial_mean, initial_covariance, &
      distribution, phi, max_iterations, convergence_tolerance, offset, &
      auxiliary, state_offset, initial_mode) result(out)
      !! Find a Laplace Gaussian approximation with the same conditional mode.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: observation_loading(:, :) !! Observation loading matrix.
      real(dp), intent(in) :: transition(:, :, :) !! State transition matrix.
      real(dp), intent(in) :: state_noise_loading(:, :, :) !! State noise loading.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      integer, intent(in) :: distribution !! Probability-distribution specification.
      integer, intent(in) :: max_iterations !! Maximum number of algorithm iterations.
      real(dp), intent(in) :: phi !! Autoregressive or model coefficient.
      real(dp), intent(in) :: convergence_tolerance !! Convergence tolerance.
      real(dp), intent(in), optional :: offset(:) !! Known additive offset.
      real(dp), intent(in), optional :: auxiliary(:) !! Auxiliary.
      real(dp), intent(in), optional :: state_offset(:, :) !! State offset.
      real(dp), intent(in), optional :: initial_mode(:) !! Initial mode.
      type(bssm_gaussian_approximation_t) :: out
      type(ssm_model_t) :: model
      type(kfs_filter_t) :: filtered
      type(kfs_smoother_t) :: smoothed
      real(dp), allocatable :: offset_work(:), auxiliary_work(:)
      real(dp), allocatable :: state_offset_work(:, :), center(:, :)
      real(dp), allocatable :: mode(:), next_mode(:), next_state(:, :)
      real(dp), allocatable :: factor(:, :)
      integer, allocatable :: status(:)
      integer :: state, noise, times, time, component, iteration, info

      state = size(initial_mean)
      noise = size(state_noise_loading, 2)
      times = size(y)
      if (state < 1 .or. noise < 1 .or. times < 1 .or. &
         max_iterations < 1 .or. convergence_tolerance < 0.0_dp .or. &
         distribution < bssm_svm .or. distribution > bssm_gamma .or. &
         any(shape(observation_loading) /= [state, times]) .or. &
         any(shape(transition) /= [state, state, times]) .or. &
         any(shape(state_noise_loading) /= [state, noise, times]) .or. &
         any(shape(initial_covariance) /= [state, state]) .or. &
         .not. all(ieee_is_finite(initial_mean)) .or. &
         .not. all(ieee_is_finite(initial_covariance))) then
         out%info = 1
         return
      end if
      allocate(offset_work(times), auxiliary_work(times))
      allocate(state_offset_work(state, times), center(state, times))
      offset_work = 0.0_dp
      auxiliary_work = 1.0_dp
      state_offset_work = 0.0_dp
      if (present(offset)) then
         if (size(offset) /= times .or. .not. all(ieee_is_finite(offset))) then
            out%info = 1
            return
         end if
         offset_work = offset
      end if
      if (present(auxiliary)) then
         if (size(auxiliary) /= times .or. &
            .not. all(ieee_is_finite(auxiliary))) then
            out%info = 1
            return
         end if
         auxiliary_work = auxiliary
      end if
      if (present(state_offset)) then
         if (any(shape(state_offset) /= [state, times]) .or. &
            .not. all(ieee_is_finite(state_offset))) then
            out%info = 1
            return
         end if
         state_offset_work = state_offset
      end if
      center(:, 1) = 0.0_dp
      do time = 1, times - 1
         center(:, time + 1) = state_offset_work(:, time) + &
            matmul(transition(:, :, time), center(:, time))
      end do
      allocate(mode(times), next_mode(times), next_state(state, times))
      mode(1) = offset_work(1) + dot_product(observation_loading(:, 1), &
         initial_mean)
      next_state(:, 1) = initial_mean
      do time = 2, times
         next_state(:, time) = state_offset_work(:, time - 1) + &
            matmul(transition(:, :, time - 1), next_state(:, time - 1))
         mode(time) = offset_work(time) + &
            dot_product(observation_loading(:, time), next_state(:, time))
      end do
      if (present(initial_mode)) then
         if (size(initial_mode) /= times .or. &
            .not. all(ieee_is_finite(initial_mode))) then
            out%info = 1
            return
         end if
         mode = initial_mode
      end if
      allocate(out%pseudo_observation(times), out%observation_variance(times))
      allocate(out%mode_state(state, times), out%mode_signal(times))
      allocate(out%scaling(times), status(times))
      allocate(model%y(times, 1), model%z(1, state, times))
      allocate(model%h(1, 1, times), model%transition(state, state, times))
      allocate(model%r(state, noise, times), model%q(noise, noise, 1))
      allocate(model%a1(state), model%p1(state, state))
      allocate(model%p1inf(state, state), model%missing(times, 1))
      model%z(1, :, :) = observation_loading
      model%transition = transition
      model%r = state_noise_loading
      model%q = 0.0_dp
      do component = 1, noise
         model%q(component, component, 1) = 1.0_dp
      end do
      model%a1 = initial_mean
      model%p1 = initial_covariance
      model%p1inf = 0.0_dp
      model%missing(:, 1) = .not. ieee_is_finite(y)
      do iteration = 1, max_iterations
         out%iterations = iteration
         call bssm_laplace_pseudo_observation(y, mode, distribution, phi, &
            auxiliary_work, out%pseudo_observation, &
            out%observation_variance, status)
         if (any(status /= 0)) then
            out%info = 2
            return
         end if
         do time = 1, times
            model%y(time, 1) = out%pseudo_observation(time) - &
               offset_work(time) - dot_product(observation_loading(:, time), &
               center(:, time))
            model%h(1, 1, time) = out%observation_variance(time)
         end do
         filtered = kfs_filter(model)
         if (filtered%info /= 0) then
            out%info = 100 + filtered%info
            return
         end if
         smoothed = kfs_smooth(model, filtered)
         if (smoothed%info /= 0) then
            out%info = 200 + smoothed%info
            return
         end if
         next_state = smoothed%state + center
         do time = 1, times
            next_mode(time) = offset_work(time) + &
               dot_product(observation_loading(:, time), next_state(:, time))
         end do
         out%difference = sum((next_mode - mode)**2)/real(times, dp)
         mode = next_mode
         out%mode_state = next_state
         if (out%difference <= convergence_tolerance) then
            out%converged = .true.
            exit
         end if
      end do
      out%mode_signal = mode
      out%gaussian_log_likelihood = filtered%log_likelihood
      out%scaling = 0.0_dp
      do time = 1, times
         if (.not. ieee_is_finite(y(time))) cycle
         out%scaling(time) = bssm_observation_log_density(y(time), mode(time), &
            distribution, phi, auxiliary_work(time)) - &
            normal_log_density(out%pseudo_observation(time), mode(time), &
            sqrt(out%observation_variance(time)))
      end do
      out%corrected_log_likelihood = out%gaussian_log_likelihood + &
         sum(out%scaling)
      allocate(out%proposal_factor(state, state, times))
      allocate(out%conditional_matrix(state, state, times))
      out%conditional_matrix = smoothed%conditional_matrix
      do time = 1, times
         call positive_semidefinite_factor( &
            smoothed%conditional_covariance(:, :, time), factor, info)
         if (info /= 0) then
            out%info = 300 + time
            return
         end if
         out%proposal_factor(:, :, time) = factor
      end do
   end function bssm_gaussian_approximation

   pure function bssm_multivariate_gaussian_approximation(y, &
      observation_loading, transition, state_noise_loading, initial_mean, &
      initial_covariance, distribution, phi, max_iterations, &
      convergence_tolerance, offset, auxiliary, state_offset, initial_mode) &
      result(out)
      !! Find a mixed-family multivariate Laplace Gaussian approximation.
      real(dp), intent(in) :: y(:, :) !! Response or time-series observations.
      real(dp), intent(in) :: observation_loading(:, :, :) !! Observation loading matrix.
      real(dp), intent(in) :: transition(:, :, :) !! State transition matrix.
      real(dp), intent(in) :: state_noise_loading(:, :, :) !! State noise loading.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      integer, intent(in) :: distribution(:) !! Probability-distribution specification.
      integer, intent(in) :: max_iterations !! Maximum number of algorithm iterations.
      real(dp), intent(in) :: phi(:) !! Autoregressive or model coefficient.
      real(dp), intent(in) :: convergence_tolerance !! Convergence tolerance.
      real(dp), intent(in), optional :: offset(:, :) !! Known additive offset.
      real(dp), intent(in), optional :: auxiliary(:, :) !! Auxiliary.
      real(dp), intent(in), optional :: state_offset(:, :) !! State offset.
      real(dp), intent(in), optional :: initial_mode(:, :) !! Initial mode.
      type(bssm_multivariate_approximation_t) :: out
      type(ssm_model_t) :: model
      type(kfs_filter_t) :: filtered
      type(kfs_smoother_t) :: smoothed
      real(dp), allocatable :: offset_work(:, :), auxiliary_work(:, :)
      real(dp), allocatable :: state_offset_work(:, :), center(:, :)
      real(dp), allocatable :: mode(:, :), next_mode(:, :)
      real(dp), allocatable :: next_state(:, :)
      real(dp), allocatable :: factor(:, :)
      integer, allocatable :: status(:, :)
      integer :: series, state, noise, times
      integer :: time, component, iteration, info

      series = size(y, 1)
      times = size(y, 2)
      state = size(initial_mean)
      noise = size(state_noise_loading, 2)
      if (series < 1 .or. state < 1 .or. noise < 1 .or. times < 1 .or. &
         max_iterations < 1 .or. convergence_tolerance < 0.0_dp .or. &
         size(distribution) /= series .or. size(phi) /= series .or. &
         any(distribution < bssm_poisson) .or. &
         any(distribution > bssm_gaussian) .or. &
         any(shape(observation_loading) /= [series, state, times]) .or. &
         any(shape(transition) /= [state, state, times]) .or. &
         any(shape(state_noise_loading) /= [state, noise, times]) .or. &
         any(shape(initial_covariance) /= [state, state]) .or. &
         .not. all(ieee_is_finite(phi)) .or. &
         .not. all(ieee_is_finite(initial_mean)) .or. &
         .not. all(ieee_is_finite(initial_covariance))) then
         out%info = 1
         return
      end if
      allocate(offset_work(series, times), auxiliary_work(series, times))
      allocate(state_offset_work(state, times), center(state, times))
      offset_work = 0.0_dp
      auxiliary_work = 1.0_dp
      state_offset_work = 0.0_dp
      if (present(offset)) then
         if (any(shape(offset) /= [series, times]) .or. &
            .not. all(ieee_is_finite(offset))) then
            out%info = 1
            return
         end if
         offset_work = offset
      end if
      if (present(auxiliary)) then
         if (any(shape(auxiliary) /= [series, times]) .or. &
            .not. all(ieee_is_finite(auxiliary))) then
            out%info = 1
            return
         end if
         auxiliary_work = auxiliary
      end if
      if (present(state_offset)) then
         if (any(shape(state_offset) /= [state, times]) .or. &
            .not. all(ieee_is_finite(state_offset))) then
            out%info = 1
            return
         end if
         state_offset_work = state_offset
      end if
      center(:, 1) = 0.0_dp
      do time = 1, times - 1
         center(:, time + 1) = state_offset_work(:, time) + &
            matmul(transition(:, :, time), center(:, time))
      end do
      allocate(mode(series, times), next_mode(series, times))
      allocate(next_state(state, times))
      next_state(:, 1) = initial_mean
      do time = 1, times
         if (time > 1) then
            next_state(:, time) = state_offset_work(:, time - 1) + &
               matmul(transition(:, :, time - 1), next_state(:, time - 1))
         end if
         mode(:, time) = offset_work(:, time) + &
            matmul(observation_loading(:, :, time), next_state(:, time))
      end do
      if (present(initial_mode)) then
         if (any(shape(initial_mode) /= [series, times]) .or. &
            .not. all(ieee_is_finite(initial_mode))) then
            out%info = 1
            return
         end if
         mode = initial_mode
      end if
      allocate(out%pseudo_observation(series, times))
      allocate(out%observation_variance(series, times))
      allocate(out%mode_state(state, times), out%mode_signal(series, times))
      allocate(out%scaling(series, times), status(series, times))
      allocate(model%y(times, series), model%z(series, state, times))
      allocate(model%h(series, series, times))
      allocate(model%transition(state, state, times))
      allocate(model%r(state, noise, times), model%q(noise, noise, 1))
      allocate(model%a1(state), model%p1(state, state))
      allocate(model%p1inf(state, state), model%missing(times, series))
      model%z = observation_loading
      model%transition = transition
      model%r = state_noise_loading
      model%q = 0.0_dp
      do component = 1, noise
         model%q(component, component, 1) = 1.0_dp
      end do
      model%a1 = initial_mean
      model%p1 = initial_covariance
      model%p1inf = 0.0_dp
      model%missing = transpose(.not. ieee_is_finite(y))
      do iteration = 1, max_iterations
         out%iterations = iteration
         do component = 1, series
            call bssm_laplace_pseudo_observation(y(component, :), &
               mode(component, :), distribution(component), phi(component), &
               auxiliary_work(component, :), &
               out%pseudo_observation(component, :), &
               out%observation_variance(component, :), status(component, :))
         end do
         if (any(status /= 0)) then
            out%info = 2
            return
         end if
         model%h = 0.0_dp
         do time = 1, times
            do component = 1, series
               model%y(time, component) = &
                  out%pseudo_observation(component, time) - &
                  offset_work(component, time) - &
                  dot_product(observation_loading(component, :, time), &
                  center(:, time))
               model%h(component, component, time) = &
                  out%observation_variance(component, time)
            end do
         end do
         filtered = kfs_filter(model)
         if (filtered%info /= 0) then
            out%info = 100 + filtered%info
            return
         end if
         smoothed = kfs_smooth(model, filtered)
         if (smoothed%info /= 0) then
            out%info = 200 + smoothed%info
            return
         end if
         next_state = smoothed%state + center
         do time = 1, times
            next_mode(:, time) = offset_work(:, time) + &
               matmul(observation_loading(:, :, time), next_state(:, time))
         end do
         out%difference = sum((next_mode - mode)**2)/real(series*times, dp)
         mode = next_mode
         out%mode_state = next_state
         if (out%difference <= convergence_tolerance) then
            out%converged = .true.
            exit
         end if
      end do
      out%mode_signal = mode
      out%gaussian_log_likelihood = filtered%log_likelihood
      out%scaling = 0.0_dp
      do time = 1, times
         do component = 1, series
            if (.not. ieee_is_finite(y(component, time))) cycle
            out%scaling(component, time) = bssm_observation_log_density( &
               y(component, time), mode(component, time), &
               distribution(component), phi(component), &
               auxiliary_work(component, time)) - normal_log_density( &
               out%pseudo_observation(component, time), mode(component, time), &
               sqrt(out%observation_variance(component, time)))
         end do
      end do
      out%corrected_log_likelihood = out%gaussian_log_likelihood + &
         sum(out%scaling)
      allocate(out%proposal_factor(state, state, times))
      allocate(out%conditional_matrix(state, state, times))
      out%conditional_matrix = smoothed%conditional_matrix
      do time = 1, times
         call positive_semidefinite_factor( &
            smoothed%conditional_covariance(:, :, time), factor, info)
         if (info /= 0) then
            out%info = 300 + time
            return
         end if
         out%proposal_factor(:, :, time) = factor
      end do
   end function bssm_multivariate_gaussian_approximation

   pure function bssm_multivariate_bootstrap_filter_draws(y, &
      observation_loading, transition, state_noise_loading, initial_mean, &
      initial_covariance, distribution, phi, initial_normals, &
      innovation_normals, resampling_uniforms, offset, auxiliary, &
      state_offset) result(out)
      !! Run a draw-driven bootstrap filter for mixed observation families.
      real(dp), intent(in) :: y(:, :) !! Response or time-series observations.
      real(dp), intent(in) :: observation_loading(:, :, :) !! Observation loading matrix.
      real(dp), intent(in) :: transition(:, :, :) !! State transition matrix.
      real(dp), intent(in) :: state_noise_loading(:, :, :) !! State noise loading.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      integer, intent(in) :: distribution(:) !! Probability-distribution specification.
      real(dp), intent(in) :: phi(:) !! Autoregressive or model coefficient.
      real(dp), intent(in) :: initial_normals(:, :) !! Initial normals.
      real(dp), intent(in) :: innovation_normals(:, :, :) !! Innovation normals.
      real(dp), intent(in) :: resampling_uniforms(:, :) !! Resampling uniforms.
      real(dp), intent(in), optional :: offset(:, :) !! Known additive offset.
      real(dp), intent(in), optional :: auxiliary(:, :) !! Auxiliary.
      real(dp), intent(in), optional :: state_offset(:, :) !! State offset.
      type(bssm_particle_filter_t) :: out
      real(dp), allocatable :: offset_work(:, :), auxiliary_work(:, :)
      real(dp), allocatable :: state_offset_work(:, :), factor(:, :)
      real(dp), allocatable :: parent(:, :), log_weight(:), probability(:)
      real(dp) :: increment, signal
      integer :: series, state, noise, particles, times
      integer :: component, particle, time, info

      series = size(y, 1)
      times = size(y, 2)
      state = size(initial_mean)
      noise = size(state_noise_loading, 2)
      particles = size(initial_normals, 2)
      if (series < 1 .or. times < 1 .or. state < 1 .or. noise < 1 .or. &
         particles < 1 .or. size(distribution) /= series .or. &
         size(phi) /= series .or. any(distribution < bssm_poisson) .or. &
         any(distribution > bssm_gaussian) .or. &
         any(shape(observation_loading) /= [series, state, times]) .or. &
         any(shape(transition) /= [state, state, times]) .or. &
         any(shape(state_noise_loading) /= [state, noise, times]) .or. &
         any(shape(initial_covariance) /= [state, state]) .or. &
         size(initial_normals, 1) /= state .or. &
         any(shape(innovation_normals) /= [noise, particles, times]) .or. &
         any(shape(resampling_uniforms) /= [particles, times]) .or. &
         .not. all(ieee_is_finite(phi)) .or. &
         .not. all(ieee_is_finite(initial_mean)) .or. &
         .not. all(ieee_is_finite(initial_covariance)) .or. &
         .not. all(ieee_is_finite(initial_normals)) .or. &
         .not. all(ieee_is_finite(innovation_normals)) .or. &
         .not. all(ieee_is_finite(resampling_uniforms))) then
         out%info = 1
         return
      end if
      allocate(offset_work(series, times), auxiliary_work(series, times))
      allocate(state_offset_work(state, times))
      offset_work = 0.0_dp
      auxiliary_work = 1.0_dp
      state_offset_work = 0.0_dp
      if (present(offset)) then
         if (any(shape(offset) /= [series, times]) .or. &
            .not. all(ieee_is_finite(offset))) then
            out%info = 1
            return
         end if
         offset_work = offset
      end if
      if (present(auxiliary)) then
         if (any(shape(auxiliary) /= [series, times]) .or. &
            .not. all(ieee_is_finite(auxiliary))) then
            out%info = 1
            return
         end if
         auxiliary_work = auxiliary
      end if
      if (present(state_offset)) then
         if (any(shape(state_offset) /= [state, times]) .or. &
            .not. all(ieee_is_finite(state_offset))) then
            out%info = 1
            return
         end if
         state_offset_work = state_offset
      end if
      call positive_semidefinite_factor(initial_covariance, factor, info)
      if (info /= 0) then
         out%info = 2
         return
      end if
      allocate(out%particles(state, particles, times + 1))
      allocate(out%weights(particles, times), out%ancestors(particles, times))
      allocate(out%predicted_mean(state, times + 1))
      allocate(out%filtered_mean(state, times))
      allocate(out%predicted_covariance(state, state, times + 1))
      allocate(out%filtered_covariance(state, state, times))
      allocate(parent(state, particles), log_weight(particles))
      allocate(probability(particles))
      do particle = 1, particles
         out%particles(:, particle, 1) = initial_mean + &
            matmul(factor, initial_normals(:, particle))
      end do
      do time = 1, times
         call particle_summary(out%particles(:, :, time), &
            out%predicted_mean(:, time), &
            out%predicted_covariance(:, :, time))
         log_weight = 0.0_dp
         do particle = 1, particles
            do component = 1, series
               if (.not. ieee_is_finite(y(component, time))) cycle
               signal = offset_work(component, time) + &
                  dot_product(observation_loading(component, :, time), &
                  out%particles(:, particle, time))
               log_weight(particle) = log_weight(particle) + &
                  bssm_observation_log_density(y(component, time), signal, &
                  distribution(component), phi(component), &
                  auxiliary_work(component, time))
            end do
         end do
         call normalize_log_weights(log_weight, probability, increment, info)
         if (info /= 0) then
            out%info = 3
            out%log_likelihood = -huge(1.0_dp)
            return
         end if
         out%log_likelihood = out%log_likelihood + increment
         out%weights(:, time) = probability
         call weighted_particle_summary(out%particles(:, :, time), &
            probability, out%filtered_mean(:, time), &
            out%filtered_covariance(:, :, time))
         out%ancestors(:, time) = bssm_stratified_resample(probability, &
            resampling_uniforms(:, time))
         do particle = 1, particles
            parent(:, particle) = out%particles(:, &
               out%ancestors(particle, time), time)
            out%particles(:, particle, time + 1) = &
               state_offset_work(:, time) + &
               matmul(transition(:, :, time), parent(:, particle)) + &
               matmul(state_noise_loading(:, :, time), &
               innovation_normals(:, particle, time))
         end do
      end do
      call particle_summary(out%particles(:, :, times + 1), &
         out%predicted_mean(:, times + 1), &
         out%predicted_covariance(:, :, times + 1))
   end function bssm_multivariate_bootstrap_filter_draws

   function bssm_multivariate_bootstrap_filter(y, observation_loading, &
      transition, state_noise_loading, initial_mean, initial_covariance, &
      distribution, phi, particles, offset, auxiliary, state_offset) &
      result(out)
      !! Run a multivariate mixed-family bootstrap filter with shared randomness.
      real(dp), intent(in) :: y(:, :) !! Response or time-series observations.
      real(dp), intent(in) :: observation_loading(:, :, :) !! Observation loading matrix.
      real(dp), intent(in) :: transition(:, :, :) !! State transition matrix.
      real(dp), intent(in) :: state_noise_loading(:, :, :) !! State noise loading.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      integer, intent(in) :: distribution(:) !! Probability-distribution specification.
      integer, intent(in) :: particles !! Number of particles.
      real(dp), intent(in) :: phi(:) !! Autoregressive or model coefficient.
      real(dp), intent(in), optional :: offset(:, :) !! Known additive offset.
      real(dp), intent(in), optional :: auxiliary(:, :) !! Auxiliary.
      real(dp), intent(in), optional :: state_offset(:, :) !! State offset.
      type(bssm_particle_filter_t) :: out
      real(dp), allocatable :: initial_normals(:, :)
      real(dp), allocatable :: innovation_normals(:, :, :), uniforms(:, :)
      integer :: time

      if (particles < 1 .or. size(y, 1) < 1 .or. size(y, 2) < 1) then
         out%info = 1
         return
      end if
      allocate(initial_normals(size(initial_mean), particles))
      allocate(innovation_normals(size(state_noise_loading, 2), particles, &
         size(y, 2)), uniforms(particles, size(y, 2)))
      call random_standard_normal_matrix(initial_normals)
      do time = 1, size(y, 2)
         call random_standard_normal_matrix(innovation_normals(:, :, time))
      end do
      call random_number(uniforms)
      out = bssm_multivariate_bootstrap_filter_draws(y, observation_loading, &
         transition, state_noise_loading, initial_mean, initial_covariance, &
         distribution, phi, initial_normals, innovation_normals, uniforms, &
         offset, auxiliary, state_offset)
   end function bssm_multivariate_bootstrap_filter

   pure function bssm_multivariate_psi_filter_draws(y, observation_loading, &
      transition, state_noise_loading, initial_mean, initial_covariance, &
      distribution, phi, approximation_iterations, convergence_tolerance, &
      proposal_normals, terminal_normals, resampling_uniforms, offset, &
      auxiliary, state_offset, initial_mode) result(out)
      !! Run a draw-driven psi filter for mixed multivariate observations.
      real(dp), intent(in) :: y(:, :) !! Response or time-series observations.
      real(dp), intent(in) :: observation_loading(:, :, :) !! Observation loading matrix.
      real(dp), intent(in) :: transition(:, :, :) !! State transition matrix.
      real(dp), intent(in) :: state_noise_loading(:, :, :) !! State noise loading.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      integer, intent(in) :: distribution(:) !! Probability-distribution specification.
      integer, intent(in) :: approximation_iterations !! Number of approximation iterations.
      real(dp), intent(in) :: phi(:) !! Autoregressive or model coefficient.
      real(dp), intent(in) :: convergence_tolerance !! Convergence tolerance.
      real(dp), intent(in) :: proposal_normals(:, :, :) !! Standard-normal proposal draws.
      real(dp), intent(in) :: terminal_normals(:, :) !! Terminal normals.
      real(dp), intent(in) :: resampling_uniforms(:, :) !! Resampling uniforms.
      real(dp), intent(in), optional :: offset(:, :) !! Known additive offset.
      real(dp), intent(in), optional :: auxiliary(:, :) !! Auxiliary.
      real(dp), intent(in), optional :: state_offset(:, :) !! State offset.
      real(dp), intent(in), optional :: initial_mode(:, :) !! Initial mode.
      type(bssm_particle_filter_t) :: out
      type(bssm_multivariate_approximation_t) :: approximation
      real(dp), allocatable :: offset_work(:, :), auxiliary_work(:, :)
      real(dp), allocatable :: state_offset_work(:, :), parent(:, :)
      real(dp), allocatable :: log_weight(:), probability(:)
      real(dp) :: increment, signal, gaussian_log_density
      integer :: series, state, noise, particles, times
      integer :: component, particle, time, info

      series = size(y, 1)
      times = size(y, 2)
      state = size(initial_mean)
      noise = size(state_noise_loading, 2)
      particles = size(proposal_normals, 2)
      if (series < 1 .or. times < 1 .or. state < 1 .or. noise < 1 .or. &
         particles < 1 .or. size(distribution) /= series .or. &
         size(phi) /= series .or. any(distribution < bssm_poisson) .or. &
         any(distribution > bssm_gaussian) .or. &
         any(shape(observation_loading) /= [series, state, times]) .or. &
         any(shape(transition) /= [state, state, times]) .or. &
         any(shape(state_noise_loading) /= [state, noise, times]) .or. &
         any(shape(initial_covariance) /= [state, state]) .or. &
         any(shape(proposal_normals) /= [state, particles, times]) .or. &
         any(shape(terminal_normals) /= [noise, particles]) .or. &
         any(shape(resampling_uniforms) /= [particles, times]) .or. &
         .not. all(ieee_is_finite(proposal_normals)) .or. &
         .not. all(ieee_is_finite(terminal_normals)) .or. &
         .not. all(ieee_is_finite(resampling_uniforms))) then
         out%info = 1
         return
      end if
      allocate(offset_work(series, times), auxiliary_work(series, times))
      allocate(state_offset_work(state, times))
      offset_work = 0.0_dp
      auxiliary_work = 1.0_dp
      state_offset_work = 0.0_dp
      if (present(offset)) then
         if (any(shape(offset) /= [series, times]) .or. &
            .not. all(ieee_is_finite(offset))) then
            out%info = 1
            return
         end if
         offset_work = offset
      end if
      if (present(auxiliary)) then
         if (any(shape(auxiliary) /= [series, times]) .or. &
            .not. all(ieee_is_finite(auxiliary))) then
            out%info = 1
            return
         end if
         auxiliary_work = auxiliary
      end if
      if (present(state_offset)) then
         if (any(shape(state_offset) /= [state, times]) .or. &
            .not. all(ieee_is_finite(state_offset))) then
            out%info = 1
            return
         end if
         state_offset_work = state_offset
      end if
      approximation = bssm_multivariate_gaussian_approximation(y, &
         observation_loading, transition, state_noise_loading, initial_mean, &
         initial_covariance, distribution, phi, approximation_iterations, &
         convergence_tolerance, offset_work, auxiliary_work, &
         state_offset_work, initial_mode)
      if (approximation%info /= 0) then
         out%info = 1000 + approximation%info
         return
      end if
      allocate(out%particles(state, particles, times + 1))
      allocate(out%weights(particles, times), out%ancestors(particles, times))
      allocate(out%predicted_mean(state, times + 1))
      allocate(out%filtered_mean(state, times))
      allocate(out%predicted_covariance(state, state, times + 1))
      allocate(out%filtered_covariance(state, state, times))
      allocate(parent(state, particles), log_weight(particles))
      allocate(probability(particles))
      do particle = 1, particles
         out%particles(:, particle, 1) = approximation%mode_state(:, 1) + &
            matmul(approximation%proposal_factor(:, :, 1), &
            proposal_normals(:, particle, 1))
      end do
      out%log_likelihood = approximation%corrected_log_likelihood
      do time = 1, times
         call particle_summary(out%particles(:, :, time), &
            out%predicted_mean(:, time), &
            out%predicted_covariance(:, :, time))
         log_weight = 0.0_dp
         do particle = 1, particles
            do component = 1, series
               if (.not. ieee_is_finite(y(component, time))) cycle
               signal = offset_work(component, time) + &
                  dot_product(observation_loading(component, :, time), &
                  out%particles(:, particle, time))
               gaussian_log_density = normal_log_density( &
                  approximation%pseudo_observation(component, time), signal, &
                  sqrt(approximation%observation_variance(component, time)))
               log_weight(particle) = log_weight(particle) + &
                  bssm_observation_log_density(y(component, time), signal, &
                  distribution(component), phi(component), &
                  auxiliary_work(component, time)) - gaussian_log_density - &
                  approximation%scaling(component, time)
            end do
         end do
         call normalize_log_weights(log_weight, probability, increment, info)
         if (info /= 0) then
            out%info = 2
            out%log_likelihood = -huge(1.0_dp)
            return
         end if
         out%log_likelihood = out%log_likelihood + increment
         out%weights(:, time) = probability
         call weighted_particle_summary(out%particles(:, :, time), &
            probability, out%filtered_mean(:, time), &
            out%filtered_covariance(:, :, time))
         out%ancestors(:, time) = bssm_stratified_resample(probability, &
            resampling_uniforms(:, time))
         do particle = 1, particles
            parent(:, particle) = out%particles(:, &
               out%ancestors(particle, time), time)
            if (time < times) then
               out%particles(:, particle, time + 1) = &
                  approximation%mode_state(:, time + 1) + &
                  matmul(approximation%conditional_matrix(:, :, time + 1), &
                  parent(:, particle) - approximation%mode_state(:, time)) + &
                  matmul(approximation%proposal_factor(:, :, time + 1), &
                  proposal_normals(:, particle, time + 1))
            else
               out%particles(:, particle, time + 1) = &
                  state_offset_work(:, time) + &
                  matmul(transition(:, :, time), parent(:, particle)) + &
                  matmul(state_noise_loading(:, :, time), &
                  terminal_normals(:, particle))
            end if
         end do
      end do
      call particle_summary(out%particles(:, :, times + 1), &
         out%predicted_mean(:, times + 1), &
         out%predicted_covariance(:, :, times + 1))
   end function bssm_multivariate_psi_filter_draws

   function bssm_multivariate_psi_filter(y, observation_loading, transition, &
      state_noise_loading, initial_mean, initial_covariance, distribution, &
      phi, particles, approximation_iterations, convergence_tolerance, &
      offset, auxiliary, state_offset, initial_mode) result(out)
      !! Run a mixed-family multivariate psi filter with shared randomness.
      real(dp), intent(in) :: y(:, :) !! Response or time-series observations.
      real(dp), intent(in) :: observation_loading(:, :, :) !! Observation loading matrix.
      real(dp), intent(in) :: transition(:, :, :) !! State transition matrix.
      real(dp), intent(in) :: state_noise_loading(:, :, :) !! State noise loading.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      integer, intent(in) :: distribution(:) !! Probability-distribution specification.
      integer, intent(in) :: particles !! Number of particles.
      integer, intent(in) :: approximation_iterations !! Number of approximation iterations.
      real(dp), intent(in) :: phi(:) !! Autoregressive or model coefficient.
      real(dp), intent(in) :: convergence_tolerance !! Convergence tolerance.
      real(dp), intent(in), optional :: offset(:, :) !! Known additive offset.
      real(dp), intent(in), optional :: auxiliary(:, :) !! Auxiliary.
      real(dp), intent(in), optional :: state_offset(:, :) !! State offset.
      real(dp), intent(in), optional :: initial_mode(:, :) !! Initial mode.
      type(bssm_particle_filter_t) :: out
      real(dp), allocatable :: proposal_normals(:, :, :)
      real(dp), allocatable :: terminal_normals(:, :), uniforms(:, :)
      integer :: time

      if (particles < 1 .or. size(y, 1) < 1 .or. size(y, 2) < 1) then
         out%info = 1
         return
      end if
      allocate(proposal_normals(size(initial_mean), particles, size(y, 2)))
      allocate(terminal_normals(size(state_noise_loading, 2), particles))
      allocate(uniforms(particles, size(y, 2)))
      do time = 1, size(y, 2)
         call random_standard_normal_matrix(proposal_normals(:, :, time))
      end do
      call random_standard_normal_matrix(terminal_normals)
      call random_number(uniforms)
      out = bssm_multivariate_psi_filter_draws(y, observation_loading, &
         transition, state_noise_loading, initial_mean, initial_covariance, &
         distribution, phi, approximation_iterations, convergence_tolerance, &
         proposal_normals, terminal_normals, uniforms, offset, auxiliary, &
         state_offset, initial_mode)
   end function bssm_multivariate_psi_filter

   pure function bssm_multivariate_spdk_importance_draws(y, &
      observation_loading, transition, state_noise_loading, initial_mean, &
      initial_covariance, distribution, phi, approximation_iterations, &
      convergence_tolerance, proposal_normals, terminal_normals, offset, &
      auxiliary, state_offset, initial_mode) result(out)
      !! Draw mixed-family multivariate trajectories for SPDK sampling.
      real(dp), intent(in) :: y(:, :) !! Response or time-series observations.
      real(dp), intent(in) :: observation_loading(:, :, :) !! Observation loading matrix.
      real(dp), intent(in) :: transition(:, :, :) !! State transition matrix.
      real(dp), intent(in) :: state_noise_loading(:, :, :) !! State noise loading.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      integer, intent(in) :: distribution(:) !! Probability-distribution specification.
      integer, intent(in) :: approximation_iterations !! Number of approximation iterations.
      real(dp), intent(in) :: phi(:) !! Autoregressive or model coefficient.
      real(dp), intent(in) :: convergence_tolerance !! Convergence tolerance.
      real(dp), intent(in) :: proposal_normals(:, :, :) !! Standard-normal proposal draws.
      real(dp), intent(in) :: terminal_normals(:, :) !! Terminal normals.
      real(dp), intent(in), optional :: offset(:, :) !! Known additive offset.
      real(dp), intent(in), optional :: auxiliary(:, :) !! Auxiliary.
      real(dp), intent(in), optional :: state_offset(:, :) !! State offset.
      real(dp), intent(in), optional :: initial_mode(:, :) !! Initial mode.
      type(bssm_importance_sample_t) :: out
      type(bssm_multivariate_approximation_t) :: approximation
      real(dp), allocatable :: offset_work(:, :), auxiliary_work(:, :)
      real(dp), allocatable :: state_offset_work(:, :)
      real(dp) :: increment, signal, gaussian_log_density
      integer :: series, state, noise, samples, times
      integer :: component, sample, time, info

      series = size(y, 1)
      times = size(y, 2)
      state = size(initial_mean)
      noise = size(state_noise_loading, 2)
      samples = size(proposal_normals, 2)
      if (series < 1 .or. times < 1 .or. state < 1 .or. noise < 1 .or. &
         samples < 1 .or. size(distribution) /= series .or. &
         size(phi) /= series .or. any(distribution < bssm_poisson) .or. &
         any(distribution > bssm_gaussian) .or. &
         any(shape(observation_loading) /= [series, state, times]) .or. &
         any(shape(transition) /= [state, state, times]) .or. &
         any(shape(state_noise_loading) /= [state, noise, times]) .or. &
         any(shape(initial_covariance) /= [state, state]) .or. &
         any(shape(proposal_normals) /= [state, samples, times]) .or. &
         any(shape(terminal_normals) /= [noise, samples]) .or. &
         .not. all(ieee_is_finite(proposal_normals)) .or. &
         .not. all(ieee_is_finite(terminal_normals))) then
         out%info = 1
         return
      end if
      allocate(offset_work(series, times), auxiliary_work(series, times))
      allocate(state_offset_work(state, times))
      offset_work = 0.0_dp
      auxiliary_work = 1.0_dp
      state_offset_work = 0.0_dp
      if (present(offset)) then
         if (any(shape(offset) /= [series, times]) .or. &
            .not. all(ieee_is_finite(offset))) then
            out%info = 1
            return
         end if
         offset_work = offset
      end if
      if (present(auxiliary)) then
         if (any(shape(auxiliary) /= [series, times]) .or. &
            .not. all(ieee_is_finite(auxiliary))) then
            out%info = 1
            return
         end if
         auxiliary_work = auxiliary
      end if
      if (present(state_offset)) then
         if (any(shape(state_offset) /= [state, times]) .or. &
            .not. all(ieee_is_finite(state_offset))) then
            out%info = 1
            return
         end if
         state_offset_work = state_offset
      end if
      approximation = bssm_multivariate_gaussian_approximation(y, &
         observation_loading, transition, state_noise_loading, initial_mean, &
         initial_covariance, distribution, phi, approximation_iterations, &
         convergence_tolerance, offset_work, auxiliary_work, &
         state_offset_work, initial_mode)
      if (approximation%info /= 0) then
         out%info = 1000 + approximation%info
         return
      end if
      allocate(out%trajectories(state, times + 1, samples))
      allocate(out%log_weight(samples), out%weight(samples))
      allocate(out%mean(state, times + 1))
      allocate(out%covariance(state, state, times + 1))
      do sample = 1, samples
         out%trajectories(:, 1, sample) = approximation%mode_state(:, 1) + &
            matmul(approximation%proposal_factor(:, :, 1), &
            proposal_normals(:, sample, 1))
         do time = 2, times
            out%trajectories(:, time, sample) = &
               approximation%mode_state(:, time) + &
               matmul(approximation%conditional_matrix(:, :, time), &
               out%trajectories(:, time - 1, sample) - &
               approximation%mode_state(:, time - 1)) + &
               matmul(approximation%proposal_factor(:, :, time), &
               proposal_normals(:, sample, time))
         end do
         out%trajectories(:, times + 1, sample) = &
            state_offset_work(:, times) + matmul(transition(:, :, times), &
            out%trajectories(:, times, sample)) + &
            matmul(state_noise_loading(:, :, times), &
            terminal_normals(:, sample))
      end do
      out%log_weight = 0.0_dp
      do sample = 1, samples
         do time = 1, times
            do component = 1, series
               if (.not. ieee_is_finite(y(component, time))) cycle
               signal = offset_work(component, time) + &
                  dot_product(observation_loading(component, :, time), &
                  out%trajectories(:, time, sample))
               gaussian_log_density = normal_log_density( &
                  approximation%pseudo_observation(component, time), signal, &
                  sqrt(approximation%observation_variance(component, time)))
               out%log_weight(sample) = out%log_weight(sample) + &
                  bssm_observation_log_density(y(component, time), signal, &
                  distribution(component), phi(component), &
                  auxiliary_work(component, time)) - gaussian_log_density - &
                  approximation%scaling(component, time)
            end do
         end do
      end do
      call normalize_log_weights(out%log_weight, out%weight, increment, info)
      if (info /= 0) then
         out%info = 2
         out%log_likelihood = -huge(1.0_dp)
         return
      end if
      out%log_likelihood = approximation%corrected_log_likelihood + increment
      out%effective_sample_size = 1.0_dp/sum(out%weight**2)
      do time = 1, times + 1
         call weighted_particle_summary(out%trajectories(:, time, :), &
            out%weight, out%mean(:, time), out%covariance(:, :, time))
      end do
   end function bssm_multivariate_spdk_importance_draws

   function bssm_multivariate_spdk_importance(y, observation_loading, &
      transition, state_noise_loading, initial_mean, initial_covariance, &
      distribution, phi, samples, approximation_iterations, &
      convergence_tolerance, use_antithetic, offset, auxiliary, state_offset, &
      initial_mode) result(out)
      !! Run multivariate SPDK importance sampling with shared randomness.
      real(dp), intent(in) :: y(:, :) !! Response or time-series observations.
      real(dp), intent(in) :: observation_loading(:, :, :) !! Observation loading matrix.
      real(dp), intent(in) :: transition(:, :, :) !! State transition matrix.
      real(dp), intent(in) :: state_noise_loading(:, :, :) !! State noise loading.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      integer, intent(in) :: distribution(:) !! Probability-distribution specification.
      integer, intent(in) :: samples !! Samples.
      integer, intent(in) :: approximation_iterations !! Number of approximation iterations.
      real(dp), intent(in) :: phi(:) !! Autoregressive or model coefficient.
      real(dp), intent(in) :: convergence_tolerance !! Convergence tolerance.
      logical, intent(in) :: use_antithetic !! Whether to use the antithetic.
      real(dp), intent(in), optional :: offset(:, :) !! Known additive offset.
      real(dp), intent(in), optional :: auxiliary(:, :) !! Auxiliary.
      real(dp), intent(in), optional :: state_offset(:, :) !! State offset.
      real(dp), intent(in), optional :: initial_mode(:, :) !! Initial mode.
      type(bssm_importance_sample_t) :: out
      real(dp), allocatable :: proposal_normals(:, :, :)
      real(dp), allocatable :: terminal_normals(:, :)
      integer :: base_samples, paired_samples, time

      if (samples < 1 .or. size(y, 1) < 1 .or. size(y, 2) < 1) then
         out%info = 1
         return
      end if
      allocate(proposal_normals(size(initial_mean), samples, size(y, 2)))
      allocate(terminal_normals(size(state_noise_loading, 2), samples))
      if (use_antithetic .and. samples > 1) then
         base_samples = samples/2
         paired_samples = 2*base_samples
         do time = 1, size(y, 2)
            call random_standard_normal_matrix( &
               proposal_normals(:, 1:base_samples, time))
            proposal_normals(:, base_samples + 1:paired_samples, time) = &
               -proposal_normals(:, 1:base_samples, time)
         end do
         call random_standard_normal_matrix( &
            terminal_normals(:, 1:base_samples))
         terminal_normals(:, base_samples + 1:paired_samples) = &
            -terminal_normals(:, 1:base_samples)
         if (paired_samples < samples) then
            do time = 1, size(y, 2)
               call random_standard_normal_matrix( &
                  proposal_normals(:, samples:samples, time))
            end do
            call random_standard_normal_matrix( &
               terminal_normals(:, samples:samples))
         end if
      else
         do time = 1, size(y, 2)
            call random_standard_normal_matrix(proposal_normals(:, :, time))
         end do
         call random_standard_normal_matrix(terminal_normals)
      end if
      out = bssm_multivariate_spdk_importance_draws(y, observation_loading, &
         transition, state_noise_loading, initial_mean, initial_covariance, &
         distribution, phi, approximation_iterations, convergence_tolerance, &
         proposal_normals, terminal_normals, offset, auxiliary, state_offset, &
         initial_mode)
   end function bssm_multivariate_spdk_importance

   pure function bssm_psi_filter_draws(y, observation_loading, transition, &
      state_noise_loading, initial_mean, initial_covariance, distribution, &
      phi, approximation_iterations, convergence_tolerance, proposal_normals, &
      terminal_normals, resampling_uniforms, offset, auxiliary, state_offset, &
      initial_mode) result(out)
      !! Run a draw-driven psi auxiliary particle filter.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: observation_loading(:, :) !! Observation loading matrix.
      real(dp), intent(in) :: transition(:, :, :) !! State transition matrix.
      real(dp), intent(in) :: state_noise_loading(:, :, :) !! State noise loading.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      integer, intent(in) :: distribution !! Probability-distribution specification.
      integer, intent(in) :: approximation_iterations !! Number of approximation iterations.
      real(dp), intent(in) :: phi !! Autoregressive or model coefficient.
      real(dp), intent(in) :: convergence_tolerance !! Convergence tolerance.
      real(dp), intent(in) :: proposal_normals(:, :, :) !! Standard-normal proposal draws.
      real(dp), intent(in) :: terminal_normals(:, :) !! Terminal normals.
      real(dp), intent(in) :: resampling_uniforms(:, :) !! Resampling uniforms.
      real(dp), intent(in), optional :: offset(:) !! Known additive offset.
      real(dp), intent(in), optional :: auxiliary(:) !! Auxiliary.
      real(dp), intent(in), optional :: state_offset(:, :) !! State offset.
      real(dp), intent(in), optional :: initial_mode(:) !! Initial mode.
      type(bssm_particle_filter_t) :: out
      type(bssm_gaussian_approximation_t) :: approximation
      real(dp), allocatable :: offset_work(:), auxiliary_work(:)
      real(dp), allocatable :: state_offset_work(:, :), parent(:, :)
      real(dp), allocatable :: log_weight(:), probability(:)
      real(dp) :: increment, signal, gaussian_log_density
      integer :: state, noise, particles, times, particle, time, info

      state = size(initial_mean)
      noise = size(state_noise_loading, 2)
      particles = size(proposal_normals, 2)
      times = size(y)
      if (state < 1 .or. noise < 1 .or. particles < 1 .or. times < 1 .or. &
         any(shape(observation_loading) /= [state, times]) .or. &
         any(shape(transition) /= [state, state, times]) .or. &
         any(shape(state_noise_loading) /= [state, noise, times]) .or. &
         any(shape(initial_covariance) /= [state, state]) .or. &
         any(shape(proposal_normals) /= [state, particles, times]) .or. &
         any(shape(terminal_normals) /= [noise, particles]) .or. &
         any(shape(resampling_uniforms) /= [particles, times]) .or. &
         .not. all(ieee_is_finite(proposal_normals)) .or. &
         .not. all(ieee_is_finite(terminal_normals)) .or. &
         .not. all(ieee_is_finite(resampling_uniforms))) then
         out%info = 1
         return
      end if
      allocate(offset_work(times), auxiliary_work(times))
      allocate(state_offset_work(state, times))
      offset_work = 0.0_dp
      auxiliary_work = 1.0_dp
      state_offset_work = 0.0_dp
      if (present(offset)) then
         if (size(offset) /= times .or. .not. all(ieee_is_finite(offset))) then
            out%info = 1
            return
         end if
         offset_work = offset
      end if
      if (present(auxiliary)) then
         if (size(auxiliary) /= times .or. &
            .not. all(ieee_is_finite(auxiliary))) then
            out%info = 1
            return
         end if
         auxiliary_work = auxiliary
      end if
      if (present(state_offset)) then
         if (any(shape(state_offset) /= [state, times]) .or. &
            .not. all(ieee_is_finite(state_offset))) then
            out%info = 1
            return
         end if
         state_offset_work = state_offset
      end if
      approximation = bssm_gaussian_approximation(y, observation_loading, &
         transition, state_noise_loading, initial_mean, initial_covariance, &
         distribution, phi, approximation_iterations, convergence_tolerance, &
         offset_work, auxiliary_work, state_offset_work, initial_mode)
      if (approximation%info /= 0) then
         out%info = 1000 + approximation%info
         return
      end if
      allocate(out%particles(state, particles, times + 1))
      allocate(out%weights(particles, times), out%ancestors(particles, times))
      allocate(out%predicted_mean(state, times + 1))
      allocate(out%filtered_mean(state, times))
      allocate(out%predicted_covariance(state, state, times + 1))
      allocate(out%filtered_covariance(state, state, times))
      allocate(parent(state, particles), log_weight(particles))
      allocate(probability(particles))
      do particle = 1, particles
         out%particles(:, particle, 1) = approximation%mode_state(:, 1) + &
            matmul(approximation%proposal_factor(:, :, 1), &
            proposal_normals(:, particle, 1))
      end do
      out%log_likelihood = approximation%corrected_log_likelihood
      do time = 1, times
         call particle_summary(out%particles(:, :, time), &
            out%predicted_mean(:, time), &
            out%predicted_covariance(:, :, time))
         if (ieee_is_finite(y(time))) then
            do particle = 1, particles
               signal = offset_work(time) + &
                  dot_product(observation_loading(:, time), &
                  out%particles(:, particle, time))
               gaussian_log_density = normal_log_density( &
                  approximation%pseudo_observation(time), signal, &
                  sqrt(approximation%observation_variance(time)))
               log_weight(particle) = bssm_observation_log_density(y(time), &
                  signal, distribution, phi, auxiliary_work(time)) - &
                  gaussian_log_density - approximation%scaling(time)
            end do
            call normalize_log_weights(log_weight, probability, increment, info)
            if (info /= 0) then
               out%info = 2
               out%log_likelihood = -huge(1.0_dp)
               return
            end if
            out%log_likelihood = out%log_likelihood + increment
         else
            probability = 1.0_dp/real(particles, dp)
         end if
         out%weights(:, time) = probability
         call weighted_particle_summary(out%particles(:, :, time), &
            probability, out%filtered_mean(:, time), &
            out%filtered_covariance(:, :, time))
         out%ancestors(:, time) = bssm_stratified_resample(probability, &
            resampling_uniforms(:, time))
         do particle = 1, particles
            parent(:, particle) = out%particles(:, &
               out%ancestors(particle, time), time)
            if (time < times) then
               out%particles(:, particle, time + 1) = &
                  approximation%mode_state(:, time + 1) + &
                  matmul(approximation%conditional_matrix(:, :, time + 1), &
                  parent(:, particle) - approximation%mode_state(:, time)) + &
                  matmul(approximation%proposal_factor(:, :, time + 1), &
                  proposal_normals(:, particle, time + 1))
            else
               out%particles(:, particle, time + 1) = &
                  state_offset_work(:, time) + &
                  matmul(transition(:, :, time), parent(:, particle)) + &
                  matmul(state_noise_loading(:, :, time), &
                  terminal_normals(:, particle))
            end if
         end do
      end do
      call particle_summary(out%particles(:, :, times + 1), &
         out%predicted_mean(:, times + 1), &
         out%predicted_covariance(:, :, times + 1))
   end function bssm_psi_filter_draws

   function bssm_psi_filter(y, observation_loading, transition, &
      state_noise_loading, initial_mean, initial_covariance, distribution, &
      phi, particles, approximation_iterations, convergence_tolerance, offset, &
      auxiliary, state_offset, initial_mode) result(out)
      !! Run a psi auxiliary particle filter using the shared random stream.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: observation_loading(:, :) !! Observation loading matrix.
      real(dp), intent(in) :: transition(:, :, :) !! State transition matrix.
      real(dp), intent(in) :: state_noise_loading(:, :, :) !! State noise loading.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      integer, intent(in) :: distribution !! Probability-distribution specification.
      integer, intent(in) :: particles !! Number of particles.
      integer, intent(in) :: approximation_iterations !! Number of approximation iterations.
      real(dp), intent(in) :: phi !! Autoregressive or model coefficient.
      real(dp), intent(in) :: convergence_tolerance !! Convergence tolerance.
      real(dp), intent(in), optional :: offset(:) !! Known additive offset.
      real(dp), intent(in), optional :: auxiliary(:) !! Auxiliary.
      real(dp), intent(in), optional :: state_offset(:, :) !! State offset.
      real(dp), intent(in), optional :: initial_mode(:) !! Initial mode.
      type(bssm_particle_filter_t) :: out
      real(dp), allocatable :: proposal_normals(:, :, :)
      real(dp), allocatable :: terminal_normals(:, :), uniforms(:, :)
      integer :: time

      if (particles < 1 .or. size(y) < 1) then
         out%info = 1
         return
      end if
      allocate(proposal_normals(size(initial_mean), particles, size(y)))
      allocate(terminal_normals(size(state_noise_loading, 2), particles))
      allocate(uniforms(particles, size(y)))
      do time = 1, size(y)
         call random_standard_normal_matrix(proposal_normals(:, :, time))
      end do
      call random_standard_normal_matrix(terminal_normals)
      call random_number(uniforms)
      out = bssm_psi_filter_draws(y, observation_loading, transition, &
         state_noise_loading, initial_mean, initial_covariance, distribution, &
         phi, approximation_iterations, convergence_tolerance, &
         proposal_normals, terminal_normals, uniforms, offset, auxiliary, &
         state_offset, initial_mode)
   end function bssm_psi_filter

   pure function bssm_spdk_importance_draws(y, observation_loading, &
      transition, state_noise_loading, initial_mean, initial_covariance, &
      distribution, phi, approximation_iterations, convergence_tolerance, &
      proposal_normals, terminal_normals, offset, auxiliary, state_offset, &
      initial_mode) result(out)
      !! Draw complete Gaussian-approximation trajectories for SPDK sampling.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: observation_loading(:, :) !! Observation loading matrix.
      real(dp), intent(in) :: transition(:, :, :) !! State transition matrix.
      real(dp), intent(in) :: state_noise_loading(:, :, :) !! State noise loading.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      integer, intent(in) :: distribution !! Probability-distribution specification.
      integer, intent(in) :: approximation_iterations !! Number of approximation iterations.
      real(dp), intent(in) :: phi !! Autoregressive or model coefficient.
      real(dp), intent(in) :: convergence_tolerance !! Convergence tolerance.
      real(dp), intent(in) :: proposal_normals(:, :, :) !! Standard-normal proposal draws.
      real(dp), intent(in) :: terminal_normals(:, :) !! Terminal normals.
      real(dp), intent(in), optional :: offset(:) !! Known additive offset.
      real(dp), intent(in), optional :: auxiliary(:) !! Auxiliary.
      real(dp), intent(in), optional :: state_offset(:, :) !! State offset.
      real(dp), intent(in), optional :: initial_mode(:) !! Initial mode.
      type(bssm_importance_sample_t) :: out
      type(bssm_gaussian_approximation_t) :: approximation
      real(dp), allocatable :: offset_work(:), auxiliary_work(:)
      real(dp), allocatable :: state_offset_work(:, :)
      real(dp) :: increment, signal, gaussian_log_density
      integer :: state, noise, samples, times, sample, time, info

      state = size(initial_mean)
      noise = size(state_noise_loading, 2)
      samples = size(proposal_normals, 2)
      times = size(y)
      if (state < 1 .or. noise < 1 .or. samples < 1 .or. times < 1 .or. &
         any(shape(observation_loading) /= [state, times]) .or. &
         any(shape(transition) /= [state, state, times]) .or. &
         any(shape(state_noise_loading) /= [state, noise, times]) .or. &
         any(shape(initial_covariance) /= [state, state]) .or. &
         any(shape(proposal_normals) /= [state, samples, times]) .or. &
         any(shape(terminal_normals) /= [noise, samples]) .or. &
         .not. all(ieee_is_finite(proposal_normals)) .or. &
         .not. all(ieee_is_finite(terminal_normals))) then
         out%info = 1
         return
      end if
      allocate(offset_work(times), auxiliary_work(times))
      allocate(state_offset_work(state, times))
      offset_work = 0.0_dp
      auxiliary_work = 1.0_dp
      state_offset_work = 0.0_dp
      if (present(offset)) then
         if (size(offset) /= times .or. .not. all(ieee_is_finite(offset))) then
            out%info = 1
            return
         end if
         offset_work = offset
      end if
      if (present(auxiliary)) then
         if (size(auxiliary) /= times .or. &
            .not. all(ieee_is_finite(auxiliary))) then
            out%info = 1
            return
         end if
         auxiliary_work = auxiliary
      end if
      if (present(state_offset)) then
         if (any(shape(state_offset) /= [state, times]) .or. &
            .not. all(ieee_is_finite(state_offset))) then
            out%info = 1
            return
         end if
         state_offset_work = state_offset
      end if
      approximation = bssm_gaussian_approximation(y, observation_loading, &
         transition, state_noise_loading, initial_mean, initial_covariance, &
         distribution, phi, approximation_iterations, convergence_tolerance, &
         offset_work, auxiliary_work, state_offset_work, initial_mode)
      if (approximation%info /= 0) then
         out%info = 1000 + approximation%info
         return
      end if
      allocate(out%trajectories(state, times + 1, samples))
      allocate(out%log_weight(samples), out%weight(samples))
      allocate(out%mean(state, times + 1))
      allocate(out%covariance(state, state, times + 1))
      do sample = 1, samples
         out%trajectories(:, 1, sample) = approximation%mode_state(:, 1) + &
            matmul(approximation%proposal_factor(:, :, 1), &
            proposal_normals(:, sample, 1))
         do time = 2, times
            out%trajectories(:, time, sample) = &
               approximation%mode_state(:, time) + &
               matmul(approximation%conditional_matrix(:, :, time), &
               out%trajectories(:, time - 1, sample) - &
               approximation%mode_state(:, time - 1)) + &
               matmul(approximation%proposal_factor(:, :, time), &
               proposal_normals(:, sample, time))
         end do
         out%trajectories(:, times + 1, sample) = &
            state_offset_work(:, times) + &
            matmul(transition(:, :, times), &
            out%trajectories(:, times, sample)) + &
            matmul(state_noise_loading(:, :, times), &
            terminal_normals(:, sample))
      end do
      out%log_weight = 0.0_dp
      do sample = 1, samples
         do time = 1, times
            if (.not. ieee_is_finite(y(time))) cycle
            signal = offset_work(time) + &
               dot_product(observation_loading(:, time), &
               out%trajectories(:, time, sample))
            gaussian_log_density = normal_log_density( &
               approximation%pseudo_observation(time), signal, &
               sqrt(approximation%observation_variance(time)))
            out%log_weight(sample) = out%log_weight(sample) + &
               bssm_observation_log_density(y(time), signal, distribution, &
               phi, auxiliary_work(time)) - gaussian_log_density - &
               approximation%scaling(time)
         end do
      end do
      call normalize_log_weights(out%log_weight, out%weight, increment, info)
      if (info /= 0) then
         out%info = 2
         out%log_likelihood = -huge(1.0_dp)
         return
      end if
      out%log_likelihood = approximation%corrected_log_likelihood + increment
      out%effective_sample_size = 1.0_dp/sum(out%weight**2)
      do time = 1, times + 1
         call weighted_particle_summary(out%trajectories(:, time, :), &
            out%weight, out%mean(:, time), out%covariance(:, :, time))
      end do
   end function bssm_spdk_importance_draws

   function bssm_spdk_importance(y, observation_loading, transition, &
      state_noise_loading, initial_mean, initial_covariance, distribution, &
      phi, samples, approximation_iterations, convergence_tolerance, &
      use_antithetic, offset, auxiliary, state_offset, initial_mode) result(out)
      !! Run SPDK importance sampling using the shared random stream.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: observation_loading(:, :) !! Observation loading matrix.
      real(dp), intent(in) :: transition(:, :, :) !! State transition matrix.
      real(dp), intent(in) :: state_noise_loading(:, :, :) !! State noise loading.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      integer, intent(in) :: distribution !! Probability-distribution specification.
      integer, intent(in) :: samples !! Samples.
      integer, intent(in) :: approximation_iterations !! Number of approximation iterations.
      real(dp), intent(in) :: phi !! Autoregressive or model coefficient.
      real(dp), intent(in) :: convergence_tolerance !! Convergence tolerance.
      logical, intent(in) :: use_antithetic !! Whether to use the antithetic.
      real(dp), intent(in), optional :: offset(:) !! Known additive offset.
      real(dp), intent(in), optional :: auxiliary(:) !! Auxiliary.
      real(dp), intent(in), optional :: state_offset(:, :) !! State offset.
      real(dp), intent(in), optional :: initial_mode(:) !! Initial mode.
      type(bssm_importance_sample_t) :: out
      real(dp), allocatable :: proposal_normals(:, :, :)
      real(dp), allocatable :: terminal_normals(:, :)
      integer :: base_samples, paired_samples, time

      if (samples < 1 .or. size(y) < 1) then
         out%info = 1
         return
      end if
      allocate(proposal_normals(size(initial_mean), samples, size(y)))
      allocate(terminal_normals(size(state_noise_loading, 2), samples))
      if (use_antithetic .and. samples > 1) then
         base_samples = samples/2
         paired_samples = 2*base_samples
         do time = 1, size(y)
            call random_standard_normal_matrix( &
               proposal_normals(:, 1:base_samples, time))
            proposal_normals(:, base_samples + 1:paired_samples, time) = &
               -proposal_normals(:, 1:base_samples, time)
         end do
         call random_standard_normal_matrix( &
            terminal_normals(:, 1:base_samples))
         terminal_normals(:, base_samples + 1:paired_samples) = &
            -terminal_normals(:, 1:base_samples)
         if (paired_samples < samples) then
            do time = 1, size(y)
               call random_standard_normal_matrix( &
                  proposal_normals(:, samples:samples, time))
            end do
            call random_standard_normal_matrix( &
               terminal_normals(:, samples:samples))
         end if
      else
         do time = 1, size(y)
            call random_standard_normal_matrix(proposal_normals(:, :, time))
         end do
         call random_standard_normal_matrix(terminal_normals)
      end if
      out = bssm_spdk_importance_draws(y, observation_loading, transition, &
         state_noise_loading, initial_mean, initial_covariance, distribution, &
         phi, approximation_iterations, convergence_tolerance, &
         proposal_normals, terminal_normals, offset, auxiliary, state_offset, &
         initial_mode)
   end function bssm_spdk_importance

   pure function bssm_ekpf_draws(y, initial_mean, initial_covariance, &
      parameters, noise_dimension, observation, transition, initial_normals, &
      proposal_normals, resampling_uniforms, iekf_iterations, &
      convergence_tolerance) result(out)
      !! Run a draw-driven extended Kalman particle filter.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: parameters(:) !! Model parameter values.
      real(dp), intent(in) :: initial_normals(:, :) !! Initial normals.
      integer, intent(in) :: noise_dimension !! Noise dimension.
      procedure(bssm_gaussian_observation_t) :: observation !! Observed value or vector.
      procedure(bssm_nonlinear_transition_t) :: transition !! State transition matrix.
      real(dp), intent(in) :: proposal_normals(:, :, :) !! Standard-normal proposal draws.
      real(dp), intent(in) :: resampling_uniforms(:, :) !! Resampling uniforms.
      integer, intent(in), optional :: iekf_iterations !! Number of iterated extended Kalman filter iterations.
      real(dp), intent(in), optional :: convergence_tolerance !! Convergence tolerance.
      type(bssm_particle_filter_t) :: out
      real(dp), allocatable :: prior_mean(:, :), proposal_mean(:, :)
      real(dp), allocatable :: prior_covariance(:, :, :)
      real(dp), allocatable :: proposal_covariance(:, :, :)
      real(dp), allocatable :: factor(:, :), log_weight(:), probability(:)
      real(dp), allocatable :: parent(:, :), noise_loading(:, :)
      real(dp) :: observation_mean, observation_sd, increment
      real(dp), allocatable :: jacobian(:)
      logical :: has_next_observation
      integer :: state, particles, times, particle, time, info, iteration_limit
      real(dp) :: tolerance

      state = size(initial_mean)
      particles = size(initial_normals, 2)
      times = size(y)
      iteration_limit = 0
      if (present(iekf_iterations)) iteration_limit = iekf_iterations
      tolerance = 1.0e-4_dp
      if (present(convergence_tolerance)) tolerance = convergence_tolerance
      if (state < 1 .or. particles < 1 .or. times < 1 .or. &
         noise_dimension < 1 .or. iteration_limit < 0 .or. tolerance < 0.0_dp &
         .or. size(initial_normals, 1) /= state .or. &
         any(shape(initial_covariance) /= [state, state]) .or. &
         any(shape(proposal_normals) /= [state, particles, times]) .or. &
         any(shape(resampling_uniforms) /= [particles, times]) .or. &
         .not. all(ieee_is_finite(initial_mean)) .or. &
         .not. all(ieee_is_finite(initial_covariance)) .or. &
         .not. all(ieee_is_finite(parameters)) .or. &
         .not. all(ieee_is_finite(initial_normals)) .or. &
         .not. all(ieee_is_finite(proposal_normals)) .or. &
         .not. all(ieee_is_finite(resampling_uniforms))) then
         out%info = 1
         return
      end if
      allocate(out%particles(state, particles, times + 1))
      allocate(out%weights(particles, times), out%ancestors(particles, times))
      allocate(out%predicted_mean(state, times + 1))
      allocate(out%filtered_mean(state, times))
      allocate(out%predicted_covariance(state, state, times + 1))
      allocate(out%filtered_covariance(state, state, times))
      allocate(prior_mean(state, particles), proposal_mean(state, particles))
      allocate(prior_covariance(state, state, particles))
      allocate(proposal_covariance(state, state, particles))
      allocate(log_weight(particles), probability(particles))
      allocate(parent(state, particles), noise_loading(state, noise_dimension))
      allocate(jacobian(state))
      do particle = 1, particles
         prior_mean(:, particle) = initial_mean
         prior_covariance(:, :, particle) = initial_covariance
         if (ieee_is_finite(y(1))) then
            call ekf_proposal_update(1, y(1), initial_mean, &
               initial_covariance, parameters, observation, &
               proposal_mean(:, particle), proposal_covariance(:, :, particle), &
               info, iteration_limit, tolerance)
            if (info /= 0) then
               out%info = 2
               return
            end if
         else
            proposal_mean(:, particle) = initial_mean
            proposal_covariance(:, :, particle) = initial_covariance
         end if
         call positive_semidefinite_factor(proposal_covariance(:, :, particle), &
            factor, info)
         if (info /= 0) then
            out%info = 2
            return
         end if
         out%particles(:, particle, 1) = proposal_mean(:, particle) + &
            matmul(factor, initial_normals(:, particle))
      end do
      do time = 1, times
         call particle_summary(out%particles(:, :, time), &
            out%predicted_mean(:, time), &
            out%predicted_covariance(:, :, time))
         if (ieee_is_finite(y(time))) then
            do particle = 1, particles
               call observation(time, out%particles(:, particle, time), &
                  parameters, observation_mean, jacobian, observation_sd)
               if (.not. ieee_is_finite(observation_mean) .or. &
                  .not. ieee_is_finite(observation_sd) .or. &
                  observation_sd <= 0.0_dp .or. &
                  .not. all(ieee_is_finite(jacobian))) then
                  out%info = 3
                  return
               end if
               log_weight(particle) = normal_log_density(y(time), &
                  observation_mean, observation_sd) + &
                  multivariate_normal_log_density( &
                  out%particles(:, particle, time), prior_mean(:, particle), &
                  prior_covariance(:, :, particle)) - &
                  multivariate_normal_log_density( &
                  out%particles(:, particle, time), proposal_mean(:, particle), &
                  proposal_covariance(:, :, particle))
            end do
            call normalize_log_weights(log_weight, probability, increment, info)
            if (info /= 0) then
               out%info = 4
               out%log_likelihood = -huge(1.0_dp)
               return
            end if
            out%log_likelihood = out%log_likelihood + increment
         else
            probability = 1.0_dp/real(particles, dp)
         end if
         out%weights(:, time) = probability
         call weighted_particle_summary(out%particles(:, :, time), probability, &
            out%filtered_mean(:, time), &
            out%filtered_covariance(:, :, time))
         out%ancestors(:, time) = bssm_stratified_resample(probability, &
            resampling_uniforms(:, time))
         do particle = 1, particles
            parent(:, particle) = &
               out%particles(:, out%ancestors(particle, time), time)
            call transition(time, parent(:, particle), parameters, &
               prior_mean(:, particle), noise_loading)
            if (.not. all(ieee_is_finite(prior_mean(:, particle))) .or. &
               .not. all(ieee_is_finite(noise_loading))) then
               out%info = 5
               return
            end if
            prior_covariance(:, :, particle) = &
               matmul(noise_loading, transpose(noise_loading))
            has_next_observation = .false.
            if (time < times) then
               has_next_observation = ieee_is_finite(y(time + 1))
            end if
            if (has_next_observation) then
               call ekf_proposal_update(time + 1, y(time + 1), &
                  prior_mean(:, particle), prior_covariance(:, :, particle), &
                  parameters, observation, proposal_mean(:, particle), &
                  proposal_covariance(:, :, particle), info, &
                  iteration_limit, tolerance)
               if (info /= 0) then
                  out%info = 6
                  return
               end if
            else
               proposal_mean(:, particle) = prior_mean(:, particle)
               proposal_covariance(:, :, particle) = &
                  prior_covariance(:, :, particle)
            end if
            call positive_semidefinite_factor( &
               proposal_covariance(:, :, particle), factor, info)
            if (info /= 0) then
               out%info = 6
               return
            end if
            out%particles(:, particle, time + 1) = &
               proposal_mean(:, particle) + &
               matmul(factor, proposal_normals(:, particle, time))
         end do
      end do
      call particle_summary(out%particles(:, :, times + 1), &
         out%predicted_mean(:, times + 1), &
         out%predicted_covariance(:, :, times + 1))
   end function bssm_ekpf_draws

   function bssm_ekpf(y, initial_mean, initial_covariance, parameters, &
      noise_dimension, observation, transition, particles, iekf_iterations, &
      convergence_tolerance) result(out)
      !! Run an extended Kalman particle filter using the shared random stream.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: parameters(:) !! Model parameter values.
      integer, intent(in) :: noise_dimension !! Noise dimension.
      integer, intent(in) :: particles !! Number of particles.
      integer, intent(in), optional :: iekf_iterations !! Number of iterated extended Kalman filter iterations.
      real(dp), intent(in), optional :: convergence_tolerance !! Convergence tolerance.
      procedure(bssm_gaussian_observation_t) :: observation !! Observed value or vector.
      procedure(bssm_nonlinear_transition_t) :: transition !! State transition matrix.
      type(bssm_particle_filter_t) :: out
      real(dp), allocatable :: initial_normals(:, :)
      real(dp), allocatable :: proposal_normals(:, :, :), uniforms(:, :)
      integer :: time, iteration_limit
      real(dp) :: tolerance

      iteration_limit = 0
      if (present(iekf_iterations)) iteration_limit = iekf_iterations
      tolerance = 1.0e-4_dp
      if (present(convergence_tolerance)) tolerance = convergence_tolerance
      if (particles < 1 .or. noise_dimension < 1 .or. size(y) < 1 .or. &
         iteration_limit < 0 .or. tolerance < 0.0_dp) then
         out%info = 1
         return
      end if
      allocate(initial_normals(size(initial_mean), particles))
      allocate(proposal_normals(size(initial_mean), particles, size(y)))
      allocate(uniforms(particles, size(y)))
      call random_standard_normal_matrix(initial_normals)
      do time = 1, size(y)
         call random_standard_normal_matrix(proposal_normals(:, :, time))
      end do
      call random_number(uniforms)
      out = bssm_ekpf_draws(y, initial_mean, initial_covariance, parameters, &
         noise_dimension, observation, transition, initial_normals, &
         proposal_normals, uniforms, iteration_limit, tolerance)
   end function bssm_ekpf

   pure function bssm_multivariate_ekpf_draws(y, initial_mean, &
      initial_covariance, parameters, observation_noise_dimension, &
      state_noise_dimension, observation, transition, initial_normals, &
      proposal_normals, resampling_uniforms, iekf_iterations, &
      convergence_tolerance) result(out)
      !! Run a supplied-draw multivariate extended Kalman particle filter.
      real(dp), intent(in) :: y(:, :) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: parameters(:) !! Model parameter values.
      integer, intent(in) :: observation_noise_dimension !! Observation noise dimension.
      integer, intent(in) :: state_noise_dimension !! State noise dimension.
      procedure(bssm_multivariate_gaussian_observation_t) :: observation !! Observed value or vector.
      procedure(bssm_nonlinear_transition_t) :: transition !! State transition matrix.
      real(dp), intent(in) :: initial_normals(:, :) !! Initial normals.
      real(dp), intent(in) :: proposal_normals(:, :, :) !! Standard-normal proposal draws.
      real(dp), intent(in) :: resampling_uniforms(:, :) !! Resampling uniforms.
      integer, intent(in), optional :: iekf_iterations !! Number of iterated extended Kalman filter iterations.
      real(dp), intent(in), optional :: convergence_tolerance !! Convergence tolerance.
      type(bssm_particle_filter_t) :: out
      real(dp), allocatable :: prior_mean(:, :), proposal_mean(:, :)
      real(dp), allocatable :: prior_covariance(:, :, :)
      real(dp), allocatable :: proposal_covariance(:, :, :)
      real(dp), allocatable :: factor(:, :), log_weight(:), probability(:)
      real(dp), allocatable :: parent(:, :), noise_loading(:, :)
      real(dp), allocatable :: innovation(:), innovation_covariance(:, :)
      real(dp) :: increment, update_log_density, tolerance
      logical :: has_next_observation
      integer :: series, state, particles, times, particle, time, info
      integer :: iteration_limit, iterations_used

      series = size(y, 1)
      times = size(y, 2)
      state = size(initial_mean)
      particles = size(initial_normals, 2)
      iteration_limit = 0
      if (present(iekf_iterations)) iteration_limit = iekf_iterations
      tolerance = 1.0e-4_dp
      if (present(convergence_tolerance)) tolerance = convergence_tolerance
      if (series < 1 .or. state < 1 .or. particles < 1 .or. times < 1 .or. &
         observation_noise_dimension < 1 .or. state_noise_dimension < 1 .or. &
         iteration_limit < 0 .or. tolerance < 0.0_dp .or. &
         size(initial_normals, 1) /= state .or. &
         any(shape(initial_covariance) /= [state, state]) .or. &
         any(shape(proposal_normals) /= [state, particles, times]) .or. &
         any(shape(resampling_uniforms) /= [particles, times]) .or. &
         .not. all(ieee_is_finite(initial_mean)) .or. &
         .not. all(ieee_is_finite(initial_covariance)) .or. &
         .not. all(ieee_is_finite(parameters)) .or. &
         .not. all(ieee_is_finite(initial_normals)) .or. &
         .not. all(ieee_is_finite(proposal_normals)) .or. &
         .not. all(ieee_is_finite(resampling_uniforms))) then
         out%info = 1
         return
      end if
      allocate(out%particles(state, particles, times + 1))
      allocate(out%weights(particles, times), out%ancestors(particles, times))
      allocate(out%predicted_mean(state, times + 1))
      allocate(out%filtered_mean(state, times))
      allocate(out%predicted_covariance(state, state, times + 1))
      allocate(out%filtered_covariance(state, state, times))
      allocate(prior_mean(state, particles), proposal_mean(state, particles))
      allocate(prior_covariance(state, state, particles))
      allocate(proposal_covariance(state, state, particles))
      allocate(log_weight(particles), probability(particles))
      allocate(parent(state, particles))
      allocate(noise_loading(state, state_noise_dimension))
      allocate(innovation(series), innovation_covariance(series, series))
      do particle = 1, particles
         prior_mean(:, particle) = initial_mean
         prior_covariance(:, :, particle) = initial_covariance
         if (any(ieee_is_finite(y(:, 1)))) then
            call multivariate_ekf_proposal_update(1, y(:, 1), initial_mean, &
               initial_covariance, parameters, observation_noise_dimension, &
               observation, proposal_mean(:, particle), &
               proposal_covariance(:, :, particle), innovation, &
               innovation_covariance, update_log_density, info, &
               iteration_limit, tolerance, iterations_used)
            if (info /= 0) then
               out%info = 2
               return
            end if
         else
            proposal_mean(:, particle) = initial_mean
            proposal_covariance(:, :, particle) = initial_covariance
         end if
         call positive_semidefinite_factor(proposal_covariance(:, :, particle), &
            factor, info)
         if (info /= 0) then
            out%info = 2
            return
         end if
         out%particles(:, particle, 1) = proposal_mean(:, particle) + &
            matmul(factor, initial_normals(:, particle))
      end do
      do time = 1, times
         call particle_summary(out%particles(:, :, time), &
            out%predicted_mean(:, time), &
            out%predicted_covariance(:, :, time))
         if (any(ieee_is_finite(y(:, time)))) then
            do particle = 1, particles
               log_weight(particle) = &
                  multivariate_gaussian_observation_log_density(time, &
                  y(:, time), out%particles(:, particle, time), parameters, &
                  observation_noise_dimension, observation) + &
                  multivariate_normal_log_density( &
                  out%particles(:, particle, time), prior_mean(:, particle), &
                  prior_covariance(:, :, particle)) - &
                  multivariate_normal_log_density( &
                  out%particles(:, particle, time), proposal_mean(:, particle), &
                  proposal_covariance(:, :, particle))
            end do
            call normalize_log_weights(log_weight, probability, increment, &
               info)
            if (info /= 0) then
               out%info = 4
               out%log_likelihood = -huge(1.0_dp)
               return
            end if
            out%log_likelihood = out%log_likelihood + increment
         else
            probability = 1.0_dp/real(particles, dp)
         end if
         out%weights(:, time) = probability
         call weighted_particle_summary(out%particles(:, :, time), &
            probability, out%filtered_mean(:, time), &
            out%filtered_covariance(:, :, time))
         out%ancestors(:, time) = bssm_stratified_resample(probability, &
            resampling_uniforms(:, time))
         do particle = 1, particles
            parent(:, particle) = &
               out%particles(:, out%ancestors(particle, time), time)
            call transition(time, parent(:, particle), parameters, &
               prior_mean(:, particle), noise_loading)
            if (.not. all(ieee_is_finite(prior_mean(:, particle))) .or. &
               .not. all(ieee_is_finite(noise_loading))) then
               out%info = 5
               return
            end if
            prior_covariance(:, :, particle) = &
               matmul(noise_loading, transpose(noise_loading))
            has_next_observation = .false.
            if (time < times) then
               has_next_observation = any(ieee_is_finite(y(:, time + 1)))
            end if
            if (has_next_observation) then
               call multivariate_ekf_proposal_update(time + 1, &
                  y(:, time + 1), prior_mean(:, particle), &
                  prior_covariance(:, :, particle), parameters, &
                  observation_noise_dimension, observation, &
                  proposal_mean(:, particle), &
                  proposal_covariance(:, :, particle), innovation, &
                  innovation_covariance, update_log_density, info, &
                  iteration_limit, tolerance, iterations_used)
               if (info /= 0) then
                  out%info = 6
                  return
               end if
            else
               proposal_mean(:, particle) = prior_mean(:, particle)
               proposal_covariance(:, :, particle) = &
                  prior_covariance(:, :, particle)
            end if
            call positive_semidefinite_factor( &
               proposal_covariance(:, :, particle), factor, info)
            if (info /= 0) then
               out%info = 6
               return
            end if
            out%particles(:, particle, time + 1) = &
               proposal_mean(:, particle) + &
               matmul(factor, proposal_normals(:, particle, time))
         end do
      end do
      call particle_summary(out%particles(:, :, times + 1), &
         out%predicted_mean(:, times + 1), &
         out%predicted_covariance(:, :, times + 1))
   end function bssm_multivariate_ekpf_draws

   function bssm_multivariate_ekpf(y, initial_mean, initial_covariance, &
      parameters, observation_noise_dimension, state_noise_dimension, &
      observation, transition, particles, iekf_iterations, &
      convergence_tolerance) result(out)
      !! Run a multivariate EKF-proposal filter from the shared random stream.
      real(dp), intent(in) :: y(:, :) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: parameters(:) !! Model parameter values.
      integer, intent(in) :: observation_noise_dimension !! Observation noise dimension.
      integer, intent(in) :: state_noise_dimension !! State noise dimension.
      integer, intent(in) :: particles !! Number of particles.
      procedure(bssm_multivariate_gaussian_observation_t) :: observation !! Observed value or vector.
      procedure(bssm_nonlinear_transition_t) :: transition !! State transition matrix.
      integer, intent(in), optional :: iekf_iterations !! Number of iterated extended Kalman filter iterations.
      real(dp), intent(in), optional :: convergence_tolerance !! Convergence tolerance.
      type(bssm_particle_filter_t) :: out
      real(dp), allocatable :: initial_normals(:, :)
      real(dp), allocatable :: proposal_normals(:, :, :), uniforms(:, :)
      integer :: time, iteration_limit
      real(dp) :: tolerance

      iteration_limit = 0
      if (present(iekf_iterations)) iteration_limit = iekf_iterations
      tolerance = 1.0e-4_dp
      if (present(convergence_tolerance)) tolerance = convergence_tolerance
      if (particles < 1 .or. observation_noise_dimension < 1 .or. &
         state_noise_dimension < 1 .or. size(y, 1) < 1 .or. &
         size(y, 2) < 1 .or. iteration_limit < 0 .or. tolerance < 0.0_dp) then
         out%info = 1
         return
      end if
      allocate(initial_normals(size(initial_mean), particles))
      allocate(proposal_normals(size(initial_mean), particles, size(y, 2)))
      allocate(uniforms(particles, size(y, 2)))
      call random_standard_normal_matrix(initial_normals)
      do time = 1, size(y, 2)
         call random_standard_normal_matrix(proposal_normals(:, :, time))
      end do
      call random_number(uniforms)
      out = bssm_multivariate_ekpf_draws(y, initial_mean, initial_covariance, &
         parameters, observation_noise_dimension, state_noise_dimension, &
         observation, transition, initial_normals, proposal_normals, uniforms, &
         iteration_limit, tolerance)
   end function bssm_multivariate_ekpf

   pure function bssm_nonlinear_ekpf_pmmh_draws(y, initial_mean, &
      initial_covariance, initial_parameters, noise_dimension, observation, &
      transition, prior, iekf_iterations, convergence_tolerance, &
      proposal_factor, parameter_normals, acceptance_uniforms, &
      initial_normals, state_normals, resampling_uniforms, &
      target_acceptance, adaptation_exponent) result(out)
      !! Run draw-driven PMMH with an EKF-proposal particle filter.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: initial_parameters(:) !! Initial parameter values.
      real(dp), intent(in) :: convergence_tolerance !! Convergence tolerance.
      real(dp), intent(in) :: proposal_factor(:, :) !! Proposal factor.
      real(dp), intent(in) :: parameter_normals(:, :) !! Parameter normals.
      real(dp), intent(in) :: acceptance_uniforms(:) !! Acceptance uniforms.
      real(dp), intent(in) :: initial_normals(:, :, :) !! Initial normals.
      real(dp), intent(in) :: state_normals(:, :, :, :) !! State normals.
      real(dp), intent(in) :: resampling_uniforms(:, :, :) !! Resampling uniforms.
      integer, intent(in) :: noise_dimension !! Noise dimension.
      integer, intent(in) :: iekf_iterations !! Number of iterated extended Kalman filter iterations.
      procedure(bssm_gaussian_observation_t) :: observation !! Observed value or vector.
      procedure(bssm_nonlinear_transition_t) :: transition !! State transition matrix.
      procedure(bssm_parameter_log_density_t) :: prior !! Prior-distribution specification.
      real(dp), intent(in), optional :: target_acceptance !! Target acceptance.
      real(dp), intent(in), optional :: adaptation_exponent !! Adaptation exponent.
      type(bssm_mcmc_t) :: out
      integer :: state, particles, times, iterations
      real(dp) :: exponent

      state = size(initial_mean)
      particles = size(initial_normals, 2)
      times = size(y)
      iterations = size(parameter_normals, 2)
      if (state < 1 .or. particles < 1 .or. times < 1 .or. &
         noise_dimension < 1 .or. iekf_iterations < 0 .or. &
         convergence_tolerance < 0.0_dp .or. iterations < 1 .or. &
         any(shape(initial_covariance) /= [state, state]) .or. &
         any(shape(initial_normals) /= &
         [state, particles, iterations + 1]) .or. &
         any(shape(state_normals) /= &
         [state, particles, times, iterations + 1]) .or. &
         any(shape(resampling_uniforms) /= &
         [particles, times, iterations + 1]) .or. &
         .not. all(ieee_is_finite(initial_normals)) .or. &
         .not. all(ieee_is_finite(state_normals)) .or. &
         .not. all(ieee_is_finite(resampling_uniforms)) .or. &
         (present(adaptation_exponent) .and. &
         .not. present(target_acceptance))) then
         out%info = 1
         return
      end if
      exponent = 0.6_dp
      if (present(adaptation_exponent)) exponent = adaptation_exponent
      if (present(target_acceptance)) then
         out = bssm_pmmh_kernel_draws(initial_parameters, prior, estimator, &
            proposal_factor, parameter_normals, acceptance_uniforms, &
            target_acceptance, exponent)
      else
         out = bssm_pmmh_kernel_draws(initial_parameters, prior, estimator, &
            proposal_factor, parameter_normals, acceptance_uniforms)
      end if

   contains

      pure function estimator(parameters, draw_index) result(estimate)
         !! Estimate the EKF-proposal particle-filter log likelihood.
         real(dp), intent(in) :: parameters(:) !! Model parameter values.
         integer, intent(in) :: draw_index !! Index of draw.
         type(bssm_likelihood_estimate_t) :: estimate
         type(bssm_particle_filter_t) :: fit

         fit = bssm_ekpf_draws(y, initial_mean, initial_covariance, parameters, &
            noise_dimension, observation, transition, &
            initial_normals(:, :, draw_index), &
            state_normals(:, :, :, draw_index), &
            resampling_uniforms(:, :, draw_index), iekf_iterations, &
            convergence_tolerance)
         estimate%value = fit%log_likelihood
         estimate%info = fit%info
      end function estimator

   end function bssm_nonlinear_ekpf_pmmh_draws

   function bssm_nonlinear_ekpf_pmmh(y, initial_mean, initial_covariance, &
      initial_parameters, noise_dimension, observation, transition, prior, &
      iekf_iterations, convergence_tolerance, proposal_factor, particles, &
      iterations, target_acceptance, adaptation_exponent) result(out)
      !! Run EKF-proposal particle-filter PMMH using shared randomness.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: initial_parameters(:) !! Initial parameter values.
      real(dp), intent(in) :: convergence_tolerance !! Convergence tolerance.
      real(dp), intent(in) :: proposal_factor(:, :) !! Proposal factor.
      integer, intent(in) :: noise_dimension !! Noise dimension.
      integer, intent(in) :: iekf_iterations !! Number of iterated extended Kalman filter iterations.
      integer, intent(in) :: particles !! Number of particles.
      integer, intent(in) :: iterations !! Number of algorithm iterations.
      procedure(bssm_gaussian_observation_t) :: observation !! Observed value or vector.
      procedure(bssm_nonlinear_transition_t) :: transition !! State transition matrix.
      procedure(bssm_parameter_log_density_t) :: prior !! Prior-distribution specification.
      real(dp), intent(in), optional :: target_acceptance !! Target acceptance.
      real(dp), intent(in), optional :: adaptation_exponent !! Adaptation exponent.
      type(bssm_mcmc_t) :: out
      real(dp), allocatable :: parameter_normals(:, :), acceptance_uniforms(:)
      real(dp), allocatable :: initial_normals(:, :, :)
      real(dp), allocatable :: state_normals(:, :, :, :)
      real(dp), allocatable :: resampling_uniforms(:, :, :)
      real(dp) :: exponent
      integer :: draw_index, time

      if (particles < 1 .or. iterations < 1 .or. noise_dimension < 1 .or. &
         iekf_iterations < 0 .or. convergence_tolerance < 0.0_dp .or. &
         size(y) < 1 .or. size(initial_mean) < 1 .or. &
         size(initial_parameters) < 1 .or. &
         (present(adaptation_exponent) .and. &
         .not. present(target_acceptance))) then
         out%info = 1
         return
      end if
      allocate(parameter_normals(size(initial_parameters), iterations))
      allocate(acceptance_uniforms(iterations))
      allocate(initial_normals(size(initial_mean), particles, iterations + 1))
      allocate(state_normals(size(initial_mean), particles, size(y), &
         iterations + 1))
      allocate(resampling_uniforms(particles, size(y), iterations + 1))
      call random_standard_normal_matrix(parameter_normals)
      call random_number(acceptance_uniforms)
      call random_number(resampling_uniforms)
      do draw_index = 1, iterations + 1
         call random_standard_normal_matrix(initial_normals(:, :, draw_index))
         do time = 1, size(y)
            call random_standard_normal_matrix( &
               state_normals(:, :, time, draw_index))
         end do
      end do
      exponent = 0.6_dp
      if (present(adaptation_exponent)) exponent = adaptation_exponent
      if (present(target_acceptance)) then
         out = bssm_nonlinear_ekpf_pmmh_draws(y, initial_mean, &
            initial_covariance, initial_parameters, noise_dimension, &
            observation, transition, prior, iekf_iterations, &
            convergence_tolerance, proposal_factor, parameter_normals, &
            acceptance_uniforms, initial_normals, state_normals, &
            resampling_uniforms, target_acceptance, exponent)
      else
         out = bssm_nonlinear_ekpf_pmmh_draws(y, initial_mean, &
            initial_covariance, initial_parameters, noise_dimension, &
            observation, transition, prior, iekf_iterations, &
            convergence_tolerance, proposal_factor, parameter_normals, &
            acceptance_uniforms, initial_normals, state_normals, &
            resampling_uniforms)
      end if
   end function bssm_nonlinear_ekpf_pmmh

   pure function bssm_multivariate_ekpf_pmmh_draws(y, initial_mean, &
      initial_covariance, initial_parameters, observation_noise_dimension, &
      state_noise_dimension, observation, transition, prior, iekf_iterations, &
      convergence_tolerance, proposal_factor, parameter_normals, &
      acceptance_uniforms, initial_normals, state_normals, &
      resampling_uniforms, target_acceptance, adaptation_exponent) result(out)
      !! Run draw-driven multivariate EKF-proposal particle PMMH.
      real(dp), intent(in) :: y(:, :) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: initial_parameters(:) !! Initial parameter values.
      real(dp), intent(in) :: convergence_tolerance !! Convergence tolerance.
      real(dp), intent(in) :: proposal_factor(:, :) !! Proposal factor.
      real(dp), intent(in) :: parameter_normals(:, :) !! Parameter normals.
      real(dp), intent(in) :: acceptance_uniforms(:) !! Acceptance uniforms.
      real(dp), intent(in) :: initial_normals(:, :, :) !! Initial normals.
      real(dp), intent(in) :: state_normals(:, :, :, :) !! State normals.
      real(dp), intent(in) :: resampling_uniforms(:, :, :) !! Resampling uniforms.
      integer, intent(in) :: observation_noise_dimension !! Observation noise dimension.
      integer, intent(in) :: state_noise_dimension !! State noise dimension.
      integer, intent(in) :: iekf_iterations !! Number of iterated extended Kalman filter iterations.
      procedure(bssm_multivariate_gaussian_observation_t) :: observation !! Observed value or vector.
      procedure(bssm_nonlinear_transition_t) :: transition !! State transition matrix.
      procedure(bssm_parameter_log_density_t) :: prior !! Prior-distribution specification.
      real(dp), intent(in), optional :: target_acceptance !! Target acceptance.
      real(dp), intent(in), optional :: adaptation_exponent !! Adaptation exponent.
      type(bssm_mcmc_t) :: out
      real(dp) :: exponent
      integer :: state, particles, times, iterations

      state = size(initial_mean)
      particles = size(initial_normals, 2)
      times = size(y, 2)
      iterations = size(parameter_normals, 2)
      if (size(y, 1) < 1 .or. state < 1 .or. particles < 1 .or. &
         times < 1 .or. iterations < 1 .or. &
         observation_noise_dimension < 1 .or. state_noise_dimension < 1 .or. &
         iekf_iterations < 0 .or. convergence_tolerance < 0.0_dp .or. &
         any(shape(initial_covariance) /= [state, state]) .or. &
         any(shape(initial_normals) /= &
         [state, particles, iterations + 1]) .or. &
         any(shape(state_normals) /= &
         [state, particles, times, iterations + 1]) .or. &
         any(shape(resampling_uniforms) /= &
         [particles, times, iterations + 1]) .or. &
         .not. all(ieee_is_finite(initial_normals)) .or. &
         .not. all(ieee_is_finite(state_normals)) .or. &
         .not. all(ieee_is_finite(resampling_uniforms)) .or. &
         (present(adaptation_exponent) .and. &
         .not. present(target_acceptance))) then
         out%info = 1
         return
      end if
      exponent = 0.6_dp
      if (present(adaptation_exponent)) exponent = adaptation_exponent
      if (present(target_acceptance)) then
         out = bssm_pmmh_kernel_draws(initial_parameters, prior, estimator, &
            proposal_factor, parameter_normals, acceptance_uniforms, &
            target_acceptance, exponent)
      else
         out = bssm_pmmh_kernel_draws(initial_parameters, prior, estimator, &
            proposal_factor, parameter_normals, acceptance_uniforms)
      end if

   contains

      pure function estimator(parameters, draw_index) result(estimate)
         !! Estimate the multivariate EKF-proposal particle log likelihood.
         real(dp), intent(in) :: parameters(:) !! Model parameter values.
         integer, intent(in) :: draw_index !! Index of draw.
         type(bssm_likelihood_estimate_t) :: estimate
         type(bssm_particle_filter_t) :: fit

         fit = bssm_multivariate_ekpf_draws(y, initial_mean, &
            initial_covariance, parameters, observation_noise_dimension, &
            state_noise_dimension, observation, transition, &
            initial_normals(:, :, draw_index), &
            state_normals(:, :, :, draw_index), &
            resampling_uniforms(:, :, draw_index), iekf_iterations, &
            convergence_tolerance)
         estimate%value = fit%log_likelihood
         estimate%info = fit%info
      end function estimator

   end function bssm_multivariate_ekpf_pmmh_draws

   function bssm_multivariate_ekpf_pmmh(y, initial_mean, initial_covariance, &
      initial_parameters, observation_noise_dimension, state_noise_dimension, &
      observation, transition, prior, iekf_iterations, &
      convergence_tolerance, proposal_factor, particles, iterations, &
      target_acceptance, adaptation_exponent) result(out)
      !! Run multivariate EKF-proposal particle PMMH with shared randomness.
      real(dp), intent(in) :: y(:, :) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: initial_parameters(:) !! Initial parameter values.
      real(dp), intent(in) :: convergence_tolerance !! Convergence tolerance.
      real(dp), intent(in) :: proposal_factor(:, :) !! Proposal factor.
      integer, intent(in) :: observation_noise_dimension !! Observation noise dimension.
      integer, intent(in) :: state_noise_dimension !! State noise dimension.
      integer, intent(in) :: iekf_iterations !! Number of iterated extended Kalman filter iterations.
      integer, intent(in) :: particles !! Number of particles.
      integer, intent(in) :: iterations !! Number of algorithm iterations.
      procedure(bssm_multivariate_gaussian_observation_t) :: observation !! Observed value or vector.
      procedure(bssm_nonlinear_transition_t) :: transition !! State transition matrix.
      procedure(bssm_parameter_log_density_t) :: prior !! Prior-distribution specification.
      real(dp), intent(in), optional :: target_acceptance !! Target acceptance.
      real(dp), intent(in), optional :: adaptation_exponent !! Adaptation exponent.
      type(bssm_mcmc_t) :: out
      real(dp), allocatable :: parameter_normals(:, :), acceptance_uniforms(:)
      real(dp), allocatable :: initial_normals(:, :, :)
      real(dp), allocatable :: state_normals(:, :, :, :)
      real(dp), allocatable :: resampling_uniforms(:, :, :)
      real(dp) :: exponent
      integer :: draw_index, time

      if (size(y, 1) < 1 .or. size(y, 2) < 1 .or. &
         size(initial_mean) < 1 .or. size(initial_parameters) < 1 .or. &
         observation_noise_dimension < 1 .or. state_noise_dimension < 1 .or. &
         iekf_iterations < 0 .or. convergence_tolerance < 0.0_dp .or. &
         particles < 1 .or. iterations < 1 .or. &
         (present(adaptation_exponent) .and. &
         .not. present(target_acceptance))) then
         out%info = 1
         return
      end if
      allocate(parameter_normals(size(initial_parameters), iterations))
      allocate(acceptance_uniforms(iterations))
      allocate(initial_normals(size(initial_mean), particles, iterations + 1))
      allocate(state_normals(size(initial_mean), particles, size(y, 2), &
         iterations + 1))
      allocate(resampling_uniforms(particles, size(y, 2), iterations + 1))
      call random_standard_normal_matrix(parameter_normals)
      call random_number(acceptance_uniforms)
      call random_number(resampling_uniforms)
      do draw_index = 1, iterations + 1
         call random_standard_normal_matrix(initial_normals(:, :, draw_index))
         do time = 1, size(y, 2)
            call random_standard_normal_matrix( &
               state_normals(:, :, time, draw_index))
         end do
      end do
      exponent = 0.6_dp
      if (present(adaptation_exponent)) exponent = adaptation_exponent
      if (present(target_acceptance)) then
         out = bssm_multivariate_ekpf_pmmh_draws(y, initial_mean, &
            initial_covariance, initial_parameters, &
            observation_noise_dimension, state_noise_dimension, observation, &
            transition, prior, iekf_iterations, convergence_tolerance, &
            proposal_factor, parameter_normals, acceptance_uniforms, &
            initial_normals, state_normals, resampling_uniforms, &
            target_acceptance, exponent)
      else
         out = bssm_multivariate_ekpf_pmmh_draws(y, initial_mean, &
            initial_covariance, initial_parameters, &
            observation_noise_dimension, state_noise_dimension, observation, &
            transition, prior, iekf_iterations, convergence_tolerance, &
            proposal_factor, parameter_normals, acceptance_uniforms, &
            initial_normals, state_normals, resampling_uniforms)
      end if
   end function bssm_multivariate_ekpf_pmmh

   pure function bssm_multivariate_ekpf_da_pmmh_draws(y, initial_mean, &
      initial_covariance, initial_parameters, observation_noise_dimension, &
      state_noise_dimension, observation, transition_model, transition, prior, &
      iekf_iterations, convergence_tolerance, proposal_factor, &
      parameter_normals, first_stage_uniforms, second_stage_uniforms, &
      initial_normals, state_normals, resampling_uniforms, target_acceptance, &
      adaptation_exponent) result(out)
      !! Run draw-driven IEKF-screened multivariate EKPF-PMMH.
      real(dp), intent(in) :: y(:, :) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: initial_parameters(:) !! Initial parameter values.
      real(dp), intent(in) :: convergence_tolerance !! Convergence tolerance.
      real(dp), intent(in) :: proposal_factor(:, :) !! Proposal factor.
      real(dp), intent(in) :: parameter_normals(:, :) !! Parameter normals.
      real(dp), intent(in) :: first_stage_uniforms(:) !! First stage uniforms.
      real(dp), intent(in) :: second_stage_uniforms(:) !! Second stage uniforms.
      real(dp), intent(in) :: initial_normals(:, :, :) !! Initial normals.
      real(dp), intent(in) :: state_normals(:, :, :, :) !! State normals.
      real(dp), intent(in) :: resampling_uniforms(:, :, :) !! Resampling uniforms.
      integer, intent(in) :: observation_noise_dimension !! Observation noise dimension.
      integer, intent(in) :: state_noise_dimension !! State noise dimension.
      integer, intent(in) :: iekf_iterations !! Number of iterated extended Kalman filter iterations.
      procedure(bssm_multivariate_gaussian_observation_t) :: observation !! Observed value or vector.
      procedure(bssm_nonlinear_transition_jacobian_t) :: transition_model !! Transition model callback procedure.
      procedure(bssm_nonlinear_transition_t) :: transition !! State transition matrix.
      procedure(bssm_parameter_log_density_t) :: prior !! Prior-distribution specification.
      real(dp), intent(in), optional :: target_acceptance !! Target acceptance.
      real(dp), intent(in), optional :: adaptation_exponent !! Adaptation exponent.
      type(bssm_da_mcmc_t) :: out
      real(dp) :: exponent
      integer :: state, particles, times, iterations

      state = size(initial_mean)
      particles = size(initial_normals, 2)
      times = size(y, 2)
      iterations = size(parameter_normals, 2)
      if (size(y, 1) < 1 .or. state < 1 .or. particles < 1 .or. &
         times < 1 .or. iterations < 1 .or. &
         observation_noise_dimension < 1 .or. state_noise_dimension < 1 .or. &
         iekf_iterations < 0 .or. convergence_tolerance < 0.0_dp .or. &
         any(shape(initial_covariance) /= [state, state]) .or. &
         any(shape(initial_normals) /= &
         [state, particles, iterations + 1]) .or. &
         any(shape(state_normals) /= &
         [state, particles, times, iterations + 1]) .or. &
         any(shape(resampling_uniforms) /= &
         [particles, times, iterations + 1]) .or. &
         .not. all(ieee_is_finite(initial_normals)) .or. &
         .not. all(ieee_is_finite(state_normals)) .or. &
         .not. all(ieee_is_finite(resampling_uniforms)) .or. &
         (present(adaptation_exponent) .and. &
         .not. present(target_acceptance))) then
         out%info = 1
         return
      end if
      exponent = 0.6_dp
      if (present(adaptation_exponent)) exponent = adaptation_exponent
      if (present(target_acceptance)) then
         out = bssm_da_pmmh_kernel_draws(initial_parameters, prior, &
            coarse_estimator, fine_estimator, proposal_factor, &
            parameter_normals, first_stage_uniforms, second_stage_uniforms, &
            target_acceptance, exponent)
      else
         out = bssm_da_pmmh_kernel_draws(initial_parameters, prior, &
            coarse_estimator, fine_estimator, proposal_factor, &
            parameter_normals, first_stage_uniforms, second_stage_uniforms)
      end if

   contains

      pure function coarse_estimator(parameters, draw_index) result(estimate)
         !! Evaluate the deterministic multivariate IEKF likelihood.
         real(dp), intent(in) :: parameters(:) !! Model parameter values.
         integer, intent(in) :: draw_index !! Index of draw.
         type(bssm_likelihood_estimate_t) :: estimate
         type(bssm_multivariate_ekf_t) :: fit

         fit = bssm_multivariate_iekf(y, initial_mean, initial_covariance, &
            parameters, observation_noise_dimension, state_noise_dimension, &
            observation, transition_model, iekf_iterations, &
            convergence_tolerance)
         estimate%value = fit%log_likelihood + 0.0_dp*real(draw_index, dp)
         estimate%info = fit%info
      end function coarse_estimator

      pure function fine_estimator(parameters, draw_index) result(estimate)
         !! Evaluate the multivariate EKF-proposal particle likelihood.
         real(dp), intent(in) :: parameters(:) !! Model parameter values.
         integer, intent(in) :: draw_index !! Index of draw.
         type(bssm_likelihood_estimate_t) :: estimate
         type(bssm_particle_filter_t) :: fit

         fit = bssm_multivariate_ekpf_draws(y, initial_mean, &
            initial_covariance, parameters, observation_noise_dimension, &
            state_noise_dimension, observation, transition, &
            initial_normals(:, :, draw_index), &
            state_normals(:, :, :, draw_index), &
            resampling_uniforms(:, :, draw_index), iekf_iterations, &
            convergence_tolerance)
         estimate%value = fit%log_likelihood
         estimate%info = fit%info
      end function fine_estimator

   end function bssm_multivariate_ekpf_da_pmmh_draws

   function bssm_multivariate_ekpf_da_pmmh(y, initial_mean, &
      initial_covariance, initial_parameters, observation_noise_dimension, &
      state_noise_dimension, observation, transition_model, transition, prior, &
      iekf_iterations, convergence_tolerance, proposal_factor, particles, &
      iterations, target_acceptance, adaptation_exponent) result(out)
      !! Run IEKF-screened multivariate EKPF-PMMH with shared randomness.
      real(dp), intent(in) :: y(:, :) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: initial_parameters(:) !! Initial parameter values.
      real(dp), intent(in) :: convergence_tolerance !! Convergence tolerance.
      real(dp), intent(in) :: proposal_factor(:, :) !! Proposal factor.
      integer, intent(in) :: observation_noise_dimension !! Observation noise dimension.
      integer, intent(in) :: state_noise_dimension !! State noise dimension.
      integer, intent(in) :: iekf_iterations !! Number of iterated extended Kalman filter iterations.
      integer, intent(in) :: particles !! Number of particles.
      integer, intent(in) :: iterations !! Number of algorithm iterations.
      procedure(bssm_multivariate_gaussian_observation_t) :: observation !! Observed value or vector.
      procedure(bssm_nonlinear_transition_jacobian_t) :: transition_model !! Transition model callback procedure.
      procedure(bssm_nonlinear_transition_t) :: transition !! State transition matrix.
      procedure(bssm_parameter_log_density_t) :: prior !! Prior-distribution specification.
      real(dp), intent(in), optional :: target_acceptance !! Target acceptance.
      real(dp), intent(in), optional :: adaptation_exponent !! Adaptation exponent.
      type(bssm_da_mcmc_t) :: out
      real(dp), allocatable :: parameter_normals(:, :)
      real(dp), allocatable :: first_stage_uniforms(:), second_stage_uniforms(:)
      real(dp), allocatable :: initial_normals(:, :, :)
      real(dp), allocatable :: state_normals(:, :, :, :)
      real(dp), allocatable :: resampling_uniforms(:, :, :)
      real(dp) :: exponent
      integer :: draw_index, time

      if (size(y, 1) < 1 .or. size(y, 2) < 1 .or. &
         size(initial_mean) < 1 .or. size(initial_parameters) < 1 .or. &
         observation_noise_dimension < 1 .or. state_noise_dimension < 1 .or. &
         iekf_iterations < 0 .or. convergence_tolerance < 0.0_dp .or. &
         particles < 1 .or. iterations < 1 .or. &
         (present(adaptation_exponent) .and. &
         .not. present(target_acceptance))) then
         out%info = 1
         return
      end if
      allocate(parameter_normals(size(initial_parameters), iterations))
      allocate(first_stage_uniforms(iterations), second_stage_uniforms(iterations))
      allocate(initial_normals(size(initial_mean), particles, iterations + 1))
      allocate(state_normals(size(initial_mean), particles, size(y, 2), &
         iterations + 1))
      allocate(resampling_uniforms(particles, size(y, 2), iterations + 1))
      call random_standard_normal_matrix(parameter_normals)
      call random_number(first_stage_uniforms)
      call random_number(second_stage_uniforms)
      call random_number(resampling_uniforms)
      do draw_index = 1, iterations + 1
         call random_standard_normal_matrix(initial_normals(:, :, draw_index))
         do time = 1, size(y, 2)
            call random_standard_normal_matrix( &
               state_normals(:, :, time, draw_index))
         end do
      end do
      exponent = 0.6_dp
      if (present(adaptation_exponent)) exponent = adaptation_exponent
      if (present(target_acceptance)) then
         out = bssm_multivariate_ekpf_da_pmmh_draws(y, initial_mean, &
            initial_covariance, initial_parameters, &
            observation_noise_dimension, state_noise_dimension, observation, &
            transition_model, transition, prior, iekf_iterations, &
            convergence_tolerance, proposal_factor, parameter_normals, &
            first_stage_uniforms, second_stage_uniforms, initial_normals, &
            state_normals, resampling_uniforms, target_acceptance, exponent)
      else
         out = bssm_multivariate_ekpf_da_pmmh_draws(y, initial_mean, &
            initial_covariance, initial_parameters, &
            observation_noise_dimension, state_noise_dimension, observation, &
            transition_model, transition, prior, iekf_iterations, &
            convergence_tolerance, proposal_factor, parameter_normals, &
            first_stage_uniforms, second_stage_uniforms, initial_normals, &
            state_normals, resampling_uniforms)
      end if
   end function bssm_multivariate_ekpf_da_pmmh

   pure function bssm_nonlinear_ekpf_da_pmmh_draws(y, initial_mean, &
      initial_covariance, initial_parameters, noise_dimension, observation, &
      transition_model, transition, prior, iekf_iterations, &
      convergence_tolerance, proposal_factor, parameter_normals, &
      first_stage_uniforms, second_stage_uniforms, initial_normals, &
      state_normals, resampling_uniforms, target_acceptance, &
      adaptation_exponent) result(out)
      !! Run delayed acceptance with EKF and EKF-proposal particle likelihoods.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: initial_parameters(:) !! Initial parameter values.
      real(dp), intent(in) :: convergence_tolerance !! Convergence tolerance.
      real(dp), intent(in) :: proposal_factor(:, :) !! Proposal factor.
      real(dp), intent(in) :: parameter_normals(:, :) !! Parameter normals.
      real(dp), intent(in) :: first_stage_uniforms(:) !! First stage uniforms.
      real(dp), intent(in) :: second_stage_uniforms(:) !! Second stage uniforms.
      real(dp), intent(in) :: initial_normals(:, :, :) !! Initial normals.
      real(dp), intent(in) :: state_normals(:, :, :, :) !! State normals.
      real(dp), intent(in) :: resampling_uniforms(:, :, :) !! Resampling uniforms.
      integer, intent(in) :: noise_dimension !! Noise dimension.
      integer, intent(in) :: iekf_iterations !! Number of iterated extended Kalman filter iterations.
      procedure(bssm_gaussian_observation_t) :: observation !! Observed value or vector.
      procedure(bssm_nonlinear_transition_jacobian_t) :: transition_model !! Transition model callback procedure.
      procedure(bssm_nonlinear_transition_t) :: transition !! State transition matrix.
      procedure(bssm_parameter_log_density_t) :: prior !! Prior-distribution specification.
      real(dp), intent(in), optional :: target_acceptance !! Target acceptance.
      real(dp), intent(in), optional :: adaptation_exponent !! Adaptation exponent.
      type(bssm_da_mcmc_t) :: out
      real(dp) :: exponent
      integer :: state, particles, times, iterations

      state = size(initial_mean)
      particles = size(initial_normals, 2)
      times = size(y)
      iterations = size(parameter_normals, 2)
      if (state < 1 .or. particles < 1 .or. times < 1 .or. iterations < 1 .or. &
         noise_dimension < 1 .or. iekf_iterations < 0 .or. &
         convergence_tolerance < 0.0_dp .or. &
         any(shape(initial_covariance) /= [state, state]) .or. &
         any(shape(initial_normals) /= &
         [state, particles, iterations + 1]) .or. &
         any(shape(state_normals) /= &
         [state, particles, times, iterations + 1]) .or. &
         any(shape(resampling_uniforms) /= &
         [particles, times, iterations + 1]) .or. &
         .not. all(ieee_is_finite(initial_normals)) .or. &
         .not. all(ieee_is_finite(state_normals)) .or. &
         .not. all(ieee_is_finite(resampling_uniforms)) .or. &
         (present(adaptation_exponent) .and. &
         .not. present(target_acceptance))) then
         out%info = 1
         return
      end if
      exponent = 0.6_dp
      if (present(adaptation_exponent)) exponent = adaptation_exponent
      if (present(target_acceptance)) then
         out = bssm_da_pmmh_kernel_draws(initial_parameters, prior, &
            coarse_estimator, fine_estimator, proposal_factor, &
            parameter_normals, first_stage_uniforms, second_stage_uniforms, &
            target_acceptance, exponent)
      else
         out = bssm_da_pmmh_kernel_draws(initial_parameters, prior, &
            coarse_estimator, fine_estimator, proposal_factor, &
            parameter_normals, first_stage_uniforms, second_stage_uniforms)
      end if

   contains

      pure function coarse_estimator(parameters, draw_index) result(estimate)
         !! Evaluate the deterministic extended-Kalman likelihood.
         real(dp), intent(in) :: parameters(:) !! Model parameter values.
         integer, intent(in) :: draw_index !! Index of draw.
         type(bssm_likelihood_estimate_t) :: estimate
         type(bssm_ekf_t) :: fit

         fit = bssm_iekf(y, initial_mean, initial_covariance, parameters, &
            noise_dimension, observation, transition_model, iekf_iterations, &
            convergence_tolerance)
         estimate%value = fit%log_likelihood + 0.0_dp*real(draw_index, dp)
         estimate%info = fit%info
      end function coarse_estimator

      pure function fine_estimator(parameters, draw_index) result(estimate)
         !! Evaluate the EKF-proposal particle likelihood.
         real(dp), intent(in) :: parameters(:) !! Model parameter values.
         integer, intent(in) :: draw_index !! Index of draw.
         type(bssm_likelihood_estimate_t) :: estimate
         type(bssm_particle_filter_t) :: fit

         fit = bssm_ekpf_draws(y, initial_mean, initial_covariance, parameters, &
            noise_dimension, observation, transition, &
            initial_normals(:, :, draw_index), &
            state_normals(:, :, :, draw_index), &
            resampling_uniforms(:, :, draw_index), iekf_iterations, &
            convergence_tolerance)
         estimate%value = fit%log_likelihood
         estimate%info = fit%info
      end function fine_estimator

   end function bssm_nonlinear_ekpf_da_pmmh_draws

   function bssm_nonlinear_ekpf_da_pmmh(y, initial_mean, &
      initial_covariance, initial_parameters, noise_dimension, observation, &
      transition_model, transition, prior, iekf_iterations, &
      convergence_tolerance, proposal_factor, particles, iterations, &
      target_acceptance, adaptation_exponent) result(out)
      !! Run EKF-screened EKF-proposal PMMH using the shared random stream.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: initial_parameters(:) !! Initial parameter values.
      real(dp), intent(in) :: convergence_tolerance !! Convergence tolerance.
      real(dp), intent(in) :: proposal_factor(:, :) !! Proposal factor.
      integer, intent(in) :: noise_dimension !! Noise dimension.
      integer, intent(in) :: iekf_iterations !! Number of iterated extended Kalman filter iterations.
      integer, intent(in) :: particles !! Number of particles.
      integer, intent(in) :: iterations !! Number of algorithm iterations.
      procedure(bssm_gaussian_observation_t) :: observation !! Observed value or vector.
      procedure(bssm_nonlinear_transition_jacobian_t) :: transition_model !! Transition model callback procedure.
      procedure(bssm_nonlinear_transition_t) :: transition !! State transition matrix.
      procedure(bssm_parameter_log_density_t) :: prior !! Prior-distribution specification.
      real(dp), intent(in), optional :: target_acceptance !! Target acceptance.
      real(dp), intent(in), optional :: adaptation_exponent !! Adaptation exponent.
      type(bssm_da_mcmc_t) :: out
      real(dp), allocatable :: parameter_normals(:, :)
      real(dp), allocatable :: first_stage_uniforms(:), second_stage_uniforms(:)
      real(dp), allocatable :: initial_normals(:, :, :)
      real(dp), allocatable :: state_normals(:, :, :, :)
      real(dp), allocatable :: resampling_uniforms(:, :, :)
      real(dp) :: exponent
      integer :: draw_index, time

      if (particles < 1 .or. iterations < 1 .or. noise_dimension < 1 .or. &
         size(y) < 1 .or. size(initial_mean) < 1 .or. &
         size(initial_parameters) < 1 .or. &
         (present(adaptation_exponent) .and. &
         .not. present(target_acceptance))) then
         out%info = 1
         return
      end if
      allocate(parameter_normals(size(initial_parameters), iterations))
      allocate(first_stage_uniforms(iterations), second_stage_uniforms(iterations))
      allocate(initial_normals(size(initial_mean), particles, iterations + 1))
      allocate(state_normals(size(initial_mean), particles, size(y), &
         iterations + 1))
      allocate(resampling_uniforms(particles, size(y), iterations + 1))
      call random_standard_normal_matrix(parameter_normals)
      call random_number(first_stage_uniforms)
      call random_number(second_stage_uniforms)
      call random_number(resampling_uniforms)
      do draw_index = 1, iterations + 1
         call random_standard_normal_matrix(initial_normals(:, :, draw_index))
         do time = 1, size(y)
            call random_standard_normal_matrix( &
               state_normals(:, :, time, draw_index))
         end do
      end do
      exponent = 0.6_dp
      if (present(adaptation_exponent)) exponent = adaptation_exponent
      if (present(target_acceptance)) then
         out = bssm_nonlinear_ekpf_da_pmmh_draws(y, initial_mean, &
            initial_covariance, initial_parameters, noise_dimension, &
            observation, transition_model, transition, prior, &
            iekf_iterations, convergence_tolerance, proposal_factor, &
            parameter_normals, first_stage_uniforms, second_stage_uniforms, &
            initial_normals, state_normals, resampling_uniforms, &
            target_acceptance, exponent)
      else
         out = bssm_nonlinear_ekpf_da_pmmh_draws(y, initial_mean, &
            initial_covariance, initial_parameters, noise_dimension, &
            observation, transition_model, transition, prior, &
            iekf_iterations, convergence_tolerance, proposal_factor, &
            parameter_normals, first_stage_uniforms, second_stage_uniforms, &
            initial_normals, state_normals, resampling_uniforms)
      end if
   end function bssm_nonlinear_ekpf_da_pmmh

   pure subroutine bssm_trace_resampled_sde_path(filter, terminal_uniform, &
      trajectory, info)
      !! Trace one SDE path from the final resampled particle generation.
      type(bssm_particle_filter_t), intent(in) :: filter !! Filter.
      real(dp), intent(in) :: terminal_uniform !! Terminal uniform.
      real(dp), allocatable, intent(out) :: trajectory(:, :) !! Trajectory.
      integer, intent(out) :: info !! Status code; zero indicates success.
      integer :: state, particles, times, current, time

      info = 1
      if (filter%info /= 0 .or. .not. allocated(filter%particles) .or. &
         .not. allocated(filter%weights) .or. &
         .not. allocated(filter%ancestors)) return
      state = size(filter%particles, 1)
      particles = size(filter%particles, 2)
      times = size(filter%weights, 2)
      if (state < 1 .or. particles < 1 .or. times < 1 .or. &
         size(filter%particles, 3) /= times + 1 .or. &
         any(shape(filter%ancestors) /= [particles, times]) .or. &
         .not. ieee_is_finite(terminal_uniform) .or. &
         terminal_uniform < 0.0_dp .or. terminal_uniform > 1.0_dp) return
      allocate(trajectory(state, times + 1))
      current = min(particles, 1 + int(terminal_uniform*real(particles, dp)))
      trajectory(:, times + 1) = filter%particles(:, current, times + 1)
      do time = times, 1, -1
         current = filter%ancestors(current, time)
         if (current < 1 .or. current > particles) then
            info = 2
            return
         end if
         trajectory(:, time) = filter%particles(:, current, time)
      end do
      info = 0
   end subroutine bssm_trace_resampled_sde_path

   pure function bssm_particle_smoother(filter) result(out)
      !! Trace particle genealogies and compute weighted smoothed moments.
      type(bssm_particle_filter_t), intent(in) :: filter !! Filter.
      type(bssm_particle_smoother_t) :: out
      real(dp), allocatable :: sample(:, :), probability(:)
      integer :: state, particles, times, endpoint, current, time

      if (filter%info /= 0 .or. .not. allocated(filter%particles) .or. &
         .not. allocated(filter%weights) .or. &
         .not. allocated(filter%ancestors)) then
         out%info = 1
         return
      end if
      state = size(filter%particles, 1)
      particles = size(filter%particles, 2)
      times = size(filter%weights, 2)
      if (size(filter%particles, 3) /= times + 1 .or. &
         any(shape(filter%ancestors) /= [particles, times])) then
         out%info = 1
         return
      end if
      allocate(out%trajectories(state, times, particles))
      allocate(out%mean(state, times), out%covariance(state, state, times))
      allocate(sample(state, particles), probability(particles))
      probability = filter%weights(:, times)
      do endpoint = 1, particles
         current = endpoint
         out%trajectories(:, times, endpoint) = &
            filter%particles(:, current, times)
         do time = times - 1, 1, -1
            current = filter%ancestors(current, time)
            if (current < 1 .or. current > particles) then
               out%info = 2
               return
            end if
            out%trajectories(:, time, endpoint) = &
               filter%particles(:, current, time)
         end do
      end do
      do time = 1, times
         sample = out%trajectories(:, time, :)
         call weighted_particle_summary(sample, probability, &
            out%mean(:, time), out%covariance(:, :, time))
      end do
   end function bssm_particle_smoother

   pure function bssm_psi_particle_smoother(filter) result(out)
      !! Trace terminal psi-filter genealogies through all n+1 states.
      type(bssm_particle_filter_t), intent(in) :: filter !! Filter.
      type(bssm_particle_smoother_t) :: out
      real(dp), allocatable :: sample(:, :)
      integer :: state, particles, times, endpoint, current, time

      if (filter%info /= 0 .or. .not. allocated(filter%particles) .or. &
         .not. allocated(filter%weights) .or. &
         .not. allocated(filter%ancestors)) then
         out%info = 1
         return
      end if
      state = size(filter%particles, 1)
      particles = size(filter%particles, 2)
      times = size(filter%weights, 2)
      if (state < 1 .or. particles < 1 .or. times < 1 .or. &
         size(filter%particles, 3) /= times + 1 .or. &
         any(shape(filter%ancestors) /= [particles, times])) then
         out%info = 1
         return
      end if
      allocate(out%trajectories(state, times + 1, particles))
      allocate(out%mean(state, times + 1))
      allocate(out%covariance(state, state, times + 1))
      allocate(sample(state, particles))
      do endpoint = 1, particles
         current = endpoint
         out%trajectories(:, times + 1, endpoint) = &
            filter%particles(:, current, times + 1)
         do time = times, 1, -1
            current = filter%ancestors(current, time)
            if (current < 1 .or. current > particles) then
               out%info = 2
               return
            end if
            out%trajectories(:, time, endpoint) = &
               filter%particles(:, current, time)
         end do
      end do
      do time = 1, times + 1
         sample = out%trajectories(:, time, :)
         call particle_summary(sample, out%mean(:, time), &
            out%covariance(:, :, time))
      end do
   end function bssm_psi_particle_smoother

   pure function bssm_post_corrected_particle_moments(filters, sample_weight) &
      result(out)
      !! Combine terminal particle smoothers over corrected chain weights.
      type(bssm_particle_filter_t), intent(in) :: filters(:) !! Filters.
      real(dp), intent(in) :: sample_weight(:) !! Sample weight.
      type(bssm_state_posterior_t) :: out
      type(bssm_particle_smoother_t) :: smoother
      real(dp), allocatable :: conditional_mean(:, :, :)
      real(dp), allocatable :: conditional_covariance(:, :, :, :)
      integer :: state, times, samples, sample

      samples = size(filters)
      if (samples < 1 .or. size(sample_weight) /= samples) then
         out%info = 1
         return
      end if
      smoother = bssm_psi_particle_smoother(filters(1))
      if (smoother%info /= 0) then
         out%info = 100 + smoother%info
         return
      end if
      state = size(smoother%mean, 1)
      times = size(smoother%mean, 2)
      allocate(conditional_mean(state, times, samples))
      allocate(conditional_covariance(state, state, times, samples))
      conditional_mean(:, :, 1) = smoother%mean
      conditional_covariance(:, :, :, 1) = smoother%covariance
      do sample = 2, samples
         smoother = bssm_psi_particle_smoother(filters(sample))
         if (smoother%info /= 0) then
            out%info = 100 + smoother%info
            return
         end if
         if (any(shape(smoother%mean) /= [state, times]) .or. &
            any(shape(smoother%covariance) /= [state, state, times])) then
            out%info = 1
            return
         end if
         conditional_mean(:, :, sample) = smoother%mean
         conditional_covariance(:, :, :, sample) = smoother%covariance
      end do
      out = bssm_corrected_state_moments(conditional_mean, &
         conditional_covariance, sample_weight)
   end function bssm_post_corrected_particle_moments

   pure function bssm_post_corrected_particle_trajectories_draws(filters, &
      sample_weight, uniforms) result(out)
      !! Resample terminal particle paths over corrected chain weights.
      type(bssm_particle_filter_t), intent(in) :: filters(:) !! Filters.
      real(dp), intent(in) :: sample_weight(:) !! Sample weight.
      real(dp), intent(in) :: uniforms(:, :) !! Uniforms.
      type(bssm_trajectory_sample_t) :: out
      type(bssm_particle_smoother_t) :: smoother
      real(dp), allocatable :: trajectories(:, :, :, :)
      real(dp), allocatable :: particle_weight(:, :)
      integer :: state, times, particles, samples, sample

      samples = size(filters)
      if (samples < 1 .or. size(sample_weight) /= samples) then
         out%info = 1
         return
      end if
      smoother = bssm_psi_particle_smoother(filters(1))
      if (smoother%info /= 0) then
         out%info = 100 + smoother%info
         return
      end if
      state = size(smoother%trajectories, 1)
      times = size(smoother%trajectories, 2)
      particles = size(smoother%trajectories, 3)
      allocate(trajectories(state, times, particles, samples))
      allocate(particle_weight(particles, samples))
      trajectories(:, :, :, 1) = smoother%trajectories
      do sample = 2, samples
         smoother = bssm_psi_particle_smoother(filters(sample))
         if (smoother%info /= 0) then
            out%info = 100 + smoother%info
            return
         end if
         if (any(shape(smoother%trajectories) /= &
            [state, times, particles])) then
            out%info = 1
            return
         end if
         trajectories(:, :, :, sample) = smoother%trajectories
      end do
      particle_weight = 1.0_dp/real(particles, dp)
      out = bssm_corrected_trajectory_draws(trajectories, particle_weight, &
         sample_weight, uniforms)
   end function bssm_post_corrected_particle_trajectories_draws

   function bssm_post_corrected_particle_trajectories(filters, sample_weight, &
      draws) result(out)
      !! Resample corrected terminal particle paths with shared randomness.
      type(bssm_particle_filter_t), intent(in) :: filters(:) !! Filters.
      real(dp), intent(in) :: sample_weight(:) !! Sample weight.
      integer, intent(in) :: draws !! Draws.
      type(bssm_trajectory_sample_t) :: out
      real(dp), allocatable :: uniforms(:, :)

      if (draws < 1) then
         out%info = 1
         return
      end if
      allocate(uniforms(2, draws))
      call random_number(uniforms)
      out = bssm_post_corrected_particle_trajectories_draws(filters, &
         sample_weight, uniforms)
   end function bssm_post_corrected_particle_trajectories

   pure subroutine form_sigma_points(mean, factor, scale, sigma)
      !! Construct symmetric scaled unscented-transform sigma points.
      real(dp), intent(in) :: mean(:) !! Mean value or vector.
      real(dp), intent(in) :: factor(:, :) !! Factor.
      real(dp), intent(in) :: scale !! Scale.
      real(dp), intent(out) :: sigma(:, :) !! Scale parameter or standard deviation.
      integer :: state

      sigma(:, 1) = mean
      do state = 1, size(mean)
         sigma(:, state + 1) = mean + scale*factor(:, state)
         sigma(:, state + 1 + size(mean)) = mean - scale*factor(:, state)
      end do
   end subroutine form_sigma_points

   pure subroutine bssm_ekf_smoothing_inputs(fit, parameters, &
      noise_dimension, transition, model, filtered, info)
      !! Adapt nonlinear EKF output to the shared linearized RTS smoother.
      type(bssm_ekf_t), intent(in) :: fit !! Previously fitted model.
      real(dp), intent(in) :: parameters(:) !! Model parameter values.
      integer, intent(in) :: noise_dimension !! Noise dimension.
      procedure(bssm_nonlinear_transition_jacobian_t) :: transition !! State transition matrix.
      type(ssm_model_t), intent(out) :: model !! Model specification.
      type(kfs_filter_t), intent(out) :: filtered !! Filtered.
      integer, intent(out) :: info !! Status code; zero indicates success.
      real(dp), allocatable :: mean(:), noise_loading(:, :)
      integer :: state, times, time

      state = size(fit%filtered_mean, 1)
      times = size(fit%filtered_mean, 2)
      allocate(model%y(times, 1), model%a1(state))
      allocate(model%transition(state, state, times))
      allocate(mean(state), noise_loading(state, noise_dimension))
      model%y = 0.0_dp
      model%a1 = fit%predicted_mean(:, 1)
      model%transition = 0.0_dp
      do time = 1, times - 1
         call transition(time, fit%filtered_mean(:, time), parameters, mean, &
            model%transition(:, :, time), noise_loading)
         if (.not. all(ieee_is_finite(mean)) .or. &
            .not. all(ieee_is_finite(model%transition(:, :, time))) .or. &
            .not. all(ieee_is_finite(noise_loading))) then
            info = time
            return
         end if
      end do
      allocate(filtered%a_pred(state, times))
      allocate(filtered%a_filt(state, times))
      allocate(filtered%p_pred(state, state, times))
      allocate(filtered%p_filt(state, state, times))
      filtered%a_pred = fit%predicted_mean(:, 1:times)
      filtered%a_filt = fit%filtered_mean
      filtered%p_pred = fit%predicted_covariance(:, :, 1:times)
      filtered%p_filt = fit%filtered_covariance
      info = 0
   end subroutine bssm_ekf_smoothing_inputs

   pure subroutine bssm_multivariate_ekf_smoothing_inputs(fit, parameters, &
      noise_dimension, transition, model, filtered, info)
      !! Adapt multivariate EKF output to the shared linearized RTS smoother.
      type(bssm_multivariate_ekf_t), intent(in) :: fit !! Previously fitted model.
      real(dp), intent(in) :: parameters(:) !! Model parameter values.
      integer, intent(in) :: noise_dimension !! Noise dimension.
      procedure(bssm_nonlinear_transition_jacobian_t) :: transition !! State transition matrix.
      type(ssm_model_t), intent(out) :: model !! Model specification.
      type(kfs_filter_t), intent(out) :: filtered !! Filtered.
      integer, intent(out) :: info !! Status code; zero indicates success.
      real(dp), allocatable :: mean(:), noise_loading(:, :)
      integer :: state, times, time

      state = size(fit%filtered_mean, 1)
      times = size(fit%filtered_mean, 2)
      allocate(model%y(times, 1), model%a1(state))
      allocate(model%transition(state, state, times))
      allocate(mean(state), noise_loading(state, noise_dimension))
      model%y = 0.0_dp
      model%a1 = fit%predicted_mean(:, 1)
      model%transition = 0.0_dp
      do time = 1, times - 1
         call transition(time, fit%filtered_mean(:, time), parameters, mean, &
            model%transition(:, :, time), noise_loading)
         if (.not. all(ieee_is_finite(mean)) .or. &
            .not. all(ieee_is_finite(model%transition(:, :, time))) .or. &
            .not. all(ieee_is_finite(noise_loading))) then
            info = time
            return
         end if
      end do
      allocate(filtered%a_pred(state, times))
      allocate(filtered%a_filt(state, times))
      allocate(filtered%p_pred(state, state, times))
      allocate(filtered%p_filt(state, state, times))
      filtered%a_pred = fit%predicted_mean(:, 1:times)
      filtered%a_filt = fit%filtered_mean
      filtered%p_pred = fit%predicted_covariance(:, :, 1:times)
      filtered%p_filt = fit%filtered_covariance
      info = 0
   end subroutine bssm_multivariate_ekf_smoothing_inputs

   pure subroutine multivariate_ekf_proposal_update(time, observation_value, &
      prior_mean, prior_covariance, parameters, noise_dimension, observation, &
      posterior_mean, posterior_covariance, final_innovation, &
      final_innovation_covariance, log_density, info, max_iterations, &
      convergence_tolerance, iterations_used)
      !! Form one partially observed multivariate extended Kalman update.
      integer, intent(in) :: time !! Observation times.
      integer, intent(in) :: noise_dimension !! Noise dimension.
      real(dp), intent(in) :: observation_value(:) !! Observation value.
      real(dp), intent(in) :: prior_mean(:) !! Prior mean.
      real(dp), intent(in) :: prior_covariance(:, :) !! Prior covariance.
      real(dp), intent(in) :: parameters(:) !! Model parameter values.
      procedure(bssm_multivariate_gaussian_observation_t) :: observation !! Observed value or vector.
      real(dp), intent(out) :: posterior_mean(:) !! Posterior mean.
      real(dp), intent(out) :: posterior_covariance(:, :) !! Posterior covariance matrix.
      real(dp), intent(out) :: final_innovation(:) !! Final innovation.
      real(dp), intent(out) :: final_innovation_covariance(:, :) !! Final innovation covariance matrix.
      real(dp), intent(out) :: log_density !! Log-density value.
      integer, intent(out) :: info !! Status code; zero indicates success.
      integer, intent(out) :: iterations_used !! Iterations used.
      integer, intent(in) :: max_iterations !! Maximum number of algorithm iterations.
      real(dp), intent(in) :: convergence_tolerance !! Convergence tolerance.
      real(dp), allocatable :: observed_mean(:), observed_jacobian(:, :)
      real(dp), allocatable :: observation_covariance(:, :)
      real(dp), allocatable :: innovation_covariance(:, :), inverse(:, :)
      real(dp), allocatable :: gain(:, :), innovation(:), estimate(:)
      real(dp), allocatable :: next_estimate(:), residual_operator(:, :)
      integer, allocatable :: observed(:)
      real(dp) :: log_determinant, difference
      integer :: series, state, observed_count, component, iteration

      series = size(observation_value)
      state = size(prior_mean)
      observed_count = count(ieee_is_finite(observation_value))
      allocate(observed(observed_count))
      observed = pack([(component, component = 1, series)], &
         ieee_is_finite(observation_value))
      allocate(observed_mean(observed_count))
      allocate(observed_jacobian(observed_count, state))
      allocate(observation_covariance(observed_count, observed_count))
      allocate(innovation_covariance(observed_count, observed_count))
      allocate(inverse(observed_count, observed_count))
      allocate(gain(state, observed_count), innovation(observed_count))
      allocate(estimate(state))
      allocate(next_estimate(state), residual_operator(state, state))
      call evaluate_observation(prior_mean, observed_mean, &
         observed_jacobian, observation_covariance, info)
      if (info /= 0) return
      innovation_covariance = matmul(matmul(observed_jacobian, &
         prior_covariance), transpose(observed_jacobian)) + &
         observation_covariance
      innovation_covariance = 0.5_dp*(innovation_covariance + &
         transpose(innovation_covariance))
      call inverse_logdet(innovation_covariance, inverse, log_determinant, &
         info, 1.0e-12_dp)
      if (info /= 0) then
         info = 2
         return
      end if
      gain = matmul(matmul(prior_covariance, transpose(observed_jacobian)), &
         inverse)
      innovation = observation_value(observed) - observed_mean
      estimate = prior_mean + matmul(gain, innovation)
      iteration = 0
      do while (iteration < max_iterations)
         call evaluate_observation(estimate, observed_mean, &
            observed_jacobian, observation_covariance, info)
         if (info /= 0) return
         innovation_covariance = matmul(matmul(observed_jacobian, &
            prior_covariance), transpose(observed_jacobian)) + &
            observation_covariance
         innovation_covariance = 0.5_dp*(innovation_covariance + &
            transpose(innovation_covariance))
         call inverse_logdet(innovation_covariance, inverse, log_determinant, &
            info, 1.0e-12_dp)
         if (info /= 0) then
            info = 2
            return
         end if
         innovation = observation_value(observed) - observed_mean - &
            matmul(observed_jacobian, prior_mean - estimate)
         gain = matmul(matmul(prior_covariance, &
            transpose(observed_jacobian)), inverse)
         next_estimate = prior_mean + matmul(gain, innovation)
         difference = sum((estimate - next_estimate)**2)/real(state, dp)
         estimate = next_estimate
         iteration = iteration + 1
         if (difference <= convergence_tolerance) exit
      end do
      posterior_mean = estimate
      residual_operator = identity_matrix(state) - &
         matmul(gain, observed_jacobian)
      posterior_covariance = &
         matmul(matmul(residual_operator, prior_covariance), &
         transpose(residual_operator)) + &
         matmul(matmul(gain, observation_covariance), transpose(gain))
      posterior_covariance = 0.5_dp*(posterior_covariance + &
         transpose(posterior_covariance))
      final_innovation = 0.0_dp
      final_innovation(observed) = innovation
      final_innovation_covariance = 0.0_dp
      do component = 1, observed_count
         final_innovation_covariance(observed, observed(component)) = &
            innovation_covariance(:, component)
      end do
      log_density = -0.5_dp*(real(observed_count, dp)* &
         log(2.0_dp*acos(-1.0_dp)) + &
         log_determinant + dot_product(innovation, matmul(inverse, innovation)))
      iterations_used = iteration
      info = 0

   contains

      pure subroutine evaluate_observation(state_value, compact_mean, &
         compact_jacobian, compact_covariance, status)
         !! Evaluate and compact the nonlinear observation moments.
         real(dp), intent(in) :: state_value(:) !! State value.
         real(dp), intent(out) :: compact_mean(:) !! Compact mean.
         real(dp), intent(out) :: compact_jacobian(:, :) !! Compact jacobian.
         real(dp), intent(out) :: compact_covariance(:, :) !! Compact covariance matrix.
         integer, intent(out) :: status !! Status.
         real(dp), allocatable :: full_mean(:), full_jacobian(:, :)
         real(dp), allocatable :: full_noise_loading(:, :)
         real(dp), allocatable :: full_covariance(:, :)
         integer :: row, column

         allocate(full_mean(series), full_jacobian(series, state))
         allocate(full_noise_loading(series, noise_dimension))
         call observation(time, state_value, parameters, full_mean, &
            full_jacobian, full_noise_loading)
         if (.not. all(ieee_is_finite(full_mean)) .or. &
            .not. all(ieee_is_finite(full_jacobian)) .or. &
            .not. all(ieee_is_finite(full_noise_loading))) then
            status = 1
            return
         end if
         full_covariance = matmul(full_noise_loading, &
            transpose(full_noise_loading))
         compact_mean = full_mean(observed)
         do row = 1, observed_count
            compact_jacobian(row, :) = full_jacobian(observed(row), :)
            do column = 1, observed_count
               compact_covariance(row, column) = &
                  full_covariance(observed(row), observed(column))
            end do
         end do
         status = 0
      end subroutine evaluate_observation

   end subroutine multivariate_ekf_proposal_update

   pure subroutine ekf_proposal_update(time, observation_value, prior_mean, &
      prior_covariance, parameters, observation, posterior_mean, &
      posterior_covariance, info, max_iterations, convergence_tolerance, &
      iterations_used, final_innovation, final_innovation_variance)
      !! Form one scalar-observation extended Kalman proposal update.
      integer, intent(in) :: time !! Observation times.
      real(dp), intent(in) :: observation_value !! Observation value.
      real(dp), intent(in) :: prior_mean(:) !! Prior mean.
      real(dp), intent(in) :: prior_covariance(:, :) !! Prior covariance.
      real(dp), intent(in) :: parameters(:) !! Model parameter values.
      procedure(bssm_gaussian_observation_t) :: observation !! Observed value or vector.
      real(dp), intent(out) :: posterior_mean(:) !! Posterior mean.
      real(dp), intent(out) :: posterior_covariance(:, :) !! Posterior covariance matrix.
      integer, intent(out) :: info !! Status code; zero indicates success.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      real(dp), intent(in), optional :: convergence_tolerance !! Convergence tolerance.
      integer, intent(out), optional :: iterations_used !! Iterations used.
      real(dp), intent(out), optional :: final_innovation !! Final innovation.
      real(dp), intent(out), optional :: final_innovation_variance !! Final innovation variance.
      real(dp), allocatable :: jacobian(:), gain(:), residual_operator(:, :)
      real(dp), allocatable :: estimate(:), next_estimate(:)
      real(dp) :: observation_mean, standard_deviation, innovation_variance
      real(dp) :: innovation, difference, tolerance
      integer :: state, iteration, limit

      state = size(prior_mean)
      allocate(jacobian(state), gain(state), residual_operator(state, state))
      allocate(estimate(state), next_estimate(state))
      limit = 0
      if (present(max_iterations)) limit = max(0, max_iterations)
      tolerance = 1.0e-4_dp
      if (present(convergence_tolerance)) then
         tolerance = max(convergence_tolerance, 0.0_dp)
      end if
      call observation(time, prior_mean, parameters, observation_mean, &
         jacobian, standard_deviation)
      if (.not. ieee_is_finite(observation_mean) .or. &
         .not. ieee_is_finite(standard_deviation) .or. &
         standard_deviation <= 0.0_dp .or. &
         .not. all(ieee_is_finite(jacobian))) then
         info = 1
         return
      end if
      innovation_variance = dot_product(jacobian, &
         matmul(prior_covariance, jacobian)) + standard_deviation**2
      if (.not. ieee_is_finite(innovation_variance) .or. &
         innovation_variance <= 0.0_dp) then
         info = 2
         return
      end if
      gain = matmul(prior_covariance, jacobian)/innovation_variance
      innovation = observation_value - observation_mean
      estimate = prior_mean + gain*innovation
      iteration = 0
      do while (iteration < limit)
         call observation(time, estimate, parameters, observation_mean, &
            jacobian, standard_deviation)
         if (.not. ieee_is_finite(observation_mean) .or. &
            .not. ieee_is_finite(standard_deviation) .or. &
            standard_deviation <= 0.0_dp .or. &
            .not. all(ieee_is_finite(jacobian))) then
            info = 1
            return
         end if
         innovation_variance = dot_product(jacobian, &
            matmul(prior_covariance, jacobian)) + standard_deviation**2
         if (.not. ieee_is_finite(innovation_variance) .or. &
            innovation_variance <= 0.0_dp) then
            info = 2
            return
         end if
         innovation = observation_value - observation_mean - &
            dot_product(jacobian, prior_mean - estimate)
         gain = matmul(prior_covariance, jacobian)/innovation_variance
         next_estimate = prior_mean + gain*innovation
         difference = sum((estimate - next_estimate)**2)/real(state, dp)
         estimate = next_estimate
         iteration = iteration + 1
         if (difference <= tolerance) exit
      end do
      posterior_mean = estimate
      residual_operator = identity_matrix(state) - outer_product(gain, jacobian)
      posterior_covariance = &
         matmul(matmul(residual_operator, prior_covariance), &
         transpose(residual_operator)) + standard_deviation**2* &
         outer_product(gain, gain)
      posterior_covariance = 0.5_dp*(posterior_covariance + &
         transpose(posterior_covariance))
      if (present(iterations_used)) iterations_used = iteration
      if (present(final_innovation)) final_innovation = innovation
      if (present(final_innovation_variance)) then
         final_innovation_variance = innovation_variance
      end if
      info = 0
   end subroutine ekf_proposal_update

   pure function multivariate_gaussian_observation_log_density(time, value, &
      state, parameters, noise_dimension, observation) result(log_density)
      !! Evaluate a partially observed nonlinear Gaussian vector density.
      integer, intent(in) :: time !! Observation times.
      integer, intent(in) :: noise_dimension !! Noise dimension.
      real(dp), intent(in) :: value(:) !! Input value.
      real(dp), intent(in) :: state(:) !! State vector or state sequence.
      real(dp), intent(in) :: parameters(:) !! Model parameter values.
      procedure(bssm_multivariate_gaussian_observation_t) :: observation !! Observed value or vector.
      real(dp) :: log_density
      real(dp), allocatable :: mean(:), jacobian(:, :), noise_loading(:, :)
      real(dp), allocatable :: covariance(:, :), compact_covariance(:, :)
      integer, allocatable :: observed(:)
      integer :: series, observed_count, component, row, column

      series = size(value)
      observed_count = count(ieee_is_finite(value))
      if (series < 1 .or. observed_count < 1 .or. noise_dimension < 1) then
         log_density = 0.0_dp
         return
      end if
      allocate(mean(series), jacobian(series, size(state)))
      allocate(noise_loading(series, noise_dimension))
      call observation(time, state, parameters, mean, jacobian, noise_loading)
      if (.not. all(ieee_is_finite(mean)) .or. &
         .not. all(ieee_is_finite(jacobian)) .or. &
         .not. all(ieee_is_finite(noise_loading))) then
         log_density = -huge(1.0_dp)
         return
      end if
      allocate(observed(observed_count))
      observed = pack([(component, component = 1, series)], &
         ieee_is_finite(value))
      covariance = matmul(noise_loading, transpose(noise_loading))
      allocate(compact_covariance(observed_count, observed_count))
      do row = 1, observed_count
         do column = 1, observed_count
            compact_covariance(row, column) = &
               covariance(observed(row), observed(column))
         end do
      end do
      log_density = multivariate_normal_log_density(value(observed), &
         mean(observed), compact_covariance)
   end function multivariate_gaussian_observation_log_density

   pure subroutine ram_proposal_update(factor, normal, probability, target, &
      iteration, exponent, info)
      !! Apply one robust adaptive Metropolis proposal-factor update.
      real(dp), allocatable, intent(inout) :: factor(:, :) !! Factor, updated in place.
      real(dp), intent(in) :: normal(:) !! Normal.
      real(dp), intent(in) :: probability !! Probability value.
      real(dp), intent(in) :: target !! Target.
      real(dp), intent(in) :: exponent !! Exponent.
      integer, intent(in) :: iteration !! Iteration.
      integer, intent(out) :: info !! Status code; zero indicates success.
      real(dp), allocatable :: direction(:), covariance(:, :), next_factor(:, :)
      real(dp) :: norm_squared, step_size

      info = 0
      norm_squared = dot_product(normal, normal)
      if (norm_squared <= tiny(1.0_dp)) return
      direction = matmul(factor, normal)
      step_size = real(iteration, dp)**(-exponent)
      covariance = matmul(factor, transpose(factor)) + &
         step_size*(probability - target)* &
         outer_product(direction, direction)/norm_squared
      call positive_semidefinite_factor(covariance, next_factor, info)
      if (info == 0) call move_alloc(next_factor, factor)
   end subroutine ram_proposal_update

   pure subroutine positive_semidefinite_factor(covariance, factor, info)
      !! Form a square root of a symmetric positive-semidefinite matrix.
      real(dp), intent(in) :: covariance(:, :) !! Covariance matrix.
      real(dp), allocatable, intent(out) :: factor(:, :) !! Factor.
      integer, intent(out) :: info !! Status code; zero indicates success.
      real(dp), allocatable :: values(:), vectors(:, :)
      real(dp) :: tolerance
      integer :: n

      n = size(covariance, 1)
      if (size(covariance, 2) /= n) then
         allocate(factor(0, 0))
         info = 1
         return
      end if
      call symmetric_eigen(0.5_dp*(covariance + transpose(covariance)), &
         values, vectors, info)
      if (info /= 0) return
      tolerance = 1.0e-10_dp*max(1.0_dp, maxval(abs(values)))
      if (minval(values) < -tolerance) then
         info = 2
         return
      end if
      allocate(factor(n, n))
      factor = vectors*spread(sqrt(max(values, 0.0_dp)), 1, n)
   end subroutine positive_semidefinite_factor

   pure subroutine bssm_form_nonlinear_prediction(state, parameters, &
      observation, response_normals, signal, mean, response, info)
      !! Evaluate nonlinear Gaussian means and add observation disturbances.
      real(dp), intent(in) :: state(:, :, :) !! State vector or state sequence.
      real(dp), intent(in) :: parameters(:, :) !! Model parameter values.
      procedure(bssm_nonlinear_prediction_observation_t) :: observation !! Observed value or vector.
      real(dp), intent(in) :: response_normals(:, :, :) !! Response normals.
      real(dp), allocatable, intent(out) :: signal(:, :, :) !! Signal.
      real(dp), allocatable, intent(out) :: mean(:, :, :) !! Mean value or vector.
      real(dp), allocatable, intent(out) :: response(:, :, :) !! Response observations.
      integer, intent(out) :: info !! Status code; zero indicates success.
      real(dp), allocatable :: observation_mean(:), noise_loading(:, :)
      integer :: observation_dimension, times, samples, sample, time

      observation_dimension = size(response_normals, 1)
      times = size(state, 2)
      samples = size(state, 3)
      allocate(signal(observation_dimension, times, samples))
      allocate(mean(observation_dimension, times, samples))
      allocate(response(observation_dimension, times, samples))
      allocate(observation_mean(observation_dimension))
      allocate(noise_loading(observation_dimension, observation_dimension))
      info = 0
      do sample = 1, samples
         do time = 1, times
            call observation(time, state(:, time, sample), &
               parameters(:, sample), observation_mean, noise_loading)
            if (.not. all(ieee_is_finite(observation_mean)) .or. &
               .not. all(ieee_is_finite(noise_loading))) then
               info = 10 + time
               return
            end if
            signal(:, time, sample) = observation_mean
            mean(:, time, sample) = observation_mean
            response(:, time, sample) = observation_mean + &
               matmul(noise_loading, response_normals(:, time, sample))
         end do
      end do
      if (.not. all(ieee_is_finite(response))) info = 1
   end subroutine bssm_form_nonlinear_prediction

   pure subroutine bssm_prediction_options(series, state, times, offset, &
      auxiliary, state_offset, offset_work, auxiliary_work, &
      state_offset_work, info)
      !! Validate and expand optional predictive offsets and exposures.
      integer, intent(in) :: series !! Time-series observations.
      integer, intent(in) :: state !! State vector or state sequence.
      integer, intent(in) :: times !! Times.
      real(dp), intent(in), optional :: offset(:, :) !! Known additive offset.
      real(dp), intent(in), optional :: auxiliary(:, :) !! Auxiliary.
      real(dp), intent(in), optional :: state_offset(:, :) !! State offset.
      real(dp), allocatable, intent(out) :: offset_work(:, :) !! Offset work.
      real(dp), allocatable, intent(out) :: auxiliary_work(:, :) !! Auxiliary work.
      real(dp), allocatable, intent(out) :: state_offset_work(:, :) !! State offset work.
      integer, intent(out) :: info !! Status code; zero indicates success.

      allocate(offset_work(series, times), auxiliary_work(series, times))
      allocate(state_offset_work(state, times))
      offset_work = 0.0_dp
      auxiliary_work = 1.0_dp
      state_offset_work = 0.0_dp
      info = 1
      if (present(offset)) then
         if (any(shape(offset) /= [series, times]) .or. &
            .not. all(ieee_is_finite(offset))) return
         offset_work = offset
      end if
      if (present(auxiliary)) then
         if (any(shape(auxiliary) /= [series, times]) .or. &
            .not. all(ieee_is_finite(auxiliary)) .or. &
            any(auxiliary < 0.0_dp)) return
         auxiliary_work = auxiliary
      end if
      if (present(state_offset)) then
         if (any(shape(state_offset) /= [state, times]) .or. &
            .not. all(ieee_is_finite(state_offset))) return
         state_offset_work = state_offset
      end if
      info = 0
   end subroutine bssm_prediction_options

   pure subroutine bssm_form_predictive_response(state, observation_loading, &
      distribution, phi, response_normals, response_uniforms, offset, &
      auxiliary, signal, mean, response, info, observation_noise_loading)
      !! Transform predictive states into means and sampled observations.
      real(dp), intent(in) :: state(:, :, :) !! State vector or state sequence.
      real(dp), intent(in) :: observation_loading(:, :, :) !! Observation loading matrix.
      integer, intent(in) :: distribution(:) !! Probability-distribution specification.
      real(dp), intent(in) :: phi(:) !! Autoregressive or model coefficient.
      real(dp), intent(in) :: response_normals(:, :, :) !! Response normals.
      real(dp), intent(in) :: response_uniforms(:, :, :) !! Response uniforms.
      real(dp), intent(in) :: offset(:, :) !! Known additive offset.
      real(dp), intent(in) :: auxiliary(:, :) !! Auxiliary.
      real(dp), intent(in), optional :: observation_noise_loading(:, :, :) !! Observation noise loading.
      real(dp), allocatable, intent(out) :: signal(:, :, :) !! Signal.
      real(dp), allocatable, intent(out) :: mean(:, :, :) !! Mean value or vector.
      real(dp), allocatable, intent(out) :: response(:, :, :) !! Response observations.
      integer, intent(out) :: info !! Status code; zero indicates success.
      real(dp), allocatable :: gaussian_noise(:)
      real(dp) :: total_mean, gamma_value
      integer :: count
      integer :: series, times, samples, sample, time, component, tz, status

      series = size(observation_loading, 1)
      times = size(state, 2)
      samples = size(state, 3)
      allocate(signal(series, times, samples), mean(series, times, samples))
      allocate(response(series, times, samples))
      allocate(gaussian_noise(series))
      info = 0
      if (.not. all(ieee_is_finite(phi))) then
         info = 1
         return
      end if
      if (present(observation_noise_loading)) then
         if (size(observation_noise_loading, 1) /= series .or. &
            size(observation_noise_loading, 2) /= series .or. &
            .not. bssm_valid_prediction_extent( &
            size(observation_noise_loading, 3), times) .or. &
            .not. all(ieee_is_finite(observation_noise_loading))) then
            info = 1
            return
         end if
      end if
      do sample = 1, samples
         do time = 1, times
            tz = min(time, size(observation_loading, 3))
            gaussian_noise = 0.0_dp
            if (present(observation_noise_loading)) then
               gaussian_noise = matmul(observation_noise_loading(:, :, &
                  min(time, size(observation_noise_loading, 3))), &
                  response_normals(:, time, sample))
            end if
            signal(:, time, sample) = offset(:, time) + &
               matmul(observation_loading(:, :, tz), state(:, time, sample))
            do component = 1, series
               select case (distribution(component))
               case (bssm_svm)
                  if (phi(component) <= 0.0_dp) then
                     info = 10 + component
                     return
                  end if
                  mean(component, time, sample) = 0.0_dp
                  response(component, time, sample) = phi(component)* &
                     sqrt(bounded_exp(signal(component, time, sample)))* &
                     response_normals(component, time, sample)
               case (bssm_poisson)
                  mean(component, time, sample) = &
                     bounded_exp(signal(component, time, sample))
                  total_mean = auxiliary(component, time)* &
                     mean(component, time, sample)
                  call bssm_poisson_quantile(total_mean, &
                     response_uniforms(component, time, sample), count, &
                     status)
                  response(component, time, sample) = real(count, dp)
                  if (status /= 0) then
                     info = 20 + component
                     return
                  end if
               case (bssm_binomial)
                  if (abs(auxiliary(component, time) - &
                     anint(auxiliary(component, time))) > 1.0e-8_dp) then
                     info = 30 + component
                     return
                  end if
                  mean(component, time, sample) = &
                     logistic(signal(component, time, sample))
                  call bssm_binomial_quantile( &
                     nint(auxiliary(component, time)), &
                     mean(component, time, sample), &
                     response_uniforms(component, time, sample), count, status)
                  response(component, time, sample) = real(count, dp)
                  if (status /= 0) then
                     info = 40 + component
                     return
                  end if
               case (bssm_negative_binomial)
                  if (phi(component) <= 0.0_dp) then
                     info = 50 + component
                     return
                  end if
                  mean(component, time, sample) = &
                     bounded_exp(signal(component, time, sample))
                  total_mean = auxiliary(component, time)* &
                     mean(component, time, sample)
                  call bssm_negative_binomial_quantile(total_mean, &
                     phi(component), response_uniforms(component, time, &
                     sample), count, status)
                  response(component, time, sample) = real(count, dp)
                  if (status /= 0) then
                     info = 60 + component
                     return
                  end if
               case (bssm_gamma)
                  if (phi(component) <= 0.0_dp) then
                     info = 70 + component
                     return
                  end if
                  mean(component, time, sample) = &
                     bounded_exp(signal(component, time, sample))
                  total_mean = auxiliary(component, time)* &
                     mean(component, time, sample)
                  call bssm_gamma_quantile(phi(component), &
                     response_uniforms(component, time, sample), gamma_value, &
                     status)
                  response(component, time, sample) = &
                     total_mean*gamma_value/phi(component)
                  if (status /= 0) then
                     info = 80 + component
                     return
                  end if
               case (bssm_gaussian)
                  if (.not. present(observation_noise_loading) .and. &
                     phi(component) <= 0.0_dp) then
                     info = 90 + component
                     return
                  end if
                  mean(component, time, sample) = &
                     signal(component, time, sample)
                  if (present(observation_noise_loading)) then
                     response(component, time, sample) = &
                        mean(component, time, sample) + &
                        gaussian_noise(component)
                  else
                     response(component, time, sample) = &
                        mean(component, time, sample) + phi(component)* &
                        response_normals(component, time, sample)
                  end if
               end select
            end do
         end do
      end do
      if (.not. all(ieee_is_finite(signal)) .or. &
         .not. all(ieee_is_finite(mean)) .or. &
         .not. all(ieee_is_finite(response))) info = 2
   end subroutine bssm_form_predictive_response

   pure subroutine bssm_poisson_quantile(mean, uniform, value, info)
      !! Invert a Poisson CDF using normalized log probabilities.
      real(dp), intent(in) :: mean !! Mean value or vector.
      real(dp), intent(in) :: uniform !! Uniform.
      integer, intent(out) :: value !! Input value.
      integer, intent(out) :: info !! Status code; zero indicates success.
      real(dp) :: log_probability, maximum, total, cumulative, target
      integer :: upper, candidate

      value = 0
      info = 0
      if (mean < 0.0_dp .or. .not. ieee_is_finite(mean)) then
         info = 1
         return
      end if
      if (mean == 0.0_dp) return
      upper = ceiling(mean + 12.0_dp*sqrt(mean + 1.0_dp) + 20.0_dp)
      if (upper > 1000000) then
         info = 2
         return
      end if
      maximum = -huge(1.0_dp)
      do candidate = 0, upper
         log_probability = -mean + real(candidate, dp)*log(mean) - &
            log_gamma(real(candidate + 1, dp))
         maximum = max(maximum, log_probability)
      end do
      total = 0.0_dp
      do candidate = 0, upper
         log_probability = -mean + real(candidate, dp)*log(mean) - &
            log_gamma(real(candidate + 1, dp))
         total = total + exp(log_probability - maximum)
      end do
      target = uniform*total
      cumulative = 0.0_dp
      do candidate = 0, upper
         log_probability = -mean + real(candidate, dp)*log(mean) - &
            log_gamma(real(candidate + 1, dp))
         cumulative = cumulative + exp(log_probability - maximum)
         if (cumulative >= target) then
            value = candidate
            return
         end if
      end do
      value = upper
   end subroutine bssm_poisson_quantile

   pure subroutine bssm_binomial_quantile(trials, probability, uniform, &
      value, info)
      !! Invert a binomial CDF using normalized log probabilities.
      integer, intent(in) :: trials !! Trials.
      real(dp), intent(in) :: probability !! Probability value.
      real(dp), intent(in) :: uniform !! Uniform.
      integer, intent(out) :: value !! Input value.
      integer, intent(out) :: info !! Status code; zero indicates success.
      real(dp) :: log_probability, maximum, total, cumulative, target
      integer :: candidate

      value = 0
      info = 0
      if (trials < 0 .or. trials > 1000000 .or. probability < 0.0_dp .or. &
         probability > 1.0_dp) then
         info = 1
         return
      end if
      if (trials == 0 .or. probability == 0.0_dp) return
      if (probability == 1.0_dp) then
         value = trials
         return
      end if
      maximum = -huge(1.0_dp)
      do candidate = 0, trials
         log_probability = log_gamma(real(trials + 1, dp)) - &
            log_gamma(real(candidate + 1, dp)) - &
            log_gamma(real(trials - candidate + 1, dp)) + &
            real(candidate, dp)*log(probability) + &
            real(trials - candidate, dp)*log(1.0_dp - probability)
         maximum = max(maximum, log_probability)
      end do
      total = 0.0_dp
      do candidate = 0, trials
         log_probability = log_gamma(real(trials + 1, dp)) - &
            log_gamma(real(candidate + 1, dp)) - &
            log_gamma(real(trials - candidate + 1, dp)) + &
            real(candidate, dp)*log(probability) + &
            real(trials - candidate, dp)*log(1.0_dp - probability)
         total = total + exp(log_probability - maximum)
      end do
      target = uniform*total
      cumulative = 0.0_dp
      do candidate = 0, trials
         log_probability = log_gamma(real(trials + 1, dp)) - &
            log_gamma(real(candidate + 1, dp)) - &
            log_gamma(real(trials - candidate + 1, dp)) + &
            real(candidate, dp)*log(probability) + &
            real(trials - candidate, dp)*log(1.0_dp - probability)
         cumulative = cumulative + exp(log_probability - maximum)
         if (cumulative >= target) then
            value = candidate
            return
         end if
      end do
      value = trials
   end subroutine bssm_binomial_quantile

   pure subroutine bssm_negative_binomial_quantile(mean, dispersion, uniform, &
      value, info)
      !! Invert a mean-dispersion negative-binomial CDF.
      real(dp), intent(in) :: mean !! Mean value or vector.
      real(dp), intent(in) :: dispersion !! Dispersion.
      real(dp), intent(in) :: uniform !! Uniform.
      integer, intent(out) :: value !! Input value.
      integer, intent(out) :: info !! Status code; zero indicates success.
      real(dp) :: variance, log_probability, maximum, total, cumulative
      real(dp) :: target, log_denominator
      integer :: upper, candidate

      value = 0
      info = 0
      if (mean < 0.0_dp .or. dispersion <= 0.0_dp .or. &
         .not. ieee_is_finite(mean)) then
         info = 1
         return
      end if
      if (mean == 0.0_dp) return
      variance = mean + mean*mean/dispersion
      upper = ceiling(mean + 20.0_dp*sqrt(variance) + 100.0_dp)
      if (upper > 1000000) then
         info = 2
         return
      end if
      log_denominator = log(dispersion + mean)
      maximum = -huge(1.0_dp)
      do candidate = 0, upper
         log_probability = log_gamma(real(candidate, dp) + dispersion) - &
            log_gamma(dispersion) - log_gamma(real(candidate + 1, dp)) + &
            dispersion*(log(dispersion) - log_denominator) + &
            real(candidate, dp)*(log(mean) - log_denominator)
         maximum = max(maximum, log_probability)
      end do
      total = 0.0_dp
      do candidate = 0, upper
         log_probability = log_gamma(real(candidate, dp) + dispersion) - &
            log_gamma(dispersion) - log_gamma(real(candidate + 1, dp)) + &
            dispersion*(log(dispersion) - log_denominator) + &
            real(candidate, dp)*(log(mean) - log_denominator)
         total = total + exp(log_probability - maximum)
      end do
      target = uniform*total
      cumulative = 0.0_dp
      do candidate = 0, upper
         log_probability = log_gamma(real(candidate, dp) + dispersion) - &
            log_gamma(dispersion) - log_gamma(real(candidate + 1, dp)) + &
            dispersion*(log(dispersion) - log_denominator) + &
            real(candidate, dp)*(log(mean) - log_denominator)
         cumulative = cumulative + exp(log_probability - maximum)
         if (cumulative >= target) then
            value = candidate
            return
         end if
      end do
      value = upper
   end subroutine bssm_negative_binomial_quantile

   pure subroutine bssm_gamma_quantile(shape, uniform, value, info)
      !! Invert a unit-scale Gamma CDF by bracketed bisection.
      real(dp), intent(in) :: shape !! Shape.
      real(dp), intent(in) :: uniform !! Uniform.
      real(dp), intent(out) :: value !! Input value.
      integer, intent(out) :: info !! Status code; zero indicates success.
      real(dp) :: lower, upper, target
      integer :: iteration

      value = 0.0_dp
      info = 0
      if (shape <= 0.0_dp .or. .not. ieee_is_finite(shape)) then
         info = 1
         return
      end if
      if (uniform <= 0.0_dp) return
      target = min(uniform, 1.0_dp - epsilon(1.0_dp))
      lower = 0.0_dp
      upper = max(1.0_dp, shape)
      do iteration = 1, 100
         if (bssm_regularized_gamma(shape, upper) >= target) exit
         upper = 2.0_dp*upper
         if (.not. ieee_is_finite(upper)) then
            info = 2
            return
         end if
      end do
      do iteration = 1, 100
         value = 0.5_dp*(lower + upper)
         if (bssm_regularized_gamma(shape, value) < target) then
            lower = value
         else
            upper = value
         end if
      end do
      value = 0.5_dp*(lower + upper)
   end subroutine bssm_gamma_quantile

   pure real(dp) function bssm_regularized_gamma(shape, value) &
      result(probability)
      !! Evaluate the regularized lower incomplete Gamma function.
      real(dp), intent(in) :: shape !! Shape.
      real(dp), intent(in) :: value !! Input value.
      real(dp) :: sum, term, ap, b, c, d, h, factor, an, delta
      integer :: iteration

      if (value <= 0.0_dp) then
         probability = 0.0_dp
         return
      end if
      factor = exp(-value + shape*log(value) - log_gamma(shape))
      if (value < shape + 1.0_dp) then
         ap = shape
         sum = 1.0_dp/shape
         term = sum
         do iteration = 1, 200
            ap = ap + 1.0_dp
            term = term*value/ap
            sum = sum + term
            if (abs(term) <= epsilon(1.0_dp)*abs(sum)) exit
         end do
         probability = factor*sum
      else
         b = value + 1.0_dp - shape
         c = 1.0_dp/tiny(1.0_dp)
         d = 1.0_dp/max(abs(b), tiny(1.0_dp))
         if (b < 0.0_dp) d = -d
         h = d
         do iteration = 1, 200
            an = -real(iteration, dp)*(real(iteration, dp) - shape)
            b = b + 2.0_dp
            d = an*d + b
            if (abs(d) < tiny(1.0_dp)) d = tiny(1.0_dp)
            c = b + an/c
            if (abs(c) < tiny(1.0_dp)) c = tiny(1.0_dp)
            d = 1.0_dp/d
            delta = d*c
            h = h*delta
            if (abs(delta - 1.0_dp) <= 10.0_dp*epsilon(1.0_dp)) exit
         end do
         probability = 1.0_dp - factor*h
      end if
      probability = min(1.0_dp, max(0.0_dp, probability))
   end function bssm_regularized_gamma

   pure logical function bssm_valid_prediction_extent(extent, times) &
      result(valid)
      !! Test a constant or fully time-varying predictive matrix extent.
      integer, intent(in) :: extent !! Extent.
      integer, intent(in) :: times !! Times.

      valid = extent == 1 .or. extent >= times
   end function bssm_valid_prediction_extent

   pure function bssm_conditional_simulation_draws(center, factor, &
      conditional_matrix, normals) result(out)
      !! Draw paths from a forward conditional Gaussian representation.
      real(dp), intent(in) :: center(:, :) !! Center.
      real(dp), intent(in) :: factor(:, :, :) !! Factor.
      real(dp), intent(in) :: conditional_matrix(:, :, :) !! Conditional matrix.
      real(dp), intent(in) :: normals(:, :, :) !! Independent standard-normal draws.
      type(bssm_simulation_smoother_t) :: out
      real(dp), allocatable :: difference(:)
      integer :: state, times, samples, sample, time

      state = size(center, 1)
      times = size(center, 2)
      samples = size(normals, 2)
      if (state < 1 .or. times < 1 .or. samples < 1 .or. &
         any(shape(factor) /= [state, state, times]) .or. &
         any(shape(conditional_matrix) /= [state, state, times]) .or. &
         any(shape(normals) /= [state, samples, times]) .or. &
         .not. all(ieee_is_finite(center)) .or. &
         .not. all(ieee_is_finite(factor)) .or. &
         .not. all(ieee_is_finite(conditional_matrix)) .or. &
         .not. all(ieee_is_finite(normals))) then
         out%info = 1
         return
      end if
      allocate(out%trajectories(state, times, samples))
      allocate(out%mean(state, times), out%covariance(state, state, times))
      allocate(difference(state))
      out%mean = 0.0_dp
      out%covariance = 0.0_dp
      do sample = 1, samples
         out%trajectories(:, 1, sample) = center(:, 1) + &
            matmul(factor(:, :, 1), normals(:, sample, 1))
         do time = 2, times
            out%trajectories(:, time, sample) = center(:, time) + &
               matmul(conditional_matrix(:, :, time), &
               out%trajectories(:, time - 1, sample) - center(:, time - 1)) + &
               matmul(factor(:, :, time), normals(:, sample, time))
         end do
         out%mean = out%mean + out%trajectories(:, :, sample)
      end do
      out%mean = out%mean/real(samples, dp)
      do sample = 1, samples
         do time = 1, times
            difference = out%trajectories(:, time, sample) - &
               out%mean(:, time)
            out%covariance(:, :, time) = out%covariance(:, :, time) + &
               outer_product(difference, difference)
         end do
      end do
      out%covariance = out%covariance/real(samples, dp)
   end function bssm_conditional_simulation_draws

   subroutine bssm_random_approximation_normals(state, samples, times, &
      use_antithetic, normals)
      !! Generate independent or paired normals for approximation paths.
      integer, intent(in) :: state !! State vector or state sequence.
      integer, intent(in) :: samples !! Samples.
      integer, intent(in) :: times !! Times.
      logical, intent(in) :: use_antithetic !! Whether to use the antithetic.
      real(dp), allocatable, intent(out) :: normals(:, :, :) !! Independent standard-normal draws.
      integer :: base_samples, paired_samples, time

      allocate(normals(state, samples, times))
      if (use_antithetic .and. samples > 1) then
         base_samples = samples/2
         paired_samples = 2*base_samples
         do time = 1, times
            call random_standard_normal_matrix( &
               normals(:, 1:base_samples, time))
            normals(:, base_samples + 1:paired_samples, time) = &
               -normals(:, 1:base_samples, time)
            if (paired_samples < samples) then
               call random_standard_normal_matrix( &
                  normals(:, samples:samples, time))
            end if
         end do
      else
         do time = 1, times
            call random_standard_normal_matrix(normals(:, :, time))
         end do
      end if
   end subroutine bssm_random_approximation_normals

   pure integer function bssm_time_index(matrix, time) result(index)
      !! Select a constant or time-varying state-space matrix slice.
      real(dp), intent(in) :: matrix(:, :, :) !! Input matrix.
      integer, intent(in) :: time !! Observation times.

      index = min(time, size(matrix, 3))
   end function bssm_time_index

   pure logical function bssm_is_observed(model, time, component) result(observed)
      !! Test whether one state-space observation is available.
      type(ssm_model_t), intent(in) :: model !! Model specification.
      integer, intent(in) :: time !! Observation times.
      integer, intent(in) :: component !! Component.

      observed = ieee_is_finite(model%y(time, component))
      if (allocated(model%missing)) then
         observed = observed .and. .not. model%missing(time, component)
      end if
   end function bssm_is_observed

   pure elemental logical function valid_log_value(value) result(valid)
      !! Test whether a log density is finite and above the rejection sentinel.
      real(dp), intent(in) :: value !! Input value.

      valid = ieee_is_finite(value) .and. &
         value > -0.5_dp*huge(1.0_dp)
   end function valid_log_value

   pure integer function discrete_index(probability, uniform) result(index)
      !! Select an index from normalized nonnegative probabilities.
      real(dp), intent(in) :: probability(:) !! Probability value.
      real(dp), intent(in) :: uniform !! Uniform.
      real(dp) :: cumulative

      cumulative = 0.0_dp
      do index = 1, size(probability) - 1
         cumulative = cumulative + probability(index)
         if (uniform < cumulative) return
      end do
      index = size(probability)
   end function discrete_index

   pure subroutine normalize_log_weights(log_weight, probability, &
      log_mean_weight, info)
      !! Normalize log weights and return their log arithmetic mean.
      real(dp), intent(in) :: log_weight(:) !! Log weight.
      real(dp), intent(out) :: probability(:) !! Probability value.
      real(dp), intent(out) :: log_mean_weight !! Log mean weight.
      integer, intent(out) :: info !! Status code; zero indicates success.
      real(dp) :: maximum, total

      maximum = maxval(log_weight)
      if (maximum <= -0.5_dp*huge(1.0_dp)) then
         probability = 0.0_dp
         log_mean_weight = -huge(1.0_dp)
         info = 1
         return
      end if
      probability = exp(log_weight - maximum)
      total = sum(probability)
      if (.not. ieee_is_finite(total) .or. total <= 0.0_dp) then
         probability = 0.0_dp
         log_mean_weight = -huge(1.0_dp)
         info = 1
         return
      end if
      probability = probability/total
      log_mean_weight = maximum + log(total/real(size(log_weight), dp))
      info = 0
   end subroutine normalize_log_weights

   pure subroutine particle_summary(particles, mean, covariance)
      !! Compute equally weighted particle means and population covariances.
      real(dp), intent(in) :: particles(:, :) !! Number of particles.
      real(dp), intent(out) :: mean(:) !! Mean value or vector.
      real(dp), intent(out) :: covariance(:, :) !! Covariance matrix.
      real(dp), allocatable :: difference(:)
      integer :: particle

      mean = sum(particles, dim=2)/real(size(particles, 2), dp)
      covariance = 0.0_dp
      allocate(difference(size(mean)))
      do particle = 1, size(particles, 2)
         difference = particles(:, particle) - mean
         covariance = covariance + outer_product(difference, difference)
      end do
      covariance = covariance/real(size(particles, 2), dp)
   end subroutine particle_summary

   pure subroutine weighted_particle_summary(particles, probability, mean, &
      covariance)
      !! Compute normalized weighted particle means and covariances.
      real(dp), intent(in) :: particles(:, :) !! Number of particles.
      real(dp), intent(in) :: probability(:) !! Probability value.
      real(dp), intent(out) :: mean(:) !! Mean value or vector.
      real(dp), intent(out) :: covariance(:, :) !! Covariance matrix.
      real(dp), allocatable :: difference(:)
      real(dp) :: total
      integer :: particle

      total = sum(probability)
      mean = matmul(particles, probability)/total
      covariance = 0.0_dp
      allocate(difference(size(mean)))
      do particle = 1, size(particles, 2)
         difference = particles(:, particle) - mean
         covariance = covariance + probability(particle)* &
            outer_product(difference, difference)
      end do
      covariance = covariance/total
   end subroutine weighted_particle_summary

   pure elemental real(dp) function softplus(x) result(value)
      !! Evaluate log(1 + exp(x)) without avoidable overflow.
      real(dp), intent(in) :: x !! Input data or predictor values.

      if (x > 0.0_dp) then
         value = x + log(1.0_dp + exp(-x))
      else
         value = log(1.0_dp + exp(x))
      end if
   end function softplus

   pure elemental real(dp) function bounded_exp(x) result(value)
      !! Exponentiate while retaining a finite positive result.
      real(dp), intent(in) :: x !! Input data or predictor values.

      value = exp(min(max(x, log(tiny(1.0_dp))), &
         log(huge(1.0_dp)) - 2.0_dp))
   end function bounded_exp

   pure elemental real(dp) function logistic(x) result(value)
      !! Evaluate the logistic inverse link without overflow.
      real(dp), intent(in) :: x !! Input data or predictor values.

      if (x >= 0.0_dp) then
         value = 1.0_dp/(1.0_dp + exp(-x))
      else
         value = exp(x)/(1.0_dp + exp(x))
      end if
      value = min(max(value, sqrt(epsilon(1.0_dp))), &
         1.0_dp - sqrt(epsilon(1.0_dp)))
   end function logistic

   pure elemental real(dp) function log_add_exp(x, y) result(value)
      !! Evaluate log(exp(x) + exp(y)) stably.
      real(dp), intent(in) :: x !! Input data or predictor values.
      real(dp), intent(in) :: y !! Response or time-series observations.
      real(dp) :: maximum

      maximum = max(x, y)
      value = maximum + log(exp(x - maximum) + exp(y - maximum))
   end function log_add_exp

end module bssm_mod
