! SPDX-License-Identifier: GPL-2.0-or-later
! SPDX-FileComment: Distinct algorithms translated from the R bvartools package.
module bvartools_mod
   !! Prior construction and posterior simulation kernels for Bayesian VAR models.
   use kind_mod, only: dp
   use linalg_mod, only: invert_matrix, kronecker_product, cholesky_lower, &
      symmetric_eigen, identity_matrix, diagonal_matrix
   use random_mod, only: multivariate_normal_from_standard, random_gamma, &
      random_standard_normal, random_uniform
   use stats_mod, only: quantile, sorted
   use mts_mod, only: mts_var_psi
   use time_series_var_utils_mod, only: build_var_data, build_varx_data
   use special_functions_mod, only: multivariate_normal_log_density
   use kfas_mod, only: ssm_model_t
   use bssm_mod, only: bssm_simulation_smoother_t, bssm_simulation_smoother
   implicit none
   private

   type, public :: bvartools_minnesota_prior_t
      !! Minnesota prior means, variances, and residual precision.
      real(dp), allocatable :: mean(:, :)
      real(dp), allocatable :: variance(:, :)
      real(dp), allocatable :: residual_precision(:, :)
      integer :: info = 0
   end type bvartools_minnesota_prior_t

   type, public :: bvartools_loglik_normal_t
      !! Observation-wise multivariate Gaussian log likelihoods.
      real(dp), allocatable :: value(:)
      integer :: info = 0
   end type bvartools_loglik_normal_t

   type, public :: bvartools_kalman_dk_t
      !! Durbin-Koopman conditional state draw and status.
      real(dp), allocatable :: state(:, :)
      integer :: info = 0
   end type bvartools_kalman_dk_t

   type, public :: bvartools_inclusion_prior_t
      !! Coefficient inclusion probabilities and selectable coefficient positions.
      real(dp), allocatable :: probability(:, :)
      integer, allocatable :: include(:)
      integer :: info = 0
   end type bvartools_inclusion_prior_t

   type, public :: bvartools_ssvs_result_t
      !! SSVS posterior probabilities, indicators, and diagonal prior precision.
      real(dp), allocatable :: posterior_probability(:)
      real(dp), allocatable :: precision(:)
      logical, allocatable :: included(:)
      integer :: info = 0
   end type bvartools_ssvs_result_t

   type, public :: bvartools_ssvs_prior_t
      !! Fixed or semiautomatic spike-and-slab scales for SSVS coefficients.
      real(dp), allocatable :: tau0(:)
      real(dp), allocatable :: tau1(:)
      real(dp), allocatable :: coefficients(:, :)
      real(dp), allocatable :: residual_covariance(:, :)
      real(dp), allocatable :: regression_standard_error(:)
      integer :: regression_parameters = 0
      integer :: covariance_parameters = 0
      logical :: semiautomatic = .false.
      integer :: info = 0
   end type bvartools_ssvs_prior_t

   type, public :: bvartools_normal_posterior_t
      !! Mean and covariance of a Gaussian coefficient posterior.
      real(dp), allocatable :: mean(:)
      real(dp), allocatable :: covariance(:, :)
      integer :: info = 0
   end type bvartools_normal_posterior_t

   type, public :: bvartools_gamma_posterior_t
      !! Independent gamma precision posterior parameters.
      real(dp), allocatable :: shape(:)
      real(dp), allocatable :: rate(:)
      integer :: info = 0
   end type bvartools_gamma_posterior_t

   type, public :: bvartools_covariance_data_t
      !! Prepared lower-triangular covariance regression data and precision.
      real(dp), allocatable :: response(:)
      real(dp), allocatable :: design(:, :)
      real(dp), allocatable :: precision(:, :)
      integer :: info = 0
   end type bvartools_covariance_data_t

   type, public :: bvartools_bvs_result_t
      !! Bayesian variable-selection indicators and posterior log odds.
      logical, allocatable :: included(:)
      real(dp), allocatable :: log_odds(:)
      integer :: info = 0
   end type bvartools_bvs_result_t

   type, public :: bvartools_bvar_draws_t
      !! Retained constant-parameter BVAR Gibbs draws.
      real(dp), allocatable :: coefficients(:, :)
      real(dp), allocatable :: structural(:, :)
      real(dp), allocatable :: structural_impact(:, :, :)
      real(dp), allocatable :: covariance(:, :, :)
      logical, allocatable :: included(:, :)
      logical, allocatable :: structural_included(:, :)
      integer :: retained_draws = 0
      integer :: info = 0
   end type bvartools_bvar_draws_t

   type, public :: bvartools_bvar_prior_bundle_t
      !! Dimensioned priors and OLS initial values for a constant BVAR.
      real(dp), allocatable :: coefficient_mean(:)
      real(dp), allocatable :: coefficient_precision(:, :)
      real(dp), allocatable :: initial_coefficients(:, :)
      real(dp), allocatable :: initial_covariance(:, :)
      real(dp), allocatable :: covariance_scale(:, :)
      real(dp), allocatable :: gamma_shape(:)
      real(dp), allocatable :: gamma_rate(:)
      real(dp), allocatable :: tau0(:)
      real(dp), allocatable :: tau1(:)
      real(dp), allocatable :: inclusion_probability(:)
      integer, allocatable :: selectable(:)
      real(dp) :: covariance_df = 0.0_dp
      character(len=8) :: covariance_prior = 'wishart'
      logical :: use_ssvs = .false.
      integer :: info = 0
   end type bvartools_bvar_prior_bundle_t

   type, public :: bvartools_bvar_fit_t
      !! Posterior draws and preparation metadata from an end-to-end BVAR fit.
      type(bvartools_bvar_draws_t) :: draws
      integer :: observations = 0
      integer :: variables = 0
      integer :: lag_order = 0
      integer :: info = 0
   end type bvartools_bvar_fit_t

   type, public :: bvartools_cointegration_draw_t
      !! One KLS posterior draw for a Bayesian cointegration model.
      real(dp), allocatable :: alpha(:, :)
      real(dp), allocatable :: beta(:, :)
      real(dp), allocatable :: pi(:, :)
      real(dp), allocatable :: gamma(:, :)
      integer :: info = 0
   end type bvartools_cointegration_draw_t

   type, public :: bvartools_bvec_draws_t
      !! Retained constant-parameter Bayesian VECM Gibbs draws.
      real(dp), allocatable :: alpha(:, :)
      real(dp), allocatable :: beta(:, :)
      real(dp), allocatable :: pi(:, :)
      real(dp), allocatable :: gamma(:, :)
      real(dp), allocatable :: structural(:, :)
      real(dp), allocatable :: structural_impact(:, :, :)
      real(dp), allocatable :: covariance(:, :, :)
      logical, allocatable :: included(:, :)
      logical, allocatable :: structural_included(:, :)
      integer :: retained_draws = 0
      integer :: info = 0
   end type bvartools_bvec_draws_t

   type, public :: bvartools_bvec_prior_bundle_t
      !! Dimensioned priors and OLS initial values for a constant Bayesian VECM.
      real(dp), allocatable :: initial_beta(:, :)
      real(dp), allocatable :: initial_gamma(:, :)
      real(dp), allocatable :: initial_covariance(:, :)
      real(dp), allocatable :: cointegration_precision(:, :)
      real(dp), allocatable :: loading_precision(:, :)
      real(dp), allocatable :: gamma_mean(:)
      real(dp), allocatable :: gamma_precision(:, :)
      real(dp), allocatable :: covariance_scale(:, :)
      real(dp), allocatable :: error_shape(:)
      real(dp), allocatable :: error_rate(:)
      real(dp), allocatable :: tau0(:)
      real(dp), allocatable :: tau1(:)
      real(dp), allocatable :: inclusion_probability(:)
      integer, allocatable :: selectable(:)
      real(dp) :: shrinkage = 0.0_dp
      real(dp) :: covariance_df = 0.0_dp
      character(len=8) :: covariance_prior = 'wishart'
      logical :: use_ssvs = .false.
      integer :: rank = 0
      integer :: info = 0
   end type bvartools_bvec_prior_bundle_t

   type, public :: bvartools_bvec_fit_t
      !! Posterior draws and preparation metadata from an end-to-end BVEC fit.
      type(bvartools_bvec_draws_t) :: draws
      type(bvartools_bvar_draws_t) :: rank_zero_draws
      integer :: observations = 0
      integer :: variables = 0
      integer :: level_lags = 0
      integer :: rank = 0
      logical :: rank_zero = .false.
      integer :: info = 0
   end type bvartools_bvec_fit_t

   type, public :: bvartools_tvp_bvar_draws_t
      !! Retained time-varying BVAR coefficient and variance draws.
      real(dp), allocatable :: states(:, :)
      real(dp), allocatable :: state_variance(:, :)
      real(dp), allocatable :: initial_state(:, :)
      real(dp), allocatable :: covariance(:, :, :)
      real(dp), allocatable :: time_covariance(:, :, :, :)
      real(dp), allocatable :: log_variance(:, :)
      real(dp), allocatable :: log_variance_state_variance(:, :)
      real(dp), allocatable :: log_variance_initial(:, :)
      logical, allocatable :: included(:, :)
      integer :: retained_draws = 0
      integer :: info = 0
   end type bvartools_tvp_bvar_draws_t

   type, public :: bvartools_tvp_bvar_prior_bundle_t
      !! Dimensioned state, observation, BVS, and volatility priors for a TVP-BVAR.
      real(dp), allocatable :: initial_state(:)
      real(dp), allocatable :: initial_covariance(:, :)
      real(dp), allocatable :: initial_state_prior_mean(:)
      real(dp), allocatable :: initial_state_prior_precision(:, :)
      real(dp), allocatable :: state_shape(:)
      real(dp), allocatable :: state_rate(:)
      real(dp), allocatable :: covariance_scale(:, :)
      real(dp), allocatable :: measurement_shape(:)
      real(dp), allocatable :: measurement_rate(:)
      real(dp), allocatable :: inclusion_probability(:)
      integer, allocatable :: selectable(:)
      real(dp), allocatable :: initial_log_variance(:, :)
      real(dp), allocatable :: initial_log_variance_level(:)
      real(dp), allocatable :: initial_log_variance_state_variance(:)
      real(dp), allocatable :: log_variance_state_shape(:)
      real(dp), allocatable :: log_variance_state_rate(:)
      real(dp), allocatable :: log_variance_initial_prior_mean(:)
      real(dp), allocatable :: log_variance_initial_prior_precision(:, :)
      real(dp), allocatable :: log_variance_offset(:)
      real(dp) :: covariance_df = 0.0_dp
      character(len=8) :: observation_prior = 'wishart'
      character(len=4) :: volatility_method = 'ksc '
      logical :: use_bvs = .false.
      integer :: info = 0
   end type bvartools_tvp_bvar_prior_bundle_t

   type, public :: bvartools_tvp_bvar_fit_t
      !! Posterior draws and preparation metadata from an end-to-end TVP-BVAR fit.
      type(bvartools_tvp_bvar_draws_t) :: draws
      integer :: observations = 0
      integer :: variables = 0
      integer :: lag_order = 0
      integer :: info = 0
   end type bvartools_tvp_bvar_fit_t

   type, public :: bvartools_tvp_covariance_draws_t
      !! Retained lower-triangular covariance-state Gibbs draws.
      real(dp), allocatable :: states(:, :)
      real(dp), allocatable :: state_variance(:, :)
      real(dp), allocatable :: initial_state(:, :)
      real(dp), allocatable :: diagonal_variance(:, :)
      real(dp), allocatable :: covariance(:, :, :, :)
      logical, allocatable :: included(:, :)
      integer :: retained_draws = 0
      integer :: info = 0
   end type bvartools_tvp_covariance_draws_t

   type, public :: bvartools_joint_tvp_bvar_draws_t
      !! Retained draws from the joint TVP-BVAR and triangular-covariance sampler.
      real(dp), allocatable :: coefficient_states(:, :)
      real(dp), allocatable :: coefficient_state_variance(:, :)
      real(dp), allocatable :: coefficient_initial_state(:, :)
      real(dp), allocatable :: covariance_states(:, :)
      real(dp), allocatable :: covariance_state_variance(:, :)
      real(dp), allocatable :: covariance_initial_state(:, :)
      real(dp), allocatable :: diagonal_variance(:, :)
      real(dp), allocatable :: covariance(:, :, :, :)
      logical, allocatable :: coefficient_included(:, :)
      logical, allocatable :: covariance_included(:, :)
      logical :: time_varying_covariance = .true.
      integer :: retained_draws = 0
      integer :: info = 0
   end type bvartools_joint_tvp_bvar_draws_t

   type, public :: bvartools_structural_tvp_bvar_prior_t
      !! Priors for a structural TVP-BVAR with triangular stochastic covariance.
      real(dp), allocatable :: initial_state(:)
      real(dp), allocatable :: initial_state_prior_mean(:)
      real(dp), allocatable :: initial_state_prior_precision(:, :)
      real(dp), allocatable :: state_shape(:)
      real(dp), allocatable :: state_rate(:)
      real(dp), allocatable :: initial_diagonal_variance(:)
      real(dp), allocatable :: initial_covariance_state(:)
      real(dp), allocatable :: covariance_initial_prior_mean(:)
      real(dp), allocatable :: covariance_initial_prior_precision(:, :)
      real(dp), allocatable :: covariance_state_shape(:)
      real(dp), allocatable :: covariance_state_rate(:)
      real(dp), allocatable :: measurement_shape(:)
      real(dp), allocatable :: measurement_rate(:)
      real(dp), allocatable :: inclusion_probability(:)
      integer, allocatable :: selectable(:)
      real(dp), allocatable :: covariance_inclusion_probability(:)
      integer, allocatable :: covariance_selectable(:)
      real(dp), allocatable :: covariance_tau0(:)
      real(dp), allocatable :: covariance_tau1(:)
      real(dp), allocatable :: covariance_ssvs_probability(:)
      integer :: reduced_form_states = 0
      integer :: structural_states = 0
      logical :: time_varying_covariance = .true.
      logical :: use_covariance_ssvs = .false.
      integer :: info = 0
   end type bvartools_structural_tvp_bvar_prior_t

   type, public :: bvartools_structural_tvp_bvar_fit_t
      !! Joint draws split into reduced-form and contemporaneous structural blocks.
      type(bvartools_joint_tvp_bvar_draws_t) :: draws
      real(dp), allocatable :: reduced_form_states(:, :)
      real(dp), allocatable :: structural_states(:, :)
      real(dp), allocatable :: structural_impact(:, :, :, :)
      logical, allocatable :: reduced_form_included(:, :)
      logical, allocatable :: structural_included(:, :)
      integer :: observations = 0
      integer :: variables = 0
      integer :: lag_order = 0
      integer :: info = 0
   end type bvartools_structural_tvp_bvar_fit_t

   type, public :: bvartools_stochastic_volatility_t
      !! One KSC random-walk log-volatility draw.
      real(dp), allocatable :: log_variance(:, :)
      integer, allocatable :: component(:, :)
      integer :: info = 0
   end type bvartools_stochastic_volatility_t

   type, public :: bvartools_tvp_bvec_draws_t
      !! Retained time-varying Bayesian VECM coefficient and variance draws.
      real(dp), allocatable :: alpha(:, :)
      real(dp), allocatable :: beta(:, :)
      real(dp), allocatable :: pi(:, :)
      real(dp), allocatable :: gamma(:, :)
      real(dp), allocatable :: structural(:, :)
      real(dp), allocatable :: structural_impact(:, :, :, :)
      real(dp), allocatable :: state_variance(:, :)
      real(dp), allocatable :: covariance(:, :, :)
      real(dp), allocatable :: time_covariance(:, :, :, :)
      real(dp), allocatable :: covariance_state(:, :)
      real(dp), allocatable :: covariance_state_variance(:, :)
      real(dp), allocatable :: covariance_initial_state(:, :)
      real(dp), allocatable :: log_variance(:, :)
      real(dp), allocatable :: log_variance_state_variance(:, :)
      real(dp), allocatable :: log_variance_initial(:, :)
      logical, allocatable :: included(:, :)
      logical, allocatable :: structural_included(:, :)
      logical, allocatable :: covariance_included(:, :)
      integer :: retained_draws = 0
      integer :: info = 0
   end type bvartools_tvp_bvec_draws_t

   type, public :: bvartools_tvp_bvec_prior_bundle_t
      !! Dimensioned state, cointegration, covariance, and volatility priors for a TVP-BVEC.
      real(dp), allocatable :: initial_alpha(:, :)
      real(dp), allocatable :: initial_beta(:, :)
      real(dp), allocatable :: initial_gamma(:, :)
      real(dp), allocatable :: initial_covariance(:, :)
      real(dp), allocatable :: initial_state_prior_mean(:)
      real(dp), allocatable :: initial_state_prior_precision(:, :)
      real(dp), allocatable :: state_shape(:)
      real(dp), allocatable :: state_rate(:)
      real(dp), allocatable :: beta_state_precision(:, :)
      real(dp), allocatable :: beta_initial_prior_mean(:)
      real(dp), allocatable :: beta_initial_prior_precision(:, :)
      real(dp), allocatable :: covariance_scale(:, :)
      real(dp), allocatable :: measurement_shape(:)
      real(dp), allocatable :: measurement_rate(:)
      real(dp), allocatable :: inclusion_probability(:)
      integer, allocatable :: selectable(:)
      real(dp), allocatable :: initial_covariance_state(:)
      real(dp), allocatable :: covariance_state_initial_prior_mean(:)
      real(dp), allocatable :: covariance_state_initial_prior_precision(:, :)
      real(dp), allocatable :: covariance_state_shape(:)
      real(dp), allocatable :: covariance_state_rate(:)
      real(dp), allocatable :: covariance_inclusion_probability(:)
      integer, allocatable :: covariance_selectable(:)
      real(dp), allocatable :: initial_log_variance(:, :)
      real(dp), allocatable :: initial_log_variance_level(:)
      real(dp), allocatable :: initial_log_variance_state_variance(:)
      real(dp), allocatable :: log_variance_state_shape(:)
      real(dp), allocatable :: log_variance_state_rate(:)
      real(dp), allocatable :: log_variance_initial_prior_mean(:)
      real(dp), allocatable :: log_variance_initial_prior_precision(:, :)
      real(dp), allocatable :: log_variance_offset(:)
      real(dp) :: beta_persistence = 0.999_dp
      real(dp) :: covariance_df = 0.0_dp
      character(len=8) :: observation_prior = 'wishart'
      character(len=4) :: volatility_method = 'ksc '
      logical :: use_bvs = .false.
      logical :: use_covariance_state = .false.
      logical :: use_covariance_bvs = .false.
      integer :: rank = 0
      integer :: info = 0
   end type bvartools_tvp_bvec_prior_bundle_t

   type, public :: bvartools_tvp_bvec_fit_t
      !! Posterior draws and preparation metadata from an end-to-end TVP-BVEC fit.
      type(bvartools_tvp_bvec_draws_t) :: draws
      integer :: observations = 0
      integer :: variables = 0
      integer :: level_lags = 0
      integer :: rank = 0
      integer :: info = 0
   end type bvartools_tvp_bvec_fit_t

   type, public :: bvartools_dfm_draws_t
      !! Retained Bayesian dynamic-factor-model Gibbs draws.
      real(dp), allocatable :: loadings(:, :)
      real(dp), allocatable :: factors(:, :)
      real(dp), allocatable :: measurement_variance(:, :)
      real(dp), allocatable :: transition(:, :)
      real(dp), allocatable :: factor_variance(:, :)
      integer :: retained_draws = 0
      integer :: info = 0
   end type bvartools_dfm_draws_t

   type, public :: bvartools_dfm_prior_t
      !! Dimensioned priors for one Bayesian dynamic factor model.
      real(dp), allocatable :: loading_precision(:, :)
      real(dp), allocatable :: measurement_shape(:)
      real(dp), allocatable :: measurement_rate(:)
      real(dp), allocatable :: transition_mean(:)
      real(dp), allocatable :: transition_precision(:, :)
      real(dp), allocatable :: factor_shape(:)
      real(dp), allocatable :: factor_rate(:)
      integer :: info = 0
   end type bvartools_dfm_prior_t

   type, public :: bvartools_dfm_grid_draws_t
      !! Priors and posterior draws for a grid of Bayesian dynamic factor models.
      type(bvartools_dfm_prior_t), allocatable :: prior(:)
      type(bvartools_dfm_draws_t), allocatable :: draws(:)
      integer, allocatable :: factor_count(:)
      integer, allocatable :: lag_order(:)
      integer :: failed_model = 0
      integer :: info = 0
   end type bvartools_dfm_grid_draws_t

   type, public :: bvartools_model_likelihood_data_t
      !! Posterior residual and covariance draws supplied for one model comparison.
      real(dp), allocatable :: residual(:, :, :)
      real(dp), allocatable :: covariance(:, :, :, :)
      integer :: parameter_count = 0
   end type bvartools_model_likelihood_data_t

   type, public :: bvartools_model_comparison_t
      !! Posterior-mean log likelihood and information criteria for one model.
      real(dp), allocatable :: observation_log_likelihood(:)
      real(dp) :: log_likelihood = 0.0_dp
      real(dp) :: aic = 0.0_dp
      real(dp) :: bic = 0.0_dp
      real(dp) :: hq = 0.0_dp
      integer :: parameter_count = 0
      integer :: observations = 0
      integer :: draws = 0
      integer :: failed_observation = 0
      integer :: failed_draw = 0
      integer :: info = 0
   end type bvartools_model_comparison_t

   type, public :: bvartools_model_comparison_set_t
      !! Information criteria and minimizing model positions for a model collection.
      type(bvartools_model_comparison_t), allocatable :: model(:)
      integer :: best_aic = 0
      integer :: best_bic = 0
      integer :: best_hq = 0
      integer :: failed_model = 0
      integer :: info = 0
   end type bvartools_model_comparison_set_t

   type, public :: bvartools_predictive_t
      !! Posterior predictive paths and equal-tail summaries for a Bayesian VAR.
      real(dp), allocatable :: paths(:, :, :)
      real(dp), allocatable :: mean(:, :)
      real(dp), allocatable :: lower(:, :)
      real(dp), allocatable :: median(:, :)
      real(dp), allocatable :: upper(:, :)
      real(dp) :: probability = 0.0_dp
      integer :: info = 0
   end type bvartools_predictive_t

   type, public :: bvartools_irf_t
      !! Draw-wise Bayesian VAR impulse responses and equal-tail summaries.
      real(dp), allocatable :: paths(:, :, :, :)
      real(dp), allocatable :: mean(:, :, :)
      real(dp), allocatable :: lower(:, :, :)
      real(dp), allocatable :: median(:, :, :)
      real(dp), allocatable :: upper(:, :, :)
      real(dp) :: probability = 0.0_dp
      character(len=4) :: response_type = 'feir'
      integer :: info = 0
   end type bvartools_irf_t

   type, public :: bvartools_fevd_t
      !! Draw-wise Bayesian VAR forecast-error variance decompositions.
      real(dp), allocatable :: paths(:, :, :, :)
      real(dp), allocatable :: mean(:, :, :)
      real(dp), allocatable :: lower(:, :, :)
      real(dp), allocatable :: median(:, :, :)
      real(dp), allocatable :: upper(:, :, :)
      real(dp) :: probability = 0.0_dp
      character(len=4) :: response_type = 'oir '
      logical :: normalized = .false.
      integer :: info = 0
   end type bvartools_fevd_t

   type, public :: bvartools_level_var_draws_t
      !! Constant BVEC posterior draws converted to a level-VAR representation.
      real(dp), allocatable :: ar(:, :, :, :)
      real(dp), allocatable :: restricted_deterministic(:, :, :)
      real(dp), allocatable :: unrestricted_deterministic(:, :, :)
      real(dp), allocatable :: exogenous(:, :, :, :)
      integer :: retained_draws = 0
      integer :: info = 0
   end type bvartools_level_var_draws_t

   type, public :: bvartools_tvp_level_var_draws_t
      !! TVP-BVEC posterior draws converted to period-specific level VARs.
      real(dp), allocatable :: ar(:, :, :, :, :)
      real(dp), allocatable :: restricted_deterministic(:, :, :, :)
      real(dp), allocatable :: unrestricted_deterministic(:, :, :, :)
      real(dp), allocatable :: exogenous(:, :, :, :, :)
      integer :: observations = 0
      integer :: retained_draws = 0
      integer :: info = 0
   end type bvartools_tvp_level_var_draws_t

   type, public :: bvartools_var_data_t
      !! Aligned response, regressor, SUR, and holdout matrices for a BVAR.
      real(dp), allocatable :: y(:, :)
      real(dp), allocatable :: x(:, :)
      real(dp), allocatable :: sur(:, :)
      real(dp), allocatable :: tvp_sur(:, :)
      real(dp), allocatable :: structural(:, :)
      real(dp), allocatable :: holdout_y(:, :)
      real(dp), allocatable :: holdout_x(:, :)
      integer :: observations = 0
      integer :: variables = 0
      integer :: lag_order = 0
      integer :: endogenous_columns = 0
      integer :: exogenous_columns = 0
      integer :: deterministic_columns = 0
      integer :: info = 0
   end type bvartools_var_data_t

   type, public :: bvartools_vecm_data_t
      !! Aligned response, error-correction, short-run, and SUR matrices for a BVEC.
      real(dp), allocatable :: y(:, :)
      real(dp), allocatable :: w(:, :)
      real(dp), allocatable :: x(:, :)
      real(dp), allocatable :: sur(:, :)
      real(dp), allocatable :: tvp_sur(:, :)
      real(dp), allocatable :: structural(:, :)
      real(dp), allocatable :: holdout_y(:, :)
      real(dp), allocatable :: holdout_w(:, :)
      real(dp), allocatable :: holdout_x(:, :)
      integer :: observations = 0
      integer :: variables = 0
      integer :: level_lags = 0
      integer :: restricted_deterministic_columns = 0
      integer :: short_run_columns = 0
      integer :: unrestricted_deterministic_columns = 0
      integer :: info = 0
   end type bvartools_vecm_data_t

   type, public :: bvartools_dfm_data_t
      !! Standardized observations and model grid prepared for Bayesian DFM sampling.
      real(dp), allocatable :: x(:, :)
      real(dp), allocatable :: mean(:)
      real(dp), allocatable :: standard_deviation(:)
      integer, allocatable :: factor_count(:)
      integer, allocatable :: lag_order(:)
      integer :: observations = 0
      integer :: variables = 0
      integer :: iterations = 50000
      integer :: burnin = 5000
      integer :: info = 0
   end type bvartools_dfm_data_t

   public :: bvartools_minnesota_prior
   public :: bvartools_loglik_normal
   public :: bvartools_kalman_dk
   public :: bvartools_inclusion_prior
   public :: bvartools_ssvs
   public :: bvartools_ssvs_prior
   public :: bvartools_normal_posterior
   public :: bvartools_normal_draw
   public :: bvartools_sur_normal_posterior
   public :: bvartools_sur_normal_draw
   public :: bvartools_measurement_variance_posterior
   public :: bvartools_state_variance_posterior
   public :: bvartools_gamma_precision_draw
   public :: bvartools_covar_prepare_data
   public :: bvartools_covar_vector_to_matrix
   public :: bvartools_sur_const_to_tvp
   public :: bvartools_covar_const_posterior
   public :: bvartools_covar_tvp_posterior
   public :: bvartools_bvs
   public :: bvartools_random_walk_posterior
   public :: bvartools_bvar_gibbs
   public :: bvartools_structural_bvar_gibbs
   public :: bvartools_structural_impacts
   public :: bvartools_prepare_bvar_prior
   public :: bvartools_fit_bvar
   public :: bvartools_fit_structural_bvar
   public :: bvartools_cointegration_draw
   public :: bvartools_cointegration_sur_draw
   public :: bvartools_bvec_gibbs
   public :: bvartools_prepare_bvec_prior
   public :: bvartools_fit_bvec
   public :: bvartools_fit_structural_bvec
   public :: bvartools_initial_state_posterior
   public :: bvartools_tvp_bvar_gibbs
   public :: bvartools_prepare_tvp_bvar_prior
   public :: bvartools_fit_tvp_bvar
   public :: bvartools_tvp_covariance_gibbs
   public :: bvartools_joint_tvp_bvar_gibbs
   public :: bvartools_prepare_structural_tvp_bvar_prior
   public :: bvartools_fit_structural_tvp_bvar
   public :: bvartools_fit_structural_tvp_bvec
   public :: bvartools_tvp_bvs
   public :: bvartools_tvp_covariance_bvs
   public :: bvartools_stochastic_volatility_draw
   public :: bvartools_stochastic_volatility
   public :: bvartools_stochastic_volatility_ocsn_draw
   public :: bvartools_stochastic_volatility_ocsn
   public :: bvartools_tvp_bvec_gibbs
   public :: bvartools_prepare_tvp_bvec_prior
   public :: bvartools_fit_tvp_bvec
   public :: bvartools_dfm_factor_posterior
   public :: bvartools_dfm_gibbs
   public :: bvartools_dfm_prior
   public :: bvartools_dfm_grid_gibbs
   public :: bvartools_bvar_predictive_from_random
   public :: bvartools_bvar_predictive
   public :: bvartools_tvp_bvar_predictive_from_random
   public :: bvartools_tvp_bvar_predictive
   public :: bvartools_bvar_irf
   public :: bvartools_tvp_bvar_irf
   public :: bvartools_bvar_fevd
   public :: bvartools_tvp_bvar_fevd
   public :: bvartools_vecm_level_ar
   public :: bvartools_vecm_level_exogenous
   public :: bvartools_bvec_to_level_var
   public :: bvartools_tvp_bvec_to_level_var
   public :: bvartools_reconstruct_levels
   public :: bvartools_prepare_var
   public :: bvartools_prepare_vecm
   public :: bvartools_prepare_dfm
   public :: bvartools_model_comparison
   public :: bvartools_compare_models

   interface bvartools_model_comparison
      module procedure bvartools_model_comparison_constant
      module procedure bvartools_model_comparison_tvp
   end interface bvartools_model_comparison

contains

   pure function bvartools_loglik_normal(residual, covariance) result(out)
      !! Evaluate observation-wise multivariate Gaussian log likelihoods.
      real(dp), intent(in) :: residual(:, :) !! Residual vectors by variable and observation.
      real(dp), intent(in) :: covariance(:, :) !! Constant or vertically stacked observation covariances.
      type(bvartools_loglik_normal_t) :: out
      real(dp), allocatable :: zero(:)
      integer :: variables, observations, time, first

      variables = size(residual, 1)
      observations = size(residual, 2)
      if (variables < 1 .or. observations < 1 .or. &
         size(covariance, 2) /= variables .or. &
         (size(covariance, 1) /= variables .and. &
         size(covariance, 1) /= variables*observations)) then
         out%info = 1
         return
      end if
      allocate(out%value(observations), zero(variables))
      zero = 0.0_dp
      if (size(covariance, 1) == variables) then
         do time = 1, observations
            out%value(time) = multivariate_normal_log_density( &
               residual(:, time), zero, covariance)
         end do
      else
         do time = 1, observations
            first = (time - 1)*variables + 1
            out%value(time) = multivariate_normal_log_density( &
               residual(:, time), zero, &
               covariance(first:first + variables - 1, :))
         end do
      end if
      if (any(out%value <= -0.5_dp*huge(1.0_dp))) out%info = 2
   end function bvartools_loglik_normal

   function bvartools_kalman_dk(y, z, observation_covariance, &
      state_covariance, transition, initial_state, initial_covariance) &
      result(out)
      !! Draw a Durbin-Koopman state path using the shared Gaussian smoother.
      real(dp), intent(in) :: y(:, :) !! Observations by variable and time.
      real(dp), intent(in) :: z(:, :) !! Vertically stacked observation loading matrices.
      real(dp), intent(in) :: observation_covariance(:, :) !! Constant or stacked observation covariances.
      real(dp), intent(in) :: state_covariance(:, :) !! Constant or stacked state innovation covariances.
      real(dp), intent(in) :: transition(:, :) !! Constant or stacked state transition matrices.
      real(dp), intent(in) :: initial_state(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance.
      type(bvartools_kalman_dk_t) :: out
      type(ssm_model_t) :: model
      type(bssm_simulation_smoother_t) :: smoothed
      integer :: variables, observations, states, time, first

      variables = size(y, 1)
      observations = size(y, 2)
      states = size(initial_state)
      if (variables < 1 .or. observations < 1 .or. states < 1 .or. &
         any(shape(z) /= [variables*observations, states]) .or. &
         size(observation_covariance, 2) /= variables .or. &
         (size(observation_covariance, 1) /= variables .and. &
         size(observation_covariance, 1) /= variables*observations) .or. &
         size(state_covariance, 2) /= states .or. &
         (size(state_covariance, 1) /= states .and. &
         size(state_covariance, 1) /= states*observations) .or. &
         size(transition, 2) /= states .or. &
         (size(transition, 1) /= states .and. &
         size(transition, 1) /= states*observations) .or. &
         any(shape(initial_covariance) /= [states, states])) then
         out%info = 1
         return
      end if

      allocate(model%y(observations + 1, variables))
      allocate(model%z(variables, states, observations + 1))
      allocate(model%h(variables, variables, observations + 1))
      allocate(model%transition(states, states, observations + 1))
      allocate(model%r(states, states, 1))
      allocate(model%q(states, states, observations + 1))
      allocate(model%a1(states), model%p1(states, states))
      allocate(model%p1inf(states, states))
      allocate(model%missing(observations + 1, variables))
      model%y = 0.0_dp
      model%y(:observations, :) = transpose(y)
      model%z = 0.0_dp
      model%h = 0.0_dp
      model%transition = 0.0_dp
      model%q = 0.0_dp
      do time = 1, observations
         first = (time - 1)*variables + 1
         model%z(:, :, time) = z(first:first + variables - 1, :)
         if (size(observation_covariance, 1) == variables) then
            model%h(:, :, time) = observation_covariance
         else
            model%h(:, :, time) = observation_covariance( &
               first:first + variables - 1, :)
         end if
         first = (time - 1)*states + 1
         if (size(state_covariance, 1) == states) then
            model%q(:, :, time) = state_covariance
         else
            model%q(:, :, time) = state_covariance( &
               first:first + states - 1, :)
         end if
         if (size(transition, 1) == states) then
            model%transition(:, :, time) = transition
         else
            model%transition(:, :, time) = transition( &
               first:first + states - 1, :)
         end if
      end do
      model%z(:, :, observations + 1) = model%z(:, :, observations)
      model%h(:, :, observations + 1) = model%h(:, :, observations)
      model%transition(:, :, observations + 1) = &
         model%transition(:, :, observations)
      model%q(:, :, observations + 1) = model%q(:, :, observations)
      model%r(:, :, 1) = identity_matrix(states)
      model%a1 = initial_state
      model%p1 = initial_covariance
      model%p1inf = 0.0_dp
      model%missing = .false.
      model%missing(observations + 1, :) = .true.

      smoothed = bssm_simulation_smoother(model, 1)
      if (smoothed%info /= 0) then
         out%info = 10 + smoothed%info
         return
      end if
      out%state = smoothed%trajectories(:, :, 1)
   end function bvartools_kalman_dk

   pure function bvartools_minnesota_prior(residual_covariance, order, kappa0, &
      kappa1, kappa2, kappa3, exogenous_variance, exogenous_lags, &
      deterministic_count, cointegrated_var, maximum_variance) result(out)
      !! Construct the bvartools Minnesota prior for a VAR coefficient matrix.
      real(dp), intent(in) :: residual_covariance(:, :) !! Unrestricted residual covariance matrix.
      integer, intent(in) :: order !! Number of endogenous lags.
      real(dp), intent(in), optional :: kappa0 !! Own-lag variance scale.
      real(dp), intent(in), optional :: kappa1 !! Other-endogenous relative variance.
      real(dp), intent(in), optional :: kappa2 !! Exogenous relative variance.
      real(dp), intent(in), optional :: kappa3 !! Deterministic relative variance.
      real(dp), intent(in), optional :: exogenous_variance(:) !! Marginal exogenous variances.
      integer, intent(in), optional :: exogenous_lags !! Largest exogenous lag, including lag zero.
      integer, intent(in), optional :: deterministic_count !! Number of deterministic regressors.
      logical, intent(in), optional :: cointegrated_var !! Set first own-lag means to one.
      real(dp), intent(in), optional :: maximum_variance !! Upper bound on prior variances.
      type(bvartools_minnesota_prior_t) :: out
      real(dp), allocatable :: endogenous_variance(:), inverse(:, :)
      real(dp) :: own_scale, cross_scale, exogenous_scale, deterministic_scale
      integer :: k, p, m, s, d, columns, equation, predictor, lag, column, status
      logical :: unit_root_mean

      k = size(residual_covariance, 1)
      p = order
      m = 0
      if (present(exogenous_variance)) m = size(exogenous_variance)
      s = 0
      if (present(exogenous_lags)) s = exogenous_lags
      d = 0
      if (present(deterministic_count)) d = deterministic_count
      own_scale = 2.0_dp
      if (present(kappa0)) own_scale = kappa0
      cross_scale = 0.5_dp
      if (present(kappa1)) cross_scale = kappa1
      exogenous_scale = 1.0_dp
      if (present(kappa2)) exogenous_scale = kappa2
      deterministic_scale = 5.0_dp
      if (present(kappa3)) deterministic_scale = kappa3
      unit_root_mean = .false.
      if (present(cointegrated_var)) unit_root_mean = cointegrated_var
      if (k < 1 .or. size(residual_covariance, 2) /= k .or. p < 0 .or. s < 0 .or. &
         d < 0 .or. min(own_scale, cross_scale, exogenous_scale, &
         deterministic_scale) <= 0.0_dp) then
         out%info = 1
         return
      end if
      if (present(exogenous_variance)) then
         if (any(exogenous_variance <= 0.0_dp)) then
            out%info = 2
            return
         end if
      end if
      endogenous_variance = [(max(residual_covariance(equation, equation), &
         tiny(1.0_dp)), equation=1, k)]
      columns = k*p + m*(s + 1) + d
      allocate(out%mean(k, columns), out%variance(k, columns))
      out%mean = 0.0_dp
      out%variance = 0.0_dp
      do lag = 1, p
         do predictor = 1, k
            column = (lag - 1)*k + predictor
            do equation = 1, k
               if (equation == predictor) then
                  out%variance(equation, column) = own_scale/real(lag*lag, dp)
               else
                  out%variance(equation, column) = own_scale*cross_scale/ &
                     real(lag*lag, dp)*endogenous_variance(equation)/ &
                     endogenous_variance(predictor)
               end if
            end do
         end do
      end do
      if (unit_root_mean .and. p > 0) then
         do equation = 1, k
            out%mean(equation, equation) = 1.0_dp
         end do
      end if
      if (m > 0) then
         do lag = 0, s
            do predictor = 1, m
               column = k*p + lag*m + predictor
               do equation = 1, k
                  out%variance(equation, column) = own_scale*exogenous_scale/ &
                     real((lag + 1)*(lag + 1), dp)*endogenous_variance(equation)/ &
                     exogenous_variance(predictor)
               end do
            end do
         end do
      end if
      do column = k*p + m*(s + 1) + 1, columns
         out%variance(:, column) = own_scale*deterministic_scale*endogenous_variance
      end do
      if (present(maximum_variance)) then
         if (maximum_variance <= 0.0_dp) then
            out%info = 3
            return
         end if
         out%variance = min(out%variance, maximum_variance)
      end if
      call invert_matrix(residual_covariance, inverse, status)
      if (status /= 0) then
         out%info = 4
         return
      end if
      out%residual_precision = inverse
   end function bvartools_minnesota_prior

   pure function bvartools_inclusion_prior(variables, order, probability, &
      minnesota_like, kappa, exogenous_count, exogenous_lags, &
      deterministic_count, exclude_deterministics) result(out)
      !! Construct VAR coefficient inclusion probabilities used by BVS and SSVS.
      integer, intent(in) :: variables !! Number of endogenous variables.
      integer, intent(in) :: order !! Number of endogenous lags.
      real(dp), intent(in), optional :: probability !! Common prior inclusion probability.
      logical, intent(in), optional :: minnesota_like !! Use lag-decaying category probabilities.
      real(dp), intent(in), optional :: kappa(4) !! Own, cross, exogenous, and deterministic probabilities.
      integer, intent(in), optional :: exogenous_count !! Number of exogenous variables.
      integer, intent(in), optional :: exogenous_lags !! Largest exogenous lag, including lag zero.
      integer, intent(in), optional :: deterministic_count !! Number of deterministic regressors.
      logical, intent(in), optional :: exclude_deterministics !! Exclude deterministic positions from selection.
      type(bvartools_inclusion_prior_t) :: out
      real(dp) :: common, category(4)
      integer :: k, p, m, s, d, columns, lag, equation, predictor, column, position, count
      logical :: use_categories, omit_deterministics

      k = variables
      p = order
      m = 0
      if (present(exogenous_count)) m = exogenous_count
      s = 0
      if (present(exogenous_lags)) s = exogenous_lags
      d = 0
      if (present(deterministic_count)) d = deterministic_count
      common = 0.5_dp
      if (present(probability)) common = probability
      category = [0.8_dp, 0.5_dp, 0.5_dp, 0.8_dp]
      if (present(kappa)) category = kappa
      use_categories = .false.
      if (present(minnesota_like)) use_categories = minnesota_like
      omit_deterministics = .true.
      if (present(exclude_deterministics)) omit_deterministics = exclude_deterministics
      if (k < 1 .or. p < 0 .or. m < 0 .or. s < 0 .or. d < 0 .or. &
         common < 0.0_dp .or. common > 1.0_dp .or. any(category < 0.0_dp) .or. &
         any(category > 1.0_dp)) then
         out%info = 1
         return
      end if
      columns = k*p + m*(s + 1) + d
      allocate(out%probability(k, columns))
      out%probability = common
      if (use_categories) then
         do lag = 1, p
            do predictor = 1, k
               column = (lag - 1)*k + predictor
               out%probability(:, column) = category(2)/real(lag, dp)
               out%probability(predictor, column) = category(1)/real(lag, dp)
            end do
         end do
         do lag = 0, s
            do predictor = 1, m
               column = k*p + lag*m + predictor
               out%probability(:, column) = category(3)/real(lag + 1, dp)
            end do
         end do
         if (d > 0) out%probability(:, columns - d + 1:) = category(4)
      end if
      count = k*columns
      if (omit_deterministics) count = count - k*d
      allocate(out%include(count))
      position = 0
      do column = 1, columns
         if (omit_deterministics .and. column > columns - d) cycle
         do equation = 1, k
            position = position + 1
            out%include(position) = equation + (column - 1)*k
         end do
      end do
   end function bvartools_inclusion_prior

   pure function bvartools_ssvs(coefficients, tau0, tau1, prior_probability, &
      uniforms, selectable) result(out)
      !! Update George-Sun-Ni SSVS indicators from supplied uniform draws.
      real(dp), intent(in) :: coefficients(:) !! Current vectorized coefficient draw.
      real(dp), intent(in) :: tau0(:) !! Excluded-state prior standard deviations.
      real(dp), intent(in) :: tau1(:) !! Included-state prior standard deviations.
      real(dp), intent(in) :: prior_probability(:) !! Prior inclusion probabilities.
      real(dp), intent(in) :: uniforms(:) !! Uniform draws for selectable coefficients.
      integer, intent(in), optional :: selectable(:) !! One-based positions to update.
      type(bvartools_ssvs_result_t) :: out
      real(dp) :: log0, log1, maximum_log
      integer :: m, item, position, selectable_count

      m = size(coefficients)
      if (size(tau0) /= m .or. size(tau1) /= m .or. &
         size(prior_probability) /= m .or. any(tau0 <= 0.0_dp) .or. &
         any(tau1 <= 0.0_dp) .or. any(prior_probability < 0.0_dp) .or. &
         any(prior_probability > 1.0_dp)) then
         out%info = 1
         return
      end if
      if (present(selectable)) then
         if (size(uniforms) /= size(selectable) .or. any(selectable < 1) .or. &
            any(selectable > m)) then
            out%info = 2
            return
         end if
      else if (size(uniforms) /= m) then
         out%info = 3
         return
      end if
      allocate(out%posterior_probability(m), out%precision(m), out%included(m))
      out%posterior_probability = 1.0_dp
      out%precision = 1.0_dp/tau1**2
      out%included = .true.
      selectable_count = m
      if (present(selectable)) selectable_count = size(selectable)
      do item = 1, selectable_count
         if (present(selectable)) then
            position = selectable(item)
         else
            position = item
         end if
         log0 = -log(tau0(position)) - 0.5_dp*(coefficients(position)/tau0(position))**2 + &
            log(max(1.0_dp - prior_probability(position), tiny(1.0_dp)))
         log1 = -log(tau1(position)) - 0.5_dp*(coefficients(position)/tau1(position))**2 + &
            log(max(prior_probability(position), tiny(1.0_dp)))
         maximum_log = max(log0, log1)
         out%posterior_probability(position) = exp(log1 - maximum_log)/ &
            (exp(log0 - maximum_log) + exp(log1 - maximum_log))
         out%included(position) = uniforms(item) <= out%posterior_probability(position)
         if (.not. out%included(position)) out%precision(position) = 1.0_dp/tau0(position)**2
      end do
   end function bvartools_ssvs

   pure function bvartools_ssvs_prior(y, x, tau, semiautomatic, &
      covariance_count) result(out)
      !! Construct fixed or OLS-scaled SSVS prior standard deviations.
      real(dp), intent(in) :: y(:, :) !! Responses by variable and observation.
      real(dp), intent(in) :: x(:, :) !! Regressors by coefficient and observation.
      real(dp), intent(in), optional :: tau(2) !! Fixed excluded and included standard deviations.
      real(dp), intent(in), optional :: semiautomatic(2) !! OLS standard-error multipliers.
      integer, intent(in), optional :: covariance_count !! Number of contemporaneous covariance coefficients.
      type(bvartools_ssvs_prior_t) :: out
      real(dp), allocatable :: cross_product(:, :), inverse(:, :), residual(:, :)
      real(dp), allocatable :: coefficient_covariance(:, :)
      real(dp) :: fixed_scale(2), automatic_scale(2)
      integer :: variables, observations, regressors, regression_count
      integer :: covariance_parameters, total, degrees_freedom, item, status

      variables = size(y, 1)
      observations = size(y, 2)
      regressors = size(x, 1)
      fixed_scale = [0.05_dp, 10.0_dp]
      if (present(tau)) fixed_scale = tau
      covariance_parameters = 0
      if (present(covariance_count)) covariance_parameters = covariance_count
      if (variables < 1 .or. observations < 1 .or. regressors < 1 .or. &
         size(x, 2) /= observations .or. any(fixed_scale <= 0.0_dp) .or. &
         covariance_parameters < 0) then
         out%info = 1
         return
      end if
      regression_count = variables*regressors
      total = regression_count + covariance_parameters
      allocate(out%tau0(total), out%tau1(total))
      out%tau0 = fixed_scale(1)
      out%tau1 = fixed_scale(2)
      out%regression_parameters = regression_count
      out%covariance_parameters = covariance_parameters
      if (.not. present(semiautomatic)) return
      automatic_scale = semiautomatic
      if (any(automatic_scale <= 0.0_dp)) then
         out%info = 2
         return
      end if
      degrees_freedom = observations - regressors
      if (degrees_freedom <= 0) then
         out%info = 3
         return
      end if
      cross_product = matmul(x, transpose(x))
      call invert_matrix(cross_product, inverse, status)
      if (status /= 0) then
         out%info = 4
         return
      end if
      out%coefficients = matmul(matmul(y, transpose(x)), inverse)
      residual = y - matmul(out%coefficients, x)
      out%residual_covariance = matmul(residual, transpose(residual))/ &
         real(degrees_freedom, dp)
      coefficient_covariance = kronecker_product(inverse, &
         out%residual_covariance)
      allocate(out%regression_standard_error(regression_count))
      do item = 1, regression_count
         if (coefficient_covariance(item, item) < 0.0_dp) then
            out%info = 5
            return
         end if
         out%regression_standard_error(item) = &
            sqrt(coefficient_covariance(item, item))
      end do
      out%tau0(:regression_count) = &
         automatic_scale(1)*out%regression_standard_error
      out%tau1(:regression_count) = &
         automatic_scale(2)*out%regression_standard_error
      out%semiautomatic = .true.
   end function bvartools_ssvs_prior

   pure function bvartools_prepare_bvar_prior(data, coefficient_precision, &
      deterministic_precision, covariance_prior, covariance_df, &
      covariance_scale, error_shape, error_rate, ssvs_tau, &
      ssvs_semiautomatic, inclusion_probability, exclude_deterministics) &
      result(out)
      !! Dimension priors and OLS initial values for a prepared constant BVAR.
      type(bvartools_var_data_t), intent(in) :: data !! Prepared BVAR response and regressor matrices.
      real(dp), intent(in), optional :: coefficient_precision !! Regular coefficient prior precision.
      real(dp), intent(in), optional :: deterministic_precision !! Deterministic-term prior precision.
      character(len=*), intent(in), optional :: covariance_prior !! Innovation prior: wishart or gamma.
      real(dp), intent(in), optional :: covariance_df !! Inverse-Wishart degrees of freedom.
      real(dp), intent(in), optional :: covariance_scale !! Diagonal inverse-Wishart scale.
      real(dp), intent(in), optional :: error_shape !! Independent gamma precision shape.
      real(dp), intent(in), optional :: error_rate !! Independent gamma precision rate.
      real(dp), intent(in), optional :: ssvs_tau(2) !! Fixed SSVS spike and slab standard deviations.
      real(dp), intent(in), optional :: ssvs_semiautomatic(2) !! SSVS OLS standard-error multipliers.
      real(dp), intent(in), optional :: inclusion_probability !! Common SSVS inclusion probability.
      logical, intent(in), optional :: exclude_deterministics !! Exclude deterministic coefficients from SSVS.
      type(bvartools_bvar_prior_bundle_t) :: out
      type(bvartools_ssvs_prior_t) :: ssvs
      real(dp) :: regular_precision, deterministic_value, scale_value
      real(dp) :: shape_value, rate_value, probability, fixed_tau(2)
      real(dp), allocatable :: residual(:, :)
      integer :: variables, regressors, coefficients, deterministic_count
      integer :: item, selectable_count, status
      logical :: omit_deterministics

      regular_precision = 1.0_dp
      if (present(coefficient_precision)) regular_precision = coefficient_precision
      deterministic_value = 0.1_dp
      if (present(deterministic_precision)) deterministic_value = deterministic_precision
      scale_value = 1.0_dp
      if (present(covariance_scale)) scale_value = covariance_scale
      shape_value = 3.0_dp
      if (present(error_shape)) shape_value = error_shape
      rate_value = 0.0001_dp
      if (present(error_rate)) rate_value = error_rate
      probability = 0.5_dp
      if (present(inclusion_probability)) probability = inclusion_probability
      fixed_tau = [0.05_dp, 10.0_dp]
      if (present(ssvs_tau)) fixed_tau = ssvs_tau
      omit_deterministics = .false.
      if (present(exclude_deterministics)) omit_deterministics = exclude_deterministics
      if (present(covariance_prior)) out%covariance_prior = trim(covariance_prior)
      variables = data%variables
      regressors = size(data%x, 1)
      coefficients = variables*regressors
      deterministic_count = data%deterministic_columns
      if (data%info /= 0 .or. .not. allocated(data%y) .or. &
         .not. allocated(data%x) .or. variables < 1 .or. data%observations < 2 .or. &
         size(data%y, 1) /= variables .or. size(data%y, 2) /= data%observations .or. &
         size(data%x, 2) /= data%observations .or. regressors < 1 .or. &
         deterministic_count < 0 .or. deterministic_count > regressors .or. &
         regular_precision < 0.0_dp .or. deterministic_value < 0.0_dp .or. &
         scale_value <= 0.0_dp .or. shape_value < 0.0_dp .or. rate_value <= 0.0_dp .or. &
         probability < 0.0_dp .or. probability > 1.0_dp .or. &
         any(fixed_tau <= 0.0_dp) .or. &
         (trim(out%covariance_prior) /= 'wishart' .and. &
         trim(out%covariance_prior) /= 'gamma')) then
         out%info = 1
         return
      end if
      call multivariate_ols_initial(data%y, data%x, out%initial_coefficients, &
         out%initial_covariance, residual, status)
      if (status /= 0) then
         out%info = 2
         return
      end if
      allocate(out%coefficient_mean(coefficients), &
         out%coefficient_precision(coefficients, coefficients), &
         out%covariance_scale(variables, variables), out%gamma_shape(variables), &
         out%gamma_rate(variables))
      out%coefficient_mean = 0.0_dp
      out%coefficient_precision = 0.0_dp
      do item = 1, coefficients
         out%coefficient_precision(item, item) = regular_precision
      end do
      if (deterministic_count > 0) then
         do item = coefficients - variables*deterministic_count + 1, coefficients
            out%coefficient_precision(item, item) = deterministic_value
         end do
      end if
      out%covariance_scale = 0.0_dp
      do item = 1, variables
         out%covariance_scale(item, item) = scale_value
      end do
      out%covariance_df = real(variables, dp)
      if (present(covariance_df)) out%covariance_df = covariance_df
      out%gamma_shape = shape_value
      out%gamma_rate = rate_value
      if (out%covariance_df <= real(variables - 1, dp)) then
         out%info = 3
         return
      end if
      if (.not. present(ssvs_tau) .and. .not. present(ssvs_semiautomatic)) return
      if (present(ssvs_semiautomatic)) then
         ssvs = bvartools_ssvs_prior(data%y, data%x, fixed_tau, &
            ssvs_semiautomatic)
      else
         ssvs = bvartools_ssvs_prior(data%y, data%x, fixed_tau)
      end if
      if (ssvs%info /= 0) then
         out%info = 4
         return
      end if
      out%tau0 = ssvs%tau0
      out%tau1 = ssvs%tau1
      allocate(out%inclusion_probability(coefficients))
      out%inclusion_probability = probability
      selectable_count = coefficients
      if (omit_deterministics) then
         selectable_count = selectable_count - variables*deterministic_count
      end if
      allocate(out%selectable(selectable_count))
      out%selectable = [(item, item=1, selectable_count)]
      out%coefficient_precision = 0.0_dp
      do item = 1, coefficients
         out%coefficient_precision(item, item) = 1.0_dp/out%tau1(item)**2
      end do
      out%use_ssvs = .true.
   end function bvartools_prepare_bvar_prior

   pure function bvartools_prepare_bvec_prior(data, rank, &
      coefficient_precision, deterministic_precision, cointegration_shrinkage, &
      cointegration_precision, loading_precision, covariance_prior, &
      covariance_df, covariance_scale, error_shape, error_rate, ssvs_tau, &
      ssvs_semiautomatic, inclusion_probability, exclude_deterministics) &
      result(out)
      !! Dimension priors and OLS initial values for a prepared constant Bayesian VECM.
      type(bvartools_vecm_data_t), intent(in) :: data !! Prepared BVEC response and regressor matrices.
      integer, intent(in) :: rank !! Positive cointegration rank.
      real(dp), intent(in), optional :: coefficient_precision !! Regular short-run coefficient precision.
      real(dp), intent(in), optional :: deterministic_precision !! Unrestricted deterministic prior precision.
      real(dp), intent(in), optional :: cointegration_shrinkage !! Cointegration-space shrinkage parameter.
      real(dp), intent(in), optional :: cointegration_precision !! Diagonal cointegration prior precision.
      real(dp), intent(in), optional :: loading_precision !! Diagonal loading prior precision.
      character(len=*), intent(in), optional :: covariance_prior !! Innovation prior: wishart or gamma.
      real(dp), intent(in), optional :: covariance_df !! Inverse-Wishart degrees of freedom.
      real(dp), intent(in), optional :: covariance_scale !! Diagonal inverse-Wishart scale metadata.
      real(dp), intent(in), optional :: error_shape !! Independent gamma precision shape.
      real(dp), intent(in), optional :: error_rate !! Independent gamma precision rate.
      real(dp), intent(in), optional :: ssvs_tau(2) !! Fixed SSVS spike and slab standard deviations.
      real(dp), intent(in), optional :: ssvs_semiautomatic(2) !! Conditional OLS standard-error multipliers.
      real(dp), intent(in), optional :: inclusion_probability !! Common SSVS inclusion probability.
      logical, intent(in), optional :: exclude_deterministics !! Exclude unrestricted deterministic coefficients.
      type(bvartools_bvec_prior_bundle_t) :: out
      type(bvartools_ssvs_prior_t) :: full_ssvs
      real(dp), allocatable :: combined(:, :), ols_coefficients(:, :), residual(:, :)
      real(dp) :: regular_precision, deterministic_value, coint_precision
      real(dp) :: loading_value, scale_value, shape_value, rate_value
      real(dp) :: probability, fixed_tau(2)
      integer :: variables, w_count, x_count, gamma_count, full_count
      integer :: deterministic_count, item, omitted, selectable_count, status
      logical :: omit_deterministics

      regular_precision = 1.0_dp
      if (present(coefficient_precision)) regular_precision = coefficient_precision
      deterministic_value = 0.1_dp
      if (present(deterministic_precision)) deterministic_value = deterministic_precision
      out%shrinkage = 0.0_dp
      if (present(cointegration_shrinkage)) out%shrinkage = cointegration_shrinkage
      coint_precision = 1.0_dp
      if (present(cointegration_precision)) coint_precision = cointegration_precision
      loading_value = 1.0_dp
      if (present(loading_precision)) loading_value = loading_precision
      scale_value = 1.0_dp
      if (present(covariance_scale)) scale_value = covariance_scale
      shape_value = 3.0_dp
      if (present(error_shape)) shape_value = error_shape
      rate_value = 0.0001_dp
      if (present(error_rate)) rate_value = error_rate
      probability = 0.5_dp
      if (present(inclusion_probability)) probability = inclusion_probability
      fixed_tau = [0.05_dp, 10.0_dp]
      if (present(ssvs_tau)) fixed_tau = ssvs_tau
      omit_deterministics = .false.
      if (present(exclude_deterministics)) omit_deterministics = exclude_deterministics
      if (present(covariance_prior)) out%covariance_prior = trim(covariance_prior)
      variables = data%variables
      w_count = size(data%w, 1)
      x_count = size(data%x, 1)
      gamma_count = variables*x_count
      deterministic_count = data%unrestricted_deterministic_columns
      if (data%info /= 0 .or. .not. allocated(data%y) .or. &
         .not. allocated(data%w) .or. .not. allocated(data%x) .or. &
         variables < 1 .or. data%observations < 2 .or. rank < 0 .or. &
         rank > min(variables, w_count) .or. size(data%y, 1) /= variables .or. &
         size(data%y, 2) /= data%observations .or. &
         size(data%w, 2) /= data%observations .or. &
         size(data%x, 2) /= data%observations .or. deterministic_count < 0 .or. &
         deterministic_count > x_count .or. regular_precision < 0.0_dp .or. &
         deterministic_value < 0.0_dp .or. out%shrinkage < 0.0_dp .or. &
         coint_precision < 0.0_dp .or. loading_value < 0.0_dp .or. &
         scale_value <= 0.0_dp .or. shape_value < 0.0_dp .or. &
         rate_value <= 0.0_dp .or. probability < 0.0_dp .or. &
         probability > 1.0_dp .or. any(fixed_tau <= 0.0_dp) .or. &
         (trim(out%covariance_prior) /= 'wishart' .and. &
         trim(out%covariance_prior) /= 'gamma')) then
         out%info = 1
         return
      end if
      if (rank == 0 .and. x_count == 0) then
         out%info = 1
         return
      end if
      if (rank > 0) then
         allocate(combined(w_count + x_count, data%observations))
         combined(:w_count, :) = data%w
         if (x_count > 0) combined(w_count + 1:, :) = data%x
      else
         allocate(combined(x_count, data%observations))
         combined = data%x
      end if
      call multivariate_ols_initial(data%y, combined, ols_coefficients, &
         out%initial_covariance, residual, status)
      if (status /= 0) then
         out%info = 2
         return
      end if
      allocate(out%initial_beta(w_count, rank), out%initial_gamma(variables, x_count), &
         out%cointegration_precision(w_count, w_count), &
         out%loading_precision(variables, variables), out%gamma_mean(gamma_count), &
         out%gamma_precision(gamma_count, gamma_count), &
         out%covariance_scale(variables, variables), out%error_shape(variables), &
         out%error_rate(variables))
      out%initial_beta = 0.0_dp
      do item = 1, rank
         out%initial_beta(item, item) = 1.0_dp
      end do
      if (x_count > 0) then
         if (rank > 0) then
            out%initial_gamma = ols_coefficients(:, w_count + 1:)
         else
            out%initial_gamma = ols_coefficients
         end if
      end if
      out%cointegration_precision = 0.0_dp
      out%loading_precision = 0.0_dp
      do item = 1, w_count
         out%cointegration_precision(item, item) = coint_precision
      end do
      do item = 1, variables
         out%loading_precision(item, item) = loading_value
      end do
      out%gamma_mean = 0.0_dp
      out%gamma_precision = 0.0_dp
      do item = 1, gamma_count
         out%gamma_precision(item, item) = regular_precision
      end do
      if (deterministic_count > 0) then
         do item = gamma_count - variables*deterministic_count + 1, gamma_count
            out%gamma_precision(item, item) = deterministic_value
         end do
      end if
      out%covariance_scale = 0.0_dp
      do item = 1, variables
         out%covariance_scale(item, item) = scale_value
      end do
      out%covariance_df = real(variables + rank, dp)
      if (present(covariance_df)) out%covariance_df = covariance_df
      out%error_shape = shape_value + real(rank, dp)
      out%error_rate = rate_value
      out%rank = rank
      if (out%covariance_df <= real(variables - 1, dp)) then
         out%info = 3
         return
      end if
      if (.not. present(ssvs_tau) .and. .not. present(ssvs_semiautomatic)) return
      if (gamma_count < 1) then
         out%info = 4
         return
      end if
      if (present(ssvs_semiautomatic)) then
         full_ssvs = bvartools_ssvs_prior(data%y, combined, fixed_tau, &
            ssvs_semiautomatic)
      else
         full_ssvs = bvartools_ssvs_prior(data%y, combined, fixed_tau)
      end if
      if (full_ssvs%info /= 0) then
         out%info = 5
         return
      end if
      omitted = 0
      if (rank > 0) omitted = variables*w_count
      full_count = omitted + gamma_count
      out%tau0 = full_ssvs%tau0(omitted + 1:full_count)
      out%tau1 = full_ssvs%tau1(omitted + 1:full_count)
      allocate(out%inclusion_probability(gamma_count))
      out%inclusion_probability = probability
      selectable_count = gamma_count
      if (omit_deterministics) then
         selectable_count = selectable_count - variables*deterministic_count
      end if
      allocate(out%selectable(selectable_count))
      out%selectable = [(item, item=1, selectable_count)]
      out%gamma_precision = 0.0_dp
      do item = 1, gamma_count
         out%gamma_precision(item, item) = 1.0_dp/out%tau1(item)**2
      end do
      out%use_ssvs = .true.
   end function bvartools_prepare_bvec_prior

   pure subroutine multivariate_ols_initial(y, x, coefficients, covariance, &
      residual, info)
      !! Compute multivariate OLS coefficients and maximum-likelihood covariance.
      real(dp), intent(in) :: y(:, :) !! Responses by variable and observation.
      real(dp), intent(in) :: x(:, :) !! Regressors by coefficient and observation.
      real(dp), allocatable, intent(out) :: coefficients(:, :) !! OLS coefficient matrix.
      real(dp), allocatable, intent(out) :: covariance(:, :) !! Residual covariance matrix.
      real(dp), allocatable, intent(out) :: residual(:, :) !! OLS residual matrix.
      integer, intent(out) :: info !! Zero for a successful full-rank regression.
      real(dp), allocatable :: cross_product(:, :), inverse(:, :)
      integer :: observations

      observations = size(y, 2)
      if (size(x, 2) /= observations .or. size(y, 1) < 1 .or. &
         size(x, 1) < 1 .or. observations <= size(x, 1)) then
         info = 1
         return
      end if
      cross_product = matmul(x, transpose(x))
      call invert_matrix(cross_product, inverse, info)
      if (info /= 0) return
      coefficients = matmul(matmul(y, transpose(x)), inverse)
      residual = y - matmul(coefficients, x)
      covariance = matmul(residual, transpose(residual))/real(observations, dp)
   end subroutine multivariate_ols_initial

   function bvartools_fit_bvar(data, prior, iterations, burnin) result(out)
      !! Fit a prepared constant BVAR using a compatible prior bundle.
      type(bvartools_var_data_t), intent(in) :: data !! Prepared BVAR response and regressors.
      type(bvartools_bvar_prior_bundle_t), intent(in) :: prior !! Dimensioned BVAR prior bundle.
      integer, intent(in) :: iterations !! Number of retained Gibbs draws.
      integer, intent(in) :: burnin !! Number of discarded Gibbs draws.
      type(bvartools_bvar_fit_t) :: out
      integer :: coefficients

      if (data%info /= 0 .or. prior%info /= 0 .or. iterations < 1 .or. burnin < 0 .or. &
         .not. allocated(data%y) .or. .not. allocated(data%x) .or. &
         .not. allocated(prior%coefficient_mean) .or. &
         .not. allocated(prior%coefficient_precision) .or. &
         .not. allocated(prior%initial_covariance)) then
         out%info = 1
         return
      end if
      coefficients = data%variables*size(data%x, 1)
      if (size(prior%coefficient_mean) /= coefficients .or. &
         any(shape(prior%coefficient_precision) /= [coefficients, coefficients]) .or. &
         any(shape(prior%initial_covariance) /= [data%variables, data%variables])) then
         out%info = 2
         return
      end if
      if (prior%use_ssvs) then
         if (.not. allocated(prior%tau0) .or. .not. allocated(prior%tau1) .or. &
            .not. allocated(prior%inclusion_probability) .or. &
            .not. allocated(prior%selectable)) then
            out%info = 3
            return
         end if
      end if
      select case (trim(prior%covariance_prior))
      case ('wishart')
         if (.not. allocated(prior%covariance_scale)) then
            out%info = 4
            return
         end if
         if (prior%use_ssvs) then
            out%draws = bvartools_bvar_gibbs(data%y, data%x, &
               prior%coefficient_mean, prior%coefficient_precision, &
               prior%initial_covariance, iterations, burnin, &
               prior%covariance_scale, prior%covariance_df, tau0=prior%tau0, &
               tau1=prior%tau1, inclusion_probability=prior%inclusion_probability, &
               selectable=prior%selectable)
         else
            out%draws = bvartools_bvar_gibbs(data%y, data%x, &
               prior%coefficient_mean, prior%coefficient_precision, &
               prior%initial_covariance, iterations, burnin, &
               prior%covariance_scale, prior%covariance_df)
         end if
      case ('gamma')
         if (.not. allocated(prior%gamma_shape) .or. &
            .not. allocated(prior%gamma_rate)) then
            out%info = 4
            return
         end if
         if (prior%use_ssvs) then
            out%draws = bvartools_bvar_gibbs(data%y, data%x, &
               prior%coefficient_mean, prior%coefficient_precision, &
               prior%initial_covariance, iterations, burnin, &
               gamma_shape=prior%gamma_shape, gamma_rate=prior%gamma_rate, &
               tau0=prior%tau0, tau1=prior%tau1, &
               inclusion_probability=prior%inclusion_probability, &
               selectable=prior%selectable)
         else
            out%draws = bvartools_bvar_gibbs(data%y, data%x, &
               prior%coefficient_mean, prior%coefficient_precision, &
               prior%initial_covariance, iterations, burnin, &
               gamma_shape=prior%gamma_shape, gamma_rate=prior%gamma_rate)
         end if
      case default
         out%info = 4
         return
      end select
      if (out%draws%info /= 0) then
         out%info = 10 + out%draws%info
         return
      end if
      out%observations = data%observations
      out%variables = data%variables
      out%lag_order = data%lag_order
   end function bvartools_fit_bvar

   function bvartools_fit_bvec(data, prior, iterations, burnin) result(out)
      !! Fit a prepared constant BVEC model, routing rank zero through BVAR sampling.
      type(bvartools_vecm_data_t), intent(in) :: data !! Prepared BVEC response and regressors.
      type(bvartools_bvec_prior_bundle_t), intent(in) :: prior !! Dimensioned BVEC prior bundle.
      integer, intent(in) :: iterations !! Number of retained Gibbs draws.
      integer, intent(in) :: burnin !! Number of discarded Gibbs draws.
      type(bvartools_bvec_fit_t) :: out
      integer :: gamma_count

      if (data%info /= 0 .or. prior%info /= 0 .or. iterations < 1 .or. burnin < 0 .or. &
         .not. allocated(data%y) .or. .not. allocated(data%w) .or. &
         .not. allocated(data%x) .or. .not. allocated(prior%initial_beta) .or. &
         .not. allocated(prior%initial_covariance) .or. &
         .not. allocated(prior%cointegration_precision) .or. &
         .not. allocated(prior%loading_precision) .or. &
         .not. allocated(prior%gamma_mean) .or. &
         .not. allocated(prior%gamma_precision)) then
         out%info = 1
         return
      end if
      gamma_count = data%variables*size(data%x, 1)
      if (prior%rank < 0 .or. &
         any(shape(prior%initial_beta) /= [size(data%w, 1), prior%rank]) .or. &
         any(shape(prior%initial_covariance) /= [data%variables, data%variables]) .or. &
         size(prior%gamma_mean) /= gamma_count .or. &
         any(shape(prior%gamma_precision) /= [gamma_count, gamma_count])) then
         out%info = 2
         return
      end if
      if (prior%use_ssvs) then
         if (.not. allocated(prior%tau0) .or. .not. allocated(prior%tau1) .or. &
            .not. allocated(prior%inclusion_probability) .or. &
            .not. allocated(prior%selectable)) then
            out%info = 3
            return
         end if
      end if
      if (prior%rank == 0) then
         call fit_rank_zero_bvec(data, prior, iterations, burnin, &
            out%rank_zero_draws, out%info)
         out%rank_zero = .true.
      else
         call fit_positive_rank_bvec(data, prior, iterations, burnin, &
            out%draws, out%info)
      end if
      if (out%info /= 0) return
      out%observations = data%observations
      out%variables = data%variables
      out%level_lags = data%level_lags
      out%rank = prior%rank
   end function bvartools_fit_bvec

   function bvartools_fit_structural_bvar(data, prior, iterations, burnin, &
      structural_prior_precision, structural_tau0, structural_tau1, &
      structural_inclusion_probability) result(out)
      !! Fit a prepared constant structural BVAR using an existing prior bundle.
      type(bvartools_var_data_t), intent(in) :: data !! Prepared structural BVAR data.
      type(bvartools_bvar_prior_bundle_t), intent(in) :: prior !! Reduced-form BVAR prior bundle.
      integer, intent(in) :: iterations !! Number of retained Gibbs draws.
      integer, intent(in) :: burnin !! Number of discarded Gibbs draws.
      real(dp), intent(in), optional :: structural_prior_precision !! Structural prior precision.
      real(dp), intent(in), optional :: structural_tau0 !! Structural SSVS spike deviation.
      real(dp), intent(in), optional :: structural_tau1 !! Structural SSVS slab deviation.
      real(dp), intent(in), optional :: structural_inclusion_probability !! Structural SSVS probability.
      type(bvartools_bvar_fit_t) :: out
      real(dp), allocatable :: reduced_tau0(:), reduced_tau1(:), reduced_probability(:)
      real(dp), allocatable :: structural_spike(:), structural_slab(:)
      real(dp), allocatable :: structural_probability(:)
      integer, allocatable :: reduced_selectable(:), structural_selectable(:)
      real(dp) :: precision_value, spike, slab, probability
      integer :: reduced, structural, item

      precision_value = 1.0_dp
      if (present(structural_prior_precision)) precision_value = &
         structural_prior_precision
      spike = 0.05_dp
      if (present(structural_tau0)) spike = structural_tau0
      slab = 10.0_dp
      if (present(structural_tau1)) slab = structural_tau1
      probability = 0.5_dp
      if (present(structural_inclusion_probability)) probability = &
         structural_inclusion_probability
      reduced = size(prior%coefficient_mean)
      structural = data%variables*(data%variables - 1)/2
      if (data%info /= 0 .or. prior%info /= 0 .or. data%variables < 2 .or. &
         .not. allocated(data%structural) .or. &
         any(shape(data%structural) /= &
         [data%variables*data%observations, structural]) .or. &
         precision_value <= 0.0_dp .or. spike <= 0.0_dp .or. slab <= 0.0_dp .or. &
         probability < 0.0_dp .or. probability > 1.0_dp) then
         out%info = 1
         return
      end if
      allocate(reduced_tau0(reduced), reduced_tau1(reduced), &
         reduced_probability(reduced), structural_spike(structural), &
         structural_slab(structural), structural_probability(structural))
      if (prior%use_ssvs) then
         reduced_tau0 = prior%tau0
         reduced_tau1 = prior%tau1
         reduced_probability = prior%inclusion_probability
         reduced_selectable = prior%selectable
      else
         do item = 1, reduced
            reduced_tau1(item) = 1.0_dp/sqrt(max( &
               prior%coefficient_precision(item, item), tiny(1.0_dp)))
         end do
         reduced_tau0 = reduced_tau1
         reduced_probability = 0.5_dp
         allocate(reduced_selectable(0))
      end if
      structural_spike = spike
      structural_slab = slab
      structural_probability = probability
      if (present(structural_inclusion_probability)) then
         structural_selectable = [(item, item=1, structural)]
      else
         structural_spike = 1.0_dp/sqrt(precision_value)
         structural_slab = structural_spike
         allocate(structural_selectable(0))
      end if
      select case (trim(prior%covariance_prior))
      case ('wishart')
         out%draws = bvartools_structural_bvar_gibbs(data%y, data%sur, &
            data%structural, prior%coefficient_mean, prior%coefficient_precision, &
            spread(0.0_dp, 1, structural), &
            precision_value*identity_matrix(structural), prior%initial_covariance, &
            iterations, burnin, prior%covariance_scale, prior%covariance_df, &
            reduced_tau0=reduced_tau0, reduced_tau1=reduced_tau1, &
            reduced_inclusion_probability=reduced_probability, &
            reduced_selectable=reduced_selectable, structural_tau0=structural_spike, &
            structural_tau1=structural_slab, &
            structural_inclusion_probability=structural_probability, &
            structural_selectable=structural_selectable)
      case ('gamma')
         out%draws = bvartools_structural_bvar_gibbs(data%y, data%sur, &
            data%structural, prior%coefficient_mean, prior%coefficient_precision, &
            spread(0.0_dp, 1, structural), &
            precision_value*identity_matrix(structural), prior%initial_covariance, &
            iterations, burnin, gamma_shape=prior%gamma_shape, &
            gamma_rate=prior%gamma_rate, reduced_tau0=reduced_tau0, &
            reduced_tau1=reduced_tau1, &
            reduced_inclusion_probability=reduced_probability, &
            reduced_selectable=reduced_selectable, structural_tau0=structural_spike, &
            structural_tau1=structural_slab, &
            structural_inclusion_probability=structural_probability, &
            structural_selectable=structural_selectable)
      case default
         out%info = 2
         return
      end select
      if (out%draws%info /= 0) then
         out%info = 10 + out%draws%info
         return
      end if
      out%observations = data%observations
      out%variables = data%variables
      out%lag_order = data%lag_order
   end function bvartools_fit_structural_bvar

   function bvartools_fit_structural_bvec(data, prior, iterations, burnin, &
      structural_prior_precision, structural_tau0, structural_tau1, &
      structural_inclusion_probability) result(out)
      !! Fit a prepared positive-rank constant structural BVEC model.
      type(bvartools_vecm_data_t), intent(in) :: data !! Prepared structural BVEC data.
      type(bvartools_bvec_prior_bundle_t), intent(in) :: prior !! Reduced-form BVEC prior bundle.
      integer, intent(in) :: iterations !! Number of retained Gibbs draws.
      integer, intent(in) :: burnin !! Number of discarded Gibbs draws.
      real(dp), intent(in), optional :: structural_prior_precision !! Structural prior precision.
      real(dp), intent(in), optional :: structural_tau0 !! Structural SSVS spike deviation.
      real(dp), intent(in), optional :: structural_tau1 !! Structural SSVS slab deviation.
      real(dp), intent(in), optional :: structural_inclusion_probability !! Structural SSVS probability.
      type(bvartools_bvec_fit_t) :: out
      real(dp), allocatable :: gamma_tau0(:), gamma_tau1(:), gamma_probability(:)
      integer, allocatable :: gamma_selectable(:), structural_selectable(:)
      real(dp) :: precision_value, spike, slab, probability
      integer :: gamma_count, structural, item

      precision_value = 1.0_dp
      if (present(structural_prior_precision)) precision_value = &
         structural_prior_precision
      spike = 0.05_dp
      if (present(structural_tau0)) spike = structural_tau0
      slab = 10.0_dp
      if (present(structural_tau1)) slab = structural_tau1
      probability = 0.5_dp
      if (present(structural_inclusion_probability)) probability = &
         structural_inclusion_probability
      gamma_count = size(prior%gamma_mean)
      structural = data%variables*(data%variables - 1)/2
      if (data%info /= 0 .or. prior%info /= 0 .or. prior%rank < 1 .or. &
         gamma_count < 1 .or. .not. allocated(data%structural) .or. &
         precision_value <= 0.0_dp .or. spike <= 0.0_dp .or. slab <= 0.0_dp .or. &
         probability < 0.0_dp .or. probability > 1.0_dp) then
         out%info = 1
         return
      end if
      allocate(gamma_tau0(gamma_count), gamma_tau1(gamma_count), &
         gamma_probability(gamma_count))
      if (prior%use_ssvs) then
         gamma_tau0 = prior%tau0
         gamma_tau1 = prior%tau1
         gamma_probability = prior%inclusion_probability
         gamma_selectable = prior%selectable
      else
         do item = 1, gamma_count
            gamma_tau1(item) = 1.0_dp/sqrt(max( &
               prior%gamma_precision(item, item), tiny(1.0_dp)))
         end do
         gamma_tau0 = gamma_tau1
         gamma_probability = 0.5_dp
         allocate(gamma_selectable(0))
      end if
      if (present(structural_inclusion_probability)) then
         structural_selectable = [(item, item=1, structural)]
      else
         spike = 1.0_dp/sqrt(precision_value)
         slab = spike
         allocate(structural_selectable(0))
      end if
      if (trim(prior%covariance_prior) == 'wishart') then
         out%draws = bvartools_bvec_gibbs(data%y, data%w, prior%initial_beta, &
            prior%initial_covariance, iterations, burnin, prior%shrinkage, &
            prior%cointegration_precision, prior%loading_precision, &
            covariance_prior_df=prior%covariance_df, x=data%x, &
            gamma_prior_mean=prior%gamma_mean, &
            gamma_prior_precision=prior%gamma_precision, tau0=gamma_tau0, &
            tau1=gamma_tau1, inclusion_probability=gamma_probability, &
            selectable=gamma_selectable, structural_design=data%structural, &
            structural_prior_mean=spread(0.0_dp, 1, structural), &
            structural_prior_precision=precision_value*identity_matrix(structural), &
            structural_tau0=spread(spike, 1, structural), &
            structural_tau1=spread(slab, 1, structural), &
            structural_inclusion_probability=spread(probability, 1, structural), &
            structural_selectable=structural_selectable)
      else
         out%draws = bvartools_bvec_gibbs(data%y, data%w, prior%initial_beta, &
            prior%initial_covariance, iterations, burnin, prior%shrinkage, &
            prior%cointegration_precision, prior%loading_precision, x=data%x, &
            gamma_prior_mean=prior%gamma_mean, &
            gamma_prior_precision=prior%gamma_precision, &
            gamma_shape=prior%error_shape, gamma_rate=prior%error_rate, &
            tau0=gamma_tau0, tau1=gamma_tau1, &
            inclusion_probability=gamma_probability, selectable=gamma_selectable, &
            structural_design=data%structural, &
            structural_prior_mean=spread(0.0_dp, 1, structural), &
            structural_prior_precision=precision_value*identity_matrix(structural), &
            structural_tau0=spread(spike, 1, structural), &
            structural_tau1=spread(slab, 1, structural), &
            structural_inclusion_probability=spread(probability, 1, structural), &
            structural_selectable=structural_selectable)
      end if
      if (out%draws%info /= 0) then
         out%info = 10 + out%draws%info
         return
      end if
      out%observations = data%observations
      out%variables = data%variables
      out%level_lags = data%level_lags
      out%rank = prior%rank
   end function bvartools_fit_structural_bvec

   subroutine fit_rank_zero_bvec(data, prior, iterations, burnin, draws, info)
      !! Sample a rank-zero BVEC specification as an unrestricted differenced BVAR.
      type(bvartools_vecm_data_t), intent(in) :: data !! Prepared BVEC response and regressors.
      type(bvartools_bvec_prior_bundle_t), intent(in) :: prior !! Rank-zero BVEC prior bundle.
      integer, intent(in) :: iterations !! Number of retained Gibbs draws.
      integer, intent(in) :: burnin !! Number of discarded Gibbs draws.
      type(bvartools_bvar_draws_t), intent(out) :: draws !! Rank-zero posterior draws.
      integer, intent(out) :: info !! Wrapper status code.

      info = 0
      select case (trim(prior%covariance_prior))
      case ('wishart')
         if (.not. allocated(prior%covariance_scale)) then
            info = 4
            return
         end if
         if (prior%use_ssvs) then
            draws = bvartools_bvar_gibbs(data%y, data%x, prior%gamma_mean, &
               prior%gamma_precision, prior%initial_covariance, iterations, burnin, &
               prior%covariance_scale, prior%covariance_df, tau0=prior%tau0, &
               tau1=prior%tau1, inclusion_probability=prior%inclusion_probability, &
               selectable=prior%selectable)
         else
            draws = bvartools_bvar_gibbs(data%y, data%x, prior%gamma_mean, &
               prior%gamma_precision, prior%initial_covariance, iterations, burnin, &
               prior%covariance_scale, prior%covariance_df)
         end if
      case ('gamma')
         if (.not. allocated(prior%error_shape) .or. &
            .not. allocated(prior%error_rate)) then
            info = 4
            return
         end if
         if (prior%use_ssvs) then
            draws = bvartools_bvar_gibbs(data%y, data%x, prior%gamma_mean, &
               prior%gamma_precision, prior%initial_covariance, iterations, burnin, &
               gamma_shape=prior%error_shape, gamma_rate=prior%error_rate, &
               tau0=prior%tau0, tau1=prior%tau1, &
               inclusion_probability=prior%inclusion_probability, &
               selectable=prior%selectable)
         else
            draws = bvartools_bvar_gibbs(data%y, data%x, prior%gamma_mean, &
               prior%gamma_precision, prior%initial_covariance, iterations, burnin, &
               gamma_shape=prior%error_shape, gamma_rate=prior%error_rate)
         end if
      case default
         info = 4
         return
      end select
      if (draws%info /= 0) info = 10 + draws%info
   end subroutine fit_rank_zero_bvec

   subroutine fit_positive_rank_bvec(data, prior, iterations, burnin, draws, info)
      !! Sample a positive-rank constant Bayesian VECM from a prior bundle.
      type(bvartools_vecm_data_t), intent(in) :: data !! Prepared BVEC response and regressors.
      type(bvartools_bvec_prior_bundle_t), intent(in) :: prior !! Positive-rank BVEC prior bundle.
      integer, intent(in) :: iterations !! Number of retained Gibbs draws.
      integer, intent(in) :: burnin !! Number of discarded Gibbs draws.
      type(bvartools_bvec_draws_t), intent(out) :: draws !! Positive-rank posterior draws.
      integer, intent(out) :: info !! Wrapper status code.
      logical :: has_x

      info = 0
      has_x = size(data%x, 1) > 0
      select case (trim(prior%covariance_prior))
      case ('wishart')
         if (has_x) then
            if (prior%use_ssvs) then
               draws = bvartools_bvec_gibbs(data%y, data%w, prior%initial_beta, &
                  prior%initial_covariance, iterations, burnin, prior%shrinkage, &
                  prior%cointegration_precision, prior%loading_precision, &
                  covariance_prior_df=prior%covariance_df, x=data%x, &
                  gamma_prior_mean=prior%gamma_mean, &
                  gamma_prior_precision=prior%gamma_precision, tau0=prior%tau0, &
                  tau1=prior%tau1, &
                  inclusion_probability=prior%inclusion_probability, &
                  selectable=prior%selectable)
            else
               draws = bvartools_bvec_gibbs(data%y, data%w, prior%initial_beta, &
                  prior%initial_covariance, iterations, burnin, prior%shrinkage, &
                  prior%cointegration_precision, prior%loading_precision, &
                  covariance_prior_df=prior%covariance_df, x=data%x, &
                  gamma_prior_mean=prior%gamma_mean, &
                  gamma_prior_precision=prior%gamma_precision)
            end if
         else
            draws = bvartools_bvec_gibbs(data%y, data%w, prior%initial_beta, &
               prior%initial_covariance, iterations, burnin, prior%shrinkage, &
               prior%cointegration_precision, prior%loading_precision, &
               covariance_prior_df=prior%covariance_df)
         end if
      case ('gamma')
         if (.not. allocated(prior%error_shape) .or. &
            .not. allocated(prior%error_rate)) then
            info = 4
            return
         end if
         if (has_x) then
            if (prior%use_ssvs) then
               draws = bvartools_bvec_gibbs(data%y, data%w, prior%initial_beta, &
                  prior%initial_covariance, iterations, burnin, prior%shrinkage, &
                  prior%cointegration_precision, prior%loading_precision, &
                  x=data%x, gamma_prior_mean=prior%gamma_mean, &
                  gamma_prior_precision=prior%gamma_precision, &
                  gamma_shape=prior%error_shape, gamma_rate=prior%error_rate, &
                  tau0=prior%tau0, tau1=prior%tau1, &
                  inclusion_probability=prior%inclusion_probability, &
                  selectable=prior%selectable)
            else
               draws = bvartools_bvec_gibbs(data%y, data%w, prior%initial_beta, &
                  prior%initial_covariance, iterations, burnin, prior%shrinkage, &
                  prior%cointegration_precision, prior%loading_precision, &
                  x=data%x, gamma_prior_mean=prior%gamma_mean, &
                  gamma_prior_precision=prior%gamma_precision, &
                  gamma_shape=prior%error_shape, gamma_rate=prior%error_rate)
            end if
         else
            draws = bvartools_bvec_gibbs(data%y, data%w, prior%initial_beta, &
               prior%initial_covariance, iterations, burnin, prior%shrinkage, &
               prior%cointegration_precision, prior%loading_precision, &
               gamma_shape=prior%error_shape, gamma_rate=prior%error_rate)
         end if
      case default
         info = 4
         return
      end select
      if (draws%info /= 0) info = 10 + draws%info
   end subroutine fit_positive_rank_bvec

   pure function bvartools_prepare_tvp_bvar_prior(data, observation_prior, &
      initial_state_precision, state_shape, state_rate, deterministic_state_rate, &
      covariance_df, covariance_scale, measurement_shape, measurement_rate, &
      bvs_probability, exclude_deterministics, volatility_state_variance, &
      volatility_shape, volatility_rate, volatility_initial_precision, &
      volatility_offset, volatility_method) result(out)
      !! Dimension state, observation, BVS, and volatility priors for a TVP-BVAR.
      type(bvartools_var_data_t), intent(in) :: data !! Prepared BVAR response, regressors, and SUR design.
      character(len=*), intent(in), optional :: observation_prior !! Observation mode: wishart, gamma, or sv.
      real(dp), intent(in), optional :: initial_state_precision !! Diagonal initial coefficient-state precision.
      real(dp), intent(in), optional :: state_shape !! Coefficient-state innovation precision shape.
      real(dp), intent(in), optional :: state_rate !! Coefficient-state innovation precision rate.
      real(dp), intent(in), optional :: deterministic_state_rate !! Deterministic coefficient-state rate.
      real(dp), intent(in), optional :: covariance_df !! Inverse-Wishart degrees of freedom.
      real(dp), intent(in), optional :: covariance_scale !! Diagonal inverse-Wishart scale.
      real(dp), intent(in), optional :: measurement_shape !! Diagonal observation precision shape.
      real(dp), intent(in), optional :: measurement_rate !! Diagonal observation precision rate.
      real(dp), intent(in), optional :: bvs_probability !! Common coefficient-trajectory inclusion probability.
      logical, intent(in), optional :: exclude_deterministics !! Exclude deterministic trajectories from BVS.
      real(dp), intent(in), optional :: volatility_state_variance !! Initial log-variance innovation variance.
      real(dp), intent(in), optional :: volatility_shape !! Log-variance innovation precision shape.
      real(dp), intent(in), optional :: volatility_rate !! Log-variance innovation precision rate.
      real(dp), intent(in), optional :: volatility_initial_precision !! Initial log-variance prior precision.
      real(dp), intent(in), optional :: volatility_offset !! Positive log-square transformation offset.
      character(len=*), intent(in), optional :: volatility_method !! Mixture method: ksc or ocsn.
      type(bvartools_tvp_bvar_prior_bundle_t) :: out
      real(dp), allocatable :: coefficients(:, :), residual(:, :)
      real(dp) :: initial_precision_value, state_shape_value, state_rate_value
      real(dp) :: deterministic_rate_value, scale_value, measurement_shape_value
      real(dp) :: measurement_rate_value, probability, volatility_variance_value
      real(dp) :: volatility_shape_value, volatility_rate_value
      real(dp) :: volatility_precision_value, offset_value
      integer :: variables, states, deterministic_states, selectable_count
      integer :: item, time, status
      logical :: omit_deterministics

      if (present(observation_prior)) out%observation_prior = trim(observation_prior)
      if (present(volatility_method)) out%volatility_method = trim(volatility_method)
      initial_precision_value = 1.0_dp
      if (present(initial_state_precision)) initial_precision_value = initial_state_precision
      state_shape_value = 3.0_dp
      if (present(state_shape)) state_shape_value = state_shape
      state_rate_value = 0.0001_dp
      if (present(state_rate)) state_rate_value = state_rate
      deterministic_rate_value = 0.01_dp
      if (present(deterministic_state_rate)) deterministic_rate_value = deterministic_state_rate
      scale_value = 1.0_dp
      if (present(covariance_scale)) scale_value = covariance_scale
      measurement_shape_value = 3.0_dp
      if (present(measurement_shape)) measurement_shape_value = measurement_shape
      measurement_rate_value = 0.0001_dp
      if (present(measurement_rate)) measurement_rate_value = measurement_rate
      probability = 0.5_dp
      if (present(bvs_probability)) probability = bvs_probability
      volatility_variance_value = 0.05_dp
      if (present(volatility_state_variance)) then
         volatility_variance_value = volatility_state_variance
      end if
      volatility_shape_value = 3.0_dp
      if (present(volatility_shape)) volatility_shape_value = volatility_shape
      volatility_rate_value = 0.0001_dp
      if (present(volatility_rate)) volatility_rate_value = volatility_rate
      volatility_precision_value = 0.01_dp
      if (present(volatility_initial_precision)) then
         volatility_precision_value = volatility_initial_precision
      end if
      offset_value = 0.0001_dp
      if (present(volatility_offset)) offset_value = volatility_offset
      omit_deterministics = .false.
      if (present(exclude_deterministics)) omit_deterministics = exclude_deterministics
      variables = data%variables
      states = variables*size(data%x, 1)
      deterministic_states = variables*data%deterministic_columns
      if (data%info /= 0 .or. .not. allocated(data%y) .or. &
         .not. allocated(data%x) .or. .not. allocated(data%sur) .or. &
         variables < 1 .or. states < 1 .or. &
         any(shape(data%sur) /= [variables*data%observations, states]) .or. &
         initial_precision_value < 0.0_dp .or. state_shape_value <= 0.0_dp .or. &
         state_rate_value <= 0.0_dp .or. deterministic_rate_value <= 0.0_dp .or. &
         scale_value <= 0.0_dp .or. measurement_shape_value <= 0.0_dp .or. &
         measurement_rate_value <= 0.0_dp .or. probability < 0.0_dp .or. &
         probability > 1.0_dp .or. volatility_variance_value <= 0.0_dp .or. &
         volatility_shape_value <= 0.0_dp .or. volatility_rate_value <= 0.0_dp .or. &
         volatility_precision_value < 0.0_dp .or. offset_value <= 0.0_dp .or. &
         (trim(out%observation_prior) /= 'wishart' .and. &
         trim(out%observation_prior) /= 'gamma' .and. &
         trim(out%observation_prior) /= 'sv') .or. &
         (trim(out%volatility_method) /= 'ksc' .and. &
         trim(out%volatility_method) /= 'ocsn')) then
         out%info = 1
         return
      end if
      call multivariate_ols_initial(data%y, data%x, coefficients, &
         out%initial_covariance, residual, status)
      if (status /= 0) then
         out%info = 2
         return
      end if
      allocate(out%initial_state(states), out%initial_state_prior_mean(states), &
         out%initial_state_prior_precision(states, states), out%state_shape(states), &
         out%state_rate(states), out%covariance_scale(variables, variables), &
         out%measurement_shape(variables), out%measurement_rate(variables))
      out%initial_state = reshape(coefficients, [states])
      out%initial_state_prior_mean = out%initial_state
      out%initial_state_prior_precision = 0.0_dp
      do item = 1, states
         out%initial_state_prior_precision(item, item) = initial_precision_value
      end do
      out%state_shape = state_shape_value
      out%state_rate = state_rate_value
      if (deterministic_states > 0) then
         out%state_rate(states - deterministic_states + 1:) = deterministic_rate_value
      end if
      out%covariance_scale = 0.0_dp
      do item = 1, variables
         out%covariance_scale(item, item) = scale_value
      end do
      out%covariance_df = real(variables, dp)
      if (present(covariance_df)) out%covariance_df = covariance_df
      out%measurement_shape = measurement_shape_value
      out%measurement_rate = measurement_rate_value
      if (out%covariance_df <= real(variables - 1, dp)) then
         out%info = 3
         return
      end if
      if (present(bvs_probability)) then
         allocate(out%inclusion_probability(states))
         out%inclusion_probability = probability
         selectable_count = states
         if (omit_deterministics) selectable_count = states - deterministic_states
         allocate(out%selectable(selectable_count))
         out%selectable = [(item, item=1, selectable_count)]
         out%use_bvs = .true.
      end if
      if (trim(out%observation_prior) == 'sv') then
         allocate(out%initial_log_variance(variables, data%observations), &
            out%initial_log_variance_level(variables), &
            out%initial_log_variance_state_variance(variables), &
            out%log_variance_state_shape(variables), &
            out%log_variance_state_rate(variables), &
            out%log_variance_initial_prior_mean(variables), &
            out%log_variance_initial_prior_precision(variables, variables), &
            out%log_variance_offset(variables))
         do item = 1, variables
            out%initial_log_variance_level(item) = &
               log(max(out%initial_covariance(item, item), tiny(1.0_dp)))
         end do
         do time = 1, data%observations
            out%initial_log_variance(:, time) = out%initial_log_variance_level
         end do
         out%initial_log_variance_state_variance = volatility_variance_value
         out%log_variance_state_shape = volatility_shape_value
         out%log_variance_state_rate = volatility_rate_value
         out%log_variance_initial_prior_mean = out%initial_log_variance_level
         out%log_variance_initial_prior_precision = 0.0_dp
         do item = 1, variables
            out%log_variance_initial_prior_precision(item, item) = &
               volatility_precision_value
         end do
         out%log_variance_offset = offset_value
      end if
   end function bvartools_prepare_tvp_bvar_prior

   function bvartools_fit_tvp_bvar(data, prior, iterations, burnin) result(out)
      !! Fit a prepared TVP-BVAR using a compatible state and observation bundle.
      type(bvartools_var_data_t), intent(in) :: data !! Prepared BVAR response and SUR design.
      type(bvartools_tvp_bvar_prior_bundle_t), intent(in) :: prior !! Dimensioned TVP-BVAR prior bundle.
      integer, intent(in) :: iterations !! Number of retained Gibbs draws.
      integer, intent(in) :: burnin !! Number of discarded Gibbs draws.
      type(bvartools_tvp_bvar_fit_t) :: out

      if (data%info /= 0 .or. prior%info /= 0 .or. iterations < 1 .or. burnin < 0) then
         out%info = 1
         return
      end if
      select case (trim(prior%observation_prior))
      case ('wishart')
         if (prior%use_bvs) then
            out%draws = bvartools_tvp_bvar_gibbs(data%y, data%sur, &
               prior%initial_state, prior%initial_covariance, iterations, burnin, &
               prior%initial_state_prior_mean, &
               prior%initial_state_prior_precision, prior%state_shape, &
               prior%state_rate, prior%covariance_scale, prior%covariance_df, &
               inclusion_probability=prior%inclusion_probability, &
               selectable=prior%selectable)
         else
            out%draws = bvartools_tvp_bvar_gibbs(data%y, data%sur, &
               prior%initial_state, prior%initial_covariance, iterations, burnin, &
               prior%initial_state_prior_mean, &
               prior%initial_state_prior_precision, prior%state_shape, &
               prior%state_rate, prior%covariance_scale, prior%covariance_df)
         end if
      case ('gamma')
         if (prior%use_bvs) then
            out%draws = bvartools_tvp_bvar_gibbs(data%y, data%sur, &
               prior%initial_state, prior%initial_covariance, iterations, burnin, &
               prior%initial_state_prior_mean, &
               prior%initial_state_prior_precision, prior%state_shape, &
               prior%state_rate, measurement_shape=prior%measurement_shape, &
               measurement_rate=prior%measurement_rate, &
               inclusion_probability=prior%inclusion_probability, &
               selectable=prior%selectable)
         else
            out%draws = bvartools_tvp_bvar_gibbs(data%y, data%sur, &
               prior%initial_state, prior%initial_covariance, iterations, burnin, &
               prior%initial_state_prior_mean, &
               prior%initial_state_prior_precision, prior%state_shape, &
               prior%state_rate, measurement_shape=prior%measurement_shape, &
               measurement_rate=prior%measurement_rate)
         end if
      case ('sv')
         if (prior%use_bvs) then
            out%draws = bvartools_tvp_bvar_gibbs(data%y, data%sur, &
               prior%initial_state, prior%initial_covariance, iterations, burnin, &
               prior%initial_state_prior_mean, &
               prior%initial_state_prior_precision, prior%state_shape, &
               prior%state_rate, inclusion_probability=prior%inclusion_probability, &
               selectable=prior%selectable, &
               initial_log_variance=prior%initial_log_variance, &
               initial_log_variance_level=prior%initial_log_variance_level, &
               initial_log_variance_state_variance= &
               prior%initial_log_variance_state_variance, &
               log_variance_state_shape=prior%log_variance_state_shape, &
               log_variance_state_rate=prior%log_variance_state_rate, &
               log_variance_initial_prior_mean= &
               prior%log_variance_initial_prior_mean, &
               log_variance_initial_prior_precision= &
               prior%log_variance_initial_prior_precision, &
               log_variance_offset=prior%log_variance_offset, &
               volatility_method=prior%volatility_method)
         else
            out%draws = bvartools_tvp_bvar_gibbs(data%y, data%sur, &
               prior%initial_state, prior%initial_covariance, iterations, burnin, &
               prior%initial_state_prior_mean, &
               prior%initial_state_prior_precision, prior%state_shape, &
               prior%state_rate, initial_log_variance=prior%initial_log_variance, &
               initial_log_variance_level=prior%initial_log_variance_level, &
               initial_log_variance_state_variance= &
               prior%initial_log_variance_state_variance, &
               log_variance_state_shape=prior%log_variance_state_shape, &
               log_variance_state_rate=prior%log_variance_state_rate, &
               log_variance_initial_prior_mean= &
               prior%log_variance_initial_prior_mean, &
               log_variance_initial_prior_precision= &
               prior%log_variance_initial_prior_precision, &
               log_variance_offset=prior%log_variance_offset, &
               volatility_method=prior%volatility_method)
         end if
      case default
         out%info = 2
         return
      end select
      if (out%draws%info /= 0) then
         out%info = 10 + out%draws%info
         return
      end if
      out%observations = data%observations
      out%variables = data%variables
      out%lag_order = data%lag_order
   end function bvartools_fit_tvp_bvar

   pure function bvartools_prepare_structural_tvp_bvar_prior(data, &
      time_varying_covariance, initial_state_precision, state_shape, state_rate, &
      structural_state_rate, covariance_initial_precision, &
      covariance_state_shape, covariance_state_rate, measurement_shape, &
      measurement_rate, coefficient_bvs_probability, structural_bvs_probability, &
      covariance_bvs_probability, covariance_ssvs_probability, covariance_tau0, &
      covariance_tau1) result(out)
      !! Construct priors for the identified structural joint TVP-BVAR sampler.
      type(bvartools_var_data_t), intent(in) :: data !! Prepared structural BVAR data.
      logical, intent(in), optional :: time_varying_covariance !! Use random-walk covariance coefficients.
      real(dp), intent(in), optional :: initial_state_precision !! Initial coefficient-state precision.
      real(dp), intent(in), optional :: state_shape !! Coefficient-state precision shape.
      real(dp), intent(in), optional :: state_rate !! Reduced-form coefficient-state precision rate.
      real(dp), intent(in), optional :: structural_state_rate !! Structural coefficient-state precision rate.
      real(dp), intent(in), optional :: covariance_initial_precision !! Initial covariance-state precision.
      real(dp), intent(in), optional :: covariance_state_shape !! Covariance-state precision shape.
      real(dp), intent(in), optional :: covariance_state_rate !! Covariance-state precision rate.
      real(dp), intent(in), optional :: measurement_shape !! Orthogonal innovation-precision shape.
      real(dp), intent(in), optional :: measurement_rate !! Orthogonal innovation-precision rate.
      real(dp), intent(in), optional :: coefficient_bvs_probability !! Reduced-form trajectory BVS probability.
      real(dp), intent(in), optional :: structural_bvs_probability !! Structural trajectory BVS probability.
      real(dp), intent(in), optional :: covariance_bvs_probability !! Covariance trajectory BVS probability.
      real(dp), intent(in), optional :: covariance_ssvs_probability !! Constant covariance SSVS probability.
      real(dp), intent(in), optional :: covariance_tau0 !! Constant covariance SSVS spike deviation.
      real(dp), intent(in), optional :: covariance_tau1 !! Constant covariance SSVS slab deviation.
      type(bvartools_structural_tvp_bvar_prior_t) :: out
      real(dp), allocatable :: coefficients(:, :), residual(:, :), covariance(:, :)
      real(dp) :: initial_precision_value, shape_value, rate_value
      real(dp) :: structural_rate_value, covariance_precision_value
      real(dp) :: covariance_shape_value, covariance_rate_value
      real(dp) :: measurement_shape_value, measurement_rate_value
      real(dp) :: probability, spike, slab
      integer :: variables, observations, reduced_states, structural_states
      integer :: total_states, covariance_count, selected, item, status
      logical :: covariance_is_tvp

      covariance_is_tvp = .true.
      if (present(time_varying_covariance)) covariance_is_tvp = &
         time_varying_covariance
      initial_precision_value = 1.0_dp
      if (present(initial_state_precision)) initial_precision_value = &
         initial_state_precision
      shape_value = 3.0_dp
      if (present(state_shape)) shape_value = state_shape
      rate_value = 0.0001_dp
      if (present(state_rate)) rate_value = state_rate
      structural_rate_value = 0.01_dp
      if (present(structural_state_rate)) structural_rate_value = &
         structural_state_rate
      covariance_precision_value = 1.0_dp
      if (present(covariance_initial_precision)) covariance_precision_value = &
         covariance_initial_precision
      covariance_shape_value = 3.0_dp
      if (present(covariance_state_shape)) covariance_shape_value = &
         covariance_state_shape
      covariance_rate_value = 0.0001_dp
      if (present(covariance_state_rate)) covariance_rate_value = &
         covariance_state_rate
      measurement_shape_value = 3.0_dp
      if (present(measurement_shape)) measurement_shape_value = measurement_shape
      measurement_rate_value = 0.0001_dp
      if (present(measurement_rate)) measurement_rate_value = measurement_rate
      variables = data%variables
      observations = data%observations
      reduced_states = 0
      if (allocated(data%sur)) reduced_states = size(data%sur, 2)
      structural_states = variables*(variables - 1)/2
      total_states = reduced_states + structural_states
      covariance_count = structural_states
      if (data%info /= 0 .or. variables < 2 .or. observations < 1 .or. &
         .not. allocated(data%y) .or. .not. allocated(data%x) .or. &
         .not. allocated(data%sur) .or. .not. allocated(data%structural) .or. &
         any(shape(data%sur) /= [variables*observations, reduced_states]) .or. &
         any(shape(data%structural) /= &
         [variables*observations, structural_states]) .or. &
         initial_precision_value < 0.0_dp .or. shape_value <= 0.0_dp .or. &
         rate_value <= 0.0_dp .or. structural_rate_value <= 0.0_dp .or. &
         covariance_precision_value < 0.0_dp .or. &
         covariance_shape_value <= 0.0_dp .or. covariance_rate_value <= 0.0_dp .or. &
         measurement_shape_value <= 0.0_dp .or. measurement_rate_value <= 0.0_dp) then
         out%info = 1
         return
      end if
      if (present(coefficient_bvs_probability)) then
         if (coefficient_bvs_probability < 0.0_dp .or. &
            coefficient_bvs_probability > 1.0_dp) then
            out%info = 2
            return
         end if
      end if
      if (present(structural_bvs_probability)) then
         if (structural_bvs_probability < 0.0_dp .or. &
            structural_bvs_probability > 1.0_dp) then
            out%info = 3
            return
         end if
      end if
      if (present(covariance_bvs_probability)) then
         if (.not. covariance_is_tvp .or. covariance_bvs_probability < 0.0_dp .or. &
            covariance_bvs_probability > 1.0_dp) then
            out%info = 4
            return
         end if
      end if
      if (present(covariance_ssvs_probability)) then
         if (covariance_is_tvp .or. covariance_ssvs_probability < 0.0_dp .or. &
            covariance_ssvs_probability > 1.0_dp) then
            out%info = 5
            return
         end if
      end if
      call multivariate_ols_initial(data%y, data%x, coefficients, covariance, &
         residual, status)
      if (status /= 0) then
         out%info = 6
         return
      end if
      allocate(out%initial_state(total_states), &
         out%initial_state_prior_mean(total_states), &
         out%initial_state_prior_precision(total_states, total_states), &
         out%state_shape(total_states), out%state_rate(total_states), &
         out%initial_diagonal_variance(variables), &
         out%initial_covariance_state(covariance_count), &
         out%covariance_initial_prior_mean(covariance_count), &
         out%covariance_initial_prior_precision(covariance_count, covariance_count), &
         out%covariance_state_shape(covariance_count), &
         out%covariance_state_rate(covariance_count), &
         out%measurement_shape(variables), out%measurement_rate(variables))
      out%initial_state = 0.0_dp
      out%initial_state(:reduced_states) = reshape(coefficients, [reduced_states])
      out%initial_state_prior_mean = out%initial_state
      out%initial_state_prior_precision = &
         initial_precision_value*identity_matrix(total_states)
      out%state_shape = shape_value
      out%state_rate(:reduced_states) = rate_value
      out%state_rate(reduced_states + 1:) = structural_rate_value
      do item = 1, variables
         out%initial_diagonal_variance(item) = &
            max(covariance(item, item), tiny(1.0_dp))
      end do
      out%initial_covariance_state = 0.0_dp
      out%covariance_initial_prior_mean = 0.0_dp
      out%covariance_initial_prior_precision = &
         covariance_precision_value*identity_matrix(covariance_count)
      out%covariance_state_shape = covariance_shape_value
      out%covariance_state_rate = covariance_rate_value
      out%measurement_shape = measurement_shape_value
      out%measurement_rate = measurement_rate_value
      out%reduced_form_states = reduced_states
      out%structural_states = structural_states
      out%time_varying_covariance = covariance_is_tvp
      selected = 0
      if (present(coefficient_bvs_probability)) selected = selected + reduced_states
      if (present(structural_bvs_probability)) selected = selected + structural_states
      allocate(out%inclusion_probability(total_states), out%selectable(selected))
      out%inclusion_probability = 0.5_dp
      selected = 0
      if (present(coefficient_bvs_probability)) then
         out%inclusion_probability(:reduced_states) = coefficient_bvs_probability
         out%selectable(1:reduced_states) = [(item, item=1, reduced_states)]
         selected = reduced_states
      end if
      if (present(structural_bvs_probability)) then
         out%inclusion_probability(reduced_states + 1:) = &
            structural_bvs_probability
         out%selectable(selected + 1:selected + structural_states) = &
            [(reduced_states + item, item=1, structural_states)]
      end if
      if (covariance_is_tvp) then
         selected = 0
         if (present(covariance_bvs_probability)) selected = covariance_count
         allocate(out%covariance_inclusion_probability(covariance_count), &
            out%covariance_selectable(selected))
         out%covariance_inclusion_probability = 0.5_dp
         if (present(covariance_bvs_probability)) then
            out%covariance_inclusion_probability = covariance_bvs_probability
            out%covariance_selectable = [(item, item=1, covariance_count)]
         end if
      else if (present(covariance_ssvs_probability)) then
         spike = 0.05_dp
         slab = 10.0_dp
         if (present(covariance_tau0)) spike = covariance_tau0
         if (present(covariance_tau1)) slab = covariance_tau1
         if (spike <= 0.0_dp .or. slab <= 0.0_dp) then
            out%info = 7
            return
         end if
         probability = covariance_ssvs_probability
         allocate(out%covariance_tau0(covariance_count), &
            out%covariance_tau1(covariance_count), &
            out%covariance_ssvs_probability(covariance_count), &
            out%covariance_selectable(covariance_count))
         out%covariance_tau0 = spike
         out%covariance_tau1 = slab
         out%covariance_ssvs_probability = probability
         out%covariance_selectable = [(item, item=1, covariance_count)]
         out%use_covariance_ssvs = .true.
      else
         allocate(out%covariance_selectable(0))
      end if
   end function bvartools_prepare_structural_tvp_bvar_prior

   function bvartools_fit_structural_tvp_bvar(data, prior, iterations, burnin) &
      result(out)
      !! Fit a prepared identified structural TVP-BVAR and split its state blocks.
      type(bvartools_var_data_t), intent(in) :: data !! Prepared structural BVAR data.
      type(bvartools_structural_tvp_bvar_prior_t), intent(in) :: prior !! Structural joint prior bundle.
      integer, intent(in) :: iterations !! Number of retained Gibbs draws.
      integer, intent(in) :: burnin !! Number of discarded Gibbs draws.
      type(bvartools_structural_tvp_bvar_fit_t) :: out
      real(dp), allocatable :: design(:, :)
      integer :: rows, reduced, structural, observations, first, last, time, status

      if (data%info /= 0 .or. prior%info /= 0 .or. iterations < 1 .or. burnin < 0 .or. &
         .not. allocated(data%sur) .or. .not. allocated(data%structural)) then
         out%info = 1
         return
      end if
      rows = size(data%sur, 1)
      reduced = prior%reduced_form_states
      structural = prior%structural_states
      observations = data%observations
      if (size(data%sur, 2) /= reduced .or. size(data%structural, 1) /= rows .or. &
         size(data%structural, 2) /= structural) then
         out%info = 2
         return
      end if
      allocate(design(rows, reduced + structural))
      design(:, :reduced) = data%sur
      design(:, reduced + 1:) = data%structural
      if (prior%time_varying_covariance) then
         out%draws = bvartools_joint_tvp_bvar_gibbs(data%y, design, &
            prior%initial_state, prior%initial_diagonal_variance, &
            prior%initial_covariance_state, iterations, burnin, &
            prior%initial_state_prior_mean, prior%initial_state_prior_precision, &
            prior%state_shape, prior%state_rate, &
            prior%covariance_initial_prior_mean, &
            prior%covariance_initial_prior_precision, &
            prior%covariance_state_shape, prior%covariance_state_rate, &
            prior%measurement_shape, prior%measurement_rate, &
            time_varying_covariance=.true., &
            inclusion_probability=prior%inclusion_probability, &
            selectable=prior%selectable, &
            covariance_inclusion_probability= &
            prior%covariance_inclusion_probability, &
            covariance_selectable=prior%covariance_selectable)
      else if (prior%use_covariance_ssvs) then
         out%draws = bvartools_joint_tvp_bvar_gibbs(data%y, design, &
            prior%initial_state, prior%initial_diagonal_variance, &
            prior%initial_covariance_state, iterations, burnin, &
            prior%initial_state_prior_mean, prior%initial_state_prior_precision, &
            prior%state_shape, prior%state_rate, &
            prior%covariance_initial_prior_mean, &
            prior%covariance_initial_prior_precision, &
            prior%covariance_state_shape, prior%covariance_state_rate, &
            prior%measurement_shape, prior%measurement_rate, &
            time_varying_covariance=.false., &
            inclusion_probability=prior%inclusion_probability, &
            selectable=prior%selectable, covariance_tau0=prior%covariance_tau0, &
            covariance_tau1=prior%covariance_tau1, &
            covariance_ssvs_probability=prior%covariance_ssvs_probability, &
            covariance_selectable=prior%covariance_selectable)
      else
         out%draws = bvartools_joint_tvp_bvar_gibbs(data%y, design, &
            prior%initial_state, prior%initial_diagonal_variance, &
            prior%initial_covariance_state, iterations, burnin, &
            prior%initial_state_prior_mean, prior%initial_state_prior_precision, &
            prior%state_shape, prior%state_rate, &
            prior%covariance_initial_prior_mean, &
            prior%covariance_initial_prior_precision, &
            prior%covariance_state_shape, prior%covariance_state_rate, &
            prior%measurement_shape, prior%measurement_rate, &
            time_varying_covariance=.false., &
            inclusion_probability=prior%inclusion_probability, &
            selectable=prior%selectable)
      end if
      if (out%draws%info /= 0) then
         out%info = 10 + out%draws%info
         return
      end if
      allocate(out%reduced_form_states(reduced*observations, iterations), &
         out%structural_states(structural*observations, iterations))
      do time = 1, observations
         first = (time - 1)*(reduced + structural) + 1
         last = first + reduced - 1
         out%reduced_form_states((time - 1)*reduced + 1:time*reduced, :) = &
            out%draws%coefficient_states(first:last, :)
         out%structural_states((time - 1)*structural + 1:time*structural, :) = &
            out%draws%coefficient_states(last + 1: &
            first + reduced + structural - 1, :)
      end do
      call bvartools_structural_impacts(out%structural_states, data%variables, &
         observations, out%structural_impact, status)
      if (status /= 0) then
         out%info = 3
         return
      end if
      if (allocated(out%draws%coefficient_included)) then
         allocate(out%reduced_form_included(reduced, iterations), &
            out%structural_included(structural, iterations))
         out%reduced_form_included = &
            out%draws%coefficient_included(:reduced, :)
         out%structural_included = &
            out%draws%coefficient_included(reduced + 1:, :)
      end if
      out%observations = observations
      out%variables = data%variables
      out%lag_order = data%lag_order
   end function bvartools_fit_structural_tvp_bvar

   pure function bvartools_prepare_tvp_bvec_prior(data, rank, &
      observation_prior, initial_state_precision, state_shape, state_rate, &
      deterministic_state_rate, beta_persistence, beta_precision, &
      covariance_df, covariance_scale, measurement_shape, measurement_rate, &
      bvs_probability, exclude_deterministics, time_varying_covariance, &
      covariance_bvs_probability, covariance_state_precision, &
      covariance_state_shape, covariance_state_rate, volatility_state_variance, &
      volatility_shape, volatility_rate, volatility_initial_precision, &
      volatility_offset, volatility_method) result(out)
      !! Dimension state, cointegration, covariance, and volatility priors for a TVP-BVEC.
      type(bvartools_vecm_data_t), intent(in) :: data !! Prepared BVEC response and regressor matrices.
      integer, intent(in) :: rank !! Positive cointegration rank.
      character(len=*), intent(in), optional :: observation_prior !! Observation mode: wishart, gamma, or sv.
      real(dp), intent(in), optional :: initial_state_precision !! Initial alpha and gamma state precision.
      real(dp), intent(in), optional :: state_shape !! Alpha and gamma innovation precision shape.
      real(dp), intent(in), optional :: state_rate !! Alpha and gamma innovation precision rate.
      real(dp), intent(in), optional :: deterministic_state_rate !! Unrestricted deterministic state rate.
      real(dp), intent(in), optional :: beta_persistence !! Cointegration-state autoregressive persistence.
      real(dp), intent(in), optional :: beta_precision !! Cointegration state and initial precision scale.
      real(dp), intent(in), optional :: covariance_df !! Inverse-Wishart degrees of freedom.
      real(dp), intent(in), optional :: covariance_scale !! Diagonal inverse-Wishart scale.
      real(dp), intent(in), optional :: measurement_shape !! Diagonal observation precision shape.
      real(dp), intent(in), optional :: measurement_rate !! Diagonal observation precision rate.
      real(dp), intent(in), optional :: bvs_probability !! Unrestricted trajectory inclusion probability.
      logical, intent(in), optional :: exclude_deterministics !! Exclude deterministic trajectories from BVS.
      logical, intent(in), optional :: time_varying_covariance !! Enable lower-triangular covariance states.
      real(dp), intent(in), optional :: covariance_bvs_probability !! Covariance-state inclusion probability.
      real(dp), intent(in), optional :: covariance_state_precision !! Initial covariance-state precision.
      real(dp), intent(in), optional :: covariance_state_shape !! Covariance-state innovation precision shape.
      real(dp), intent(in), optional :: covariance_state_rate !! Covariance-state innovation precision rate.
      real(dp), intent(in), optional :: volatility_state_variance !! Initial log-variance innovation variance.
      real(dp), intent(in), optional :: volatility_shape !! Log-variance innovation precision shape.
      real(dp), intent(in), optional :: volatility_rate !! Log-variance innovation precision rate.
      real(dp), intent(in), optional :: volatility_initial_precision !! Initial log-variance prior precision.
      real(dp), intent(in), optional :: volatility_offset !! Positive log-square transformation offset.
      character(len=*), intent(in), optional :: volatility_method !! Mixture method: ksc or ocsn.
      type(bvartools_tvp_bvec_prior_bundle_t) :: out
      real(dp), allocatable :: combined(:, :), coefficients(:, :), residual(:, :)
      real(dp) :: initial_precision_value, state_shape_value, state_rate_value
      real(dp) :: deterministic_rate_value, beta_precision_value, scale_value
      real(dp) :: measurement_shape_value, measurement_rate_value, probability
      real(dp) :: covariance_probability, covariance_precision_value
      real(dp) :: covariance_shape_value, covariance_rate_value
      real(dp) :: volatility_variance_value, volatility_shape_value
      real(dp) :: volatility_rate_value, volatility_precision_value, offset_value
      integer :: variables, w_count, x_count, alpha_count, gamma_count, states
      integer :: beta_count, covariance_count, deterministic_states
      integer :: selectable_count, item, time, status
      logical :: omit_deterministics, use_covariance

      if (present(observation_prior)) out%observation_prior = trim(observation_prior)
      if (present(volatility_method)) out%volatility_method = trim(volatility_method)
      initial_precision_value = 1.0_dp
      if (present(initial_state_precision)) initial_precision_value = initial_state_precision
      state_shape_value = 3.0_dp
      if (present(state_shape)) state_shape_value = state_shape
      state_rate_value = 0.0001_dp
      if (present(state_rate)) state_rate_value = state_rate
      deterministic_rate_value = 0.01_dp
      if (present(deterministic_state_rate)) deterministic_rate_value = deterministic_state_rate
      if (present(beta_persistence)) out%beta_persistence = beta_persistence
      beta_precision_value = 1.0_dp - out%beta_persistence**2
      if (present(beta_precision)) beta_precision_value = beta_precision
      scale_value = 1.0_dp
      if (present(covariance_scale)) scale_value = covariance_scale
      measurement_shape_value = 3.0_dp
      if (present(measurement_shape)) measurement_shape_value = measurement_shape
      measurement_rate_value = 0.0001_dp
      if (present(measurement_rate)) measurement_rate_value = measurement_rate
      probability = 0.5_dp
      if (present(bvs_probability)) probability = bvs_probability
      covariance_probability = 0.5_dp
      if (present(covariance_bvs_probability)) then
         covariance_probability = covariance_bvs_probability
      end if
      covariance_precision_value = 1.0_dp
      if (present(covariance_state_precision)) then
         covariance_precision_value = covariance_state_precision
      end if
      covariance_shape_value = 3.0_dp
      if (present(covariance_state_shape)) covariance_shape_value = covariance_state_shape
      covariance_rate_value = 0.0001_dp
      if (present(covariance_state_rate)) covariance_rate_value = covariance_state_rate
      volatility_variance_value = 0.05_dp
      if (present(volatility_state_variance)) then
         volatility_variance_value = volatility_state_variance
      end if
      volatility_shape_value = 3.0_dp
      if (present(volatility_shape)) volatility_shape_value = volatility_shape
      volatility_rate_value = 0.0001_dp
      if (present(volatility_rate)) volatility_rate_value = volatility_rate
      volatility_precision_value = 0.01_dp
      if (present(volatility_initial_precision)) then
         volatility_precision_value = volatility_initial_precision
      end if
      offset_value = 0.0001_dp
      if (present(volatility_offset)) offset_value = volatility_offset
      omit_deterministics = .false.
      if (present(exclude_deterministics)) omit_deterministics = exclude_deterministics
      use_covariance = .false.
      if (present(time_varying_covariance)) use_covariance = time_varying_covariance
      variables = data%variables
      w_count = size(data%w, 1)
      x_count = size(data%x, 1)
      alpha_count = variables*rank
      gamma_count = variables*x_count
      states = alpha_count + gamma_count
      beta_count = w_count*rank
      covariance_count = variables*(variables - 1)/2
      deterministic_states = variables*data%unrestricted_deterministic_columns
      if (data%info /= 0 .or. .not. allocated(data%y) .or. &
         .not. allocated(data%w) .or. .not. allocated(data%x) .or. &
         variables < 1 .or. rank < 1 .or. rank > min(variables, w_count) .or. &
         states < 1 .or. abs(out%beta_persistence) >= 1.0_dp .or. &
         beta_precision_value <= 0.0_dp .or. initial_precision_value < 0.0_dp .or. &
         state_shape_value <= 0.0_dp .or. state_rate_value <= 0.0_dp .or. &
         deterministic_rate_value <= 0.0_dp .or. scale_value <= 0.0_dp .or. &
         measurement_shape_value <= 0.0_dp .or. measurement_rate_value <= 0.0_dp .or. &
         probability < 0.0_dp .or. probability > 1.0_dp .or. &
         covariance_probability < 0.0_dp .or. covariance_probability > 1.0_dp .or. &
         covariance_precision_value < 0.0_dp .or. covariance_shape_value <= 0.0_dp .or. &
         covariance_rate_value <= 0.0_dp .or. volatility_variance_value <= 0.0_dp .or. &
         volatility_shape_value <= 0.0_dp .or. volatility_rate_value <= 0.0_dp .or. &
         volatility_precision_value < 0.0_dp .or. offset_value <= 0.0_dp .or. &
         (trim(out%observation_prior) /= 'wishart' .and. &
         trim(out%observation_prior) /= 'gamma' .and. &
         trim(out%observation_prior) /= 'sv') .or. &
         (trim(out%volatility_method) /= 'ksc' .and. &
         trim(out%volatility_method) /= 'ocsn') .or. &
         (use_covariance .and. trim(out%observation_prior) /= 'gamma') .or. &
         (present(covariance_bvs_probability) .and. .not. use_covariance)) then
         out%info = 1
         return
      end if
      allocate(combined(w_count + x_count, data%observations))
      combined(:w_count, :) = data%w
      if (x_count > 0) combined(w_count + 1:, :) = data%x
      call multivariate_ols_initial(data%y, combined, coefficients, &
         out%initial_covariance, residual, status)
      if (status /= 0) then
         out%info = 2
         return
      end if
      allocate(out%initial_alpha(variables, rank), out%initial_beta(w_count, rank), &
         out%initial_gamma(variables, x_count), &
         out%initial_state_prior_mean(states), &
         out%initial_state_prior_precision(states, states), out%state_shape(states), &
         out%state_rate(states), out%beta_state_precision(beta_count, beta_count), &
         out%beta_initial_prior_mean(beta_count), &
         out%beta_initial_prior_precision(beta_count, beta_count), &
         out%covariance_scale(variables, variables), &
         out%measurement_shape(variables), out%measurement_rate(variables))
      out%initial_beta = 0.0_dp
      do item = 1, rank
         out%initial_beta(item, item) = 1.0_dp
      end do
      out%initial_alpha = coefficients(:, :rank)
      if (x_count > 0) out%initial_gamma = coefficients(:, w_count + 1:)
      out%initial_state_prior_mean(:alpha_count) = reshape(out%initial_alpha, &
         [alpha_count])
      if (gamma_count > 0) then
         out%initial_state_prior_mean(alpha_count + 1:) = &
            reshape(out%initial_gamma, [gamma_count])
      end if
      out%initial_state_prior_precision = 0.0_dp
      do item = 1, states
         out%initial_state_prior_precision(item, item) = initial_precision_value
      end do
      out%state_shape = state_shape_value
      out%state_rate = state_rate_value
      if (deterministic_states > 0) then
         out%state_rate(states - deterministic_states + 1:) = deterministic_rate_value
      end if
      out%beta_state_precision = 0.0_dp
      out%beta_initial_prior_precision = 0.0_dp
      do item = 1, beta_count
         out%beta_state_precision(item, item) = beta_precision_value
         out%beta_initial_prior_precision(item, item) = beta_precision_value
      end do
      out%beta_initial_prior_mean = reshape(out%initial_beta, [beta_count])
      out%covariance_scale = 0.0_dp
      do item = 1, variables
         out%covariance_scale(item, item) = scale_value
      end do
      out%covariance_df = real(variables + rank, dp)
      if (present(covariance_df)) out%covariance_df = covariance_df
      out%measurement_shape = measurement_shape_value + real(rank, dp)
      out%measurement_rate = measurement_rate_value
      out%rank = rank
      if (out%covariance_df <= real(variables - 1, dp)) then
         out%info = 3
         return
      end if
      if (present(bvs_probability)) then
         if (gamma_count < 1) then
            out%info = 4
            return
         end if
         allocate(out%inclusion_probability(gamma_count))
         out%inclusion_probability = probability
         selectable_count = gamma_count
         if (omit_deterministics) selectable_count = gamma_count - deterministic_states
         allocate(out%selectable(selectable_count))
         out%selectable = [(item, item=1, selectable_count)]
         out%use_bvs = .true.
      end if
      if (use_covariance) then
         allocate(out%initial_covariance_state(covariance_count), &
            out%covariance_state_initial_prior_mean(covariance_count), &
            out%covariance_state_initial_prior_precision(covariance_count, covariance_count), &
            out%covariance_state_shape(covariance_count), &
            out%covariance_state_rate(covariance_count))
         out%initial_covariance_state = 0.0_dp
         out%covariance_state_initial_prior_mean = 0.0_dp
         out%covariance_state_initial_prior_precision = 0.0_dp
         do item = 1, covariance_count
            out%covariance_state_initial_prior_precision(item, item) = &
               covariance_precision_value
         end do
         out%covariance_state_shape = covariance_shape_value
         out%covariance_state_rate = covariance_rate_value
         out%use_covariance_state = .true.
         if (present(covariance_bvs_probability)) then
            allocate(out%covariance_inclusion_probability(covariance_count), &
               out%covariance_selectable(covariance_count))
            out%covariance_inclusion_probability = covariance_probability
            out%covariance_selectable = [(item, item=1, covariance_count)]
            out%use_covariance_bvs = .true.
         end if
      end if
      if (trim(out%observation_prior) == 'sv') then
         allocate(out%initial_log_variance(variables, data%observations), &
            out%initial_log_variance_level(variables), &
            out%initial_log_variance_state_variance(variables), &
            out%log_variance_state_shape(variables), &
            out%log_variance_state_rate(variables), &
            out%log_variance_initial_prior_mean(variables), &
            out%log_variance_initial_prior_precision(variables, variables), &
            out%log_variance_offset(variables))
         do item = 1, variables
            out%initial_log_variance_level(item) = &
               log(max(out%initial_covariance(item, item), tiny(1.0_dp)))
         end do
         do time = 1, data%observations
            out%initial_log_variance(:, time) = out%initial_log_variance_level
         end do
         out%initial_log_variance_state_variance = volatility_variance_value
         out%log_variance_state_shape = volatility_shape_value
         out%log_variance_state_rate = volatility_rate_value
         out%log_variance_initial_prior_mean = out%initial_log_variance_level
         out%log_variance_initial_prior_precision = 0.0_dp
         do item = 1, variables
            out%log_variance_initial_prior_precision(item, item) = &
               volatility_precision_value
         end do
         out%log_variance_offset = offset_value
      end if
   end function bvartools_prepare_tvp_bvec_prior

   function bvartools_fit_tvp_bvec(data, prior, iterations, burnin) result(out)
      !! Fit a prepared TVP-BVEC using a compatible state and observation bundle.
      type(bvartools_vecm_data_t), intent(in) :: data !! Prepared BVEC response and regressors.
      type(bvartools_tvp_bvec_prior_bundle_t), intent(in) :: prior !! Dimensioned TVP-BVEC prior bundle.
      integer, intent(in) :: iterations !! Number of retained Gibbs draws.
      integer, intent(in) :: burnin !! Number of discarded Gibbs draws.
      type(bvartools_tvp_bvec_fit_t) :: out

      if (data%info /= 0 .or. prior%info /= 0 .or. prior%rank < 1 .or. &
         iterations < 1 .or. burnin < 0) then
         out%info = 1
         return
      end if
      select case (trim(prior%observation_prior))
      case ('wishart')
         call fit_tvp_bvec_wishart(data, prior, iterations, burnin, &
            out%draws, out%info)
      case ('gamma')
         call fit_tvp_bvec_gamma(data, prior, iterations, burnin, &
            out%draws, out%info)
      case ('sv')
         call fit_tvp_bvec_sv(data, prior, iterations, burnin, &
            out%draws, out%info)
      case default
         out%info = 2
      end select
      if (out%info /= 0) return
      out%observations = data%observations
      out%variables = data%variables
      out%level_lags = data%level_lags
      out%rank = prior%rank
   end function bvartools_fit_tvp_bvec

   function bvartools_fit_structural_tvp_bvec(data, prior, iterations, burnin, &
      structural_initial_precision, structural_state_shape, structural_state_rate, &
      structural_bvs_probability) result(out)
      !! Fit an identified structural TVP-BVEC under Wishart, gamma, or SV errors.
      type(bvartools_vecm_data_t), intent(in) :: data !! Prepared structural BVEC data.
      type(bvartools_tvp_bvec_prior_bundle_t), intent(in) :: prior !! Dimensioned TVP-BVEC prior bundle.
      integer, intent(in) :: iterations !! Number of retained Gibbs draws.
      integer, intent(in) :: burnin !! Number of discarded Gibbs draws.
      real(dp), intent(in), optional :: structural_initial_precision !! Structural initial-state precision.
      real(dp), intent(in), optional :: structural_state_shape !! Structural-state precision shape.
      real(dp), intent(in), optional :: structural_state_rate !! Structural-state precision rate.
      real(dp), intent(in), optional :: structural_bvs_probability !! Structural trajectory BVS probability.
      type(bvartools_tvp_bvec_fit_t) :: out
      real(dp), allocatable :: prior_mean(:), prior_precision(:, :)
      real(dp), allocatable :: shapes(:), rates(:), structural_probability(:)
      integer, allocatable :: structural_selectable(:)
      real(dp), allocatable :: gamma_probability(:), covariance_probability(:)
      integer, allocatable :: gamma_selectable(:), covariance_selectable(:)
      real(dp) :: initial_precision, shape_value, rate_value, probability
      integer :: old_states, structural_count, total_states, item, gamma_count, status

      initial_precision = 1.0_dp
      if (present(structural_initial_precision)) initial_precision = &
         structural_initial_precision
      shape_value = 3.0_dp
      if (present(structural_state_shape)) shape_value = structural_state_shape
      rate_value = 0.01_dp
      if (present(structural_state_rate)) rate_value = structural_state_rate
      probability = 0.5_dp
      if (present(structural_bvs_probability)) probability = &
         structural_bvs_probability
      structural_count = data%variables*(data%variables - 1)/2
      gamma_count = size(prior%initial_gamma)
      old_states = size(prior%initial_state_prior_mean)
      total_states = old_states + structural_count
      if (data%info /= 0 .or. prior%info /= 0 .or. data%variables < 2 .or. &
         .not. allocated(data%structural) .or. .not. allocated(data%x) .or. &
         .not. allocated(prior%initial_gamma) .or. &
         any(shape(data%structural) /= &
         [data%variables*data%observations, structural_count]) .or. &
         initial_precision < 0.0_dp .or. shape_value <= 0.0_dp .or. &
         rate_value <= 0.0_dp .or. probability < 0.0_dp .or. probability > 1.0_dp .or. &
         iterations < 1 .or. burnin < 0) then
         out%info = 1
         return
      end if
      allocate(prior_mean(total_states), prior_precision(total_states, total_states), &
         shapes(total_states), rates(total_states), &
         structural_probability(structural_count))
      prior_mean = 0.0_dp
      prior_mean(:old_states) = prior%initial_state_prior_mean
      prior_precision = 0.0_dp
      prior_precision(:old_states, :old_states) = prior%initial_state_prior_precision
      do item = old_states + 1, total_states
         prior_precision(item, item) = initial_precision
      end do
      shapes(:old_states) = prior%state_shape
      shapes(old_states + 1:) = shape_value
      rates(:old_states) = prior%state_rate
      rates(old_states + 1:) = rate_value
      structural_probability = probability
      if (present(structural_bvs_probability)) then
         allocate(structural_selectable(structural_count))
         structural_selectable = [(item, item=1, structural_count)]
      else
         allocate(structural_selectable(0))
      end if
      allocate(gamma_probability(gamma_count))
      gamma_probability = 0.5_dp
      if (prior%use_bvs) then
         gamma_probability = prior%inclusion_probability
         gamma_selectable = prior%selectable
      else
         allocate(gamma_selectable(0))
      end if
      if (prior%use_covariance_state) then
         allocate(covariance_probability(structural_count))
         covariance_probability = 0.5_dp
         if (prior%use_covariance_bvs) then
            covariance_probability = prior%covariance_inclusion_probability
            covariance_selectable = prior%covariance_selectable
         else
            allocate(covariance_selectable(0))
         end if
      end if
      select case (trim(prior%observation_prior))
      case ('wishart')
         out%draws = bvartools_tvp_bvec_gibbs(data%y, data%w, prior%initial_alpha, &
            prior%initial_beta, prior%initial_covariance, iterations, burnin, &
            prior_mean, prior_precision, shapes, rates, prior%beta_persistence, &
            prior%beta_state_precision, prior%beta_initial_prior_mean, &
            prior%beta_initial_prior_precision, prior%covariance_scale, &
            prior%covariance_df, x=data%x, initial_gamma=prior%initial_gamma, &
            inclusion_probability=gamma_probability, selectable=gamma_selectable, &
            structural_design=data%structural, &
            initial_structural_state=spread(0.0_dp, 1, structural_count), &
            structural_inclusion_probability=structural_probability, &
            structural_selectable=structural_selectable)
      case ('sv')
         out%draws = bvartools_tvp_bvec_gibbs(data%y, data%w, prior%initial_alpha, &
            prior%initial_beta, prior%initial_covariance, iterations, burnin, &
            prior_mean, prior_precision, shapes, rates, prior%beta_persistence, &
            prior%beta_state_precision, prior%beta_initial_prior_mean, &
            prior%beta_initial_prior_precision, x=data%x, &
            initial_gamma=prior%initial_gamma, &
            inclusion_probability=gamma_probability, selectable=gamma_selectable, &
            initial_log_variance=prior%initial_log_variance, &
            initial_log_variance_level=prior%initial_log_variance_level, &
            initial_log_variance_state_variance= &
            prior%initial_log_variance_state_variance, &
            log_variance_state_shape=prior%log_variance_state_shape, &
            log_variance_state_rate=prior%log_variance_state_rate, &
            log_variance_initial_prior_mean= &
            prior%log_variance_initial_prior_mean, &
            log_variance_initial_prior_precision= &
            prior%log_variance_initial_prior_precision, &
            log_variance_offset=prior%log_variance_offset, &
            volatility_method=prior%volatility_method, &
            structural_design=data%structural, &
            initial_structural_state=spread(0.0_dp, 1, structural_count), &
            structural_inclusion_probability=structural_probability, &
            structural_selectable=structural_selectable)
      case ('gamma')
         if (prior%use_covariance_state) then
            out%draws = bvartools_tvp_bvec_gibbs(data%y, data%w, &
               prior%initial_alpha, prior%initial_beta, prior%initial_covariance, &
               iterations, burnin, prior_mean, prior_precision, shapes, rates, &
               prior%beta_persistence, prior%beta_state_precision, &
               prior%beta_initial_prior_mean, prior%beta_initial_prior_precision, &
               measurement_shape=prior%measurement_shape, &
               measurement_rate=prior%measurement_rate, x=data%x, &
               initial_gamma=prior%initial_gamma, &
               inclusion_probability=gamma_probability, &
               selectable=gamma_selectable, &
               initial_covariance_state=prior%initial_covariance_state, &
               covariance_state_initial_prior_mean= &
               prior%covariance_state_initial_prior_mean, &
               covariance_state_initial_prior_precision= &
               prior%covariance_state_initial_prior_precision, &
               covariance_state_shape=prior%covariance_state_shape, &
               covariance_state_rate=prior%covariance_state_rate, &
               covariance_inclusion_probability=covariance_probability, &
               covariance_selectable=covariance_selectable, &
               structural_design=data%structural, &
               initial_structural_state=spread(0.0_dp, 1, structural_count), &
               structural_inclusion_probability=structural_probability, &
               structural_selectable=structural_selectable)
         else
            out%draws = bvartools_tvp_bvec_gibbs(data%y, data%w, &
               prior%initial_alpha, prior%initial_beta, prior%initial_covariance, &
               iterations, burnin, prior_mean, prior_precision, shapes, rates, &
               prior%beta_persistence, prior%beta_state_precision, &
               prior%beta_initial_prior_mean, prior%beta_initial_prior_precision, &
               measurement_shape=prior%measurement_shape, &
               measurement_rate=prior%measurement_rate, x=data%x, &
               initial_gamma=prior%initial_gamma, &
               inclusion_probability=gamma_probability, &
               selectable=gamma_selectable, structural_design=data%structural, &
               initial_structural_state=spread(0.0_dp, 1, structural_count), &
               structural_inclusion_probability=structural_probability, &
               structural_selectable=structural_selectable)
         end if
      case default
         out%info = 2
         return
      end select
      if (out%draws%info /= 0) then
         out%info = 10 + out%draws%info
         return
      end if
      call bvartools_structural_impacts(out%draws%structural, data%variables, &
         data%observations, out%draws%structural_impact, status)
      if (status /= 0) then
         out%info = 3
         return
      end if
      out%observations = data%observations
      out%variables = data%variables
      out%level_lags = data%level_lags
      out%rank = prior%rank
   end function bvartools_fit_structural_tvp_bvec

   subroutine fit_tvp_bvec_wishart(data, prior, iterations, burnin, draws, info)
      !! Dispatch a Wishart TVP-BVEC fit with optional coefficient BVS.
      type(bvartools_vecm_data_t), intent(in) :: data !! Prepared BVEC response and regressors.
      type(bvartools_tvp_bvec_prior_bundle_t), intent(in) :: prior !! TVP-BVEC prior bundle.
      integer, intent(in) :: iterations !! Number of retained Gibbs draws.
      integer, intent(in) :: burnin !! Number of discarded Gibbs draws.
      type(bvartools_tvp_bvec_draws_t), intent(out) :: draws !! Posterior draws.
      integer, intent(out) :: info !! Wrapper status code.

      if (prior%use_covariance_state) then
         info = 3
         return
      end if
      if (prior%use_bvs) then
         draws = bvartools_tvp_bvec_gibbs(data%y, data%w, prior%initial_alpha, &
            prior%initial_beta, prior%initial_covariance, iterations, burnin, &
            prior%initial_state_prior_mean, prior%initial_state_prior_precision, &
            prior%state_shape, prior%state_rate, prior%beta_persistence, &
            prior%beta_state_precision, prior%beta_initial_prior_mean, &
            prior%beta_initial_prior_precision, prior%covariance_scale, &
            prior%covariance_df, x=data%x, initial_gamma=prior%initial_gamma, &
            inclusion_probability=prior%inclusion_probability, &
            selectable=prior%selectable)
      else
         draws = bvartools_tvp_bvec_gibbs(data%y, data%w, prior%initial_alpha, &
            prior%initial_beta, prior%initial_covariance, iterations, burnin, &
            prior%initial_state_prior_mean, prior%initial_state_prior_precision, &
            prior%state_shape, prior%state_rate, prior%beta_persistence, &
            prior%beta_state_precision, prior%beta_initial_prior_mean, &
            prior%beta_initial_prior_precision, prior%covariance_scale, &
            prior%covariance_df, x=data%x, initial_gamma=prior%initial_gamma)
      end if
      info = 0
      if (draws%info /= 0) info = 10 + draws%info
   end subroutine fit_tvp_bvec_wishart

   subroutine fit_tvp_bvec_sv(data, prior, iterations, burnin, draws, info)
      !! Dispatch a stochastic-volatility TVP-BVEC fit with optional coefficient BVS.
      type(bvartools_vecm_data_t), intent(in) :: data !! Prepared BVEC response and regressors.
      type(bvartools_tvp_bvec_prior_bundle_t), intent(in) :: prior !! TVP-BVEC prior bundle.
      integer, intent(in) :: iterations !! Number of retained Gibbs draws.
      integer, intent(in) :: burnin !! Number of discarded Gibbs draws.
      type(bvartools_tvp_bvec_draws_t), intent(out) :: draws !! Posterior draws.
      integer, intent(out) :: info !! Wrapper status code.

      if (prior%use_covariance_state) then
         info = 3
         return
      end if
      if (prior%use_bvs) then
         draws = bvartools_tvp_bvec_gibbs(data%y, data%w, prior%initial_alpha, &
            prior%initial_beta, prior%initial_covariance, iterations, burnin, &
            prior%initial_state_prior_mean, prior%initial_state_prior_precision, &
            prior%state_shape, prior%state_rate, prior%beta_persistence, &
            prior%beta_state_precision, prior%beta_initial_prior_mean, &
            prior%beta_initial_prior_precision, x=data%x, &
            initial_gamma=prior%initial_gamma, &
            inclusion_probability=prior%inclusion_probability, &
            selectable=prior%selectable, &
            initial_log_variance=prior%initial_log_variance, &
            initial_log_variance_level=prior%initial_log_variance_level, &
            initial_log_variance_state_variance= &
            prior%initial_log_variance_state_variance, &
            log_variance_state_shape=prior%log_variance_state_shape, &
            log_variance_state_rate=prior%log_variance_state_rate, &
            log_variance_initial_prior_mean=prior%log_variance_initial_prior_mean, &
            log_variance_initial_prior_precision= &
            prior%log_variance_initial_prior_precision, &
            log_variance_offset=prior%log_variance_offset, &
            volatility_method=prior%volatility_method)
      else
         draws = bvartools_tvp_bvec_gibbs(data%y, data%w, prior%initial_alpha, &
            prior%initial_beta, prior%initial_covariance, iterations, burnin, &
            prior%initial_state_prior_mean, prior%initial_state_prior_precision, &
            prior%state_shape, prior%state_rate, prior%beta_persistence, &
            prior%beta_state_precision, prior%beta_initial_prior_mean, &
            prior%beta_initial_prior_precision, x=data%x, &
            initial_gamma=prior%initial_gamma, &
            initial_log_variance=prior%initial_log_variance, &
            initial_log_variance_level=prior%initial_log_variance_level, &
            initial_log_variance_state_variance= &
            prior%initial_log_variance_state_variance, &
            log_variance_state_shape=prior%log_variance_state_shape, &
            log_variance_state_rate=prior%log_variance_state_rate, &
            log_variance_initial_prior_mean=prior%log_variance_initial_prior_mean, &
            log_variance_initial_prior_precision= &
            prior%log_variance_initial_prior_precision, &
            log_variance_offset=prior%log_variance_offset, &
            volatility_method=prior%volatility_method)
      end if
      info = 0
      if (draws%info /= 0) info = 10 + draws%info
   end subroutine fit_tvp_bvec_sv

   subroutine fit_tvp_bvec_gamma(data, prior, iterations, burnin, draws, info)
      !! Dispatch a diagonal-gamma TVP-BVEC fit with optional state selections.
      type(bvartools_vecm_data_t), intent(in) :: data !! Prepared BVEC response and regressors.
      type(bvartools_tvp_bvec_prior_bundle_t), intent(in) :: prior !! TVP-BVEC prior bundle.
      integer, intent(in) :: iterations !! Number of retained Gibbs draws.
      integer, intent(in) :: burnin !! Number of discarded Gibbs draws.
      type(bvartools_tvp_bvec_draws_t), intent(out) :: draws !! Posterior draws.
      integer, intent(out) :: info !! Wrapper status code.

      if (.not. prior%use_covariance_state) then
         if (prior%use_bvs) then
            draws = bvartools_tvp_bvec_gibbs(data%y, data%w, prior%initial_alpha, &
               prior%initial_beta, prior%initial_covariance, iterations, burnin, &
               prior%initial_state_prior_mean, prior%initial_state_prior_precision, &
               prior%state_shape, prior%state_rate, prior%beta_persistence, &
               prior%beta_state_precision, prior%beta_initial_prior_mean, &
               prior%beta_initial_prior_precision, &
               measurement_shape=prior%measurement_shape, &
               measurement_rate=prior%measurement_rate, x=data%x, &
               initial_gamma=prior%initial_gamma, &
               inclusion_probability=prior%inclusion_probability, &
               selectable=prior%selectable)
         else
            draws = bvartools_tvp_bvec_gibbs(data%y, data%w, prior%initial_alpha, &
               prior%initial_beta, prior%initial_covariance, iterations, burnin, &
               prior%initial_state_prior_mean, prior%initial_state_prior_precision, &
               prior%state_shape, prior%state_rate, prior%beta_persistence, &
               prior%beta_state_precision, prior%beta_initial_prior_mean, &
               prior%beta_initial_prior_precision, &
               measurement_shape=prior%measurement_shape, &
               measurement_rate=prior%measurement_rate, x=data%x, &
               initial_gamma=prior%initial_gamma)
         end if
      else if (prior%use_bvs .and. prior%use_covariance_bvs) then
         draws = bvartools_tvp_bvec_gibbs(data%y, data%w, prior%initial_alpha, &
            prior%initial_beta, prior%initial_covariance, iterations, burnin, &
            prior%initial_state_prior_mean, prior%initial_state_prior_precision, &
            prior%state_shape, prior%state_rate, prior%beta_persistence, &
            prior%beta_state_precision, prior%beta_initial_prior_mean, &
            prior%beta_initial_prior_precision, &
            measurement_shape=prior%measurement_shape, &
            measurement_rate=prior%measurement_rate, x=data%x, &
            initial_gamma=prior%initial_gamma, &
            inclusion_probability=prior%inclusion_probability, &
            selectable=prior%selectable, &
            initial_covariance_state=prior%initial_covariance_state, &
            covariance_state_initial_prior_mean= &
            prior%covariance_state_initial_prior_mean, &
            covariance_state_initial_prior_precision= &
            prior%covariance_state_initial_prior_precision, &
            covariance_state_shape=prior%covariance_state_shape, &
            covariance_state_rate=prior%covariance_state_rate, &
            covariance_inclusion_probability= &
            prior%covariance_inclusion_probability, &
            covariance_selectable=prior%covariance_selectable)
      else if (prior%use_bvs) then
         draws = bvartools_tvp_bvec_gibbs(data%y, data%w, prior%initial_alpha, &
            prior%initial_beta, prior%initial_covariance, iterations, burnin, &
            prior%initial_state_prior_mean, prior%initial_state_prior_precision, &
            prior%state_shape, prior%state_rate, prior%beta_persistence, &
            prior%beta_state_precision, prior%beta_initial_prior_mean, &
            prior%beta_initial_prior_precision, &
            measurement_shape=prior%measurement_shape, &
            measurement_rate=prior%measurement_rate, x=data%x, &
            initial_gamma=prior%initial_gamma, &
            inclusion_probability=prior%inclusion_probability, &
            selectable=prior%selectable, &
            initial_covariance_state=prior%initial_covariance_state, &
            covariance_state_initial_prior_mean= &
            prior%covariance_state_initial_prior_mean, &
            covariance_state_initial_prior_precision= &
            prior%covariance_state_initial_prior_precision, &
            covariance_state_shape=prior%covariance_state_shape, &
            covariance_state_rate=prior%covariance_state_rate)
      else if (prior%use_covariance_bvs) then
         draws = bvartools_tvp_bvec_gibbs(data%y, data%w, prior%initial_alpha, &
            prior%initial_beta, prior%initial_covariance, iterations, burnin, &
            prior%initial_state_prior_mean, prior%initial_state_prior_precision, &
            prior%state_shape, prior%state_rate, prior%beta_persistence, &
            prior%beta_state_precision, prior%beta_initial_prior_mean, &
            prior%beta_initial_prior_precision, &
            measurement_shape=prior%measurement_shape, &
            measurement_rate=prior%measurement_rate, x=data%x, &
            initial_gamma=prior%initial_gamma, &
            initial_covariance_state=prior%initial_covariance_state, &
            covariance_state_initial_prior_mean= &
            prior%covariance_state_initial_prior_mean, &
            covariance_state_initial_prior_precision= &
            prior%covariance_state_initial_prior_precision, &
            covariance_state_shape=prior%covariance_state_shape, &
            covariance_state_rate=prior%covariance_state_rate, &
            covariance_inclusion_probability= &
            prior%covariance_inclusion_probability, &
            covariance_selectable=prior%covariance_selectable)
      else
         draws = bvartools_tvp_bvec_gibbs(data%y, data%w, prior%initial_alpha, &
            prior%initial_beta, prior%initial_covariance, iterations, burnin, &
            prior%initial_state_prior_mean, prior%initial_state_prior_precision, &
            prior%state_shape, prior%state_rate, prior%beta_persistence, &
            prior%beta_state_precision, prior%beta_initial_prior_mean, &
            prior%beta_initial_prior_precision, &
            measurement_shape=prior%measurement_shape, &
            measurement_rate=prior%measurement_rate, x=data%x, &
            initial_gamma=prior%initial_gamma, &
            initial_covariance_state=prior%initial_covariance_state, &
            covariance_state_initial_prior_mean= &
            prior%covariance_state_initial_prior_mean, &
            covariance_state_initial_prior_precision= &
            prior%covariance_state_initial_prior_precision, &
            covariance_state_shape=prior%covariance_state_shape, &
            covariance_state_rate=prior%covariance_state_rate)
      end if
      info = 0
      if (draws%info /= 0) info = 10 + draws%info
   end subroutine fit_tvp_bvec_gamma

   pure function bvartools_normal_posterior(y, x, residual_precision, &
      prior_mean, prior_precision) result(out)
      !! Compute posterior moments for a vectorized multivariate regression.
      real(dp), intent(in) :: y(:, :) !! Endogenous variables by observation.
      real(dp), intent(in) :: x(:, :) !! Regressors by observation.
      real(dp), intent(in) :: residual_precision(:, :) !! Residual precision matrix.
      real(dp), intent(in) :: prior_mean(:) !! Coefficient prior mean.
      real(dp), intent(in) :: prior_precision(:, :) !! Coefficient prior precision.
      type(bvartools_normal_posterior_t) :: out
      real(dp), allocatable :: posterior_precision(:, :), inverse(:, :), score(:)
      integer :: k, m, status

      k = size(y, 1)
      m = size(x, 1)
      if (size(y, 2) /= size(x, 2) .or. any(shape(residual_precision) /= [k, k]) .or. &
         size(prior_mean) /= k*m .or. any(shape(prior_precision) /= [k*m, k*m])) then
         out%info = 1
         return
      end if
      posterior_precision = prior_precision + kronecker_product( &
         matmul(x, transpose(x)), residual_precision)
      call invert_matrix(posterior_precision, inverse, status)
      if (status /= 0) then
         out%info = 2
         return
      end if
      score = matmul(prior_precision, prior_mean) + reshape( &
         matmul(residual_precision, matmul(y, transpose(x))), [k*m])
      out%covariance = inverse
      out%mean = matmul(inverse, score)
   end function bvartools_normal_posterior

   pure function bvartools_normal_draw(posterior, standard_normal) result(draw)
      !! Transform standard normals into a draw from a Gaussian posterior.
      type(bvartools_normal_posterior_t), intent(in) :: posterior !! Gaussian posterior moments.
      real(dp), intent(in) :: standard_normal(:) !! Independent standard-normal draws.
      real(dp), allocatable :: draw(:)
      integer :: status

      allocate(draw(size(posterior%mean)))
      if (posterior%info /= 0 .or. size(standard_normal) /= size(draw)) then
         draw = 0.0_dp
         return
      end if
      call multivariate_normal_from_standard(posterior%mean, posterior%covariance, &
         standard_normal, draw, status)
      if (status /= 0) draw = posterior%mean
   end function bvartools_normal_draw

   pure function bvartools_sur_normal_posterior(y, z, residual_precision, &
      prior_mean, prior_precision) result(out)
      !! Compute Gaussian posterior moments for a SUR design and time-varying precision.
      real(dp), intent(in) :: y(:, :) !! Endogenous variables by observation.
      real(dp), intent(in) :: z(:, :) !! Stacked observation-major SUR design.
      real(dp), intent(in) :: residual_precision(:, :, :) !! Precision matrix by observation.
      real(dp), intent(in) :: prior_mean(:) !! Coefficient prior mean.
      real(dp), intent(in) :: prior_precision(:, :) !! Coefficient prior precision.
      type(bvartools_normal_posterior_t) :: out
      real(dp), allocatable :: posterior_precision(:, :), score(:), block(:, :), inverse(:, :)
      integer :: k, observations, parameters, time, first, last, status

      k = size(y, 1)
      observations = size(y, 2)
      parameters = size(z, 2)
      if (size(z, 1) /= k*observations .or. &
         any(shape(residual_precision) /= [k, k, observations]) .or. &
         size(prior_mean) /= parameters .or. &
         any(shape(prior_precision) /= [parameters, parameters])) then
         out%info = 1
         return
      end if
      posterior_precision = prior_precision
      score = matmul(prior_precision, prior_mean)
      do time = 1, observations
         first = (time - 1)*k + 1
         last = time*k
         block = z(first:last, :)
         posterior_precision = posterior_precision + matmul(transpose(block), &
            matmul(residual_precision(:, :, time), block))
         score = score + matmul(transpose(block), &
            matmul(residual_precision(:, :, time), y(:, time)))
      end do
      call invert_matrix(posterior_precision, inverse, status)
      if (status /= 0) then
         out%info = 2
         return
      end if
      out%covariance = inverse
      out%mean = matmul(inverse, score)
   end function bvartools_sur_normal_posterior

   pure function bvartools_sur_normal_draw(posterior, standard_normal) result(draw)
      !! Draw SUR coefficients from precomputed Gaussian posterior moments.
      type(bvartools_normal_posterior_t), intent(in) :: posterior !! Gaussian posterior moments.
      real(dp), intent(in) :: standard_normal(:) !! Independent standard-normal draws.
      real(dp), allocatable :: draw(:)

      draw = bvartools_normal_draw(posterior, standard_normal)
   end function bvartools_sur_normal_draw

   pure function bvartools_measurement_variance_posterior(errors, shape_prior, &
      rate_prior) result(out)
      !! Compute independent gamma precision posteriors for measurement errors.
      real(dp), intent(in) :: errors(:, :) !! Measurement errors by variable and observation.
      real(dp), intent(in) :: shape_prior(:) !! Prior gamma shapes.
      real(dp), intent(in) :: rate_prior(:) !! Prior gamma rates.
      type(bvartools_gamma_posterior_t) :: out
      integer :: k

      k = size(errors, 1)
      if (size(shape_prior) /= k .or. size(rate_prior) /= k .or. &
         any(shape_prior <= 0.0_dp) .or. any(rate_prior <= 0.0_dp)) then
         out%info = 1
         return
      end if
      out%shape = shape_prior + 0.5_dp*real(size(errors, 2), dp)
      out%rate = rate_prior + 0.5_dp*sum(errors**2, dim=2)
   end function bvartools_measurement_variance_posterior

   pure function bvartools_state_variance_posterior(states, initial_state, &
      shape_prior, rate_prior) result(out)
      !! Compute gamma precision posteriors for random-walk state innovations.
      real(dp), intent(in) :: states(:, :) !! State draws by coefficient and observation.
      real(dp), intent(in) :: initial_state(:) !! State preceding the first observation.
      real(dp), intent(in) :: shape_prior(:) !! Prior gamma shapes.
      real(dp), intent(in) :: rate_prior(:) !! Prior gamma rates.
      type(bvartools_gamma_posterior_t) :: out
      real(dp), allocatable :: innovations(:, :)
      integer :: k

      k = size(states, 1)
      if (size(initial_state) /= k .or. size(shape_prior) /= k .or. &
         size(rate_prior) /= k .or. any(shape_prior <= 0.0_dp) .or. &
         any(rate_prior <= 0.0_dp)) then
         out%info = 1
         return
      end if
      allocate(innovations(k, size(states, 2)))
      innovations(:, 1) = states(:, 1) - initial_state
      if (size(states, 2) > 1) innovations(:, 2:) = states(:, 2:) - states(:, :size(states, 2) - 1)
      out%shape = shape_prior + 0.5_dp*real(size(states, 2), dp)
      out%rate = rate_prior + 0.5_dp*sum(innovations**2, dim=2)
   end function bvartools_state_variance_posterior

   function bvartools_gamma_precision_draw(posterior, inverse) result(draw)
      !! Draw independent precisions or variances from gamma posterior parameters.
      type(bvartools_gamma_posterior_t), intent(in) :: posterior !! Gamma posterior parameters.
      logical, intent(in), optional :: inverse !! Return variances instead of precisions.
      real(dp), allocatable :: draw(:)
      logical :: return_variance
      integer :: item

      return_variance = .false.
      if (present(inverse)) return_variance = inverse
      allocate(draw(size(posterior%shape)))
      if (posterior%info /= 0) then
         draw = 0.0_dp
         return
      end if
      do item = 1, size(draw)
         draw(item) = random_gamma(posterior%shape(item), 1.0_dp/posterior%rate(item))
      end do
      if (return_variance) draw = 1.0_dp/max(draw, tiny(1.0_dp))
   end function bvartools_gamma_precision_draw

   pure function bvartools_covar_prepare_data(y, precision, time_varying) result(out)
      !! Prepare Primiceri lower-triangular covariance regressions.
      real(dp), intent(in) :: y(:, :) !! Reduced-form residuals by variable and observation.
      real(dp), intent(in) :: precision(:, :, :) !! Diagonal innovation precision by observation.
      logical, intent(in), optional :: time_varying !! Create observation-specific coefficient blocks.
      type(bvartools_covariance_data_t) :: out
      integer :: k, observations, covariance_count, columns
      integer :: time, equation, first, last, row
      logical :: tvp

      k = size(y, 1)
      observations = size(y, 2)
      covariance_count = k*(k - 1)/2
      tvp = .false.
      if (present(time_varying)) tvp = time_varying
      if (k < 2 .or. any(shape(precision) /= [k, k, observations])) then
         out%info = 1
         return
      end if
      columns = covariance_count
      if (tvp) columns = covariance_count*observations
      allocate(out%response((k - 1)*observations), &
         out%design((k - 1)*observations, columns), &
         out%precision((k - 1)*observations, (k - 1)*observations))
      out%design = 0.0_dp
      out%precision = 0.0_dp
      do time = 1, observations
         do equation = 2, k
            row = (time - 1)*(k - 1) + equation - 1
            out%response(row) = y(equation, time)
            first = (equation - 1)*(equation - 2)/2 + 1
            last = equation*(equation - 1)/2
            if (tvp) then
               first = first + (time - 1)*covariance_count
               last = last + (time - 1)*covariance_count
            end if
            out%design(row, first:last) = -y(:equation - 1, time)
            out%precision(row, row) = precision(equation, equation, time)
         end do
      end do
   end function bvartools_covar_prepare_data

   pure function bvartools_covar_vector_to_matrix(psi, variables, observations) &
      result(matrix)
      !! Convert constant or time-varying lower-triangular coefficients to blocks.
      real(dp), intent(in) :: psi(:) !! Packed strict-lower-triangular coefficients.
      integer, intent(in) :: variables !! Number of endogenous variables.
      integer, intent(in), optional :: observations !! Number of time-varying blocks.
      real(dp), allocatable :: matrix(:, :)
      integer :: k, count, blocks, block, equation, first, last, offset, index

      k = variables
      count = k*(k - 1)/2
      blocks = 1
      if (present(observations)) blocks = observations
      if (k < 2 .or. size(psi) /= count*blocks) then
         allocate(matrix(0, 0))
         return
      end if
      allocate(matrix(k*blocks, k*blocks))
      matrix = 0.0_dp
      do index = 1, k*blocks
         matrix(index, index) = 1.0_dp
      end do
      do block = 1, blocks
         offset = (block - 1)*k
         do equation = 2, k
            first = (block - 1)*count + (equation - 1)*(equation - 2)/2 + 1
            last = (block - 1)*count + equation*(equation - 1)/2
            matrix(offset + equation, offset + 1:offset + equation - 1) = psi(first:last)
         end do
      end do
   end function bvartools_covar_vector_to_matrix

   pure function bvartools_sur_const_to_tvp(design, variables, observations) result(expanded)
      !! Expand a constant SUR design into observation-specific coefficient blocks.
      real(dp), intent(in) :: design(:, :) !! Observation-major constant SUR design.
      integer, intent(in) :: variables !! Number of endogenous variables.
      integer, intent(in) :: observations !! Number of observations.
      real(dp), allocatable :: expanded(:, :)
      integer :: parameters, time, row_first, row_last, column_first, column_last

      if (variables < 1 .or. observations < 1 .or. &
         size(design, 1) /= variables*observations) then
         allocate(expanded(0, 0))
         return
      end if
      parameters = size(design, 2)
      allocate(expanded(variables*observations, parameters*observations))
      expanded = 0.0_dp
      do time = 1, observations
         row_first = (time - 1)*variables + 1
         row_last = time*variables
         column_first = (time - 1)*parameters + 1
         column_last = time*parameters
         expanded(row_first:row_last, column_first:column_last) = &
            design(row_first:row_last, :)
      end do
   end function bvartools_sur_const_to_tvp

   pure function bvartools_covar_const_posterior(y, precision, prior_mean, &
      prior_precision) result(out)
      !! Compute posterior moments of constant lower-triangular covariance coefficients.
      real(dp), intent(in) :: y(:, :) !! Reduced-form residuals by variable and observation.
      real(dp), intent(in) :: precision(:, :, :) !! Diagonal innovation precision by observation.
      real(dp), intent(in) :: prior_mean(:) !! Covariance-coefficient prior mean.
      real(dp), intent(in) :: prior_precision(:, :) !! Covariance-coefficient prior precision.
      type(bvartools_normal_posterior_t) :: out
      type(bvartools_covariance_data_t) :: prepared
      real(dp), allocatable :: posterior_precision(:, :), inverse(:, :), score(:)
      integer :: status

      prepared = bvartools_covar_prepare_data(y, precision)
      if (prepared%info /= 0 .or. size(prior_mean) /= size(prepared%design, 2) .or. &
         any(shape(prior_precision) /= [size(prior_mean), size(prior_mean)])) then
         out%info = 1
         return
      end if
      posterior_precision = prior_precision + matmul(transpose(prepared%design), &
         matmul(prepared%precision, prepared%design))
      score = matmul(prior_precision, prior_mean) + matmul(transpose(prepared%design), &
         matmul(prepared%precision, prepared%response))
      call invert_matrix(posterior_precision, inverse, status)
      if (status /= 0) then
         out%info = 2
         return
      end if
      out%covariance = inverse
      out%mean = matmul(inverse, score)
   end function bvartools_covar_const_posterior

   pure function bvartools_covar_tvp_posterior(y, precision, state_precision, &
      initial_state, included) result(out)
      !! Compute Chan-Jeliazkov posterior moments of random-walk covariance states.
      real(dp), intent(in) :: y(:, :) !! Reduced-form residuals by variable and observation.
      real(dp), intent(in) :: precision(:, :, :) !! Diagonal innovation precision by observation.
      real(dp), intent(in) :: state_precision(:, :, :) !! State innovation precision by observation.
      real(dp), intent(in) :: initial_state(:) !! Initial covariance coefficient state.
      logical, intent(in), optional :: included(:) !! Included covariance-state trajectories.
      type(bvartools_normal_posterior_t) :: out
      type(bvartools_covariance_data_t) :: prepared
      real(dp), allocatable :: difference(:, :), block_precision(:, :), prior_path(:)
      real(dp), allocatable :: posterior_precision(:, :), inverse(:, :), score(:)
      integer :: count, observations, time, first, last, status, index

      observations = size(y, 2)
      count = size(initial_state)
      prepared = bvartools_covar_prepare_data(y, precision, .true.)
      if (prepared%info /= 0 .or. &
         any(shape(state_precision) /= [count, count, observations]) .or. &
         size(prepared%design, 2) /= count*observations) then
         out%info = 1
         return
      end if
      if (present(included)) then
         if (size(included) /= count) then
            out%info = 1
            return
         end if
         do time = 1, observations
            first = (time - 1)*count + 1
            last = time*count
            do index = 1, count
               if (.not. included(index)) prepared%design(:, first + index - 1) = 0.0_dp
            end do
         end do
      end if
      allocate(difference(count*observations, count*observations), &
         block_precision(count*observations, count*observations), &
         prior_path(count*observations))
      difference = 0.0_dp
      block_precision = 0.0_dp
      do index = 1, count*observations
         difference(index, index) = 1.0_dp
      end do
      do time = 2, observations
         first = (time - 1)*count + 1
         last = time*count
         difference(first:last, first - count:last - count) = -identity_matrix(count)
      end do
      do time = 1, observations
         first = (time - 1)*count + 1
         last = time*count
         block_precision(first:last, first:last) = state_precision(:, :, time)
         prior_path(first:last) = initial_state
      end do
      posterior_precision = matmul(transpose(difference), &
         matmul(block_precision, difference)) + matmul(transpose(prepared%design), &
         matmul(prepared%precision, prepared%design))
      score = matmul(transpose(difference), matmul(block_precision, &
         matmul(difference, prior_path))) + matmul(transpose(prepared%design), &
         matmul(prepared%precision, prepared%response))
      call invert_matrix(posterior_precision, inverse, status)
      if (status /= 0) then
         out%info = 2
         return
      end if
      out%covariance = inverse
      out%mean = matmul(inverse, score)
   end function bvartools_covar_tvp_posterior

   pure function bvartools_bvs(y, z, coefficients, current_inclusion, precision, &
      prior_probability, proposal_uniforms, acceptance_uniforms, selectable) result(out)
      !! Update Korobilis Bayesian variable-selection indicators in supplied order.
      real(dp), intent(in) :: y(:, :) !! Endogenous variables by observation.
      real(dp), intent(in) :: z(:, :) !! Observation-major SUR design.
      real(dp), intent(in) :: coefficients(:) !! Current constant coefficient draw.
      logical, intent(in) :: current_inclusion(:) !! Current inclusion indicators.
      real(dp), intent(in) :: precision(:, :, :) !! Residual precision by observation.
      real(dp), intent(in) :: prior_probability(:) !! Prior inclusion probabilities.
      real(dp), intent(in) :: proposal_uniforms(:) !! Uniforms for prior-gated proposals.
      real(dp), intent(in) :: acceptance_uniforms(:) !! Uniforms for likelihood acceptance.
      integer, intent(in), optional :: selectable(:) !! One-based update positions and order.
      type(bvartools_bvs_result_t) :: out
      real(dp), allocatable :: theta(:), residual0(:), residual1(:)
      real(dp) :: gate, log0, log1
      integer :: k, observations, m, updates, item, position

      k = size(y, 1)
      observations = size(y, 2)
      m = size(coefficients)
      updates = m
      if (present(selectable)) updates = size(selectable)
      if (size(z, 1) /= k*observations .or. size(z, 2) /= m .or. &
         size(current_inclusion) /= m .or. size(prior_probability) /= m .or. &
         any(shape(precision) /= [k, k, observations]) .or. &
         size(proposal_uniforms) /= updates .or. size(acceptance_uniforms) /= updates .or. &
         any(prior_probability < 0.0_dp) .or. any(prior_probability > 1.0_dp)) then
         out%info = 1
         return
      end if
      if (present(selectable)) then
         if (any(selectable < 1) .or. any(selectable > m)) then
            out%info = 2
            return
         end if
      end if
      out%included = current_inclusion
      allocate(out%log_odds(m))
      out%log_odds = 0.0_dp
      do item = 1, updates
         position = item
         if (present(selectable)) position = selectable(item)
         if (out%included(position)) then
            gate = prior_probability(position)
         else
            gate = 1.0_dp - prior_probability(position)
         end if
         if (proposal_uniforms(item) >= gate) cycle
         theta = coefficients
         where (.not. out%included) theta = 0.0_dp
         theta(position) = 0.0_dp
         residual0 = reshape(y, [k*observations]) - matmul(z, theta)
         theta(position) = coefficients(position)
         residual1 = reshape(y, [k*observations]) - matmul(z, theta)
         log0 = -0.5_dp*block_quadratic(residual0, precision) + &
            log(max(1.0_dp - prior_probability(position), tiny(1.0_dp)))
         log1 = -0.5_dp*block_quadratic(residual1, precision) + &
            log(max(prior_probability(position), tiny(1.0_dp)))
         out%log_odds(position) = log1 - log0
         out%included(position) = log1 - log0 >= &
            log(max(acceptance_uniforms(item), tiny(1.0_dp)))
      end do
   end function bvartools_bvs

   pure function bvartools_tvp_bvs(y, design, states, current_inclusion, &
      precision, prior_probability, proposal_uniforms, acceptance_uniforms, &
      selectable) result(out)
      !! Update inclusion indicators shared across coefficient-state trajectories.
      real(dp), intent(in) :: y(:, :) !! Endogenous variables by observation.
      real(dp), intent(in) :: design(:, :) !! Observation-major constant SUR design.
      real(dp), intent(in) :: states(:, :) !! Coefficient states by trajectory and observation.
      logical, intent(in) :: current_inclusion(:) !! Current trajectory inclusion indicators.
      real(dp), intent(in) :: precision(:, :, :) !! Measurement precision by observation.
      real(dp), intent(in) :: prior_probability(:) !! Prior trajectory inclusion probabilities.
      real(dp), intent(in) :: proposal_uniforms(:) !! Uniforms for prior-gated proposals.
      real(dp), intent(in) :: acceptance_uniforms(:) !! Uniforms for likelihood acceptance.
      integer, intent(in), optional :: selectable(:) !! One-based trajectory update order.
      type(bvartools_bvs_result_t) :: out
      real(dp), allocatable :: expanded(:, :), theta0(:, :), theta1(:, :)
      real(dp), allocatable :: residual0(:), residual1(:)
      real(dp) :: gate, log0, log1
      integer :: k, observations, trajectories, updates, item, position

      k = size(y, 1)
      observations = size(y, 2)
      trajectories = size(states, 1)
      updates = trajectories
      if (present(selectable)) updates = size(selectable)
      if (size(states, 2) /= observations .or. size(design, 1) /= k*observations .or. &
         size(design, 2) /= trajectories .or. size(current_inclusion) /= trajectories .or. &
         size(prior_probability) /= trajectories .or. &
         any(shape(precision) /= [k, k, observations]) .or. &
         size(proposal_uniforms) /= updates .or. &
         size(acceptance_uniforms) /= updates .or. &
         any(prior_probability < 0.0_dp) .or. any(prior_probability > 1.0_dp)) then
         out%info = 1
         return
      end if
      if (present(selectable)) then
         if (any(selectable < 1) .or. any(selectable > trajectories)) then
            out%info = 2
            return
         end if
      end if
      expanded = bvartools_sur_const_to_tvp(design, k, observations)
      out%included = current_inclusion
      allocate(out%log_odds(trajectories))
      out%log_odds = 0.0_dp
      do item = 1, updates
         position = item
         if (present(selectable)) position = selectable(item)
         if (out%included(position)) then
            gate = prior_probability(position)
         else
            gate = 1.0_dp - prior_probability(position)
         end if
         if (proposal_uniforms(item) >= gate) cycle
         theta0 = states
         theta1 = states
         where (.not. spread(out%included, 2, observations))
            theta0 = 0.0_dp
            theta1 = 0.0_dp
         end where
         theta0(position, :) = 0.0_dp
         theta1(position, :) = states(position, :)
         residual0 = reshape(y, [k*observations]) - &
            matmul(expanded, reshape(theta0, [trajectories*observations]))
         residual1 = reshape(y, [k*observations]) - &
            matmul(expanded, reshape(theta1, [trajectories*observations]))
         log0 = -0.5_dp*block_quadratic(residual0, precision) + &
            log(max(1.0_dp - prior_probability(position), tiny(1.0_dp)))
         log1 = -0.5_dp*block_quadratic(residual1, precision) + &
            log(max(prior_probability(position), tiny(1.0_dp)))
         out%log_odds(position) = log1 - log0
         out%included(position) = log1 - log0 >= &
            log(max(acceptance_uniforms(item), tiny(1.0_dp)))
      end do
   end function bvartools_tvp_bvs

   pure function bvartools_tvp_covariance_bvs(residuals, states, &
      current_inclusion, precision, prior_probability, proposal_uniforms, &
      acceptance_uniforms, selectable) result(out)
      !! Update inclusion indicators shared across covariance-state trajectories.
      real(dp), intent(in) :: residuals(:, :) !! Reduced-form residuals by observation.
      real(dp), intent(in) :: states(:, :) !! Covariance states by trajectory and observation.
      logical, intent(in) :: current_inclusion(:) !! Current trajectory inclusion indicators.
      real(dp), intent(in) :: precision(:, :, :) !! Orthogonal innovation precision by observation.
      real(dp), intent(in) :: prior_probability(:) !! Prior trajectory inclusion probabilities.
      real(dp), intent(in) :: proposal_uniforms(:) !! Uniforms for prior-gated proposals.
      real(dp), intent(in) :: acceptance_uniforms(:) !! Uniforms for likelihood acceptance.
      integer, intent(in), optional :: selectable(:) !! One-based trajectory update order.
      type(bvartools_bvs_result_t) :: out
      type(bvartools_covariance_data_t) :: prepared
      real(dp), allocatable :: theta0(:, :), theta1(:, :), residual0(:), residual1(:)
      real(dp) :: gate, log0, log1
      integer :: observations, trajectories, updates, item, position

      observations = size(residuals, 2)
      trajectories = size(states, 1)
      updates = trajectories
      if (present(selectable)) updates = size(selectable)
      prepared = bvartools_covar_prepare_data(residuals, precision, .true.)
      if (prepared%info /= 0 .or. size(states, 2) /= observations .or. &
         size(prepared%design, 2) /= trajectories*observations .or. &
         size(current_inclusion) /= trajectories .or. &
         size(prior_probability) /= trajectories .or. &
         size(proposal_uniforms) /= updates .or. &
         size(acceptance_uniforms) /= updates .or. &
         any(prior_probability < 0.0_dp) .or. any(prior_probability > 1.0_dp)) then
         out%info = 1
         return
      end if
      if (present(selectable)) then
         if (any(selectable < 1) .or. any(selectable > trajectories)) then
            out%info = 2
            return
         end if
      end if
      out%included = current_inclusion
      allocate(out%log_odds(trajectories))
      out%log_odds = 0.0_dp
      do item = 1, updates
         position = item
         if (present(selectable)) position = selectable(item)
         if (out%included(position)) then
            gate = prior_probability(position)
         else
            gate = 1.0_dp - prior_probability(position)
         end if
         if (proposal_uniforms(item) >= gate) cycle
         theta0 = states
         theta1 = states
         where (.not. spread(out%included, 2, observations))
            theta0 = 0.0_dp
            theta1 = 0.0_dp
         end where
         theta0(position, :) = 0.0_dp
         theta1(position, :) = states(position, :)
         residual0 = prepared%response - matmul(prepared%design, &
            reshape(theta0, [trajectories*observations]))
         residual1 = prepared%response - matmul(prepared%design, &
            reshape(theta1, [trajectories*observations]))
         log0 = -0.5_dp*dot_product(residual0, matmul(prepared%precision, residual0)) + &
            log(max(1.0_dp - prior_probability(position), tiny(1.0_dp)))
         log1 = -0.5_dp*dot_product(residual1, matmul(prepared%precision, residual1)) + &
            log(max(prior_probability(position), tiny(1.0_dp)))
         out%log_odds(position) = log1 - log0
         out%included(position) = log1 - log0 >= &
            log(max(acceptance_uniforms(item), tiny(1.0_dp)))
      end do
   end function bvartools_tvp_covariance_bvs

   pure function bvartools_stochastic_volatility_draw(residuals, current_log_variance, &
      state_variance, initial_log_variance, offset, mixture_uniforms, &
      standard_normals) result(out)
      !! Draw random-walk log variances using the KSC seven-normal mixture.
      real(dp), intent(in) :: residuals(:, :) !! Innovations by variable and observation.
      real(dp), intent(in) :: current_log_variance(:, :) !! Current log-variance paths.
      real(dp), intent(in) :: state_variance(:) !! Log-variance innovation variances.
      real(dp), intent(in) :: initial_log_variance(:) !! Pre-sample log variances.
      real(dp), intent(in) :: offset(:) !! Positive constants added before taking logarithms.
      real(dp), intent(in) :: mixture_uniforms(:, :) !! Uniforms selecting mixture components.
      real(dp), intent(in) :: standard_normals(:, :) !! Standard normals for path draws.
      type(bvartools_stochastic_volatility_t) :: out
      real(dp), parameter :: probability(7) = [0.00730_dp, 0.10556_dp, &
         0.00002_dp, 0.04395_dp, 0.34001_dp, 0.24566_dp, 0.25750_dp]
      real(dp), parameter :: mixture_mean(7) = [-11.40039_dp, -5.24321_dp, &
         -9.83726_dp, 1.50746_dp, -0.65098_dp, 0.52478_dp, -2.35859_dp]
      real(dp), parameter :: mixture_variance(7) = [5.79596_dp, 2.61369_dp, &
         5.17950_dp, 0.16735_dp, 0.64009_dp, 0.34023_dp, 1.26261_dp]

      out = stochastic_volatility_mixture_draw(residuals, current_log_variance, &
         state_variance, initial_log_variance, offset, mixture_uniforms, &
         standard_normals, probability, mixture_mean, mixture_variance)
   end function bvartools_stochastic_volatility_draw

   function bvartools_stochastic_volatility(residuals, current_log_variance, &
      state_variance, initial_log_variance, offset) result(out)
      !! Draw KSC random-walk log variances using the shared random stream.
      real(dp), intent(in) :: residuals(:, :) !! Innovations by variable and observation.
      real(dp), intent(in) :: current_log_variance(:, :) !! Current log-variance paths.
      real(dp), intent(in) :: state_variance(:) !! Log-variance innovation variances.
      real(dp), intent(in) :: initial_log_variance(:) !! Pre-sample log variances.
      real(dp), intent(in) :: offset(:) !! Positive constants added before taking logarithms.
      type(bvartools_stochastic_volatility_t) :: out
      real(dp), allocatable :: uniforms(:, :), normals(:, :)
      integer :: variable, time

      allocate(uniforms(size(residuals, 1), size(residuals, 2)), &
         normals(size(residuals, 1), size(residuals, 2)))
      do time = 1, size(residuals, 2)
         do variable = 1, size(residuals, 1)
            uniforms(variable, time) = random_uniform()
            normals(variable, time) = random_standard_normal()
         end do
      end do
      out = bvartools_stochastic_volatility_draw(residuals, current_log_variance, &
         state_variance, initial_log_variance, offset, uniforms, normals)
   end function bvartools_stochastic_volatility

   pure function bvartools_stochastic_volatility_ocsn_draw(residuals, &
      current_log_variance, state_variance, initial_log_variance, offset, &
      mixture_uniforms, standard_normals) result(out)
      !! Draw random-walk log variances using the OCSN ten-normal mixture.
      real(dp), intent(in) :: residuals(:, :) !! Innovations by variable and observation.
      real(dp), intent(in) :: current_log_variance(:, :) !! Current log-variance paths.
      real(dp), intent(in) :: state_variance(:) !! Log-variance innovation variances.
      real(dp), intent(in) :: initial_log_variance(:) !! Pre-sample log variances.
      real(dp), intent(in) :: offset(:) !! Positive constants added before taking logarithms.
      real(dp), intent(in) :: mixture_uniforms(:, :) !! Uniforms selecting mixture components.
      real(dp), intent(in) :: standard_normals(:, :) !! Standard normals for path draws.
      type(bvartools_stochastic_volatility_t) :: out
      real(dp), parameter :: probability(10) = [0.00609_dp, 0.04775_dp, &
         0.13057_dp, 0.20674_dp, 0.22715_dp, 0.18842_dp, 0.12047_dp, &
         0.05591_dp, 0.01575_dp, 0.00115_dp]
      real(dp), parameter :: mixture_mean(10) = [1.92677_dp, 1.34744_dp, &
         0.73504_dp, 0.02266_dp, -0.85173_dp, -1.97278_dp, -3.46788_dp, &
         -5.55246_dp, -8.68384_dp, -14.65000_dp]
      real(dp), parameter :: mixture_variance(10) = [0.11265_dp, 0.17788_dp, &
         0.26768_dp, 0.40611_dp, 0.62699_dp, 0.98583_dp, 1.57469_dp, &
         2.54498_dp, 4.16591_dp, 7.33342_dp]

      out = stochastic_volatility_mixture_draw(residuals, current_log_variance, &
         state_variance, initial_log_variance, offset, mixture_uniforms, &
         standard_normals, probability, mixture_mean, mixture_variance)
   end function bvartools_stochastic_volatility_ocsn_draw

   function bvartools_stochastic_volatility_ocsn(residuals, current_log_variance, &
      state_variance, initial_log_variance, offset) result(out)
      !! Draw OCSN random-walk log variances using the shared random stream.
      real(dp), intent(in) :: residuals(:, :) !! Innovations by variable and observation.
      real(dp), intent(in) :: current_log_variance(:, :) !! Current log-variance paths.
      real(dp), intent(in) :: state_variance(:) !! Log-variance innovation variances.
      real(dp), intent(in) :: initial_log_variance(:) !! Pre-sample log variances.
      real(dp), intent(in) :: offset(:) !! Positive constants added before taking logarithms.
      type(bvartools_stochastic_volatility_t) :: out
      real(dp), allocatable :: uniforms(:, :), normals(:, :)
      integer :: variable, time

      allocate(uniforms(size(residuals, 1), size(residuals, 2)), &
         normals(size(residuals, 1), size(residuals, 2)))
      do time = 1, size(residuals, 2)
         do variable = 1, size(residuals, 1)
            uniforms(variable, time) = random_uniform()
            normals(variable, time) = random_standard_normal()
         end do
      end do
      out = bvartools_stochastic_volatility_ocsn_draw(residuals, &
         current_log_variance, state_variance, initial_log_variance, offset, &
         uniforms, normals)
   end function bvartools_stochastic_volatility_ocsn

   pure function stochastic_volatility_mixture_draw(residuals, &
      current_log_variance, state_variance, initial_log_variance, offset, &
      mixture_uniforms, standard_normals, probability, mixture_mean, &
      mixture_variance) result(out)
      !! Draw log-variance paths for a supplied finite normal mixture.
      real(dp), intent(in) :: residuals(:, :) !! Innovations by variable and observation.
      real(dp), intent(in) :: current_log_variance(:, :) !! Current log-variance paths.
      real(dp), intent(in) :: state_variance(:) !! Log-variance innovation variances.
      real(dp), intent(in) :: initial_log_variance(:) !! Pre-sample log variances.
      real(dp), intent(in) :: offset(:) !! Positive constants added before taking logarithms.
      real(dp), intent(in) :: mixture_uniforms(:, :) !! Uniforms selecting mixture components.
      real(dp), intent(in) :: standard_normals(:, :) !! Standard normals for path draws.
      real(dp), intent(in) :: probability(:) !! Mixture component probabilities.
      real(dp), intent(in) :: mixture_mean(:) !! Mixture component means.
      real(dp), intent(in) :: mixture_variance(:) !! Mixture component variances.
      type(bvartools_stochastic_volatility_t) :: out
      real(dp), allocatable :: difference(:, :), random_walk_precision(:, :)
      real(dp), allocatable :: posterior_precision(:, :), posterior_covariance(:, :)
      real(dp), allocatable :: score(:), mean(:), path(:), transformed(:)
      real(dp), allocatable :: log_weight(:), weight(:)
      real(dp) :: maximum, cumulative
      integer :: k, observations, components, variable, time, component, status

      k = size(residuals, 1)
      observations = size(residuals, 2)
      components = size(probability)
      if (k < 1 .or. observations < 1 .or. components < 1 .or. &
         any(shape(current_log_variance) /= [k, observations]) .or. &
         size(state_variance) /= k .or. size(initial_log_variance) /= k .or. &
         size(offset) /= k .or. any(state_variance <= 0.0_dp) .or. &
         any(offset <= 0.0_dp) .or. &
         any(shape(mixture_uniforms) /= [k, observations]) .or. &
         any(shape(standard_normals) /= [k, observations]) .or. &
         size(mixture_mean) /= components .or. &
         size(mixture_variance) /= components .or. any(probability <= 0.0_dp) .or. &
         any(mixture_variance <= 0.0_dp) .or. &
         any(mixture_uniforms < 0.0_dp) .or. any(mixture_uniforms > 1.0_dp)) then
         out%info = 1
         return
      end if
      allocate(out%log_variance(k, observations), out%component(k, observations), &
         difference(observations, observations), log_weight(components), &
         weight(components), path(observations))
      difference = 0.0_dp
      do time = 1, observations
         difference(time, time) = 1.0_dp
         if (time > 1) difference(time, time - 1) = -1.0_dp
      end do
      random_walk_precision = matmul(transpose(difference), difference)
      do variable = 1, k
         transformed = log(residuals(variable, :)**2 + offset(variable))
         do time = 1, observations
            do component = 1, components
               log_weight(component) = log(probability(component)) - &
                  0.5_dp*(log(2.0_dp*acos(-1.0_dp)*mixture_variance(component)) + &
                  (transformed(time) - current_log_variance(variable, time) - &
                  mixture_mean(component))**2/mixture_variance(component))
            end do
            maximum = maxval(log_weight)
            weight = exp(log_weight - maximum)
            weight = weight/sum(weight)
            cumulative = 0.0_dp
            out%component(variable, time) = components
            do component = 1, components
               cumulative = cumulative + weight(component)
               if (mixture_uniforms(variable, time) <= cumulative) then
                  out%component(variable, time) = component
                  exit
               end if
            end do
         end do
         posterior_precision = random_walk_precision/state_variance(variable)
         do time = 1, observations
            component = out%component(variable, time)
            posterior_precision(time, time) = posterior_precision(time, time) + &
               1.0_dp/mixture_variance(component)
         end do
         score = matmul(random_walk_precision/state_variance(variable), &
            spread(initial_log_variance(variable), 1, observations))
         do time = 1, observations
            component = out%component(variable, time)
            score(time) = score(time) + (transformed(time) - &
               mixture_mean(component))/mixture_variance(component)
         end do
         call invert_matrix(posterior_precision, posterior_covariance, status)
         if (status /= 0) then
            out%info = 2
            return
         end if
         mean = matmul(posterior_covariance, score)
         call multivariate_normal_from_standard(mean, posterior_covariance, &
            standard_normals(variable, :), path, status)
         if (status /= 0) then
            out%info = 3
            return
         end if
         out%log_variance(variable, :) = path
      end do
   end function stochastic_volatility_mixture_draw

   pure function bvartools_random_walk_posterior(y, design, measurement_precision, &
      state_precision, initial_state, persistence) result(out)
      !! Compute the joint Gaussian posterior of random-walk regression states.
      real(dp), intent(in) :: y(:, :) !! Endogenous variables by observation.
      real(dp), intent(in) :: design(:, :) !! Constant observation-major SUR design.
      real(dp), intent(in) :: measurement_precision(:, :, :) !! Measurement precision by observation.
      real(dp), intent(in) :: state_precision(:, :, :) !! State innovation precision by observation.
      real(dp), intent(in) :: initial_state(:) !! State before the first observation.
      real(dp), intent(in), optional :: persistence !! Scalar state-transition persistence.
      type(bvartools_normal_posterior_t) :: out
      real(dp), allocatable :: expanded(:, :), difference(:, :), block_precision(:, :)
      real(dp), allocatable :: measurement_blocks(:, :), innovation_mean(:), response(:)
      real(dp), allocatable :: posterior_precision(:, :), inverse(:, :), score(:)
      integer :: k, observations, states, time, first, last, index, status
      real(dp) :: rho

      k = size(y, 1)
      observations = size(y, 2)
      states = size(initial_state)
      rho = 1.0_dp
      if (present(persistence)) rho = persistence
      if (size(design, 1) /= k*observations .or. size(design, 2) /= states .or. &
         any(shape(measurement_precision) /= [k, k, observations]) .or. &
         any(shape(state_precision) /= [states, states, observations])) then
         out%info = 1
         return
      end if
      expanded = bvartools_sur_const_to_tvp(design, k, observations)
      allocate(difference(states*observations, states*observations), &
         block_precision(states*observations, states*observations), &
         measurement_blocks(k*observations, k*observations), &
         innovation_mean(states*observations))
      difference = 0.0_dp
      block_precision = 0.0_dp
      measurement_blocks = 0.0_dp
      innovation_mean = 0.0_dp
      innovation_mean(1:states) = rho*initial_state
      do index = 1, states*observations
         difference(index, index) = 1.0_dp
      end do
      do time = 1, observations
         first = (time - 1)*states + 1
         last = time*states
         block_precision(first:last, first:last) = state_precision(:, :, time)
         if (time > 1) difference(first:last, first - states:last - states) = &
            -rho*identity_matrix(states)
         first = (time - 1)*k + 1
         last = time*k
         measurement_blocks(first:last, first:last) = measurement_precision(:, :, time)
      end do
      response = reshape(y, [k*observations])
      posterior_precision = matmul(transpose(difference), &
         matmul(block_precision, difference)) + matmul(transpose(expanded), &
         matmul(measurement_blocks, expanded))
      score = matmul(transpose(difference), matmul(block_precision, &
         innovation_mean)) + matmul(transpose(expanded), &
         matmul(measurement_blocks, response))
      call invert_matrix(posterior_precision, inverse, status)
      if (status /= 0) then
         out%info = 2
         return
      end if
      out%covariance = inverse
      out%mean = matmul(inverse, score)
   end function bvartools_random_walk_posterior

   pure function bvartools_cointegration_draw(y, beta, w, residual_precision, &
      shrinkage, cointegration_precision, loading_precision, loading_normals, &
      beta_normals, x, gamma_prior_mean, gamma_prior_precision) result(out)
      !! Draw KLS-normalized coefficients for a Bayesian cointegration model.
      real(dp), intent(in) :: y(:, :) !! Differenced endogenous variables by observation.
      real(dp), intent(in) :: beta(:, :) !! Current normalized cointegration matrix.
      real(dp), intent(in) :: w(:, :) !! Regressors entering the cointegration term.
      real(dp), intent(in) :: residual_precision(:, :) !! Innovation precision matrix.
      real(dp), intent(in) :: shrinkage !! Cointegration-space prior shrinkage.
      real(dp), intent(in) :: cointegration_precision(:, :) !! Cointegration-space prior precision.
      real(dp), intent(in) :: loading_precision(:, :) !! Loading-space prior precision.
      real(dp), intent(in) :: loading_normals(:) !! Standard normals for loading and gamma coefficients.
      real(dp), intent(in) :: beta_normals(:) !! Standard normals for unrestricted beta coefficients.
      real(dp), intent(in), optional :: x(:, :) !! Regressors outside the cointegration term.
      real(dp), intent(in), optional :: gamma_prior_mean(:) !! Prior mean of unrestricted coefficients.
      real(dp), intent(in), optional :: gamma_prior_precision(:, :) !! Prior precision of unrestricted coefficients.
      type(bvartools_cointegration_draw_t) :: out
      real(dp), allocatable :: z(:, :), prior_mean(:), prior_precision(:, :)
      real(dp), allocatable :: posterior_precision(:, :), posterior_covariance(:, :)
      real(dp), allocatable :: posterior_mean(:), score(:), coefficients(:)
      real(dp), allocatable :: y_adjusted(:, :), normalized_alpha(:, :)
      real(dp), allocatable :: alpha_root(:, :), alpha_inverse_root(:, :)
      real(dp), allocatable :: beta_precision(:, :), beta_covariance(:, :)
      real(dp), allocatable :: beta_mean(:), beta_vector(:), unrestricted_beta(:, :)
      real(dp), allocatable :: beta_root(:, :), beta_inverse_root(:, :)
      real(dp), allocatable :: beta_cross(:, :), alpha_cross(:, :)
      integer :: k, observations, m, rank, regressors, loading_count
      integer :: gamma_count, beta_count, total_count, status

      k = size(y, 1)
      observations = size(y, 2)
      m = size(w, 1)
      rank = size(beta, 2)
      regressors = 0
      if (present(x)) regressors = size(x, 1)
      loading_count = k*rank
      gamma_count = k*regressors
      beta_count = m*rank
      total_count = loading_count + gamma_count
      if (k < 1 .or. observations < 1 .or. m < 1 .or. rank < 1 .or. &
         size(beta, 1) /= m .or. size(w, 2) /= observations .or. &
         any(shape(residual_precision) /= [k, k]) .or. &
         any(shape(cointegration_precision) /= [m, m]) .or. &
         any(shape(loading_precision) /= [k, k]) .or. &
         size(loading_normals) /= total_count .or. &
         size(beta_normals) /= beta_count .or. shrinkage < 0.0_dp) then
         out%info = 1
         return
      end if
      if (present(x)) then
         if (size(x, 2) /= observations .or. .not. present(gamma_prior_mean) .or. &
            .not. present(gamma_prior_precision)) then
            out%info = 2
            return
         end if
         if (size(gamma_prior_mean) /= gamma_count .or. &
            any(shape(gamma_prior_precision) /= [gamma_count, gamma_count])) then
            out%info = 3
            return
         end if
      else if (present(gamma_prior_mean) .or. present(gamma_prior_precision)) then
         out%info = 4
         return
      end if
      allocate(z(rank + regressors, observations), prior_mean(total_count), &
         prior_precision(total_count, total_count))
      z(1:rank, :) = matmul(transpose(beta), w)
      if (present(x)) z(rank + 1:, :) = x
      prior_mean = 0.0_dp
      prior_precision = 0.0_dp
      beta_cross = matmul(transpose(beta), &
         matmul(cointegration_precision, beta))
      prior_precision(1:loading_count, 1:loading_count) = &
         kronecker_product(shrinkage*beta_cross, loading_precision)
      if (present(x)) then
         prior_mean(loading_count + 1:) = gamma_prior_mean
         prior_precision(loading_count + 1:, loading_count + 1:) = &
            gamma_prior_precision
      end if
      posterior_precision = prior_precision + kronecker_product( &
         matmul(z, transpose(z)), residual_precision)
      call invert_matrix(posterior_precision, posterior_covariance, status)
      if (status /= 0) then
         out%info = 5
         return
      end if
      score = matmul(prior_precision, prior_mean) + reshape(matmul( &
         residual_precision, matmul(y, transpose(z))), [total_count])
      posterior_mean = matmul(posterior_covariance, score)
      allocate(coefficients(total_count))
      call multivariate_normal_from_standard(posterior_mean, posterior_covariance, &
         loading_normals, coefficients, status)
      if (status /= 0) then
         out%info = 6
         return
      end if
      out%alpha = reshape(coefficients(1:loading_count), [k, rank])
      allocate(out%gamma(k, regressors))
      y_adjusted = y
      if (present(x)) then
         out%gamma = reshape(coefficients(loading_count + 1:), [k, regressors])
         y_adjusted = y - matmul(out%gamma, x)
      end if

      alpha_cross = matmul(transpose(out%alpha), out%alpha)
      call symmetric_roots(alpha_cross, alpha_root, alpha_inverse_root, status)
      if (status /= 0) then
         out%info = 7
         return
      end if
      normalized_alpha = matmul(out%alpha, alpha_inverse_root)
      beta_precision = kronecker_product(matmul(transpose(normalized_alpha), &
         matmul(residual_precision, normalized_alpha)), &
         matmul(w, transpose(w))) + kronecker_product( &
         matmul(transpose(normalized_alpha), matmul(loading_precision, &
         normalized_alpha)), shrinkage*cointegration_precision)
      call invert_matrix(beta_precision, beta_covariance, status)
      if (status /= 0) then
         out%info = 8
         return
      end if
      beta_mean = matmul(beta_covariance, reshape(matmul(w, matmul( &
         transpose(y_adjusted), matmul(residual_precision, normalized_alpha))), &
         [beta_count]))
      allocate(beta_vector(beta_count))
      call multivariate_normal_from_standard(beta_mean, beta_covariance, &
         beta_normals, beta_vector, status)
      if (status /= 0) then
         out%info = 9
         return
      end if
      unrestricted_beta = reshape(beta_vector, [m, rank])
      beta_cross = matmul(transpose(unrestricted_beta), unrestricted_beta)
      call symmetric_roots(beta_cross, beta_root, beta_inverse_root, status)
      if (status /= 0) then
         out%info = 10
         return
      end if
      out%alpha = matmul(normalized_alpha, beta_root)
      out%beta = matmul(unrestricted_beta, beta_inverse_root)
      out%pi = matmul(out%alpha, transpose(out%beta))
   end function bvartools_cointegration_draw

   pure function bvartools_cointegration_sur_draw(y, beta, w, &
      residual_precision, shrinkage, cointegration_precision, loading_precision, &
      loading_normals, beta_normals, x, gamma_prior_mean, &
      gamma_prior_precision) result(out)
      !! Draw KLS-normalized cointegration coefficients with time-varying precision.
      real(dp), intent(in) :: y(:, :) !! Differenced endogenous variables by observation.
      real(dp), intent(in) :: beta(:, :) !! Current normalized cointegration matrix.
      real(dp), intent(in) :: w(:, :) !! Regressors entering the cointegration term.
      real(dp), intent(in) :: residual_precision(:, :, :) !! Innovation precision by observation.
      real(dp), intent(in) :: shrinkage !! Cointegration-space prior shrinkage.
      real(dp), intent(in) :: cointegration_precision(:, :) !! Cointegration-space prior precision.
      real(dp), intent(in) :: loading_precision(:, :) !! Loading-space prior precision.
      real(dp), intent(in) :: loading_normals(:) !! Standard normals for loading and gamma coefficients.
      real(dp), intent(in) :: beta_normals(:) !! Standard normals for unrestricted beta coefficients.
      real(dp), intent(in), optional :: x(:, :) !! Regressors outside the cointegration term.
      real(dp), intent(in), optional :: gamma_prior_mean(:) !! Prior mean of unrestricted coefficients.
      real(dp), intent(in), optional :: gamma_prior_precision(:, :) !! Prior precision of unrestricted coefficients.
      type(bvartools_cointegration_draw_t) :: out
      type(bvartools_normal_posterior_t) :: posterior
      real(dp), allocatable :: design(:, :), beta_design(:, :), prior_mean(:)
      real(dp), allocatable :: prior_precision(:, :), coefficients(:)
      real(dp), allocatable :: y_adjusted(:, :), normalized_alpha(:, :)
      real(dp), allocatable :: alpha_root(:, :), alpha_inverse_root(:, :)
      real(dp), allocatable :: beta_root(:, :), beta_inverse_root(:, :)
      real(dp), allocatable :: beta_cross(:, :), alpha_cross(:, :)
      real(dp), allocatable :: beta_vector(:), unrestricted_beta(:, :)
      real(dp), allocatable :: cointegration_regressor(:, :)
      integer :: k, observations, m, rank, regressors, loading_count
      integer :: gamma_count, beta_count, total_count, time, factor, equation
      integer :: regressor, column, first, last, status

      k = size(y, 1)
      observations = size(y, 2)
      m = size(w, 1)
      rank = size(beta, 2)
      regressors = 0
      if (present(x)) regressors = size(x, 1)
      loading_count = k*rank
      gamma_count = k*regressors
      beta_count = m*rank
      total_count = loading_count + gamma_count
      if (k < 1 .or. observations < 1 .or. m < 1 .or. rank < 1 .or. &
         size(beta, 1) /= m .or. size(w, 2) /= observations .or. &
         any(shape(residual_precision) /= [k, k, observations]) .or. &
         any(shape(cointegration_precision) /= [m, m]) .or. &
         any(shape(loading_precision) /= [k, k]) .or. &
         size(loading_normals) /= total_count .or. &
         size(beta_normals) /= beta_count .or. shrinkage < 0.0_dp) then
         out%info = 1
         return
      end if
      if (present(x)) then
         if (size(x, 2) /= observations .or. .not. present(gamma_prior_mean) .or. &
            .not. present(gamma_prior_precision)) then
            out%info = 2
            return
         end if
         if (size(gamma_prior_mean) /= gamma_count .or. &
            any(shape(gamma_prior_precision) /= [gamma_count, gamma_count])) then
            out%info = 3
            return
         end if
      else if (present(gamma_prior_mean) .or. present(gamma_prior_precision)) then
         out%info = 4
         return
      end if

      cointegration_regressor = matmul(transpose(beta), w)
      allocate(design(k*observations, total_count), prior_mean(total_count), &
         prior_precision(total_count, total_count))
      design = 0.0_dp
      do time = 1, observations
         first = (time - 1)*k + 1
         do factor = 1, rank
            do equation = 1, k
               column = (factor - 1)*k + equation
               design(first + equation - 1, column) = &
                  cointegration_regressor(factor, time)
            end do
         end do
         do regressor = 1, regressors
            do equation = 1, k
               column = loading_count + (regressor - 1)*k + equation
               design(first + equation - 1, column) = x(regressor, time)
            end do
         end do
      end do
      prior_mean = 0.0_dp
      prior_precision = 0.0_dp
      beta_cross = matmul(transpose(beta), &
         matmul(cointegration_precision, beta))
      prior_precision(1:loading_count, 1:loading_count) = &
         kronecker_product(shrinkage*beta_cross, loading_precision)
      if (present(x)) then
         prior_mean(loading_count + 1:) = gamma_prior_mean
         prior_precision(loading_count + 1:, loading_count + 1:) = &
            gamma_prior_precision
      end if
      posterior = bvartools_sur_normal_posterior(y, design, residual_precision, &
         prior_mean, prior_precision)
      if (posterior%info /= 0) then
         out%info = 5
         return
      end if
      allocate(coefficients(total_count))
      call multivariate_normal_from_standard(posterior%mean, posterior%covariance, &
         loading_normals, coefficients, status)
      if (status /= 0) then
         out%info = 6
         return
      end if
      out%alpha = reshape(coefficients(1:loading_count), [k, rank])
      allocate(out%gamma(k, regressors))
      y_adjusted = y
      if (present(x)) then
         out%gamma = reshape(coefficients(loading_count + 1:), [k, regressors])
         y_adjusted = y - matmul(out%gamma, x)
      end if

      alpha_cross = matmul(transpose(out%alpha), out%alpha)
      call symmetric_roots(alpha_cross, alpha_root, alpha_inverse_root, status)
      if (status /= 0) then
         out%info = 7
         return
      end if
      normalized_alpha = matmul(out%alpha, alpha_inverse_root)
      deallocate(prior_mean, prior_precision)
      allocate(beta_design(k*observations, beta_count), prior_mean(beta_count), &
         prior_precision(beta_count, beta_count))
      beta_design = 0.0_dp
      do time = 1, observations
         first = (time - 1)*k + 1
         last = time*k
         do factor = 1, rank
            do regressor = 1, m
               column = (factor - 1)*m + regressor
               beta_design(first:last, column) = &
                  normalized_alpha(:, factor)*w(regressor, time)
            end do
         end do
      end do
      prior_mean = 0.0_dp
      prior_precision = kronecker_product( &
         matmul(transpose(normalized_alpha), matmul(loading_precision, &
         normalized_alpha)), shrinkage*cointegration_precision)
      posterior = bvartools_sur_normal_posterior(y_adjusted, beta_design, &
         residual_precision, prior_mean, prior_precision)
      if (posterior%info /= 0) then
         out%info = 8
         return
      end if
      allocate(beta_vector(beta_count))
      call multivariate_normal_from_standard(posterior%mean, posterior%covariance, &
         beta_normals, beta_vector, status)
      if (status /= 0) then
         out%info = 9
         return
      end if
      unrestricted_beta = reshape(beta_vector, [m, rank])
      beta_cross = matmul(transpose(unrestricted_beta), unrestricted_beta)
      call symmetric_roots(beta_cross, beta_root, beta_inverse_root, status)
      if (status /= 0) then
         out%info = 10
         return
      end if
      out%alpha = matmul(normalized_alpha, beta_root)
      out%beta = matmul(unrestricted_beta, beta_inverse_root)
      out%pi = matmul(out%alpha, transpose(out%beta))
   end function bvartools_cointegration_sur_draw

   function bvartools_bvar_gibbs(y, x, prior_mean, prior_precision, &
      initial_covariance, iterations, burnin, covariance_prior_scale, &
      covariance_prior_df, gamma_shape, gamma_rate, tau0, tau1, &
      inclusion_probability, selectable) result(out)
      !! Run the bvartools constant-parameter BVAR Gibbs sampler.
      real(dp), intent(in) :: y(:, :) !! Endogenous variables by observation.
      real(dp), intent(in) :: x(:, :) !! Regressors by observation.
      real(dp), intent(in) :: prior_mean(:) !! Coefficient prior mean.
      real(dp), intent(in) :: prior_precision(:, :) !! Coefficient prior precision.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial residual covariance.
      integer, intent(in) :: iterations !! Number of retained Gibbs draws.
      integer, intent(in) :: burnin !! Number of discarded Gibbs draws.
      real(dp), intent(in), optional :: covariance_prior_scale(:, :) !! Inverse-Wishart scale matrix.
      real(dp), intent(in), optional :: covariance_prior_df !! Inverse-Wishart prior degrees of freedom.
      real(dp), intent(in), optional :: gamma_shape(:) !! Diagonal precision prior shapes.
      real(dp), intent(in), optional :: gamma_rate(:) !! Diagonal precision prior rates.
      real(dp), intent(in), optional :: tau0(:) !! SSVS excluded-state standard deviations.
      real(dp), intent(in), optional :: tau1(:) !! SSVS included-state standard deviations.
      real(dp), intent(in), optional :: inclusion_probability(:) !! SSVS inclusion probabilities.
      integer, intent(in), optional :: selectable(:) !! SSVS coefficient positions.
      type(bvartools_bvar_draws_t) :: out
      type(bvartools_normal_posterior_t) :: posterior
      type(bvartools_ssvs_result_t) :: selection
      real(dp), allocatable :: covariance(:, :), precision(:, :), coefficient(:)
      real(dp), allocatable :: current_precision(:, :), normal(:), uniforms(:)
      real(dp), allocatable :: residual(:, :), posterior_scale(:, :), diagonal_precision(:)
      real(dp) :: posterior_df
      integer :: k, regressors, coefficient_count, total, draw, retained, status, item
      logical :: full_covariance, use_ssvs

      k = size(y, 1)
      regressors = size(x, 1)
      coefficient_count = k*regressors
      total = iterations + burnin
      full_covariance = present(covariance_prior_scale)
      use_ssvs = present(tau0)
      if (size(y, 2) /= size(x, 2) .or. size(prior_mean) /= coefficient_count .or. &
         any(shape(prior_precision) /= [coefficient_count, coefficient_count]) .or. &
         any(shape(initial_covariance) /= [k, k]) .or. iterations < 1 .or. burnin < 0) then
         out%info = 1
         return
      end if
      if (full_covariance) then
         if (.not. present(covariance_prior_df)) then
            out%info = 2
            return
         end if
         if (any(shape(covariance_prior_scale) /= [k, k]) .or. covariance_prior_df <= &
            real(k - 1, dp)) then
            out%info = 3
            return
         end if
      else
         if (.not. present(gamma_shape) .or. .not. present(gamma_rate)) then
            out%info = 4
            return
         end if
         if (size(gamma_shape) /= k .or. size(gamma_rate) /= k .or. &
            any(gamma_shape <= 0.0_dp) .or. any(gamma_rate <= 0.0_dp)) then
            out%info = 5
            return
         end if
      end if
      if (use_ssvs) then
         if (.not. present(tau1) .or. .not. present(inclusion_probability)) then
            out%info = 6
            return
         end if
         if (size(tau0) /= coefficient_count .or. size(tau1) /= coefficient_count .or. &
            size(inclusion_probability) /= coefficient_count) then
            out%info = 7
            return
         end if
      end if
      covariance = initial_covariance
      call invert_matrix(covariance, precision, status)
      if (status /= 0) then
         out%info = 8
         return
      end if
      current_precision = prior_precision
      coefficient = prior_mean
      allocate(out%coefficients(coefficient_count, iterations), &
         out%covariance(k, k, iterations), normal(coefficient_count))
      if (use_ssvs) then
         allocate(out%included(coefficient_count, iterations))
         if (present(selectable)) then
            allocate(uniforms(size(selectable)))
         else
            allocate(uniforms(coefficient_count))
         end if
      end if
      retained = 0
      do draw = 1, total
         posterior = bvartools_normal_posterior(y, x, precision, prior_mean, current_precision)
         if (posterior%info /= 0) then
            out%info = 9
            return
         end if
         do item = 1, coefficient_count
            normal(item) = random_standard_normal()
         end do
         coefficient = bvartools_normal_draw(posterior, normal)
         if (use_ssvs) then
            do item = 1, size(uniforms)
               uniforms(item) = random_uniform()
            end do
            selection = bvartools_ssvs(coefficient, tau0, tau1, &
               inclusion_probability, uniforms, selectable)
            if (selection%info /= 0) then
               out%info = 10
               return
            end if
            current_precision = 0.0_dp
            do item = 1, coefficient_count
               current_precision(item, item) = selection%precision(item)
            end do
         end if
         residual = y - reshape(matmul(reshape(coefficient, [k, regressors]), x), &
            shape(y))
         if (full_covariance) then
            posterior_scale = covariance_prior_scale + matmul(residual, transpose(residual))
            posterior_df = covariance_prior_df + real(size(y, 2), dp)
            call inverse_wishart_draw(posterior_scale, posterior_df, covariance, precision, status)
            if (status /= 0) then
               out%info = 11
               return
            end if
         else
            allocate(diagonal_precision(k))
            do item = 1, k
               diagonal_precision(item) = random_gamma(gamma_shape(item) + &
                  0.5_dp*real(size(y, 2), dp), 1.0_dp/(gamma_rate(item) + &
                  0.5_dp*sum(residual(item, :)**2)))
            end do
            covariance = 0.0_dp
            precision = 0.0_dp
            do item = 1, k
               precision(item, item) = diagonal_precision(item)
               covariance(item, item) = 1.0_dp/diagonal_precision(item)
            end do
            deallocate(diagonal_precision)
         end if
         if (draw > burnin) then
            retained = retained + 1
            out%coefficients(:, retained) = coefficient
            out%covariance(:, :, retained) = covariance
            if (use_ssvs) out%included(:, retained) = selection%included
         end if
      end do
      out%retained_draws = retained
   end function bvartools_bvar_gibbs

   function bvartools_structural_bvar_gibbs(y, reduced_design, &
      structural_design, reduced_prior_mean, reduced_prior_precision, &
      structural_prior_mean, structural_prior_precision, initial_covariance, &
      iterations, burnin, covariance_prior_scale, covariance_prior_df, &
      gamma_shape, gamma_rate, reduced_tau0, reduced_tau1, &
      reduced_inclusion_probability, reduced_selectable, structural_tau0, &
      structural_tau1, structural_inclusion_probability, &
      structural_selectable) result(out)
      !! Sample a constant identified structural BVAR in observation-major SUR form.
      real(dp), intent(in) :: y(:, :) !! Endogenous variables by observation.
      real(dp), intent(in) :: reduced_design(:, :) !! Reduced-form SUR design.
      real(dp), intent(in) :: structural_design(:, :) !! Contemporaneous structural SUR design.
      real(dp), intent(in) :: reduced_prior_mean(:) !! Reduced-form coefficient prior mean.
      real(dp), intent(in) :: reduced_prior_precision(:, :) !! Reduced-form coefficient prior precision.
      real(dp), intent(in) :: structural_prior_mean(:) !! Structural coefficient prior mean.
      real(dp), intent(in) :: structural_prior_precision(:, :) !! Structural coefficient prior precision.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial innovation covariance.
      integer, intent(in) :: iterations !! Number of retained Gibbs draws.
      integer, intent(in) :: burnin !! Number of discarded Gibbs draws.
      real(dp), intent(in), optional :: covariance_prior_scale(:, :) !! Inverse-Wishart scale.
      real(dp), intent(in), optional :: covariance_prior_df !! Inverse-Wishart degrees of freedom.
      real(dp), intent(in), optional :: gamma_shape(:) !! Diagonal precision prior shapes.
      real(dp), intent(in), optional :: gamma_rate(:) !! Diagonal precision prior rates.
      real(dp), intent(in), optional :: reduced_tau0(:) !! Reduced-form SSVS spike deviations.
      real(dp), intent(in), optional :: reduced_tau1(:) !! Reduced-form SSVS slab deviations.
      real(dp), intent(in), optional :: reduced_inclusion_probability(:) !! Reduced-form SSVS probabilities.
      integer, intent(in), optional :: reduced_selectable(:) !! Selectable reduced-form coefficients.
      real(dp), intent(in), optional :: structural_tau0(:) !! Structural SSVS spike deviations.
      real(dp), intent(in), optional :: structural_tau1(:) !! Structural SSVS slab deviations.
      real(dp), intent(in), optional :: structural_inclusion_probability(:) !! Structural SSVS probabilities.
      integer, intent(in), optional :: structural_selectable(:) !! Selectable structural coefficients.
      type(bvartools_bvar_draws_t) :: out
      type(bvartools_normal_posterior_t) :: posterior
      type(bvartools_ssvs_result_t) :: reduced_selection, structural_selection
      real(dp), allocatable :: design(:, :), prior_mean(:), prior_precision(:, :)
      real(dp), allocatable :: current_precision(:, :), coefficient(:), normal(:)
      real(dp), allocatable :: covariance(:, :), precision(:, :), residual(:, :)
      real(dp), allocatable :: residual_vector(:), posterior_scale(:, :)
      real(dp), allocatable :: diagonal_precision(:), reduced_uniforms(:)
      real(dp), allocatable :: structural_uniforms(:), contemporaneous(:, :), impact(:, :)
      real(dp) :: posterior_df
      integer :: k, observations, reduced_count, structural_count, total_count
      integer :: total, draw, retained, item, status, equation, source, position
      logical :: full_covariance, use_reduced_ssvs, use_structural_ssvs

      k = size(y, 1)
      observations = size(y, 2)
      reduced_count = size(reduced_design, 2)
      structural_count = size(structural_design, 2)
      total_count = reduced_count + structural_count
      total = iterations + burnin
      full_covariance = present(covariance_prior_scale)
      use_reduced_ssvs = present(reduced_tau0)
      use_structural_ssvs = present(structural_tau0)
      if (k < 2 .or. observations < 1 .or. structural_count /= k*(k - 1)/2 .or. &
         size(reduced_design, 1) /= k*observations .or. &
         size(structural_design, 1) /= k*observations .or. &
         size(reduced_prior_mean) /= reduced_count .or. &
         any(shape(reduced_prior_precision) /= [reduced_count, reduced_count]) .or. &
         size(structural_prior_mean) /= structural_count .or. &
         any(shape(structural_prior_precision) /= &
         [structural_count, structural_count]) .or. &
         any(shape(initial_covariance) /= [k, k]) .or. iterations < 1 .or. burnin < 0) then
         out%info = 1
         return
      end if
      if (full_covariance) then
         if (.not. present(covariance_prior_df) .or. &
            any(shape(covariance_prior_scale) /= [k, k]) .or. &
            covariance_prior_df <= real(k - 1, dp)) then
            out%info = 2
            return
         end if
      else if (.not. present(gamma_shape) .or. .not. present(gamma_rate)) then
         out%info = 3
         return
      else if (size(gamma_shape) /= k .or. size(gamma_rate) /= k .or. &
         any(gamma_shape <= 0.0_dp) .or. any(gamma_rate <= 0.0_dp)) then
         out%info = 4
         return
      end if
      if (use_reduced_ssvs) then
         if (.not. present(reduced_tau1) .or. &
            .not. present(reduced_inclusion_probability) .or. &
            size(reduced_tau0) /= reduced_count .or. &
            size(reduced_tau1) /= reduced_count .or. &
            size(reduced_inclusion_probability) /= reduced_count) then
            out%info = 5
            return
         end if
      end if
      if (use_structural_ssvs) then
         if (.not. present(structural_tau1) .or. &
            .not. present(structural_inclusion_probability) .or. &
            size(structural_tau0) /= structural_count .or. &
            size(structural_tau1) /= structural_count .or. &
            size(structural_inclusion_probability) /= structural_count) then
            out%info = 6
            return
         end if
      end if
      allocate(design(k*observations, total_count), prior_mean(total_count), &
         prior_precision(total_count, total_count), current_precision(total_count, total_count))
      design(:, :reduced_count) = reduced_design
      design(:, reduced_count + 1:) = structural_design
      prior_mean = [reduced_prior_mean, structural_prior_mean]
      prior_precision = 0.0_dp
      prior_precision(:reduced_count, :reduced_count) = reduced_prior_precision
      prior_precision(reduced_count + 1:, reduced_count + 1:) = &
         structural_prior_precision
      current_precision = prior_precision
      covariance = initial_covariance
      call invert_matrix(covariance, precision, status)
      if (status /= 0) then
         out%info = 7
         return
      end if
      allocate(normal(total_count), out%coefficients(reduced_count, iterations), &
         out%structural(structural_count, iterations), &
         out%structural_impact(k, k, iterations), out%covariance(k, k, iterations))
      if (use_reduced_ssvs) then
         allocate(out%included(reduced_count, iterations))
         if (present(reduced_selectable)) then
            allocate(reduced_uniforms(size(reduced_selectable)))
         else
            allocate(reduced_uniforms(reduced_count))
         end if
      end if
      if (use_structural_ssvs) then
         allocate(out%structural_included(structural_count, iterations))
         if (present(structural_selectable)) then
            allocate(structural_uniforms(size(structural_selectable)))
         else
            allocate(structural_uniforms(structural_count))
         end if
      end if
      retained = 0
      do draw = 1, total
         posterior = bvartools_sur_normal_posterior(y, design, &
            spread(precision, 3, observations), prior_mean, current_precision)
         if (posterior%info /= 0) then
            out%info = 8
            return
         end if
         do item = 1, total_count
            normal(item) = random_standard_normal()
         end do
         coefficient = bvartools_normal_draw(posterior, normal)
         if (use_reduced_ssvs) then
            do item = 1, size(reduced_uniforms)
               reduced_uniforms(item) = random_uniform()
            end do
            reduced_selection = bvartools_ssvs(coefficient(:reduced_count), &
               reduced_tau0, reduced_tau1, reduced_inclusion_probability, &
               reduced_uniforms, reduced_selectable)
            do item = 1, reduced_count
               current_precision(item, item) = reduced_selection%precision(item)
            end do
         end if
         if (use_structural_ssvs) then
            do item = 1, size(structural_uniforms)
               structural_uniforms(item) = random_uniform()
            end do
            structural_selection = bvartools_ssvs(coefficient(reduced_count + 1:), &
               structural_tau0, structural_tau1, structural_inclusion_probability, &
               structural_uniforms, structural_selectable)
            do item = 1, structural_count
               current_precision(reduced_count + item, reduced_count + item) = &
                  structural_selection%precision(item)
            end do
         end if
         residual_vector = reshape(y, [k*observations]) - matmul(design, coefficient)
         residual = reshape(residual_vector, [k, observations])
         if (full_covariance) then
            posterior_scale = covariance_prior_scale + matmul(residual, transpose(residual))
            posterior_df = covariance_prior_df + real(observations, dp)
            call inverse_wishart_draw(posterior_scale, posterior_df, covariance, &
               precision, status)
         else
            allocate(diagonal_precision(k))
            do item = 1, k
               diagonal_precision(item) = random_gamma(gamma_shape(item) + &
                  0.5_dp*real(observations, dp), 1.0_dp/(gamma_rate(item) + &
                  0.5_dp*sum(residual(item, :)**2)))
            end do
            precision = diagonal_matrix(diagonal_precision)
            covariance = diagonal_matrix(1.0_dp/diagonal_precision)
            deallocate(diagonal_precision)
            status = 0
         end if
         if (status /= 0) then
            out%info = 9
            return
         end if
         if (draw > burnin) then
            retained = retained + 1
            out%coefficients(:, retained) = coefficient(:reduced_count)
            out%structural(:, retained) = coefficient(reduced_count + 1:)
            out%covariance(:, :, retained) = covariance
            contemporaneous = identity_matrix(k)
            position = 0
            do equation = 2, k
               do source = 1, equation - 1
                  position = position + 1
                  contemporaneous(equation, source) = &
                     coefficient(reduced_count + position)
               end do
            end do
            call invert_matrix(contemporaneous, impact, status)
            if (status /= 0) then
               out%info = 10
               return
            end if
            out%structural_impact(:, :, retained) = impact
            if (use_reduced_ssvs) out%included(:, retained) = &
               reduced_selection%included
            if (use_structural_ssvs) out%structural_included(:, retained) = &
               structural_selection%included
         end if
      end do
      out%retained_draws = retained
   end function bvartools_structural_bvar_gibbs

   pure function bvartools_bvar_predictive_from_random(draws, history, lags, &
      standard_normals, future_regressors, structural_impact, probability) result(out)
      !! Simulate constant-parameter BVAR forecasts from supplied normal variates.
      type(bvartools_bvar_draws_t), intent(in) :: draws !! Retained coefficient and covariance draws.
      real(dp), intent(in) :: history(:, :) !! Endogenous history, ordered from oldest to newest observation.
      integer, intent(in) :: lags !! Number of endogenous lags in each coefficient draw.
      real(dp), intent(in) :: standard_normals(:, :, :) !! Independent normals by variable, horizon, and draw.
      real(dp), intent(in), optional :: future_regressors(:, :) !! Future exogenous and deterministic regressors.
      real(dp), intent(in), optional :: structural_impact(:, :, :) !! Draw-specific inverse contemporaneous matrices.
      real(dp), intent(in), optional :: probability !! Equal-tail posterior interval probability.
      type(bvartools_predictive_t) :: out

      if (.not. allocated(draws%coefficients) .or. &
         .not. allocated(draws%covariance)) then
         out%info = 1
         return
      end if
      out = predictive_kernel(draws%coefficients, draws%covariance, history, lags, &
         standard_normals, future_regressors, structural_impact, probability)
      if (draws%info /= 0 .and. out%info == 0) out%info = 10
   end function bvartools_bvar_predictive_from_random

   function bvartools_bvar_predictive(draws, history, lags, horizon, &
      future_regressors, structural_impact, probability) result(out)
      !! Simulate constant-parameter BVAR forecasts using the shared RNG.
      type(bvartools_bvar_draws_t), intent(in) :: draws !! Retained coefficient and covariance draws.
      real(dp), intent(in) :: history(:, :) !! Endogenous history, ordered from oldest to newest observation.
      integer, intent(in) :: lags !! Number of endogenous lags in each coefficient draw.
      integer, intent(in) :: horizon !! Number of forecast periods.
      real(dp), intent(in), optional :: future_regressors(:, :) !! Future exogenous and deterministic regressors.
      real(dp), intent(in), optional :: structural_impact(:, :, :) !! Draw-specific inverse contemporaneous matrices.
      real(dp), intent(in), optional :: probability !! Equal-tail posterior interval probability.
      type(bvartools_predictive_t) :: out
      real(dp), allocatable :: standard_normals(:, :, :)
      integer :: variable, step, draw

      if (horizon < 1 .or. draws%retained_draws < 1) then
         out%info = 1
         return
      end if
      allocate(standard_normals(size(history, 1), horizon, draws%retained_draws))
      do draw = 1, draws%retained_draws
         do step = 1, horizon
            do variable = 1, size(history, 1)
               standard_normals(variable, step, draw) = random_standard_normal()
            end do
         end do
      end do
      out = bvartools_bvar_predictive_from_random(draws, history, lags, &
         standard_normals, future_regressors, structural_impact, probability)
   end function bvartools_bvar_predictive

   pure function bvartools_tvp_bvar_predictive_from_random(draws, history, &
      lags, observations, standard_normals, future_regressors, &
      structural_impact, probability) result(out)
      !! Simulate TVP-BVAR forecasts using each draw's terminal parameter state.
      type(bvartools_tvp_bvar_draws_t), intent(in) :: draws !! Retained TVP coefficient and covariance draws.
      real(dp), intent(in) :: history(:, :) !! Endogenous history, ordered from oldest to newest observation.
      integer, intent(in) :: lags !! Number of endogenous lags in each coefficient state.
      integer, intent(in) :: observations !! Number of coefficient-state observations stored per draw.
      real(dp), intent(in) :: standard_normals(:, :, :) !! Independent normals by variable, horizon, and draw.
      real(dp), intent(in), optional :: future_regressors(:, :) !! Future exogenous and deterministic regressors.
      real(dp), intent(in), optional :: structural_impact(:, :, :) !! Draw-specific inverse contemporaneous matrices.
      real(dp), intent(in), optional :: probability !! Equal-tail posterior interval probability.
      type(bvartools_predictive_t) :: out
      real(dp), allocatable :: terminal_coefficients(:, :), terminal_covariance(:, :, :)
      integer :: states, draws_count, first, draw

      draws_count = draws%retained_draws
      if (.not. allocated(draws%states)) then
         out%info = 1
         return
      end if
      if (observations < 1 .or. draws_count < 1 .or. &
         size(draws%states, 2) /= draws_count .or. &
         mod(size(draws%states, 1), observations) /= 0) then
         out%info = 1
         return
      end if
      states = size(draws%states, 1)/observations
      first = (observations - 1)*states + 1
      terminal_coefficients = draws%states(first:first + states - 1, :)
      allocate(terminal_covariance(size(history, 1), size(history, 1), draws_count))
      if (allocated(draws%time_covariance)) then
         if (size(draws%time_covariance, 1) /= size(history, 1) .or. &
            size(draws%time_covariance, 2) /= size(history, 1) .or. &
            size(draws%time_covariance, 3) < observations .or. &
            size(draws%time_covariance, 4) /= draws_count) then
            out%info = 2
            return
         end if
         terminal_covariance = draws%time_covariance(:, :, observations, :)
      else if (allocated(draws%covariance)) then
         if (any(shape(draws%covariance) /= &
            [size(history, 1), size(history, 1), draws_count])) then
            out%info = 3
            return
         end if
         terminal_covariance = draws%covariance
      else
         out%info = 4
         return
      end if
      out = predictive_kernel(terminal_coefficients, terminal_covariance, &
         history, lags, standard_normals, future_regressors, structural_impact, &
         probability)
      if (draws%info /= 0 .and. out%info == 0) out%info = 10
   end function bvartools_tvp_bvar_predictive_from_random

   function bvartools_tvp_bvar_predictive(draws, history, lags, observations, &
      horizon, future_regressors, structural_impact, probability) result(out)
      !! Simulate terminal-state TVP-BVAR forecasts using the shared RNG.
      type(bvartools_tvp_bvar_draws_t), intent(in) :: draws !! Retained TVP coefficient and covariance draws.
      real(dp), intent(in) :: history(:, :) !! Endogenous history, ordered from oldest to newest observation.
      integer, intent(in) :: lags !! Number of endogenous lags in each coefficient state.
      integer, intent(in) :: observations !! Number of coefficient-state observations stored per draw.
      integer, intent(in) :: horizon !! Number of forecast periods.
      real(dp), intent(in), optional :: future_regressors(:, :) !! Future exogenous and deterministic regressors.
      real(dp), intent(in), optional :: structural_impact(:, :, :) !! Draw-specific inverse contemporaneous matrices.
      real(dp), intent(in), optional :: probability !! Equal-tail posterior interval probability.
      type(bvartools_predictive_t) :: out
      real(dp), allocatable :: standard_normals(:, :, :)
      integer :: variable, step, draw

      if (horizon < 1 .or. draws%retained_draws < 1) then
         out%info = 1
         return
      end if
      allocate(standard_normals(size(history, 1), horizon, draws%retained_draws))
      do draw = 1, draws%retained_draws
         do step = 1, horizon
            do variable = 1, size(history, 1)
               standard_normals(variable, step, draw) = random_standard_normal()
            end do
         end do
      end do
      out = bvartools_tvp_bvar_predictive_from_random(draws, history, lags, &
         observations, standard_normals, future_regressors, structural_impact, &
         probability)
   end function bvartools_tvp_bvar_predictive

   pure function bvartools_bvar_irf(draws, lags, horizon, response_type, shock, &
      scale_by_standard_deviation, structural_impact, cumulative, probability) &
      result(out)
      !! Compute posterior impulse responses for a constant-parameter BVAR.
      type(bvartools_bvar_draws_t), intent(in) :: draws !! Retained coefficient and covariance draws.
      integer, intent(in) :: lags !! Number of endogenous lags in each coefficient draw.
      integer, intent(in) :: horizon !! Largest impulse-response horizon, including responses from lag zero.
      character(len=*), intent(in) :: response_type !! Response definition: feir, oir, gir, sir, or sgir.
      real(dp), intent(in), optional :: shock(:) !! Signed shock magnitude for each impulse variable.
      logical, intent(in), optional :: scale_by_standard_deviation !! Scale shocks by draw-specific innovation deviations.
      real(dp), intent(in), optional :: structural_impact(:, :, :) !! Draw-specific inverse contemporaneous matrices.
      logical, intent(in), optional :: cumulative !! Accumulate responses from lag zero through each horizon.
      real(dp), intent(in), optional :: probability !! Equal-tail posterior interval probability.
      type(bvartools_irf_t) :: out

      if (.not. allocated(draws%coefficients) .or. &
         .not. allocated(draws%covariance)) then
         out%info = 1
         return
      end if
      out = irf_kernel(draws%coefficients, draws%covariance, lags, horizon, &
         response_type, shock, scale_by_standard_deviation, structural_impact, &
         cumulative, probability)
      if (draws%info /= 0 .and. out%info == 0) out%info = 10
   end function bvartools_bvar_irf

   pure function bvartools_tvp_bvar_irf(draws, variables, lags, observations, &
      period, horizon, response_type, shock, scale_by_standard_deviation, &
      structural_impact, cumulative, probability) result(out)
      !! Compute posterior impulse responses from one stored TVP-BVAR period.
      type(bvartools_tvp_bvar_draws_t), intent(in) :: draws !! Retained TVP coefficient and covariance draws.
      integer, intent(in) :: variables !! Number of endogenous variables.
      integer, intent(in) :: lags !! Number of endogenous lags in each coefficient state.
      integer, intent(in) :: observations !! Number of coefficient-state observations stored per draw.
      integer, intent(in) :: period !! Stored observation used for the impulse responses.
      integer, intent(in) :: horizon !! Largest impulse-response horizon, including responses from lag zero.
      character(len=*), intent(in) :: response_type !! Response definition: feir, oir, gir, sir, or sgir.
      real(dp), intent(in), optional :: shock(:) !! Signed shock magnitude for each impulse variable.
      logical, intent(in), optional :: scale_by_standard_deviation !! Scale shocks by draw-specific innovation deviations.
      real(dp), intent(in), optional :: structural_impact(:, :, :) !! Draw-specific inverse contemporaneous matrices.
      logical, intent(in), optional :: cumulative !! Accumulate responses from lag zero through each horizon.
      real(dp), intent(in), optional :: probability !! Equal-tail posterior interval probability.
      type(bvartools_irf_t) :: out
      real(dp), allocatable :: period_coefficients(:, :), period_covariance(:, :, :)
      integer :: states, draws_count, first

      if (.not. allocated(draws%states)) then
         out%info = 1
         return
      end if
      draws_count = draws%retained_draws
      if (variables < 1 .or. observations < 1 .or. period < 1 .or. &
         period > observations .or. draws_count < 1 .or. &
         size(draws%states, 2) /= draws_count .or. &
         mod(size(draws%states, 1), observations) /= 0) then
         out%info = 2
         return
      end if
      states = size(draws%states, 1)/observations
      first = (period - 1)*states + 1
      period_coefficients = draws%states(first:first + states - 1, :)
      allocate(period_covariance(variables, variables, draws_count))
      if (allocated(draws%time_covariance)) then
         if (any(shape(draws%time_covariance) /= &
            [variables, variables, observations, draws_count])) then
            out%info = 3
            return
         end if
         period_covariance = draws%time_covariance(:, :, period, :)
      else if (allocated(draws%covariance)) then
         if (any(shape(draws%covariance) /= &
            [variables, variables, draws_count])) then
            out%info = 4
            return
         end if
         period_covariance = draws%covariance
      else
         out%info = 5
         return
      end if
      out = irf_kernel(period_coefficients, period_covariance, lags, horizon, &
         response_type, shock, scale_by_standard_deviation, structural_impact, &
         cumulative, probability)
      if (draws%info /= 0 .and. out%info == 0) out%info = 10
   end function bvartools_tvp_bvar_irf

   pure function bvartools_bvar_fevd(draws, lags, horizon, response_type, &
      normalize_generalized, structural_impact, probability) result(out)
      !! Compute posterior FEVDs for a constant-parameter BVAR.
      type(bvartools_bvar_draws_t), intent(in) :: draws !! Retained coefficient and covariance draws.
      integer, intent(in) :: lags !! Number of endogenous lags in each coefficient draw.
      integer, intent(in) :: horizon !! Largest forecast-error horizon, including horizon zero.
      character(len=*), intent(in) :: response_type !! Decomposition definition: oir, gir, sir, or sgir.
      logical, intent(in), optional :: normalize_generalized !! Normalize generalized shares to sum to one.
      real(dp), intent(in), optional :: structural_impact(:, :, :) !! Draw-specific inverse contemporaneous matrices.
      real(dp), intent(in), optional :: probability !! Equal-tail posterior interval probability.
      type(bvartools_fevd_t) :: out

      if (.not. allocated(draws%coefficients) .or. &
         .not. allocated(draws%covariance)) then
         out%info = 1
         return
      end if
      out = fevd_kernel(draws%coefficients, draws%covariance, lags, horizon, &
         response_type, normalize_generalized, structural_impact, probability)
      if (draws%info /= 0 .and. out%info == 0) out%info = 10
   end function bvartools_bvar_fevd

   pure function bvartools_tvp_bvar_fevd(draws, variables, lags, observations, &
      period, horizon, response_type, normalize_generalized, structural_impact, &
      probability) result(out)
      !! Compute posterior FEVDs from one stored TVP-BVAR period.
      type(bvartools_tvp_bvar_draws_t), intent(in) :: draws !! Retained TVP coefficient and covariance draws.
      integer, intent(in) :: variables !! Number of endogenous variables.
      integer, intent(in) :: lags !! Number of endogenous lags in each coefficient state.
      integer, intent(in) :: observations !! Number of coefficient-state observations stored per draw.
      integer, intent(in) :: period !! Stored observation used for the decomposition.
      integer, intent(in) :: horizon !! Largest forecast-error horizon, including horizon zero.
      character(len=*), intent(in) :: response_type !! Decomposition definition: oir, gir, sir, or sgir.
      logical, intent(in), optional :: normalize_generalized !! Normalize generalized shares to sum to one.
      real(dp), intent(in), optional :: structural_impact(:, :, :) !! Draw-specific inverse contemporaneous matrices.
      real(dp), intent(in), optional :: probability !! Equal-tail posterior interval probability.
      type(bvartools_fevd_t) :: out
      real(dp), allocatable :: period_coefficients(:, :), period_covariance(:, :, :)
      integer :: states, draws_count, first

      if (.not. allocated(draws%states)) then
         out%info = 1
         return
      end if
      draws_count = draws%retained_draws
      if (variables < 1 .or. observations < 1 .or. period < 1 .or. &
         period > observations .or. draws_count < 1 .or. &
         size(draws%states, 2) /= draws_count .or. &
         mod(size(draws%states, 1), observations) /= 0) then
         out%info = 2
         return
      end if
      states = size(draws%states, 1)/observations
      first = (period - 1)*states + 1
      period_coefficients = draws%states(first:first + states - 1, :)
      allocate(period_covariance(variables, variables, draws_count))
      if (allocated(draws%time_covariance)) then
         if (any(shape(draws%time_covariance) /= &
            [variables, variables, observations, draws_count])) then
            out%info = 3
            return
         end if
         period_covariance = draws%time_covariance(:, :, period, :)
      else if (allocated(draws%covariance)) then
         if (any(shape(draws%covariance) /= &
            [variables, variables, draws_count])) then
            out%info = 4
            return
         end if
         period_covariance = draws%covariance
      else
         out%info = 5
         return
      end if
      out = fevd_kernel(period_coefficients, period_covariance, lags, horizon, &
         response_type, normalize_generalized, structural_impact, probability)
      if (draws%info /= 0 .and. out%info == 0) out%info = 10
   end function bvartools_tvp_bvar_fevd

   pure function bvartools_vecm_level_ar(pi, gamma) result(ar)
      !! Convert VECM error-correction and differenced-lag blocks to level-VAR blocks.
      real(dp), intent(in) :: pi(:, :) !! Endogenous error-correction coefficient matrix.
      real(dp), intent(in) :: gamma(:, :, :) !! Differenced endogenous lag coefficients.
      real(dp), allocatable :: ar(:, :, :)
      integer :: variables, level_lags, lag

      variables = size(pi, 1)
      level_lags = size(gamma, 3) + 1
      if (variables < 1 .or. size(pi, 2) /= variables .or. &
         size(gamma, 1) /= variables .or. size(gamma, 2) /= variables) then
         allocate(ar(0, 0, 0))
         return
      end if
      allocate(ar(variables, variables, level_lags))
      if (level_lags == 1) then
         ar(:, :, 1) = identity_matrix(variables) + pi
         return
      end if
      ar(:, :, 1) = identity_matrix(variables) + pi + gamma(:, :, 1)
      do lag = 2, level_lags - 1
         ar(:, :, lag) = gamma(:, :, lag) - gamma(:, :, lag - 1)
      end do
      ar(:, :, level_lags) = -gamma(:, :, level_lags - 1)
   end function bvartools_vecm_level_ar

   pure function bvartools_vecm_level_exogenous(pi_x, differences) result(level)
      !! Convert exogenous error-correction and differenced terms to level lag blocks.
      real(dp), intent(in) :: pi_x(:, :) !! Exogenous error-correction coefficients on lag-one levels.
      real(dp), intent(in) :: differences(:, :, :) !! Coefficients on current and lagged exogenous differences.
      real(dp), allocatable :: level(:, :, :)
      integer :: equations, variables, difference_blocks, block

      equations = size(pi_x, 1)
      variables = size(pi_x, 2)
      difference_blocks = size(differences, 3)
      if (equations < 1 .or. variables < 1 .or. difference_blocks < 1 .or. &
         size(differences, 1) /= equations .or. &
         size(differences, 2) /= variables) then
         allocate(level(0, 0, 0))
         return
      end if
      allocate(level(equations, variables, difference_blocks + 1))
      level(:, :, 1) = differences(:, :, 1)
      if (difference_blocks == 1) then
         level(:, :, 2) = pi_x - differences(:, :, 1)
         return
      end if
      level(:, :, 2) = pi_x - differences(:, :, 1) + differences(:, :, 2)
      do block = 3, difference_blocks
         level(:, :, block) = differences(:, :, block) - &
            differences(:, :, block - 1)
      end do
      level(:, :, difference_blocks + 1) = -differences(:, :, difference_blocks)
   end function bvartools_vecm_level_exogenous

   pure function bvartools_bvec_to_level_var(draws, variables, level_lags, &
      exogenous_error_correction, exogenous_differences) result(out)
      !! Convert constant BVEC posterior draws to level-VAR coefficient draws.
      type(bvartools_bvec_draws_t), intent(in) :: draws !! Retained BVEC coefficient draws.
      integer, intent(in) :: variables !! Number of endogenous variables.
      integer, intent(in) :: level_lags !! Lag order of the resulting level VAR.
      real(dp), intent(in), optional :: exogenous_error_correction(:, :, :) !! Exogenous level coefficients by draw.
      real(dp), intent(in), optional :: exogenous_differences(:, :, :, :) !! Exogenous difference coefficients by lag and draw.
      type(bvartools_level_var_draws_t) :: out
      real(dp), allocatable :: pi_matrix(:, :), gamma_matrix(:, :)
      real(dp), allocatable :: gamma_endogenous(:, :, :)
      integer :: draws_count, pi_columns, gamma_columns, endogenous_columns
      integer :: restricted_count, unrestricted_count, exogenous_count
      integer :: exogenous_difference_count, draw

      if (.not. allocated(draws%pi)) then
         out%info = 1
         return
      end if
      draws_count = draws%retained_draws
      if (variables < 1 .or. level_lags < 1 .or. draws_count < 1 .or. &
         size(draws%pi, 2) /= draws_count .or. &
         mod(size(draws%pi, 1), variables) /= 0) then
         out%info = 2
         return
      end if
      pi_columns = size(draws%pi, 1)/variables
      restricted_count = pi_columns - variables
      if (restricted_count < 0) then
         out%info = 3
         return
      end if
      endogenous_columns = variables*(level_lags - 1)
      gamma_columns = 0
      if (allocated(draws%gamma)) then
         if (size(draws%gamma, 2) /= draws_count .or. &
            mod(size(draws%gamma, 1), variables) /= 0) then
            out%info = 4
            return
         end if
         gamma_columns = size(draws%gamma, 1)/variables
      end if
      unrestricted_count = gamma_columns - endogenous_columns
      if (unrestricted_count < 0) then
         out%info = 5
         return
      end if
      if (present(exogenous_error_correction) .neqv. &
         present(exogenous_differences)) then
         out%info = 6
         return
      end if
      exogenous_count = 0
      exogenous_difference_count = 0
      if (present(exogenous_error_correction)) then
         exogenous_count = size(exogenous_error_correction, 2)
         exogenous_difference_count = size(exogenous_differences, 3)
         if (any(shape(exogenous_error_correction) /= &
            [variables, exogenous_count, draws_count]) .or. &
            size(exogenous_differences, 1) /= variables .or. &
            size(exogenous_differences, 2) /= exogenous_count .or. &
            size(exogenous_differences, 4) /= draws_count .or. &
            exogenous_difference_count < 1) then
            out%info = 7
            return
         end if
      end if

      allocate(out%ar(variables, variables, level_lags, draws_count))
      allocate(out%restricted_deterministic(variables, restricted_count, draws_count))
      allocate(out%unrestricted_deterministic(variables, unrestricted_count, draws_count))
      if (exogenous_count > 0) allocate(out%exogenous(variables, exogenous_count, &
         exogenous_difference_count + 1, draws_count))
      allocate(pi_matrix(variables, pi_columns))
      allocate(gamma_matrix(variables, gamma_columns))
      allocate(gamma_endogenous(variables, variables, level_lags - 1))
      do draw = 1, draws_count
         pi_matrix = reshape(draws%pi(:, draw), [variables, pi_columns])
         if (gamma_columns > 0) gamma_matrix = reshape(draws%gamma(:, draw), &
            [variables, gamma_columns])
         if (level_lags > 1) gamma_endogenous = reshape( &
            gamma_matrix(:, :endogenous_columns), &
            [variables, variables, level_lags - 1])
         out%ar(:, :, :, draw) = bvartools_vecm_level_ar( &
            pi_matrix(:, :variables), gamma_endogenous)
         if (restricted_count > 0) out%restricted_deterministic(:, :, draw) = &
            pi_matrix(:, variables + 1:)
         if (unrestricted_count > 0) out%unrestricted_deterministic(:, :, draw) = &
            gamma_matrix(:, endogenous_columns + 1:)
         if (exogenous_count > 0) out%exogenous(:, :, :, draw) = &
            bvartools_vecm_level_exogenous( &
            exogenous_error_correction(:, :, draw), &
            exogenous_differences(:, :, :, draw))
      end do
      out%retained_draws = draws_count
      if (draws%info /= 0) out%info = 10
   end function bvartools_bvec_to_level_var

   pure function bvartools_tvp_bvec_to_level_var(draws, variables, level_lags, &
      observations, exogenous_error_correction, exogenous_differences) result(out)
      !! Convert all TVP-BVEC periods and posterior draws to level-VAR coefficients.
      type(bvartools_tvp_bvec_draws_t), intent(in) :: draws !! Retained TVP-BVEC coefficient draws.
      integer, intent(in) :: variables !! Number of endogenous variables.
      integer, intent(in) :: level_lags !! Lag order of the resulting level VAR.
      integer, intent(in) :: observations !! Number of stored coefficient-state observations.
      real(dp), intent(in), optional :: exogenous_error_correction(:, :, :, :) !! Exogenous level coefficients by period and draw.
      real(dp), intent(in), optional :: exogenous_differences(:, :, :, :, :) !! Exogenous differences by lag, period, and draw.
      type(bvartools_tvp_level_var_draws_t) :: out
      real(dp), allocatable :: pi_matrix(:, :), gamma_matrix(:, :)
      real(dp), allocatable :: gamma_endogenous(:, :, :)
      integer :: draws_count, pi_columns, gamma_columns, endogenous_columns
      integer :: restricted_count, unrestricted_count, exogenous_count
      integer :: exogenous_difference_count, draw, time, first, last

      if (.not. allocated(draws%pi)) then
         out%info = 1
         return
      end if
      draws_count = draws%retained_draws
      if (variables < 1 .or. level_lags < 1 .or. observations < 1 .or. &
         draws_count < 1 .or. size(draws%pi, 2) /= draws_count .or. &
         mod(size(draws%pi, 1), variables*observations) /= 0) then
         out%info = 2
         return
      end if
      pi_columns = size(draws%pi, 1)/(variables*observations)
      restricted_count = pi_columns - variables
      if (restricted_count < 0) then
         out%info = 3
         return
      end if
      endogenous_columns = variables*(level_lags - 1)
      gamma_columns = 0
      if (allocated(draws%gamma)) then
         if (size(draws%gamma, 2) /= draws_count .or. &
            mod(size(draws%gamma, 1), variables*observations) /= 0) then
            out%info = 4
            return
         end if
         gamma_columns = size(draws%gamma, 1)/(variables*observations)
      end if
      unrestricted_count = gamma_columns - endogenous_columns
      if (unrestricted_count < 0) then
         out%info = 5
         return
      end if
      if (present(exogenous_error_correction) .neqv. &
         present(exogenous_differences)) then
         out%info = 6
         return
      end if
      exogenous_count = 0
      exogenous_difference_count = 0
      if (present(exogenous_error_correction)) then
         exogenous_count = size(exogenous_error_correction, 2)
         exogenous_difference_count = size(exogenous_differences, 3)
         if (any(shape(exogenous_error_correction) /= &
            [variables, exogenous_count, observations, draws_count]) .or. &
            size(exogenous_differences, 1) /= variables .or. &
            size(exogenous_differences, 2) /= exogenous_count .or. &
            size(exogenous_differences, 4) /= observations .or. &
            size(exogenous_differences, 5) /= draws_count .or. &
            exogenous_difference_count < 1) then
            out%info = 7
            return
         end if
      end if

      allocate(out%ar(variables, variables, level_lags, observations, draws_count))
      allocate(out%restricted_deterministic(variables, restricted_count, &
         observations, draws_count))
      allocate(out%unrestricted_deterministic(variables, unrestricted_count, &
         observations, draws_count))
      if (exogenous_count > 0) allocate(out%exogenous(variables, exogenous_count, &
         exogenous_difference_count + 1, observations, draws_count))
      allocate(pi_matrix(variables, pi_columns))
      allocate(gamma_matrix(variables, gamma_columns))
      allocate(gamma_endogenous(variables, variables, level_lags - 1))
      do draw = 1, draws_count
         do time = 1, observations
            first = (time - 1)*variables*pi_columns + 1
            last = time*variables*pi_columns
            pi_matrix = reshape(draws%pi(first:last, draw), &
               [variables, pi_columns])
            if (gamma_columns > 0) then
               first = (time - 1)*variables*gamma_columns + 1
               last = time*variables*gamma_columns
               gamma_matrix = reshape(draws%gamma(first:last, draw), &
                  [variables, gamma_columns])
            end if
            if (level_lags > 1) gamma_endogenous = reshape( &
               gamma_matrix(:, :endogenous_columns), &
               [variables, variables, level_lags - 1])
            out%ar(:, :, :, time, draw) = bvartools_vecm_level_ar( &
               pi_matrix(:, :variables), gamma_endogenous)
            if (restricted_count > 0) &
               out%restricted_deterministic(:, :, time, draw) = &
               pi_matrix(:, variables + 1:)
            if (unrestricted_count > 0) &
               out%unrestricted_deterministic(:, :, time, draw) = &
               gamma_matrix(:, endogenous_columns + 1:)
            if (exogenous_count > 0) out%exogenous(:, :, :, time, draw) = &
               bvartools_vecm_level_exogenous( &
               exogenous_error_correction(:, :, time, draw), &
               exogenous_differences(:, :, :, time, draw))
         end do
      end do
      out%observations = observations
      out%retained_draws = draws_count
      if (draws%info /= 0) out%info = 10
   end function bvartools_tvp_bvec_to_level_var

   pure function bvartools_reconstruct_levels(initial_level, differences) result(levels)
      !! Reconstruct level observations from an initial value and first differences.
      real(dp), intent(in) :: initial_level(:) !! Level immediately preceding the first difference.
      real(dp), intent(in) :: differences(:, :) !! First differences by variable and observation.
      real(dp), allocatable :: levels(:, :)
      integer :: variables, observations, time

      variables = size(initial_level)
      observations = size(differences, 2)
      if (variables < 1 .or. size(differences, 1) /= variables) then
         allocate(levels(0, 0))
         return
      end if
      allocate(levels(variables, observations + 1))
      levels(:, 1) = initial_level
      do time = 1, observations
         levels(:, time + 1) = levels(:, time) + differences(:, time)
      end do
   end function bvartools_reconstruct_levels

   pure function bvartools_prepare_var(series, lag_order, exogenous, &
      exogenous_order, include_constant, include_trend, seasonal_period, &
      holdout, time_varying, structural) result(out)
      !! Construct aligned BVAR response, regressor, and SUR matrices from raw series.
      real(dp), intent(in) :: series(:, :) !! Endogenous variables by chronological observation.
      integer, intent(in) :: lag_order !! Number of endogenous lags.
      real(dp), intent(in), optional :: exogenous(:, :) !! Exogenous variables by chronological observation.
      integer, intent(in), optional :: exogenous_order !! Largest exogenous lag, including current values.
      logical, intent(in), optional :: include_constant !! Include an unrestricted constant regressor.
      logical, intent(in), optional :: include_trend !! Include an unrestricted linear trend regressor.
      integer, intent(in), optional :: seasonal_period !! Period used for seasonal indicator regressors.
      integer, intent(in), optional :: holdout !! Number of aligned observations reserved for testing.
      logical, intent(in), optional :: time_varying !! Construct an observation-specific SUR design.
      logical, intent(in), optional :: structural !! Construct the contemporaneous structural design.
      type(bvartools_var_data_t) :: out
      real(dp), allocatable :: response_rows(:, :), ar_rows(:, :), exogenous_rows(:, :)
      real(dp), allocatable :: full_x(:, :)
      integer :: variables, observations, x_variables, x_order, leading
      integer :: usable, training, test_count, deterministic_count, period
      integer :: row, time, phase
      logical :: constant, trend, tvp, use_structural

      variables = size(series, 1)
      observations = size(series, 2)
      x_order = 0
      if (present(exogenous_order)) x_order = exogenous_order
      constant = .true.
      if (present(include_constant)) constant = include_constant
      trend = .false.
      if (present(include_trend)) trend = include_trend
      period = 0
      if (present(seasonal_period)) period = seasonal_period
      test_count = 0
      if (present(holdout)) test_count = holdout
      tvp = .false.
      if (present(time_varying)) tvp = time_varying
      use_structural = .false.
      if (present(structural)) use_structural = structural
      if (variables < 1 .or. lag_order < 1 .or. observations <= lag_order .or. &
         x_order < 0 .or. period < 0 .or. period == 1 .or. test_count < 0) then
         out%info = 1
         return
      end if
      if (present(exogenous)) then
         if (size(exogenous, 2) /= observations .or. size(exogenous, 1) < 1 .or. &
            observations <= max(lag_order, x_order)) then
            out%info = 2
            return
         end if
         call build_varx_data(transpose(series), transpose(exogenous), lag_order, &
            x_order, response_rows, ar_rows, exogenous_rows, .true.)
         x_variables = size(exogenous, 1)*(x_order + 1)
         leading = max(lag_order, x_order)
      else
         if (x_order /= 0) then
            out%info = 4
            return
         end if
         call build_var_data(transpose(series), lag_order, response_rows, ar_rows)
         allocate(exogenous_rows(size(response_rows, 1), 0))
         x_variables = 0
         leading = lag_order
      end if
      usable = size(response_rows, 1)
      if (test_count >= usable) then
         out%info = 3
         return
      end if
      deterministic_count = merge(1, 0, constant) + merge(1, 0, trend)
      if (period >= 2) deterministic_count = deterministic_count + period - 1
      allocate(full_x(variables*lag_order + x_variables + deterministic_count, usable))
      full_x(:variables*lag_order, :) = transpose(ar_rows)
      if (x_variables > 0) full_x(variables*lag_order + 1: &
         variables*lag_order + x_variables, :) = transpose(exogenous_rows)
      row = variables*lag_order + x_variables
      if (constant) then
         row = row + 1
         full_x(row, :) = 1.0_dp
      end if
      if (trend) then
         row = row + 1
         full_x(row, :) = [(real(time, dp), time=1, usable)]
      end if
      if (period >= 2) then
         full_x(row + 1:row + period - 1, :) = 0.0_dp
         do time = 1, usable
            phase = modulo(leading + time - 1, period) + 1
            if (phase < period) full_x(row + phase, time) = 1.0_dp
         end do
      end if
      training = usable - test_count
      allocate(out%y(variables, training), out%x(size(full_x, 1), training))
      allocate(out%holdout_y(variables, test_count), &
         out%holdout_x(size(full_x, 1), test_count))
      out%y = transpose(response_rows(:training, :))
      out%x = full_x(:, :training)
      if (test_count > 0) then
         out%holdout_y = transpose(response_rows(training + 1:, :))
         out%holdout_x = full_x(:, training + 1:)
      end if
      out%sur = make_sur_design(out%x, variables)
      if (tvp) out%tvp_sur = bvartools_sur_const_to_tvp(out%sur, variables, training)
      if (use_structural) out%structural = make_structural_design(out%y)
      out%observations = training
      out%variables = variables
      out%lag_order = lag_order
      out%endogenous_columns = variables*lag_order
      out%exogenous_columns = x_variables
      out%deterministic_columns = deterministic_count
   end function bvartools_prepare_var

   pure function bvartools_prepare_vecm(series, level_lags, exogenous, &
      exogenous_order, constant, trend, seasonal, seasonal_period, holdout, &
      time_varying, structural) result(out)
      !! Construct aligned BVEC response, error-correction, and short-run matrices.
      real(dp), intent(in) :: series(:, :) !! Endogenous levels by chronological observation.
      integer, intent(in) :: level_lags !! Lag order of the corresponding level VAR.
      real(dp), intent(in), optional :: exogenous(:, :) !! Exogenous levels by chronological observation.
      integer, intent(in), optional :: exogenous_order !! Number of exogenous difference blocks.
      character(len=*), intent(in), optional :: constant !! Constant placement: none, restricted, or unrestricted.
      character(len=*), intent(in), optional :: trend !! Trend placement: none, restricted, or unrestricted.
      character(len=*), intent(in), optional :: seasonal !! Seasonal placement: none, restricted, or unrestricted.
      integer, intent(in), optional :: seasonal_period !! Period used for seasonal indicator regressors.
      integer, intent(in), optional :: holdout !! Number of aligned observations reserved for testing.
      logical, intent(in), optional :: time_varying !! Construct an observation-specific short-run SUR design.
      logical, intent(in), optional :: structural !! Construct the contemporaneous structural design.
      type(bvartools_vecm_data_t) :: out
      real(dp), allocatable :: full_y(:, :), full_w(:, :), full_x(:, :)
      character(len=12) :: constant_location, trend_location, seasonal_location
      integer :: variables, observations, x_variables, x_order, leading, usable
      integer :: restricted_count, unrestricted_count, short_run_count
      integer :: training, test_count, period, time, lag, row, phase, original_time
      logical :: tvp, use_structural

      variables = size(series, 1)
      observations = size(series, 2)
      x_variables = 0
      if (present(exogenous)) x_variables = size(exogenous, 1)
      x_order = 0
      if (present(exogenous)) x_order = 1
      if (present(exogenous_order)) x_order = exogenous_order
      constant_location = 'none'
      if (present(constant)) constant_location = trim(constant)
      trend_location = 'none'
      if (present(trend)) trend_location = trim(trend)
      seasonal_location = 'none'
      if (present(seasonal)) seasonal_location = trim(seasonal)
      period = 0
      if (present(seasonal_period)) period = seasonal_period
      test_count = 0
      if (present(holdout)) test_count = holdout
      tvp = .false.
      if (present(time_varying)) tvp = time_varying
      use_structural = .false.
      if (present(structural)) use_structural = structural
      if (variables < 1 .or. level_lags < 1 .or. observations <= level_lags .or. &
         x_order < 0 .or. test_count < 0 .or. period < 0 .or. period == 1 .or. &
         .not. valid_deterministic_location(constant_location) .or. &
         .not. valid_deterministic_location(trend_location) .or. &
         .not. valid_deterministic_location(seasonal_location)) then
         out%info = 1
         return
      end if
      if (present(exogenous)) then
         if (size(exogenous, 2) /= observations .or. x_variables < 1 .or. &
            x_order < 1) then
            out%info = 2
            return
         end if
      else if (x_order /= 0) then
         out%info = 3
         return
      end if
      if (trim(seasonal_location) /= 'none' .and. period < 2) then
         out%info = 4
         return
      end if
      leading = max(level_lags, x_order)
      usable = observations - leading
      if (test_count >= usable) then
         out%info = 5
         return
      end if
      restricted_count = deterministic_count(constant_location, trend_location, &
         seasonal_location, period, 'restricted')
      unrestricted_count = deterministic_count(constant_location, trend_location, &
         seasonal_location, period, 'unrestricted')
      short_run_count = variables*(level_lags - 1) + x_variables*x_order
      allocate(full_y(variables, usable))
      allocate(full_w(variables + x_variables + restricted_count, usable))
      allocate(full_x(short_run_count + unrestricted_count, usable))
      full_w = 0.0_dp
      full_x = 0.0_dp
      do time = 1, usable
         original_time = leading + time
         full_y(:, time) = series(:, original_time) - series(:, original_time - 1)
         full_w(:variables, time) = series(:, original_time - 1)
         row = variables
         if (x_variables > 0) then
            full_w(row + 1:row + x_variables, time) = &
               exogenous(:, original_time - 1)
            row = row + x_variables
         end if
         call fill_deterministic_column(full_w(:, time), row, restricted_count, &
            constant_location, trend_location, seasonal_location, period, &
            time, original_time, 'restricted')
         row = 0
         do lag = 1, level_lags - 1
            full_x(row + 1:row + variables, time) = &
               series(:, original_time - lag) - &
               series(:, original_time - lag - 1)
            row = row + variables
         end do
         if (x_variables > 0) then
            do lag = 0, x_order - 1
               full_x(row + 1:row + x_variables, time) = &
                  exogenous(:, original_time - lag) - &
                  exogenous(:, original_time - lag - 1)
               row = row + x_variables
            end do
         end if
         call fill_deterministic_column(full_x(:, time), row, unrestricted_count, &
            constant_location, trend_location, seasonal_location, period, &
            time, original_time, 'unrestricted')
      end do
      training = usable - test_count
      allocate(out%y(variables, training), out%w(size(full_w, 1), training))
      allocate(out%x(size(full_x, 1), training))
      allocate(out%holdout_y(variables, test_count), &
         out%holdout_w(size(full_w, 1), test_count), &
         out%holdout_x(size(full_x, 1), test_count))
      out%y = full_y(:, :training)
      out%w = full_w(:, :training)
      out%x = full_x(:, :training)
      if (test_count > 0) then
         out%holdout_y = full_y(:, training + 1:)
         out%holdout_w = full_w(:, training + 1:)
         out%holdout_x = full_x(:, training + 1:)
      end if
      out%sur = make_sur_design(out%x, variables)
      if (tvp) out%tvp_sur = bvartools_sur_const_to_tvp(out%sur, variables, training)
      if (use_structural) out%structural = make_structural_design(out%y)
      out%observations = training
      out%variables = variables
      out%level_lags = level_lags
      out%restricted_deterministic_columns = restricted_count
      out%short_run_columns = short_run_count
      out%unrestricted_deterministic_columns = unrestricted_count
   end function bvartools_prepare_vecm

   pure function bvartools_prepare_dfm(series, lag_orders, factor_counts, &
      iterations, burnin) result(out)
      !! Standardize observations and enumerate the bvartools DFM model grid.
      real(dp), intent(in) :: series(:, :) !! Observed variables by chronological observation.
      integer, intent(in), optional :: lag_orders(:) !! Candidate factor-VAR lag orders.
      integer, intent(in), optional :: factor_counts(:) !! Candidate numbers of latent factors.
      integer, intent(in), optional :: iterations !! Number of retained Gibbs draws.
      integer, intent(in), optional :: burnin !! Number of discarded Gibbs draws.
      type(bvartools_dfm_data_t) :: out
      integer, allocatable :: orders(:), counts(:)
      real(dp) :: centered_sum
      integer :: variables, observations, models, variable, order, factor, model

      variables = size(series, 1)
      observations = size(series, 2)
      orders = [2]
      if (present(lag_orders)) orders = lag_orders
      counts = [1]
      if (present(factor_counts)) counts = factor_counts
      if (present(iterations)) out%iterations = iterations
      if (present(burnin)) out%burnin = burnin
      if (variables < 1 .or. observations < 2 .or. size(orders) < 1 .or. &
         size(counts) < 1 .or. any(orders < 0) .or. any(counts < 1) .or. &
         any(counts > variables) .or. out%iterations < 1 .or. out%burnin < 0) then
         out%info = 1
         return
      end if
      allocate(out%x(variables, observations), out%mean(variables), &
         out%standard_deviation(variables))
      do variable = 1, variables
         out%mean(variable) = sum(series(variable, :))/real(observations, dp)
         centered_sum = sum((series(variable, :) - out%mean(variable))**2)
         out%standard_deviation(variable) = &
            sqrt(centered_sum/real(observations - 1, dp))
         if (out%standard_deviation(variable) <= sqrt(tiny(1.0_dp))) then
            out%info = 2
            return
         end if
         out%x(variable, :) = (series(variable, :) - out%mean(variable))/ &
            out%standard_deviation(variable)
      end do
      models = size(orders)*size(counts)
      allocate(out%factor_count(models), out%lag_order(models))
      model = 0
      do factor = 1, size(counts)
         do order = 1, size(orders)
            model = model + 1
            out%factor_count(model) = counts(factor)
            out%lag_order(model) = orders(order)
         end do
      end do
      out%observations = observations
      out%variables = variables
   end function bvartools_prepare_dfm

   pure function bvartools_dfm_prior(variables, factors, lags, &
      loading_precision, measurement_shape, measurement_rate, &
      transition_precision, factor_shape, factor_rate) result(out)
      !! Construct the default or caller-scaled priors used by add_priors.dfmodel.
      integer, intent(in) :: variables !! Number of observed variables.
      integer, intent(in) :: factors !! Number of latent factors.
      integer, intent(in) :: lags !! Factor-VAR lag order.
      real(dp), intent(in), optional :: loading_precision !! Diagonal prior precision for free loadings.
      real(dp), intent(in), optional :: measurement_shape !! Measurement-variance inverse-gamma shape.
      real(dp), intent(in), optional :: measurement_rate !! Measurement-variance inverse-gamma rate.
      real(dp), intent(in), optional :: transition_precision !! Diagonal factor-transition prior precision.
      real(dp), intent(in), optional :: factor_shape !! Factor-variance inverse-gamma shape.
      real(dp), intent(in), optional :: factor_rate !! Factor-variance inverse-gamma rate.
      type(bvartools_dfm_prior_t) :: out
      real(dp) :: lambda_value, measurement_shape_value, measurement_rate_value
      real(dp) :: transition_value, factor_shape_value, factor_rate_value
      integer :: loading_count, transition_count, item

      lambda_value = 0.01_dp
      if (present(loading_precision)) lambda_value = loading_precision
      measurement_shape_value = 5.0_dp
      if (present(measurement_shape)) measurement_shape_value = measurement_shape
      measurement_rate_value = 4.0_dp
      if (present(measurement_rate)) measurement_rate_value = measurement_rate
      transition_value = 0.01_dp
      if (present(transition_precision)) transition_value = transition_precision
      factor_shape_value = 5.0_dp
      if (present(factor_shape)) factor_shape_value = factor_shape
      factor_rate_value = 4.0_dp
      if (present(factor_rate)) factor_rate_value = factor_rate
      if (variables < 1 .or. factors < 1 .or. factors > variables .or. lags < 0 .or. &
         lambda_value < 0.0_dp .or. measurement_shape_value < 0.0_dp .or. &
         measurement_rate_value <= 0.0_dp .or. transition_value < 0.0_dp .or. &
         factor_shape_value < 0.0_dp .or. factor_rate_value <= 0.0_dp) then
         out%info = 1
         return
      end if
      loading_count = (2*variables - factors - 1)*factors/2
      transition_count = factors*factors*lags
      allocate(out%loading_precision(loading_count, loading_count), &
         out%measurement_shape(variables), out%measurement_rate(variables), &
         out%transition_mean(transition_count), &
         out%transition_precision(transition_count, transition_count), &
         out%factor_shape(factors), out%factor_rate(factors))
      out%loading_precision = 0.0_dp
      do item = 1, loading_count
         out%loading_precision(item, item) = lambda_value
      end do
      out%measurement_shape = measurement_shape_value
      out%measurement_rate = measurement_rate_value
      out%transition_mean = 0.0_dp
      out%transition_precision = 0.0_dp
      do item = 1, transition_count
         out%transition_precision(item, item) = transition_value
      end do
      out%factor_shape = factor_shape_value
      out%factor_rate = factor_rate_value
   end function bvartools_dfm_prior

   function bvartools_dfm_grid_gibbs(data, loading_precision, &
      measurement_shape, measurement_rate, transition_precision, &
      factor_shape, factor_rate) result(out)
      !! Fit every dynamic factor specification prepared by bvartools_prepare_dfm.
      type(bvartools_dfm_data_t), intent(in) :: data !! Standardized observations and DFM model grid.
      real(dp), intent(in), optional :: loading_precision !! Diagonal prior precision for free loadings.
      real(dp), intent(in), optional :: measurement_shape !! Measurement-variance inverse-gamma shape.
      real(dp), intent(in), optional :: measurement_rate !! Measurement-variance inverse-gamma rate.
      real(dp), intent(in), optional :: transition_precision !! Diagonal factor-transition prior precision.
      real(dp), intent(in), optional :: factor_shape !! Factor-variance inverse-gamma shape.
      real(dp), intent(in), optional :: factor_rate !! Factor-variance inverse-gamma rate.
      type(bvartools_dfm_grid_draws_t) :: out
      real(dp) :: lambda_value, measurement_shape_value, measurement_rate_value
      real(dp) :: transition_value, factor_shape_value, factor_rate_value
      integer :: models, model

      lambda_value = 0.01_dp
      if (present(loading_precision)) lambda_value = loading_precision
      measurement_shape_value = 5.0_dp
      if (present(measurement_shape)) measurement_shape_value = measurement_shape
      measurement_rate_value = 4.0_dp
      if (present(measurement_rate)) measurement_rate_value = measurement_rate
      transition_value = 0.01_dp
      if (present(transition_precision)) transition_value = transition_precision
      factor_shape_value = 5.0_dp
      if (present(factor_shape)) factor_shape_value = factor_shape
      factor_rate_value = 4.0_dp
      if (present(factor_rate)) factor_rate_value = factor_rate
      if (data%info /= 0 .or. .not. allocated(data%x) .or. &
         .not. allocated(data%factor_count) .or. .not. allocated(data%lag_order) .or. &
         size(data%x, 1) /= data%variables .or. &
         size(data%x, 2) /= data%observations .or. &
         size(data%factor_count) /= size(data%lag_order) .or. &
         size(data%factor_count) < 1 .or. data%iterations < 1 .or. data%burnin < 0) then
         out%info = 1
         return
      end if
      models = size(data%factor_count)
      allocate(out%prior(models), out%draws(models), out%factor_count(models), &
         out%lag_order(models))
      out%factor_count = data%factor_count
      out%lag_order = data%lag_order
      do model = 1, models
         out%prior(model) = bvartools_dfm_prior(data%variables, &
            data%factor_count(model), data%lag_order(model), lambda_value, &
            measurement_shape_value, measurement_rate_value, transition_value, &
            factor_shape_value, factor_rate_value)
         if (out%prior(model)%info /= 0) then
            out%failed_model = model
            out%info = 2
            return
         end if
         out%draws(model) = bvartools_dfm_gibbs(data%x, &
            data%factor_count(model), data%lag_order(model), data%iterations, &
            data%burnin, out%prior(model)%loading_precision, &
            out%prior(model)%measurement_shape, out%prior(model)%measurement_rate, &
            out%prior(model)%transition_mean, &
            out%prior(model)%transition_precision, out%prior(model)%factor_shape, &
            out%prior(model)%factor_rate)
         if (out%draws(model)%info /= 0) then
            out%failed_model = model
            out%info = 3
            return
         end if
      end do
   end function bvartools_dfm_grid_gibbs

   pure function bvartools_model_comparison_constant(residual, covariance, &
      parameter_count) result(out)
      !! Compute bvartools posterior information criteria with constant covariances.
      real(dp), intent(in) :: residual(:, :, :) !! Residuals by variable, observation, and posterior draw.
      real(dp), intent(in) :: covariance(:, :, :) !! Constant covariance matrix by posterior draw.
      integer, intent(in) :: parameter_count !! Effective regression parameter count.
      type(bvartools_model_comparison_t) :: out
      real(dp), allocatable :: inverse(:, :), lower(:, :)
      real(dp) :: log_determinant, contribution, constant
      integer :: variables, observations, draws, draw, time, status

      variables = size(residual, 1)
      observations = size(residual, 2)
      draws = size(residual, 3)
      if (variables < 1 .or. observations < 2 .or. draws < 1 .or. &
         any(shape(covariance) /= [variables, variables, draws]) .or. &
         parameter_count < 0) then
         out%info = 1
         return
      end if
      allocate(out%observation_log_likelihood(observations))
      out%observation_log_likelihood = 0.0_dp
      constant = -0.5_dp*real(variables, dp)*log(2.0_dp*acos(-1.0_dp))
      do draw = 1, draws
         call cholesky_lower(covariance(:, :, draw), lower, status)
         if (status /= 0) then
            out%failed_draw = draw
            out%info = 2
            return
         end if
         log_determinant = 2.0_dp*sum(log([(lower(time, time), &
            time=1, variables)]))
         call invert_matrix(covariance(:, :, draw), inverse, status)
         if (status /= 0) then
            out%failed_draw = draw
            out%info = 2
            return
         end if
         do time = 1, observations
            contribution = constant - 0.5_dp*log_determinant - &
               0.5_dp*dot_product(residual(:, time, draw), &
               matmul(inverse, residual(:, time, draw)))
            out%observation_log_likelihood(time) = &
               out%observation_log_likelihood(time) + contribution/real(draws, dp)
         end do
      end do
      call finish_model_comparison(out, parameter_count, observations, draws)
   end function bvartools_model_comparison_constant

   pure function bvartools_model_comparison_tvp(residual, covariance, &
      parameter_count) result(out)
      !! Compute bvartools posterior information criteria with time-varying covariances.
      real(dp), intent(in) :: residual(:, :, :) !! Residuals by variable, observation, and posterior draw.
      real(dp), intent(in) :: covariance(:, :, :, :) !! Covariance by variable, observation, and posterior draw.
      integer, intent(in) :: parameter_count !! Effective regression parameter count.
      type(bvartools_model_comparison_t) :: out
      real(dp), allocatable :: inverse(:, :), lower(:, :)
      real(dp) :: log_determinant, contribution, constant
      integer :: variables, observations, draws, draw, time, diagonal, status

      variables = size(residual, 1)
      observations = size(residual, 2)
      draws = size(residual, 3)
      if (variables < 1 .or. observations < 2 .or. draws < 1 .or. &
         any(shape(covariance) /= [variables, variables, observations, draws]) .or. &
         parameter_count < 0) then
         out%info = 1
         return
      end if
      allocate(out%observation_log_likelihood(observations))
      out%observation_log_likelihood = 0.0_dp
      constant = -0.5_dp*real(variables, dp)*log(2.0_dp*acos(-1.0_dp))
      do draw = 1, draws
         do time = 1, observations
            call cholesky_lower(covariance(:, :, time, draw), lower, status)
            if (status /= 0) then
               out%failed_observation = time
               out%failed_draw = draw
               out%info = 2
               return
            end if
            log_determinant = 0.0_dp
            do diagonal = 1, variables
               log_determinant = log_determinant + 2.0_dp*log(lower(diagonal, diagonal))
            end do
            call invert_matrix(covariance(:, :, time, draw), inverse, status)
            if (status /= 0) then
               out%failed_observation = time
               out%failed_draw = draw
               out%info = 2
               return
            end if
            contribution = constant - 0.5_dp*log_determinant - &
               0.5_dp*dot_product(residual(:, time, draw), &
               matmul(inverse, residual(:, time, draw)))
            out%observation_log_likelihood(time) = &
               out%observation_log_likelihood(time) + contribution/real(draws, dp)
         end do
      end do
      call finish_model_comparison(out, parameter_count, observations, draws)
   end function bvartools_model_comparison_tvp

   pure function bvartools_compare_models(data) result(out)
      !! Compare a collection of posterior BVAR or BVEC likelihood data.
      type(bvartools_model_likelihood_data_t), intent(in) :: data(:) !! Residual, covariance, and parameter data by model.
      type(bvartools_model_comparison_set_t) :: out
      integer :: models, model, covariance_periods

      models = size(data)
      if (models < 1) then
         out%info = 1
         return
      end if
      allocate(out%model(models))
      do model = 1, models
         if (.not. allocated(data(model)%residual) .or. &
            .not. allocated(data(model)%covariance)) then
            out%failed_model = model
            out%info = 1
            return
         end if
         covariance_periods = size(data(model)%covariance, 3)
         if (covariance_periods == 1) then
            out%model(model) = bvartools_model_comparison( &
               data(model)%residual, data(model)%covariance(:, :, 1, :), &
               data(model)%parameter_count)
         else if (covariance_periods == size(data(model)%residual, 2)) then
            out%model(model) = bvartools_model_comparison( &
               data(model)%residual, data(model)%covariance, &
               data(model)%parameter_count)
         else
            out%failed_model = model
            out%info = 1
            return
         end if
         if (out%model(model)%info /= 0) then
            out%failed_model = model
            out%info = 2
            return
         end if
      end do
      out%best_aic = minloc([(out%model(model)%aic, model=1, models)], dim=1)
      out%best_bic = minloc([(out%model(model)%bic, model=1, models)], dim=1)
      out%best_hq = minloc([(out%model(model)%hq, model=1, models)], dim=1)
   end function bvartools_compare_models

   pure subroutine finish_model_comparison(out, parameter_count, observations, draws)
      !! Finish posterior likelihood aggregation and information criteria.
      type(bvartools_model_comparison_t), intent(inout) :: out !! Partially accumulated comparison result.
      integer, intent(in) :: parameter_count !! Effective regression parameter count.
      integer, intent(in) :: observations !! Number of likelihood observations.
      integer, intent(in) :: draws !! Number of posterior draws.

      out%log_likelihood = sum(out%observation_log_likelihood)
      out%parameter_count = parameter_count
      out%observations = observations
      out%draws = draws
      out%aic = 2.0_dp*real(parameter_count, dp) - 2.0_dp*out%log_likelihood
      out%bic = log(real(observations, dp))*real(parameter_count, dp) - &
         2.0_dp*out%log_likelihood
      out%hq = 2.0_dp*log(log(real(observations, dp)))* &
         real(parameter_count, dp) - 2.0_dp*out%log_likelihood
   end subroutine finish_model_comparison

   pure function bvartools_initial_state_posterior(first_state, state_precision, &
      prior_mean, prior_precision, persistence) result(out)
      !! Compute the Gaussian posterior of the state preceding a random walk.
      real(dp), intent(in) :: first_state(:) !! State at the first observation.
      real(dp), intent(in) :: state_precision(:, :) !! First state-innovation precision.
      real(dp), intent(in) :: prior_mean(:) !! Prior mean of the pre-sample state.
      real(dp), intent(in) :: prior_precision(:, :) !! Prior precision of the pre-sample state.
      real(dp), intent(in), optional :: persistence !! Scalar transition persistence.
      type(bvartools_normal_posterior_t) :: out
      real(dp), allocatable :: posterior_precision(:, :), inverse(:, :), score(:)
      integer :: states, status
      real(dp) :: rho

      states = size(first_state)
      rho = 1.0_dp
      if (present(persistence)) rho = persistence
      if (size(prior_mean) /= states .or. &
         any(shape(state_precision) /= [states, states]) .or. &
         any(shape(prior_precision) /= [states, states])) then
         out%info = 1
         return
      end if
      posterior_precision = prior_precision + rho*rho*state_precision
      score = matmul(prior_precision, prior_mean) + &
         rho*matmul(state_precision, first_state)
      call invert_matrix(posterior_precision, inverse, status)
      if (status /= 0) then
         out%info = 2
         return
      end if
      out%covariance = inverse
      out%mean = matmul(inverse, score)
   end function bvartools_initial_state_posterior

   function bvartools_tvp_bvar_gibbs(y, design, initial_state, &
      initial_covariance, iterations, burnin, initial_state_prior_mean, &
      initial_state_prior_precision, state_shape, state_rate, &
      covariance_prior_scale, covariance_prior_df, measurement_shape, &
      measurement_rate, inclusion_probability, selectable, &
      initial_log_variance, initial_log_variance_level, &
      initial_log_variance_state_variance, &
      log_variance_state_shape, log_variance_state_rate, &
      log_variance_initial_prior_mean, log_variance_initial_prior_precision, &
      log_variance_offset, volatility_method) result(out)
      !! Run the random-walk coefficient core of the bvartools TVP-BVAR sampler.
      real(dp), intent(in) :: y(:, :) !! Endogenous variables by observation.
      real(dp), intent(in) :: design(:, :) !! Observation-major constant SUR design.
      real(dp), intent(in) :: initial_state(:) !! Initial pre-sample coefficient state.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial measurement covariance.
      integer, intent(in) :: iterations !! Number of retained Gibbs draws.
      integer, intent(in) :: burnin !! Number of discarded Gibbs draws.
      real(dp), intent(in) :: initial_state_prior_mean(:) !! Prior mean of the pre-sample state.
      real(dp), intent(in) :: initial_state_prior_precision(:, :) !! Prior precision of the pre-sample state.
      real(dp), intent(in) :: state_shape(:) !! State-innovation gamma prior shapes.
      real(dp), intent(in) :: state_rate(:) !! State-innovation gamma prior rates.
      real(dp), intent(in), optional :: covariance_prior_scale(:, :) !! Inverse-Wishart measurement scale.
      real(dp), intent(in), optional :: covariance_prior_df !! Inverse-Wishart prior degrees of freedom.
      real(dp), intent(in), optional :: measurement_shape(:) !! Diagonal measurement-precision prior shapes.
      real(dp), intent(in), optional :: measurement_rate(:) !! Diagonal measurement-precision prior rates.
      real(dp), intent(in), optional :: inclusion_probability(:) !! BVS trajectory inclusion probabilities.
      integer, intent(in), optional :: selectable(:) !! One-based BVS trajectory positions.
      real(dp), intent(in), optional :: initial_log_variance(:, :) !! Initial log-variance paths.
      real(dp), intent(in), optional :: initial_log_variance_level(:) !! Initial pre-sample log variances.
      real(dp), intent(in), optional :: initial_log_variance_state_variance(:) !! Initial log-variance innovation variances.
      real(dp), intent(in), optional :: log_variance_state_shape(:) !! Log-variance innovation precision shapes.
      real(dp), intent(in), optional :: log_variance_state_rate(:) !! Log-variance innovation precision rates.
      real(dp), intent(in), optional :: log_variance_initial_prior_mean(:) !! Prior mean of pre-sample log variances.
      real(dp), intent(in), optional :: log_variance_initial_prior_precision(:, :) !! Prior precision of pre-sample log variances.
      real(dp), intent(in), optional :: log_variance_offset(:) !! Positive log-square transformation offsets.
      character(len=*), intent(in), optional :: volatility_method !! Log-variance mixture method: ksc or ocsn.
      type(bvartools_tvp_bvar_draws_t) :: out
      type(bvartools_normal_posterior_t) :: state_posterior, initial_posterior
      type(bvartools_gamma_posterior_t) :: variance_posterior
      type(bvartools_bvs_result_t) :: selection
      type(bvartools_stochastic_volatility_t) :: volatility_draw
      real(dp), allocatable :: state(:, :), state_vector(:), current_initial(:)
      real(dp), allocatable :: state_precision(:, :, :), measurement_precision(:, :, :)
      real(dp), allocatable :: precision(:, :), covariance(:, :), normal(:)
      real(dp), allocatable :: state_precision_diagonal(:), residual(:, :)
      real(dp), allocatable :: posterior_scale(:, :), measurement_diagonal(:)
      real(dp), allocatable :: active_design(:, :), proposal_uniforms(:)
      real(dp), allocatable :: acceptance_uniforms(:)
      logical, allocatable :: current_inclusion(:)
      real(dp), allocatable :: log_variance(:, :), log_variance_initial(:)
      real(dp), allocatable :: log_variance_state_variance(:)
      real(dp), allocatable :: log_variance_state_precision(:, :)
      real(dp), allocatable :: log_variance_normal(:)
      real(dp) :: posterior_df
      integer :: k, observations, states, total, draw, retained, time, item, status
      integer :: first, last
      logical :: full_covariance, use_bvs, use_sv, use_ocsn

      k = size(y, 1)
      observations = size(y, 2)
      states = size(initial_state)
      total = iterations + burnin
      full_covariance = present(covariance_prior_scale)
      use_bvs = present(inclusion_probability)
      use_sv = present(initial_log_variance)
      use_ocsn = .false.
      if (present(volatility_method)) then
         select case (trim(volatility_method))
         case ('ksc', 'KSC')
         case ('ocsn', 'OCSN')
            use_ocsn = .true.
         case default
            out%info = 19
            return
         end select
      end if
      if (k < 1 .or. observations < 1 .or. states < 1 .or. &
         size(design, 1) /= k*observations .or. size(design, 2) /= states .or. &
         any(shape(initial_covariance) /= [k, k]) .or. &
         size(initial_state_prior_mean) /= states .or. &
         any(shape(initial_state_prior_precision) /= [states, states]) .or. &
         size(state_shape) /= states .or. size(state_rate) /= states .or. &
         any(state_shape <= 0.0_dp) .or. any(state_rate <= 0.0_dp) .or. &
         iterations < 1 .or. burnin < 0) then
         out%info = 1
         return
      end if
      if (full_covariance) then
         if (.not. present(covariance_prior_df)) then
            out%info = 2
            return
         end if
         if (any(shape(covariance_prior_scale) /= [k, k]) .or. &
            covariance_prior_df <= real(k - 1, dp)) then
            out%info = 3
            return
         end if
      else if (.not. use_sv) then
         if (.not. present(measurement_shape) .or. &
            .not. present(measurement_rate)) then
            out%info = 4
            return
         end if
         if (size(measurement_shape) /= k .or. size(measurement_rate) /= k .or. &
            any(measurement_shape <= 0.0_dp) .or. any(measurement_rate <= 0.0_dp)) then
            out%info = 5
            return
         end if
      end if
      if (use_sv) then
         if (full_covariance .or. .not. present(initial_log_variance_level) .or. &
            .not. present(initial_log_variance_state_variance) .or. &
            .not. present(log_variance_state_shape) .or. &
            .not. present(log_variance_state_rate) .or. &
            .not. present(log_variance_initial_prior_mean) .or. &
            .not. present(log_variance_initial_prior_precision) .or. &
            .not. present(log_variance_offset)) then
            out%info = 14
            return
         end if
         if (any(shape(initial_log_variance) /= [k, observations]) .or. &
            size(initial_log_variance_level) /= k .or. &
            size(initial_log_variance_state_variance) /= k .or. &
            size(log_variance_state_shape) /= k .or. &
            size(log_variance_state_rate) /= k .or. &
            size(log_variance_initial_prior_mean) /= k .or. &
            any(shape(log_variance_initial_prior_precision) /= [k, k]) .or. &
            size(log_variance_offset) /= k .or. &
            any(log_variance_state_shape <= 0.0_dp) .or. &
            any(log_variance_state_rate <= 0.0_dp) .or. &
            any(initial_log_variance_state_variance <= 0.0_dp) .or. &
            any(log_variance_offset <= 0.0_dp)) then
            out%info = 15
            return
         end if
      else if (present(initial_log_variance_level) .or. &
         present(initial_log_variance_state_variance) .or. &
         present(log_variance_state_shape) .or. present(log_variance_state_rate) .or. &
         present(log_variance_initial_prior_mean) .or. &
         present(log_variance_initial_prior_precision) .or. &
         present(log_variance_offset) .or. present(volatility_method)) then
         out%info = 16
         return
      end if
      if (use_bvs) then
         if (size(inclusion_probability) /= states .or. &
            any(inclusion_probability < 0.0_dp) .or. &
            any(inclusion_probability > 1.0_dp)) then
            out%info = 10
            return
         end if
         if (present(selectable)) then
            if (any(selectable < 1) .or. any(selectable > states)) then
               out%info = 11
               return
            end if
         end if
      else if (present(selectable)) then
         out%info = 12
         return
      end if

      covariance = initial_covariance
      call invert_matrix(covariance, precision, status)
      if (status /= 0) then
         out%info = 6
         return
      end if
      current_initial = initial_state
      allocate(state(states, observations), state_precision(states, states, observations), &
         measurement_precision(k, k, observations), state_precision_diagonal(states), &
         normal(states*observations), out%states(states*observations, iterations), &
         out%state_variance(states, iterations), &
         out%initial_state(states, iterations), out%covariance(k, k, iterations))
      allocate(out%time_covariance(k, k, observations, iterations))
      active_design = design
      if (use_bvs) then
         allocate(current_inclusion(states), out%included(states, iterations))
         current_inclusion = .true.
         if (present(selectable)) then
            allocate(proposal_uniforms(size(selectable)), &
               acceptance_uniforms(size(selectable)))
         else
            allocate(proposal_uniforms(states), acceptance_uniforms(states))
         end if
      end if
      if (use_sv) then
         log_variance = initial_log_variance
         log_variance_initial = initial_log_variance_level
         log_variance_state_variance = initial_log_variance_state_variance
         allocate(log_variance_state_precision(k, k), log_variance_normal(k), &
            out%log_variance(k*observations, iterations), &
            out%log_variance_state_variance(k, iterations), &
            out%log_variance_initial(k, iterations))
         log_variance_state_precision = 0.0_dp
      end if
      state_precision_diagonal = 1.0_dp/state_rate
      state_precision = 0.0_dp
      measurement_precision = 0.0_dp
      do time = 1, observations
         do item = 1, states
            state_precision(item, item, time) = state_precision_diagonal(item)
         end do
         if (use_sv) then
            do item = 1, k
               measurement_precision(item, item, time) = &
                  exp(-log_variance(item, time))
            end do
         else
            measurement_precision(:, :, time) = precision
         end if
      end do
      posterior_df = 0.0_dp
      if (full_covariance) posterior_df = covariance_prior_df + real(observations, dp)
      retained = 0
      do draw = 1, total
         active_design = design
         if (use_bvs) then
            do item = 1, states
               if (.not. current_inclusion(item)) active_design(:, item) = 0.0_dp
            end do
         end if
         state_posterior = bvartools_random_walk_posterior(y, active_design, &
            measurement_precision, state_precision, current_initial)
         if (state_posterior%info /= 0) then
            out%info = 7
            return
         end if
         do item = 1, size(normal)
            normal(item) = random_standard_normal()
         end do
         state_vector = bvartools_normal_draw(state_posterior, normal)
         state = reshape(state_vector, [states, observations])
         if (use_bvs) then
            do item = 1, size(proposal_uniforms)
               proposal_uniforms(item) = random_uniform()
               acceptance_uniforms(item) = random_uniform()
            end do
            selection = bvartools_tvp_bvs(y, design, state, current_inclusion, &
               measurement_precision, inclusion_probability, proposal_uniforms, &
               acceptance_uniforms, selectable)
            if (selection%info /= 0) then
               out%info = 13
               return
            end if
            current_inclusion = selection%included
            do item = 1, states
               if (.not. current_inclusion(item)) state(item, :) = 0.0_dp
            end do
            state_vector = reshape(state, [states*observations])
         end if
         allocate(residual(k, observations))
         do time = 1, observations
            first = (time - 1)*k + 1
            last = time*k
            residual(:, time) = y(:, time) - &
               matmul(design(first:last, :), state(:, time))
         end do
         if (use_sv) then
            if (use_ocsn) then
               volatility_draw = bvartools_stochastic_volatility_ocsn(residual, &
                  log_variance, log_variance_state_variance, &
                  log_variance_initial, log_variance_offset)
            else
               volatility_draw = bvartools_stochastic_volatility(residual, &
                  log_variance, log_variance_state_variance, &
                  log_variance_initial, log_variance_offset)
            end if
            if (volatility_draw%info /= 0) then
               out%info = 17
               return
            end if
            log_variance = volatility_draw%log_variance
            variance_posterior = bvartools_state_variance_posterior(log_variance, &
               log_variance_initial, log_variance_state_shape, log_variance_state_rate)
            measurement_diagonal = bvartools_gamma_precision_draw(variance_posterior)
            log_variance_state_variance = &
               1.0_dp/max(measurement_diagonal, tiny(1.0_dp))
            log_variance_state_precision = 0.0_dp
            do item = 1, k
               log_variance_state_precision(item, item) = measurement_diagonal(item)
            end do
            initial_posterior = bvartools_initial_state_posterior( &
               log_variance(:, 1), log_variance_state_precision, &
               log_variance_initial_prior_mean, &
               log_variance_initial_prior_precision)
            if (initial_posterior%info /= 0) then
               out%info = 18
               return
            end if
            do item = 1, k
               log_variance_normal(item) = random_standard_normal()
            end do
            log_variance_initial = bvartools_normal_draw(initial_posterior, &
               log_variance_normal)
            covariance = 0.0_dp
            precision = 0.0_dp
            do item = 1, k
               covariance(item, item) = exp(log_variance(item, observations))
               precision(item, item) = exp(-log_variance(item, observations))
            end do
         else if (full_covariance) then
            posterior_scale = covariance_prior_scale + &
               matmul(residual, transpose(residual))
            call inverse_wishart_draw(posterior_scale, posterior_df, covariance, &
               precision, status)
            if (status /= 0) then
               out%info = 8
               return
            end if
         else
            variance_posterior = bvartools_measurement_variance_posterior( &
               residual, measurement_shape, measurement_rate)
            measurement_diagonal = bvartools_gamma_precision_draw(variance_posterior)
            covariance = 0.0_dp
            precision = 0.0_dp
            do item = 1, k
               precision(item, item) = measurement_diagonal(item)
               covariance(item, item) = 1.0_dp/measurement_diagonal(item)
            end do
         end if
         deallocate(residual)
         variance_posterior = bvartools_state_variance_posterior(state, &
            current_initial, state_shape, state_rate)
         state_precision_diagonal = bvartools_gamma_precision_draw(variance_posterior)
         state_precision(:, :, 1) = 0.0_dp
         do item = 1, states
            state_precision(item, item, 1) = state_precision_diagonal(item)
         end do
         initial_posterior = bvartools_initial_state_posterior(state(:, 1), &
            state_precision(:, :, 1), initial_state_prior_mean, &
            initial_state_prior_precision)
         if (initial_posterior%info /= 0) then
            out%info = 9
            return
         end if
         do item = 1, states
            normal(item) = random_standard_normal()
         end do
         current_initial = bvartools_normal_draw(initial_posterior, normal(1:states))
         do time = 1, observations
            state_precision(:, :, time) = 0.0_dp
            do item = 1, states
               state_precision(item, item, time) = state_precision_diagonal(item)
            end do
            measurement_precision(:, :, time) = 0.0_dp
            if (use_sv) then
               do item = 1, k
                  measurement_precision(item, item, time) = &
                     exp(-log_variance(item, time))
               end do
            else
               measurement_precision(:, :, time) = precision
            end if
         end do
         if (draw > burnin) then
            retained = retained + 1
            out%states(:, retained) = state_vector
            out%state_variance(:, retained) = &
               1.0_dp/max(state_precision_diagonal, tiny(1.0_dp))
            out%initial_state(:, retained) = current_initial
            out%covariance(:, :, retained) = covariance
            do time = 1, observations
               out%time_covariance(:, :, time, retained) = 0.0_dp
               if (use_sv) then
                  do item = 1, k
                     out%time_covariance(item, item, time, retained) = &
                        exp(log_variance(item, time))
                  end do
               else
                  out%time_covariance(:, :, time, retained) = covariance
               end if
            end do
            if (use_sv) then
               out%log_variance(:, retained) = &
                  reshape(log_variance, [k*observations])
               out%log_variance_state_variance(:, retained) = &
                  log_variance_state_variance
               out%log_variance_initial(:, retained) = log_variance_initial
            end if
            if (use_bvs) out%included(:, retained) = current_inclusion
         end if
      end do
      out%retained_draws = retained
   end function bvartools_tvp_bvar_gibbs

   function bvartools_tvp_bvec_gibbs(y, w, initial_alpha, initial_beta, &
      initial_covariance, iterations, burnin, initial_state_prior_mean, &
      initial_state_prior_precision, state_shape, state_rate, beta_persistence, &
      beta_state_precision, beta_initial_prior_mean, beta_initial_prior_precision, &
      covariance_prior_scale, covariance_prior_df, measurement_shape, &
      measurement_rate, x, initial_gamma, inclusion_probability, &
      selectable, initial_covariance_state, &
      covariance_state_initial_prior_mean, &
       covariance_state_initial_prior_precision, covariance_state_shape, &
       covariance_state_rate, covariance_inclusion_probability, &
       covariance_selectable, initial_log_variance, initial_log_variance_level, &
       initial_log_variance_state_variance, log_variance_state_shape, &
       log_variance_state_rate, log_variance_initial_prior_mean, &
       log_variance_initial_prior_precision, log_variance_offset, &
       volatility_method, structural_design, initial_structural_state, &
       structural_inclusion_probability, structural_selectable) result(out)
      !! Run the time-varying loading and cointegration core of bvectvpalg.
      real(dp), intent(in) :: y(:, :) !! Differenced endogenous variables by observation.
      real(dp), intent(in) :: w(:, :) !! Regressors entering the cointegration term.
      real(dp), intent(in) :: initial_alpha(:, :) !! Initial loading matrix.
      real(dp), intent(in) :: initial_beta(:, :) !! Initial cointegration matrix.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial innovation covariance.
      integer, intent(in) :: iterations !! Number of retained Gibbs draws.
      integer, intent(in) :: burnin !! Number of discarded Gibbs draws.
      real(dp), intent(in) :: initial_state_prior_mean(:) !! Prior mean of the pre-sample alpha and gamma state.
      real(dp), intent(in) :: initial_state_prior_precision(:, :) !! Prior precision of the pre-sample alpha and gamma state.
      real(dp), intent(in) :: state_shape(:) !! Alpha and gamma innovation precision shapes.
      real(dp), intent(in) :: state_rate(:) !! Alpha and gamma innovation precision rates.
      real(dp), intent(in) :: beta_persistence !! AR persistence of cointegration states.
      real(dp), intent(in) :: beta_state_precision(:, :) !! Cointegration-state innovation precision.
      real(dp), intent(in) :: beta_initial_prior_mean(:) !! Prior mean of the pre-sample beta state.
      real(dp), intent(in) :: beta_initial_prior_precision(:, :) !! Prior precision of the pre-sample beta state.
      real(dp), intent(in), optional :: covariance_prior_scale(:, :) !! Inverse-Wishart innovation scale.
      real(dp), intent(in), optional :: covariance_prior_df !! Inverse-Wishart prior degrees of freedom.
      real(dp), intent(in), optional :: measurement_shape(:) !! Diagonal innovation-precision prior shapes.
      real(dp), intent(in), optional :: measurement_rate(:) !! Diagonal innovation-precision prior rates.
      real(dp), intent(in), optional :: x(:, :) !! Regressors outside the cointegration term.
      real(dp), intent(in), optional :: initial_gamma(:, :) !! Initial unrestricted coefficient matrix.
      real(dp), intent(in), optional :: inclusion_probability(:) !! BVS probabilities for unrestricted trajectories.
      integer, intent(in), optional :: selectable(:) !! One-based unrestricted BVS positions.
      real(dp), intent(in), optional :: initial_covariance_state(:) !! Initial packed covariance state.
      real(dp), intent(in), optional :: covariance_state_initial_prior_mean(:) !! Prior mean of the pre-sample covariance state.
      real(dp), intent(in), optional :: &
         covariance_state_initial_prior_precision(:, :) !! Prior precision of the pre-sample covariance state.
      real(dp), intent(in), optional :: covariance_state_shape(:) !! Covariance-state precision prior shapes.
      real(dp), intent(in), optional :: covariance_state_rate(:) !! Covariance-state precision prior rates.
      real(dp), intent(in), optional :: covariance_inclusion_probability(:) !! Covariance-state BVS probabilities.
      integer, intent(in), optional :: covariance_selectable(:) !! One-based covariance-state BVS positions.
      real(dp), intent(in), optional :: initial_log_variance(:, :) !! Initial log-variance paths.
      real(dp), intent(in), optional :: initial_log_variance_level(:) !! Initial pre-sample log variances.
      real(dp), intent(in), optional :: initial_log_variance_state_variance(:) !! Initial log-variance innovation variances.
      real(dp), intent(in), optional :: log_variance_state_shape(:) !! Log-variance innovation precision shapes.
      real(dp), intent(in), optional :: log_variance_state_rate(:) !! Log-variance innovation precision rates.
      real(dp), intent(in), optional :: log_variance_initial_prior_mean(:) !! Prior mean of pre-sample log variances.
      real(dp), intent(in), optional :: &
         log_variance_initial_prior_precision(:, :) !! Prior precision of pre-sample log variances.
      real(dp), intent(in), optional :: log_variance_offset(:) !! Positive log-square transformation offsets.
      character(len=*), intent(in), optional :: volatility_method !! Log-variance mixture method: ksc or ocsn.
      real(dp), intent(in), optional :: structural_design(:, :) !! Observation-major contemporaneous SUR design.
      real(dp), intent(in), optional :: initial_structural_state(:) !! Initial structural coefficient state.
      real(dp), intent(in), optional :: structural_inclusion_probability(:) !! Structural trajectory BVS probabilities.
      integer, intent(in), optional :: structural_selectable(:) !! One-based structural BVS positions.
      type(bvartools_tvp_bvec_draws_t) :: out
      type(bvartools_normal_posterior_t) :: posterior, initial_posterior
      type(bvartools_gamma_posterior_t) :: gamma_posterior
      type(bvartools_bvs_result_t) :: selection
      type(bvartools_bvs_result_t) :: covariance_selection
      type(bvartools_stochastic_volatility_t) :: volatility_draw
      real(dp), allocatable :: alpha(:, :, :), beta(:, :, :), gamma(:, :, :)
      real(dp), allocatable :: pi(:, :, :), state(:, :), beta_state(:, :)
      real(dp), allocatable :: state_initial(:), beta_initial(:), design(:, :)
      real(dp), allocatable :: beta_design(:, :), y_adjusted(:, :), residual(:, :)
      real(dp), allocatable :: state_precision(:, :, :), beta_precision(:, :, :)
      real(dp), allocatable :: measurement_precision(:, :, :), covariance(:, :)
      real(dp), allocatable :: precision(:, :), state_precision_diagonal(:)
      real(dp), allocatable :: normal(:), beta_normal(:), state_vector(:)
      real(dp), allocatable :: beta_vector(:), posterior_scale(:, :)
      real(dp), allocatable :: measurement_diagonal(:)
      real(dp), allocatable :: selection_design(:, :), full_probability(:)
      real(dp), allocatable :: proposal_uniforms(:), acceptance_uniforms(:)
      logical, allocatable :: current_inclusion(:)
      integer, allocatable :: selection_positions(:)
      real(dp), allocatable :: covariance_state(:, :), covariance_state_vector(:)
      real(dp), allocatable :: covariance_state_initial(:)
      real(dp), allocatable :: covariance_state_precision(:, :, :)
      real(dp), allocatable :: covariance_state_precision_diagonal(:)
      real(dp), allocatable :: omega_precision(:, :, :), transformed_residual(:, :)
      real(dp), allocatable :: covariance_normal(:), transform_blocks(:, :)
      real(dp), allocatable :: transform(:, :), full_precision(:, :)
      real(dp), allocatable :: inverse_precision(:, :)
      real(dp), allocatable :: current_time_covariance(:, :, :)
      real(dp), allocatable :: covariance_proposal_uniforms(:)
      real(dp), allocatable :: covariance_acceptance_uniforms(:)
      logical, allocatable :: current_covariance_inclusion(:)
      integer, allocatable :: covariance_selection_positions(:)
      real(dp), allocatable :: log_variance(:, :), log_variance_initial(:)
      real(dp), allocatable :: log_variance_state_variance(:)
      real(dp), allocatable :: log_variance_state_precision(:, :)
      real(dp), allocatable :: log_variance_normal(:)
      real(dp) :: posterior_df
      integer :: k, observations, m, rank, regressors, alpha_count, gamma_count
      integer :: states, beta_count, total, draw, retained, time, item, equation
      integer :: factor, regressor, column, first, last, status
      integer :: covariance_count, structural_count, unrestricted_count
      logical :: full_covariance, use_bvs, use_covariance, use_covariance_bvs
      logical :: use_gamma_bvs, use_structural_bvs
      logical :: use_sv, use_ocsn

      k = size(y, 1)
      observations = size(y, 2)
      m = size(w, 1)
      rank = size(initial_beta, 2)
      regressors = 0
      if (present(x)) regressors = size(x, 1)
      structural_count = 0
      if (present(structural_design)) structural_count = size(structural_design, 2)
      alpha_count = k*rank
      gamma_count = k*regressors
      unrestricted_count = gamma_count + structural_count
      states = alpha_count + unrestricted_count
      beta_count = m*rank
      total = iterations + burnin
      full_covariance = present(covariance_prior_scale)
      use_gamma_bvs = present(inclusion_probability)
      use_structural_bvs = present(structural_inclusion_probability)
      use_bvs = use_gamma_bvs .or. use_structural_bvs
      covariance_count = k*(k - 1)/2
      use_covariance = present(initial_covariance_state)
      use_covariance_bvs = present(covariance_inclusion_probability)
      use_sv = present(initial_log_variance)
      use_ocsn = .false.
      if (present(volatility_method)) then
         select case (trim(volatility_method))
         case ('ksc', 'KSC')
         case ('ocsn', 'OCSN')
            use_ocsn = .true.
         case default
            out%info = 35
            return
         end select
      end if
      if (k < 1 .or. observations < 1 .or. m < 1 .or. rank < 1 .or. &
         size(w, 2) /= observations .or. &
         any(shape(initial_alpha) /= [k, rank]) .or. &
         size(initial_beta, 1) /= m .or. &
         any(shape(initial_covariance) /= [k, k]) .or. &
         size(initial_state_prior_mean) /= states .or. &
         any(shape(initial_state_prior_precision) /= [states, states]) .or. &
         size(state_shape) /= states .or. size(state_rate) /= states .or. &
         any(state_shape <= 0.0_dp) .or. any(state_rate <= 0.0_dp) .or. &
         any(shape(beta_state_precision) /= [beta_count, beta_count]) .or. &
         size(beta_initial_prior_mean) /= beta_count .or. &
         any(shape(beta_initial_prior_precision) /= [beta_count, beta_count]) .or. &
         abs(beta_persistence) >= 1.0_dp .or. iterations < 1 .or. burnin < 0) then
         out%info = 1
         return
      end if
      if (present(x)) then
         if (size(x, 2) /= observations .or. .not. present(initial_gamma)) then
            out%info = 2
            return
         end if
         if (any(shape(initial_gamma) /= [k, regressors])) then
            out%info = 3
            return
         end if
      else if (present(initial_gamma)) then
         out%info = 4
         return
      end if
      if (present(structural_design)) then
         if (size(structural_design, 1) /= k*observations .or. &
            structural_count < 1 .or. .not. present(initial_structural_state)) then
            out%info = 36
            return
         end if
         if (size(initial_structural_state) /= structural_count) then
            out%info = 37
            return
         end if
      else if (present(initial_structural_state)) then
         out%info = 38
         return
      end if
      if (use_gamma_bvs) then
         if (gamma_count < 1 .or. size(inclusion_probability) /= gamma_count .or. &
            any(inclusion_probability < 0.0_dp) .or. &
            any(inclusion_probability > 1.0_dp)) then
            out%info = 15
            return
         end if
         if (present(selectable)) then
            if (any(selectable < 1) .or. any(selectable > gamma_count)) then
               out%info = 16
               return
            end if
         end if
      else if (present(selectable)) then
         out%info = 17
         return
      end if
      if (use_structural_bvs) then
         if (structural_count < 1 .or. &
            size(structural_inclusion_probability) /= structural_count .or. &
            any(structural_inclusion_probability < 0.0_dp) .or. &
            any(structural_inclusion_probability > 1.0_dp)) then
            out%info = 39
            return
         end if
         if (present(structural_selectable)) then
            if (any(structural_selectable < 1) .or. &
               any(structural_selectable > structural_count)) then
               out%info = 40
               return
            end if
         end if
      else if (present(structural_selectable)) then
         out%info = 41
         return
      end if
      if (use_covariance) then
         if (k < 2 .or. full_covariance .or. use_sv .or. &
            .not. present(covariance_state_initial_prior_mean) .or. &
            .not. present(covariance_state_initial_prior_precision) .or. &
            .not. present(covariance_state_shape) .or. &
            .not. present(covariance_state_rate)) then
            out%info = 19
            return
         end if
         if (size(initial_covariance_state) /= covariance_count .or. &
            size(covariance_state_initial_prior_mean) /= covariance_count .or. &
            any(shape(covariance_state_initial_prior_precision) /= &
            [covariance_count, covariance_count]) .or. &
            size(covariance_state_shape) /= covariance_count .or. &
            size(covariance_state_rate) /= covariance_count .or. &
            any(covariance_state_shape <= 0.0_dp) .or. &
            any(covariance_state_rate <= 0.0_dp)) then
            out%info = 20
            return
         end if
      else if (present(covariance_state_initial_prior_mean) .or. &
         present(covariance_state_initial_prior_precision) .or. &
         present(covariance_state_shape) .or. present(covariance_state_rate)) then
         out%info = 21
         return
      end if
      if (use_covariance_bvs) then
         if (.not. use_covariance .or. &
            size(covariance_inclusion_probability) /= covariance_count .or. &
            any(covariance_inclusion_probability < 0.0_dp) .or. &
            any(covariance_inclusion_probability > 1.0_dp)) then
            out%info = 22
            return
         end if
         if (present(covariance_selectable)) then
            if (any(covariance_selectable < 1) .or. &
               any(covariance_selectable > covariance_count)) then
               out%info = 23
               return
            end if
         end if
      else if (present(covariance_selectable)) then
         out%info = 24
         return
      end if
      if (use_sv) then
         if (full_covariance .or. use_covariance .or. &
            .not. present(initial_log_variance_level) .or. &
            .not. present(initial_log_variance_state_variance) .or. &
            .not. present(log_variance_state_shape) .or. &
            .not. present(log_variance_state_rate) .or. &
            .not. present(log_variance_initial_prior_mean) .or. &
            .not. present(log_variance_initial_prior_precision) .or. &
            .not. present(log_variance_offset)) then
            out%info = 30
            return
         end if
         if (any(shape(initial_log_variance) /= [k, observations]) .or. &
            size(initial_log_variance_level) /= k .or. &
            size(initial_log_variance_state_variance) /= k .or. &
            size(log_variance_state_shape) /= k .or. &
            size(log_variance_state_rate) /= k .or. &
            size(log_variance_initial_prior_mean) /= k .or. &
            any(shape(log_variance_initial_prior_precision) /= [k, k]) .or. &
            size(log_variance_offset) /= k .or. &
            any(initial_log_variance_state_variance <= 0.0_dp) .or. &
            any(log_variance_state_shape <= 0.0_dp) .or. &
            any(log_variance_state_rate <= 0.0_dp) .or. &
            any(log_variance_offset <= 0.0_dp)) then
            out%info = 31
            return
         end if
      else if (present(initial_log_variance_level) .or. &
         present(initial_log_variance_state_variance) .or. &
         present(log_variance_state_shape) .or. present(log_variance_state_rate) .or. &
         present(log_variance_initial_prior_mean) .or. &
         present(log_variance_initial_prior_precision) .or. &
         present(log_variance_offset) .or. present(volatility_method)) then
         out%info = 32
         return
      end if
      if (full_covariance) then
         if (.not. present(covariance_prior_df)) then
            out%info = 5
            return
         end if
         if (any(shape(covariance_prior_scale) /= [k, k]) .or. &
            covariance_prior_df <= real(k - 1, dp)) then
            out%info = 6
            return
         end if
      else if (.not. use_sv) then
         if (.not. present(measurement_shape) .or. &
            .not. present(measurement_rate)) then
            out%info = 7
            return
         end if
         if (size(measurement_shape) /= k .or. size(measurement_rate) /= k .or. &
            any(measurement_shape <= 0.0_dp) .or. any(measurement_rate <= 0.0_dp)) then
            out%info = 8
            return
         end if
      end if

      covariance = initial_covariance
      call invert_matrix(covariance, precision, status)
      if (status /= 0) then
         out%info = 9
         return
      end if
      allocate(alpha(k, rank, observations), beta(m, rank, observations), &
         gamma(k, regressors, observations), pi(k, m, observations), &
         state(states, observations), beta_state(beta_count, observations), &
         state_initial(states), beta_initial(beta_count), design(k*observations, states), &
         beta_design(k*observations, beta_count), y_adjusted(k, observations), &
         residual(k, observations), state_precision(states, states, observations), &
         beta_precision(beta_count, beta_count, observations), &
         measurement_precision(k, k, observations), &
         state_precision_diagonal(states), normal(states*observations), &
         beta_normal(beta_count*observations), &
         out%alpha(alpha_count*observations, iterations), &
         out%beta(beta_count*observations, iterations), &
          out%pi(k*m*observations, iterations), &
          out%gamma(gamma_count*observations, iterations), &
          out%structural(structural_count*observations, iterations), &
          out%state_variance(states, iterations), out%covariance(k, k, iterations))
      allocate(out%time_covariance(k, k, observations, iterations), &
         current_time_covariance(k, k, observations))
      if (use_bvs) then
         allocate(current_inclusion(states), full_probability(states))
         current_inclusion = .true.
         full_probability = 0.5_dp
         item = 0
         if (use_gamma_bvs) then
            full_probability(alpha_count + 1:alpha_count + gamma_count) = &
               inclusion_probability
            allocate(out%included(gamma_count, iterations))
            if (present(selectable)) then
               item = size(selectable)
            else
               item = gamma_count
            end if
         end if
         if (use_structural_bvs) then
            full_probability(alpha_count + gamma_count + 1:) = &
               structural_inclusion_probability
            allocate(out%structural_included(structural_count, iterations))
            if (present(structural_selectable)) then
               item = item + size(structural_selectable)
            else
               item = item + structural_count
            end if
         end if
         allocate(selection_positions(item), proposal_uniforms(item), &
            acceptance_uniforms(item))
         item = 0
         if (use_gamma_bvs) then
            if (present(selectable)) then
               selection_positions(1:size(selectable)) = alpha_count + selectable
               item = size(selectable)
            else
               selection_positions(1:gamma_count) = &
                  [(alpha_count + factor, factor=1, gamma_count)]
               item = gamma_count
            end if
         end if
         if (use_structural_bvs) then
            if (present(structural_selectable)) then
               selection_positions(item + 1:) = alpha_count + gamma_count + &
                  structural_selectable
            else
               selection_positions(item + 1:) = &
                  [(alpha_count + gamma_count + factor, &
                  factor=1, structural_count)]
            end if
         end if
      end if
      if (use_covariance) then
         allocate(covariance_state(covariance_count, observations), &
            covariance_state_initial(covariance_count), &
            covariance_state_precision(covariance_count, covariance_count, observations), &
            covariance_state_precision_diagonal(covariance_count), &
            omega_precision(k, k, observations), &
            transformed_residual(k, observations), &
            covariance_normal(covariance_count*observations), &
            out%covariance_state(covariance_count*observations, iterations), &
            out%covariance_state_variance(covariance_count, iterations), &
            out%covariance_initial_state(covariance_count, iterations))
         covariance_state = spread(initial_covariance_state, 2, observations)
         covariance_state_initial = initial_covariance_state
         covariance_state_precision_diagonal = 1.0_dp/covariance_state_rate
         covariance_state_precision = 0.0_dp
         omega_precision = 0.0_dp
         do time = 1, observations
            do item = 1, covariance_count
               covariance_state_precision(item, item, time) = &
                  covariance_state_precision_diagonal(item)
            end do
            do item = 1, k
               omega_precision(item, item, time) = precision(item, item)
            end do
         end do
         if (use_covariance_bvs) then
            allocate(current_covariance_inclusion(covariance_count), &
               out%covariance_included(covariance_count, iterations))
            current_covariance_inclusion = .true.
            if (present(covariance_selectable)) then
               allocate(covariance_selection_positions(size(covariance_selectable)), &
                  covariance_proposal_uniforms(size(covariance_selectable)), &
                  covariance_acceptance_uniforms(size(covariance_selectable)))
               covariance_selection_positions = covariance_selectable
            else
               allocate(covariance_selection_positions(covariance_count), &
                  covariance_proposal_uniforms(covariance_count), &
                  covariance_acceptance_uniforms(covariance_count))
               covariance_selection_positions = [(item, item=1, covariance_count)]
            end if
         end if
      end if
      if (use_sv) then
         log_variance = initial_log_variance
         log_variance_initial = initial_log_variance_level
         log_variance_state_variance = initial_log_variance_state_variance
         allocate(log_variance_state_precision(k, k), log_variance_normal(k), &
            out%log_variance(k*observations, iterations), &
            out%log_variance_state_variance(k, iterations), &
            out%log_variance_initial(k, iterations))
         log_variance_state_precision = 0.0_dp
      end if
      alpha = spread(initial_alpha, 3, observations)
      beta = spread(initial_beta, 3, observations)
      gamma = 0.0_dp
      if (present(initial_gamma)) gamma = spread(initial_gamma, 3, observations)
      state_initial(1:alpha_count) = reshape(initial_alpha, [alpha_count])
      if (present(initial_gamma)) then
         state_initial(alpha_count + 1:alpha_count + gamma_count) = &
            reshape(initial_gamma, [gamma_count])
      end if
      if (present(initial_structural_state)) then
         state_initial(alpha_count + gamma_count + 1:) = initial_structural_state
      end if
      beta_initial = reshape(initial_beta, [beta_count])
      state_precision_diagonal = 1.0_dp/state_rate
      state_precision = 0.0_dp
      beta_precision = 0.0_dp
      if (use_covariance) then
         transform_blocks = bvartools_covar_vector_to_matrix( &
            reshape(covariance_state, [covariance_count*observations]), k, observations)
      end if
      do time = 1, observations
         do item = 1, states
            state_precision(item, item, time) = state_precision_diagonal(item)
         end do
         beta_precision(:, :, time) = beta_state_precision
         if (use_covariance) then
            first = (time - 1)*k + 1
            last = time*k
            transform = transform_blocks(first:last, first:last)
            measurement_precision(:, :, time) = matmul(transpose(transform), &
               matmul(omega_precision(:, :, time), transform))
            call invert_matrix(measurement_precision(:, :, time), inverse_precision, &
               status)
            if (status /= 0) then
               out%info = 25
               return
            end if
            current_time_covariance(:, :, time) = inverse_precision
         else if (use_sv) then
            measurement_precision(:, :, time) = 0.0_dp
            current_time_covariance(:, :, time) = 0.0_dp
            do item = 1, k
               measurement_precision(item, item, time) = &
                  exp(-log_variance(item, time))
               current_time_covariance(item, item, time) = &
                  exp(log_variance(item, time))
            end do
         else
            measurement_precision(:, :, time) = precision
            current_time_covariance(:, :, time) = covariance
         end if
      end do
      posterior_df = 0.0_dp
      if (full_covariance) posterior_df = covariance_prior_df + real(observations, dp)
      retained = 0
      do draw = 1, total
         design = 0.0_dp
         do time = 1, observations
            first = (time - 1)*k + 1
            last = time*k
            do factor = 1, rank
               do equation = 1, k
                  column = (factor - 1)*k + equation
                  design(first + equation - 1, column) = &
                     dot_product(beta(:, factor, time), w(:, time))
               end do
            end do
            do regressor = 1, regressors
               do equation = 1, k
                  column = alpha_count + (regressor - 1)*k + equation
                  design(first + equation - 1, column) = x(regressor, time)
               end do
            end do
         end do
         if (present(structural_design)) then
            design(:, alpha_count + gamma_count + 1:) = structural_design
         end if
         selection_design = design
         if (use_bvs) then
            do item = alpha_count + 1, states
               if (.not. current_inclusion(item)) design(:, item) = 0.0_dp
            end do
         end if
         posterior = bvartools_random_walk_posterior(y, design, measurement_precision, &
            state_precision, state_initial)
         if (posterior%info /= 0) then
            out%info = 10
            return
         end if
         do item = 1, size(normal)
            normal(item) = random_standard_normal()
         end do
         state_vector = bvartools_normal_draw(posterior, normal)
         state = reshape(state_vector, [states, observations])
         if (use_bvs) then
            do item = 1, size(proposal_uniforms)
               proposal_uniforms(item) = random_uniform()
               acceptance_uniforms(item) = random_uniform()
            end do
            selection = bvartools_tvp_bvs(y, selection_design, state, &
               current_inclusion, measurement_precision, full_probability, &
               proposal_uniforms, acceptance_uniforms, selection_positions)
            if (selection%info /= 0) then
               out%info = 18
               return
            end if
            current_inclusion = selection%included
            do item = alpha_count + 1, states
               if (.not. current_inclusion(item)) state(item, :) = 0.0_dp
            end do
            state_vector = reshape(state, [states*observations])
         end if
         do time = 1, observations
            alpha(:, :, time) = reshape(state(1:alpha_count, time), [k, rank])
            if (regressors > 0) then
               gamma(:, :, time) = reshape(state(alpha_count + 1: &
                  alpha_count + gamma_count, time), &
                  [k, regressors])
               y_adjusted(:, time) = y(:, time) - matmul(gamma(:, :, time), x(:, time))
            else
               y_adjusted(:, time) = y(:, time)
            end if
            if (structural_count > 0) then
               first = (time - 1)*k + 1
               last = time*k
               y_adjusted(:, time) = y_adjusted(:, time) - &
                  matmul(structural_design(first:last, :), &
                  state(alpha_count + gamma_count + 1:, time))
            end if
         end do
         beta_design = 0.0_dp
         do time = 1, observations
            first = (time - 1)*k + 1
            do factor = 1, rank
               do regressor = 1, m
                  column = (factor - 1)*m + regressor
                  beta_design(first:first + k - 1, column) = &
                     alpha(:, factor, time)*w(regressor, time)
               end do
            end do
         end do
         posterior = bvartools_random_walk_posterior(y_adjusted, beta_design, &
            measurement_precision, beta_precision, beta_initial, beta_persistence)
         if (posterior%info /= 0) then
            out%info = 11
            return
         end if
         do item = 1, size(beta_normal)
            beta_normal(item) = random_standard_normal()
         end do
         beta_vector = bvartools_normal_draw(posterior, beta_normal)
         beta_state = reshape(beta_vector, [beta_count, observations])
         do time = 1, observations
            beta(:, :, time) = reshape(beta_state(:, time), [m, rank])
            pi(:, :, time) = matmul(alpha(:, :, time), transpose(beta(:, :, time)))
            residual(:, time) = y_adjusted(:, time) - matmul(pi(:, :, time), w(:, time))
         end do
         if (use_covariance) then
            if (use_covariance_bvs) then
               posterior = bvartools_covar_tvp_posterior(residual, omega_precision, &
                  covariance_state_precision, covariance_state_initial, &
                  current_covariance_inclusion)
            else
               posterior = bvartools_covar_tvp_posterior(residual, omega_precision, &
                  covariance_state_precision, covariance_state_initial)
            end if
            if (posterior%info /= 0) then
               out%info = 26
               return
            end if
            do item = 1, size(covariance_normal)
               covariance_normal(item) = random_standard_normal()
            end do
            covariance_state_vector = bvartools_normal_draw(posterior, covariance_normal)
            covariance_state = reshape(covariance_state_vector, &
               [covariance_count, observations])
            if (use_covariance_bvs) then
               do item = 1, size(covariance_proposal_uniforms)
                  covariance_proposal_uniforms(item) = random_uniform()
                  covariance_acceptance_uniforms(item) = random_uniform()
               end do
               covariance_selection = bvartools_tvp_covariance_bvs(residual, &
                  covariance_state, current_covariance_inclusion, omega_precision, &
                  covariance_inclusion_probability, covariance_proposal_uniforms, &
                  covariance_acceptance_uniforms, covariance_selection_positions)
               if (covariance_selection%info /= 0) then
                  out%info = 27
                  return
               end if
               current_covariance_inclusion = covariance_selection%included
               do item = 1, covariance_count
                  if (.not. current_covariance_inclusion(item)) then
                     covariance_state(item, :) = 0.0_dp
                  end if
               end do
               covariance_state_vector = reshape(covariance_state, &
                  [covariance_count*observations])
            end if
            transform_blocks = bvartools_covar_vector_to_matrix( &
               covariance_state_vector, k, observations)
            do time = 1, observations
               first = (time - 1)*k + 1
               last = time*k
               transform = transform_blocks(first:last, first:last)
               transformed_residual(:, time) = matmul(transform, residual(:, time))
            end do
            gamma_posterior = bvartools_measurement_variance_posterior( &
               transformed_residual, measurement_shape, measurement_rate)
            measurement_diagonal = bvartools_gamma_precision_draw(gamma_posterior)
            gamma_posterior = bvartools_state_variance_posterior(covariance_state, &
               covariance_state_initial, covariance_state_shape, covariance_state_rate)
            covariance_state_precision_diagonal = &
               bvartools_gamma_precision_draw(gamma_posterior)
            covariance_state_precision(:, :, 1) = 0.0_dp
            do item = 1, covariance_count
               covariance_state_precision(item, item, 1) = &
                  covariance_state_precision_diagonal(item)
            end do
            initial_posterior = bvartools_initial_state_posterior( &
               covariance_state(:, 1), covariance_state_precision(:, :, 1), &
               covariance_state_initial_prior_mean, &
               covariance_state_initial_prior_precision)
            if (initial_posterior%info /= 0) then
               out%info = 28
               return
            end if
            do item = 1, covariance_count
               covariance_normal(item) = random_standard_normal()
            end do
            covariance_state_initial = bvartools_normal_draw(initial_posterior, &
               covariance_normal(1:covariance_count))
            do time = 1, observations
               covariance_state_precision(:, :, time) = 0.0_dp
               omega_precision(:, :, time) = 0.0_dp
               do item = 1, covariance_count
                  covariance_state_precision(item, item, time) = &
                     covariance_state_precision_diagonal(item)
               end do
               do item = 1, k
                  omega_precision(item, item, time) = measurement_diagonal(item)
               end do
               first = (time - 1)*k + 1
               last = time*k
               transform = transform_blocks(first:last, first:last)
               full_precision = matmul(transpose(transform), &
                  matmul(omega_precision(:, :, time), transform))
               measurement_precision(:, :, time) = full_precision
               call invert_matrix(full_precision, inverse_precision, status)
               if (status /= 0) then
                  out%info = 29
                  return
               end if
               current_time_covariance(:, :, time) = inverse_precision
            end do
            covariance = current_time_covariance(:, :, observations)
            precision = measurement_precision(:, :, observations)
         else if (use_sv) then
            if (use_ocsn) then
               volatility_draw = bvartools_stochastic_volatility_ocsn(residual, &
                  log_variance, log_variance_state_variance, &
                  log_variance_initial, log_variance_offset)
            else
               volatility_draw = bvartools_stochastic_volatility(residual, &
                  log_variance, log_variance_state_variance, &
                  log_variance_initial, log_variance_offset)
            end if
            if (volatility_draw%info /= 0) then
               out%info = 33
               return
            end if
            log_variance = volatility_draw%log_variance
            gamma_posterior = bvartools_state_variance_posterior(log_variance, &
               log_variance_initial, log_variance_state_shape, log_variance_state_rate)
            measurement_diagonal = bvartools_gamma_precision_draw(gamma_posterior)
            log_variance_state_variance = &
               1.0_dp/max(measurement_diagonal, tiny(1.0_dp))
            log_variance_state_precision = 0.0_dp
            do item = 1, k
               log_variance_state_precision(item, item) = measurement_diagonal(item)
            end do
            initial_posterior = bvartools_initial_state_posterior( &
               log_variance(:, 1), log_variance_state_precision, &
               log_variance_initial_prior_mean, &
               log_variance_initial_prior_precision)
            if (initial_posterior%info /= 0) then
               out%info = 34
               return
            end if
            do item = 1, k
               log_variance_normal(item) = random_standard_normal()
            end do
            log_variance_initial = bvartools_normal_draw(initial_posterior, &
               log_variance_normal)
            covariance = 0.0_dp
            precision = 0.0_dp
            do item = 1, k
               covariance(item, item) = exp(log_variance(item, observations))
               precision(item, item) = exp(-log_variance(item, observations))
            end do
         else if (full_covariance) then
            posterior_scale = covariance_prior_scale + matmul(residual, transpose(residual))
            call inverse_wishart_draw(posterior_scale, posterior_df, covariance, &
               precision, status)
            if (status /= 0) then
               out%info = 12
               return
            end if
         else
            gamma_posterior = bvartools_measurement_variance_posterior(residual, &
               measurement_shape, measurement_rate)
            measurement_diagonal = bvartools_gamma_precision_draw(gamma_posterior)
            covariance = 0.0_dp
            precision = 0.0_dp
            do item = 1, k
               precision(item, item) = measurement_diagonal(item)
               covariance(item, item) = 1.0_dp/measurement_diagonal(item)
            end do
         end if
         gamma_posterior = bvartools_state_variance_posterior(state, state_initial, &
            state_shape, state_rate)
         state_precision_diagonal = bvartools_gamma_precision_draw(gamma_posterior)
         state_precision(:, :, 1) = 0.0_dp
         do item = 1, states
            state_precision(item, item, 1) = state_precision_diagonal(item)
         end do
         initial_posterior = bvartools_initial_state_posterior(state(:, 1), &
            state_precision(:, :, 1), initial_state_prior_mean, &
            initial_state_prior_precision)
         if (initial_posterior%info /= 0) then
            out%info = 13
            return
         end if
         do item = 1, states
            normal(item) = random_standard_normal()
         end do
         state_initial = bvartools_normal_draw(initial_posterior, normal(1:states))
         initial_posterior = bvartools_initial_state_posterior(beta_state(:, 1), &
            beta_state_precision, beta_initial_prior_mean, &
            beta_initial_prior_precision, beta_persistence)
         if (initial_posterior%info /= 0) then
            out%info = 14
            return
         end if
         do item = 1, beta_count
            beta_normal(item) = random_standard_normal()
         end do
         beta_initial = bvartools_normal_draw(initial_posterior, &
            beta_normal(1:beta_count))
         do time = 1, observations
            state_precision(:, :, time) = 0.0_dp
            do item = 1, states
               state_precision(item, item, time) = state_precision_diagonal(item)
            end do
            if (use_sv) then
               measurement_precision(:, :, time) = 0.0_dp
               current_time_covariance(:, :, time) = 0.0_dp
               do item = 1, k
                  measurement_precision(item, item, time) = &
                     exp(-log_variance(item, time))
                  current_time_covariance(item, item, time) = &
                     exp(log_variance(item, time))
               end do
            else if (.not. use_covariance) then
               measurement_precision(:, :, time) = precision
            end if
         end do
         if (draw > burnin) then
            retained = retained + 1
            out%alpha(:, retained) = reshape(alpha, [alpha_count*observations])
            out%beta(:, retained) = reshape(beta, [beta_count*observations])
            out%pi(:, retained) = reshape(pi, [k*m*observations])
            if (gamma_count > 0) then
               out%gamma(:, retained) = reshape(gamma, [gamma_count*observations])
            end if
            if (structural_count > 0) then
               out%structural(:, retained) = reshape( &
                  state(alpha_count + gamma_count + 1:, :), &
                  [structural_count*observations])
            end if
            out%state_variance(:, retained) = &
               1.0_dp/max(state_precision_diagonal, tiny(1.0_dp))
            out%covariance(:, :, retained) = covariance
            if (use_covariance) then
               out%time_covariance(:, :, :, retained) = current_time_covariance
               out%covariance_state(:, retained) = covariance_state_vector
               out%covariance_state_variance(:, retained) = &
                  1.0_dp/max(covariance_state_precision_diagonal, tiny(1.0_dp))
               out%covariance_initial_state(:, retained) = covariance_state_initial
               if (use_covariance_bvs) out%covariance_included(:, retained) = &
                  current_covariance_inclusion
            else if (use_sv) then
               out%time_covariance(:, :, :, retained) = current_time_covariance
               out%log_variance(:, retained) = &
                  reshape(log_variance, [k*observations])
               out%log_variance_state_variance(:, retained) = &
                  log_variance_state_variance
               out%log_variance_initial(:, retained) = log_variance_initial
            else
               do time = 1, observations
                  out%time_covariance(:, :, time, retained) = covariance
               end do
            end if
            if (use_gamma_bvs) out%included(:, retained) = current_inclusion( &
               alpha_count + 1:alpha_count + gamma_count)
            if (use_structural_bvs) out%structural_included(:, retained) = &
               current_inclusion(alpha_count + gamma_count + 1:)
         end if
      end do
      out%retained_draws = retained
   end function bvartools_tvp_bvec_gibbs

   function bvartools_tvp_covariance_gibbs(residuals, initial_diagonal_variance, &
      initial_state, iterations, burnin, initial_state_prior_mean, &
      initial_state_prior_precision, state_shape, state_rate, &
      measurement_shape, measurement_rate, inclusion_probability, &
      selectable) result(out)
      !! Sample time-varying lower-triangular covariance states and variances.
      real(dp), intent(in) :: residuals(:, :) !! Reduced-form residuals by observation.
      real(dp), intent(in) :: initial_diagonal_variance(:) !! Initial orthogonal innovation variances.
      real(dp), intent(in) :: initial_state(:) !! Initial packed covariance state.
      integer, intent(in) :: iterations !! Number of retained Gibbs draws.
      integer, intent(in) :: burnin !! Number of discarded Gibbs draws.
      real(dp), intent(in) :: initial_state_prior_mean(:) !! Prior mean of the pre-sample covariance state.
      real(dp), intent(in) :: initial_state_prior_precision(:, :) !! Prior precision of the pre-sample covariance state.
      real(dp), intent(in) :: state_shape(:) !! Covariance-state precision prior shapes.
      real(dp), intent(in) :: state_rate(:) !! Covariance-state precision prior rates.
      real(dp), intent(in) :: measurement_shape(:) !! Orthogonal innovation-precision prior shapes.
      real(dp), intent(in) :: measurement_rate(:) !! Orthogonal innovation-precision prior rates.
      real(dp), intent(in), optional :: inclusion_probability(:) !! BVS trajectory inclusion probabilities.
      integer, intent(in), optional :: selectable(:) !! One-based BVS trajectory positions.
      type(bvartools_tvp_covariance_draws_t) :: out
      type(bvartools_normal_posterior_t) :: posterior, initial_posterior
      type(bvartools_gamma_posterior_t) :: gamma_posterior
      type(bvartools_bvs_result_t) :: selection
      real(dp), allocatable :: state(:, :), state_vector(:), current_initial(:)
      real(dp), allocatable :: state_precision(:, :, :), diagonal_precision(:, :, :)
      real(dp), allocatable :: state_precision_diagonal(:), measurement_precision(:)
      real(dp), allocatable :: normal(:), transform_blocks(:, :), transform(:, :)
      real(dp), allocatable :: transformed_residuals(:, :), full_precision(:, :)
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: proposal_uniforms(:), acceptance_uniforms(:)
      logical, allocatable :: current_inclusion(:)
      integer :: k, observations, count, total, draw, retained, time, item, status
      integer :: first, last
      logical :: use_bvs

      k = size(residuals, 1)
      observations = size(residuals, 2)
      count = k*(k - 1)/2
      total = iterations + burnin
      use_bvs = present(inclusion_probability)
      if (k < 2 .or. observations < 1 .or. size(initial_diagonal_variance) /= k .or. &
         any(initial_diagonal_variance <= 0.0_dp) .or. size(initial_state) /= count .or. &
         size(initial_state_prior_mean) /= count .or. &
         any(shape(initial_state_prior_precision) /= [count, count]) .or. &
         size(state_shape) /= count .or. size(state_rate) /= count .or. &
         any(state_shape <= 0.0_dp) .or. any(state_rate <= 0.0_dp) .or. &
         size(measurement_shape) /= k .or. size(measurement_rate) /= k .or. &
         any(measurement_shape <= 0.0_dp) .or. any(measurement_rate <= 0.0_dp) .or. &
         iterations < 1 .or. burnin < 0) then
         out%info = 1
         return
      end if
      if (use_bvs) then
         if (size(inclusion_probability) /= count .or. &
            any(inclusion_probability < 0.0_dp) .or. &
            any(inclusion_probability > 1.0_dp)) then
            out%info = 5
            return
         end if
         if (present(selectable)) then
            if (any(selectable < 1) .or. any(selectable > count)) then
               out%info = 6
               return
            end if
         end if
      else if (present(selectable)) then
         out%info = 7
         return
      end if

      current_initial = initial_state
      allocate(state(count, observations), state_precision(count, count, observations), &
         diagonal_precision(k, k, observations), state_precision_diagonal(count), &
         measurement_precision(k), normal(count*observations), &
         transformed_residuals(k, observations), &
         out%states(count*observations, iterations), &
         out%state_variance(count, iterations), &
         out%initial_state(count, iterations), &
         out%diagonal_variance(k, iterations), &
         out%covariance(k, k, observations, iterations))
      if (use_bvs) then
         allocate(current_inclusion(count), out%included(count, iterations))
         current_inclusion = .true.
         if (present(selectable)) then
            allocate(proposal_uniforms(size(selectable)), &
               acceptance_uniforms(size(selectable)))
         else
            allocate(proposal_uniforms(count), acceptance_uniforms(count))
         end if
      end if
      state_precision_diagonal = 1.0_dp/state_rate
      measurement_precision = 1.0_dp/initial_diagonal_variance
      state_precision = 0.0_dp
      diagonal_precision = 0.0_dp
      do time = 1, observations
         do item = 1, count
            state_precision(item, item, time) = state_precision_diagonal(item)
         end do
         do item = 1, k
            diagonal_precision(item, item, time) = measurement_precision(item)
         end do
      end do
      retained = 0
      do draw = 1, total
         if (use_bvs) then
            posterior = bvartools_covar_tvp_posterior(residuals, diagonal_precision, &
               state_precision, current_initial, current_inclusion)
         else
            posterior = bvartools_covar_tvp_posterior(residuals, diagonal_precision, &
               state_precision, current_initial)
         end if
         if (posterior%info /= 0) then
            out%info = 2
            return
         end if
         do item = 1, size(normal)
            normal(item) = random_standard_normal()
         end do
         state_vector = bvartools_normal_draw(posterior, normal)
         state = reshape(state_vector, [count, observations])
         if (use_bvs) then
            do item = 1, size(proposal_uniforms)
               proposal_uniforms(item) = random_uniform()
               acceptance_uniforms(item) = random_uniform()
            end do
            selection = bvartools_tvp_covariance_bvs(residuals, state, &
               current_inclusion, diagonal_precision, inclusion_probability, &
               proposal_uniforms, acceptance_uniforms, selectable)
            if (selection%info /= 0) then
               out%info = 8
               return
            end if
            current_inclusion = selection%included
            do item = 1, count
               if (.not. current_inclusion(item)) state(item, :) = 0.0_dp
            end do
            state_vector = reshape(state, [count*observations])
         end if
         transform_blocks = bvartools_covar_vector_to_matrix(state_vector, k, observations)
         do time = 1, observations
            first = (time - 1)*k + 1
            last = time*k
            transform = transform_blocks(first:last, first:last)
            transformed_residuals(:, time) = matmul(transform, residuals(:, time))
         end do
         gamma_posterior = bvartools_measurement_variance_posterior( &
            transformed_residuals, measurement_shape, measurement_rate)
         measurement_precision = bvartools_gamma_precision_draw(gamma_posterior)
         gamma_posterior = bvartools_state_variance_posterior(state, current_initial, &
            state_shape, state_rate)
         state_precision_diagonal = bvartools_gamma_precision_draw(gamma_posterior)
         state_precision(:, :, 1) = 0.0_dp
         do item = 1, count
            state_precision(item, item, 1) = state_precision_diagonal(item)
         end do
         initial_posterior = bvartools_initial_state_posterior(state(:, 1), &
            state_precision(:, :, 1), initial_state_prior_mean, &
            initial_state_prior_precision)
         if (initial_posterior%info /= 0) then
            out%info = 3
            return
         end if
         do item = 1, count
            normal(item) = random_standard_normal()
         end do
         current_initial = bvartools_normal_draw(initial_posterior, normal(1:count))
         do time = 1, observations
            state_precision(:, :, time) = 0.0_dp
            diagonal_precision(:, :, time) = 0.0_dp
            do item = 1, count
               state_precision(item, item, time) = state_precision_diagonal(item)
            end do
            do item = 1, k
               diagonal_precision(item, item, time) = measurement_precision(item)
            end do
         end do
         if (draw > burnin) then
            retained = retained + 1
            out%states(:, retained) = state_vector
            out%state_variance(:, retained) = &
               1.0_dp/max(state_precision_diagonal, tiny(1.0_dp))
            out%initial_state(:, retained) = current_initial
            out%diagonal_variance(:, retained) = &
               1.0_dp/max(measurement_precision, tiny(1.0_dp))
            if (use_bvs) out%included(:, retained) = current_inclusion
            do time = 1, observations
               first = (time - 1)*k + 1
               last = time*k
               transform = transform_blocks(first:last, first:last)
               full_precision = matmul(transpose(transform), &
                  matmul(diagonal_precision(:, :, time), transform))
               call invert_matrix(full_precision, covariance, status)
               if (status /= 0) then
                  out%info = 4
                  return
               end if
               out%covariance(:, :, time, retained) = covariance
            end do
         end if
      end do
      out%retained_draws = retained
   end function bvartools_tvp_covariance_gibbs

   function bvartools_joint_tvp_bvar_gibbs(y, design, initial_state, &
      initial_diagonal_variance, initial_covariance_state, iterations, burnin, &
      initial_state_prior_mean, initial_state_prior_precision, state_shape, &
      state_rate, covariance_initial_prior_mean, &
      covariance_initial_prior_precision, covariance_state_shape, &
      covariance_state_rate, measurement_shape, measurement_rate, &
      time_varying_covariance, inclusion_probability, selectable, &
      covariance_inclusion_probability, covariance_selectable, covariance_tau0, &
      covariance_tau1, covariance_ssvs_probability) result(out)
      !! Jointly sample TVP-BVAR coefficients and lower-triangular covariance terms.
      real(dp), intent(in) :: y(:, :) !! Endogenous variables by observation.
      real(dp), intent(in) :: design(:, :) !! Observation-major constant SUR design.
      real(dp), intent(in) :: initial_state(:) !! Initial pre-sample coefficient state.
      real(dp), intent(in) :: initial_diagonal_variance(:) !! Initial orthogonal innovation variances.
      real(dp), intent(in) :: initial_covariance_state(:) !! Initial packed covariance coefficients.
      integer, intent(in) :: iterations !! Number of retained Gibbs draws.
      integer, intent(in) :: burnin !! Number of discarded Gibbs draws.
      real(dp), intent(in) :: initial_state_prior_mean(:) !! Prior mean of the coefficient initial state.
      real(dp), intent(in) :: initial_state_prior_precision(:, :) !! Prior precision of the coefficient initial state.
      real(dp), intent(in) :: state_shape(:) !! Coefficient-state precision prior shapes.
      real(dp), intent(in) :: state_rate(:) !! Coefficient-state precision prior rates.
      real(dp), intent(in) :: covariance_initial_prior_mean(:) !! Prior mean of covariance coefficients or their initial state.
      real(dp), intent(in) :: covariance_initial_prior_precision(:, :) !! Covariance initial-state prior precision.
      real(dp), intent(in) :: covariance_state_shape(:) !! Covariance-state precision prior shapes.
      real(dp), intent(in) :: covariance_state_rate(:) !! Covariance-state precision prior rates.
      real(dp), intent(in) :: measurement_shape(:) !! Orthogonal innovation-precision prior shapes.
      real(dp), intent(in) :: measurement_rate(:) !! Orthogonal innovation-precision prior rates.
      logical, intent(in), optional :: time_varying_covariance !! Use random-walk rather than constant covariance coefficients.
      real(dp), intent(in), optional :: inclusion_probability(:) !! Coefficient trajectory BVS probabilities.
      integer, intent(in), optional :: selectable(:) !! Selectable coefficient trajectory positions.
      real(dp), intent(in), optional :: covariance_inclusion_probability(:) !! Covariance trajectory BVS probabilities.
      integer, intent(in), optional :: covariance_selectable(:) !! Selectable covariance coefficient positions.
      real(dp), intent(in), optional :: covariance_tau0(:) !! Excluded-state covariance SSVS deviations.
      real(dp), intent(in), optional :: covariance_tau1(:) !! Included-state covariance SSVS deviations.
      real(dp), intent(in), optional :: covariance_ssvs_probability(:) !! Covariance SSVS inclusion probabilities.
      type(bvartools_joint_tvp_bvar_draws_t) :: out
      type(bvartools_normal_posterior_t) :: coefficient_posterior
      type(bvartools_normal_posterior_t) :: covariance_posterior
      type(bvartools_normal_posterior_t) :: initial_posterior
      type(bvartools_gamma_posterior_t) :: gamma_posterior
      type(bvartools_bvs_result_t) :: selection
      type(bvartools_ssvs_result_t) :: ssvs_selection
      real(dp), allocatable :: coefficient_state(:, :), coefficient_vector(:)
      real(dp), allocatable :: coefficient_initial(:), coefficient_precision(:, :, :)
      real(dp), allocatable :: coefficient_precision_diagonal(:), coefficient_normal(:)
      real(dp), allocatable :: covariance_state(:, :), covariance_vector(:)
      real(dp), allocatable :: covariance_initial(:), covariance_precision(:, :, :)
      real(dp), allocatable :: covariance_precision_diagonal(:), covariance_normal(:)
      real(dp), allocatable :: covariance_prior_precision(:, :)
      real(dp), allocatable :: diagonal_precision(:, :, :), measurement_precision(:)
      real(dp), allocatable :: full_precision(:, :, :), transform_blocks(:, :)
      real(dp), allocatable :: transform(:, :), residual(:, :), transformed_residual(:, :)
      real(dp), allocatable :: active_design(:, :), proposal_uniforms(:)
      real(dp), allocatable :: acceptance_uniforms(:), covariance_proposal_uniforms(:)
      real(dp), allocatable :: covariance_acceptance_uniforms(:), ssvs_uniforms(:)
      real(dp), allocatable :: covariance_matrix(:, :)
      logical, allocatable :: coefficient_included(:), covariance_included(:)
      logical :: covariance_is_tvp, use_coefficient_bvs, use_covariance_bvs
      logical :: use_covariance_ssvs
      integer :: k, observations, states, covariance_count, total
      integer :: draw, retained, time, item, first, last, status

      k = size(y, 1)
      observations = size(y, 2)
      states = size(initial_state)
      covariance_count = k*(k - 1)/2
      total = iterations + burnin
      covariance_is_tvp = .true.
      if (present(time_varying_covariance)) covariance_is_tvp = &
         time_varying_covariance
      use_coefficient_bvs = present(inclusion_probability)
      use_covariance_bvs = present(covariance_inclusion_probability)
      use_covariance_ssvs = present(covariance_tau0) .or. &
         present(covariance_tau1) .or. present(covariance_ssvs_probability)
      if (k < 2 .or. observations < 1 .or. states < 1 .or. &
         size(design, 1) /= k*observations .or. size(design, 2) /= states .or. &
         size(initial_diagonal_variance) /= k .or. &
         any(initial_diagonal_variance <= 0.0_dp) .or. &
         size(initial_covariance_state) /= covariance_count .or. &
         size(initial_state_prior_mean) /= states .or. &
         any(shape(initial_state_prior_precision) /= [states, states]) .or. &
         size(state_shape) /= states .or. size(state_rate) /= states .or. &
         any(state_shape <= 0.0_dp) .or. any(state_rate <= 0.0_dp) .or. &
         size(covariance_initial_prior_mean) /= covariance_count .or. &
         any(shape(covariance_initial_prior_precision) /= &
         [covariance_count, covariance_count]) .or. &
         size(covariance_state_shape) /= covariance_count .or. &
         size(covariance_state_rate) /= covariance_count .or. &
         any(covariance_state_shape <= 0.0_dp) .or. &
         any(covariance_state_rate <= 0.0_dp) .or. &
         size(measurement_shape) /= k .or. size(measurement_rate) /= k .or. &
         any(measurement_shape <= 0.0_dp) .or. any(measurement_rate <= 0.0_dp) .or. &
         iterations < 1 .or. burnin < 0) then
         out%info = 1
         return
      end if
      if (use_coefficient_bvs) then
         if (size(inclusion_probability) /= states .or. &
            any(inclusion_probability < 0.0_dp) .or. &
            any(inclusion_probability > 1.0_dp)) then
            out%info = 2
            return
         end if
         if (present(selectable)) then
            if (any(selectable < 1) .or. any(selectable > states)) then
               out%info = 3
               return
            end if
         end if
      else if (present(selectable)) then
         out%info = 4
         return
      end if
      if (use_covariance_bvs) then
         if (.not. covariance_is_tvp .or. &
            size(covariance_inclusion_probability) /= covariance_count .or. &
            any(covariance_inclusion_probability < 0.0_dp) .or. &
            any(covariance_inclusion_probability > 1.0_dp)) then
            out%info = 5
            return
         end if
         if (present(covariance_selectable)) then
            if (any(covariance_selectable < 1) .or. &
               any(covariance_selectable > covariance_count)) then
               out%info = 6
               return
            end if
         end if
      else if (present(covariance_selectable) .and. .not. use_covariance_ssvs) then
         out%info = 7
         return
      end if
      if (use_covariance_ssvs) then
         if (covariance_is_tvp .or. .not. present(covariance_tau0) .or. &
            .not. present(covariance_tau1) .or. &
            .not. present(covariance_ssvs_probability)) then
            out%info = 8
            return
         end if
         if (size(covariance_tau0) /= covariance_count .or. &
            size(covariance_tau1) /= covariance_count .or. &
            size(covariance_ssvs_probability) /= covariance_count .or. &
            any(covariance_tau0 <= 0.0_dp) .or. any(covariance_tau1 <= 0.0_dp) .or. &
            any(covariance_ssvs_probability < 0.0_dp) .or. &
            any(covariance_ssvs_probability > 1.0_dp)) then
            out%info = 9
            return
         end if
      end if

      out%time_varying_covariance = covariance_is_tvp
      coefficient_initial = initial_state
      covariance_initial = initial_covariance_state
      allocate(coefficient_state(states, observations), &
         coefficient_precision(states, states, observations), &
         coefficient_precision_diagonal(states), &
         coefficient_normal(states*observations), &
         covariance_state(covariance_count, observations), &
         covariance_precision(covariance_count, covariance_count, observations), &
         covariance_precision_diagonal(covariance_count), &
         covariance_normal(covariance_count*observations), &
         covariance_prior_precision(covariance_count, covariance_count), &
         diagonal_precision(k, k, observations), measurement_precision(k), &
         full_precision(k, k, observations), residual(k, observations), &
         transformed_residual(k, observations), &
         out%coefficient_states(states*observations, iterations), &
         out%coefficient_state_variance(states, iterations), &
         out%coefficient_initial_state(states, iterations), &
         out%covariance_states(covariance_count*observations, iterations), &
         out%covariance_initial_state(covariance_count, iterations), &
         out%diagonal_variance(k, iterations), &
         out%covariance(k, k, observations, iterations))
      if (covariance_is_tvp) then
         allocate(out%covariance_state_variance(covariance_count, iterations))
      end if
      active_design = design
      covariance_state = spread(initial_covariance_state, 2, observations)
      coefficient_precision_diagonal = 1.0_dp/state_rate
      covariance_precision_diagonal = 1.0_dp/covariance_state_rate
      measurement_precision = 1.0_dp/initial_diagonal_variance
      covariance_prior_precision = covariance_initial_prior_precision
      coefficient_precision = 0.0_dp
      covariance_precision = 0.0_dp
      diagonal_precision = 0.0_dp
      do time = 1, observations
         do item = 1, states
            coefficient_precision(item, item, time) = &
               coefficient_precision_diagonal(item)
         end do
         do item = 1, covariance_count
            covariance_precision(item, item, time) = &
               covariance_precision_diagonal(item)
         end do
         do item = 1, k
            diagonal_precision(item, item, time) = measurement_precision(item)
         end do
      end do
      if (use_coefficient_bvs) then
         allocate(coefficient_included(states), &
            out%coefficient_included(states, iterations))
         coefficient_included = .true.
         if (present(selectable)) then
            allocate(proposal_uniforms(size(selectable)), &
               acceptance_uniforms(size(selectable)))
         else
            allocate(proposal_uniforms(states), acceptance_uniforms(states))
         end if
      end if
      if (use_covariance_bvs .or. use_covariance_ssvs) then
         allocate(covariance_included(covariance_count), &
            out%covariance_included(covariance_count, iterations))
         covariance_included = .true.
      end if
      if (use_covariance_bvs) then
         if (present(covariance_selectable)) then
            allocate(covariance_proposal_uniforms(size(covariance_selectable)), &
               covariance_acceptance_uniforms(size(covariance_selectable)))
         else
            allocate(covariance_proposal_uniforms(covariance_count), &
               covariance_acceptance_uniforms(covariance_count))
         end if
      end if
      if (use_covariance_ssvs) then
         if (present(covariance_selectable)) then
            allocate(ssvs_uniforms(size(covariance_selectable)))
         else
            allocate(ssvs_uniforms(covariance_count))
         end if
      end if
      call update_joint_precision(covariance_state, diagonal_precision, &
         full_precision, status)
      if (status /= 0) then
         out%info = 10
         return
      end if

      retained = 0
      do draw = 1, total
         active_design = design
         if (use_coefficient_bvs) then
            do item = 1, states
               if (.not. coefficient_included(item)) active_design(:, item) = 0.0_dp
            end do
         end if
         coefficient_posterior = bvartools_random_walk_posterior(y, active_design, &
            full_precision, coefficient_precision, coefficient_initial)
         if (coefficient_posterior%info /= 0) then
            out%info = 11
            return
         end if
         do item = 1, size(coefficient_normal)
            coefficient_normal(item) = random_standard_normal()
         end do
         coefficient_vector = bvartools_normal_draw(coefficient_posterior, &
            coefficient_normal)
         coefficient_state = reshape(coefficient_vector, [states, observations])
         if (use_coefficient_bvs) then
            do item = 1, size(proposal_uniforms)
               proposal_uniforms(item) = random_uniform()
               acceptance_uniforms(item) = random_uniform()
            end do
            selection = bvartools_tvp_bvs(y, design, coefficient_state, &
               coefficient_included, full_precision, inclusion_probability, &
               proposal_uniforms, acceptance_uniforms, selectable)
            if (selection%info /= 0) then
               out%info = 12
               return
            end if
            coefficient_included = selection%included
            do item = 1, states
               if (.not. coefficient_included(item)) coefficient_state(item, :) = 0.0_dp
            end do
            coefficient_vector = reshape(coefficient_state, [states*observations])
         end if
         do time = 1, observations
            first = (time - 1)*k + 1
            last = time*k
            residual(:, time) = y(:, time) - &
               matmul(design(first:last, :), coefficient_state(:, time))
         end do

         if (covariance_is_tvp) then
            if (use_covariance_bvs) then
               covariance_posterior = bvartools_covar_tvp_posterior(residual, &
                  diagonal_precision, covariance_precision, covariance_initial, &
                  covariance_included)
            else
               covariance_posterior = bvartools_covar_tvp_posterior(residual, &
                  diagonal_precision, covariance_precision, covariance_initial)
            end if
            if (covariance_posterior%info /= 0) then
               out%info = 13
               return
            end if
            do item = 1, covariance_count*observations
               covariance_normal(item) = random_standard_normal()
            end do
            covariance_vector = bvartools_normal_draw(covariance_posterior, &
               covariance_normal)
            covariance_state = reshape(covariance_vector, &
               [covariance_count, observations])
            if (use_covariance_bvs) then
               do item = 1, size(covariance_proposal_uniforms)
                  covariance_proposal_uniforms(item) = random_uniform()
                  covariance_acceptance_uniforms(item) = random_uniform()
               end do
               selection = bvartools_tvp_covariance_bvs(residual, covariance_state, &
                  covariance_included, diagonal_precision, &
                  covariance_inclusion_probability, covariance_proposal_uniforms, &
                  covariance_acceptance_uniforms, covariance_selectable)
               if (selection%info /= 0) then
                  out%info = 14
                  return
               end if
               covariance_included = selection%included
               do item = 1, covariance_count
                  if (.not. covariance_included(item)) &
                     covariance_state(item, :) = 0.0_dp
               end do
            end if
         else
            covariance_posterior = bvartools_covar_const_posterior(residual, &
               diagonal_precision, covariance_initial_prior_mean, &
               covariance_prior_precision)
            if (covariance_posterior%info /= 0) then
               out%info = 15
               return
            end if
            do item = 1, covariance_count
               covariance_normal(item) = random_standard_normal()
            end do
            covariance_vector = bvartools_normal_draw(covariance_posterior, &
               covariance_normal(1:covariance_count))
            covariance_state = spread(covariance_vector, 2, observations)
            covariance_initial = covariance_vector
            if (use_covariance_ssvs) then
               do item = 1, size(ssvs_uniforms)
                  ssvs_uniforms(item) = random_uniform()
               end do
               ssvs_selection = bvartools_ssvs(covariance_vector, covariance_tau0, &
                  covariance_tau1, covariance_ssvs_probability, ssvs_uniforms, &
                  covariance_selectable)
               if (ssvs_selection%info /= 0) then
                  out%info = 16
                  return
               end if
               covariance_included = ssvs_selection%included
               covariance_prior_precision = covariance_initial_prior_precision
               do item = 1, covariance_count
                  covariance_prior_precision(item, item) = &
                     ssvs_selection%precision(item)
               end do
            end if
         end if

         call transform_joint_residual(residual, covariance_state, transformed_residual)
         gamma_posterior = bvartools_measurement_variance_posterior( &
            transformed_residual, measurement_shape, measurement_rate)
         measurement_precision = bvartools_gamma_precision_draw(gamma_posterior)
         gamma_posterior = bvartools_state_variance_posterior(coefficient_state, &
            coefficient_initial, state_shape, state_rate)
         coefficient_precision_diagonal = &
            bvartools_gamma_precision_draw(gamma_posterior)
         initial_posterior = bvartools_initial_state_posterior( &
            coefficient_state(:, 1), diagonal_matrix(coefficient_precision_diagonal), &
            initial_state_prior_mean, initial_state_prior_precision)
         if (initial_posterior%info /= 0) then
            out%info = 17
            return
         end if
         do item = 1, states
            coefficient_normal(item) = random_standard_normal()
         end do
         coefficient_initial = bvartools_normal_draw(initial_posterior, &
            coefficient_normal(1:states))
         if (covariance_is_tvp) then
            gamma_posterior = bvartools_state_variance_posterior(covariance_state, &
               covariance_initial, covariance_state_shape, covariance_state_rate)
            covariance_precision_diagonal = &
               bvartools_gamma_precision_draw(gamma_posterior)
            initial_posterior = bvartools_initial_state_posterior( &
               covariance_state(:, 1), diagonal_matrix(covariance_precision_diagonal), &
               covariance_initial_prior_mean, covariance_initial_prior_precision)
            if (initial_posterior%info /= 0) then
               out%info = 18
               return
            end if
            do item = 1, covariance_count
               covariance_normal(item) = random_standard_normal()
            end do
            covariance_initial = bvartools_normal_draw(initial_posterior, &
               covariance_normal(1:covariance_count))
         end if
         do time = 1, observations
            coefficient_precision(:, :, time) = &
               diagonal_matrix(coefficient_precision_diagonal)
            covariance_precision(:, :, time) = &
               diagonal_matrix(covariance_precision_diagonal)
            diagonal_precision(:, :, time) = diagonal_matrix(measurement_precision)
         end do
         call update_joint_precision(covariance_state, diagonal_precision, &
            full_precision, status)
         if (status /= 0) then
            out%info = 19
            return
         end if
         if (draw > burnin) then
            retained = retained + 1
            out%coefficient_states(:, retained) = coefficient_vector
            out%coefficient_state_variance(:, retained) = &
               1.0_dp/max(coefficient_precision_diagonal, tiny(1.0_dp))
            out%coefficient_initial_state(:, retained) = coefficient_initial
            out%covariance_states(:, retained) = reshape(covariance_state, &
               [covariance_count*observations])
            out%covariance_initial_state(:, retained) = covariance_initial
            out%diagonal_variance(:, retained) = &
               1.0_dp/max(measurement_precision, tiny(1.0_dp))
            if (covariance_is_tvp) out%covariance_state_variance(:, retained) = &
               1.0_dp/max(covariance_precision_diagonal, tiny(1.0_dp))
            if (use_coefficient_bvs) out%coefficient_included(:, retained) = &
               coefficient_included
            if (use_covariance_bvs .or. use_covariance_ssvs) &
               out%covariance_included(:, retained) = covariance_included
            do time = 1, observations
               call invert_matrix(full_precision(:, :, time), covariance_matrix, status)
               if (status /= 0) then
                  out%info = 20
                  return
               end if
               out%covariance(:, :, time, retained) = covariance_matrix
            end do
         end if
      end do
      out%retained_draws = retained
   end function bvartools_joint_tvp_bvar_gibbs

   function bvartools_bvec_gibbs(y, w, initial_beta, initial_covariance, &
      iterations, burnin, shrinkage, cointegration_precision, &
      loading_precision, covariance_prior_df, x, gamma_prior_mean, &
      gamma_prior_precision, gamma_shape, gamma_rate, tau0, tau1, &
      inclusion_probability, selectable, structural_design, &
      structural_prior_mean, structural_prior_precision, structural_tau0, &
      structural_tau1, structural_inclusion_probability, &
      structural_selectable) result(out)
      !! Run the constant-parameter bvartools Bayesian VECM Gibbs sampler.
      real(dp), intent(in) :: y(:, :) !! Differenced endogenous variables by observation.
      real(dp), intent(in) :: w(:, :) !! Regressors entering the cointegration term.
      real(dp), intent(in) :: initial_beta(:, :) !! Initial normalized cointegration matrix.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial innovation covariance matrix.
      integer, intent(in) :: iterations !! Number of retained Gibbs draws.
      integer, intent(in) :: burnin !! Number of discarded Gibbs draws.
      real(dp), intent(in) :: shrinkage !! Cointegration-space prior shrinkage.
      real(dp), intent(in) :: cointegration_precision(:, :) !! Cointegration-space prior precision.
      real(dp), intent(in) :: loading_precision(:, :) !! Initial loading-space prior precision.
      real(dp), intent(in), optional :: covariance_prior_df !! Inverse-Wishart prior degrees of freedom.
      real(dp), intent(in), optional :: x(:, :) !! Regressors outside the cointegration term.
      real(dp), intent(in), optional :: gamma_prior_mean(:) !! Prior mean of unrestricted coefficients.
      real(dp), intent(in), optional :: gamma_prior_precision(:, :) !! Prior precision of unrestricted coefficients.
      real(dp), intent(in), optional :: gamma_shape(:) !! Diagonal innovation-precision prior shapes.
      real(dp), intent(in), optional :: gamma_rate(:) !! Diagonal innovation-precision prior rates.
      real(dp), intent(in), optional :: tau0(:) !! SSVS excluded-state standard deviations.
      real(dp), intent(in), optional :: tau1(:) !! SSVS included-state standard deviations.
      real(dp), intent(in), optional :: inclusion_probability(:) !! SSVS inclusion probabilities.
      integer, intent(in), optional :: selectable(:) !! SSVS unrestricted-coefficient positions.
      real(dp), intent(in), optional :: structural_design(:, :) !! Contemporaneous structural SUR design.
      real(dp), intent(in), optional :: structural_prior_mean(:) !! Structural coefficient prior mean.
      real(dp), intent(in), optional :: structural_prior_precision(:, :) !! Structural coefficient prior precision.
      real(dp), intent(in), optional :: structural_tau0(:) !! Structural SSVS spike deviations.
      real(dp), intent(in), optional :: structural_tau1(:) !! Structural SSVS slab deviations.
      real(dp), intent(in), optional :: structural_inclusion_probability(:) !! Structural SSVS probabilities.
      integer, intent(in), optional :: structural_selectable(:) !! Selectable structural coefficients.
      type(bvartools_bvec_draws_t) :: out
      type(bvartools_cointegration_draw_t) :: coefficient_draw
      type(bvartools_ssvs_result_t) :: selection
      type(bvartools_ssvs_result_t) :: structural_selection
      type(bvartools_normal_posterior_t) :: structural_posterior
      real(dp), allocatable :: beta(:, :), covariance(:, :), precision(:, :)
      real(dp), allocatable :: current_loading_precision(:, :)
      real(dp), allocatable :: current_gamma_precision(:, :), loading_normals(:)
      real(dp), allocatable :: beta_normals(:), uniforms(:), residual(:, :)
      real(dp), allocatable :: posterior_scale(:, :), diagonal_precision(:)
      real(dp), allocatable :: gamma_vector(:), beta_cross(:, :)
      real(dp), allocatable :: structural_coefficient(:), structural_precision(:, :)
      real(dp), allocatable :: structural_normal(:), structural_uniforms(:)
      real(dp), allocatable :: adjusted_y(:, :), structural_residual(:)
      real(dp), allocatable :: structural_response(:, :)
      real(dp), allocatable :: contemporaneous(:, :), impact(:, :)
      real(dp) :: posterior_df
      integer :: k, observations, m, rank, regressors, gamma_count
      integer :: total, draw, retained, item, status, structural_count
      integer :: equation, source, position
      logical :: full_covariance, use_ssvs, use_structural, structural_ssvs

      k = size(y, 1)
      observations = size(y, 2)
      m = size(w, 1)
      rank = size(initial_beta, 2)
      regressors = 0
      if (present(x)) regressors = size(x, 1)
      gamma_count = k*regressors
      total = iterations + burnin
      full_covariance = present(covariance_prior_df)
      use_ssvs = present(tau0)
      use_structural = present(structural_design)
      structural_count = 0
      if (use_structural) structural_count = size(structural_design, 2)
      structural_ssvs = present(structural_tau0)
      if (k < 1 .or. observations < 1 .or. m < 1 .or. rank < 1 .or. &
         size(w, 2) /= observations .or. size(initial_beta, 1) /= m .or. &
         any(shape(initial_covariance) /= [k, k]) .or. &
         any(shape(cointegration_precision) /= [m, m]) .or. &
         any(shape(loading_precision) /= [k, k]) .or. iterations < 1 .or. &
         burnin < 0 .or. shrinkage < 0.0_dp) then
         out%info = 1
         return
      end if
      if (present(x)) then
         if (size(x, 2) /= observations .or. .not. present(gamma_prior_mean) .or. &
            .not. present(gamma_prior_precision)) then
            out%info = 2
            return
         end if
         if (size(gamma_prior_mean) /= gamma_count .or. &
            any(shape(gamma_prior_precision) /= [gamma_count, gamma_count])) then
            out%info = 3
            return
         end if
      else if (present(gamma_prior_mean) .or. present(gamma_prior_precision)) then
         out%info = 4
         return
      end if
      if (use_structural) then
         if (size(structural_design, 1) /= k*observations .or. &
            structural_count < 1 .or. .not. present(structural_prior_mean) .or. &
            .not. present(structural_prior_precision)) then
            out%info = 14
            return
         end if
         if (size(structural_prior_mean) /= structural_count .or. &
            any(shape(structural_prior_precision) /= &
            [structural_count, structural_count])) then
            out%info = 15
            return
         end if
      end if
      if (structural_ssvs) then
         if (.not. use_structural .or. .not. present(structural_tau1) .or. &
            .not. present(structural_inclusion_probability) .or. &
            size(structural_tau0) /= structural_count .or. &
            size(structural_tau1) /= structural_count .or. &
            size(structural_inclusion_probability) /= structural_count) then
            out%info = 16
            return
         end if
      end if
      if (full_covariance) then
         if (covariance_prior_df <= real(k - 1, dp)) then
            out%info = 5
            return
         end if
      else
         if (.not. present(gamma_shape) .or. .not. present(gamma_rate)) then
            out%info = 6
            return
         end if
         if (size(gamma_shape) /= k .or. size(gamma_rate) /= k .or. &
            any(gamma_shape <= 0.0_dp) .or. any(gamma_rate <= 0.0_dp)) then
            out%info = 7
            return
         end if
      end if
      if (use_ssvs) then
         if (gamma_count < 1 .or. .not. present(tau1) .or. &
            .not. present(inclusion_probability)) then
            out%info = 8
            return
         end if
         if (size(tau0) /= gamma_count .or. size(tau1) /= gamma_count .or. &
            size(inclusion_probability) /= gamma_count) then
            out%info = 9
            return
         end if
      end if

      beta = initial_beta
      covariance = initial_covariance
      current_loading_precision = loading_precision
      call invert_matrix(covariance, precision, status)
      if (status /= 0) then
         out%info = 10
         return
      end if
      if (present(x)) current_gamma_precision = gamma_prior_precision
      adjusted_y = y
      structural_response = y
      if (use_structural) then
         structural_coefficient = structural_prior_mean
         structural_precision = structural_prior_precision
         allocate(structural_normal(structural_count), &
            out%structural(structural_count, iterations), &
            out%structural_impact(k, k, iterations))
         if (structural_ssvs) then
            allocate(out%structural_included(structural_count, iterations))
            if (present(structural_selectable)) then
               allocate(structural_uniforms(size(structural_selectable)))
            else
               allocate(structural_uniforms(structural_count))
            end if
         end if
      end if
      allocate(loading_normals(k*rank + gamma_count), beta_normals(m*rank), &
         out%alpha(k*rank, iterations), out%beta(m*rank, iterations), &
         out%pi(k*m, iterations), out%gamma(gamma_count, iterations), &
         out%covariance(k, k, iterations))
      if (use_ssvs) then
         allocate(out%included(gamma_count, iterations))
         if (present(selectable)) then
            allocate(uniforms(size(selectable)))
         else
            allocate(uniforms(gamma_count))
         end if
      end if
      posterior_df = 0.0_dp
      if (full_covariance) posterior_df = covariance_prior_df + real(observations, dp)
      retained = 0
      do draw = 1, total
         if (use_structural) then
            structural_posterior = bvartools_sur_normal_posterior( &
               structural_response, &
               structural_design, spread(precision, 3, observations), &
               structural_prior_mean, structural_precision)
            if (structural_posterior%info /= 0) then
               out%info = 17
               return
            end if
            do item = 1, structural_count
               structural_normal(item) = random_standard_normal()
            end do
            structural_coefficient = bvartools_normal_draw(structural_posterior, &
               structural_normal)
            if (structural_ssvs) then
               do item = 1, size(structural_uniforms)
                  structural_uniforms(item) = random_uniform()
               end do
               structural_selection = bvartools_ssvs(structural_coefficient, &
                  structural_tau0, structural_tau1, &
                  structural_inclusion_probability, structural_uniforms, &
                  structural_selectable)
               structural_precision = 0.0_dp
               do item = 1, structural_count
                  structural_precision(item, item) = &
                     structural_selection%precision(item)
               end do
            end if
            structural_residual = matmul(structural_design, structural_coefficient)
            adjusted_y = y - reshape(structural_residual, [k, observations])
         end if
         do item = 1, size(loading_normals)
            loading_normals(item) = random_standard_normal()
         end do
         do item = 1, size(beta_normals)
            beta_normals(item) = random_standard_normal()
         end do
         if (present(x)) then
            coefficient_draw = bvartools_cointegration_draw(adjusted_y, beta, w, &
               precision, shrinkage, cointegration_precision, &
               current_loading_precision, &
               loading_normals, beta_normals, x, gamma_prior_mean, &
               current_gamma_precision)
         else
            coefficient_draw = bvartools_cointegration_draw(adjusted_y, beta, w, &
               precision, shrinkage, cointegration_precision, &
               current_loading_precision, &
               loading_normals, beta_normals)
         end if
         if (coefficient_draw%info /= 0) then
            out%info = 11
            return
         end if
         beta = coefficient_draw%beta
         if (use_ssvs) then
            gamma_vector = reshape(coefficient_draw%gamma, [gamma_count])
            do item = 1, size(uniforms)
               uniforms(item) = random_uniform()
            end do
            selection = bvartools_ssvs(gamma_vector, tau0, tau1, &
               inclusion_probability, uniforms, selectable)
            if (selection%info /= 0) then
               out%info = 12
               return
            end if
            do item = 1, gamma_count
               current_gamma_precision(item, item) = selection%precision(item)
            end do
         end if
         structural_response = y - matmul(coefficient_draw%pi, w)
         if (present(x)) structural_response = structural_response - &
            matmul(coefficient_draw%gamma, x)
         residual = structural_response
         if (use_structural) residual = residual - &
            reshape(structural_residual, [k, observations])
         if (full_covariance) then
            beta_cross = matmul(transpose(beta), &
               matmul(cointegration_precision, beta))
            posterior_scale = matmul(coefficient_draw%alpha, matmul( &
               shrinkage*beta_cross, transpose(coefficient_draw%alpha))) + &
               matmul(residual, transpose(residual))
            call inverse_wishart_draw(posterior_scale, posterior_df, covariance, &
               precision, status)
            if (status /= 0) then
               out%info = 13
               return
            end if
         else
            allocate(diagonal_precision(k))
            do item = 1, k
               diagonal_precision(item) = random_gamma(gamma_shape(item) + &
                  0.5_dp*real(observations, dp), 1.0_dp/(gamma_rate(item) + &
                  0.5_dp*sum(residual(item, :)**2)))
            end do
            covariance = 0.0_dp
            precision = 0.0_dp
            do item = 1, k
               precision(item, item) = diagonal_precision(item)
               covariance(item, item) = 1.0_dp/diagonal_precision(item)
            end do
            deallocate(diagonal_precision)
         end if
         current_loading_precision = precision
         if (draw > burnin) then
            retained = retained + 1
            out%alpha(:, retained) = reshape(coefficient_draw%alpha, [k*rank])
            out%beta(:, retained) = reshape(coefficient_draw%beta, [m*rank])
            out%pi(:, retained) = reshape(coefficient_draw%pi, [k*m])
            if (gamma_count > 0) then
               out%gamma(:, retained) = reshape(coefficient_draw%gamma, [gamma_count])
            end if
            if (use_structural) out%structural(:, retained) = structural_coefficient
            if (use_structural) then
               contemporaneous = identity_matrix(k)
               position = 0
               do equation = 2, k
                  do source = 1, equation - 1
                     position = position + 1
                     contemporaneous(equation, source) = &
                        structural_coefficient(position)
                  end do
               end do
               call invert_matrix(contemporaneous, impact, status)
               if (status /= 0) then
                  out%info = 18
                  return
               end if
               out%structural_impact(:, :, retained) = impact
            end if
            out%covariance(:, :, retained) = covariance
            if (use_ssvs) out%included(:, retained) = selection%included
            if (structural_ssvs) out%structural_included(:, retained) = &
               structural_selection%included
         end if
      end do
      out%retained_draws = retained
   end function bvartools_bvec_gibbs

   pure function bvartools_dfm_factor_posterior(x, loadings, &
      measurement_precision, transition, factor_precision) result(out)
      !! Compute the joint Gaussian posterior of dynamic-factor paths.
      real(dp), intent(in) :: x(:, :) !! Standardized observed series by observation.
      real(dp), intent(in) :: loadings(:, :) !! Identified measurement loading matrix.
      real(dp), intent(in) :: measurement_precision(:) !! Diagonal measurement-error precisions.
      real(dp), intent(in) :: transition(:, :) !! Factor VAR coefficient matrices by lag.
      real(dp), intent(in) :: factor_precision(:) !! Diagonal factor-innovation precisions.
      type(bvartools_normal_posterior_t) :: out
      real(dp), allocatable :: difference(:, :), innovation_precision(:, :)
      real(dp), allocatable :: posterior_precision(:, :), inverse(:, :)
      real(dp), allocatable :: weighted_loadings(:, :), observation_precision(:, :)
      real(dp), allocatable :: score(:)
      integer :: m, n, observations, lags, dimension, time, lag, item
      integer :: first, last, source_first, source_last, status

      m = size(x, 1)
      observations = size(x, 2)
      n = size(loadings, 2)
      if (m < 1 .or. n < 1 .or. observations < 1 .or. &
         size(loadings, 1) /= m .or. size(measurement_precision) /= m .or. &
         size(factor_precision) /= n .or. size(transition, 1) /= n .or. &
         mod(size(transition, 2), n) /= 0 .or. &
         any(measurement_precision <= 0.0_dp) .or. &
         any(factor_precision <= 0.0_dp)) then
         out%info = 1
         return
      end if
      lags = size(transition, 2)/n
      dimension = n*observations
      difference = identity_matrix(dimension)
      do time = 1, observations
         first = (time - 1)*n + 1
         last = time*n
         do lag = 1, min(lags, time - 1)
            source_first = (time - lag - 1)*n + 1
            source_last = source_first + n - 1
            difference(first:last, source_first:source_last) = &
               -transition(:, (lag - 1)*n + 1:lag*n)
         end do
      end do
      allocate(innovation_precision(dimension, dimension))
      innovation_precision = 0.0_dp
      do time = 1, observations
         first = (time - 1)*n + 1
         do item = 1, n
            innovation_precision(first + item - 1, first + item - 1) = &
               factor_precision(item)
         end do
      end do
      weighted_loadings = loadings
      do item = 1, m
         weighted_loadings(item, :) = measurement_precision(item)*loadings(item, :)
      end do
      observation_precision = matmul(transpose(loadings), weighted_loadings)
      posterior_precision = matmul(transpose(difference), &
         matmul(innovation_precision, difference))
      allocate(score(dimension))
      score = 0.0_dp
      do time = 1, observations
         first = (time - 1)*n + 1
         last = time*n
         posterior_precision(first:last, first:last) = &
            posterior_precision(first:last, first:last) + observation_precision
         score(first:last) = matmul(transpose(weighted_loadings), x(:, time))
      end do
      call invert_matrix(posterior_precision, inverse, status)
      if (status /= 0) then
         out%info = 2
         return
      end if
      out%covariance = inverse
      out%mean = matmul(inverse, score)
   end function bvartools_dfm_factor_posterior

   function bvartools_dfm_gibbs(x, factors, lags, iterations, burnin, &
      loading_prior_precision, measurement_shape, measurement_rate, &
      transition_prior_mean, transition_prior_precision, factor_shape, &
      factor_rate) result(out)
      !! Sample an identified Bayesian dynamic factor model.
      real(dp), intent(in) :: x(:, :) !! Standardized observed series by observation.
      integer, intent(in) :: factors !! Number of latent factors.
      integer, intent(in) :: lags !! Factor VAR lag order.
      integer, intent(in) :: iterations !! Number of retained Gibbs draws.
      integer, intent(in) :: burnin !! Number of discarded Gibbs draws.
      real(dp), intent(in) :: loading_prior_precision(:, :) !! Prior precision of free loadings.
      real(dp), intent(in) :: measurement_shape(:) !! Measurement-precision gamma prior shapes.
      real(dp), intent(in) :: measurement_rate(:) !! Measurement-precision gamma prior rates.
      real(dp), intent(in) :: transition_prior_mean(:) !! Prior mean of vectorized factor VAR coefficients.
      real(dp), intent(in) :: transition_prior_precision(:, :) !! Prior precision of factor VAR coefficients.
      real(dp), intent(in) :: factor_shape(:) !! Factor-precision gamma prior shapes.
      real(dp), intent(in) :: factor_rate(:) !! Factor-precision gamma prior rates.
      type(bvartools_dfm_draws_t) :: out
      type(bvartools_normal_posterior_t) :: posterior, loading_posterior
      real(dp), allocatable :: loadings(:, :), factor_path(:, :), transition(:, :)
      real(dp), allocatable :: measurement_precision(:), factor_precision(:)
      real(dp), allocatable :: standard_normal(:), loading_normal(:)
      real(dp), allocatable :: transition_normal(:), vector_draw(:)
      real(dp), allocatable :: loading_precision(:, :), inverse(:, :), score(:)
      real(dp), allocatable :: response(:), factor_subset(:, :), residual(:, :)
      real(dp), allocatable :: factor_residual(:, :), lagged_factor(:, :)
      real(dp) :: mean_value, sample_var
      integer :: m, observations, n, loading_count, transition_count
      integer :: total, draw, retained, row, item, time, lag, free_count
      integer :: first, last, status

      m = size(x, 1)
      observations = size(x, 2)
      n = factors
      loading_count = (2*m - n - 1)*n/2
      transition_count = n*n*lags
      if (m < 1 .or. observations < 2 .or. n < 1 .or. n > m .or. lags < 0 .or. &
         iterations < 1 .or. burnin < 0 .or. loading_count < 0 .or. &
         any(shape(loading_prior_precision) /= [loading_count, loading_count]) .or. &
         size(measurement_shape) /= m .or. size(measurement_rate) /= m .or. &
         size(transition_prior_mean) /= transition_count .or. &
         any(shape(transition_prior_precision) /= &
         [transition_count, transition_count]) .or. &
         size(factor_shape) /= n .or. size(factor_rate) /= n .or. &
         any(measurement_shape < 0.0_dp) .or. any(measurement_rate <= 0.0_dp) .or. &
         any(factor_shape < 0.0_dp) .or. any(factor_rate <= 0.0_dp)) then
         out%info = 1
         return
      end if
      total = iterations + burnin
      allocate(loadings(m, n), factor_path(n, observations), &
         transition(n, n*lags), measurement_precision(m), factor_precision(n), &
         standard_normal(n*observations), loading_normal(n), &
         transition_normal(transition_count), residual(m, observations), &
         factor_residual(n, observations), lagged_factor(n*lags, observations), &
         out%loadings(m*n, iterations), out%factors(n*observations, iterations), &
         out%measurement_variance(m, iterations), &
         out%transition(transition_count, iterations), &
         out%factor_variance(n, iterations))
      loadings = 0.0_dp
      do item = 1, n
         loadings(item, item) = 1.0_dp
      end do
      factor_path = 0.0_dp
      transition = 0.0_dp
      factor_precision = 1.0_dp
      do row = 1, m
         mean_value = sum(x(row, :))/real(observations, dp)
         sample_var = sum((x(row, :) - mean_value)**2)/real(observations - 1, dp)
         if (sample_var <= tiny(1.0_dp)) then
            out%info = 2
            return
         end if
         measurement_precision(row) = 1.0_dp/sample_var
      end do
      retained = 0
      do draw = 1, total
         posterior = bvartools_dfm_factor_posterior(x, loadings, &
            measurement_precision, transition, factor_precision)
         if (posterior%info /= 0) then
            out%info = 3
            return
         end if
         do item = 1, size(standard_normal)
            standard_normal(item) = random_standard_normal()
         end do
         vector_draw = bvartools_normal_draw(posterior, standard_normal)
         factor_path = reshape(vector_draw, [n, observations])
         first = 1
         do row = 2, m
            free_count = min(row - 1, n)
            last = first + free_count - 1
            factor_subset = factor_path(1:free_count, :)
            if (row <= n) then
               response = x(row, :) - factor_path(row, :)
            else
               response = x(row, :)
            end if
            loading_precision = loading_prior_precision(first:last, first:last) + &
               measurement_precision(row)*matmul(factor_subset, &
               transpose(factor_subset))
            call invert_matrix(loading_precision, inverse, status)
            if (status /= 0) then
               out%info = 4
               return
            end if
            score = measurement_precision(row)*matmul(factor_subset, response)
            loading_posterior%mean = matmul(inverse, score)
            loading_posterior%covariance = inverse
            loading_posterior%info = 0
            do item = 1, free_count
               loading_normal(item) = random_standard_normal()
            end do
            vector_draw = bvartools_normal_draw(loading_posterior, &
               loading_normal(1:free_count))
            loadings(row, 1:free_count) = vector_draw
            first = last + 1
         end do
         residual = x - matmul(loadings, factor_path)
         do row = 1, m
            measurement_precision(row) = random_gamma( &
               measurement_shape(row) + 0.5_dp*real(observations, dp), &
               1.0_dp/(measurement_rate(row) + 0.5_dp*sum(residual(row, :)**2)))
         end do
         factor_residual = factor_path
         do time = 1, observations
            do lag = 1, min(lags, time - 1)
               factor_residual(:, time) = factor_residual(:, time) - &
                  matmul(transition(:, (lag - 1)*n + 1:lag*n), &
                  factor_path(:, time - lag))
            end do
         end do
         do item = 1, n
            factor_precision(item) = random_gamma( &
               factor_shape(item) + 0.5_dp*real(observations, dp), &
               1.0_dp/(factor_rate(item) + &
               0.5_dp*sum(factor_residual(item, :)**2)))
         end do
         if (lags > 0) then
            lagged_factor = 0.0_dp
            do lag = 1, lags
               do time = lag + 1, observations
                  lagged_factor((lag - 1)*n + 1:lag*n, time) = &
                     factor_path(:, time - lag)
               end do
            end do
            loading_precision = 0.0_dp
            do item = 1, n
               loading_precision(item, item) = factor_precision(item)
            end do
            posterior = bvartools_normal_posterior(factor_path, lagged_factor, &
               loading_precision, transition_prior_mean, &
               transition_prior_precision)
            if (posterior%info /= 0) then
               out%info = 5
               return
            end if
            do item = 1, transition_count
               transition_normal(item) = random_standard_normal()
            end do
            vector_draw = bvartools_normal_draw(posterior, transition_normal)
            transition = reshape(vector_draw, [n, n*lags])
         end if
         if (draw > burnin) then
            retained = retained + 1
            out%loadings(:, retained) = reshape(loadings, [m*n])
            out%factors(:, retained) = reshape(factor_path, [n*observations])
            out%measurement_variance(:, retained) = &
               1.0_dp/max(measurement_precision, tiny(1.0_dp))
            if (lags > 0) out%transition(:, retained) = reshape(transition, &
               [transition_count])
            out%factor_variance(:, retained) = &
               1.0_dp/max(factor_precision, tiny(1.0_dp))
         end if
      end do
      out%retained_draws = retained
   end function bvartools_dfm_gibbs

   pure function fevd_kernel(coefficients, covariance_draws, lags, horizon, &
      response_type, normalize_generalized, structural_impact, probability) &
      result(out)
      !! Accumulate identified responses into posterior forecast-error shares.
      real(dp), intent(in) :: coefficients(:, :) !! Vectorized coefficient matrix by posterior draw.
      real(dp), intent(in) :: covariance_draws(:, :, :) !! Innovation covariance by posterior draw.
      integer, intent(in) :: lags !! Number of endogenous lags in each coefficient draw.
      integer, intent(in) :: horizon !! Largest forecast-error horizon.
      character(len=*), intent(in) :: response_type !! Decomposition definition: oir, gir, sir, or sgir.
      logical, intent(in), optional :: normalize_generalized !! Normalize generalized shares to sum to one.
      real(dp), intent(in), optional :: structural_impact(:, :, :) !! Draw-specific inverse contemporaneous matrices.
      real(dp), intent(in), optional :: probability !! Equal-tail posterior interval probability.
      type(bvartools_fevd_t) :: out
      type(bvartools_irf_t) :: identified, base
      real(dp), allocatable :: numerator(:, :), mse(:), phi_row(:), transformed(:)
      real(dp), allocatable :: values(:), ordered(:)
      real(dp) :: interval_probability, tail, total
      integer :: variables, draws_count, draw, step, response, impulse
      logical :: normalize, generalized

      variables = size(covariance_draws, 1)
      draws_count = size(coefficients, 2)
      interval_probability = 0.95_dp
      if (present(probability)) interval_probability = probability
      normalize = .false.
      if (present(normalize_generalized)) normalize = normalize_generalized
      generalized = trim(response_type) == 'gir' .or. trim(response_type) == 'sgir'
      select case (trim(response_type))
      case ('oir')
         identified = irf_kernel(coefficients, covariance_draws, lags, horizon, &
            response_type, scale_by_standard_deviation=.true.)
      case ('gir')
         identified = irf_kernel(coefficients, covariance_draws, lags, horizon, &
            response_type, scale_by_standard_deviation=.true.)
      case ('sir')
         identified = irf_kernel(coefficients, covariance_draws, lags, horizon, &
            response_type, structural_impact=structural_impact)
      case ('sgir')
         identified = irf_kernel(coefficients, covariance_draws, lags, horizon, &
            response_type, scale_by_standard_deviation=.true., &
            structural_impact=structural_impact)
      case default
         out%info = 1
         return
      end select
      if (identified%info /= 0) then
         out%info = 10 + identified%info
         return
      end if
      base = irf_kernel(coefficients, covariance_draws, lags, horizon, 'feir')
      if (base%info /= 0) then
         out%info = 30 + base%info
         return
      end if
      if (interval_probability <= 0.0_dp .or. interval_probability >= 1.0_dp) then
         out%info = 2
         return
      end if

      allocate(out%paths(variables, variables, horizon + 1, draws_count))
      allocate(out%mean(variables, variables, horizon + 1))
      allocate(out%lower(variables, variables, horizon + 1))
      allocate(out%median(variables, variables, horizon + 1))
      allocate(out%upper(variables, variables, horizon + 1))
      allocate(numerator(variables, variables), mse(variables))
      allocate(phi_row(variables), transformed(variables), values(draws_count))
      do draw = 1, draws_count
         numerator = 0.0_dp
         mse = 0.0_dp
         do step = 1, horizon + 1
            numerator = numerator + identified%paths(:, :, step, draw)**2
            if (generalized) then
               do response = 1, variables
                  phi_row = base%paths(response, :, step, draw)
                  if (trim(response_type) == 'sgir') then
                     transformed = matmul(transpose(structural_impact(:, :, draw)), &
                        phi_row)
                  else
                     transformed = phi_row
                  end if
                  mse(response) = mse(response) + dot_product(transformed, &
                     matmul(covariance_draws(:, :, draw), transformed))
               end do
            else
               do response = 1, variables
                  mse(response) = sum(numerator(response, :))
               end do
            end if
            do response = 1, variables
               if (mse(response) <= tiny(1.0_dp)) then
                  out%info = 3
                  return
               end if
               out%paths(response, :, step, draw) = &
                  numerator(response, :)/mse(response)
               if (normalize .and. generalized) then
                  total = sum(out%paths(response, :, step, draw))
                  if (total <= tiny(1.0_dp)) then
                     out%info = 4
                     return
                  end if
                  out%paths(response, :, step, draw) = &
                     out%paths(response, :, step, draw)/total
               end if
            end do
         end do
      end do
      tail = 0.5_dp*(1.0_dp - interval_probability)
      do step = 1, horizon + 1
         do impulse = 1, variables
            do response = 1, variables
               values = out%paths(response, impulse, step, :)
               ordered = sorted(values)
               out%mean(response, impulse, step) = sum(values)/real(draws_count, dp)
               out%lower(response, impulse, step) = quantile(ordered, tail)
               out%median(response, impulse, step) = quantile(ordered, 0.5_dp)
               out%upper(response, impulse, step) = quantile(ordered, 1.0_dp - tail)
            end do
         end do
      end do
      out%probability = interval_probability
      out%response_type = trim(response_type)
      out%normalized = normalize .and. generalized
   end function fevd_kernel

   pure function irf_kernel(coefficients, covariance_draws, lags, horizon, &
      response_type, shock, scale_by_standard_deviation, structural_impact, &
      cumulative, probability) result(out)
      !! Transform posterior VAR moving-average matrices into impulse responses.
      real(dp), intent(in) :: coefficients(:, :) !! Vectorized coefficient matrix by posterior draw.
      real(dp), intent(in) :: covariance_draws(:, :, :) !! Innovation covariance by posterior draw.
      integer, intent(in) :: lags !! Number of endogenous lags in each coefficient draw.
      integer, intent(in) :: horizon !! Largest impulse-response horizon.
      character(len=*), intent(in) :: response_type !! Response definition: feir, oir, gir, sir, or sgir.
      real(dp), intent(in), optional :: shock(:) !! Signed shock magnitude for each impulse variable.
      logical, intent(in), optional :: scale_by_standard_deviation !! Scale shocks by innovation deviations.
      real(dp), intent(in), optional :: structural_impact(:, :, :) !! Draw-specific inverse contemporaneous matrices.
      logical, intent(in), optional :: cumulative !! Accumulate responses through each horizon.
      real(dp), intent(in), optional :: probability !! Equal-tail posterior interval probability.
      type(bvartools_irf_t) :: out
      real(dp), allocatable :: ar(:, :, :), psi(:, :, :), factor(:, :), lower_factor(:, :)
      real(dp), allocatable :: beta(:, :), magnitude(:), values(:), ordered(:)
      real(dp) :: interval_probability, tail, deviation
      integer :: variables, draws_count, regressors, draw, lag, impulse, response
      integer :: step, first, last, status
      logical :: use_deviation, accumulate

      variables = size(covariance_draws, 1)
      draws_count = size(coefficients, 2)
      interval_probability = 0.95_dp
      if (present(probability)) interval_probability = probability
      use_deviation = .false.
      if (present(scale_by_standard_deviation)) use_deviation = &
         scale_by_standard_deviation
      accumulate = .false.
      if (present(cumulative)) accumulate = cumulative
      if (variables < 1 .or. draws_count < 1 .or. lags < 1 .or. horizon < 0 .or. &
         size(covariance_draws, 2) /= variables .or. &
         size(covariance_draws, 3) /= draws_count .or. &
         mod(size(coefficients, 1), variables) /= 0 .or. &
         interval_probability <= 0.0_dp .or. interval_probability >= 1.0_dp) then
         out%info = 1
         return
      end if
      regressors = size(coefficients, 1)/variables
      if (regressors < variables*lags) then
         out%info = 2
         return
      end if
      select case (trim(response_type))
      case ('feir', 'oir', 'gir', 'sir', 'sgir')
      case default
         out%info = 3
         return
      end select
      if (present(shock)) then
         if (size(shock) /= variables) then
            out%info = 4
            return
         end if
      end if
      if (trim(response_type) == 'sir' .or. trim(response_type) == 'sgir') then
         if (.not. present(structural_impact)) then
            out%info = 5
            return
         end if
      end if
      if (present(structural_impact)) then
         if (any(shape(structural_impact) /= [variables, variables, draws_count])) then
            out%info = 6
            return
         end if
      end if

      allocate(out%paths(variables, variables, horizon + 1, draws_count))
      allocate(out%mean(variables, variables, horizon + 1))
      allocate(out%lower(variables, variables, horizon + 1))
      allocate(out%median(variables, variables, horizon + 1))
      allocate(out%upper(variables, variables, horizon + 1))
      allocate(ar(variables, variables, lags), factor(variables, variables))
      allocate(beta(variables, regressors), magnitude(variables), values(draws_count))
      do draw = 1, draws_count
         beta = reshape(coefficients(:, draw), [variables, regressors])
         do lag = 1, lags
            first = (lag - 1)*variables + 1
            last = lag*variables
            ar(:, :, lag) = beta(:, first:last)
         end do
         psi = mts_var_psi(ar, horizon)
         factor = 0.0_dp
         select case (trim(response_type))
         case ('feir')
            factor = identity_matrix(variables)
         case ('oir')
            call cholesky_lower(covariance_draws(:, :, draw), lower_factor, status)
            if (status /= 0) then
               out%info = 7
               return
            end if
            factor = lower_factor
            do impulse = 1, variables
               if (abs(factor(impulse, impulse)) <= tiny(1.0_dp)) then
                  out%info = 8
                  return
               end if
               factor(:, impulse) = factor(:, impulse)/factor(impulse, impulse)
            end do
         case ('gir')
            factor = covariance_draws(:, :, draw)
            do impulse = 1, variables
               if (covariance_draws(impulse, impulse, draw) <= 0.0_dp) then
                  out%info = 9
                  return
               end if
               factor(:, impulse) = factor(:, impulse)/ &
                  covariance_draws(impulse, impulse, draw)
            end do
         case ('sir')
            factor = structural_impact(:, :, draw)
         case ('sgir')
            factor = matmul(structural_impact(:, :, draw), &
               covariance_draws(:, :, draw))
            do impulse = 1, variables
               if (covariance_draws(impulse, impulse, draw) <= 0.0_dp) then
                  out%info = 9
                  return
               end if
               factor(:, impulse) = factor(:, impulse)/ &
                  covariance_draws(impulse, impulse, draw)
            end do
         end select
         magnitude = 1.0_dp
         if (present(shock)) magnitude = shock
         if (use_deviation) then
            do impulse = 1, variables
               if (trim(response_type) == 'oir') then
                  deviation = lower_factor(impulse, impulse)
               else
                  deviation = sqrt(max(covariance_draws(impulse, impulse, draw), &
                     0.0_dp))
               end if
               magnitude(impulse) = magnitude(impulse)*deviation
            end do
         end if
         factor = factor*spread(magnitude, 1, variables)
         do step = 1, horizon + 1
            out%paths(:, :, step, draw) = matmul(psi(:, :, step), factor)
            if (accumulate .and. step > 1) out%paths(:, :, step, draw) = &
               out%paths(:, :, step - 1, draw) + out%paths(:, :, step, draw)
         end do
      end do
      tail = 0.5_dp*(1.0_dp - interval_probability)
      do step = 1, horizon + 1
         do impulse = 1, variables
            do response = 1, variables
               values = out%paths(response, impulse, step, :)
               ordered = sorted(values)
               out%mean(response, impulse, step) = sum(values)/real(draws_count, dp)
               out%lower(response, impulse, step) = quantile(ordered, tail)
               out%median(response, impulse, step) = quantile(ordered, 0.5_dp)
               out%upper(response, impulse, step) = quantile(ordered, 1.0_dp - tail)
            end do
         end do
      end do
      out%probability = interval_probability
      out%response_type = trim(response_type)
   end function irf_kernel

   pure function predictive_kernel(coefficients, covariance_draws, history, &
      lags, standard_normals, future_regressors, structural_impact, probability) &
      result(out)
      !! Apply BVAR coefficient draws recursively and summarize predictive paths.
      real(dp), intent(in) :: coefficients(:, :) !! Vectorized coefficient matrix by posterior draw.
      real(dp), intent(in) :: covariance_draws(:, :, :) !! Innovation covariance by posterior draw.
      real(dp), intent(in) :: history(:, :) !! Endogenous history, ordered from oldest to newest observation.
      integer, intent(in) :: lags !! Number of endogenous lags in each coefficient draw.
      real(dp), intent(in) :: standard_normals(:, :, :) !! Independent normals by variable, horizon, and draw.
      real(dp), intent(in), optional :: future_regressors(:, :) !! Future exogenous and deterministic regressors.
      real(dp), intent(in), optional :: structural_impact(:, :, :) !! Draw-specific inverse contemporaneous matrices.
      real(dp), intent(in), optional :: probability !! Equal-tail posterior interval probability.
      type(bvartools_predictive_t) :: out
      real(dp), allocatable :: beta(:, :), regressor(:), trajectory(:, :)
      real(dp), allocatable :: innovation(:), ordered(:), values(:), linear_mean(:)
      real(dp) :: interval_probability, tail
      integer :: variables, horizon, draws_count, regressors, extra_count
      integer :: draw, step, lag, first, last, variable, status

      variables = size(history, 1)
      horizon = size(standard_normals, 2)
      draws_count = size(coefficients, 2)
      interval_probability = 0.95_dp
      if (present(probability)) interval_probability = probability
      if (variables < 1 .or. size(history, 2) < max(1, lags) .or. &
         lags < 0 .or. horizon < 1 .or. draws_count < 1 .or. &
         mod(size(coefficients, 1), variables) /= 0 .or. &
         any(shape(covariance_draws) /= [variables, variables, draws_count]) .or. &
         any(shape(standard_normals) /= [variables, horizon, draws_count]) .or. &
         interval_probability <= 0.0_dp .or. interval_probability >= 1.0_dp) then
         out%info = 1
         return
      end if
      regressors = size(coefficients, 1)/variables
      extra_count = regressors - variables*lags
      if (extra_count < 0) then
         out%info = 2
         return
      end if
      if (extra_count > 0) then
         if (.not. present(future_regressors)) then
            out%info = 3
            return
         end if
         if (any(shape(future_regressors) /= [extra_count, horizon])) then
            out%info = 4
            return
         end if
      else if (present(future_regressors)) then
         if (size(future_regressors, 1) /= 0 .or. &
            size(future_regressors, 2) /= horizon) then
            out%info = 5
            return
         end if
      end if
      if (present(structural_impact)) then
         if (any(shape(structural_impact) /= [variables, variables, draws_count])) then
            out%info = 6
            return
         end if
      end if

      allocate(out%paths(variables, horizon, draws_count))
      allocate(out%mean(variables, horizon), out%lower(variables, horizon))
      allocate(out%median(variables, horizon), out%upper(variables, horizon))
      allocate(beta(variables, regressors), regressor(regressors))
      allocate(trajectory(variables, max(1, lags) + horizon))
      allocate(innovation(variables), linear_mean(variables), values(draws_count))
      do draw = 1, draws_count
         if (lags > 0) then
            trajectory(:, :lags) = history(:, size(history, 2) - lags + 1:)
         else
            trajectory(:, 1) = history(:, size(history, 2))
         end if
         beta = reshape(coefficients(:, draw), [variables, regressors])
         do step = 1, horizon
            do lag = 1, lags
               first = (lag - 1)*variables + 1
               last = lag*variables
               regressor(first:last) = trajectory(:, lags + step - lag)
            end do
            if (extra_count > 0) then
               regressor(variables*lags + 1:) = future_regressors(:, step)
            end if
            call multivariate_normal_from_standard( &
               [(0.0_dp, variable=1, variables)], covariance_draws(:, :, draw), &
               standard_normals(:, step, draw), innovation, status)
            if (status /= 0) then
               out%info = 7
               return
            end if
            if (regressors > 0) then
               linear_mean = matmul(beta, regressor)
            else
               linear_mean = trajectory(:, step)
            end if
            if (present(structural_impact)) then
               out%paths(:, step, draw) = &
                  matmul(structural_impact(:, :, draw), linear_mean + innovation)
            else
               out%paths(:, step, draw) = linear_mean + innovation
            end if
            trajectory(:, max(1, lags) + step) = out%paths(:, step, draw)
         end do
      end do
      tail = 0.5_dp*(1.0_dp - interval_probability)
      do step = 1, horizon
         do variable = 1, variables
            values = out%paths(variable, step, :)
            ordered = sorted(values)
            out%mean(variable, step) = sum(values)/real(draws_count, dp)
            out%lower(variable, step) = quantile(ordered, tail)
            out%median(variable, step) = quantile(ordered, 0.5_dp)
            out%upper(variable, step) = quantile(ordered, 1.0_dp - tail)
         end do
      end do
      out%probability = interval_probability
   end function predictive_kernel

   pure function make_sur_design(design, variables) result(sur)
      !! Expand regressor observations into the constant-coefficient SUR layout.
      real(dp), intent(in) :: design(:, :) !! Regressors by observation.
      integer, intent(in) :: variables !! Number of response equations.
      real(dp), allocatable :: sur(:, :)
      real(dp), allocatable :: identity(:, :)
      integer :: observations, regressors, time, regressor
      integer :: row_first, row_last, column_first, column_last

      observations = size(design, 2)
      regressors = size(design, 1)
      if (variables < 1) then
         allocate(sur(0, 0))
         return
      end if
      allocate(sur(variables*observations, variables*regressors))
      sur = 0.0_dp
      identity = identity_matrix(variables)
      do time = 1, observations
         row_first = (time - 1)*variables + 1
         row_last = time*variables
         do regressor = 1, regressors
            column_first = (regressor - 1)*variables + 1
            column_last = regressor*variables
            sur(row_first:row_last, column_first:column_last) = &
               design(regressor, time)*identity
         end do
      end do
   end function make_sur_design

   pure function make_structural_design(response) result(design)
      !! Form the off-diagonal contemporaneous-coefficient regression design.
      real(dp), intent(in) :: response(:, :) !! Response variables by observation.
      real(dp), allocatable :: design(:, :)
      integer :: variables, observations, time, source, equation, column, row

      variables = size(response, 1)
      observations = size(response, 2)
      if (variables < 2) then
         allocate(design(0, 0))
         return
      end if
      allocate(design(variables*observations, variables*(variables - 1)/2))
      design = 0.0_dp
      column = 0
      do equation = 2, variables
         do source = 1, equation - 1
            column = column + 1
            do time = 1, observations
               row = (time - 1)*variables + equation
               design(row, column) = -response(source, time)
            end do
         end do
      end do
   end function make_structural_design

   pure logical function valid_deterministic_location(location) result(valid)
      !! Test a deterministic-term placement option.
      character(len=*), intent(in) :: location !! Placement option.

      valid = trim(location) == 'none' .or. trim(location) == 'restricted' .or. &
         trim(location) == 'unrestricted'
   end function valid_deterministic_location

   pure integer function deterministic_count(constant, trend, seasonal, period, &
      target) result(count)
      !! Count deterministic regressors assigned to one VECM coefficient block.
      character(len=*), intent(in) :: constant !! Constant placement.
      character(len=*), intent(in) :: trend !! Trend placement.
      character(len=*), intent(in) :: seasonal !! Seasonal placement.
      integer, intent(in) :: period !! Seasonal period.
      character(len=*), intent(in) :: target !! Placement block being counted.

      count = 0
      if (trim(constant) == trim(target)) count = count + 1
      if (trim(trend) == trim(target)) count = count + 1
      if (trim(seasonal) == trim(target) .and. period >= 2) count = count + period - 1
   end function deterministic_count

   pure subroutine fill_deterministic_column(vector, leading, expected, constant, &
      trend, seasonal, period, aligned_time, original_time, target)
      !! Append one observation of selected deterministic regressors to a vector.
      real(dp), intent(inout) :: vector(:) !! Regressor vector being filled.
      integer, intent(in) :: leading !! Number of preceding nondeterministic regressors.
      integer, intent(in) :: expected !! Number of deterministic regressors to append.
      character(len=*), intent(in) :: constant !! Constant placement.
      character(len=*), intent(in) :: trend !! Trend placement.
      character(len=*), intent(in) :: seasonal !! Seasonal placement.
      integer, intent(in) :: period !! Seasonal period.
      integer, intent(in) :: aligned_time !! Index within the aligned estimation sample.
      integer, intent(in) :: original_time !! Index within the original input series.
      character(len=*), intent(in) :: target !! Placement block being filled.
      integer :: position, phase

      if (expected == 0) return
      position = leading
      if (trim(constant) == trim(target)) then
         position = position + 1
         vector(position) = 1.0_dp
      end if
      if (trim(trend) == trim(target)) then
         position = position + 1
         vector(position) = real(aligned_time, dp)
      end if
      if (trim(seasonal) == trim(target) .and. period >= 2) then
         vector(position + 1:position + period - 1) = 0.0_dp
         phase = modulo(original_time - 1, period) + 1
         if (phase < period) vector(position + phase) = 1.0_dp
      end if
   end subroutine fill_deterministic_column

   pure subroutine symmetric_roots(matrix, root, inverse_root, info)
      !! Compute the square root and inverse square root of a positive-definite matrix.
      real(dp), intent(in) :: matrix(:, :) !! Symmetric positive-definite matrix.
      real(dp), allocatable, intent(out) :: root(:, :) !! Symmetric matrix square root.
      real(dp), allocatable, intent(out) :: inverse_root(:, :) !! Symmetric inverse square root.
      integer, intent(out) :: info !! Zero on success.
      real(dp), allocatable :: values(:), vectors(:, :), diagonal(:, :)
      real(dp) :: tolerance
      integer :: n, item

      n = size(matrix, 1)
      if (size(matrix, 2) /= n) then
         info = 1
         return
      end if
      call symmetric_eigen(0.5_dp*(matrix + transpose(matrix)), values, vectors, info)
      if (info /= 0) return
      tolerance = 100.0_dp*epsilon(1.0_dp)*max(1.0_dp, maxval(abs(values)))
      if (any(values <= tolerance)) then
         info = 2
         return
      end if
      allocate(diagonal(n, n))
      diagonal = 0.0_dp
      do item = 1, n
         diagonal(item, item) = sqrt(values(item))
      end do
      root = matmul(vectors, matmul(diagonal, transpose(vectors)))
      do item = 1, n
         diagonal(item, item) = 1.0_dp/sqrt(values(item))
      end do
      inverse_root = matmul(vectors, matmul(diagonal, transpose(vectors)))
   end subroutine symmetric_roots

   subroutine inverse_wishart_draw(scale, degrees_of_freedom, covariance, precision, info)
      !! Draw an inverse-Wishart covariance by the Bartlett decomposition.
      real(dp), intent(in) :: scale(:, :) !! Inverse-Wishart scale matrix.
      real(dp), intent(in) :: degrees_of_freedom !! Posterior degrees of freedom.
      real(dp), allocatable, intent(out) :: covariance(:, :) !! Drawn covariance matrix.
      real(dp), allocatable, intent(out) :: precision(:, :) !! Drawn precision matrix.
      integer, intent(out) :: info !! Zero on success.
      real(dp), allocatable :: inverse_scale(:, :), lower(:, :), bartlett(:, :), factor(:, :)
      integer :: k, row, column

      k = size(scale, 1)
      call invert_matrix(scale, inverse_scale, info)
      if (info /= 0) return
      call cholesky_lower(inverse_scale, lower, info)
      if (info /= 0) return
      allocate(bartlett(k, k))
      bartlett = 0.0_dp
      do row = 1, k
         bartlett(row, row) = sqrt(2.0_dp*random_gamma( &
            0.5_dp*(degrees_of_freedom - real(row, dp) + 1.0_dp)))
         do column = 1, row - 1
            bartlett(row, column) = random_standard_normal()
         end do
      end do
      factor = matmul(lower, bartlett)
      precision = matmul(factor, transpose(factor))
      call invert_matrix(precision, covariance, info)
   end subroutine inverse_wishart_draw

   pure subroutine update_joint_precision(covariance_state, diagonal_precision, &
      full_precision, info)
      !! Form reduced-form precision matrices from triangular covariance states.
      real(dp), intent(in) :: covariance_state(:, :) !! Packed covariance states by observation.
      real(dp), intent(in) :: diagonal_precision(:, :, :) !! Orthogonal innovation precision by observation.
      real(dp), intent(out) :: full_precision(:, :, :) !! Reduced-form precision by observation.
      integer, intent(out) :: info !! Zero on success.
      real(dp), allocatable :: blocks(:, :), transform(:, :)
      integer :: k, observations, count, time, first, last

      k = size(diagonal_precision, 1)
      observations = size(diagonal_precision, 3)
      count = k*(k - 1)/2
      if (size(diagonal_precision, 2) /= k .or. &
         any(shape(covariance_state) /= [count, observations]) .or. &
         any(shape(full_precision) /= [k, k, observations])) then
         info = 1
         return
      end if
      blocks = bvartools_covar_vector_to_matrix(reshape(covariance_state, &
         [count*observations]), k, observations)
      do time = 1, observations
         first = (time - 1)*k + 1
         last = time*k
         transform = blocks(first:last, first:last)
         full_precision(:, :, time) = matmul(transpose(transform), &
            matmul(diagonal_precision(:, :, time), transform))
      end do
      info = 0
   end subroutine update_joint_precision

   pure subroutine bvartools_structural_impacts(states, variables, observations, &
      impacts, info)
      !! Reconstruct inverse lower-triangular contemporaneous matrices by draw.
      real(dp), intent(in) :: states(:, :) !! Packed structural states by time and draw.
      integer, intent(in) :: variables !! Number of endogenous variables.
      integer, intent(in) :: observations !! Number of state periods per draw.
      real(dp), allocatable, intent(out) :: impacts(:, :, :, :) !! Impact matrices by time and draw.
      integer, intent(out) :: info !! Zero on success.
      real(dp), allocatable :: contemporaneous(:, :), inverse(:, :)
      integer :: count, draws, draw, time, equation, source, position, status

      count = variables*(variables - 1)/2
      draws = size(states, 2)
      if (variables < 2 .or. observations < 1 .or. &
         size(states, 1) /= count*observations) then
         info = 1
         return
      end if
      allocate(impacts(variables, variables, observations, draws))
      do draw = 1, draws
         do time = 1, observations
            contemporaneous = identity_matrix(variables)
            position = (time - 1)*count
            do equation = 2, variables
               do source = 1, equation - 1
                  position = position + 1
                  contemporaneous(equation, source) = states(position, draw)
               end do
            end do
            call invert_matrix(contemporaneous, inverse, status)
            if (status /= 0) then
               info = 2
               return
            end if
            impacts(:, :, time, draw) = inverse
         end do
      end do
      info = 0
   end subroutine bvartools_structural_impacts

   pure subroutine transform_joint_residual(residual, covariance_state, transformed)
      !! Apply lower-triangular covariance transforms to reduced-form residuals.
      real(dp), intent(in) :: residual(:, :) !! Reduced-form residuals by observation.
      real(dp), intent(in) :: covariance_state(:, :) !! Packed covariance states by observation.
      real(dp), intent(out) :: transformed(:, :) !! Orthogonalized residuals by observation.
      real(dp), allocatable :: blocks(:, :), transform(:, :)
      integer :: k, observations, count, time, first, last

      k = size(residual, 1)
      observations = size(residual, 2)
      count = size(covariance_state, 1)
      blocks = bvartools_covar_vector_to_matrix(reshape(covariance_state, &
         [count*observations]), k, observations)
      do time = 1, observations
         first = (time - 1)*k + 1
         last = time*k
         transform = blocks(first:last, first:last)
         transformed(:, time) = matmul(transform, residual(:, time))
      end do
   end subroutine transform_joint_residual

   pure real(dp) function block_quadratic(vector, precision) result(value)
      !! Evaluate a quadratic form using observation-specific precision blocks.
      real(dp), intent(in) :: vector(:) !! Observation-major stacked residual vector.
      real(dp), intent(in) :: precision(:, :, :) !! Precision block by observation.
      integer :: k, observations, time, first, last

      k = size(precision, 1)
      observations = size(precision, 3)
      value = 0.0_dp
      do time = 1, observations
         first = (time - 1)*k + 1
         last = time*k
         value = value + dot_product(vector(first:last), &
            matmul(precision(:, :, time), vector(first:last)))
      end do
   end function block_quadratic

end module bvartools_mod
