! SPDX-License-Identifier: GPL-2.0-or-later
! SPDX-FileComment: Smooth-transition autoregression translated from R tsDyn.
module tsdyn_mod
   !! Nonlinear time-series algorithms translated from the CRAN tsDyn package.
   use kind_mod, only: dp
   use stats_mod, only: ols_fit, quantile, sorted
   use random_mod, only: random_standard_normal
   use linalg_mod, only: invert_matrix, inverse_logdet, symmetric_eigen, &
      cholesky_lower, kronecker_product, outer_product, identity_matrix
   use time_series_stats_mod, only: yule_walker_result_t, yule_walker_fit
   use nts_mod, only: nts_tar_model_t, nts_tar_simulation_t, &
      nts_tar_forecast_t, nts_mtar_model_t, nts_mtar_simulation_t, &
      nts_mtar_forecast_t, nts_utar_estimate, &
      nts_utar_simulate_from_innovations, nts_utar_forecast_draws, &
      nts_mtar_estimate, nts_mtar_simulate_from_standard, &
      nts_mtar_forecast_draws
   use neural_network_mod, only: neural_network_t, neural_network_fit, &
      neural_network_predict, neural_network_parameter_count
   use spline_mod, only: spline_basis_t, spline_basis_create, &
      spline_basis_matrix, spline_basis_values, spline_difference_penalty, &
      penalized_regression_t, penalized_regression_fit
   use special_functions_mod, only: regularized_beta, regularized_gamma_q
   use polynomial_mod, only: polynomial_roots_t, polynomial_roots
   use rolling_forecast_mod, only: rolling_forecast_result_t
   use utils_mod, only: lowercase
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, &
      ieee_quiet_nan
   implicit none
   private

   type, public :: tsdyn_lstar_simulation_t
      !! Simulated logistic smooth-transition autoregression and parameters.
      real(dp), allocatable :: series(:)
      real(dp), allocatable :: innovations(:)
      real(dp), allocatable :: low_coefficients(:)
      real(dp), allocatable :: transition_coefficients(:)
      real(dp), allocatable :: transition_weight(:)
      real(dp) :: gamma = 0.0_dp
      real(dp) :: threshold = 0.0_dp
      integer :: threshold_lag = 1
      integer :: burnin = 0
      integer :: info = 0
   end type tsdyn_lstar_simulation_t

   type, public :: tsdyn_lstar_model_t
      !! Variable-projection least-squares fit of a logistic STAR model.
      real(dp), allocatable :: data(:)
      real(dp), allocatable :: low_coefficients(:)
      real(dp), allocatable :: transition_coefficients(:)
      real(dp), allocatable :: fitted(:)
      real(dp), allocatable :: residuals(:)
      real(dp), allocatable :: transition_weight(:)
      real(dp), allocatable :: threshold_data(:)
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: standard_error(:)
      real(dp), allocatable :: confidence_lower(:)
      real(dp), allocatable :: confidence_upper(:)
      real(dp) :: gamma = 0.0_dp
      real(dp) :: threshold = 0.0_dp
      real(dp) :: residual_standard_deviation = 0.0_dp
      real(dp) :: rss = huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      real(dp) :: bic = huge(1.0_dp)
      integer :: low_order = 0
      integer :: high_order = 0
      integer :: threshold_lag = 1
      integer :: first_fitted = 0
      logical :: include_constant = .true.
      logical :: include_trend = .false.
      logical :: external_threshold = .false.
      logical :: converged = .false.
      integer :: iterations = 0
      integer :: info = 0
   end type tsdyn_lstar_model_t

   type, public :: tsdyn_lstar_forecast_t
      !! Recursive point forecasts and transition weights from an LSTAR model.
      real(dp), allocatable :: mean(:)
      real(dp), allocatable :: transition_weight(:)
      integer :: info = 0
   end type tsdyn_lstar_forecast_t

   type, public :: tsdyn_lstar_selection_t
      !! Candidate information criteria and selected LSTAR model.
      type(tsdyn_lstar_model_t) :: model
      integer, allocatable :: low_order(:)
      integer, allocatable :: high_order(:)
      integer, allocatable :: threshold_lag(:)
      real(dp), allocatable :: aic(:)
      real(dp), allocatable :: bic(:)
      real(dp), allocatable :: gamma(:)
      real(dp), allocatable :: threshold(:)
      integer :: selected = 0
      integer :: info = 0
   end type tsdyn_lstar_selection_t

   type, public :: tsdyn_regime_test_t
      !! Auxiliary-regression F test for an additional STAR regime.
      real(dp) :: statistic = 0.0_dp
      real(dp) :: p_value = 1.0_dp
      integer :: numerator_df = 0
      integer :: denominator_df = 0
      integer :: info = 0
   end type tsdyn_regime_test_t

   type, public :: tsdyn_star_model_t
      !! Additive multi-regime logistic smooth-transition autoregression.
      real(dp), allocatable :: data(:)
      real(dp), allocatable :: coefficients(:, :)
      real(dp), allocatable :: gamma(:)
      real(dp), allocatable :: threshold(:)
      real(dp), allocatable :: fitted(:)
      real(dp), allocatable :: residuals(:)
      real(dp), allocatable :: transition_weight(:, :)
      real(dp), allocatable :: threshold_data(:)
      real(dp) :: residual_standard_deviation = 0.0_dp
      real(dp) :: rss = huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      real(dp) :: bic = huge(1.0_dp)
      integer :: order = 0
      integer :: regimes = 0
      integer :: threshold_lag = 1
      integer :: first_fitted = 0
      logical :: include_constant = .true.
      logical :: include_trend = .false.
      logical :: converged = .false.
      integer :: iterations = 0
      integer :: info = 0
   end type tsdyn_star_model_t

   type, public :: tsdyn_star_forecast_t
      !! Recursive forecasts and transition weights from a multi-regime STAR.
      real(dp), allocatable :: mean(:)
      real(dp), allocatable :: transition_weight(:, :)
      integer :: info = 0
   end type tsdyn_star_forecast_t

   type, public :: tsdyn_llar_model_t
      !! Epsilon diagnostics, fitted values, and selected local-linear AR model.
      real(dp), allocatable :: data(:)
      real(dp), allocatable :: epsilon(:)
      real(dp), allocatable :: normalized_rmse(:)
      real(dp), allocatable :: usable_fraction(:)
      real(dp), allocatable :: average_neighbors(:)
      real(dp), allocatable :: fitted(:)
      integer, allocatable :: neighbor_count(:)
      real(dp) :: selected_epsilon = 0.0_dp
      integer :: selected = 0
      integer :: order = 0
      integer :: delay = 1
      integer :: steps = 1
      integer :: first_fitted = 0
      integer :: info = 0
   end type tsdyn_llar_model_t

   type, public :: tsdyn_llar_forecast_t
      !! Recursive local-linear forecasts and neighborhood diagnostics.
      real(dp), allocatable :: mean(:)
      real(dp), allocatable :: epsilon(:)
      integer, allocatable :: neighbor_count(:)
      integer :: info = 0
   end type tsdyn_llar_forecast_t

   type, public :: tsdyn_aar_model_t
      !! Penalized additive autoregression with one spline component per lag.
      real(dp), allocatable :: data(:)
      type(spline_basis_t), allocatable :: basis(:)
      real(dp), allocatable :: basis_center(:, :)
      real(dp), allocatable :: coefficients(:, :)
      real(dp), allocatable :: lambda(:)
      real(dp), allocatable :: component_df(:)
      real(dp), allocatable :: fitted(:)
      real(dp), allocatable :: residuals(:)
      real(dp) :: intercept = 0.0_dp
      real(dp) :: rss = huge(1.0_dp)
      real(dp) :: gcv = huge(1.0_dp)
      real(dp) :: effective_df = 0.0_dp
      integer :: order = 0
      integer :: delay = 1
      integer :: steps = 1
      integer :: basis_count = 0
      integer :: first_fitted = 0
      integer :: info = 0
   end type tsdyn_aar_model_t

   type, public :: tsdyn_aar_forecast_t
      !! Recursive point forecasts from an additive autoregression.
      real(dp), allocatable :: mean(:)
      integer :: info = 0
   end type tsdyn_aar_forecast_t

   type, public :: tsdyn_tvecm_model_t
      !! Conditional least-squares threshold vector error-correction model.
      real(dp), allocatable :: data(:, :)
      real(dp), allocatable :: coefficients(:, :)
      real(dp), allocatable :: cointegration(:)
      real(dp), allocatable :: cointegration_offset(:)
      real(dp), allocatable :: thresholds(:)
      real(dp), allocatable :: error_correction(:)
      real(dp), allocatable :: fitted(:, :)
      real(dp), allocatable :: residuals(:, :)
      real(dp), allocatable :: covariance(:, :)
      integer, allocatable :: regime(:)
      real(dp) :: beta = 0.0_dp
      real(dp) :: rss = huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      real(dp) :: bic = huge(1.0_dp)
      real(dp) :: trim = 0.05_dp
      integer :: lag = 1
      integer :: threshold_count = 1
      integer :: first_fitted = 0
      logical :: only_error_correction = .false.
      logical :: shared_deterministic = .false.
      logical :: include_constant = .true.
      logical :: include_trend = .false.
      integer :: info = 0
   end type tsdyn_tvecm_model_t

   type, public :: tsdyn_tvecm_forecast_t
      !! Recursive level forecasts and selected regimes from a TVECM.
      real(dp), allocatable :: mean(:, :)
      integer, allocatable :: regime(:)
      integer :: info = 0
   end type tsdyn_tvecm_forecast_t

   type, public :: tsdyn_tvecm_simulation_t
      !! TVECM simulation with supplied innovations and selected regimes.
      real(dp), allocatable :: series(:, :)
      real(dp), allocatable :: innovations(:, :)
      integer, allocatable :: regime(:)
      integer :: info = 0
   end type tsdyn_tvecm_simulation_t

   type, public :: tsdyn_setar_bootstrap_t
      !! Residual-bootstrap SETAR path and optional threshold-model refit.
      real(dp), allocatable :: series(:)
      real(dp), allocatable :: innovations(:)
      type(nts_tar_model_t) :: fitted_model
      logical :: refitted = .false.
      integer :: info = 0
   end type tsdyn_setar_bootstrap_t

   type, public :: tsdyn_tvar_bootstrap_t
      !! Residual-bootstrap TVAR path and optional threshold-model refit.
      real(dp), allocatable :: series(:, :)
      real(dp), allocatable :: innovations(:, :)
      type(nts_mtar_model_t) :: fitted_model
      logical :: refitted = .false.
      integer :: info = 0
   end type tsdyn_tvar_bootstrap_t

   type, public :: tsdyn_tvecm_bootstrap_t
      !! Residual-bootstrap TVECM path and optional threshold-model refit.
      real(dp), allocatable :: series(:, :)
      real(dp), allocatable :: innovations(:, :)
      type(tsdyn_tvecm_model_t) :: fitted_model
      logical :: refitted = .false.
      integer :: info = 0
   end type tsdyn_tvecm_bootstrap_t

   type, public :: tsdyn_threshold_test_t
      !! Supremum threshold-test path and optional bootstrap distribution.
      real(dp), allocatable :: thresholds(:, :)
      real(dp), allocatable :: statistic_path(:)
      real(dp), allocatable :: bootstrap_statistic(:)
      real(dp), allocatable :: critical_values(:)
      real(dp) :: statistic = 0.0_dp
      real(dp) :: p_value = 1.0_dp
      real(dp) :: beta = 0.0_dp
      integer :: selected = 0
      integer :: info = 0
   end type tsdyn_threshold_test_t

   type, public :: tsdyn_regime_count_test_t
      !! Tests comparing linear, two-regime, and three-regime threshold models.
      real(dp) :: fit_measure(3) = huge(1.0_dp)
      real(dp) :: statistic(3) = 0.0_dp
      real(dp), allocatable :: bootstrap_statistic(:, :)
      real(dp) :: p_value(3) = 1.0_dp
      integer :: info = 0
   end type tsdyn_regime_count_test_t

   type, public :: tsdyn_bbc_test_t
      !! Bec-Ben Salem-Carrasco unit-root test over symmetric SETAR thresholds.
      real(dp), allocatable :: threshold(:)
      real(dp), allocatable :: statistic_path(:)
      real(dp), allocatable :: bootstrap_statistic(:)
      real(dp) :: critical_values(3) = 0.0_dp
      real(dp) :: statistic = 0.0_dp
      real(dp) :: p_value = 1.0_dp
      character(len=5) :: method = 'LR'
      integer :: selected = 0
      integer :: info = 0
   end type tsdyn_bbc_test_t

   type, public :: tsdyn_kapshin_test_t
      !! Kapetanios-Shin nonlinear unit-root test over asymmetric thresholds.
      real(dp), allocatable :: thresholds(:, :)
      real(dp), allocatable :: wald_path(:)
      real(dp), allocatable :: bootstrap_statistic(:, :)
      real(dp) :: statistic(3) = 0.0_dp
      real(dp) :: p_value(3) = 1.0_dp
      character(len=5) :: deterministic = 'none'
      integer :: selected = 0
      integer :: info = 0
   end type tsdyn_kapshin_test_t

   type, public :: tsdyn_nnet_model_t
      !! Delay-embedded neural-network autoregression translated from tsDyn.
      type(neural_network_t) :: network
      real(dp), allocatable :: data(:)
      real(dp), allocatable :: fitted(:)
      real(dp), allocatable :: residuals(:)
      real(dp) :: rss = huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      real(dp) :: bic = huge(1.0_dp)
      integer :: order = 0
      integer :: delay = 1
      integer :: steps = 1
      integer :: hidden_count = 0
      integer :: first_fitted = 0
      integer :: info = 0
   end type tsdyn_nnet_model_t

   type, public :: tsdyn_nnet_forecast_t
      !! Recursive point forecasts from a neural-network autoregression.
      real(dp), allocatable :: mean(:)
      integer :: info = 0
   end type tsdyn_nnet_forecast_t

   type, public :: tsdyn_nnet_selection_t
      !! AIC/BIC comparison and selected neural-network autoregression.
      type(tsdyn_nnet_model_t) :: model
      integer, allocatable :: hidden_count(:)
      real(dp), allocatable :: aic(:)
      real(dp), allocatable :: bic(:)
      integer :: selected = 0
      integer :: info = 0
   end type tsdyn_nnet_selection_t

   type, public :: tsdyn_forecast_distribution_t
      !! Simulated nonlinear forecast paths and pointwise distribution summaries.
      real(dp), allocatable :: paths(:, :)
      real(dp), allocatable :: mean(:)
      real(dp), allocatable :: standard_error(:)
      real(dp), allocatable :: lower(:)
      real(dp), allocatable :: upper(:)
      real(dp) :: level = 0.95_dp
      logical :: first_innovation_zero = .false.
      integer :: info = 0
   end type tsdyn_forecast_distribution_t

   type, public :: tsdyn_setar_selection_t
      !! Model-selection results across SETAR specifications.
      type(nts_tar_model_t) :: model
      integer, allocatable :: delay(:)
      integer, allocatable :: threshold_count(:)
      integer, allocatable :: ar_order(:, :)
      integer, allocatable :: threshold_variable_index(:)
      real(dp), allocatable :: thresholds(:, :)
      real(dp), allocatable :: score(:)
      character(len=16) :: criterion = 'AIC'
      character(len=8) :: common = 'none'
      logical :: momentum = .false.
      logical :: same_sample = .false.
      integer :: selected = 0
      integer :: info = 0
   end type tsdyn_setar_selection_t

   type, public :: tsdyn_tvar_selection_t
      !! Model-selection results across TVAR specifications.
      type(nts_mtar_model_t) :: model
      integer, allocatable :: threshold_component(:)
      integer, allocatable :: threshold_variable_index(:)
      integer, allocatable :: delay(:)
      integer, allocatable :: threshold_count(:)
      integer, allocatable :: ar_order(:, :)
      real(dp), allocatable :: thresholds(:, :)
      real(dp), allocatable :: score(:)
      character(len=8) :: criterion = 'AIC'
      logical :: momentum = .false.
      logical :: common_intercept = .false.
      logical :: same_sample = .false.
      integer :: selected = 0
      integer :: info = 0
   end type tsdyn_tvar_selection_t

   type, public :: tsdyn_girf_t
      !! Paired-simulation generalized impulse responses for a nonlinear model.
      real(dp), allocatable :: response(:, :, :)
      integer :: simulations = 0
      integer :: info = 0
   end type tsdyn_girf_t

   type, public :: tsdyn_irf_t
      !! Regime-conditional impulse responses and optional bootstrap intervals.
      real(dp), allocatable :: response(:, :, :)
      real(dp), allocatable :: lower(:, :, :)
      real(dp), allocatable :: upper(:, :, :)
      real(dp) :: level = 0.95_dp
      integer :: regime = 0
      integer :: bootstrap_count = 0
      logical :: cumulative = .false.
      logical :: orthogonalized = .false.
      integer :: info = 0
   end type tsdyn_irf_t

   type, public :: tsdyn_fevd_t
      !! Normalized forecast-error variance shares formed from nonlinear GIRFs.
      real(dp), allocatable :: decomposition(:, :, :)
      integer :: info = 0
   end type tsdyn_fevd_t

   type, public :: tsdyn_delta_test_t
      !! Manzan conditional-independence statistics and permutation inference.
      integer, allocatable :: embedding_dimension(:)
      real(dp), allocatable :: epsilon(:)
      real(dp), allocatable :: statistic(:, :)
      real(dp), allocatable :: permutation_statistic(:, :, :)
      real(dp), allocatable :: p_value(:, :)
      integer :: info = 0
   end type tsdyn_delta_test_t

   type, public :: tsdyn_delta_linear_test_t
      !! Manzan linearity statistics and fitted-AR bootstrap inference.
      integer, allocatable :: embedding_dimension(:)
      real(dp), allocatable :: epsilon(:)
      real(dp), allocatable :: delta(:, :)
      real(dp), allocatable :: linear_delta(:)
      real(dp), allocatable :: statistic(:, :)
      real(dp), allocatable :: bootstrap_statistic(:, :, :)
      real(dp), allocatable :: p_value(:, :)
      real(dp), allocatable :: ar_coefficients(:)
      real(dp) :: innovation_variance = 0.0_dp
      integer :: ar_order = 0
      integer :: info = 0
   end type tsdyn_delta_linear_test_t

   type, public :: tsdyn_regime_diagnostics_t
      !! Regime equilibrium means, characteristic roots, and stability flags.
      character(len=16), allocatable :: regime_label(:)
      real(dp), allocatable :: equilibrium_mean(:)
      complex(dp), allocatable :: roots(:, :)
      real(dp), allocatable :: root_modulus(:, :)
      integer, allocatable :: ar_order(:)
      logical, allocatable :: stable(:)
      logical, allocatable :: unit_root(:)
      logical, allocatable :: mean_defined(:)
      integer :: info = 0
   end type tsdyn_regime_diagnostics_t

   type, public :: tsdyn_rank_test_t
      !! Johansen rank statistics and Doornik gamma-approximation p-values.
      real(dp), allocatable :: eigenvalues(:)
      real(dp), allocatable :: trace_statistic(:)
      real(dp), allocatable :: trace_p_value(:)
      real(dp), allocatable :: adjusted_trace_p_value(:)
      real(dp), allocatable :: eigen_statistic(:)
      real(dp), allocatable :: eigen_p_value(:)
      real(dp) :: selected_p_value = 1.0_dp
      integer :: selected_rank = 0
      integer :: tested_rank = -1
      integer :: sample_size = 0
      character(len=4) :: deterministic = 'H_lc'
      character(len=5) :: test_type = 'eigen'
      integer :: info = 0
   end type tsdyn_rank_test_t

   type, public :: tsdyn_rank_selection_t
      !! Information-criterion grid and selected cointegration ranks and lags.
      real(dp), allocatable :: aic(:, :)
      real(dp), allocatable :: bic(:, :)
      real(dp), allocatable :: hq(:, :)
      real(dp), allocatable :: log_likelihood(:, :)
      integer :: selected_rank(3) = 0
      integer :: selected_lag(3) = 0
      integer :: rank_max = 0
      integer :: lag_max = 0
      logical :: same_sample = .true.
      character(len=5) :: deterministic = 'const'
      character(len=3) :: fit_measure = 'SSR'
      integer :: info = 0
   end type tsdyn_rank_selection_t

   type, public :: tsdyn_setar_variance_t
      !! Pooled and regime-specific SETAR residual variances.
      real(dp), allocatable :: regime(:)
      integer, allocatable :: degrees_of_freedom(:)
      real(dp) :: pooled = 0.0_dp
      integer :: pooled_degrees_of_freedom = 0
      character(len=3) :: adjustment = 'OLS'
      integer :: info = 0
   end type tsdyn_setar_variance_t

   type, public :: tsdyn_setar_inference_t
      !! Conditional SETAR coefficient covariance blocks and standard errors.
      real(dp), allocatable :: covariance(:, :, :)
      real(dp), allocatable :: standard_error(:, :)
      logical, allocatable :: active(:, :)
      integer, allocatable :: parameter_count(:)
      type(tsdyn_setar_variance_t) :: residual_variance
      character(len=6) :: variance_mode = 'pooled'
      integer :: info = 0
   end type tsdyn_setar_inference_t

   type, public :: tsdyn_tvar_inference_t
      !! Conditional TVAR Kronecker covariance blocks and standard errors.
      real(dp), allocatable :: coefficient_covariance(:, :, :)
      real(dp), allocatable :: standard_error(:, :, :)
      logical, allocatable :: active(:, :)
      integer, allocatable :: parameter_count(:)
      real(dp), allocatable :: pooled_residual_covariance(:, :)
      real(dp), allocatable :: regime_residual_covariance(:, :, :)
      integer, allocatable :: degrees_of_freedom(:)
      integer :: pooled_degrees_of_freedom = 0
      character(len=6) :: variance_mode = 'pooled'
      character(len=3) :: adjustment = 'OLS'
      integer :: info = 0
   end type tsdyn_tvar_inference_t

   type, public :: tsdyn_tvecm_inference_t
      !! Conditional TVECM coefficient inference after threshold selection.
      real(dp), allocatable :: coefficient_covariance(:, :)
      real(dp), allocatable :: standard_error(:, :)
      real(dp), allocatable :: t_statistic(:, :)
      real(dp), allocatable :: p_value(:, :)
      real(dp), allocatable :: pooled_residual_covariance(:, :)
      real(dp), allocatable :: regime_residual_covariance(:, :, :)
      integer, allocatable :: regime_observations(:)
      integer, allocatable :: regime_degrees_of_freedom(:)
      integer :: pooled_degrees_of_freedom = 0
      integer :: parameter_count = 0
      character(len=6) :: variance_mode = 'pooled'
      character(len=3) :: adjustment = 'OLS'
      integer :: info = 0
   end type tsdyn_tvecm_inference_t

   type, public :: tsdyn_regime_path_t
      !! Full-series regime labels, transition weights, and alignment validity.
      integer, allocatable :: regime(:)
      real(dp), allocatable :: transition_weight(:, :)
      logical, allocatable :: valid(:)
      integer :: first_valid = 0
      integer :: regime_count = 0
      integer :: info = 0
   end type tsdyn_regime_path_t

   type :: tsdyn_johansen_kernel_t
      !! Internal reduced-rank eigenvalues and unrestricted covariance determinant.
      real(dp), allocatable :: eigenvalues(:)
      real(dp) :: log_determinant_s00 = 0.0_dp
      integer :: sample_size = 0
      integer :: unrestricted_count = 0
      integer :: restricted_count = 0
      integer :: info = 0
   end type tsdyn_johansen_kernel_t

   public :: tsdyn_logistic_transition
   public :: tsdyn_lstar_simulate_from_innovations, tsdyn_lstar_simulate
   public :: tsdyn_lstar_fit, tsdyn_lstar_forecast
   public :: tsdyn_lstar_select, tsdyn_lstar_gradient, tsdyn_lstar_regime_test
   public :: tsdyn_star_fit, tsdyn_star_forecast, tsdyn_star_gradient
   public :: tsdyn_star_regime_test, tsdyn_star_add_regime
   public :: tsdyn_llar_fit, tsdyn_llar_forecast
   public :: tsdyn_aar_fit, tsdyn_aar_component, tsdyn_aar_forecast
   public :: tsdyn_tvecm_fit, tsdyn_tvecm_forecast
   public :: tsdyn_tvecm_simulate_from_innovations
   public :: tsdyn_setar_bootstrap, tsdyn_tvar_bootstrap
   public :: tsdyn_tvecm_bootstrap
   public :: nts_tar_model_t, nts_tar_simulation_t, nts_tar_forecast_t
   public :: nts_mtar_model_t, nts_mtar_simulation_t, nts_mtar_forecast_t
   public :: tsdyn_setar_fit, tsdyn_setar_simulate_from_innovations
   public :: tsdyn_setar_forecast_draws
   public :: tsdyn_setar_select
   public :: tsdyn_setar_restricted_fit
   public :: tsdyn_tvar_fit, tsdyn_tvar_simulate_from_standard
   public :: tsdyn_tvar_forecast_draws
   public :: tsdyn_tvar_select
   public :: tsdyn_tvar_restricted_fit
   public :: tsdyn_hansen_seo_test, tsdyn_seo_test
   public :: tsdyn_setar_linearity_test, tsdyn_tvar_lr_test
   public :: tsdyn_bbc_test, tsdyn_kapshin_test
   public :: tsdyn_nnet_fit, tsdyn_nnet_predict
   public :: tsdyn_nnet_forecast, tsdyn_nnet_select
   public :: tsdyn_setar_forecast_distribution
   public :: tsdyn_lstar_forecast_distribution
   public :: tsdyn_star_forecast_distribution
   public :: tsdyn_llar_forecast_distribution
   public :: tsdyn_aar_forecast_distribution
   public :: tsdyn_nnet_forecast_distribution
   public :: rolling_forecast_result_t
   public :: tsdyn_setar_rolling_forecast, tsdyn_lstar_rolling_forecast
   public :: tsdyn_star_rolling_forecast, tsdyn_llar_rolling_forecast
   public :: tsdyn_aar_rolling_forecast, tsdyn_nnet_rolling_forecast
   public :: tsdyn_tvar_girf_from_innovations
   public :: tsdyn_tvecm_girf_from_innovations, tsdyn_nonlinear_fevd
   public :: tsdyn_setar_girf_from_innovations
   public :: tsdyn_setar_regime_irf, tsdyn_tvar_regime_irf
   public :: tsdyn_tvecm_regime_irf
   public :: tsdyn_setar_irf_bootstrap, tsdyn_tvar_irf_bootstrap
   public :: tsdyn_tvecm_irf_bootstrap
   public :: tsdyn_delta, tsdyn_delta_test
   public :: tsdyn_delta_linear, tsdyn_delta_linear_test
   public :: tsdyn_ar_regime_diagnostics, tsdyn_setar_regime_diagnostics
   public :: tsdyn_lstar_regime_diagnostics, tsdyn_star_regime_diagnostics
   public :: tsdyn_rank_test, tsdyn_rank_select
   public :: tsdyn_setar_residual_variance, tsdyn_setar_inference
   public :: tsdyn_tvar_inference, tsdyn_tvecm_inference
   public :: tsdyn_setar_regimes, tsdyn_lstar_regimes, tsdyn_star_regimes
   public :: tsdyn_tvar_regimes, tsdyn_tvecm_regimes

contains

   pure function tsdyn_setar_regimes(model, series, threshold_variable) &
      result(out)
      !! Classify SETAR response times by their lagged threshold driver.
      type(nts_tar_model_t), intent(in) :: model !! Fitted SETAR model.
      real(dp), intent(in), optional :: series(:) !! Replacement series to classify.
      real(dp), intent(in), optional :: threshold_variable(:) !! External threshold-driver series.
      type(tsdyn_regime_path_t) :: out
      real(dp), allocatable :: data(:)
      real(dp) :: threshold_value, nan_value
      integer :: observations, maximum_order, time

      if (present(series)) then
         data = series
      else if (allocated(model%data)) then
         data = model%data
      else
         out%info = 1
         return
      end if
      observations = size(data)
      if (.not. allocated(model%ar_order) .or. &
         .not. allocated(model%thresholds) .or. model%delay < 1 .or. &
         observations < 1 .or. .not. all(ieee_is_finite(data))) then
         out%info = 1
         return
      end if
      if (present(threshold_variable)) then
         if (size(threshold_variable) /= observations .or. &
            .not. all(ieee_is_finite(threshold_variable))) then
            out%info = 2
            return
         end if
      end if
      maximum_order = max(1, maxval(model%ar_order))
      out%first_valid = max(maximum_order, model%delay) + 1
      if (observations < out%first_valid) then
         out%info = 1
         return
      end if
      nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
      allocate(out%regime(observations), out%valid(observations), &
         out%transition_weight(observations, 0))
      out%regime = 0
      out%valid = .false.
      out%transition_weight = nan_value
      do time = out%first_valid, observations
         if (present(threshold_variable)) then
            threshold_value = threshold_variable(time - model%delay)
         else
            threshold_value = data(time - model%delay)
         end if
         out%regime(time) = inference_threshold_regime(threshold_value, &
            model%thresholds)
         out%valid(time) = .true.
      end do
      out%regime_count = size(model%thresholds) + 1
      out%info = 0
   end function tsdyn_setar_regimes

   pure function tsdyn_lstar_regimes(model, series, threshold_variable) &
      result(out)
      !! Return LSTAR transition probabilities and low/high regime labels.
      type(tsdyn_lstar_model_t), intent(in) :: model !! Fitted LSTAR model.
      real(dp), intent(in), optional :: series(:) !! Replacement series to classify.
      real(dp), intent(in), optional :: threshold_variable(:) !! External contemporaneous threshold driver.
      type(tsdyn_regime_path_t) :: out
      real(dp), allocatable :: data(:)
      real(dp) :: threshold_value, weight, nan_value
      integer :: observations, maximum_order, time

      if (present(series)) then
         data = series
      else if (allocated(model%data)) then
         data = model%data
      else
         out%info = 1
         return
      end if
      observations = size(data)
      if (model%threshold_lag < 1 .or. model%gamma <= 0.0_dp .or. &
         observations < 1 .or. .not. all(ieee_is_finite(data))) then
         out%info = 1
         return
      end if
      if (model%external_threshold .and. .not. present(threshold_variable)) then
         out%info = 2
         return
      end if
      if (present(threshold_variable)) then
         if (size(threshold_variable) /= observations .or. &
            .not. all(ieee_is_finite(threshold_variable))) then
            out%info = 2
            return
         end if
      end if
      maximum_order = max(model%low_order, model%high_order)
      out%first_valid = max(maximum_order, model%threshold_lag) + 1
      if (observations < out%first_valid) then
         out%info = 1
         return
      end if
      nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
      allocate(out%regime(observations), out%valid(observations), &
         out%transition_weight(observations, 1))
      out%regime = 0
      out%valid = .false.
      out%transition_weight = nan_value
      do time = out%first_valid, observations
         if (present(threshold_variable)) then
            threshold_value = threshold_variable(time)
         else
            threshold_value = data(time - model%threshold_lag)
         end if
         weight = tsdyn_logistic_transition(threshold_value, model%gamma, &
            model%threshold)
         out%transition_weight(time, 1) = weight
         out%regime(time) = merge(1, 2, weight <= 0.5_dp)
         out%valid(time) = .true.
      end do
      out%regime_count = 2
      out%info = 0
   end function tsdyn_lstar_regimes

   pure function tsdyn_star_regimes(model, series) result(out)
      !! Return additive STAR transition weights and dominant asymptotic regimes.
      type(tsdyn_star_model_t), intent(in) :: model !! Fitted additive STAR model.
      real(dp), intent(in), optional :: series(:) !! Replacement series to classify.
      type(tsdyn_regime_path_t) :: out
      real(dp), allocatable :: data(:)
      real(dp) :: threshold_value, nan_value
      integer :: observations, time, component

      if (present(series)) then
         data = series
      else if (allocated(model%data)) then
         data = model%data
      else
         out%info = 1
         return
      end if
      observations = size(data)
      if (model%regimes < 2 .or. model%order < 1 .or. &
         model%threshold_lag < 1 .or. .not. allocated(model%gamma) .or. &
         .not. allocated(model%threshold) .or. &
         size(model%gamma) /= model%regimes - 1 .or. &
         size(model%threshold) /= model%regimes - 1 .or. &
         .not. all(ieee_is_finite(data))) then
         out%info = 1
         return
      end if
      out%first_valid = max(model%order, model%threshold_lag) + 1
      if (observations < out%first_valid) then
         out%info = 1
         return
      end if
      nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
      allocate(out%regime(observations), out%valid(observations), &
         out%transition_weight(observations, model%regimes - 1))
      out%regime = 0
      out%valid = .false.
      out%transition_weight = nan_value
      do time = out%first_valid, observations
         threshold_value = data(time - model%threshold_lag)
         do component = 1, model%regimes - 1
            out%transition_weight(time, component) = tsdyn_logistic_transition( &
               threshold_value, model%gamma(component), &
               model%threshold(component))
         end do
         out%regime(time) = 1 + count(out%transition_weight(time, :) > 0.5_dp)
         out%valid(time) = .true.
      end do
      out%regime_count = model%regimes
      out%info = 0
   end function tsdyn_star_regimes

   pure function tsdyn_tvar_regimes(model, series, threshold_variable) &
      result(out)
      !! Classify TVAR response times by their lagged scalar threshold driver.
      type(nts_mtar_model_t), intent(in) :: model !! Fitted TVAR model.
      real(dp), intent(in), optional :: series(:, :) !! Replacement multivariate series.
      real(dp), intent(in), optional :: threshold_variable(:) !! External threshold-driver series.
      type(tsdyn_regime_path_t) :: out
      real(dp), allocatable :: data(:, :)
      real(dp) :: threshold_value, nan_value
      integer :: observations, maximum_order, time

      if (present(series)) then
         data = series
      else if (allocated(model%data)) then
         data = model%data
      else
         out%info = 1
         return
      end if
      observations = size(data, 1)
      if (.not. allocated(model%ar_order) .or. &
         .not. allocated(model%thresholds) .or. model%delay < 1 .or. &
         model%threshold_component < 1 .or. &
         model%threshold_component > size(data, 2) .or. &
         .not. all(ieee_is_finite(data))) then
         out%info = 1
         return
      end if
      if (present(threshold_variable)) then
         if (size(threshold_variable) /= observations .or. &
            .not. all(ieee_is_finite(threshold_variable))) then
            out%info = 2
            return
         end if
      end if
      maximum_order = max(1, maxval(model%ar_order))
      out%first_valid = max(maximum_order, model%delay) + 1
      if (observations < out%first_valid) then
         out%info = 1
         return
      end if
      nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
      allocate(out%regime(observations), out%valid(observations), &
         out%transition_weight(observations, 0))
      out%regime = 0
      out%valid = .false.
      out%transition_weight = nan_value
      do time = out%first_valid, observations
         if (present(threshold_variable)) then
            threshold_value = threshold_variable(time - model%delay)
         else
            threshold_value = data(time - model%delay, &
               model%threshold_component)
         end if
         out%regime(time) = inference_threshold_regime(threshold_value, &
            model%thresholds)
         out%valid(time) = .true.
      end do
      out%regime_count = size(model%thresholds) + 1
      out%info = 0
   end function tsdyn_tvar_regimes

   pure function tsdyn_tvecm_regimes(model, series) result(out)
      !! Classify TVECM response times by the lagged error-correction term.
      type(tsdyn_tvecm_model_t), intent(in) :: model !! Fitted TVECM model.
      real(dp), intent(in), optional :: series(:, :) !! Replacement multivariate level series.
      type(tsdyn_regime_path_t) :: out
      real(dp), allocatable :: data(:, :), cointegration(:)
      real(dp) :: error_correction, nan_value
      integer :: observations, time

      if (present(series)) then
         data = series
      else if (allocated(model%data)) then
         data = model%data
      else
         out%info = 1
         return
      end if
      observations = size(data, 1)
      if (size(data, 2) < 2 .or. model%lag < 0 .or. &
         .not. allocated(model%thresholds) .or. &
         .not. all(ieee_is_finite(data))) then
         out%info = 1
         return
      end if
      out%first_valid = model%lag + 2
      if (observations < out%first_valid) then
         out%info = 1
         return
      end if
      cointegration = tvecm_cointegration_vector(model, size(data, 2))
      if (size(cointegration) /= size(data, 2)) then
         out%info = 1
         return
      end if
      nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
      allocate(out%regime(observations), out%valid(observations), &
         out%transition_weight(observations, 0))
      out%regime = 0
      out%valid = .false.
      out%transition_weight = nan_value
      do time = out%first_valid, observations
         error_correction = dot_product(cointegration, data(time - 1, :))
         error_correction = error_correction - &
            tvecm_cointegration_offset_value(model, time - 1)
         out%regime(time) = inference_threshold_regime(error_correction, &
            model%thresholds)
         out%valid(time) = .true.
      end do
      out%regime_count = size(model%thresholds) + 1
      out%info = 0
   end function tsdyn_tvecm_regimes

   pure function tsdyn_setar_residual_variance(model, adjustment, &
      threshold_variable) result(out)
      !! Estimate pooled and regime-specific SETAR innovation variances.
      type(nts_tar_model_t), intent(in) :: model !! Fitted SETAR model.
      character(len=*), intent(in), optional :: adjustment !! Denominator choice: OLS or ML.
      real(dp), intent(in), optional :: threshold_variable(:) !! External threshold series used by the fit.
      type(tsdyn_setar_variance_t) :: out
      character(len=8) :: selected_adjustment
      real(dp), allocatable :: sum_squares(:)
      real(dp) :: threshold_value
      integer :: regimes, maximum_order, start, usable, row, time, regime
      integer :: parameter_total, denominator

      if (model%representation /= 'level') then
         out%info = 3
         return
      end if
      if (allocated(model%include_trend)) then
         if (any(model%include_trend)) then
            out%info = 3
            return
         end if
      end if
      selected_adjustment = 'ols'
      if (present(adjustment)) selected_adjustment = &
         lowercase(adjustl(adjustment))
      if (.not. allocated(model%data) .or. .not. allocated(model%residuals) .or. &
         .not. allocated(model%ar_order) .or. &
         .not. allocated(model%include_mean) .or. &
         .not. allocated(model%thresholds) .or. &
         (trim(selected_adjustment) /= 'ols' .and. &
         trim(selected_adjustment) /= 'ml')) then
         out%info = 1
         return
      end if
      regimes = size(model%ar_order)
      maximum_order = max(1, maxval(model%ar_order))
      start = max(maximum_order, model%delay) + 1
      usable = size(model%data) - start + 1
      if (regimes < 2 .or. size(model%include_mean) /= regimes .or. &
         size(model%thresholds) /= regimes - 1 .or. usable < 1 .or. &
         size(model%residuals) /= usable) then
         out%info = 1
         return
      end if
      if (present(threshold_variable)) then
         if (size(threshold_variable) /= size(model%data) .or. &
            .not. all(ieee_is_finite(threshold_variable))) then
            out%info = 2
            return
         end if
      end if
      allocate(sum_squares(regimes), out%regime(regimes), &
         out%degrees_of_freedom(regimes))
      sum_squares = 0.0_dp
      out%degrees_of_freedom = 0
      do row = 1, usable
         time = start + row - 1
         if (present(threshold_variable)) then
            threshold_value = threshold_variable(time - model%delay)
         else
            threshold_value = model%data(time - model%delay)
         end if
         regime = inference_threshold_regime(threshold_value, model%thresholds)
         sum_squares(regime) = sum_squares(regime) + model%residuals(row)**2
         out%degrees_of_freedom(regime) = &
            out%degrees_of_freedom(regime) + 1
      end do
      parameter_total = sum(model%ar_order) + count(model%include_mean)
      do regime = 1, regimes
         denominator = out%degrees_of_freedom(regime)
         if (trim(selected_adjustment) == 'ols') denominator = denominator - &
            model%ar_order(regime) - merge(1, 0, model%include_mean(regime))
         if (denominator <= 0) then
            out%info = 3
            return
         end if
         out%regime(regime) = sum_squares(regime)/real(denominator, dp)
         out%degrees_of_freedom(regime) = denominator
      end do
      denominator = usable
      if (trim(selected_adjustment) == 'ols') denominator = denominator - &
         parameter_total
      if (denominator <= 0) then
         out%info = 3
         return
      end if
      out%pooled = sum(sum_squares)/real(denominator, dp)
      out%pooled_degrees_of_freedom = denominator
      if (trim(selected_adjustment) == 'ols') then
         out%adjustment = 'OLS'
      else
         out%adjustment = 'ML'
      end if
      out%info = 0
   end function tsdyn_setar_residual_variance

   pure function tsdyn_setar_inference(model, variance_mode, adjustment, &
      threshold_variable) result(out)
      !! Compute conditional SETAR coefficient covariance and standard errors.
      type(nts_tar_model_t), intent(in) :: model !! Fitted SETAR model.
      character(len=*), intent(in), optional :: variance_mode !! Use pooled or regime residual variance.
      character(len=*), intent(in), optional :: adjustment !! Denominator choice: OLS or ML.
      real(dp), intent(in), optional :: threshold_variable(:) !! External threshold series used by the fit.
      type(tsdyn_setar_inference_t) :: out
      real(dp), allocatable :: gram(:, :, :), inverse(:, :), regressor(:)
      real(dp) :: threshold_value, selected_variance
      character(len=8) :: selected_mode
      integer :: regimes, maximum_order, maximum_parameters, start, usable
      integer :: row, time, regime, parameter_count, column, lag, status

      if (model%representation /= 'level') then
         out%info = 3
         return
      end if
      if (allocated(model%include_trend)) then
         if (any(model%include_trend)) then
            out%info = 3
            return
         end if
      end if
      selected_mode = 'pooled'
      if (present(variance_mode)) selected_mode = &
         lowercase(adjustl(variance_mode))
      if (trim(selected_mode) /= 'pooled' .and. trim(selected_mode) /= 'regime') then
         out%info = 1
         return
      end if
      out%residual_variance = tsdyn_setar_residual_variance(model, adjustment, &
         threshold_variable)
      if (out%residual_variance%info /= 0) then
         out%info = 10 + out%residual_variance%info
         return
      end if
      regimes = size(model%ar_order)
      maximum_order = max(1, maxval(model%ar_order))
      maximum_parameters = maxval(model%ar_order + &
         merge(1, 0, model%include_mean))
      start = max(maximum_order, model%delay) + 1
      usable = size(model%data) - start + 1
      allocate(gram(maximum_parameters, maximum_parameters, regimes), &
         regressor(maximum_parameters), &
         out%covariance(maximum_parameters, maximum_parameters, regimes), &
         out%standard_error(regimes, size(model%coefficients, 2)), &
         out%active(regimes, size(model%coefficients, 2)), &
         out%parameter_count(regimes))
      gram = 0.0_dp
      out%covariance = 0.0_dp
      out%standard_error = 0.0_dp
      out%active = .false.
      out%parameter_count = model%ar_order + merge(1, 0, model%include_mean)
      do row = 1, usable
         time = start + row - 1
         if (present(threshold_variable)) then
            threshold_value = threshold_variable(time - model%delay)
         else
            threshold_value = model%data(time - model%delay)
         end if
         regime = inference_threshold_regime(threshold_value, model%thresholds)
         parameter_count = model%ar_order(regime) + &
            merge(1, 0, model%include_mean(regime))
         regressor = 0.0_dp
         column = 0
         if (model%include_mean(regime)) then
            column = 1
            regressor(1) = 1.0_dp
         end if
         do lag = 1, model%ar_order(regime)
            regressor(column + lag) = model%data(time - lag)
         end do
         gram(:parameter_count, :parameter_count, regime) = &
            gram(:parameter_count, :parameter_count, regime) + &
            outer_product(regressor(:parameter_count), regressor(:parameter_count))
      end do
      do regime = 1, regimes
         parameter_count = model%ar_order(regime) + &
            merge(1, 0, model%include_mean(regime))
         call invert_matrix(gram(:parameter_count, :parameter_count, regime), &
            inverse, status)
         if (status /= 0) then
            out%info = 20 + regime
            return
         end if
         if (trim(selected_mode) == 'pooled') then
            selected_variance = out%residual_variance%pooled
         else
            selected_variance = out%residual_variance%regime(regime)
         end if
         out%covariance(:parameter_count, :parameter_count, regime) = &
            selected_variance*inverse
         column = 0
         if (model%include_mean(regime)) then
            column = 1
            out%active(regime, 1) = .true.
            out%standard_error(regime, 1) = sqrt(max(0.0_dp, &
               out%covariance(1, 1, regime)))
         end if
         do lag = 1, model%ar_order(regime)
            out%active(regime, lag + 1) = .true.
            out%standard_error(regime, lag + 1) = sqrt(max(0.0_dp, &
               out%covariance(column + lag, column + lag, regime)))
         end do
         deallocate(inverse)
      end do
      out%variance_mode = selected_mode
      out%info = 0
   end function tsdyn_setar_inference

   pure function tsdyn_tvar_inference(model, variance_mode, adjustment, &
      threshold_variable) result(out)
      !! Compute conditional TVAR Kronecker covariance blocks and standard errors.
      type(nts_mtar_model_t), intent(in) :: model !! Fitted TVAR model.
      character(len=*), intent(in), optional :: variance_mode !! Use pooled or regime residual covariance.
      character(len=*), intent(in), optional :: adjustment !! Denominator choice: OLS or ML.
      real(dp), intent(in), optional :: threshold_variable(:) !! External threshold series used by the fit.
      type(tsdyn_tvar_inference_t) :: out
      real(dp), allocatable :: gram(:, :, :), inverse(:, :), regressor(:)
      real(dp), allocatable :: residual_sum(:, :, :), selected_covariance(:, :)
      real(dp), allocatable :: compact_covariance(:, :)
      real(dp) :: threshold_value
      character(len=8) :: selected_mode, selected_adjustment
      integer :: observations, variables, regimes, maximum_order
      integer :: maximum_parameters, stored_parameters, start, usable
      integer :: row, time, regime, parameter_count, column, lag, variable
      integer :: equation, compact_index, denominator, parameter_total, status

      if (allocated(model%include_trend)) then
         if (any(model%include_trend)) then
            out%info = 3
            return
         end if
      end if
      selected_mode = 'pooled'
      if (present(variance_mode)) selected_mode = lowercase(adjustl(variance_mode))
      selected_adjustment = 'ols'
      if (present(adjustment)) selected_adjustment = &
         lowercase(adjustl(adjustment))
      if (.not. allocated(model%data) .or. .not. allocated(model%residuals) .or. &
         .not. allocated(model%ar_order) .or. &
         .not. allocated(model%include_mean) .or. &
         .not. allocated(model%thresholds) .or. &
         (trim(selected_mode) /= 'pooled' .and. trim(selected_mode) /= 'regime') .or. &
         (trim(selected_adjustment) /= 'ols' .and. &
         trim(selected_adjustment) /= 'ml')) then
         out%info = 1
         return
      end if
      observations = size(model%data, 1)
      variables = size(model%data, 2)
      regimes = size(model%ar_order)
      maximum_order = max(1, maxval(model%ar_order))
      maximum_parameters = maxval(variables*model%ar_order + &
         merge(1, 0, model%include_mean))
      stored_parameters = 1 + variables*maximum_order
      start = max(maximum_order, model%delay) + 1
      usable = observations - start + 1
      if (variables < 1 .or. regimes < 2 .or. &
         size(model%residuals, 1) /= usable .or. &
         size(model%residuals, 2) /= variables .or. &
         size(model%thresholds) /= regimes - 1) then
         out%info = 1
         return
      end if
      if (present(threshold_variable)) then
         if (size(threshold_variable) /= observations .or. &
            .not. all(ieee_is_finite(threshold_variable))) then
            out%info = 2
            return
         end if
      end if
      allocate(gram(maximum_parameters, maximum_parameters, regimes), &
         residual_sum(variables, variables, regimes), &
         regressor(maximum_parameters), &
         out%coefficient_covariance(variables*maximum_parameters, &
         variables*maximum_parameters, regimes), &
         out%standard_error(variables, stored_parameters, regimes), &
         out%active(stored_parameters, regimes), &
         out%parameter_count(regimes), &
         out%pooled_residual_covariance(variables, variables), &
         out%regime_residual_covariance(variables, variables, regimes), &
         out%degrees_of_freedom(regimes))
      gram = 0.0_dp
      residual_sum = 0.0_dp
      out%coefficient_covariance = 0.0_dp
      out%standard_error = 0.0_dp
      out%active = .false.
      out%parameter_count = variables*model%ar_order + &
         merge(1, 0, model%include_mean)
      out%degrees_of_freedom = 0
      do row = 1, usable
         time = start + row - 1
         if (present(threshold_variable)) then
            threshold_value = threshold_variable(time - model%delay)
         else
            threshold_value = model%data(time - model%delay, &
               model%threshold_component)
         end if
         regime = inference_threshold_regime(threshold_value, model%thresholds)
         out%degrees_of_freedom(regime) = out%degrees_of_freedom(regime) + 1
         residual_sum(:, :, regime) = residual_sum(:, :, regime) + &
            outer_product(model%residuals(row, :), model%residuals(row, :))
         parameter_count = variables*model%ar_order(regime) + &
            merge(1, 0, model%include_mean(regime))
         regressor = 0.0_dp
         column = 0
         if (model%include_mean(regime)) then
            column = 1
            regressor(1) = 1.0_dp
         end if
         do lag = 1, model%ar_order(regime)
            do variable = 1, variables
               regressor(column + (lag - 1)*variables + variable) = &
                  model%data(time - lag, variable)
            end do
         end do
         gram(:parameter_count, :parameter_count, regime) = &
            gram(:parameter_count, :parameter_count, regime) + &
            outer_product(regressor(:parameter_count), regressor(:parameter_count))
      end do
      parameter_total = sum(variables*model%ar_order + &
         merge(1, 0, model%include_mean))
      denominator = usable
      if (trim(selected_adjustment) == 'ols') denominator = denominator - &
         parameter_total
      if (denominator <= 0) then
         out%info = 3
         return
      end if
      out%pooled_degrees_of_freedom = denominator
      out%pooled_residual_covariance = sum(residual_sum, dim=3)/ &
         real(denominator, dp)
      do regime = 1, regimes
         parameter_count = variables*model%ar_order(regime) + &
            merge(1, 0, model%include_mean(regime))
         denominator = out%degrees_of_freedom(regime)
         if (trim(selected_adjustment) == 'ols') denominator = denominator - &
            parameter_count
         if (denominator <= 0) then
            out%info = 3
            return
         end if
         out%degrees_of_freedom(regime) = denominator
         out%regime_residual_covariance(:, :, regime) = &
            residual_sum(:, :, regime)/real(denominator, dp)
         call invert_matrix(gram(:parameter_count, :parameter_count, regime), &
            inverse, status)
         if (status /= 0) then
            out%info = 20 + regime
            return
         end if
         allocate(selected_covariance(variables, variables))
         if (trim(selected_mode) == 'pooled') then
            selected_covariance = out%pooled_residual_covariance
         else
            selected_covariance = out%regime_residual_covariance(:, :, regime)
         end if
         compact_covariance = kronecker_product(selected_covariance, inverse)
         out%coefficient_covariance(:variables*parameter_count, &
            :variables*parameter_count, regime) = compact_covariance
         column = 0
         if (model%include_mean(regime)) then
            column = 1
            out%active(1, regime) = .true.
         end if
         do lag = 1, model%ar_order(regime)
            out%active(2 + (lag - 1)*variables:1 + lag*variables, regime) = .true.
         end do
         do equation = 1, variables
            if (model%include_mean(regime)) then
               compact_index = (equation - 1)*parameter_count + 1
               out%standard_error(equation, 1, regime) = sqrt(max(0.0_dp, &
                  compact_covariance(compact_index, compact_index)))
            end if
            do lag = 1, model%ar_order(regime)
               do variable = 1, variables
                  compact_index = (equation - 1)*parameter_count + column + &
                     (lag - 1)*variables + variable
                  out%standard_error(equation, &
                     1 + (lag - 1)*variables + variable, regime) = &
                     sqrt(max(0.0_dp, compact_covariance(compact_index, &
                     compact_index)))
               end do
            end do
         end do
         deallocate(inverse, selected_covariance, compact_covariance)
      end do
      out%variance_mode = selected_mode
      if (trim(selected_adjustment) == 'ols') then
         out%adjustment = 'OLS'
      else
         out%adjustment = 'ML'
      end if
      out%info = 0
   end function tsdyn_tvar_inference

   pure function tsdyn_tvecm_inference(model, variance_mode, adjustment) &
      result(out)
      !! Compute conditional TVECM covariance, standard errors, and tests.
      type(tsdyn_tvecm_model_t), intent(in) :: model !! Fitted TVECM model.
      character(len=*), intent(in), optional :: variance_mode !! Use pooled or regime residual covariance.
      character(len=*), intent(in), optional :: adjustment !! Denominator choice: OLS or ML.
      type(tsdyn_tvecm_inference_t) :: out
      real(dp), allocatable :: response(:, :), common(:, :), error_correction(:)
      real(dp), allocatable :: design(:, :), gram_inverse(:, :), bread(:, :)
      real(dp), allocatable :: residual_sum(:, :, :), selected_covariance(:, :)
      real(dp), allocatable :: meat(:, :), contribution(:, :), regressor_outer(:, :)
      real(dp), allocatable :: coefficient_vector(:)
      real(dp), allocatable :: cointegration(:)
      integer, allocatable :: regime(:)
      character(len=8) :: selected_mode, selected_adjustment
      real(dp) :: statistic, argument, nan_value
      integer :: observations, variables, regressors, regimes, rows
      integer :: deterministic_count, shared_columns
      integer :: row, selected_regime, denominator, regime_parameters
      integer :: equation, regressor, index, status

      selected_mode = 'pooled'
      if (present(variance_mode)) selected_mode = lowercase(adjustl(variance_mode))
      selected_adjustment = 'ols'
      if (present(adjustment)) selected_adjustment = &
         lowercase(adjustl(adjustment))
      if (model%info /= 0 .or. .not. allocated(model%data) .or. &
         .not. allocated(model%coefficients) .or. &
         .not. allocated(model%thresholds) .or. &
         .not. allocated(model%residuals) .or. &
         (trim(selected_mode) /= 'pooled' .and. &
         trim(selected_mode) /= 'regime') .or. &
         (trim(selected_adjustment) /= 'ols' .and. &
         trim(selected_adjustment) /= 'ml')) then
         out%info = 1
         return
      end if
      observations = size(model%data, 1)
      variables = size(model%data, 2)
      regimes = size(model%thresholds) + 1
      if (variables < 2 .or. model%lag < 1 .or. &
         (regimes /= 2 .and. regimes /= 3) .or. &
         observations <= model%lag + 1) then
         out%info = 1
         return
      end if
      cointegration = tvecm_cointegration_vector(model, variables)
      if (size(cointegration) /= variables) then
         out%info = 1
         return
      end if
      if (allocated(model%cointegration_offset)) then
         call tvecm_regression_data(model%data, model%lag, cointegration, &
            model%include_constant, model%include_trend, response, common, &
            error_correction, model%cointegration_offset)
      else
         call tvecm_regression_data(model%data, model%lag, cointegration, &
            model%include_constant, model%include_trend, response, common, &
            error_correction)
      end if
      deterministic_count = merge(1, 0, model%include_constant) + &
         merge(1, 0, model%include_trend)
      shared_columns = merge(deterministic_count, 0, &
         model%shared_deterministic .and. .not. model%only_error_correction)
      call tvecm_design(error_correction, common, model%thresholds, &
         model%only_error_correction, model%trim, design, regime, status, &
         shared_columns)
      if (status /= 0) then
         out%info = 2
         return
      end if
      rows = size(design, 1)
      regressors = size(design, 2)
      if (any(shape(model%coefficients) /= [variables, regressors]) .or. &
         any(shape(model%residuals) /= [rows, variables])) then
         out%info = 2
         return
      end if
      call invert_matrix(matmul(transpose(design), design), gram_inverse, status)
      if (status /= 0) then
         out%info = 3
         return
      end if
      allocate(residual_sum(variables, variables, regimes), &
         out%regime_observations(regimes), &
         out%regime_degrees_of_freedom(regimes), &
         out%pooled_residual_covariance(variables, variables), &
         out%regime_residual_covariance(variables, variables, regimes))
      residual_sum = 0.0_dp
      out%regime_observations = 0
      do row = 1, rows
         selected_regime = regime(row)
         out%regime_observations(selected_regime) = &
            out%regime_observations(selected_regime) + 1
         residual_sum(:, :, selected_regime) = &
            residual_sum(:, :, selected_regime) + &
            outer_product(model%residuals(row, :), model%residuals(row, :))
      end do
      denominator = rows
      if (trim(selected_adjustment) == 'ols') denominator = denominator - regressors
      if (denominator <= 0) then
         out%info = 4
         return
      end if
      out%pooled_degrees_of_freedom = denominator
      out%pooled_residual_covariance = sum(residual_sum, dim=3)/ &
         real(denominator, dp)
      regime_parameters = 1 + size(common, 2)
      do selected_regime = 1, regimes
         denominator = out%regime_observations(selected_regime)
         if (trim(selected_adjustment) == 'ols') denominator = &
            denominator - regime_parameters
         if (denominator <= 0) then
            out%info = 4
            return
         end if
         out%regime_degrees_of_freedom(selected_regime) = denominator
         out%regime_residual_covariance(:, :, selected_regime) = &
            residual_sum(:, :, selected_regime)/real(denominator, dp)
      end do
      if (trim(selected_mode) == 'pooled') then
         out%coefficient_covariance = kronecker_product( &
            out%pooled_residual_covariance, gram_inverse)
      else
         allocate(meat(variables*regressors, variables*regressors))
         meat = 0.0_dp
         do row = 1, rows
            regressor_outer = outer_product(design(row, :), design(row, :))
            selected_covariance = &
               out%regime_residual_covariance(:, :, regime(row))
            contribution = kronecker_product(selected_covariance, regressor_outer)
            meat = meat + contribution
         end do
         bread = kronecker_product(identity_matrix(variables), gram_inverse)
         out%coefficient_covariance = matmul(bread, matmul(meat, bread))
         out%coefficient_covariance = 0.5_dp*(out%coefficient_covariance + &
            transpose(out%coefficient_covariance))
      end if
      allocate(out%standard_error(variables, regressors), &
         out%t_statistic(variables, regressors), &
         out%p_value(variables, regressors), &
         coefficient_vector(variables*regressors))
      coefficient_vector = reshape(transpose(model%coefficients), &
         [variables*regressors])
      nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
      do equation = 1, variables
         do regressor = 1, regressors
            index = (equation - 1)*regressors + regressor
            out%standard_error(equation, regressor) = sqrt(max(0.0_dp, &
               out%coefficient_covariance(index, index)))
            if (out%standard_error(equation, regressor) <= tiny(1.0_dp)) then
               out%t_statistic(equation, regressor) = nan_value
               out%p_value(equation, regressor) = nan_value
            else
               statistic = coefficient_vector(index)/ &
                  out%standard_error(equation, regressor)
               out%t_statistic(equation, regressor) = statistic
               argument = real(out%pooled_degrees_of_freedom, dp)/ &
                  (real(out%pooled_degrees_of_freedom, dp) + statistic*statistic)
               out%p_value(equation, regressor) = regularized_beta(argument, &
                  0.5_dp*real(out%pooled_degrees_of_freedom, dp), 0.5_dp)
            end if
         end do
      end do
      out%parameter_count = variables*regressors
      out%variance_mode = selected_mode
      if (trim(selected_adjustment) == 'ols') then
         out%adjustment = 'OLS'
      else
         out%adjustment = 'ML'
      end if
      out%info = 0
   end function tsdyn_tvecm_inference

   pure integer function inference_threshold_regime(value, thresholds) &
      result(regime)
      !! Return the ordered threshold regime containing a scalar value.
      real(dp), intent(in) :: value !! Threshold-variable value.
      real(dp), intent(in) :: thresholds(:) !! Ordered regime thresholds.

      regime = 1
      do while (regime <= size(thresholds))
         if (value <= thresholds(regime)) exit
         regime = regime + 1
      end do
   end function inference_threshold_regime

   pure function tsdyn_rank_test(data, lag, deterministic, test_type, &
      null_rank, significance, fitted_rank) result(out)
      !! Test Johansen cointegration rank using Doornik gamma p-value approximations.
      real(dp), intent(in) :: data(:, :) !! Multivariate level observations.
      integer, intent(in) :: lag !! Number of lagged differences in the VECM.
      character(len=*), intent(in), optional :: deterministic !! Doornik case H_z, H_c, H_lc, H_l, or H_ql.
      character(len=*), intent(in), optional :: test_type !! Selection test: eigen or trace.
      integer, intent(in), optional :: null_rank !! Specific null rank; omit for automatic selection.
      real(dp), intent(in), optional :: significance !! Rejection significance level.
      integer, intent(in), optional :: fitted_rank !! Rank used for the small-sample parameter correction.
      type(tsdyn_rank_test_t) :: out
      type(tsdyn_johansen_kernel_t) :: kernel
      character(len=8) :: category, selected_test
      real(dp) :: selected_significance, lambda
      integer :: variables, rank, category_index, fitted, parameters
      integer :: adjusted_sample, selected_index
      logical :: unrestricted_constant, unrestricted_trend
      logical :: restricted_constant, restricted_trend

      category = 'h_lc'
      if (present(deterministic)) category = lowercase(adjustl(deterministic))
      selected_test = 'eigen'
      if (present(test_type)) selected_test = lowercase(adjustl(test_type))
      selected_significance = 0.05_dp
      if (present(significance)) selected_significance = significance
      call rank_category_configuration(category, category_index, &
         unrestricted_constant, unrestricted_trend, restricted_constant, &
         restricted_trend)
      variables = size(data, 2)
      fitted = min(1, max(0, variables - 1))
      if (present(fitted_rank)) fitted = fitted_rank
      if (category_index == 0 .or. &
         (trim(selected_test) /= 'eigen' .and. trim(selected_test) /= 'trace') .or. &
         selected_significance <= 0.0_dp .or. selected_significance >= 1.0_dp .or. &
         fitted < 0 .or. fitted > variables) then
         out%info = 1
         return
      end if
      kernel = johansen_rank_kernel(data, lag, unrestricted_constant, &
         unrestricted_trend, restricted_constant, restricted_trend)
      if (kernel%info /= 0) then
         out%info = 10 + kernel%info
         return
      end if
      allocate(out%eigenvalues(variables), out%trace_statistic(0:variables - 1), &
         out%trace_p_value(0:variables - 1), &
         out%adjusted_trace_p_value(0:variables - 1), &
         out%eigen_statistic(0:variables - 1), &
         out%eigen_p_value(0:variables - 1))
      out%eigenvalues = kernel%eigenvalues
      parameters = variables*(lag*variables + kernel%unrestricted_count) + &
         (2*variables + kernel%restricted_count)*fitted - fitted*fitted
      adjusted_sample = max(1, kernel%sample_size - parameters/variables)
      do rank = 0, variables - 1
         out%trace_statistic(rank) = 0.0_dp
         do selected_index = rank + 1, variables
            lambda = max(0.0_dp, min(1.0_dp - epsilon(1.0_dp), &
               kernel%eigenvalues(selected_index)))
            out%trace_statistic(rank) = out%trace_statistic(rank) - &
               real(kernel%sample_size, dp)*log(1.0_dp - lambda)
         end do
         lambda = max(0.0_dp, min(1.0_dp - epsilon(1.0_dp), &
            kernel%eigenvalues(rank + 1)))
         out%eigen_statistic(rank) = -real(kernel%sample_size, dp)* &
            log(1.0_dp - lambda)
         out%trace_p_value(rank) = doornik_rank_p_value( &
            out%trace_statistic(rank), variables - rank, category_index, &
            .true., .false., kernel%sample_size)
         out%adjusted_trace_p_value(rank) = doornik_rank_p_value( &
            out%trace_statistic(rank), variables - rank, category_index, &
            .true., .true., adjusted_sample)
         out%eigen_p_value(rank) = doornik_rank_p_value( &
            out%eigen_statistic(rank), variables - rank, category_index, &
            .false., .false., kernel%sample_size)
      end do
      out%selected_rank = variables
      if (present(null_rank)) then
         if (null_rank < 0 .or. null_rank >= variables) then
            out%info = 2
            return
         end if
         out%tested_rank = null_rank
         out%selected_rank = null_rank
         if (trim(selected_test) == 'trace') then
            out%selected_p_value = out%trace_p_value(null_rank)
         else
            out%selected_p_value = out%eigen_p_value(null_rank)
         end if
      else
         do rank = 0, variables - 1
            if (trim(selected_test) == 'trace') then
               out%selected_p_value = out%trace_p_value(rank)
            else
               out%selected_p_value = out%eigen_p_value(rank)
            end if
            if (out%selected_p_value > selected_significance) then
               out%selected_rank = rank
               exit
            end if
         end do
      end if
      out%sample_size = kernel%sample_size
      select case (category_index)
      case (1)
         out%deterministic = 'H_z'
      case (2)
         out%deterministic = 'H_c'
      case (3)
         out%deterministic = 'H_lc'
      case (4)
         out%deterministic = 'H_l'
      case (5)
         out%deterministic = 'H_ql'
      end select
      out%test_type = selected_test
      out%info = 0
   end function tsdyn_rank_test

   pure function tsdyn_rank_select(data, lag_max, rank_max, deterministic, &
      fit_measure, same_sample) result(out)
      !! Select cointegration rank and VECM lag jointly by information criteria.
      real(dp), intent(in) :: data(:, :) !! Multivariate level observations.
      integer, intent(in) :: lag_max !! Largest number of lagged differences.
      integer, intent(in), optional :: rank_max !! Largest cointegration rank to evaluate.
      character(len=*), intent(in), optional :: deterministic !! Unrestricted terms: none, const, trend, or both.
      character(len=*), intent(in), optional :: fit_measure !! Scoring base: SSR or LL.
      logical, intent(in), optional :: same_sample !! Whether every lag uses the same terminal sample length.
      type(tsdyn_rank_selection_t) :: out
      type(tsdyn_johansen_kernel_t) :: kernel
      real(dp), allocatable :: selected_data(:, :)
      character(len=8) :: terms, measure
      real(dp) :: fit, log_likelihood, log_covariance, lambda
      real(dp) :: penalty(3), best(3)
      integer :: variables, maximum_rank, lag, rank, start, criterion
      integer :: parameters, location_rank(3), location_lag(3)
      logical :: common_sample, unrestricted_constant, unrestricted_trend

      variables = size(data, 2)
      maximum_rank = variables
      if (present(rank_max)) maximum_rank = rank_max
      terms = 'const'
      if (present(deterministic)) terms = lowercase(adjustl(deterministic))
      measure = 'ssr'
      if (present(fit_measure)) measure = lowercase(adjustl(fit_measure))
      common_sample = .true.
      if (present(same_sample)) common_sample = same_sample
      unrestricted_constant = trim(terms) == 'const' .or. trim(terms) == 'both'
      unrestricted_trend = trim(terms) == 'trend' .or. trim(terms) == 'both'
      if (variables < 2 .or. lag_max < 1 .or. maximum_rank < 0 .or. &
         maximum_rank > variables .or. size(data, 1) <= lag_max + variables + 2 .or. &
         (trim(terms) /= 'none' .and. trim(terms) /= 'const' .and. &
         trim(terms) /= 'trend' .and. trim(terms) /= 'both') .or. &
         (trim(measure) /= 'ssr' .and. trim(measure) /= 'll') .or. &
         .not. all(ieee_is_finite(data))) then
         out%info = 1
         return
      end if
      allocate(out%aic(0:maximum_rank, 1:lag_max), &
         out%bic(0:maximum_rank, 1:lag_max), &
         out%hq(0:maximum_rank, 1:lag_max), &
         out%log_likelihood(0:maximum_rank, 1:lag_max))
      out%aic = huge(1.0_dp)
      out%bic = huge(1.0_dp)
      out%hq = huge(1.0_dp)
      out%log_likelihood = -huge(1.0_dp)
      do lag = 1, lag_max
         start = 1
         if (common_sample) start = lag_max - lag + 1
         selected_data = data(start:, :)
         kernel = johansen_rank_kernel(selected_data, lag, &
            unrestricted_constant, unrestricted_trend, .false., .false.)
         if (kernel%info /= 0) cycle
         penalty = [2.0_dp, log(real(kernel%sample_size, dp)), &
            2.0_dp*log(log(real(kernel%sample_size, dp)))]
         log_covariance = kernel%log_determinant_s00
         do rank = 0, maximum_rank
            if (rank > 0) then
               lambda = max(0.0_dp, min(1.0_dp - epsilon(1.0_dp), &
                  kernel%eigenvalues(rank)))
               log_covariance = log_covariance + log(1.0_dp - lambda)
            end if
            log_likelihood = -0.5_dp*real(kernel%sample_size*variables, dp)* &
               (log(2.0_dp*acos(-1.0_dp)) + 1.0_dp) - &
               0.5_dp*real(kernel%sample_size, dp)*log_covariance
            parameters = variables*(lag*variables + kernel%unrestricted_count) + &
               2*variables*rank - rank*rank
            if (trim(measure) == 'll') then
               fit = -2.0_dp*log_likelihood
            else
               fit = real(kernel%sample_size, dp)*log_covariance
            end if
            out%aic(rank, lag) = fit + penalty(1)*real(parameters, dp)
            out%bic(rank, lag) = fit + penalty(2)*real(parameters, dp)
            out%hq(rank, lag) = fit + penalty(3)*real(parameters, dp)
            out%log_likelihood(rank, lag) = log_likelihood
         end do
      end do
      best = huge(1.0_dp)
      location_rank = 0
      location_lag = 0
      do lag = 1, lag_max
         do rank = 0, maximum_rank
            do criterion = 1, 3
               select case (criterion)
               case (1)
                  fit = out%aic(rank, lag)
               case (2)
                  fit = out%bic(rank, lag)
               case default
                  fit = out%hq(rank, lag)
               end select
               if (fit < best(criterion)) then
                  best(criterion) = fit
                  location_rank(criterion) = rank
                  location_lag(criterion) = lag
               end if
            end do
         end do
      end do
      if (any(best >= huge(1.0_dp))) then
         out%info = 2
         return
      end if
      out%selected_rank = location_rank
      out%selected_lag = location_lag
      out%rank_max = maximum_rank
      out%lag_max = lag_max
      out%same_sample = common_sample
      out%deterministic = terms
      if (trim(measure) == 'll') then
         out%fit_measure = 'LL'
      else
         out%fit_measure = 'SSR'
      end if
      out%info = 0
   end function tsdyn_rank_select

   pure function johansen_rank_kernel(data, lag, unrestricted_constant, &
      unrestricted_trend, restricted_constant, restricted_trend) result(out)
      !! Compute the reduced-rank eigenvalues for a specified VECM design.
      real(dp), intent(in) :: data(:, :) !! Multivariate level observations.
      integer, intent(in) :: lag !! Number of lagged differences.
      logical, intent(in) :: unrestricted_constant !! Include a short-run constant.
      logical, intent(in) :: unrestricted_trend !! Include a short-run trend.
      logical, intent(in) :: restricted_constant !! Include a cointegrating-relation constant.
      logical, intent(in) :: restricted_trend !! Include a cointegrating-relation trend.
      type(tsdyn_johansen_kernel_t) :: out
      real(dp), allocatable :: differences(:, :), response(:, :), level(:, :)
      real(dp), allocatable :: regressors(:, :), r0(:, :), rk(:, :)
      real(dp), allocatable :: s00(:, :), s0k(:, :), skk(:, :)
      real(dp), allocatable :: s00_inverse(:, :), chol(:, :), chol_inverse(:, :)
      real(dp), allocatable :: eigenmatrix(:, :), eigenvectors(:, :)
      integer :: observations, variables, rows, unrestricted, restricted
      integer :: row, time, column, current_lag, status

      observations = size(data, 1)
      variables = size(data, 2)
      unrestricted = merge(1, 0, unrestricted_constant) + &
         merge(1, 0, unrestricted_trend)
      restricted = merge(1, 0, restricted_constant) + &
         merge(1, 0, restricted_trend)
      rows = observations - lag - 1
      if (variables < 2 .or. lag < 0 .or. rows <= variables + 1 .or. &
         .not. all(ieee_is_finite(data))) then
         out%info = 1
         return
      end if
      allocate(differences(observations - 1, variables), &
         response(rows, variables), level(rows, variables + restricted), &
         regressors(rows, lag*variables + unrestricted))
      differences = data(2:, :) - data(:observations - 1, :)
      regressors = 0.0_dp
      do row = 1, rows
         time = lag + row + 1
         response(row, :) = differences(time - 1, :)
         level(row, :variables) = data(time - 1, :)
         column = 0
         if (unrestricted_constant) then
            column = column + 1
            regressors(row, column) = 1.0_dp
         end if
         if (unrestricted_trend) then
            column = column + 1
            regressors(row, column) = real(row, dp)
         end if
         do current_lag = 1, lag
            regressors(row, column + 1:column + variables) = &
               differences(time - 1 - current_lag, :)
            column = column + variables
         end do
         column = variables
         if (restricted_constant) then
            column = column + 1
            level(row, column) = 1.0_dp
         end if
         if (restricted_trend) then
            column = column + 1
            level(row, column) = real(row, dp)
         end if
      end do
      call residualize_rank_data(response, regressors, r0, status)
      if (status /= 0) then
         out%info = 10 + status
         return
      end if
      call residualize_rank_data(level, regressors, rk, status)
      if (status /= 0) then
         out%info = 20 + status
         return
      end if
      s00 = matmul(transpose(r0), r0)/real(rows, dp)
      s0k = matmul(transpose(r0), rk)/real(rows, dp)
      skk = matmul(transpose(rk), rk)/real(rows, dp)
      allocate(s00_inverse(variables, variables))
      call inverse_logdet(s00, s00_inverse, out%log_determinant_s00, &
         status, 100.0_dp*epsilon(1.0_dp))
      if (status /= 0) then
         out%info = 30 + status
         return
      end if
      call cholesky_lower(skk, chol, status)
      if (status /= 0) then
         out%info = 40 + status
         return
      end if
      call invert_matrix(chol, chol_inverse, status)
      if (status /= 0) then
         out%info = 50 + status
         return
      end if
      eigenmatrix = matmul(chol_inverse, matmul(transpose(s0k), &
         matmul(s00_inverse, matmul(s0k, transpose(chol_inverse)))))
      eigenmatrix = 0.5_dp*(eigenmatrix + transpose(eigenmatrix))
      call symmetric_eigen(eigenmatrix, out%eigenvalues, eigenvectors, status)
      if (status /= 0) then
         out%info = 60 + status
         return
      end if
      out%eigenvalues = max(0.0_dp, min(1.0_dp - epsilon(1.0_dp), &
         out%eigenvalues(:variables)))
      out%sample_size = rows
      out%unrestricted_count = unrestricted
      out%restricted_count = restricted
      out%info = 0
   end function johansen_rank_kernel

   pure subroutine residualize_rank_data(values, regressors, residuals, info)
      !! Residualize each column against the common short-run regressors.
      real(dp), intent(in) :: values(:, :) !! Variables to residualize.
      real(dp), intent(in) :: regressors(:, :) !! Common regression design.
      real(dp), allocatable, intent(out) :: residuals(:, :) !! Residualized variables.
      integer, intent(out) :: info !! Status code; zero indicates success.
      real(dp), allocatable :: inverse(:, :), coefficients(:, :)

      if (size(values, 1) /= size(regressors, 1)) then
         info = 1
         return
      end if
      if (size(regressors, 2) == 0) then
         residuals = values
         info = 0
         return
      end if
      call invert_matrix(matmul(transpose(regressors), regressors), inverse, info)
      if (info /= 0) return
      coefficients = matmul(inverse, matmul(transpose(regressors), values))
      residuals = values - matmul(regressors, coefficients)
   end subroutine residualize_rank_data

   pure subroutine rank_category_configuration(category, category_index, &
      unrestricted_constant, unrestricted_trend, restricted_constant, &
      restricted_trend)
      !! Decode the five Johansen deterministic-term cases used by tsDyn.
      character(len=*), intent(in) :: category !! Doornik category name.
      integer, intent(out) :: category_index !! Table column index.
      logical, intent(out) :: unrestricted_constant !! Short-run constant flag.
      logical, intent(out) :: unrestricted_trend !! Short-run trend flag.
      logical, intent(out) :: restricted_constant !! Long-run constant flag.
      logical, intent(out) :: restricted_trend !! Long-run trend flag.

      category_index = 0
      unrestricted_constant = .false.
      unrestricted_trend = .false.
      restricted_constant = .false.
      restricted_trend = .false.
      select case (trim(category))
      case ('h_z', 'none')
         category_index = 1
      case ('h_c', 'restricted_constant')
         category_index = 2
         restricted_constant = .true.
      case ('h_lc', 'constant', 'const')
         category_index = 3
         unrestricted_constant = .true.
      case ('h_l', 'restricted_trend')
         category_index = 4
         unrestricted_constant = .true.
         restricted_trend = .true.
      case ('h_ql', 'both')
         category_index = 5
         unrestricted_constant = .true.
         unrestricted_trend = .true.
      end select
   end subroutine rank_category_configuration

   pure real(dp) function doornik_rank_p_value(statistic, nmp, category, &
      trace_test, small_sample, sample_size) result(probability)
      !! Evaluate Doornik's gamma approximation to a Johansen-test tail probability.
      real(dp), intent(in) :: statistic !! Johansen trace or maximum-eigenvalue statistic.
      integer, intent(in) :: nmp !! Number of remaining stochastic trends.
      integer, intent(in) :: category !! Deterministic-case table column.
      logical, intent(in) :: trace_test !! Whether the statistic is a trace statistic.
      logical, intent(in) :: small_sample !! Whether to apply the trace finite-sample correction.
      integer, intent(in) :: sample_size !! Effective sample size for the correction.
      real(dp), parameter :: trace_mean(6, 5) = reshape([ &
         2.0_dp, -1.0_dp, 0.07_dp, 0.07_dp, 0.0_dp, 0.0_dp, &
         2.0_dp, 2.01_dp, 0.0_dp, 0.06_dp, 0.05_dp, 0.0_dp, &
         2.0_dp, 1.05_dp, -1.55_dp, -0.5_dp, -0.23_dp, 0.0_dp, &
         2.0_dp, 4.05_dp, 0.5_dp, -0.23_dp, -0.07_dp, 0.0_dp, &
         2.0_dp, 2.85_dp, -5.1_dp, -0.1_dp, -0.06_dp, 1.35_dp], [6, 5])
      real(dp), parameter :: trace_variance(5, 5) = reshape([ &
         3.0_dp, -0.33_dp, -0.55_dp, 0.0_dp, 0.0_dp, &
         3.0_dp, 3.6_dp, 0.75_dp, -0.4_dp, -0.3_dp, &
         3.0_dp, 1.8_dp, 0.0_dp, -2.8_dp, -1.1_dp, &
         3.0_dp, 5.7_dp, 3.2_dp, -1.3_dp, -0.5_dp, &
         3.0_dp, 4.0_dp, 0.8_dp, -5.8_dp, -2.66_dp], [5, 5])
      real(dp), parameter :: eigen_mean(5, 5) = reshape([ &
         6.0019_dp, -2.7558_dp, 0.67185_dp, 0.1149_dp, -2.7764_dp, &
         5.9498_dp, 0.43402_dp, 0.04836_dp, 0.018198_dp, -2.3669_dp, &
         5.8271_dp, -1.6487_dp, -1.6118_dp, -0.25949_dp, -1.5666_dp, &
         5.8658_dp, 2.5595_dp, -0.34443_dp, -0.077991_dp, -1.7552_dp, &
         5.6364_dp, -0.90531_dp, -3.5166_dp, -0.47966_dp, -0.21447_dp], [5, 5])
      real(dp), parameter :: eigen_variance(5, 5) = reshape([ &
         1.8806_dp, -15.499_dp, 1.1136_dp, 0.070508_dp, 14.714_dp, &
         2.2231_dp, -7.9064_dp, 0.58592_dp, -0.034324_dp, 12.058_dp, &
         2.0785_dp, -9.7846_dp, -3.368_dp, -0.24528_dp, 13.074_dp, &
         1.9955_dp, -5.5428_dp, 1.2425_dp, 0.41949_dp, 12.841_dp, &
         2.0899_dp, -5.3303_dp, -7.1523_dp, -0.2526_dp, 12.393_dp], [5, 5])
      real(dp), parameter :: correction_mean(7, 5) = reshape([ &
         -0.101_dp, 0.499_dp, 0.896_dp, -0.562_dp, 0.00229_dp, 0.00662_dp, 0.0_dp, &
         0.0_dp, 0.465_dp, 0.984_dp, -0.273_dp, -244.0_dp, 0.0_dp, 0.0_dp, &
         0.134_dp, 0.422_dp, 1.02_dp, 2.17_dp, -0.00182_dp, 0.0_dp, -0.00321_dp, &
         0.0252_dp, 0.448_dp, 1.09_dp, -0.353_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
         -0.819_dp, 0.615_dp, 0.896_dp, 2.43_dp, 0.00149_dp, 0.0_dp, 0.0_dp], [7, 5])
      real(dp), parameter :: correction_variance(7, 5) = reshape([ &
         -0.204_dp, 0.98_dp, 3.11_dp, -2.14_dp, 0.0499_dp, -0.0103_dp, -0.00902_dp, &
         0.224_dp, 0.863_dp, 3.38_dp, -0.807_dp, 0.0_dp, 0.0_dp, -0.0091_dp, &
         0.422_dp, 0.734_dp, 3.76_dp, 4.32_dp, -0.00606_dp, 0.0_dp, -0.00718_dp, &
         0.0_dp, 0.836_dp, 3.99_dp, -1.33_dp, -0.00298_dp, -0.00139_dp, -0.00268_dp, &
         -1.29_dp, 1.01_dp, 3.92_dp, 4.67_dp, 0.00484_dp, -0.00127_dp, -0.0199_dp], [7, 5])
      real(dp) :: mean_value, variance_value, shape, rate, sample
      real(dp) :: values_mean(6), values_variance(5), correction_values(7)

      if (nmp < 1 .or. category < 1 .or. category > 5 .or. &
         sample_size < 1 .or. statistic < 0.0_dp) then
         probability = 0.0_dp
         return
      end if
      if (trace_test) then
         values_mean = [real(nmp*nmp, dp), real(nmp, dp), 1.0_dp, &
            merge(1.0_dp, 0.0_dp, nmp == 1), &
            merge(1.0_dp, 0.0_dp, nmp == 2), sqrt(real(nmp, dp))]
         values_variance = [real(nmp*nmp, dp), real(nmp, dp), 1.0_dp, &
            merge(1.0_dp, 0.0_dp, nmp == 1), &
            merge(1.0_dp, 0.0_dp, nmp == 2)]
         mean_value = dot_product(values_mean, trace_mean(:, category))
         variance_value = dot_product(values_variance, &
            trace_variance(:, category))
         if (small_sample) then
            sample = real(sample_size, dp)
            correction_values = [sqrt(real(nmp, dp))/sample, &
               real(nmp, dp)/sample, real(nmp*nmp, dp)/(sample*sample), &
               merge(1.0_dp/sample, 0.0_dp, nmp == 1), &
               merge(1.0_dp, 0.0_dp, nmp == 1), &
               merge(1.0_dp, 0.0_dp, nmp == 2), &
               merge(1.0_dp, 0.0_dp, nmp == 3)]
            mean_value = exp(log(mean_value) + dot_product(correction_values, &
               correction_mean(:, category)))
            variance_value = exp(log(variance_value) + &
               dot_product(correction_values, correction_variance(:, category)))
         end if
      else
         values_variance = [real(nmp, dp), 1.0_dp, &
            merge(1.0_dp, 0.0_dp, nmp == 1), &
            merge(1.0_dp, 0.0_dp, nmp == 2), sqrt(real(nmp, dp))]
         mean_value = dot_product(values_variance, eigen_mean(:, category))
         variance_value = dot_product(values_variance, &
            eigen_variance(:, category))
      end if
      if (mean_value <= 0.0_dp .or. variance_value <= 0.0_dp) then
         probability = 0.0_dp
         return
      end if
      shape = mean_value*mean_value/variance_value
      rate = mean_value/variance_value
      probability = regularized_gamma_q(shape, rate*statistic)
   end function doornik_rank_p_value

   pure function tsdyn_ar_regime_diagnostics(intercept, ar, &
      include_constant, include_trend, tolerance) result(out)
      !! Diagnose the equilibrium mean and characteristic roots of a linear AR model.
      real(dp), intent(in) :: intercept !! Model intercept when a constant is included.
      real(dp), intent(in) :: ar(:) !! Autoregressive coefficients in increasing lag order.
      logical, intent(in), optional :: include_constant !! Whether the model includes a constant.
      logical, intent(in), optional :: include_trend !! Whether the model includes a time trend.
      real(dp), intent(in), optional :: tolerance !! Unit-root and zero-denominator tolerance.
      type(tsdyn_regime_diagnostics_t) :: out
      real(dp), allocatable :: coefficients(:, :), constants(:)
      integer :: orders(1)
      logical :: has_constant(1), has_trend(1)
      character(len=16) :: labels(1)

      allocate(coefficients(1, size(ar)), constants(1))
      coefficients(1, :) = ar
      constants(1) = intercept
      orders = size(ar)
      has_constant = .true.
      if (present(include_constant)) has_constant = include_constant
      has_trend = .false.
      if (present(include_trend)) has_trend = include_trend
      labels = 'all'
      out = regime_diagnostics(constants, coefficients, orders, has_constant, &
         has_trend, labels, tolerance)
   end function tsdyn_ar_regime_diagnostics

   pure function tsdyn_setar_regime_diagnostics(model, tolerance) result(out)
      !! Diagnose regime equilibrium means and characteristic roots of a SETAR model.
      type(nts_tar_model_t), intent(in) :: model !! Fitted SETAR model.
      real(dp), intent(in), optional :: tolerance !! Unit-root and zero-denominator tolerance.
      type(tsdyn_regime_diagnostics_t) :: out
      real(dp), allocatable :: ar(:, :), intercept(:)
      logical, allocatable :: has_trend(:)
      character(len=16), allocatable :: labels(:)
      integer :: regimes, regime, order

      if (.not. allocated(model%coefficients) .or. &
         .not. allocated(model%ar_order) .or. &
         .not. allocated(model%include_mean)) then
         out%info = 1
         return
      end if
      regimes = size(model%coefficients, 1)
      if (regimes < 1 .or. size(model%ar_order) /= regimes .or. &
         size(model%include_mean) /= regimes .or. &
         any(model%ar_order < 1) .or. &
         size(model%coefficients, 2) < maxval(model%ar_order) + 1) then
         out%info = 1
         return
      end if
      allocate(ar(regimes, maxval(model%ar_order)), intercept(regimes), &
         has_trend(regimes), labels(regimes))
      ar = 0.0_dp
      intercept = model%coefficients(:, 1)
      has_trend = .false.
      call regime_labels(labels)
      do regime = 1, regimes
         order = model%ar_order(regime)
         ar(regime, :order) = model%coefficients(regime, 2:order + 1)
      end do
      out = regime_diagnostics(intercept, ar, model%ar_order, &
         model%include_mean, has_trend, labels, tolerance)
   end function tsdyn_setar_regime_diagnostics

   pure function tsdyn_lstar_regime_diagnostics(model, tolerance) result(out)
      !! Diagnose low- and high-regime equilibrium means and roots of an LSTAR model.
      type(tsdyn_lstar_model_t), intent(in) :: model !! Fitted LSTAR model.
      real(dp), intent(in), optional :: tolerance !! Unit-root and zero-denominator tolerance.
      type(tsdyn_regime_diagnostics_t) :: out
      real(dp), allocatable :: ar(:, :), intercept(:)
      integer :: orders(2), deterministic, maximum_order
      logical :: has_constant(2), has_trend(2)
      character(len=16) :: labels(2)

      if (.not. allocated(model%low_coefficients) .or. &
         .not. allocated(model%transition_coefficients) .or. &
         model%low_order < 1 .or. model%high_order < 1) then
         out%info = 1
         return
      end if
      deterministic = merge(1, 0, model%include_constant) + &
         merge(1, 0, model%include_trend)
      if (size(model%low_coefficients) < deterministic + model%low_order .or. &
         size(model%transition_coefficients) < deterministic + model%high_order) then
         out%info = 1
         return
      end if
      maximum_order = max(model%low_order, model%high_order)
      allocate(ar(2, maximum_order), intercept(2))
      ar = 0.0_dp
      intercept = 0.0_dp
      ar(1, :model%low_order) = &
         model%low_coefficients(deterministic + 1:deterministic + model%low_order)
      ar(2, :) = ar(1, :)
      ar(2, :model%high_order) = ar(2, :model%high_order) + &
         model%transition_coefficients(deterministic + 1:deterministic + model%high_order)
      if (model%include_constant) then
         intercept(1) = model%low_coefficients(1)
         intercept(2) = intercept(1) + model%transition_coefficients(1)
      end if
      orders = [model%low_order, maximum_order]
      has_constant = model%include_constant
      has_trend = model%include_trend
      labels = ['low             ', 'high            ']
      out = regime_diagnostics(intercept, ar, orders, has_constant, &
         has_trend, labels, tolerance)
   end function tsdyn_lstar_regime_diagnostics

   pure function tsdyn_star_regime_diagnostics(model, tolerance) result(out)
      !! Diagnose cumulative asymptotic regimes of an additive STAR model.
      type(tsdyn_star_model_t), intent(in) :: model !! Fitted additive STAR model.
      real(dp), intent(in), optional :: tolerance !! Unit-root and zero-denominator tolerance.
      type(tsdyn_regime_diagnostics_t) :: out
      real(dp), allocatable :: ar(:, :), intercept(:)
      integer, allocatable :: orders(:)
      logical, allocatable :: has_constant(:), has_trend(:)
      character(len=16), allocatable :: labels(:)
      integer :: deterministic, regime

      if (.not. allocated(model%coefficients) .or. model%regimes < 1 .or. &
         model%order < 1 .or. size(model%coefficients, 1) /= model%regimes) then
         out%info = 1
         return
      end if
      deterministic = merge(1, 0, model%include_constant) + &
         merge(1, 0, model%include_trend)
      if (size(model%coefficients, 2) < deterministic + model%order) then
         out%info = 1
         return
      end if
      allocate(ar(model%regimes, model%order), intercept(model%regimes), &
         orders(model%regimes), has_constant(model%regimes), &
         has_trend(model%regimes), labels(model%regimes))
      ar = 0.0_dp
      intercept = 0.0_dp
      call regime_labels(labels)
      do regime = 1, model%regimes
         ar(regime, :) = sum(model%coefficients(:regime, &
            deterministic + 1:deterministic + model%order), dim=1)
         if (model%include_constant) &
            intercept(regime) = sum(model%coefficients(:regime, 1))
      end do
      orders = model%order
      has_constant = model%include_constant
      has_trend = model%include_trend
      out = regime_diagnostics(intercept, ar, orders, has_constant, &
         has_trend, labels, tolerance)
   end function tsdyn_star_regime_diagnostics

   pure function regime_diagnostics(intercept, ar, orders, include_constant, &
      include_trend, labels, tolerance) result(out)
      !! Form common equilibrium and root diagnostics for AR coefficient regimes.
      real(dp), intent(in) :: intercept(:) !! Regime intercepts.
      real(dp), intent(in) :: ar(:, :) !! Regime-by-lag autoregressive coefficients.
      integer, intent(in) :: orders(:) !! Effective AR order of each regime.
      logical, intent(in) :: include_constant(:) !! Constant-inclusion flags.
      logical, intent(in) :: include_trend(:) !! Trend-inclusion flags.
      character(len=*), intent(in) :: labels(:) !! Regime labels.
      real(dp), intent(in), optional :: tolerance !! Unit-root and zero-denominator tolerance.
      type(tsdyn_regime_diagnostics_t) :: out
      type(polynomial_roots_t) :: root_result
      real(dp), allocatable :: polynomial(:)
      real(dp) :: selected_tolerance, denominator, nan_value
      integer :: regimes, maximum_order, regime, order

      regimes = size(ar, 1)
      selected_tolerance = sqrt(epsilon(1.0_dp))
      if (present(tolerance)) selected_tolerance = tolerance
      if (regimes < 1 .or. size(intercept) /= regimes .or. &
         size(orders) /= regimes .or. size(include_constant) /= regimes .or. &
         size(include_trend) /= regimes .or. size(labels) /= regimes .or. &
         any(orders < 1) .or. any(orders > size(ar, 2)) .or. &
         selected_tolerance < 0.0_dp) then
         out%info = 1
         return
      end if
      maximum_order = maxval(orders)
      nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
      allocate(out%regime_label(regimes), out%equilibrium_mean(regimes), &
         out%roots(maximum_order, regimes), &
         out%root_modulus(maximum_order, regimes), out%ar_order(regimes), &
         out%stable(regimes), out%unit_root(regimes), &
         out%mean_defined(regimes))
      out%regime_label = labels
      out%equilibrium_mean = nan_value
      out%roots = cmplx(nan_value, nan_value, dp)
      out%root_modulus = nan_value
      out%ar_order = orders
      out%stable = .false.
      out%unit_root = .false.
      out%mean_defined = .false.
      do regime = 1, regimes
         order = orders(regime)
         do while (order > 0)
            if (abs(ar(regime, order)) > selected_tolerance) exit
            order = order - 1
         end do
         out%ar_order(regime) = order
         allocate(polynomial(0:order))
         polynomial(0) = 1.0_dp
         if (order > 0) polynomial(1:) = -ar(regime, :order)
         root_result = polynomial_roots(polynomial)
         deallocate(polynomial)
         if (root_result%info /= 0 .or. size(root_result%roots) /= order) then
            out%info = 10 + root_result%info
            return
         end if
         out%roots(:order, regime) = root_result%roots
         out%root_modulus(:order, regime) = abs(root_result%roots)
         if (order > 0) then
            out%unit_root(regime) = any(abs(abs(root_result%roots) - 1.0_dp) <= &
               selected_tolerance)
            out%stable(regime) = all(abs(root_result%roots) > &
               1.0_dp + selected_tolerance)
         else
            out%stable(regime) = .true.
         end if
         if (.not. include_trend(regime)) then
            if (.not. include_constant(regime)) then
               out%equilibrium_mean(regime) = 0.0_dp
               out%mean_defined(regime) = .true.
            else
               denominator = 1.0_dp - sum(ar(regime, :orders(regime)))
               if (abs(denominator) > selected_tolerance) then
                  out%equilibrium_mean(regime) = intercept(regime)/denominator
                  out%mean_defined(regime) = .true.
               end if
            end if
         end if
      end do
      out%info = 0
   end function regime_diagnostics

   pure subroutine regime_labels(labels)
      !! Assign compact ordered labels to nonlinear autoregressive regimes.
      character(len=*), intent(out) :: labels(:) !! Labels to populate.
      character(len=16), parameter :: ordered(9) = [character(len=16) :: &
         'regime 1', 'regime 2', 'regime 3', 'regime 4', 'regime 5', &
         'regime 6', 'regime 7', 'regime 8', 'regime 9']
      integer :: regime

      do regime = 1, size(labels)
         if (regime <= size(ordered)) then
            labels(regime) = ordered(regime)
         else
            labels(regime) = 'regime'
         end if
      end do
      if (size(labels) == 2) labels = ['low             ', 'high            ']
      if (size(labels) == 3) labels = ['low             ', 'medium          ', &
         'high            ']
   end subroutine regime_labels

   pure real(dp) function tsdyn_delta(series, embedding_dimension, delay, &
      epsilon_value) result(value)
      !! Compute Manzan's correlation-integral conditional-independence statistic.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      integer, intent(in) :: embedding_dimension !! Conditioning embedding dimension.
      integer, intent(in) :: delay !! Delay between embedded coordinates.
      real(dp), intent(in) :: epsilon_value !! Correlation-integral length scale.
      real(dp) :: lower_integral, middle_integral, upper_integral

      value = ieee_value(0.0_dp, ieee_quiet_nan)
      if (embedding_dimension < 2 .or. delay < 1 .or. &
         epsilon_value <= 0.0_dp) return
      if (size(series) <= embedding_dimension*delay + 1) return
      lower_integral = delta_correlation_integral(series, &
         embedding_dimension - 1, delay, epsilon_value)
      middle_integral = delta_correlation_integral(series, &
         embedding_dimension, delay, epsilon_value)
      upper_integral = delta_correlation_integral(series, &
         embedding_dimension + 1, delay, epsilon_value)
      if (lower_integral <= 0.0_dp .or. upper_integral <= 0.0_dp) return
      value = 1.0_dp - middle_integral**2/(lower_integral*upper_integral)
   end function tsdyn_delta

   pure function tsdyn_delta_test(series, embedding_dimensions, delay, &
      epsilon_values, permutation_indices) result(out)
      !! Test conditional independence using caller-supplied random permutations.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      integer, intent(in) :: embedding_dimensions(:) !! Embedding dimensions to test.
      integer, intent(in) :: delay !! Delay between embedded coordinates.
      real(dp), intent(in) :: epsilon_values(:) !! Correlation-integral length scales.
      integer, intent(in) :: permutation_indices(:, :) !! Permuted observation indices by column.
      type(tsdyn_delta_test_t) :: out
      real(dp), allocatable :: permuted(:)
      integer :: dimension_index, epsilon_index, replication, exceedances
      integer :: first_index

      out%embedding_dimension = embedding_dimensions
      out%epsilon = epsilon_values
      allocate(out%statistic(size(embedding_dimensions), size(epsilon_values)))
      allocate(out%p_value(size(embedding_dimensions), size(epsilon_values)))
      allocate(out%permutation_statistic(size(permutation_indices, 2), &
         size(embedding_dimensions), size(epsilon_values)))
      out%statistic = ieee_value(0.0_dp, ieee_quiet_nan)
      out%p_value = 1.0_dp
      out%permutation_statistic = ieee_value(0.0_dp, ieee_quiet_nan)
      if (size(series) < 4 .or. size(embedding_dimensions) == 0 .or. &
         size(epsilon_values) == 0 .or. size(permutation_indices, 1) /= &
         size(series) .or. size(permutation_indices, 2) == 0 .or. &
         any(embedding_dimensions < 2) .or. delay < 1 .or. &
         any(epsilon_values <= 0.0_dp) .or. &
         any(permutation_indices < 1) .or. &
         any(permutation_indices > size(series))) then
         out%info = 1
         return
      end if
      do replication = 1, size(permutation_indices, 2)
         do first_index = 1, size(series) - 1
            if (any(permutation_indices(first_index + 1:, replication) == &
               permutation_indices(first_index, replication))) then
               out%info = 1
               return
            end if
         end do
      end do
      allocate(permuted(size(series)))
      do dimension_index = 1, size(embedding_dimensions)
         do epsilon_index = 1, size(epsilon_values)
            out%statistic(dimension_index, epsilon_index) = tsdyn_delta(series, &
               embedding_dimensions(dimension_index), delay, &
               epsilon_values(epsilon_index))
            do replication = 1, size(permutation_indices, 2)
               permuted = series(permutation_indices(:, replication))
               out%permutation_statistic(replication, dimension_index, &
                  epsilon_index) = tsdyn_delta(permuted, &
                  embedding_dimensions(dimension_index), delay, &
                  epsilon_values(epsilon_index))
            end do
            if (.not. ieee_is_finite(out%statistic(dimension_index, &
               epsilon_index)) .or. any(.not. ieee_is_finite( &
               out%permutation_statistic(:, dimension_index, &
               epsilon_index)))) then
               out%info = 2
               cycle
            end if
            exceedances = count(out%permutation_statistic(:, dimension_index, &
               epsilon_index) >= out%statistic(dimension_index, epsilon_index))
            out%p_value(dimension_index, epsilon_index) = &
               real(1 + exceedances, dp)/real(1 + size(permutation_indices, 2), dp)
         end do
      end do
   end function tsdyn_delta_test

   pure real(dp) function tsdyn_delta_linear(series, embedding_dimension, &
      delay) result(value)
      !! Compute the Gaussian linear benchmark for Manzan's delta statistic.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      integer, intent(in) :: embedding_dimension !! Conditioning embedding dimension.
      integer, intent(in) :: delay !! Delay between embedded coordinates.
      real(dp) :: first_eigenvalue, second_eigenvalue
      integer :: status

      value = ieee_value(0.0_dp, ieee_quiet_nan)
      if (embedding_dimension < 2 .or. delay < 1 .or. &
         size(series) <= embedding_dimension*delay) return
      call delta_largest_covariance_eigenvalue(series, embedding_dimension + 1, &
         delay, first_eigenvalue, status)
      if (status /= 0) return
      call delta_largest_covariance_eigenvalue(series, embedding_dimension, &
         delay, second_eigenvalue, status)
      if (status /= 0 .or. second_eigenvalue <= tiny(1.0_dp)) return
      value = 1.0_dp - first_eigenvalue/second_eigenvalue
   end function tsdyn_delta_linear

   pure function tsdyn_delta_linear_test(series, embedding_dimensions, delay, &
      epsilon_values, standard_normal_draws, maximum_ar_order) result(out)
      !! Test nonlinearity against an AIC-selected AR bootstrap null model.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      integer, intent(in) :: embedding_dimensions(:) !! Embedding dimensions to test.
      integer, intent(in) :: delay !! Delay between embedded coordinates.
      real(dp), intent(in) :: epsilon_values(:) !! Correlation-integral length scales.
      real(dp), intent(in) :: standard_normal_draws(:, :) !! Standard-normal bootstrap draws by column.
      integer, intent(in), optional :: maximum_ar_order !! Largest AR order considered under the null.
      type(tsdyn_delta_linear_test_t) :: out
      type(yule_walker_result_t) :: candidate, selected
      real(dp), allocatable :: centered(:), simulated(:), full_simulation(:)
      real(dp) :: candidate_delta, candidate_linear
      integer :: order_limit, order, dimension_index, epsilon_index
      integer :: replication, time, lag, burnin, exceedances

      out%embedding_dimension = embedding_dimensions
      out%epsilon = epsilon_values
      allocate(out%delta(size(embedding_dimensions), size(epsilon_values)))
      allocate(out%linear_delta(size(embedding_dimensions)))
      allocate(out%statistic(size(embedding_dimensions), size(epsilon_values)))
      allocate(out%p_value(size(embedding_dimensions), size(epsilon_values)))
      allocate(out%bootstrap_statistic(size(standard_normal_draws, 2), &
         size(embedding_dimensions), size(epsilon_values)))
      out%delta = ieee_value(0.0_dp, ieee_quiet_nan)
      out%linear_delta = ieee_value(0.0_dp, ieee_quiet_nan)
      out%statistic = ieee_value(0.0_dp, ieee_quiet_nan)
      out%bootstrap_statistic = ieee_value(0.0_dp, ieee_quiet_nan)
      out%p_value = 1.0_dp
      if (size(series) < 4 .or. size(embedding_dimensions) == 0 .or. &
         size(epsilon_values) == 0 .or. size(standard_normal_draws, 1) < &
         size(series) .or. size(standard_normal_draws, 2) == 0 .or. &
         any(embedding_dimensions < 2) .or. delay < 1 .or. &
         any(epsilon_values <= 0.0_dp)) then
         out%info = 1
         return
      end if
      order_limit = min(size(series) - 2, &
         max(0, int(10.0_dp*log10(real(size(series), dp)))))
      if (present(maximum_ar_order)) order_limit = min(order_limit, &
         max(0, maximum_ar_order))
      centered = series - sum(series)/real(size(series), dp)
      selected = yule_walker_fit(centered, 0)
      do order = 1, order_limit
         candidate = yule_walker_fit(centered, order)
         if (candidate%info == 0 .and. candidate%criterion < &
            selected%criterion) selected = candidate
      end do
      if (selected%info /= 0) then
         out%info = 2
         return
      end if
      out%ar_order = size(selected%coefficients)
      out%ar_coefficients = selected%coefficients
      out%innovation_variance = selected%variance
      do dimension_index = 1, size(embedding_dimensions)
         out%linear_delta(dimension_index) = tsdyn_delta_linear(series, &
            embedding_dimensions(dimension_index), delay)
         do epsilon_index = 1, size(epsilon_values)
            out%delta(dimension_index, epsilon_index) = tsdyn_delta(series, &
               embedding_dimensions(dimension_index), delay, &
               epsilon_values(epsilon_index))
            out%statistic(dimension_index, epsilon_index) = &
               out%delta(dimension_index, epsilon_index) - &
               out%linear_delta(dimension_index)
         end do
      end do
      burnin = size(standard_normal_draws, 1) - size(series)
      allocate(full_simulation(size(standard_normal_draws, 1)))
      allocate(simulated(size(series)))
      do replication = 1, size(standard_normal_draws, 2)
         full_simulation = 0.0_dp
         do time = 1, size(full_simulation)
            full_simulation(time) = sqrt(selected%variance)* &
               standard_normal_draws(time, replication)
            do lag = 1, min(out%ar_order, time - 1)
               full_simulation(time) = full_simulation(time) + &
                  selected%coefficients(lag)*full_simulation(time - lag)
            end do
         end do
         simulated = full_simulation(burnin + 1:)
         do dimension_index = 1, size(embedding_dimensions)
            candidate_linear = tsdyn_delta_linear(simulated, &
               embedding_dimensions(dimension_index), delay)
            do epsilon_index = 1, size(epsilon_values)
               candidate_delta = tsdyn_delta(simulated, &
                  embedding_dimensions(dimension_index), delay, &
                  epsilon_values(epsilon_index))
               out%bootstrap_statistic(replication, dimension_index, &
                  epsilon_index) = candidate_delta - candidate_linear
            end do
         end do
      end do
      if (any(.not. ieee_is_finite(out%statistic)) .or. &
         any(.not. ieee_is_finite(out%bootstrap_statistic))) then
         out%info = 3
         return
      end if
      do dimension_index = 1, size(embedding_dimensions)
         do epsilon_index = 1, size(epsilon_values)
            exceedances = count(out%bootstrap_statistic(:, dimension_index, &
               epsilon_index) >= out%statistic(dimension_index, epsilon_index))
            out%p_value(dimension_index, epsilon_index) = &
               real(1 + exceedances, dp)/ &
               real(1 + size(standard_normal_draws, 2), dp)
         end do
      end do
   end function tsdyn_delta_linear_test

   pure elemental real(dp) function tsdyn_logistic_transition(value, gamma, &
      threshold) result(weight)
      !! Evaluate the logistic smooth-transition function used by tsDyn.
      real(dp), intent(in) :: value !! Transition-variable value.
      real(dp), intent(in) :: gamma !! Positive transition slope.
      real(dp), intent(in) :: threshold !! Transition midpoint.
      real(dp) :: argument

      argument = max(-700.0_dp, min(700.0_dp, gamma*(value - threshold)))
      weight = 1.0_dp/(1.0_dp + exp(-argument))
   end function tsdyn_logistic_transition

   pure function tsdyn_setar_forecast_distribution(model, innovations, level, &
      zero_first) result(out)
      !! Simulate SETAR forecasts from supplied future innovation paths.
      type(nts_tar_model_t), intent(in) :: model !! Fitted SETAR model.
      real(dp), intent(in) :: innovations(:, :) !! Horizon-by-simulation innovations.
      real(dp), intent(in), optional :: level !! Central empirical interval coverage.
      logical, intent(in), optional :: zero_first !! Set each path's first innovation to zero.
      type(tsdyn_forecast_distribution_t) :: out
      real(dp), allocatable :: history(:)
      real(dp) :: selected_level, innovation
      integer :: observations, horizon, simulations, simulation, step, time
      integer :: regime, lag
      logical :: remove_first

      selected_level = 0.95_dp
      if (present(level)) selected_level = level
      remove_first = .false.
      if (present(zero_first)) remove_first = zero_first
      horizon = size(innovations, 1)
      simulations = size(innovations, 2)
      if (model%info /= 0 .or. .not. allocated(model%data) .or. horizon < 1 .or. &
         simulations < 1 .or. selected_level <= 0.0_dp .or. &
         selected_level >= 1.0_dp .or. .not. all(ieee_is_finite(innovations))) then
         out%info = 1
         return
      end if
      observations = size(model%data)
      allocate(out%paths(horizon, simulations), history(observations + horizon))
      do simulation = 1, simulations
         history(:observations) = model%data
         do step = 1, horizon
            time = observations + step
            regime = tvecm_regime(history(time - model%delay), model%thresholds)
            innovation = innovations(step, simulation)
            if (remove_first .and. step == 1) innovation = 0.0_dp
            history(time) = model%coefficients(regime, 1) + innovation
            do lag = 1, model%ar_order(regime)
               history(time) = history(time) + &
                  model%coefficients(regime, lag + 1)*history(time - lag)
            end do
            out%paths(step, simulation) = history(time)
         end do
      end do
      call summarize_forecast_distribution(out%paths, selected_level, &
         out%mean, out%standard_error, out%lower, out%upper)
      out%level = selected_level
      out%first_innovation_zero = remove_first
   end function tsdyn_setar_forecast_distribution

   pure function tsdyn_lstar_forecast_distribution(model, innovations, level, &
      zero_first, future_threshold) result(out)
      !! Simulate LSTAR forecasts from supplied future innovation paths.
      type(tsdyn_lstar_model_t), intent(in) :: model !! Fitted LSTAR model.
      real(dp), intent(in) :: innovations(:, :) !! Horizon-by-simulation innovations.
      real(dp), intent(in), optional :: level !! Central empirical interval coverage.
      logical, intent(in), optional :: zero_first !! Set each path's first innovation to zero.
      real(dp), intent(in), optional :: future_threshold(:) !! Future external transition values.
      type(tsdyn_forecast_distribution_t) :: out
      real(dp), allocatable :: history(:)
      real(dp) :: selected_level, innovation, transition_value, weight
      real(dp) :: low_mean, transition_mean, trend
      integer :: observations, horizon, simulations, simulation, step, time
      integer :: lag, position
      logical :: remove_first

      selected_level = 0.95_dp
      if (present(level)) selected_level = level
      remove_first = .false.
      if (present(zero_first)) remove_first = zero_first
      horizon = size(innovations, 1)
      simulations = size(innovations, 2)
      if (model%info /= 0 .or. .not. allocated(model%data) .or. horizon < 1 .or. &
         simulations < 1 .or. selected_level <= 0.0_dp .or. &
         selected_level >= 1.0_dp .or. .not. all(ieee_is_finite(innovations))) then
         out%info = 1
         return
      end if
      if (model%external_threshold .and. .not. present(future_threshold)) then
         out%info = 2
         return
      end if
      if (present(future_threshold)) then
         if (size(future_threshold) < horizon .or. &
            .not. all(ieee_is_finite(future_threshold(:horizon)))) then
            out%info = 2
            return
         end if
      end if
      observations = size(model%data)
      allocate(out%paths(horizon, simulations), history(observations + horizon))
      do simulation = 1, simulations
         history(:observations) = model%data
         do step = 1, horizon
            time = observations + step
            if (model%external_threshold) then
               transition_value = future_threshold(step)
            else
               transition_value = history(time - model%threshold_lag)
            end if
            weight = tsdyn_logistic_transition(transition_value, model%gamma, &
               model%threshold)
            position = 0
            low_mean = 0.0_dp
            transition_mean = 0.0_dp
            if (model%include_constant) then
               position = position + 1
               low_mean = model%low_coefficients(position)
               transition_mean = model%transition_coefficients(position)
            end if
            if (model%include_trend) then
               position = position + 1
               trend = real(time, dp)
               low_mean = low_mean + model%low_coefficients(position)*trend
               transition_mean = transition_mean + &
                  model%transition_coefficients(position)*trend
            end if
            do lag = 1, model%low_order
               low_mean = low_mean + model%low_coefficients(position + lag)* &
                  history(time - lag)
            end do
            do lag = 1, model%high_order
               transition_mean = transition_mean + &
                  model%transition_coefficients(position + lag)*history(time - lag)
            end do
            innovation = innovations(step, simulation)
            if (remove_first .and. step == 1) innovation = 0.0_dp
            history(time) = low_mean + weight*transition_mean + innovation
            out%paths(step, simulation) = history(time)
         end do
      end do
      call summarize_forecast_distribution(out%paths, selected_level, &
         out%mean, out%standard_error, out%lower, out%upper)
      out%level = selected_level
      out%first_innovation_zero = remove_first
   end function tsdyn_lstar_forecast_distribution

   pure function tsdyn_star_forecast_distribution(model, innovations, level, &
      zero_first) result(out)
      !! Simulate multi-regime STAR forecasts from supplied innovations.
      type(tsdyn_star_model_t), intent(in) :: model !! Fitted STAR model.
      real(dp), intent(in) :: innovations(:, :) !! Horizon-by-simulation innovations.
      real(dp), intent(in), optional :: level !! Central empirical interval coverage.
      logical, intent(in), optional :: zero_first !! Set each path's first innovation to zero.
      type(tsdyn_forecast_distribution_t) :: out
      real(dp), allocatable :: history(:), regressor(:)
      real(dp) :: selected_level, innovation, transition_value, prediction
      integer :: observations, horizon, simulations, simulation, step, time
      integer :: component, column, lag
      logical :: remove_first

      selected_level = 0.95_dp
      if (present(level)) selected_level = level
      remove_first = .false.
      if (present(zero_first)) remove_first = zero_first
      horizon = size(innovations, 1)
      simulations = size(innovations, 2)
      if (model%info /= 0 .or. .not. allocated(model%data) .or. horizon < 1 .or. &
         simulations < 1 .or. selected_level <= 0.0_dp .or. &
         selected_level >= 1.0_dp .or. .not. all(ieee_is_finite(innovations))) then
         out%info = 1
         return
      end if
      observations = size(model%data)
      allocate(out%paths(horizon, simulations), history(observations + horizon))
      allocate(regressor(size(model%coefficients, 2)))
      do simulation = 1, simulations
         history(:observations) = model%data
         do step = 1, horizon
            time = observations + step
            column = 0
            if (model%include_constant) then
               column = column + 1
               regressor(column) = 1.0_dp
            end if
            if (model%include_trend) then
               column = column + 1
               regressor(column) = real(time, dp)
            end if
            do lag = 1, model%order
               regressor(column + lag) = history(time - lag)
            end do
            transition_value = history(time - model%threshold_lag)
            prediction = dot_product(model%coefficients(1, :), regressor)
            do component = 1, model%regimes - 1
               prediction = prediction + tsdyn_logistic_transition( &
                  transition_value, model%gamma(component), &
                  model%threshold(component))*dot_product( &
                  model%coefficients(component + 1, :), regressor)
            end do
            innovation = innovations(step, simulation)
            if (remove_first .and. step == 1) innovation = 0.0_dp
            history(time) = prediction + innovation
            out%paths(step, simulation) = history(time)
         end do
      end do
      call summarize_forecast_distribution(out%paths, selected_level, &
         out%mean, out%standard_error, out%lower, out%upper)
      out%level = selected_level
      out%first_innovation_zero = remove_first
   end function tsdyn_star_forecast_distribution

   pure function tsdyn_llar_forecast_distribution(model, innovations, level, &
      zero_first, enlarge, expansion_percent) result(out)
      !! Simulate local-linear AR forecasts from supplied innovations.
      type(tsdyn_llar_model_t), intent(in) :: model !! Fitted local-linear AR model.
      real(dp), intent(in) :: innovations(:, :) !! Horizon-by-simulation innovations.
      real(dp), intent(in), optional :: level !! Central empirical interval coverage.
      logical, intent(in), optional :: zero_first !! Set each path's first innovation to zero.
      logical, intent(in), optional :: enlarge !! Enlarge empty neighborhoods when true.
      real(dp), intent(in), optional :: expansion_percent !! Radius increase percentage.
      type(tsdyn_forecast_distribution_t) :: out
      real(dp), allocatable :: history(:)
      real(dp) :: selected_level, innovation, epsilon_value, expansion, prediction
      integer :: observations, horizon, simulations, simulation, step, time
      integer :: neighbors, status
      logical :: remove_first, allow_enlargement

      selected_level = 0.95_dp
      if (present(level)) selected_level = level
      remove_first = .false.
      if (present(zero_first)) remove_first = zero_first
      allow_enlargement = .true.
      if (present(enlarge)) allow_enlargement = enlarge
      expansion = 20.0_dp
      if (present(expansion_percent)) expansion = expansion_percent
      horizon = size(innovations, 1)
      simulations = size(innovations, 2)
      if (model%info /= 0 .or. .not. allocated(model%data) .or. horizon < 1 .or. &
         simulations < 1 .or. selected_level <= 0.0_dp .or. &
         selected_level >= 1.0_dp .or. expansion <= 0.0_dp .or. &
         .not. all(ieee_is_finite(innovations))) then
         out%info = 1
         return
      end if
      observations = size(model%data)
      allocate(out%paths(horizon, simulations), history(observations + horizon))
      do simulation = 1, simulations
         history(:observations) = model%data
         do step = 1, horizon
            time = observations + step
            epsilon_value = model%selected_epsilon
            do
               call llar_next_prediction(history(:time - 1), model%order, &
                  model%delay, model%steps, epsilon_value, prediction, &
                  neighbors, status)
               if (status == 0 .or. .not. allow_enlargement) exit
               epsilon_value = epsilon_value*(1.0_dp + expansion/100.0_dp)
               if (epsilon_value > 1000.0_dp*(maxval(history(:time - 1)) - &
                  minval(history(:time - 1)))) exit
            end do
            if (status /= 0) then
               out%info = 2
               return
            end if
            innovation = innovations(step, simulation)
            if (remove_first .and. step == 1) innovation = 0.0_dp
            history(time) = prediction + innovation
            out%paths(step, simulation) = history(time)
         end do
      end do
      call summarize_forecast_distribution(out%paths, selected_level, &
         out%mean, out%standard_error, out%lower, out%upper)
      out%level = selected_level
      out%first_innovation_zero = remove_first
   end function tsdyn_llar_forecast_distribution

   pure function tsdyn_aar_forecast_distribution(model, innovations, level, &
      zero_first) result(out)
      !! Simulate additive-AR forecasts from supplied innovations.
      type(tsdyn_aar_model_t), intent(in) :: model !! Fitted additive AR model.
      real(dp), intent(in) :: innovations(:, :) !! Horizon-by-simulation innovations.
      real(dp), intent(in), optional :: level !! Central empirical interval coverage.
      logical, intent(in), optional :: zero_first !! Set each path's first innovation to zero.
      type(tsdyn_forecast_distribution_t) :: out
      real(dp), allocatable :: history(:), basis_value(:)
      real(dp) :: selected_level, innovation, prediction
      integer :: observations, horizon, simulations, simulation, step, time
      integer :: component, position
      logical :: remove_first

      selected_level = 0.95_dp
      if (present(level)) selected_level = level
      remove_first = .false.
      if (present(zero_first)) remove_first = zero_first
      horizon = size(innovations, 1)
      simulations = size(innovations, 2)
      if (model%info /= 0 .or. .not. allocated(model%data) .or. horizon < 1 .or. &
         simulations < 1 .or. selected_level <= 0.0_dp .or. &
         selected_level >= 1.0_dp .or. .not. all(ieee_is_finite(innovations))) then
         out%info = 1
         return
      end if
      observations = size(model%data)
      allocate(out%paths(horizon, simulations), history(observations + horizon))
      do simulation = 1, simulations
         history(:observations) = model%data
         do step = 1, horizon
            time = observations + step
            prediction = model%intercept
            do component = 1, model%order
               position = time - 1 - (component - 1)*model%delay
               basis_value = spline_basis_values(model%basis(component), &
                  history(position))
               prediction = prediction + dot_product(basis_value - &
                  model%basis_center(:, component), &
                  model%coefficients(:, component))
            end do
            innovation = innovations(step, simulation)
            if (remove_first .and. step == 1) innovation = 0.0_dp
            history(time) = prediction + innovation
            out%paths(step, simulation) = history(time)
         end do
      end do
      call summarize_forecast_distribution(out%paths, selected_level, &
         out%mean, out%standard_error, out%lower, out%upper)
      out%level = selected_level
      out%first_innovation_zero = remove_first
   end function tsdyn_aar_forecast_distribution

   pure function tsdyn_nnet_forecast_distribution(model, innovations, level, &
      zero_first) result(out)
      !! Simulate neural autoregression forecasts from supplied innovations.
      type(tsdyn_nnet_model_t), intent(in) :: model !! Fitted neural autoregression.
      real(dp), intent(in) :: innovations(:, :) !! Horizon-by-simulation innovations.
      real(dp), intent(in), optional :: level !! Central empirical interval coverage.
      logical, intent(in), optional :: zero_first !! Set each path's first innovation to zero.
      type(tsdyn_forecast_distribution_t) :: out
      real(dp), allocatable :: history(:), predictors(:, :), prediction(:)
      real(dp) :: selected_level, innovation
      integer :: observations, horizon, simulations, simulation, step, time
      integer :: component
      logical :: remove_first

      selected_level = 0.95_dp
      if (present(level)) selected_level = level
      remove_first = .false.
      if (present(zero_first)) remove_first = zero_first
      horizon = size(innovations, 1)
      simulations = size(innovations, 2)
      if (model%info /= 0 .or. model%steps /= 1 .or. &
         .not. allocated(model%data) .or. horizon < 1 .or. simulations < 1 .or. &
         selected_level <= 0.0_dp .or. selected_level >= 1.0_dp .or. &
         .not. all(ieee_is_finite(innovations))) then
         out%info = 1
         return
      end if
      observations = size(model%data)
      allocate(out%paths(horizon, simulations), history(observations + horizon))
      allocate(predictors(1, model%order))
      do simulation = 1, simulations
         history(:observations) = model%data
         do step = 1, horizon
            time = observations + step
            do component = 1, model%order
               predictors(1, component) = history(time - 1 - &
                  (component - 1)*model%delay)
            end do
            prediction = tsdyn_nnet_predict(model, predictors)
            if (size(prediction) /= 1 .or. .not. ieee_is_finite(prediction(1))) then
               out%info = 2
               return
            end if
            innovation = innovations(step, simulation)
            if (remove_first .and. step == 1) innovation = 0.0_dp
            history(time) = prediction(1) + innovation
            out%paths(step, simulation) = history(time)
         end do
      end do
      call summarize_forecast_distribution(out%paths, selected_level, &
         out%mean, out%standard_error, out%lower, out%upper)
      out%level = selected_level
      out%first_innovation_zero = remove_first
   end function tsdyn_nnet_forecast_distribution

   pure function tsdyn_setar_rolling_forecast(model, new_data, horizons, &
      refit_every) result(out)
      !! Evaluate rolling SETAR forecasts with optional periodic refitting.
      type(nts_tar_model_t), intent(in) :: model !! Fitted SETAR model.
      real(dp), intent(in) :: new_data(:) !! Sequential realized observations.
      integer, intent(in) :: horizons(:) !! Positive forecast horizons.
      integer, intent(in), optional :: refit_every !! Refit interval; zero disables refitting.
      type(rolling_forecast_result_t) :: out
      type(nts_tar_model_t) :: working, candidate, forecast_model
      type(tsdyn_forecast_distribution_t) :: distribution
      real(dp), allocatable :: full(:), zero_innovations(:, :)
      integer :: interval, horizon_index, target_index, horizon, origin, target

      interval = 0
      if (present(refit_every)) interval = refit_every
      if (.not. allocated(model%data)) then
         out%info = 2
         return
      end if
      call initialize_rolling_result(model%data, new_data, horizons, interval, out)
      if (out%info /= 0 .or. model%info /= 0) then
         if (out%info == 0) out%info = 2
         return
      end if
      full = [model%data, new_data]
      allocate(zero_innovations(maxval(horizons), 1))
      zero_innovations = 0.0_dp
      do horizon_index = 1, size(horizons)
         horizon = horizons(horizon_index)
         working = model
         do target_index = 1, size(new_data)
            origin = out%origin(target_index, horizon_index)
            target = out%target(target_index)
            if (origin <= max(maxval(model%ar_order), model%delay)) cycle
            if (rolling_refit_due(origin, out%training_size, interval)) then
               candidate = tsdyn_setar_fit(full(:origin), model%ar_order, &
                  size(model%thresholds), model%delay, &
                  include_mean=model%include_mean)
               if (candidate%info == 0) then
                  working = candidate
                  out%refitted(target_index, horizon_index) = .true.
               end if
            end if
            forecast_model = working
            forecast_model%data = full(:origin)
            distribution = tsdyn_setar_forecast_distribution(forecast_model, &
               zero_innovations(:horizon, :))
            if (distribution%info /= 0) cycle
            call record_rolling_forecast(out, target_index, horizon_index, &
               distribution%mean(horizon), full(target))
         end do
      end do
   end function tsdyn_setar_rolling_forecast

   pure function tsdyn_lstar_rolling_forecast(model, new_data, horizons, &
      refit_every, external_threshold) result(out)
      !! Evaluate rolling LSTAR forecasts with optional periodic refitting.
      type(tsdyn_lstar_model_t), intent(in) :: model !! Fitted LSTAR model.
      real(dp), intent(in) :: new_data(:) !! Sequential realized observations.
      integer, intent(in) :: horizons(:) !! Positive forecast horizons.
      integer, intent(in), optional :: refit_every !! Refit interval; zero disables refitting.
      real(dp), intent(in), optional :: external_threshold(:) !! Transition series aligned with training and new data.
      type(rolling_forecast_result_t) :: out
      type(tsdyn_lstar_model_t) :: working, candidate, forecast_model
      type(tsdyn_lstar_forecast_t) :: forecast_value
      real(dp), allocatable :: full(:)
      integer :: interval, horizon_index, target_index, horizon, origin, target

      interval = 0
      if (present(refit_every)) interval = refit_every
      if (.not. allocated(model%data)) then
         out%info = 2
         return
      end if
      call initialize_rolling_result(model%data, new_data, horizons, interval, out)
      if (out%info /= 0 .or. model%info /= 0) then
         if (out%info == 0) out%info = 2
         return
      end if
      full = [model%data, new_data]
      if (model%external_threshold) then
         if (.not. present(external_threshold)) then
            out%info = 3
            return
         end if
         if (size(external_threshold) /= size(full) .or. &
            .not. all(ieee_is_finite(external_threshold))) then
            out%info = 3
            return
         end if
      end if
      do horizon_index = 1, size(horizons)
         horizon = horizons(horizon_index)
         working = model
         do target_index = 1, size(new_data)
            origin = out%origin(target_index, horizon_index)
            target = out%target(target_index)
            if (origin <= max(model%threshold_lag, &
               max(model%low_order, model%high_order))) cycle
            if (rolling_refit_due(origin, out%training_size, interval)) then
               if (model%external_threshold) then
                  candidate = tsdyn_lstar_fit(full(:origin), model%low_order, &
                     model%high_order, model%threshold_lag, &
                     model%include_constant, model%include_trend, &
                     external_threshold(:origin))
               else
                  candidate = tsdyn_lstar_fit(full(:origin), model%low_order, &
                     model%high_order, model%threshold_lag, &
                     model%include_constant, model%include_trend)
               end if
               if (candidate%info == 0) then
                  working = candidate
                  out%refitted(target_index, horizon_index) = .true.
               end if
            end if
            forecast_model = working
            forecast_model%data = full(:origin)
            if (model%external_threshold) then
               forecast_value = tsdyn_lstar_forecast(forecast_model, horizon, &
                  external_threshold(origin + 1:target))
            else
               forecast_value = tsdyn_lstar_forecast(forecast_model, horizon)
            end if
            if (forecast_value%info /= 0) cycle
            call record_rolling_forecast(out, target_index, horizon_index, &
               forecast_value%mean(horizon), full(target))
         end do
      end do
   end function tsdyn_lstar_rolling_forecast

   pure function tsdyn_star_rolling_forecast(model, new_data, horizons, &
      refit_every) result(out)
      !! Evaluate rolling general-STAR forecasts with periodic refitting.
      type(tsdyn_star_model_t), intent(in) :: model !! Fitted STAR model.
      real(dp), intent(in) :: new_data(:) !! Sequential realized observations.
      integer, intent(in) :: horizons(:) !! Positive forecast horizons.
      integer, intent(in), optional :: refit_every !! Refit interval; zero disables refitting.
      type(rolling_forecast_result_t) :: out
      type(tsdyn_star_model_t) :: working, candidate, forecast_model
      type(tsdyn_star_forecast_t) :: forecast_value
      real(dp), allocatable :: full(:)
      integer :: interval, horizon_index, target_index, horizon, origin, target

      interval = 0
      if (present(refit_every)) interval = refit_every
      if (.not. allocated(model%data)) then
         out%info = 2
         return
      end if
      call initialize_rolling_result(model%data, new_data, horizons, interval, out)
      if (out%info /= 0 .or. model%info /= 0) then
         if (out%info == 0) out%info = 2
         return
      end if
      full = [model%data, new_data]
      do horizon_index = 1, size(horizons)
         horizon = horizons(horizon_index)
         working = model
         do target_index = 1, size(new_data)
            origin = out%origin(target_index, horizon_index)
            target = out%target(target_index)
            if (origin <= max(model%order, model%threshold_lag)) cycle
            if (rolling_refit_due(origin, out%training_size, interval)) then
               candidate = tsdyn_star_fit(full(:origin), model%order, &
                  model%regimes, model%threshold_lag, model%include_constant, &
                  model%include_trend)
               if (candidate%info == 0) then
                  working = candidate
                  out%refitted(target_index, horizon_index) = .true.
               end if
            end if
            forecast_model = working
            forecast_model%data = full(:origin)
            forecast_value = tsdyn_star_forecast(forecast_model, horizon)
            if (forecast_value%info /= 0) cycle
            call record_rolling_forecast(out, target_index, horizon_index, &
               forecast_value%mean(horizon), full(target))
         end do
      end do
   end function tsdyn_star_rolling_forecast

   pure function tsdyn_llar_rolling_forecast(model, new_data, horizons, &
      refit_every) result(out)
      !! Evaluate rolling local-linear AR forecasts with periodic refitting.
      type(tsdyn_llar_model_t), intent(in) :: model !! Fitted local-linear AR model.
      real(dp), intent(in) :: new_data(:) !! Sequential realized observations.
      integer, intent(in) :: horizons(:) !! Positive forecast horizons.
      integer, intent(in), optional :: refit_every !! Refit interval; zero disables refitting.
      type(rolling_forecast_result_t) :: out
      type(tsdyn_llar_model_t) :: working, candidate, forecast_model
      type(tsdyn_llar_forecast_t) :: forecast_value
      real(dp), allocatable :: full(:)
      integer :: interval, horizon_index, target_index, horizon, origin, target

      interval = 0
      if (present(refit_every)) interval = refit_every
      if (.not. allocated(model%data)) then
         out%info = 2
         return
      end if
      call initialize_rolling_result(model%data, new_data, horizons, interval, out)
      if (out%info /= 0 .or. model%info /= 0) then
         if (out%info == 0) out%info = 2
         return
      end if
      full = [model%data, new_data]
      do horizon_index = 1, size(horizons)
         horizon = horizons(horizon_index)
         working = model
         do target_index = 1, size(new_data)
            origin = out%origin(target_index, horizon_index)
            target = out%target(target_index)
            if (origin <= (model%order - 1)*model%delay + model%steps) cycle
            if (rolling_refit_due(origin, out%training_size, interval)) then
               candidate = tsdyn_llar_fit(full(:origin), model%order, &
                  model%delay, model%steps, model%epsilon)
               if (candidate%info == 0) then
                  working = candidate
                  out%refitted(target_index, horizon_index) = .true.
               end if
            end if
            forecast_model = working
            forecast_model%data = full(:origin)
            forecast_value = tsdyn_llar_forecast(forecast_model, horizon)
            if (forecast_value%info /= 0) cycle
            call record_rolling_forecast(out, target_index, horizon_index, &
               forecast_value%mean(horizon), full(target))
         end do
      end do
   end function tsdyn_llar_rolling_forecast

   pure function tsdyn_aar_rolling_forecast(model, new_data, horizons, &
      refit_every) result(out)
      !! Evaluate rolling additive-AR forecasts with periodic refitting.
      type(tsdyn_aar_model_t), intent(in) :: model !! Fitted additive AR model.
      real(dp), intent(in) :: new_data(:) !! Sequential realized observations.
      integer, intent(in) :: horizons(:) !! Positive forecast horizons.
      integer, intent(in), optional :: refit_every !! Refit interval; zero disables refitting.
      type(rolling_forecast_result_t) :: out
      type(tsdyn_aar_model_t) :: working, candidate, forecast_model
      type(tsdyn_aar_forecast_t) :: forecast_value
      real(dp), allocatable :: full(:)
      integer :: interval, horizon_index, target_index, horizon, origin, target

      interval = 0
      if (present(refit_every)) interval = refit_every
      if (.not. allocated(model%data)) then
         out%info = 2
         return
      end if
      call initialize_rolling_result(model%data, new_data, horizons, interval, out)
      if (out%info /= 0 .or. model%info /= 0) then
         if (out%info == 0) out%info = 2
         return
      end if
      full = [model%data, new_data]
      do horizon_index = 1, size(horizons)
         horizon = horizons(horizon_index)
         working = model
         do target_index = 1, size(new_data)
            origin = out%origin(target_index, horizon_index)
            target = out%target(target_index)
            if (origin <= (model%order - 1)*model%delay + model%steps) cycle
            if (rolling_refit_due(origin, out%training_size, interval)) then
               candidate = tsdyn_aar_fit(full(:origin), model%order, model%delay, &
                  model%steps, model%basis_count, model%lambda)
               if (candidate%info == 0) then
                  working = candidate
                  out%refitted(target_index, horizon_index) = .true.
               end if
            end if
            forecast_model = working
            forecast_model%data = full(:origin)
            forecast_value = tsdyn_aar_forecast(forecast_model, horizon)
            if (forecast_value%info /= 0) cycle
            call record_rolling_forecast(out, target_index, horizon_index, &
               forecast_value%mean(horizon), full(target))
         end do
      end do
   end function tsdyn_aar_rolling_forecast

   pure function tsdyn_nnet_rolling_forecast(model, new_data, horizons, &
      refit_every) result(out)
      !! Evaluate rolling neural-AR forecasts with periodic refitting.
      type(tsdyn_nnet_model_t), intent(in) :: model !! Fitted neural autoregression.
      real(dp), intent(in) :: new_data(:) !! Sequential realized observations.
      integer, intent(in) :: horizons(:) !! Positive forecast horizons.
      integer, intent(in), optional :: refit_every !! Refit interval; zero disables refitting.
      type(rolling_forecast_result_t) :: out
      type(tsdyn_nnet_model_t) :: working, candidate, forecast_model
      type(tsdyn_nnet_forecast_t) :: forecast_value
      real(dp), allocatable :: full(:)
      integer :: interval, horizon_index, target_index, horizon, origin, target

      interval = 0
      if (present(refit_every)) interval = refit_every
      if (.not. allocated(model%data)) then
         out%info = 2
         return
      end if
      call initialize_rolling_result(model%data, new_data, horizons, interval, out)
      if (out%info /= 0 .or. model%info /= 0 .or. model%steps /= 1) then
         if (out%info == 0) out%info = 2
         return
      end if
      full = [model%data, new_data]
      do horizon_index = 1, size(horizons)
         horizon = horizons(horizon_index)
         working = model
         do target_index = 1, size(new_data)
            origin = out%origin(target_index, horizon_index)
            target = out%target(target_index)
            if (origin <= (model%order - 1)*model%delay + 1) cycle
            if (rolling_refit_due(origin, out%training_size, interval)) then
               candidate = tsdyn_nnet_fit(full(:origin), model%order, &
                  model%hidden_count, model%delay, model%steps, decay= &
                  model%network%decay)
               if (candidate%info == 0) then
                  working = candidate
                  out%refitted(target_index, horizon_index) = .true.
               end if
            end if
            forecast_model = working
            forecast_model%data = full(:origin)
            forecast_value = tsdyn_nnet_forecast(forecast_model, horizon)
            if (forecast_value%info /= 0) cycle
            call record_rolling_forecast(out, target_index, horizon_index, &
               forecast_value%mean(horizon), full(target))
         end do
      end do
   end function tsdyn_nnet_rolling_forecast

   pure function tsdyn_lstar_simulate_from_innovations(low_intercept, low_ar, &
      transition_intercept, transition_ar, gamma, threshold, threshold_lag, &
      innovations, burnin) result(out)
      !! Simulate an LSTAR process from supplied innovations.
      real(dp), intent(in) :: low_intercept !! Baseline-regime intercept.
      real(dp), intent(in) :: low_ar(:) !! Baseline autoregressive coefficients.
      real(dp), intent(in) :: transition_intercept !! Smooth-transition intercept increment.
      real(dp), intent(in) :: transition_ar(:) !! Smooth-transition AR increments.
      real(dp), intent(in) :: gamma !! Positive logistic transition slope.
      real(dp), intent(in) :: threshold !! Logistic transition midpoint.
      integer, intent(in) :: threshold_lag !! Positive lag of the transition variable.
      real(dp), intent(in) :: innovations(:) !! Supplied process innovations.
      integer, intent(in), optional :: burnin !! Number of leading values to discard.
      type(tsdyn_lstar_simulation_t) :: out
      real(dp), allocatable :: complete(:), weight(:)
      real(dp) :: low_mean, transition_mean
      integer :: total, memory, selected_burnin, time, lag

      total = size(innovations)
      selected_burnin = 0
      if (present(burnin)) selected_burnin = burnin
      memory = max(threshold_lag, max(size(low_ar), size(transition_ar)))
      if (total <= memory .or. selected_burnin < 0 .or. selected_burnin >= total .or. &
         threshold_lag < 1 .or. gamma <= 0.0_dp .or. &
         .not. all(ieee_is_finite(low_ar)) .or. &
         .not. all(ieee_is_finite(transition_ar)) .or. &
         .not. all(ieee_is_finite(innovations)) .or. &
         .not. ieee_is_finite(low_intercept) .or. &
         .not. ieee_is_finite(transition_intercept) .or. &
         .not. ieee_is_finite(gamma) .or. .not. ieee_is_finite(threshold)) then
         out%info = 1
         return
      end if
      allocate(complete(total), weight(total))
      complete(:memory) = innovations(:memory)
      weight(:memory) = tsdyn_logistic_transition(complete(:memory), gamma, threshold)
      do time = memory + 1, total
         low_mean = low_intercept
         do lag = 1, size(low_ar)
            low_mean = low_mean + low_ar(lag)*complete(time - lag)
         end do
         transition_mean = transition_intercept
         do lag = 1, size(transition_ar)
            transition_mean = transition_mean + &
               transition_ar(lag)*complete(time - lag)
         end do
         weight(time) = tsdyn_logistic_transition( &
            complete(time - threshold_lag), gamma, threshold)
         complete(time) = low_mean + weight(time)*transition_mean + innovations(time)
      end do
      out%series = complete(selected_burnin + 1:)
      out%innovations = innovations(selected_burnin + 1:)
      out%transition_weight = weight(selected_burnin + 1:)
      out%low_coefficients = [low_intercept, low_ar]
      out%transition_coefficients = [transition_intercept, transition_ar]
      out%gamma = gamma
      out%threshold = threshold
      out%threshold_lag = threshold_lag
      out%burnin = selected_burnin
   end function tsdyn_lstar_simulate_from_innovations

   function tsdyn_lstar_simulate(low_intercept, low_ar, transition_intercept, &
      transition_ar, gamma, threshold, threshold_lag, innovation_sd, &
      observations, burnin) result(out)
      !! Simulate an LSTAR process using the shared normal random stream.
      real(dp), intent(in) :: low_intercept !! Baseline-regime intercept.
      real(dp), intent(in) :: low_ar(:) !! Baseline autoregressive coefficients.
      real(dp), intent(in) :: transition_intercept !! Smooth-transition intercept increment.
      real(dp), intent(in) :: transition_ar(:) !! Smooth-transition AR increments.
      real(dp), intent(in) :: gamma !! Positive logistic transition slope.
      real(dp), intent(in) :: threshold !! Logistic transition midpoint.
      integer, intent(in) :: threshold_lag !! Positive lag of the transition variable.
      real(dp), intent(in) :: innovation_sd !! Positive innovation standard deviation.
      integer, intent(in) :: observations !! Number of retained observations.
      integer, intent(in), optional :: burnin !! Number of leading values to discard.
      type(tsdyn_lstar_simulation_t) :: out
      real(dp), allocatable :: innovations(:)
      integer :: selected_burnin, time

      selected_burnin = 500
      if (present(burnin)) selected_burnin = burnin
      if (observations < 1 .or. selected_burnin < 0 .or. innovation_sd <= 0.0_dp) then
         out%info = 1
         return
      end if
      allocate(innovations(observations + selected_burnin))
      do time = 1, size(innovations)
         innovations(time) = innovation_sd*random_standard_normal()
      end do
      out = tsdyn_lstar_simulate_from_innovations(low_intercept, low_ar, &
         transition_intercept, transition_ar, gamma, threshold, threshold_lag, &
         innovations, selected_burnin)
   end function tsdyn_lstar_simulate

   pure function tsdyn_lstar_fit(series, low_order, high_order, threshold_lag, &
      include_constant, include_trend, threshold_variable) result(model)
      !! Fit a two-regime logistic STAR model by variable-projection least squares.
      real(dp), intent(in) :: series(:) !! Univariate time series.
      integer, intent(in) :: low_order !! Baseline-regime AR order.
      integer, intent(in) :: high_order !! Smooth-transition AR order.
      integer, intent(in), optional :: threshold_lag !! Positive internal transition lag.
      logical, intent(in), optional :: include_constant !! Whether each component has an intercept.
      logical, intent(in), optional :: include_trend !! Whether each component has a trend.
      real(dp), intent(in), optional :: threshold_variable(:) !! External transition variable.
      type(tsdyn_lstar_model_t) :: model
      real(dp), allocatable :: response(:), low_design(:, :), high_design(:, :)
      real(dp), allocatable :: transition_data(:), ordered(:), coefficient(:)
      real(dp), allocatable :: fitted(:), residual(:), weight(:)
      real(dp) :: lower_threshold, upper_threshold, gamma, threshold
      real(dp) :: best_gamma, best_threshold, rss, best_rss
      real(dp) :: log_gamma_step, threshold_step, candidate_log_gamma
      integer :: selected_lag, start, rows, gamma_index, threshold_index
      integer :: iteration, direction_gamma, direction_threshold, info
      logical :: selected_constant, selected_trend, improved

      selected_lag = 1
      selected_constant = .true.
      selected_trend = .false.
      if (present(threshold_lag)) selected_lag = threshold_lag
      if (present(include_constant)) selected_constant = include_constant
      if (present(include_trend)) selected_trend = include_trend
      start = max(selected_lag, max(low_order, high_order)) + 1
      if (low_order < 0 .or. high_order < 0 .or. selected_lag < 1 .or. &
         start >= size(series) .or. .not. all(ieee_is_finite(series))) then
         model%info = 1
         return
      end if
      if (present(threshold_variable)) then
         if (size(threshold_variable) /= size(series) .or. &
            .not. all(ieee_is_finite(threshold_variable))) then
            model%info = 1
            return
         end if
      end if
      call lstar_regression_data(series, low_order, high_order, selected_lag, &
         selected_constant, selected_trend, threshold_variable, response, &
         low_design, high_design, transition_data)
      rows = size(response)
      if (rows <= size(low_design, 2) + size(high_design, 2) + 2) then
         model%info = 1
         return
      end if
      ordered = sorted(transition_data)
      lower_threshold = quantile(ordered, 0.10_dp)
      upper_threshold = quantile(ordered, 0.90_dp)
      if (upper_threshold <= lower_threshold) then
         model%info = 2
         return
      end if
      best_rss = huge(1.0_dp)
      best_gamma = 1.0_dp
      best_threshold = 0.5_dp*(lower_threshold + upper_threshold)
      do gamma_index = 0, 19
         gamma = exp(log(0.1_dp) + real(gamma_index, dp)* &
            (log(100.0_dp) - log(0.1_dp))/19.0_dp)
         do threshold_index = 0, 59
            threshold = lower_threshold + real(threshold_index, dp)* &
               (upper_threshold - lower_threshold)/59.0_dp
            call lstar_profile_fit(response, low_design, high_design, &
               transition_data, gamma, threshold, coefficient, fitted, residual, &
               weight, rss, info)
            if (info == 0 .and. rss < best_rss) then
               best_rss = rss
               best_gamma = gamma
               best_threshold = threshold
            end if
         end do
      end do
      log_gamma_step = 0.25_dp*(log(100.0_dp) - log(0.1_dp))
      threshold_step = 0.10_dp*(upper_threshold - lower_threshold)
      do iteration = 1, 100
         improved = .false.
         do direction_gamma = -1, 1
            do direction_threshold = -1, 1
               if (direction_gamma == 0 .and. direction_threshold == 0) cycle
               candidate_log_gamma = log(best_gamma) + &
                  real(direction_gamma, dp)*log_gamma_step
               gamma = exp(max(log(1.0e-4_dp), min(log(1.0e4_dp), &
                  candidate_log_gamma)))
               threshold = best_threshold + &
                  real(direction_threshold, dp)*threshold_step
               call lstar_profile_fit(response, low_design, high_design, &
                  transition_data, gamma, threshold, coefficient, fitted, &
                  residual, weight, rss, info)
               if (info == 0 .and. rss < best_rss) then
                  best_rss = rss
                  best_gamma = gamma
                  best_threshold = threshold
                  improved = .true.
               end if
            end do
         end do
         if (.not. improved) then
            log_gamma_step = 0.5_dp*log_gamma_step
            threshold_step = 0.5_dp*threshold_step
         end if
         if (max(log_gamma_step, threshold_step) < 1.0e-6_dp) exit
      end do
      call lstar_profile_fit(response, low_design, high_design, transition_data, &
         best_gamma, best_threshold, coefficient, fitted, residual, weight, &
         rss, info)
      if (info /= 0) then
         model%info = info
         return
      end if
      model%data = series
      model%low_coefficients = coefficient(:size(low_design, 2))
      model%transition_coefficients = coefficient(size(low_design, 2) + 1:)
      model%fitted = fitted
      model%residuals = residual
      model%transition_weight = weight
      model%threshold_data = transition_data
      model%gamma = best_gamma
      model%threshold = best_threshold
      model%rss = rss
      model%low_order = low_order
      model%high_order = high_order
      model%threshold_lag = selected_lag
      model%first_fitted = start
      model%include_constant = selected_constant
      model%include_trend = selected_trend
      model%external_threshold = present(threshold_variable)
      model%iterations = iteration
      model%converged = max(log_gamma_step, threshold_step) < 1.0e-6_dp
      model%residual_standard_deviation = sqrt(rss/real(rows - size(coefficient), dp))
      model%aic = real(rows, dp)*log(rss/real(rows, dp)) + &
         2.0_dp*real(size(coefficient) + 2, dp)
      model%bic = real(rows, dp)*log(rss/real(rows, dp)) + &
         log(real(rows, dp))*real(size(coefficient) + 2, dp)
      call lstar_inference(response, low_design, high_design, transition_data, &
         coefficient, best_gamma, best_threshold, &
         model%residual_standard_deviation**2, model%covariance, &
         model%standard_error, model%confidence_lower, model%confidence_upper)
   end function tsdyn_lstar_fit

   pure function tsdyn_lstar_forecast(model, horizon, future_threshold) result(out)
      !! Compute recursive naive forecasts from a fitted LSTAR model.
      type(tsdyn_lstar_model_t), intent(in) :: model !! Fitted LSTAR model.
      integer, intent(in) :: horizon !! Positive forecast horizon.
      real(dp), intent(in), optional :: future_threshold(:) !! Future external transition values.
      type(tsdyn_lstar_forecast_t) :: out
      real(dp), allocatable :: history(:)
      real(dp) :: transition_value, low_mean, transition_mean, trend
      integer :: observations, time, step, lag, position

      if (model%info /= 0 .or. horizon < 1 .or. .not. allocated(model%data)) then
         out%info = 1
         return
      end if
      if (model%external_threshold .and. (.not. present(future_threshold))) then
         out%info = 2
         return
      end if
      if (present(future_threshold)) then
         if (size(future_threshold) < horizon .or. &
            .not. all(ieee_is_finite(future_threshold(:horizon)))) then
            out%info = 2
            return
         end if
      end if
      observations = size(model%data)
      allocate(history(observations + horizon), out%mean(horizon))
      allocate(out%transition_weight(horizon))
      history(:observations) = model%data
      do step = 1, horizon
         time = observations + step
         if (model%external_threshold) then
            transition_value = future_threshold(step)
         else
            transition_value = history(time - model%threshold_lag)
         end if
         out%transition_weight(step) = tsdyn_logistic_transition(transition_value, &
            model%gamma, model%threshold)
         position = 0
         low_mean = 0.0_dp
         transition_mean = 0.0_dp
         if (model%include_constant) then
            position = position + 1
            low_mean = model%low_coefficients(position)
            transition_mean = model%transition_coefficients(position)
         end if
         if (model%include_trend) then
            position = position + 1
            trend = real(time, dp)
            low_mean = low_mean + model%low_coefficients(position)*trend
            transition_mean = transition_mean + &
               model%transition_coefficients(position)*trend
         end if
         do lag = 1, model%low_order
            low_mean = low_mean + model%low_coefficients(position + lag)* &
               history(time - lag)
         end do
         do lag = 1, model%high_order
            transition_mean = transition_mean + &
               model%transition_coefficients(position + lag)*history(time - lag)
         end do
         history(time) = low_mean + out%transition_weight(step)*transition_mean
         out%mean(step) = history(time)
      end do
   end function tsdyn_lstar_forecast

   pure function tsdyn_lstar_select(series, maximum_order, criterion, &
      include_constant, include_trend) result(out)
      !! Select LSTAR low/high orders and transition lag by AIC or BIC.
      real(dp), intent(in) :: series(:) !! Univariate time series.
      integer, intent(in) :: maximum_order !! Largest AR order and transition lag.
      character(len=*), intent(in), optional :: criterion !! Selection criterion, AIC or BIC.
      logical, intent(in), optional :: include_constant !! Whether each component has an intercept.
      logical, intent(in), optional :: include_trend !! Whether each component has a trend.
      type(tsdyn_lstar_selection_t) :: out
      type(tsdyn_lstar_model_t) :: candidate_model
      character(len=:), allocatable :: selected_criterion
      real(dp) :: score, best_score
      integer :: candidates, candidate, low_order, high_order, threshold_lag
      logical :: selected_constant, selected_trend

      selected_criterion = 'aic'
      selected_constant = .true.
      selected_trend = .false.
      if (present(criterion)) selected_criterion = lowercase(trim(criterion))
      if (present(include_constant)) selected_constant = include_constant
      if (present(include_trend)) selected_trend = include_trend
      if (maximum_order < 1 .or. &
         (selected_criterion /= 'aic' .and. selected_criterion /= 'bic')) then
         out%info = 1
         return
      end if
      candidates = maximum_order**3
      allocate(out%low_order(candidates), out%high_order(candidates))
      allocate(out%threshold_lag(candidates), out%aic(candidates))
      allocate(out%bic(candidates), out%gamma(candidates), out%threshold(candidates))
      out%aic = huge(1.0_dp)
      out%bic = huge(1.0_dp)
      out%gamma = 0.0_dp
      out%threshold = 0.0_dp
      best_score = huge(1.0_dp)
      candidate = 0
      do threshold_lag = 1, maximum_order
         do low_order = 1, maximum_order
            do high_order = 1, maximum_order
               candidate = candidate + 1
               out%low_order(candidate) = low_order
               out%high_order(candidate) = high_order
               out%threshold_lag(candidate) = threshold_lag
               candidate_model = tsdyn_lstar_fit(series, low_order, high_order, &
                  threshold_lag, selected_constant, selected_trend)
               if (candidate_model%info /= 0) cycle
               out%aic(candidate) = candidate_model%aic
               out%bic(candidate) = candidate_model%bic
               out%gamma(candidate) = candidate_model%gamma
               out%threshold(candidate) = candidate_model%threshold
               if (selected_criterion == 'aic') then
                  score = candidate_model%aic
               else
                  score = candidate_model%bic
               end if
               if (score < best_score) then
                  best_score = score
                  out%selected = candidate
                  out%model = candidate_model
               end if
            end do
         end do
      end do
      if (out%selected == 0) out%info = 2
   end function tsdyn_lstar_select

   pure function tsdyn_lstar_gradient(model) result(gradient)
      !! Evaluate the LSTAR fitted-value gradient for all model parameters.
      type(tsdyn_lstar_model_t), intent(in) :: model !! Fitted LSTAR model.
      real(dp), allocatable :: gradient(:, :)
      real(dp), allocatable :: response(:), low_design(:, :), high_design(:, :)
      real(dp), allocatable :: internal_threshold(:), transition_mean(:)
      real(dp) :: derivative
      integer :: rows, low_columns, high_columns, row

      if (model%info /= 0 .or. .not. allocated(model%threshold_data)) then
         allocate(gradient(0, 0))
         return
      end if
      call lstar_regression_data(model%data, model%low_order, model%high_order, &
         model%threshold_lag, model%include_constant, model%include_trend, &
         response=response, low_design=low_design, high_design=high_design, &
         transition_data=internal_threshold)
      rows = size(response)
      low_columns = size(low_design, 2)
      high_columns = size(high_design, 2)
      allocate(gradient(rows, low_columns + high_columns + 2))
      gradient(:, :low_columns) = low_design
      do row = 1, rows
         gradient(row, low_columns + 1:low_columns + high_columns) = &
            model%transition_weight(row)*high_design(row, :)
      end do
      transition_mean = matmul(high_design, model%transition_coefficients)
      do row = 1, rows
         derivative = model%transition_weight(row)*(1.0_dp - &
            model%transition_weight(row))
         gradient(row, low_columns + high_columns + 1) = transition_mean(row)* &
            derivative*(model%threshold_data(row) - model%threshold)
         gradient(row, low_columns + high_columns + 2) = -transition_mean(row)* &
            model%gamma*derivative
      end do
   end function tsdyn_lstar_gradient

   pure function tsdyn_lstar_regime_test(model) result(out)
      !! Test an LSTAR fit for an omitted additional smooth-transition regime.
      type(tsdyn_lstar_model_t), intent(in) :: model !! Fitted LSTAR model.
      type(tsdyn_regime_test_t) :: out
      real(dp), allocatable :: gradient(:, :), response(:), base(:, :), unused(:, :)
      real(dp), allocatable :: internal_threshold(:), alternative(:, :)
      integer :: row, columns

      if (model%info /= 0 .or. .not. allocated(model%threshold_data)) then
         out%info = 1
         return
      end if
      gradient = tsdyn_lstar_gradient(model)
      call lstar_regression_data(model%data, max(model%low_order, model%high_order), &
         max(model%low_order, model%high_order), model%threshold_lag, &
         model%include_constant, model%include_trend, response=response, &
         low_design=base, high_design=unused, transition_data=internal_threshold)
      columns = 3*size(base, 2)
      allocate(alternative(size(base, 1), columns))
      do row = 1, size(base, 1)
         alternative(row, 1:size(base, 2)) = base(row, :)*model%threshold_data(row)
         alternative(row, size(base, 2) + 1:2*size(base, 2)) = base(row, :)* &
            model%threshold_data(row)**2
         alternative(row, 2*size(base, 2) + 1:) = base(row, :)* &
            model%threshold_data(row)**3
      end do
      out = star_auxiliary_test(model%residuals, gradient, alternative)
   end function tsdyn_lstar_regime_test

   pure function tsdyn_star_fit(series, order, regimes, threshold_lag, &
      include_constant, include_trend) result(model)
      !! Fit an additive multi-regime logistic STAR model.
      real(dp), intent(in) :: series(:) !! Univariate time series.
      integer, intent(in) :: order !! Common autoregressive order.
      integer, intent(in) :: regimes !! Number of additive STAR regimes.
      integer, intent(in), optional :: threshold_lag !! Positive transition-variable lag.
      logical, intent(in), optional :: include_constant !! Whether each component has an intercept.
      logical, intent(in), optional :: include_trend !! Whether each component has a trend.
      type(tsdyn_star_model_t) :: model
      real(dp), allocatable :: response(:), base(:, :), unused(:, :), transition_data(:)
      real(dp), allocatable :: nonlinear(:), coefficient(:), fitted(:), residual(:)
      real(dp), allocatable :: weight(:, :), candidate_nonlinear(:)
      real(dp) :: rss, candidate_rss, step, scale, ordered_range
      integer :: selected_lag, nonlinear_count, component, iteration, direction, info
      logical :: selected_constant, selected_trend, improved

      selected_lag = 1
      selected_constant = .true.
      selected_trend = .false.
      if (present(threshold_lag)) selected_lag = threshold_lag
      if (present(include_constant)) selected_constant = include_constant
      if (present(include_trend)) selected_trend = include_trend
      if (order < 1 .or. regimes < 2 .or. selected_lag < 1 .or. &
         size(series) <= max(order, selected_lag) + regimes*(order + 2) .or. &
         .not. all(ieee_is_finite(series))) then
         model%info = 1
         return
      end if
      call lstar_regression_data(series, order, order, selected_lag, &
         selected_constant, selected_trend, response=response, low_design=base, &
         high_design=unused, transition_data=transition_data)
      nonlinear_count = 2*(regimes - 1)
      allocate(nonlinear(nonlinear_count), candidate_nonlinear(nonlinear_count))
      scale = sqrt(sum((transition_data - sum(transition_data)/ &
         real(size(transition_data), dp))**2)/real(size(transition_data) - 1, dp))
      scale = max(scale, 1.0e-4_dp)
      do component = 1, regimes - 1
         nonlinear(2*component - 1) = log(5.0_dp/scale)
         nonlinear(2*component) = quantile(sorted(transition_data), &
            real(component, dp)/real(regimes, dp))
      end do
      call star_profile_fit(response, base, transition_data, nonlinear, &
         coefficient, fitted, residual, weight, rss, info)
      if (info /= 0) then
         model%info = 2
         return
      end if
      ordered_range = maxval(transition_data) - minval(transition_data)
      step = max(0.25_dp, 0.10_dp*ordered_range)
      do iteration = 1, 160
         improved = .false.
         do component = 1, nonlinear_count
            do direction = -1, 1, 2
               candidate_nonlinear = nonlinear
               candidate_nonlinear(component) = candidate_nonlinear(component) + &
                  real(direction, dp)*step
               call star_profile_fit(response, base, transition_data, &
                  candidate_nonlinear, coefficient, fitted, residual, weight, &
                  candidate_rss, info)
               if (info == 0 .and. candidate_rss < rss) then
                  nonlinear = candidate_nonlinear
                  rss = candidate_rss
                  improved = .true.
               end if
            end do
         end do
         if (.not. improved) step = 0.5_dp*step
         if (step < 1.0e-6_dp) exit
      end do
      call star_profile_fit(response, base, transition_data, nonlinear, &
         coefficient, fitted, residual, weight, rss, info)
      if (info /= 0) then
         model%info = info
         return
      end if
      model%data = series
      model%coefficients = transpose(reshape(coefficient, &
         [size(base, 2), regimes]))
      allocate(model%gamma(regimes - 1), model%threshold(regimes - 1))
      do component = 1, regimes - 1
         model%gamma(component) = exp(nonlinear(2*component - 1))
         model%threshold(component) = nonlinear(2*component)
      end do
      model%fitted = fitted
      model%residuals = residual
      model%transition_weight = weight
      model%threshold_data = transition_data
      model%rss = rss
      model%order = order
      model%regimes = regimes
      model%threshold_lag = selected_lag
      model%first_fitted = max(order, selected_lag) + 1
      model%include_constant = selected_constant
      model%include_trend = selected_trend
      model%iterations = iteration
      model%converged = step < 1.0e-6_dp
      model%residual_standard_deviation = sqrt(rss/real(size(response) - &
         size(coefficient), dp))
      model%aic = real(size(response), dp)*log(rss/real(size(response), dp)) + &
         2.0_dp*real(size(coefficient) + nonlinear_count, dp)
      model%bic = real(size(response), dp)*log(rss/real(size(response), dp)) + &
         log(real(size(response), dp))*real(size(coefficient) + nonlinear_count, dp)
   end function tsdyn_star_fit

   pure function tsdyn_star_forecast(model, horizon) result(out)
      !! Compute recursive point forecasts from a multi-regime STAR model.
      type(tsdyn_star_model_t), intent(in) :: model !! Fitted multi-regime STAR model.
      integer, intent(in) :: horizon !! Positive forecast horizon.
      type(tsdyn_star_forecast_t) :: out
      real(dp), allocatable :: history(:), regressor(:)
      real(dp) :: transition_value, prediction
      integer :: observations, step_index, time, component, column, lag

      if (model%info /= 0 .or. horizon < 1 .or. .not. allocated(model%data)) then
         out%info = 1
         return
      end if
      observations = size(model%data)
      allocate(history(observations + horizon), out%mean(horizon))
      allocate(out%transition_weight(horizon, model%regimes - 1))
      allocate(regressor(size(model%coefficients, 2)))
      history(:observations) = model%data
      do step_index = 1, horizon
         time = observations + step_index
         column = 0
         if (model%include_constant) then
            column = column + 1
            regressor(column) = 1.0_dp
         end if
         if (model%include_trend) then
            column = column + 1
            regressor(column) = real(time, dp)
         end if
         do lag = 1, model%order
            regressor(column + lag) = history(time - lag)
         end do
         transition_value = history(time - model%threshold_lag)
         prediction = dot_product(model%coefficients(1, :), regressor)
         do component = 1, model%regimes - 1
            out%transition_weight(step_index, component) = &
               tsdyn_logistic_transition(transition_value, model%gamma(component), &
               model%threshold(component))
            prediction = prediction + out%transition_weight(step_index, component)* &
               dot_product(model%coefficients(component + 1, :), regressor)
         end do
         history(time) = prediction
         out%mean(step_index) = prediction
      end do
   end function tsdyn_star_forecast

   pure function tsdyn_star_gradient(model) result(gradient)
      !! Evaluate the fitted-value gradient of a multi-regime STAR model.
      type(tsdyn_star_model_t), intent(in) :: model !! Fitted multi-regime STAR model.
      real(dp), allocatable :: gradient(:, :)
      real(dp), allocatable :: response(:), base(:, :), unused(:, :), threshold_data(:)
      real(dp) :: derivative, component_mean
      integer :: rows, columns, component, row, offset

      if (model%info /= 0 .or. .not. allocated(model%transition_weight)) then
         allocate(gradient(0, 0))
         return
      end if
      call lstar_regression_data(model%data, model%order, model%order, &
         model%threshold_lag, model%include_constant, model%include_trend, &
         response=response, low_design=base, high_design=unused, &
         transition_data=threshold_data)
      rows = size(base, 1)
      columns = size(base, 2)
      allocate(gradient(rows, model%regimes*columns + 2*(model%regimes - 1)))
      gradient(:, :columns) = base
      do component = 1, model%regimes - 1
         offset = component*columns
         do row = 1, rows
            gradient(row, offset + 1:offset + columns) = &
               model%transition_weight(row, component)*base(row, :)
            derivative = model%transition_weight(row, component)*(1.0_dp - &
               model%transition_weight(row, component))
            component_mean = dot_product(model%coefficients(component + 1, :), &
               base(row, :))
            gradient(row, model%regimes*columns + 2*component - 1) = &
               component_mean*derivative*(model%threshold_data(row) - &
               model%threshold(component))
            gradient(row, model%regimes*columns + 2*component) = &
               -component_mean*derivative*model%gamma(component)
         end do
      end do
   end function tsdyn_star_gradient

   pure function tsdyn_star_regime_test(model) result(out)
      !! Test a multi-regime STAR fit for another smooth-transition component.
      type(tsdyn_star_model_t), intent(in) :: model !! Fitted multi-regime STAR model.
      type(tsdyn_regime_test_t) :: out
      real(dp), allocatable :: gradient(:, :), response(:), base(:, :), unused(:, :)
      real(dp), allocatable :: threshold_data(:), alternative(:, :)
      integer :: row, columns

      if (model%info /= 0) then
         out%info = 1
         return
      end if
      gradient = tsdyn_star_gradient(model)
      call lstar_regression_data(model%data, model%order, model%order, &
         model%threshold_lag, model%include_constant, model%include_trend, &
         response=response, low_design=base, high_design=unused, &
         transition_data=threshold_data)
      columns = size(base, 2)
      allocate(alternative(size(base, 1), 3*columns))
      do row = 1, size(base, 1)
         alternative(row, :columns) = base(row, :)*model%threshold_data(row)
         alternative(row, columns + 1:2*columns) = base(row, :)* &
            model%threshold_data(row)**2
         alternative(row, 2*columns + 1:) = base(row, :)*model%threshold_data(row)**3
      end do
      out = star_auxiliary_test(model%residuals, gradient, alternative)
   end function tsdyn_star_regime_test

   pure function tsdyn_star_add_regime(model) result(expanded)
      !! Refit a STAR model after adding one logistic transition component.
      type(tsdyn_star_model_t), intent(in) :: model !! Existing multi-regime STAR model.
      type(tsdyn_star_model_t) :: expanded

      if (model%info /= 0) then
         expanded%info = 1
         return
      end if
      expanded = tsdyn_star_fit(model%data, model%order, model%regimes + 1, &
         model%threshold_lag, model%include_constant, model%include_trend)
   end function tsdyn_star_add_regime

   pure function tsdyn_llar_fit(series, order, delay, steps, epsilon_values) &
      result(model)
      !! Select and fit a local-linear autoregression over epsilon neighborhoods.
      real(dp), intent(in) :: series(:) !! Univariate time series.
      integer, intent(in) :: order !! Delay-embedding dimension.
      integer, intent(in), optional :: delay !! Positive embedding delay.
      integer, intent(in), optional :: steps !! Positive forecast step and Theiler window.
      real(dp), intent(in), optional :: epsilon_values(:) !! Positive neighborhood radii.
      type(tsdyn_llar_model_t) :: model
      real(dp), allocatable :: candidate_epsilon(:), target(:)
      real(dp) :: minimum_epsilon, maximum_epsilon, prediction, error_sum
      real(dp) :: target_sum, target_square_sum, target_sd, nan
      integer :: selected_delay, selected_steps, embedding_count
      integer :: epsilon_index, query, neighbors, usable, total_neighbors, info

      selected_delay = 1
      selected_steps = selected_delay
      if (present(delay)) selected_delay = delay
      if (present(steps)) selected_steps = steps
      embedding_count = size(series) - (order - 1)*selected_delay - selected_steps
      if (order < 1 .or. selected_delay < 1 .or. selected_steps < 1 .or. &
         embedding_count < 2*(order + 1) + 2 .or. &
         .not. all(ieee_is_finite(series)) .or. &
         maxval(series) <= minval(series)) then
         model%info = 1
         return
      end if
      if (present(epsilon_values)) then
         if (size(epsilon_values) < 1 .or. any(epsilon_values <= 0.0_dp) .or. &
            .not. all(ieee_is_finite(epsilon_values))) then
            model%info = 1
            return
         end if
         candidate_epsilon = epsilon_values
      else
         minimum_epsilon = 0.5_dp*sqrt(sum((series - sum(series)/ &
            real(size(series), dp))**2)/real(size(series) - 1, dp))
         maximum_epsilon = maxval(series) - minval(series)
         allocate(candidate_epsilon(30))
         do epsilon_index = 1, 30
            candidate_epsilon(epsilon_index) = exp(log(minimum_epsilon) + &
               real(epsilon_index - 1, dp)*(log(maximum_epsilon) - &
               log(minimum_epsilon))/29.0_dp)
         end do
      end if
      model%data = series
      model%epsilon = candidate_epsilon
      allocate(model%normalized_rmse(size(candidate_epsilon)))
      allocate(model%usable_fraction(size(candidate_epsilon)))
      allocate(model%average_neighbors(size(candidate_epsilon)))
      nan = ieee_value(0.0_dp, ieee_quiet_nan)
      model%normalized_rmse = nan
      model%usable_fraction = 0.0_dp
      model%average_neighbors = nan
      allocate(target(embedding_count))
      do query = 1, embedding_count
         target(query) = series(query + (order - 1)*selected_delay + selected_steps)
      end do
      do epsilon_index = 1, size(candidate_epsilon)
         error_sum = 0.0_dp
         target_sum = 0.0_dp
         target_square_sum = 0.0_dp
         usable = 0
         total_neighbors = 0
         do query = 1, embedding_count
            call llar_embedded_prediction(series, order, selected_delay, &
               selected_steps, embedding_count, query, candidate_epsilon(epsilon_index), &
               prediction, neighbors, info)
            if (info /= 0) cycle
            error_sum = error_sum + (prediction - target(query))**2
            target_sum = target_sum + target(query)
            target_square_sum = target_square_sum + target(query)**2
            usable = usable + 1
            total_neighbors = total_neighbors + neighbors
         end do
         model%usable_fraction(epsilon_index) = real(usable, dp)/ &
            real(embedding_count, dp)
         if (usable > 1) then
            target_sd = sqrt(max(0.0_dp, (target_square_sum - &
               target_sum*target_sum/real(usable, dp))/real(usable - 1, dp)))
            if (target_sd > tiny(1.0_dp)) then
               model%normalized_rmse(epsilon_index) = &
                  sqrt(error_sum/real(usable, dp))/target_sd
            end if
            model%average_neighbors(epsilon_index) = &
               real(total_neighbors, dp)/real(usable, dp)
         end if
      end do
      model%selected = llar_minimum_finite_index(model%normalized_rmse)
      if (model%selected == 0) then
         model%info = 2
         return
      end if
      model%selected_epsilon = model%epsilon(model%selected)
      allocate(model%fitted(embedding_count), model%neighbor_count(embedding_count))
      model%fitted = nan
      model%neighbor_count = 0
      do query = 1, embedding_count
         call llar_embedded_prediction(series, order, selected_delay, &
            selected_steps, embedding_count, query, model%selected_epsilon, &
            prediction, neighbors, info)
         model%neighbor_count(query) = neighbors
         if (info == 0) model%fitted(query) = prediction
      end do
      model%order = order
      model%delay = selected_delay
      model%steps = selected_steps
      model%first_fitted = (order - 1)*selected_delay + selected_steps + 1
   end function tsdyn_llar_fit

   pure function tsdyn_llar_forecast(model, horizon, enlarge, &
      expansion_percent) result(out)
      !! Recursively forecast a fitted local-linear autoregression.
      type(tsdyn_llar_model_t), intent(in) :: model !! Fitted local-linear AR model.
      integer, intent(in) :: horizon !! Positive number of recursive forecasts.
      logical, intent(in), optional :: enlarge !! Whether to enlarge empty neighborhoods.
      real(dp), intent(in), optional :: expansion_percent !! Radius increase percentage.
      type(tsdyn_llar_forecast_t) :: out
      real(dp), allocatable :: history(:)
      real(dp) :: epsilon, expansion, prediction
      integer :: observations, step_index, neighbors, info
      logical :: allow_enlargement

      allow_enlargement = .true.
      expansion = 20.0_dp
      if (present(enlarge)) allow_enlargement = enlarge
      if (present(expansion_percent)) expansion = expansion_percent
      if (model%info /= 0 .or. horizon < 1 .or. expansion <= 0.0_dp .or. &
         .not. allocated(model%data)) then
         out%info = 1
         return
      end if
      observations = size(model%data)
      allocate(history(observations + horizon), out%mean(horizon))
      allocate(out%epsilon(horizon), out%neighbor_count(horizon))
      history(:observations) = model%data
      do step_index = 1, horizon
         epsilon = model%selected_epsilon
         do
            call llar_next_prediction(history(:observations + step_index - 1), &
               model%order, model%delay, model%steps, epsilon, prediction, &
               neighbors, info)
            if (info == 0 .or. .not. allow_enlargement) exit
            epsilon = epsilon*(1.0_dp + expansion/100.0_dp)
            if (epsilon > 1000.0_dp*(maxval(history(:observations + &
               step_index - 1)) - minval(history(:observations + &
               step_index - 1)))) exit
         end do
         if (info /= 0) then
            out%info = 2
            return
         end if
         history(observations + step_index) = prediction
         out%mean(step_index) = prediction
         out%epsilon(step_index) = epsilon
         out%neighbor_count(step_index) = neighbors
      end do
   end function tsdyn_llar_forecast

   pure function tsdyn_aar_fit(series, order, delay, steps, basis_count, &
      lambda_values) result(model)
      !! Fit an additive autoregression with GCV-selected lag-spline penalties.
      real(dp), intent(in) :: series(:) !! Univariate time series.
      integer, intent(in) :: order !! Number of additive lag components.
      integer, intent(in), optional :: delay !! Positive delay between lag components.
      integer, intent(in), optional :: steps !! Positive forecast step.
      integer, intent(in), optional :: basis_count !! Cubic B-spline functions per lag.
      real(dp), intent(in), optional :: lambda_values(:) !! Candidate nonnegative penalties.
      type(tsdyn_aar_model_t) :: model
      type(penalized_regression_t) :: regression, candidate_regression
      real(dp), allocatable :: response(:), lag_data(:, :), design(:, :)
      real(dp), allocatable :: base_penalty(:, :), full_penalty(:, :)
      real(dp), allocatable :: candidates(:), lambda(:)
      real(dp), allocatable :: basis_matrix(:, :)
      real(dp) :: best_gcv, best_lambda
      integer :: selected_delay, selected_steps, selected_basis, free_count, rows
      integer :: component, candidate, iteration, column_start, index

      selected_delay = 1
      selected_steps = selected_delay
      selected_basis = 8
      if (present(delay)) selected_delay = delay
      if (present(steps)) selected_steps = steps
      if (present(basis_count)) selected_basis = basis_count
      free_count = selected_basis - 1
      rows = size(series) - (order - 1)*selected_delay - selected_steps
      if (order < 1 .or. selected_delay < 1 .or. selected_steps < 1 .or. &
         selected_basis < 4 .or. rows <= 1 + order*free_count .or. &
         .not. all(ieee_is_finite(series))) then
         model%info = 1
         return
      end if
      if (present(lambda_values)) then
         if (size(lambda_values) < 1 .or. any(lambda_values < 0.0_dp) .or. &
            .not. all(ieee_is_finite(lambda_values))) then
            model%info = 1
            return
         end if
         candidates = lambda_values
      else
         allocate(candidates(17))
         do candidate = 1, size(candidates)
            candidates(candidate) = exp(log(1.0e-4_dp) + real(candidate - 1, dp)* &
               (log(1.0e4_dp) - log(1.0e-4_dp))/real(size(candidates) - 1, dp))
         end do
      end if
      call aar_embedding_data(series, order, selected_delay, selected_steps, &
         response, lag_data)
      allocate(model%basis(order), model%basis_center(selected_basis, order))
      allocate(design(rows, 1 + order*free_count))
      design(:, 1) = 1.0_dp
      do component = 1, order
         model%basis(component) = spline_basis_create(lag_data(:, component), &
            selected_basis)
         if (model%basis(component)%info /= 0) then
            model%info = 2
            return
         end if
         basis_matrix = spline_basis_matrix(model%basis(component), &
            lag_data(:, component))
         model%basis_center(:, component) = sum(basis_matrix, dim=1)/real(rows, dp)
         column_start = 2 + (component - 1)*free_count
         do index = 1, rows
            design(index, column_start:column_start + free_count - 1) = &
               basis_matrix(index, :free_count) - &
               model%basis_center(:free_count, component)
         end do
      end do
      full_penalty = spline_difference_penalty(selected_basis, 2)
      base_penalty = full_penalty(:free_count, :free_count)
      allocate(lambda(order))
      lambda = 1.0_dp
      regression = aar_penalized_fit(design, response, base_penalty, lambda, &
         free_count)
      if (regression%info /= 0) then
         model%info = 3
         return
      end if
      do iteration = 1, 5
         do component = 1, order
            best_gcv = regression%gcv
            best_lambda = lambda(component)
            do candidate = 1, size(candidates)
               lambda(component) = candidates(candidate)
               candidate_regression = aar_penalized_fit(design, response, &
                  base_penalty, lambda, free_count)
               if (candidate_regression%info == 0 .and. &
                  candidate_regression%gcv < best_gcv) then
                  best_gcv = candidate_regression%gcv
                  best_lambda = candidates(candidate)
                  regression = candidate_regression
               end if
            end do
            lambda(component) = best_lambda
         end do
      end do
      model%data = series
      model%lambda = lambda
      model%intercept = regression%coefficients(1)
      allocate(model%coefficients(selected_basis, order))
      model%coefficients = 0.0_dp
      do component = 1, order
         column_start = 2 + (component - 1)*free_count
         model%coefficients(:free_count, component) = regression%coefficients( &
            column_start:column_start + free_count - 1)
      end do
      model%fitted = regression%fitted
      model%residuals = regression%residuals
      model%rss = regression%rss
      model%gcv = regression%gcv
      model%effective_df = regression%effective_df
      call aar_component_degrees(design, base_penalty, lambda, free_count, &
         model%component_df)
      model%order = order
      model%delay = selected_delay
      model%steps = selected_steps
      model%basis_count = selected_basis
      model%first_fitted = (order - 1)*selected_delay + selected_steps + 1
   end function tsdyn_aar_fit

   pure function tsdyn_aar_component(model, component, values) result(effect)
      !! Evaluate one centered additive lag component at supplied values.
      type(tsdyn_aar_model_t), intent(in) :: model !! Fitted additive AR model.
      integer, intent(in) :: component !! One-based lag-component index.
      real(dp), intent(in) :: values(:) !! Values at which to evaluate the component.
      real(dp), allocatable :: effect(:)
      real(dp), allocatable :: matrix(:, :)
      integer :: row

      if (model%info /= 0 .or. component < 1 .or. component > model%order) then
         allocate(effect(0))
         return
      end if
      matrix = spline_basis_matrix(model%basis(component), values)
      do row = 1, size(matrix, 1)
         matrix(row, :) = matrix(row, :) - model%basis_center(:, component)
      end do
      effect = matmul(matrix, model%coefficients(:, component))
   end function tsdyn_aar_component

   pure function tsdyn_aar_forecast(model, horizon) result(out)
      !! Recursively forecast an additive autoregression.
      type(tsdyn_aar_model_t), intent(in) :: model !! Fitted additive AR model.
      integer, intent(in) :: horizon !! Positive forecast horizon.
      type(tsdyn_aar_forecast_t) :: out
      real(dp), allocatable :: history(:), basis_value(:)
      real(dp) :: prediction
      integer :: observations, step_index, time, component, position

      if (model%info /= 0 .or. horizon < 1 .or. .not. allocated(model%data)) then
         out%info = 1
         return
      end if
      observations = size(model%data)
      allocate(history(observations + horizon), out%mean(horizon))
      history(:observations) = model%data
      do step_index = 1, horizon
         time = observations + step_index
         prediction = model%intercept
         do component = 1, model%order
            position = time - 1 - (component - 1)*model%delay
            basis_value = spline_basis_values(model%basis(component), history(position))
            prediction = prediction + dot_product(basis_value - &
               model%basis_center(:, component), model%coefficients(:, component))
         end do
         history(time) = prediction
         out%mean(step_index) = prediction
      end do
   end function tsdyn_aar_forecast

   pure function tsdyn_nnet_fit(series, order, hidden_count, delay, steps, &
      max_iterations, tolerance, decay, initial_parameters) result(model)
      !! Fit a delay-embedded neural-network autoregression with linear output.
      real(dp), intent(in) :: series(:) !! Observed univariate series.
      integer, intent(in) :: order !! Delay-embedding dimension.
      integer, intent(in) :: hidden_count !! Number of sigmoid hidden units.
      integer, intent(in), optional :: delay !! Delay between embedding components.
      integer, intent(in), optional :: steps !! Direct forecast step.
      integer, intent(in), optional :: max_iterations !! Maximum analytic BFGS iterations.
      real(dp), intent(in), optional :: tolerance !! Gradient convergence tolerance.
      real(dp), intent(in), optional :: decay !! L2 network weight decay.
      real(dp), intent(in), optional :: initial_parameters(:) !! Initial packed network weights.
      type(tsdyn_nnet_model_t) :: model
      real(dp), allocatable :: response(:), predictors(:, :), response_matrix(:, :)
      real(dp), allocatable :: prediction(:, :)
      real(dp) :: selected_tolerance, selected_decay
      integer :: selected_delay, selected_steps, selected_iterations, rows, parameters

      selected_delay = 1
      if (present(delay)) selected_delay = delay
      selected_steps = selected_delay
      if (present(steps)) selected_steps = steps
      selected_iterations = 1000
      if (present(max_iterations)) selected_iterations = max_iterations
      selected_tolerance = 1.0e-6_dp
      if (present(tolerance)) selected_tolerance = tolerance
      selected_decay = 0.0_dp
      if (present(decay)) selected_decay = decay
      rows = size(series) - (order - 1)*selected_delay - selected_steps
      if (order < 1 .or. hidden_count < 1 .or. selected_delay < 1 .or. &
         selected_steps < 1 .or. rows < 3 .or. selected_iterations < 1 .or. &
         selected_tolerance <= 0.0_dp .or. selected_decay < 0.0_dp .or. &
         .not. all(ieee_is_finite(series))) then
         model%info = 1
         return
      end if
      call aar_embedding_data(series, order, selected_delay, selected_steps, &
         response, predictors)
      allocate(response_matrix(rows, 1))
      response_matrix(:, 1) = response
      model%network = neural_network_fit(predictors, response_matrix, hidden_count, &
         selected_iterations, selected_tolerance, selected_decay, initial_parameters)
      if (model%network%info /= 0) then
         model%info = 2
         return
      end if
      prediction = neural_network_predict(model%network, predictors)
      model%data = series
      model%fitted = prediction(:, 1)
      model%residuals = response - model%fitted
      model%rss = sum(model%residuals**2)
      parameters = neural_network_parameter_count(order, hidden_count, 1)
      model%aic = real(rows, dp)*log(max(model%rss/real(rows, dp), &
         tiny(1.0_dp))) + 2.0_dp*real(parameters, dp)
      model%bic = real(rows, dp)*log(max(model%rss/real(rows, dp), &
         tiny(1.0_dp))) + log(real(rows, dp))*real(parameters, dp)
      model%order = order
      model%delay = selected_delay
      model%steps = selected_steps
      model%hidden_count = hidden_count
      model%first_fitted = (order - 1)*selected_delay + selected_steps + 1
   end function tsdyn_nnet_fit

   pure function tsdyn_nnet_predict(model, predictors) result(prediction)
      !! Evaluate a fitted neural autoregression for supplied embedding rows.
      type(tsdyn_nnet_model_t), intent(in) :: model !! Fitted neural autoregression.
      real(dp), intent(in) :: predictors(:, :) !! Embedding rows ordered newest to oldest.
      real(dp), allocatable :: prediction(:)
      real(dp), allocatable :: matrix(:, :)

      if (model%info /= 0 .or. size(predictors, 2) /= model%order) then
         allocate(prediction(0))
         return
      end if
      matrix = neural_network_predict(model%network, predictors)
      prediction = matrix(:, 1)
   end function tsdyn_nnet_predict

   pure function tsdyn_nnet_forecast(model, horizon) result(out)
      !! Recursively forecast a one-step neural-network autoregression.
      type(tsdyn_nnet_model_t), intent(in) :: model !! Fitted neural autoregression.
      integer, intent(in) :: horizon !! Positive forecast horizon.
      type(tsdyn_nnet_forecast_t) :: out
      real(dp), allocatable :: history(:), predictors(:, :), prediction(:)
      integer :: origin, step, component

      if (model%info /= 0 .or. horizon < 1 .or. model%steps /= 1) then
         out%info = 1
         return
      end if
      origin = size(model%data)
      allocate(history(origin + horizon), out%mean(horizon))
      allocate(predictors(1, model%order))
      history(:origin) = model%data
      do step = 1, horizon
         do component = 1, model%order
            predictors(1, component) = history(origin + step - 1 - &
               (component - 1)*model%delay)
         end do
         prediction = tsdyn_nnet_predict(model, predictors)
         if (size(prediction) /= 1 .or. .not. ieee_is_finite(prediction(1))) then
            out%info = 2
            return
         end if
         out%mean(step) = prediction(1)
         history(origin + step) = prediction(1)
      end do
   end function tsdyn_nnet_forecast

   pure function tsdyn_nnet_select(series, order, hidden_counts, criterion, &
      delay, steps, max_iterations, tolerance, decay) result(out)
      !! Select neural-network hidden size using AIC or BIC.
      real(dp), intent(in) :: series(:) !! Observed univariate series.
      integer, intent(in) :: order !! Delay-embedding dimension.
      integer, intent(in) :: hidden_counts(:) !! Candidate positive hidden-unit counts.
      character(len=*), intent(in), optional :: criterion !! `aic` or `bic`.
      integer, intent(in), optional :: delay !! Delay between embedding components.
      integer, intent(in), optional :: steps !! Direct forecast step.
      integer, intent(in), optional :: max_iterations !! Maximum BFGS iterations per fit.
      real(dp), intent(in), optional :: tolerance !! Gradient convergence tolerance.
      real(dp), intent(in), optional :: decay !! L2 network weight decay.
      type(tsdyn_nnet_selection_t) :: out
      type(tsdyn_nnet_model_t) :: candidate
      character(len=:), allocatable :: selected_criterion
      real(dp) :: score, best_score
      integer :: candidate_index

      selected_criterion = 'aic'
      if (present(criterion)) selected_criterion = lowercase(criterion)
      if (size(hidden_counts) < 1 .or. any(hidden_counts < 1) .or. &
         (selected_criterion /= 'aic' .and. selected_criterion /= 'bic')) then
         out%info = 1
         return
      end if
      out%hidden_count = hidden_counts
      allocate(out%aic(size(hidden_counts)), out%bic(size(hidden_counts)))
      out%aic = huge(1.0_dp)
      out%bic = huge(1.0_dp)
      best_score = huge(1.0_dp)
      do candidate_index = 1, size(hidden_counts)
         candidate = tsdyn_nnet_fit(series, order, hidden_counts(candidate_index), &
            delay, steps, max_iterations, tolerance, decay)
         if (candidate%info /= 0) cycle
         out%aic(candidate_index) = candidate%aic
         out%bic(candidate_index) = candidate%bic
         if (selected_criterion == 'aic') then
            score = candidate%aic
         else
            score = candidate%bic
         end if
         if (score < best_score) then
            best_score = score
            out%selected = candidate_index
            out%model = candidate
         end if
      end do
      if (out%selected == 0) out%info = 2
   end function tsdyn_nnet_select

   pure function tsdyn_tvecm_fit(data, lag, threshold_count, trim, beta_values, &
      threshold_grid_size, only_error_correction, include_constant, &
      include_trend, cointegrating_coefficients, fixed_thresholds, &
      threshold_bounds, cointegration_offset, shared_deterministic) result(model)
      !! Estimate a threshold VECM by cointegration and threshold grid search.
      real(dp), intent(in) :: data(:, :) !! Time-by-variable level observations.
      integer, intent(in) :: lag !! Number of lagged differences.
      integer, intent(in), optional :: threshold_count !! One or two thresholds.
      real(dp), intent(in), optional :: trim !! Minimum regime fraction.
      real(dp), intent(in), optional :: beta_values(:) !! Cointegrating coefficients to search.
      integer, intent(in), optional :: threshold_grid_size !! Threshold candidates per beta.
      logical, intent(in), optional :: only_error_correction !! Switch only ECT coefficients.
      logical, intent(in), optional :: include_constant !! Include deterministic intercepts.
      logical, intent(in), optional :: include_trend !! Include deterministic trends.
      real(dp), intent(in), optional :: cointegrating_coefficients(:) !! Fixed coefficients after normalized y1.
      real(dp), intent(in), optional :: fixed_thresholds(:) !! Exact ordered thresholds, bypassing grid search.
      real(dp), intent(in), optional :: threshold_bounds(:, :) !! Threshold-by-lower/upper search bounds.
      real(dp), intent(in), optional :: cointegration_offset(:) !! Supplied long-run deterministic or exogenous contribution.
      logical, intent(in), optional :: shared_deterministic !! Share intercept and trend across regimes.
      type(tsdyn_tvecm_model_t) :: model
      real(dp), allocatable :: betas(:), candidates(:, :), response(:, :)
      real(dp), allocatable :: common(:, :), ect(:)
      real(dp), allocatable :: ordered(:), threshold_grid(:), thresholds(:)
      real(dp), allocatable :: design(:, :), coefficient(:, :), fitted(:, :)
      real(dp), allocatable :: residual(:, :)
      integer, allocatable :: regime(:)
      real(dp) :: selected_trim, beta_ols, beta_se, denominator, rss, best_rss
      integer :: selected_thresholds, selected_grid, beta_index, first, second
      integer :: rows, variables, deterministic_count, shared_columns
      integer :: info, regime_index
      logical :: selected_only_ect, selected_constant, selected_trend
      logical :: selected_shared_deterministic

      selected_thresholds = 1
      selected_trim = 0.05_dp
      selected_grid = 25
      selected_only_ect = .false.
      selected_constant = .true.
      selected_trend = .false.
      selected_shared_deterministic = .false.
      if (present(threshold_count)) selected_thresholds = threshold_count
      if (present(trim)) selected_trim = trim
      if (present(threshold_grid_size)) selected_grid = threshold_grid_size
      if (present(only_error_correction)) selected_only_ect = only_error_correction
      if (present(include_constant)) selected_constant = include_constant
      if (present(include_trend)) selected_trend = include_trend
      if (present(shared_deterministic)) selected_shared_deterministic = &
         shared_deterministic
      variables = size(data, 2)
      deterministic_count = merge(1, 0, selected_constant) + &
         merge(1, 0, selected_trend)
      shared_columns = merge(deterministic_count, 0, &
         selected_shared_deterministic .and. .not. selected_only_ect)
      rows = size(data, 1) - lag - 1
      if (variables < 2 .or. lag < 1 .or. &
         (selected_thresholds /= 1 .and. selected_thresholds /= 2) .or. &
         selected_trim <= 0.0_dp .or. selected_trim >= 1.0_dp/ &
         real(selected_thresholds + 1, dp) .or. selected_grid < 3 .or. &
         rows < 20 .or. .not. all(ieee_is_finite(data))) then
         model%info = 1
         return
      end if
      if (present(cointegration_offset)) then
         if (size(cointegration_offset) /= size(data, 1) .or. &
            .not. all(ieee_is_finite(cointegration_offset))) then
            model%info = 1
            return
         end if
      end if
      if (present(fixed_thresholds)) then
         if (size(fixed_thresholds) /= selected_thresholds .or. &
            .not. all(ieee_is_finite(fixed_thresholds))) then
            model%info = 1
            return
         end if
         if (selected_thresholds == 2) then
            if (fixed_thresholds(2) <= fixed_thresholds(1)) then
               model%info = 1
               return
            end if
         end if
      end if
      if (present(threshold_bounds)) then
         if (any(shape(threshold_bounds) /= [selected_thresholds, 2])) then
            model%info = 1
            return
         end if
         if (.not. all(ieee_is_finite(threshold_bounds)) .or. &
            any(threshold_bounds(:, 2) < threshold_bounds(:, 1))) then
            model%info = 1
            return
         end if
         if (present(fixed_thresholds)) then
            if (any(fixed_thresholds < threshold_bounds(:, 1)) .or. &
               any(fixed_thresholds > threshold_bounds(:, 2))) then
               model%info = 1
               return
            end if
         end if
      end if
      if (present(cointegrating_coefficients)) then
         if (size(cointegrating_coefficients) /= variables - 1 .or. &
            .not. all(ieee_is_finite(cointegrating_coefficients)) .or. &
            present(beta_values)) then
            model%info = 1
            return
         end if
         allocate(candidates(variables, 1))
         candidates(:, 1) = [1.0_dp, -cointegrating_coefficients]
      else if (variables > 2) then
         model%info = 2
         return
      else if (present(beta_values)) then
         if (size(beta_values) < 1 .or. .not. all(ieee_is_finite(beta_values))) then
            model%info = 1
            return
         end if
         betas = beta_values
      else
         denominator = sum(data(:, 2)*data(:, 2))
         if (denominator <= tiny(1.0_dp)) then
            model%info = 2
            return
         end if
         beta_ols = sum(data(:, 1)*data(:, 2))/denominator
         beta_se = sqrt(sum((data(:, 1) - beta_ols*data(:, 2))**2)/ &
            real(size(data, 1) - 1, dp)/denominator)
         beta_se = max(beta_se, 0.02_dp*max(1.0_dp, abs(beta_ols)))
         allocate(betas(21))
         do beta_index = 1, size(betas)
            betas(beta_index) = beta_ols - 2.0_dp*beta_se + &
               4.0_dp*beta_se*real(beta_index - 1, dp)/real(size(betas) - 1, dp)
         end do
      end if
      if (.not. allocated(candidates)) then
         allocate(candidates(2, size(betas)))
         candidates(1, :) = 1.0_dp
         candidates(2, :) = -betas
      end if
      best_rss = huge(1.0_dp)
      do beta_index = 1, size(candidates, 2)
         if (present(cointegration_offset)) then
            call tvecm_regression_data(data, lag, candidates(:, beta_index), &
               selected_constant, selected_trend, response, common, ect, &
               cointegration_offset)
         else
            call tvecm_regression_data(data, lag, candidates(:, beta_index), &
               selected_constant, selected_trend, response, common, ect)
         end if
         if (present(fixed_thresholds)) then
            thresholds = fixed_thresholds
            call tvecm_design(ect, common, thresholds, selected_only_ect, &
               selected_trim, design, regime, info, shared_columns)
            if (info == 0) then
               call tvecm_multivariate_regression(design, response, coefficient, &
                  fitted, residual, rss, info)
            end if
            if (info == 0 .and. rss < best_rss) then
               best_rss = rss
               model%cointegration = candidates(:, beta_index)
               if (variables == 2) model%beta = -candidates(2, beta_index)
               model%thresholds = thresholds
            end if
            cycle
         end if
         ordered = sorted(ect)
         allocate(threshold_grid(selected_grid))
         do first = 1, selected_grid
            threshold_grid(first) = quantile(ordered, selected_trim + &
               (1.0_dp - 2.0_dp*selected_trim)*real(first - 1, dp)/ &
               real(selected_grid - 1, dp))
         end do
         if (selected_thresholds == 1) then
            do first = 1, selected_grid
               thresholds = [threshold_grid(first)]
               if (present(threshold_bounds)) then
                  if (thresholds(1) < threshold_bounds(1, 1) .or. &
                     thresholds(1) > threshold_bounds(1, 2)) cycle
               end if
               call tvecm_design(ect, common, thresholds, selected_only_ect, &
                  selected_trim, design, regime, info, shared_columns)
               if (info /= 0) cycle
               call tvecm_multivariate_regression(design, response, coefficient, &
                  fitted, residual, rss, info)
               if (info == 0 .and. rss < best_rss) then
                  best_rss = rss
                  model%cointegration = candidates(:, beta_index)
                  if (variables == 2) model%beta = -candidates(2, beta_index)
                  model%thresholds = thresholds
               end if
            end do
         else
            do first = 1, selected_grid - 1
               do second = first + 1, selected_grid
                  thresholds = [threshold_grid(first), threshold_grid(second)]
                  if (present(threshold_bounds)) then
                     if (thresholds(1) < threshold_bounds(1, 1) .or. &
                        thresholds(1) > threshold_bounds(1, 2) .or. &
                        thresholds(2) < threshold_bounds(2, 1) .or. &
                        thresholds(2) > threshold_bounds(2, 2)) cycle
                  end if
                  call tvecm_design(ect, common, thresholds, selected_only_ect, &
                     selected_trim, design, regime, info, shared_columns)
                  if (info /= 0) cycle
                  call tvecm_multivariate_regression(design, response, coefficient, &
                     fitted, residual, rss, info)
                  if (info == 0 .and. rss < best_rss) then
                     best_rss = rss
                     model%cointegration = candidates(:, beta_index)
                     if (variables == 2) model%beta = -candidates(2, beta_index)
                     model%thresholds = thresholds
                  end if
               end do
            end do
         end if
         deallocate(threshold_grid)
      end do
      if (.not. allocated(model%thresholds)) then
         model%info = 3
         return
      end if
      if (present(cointegration_offset)) then
         call tvecm_regression_data(data, lag, model%cointegration, &
            selected_constant, selected_trend, response, common, ect, &
            cointegration_offset)
      else
         call tvecm_regression_data(data, lag, model%cointegration, &
            selected_constant, selected_trend, response, common, ect)
      end if
      call tvecm_design(ect, common, model%thresholds, selected_only_ect, &
         selected_trim, design, regime, info, shared_columns)
      call tvecm_multivariate_regression(design, response, coefficient, fitted, &
         residual, rss, info)
      if (info /= 0) then
         model%info = info
         return
      end if
      model%data = data
      model%coefficients = coefficient
      model%error_correction = ect
      model%fitted = fitted
      model%residuals = residual
      model%regime = regime
      model%rss = rss
      model%covariance = matmul(transpose(residual), residual)/real(rows, dp)
      model%lag = lag
      model%threshold_count = selected_thresholds
      model%trim = selected_trim
      model%first_fitted = lag + 2
      model%only_error_correction = selected_only_ect
      model%shared_deterministic = selected_shared_deterministic
      model%include_constant = selected_constant
      model%include_trend = selected_trend
      if (present(cointegration_offset)) then
         model%cointegration_offset = cointegration_offset
      else
         allocate(model%cointegration_offset(size(data, 1)))
         model%cointegration_offset = 0.0_dp
      end if
      model%aic = real(rows, dp)*log(rss/real(variables*rows, dp)) + &
         2.0_dp*real(size(coefficient) + selected_thresholds + variables - 1, dp)
      model%bic = real(rows, dp)*log(rss/real(variables*rows, dp)) + &
         log(real(rows, dp))*real(size(coefficient) + selected_thresholds + &
         variables - 1, dp)
      do regime_index = 1, selected_thresholds + 1
         if (count(regime == regime_index) == 0) then
            model%info = 4
            return
         end if
      end do
   end function tsdyn_tvecm_fit

   pure function tsdyn_tvecm_forecast(model, horizon) result(out)
      !! Recursively forecast levels from a fitted threshold VECM.
      type(tsdyn_tvecm_model_t), intent(in) :: model !! Fitted threshold VECM.
      integer, intent(in) :: horizon !! Positive forecast horizon.
      type(tsdyn_tvecm_forecast_t) :: out
      real(dp), allocatable :: history(:, :), change(:)
      integer :: observations, variables, step_index, time, regime, info

      if (model%info /= 0 .or. horizon < 1 .or. .not. allocated(model%data)) then
         out%info = 1
         return
      end if
      observations = size(model%data, 1)
      variables = size(model%data, 2)
      allocate(history(observations + horizon, variables), &
         out%mean(horizon, variables))
      allocate(out%regime(horizon), change(variables))
      history(:observations, :) = model%data
      do step_index = 1, horizon
         time = observations + step_index
         call tvecm_next_change(model, history, time, &
            size(model%residuals, 1) + step_index, change, regime, info)
         if (info /= 0) then
            out%info = info
            return
         end if
         history(time, :) = history(time - 1, :) + change
         out%mean(step_index, :) = history(time, :)
         out%regime(step_index) = regime
      end do
   end function tsdyn_tvecm_forecast

   pure function tsdyn_tvecm_simulate_from_innovations(model, starting, &
      innovations) result(out)
      !! Simulate a fitted TVECM from starting levels and supplied innovations.
      type(tsdyn_tvecm_model_t), intent(in) :: model !! TVECM coefficient specification.
      real(dp), intent(in) :: starting(:, :) !! At least lag-plus-one starting level rows.
      real(dp), intent(in) :: innovations(:, :) !! Time-by-variable supplied innovations.
      type(tsdyn_tvecm_simulation_t) :: out
      real(dp), allocatable :: history(:, :), change(:)
      integer :: initial_rows, variables, step_index, time, regime, info
      integer :: offset_time, trend_index

      initial_rows = size(starting, 1)
      variables = size(starting, 2)
      if (model%info /= 0 .or. variables < 2 .or. &
         initial_rows < model%lag + 1 .or. &
         size(innovations, 2) /= variables .or. &
         size(innovations, 1) < 1 .or. .not. all(ieee_is_finite(starting)) .or. &
         .not. all(ieee_is_finite(innovations))) then
         out%info = 1
         return
      end if
      allocate(history(initial_rows + size(innovations, 1), variables), &
         change(variables))
      allocate(out%regime(size(innovations, 1)))
      history(:initial_rows, :) = starting
      do step_index = 1, size(innovations, 1)
         time = initial_rows + step_index
         offset_time = time - 1
         if (allocated(model%data)) offset_time = &
            size(model%data, 1) + step_index - 1
         trend_index = step_index
         if (allocated(model%residuals)) trend_index = &
            size(model%residuals, 1) + step_index
         call tvecm_next_change(model, history, time, trend_index, change, &
            regime, info, offset_time)
         if (info /= 0) then
            out%info = info
            return
         end if
         history(time, :) = history(time - 1, :) + change + innovations(step_index, :)
         out%regime(step_index) = regime
      end do
      out%series = history(initial_rows + 1:, :)
      out%innovations = innovations
   end function tsdyn_tvecm_simulate_from_innovations

   pure function tsdyn_setar_bootstrap(model, innovations, refit, trim, &
      grid_count) result(out)
      !! Generate and optionally refit a SETAR residual-bootstrap series.
      type(nts_tar_model_t), intent(in) :: model !! Fitted level-threshold SETAR model.
      real(dp), intent(in) :: innovations(:) !! Resampled fitted residuals in time order.
      logical, intent(in), optional :: refit !! Re-estimate thresholds and coefficients when true.
      real(dp), intent(in), optional :: trim !! Minimum regime fraction for refitting.
      integer, intent(in), optional :: grid_count !! Maximum threshold grid size for refitting.
      type(tsdyn_setar_bootstrap_t) :: out
      real(dp) :: selected_trim
      integer :: selected_grid, start, row, time, regime, lag
      logical :: perform_refit

      perform_refit = .false.
      if (present(refit)) perform_refit = refit
      selected_trim = 0.1_dp
      if (present(trim)) selected_trim = trim
      selected_grid = 40
      if (present(grid_count)) selected_grid = grid_count
      if (model%info /= 0 .or. .not. allocated(model%data) .or. &
         .not. allocated(model%residuals) .or. size(innovations) /= &
         size(model%residuals) .or. .not. all(ieee_is_finite(innovations))) then
         out%info = 1
         return
      end if
      start = size(model%data) - size(model%residuals) + 1
      if (start <= max(maxval(model%ar_order), model%delay)) then
         out%info = 2
         return
      end if
      allocate(out%series(size(model%data)))
      out%series(:start - 1) = model%data(:start - 1)
      do row = 1, size(innovations)
         time = start + row - 1
         regime = tvecm_regime(out%series(time - model%delay), &
            model%thresholds)
         out%series(time) = model%coefficients(regime, 1) + innovations(row)
         do lag = 1, model%ar_order(regime)
            out%series(time) = out%series(time) + &
               model%coefficients(regime, lag + 1)*out%series(time - lag)
         end do
      end do
      out%innovations = innovations
      if (perform_refit) then
         out%fitted_model = tsdyn_setar_fit(out%series, model%ar_order, &
            size(model%thresholds), model%delay, selected_trim, selected_grid, &
            include_mean=model%include_mean)
         if (out%fitted_model%info /= 0) then
            out%info = 10 + out%fitted_model%info
            return
         end if
         out%refitted = .true.
      end if
   end function tsdyn_setar_bootstrap

   pure function tsdyn_tvar_bootstrap(model, innovations, refit, trim, &
      grid_count) result(out)
      !! Generate and optionally refit a TVAR residual-bootstrap series.
      type(nts_mtar_model_t), intent(in) :: model !! Fitted level-threshold TVAR model.
      real(dp), intent(in) :: innovations(:, :) !! Resampled residual rows in time order.
      logical, intent(in), optional :: refit !! Re-estimate thresholds and coefficients when true.
      real(dp), intent(in), optional :: trim !! Minimum regime fraction for refitting.
      integer, intent(in), optional :: grid_count !! Maximum threshold grid size for refitting.
      type(tsdyn_tvar_bootstrap_t) :: out
      real(dp) :: selected_trim
      integer :: selected_grid, start, row, time, regime, lag
      logical :: perform_refit

      perform_refit = .false.
      if (present(refit)) perform_refit = refit
      selected_trim = 0.1_dp
      if (present(trim)) selected_trim = trim
      selected_grid = 30
      if (present(grid_count)) selected_grid = grid_count
      if (model%info /= 0 .or. .not. allocated(model%data) .or. &
         .not. allocated(model%residuals) .or. size(innovations, 1) /= &
         size(model%residuals, 1) .or. size(innovations, 2) /= &
         size(model%data, 2) .or. .not. all(ieee_is_finite(innovations))) then
         out%info = 1
         return
      end if
      start = size(model%data, 1) - size(model%residuals, 1) + 1
      if (start <= max(maxval(model%ar_order), model%delay)) then
         out%info = 2
         return
      end if
      allocate(out%series(size(model%data, 1), size(model%data, 2)))
      out%series(:start - 1, :) = model%data(:start - 1, :)
      do row = 1, size(innovations, 1)
         time = start + row - 1
         regime = tvecm_regime(out%series(time - model%delay, &
            model%threshold_component), model%thresholds)
         out%series(time, :) = model%intercept(:, regime) + innovations(row, :)
         do lag = 1, model%ar_order(regime)
            out%series(time, :) = out%series(time, :) + &
               matmul(model%ar(:, :, lag, regime), out%series(time - lag, :))
         end do
      end do
      out%innovations = innovations
      if (perform_refit) then
         out%fitted_model = tsdyn_tvar_fit(out%series, model%ar_order, &
            size(model%thresholds), model%threshold_component, model%delay, &
            selected_trim, selected_grid, include_mean=model%include_mean)
         if (out%fitted_model%info /= 0) then
            out%info = 10 + out%fitted_model%info
            return
         end if
         out%refitted = .true.
      end if
   end function tsdyn_tvar_bootstrap

   pure function tsdyn_tvecm_bootstrap(model, innovations, refit, beta_values, &
      threshold_grid_size, cointegrating_coefficients, fixed_thresholds, &
      threshold_bounds) result(out)
      !! Generate and optionally refit a TVECM residual-bootstrap series.
      type(tsdyn_tvecm_model_t), intent(in) :: model !! Fitted threshold VECM.
      real(dp), intent(in) :: innovations(:, :) !! Resampled multivariate residual rows.
      logical, intent(in), optional :: refit !! Re-estimate thresholds and coefficients when true.
      real(dp), intent(in), optional :: beta_values(:) !! Cointegrating coefficients searched during refitting.
      integer, intent(in), optional :: threshold_grid_size !! Threshold candidates per beta during refitting.
      real(dp), intent(in), optional :: cointegrating_coefficients(:) !! Fixed multivariate refit coefficients after y1.
      real(dp), intent(in), optional :: fixed_thresholds(:) !! Exact thresholds retained during refitting.
      real(dp), intent(in), optional :: threshold_bounds(:, :) !! Threshold search bounds during refitting.
      type(tsdyn_tvecm_bootstrap_t) :: out
      real(dp), allocatable :: change(:), refit_cointegration(:)
      integer :: selected_grid, start, variables, row, time, regime, status
      logical :: perform_refit

      perform_refit = .false.
      if (present(refit)) perform_refit = refit
      selected_grid = 25
      if (present(threshold_grid_size)) selected_grid = threshold_grid_size
      if (model%info /= 0 .or. .not. allocated(model%data) .or. &
         .not. allocated(model%residuals)) then
         out%info = 1
         return
      end if
      variables = size(model%data, 2)
      if (size(innovations, 1) /= &
         size(model%residuals, 1) .or. size(innovations, 2) /= variables .or. &
         .not. all(ieee_is_finite(innovations)) .or. &
         (variables > 2 .and. present(beta_values))) then
         out%info = 1
         return
      end if
      start = size(model%data, 1) - size(model%residuals, 1) + 1
      if (start <= model%lag + 1) then
         out%info = 2
         return
      end if
      allocate(out%series(size(model%data, 1), variables), change(variables))
      out%series(:start - 1, :) = model%data(:start - 1, :)
      do row = 1, size(innovations, 1)
         time = start + row - 1
         call tvecm_next_change(model, out%series, time, row, change, regime, &
            status)
         if (status /= 0) then
            out%info = 3
            return
         end if
         out%series(time, :) = out%series(time - 1, :) + change + &
            innovations(row, :)
      end do
      out%innovations = innovations
      if (perform_refit) then
         if (variables > 2) then
            if (present(cointegrating_coefficients)) then
               refit_cointegration = cointegrating_coefficients
            else if (allocated(model%cointegration)) then
               refit_cointegration = -model%cointegration(2:)
            else
               out%info = 4
               return
            end if
            out%fitted_model = tsdyn_tvecm_fit(out%series, model%lag, &
               model%threshold_count, model%trim, &
               threshold_grid_size=selected_grid, &
               only_error_correction=model%only_error_correction, &
               include_constant=model%include_constant, &
               include_trend=model%include_trend, &
               cointegrating_coefficients=refit_cointegration, &
               fixed_thresholds=fixed_thresholds, &
               threshold_bounds=threshold_bounds, &
               cointegration_offset=model%cointegration_offset, &
               shared_deterministic=model%shared_deterministic)
         else if (present(cointegrating_coefficients)) then
            out%fitted_model = tsdyn_tvecm_fit(out%series, model%lag, &
               model%threshold_count, model%trim, &
               threshold_grid_size=selected_grid, &
               only_error_correction=model%only_error_correction, &
               include_constant=model%include_constant, &
               include_trend=model%include_trend, &
               cointegrating_coefficients=cointegrating_coefficients, &
               fixed_thresholds=fixed_thresholds, &
               threshold_bounds=threshold_bounds, &
               cointegration_offset=model%cointegration_offset, &
               shared_deterministic=model%shared_deterministic)
         else if (present(beta_values)) then
            out%fitted_model = tsdyn_tvecm_fit(out%series, model%lag, &
               model%threshold_count, model%trim, beta_values, selected_grid, &
               model%only_error_correction, model%include_constant, &
               model%include_trend, fixed_thresholds=fixed_thresholds, &
               threshold_bounds=threshold_bounds, &
               cointegration_offset=model%cointegration_offset, &
               shared_deterministic=model%shared_deterministic)
         else
            out%fitted_model = tsdyn_tvecm_fit(out%series, model%lag, &
               model%threshold_count, model%trim, &
               threshold_grid_size=selected_grid, &
               only_error_correction=model%only_error_correction, &
               include_constant=model%include_constant, &
               include_trend=model%include_trend, &
               fixed_thresholds=fixed_thresholds, &
               threshold_bounds=threshold_bounds, &
               cointegration_offset=model%cointegration_offset, &
               shared_deterministic=model%shared_deterministic)
         end if
         if (out%fitted_model%info /= 0) then
            out%info = 10 + out%fitted_model%info
            return
         end if
         out%refitted = .true.
      end if
   end function tsdyn_tvecm_bootstrap

   pure function tsdyn_setar_fit(series, ar_order, threshold_count, delay, &
      trim, grid_count, momentum, threshold_variable, include_mean, &
      include_trend, representation, lag_spacing, forecast_step, &
      lag_indices) result(model)
      !! Fit a one- or two-threshold SETAR model by trimmed grid search.
      real(dp), intent(in) :: series(:) !! Observed univariate series.
      integer, intent(in) :: ar_order(:) !! Regime-specific autoregressive orders.
      integer, intent(in) :: threshold_count !! Number of thresholds, one or two.
      integer, intent(in) :: delay !! Delay applied to the transition variable.
      real(dp), intent(in), optional :: trim !! Minimum fraction in each outer regime.
      integer, intent(in), optional :: grid_count !! Maximum threshold grid size.
      logical, intent(in), optional :: momentum !! Use delayed first differences for switching.
      real(dp), intent(in), optional :: threshold_variable(:) !! External transition series.
      logical, intent(in), optional :: include_mean(:) !! Regime intercept flags.
      logical, intent(in), optional :: include_trend(:) !! Regime trend flags.
      character(len=*), intent(in), optional :: representation !! `level`, `diff`, or `ADF`.
      integer, intent(in), optional :: lag_spacing !! Spacing between candidate lag positions.
      integer, intent(in), optional :: forecast_step !! Lead from the first predictor to response.
      integer, intent(in), optional :: lag_indices(:, :) !! Selected candidate-lag indices by regime.
      type(nts_tar_model_t) :: model
      type(nts_tar_model_t) :: candidate_model
      real(dp), allocatable :: transition(:), ordered(:), grid(:)
      real(dp) :: selected_trim, score, best_score
      integer :: selected_grid, start, usable, i, j, maximum_lag
      logical :: use_momentum
      character(len=:), allocatable :: selected_representation
      integer :: selected_spacing, selected_step

      selected_trim = 0.1_dp
      if (present(trim)) selected_trim = trim
      selected_grid = 40
      if (present(grid_count)) selected_grid = grid_count
      use_momentum = .false.
      if (present(momentum)) use_momentum = momentum
      selected_representation = 'level'
      if (present(representation)) selected_representation = &
         lowercase(adjustl(representation))
      selected_spacing = 1
      if (present(lag_spacing)) selected_spacing = lag_spacing
      selected_step = 1
      if (present(forecast_step)) selected_step = forecast_step
      if (threshold_count < 1 .or. threshold_count > 2 .or. &
         size(ar_order) /= threshold_count + 1 .or. selected_trim <= 0.0_dp .or. &
         selected_trim >= 1.0_dp/real(threshold_count + 1, dp) .or. &
         selected_grid < 2 .or. delay < 1 .or. selected_spacing < 1 .or. &
         selected_step < 1 .or. (selected_representation /= 'level' .and. &
         selected_representation /= 'diff' .and. &
         selected_representation /= 'adf')) then
         model%info = 1
         return
      end if
      maximum_lag = setar_maximum_lag(ar_order, selected_spacing, selected_step, &
         lag_indices)
      if (maximum_lag < 0) then
         model%info = 1
         return
      end if
      call setar_transition_series(series, delay, use_momentum, &
         threshold_variable, transition, start)
      if (.not. allocated(transition)) then
         model%info = 2
         return
      end if
      usable = size(series) - max(start, maximum_lag + 2) + 1
      if (selected_representation == 'level') usable = &
         size(series) - max(start, maximum_lag + 1) + 1
      if (usable < 10) then
         model%info = 3
         return
      end if
      ordered = sorted(transition(size(transition) - usable + 1:))
      call threshold_candidate_grid(ordered, selected_trim, selected_grid, grid)
      if (size(grid) < threshold_count) then
         model%info = 4
         return
      end if
      best_score = huge(1.0_dp)
      if (threshold_count == 1) then
         do i = 1, size(grid)
            candidate_model = setar_design_estimate(series, ar_order, grid(i:i), &
               delay, transition, start, include_mean, include_trend, &
               selected_representation, selected_spacing, selected_step, lag_indices)
            if (candidate_model%info /= 0) cycle
            score = sum(candidate_model%residuals**2)
            if (score < best_score) then
               best_score = score
               model = candidate_model
            end if
         end do
      else
         do i = 1, size(grid) - 1
            do j = i + 1, size(grid)
            candidate_model = setar_design_estimate(series, ar_order, &
               [grid(i), grid(j)], delay, transition, start, include_mean, include_trend, &
               selected_representation, selected_spacing, selected_step, lag_indices)
               if (candidate_model%info /= 0) cycle
               if (minval(real(candidate_model%regime_observations, dp)) < &
                  selected_trim*real(size(candidate_model%residuals), dp)) cycle
               score = sum(candidate_model%residuals**2)
               if (score < best_score) then
                  best_score = score
                  model = candidate_model
               end if
            end do
         end do
      end if
      if (best_score == huge(1.0_dp)) model%info = 5
   end function tsdyn_setar_fit

   pure integer function setar_maximum_lag(ar_order, lag_spacing, forecast_step, &
      lag_indices) result(maximum_lag)
      !! Return the largest actual lag implied by a generalized SETAR design.
      integer, intent(in) :: ar_order(:) !! Numbers of selected lags by regime.
      integer, intent(in) :: lag_spacing !! Spacing between candidate lag positions.
      integer, intent(in) :: forecast_step !! Lead from first predictor to response.
      integer, intent(in), optional :: lag_indices(:, :) !! Candidate-lag indices by regime.
      integer :: regime, term, candidate

      maximum_lag = -1
      if (size(ar_order) < 1 .or. any(ar_order < 0) .or. lag_spacing < 1 .or. &
         forecast_step < 1) return
      maximum_lag = 0
      do regime = 1, size(ar_order)
         do term = 1, ar_order(regime)
            candidate = term
            if (present(lag_indices)) then
               if (size(lag_indices, 1) < size(ar_order) .or. &
                  size(lag_indices, 2) < ar_order(regime)) then
                  maximum_lag = -1
                  return
               end if
               candidate = lag_indices(regime, term)
               if (candidate < 1) then
                  maximum_lag = -1
                  return
               end if
            end if
            maximum_lag = max(maximum_lag, forecast_step + &
               (candidate - 1)*lag_spacing)
         end do
      end do
   end function setar_maximum_lag

   pure function setar_design_estimate(series, ar_order, thresholds, delay, &
      transition, first_transition, include_mean, include_trend, representation, &
      lag_spacing, forecast_step, lag_indices) result(model)
      !! Fit a fixed-threshold SETAR with generalized lags and representations.
      real(dp), intent(in) :: series(:) !! Observed univariate series.
      integer, intent(in) :: ar_order(:) !! Numbers of selected lags by regime.
      real(dp), intent(in) :: thresholds(:) !! Ordered fixed thresholds.
      integer, intent(in) :: delay !! Delay applied to the transition variable.
      real(dp), intent(in) :: transition(:) !! Complete transition series.
      integer, intent(in) :: first_transition !! First response with valid transition.
      logical, intent(in), optional :: include_mean(:) !! Regime intercept flags.
      logical, intent(in), optional :: include_trend(:) !! Regime trend flags.
      character(len=*), intent(in) :: representation !! `level`, `diff`, or `adf`.
      integer, intent(in) :: lag_spacing !! Spacing between candidate lag positions.
      integer, intent(in) :: forecast_step !! Lead from first predictor to response.
      integer, intent(in), optional :: lag_indices(:, :) !! Candidate-lag indices by regime.
      type(nts_tar_model_t) :: model
      real(dp), allocatable :: design(:, :), response(:), coefficient(:)
      real(dp), allocatable :: standard_error(:), residual(:), regime_rss(:)
      logical, allocatable :: means(:), trends(:), selected_rows(:)
      integer, allocatable :: actual_lag(:, :), row_regime(:)
      real(dp) :: rss
      integer :: regimes, maximum_lag, start, rows, regime, row, time, term
      integer :: columns, column, status, candidate, count_rows
      logical :: differenced, adf

      regimes = size(ar_order)
      maximum_lag = setar_maximum_lag(ar_order, lag_spacing, forecast_step, &
         lag_indices)
      differenced = representation == 'diff' .or. representation == 'adf'
      adf = representation == 'adf'
      if (regimes /= size(thresholds) + 1 .or. maximum_lag < 0 .or. &
         size(transition) /= size(series)) then
         model%info = 1
         return
      end if
      allocate(means(regimes), trends(regimes))
      means = .true.
      trends = .false.
      if (present(include_mean)) then
         if (size(include_mean) /= regimes) then
            model%info = 1
            return
         end if
         means = include_mean
      end if
      if (present(include_trend)) then
         if (size(include_trend) /= regimes) then
            model%info = 1
            return
         end if
         trends = include_trend
      end if
      allocate(actual_lag(regimes, max(1, maxval(ar_order))))
      actual_lag = 0
      do regime = 1, regimes
         do term = 1, ar_order(regime)
            candidate = term
            if (present(lag_indices)) candidate = lag_indices(regime, term)
            actual_lag(regime, term) = forecast_step + &
               (candidate - 1)*lag_spacing
         end do
      end do
      start = max(first_transition, maximum_lag + 1)
      if (differenced) start = max(start, maximum_lag + 2)
      rows = size(series) - start + 1
      if (rows < 10) then
         model%info = 2
         return
      end if
      allocate(model%coefficients(regimes, maximum_lag + 1))
      allocate(model%trend_coefficient(regimes), model%adf_coefficient(regimes))
      allocate(model%lag_active(regimes, max(1, maximum_lag)))
      allocate(model%ar_order(regimes), model%lag_count(regimes))
      allocate(model%regime_observations(regimes), model%innovation_sd(regimes))
      allocate(model%include_mean(regimes), model%include_trend(regimes))
      allocate(model%residuals(rows), model%standardized_residuals(rows))
      allocate(row_regime(rows), regime_rss(regimes))
      model%coefficients = 0.0_dp
      model%trend_coefficient = 0.0_dp
      model%adf_coefficient = 0.0_dp
      model%lag_active = .false.
      model%ar_order = 0
      model%lag_count = ar_order
      model%include_mean = means
      model%include_trend = trends
      model%representation = representation
      model%forecast_step = forecast_step
      model%residuals = 0.0_dp
      do row = 1, rows
         time = start + row - 1
         row_regime(row) = tvecm_regime(transition(time - delay), thresholds)
      end do
      model%aic = 0.0_dp
      do regime = 1, regimes
         selected_rows = row_regime == regime
         count_rows = count(selected_rows)
         columns = merge(1, 0, means(regime)) + merge(1, 0, trends(regime)) + &
            merge(1, 0, adf) + ar_order(regime)
         if (columns < 1 .or. count_rows <= columns) then
            model%info = 3
            return
         end if
         allocate(design(count_rows, columns), response(count_rows))
         row = 0
         do time = start, size(series)
            if (row_regime(time - start + 1) /= regime) cycle
            row = row + 1
            column = 0
            if (means(regime)) then
               column = column + 1
               design(row, column) = 1.0_dp
            end if
            if (trends(regime)) then
               column = column + 1
               design(row, column) = real(time - start + 1, dp)
            end if
            if (adf) then
               column = column + 1
               design(row, column) = series(time - forecast_step)
            end if
            do term = 1, ar_order(regime)
               column = column + 1
               if (differenced) then
                  design(row, column) = series(time - actual_lag(regime, term)) - &
                     series(time - actual_lag(regime, term) - 1)
               else
                  design(row, column) = series(time - actual_lag(regime, term))
               end if
            end do
            if (differenced) then
               response(row) = series(time) - series(time - 1)
            else
               response(row) = series(time)
            end if
         end do
         call ols_fit(design, response, coefficient, standard_error, residual, &
            rss, status)
         if (status /= 0 .or. rss <= 0.0_dp) then
            model%info = 10 + status
            return
         end if
         column = 0
         if (means(regime)) then
            column = column + 1
            model%coefficients(regime, 1) = coefficient(column)
         end if
         if (trends(regime)) then
            column = column + 1
            model%trend_coefficient(regime) = coefficient(column)
         end if
         if (adf) then
            column = column + 1
            model%adf_coefficient(regime) = coefficient(column)
         end if
         do term = 1, ar_order(regime)
            column = column + 1
            model%coefficients(regime, actual_lag(regime, term) + 1) = &
               coefficient(column)
            model%lag_active(regime, actual_lag(regime, term)) = .true.
            model%ar_order(regime) = max(model%ar_order(regime), &
               actual_lag(regime, term))
         end do
         model%residuals = unpack(residual, selected_rows, model%residuals)
         regime_rss(regime) = rss
         model%regime_observations(regime) = count_rows
         model%innovation_sd(regime) = sqrt(rss/real(count_rows - columns, dp))
         model%aic = model%aic + real(count_rows, dp)* &
            log(rss/real(count_rows, dp)) + 2.0_dp*real(columns, dp)
         deallocate(design, response, coefficient, standard_error, residual)
      end do
      do row = 1, rows
         model%standardized_residuals(row) = model%residuals(row)/ &
            model%innovation_sd(row_regime(row))
      end do
      model%data = series
      model%thresholds = thresholds
      model%delay = delay
   end function setar_design_estimate

   pure function tsdyn_setar_restricted_fit(series, ar_order, thresholds, &
      delay, common, outer_symmetric, momentum, threshold_variable, &
      include_mean) result(model)
      !! Fit a fixed-threshold SETAR with shared or outer-symmetric parameters.
      real(dp), intent(in) :: series(:) !! Observed univariate series.
      integer, intent(in) :: ar_order(:) !! Regime-specific autoregressive orders.
      real(dp), intent(in) :: thresholds(:) !! Ordered fixed thresholds.
      integer, intent(in) :: delay !! Delay applied to the transition variable.
      character(len=*), intent(in), optional :: common !! `none`, `include`, `lags`, or `both`.
      logical, intent(in), optional :: outer_symmetric !! Share low and high regime parameters.
      logical, intent(in), optional :: momentum !! Use delayed first differences for switching.
      real(dp), intent(in), optional :: threshold_variable(:) !! External transition series.
      logical, intent(in), optional :: include_mean(:) !! Regime intercept flags.
      type(nts_tar_model_t) :: model
      real(dp), allocatable :: transition(:), design(:, :), response(:), beta(:)
      real(dp), allocatable :: standard_error(:), residual(:)
      real(dp), allocatable :: regime_rss(:)
      integer, allocatable :: map(:, :), row_regime(:)
      logical, allocatable :: means(:)
      character(len=:), allocatable :: common_mode
      real(dp) :: rss
      integer :: regimes, maximum_order, start, first_transition, rows
      integer :: row, time, regime, lag_index, column, status, count_regime_rows
      logical :: share_include, share_lags, symmetric, use_momentum

      regimes = size(thresholds) + 1
      maximum_order = maxval(ar_order)
      common_mode = 'none'
      if (present(common)) common_mode = lowercase(trim(common))
      share_include = common_mode == 'include' .or. common_mode == 'both'
      share_lags = common_mode == 'lags' .or. common_mode == 'both'
      symmetric = .false.
      if (present(outer_symmetric)) symmetric = outer_symmetric
      use_momentum = .false.
      if (present(momentum)) use_momentum = momentum
      if (regimes < 2 .or. size(ar_order) /= regimes .or. any(ar_order < 0) .or. &
         delay < 1 .or. (common_mode /= 'none' .and. &
         common_mode /= 'include' .and. common_mode /= 'lags' .and. &
         common_mode /= 'both') .or. (symmetric .and. regimes /= 3)) then
         model%info = 1
         return
      end if
      if (size(thresholds) > 1) then
         if (any(thresholds(2:) <= thresholds(:size(thresholds) - 1))) then
            model%info = 1
            return
         end if
      end if
      allocate(means(regimes))
      means = .true.
      if (present(include_mean)) then
         if (size(include_mean) /= regimes) then
            model%info = 2
            return
         end if
         means = include_mean
      end if
      if (share_include .and. any(means .neqv. means(1))) then
         model%info = 2
         return
      end if
      if (symmetric) then
         if (means(1) .neqv. means(regimes)) then
            model%info = 2
            return
         end if
      end if
      call setar_transition_series(series, delay, use_momentum, &
         threshold_variable, transition, first_transition)
      if (.not. allocated(transition)) then
         model%info = 3
         return
      end if
      start = max(maximum_order + 1, first_transition)
      rows = size(series) - start + 1
      allocate(map(regimes, maximum_order + 1))
      map = 0
      column = 0
      if (share_include .and. means(1)) then
         column = column + 1
         map(:, 1) = column
      else
         do regime = 1, regimes
            if (.not. means(regime)) cycle
            if (symmetric .and. regime == regimes) then
               map(regime, 1) = map(1, 1)
            else
               column = column + 1
               map(regime, 1) = column
            end if
         end do
      end if
      do lag_index = 1, maximum_order
         if (share_lags) then
            if (any(ar_order >= lag_index)) then
               column = column + 1
               where (ar_order >= lag_index) map(:, lag_index + 1) = column
            end if
         else
            do regime = 1, regimes
               if (ar_order(regime) < lag_index) cycle
               if (symmetric .and. regime == regimes .and. &
                  ar_order(1) >= lag_index) then
                  map(regime, lag_index + 1) = map(1, lag_index + 1)
               else
                  column = column + 1
                  map(regime, lag_index + 1) = column
               end if
            end do
         end if
      end do
      if (rows <= column) then
         model%info = 4
         return
      end if
      allocate(design(rows, column), response(rows), row_regime(rows))
      design = 0.0_dp
      do row = 1, rows
         time = start + row - 1
         row_regime(row) = tvecm_regime(transition(time - delay), thresholds)
         regime = row_regime(row)
         response(row) = series(time)
         if (map(regime, 1) > 0) design(row, map(regime, 1)) = 1.0_dp
         do lag_index = 1, ar_order(regime)
            design(row, map(regime, lag_index + 1)) = series(time - lag_index)
         end do
      end do
      call ols_fit(design, response, beta, standard_error, residual, rss, status)
      if (status /= 0 .or. rss <= 0.0_dp) then
         model%info = 10 + status
         return
      end if
      allocate(model%coefficients(regimes, maximum_order + 1))
      allocate(model%innovation_sd(regimes), model%regime_observations(regimes))
      allocate(model%residuals(rows), model%standardized_residuals(rows))
      allocate(regime_rss(regimes))
      model%coefficients = 0.0_dp
      do regime = 1, regimes
         do lag_index = 0, maximum_order
            if (map(regime, lag_index + 1) > 0) then
               model%coefficients(regime, lag_index + 1) = &
                  beta(map(regime, lag_index + 1))
            end if
         end do
         count_regime_rows = count(row_regime == regime)
         model%regime_observations(regime) = count_regime_rows
         regime_rss(regime) = sum(merge(residual**2, 0.0_dp, row_regime == regime))
         if (count_regime_rows < 2 .or. regime_rss(regime) <= 0.0_dp) then
            model%info = 5
            return
         end if
         model%innovation_sd(regime) = sqrt(regime_rss(regime)/ &
            real(count_regime_rows, dp))
      end do
      model%residuals = residual
      do row = 1, rows
         model%standardized_residuals(row) = residual(row)/ &
            model%innovation_sd(row_regime(row))
      end do
      model%aic = 2.0_dp*real(column, dp)
      do regime = 1, regimes
         model%aic = model%aic + real(model%regime_observations(regime), dp)* &
            log(regime_rss(regime)/real(model%regime_observations(regime), dp))
      end do
      model%data = series
      model%ar_order = ar_order
      model%thresholds = thresholds
      model%include_mean = means
      model%delay = delay
   end function tsdyn_setar_restricted_fit

   pure function tsdyn_setar_simulate_from_innovations(model, innovations, &
      burnin) result(out)
      !! Simulate a fitted SETAR specification from supplied innovations.
      type(nts_tar_model_t), intent(in) :: model !! SETAR parameter specification.
      real(dp), intent(in) :: innovations(:) !! Innovation sequence including burn-in.
      integer, intent(in), optional :: burnin !! Initial observations to discard.
      type(nts_tar_simulation_t) :: out
      real(dp), allocatable :: ar(:, :), intercept(:)
      integer :: regimes, regime, discard

      discard = 0
      if (present(burnin)) discard = burnin
      if (model%info /= 0 .or. .not. allocated(model%coefficients) .or. &
         .not. allocated(model%ar_order) .or. .not. allocated(model%thresholds)) then
         out%info = 1
         return
      end if
      if (allocated(model%lag_active)) then
         out = setar_general_simulate_from_innovations(model, innovations, discard)
         return
      end if
      regimes = size(model%ar_order)
      allocate(ar(regimes, size(model%coefficients, 2) - 1), intercept(regimes))
      ar = model%coefficients(:, 2:)
      intercept = model%coefficients(:, 1)
      out = nts_utar_simulate_from_innovations(intercept, ar, model%ar_order, &
         model%thresholds, model%innovation_sd, model%delay, innovations, discard)
   end function tsdyn_setar_simulate_from_innovations

   pure function setar_general_simulate_from_innovations(model, innovations, &
      burnin) result(out)
      !! Simulate a generalized SETAR model from supplied standard innovations.
      type(nts_tar_model_t), intent(in) :: model !! Generalized SETAR parameter specification.
      real(dp), intent(in) :: innovations(:) !! Standard-normal innovations including burn-in.
      integer, intent(in) :: burnin !! Number of initial observations to discard.
      type(nts_tar_simulation_t) :: out
      real(dp), allocatable :: work(:), errors(:)
      real(dp) :: prediction
      integer :: total, maximum_lag, start, time, regime, lag

      total = size(innovations)
      maximum_lag = max(1, maxval(model%ar_order))
      start = max(maximum_lag + 1, model%delay + 1)
      if (model%representation /= 'level') start = max(start, maximum_lag + 2)
      if (burnin < 0 .or. burnin >= total .or. total < start .or. &
         .not. allocated(model%innovation_sd) .or. &
         .not. allocated(model%include_trend) .or. &
         .not. all(ieee_is_finite(innovations))) then
         out%info = 1
         return
      end if
      allocate(work(total), errors(total))
      work = 0.0_dp
      errors = 0.0_dp
      work(:start - 1) = model%innovation_sd(1)*innovations(:start - 1)
      errors(:start - 1) = work(:start - 1)
      do time = start, total
         regime = tvecm_regime(work(time - model%delay), model%thresholds)
         prediction = model%coefficients(regime, 1)
         if (model%include_trend(regime)) prediction = prediction + &
            model%trend_coefficient(regime)*real(time - start + 1, dp)
         if (model%representation == 'adf') prediction = prediction + &
            model%adf_coefficient(regime)*work(time - model%forecast_step)
         do lag = 1, model%ar_order(regime)
            if (.not. model%lag_active(regime, lag)) cycle
            if (model%representation == 'level') then
               prediction = prediction + model%coefficients(regime, lag + 1)* &
                  work(time - lag)
            else
               prediction = prediction + model%coefficients(regime, lag + 1)* &
                  (work(time - lag) - work(time - lag - 1))
            end if
         end do
         errors(time) = model%innovation_sd(regime)*innovations(time)
         prediction = prediction + errors(time)
         if (model%representation == 'level') then
            work(time) = prediction
         else
            work(time) = work(time - 1) + prediction
         end if
      end do
      out%series = work(burnin + 1:)
      out%innovations = errors(burnin + 1:)
      out%intercept = model%coefficients(:, 1)
      out%trend_coefficient = model%trend_coefficient
      out%adf_coefficient = model%adf_coefficient
      out%ar = model%coefficients(:, 2:)
      out%ar_order = model%ar_order
      out%lag_active = model%lag_active
      out%include_trend = model%include_trend
      out%thresholds = model%thresholds
      out%innovation_sd = model%innovation_sd
      out%representation = model%representation
      out%forecast_step = model%forecast_step
      out%delay = model%delay
      out%burnin = burnin
   end function setar_general_simulate_from_innovations

   pure function tsdyn_setar_select(series, ar_order, threshold_count, &
      maximum_delay, trim, grid_count, momentum, threshold_variable, &
      include_mean, criterion, order_candidates, threshold_counts, &
      same_lags, common, outer_symmetric, same_sample, &
      threshold_variables, include_trend, representation, lag_spacing, &
      forecast_step, lag_indices) result(out)
      !! Select a SETAR specification over thresholds, delays, and lag orders.
      real(dp), intent(in) :: series(:) !! Observed univariate series.
      integer, intent(in) :: ar_order(:) !! Regime-specific autoregressive orders.
      integer, intent(in) :: threshold_count !! Number of thresholds.
      integer, intent(in) :: maximum_delay !! Largest transition delay considered.
      real(dp), intent(in), optional :: trim !! Minimum outer-regime fraction.
      integer, intent(in), optional :: grid_count !! Maximum threshold candidates.
      logical, intent(in), optional :: momentum !! Use first-difference transitions.
      real(dp), intent(in), optional :: threshold_variable(:) !! External transition series.
      logical, intent(in), optional :: include_mean(:) !! Regime intercept flags.
      character(len=*), intent(in), optional :: criterion !! `SSR`, `AIC`, `BIC`, or `pooled-AIC`.
      integer, intent(in), optional :: order_candidates(:) !! Permitted AR orders in each regime.
      integer, intent(in), optional :: threshold_counts(:) !! Permitted threshold counts, one or two.
      logical, intent(in), optional :: same_lags !! Require the same AR order in every regime.
      character(len=*), intent(in), optional :: common !! `none`, `include`, `lags`, or `both`.
      logical, intent(in), optional :: outer_symmetric !! Share outer-regime coefficients.
      logical, intent(in), optional :: same_sample !! Score every model on a common sample.
      real(dp), intent(in), optional :: threshold_variables(:, :) !! Candidate external transition series.
      logical, intent(in), optional :: include_trend(:) !! Regime trend flags.
      character(len=*), intent(in), optional :: representation !! `level`, `diff`, or `ADF`.
      integer, intent(in), optional :: lag_spacing !! Spacing between candidate lag positions.
      integer, intent(in), optional :: forecast_step !! Lead from first predictor to response.
      integer, intent(in), optional :: lag_indices(:, :) !! Candidate-lag indices by regime.
      type(tsdyn_setar_selection_t) :: out
      type(nts_tar_model_t) :: candidate
      integer, allocatable :: counts(:), orders(:)
      logical, allocatable :: means(:)
      character(len=:), allocatable :: selected_criterion, common_mode
      real(dp) :: candidate_score
      integer :: candidate_count, candidate_index, count_index, delay
      integer :: regimes, combinations, combination, variable_index
      integer :: variable_count, order_count, common_start, code, regime
      logical :: equal_orders, symmetric, use_common_sample, use_momentum
      logical :: means_valid

      selected_criterion = 'aic'
      if (present(criterion)) selected_criterion = lowercase(adjustl(criterion))
      common_mode = 'none'
      if (present(common)) common_mode = lowercase(adjustl(common))
      equal_orders = .false.
      if (present(same_lags)) equal_orders = same_lags
      symmetric = .false.
      if (present(outer_symmetric)) symmetric = outer_symmetric
      use_common_sample = .false.
      if (present(same_sample)) use_common_sample = same_sample
      use_momentum = .false.
      if (present(momentum)) use_momentum = momentum
      if (present(threshold_counts)) then
         counts = threshold_counts
      else
         counts = [threshold_count]
      end if
      if (present(order_candidates)) then
         order_count = size(order_candidates)
      else
         order_count = 1
      end if
      variable_count = 1
      if (present(threshold_variables)) variable_count = size(threshold_variables, 2)
      if (maximum_delay < 1 .or. size(counts) < 1 .or. &
         any(counts < 1) .or. any(counts > 2) .or. order_count < 1 .or. &
         variable_count < 1 .or. &
         (present(threshold_variable) .and. present(threshold_variables)) .or. &
         (selected_criterion /= 'ssr' .and. selected_criterion /= 'aic' .and. &
         selected_criterion /= 'bic' .and. &
         selected_criterion /= 'pooled-aic') .or. &
         (common_mode /= 'none' .and. common_mode /= 'include' .and. &
         common_mode /= 'lags' .and. common_mode /= 'both')) then
         out%info = 1
         return
      end if
      if (present(order_candidates)) then
         if (any(order_candidates < 0)) then
            out%info = 1
            return
         end if
      end if
      if (present(threshold_variables)) then
         if (size(threshold_variables, 1) /= size(series)) then
            out%info = 1
            return
         end if
      end if
      candidate_count = 0
      do count_index = 1, size(counts)
         regimes = counts(count_index) + 1
         combinations = 1
         if (present(order_candidates)) then
            combinations = order_count
            if (.not. equal_orders) combinations = order_count**regimes
         end if
         candidate_count = candidate_count + maximum_delay*variable_count*combinations
      end do
      allocate(out%delay(candidate_count), out%threshold_count(candidate_count))
      allocate(out%ar_order(3, candidate_count))
      allocate(out%threshold_variable_index(candidate_count))
      allocate(out%thresholds(2, candidate_count), out%score(candidate_count))
      out%score = huge(1.0_dp)
      out%ar_order = 0
      out%thresholds = 0.0_dp
      out%criterion = selected_criterion
      out%common = common_mode
      out%momentum = use_momentum
      out%same_sample = use_common_sample
      common_start = setar_selection_common_start(ar_order, maximum_delay, &
         use_momentum, present(threshold_variable) .or. &
         present(threshold_variables), order_candidates)
      candidate_index = 0
      do count_index = 1, size(counts)
         regimes = counts(count_index) + 1
         combinations = 1
         if (present(order_candidates)) then
            combinations = order_count
            if (.not. equal_orders) combinations = order_count**regimes
         end if
         allocate(orders(regimes), means(regimes))
         means = .true.
         means_valid = .true.
         if (present(include_mean)) then
            means_valid = size(include_mean) >= regimes
            if (means_valid) means = include_mean(:regimes)
         end if
         do combination = 1, combinations
            if (present(order_candidates)) then
               code = combination - 1
               do regime = 1, regimes
                  orders(regime) = order_candidates(mod(code, order_count) + 1)
                  if (equal_orders) then
                     orders = orders(regime)
                     exit
                  end if
                  code = code/order_count
               end do
            else if (size(ar_order) == regimes) then
               orders = ar_order
            else
               orders = -1
            end if
            do variable_index = 1, variable_count
               do delay = 1, maximum_delay
                  candidate_index = candidate_index + 1
                  out%delay(candidate_index) = delay
                  out%threshold_count(candidate_index) = counts(count_index)
                  out%ar_order(:regimes, candidate_index) = orders
                  out%threshold_variable_index(candidate_index) = 0
                  if (present(threshold_variables)) then
                     out%threshold_variable_index(candidate_index) = variable_index
                  else if (present(threshold_variable)) then
                     out%threshold_variable_index(candidate_index) = 1
                  end if
                  if (any(orders < 0) .or. .not. means_valid .or. &
                     (symmetric .and. regimes /= 3)) cycle
                  if (present(threshold_variables)) then
                     candidate = setar_selection_fit(series, orders, counts(count_index), &
                        delay, trim, grid_count, use_momentum, means, common_mode, &
                        symmetric, threshold_variables(:, variable_index), &
                        include_trend, representation, lag_spacing, forecast_step, &
                        lag_indices)
                  else if (present(threshold_variable)) then
                     candidate = setar_selection_fit(series, orders, counts(count_index), &
                        delay, trim, grid_count, use_momentum, means, common_mode, &
                        symmetric, threshold_variable, include_trend, representation, &
                        lag_spacing, forecast_step, lag_indices)
                  else
                     candidate = setar_selection_fit(series, orders, counts(count_index), &
                        delay, trim, grid_count, use_momentum, means, common_mode, &
                        symmetric, include_trend=include_trend, &
                        representation=representation, lag_spacing=lag_spacing, &
                        forecast_step=forecast_step, lag_indices=lag_indices)
                  end if
                  if (candidate%info /= 0) cycle
                  out%thresholds(:counts(count_index), candidate_index) = &
                     candidate%thresholds
                  if (present(threshold_variables)) then
                     candidate_score = setar_selection_score(candidate, selected_criterion, &
                        common_mode, symmetric, use_common_sample, common_start, &
                        use_momentum, threshold_variables(:, variable_index))
                  else if (present(threshold_variable)) then
                     candidate_score = setar_selection_score(candidate, selected_criterion, &
                        common_mode, symmetric, use_common_sample, common_start, &
                        use_momentum, threshold_variable)
                  else
                     candidate_score = setar_selection_score(candidate, selected_criterion, &
                        common_mode, symmetric, use_common_sample, common_start, use_momentum)
                  end if
                  out%score(candidate_index) = candidate_score
                  if (out%selected == 0 .or. &
                     candidate_score < out%score(out%selected)) then
                     out%selected = candidate_index
                     out%model = candidate
                  end if
               end do
            end do
         end do
         deallocate(orders, means)
      end do
      if (out%selected == 0) out%info = 2
   end function tsdyn_setar_select

   pure function setar_selection_fit(series, ar_order, threshold_count, delay, &
      trim, grid_count, momentum, include_mean, common, outer_symmetric, &
      threshold_variable, include_trend, representation, lag_spacing, &
      forecast_step, lag_indices) result(model)
      !! Fit one candidate SETAR specification and apply requested restrictions.
      real(dp), intent(in) :: series(:) !! Observed univariate series.
      integer, intent(in) :: ar_order(:) !! Regime-specific autoregressive orders.
      integer, intent(in) :: threshold_count !! Number of thresholds.
      integer, intent(in) :: delay !! Transition delay.
      real(dp), intent(in), optional :: trim !! Minimum outer-regime fraction.
      integer, intent(in), optional :: grid_count !! Maximum threshold candidates.
      logical, intent(in) :: momentum !! Use first-difference transitions.
      logical, intent(in) :: include_mean(:) !! Regime intercept flags.
      character(len=*), intent(in) :: common !! Parameter-sharing mode.
      logical, intent(in) :: outer_symmetric !! Share outer-regime coefficients.
      real(dp), intent(in), optional :: threshold_variable(:) !! External transition series.
      logical, intent(in), optional :: include_trend(:) !! Regime trend flags.
      character(len=*), intent(in), optional :: representation !! Model representation.
      integer, intent(in), optional :: lag_spacing !! Candidate-lag spacing.
      integer, intent(in), optional :: forecast_step !! Predictor-to-response lead.
      integer, intent(in), optional :: lag_indices(:, :) !! Candidate-lag indices by regime.
      type(nts_tar_model_t) :: model

      model = tsdyn_setar_fit(series, ar_order, threshold_count, delay, trim, &
         grid_count, momentum, threshold_variable, include_mean, include_trend, &
         representation, lag_spacing, forecast_step, lag_indices)
      if (model%info /= 0) return
      if (common /= 'none' .or. outer_symmetric) then
         if (present(include_trend) .or. present(representation) .or. &
            present(lag_spacing) .or. present(forecast_step) .or. &
            present(lag_indices)) then
            model%info = 90
            return
         end if
         model = tsdyn_setar_restricted_fit(series, ar_order, model%thresholds, &
            delay, common, outer_symmetric, momentum, threshold_variable, include_mean)
      end if
   end function setar_selection_fit

   pure integer function setar_selection_common_start(ar_order, maximum_delay, &
      momentum, external, order_candidates) result(first)
      !! Find the first response shared by every candidate specification.
      integer, intent(in) :: ar_order(:) !! Baseline regime-specific AR orders.
      integer, intent(in) :: maximum_delay !! Largest transition delay.
      logical, intent(in) :: momentum !! Use first-difference transitions.
      logical, intent(in) :: external !! External transition variables are supplied.
      integer, intent(in), optional :: order_candidates(:) !! Permitted AR orders.
      integer :: maximum_order, transition_start

      maximum_order = maxval(ar_order)
      if (present(order_candidates)) maximum_order = maxval(order_candidates)
      transition_start = maximum_delay + 1
      if (momentum .and. .not. external) transition_start = transition_start + 1
      first = max(maximum_order + 1, transition_start)
   end function setar_selection_common_start

   pure real(dp) function setar_selection_score(model, criterion, common, &
      outer_symmetric, same_sample, common_start, momentum, threshold_variable) &
      result(score)
      !! Score a fitted SETAR on its estimation sample or a common sample.
      type(nts_tar_model_t), intent(in) :: model !! Fitted SETAR candidate.
      character(len=*), intent(in) :: criterion !! Requested selection criterion.
      character(len=*), intent(in) :: common !! Parameter-sharing mode.
      logical, intent(in) :: outer_symmetric !! Share outer-regime coefficients.
      logical, intent(in) :: same_sample !! Use the common comparison sample.
      integer, intent(in) :: common_start !! First common response index.
      logical, intent(in) :: momentum !! Use first-difference transitions.
      real(dp), intent(in), optional :: threshold_variable(:) !! External transition series.
      real(dp), allocatable :: transition(:), regime_rss(:)
      integer, allocatable :: regime_count(:)
      real(dp) :: prediction, residual, total_rss
      integer :: first_transition, first, time, regime, lag_index
      integer :: observation_count, parameter_count

      score = huge(1.0_dp)
      call setar_transition_series(model%data, model%delay, momentum, &
         threshold_variable, transition, first_transition)
      if (.not. allocated(transition)) return
      first = max(first_transition, maxval(model%ar_order) + 1)
      if (same_sample) first = max(first, common_start)
      observation_count = size(model%data) - first + 1
      if (observation_count < 1) return
      allocate(regime_rss(size(model%ar_order)), regime_count(size(model%ar_order)))
      regime_rss = 0.0_dp
      regime_count = 0
      do time = first, size(model%data)
         regime = tvecm_regime(transition(time - model%delay), model%thresholds)
         prediction = model%coefficients(regime, 1)
         if (allocated(model%include_trend)) then
            if (model%include_trend(regime)) prediction = prediction + &
               model%trend_coefficient(regime)*real(time - first + 1, dp)
         end if
         if (model%representation == 'adf') prediction = prediction + &
            model%adf_coefficient(regime)*model%data(time - model%forecast_step)
         do lag_index = 1, model%ar_order(regime)
            if (allocated(model%lag_active)) then
               if (.not. model%lag_active(regime, lag_index)) cycle
            end if
            if (model%representation == 'level') then
               prediction = prediction + &
                  model%coefficients(regime, lag_index + 1)* &
                  model%data(time - lag_index)
            else
               prediction = prediction + &
                  model%coefficients(regime, lag_index + 1)* &
                  (model%data(time - lag_index) - &
                  model%data(time - lag_index - 1))
            end if
         end do
         if (model%representation == 'level') then
            residual = model%data(time) - prediction
         else
            residual = model%data(time) - model%data(time - 1) - prediction
         end if
         regime_rss(regime) = regime_rss(regime) + residual**2
         regime_count(regime) = regime_count(regime) + 1
      end do
      if (any(regime_count < 1) .or. any(regime_rss <= 0.0_dp)) return
      total_rss = sum(regime_rss)
      parameter_count = setar_model_parameter_count(model, common, outer_symmetric)
      select case (criterion)
      case ('ssr')
         score = total_rss
      case ('aic', 'bic')
         score = real(observation_count, dp)*log(total_rss/ &
            real(observation_count, dp))
         parameter_count = parameter_count + size(model%thresholds)
         if (criterion == 'aic') then
            score = score + 2.0_dp*real(parameter_count, dp)
         else
            score = score + log(real(observation_count, dp))* &
               real(parameter_count, dp)
         end if
      case ('pooled-aic')
         score = 0.0_dp
         do regime = 1, size(regime_count)
            score = score + real(regime_count(regime), dp)* &
               log(regime_rss(regime)/real(regime_count(regime), dp))
         end do
         score = score + 2.0_dp*real(parameter_count + size(regime_count), dp)
      end select
   end function setar_selection_score

   pure integer function setar_selection_parameter_count(ar_order, include_mean, &
      common, outer_symmetric) result(count_parameters)
      !! Count free SETAR regression coefficients under parameter restrictions.
      integer, intent(in) :: ar_order(:) !! Regime-specific autoregressive orders.
      logical, intent(in) :: include_mean(:) !! Regime intercept flags.
      character(len=*), intent(in) :: common !! Parameter-sharing mode.
      logical, intent(in) :: outer_symmetric !! Share outer-regime coefficients.
      logical :: share_include, share_lags

      share_include = common == 'include' .or. common == 'both'
      share_lags = common == 'lags' .or. common == 'both'
      if (share_include) then
         count_parameters = merge(1, 0, include_mean(1))
      else
         count_parameters = count(include_mean)
         if (outer_symmetric .and. include_mean(1) .and. &
            include_mean(size(include_mean))) count_parameters = count_parameters - 1
      end if
      if (share_lags) then
         count_parameters = count_parameters + maxval(ar_order)
      else
         count_parameters = count_parameters + sum(ar_order)
         if (outer_symmetric) count_parameters = count_parameters - &
            min(ar_order(1), ar_order(size(ar_order)))
      end if
   end function setar_selection_parameter_count

   pure integer function setar_model_parameter_count(model, common, &
      outer_symmetric) result(count_parameters)
      !! Count free coefficients stored by a fitted generalized SETAR model.
      type(nts_tar_model_t), intent(in) :: model !! Fitted SETAR model.
      character(len=*), intent(in) :: common !! Parameter-sharing mode.
      logical, intent(in) :: outer_symmetric !! Share outer-regime coefficients.

      if (.not. allocated(model%lag_count)) then
         count_parameters = setar_selection_parameter_count(model%ar_order, &
            model%include_mean, common, outer_symmetric)
         return
      end if
      count_parameters = sum(model%lag_count) + count(model%include_mean)
      if (allocated(model%include_trend)) count_parameters = count_parameters + &
         count(model%include_trend)
      if (model%representation == 'adf') count_parameters = count_parameters + &
         size(model%ar_order)
   end function setar_model_parameter_count

   pure function tsdyn_setar_forecast_draws(model, origin, normal_draws, &
      level) result(out)
      !! Produce SETAR forecasts from supplied standard-normal draws.
      type(nts_tar_model_t), intent(in) :: model !! Fitted SETAR model.
      integer, intent(in) :: origin !! Forecast origin in the fitted series.
      real(dp), intent(in) :: normal_draws(:, :) !! Horizon-by-simulation normal draws.
      real(dp), intent(in), optional :: level !! Central interval coverage.
      type(nts_tar_forecast_t) :: out

      if (allocated(model%lag_active)) then
         out = setar_general_forecast_draws(model, origin, normal_draws, level)
      else
         out = nts_utar_forecast_draws(model, origin, normal_draws, level)
      end if
   end function tsdyn_setar_forecast_draws

   pure function setar_general_forecast_draws(model, origin, normal_draws, &
      level) result(out)
      !! Forecast a generalized SETAR from supplied standard-normal draws.
      type(nts_tar_model_t), intent(in) :: model !! Fitted generalized SETAR model.
      integer, intent(in) :: origin !! Forecast origin in the fitted series.
      real(dp), intent(in) :: normal_draws(:, :) !! Horizon-by-simulation normal draws.
      real(dp), intent(in), optional :: level !! Central interval coverage.
      type(nts_tar_forecast_t) :: out
      real(dp), allocatable :: work(:), ordered(:)
      real(dp) :: selected_level, prediction, transition_value
      integer :: horizon, simulations, simulation, step, time, regime, lag
      integer :: first_fitted

      selected_level = 0.95_dp
      if (present(level)) selected_level = level
      horizon = size(normal_draws, 1)
      simulations = size(normal_draws, 2)
      if (model%info /= 0 .or. origin < maxval(model%ar_order) + 1 .or. &
         origin > size(model%data) .or. horizon < 1 .or. simulations < 1 .or. &
         selected_level <= 0.0_dp .or. selected_level >= 1.0_dp) then
         out%info = 1
         return
      end if
      allocate(out%simulations(horizon, simulations), out%mean(horizon))
      allocate(out%lower(horizon), out%upper(horizon))
      first_fitted = size(model%data) - size(model%residuals) + 1
      do simulation = 1, simulations
         allocate(work(origin + horizon))
         work(:origin) = model%data(:origin)
         do step = 1, horizon
            time = origin + step
            transition_value = work(time - model%delay)
            regime = tvecm_regime(transition_value, model%thresholds)
            prediction = model%coefficients(regime, 1)
            if (model%include_trend(regime)) prediction = prediction + &
               model%trend_coefficient(regime)*real(time - first_fitted + 1, dp)
            if (model%representation == 'adf') prediction = prediction + &
               model%adf_coefficient(regime)*work(time - model%forecast_step)
            do lag = 1, model%ar_order(regime)
               if (.not. model%lag_active(regime, lag)) cycle
               if (model%representation == 'level') then
                  prediction = prediction + model%coefficients(regime, lag + 1)* &
                     work(time - lag)
               else
                  prediction = prediction + model%coefficients(regime, lag + 1)* &
                     (work(time - lag) - work(time - lag - 1))
               end if
            end do
            prediction = prediction + model%innovation_sd(regime)* &
               normal_draws(step, simulation)
            if (model%representation == 'level') then
               work(time) = prediction
            else
               work(time) = work(time - 1) + prediction
            end if
            out%simulations(step, simulation) = work(time)
         end do
         deallocate(work)
      end do
      do step = 1, horizon
         ordered = sorted(out%simulations(step, :))
         out%mean(step) = sum(out%simulations(step, :))/real(simulations, dp)
         out%lower(step) = quantile(ordered, 0.5_dp*(1.0_dp - selected_level))
         out%upper(step) = quantile(ordered, 0.5_dp*(1.0_dp + selected_level))
      end do
      out%level = selected_level
      out%origin = origin
   end function setar_general_forecast_draws

   pure function tsdyn_tvar_fit(series, ar_order, threshold_count, &
      threshold_component, delay, trim, grid_count, momentum, &
      threshold_variable, include_mean, include_trend) result(model)
      !! Fit a one- or two-threshold vector autoregression by grid search.
      real(dp), intent(in) :: series(:, :) !! Time-by-variable observations.
      integer, intent(in) :: ar_order(:) !! Regime-specific VAR orders.
      integer, intent(in) :: threshold_count !! Number of thresholds, one or two.
      integer, intent(in) :: threshold_component !! Variable driving regime selection.
      integer, intent(in) :: delay !! Delay applied to the transition variable.
      real(dp), intent(in), optional :: trim !! Minimum outer-regime fraction.
      integer, intent(in), optional :: grid_count !! Maximum threshold grid size.
      logical, intent(in), optional :: momentum !! Use delayed first differences.
      real(dp), intent(in), optional :: threshold_variable(:) !! External transition series.
      logical, intent(in), optional :: include_mean(:) !! Regime intercept flags.
      logical, intent(in), optional :: include_trend(:) !! Regime trend flags.
      type(nts_mtar_model_t) :: model
      type(nts_mtar_model_t) :: candidate_model
      real(dp), allocatable :: transition(:), ordered(:), grid(:), inverse(:, :)
      real(dp) :: selected_trim, score, best_score, logdet
      integer :: selected_grid, start, usable, i, j, regime, status
      logical :: use_momentum

      selected_trim = 0.1_dp
      if (present(trim)) selected_trim = trim
      selected_grid = 30
      if (present(grid_count)) selected_grid = grid_count
      use_momentum = .false.
      if (present(momentum)) use_momentum = momentum
      if (threshold_count < 1 .or. threshold_count > 2 .or. &
         size(ar_order) /= threshold_count + 1 .or. threshold_component < 1 .or. &
         threshold_component > size(series, 2) .or. selected_trim <= 0.0_dp .or. &
         selected_trim >= 1.0_dp/real(threshold_count + 1, dp) .or. &
         selected_grid < 2 .or. delay < 1) then
         model%info = 1
         return
      end if
      call setar_transition_series(series(:, threshold_component), delay, &
         use_momentum, threshold_variable, transition, start)
      if (.not. allocated(transition)) then
         model%info = 2
         return
      end if
      usable = size(series, 1) - max(start, maxval(ar_order) + 1) + 1
      if (usable < 10) then
         model%info = 3
         return
      end if
      ordered = sorted(transition(size(transition) - usable + 1:))
      call threshold_candidate_grid(ordered, selected_trim, selected_grid, grid)
      best_score = huge(1.0_dp)
      do i = 1, size(grid) - threshold_count + 1
         do j = i + threshold_count - 1, merge(i, size(grid), threshold_count == 1)
            if (threshold_count == 1 .and. j /= i) cycle
            if (threshold_count == 1) then
               candidate_model = tsdyn_tvar_restricted_fit(series, ar_order, &
                  grid(i:i), threshold_component, delay, momentum=use_momentum, &
                  threshold_variable=threshold_variable, include_mean=include_mean, &
                  include_trend=include_trend)
            else
               candidate_model = tsdyn_tvar_restricted_fit(series, ar_order, &
                  [grid(i), grid(j)], threshold_component, delay, &
                  momentum=use_momentum, threshold_variable=threshold_variable, &
                  include_mean=include_mean, include_trend=include_trend)
            end if
            if (candidate_model%info /= 0) cycle
            if (minval(real(candidate_model%regime_observations, dp)) < &
               selected_trim*real(size(candidate_model%residuals, 1), dp)) cycle
            score = 0.0_dp
            allocate(inverse(size(series, 2), size(series, 2)))
            do regime = 1, size(candidate_model%ar_order)
               call inverse_logdet(candidate_model%covariance(:, :, regime), &
                  inverse, logdet, status, 1.0e-12_dp)
               if (status /= 0) then
                  score = huge(1.0_dp)
                  exit
               end if
               score = score + real(candidate_model%regime_observations(regime), dp)*logdet
            end do
            deallocate(inverse)
            if (score < best_score) then
               best_score = score
               model = candidate_model
            end if
         end do
      end do
      if (best_score == huge(1.0_dp)) model%info = 5
   end function tsdyn_tvar_fit

   pure function tsdyn_tvar_restricted_fit(series, ar_order, thresholds, &
      threshold_component, delay, common_intercept, momentum, &
      threshold_variable, include_mean, include_trend) result(model)
      !! Fit a fixed-threshold TVAR with an optional common intercept.
      real(dp), intent(in) :: series(:, :) !! Time-by-variable observations.
      integer, intent(in) :: ar_order(:) !! Regime-specific VAR orders.
      real(dp), intent(in) :: thresholds(:) !! Ordered fixed thresholds.
      integer, intent(in) :: threshold_component !! Variable driving regime selection.
      integer, intent(in) :: delay !! Delay applied to the transition variable.
      logical, intent(in), optional :: common_intercept !! Share one intercept across regimes.
      logical, intent(in), optional :: momentum !! Use delayed first differences for switching.
      real(dp), intent(in), optional :: threshold_variable(:) !! External transition series.
      logical, intent(in), optional :: include_mean(:) !! Regime intercept flags.
      logical, intent(in), optional :: include_trend(:) !! Regime trend flags.
      type(nts_mtar_model_t) :: model
      real(dp), allocatable :: transition(:), design(:, :), response(:, :)
      real(dp), allocatable :: coefficient(:, :), fitted(:, :), residual(:, :)
      real(dp), allocatable :: inverse(:, :), regime_residual(:, :)
      integer, allocatable :: map(:, :), row_regime(:)
      logical, allocatable :: means(:), trends(:)
      real(dp) :: rss, logdet
      integer :: variables, regimes, maximum_order, first_transition, start, rows
      integer :: row, time, regime, lag_index, variable, component, column, status
      logical :: shared, use_momentum

      variables = size(series, 2)
      regimes = size(thresholds) + 1
      maximum_order = maxval(ar_order)
      shared = .false.
      if (present(common_intercept)) shared = common_intercept
      use_momentum = .false.
      if (present(momentum)) use_momentum = momentum
      if (variables < 1 .or. regimes < 2 .or. size(ar_order) /= regimes .or. &
         any(ar_order < 0) .or. threshold_component < 1 .or. &
         threshold_component > variables .or. delay < 1) then
         model%info = 1
         return
      end if
      if (size(thresholds) > 1) then
         if (any(thresholds(2:) <= thresholds(:size(thresholds) - 1))) then
            model%info = 1
            return
         end if
      end if
      allocate(means(regimes), trends(regimes))
      means = .true.
      trends = .false.
      if (present(include_mean)) then
         if (size(include_mean) /= regimes) then
            model%info = 2
            return
         end if
         means = include_mean
      end if
      if (present(include_trend)) then
         if (size(include_trend) /= regimes) then
            model%info = 2
            return
         end if
         trends = include_trend
      end if
      if (shared .and. any(means .neqv. means(1))) then
         model%info = 2
         return
      end if
      call setar_transition_series(series(:, threshold_component), delay, &
         use_momentum, threshold_variable, transition, first_transition)
      if (.not. allocated(transition)) then
         model%info = 3
         return
      end if
      start = max(maximum_order + 1, first_transition)
      rows = size(series, 1) - start + 1
      allocate(map(regimes, 2 + variables*maximum_order))
      map = 0
      column = 0
      if (shared .and. means(1)) then
         column = 1
         map(:, 1) = 1
      else
         do regime = 1, regimes
            if (.not. means(regime)) cycle
            column = column + 1
            map(regime, 1) = column
         end do
      end if
      do regime = 1, regimes
         if (trends(regime)) then
            column = column + 1
            map(regime, 2) = column
         end if
      end do
      do regime = 1, regimes
         do lag_index = 1, ar_order(regime)
            do variable = 1, variables
               component = 2 + (lag_index - 1)*variables + variable
               column = column + 1
               map(regime, component) = column
            end do
         end do
      end do
      if (rows <= column) then
         model%info = 4
         return
      end if
      allocate(design(rows, column), response(rows, variables), row_regime(rows))
      design = 0.0_dp
      do row = 1, rows
         time = start + row - 1
         regime = tvecm_regime(transition(time - delay), thresholds)
         row_regime(row) = regime
         response(row, :) = series(time, :)
         if (map(regime, 1) > 0) design(row, map(regime, 1)) = 1.0_dp
         if (map(regime, 2) > 0) design(row, map(regime, 2)) = real(row, dp)
         do lag_index = 1, ar_order(regime)
            do variable = 1, variables
               component = 2 + (lag_index - 1)*variables + variable
               design(row, map(regime, component)) = series(time - lag_index, variable)
            end do
         end do
      end do
      call tvecm_multivariate_regression(design, response, coefficient, fitted, &
         residual, rss, status)
      if (status /= 0) then
         model%info = 10 + status
         return
      end if
      allocate(model%intercept(variables, regimes))
      allocate(model%trend(variables, regimes))
      allocate(model%ar(variables, variables, maximum_order, regimes))
      allocate(model%covariance(variables, variables, regimes))
      allocate(model%regime_observations(regimes))
      model%intercept = 0.0_dp
      model%trend = 0.0_dp
      model%ar = 0.0_dp
      model%aic = 2.0_dp*real(variables*column, dp)
      do regime = 1, regimes
         if (map(regime, 1) > 0) then
            model%intercept(:, regime) = coefficient(:, map(regime, 1))
         end if
         if (map(regime, 2) > 0) then
            model%trend(:, regime) = coefficient(:, map(regime, 2))
         end if
         do lag_index = 1, ar_order(regime)
            do variable = 1, variables
               component = 2 + (lag_index - 1)*variables + variable
               model%ar(:, variable, lag_index, regime) = &
                  coefficient(:, map(regime, component))
            end do
         end do
         model%regime_observations(regime) = count(row_regime == regime)
         if (model%regime_observations(regime) <= variables) then
            model%info = 5
            return
         end if
         regime_residual = reshape(pack(residual, spread(row_regime == regime, &
            2, variables)), [model%regime_observations(regime), variables])
         model%covariance(:, :, regime) = matmul(transpose(regime_residual), &
            regime_residual)/real(model%regime_observations(regime), dp)
         allocate(inverse(variables, variables))
         call inverse_logdet(model%covariance(:, :, regime), inverse, logdet, &
            status, 1.0e-12_dp)
         deallocate(inverse)
         if (status /= 0) then
            model%info = 20 + regime
            return
         end if
         model%aic = model%aic + real(model%regime_observations(regime), dp)*logdet
      end do
      model%data = series
      model%residuals = residual
      allocate(model%standardized_residuals(rows, variables))
      do row = 1, rows
         regime = row_regime(row)
         do variable = 1, variables
            model%standardized_residuals(row, variable) = residual(row, variable)/ &
               sqrt(model%covariance(variable, variable, regime))
         end do
      end do
      model%ar_order = ar_order
      model%thresholds = thresholds
      model%include_mean = means
      model%include_trend = trends
      model%threshold_component = threshold_component
      model%delay = delay
   end function tsdyn_tvar_restricted_fit

   pure function tsdyn_tvar_simulate_from_standard(model, standard_draws, &
      burnin) result(out)
      !! Simulate a fitted TVAR specification from supplied standard-normal draws.
      type(nts_mtar_model_t), intent(in) :: model !! TVAR parameter specification.
      real(dp), intent(in) :: standard_draws(:, :) !! Time-by-variable normal draws.
      integer, intent(in), optional :: burnin !! Initial observations to discard.
      type(nts_mtar_simulation_t) :: out
      integer :: discard

      discard = 0
      if (present(burnin)) discard = burnin
      if (model%info /= 0) then
         out%info = 1
         return
      end if
      if (allocated(model%include_trend)) then
         if (any(model%include_trend)) then
            out = tvar_trend_simulate_from_standard(model, standard_draws, discard)
            return
         end if
      end if
      out = nts_mtar_simulate_from_standard(model%intercept, model%ar, &
         model%ar_order, model%covariance, model%thresholds, &
         model%threshold_component, model%delay, standard_draws, discard)
   end function tsdyn_tvar_simulate_from_standard

   pure function tvar_trend_simulate_from_standard(model, standard_draws, &
      burnin) result(out)
      !! Simulate a trend-enabled TVAR from supplied standard-normal draws.
      type(nts_mtar_model_t), intent(in) :: model !! Trend-enabled TVAR specification.
      real(dp), intent(in) :: standard_draws(:, :) !! Time-by-variable standard normals.
      integer, intent(in) :: burnin !! Number of initial observations to discard.
      type(nts_mtar_simulation_t) :: out
      real(dp), allocatable :: work(:, :), errors(:, :), lower(:, :, :)
      real(dp), allocatable :: factor(:, :)
      integer :: variables, regimes, total, maximum_order, start
      integer :: time, regime, lag, status

      variables = size(model%intercept, 1)
      regimes = size(model%intercept, 2)
      total = size(standard_draws, 1)
      maximum_order = max(1, maxval(model%ar_order))
      start = max(maximum_order, model%delay) + 1
      if (burnin < 0 .or. burnin >= total .or. total < start .or. &
         size(standard_draws, 2) /= variables .or. &
         .not. all(ieee_is_finite(standard_draws))) then
         out%info = 1
         return
      end if
      allocate(lower(variables, variables, regimes))
      do regime = 1, regimes
         call cholesky_lower(model%covariance(:, :, regime), factor, status)
         if (status /= 0) then
            out%info = 10 + regime
            return
         end if
         lower(:, :, regime) = factor
         deallocate(factor)
      end do
      allocate(work(total, variables), errors(total, variables))
      work = 0.0_dp
      errors = 0.0_dp
      do time = 1, start - 1
         errors(time, :) = matmul(lower(:, :, 1), standard_draws(time, :))
         work(time, :) = model%intercept(:, 1) + errors(time, :)
      end do
      do time = start, total
         regime = tvecm_regime(work(time - model%delay, &
            model%threshold_component), model%thresholds)
         errors(time, :) = matmul(lower(:, :, regime), standard_draws(time, :))
         work(time, :) = model%intercept(:, regime) + errors(time, :)
         if (model%include_trend(regime)) work(time, :) = work(time, :) + &
            model%trend(:, regime)*real(time - start + 1, dp)
         do lag = 1, model%ar_order(regime)
            work(time, :) = work(time, :) + &
               matmul(model%ar(:, :, lag, regime), work(time - lag, :))
         end do
      end do
      out%series = work(burnin + 1:, :)
      out%innovations = errors(burnin + 1:, :)
      out%intercept = model%intercept
      out%trend = model%trend
      out%ar = model%ar
      out%ar_order = model%ar_order
      out%covariance = model%covariance
      out%thresholds = model%thresholds
      allocate(out%regime_observations(regimes))
      out%regime_observations = 0
      do time = burnin + 1, total
         regime = 1
         if (time > model%delay) regime = tvecm_regime( &
            work(time - model%delay, model%threshold_component), &
            model%thresholds)
         out%regime_observations(regime) = out%regime_observations(regime) + 1
      end do
      out%threshold_component = model%threshold_component
      out%delay = model%delay
      out%burnin = burnin
   end function tvar_trend_simulate_from_standard

   pure function tsdyn_tvar_select(series, ar_order, threshold_count, &
      maximum_delay, threshold_components, trim, grid_count, momentum, &
      include_mean, criterion, order_candidates, threshold_counts, &
      same_lags, common_intercept, same_sample, threshold_variables, &
      include_trend) result(out)
      !! Select a TVAR specification over thresholds, variables, delays, and orders.
      real(dp), intent(in) :: series(:, :) !! Time-by-variable observations.
      integer, intent(in) :: ar_order(:) !! Regime-specific VAR orders.
      integer, intent(in) :: threshold_count !! Number of thresholds.
      integer, intent(in) :: maximum_delay !! Largest transition delay considered.
      integer, intent(in), optional :: threshold_components(:) !! Candidate transition variables.
      real(dp), intent(in), optional :: trim !! Minimum outer-regime fraction.
      integer, intent(in), optional :: grid_count !! Maximum threshold candidates.
      logical, intent(in), optional :: momentum !! Use first-difference transitions.
      logical, intent(in), optional :: include_mean(:) !! Regime intercept flags.
      character(len=*), intent(in), optional :: criterion !! `SSR`, `AIC`, or `BIC`.
      integer, intent(in), optional :: order_candidates(:) !! Permitted VAR orders in each regime.
      integer, intent(in), optional :: threshold_counts(:) !! Permitted threshold counts, one or two.
      logical, intent(in), optional :: same_lags !! Require the same VAR order in every regime.
      logical, intent(in), optional :: common_intercept !! Share one intercept across regimes.
      logical, intent(in), optional :: same_sample !! Score every model on a common sample.
      real(dp), intent(in), optional :: threshold_variables(:, :) !! Candidate external transition series.
      logical, intent(in), optional :: include_trend(:) !! Regime trend flags.
      type(tsdyn_tvar_selection_t) :: out
      type(nts_mtar_model_t) :: candidate
      integer, allocatable :: components(:), counts(:), orders(:)
      logical, allocatable :: means(:)
      character(len=:), allocatable :: selected_criterion
      real(dp) :: candidate_score
      integer :: component_index, delay, candidate_index, candidates
      integer :: count_index, regimes, combinations, combination, code, regime
      integer :: order_count, external_count, transition_index, common_start
      logical :: equal_orders, shared_intercept, use_common_sample
      logical :: use_momentum, means_valid, has_internal

      selected_criterion = 'aic'
      if (present(criterion)) selected_criterion = lowercase(adjustl(criterion))
      equal_orders = .false.
      if (present(same_lags)) equal_orders = same_lags
      shared_intercept = .false.
      if (present(common_intercept)) shared_intercept = common_intercept
      use_common_sample = .false.
      if (present(same_sample)) use_common_sample = same_sample
      use_momentum = .false.
      if (present(momentum)) use_momentum = momentum
      if (present(threshold_counts)) then
         counts = threshold_counts
      else
         counts = [threshold_count]
      end if
      if (present(order_candidates)) then
         order_count = size(order_candidates)
      else
         order_count = 1
      end if
      external_count = 0
      if (present(threshold_variables)) external_count = size(threshold_variables, 2)
      if (maximum_delay < 1 .or. size(series, 2) < 1 .or. &
         size(counts) < 1 .or. any(counts < 1) .or. any(counts > 2) .or. &
         order_count < 1 .or. external_count < 0 .or. &
         (selected_criterion /= 'ssr' .and. selected_criterion /= 'aic' .and. &
         selected_criterion /= 'bic')) then
         out%info = 1
         return
      end if
      if (present(order_candidates)) then
         if (any(order_candidates < 0)) then
            out%info = 1
            return
         end if
      end if
      if (present(threshold_variables)) then
         if (size(threshold_variables, 1) /= size(series, 1)) then
            out%info = 1
            return
         end if
      end if
      if (present(threshold_components)) then
         if (size(threshold_components) > 0) then
            if (minval(threshold_components) < 1 .or. &
               maxval(threshold_components) > size(series, 2)) then
               out%info = 1
               return
            end if
         end if
         components = threshold_components
      else
         allocate(components(size(series, 2)))
         components = [(component_index, component_index=1, size(series, 2))]
      end if
      if (size(components) + external_count < 1) then
         out%info = 1
         return
      end if
      candidates = 0
      do count_index = 1, size(counts)
         regimes = counts(count_index) + 1
         combinations = 1
         if (present(order_candidates)) then
            combinations = order_count
            if (.not. equal_orders) combinations = order_count**regimes
         end if
         candidates = candidates + maximum_delay*combinations* &
            (size(components) + external_count)
      end do
      allocate(out%threshold_component(candidates))
      allocate(out%threshold_variable_index(candidates), out%delay(candidates))
      allocate(out%threshold_count(candidates), out%ar_order(3, candidates))
      allocate(out%thresholds(2, candidates), out%score(candidates))
      out%score = huge(1.0_dp)
      out%ar_order = 0
      out%thresholds = 0.0_dp
      out%criterion = selected_criterion
      out%momentum = use_momentum
      out%common_intercept = shared_intercept
      out%same_sample = use_common_sample
      has_internal = size(components) > 0
      common_start = tvar_selection_common_start(ar_order, maximum_delay, &
         use_momentum, has_internal, order_candidates)
      candidate_index = 0
      do count_index = 1, size(counts)
         regimes = counts(count_index) + 1
         combinations = 1
         if (present(order_candidates)) then
            combinations = order_count
            if (.not. equal_orders) combinations = order_count**regimes
         end if
         allocate(orders(regimes), means(regimes))
         means = .true.
         means_valid = .true.
         if (present(include_mean)) then
            means_valid = size(include_mean) >= regimes
            if (means_valid) means = include_mean(:regimes)
         end if
         do combination = 1, combinations
            if (present(order_candidates)) then
               code = combination - 1
               do regime = 1, regimes
                  orders(regime) = order_candidates(mod(code, order_count) + 1)
                  if (equal_orders) then
                     orders = orders(regime)
                     exit
                  end if
                  code = code/order_count
               end do
            else if (size(ar_order) == regimes) then
               orders = ar_order
            else
               orders = -1
            end if
            do transition_index = 1, size(components) + external_count
               do delay = 1, maximum_delay
                  candidate_index = candidate_index + 1
                  out%delay(candidate_index) = delay
                  out%threshold_count(candidate_index) = counts(count_index)
                  out%ar_order(:regimes, candidate_index) = orders
                  out%threshold_component(candidate_index) = 0
                  out%threshold_variable_index(candidate_index) = 0
                  if (transition_index <= size(components)) then
                     component_index = components(transition_index)
                     out%threshold_component(candidate_index) = component_index
                  else
                     component_index = 1
                     out%threshold_variable_index(candidate_index) = &
                        transition_index - size(components)
                  end if
                  if (any(orders < 0) .or. .not. means_valid) cycle
                  if (transition_index <= size(components)) then
                     candidate = tvar_selection_fit(series, orders, &
                        counts(count_index), component_index, delay, trim, &
                        grid_count, use_momentum, means, shared_intercept, &
                        include_trend=include_trend)
                  else
                     candidate = tvar_selection_fit(series, orders, &
                        counts(count_index), component_index, delay, trim, &
                        grid_count, use_momentum, means, shared_intercept, &
                        threshold_variables(:, out%threshold_variable_index(candidate_index)), &
                        include_trend)
                  end if
                  if (candidate%info /= 0) cycle
                  out%thresholds(:counts(count_index), candidate_index) = &
                     candidate%thresholds
                  if (transition_index <= size(components)) then
                     candidate_score = tvar_selection_score(candidate, &
                        selected_criterion, shared_intercept, use_common_sample, &
                        common_start, use_momentum)
                  else
                     candidate_score = tvar_selection_score(candidate, &
                        selected_criterion, shared_intercept, use_common_sample, &
                        common_start, use_momentum, threshold_variables(:, &
                        out%threshold_variable_index(candidate_index)))
                  end if
                  if (candidate_score >= huge(1.0_dp)) cycle
                  out%score(candidate_index) = candidate_score
                  if (out%selected == 0 .or. &
                     candidate_score < out%score(out%selected)) then
                     out%selected = candidate_index
                     out%model = candidate
                  end if
               end do
            end do
         end do
         deallocate(orders, means)
      end do
      if (out%selected == 0) out%info = 2
   end function tsdyn_tvar_select

   pure function tvar_selection_fit(series, ar_order, threshold_count, &
      threshold_component, delay, trim, grid_count, momentum, include_mean, &
      common_intercept, threshold_variable, include_trend) result(model)
      !! Fit one candidate TVAR specification and apply its intercept restriction.
      real(dp), intent(in) :: series(:, :) !! Time-by-variable observations.
      integer, intent(in) :: ar_order(:) !! Regime-specific VAR orders.
      integer, intent(in) :: threshold_count !! Number of thresholds.
      integer, intent(in) :: threshold_component !! Internal transition component.
      integer, intent(in) :: delay !! Transition delay.
      real(dp), intent(in), optional :: trim !! Minimum outer-regime fraction.
      integer, intent(in), optional :: grid_count !! Maximum threshold candidates.
      logical, intent(in) :: momentum !! Use first-difference transitions.
      logical, intent(in) :: include_mean(:) !! Regime intercept flags.
      logical, intent(in) :: common_intercept !! Share one intercept across regimes.
      real(dp), intent(in), optional :: threshold_variable(:) !! External transition series.
      logical, intent(in), optional :: include_trend(:) !! Regime trend flags.
      type(nts_mtar_model_t) :: model

      model = tsdyn_tvar_fit(series, ar_order, threshold_count, &
         threshold_component, delay, trim, grid_count, momentum, &
         threshold_variable, include_mean, include_trend)
      if (model%info /= 0) return
      if (common_intercept) then
         model = tsdyn_tvar_restricted_fit(series, ar_order, model%thresholds, &
            threshold_component, delay, common_intercept, momentum, &
            threshold_variable, include_mean, include_trend)
      end if
   end function tvar_selection_fit

   pure integer function tvar_selection_common_start(ar_order, maximum_delay, &
      momentum, has_internal, order_candidates) result(first)
      !! Find the first response shared by every TVAR candidate specification.
      integer, intent(in) :: ar_order(:) !! Baseline regime-specific VAR orders.
      integer, intent(in) :: maximum_delay !! Largest transition delay.
      logical, intent(in) :: momentum !! Use first-difference transitions.
      logical, intent(in) :: has_internal !! Internal transition candidates are included.
      integer, intent(in), optional :: order_candidates(:) !! Permitted VAR orders.
      integer :: maximum_order, transition_start

      maximum_order = maxval(ar_order)
      if (present(order_candidates)) maximum_order = maxval(order_candidates)
      transition_start = maximum_delay + 1
      if (momentum .and. has_internal) transition_start = transition_start + 1
      first = max(maximum_order + 1, transition_start)
   end function tvar_selection_common_start

   pure real(dp) function tvar_selection_score(model, criterion, &
      common_intercept, same_sample, common_start, momentum, &
      threshold_variable) result(score)
      !! Score a fitted TVAR using pooled multivariate residual covariance.
      type(nts_mtar_model_t), intent(in) :: model !! Fitted TVAR candidate.
      character(len=*), intent(in) :: criterion !! Requested selection criterion.
      logical, intent(in) :: common_intercept !! Share one intercept across regimes.
      logical, intent(in) :: same_sample !! Use the common comparison sample.
      integer, intent(in) :: common_start !! First common response index.
      logical, intent(in) :: momentum !! Use first-difference transitions.
      real(dp), intent(in), optional :: threshold_variable(:) !! External transition series.
      real(dp), allocatable :: transition(:), covariance(:, :), inverse(:, :)
      real(dp), allocatable :: residual(:), prediction(:)
      real(dp) :: log_determinant, residual_sum_squares
      integer :: first_transition, first, time, regime, lag_index
      integer :: observations, variables, parameter_count, status

      score = huge(1.0_dp)
      variables = size(model%data, 2)
      call setar_transition_series(model%data(:, model%threshold_component), &
         model%delay, momentum, threshold_variable, transition, first_transition)
      if (.not. allocated(transition)) return
      first = max(first_transition, maxval(model%ar_order) + 1)
      if (same_sample) first = max(first, common_start)
      observations = size(model%data, 1) - first + 1
      if (observations <= variables) return
      allocate(covariance(variables, variables), inverse(variables, variables))
      allocate(residual(variables), prediction(variables))
      covariance = 0.0_dp
      residual_sum_squares = 0.0_dp
      do time = first, size(model%data, 1)
         regime = tvecm_regime(transition(time - model%delay), model%thresholds)
         prediction = model%intercept(:, regime)
         if (allocated(model%trend)) then
            if (model%include_trend(regime)) prediction = prediction + &
               model%trend(:, regime)*real(time - first + 1, dp)
         end if
         do lag_index = 1, model%ar_order(regime)
            prediction = prediction + matmul(model%ar(:, :, lag_index, regime), &
               model%data(time - lag_index, :))
         end do
         residual = model%data(time, :) - prediction
         covariance = covariance + spread(residual, 2, variables)* &
            spread(residual, 1, variables)
         residual_sum_squares = residual_sum_squares + sum(residual**2)
      end do
      covariance = covariance/real(observations, dp)
      if (criterion == 'ssr') then
         score = residual_sum_squares
         return
      end if
      call inverse_logdet(covariance, inverse, log_determinant, status, 1.0e-12_dp)
      if (status /= 0) return
      parameter_count = tvar_selection_parameter_count(model%ar_order, &
         model%include_mean, variables, common_intercept) + size(model%thresholds)
      score = real(observations, dp)*log_determinant
      if (criterion == 'aic') then
         score = score + 2.0_dp*real(parameter_count, dp)
      else
         score = score + log(real(observations, dp))*real(parameter_count, dp)
      end if
   end function tvar_selection_score

   pure integer function tvar_selection_parameter_count(ar_order, include_mean, &
      variables, common_intercept) result(count_parameters)
      !! Count free TVAR regression coefficients under an intercept restriction.
      integer, intent(in) :: ar_order(:) !! Regime-specific VAR orders.
      logical, intent(in) :: include_mean(:) !! Regime intercept flags.
      integer, intent(in) :: variables !! Number of response variables.
      logical, intent(in) :: common_intercept !! Share one intercept across regimes.
      integer :: regressor_count

      if (common_intercept) then
         regressor_count = merge(1, 0, include_mean(1))
      else
         regressor_count = count(include_mean)
      end if
      regressor_count = regressor_count + variables*sum(ar_order)
      count_parameters = variables*regressor_count
   end function tvar_selection_parameter_count

   pure function tsdyn_tvar_forecast_draws(model, origin, standard_draws, &
      level) result(out)
      !! Produce TVAR forecasts from supplied standard-normal draws.
      type(nts_mtar_model_t), intent(in) :: model !! Fitted TVAR model.
      integer, intent(in) :: origin !! Forecast origin in the fitted series.
      real(dp), intent(in) :: standard_draws(:, :, :) !! Horizon-by-variable draws.
      real(dp), intent(in), optional :: level !! Central interval coverage.
      type(nts_mtar_forecast_t) :: out
      real(dp), allocatable :: transposed(:, :, :)
      integer :: simulation, step

      if (allocated(model%trend)) then
         out = tvar_general_forecast_draws(model, origin, standard_draws, level)
         return
      end if
      allocate(transposed(size(standard_draws, 2), size(standard_draws, 1), &
         size(standard_draws, 3)))
      do simulation = 1, size(standard_draws, 3)
         do step = 1, size(standard_draws, 1)
            transposed(:, step, simulation) = standard_draws(step, :, simulation)
         end do
      end do
      out = nts_mtar_forecast_draws(model, origin, transposed, level)
   end function tsdyn_tvar_forecast_draws

   pure function tvar_general_forecast_draws(model, origin, standard_draws, &
      level) result(out)
      !! Forecast a deterministic TVAR from supplied multivariate normal draws.
      type(nts_mtar_model_t), intent(in) :: model !! Fitted TVAR model.
      integer, intent(in) :: origin !! Forecast origin in the fitted series.
      real(dp), intent(in) :: standard_draws(:, :, :) !! Horizon-variable-simulation draws.
      real(dp), intent(in), optional :: level !! Central interval coverage.
      type(nts_mtar_forecast_t) :: out
      real(dp), allocatable :: work(:, :), lower_factor(:, :, :), factor(:, :)
      real(dp), allocatable :: prediction(:), ordered(:)
      real(dp) :: selected_level
      integer :: horizon, variables, simulations, simulation, step, time
      integer :: regime, lag_index, status, first_fitted

      selected_level = 0.95_dp
      if (present(level)) selected_level = level
      horizon = size(standard_draws, 1)
      variables = size(standard_draws, 2)
      simulations = size(standard_draws, 3)
      if (model%info /= 0 .or. variables /= size(model%data, 2) .or. &
         origin < max(maxval(model%ar_order), model%delay) .or. &
         origin > size(model%data, 1) .or. horizon < 1 .or. simulations < 1 .or. &
         selected_level <= 0.0_dp .or. selected_level >= 1.0_dp) then
         out%info = 1
         return
      end if
      allocate(lower_factor(variables, variables, size(model%ar_order)))
      do regime = 1, size(model%ar_order)
         call cholesky_lower(model%covariance(:, :, regime), factor, status)
         if (status /= 0) then
            out%info = 2
            return
         end if
         lower_factor(:, :, regime) = factor
      end do
      allocate(out%simulations(horizon, variables, simulations))
      allocate(out%mean(horizon, variables), out%lower(horizon, variables))
      allocate(out%upper(horizon, variables), prediction(variables))
      first_fitted = size(model%data, 1) - size(model%residuals, 1) + 1
      do simulation = 1, simulations
         allocate(work(origin + horizon, variables))
         work(:origin, :) = model%data(:origin, :)
         do step = 1, horizon
            time = origin + step
            regime = tvecm_regime(work(time - model%delay, &
               model%threshold_component), model%thresholds)
            prediction = model%intercept(:, regime)
            if (model%include_trend(regime)) prediction = prediction + &
               model%trend(:, regime)*real(time - first_fitted + 1, dp)
            do lag_index = 1, model%ar_order(regime)
               prediction = prediction + matmul(model%ar(:, :, lag_index, regime), &
                  work(time - lag_index, :))
            end do
            prediction = prediction + matmul(lower_factor(:, :, regime), &
               standard_draws(step, :, simulation))
            work(time, :) = prediction
            out%simulations(step, :, simulation) = prediction
         end do
         deallocate(work)
      end do
      do step = 1, horizon
         do regime = 1, variables
            ordered = sorted(out%simulations(step, regime, :))
            out%mean(step, regime) = sum(out%simulations(step, regime, :))/ &
               real(simulations, dp)
            out%lower(step, regime) = quantile(ordered, &
               0.5_dp*(1.0_dp - selected_level))
            out%upper(step, regime) = quantile(ordered, &
               0.5_dp*(1.0_dp + selected_level))
         end do
      end do
      out%level = selected_level
      out%origin = origin
   end function tvar_general_forecast_draws

   pure function tsdyn_setar_linearity_test(series, order, delays, trim, &
      grid_count, momentum, bootstrap_series) result(out)
      !! Compare linear AR, two-regime SETAR, and three-regime SETAR fits.
      real(dp), intent(in) :: series(:) !! Observed univariate series.
      integer, intent(in) :: order !! Common autoregressive order.
      integer, intent(in) :: delays(:) !! Transition delays searched under each alternative.
      real(dp), intent(in), optional :: trim !! Minimum threshold-regime fraction.
      integer, intent(in), optional :: grid_count !! Maximum threshold candidates.
      logical, intent(in), optional :: momentum !! Use first-difference transitions.
      real(dp), intent(in), optional :: bootstrap_series(:, :) !! Null-bootstrap series by column.
      type(tsdyn_regime_count_test_t) :: out
      real(dp) :: measure(3), bootstrap_measure(3), bootstrap_statistics(3)
      integer :: replication, status

      call setar_regime_measures(series, order, delays, trim, grid_count, &
         momentum, measure, status)
      if (status /= 0) then
         out%info = status
         return
      end if
      out%fit_measure = measure
      call regime_count_statistics(measure, real(size(series) - order, dp), &
         .false., out%statistic)
      if (present(bootstrap_series)) then
         if (size(bootstrap_series, 1) /= size(series)) then
            out%info = 4
            return
         end if
         allocate(out%bootstrap_statistic(size(bootstrap_series, 2), 3))
         do replication = 1, size(bootstrap_series, 2)
            call setar_regime_measures(bootstrap_series(:, replication), order, &
               delays, trim, grid_count, momentum, bootstrap_measure, status)
            if (status /= 0) then
               out%bootstrap_statistic(replication, :) = 0.0_dp
            else
               call regime_count_statistics(bootstrap_measure, &
                  real(size(series) - order, dp), .false., &
                  bootstrap_statistics)
               out%bootstrap_statistic(replication, :) = bootstrap_statistics
            end if
         end do
         do status = 1, 3
            out%p_value(status) = real(count(out%bootstrap_statistic(:, status) > &
               out%statistic(status)), dp)/real(size(bootstrap_series, 2), dp)
         end do
      end if
   end function tsdyn_setar_linearity_test

   pure function tsdyn_tvar_lr_test(series, order, threshold_component, &
      delays, trim, grid_count, momentum, bootstrap_series) result(out)
      !! Compare linear VAR, two-regime TVAR, and three-regime TVAR fits.
      real(dp), intent(in) :: series(:, :) !! Time-by-variable observations.
      integer, intent(in) :: order !! Common VAR order.
      integer, intent(in) :: threshold_component !! Variable driving regime selection.
      integer, intent(in) :: delays(:) !! Transition delays searched under each alternative.
      real(dp), intent(in), optional :: trim !! Minimum threshold-regime fraction.
      integer, intent(in), optional :: grid_count !! Maximum threshold candidates.
      logical, intent(in), optional :: momentum !! Use first-difference transitions.
      real(dp), intent(in), optional :: bootstrap_series(:, :, :) !! Null-bootstrap series by slice.
      type(tsdyn_regime_count_test_t) :: out
      real(dp) :: measure(3), bootstrap_measure(3), bootstrap_statistics(3)
      integer :: replication, status

      call tvar_regime_measures(series, order, threshold_component, delays, trim, &
         grid_count, momentum, measure, status)
      if (status /= 0) then
         out%info = status
         return
      end if
      out%fit_measure = measure
      call regime_count_statistics(measure, real(size(series, 1) - order, dp), &
         .true., out%statistic)
      if (present(bootstrap_series)) then
         if (size(bootstrap_series, 1) /= size(series, 1) .or. &
            size(bootstrap_series, 2) /= size(series, 2)) then
            out%info = 4
            return
         end if
         allocate(out%bootstrap_statistic(size(bootstrap_series, 3), 3))
         do replication = 1, size(bootstrap_series, 3)
            call tvar_regime_measures(bootstrap_series(:, :, replication), order, &
               threshold_component, delays, trim, grid_count, momentum, &
               bootstrap_measure, status)
            if (status /= 0) then
               out%bootstrap_statistic(replication, :) = 0.0_dp
            else
               call regime_count_statistics(bootstrap_measure, &
                  real(size(series, 1) - order, dp), .true., &
                  bootstrap_statistics)
               out%bootstrap_statistic(replication, :) = bootstrap_statistics
            end if
         end do
         do status = 1, 3
            out%p_value(status) = real(count(out%bootstrap_statistic(:, status) > &
               out%statistic(status)), dp)/real(size(bootstrap_series, 3), dp)
         end do
      end if
   end function tsdyn_tvar_lr_test

   pure function tsdyn_bbc_test(series, order, method, trim, &
      minimum_observations, bootstrap_series) result(out)
      !! Test a unit root against a stationary symmetric three-regime SETAR.
      real(dp), intent(in) :: series(:) !! Observed univariate series.
      integer, intent(in) :: order !! Number of lagged differences.
      character(len=*), intent(in), optional :: method !! `LR`, `Wald`, or `LM` statistic.
      real(dp), intent(in), optional :: trim !! Minimum fraction in each regime.
      logical, intent(in), optional :: minimum_observations !! Enforce counts instead of fractions.
      real(dp), intent(in), optional :: bootstrap_series(:, :) !! Null-bootstrap series by column.
      type(tsdyn_bbc_test_t) :: out
      real(dp), allocatable :: centered(:), bootstrap_centered(:)
      real(dp) :: selected_trim
      integer :: replication, status
      logical :: use_counts
      character(len=:), allocatable :: selected_method

      selected_method = 'lr'
      if (present(method)) selected_method = lowercase(method)
      selected_trim = 0.1_dp
      if (present(trim)) selected_trim = trim
      use_counts = .false.
      if (present(minimum_observations)) use_counts = minimum_observations
      if (selected_method /= 'lr' .and. selected_method /= 'wald' .and. &
         selected_method /= 'lm') then
         out%info = 1
         return
      end if
      out%method = selected_method
      centered = series
      if (abs(sum(centered)/real(size(centered), dp)) > 0.01_dp) then
         centered = centered - sum(centered)/real(size(centered), dp)
      end if
      call bbc_statistic_path(centered, order, selected_method, selected_trim, &
         use_counts, out%threshold, out%statistic_path, out%selected, &
         out%statistic, status)
      if (status /= 0) then
         out%info = status
         return
      end if
      select case (selected_method)
      case ('lr')
         out%critical_values = [15.772_dp, 17.898_dp, 22.232_dp]
      case ('wald')
         out%critical_values = [16.181_dp, 18.400_dp, 23.010_dp]
      case ('lm')
         out%critical_values = [15.587_dp, 17.630_dp, 21.756_dp]
      end select
      if (present(bootstrap_series)) then
         if (size(bootstrap_series, 1) /= size(series)) then
            out%info = 4
            return
         end if
         allocate(out%bootstrap_statistic(size(bootstrap_series, 2)))
         do replication = 1, size(bootstrap_series, 2)
            bootstrap_centered = bootstrap_series(:, replication)
            if (abs(sum(bootstrap_centered)/real(size(bootstrap_centered), dp)) > &
               0.01_dp) then
               bootstrap_centered = bootstrap_centered - &
                  sum(bootstrap_centered)/real(size(bootstrap_centered), dp)
            end if
            call bbc_statistic_value(bootstrap_centered, order, selected_method, &
               selected_trim, use_counts, out%bootstrap_statistic(replication), &
               status)
            if (status /= 0) out%bootstrap_statistic(replication) = 0.0_dp
         end do
         out%p_value = real(count(out%bootstrap_statistic > out%statistic), dp)/ &
            real(size(out%bootstrap_statistic), dp)
      end if
   end function tsdyn_bbc_test

   pure function tsdyn_kapshin_test(series, order, deterministic, grid_scale, &
      grid_exponent, points, minimum_middle, bootstrap_series) result(out)
      !! Test a unit root using Kapetanios-Shin outer-regime Wald statistics.
      real(dp), intent(in) :: series(:) !! Observed univariate series.
      integer, intent(in) :: order !! Number of lagged differences.
      character(len=*), intent(in), optional :: deterministic !! `none`, `const`, `trend`, or `both`.
      real(dp), intent(in), optional :: grid_scale !! Central-grid scale parameter.
      real(dp), intent(in), optional :: grid_exponent !! Central-grid rate exponent.
      integer, intent(in), optional :: points !! Number of marginal threshold points.
      integer, intent(in), optional :: minimum_middle !! Minimum middle-regime observations.
      real(dp), intent(in), optional :: bootstrap_series(:, :) !! Null-bootstrap series by column.
      type(tsdyn_kapshin_test_t) :: out
      real(dp), allocatable :: adjusted(:), bootstrap_adjusted(:)
      real(dp) :: selected_scale, selected_exponent, bootstrap_statistics(3)
      integer :: selected_points, selected_middle, replication, status, statistic_index
      character(len=:), allocatable :: selected_deterministic

      selected_deterministic = 'none'
      if (present(deterministic)) then
         selected_deterministic = lowercase(trim(deterministic))
      end if
      selected_scale = 3.0_dp
      if (present(grid_scale)) selected_scale = grid_scale
      selected_exponent = 0.5_dp
      if (present(grid_exponent)) selected_exponent = grid_exponent
      selected_points = 0
      if (present(points)) selected_points = points
      selected_middle = 10
      if (present(minimum_middle)) selected_middle = minimum_middle
      if (selected_deterministic /= 'none' .and. &
         selected_deterministic /= 'const' .and. &
         selected_deterministic /= 'trend' .and. &
         selected_deterministic /= 'both') then
         out%info = 1
         return
      end if
      out%deterministic = selected_deterministic
      call remove_deterministic(series, selected_deterministic, adjusted, status)
      if (status /= 0) then
         out%info = 2
         return
      end if
      call kapshin_statistic_path(adjusted, order, selected_scale, &
         selected_exponent, selected_points, selected_middle, out%thresholds, &
         out%wald_path, out%selected, out%statistic, status)
      if (status /= 0) then
         out%info = status
         return
      end if
      if (present(bootstrap_series)) then
         if (size(bootstrap_series, 1) /= size(series)) then
            out%info = 4
            return
         end if
         allocate(out%bootstrap_statistic(size(bootstrap_series, 2), 3))
         do replication = 1, size(bootstrap_series, 2)
            call remove_deterministic(bootstrap_series(:, replication), &
               selected_deterministic, bootstrap_adjusted, status)
            if (status == 0) then
               call kapshin_statistic_value(bootstrap_adjusted, order, &
                  selected_scale, selected_exponent, selected_points, &
                  selected_middle, bootstrap_statistics, status)
               out%bootstrap_statistic(replication, :) = bootstrap_statistics
            end if
            if (status /= 0) out%bootstrap_statistic(replication, :) = 0.0_dp
         end do
         do statistic_index = 1, 3
            out%p_value(statistic_index) = real(count( &
               out%bootstrap_statistic(:, statistic_index) > &
               out%statistic(statistic_index)), dp)/ &
               real(size(out%bootstrap_statistic, 1), dp)
         end do
      end if
   end function tsdyn_kapshin_test

   pure function tsdyn_hansen_seo_test(data, lag, beta, trim, grid_count, &
      wild_multipliers) result(out)
      !! Test linear against threshold cointegration using a supremum score test.
      real(dp), intent(in) :: data(:, :) !! Time-by-two level observations.
      integer, intent(in) :: lag !! Number of lagged differences.
      real(dp), intent(in) :: beta !! Fixed cointegrating coefficient.
      real(dp), intent(in), optional :: trim !! Minimum threshold-regime fraction.
      integer, intent(in), optional :: grid_count !! Maximum threshold candidates.
      real(dp), intent(in), optional :: wild_multipliers(:, :) !! Row-by-bootstrap multipliers.
      type(tsdyn_threshold_test_t) :: out
      real(dp), allocatable :: response(:, :), common(:, :), ect(:), design(:, :)
      real(dp), allocatable :: coefficient(:, :), fitted(:, :), residual(:, :)
      real(dp), allocatable :: ordered(:), grid(:), bootstrap_response(:, :)
      real(dp), allocatable :: ordered_bootstrap(:)
      real(dp) :: selected_trim, rss
      integer :: selected_grid, candidate, replicate, status

      selected_trim = 0.05_dp
      if (present(trim)) selected_trim = trim
      selected_grid = 100
      if (present(grid_count)) selected_grid = grid_count
      out%beta = beta
      if (size(data, 2) /= 2 .or. lag < 0 .or. &
         selected_trim <= 0.0_dp .or. selected_trim >= 0.5_dp) then
         out%info = 1
         return
      end if
      call tvecm_regression_data(data, lag, [1.0_dp, -beta], .true., .false., &
         response, common, ect)
      allocate(design(size(ect), 1 + size(common, 2)))
      design(:, 1) = ect
      design(:, 2:) = common
      call tvecm_multivariate_regression(design, response, coefficient, fitted, &
         residual, rss, status)
      if (status /= 0) then
         out%info = 2
         return
      end if
      ordered = sorted(ect)
      call threshold_candidate_grid(ordered, selected_trim, selected_grid, grid)
      allocate(out%thresholds(size(grid), 1), out%statistic_path(size(grid)))
      out%thresholds(:, 1) = grid
      do candidate = 1, size(grid)
         out%statistic_path(candidate) = hansen_seo_score(response, design, ect, &
            grid(candidate), selected_trim)
      end do
      out%selected = maxloc(out%statistic_path, dim=1)
      out%statistic = out%statistic_path(out%selected)
      if (present(wild_multipliers)) then
         if (size(wild_multipliers, 1) /= size(response, 1)) then
            out%info = 3
            return
         end if
         allocate(out%bootstrap_statistic(size(wild_multipliers, 2)))
         allocate(bootstrap_response(size(response, 1), 2))
         do replicate = 1, size(wild_multipliers, 2)
            bootstrap_response = residual*spread(wild_multipliers(:, replicate), 2, 2)
            out%bootstrap_statistic(replicate) = 0.0_dp
            do candidate = 1, size(grid)
               out%bootstrap_statistic(replicate) = max( &
                  out%bootstrap_statistic(replicate), &
                  hansen_seo_score(bootstrap_response, design, ect, &
                  grid(candidate), selected_trim))
            end do
         end do
         out%p_value = real(count(out%bootstrap_statistic > out%statistic), dp)/ &
            real(size(out%bootstrap_statistic), dp)
         ordered_bootstrap = sorted(out%bootstrap_statistic)
         out%critical_values = [quantile(ordered_bootstrap, 0.90_dp), &
            quantile(ordered_bootstrap, 0.95_dp), &
            quantile(ordered_bootstrap, 0.99_dp)]
      end if
   end function tsdyn_hansen_seo_test

   pure function tsdyn_seo_test(data, lag, beta, trim, grid_count, &
      residual_indices) result(out)
      !! Test no cointegration against outer-regime threshold cointegration.
      real(dp), intent(in) :: data(:, :) !! Time-by-two level observations.
      integer, intent(in) :: lag !! Number of lagged differences.
      real(dp), intent(in) :: beta !! Fixed cointegrating coefficient.
      real(dp), intent(in), optional :: trim !! Minimum fraction in each regime.
      integer, intent(in), optional :: grid_count !! Maximum one-dimensional grid size.
      integer, intent(in), optional :: residual_indices(:, :) !! Row-by-bootstrap residual indices.
      type(tsdyn_threshold_test_t) :: out
      real(dp), allocatable :: response(:, :), common(:, :), ect(:)
      real(dp), allocatable :: coefficient(:, :), fitted(:, :), residual(:, :)
      real(dp), allocatable :: bootstrap_data(:, :), bootstrap_response(:, :)
      real(dp), allocatable :: bootstrap_common(:, :), bootstrap_ect(:)
      real(dp), allocatable :: ordered_bootstrap(:), ignored_thresholds(:, :)
      real(dp), allocatable :: ignored_path(:)
      real(dp) :: selected_trim, rss
      integer :: selected_grid, replicate, ignored_selected, status

      selected_trim = 0.1_dp
      if (present(trim)) selected_trim = trim
      selected_grid = 30
      if (present(grid_count)) selected_grid = grid_count
      out%beta = beta
      if (size(data, 2) /= 2 .or. lag < 0 .or. &
         selected_trim <= 0.0_dp .or. selected_trim >= 1.0_dp/3.0_dp) then
         out%info = 1
         return
      end if
      call tvecm_regression_data(data, lag, [1.0_dp, -beta], .true., .false., &
         response, common, ect)
      call seo_supremum_wald(response, common, ect, selected_trim, selected_grid, &
         out%thresholds, out%statistic_path, out%selected, out%statistic, status)
      if (status /= 0) then
         out%info = 2
         return
      end if
      call tvecm_multivariate_regression(common, response, coefficient, fitted, &
         residual, rss, status)
      if (status /= 0) then
         out%info = 3
         return
      end if
      if (present(residual_indices)) then
         if (size(residual_indices, 1) /= size(residual, 1) .or. &
            minval(residual_indices) < 1 .or. &
            maxval(residual_indices) > size(residual, 1)) then
            out%info = 4
            return
         end if
         allocate(out%bootstrap_statistic(size(residual_indices, 2)))
         do replicate = 1, size(residual_indices, 2)
            call simulate_difference_null(data, lag, coefficient, residual, &
               residual_indices(:, replicate), bootstrap_data)
            call tvecm_regression_data(bootstrap_data, lag, [1.0_dp, -beta], &
               .true., .false., &
               bootstrap_response, bootstrap_common, bootstrap_ect)
            call seo_supremum_wald(bootstrap_response, bootstrap_common, &
               bootstrap_ect, selected_trim, selected_grid, ignored_thresholds, &
               ignored_path, ignored_selected, out%bootstrap_statistic(replicate), &
               status)
            if (status /= 0) out%bootstrap_statistic(replicate) = 0.0_dp
         end do
         out%p_value = real(count(out%bootstrap_statistic > out%statistic), dp)/ &
            real(size(out%bootstrap_statistic), dp)
         ordered_bootstrap = sorted(out%bootstrap_statistic)
         out%critical_values = [quantile(ordered_bootstrap, 0.90_dp), &
            quantile(ordered_bootstrap, 0.95_dp), &
            quantile(ordered_bootstrap, 0.975_dp), &
            quantile(ordered_bootstrap, 0.99_dp)]
      end if
   end function tsdyn_seo_test

   pure function tsdyn_setar_regime_irf(model, regime, horizon, cumulative) &
      result(out)
      !! Compute the linear impulse response conditional on one SETAR regime.
      type(nts_tar_model_t), intent(in) :: model !! Fitted SETAR model.
      integer, intent(in) :: regime !! Regime held fixed for propagation.
      integer, intent(in) :: horizon !! Largest response horizon, including zero.
      logical, intent(in), optional :: cumulative !! Accumulate responses over horizons.
      type(tsdyn_irf_t) :: out
      real(dp), allocatable :: response(:)
      integer :: step, lag_index
      logical :: accumulate

      accumulate = .false.
      if (present(cumulative)) accumulate = cumulative
      if (model%info /= 0 .or. .not. allocated(model%ar_order) .or. &
         .not. allocated(model%coefficients) .or. horizon < 0 .or. &
         regime < 1 .or. regime > size(model%ar_order)) then
         out%info = 1
         return
      end if
      allocate(response(horizon + 1), out%response(horizon + 1, 1, 1))
      response = 0.0_dp
      response(1) = 1.0_dp
      do step = 1, horizon
         do lag_index = 1, min(model%ar_order(regime), step)
            response(step + 1) = response(step + 1) + &
               model%coefficients(regime, lag_index + 1)* &
               response(step - lag_index + 1)
         end do
      end do
      if (accumulate) call cumulative_irf(response)
      out%response(:, 1, 1) = response
      out%regime = regime
      out%cumulative = accumulate
   end function tsdyn_setar_regime_irf

   pure function tsdyn_tvar_regime_irf(model, regime, horizon, cumulative, &
      orthogonalized) result(out)
      !! Compute fixed-regime TVAR responses to unit or orthogonalized shocks.
      type(nts_mtar_model_t), intent(in) :: model !! Fitted TVAR model.
      integer, intent(in) :: regime !! Regime held fixed for propagation.
      integer, intent(in) :: horizon !! Largest response horizon, including zero.
      logical, intent(in), optional :: cumulative !! Accumulate responses over horizons.
      logical, intent(in), optional :: orthogonalized !! Use Cholesky impact shocks.
      type(tsdyn_irf_t) :: out
      real(dp), allocatable :: impact(:, :), response(:, :, :)
      integer :: status
      logical :: accumulate, orthogonal

      accumulate = .false.
      if (present(cumulative)) accumulate = cumulative
      orthogonal = .false.
      if (present(orthogonalized)) orthogonal = orthogonalized
      if (model%info /= 0 .or. .not. allocated(model%ar_order) .or. &
         horizon < 0 .or. regime < 1 .or. regime > size(model%ar_order)) then
         out%info = 1
         return
      end if
      call tvar_impact_matrix(model, orthogonal, impact, status)
      if (status /= 0) then
         out%info = 2
         return
      end if
      call tvar_fixed_regime_response(model, regime, horizon, impact, response, status)
      if (status /= 0) then
         out%info = 3
         return
      end if
      if (accumulate) call cumulative_irf_array(response)
      out%response = response
      out%regime = regime
      out%cumulative = accumulate
      out%orthogonalized = orthogonal
   end function tsdyn_tvar_regime_irf

   pure function tsdyn_tvecm_regime_irf(model, regime, horizon, cumulative, &
      orthogonalized) result(out)
      !! Compute fixed-regime TVECM level responses to supplied shock conventions.
      type(tsdyn_tvecm_model_t), intent(in) :: model !! Fitted TVECM model.
      integer, intent(in) :: regime !! Regime held fixed for propagation.
      integer, intent(in) :: horizon !! Largest response horizon, including zero.
      logical, intent(in), optional :: cumulative !! Accumulate level responses.
      logical, intent(in), optional :: orthogonalized !! Use Cholesky impact shocks.
      type(tsdyn_irf_t) :: out
      real(dp), allocatable :: impact(:, :), response(:, :, :)
      integer :: status
      logical :: accumulate, orthogonal

      accumulate = .false.
      if (present(cumulative)) accumulate = cumulative
      orthogonal = .false.
      if (present(orthogonalized)) orthogonal = orthogonalized
      if (model%info /= 0 .or. .not. allocated(model%thresholds) .or. &
         horizon < 0 .or. regime < 1 .or. &
         regime > size(model%thresholds) + 1) then
         out%info = 1
         return
      end if
      call tvecm_impact_matrix(model, orthogonal, impact, status)
      if (status /= 0) then
         out%info = 2
         return
      end if
      call tvecm_fixed_regime_response(model, regime, horizon, impact, &
         response, status)
      if (status /= 0) then
         out%info = 3
         return
      end if
      if (accumulate) call cumulative_irf_array(response)
      out%response = response
      out%regime = regime
      out%cumulative = accumulate
      out%orthogonalized = orthogonal
   end function tsdyn_tvecm_regime_irf

   pure function tsdyn_setar_girf_from_innovations(model, history, innovations, &
      shocks) result(out)
      !! Compute nonlinear SETAR GIRFs for supplied histories and innovations.
      type(nts_tar_model_t), intent(in) :: model !! Fitted SETAR model.
      real(dp), intent(in) :: history(:) !! Presample level history.
      real(dp), intent(in) :: innovations(:, :) !! Horizon-by-replication innovations.
      real(dp), intent(in) :: shocks(:) !! Alternative first-period shock sizes.
      type(tsdyn_girf_t) :: out
      real(dp), allocatable :: baseline(:), shocked(:)
      integer :: replication, shock, status

      if (model%info /= 0 .or. size(innovations, 1) < 1 .or. &
         size(innovations, 2) < 1 .or. size(shocks) < 1) then
         out%info = 1
         return
      end if
      allocate(out%response(size(innovations, 1), 1, size(shocks)))
      out%response = 0.0_dp
      do replication = 1, size(innovations, 2)
         call setar_impulse_path(model, history, innovations(:, replication), &
            0.0_dp, .false., baseline, status)
         if (status /= 0) then
            out%info = 2
            return
         end if
         do shock = 1, size(shocks)
            call setar_impulse_path(model, history, innovations(:, replication), &
               shocks(shock), .true., shocked, status)
            if (status /= 0) then
               out%info = 2
               return
            end if
            out%response(:, 1, shock) = out%response(:, 1, shock) + &
               shocked - baseline
         end do
      end do
      out%simulations = size(innovations, 2)
      out%response = out%response/real(out%simulations, dp)
   end function tsdyn_setar_girf_from_innovations

   pure function tsdyn_setar_irf_bootstrap(model, bootstrap_models, regime, &
      horizon, level, cumulative) result(out)
      !! Form SETAR regime-IRF confidence intervals from supplied bootstrap refits.
      type(nts_tar_model_t), intent(in) :: model !! Original fitted SETAR model.
      type(nts_tar_model_t), intent(in) :: bootstrap_models(:) !! Bootstrap refitted models.
      integer, intent(in) :: regime !! Regime held fixed for propagation.
      integer, intent(in) :: horizon !! Largest response horizon, including zero.
      real(dp), intent(in), optional :: level !! Central confidence level.
      logical, intent(in), optional :: cumulative !! Accumulate responses over horizons.
      type(tsdyn_irf_t) :: out
      type(tsdyn_irf_t) :: candidate
      real(dp), allocatable :: draws(:, :, :, :)
      real(dp) :: selected_level
      integer :: bootstrap

      selected_level = 0.95_dp
      if (present(level)) selected_level = level
      out = tsdyn_setar_regime_irf(model, regime, horizon, cumulative)
      if (out%info /= 0 .or. size(bootstrap_models) < 1 .or. &
         selected_level <= 0.0_dp .or. selected_level >= 1.0_dp) then
         out%info = 1
         return
      end if
      allocate(draws(horizon + 1, 1, 1, size(bootstrap_models)))
      do bootstrap = 1, size(bootstrap_models)
         candidate = tsdyn_setar_regime_irf(bootstrap_models(bootstrap), &
            regime, horizon, cumulative)
         if (candidate%info /= 0) then
            out%info = 2
            return
         end if
         draws(:, :, :, bootstrap) = candidate%response
      end do
      call irf_bootstrap_intervals(draws, selected_level, out%lower, out%upper)
      out%level = selected_level
      out%bootstrap_count = size(bootstrap_models)
   end function tsdyn_setar_irf_bootstrap

   pure function tsdyn_tvar_irf_bootstrap(model, bootstrap_models, regime, &
      horizon, level, cumulative, orthogonalized) result(out)
      !! Form TVAR regime-IRF confidence intervals from supplied bootstrap refits.
      type(nts_mtar_model_t), intent(in) :: model !! Original fitted TVAR model.
      type(nts_mtar_model_t), intent(in) :: bootstrap_models(:) !! Bootstrap refitted models.
      integer, intent(in) :: regime !! Regime held fixed for propagation.
      integer, intent(in) :: horizon !! Largest response horizon, including zero.
      real(dp), intent(in), optional :: level !! Central confidence level.
      logical, intent(in), optional :: cumulative !! Accumulate responses over horizons.
      logical, intent(in), optional :: orthogonalized !! Use Cholesky impact shocks.
      type(tsdyn_irf_t) :: out
      type(tsdyn_irf_t) :: candidate
      real(dp), allocatable :: draws(:, :, :, :)
      real(dp) :: selected_level
      integer :: bootstrap

      selected_level = 0.95_dp
      if (present(level)) selected_level = level
      out = tsdyn_tvar_regime_irf(model, regime, horizon, cumulative, orthogonalized)
      if (out%info /= 0 .or. size(bootstrap_models) < 1 .or. &
         selected_level <= 0.0_dp .or. selected_level >= 1.0_dp) then
         out%info = 1
         return
      end if
      allocate(draws, mold=spread(out%response, 4, size(bootstrap_models)))
      do bootstrap = 1, size(bootstrap_models)
         candidate = tsdyn_tvar_regime_irf(bootstrap_models(bootstrap), regime, &
            horizon, cumulative, orthogonalized)
         if (candidate%info /= 0) then
            out%info = 2
            return
         end if
         draws(:, :, :, bootstrap) = candidate%response
      end do
      call irf_bootstrap_intervals(draws, selected_level, out%lower, out%upper)
      out%level = selected_level
      out%bootstrap_count = size(bootstrap_models)
   end function tsdyn_tvar_irf_bootstrap

   pure function tsdyn_tvecm_irf_bootstrap(model, bootstrap_models, regime, &
      horizon, level, cumulative, orthogonalized) result(out)
      !! Form TVECM regime-IRF confidence intervals from supplied bootstrap refits.
      type(tsdyn_tvecm_model_t), intent(in) :: model !! Original fitted TVECM model.
      type(tsdyn_tvecm_model_t), intent(in) :: bootstrap_models(:) !! Bootstrap refitted models.
      integer, intent(in) :: regime !! Regime held fixed for propagation.
      integer, intent(in) :: horizon !! Largest response horizon, including zero.
      real(dp), intent(in), optional :: level !! Central confidence level.
      logical, intent(in), optional :: cumulative !! Accumulate responses over horizons.
      logical, intent(in), optional :: orthogonalized !! Use Cholesky impact shocks.
      type(tsdyn_irf_t) :: out
      type(tsdyn_irf_t) :: candidate
      real(dp), allocatable :: draws(:, :, :, :)
      real(dp) :: selected_level
      integer :: bootstrap

      selected_level = 0.95_dp
      if (present(level)) selected_level = level
      out = tsdyn_tvecm_regime_irf(model, regime, horizon, cumulative, &
         orthogonalized)
      if (out%info /= 0 .or. size(bootstrap_models) < 1 .or. &
         selected_level <= 0.0_dp .or. selected_level >= 1.0_dp) then
         out%info = 1
         return
      end if
      allocate(draws, mold=spread(out%response, 4, size(bootstrap_models)))
      do bootstrap = 1, size(bootstrap_models)
         candidate = tsdyn_tvecm_regime_irf(bootstrap_models(bootstrap), &
            regime, horizon, cumulative, orthogonalized)
         if (candidate%info /= 0) then
            out%info = 2
            return
         end if
         draws(:, :, :, bootstrap) = candidate%response
      end do
      call irf_bootstrap_intervals(draws, selected_level, out%lower, out%upper)
      out%level = selected_level
      out%bootstrap_count = size(bootstrap_models)
   end function tsdyn_tvecm_irf_bootstrap

   pure function tsdyn_tvar_girf_from_innovations(model, history, innovations, &
      shocks) result(out)
      !! Compute paired-simulation generalized impulse responses for a TVAR.
      type(nts_mtar_model_t), intent(in) :: model !! Fitted threshold VAR model.
      real(dp), intent(in) :: history(:, :) !! Presample level history.
      real(dp), intent(in) :: innovations(:, :, :) !! Horizon-by-variable-by-replication innovations.
      real(dp), intent(in) :: shocks(:, :) !! Variable-by-shock impact vectors.
      type(tsdyn_girf_t) :: out
      real(dp), allocatable :: baseline(:, :), shocked(:, :)
      integer :: replication, shock, horizon, variables, status

      horizon = size(innovations, 1)
      variables = size(innovations, 2)
      if (model%info /= 0 .or. size(history, 2) /= variables .or. &
         size(shocks, 1) /= variables .or. size(innovations, 3) < 1) then
         out%info = 1
         return
      end if
      allocate(out%response(horizon, variables, size(shocks, 2)))
      out%response = 0.0_dp
      do replication = 1, size(innovations, 3)
         call tvar_impulse_path(model, history, innovations(:, :, replication), &
            [real(dp) ::], baseline, status)
         if (status /= 0) then
            out%info = 2
            return
         end if
         do shock = 1, size(shocks, 2)
            call tvar_impulse_path(model, history, innovations(:, :, replication), &
               shocks(:, shock), shocked, status)
            if (status /= 0) then
               out%info = 2
               return
            end if
            out%response(:, :, shock) = out%response(:, :, shock) + &
               shocked - baseline
         end do
      end do
      out%simulations = size(innovations, 3)
      out%response = out%response/real(out%simulations, dp)
   end function tsdyn_tvar_girf_from_innovations

   pure function tsdyn_tvecm_girf_from_innovations(model, history, &
      innovations, shocks) result(out)
      !! Compute paired-simulation generalized impulse responses for a TVECM.
      type(tsdyn_tvecm_model_t), intent(in) :: model !! Fitted threshold VECM.
      real(dp), intent(in) :: history(:, :) !! Presample level history.
      real(dp), intent(in) :: innovations(:, :, :) !! Horizon-by-variable-by-replication innovations.
      real(dp), intent(in) :: shocks(:, :) !! Variable-by-shock impact vectors.
      type(tsdyn_girf_t) :: out
      real(dp), allocatable :: baseline(:, :), shocked(:, :)
      integer :: replication, shock, variables, status

      variables = size(innovations, 2)
      if (model%info /= 0 .or. variables < 2 .or. &
         size(history, 2) /= variables .or. &
         size(shocks, 1) /= variables .or. &
         size(innovations, 3) < 1) then
         out%info = 1
         return
      end if
      allocate(out%response(size(innovations, 1), variables, size(shocks, 2)))
      out%response = 0.0_dp
      do replication = 1, size(innovations, 3)
         call tvecm_impulse_path(model, history, innovations(:, :, replication), &
            [real(dp) ::], baseline, status)
         if (status /= 0) then
            out%info = 2
            return
         end if
         do shock = 1, size(shocks, 2)
            call tvecm_impulse_path(model, history, &
               innovations(:, :, replication), shocks(:, shock), shocked, status)
            if (status /= 0) then
               out%info = 2
               return
            end if
            out%response(:, :, shock) = out%response(:, :, shock) + &
               shocked - baseline
         end do
      end do
      out%simulations = size(innovations, 3)
      out%response = out%response/real(out%simulations, dp)
   end function tsdyn_tvecm_girf_from_innovations

   pure function tsdyn_nonlinear_fevd(girf) result(out)
      !! Normalize cumulative squared generalized responses into variance shares.
      type(tsdyn_girf_t), intent(in) :: girf !! Nonlinear generalized responses.
      type(tsdyn_fevd_t) :: out
      real(dp) :: denominator
      integer :: horizon, response

      if (girf%info /= 0 .or. .not. allocated(girf%response)) then
         out%info = 1
         return
      end if
      allocate(out%decomposition, mold=girf%response)
      out%decomposition = 0.0_dp
      do horizon = 1, size(girf%response, 1)
         do response = 1, size(girf%response, 2)
            out%decomposition(horizon, response, :) = sum( &
               girf%response(:horizon, response, :)**2, dim=1)
            denominator = sum(out%decomposition(horizon, response, :))
            if (denominator > tiny(1.0_dp)) then
               out%decomposition(horizon, response, :) = &
                  out%decomposition(horizon, response, :)/denominator
            end if
         end do
      end do
   end function tsdyn_nonlinear_fevd

   pure real(dp) function hansen_seo_score(response, restricted_design, ect, &
      threshold, trim) result(statistic)
      !! Evaluate the heteroskedasticity-robust threshold score at one split.
      real(dp), intent(in) :: response(:, :) !! Aligned dependent changes.
      real(dp), intent(in) :: restricted_design(:, :) !! Linear VECM design matrix.
      real(dp), intent(in) :: ect(:) !! Lagged error-correction values.
      real(dp), intent(in) :: threshold !! Candidate threshold.
      real(dp), intent(in) :: trim !! Minimum regime fraction.
      real(dp), allocatable :: coefficient(:, :), fitted(:, :), residual(:, :)
      real(dp), allocatable :: cross_product(:, :), inverse(:, :), interacted(:, :)
      real(dp), allocatable :: residualized(:, :), score(:), meat(:, :), vector(:)
      real(dp) :: rss
      integer :: rows, columns, variables, row, equation, column, index, status

      statistic = 0.0_dp
      rows = size(response, 1)
      columns = size(restricted_design, 2)
      variables = size(response, 2)
      if (count(ect <= threshold) <= int(trim*real(rows, dp)) .or. &
         count(ect > threshold) <= int(trim*real(rows, dp))) return
      call tvecm_multivariate_regression(restricted_design, response, coefficient, &
         fitted, residual, rss, status)
      if (status /= 0) return
      cross_product = matmul(transpose(restricted_design), restricted_design)
      call invert_matrix(cross_product, inverse, status)
      if (status /= 0) return
      allocate(interacted(rows, columns))
      interacted = restricted_design*spread(merge(1.0_dp, 0.0_dp, &
         ect <= threshold), 2, columns)
      residualized = interacted - matmul(restricted_design, matmul(inverse, &
         matmul(transpose(restricted_design), interacted)))
      allocate(score(columns*variables), meat(columns*variables, columns*variables))
      allocate(vector(columns*variables))
      score = 0.0_dp
      meat = 0.0_dp
      do row = 1, rows
         index = 0
         do equation = 1, variables
            do column = 1, columns
               index = index + 1
               vector(index) = residualized(row, column)*residual(row, equation)
            end do
         end do
         score = score + vector
         meat = meat + spread(vector, 2, size(vector))* &
            spread(vector, 1, size(vector))
      end do
      do index = 1, size(score)
         meat(index, index) = meat(index, index) + 1.0e-10_dp
      end do
      call invert_matrix(meat, inverse, status)
      if (status /= 0) return
      statistic = max(0.0_dp, dot_product(score, matmul(inverse, score)))
   end function hansen_seo_score

   pure subroutine seo_supremum_wald(response, common, ect, trim, grid_count, &
      thresholds, path, selected, statistic, info)
      !! Search the two-threshold outer-regime Seo supremum-Wald statistic.
      real(dp), intent(in) :: response(:, :) !! Aligned dependent changes.
      real(dp), intent(in) :: common(:, :) !! Deterministic and lagged-change regressors.
      real(dp), intent(in) :: ect(:) !! Lagged error-correction values.
      real(dp), intent(in) :: trim !! Minimum fraction in every regime.
      integer, intent(in) :: grid_count !! Maximum marginal threshold-grid size.
      real(dp), allocatable, intent(out) :: thresholds(:, :) !! Evaluated threshold pairs.
      real(dp), allocatable, intent(out) :: path(:) !! Wald statistic for each pair.
      integer, intent(out) :: selected !! Index of the supremum pair.
      real(dp), intent(out) :: statistic !! Supremum Wald statistic.
      integer, intent(out) :: info !! Zero on success.
      real(dp), allocatable :: ordered(:), grid(:), threshold_work(:, :), path_work(:)
      integer :: first, second, pair_count, minimum

      info = 1
      selected = 0
      statistic = 0.0_dp
      ordered = sorted(ect)
      call threshold_candidate_grid(ordered, trim, grid_count, grid)
      allocate(threshold_work(max(1, size(grid)*(size(grid) - 1)/2), 2))
      allocate(path_work(size(threshold_work, 1)))
      pair_count = 0
      minimum = ceiling(trim*real(size(ect), dp))
      do first = 1, size(grid) - 1
         do second = first + 1, size(grid)
            if (count(ect < grid(first)) < minimum .or. &
               count(ect >= grid(first) .and. ect <= grid(second)) < minimum .or. &
               count(ect > grid(second)) < minimum) cycle
            pair_count = pair_count + 1
            threshold_work(pair_count, :) = [grid(first), grid(second)]
            path_work(pair_count) = seo_wald_at(response, common, ect, &
               grid(first), grid(second))
         end do
      end do
      if (pair_count == 0) then
         allocate(thresholds(0, 2), path(0))
         return
      end if
      thresholds = threshold_work(:pair_count, :)
      path = path_work(:pair_count)
      selected = maxloc(path, dim=1)
      statistic = path(selected)
      info = 0
   end subroutine seo_supremum_wald

   pure real(dp) function seo_wald_at(response, common, ect, lower, upper) &
      result(statistic)
      !! Evaluate the Seo outer-regime Wald statistic for one threshold pair.
      real(dp), intent(in) :: response(:, :) !! Aligned dependent changes.
      real(dp), intent(in) :: common(:, :) !! Regressors under no cointegration.
      real(dp), intent(in) :: ect(:) !! Lagged error-correction values.
      real(dp), intent(in) :: lower !! Lower threshold.
      real(dp), intent(in) :: upper !! Upper threshold.
      real(dp), allocatable :: z(:, :), common_cross(:, :), common_inverse(:, :)
      real(dp), allocatable :: residualized(:, :), a(:, :), a_inverse(:, :)
      real(dp), allocatable :: alpha(:, :), design(:, :), coefficient(:, :)
      real(dp), allocatable :: fitted(:, :), residual(:, :), sigma(:, :), sigma_inverse(:, :)
      real(dp) :: rss
      integer :: first, second, equation, other_equation, status

      statistic = 0.0_dp
      allocate(z(size(ect), 2))
      z(:, 1) = merge(ect, 0.0_dp, ect < lower)
      z(:, 2) = merge(ect, 0.0_dp, ect > upper)
      common_cross = matmul(transpose(common), common)
      call invert_matrix(common_cross, common_inverse, status)
      if (status /= 0) return
      residualized = z - matmul(common, matmul(common_inverse, &
         matmul(transpose(common), z)))
      a = matmul(transpose(z), residualized)
      call invert_matrix(a, a_inverse, status)
      if (status /= 0) return
      alpha = matmul(a_inverse, matmul(transpose(residualized), response))
      allocate(design(size(ect), 2 + size(common, 2)))
      design(:, :2) = z
      design(:, 3:) = common
      call tvecm_multivariate_regression(design, response, coefficient, fitted, &
         residual, rss, status)
      if (status /= 0) return
      sigma = matmul(transpose(residual), residual)/real(size(residual, 1), dp)
      call invert_matrix(sigma, sigma_inverse, status)
      if (status /= 0) return
      do first = 1, 2
         do second = 1, 2
            do equation = 1, size(response, 2)
               do other_equation = 1, size(response, 2)
                  statistic = statistic + alpha(first, equation)*a(first, second)* &
                     sigma_inverse(equation, other_equation)* &
                     alpha(second, other_equation)
               end do
            end do
         end do
      end do
      statistic = max(0.0_dp, statistic)
   end function seo_wald_at

   pure subroutine simulate_difference_null(data, lag, coefficient, residual, &
      indices, simulated)
      !! Simulate integrated levels under a VAR in differences with resampled errors.
      real(dp), intent(in) :: data(:, :) !! Original level observations.
      integer, intent(in) :: lag !! Difference autoregressive order.
      real(dp), intent(in) :: coefficient(:, :) !! Equation-by-regressor null coefficients.
      real(dp), intent(in) :: residual(:, :) !! Centered null residuals.
      integer, intent(in) :: indices(:) !! Resampled residual row indices.
      real(dp), allocatable, intent(out) :: simulated(:, :) !! Bootstrap level series.
      real(dp), allocatable :: regressor(:), change(:)
      integer :: time, row, lag_index, column

      allocate(simulated, mold=data)
      simulated = 0.0_dp
      simulated(:lag + 1, :) = data(:lag + 1, :)
      allocate(regressor(1 + 2*lag), change(2))
      do time = lag + 2, size(data, 1)
         row = time - lag - 1
         regressor(1) = 1.0_dp
         column = 1
         do lag_index = 1, lag
            regressor(column + 1:column + 2) = &
               simulated(time - lag_index, :) - &
               simulated(time - lag_index - 1, :)
            column = column + 2
         end do
         change = matmul(coefficient, regressor) + residual(indices(row), :)
         simulated(time, :) = simulated(time - 1, :) + change
      end do
   end subroutine simulate_difference_null

   pure subroutine cumulative_irf(response)
      !! Replace a scalar impulse response by its running sum.
      real(dp), intent(inout) :: response(:) !! Response sequence to accumulate.
      integer :: step

      do step = 2, size(response)
         response(step) = response(step) + response(step - 1)
      end do
   end subroutine cumulative_irf

   pure subroutine cumulative_irf_array(response)
      !! Replace multivariate impulse responses by horizon-wise running sums.
      real(dp), intent(inout) :: response(:, :, :) !! Horizon-response-shock array.
      integer :: step

      do step = 2, size(response, 1)
         response(step, :, :) = response(step, :, :) + response(step - 1, :, :)
      end do
   end subroutine cumulative_irf_array

   pure subroutine irf_bootstrap_intervals(draws, level, lower, upper)
      !! Reduce bootstrap response draws to equal-tailed pointwise intervals.
      real(dp), intent(in) :: draws(:, :, :, :) !! Horizon-response-shock-bootstrap draws.
      real(dp), intent(in) :: level !! Central confidence level.
      real(dp), allocatable, intent(out) :: lower(:, :, :) !! Pointwise lower endpoints.
      real(dp), allocatable, intent(out) :: upper(:, :, :) !! Pointwise upper endpoints.
      real(dp), allocatable :: ordered(:)
      real(dp) :: tail
      integer :: step, response_index, shock

      allocate(lower(size(draws, 1), size(draws, 2), size(draws, 3)))
      allocate(upper(size(draws, 1), size(draws, 2), size(draws, 3)))
      tail = 0.5_dp*(1.0_dp - level)
      do shock = 1, size(draws, 3)
         do response_index = 1, size(draws, 2)
            do step = 1, size(draws, 1)
               ordered = sorted(draws(step, response_index, shock, :))
               lower(step, response_index, shock) = quantile(ordered, tail)
               upper(step, response_index, shock) = quantile(ordered, 1.0_dp - tail)
            end do
         end do
      end do
   end subroutine irf_bootstrap_intervals

   pure subroutine setar_impulse_path(model, history, innovations, shock, &
      apply_shock, path, info)
      !! Simulate one SETAR path, optionally perturbing its first innovation.
      type(nts_tar_model_t), intent(in) :: model !! Fitted SETAR model.
      real(dp), intent(in) :: history(:) !! Presample level history.
      real(dp), intent(in) :: innovations(:) !! Future innovation sequence.
      real(dp), intent(in) :: shock !! First-period additive shock.
      logical, intent(in) :: apply_shock !! Whether to add the shock.
      real(dp), allocatable, intent(out) :: path(:) !! Simulated future levels.
      integer, intent(out) :: info !! Zero on success.
      real(dp), allocatable :: work(:)
      real(dp) :: value
      integer :: start, step, time, regime, lag_index

      info = 1
      start = size(history)
      if (model%info /= 0 .or. .not. allocated(model%ar_order) .or. &
         start < max(maxval(model%ar_order), model%delay) .or. &
         size(innovations) < 1) return
      allocate(work(start + size(innovations)), path(size(innovations)))
      work(:start) = history
      do step = 1, size(innovations)
         time = start + step
         regime = tvecm_regime(work(time - model%delay), model%thresholds)
         value = model%coefficients(regime, 1)
         do lag_index = 1, model%ar_order(regime)
            value = value + model%coefficients(regime, lag_index + 1)* &
               work(time - lag_index)
         end do
         value = value + innovations(step)
         if (step == 1 .and. apply_shock) value = value + shock
         work(time) = value
         path(step) = value
      end do
      info = 0
   end subroutine setar_impulse_path

   pure subroutine tvar_impact_matrix(model, orthogonalized, impact, info)
      !! Construct identity or pooled-covariance Cholesky TVAR impact shocks.
      type(nts_mtar_model_t), intent(in) :: model !! Fitted TVAR model.
      logical, intent(in) :: orthogonalized !! Whether to orthogonalize shocks.
      real(dp), allocatable, intent(out) :: impact(:, :) !! Variable-by-shock impacts.
      integer, intent(out) :: info !! Zero on success.
      real(dp), allocatable :: covariance(:, :)
      integer :: variables

      info = 1
      variables = size(model%data, 2)
      if (variables < 1) return
      if (.not. orthogonalized) then
         impact = identity_matrix(variables)
         info = 0
         return
      end if
      if (.not. allocated(model%residuals) .or. size(model%residuals, 1) < 1) return
      covariance = matmul(transpose(model%residuals), model%residuals)/ &
         real(size(model%residuals, 1), dp)
      call cholesky_lower(covariance, impact, info)
   end subroutine tvar_impact_matrix

   pure subroutine tvecm_impact_matrix(model, orthogonalized, impact, info)
      !! Construct identity or covariance-Cholesky TVECM impact shocks.
      type(tsdyn_tvecm_model_t), intent(in) :: model !! Fitted TVECM model.
      logical, intent(in) :: orthogonalized !! Whether to orthogonalize shocks.
      real(dp), allocatable, intent(out) :: impact(:, :) !! Variable-by-shock impacts.
      integer, intent(out) :: info !! Zero on success.
      integer :: variables

      info = 1
      if (.not. allocated(model%covariance)) return
      variables = size(model%covariance, 1)
      if (variables < 1 .or. size(model%covariance, 2) /= variables) return
      if (.not. orthogonalized) then
         impact = identity_matrix(variables)
         info = 0
         return
      end if
      call cholesky_lower(model%covariance, impact, info)
   end subroutine tvecm_impact_matrix

   pure subroutine tvar_fixed_regime_response(model, regime, horizon, impact, &
      response, info)
      !! Compute paired TVAR paths while holding the selected regime fixed.
      type(nts_mtar_model_t), intent(in) :: model !! Fitted TVAR model.
      integer, intent(in) :: regime !! Regime held fixed.
      integer, intent(in) :: horizon !! Largest response horizon, including zero.
      real(dp), intent(in) :: impact(:, :) !! Variable-by-shock impact matrix.
      real(dp), allocatable, intent(out) :: response(:, :, :) !! Conditional responses.
      integer, intent(out) :: info !! Zero on success.
      real(dp), allocatable :: history(:, :), innovations(:, :)
      real(dp), allocatable :: baseline(:, :), shocked(:, :)
      integer :: start, shock, status, variables

      info = 1
      variables = size(impact, 1)
      if (size(impact, 2) < 1 .or. variables /= size(model%data, 2)) return
      start = max(maxval(model%ar_order), model%delay)
      allocate(history(start, variables), innovations(horizon + 1, variables))
      allocate(response(horizon + 1, variables, size(impact, 2)))
      history = 0.0_dp
      innovations = 0.0_dp
      call tvar_impulse_path(model, history, innovations, [real(dp) ::], &
         baseline, status, regime)
      if (status /= 0) return
      do shock = 1, size(impact, 2)
         call tvar_impulse_path(model, history, innovations, impact(:, shock), &
            shocked, status, regime)
         if (status /= 0) return
         response(:, :, shock) = shocked - baseline
      end do
      info = 0
   end subroutine tvar_fixed_regime_response

   pure subroutine tvecm_fixed_regime_response(model, regime, horizon, impact, &
      response, info)
      !! Compute paired TVECM paths while holding the selected regime fixed.
      type(tsdyn_tvecm_model_t), intent(in) :: model !! Fitted TVECM model.
      integer, intent(in) :: regime !! Regime held fixed.
      integer, intent(in) :: horizon !! Largest response horizon, including zero.
      real(dp), intent(in) :: impact(:, :) !! Variable-by-shock impact matrix.
      real(dp), allocatable, intent(out) :: response(:, :, :) !! Conditional responses.
      integer, intent(out) :: info !! Zero on success.
      real(dp), allocatable :: history(:, :), innovations(:, :)
      real(dp), allocatable :: baseline(:, :), shocked(:, :)
      integer :: shock, status, variables

      info = 1
      variables = size(impact, 1)
      if (size(impact, 2) < 1 .or. variables /= size(model%data, 2)) return
      allocate(history(model%lag + 1, variables))
      allocate(innovations(horizon + 1, variables))
      allocate(response(horizon + 1, variables, size(impact, 2)))
      history = 0.0_dp
      innovations = 0.0_dp
      call tvecm_impulse_path(model, history, innovations, [real(dp) ::], &
         baseline, status, regime)
      if (status /= 0) return
      do shock = 1, size(impact, 2)
         call tvecm_impulse_path(model, history, innovations, impact(:, shock), &
            shocked, status, regime)
         if (status /= 0) return
         response(:, :, shock) = shocked - baseline
      end do
      info = 0
   end subroutine tvecm_fixed_regime_response

   pure subroutine tvar_impulse_path(model, history, innovations, shock, path, &
      info, fixed_regime)
      !! Simulate one TVAR path, optionally perturbing its first innovation.
      type(nts_mtar_model_t), intent(in) :: model !! Fitted threshold VAR model.
      real(dp), intent(in) :: history(:, :) !! Presample level history.
      real(dp), intent(in) :: innovations(:, :) !! Horizon-by-variable innovations.
      real(dp), intent(in) :: shock(:) !! First-period impact, or an empty vector.
      real(dp), allocatable, intent(out) :: path(:, :) !! Simulated future levels.
      integer, intent(out) :: info !! Zero on success.
      integer, intent(in), optional :: fixed_regime !! Regime override for conditional IRFs.
      real(dp), allocatable :: work(:, :), value(:)
      integer :: start, time, step, regime, lag_index

      info = 1
      start = size(history, 1)
      if (size(history, 2) /= size(innovations, 2) .or. &
         start < max(maxval(model%ar_order), model%delay) .or. &
         (size(shock) /= 0 .and. size(shock) /= size(innovations, 2))) return
      allocate(work(start + size(innovations, 1), size(innovations, 2)))
      allocate(path(size(innovations, 1), size(innovations, 2)))
      allocate(value(size(innovations, 2)))
      work(:start, :) = history
      do step = 1, size(innovations, 1)
         time = start + step
         regime = tvecm_regime(work(time - model%delay, &
            model%threshold_component), model%thresholds)
         if (present(fixed_regime)) regime = fixed_regime
         value = model%intercept(:, regime)
         if (allocated(model%trend)) then
            if (model%include_trend(regime)) value = value + &
               model%trend(:, regime)*real(size(model%residuals, 1) + step, dp)
         end if
         do lag_index = 1, model%ar_order(regime)
            value = value + matmul(model%ar(:, :, lag_index, regime), &
               work(time - lag_index, :))
         end do
         value = value + innovations(step, :)
         if (step == 1 .and. size(shock) > 0) value = value + shock
         work(time, :) = value
         path(step, :) = value
      end do
      info = 0
   end subroutine tvar_impulse_path

   pure subroutine tvecm_impulse_path(model, history, innovations, shock, path, &
      info, fixed_regime)
      !! Simulate one TVECM path, optionally perturbing its first innovation.
      type(tsdyn_tvecm_model_t), intent(in) :: model !! Fitted threshold VECM.
      real(dp), intent(in) :: history(:, :) !! Presample level history.
      real(dp), intent(in) :: innovations(:, :) !! Horizon-by-variable innovations.
      real(dp), intent(in) :: shock(:) !! First-period impact, or an empty vector.
      real(dp), allocatable, intent(out) :: path(:, :) !! Simulated future levels.
      integer, intent(out) :: info !! Zero on success.
      integer, intent(in), optional :: fixed_regime !! Regime override for conditional IRFs.
      real(dp), allocatable :: work(:, :), change(:)
      integer :: start, variables, time, step, regime, status, offset_time
      integer :: trend_index

      info = 1
      start = size(history, 1)
      variables = size(innovations, 2)
      if (variables < 2 .or. size(history, 2) /= variables .or. &
         start < model%lag + 1 .or. &
         (size(shock) /= 0 .and. size(shock) /= variables)) return
      allocate(work(start + size(innovations, 1), variables))
      allocate(path(size(innovations, 1), variables), change(variables))
      work(:start, :) = history
      do step = 1, size(innovations, 1)
         time = start + step
         offset_time = time - 1
         if (allocated(model%data)) offset_time = &
            size(model%data, 1) + step - 1
         trend_index = step
         if (allocated(model%residuals)) trend_index = &
            size(model%residuals, 1) + step
          call tvecm_next_change(model, work, time, trend_index, change, regime, &
             status, offset_time, fixed_regime)
         if (status /= 0) return
         change = change + innovations(step, :)
         if (step == 1 .and. size(shock) > 0) change = change + shock
         work(time, :) = work(time - 1, :) + change
         path(step, :) = work(time, :)
      end do
      info = 0
   end subroutine tvecm_impulse_path

   pure subroutine setar_regime_measures(series, order, delays, trim, &
      grid_count, momentum, measure, info)
      !! Compute residual sums of squares for linear, two-, and three-regime AR fits.
      real(dp), intent(in) :: series(:) !! Observed univariate series.
      integer, intent(in) :: order !! Common autoregressive order.
      integer, intent(in) :: delays(:) !! Candidate transition delays.
      real(dp), intent(in), optional :: trim !! Minimum threshold-regime fraction.
      integer, intent(in), optional :: grid_count !! Maximum threshold candidates.
      logical, intent(in), optional :: momentum !! Use first-difference transitions.
      real(dp), intent(out) :: measure(3) !! AR, two-regime, and three-regime RSS.
      integer, intent(out) :: info !! Zero on success.
      type(nts_tar_model_t) :: model
      real(dp), allocatable :: design(:, :), response(:), beta(:), standard_error(:)
      real(dp), allocatable :: residual(:)
      real(dp) :: rss
      integer :: row, lag_index, delay, status

      info = 1
      measure = huge(1.0_dp)
      if (order < 1 .or. size(series) <= 3*order + 5 .or. size(delays) < 1 .or. &
         minval(delays) < 1) return
      allocate(design(size(series) - order, order + 1))
      allocate(response(size(series) - order))
      design(:, 1) = 1.0_dp
      do row = 1, size(response)
         response(row) = series(order + row)
         do lag_index = 1, order
            design(row, lag_index + 1) = series(order + row - lag_index)
         end do
      end do
      call ols_fit(design, response, beta, standard_error, residual, rss, status)
      if (status /= 0 .or. rss <= 0.0_dp) return
      measure(1) = rss
      do delay = 1, size(delays)
         model = tsdyn_setar_fit(series, [order, order], 1, delays(delay), trim, &
            grid_count, momentum)
         if (model%info == 0) measure(2) = min(measure(2), sum(model%residuals**2))
         model = tsdyn_setar_fit(series, [order, order, order], 2, delays(delay), &
            trim, grid_count, momentum)
         if (model%info == 0) measure(3) = min(measure(3), sum(model%residuals**2))
      end do
      if (any(measure == huge(1.0_dp))) then
         info = 2
         return
      end if
      info = 0
   end subroutine setar_regime_measures

   pure subroutine tvar_regime_measures(series, order, threshold_component, &
      delays, trim, grid_count, momentum, measure, info)
      !! Compute covariance log determinants for linear, two-, and three-regime VAR fits.
      real(dp), intent(in) :: series(:, :) !! Time-by-variable observations.
      integer, intent(in) :: order !! Common VAR order.
      integer, intent(in) :: threshold_component !! Transition variable index.
      integer, intent(in) :: delays(:) !! Candidate transition delays.
      real(dp), intent(in), optional :: trim !! Minimum threshold-regime fraction.
      integer, intent(in), optional :: grid_count !! Maximum threshold candidates.
      logical, intent(in), optional :: momentum !! Use first-difference transitions.
      real(dp), intent(out) :: measure(3) !! VAR, two-regime, and three-regime log determinants.
      integer, intent(out) :: info !! Zero on success.
      type(nts_mtar_model_t) :: model
      real(dp), allocatable :: design(:, :), response(:, :), coefficient(:, :)
      real(dp), allocatable :: fitted(:, :), residual(:, :)
      real(dp) :: rss, value
      integer :: variables, row, lag_index, offset, delay, status

      info = 1
      measure = huge(1.0_dp)
      variables = size(series, 2)
      if (order < 1 .or. size(series, 1) <= 3*variables*order + 8 .or. &
         threshold_component < 1 .or. threshold_component > variables .or. &
         size(delays) < 1 .or. minval(delays) < 1) return
      allocate(design(size(series, 1) - order, 1 + variables*order))
      allocate(response(size(series, 1) - order, variables))
      design(:, 1) = 1.0_dp
      do row = 1, size(response, 1)
         response(row, :) = series(order + row, :)
         do lag_index = 1, order
            offset = 1 + (lag_index - 1)*variables
            design(row, offset + 1:offset + variables) = &
               series(order + row - lag_index, :)
         end do
      end do
      call tvecm_multivariate_regression(design, response, coefficient, fitted, &
         residual, rss, status)
      if (status /= 0) return
      call residual_logdet(residual, measure(1), status)
      if (status /= 0) return
      do delay = 1, size(delays)
         model = tsdyn_tvar_fit(series, [order, order], 1, threshold_component, &
            delays(delay), trim, grid_count, momentum)
         if (model%info == 0) then
            call residual_logdet(model%residuals, value, status)
            if (status == 0) measure(2) = min(measure(2), value)
         end if
         model = tsdyn_tvar_fit(series, [order, order, order], 2, &
            threshold_component, delays(delay), trim, grid_count, momentum)
         if (model%info == 0) then
            call residual_logdet(model%residuals, value, status)
            if (status == 0) measure(3) = min(measure(3), value)
         end if
      end do
      if (any(measure == huge(1.0_dp))) then
         info = 2
         return
      end if
      info = 0
   end subroutine tvar_regime_measures

   pure subroutine residual_logdet(residual, logdet, info)
      !! Compute the log determinant of a residual covariance matrix.
      real(dp), intent(in) :: residual(:, :) !! Time-by-equation residual matrix.
      real(dp), intent(out) :: logdet !! Covariance log determinant.
      integer, intent(out) :: info !! Zero on success.
      real(dp), allocatable :: covariance(:, :), inverse(:, :)

      covariance = matmul(transpose(residual), residual)/real(size(residual, 1), dp)
      allocate(inverse(size(residual, 2), size(residual, 2)))
      call inverse_logdet(covariance, inverse, logdet, info, 1.0e-12_dp)
   end subroutine residual_logdet

   pure subroutine regime_count_statistics(measure, scale, logarithmic, statistic)
      !! Form 1-vs-2, 1-vs-3, and 2-vs-3 threshold-model statistics.
      real(dp), intent(in) :: measure(3) !! Ordered model fit measures.
      real(dp), intent(in) :: scale !! Effective observation count.
      logical, intent(in) :: logarithmic !! Measures are covariance log determinants.
      real(dp), intent(out) :: statistic(3) !! Three nested-model statistics.

      if (logarithmic) then
         statistic = scale*[measure(1) - measure(2), measure(1) - measure(3), &
            measure(2) - measure(3)]
      else
         statistic = scale*[(measure(1) - measure(2))/measure(2), &
            (measure(1) - measure(3))/measure(3), &
            (measure(2) - measure(3))/measure(3)]
      end if
      statistic = max(0.0_dp, statistic)
   end subroutine regime_count_statistics

   pure subroutine bbc_statistic_path(series, order, method, trim, use_counts, &
      thresholds, path, selected, statistic, info)
      !! Evaluate a BBC statistic across admissible symmetric threshold pairs.
      real(dp), intent(in) :: series(:) !! Centered univariate series.
      integer, intent(in) :: order !! Number of lagged differences.
      character(len=*), intent(in) :: method !! Selected BBC statistic.
      real(dp), intent(in) :: trim !! Minimum regime fraction.
      logical, intent(in) :: use_counts !! Enforce an absolute count criterion.
      real(dp), allocatable, intent(out) :: thresholds(:) !! Evaluated positive thresholds.
      real(dp), allocatable, intent(out) :: path(:) !! Statistic at each threshold.
      integer, intent(out) :: selected !! Maximizing threshold index.
      real(dp), intent(out) :: statistic !! Supremum statistic.
      integer, intent(out) :: info !! Zero on success.
      real(dp), allocatable :: response(:), level(:), differences(:, :)
      real(dp), allocatable :: ordered(:), threshold_work(:), path_work(:)
      real(dp) :: value
      integer :: candidate, path_count, minimum, status

      info = 1
      selected = 0
      statistic = 0.0_dp
      call unit_root_regression_data(series, order, response, level, differences, status)
      if (status /= 0) return
      ordered = sorted(abs(level))
      allocate(threshold_work(size(ordered)), path_work(size(ordered)))
      path_count = 0
      minimum = max(order + 2, ceiling(trim*real(size(level), dp)))
      if (.not. use_counts) minimum = ceiling(trim*real(size(level), dp))
      do candidate = 1, size(ordered)
         if (candidate > 1) then
            if (ordered(candidate) <= ordered(candidate - 1)) cycle
         end if
         value = ordered(candidate)
         if (count(level <= -value) < minimum .or. &
            count(level > value) < minimum .or. &
            count(level > -value .and. level <= value) < minimum) cycle
         call bbc_statistic_at(response, level, differences, value, method, &
            path_work(path_count + 1), status)
         if (status /= 0) cycle
         path_count = path_count + 1
         threshold_work(path_count) = value
      end do
      if (path_count == 0) then
         allocate(thresholds(0), path(0))
         info = 2
         return
      end if
      thresholds = threshold_work(:path_count)
      path = path_work(:path_count)
      selected = maxloc(path, dim=1)
      statistic = path(selected)
      info = 0
   end subroutine bbc_statistic_path

   pure subroutine bbc_statistic_value(series, order, method, trim, use_counts, &
      statistic, info)
      !! Return only the supremum BBC statistic for bootstrap evaluation.
      real(dp), intent(in) :: series(:) !! Centered bootstrap series.
      integer, intent(in) :: order !! Number of lagged differences.
      character(len=*), intent(in) :: method !! Selected BBC statistic.
      real(dp), intent(in) :: trim !! Minimum regime fraction.
      logical, intent(in) :: use_counts !! Enforce an absolute count criterion.
      real(dp), intent(out) :: statistic !! Supremum statistic.
      integer, intent(out) :: info !! Zero on success.
      real(dp), allocatable :: thresholds(:), path(:)
      integer :: selected

      call bbc_statistic_path(series, order, method, trim, use_counts, thresholds, &
         path, selected, statistic, info)
   end subroutine bbc_statistic_value

   pure subroutine bbc_statistic_at(response, level, differences, threshold, &
      method, statistic, info)
      !! Evaluate one BBC LR, Wald, or LM statistic.
      real(dp), intent(in) :: response(:) !! Aligned first differences.
      real(dp), intent(in) :: level(:) !! Lagged levels.
      real(dp), intent(in) :: differences(:, :) !! Lagged first differences.
      real(dp), intent(in) :: threshold !! Positive symmetric threshold.
      character(len=*), intent(in) :: method !! `lr`, `wald`, or `lm`.
      real(dp), intent(out) :: statistic !! Candidate statistic.
      integer, intent(out) :: info !! Zero on success.
      real(dp), allocatable :: unrestricted(:, :), restricted(:, :)
      real(dp), allocatable :: beta_u(:), beta_r(:), se(:), residual_u(:), residual_r(:)
      real(dp), allocatable :: inverse(:, :), restriction_covariance(:, :)
      real(dp), allocatable :: restriction_inverse(:, :), restricted_level(:), score(:)
      real(dp) :: rss_u, rss_r, sigma2
      integer :: rows, order, row, regime, block_u, block_r, lag_index, status

      info = 1
      statistic = 0.0_dp
      rows = size(response)
      order = size(differences, 2)
      allocate(unrestricted(rows, 3*(order + 2)))
      allocate(restricted(rows, 3*(order + 1)))
      unrestricted = 0.0_dp
      restricted = 0.0_dp
      do row = 1, rows
         regime = tvecm_regime(level(row), [-threshold, threshold])
         block_u = (regime - 1)*(order + 2)
         block_r = (regime - 1)*(order + 1)
         unrestricted(row, block_u + 1) = 1.0_dp
         unrestricted(row, block_u + 2) = level(row)
         restricted(row, block_r + 1) = 1.0_dp
         do lag_index = 1, order
            unrestricted(row, block_u + 2 + lag_index) = &
               differences(row, lag_index)
            restricted(row, block_r + 1 + lag_index) = &
               differences(row, lag_index)
         end do
      end do
      call ols_fit(unrestricted, response, beta_u, se, residual_u, rss_u, status)
      if (status /= 0 .or. rss_u <= 0.0_dp) return
      call ols_fit(restricted, response, beta_r, se, residual_r, rss_r, status)
      if (status /= 0 .or. rss_r <= 0.0_dp) return
      select case (method)
      case ('lr')
         statistic = real(rows, dp)*log(rss_r/rss_u)
      case ('wald')
         call invert_matrix(matmul(transpose(unrestricted), unrestricted), &
            inverse, status)
         if (status /= 0) return
         allocate(restriction_covariance(3, 3), restricted_level(3))
         do regime = 1, 3
            restricted_level(regime) = beta_u((regime - 1)*(order + 2) + 2)
            do lag_index = 1, 3
               restriction_covariance(regime, lag_index) = inverse( &
                  (regime - 1)*(order + 2) + 2, &
                  (lag_index - 1)*(order + 2) + 2)
            end do
         end do
         call invert_matrix(restriction_covariance, restriction_inverse, status)
         if (status /= 0) return
         sigma2 = rss_u/real(rows, dp)
         statistic = dot_product(restricted_level, &
            matmul(restriction_inverse, restricted_level))/sigma2
      case ('lm')
         call invert_matrix(matmul(transpose(unrestricted), unrestricted), &
            inverse, status)
         if (status /= 0) return
         score = matmul(transpose(unrestricted), residual_r)
         sigma2 = rss_r/real(rows, dp)
         statistic = dot_product(score, matmul(inverse, score))/sigma2
      end select
      statistic = max(0.0_dp, statistic)
      info = 0
   end subroutine bbc_statistic_at

   pure subroutine kapshin_statistic_path(series, order, grid_scale, &
      grid_exponent, requested_points, minimum_middle, thresholds, path, &
      selected, statistic, info)
      !! Evaluate Kapetanios-Shin Wald statistics across an asymmetric grid.
      real(dp), intent(in) :: series(:) !! Deterministically adjusted series.
      integer, intent(in) :: order !! Number of lagged differences.
      real(dp), intent(in) :: grid_scale !! Central-grid scale parameter.
      real(dp), intent(in) :: grid_exponent !! Central-grid rate exponent.
      integer, intent(in) :: requested_points !! Requested marginal point count or zero.
      integer, intent(in) :: minimum_middle !! Minimum middle-regime observations.
      real(dp), allocatable, intent(out) :: thresholds(:, :) !! Evaluated lower-upper pairs.
      real(dp), allocatable, intent(out) :: path(:) !! Wald statistic path.
      integer, intent(out) :: selected !! Supremum pair index.
      real(dp), intent(out) :: statistic(3) !! Supremum, average, and exponential average.
      integer, intent(out) :: info !! Zero on success.
      real(dp), allocatable :: response(:), level(:), differences(:, :), ordered(:)
      real(dp), allocatable :: grid(:), threshold_work(:, :), path_work(:)
      real(dp) :: p1, p2, probability
      integer :: points, middle, lower_index, upper_index, pair_count, status

      info = 1
      selected = 0
      statistic = 0.0_dp
      call unit_root_regression_data(series, order, response, level, differences, status)
      if (status /= 0 .or. grid_scale <= 0.0_dp .or. grid_exponent <= 0.0_dp .or. &
         minimum_middle < 1) return
      p1 = 0.5_dp - grid_scale/real(size(level), dp)**grid_exponent
      p2 = 0.5_dp + grid_scale/real(size(level), dp)**grid_exponent
      p1 = max(0.0_dp, p1)
      p2 = min(1.0_dp, p2)
      points = requested_points
      if (points <= 0) points = max(3, nint((p2 - p1)*real(size(level), dp)))
      if (points < 2 .or. p2 <= p1) then
         info = 2
         return
      end if
      ordered = sorted(level)
      allocate(grid(points))
      do lower_index = 1, points
         probability = p1 + (p2 - p1)*real(lower_index - 1, dp)/ &
            real(points - 1, dp)
         grid(lower_index) = quantile(ordered, probability)
      end do
      middle = (points + 1)/2
      allocate(threshold_work(points*points, 2), path_work(points*points))
      pair_count = 0
      do lower_index = 1, middle
         do upper_index = middle, points
            if (grid(upper_index) <= grid(lower_index)) cycle
            if (count(level > grid(lower_index) .and. &
               level <= grid(upper_index)) < minimum_middle) cycle
            call kapshin_wald_at(response, level, differences, grid(lower_index), &
               grid(upper_index), path_work(pair_count + 1), status)
            if (status /= 0) cycle
            pair_count = pair_count + 1
            threshold_work(pair_count, :) = &
               [grid(lower_index), grid(upper_index)]
         end do
      end do
      if (pair_count == 0) then
         allocate(thresholds(0, 2), path(0))
         info = 3
         return
      end if
      thresholds = threshold_work(:pair_count, :)
      path = path_work(:pair_count)
      selected = maxloc(path, dim=1)
      statistic(1) = path(selected)
      statistic(2) = sum(path)/real(pair_count, dp)
      statistic(3) = sum(exp(min(700.0_dp, path/2.0_dp)))/real(pair_count, dp)
      info = 0
   end subroutine kapshin_statistic_path

   pure subroutine kapshin_statistic_value(series, order, grid_scale, &
      grid_exponent, points, minimum_middle, statistic, info)
      !! Return only the three Kapetanios-Shin aggregate statistics.
      real(dp), intent(in) :: series(:) !! Adjusted bootstrap series.
      integer, intent(in) :: order !! Number of lagged differences.
      real(dp), intent(in) :: grid_scale !! Central-grid scale parameter.
      real(dp), intent(in) :: grid_exponent !! Central-grid rate exponent.
      integer, intent(in) :: points !! Marginal threshold point count or zero.
      integer, intent(in) :: minimum_middle !! Minimum middle-regime observations.
      real(dp), intent(out) :: statistic(3) !! Supremum, average, and exponential average.
      integer, intent(out) :: info !! Zero on success.
      real(dp), allocatable :: thresholds(:, :), path(:)
      integer :: selected

      call kapshin_statistic_path(series, order, grid_scale, grid_exponent, &
         points, minimum_middle, thresholds, path, selected, statistic, info)
   end subroutine kapshin_statistic_value

   pure subroutine kapshin_wald_at(response, level, differences, lower, upper, &
      statistic, info)
      !! Evaluate the outer-regime level-coefficient Wald statistic.
      real(dp), intent(in) :: response(:) !! Aligned first differences.
      real(dp), intent(in) :: level(:) !! Lagged levels.
      real(dp), intent(in) :: differences(:, :) !! Lagged first differences.
      real(dp), intent(in) :: lower !! Lower threshold.
      real(dp), intent(in) :: upper !! Upper threshold.
      real(dp), intent(out) :: statistic !! Wald statistic.
      integer, intent(out) :: info !! Zero on success.
      real(dp), allocatable :: design(:, :), beta(:), standard_error(:), residual(:)
      real(dp), allocatable :: inverse(:, :), covariance(:, :), covariance_inverse(:, :)
      real(dp) :: rss, sigma2
      integer :: status

      info = 1
      statistic = 0.0_dp
      allocate(design(size(response), 2 + size(differences, 2)))
      design(:, 1) = merge(level, 0.0_dp, level <= lower)
      design(:, 2) = merge(level, 0.0_dp, level > upper)
      if (size(differences, 2) > 0) design(:, 3:) = differences
      call ols_fit(design, response, beta, standard_error, residual, rss, status)
      if (status /= 0 .or. rss <= 0.0_dp) return
      call invert_matrix(matmul(transpose(design), design), inverse, status)
      if (status /= 0) return
      covariance = inverse(:2, :2)
      call invert_matrix(covariance, covariance_inverse, status)
      if (status /= 0) return
      sigma2 = rss/real(size(beta), dp)
      statistic = dot_product(beta(:2), matmul(covariance_inverse, beta(:2)))/sigma2
      statistic = max(0.0_dp, statistic)
      info = 0
   end subroutine kapshin_wald_at

   pure subroutine unit_root_regression_data(series, order, response, level, &
      differences, info)
      !! Align first differences, lagged levels, and lagged differences.
      real(dp), intent(in) :: series(:) !! Univariate level series.
      integer, intent(in) :: order !! Number of lagged differences.
      real(dp), allocatable, intent(out) :: response(:) !! Current first differences.
      real(dp), allocatable, intent(out) :: level(:) !! Lagged levels.
      real(dp), allocatable, intent(out) :: differences(:, :) !! Lagged first differences.
      integer, intent(out) :: info !! Zero on success.
      integer :: rows, row, time, lag_index

      info = 1
      if (order < 0 .or. size(series) <= order + 4 .or. &
         .not. all(ieee_is_finite(series))) return
      rows = size(series) - order - 1
      allocate(response(rows), level(rows), differences(rows, order))
      do row = 1, rows
         time = order + 1 + row
         response(row) = series(time) - series(time - 1)
         level(row) = series(time - 1)
         do lag_index = 1, order
            differences(row, lag_index) = series(time - lag_index) - &
               series(time - lag_index - 1)
         end do
      end do
      info = 0
   end subroutine unit_root_regression_data

   pure subroutine remove_deterministic(series, deterministic, adjusted, info)
      !! Remove the deterministic terms used by the Kapetanios-Shin cases.
      real(dp), intent(in) :: series(:) !! Original series.
      character(len=*), intent(in) :: deterministic !! Deterministic-case name.
      real(dp), allocatable, intent(out) :: adjusted(:) !! Residualized series.
      integer, intent(out) :: info !! Zero on success.
      real(dp), allocatable :: design(:, :), beta(:), standard_error(:), residual(:)
      real(dp) :: rss
      integer :: time, status

      info = 0
      if (deterministic == 'none') then
         adjusted = series
         return
      end if
      if (deterministic == 'const') then
         adjusted = series - sum(series)/real(size(series), dp)
         return
      end if
      if (deterministic == 'trend') then
         allocate(design(size(series), 1))
         design(:, 1) = [(real(time, dp), time=1, size(series))]
      else if (deterministic == 'both') then
         allocate(design(size(series), 2))
         design(:, 1) = 1.0_dp
         design(:, 2) = [(real(time, dp), time=1, size(series))]
      else
         info = 1
         return
      end if
      call ols_fit(design, series, beta, standard_error, residual, rss, status)
      if (status /= 0) then
         info = 2
         return
      end if
      adjusted = residual
   end subroutine remove_deterministic

   pure subroutine tvecm_regression_data(data, lag, cointegration, include_constant, &
      include_trend, response, common, error_correction, cointegration_offset)
      !! Construct aligned TVECM changes, common regressors, and lagged ECT.
      real(dp), intent(in) :: data(:, :) !! Time-by-variable level observations.
      integer, intent(in) :: lag !! Number of lagged differences.
      real(dp), intent(in) :: cointegration(:) !! Normalized cointegrating vector.
      logical, intent(in) :: include_constant !! Whether to include an intercept.
      logical, intent(in) :: include_trend !! Whether to include a trend.
      real(dp), allocatable, intent(out) :: response(:, :) !! Aligned level changes.
      real(dp), allocatable, intent(out) :: common(:, :) !! Deterministic and lagged changes.
      real(dp), allocatable, intent(out) :: error_correction(:) !! Lagged ECT values.
      real(dp), intent(in), optional :: cointegration_offset(:) !! Long-run contribution by input time.
      integer :: rows, variables, deterministic, row, time, column, lag_index

      rows = size(data, 1) - lag - 1
      variables = size(data, 2)
      deterministic = merge(1, 0, include_constant) + merge(1, 0, include_trend)
      allocate(response(rows, variables), &
         common(rows, deterministic + variables*lag))
      allocate(error_correction(rows))
      do row = 1, rows
         time = lag + 1 + row
         response(row, :) = data(time, :) - data(time - 1, :)
         error_correction(row) = dot_product(cointegration, data(time - 1, :))
         if (present(cointegration_offset)) error_correction(row) = &
            error_correction(row) - cointegration_offset(time - 1)
         column = 0
         if (include_constant) then
            column = column + 1
            common(row, column) = 1.0_dp
         end if
         if (include_trend) then
            column = column + 1
            common(row, column) = real(row, dp)
         end if
         do lag_index = 1, lag
            common(row, column + variables*(lag_index - 1) + 1: &
               column + variables*lag_index) = &
               data(time - lag_index, :) - data(time - lag_index - 1, :)
         end do
      end do
   end subroutine tvecm_regression_data

   pure subroutine tvecm_design(error_correction, common, thresholds, only_ect, &
      trim, design, regime, info, shared_columns)
      !! Build a trimmed regime-interacted TVECM regression matrix.
      real(dp), intent(in) :: error_correction(:) !! Lagged error-correction values.
      real(dp), intent(in) :: common(:, :) !! Deterministic and lagged-change regressors.
      real(dp), intent(in) :: thresholds(:) !! One or two ordered thresholds.
      logical, intent(in) :: only_ect !! Whether only ECT adjustment switches.
      real(dp), intent(in) :: trim !! Minimum regime fraction.
      real(dp), allocatable, intent(out) :: design(:, :) !! Regime-interacted design.
      integer, allocatable, intent(out) :: regime(:) !! Regime assignment by row.
      integer, intent(out) :: info !! Zero on success.
      integer, intent(in), optional :: shared_columns !! Leading common columns shared across regimes.
      integer :: regimes, shared, row, selected, block, columns, switching_common

      info = 1
      regimes = size(thresholds) + 1
      shared = 0
      if (present(shared_columns)) shared = shared_columns
      if (size(common, 1) /= size(error_correction) .or. &
         (size(thresholds) /= 1 .and. size(thresholds) /= 2) .or. &
         shared < 0 .or. shared > size(common, 2)) return
      if (size(thresholds) == 2) then
         if (thresholds(2) <= thresholds(1)) return
      end if
      allocate(regime(size(error_correction)))
      do row = 1, size(error_correction)
         regime(row) = tvecm_regime(error_correction(row), thresholds)
      end do
      do selected = 1, regimes
         if (real(count(regime == selected), dp)/real(size(regime), dp) <= trim) return
      end do
      if (only_ect) then
         columns = regimes + size(common, 2)
         allocate(design(size(error_correction), columns))
         design = 0.0_dp
         do row = 1, size(error_correction)
            design(row, regime(row)) = error_correction(row)
            design(row, regimes + 1:) = common(row, :)
         end do
      else if (shared > 0) then
         switching_common = size(common, 2) - shared
         columns = shared + regimes*(1 + switching_common)
         allocate(design(size(error_correction), columns))
         design = 0.0_dp
         design(:, :shared) = common(:, :shared)
         do row = 1, size(error_correction)
            block = shared + (regime(row) - 1)*(1 + switching_common)
            design(row, block + 1) = error_correction(row)
            if (switching_common > 0) then
               design(row, block + 2:block + 1 + switching_common) = &
                  common(row, shared + 1:)
            end if
         end do
      else
         columns = regimes*(1 + size(common, 2))
         allocate(design(size(error_correction), columns))
         design = 0.0_dp
         do row = 1, size(error_correction)
            block = (regime(row) - 1)*(1 + size(common, 2))
            design(row, block + 1) = error_correction(row)
            design(row, block + 2:block + 1 + size(common, 2)) = common(row, :)
         end do
      end if
      info = 0
   end subroutine tvecm_design

   pure subroutine tvecm_multivariate_regression(design, response, coefficient, &
      fitted, residual, rss, info)
      !! Fit all TVECM change equations by common-design least squares.
      real(dp), intent(in) :: design(:, :) !! TVECM regression matrix.
      real(dp), intent(in) :: response(:, :) !! Aligned multivariate changes.
      real(dp), allocatable, intent(out) :: coefficient(:, :) !! Equation-by-regressor coefficients.
      real(dp), allocatable, intent(out) :: fitted(:, :) !! Fitted multivariate changes.
      real(dp), allocatable, intent(out) :: residual(:, :) !! Multivariate residuals.
      real(dp), intent(out) :: rss !! Sum of squared residuals over equations.
      integer, intent(out) :: info !! Zero on success.
      real(dp), allocatable :: cross_product(:, :), inverse(:, :)
      integer :: status

      info = 1
      rss = huge(1.0_dp)
      if (size(design, 1) /= size(response, 1) .or. &
         size(design, 1) <= size(design, 2)) return
      cross_product = matmul(transpose(design), design)
      call invert_matrix(cross_product, inverse, status)
      if (status /= 0) return
      coefficient = transpose(matmul(inverse, matmul(transpose(design), response)))
      fitted = matmul(design, transpose(coefficient))
      residual = response - fitted
      rss = sum(residual*residual)
      if (.not. ieee_is_finite(rss)) return
      info = 0
   end subroutine tvecm_multivariate_regression

   pure subroutine tvecm_next_change(model, history, time, trend, change, &
      regime, info, offset_time, fixed_regime)
      !! Evaluate one conditional TVECM level change from available history.
      type(tsdyn_tvecm_model_t), intent(in) :: model !! Fitted TVECM model.
      real(dp), intent(in) :: history(:, :) !! Level history including output capacity.
      integer, intent(in) :: time !! One-based row to be generated.
      integer, intent(in) :: trend !! Deterministic trend index.
      real(dp), intent(out) :: change(:) !! Conditional multivariate change.
      integer, intent(out) :: regime !! Selected threshold regime.
      integer, intent(out) :: info !! Zero on success.
      integer, intent(in), optional :: offset_time !! Time index used for the lagged long-run offset.
      integer, intent(in), optional :: fixed_regime !! Regime override for conditional IRFs.
      real(dp), allocatable :: common(:), regressor(:)
      real(dp), allocatable :: cointegration(:)
      real(dp) :: ect
      integer :: variables, deterministic, switching_common
      integer :: column, lag_index, block, width, selected_offset_time

      info = 1
      change = 0.0_dp
      variables = size(history, 2)
      if (size(change) /= variables .or. variables < 2 .or. &
         time <= model%lag + 1 .or. &
         time > size(history, 1)) return
      cointegration = tvecm_cointegration_vector(model, variables)
      if (size(cointegration) /= variables) return
      selected_offset_time = time - 1
      if (present(offset_time)) selected_offset_time = offset_time
      ect = dot_product(cointegration, history(time - 1, :)) - &
         tvecm_cointegration_offset_value(model, selected_offset_time)
      regime = tvecm_regime(ect, model%thresholds)
      if (present(fixed_regime)) regime = fixed_regime
      deterministic = merge(1, 0, model%include_constant) + &
         merge(1, 0, model%include_trend)
      allocate(common(deterministic + variables*model%lag))
      column = 0
      if (model%include_constant) then
         column = column + 1
         common(column) = 1.0_dp
      end if
      if (model%include_trend) then
         column = column + 1
         common(column) = real(trend, dp)
      end if
      do lag_index = 1, model%lag
         common(column + variables*(lag_index - 1) + 1: &
            column + variables*lag_index) = &
            history(time - lag_index, :) - history(time - lag_index - 1, :)
      end do
      if (model%only_error_correction) then
         change = model%coefficients(:, regime)*ect + matmul( &
            model%coefficients(:, size(model%thresholds) + 2:), common)
      else if (model%shared_deterministic .and. deterministic > 0) then
         switching_common = size(common) - deterministic
         width = 1 + switching_common
         block = deterministic + (regime - 1)*width
         allocate(regressor(width))
         regressor = [ect, common(deterministic + 1:)]
         change = matmul(model%coefficients(:, :deterministic), &
            common(:deterministic)) + matmul( &
            model%coefficients(:, block + 1:block + width), regressor)
      else
         width = 1 + size(common)
         block = (regime - 1)*width
         allocate(regressor(width))
         regressor = [ect, common]
         change = matmul(model%coefficients(:, block + 1:block + width), regressor)
      end if
      info = 0
   end subroutine tvecm_next_change

   pure function tvecm_cointegration_vector(model, variables) result(vector)
      !! Return the normalized cointegrating vector stored by a TVECM model.
      type(tsdyn_tvecm_model_t), intent(in) :: model !! TVECM model specification.
      integer, intent(in) :: variables !! Required number of vector elements.
      real(dp), allocatable :: vector(:)

      if (allocated(model%cointegration)) then
         if (size(model%cointegration) == variables) then
            vector = model%cointegration
         else
            allocate(vector(0))
         end if
      else if (variables == 2) then
         vector = [1.0_dp, -model%beta]
      else
         allocate(vector(0))
      end if
   end function tvecm_cointegration_vector

   pure real(dp) function tvecm_cointegration_offset_value(model, time) &
      result(value)
      !! Return a stored long-run offset, holding its final value after sample.
      type(tsdyn_tvecm_model_t), intent(in) :: model !! TVECM model specification.
      integer, intent(in) :: time !! Requested one-based offset time.

      value = 0.0_dp
      if (.not. allocated(model%cointegration_offset)) return
      if (size(model%cointegration_offset) < 1 .or. time < 1) return
      value = model%cointegration_offset(min(time, &
         size(model%cointegration_offset)))
   end function tvecm_cointegration_offset_value

   pure integer function tvecm_regime(error_correction, thresholds) result(regime)
      !! Map one error-correction value to its threshold regime.
      real(dp), intent(in) :: error_correction !! Error-correction value.
      real(dp), intent(in) :: thresholds(:) !! Ordered threshold values.

      regime = 1
      if (size(thresholds) >= 1) then
         if (error_correction > thresholds(1)) regime = 2
      end if
      if (size(thresholds) >= 2) then
         if (error_correction > thresholds(2)) regime = 3
      end if
   end function tvecm_regime

   pure subroutine setar_transition_series(series, delay, momentum, external, &
      transition, first_usable)
      !! Build the level or momentum transition series used by threshold models.
      real(dp), intent(in) :: series(:) !! Series supplying the default transition variable.
      integer, intent(in) :: delay !! Transition delay.
      logical, intent(in) :: momentum !! Whether to use first differences.
      real(dp), intent(in), optional :: external(:) !! External transition variable.
      real(dp), allocatable, intent(out) :: transition(:) !! Complete transition series.
      integer, intent(out) :: first_usable !! First response index with a valid transition.
      integer :: time

      first_usable = delay + 1
      if (delay < 1 .or. size(series) < delay + 2) return
      if (present(external)) then
         if (size(external) /= size(series) .or. &
            .not. all(ieee_is_finite(external))) return
         transition = external
      else if (momentum) then
         allocate(transition(size(series)))
         transition(1) = 0.0_dp
         do time = 2, size(series)
            transition(time) = series(time) - series(time - 1)
         end do
         first_usable = delay + 2
      else
         transition = series
      end if
   end subroutine setar_transition_series

   pure subroutine threshold_candidate_grid(ordered, trim, maximum, grid)
      !! Select an approximately uniform unique grid from trimmed order statistics.
      real(dp), intent(in) :: ordered(:) !! Sorted transition values.
      real(dp), intent(in) :: trim !! Fraction removed from each tail.
      integer, intent(in) :: maximum !! Maximum candidate count.
      real(dp), allocatable, intent(out) :: grid(:) !! Unique threshold candidates.
      real(dp), allocatable :: work(:)
      real(dp) :: probability, value
      integer :: candidate, count

      if (size(ordered) < 2 .or. maximum < 1) then
         allocate(grid(0))
         return
      end if
      allocate(work(maximum))
      count = 0
      do candidate = 1, maximum
         if (maximum == 1) then
            probability = 0.5_dp
         else
            probability = trim + (1.0_dp - 2.0_dp*trim)* &
               real(candidate - 1, dp)/real(maximum - 1, dp)
         end if
         value = quantile(ordered, probability)
         if (count == 0) then
            count = 1
            work(count) = value
         else if (value > work(count)) then
            count = count + 1
            work(count) = value
         end if
      end do
      allocate(grid(count))
      grid = work(:count)
   end subroutine threshold_candidate_grid

   pure subroutine aar_embedding_data(series, order, delay, steps, response, &
      lag_data)
      !! Construct aligned responses and delay-embedded AAR lag predictors.
      real(dp), intent(in) :: series(:) !! Complete univariate series.
      integer, intent(in) :: order !! Number of lag components.
      integer, intent(in) :: delay !! Delay between components.
      integer, intent(in) :: steps !! Forecast step.
      real(dp), allocatable, intent(out) :: response(:) !! Aligned future responses.
      real(dp), allocatable, intent(out) :: lag_data(:, :) !! Aligned lag predictors.
      integer :: rows, row, component, start

      rows = size(series) - (order - 1)*delay - steps
      allocate(response(rows), lag_data(rows, order))
      do row = 1, rows
         start = row
         response(row) = series(start + (order - 1)*delay + steps)
         do component = 1, order
            lag_data(row, component) = series(start + (order - component)*delay)
         end do
      end do
   end subroutine aar_embedding_data

   pure function aar_penalized_fit(design, response, base_penalty, lambda, &
      basis_count) result(fit)
      !! Fit an additive spline design for supplied component penalties.
      real(dp), intent(in) :: design(:, :) !! Intercept and centered spline design.
      real(dp), intent(in) :: response(:) !! Aligned response vector.
      real(dp), intent(in) :: base_penalty(:, :) !! One-component spline penalty.
      real(dp), intent(in) :: lambda(:) !! Component smoothing parameters.
      integer, intent(in) :: basis_count !! Basis functions per component.
      type(penalized_regression_t) :: fit
      real(dp), allocatable :: penalty(:, :)
      integer :: component, first, last

      allocate(penalty(size(design, 2), size(design, 2)))
      penalty = 0.0_dp
      do component = 1, size(lambda)
         first = 2 + (component - 1)*basis_count
         last = first + basis_count - 1
         penalty(first:last, first:last) = lambda(component)*base_penalty
      end do
      fit = penalized_regression_fit(design, response, penalty, 1.0_dp)
   end function aar_penalized_fit

   pure subroutine aar_component_degrees(design, base_penalty, lambda, &
      basis_count, component_df)
      !! Decompose additive-model effective degrees of freedom by lag component.
      real(dp), intent(in) :: design(:, :) !! Intercept and centered spline design.
      real(dp), intent(in) :: base_penalty(:, :) !! One-component spline penalty.
      real(dp), intent(in) :: lambda(:) !! Component smoothing parameters.
      integer, intent(in) :: basis_count !! Basis functions per component.
      real(dp), allocatable, intent(out) :: component_df(:) !! Component effective degrees.
      real(dp), allocatable :: penalty(:, :), cross_product(:, :), inverse(:, :)
      real(dp), allocatable :: influence(:, :)
      integer :: component, first, last, index, info

      allocate(component_df(size(lambda)))
      component_df = 0.0_dp
      allocate(penalty(size(design, 2), size(design, 2)))
      penalty = 0.0_dp
      do component = 1, size(lambda)
         first = 2 + (component - 1)*basis_count
         last = first + basis_count - 1
         penalty(first:last, first:last) = lambda(component)*base_penalty
      end do
      cross_product = matmul(transpose(design), design)
      call invert_matrix(cross_product + penalty, inverse, info)
      if (info /= 0) return
      influence = matmul(inverse, cross_product)
      do component = 1, size(lambda)
         first = 2 + (component - 1)*basis_count
         last = first + basis_count - 1
         do index = first, last
            component_df(component) = component_df(component) + influence(index, index)
         end do
      end do
   end subroutine aar_component_degrees

   pure subroutine llar_embedded_prediction(series, order, delay, steps, &
      embedding_count, query, epsilon, prediction, neighbors, info)
      !! Fit and evaluate one local-linear model in delay-embedding space.
      real(dp), intent(in) :: series(:) !! Complete univariate series.
      integer, intent(in) :: order !! Delay-embedding dimension.
      integer, intent(in) :: delay !! Positive embedding delay.
      integer, intent(in) :: steps !! Forecast step and Theiler window.
      integer, intent(in) :: embedding_count !! Number of training embeddings.
      integer, intent(in) :: query !! One-based query embedding index.
      real(dp), intent(in) :: epsilon !! Neighborhood radius.
      real(dp), intent(out) :: prediction !! Local-linear prediction.
      integer, intent(out) :: neighbors !! Number of eligible neighbors.
      integer, intent(out) :: info !! Zero on success.
      real(dp), allocatable :: design(:, :), response(:), query_vector(:)
      real(dp) :: distance
      integer :: candidate, component, row

      info = 1
      prediction = 0.0_dp
      neighbors = 0
      if (embedding_count < 1 .or. query < 1 .or. epsilon <= 0.0_dp .or. &
         query + (order - 1)*delay > size(series)) return
      do candidate = 1, embedding_count
         if (abs(candidate - query) <= steps) cycle
         distance = llar_embedding_distance(series, order, delay, query, candidate)
         if (distance < epsilon) neighbors = neighbors + 1
      end do
      if (neighbors <= 2*(order + 1)) return
      allocate(design(neighbors, order + 1), response(neighbors))
      allocate(query_vector(order + 1))
      row = 0
      do candidate = 1, embedding_count
         if (abs(candidate - query) <= steps) cycle
         distance = llar_embedding_distance(series, order, delay, query, candidate)
         if (distance >= epsilon) cycle
         row = row + 1
         design(row, 1) = 1.0_dp
         do component = 1, order
            design(row, component + 1) = &
               series(candidate + (order - component)*delay)
         end do
         response(row) = series(candidate + (order - 1)*delay + steps)
      end do
      query_vector(1) = 1.0_dp
      do component = 1, order
         query_vector(component + 1) = series(query + (order - component)*delay)
      end do
      call llar_regression_prediction(design, response, query_vector, &
         prediction, info)
   end subroutine llar_embedded_prediction

   pure subroutine llar_next_prediction(series, order, delay, steps, epsilon, &
      prediction, neighbors, info)
      !! Forecast from the terminal delay embedding of a series.
      real(dp), intent(in) :: series(:) !! Available series history.
      integer, intent(in) :: order !! Delay-embedding dimension.
      integer, intent(in) :: delay !! Positive embedding delay.
      integer, intent(in) :: steps !! Forecast step and Theiler window.
      real(dp), intent(in) :: epsilon !! Neighborhood radius.
      real(dp), intent(out) :: prediction !! Local-linear forecast.
      integer, intent(out) :: neighbors !! Number of eligible neighbors.
      integer, intent(out) :: info !! Zero on success.
      integer :: embedding_count, query

      embedding_count = size(series) - (order - 1)*delay - steps
      query = embedding_count + steps
      call llar_embedded_prediction(series, order, delay, steps, embedding_count, &
         query, epsilon, prediction, neighbors, info)
   end subroutine llar_next_prediction

   pure real(dp) function llar_embedding_distance(series, order, delay, left, &
      right) result(distance)
      !! Return Euclidean distance between two delay embeddings.
      real(dp), intent(in) :: series(:) !! Complete univariate series.
      integer, intent(in) :: order !! Delay-embedding dimension.
      integer, intent(in) :: delay !! Positive embedding delay.
      integer, intent(in) :: left !! First embedding start index.
      integer, intent(in) :: right !! Second embedding start index.
      integer :: component

      distance = 0.0_dp
      do component = 0, order - 1
         distance = distance + (series(left + component*delay) - &
            series(right + component*delay))**2
      end do
      distance = sqrt(distance)
   end function llar_embedding_distance

   pure subroutine llar_regression_prediction(design, response, query, &
      prediction, info)
      !! Fit a local regression and evaluate it at one query vector.
      real(dp), intent(in) :: design(:, :) !! Local intercept and lag regressors.
      real(dp), intent(in) :: response(:) !! Local future observations.
      real(dp), intent(in) :: query(:) !! Query intercept and lag regressors.
      real(dp), intent(out) :: prediction !! Evaluated local-linear prediction.
      integer, intent(out) :: info !! Zero on success.
      real(dp), allocatable :: coefficient(:), standard_error(:), residual(:)
      real(dp), allocatable :: scaled(:, :), scaled_query(:)
      real(dp), allocatable :: cross_product(:, :), inverse(:, :)
      real(dp) :: rss, norm
      integer :: column, status

      call ols_fit(design, response, coefficient, standard_error, residual, &
         rss, status)
      if (status == 0 .and. all(ieee_is_finite(coefficient))) then
         prediction = dot_product(query, coefficient)
         info = 0
         return
      end if
      allocate(scaled(size(design, 1), size(design, 2)))
      allocate(scaled_query(size(query)))
      do column = 1, size(design, 2)
         norm = sqrt(sum(design(:, column)*design(:, column)))
         if (norm > tiny(1.0_dp)) then
            scaled(:, column) = design(:, column)/norm
            scaled_query(column) = query(column)/norm
         else
            scaled(:, column) = 0.0_dp
            scaled_query(column) = 0.0_dp
         end if
      end do
      cross_product = matmul(transpose(scaled), scaled)
      do column = 1, size(cross_product, 1)
         cross_product(column, column) = cross_product(column, column) + 1.0e-10_dp
      end do
      call invert_matrix(cross_product, inverse, status)
      if (status /= 0) then
         info = 1
         prediction = 0.0_dp
         return
      end if
      coefficient = matmul(inverse, matmul(transpose(scaled), response))
      prediction = dot_product(scaled_query, coefficient)
      info = merge(0, 1, ieee_is_finite(prediction))
   end subroutine llar_regression_prediction

   pure integer function llar_minimum_finite_index(values) result(index)
      !! Return the index of the smallest finite value, or zero if absent.
      real(dp), intent(in) :: values(:) !! Candidate objective values.
      real(dp) :: best
      integer :: candidate

      index = 0
      best = huge(1.0_dp)
      do candidate = 1, size(values)
         if (ieee_is_finite(values(candidate)) .and. values(candidate) < best) then
            best = values(candidate)
            index = candidate
         end if
      end do
   end function llar_minimum_finite_index

   pure subroutine lstar_inference(response, low_design, high_design, &
      transition_data, coefficient, gamma, threshold, residual_variance, &
      covariance, standard_error, confidence_lower, confidence_upper)
      !! Compute numerical-Hessian covariance and approximate LSTAR intervals.
      real(dp), intent(in) :: response(:) !! Aligned response vector.
      real(dp), intent(in) :: low_design(:, :) !! Baseline regression matrix.
      real(dp), intent(in) :: high_design(:, :) !! Transition regression matrix.
      real(dp), intent(in) :: transition_data(:) !! Transition-variable values.
      real(dp), intent(in) :: coefficient(:) !! Profiled linear coefficients.
      real(dp), intent(in) :: gamma !! Estimated logistic slope.
      real(dp), intent(in) :: threshold !! Estimated logistic midpoint.
      real(dp), intent(in) :: residual_variance !! Estimated innovation variance.
      real(dp), allocatable, intent(out) :: covariance(:, :) !! Parameter covariance matrix.
      real(dp), allocatable, intent(out) :: standard_error(:) !! Parameter standard errors.
      real(dp), allocatable, intent(out) :: confidence_lower(:) !! Approximate 95 percent lower limits.
      real(dp), allocatable, intent(out) :: confidence_upper(:) !! Approximate 95 percent upper limits.
      real(dp), allocatable :: parameter(:), step(:), hessian(:, :), inverse(:, :)
      real(dp), allocatable :: pp(:), pm(:), mp(:), mm(:)
      real(dp) :: center, plus, minus, nan
      integer :: count, left, right, info

      count = size(coefficient) + 2
      allocate(parameter(count), step(count), hessian(count, count))
      parameter(:size(coefficient)) = coefficient
      parameter(count - 1) = gamma
      parameter(count) = threshold
      step = 1.0e-4_dp*(abs(parameter) + 1.0_dp)
      center = lstar_full_rss(parameter, response, low_design, high_design, &
         transition_data)
      hessian = 0.0_dp
      do left = 1, count
         pp = parameter
         pm = parameter
         pp(left) = pp(left) + step(left)
         pm(left) = pm(left) - step(left)
         plus = lstar_full_rss(pp, response, low_design, high_design, transition_data)
         minus = lstar_full_rss(pm, response, low_design, high_design, transition_data)
         hessian(left, left) = (plus - 2.0_dp*center + minus)/(step(left)**2)
         do right = left + 1, count
            pp = parameter
            pm = parameter
            mp = parameter
            mm = parameter
            pp(left) = pp(left) + step(left)
            pp(right) = pp(right) + step(right)
            pm(left) = pm(left) + step(left)
            pm(right) = pm(right) - step(right)
            mp(left) = mp(left) - step(left)
            mp(right) = mp(right) + step(right)
            mm(left) = mm(left) - step(left)
            mm(right) = mm(right) - step(right)
            hessian(left, right) = (lstar_full_rss(pp, response, low_design, &
               high_design, transition_data) - lstar_full_rss(pm, response, &
               low_design, high_design, transition_data) - lstar_full_rss(mp, &
               response, low_design, high_design, transition_data) + &
               lstar_full_rss(mm, response, low_design, high_design, &
               transition_data))/(4.0_dp*step(left)*step(right))
            hessian(right, left) = hessian(left, right)
         end do
      end do
      call invert_matrix(hessian, inverse, info)
      allocate(covariance(count, count), standard_error(count))
      allocate(confidence_lower(count), confidence_upper(count))
      if (info /= 0 .or. .not. all(ieee_is_finite(inverse))) then
         nan = ieee_value(0.0_dp, ieee_quiet_nan)
         covariance = nan
         standard_error = nan
         confidence_lower = nan
         confidence_upper = nan
         return
      end if
      covariance = 2.0_dp*residual_variance*inverse
      do left = 1, count
         standard_error(left) = sqrt(max(0.0_dp, covariance(left, left)))
      end do
      confidence_lower = parameter - 1.959963984540054_dp*standard_error
      confidence_upper = parameter + 1.959963984540054_dp*standard_error
   end subroutine lstar_inference

   pure real(dp) function lstar_full_rss(parameter, response, low_design, &
      high_design, transition_data) result(rss)
      !! Evaluate LSTAR RSS for a complete linear and nonlinear parameter vector.
      real(dp), intent(in) :: parameter(:) !! Linear coefficients, gamma, and threshold.
      real(dp), intent(in) :: response(:) !! Aligned response vector.
      real(dp), intent(in) :: low_design(:, :) !! Baseline regression matrix.
      real(dp), intent(in) :: high_design(:, :) !! Transition regression matrix.
      real(dp), intent(in) :: transition_data(:) !! Transition-variable values.
      real(dp), allocatable :: weight(:), fitted(:)
      integer :: low_columns, high_columns, count

      low_columns = size(low_design, 2)
      high_columns = size(high_design, 2)
      count = size(parameter)
      if (count /= low_columns + high_columns + 2 .or. &
         parameter(count - 1) <= 0.0_dp) then
         rss = huge(1.0_dp)
         return
      end if
      weight = tsdyn_logistic_transition(transition_data, parameter(count - 1), &
         parameter(count))
      fitted = matmul(low_design, parameter(:low_columns)) + weight* &
         matmul(high_design, parameter(low_columns + 1:low_columns + high_columns))
      rss = sum((response - fitted)**2)
   end function lstar_full_rss

   pure subroutine star_profile_fit(response, base, transition_data, nonlinear, &
      coefficient, fitted, residual, weight, rss, info)
      !! Profile multi-regime STAR linear coefficients for fixed transitions.
      real(dp), intent(in) :: response(:) !! Aligned response vector.
      real(dp), intent(in) :: base(:, :) !! Common linear regressors.
      real(dp), intent(in) :: transition_data(:) !! Transition-variable values.
      real(dp), intent(in) :: nonlinear(:) !! Alternating log-gamma and threshold values.
      real(dp), allocatable, intent(out) :: coefficient(:) !! Profiled linear coefficients.
      real(dp), allocatable, intent(out) :: fitted(:) !! Profiled fitted values.
      real(dp), allocatable, intent(out) :: residual(:) !! Profiled residuals.
      real(dp), allocatable, intent(out) :: weight(:, :) !! Transition weights by component.
      real(dp), intent(out) :: rss !! Residual sum of squares.
      integer, intent(out) :: info !! Zero on success.
      real(dp), allocatable :: design(:, :), standard_error(:)
      integer :: components, columns, component, row, offset

      components = size(nonlinear)/2
      columns = size(base, 2)
      info = 1
      rss = huge(1.0_dp)
      if (components < 1 .or. size(nonlinear) /= 2*components .or. &
         any(nonlinear(1::2) > log(1.0e4_dp)) .or. &
         any(nonlinear(1::2) < log(1.0e-4_dp))) return
      if (components > 1) then
         do component = 2, components
            if (nonlinear(2*component) <= nonlinear(2*component - 2)) return
         end do
      end if
      allocate(weight(size(response), components))
      allocate(design(size(response), columns*(components + 1)))
      design(:, :columns) = base
      do component = 1, components
         weight(:, component) = tsdyn_logistic_transition(transition_data, &
            exp(nonlinear(2*component - 1)), nonlinear(2*component))
         offset = component*columns
         do row = 1, size(response)
            design(row, offset + 1:offset + columns) = &
               weight(row, component)*base(row, :)
         end do
      end do
      call ols_fit(design, response, coefficient, standard_error, residual, rss, info)
      if (info == 0) fitted = response - residual
   end subroutine star_profile_fit

   pure function star_auxiliary_test(residual, null_design, alternative_design) &
      result(out)
      !! Compare nested auxiliary regressions for an omitted STAR component.
      real(dp), intent(in) :: residual(:) !! Fitted-model residuals.
      real(dp), intent(in) :: null_design(:, :) !! Null gradient regressors.
      real(dp), intent(in) :: alternative_design(:, :) !! Additional nonlinear regressors.
      type(tsdyn_regime_test_t) :: out
      real(dp), allocatable :: full_design(:, :)
      real(dp) :: null_rss, full_rss
      integer :: info, null_columns, full_columns

      if (size(null_design, 1) /= size(residual) .or. &
         size(alternative_design, 1) /= size(residual)) then
         out%info = 1
         return
      end if
      null_columns = size(null_design, 2)
      full_columns = null_columns + size(alternative_design, 2)
      call star_stable_regression_rss(null_design, residual, null_rss, info)
      if (info /= 0) then
         null_rss = sum(residual*residual)
         null_columns = 0
         full_columns = size(alternative_design, 2)
         full_design = alternative_design
      else
         allocate(full_design(size(residual), full_columns))
         full_design(:, :null_columns) = null_design
         full_design(:, null_columns + 1:) = alternative_design
      end if
      call star_stable_regression_rss(full_design, residual, full_rss, info)
      out%numerator_df = size(alternative_design, 2)
      out%denominator_df = size(residual) - full_columns
      if (info /= 0 .or. full_rss <= 0.0_dp .or. out%denominator_df <= 0) then
         out%info = 2
         return
      end if
      out%statistic = max(0.0_dp, null_rss - full_rss)* &
         real(out%denominator_df, dp)/(real(out%numerator_df, dp)*full_rss)
      out%p_value = star_f_upper_probability(out%statistic, out%numerator_df, &
         out%denominator_df)
   end function star_auxiliary_test

   pure subroutine star_stable_regression_rss(design, response, rss, info)
      !! Compute auxiliary-regression RSS after column scaling and light regularization.
      real(dp), intent(in) :: design(:, :) !! Auxiliary regression matrix.
      real(dp), intent(in) :: response(:) !! Auxiliary response vector.
      real(dp), intent(out) :: rss !! Residual sum of squares.
      integer, intent(out) :: info !! Zero on success.
      real(dp), allocatable :: scaled(:, :), cross_product(:, :), inverse(:, :)
      real(dp), allocatable :: coefficient(:), residual(:)
      real(dp) :: norm
      integer :: column, status

      info = 1
      rss = huge(1.0_dp)
      if (size(design, 1) /= size(response) .or. size(design, 2) < 1) return
      allocate(scaled(size(design, 1), size(design, 2)))
      do column = 1, size(design, 2)
         norm = sqrt(sum(design(:, column)*design(:, column)))
         if (norm > tiny(1.0_dp)) then
            scaled(:, column) = design(:, column)/norm
         else
            scaled(:, column) = 0.0_dp
         end if
      end do
      cross_product = matmul(transpose(scaled), scaled)
      do column = 1, size(cross_product, 1)
         cross_product(column, column) = cross_product(column, column) + 1.0e-10_dp
      end do
      call invert_matrix(cross_product, inverse, status)
      if (status /= 0) return
      coefficient = matmul(inverse, matmul(transpose(scaled), response))
      residual = response - matmul(scaled, coefficient)
      rss = sum(residual*residual)
      if (.not. ieee_is_finite(rss)) return
      info = 0
   end subroutine star_stable_regression_rss

   pure real(dp) function star_f_upper_probability(statistic, numerator_df, &
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
   end function star_f_upper_probability

   pure subroutine lstar_regression_data(series, low_order, high_order, &
      threshold_lag, include_constant, include_trend, threshold_variable, &
      response, low_design, high_design, transition_data)
      !! Build aligned LSTAR responses, regressors, and transition data.
      real(dp), intent(in) :: series(:) !! Complete univariate series.
      integer, intent(in) :: low_order !! Baseline AR order.
      integer, intent(in) :: high_order !! Transition AR order.
      integer, intent(in) :: threshold_lag !! Internal transition lag.
      logical, intent(in) :: include_constant !! Whether to include intercepts.
      logical, intent(in) :: include_trend !! Whether to include trends.
      real(dp), intent(in), optional :: threshold_variable(:) !! External transition series.
      real(dp), allocatable, intent(out) :: response(:) !! Aligned response.
      real(dp), allocatable, intent(out) :: low_design(:, :) !! Baseline regressors.
      real(dp), allocatable, intent(out) :: high_design(:, :) !! Transition regressors.
      real(dp), allocatable, intent(out) :: transition_data(:) !! Aligned transition variable.
      integer :: start, rows, constants, row, time, column, lag

      start = max(threshold_lag, max(low_order, high_order)) + 1
      rows = size(series) - start + 1
      constants = merge(1, 0, include_constant) + merge(1, 0, include_trend)
      allocate(response(rows), transition_data(rows))
      allocate(low_design(rows, constants + low_order))
      allocate(high_design(rows, constants + high_order))
      low_design = 0.0_dp
      high_design = 0.0_dp
      do row = 1, rows
         time = start + row - 1
         response(row) = series(time)
         if (present(threshold_variable)) then
            transition_data(row) = threshold_variable(time)
         else
            transition_data(row) = series(time - threshold_lag)
         end if
         column = 0
         if (include_constant) then
            column = column + 1
            low_design(row, column) = 1.0_dp
            high_design(row, column) = 1.0_dp
         end if
         if (include_trend) then
            column = column + 1
            low_design(row, column) = real(time, dp)
            high_design(row, column) = real(time, dp)
         end if
         do lag = 1, low_order
            low_design(row, column + lag) = series(time - lag)
         end do
         do lag = 1, high_order
            high_design(row, column + lag) = series(time - lag)
         end do
      end do
   end subroutine lstar_regression_data

   pure subroutine lstar_profile_fit(response, low_design, high_design, &
      transition_data, gamma, threshold, coefficient, fitted, residual, &
      weight, rss, info)
      !! Profile LSTAR linear coefficients for fixed transition parameters.
      real(dp), intent(in) :: response(:) !! Aligned response vector.
      real(dp), intent(in) :: low_design(:, :) !! Baseline regression matrix.
      real(dp), intent(in) :: high_design(:, :) !! Transition regression matrix.
      real(dp), intent(in) :: transition_data(:) !! Transition-variable values.
      real(dp), intent(in) :: gamma !! Positive logistic slope.
      real(dp), intent(in) :: threshold !! Logistic midpoint.
      real(dp), allocatable, intent(out) :: coefficient(:) !! Profiled coefficients.
      real(dp), allocatable, intent(out) :: fitted(:) !! Profiled fitted values.
      real(dp), allocatable, intent(out) :: residual(:) !! Profiled residuals.
      real(dp), allocatable, intent(out) :: weight(:) !! Logistic transition weights.
      real(dp), intent(out) :: rss !! Residual sum of squares.
      integer, intent(out) :: info !! Zero on success.
      real(dp), allocatable :: design(:, :), standard_error(:)
      integer :: low_columns, row

      low_columns = size(low_design, 2)
      allocate(weight(size(response)))
      weight = tsdyn_logistic_transition(transition_data, gamma, threshold)
      allocate(design(size(response), low_columns + size(high_design, 2)))
      design(:, :low_columns) = low_design
      do row = 1, size(response)
         design(row, low_columns + 1:) = weight(row)*high_design(row, :)
      end do
      call ols_fit(design, response, coefficient, standard_error, residual, rss, info)
      if (info /= 0) return
      fitted = response - residual
   end subroutine lstar_profile_fit

   pure real(dp) function delta_correlation_integral(series, &
      embedding_dimension, delay, epsilon_value) result(value)
      !! Estimate a correlation integral with the maximum embedding norm.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      integer, intent(in) :: embedding_dimension !! Number of embedded coordinates.
      integer, intent(in) :: delay !! Delay between embedded coordinates.
      real(dp), intent(in) :: epsilon_value !! Maximum-norm distance threshold.
      integer :: rows, first, second, coordinate, first_time, second_time
      integer :: close_pairs, eligible_pairs
      real(dp) :: distance

      value = 0.0_dp
      rows = size(series) - (embedding_dimension - 1)*delay
      if (embedding_dimension < 1 .or. delay < 1 .or. rows < 3 .or. &
         epsilon_value <= 0.0_dp) return
      close_pairs = 0
      eligible_pairs = 0
      do first = 1, rows - 2
         do second = first + 2, rows
            eligible_pairs = eligible_pairs + 1
            distance = 0.0_dp
            do coordinate = 0, embedding_dimension - 1
               first_time = first + coordinate*delay
               second_time = second + coordinate*delay
               distance = max(distance, abs(series(first_time) - &
                  series(second_time)))
            end do
            if (distance <= epsilon_value) close_pairs = close_pairs + 1
         end do
      end do
      if (eligible_pairs > 0) value = real(close_pairs, dp)/ &
         real(eligible_pairs, dp)
   end function delta_correlation_integral

   pure subroutine summarize_forecast_distribution(paths, level, mean, &
      standard_error, lower, upper)
      !! Summarize simulated paths by horizon using empirical type-7 quantiles.
      real(dp), intent(in) :: paths(:, :) !! Horizon-by-simulation forecast paths.
      real(dp), intent(in) :: level !! Central interval coverage.
      real(dp), allocatable, intent(out) :: mean(:) !! Pointwise path means.
      real(dp), allocatable, intent(out) :: standard_error(:) !! Pointwise sample standard deviations.
      real(dp), allocatable, intent(out) :: lower(:) !! Lower empirical interval limits.
      real(dp), allocatable, intent(out) :: upper(:) !! Upper empirical interval limits.
      real(dp), allocatable :: ordered(:)
      real(dp) :: tail
      integer :: step, simulations

      simulations = size(paths, 2)
      tail = 0.5_dp*(1.0_dp - level)
      allocate(mean(size(paths, 1)), standard_error(size(paths, 1)))
      allocate(lower(size(paths, 1)), upper(size(paths, 1)))
      do step = 1, size(paths, 1)
         mean(step) = sum(paths(step, :))/real(simulations, dp)
         if (simulations > 1) then
            standard_error(step) = sqrt(sum((paths(step, :) - &
               mean(step))**2)/real(simulations - 1, dp))
         else
            standard_error(step) = 0.0_dp
         end if
         ordered = sorted(paths(step, :))
         lower(step) = quantile(ordered, tail)
         upper(step) = quantile(ordered, 1.0_dp - tail)
      end do
   end subroutine summarize_forecast_distribution

   pure subroutine initialize_rolling_result(training, new_data, horizons, &
      refit_every, out)
      !! Allocate and align a rolling-origin result table.
      real(dp), intent(in) :: training(:) !! Original estimation sample.
      real(dp), intent(in) :: new_data(:) !! Sequential realized observations.
      integer, intent(in) :: horizons(:) !! Positive forecast horizons.
      integer, intent(in) :: refit_every !! Refit interval; zero disables refitting.
      type(rolling_forecast_result_t), intent(out) :: out !! Initialized rolling result.
      real(dp) :: nan
      integer :: target_index, horizon_index

      if (size(training) < 2 .or. size(new_data) < 1 .or. &
         size(horizons) < 1 .or. any(horizons < 1) .or. refit_every < 0 .or. &
         .not. all(ieee_is_finite(training)) .or. &
         .not. all(ieee_is_finite(new_data))) then
         out%info = 1
         return
      end if
      nan = ieee_value(0.0_dp, ieee_quiet_nan)
      out%training_size = size(training)
      out%horizon = horizons
      allocate(out%target(size(new_data)))
      allocate(out%origin(size(new_data), size(horizons)))
      allocate(out%forecast(size(new_data), size(horizons)))
      allocate(out%actual(size(new_data), size(horizons)))
      allocate(out%error(size(new_data), size(horizons)))
      allocate(out%valid(size(new_data), size(horizons)))
      allocate(out%refitted(size(new_data), size(horizons)))
      out%forecast = nan
      out%error = nan
      out%valid = .false.
      out%refitted = .false.
      do target_index = 1, size(new_data)
         out%target(target_index) = size(training) + target_index
         out%actual(target_index, :) = new_data(target_index)
         do horizon_index = 1, size(horizons)
            out%origin(target_index, horizon_index) = &
               out%target(target_index) - horizons(horizon_index)
         end do
      end do
   end subroutine initialize_rolling_result

   pure elemental logical function rolling_refit_due(origin, training_size, &
      refit_every) result(due)
      !! Test whether a rolling origin is on the requested refit schedule.
      integer, intent(in) :: origin !! Absolute forecast origin.
      integer, intent(in) :: training_size !! Original estimation-sample size.
      integer, intent(in) :: refit_every !! Positive refit interval or zero.

      due = refit_every > 0 .and. origin > training_size
      if (due) due = modulo(origin - training_size, refit_every) == 0
   end function rolling_refit_due

   pure subroutine record_rolling_forecast(out, target_index, horizon_index, &
      forecast_value, actual_value)
      !! Store one valid rolling forecast and its signed error.
      type(rolling_forecast_result_t), intent(inout) :: out !! Rolling result table.
      integer, intent(in) :: target_index !! Target row to populate.
      integer, intent(in) :: horizon_index !! Horizon column to populate.
      real(dp), intent(in) :: forecast_value !! Forecast at the requested target.
      real(dp), intent(in) :: actual_value !! Realized target observation.

      if (.not. ieee_is_finite(forecast_value) .or. &
         .not. ieee_is_finite(actual_value)) return
      out%forecast(target_index, horizon_index) = forecast_value
      out%actual(target_index, horizon_index) = actual_value
      out%error(target_index, horizon_index) = actual_value - forecast_value
      out%valid(target_index, horizon_index) = .true.
   end subroutine record_rolling_forecast

   pure subroutine delta_largest_covariance_eigenvalue(series, &
      embedding_dimension, delay, eigenvalue, info)
      !! Return the largest covariance eigenvalue of a delay embedding.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      integer, intent(in) :: embedding_dimension !! Number of embedded coordinates.
      integer, intent(in) :: delay !! Delay between embedded coordinates.
      real(dp), intent(out) :: eigenvalue !! Largest sample-covariance eigenvalue.
      integer, intent(out) :: info !! Zero on success.
      real(dp), allocatable :: embedding(:, :), covariance_matrix(:, :)
      real(dp), allocatable :: eigenvalues(:), eigenvectors(:, :)
      integer :: rows, coordinate, row

      eigenvalue = 0.0_dp
      info = 1
      rows = size(series) - (embedding_dimension - 1)*delay
      if (embedding_dimension < 1 .or. delay < 1 .or. rows < 2) return
      allocate(embedding(rows, embedding_dimension))
      do coordinate = 1, embedding_dimension
         do row = 1, rows
            embedding(row, coordinate) = series(row + &
               (coordinate - 1)*delay)
         end do
         embedding(:, coordinate) = embedding(:, coordinate) - &
            sum(embedding(:, coordinate))/real(rows, dp)
      end do
      covariance_matrix = matmul(transpose(embedding), embedding)/ &
         real(rows - 1, dp)
      call symmetric_eigen(covariance_matrix, eigenvalues, eigenvectors, info)
      if (info == 0) eigenvalue = eigenvalues(1)
   end subroutine delta_largest_covariance_eigenvalue

end module tsdyn_mod
