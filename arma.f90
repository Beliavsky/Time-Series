! SPDX-License-Identifier: MIT
! SPDX-FileComment: Topic facade for autoregressive and moving-average models.
module arma_mod
   !! Re-export principal AR, MA, ARMA, ARIMA, and SARIMA algorithms.
   use utils_mod, only: real_vector_t
   use arima2_mod, only: arima2_fit_t, arima2_aic_table_t, arima2_profile_t, &
      arma_polynomial_roots, durbin_levinson_coefficients, arima2_aic_table, &
      arima2_profile, arima2_fit, sample_arma_coefficients
   use astsa_mod, only: astsa_sarima_fit_t, astsa_sarima_forecast_t, &
      astsa_sarima_diagnostics_t, sarima_sim, sarima_sim_from_innovations, &
      sarima_likelihood, sarima_exact_likelihood, sarima_fit, sarima_forecast, &
      sarima_diagnostics, arma_spectrum, arma_check, arma_to_ar
   use itsmr_mod, only: itsmr_arma_model_t, itsmr_forecast_t, arma_acvf, &
      innovations_algorithm, innovations_ma_fit, hannan_rissanen_fit, &
      arma_mle_fit, arma_autofit, arma_infinite_ma, arma_infinite_ar, &
      arma_residuals, arma_forecast, burg_ar_fit
   use tfarima_mod, only: tfarima_forecast_t, tfarima_arima_forecast, &
      tfarima_arima_backcast, tfarima_autocovariance, &
      tfarima_partial_autocorrelation, tfarima_autocovariance_to_ma, &
      tfarima_combine_arima, tfarima_psi_weights, tfarima_pi_weights
   use robustarima_mod, only: robustarima_fit_t, robustarima_forecast_t, &
      robustarima_fit, robustarima_select_ar_order, robustarima_forecast
   use expar_mod, only: expar_fit_t, expar_selection_t, expar_forecast_t, &
      expar_evaluate, expar_initial_parameters, expar_fit, expar_select, &
      expar_forecast
   use exparma_mod, only: exparma_fit_t, exparma_selection_t, &
      exparma_evaluate, exparma_initial_parameters, exparma_fit, exparma_select
   use baystar_mod, only: baystar_prior_t, baystar_fit_t, baystar_simulation_t, &
      baystar_default_prior, baystar_log_likelihood, &
      baystar_simulate_from_innovations, baystar_simulate, baystar_fit
   use mixar_mod, only: mixar_model_t, mixar_filter_t, mixar_fit_t, &
      mixar_simulation_t, mixar_predictive_distribution_t, &
      mixar_exact_forecast_t, mixar_forecast_sample_t, &
      mixar_regression_model_t, mixar_regression_fit_t, &
      mixar_regression_simulation_t, &
      mixar_seasonal_model_t, mixar_seasonal_fit_t, &
      mixar_seasonal_simulation_t, &
      mixar_inference_t, mixar_bic_selection_t, &
      mixar_model, mixar_component_locations, &
      mixar_conditional_density, mixar_conditional_cdf, mixar_filter, &
      mixar_log_likelihood, mixar_simulate_from_draws, mixar_simulate, &
      mixar_is_stable, mixar_fit, mixar_general_fit, mixar_standard_density, &
      mixar_standard_cdf, mixar_exact_forecast, mixar_predictive_density, &
      mixar_predictive_cdf, mixar_predictive_quantile, &
      mixar_forecast_from_draws, mixar_forecast
   use mixar_mod, only: mixar_regression_model, mixar_regression_filter, &
      mixar_regression_log_likelihood, mixar_regression_fit, &
      mixar_regression_simulate_from_draws, mixar_regression_simulate
   use mixar_mod, only: mixar_seasonal_model, mixar_seasonal_expanded_model, &
      mixar_seasonal_filter, mixar_seasonal_log_likelihood, &
      mixar_seasonal_fit, mixar_seasonal_is_stable, &
      mixar_seasonal_simulate_from_draws, mixar_seasonal_simulate
   use mixar_mod, only: mixar_observed_inference, &
      mixar_seasonal_observed_inference, mixar_bic, mixar_select_bic
   use mixar_mod, only: mixar_bayesian_random_t, mixar_bayesian_draws_t, &
      mixar_marginal_likelihood_t, mixar_bayesian_random, &
      mixar_bayesian_sample_from_random, mixar_bayesian_sample, &
      mixar_bayesian_relabel, mixar_marginal_likelihood
   use mixar_mod, only: mixar_order_random_t, mixar_order_selection_t, &
      mixar_order_prior_flat, mixar_order_prior_ratio, &
      mixar_order_prior_poisson, mixar_order_birth_death_probability, &
      mixar_order_random, mixar_order_select_from_random, mixar_order_select
   use mixar_mod, only: mixar_initialization_t, mixar_multistart_t, &
      mixar_initialize_from_indices, mixar_random_subsample_indices, &
      mixar_random_initialize, mixar_multistart_fit, &
      mixar_random_multistart_fit
   use mixar_mod, only: mixar_moment_t, mixar_stationary_moments_t, &
      mixar_standard_moment, mixar_standard_absolute_moment, &
      mixar_innovation_moment, mixar_conditional_moment, &
      mixar_conditional_central_moment, mixar_conditional_kurtosis, &
      mixar_conditional_excess_kurtosis, mixar_stationary_moments
   use nts_mod, only: nts_tar_simulation_t, nts_tar_model_t, &
      nts_tar_search_t, nts_tar_forecast_t, nts_msar_simulation_t, &
      nts_tar_backtest_t, nts_tvar_fit_t, nts_rcar_fit_t, &
      nts_utar_simulate_from_innovations, nts_utar_simulate, &
      nts_utar_estimate, nts_utar_threshold_search, &
      nts_utar_forecast_draws, nts_utar_forecast, &
      nts_msar_simulate_from_draws, nts_msar_simulate, &
      nts_tar_backtest_draws, nts_tar_backtest, &
      nts_tvar_filter_smooth, nts_tvar_fit, nts_rcar_fit
   use tsdyn_mod, only: tsdyn_lstar_simulation_t, tsdyn_lstar_model_t, &
      tsdyn_lstar_forecast_t, tsdyn_lstar_selection_t, tsdyn_regime_test_t, &
      tsdyn_star_model_t, tsdyn_star_forecast_t, tsdyn_llar_model_t, &
      tsdyn_llar_forecast_t, tsdyn_aar_model_t, tsdyn_aar_forecast_t, &
      tsdyn_forecast_distribution_t, tsdyn_regime_diagnostics_t, &
      tsdyn_setar_variance_t, tsdyn_setar_inference_t, &
      tsdyn_regime_path_t, tsdyn_irf_t, tsdyn_girf_t, &
      tsdyn_logistic_transition, &
      tsdyn_lstar_simulate_from_innovations, tsdyn_lstar_simulate, &
      tsdyn_lstar_fit, tsdyn_lstar_forecast, tsdyn_lstar_select, &
      tsdyn_lstar_gradient, tsdyn_lstar_regime_test, tsdyn_star_fit, &
      tsdyn_star_forecast, tsdyn_star_gradient, tsdyn_star_regime_test, &
      tsdyn_star_add_regime, tsdyn_llar_fit, tsdyn_llar_forecast, &
      tsdyn_aar_fit, tsdyn_aar_component, tsdyn_aar_forecast, &
      tsdyn_setar_fit, tsdyn_setar_simulate_from_innovations, &
      tsdyn_setar_forecast_draws, tsdyn_setar_selection_t, &
      tsdyn_setar_select, tsdyn_regime_count_test_t, &
      tsdyn_setar_bootstrap_t, tsdyn_setar_bootstrap, &
      tsdyn_setar_restricted_fit, tsdyn_setar_linearity_test, &
      tsdyn_nnet_model_t, tsdyn_nnet_forecast_t, tsdyn_nnet_selection_t, &
      tsdyn_nnet_fit, tsdyn_nnet_predict, tsdyn_nnet_forecast, &
      tsdyn_nnet_select, tsdyn_setar_forecast_distribution, &
      tsdyn_lstar_forecast_distribution, tsdyn_star_forecast_distribution, &
      tsdyn_llar_forecast_distribution, tsdyn_aar_forecast_distribution, &
      tsdyn_nnet_forecast_distribution, tsdyn_ar_regime_diagnostics, &
      tsdyn_setar_regime_diagnostics, tsdyn_lstar_regime_diagnostics, &
      tsdyn_star_regime_diagnostics, tsdyn_setar_regime_irf, &
      tsdyn_setar_girf_from_innovations, tsdyn_setar_irf_bootstrap
   use tsdyn_mod, only: tsdyn_setar_residual_variance, tsdyn_setar_inference
   use tsdyn_mod, only: tsdyn_setar_regimes, tsdyn_lstar_regimes, &
      tsdyn_star_regimes
   use ugmar_mod, only: ugmar_model_t, ugmar_evaluation_t, &
      ugmar_stationary_moments_t, ugmar_fit_t, ugmar_genetic_fit_t, &
      ugmar_multistart_fit_t, ugmar_constraints_t, &
      ugmar_likelihood_profile_t, ugmar_regime_conversion_t, &
      ugmar_ar_roots_t, ugmar_model, ugmar_random_model, ugmar_smart_model, &
      ugmar_ar_roots, &
      ugmar_is_stationary, ugmar_evaluate, ugmar_log_likelihood, &
      ugmar_mixing_weights, ugmar_regime_means, &
      ugmar_model_from_regime_means, ugmar_stationary_moments, &
      ugmar_estimate, ugmar_genetic_estimate, ugmar_multistart_estimate, &
      ugmar_estimate_constrained, ugmar_profile_likelihood, &
      ugmar_convert_student_regimes
   use setartree_mod, only: setartree_category_levels_t, setartree_model_t, &
      setarforest_model_t, &
      setartree_stop_both, setartree_stop_linearity, setartree_stop_error, &
      setartree_fit, setartree_fit_categorical, setartree_fit_series, &
      setarforest_fit, setarforest_fit_categorical, setarforest_fit_series
   use tseriestarma_mod, only: tseriestarma_model_t, &
      tseriestarma_evaluation_t, tseriestarma_simulation_t, &
      tseriestarma_fit_t, tseriestarma_fit2_t, tseriestarma_derivatives_t, &
      tseriestarma_innovation_normal, tseriestarma_innovation_student, &
      tseriestarma_model, tseriestarma_evaluate, &
      tseriestarma_simulate_from_standard, tseriestarma_simulate, &
      tseriestarma_initial_parameters, tseriestarma_fit, &
      tseriestarma_residual_derivatives, &
      tseriestarma_least_squares_gradient, tseriestarma_robust_weights, &
      tseriestarma_robust_fit, tseriestarma_fit2
   implicit none
   private

   public :: arima2_fit_t, arima2_aic_table_t, arima2_profile_t
   public :: arma_polynomial_roots, durbin_levinson_coefficients
   public :: arima2_aic_table, arima2_profile, arima2_fit, sample_arma_coefficients
   public :: astsa_sarima_fit_t, astsa_sarima_forecast_t, astsa_sarima_diagnostics_t
   public :: sarima_sim, sarima_sim_from_innovations
   public :: sarima_likelihood, sarima_exact_likelihood, sarima_fit
   public :: sarima_forecast, sarima_diagnostics, arma_spectrum, arma_check, arma_to_ar
   public :: itsmr_arma_model_t, itsmr_forecast_t, arma_acvf
   public :: innovations_algorithm, innovations_ma_fit, hannan_rissanen_fit
   public :: arma_mle_fit, arma_autofit, arma_infinite_ma, arma_infinite_ar
   public :: arma_residuals, arma_forecast, burg_ar_fit
   public :: tfarima_forecast_t, tfarima_arima_forecast, tfarima_arima_backcast
   public :: tfarima_autocovariance, tfarima_partial_autocorrelation
   public :: tfarima_autocovariance_to_ma, tfarima_combine_arima
   public :: tfarima_psi_weights, tfarima_pi_weights
   public :: robustarima_fit_t, robustarima_forecast_t
   public :: robustarima_fit, robustarima_select_ar_order, robustarima_forecast
   public :: expar_fit_t, expar_selection_t, expar_forecast_t
   public :: expar_evaluate, expar_initial_parameters, expar_fit
   public :: expar_select, expar_forecast
   public :: exparma_fit_t, exparma_selection_t
   public :: exparma_evaluate, exparma_initial_parameters
   public :: exparma_fit, exparma_select
   public :: baystar_prior_t, baystar_fit_t, baystar_simulation_t
   public :: baystar_default_prior, baystar_log_likelihood
   public :: baystar_simulate_from_innovations, baystar_simulate, baystar_fit
   public :: mixar_model_t, mixar_filter_t, mixar_fit_t, mixar_simulation_t
   public :: mixar_predictive_distribution_t, mixar_exact_forecast_t
   public :: mixar_forecast_sample_t
   public :: mixar_regression_model_t, mixar_regression_fit_t
   public :: mixar_regression_simulation_t
   public :: mixar_seasonal_model_t, mixar_seasonal_fit_t
   public :: mixar_seasonal_simulation_t
   public :: mixar_inference_t, mixar_bic_selection_t
   public :: mixar_model, mixar_component_locations, mixar_conditional_density
   public :: mixar_conditional_cdf, mixar_filter, mixar_log_likelihood
   public :: mixar_simulate_from_draws, mixar_simulate, mixar_is_stable, mixar_fit
   public :: mixar_general_fit, mixar_standard_density, mixar_standard_cdf
   public :: mixar_exact_forecast, mixar_predictive_density
   public :: mixar_predictive_cdf, mixar_predictive_quantile
   public :: mixar_forecast_from_draws, mixar_forecast
   public :: mixar_regression_model, mixar_regression_filter
   public :: mixar_regression_log_likelihood, mixar_regression_fit
   public :: mixar_regression_simulate_from_draws, mixar_regression_simulate
   public :: mixar_seasonal_model, mixar_seasonal_expanded_model
   public :: mixar_seasonal_filter, mixar_seasonal_log_likelihood
   public :: mixar_seasonal_fit, mixar_seasonal_is_stable
   public :: mixar_seasonal_simulate_from_draws, mixar_seasonal_simulate
   public :: mixar_observed_inference, mixar_seasonal_observed_inference
   public :: mixar_bic, mixar_select_bic
   public :: mixar_bayesian_random_t, mixar_bayesian_draws_t
   public :: mixar_marginal_likelihood_t, mixar_bayesian_random
   public :: mixar_bayesian_sample_from_random, mixar_bayesian_sample
   public :: mixar_bayesian_relabel, mixar_marginal_likelihood
   public :: mixar_order_random_t, mixar_order_selection_t
   public :: mixar_order_prior_flat, mixar_order_prior_ratio
   public :: mixar_order_prior_poisson, mixar_order_birth_death_probability
   public :: mixar_order_random, mixar_order_select_from_random
   public :: mixar_order_select
   public :: mixar_initialization_t, mixar_multistart_t
   public :: mixar_initialize_from_indices, mixar_random_subsample_indices
   public :: mixar_random_initialize, mixar_multistart_fit
   public :: mixar_random_multistart_fit
   public :: mixar_moment_t, mixar_stationary_moments_t
   public :: mixar_standard_moment, mixar_standard_absolute_moment
   public :: mixar_innovation_moment, mixar_conditional_moment
   public :: mixar_conditional_central_moment, mixar_conditional_kurtosis
   public :: mixar_conditional_excess_kurtosis, mixar_stationary_moments
   public :: nts_tar_simulation_t, nts_tar_model_t
   public :: nts_tar_search_t, nts_tar_forecast_t
   public :: nts_msar_simulation_t
   public :: nts_tar_backtest_t, nts_tvar_fit_t, nts_rcar_fit_t
   public :: nts_utar_simulate_from_innovations, nts_utar_simulate
   public :: nts_utar_estimate, nts_utar_threshold_search
   public :: nts_utar_forecast_draws, nts_utar_forecast
   public :: nts_msar_simulate_from_draws, nts_msar_simulate
   public :: nts_tar_backtest_draws, nts_tar_backtest
   public :: nts_tvar_filter_smooth, nts_tvar_fit, nts_rcar_fit
   public :: tsdyn_lstar_simulation_t, tsdyn_lstar_model_t
   public :: tsdyn_lstar_forecast_t, tsdyn_lstar_selection_t
   public :: tsdyn_regime_test_t, tsdyn_star_model_t, tsdyn_star_forecast_t
   public :: tsdyn_llar_model_t, tsdyn_llar_forecast_t
   public :: tsdyn_aar_model_t, tsdyn_aar_forecast_t
   public :: tsdyn_forecast_distribution_t
   public :: tsdyn_regime_diagnostics_t
   public :: tsdyn_setar_variance_t, tsdyn_setar_inference_t
   public :: tsdyn_regime_path_t
   public :: tsdyn_irf_t, tsdyn_girf_t
   public :: tsdyn_logistic_transition
   public :: tsdyn_lstar_simulate_from_innovations, tsdyn_lstar_simulate
   public :: tsdyn_lstar_fit, tsdyn_lstar_forecast
   public :: tsdyn_lstar_select, tsdyn_lstar_gradient, tsdyn_lstar_regime_test
   public :: tsdyn_star_fit, tsdyn_star_forecast, tsdyn_star_gradient
   public :: tsdyn_star_regime_test, tsdyn_star_add_regime
   public :: tsdyn_llar_fit, tsdyn_llar_forecast
   public :: tsdyn_aar_fit, tsdyn_aar_component, tsdyn_aar_forecast
   public :: tsdyn_setar_fit, tsdyn_setar_simulate_from_innovations
   public :: tsdyn_setar_forecast_draws, tsdyn_setar_selection_t
   public :: tsdyn_setar_select
   public :: tsdyn_setar_bootstrap_t, tsdyn_setar_bootstrap
   public :: tsdyn_regime_count_test_t
   public :: tsdyn_setar_restricted_fit, tsdyn_setar_linearity_test
   public :: tsdyn_nnet_model_t, tsdyn_nnet_forecast_t, tsdyn_nnet_selection_t
   public :: tsdyn_nnet_fit, tsdyn_nnet_predict, tsdyn_nnet_forecast
   public :: tsdyn_nnet_select
   public :: tsdyn_setar_forecast_distribution
   public :: tsdyn_lstar_forecast_distribution, tsdyn_star_forecast_distribution
   public :: tsdyn_llar_forecast_distribution, tsdyn_aar_forecast_distribution
   public :: tsdyn_nnet_forecast_distribution
   public :: tsdyn_ar_regime_diagnostics, tsdyn_setar_regime_diagnostics
   public :: tsdyn_lstar_regime_diagnostics, tsdyn_star_regime_diagnostics
   public :: tsdyn_setar_residual_variance, tsdyn_setar_inference
   public :: tsdyn_setar_regimes, tsdyn_lstar_regimes, tsdyn_star_regimes
   public :: tsdyn_setar_regime_irf, tsdyn_setar_girf_from_innovations
   public :: tsdyn_setar_irf_bootstrap
   public :: ugmar_model_t, ugmar_evaluation_t, ugmar_stationary_moments_t
   public :: ugmar_fit_t, ugmar_genetic_fit_t, ugmar_multistart_fit_t
   public :: ugmar_constraints_t, ugmar_likelihood_profile_t
   public :: ugmar_regime_conversion_t, ugmar_ar_roots_t
   public :: ugmar_model, ugmar_is_stationary, ugmar_evaluate
   public :: ugmar_random_model, ugmar_smart_model, ugmar_ar_roots
   public :: ugmar_log_likelihood, ugmar_mixing_weights
   public :: ugmar_regime_means, ugmar_model_from_regime_means
   public :: ugmar_stationary_moments, ugmar_estimate
   public :: ugmar_genetic_estimate, ugmar_multistart_estimate
   public :: ugmar_estimate_constrained, ugmar_profile_likelihood
   public :: ugmar_convert_student_regimes
   public :: setartree_category_levels_t, setartree_model_t, setarforest_model_t
   public :: real_vector_t
   public :: setartree_stop_both, setartree_stop_linearity, setartree_stop_error
   public :: setartree_fit, setartree_fit_categorical, setartree_fit_series
   public :: setarforest_fit, setarforest_fit_categorical, setarforest_fit_series
   public :: tseriestarma_model_t, tseriestarma_evaluation_t
   public :: tseriestarma_simulation_t, tseriestarma_fit_t
   public :: tseriestarma_fit2_t
   public :: tseriestarma_derivatives_t
   public :: tseriestarma_innovation_normal, tseriestarma_innovation_student
   public :: tseriestarma_model, tseriestarma_evaluate
   public :: tseriestarma_simulate_from_standard, tseriestarma_simulate
   public :: tseriestarma_initial_parameters, tseriestarma_fit
   public :: tseriestarma_residual_derivatives
   public :: tseriestarma_least_squares_gradient
   public :: tseriestarma_robust_weights, tseriestarma_robust_fit
   public :: tseriestarma_fit2
end module arma_mod
