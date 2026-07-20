! SPDX-License-Identifier: MIT
! SPDX-FileComment: Topic facade for Bayesian time-series algorithms.
module bayesian_time_series_mod
   !! Re-export Bayesian state-space, regression, and simulation algorithms.
   use bayesforecast_mod, only: bf_mcmc_t, bf_information_criteria_t, &
      bf_log_prior, bf_log_posterior, bf_metropolis_sample, &
      bf_posterior_interval, bf_information_criteria, bf_bayes_factor
   use bsts_mod, only: bsts_mcmc_t, bsts_prediction_t, bsts_local_level, &
      bsts_local_linear_trend, bsts_seasonal, bsts_spike_slab, &
      bsts_dynamic_regression, bsts_predict, bsts_mbsts
   use bssm_mod, only: bssm_mcmc_t, bssm_da_mcmc_t, bssm_approximate_mcmc, &
      bssm_nonlinear_approximate_mcmc, bssm_nonlinear_pmmh, &
      bssm_sde_pmmh, bssm_importance_post_correction, bssm_mcmc_diagnostics
   use astsa_mod, only: astsa_ar_mcmc_t, astsa_sv_mcmc_t, ar_mcmc, &
      ar_mcmc_draws, sv_mcmc, astsa_ffbs, astsa_ffbs_draws
   use baystar_mod, only: baystar_prior_t, baystar_coefficient_posterior_t, &
      baystar_variance_posterior_t, baystar_summary_t, baystar_fit_t, &
      baystar_simulation_t, baystar_default_prior, baystar_log_likelihood, &
      baystar_coefficient_posterior, baystar_variance_posterior, &
      baystar_delay_probabilities, baystar_summary, &
      baystar_simulate_from_innovations, baystar_simulate, baystar_fit
   use mts_mod, only: mts_bvar_fit_t, mts_bvar_prior_t, mts_bvar_fit, &
      mts_minnesota_prior
   use bvartools_mod, only: bvartools_minnesota_prior_t, &
      bvartools_loglik_normal_t, bvartools_kalman_dk_t, &
      bvartools_inclusion_prior_t, bvartools_ssvs_result_t, &
      bvartools_ssvs_prior_t, &
      bvartools_normal_posterior_t, bvartools_gamma_posterior_t, &
      bvartools_covariance_data_t, bvartools_bvs_result_t, &
      bvartools_bvar_draws_t, bvartools_cointegration_draw_t, &
      bvartools_bvar_prior_bundle_t, bvartools_bvar_fit_t, &
      bvartools_bvec_draws_t, bvartools_bvec_prior_bundle_t, &
      bvartools_bvec_fit_t, &
      bvartools_tvp_bvar_draws_t, &
      bvartools_tvp_bvar_prior_bundle_t, bvartools_tvp_bvar_fit_t, &
      bvartools_tvp_covariance_draws_t, bvartools_joint_tvp_bvar_draws_t, &
      bvartools_structural_tvp_bvar_prior_t, &
      bvartools_structural_tvp_bvar_fit_t, &
      bvartools_stochastic_volatility_t, &
      bvartools_tvp_bvec_draws_t, &
      bvartools_tvp_bvec_prior_bundle_t, bvartools_tvp_bvec_fit_t, &
      bvartools_dfm_draws_t, bvartools_dfm_prior_t, &
      bvartools_dfm_grid_draws_t, &
      bvartools_model_likelihood_data_t, bvartools_model_comparison_t, &
      bvartools_model_comparison_set_t, &
      bvartools_predictive_t, bvartools_irf_t, bvartools_fevd_t, &
      bvartools_level_var_draws_t, bvartools_tvp_level_var_draws_t, &
      bvartools_var_data_t, bvartools_vecm_data_t, bvartools_dfm_data_t, &
      bvartools_loglik_normal, bvartools_kalman_dk, &
      bvartools_minnesota_prior, bvartools_inclusion_prior, bvartools_ssvs, &
      bvartools_ssvs_prior, &
      bvartools_normal_posterior, bvartools_normal_draw, &
      bvartools_sur_normal_posterior, bvartools_sur_normal_draw, &
      bvartools_measurement_variance_posterior, &
      bvartools_state_variance_posterior, bvartools_gamma_precision_draw, &
      bvartools_covar_prepare_data, bvartools_covar_vector_to_matrix, &
      bvartools_sur_const_to_tvp, bvartools_covar_const_posterior, &
      bvartools_covar_tvp_posterior, bvartools_bvs, &
      bvartools_random_walk_posterior, bvartools_bvar_gibbs, &
      bvartools_structural_bvar_gibbs, &
      bvartools_prepare_bvar_prior, bvartools_fit_bvar, &
      bvartools_fit_structural_bvar, &
      bvartools_cointegration_draw, bvartools_cointegration_sur_draw, &
      bvartools_bvec_gibbs, bvartools_prepare_bvec_prior, bvartools_fit_bvec, &
      bvartools_fit_structural_bvec, &
      bvartools_initial_state_posterior, bvartools_tvp_bvar_gibbs, &
      bvartools_prepare_tvp_bvar_prior, bvartools_fit_tvp_bvar, &
      bvartools_tvp_covariance_gibbs, bvartools_joint_tvp_bvar_gibbs, &
      bvartools_prepare_structural_tvp_bvar_prior, &
      bvartools_fit_structural_tvp_bvar
   use bvartools_mod, only: bvartools_fit_structural_tvp_bvec
   use bvartools_mod, only: bvartools_tvp_bvs, bvartools_tvp_covariance_bvs
   use bvartools_mod, only: bvartools_structural_impacts
   use bvartools_mod, only: bvartools_stochastic_volatility_draw, &
      bvartools_stochastic_volatility, bvartools_tvp_bvec_gibbs, &
      bvartools_prepare_tvp_bvec_prior, bvartools_fit_tvp_bvec, &
      bvartools_stochastic_volatility_ocsn_draw, &
      bvartools_stochastic_volatility_ocsn, &
      bvartools_dfm_factor_posterior, bvartools_dfm_gibbs, &
      bvartools_dfm_prior, bvartools_dfm_grid_gibbs
   use bvartools_mod, only: bvartools_model_comparison, bvartools_compare_models
   use bvartools_mod, only: bvartools_bvar_predictive_from_random, &
      bvartools_bvar_predictive, bvartools_tvp_bvar_predictive_from_random, &
      bvartools_tvp_bvar_predictive
   use bvartools_mod, only: bvartools_bvar_irf, bvartools_tvp_bvar_irf, &
      bvartools_bvar_fevd, bvartools_tvp_bvar_fevd
   use bvartools_mod, only: bvartools_vecm_level_ar, &
      bvartools_vecm_level_exogenous, bvartools_bvec_to_level_var, &
      bvartools_tvp_bvec_to_level_var, bvartools_reconstruct_levels
   use bvartools_mod, only: bvartools_prepare_var, bvartools_prepare_vecm, &
      bvartools_prepare_dfm
   use bvar_mod, only: bvar_dummy_observations_t, bvar_evidence_t, &
      bvar_hyper_evaluation_t, bvar_metropolis_result_t, &
      bvar_minnesota_variance, bvar_soc_dummy, bvar_sur_dummy, &
      bvar_conjugate_evidence, bvar_gamma_log_density, &
      bvar_inverse_gamma_log_density, bvar_hyper_log_posterior, &
      bvar_hierarchical_metropolis_from_random, bvar_hierarchical_metropolis
   use bvars_mod, only: bvars_data_t, bvars_prior_t, bvars_draws_t, &
      bvars_common_sv_t, &
      bvars_predictive_t, bvars_fitted_t, bvars_prepare, bvars_default_prior, &
      bvars_matrix_normal_from_standard, bvars_matrix_normal, &
      bvars_conjugate_draws, bvars_student_log_kernel, &
      bvars_student_t_draws, bvars_student_scale_forecast, &
      bvars_sv_auxiliary_mixture, &
      bvars_centered_sv_update, bvars_centered_sv_draws, &
      bvars_noncentered_sv_update, bvars_noncentered_sv_draws, &
      bvars_student_sv_draws, &
      bvars_common_variance_forecast, &
      bvars_shocks, bvars_fitted_from_random, &
      bvars_fitted, bvars_forecast_from_random, bvars_forecast, &
      bvars_conditional_forecast_from_random, bvars_conditional_forecast, &
      bvars_fevd
   use bvarsv_mod, only: bvarsv_data_t, bvarsv_state_draw_t, &
      bvarsv_volatility_draw_t, &
      bvarsv_contemporaneous_draw_t, &
      bvarsv_covariance_draw_t, &
      bvarsv_prior_t, &
      bvarsv_draws_t, &
      bvarsv_predictive_t, &
      bvarsv_irf_t, &
      bvarsv_simulation_t, &
      bvarsv_predictive_draws_t, &
      bvarsv_parameter_path_t, &
      bvarsv_carter_kohn, bvarsv_contemporaneous_matrix, &
      bvarsv_covariance_from_state, bvarsv_log_volatility_update, &
      bvarsv_contemporaneous_update, bvarsv_random_walk_covariance, &
      bvarsv_contemporaneous_covariance_update
   use bvarsv_mod, only: bvarsv_prepare, bvarsv_coefficient_update
   use bvarsv_mod, only: bvarsv_ols_prior
   use bvarsv_mod, only: bvarsv_gibbs
   use bvarsv_mod, only: bvarsv_forecast
   use bvarsv_mod, only: bvarsv_irf
   use bvarsv_mod, only: bvarsv_predictive_density
   use bvarsv_mod, only: bvarsv_simulate_var1
   use bvarsv_mod, only: bvarsv_predictive_draws, bvarsv_parameter_draws
   use mixar_mod, only: mixar_bayesian_random_t, mixar_bayesian_draws_t, &
      mixar_marginal_likelihood_t, mixar_bayesian_random, &
      mixar_bayesian_sample_from_random, mixar_bayesian_sample, &
      mixar_bayesian_relabel, mixar_marginal_likelihood
   use mixar_mod, only: mixar_order_random_t, mixar_order_selection_t, &
      mixar_order_prior_flat, mixar_order_prior_ratio, &
      mixar_order_prior_poisson, mixar_order_birth_death_probability, &
      mixar_order_random, mixar_order_select_from_random, mixar_order_select
   implicit none
   private

   public :: bf_mcmc_t, bf_information_criteria_t
   public :: bf_log_prior, bf_log_posterior, bf_metropolis_sample
   public :: bf_posterior_interval, bf_information_criteria, bf_bayes_factor
   public :: bsts_mcmc_t, bsts_prediction_t
   public :: bsts_local_level, bsts_local_linear_trend, bsts_seasonal
   public :: bsts_spike_slab, bsts_dynamic_regression, bsts_predict, bsts_mbsts
   public :: bssm_mcmc_t, bssm_da_mcmc_t, bssm_approximate_mcmc
   public :: bssm_nonlinear_approximate_mcmc, bssm_nonlinear_pmmh
   public :: bssm_sde_pmmh, bssm_importance_post_correction, bssm_mcmc_diagnostics
   public :: astsa_ar_mcmc_t, astsa_sv_mcmc_t, ar_mcmc, ar_mcmc_draws
   public :: baystar_prior_t, baystar_coefficient_posterior_t
   public :: baystar_variance_posterior_t, baystar_summary_t
   public :: baystar_fit_t, baystar_simulation_t
   public :: baystar_default_prior, baystar_log_likelihood
   public :: baystar_coefficient_posterior, baystar_variance_posterior
   public :: baystar_delay_probabilities, baystar_summary
   public :: baystar_simulate_from_innovations, baystar_simulate, baystar_fit
   public :: sv_mcmc, astsa_ffbs, astsa_ffbs_draws
   public :: mts_bvar_fit_t, mts_bvar_prior_t, mts_bvar_fit, mts_minnesota_prior
   public :: bvartools_minnesota_prior_t, bvartools_inclusion_prior_t
   public :: bvartools_ssvs_result_t, bvartools_normal_posterior_t
   public :: bvartools_ssvs_prior_t
   public :: bvartools_gamma_posterior_t, bvartools_covariance_data_t
   public :: bvartools_bvs_result_t, bvartools_minnesota_prior
   public :: bvartools_bvar_draws_t, bvartools_cointegration_draw_t
   public :: bvartools_bvar_prior_bundle_t, bvartools_bvar_fit_t
   public :: bvartools_bvec_draws_t, bvartools_bvec_prior_bundle_t
   public :: bvartools_bvec_fit_t
   public :: bvartools_tvp_bvar_draws_t
   public :: bvartools_tvp_bvar_prior_bundle_t, bvartools_tvp_bvar_fit_t
   public :: bvartools_tvp_covariance_draws_t
   public :: bvartools_joint_tvp_bvar_draws_t
   public :: bvartools_structural_tvp_bvar_prior_t
   public :: bvartools_structural_tvp_bvar_fit_t
   public :: bvartools_stochastic_volatility_t
   public :: bvartools_tvp_bvec_draws_t
   public :: bvartools_tvp_bvec_prior_bundle_t, bvartools_tvp_bvec_fit_t
   public :: bvartools_dfm_draws_t
   public :: bvartools_dfm_prior_t, bvartools_dfm_grid_draws_t
   public :: bvartools_model_likelihood_data_t, bvartools_model_comparison_t
   public :: bvartools_model_comparison_set_t
   public :: bvartools_predictive_t
   public :: bvartools_irf_t
   public :: bvartools_fevd_t
   public :: bvartools_level_var_draws_t, bvartools_tvp_level_var_draws_t
   public :: bvartools_var_data_t, bvartools_vecm_data_t, bvartools_dfm_data_t
   public :: bvartools_loglik_normal_t, bvartools_kalman_dk_t
   public :: bvartools_loglik_normal, bvartools_kalman_dk
   public :: bvartools_inclusion_prior, bvartools_ssvs
   public :: bvartools_ssvs_prior
   public :: bvartools_normal_posterior, bvartools_normal_draw
   public :: bvartools_sur_normal_posterior, bvartools_sur_normal_draw
   public :: bvartools_measurement_variance_posterior
   public :: bvartools_state_variance_posterior, bvartools_gamma_precision_draw
   public :: bvartools_covar_prepare_data, bvartools_covar_vector_to_matrix
   public :: bvartools_sur_const_to_tvp, bvartools_covar_const_posterior
   public :: bvartools_covar_tvp_posterior, bvartools_bvs
   public :: bvartools_random_walk_posterior, bvartools_bvar_gibbs
   public :: bvartools_structural_bvar_gibbs
   public :: bvartools_prepare_bvar_prior, bvartools_fit_bvar
   public :: bvartools_fit_structural_bvar
   public :: bvartools_cointegration_draw
   public :: bvartools_cointegration_sur_draw
   public :: bvartools_bvec_gibbs
   public :: bvartools_prepare_bvec_prior, bvartools_fit_bvec
   public :: bvartools_fit_structural_bvec
   public :: bvartools_initial_state_posterior, bvartools_tvp_bvar_gibbs
   public :: bvartools_prepare_tvp_bvar_prior, bvartools_fit_tvp_bvar
   public :: bvartools_tvp_covariance_gibbs
   public :: bvartools_joint_tvp_bvar_gibbs
   public :: bvartools_prepare_structural_tvp_bvar_prior
   public :: bvartools_fit_structural_tvp_bvar
   public :: bvartools_fit_structural_tvp_bvec
   public :: bvartools_tvp_bvs, bvartools_tvp_covariance_bvs
   public :: bvartools_structural_impacts
   public :: bvartools_stochastic_volatility_draw
   public :: bvartools_stochastic_volatility
   public :: bvartools_stochastic_volatility_ocsn_draw
   public :: bvartools_stochastic_volatility_ocsn
   public :: bvartools_tvp_bvec_gibbs
   public :: bvartools_prepare_tvp_bvec_prior, bvartools_fit_tvp_bvec
   public :: bvartools_dfm_factor_posterior, bvartools_dfm_gibbs
   public :: bvartools_dfm_prior, bvartools_dfm_grid_gibbs
   public :: bvartools_model_comparison, bvartools_compare_models
   public :: bvartools_bvar_predictive_from_random, bvartools_bvar_predictive
   public :: bvartools_tvp_bvar_predictive_from_random
   public :: bvartools_tvp_bvar_predictive
   public :: bvartools_bvar_irf, bvartools_tvp_bvar_irf
   public :: bvartools_bvar_fevd, bvartools_tvp_bvar_fevd
   public :: bvartools_vecm_level_ar, bvartools_vecm_level_exogenous
   public :: bvartools_bvec_to_level_var, bvartools_tvp_bvec_to_level_var
   public :: bvartools_reconstruct_levels
   public :: bvartools_prepare_var, bvartools_prepare_vecm, bvartools_prepare_dfm
   public :: bvar_dummy_observations_t, bvar_evidence_t
   public :: bvar_hyper_evaluation_t, bvar_metropolis_result_t
   public :: bvar_minnesota_variance, bvar_soc_dummy, bvar_sur_dummy
   public :: bvar_conjugate_evidence
   public :: bvar_gamma_log_density, bvar_inverse_gamma_log_density
   public :: bvar_hyper_log_posterior
   public :: bvar_hierarchical_metropolis_from_random
   public :: bvar_hierarchical_metropolis
   public :: bvars_data_t, bvars_prior_t, bvars_draws_t
   public :: bvars_common_sv_t
   public :: bvars_predictive_t, bvars_fitted_t
   public :: bvars_prepare, bvars_default_prior
   public :: bvars_matrix_normal_from_standard, bvars_matrix_normal
   public :: bvars_conjugate_draws, bvars_shocks
   public :: bvars_student_log_kernel, bvars_student_t_draws
   public :: bvars_student_scale_forecast
   public :: bvars_sv_auxiliary_mixture
   public :: bvars_centered_sv_update, bvars_centered_sv_draws
   public :: bvars_noncentered_sv_update, bvars_noncentered_sv_draws
   public :: bvars_student_sv_draws
   public :: bvars_common_variance_forecast
   public :: bvars_fitted_from_random, bvars_fitted
   public :: bvars_forecast_from_random, bvars_forecast
   public :: bvars_conditional_forecast_from_random, bvars_conditional_forecast
   public :: bvars_fevd
   public :: bvarsv_data_t, bvarsv_state_draw_t, bvarsv_volatility_draw_t
   public :: bvarsv_contemporaneous_draw_t
   public :: bvarsv_covariance_draw_t
   public :: bvarsv_prior_t
   public :: bvarsv_draws_t
   public :: bvarsv_predictive_t
   public :: bvarsv_irf_t
   public :: bvarsv_simulation_t
   public :: bvarsv_predictive_draws_t, bvarsv_parameter_path_t
   public :: bvarsv_carter_kohn, bvarsv_contemporaneous_matrix
   public :: bvarsv_covariance_from_state, bvarsv_log_volatility_update
   public :: bvarsv_contemporaneous_update
   public :: bvarsv_random_walk_covariance
   public :: bvarsv_contemporaneous_covariance_update
   public :: bvarsv_prepare, bvarsv_coefficient_update
   public :: bvarsv_ols_prior
   public :: bvarsv_gibbs
   public :: bvarsv_forecast
   public :: bvarsv_irf
   public :: bvarsv_predictive_density
   public :: bvarsv_simulate_var1
   public :: bvarsv_predictive_draws, bvarsv_parameter_draws
   public :: mixar_bayesian_random_t, mixar_bayesian_draws_t
   public :: mixar_marginal_likelihood_t, mixar_bayesian_random
   public :: mixar_bayesian_sample_from_random, mixar_bayesian_sample
   public :: mixar_bayesian_relabel, mixar_marginal_likelihood
   public :: mixar_order_random_t, mixar_order_selection_t
   public :: mixar_order_prior_flat, mixar_order_prior_ratio
   public :: mixar_order_prior_poisson, mixar_order_birth_death_probability
   public :: mixar_order_random, mixar_order_select_from_random
   public :: mixar_order_select
end module bayesian_time_series_mod
