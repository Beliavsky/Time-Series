! SPDX-License-Identifier: MIT
! SPDX-FileComment: Topic facade for multivariate time-series models.
module multivariate_mod
   !! Re-export VAR, VARMA, VECM, cointegration, and factor-model algorithms.
   use mts_mod, only: mts_var_fit_t, mts_var_forecast_t, mts_vecm_fit_t, &
      mts_var, mts_vars, mts_var_forecast, mts_var_order, mts_var_irf, mts_fevd, &
      mts_varx, mts_varx_forecast, mts_varma_fit, mts_varma_forecast, &
      mts_vecm_fit, mts_vecm_forecast, mts_factor_fit, mts_factor_forecast
   use urca_mod, only: johansen_result_t, johansen_test
   use fcvar_mod, only: fcvar_fit_t, fcvar_estimation_t, fcvar_rank_tests_t, &
      fcvar_estimate, fcvar_forecast, fcvar_rank_tests, fcvar_lag_select, &
      fcvar_likelihood_ratio, fcvar_simulate
   use mar_mod, only: mar_fit_t, mar_modes_t, mar_simulation_t, mar_estimate, &
      mar_eigenmodes, mar_pca, mar_simulate, display_mar_simulation
   use nts_mod, only: nts_mtar_simulation_t, nts_mtar_model_t, &
      nts_mtar_search_t, nts_mtar_forecast_t, nts_mtar_refinement_t, &
      nts_mtar_simulate_from_standard, nts_mtar_simulate, &
      nts_mtar_estimate, nts_mtar_threshold_search, &
      nts_mtar_forecast_draws, nts_mtar_forecast, nts_mtar_refine
   use tsdyn_mod, only: tsdyn_tvecm_model_t, tsdyn_tvecm_forecast_t, &
      tsdyn_tvecm_simulation_t, tsdyn_threshold_test_t, tsdyn_girf_t, &
      tsdyn_irf_t, &
      tsdyn_fevd_t, tsdyn_tvar_selection_t, tsdyn_regime_count_test_t, &
      tsdyn_rank_test_t, tsdyn_rank_selection_t, &
      tsdyn_tvar_inference_t, tsdyn_tvecm_inference_t, &
      tsdyn_regime_path_t, &
      tsdyn_tvar_bootstrap_t, tsdyn_tvecm_bootstrap_t, &
      tsdyn_tvecm_fit, &
      tsdyn_tvecm_forecast, &
      tsdyn_tvecm_simulate_from_innovations, tsdyn_tvar_fit, &
      tsdyn_tvar_simulate_from_standard, tsdyn_tvar_forecast_draws, &
      tsdyn_tvar_select, tsdyn_tvar_restricted_fit, tsdyn_tvar_lr_test, &
      tsdyn_tvar_bootstrap, tsdyn_tvecm_bootstrap, &
      tsdyn_hansen_seo_test, tsdyn_seo_test, &
      tsdyn_rank_test, tsdyn_rank_select, &
      tsdyn_tvar_inference, tsdyn_tvecm_inference, &
      tsdyn_tvar_regimes, tsdyn_tvecm_regimes, &
      tsdyn_tvar_girf_from_innovations, &
      tsdyn_tvecm_girf_from_innovations, tsdyn_nonlinear_fevd, &
      tsdyn_tvar_regime_irf, tsdyn_tvecm_regime_irf, &
      tsdyn_tvar_irf_bootstrap, tsdyn_tvecm_irf_bootstrap
   use bigtime_mod, only: bigtime_var_fit_t, bigtime_varx_fit_t, &
      bigtime_sparse_var, bigtime_sparse_varx, bigtime_var_cv, &
      bigtime_var_forecast, bigtime_varx_forecast, bigtime_var_simulate
   use bigvar_mod, only: bigvar_fit_t, bigvar_varx_fit_t, &
      bigvar_structured_var, bigvar_structured_varx, bigvar_forecast, &
      bigvar_varx_forecast, bigvar_var_validate, bigvar_var_simulate
   use varshrink_mod, only: varshrink_var_fit_t, varshrink_var_ridge, &
      varshrink_semibayes, varshrink_covariance_shrink
   use vars_mod, only: vars_fit_t, vars_selection_t, vars_bq_t, vars_svar_t, &
      vars_svec_t, vars_vec2var_t, vars_irf_bootstrap_t, &
      vars_structural_bootstrap_t, vars_svec_bootstrap_t, &
      vars_structural_irf_t, vars_fevd_t, vars_fit, &
      vars_select, vars_restrict, vars_restrict_ser, vars_phi, vars_psi, &
      vars_roots, vars_bq, vars_svar, vars_svar_scoring, vars_svec, &
      vars_vec2var, vars_irf_bootstrap, vars_svar_bootstrap, vars_svec_bootstrap, &
      vars_svar_irf, vars_svec_irf, vars_svar_fevd, vars_svec_fevd, &
      display_vars_fit
   use var_etp_mod, only: var_etp_bias_result_t, var_etp_test_result_t, &
      var_etp_pope, var_etp_bootstrap_bias, var_etp_restrict, &
      var_etp_wald_test, var_etp_lr_test
   use bvartools_mod, only: bvartools_irf_t, bvartools_fevd_t, &
      bvartools_loglik_normal_t, bvartools_loglik_normal, &
      bvartools_level_var_draws_t, bvartools_tvp_level_var_draws_t, &
      bvartools_var_data_t, bvartools_vecm_data_t, bvartools_dfm_data_t, &
      bvartools_dfm_prior_t, bvartools_dfm_grid_draws_t, &
      bvartools_model_likelihood_data_t, bvartools_model_comparison_t, &
      bvartools_model_comparison_set_t, &
      bvartools_ssvs_prior_t, bvartools_ssvs_prior, &
      bvartools_bvar_prior_bundle_t, bvartools_bvar_fit_t, &
      bvartools_bvec_prior_bundle_t, bvartools_bvec_fit_t, &
      bvartools_tvp_bvar_prior_bundle_t, bvartools_tvp_bvar_fit_t, &
      bvartools_joint_tvp_bvar_draws_t, bvartools_joint_tvp_bvar_gibbs, &
      bvartools_structural_tvp_bvar_prior_t, &
      bvartools_structural_tvp_bvar_fit_t, &
      bvartools_prepare_structural_tvp_bvar_prior, &
      bvartools_fit_structural_tvp_bvar, &
      bvartools_fit_structural_tvp_bvec, &
      bvartools_structural_bvar_gibbs, &
      bvartools_fit_structural_bvar, bvartools_fit_structural_bvec, &
      bvartools_structural_impacts, &
      bvartools_tvp_bvec_prior_bundle_t, bvartools_tvp_bvec_fit_t, &
      bvartools_bvar_irf, bvartools_tvp_bvar_irf, &
      bvartools_bvar_fevd, bvartools_tvp_bvar_fevd, &
      bvartools_vecm_level_ar, bvartools_vecm_level_exogenous, &
      bvartools_bvec_to_level_var, bvartools_tvp_bvec_to_level_var, &
      bvartools_reconstruct_levels, bvartools_prepare_var, &
      bvartools_prepare_vecm, bvartools_prepare_dfm, bvartools_dfm_prior, &
      bvartools_dfm_grid_gibbs
   use bvartools_mod, only: bvartools_prepare_bvar_prior, bvartools_fit_bvar, &
      bvartools_prepare_bvec_prior, bvartools_fit_bvec
   use bvartools_mod, only: bvartools_prepare_tvp_bvar_prior, &
      bvartools_fit_tvp_bvar, bvartools_prepare_tvp_bvec_prior, &
      bvartools_fit_tvp_bvec
   use bvartools_mod, only: bvartools_model_comparison, bvartools_compare_models
   use bvars_mod, only: bvars_data_t, bvars_prior_t, bvars_draws_t, &
      bvars_common_sv_t, &
      bvars_predictive_t, bvars_fitted_t, bvars_prepare, bvars_default_prior, &
      bvars_conjugate_draws, bvars_student_t_draws, &
      bvars_student_scale_forecast, bvars_sv_auxiliary_mixture, &
      bvars_shocks, bvars_fitted, bvars_forecast, &
      bvars_conditional_forecast, &
      bvars_fevd, bvars_centered_sv_draws, bvars_noncentered_sv_draws, &
      bvars_student_sv_draws, &
      bvars_common_variance_forecast
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
   use gmvarkit_mod, only: gmvarkit_model_t, gmvarkit_regime_moments_t, &
      gmvarkit_location_parameters_t, &
      gmvarkit_evaluation_t, gmvarkit_simulation_t, gmvarkit_forecast_t, &
      gmvarkit_girf_t, gmvarkit_gfevd_t, &
      gmvarkit_girf_inference_t, gmvarkit_gfevd_inference_t, &
      gmvarkit_linear_irf_t, &
      gmvarkit_unconditional_moments_t, gmvarkit_pearson_residuals_t, &
      gmvarkit_information_criteria_t, &
      gmvarkit_structural_t, &
      gmvarkit_linear_irf_bootstrap_t, &
      gmvarkit_structural_constraints_t, gmvarkit_structural_fit_t, &
      gmvarkit_structural_multistart_fit_t, &
      gmvarkit_likelihood_profile_t, &
      gmvarkit_regime_conversion_t, &
      gmvarkit_companion_eigen_t, gmvarkit_covariance_eigen_t, &
      gmvarkit_fit_t, &
      gmvarkit_genetic_fit_t, &
      gmvarkit_multistart_fit_t, &
      gmvarkit_constraints_t, &
      gmvarkit_inference_t, gmvarkit_hypothesis_test_t, &
      gmvarkit_score_t, gmvarkit_quantile_residuals_t, &
      gmvarkit_residual_tests_t, &
      gmvarkit_regime_moments, gmvarkit_evaluate, gmvarkit_simulate, &
      gmvarkit_location_parameters, gmvarkit_model_from_regime_means, &
      gmvarkit_forecast, gmvarkit_estimate, gmvarkit_genetic_estimate, &
      gmvarkit_multistart_estimate, &
      gmvarkit_estimate_constrained, gmvarkit_inference, &
      gmvarkit_wald_test, gmvarkit_likelihood_ratio, &
      gmvarkit_score_matrix, gmvarkit_rao_test, &
      gmvarkit_quantile_residuals, gmvarkit_quantile_residual_tests, &
      gmvarkit_girf, gmvarkit_gfevd, gmvarkit_linear_irf, &
      gmvarkit_girf_inference, gmvarkit_gfevd_inference, &
      gmvarkit_unconditional_moments, gmvarkit_pearson_residuals, &
      gmvarkit_information_criteria
   use gmvarkit_mod, only: gmvarkit_structural_covariances, &
      gmvarkit_identify_structural, gmvarkit_reorder_structural, &
      gmvarkit_swap_structural_signs, gmvarkit_rebase_structural, &
      gmvarkit_linear_irf_bootstrap
   use gmvarkit_mod, only: gmvarkit_estimate_structural, &
      gmvarkit_profile_likelihood
   use gmvarkit_mod, only: gmvarkit_structural_multistart_estimate
   use gmvarkit_mod, only: gmvarkit_convert_student_regimes
   use gmvarkit_mod, only: gmvarkit_companion_eigenvalues, &
      gmvarkit_covariance_eigenvalues
   use mixar_mod, only: mixar_var_model_t, mixar_var_filter_t, &
      mixar_var_fit_t, mixar_var_simulation_t, mixar_var_forecast_t, &
      mixar_var_diagnostics_t, mixar_var_model, &
      mixar_var_component_locations, mixar_var_filter, &
      mixar_var_log_likelihood, mixar_var_fit, mixar_var_is_stable, &
      mixar_var_simulate_from_draws, mixar_var_simulate, &
      mixar_var_forecast_from_draws, mixar_var_forecast, mixar_var_diagnose
   use nlints_mod, only: nlints_varmlp_t, nlints_varmlp_fit, &
      nlints_varmlp_update, nlints_varmlp_predict, nlints_varmlp_forecast, &
      nlints_varmlp_forecast_table, nlints_varmlp_save, nlints_varmlp_load
   use starvars_mod, only: starvars_start_t, starvars_fit_t, &
      starvars_method_nls, starvars_method_ml, starvars_logistic, &
      starvars_starting, starvars_fit
   implicit none
   private

   interface display
      procedure :: display_vars_fit
      procedure :: display_mar_simulation
   end interface display

   public :: mts_var_fit_t, mts_var_forecast_t, mts_vecm_fit_t
   public :: mts_var, mts_vars, mts_var_forecast, mts_var_order, mts_var_irf, mts_fevd
   public :: mts_varx, mts_varx_forecast, mts_varma_fit, mts_varma_forecast
   public :: mts_vecm_fit, mts_vecm_forecast, mts_factor_fit, mts_factor_forecast
   public :: johansen_result_t, johansen_test
   public :: fcvar_fit_t, fcvar_estimation_t, fcvar_rank_tests_t
   public :: fcvar_estimate, fcvar_forecast, fcvar_rank_tests, fcvar_lag_select
   public :: fcvar_likelihood_ratio, fcvar_simulate
   public :: mar_fit_t, mar_modes_t, mar_simulation_t
   public :: mar_estimate, mar_eigenmodes, mar_pca
   public :: mar_simulate
   public :: nts_mtar_simulation_t, nts_mtar_model_t
   public :: nts_mtar_search_t, nts_mtar_forecast_t
   public :: nts_mtar_refinement_t
   public :: nts_mtar_simulate_from_standard, nts_mtar_simulate
   public :: nts_mtar_estimate, nts_mtar_threshold_search
   public :: nts_mtar_forecast_draws, nts_mtar_forecast
   public :: nts_mtar_refine
   public :: tsdyn_tvecm_model_t, tsdyn_tvecm_forecast_t
   public :: tsdyn_tvecm_simulation_t
   public :: tsdyn_tvar_bootstrap_t, tsdyn_tvecm_bootstrap_t
   public :: tsdyn_tvecm_fit, tsdyn_tvecm_forecast
   public :: tsdyn_tvecm_simulate_from_innovations
   public :: tsdyn_threshold_test_t, tsdyn_girf_t, tsdyn_irf_t, tsdyn_fevd_t
   public :: tsdyn_rank_test_t, tsdyn_rank_selection_t
   public :: tsdyn_rank_test, tsdyn_rank_select
   public :: tsdyn_tvar_inference_t, tsdyn_tvar_inference
   public :: tsdyn_tvecm_inference_t, tsdyn_tvecm_inference
   public :: tsdyn_regime_path_t, tsdyn_tvar_regimes, tsdyn_tvecm_regimes
   public :: tsdyn_tvar_fit, tsdyn_tvar_simulate_from_standard
   public :: tsdyn_tvar_forecast_draws
   public :: tsdyn_tvar_bootstrap, tsdyn_tvecm_bootstrap
   public :: tsdyn_tvar_selection_t, tsdyn_tvar_select
   public :: tsdyn_regime_count_test_t
   public :: tsdyn_tvar_restricted_fit, tsdyn_tvar_lr_test
   public :: tsdyn_hansen_seo_test, tsdyn_seo_test
   public :: tsdyn_tvar_girf_from_innovations
   public :: tsdyn_tvecm_girf_from_innovations, tsdyn_nonlinear_fevd
   public :: tsdyn_tvar_regime_irf, tsdyn_tvecm_regime_irf
   public :: tsdyn_tvar_irf_bootstrap, tsdyn_tvecm_irf_bootstrap
   public :: bigtime_var_fit_t, bigtime_varx_fit_t
   public :: bigtime_sparse_var, bigtime_sparse_varx, bigtime_var_cv
   public :: bigtime_var_forecast, bigtime_varx_forecast, bigtime_var_simulate
   public :: bigvar_fit_t, bigvar_varx_fit_t
   public :: bigvar_structured_var, bigvar_structured_varx, bigvar_forecast
   public :: bigvar_varx_forecast, bigvar_var_validate, bigvar_var_simulate
   public :: varshrink_var_fit_t, varshrink_var_ridge
   public :: varshrink_semibayes, varshrink_covariance_shrink
   public :: vars_fit_t, vars_selection_t, vars_bq_t, vars_svar_t, vars_svec_t
   public :: vars_vec2var_t, vars_irf_bootstrap_t, vars_structural_bootstrap_t
   public :: vars_svec_bootstrap_t, vars_structural_irf_t, vars_fevd_t
   public :: vars_fit, vars_select, vars_restrict, vars_restrict_ser
   public :: vars_phi, vars_psi, vars_roots, vars_bq, vars_svar, vars_svar_scoring
   public :: vars_svec, vars_vec2var, vars_irf_bootstrap, vars_svar_bootstrap
   public :: vars_svec_bootstrap, vars_svar_irf, vars_svec_irf
   public :: vars_svar_fevd, vars_svec_fevd
   public :: display
   public :: var_etp_bias_result_t, var_etp_test_result_t
   public :: var_etp_pope, var_etp_bootstrap_bias, var_etp_restrict
   public :: var_etp_wald_test, var_etp_lr_test
   public :: bvartools_irf_t, bvartools_bvar_irf, bvartools_tvp_bvar_irf
   public :: bvartools_loglik_normal_t, bvartools_loglik_normal
   public :: bvartools_fevd_t, bvartools_bvar_fevd, bvartools_tvp_bvar_fevd
   public :: bvartools_level_var_draws_t, bvartools_tvp_level_var_draws_t
   public :: bvartools_vecm_level_ar, bvartools_vecm_level_exogenous
   public :: bvartools_bvec_to_level_var, bvartools_tvp_bvec_to_level_var
   public :: bvartools_reconstruct_levels
   public :: bvartools_var_data_t, bvartools_vecm_data_t, bvartools_dfm_data_t
   public :: bvartools_prepare_var, bvartools_prepare_vecm, bvartools_prepare_dfm
   public :: bvartools_dfm_prior_t, bvartools_dfm_grid_draws_t
   public :: bvartools_dfm_prior, bvartools_dfm_grid_gibbs
   public :: bvartools_model_likelihood_data_t, bvartools_model_comparison_t
   public :: bvartools_model_comparison_set_t
   public :: bvartools_model_comparison, bvartools_compare_models
   public :: bvartools_ssvs_prior_t, bvartools_ssvs_prior
   public :: bvartools_bvar_prior_bundle_t, bvartools_bvar_fit_t
   public :: bvartools_bvec_prior_bundle_t, bvartools_bvec_fit_t
   public :: bvartools_prepare_bvar_prior, bvartools_fit_bvar
   public :: bvartools_prepare_bvec_prior, bvartools_fit_bvec
   public :: bvartools_tvp_bvar_prior_bundle_t, bvartools_tvp_bvar_fit_t
   public :: bvartools_joint_tvp_bvar_draws_t, bvartools_joint_tvp_bvar_gibbs
   public :: bvartools_structural_tvp_bvar_prior_t
   public :: bvartools_structural_tvp_bvar_fit_t
   public :: bvartools_prepare_structural_tvp_bvar_prior
   public :: bvartools_fit_structural_tvp_bvar
   public :: bvartools_fit_structural_tvp_bvec
   public :: bvartools_structural_bvar_gibbs
   public :: bvartools_fit_structural_bvar, bvartools_fit_structural_bvec
   public :: bvartools_structural_impacts
   public :: bvartools_tvp_bvec_prior_bundle_t, bvartools_tvp_bvec_fit_t
   public :: bvartools_prepare_tvp_bvar_prior, bvartools_fit_tvp_bvar
   public :: bvartools_prepare_tvp_bvec_prior, bvartools_fit_tvp_bvec
   public :: bvars_data_t, bvars_prior_t, bvars_draws_t
   public :: bvars_common_sv_t
   public :: bvars_predictive_t, bvars_fitted_t
   public :: bvars_prepare, bvars_default_prior, bvars_conjugate_draws
   public :: bvars_student_t_draws, bvars_student_scale_forecast
   public :: bvars_sv_auxiliary_mixture
   public :: bvars_centered_sv_draws, bvars_common_variance_forecast
   public :: bvars_noncentered_sv_draws
   public :: bvars_student_sv_draws
   public :: bvars_shocks, bvars_fitted, bvars_forecast
   public :: bvars_conditional_forecast, bvars_fevd
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
   public :: gmvarkit_model_t, gmvarkit_regime_moments_t
   public :: gmvarkit_location_parameters_t
   public :: gmvarkit_evaluation_t
   public :: gmvarkit_simulation_t, gmvarkit_forecast_t
   public :: gmvarkit_girf_t, gmvarkit_gfevd_t
   public :: gmvarkit_girf_inference_t, gmvarkit_gfevd_inference_t
   public :: gmvarkit_linear_irf_t
   public :: gmvarkit_unconditional_moments_t, gmvarkit_pearson_residuals_t
   public :: gmvarkit_information_criteria_t
   public :: gmvarkit_structural_t
   public :: gmvarkit_linear_irf_bootstrap_t
   public :: gmvarkit_structural_constraints_t, gmvarkit_structural_fit_t
   public :: gmvarkit_structural_multistart_fit_t
   public :: gmvarkit_likelihood_profile_t
   public :: gmvarkit_regime_conversion_t
   public :: gmvarkit_companion_eigen_t, gmvarkit_covariance_eigen_t
   public :: gmvarkit_fit_t
   public :: gmvarkit_genetic_fit_t
   public :: gmvarkit_multistart_fit_t
   public :: gmvarkit_constraints_t
   public :: gmvarkit_inference_t, gmvarkit_hypothesis_test_t
   public :: gmvarkit_score_t, gmvarkit_quantile_residuals_t
   public :: gmvarkit_residual_tests_t
   public :: gmvarkit_regime_moments, gmvarkit_evaluate
   public :: gmvarkit_location_parameters, gmvarkit_model_from_regime_means
   public :: gmvarkit_simulate, gmvarkit_forecast
   public :: gmvarkit_girf, gmvarkit_gfevd
   public :: gmvarkit_girf_inference, gmvarkit_gfevd_inference
   public :: gmvarkit_linear_irf
   public :: gmvarkit_unconditional_moments, gmvarkit_pearson_residuals
   public :: gmvarkit_information_criteria
   public :: gmvarkit_structural_covariances, gmvarkit_identify_structural
   public :: gmvarkit_reorder_structural, gmvarkit_swap_structural_signs
   public :: gmvarkit_rebase_structural
   public :: gmvarkit_linear_irf_bootstrap
   public :: gmvarkit_estimate_structural
   public :: gmvarkit_structural_multistart_estimate
   public :: gmvarkit_profile_likelihood
   public :: gmvarkit_convert_student_regimes
   public :: gmvarkit_companion_eigenvalues, gmvarkit_covariance_eigenvalues
   public :: gmvarkit_estimate
   public :: gmvarkit_genetic_estimate
   public :: gmvarkit_multistart_estimate
   public :: gmvarkit_estimate_constrained
   public :: gmvarkit_inference, gmvarkit_wald_test
   public :: gmvarkit_likelihood_ratio
   public :: gmvarkit_score_matrix, gmvarkit_rao_test
   public :: gmvarkit_quantile_residuals, gmvarkit_quantile_residual_tests
   public :: mixar_var_model_t, mixar_var_filter_t, mixar_var_fit_t
   public :: mixar_var_simulation_t, mixar_var_forecast_t
   public :: mixar_var_diagnostics_t, mixar_var_model
   public :: mixar_var_component_locations, mixar_var_filter
   public :: mixar_var_log_likelihood, mixar_var_fit, mixar_var_is_stable
   public :: mixar_var_simulate_from_draws, mixar_var_simulate
   public :: mixar_var_forecast_from_draws, mixar_var_forecast
   public :: mixar_var_diagnose
   public :: nlints_varmlp_t, nlints_varmlp_fit, nlints_varmlp_update
   public :: nlints_varmlp_predict, nlints_varmlp_forecast
   public :: nlints_varmlp_forecast_table
   public :: nlints_varmlp_save, nlints_varmlp_load
   public :: starvars_start_t, starvars_fit_t
   public :: starvars_method_nls, starvars_method_ml
   public :: starvars_logistic, starvars_starting, starvars_fit
end module multivariate_mod
