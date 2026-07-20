! SPDX-License-Identifier: MIT
! SPDX-FileComment: Topic facade for time-series diagnostics.
module diagnostics_mod
   !! Re-export residual, portmanteau, accuracy, and model diagnostic algorithms.
   use time_series_diagnostics_mod, only: weighted_box_test_t, &
      multivariate_white_noise_test_t, weighted_box_test, &
      multivariate_white_noise_test, box_test_box_pierce, box_test_ljung_box, &
      box_test_monti
   use forecast_mod, only: accuracy_result_t, dm_result_t, forecast_accuracy, &
      rolling_forecast_accuracy, dm_test
   use itsmr_mod, only: itsmr_randomness_tests_t, residual_randomness_tests
   use astsa_mod, only: astsa_sarima_diagnostics_t, sarima_diagnostics, &
      astsa_linearity_test_t, test_linearity
   use arfima_mod, only: arfima_diagnostics_t, arfima_diagnostics, &
      arfima_mode_diagnostics
   use garma_mod, only: garma_gof_t, garma_goodness_of_fit, garma_accuracy
   use mts_mod, only: mts_diagnostic_t, mts_arch_test, mts_march_test, &
      mts_granger_test, mts_diagnostic, mts_var_backtest
   use tsissm_mod, only: tsissm_diagnostics_t, tsissm_diagnose, &
      tsissm_diagnose_structural, tsissm_weighted_box_test
   use vars_mod, only: vars_serial_tests_t, vars_normality_tests_t, &
      vars_arch_tests_t, vars_stability_t, vars_serial_test, &
      vars_normality_test, vars_arch_test, vars_instantaneous_causality, &
      vars_ols_cusum, vars_granger_test, vars_granger_bootstrap
   use nts_mod, only: nts_nonlinearity_test_t, nts_rank_portmanteau_t, &
      nts_prnd_test_t, nts_threshold_test, nts_tsay_test, &
      nts_rank_portmanteau, nts_quadratic_f_test, nts_prnd_test
   use tsdyn_mod, only: tsdyn_bbc_test_t, tsdyn_kapshin_test_t, &
      tsdyn_delta_test_t, tsdyn_delta_linear_test_t, tsdyn_bbc_test, &
      tsdyn_kapshin_test, tsdyn_delta, tsdyn_delta_test, &
      tsdyn_delta_linear, tsdyn_delta_linear_test, &
      tsdyn_regime_diagnostics_t, tsdyn_ar_regime_diagnostics, &
      tsdyn_setar_regime_diagnostics, tsdyn_lstar_regime_diagnostics, &
      tsdyn_star_regime_diagnostics
   use tsdyn_mod, only: tsdyn_rank_test_t, tsdyn_rank_test
   use tsdyn_mod, only: tsdyn_setar_variance_t, &
      tsdyn_setar_residual_variance
   use mixar_mod, only: mixar_inference_t, mixar_diagnostics_t, &
      mixar_bic_selection_t, mixar_observed_inference, &
      mixar_seasonal_observed_inference, mixar_regression_observed_inference, &
      mixar_diagnose, mixar_bic, mixar_select_bic
   use mixar_mod, only: mixar_var_diagnostics_t, mixar_var_diagnose
   use ugmar_mod, only: ugmar_quantile_residuals_t, &
      ugmar_hypothesis_test_t, ugmar_residual_tests_t, &
      ugmar_inference_t, ugmar_quantile_residuals, &
      ugmar_quantile_residual_tests, ugmar_inference, ugmar_wald_test, &
      ugmar_likelihood_ratio
   use tseriestarma_mod, only: tseriestarma_tar_test_t, &
      tseriestarma_tarma_test_t, tseriestarma_tar_test, &
      tseriestarma_tarma_test, tseriestarma_tar_bootstrap_t, &
      tseriestarma_tar_bootstrap_from_random, tseriestarma_tar_bootstrap, &
      tseriestarma_bootstrap_iid, tseriestarma_bootstrap_hansen, &
      tseriestarma_bootstrap_rademacher, tseriestarma_bootstrap_normal, &
      tseriestarma_unit_root_test_t, tseriestarma_unit_root_bootstrap_t, &
      tseriestarma_unit_root_test, &
      tseriestarma_unit_root_bootstrap_from_random, &
      tseriestarma_unit_root_bootstrap, &
      tseriestarma_unit_root_critical_values, &
      tseriestarma_unit_root_critical_probabilities, &
      tseriestarma_arma_garch_fit_t, tseriestarma_garch_test_t, &
      tseriestarma_arma_garch_fit, tseriestarma_garch_test, &
      tseriestarma_garch_critical_values, &
      tseriestarma_garch_critical_probabilities
   use nlints_mod, only: nlints_causality_test_t, &
      nlints_neural_causality_t, nlints_entropy_discrete, &
      nlints_joint_entropy_discrete, nlints_mutual_information_discrete, &
      nlints_multivariate_information_discrete, &
      nlints_transfer_entropy_discrete, nlints_entropy_continuous, &
      nlints_joint_entropy_continuous, &
      nlints_mutual_information_continuous, &
      nlints_transfer_entropy_continuous, nlints_granger_test, &
      nlints_neural_granger_test, nlints_df_test
   use echos_mod, only: echos_kpss_t, echos_kpss, &
      echos_estimate_differences
   use starvars_mod, only: starvars_joint_test_t, &
      starvars_long_run_variance_t, starvars_cumsum_t, &
      starvars_joint_linearity_test, starvars_long_run_variance, &
      starvars_multiple_cumsum
   use rugarch_mod, only: rugarch_berkowitz_test_t, &
      rugarch_direction_test_t, rugarch_var_test_t, rugarch_es_test_t, &
      rugarch_direction_pt, rugarch_direction_ag, rugarch_berkowitz_test, &
      rugarch_direction_test, rugarch_var_test, rugarch_es_test
   use rugarch_diagnostics_mod, only: rugarch_nyblom_test_t, &
      rugarch_gof_test_t, rugarch_var_duration_test_t, rugarch_gmm_test_t, &
      rugarch_hong_li_test_t, rugarch_nyblom_test, rugarch_gof_test, &
      rugarch_var_duration_test, rugarch_gmm_test, rugarch_hong_li_test
   implicit none
   private

   public :: weighted_box_test_t, multivariate_white_noise_test_t
   public :: weighted_box_test, multivariate_white_noise_test
   public :: box_test_box_pierce, box_test_ljung_box, box_test_monti
   public :: accuracy_result_t, dm_result_t, forecast_accuracy, dm_test
   public :: rolling_forecast_accuracy
   public :: itsmr_randomness_tests_t, residual_randomness_tests
   public :: astsa_sarima_diagnostics_t, sarima_diagnostics
   public :: astsa_linearity_test_t, test_linearity
   public :: arfima_diagnostics_t, arfima_diagnostics, arfima_mode_diagnostics
   public :: garma_gof_t, garma_goodness_of_fit, garma_accuracy
   public :: mts_diagnostic_t, mts_arch_test, mts_march_test
   public :: mts_granger_test, mts_diagnostic, mts_var_backtest
   public :: tsissm_diagnostics_t, tsissm_diagnose
   public :: tsissm_diagnose_structural, tsissm_weighted_box_test
   public :: vars_serial_tests_t, vars_normality_tests_t, vars_arch_tests_t
   public :: vars_stability_t, vars_serial_test, vars_normality_test
   public :: vars_arch_test, vars_instantaneous_causality, vars_ols_cusum
   public :: vars_granger_test, vars_granger_bootstrap
   public :: nts_nonlinearity_test_t, nts_threshold_test, nts_tsay_test
   public :: nts_rank_portmanteau_t, nts_prnd_test_t
   public :: nts_rank_portmanteau, nts_quadratic_f_test, nts_prnd_test
   public :: tsdyn_bbc_test_t, tsdyn_kapshin_test_t
   public :: tsdyn_bbc_test, tsdyn_kapshin_test
   public :: tsdyn_delta_test_t, tsdyn_delta_linear_test_t
   public :: tsdyn_delta, tsdyn_delta_test
   public :: tsdyn_delta_linear, tsdyn_delta_linear_test
   public :: tsdyn_regime_diagnostics_t, tsdyn_ar_regime_diagnostics
   public :: tsdyn_setar_regime_diagnostics, tsdyn_lstar_regime_diagnostics
   public :: tsdyn_star_regime_diagnostics
   public :: tsdyn_rank_test_t, tsdyn_rank_test
   public :: tsdyn_setar_variance_t, tsdyn_setar_residual_variance
   public :: mixar_inference_t, mixar_diagnostics_t, mixar_bic_selection_t
   public :: mixar_observed_inference, mixar_seasonal_observed_inference
   public :: mixar_regression_observed_inference, mixar_diagnose
   public :: mixar_bic, mixar_select_bic
   public :: mixar_var_diagnostics_t, mixar_var_diagnose
   public :: ugmar_quantile_residuals_t, ugmar_hypothesis_test_t
   public :: ugmar_residual_tests_t, ugmar_inference_t
   public :: ugmar_quantile_residuals
   public :: ugmar_quantile_residual_tests
   public :: ugmar_inference, ugmar_wald_test, ugmar_likelihood_ratio
   public :: tseriestarma_tar_test_t, tseriestarma_tar_test
   public :: tseriestarma_tarma_test_t, tseriestarma_tarma_test
   public :: tseriestarma_tar_bootstrap_t
   public :: tseriestarma_tar_bootstrap_from_random
   public :: tseriestarma_tar_bootstrap
   public :: tseriestarma_bootstrap_iid, tseriestarma_bootstrap_hansen
   public :: tseriestarma_bootstrap_rademacher, tseriestarma_bootstrap_normal
   public :: tseriestarma_unit_root_test_t, tseriestarma_unit_root_bootstrap_t
   public :: tseriestarma_unit_root_test
   public :: tseriestarma_unit_root_bootstrap_from_random
   public :: tseriestarma_unit_root_bootstrap
   public :: tseriestarma_unit_root_critical_values
   public :: tseriestarma_unit_root_critical_probabilities
   public :: tseriestarma_arma_garch_fit_t, tseriestarma_garch_test_t
   public :: tseriestarma_arma_garch_fit, tseriestarma_garch_test
   public :: tseriestarma_garch_critical_values
   public :: tseriestarma_garch_critical_probabilities
   public :: nlints_causality_test_t, nlints_neural_causality_t
   public :: nlints_entropy_discrete, nlints_joint_entropy_discrete
   public :: nlints_mutual_information_discrete
   public :: nlints_multivariate_information_discrete
   public :: nlints_transfer_entropy_discrete
   public :: nlints_entropy_continuous, nlints_joint_entropy_continuous
   public :: nlints_mutual_information_continuous
   public :: nlints_transfer_entropy_continuous
   public :: nlints_granger_test, nlints_neural_granger_test, nlints_df_test
   public :: echos_kpss_t, echos_kpss, echos_estimate_differences
   public :: starvars_joint_test_t, starvars_long_run_variance_t
   public :: starvars_cumsum_t, starvars_joint_linearity_test
   public :: starvars_long_run_variance, starvars_multiple_cumsum
   public :: rugarch_berkowitz_test_t, rugarch_direction_test_t
   public :: rugarch_var_test_t, rugarch_es_test_t
   public :: rugarch_direction_pt, rugarch_direction_ag
   public :: rugarch_berkowitz_test, rugarch_direction_test
   public :: rugarch_var_test, rugarch_es_test
   public :: rugarch_nyblom_test_t, rugarch_gof_test_t
   public :: rugarch_var_duration_test_t, rugarch_gmm_test_t
   public :: rugarch_hong_li_test_t, rugarch_nyblom_test
   public :: rugarch_gof_test, rugarch_var_duration_test
   public :: rugarch_gmm_test, rugarch_hong_li_test
end module diagnostics_mod
