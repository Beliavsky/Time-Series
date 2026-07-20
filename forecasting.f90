! SPDX-License-Identifier: MIT
! SPDX-FileComment: Topic facade for forecasting algorithms.
module forecasting_mod
   !! Re-export principal univariate, multivariate, and structural forecasts.
   use utils_mod, only: real_vector_t
   use rolling_forecast_mod, only: rolling_forecast_result_t
   use forecast_mod, only: forecast_result_t, accuracy_result_t, meanf, naive, &
      rwf, snaive, ses, holt, holt_winters, croston, thetaf, &
      rolling_forecast_accuracy
   use itsmr_mod, only: itsmr_forecast_t, arma_forecast, arar_forecast, &
      transformed_arma_forecast
   use astsa_mod, only: astsa_sarima_forecast_t, sarima_forecast
   use arfima_mod, only: arfima_forecast_t, arfima_forecast, &
      arfima_average_forecast, arfima_forecast_paths
   use garma_mod, only: garma_forecast_t, garma_forecast, garma_regression_forecast
   use tfarima_mod, only: tfarima_forecast_t, tfarima_arima_forecast, &
      tfarima_exact_transfer_forecast, tfarima_structural_forecast
   use robustarima_mod, only: robustarima_forecast_t, robustarima_forecast
   use expar_mod, only: expar_forecast_t, expar_forecast
   use mts_mod, only: mts_var_forecast_t, mts_var_forecast, mts_varx_forecast, &
      mts_varma_forecast, mts_vecm_forecast, mts_factor_forecast
   use bigtime_mod, only: bigtime_forecast_t, bigtime_var_forecast, &
      bigtime_varx_forecast, bigtime_varma_forecast
   use bigvar_mod, only: bigvar_forecast, bigvar_direct_forecast, &
      bigvar_varx_forecast, bigvar_var_forecast_interval
   use bsts_mod, only: bsts_prediction_t, bsts_predict, bsts_predict_draws
   use bssm_mod, only: bssm_prediction_t, bssm_predictive_sample, &
      bssm_predictive_draws
   use var_etp_mod, only: var_etp_forecast_result_t, var_etp_forecast, &
      var_etp_bootstrap_prediction
   use bvartools_mod, only: bvartools_predictive_t, &
      bvartools_bvar_predictive_from_random, bvartools_bvar_predictive, &
      bvartools_tvp_bvar_predictive_from_random, bvartools_tvp_bvar_predictive
   use tsdyn_mod, only: tsdyn_setar_rolling_forecast, &
      tsdyn_lstar_rolling_forecast, tsdyn_star_rolling_forecast, &
      tsdyn_llar_rolling_forecast, tsdyn_aar_rolling_forecast, &
      tsdyn_nnet_rolling_forecast
   use mixar_mod, only: mixar_predictive_distribution_t, &
      mixar_exact_forecast_t, mixar_forecast_sample_t, mixar_exact_forecast, &
      mixar_predictive_density, mixar_predictive_cdf, &
      mixar_predictive_quantile, mixar_forecast_from_draws, mixar_forecast
   use mixar_mod, only: mixar_regression_exact_forecast, &
      mixar_regression_forecast_from_draws, mixar_regression_forecast
   use mixar_mod, only: mixar_seasonal_exact_forecast, &
      mixar_seasonal_forecast_from_draws, mixar_seasonal_forecast
   use mixar_mod, only: mixar_var_forecast_t, &
      mixar_var_forecast_from_draws, mixar_var_forecast
   use ugmar_mod, only: ugmar_simulation_t, ugmar_forecast_t, &
      ugmar_simulate, ugmar_forecast
   use setartree_mod, only: setartree_prediction_t, setartree_forecast_t, &
      setartree_predict, setartree_predict_categorical, setartree_forecast, &
      setarforest_predict, setarforest_predict_categorical, setarforest_forecast
   use tseriestarma_mod, only: tseriestarma_forecast_t, &
      tseriestarma_forecast_from_standard, tseriestarma_forecast
   use gmdh_mod, only: gmdh_model_t, gmdh_forecast_t, gmdh_fit, &
      gmdh_predict, gmdh_forecast
   use nnfor_mod, only: nnfor_elm_model_t, nnfor_mlp_model_t, &
      nnfor_forecast_t, nnfor_hidden_selection_t, &
      nnfor_elm_hidden_selection_t, nnfor_seasonality_t, nnfor_trend_t, &
      nnfor_canova_hansen_t, nnfor_mseason_t, nnfor_difference_selection_t, &
      nnfor_elm_auto_t, nnfor_mlp_auto_t, &
      nnfor_ch_trigonometric, nnfor_ch_dummy, &
      nnfor_elm, nnfor_elm_layers, &
      nnfor_elm_forecast, nnfor_mlp, &
      nnfor_mlp_layers, nnfor_mlp_forecast, nnfor_elm_preprocessed, &
      nnfor_elm_preprocessed_forecast, nnfor_elm_preprocessed_layers, &
      nnfor_mlp_preprocessed, &
      nnfor_mlp_preprocessed_layers, nnfor_mlp_preprocessed_forecast, &
      nnfor_select_hidden_count, &
      nnfor_select_hidden_count_folds, nnfor_select_hidden_count_random, &
      nnfor_select_elm_hidden, nnfor_season_check, nnfor_mseason_test, &
      nnfor_canova_hansen, &
      nnfor_trend_check, &
      nnfor_select_differences, nnfor_elm_auto, nnfor_elm_auto_forecast, &
      nnfor_mlp_auto, nnfor_mlp_auto_forecast, &
      nnfor_elm_thief, nnfor_mlp_thief, &
      nnfor_elm_refit, nnfor_elm_retrain, nnfor_mlp_refit, nnfor_mlp_retrain
   use narfima_mod, only: narfima_model_t, narfima_forecast_t, &
      narfima_fit, narfima_auto_narfima, narfima_auto_narima, &
      narfima_auto_nbsts, narfima_auto_nnaive, narfima_forecast, &
      narfima_forecast_from_innovations
   use nlints_mod, only: nlints_varmlp_t, nlints_varmlp_fit, &
      nlints_varmlp_update, nlints_varmlp_predict, nlints_varmlp_forecast, &
      nlints_varmlp_forecast_table, nlints_varmlp_save, nlints_varmlp_load
   use tslstmplus_mod, only: tslstmplus_model_t, tslstmplus_prepared_t, &
      tslstmplus_prepare_data, tslstmplus_prepare_sequences, &
      tslstmplus_fit, tslstmplus_predict, tslstmplus_forecast
   use tsann_mod, only: tsann_fit_t, tsann_trace_t, &
      tsann_select_validation, tsann_select_test, tsann_maximum_lag, &
      tsann_auto_fit
   use echos_mod, only: echos_model_t, echos_forecast_t, echos_tuning_t, &
      echos_ridge_t, echos_ic_aic, echos_ic_aicc, echos_ic_bic, &
      echos_ic_hqc, echos_run_reservoir, echos_fit_ridge, echos_fit, &
      echos_forecast, echos_tune
   use starvars_mod, only: starvars_forecast_t, starvars_forecast_naive, &
      starvars_forecast_monte_carlo, starvars_forecast_bootstrap, &
      starvars_forecast
   use rugarch_mod, only: rugarch_fit_t, rugarch_forecast_t, rugarch_forecast
   implicit none
   private

   public :: forecast_result_t, meanf, naive, rwf, snaive, ses
   public :: holt, holt_winters, croston, thetaf
   public :: itsmr_forecast_t, arma_forecast, arar_forecast
   public :: transformed_arma_forecast
   public :: astsa_sarima_forecast_t, sarima_forecast
   public :: arfima_forecast_t, arfima_forecast
   public :: arfima_average_forecast, arfima_forecast_paths
   public :: garma_forecast_t, garma_forecast, garma_regression_forecast
   public :: tfarima_forecast_t, tfarima_arima_forecast
   public :: tfarima_exact_transfer_forecast, tfarima_structural_forecast
   public :: robustarima_forecast_t, robustarima_forecast
   public :: expar_forecast_t, expar_forecast
   public :: mts_var_forecast_t, mts_var_forecast, mts_varx_forecast
   public :: mts_varma_forecast, mts_vecm_forecast, mts_factor_forecast
   public :: bigtime_forecast_t, bigtime_var_forecast
   public :: bigtime_varx_forecast, bigtime_varma_forecast
   public :: bigvar_forecast, bigvar_direct_forecast
   public :: bigvar_varx_forecast, bigvar_var_forecast_interval
   public :: bsts_prediction_t, bsts_predict, bsts_predict_draws
   public :: bssm_prediction_t, bssm_predictive_sample, bssm_predictive_draws
   public :: var_etp_forecast_result_t, var_etp_forecast
   public :: var_etp_bootstrap_prediction
   public :: bvartools_predictive_t
   public :: bvartools_bvar_predictive_from_random, bvartools_bvar_predictive
   public :: bvartools_tvp_bvar_predictive_from_random
   public :: bvartools_tvp_bvar_predictive
   public :: rolling_forecast_result_t
   public :: accuracy_result_t, rolling_forecast_accuracy
   public :: tsdyn_setar_rolling_forecast, tsdyn_lstar_rolling_forecast
   public :: tsdyn_star_rolling_forecast, tsdyn_llar_rolling_forecast
   public :: tsdyn_aar_rolling_forecast, tsdyn_nnet_rolling_forecast
   public :: mixar_predictive_distribution_t, mixar_exact_forecast_t
   public :: mixar_forecast_sample_t, mixar_exact_forecast
   public :: mixar_predictive_density, mixar_predictive_cdf
   public :: mixar_predictive_quantile, mixar_forecast_from_draws
   public :: mixar_forecast
   public :: mixar_regression_exact_forecast
   public :: mixar_regression_forecast_from_draws, mixar_regression_forecast
   public :: mixar_seasonal_exact_forecast
   public :: mixar_seasonal_forecast_from_draws, mixar_seasonal_forecast
   public :: mixar_var_forecast_t, mixar_var_forecast_from_draws
   public :: mixar_var_forecast
   public :: ugmar_simulation_t, ugmar_forecast_t
   public :: ugmar_simulate, ugmar_forecast
   public :: setartree_prediction_t, setartree_forecast_t
   public :: real_vector_t
   public :: setartree_predict, setartree_predict_categorical, setartree_forecast
   public :: setarforest_predict, setarforest_predict_categorical
   public :: setarforest_forecast
   public :: tseriestarma_forecast_t
   public :: tseriestarma_forecast_from_standard, tseriestarma_forecast
   public :: gmdh_model_t, gmdh_forecast_t
   public :: gmdh_fit, gmdh_predict, gmdh_forecast
   public :: nnfor_elm_model_t, nnfor_mlp_model_t, nnfor_forecast_t
   public :: nnfor_hidden_selection_t, nnfor_seasonality_t, nnfor_trend_t
   public :: nnfor_elm_hidden_selection_t
   public :: nnfor_canova_hansen_t
   public :: nnfor_mseason_t
   public :: nnfor_difference_selection_t
   public :: nnfor_ch_trigonometric, nnfor_ch_dummy
   public :: nnfor_elm_auto_t, nnfor_mlp_auto_t
   public :: nnfor_elm, nnfor_elm_forecast
   public :: nnfor_elm_layers
   public :: nnfor_mlp, nnfor_mlp_forecast
   public :: nnfor_mlp_layers
   public :: nnfor_elm_preprocessed, nnfor_elm_preprocessed_forecast
   public :: nnfor_elm_preprocessed_layers
   public :: nnfor_mlp_preprocessed, nnfor_mlp_preprocessed_forecast
   public :: nnfor_mlp_preprocessed_layers
   public :: nnfor_select_hidden_count, nnfor_select_hidden_count_folds
   public :: nnfor_select_hidden_count_random
   public :: nnfor_select_elm_hidden
   public :: nnfor_season_check, nnfor_mseason_test, nnfor_canova_hansen
   public :: nnfor_trend_check, nnfor_select_differences
   public :: nnfor_elm_auto, nnfor_elm_auto_forecast
   public :: nnfor_mlp_auto, nnfor_mlp_auto_forecast
   public :: nnfor_elm_thief, nnfor_mlp_thief
   public :: nnfor_elm_refit, nnfor_elm_retrain
   public :: nnfor_mlp_refit, nnfor_mlp_retrain
   public :: narfima_model_t, narfima_forecast_t, narfima_fit
   public :: narfima_auto_narfima, narfima_auto_narima
   public :: narfima_auto_nbsts, narfima_auto_nnaive
   public :: narfima_forecast, narfima_forecast_from_innovations
   public :: nlints_varmlp_t, nlints_varmlp_fit, nlints_varmlp_update
   public :: nlints_varmlp_predict, nlints_varmlp_forecast
   public :: nlints_varmlp_forecast_table
   public :: nlints_varmlp_save, nlints_varmlp_load
   public :: tslstmplus_model_t, tslstmplus_prepared_t
   public :: tslstmplus_prepare_data, tslstmplus_prepare_sequences
   public :: tslstmplus_fit, tslstmplus_predict, tslstmplus_forecast
   public :: tsann_fit_t, tsann_trace_t
   public :: tsann_select_validation, tsann_select_test
   public :: tsann_maximum_lag, tsann_auto_fit
   public :: echos_model_t, echos_forecast_t, echos_tuning_t, echos_ridge_t
   public :: echos_ic_aic, echos_ic_aicc, echos_ic_bic, echos_ic_hqc
   public :: echos_run_reservoir, echos_fit_ridge, echos_fit
   public :: echos_forecast, echos_tune
   public :: starvars_forecast_t, starvars_forecast_naive
   public :: starvars_forecast_monte_carlo, starvars_forecast_bootstrap
   public :: starvars_forecast
   public :: rugarch_fit_t, rugarch_forecast_t, rugarch_forecast
end module forecasting_mod
