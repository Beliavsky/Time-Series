! SPDX-License-Identifier: MIT
! SPDX-FileComment: Topic facade for regression with time-series errors.
module regression_time_series_mod
   !! Re-export dynamic regression, transfer-function, and intervention algorithms.
   use utils_mod, only: real_vector_t
   use astsa_mod, only: astsa_lag_regression_t, astsa_stochastic_regression_t, &
      lag_reg, stochastic_regression, pre_white
   use arfima_mod, only: arfima_regression_fit_t, arfima_regression_fit, &
      arfima_regression_fit_fixed, arfima_regression_forecast, &
      arfima_transfer, arfima_transfer_response
   use tfarima_mod, only: tfarima_transfer_fit_t, tfarima_transfer_spec_t, &
      tfarima_transfer, tfarima_transfer_fit, tfarima_exact_transfer_fit, &
      tfarima_select_transfer, tfarima_identify_transfer, tfarima_intervention, &
      tfarima_calendar_regressors, tfarima_seasonal_dummies, &
      tfarima_harmonic_regressors
   use mts_mod, only: mts_transfer_fit_t, mts_regts_fit_t, mts_transfer_fit, &
      mts_transfer2_fit, mts_transfer2_forecast, mts_regts_fit, &
      mts_regts_refine, mts_regts_forecast, mts_multivariate_regression
   use bsts_mod, only: bsts_dynamic_regression_t, bsts_dynamic_regression, &
      bsts_dynamic_regression_predict, bsts_dynamic_regression_ar, &
      bsts_regression_holiday, bsts_regression_holiday_predict
   use var_etp_mod, only: var_etp_predictive_result_t, var_etp_order_result_t, &
      var_etp_predictive_regression, var_etp_predictive_order, &
      var_etp_predictive_forecast
   use bentcablear_mod, only: bentcable_residuals_t, bentcable_fit_t, &
      bentcable_profile_t, bentcable_change_t, bentcable_basis, &
      bentcable_value, bentcable_design_matrix, bentcable_residuals, &
      bentcable_sse, bentcable_stationary, bentcable_profile, bentcable_fit, &
      bentcable_fisher_information, bentcable_change_interval, &
      bentcable_ar_covariance_t, bentcable_ar_covariance, &
      bentcable_fit_iterative_yw
   use mixar_mod, only: mixar_regression_model_t, mixar_regression_fit_t, &
      mixar_regression_simulation_t, mixar_regression_model, &
      mixar_regression_filter, mixar_regression_log_likelihood, &
      mixar_regression_fit, mixar_regression_simulate_from_draws, &
      mixar_regression_simulate, mixar_regression_exact_forecast, &
      mixar_regression_forecast_from_draws, mixar_regression_forecast, &
      mixar_inference_t, mixar_regression_observed_inference
   use setartree_mod, only: setartree_category_levels_t, setartree_model_t, &
      setarforest_model_t, setartree_prediction_t, setartree_fit, &
      setartree_fit_categorical, setartree_predict, &
      setartree_predict_categorical, setarforest_fit, &
      setarforest_fit_categorical, setarforest_predict, &
      setarforest_predict_categorical
   use tseriestarma_mod, only: tseriestarma_fit2_t, tseriestarma_fit2
   use tslstmplus_mod, only: tslstmplus_model_t, tslstmplus_prepared_t, &
      tslstmplus_prepare_data, tslstmplus_prepare_sequences, &
      tslstmplus_fit, tslstmplus_predict, tslstmplus_forecast
   implicit none
   private

   public :: astsa_lag_regression_t, astsa_stochastic_regression_t
   public :: lag_reg, stochastic_regression, pre_white
   public :: arfima_regression_fit_t, arfima_regression_fit
   public :: arfima_regression_fit_fixed, arfima_regression_forecast
   public :: arfima_transfer, arfima_transfer_response
   public :: tfarima_transfer_fit_t, tfarima_transfer_spec_t
   public :: tfarima_transfer, tfarima_transfer_fit, tfarima_exact_transfer_fit
   public :: tfarima_select_transfer, tfarima_identify_transfer
   public :: tfarima_intervention, tfarima_calendar_regressors
   public :: tfarima_seasonal_dummies, tfarima_harmonic_regressors
   public :: mts_transfer_fit_t, mts_regts_fit_t, mts_transfer_fit
   public :: mts_transfer2_fit, mts_transfer2_forecast
   public :: mts_regts_fit, mts_regts_refine, mts_regts_forecast
   public :: mts_multivariate_regression
   public :: bsts_dynamic_regression_t, bsts_dynamic_regression
   public :: bsts_dynamic_regression_predict, bsts_dynamic_regression_ar
   public :: bsts_regression_holiday, bsts_regression_holiday_predict
   public :: var_etp_predictive_result_t, var_etp_order_result_t
   public :: var_etp_predictive_regression, var_etp_predictive_order
   public :: var_etp_predictive_forecast
   public :: bentcable_residuals_t, bentcable_fit_t
   public :: bentcable_profile_t, bentcable_change_t
   public :: bentcable_basis, bentcable_value, bentcable_design_matrix
   public :: bentcable_residuals, bentcable_sse, bentcable_stationary
   public :: bentcable_profile, bentcable_fit
   public :: bentcable_fisher_information, bentcable_change_interval
   public :: bentcable_ar_covariance_t, bentcable_ar_covariance
   public :: bentcable_fit_iterative_yw
   public :: mixar_regression_model_t, mixar_regression_fit_t
   public :: mixar_regression_simulation_t, mixar_regression_model
   public :: mixar_regression_filter, mixar_regression_log_likelihood
   public :: mixar_regression_fit, mixar_regression_simulate_from_draws
   public :: mixar_regression_simulate, mixar_regression_exact_forecast
   public :: mixar_regression_forecast_from_draws, mixar_regression_forecast
   public :: mixar_inference_t, mixar_regression_observed_inference
   public :: setartree_category_levels_t, setartree_model_t
   public :: setarforest_model_t, setartree_prediction_t
   public :: real_vector_t
   public :: setartree_fit, setartree_fit_categorical
   public :: setartree_predict, setartree_predict_categorical
   public :: setarforest_fit, setarforest_fit_categorical
   public :: setarforest_predict, setarforest_predict_categorical
   public :: tseriestarma_fit2_t, tseriestarma_fit2
   public :: tslstmplus_model_t, tslstmplus_prepared_t
   public :: tslstmplus_prepare_data, tslstmplus_prepare_sequences
   public :: tslstmplus_fit, tslstmplus_predict, tslstmplus_forecast
end module regression_time_series_mod
