! SPDX-License-Identifier: MIT
! SPDX-FileComment: Compilation test for topic facade modules.
program test_facades
   !! Verify that every topic facade and a representative export compile.
   use arma_mod, only: arma_acvf, tsdyn_lstar_fit, tsdyn_aar_fit, &
      tsdyn_setar_fit, tsdyn_setar_linearity_test
   use arma_mod, only: tsdyn_nnet_fit, tsdyn_setar_bootstrap, &
      tsdyn_setar_forecast_distribution, tsdyn_lstar_forecast_distribution, &
      tsdyn_ar_regime_diagnostics
   use long_memory_mod, only: fracdiff_difference, narfima_auto_narfima
   use state_space_mod, only: make_local_level, nts_smc_filter
   use multivariate_mod, only: johansen_test, vars_fit, tsdyn_tvecm_fit, &
      tsdyn_tvar_fit, tsdyn_hansen_seo_test, tsdyn_nonlinear_fevd
   use multivariate_mod, only: tsdyn_tvar_lr_test, tsdyn_tvar_bootstrap, &
      tsdyn_tvecm_bootstrap, tsdyn_rank_select
   use multivariate_mod, only: tsdyn_tvar_inference, tsdyn_tvecm_inference
   use multivariate_mod, only: tsdyn_tvar_regimes, tsdyn_tvecm_regimes
   use multivariate_mod, only: starvars_fit
   use volatility_mod, only: sv_filter
   use volatility_mod, only: starvars_realized_covariance
   use volatility_mod, only: rugarch_spec, rugarch_fit
   use spectral_mod, only: real_dft
   use diagnostics_mod, only: weighted_box_test, vars_normality_test, &
      tsdyn_bbc_test, tsdyn_kapshin_test, tsdyn_delta, tsdyn_delta_test, &
      tsdyn_delta_linear, tsdyn_delta_linear_test
   use diagnostics_mod, only: tsdyn_star_regime_diagnostics
   use diagnostics_mod, only: tsdyn_rank_test
   use diagnostics_mod, only: tsdyn_setar_residual_variance
   use diagnostics_mod, only: echos_kpss
   use diagnostics_mod, only: starvars_joint_linearity_test
   use diagnostics_mod, only: rugarch_var_test
   use arma_mod, only: tsdyn_setar_inference
   use arma_mod, only: tsdyn_setar_regimes, tsdyn_lstar_regimes, &
      tsdyn_star_regimes
   use forecasting_mod, only: meanf, tsdyn_setar_rolling_forecast, &
      tsdyn_lstar_rolling_forecast, rolling_forecast_accuracy, gmdh_forecast, &
      nnfor_elm_forecast
   use forecasting_mod, only: narfima_forecast
   use forecasting_mod, only: tsann_auto_fit
   use forecasting_mod, only: echos_fit, echos_forecast
   use forecasting_mod, only: starvars_forecast
   use forecasting_mod, only: rugarch_forecast
   use regression_time_series_mod, only: lag_reg
   use bayesian_time_series_mod, only: bf_log_prior
   use markov_switching_mod, only: mswm_gaussian_fit
   use count_time_series_mod, only: nts_acmx_fit
   use functional_time_series_mod, only: nts_cfar_estimate
   implicit none

   print '(a)', 'All topic facade modules compiled successfully.'
end program test_facades
