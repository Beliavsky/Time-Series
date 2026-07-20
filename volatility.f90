! SPDX-License-Identifier: MIT
! SPDX-FileComment: Topic facade for conditional and stochastic volatility.
module volatility_mod
   !! Re-export ARCH, multivariate GARCH, and stochastic-volatility algorithms.
   use astsa_mod, only: astsa_sv_filter_t, astsa_sv_fit_t, astsa_sv_mcmc_t, &
      sv_filter, sv_mle, sv_particle_filter, sv_particle_filter_draws, sv_mcmc
   use bayesforecast_mod, only: bf_sv_filter_t, bf_garch_filter, &
      bf_stochastic_volatility_filter, bf_garch_predict, &
      bf_stochastic_volatility_predict
   use mts_mod, only: mts_bekk_fit_t, mts_dcc_fit_t, mts_adcc_fit_t, &
      mts_bekk_fit, mts_bekk_forecast, mts_dcc_fit, mts_dcc_forecast, &
      mts_adcc_fit, mts_adcc_t_fit, mts_adcc_forecast, mts_arch_test, &
      mts_march_test, mts_common_volatility
   use tsissm_mod, only: tsissm_garch_t, tsissm_garch_recursion, &
      tsissm_initialize_variance
   use bvartools_mod, only: bvartools_stochastic_volatility_t, &
      bvartools_stochastic_volatility_draw, bvartools_stochastic_volatility, &
      bvartools_stochastic_volatility_ocsn_draw, &
      bvartools_stochastic_volatility_ocsn
   use tseriestarma_mod, only: tseriestarma_arma_garch_fit_t, &
      tseriestarma_arma_garch_fit
   use starvars_mod, only: starvars_realized_covariance_t, &
      starvars_frequency_daily, starvars_frequency_monthly, &
      starvars_frequency_quarterly, starvars_frequency_yearly, &
      starvars_realized_covariance
   use rugarch_mod, only: rugarch_spec_t, rugarch_parameters_t, &
      rugarch_filter_t, rugarch_fit_t, rugarch_simulation_t, &
      rugarch_model_sgarch, rugarch_model_igarch, rugarch_model_egarch, &
      rugarch_model_gjrgarch, rugarch_model_aparch, rugarch_model_figarch, &
      rugarch_model_csgarch, rugarch_model_realgarch, &
      rugarch_model_fgarch, rugarch_fgarch_garch, rugarch_fgarch_tgarch, &
      rugarch_fgarch_avgarch, rugarch_fgarch_ngarch, &
      rugarch_fgarch_nagarch, rugarch_fgarch_aparch, &
      rugarch_fgarch_allgarch, rugarch_fgarch_gjrgarch, &
      rugarch_distribution_normal, rugarch_distribution_student, &
      rugarch_distribution_ged, rugarch_distribution_skew_normal, &
      rugarch_distribution_skew_student, rugarch_distribution_skew_ged, &
      rugarch_distribution_johnson_su, rugarch_distribution_nig, &
      rugarch_distribution_ghyp, rugarch_distribution_gh_skew_student, &
      rugarch_spec, rugarch_filter, rugarch_fit, &
      rugarch_simulate, rugarch_coefficients, rugarch_persistence, &
      rugarch_unconditional_variance, rugarch_half_life, rugarch_log_density
   use rugarch_extensions_mod, only: rugarch_bootstrap_t, rugarch_roll_t, &
      rugarch_multifit_t, rugarch_parameter_distribution_t, &
      rugarch_multifilter_t, rugarch_multiforecast_t, &
      rugarch_model_confidence_t, rugarch_news_impact_t, &
      rugarch_sign_bias_test_t, rugarch_bootstrap_forecast, rugarch_roll, &
      rugarch_multifit, rugarch_parameter_distribution, &
      rugarch_multifilter, rugarch_multiforecast, &
      rugarch_model_confidence, rugarch_news_impact, rugarch_sign_bias_test
   use distribution_fit_mod, only: distribution_fit_t, fit_distribution
   use distribution_mod, only: distribution_log_density, &
      distribution_density, distribution_cdf, distribution_quantile, &
      random_distribution, standardized_skewness, &
      standardized_excess_kurtosis
   implicit none
   private

   public :: astsa_sv_filter_t, astsa_sv_fit_t, astsa_sv_mcmc_t
   public :: sv_filter, sv_mle, sv_particle_filter, sv_particle_filter_draws
   public :: sv_mcmc
   public :: bf_sv_filter_t, bf_garch_filter, bf_stochastic_volatility_filter
   public :: bf_garch_predict, bf_stochastic_volatility_predict
   public :: mts_bekk_fit_t, mts_dcc_fit_t, mts_adcc_fit_t
   public :: mts_bekk_fit, mts_bekk_forecast, mts_dcc_fit, mts_dcc_forecast
   public :: mts_adcc_fit, mts_adcc_t_fit, mts_adcc_forecast
   public :: mts_arch_test, mts_march_test, mts_common_volatility
   public :: tsissm_garch_t, tsissm_garch_recursion, tsissm_initialize_variance
   public :: bvartools_stochastic_volatility_t
   public :: bvartools_stochastic_volatility_draw
   public :: bvartools_stochastic_volatility
   public :: bvartools_stochastic_volatility_ocsn_draw
   public :: bvartools_stochastic_volatility_ocsn
   public :: tseriestarma_arma_garch_fit_t, tseriestarma_arma_garch_fit
   public :: starvars_realized_covariance_t, starvars_realized_covariance
   public :: starvars_frequency_daily, starvars_frequency_monthly
   public :: starvars_frequency_quarterly, starvars_frequency_yearly
   public :: rugarch_spec_t, rugarch_parameters_t, rugarch_filter_t
   public :: rugarch_fit_t, rugarch_simulation_t
   public :: rugarch_model_sgarch, rugarch_model_igarch, rugarch_model_egarch
   public :: rugarch_model_gjrgarch, rugarch_model_aparch
   public :: rugarch_model_figarch
   public :: rugarch_model_csgarch, rugarch_model_realgarch
   public :: rugarch_model_fgarch
   public :: rugarch_fgarch_garch, rugarch_fgarch_tgarch
   public :: rugarch_fgarch_avgarch, rugarch_fgarch_ngarch
   public :: rugarch_fgarch_nagarch, rugarch_fgarch_aparch
   public :: rugarch_fgarch_allgarch, rugarch_fgarch_gjrgarch
   public :: rugarch_distribution_normal, rugarch_distribution_student
   public :: rugarch_distribution_ged
   public :: rugarch_distribution_skew_normal, rugarch_distribution_skew_student
   public :: rugarch_distribution_skew_ged, rugarch_distribution_johnson_su
   public :: rugarch_distribution_nig, rugarch_distribution_ghyp
   public :: rugarch_distribution_gh_skew_student
   public :: rugarch_spec, rugarch_filter, rugarch_fit, rugarch_simulate
   public :: rugarch_coefficients, rugarch_persistence
   public :: rugarch_unconditional_variance, rugarch_half_life
   public :: rugarch_log_density
   public :: rugarch_bootstrap_t, rugarch_roll_t, rugarch_multifit_t
   public :: rugarch_multifilter_t, rugarch_multiforecast_t
   public :: rugarch_parameter_distribution_t, rugarch_model_confidence_t
   public :: rugarch_news_impact_t, rugarch_sign_bias_test_t
   public :: rugarch_bootstrap_forecast, rugarch_roll, rugarch_multifit
   public :: rugarch_multifilter, rugarch_multiforecast
   public :: rugarch_parameter_distribution, rugarch_model_confidence
   public :: rugarch_news_impact, rugarch_sign_bias_test
   public :: distribution_fit_t, fit_distribution
   public :: distribution_log_density, distribution_density
   public :: distribution_cdf, distribution_quantile, random_distribution
   public :: standardized_skewness, standardized_excess_kurtosis
end module volatility_mod
