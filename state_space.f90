! SPDX-License-Identifier: MIT
! SPDX-FileComment: Topic facade for state-space models.
module state_space_mod
   !! Re-export linear, nonlinear, Gaussian, and particle state-space algorithms.
   use kfas_mod, only: ssm_model_t, kfs_filter_t, kfs_smoother_t, &
      kfs_prediction_t, kfs_filter, kfs_filter_diffuse, kfs_smooth, kfs_loglik, &
      kfs_predict, kfs_fast_smooth, kfs_standardized_innovations, &
      kfs_disturbance_smooth, make_local_level, make_local_linear_trend
   use astsa_mod, only: astsa_filter_t, astsa_smoother_t, astsa_em_result_t, &
      astsa_kfilter, astsa_ksmooth, astsa_kfilter_correlated, &
      astsa_ksmooth_correlated, astsa_em, astsa_ffbs, astsa_ffbs_draws, ssm_fit
   use bssm_mod, only: bssm_particle_filter_t, bssm_particle_smoother_t, &
      bssm_ekf_t, bssm_ekf_smoother_t, bssm_bootstrap_filter, &
      bssm_nonlinear_bootstrap_filter, bssm_iekf, bssm_ekf_smoother, &
      bssm_ukf, bssm_particle_smoother, bssm_simulation_smoother
   use tsss_mod, only: tsss_particle_filter_t, tsss_structural_model_t, &
      tsss_particle_filter, tsss_nonlinear_particle_filter, &
      tsss_structural_model, tsss_tsmooth, tsss_simssm, tsss_ngsim
   use tfarima_mod, only: tfarima_structural_ssm_t, &
      tfarima_structural_filter_t, tfarima_structural_smoother_t, &
      tfarima_arima_to_structural_ssm, tfarima_structural_initialize, &
      tfarima_structural_filter, tfarima_structural_smooth, &
      tfarima_structural_forecast
   use bvartools_mod, only: bvartools_kalman_dk_t, bvartools_kalman_dk
   use nts_mod, only: nts_tvar_fit_t, nts_tvar_filter_smooth, nts_tvar_fit, &
      nts_smc_filter_t, nts_smc_smoother_t, nts_smc_marginal_smoother_t, &
      nts_rb_smc_filter_t, nts_smc_filter_draws, nts_smc_filter, &
      nts_smc_smooth, nts_smc_marginal_smooth, &
      nts_rb_smc_filter_draws, nts_rb_smc_filter
   implicit none
   private

   public :: ssm_model_t, kfs_filter_t, kfs_smoother_t, kfs_prediction_t
   public :: kfs_filter, kfs_filter_diffuse, kfs_smooth, kfs_loglik, kfs_predict
   public :: kfs_fast_smooth, kfs_standardized_innovations
   public :: kfs_disturbance_smooth, make_local_level, make_local_linear_trend
   public :: astsa_filter_t, astsa_smoother_t, astsa_em_result_t
   public :: astsa_kfilter, astsa_ksmooth, astsa_kfilter_correlated
   public :: astsa_ksmooth_correlated, astsa_em, astsa_ffbs, astsa_ffbs_draws
   public :: ssm_fit
   public :: bssm_particle_filter_t, bssm_particle_smoother_t
   public :: bssm_ekf_t, bssm_ekf_smoother_t, bssm_bootstrap_filter
   public :: bssm_nonlinear_bootstrap_filter, bssm_iekf, bssm_ekf_smoother
   public :: bssm_ukf, bssm_particle_smoother, bssm_simulation_smoother
   public :: tsss_particle_filter_t, tsss_structural_model_t
   public :: tsss_particle_filter, tsss_nonlinear_particle_filter
   public :: tsss_structural_model, tsss_tsmooth, tsss_simssm, tsss_ngsim
   public :: tfarima_structural_ssm_t, tfarima_structural_filter_t
   public :: tfarima_structural_smoother_t, tfarima_arima_to_structural_ssm
   public :: tfarima_structural_initialize, tfarima_structural_filter
   public :: tfarima_structural_smooth, tfarima_structural_forecast
   public :: bvartools_kalman_dk_t, bvartools_kalman_dk
   public :: nts_tvar_fit_t, nts_tvar_filter_smooth, nts_tvar_fit
   public :: nts_smc_filter_t, nts_smc_smoother_t
   public :: nts_smc_marginal_smoother_t, nts_rb_smc_filter_t
   public :: nts_smc_filter_draws, nts_smc_filter, nts_smc_smooth
   public :: nts_smc_marginal_smooth
   public :: nts_rb_smc_filter_draws, nts_rb_smc_filter
end module state_space_mod
