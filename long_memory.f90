! SPDX-License-Identifier: MIT
! SPDX-FileComment: Topic facade for fractional and long-memory models.
module long_memory_mod
   !! Re-export principal fractional differencing and long-memory algorithms.
   use arfima_mod, only: arfima_model_t, arfima_fit_t, arfima_forecast_t, &
      arfima_diagnostics_t, arfima_model, arfima_fractional_weights, &
      arfima_fdwn_acvf, arfima_fgn_acvf, arfima_hd_acvf, arfima_farma_acvf, &
      arfima_model_acvf, arfima_fit, arfima_fit_fixed, arfima_forecast, &
      arfima_simulate, arfima_simulate_from_innovations, arfima_diagnostics
   use fracdiff_mod, only: fracdiff_semiparametric_t, fracdiff_fit_t, &
      fracdiff_simulation_t, fracdiff_difference, fracdiff_gph, fracdiff_sperio, &
      fracdiff_hr_filter, fracdiff_fit, fracdiff_simulate, &
      fracdiff_simulate_from_innovations
   use fracdist_mod, only: fracdist_probability_t, fracdist_critical_values_t, &
      fracdist_p_value, fracdist_critical_values
   use nsarfima_mod, only: nsarfima_fit_t, nsarfima_simulation_t, &
      nsarfima_residuals, nsarfima_residual_acf, nsarfima_mde, nsarfima_pml, &
      nsarfima_simulate, nsarfima_simulate_from_innovations
   use esemifar_mod, only: esemifar_model_t, esemifar_forecast_t, &
      esemifar_farima_to_ar, esemifar_farima_to_ma, esemifar_order_selection, &
      esemifar_forecast_normal, esemifar_forecast_bootstrap
   use narfima_mod, only: narfima_model_t, narfima_forecast_t, narfima_fit, &
      narfima_auto_narfima, narfima_forecast, &
      narfima_forecast_from_innovations
   implicit none
   private

   public :: arfima_model_t, arfima_fit_t, arfima_forecast_t, arfima_diagnostics_t
   public :: arfima_model, arfima_fractional_weights
   public :: arfima_fdwn_acvf, arfima_fgn_acvf, arfima_hd_acvf
   public :: arfima_farma_acvf, arfima_model_acvf
   public :: arfima_fit, arfima_fit_fixed, arfima_forecast
   public :: arfima_simulate, arfima_simulate_from_innovations, arfima_diagnostics
   public :: fracdiff_semiparametric_t, fracdiff_fit_t, fracdiff_simulation_t
   public :: fracdiff_difference, fracdiff_gph, fracdiff_sperio, fracdiff_hr_filter
   public :: fracdiff_fit, fracdiff_simulate, fracdiff_simulate_from_innovations
   public :: fracdist_probability_t, fracdist_critical_values_t
   public :: fracdist_p_value, fracdist_critical_values
   public :: nsarfima_fit_t, nsarfima_simulation_t
   public :: nsarfima_residuals, nsarfima_residual_acf, nsarfima_mde, nsarfima_pml
   public :: nsarfima_simulate, nsarfima_simulate_from_innovations
   public :: esemifar_model_t, esemifar_forecast_t
   public :: esemifar_farima_to_ar, esemifar_farima_to_ma, esemifar_order_selection
   public :: esemifar_forecast_normal, esemifar_forecast_bootstrap
   public :: narfima_model_t, narfima_forecast_t
   public :: narfima_fit, narfima_auto_narfima
   public :: narfima_forecast, narfima_forecast_from_innovations
end module long_memory_mod
