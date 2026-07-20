! SPDX-License-Identifier: MIT
! SPDX-FileComment: Topic facade for functional time-series models.
module functional_time_series_mod
   !! Re-export continuous functional autoregression algorithms and result types.
   use nts_mod, only: nts_cfar_simulation_t, nts_cfar_irregular_simulation_t, &
      nts_cfar_model_t, &
      nts_cfar_order_tests_t, nts_cfar_forecast_t, &
      nts_cfar_simulate_from_standard, nts_cfar_simulate, &
      nts_cfar1_simulate, nts_cfar2_simulate, nts_cfar_estimate, &
      nts_cfar_order_tests, nts_cfar_forecast, nts_cfar_partial_forecast, &
      nts_cfar_irregular_simulate, nts_cfar2_irregular_simulate, &
      nts_cfar_irregular_estimate, nts_cfar_irregular_order_tests
   implicit none
   private

   public :: nts_cfar_simulation_t, nts_cfar_model_t
   public :: nts_cfar_irregular_simulation_t
   public :: nts_cfar_order_tests_t, nts_cfar_forecast_t
   public :: nts_cfar_simulate_from_standard, nts_cfar_simulate
   public :: nts_cfar1_simulate, nts_cfar2_simulate, nts_cfar_estimate
   public :: nts_cfar_order_tests, nts_cfar_forecast, nts_cfar_partial_forecast
   public :: nts_cfar_irregular_simulate, nts_cfar2_irregular_simulate
   public :: nts_cfar_irregular_estimate, nts_cfar_irregular_order_tests
end module functional_time_series_mod
