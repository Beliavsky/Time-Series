! SPDX-License-Identifier: MIT
! SPDX-FileComment: Topic facade for Markov-switching models.
module markov_switching_mod
   !! Re-export filtering, smoothing, estimation, and diagnostics for switching models.
   use mswm_mod, only: mswm_filter_t, mswm_fit_t, mswm_glm_fit_t, &
      mswm_ar_data_t, mswm_inference_t, mswm_intervals_t, &
      mswm_multistart_t, mswm_glm_multistart_t, &
      mswm_gaussian_filter, mswm_gaussian_fit, mswm_glm_filter, &
      mswm_glm_fit, mswm_residuals, mswm_states, &
      mswm_ar_data, mswm_gaussian_ar_fit, mswm_glm_ar_fit, &
      mswm_gaussian_inference, mswm_gaussian_intervals, &
      mswm_glm_inference, mswm_glm_intervals, &
      mswm_gaussian_multistart, mswm_gaussian_random_starts, &
      mswm_glm_multistart, mswm_glm_random_starts, display
   use nts_mod, only: nts_msar_simulation_t, nts_msar_simulate_from_draws, &
      nts_msar_simulate
   implicit none
   private

   public :: mswm_filter_t, mswm_fit_t, mswm_glm_fit_t
   public :: mswm_ar_data_t, mswm_inference_t, mswm_intervals_t
   public :: mswm_multistart_t, mswm_glm_multistart_t
   public :: mswm_gaussian_filter, mswm_gaussian_fit
   public :: mswm_glm_filter, mswm_glm_fit
   public :: mswm_ar_data, mswm_gaussian_ar_fit, mswm_glm_ar_fit
   public :: mswm_residuals, mswm_states
   public :: mswm_gaussian_inference, mswm_gaussian_intervals
   public :: mswm_glm_inference, mswm_glm_intervals
   public :: mswm_gaussian_multistart, mswm_gaussian_random_starts
   public :: mswm_glm_multistart, mswm_glm_random_starts
   public :: display
   public :: nts_msar_simulation_t, nts_msar_simulate_from_draws
   public :: nts_msar_simulate
end module markov_switching_mod
