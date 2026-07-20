! SPDX-License-Identifier: MIT
! SPDX-FileComment: Topic facade for discrete-valued time-series models.
module count_time_series_mod
   !! Re-export conditional count time-series models and their result types.
   use nts_mod, only: nts_acmx_fit_t, nts_acmx_fit
   implicit none
   private

   public :: nts_acmx_fit_t, nts_acmx_fit
end module count_time_series_mod
