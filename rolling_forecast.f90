! SPDX-License-Identifier: MIT
! SPDX-FileComment: Generic rolling-origin forecast result types.
module rolling_forecast_mod
   !! Define shared containers for rolling-origin forecast evaluation.
   use kind_mod, only: dp
   implicit none
   private

   type, public :: rolling_forecast_result_t
      !! Forecasts and errors aligned by evaluation target and horizon.
      real(dp), allocatable :: forecast(:, :)
      real(dp), allocatable :: actual(:, :)
      real(dp), allocatable :: error(:, :)
      integer, allocatable :: origin(:, :)
      integer, allocatable :: target(:)
      integer, allocatable :: horizon(:)
      logical, allocatable :: valid(:, :)
      logical, allocatable :: refitted(:, :)
      integer :: training_size = 0
      integer :: info = 0
   end type rolling_forecast_result_t

end module rolling_forecast_mod
