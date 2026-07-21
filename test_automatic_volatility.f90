! SPDX-License-Identifier: MIT
! SPDX-FileComment: Tests for automatic conditional-volatility model selection.
program test_automatic_volatility
   !! Exercise volatility candidate fitting, ranking, and sigma forecasts.
   use kind_mod, only: dp
   use automatic_volatility_mod, only: automatic_volatility_options_t, &
      automatic_volatility_result_t, automatic_volatility
   implicit none
   type(automatic_volatility_options_t) :: options
   type(automatic_volatility_result_t) :: result
   real(dp) :: series(80), scale
   integer :: i

   scale = 0.2_dp
   do i = 1, size(series)
      scale = 0.15_dp + 0.75_dp*scale + &
         0.12_dp*merge(1.0_dp, 0.0_dp, modulo(i, 9) == 0)
      series(i) = sqrt(scale)*sin(1.73_dp*real(i, dp))
   end do
   options%selection = "bic"
   options%horizon = 4
   result = automatic_volatility(series, options)
   if (result%info /= 0) error stop "automatic volatility status"
   if (size(result%candidates) /= 4) then
      error stop "automatic volatility candidate count"
   end if
   if (result%selected /= 1 .or. &
      .not. result%candidates(1)%selected) then
      error stop "automatic volatility selection"
   end if
   if (size(result%sigma_forecast) /= options%horizon) then
      error stop "automatic volatility horizon"
   end if
   if (any(result%sigma_forecast <= 0.0_dp) .or. &
      any(result%sigma_forecast /= result%sigma_forecast)) then
      error stop "automatic volatility forecast"
   end if
   do i = 2, size(result%candidates)
      if (result%candidates(i)%bic < result%candidates(i - 1)%bic) then
         error stop "automatic volatility ranking"
      end if
   end do

   print '(a)', "automatic volatility tests passed"
end program test_automatic_volatility
