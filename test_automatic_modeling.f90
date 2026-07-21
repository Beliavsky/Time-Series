! SPDX-License-Identifier: MIT
! SPDX-FileComment: Tests for automatic univariate time-series modeling.
program test_automatic_modeling
   !! Exercise profiling, candidate gating, ranking, and forecasts.
   use kind_mod, only: dp
   use automatic_modeling_mod, only: automatic_model_options_t, &
      automatic_model_result_t, series_profile_t, profile_series, automatic_model
   implicit none
   real(dp) :: autoregressive(120), moderate(120), seasonal(120), innovation
   real(dp) :: large_stationary(2101)
   type(series_profile_t) :: profile
   type(automatic_model_options_t) :: options
   type(automatic_model_result_t) :: result
   integer :: i

   autoregressive(1) = 0.5_dp
   do i = 2, size(autoregressive)
      innovation = 0.08_dp*sin(0.73_dp*real(i, dp))
      autoregressive(i) = 0.85_dp*autoregressive(i - 1) + innovation
   end do
   profile = profile_series(autoregressive, max_lag=12)
   if (profile%info /= 0) error stop "automatic profile status"
   if (.not. profile%autocorrelation_detected) error stop "automatic AR detection"
   if (size(profile%acf) /= 13) error stop "automatic ACF size"
   if (size(profile%squared_acf) /= 13) error stop "automatic squared ACF size"

   options%horizon = 6
   options%max_lag = 12
   options%max_ar_order = 6
   options%validation_size = 18
   result = automatic_model(autoregressive, options)
   if (result%info /= 0) error stop "automatic model status"
   if (size(result%candidates) < 6) error stop "automatic candidate count"
   if (result%selected /= 1) error stop "automatic ranking"
   if (.not. result%candidates(1)%selected) error stop "automatic selection flag"
   if (size(result%forecast) /= options%horizon) error stop "automatic horizon"
   if (any(result%forecast /= result%forecast)) error stop "automatic finite forecast"
   do i = 2, size(result%candidates)
      if (result%candidates(i)%rmse < result%candidates(i - 1)%rmse) then
         error stop "automatic candidate order"
      end if
   end do

   moderate(1) = 0.25_dp
   do i = 2, size(moderate)
      innovation = 0.08_dp*sin(0.73_dp*real(i, dp)) + &
         0.05_dp*sin(2.11_dp*real(i, dp)) + &
         0.03_dp*cos(1.37_dp*real(i, dp))
      moderate(i) = 0.55_dp*moderate(i - 1) + innovation
   end do
   profile = profile_series(moderate, max_lag=12)
   if (.not. profile%autocorrelation_detected .or. &
      profile%differencing_suggested) then
      error stop "automatic stationary ARMA gating profile"
   end if
   result = automatic_model(moderate, options)
   if (result%info /= 0) error stop "automatic ARMA status"
   do i = 1, size(result%candidates)
      if (index(result%candidates(i)%name, "ARMA") > 0) exit
   end do
   if (i > size(result%candidates)) error stop "automatic ARMA candidate"
   if (result%candidates(i)%ma_order < 1 .or. &
      .not. allocated(result%candidates(i)%ma_coefficients)) then
      error stop "automatic ARMA parameters"
   end if
   if (.not. allocated(result%candidates(i)%order_fits)) then
      error stop "automatic ARMA order diagnostics allocation"
   end if
   if (size(result%candidates(i)%order_fits) /= 9) then
      error stop "automatic ARMA order diagnostics size"
   end if
   if (count(result%candidates(i)%order_fits%converged) < 1) then
      error stop "automatic ARMA order diagnostics convergence"
   end if
   do i = 1, size(result%candidates)
      if (index(result%candidates(i)%name, "Yule-Walker") > 0) exit
   end do
   if (i > size(result%candidates)) error stop "automatic AR candidate"
   if (.not. allocated(result%candidates(i)%order_fits)) then
      error stop "automatic AR order diagnostics allocation"
   end if
   if (size(result%candidates(i)%order_fits) /= options%max_ar_order) then
      error stop "automatic AR order diagnostics size"
   end if

   do i = 1, size(seasonal)
      seasonal(i) = 2.0_dp + 0.1_dp*real(i, dp) + &
         1.5_dp*sin(2.0_dp*acos(-1.0_dp)*real(i, dp)/12.0_dp)
   end do
   options%frequency = 12
   options%validation_size = 12
   result = automatic_model(seasonal, options)
   if (result%info /= 0) error stop "automatic seasonal status"
   if (.not. result%profile%seasonality_detected) then
      error stop "automatic seasonal detection"
   end if
   if (.not. any([(index(result%candidates(i)%name, "Seasonal") > 0 .or. &
      index(result%candidates(i)%name, "Holt-Winters") > 0, &
      i=1, size(result%candidates))])) then
      error stop "automatic seasonal candidates"
   end if

   options%frequency = 1
   options%selection = "aicc"
   result = automatic_model(autoregressive, options)
   if (result%info /= 0) error stop "automatic AICc status"
   if (result%validation_size /= 0) error stop "automatic AICc full data"
   if (trim(result%selection_criterion) /= "aicc") then
      error stop "automatic AICc criterion"
   end if
   if (result%candidates(1)%effective_observations <= 0) then
      error stop "automatic AICc observations"
   end if
   do i = 2, size(result%candidates)
      if (result%candidates(i)%aicc < result%candidates(i - 1)%aicc) then
         error stop "automatic AICc candidate order"
      end if
   end do

   options%selection = "BIC"
   result = automatic_model(autoregressive, options)
   if (result%info /= 0) error stop "automatic BIC status"
   if (trim(result%selection_criterion) /= "bic") then
      error stop "automatic BIC criterion"
   end if
   do i = 2, size(result%candidates)
      if (result%candidates(i)%bic < result%candidates(i - 1)%bic) then
         error stop "automatic BIC candidate order"
      end if
   end do

   large_stationary(1) = 0.25_dp
   do i = 2, size(large_stationary)
      innovation = 0.12_dp*(2.0_dp*modulo( &
         sin(12.9898_dp*real(i, dp))*43758.5453_dp, 1.0_dp) - 1.0_dp)
      large_stationary(i) = 0.55_dp*large_stationary(i - 1) + innovation
   end do
   options%selection = "aicc"
   result = automatic_model(large_stationary, options)
   if (result%info /= 0) error stop "automatic large ARMA status"
   do i = 1, size(result%candidates)
      if (index(result%candidates(i)%name, "ARMA") > 0) exit
   end do
   if (i > size(result%candidates) .or. &
      .not. result%candidates(i)%converged) then
      error stop "automatic large ARMA screening"
   end if

   print '(a)', "automatic modeling tests passed"
end program test_automatic_modeling
