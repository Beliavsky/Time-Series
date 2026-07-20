! SPDX-License-Identifier: BSD-2-Clause
! SPDX-FileComment: Regression tests for the itsmr translation.
program test_itsmr
   use kind_mod, only: dp
   use itsmr_mod
   implicit none
   type(itsmr_arma_model_t) :: model, fitted, hannan_fit, mle_fit, selected_fit, forecast_model, burg_model
   type(itsmr_arma_model_t) :: transform_model
   type(itsmr_innovations_t) :: innovations
   type(itsmr_forecast_t) :: forecast_result
   type(itsmr_arar_t) :: arar_result
   type(itsmr_residuals_t) :: residual_result
   type(harmonic_regression_result_t) :: harmonic_fit
   type(itsmr_transformed_forecast_t) :: transformed_forecast
   type(itsmr_rank_filter_t) :: rank_filter
   type(itsmr_randomness_tests_t) :: randomness
   type(itsmr_transform_t) :: transform_steps(2), seasonal_step(1)
   real(dp), allocatable :: covariance(:)
   real(dp), allocatable :: smoothed(:)
   real(dp) :: series(40), noise(40)
   real(dp) :: harmonic_series(24), expected_harmonic(6), angle
   real(dp) :: transform_series(40), seasonal_series(40)
   real(dp) :: spectral_series(16), expected_spectral(16)
   real(dp), parameter :: seasonal_pattern(4) = [1.0_dp, -1.0_dp, 2.0_dp, -2.0_dp]
   real(dp), parameter :: deaths(72) = [9007.0_dp, 8106.0_dp, 8928.0_dp, 9137.0_dp, &
      10017.0_dp, 10826.0_dp, 11317.0_dp, 10744.0_dp, 9713.0_dp, 9938.0_dp, &
      9161.0_dp, 8927.0_dp, 7750.0_dp, 6981.0_dp, 8038.0_dp, 8422.0_dp, &
      8714.0_dp, 9512.0_dp, 10120.0_dp, 9823.0_dp, 8743.0_dp, 9129.0_dp, &
      8710.0_dp, 8680.0_dp, 8162.0_dp, 7306.0_dp, 8124.0_dp, 7870.0_dp, &
      9387.0_dp, 9556.0_dp, 10093.0_dp, 9620.0_dp, 8285.0_dp, 8433.0_dp, &
      8160.0_dp, 8034.0_dp, 7717.0_dp, 7461.0_dp, 7776.0_dp, 7925.0_dp, &
      8634.0_dp, 8945.0_dp, 10078.0_dp, 9179.0_dp, 8037.0_dp, 8488.0_dp, &
      7874.0_dp, 8647.0_dp, 7792.0_dp, 6957.0_dp, 7726.0_dp, 8106.0_dp, &
      8890.0_dp, 9299.0_dp, 10625.0_dp, 9302.0_dp, 8314.0_dp, 8850.0_dp, &
      8265.0_dp, 8796.0_dp, 7836.0_dp, 6892.0_dp, 7791.0_dp, 8129.0_dp, &
      9115.0_dp, 9434.0_dp, 10484.0_dp, 9827.0_dp, 9110.0_dp, 9070.0_dp, &
      8633.0_dp, 9240.0_dp]
   integer :: i

   allocate(model%ar(1), model%ma(0))
   model%ar = 0.5_dp
   model%innovation_variance = 2.0_dp
   covariance = arma_acvf(model, 4)
   call check(maxval(abs(covariance - [2.0_dp/0.75_dp, 1.0_dp/0.75_dp, &
      0.5_dp/0.75_dp, 0.25_dp/0.75_dp, 0.125_dp/0.75_dp])) < 1.0e-11_dp, &
      'AR(1) theoretical ACVF')

   deallocate(model%ar, model%ma)
   allocate(model%ar(0), model%ma(1))
   model%ma = 0.5_dp
   model%innovation_variance = 2.0_dp
   covariance = arma_acvf(model, 3)
   call check(maxval(abs(covariance - [2.5_dp, 1.0_dp, 0.0_dp, 0.0_dp])) < 1.0e-11_dp, &
      'MA(1) theoretical ACVF')
   innovations = innovations_algorithm(covariance)
   call check(innovations%info == 0 .and. all(innovations%variance > 0.0_dp), &
      'innovations recursion')

   do i = 1, 40
      noise(i) = sin(1.7_dp*real(i, dp)) + 0.3_dp*cos(0.43_dp*real(i, dp))
   end do
   series(1) = noise(1)
   do i = 2, 40
      series(i) = noise(i) + 0.4_dp*noise(i - 1)
   end do
   fitted = innovations_ma_fit(series, 1, 12)
   call check(fitted%info == 0 .and. size(fitted%ma) == 1 .and. &
      fitted%innovation_variance > 0.0_dp, 'innovations MA fit')
   call check(fitted%ma_standard_error(1) > 0.0_dp .and. &
      fitted%aicc < huge(1.0_dp), 'innovations MA inference')

   series = 0.0_dp
   do i = 2, 40
      series(i) = 0.55_dp*series(i - 1) + noise(i) + 0.25_dp*noise(i - 1)
   end do
   hannan_fit = hannan_rissanen_fit(series, 1, 1)
   call check(hannan_fit%info == 0 .and. size(hannan_fit%ar) == 1 .and. &
      size(hannan_fit%ma) == 1, 'Hannan-Rissanen ARMA fit')
   call check(hannan_fit%innovation_variance > 0.0_dp .and. &
      all(hannan_fit%ar_standard_error > 0.0_dp) .and. &
      all(hannan_fit%ma_standard_error > 0.0_dp), 'Hannan-Rissanen inference')

   mle_fit = arma_mle_fit(series, 1, 1, 80, 1.0e-5_dp)
   call check(mle_fit%info == 0 .and. mle_fit%converged, 'exact ARMA maximum likelihood')
   call check(abs(mle_fit%ar(1)) < 1.0_dp .and. abs(mle_fit%ma(1)) < 1.0_dp .and. &
      all(mle_fit%ar_standard_error > 0.0_dp) .and. &
      all(mle_fit%ma_standard_error > 0.0_dp), 'stable ARMA inference')
   selected_fit = arma_autofit(series, [0, 1], [0], 80, 1.0e-5_dp)
   call check(selected_fit%info == 0 .and. selected_fit%aicc < huge(1.0_dp), &
      'automatic ARMA order selection')

   allocate(forecast_model%ar(1), forecast_model%ma(1))
   forecast_model%ar = 0.5_dp
   forecast_model%ma = 0.2_dp
   forecast_model%innovation_variance = 2.0_dp
   covariance = arma_infinite_ma(forecast_model, 3)
   call check(maxval(abs(covariance - [1.0_dp, 0.7_dp, 0.35_dp, 0.175_dp])) < 1.0e-12_dp, &
      'ARMA infinite-MA coefficients')
   covariance = arma_infinite_ar(forecast_model, 3)
   call check(maxval(abs(covariance - [1.0_dp, -0.7_dp, 0.14_dp, -0.028_dp])) < 1.0e-12_dp, &
      'ARMA infinite-AR coefficients')
   residual_result = arma_residuals(series, forecast_model)
   call check(residual_result%info == 0 .and. &
      maxval(abs(residual_result%fitted + residual_result%innovations - series)) < 1.0e-12_dp, &
      'exact ARMA residual decomposition')
   call check(all(residual_result%variance > 0.0_dp) .and. &
      maxval(abs(residual_result%standardized - residual_result%innovations/ &
      sqrt(residual_result%variance))) < 1.0e-12_dp, 'standardized ARMA innovations')
   deallocate(forecast_model%ma)
   allocate(forecast_model%ma(0))
   forecast_result = arma_forecast(series, forecast_model, 4, 0.90_dp)
   call check(forecast_result%info == 0 .and. &
      abs(forecast_result%mean(1) - (sum(series)/40.0_dp + &
      0.5_dp*(series(40) - sum(series)/40.0_dp))) < 1.0e-11_dp, 'AR forecast recursion')
   call check(maxval(abs(forecast_result%standard_error - &
      sqrt(2.0_dp*[1.0_dp, 1.25_dp, 1.3125_dp, 1.328125_dp]))) < 1.0e-11_dp .and. &
      all(forecast_result%lower < forecast_result%mean) .and. &
      all(forecast_result%upper > forecast_result%mean), 'AR forecast uncertainty')

   arar_result = arar_forecast(deaths, 10)
   call check(arar_result%info == 0 .and. all(arar_result%lags(2:) > arar_result%lags(:3)), &
      'ARAR sparse lag selection')
   call check(maxval(abs(arar_result%mean/[8167.8_dp, 7195.8_dp, 7982.0_dp, 8283.5_dp, &
      9144.1_dp, 9464.9_dp, 10541.0_dp, 9640.8_dp, 8902.7_dp, 9096.7_dp] - 1.0_dp)) < 1.0e-4_dp, &
      'ARAR reference forecasts')
   call check(maxval(abs((arar_result%standard_error/arar_result%standard_error(1))/ &
      ([323.35_dp, 375.68_dp, 392.34_dp, 414.78_dp, 431.69_dp, 441.90_dp, &
      449.82_dp, 455.98_dp, 460.45_dp, 463.77_dp]/323.35_dp) - 1.0_dp)) < 1.0e-4_dp, &
      'ARAR reference standard errors')

   burg_model = burg_ar_fit(deaths, 4)
   call check(burg_model%info == 0 .and. maxval(abs(burg_model%ar - &
      [0.7539377705729109_dp, 0.0328283459223125_dp, -0.1383061772755615_dp, &
      -0.1579254556424599_dp])) < 1.0e-12_dp, 'Burg AR reference coefficients')
   call check(maxval(abs(burg_model%ar_standard_error - &
      [0.1160858748708147_dp, 0.1443546690010265_dp, 0.1443546690010264_dp, &
      0.1160858748708147_dp])) < 1.0e-12_dp .and. burg_model%innovation_variance > 0.0_dp, &
      'Burg AR inference')

   do i = 1, 24
      angle = 2.0_dp*acos(-1.0_dp)*real(i - 1, dp)/4.0_dp
      harmonic_series(i) = 3.0_dp + 0.2_dp*real(i, dp) + 2.0_dp*cos(angle) - 0.5_dp*sin(angle)
   end do
   harmonic_fit = harmonic_regression(harmonic_series, [4.0_dp], 6, 1)
   call check(harmonic_fit%info == 0 .and. maxval(abs(harmonic_fit%coefficients - &
      [3.0_dp, 0.2_dp, 2.0_dp, -0.5_dp])) < 1.0e-12_dp .and. &
      maxval(abs(harmonic_fit%residuals)) < 1.0e-12_dp, 'harmonic regression decomposition')
   do i = 1, 6
      angle = 2.0_dp*acos(-1.0_dp)*real(23 + i, dp)/4.0_dp
      expected_harmonic(i) = 3.0_dp + 0.2_dp*real(24 + i, dp) + &
         2.0_dp*cos(angle) - 0.5_dp*sin(angle)
   end do
   call check(maxval(abs(harmonic_fit%forecast - expected_harmonic)) < 1.0e-12_dp .and. &
      harmonic_fit%polynomial_order == 1, 'harmonic regression extrapolation')

   allocate(transform_model%ar(0), transform_model%ma(0))
   transform_model%innovation_variance = 2.0_dp
   transform_series = [(real(i, dp), i=1, 40)]
   transform_steps(1)%kind = itsmr_transform_difference
   transform_steps(1)%lag = 1
   transformed_forecast = transformed_arma_forecast(transform_series, transform_model, &
      transform_steps(:1), 4, 0.90_dp)
   call check(transformed_forecast%info == 0 .and. &
      maxval(abs(transformed_forecast%mean - [41.0_dp, 42.0_dp, 43.0_dp, 44.0_dp])) < 1.0e-12_dp, &
      'inverse-difference forecasts')
   call check(maxval(abs(transformed_forecast%standard_error - &
      sqrt(2.0_dp*[1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]))) < 1.0e-12_dp, &
      'integrated forecast uncertainty')

   transform_series = exp(0.1_dp*[(real(i, dp), i=1, 40)])
   transform_steps(1)%kind = itsmr_transform_log
   transform_steps(2)%kind = itsmr_transform_difference
   transform_steps(2)%lag = 1
   transformed_forecast = transformed_arma_forecast(transform_series, transform_model, transform_steps, 4)
   call check(transformed_forecast%info == 0 .and. maxval(abs(transformed_forecast%mean - &
      exp(0.1_dp*[41.0_dp, 42.0_dp, 43.0_dp, 44.0_dp]))) < 1.0e-10_dp .and. &
      all(transformed_forecast%lower < transformed_forecast%mean) .and. &
      all(transformed_forecast%upper > transformed_forecast%mean), 'log-difference forecast reversal')

   do i = 1, 40
      seasonal_series(i) = 10.0_dp + seasonal_pattern(modulo(i - 1, 4) + 1)
   end do
   seasonal_step(1)%kind = itsmr_transform_seasonal
   seasonal_step(1)%lag = 4
   transformed_forecast = transformed_arma_forecast(seasonal_series, transform_model, seasonal_step, 4)
   call check(transformed_forecast%info == 0 .and. &
      maxval(abs(transformed_forecast%mean - [11.0_dp, 9.0_dp, 12.0_dp, 8.0_dp])) < 1.0e-12_dp, &
      'seasonal forecast reversal')

   do i = 1, 16
      spectral_series(i) = 3.0_dp + 2.0_dp*cos(2.0_dp*acos(-1.0_dp)*2.0_dp* &
         real(i - 1, dp)/16.0_dp) + cos(2.0_dp*acos(-1.0_dp)*5.0_dp*real(i - 1, dp)/16.0_dp)
      expected_spectral(i) = 3.0_dp + 2.0_dp*cos(2.0_dp*acos(-1.0_dp)*2.0_dp* &
         real(i - 1, dp)/16.0_dp)
   end do
   rank_filter = spectral_rank_filter(spectral_series, 1)
   call check(rank_filter%info == 0 .and. all(rank_filter%retained_bins == [2]) .and. &
      maxval(abs(rank_filter%filtered - expected_spectral)) < 1.0e-12_dp, &
      'spectral rank reconstruction')
   call check(abs(rank_filter%frequency(1) - 0.125_dp) < 1.0e-12_dp .and. &
      abs(rank_filter%amplitude(1) - 2.0_dp) < 1.0e-12_dp, 'spectral rank diagnostics')
   rank_filter = spectral_rank_filter(spectral_series, 0)
   call check(maxval(abs(rank_filter%filtered - 3.0_dp)) < 1.0e-12_dp, &
      'spectral rank mean preservation')

   randomness = residual_randomness_tests(noise, 10)
   call check(randomness%info == 0 .and. abs(randomness%ljung_box - &
      138.08060630910168_dp) < 1.0e-11_dp .and. abs(randomness%mcleod_li - &
      57.83028377675973_dp) < 1.0e-11_dp, 'portmanteau randomness tests')
   call check(abs(randomness%ljung_box_p_value/1.0420756605933216e-24_dp - 1.0_dp) < 1.0e-10_dp .and. &
      abs(randomness%mcleod_li_p_value/9.304696888919567e-9_dp - 1.0_dp) < 1.0e-10_dp, &
      'portmanteau p-values')
   call check(randomness%turning_points == 20.0_dp .and. randomness%difference_signs == 18.0_dp .and. &
      randomness%rank_statistic == 365.0_dp .and. &
      abs(randomness%turning_point_p_value - 0.040666630545760925_dp) < 1.0e-12_dp .and. &
      abs(randomness%difference_sign_p_value - 0.41707705952056467_dp) < 1.0e-12_dp .and. &
      abs(randomness%rank_p_value - 0.5601951093846056_dp) < 1.0e-12_dp, &
      'order-based randomness tests')

   smoothed = itsmr_moving_average([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp], 1)
   call check(maxval(abs(smoothed - [4.0_dp/3.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, &
      14.0_dp/3.0_dp])) < 1.0e-12_dp, 'endpoint-replicated moving average')
   smoothed = itsmr_exponential_smooth([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp], 0.5_dp)
   call check(maxval(abs(smoothed - [1.0_dp, 1.5_dp, 2.25_dp, 3.125_dp, 4.0625_dp])) < 1.0e-12_dp, &
      'ITSMR exponential smoothing')
   smoothed = itsmr_fft_smooth(spectral_series, 0.25_dp)
   call check(maxval(abs(smoothed - expected_spectral)) < 1.0e-12_dp, 'ITSMR low-pass FFT smoothing')
   smoothed = itsmr_seasonal_component(seasonal_series, 4)
   call check(maxval(abs(smoothed - &
      [(seasonal_pattern(modulo(i - 1, 4) + 1), i=1, 40)])) < 1.0e-12_dp, &
      'ITSMR seasonal component')
   print '(a)', 'All itsmr_mod tests passed.'

contains

   subroutine check(ok, name)
      !! Stop the test program when a named assertion fails.
      logical, intent(in) :: ok !! Flag controlling ok.
      character(len=*), intent(in) :: name !! Name.
      if (.not. ok) then
         print '(a)', 'FAILED: '//name
         error stop 1
      end if
   end subroutine check
end program test_itsmr
