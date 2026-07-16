! SPDX-License-Identifier: GPL-2.0-or-later
! SPDX-FileComment: Regression tests for the TSSS translation.
program test_tsss
   use kind_mod, only: dp
   use tsss_mod
   use time_series_random_mod, only: set_random_seed
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   type(tsss_particle_filter_t) :: filtered, nonlinear
   type(tsss_tvvar_t) :: changing_variance
   type(tsss_tvar_t) :: changing_ar
   type(tsss_tvspc_t) :: spectrum
   type(tsss_ngsmooth_t) :: grid_smoother
   type(tsss_simulation_t) :: gaussian_simulation, non_gaussian_simulation
   type(tsss_structural_model_t) :: structural
   type(tsss_tsmooth_t) :: general_smoother
   type(tsss_density_t) :: density
   type(tsss_kl_t) :: divergence
   type(tsss_boxcox_t) :: transformed
   real(dp) :: observations(12), expected, series(120), coefficients(1, 3)
   real(dp) :: transition(1, 1), loading(1, 1), observation_loading(1)
   real(dp) :: covariance(1, 1), initial_state(1), spacing, integral
   integer :: i

   observations = 0.0_dp
   call set_random_seed(12345)
   filtered = tsss_particle_filter(observations, 500, tsss_system_gaussian, 3, &
      tsss_initial_zero, 1.0_dp, 0.0_dp, 0.99_dp, 4.0_dp, 1.0_dp, -5.0_dp, 5.0_dp)
   expected = -0.5_dp*real(size(observations), dp)*log(2.0_dp*acos(-1.0_dp))
   call check(filtered%info == 0 .and. all(shape(filtered%quantile) == [12, 7]) .and. &
      maxval(abs(filtered%quantile)) == 0.0_dp .and. maxval(abs(filtered%mean)) == 0.0_dp .and. &
      abs(filtered%log_likelihood - expected) < 1.0e-11_dp, &
      'linear fixed-lag particle smoothing')

   call set_random_seed(54321)
   nonlinear = tsss_nonlinear_particle_filter(observations, 600, 2, 1.0_dp, &
      0.2_dp, -20.0_dp, 20.0_dp)
   call check(nonlinear%info == 0 .and. all(shape(nonlinear%quantile) == [12, 7]) .and. &
      all(nonlinear%quantile(:, 2:) >= nonlinear%quantile(:, :6)) .and. &
      ieee_is_finite(nonlinear%log_likelihood), 'nonlinear fixed-lag particle smoothing')

   do i = 1, size(series)
      series(i) = merge(1.0_dp, -1.0_dp, mod(i, 2) == 0)
   end do
   changing_variance = tsss_tvvar(series, 1)
   call check(changing_variance%info == 0 .and. size(changing_variance%variance) == 60 .and. &
      all(changing_variance%variance > 0.0_dp) .and. &
      maxval(changing_variance%variance) - minval(changing_variance%variance) < 1.0e-8_dp, &
      'time-varying variance')

   series(1) = 0.25_dp
   do i = 2, size(series)
      series(i) = 0.7_dp*series(i - 1) + 0.1_dp*sin(0.37_dp*real(i, dp))
   end do
   changing_ar = tsss_tvar(series, 1, 1, 10)
   call check(changing_ar%info == 0 .and. all(shape(changing_ar%ar) == [1, 12]) .and. &
      ieee_is_finite(changing_ar%log_likelihood), 'time-varying autoregression')

   coefficients(1, :) = 0.5_dp
   spectrum = tsss_tvspc(coefficients, 2.0_dp, 4, 10)
   call check(spectrum%info == 0 .and. all(shape(spectrum%log_spectrum) == [11, 3]) .and. &
      abs(spectrum%log_spectrum(0, 1) - log10(8.0_dp)) < 1.0e-12_dp .and. &
      abs(spectrum%frequency(10) - 0.5_dp) < 1.0e-12_dp, 'evolutionary spectrum')

   observations = [(0.2_dp*sin(0.4_dp*real(i, dp)), i = 1, 12)]
   grid_smoother = tsss_ngsmth(observations, tsss_noise_laplace, 0.05_dp, 1.0_dp, &
      tsss_noise_pearson, 0.2_dp, 2.5_dp, 80)
   spacing = grid_smoother%grid(1) - grid_smoother%grid(0)
   integral = spacing*(sum(grid_smoother%density(:, 6)) - &
      0.5_dp*(grid_smoother%density(0, 6) + grid_smoother%density(80, 6)))
   call check(grid_smoother%info == 0 .and. all(shape(grid_smoother%quantile) == [12, 7]) .and. &
      abs(integral - 1.0_dp) < 1.0e-10_dp .and. &
      all(grid_smoother%quantile(:, 2:) >= grid_smoother%quantile(:, :6)), &
      'non-Gaussian grid smoother')

   transition = 1.0_dp
   loading = 1.0_dp
   observation_loading = 1.0_dp
   covariance = 0.0_dp
   initial_state = 2.0_dp
   gaussian_simulation = tsss_simssm(transition, loading, observation_loading, covariance, &
      0.0_dp, initial_state, 8)
   call check(gaussian_simulation%info == 0 .and. all(gaussian_simulation%state == 2.0_dp) .and. &
      all(gaussian_simulation%observation == 2.0_dp), 'Gaussian state-space simulation')

   covariance = 0.1_dp
   call set_random_seed(13579)
   non_gaussian_simulation = tsss_ngsim(transition, loading, observation_loading, covariance, &
      0.2_dp, initial_state, 20, tsss_noise_pearson, 3.0_dp, tsss_noise_laplace, 1.0_dp)
   call check(non_gaussian_simulation%info == 0 .and. &
      all(shape(non_gaussian_simulation%state) == [1, 20]) .and. &
      all(ieee_is_finite(non_gaussian_simulation%observation)), &
      'non-Gaussian state-space simulation')

   structural = tsss_structural_model(trend=[1.0_dp, 0.5_dp], seasonal_order=1, &
      seasonal=[0.2_dp, -0.1_dp, 0.1_dp], ar_coefficients=[0.5_dp, -0.2_dp], &
      ar_initial=[0.3_dp, 0.1_dp], trend_variance=0.0_dp, seasonal_variance=0.0_dp, &
      ar_variance=0.0_dp, observation_variance=0.0_dp)
   call check(structural%info == 0 .and. all(shape(structural%transition) == [7, 7]) .and. &
      all(shape(structural%system_loading) == [7, 3]) .and. &
      all(structural%observation_loading == [1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
         0.0_dp, 1.0_dp, 0.0_dp]), 'structural state-space model builder')
   gaussian_simulation = tsss_simssm(structural%transition, structural%system_loading, &
      structural%observation_loading, structural%system_covariance, &
      structural%observation_variance, structural%initial_state, 1)
   call check(gaussian_simulation%info == 0 .and. &
      abs(gaussian_simulation%observation(1) - 1.43_dp) < 1.0e-12_dp, &
      'structural model simulation')

   observations = [(0.1_dp*real(i, dp), i = 1, 12)]
   transition = 1.0_dp
   loading = 1.0_dp
   observation_loading = 1.0_dp
   covariance = 0.02_dp
   initial_state = 0.0_dp
   general_smoother = tsss_tsmooth(observations, transition, loading, observation_loading, &
      covariance, 0.1_dp, initial_state, reshape([1.0_dp], [1, 1]), filter_end=10, &
      predict_end=15, missing_start=[5], missing_count=[2])
   call check(general_smoother%info == 0 .and. &
      all(shape(general_smoother%smoothed_state) == [1, 15]) .and. &
      all(shape(general_smoother%smoothed_covariance) == [1, 1, 15]) .and. &
      general_smoother%observations == 8 .and. &
      ieee_is_finite(general_smoother%log_likelihood) .and. &
      abs(general_smoother%prediction_error(5)) > 0.0_dp, &
      'general TSSS smoothing compatibility')

   density = tsss_pdfunc(tsss_density_gaussian, [0.0_dp, 1.0_dp, 0.0_dp], -1.0_dp, 1.0_dp, 3)
   call check(density%info == 0 .and. all(density%grid == [-1.0_dp, 0.0_dp, 1.0_dp]) .and. &
      abs(density%density(2) - 1.0_dp/sqrt(2.0_dp*acos(-1.0_dp))) < 1.0e-12_dp, &
      'TSSS probability density utility')

   divergence = tsss_klinfo(tsss_density_gaussian, [0.0_dp, 1.0_dp, 0.0_dp], &
      tsss_density_gaussian, [0.0_dp, 1.0_dp, 0.0_dp], -8.0_dp, 8.0_dp)
   call check(divergence%info == 0 .and. maxval(abs(divergence%information)) < 1.0e-14_dp .and. &
      abs(divergence%reference_mass(4) - 1.0_dp) < 1.0e-10_dp, &
      'Kullback-Leibler information')

   transformed = tsss_boxcox([1.0_dp, 2.0_dp, 4.0_dp, 8.0_dp, 16.0_dp])
   call check(transformed%info == 0 .and. abs(transformed%lambda(1) - 1.0_dp) < 1.0e-12_dp .and. &
      abs(transformed%lambda(21) + 1.0_dp) < 1.0e-12_dp .and. &
      minval(transformed%adjusted_aic) == &
      minval(transformed%adjusted_aic, mask=transformed%lambda == transformed%best_lambda), &
      'TSSS Box-Cox likelihood selection')
   print '(a)', 'All tsss_mod tests passed.'

contains

   subroutine check(ok, name)
      ! Stop the test program when a named assertion fails.
      logical, intent(in) :: ok
      character(len=*), intent(in) :: name

      if (.not. ok) then
         print '(a)', 'FAILED: '//name
         error stop 1
      end if
   end subroutine check
end program test_tsss
