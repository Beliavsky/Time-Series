! SPDX-License-Identifier: GPL-3.0-only
! SPDX-FileComment: Regression tests for the esemifar translation.
program test_esemifar
   use kind_mod, only: dp
   use esemifar_mod
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   type(esemifar_smooth_t) :: smooth, derivative
   type(esemifar_order_selection_t) :: selection
   type(esemifar_model_t) :: model, derivative_model
   type(esemifar_forecast_t) :: normal_forecast, first_bootstrap, second_bootstrap
   type(esemifar_forecast_t) :: first_advanced, second_advanced
   real(dp), allocatable :: coefficients(:), series(:), trend(:)
   integer :: i

   coefficients = esemifar_arma_to_ma([0.75_dp], [0.5_dp], 4)
   call check(maxval(abs(coefficients - [1.0_dp, 1.25_dp, 0.9375_dp, &
      0.703125_dp, 0.52734375_dp])) < 1.0e-14_dp, 'ARMA infinite MA filter')
   coefficients = esemifar_arma_to_ar([0.75_dp], [0.5_dp], 4)
   call check(maxval(abs(coefficients - [-1.0_dp, 1.25_dp, -0.625_dp, &
      0.3125_dp, -0.15625_dp])) < 1.0e-14_dp, 'ARMA infinite AR filter')
   coefficients = esemifar_d_coefficients(0.3_dp, 6)
   call check(maxval(abs(coefficients - [1.0_dp, -0.3_dp, -0.105_dp, &
      -0.0595_dp, -0.0401625_dp, -0.02972025_dp, -0.0232808625_dp])) < &
      2.0e-15_dp, 'fractional differencing coefficients against R')
   call check(abs(esemifar_kdf(0, 0, 0.2_dp) - 9.425056505520672_dp) < &
      2.0e-13_dp .and. abs(esemifar_kdf(2, 2, 0.2_dp) - &
      1.279945945194150_dp) < 2.0e-12_dp, 'kernel constants against R')

   allocate(series(80), trend(80))
   do i = 1, 80
      trend(i) = 1.0_dp + 0.02_dp*real(i, dp) + &
         0.0005_dp*real(i, dp)**2
      series(i) = trend(i) + 0.12_dp*sin(1.31_dp*real(i, dp))
   end do
   smooth = esemifar_smooth(series, derivative_order=0, polynomial_order=1, &
      kernel_smoothness=1, bandwidth=0.15_dp)
   derivative = esemifar_derivative_fit(trend, 1, 0.15_dp, &
      polynomial_order=2, kernel_smoothness=1)
   call check(smooth%info == 0 .and. size(smooth%estimate) == 80 .and. &
      size(smooth%weights, 1) == 80 .and. size(smooth%residuals) == 80, &
      'boundary-aware local polynomial trend')
   call check(derivative%info == 0 .and. maxval(abs(derivative%estimate - &
      [(80.0_dp*(0.02_dp + 0.001_dp*real(i, dp)), i=1,80)])) < 2.0e-3_dp, &
      'local polynomial derivative reproduces a quadratic')

   selection = esemifar_order_selection(smooth%residuals, 0, 0, &
      truncation=40, max_iterations=100)
   call check(selection%info == 0 .and. selection%ar_order == 0 .and. &
      selection%ma_order == 0 .and. all(shape(selection%criterion) == [1, 1]), &
      'FARIMA information-criterion grid')
   model = esemifar_trend_fit(series, polynomial_order=1, kernel_smoothness=1, &
      initial_bandwidth=0.15_dp, p_max=0, q_max=0, max_iterations=6)
   call check(model%info == 0 .and. model%smoother%bandwidth > 0.0_dp .and. &
      model%smoother%bandwidth < 0.5_dp .and. model%iterations >= 1 .and. &
      size(model%bandwidth_history) == model%iterations .and. &
      all(ieee_is_finite(model%smoother%estimate)), &
      'iterative plug-in trend bandwidth selection')
   derivative_model = esemifar_derivative_ipi(series, 1, pilot_order=1, &
      kernel_smoothness=1, pilot_kernel_smoothness=1, initial_bandwidth=0.15_dp, &
      p_max=0, q_max=0, max_iterations=5)
   call check(derivative_model%info == 0 .and. &
      derivative_model%smoother%derivative_order == 1 .and. &
      derivative_model%smoother%bandwidth > 0.0_dp .and. &
      derivative_model%smoother%bandwidth < 0.5_dp .and. &
      all(ieee_is_finite(derivative_model%smoother%estimate)), &
      'iterative plug-in derivative bandwidth selection')
   smooth = model%smoother
   selection%fit = model%farima
   model%smoother = smooth
   model%farima = selection%fit
   model%variance_factor = 1.0_dp
   normal_forecast = esemifar_forecast_normal(model, 5, [0.8_dp, 0.95_dp])
   call check(normal_forecast%info == 0 .and. size(normal_forecast%mean) == 5 .and. &
      all(shape(normal_forecast%lower) == [5, 2]) .and. &
      all(normal_forecast%lower <= spread(normal_forecast%mean, 2, 2)) .and. &
      all(normal_forecast%upper >= spread(normal_forecast%mean, 2, 2)), &
      'analytic ESEMIFAR forecasts')
   first_bootstrap = esemifar_forecast_bootstrap(model, 4, [0.8_dp], 100, seed=81)
   second_bootstrap = esemifar_forecast_bootstrap(model, 4, [0.8_dp], 100, seed=81)
   call check(first_bootstrap%info == 0 .and. &
      maxval(abs(first_bootstrap%mean - second_bootstrap%mean)) < 1.0e-14_dp .and. &
      maxval(abs(first_bootstrap%lower - second_bootstrap%lower)) < 1.0e-14_dp, &
      'seeded residual-bootstrap forecasts')
   first_advanced = esemifar_forecast_bootstrap_advanced(model, 2, [0.8_dp], 2, &
      seed=19, burn_in=20, max_iterations=40)
   second_advanced = esemifar_forecast_bootstrap_advanced(model, 2, [0.8_dp], 2, &
      seed=19, burn_in=20, max_iterations=40)
   call check(first_advanced%info == 0 .and. second_advanced%info == 0 .and. &
      maxval(abs(first_advanced%mean - second_advanced%mean)) < 1.0e-14_dp .and. &
      maxval(abs(first_advanced%lower - second_advanced%lower)) < 1.0e-14_dp, &
      'seeded FARIMA-refitted predictive-root forecasts')

contains

   subroutine check(condition, label)
      !! Stop the test program when an assertion fails.
      logical, intent(in) :: condition !! Flag controlling condition.
      character(*), intent(in) :: label !! Label.

      if (.not. condition) then
         write (*, '(a)') 'FAILED: '//label
         error stop 1
      end if
   end subroutine check

end program test_esemifar
