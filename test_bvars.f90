! SPDX-License-Identifier: GPL-3.0-or-later
! SPDX-FileComment: Numerical tests for algorithms translated from R bvars.
program test_bvars
   use kind_mod, only: dp
   use bvars_mod
   use bvartools_mod, only: bvartools_fevd_t
   use random_mod, only: set_random_seed
   use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
   implicit none

   real(dp) :: series(12, 2), exogenous(12, 1)
   real(dp) :: zero_normal(2, 2), fitted_normal(2, 10, 3)
   real(dp) :: forecast_normal(2, 2, 3), future_exogenous(2, 1)
   real(dp) :: conditional_path(2, 2)
   real(dp), allocatable :: matrix_draw(:, :), shocks(:, :, :)
   type(bvars_data_t) :: data
   type(bvars_prior_t) :: prior
   type(bvars_draws_t) :: draws, student_draws, centered_sv_draws
   type(bvars_draws_t) :: noncentered_sv_draws
   type(bvars_draws_t) :: centered_student_sv_draws
   type(bvars_draws_t) :: noncentered_student_sv_draws
   type(bvars_fitted_t) :: fitted
   type(bvars_predictive_t) :: prediction
   type(bvartools_fevd_t) :: fevd
   integer :: time
   real(dp), allocatable :: future_scale(:, :)
   real(dp) :: sv_mixture(3, 1)
   real(dp), allocatable :: generated_sv_mixture(:, :)
   real(dp), allocatable :: future_variance(:, :)

   do time = 1, 12
      series(time, 1) = 0.1_dp*real(time, dp) + &
         0.05_dp*sin(real(time, dp))
      series(time, 2) = 0.2_dp + 0.7_dp*series(time, 1) + &
         0.03_dp*cos(real(time, dp))
      exogenous(time, 1) = sin(0.2_dp*real(time, dp))
   end do
   data = bvars_prepare(series, 2, exogenous)
   call assert_true(data%info == 0 .and. data%variables == 2 .and. &
      data%observations == 10 .and. data%regressors == 6, &
      'prepared BVAR dimensions')

   prior = bvars_default_prior(data, stationary=[.false., .true.], &
      innovation_variance=[1.0_dp, 2.0_dp])
   call assert_true(prior%info == 0 .and. &
      all(shape(prior%coefficient_mean) == [2, 6]), &
      'default prior dimensions')
   call assert_close(prior%coefficient_mean(1, 1), 1.0_dp, 1.0e-14_dp, &
      'random-walk prior mean')
   call assert_close(prior%coefficient_mean(2, 2), 0.0_dp, 1.0e-14_dp, &
      'stationary prior mean')

   zero_normal = 0.0_dp
   matrix_draw = bvars_matrix_normal_from_standard(reshape( &
      [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], [2, 2]), &
      reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2]), &
      reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2]), zero_normal)
   call assert_close(matrix_draw(2, 2), 4.0_dp, 1.0e-14_dp, &
      'matrix-normal supplied draw')

   call set_random_seed(810)
   draws = bvars_conjugate_draws(data, prior, 3, burnin=1)
   call assert_true(draws%info == 0 .and. draws%retained_draws == 3 .and. &
      all(shape(draws%coefficient) == [2, 6, 3]) .and. &
      all(shape(draws%covariance) == [2, 2, 3]), &
      'conjugate posterior draws')
   call assert_true(all(draws%covariance(1, 1, :) > 0.0_dp), &
      'positive covariance draws')

   call assert_true(bvars_student_log_kernel(8.0_dp, [1.0_dp, 1.5_dp]) > &
      -huge(1.0_dp), 'finite Student-t degrees kernel')
   call set_random_seed(811)
   student_draws = bvars_student_t_draws(data, prior, 3, burnin=1, &
      initial_df=10.0_dp)
   call assert_true(student_draws%info == 0 .and. &
      student_draws%retained_draws == 3 .and. &
      all(student_draws%degrees_of_freedom > 2.0_dp) .and. &
      all(student_draws%scale_mixture > 0.0_dp), &
      'Student-t posterior sampler')
   future_scale = bvars_student_scale_forecast( &
      student_draws%degrees_of_freedom, 2)
   call assert_true(all(shape(future_scale) == [2, 3]) .and. &
      all(future_scale > 0.0_dp), 'Student-t forecast scales')

   sv_mixture(:, 1) = [0.0_dp, -1.2704_dp, 4.9348_dp]
   call set_random_seed(809)
   generated_sv_mixture = bvars_sv_auxiliary_mixture(2, sample_count=2000, &
      components=3, max_iterations=100)
   call assert_true(all(shape(generated_sv_mixture) == [3, 3]) .and. &
      abs(sum(generated_sv_mixture(1, :)) - 1.0_dp) < 1.0e-12_dp .and. &
      all(generated_sv_mixture(3, :) > 0.0_dp), &
      'generated stochastic-volatility auxiliary mixture')
   call set_random_seed(812)
   centered_sv_draws = bvars_centered_sv_draws(data, prior, sv_mixture, 3, &
      burnin=1)
   call assert_true(centered_sv_draws%info == 0 .and. &
      centered_sv_draws%retained_draws == 3 .and. &
      all(shape(centered_sv_draws%log_variance) == [10, 3]) .and. &
      all(centered_sv_draws%common_variance > 0.0_dp) .and. &
      all(abs(centered_sv_draws%persistence) < 1.0_dp), &
      'centred common stochastic-volatility sampler')
   future_variance = bvars_common_variance_forecast( &
      centered_sv_draws%log_variance(10, :), &
      centered_sv_draws%persistence, &
      sqrt(centered_sv_draws%volatility_innovation_variance), 2)
   call assert_true(all(shape(future_variance) == [2, 3]) .and. &
      all(future_variance > 0.0_dp), 'common stochastic-volatility forecasts')

   call set_random_seed(813)
   noncentered_sv_draws = bvars_noncentered_sv_draws(data, prior, &
      sv_mixture, 3, burnin=1)
   call assert_true(noncentered_sv_draws%info == 0 .and. &
      noncentered_sv_draws%retained_draws == 3 .and. &
      all(shape(noncentered_sv_draws%log_variance) == [10, 3]) .and. &
      all(noncentered_sv_draws%common_variance > 0.0_dp) .and. &
      all(abs(noncentered_sv_draws%persistence) < 1.0_dp) .and. &
      all(abs(noncentered_sv_draws%volatility_loading) > 0.0_dp), &
      'non-centred common stochastic-volatility sampler')

   call set_random_seed(814)
   centered_student_sv_draws = bvars_student_sv_draws(data, prior, &
      sv_mixture, .true., 3, burnin=1, initial_df=10.0_dp)
   call assert_true(centered_student_sv_draws%info == 0 .and. &
      centered_student_sv_draws%retained_draws == 3 .and. &
      all(centered_student_sv_draws%degrees_of_freedom > 2.0_dp) .and. &
      all(centered_student_sv_draws%scale_mixture > 0.0_dp) .and. &
      all(centered_student_sv_draws%common_variance > 0.0_dp), &
      'centred Student-t stochastic-volatility sampler')

   call set_random_seed(815)
   noncentered_student_sv_draws = bvars_student_sv_draws(data, prior, &
      sv_mixture, .false., 3, burnin=1, initial_df=10.0_dp)
   call assert_true(noncentered_student_sv_draws%info == 0 .and. &
      noncentered_student_sv_draws%retained_draws == 3 .and. &
      all(noncentered_student_sv_draws%degrees_of_freedom > 2.0_dp) .and. &
      all(noncentered_student_sv_draws%scale_mixture > 0.0_dp) .and. &
      all(noncentered_student_sv_draws%common_variance > 0.0_dp), &
      'non-centred Student-t stochastic-volatility sampler')

   shocks = bvars_shocks(draws, data)
   call assert_true(all(shape(shocks) == [2, 10, 3]), &
      'posterior shock dimensions')
   fitted_normal = 0.0_dp
   fitted = bvars_fitted_from_random(draws, data, fitted_normal)
   call assert_true(fitted%info == 0 .and. &
      all(shape(fitted%value) == [2, 10, 3]), &
      'posterior fitted-density draws')
   call assert_close(fitted%value(1, 1, 1), &
      dot_product(draws%coefficient(1, :, 1), data%x(:, 1)), &
      1.0e-12_dp, 'zero-innovation fitted value')

   forecast_normal = 0.0_dp
   future_exogenous(:, 1) = [0.1_dp, 0.2_dp]
   prediction = bvars_forecast_from_random(draws, data, 2, &
      forecast_normal, future_exogenous)
   call assert_true(prediction%info == 0 .and. &
      all(shape(prediction%path) == [2, 2, 3]) .and. &
      all(abs(prediction%path - prediction%mean) < 1.0e-12_dp), &
      'recursive posterior forecasts')

   conditional_path = ieee_value(0.0_dp, ieee_quiet_nan)
   conditional_path(1, 1) = 0.75_dp
   conditional_path(2, 2) = 0.60_dp
   prediction = bvars_conditional_forecast_from_random(draws, data, 2, &
      forecast_normal, conditional_path, future_exogenous)
   call assert_true(prediction%info == 0 .and. &
      all(abs(prediction%path(1, 1, :) - 0.75_dp) < 1.0e-14_dp) .and. &
      all(abs(prediction%path(2, 2, :) - 0.60_dp) < 1.0e-14_dp), &
      'conditional recursive forecasts')

   fevd = bvars_fevd(draws, 2, 2)
   call assert_true(fevd%info == 0 .and. &
      all(shape(fevd%paths) == [2, 2, 3, 3]), &
      'orthogonalized posterior FEVD')

   print '(a)', 'bvars tests passed'

contains

   subroutine assert_true(condition, message)
      !! Stop the test program when a logical assertion fails.
      logical, intent(in) :: condition !! Assertion condition.
      character(len=*), intent(in) :: message !! Failure description.

      if (.not. condition) error stop 'FAIL: '//trim(message)
   end subroutine assert_true

   subroutine assert_close(actual, expected, tolerance, message)
      !! Stop the test program when two real values differ materially.
      real(dp), intent(in) :: actual !! Computed value.
      real(dp), intent(in) :: expected !! Reference value.
      real(dp), intent(in) :: tolerance !! Absolute tolerance.
      character(len=*), intent(in) :: message !! Failure description.

      if (abs(actual - expected) > tolerance) then
         print '(a,2es24.15)', trim(message)//': ', actual, expected
         error stop 'FAIL'
      end if
   end subroutine assert_close

end program test_bvars
