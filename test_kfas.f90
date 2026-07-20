! SPDX-License-Identifier: GPL-2.0-or-later
! SPDX-FileComment: Regression tests for the KFAS translation.
program test_kfas
   use kind_mod, only: dp
   use kfas_mod
   use linalg_mod, only: inverse_logdet
   implicit none
   real(dp) :: y(5)
   type(ssm_model_t) :: model
   type(kfs_filter_t) :: filtered
   type(kfs_smoother_t) :: smoothed
   type(kfs_smoother_t) :: fast_smoothed
   type(kfs_conditional_t) :: conditional
   type(kfs_prediction_t) :: predicted
   type(kfs_filter_t) :: diffuse_filtered
   type(kfs_disturbance_t) :: disturbances
   type(kfs_loglik_result_t) :: loglik
   real(dp), allocatable :: standardized(:, :)
   real(dp) :: covariance(2, 2), covariance_inverse(2, 2), log_determinant
   real(dp) :: singular_covariance(2, 2, 2)
   real(dp) :: singular_lag_covariance(2, 2, 1)
   integer :: status

   covariance = reshape([1.0_dp, -2.0_dp, -2.0_dp, 5.0_dp], [2, 2])
   call inverse_logdet(covariance, covariance_inverse, log_determinant, &
      status, 100.0_dp*epsilon(1.0_dp))
   call check(status == 0 .and. abs(log_determinant) < 1.0e-12_dp .and. &
      maxval(abs(covariance_inverse - &
      reshape([5.0_dp, 2.0_dp, 2.0_dp, 1.0_dp], [2, 2]))) < 1.0e-12_dp, &
      'pivoted covariance inverse')

   y = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
   model = make_local_level(y, 1.0_dp, 0.5_dp, 0.0_dp, 2.0_dp)
   filtered = kfs_filter(model)
   loglik = kfs_loglik(model)
   call check(filtered%info == 0, 'filter status')
   call check(loglik%info == 0 .and. abs(loglik%value - filtered%log_likelihood) < 1.0e-12_dp, &
      'log likelihood result')
   call check(filtered%observations == 5, 'observation count')
   call check(abs(filtered%a_filt(1, 1) - 2.0_dp/3.0_dp) < 1.0e-12_dp, 'first update')
   smoothed = kfs_smooth(model, filtered)
   call check(smoothed%info == 0, 'smoother status')
   call check(smoothed%state(1, 1) > filtered%a_filt(1, 1), 'backward smoothing')
   call check(all(shape(smoothed%lag_one_covariance) == [1, 1, 4]) .and. &
      all(shape(smoothed%conditional_matrix) == [1, 1, 5]) .and. &
      all(smoothed%conditional_covariance >= 0.0_dp), &
      'smoother covariance outputs')
   fast_smoothed = kfs_fast_smooth(model, filtered)
   call check(fast_smoothed%info == 0 .and. &
      maxval(abs(fast_smoothed%state - smoothed%state)) < 1.0e-12_dp .and. &
      .not. allocated(fast_smoothed%covariance), 'fast smoother means')
   predicted = kfs_predict(model, filtered, 3)
   call check(predicted%info == 0, 'prediction status')
   call check(all(abs(predicted%mean(:, 1) - filtered%a_filt(1, 5)) < 1.0e-12_dp), 'local-level mean')
   call check(predicted%covariance(1, 1, 3) > predicted%covariance(1, 1, 1), 'prediction variance')
   standardized = kfs_standardized_innovations(filtered)
   call check(abs(standardized(1, 1) - 1.0_dp/sqrt(3.0_dp)) < 1.0e-12_dp, &
      'standardized innovation')
   disturbances = kfs_disturbance_smooth(model, filtered)
   call check(disturbances%info == 0, 'disturbance smoother status')
   call check(all(abs(disturbances%state(:, 1)) < huge(1.0_dp)), 'state disturbances')
   call check(all(abs(disturbances%observation(:, 1)) < huge(1.0_dp)), 'observation disturbances')

   model = make_local_level([1.0_dp, 2.0_dp], 1.0_dp, 1.0_dp, &
      0.0_dp, 1.0_dp)
   filtered = kfs_filter(model)
   smoothed = kfs_smooth(model, filtered)
   call check(smoothed%info == 0 .and. &
      maxval(abs(smoothed%state(1, :) - [0.8_dp, 1.4_dp])) < 1.0e-12_dp .and. &
      maxval(abs(smoothed%covariance(1, 1, :) - &
      [0.4_dp, 0.6_dp])) < 1.0e-12_dp .and. &
      abs(smoothed%lag_one_covariance(1, 1, 1) - 0.2_dp) < 1.0e-12_dp .and. &
      abs(smoothed%conditional_matrix(1, 1, 2) - 0.5_dp) < 1.0e-12_dp .and. &
      abs(smoothed%conditional_covariance(1, 1, 2) - 0.5_dp) < 1.0e-12_dp, &
      'analytic lag-one and conditional covariance smoother')
   singular_covariance = 0.0_dp
   singular_covariance(1, 1, :) = [1.0_dp, 2.0_dp]
   singular_lag_covariance = 0.0_dp
   singular_lag_covariance(1, 1, 1) = 0.5_dp
   conditional = kfs_conditional_covariance(singular_covariance, &
      singular_lag_covariance)
   call check(conditional%info == 0 .and. &
      abs(conditional%matrix(1, 1, 2) - 0.5_dp) < 1.0e-12_dp .and. &
      abs(conditional%covariance(1, 1, 2) - 1.75_dp) < 1.0e-12_dp .and. &
      maxval(abs(conditional%matrix(2, :, :))) < 1.0e-14_dp .and. &
      maxval(abs(conditional%covariance(2, :, :))) < 1.0e-14_dp, &
      'singular conditional covariance pseudoinverse')

   model%a1 = 0.0_dp
   model%p1 = 0.0_dp
   model%p1inf = 1.0_dp
   diffuse_filtered = kfs_filter_diffuse(model)
   call check(diffuse_filtered%info == 0, 'diffuse filter status')
   call check(abs(diffuse_filtered%a_filt(1, 1) - y(1)) < 1.0e-12_dp, 'diffuse first state')
   call check(abs(diffuse_filtered%p_filt(1, 1, 1) - 1.0_dp) < 1.0e-12_dp, &
      'diffuse finite covariance')
   print '(a)', 'All kfas_mod tests passed.'

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
end program
