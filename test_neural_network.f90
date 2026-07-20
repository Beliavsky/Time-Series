! SPDX-License-Identifier: MIT
! SPDX-FileComment: Numerical tests for reusable neural-network regression.
program test_neural_network
   use kind_mod, only: dp
   use neural_network_mod, only: neural_network_t, neural_network_fit, &
      neural_network_predict, neural_network_parameter_count
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none

   type(neural_network_t) :: network, multilayer_network, direct_network
   real(dp) :: predictors(81, 1), response(81, 1), new_predictors(3, 1)
   real(dp) :: linear_response(81, 1)
   real(dp), allocatable :: prediction(:, :)
   integer :: observation

   do observation = 1, size(predictors, 1)
      predictors(observation, 1) = -1.0_dp + &
         2.0_dp*real(observation - 1, dp)/real(size(predictors, 1) - 1, dp)
      response(observation, 1) = predictors(observation, 1)**2
      linear_response(observation, 1) = 1.5_dp + &
         2.25_dp*predictors(observation, 1)
   end do
   network = neural_network_fit(predictors, response, 4, 800, 1.0e-7_dp, &
      1.0e-6_dp)
   call check(network%info == 0 .and. ieee_is_finite(network%rss), &
      'neural regression training status')
   call check(network%rss < 0.05_dp*sum((response - &
      sum(response)/real(size(response), dp))**2), &
      'neural regression nonlinear fit')
   call check(neural_network_parameter_count(1, 4, 1) == 13, &
      'network parameter count')
   call check(neural_network_parameter_count(1, 4, 1, .true.) == 14, &
      'direct network parameter count')
   call check(neural_network_parameter_count(1, [4, 3], 1) == 27, &
      'multilayer network parameter count')

   new_predictors(:, 1) = [-0.75_dp, 0.0_dp, 0.75_dp]
   prediction = neural_network_predict(network, new_predictors)
   call check(all(shape(prediction) == [3, 1]) .and. &
      all(ieee_is_finite(prediction)), 'neural regression prediction')
   call check(abs(prediction(1, 1) - prediction(3, 1)) < 0.08_dp .and. &
      prediction(2, 1) < prediction(1, 1), 'neural regression learned curvature')

   multilayer_network = neural_network_fit(predictors, response, [4, 3], &
      1200, 1.0e-7_dp, 1.0e-6_dp)
   prediction = neural_network_predict(multilayer_network, new_predictors)
   call check(multilayer_network%info == 0 .and. &
      multilayer_network%layer_count == 2 .and. &
      all(multilayer_network%hidden_counts == [4, 3]) .and. &
      multilayer_network%rss < 0.08_dp*sum((response - &
      sum(response)/real(size(response), dp))**2), &
      'multilayer neural regression fit')
   call check(all(shape(prediction) == [3, 1]) .and. &
      all(ieee_is_finite(prediction)), 'multilayer neural regression prediction')

   direct_network = neural_network_fit(predictors, linear_response, 1, 500, &
      1.0e-8_dp, 0.0_dp, direct=.true.)
   prediction = neural_network_predict(direct_network, predictors)
   call check(direct_network%info == 0 .and. direct_network%direct .and. &
      all(shape(direct_network%direct_weights) == [1, 1]) .and. &
      maxval(abs(prediction - linear_response)) < 1.0e-4_dp, &
      'jointly optimized direct neural connections')

   print '(a)', 'Neural-network tests passed'

contains

   subroutine check(condition, message)
      !! Stop the test program when a condition fails.
      logical, intent(in) :: condition !! Condition expected to be true.
      character(len=*), intent(in) :: message !! Failure message.

      if (.not. condition) then
         print '(a)', 'FAILED: '//trim(message)
         error stop 1
      end if
   end subroutine check

end program test_neural_network
