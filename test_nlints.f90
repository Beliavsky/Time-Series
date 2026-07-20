! SPDX-License-Identifier: GPL-2.0-or-later
! SPDX-FileComment: Regression tests for the R NlinTS package translation.
program test_nlints
   use kind_mod, only: dp
   use nlints_mod
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none

   type(nlints_varmlp_t) :: network, updated, loaded
   type(nlints_causality_test_t) :: linear_test
   type(nlints_neural_causality_t) :: neural_test
   integer :: categories(5), paired(5), matrix_categories(5, 2)
   real(dp) :: data(120, 2), target(120), source(120), forecast(2)
   real(dp) :: integrated_target(120), integrated_source(120)
   real(dp), allocatable :: forecast_table(:, :), loaded_table(:, :)
   real(dp) :: continuous(60), shifted(60)
   integer :: time, scratch_unit, display_size, save_info, delete_unit

   categories = [3, 2, 4, 4, 3]
   paired = [1, 4, 4, 3, 3]
   matrix_categories(:, 1) = categories
   matrix_categories(:, 2) = paired
   call check(abs(nlints_entropy_discrete(categories) - &
      1.521928094887362_dp) < 1.0e-12_dp, 'discrete Shannon entropy')
   call check(abs(nlints_mutual_information_discrete(categories, categories) - &
      nlints_entropy_discrete(categories)) < 1.0e-12_dp, &
      'identical-variable mutual information')
   call check(abs(nlints_multivariate_information_discrete(matrix_categories) - &
      nlints_mutual_information_discrete(categories, paired)) < 1.0e-12_dp, &
      'two-variable interaction information')
   call check(ieee_is_finite(nlints_transfer_entropy_discrete( &
      categories, paired)), 'discrete transfer entropy')

   do time = 1, size(continuous)
      continuous(time) = sin(0.17_dp*real(time, dp)) + &
         0.013_dp*real(time, dp)
      shifted(time) = continuous(time) + 0.03_dp*cos(0.71_dp*real(time, dp))
   end do
   call check(ieee_is_finite(nlints_entropy_continuous(continuous, 3)) .and. &
      ieee_is_finite(nlints_mutual_information_continuous(continuous, &
      shifted, 3, 1)) .and. &
      ieee_is_finite(nlints_mutual_information_continuous(continuous, &
      shifted, 3, 2)) .and. &
      ieee_is_finite(nlints_transfer_entropy_continuous(continuous, &
      shifted, 1, 1, 3)), 'continuous KNN information estimators')

   source(1) = 0.2_dp
   target(1) = -0.1_dp
   do time = 2, size(target)
      source(time) = 0.72_dp*source(time - 1) + &
         0.12_dp*sin(0.83_dp*real(time, dp))
      target(time) = 0.35_dp*target(time - 1) + &
         0.75_dp*source(time - 1) + 0.08_dp*cos(1.17_dp*real(time, dp))
   end do
   data(:, 1) = target
   data(:, 2) = source
   network = nlints_varmlp_fit(data, 2, [6], iterations=250, &
      learning_rate=0.01_dp, optimizer=nlints_optimizer_adam, batch_size=12, &
      activations=[nlints_activation_tanh, nlints_activation_linear])
   call check(network%info == 0 .and. all(ieee_is_finite(network%rss)) .and. &
      all(shape(network%fitted) == [120, 2]), 'Adam VARNN fit')
   forecast = nlints_varmlp_forecast(network, data)
   call check(size(forecast) == 2 .and. all(ieee_is_finite(forecast)), &
      'VARNN next-step forecast')
   forecast_table = nlints_varmlp_forecast_table(network, data)
   call check(all(shape(forecast_table) == [119, 2]) .and. &
      maxval(abs(forecast_table(:118, :) - network%fitted(3:, :))) < &
      1.0e-12_dp .and. maxval(abs(forecast_table(119, :) - forecast)) < &
      1.0e-12_dp, 'package-compatible rolling forecast table')
   updated = nlints_varmlp_update(network, data, iterations=5, batch_size=20)
   call check(updated%info == 0 .and. updated%iterations == 255, &
      'incremental VARNN update')
   call nlints_varmlp_save(network, 'nlints_test_model.tmp', save_info)
   loaded = nlints_varmlp_load('nlints_test_model.tmp')
   loaded_table = nlints_varmlp_forecast_table(loaded, data)
   call check(save_info == 0 .and. loaded%info == 0 .and. &
      all(loaded%hidden_counts == network%hidden_counts) .and. &
      maxval(abs(loaded_table - forecast_table)) < 1.0e-12_dp, &
      'versioned VARNN save and load round trip')
   open(newunit=delete_unit, file='nlints_test_model.tmp', status='old')
   close(delete_unit, status='delete')

   linear_test = nlints_granger_test(target, source, 1)
   call check(linear_test%info == 0 .and. linear_test%gci > 0.0_dp .and. &
      linear_test%p_value < 0.05_dp, 'linear Granger causality test')
   integrated_source(1) = source(1)
   integrated_target(1) = target(1)
   do time = 2, size(target)
      integrated_source(time) = integrated_source(time - 1) + source(time)
      integrated_target(time) = integrated_target(time - 1) + target(time)
   end do
   linear_test = nlints_granger_test(integrated_target, integrated_source, 1, &
      difference=.true.)
   call check(linear_test%info == 0 .and. &
      ieee_is_finite(linear_test%statistic) .and. &
      ieee_is_finite(linear_test%p_value), &
      'ADF-selected differenced Granger causality test')
   neural_test = nlints_neural_granger_test(target, source, 1, [2], [3], &
      iterations=100, learning_rate=0.01_dp, optimizer=nlints_optimizer_adam, &
      batch_size=12, seed=17)
   call check(neural_test%test%info == 0 .and. &
      ieee_is_finite(neural_test%test%gci), 'neural Granger causality test')

   open(newunit=scratch_unit, status='scratch', action='readwrite')
   call display(network, scratch_unit)
   call display(linear_test, scratch_unit)
   call display(neural_test, scratch_unit)
   inquire(unit=scratch_unit, size=display_size)
   close(scratch_unit)
   call check(display_size > 0, 'NlinTS display methods')

   print '(a)', 'NlinTS tests passed'

contains

   subroutine check(condition, message)
      !! Stop the test program when a condition fails.
      logical, intent(in) :: condition !! Condition expected to be true.
      character(len=*), intent(in) :: message !! Failure message.

      if (.not. condition) then
         print '(a)', 'FAILED: '//message
         error stop 1
      end if
   end subroutine check

end program test_nlints
