! SPDX-License-Identifier: MIT
! SPDX-FileComment: Numerical tests for algorithms translated from R setartree.
program test_setartree
   use kind_mod, only: dp
   use random_mod, only: set_random_seed
   use utils_mod, only: real_vector_t
   use setartree_mod, only: setartree_model_t, setarforest_model_t, &
      setartree_prediction_t, setartree_forecast_t, setartree_stop_error, &
      setartree_fit, setartree_fit_series, setartree_predict, &
      setartree_forecast, setarforest_fit, setarforest_predict, &
      setarforest_fit_series, setarforest_forecast, &
      setartree_fit_categorical, setartree_predict_categorical, &
      setarforest_fit_categorical, setarforest_predict_categorical
   implicit none

   type(setartree_model_t) :: tree, series_tree, ragged_tree, weighted_tree
   type(setartree_model_t) :: categorical_tree
   type(setarforest_model_t) :: forest, ragged_forest, categorical_forest
   type(setartree_prediction_t) :: prediction, forest_prediction
   type(setartree_prediction_t) :: categorical_prediction
   type(setartree_forecast_t) :: forecast, forest_forecast
   real(dp) :: x(120, 2), y(120), weighted_y(120), weights(120)
   real(dp) :: new_x(2, 2)
   real(dp) :: series(40, 2), history(4, 2)
   real(dp) :: categorical_x(120, 1), categorical_y(120)
   real(dp) :: categorical_new_x(3, 1)
   type(real_vector_t) :: ragged(3)
   integer :: categories(120, 1), new_categories(3, 1)
   integer :: observation, tree_index, failures

   failures = 0
   do observation = 1, 120
      x(observation, 1) = -1.0_dp + &
         2.0_dp*real(observation - 1, dp)/119.0_dp
      x(observation, 2) = sin(0.37_dp*real(observation, dp))
      if (x(observation, 1) < 0.0_dp) then
         y(observation) = 0.8_dp*x(observation, 1) + &
            0.1_dp*x(observation, 2)
      else
         y(observation) = -0.5_dp*x(observation, 1) + &
            0.1_dp*x(observation, 2)
      end if
   end do
   tree = setartree_fit(x, y, maximum_depth=3, threshold_count=15, &
      error_threshold=0.02_dp, stopping_criterion=setartree_stop_error)
   call check(tree%info == 0, "tree fit status", failures)
   call check(tree%node_count >= 3 .and. tree%depth >= 1, &
      "tree split", failures)

   new_x = reshape([-0.8_dp, 0.8_dp, 0.2_dp, 0.2_dp], [2, 2])
   prediction = setartree_predict(tree, new_x, [0.8_dp, 0.95_dp])
   call check(prediction%info == 0, "tree prediction status", failures)
   call check(abs(prediction%value(1) + 0.62_dp) < 0.08_dp, &
      "left leaf prediction", failures)
   call check(abs(prediction%value(2) + 0.38_dp) < 0.08_dp, &
      "right leaf prediction", failures)
   call check(all(prediction%lower <= spread(prediction%value, 2, 2)) .and. &
      all(prediction%upper >= spread(prediction%value, 2, 2)), &
      "tree prediction intervals", failures)

   weighted_y = y
   weighted_y(120) = weighted_y(120) + 10.0_dp
   weights = 1.0_dp
   weights(120) = 0.0_dp
   weighted_tree = setartree_fit(x, weighted_y, maximum_depth=3, &
      threshold_count=15, error_threshold=0.02_dp, &
      stopping_criterion=setartree_stop_error, weights=weights)
   call check(weighted_tree%info == 0 .and. &
      abs(weighted_tree%effective_sample_size(1) - 119.0_dp) < 1.0e-12_dp, &
      "weighted effective sample size", failures)
   prediction = setartree_predict(weighted_tree, new_x)
   call check(all(abs(prediction%value - [-0.62_dp, -0.38_dp]) < 0.08_dp), &
      "zero-weight outlier", failures)

   series(1, :) = [1.0_dp, -1.0_dp]
   do observation = 2, 40
      series(observation, 1) = 0.9_dp*series(observation - 1, 1)
      series(observation, 2) = 0.8_dp*series(observation - 1, 2)
   end do
   series_tree = setartree_fit_series(series, 2, maximum_depth=2, &
      stopping_criterion=setartree_stop_error, error_threshold=0.01_dp)
   call check(series_tree%info == 0, "pooled series fit", failures)
   history = series(37:40, :)
   forecast = setartree_forecast(series_tree, history, 3)
   call check(forecast%info == 0, "tree forecast status", failures)
   call check(all(shape(forecast%mean) == [2, 3]), &
      "tree forecast shape", failures)

   ragged(1)%values = series(:, 1)
   ragged(2)%values = series(5:, 2)
   ragged(3)%values = 0.6_dp*series(:31, 1)
   ragged_tree = setartree_fit_series(ragged, 2, maximum_depth=2, &
      stopping_criterion=setartree_stop_error, error_threshold=0.01_dp)
   call check(ragged_tree%info == 0, "ragged pooled series fit", failures)
   forecast = setartree_forecast(ragged_tree, ragged, 2)
   call check(forecast%info == 0 .and. &
      all(shape(forecast%mean) == [3, 2]), "ragged tree forecast", failures)

   call set_random_seed(9182)
   weights = 0.5_dp
   forest = setarforest_fit(x, y, tree_count=4, bagging_fraction=0.8_dp, &
      randomized_parameters=.true., maximum_depth=3, threshold_count=15, &
      significance=0.04_dp, significance_divider=3.0_dp, &
      error_threshold=0.02_dp, stopping_criterion=setartree_stop_error, &
      randomize_significance=.false., &
      randomize_significance_divider=.false., &
      randomize_error_threshold=.true., weights=weights)
   call check(forest%info == 0 .and. size(forest%tree) == 4, &
      "forest fit", failures)
   call check(all(abs(forest%tree%significance - 0.04_dp) < 1.0e-12_dp), &
      "fixed forest significance", failures)
   call check(all(abs(forest%tree%significance_divider - 3.0_dp) < &
      1.0e-12_dp), "fixed forest significance divider", failures)
   call check(all(forest%tree%error_threshold >= 0.001_dp .and. &
      forest%tree%error_threshold < 0.05_dp), &
      "random forest error thresholds", failures)
   do tree_index = 1, size(forest%tree)
      call check(abs(forest%tree(tree_index)%effective_sample_size(1) - &
         0.5_dp*real(forest%tree(tree_index)%sample_count(1), dp)) < &
         1.0e-12_dp, "forest subsample weights", failures)
   end do
   forest_prediction = setarforest_predict(forest, new_x)
   call check(forest_prediction%info == 0, "forest prediction status", failures)
   call check(all(abs(forest_prediction%value - prediction%value) < 0.15_dp), &
      "forest predictions", failures)
   forest_forecast = setarforest_forecast(forest, history, 2)
   call check(forest_forecast%info == 0, "forest forecast status", failures)
   call check(all(shape(forest_forecast%mean) == [2, 2]), &
      "forest forecast shape", failures)
   call set_random_seed(9183)
   ragged_forest = setarforest_fit_series(ragged, 2, tree_count=3, &
      randomized_parameters=.false., maximum_depth=2, &
      stopping_criterion=setartree_stop_error)
   call check(ragged_forest%info == 0, "ragged forest fit", failures)
   forest_forecast = setarforest_forecast(ragged_forest, ragged, 2)
   call check(forest_forecast%info == 0 .and. &
      all(shape(forest_forecast%mean) == [3, 2]), &
      "ragged forest forecast", failures)

   do observation = 1, 120
      categorical_x(observation, 1) = -1.0_dp + &
         2.0_dp*real(observation - 1, dp)/119.0_dp
      categories(observation, 1) = 10 + 10*mod(observation - 1, 3)
      select case (categories(observation, 1))
      case (10)
         categorical_y(observation) = &
            0.5_dp*categorical_x(observation, 1) + 1.0_dp
      case (20)
         categorical_y(observation) = &
            0.5_dp*categorical_x(observation, 1) - 0.5_dp
      case default
         categorical_y(observation) = 0.5_dp*categorical_x(observation, 1)
      end select
   end do
   categorical_tree = setartree_fit_categorical(categorical_x, categories, &
      categorical_y, maximum_depth=0, error_threshold=0.01_dp, &
      stopping_criterion=setartree_stop_error)
   call check(categorical_tree%info == 0 .and. &
      categorical_tree%numerical_feature_count == 1 .and. &
      categorical_tree%feature_count == 3, "categorical tree fit", failures)
   call check(all(categorical_tree%categorical_levels(1)%value == &
      [10, 20, 30]), "first-seen category levels", failures)
   categorical_new_x(:, 1) = 0.4_dp
   new_categories(:, 1) = [10, 20, 30]
   categorical_prediction = setartree_predict_categorical(categorical_tree, &
      categorical_new_x, new_categories)
   call check(categorical_prediction%info == 0 .and. &
      all(abs(categorical_prediction%value - [1.2_dp, -0.3_dp, 0.2_dp]) < &
      1.0e-8_dp), "categorical tree prediction", failures)
   new_categories(3, 1) = 99
   categorical_prediction = setartree_predict_categorical(categorical_tree, &
      categorical_new_x, new_categories)
   call check(categorical_prediction%info == 2, &
      "unseen category rejection", failures)
   call set_random_seed(9184)
   categorical_forest = setarforest_fit_categorical(categorical_x, categories, &
      categorical_y, tree_count=3, randomized_parameters=.false., &
      maximum_depth=0, stopping_criterion=setartree_stop_error, &
      error_threshold=0.01_dp)
   new_categories(:, 1) = [10, 20, 30]
   categorical_prediction = setarforest_predict_categorical( &
      categorical_forest, categorical_new_x, new_categories)
   call check(categorical_forest%info == 0 .and. &
      categorical_prediction%info == 0, "categorical forest", failures)

   if (failures /= 0) error stop 1
   print *, "setartree tests passed"

contains

   subroutine check(condition, label, failure_count)
      !! Record a failed numerical or structural assertion.
      logical, intent(in) :: condition !! Assertion outcome.
      character(len=*), intent(in) :: label !! Assertion description.
      integer, intent(inout) :: failure_count !! Accumulated failure count.

      if (.not. condition) then
         print *, "FAILED: ", trim(label)
         failure_count = failure_count + 1
      end if
   end subroutine check

end program test_setartree
