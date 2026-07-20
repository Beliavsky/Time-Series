! SPDX-License-Identifier: GPL-3.0-only
! SPDX-FileComment: Regression tests for the R nnfor package translation.
program test_nnfor
   use kind_mod, only: dp
   use nnfor_mod
   use random_mod, only: set_random_seed
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan
   implicit none

   type(nnfor_scaled_t) :: scaled
   type(nnfor_elm_fast_t) :: regression
   type(nnfor_elm_fast_t) :: multilayer_regression
   type(nnfor_elm_fast_t) :: orthogonal_regression, skipped_orthogonal_regression
   type(nnfor_elm_fast_t) :: lasso_regression, stepwise_regression
   type(nnfor_output_fit_t) :: lasso_output, stepwise_output
   type(nnfor_elm_model_t) :: elm_model
   type(nnfor_elm_model_t) :: multilayer_elm
   type(nnfor_elm_model_t) :: multilayer_preprocessed_elm
   type(nnfor_elm_model_t) :: random_multilayer_elm
   type(nnfor_elm_model_t) :: random_multilayer_preprocessed_elm
   type(nnfor_elm_model_t) :: orthogonal_multilayer_elm
   type(nnfor_elm_model_t) :: orthogonal_preprocessed_elm
   type(nnfor_elm_model_t) :: refitted_elm, retrained_elm
   type(nnfor_mlp_model_t) :: mlp_model
   type(nnfor_mlp_model_t) :: multilayer_mlp
   type(nnfor_mlp_model_t) :: multilayer_preprocessed_mlp
   type(nnfor_mlp_model_t) :: refitted_mlp, retrained_mlp
   type(nnfor_elm_model_t) :: preprocessed_model
   type(nnfor_preprocessing_t) :: preprocessing
   type(nnfor_lag_selection_t) :: lag_selection
   type(nnfor_hidden_selection_t) :: hidden_selection
   type(nnfor_elm_hidden_selection_t) :: elm_hidden_selection
   type(nnfor_seasonality_t) :: seasonality
   type(nnfor_canova_hansen_t) :: canova_hansen, dummy_canova_hansen
   type(nnfor_canova_hansen_t) :: rejecting_canova_hansen
   type(nnfor_mseason_t) :: multiplicative_test, additive_test
   type(nnfor_trend_t) :: trend
   type(nnfor_difference_selection_t) :: difference_selection
   type(nnfor_elm_auto_t) :: automatic_elm
   type(nnfor_elm_auto_t) :: selected_automatic_elm
   type(nnfor_elm_auto_t) :: multiple_automatic_elm
   type(nnfor_mlp_auto_t) :: automatic_mlp
   type(nnfor_mlp_auto_t) :: multiple_automatic_mlp
   type(nnfor_elm_model_t) :: multiple_elm
   type(nnfor_mlp_model_t) :: multiple_mlp
   type(integer_vector_t) :: exogenous_lags(1)
   type(nnfor_elm_layer_t) :: multilayer_elm_weights(2, 2)
   type(nnfor_elm_layer_t) :: multilayer_preprocessed_elm_weights(2, 1)
   type(nnfor_forecast_t) :: prediction, elm_forecast, mlp_forecast
   type(nnfor_forecast_t) :: elm_thief_forecast, mlp_thief_forecast
   real(dp) :: predictors(60, 2), response(60), weights(3, 4, 3)
   real(dp) :: output_design(60, 4)
   real(dp) :: series(90), elm_weights(3, 5, 2), initial(13, 2)
   real(dp) :: multilayer_initial(31, 1), multilayer_preprocessed_initial(32, 1)
   real(dp) :: seasonal_series(120), exogenous(124, 1)
   real(dp) :: preprocessed_weights(7, 4, 2), expected_future(4)
   real(dp) :: preprocessed_initial(17, 2)
   real(dp) :: multiple_series(124), multiple_weights(7, 3, 2)
   real(dp) :: ch_series(180), rejecting_ch_series(240)
   real(dp) :: multiplicative_series(120), additive_series(120), level(120)
   real(dp) :: seasonal_pattern(12)
   real(dp), allocatable :: restored(:)
   integer :: fold_ids(60)
   integer :: no_differences(0)
   integer :: row, column, member, scratch_unit, display_size

   scaled = nnfor_linscale([-3.0_dp, 1.0_dp, 5.0_dp], -0.8_dp, 0.8_dp)
   restored = nnfor_apply_scale(scaled%values, scaled%scale, .true.)
   call check(scaled%info == 0 .and. &
      maxval(abs(restored - [-3.0_dp, 1.0_dp, 5.0_dp])) < 1.0e-13_dp, &
      "reversible linear scaling")
   call check(abs(nnfor_fast_sigmoid(2.0_dp) - 2.0_dp/3.0_dp) < &
      1.0e-14_dp, "fast sigmoid")

   do row = 1, size(response)
      predictors(row, 1) = sin(0.13_dp*real(row, dp))
      predictors(row, 2) = cos(0.21_dp*real(row, dp))
      response(row) = 1.2_dp + 2.0_dp*predictors(row, 1) - &
         0.7_dp*predictors(row, 2)
   end do
   do member = 1, 3
      do column = 1, 4
         weights(:, column, member) = [0.03_dp*real(column + member, dp), &
            0.07_dp*real(2*column - member, dp), &
            -0.05_dp*real(column + 2*member, dp)]
      end do
   end do
   regression = nnfor_elm_fast_from_weights(response, predictors, weights, &
      estimator=nnfor_estimator_least_squares, &
      combination=nnfor_combine_mean, direct=.true.)
   prediction = nnfor_elm_fast_predict(regression, predictors)
   call check(regression%info == 0 .and. regression%mse < 1.0e-18_dp .and. &
      maxval(abs(prediction%mean - response)) < 1.0e-8_dp .and. &
      all(shape(prediction%all_mean) == [60, 3]), &
      "fast ELM ensemble with direct connections")

   do member = 1, 2
      allocate(multilayer_elm_weights(1, member)%input_weights(3, 4))
      allocate(multilayer_elm_weights(2, member)%input_weights(5, 3))
      do column = 1, 4
         multilayer_elm_weights(1, member)%input_weights(:, column) = &
            0.05_dp*sin([(real(row + column + member, dp), row=1, 3)])
      end do
      do column = 1, 3
         multilayer_elm_weights(2, member)%input_weights(:, column) = &
            0.06_dp*cos([(real(row + 2*column + member, dp), row=1, 5)])
      end do
   end do
   multilayer_regression = nnfor_elm_fast_layers_from_weights(response, &
      predictors, multilayer_elm_weights, &
      estimator=nnfor_estimator_least_squares, &
      combination=nnfor_combine_mean, direct=.true.)
   prediction = nnfor_elm_fast_predict(multilayer_regression, predictors)
   call check(multilayer_regression%info == 0 .and. &
      all(multilayer_regression%members(1)%hidden_counts == [4, 3]) .and. &
      maxval(abs(prediction%mean - response)) < 1.0e-7_dp, &
      "multilayer fast ELM with direct connections")
   call set_random_seed(90127)
   orthogonal_regression = nnfor_elm_fast(response, predictors, hidden_count=3, &
      repetitions=2, estimator=nnfor_estimator_least_squares, &
      combination=nnfor_combine_mean, orthogonal=.true.)
   call check(orthogonal_regression%info == 0 .and. &
      orthogonal_regression%orthogonal .and. &
      orthogonal_regression%orthogonalized_layer_count == 1 .and. &
      columns_orthonormal(orthogonal_regression%members(1)%input_weights, &
      1.0e-12_dp), "orthogonal fast ELM initialization")
   call set_random_seed(90127)
   skipped_orthogonal_regression = nnfor_elm_fast(response, predictors, &
      hidden_count=4, repetitions=1, estimator=nnfor_estimator_least_squares, &
      orthogonal=.true.)
   call check(skipped_orthogonal_regression%info == 0 .and. &
      .not. skipped_orthogonal_regression%orthogonal .and. &
      skipped_orthogonal_regression%orthogonalized_layer_count == 0, &
      "infeasible orthogonal ELM initialization is skipped")

   elm_hidden_selection = nnfor_select_elm_hidden_from_weights(response, &
      predictors, weights)
   call check(elm_hidden_selection%info == 0 .and. &
      elm_hidden_selection%selected >= 1 .and. &
      elm_hidden_selection%selected <= 4 .and. &
      elm_hidden_selection%candidate_hidden == 4 .and. &
      elm_hidden_selection%repetitions == 3 .and. &
      all(elm_hidden_selection%significant_count >= 0) .and. &
      all(elm_hidden_selection%significant_count <= 4), &
      "supplied-weight automatic ELM hidden-size heuristic")
   call set_random_seed(42017)
   elm_hidden_selection = nnfor_select_elm_hidden(response, predictors, &
      repetitions=3, maximum_hidden=4)
   call check(elm_hidden_selection%info == 0 .and. &
      elm_hidden_selection%selected >= 1 .and. &
      elm_hidden_selection%selected <= 4 .and. &
      size(elm_hidden_selection%significant_count) == 3, &
      "shared-RNG automatic ELM hidden-size heuristic")

   output_design(:, 1) = 1.0_dp
   output_design(:, 2) = predictors(:, 1)
   output_design(:, 3) = predictors(:, 2)
   output_design(:, 4) = sin(0.77_dp*[(real(row, dp), row=1, 60)])
   lasso_output = nnfor_lasso_output_fit(output_design, response, [0.05_dp], &
      0.7_dp)
   stepwise_output = nnfor_stepwise_output_fit(output_design, response)
   call check(lasso_output%info == 0 .and. &
      abs(lasso_output%coefficients(4)) < 1.0e-12_dp .and. &
      stepwise_output%info == 0 .and. .not. stepwise_output%active(4), &
      "lasso and backward-AIC output selection")
   restored = nnfor_combine(reshape([0.9_dp, 0.95_dp, 1.0_dp, 4.0_dp], &
      [1, 4]), nnfor_combine_mode)
   call check(abs(nnfor_kde_mode([0.9_dp, 0.95_dp, 1.0_dp, 4.0_dp]) - &
      0.95_dp) < 0.2_dp .and. &
      abs(restored(1) - 0.95_dp) < 0.2_dp, &
      "KDE-mode ensemble combination")

   lasso_regression = nnfor_elm_fast_from_weights(response, predictors, &
      weights, estimator=nnfor_estimator_lasso, &
      combination=nnfor_combine_mode, direct=.true., lambdas=[0.001_dp])
   stepwise_regression = nnfor_elm_fast_from_weights(response, predictors, &
      weights, estimator=nnfor_estimator_stepwise, &
      combination=nnfor_combine_mean, direct=.true.)
   call check(lasso_regression%info == 0 .and. &
      stepwise_regression%info == 0 .and. &
      lasso_regression%mse >= 0.0_dp .and. stepwise_regression%mse >= 0.0_dp, &
      "ELM lasso and stepwise estimator dispatch")

   series(1:2) = [0.2_dp, 0.25_dp]
   do row = 3, size(series)
      series(row) = 0.15_dp + 0.65_dp*series(row - 1) - &
         0.2_dp*series(row - 2) + 0.01_dp*sin(0.3_dp*real(row, dp))
   end do
   do member = 1, 2
      do column = 1, 5
         elm_weights(:, column, member) = &
            [0.02_dp*real(column, dp), 0.04_dp*real(column + member, dp), &
            -0.03_dp*real(2*column + member, dp)]
      end do
   end do
   elm_model = nnfor_elm_from_weights(series, [1, 2], elm_weights, &
      estimator=nnfor_estimator_ridge, combination=nnfor_combine_median, &
      direct=.true., lambdas=[0.0_dp, 0.001_dp, 0.01_dp])
   elm_forecast = nnfor_elm_forecast(elm_model, 4)
   call check(elm_model%info == 0 .and. size(elm_model%fitted) == 88 .and. &
      elm_model%mse >= 0.0_dp .and. size(elm_forecast%mean) == 4 .and. &
      all(shape(elm_forecast%all_mean) == [4, 2]), &
      "lagged ELM fitting and recursive forecast")
   multilayer_elm = nnfor_elm_layers_from_weights(series, [1, 2], &
      multilayer_elm_weights, estimator=nnfor_estimator_ridge, &
      combination=nnfor_combine_mean, direct=.true., lambdas=[0.0_dp, 0.01_dp])
   elm_forecast = nnfor_elm_forecast(multilayer_elm, 4)
   call check(multilayer_elm%info == 0 .and. elm_forecast%info == 0 .and. &
      size(elm_forecast%mean) == 4 .and. &
      all(abs(elm_forecast%mean) < huge(1.0_dp)), &
      "multilayer lagged ELM fit and recursive forecast")
   call set_random_seed(51803)
   random_multilayer_elm = nnfor_elm_layers(series, [1, 2], [3, 2], &
      repetitions=1, estimator=nnfor_estimator_least_squares, &
      combination=nnfor_combine_mean, direct=.true.)
   call check(random_multilayer_elm%info == 0 .and. &
      all(random_multilayer_elm%network%members(1)%hidden_counts == [3, 2]), &
      "shared-RNG multilayer lagged ELM")
   call set_random_seed(90127)
   orthogonal_multilayer_elm = nnfor_elm_layers(series, [1, 2], [3, 2], &
      repetitions=1, estimator=nnfor_estimator_least_squares, &
      combination=nnfor_combine_mean, orthogonal=.true.)
   call check(orthogonal_multilayer_elm%info == 0 .and. &
      orthogonal_multilayer_elm%network%orthogonal .and. &
      orthogonal_multilayer_elm%network%orthogonalized_layer_count == 2 .and. &
      columns_orthonormal(orthogonal_multilayer_elm%network%members(1)% &
      layers(1)%input_weights, 1.0e-12_dp) .and. &
      columns_orthonormal(orthogonal_multilayer_elm%network%members(1)% &
      layers(2)%input_weights, 1.0e-12_dp), &
      "orthogonal multilayer ELM initialization")
   preprocessing = nnfor_preprocess(series, [1, 2, 3, 4, 5])
   lag_selection = nnfor_select_lags(preprocessing, &
      [.true., .false., .false., .false., .false.])
   call check(lag_selection%info == 0 .and. &
      any(lag_selection%selected_lags == 1) .and. &
      size(lag_selection%selected_lags) >= 1 .and. &
      size(lag_selection%selected_lags) <= 5, &
      "backward-AIC lag selection with forced lag")

   do member = 1, 2
      do row = 1, 13
         initial(row, member) = 0.12_dp*sin(0.37_dp*real(row + 3*member, dp))
      end do
   end do
   mlp_model = nnfor_mlp_from_initial(series, [1, 2], 3, initial, &
      combination=nnfor_combine_mean, max_iterations=300, tolerance=1.0e-6_dp)
   mlp_forecast = nnfor_mlp_forecast(mlp_model, 3)
   call check(mlp_model%info == 0 .and. mlp_model%repetitions == 2 .and. &
      size(mlp_model%fitted) == 88 .and. size(mlp_forecast%mean) == 3 .and. &
      all(abs(mlp_forecast%mean) < huge(1.0_dp)), &
      "lagged MLP ensemble and recursive forecast")

   do row = 1, size(multilayer_initial, 1)
      multilayer_initial(row, 1) = 0.1_dp*sin(0.31_dp*real(row, dp))
   end do
   multilayer_mlp = nnfor_mlp_layers_from_initial(series, [1, 2], [4, 3], &
      multilayer_initial, combination=nnfor_combine_mean, max_iterations=500, &
      tolerance=1.0e-6_dp)
   mlp_forecast = nnfor_mlp_forecast(multilayer_mlp, 3)
   call check(multilayer_mlp%info == 0 .and. &
      all(multilayer_mlp%hidden_counts == [4, 3]) .and. &
      multilayer_mlp%members(1)%layer_count == 2 .and. &
      mlp_forecast%info == 0 .and. all(abs(mlp_forecast%mean) < huge(1.0_dp)), &
      "multilayer lagged MLP fit and recursive forecast")

   do row = 1, size(exogenous, 1)
      exogenous(row, 1) = cos(0.17_dp*real(row, dp))
      if (row <= size(seasonal_series)) then
         seasonal_series(row) = 5.0_dp + 0.3_dp*real(row, dp) + &
            2.0_dp*sin(2.0_dp*acos(-1.0_dp)*real(row, dp)/12.0_dp)
      else
         expected_future(row - size(seasonal_series)) = 5.0_dp + &
            0.3_dp*real(row, dp) + 2.0_dp*sin(2.0_dp*acos(-1.0_dp)* &
            real(row, dp)/12.0_dp)
      end if
   end do
   exogenous_lags(1)%values = [0, 1]
   preprocessing = nnfor_preprocess(seasonal_series, [1, 2], [1, 12], &
      period=12, seasonal_type=nnfor_seasonal_trigonometric, &
      exogenous=exogenous(:120, :), exogenous_lags=exogenous_lags)
   call check(preprocessing%info == 0 .and. &
      all(shape(preprocessing%predictors) == [105, 6]) .and. &
      preprocessing%start_index == 16 .and. &
      maxval(abs(preprocessing%response)) < 1.0e-12_dp .and. &
      size(nnfor_difference(seasonal_series, [1, 12])) == 107, &
      "differencing and seasonal-exogenous preprocessing")
   do member = 1, 2
      do column = 1, 4
         preprocessed_weights(:, column, member) = &
            0.04_dp*sin([(real(row + 2*column + member, dp), row=1, 7)])
      end do
   end do
   preprocessed_model = nnfor_elm_preprocessed_from_weights( &
      seasonal_series, [1, 2], preprocessed_weights, &
      difference_lags=[1, 12], period=12, &
      seasonal_type=nnfor_seasonal_trigonometric, &
      exogenous=exogenous(:120, :), exogenous_lags=exogenous_lags, &
      estimator=nnfor_estimator_least_squares, direct=.true.)
   elm_forecast = nnfor_elm_preprocessed_forecast(preprocessed_model, 4, &
      exogenous)
   call check(preprocessed_model%info == 0 .and. &
      preprocessed_model%extended_preprocessing .and. &
      size(preprocessed_model%fitted) == 105 .and. &
      maxval(abs(preprocessed_model%residuals)) < 1.0e-8_dp .and. &
      elm_forecast%info == 0 .and. &
      maxval(abs(elm_forecast%mean - expected_future)) < 1.0e-7_dp, &
      "preprocessed ELM fit and inverse-differenced forecast")
   refitted_elm = nnfor_elm_refit(preprocessed_model, seasonal_series, &
      exogenous(:120, :))
   call check(refitted_elm%info == 0 .and. &
      maxval(abs(refitted_elm%fitted - preprocessed_model%fitted)) < &
      1.0e-12_dp, "exogenous ELM fixed-weight refit")
   refitted_elm = nnfor_elm_refit(preprocessed_model, seasonal_series)
   call check(refitted_elm%info /= 0, &
      "ELM refit rejects missing required exogenous regressors")

   do member = 1, 2
      do row = 1, 17
         preprocessed_initial(row, member) = &
            0.1_dp*sin(0.29_dp*real(row + 5*member, dp))
      end do
   end do
   mlp_model = nnfor_mlp_preprocessed_from_initial(seasonal_series, [1, 2], &
      2, preprocessed_initial, difference_lags=[1, 12], period=12, &
      seasonal_type=nnfor_seasonal_trigonometric, &
      exogenous=exogenous(:120, :), exogenous_lags=exogenous_lags, &
      combination=nnfor_combine_mean, max_iterations=400, tolerance=1.0e-7_dp)
   mlp_forecast = nnfor_mlp_preprocessed_forecast(mlp_model, 4, exogenous)
   call check(mlp_model%info == 0 .and. mlp_model%extended_preprocessing .and. &
      size(mlp_model%fitted) == 105 .and. mlp_forecast%info == 0 .and. &
      maxval(abs(mlp_forecast%mean - expected_future)) < 1.0e-5_dp, &
      "preprocessed MLP fit and inverse-differenced forecast")

   hidden_selection = nnfor_select_hidden_count(response, predictors, &
      maximum_hidden=3, repetitions=2, validation_weight=0.75_dp, &
      max_iterations=250, tolerance=1.0e-6_dp, &
      combination=nnfor_combine_mean)
   call check(hidden_selection%info == 0 .and. &
      hidden_selection%selected >= 1 .and. &
      hidden_selection%selected <= 3 .and. &
      all(hidden_selection%mse >= 0.0_dp), &
      "validation-based hidden-size selection")

   fold_ids = [(1 + mod(row - 1, 5), row=1, size(fold_ids))]
   hidden_selection = nnfor_select_hidden_count_folds(response, predictors, &
      fold_ids, maximum_hidden=1, repetitions=1, max_iterations=200, &
      tolerance=1.0e-6_dp, combination=nnfor_combine_mean)
   call check(hidden_selection%info == 0 .and. &
      hidden_selection%selected == 1 .and. &
      hidden_selection%mse(1) >= 0.0_dp, &
      "supplied five-fold hidden-size selection")
   call set_random_seed(42017)
   hidden_selection = nnfor_select_hidden_count_random(response, predictors, &
      method=nnfor_hidden_holdout, maximum_hidden=1, repetitions=1, &
      validation_fraction=0.2_dp, max_iterations=200, tolerance=1.0e-6_dp, &
      combination=nnfor_combine_mean)
   call check(hidden_selection%info == 0 .and. &
      hidden_selection%selected == 1, &
      "shared-RNG random-holdout hidden-size selection")
   call set_random_seed(42017)
   hidden_selection = nnfor_select_hidden_count_random(response, predictors, &
      method=nnfor_hidden_cross_validation, maximum_hidden=1, repetitions=1, &
      fold_count=5, max_iterations=200, tolerance=1.0e-6_dp, &
      combination=nnfor_combine_mean)
   call check(hidden_selection%info == 0 .and. &
      hidden_selection%selected == 1, &
      "shared-RNG five-fold hidden-size selection")

   seasonality = nnfor_season_check(seasonal_series, 12)
   ch_series(:12) = sin(0.31_dp*[(real(row, dp), row=1, 12)])
   do row = 13, size(ch_series)
      ch_series(row) = ch_series(row - 12) + &
         0.3_dp*sin(0.17_dp*real(row, dp)) + &
         0.1_dp*cos(0.43_dp*real(row, dp))
   end do
   canova_hansen = nnfor_canova_hansen(ch_series, 12, &
      nnfor_ch_trigonometric)
   dummy_canova_hansen = nnfor_canova_hansen(ch_series, 12, nnfor_ch_dummy)
   call check(canova_hansen%info == 0 .and. &
      abs(canova_hansen%statistic - 0.7450012_dp) < 2.0e-6_dp .and. &
      abs(canova_hansen%critical_value - 2.218139_dp) < 2.0e-6_dp .and. &
      canova_hansen%newey_west_order == 14 .and. &
      .not. canova_hansen%difference_required .and. &
      dummy_canova_hansen%info == 0 .and. &
      abs(dummy_canova_hansen%statistic - 0.7725427_dp) < 2.0e-6_dp, &
      "Canova-Hansen statistics agree with uroot")
   do row = 1, 12
      rejecting_ch_series(row) = sin(0.7_dp*real(row, dp)**2) + &
         0.3_dp*cos(0.87_dp*real(row, dp)**2)
   end do
   do row = 13, size(rejecting_ch_series)
      rejecting_ch_series(row) = rejecting_ch_series(row - 12) + &
         sin(0.7_dp*real(row, dp)**2) + &
         0.3_dp*cos(0.87_dp*real(row, dp)**2)
   end do
   rejecting_canova_hansen = nnfor_canova_hansen(rejecting_ch_series, 12, &
      nnfor_ch_trigonometric)
   call check(rejecting_canova_hansen%info == 0 .and. &
      abs(rejecting_canova_hansen%statistic - 2.478181_dp) < 2.0e-6_dp .and. &
      rejecting_canova_hansen%difference_required, &
      "Canova-Hansen rejection agrees with uroot")
   seasonal_pattern = [-0.18_dp, -0.12_dp, -0.08_dp, -0.04_dp, &
      -0.01_dp, 0.02_dp, 0.04_dp, 0.06_dp, 0.08_dp, 0.10_dp, &
      0.12_dp, 0.16_dp]
   do row = 1, size(level)
      level(row) = 20.0_dp + 0.08_dp*real(row, dp)
      multiplicative_series(row) = level(row)*(1.0_dp + &
         seasonal_pattern(1 + mod(row - 1, 12))) + &
         0.03_dp*sin(0.37_dp*real(row, dp))
      additive_series(row) = level(row) + &
         5.0_dp*seasonal_pattern(1 + mod(row - 1, 12)) + &
         0.2_dp*sin(0.37_dp*real(row, dp))
   end do
   multiplicative_test = nnfor_mseason_test(multiplicative_series, 12, level)
   additive_test = nnfor_mseason_test(additive_series, 12, level)
   call check(multiplicative_test%info == 0 .and. &
      abs(multiplicative_test%statistic - 0.9991259_dp) < 2.0e-7_dp .and. &
      abs(multiplicative_test%p_value - 1.275541e-12_dp) < 1.0e-16_dp .and. &
      multiplicative_test%multiplicative .and. additive_test%info == 0 .and. &
      abs(additive_test%statistic + 0.022075_dp) < 2.0e-6_dp .and. &
      additive_test%p_value == 1.0_dp .and. &
      .not. additive_test%multiplicative, &
      "multiplicative seasonality test agrees with tsutils")
   trend = nnfor_trend_check([(2.0_dp + 0.4_dp*real(row, dp), &
      row=1, 60)])
   difference_selection = nnfor_select_differences(seasonal_series, 12)
   call check(seasonality%info == 0 .and. seasonality%seasonal .and. &
      seasonality%p_value < 0.05_dp .and. trend%info == 0 .and. &
      trend%trending .and. difference_selection%info == 0 .and. &
      all(difference_selection%difference_lags == [1]) .and. &
      .not. difference_selection%canova_hansen%difference_required, &
      "automatic trend and seasonal differencing selection")

   call set_random_seed(73129)
   automatic_elm = nnfor_elm_auto(series, candidate_lags=[1, 2, 3, 4], &
      keep=[.true., .false., .false., .false.], &
      difference_lags=no_differences, repetitions=2, &
      selection_repetitions=3, estimator=nnfor_estimator_least_squares, &
      combination=nnfor_combine_mean, direct=.true.)
   elm_forecast = nnfor_elm_auto_forecast(automatic_elm, 3)
   call check(automatic_elm%info == 0 .and. &
      automatic_elm%model%info == 0 .and. &
      .not. automatic_elm%automatic_differences .and. &
      any(automatic_elm%lag_selection%selected_lags == 1) .and. &
      automatic_elm%hidden_selection%selected >= 1 .and. &
      elm_forecast%info == 0 .and. size(elm_forecast%mean) == 3, &
      "end-to-end automatic ELM specification and forecast")

   call set_random_seed(73129)
   selected_automatic_elm = nnfor_elm_auto(seasonal_series, period=12, &
      candidate_lags=[1, 2], repetitions=1, selection_repetitions=2, &
      estimator=nnfor_estimator_least_squares, &
      combination=nnfor_combine_mean, direct=.true.)
   elm_forecast = nnfor_elm_auto_forecast(selected_automatic_elm, 4)
   call check(selected_automatic_elm%info == 0 .and. &
      selected_automatic_elm%automatic_differences .and. &
      all(selected_automatic_elm%difference_lags == [1]) .and. &
      elm_forecast%info == 0 .and. &
      maxval(abs(elm_forecast%mean - expected_future)) < 1.0e-7_dp, &
      "automatic differencing carried through ELM fitting and forecast")

   call set_random_seed(73129)
   automatic_mlp = nnfor_mlp_auto(series, candidate_lags=[1, 2, 3, 4], &
      keep=[.true., .false., .false., .false.], &
      difference_lags=no_differences, repetitions=1, &
      hidden_method=nnfor_hidden_terminal, maximum_hidden=1, &
      selection_repetitions=1, combination=nnfor_combine_mean, &
      max_iterations=300, tolerance=1.0e-6_dp)
   mlp_forecast = nnfor_mlp_auto_forecast(automatic_mlp, 3)
   call check(automatic_mlp%info == 0 .and. &
      automatic_mlp%model%info == 0 .and. &
      .not. automatic_mlp%automatic_differences .and. &
      any(automatic_mlp%lag_selection%selected_lags == 1) .and. &
      automatic_mlp%hidden_selection%selected == 1 .and. &
      mlp_forecast%info == 0 .and. size(mlp_forecast%mean) == 3, &
      "end-to-end automatic MLP specification and forecast")

   call set_random_seed(73129)
   elm_thief_forecast = nnfor_elm_thief(series, horizon=3, &
      candidate_lags=[1, 2, 3, 4], &
      keep=[.true., .false., .false., .false.], &
      difference_lags=no_differences, repetitions=1, &
      selection_repetitions=1, estimator=nnfor_estimator_least_squares, &
      combination=nnfor_combine_mean, direct=.true.)
   call set_random_seed(73129)
   mlp_thief_forecast = nnfor_mlp_thief(series, horizon=3, &
      candidate_lags=[1, 2, 3, 4], &
      keep=[.true., .false., .false., .false.], &
      difference_lags=no_differences, repetitions=1, &
      hidden_method=nnfor_hidden_terminal, maximum_hidden=1, &
      selection_repetitions=1, combination=nnfor_combine_mean, &
      max_iterations=300, tolerance=1.0e-6_dp)
   call check(elm_thief_forecast%info == 0 .and. &
      mlp_thief_forecast%info == 0 .and. &
      size(elm_thief_forecast%mean) == 3 .and. &
      size(mlp_thief_forecast%mean) == 3 .and. &
      size(elm_thief_forecast%fitted) == size(series) .and. &
      size(mlp_thief_forecast%fitted) == size(series) .and. &
      ieee_is_nan(elm_thief_forecast%fitted(1)) .and. &
      ieee_is_nan(mlp_thief_forecast%fitted(1)) .and. &
      ieee_is_finite(elm_thief_forecast%fitted(size(series))) .and. &
      ieee_is_finite(mlp_thief_forecast%fitted(size(series))), &
      "ELM and MLP temporal-hierarchy callback adapters")

   do row = 1, size(multiple_series)
      multiple_series(row) = 4.0_dp + 0.08_dp*real(row, dp) + &
         0.7_dp*sin(2.0_dp*acos(-1.0_dp)*real(row, dp)/4.0_dp) + &
         1.3_dp*cos(2.0_dp*acos(-1.0_dp)*real(row, dp)/12.0_dp)
   end do
   preprocessing = nnfor_preprocess(multiple_series(:120), [1, 2], [1], &
      periods=[4, 12], seasonal_types=[nnfor_seasonal_trigonometric, &
      nnfor_seasonal_trigonometric])
   call check(preprocessing%info == 0 .and. &
      all(preprocessing%periods == [4, 12]) .and. &
      all(preprocessing%seasonal_types == nnfor_seasonal_trigonometric) .and. &
      all(shape(preprocessing%predictors) == [117, 6]), &
      "multiple-season preprocessing matrix")
   do member = 1, 2
      do column = 1, 3
         multiple_weights(:, column, member) = &
            0.04_dp*sin([(real(row + 2*column + member, dp), row=1, 7)])
      end do
   end do
   multiple_elm = nnfor_elm_preprocessed_from_weights( &
      multiple_series(:120), [1, 2], multiple_weights, difference_lags=[1], &
      estimator=nnfor_estimator_least_squares, &
      combination=nnfor_combine_mean, direct=.true., periods=[4, 12], &
      seasonal_types=[nnfor_seasonal_trigonometric, &
      nnfor_seasonal_trigonometric])
   elm_forecast = nnfor_elm_preprocessed_forecast(multiple_elm, 4)
   call check(multiple_elm%info == 0 .and. elm_forecast%info == 0 .and. &
      maxval(abs(elm_forecast%mean - multiple_series(121:))) < 1.0e-7_dp, &
      "multiple-season ELM fit and recursive forecast")
   allocate(multilayer_preprocessed_elm_weights(1, 1)%input_weights(7, 3))
   allocate(multilayer_preprocessed_elm_weights(2, 1)%input_weights(4, 2))
   do column = 1, 3
      multilayer_preprocessed_elm_weights(1, 1)%input_weights(:, column) = &
         0.04_dp*sin([(real(row + column, dp), row=1, 7)])
   end do
   do column = 1, 2
      multilayer_preprocessed_elm_weights(2, 1)%input_weights(:, column) = &
         0.05_dp*cos([(real(row + 2*column, dp), row=1, 4)])
   end do
   multilayer_preprocessed_elm = &
      nnfor_elm_preprocessed_layers_from_weights(multiple_series(:120), &
      [1, 2], multilayer_preprocessed_elm_weights, difference_lags=[1], &
      estimator=nnfor_estimator_least_squares, &
      combination=nnfor_combine_mean, direct=.true., periods=[4, 12], &
      seasonal_types=[nnfor_seasonal_trigonometric, &
      nnfor_seasonal_trigonometric])
   elm_forecast = nnfor_elm_preprocessed_forecast(multilayer_preprocessed_elm, 4)
   call check(multilayer_preprocessed_elm%info == 0 .and. &
      multilayer_preprocessed_elm%extended_preprocessing .and. &
      elm_forecast%info == 0 .and. &
      maxval(abs(elm_forecast%mean - multiple_series(121:))) < 1.0e-7_dp, &
      "preprocessed multilayer ELM fit and recursive forecast")
   refitted_elm = nnfor_elm_refit(multilayer_preprocessed_elm, &
      multiple_series(:120))
   call check(refitted_elm%info == 0 .and. &
      maxval(abs(refitted_elm%fitted - &
      multilayer_preprocessed_elm%fitted)) < 1.0e-12_dp .and. &
      all(refitted_elm%network%members(1)%hidden_counts == [3, 2]), &
      "fixed-weight multilayer ELM refit")
   call set_random_seed(81239)
   retrained_elm = nnfor_elm_retrain(multilayer_preprocessed_elm, &
      multiple_series(:120))
   call check(retrained_elm%info == 0 .and. &
      retrained_elm%extended_preprocessing .and. &
      all(retrained_elm%preprocessing%periods == [4, 12]) .and. &
      all(retrained_elm%network%members(1)%hidden_counts == [3, 2]), &
      "multilayer ELM specification-preserving retraining")
   call set_random_seed(51803)
   random_multilayer_preprocessed_elm = nnfor_elm_preprocessed_layers( &
      multiple_series(:120), [1, 2], [3, 2], repetitions=1, &
      difference_lags=[1], estimator=nnfor_estimator_least_squares, &
      combination=nnfor_combine_mean, direct=.true., periods=[4, 12], &
      seasonal_types=[nnfor_seasonal_trigonometric, &
      nnfor_seasonal_trigonometric])
   call check(random_multilayer_preprocessed_elm%info == 0 .and. &
      random_multilayer_preprocessed_elm%extended_preprocessing, &
      "shared-RNG preprocessed multilayer ELM")
   call set_random_seed(90127)
   orthogonal_preprocessed_elm = nnfor_elm_preprocessed( &
      multiple_series(:120), [1, 2], hidden_count=3, repetitions=1, &
      difference_lags=[1], estimator=nnfor_estimator_least_squares, &
      combination=nnfor_combine_mean, direct=.true., periods=[4, 12], &
      seasonal_types=[nnfor_seasonal_trigonometric, &
      nnfor_seasonal_trigonometric], orthogonal=.true.)
   call check(orthogonal_preprocessed_elm%info == 0 .and. &
      orthogonal_preprocessed_elm%network%orthogonal .and. &
      columns_orthonormal(orthogonal_preprocessed_elm%network%members(1)% &
      input_weights, 1.0e-12_dp), &
      "orthogonal preprocessed ELM initialization")
   multiple_mlp = nnfor_mlp_preprocessed_from_initial( &
      multiple_series(:120), [1, 2], 2, preprocessed_initial, &
      difference_lags=[1], combination=nnfor_combine_mean, &
      max_iterations=400, tolerance=1.0e-7_dp, periods=[4, 12], &
      seasonal_types=[nnfor_seasonal_trigonometric, &
      nnfor_seasonal_trigonometric])
   mlp_forecast = nnfor_mlp_preprocessed_forecast(multiple_mlp, 4)
   call check(multiple_mlp%info == 0 .and. mlp_forecast%info == 0 .and. &
      all(abs(mlp_forecast%mean) < huge(1.0_dp)), &
      "multiple-season MLP fit and recursive forecast")
   do row = 1, size(multilayer_preprocessed_initial, 1)
      multilayer_preprocessed_initial(row, 1) = &
         0.08_dp*cos(0.23_dp*real(row, dp))
   end do
   multilayer_preprocessed_mlp = &
      nnfor_mlp_preprocessed_layers_from_initial(multiple_series(:120), &
      [1, 2], [3, 2], multilayer_preprocessed_initial, difference_lags=[1], &
      combination=nnfor_combine_mean, max_iterations=500, tolerance=1.0e-6_dp, &
      periods=[4, 12], seasonal_types=[nnfor_seasonal_trigonometric, &
      nnfor_seasonal_trigonometric])
   mlp_forecast = nnfor_mlp_preprocessed_forecast(multilayer_preprocessed_mlp, 4)
   call check(multilayer_preprocessed_mlp%info == 0 .and. &
      all(multilayer_preprocessed_mlp%hidden_counts == [3, 2]) .and. &
      multilayer_preprocessed_mlp%extended_preprocessing .and. &
      mlp_forecast%info == 0 .and. all(abs(mlp_forecast%mean) < huge(1.0_dp)), &
      "preprocessed multilayer MLP fit and recursive forecast")
   refitted_mlp = nnfor_mlp_refit(multilayer_preprocessed_mlp, &
      multiple_series(:120))
   call check(refitted_mlp%info == 0 .and. &
      maxval(abs(refitted_mlp%fitted - &
      multilayer_preprocessed_mlp%fitted)) < 1.0e-10_dp .and. &
      all(refitted_mlp%hidden_counts == [3, 2]), &
      "fixed-weight multilayer MLP refit")
   call set_random_seed(81239)
   retrained_mlp = nnfor_mlp_retrain(multilayer_preprocessed_mlp, &
      multiple_series(:120), max_iterations=500, tolerance=1.0e-6_dp)
   call check(retrained_mlp%info == 0 .and. &
      retrained_mlp%extended_preprocessing .and. &
      all(retrained_mlp%preprocessing%periods == [4, 12]) .and. &
      all(retrained_mlp%hidden_counts == [3, 2]), &
      "multilayer MLP specification-preserving retraining")
   call set_random_seed(73129)
   multiple_automatic_elm = nnfor_elm_auto(multiple_series(:120), &
      candidate_lags=[1, 2], difference_lags=[1], repetitions=1, &
      selection_repetitions=2, estimator=nnfor_estimator_least_squares, &
      combination=nnfor_combine_mean, direct=.true., periods=[4, 12])
   call check(multiple_automatic_elm%info == 0 .and. &
      all(multiple_automatic_elm%periods == [4, 12]) .and. &
      all(multiple_automatic_elm%seasonal_types == &
      nnfor_seasonal_trigonometric), &
      "automatic multiple-season ELM specification")
   call set_random_seed(73129)
   multiple_automatic_mlp = nnfor_mlp_auto(multiple_series(:120), &
      candidate_lags=[1, 2], difference_lags=[1], repetitions=1, &
      hidden_method=nnfor_hidden_terminal, maximum_hidden=1, &
      selection_repetitions=1, combination=nnfor_combine_mean, &
      max_iterations=300, tolerance=1.0e-6_dp, periods=[4, 12])
   mlp_forecast = nnfor_mlp_auto_forecast(multiple_automatic_mlp, 4)
   call check(multiple_automatic_mlp%info == 0 .and. &
      all(multiple_automatic_mlp%seasonal_types == &
      nnfor_seasonal_trigonometric) .and. mlp_forecast%info == 0 .and. &
      all(abs(mlp_forecast%mean) < huge(1.0_dp)), &
      "automatic multiple-season MLP specification and forecast")

   open(newunit=scratch_unit, status='scratch', action='write')
   call display(regression, scratch_unit)
   call display(elm_model, scratch_unit)
   call display(mlp_model, scratch_unit)
   call display(mlp_forecast, scratch_unit)
   call display(automatic_elm, scratch_unit)
   call display(selected_automatic_elm, scratch_unit)
   call display(automatic_mlp, scratch_unit)
   call display(multiple_elm, scratch_unit)
   call display(multilayer_elm, scratch_unit)
   call display(multilayer_preprocessed_elm, scratch_unit)
   call display(multiple_mlp, scratch_unit)
   call display(multilayer_mlp, scratch_unit)
   call display(multilayer_preprocessed_mlp, scratch_unit)
   call display(multiple_automatic_elm, scratch_unit)
   call display(multiple_automatic_mlp, scratch_unit)
   inquire(unit=scratch_unit, size=display_size)
   close(scratch_unit)
   call check(display_size > 0, "nnfor display methods")

   print '(a)', "nnfor tests passed"

contains

   pure logical function columns_orthonormal(matrix, tolerance) result(valid)
      !! Test whether matrix columns have an identity cross-product.
      real(dp), intent(in) :: matrix(:, :) !! Matrix whose columns are tested.
      real(dp), intent(in) :: tolerance !! Maximum absolute cross-product error.
      real(dp), allocatable :: cross_product(:, :)
      integer :: column

      cross_product = matmul(transpose(matrix), matrix)
      do column = 1, size(cross_product, 1)
         cross_product(column, column) = cross_product(column, column) - 1.0_dp
      end do
      valid = maxval(abs(cross_product)) <= tolerance
   end function columns_orthonormal

   subroutine check(condition, message)
      !! Stop the test program when a condition fails.
      logical, intent(in) :: condition !! Test condition.
      character(len=*), intent(in) :: message !! Failure message.

      if (.not. condition) then
         print '(a)', "FAILED: "//message
         error stop 1
      end if
   end subroutine check

end program test_nnfor
