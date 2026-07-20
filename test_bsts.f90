! SPDX-License-Identifier: MIT
! SPDX-FileComment: Regression tests for the bsts translation.
program test_bsts
   use kind_mod, only: dp
   use bsts_mod
   use calendar_mod, only: date_t
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, &
      ieee_quiet_nan
   implicit none

   type(bsts_mcmc_t) :: level_fit, trend_fit, semilocal_fit
   type(bsts_mcmc_t) :: seasonal_fit, trig_fit, random_fit
   type(bsts_mcmc_t) :: ar_fit
   type(bsts_mcmc_t) :: student_fit
   type(bsts_mcmc_t) :: monthly_fit
   type(bsts_mcmc_t) :: holiday_fit
   type(bsts_holiday_t) :: holiday, fixed_holiday, weekday_holiday
   type(bsts_holiday_t), allocatable :: regression_holidays(:)
   type(bsts_holiday_t), allocatable :: hierarchical_holidays(:)
   type(bsts_holiday_regression_t) :: holiday_regression_fit
   type(bsts_holiday_regression_t) :: hierarchical_holiday_fit
   type(bsts_shared_local_level_t) :: shared_fit
   type(bsts_non_gaussian_t) :: count_fit
   type(bsts_multivariate_prediction_t) :: multivariate_prediction
   type(bsts_prediction_t) :: prediction
   type(bsts_spike_slab_t) :: regression_fit, random_regression_fit
   type(bsts_dirm_t) :: dirm_fit
   type(bsts_mixed_t) :: mixed_fit
   type(bsts_mixed_prediction_t) :: mixed_prediction
   type(bsts_prediction_errors_t) :: prediction_errors, scaled_errors
   type(bsts_prediction_errors_t) :: holdout_errors
   type(bsts_prediction_errors_t), allocatable :: comparison_inputs(:)
   type(bsts_model_comparison_t) :: model_comparison
   type(bsts_static_intercept_t) :: static_intercept_fit
   type(bsts_numeric_timestamps_t) :: numeric_timestamps
   type(bsts_date_timestamps_t) :: date_timestamps
   type(bsts_wide_series_t) :: wide_series
   type(bsts_long_series_t) :: long_series
   type(bsts_monthly_series_t) :: monthly_series
   type(bsts_dynamic_regression_t) :: dynamic_fit, hierarchical_fit
   real(dp), allocatable :: sequence(:), aggregate(:), y(:)
   real(dp), allocatable :: count_response(:), count_scale(:)
   real(dp), allocatable :: count_normals(:, :), count_uniforms(:, :)
   real(dp), allocatable :: count_gamma(:), count_forecast_normals(:, :)
   real(dp), allocatable :: count_forecast_uniforms(:, :)
   real(dp), allocatable :: count_predictors(:, :), coefficient_normals(:, :)
   real(dp), allocatable :: coefficient_uniforms(:, :), birth_normals(:, :)
   real(dp), allocatable :: birth_uniforms(:, :)
   real(dp), allocatable :: structural_normals(:, :, :)
   real(dp), allocatable :: structural_uniforms(:), structural_gamma(:, :)
   real(dp), allocatable :: structural_forecast_normals(:, :, :)
   real(dp), allocatable :: level_normals(:, :, :), level_gamma(:, :)
   real(dp), allocatable :: trend_normals(:, :, :), trend_gamma(:, :)
   real(dp), allocatable :: forecast_state(:, :, :), forecast_observation(:, :)
   real(dp), allocatable :: seasonal_y(:), seasonal_normals(:, :, :)
   real(dp), allocatable :: seasonal_gamma(:, :), trig_y(:)
   real(dp), allocatable :: trig_normals(:, :, :), trig_gamma(:, :)
   real(dp), allocatable :: semilocal_y(:), semilocal_normals(:, :, :)
   real(dp), allocatable :: semilocal_gamma(:, :), parameter_normals(:, :)
   real(dp), allocatable :: ar_uniforms(:)
   real(dp), allocatable :: semi_forecast_state(:, :, :)
   real(dp), allocatable :: semi_forecast_observation(:, :)
   real(dp), allocatable :: regression_y(:), regression_x(:, :)
   real(dp), allocatable :: regression_normals(:, :), inclusion_uniforms(:, :)
   real(dp), allocatable :: regression_gamma(:), regression_offset(:, :)
   real(dp), allocatable :: regression_prediction_normals(:, :)
   real(dp), allocatable :: dirm_y(:), dirm_x(:, :)
   real(dp), allocatable :: dirm_state_normals(:, :, :)
   real(dp), allocatable :: dirm_coefficient_normals(:, :)
   real(dp), allocatable :: dirm_uniforms(:, :), dirm_gamma(:, :)
   real(dp), allocatable :: dirm_forecast_state(:, :)
   real(dp), allocatable :: dirm_forecast_observation(:, :)
   integer, allocatable :: dirm_time(:)
   real(dp), allocatable :: mixed_coarse(:), mixed_x(:, :)
   real(dp), allocatable :: mixed_fraction(:), mixed_latent_normals(:, :)
   real(dp), allocatable :: mixed_state_normals(:, :, :)
   real(dp), allocatable :: mixed_coefficient_normals(:, :)
   real(dp), allocatable :: mixed_uniforms(:, :), mixed_gamma(:, :)
   real(dp), allocatable :: mixed_forecast_state(:, :)
   real(dp), allocatable :: mixed_forecast_observation(:, :)
   real(dp), allocatable :: mixed_structural_normals(:, :, :)
   real(dp), allocatable :: mixed_component_gamma(:, :)
   real(dp), allocatable :: mixed_structural_forecast(:, :, :)
   integer, allocatable :: mixed_coarse_index(:)
   logical, allocatable :: mixed_ends(:)
   real(dp), allocatable :: static_response(:), static_offset(:, :)
   real(dp), allocatable :: static_normals(:), static_gamma(:)
   real(dp), allocatable :: static_forecast_normals(:, :)
   real(dp), allocatable :: dynamic_y(:), dynamic_x(:, :)
   real(dp), allocatable :: dynamic_normals(:, :, :), dynamic_gamma(:, :)
   real(dp), allocatable :: dynamic_forecast_state(:, :, :)
   real(dp), allocatable :: dynamic_forecast_observation(:, :)
   real(dp), allocatable :: ar_dynamic_y(:), ar_dynamic_normals(:, :, :)
   real(dp), allocatable :: ar_parameter_normals(:, :, :, :)
   real(dp), allocatable :: ar_dynamic_gamma(:, :)
   real(dp), allocatable :: ar_state_y(:), ar_state_normals(:, :, :)
   real(dp), allocatable :: ar_state_parameter_normals(:, :, :)
   real(dp), allocatable :: ar_state_gamma(:, :), ar_state_forecast(:, :)
   real(dp), allocatable :: auto_ar_uniforms(:, :)
   real(dp), allocatable :: student_y(:), student_state_normals(:, :, :)
   real(dp), allocatable :: student_weight_normals(:, :, :, :)
   real(dp), allocatable :: student_weight_uniforms(:, :, :, :)
   real(dp), allocatable :: student_variance_gamma(:, :)
   real(dp), allocatable :: student_degrees_normals(:, :)
   real(dp), allocatable :: student_degrees_uniforms(:, :)
   real(dp), allocatable :: student_forecast_state(:, :, :)
   real(dp), allocatable :: student_forecast_observation(:, :)
   real(dp), allocatable :: monthly_y(:), monthly_normals(:, :, :)
   real(dp), allocatable :: monthly_gamma(:, :), monthly_forecast(:, :)
   real(dp), allocatable :: holiday_y(:), holiday_normals(:, :, :)
   real(dp), allocatable :: holiday_gamma(:, :), holiday_forecast(:, :)
   real(dp), allocatable :: holiday_design(:, :), holiday_regression_y(:)
   real(dp), allocatable :: holiday_regression_normals(:, :)
   real(dp), allocatable :: holiday_regression_gamma(:)
   real(dp), allocatable :: holiday_regression_offset(:, :)
   real(dp), allocatable :: holiday_regression_forecast(:, :)
   real(dp), allocatable :: hierarchical_design(:, :), hierarchical_y(:)
   real(dp), allocatable :: hierarchical_coefficient_normals(:, :)
   real(dp), allocatable :: hierarchical_mean_normals(:, :)
   real(dp), allocatable :: hierarchical_wishart_normals(:, :, :)
   real(dp), allocatable :: hierarchical_wishart_gamma(:, :)
   real(dp), allocatable :: hierarchical_residual_gamma(:)
   real(dp), allocatable :: hierarchical_offset(:, :)
   real(dp), allocatable :: hierarchical_forecast(:, :)
   real(dp), allocatable :: shared_response(:, :), shared_factors(:, :)
   real(dp), allocatable :: shared_state_normals(:, :, :)
   real(dp), allocatable :: shared_loading_normals(:, :, :)
   real(dp), allocatable :: shared_gamma(:, :), shared_offset(:, :, :)
   real(dp), allocatable :: shared_forecast_state(:, :, :)
   real(dp), allocatable :: shared_forecast_observation(:, :, :)
   real(dp), allocatable :: shared_inclusion_uniforms(:, :, :)
   real(dp) :: shared_inclusion_prior(3, 2)
   logical :: shared_initial_inclusion(3, 2)
   real(dp) :: shared_loadings(3, 2)
   real(dp), allocatable :: mbsts_response(:, :), mbsts_predictors(:, :, :)
   real(dp), allocatable :: mbsts_state_normals(:, :, :)
   real(dp), allocatable :: mbsts_loading_normals(:, :, :)
   real(dp), allocatable :: mbsts_gamma(:, :), mbsts_offset(:, :, :)
   real(dp), allocatable :: mbsts_regression_normals(:, :, :)
   real(dp), allocatable :: mbsts_regression_uniforms(:, :, :)
   real(dp), allocatable :: mbsts_series_normals(:, :, :)
   real(dp), allocatable :: mbsts_series_gamma(:, :)
   real(dp), allocatable :: mbsts_series_forecast(:, :, :)
   real(dp), allocatable :: mbsts_trend_normals(:, :, :, :)
   real(dp), allocatable :: mbsts_trend_gamma(:, :, :)
   real(dp), allocatable :: mbsts_trend_forecast(:, :, :, :)
   real(dp), allocatable :: mbsts_seasonal_normals(:, :, :, :)
   real(dp), allocatable :: mbsts_seasonal_gamma(:, :)
   real(dp), allocatable :: mbsts_seasonal_forecast(:, :, :, :)
   real(dp), allocatable :: mbsts_future_predictors(:, :, :)
   real(dp), allocatable :: mbsts_forecast_state(:, :, :)
   real(dp), allocatable :: mbsts_forecast_observation(:, :, :)
   real(dp) :: mbsts_beta(2, 3)
   real(dp) :: expected_trend(2), expected_season(3), next_season(3)
   real(dp) :: mbsts_inclusion_prior(2, 3)
   logical :: mbsts_initial_inclusion(2, 3)
   real(dp) :: beta_first, beta_second, previous_first, previous_second
   real(dp) :: seasonal_pattern(8)
   real(dp), allocatable :: aggregate_matrix(:, :)
   type(date_t) :: monthly_dates_input(3), weekly_dates(9)
   integer :: matched_month(1)
   real(dp) :: slope
   logical :: boundaries(4)
   integer :: i, j, k

   sequence = bsts_geometric_sequence(4, 2.0_dp, 0.5_dp)
   call check(maxval(abs(sequence - [2.0_dp, 1.0_dp, 0.5_dp, &
      0.25_dp])) < 1.0e-14_dp, 'geometric sequence')
   boundaries = [.false., .true., .false., .true.]
   aggregate = bsts_harvey_cumulator([10.0_dp, 20.0_dp, 30.0_dp, 40.0_dp], &
      boundaries, [0.5_dp])
   call check(maxval(abs(aggregate - [10.0_dp, 20.0_dp, 40.0_dp, &
      60.0_dp])) < 1.0e-14_dp, 'Harvey cumulator')
   aggregate = bsts_harvey_cumulator([7.0_dp], [.true.], [0.25_dp])
   call check(size(aggregate) == 1 .and. aggregate(1) == 7.0_dp, &
      'Harvey one-observation case')
   numeric_timestamps = bsts_regularize_numeric_timestamps( &
      [1.0_dp, 2.0_dp, 2.0_dp, 4.0_dp])
   call check(numeric_timestamps%info == 0 .and. &
      all(abs(numeric_timestamps%grid - [1.0_dp, 2.0_dp, 3.0_dp, &
      4.0_dp]) < 1.0e-14_dp) .and. &
      all(numeric_timestamps%mapping == [1, 2, 2, 4]) .and. &
      .not. bsts_no_duplicates_numeric([1.0_dp, 2.0_dp, 2.0_dp]) .and. &
      .not. bsts_no_duplicates_numeric([1.0_dp, &
      ieee_value(0.0_dp, ieee_quiet_nan)]) .and. &
      .not. bsts_no_gaps_numeric([1.0_dp, 2.0_dp, 4.0_dp]) .and. &
      bsts_is_regular_numeric([1.0_dp, 2.0_dp, 3.0_dp]), &
      'numeric timestamp regularization')
   monthly_dates_input = [date_t(2024, 1, 15), date_t(2024, 2, 15), &
      date_t(2024, 4, 15)]
   date_timestamps = bsts_regularize_date_timestamps(monthly_dates_input)
   call check(date_timestamps%info == 0 .and. size(date_timestamps%grid) == 4 .and. &
      date_timestamps%grid(3)%year == 2024 .and. &
      date_timestamps%grid(3)%month == 3 .and. &
      date_timestamps%grid(3)%day == 15 .and. &
      all(date_timestamps%mapping == [1, 2, 4]) .and. &
      .not. bsts_no_gaps_date(monthly_dates_input), &
      'date timestamp regularization')
   wide_series = bsts_long_to_wide([1.0_dp, 2.0_dp, 3.0_dp], &
      [10, 20, 10], [1, 1, 2])
   call check(wide_series%info == 0 .and. &
      all(shape(wide_series%values) == [2, 2]) .and. &
      all(wide_series%timestamps == [1, 2]) .and. &
      all(wide_series%series_id == [10, 20]) .and. &
      wide_series%values(1, 1) == 1.0_dp .and. &
      wide_series%values(1, 2) == 2.0_dp .and. &
      wide_series%values(2, 1) == 3.0_dp .and. &
      .not. ieee_is_finite(wide_series%values(2, 2)), &
      'long-to-wide conversion')
   long_series = bsts_wide_to_long(wide_series%values, &
      wide_series%timestamps, wide_series%series_id)
   call check(long_series%info == 0 .and. &
      all(long_series%values == [1.0_dp, 2.0_dp, 3.0_dp]) .and. &
      all(long_series%timestamps == [1, 1, 2]) .and. &
      all(long_series%series_id == [10, 20, 10]), &
      'wide-to-long conversion')
   wide_series = bsts_long_to_wide([ieee_value(0.0_dp, ieee_quiet_nan), &
      2.0_dp], [10, 10], [1, 1])
   call check(wide_series%info == 2, 'duplicate long-format cell')
   aggregate = bsts_aggregate_time_series( &
      [10.0_dp, 20.0_dp, 30.0_dp, 40.0_dp], &
      [.false., .true., .false., .true.], &
      [0.5_dp, 0.5_dp, 0.5_dp, 0.5_dp], .false., .false.)
   call check(all(abs(aggregate - [20.0_dp, 60.0_dp, 20.0_dp]) < &
      1.0e-14_dp), 'general vector aggregation')
   aggregate_matrix = bsts_aggregate_time_series(reshape([10.0_dp, &
      20.0_dp, 30.0_dp, 40.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], &
      [4, 2]), [.false., .true., .false., .true.], &
      [0.5_dp, 0.5_dp, 0.5_dp, 0.5_dp], .true., .true.)
   call check(all(shape(aggregate_matrix) == [1, 2]) .and. &
      maxval(abs(aggregate_matrix(1, :) - [60.0_dp, 6.0_dp])) < 1.0e-14_dp, &
      'general matrix aggregation')
   weekly_dates = [date_t(2024, 1, 6), date_t(2024, 1, 13), &
      date_t(2024, 1, 20), date_t(2024, 1, 27), date_t(2024, 2, 3), &
      date_t(2024, 2, 10), date_t(2024, 2, 17), date_t(2024, 2, 24), &
      date_t(2024, 3, 2)]
   matched_month = bsts_match_week_to_month([weekly_dates(5)], &
      date_t(2024, 1, 1))
   call check(bsts_week_ends_month(weekly_dates(5)) .and. &
      abs(bsts_fraction_initial_month(weekly_dates(5)) - 4.0_dp/7.0_dp) < &
      1.0e-14_dp .and. matched_month(1) == 1, &
      'weekly month boundary helpers')
   monthly_series = bsts_aggregate_weeks_to_months( &
      reshape([(1.0_dp, i=1, 9)], [9, 1]), weekly_dates)
   call check(monthly_series%info == 0 .and. &
      all(shape(monthly_series%values) == [2, 1]) .and. &
      maxval(abs(monthly_series%values(:, 1) - &
      [31.0_dp/7.0_dp, 29.0_dp/7.0_dp])) < 1.0e-14_dp .and. &
      monthly_series%dates(1)%year == 2024 .and. &
      monthly_series%dates(1)%month == 1 .and. &
      monthly_series%dates(1)%day == 31 .and. &
      monthly_series%dates(2)%year == 2024 .and. &
      monthly_series%dates(2)%month == 2 .and. &
      monthly_series%dates(2)%day == 29, 'weekly-to-monthly aggregation')

   allocate(static_response(4), static_offset(4, 6), static_normals(6), &
      static_gamma(6), static_forecast_normals(2, 4))
   static_response = [2.0_dp, 3.0_dp, &
      ieee_value(0.0_dp, ieee_quiet_nan), 4.0_dp]
   static_offset = 1.0_dp
   static_normals = 0.0_dp
   static_gamma = 4.0_dp
   static_forecast_normals = 0.0_dp
   static_intercept_fit = bsts_static_intercept_draws(static_response, &
      1.5_dp, 2.0_dp, 0.5_dp, 2.0_dp, 0.1_dp, 2, static_normals, &
      static_gamma, static_offset)
   call check(static_intercept_fit%info == 0 .and. &
      all(shape(static_intercept_fit%contribution) == [4, 6]) .and. &
      maxval(abs(static_intercept_fit%contribution - &
      spread(static_intercept_fit%intercept, 1, 4))) < 1.0e-14_dp .and. &
      maxval(abs(static_intercept_fit%fitted - &
      static_intercept_fit%contribution - static_offset)) < 1.0e-14_dp .and. &
      all(static_intercept_fit%observation_variance > 0.0_dp), &
      'static-intercept supplied draws and offsets')
   prediction = bsts_static_intercept_predict_draws(static_intercept_fit, &
      [0.5_dp, -0.5_dp], static_forecast_normals)
   call check(prediction%info == 0 .and. &
      all(shape(prediction%draws) == [2, 4]) .and. &
      maxval(abs(prediction%draws(1, :) - prediction%draws(2, :) - &
      1.0_dp)) < 1.0e-14_dp, 'static-intercept supplied forecast')
   static_intercept_fit = bsts_static_intercept(static_response, 4, burn=1, &
      offset_draws=static_offset(:, :4))
   call check(static_intercept_fit%info == 0, &
      'random static-intercept wrapper')
   prediction = bsts_static_intercept_predict(static_intercept_fit, &
      [0.0_dp, 1.0_dp])
   call check(prediction%info == 0 .and. &
      all(shape(prediction%draws) == [2, 3]), &
      'random static-intercept forecast wrapper')

   allocate(count_response(12), count_scale(12), count_normals(12, 6), &
      count_uniforms(12, 6), count_gamma(6), &
      count_forecast_normals(3, 4), count_forecast_uniforms(3, 4), &
      count_predictors(12, 2), coefficient_normals(2, 6), &
      coefficient_uniforms(2, 6), birth_normals(2, 6), &
      birth_uniforms(2, 6), structural_normals(5, 12, 6), &
      structural_uniforms(6), structural_gamma(3, 6), &
      structural_forecast_normals(5, 3, 4))
   count_response = [2.0_dp, 3.0_dp, 1.0_dp, 4.0_dp, 2.0_dp, 3.0_dp, &
      2.0_dp, 1.0_dp, 4.0_dp, 3.0_dp, 2.0_dp, 1.0_dp]
   count_scale = 5.0_dp
   count_response(7) = ieee_value(0.0_dp, ieee_quiet_nan)
   count_normals = 0.0_dp
   count_uniforms = 0.5_dp
   count_gamma = 7.0_dp
   count_predictors(:, 1) = 1.0_dp
   count_predictors(:, 2) = [(sin(0.3_dp*real(i, dp)), i=1, 12)]
   coefficient_normals = 0.0_dp
   coefficient_uniforms = 0.5_dp
   birth_normals = 0.0_dp
   birth_uniforms = 0.5_dp
   structural_normals = 0.0_dp
   structural_uniforms = 0.5_dp
   structural_gamma = 7.0_dp
   structural_forecast_normals = 0.0_dp
   count_fit = bsts_logit_local_level_draws(count_response, count_scale, &
      0.0_dp, 10.0_dp, 0.05_dp, 2.0_dp, 0.05_dp, 0.3_dp, 2, &
      count_normals, count_uniforms, count_gamma)
   call check(count_fit%info == 0 .and. count_fit%family == 1 .and. &
      all(shape(count_fit%state) == [12, 6]) .and. &
      maxval(abs(count_fit%state)) == 0.0_dp .and. &
      maxval(abs(count_fit%fitted_mean - 2.5_dp)) < 1.0e-14_dp .and. &
      all(count_fit%state_variance > 0.0_dp), &
      'binomial-logit local-level supplied draws')
   count_forecast_normals = 0.0_dp
   count_forecast_uniforms = 0.5_dp
   prediction = bsts_logit_predict_draws(count_fit, [5.0_dp, 10.0_dp, &
      20.0_dp], count_forecast_normals, count_forecast_uniforms)
   call check(prediction%info == 0 .and. &
      all(shape(prediction%draws) == [3, 4]) .and. &
      all(prediction%draws >= 0.0_dp) .and. &
      all(prediction%draws <= spread([5.0_dp, 10.0_dp, 20.0_dp], 2, 4)), &
      'binomial-logit supplied forecast')
   count_fit = bsts_logit_trend_seasonal_draws(count_response, count_scale, &
      4, 2, [0.0_dp, 0.0_dp], identity2(), 1.0_dp, &
      [0.05_dp, 0.01_dp, 0.01_dp], [2.0_dp, 2.0_dp, 2.0_dp], &
      [0.05_dp, 0.01_dp, 0.01_dp], 2, structural_normals, &
      structural_uniforms, structural_gamma)
   call check(count_fit%info == 0 .and. count_fit%nseasons == 4 .and. &
      count_fit%season_duration == 2 .and. &
      all(shape(count_fit%structural_state) == [5, 12, 6]) .and. &
      maxval(abs(count_fit%state - count_fit%structural_state(1, :, :) - &
      count_fit%structural_state(3, :, :))) < 1.0e-14_dp .and. &
      all(count_fit%component_variance > 0.0_dp), &
      'binomial-logit trend-seasonal supplied draws')
   prediction = bsts_logit_trend_seasonal_predict_draws(count_fit, &
      [5.0_dp, 10.0_dp, 20.0_dp], structural_forecast_normals, &
      count_forecast_uniforms)
   call check(prediction%info == 0 .and. &
      all(shape(prediction%draws) == [3, 4]), &
      'binomial-logit trend-seasonal forecast')
   count_fit = bsts_logit_regression_draws(count_response, count_scale, &
      count_predictors, [0.2_dp, -0.1_dp], [1.0_dp, 1.0_dp], &
      [1.0_dp, 0.0_dp], [.false., .true.], 1, -1, 0.0_dp, 10.0_dp, &
      0.05_dp, 2.0_dp, 0.05_dp, 0.3_dp, 0.2_dp, 2, count_normals, &
      count_uniforms, coefficient_normals, coefficient_uniforms, &
      birth_normals, birth_uniforms, count_gamma)
   call check(count_fit%info == 0 .and. &
      all(count_fit%included(1, :)) .and. &
      all(.not. count_fit%included(2, :)) .and. &
      all(count_fit%coefficients(1, :) == 0.2_dp) .and. &
      all(count_fit%coefficients(2, :) == 0.0_dp) .and. &
      count_fit%inclusion_probability(1) == 1.0_dp .and. &
      count_fit%inclusion_probability(2) == 0.0_dp .and. &
      maxval(abs(count_fit%regression_contribution - 0.2_dp)) < 1.0e-14_dp, &
      'binomial-logit forced sparse regression')
   prediction = bsts_logit_regression_predict_draws(count_fit, &
      [5.0_dp, 10.0_dp, 20.0_dp], reshape([1.0_dp, 0.0_dp, 1.0_dp, &
      0.5_dp, 1.0_dp, 1.0_dp], [3, 2]), count_forecast_normals, &
      count_forecast_uniforms)
   call check(prediction%info == 0 .and. &
      all(shape(prediction%draws) == [3, 4]), &
      'binomial-logit sparse regression forecast')
   count_response = [1.0_dp, 2.0_dp, 3.0_dp, 2.0_dp, 1.0_dp, 4.0_dp, &
      2.0_dp, 3.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 2.0_dp]
   count_response(5) = ieee_value(0.0_dp, ieee_quiet_nan)
   count_scale = [(0.5_dp + 0.1_dp*real(i, dp), i=1, 12)]
   count_fit = bsts_poisson_local_level_draws(count_response, count_scale, &
      log(2.0_dp), 10.0_dp, 0.05_dp, 2.0_dp, 0.05_dp, 0.2_dp, 2, &
      count_normals, count_uniforms, count_gamma)
   call check(count_fit%info == 0 .and. count_fit%family == 2 .and. &
      maxval(abs(count_fit%state - log(2.0_dp))) < 1.0e-14_dp .and. &
      maxval(abs(count_fit%fitted_mean(:, 1) - 2.0_dp*count_scale)) < &
      1.0e-13_dp, 'Poisson exposure-offset supplied draws')
   count_fit = bsts_poisson_trend_seasonal_draws(count_response, count_scale, &
      4, 2, [log(2.0_dp), 0.0_dp], identity2(), 1.0_dp, &
      [0.05_dp, 0.01_dp, 0.01_dp], [2.0_dp, 2.0_dp, 2.0_dp], &
      [0.05_dp, 0.01_dp, 0.01_dp], 2, structural_normals, &
      structural_uniforms, structural_gamma)
   call check(count_fit%info == 0 .and. count_fit%family == 2 .and. &
      all(shape(count_fit%component_variance) == [3, 6]), &
      'Poisson exposure trend-seasonal supplied draws')
   prediction = bsts_poisson_predict_draws(count_fit, [0.5_dp, 1.0_dp, &
      2.0_dp], count_forecast_normals, count_forecast_uniforms)
   call check(prediction%info == 0 .and. &
      all(prediction%draws >= 0.0_dp) .and. &
      all(prediction%draws == real(nint(prediction%draws), dp)), &
      'Poisson exposure supplied forecast')
   count_fit = bsts_logit_local_level([1.0_dp, 2.0_dp, 3.0_dp, 2.0_dp], &
      [4.0_dp, 4.0_dp, 4.0_dp, 4.0_dp], 4, burn=1)
   call check(count_fit%info == 0, 'random binomial-logit wrapper')
   count_fit = bsts_logit_trend_seasonal([1.0_dp, 2.0_dp, 3.0_dp, 2.0_dp], &
      [4.0_dp, 4.0_dp, 4.0_dp, 4.0_dp], 4, 4, burn=1)
   call check(count_fit%info == 0, 'random binomial trend-seasonal wrapper')
   count_fit = bsts_poisson_local_level([1.0_dp, 2.0_dp, 1.0_dp, 3.0_dp], &
      [0.5_dp, 1.0_dp, 1.5_dp, 2.0_dp], 4, burn=1)
   call check(count_fit%info == 0, 'random Poisson wrapper')
   count_fit = bsts_poisson_trend_seasonal( &
      [1.0_dp, 2.0_dp, 1.0_dp, 3.0_dp], &
      [0.5_dp, 1.0_dp, 1.5_dp, 2.0_dp], 4, 4, burn=1)
   call check(count_fit%info == 0, 'random Poisson trend-seasonal wrapper')
   prediction = bsts_poisson_predict(count_fit, [0.5_dp, 1.0_dp])
   call check(prediction%info == 0 .and. &
      all(shape(prediction%draws) == [2, 3]), 'random Poisson forecast wrapper')
   count_fit = bsts_poisson_regression([1.0_dp, 2.0_dp, 1.0_dp, 3.0_dp], &
      [0.5_dp, 1.0_dp, 1.5_dp, 2.0_dp], reshape([1.0_dp, 0.0_dp, &
      1.0_dp, 0.5_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.5_dp], [4, 2]), 4, &
      burn=1, prior_inclusion_probability=[1.0_dp, 0.0_dp])
   call check(count_fit%info == 0 .and. all(count_fit%included(1, :)) .and. &
      all(.not. count_fit%included(2, :)), 'random Poisson sparse regression')
   prediction = bsts_poisson_regression_predict_draws(count_fit, &
      [0.5_dp, 1.0_dp], reshape([1.0_dp, 0.0_dp, 1.0_dp, 1.0_dp], &
      [2, 2]), count_forecast_normals(:2, :3), &
      count_forecast_uniforms(:2, :3))
   call check(prediction%info == 0 .and. &
      all(shape(prediction%draws) == [2, 3]), &
      'Poisson sparse regression forecast')

   allocate(dirm_y(12), dirm_x(12, 2), dirm_time(12), &
      dirm_state_normals(1, 4, 6), dirm_coefficient_normals(2, 6), &
      dirm_uniforms(2, 6), dirm_gamma(2, 6), &
      dirm_forecast_state(2, 4), dirm_forecast_observation(4, 4))
   dirm_time = [1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4]
   dirm_x(:, 1) = [(0.2_dp*real(i, dp), i=1, 12)]
   dirm_x(:, 2) = [(cos(0.4_dp*real(i, dp)), i=1, 12)]
   dirm_y = 1.5_dp + 0.7_dp*dirm_x(:, 1) + &
      [(0.1_dp*real(dirm_time(i), dp), i=1, 12)]
   dirm_y(7:9) = ieee_value(0.0_dp, ieee_quiet_nan)
   dirm_state_normals = 0.0_dp
   dirm_coefficient_normals = 0.0_dp
   dirm_uniforms = 0.5_dp
   do j = 1, 6
      dirm_gamma(:, j) = [8.0_dp, 4.0_dp]
   end do
   dirm_fit = bsts_dirm_draws(dirm_y, dirm_x, dirm_time, 1.0_dp, &
      2.0_dp, 0.3_dp, 0.02_dp, [0.0_dp, 0.0_dp], identity2(), &
      [1.0_dp, 0.0_dp], 2.0_dp, 0.2_dp, 2.0_dp, 0.02_dp, &
      [.true., .false.], 1, -1, 2, dirm_state_normals, &
      dirm_coefficient_normals, dirm_uniforms, dirm_gamma)
   call check(dirm_fit%info == 0 .and. dirm_fit%time_points == 4 .and. &
      all(shape(dirm_fit%state) == [4, 6]) .and. &
      all(dirm_fit%included(1, :)) .and. &
      all(.not. dirm_fit%included(2, :)) .and. &
      all(ieee_is_finite(dirm_fit%state)) .and. &
      maxval(abs(dirm_fit%contribution(1, :) - &
      dirm_fit%contribution(2, :))) < 1.0e-14_dp .and. &
      maxval(abs(dirm_fit%contribution(2, :) - &
      dirm_fit%contribution(3, :))) < 1.0e-14_dp .and. &
      all(dirm_fit%observation_variance > 0.0_dp) .and. &
      all(dirm_fit%level_variance > 0.0_dp), &
      'dynamic-intercept supplied draws')
   dirm_forecast_state = 0.0_dp
   dirm_forecast_observation = 0.0_dp
   prediction = bsts_dirm_predict_draws(dirm_fit, &
      reshape([1.0_dp, 1.0_dp, 2.0_dp, 2.0_dp, &
      0.0_dp, 0.0_dp, 0.5_dp, 0.5_dp], [4, 2]), [1, 1, 2, 2], &
      dirm_forecast_state, dirm_forecast_observation)
   call check(prediction%info == 0 .and. &
      all(shape(prediction%draws) == [4, 4]) .and. &
      maxval(abs(prediction%draws(1, :) - prediction%draws(2, :))) < &
      1.0e-14_dp .and. maxval(abs(prediction%draws(3, :) - &
      prediction%draws(4, :))) < 1.0e-14_dp, &
      'dynamic-intercept grouped forecast')
   dirm_fit = bsts_dirm(dirm_y, dirm_x, dirm_time, 5, burn=1, &
      prior_inclusion_probability=[1.0_dp, 0.0_dp])
   call check(dirm_fit%info == 0 .and. all(dirm_fit%included(1, :)) .and. &
      all(.not. dirm_fit%included(2, :)), &
      'random dynamic-intercept wrapper')
   prediction = bsts_dirm_predict(dirm_fit, &
      reshape([1.0_dp, 1.0_dp, 2.0_dp, 2.0_dp, &
      0.0_dp, 0.0_dp, 0.5_dp, 0.5_dp], [4, 2]), [1, 1, 2, 2])
   call check(prediction%info == 0 .and. &
      all(shape(prediction%draws) == [4, 4]), &
      'random dynamic-intercept forecast wrapper')

   allocate(mixed_coarse(3), mixed_x(12, 2), mixed_coarse_index(12), &
      mixed_fraction(12), mixed_ends(12), mixed_latent_normals(12, 6), &
      mixed_state_normals(1, 12, 6), mixed_coefficient_normals(2, 6), &
      mixed_uniforms(2, 6), mixed_gamma(2, 6), &
      mixed_forecast_state(8, 4), mixed_forecast_observation(8, 4), &
      mixed_structural_normals(5, 12, 6), mixed_component_gamma(4, 6), &
      mixed_structural_forecast(3, 8, 4))
   mixed_coarse = [10.0_dp, ieee_value(0.0_dp, ieee_quiet_nan), 30.0_dp]
   mixed_coarse_index = [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3]
   mixed_fraction = 1.0_dp
   mixed_fraction(4) = 0.5_dp
   mixed_ends = .false.
   mixed_ends([4, 8, 12]) = .true.
   mixed_x(:, 1) = [(sin(0.2_dp*real(i, dp)), i=1, 12)]
   mixed_x(:, 2) = [(cos(0.3_dp*real(i, dp)), i=1, 12)]
   mixed_latent_normals = 0.0_dp
   mixed_state_normals = 0.0_dp
   mixed_coefficient_normals = 0.0_dp
   mixed_uniforms = 0.5_dp
   do j = 1, 6
      mixed_gamma(:, j) = [8.0_dp, 6.0_dp]
   end do
   mixed_fit = bsts_mixed_draws(mixed_coarse, mixed_x, mixed_coarse_index, &
      mixed_fraction, mixed_ends, 2.0_dp, 4.0_dp, 0.5_dp, 0.02_dp, &
      [0.0_dp, 0.0_dp], identity2(), [1.0_dp, 0.0_dp], &
      2.0_dp, 0.2_dp, 2.0_dp, 0.02_dp, [.true., .false.], 1, -1, 2, &
      mixed_latent_normals, mixed_state_normals, mixed_coefficient_normals, &
      mixed_uniforms, mixed_gamma)
   call check(mixed_fit%info == 0 .and. &
      all(shape(mixed_fit%latent_fine) == [12, 6]) .and. &
      all(mixed_fit%included(1, :)) .and. &
      all(.not. mixed_fit%included(2, :)) .and. &
      maxval(abs(mixed_fit%coarse_fitted(1, :) - 10.0_dp)) < 1.0e-12_dp .and. &
      maxval(abs(mixed_fit%coarse_fitted(3, :) - 30.0_dp)) < 1.0e-12_dp .and. &
      all(ieee_is_finite(mixed_fit%coarse_fitted(2, :))) .and. &
      maxval(abs(mixed_fit%cumulator(1, :))) < 1.0e-14_dp .and. &
      maxval(abs(mixed_fit%cumulator(5, :) - &
      0.5_dp*mixed_fit%latent_fine(4, :))) < 1.0e-12_dp .and. &
      all(mixed_fit%observation_variance > 0.0_dp) .and. &
      all(mixed_fit%level_variance > 0.0_dp), &
      'mixed-frequency supplied draws and constraints')
   mixed_forecast_state = 0.0_dp
   mixed_forecast_observation = 0.0_dp
   mixed_structural_normals = 0.0_dp
   mixed_structural_forecast = 0.0_dp
   do j = 1, 6
      mixed_component_gamma(:, j) = [8.0_dp, 6.0_dp, 6.0_dp, 4.0_dp]
   end do
   mixed_prediction = bsts_mixed_predict_draws(mixed_fit, mixed_x(:8, :), &
      [1, 1, 1, 1, 2, 2, 2, 2], &
      [1.0_dp, 1.0_dp, 1.0_dp, 0.5_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], &
      [.false., .false., .false., .true., .false., .false., .false., .true.], &
      mixed_forecast_state, mixed_forecast_observation)
   call check(mixed_prediction%info == 0 .and. &
      all(shape(mixed_prediction%fine%draws) == [8, 4]) .and. &
      all(shape(mixed_prediction%coarse%draws) == [2, 4]) .and. &
      maxval(abs(mixed_prediction%coarse%draws(1, :) - &
      sum(mixed_prediction%fine%draws(1:3, :), dim=1) - &
      0.5_dp*mixed_prediction%fine%draws(4, :))) < 1.0e-12_dp .and. &
      maxval(abs(mixed_prediction%coarse%draws(2, :) - &
      0.5_dp*mixed_prediction%fine%draws(4, :) - &
      sum(mixed_prediction%fine%draws(5:8, :), dim=1))) < 1.0e-12_dp, &
      'mixed-frequency split-boundary forecast')
   mixed_fit = bsts_mixed_trend_seasonal_draws(mixed_coarse, mixed_x, &
      mixed_coarse_index, mixed_fraction, mixed_ends, 4, 2, &
      [2.0_dp, 0.0_dp], identity2(), 1.0_dp, 0.5_dp, &
      [0.02_dp, 0.01_dp, 0.01_dp], [0.0_dp, 0.0_dp], identity2(), &
      [1.0_dp, 0.0_dp], 2.0_dp, 0.2_dp, &
      [2.0_dp, 2.0_dp, 2.0_dp], [0.02_dp, 0.01_dp, 0.01_dp], &
      [.true., .false.], 1, -1, 2, mixed_latent_normals, &
      mixed_structural_normals, mixed_coefficient_normals, mixed_uniforms, &
      mixed_component_gamma)
   call check(mixed_fit%info == 0 .and. mixed_fit%nseasons == 4 .and. &
      mixed_fit%season_duration == 2 .and. &
      all(shape(mixed_fit%structural_state) == [5, 12, 6]) .and. &
      maxval(abs(mixed_fit%state - mixed_fit%structural_state(1, :, :) - &
      mixed_fit%structural_state(3, :, :))) < 1.0e-13_dp .and. &
      maxval(abs(mixed_fit%coarse_fitted(1, :) - 10.0_dp)) < 1.0e-12_dp .and. &
      maxval(abs(mixed_fit%coarse_fitted(3, :) - 30.0_dp)) < 1.0e-12_dp .and. &
      all(mixed_fit%component_variance > 0.0_dp), &
      'mixed-frequency trend-seasonal supplied draws')
   mixed_prediction = bsts_mixed_trend_seasonal_predict_draws(mixed_fit, &
      mixed_x(:8, :), [1, 1, 1, 1, 2, 2, 2, 2], &
      [1.0_dp, 1.0_dp, 1.0_dp, 0.5_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], &
      [.false., .false., .false., .true., .false., .false., .false., .true.], &
      mixed_structural_forecast, mixed_forecast_observation)
   call check(mixed_prediction%info == 0 .and. &
      all(shape(mixed_prediction%fine%draws) == [8, 4]) .and. &
      all(shape(mixed_prediction%coarse%draws) == [2, 4]), &
      'mixed-frequency trend-seasonal forecast')
   mixed_fit = bsts_mixed(mixed_coarse, mixed_x, mixed_coarse_index, &
      mixed_fraction, mixed_ends, 5, burn=1, &
      prior_inclusion_probability=[1.0_dp, 0.0_dp])
   call check(mixed_fit%info == 0 .and. all(mixed_fit%included(1, :)) .and. &
      all(.not. mixed_fit%included(2, :)), 'random mixed-frequency wrapper')
   mixed_prediction = bsts_mixed_predict(mixed_fit, mixed_x(:4, :), &
      [1, 1, 1, 1], [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], &
      [.false., .false., .false., .true.])
   call check(mixed_prediction%info == 0 .and. &
      all(shape(mixed_prediction%fine%draws) == [4, 4]), &
      'random mixed-frequency forecast wrapper')
   mixed_fit = bsts_mixed_trend_seasonal(mixed_coarse, mixed_x, &
      mixed_coarse_index, mixed_fraction, mixed_ends, 4, 2, 5, burn=1, &
      prior_inclusion_probability=[1.0_dp, 0.0_dp])
   call check(mixed_fit%info == 0 .and. &
      all(shape(mixed_fit%structural_state) == [5, 12, 5]), &
      'random mixed-frequency trend-seasonal wrapper')
   mixed_prediction = bsts_mixed_trend_seasonal_predict(mixed_fit, &
      mixed_x(:4, :), [1, 1, 1, 1], &
      [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], &
      [.false., .false., .false., .true.])
   call check(mixed_prediction%info == 0 .and. &
      all(shape(mixed_prediction%fine%draws) == [4, 4]), &
      'random mixed-frequency trend-seasonal forecast wrapper')

   allocate(y(24), level_normals(1, 24, 16), level_gamma(2, 16))
   do i = 1, size(y)
      y(i) = 1.0_dp + 0.04_dp*real(i, dp) + &
         0.08_dp*sin(0.7_dp*real(i, dp))
   end do
   do j = 1, size(level_normals, 3)
      do i = 1, size(y)
         level_normals(1, i, j) = 0.2_dp*sin(real(3*i + j, dp))
      end do
      level_gamma(:, j) = [14.0_dp + 0.1_dp*real(j, dp), &
         12.0_dp + 0.1_dp*real(j, dp)]
   end do
   level_fit = bsts_local_level_draws(y, y(1), 1.0_dp, 0.2_dp, &
      0.02_dp, 2.0_dp, 0.2_dp, 2.0_dp, 0.02_dp, 4, &
      level_normals, level_gamma)
   call check(level_fit%info == 0 .and. &
      all(shape(level_fit%state) == [1, 24, 16]), 'local-level Gibbs status')
   call check(all(level_fit%observation_variance > 0.0_dp) .and. &
      all(level_fit%state_variance > 0.0_dp) .and. &
      all(ieee_is_finite(level_fit%state)), 'local-level Gibbs draws')
   prediction_errors = bsts_structural_prediction_errors(level_fit, y, &
      [y(1)], reshape([1.0_dp], [1, 1]))
   scaled_errors = bsts_structural_prediction_errors(level_fit, y, &
      [y(1)], reshape([1.0_dp], [1, 1]), .true.)
   call check(prediction_errors%info == 0 .and. scaled_errors%info == 0 .and. &
      all(shape(prediction_errors%draws) == [24, 12]) .and. &
      maxval(abs(scaled_errors%draws*sqrt( &
      prediction_errors%forecast_variance) - prediction_errors%draws)) < &
      1.0e-12_dp .and. prediction_errors%rmse >= 0.0_dp .and. &
      prediction_errors%mae >= 0.0_dp, &
      'structural one-step prediction errors')
   holdout_errors = bsts_local_level_holdout_errors(y, 12, 6, burn=2, &
      standardize=.true.)
   call check(holdout_errors%info == 0 .and. &
      holdout_errors%training_size == 12 .and. &
      holdout_errors%standardized .and. &
      all(shape(holdout_errors%draws) == [24, 4]), &
      'local-level holdout prediction errors')
   allocate(comparison_inputs(2))
   comparison_inputs(1) = prediction_errors
   comparison_inputs(2) = prediction_errors
   comparison_inputs(2)%draws = 2.0_dp*comparison_inputs(2)%draws
   comparison_inputs(2)%mean = 2.0_dp*comparison_inputs(2)%mean
   model_comparison = bsts_compare_prediction_errors(comparison_inputs, 5)
   call check(model_comparison%info == 0 .and. &
      model_comparison%start_index == 5 .and. &
      all(shape(model_comparison%cumulative_absolute_error) == [20, 2]) .and. &
      model_comparison%best_rmse == 1 .and. model_comparison%best_mae == 1 .and. &
      all(model_comparison%rmse_rank == [1, 2]) .and. &
      all(model_comparison%mae_rank == [1, 2]) .and. &
      abs(model_comparison%final_absolute_error(2) - &
      2.0_dp*model_comparison%final_absolute_error(1)) < 1.0e-12_dp, &
      'BSTS model prediction-error comparison')
   comparison_inputs(2)%standardized = .true.
   model_comparison = bsts_compare_prediction_errors(comparison_inputs)
   call check(model_comparison%info /= 0, &
      'BSTS comparison rejects inconsistent scaling')
   deallocate(comparison_inputs)

   allocate(trend_normals(2, 24, 12), trend_gamma(3, 12))
   do j = 1, size(trend_normals, 3)
      do i = 1, size(y)
         trend_normals(1, i, j) = 0.15_dp*sin(real(i + 2*j, dp))
         trend_normals(2, i, j) = 0.15_dp*cos(real(2*i + j, dp))
      end do
      trend_gamma(:, j) = [14.0_dp, 12.0_dp, 12.0_dp] + &
         0.1_dp*real(j, dp)
   end do
   trend_fit = bsts_local_linear_trend_draws(y, [y(1), 0.04_dp], &
      reshape([1.0_dp, 0.0_dp, 0.0_dp, 0.2_dp], [2, 2]), 0.2_dp, &
      [0.02_dp, 0.002_dp], 2.0_dp, 0.2_dp, [2.0_dp, 2.0_dp], &
      [0.02_dp, 0.002_dp], 3, trend_normals, trend_gamma)
   call check(trend_fit%info == 0 .and. &
      all(shape(trend_fit%state) == [2, 24, 12]), &
      'local-linear-trend Gibbs status')
   call check(all(trend_fit%observation_variance > 0.0_dp) .and. &
      all(trend_fit%state_variance > 0.0_dp), &
      'local-linear-trend variance draws')
   prediction_errors = bsts_structural_prediction_errors(trend_fit, y, &
      [y(1), 0.04_dp], &
      reshape([1.0_dp, 0.0_dp, 0.0_dp, 0.2_dp], [2, 2]), .true.)
   call check(prediction_errors%info == 0 .and. &
      all(shape(prediction_errors%draws) == [24, 9]) .and. &
      prediction_errors%standardized, &
      'trend standardized prediction errors')

   allocate(student_y(36), student_state_normals(2, 36, 12), &
      student_weight_normals(2, 35, 8, 12), &
      student_weight_uniforms(2, 35, 8, 12), &
      student_variance_gamma(3, 12), student_degrees_normals(2, 12), &
      student_degrees_uniforms(2, 12))
   do i = 1, 36
      student_y(i) = 1.0_dp + 0.035_dp*real(i, dp) + &
         0.02_dp*sin(0.6_dp*real(i, dp))
   end do
   student_y(18) = student_y(18) + 0.8_dp
   do j = 1, 12
      do i = 1, 36
         student_state_normals(:, i, j) = 0.06_dp*[ &
            sin(real(i + j, dp)), cos(real(2*i + j, dp))]
      end do
      student_weight_normals(:, :, :, j) = 0.05_dp*sin(real(j, dp))
      student_weight_uniforms(:, :, :, j) = 0.5_dp
      student_variance_gamma(:, j) = [20.0_dp, 19.5_dp, 19.5_dp] + &
         0.02_dp*real(j, dp)
      student_degrees_normals(:, j) = 0.1_dp*[sin(real(j, dp)), &
         cos(real(j, dp))]
      student_degrees_uniforms(:, j) = 0.5_dp
   end do
   student_fit = bsts_student_local_linear_trend_draws(student_y, &
      [student_y(1), 0.03_dp], &
      reshape([1.0_dp, 0.0_dp, 0.0_dp, 0.2_dp], [2, 2]), &
      0.05_dp, [0.01_dp, 0.002_dp], [7.0_dp, 9.0_dp], 2.0_dp, &
      0.02_dp, [2.0_dp, 2.0_dp], [0.01_dp, 0.002_dp], &
      [3.0_dp, 3.0_dp], [30.0_dp, 30.0_dp], [0.5_dp, 0.5_dp], &
      3, .true., student_state_normals, student_weight_normals, &
      student_weight_uniforms, student_variance_gamma, &
      student_degrees_normals, student_degrees_uniforms)
   call check(student_fit%info == 0 .and. &
      all(shape(student_fit%degrees_of_freedom) == [2, 12]) .and. &
      all(shape(student_fit%state_weights) == [2, 35, 12]) .and. &
      all(student_fit%state_weights > 0.0_dp) .and. &
      all(student_fit%degrees_of_freedom >= 3.0_dp) .and. &
      all(student_fit%degrees_of_freedom <= 30.0_dp), &
      'Student local-linear-trend posterior')
   allocate(student_forecast_state(2, 2, 9), &
      student_forecast_observation(2, 9))
   student_forecast_state = 0.0_dp
   student_forecast_observation = 0.0_dp
   prediction = bsts_student_trend_predict_draws(student_fit, 2, &
      student_forecast_state, student_forecast_observation)
   call check(prediction%info == 0 .and. &
      all(shape(prediction%draws) == [2, 9]), &
      'Student trend supplied-draw forecast')
   do j = 1, 9
      call check(abs(prediction%draws(1, j) - &
         sum(student_fit%state(:, 36, j + 3))) < 1.0e-12_dp, &
         'Student trend zero-shock forecast')
   end do
   student_fit = bsts_student_local_linear_trend(student_y, 7, burn=2, &
      save_weights=.false., proposal_attempts=16)
   call check(student_fit%info == 0 .and. &
      .not. allocated(student_fit%state_weights), &
      'random Student trend wrapper')
   prediction = bsts_student_trend_predict(student_fit, 2)
   call check(prediction%info == 0 .and. &
      all(shape(prediction%draws) == [2, 5]), &
      'random Student trend forecast wrapper')

   allocate(semilocal_y(48), semilocal_normals(2, 48, 14), &
      semilocal_gamma(3, 14), parameter_normals(2, 14), ar_uniforms(14))
   semilocal_y(1) = 2.0_dp
   slope = 0.0_dp
   do i = 2, size(semilocal_y)
      slope = 0.04_dp + 0.65_dp*(slope - 0.04_dp) + &
         0.01_dp*sin(0.8_dp*real(i, dp))
      semilocal_y(i) = semilocal_y(i - 1) + slope + &
         0.015_dp*cos(0.5_dp*real(i, dp))
   end do
   do j = 1, size(semilocal_normals, 3)
      do i = 1, size(semilocal_y)
         semilocal_normals(:, i, j) = 0.12_dp*[ &
            sin(real(i + j, dp)), cos(real(2*i + j, dp))]
      end do
      semilocal_gamma(:, j) = [26.0_dp, 24.0_dp, 24.0_dp] + &
         0.1_dp*real(j, dp)
      parameter_normals(:, j) = 0.1_dp*[sin(real(j, dp)), &
         cos(real(j, dp))]
      ar_uniforms(j) = 0.2_dp + 0.6_dp*real(j - 1, dp)/13.0_dp
   end do
   semilocal_fit = bsts_semilocal_trend_draws(semilocal_y, &
      [2.0_dp, 0.0_dp], reshape([1.0_dp, 0.0_dp, 0.0_dp, 0.2_dp], &
      [2, 2]), 0.1_dp, [0.01_dp, 0.002_dp], 0.02_dp, 0.6_dp, &
      2.0_dp, 0.1_dp, [2.0_dp, 2.0_dp], [0.01_dp, 0.002_dp], &
      0.0_dp, 0.2_dp, 0.7_dp, 0.1_dp, .true., .true., 3, &
      semilocal_normals, semilocal_gamma, parameter_normals, ar_uniforms)
   call check(semilocal_fit%info == 0 .and. semilocal_fit%is_semilocal .and. &
      all(shape(semilocal_fit%state) == [2, 48, 14]), &
      'semilocal trend Gibbs status')
   call check(all(semilocal_fit%slope_ar > 0.0_dp) .and. &
      all(semilocal_fit%slope_ar < 1.0_dp) .and. &
      all(ieee_is_finite(semilocal_fit%slope_mean)) .and. &
      all(semilocal_fit%component_variance > 0.0_dp), &
      'semilocal constrained parameter draws')
   prediction = bsts_predict(semilocal_fit, 5)
   call check(prediction%info == 0 .and. &
      all(shape(prediction%draws) == [5, 11]) .and. &
      all(ieee_is_finite(prediction%draws)), 'semilocal posterior forecast')
   allocate(semi_forecast_state(2, 1, 11), semi_forecast_observation(1, 11))
   semi_forecast_state = 0.0_dp
   semi_forecast_observation = 0.0_dp
   prediction = bsts_predict_draws(semilocal_fit, 1, semi_forecast_state, &
      semi_forecast_observation)
   do j = 1, 11
      call check(abs(prediction%draws(1, j) - &
         sum(semilocal_fit%state(:, 48, j + 3))) < 1.0e-12_dp, &
         'semilocal affine level forecast')
   end do

   allocate(seasonal_y(32), seasonal_normals(3, 32, 12), &
      seasonal_gamma(2, 12))
   seasonal_pattern = [1.0_dp, 1.0_dp, -1.0_dp, -1.0_dp, &
      0.5_dp, 0.5_dp, -0.5_dp, -0.5_dp]
   do i = 1, size(seasonal_y)
      seasonal_y(i) = seasonal_pattern(1 + modulo(i - 1, 8)) + &
         0.03_dp*sin(real(i, dp))
   end do
   do j = 1, size(seasonal_normals, 3)
      do i = 1, size(seasonal_y)
         seasonal_normals(:, i, j) = 0.1_dp*[ &
            sin(real(i + j, dp)), cos(real(2*i + j, dp)), &
            sin(real(3*i + 2*j, dp))]
      end do
      seasonal_gamma(:, j) = [18.0_dp, 10.0_dp] + 0.1_dp*real(j, dp)
   end do
   seasonal_fit = bsts_seasonal_draws(seasonal_y, 4, 2, 1.0_dp, &
      0.2_dp, 0.02_dp, 2.0_dp, 0.2_dp, 2.0_dp, 0.02_dp, 2, &
      seasonal_normals, seasonal_gamma)
   call check(seasonal_fit%info == 0 .and. &
      all(shape(seasonal_fit%state) == [3, 32, 12]) .and. &
      all(shape(seasonal_fit%state_variance) == [1, 12]), &
      'duration-aware seasonal Gibbs status')
   call check(maxval(abs(seasonal_fit%transition_schedule(:, :, 1) - &
      reshape([1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 1.0_dp], [3, 3]))) < 1.0e-14_dp .and. &
      all(seasonal_fit%transition_schedule(1, :, 2) == -1.0_dp) .and. &
      seasonal_fit%state_loading_schedule(1, 1, 1) == 0.0_dp .and. &
      seasonal_fit%state_loading_schedule(1, 1, 2) == 1.0_dp, &
      'seasonal transition and innovation schedule')

   allocate(trig_y(30), trig_normals(4, 30, 10), trig_gamma(2, 10))
   do i = 1, size(trig_y)
      trig_y(i) = 1.2_dp*cos(2.0_dp*acos(-1.0_dp)*real(i, dp)/12.0_dp) + &
         0.4_dp*sin(4.0_dp*acos(-1.0_dp)*real(i, dp)/12.0_dp)
   end do
   do j = 1, size(trig_normals, 3)
      do i = 1, size(trig_y)
         trig_normals(:, i, j) = 0.1_dp*[sin(real(i + j, dp)), &
            cos(real(i + 2*j, dp)), sin(real(2*i + j, dp)), &
            cos(real(3*i + j, dp))]
      end do
      trig_gamma(:, j) = [17.0_dp, 60.0_dp] + 0.1_dp*real(j, dp)
   end do
   trig_fit = bsts_trig_draws(trig_y, 12.0_dp, [1.0_dp, 2.0_dp], &
      [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], &
      reshape([1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [4, 4]), &
      0.2_dp, 0.01_dp, 2.0_dp, 0.2_dp, 2.0_dp, 0.01_dp, 2, &
      trig_normals, trig_gamma)
   call check(trig_fit%info == 0 .and. &
      all(shape(trig_fit%state) == [4, 30, 10]) .and. &
      all(shape(trig_fit%component_variance) == [1, 10]), &
      'harmonic trigonometric Gibbs status')
   call check(maxval(abs(matmul(trig_fit%transition, &
      transpose(trig_fit%transition)) - reshape([1.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 1.0_dp], [4, 4]))) < 1.0e-12_dp .and. &
      all(trig_fit%state_variance(1, :) == &
      trig_fit%state_variance(4, :)), 'harmonic rotation and shared variance')

   allocate(monthly_y(75), monthly_normals(11, 75, 10), &
      monthly_gamma(2, 10))
   do i = 1, 75
      if (i <= 17) then
         monthly_y(i) = 1.0_dp
      else if (i <= 46) then
         monthly_y(i) = -0.5_dp
      else
         monthly_y(i) = 0.2_dp
      end if
      monthly_y(i) = monthly_y(i) + 0.01_dp*sin(real(i, dp))
   end do
   do j = 1, 10
      do i = 1, 75
         monthly_normals(:, i, j) = 0.03_dp*sin(real(i + j, dp))
      end do
      monthly_gamma(:, j) = [40.0_dp, 3.0_dp] + 0.02_dp*real(j, dp)
   end do
   monthly_fit = bsts_monthly_annual_cycle_draws(monthly_y, &
      date_t(2024, 1, 15), [(0.0_dp, i=1, 11)], &
      reshape([((merge(1.0_dp, 0.0_dp, i == j), i=1, 11), &
      j=1, 11)], [11, 11]), &
      0.1_dp, 0.01_dp, 2.0_dp, 0.1_dp, 2.0_dp, 0.01_dp, 2, &
      monthly_normals, monthly_gamma)
   call check(monthly_fit%info == 0 .and. monthly_fit%is_monthly .and. &
      all(shape(monthly_fit%state) == [11, 75, 10]) .and. &
      monthly_fit%last_date%year == 2024 .and. &
      monthly_fit%last_date%month == 3 .and. monthly_fit%last_date%day == 29, &
      'monthly annual cycle posterior')
   call check(monthly_fit%state_loading_schedule(1, 1, 17) == 1.0_dp .and. &
      monthly_fit%state_loading_schedule(1, 1, 46) == 1.0_dp .and. &
      monthly_fit%state_loading_schedule(1, 1, 16) == 0.0_dp .and. &
      all(monthly_fit%transition_schedule(1, :, 17) == -1.0_dp), &
      'monthly leap-year transition schedule')
   allocate(monthly_forecast(4, 8))
   monthly_forecast = 0.0_dp
   prediction = bsts_monthly_predict_draws(monthly_fit, 4, monthly_forecast, &
      monthly_forecast)
   call check(prediction%info == 0 .and. &
      all(shape(prediction%draws) == [4, 8]), &
      'monthly supplied-draw forecast')
   do j = 1, 8
      call check(abs(prediction%draws(1, j) - &
         monthly_fit%state(1, 75, j + 2)) < 1.0e-12_dp .and. &
         abs(prediction%draws(2, j) - prediction%draws(1, j)) < 1.0e-12_dp .and. &
         abs(prediction%draws(3, j) + &
         sum(monthly_fit%state(:, 75, j + 2))) < 1.0e-12_dp .and. &
         abs(prediction%draws(4, j) - prediction%draws(3, j)) < 1.0e-12_dp, &
         'monthly boundary forecast')
   end do
   monthly_fit = bsts_monthly_annual_cycle(monthly_y, date_t(2024, 1, 15), &
      7, burn=2)
   call check(monthly_fit%info == 0 .and. monthly_fit%is_monthly, &
      'random monthly annual cycle wrapper')
   prediction = bsts_monthly_predict(monthly_fit, 4)
   call check(prediction%info == 0 .and. &
      all(shape(prediction%draws) == [4, 5]), &
      'random monthly forecast wrapper')

   fixed_holiday = bsts_fixed_date_holiday(1, 1, 1, 1, 'New Year')
   call check(fixed_holiday%info == 0 .and. &
      bsts_holiday_position(fixed_holiday, date_t(2023, 12, 31)) == 1 .and. &
      bsts_holiday_position(fixed_holiday, date_t(2024, 1, 1)) == 2 .and. &
      bsts_holiday_position(fixed_holiday, date_t(2024, 1, 3)) == 0, &
      'fixed-date holiday window')
   weekday_holiday = bsts_nth_weekday_holiday(9, 1, 1, 0, 0, 'Labor Day')
   call check(bsts_holiday_position(weekday_holiday, &
      date_t(2024, 9, 2)) == 1, 'nth-weekday holiday')
   weekday_holiday = bsts_last_weekday_holiday(5, 1, 0, 0, 'Memorial Day')
   call check(bsts_holiday_position(weekday_holiday, &
      date_t(2024, 5, 27)) == 1, 'last-weekday holiday')
   weekday_holiday = bsts_named_holiday('EasterSunday', 1, 1)
   call check(bsts_holiday_position(weekday_holiday, &
      date_t(2024, 3, 31)) == 2, 'named holiday')
   holiday = bsts_date_range_holiday([date_t(2024, 1, 10), &
      date_t(2024, 2, 10), date_t(2024, 3, 1)], &
      [date_t(2024, 1, 12), date_t(2024, 2, 12), &
      date_t(2024, 3, 3)], 'Campaign')
   call check(holiday%info == 0 .and. bsts_holiday_width(holiday) == 3 .and. &
      bsts_holiday_position(holiday, date_t(2024, 2, 11)) == 2, &
      'date-range holiday')
   allocate(holiday_y(50), holiday_normals(3, 50, 10), holiday_gamma(2, 10))
   holiday_y = 0.0_dp
   holiday_y(10:12) = [1.0_dp, 2.0_dp, 1.0_dp]
   holiday_y(41:43) = [1.2_dp, 2.1_dp, 0.9_dp]
   do j = 1, 10
      do i = 1, 50
         holiday_normals(:, i, j) = 0.04_dp*[sin(real(i + j, dp)), &
            cos(real(2*i + j, dp)), sin(real(3*i + j, dp))]
      end do
      holiday_gamma(:, j) = [27.0_dp, 5.0_dp] + 0.02_dp*real(j, dp)
   end do
   holiday_fit = bsts_random_walk_holiday_draws(holiday_y, &
      date_t(2024, 1, 1), holiday, [0.0_dp, 0.0_dp, 0.0_dp], &
      reshape([1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 1.0_dp], [3, 3]), 0.1_dp, 0.01_dp, &
      2.0_dp, 0.1_dp, 2.0_dp, 0.01_dp, 2, holiday_normals, holiday_gamma)
   call check(holiday_fit%info == 0 .and. holiday_fit%is_holiday .and. &
      all(shape(holiday_fit%state) == [3, 50, 10]) .and. &
      holiday_fit%observation_schedule(1, 10) == 1.0_dp .and. &
      holiday_fit%observation_schedule(2, 11) == 1.0_dp .and. &
      holiday_fit%state_loading_schedule(1, 1, 9) == 1.0_dp .and. &
      holiday_fit%state_loading_schedule(1, 1, 40) == 1.0_dp, &
      'random-walk holiday posterior')
   allocate(holiday_forecast(15, 8))
   holiday_forecast = 0.0_dp
   prediction = bsts_holiday_predict_draws(holiday_fit, 15, holiday_forecast, &
      holiday_forecast)
   call check(prediction%info == 0 .and. &
      all(shape(prediction%draws) == [15, 8]) .and. &
      all(prediction%draws(:10, :) == 0.0_dp), &
      'holiday supplied-draw forecast')
   do j = 1, 8
      call check(abs(prediction%draws(11, j) - &
         holiday_fit%state(1, 50, j + 2)) < 1.0e-12_dp .and. &
         abs(prediction%draws(12, j) - &
         holiday_fit%state(2, 50, j + 2)) < 1.0e-12_dp .and. &
         abs(prediction%draws(13, j) - &
         holiday_fit%state(3, 50, j + 2)) < 1.0e-12_dp, &
         'holiday future influence window')
   end do
   holiday_fit = bsts_random_walk_holiday(holiday_y, date_t(2024, 1, 1), &
      holiday, 7, burn=2)
   call check(holiday_fit%info == 0 .and. holiday_fit%is_holiday, &
      'random holiday wrapper')
   prediction = bsts_holiday_predict(holiday_fit, 15)
   call check(prediction%info == 0 .and. &
      all(shape(prediction%draws) == [15, 5]), &
      'random holiday forecast wrapper')

   allocate(regression_holidays(2))
   regression_holidays(1) = holiday
   regression_holidays(2) = bsts_fixed_date_holiday(2, 14, 0, 0, &
      'Valentine Day')
   holiday_design = bsts_holiday_design(date_t(2024, 1, 1), 50, &
      regression_holidays)
   call check(all(shape(holiday_design) == [50, 4]) .and. &
      all(holiday_design(10, :) == [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]) .and. &
      all(holiday_design(45, :) == [0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp]), &
      'fixed holiday regression design')
   allocate(holiday_regression_y(50), holiday_regression_normals(4, 10), &
      holiday_regression_gamma(10), holiday_regression_offset(50, 10))
   holiday_regression_offset = 0.7_dp
   holiday_regression_y = 0.7_dp + matmul(holiday_design, &
      [1.0_dp, 2.0_dp, 1.0_dp, -1.5_dp])
   holiday_regression_y(11) = ieee_value(0.0_dp, ieee_quiet_nan)
   holiday_regression_normals = 0.0_dp
   holiday_regression_gamma = 27.0_dp
   holiday_regression_fit = bsts_regression_holiday_draws( &
      holiday_regression_y, date_t(2024, 1, 1), regression_holidays, &
      0.0_dp, 100.0_dp, 0.01_dp, 2.0_dp, 0.01_dp, 2, &
      holiday_regression_normals, holiday_regression_gamma, &
      holiday_regression_offset)
   call check(holiday_regression_fit%info == 0 .and. &
      all(shape(holiday_regression_fit%coefficients) == [4, 10]) .and. &
      all(holiday_regression_fit%coefficient_offset == [0, 3, 4]) .and. &
      holiday_regression_fit%coefficients(1, 1) > 0.9_dp .and. &
      holiday_regression_fit%coefficients(2, 1) > 1.8_dp .and. &
      holiday_regression_fit%coefficients(4, 1) < -1.3_dp, &
      'fixed holiday regression posterior')
   allocate(holiday_regression_forecast(15, 8))
   holiday_regression_forecast = 0.0_dp
   prediction = bsts_regression_holiday_predict_draws( &
      holiday_regression_fit, 15, holiday_regression_forecast)
   call check(prediction%info == 0 .and. &
      all(shape(prediction%draws) == [15, 8]) .and. &
      all(prediction%draws(:10, :) == 0.0_dp), &
      'fixed holiday regression supplied-draw forecast')
   do j = 1, 8
      call check(abs(prediction%draws(11, j) - &
         holiday_regression_fit%coefficients(1, j + 2)) < 1.0e-12_dp .and. &
         abs(prediction%draws(12, j) - &
         holiday_regression_fit%coefficients(2, j + 2)) < 1.0e-12_dp .and. &
         abs(prediction%draws(13, j) - &
         holiday_regression_fit%coefficients(3, j + 2)) < 1.0e-12_dp, &
         'fixed holiday regression future influence window')
   end do
   holiday_regression_fit = bsts_regression_holiday(holiday_regression_y, &
      date_t(2024, 1, 1), regression_holidays, 7, burn=2, &
      offset_draws=holiday_regression_offset(:, :7))
   call check(holiday_regression_fit%info == 0, &
      'random fixed holiday regression wrapper')
   prediction = bsts_regression_holiday_predict(holiday_regression_fit, 15)
   call check(prediction%info == 0 .and. &
      all(shape(prediction%draws) == [15, 5]), &
      'random fixed holiday regression forecast wrapper')

   allocate(hierarchical_holidays(3))
   hierarchical_holidays(1) = holiday
   hierarchical_holidays(2) = bsts_date_range_holiday( &
      [date_t(2024, 1, 20), date_t(2024, 3, 5)], &
      [date_t(2024, 1, 22), date_t(2024, 3, 7)], 'Second campaign')
   hierarchical_holidays(3) = bsts_date_range_holiday( &
      [date_t(2024, 1, 30), date_t(2024, 3, 10)], &
      [date_t(2024, 2, 1), date_t(2024, 3, 12)], 'Third campaign')
   hierarchical_design = bsts_holiday_design(date_t(2024, 1, 1), 50, &
      hierarchical_holidays)
   allocate(hierarchical_y(50), &
      hierarchical_coefficient_normals(9, 10), &
      hierarchical_mean_normals(3, 10), &
      hierarchical_wishart_normals(3, 3, 10), &
      hierarchical_wishart_gamma(3, 10), &
      hierarchical_residual_gamma(10), hierarchical_offset(50, 10))
   hierarchical_offset = 0.3_dp
   hierarchical_y = 0.3_dp + matmul(hierarchical_design, &
      [1.0_dp, 2.0_dp, 1.0_dp, 1.1_dp, 1.8_dp, 0.9_dp, &
      0.8_dp, 2.2_dp, 1.1_dp])
   hierarchical_y(21) = ieee_value(0.0_dp, ieee_quiet_nan)
   hierarchical_coefficient_normals = 0.0_dp
   hierarchical_mean_normals = 0.0_dp
   hierarchical_wishart_normals = 0.0_dp
   hierarchical_wishart_gamma = 5.0_dp
   hierarchical_residual_gamma = 27.0_dp
   hierarchical_holiday_fit = bsts_hierarchical_regression_holiday_draws( &
      hierarchical_y, date_t(2024, 1, 1), hierarchical_holidays, &
      [0.0_dp, 0.0_dp, 0.0_dp], 100.0_dp*identity3(), 5.0_dp, &
      0.2_dp*identity3(), 0.01_dp, 2.0_dp, 0.01_dp, 2, &
      hierarchical_coefficient_normals, hierarchical_mean_normals, &
      hierarchical_wishart_normals, hierarchical_wishart_gamma, &
      hierarchical_residual_gamma, hierarchical_offset)
   call check(hierarchical_holiday_fit%info == 0 .and. &
      hierarchical_holiday_fit%hierarchical .and. &
      all(shape(hierarchical_holiday_fit%coefficients) == [9, 10]) .and. &
      all(shape(hierarchical_holiday_fit%coefficient_mean) == [3, 10]) .and. &
      all(shape(hierarchical_holiday_fit%coefficient_variance) == [3, 3, 10]) .and. &
      all(hierarchical_holiday_fit%coefficient_offset == [0, 3, 6, 9]) .and. &
      all(hierarchical_holiday_fit%coefficients([1, 4, 7], 1) > 0.5_dp) .and. &
      all(hierarchical_holiday_fit%coefficients([2, 8], 1) > 1.0_dp), &
      'hierarchical holiday regression posterior')
   do j = 1, 10
      call check(all([(hierarchical_holiday_fit% &
         coefficient_variance(i, i, j) > 0.0_dp, i=1, 3)]) .and. &
         maxval(abs(hierarchical_holiday_fit%coefficient_variance(:, :, j) - &
         transpose(hierarchical_holiday_fit% &
         coefficient_variance(:, :, j)))) < 1.0e-12_dp, &
         'hierarchical holiday covariance draw')
   end do
   allocate(hierarchical_forecast(25, 8))
   hierarchical_forecast = 0.0_dp
   prediction = bsts_regression_holiday_predict_draws( &
      hierarchical_holiday_fit, 25, hierarchical_forecast)
   call check(prediction%info == 0 .and. &
      all(shape(prediction%draws) == [25, 8]) .and. &
      all(prediction%draws(:10, :) == 0.0_dp), &
      'hierarchical holiday supplied-draw forecast')
   do j = 1, 8
      call check(abs(prediction%draws(11, j) - &
         hierarchical_holiday_fit%coefficients(1, j + 2)) < 1.0e-12_dp .and. &
         abs(prediction%draws(15, j) - &
         hierarchical_holiday_fit%coefficients(4, j + 2)) < 1.0e-12_dp .and. &
         abs(prediction%draws(20, j) - &
         hierarchical_holiday_fit%coefficients(7, j + 2)) < 1.0e-12_dp, &
         'hierarchical holiday future influence windows')
   end do
   hierarchical_holiday_fit = bsts_hierarchical_regression_holiday( &
      hierarchical_y, date_t(2024, 1, 1), hierarchical_holidays, 7, burn=2, &
      offset_draws=hierarchical_offset(:, :7))
   call check(hierarchical_holiday_fit%info == 0 .and. &
      hierarchical_holiday_fit%hierarchical, &
      'random hierarchical holiday regression wrapper')

   allocate(shared_response(30, 3), shared_factors(2, 30), &
      shared_state_normals(2, 30, 8), shared_loading_normals(2, 3, 8), &
      shared_gamma(5, 8), shared_offset(30, 3, 8))
   shared_loadings = reshape([1.0_dp, 0.5_dp, -0.2_dp, &
      0.0_dp, 1.0_dp, 0.8_dp], [3, 2])
   do i = 1, 30
      shared_factors(:, i) = [0.1_dp*real(i, dp), &
         0.4_dp*sin(0.2_dp*real(i, dp))]
      shared_response(i, :) = matmul(shared_loadings, shared_factors(:, i)) + &
         0.2_dp + 0.01_dp*[sin(real(i, dp)), cos(real(i, dp)), &
         sin(2.0_dp*real(i, dp))]
   end do
   shared_response(12, 2) = ieee_value(0.0_dp, ieee_quiet_nan)
   shared_state_normals = 0.0_dp
   shared_loading_normals = 0.0_dp
   shared_gamma = 17.0_dp
   shared_offset = 0.2_dp
   shared_fit = bsts_shared_local_level_draws(shared_response, &
      shared_factors(:, 1), identity2(), shared_loadings, shared_loadings, &
      10.0_dp*spread([1.0_dp, 1.0_dp], 1, 3), &
      [0.01_dp, 0.01_dp, 0.01_dp], [0.02_dp, 0.02_dp], &
      [2.0_dp, 2.0_dp, 2.0_dp], [0.01_dp, 0.01_dp, 0.01_dp], &
      [2.0_dp, 2.0_dp], [0.02_dp, 0.02_dp], 2, shared_state_normals, &
      shared_loading_normals, shared_gamma, shared_offset)
   call check(shared_fit%info == 0 .and. &
      all(shape(shared_fit%state) == [2, 30, 8]) .and. &
      all(shape(shared_fit%loadings) == [3, 2, 8]) .and. &
      all(shape(shared_fit%contribution) == [3, 30, 8]) .and. &
      all(shared_fit%observation_variance > 0.0_dp) .and. &
      all(shared_fit%factor_variance > 0.0_dp), &
      'shared local-level posterior')
   do j = 1, 8
      call check(shared_fit%loadings(1, 1, j) == 1.0_dp .and. &
         shared_fit%loadings(1, 2, j) == 0.0_dp .and. &
         shared_fit%loadings(2, 2, j) == 1.0_dp, &
         'shared local-level loading identification')
   end do
   allocate(shared_forecast_state(2, 4, 6), &
      shared_forecast_observation(3, 4, 6))
   shared_forecast_state = 0.0_dp
   shared_forecast_observation = 0.0_dp
   multivariate_prediction = bsts_shared_local_level_predict_draws( &
      shared_fit, 4, shared_forecast_state, shared_forecast_observation)
   call check(multivariate_prediction%info == 0 .and. &
      all(shape(multivariate_prediction%draws) == [3, 4, 6]), &
      'shared local-level supplied-draw forecast')
   do j = 1, 6
      do i = 1, 4
         call check(maxval(abs(multivariate_prediction%draws(:, i, j) - &
            matmul(shared_fit%loadings(:, :, j + 2), &
            shared_fit%state(:, 30, j + 2)))) < 1.0e-12_dp, &
            'shared local-level constant forecast')
      end do
   end do
   allocate(shared_inclusion_uniforms(2, 3, 8))
   shared_inclusion_prior = reshape([1.0_dp, 0.0_dp, 1.0_dp, &
      0.0_dp, 1.0_dp, 0.5_dp], [3, 2])
   shared_initial_inclusion = .false.
   shared_inclusion_uniforms = tiny(1.0_dp)
   shared_fit = bsts_shared_local_level_draws(shared_response, &
      shared_factors(:, 1), identity2(), shared_loadings, shared_loadings, &
      10.0_dp*spread([1.0_dp, 1.0_dp], 1, 3), &
      [0.01_dp, 0.01_dp, 0.01_dp], [0.02_dp, 0.02_dp], &
      [2.0_dp, 2.0_dp, 2.0_dp], [0.01_dp, 0.01_dp, 0.01_dp], &
      [2.0_dp, 2.0_dp], [0.02_dp, 0.02_dp], 2, shared_state_normals, &
      shared_loading_normals, shared_gamma, shared_offset, &
      prior_inclusion_probability=shared_inclusion_prior, &
      initial_inclusion=shared_initial_inclusion, &
      inclusion_uniform_draws=shared_inclusion_uniforms, maximum_flips=-1)
   call check(shared_fit%info == 0 .and. shared_fit%spike_slab .and. &
      all(shape(shared_fit%included) == [3, 2, 8]) .and. &
      all(.not. shared_fit%included(2, 1, :)) .and. &
      all(shared_fit%included(3, 1, :)) .and. &
      all(shared_fit%included(3, 2, :)) .and. &
      all(shared_fit%loadings(2, 1, :) == 0.0_dp) .and. &
      shared_fit%inclusion_probability(2, 1) == 0.0_dp .and. &
      shared_fit%inclusion_probability(3, 1) == 1.0_dp .and. &
      shared_fit%inclusion_probability(3, 2) == 1.0_dp, &
      'shared local-level loading selection')
   shared_fit = bsts_shared_local_level_draws(shared_response, &
      shared_factors(:, 1), identity2(), shared_loadings, shared_loadings, &
      10.0_dp*spread([1.0_dp, 1.0_dp], 1, 3), &
      [0.01_dp, 0.01_dp, 0.01_dp], [0.02_dp, 0.02_dp], &
      [2.0_dp, 2.0_dp, 2.0_dp], [0.01_dp, 0.01_dp, 0.01_dp], &
      [2.0_dp, 2.0_dp], [0.02_dp, 0.02_dp], 2, shared_state_normals, &
      shared_loading_normals, shared_gamma, shared_offset, &
      prior_inclusion_probability=shared_inclusion_prior, &
      initial_inclusion=shared_initial_inclusion, &
      inclusion_uniform_draws=shared_inclusion_uniforms, maximum_flips=0)
   call check(shared_fit%info == 0 .and. &
      all(.not. shared_fit%included(3, 2, :)) .and. &
      all(shared_fit%loadings(3, 2, :) == 0.0_dp), &
      'shared local-level loading flip limit')
   shared_fit = bsts_shared_local_level(shared_response, 2, 7, burn=2, &
      offset_draws=shared_offset(:, :, :7))
   call check(shared_fit%info == 0 .and. shared_fit%spike_slab .and. &
      allocated(shared_fit%included), 'random shared local-level wrapper')
   multivariate_prediction = bsts_shared_local_level_predict(shared_fit, 3)
   call check(multivariate_prediction%info == 0 .and. &
      all(shape(multivariate_prediction%draws) == [3, 3, 5]), &
      'random shared local-level forecast wrapper')

   allocate(mbsts_response(30, 3), mbsts_predictors(30, 2, 3), &
      mbsts_state_normals(8, 30, 8), mbsts_loading_normals(2, 3, 8), &
      mbsts_gamma(5, 8), mbsts_offset(30, 3, 8), &
      mbsts_regression_normals(2, 3, 8), &
      mbsts_regression_uniforms(2, 3, 8), &
      mbsts_series_normals(3, 30, 8), mbsts_series_gamma(3, 8), &
      mbsts_trend_normals(2, 3, 30, 8), mbsts_trend_gamma(2, 3, 8), &
      mbsts_seasonal_normals(3, 3, 30, 8), &
      mbsts_seasonal_gamma(3, 8))
   mbsts_beta = reshape([0.2_dp, 0.4_dp, -0.1_dp, -0.3_dp, &
      0.3_dp, 0.2_dp], [2, 3])
   do j = 1, 3
      do i = 1, 30
         mbsts_predictors(i, :, j) = [1.0_dp, &
            cos(0.37_dp*real(i + j, dp))]
         mbsts_response(i, j) = dot_product(shared_loadings(j, :), &
            shared_factors(:, i)) + dot_product(mbsts_predictors(i, :, j), &
            mbsts_beta(:, j)) + 0.1_dp + 0.01_dp*sin(real(i + j, dp))
      end do
   end do
   mbsts_response(14, 3) = ieee_value(0.0_dp, ieee_quiet_nan)
   mbsts_state_normals = 0.0_dp
   mbsts_loading_normals = 0.0_dp
   mbsts_gamma = 17.0_dp
   mbsts_offset = 0.1_dp
   mbsts_series_normals = 0.0_dp
   mbsts_series_gamma = 15.0_dp
   mbsts_trend_normals = 0.0_dp
   mbsts_trend_gamma = 15.0_dp
   mbsts_seasonal_normals = 0.0_dp
   mbsts_seasonal_gamma = 15.0_dp
   shared_fit = bsts_mbsts_draws(mbsts_response, mbsts_predictors, &
      reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], &
      [2, 3]), spread(10.0_dp*identity2(), 3, 3), shared_factors(:, 1), &
      identity2(), shared_loadings, shared_loadings, &
      10.0_dp*spread([1.0_dp, 1.0_dp], 1, 3), &
      [0.01_dp, 0.01_dp, 0.01_dp], [0.02_dp, 0.02_dp], &
      [2.0_dp, 2.0_dp, 2.0_dp], [0.01_dp, 0.01_dp, 0.01_dp], &
      [2.0_dp, 2.0_dp], [0.02_dp, 0.02_dp], 2, mbsts_state_normals, &
      mbsts_loading_normals, mbsts_gamma, mbsts_offset)
   call check(shared_fit%info == 0 .and. shared_fit%is_mbsts .and. &
      all(shape(shared_fit%regression_coefficients) == [2, 3, 8]) .and. &
      all(shape(shared_fit%regression_contribution) == [3, 30, 8]) .and. &
      all(shape(shared_fit%fitted) == [3, 30, 8]) .and. &
      all(shape(shared_fit%residuals) == [3, 30, 8]), &
      'multivariate BSTS regression posterior')
   do j = 1, 8
      call check(maxval(abs(shared_fit%residuals(:, :, j) - &
         (transpose(mbsts_response) - shared_fit%fitted(:, :, j))), &
         mask=ieee_is_finite(transpose(mbsts_response))) < 1.0e-12_dp, &
         'multivariate BSTS fitted residual identity')
   end do
   allocate(mbsts_future_predictors(4, 2, 3), &
      mbsts_forecast_state(2, 4, 6), &
      mbsts_forecast_observation(3, 4, 6))
   do j = 1, 3
      do i = 1, 4
         mbsts_future_predictors(i, :, j) = [1.0_dp, &
            cos(0.37_dp*real(30 + i + j, dp))]
      end do
   end do
   mbsts_forecast_state = 0.0_dp
   mbsts_forecast_observation = 0.0_dp
   multivariate_prediction = bsts_mbsts_predict_draws(shared_fit, &
      mbsts_future_predictors, mbsts_forecast_state, &
      mbsts_forecast_observation)
   call check(multivariate_prediction%info == 0 .and. &
      all(shape(multivariate_prediction%draws) == [3, 4, 6]), &
      'multivariate BSTS supplied-draw forecast')
   do j = 1, 6
      do i = 1, 4
         call check(maxval(abs(multivariate_prediction%draws(:, i, j) - &
            matmul(shared_fit%loadings(:, :, j + 2), &
            shared_fit%state(:, 30, j + 2)) - &
            [(dot_product(mbsts_future_predictors(i, :, k), &
            shared_fit%regression_coefficients(:, k, j + 2)), k=1, 3)])) < &
            1.0e-11_dp, &
            'multivariate BSTS regression forecast identity')
      end do
   end do
   allocate(mbsts_series_forecast(3, 4, 6), &
      mbsts_trend_forecast(2, 3, 4, 6), &
      mbsts_seasonal_forecast(3, 3, 4, 6))
   mbsts_series_forecast = 0.0_dp
   mbsts_trend_forecast = 0.0_dp
   mbsts_seasonal_forecast = 0.0_dp
   shared_fit = bsts_mbsts_draws(mbsts_response, mbsts_predictors, &
      reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], &
      [2, 3]), spread(10.0_dp*identity2(), 3, 3), shared_factors(:, 1), &
      identity2(), shared_loadings, shared_loadings, &
      10.0_dp*spread([1.0_dp, 1.0_dp], 1, 3), &
      [0.01_dp, 0.01_dp, 0.01_dp], [0.02_dp, 0.02_dp], &
      [2.0_dp, 2.0_dp, 2.0_dp], [0.01_dp, 0.01_dp, 0.01_dp], &
      [2.0_dp, 2.0_dp], [0.02_dp, 0.02_dp], 2, mbsts_state_normals, &
      mbsts_loading_normals, mbsts_gamma, mbsts_offset, &
      series_initial_mean=[0.1_dp, -0.2_dp, 0.3_dp], &
      series_initial_variance=[1.0_dp, 1.0_dp, 1.0_dp], &
      series_variance=[0.01_dp, 0.02_dp, 0.03_dp], &
      series_prior_shape=[2.0_dp, 2.0_dp, 2.0_dp], &
      series_prior_rate=[0.01_dp, 0.02_dp, 0.03_dp], &
      series_state_normal_draws=mbsts_series_normals, &
      series_gamma_draws=mbsts_series_gamma, &
      trend_initial_mean=reshape([0.0_dp, 0.01_dp, 0.0_dp, -0.01_dp, &
      0.0_dp, 0.02_dp], [2, 3]), &
      trend_initial_covariance=spread(identity2(), 3, 3), &
      trend_variance=reshape([0.01_dp, 0.005_dp, 0.01_dp, 0.005_dp, &
      0.01_dp, 0.005_dp], [2, 3]), &
      trend_prior_shape=reshape([2.0_dp, 2.0_dp, 2.0_dp, 2.0_dp, &
      2.0_dp, 2.0_dp], [2, 3]), &
      trend_prior_rate=reshape([0.01_dp, 0.005_dp, 0.01_dp, 0.005_dp, &
      0.01_dp, 0.005_dp], [2, 3]), &
      trend_state_normal_draws=mbsts_trend_normals, &
      trend_gamma_draws=mbsts_trend_gamma, seasonal_nseasons=4, &
      seasonal_duration=2, &
      seasonal_initial_variance=[1.0_dp, 1.0_dp, 1.0_dp], &
      seasonal_variance=[0.01_dp, 0.02_dp, 0.03_dp], &
      seasonal_prior_shape=[2.0_dp, 2.0_dp, 2.0_dp], &
      seasonal_prior_rate=[0.01_dp, 0.02_dp, 0.03_dp], &
      seasonal_state_normal_draws=mbsts_seasonal_normals, &
      seasonal_gamma_draws=mbsts_seasonal_gamma)
   call check(shared_fit%info == 0 .and. &
      shared_fit%has_series_local_level .and. &
      all(shape(shared_fit%series_state) == [3, 30, 8]) .and. &
      all(shape(shared_fit%series_variance) == [3, 8]) .and. &
      all(shared_fit%series_variance > 0.0_dp) .and. &
      maxval(abs(shared_fit%series_contribution - &
      shared_fit%series_state)) < 1.0e-12_dp, &
      'multivariate BSTS series local-level posterior')
   call check(shared_fit%has_series_local_linear_trend .and. &
      shared_fit%has_series_seasonal .and. &
      all(shape(shared_fit%series_trend_state) == [2, 3, 30, 8]) .and. &
      all(shape(shared_fit%series_trend_variance) == [2, 3, 8]) .and. &
      all(shape(shared_fit%series_seasonal_state) == [3, 3, 30, 8]) .and. &
      all(shape(shared_fit%series_seasonal_variance) == [3, 8]) .and. &
      all(shared_fit%series_trend_variance > 0.0_dp) .and. &
      all(shared_fit%series_seasonal_variance > 0.0_dp), &
      'multivariate BSTS series trend and seasonal posterior')
   do j = 1, 8
      call check(maxval(abs(shared_fit%fitted(:, :, j) - &
         shared_fit%contribution(:, :, j) - &
         shared_fit%regression_contribution(:, :, j) - &
         shared_fit%series_contribution(:, :, j) - &
         shared_fit%series_trend_contribution(:, :, j) - &
         shared_fit%series_seasonal_contribution(:, :, j) - &
         transpose(mbsts_offset(:, :, j)))) < 1.0e-12_dp, &
         'multivariate BSTS series decomposition')
   end do
   multivariate_prediction = bsts_mbsts_predict_draws(shared_fit, &
      mbsts_future_predictors, mbsts_forecast_state, &
      mbsts_forecast_observation, mbsts_series_forecast, &
      mbsts_trend_forecast, mbsts_seasonal_forecast)
   call check(multivariate_prediction%info == 0 .and. &
      all(shape(multivariate_prediction%draws) == [3, 4, 6]), &
      'multivariate BSTS series local-level forecast')
   do j = 1, 6
      do k = 1, 3
         expected_trend = shared_fit%series_trend_state(:, k, 30, j + 2)
         expected_season = shared_fit%series_seasonal_state(:, k, 30, j + 2)
         do i = 1, 4
            expected_trend(1) = expected_trend(1) + expected_trend(2)
            if (modulo(30 + i - 1, 2) == 0) then
               next_season = [-sum(expected_season), expected_season(1), &
                  expected_season(2)]
               expected_season = next_season
            end if
            call check(abs(multivariate_prediction%draws(k, i, j) - &
               dot_product(shared_fit%loadings(k, :, j + 2), &
               shared_fit%state(:, 30, j + 2)) - &
               dot_product(mbsts_future_predictors(i, :, k), &
               shared_fit%regression_coefficients(:, k, j + 2)) - &
               shared_fit%series_state(k, 30, j + 2) - expected_trend(1) - &
               expected_season(1)) < 1.0e-11_dp, &
               'multivariate BSTS series forecast identity')
         end do
      end do
   end do
   mbsts_inclusion_prior = reshape([1.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, &
      1.0_dp, 0.5_dp], [2, 3])
   mbsts_initial_inclusion = .false.
   mbsts_regression_normals = 0.0_dp
   mbsts_regression_uniforms = tiny(1.0_dp)
   shared_fit = bsts_mbsts_draws(mbsts_response, mbsts_predictors, &
      reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], &
      [2, 3]), spread(10.0_dp*identity2(), 3, 3), shared_factors(:, 1), &
      identity2(), shared_loadings, shared_loadings, &
      10.0_dp*spread([1.0_dp, 1.0_dp], 1, 3), &
      [0.01_dp, 0.01_dp, 0.01_dp], [0.02_dp, 0.02_dp], &
      [2.0_dp, 2.0_dp, 2.0_dp], [0.01_dp, 0.01_dp, 0.01_dp], &
      [2.0_dp, 2.0_dp], [0.02_dp, 0.02_dp], 2, mbsts_state_normals, &
      mbsts_loading_normals, mbsts_gamma, mbsts_offset, &
      mbsts_inclusion_prior, mbsts_initial_inclusion, &
      mbsts_regression_normals, mbsts_regression_uniforms, 2, -1)
   call check(shared_fit%info == 0 .and. &
      shared_fit%regression_spike_slab .and. &
      all(shape(shared_fit%regression_included) == [2, 3, 8]) .and. &
      all(shared_fit%regression_included(1, :, :)) .and. &
      all(.not. shared_fit%regression_included(2, 1, :)) .and. &
      all(shared_fit%regression_included(2, 2:, :)) .and. &
      all(shared_fit%regression_coefficients(2, 1, :) == 0.0_dp) .and. &
      shared_fit%regression_inclusion_probability(2, 1) == 0.0_dp .and. &
      shared_fit%regression_inclusion_probability(2, 2) == 1.0_dp .and. &
      shared_fit%regression_inclusion_probability(2, 3) == 1.0_dp, &
      'multivariate BSTS regression selection')
   shared_fit = bsts_mbsts_draws(mbsts_response, mbsts_predictors, &
      reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], &
      [2, 3]), spread(10.0_dp*identity2(), 3, 3), shared_factors(:, 1), &
      identity2(), shared_loadings, shared_loadings, &
      10.0_dp*spread([1.0_dp, 1.0_dp], 1, 3), &
      [0.01_dp, 0.01_dp, 0.01_dp], [0.02_dp, 0.02_dp], &
      [2.0_dp, 2.0_dp, 2.0_dp], [0.01_dp, 0.01_dp, 0.01_dp], &
      [2.0_dp, 2.0_dp], [0.02_dp, 0.02_dp], 2, mbsts_state_normals, &
      mbsts_loading_normals, mbsts_gamma, mbsts_offset, &
      mbsts_inclusion_prior, mbsts_initial_inclusion, &
      mbsts_regression_normals, mbsts_regression_uniforms, 2, 0)
   call check(shared_fit%info == 0 .and. &
      all(.not. shared_fit%regression_included(2, 3, :)) .and. &
      all(shared_fit%regression_coefficients(2, 3, :) == 0.0_dp), &
      'multivariate BSTS regression flip limit')
   mbsts_inclusion_prior = reshape([1.0_dp, 0.5_dp, 1.0_dp, 0.5_dp, &
      1.0_dp, 0.5_dp], [2, 3])
   shared_fit = bsts_mbsts_draws(mbsts_response, mbsts_predictors, &
      reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], &
      [2, 3]), spread(10.0_dp*identity2(), 3, 3), shared_factors(:, 1), &
      identity2(), shared_loadings, shared_loadings, &
      10.0_dp*spread([1.0_dp, 1.0_dp], 1, 3), &
      [0.01_dp, 0.01_dp, 0.01_dp], [0.02_dp, 0.02_dp], &
      [2.0_dp, 2.0_dp, 2.0_dp], [0.01_dp, 0.01_dp, 0.01_dp], &
      [2.0_dp, 2.0_dp], [0.02_dp, 0.02_dp], 2, mbsts_state_normals, &
      mbsts_loading_normals, mbsts_gamma, mbsts_offset, &
      mbsts_inclusion_prior, mbsts_initial_inclusion, &
      mbsts_regression_normals, mbsts_regression_uniforms, 1, -1)
   call check(shared_fit%info == 0, &
      'multivariate BSTS regression model-size status')
   call check(all(shared_fit%regression_included(1, :, :)) .and. &
      all(.not. shared_fit%regression_included(2, :, :)), &
      'multivariate BSTS regression model-size indicators')
   call check(all(shared_fit%regression_coefficients(2, :, :) == 0.0_dp), &
      'multivariate BSTS regression model-size coefficients')
   shared_fit = bsts_mbsts(mbsts_response, mbsts_predictors, 2, 7, burn=2, &
      offset_draws=mbsts_offset(:, :, :7), series_local_level=.true., &
      series_local_linear_trend=.true., series_nseasons=4, &
      series_season_duration=2)
   call check(shared_fit%info == 0 .and. shared_fit%is_mbsts .and. &
      shared_fit%regression_spike_slab .and. &
      allocated(shared_fit%regression_included) .and. &
      shared_fit%has_series_local_level .and. &
      allocated(shared_fit%series_state) .and. &
      shared_fit%has_series_local_linear_trend .and. &
      shared_fit%has_series_seasonal, &
      'random multivariate BSTS wrapper')
   multivariate_prediction = bsts_mbsts_predict(shared_fit, &
      mbsts_future_predictors)
   call check(multivariate_prediction%info == 0 .and. &
      all(shape(multivariate_prediction%draws) == [3, 4, 5]), &
      'random multivariate BSTS forecast wrapper')
   shared_fit = bsts_mbsts(mbsts_response, mbsts_predictors, 2, 3, burn=1, &
      select_regression=.false.)
   call check(shared_fit%info == 0 .and. shared_fit%is_mbsts .and. &
      .not. shared_fit%regression_spike_slab .and. &
      .not. allocated(shared_fit%regression_included), &
      'random dense multivariate BSTS wrapper')

   allocate(forecast_state(2, 3, 9), forecast_observation(3, 9))
   do j = 1, 9
      do i = 1, 3
         forecast_state(1, i, j) = 0.2_dp*sin(real(i + j, dp))
         forecast_state(2, i, j) = 0.2_dp*cos(real(2*i + j, dp))
         forecast_observation(i, j) = 0.2_dp*sin(real(3*i + j, dp))
      end do
   end do
   prediction = bsts_predict_draws(trend_fit, 3, forecast_state, &
      forecast_observation)
   call check(prediction%info == 0 .and. &
      all(shape(prediction%draws) == [3, 9]), 'posterior forecast status')
   call check(all(prediction%standard_deviation > 0.0_dp) .and. &
      all(prediction%lower <= prediction%mean) .and. &
      all(prediction%upper >= prediction%mean), 'posterior forecast summaries')
   prediction = bsts_predict(seasonal_fit, 4)
   call check(prediction%info == 0 .and. &
      all(shape(prediction%draws) == [4, 10]), &
      'duration-aware seasonal forecast')

   random_fit = bsts_local_level(y, 10, burn=2)
   call check(random_fit%info == 0 .and. &
      all(random_fit%observation_variance > 0.0_dp), &
      'random local-level wrapper')
   random_fit = bsts_seasonal(seasonal_y, 4, 8, season_duration=2, burn=2)
   call check(random_fit%info == 0 .and. &
      all(random_fit%component_variance > 0.0_dp), &
      'random seasonal wrapper')
   random_fit = bsts_trig(trig_y, 12.0_dp, [1.0_dp, 2.0_dp], 8, burn=2)
   call check(random_fit%info == 0 .and. &
      all(random_fit%component_variance > 0.0_dp), &
      'random trigonometric wrapper')
   random_fit = bsts_semilocal_trend(semilocal_y, 8, burn=2, &
      force_stationary=.true., force_positive=.true.)
   call check(random_fit%info == 0 .and. &
      all(random_fit%slope_ar > 0.0_dp) .and. &
      all(random_fit%slope_ar < 1.0_dp), 'random semilocal wrapper')

   allocate(regression_y(60), regression_x(60, 3), &
      regression_normals(3, 30), inclusion_uniforms(3, 30), &
      regression_gamma(30), regression_offset(60, 30))
   do i = 1, 60
      regression_x(i, 1) = sin(0.31_dp*real(i, dp))
      regression_x(i, 2) = cos(0.73_dp*real(i, dp))
      regression_x(i, 3) = (real(i, dp) - 30.5_dp)/30.0_dp
      regression_y(i) = 2.0_dp + 1.5_dp*regression_x(i, 1) - &
         0.9_dp*regression_x(i, 3) + 0.03_dp*sin(1.7_dp*real(i, dp))
   end do
   regression_offset = 2.0_dp
   do j = 1, 30
      regression_normals(:, j) = 0.1_dp*[sin(real(j, dp)), &
         cos(real(2*j, dp)), sin(real(3*j, dp))]
      inclusion_uniforms(:, j) = [0.35_dp, 0.65_dp, 0.4_dp]
      regression_gamma(j) = 31.0_dp + 0.05_dp*real(j, dp)
   end do
   regression_fit = bsts_spike_slab_draws(regression_y, regression_x, &
      [0.0_dp, 0.0_dp, 0.0_dp], reshape([10.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 10.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 10.0_dp], [3, 3]), &
      [0.5_dp, 0.5_dp, 0.5_dp], 2.0_dp, 0.02_dp, 0.2_dp, &
      [.false., .false., .false.], 2, 3, 10, regression_normals, &
      inclusion_uniforms, regression_gamma, regression_offset)
   call check(regression_fit%info == 0 .and. &
      all(shape(regression_fit%coefficients) == [3, 30]) .and. &
      all(regression_fit%residual_variance > 0.0_dp), &
      'static spike-and-slab status')
   call check(regression_fit%inclusion_probability(1) > 0.8_dp .and. &
      regression_fit%inclusion_probability(3) > 0.8_dp .and. &
      maxval([(count(regression_fit%included(:, j)), j=1, 30)]) <= 2 .and. &
      sum(regression_fit%coefficients(1, 11:))/20.0_dp > 1.0_dp .and. &
      sum(regression_fit%coefficients(3, 11:))/20.0_dp < -0.5_dp, &
      'static spike-and-slab signal selection')
   prediction_errors = bsts_regression_prediction_errors(regression_fit, &
      regression_y - 2.0_dp, regression_x)
   scaled_errors = bsts_regression_prediction_errors(regression_fit, &
      regression_y - 2.0_dp, regression_x, .true.)
   call check(prediction_errors%info == 0 .and. scaled_errors%info == 0 .and. &
      all(shape(prediction_errors%draws) == [60, 20]) .and. &
      maxval(abs(scaled_errors%draws*sqrt( &
      prediction_errors%forecast_variance) - prediction_errors%draws)) < &
      1.0e-12_dp, 'static regression prediction errors')
   allocate(regression_prediction_normals(2, 20))
   regression_prediction_normals = 0.0_dp
   prediction = bsts_regression_predict_draws(regression_fit, reshape([ &
      1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 3]), &
      regression_prediction_normals)
   call check(prediction%info == 0 .and. &
      all(shape(prediction%draws) == [2, 20]) .and. &
      prediction%mean(1) > 1.0_dp .and. prediction%mean(2) < -0.5_dp, &
      'static regression posterior prediction')
   random_regression_fit = bsts_spike_slab(regression_y - 2.0_dp, &
      regression_x, 12, burn=2, expected_model_size=2.0_dp, &
      maximum_model_size=2)
   call check(random_regression_fit%info == 0 .and. &
      all(random_regression_fit%residual_variance > 0.0_dp), &
      'random spike-and-slab wrapper')
   prediction = bsts_regression_predict(random_regression_fit, &
      regression_x(:2, :))
   call check(prediction%info == 0 .and. &
      all(shape(prediction%draws) == [2, 10]), &
      'random regression prediction wrapper')

   allocate(dynamic_y(50), dynamic_x(50, 2), &
      dynamic_normals(2, 50, 16), dynamic_gamma(4, 16))
   do i = 1, 50
      dynamic_x(i, 1) = sin(0.23_dp*real(i, dp))
      dynamic_x(i, 2) = cos(0.37_dp*real(i, dp))
      dynamic_y(i) = (0.6_dp + 0.006_dp*real(i, dp))*dynamic_x(i, 1) + &
         (-0.5_dp + 0.003_dp*real(i, dp))*dynamic_x(i, 2) + &
         0.015_dp*sin(1.3_dp*real(i, dp))
   end do
   do j = 1, 16
      do i = 1, 50
         dynamic_normals(:, i, j) = 0.08_dp*[sin(real(i + j, dp)), &
            cos(real(2*i + j, dp))]
      end do
      dynamic_gamma(:, j) = [27.0_dp, 26.0_dp, 26.0_dp, 21.0_dp] + &
         0.02_dp*real(j, dp)
   end do
   dynamic_fit = bsts_dynamic_regression_draws(dynamic_y, dynamic_x, &
      [0.0_dp, 0.0_dp], reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2]), &
      0.05_dp, [0.001_dp, 0.001_dp], 2.0_dp, 0.02_dp, [2.0_dp, 2.0_dp], &
      [0.001_dp, 0.001_dp], 4, dynamic_normals, dynamic_gamma(:3, :))
   call check(dynamic_fit%info == 0 .and. &
      all(shape(dynamic_fit%coefficients) == [2, 50, 16]) .and. &
      all(dynamic_fit%innovation_variance > 0.0_dp) .and. &
      all(dynamic_fit%residual_variance > 0.0_dp), &
      'independent dynamic regression status')
   call check(sum(dynamic_fit%coefficients(1, 50, 5:))/12.0_dp > 0.5_dp .and. &
      sum(dynamic_fit%coefficients(2, 50, 5:))/12.0_dp < -0.2_dp, &
      'dynamic regression coefficient paths')
   allocate(dynamic_forecast_state(2, 2, 12), &
      dynamic_forecast_observation(2, 12))
   dynamic_forecast_state = 0.0_dp
   dynamic_forecast_observation = 0.0_dp
   prediction = bsts_dynamic_regression_predict_draws(dynamic_fit, &
      dynamic_x(:2, :), dynamic_forecast_state, dynamic_forecast_observation)
   call check(prediction%info == 0 .and. &
      all(shape(prediction%draws) == [2, 12]), &
      'dynamic regression supplied-draw forecast')
   do j = 1, 12
      call check(maxval(abs(prediction%draws(:, j) - &
         matmul(dynamic_x(:2, :), dynamic_fit%coefficients(:, 50, j + 4)))) < &
         1.0e-12_dp, 'dynamic regression zero-shock forecast')
   end do
   hierarchical_fit = bsts_dynamic_regression_hierarchical_draws(dynamic_y, &
      dynamic_x, [0.0_dp, 0.0_dp], &
      reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2]), 0.05_dp, &
      [0.001_dp, 0.001_dp], [0.5_dp, 0.5_dp], 2.0_dp, 0.02_dp, 10.0_dp, &
      0.009_dp, 1.0_dp, 1.0_dp, 4, dynamic_normals, dynamic_gamma)
   call check(hierarchical_fit%info == 0 .and. hierarchical_fit%hierarchical .and. &
      all(hierarchical_fit%hierarchy_rate > 0.0_dp) .and. &
      maxval(abs(hierarchical_fit%scaled_innovation_variance - &
      0.5_dp*hierarchical_fit%innovation_variance)) < 1.0e-12_dp, &
      'hierarchical dynamic regression scaling')
   dynamic_fit = bsts_dynamic_regression(dynamic_y, dynamic_x, 8, burn=2)
   call check(dynamic_fit%info == 0 .and. &
      all(dynamic_fit%innovation_variance > 0.0_dp), &
      'random independent dynamic regression wrapper')
   hierarchical_fit = bsts_dynamic_regression_hierarchical(dynamic_y, &
      dynamic_x, 8, burn=2)
   call check(hierarchical_fit%info == 0 .and. &
      all(hierarchical_fit%hierarchy_rate > 0.0_dp), &
      'random hierarchical dynamic regression wrapper')
   prediction = bsts_dynamic_regression_predict(hierarchical_fit, &
      dynamic_x(:2, :))
   call check(prediction%info == 0 .and. &
      all(shape(prediction%draws) == [2, 6]), &
      'random dynamic regression forecast wrapper')

   allocate(ar_state_y(60), ar_state_normals(2, 60, 14), &
      ar_state_parameter_normals(2, 30, 14), ar_state_gamma(2, 14))
   beta_first = 0.7_dp
   previous_first = 0.4_dp
   do i = 1, 60
      slope = beta_first
      beta_first = 0.68_dp*beta_first - 0.16_dp*previous_first + &
         0.05_dp*sin(0.37_dp*real(i, dp))
      previous_first = slope
      ar_state_y(i) = beta_first + 0.01_dp*cos(1.2_dp*real(i, dp))
   end do
   do j = 1, 14
      do i = 1, 60
         ar_state_normals(:, i, j) = 0.05_dp*[ &
            sin(real(i + j, dp)), cos(real(2*i + j, dp))]
      end do
      do i = 1, 30
         ar_state_parameter_normals(:, i, j) = 0.03_dp*[ &
            sin(real(i + j, dp)), cos(real(2*i + j, dp))]
      end do
      ar_state_gamma(:, j) = [32.0_dp, 31.0_dp] + 0.02_dp*real(j, dp)
   end do
   ar_fit = bsts_ar_draws(ar_state_y, [0.0_dp, 0.0_dp], &
      reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2]), &
      [0.5_dp, 0.0_dp], 0.03_dp, 0.003_dp, 2.0_dp, 0.02_dp, &
      2.0_dp, 0.002_dp, 4, ar_state_normals, &
      ar_state_parameter_normals, ar_state_gamma)
   call check(ar_fit%info == 0 .and. ar_fit%ar_order == 2 .and. &
      all(shape(ar_fit%state) == [2, 60, 14]) .and. &
      all(shape(ar_fit%ar_coefficients) == [2, 14]) .and. &
      all(ar_fit%component_variance > 0.0_dp) .and. &
      all(ieee_is_finite(ar_fit%state)), 'Bayesian AR(2) state status')
   allocate(ar_state_forecast(2, 10))
   ar_state_forecast = 0.0_dp
   prediction = bsts_ar_predict_draws(ar_fit, 2, ar_state_forecast, &
      ar_state_forecast)
   call check(prediction%info == 0 .and. &
      all(shape(prediction%draws) == [2, 10]), &
      'Bayesian AR supplied-draw forecast')
   do j = 1, 10
      call check(abs(prediction%draws(1, j) - dot_product( &
         ar_fit%ar_coefficients(:, j + 4), ar_fit%state(:, 60, j + 4))) < &
         1.0e-12_dp, 'Bayesian AR zero-shock forecast')
   end do
   ar_fit = bsts_ar(ar_state_y, 2, 8, burn=2, proposal_attempts=30)
   call check(ar_fit%info == 0 .and. ar_fit%ar_order == 2, &
      'random Bayesian AR wrapper')
   prediction = bsts_ar_predict(ar_fit, 3)
   call check(prediction%info == 0 .and. &
      all(shape(prediction%draws) == [3, 6]), &
      'random Bayesian AR forecast wrapper')
   allocate(auto_ar_uniforms(2, 14))
   auto_ar_uniforms = 0.5_dp
   ar_fit = bsts_auto_ar_draws(ar_state_y, [0.0_dp, 0.0_dp], &
      reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2]), &
      [0.0_dp, 0.0_dp], [.false., .false.], 0.03_dp, 0.003_dp, &
      [0.0_dp, 0.0_dp], reshape([0.25_dp, 0.0_dp, 0.0_dp, 0.16_dp], &
      [2, 2]), [1.0_dp, 0.0_dp], 2.0_dp, 0.02_dp, 2.0_dp, 0.002_dp, &
      -1, .true., 4, ar_state_normals, ar_state_parameter_normals, &
      auto_ar_uniforms, ar_state_gamma)
   call check(ar_fit%info == 0 .and. &
      all(ar_fit%ar_included(1, :)) .and. &
      .not. any(ar_fit%ar_included(2, :)) .and. &
      ar_fit%ar_inclusion_probability(1) == 1.0_dp .and. &
      ar_fit%ar_inclusion_probability(2) == 0.0_dp .and. &
      all(ar_fit%ar_coefficients(2, :) == 0.0_dp), &
      'automatic AR forced lag selection')
   prediction = bsts_ar_predict(ar_fit, 2)
   call check(prediction%info == 0 .and. &
      all(shape(prediction%draws) == [2, 10]), &
      'automatic AR posterior forecast')
   ar_fit = bsts_auto_ar(ar_state_y, 4, 8, burn=2, maximum_flips=2, &
      proposal_attempts=100)
   call check(ar_fit%info == 0 .and. ar_fit%ar_order == 4 .and. &
      all(ar_fit%ar_inclusion_probability >= 0.0_dp) .and. &
      all(ar_fit%ar_inclusion_probability <= 1.0_dp), &
      'random automatic AR wrapper')

   allocate(ar_dynamic_y(50), ar_dynamic_normals(4, 50, 14), &
      ar_parameter_normals(2, 2, 30, 14), ar_dynamic_gamma(3, 14))
   beta_first = 0.8_dp
   previous_first = 0.5_dp
   beta_second = -0.6_dp
   previous_second = -0.3_dp
   do i = 1, 50
      slope = beta_first
      beta_first = 0.72_dp*beta_first - 0.18_dp*previous_first + &
         0.04_dp*sin(0.41_dp*real(i, dp))
      previous_first = slope
      slope = beta_second
      beta_second = 0.62_dp*beta_second - 0.12_dp*previous_second + &
         0.035_dp*cos(0.29_dp*real(i, dp))
      previous_second = slope
      ar_dynamic_y(i) = beta_first*dynamic_x(i, 1) + &
         beta_second*dynamic_x(i, 2) + 0.01_dp*sin(1.1_dp*real(i, dp))
   end do
   do j = 1, 14
      do i = 1, 50
         ar_dynamic_normals(:, i, j) = 0.04_dp*[ &
            sin(real(i + j, dp)), cos(real(2*i + j, dp)), &
            sin(real(3*i + j, dp)), cos(real(i + 2*j, dp))]
      end do
      do i = 1, 30
         ar_parameter_normals(:, 1, i, j) = 0.03_dp*[ &
            sin(real(i + j, dp)), cos(real(2*i + j, dp))]
         ar_parameter_normals(:, 2, i, j) = 0.03_dp*[ &
            cos(real(i + 2*j, dp)), sin(real(3*i + j, dp))]
      end do
      ar_dynamic_gamma(:, j) = [27.0_dp, 26.0_dp, 26.0_dp] + &
         0.02_dp*real(j, dp)
   end do
   dynamic_fit = bsts_dynamic_regression_ar_draws(ar_dynamic_y, dynamic_x, 2, &
      [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], &
      reshape([1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
      0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [4, 4]), &
      reshape([0.5_dp, 0.0_dp, 0.5_dp, 0.0_dp], [2, 2]), &
      0.03_dp, [0.001_dp, 0.001_dp], [0.5_dp, 0.5_dp], 2.0_dp, &
      0.02_dp, [2.0_dp, 2.0_dp], [0.001_dp, 0.001_dp], 4, &
      ar_dynamic_normals, ar_parameter_normals, ar_dynamic_gamma)
   call check(dynamic_fit%info == 0 .and. dynamic_fit%ar_order == 2 .and. &
      all(shape(dynamic_fit%ar_coefficients) == [2, 2, 14]) .and. &
      all(dynamic_fit%innovation_variance > 0.0_dp) .and. &
      all(ieee_is_finite(dynamic_fit%coefficients)), &
      'AR(2) dynamic regression status')
   dynamic_forecast_state = 0.0_dp
   dynamic_forecast_observation = 0.0_dp
   prediction = bsts_dynamic_regression_ar_predict_draws(dynamic_fit, &
      dynamic_x(:2, :), dynamic_forecast_state(:, :, :10), &
      dynamic_forecast_observation(:, :10))
   call check(prediction%info == 0 .and. &
      all(shape(prediction%draws) == [2, 10]), &
      'AR dynamic regression supplied-draw forecast')
   do j = 1, 10
      beta_first = dot_product(dynamic_fit%ar_coefficients(:, 1, j + 4), &
         dynamic_fit%coefficients(1, 50:49:-1, j + 4))
      beta_second = dot_product(dynamic_fit%ar_coefficients(:, 2, j + 4), &
         dynamic_fit%coefficients(2, 50:49:-1, j + 4))
      call check(abs(prediction%draws(1, j) - &
         dot_product(dynamic_x(1, :), [beta_first, beta_second])) < 1.0e-12_dp, &
         'AR dynamic regression zero-shock forecast')
   end do
   dynamic_fit = bsts_dynamic_regression_ar(ar_dynamic_y, dynamic_x, 2, 7, &
      burn=2, proposal_attempts=30)
   call check(dynamic_fit%info == 0 .and. dynamic_fit%ar_order == 2, &
      'random AR dynamic regression wrapper')
   prediction = bsts_dynamic_regression_ar_predict(dynamic_fit, dynamic_x(:2, :))
   call check(prediction%info == 0 .and. &
      all(shape(prediction%draws) == [2, 5]), &
      'random AR dynamic regression forecast wrapper')

contains

   pure function identity2() result(matrix)
      !! Return a two-dimensional identity matrix for factor tests.
      real(dp) :: matrix(2, 2)

      matrix = reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2])
   end function identity2

   pure function identity3() result(matrix)
      !! Return a three-dimensional identity matrix for prior tests.
      real(dp) :: matrix(3, 3)

      matrix = reshape([1.0_dp, 0.0_dp, 0.0_dp, &
         0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [3, 3])
   end function identity3

   subroutine check(condition, label)
      !! Stop the test when an assertion fails.
      logical, intent(in) :: condition !! Flag controlling condition.
      character(*), intent(in) :: label !! Label.

      if (.not. condition) then
         write (*, '(a)') 'FAILED: '//label
         error stop 1
      end if
   end subroutine check

end program test_bsts
