! SPDX-License-Identifier: MIT
! SPDX-FileComment: Algorithms translated from the R bsts package under its MIT option.
! Bayesian structural time-series algorithms translated from bsts.
module bsts_mod
   use kind_mod, only: dp
   use kfas_mod, only: ssm_model_t, kfs_filter_t, kfs_filter
   use linalg_mod, only: invert_matrix, inverse_logdet, symmetrize, &
      symmetric_eigen, cholesky_lower
   use random_mod, only: random_gamma, &
      random_standard_normal, random_uniform, multivariate_normal_from_standard
   use stats_mod, only: normal_quantile, quantile
   use calendar_mod, only: date_t, date_valid, date_day_number, &
      date_from_day_number, date_day_of_week, date_days_in_month, date_easter
   use calendar_mod, only: operator(+), operator(-), operator(==), &
      operator(<)
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, &
      ieee_quiet_nan
   implicit none
   private

   integer, parameter, public :: bsts_holiday_fixed = 1
   integer, parameter, public :: bsts_holiday_nth_weekday = 2
   integer, parameter, public :: bsts_holiday_last_weekday = 3
   integer, parameter, public :: bsts_holiday_date_range = 4
   integer, parameter, public :: bsts_holiday_named = 5

   type, public :: bsts_holiday_t
      ! Calendar definition and influence window for one holiday.
      integer :: kind = 0
      integer :: month = 0, day = 0, weekday = 0, week_number = 0
      integer :: days_before = 0, days_after = 0
      character(len=40) :: name = ''
      type(date_t), allocatable :: range_start(:), range_end(:)
      integer :: info = 0
   end type bsts_holiday_t

   type, public :: bsts_mcmc_t
      ! Posterior draws from a Gaussian structural time-series model.
      real(dp), allocatable :: state(:, :, :)
      real(dp), allocatable :: observation_variance(:)
      real(dp), allocatable :: state_variance(:, :)
      real(dp), allocatable :: component_variance(:, :)
      real(dp), allocatable :: transition(:, :), observation(:)
      real(dp), allocatable :: transition_schedule(:, :, :)
      real(dp), allocatable :: state_loading_schedule(:, :, :)
      real(dp), allocatable :: observation_schedule(:, :)
      real(dp), allocatable :: slope_mean(:), slope_ar(:)
      real(dp), allocatable :: ar_coefficients(:, :)
      real(dp), allocatable :: ar_inclusion_probability(:)
      logical, allocatable :: ar_included(:, :)
      real(dp), allocatable :: degrees_of_freedom(:, :)
      real(dp), allocatable :: state_weights(:, :, :)
      integer :: forecast_phase = 1
      integer :: ar_order = 0
      integer :: burn = 0
      integer :: info = 0
      logical :: is_semilocal = .false.
      logical :: is_monthly = .false.
      type(date_t) :: last_date
      type(bsts_holiday_t) :: holiday
      logical :: is_holiday = .false.
   end type bsts_mcmc_t

   type, public :: bsts_prediction_t
      ! Posterior predictive draws and pointwise summaries.
      real(dp), allocatable :: draws(:, :), mean(:), standard_deviation(:)
      real(dp), allocatable :: lower(:), upper(:)
      integer :: info = 0
   end type bsts_prediction_t

   type, public :: bsts_multivariate_prediction_t
      ! Multivariate posterior predictive draws and pointwise summaries.
      real(dp), allocatable :: draws(:, :, :), mean(:, :)
      real(dp), allocatable :: standard_deviation(:, :)
      real(dp), allocatable :: lower(:, :), upper(:, :)
      integer :: info = 0
   end type bsts_multivariate_prediction_t

   type, public :: bsts_spike_slab_t
      ! Static Gaussian spike-and-slab regression posterior draws.
      real(dp), allocatable :: coefficients(:, :), residual_variance(:)
      real(dp), allocatable :: inclusion_probability(:)
      real(dp), allocatable :: log_model_probability(:)
      logical, allocatable :: included(:, :)
      integer :: burn = 0
      integer :: info = 0
   end type bsts_spike_slab_t

   type, public :: bsts_dirm_t
      ! Posterior draws from a Gaussian dynamic-intercept regression model.
      real(dp), allocatable :: state(:, :), coefficients(:, :)
      real(dp), allocatable :: observation_variance(:), level_variance(:)
      real(dp), allocatable :: contribution(:, :), fitted(:, :), residuals(:, :)
      real(dp), allocatable :: inclusion_probability(:)
      logical, allocatable :: included(:, :)
      integer :: time_points = 0
      integer :: burn = 0
      integer :: info = 0
   end type bsts_dirm_t

   type, public :: bsts_mixed_t
      ! Posterior draws from a mixed-frequency Gaussian regression model.
      real(dp), allocatable :: latent_fine(:, :), cumulator(:, :)
      real(dp), allocatable :: state(:, :), coefficients(:, :)
      real(dp), allocatable :: observation_variance(:), level_variance(:)
      real(dp), allocatable :: state_contribution(:, :)
      real(dp), allocatable :: regression_contribution(:, :)
      real(dp), allocatable :: coarse_fitted(:, :)
      real(dp), allocatable :: structural_state(:, :, :)
      real(dp), allocatable :: component_variance(:, :)
      real(dp), allocatable :: inclusion_probability(:)
      logical, allocatable :: included(:, :)
      integer :: nseasons = 0
      integer :: season_duration = 1
      integer :: burn = 0
      integer :: info = 0
   end type bsts_mixed_t

   type, public :: bsts_mixed_prediction_t
      ! Fine- and coarse-frequency posterior predictive distributions.
      type(bsts_prediction_t) :: fine
      type(bsts_prediction_t) :: coarse
      integer :: info = 0
   end type bsts_mixed_prediction_t

   type, public :: bsts_prediction_errors_t
      ! Posterior one-step-ahead Gaussian prediction errors.
      real(dp), allocatable :: draws(:, :), forecast_variance(:, :)
      real(dp), allocatable :: mean(:)
      real(dp) :: rmse = 0.0_dp
      real(dp) :: mae = 0.0_dp
      integer :: training_size = 0
      integer :: info = 0
      logical :: standardized = .false.
   end type bsts_prediction_errors_t

   type, public :: bsts_model_comparison_t
      ! Prediction-error scores for a collection of compatible models.
      real(dp), allocatable :: cumulative_absolute_error(:, :)
      real(dp), allocatable :: rmse(:), mae(:), final_absolute_error(:)
      integer, allocatable :: rmse_rank(:), mae_rank(:)
      integer :: best_rmse = 0
      integer :: best_mae = 0
      integer :: start_index = 1
      integer :: info = 0
      logical :: standardized = .false.
   end type bsts_model_comparison_t

   type, public :: bsts_static_intercept_t
      ! Posterior draws from a constant Gaussian state component.
      real(dp), allocatable :: intercept(:), observation_variance(:)
      real(dp), allocatable :: contribution(:, :), fitted(:, :), residuals(:, :)
      integer :: burn = 0
      integer :: info = 0
   end type bsts_static_intercept_t

   type, public :: bsts_numeric_timestamps_t
      ! Regular numeric timestamp grid and observation mapping.
      real(dp), allocatable :: grid(:)
      integer, allocatable :: mapping(:)
      logical :: trivial = .false.
      integer :: info = 0
   end type bsts_numeric_timestamps_t

   type, public :: bsts_date_timestamps_t
      ! Regular Gregorian date grid and observation mapping.
      type(date_t), allocatable :: grid(:)
      integer, allocatable :: mapping(:)
      logical :: trivial = .false.
      integer :: info = 0
   end type bsts_date_timestamps_t

   type, public :: bsts_wide_series_t
      ! Wide time-by-series representation with integer labels.
      real(dp), allocatable :: values(:, :)
      integer, allocatable :: timestamps(:), series_id(:)
      integer :: info = 0
   end type bsts_wide_series_t

   type, public :: bsts_long_series_t
      ! Long observation representation with integer labels.
      real(dp), allocatable :: values(:)
      integer, allocatable :: timestamps(:), series_id(:)
      integer :: info = 0
   end type bsts_long_series_t

   type, public :: bsts_monthly_series_t
      ! Monthly aggregates and their month-end dates.
      real(dp), allocatable :: values(:, :)
      type(date_t), allocatable :: dates(:)
      integer :: info = 0
   end type bsts_monthly_series_t

   type, public :: bsts_dynamic_regression_t
      ! Random-walk dynamic-regression posterior draws.
      real(dp), allocatable :: coefficients(:, :, :)
      real(dp), allocatable :: innovation_variance(:, :)
      real(dp), allocatable :: scaled_innovation_variance(:, :)
      real(dp), allocatable :: residual_variance(:), hierarchy_rate(:)
      real(dp), allocatable :: predictor_variance(:)
      real(dp), allocatable :: ar_coefficients(:, :, :)
      integer :: burn = 0
      integer :: ar_order = 0
      integer :: info = 0
      logical :: hierarchical = .false.
   end type bsts_dynamic_regression_t

   type, public :: bsts_holiday_regression_t
      ! Posterior draws from a fixed-pattern regression holiday model.
      real(dp), allocatable :: coefficients(:, :), residual_variance(:)
      real(dp), allocatable :: contribution(:, :)
      real(dp), allocatable :: coefficient_mean(:, :)
      real(dp), allocatable :: coefficient_variance(:, :, :)
      integer, allocatable :: coefficient_offset(:)
      type(bsts_holiday_t), allocatable :: holidays(:)
      type(date_t) :: first_date, last_date
      integer :: burn = 0
      integer :: info = 0
      logical :: hierarchical = .false.
   end type bsts_holiday_regression_t

   type, public :: bsts_shared_local_level_t
      ! Posterior draws from a shared multivariate local-level model.
      real(dp), allocatable :: state(:, :, :), loadings(:, :, :)
      real(dp), allocatable :: observation_variance(:, :)
      real(dp), allocatable :: factor_variance(:, :)
      real(dp), allocatable :: contribution(:, :, :)
      real(dp), allocatable :: inclusion_probability(:, :)
      real(dp), allocatable :: regression_coefficients(:, :, :)
      real(dp), allocatable :: regression_contribution(:, :, :)
      real(dp), allocatable :: regression_inclusion_probability(:, :)
      real(dp), allocatable :: series_state(:, :, :)
      real(dp), allocatable :: series_variance(:, :)
      real(dp), allocatable :: series_contribution(:, :, :)
      real(dp), allocatable :: series_trend_state(:, :, :, :)
      real(dp), allocatable :: series_trend_variance(:, :, :)
      real(dp), allocatable :: series_trend_contribution(:, :, :)
      real(dp), allocatable :: series_seasonal_state(:, :, :, :)
      real(dp), allocatable :: series_seasonal_variance(:, :)
      real(dp), allocatable :: series_seasonal_contribution(:, :, :)
      real(dp), allocatable :: fitted(:, :, :), residuals(:, :, :)
      logical, allocatable :: included(:, :, :)
      logical, allocatable :: regression_included(:, :, :)
      integer :: burn = 0
      integer :: info = 0
      logical :: spike_slab = .false.
      logical :: regression_spike_slab = .false.
      logical :: is_mbsts = .false.
      logical :: has_series_local_level = .false.
      logical :: has_series_local_linear_trend = .false.
      logical :: has_series_seasonal = .false.
      integer :: series_seasons = 0
      integer :: series_season_duration = 1
   end type bsts_shared_local_level_t

   type, public :: bsts_non_gaussian_t
      ! Posterior draws from a non-Gaussian local-level model.
      real(dp), allocatable :: state(:, :), fitted_mean(:, :)
      real(dp), allocatable :: state_variance(:)
      real(dp), allocatable :: coefficients(:, :)
      real(dp), allocatable :: regression_contribution(:, :)
      real(dp), allocatable :: inclusion_probability(:)
      logical, allocatable :: included(:, :)
      real(dp), allocatable :: structural_state(:, :, :)
      real(dp), allocatable :: component_variance(:, :)
      integer :: nseasons = 0
      integer :: season_duration = 1
      integer :: burn = 0
      integer :: info = 0
      integer :: family = 0
   end type bsts_non_gaussian_t

   public :: bsts_local_level_draws, bsts_local_level
   public :: bsts_local_linear_trend_draws, bsts_local_linear_trend
   public :: bsts_seasonal_draws, bsts_seasonal
   public :: bsts_trig_draws, bsts_trig
   public :: bsts_semilocal_trend_draws, bsts_semilocal_trend
   public :: bsts_spike_slab_draws, bsts_spike_slab
   public :: bsts_regression_predict_draws, bsts_regression_predict
   public :: bsts_dirm_draws, bsts_dirm
   public :: bsts_dirm_predict_draws, bsts_dirm_predict
   public :: bsts_mixed_draws, bsts_mixed
   public :: bsts_mixed_predict_draws, bsts_mixed_predict
   public :: bsts_mixed_trend_seasonal_draws
   public :: bsts_mixed_trend_seasonal
   public :: bsts_mixed_trend_seasonal_predict_draws
   public :: bsts_mixed_trend_seasonal_predict
   public :: bsts_structural_prediction_errors
   public :: bsts_regression_prediction_errors
   public :: bsts_local_level_holdout_errors
   public :: bsts_compare_prediction_errors
   public :: bsts_static_intercept_draws, bsts_static_intercept
   public :: bsts_static_intercept_predict_draws
   public :: bsts_static_intercept_predict
   public :: bsts_no_duplicates_numeric, bsts_no_duplicates_date
   public :: bsts_no_gaps_numeric, bsts_no_gaps_date
   public :: bsts_is_regular_numeric, bsts_is_regular_date
   public :: bsts_regularize_numeric_timestamps
   public :: bsts_regularize_date_timestamps
   public :: bsts_long_to_wide, bsts_wide_to_long
   public :: bsts_aggregate_time_series
   public :: bsts_month_distance, bsts_week_ends_month
   public :: bsts_week_ends_quarter, bsts_fraction_initial_month
   public :: bsts_match_week_to_month, bsts_aggregate_weeks_to_months

   interface bsts_aggregate_time_series
      module procedure aggregate_time_series_vector
      module procedure aggregate_time_series_matrix
   end interface bsts_aggregate_time_series
   public :: bsts_dynamic_regression_draws, bsts_dynamic_regression
   public :: bsts_dynamic_regression_hierarchical_draws
   public :: bsts_dynamic_regression_hierarchical
   public :: bsts_dynamic_regression_predict_draws
   public :: bsts_dynamic_regression_predict
   public :: bsts_dynamic_regression_ar_draws, bsts_dynamic_regression_ar
   public :: bsts_dynamic_regression_ar_predict_draws
   public :: bsts_dynamic_regression_ar_predict
   public :: bsts_predict_draws, bsts_predict
   public :: bsts_geometric_sequence, bsts_harvey_cumulator
   public :: bsts_ar_draws, bsts_ar, bsts_ar_predict_draws, bsts_ar_predict
   public :: bsts_auto_ar_draws, bsts_auto_ar
   public :: bsts_student_local_linear_trend_draws
   public :: bsts_student_local_linear_trend
   public :: bsts_student_trend_predict_draws, bsts_student_trend_predict
   public :: bsts_monthly_annual_cycle_draws, bsts_monthly_annual_cycle
   public :: bsts_monthly_predict_draws, bsts_monthly_predict
   public :: bsts_fixed_date_holiday, bsts_nth_weekday_holiday
   public :: bsts_last_weekday_holiday, bsts_date_range_holiday
   public :: bsts_named_holiday, bsts_holiday_position, bsts_holiday_width
   public :: bsts_random_walk_holiday_draws, bsts_random_walk_holiday
   public :: bsts_holiday_predict_draws, bsts_holiday_predict
   public :: bsts_holiday_design
   public :: bsts_regression_holiday_draws, bsts_regression_holiday
   public :: bsts_hierarchical_regression_holiday_draws
   public :: bsts_hierarchical_regression_holiday
   public :: bsts_regression_holiday_predict_draws
   public :: bsts_regression_holiday_predict
   public :: bsts_shared_local_level_draws, bsts_shared_local_level
   public :: bsts_shared_local_level_predict_draws
   public :: bsts_shared_local_level_predict
   public :: bsts_mbsts_draws, bsts_mbsts
   public :: bsts_mbsts_predict_draws, bsts_mbsts_predict
   public :: bsts_logit_local_level_draws, bsts_logit_local_level
   public :: bsts_poisson_local_level_draws, bsts_poisson_local_level
   public :: bsts_logit_predict_draws, bsts_logit_predict
   public :: bsts_poisson_predict_draws, bsts_poisson_predict
   public :: bsts_logit_regression_draws, bsts_logit_regression
   public :: bsts_poisson_regression_draws, bsts_poisson_regression
   public :: bsts_logit_regression_predict_draws
   public :: bsts_poisson_regression_predict_draws
   public :: bsts_logit_trend_seasonal_draws
   public :: bsts_poisson_trend_seasonal_draws
   public :: bsts_logit_trend_seasonal, bsts_poisson_trend_seasonal
   public :: bsts_logit_trend_seasonal_predict_draws
   public :: bsts_poisson_trend_seasonal_predict_draws

contains

   pure logical function bsts_no_duplicates_numeric(timestamps) result(valid)
      !! Report whether numeric timestamps contain no repeated values.
      real(dp), intent(in) :: timestamps(:) !! Timestamps.
      integer :: left, right

      valid = all(ieee_is_finite(timestamps))
      if (.not. valid) return
      do left = 1, size(timestamps) - 1
         do right = left + 1, size(timestamps)
            if (timestamps(left) == timestamps(right)) then
               valid = .false.
               return
            end if
         end do
      end do
   end function bsts_no_duplicates_numeric

   pure logical function bsts_no_duplicates_date(timestamps) result(valid)
      !! Report whether date timestamps contain no repeated values.
      type(date_t), intent(in) :: timestamps(:) !! Timestamps.
      integer :: left, right

      valid = all(date_valid(timestamps))
      if (.not. valid) return
      do left = 1, size(timestamps) - 1
         do right = left + 1, size(timestamps)
            if (timestamps(left) == timestamps(right)) then
               valid = .false.
               return
            end if
         end do
      end do
   end function bsts_no_duplicates_date

   pure logical function bsts_no_gaps_numeric(timestamps) result(valid)
      !! Detect numeric gaps at least 1.8 times the smallest increment.
      real(dp), intent(in) :: timestamps(:) !! Timestamps.
      real(dp), allocatable :: unique(:)
      real(dp) :: minimum

      if (.not. all(ieee_is_finite(timestamps))) then
         valid = .false.
         return
      end if
      unique = sorted_unique_real(timestamps)
      if (size(unique) < 3) then
         valid = size(unique) > 0
         return
      end if
      minimum = minval(unique(2:) - unique(:size(unique) - 1))
      valid = minimum > 0.0_dp .and. all(unique(2:) - &
         unique(:size(unique) - 1) < 1.8_dp*minimum)
   end function bsts_no_gaps_numeric

   pure logical function bsts_no_gaps_date(timestamps) result(valid)
      !! Detect date gaps at least 1.8 times the smallest increment.
      type(date_t), intent(in) :: timestamps(:) !! Timestamps.
      integer, allocatable :: days(:)
      real(dp), allocatable :: numeric(:)
      integer :: index

      if (.not. all(date_valid(timestamps))) then
         valid = .false.
         return
      end if
      allocate(days(size(timestamps)), numeric(size(timestamps)))
      do index = 1, size(timestamps)
         days(index) = date_day_number(timestamps(index))
      end do
      numeric = real(days, dp)
      valid = bsts_no_gaps_numeric(numeric)
   end function bsts_no_gaps_date

   pure logical function bsts_is_regular_numeric(timestamps) result(valid)
      !! Report whether numeric timestamps have no duplicates or gaps.
      real(dp), intent(in) :: timestamps(:) !! Timestamps.

      valid = bsts_no_duplicates_numeric(timestamps) .and. &
         bsts_no_gaps_numeric(timestamps)
   end function bsts_is_regular_numeric

   pure logical function bsts_is_regular_date(timestamps) result(valid)
      !! Report whether date timestamps have no duplicates or gaps.
      type(date_t), intent(in) :: timestamps(:) !! Timestamps.

      valid = bsts_no_duplicates_date(timestamps) .and. &
         bsts_no_gaps_date(timestamps)
   end function bsts_is_regular_date

   pure function bsts_regularize_numeric_timestamps(timestamps) result(out)
      !! Expand numeric timestamps to their smallest regular grid.
      real(dp), intent(in) :: timestamps(:) !! Timestamps.
      type(bsts_numeric_timestamps_t) :: out
      real(dp), allocatable :: unique(:)
      real(dp) :: increment, position, tolerance
      integer :: points, observation, index

      if (size(timestamps) < 1 .or. .not. all(ieee_is_finite(timestamps))) then
         out%info = 1
         return
      end if
      unique = sorted_unique_real(timestamps)
      if (size(unique) == 1) then
         out%grid = unique
         allocate(out%mapping(size(timestamps)))
         out%mapping = 1
         out%trivial = size(timestamps) == 1
         return
      end if
      increment = minval(unique(2:) - unique(:size(unique) - 1))
      points = nint((unique(size(unique)) - unique(1))/increment) + 1
      allocate(out%grid(points), out%mapping(size(timestamps)))
      out%grid = [(unique(1) + increment*real(index - 1, dp), index=1, points)]
      tolerance = 1.0e-8_dp*max(1.0_dp, maxval(abs(timestamps)))
      do observation = 1, size(timestamps)
         position = (timestamps(observation) - unique(1))/increment
         index = nint(position) + 1
         if (index < 1 .or. index > points .or. &
            abs(timestamps(observation) - out%grid(index)) > tolerance) then
            out%info = 2
            return
         end if
         out%mapping(observation) = index
      end do
      out%trivial = points == size(timestamps) .and. &
         all(out%mapping == [(index, index=1, size(timestamps))])
   end function bsts_regularize_numeric_timestamps

   pure function bsts_regularize_date_timestamps(timestamps) result(out)
      !! Expand Gregorian dates on a daily, weekly, monthly, quarterly, or yearly grid.
      type(date_t), intent(in) :: timestamps(:) !! Timestamps.
      type(bsts_date_timestamps_t) :: out
      integer, allocatable :: unique_days(:), deltas(:)
      type(date_t) :: first, candidate
      integer :: minimum, points, index, observation, scale, step

      if (size(timestamps) < 1 .or. .not. all(date_valid(timestamps))) then
         out%info = 1
         return
      end if
      unique_days = sorted_unique_integer(date_day_number(timestamps))
      if (size(unique_days) == 1) then
         allocate(out%grid(1), out%mapping(size(timestamps)))
         out%grid(1) = date_from_day_number(unique_days(1))
         out%mapping = 1
         out%trivial = size(timestamps) == 1
         return
      end if
      deltas = unique_days(2:) - unique_days(:size(unique_days) - 1)
      minimum = minval(deltas)
      scale = 0
      step = 0
      if (minimum == 1 .or. minimum == 7) then
         scale = 1
         step = minimum
      else if (minimum >= 28 .and. minimum <= 31 .and. &
         any(deltas == 30 .or. deltas == 31)) then
         scale = 2
      else if (minimum >= 89 .and. minimum <= 93) then
         scale = 3
      else if (minimum >= 365 .and. minimum <= 366) then
         scale = 4
      else
         out%info = 2
         return
      end if
      first = date_from_day_number(unique_days(1))
      if (scale == 1) then
         points = (unique_days(size(unique_days)) - unique_days(1))/step + 1
      else
         points = 1
         candidate = first
         do while (candidate < date_from_day_number(unique_days(size(unique_days))))
            points = points + 1
            candidate = add_calendar_months(first, &
               merge(points - 1, merge(3*(points - 1), 12*(points - 1), &
               scale == 3), scale == 2))
         end do
      end if
      allocate(out%grid(points), out%mapping(size(timestamps)))
      do index = 1, points
         if (scale == 1) then
            out%grid(index) = first + step*(index - 1)
         else if (scale == 2) then
            out%grid(index) = add_calendar_months(first, index - 1)
         else if (scale == 3) then
            out%grid(index) = add_calendar_months(first, 3*(index - 1))
         else
            out%grid(index) = add_calendar_months(first, 12*(index - 1))
         end if
      end do
      do observation = 1, size(timestamps)
         out%mapping(observation) = 0
         do index = 1, points
            if (timestamps(observation) == out%grid(index)) then
               out%mapping(observation) = index
               exit
            end if
         end do
         if (out%mapping(observation) == 0) then
            out%info = 3
            return
         end if
      end do
      out%trivial = points == size(timestamps) .and. &
         all(out%mapping == [(index, index=1, size(timestamps))])
   end function bsts_regularize_date_timestamps

   pure function bsts_long_to_wide(response, series_id, timestamps) result(out)
      !! Reshape integer-labelled long observations into a wide matrix.
      real(dp), intent(in) :: response(:) !! Response observations.
      integer, intent(in) :: series_id(:) !! Series identifier.
      integer, intent(in) :: timestamps(:) !! Timestamps.
      type(bsts_wide_series_t) :: out
      logical, allocatable :: seen(:, :)
      integer :: observation, time, series

      if (size(response) /= size(series_id) .or. &
         size(response) /= size(timestamps)) then
         out%info = 1
         return
      end if
      out%timestamps = sorted_unique_integer(timestamps)
      out%series_id = sorted_unique_integer(series_id)
      allocate(out%values(size(out%timestamps), size(out%series_id)))
      allocate(seen(size(out%timestamps), size(out%series_id)))
      out%values = ieee_value(0.0_dp, ieee_quiet_nan)
      seen = .false.
      do observation = 1, size(response)
         time = integer_position(out%timestamps, timestamps(observation))
         series = integer_position(out%series_id, series_id(observation))
         if (seen(time, series)) then
            out%info = 2
            return
         end if
         out%values(time, series) = response(observation)
         seen(time, series) = .true.
      end do
   end function bsts_long_to_wide

   pure function bsts_wide_to_long(values, timestamps, series_id, &
      remove_missing) result(out)
      !! Reshape a wide matrix into integer-labelled long observations.
      real(dp), intent(in) :: values(:, :) !! Input values.
      integer, intent(in) :: timestamps(:) !! Timestamps.
      integer, intent(in) :: series_id(:) !! Series identifier.
      logical, intent(in), optional :: remove_missing !! Flag controlling remove missing.
      type(bsts_long_series_t) :: out
      logical :: omit
      integer :: total, time, series, destination

      if (size(timestamps) /= size(values, 1) .or. &
         size(series_id) /= size(values, 2)) then
         out%info = 1
         return
      end if
      omit = .true.
      if (present(remove_missing)) omit = remove_missing
      total = size(values)
      if (omit) total = count(ieee_is_finite(values))
      allocate(out%values(total), out%timestamps(total), out%series_id(total))
      destination = 0
      do time = 1, size(values, 1)
         do series = 1, size(values, 2)
            if (omit .and. .not. ieee_is_finite(values(time, series))) cycle
            destination = destination + 1
            out%values(destination) = values(time, series)
            out%timestamps(destination) = timestamps(time)
            out%series_id(destination) = series_id(series)
         end do
      end do
   end function bsts_wide_to_long

   pure function aggregate_time_series_vector(fine_series, contains_end, &
      membership_fraction, trim_left, trim_right) result(aggregate)
      !! Aggregate a fine-scale vector using Harvey boundary fractions.
      real(dp), intent(in) :: fine_series(:) !! Fine series.
      real(dp), intent(in) :: membership_fraction(:) !! Membership fraction.
      logical, intent(in) :: contains_end(:) !! Flag controlling contains end.
      logical, intent(in), optional :: trim_left !! Flag controlling trim left.
      logical, intent(in), optional :: trim_right !! Flag controlling trim right.
      real(dp), allocatable :: aggregate(:), work(:)
      real(dp) :: accumulated
      logical :: left, right, has_remainder
      integer :: time, destination, first, last, maximum

      if (size(fine_series) < 1 .or. &
         size(contains_end) /= size(fine_series) .or. &
         size(membership_fraction) /= size(fine_series) .or. &
         any(membership_fraction <= 0.0_dp) .or. &
         any(membership_fraction > 1.0_dp)) then
         allocate(aggregate(0))
         return
      end if
      left = any(membership_fraction < 1.0_dp)
      if (present(trim_left)) left = trim_left
      has_remainder = .not. contains_end(size(contains_end)) .or. &
         membership_fraction(size(membership_fraction)) < 0.9999_dp
      right = has_remainder
      if (present(trim_right)) right = trim_right
      maximum = count(contains_end) + merge(1, 0, has_remainder)
      allocate(work(maximum))
      accumulated = 0.0_dp
      destination = 0
      do time = 1, size(fine_series)
         if (contains_end(time)) then
            destination = destination + 1
            work(destination) = accumulated + &
               membership_fraction(time)*fine_series(time)
            accumulated = (1.0_dp - membership_fraction(time))*fine_series(time)
         else
            accumulated = accumulated + fine_series(time)
         end if
      end do
      if (has_remainder) then
         destination = destination + 1
         work(destination) = accumulated
      end if
      first = merge(2, 1, left)
      last = destination - merge(1, 0, right)
      if (last < first) then
         allocate(aggregate(0))
      else
         aggregate = work(first:last)
      end if
   end function aggregate_time_series_vector

   pure function aggregate_time_series_matrix(fine_series, contains_end, &
      membership_fraction, trim_left, trim_right) result(aggregate)
      !! Aggregate each column of a fine-scale time-by-series matrix.
      real(dp), intent(in) :: fine_series(:, :) !! Fine series.
      real(dp), intent(in) :: membership_fraction(:) !! Membership fraction.
      logical, intent(in) :: contains_end(:) !! Flag controlling contains end.
      logical, intent(in), optional :: trim_left !! Flag controlling trim left.
      logical, intent(in), optional :: trim_right !! Flag controlling trim right.
      real(dp), allocatable :: aggregate(:, :), column(:)
      integer :: series

      if (size(fine_series, 2) < 1) then
         allocate(aggregate(0, 0))
         return
      end if
      if (present(trim_left) .and. present(trim_right)) then
         column = aggregate_time_series_vector(fine_series(:, 1), contains_end, &
            membership_fraction, trim_left, trim_right)
      else if (present(trim_left)) then
         column = aggregate_time_series_vector(fine_series(:, 1), contains_end, &
            membership_fraction, trim_left=trim_left)
      else if (present(trim_right)) then
         column = aggregate_time_series_vector(fine_series(:, 1), contains_end, &
            membership_fraction, trim_right=trim_right)
      else
         column = aggregate_time_series_vector(fine_series(:, 1), contains_end, &
            membership_fraction)
      end if
      allocate(aggregate(size(column), size(fine_series, 2)))
      aggregate(:, 1) = column
      do series = 2, size(fine_series, 2)
         if (present(trim_left) .and. present(trim_right)) then
            column = aggregate_time_series_vector(fine_series(:, series), &
               contains_end, membership_fraction, trim_left, trim_right)
         else if (present(trim_left)) then
            column = aggregate_time_series_vector(fine_series(:, series), &
               contains_end, membership_fraction, trim_left=trim_left)
         else if (present(trim_right)) then
            column = aggregate_time_series_vector(fine_series(:, series), &
               contains_end, membership_fraction, trim_right=trim_right)
         else
            column = aggregate_time_series_vector(fine_series(:, series), &
               contains_end, membership_fraction)
         end if
         aggregate(:, series) = column
      end do
   end function aggregate_time_series_matrix

   pure elemental integer function bsts_month_distance(date, origin) &
      result(distance)
      !! Return the signed calendar-month distance from an origin date.
      type(date_t), intent(in) :: date !! Date.
      type(date_t), intent(in) :: origin !! Origin.

      distance = 12*(date%year - origin%year) + date%month - origin%month
   end function bsts_month_distance

   pure elemental logical function bsts_week_ends_month(week_ending) &
      result(ends_month)
      !! Report whether a seven-day interval contains a month end.
      type(date_t), intent(in) :: week_ending !! Week ending.
      type(date_t) :: first_day, following_day

      first_day = week_ending - 6
      following_day = week_ending + 1
      ends_month = first_day%year /= following_day%year .or. &
         first_day%month /= following_day%month
   end function bsts_week_ends_month

   pure elemental logical function bsts_week_ends_quarter(week_ending) &
      result(ends_quarter)
      !! Report whether a seven-day interval contains a quarter end.
      type(date_t), intent(in) :: week_ending !! Week ending.
      type(date_t) :: first_day, following_day
      integer :: first_quarter, following_quarter

      first_day = week_ending - 6
      following_day = week_ending + 1
      first_quarter = 4*first_day%year + (first_day%month - 1)/3
      following_quarter = 4*following_day%year + (following_day%month - 1)/3
      ends_quarter = first_quarter /= following_quarter
   end function bsts_week_ends_quarter

   pure elemental real(dp) function bsts_fraction_initial_month(week_ending) &
      result(fraction)
      !! Return the fraction of a week in the month containing its first day.
      type(date_t), intent(in) :: week_ending !! Week ending.
      type(date_t) :: first_day

      first_day = week_ending - 6
      if (first_day%year == week_ending%year .and. &
         first_day%month == week_ending%month) then
         fraction = 1.0_dp
      else
         fraction = 1.0_dp - real(week_ending%day, dp)/7.0_dp
      end if
   end function bsts_fraction_initial_month

   pure function bsts_match_week_to_month(week_ending, origin_month) &
      result(index)
      !! Map week-ending dates to one-based months containing each first day.
      type(date_t), intent(in) :: week_ending(:) !! Week ending.
      type(date_t), intent(in) :: origin_month !! Origin month.
      integer :: index(size(week_ending))
      integer :: week

      do week = 1, size(week_ending)
         index(week) = 1 + bsts_month_distance(week_ending(week) - 6, &
            origin_month)
      end do
   end function bsts_match_week_to_month

   pure function bsts_aggregate_weeks_to_months(weekly, week_ending, &
      membership_fraction, trim_left, trim_right) result(out)
      !! Aggregate weekly time-by-series observations to calendar months.
      real(dp), intent(in) :: weekly(:, :) !! Weekly.
      type(date_t), intent(in) :: week_ending(:) !! Week ending.
      real(dp), intent(in), optional :: membership_fraction(:) !! Membership fraction.
      logical, intent(in), optional :: trim_left !! Flag controlling trim left.
      logical, intent(in), optional :: trim_right !! Flag controlling trim right.
      type(bsts_monthly_series_t) :: out
      real(dp), allocatable :: fraction(:)
      logical, allocatable :: ends(:)
      type(date_t), allocatable :: labels(:), work_labels(:)
      logical :: left, right, has_remainder
      integer :: week, destination, first, last

      if (size(weekly, 1) < 1 .or. size(weekly, 1) /= size(week_ending) .or. &
         .not. all(date_valid(week_ending))) then
         out%info = 1
         return
      end if
      allocate(fraction(size(week_ending)), ends(size(week_ending)))
      fraction = bsts_fraction_initial_month(week_ending)
      if (present(membership_fraction)) then
         if (size(membership_fraction) /= size(week_ending)) then
            out%info = 1
            return
         end if
         fraction = membership_fraction
      end if
      ends = bsts_week_ends_month(week_ending)
      left = .true.
      if (present(trim_left)) left = trim_left
      has_remainder = .not. ends(size(ends)) .or. &
         fraction(size(fraction)) < 0.9999_dp
      right = has_remainder
      if (present(trim_right)) right = trim_right
      out%values = aggregate_time_series_matrix(weekly, ends, fraction, &
         left, right)
      allocate(work_labels(count(ends) + merge(1, 0, has_remainder)))
      destination = 0
      do week = 1, size(week_ending)
         if (.not. ends(week)) cycle
         destination = destination + 1
         work_labels(destination) = last_day_in_month(week_ending(week) - 6)
      end do
      if (has_remainder) then
         destination = destination + 1
         work_labels(destination) = last_day_in_month(week_ending(size(week_ending)))
      end if
      first = merge(2, 1, left)
      last = destination - merge(1, 0, right)
      if (last < first) then
         allocate(out%dates(0))
      else
         labels = work_labels(first:last)
         out%dates = labels
      end if
   end function bsts_aggregate_weeks_to_months

   pure function bsts_fixed_date_holiday(month, day, days_before, days_after, &
      name) result(holiday)
      !! Define a holiday occurring on a fixed month and day each year.
      integer, intent(in) :: month !! Month.
      integer, intent(in) :: day !! Day.
      integer, intent(in) :: days_before !! Days before.
      integer, intent(in) :: days_after !! Days after.
      character(*), intent(in), optional :: name !! Name.
      type(bsts_holiday_t) :: holiday

      holiday%kind = bsts_holiday_fixed
      holiday%month = month
      holiday%day = day
      holiday%days_before = days_before
      holiday%days_after = days_after
      if (present(name)) holiday%name = name
      if (days_before < 0 .or. days_after < 0 .or. &
         .not. date_valid(date_t(2000, month, day))) holiday%info = 1
   end function bsts_fixed_date_holiday

   pure function bsts_nth_weekday_holiday(month, weekday, week_number, &
      days_before, days_after, name) result(holiday)
      !! Define the nth ISO weekday in a month as a holiday anchor.
      integer, intent(in) :: month !! Month.
      integer, intent(in) :: weekday !! Weekday.
      integer, intent(in) :: week_number !! Week number.
      integer, intent(in) :: days_before !! Days before.
      integer, intent(in) :: days_after !! Days after.
      character(*), intent(in), optional :: name !! Name.
      type(bsts_holiday_t) :: holiday
      type(date_t) :: anchor

      holiday%kind = bsts_holiday_nth_weekday
      holiday%month = month
      holiday%weekday = weekday
      holiday%week_number = week_number
      holiday%days_before = days_before
      holiday%days_after = days_after
      if (present(name)) holiday%name = name
      anchor = holiday_anchor_date(holiday, 2024)
      if (month < 1 .or. month > 12 .or. weekday < 1 .or. weekday > 7 .or. &
         week_number < 1 .or. days_before < 0 .or. days_after < 0 .or. &
         .not. date_valid(anchor)) holiday%info = 1
   end function bsts_nth_weekday_holiday

   pure function bsts_last_weekday_holiday(month, weekday, days_before, &
      days_after, name) result(holiday)
      !! Define the final ISO weekday in a month as a holiday anchor.
      integer, intent(in) :: month !! Month.
      integer, intent(in) :: weekday !! Weekday.
      integer, intent(in) :: days_before !! Days before.
      integer, intent(in) :: days_after !! Days after.
      character(*), intent(in), optional :: name !! Name.
      type(bsts_holiday_t) :: holiday

      holiday%kind = bsts_holiday_last_weekday
      holiday%month = month
      holiday%weekday = weekday
      holiday%days_before = days_before
      holiday%days_after = days_after
      if (present(name)) holiday%name = name
      if (month < 1 .or. month > 12 .or. weekday < 1 .or. weekday > 7 .or. &
         days_before < 0 .or. days_after < 0) holiday%info = 1
   end function bsts_last_weekday_holiday

   pure function bsts_date_range_holiday(range_start, range_end, name) &
      result(holiday)
      !! Define irregular nonoverlapping holiday influence intervals.
      type(date_t), intent(in) :: range_start(:) !! Range start.
      type(date_t), intent(in) :: range_end(:) !! Range end.
      character(*), intent(in), optional :: name !! Name.
      type(bsts_holiday_t) :: holiday
      integer :: occurrence

      holiday%kind = bsts_holiday_date_range
      if (present(name)) holiday%name = name
      if (size(range_start) < 1 .or. size(range_end) /= size(range_start)) then
         holiday%info = 1
         return
      end if
      allocate(holiday%range_start(size(range_start)), &
         holiday%range_end(size(range_end)))
      holiday%range_start = range_start
      holiday%range_end = range_end
      do occurrence = 1, size(range_start)
         if (.not. date_valid(range_start(occurrence)) .or. &
            .not. date_valid(range_end(occurrence)) .or. &
            date_day_number(range_end(occurrence)) < &
            date_day_number(range_start(occurrence))) holiday%info = 1
         if (occurrence > 1) then
            if (date_day_number(range_start(occurrence)) <= &
               date_day_number(range_end(occurrence - 1))) holiday%info = 1
         end if
      end do
   end function bsts_date_range_holiday

   pure function bsts_named_holiday(name, days_before, days_after) &
      result(holiday)
      !! Define one of the recurring named holidays supported by bsts.
      character(*), intent(in) :: name !! Name.
      integer, intent(in) :: days_before !! Days before.
      integer, intent(in) :: days_after !! Days after.
      type(bsts_holiday_t) :: holiday

      holiday%kind = bsts_holiday_named
      holiday%name = name
      holiday%days_before = days_before
      holiday%days_after = days_after
      if (days_before < 0 .or. days_after < 0 .or. &
         .not. date_valid(holiday_anchor_date(holiday, 2024))) holiday%info = 1
   end function bsts_named_holiday

   pure function bsts_holiday_width(holiday) result(width)
      !! Return the maximum number of relative days in a holiday window.
      type(bsts_holiday_t), intent(in) :: holiday !! Holiday.
      integer :: width, occurrence

      width = 0
      if (holiday%info /= 0) return
      if (holiday%kind == bsts_holiday_date_range) then
         if (.not. allocated(holiday%range_start) .or. &
            .not. allocated(holiday%range_end)) return
         do occurrence = 1, size(holiday%range_start)
            width = max(width, 1 + date_day_number(holiday%range_end(occurrence)) - &
               date_day_number(holiday%range_start(occurrence)))
         end do
      else
         width = 1 + holiday%days_before + holiday%days_after
      end if
   end function bsts_holiday_width

   pure elemental function bsts_holiday_position(holiday, date) result(position)
      !! Return the one-based relative holiday day, or zero outside its window.
      type(bsts_holiday_t), intent(in) :: holiday !! Holiday.
      type(date_t), intent(in) :: date !! Date.
      integer :: position
      type(date_t) :: anchor
      integer :: year, difference, occurrence

      position = 0
      if (holiday%info /= 0 .or. .not. date_valid(date)) return
      if (holiday%kind == bsts_holiday_date_range) then
         if (.not. allocated(holiday%range_start) .or. &
            .not. allocated(holiday%range_end)) return
         do occurrence = 1, size(holiday%range_start)
            if (date_day_number(date) >= &
               date_day_number(holiday%range_start(occurrence)) .and. &
               date_day_number(date) <= &
               date_day_number(holiday%range_end(occurrence))) then
               position = 1 + date_day_number(date) - &
                  date_day_number(holiday%range_start(occurrence))
               return
            end if
         end do
         return
      end if
      do year = date%year - 1, date%year + 1
         anchor = holiday_anchor_date(holiday, year)
         if (.not. date_valid(anchor)) cycle
         difference = date_day_number(date) - date_day_number(anchor)
         if (difference >= -holiday%days_before .and. &
            difference <= holiday%days_after) then
            position = difference + holiday%days_before + 1
            return
         end if
      end do
   end function bsts_holiday_position

   pure function bsts_holiday_design(first_date, observations, holidays) &
      result(design)
      !! Build a sparse daily design matrix for concatenated holiday windows.
      type(date_t), intent(in) :: first_date !! First date.
      integer, intent(in) :: observations !! Observed time-series values.
      type(bsts_holiday_t), intent(in) :: holidays(:) !! Holidays.
      real(dp), allocatable :: design(:, :)
      type(date_t) :: date
      integer, allocatable :: offset(:)
      integer :: total, holiday_index, time, position

      total = 0
      if (observations < 0 .or. .not. date_valid(first_date)) then
         allocate(design(0, 0))
         return
      end if
      allocate(offset(size(holidays) + 1))
      offset(1) = 0
      do holiday_index = 1, size(holidays)
         if (holidays(holiday_index)%info /= 0 .or. &
            bsts_holiday_width(holidays(holiday_index)) < 1) then
            allocate(design(0, 0))
            return
         end if
         total = total + bsts_holiday_width(holidays(holiday_index))
         offset(holiday_index + 1) = total
      end do
      allocate(design(observations, total))
      design = 0.0_dp
      date = first_date
      do time = 1, observations
         do holiday_index = 1, size(holidays)
            position = bsts_holiday_position(holidays(holiday_index), date)
            if (position > 0) design(time, offset(holiday_index) + position) = &
               1.0_dp
         end do
         date = date_from_day_number(date_day_number(date) + 1)
      end do
   end function bsts_holiday_design

   pure function holiday_anchor_date(holiday, year) result(anchor)
      !! Resolve a recurring holiday definition to its anchor date in one year.
      type(bsts_holiday_t), intent(in) :: holiday !! Holiday.
      integer, intent(in) :: year !! Year.
      type(date_t) :: anchor, first, last
      integer :: day, offset

      anchor = date_t(0, 0, 0)
      select case (holiday%kind)
      case (bsts_holiday_fixed)
         anchor = date_t(year, holiday%month, holiday%day)
      case (bsts_holiday_nth_weekday)
         first = date_t(year, holiday%month, 1)
         if (.not. date_valid(first)) return
         offset = modulo(holiday%weekday - date_day_of_week(first), 7)
         day = 1 + offset + 7*(holiday%week_number - 1)
         anchor = date_t(year, holiday%month, day)
      case (bsts_holiday_last_weekday)
         day = date_days_in_month(year, holiday%month)
         last = date_t(year, holiday%month, day)
         if (.not. date_valid(last)) return
         offset = modulo(date_day_of_week(last) - holiday%weekday, 7)
         anchor = date_t(year, holiday%month, day - offset)
      case (bsts_holiday_named)
         anchor = named_holiday_anchor(trim(holiday%name), year)
      end select
      if (.not. date_valid(anchor)) anchor = date_t(0, 0, 0)
   end function holiday_anchor_date

   pure function named_holiday_anchor(name, year) result(anchor)
      !! Resolve a canonical bsts named holiday for one calendar year.
      character(*), intent(in) :: name !! Name.
      integer, intent(in) :: year !! Year.
      type(date_t) :: anchor
      type(bsts_holiday_t) :: rule

      anchor = date_t(0, 0, 0)
      select case (name)
      case ('NewYearsDay')
         anchor = date_t(year, 1, 1)
      case ('SuperBowlSunday')
         rule = bsts_nth_weekday_holiday(2, 7, 2, 0, 0)
         anchor = holiday_anchor_date(rule, year)
      case ('MartinLutherKingDay')
         rule = bsts_nth_weekday_holiday(1, 1, 3, 0, 0)
         anchor = holiday_anchor_date(rule, year)
      case ('PresidentsDay')
         rule = bsts_nth_weekday_holiday(2, 1, 3, 0, 0)
         anchor = holiday_anchor_date(rule, year)
      case ('ValentinesDay')
         anchor = date_t(year, 2, 14)
      case ('SaintPatricksDay')
         anchor = date_t(year, 3, 17)
      case ('USDaylightSavingsTimeBegins')
         rule = bsts_nth_weekday_holiday(3, 7, 2, 0, 0)
         anchor = holiday_anchor_date(rule, year)
      case ('USDaylightSavingsTimeEnds')
         rule = bsts_nth_weekday_holiday(11, 7, 1, 0, 0)
         anchor = holiday_anchor_date(rule, year)
      case ('EasterSunday')
         anchor = date_easter(year)
      case ('USMothersDay')
         rule = bsts_nth_weekday_holiday(5, 7, 2, 0, 0)
         anchor = holiday_anchor_date(rule, year)
      case ('IndependenceDay')
         anchor = date_t(year, 7, 4)
      case ('LaborDay')
         rule = bsts_nth_weekday_holiday(9, 1, 1, 0, 0)
         anchor = holiday_anchor_date(rule, year)
      case ('ColumbusDay')
         rule = bsts_nth_weekday_holiday(10, 1, 2, 0, 0)
         anchor = holiday_anchor_date(rule, year)
      case ('Halloween')
         anchor = date_t(year, 10, 31)
      case ('Thanksgiving')
         rule = bsts_nth_weekday_holiday(11, 4, 4, 0, 0)
         anchor = holiday_anchor_date(rule, year)
      case ('MemorialDay')
         rule = bsts_last_weekday_holiday(5, 1, 0, 0)
         anchor = holiday_anchor_date(rule, year)
      case ('VeteransDay')
         anchor = date_t(year, 11, 11)
      case ('Christmas')
         anchor = date_t(year, 12, 25)
      end select
   end function named_holiday_anchor

   pure function bsts_local_level_draws(y, initial_mean, initial_variance, &
      observation_variance, level_variance, observation_prior_shape, &
      observation_prior_rate, level_prior_shape, level_prior_rate, burn, &
      state_normal_draws, gamma_draws) result(out)
      !! Draw a local-level posterior using supplied normal and gamma variates.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean !! Initial state mean.
      real(dp), intent(in) :: initial_variance !! Initial variance.
      real(dp), intent(in) :: observation_variance !! Observation-error variance.
      real(dp), intent(in) :: level_variance !! Level variance.
      real(dp), intent(in) :: observation_prior_shape !! Observation prior shape.
      real(dp), intent(in) :: observation_prior_rate !! Observation prior rate.
      real(dp), intent(in) :: level_prior_shape !! Level prior shape.
      real(dp), intent(in) :: level_prior_rate !! Level prior rate.
      integer, intent(in) :: burn !! Number of initial simulation draws to discard.
      real(dp), intent(in) :: state_normal_draws(:, :, :) !! Independent standard-normal state draws.
      real(dp), intent(in) :: gamma_draws(:, :) !! Gamma simulation draws.
      type(bsts_mcmc_t) :: out
      real(dp) :: transition(1, 1), observation(1)
      real(dp) :: state_variances(1), state_shapes(1), state_rates(1)

      transition = 1.0_dp
      observation = 1.0_dp
      state_variances = level_variance
      state_shapes = level_prior_shape
      state_rates = level_prior_rate
      out = structural_gibbs_draws(y, [initial_mean], &
         reshape([initial_variance], [1, 1]), transition, observation, &
         observation_variance, state_variances, observation_prior_shape, &
         observation_prior_rate, state_shapes, state_rates, burn, &
         state_normal_draws, gamma_draws)
   end function bsts_local_level_draws

   function bsts_local_level(y, iterations, burn, initial_mean, &
      initial_variance, observation_variance, level_variance, &
      observation_prior_shape, observation_prior_rate, level_prior_shape, &
      level_prior_rate) result(out)
      !! Sample the local-level posterior using the shared random stream.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      integer, intent(in) :: iterations !! Number of algorithm iterations.
      integer, intent(in), optional :: burn !! Number of initial simulation draws to discard.
      real(dp), intent(in), optional :: initial_mean !! Initial state mean.
      real(dp), intent(in), optional :: initial_variance !! Initial variance.
      real(dp), intent(in), optional :: observation_variance !! Observation-error variance.
      real(dp), intent(in), optional :: level_variance !! Level variance.
      real(dp), intent(in), optional :: observation_prior_shape !! Observation prior shape.
      real(dp), intent(in), optional :: observation_prior_rate !! Observation prior rate.
      real(dp), intent(in), optional :: level_prior_shape !! Level prior shape.
      real(dp), intent(in), optional :: level_prior_rate !! Level prior rate.
      type(bsts_mcmc_t) :: out
      real(dp), allocatable :: normal_draws(:, :, :), gamma_draws(:, :)
      real(dp) :: mean0, variance0, observation0, level0
      real(dp) :: observation_shape, observation_rate
      real(dp) :: level_shape, level_rate, series_variance
      integer :: discarded, iteration, observed

      if (size(y) < 2 .or. iterations < 1) then
         out%info = 1
         return
      end if
      discarded = 0
      if (present(burn)) discarded = burn
      series_variance = finite_sample_variance(y)
      mean0 = first_finite_value(y)
      if (present(initial_mean)) mean0 = initial_mean
      variance0 = max(series_variance, 1.0_dp)
      if (present(initial_variance)) variance0 = initial_variance
      observation0 = max(0.5_dp*series_variance, 1.0e-6_dp)
      if (present(observation_variance)) observation0 = observation_variance
      level0 = max(0.01_dp*series_variance, 1.0e-8_dp)
      if (present(level_variance)) level0 = level_variance
      observation_shape = 0.01_dp
      if (present(observation_prior_shape)) &
         observation_shape = observation_prior_shape
      observation_rate = max(0.01_dp*series_variance, 1.0e-8_dp)
      if (present(observation_prior_rate)) observation_rate = &
         observation_prior_rate
      level_shape = 0.01_dp
      if (present(level_prior_shape)) level_shape = level_prior_shape
      level_rate = max(0.0001_dp*series_variance, 1.0e-10_dp)
      if (present(level_prior_rate)) level_rate = level_prior_rate
      if (discarded < 0 .or. discarded >= iterations .or. &
         variance0 <= 0.0_dp .or. observation0 <= 0.0_dp .or. &
         level0 <= 0.0_dp .or. observation_shape <= 0.0_dp .or. &
         observation_rate <= 0.0_dp .or. level_shape <= 0.0_dp .or. &
         level_rate <= 0.0_dp) then
         out%info = 1
         return
      end if
      allocate(normal_draws(1, size(y), iterations), &
         gamma_draws(2, iterations))
      observed = count(ieee_is_finite(y))
      do iteration = 1, iterations
         call fill_standard_normals(normal_draws(:, :, iteration))
         gamma_draws(1, iteration) = random_gamma( &
            observation_shape + 0.5_dp*real(observed, dp))
         gamma_draws(2, iteration) = random_gamma( &
            level_shape + 0.5_dp*real(size(y) - 1, dp))
      end do
      out = bsts_local_level_draws(y, mean0, variance0, observation0, &
         level0, observation_shape, observation_rate, level_shape, &
         level_rate, discarded, normal_draws, gamma_draws)
   end function bsts_local_level

   pure function bsts_local_linear_trend_draws(y, initial_mean, &
      initial_covariance, observation_variance, component_variance, &
      observation_prior_shape, observation_prior_rate, component_prior_shape, &
      component_prior_rate, burn, state_normal_draws, gamma_draws) result(out)
      !! Draw a local-linear-trend posterior from supplied random variates.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(2) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(2, 2) !! Initial state covariance matrix.
      real(dp), intent(in) :: observation_variance !! Observation-error variance.
      real(dp), intent(in) :: component_variance(2) !! Component variance.
      real(dp), intent(in) :: observation_prior_shape !! Observation prior shape.
      real(dp), intent(in) :: observation_prior_rate !! Observation prior rate.
      real(dp), intent(in) :: component_prior_shape(2) !! Component prior shape.
      real(dp), intent(in) :: component_prior_rate(2) !! Component prior rate.
      integer, intent(in) :: burn !! Number of initial simulation draws to discard.
      real(dp), intent(in) :: state_normal_draws(:, :, :) !! Independent standard-normal state draws.
      real(dp), intent(in) :: gamma_draws(:, :) !! Gamma simulation draws.
      type(bsts_mcmc_t) :: out
      real(dp) :: transition(2, 2), observation(2)

      transition = reshape([1.0_dp, 0.0_dp, 1.0_dp, 1.0_dp], [2, 2])
      observation = [1.0_dp, 0.0_dp]
      out = structural_gibbs_draws(y, initial_mean, initial_covariance, &
         transition, observation, observation_variance, component_variance, &
         observation_prior_shape, observation_prior_rate, &
         component_prior_shape, component_prior_rate, burn, &
         state_normal_draws, gamma_draws)
   end function bsts_local_linear_trend_draws

   function bsts_local_linear_trend(y, iterations, burn, initial_mean, &
      initial_covariance, observation_variance, component_variance, &
      observation_prior_shape, observation_prior_rate, component_prior_shape, &
      component_prior_rate) result(out)
      !! Sample a local-linear-trend posterior using the shared random stream.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      integer, intent(in) :: iterations !! Number of algorithm iterations.
      integer, intent(in), optional :: burn !! Number of initial simulation draws to discard.
      real(dp), intent(in), optional :: initial_mean(2) !! Initial state mean.
      real(dp), intent(in), optional :: initial_covariance(2, 2) !! Initial state covariance matrix.
      real(dp), intent(in), optional :: observation_variance !! Observation-error variance.
      real(dp), intent(in), optional :: component_variance(2) !! Component variance.
      real(dp), intent(in), optional :: observation_prior_shape !! Observation prior shape.
      real(dp), intent(in), optional :: observation_prior_rate !! Observation prior rate.
      real(dp), intent(in), optional :: component_prior_shape(2) !! Component prior shape.
      real(dp), intent(in), optional :: component_prior_rate(2) !! Component prior rate.
      type(bsts_mcmc_t) :: out
      real(dp), allocatable :: normal_draws(:, :, :), gamma_draws(:, :)
      real(dp) :: mean0(2), covariance0(2, 2), observation0
      real(dp) :: state0(2), observation_shape, observation_rate
      real(dp) :: state_shape(2), state_rate(2), series_variance
      integer :: discarded, iteration, component, observed

      if (size(y) < 3 .or. iterations < 1) then
         out%info = 1
         return
      end if
      discarded = 0
      if (present(burn)) discarded = burn
      series_variance = finite_sample_variance(y)
      mean0 = [first_finite_value(y), finite_endpoint_slope(y)]
      if (present(initial_mean)) mean0 = initial_mean
      covariance0 = 0.0_dp
      covariance0(1, 1) = max(series_variance, 1.0_dp)
      covariance0(2, 2) = max(series_variance, 1.0_dp)
      if (present(initial_covariance)) covariance0 = initial_covariance
      observation0 = max(0.5_dp*series_variance, 1.0e-6_dp)
      if (present(observation_variance)) observation0 = observation_variance
      state0 = max(0.01_dp*series_variance, 1.0e-8_dp)
      if (present(component_variance)) state0 = component_variance
      observation_shape = 0.01_dp
      if (present(observation_prior_shape)) &
         observation_shape = observation_prior_shape
      observation_rate = max(0.01_dp*series_variance, 1.0e-8_dp)
      if (present(observation_prior_rate)) observation_rate = &
         observation_prior_rate
      state_shape = 0.01_dp
      if (present(component_prior_shape)) state_shape = component_prior_shape
      state_rate = max(0.0001_dp*series_variance, 1.0e-10_dp)
      if (present(component_prior_rate)) state_rate = component_prior_rate
      if (discarded < 0 .or. discarded >= iterations .or. &
         observation0 <= 0.0_dp .or. any(state0 <= 0.0_dp) .or. &
         observation_shape <= 0.0_dp .or. observation_rate <= 0.0_dp .or. &
         any(state_shape <= 0.0_dp) .or. any(state_rate <= 0.0_dp) .or. &
         any([covariance0(1, 1), covariance0(2, 2)] <= 0.0_dp)) then
         out%info = 1
         return
      end if
      allocate(normal_draws(2, size(y), iterations), &
         gamma_draws(3, iterations))
      observed = count(ieee_is_finite(y))
      do iteration = 1, iterations
         call fill_standard_normals(normal_draws(:, :, iteration))
         gamma_draws(1, iteration) = random_gamma( &
            observation_shape + 0.5_dp*real(observed, dp))
         do component = 1, 2
            gamma_draws(component + 1, iteration) = random_gamma( &
               state_shape(component) + 0.5_dp*real(size(y) - 1, dp))
         end do
      end do
      out = bsts_local_linear_trend_draws(y, mean0, covariance0, &
         observation0, state0, observation_shape, observation_rate, &
         state_shape, state_rate, discarded, normal_draws, gamma_draws)
   end function bsts_local_linear_trend

   pure function bsts_student_local_linear_trend_draws(y, initial_mean, &
      initial_covariance, observation_variance, component_variance, &
      initial_degrees_of_freedom, observation_prior_shape, &
      observation_prior_rate, component_prior_shape, component_prior_rate, &
      degrees_lower, degrees_upper, degrees_proposal_sd, burn, save_weights, &
      state_normal_draws, weight_normal_draws, weight_uniform_draws, &
      variance_gamma_draws, degrees_normal_draws, degrees_uniform_draws) &
      result(out)
      !! Sample a Student-t local-linear trend from supplied random streams.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(2) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(2, 2) !! Initial state covariance matrix.
      real(dp), intent(in) :: observation_variance !! Observation-error variance.
      real(dp), intent(in) :: component_variance(2) !! Component variance.
      real(dp), intent(in) :: initial_degrees_of_freedom(2) !! Initial degrees of freedom.
      real(dp), intent(in) :: observation_prior_shape !! Observation prior shape.
      real(dp), intent(in) :: observation_prior_rate !! Observation prior rate.
      real(dp), intent(in) :: component_prior_shape(2) !! Component prior shape.
      real(dp), intent(in) :: component_prior_rate(2) !! Component prior rate.
      real(dp), intent(in) :: degrees_lower(2) !! Degrees lower.
      real(dp), intent(in) :: degrees_upper(2) !! Degrees upper.
      real(dp), intent(in) :: degrees_proposal_sd(2) !! Degrees proposal standard deviation.
      integer, intent(in) :: burn !! Number of initial simulation draws to discard.
      logical, intent(in) :: save_weights !! Flag controlling save weights.
      real(dp), intent(in) :: state_normal_draws(:, :, :) !! Independent standard-normal state draws.
      real(dp), intent(in) :: weight_normal_draws(:, :, :, :) !! Weight normal simulation draws.
      real(dp), intent(in) :: weight_uniform_draws(:, :, :, :) !! Weight uniform simulation draws.
      real(dp), intent(in) :: variance_gamma_draws(:, :) !! Variance gamma simulation draws.
      real(dp), intent(in) :: degrees_normal_draws(:, :) !! Degrees normal simulation draws.
      real(dp), intent(in) :: degrees_uniform_draws(:, :) !! Degrees uniform simulation draws.
      type(bsts_mcmc_t) :: out
      type(ssm_model_t) :: model
      real(dp), allocatable :: path(:, :), disturbance(:, :)
      real(dp) :: weights(2, size(y) - 1), scale(2), degrees(2)
      real(dp) :: observation_scale, observation_sum, component_sum
      real(dp) :: rate, gamma_value, proposal, log_ratio
      integer :: iterations, attempts, iteration, time, component, status

      iterations = size(state_normal_draws, 3)
      attempts = size(weight_normal_draws, 3)
      if (size(y) < 2 .or. count(ieee_is_finite(y)) < 1 .or. &
         observation_variance <= 0.0_dp .or. any(component_variance <= 0.0_dp) .or. &
         any(initial_degrees_of_freedom <= 0.0_dp) .or. &
         observation_prior_shape <= 0.0_dp .or. observation_prior_rate <= 0.0_dp .or. &
         any(component_prior_shape <= 0.0_dp) .or. &
         any(component_prior_rate <= 0.0_dp) .or. &
         any(degrees_lower <= 0.0_dp) .or. any(degrees_upper <= degrees_lower) .or. &
         any(initial_degrees_of_freedom < degrees_lower) .or. &
         any(initial_degrees_of_freedom > degrees_upper) .or. &
         any(degrees_proposal_sd <= 0.0_dp) .or. iterations < 1 .or. &
         attempts < 1 .or. burn < 0 .or. burn >= iterations .or. &
         any(shape(state_normal_draws) /= [2, size(y), iterations]) .or. &
         any(shape(weight_normal_draws) /= [2, size(y) - 1, attempts, iterations]) .or. &
         any(shape(weight_uniform_draws) /= &
         [2, size(y) - 1, attempts, iterations]) .or. &
         any(shape(variance_gamma_draws) /= [3, iterations]) .or. &
         any(shape(degrees_normal_draws) /= [2, iterations]) .or. &
         any(shape(degrees_uniform_draws) /= [2, iterations]) .or. &
         any(weight_uniform_draws <= 0.0_dp) .or. &
         any(weight_uniform_draws >= 1.0_dp) .or. &
         any(variance_gamma_draws <= 0.0_dp) .or. &
         any(degrees_uniform_draws <= 0.0_dp) .or. &
         any(degrees_uniform_draws >= 1.0_dp) .or. &
         .not. all(ieee_is_finite(state_normal_draws)) .or. &
         .not. all(ieee_is_finite(weight_normal_draws)) .or. &
         .not. all(ieee_is_finite(degrees_normal_draws))) then
         out%info = 1
         return
      end if
      allocate(model%y(size(y), 1), model%z(1, 2, 1), &
         model%h(1, 1, 1), model%transition(2, 2, 1), &
         model%r(2, 2, size(y)), model%q(2, 2, size(y)), &
         model%a1(2), model%p1(2, 2), model%p1inf(2, 2), &
         model%missing(size(y), 1), disturbance(2, size(y) - 1), &
         out%state(2, size(y), iterations), &
         out%observation_variance(iterations), &
         out%state_variance(2, iterations), &
         out%component_variance(2, iterations), &
         out%degrees_of_freedom(2, iterations), &
         out%transition(2, 2), out%observation(2), &
         out%transition_schedule(2, 2, 1), &
         out%state_loading_schedule(2, 2, 1))
      if (save_weights) allocate(out%state_weights(2, size(y) - 1, iterations))
      model%y(:, 1) = y
      where (.not. ieee_is_finite(model%y)) model%y = 0.0_dp
      model%z(1, :, 1) = [1.0_dp, 0.0_dp]
      model%h(1, 1, 1) = observation_variance
      model%transition(:, :, 1) = reshape( &
         [1.0_dp, 0.0_dp, 1.0_dp, 1.0_dp], [2, 2])
      model%r = 0.0_dp
      do time = 1, size(y)
         model%r(1, 1, time) = 1.0_dp
         model%r(2, 2, time) = 1.0_dp
      end do
      model%a1 = initial_mean
      model%p1 = initial_covariance
      model%p1inf = 0.0_dp
      model%missing(:, 1) = .not. ieee_is_finite(y)
      observation_scale = observation_variance
      scale = component_variance
      degrees = initial_degrees_of_freedom
      weights = 1.0_dp
      out%transition = model%transition(:, :, 1)
      out%transition_schedule(:, :, 1) = out%transition
      out%observation = [1.0_dp, 0.0_dp]
      out%state_loading_schedule = 0.0_dp
      out%state_loading_schedule(1, 1, 1) = 1.0_dp
      out%state_loading_schedule(2, 2, 1) = 1.0_dp
      out%burn = burn
      do iteration = 1, iterations
         model%h(1, 1, 1) = observation_scale
         model%q = 0.0_dp
         do time = 1, size(y) - 1
            model%q(1, 1, time) = scale(1)/weights(1, time)
            model%q(2, 2, time) = scale(2)/weights(2, time)
         end do
         model%q(:, :, size(y)) = model%q(:, :, size(y) - 1)
         call state_path_from_draws(model, state_normal_draws(:, :, iteration), &
            path, status)
         if (status /= 0) then
            out%info = 10 + status
            return
         end if
         out%state(:, :, iteration) = path
         observation_sum = 0.0_dp
         do time = 1, size(y)
            if (ieee_is_finite(y(time))) observation_sum = observation_sum + &
               (y(time) - path(1, time))**2
         end do
         observation_scale = (observation_prior_rate + &
            0.5_dp*observation_sum)/variance_gamma_draws(1, iteration)
         do time = 1, size(y) - 1
            disturbance(:, time) = path(:, time + 1) - &
               matmul(model%transition(:, :, 1), path(:, time))
         end do
         do component = 1, 2
            do time = 1, size(y) - 1
               rate = 0.5_dp*(degrees(component) + &
                  disturbance(component, time)**2/scale(component))
               call gamma_from_proposals(0.5_dp*(degrees(component) + 1.0_dp), &
                  weight_normal_draws(component, time, :, iteration), &
                  weight_uniform_draws(component, time, :, iteration), &
                  gamma_value, status)
               if (status /= 0) then
                  out%info = 50 + status
                  return
               end if
               weights(component, time) = gamma_value/rate
            end do
            component_sum = sum(weights(component, :)* &
               disturbance(component, :)**2)
            scale(component) = (component_prior_rate(component) + &
               0.5_dp*component_sum)/variance_gamma_draws(component + 1, iteration)
            proposal = reflected_value(degrees(component) + &
               degrees_proposal_sd(component)*degrees_normal_draws(component, iteration), &
               degrees_lower(component), degrees_upper(component))
            log_ratio = student_weight_log_density(weights(component, :), proposal) - &
               student_weight_log_density(weights(component, :), degrees(component))
            if (log(degrees_uniform_draws(component, iteration)) < &
               min(0.0_dp, log_ratio)) degrees(component) = proposal
         end do
         out%observation_variance(iteration) = observation_scale
         out%state_variance(:, iteration) = scale
         out%component_variance(:, iteration) = scale
         out%degrees_of_freedom(:, iteration) = degrees
         if (save_weights) out%state_weights(:, :, iteration) = weights
      end do
   end function bsts_student_local_linear_trend_draws

   function bsts_student_local_linear_trend(y, iterations, burn, save_weights, &
      proposal_attempts, initial_degrees_of_freedom, degrees_lower, &
      degrees_upper, degrees_proposal_sd) result(out)
      !! Sample a Student-t local-linear trend using the random stream.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      integer, intent(in) :: iterations !! Number of algorithm iterations.
      integer, intent(in), optional :: burn !! Number of initial simulation draws to discard.
      integer, intent(in), optional :: proposal_attempts !! Proposal attempts.
      logical, intent(in), optional :: save_weights !! Flag controlling save weights.
      real(dp), intent(in), optional :: initial_degrees_of_freedom(2) !! Initial degrees of freedom.
      real(dp), intent(in), optional :: degrees_lower(2) !! Degrees lower.
      real(dp), intent(in), optional :: degrees_upper(2) !! Degrees upper.
      real(dp), intent(in), optional :: degrees_proposal_sd(2) !! Degrees proposal standard deviation.
      type(bsts_mcmc_t) :: out
      real(dp), allocatable :: state_normals(:, :, :)
      real(dp), allocatable :: weight_normals(:, :, :, :)
      real(dp), allocatable :: weight_uniforms(:, :, :, :)
      real(dp), allocatable :: variance_gammas(:, :), degrees_normals(:, :)
      real(dp), allocatable :: degrees_uniforms(:, :)
      real(dp) :: series_variance, mean0(2), covariance0(2, 2)
      real(dp) :: observation0, component0(2), degrees0(2)
      real(dp) :: lower0(2), upper0(2), proposal_sd0(2)
      integer :: discarded, attempts, iteration, time, component, proposal
      integer :: observed
      logical :: retain_weights

      if (size(y) < 2 .or. count(ieee_is_finite(y)) < 1 .or. iterations < 1) then
         out%info = 1
         return
      end if
      discarded = 0
      if (present(burn)) discarded = burn
      attempts = 32
      if (present(proposal_attempts)) attempts = proposal_attempts
      retain_weights = .false.
      if (present(save_weights)) retain_weights = save_weights
      if (discarded < 0 .or. discarded >= iterations .or. attempts < 1) then
         out%info = 1
         return
      end if
      series_variance = max(finite_sample_variance(y), 1.0e-6_dp)
      mean0 = [first_finite_value(y), 0.0_dp]
      covariance0 = reshape([series_variance, 0.0_dp, 0.0_dp, &
         0.1_dp*series_variance], [2, 2])
      observation0 = max(0.1_dp*series_variance, 1.0e-8_dp)
      component0 = [max(0.01_dp*series_variance, 1.0e-10_dp), &
         max(0.001_dp*series_variance, 1.0e-12_dp)]
      degrees0 = 10.0_dp
      if (present(initial_degrees_of_freedom)) degrees0 = &
         initial_degrees_of_freedom
      lower0 = 1.0_dp
      if (present(degrees_lower)) lower0 = degrees_lower
      upper0 = 500.0_dp
      if (present(degrees_upper)) upper0 = degrees_upper
      proposal_sd0 = 1.0_dp
      if (present(degrees_proposal_sd)) proposal_sd0 = degrees_proposal_sd
      allocate(state_normals(2, size(y), iterations), &
         weight_normals(2, size(y) - 1, attempts, iterations), &
         weight_uniforms(2, size(y) - 1, attempts, iterations), &
         variance_gammas(3, iterations), degrees_normals(2, iterations), &
         degrees_uniforms(2, iterations))
      observed = count(ieee_is_finite(y))
      do iteration = 1, iterations
         call fill_standard_normals(state_normals(:, :, iteration))
         do proposal = 1, attempts
            do time = 1, size(y) - 1
               do component = 1, 2
                  weight_normals(component, time, proposal, iteration) = &
                     random_standard_normal()
                  weight_uniforms(component, time, proposal, iteration) = &
                     max(tiny(1.0_dp), min(1.0_dp - epsilon(1.0_dp), &
                     random_uniform()))
               end do
            end do
         end do
         variance_gammas(1, iteration) = random_gamma(2.0_dp + &
            0.5_dp*real(observed, dp))
         do component = 1, 2
            variance_gammas(component + 1, iteration) = random_gamma(2.0_dp + &
               0.5_dp*real(size(y) - 1, dp))
            degrees_normals(component, iteration) = random_standard_normal()
            degrees_uniforms(component, iteration) = max(tiny(1.0_dp), &
               min(1.0_dp - epsilon(1.0_dp), random_uniform()))
         end do
      end do
      out = bsts_student_local_linear_trend_draws(y, mean0, covariance0, &
         observation0, component0, degrees0, 2.0_dp, observation0, &
         [2.0_dp, 2.0_dp], component0, lower0, upper0, proposal_sd0, &
         discarded, retain_weights, state_normals, weight_normals, &
         weight_uniforms, variance_gammas, degrees_normals, degrees_uniforms)
   end function bsts_student_local_linear_trend

   pure function bsts_student_trend_predict_draws(fit, horizon, &
      state_student_draws, observation_normal_draws) result(out)
      !! Forecast a Student local-linear trend from supplied standardized draws.
      type(bsts_mcmc_t), intent(in) :: fit !! Previously fitted model.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      real(dp), intent(in) :: state_student_draws(:, :, :) !! State student draws.
      real(dp), intent(in) :: observation_normal_draws(:, :) !! Observation normal draws.
      type(bsts_prediction_t) :: out
      real(dp), allocatable :: state(:), sorted(:)
      integer :: retained, draw, source, step

      if (fit%info /= 0 .or. .not. allocated(fit%state) .or. &
         .not. allocated(fit%degrees_of_freedom) .or. &
         .not. allocated(fit%component_variance) .or. &
         .not. allocated(fit%observation_variance)) then
         out%info = 1
         return
      end if
      retained = size(fit%state, 3) - fit%burn
      if (horizon < 1 .or. retained < 1 .or. &
         any(shape(state_student_draws) /= [2, horizon, retained]) .or. &
         any(shape(observation_normal_draws) /= [horizon, retained]) .or. &
         .not. all(ieee_is_finite(state_student_draws)) .or. &
         .not. all(ieee_is_finite(observation_normal_draws))) then
         out%info = 1
         return
      end if
      allocate(out%draws(horizon, retained), out%mean(horizon), &
         out%standard_deviation(horizon), out%lower(horizon), &
         out%upper(horizon), state(2))
      do draw = 1, retained
         source = fit%burn + draw
         state = fit%state(:, size(fit%state, 2), source)
         do step = 1, horizon
            state = matmul(reshape([1.0_dp, 0.0_dp, 1.0_dp, 1.0_dp], &
               [2, 2]), state) + sqrt(fit%component_variance(:, source))* &
               state_student_draws(:, step, draw)
            out%draws(step, draw) = state(1) + &
               sqrt(fit%observation_variance(source))* &
               observation_normal_draws(step, draw)
         end do
      end do
      do step = 1, horizon
         out%mean(step) = sum(out%draws(step, :))/real(retained, dp)
         if (retained > 1) then
            out%standard_deviation(step) = sqrt(sum((out%draws(step, :) - &
               out%mean(step))**2)/real(retained - 1, dp))
         else
            out%standard_deviation(step) = 0.0_dp
         end if
         sorted = out%draws(step, :)
         call insertion_sort(sorted)
         out%lower(step) = quantile(sorted, 0.025_dp)
         out%upper(step) = quantile(sorted, 0.975_dp)
      end do
   end function bsts_student_trend_predict_draws

   function bsts_student_trend_predict(fit, horizon) result(out)
      !! Forecast a Student local-linear trend using the random stream.
      type(bsts_mcmc_t), intent(in) :: fit !! Previously fitted model.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      type(bsts_prediction_t) :: out
      real(dp), allocatable :: state_draws(:, :, :), observation_draws(:, :)
      real(dp) :: chi_square
      integer :: retained, draw, source, step, component

      if (.not. allocated(fit%state) .or. &
         .not. allocated(fit%degrees_of_freedom)) then
         out%info = 1
         return
      end if
      retained = size(fit%state, 3) - fit%burn
      if (horizon < 1 .or. retained < 1) then
         out%info = 1
         return
      end if
      allocate(state_draws(2, horizon, retained), &
         observation_draws(horizon, retained))
      do draw = 1, retained
         source = fit%burn + draw
         do step = 1, horizon
            do component = 1, 2
               chi_square = 2.0_dp*random_gamma( &
                  0.5_dp*fit%degrees_of_freedom(component, source))
               state_draws(component, step, draw) = random_standard_normal()/ &
                  sqrt(chi_square/fit%degrees_of_freedom(component, source))
            end do
            observation_draws(step, draw) = random_standard_normal()
         end do
      end do
      out = bsts_student_trend_predict_draws(fit, horizon, state_draws, &
         observation_draws)
   end function bsts_student_trend_predict

   pure function bsts_semilocal_trend_draws(y, initial_mean, &
      initial_covariance, observation_variance, component_variance, &
      slope_mean, slope_ar, observation_prior_shape, observation_prior_rate, &
      component_prior_shape, component_prior_rate, slope_mean_prior_mean, &
      slope_mean_prior_variance, slope_ar_prior_mean, slope_ar_prior_variance, &
      force_stationary, force_positive, burn, state_normal_draws, gamma_draws, &
      parameter_normal_draws, slope_ar_uniform_draws) result(out)
      !! Draw a semilocal trend posterior from supplied independent variates.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(2) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(2, 2) !! Initial state covariance matrix.
      real(dp), intent(in) :: observation_variance !! Observation-error variance.
      real(dp), intent(in) :: component_variance(2) !! Component variance.
      real(dp), intent(in) :: slope_mean !! Slope mean.
      real(dp), intent(in) :: slope_ar !! Slope autoregressive.
      real(dp), intent(in) :: observation_prior_shape !! Observation prior shape.
      real(dp), intent(in) :: observation_prior_rate !! Observation prior rate.
      real(dp), intent(in) :: component_prior_shape(2) !! Component prior shape.
      real(dp), intent(in) :: component_prior_rate(2) !! Component prior rate.
      real(dp), intent(in) :: slope_mean_prior_mean !! Slope mean prior mean.
      real(dp), intent(in) :: slope_mean_prior_variance !! Slope mean prior variance.
      real(dp), intent(in) :: slope_ar_prior_mean !! Slope autoregressive prior mean.
      real(dp), intent(in) :: slope_ar_prior_variance !! Slope autoregressive prior variance.
      logical, intent(in) :: force_stationary !! Flag controlling force stationary.
      logical, intent(in) :: force_positive !! Flag controlling force positive.
      integer, intent(in) :: burn !! Number of initial simulation draws to discard.
      real(dp), intent(in) :: state_normal_draws(:, :, :) !! Independent standard-normal state draws.
      real(dp), intent(in) :: gamma_draws(:, :) !! Gamma simulation draws.
      real(dp), intent(in) :: parameter_normal_draws(:, :) !! Parameter normal simulation draws.
      real(dp), intent(in) :: slope_ar_uniform_draws(:) !! Slope autoregressive uniform simulation draws.
      type(bsts_mcmc_t) :: out
      type(ssm_model_t) :: model
      real(dp), allocatable :: adjusted_y(:), centered_path(:, :), path(:, :)
      real(dp) :: centered_mean(2), transition(2, 2), state_scale(2)
      real(dp) :: observation_scale, long_run_slope, ar_coefficient
      real(dp) :: observation_sum, level_sum, slope_sum, residual
      real(dp) :: posterior_variance, posterior_mean, precision, coefficient
      real(dp) :: lower_bound, upper_bound
      integer :: iterations, iteration, time, component, status, observed

      iterations = size(state_normal_draws, 3)
      if (size(y) < 3 .or. iterations < 1 .or. burn < 0 .or. &
         burn >= iterations .or. observation_variance <= 0.0_dp .or. &
         any(component_variance <= 0.0_dp) .or. &
         observation_prior_shape <= 0.0_dp .or. observation_prior_rate <= 0.0_dp .or. &
         any(component_prior_shape <= 0.0_dp) .or. &
         any(component_prior_rate <= 0.0_dp) .or. &
         slope_mean_prior_variance <= 0.0_dp .or. &
         slope_ar_prior_variance <= 0.0_dp .or. &
         size(state_normal_draws, 1) /= 2 .or. &
         size(state_normal_draws, 2) /= size(y) .or. &
         any(shape(gamma_draws) /= [3, iterations]) .or. &
         any(shape(parameter_normal_draws) /= [2, iterations]) .or. &
         size(slope_ar_uniform_draws) /= iterations .or. &
         any(gamma_draws <= 0.0_dp) .or. &
         any(slope_ar_uniform_draws <= 0.0_dp) .or. &
         any(slope_ar_uniform_draws >= 1.0_dp) .or. &
         .not. all(ieee_is_finite(state_normal_draws)) .or. &
         .not. all(ieee_is_finite(parameter_normal_draws))) then
         out%info = 1
         return
      end if
      if (force_stationary .and. abs(slope_ar) >= 1.0_dp) then
         out%info = 1
         return
      end if
      if (force_positive .and. slope_ar <= 0.0_dp) then
         out%info = 1
         return
      end if
      allocate(adjusted_y(size(y)), centered_path(2, size(y)), &
         path(2, size(y)), out%state(2, size(y), iterations), &
         out%observation_variance(iterations), &
         out%state_variance(2, iterations), &
         out%component_variance(2, iterations), out%slope_mean(iterations), &
         out%slope_ar(iterations), out%transition(2, 2), &
         out%observation(2), out%transition_schedule(2, 2, 1), &
         out%state_loading_schedule(2, 2, 1))
      out%observation = [1.0_dp, 0.0_dp]
      out%state_loading_schedule = 0.0_dp
      out%state_loading_schedule(1, 1, 1) = 1.0_dp
      out%state_loading_schedule(2, 2, 1) = 1.0_dp
      out%burn = burn
      out%is_semilocal = .true.
      observed = count(ieee_is_finite(y))
      observation_scale = observation_variance
      state_scale = component_variance
      long_run_slope = slope_mean
      ar_coefficient = slope_ar
      lower_bound = -huge(1.0_dp)
      upper_bound = huge(1.0_dp)
      if (force_stationary) then
         lower_bound = -1.0_dp
         upper_bound = 1.0_dp
      end if
      if (force_positive) lower_bound = max(0.0_dp, lower_bound)
      do iteration = 1, iterations
         do time = 1, size(y)
            adjusted_y(time) = y(time) - real(time - 1, dp)*long_run_slope
         end do
         centered_mean = [initial_mean(1), initial_mean(2) - long_run_slope]
         transition = reshape([1.0_dp, 0.0_dp, 1.0_dp, ar_coefficient], &
            [2, 2])
         model = structural_model(adjusted_y, centered_mean, &
            initial_covariance, transition, out%observation, &
            observation_scale, state_scale)
         call state_path_from_draws(model, state_normal_draws(:, :, iteration), &
            centered_path, status)
         if (status /= 0) then
            out%info = 10 + status
            return
         end if
         do time = 1, size(y)
            path(1, time) = centered_path(1, time) + &
               real(time - 1, dp)*long_run_slope
            path(2, time) = centered_path(2, time) + long_run_slope
         end do
         out%state(:, :, iteration) = path
         observation_sum = 0.0_dp
         do time = 1, size(y)
            if (ieee_is_finite(y(time))) observation_sum = observation_sum + &
               (y(time) - path(1, time))**2
         end do
         observation_scale = (observation_prior_rate + &
            0.5_dp*observation_sum)/gamma_draws(1, iteration)
         level_sum = 0.0_dp
         slope_sum = 0.0_dp
         do time = 2, size(y)
            residual = path(1, time) - path(1, time - 1) - path(2, time - 1)
            level_sum = level_sum + residual**2
            residual = path(2, time) - long_run_slope - &
               ar_coefficient*(path(2, time - 1) - long_run_slope)
            slope_sum = slope_sum + residual**2
         end do
         state_scale(1) = (component_prior_rate(1) + 0.5_dp*level_sum)/ &
            gamma_draws(2, iteration)
         state_scale(2) = (component_prior_rate(2) + 0.5_dp*slope_sum)/ &
            gamma_draws(3, iteration)
         coefficient = 1.0_dp - ar_coefficient
         precision = 1.0_dp/slope_mean_prior_variance + &
            real(size(y) - 1, dp)*coefficient**2/state_scale(2)
         posterior_variance = 1.0_dp/precision
         posterior_mean = posterior_variance*(slope_mean_prior_mean/ &
            slope_mean_prior_variance + coefficient*sum(path(2, 2:) - &
            ar_coefficient*path(2, :size(y) - 1))/state_scale(2))
         long_run_slope = posterior_mean + sqrt(posterior_variance)* &
            parameter_normal_draws(1, iteration)
         precision = 1.0_dp/slope_ar_prior_variance + &
            sum((path(2, :size(y) - 1) - long_run_slope)**2)/state_scale(2)
         posterior_variance = 1.0_dp/precision
         posterior_mean = posterior_variance*(slope_ar_prior_mean/ &
            slope_ar_prior_variance + sum((path(2, :size(y) - 1) - &
            long_run_slope)*(path(2, 2:) - long_run_slope))/state_scale(2))
         if (force_stationary .or. force_positive) then
            ar_coefficient = truncated_normal_draw(posterior_mean, &
               sqrt(posterior_variance), lower_bound, upper_bound, &
               slope_ar_uniform_draws(iteration))
         else
            ar_coefficient = posterior_mean + sqrt(posterior_variance)* &
               parameter_normal_draws(2, iteration)
         end if
         out%observation_variance(iteration) = observation_scale
         out%state_variance(:, iteration) = state_scale
         out%component_variance(:, iteration) = state_scale
         out%slope_mean(iteration) = long_run_slope
         out%slope_ar(iteration) = ar_coefficient
      end do
      out%transition = reshape([1.0_dp, 0.0_dp, 1.0_dp, ar_coefficient], &
         [2, 2])
      out%transition_schedule(:, :, 1) = out%transition
   end function bsts_semilocal_trend_draws

   function bsts_semilocal_trend(y, iterations, burn, initial_mean, &
      initial_covariance, observation_variance, component_variance, &
      slope_mean, slope_ar, observation_prior_shape, observation_prior_rate, &
      component_prior_shape, component_prior_rate, slope_mean_prior_mean, &
      slope_mean_prior_variance, slope_ar_prior_mean, slope_ar_prior_variance, &
      force_stationary, force_positive) result(out)
      !! Sample a semilocal trend posterior using the shared random stream.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      integer, intent(in) :: iterations !! Number of algorithm iterations.
      integer, intent(in), optional :: burn !! Number of initial simulation draws to discard.
      real(dp), intent(in), optional :: initial_mean(2) !! Initial state mean.
      real(dp), intent(in), optional :: initial_covariance(2, 2) !! Initial state covariance matrix.
      real(dp), intent(in), optional :: observation_variance !! Observation-error variance.
      real(dp), intent(in), optional :: component_variance(2) !! Component variance.
      real(dp), intent(in), optional :: slope_mean !! Slope mean.
      real(dp), intent(in), optional :: slope_ar !! Slope autoregressive.
      real(dp), intent(in), optional :: observation_prior_shape !! Observation prior shape.
      real(dp), intent(in), optional :: observation_prior_rate !! Observation prior rate.
      real(dp), intent(in), optional :: component_prior_shape(2) !! Component prior shape.
      real(dp), intent(in), optional :: component_prior_rate(2) !! Component prior rate.
      real(dp), intent(in), optional :: slope_mean_prior_mean !! Slope mean prior mean.
      real(dp), intent(in), optional :: slope_mean_prior_variance !! Slope mean prior variance.
      real(dp), intent(in), optional :: slope_ar_prior_mean !! Slope autoregressive prior mean.
      real(dp), intent(in), optional :: slope_ar_prior_variance !! Slope autoregressive prior variance.
      logical, intent(in), optional :: force_stationary !! Flag controlling force stationary.
      logical, intent(in), optional :: force_positive !! Flag controlling force positive.
      type(bsts_mcmc_t) :: out
      real(dp), allocatable :: state_draws(:, :, :), gamma_draws(:, :)
      real(dp), allocatable :: parameter_draws(:, :), uniform_draws(:)
      real(dp) :: mean0(2), covariance0(2, 2), observation0, state0(2)
      real(dp) :: mean_slope0, ar0, observation_shape, observation_rate
      real(dp) :: state_shape(2), state_rate(2), mean_prior, mean_prior_variance
      real(dp) :: ar_prior, ar_prior_variance, series_variance
      integer :: discarded, iteration, component, observed
      logical :: stationary, positive

      if (size(y) < 3 .or. iterations < 1) then
         out%info = 1
         return
      end if
      discarded = 0
      if (present(burn)) discarded = burn
      series_variance = finite_sample_variance(y)
      mean0 = [first_finite_value(y), 0.0_dp]
      if (present(initial_mean)) mean0 = initial_mean
      covariance0 = 0.0_dp
      covariance0(1, 1) = max(series_variance, 1.0_dp)
      covariance0(2, 2) = max(series_variance, 1.0_dp)
      if (present(initial_covariance)) covariance0 = initial_covariance
      observation0 = max(0.5_dp*series_variance, 1.0e-6_dp)
      if (present(observation_variance)) observation0 = observation_variance
      state0 = max(0.0001_dp*series_variance, 1.0e-10_dp)
      if (present(component_variance)) state0 = component_variance
      mean_slope0 = 0.0_dp
      if (present(slope_mean)) mean_slope0 = slope_mean
      ar0 = 0.8_dp
      if (present(slope_ar)) ar0 = slope_ar
      observation_shape = 0.01_dp
      if (present(observation_prior_shape)) &
         observation_shape = observation_prior_shape
      observation_rate = max(0.01_dp*series_variance, 1.0e-8_dp)
      if (present(observation_prior_rate)) observation_rate = &
         observation_prior_rate
      state_shape = 0.01_dp
      if (present(component_prior_shape)) state_shape = component_prior_shape
      state_rate = max(0.0001_dp*series_variance, 1.0e-10_dp)
      if (present(component_prior_rate)) state_rate = component_prior_rate
      mean_prior = 0.0_dp
      if (present(slope_mean_prior_mean)) mean_prior = slope_mean_prior_mean
      mean_prior_variance = max(series_variance, 1.0e-6_dp)
      if (present(slope_mean_prior_variance)) &
         mean_prior_variance = slope_mean_prior_variance
      ar_prior = 0.8_dp
      if (present(slope_ar_prior_mean)) ar_prior = slope_ar_prior_mean
      ar_prior_variance = 0.15_dp**2
      if (present(slope_ar_prior_variance)) &
         ar_prior_variance = slope_ar_prior_variance
      stationary = .true.
      if (present(force_stationary)) stationary = force_stationary
      positive = .false.
      if (present(force_positive)) positive = force_positive
      if (discarded < 0 .or. discarded >= iterations .or. &
         observation0 <= 0.0_dp .or. any(state0 <= 0.0_dp) .or. &
         observation_shape <= 0.0_dp .or. observation_rate <= 0.0_dp .or. &
         any(state_shape <= 0.0_dp) .or. any(state_rate <= 0.0_dp) .or. &
         mean_prior_variance <= 0.0_dp .or. ar_prior_variance <= 0.0_dp) then
         out%info = 1
         return
      end if
      allocate(state_draws(2, size(y), iterations), gamma_draws(3, iterations), &
         parameter_draws(2, iterations), uniform_draws(iterations))
      observed = count(ieee_is_finite(y))
      do iteration = 1, iterations
         call fill_standard_normals(state_draws(:, :, iteration))
         do component = 1, 2
            parameter_draws(component, iteration) = random_standard_normal()
         end do
         uniform_draws(iteration) = max(tiny(1.0_dp), min(1.0_dp - &
            epsilon(1.0_dp), random_uniform()))
         gamma_draws(1, iteration) = random_gamma( &
            observation_shape + 0.5_dp*real(observed, dp))
         do component = 1, 2
            gamma_draws(component + 1, iteration) = random_gamma( &
               state_shape(component) + 0.5_dp*real(size(y) - 1, dp))
         end do
      end do
      out = bsts_semilocal_trend_draws(y, mean0, covariance0, observation0, &
         state0, mean_slope0, ar0, observation_shape, observation_rate, &
         state_shape, state_rate, mean_prior, mean_prior_variance, ar_prior, &
         ar_prior_variance, stationary, positive, discarded, state_draws, &
         gamma_draws, parameter_draws, uniform_draws)
   end function bsts_semilocal_trend

   pure function bsts_seasonal_draws(y, nseasons, season_duration, &
      initial_variance, observation_variance, seasonal_variance, &
      observation_prior_shape, observation_prior_rate, seasonal_prior_shape, &
      seasonal_prior_rate, burn, state_normal_draws, gamma_draws) result(out)
      !! Draw a sum-to-zero seasonal posterior from supplied random variates.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      integer, intent(in) :: nseasons !! Nseasons.
      integer, intent(in) :: season_duration !! Season duration.
      integer, intent(in) :: burn !! Number of initial simulation draws to discard.
      real(dp), intent(in) :: initial_variance !! Initial variance.
      real(dp), intent(in) :: observation_variance !! Observation-error variance.
      real(dp), intent(in) :: seasonal_variance !! Seasonal variance.
      real(dp), intent(in) :: observation_prior_shape !! Observation prior shape.
      real(dp), intent(in) :: observation_prior_rate !! Observation prior rate.
      real(dp), intent(in) :: seasonal_prior_shape !! Seasonal prior shape.
      real(dp), intent(in) :: seasonal_prior_rate !! Seasonal prior rate.
      real(dp), intent(in) :: state_normal_draws(:, :, :) !! Independent standard-normal state draws.
      real(dp), intent(in) :: gamma_draws(:, :) !! Gamma simulation draws.
      type(bsts_mcmc_t) :: out
      type(ssm_model_t) :: model
      real(dp), allocatable :: path(:, :), initial_mean(:)
      real(dp), allocatable :: initial_covariance(:, :), seasonal_transition(:, :)
      real(dp) :: observation_scale, seasonal_scale
      real(dp) :: observation_sum, seasonal_sum, residual
      integer :: dimension, iterations, iteration, time, status, boundaries
      integer :: component

      dimension = nseasons - 1
      iterations = size(state_normal_draws, 3)
      boundaries = (size(y) - 1)/max(1, season_duration)
      if (size(y) < 2 .or. dimension < 1 .or. season_duration < 1 .or. &
         iterations < 1 .or. burn < 0 .or. burn >= iterations .or. &
         initial_variance <= 0.0_dp .or. observation_variance <= 0.0_dp .or. &
         seasonal_variance <= 0.0_dp .or. observation_prior_shape <= 0.0_dp .or. &
         observation_prior_rate <= 0.0_dp .or. seasonal_prior_shape <= 0.0_dp .or. &
         seasonal_prior_rate <= 0.0_dp .or. &
         size(state_normal_draws, 1) /= dimension .or. &
         size(state_normal_draws, 2) /= size(y) .or. &
         any(shape(gamma_draws) /= [2, iterations]) .or. &
         any(gamma_draws <= 0.0_dp) .or. &
         .not. all(ieee_is_finite(state_normal_draws))) then
         out%info = 1
         return
      end if
      allocate(initial_mean(dimension), initial_covariance(dimension, dimension), &
         seasonal_transition(dimension, dimension))
      initial_mean = 0.0_dp
      initial_covariance = 0.0_dp
      seasonal_transition = 0.0_dp
      seasonal_transition(1, :) = -1.0_dp
      do component = 1, dimension
         initial_covariance(component, component) = initial_variance
      end do
      do component = 2, dimension
         seasonal_transition(component, component - 1) = 1.0_dp
      end do
      allocate(out%state(dimension, size(y), iterations), &
         out%observation_variance(iterations), out%state_variance(1, iterations), &
         out%component_variance(1, iterations), &
         out%transition(dimension, dimension), out%observation(dimension), &
         out%transition_schedule(dimension, dimension, season_duration), &
         out%state_loading_schedule(dimension, 1, season_duration))
      out%transition = seasonal_transition
      out%observation = 0.0_dp
      out%observation(1) = 1.0_dp
      out%transition_schedule = 0.0_dp
      out%state_loading_schedule = 0.0_dp
      do time = 1, season_duration - 1
         do component = 1, dimension
            out%transition_schedule(component, component, time) = 1.0_dp
         end do
      end do
      out%transition_schedule(:, :, season_duration) = seasonal_transition
      out%state_loading_schedule(1, 1, season_duration) = 1.0_dp
      out%forecast_phase = 1 + modulo(size(y) - 1, season_duration)
      out%burn = burn
      observation_scale = observation_variance
      seasonal_scale = seasonal_variance
      do iteration = 1, iterations
         model = seasonal_model(y, nseasons, season_duration, initial_mean, &
            initial_covariance, observation_scale, seasonal_scale)
         call state_path_from_draws(model, state_normal_draws(:, :, iteration), &
            path, status)
         if (status /= 0) then
            out%info = 10 + status
            return
         end if
         out%state(:, :, iteration) = path
         observation_sum = 0.0_dp
         do time = 1, size(y)
            if (ieee_is_finite(y(time))) observation_sum = observation_sum + &
               (y(time) - path(1, time))**2
         end do
         observation_scale = (observation_prior_rate + &
            0.5_dp*observation_sum)/gamma_draws(1, iteration)
         seasonal_sum = 0.0_dp
         do time = 2, size(y)
            if (modulo(time - 1, season_duration) == 0) then
               residual = path(1, time) + sum(path(:, time - 1))
               seasonal_sum = seasonal_sum + residual**2
            end if
         end do
         seasonal_scale = (seasonal_prior_rate + 0.5_dp*seasonal_sum)/ &
            gamma_draws(2, iteration)
         out%observation_variance(iteration) = observation_scale
         out%state_variance(1, iteration) = seasonal_scale
         out%component_variance(1, iteration) = seasonal_scale
      end do
   end function bsts_seasonal_draws

   function bsts_seasonal(y, nseasons, iterations, season_duration, burn, &
      initial_variance, observation_variance, seasonal_variance, &
      observation_prior_shape, observation_prior_rate, seasonal_prior_shape, &
      seasonal_prior_rate) result(out)
      !! Sample a bsts seasonal posterior using the shared random stream.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      integer, intent(in) :: nseasons !! Nseasons.
      integer, intent(in) :: iterations !! Number of algorithm iterations.
      integer, intent(in), optional :: season_duration !! Season duration.
      integer, intent(in), optional :: burn !! Number of initial simulation draws to discard.
      real(dp), intent(in), optional :: initial_variance !! Initial variance.
      real(dp), intent(in), optional :: observation_variance !! Observation-error variance.
      real(dp), intent(in), optional :: seasonal_variance !! Seasonal variance.
      real(dp), intent(in), optional :: observation_prior_shape !! Observation prior shape.
      real(dp), intent(in), optional :: observation_prior_rate !! Observation prior rate.
      real(dp), intent(in), optional :: seasonal_prior_shape !! Seasonal prior shape.
      real(dp), intent(in), optional :: seasonal_prior_rate !! Seasonal prior rate.
      type(bsts_mcmc_t) :: out
      real(dp), allocatable :: normal_draws(:, :, :), gamma_draws(:, :)
      real(dp) :: initial0, observation0, seasonal0, observation_shape
      real(dp) :: observation_rate, seasonal_shape, seasonal_rate, series_variance
      integer :: duration, discarded, iteration, observed, boundaries

      duration = 1
      if (present(season_duration)) duration = season_duration
      discarded = 0
      if (present(burn)) discarded = burn
      if (size(y) < 2 .or. nseasons < 2 .or. duration < 1 .or. &
         iterations < 1 .or. discarded < 0 .or. discarded >= iterations) then
         out%info = 1
         return
      end if
      series_variance = finite_sample_variance(y)
      initial0 = max(series_variance, 1.0_dp)
      if (present(initial_variance)) initial0 = initial_variance
      observation0 = max(0.5_dp*series_variance, 1.0e-6_dp)
      if (present(observation_variance)) observation0 = observation_variance
      seasonal0 = max(0.0001_dp*series_variance, 1.0e-10_dp)
      if (present(seasonal_variance)) seasonal0 = seasonal_variance
      observation_shape = 0.01_dp
      if (present(observation_prior_shape)) &
         observation_shape = observation_prior_shape
      observation_rate = max(0.01_dp*series_variance, 1.0e-8_dp)
      if (present(observation_prior_rate)) observation_rate = &
         observation_prior_rate
      seasonal_shape = 0.01_dp
      if (present(seasonal_prior_shape)) seasonal_shape = seasonal_prior_shape
      seasonal_rate = max(0.0001_dp*series_variance, 1.0e-10_dp)
      if (present(seasonal_prior_rate)) seasonal_rate = seasonal_prior_rate
      if (initial0 <= 0.0_dp .or. observation0 <= 0.0_dp .or. &
         seasonal0 <= 0.0_dp .or. observation_shape <= 0.0_dp .or. &
         observation_rate <= 0.0_dp .or. seasonal_shape <= 0.0_dp .or. &
         seasonal_rate <= 0.0_dp) then
         out%info = 1
         return
      end if
      allocate(normal_draws(nseasons - 1, size(y), iterations), &
         gamma_draws(2, iterations))
      observed = count(ieee_is_finite(y))
      boundaries = (size(y) - 1)/duration
      do iteration = 1, iterations
         call fill_standard_normals(normal_draws(:, :, iteration))
         gamma_draws(1, iteration) = random_gamma( &
            observation_shape + 0.5_dp*real(observed, dp))
         gamma_draws(2, iteration) = random_gamma( &
            seasonal_shape + 0.5_dp*real(boundaries, dp))
      end do
      out = bsts_seasonal_draws(y, nseasons, duration, initial0, observation0, &
         seasonal0, observation_shape, observation_rate, seasonal_shape, &
         seasonal_rate, discarded, normal_draws, gamma_draws)
   end function bsts_seasonal

   pure function bsts_monthly_annual_cycle_draws(y, first_date, initial_mean, &
      initial_covariance, observation_variance, monthly_variance, &
      observation_prior_shape, observation_prior_rate, monthly_prior_shape, &
      monthly_prior_rate, burn, state_normal_draws, gamma_draws) result(out)
      !! Sample a calendar-driven monthly annual cycle from supplied draws.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(11) !! Initial state mean.
      type(date_t), intent(in) :: first_date !! First date.
      real(dp), intent(in) :: initial_covariance(11, 11) !! Initial state covariance matrix.
      real(dp), intent(in) :: observation_variance !! Observation-error variance.
      real(dp), intent(in) :: monthly_variance !! Monthly variance.
      real(dp), intent(in) :: observation_prior_shape !! Observation prior shape.
      real(dp), intent(in) :: observation_prior_rate !! Observation prior rate.
      real(dp), intent(in) :: monthly_prior_shape !! Monthly prior shape.
      real(dp), intent(in) :: monthly_prior_rate !! Monthly prior rate.
      integer, intent(in) :: burn !! Number of initial simulation draws to discard.
      real(dp), intent(in) :: state_normal_draws(:, :, :) !! Independent standard-normal state draws.
      real(dp), intent(in) :: gamma_draws(:, :) !! Gamma simulation draws.
      type(bsts_mcmc_t) :: out
      type(ssm_model_t) :: model
      type(date_t) :: current_date, next_date
      real(dp), allocatable :: path(:, :)
      real(dp) :: seasonal_transition(11, 11), observation_scale, monthly_scale
      real(dp) :: observation_sum, monthly_sum, residual
      integer :: iterations, iteration, time, component, status

      iterations = size(state_normal_draws, 3)
      if (size(y) < 2 .or. count(ieee_is_finite(y)) < 1 .or. &
         .not. date_valid(first_date) .or. observation_variance <= 0.0_dp .or. &
         monthly_variance <= 0.0_dp .or. observation_prior_shape <= 0.0_dp .or. &
         observation_prior_rate <= 0.0_dp .or. monthly_prior_shape <= 0.0_dp .or. &
         monthly_prior_rate <= 0.0_dp .or. iterations < 1 .or. &
         burn < 0 .or. burn >= iterations .or. &
         any(shape(state_normal_draws) /= [11, size(y), iterations]) .or. &
         any(shape(gamma_draws) /= [2, iterations]) .or. &
         any(gamma_draws <= 0.0_dp) .or. &
         .not. all(ieee_is_finite(state_normal_draws))) then
         out%info = 1
         return
      end if
      seasonal_transition = 0.0_dp
      seasonal_transition(1, :) = -1.0_dp
      do component = 2, 11
         seasonal_transition(component, component - 1) = 1.0_dp
      end do
      allocate(model%y(size(y), 1), model%z(1, 11, 1), &
         model%h(1, 1, 1), model%transition(11, 11, size(y)), &
         model%r(11, 1, size(y)), model%q(1, 1, size(y)), &
         model%a1(11), model%p1(11, 11), model%p1inf(11, 11), &
         model%missing(size(y), 1), out%state(11, size(y), iterations), &
         out%observation_variance(iterations), &
         out%state_variance(1, iterations), &
         out%component_variance(1, iterations), &
         out%transition(11, 11), out%observation(11), &
         out%transition_schedule(11, 11, size(y)), &
         out%state_loading_schedule(11, 1, size(y)))
      model%y(:, 1) = y
      where (.not. ieee_is_finite(model%y)) model%y = 0.0_dp
      model%z = 0.0_dp
      model%z(1, 1, 1) = 1.0_dp
      model%h(1, 1, 1) = observation_variance
      model%transition = 0.0_dp
      model%r = 0.0_dp
      model%q(1, 1, :) = monthly_variance
      current_date = first_date
      do time = 1, size(y)
         do component = 1, 11
            model%transition(component, component, time) = 1.0_dp
         end do
         if (time < size(y)) then
            next_date = date_from_day_number(date_day_number(current_date) + 1)
            if (next_date%month /= current_date%month) then
               model%transition(:, :, time) = seasonal_transition
               model%r(1, 1, time) = 1.0_dp
            end if
            current_date = next_date
         end if
      end do
      model%a1 = initial_mean
      model%p1 = initial_covariance
      model%p1inf = 0.0_dp
      model%missing(:, 1) = .not. ieee_is_finite(y)
      out%transition = seasonal_transition
      out%observation = 0.0_dp
      out%observation(1) = 1.0_dp
      out%transition_schedule = model%transition
      out%state_loading_schedule = model%r
      out%burn = burn
      out%is_monthly = .true.
      out%last_date = current_date
      observation_scale = observation_variance
      monthly_scale = monthly_variance
      do iteration = 1, iterations
         model%h(1, 1, 1) = observation_scale
         model%q(1, 1, :) = monthly_scale
         call state_path_from_draws(model, state_normal_draws(:, :, iteration), &
            path, status)
         if (status /= 0) then
            out%info = 10 + status
            return
         end if
         out%state(:, :, iteration) = path
         observation_sum = 0.0_dp
         do time = 1, size(y)
            if (ieee_is_finite(y(time))) observation_sum = observation_sum + &
               (y(time) - path(1, time))**2
         end do
         observation_scale = (observation_prior_rate + &
            0.5_dp*observation_sum)/gamma_draws(1, iteration)
         monthly_sum = 0.0_dp
         do time = 2, size(y)
            if (model%r(1, 1, time - 1) /= 0.0_dp) then
               residual = path(1, time) + sum(path(:, time - 1))
               monthly_sum = monthly_sum + residual**2
            end if
         end do
         monthly_scale = (monthly_prior_rate + 0.5_dp*monthly_sum)/ &
            gamma_draws(2, iteration)
         out%observation_variance(iteration) = observation_scale
         out%state_variance(1, iteration) = monthly_scale
         out%component_variance(1, iteration) = monthly_scale
      end do
   end function bsts_monthly_annual_cycle_draws

   function bsts_monthly_annual_cycle(y, first_date, iterations, burn, &
      initial_variance, observation_variance, monthly_variance, &
      observation_prior_shape, observation_prior_rate, monthly_prior_shape, &
      monthly_prior_rate) result(out)
      !! Sample a monthly annual cycle using the shared random stream.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      type(date_t), intent(in) :: first_date !! First date.
      integer, intent(in) :: iterations !! Number of algorithm iterations.
      integer, intent(in), optional :: burn !! Number of initial simulation draws to discard.
      real(dp), intent(in), optional :: initial_variance !! Initial variance.
      real(dp), intent(in), optional :: observation_variance !! Observation-error variance.
      real(dp), intent(in), optional :: monthly_variance !! Monthly variance.
      real(dp), intent(in), optional :: observation_prior_shape !! Observation prior shape.
      real(dp), intent(in), optional :: observation_prior_rate !! Observation prior rate.
      real(dp), intent(in), optional :: monthly_prior_shape !! Monthly prior shape.
      real(dp), intent(in), optional :: monthly_prior_rate !! Monthly prior rate.
      type(bsts_mcmc_t) :: out
      real(dp), allocatable :: normal_draws(:, :, :), gamma_draws(:, :)
      real(dp) :: mean0(11), covariance0(11, 11), series_variance
      real(dp) :: initial0, observation0, monthly0
      real(dp) :: observation_shape, observation_rate, monthly_shape, monthly_rate
      type(date_t) :: current_date, next_date
      integer :: discarded, iteration, component, time, observed, boundaries

      if (size(y) < 2 .or. count(ieee_is_finite(y)) < 1 .or. &
         .not. date_valid(first_date) .or. iterations < 1) then
         out%info = 1
         return
      end if
      discarded = 0
      if (present(burn)) discarded = burn
      if (discarded < 0 .or. discarded >= iterations) then
         out%info = 1
         return
      end if
      series_variance = max(finite_sample_variance(y), 1.0e-6_dp)
      initial0 = max(series_variance, 1.0_dp)
      if (present(initial_variance)) initial0 = initial_variance
      observation0 = max(0.5_dp*series_variance, 1.0e-6_dp)
      if (present(observation_variance)) observation0 = observation_variance
      monthly0 = max(0.0001_dp*series_variance, 1.0e-10_dp)
      if (present(monthly_variance)) monthly0 = monthly_variance
      observation_shape = 0.01_dp
      if (present(observation_prior_shape)) &
         observation_shape = observation_prior_shape
      observation_rate = max(0.01_dp*series_variance, 1.0e-8_dp)
      if (present(observation_prior_rate)) observation_rate = &
         observation_prior_rate
      monthly_shape = 0.01_dp
      if (present(monthly_prior_shape)) monthly_shape = monthly_prior_shape
      monthly_rate = max(0.0001_dp*series_variance, 1.0e-10_dp)
      if (present(monthly_prior_rate)) monthly_rate = monthly_prior_rate
      if (initial0 <= 0.0_dp .or. observation0 <= 0.0_dp .or. &
         monthly0 <= 0.0_dp .or. observation_shape <= 0.0_dp .or. &
         observation_rate <= 0.0_dp .or. monthly_shape <= 0.0_dp .or. &
         monthly_rate <= 0.0_dp) then
         out%info = 1
         return
      end if
      mean0 = 0.0_dp
      covariance0 = 0.0_dp
      do component = 1, 11
         covariance0(component, component) = initial0
      end do
      boundaries = 0
      current_date = first_date
      do time = 1, size(y) - 1
         next_date = date_from_day_number(date_day_number(current_date) + 1)
         if (next_date%month /= current_date%month) boundaries = boundaries + 1
         current_date = next_date
      end do
      allocate(normal_draws(11, size(y), iterations), gamma_draws(2, iterations))
      observed = count(ieee_is_finite(y))
      do iteration = 1, iterations
         call fill_standard_normals(normal_draws(:, :, iteration))
         gamma_draws(1, iteration) = random_gamma(observation_shape + &
            0.5_dp*real(observed, dp))
         gamma_draws(2, iteration) = random_gamma(monthly_shape + &
            0.5_dp*real(boundaries, dp))
      end do
      out = bsts_monthly_annual_cycle_draws(y, first_date, mean0, covariance0, &
         observation0, monthly0, observation_shape, observation_rate, &
         monthly_shape, monthly_rate, discarded, normal_draws, gamma_draws)
   end function bsts_monthly_annual_cycle

   pure function bsts_monthly_predict_draws(fit, horizon, state_normal_draws, &
      observation_normal_draws) result(out)
      !! Forecast a monthly annual cycle from supplied normal draws.
      type(bsts_mcmc_t), intent(in) :: fit !! Previously fitted model.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      real(dp), intent(in) :: state_normal_draws(:, :) !! Independent standard-normal state draws.
      real(dp), intent(in) :: observation_normal_draws(:, :) !! Observation normal draws.
      type(bsts_prediction_t) :: out
      real(dp), allocatable :: state(:), sorted(:)
      real(dp) :: transition(11, 11)
      type(date_t) :: current_date, next_date
      integer :: retained, draw, source, step, component

      if (fit%info /= 0 .or. .not. fit%is_monthly .or. &
         .not. date_valid(fit%last_date) .or. .not. allocated(fit%state) .or. &
         .not. allocated(fit%component_variance) .or. &
         .not. allocated(fit%observation_variance)) then
         out%info = 1
         return
      end if
      retained = size(fit%state, 3) - fit%burn
      if (horizon < 1 .or. retained < 1 .or. &
         any(shape(state_normal_draws) /= [horizon, retained]) .or. &
         any(shape(observation_normal_draws) /= [horizon, retained]) .or. &
         .not. all(ieee_is_finite(state_normal_draws)) .or. &
         .not. all(ieee_is_finite(observation_normal_draws))) then
         out%info = 1
         return
      end if
      transition = 0.0_dp
      transition(1, :) = -1.0_dp
      do component = 2, 11
         transition(component, component - 1) = 1.0_dp
      end do
      allocate(out%draws(horizon, retained), out%mean(horizon), &
         out%standard_deviation(horizon), out%lower(horizon), &
         out%upper(horizon), state(11))
      do draw = 1, retained
         source = fit%burn + draw
         state = fit%state(:, size(fit%state, 2), source)
         current_date = fit%last_date
         do step = 1, horizon
            next_date = date_from_day_number(date_day_number(current_date) + 1)
            if (next_date%month /= current_date%month) then
               state = matmul(transition, state)
               state(1) = state(1) + &
                  sqrt(fit%component_variance(1, source))* &
                  state_normal_draws(step, draw)
            end if
            out%draws(step, draw) = state(1) + &
               sqrt(fit%observation_variance(source))* &
               observation_normal_draws(step, draw)
            current_date = next_date
         end do
      end do
      do step = 1, horizon
         out%mean(step) = sum(out%draws(step, :))/real(retained, dp)
         if (retained > 1) then
            out%standard_deviation(step) = sqrt(sum((out%draws(step, :) - &
               out%mean(step))**2)/real(retained - 1, dp))
         else
            out%standard_deviation(step) = 0.0_dp
         end if
         sorted = out%draws(step, :)
         call insertion_sort(sorted)
         out%lower(step) = quantile(sorted, 0.025_dp)
         out%upper(step) = quantile(sorted, 0.975_dp)
      end do
   end function bsts_monthly_predict_draws

   function bsts_monthly_predict(fit, horizon) result(out)
      !! Forecast a monthly annual cycle using the shared random stream.
      type(bsts_mcmc_t), intent(in) :: fit !! Previously fitted model.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      type(bsts_prediction_t) :: out
      real(dp), allocatable :: state_normals(:, :), observation_normals(:, :)
      integer :: retained, draw, step

      if (.not. allocated(fit%state)) then
         out%info = 1
         return
      end if
      retained = size(fit%state, 3) - fit%burn
      if (horizon < 1 .or. retained < 1) then
         out%info = 1
         return
      end if
      allocate(state_normals(horizon, retained), &
         observation_normals(horizon, retained))
      do draw = 1, retained
         do step = 1, horizon
            state_normals(step, draw) = random_standard_normal()
            observation_normals(step, draw) = random_standard_normal()
         end do
      end do
      out = bsts_monthly_predict_draws(fit, horizon, state_normals, &
         observation_normals)
   end function bsts_monthly_predict

   pure function bsts_random_walk_holiday_draws(y, first_date, holiday, &
      initial_mean, initial_covariance, observation_variance, holiday_variance, &
      observation_prior_shape, observation_prior_rate, holiday_prior_shape, &
      holiday_prior_rate, burn, state_normal_draws, gamma_draws) result(out)
      !! Sample a random-walk holiday effect from supplied random draws.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      type(date_t), intent(in) :: first_date !! First date.
      type(bsts_holiday_t), intent(in) :: holiday !! Holiday.
      real(dp), intent(in) :: observation_variance !! Observation-error variance.
      real(dp), intent(in) :: holiday_variance !! Holiday variance.
      real(dp), intent(in) :: observation_prior_shape !! Observation prior shape.
      real(dp), intent(in) :: observation_prior_rate !! Observation prior rate.
      real(dp), intent(in) :: holiday_prior_shape !! Holiday prior shape.
      real(dp), intent(in) :: holiday_prior_rate !! Holiday prior rate.
      integer, intent(in) :: burn !! Number of initial simulation draws to discard.
      real(dp), intent(in) :: state_normal_draws(:, :, :) !! Independent standard-normal state draws.
      real(dp), intent(in) :: gamma_draws(:, :) !! Gamma simulation draws.
      type(bsts_mcmc_t) :: out
      type(ssm_model_t) :: model
      type(date_t) :: current_date, next_date
      real(dp), allocatable :: path(:, :)
      real(dp) :: observation_scale, holiday_scale, observation_sum
      real(dp) :: holiday_sum, residual
      integer :: width, iterations, iteration, time, component, position, status

      width = bsts_holiday_width(holiday)
      iterations = size(state_normal_draws, 3)
      if (size(y) < 2 .or. count(ieee_is_finite(y)) < 1 .or. &
         .not. date_valid(first_date) .or. holiday%info /= 0 .or. width < 1 .or. &
         size(initial_mean) /= width .or. &
         any(shape(initial_covariance) /= [width, width]) .or. &
         observation_variance <= 0.0_dp .or. holiday_variance <= 0.0_dp .or. &
         observation_prior_shape <= 0.0_dp .or. observation_prior_rate <= 0.0_dp .or. &
         holiday_prior_shape <= 0.0_dp .or. holiday_prior_rate <= 0.0_dp .or. &
         iterations < 1 .or. burn < 0 .or. burn >= iterations .or. &
         any(shape(state_normal_draws) /= [width, size(y), iterations]) .or. &
         any(shape(gamma_draws) /= [2, iterations]) .or. &
         any(gamma_draws <= 0.0_dp) .or. &
         .not. all(ieee_is_finite(state_normal_draws))) then
         out%info = 1
         return
      end if
      allocate(model%y(size(y), 1), model%z(1, width, size(y)), &
         model%h(1, 1, 1), model%transition(width, width, 1), &
         model%r(width, 1, size(y)), model%q(1, 1, size(y)), &
         model%a1(width), model%p1(width, width), &
         model%p1inf(width, width), model%missing(size(y), 1), &
         out%state(width, size(y), iterations), &
         out%observation_variance(iterations), &
         out%state_variance(1, iterations), &
         out%component_variance(1, iterations), &
         out%transition(width, width), out%observation(width), &
         out%transition_schedule(width, width, 1), &
         out%state_loading_schedule(width, 1, size(y)), &
         out%observation_schedule(width, size(y)))
      model%y(:, 1) = y
      where (.not. ieee_is_finite(model%y)) model%y = 0.0_dp
      model%z = 0.0_dp
      model%h(1, 1, 1) = observation_variance
      model%transition = 0.0_dp
      do component = 1, width
         model%transition(component, component, 1) = 1.0_dp
      end do
      model%r = 0.0_dp
      model%q(1, 1, :) = holiday_variance
      current_date = first_date
      do time = 1, size(y)
         position = bsts_holiday_position(holiday, current_date)
         if (position > 0) model%z(1, position, time) = 1.0_dp
         if (time < size(y)) then
            next_date = date_from_day_number(date_day_number(current_date) + 1)
            position = bsts_holiday_position(holiday, next_date)
            if (position > 0) model%r(position, 1, time) = 1.0_dp
            current_date = next_date
         end if
      end do
      model%a1 = initial_mean
      model%p1 = initial_covariance
      model%p1inf = 0.0_dp
      model%missing(:, 1) = .not. ieee_is_finite(y)
      out%transition = model%transition(:, :, 1)
      out%transition_schedule(:, :, 1) = out%transition
      out%observation = 0.0_dp
      out%state_loading_schedule = model%r
      out%observation_schedule = model%z(1, :, :)
      out%burn = burn
      out%last_date = current_date
      out%holiday = holiday
      out%is_holiday = .true.
      observation_scale = observation_variance
      holiday_scale = holiday_variance
      do iteration = 1, iterations
         model%h(1, 1, 1) = observation_scale
         model%q(1, 1, :) = holiday_scale
         call state_path_from_draws(model, state_normal_draws(:, :, iteration), &
            path, status)
         if (status /= 0) then
            out%info = 10 + status
            return
         end if
         out%state(:, :, iteration) = path
         observation_sum = 0.0_dp
         do time = 1, size(y)
            position = maxloc(model%z(1, :, time), dim=1)
            if (maxval(model%z(1, :, time)) == 0.0_dp) then
               residual = y(time)
            else
               residual = y(time) - path(position, time)
            end if
            if (ieee_is_finite(y(time))) observation_sum = observation_sum + &
               residual**2
         end do
         observation_scale = (observation_prior_rate + &
            0.5_dp*observation_sum)/gamma_draws(1, iteration)
         holiday_sum = 0.0_dp
         do time = 1, size(y) - 1
            position = maxloc(model%r(:, 1, time), dim=1)
            if (maxval(model%r(:, 1, time)) /= 0.0_dp) then
               residual = path(position, time + 1) - path(position, time)
               holiday_sum = holiday_sum + residual**2
            end if
         end do
         holiday_scale = (holiday_prior_rate + 0.5_dp*holiday_sum)/ &
            gamma_draws(2, iteration)
         out%observation_variance(iteration) = observation_scale
         out%state_variance(1, iteration) = holiday_scale
         out%component_variance(1, iteration) = holiday_scale
      end do
   end function bsts_random_walk_holiday_draws

   function bsts_random_walk_holiday(y, first_date, holiday, iterations, burn, &
      initial_variance, observation_variance, holiday_variance, &
      observation_prior_shape, observation_prior_rate, holiday_prior_shape, &
      holiday_prior_rate) result(out)
      !! Sample a random-walk holiday effect using the random stream.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      type(date_t), intent(in) :: first_date !! First date.
      type(bsts_holiday_t), intent(in) :: holiday !! Holiday.
      integer, intent(in) :: iterations !! Number of algorithm iterations.
      integer, intent(in), optional :: burn !! Number of initial simulation draws to discard.
      real(dp), intent(in), optional :: initial_variance !! Initial variance.
      real(dp), intent(in), optional :: observation_variance !! Observation-error variance.
      real(dp), intent(in), optional :: holiday_variance !! Holiday variance.
      real(dp), intent(in), optional :: observation_prior_shape !! Observation prior shape.
      real(dp), intent(in), optional :: observation_prior_rate !! Observation prior rate.
      real(dp), intent(in), optional :: holiday_prior_shape !! Holiday prior shape.
      real(dp), intent(in), optional :: holiday_prior_rate !! Holiday prior rate.
      type(bsts_mcmc_t) :: out
      real(dp), allocatable :: mean0(:), covariance0(:, :)
      real(dp), allocatable :: normal_draws(:, :, :), gamma_draws(:, :)
      real(dp) :: series_variance, initial0, observation0, holiday0
      real(dp) :: observation_shape, observation_rate, holiday_shape, holiday_rate
      type(date_t) :: current_date, next_date
      integer :: width, discarded, iteration, component, time, observed, updates

      width = bsts_holiday_width(holiday)
      if (size(y) < 2 .or. count(ieee_is_finite(y)) < 1 .or. &
         .not. date_valid(first_date) .or. holiday%info /= 0 .or. width < 1 .or. &
         iterations < 1) then
         out%info = 1
         return
      end if
      discarded = 0
      if (present(burn)) discarded = burn
      if (discarded < 0 .or. discarded >= iterations) then
         out%info = 1
         return
      end if
      series_variance = max(finite_sample_variance(y), 1.0e-6_dp)
      initial0 = max(series_variance, 1.0_dp)
      if (present(initial_variance)) initial0 = initial_variance
      observation0 = max(0.5_dp*series_variance, 1.0e-6_dp)
      if (present(observation_variance)) observation0 = observation_variance
      holiday0 = max(0.0001_dp*series_variance, 1.0e-10_dp)
      if (present(holiday_variance)) holiday0 = holiday_variance
      observation_shape = 0.01_dp
      if (present(observation_prior_shape)) &
         observation_shape = observation_prior_shape
      observation_rate = max(0.01_dp*series_variance, 1.0e-8_dp)
      if (present(observation_prior_rate)) observation_rate = &
         observation_prior_rate
      holiday_shape = 0.01_dp
      if (present(holiday_prior_shape)) holiday_shape = holiday_prior_shape
      holiday_rate = max(0.0001_dp*series_variance, 1.0e-10_dp)
      if (present(holiday_prior_rate)) holiday_rate = holiday_prior_rate
      if (initial0 <= 0.0_dp .or. observation0 <= 0.0_dp .or. &
         holiday0 <= 0.0_dp .or. observation_shape <= 0.0_dp .or. &
         observation_rate <= 0.0_dp .or. holiday_shape <= 0.0_dp .or. &
         holiday_rate <= 0.0_dp) then
         out%info = 1
         return
      end if
      allocate(mean0(width), covariance0(width, width), &
         normal_draws(width, size(y), iterations), gamma_draws(2, iterations))
      mean0 = 0.0_dp
      covariance0 = 0.0_dp
      do component = 1, width
         covariance0(component, component) = initial0
      end do
      updates = 0
      current_date = first_date
      do time = 1, size(y) - 1
         next_date = date_from_day_number(date_day_number(current_date) + 1)
         if (bsts_holiday_position(holiday, next_date) > 0) updates = updates + 1
         current_date = next_date
      end do
      observed = count(ieee_is_finite(y))
      do iteration = 1, iterations
         call fill_standard_normals(normal_draws(:, :, iteration))
         gamma_draws(1, iteration) = random_gamma(observation_shape + &
            0.5_dp*real(observed, dp))
         gamma_draws(2, iteration) = random_gamma(holiday_shape + &
            0.5_dp*real(updates, dp))
      end do
      out = bsts_random_walk_holiday_draws(y, first_date, holiday, mean0, &
         covariance0, observation0, holiday0, observation_shape, &
         observation_rate, holiday_shape, holiday_rate, discarded, &
         normal_draws, gamma_draws)
   end function bsts_random_walk_holiday

   pure function bsts_holiday_predict_draws(fit, horizon, state_normal_draws, &
      observation_normal_draws) result(out)
      !! Forecast a random-walk holiday effect from supplied normal draws.
      type(bsts_mcmc_t), intent(in) :: fit !! Previously fitted model.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      real(dp), intent(in) :: state_normal_draws(:, :) !! Independent standard-normal state draws.
      real(dp), intent(in) :: observation_normal_draws(:, :) !! Observation normal draws.
      type(bsts_prediction_t) :: out
      real(dp), allocatable :: state(:), sorted(:)
      type(date_t) :: current_date, next_date
      integer :: retained, width, draw, source, step, position

      if (fit%info /= 0 .or. .not. fit%is_holiday .or. &
         .not. date_valid(fit%last_date) .or. fit%holiday%info /= 0 .or. &
         .not. allocated(fit%state) .or. &
         .not. allocated(fit%component_variance) .or. &
         .not. allocated(fit%observation_variance)) then
         out%info = 1
         return
      end if
      width = size(fit%state, 1)
      retained = size(fit%state, 3) - fit%burn
      if (horizon < 1 .or. retained < 1 .or. &
         any(shape(state_normal_draws) /= [horizon, retained]) .or. &
         any(shape(observation_normal_draws) /= [horizon, retained]) .or. &
         .not. all(ieee_is_finite(state_normal_draws)) .or. &
         .not. all(ieee_is_finite(observation_normal_draws))) then
         out%info = 1
         return
      end if
      allocate(out%draws(horizon, retained), out%mean(horizon), &
         out%standard_deviation(horizon), out%lower(horizon), &
         out%upper(horizon), state(width))
      do draw = 1, retained
         source = fit%burn + draw
         state = fit%state(:, size(fit%state, 2), source)
         current_date = fit%last_date
         do step = 1, horizon
            next_date = date_from_day_number(date_day_number(current_date) + 1)
            position = bsts_holiday_position(fit%holiday, next_date)
            if (position > 0) state(position) = state(position) + &
               sqrt(fit%component_variance(1, source))* &
               state_normal_draws(step, draw)
            if (position > 0) then
               out%draws(step, draw) = state(position)
            else
               out%draws(step, draw) = 0.0_dp
            end if
            out%draws(step, draw) = out%draws(step, draw) + &
               sqrt(fit%observation_variance(source))* &
               observation_normal_draws(step, draw)
            current_date = next_date
         end do
      end do
      do step = 1, horizon
         out%mean(step) = sum(out%draws(step, :))/real(retained, dp)
         if (retained > 1) then
            out%standard_deviation(step) = sqrt(sum((out%draws(step, :) - &
               out%mean(step))**2)/real(retained - 1, dp))
         else
            out%standard_deviation(step) = 0.0_dp
         end if
         sorted = out%draws(step, :)
         call insertion_sort(sorted)
         out%lower(step) = quantile(sorted, 0.025_dp)
         out%upper(step) = quantile(sorted, 0.975_dp)
      end do
   end function bsts_holiday_predict_draws

   function bsts_holiday_predict(fit, horizon) result(out)
      !! Forecast a random-walk holiday effect using the random stream.
      type(bsts_mcmc_t), intent(in) :: fit !! Previously fitted model.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      type(bsts_prediction_t) :: out
      real(dp), allocatable :: state_normals(:, :), observation_normals(:, :)
      integer :: retained, draw, step

      if (.not. allocated(fit%state)) then
         out%info = 1
         return
      end if
      retained = size(fit%state, 3) - fit%burn
      if (horizon < 1 .or. retained < 1) then
         out%info = 1
         return
      end if
      allocate(state_normals(horizon, retained), &
         observation_normals(horizon, retained))
      do draw = 1, retained
         do step = 1, horizon
            state_normals(step, draw) = random_standard_normal()
            observation_normals(step, draw) = random_standard_normal()
         end do
      end do
      out = bsts_holiday_predict_draws(fit, horizon, state_normals, &
         observation_normals)
   end function bsts_holiday_predict

   pure function bsts_regression_holiday_draws(y, first_date, holidays, &
      coefficient_prior_mean, coefficient_prior_variance, residual_variance, &
      residual_prior_shape, residual_prior_rate, burn, coefficient_normal_draws, &
      gamma_draws, offset_draws) result(out)
      !! Sample fixed holiday-pattern regression from supplied random draws.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: coefficient_prior_mean !! Coefficient prior mean.
      type(date_t), intent(in) :: first_date !! First date.
      type(bsts_holiday_t), intent(in) :: holidays(:) !! Holidays.
      real(dp), intent(in) :: coefficient_prior_variance !! Coefficient prior variance.
      real(dp), intent(in) :: residual_variance !! Residual variance.
      real(dp), intent(in) :: residual_prior_shape !! Residual prior shape.
      real(dp), intent(in) :: residual_prior_rate !! Residual prior rate.
      integer, intent(in) :: burn !! Number of initial simulation draws to discard.
      real(dp), intent(in) :: coefficient_normal_draws(:, :) !! Coefficient normal simulation draws.
      real(dp), intent(in) :: gamma_draws(:) !! Gamma simulation draws.
      real(dp), intent(in), optional :: offset_draws(:, :) !! Offset simulation draws.
      type(bsts_holiday_regression_t) :: out
      real(dp), allocatable :: design(:, :), observed_design(:, :)
      real(dp), allocatable :: observed_y(:), working_y(:), precision(:, :)
      real(dp), allocatable :: covariance(:, :), score(:), posterior_mean(:)
      real(dp), allocatable :: coefficient(:)
      integer, allocatable :: observed_index(:)
      real(dp) :: variance, residual_sum
      integer :: observations, variables, iterations, observed
      integer :: iteration, holiday_index, component, status

      observations = size(y)
      iterations = size(coefficient_normal_draws, 2)
      design = bsts_holiday_design(first_date, observations, holidays)
      variables = size(design, 2)
      observed = count(ieee_is_finite(y))
      if (observations < 1 .or. observed < 1 .or. size(holidays) < 1 .or. &
         variables < 1 .or. size(design, 1) /= observations .or. &
         coefficient_prior_variance <= 0.0_dp .or. residual_variance <= 0.0_dp .or. &
         residual_prior_shape <= 0.0_dp .or. residual_prior_rate <= 0.0_dp .or. &
         iterations < 1 .or. burn < 0 .or. burn >= iterations .or. &
         any(shape(coefficient_normal_draws) /= [variables, iterations]) .or. &
         size(gamma_draws) /= iterations .or. any(gamma_draws <= 0.0_dp) .or. &
         .not. all(ieee_is_finite(coefficient_normal_draws))) then
         out%info = 1
         return
      end if
      if (present(offset_draws)) then
         if (any(shape(offset_draws) /= [observations, iterations]) .or. &
            .not. all(ieee_is_finite(offset_draws))) then
            out%info = 1
            return
         end if
      end if
      observed_index = pack([(component, component=1, observations)], &
         ieee_is_finite(y))
      allocate(observed_design(observed, variables), observed_y(observed), &
         working_y(observed), precision(variables, variables), &
         score(variables), posterior_mean(variables), coefficient(variables), &
         out%coefficients(variables, iterations), &
         out%residual_variance(iterations), &
         out%contribution(observations, iterations), &
         out%coefficient_offset(size(holidays) + 1), &
         out%holidays(size(holidays)))
      observed_design = design(observed_index, :)
      observed_y = y(observed_index)
      out%holidays = holidays
      out%coefficient_offset(1) = 0
      do holiday_index = 1, size(holidays)
         out%coefficient_offset(holiday_index + 1) = &
            out%coefficient_offset(holiday_index) + &
            bsts_holiday_width(holidays(holiday_index))
      end do
      out%first_date = first_date
      out%last_date = date_from_day_number(date_day_number(first_date) + &
         observations - 1)
      out%burn = burn
      variance = residual_variance
      do iteration = 1, iterations
         working_y = observed_y
         if (present(offset_draws)) working_y = working_y - &
            offset_draws(observed_index, iteration)
         precision = matmul(transpose(observed_design), observed_design)/variance
         do component = 1, variables
            precision(component, component) = precision(component, component) + &
               1.0_dp/coefficient_prior_variance
         end do
         call invert_matrix(precision, covariance, status)
         if (status /= 0) then
            out%info = 10 + status
            return
         end if
         score = matmul(transpose(observed_design), working_y)/variance + &
            coefficient_prior_mean/coefficient_prior_variance
         posterior_mean = matmul(covariance, score)
         call multivariate_normal_from_standard(posterior_mean, covariance, &
            coefficient_normal_draws(:, iteration), coefficient, status)
         if (status /= 0) then
            out%info = 20 + status
            return
         end if
         residual_sum = sum((working_y - matmul(observed_design, coefficient))**2)
         variance = (residual_prior_rate + 0.5_dp*residual_sum)/ &
            gamma_draws(iteration)
         out%coefficients(:, iteration) = coefficient
         out%residual_variance(iteration) = variance
         out%contribution(:, iteration) = matmul(design, coefficient)
      end do
   end function bsts_regression_holiday_draws

   function bsts_regression_holiday(y, first_date, holidays, iterations, burn, &
      coefficient_prior_mean, coefficient_prior_variance, residual_variance, &
      residual_prior_shape, residual_prior_rate, offset_draws) result(out)
      !! Sample fixed holiday-pattern regression using the random stream.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      type(date_t), intent(in) :: first_date !! First date.
      type(bsts_holiday_t), intent(in) :: holidays(:) !! Holidays.
      integer, intent(in) :: iterations !! Number of algorithm iterations.
      integer, intent(in), optional :: burn !! Number of initial simulation draws to discard.
      real(dp), intent(in), optional :: coefficient_prior_mean !! Coefficient prior mean.
      real(dp), intent(in), optional :: coefficient_prior_variance !! Coefficient prior variance.
      real(dp), intent(in), optional :: residual_variance !! Residual variance.
      real(dp), intent(in), optional :: residual_prior_shape !! Residual prior shape.
      real(dp), intent(in), optional :: residual_prior_rate !! Residual prior rate.
      real(dp), intent(in), optional :: offset_draws(:, :) !! Offset simulation draws.
      type(bsts_holiday_regression_t) :: out
      real(dp), allocatable :: design(:, :), normals(:, :), gammas(:)
      real(dp) :: series_variance, mean0, coefficient_variance0
      real(dp) :: residual0, shape0, rate0
      integer :: variables, discarded, iteration, component, observed

      design = bsts_holiday_design(first_date, size(y), holidays)
      variables = size(design, 2)
      if (size(y) < 1 .or. count(ieee_is_finite(y)) < 1 .or. &
         variables < 1 .or. iterations < 1) then
         out%info = 1
         return
      end if
      discarded = 0
      if (present(burn)) discarded = burn
      if (discarded < 0 .or. discarded >= iterations) then
         out%info = 1
         return
      end if
      series_variance = max(finite_sample_variance(y), 1.0e-6_dp)
      mean0 = 0.0_dp
      if (present(coefficient_prior_mean)) mean0 = coefficient_prior_mean
      coefficient_variance0 = series_variance
      if (present(coefficient_prior_variance)) &
         coefficient_variance0 = coefficient_prior_variance
      residual0 = max(0.5_dp*series_variance, 1.0e-8_dp)
      if (present(residual_variance)) residual0 = residual_variance
      shape0 = 2.0_dp
      if (present(residual_prior_shape)) shape0 = residual_prior_shape
      rate0 = residual0
      if (present(residual_prior_rate)) rate0 = residual_prior_rate
      if (coefficient_variance0 <= 0.0_dp .or. residual0 <= 0.0_dp .or. &
         shape0 <= 0.0_dp .or. rate0 <= 0.0_dp) then
         out%info = 1
         return
      end if
      allocate(normals(variables, iterations), gammas(iterations))
      observed = count(ieee_is_finite(y))
      do iteration = 1, iterations
         do component = 1, variables
            normals(component, iteration) = random_standard_normal()
         end do
         gammas(iteration) = random_gamma(shape0 + &
            0.5_dp*real(observed, dp))
      end do
      if (present(offset_draws)) then
         out = bsts_regression_holiday_draws(y, first_date, holidays, mean0, &
            coefficient_variance0, residual0, shape0, rate0, discarded, &
            normals, gammas, offset_draws)
      else
         out = bsts_regression_holiday_draws(y, first_date, holidays, mean0, &
            coefficient_variance0, residual0, shape0, rate0, discarded, &
            normals, gammas)
      end if
   end function bsts_regression_holiday

   pure function bsts_hierarchical_regression_holiday_draws(y, first_date, &
      holidays, coefficient_mean_prior, coefficient_mean_prior_covariance, &
      coefficient_variance_prior_df, coefficient_variance_prior_scale, &
      residual_variance, residual_prior_shape, residual_prior_rate, burn, &
      coefficient_normal_draws, mean_normal_draws, wishart_normal_draws, &
      wishart_gamma_draws, residual_gamma_draws, offset_draws) result(out)
      !! Sample hierarchical holiday regression from supplied random draws.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: coefficient_mean_prior(:) !! Coefficient mean prior.
      type(date_t), intent(in) :: first_date !! First date.
      type(bsts_holiday_t), intent(in) :: holidays(:) !! Holidays.
      real(dp), intent(in) :: coefficient_mean_prior_covariance(:, :) !! Coefficient mean prior covariance matrix.
      real(dp), intent(in) :: coefficient_variance_prior_df !! Coefficient variance prior degrees of freedom.
      real(dp), intent(in) :: coefficient_variance_prior_scale(:, :) !! Coefficient variance prior scale.
      real(dp), intent(in) :: residual_variance !! Residual variance.
      real(dp), intent(in) :: residual_prior_shape !! Residual prior shape.
      real(dp), intent(in) :: residual_prior_rate !! Residual prior rate.
      integer, intent(in) :: burn !! Number of initial simulation draws to discard.
      real(dp), intent(in) :: coefficient_normal_draws(:, :) !! Coefficient normal simulation draws.
      real(dp), intent(in) :: mean_normal_draws(:, :) !! Mean normal simulation draws.
      real(dp), intent(in) :: wishart_normal_draws(:, :, :) !! Wishart normal simulation draws.
      real(dp), intent(in) :: wishart_gamma_draws(:, :) !! Wishart gamma simulation draws.
      real(dp), intent(in) :: residual_gamma_draws(:) !! Residual gamma simulation draws.
      real(dp), intent(in), optional :: offset_draws(:, :) !! Offset simulation draws.
      type(bsts_holiday_regression_t) :: out
      real(dp), allocatable :: design(:, :), observed_design(:, :)
      real(dp), allocatable :: observed_y(:), working_y(:)
      real(dp), allocatable :: prior_precision(:, :), group_precision(:, :)
      real(dp), allocatable :: precision(:, :), covariance(:, :), score(:)
      real(dp), allocatable :: posterior_mean(:), coefficient(:)
      real(dp), allocatable :: mean_precision(:, :), mean_covariance(:, :)
      real(dp), allocatable :: mean_score(:), mean_current(:)
      real(dp), allocatable :: variance_current(:, :), posterior_scale(:, :)
      real(dp), allocatable :: difference(:)
      integer, allocatable :: observed_index(:)
      real(dp) :: variance, residual_sum
      integer :: observations, width, groups, variables, iterations, observed
      integer :: iteration, group, first, last, component, status

      observations = size(y)
      groups = size(holidays)
      iterations = size(coefficient_normal_draws, 2)
      width = 0
      if (groups > 0) width = bsts_holiday_width(holidays(1))
      design = bsts_holiday_design(first_date, observations, holidays)
      variables = size(design, 2)
      observed = count(ieee_is_finite(y))
      if (observations < 1 .or. observed < 1 .or. groups < 3 .or. width < 1 .or. &
         any([(bsts_holiday_width(holidays(group)) /= width, &
         group=1, groups)]) .or. variables /= groups*width .or. &
         size(design, 1) /= observations .or. &
         size(coefficient_mean_prior) /= width .or. &
         any(shape(coefficient_mean_prior_covariance) /= [width, width]) .or. &
         coefficient_variance_prior_df <= real(width - 1, dp) .or. &
         any(shape(coefficient_variance_prior_scale) /= [width, width]) .or. &
         residual_variance <= 0.0_dp .or. residual_prior_shape <= 0.0_dp .or. &
         residual_prior_rate <= 0.0_dp .or. iterations < 1 .or. &
         burn < 0 .or. burn >= iterations .or. &
         any(shape(coefficient_normal_draws) /= [variables, iterations]) .or. &
         any(shape(mean_normal_draws) /= [width, iterations]) .or. &
         any(shape(wishart_normal_draws) /= [width, width, iterations]) .or. &
         any(shape(wishart_gamma_draws) /= [width, iterations]) .or. &
         size(residual_gamma_draws) /= iterations .or. &
         any(wishart_gamma_draws <= 0.0_dp) .or. &
         any(residual_gamma_draws <= 0.0_dp) .or. &
         .not. all(ieee_is_finite(coefficient_normal_draws)) .or. &
         .not. all(ieee_is_finite(mean_normal_draws)) .or. &
         .not. all(ieee_is_finite(wishart_normal_draws))) then
         out%info = 1
         return
      end if
      if (present(offset_draws)) then
         if (any(shape(offset_draws) /= [observations, iterations]) .or. &
            .not. all(ieee_is_finite(offset_draws))) then
            out%info = 1
            return
         end if
      end if
      call invert_matrix(coefficient_mean_prior_covariance, prior_precision, &
         status)
      if (status /= 0) then
         out%info = 2
         return
      end if
      call invert_matrix(coefficient_variance_prior_scale, group_precision, &
         status)
      if (status /= 0) then
         out%info = 3
         return
      end if
      observed_index = pack([(component, component=1, observations)], &
         ieee_is_finite(y))
      allocate(observed_design(observed, variables), observed_y(observed), &
         working_y(observed), precision(variables, variables), &
         score(variables), posterior_mean(variables), coefficient(variables), &
         mean_precision(width, width), mean_score(width), mean_current(width), &
         variance_current(width, width), posterior_scale(width, width), &
         difference(width), out%coefficients(variables, iterations), &
         out%residual_variance(iterations), &
         out%contribution(observations, iterations), &
         out%coefficient_mean(width, iterations), &
         out%coefficient_variance(width, width, iterations), &
         out%coefficient_offset(groups + 1), out%holidays(groups))
      observed_design = design(observed_index, :)
      observed_y = y(observed_index)
      out%holidays = holidays
      do group = 1, groups + 1
         out%coefficient_offset(group) = (group - 1)*width
      end do
      out%first_date = first_date
      out%last_date = date_from_day_number(date_day_number(first_date) + &
         observations - 1)
      out%burn = burn
      out%hierarchical = .true.
      mean_current = coefficient_mean_prior
      variance_current = coefficient_variance_prior_scale/ &
         max(coefficient_variance_prior_df - real(width + 1, dp), 1.0_dp)
      variance = residual_variance
      do iteration = 1, iterations
         working_y = observed_y
         if (present(offset_draws)) working_y = working_y - &
            offset_draws(observed_index, iteration)
         call invert_matrix(variance_current, group_precision, status)
         if (status /= 0) then
            out%info = 10 + status
            return
         end if
         precision = matmul(transpose(observed_design), observed_design)/variance
         score = matmul(transpose(observed_design), working_y)/variance
         do group = 1, groups
            first = (group - 1)*width + 1
            last = group*width
            precision(first:last, first:last) = &
               precision(first:last, first:last) + group_precision
            score(first:last) = score(first:last) + &
               matmul(group_precision, mean_current)
         end do
         call invert_matrix(precision, covariance, status)
         if (status /= 0) then
            out%info = 20 + status
            return
         end if
         posterior_mean = matmul(covariance, score)
         call multivariate_normal_from_standard(posterior_mean, covariance, &
            coefficient_normal_draws(:, iteration), coefficient, status)
         if (status /= 0) then
            out%info = 30 + status
            return
         end if
         mean_precision = prior_precision + real(groups, dp)*group_precision
         call invert_matrix(mean_precision, mean_covariance, status)
         if (status /= 0) then
            out%info = 40 + status
            return
         end if
         mean_score = matmul(prior_precision, coefficient_mean_prior)
         do group = 1, groups
            first = (group - 1)*width + 1
            last = group*width
            mean_score = mean_score + &
               matmul(group_precision, coefficient(first:last))
         end do
         posterior_mean(:width) = matmul(mean_covariance, mean_score)
         call multivariate_normal_from_standard(posterior_mean(:width), &
            mean_covariance, mean_normal_draws(:, iteration), mean_current, &
            status)
         if (status /= 0) then
            out%info = 50 + status
            return
         end if
         posterior_scale = coefficient_variance_prior_scale
         do group = 1, groups
            first = (group - 1)*width + 1
            last = group*width
            difference = coefficient(first:last) - mean_current
            posterior_scale = posterior_scale + &
               spread(difference, 2, width)*spread(difference, 1, width)
         end do
         call inverse_wishart_from_draws( &
            coefficient_variance_prior_df + real(groups, dp), posterior_scale, &
            wishart_normal_draws(:, :, iteration), &
            wishart_gamma_draws(:, iteration), variance_current, status)
         if (status /= 0) then
            out%info = 60 + status
            return
         end if
         residual_sum = sum((working_y - matmul(observed_design, coefficient))**2)
         variance = (residual_prior_rate + 0.5_dp*residual_sum)/ &
            residual_gamma_draws(iteration)
         out%coefficients(:, iteration) = coefficient
         out%coefficient_mean(:, iteration) = mean_current
         out%coefficient_variance(:, :, iteration) = variance_current
         out%residual_variance(iteration) = variance
         out%contribution(:, iteration) = matmul(design, coefficient)
      end do
   end function bsts_hierarchical_regression_holiday_draws

   function bsts_hierarchical_regression_holiday(y, first_date, holidays, &
      iterations, burn, coefficient_mean_prior, &
      coefficient_mean_prior_covariance, coefficient_variance_prior_df, &
      coefficient_variance_prior_scale, residual_variance, &
      residual_prior_shape, residual_prior_rate, offset_draws) result(out)
      !! Sample hierarchical holiday regression using the random stream.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      type(date_t), intent(in) :: first_date !! First date.
      type(bsts_holiday_t), intent(in) :: holidays(:) !! Holidays.
      integer, intent(in) :: iterations !! Number of algorithm iterations.
      integer, intent(in), optional :: burn !! Number of initial simulation draws to discard.
      real(dp), intent(in), optional :: coefficient_mean_prior(:) !! Coefficient mean prior.
      real(dp), intent(in), optional :: coefficient_mean_prior_covariance(:, :) !! Coefficient mean prior covariance matrix.
      real(dp), intent(in), optional :: coefficient_variance_prior_df !! Coefficient variance prior degrees of freedom.
      real(dp), intent(in), optional :: coefficient_variance_prior_scale(:, :) !! Coefficient variance prior scale.
      real(dp), intent(in), optional :: residual_variance !! Residual variance.
      real(dp), intent(in), optional :: residual_prior_shape !! Residual prior shape.
      real(dp), intent(in), optional :: residual_prior_rate !! Residual prior rate.
      real(dp), intent(in), optional :: offset_draws(:, :) !! Offset simulation draws.
      type(bsts_holiday_regression_t) :: out
      real(dp), allocatable :: mean0(:), mean_covariance0(:, :), scale0(:, :)
      real(dp), allocatable :: coefficient_normals(:, :), mean_normals(:, :)
      real(dp), allocatable :: wishart_normals(:, :, :), wishart_gammas(:, :)
      real(dp), allocatable :: residual_gammas(:)
      real(dp) :: series_variance, degrees0, residual0, shape0, rate0
      real(dp) :: posterior_degrees
      integer :: width, groups, variables, discarded, observed
      integer :: iteration, component, row, column

      groups = size(holidays)
      width = 0
      if (groups > 0) width = bsts_holiday_width(holidays(1))
      variables = groups*width
      if (size(y) < 1 .or. count(ieee_is_finite(y)) < 1 .or. groups < 3 .or. &
         width < 1 .or. any([(bsts_holiday_width(holidays(component)) /= width, &
         component=1, groups)]) .or. iterations < 1) then
         out%info = 1
         return
      end if
      discarded = 0
      if (present(burn)) discarded = burn
      if (discarded < 0 .or. discarded >= iterations) then
         out%info = 1
         return
      end if
      series_variance = max(finite_sample_variance(y), 1.0e-6_dp)
      allocate(mean0(width), mean_covariance0(width, width), &
         scale0(width, width))
      mean0 = 0.0_dp
      if (present(coefficient_mean_prior)) then
         if (size(coefficient_mean_prior) /= width) then
            out%info = 1
            return
         end if
         mean0 = coefficient_mean_prior
      end if
      mean_covariance0 = identity_block(width, series_variance)
      if (present(coefficient_mean_prior_covariance)) then
         if (any(shape(coefficient_mean_prior_covariance) /= [width, width])) then
            out%info = 1
            return
         end if
         mean_covariance0 = coefficient_mean_prior_covariance
      end if
      degrees0 = real(width + 2, dp)
      if (present(coefficient_variance_prior_df)) &
         degrees0 = coefficient_variance_prior_df
      scale0 = identity_block(width, 0.1_dp*series_variance)
      if (present(coefficient_variance_prior_scale)) then
         if (any(shape(coefficient_variance_prior_scale) /= [width, width])) then
            out%info = 1
            return
         end if
         scale0 = coefficient_variance_prior_scale
      end if
      residual0 = max(0.5_dp*series_variance, 1.0e-8_dp)
      if (present(residual_variance)) residual0 = residual_variance
      shape0 = 2.0_dp
      if (present(residual_prior_shape)) shape0 = residual_prior_shape
      rate0 = residual0
      if (present(residual_prior_rate)) rate0 = residual_prior_rate
      if (degrees0 <= real(width - 1, dp) .or. residual0 <= 0.0_dp .or. &
         shape0 <= 0.0_dp .or. rate0 <= 0.0_dp) then
         out%info = 1
         return
      end if
      allocate(coefficient_normals(variables, iterations), &
         mean_normals(width, iterations), &
         wishart_normals(width, width, iterations), &
         wishart_gammas(width, iterations), residual_gammas(iterations))
      observed = count(ieee_is_finite(y))
      posterior_degrees = degrees0 + real(groups, dp)
      do iteration = 1, iterations
         do component = 1, variables
            coefficient_normals(component, iteration) = random_standard_normal()
         end do
         do row = 1, width
            mean_normals(row, iteration) = random_standard_normal()
            wishart_gammas(row, iteration) = random_gamma( &
               0.5_dp*(posterior_degrees - real(row, dp) + 1.0_dp))
            do column = 1, width
               wishart_normals(row, column, iteration) = &
                  random_standard_normal()
            end do
         end do
         residual_gammas(iteration) = random_gamma(shape0 + &
            0.5_dp*real(observed, dp))
      end do
      if (present(offset_draws)) then
         out = bsts_hierarchical_regression_holiday_draws(y, first_date, &
            holidays, mean0, mean_covariance0, degrees0, scale0, residual0, &
            shape0, rate0, discarded, coefficient_normals, mean_normals, &
            wishart_normals, wishart_gammas, residual_gammas, offset_draws)
      else
         out = bsts_hierarchical_regression_holiday_draws(y, first_date, &
            holidays, mean0, mean_covariance0, degrees0, scale0, residual0, &
            shape0, rate0, discarded, coefficient_normals, mean_normals, &
            wishart_normals, wishart_gammas, residual_gammas)
      end if
   end function bsts_hierarchical_regression_holiday

   pure function bsts_regression_holiday_predict_draws(fit, horizon, &
      observation_normal_draws) result(out)
      !! Forecast fixed holiday effects from supplied observation normals.
      type(bsts_holiday_regression_t), intent(in) :: fit !! Previously fitted model.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      real(dp), intent(in) :: observation_normal_draws(:, :) !! Observation normal draws.
      type(bsts_prediction_t) :: out
      real(dp), allocatable :: design(:, :), sorted(:)
      type(date_t) :: first_future
      integer :: retained, draw, source, step

      if (fit%info /= 0 .or. .not. allocated(fit%coefficients) .or. &
         .not. allocated(fit%residual_variance) .or. &
         .not. allocated(fit%holidays) .or. .not. date_valid(fit%last_date)) then
         out%info = 1
         return
      end if
      retained = size(fit%coefficients, 2) - fit%burn
      if (horizon < 1 .or. retained < 1 .or. &
         any(shape(observation_normal_draws) /= [horizon, retained]) .or. &
         .not. all(ieee_is_finite(observation_normal_draws))) then
         out%info = 1
         return
      end if
      first_future = date_from_day_number(date_day_number(fit%last_date) + 1)
      design = bsts_holiday_design(first_future, horizon, fit%holidays)
      if (size(design, 2) /= size(fit%coefficients, 1)) then
         out%info = 1
         return
      end if
      allocate(out%draws(horizon, retained), out%mean(horizon), &
         out%standard_deviation(horizon), out%lower(horizon), out%upper(horizon))
      do draw = 1, retained
         source = fit%burn + draw
         out%draws(:, draw) = matmul(design, fit%coefficients(:, source)) + &
            sqrt(fit%residual_variance(source))*observation_normal_draws(:, draw)
      end do
      do step = 1, horizon
         out%mean(step) = sum(out%draws(step, :))/real(retained, dp)
         if (retained > 1) then
            out%standard_deviation(step) = sqrt(sum((out%draws(step, :) - &
               out%mean(step))**2)/real(retained - 1, dp))
         else
            out%standard_deviation(step) = 0.0_dp
         end if
         sorted = out%draws(step, :)
         call insertion_sort(sorted)
         out%lower(step) = quantile(sorted, 0.025_dp)
         out%upper(step) = quantile(sorted, 0.975_dp)
      end do
   end function bsts_regression_holiday_predict_draws

   function bsts_regression_holiday_predict(fit, horizon) result(out)
      !! Forecast fixed holiday effects using the random stream.
      type(bsts_holiday_regression_t), intent(in) :: fit !! Previously fitted model.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      type(bsts_prediction_t) :: out
      real(dp), allocatable :: normals(:, :)
      integer :: retained, draw, step

      if (.not. allocated(fit%coefficients)) then
         out%info = 1
         return
      end if
      retained = size(fit%coefficients, 2) - fit%burn
      if (horizon < 1 .or. retained < 1) then
         out%info = 1
         return
      end if
      allocate(normals(horizon, retained))
      do draw = 1, retained
         do step = 1, horizon
            normals(step, draw) = random_standard_normal()
         end do
      end do
      out = bsts_regression_holiday_predict_draws(fit, horizon, normals)
   end function bsts_regression_holiday_predict

   pure function bsts_shared_local_level_draws(response, initial_mean, &
      initial_covariance, initial_loadings, loading_prior_mean, &
      loading_prior_variance, observation_variance, factor_variance, &
      observation_prior_shape, observation_prior_rate, factor_prior_shape, &
      factor_prior_rate, burn, state_normal_draws, loading_normal_draws, &
      gamma_draws, offset_draws, prior_inclusion_probability, &
      initial_inclusion, inclusion_uniform_draws, maximum_flips, &
      regression_predictors, regression_prior_mean, &
      regression_prior_covariance, regression_prior_inclusion_probability, &
      initial_regression_inclusion, regression_normal_draws, &
      regression_uniform_draws, regression_maximum_model_size, &
      regression_maximum_flips, series_initial_mean, &
      series_initial_variance, series_variance, series_prior_shape, &
      series_prior_rate, series_state_normal_draws, &
      series_gamma_draws, trend_initial_mean, trend_initial_covariance, &
      trend_variance, trend_prior_shape, trend_prior_rate, &
      trend_state_normal_draws, trend_gamma_draws, seasonal_nseasons, &
      seasonal_duration, seasonal_initial_variance, seasonal_variance, &
      seasonal_prior_shape, seasonal_prior_rate, &
      seasonal_state_normal_draws, seasonal_gamma_draws, &
      enable_series_state, enable_trend_state, &
      enable_seasonal_state) result(out)
      !! Sample a shared multivariate local-level model from supplied draws.
      real(dp), intent(in) :: response(:, :) !! Response observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: initial_loadings(:, :) !! Initial loadings.
      real(dp), intent(in) :: loading_prior_mean(:, :) !! Loading prior mean.
      real(dp), intent(in) :: loading_prior_variance(:, :) !! Loading prior variance.
      real(dp), intent(in) :: observation_variance(:) !! Observation-error variance.
      real(dp), intent(in) :: factor_variance(:) !! Factor variance.
      real(dp), intent(in) :: observation_prior_shape(:) !! Observation prior shape.
      real(dp), intent(in) :: observation_prior_rate(:) !! Observation prior rate.
      real(dp), intent(in) :: factor_prior_shape(:) !! Factor prior shape.
      real(dp), intent(in) :: factor_prior_rate(:) !! Factor prior rate.
      integer, intent(in) :: burn !! Number of initial simulation draws to discard.
      real(dp), intent(in) :: state_normal_draws(:, :, :) !! Independent standard-normal state draws.
      real(dp), intent(in) :: loading_normal_draws(:, :, :) !! Loading normal simulation draws.
      real(dp), intent(in) :: gamma_draws(:, :) !! Gamma simulation draws.
      real(dp), intent(in), optional :: offset_draws(:, :, :) !! Offset simulation draws.
      real(dp), intent(in), optional :: prior_inclusion_probability(:, :) !! Prior inclusion probability.
      logical, intent(in), optional :: initial_inclusion(:, :) !! Initial inclusion.
      real(dp), intent(in), optional :: inclusion_uniform_draws(:, :, :) !! Inclusion uniform simulation draws.
      integer, intent(in), optional :: maximum_flips !! Maximum flips.
      real(dp), intent(in), optional :: regression_predictors(:, :, :) !! Regression predictors.
      real(dp), intent(in), optional :: regression_prior_mean(:, :) !! Regression prior mean.
      real(dp), intent(in), optional :: regression_prior_covariance(:, :, :) !! Regression prior covariance matrix.
      real(dp), intent(in), optional :: regression_prior_inclusion_probability(:, :) !! Regression prior inclusion probability.
      logical, intent(in), optional :: initial_regression_inclusion(:, :) !! Initial regression inclusion.
      real(dp), intent(in), optional :: regression_normal_draws(:, :, :) !! Regression normal simulation draws.
      real(dp), intent(in), optional :: regression_uniform_draws(:, :, :) !! Regression uniform simulation draws.
      integer, intent(in), optional :: regression_maximum_model_size !! Regression maximum model size.
      integer, intent(in), optional :: regression_maximum_flips !! Regression maximum flips.
      real(dp), intent(in), optional :: series_initial_mean(:) !! Series initial mean.
      real(dp), intent(in), optional :: series_initial_variance(:) !! Series initial variance.
      real(dp), intent(in), optional :: series_variance(:) !! Series variance.
      real(dp), intent(in), optional :: series_prior_shape(:) !! Series prior shape.
      real(dp), intent(in), optional :: series_prior_rate(:) !! Series prior rate.
      real(dp), intent(in), optional :: series_state_normal_draws(:, :, :) !! Series state normal simulation draws.
      real(dp), intent(in), optional :: series_gamma_draws(:, :) !! Series gamma simulation draws.
      real(dp), intent(in), optional :: trend_initial_mean(:, :) !! Trend initial mean.
      real(dp), intent(in), optional :: trend_initial_covariance(:, :, :) !! Trend initial covariance matrix.
      real(dp), intent(in), optional :: trend_variance(:, :) !! Trend variance.
      real(dp), intent(in), optional :: trend_prior_shape(:, :) !! Trend prior shape.
      real(dp), intent(in), optional :: trend_prior_rate(:, :) !! Trend prior rate.
      real(dp), intent(in), optional :: trend_state_normal_draws(:, :, :, :) !! Trend state normal simulation draws.
      real(dp), intent(in), optional :: trend_gamma_draws(:, :, :) !! Trend gamma simulation draws.
      integer, intent(in), optional :: seasonal_nseasons !! Seasonal nseasons.
      integer, intent(in), optional :: seasonal_duration !! Seasonal duration.
      real(dp), intent(in), optional :: seasonal_initial_variance(:) !! Seasonal initial variance.
      real(dp), intent(in), optional :: seasonal_variance(:) !! Seasonal variance.
      real(dp), intent(in), optional :: seasonal_prior_shape(:) !! Seasonal prior shape.
      real(dp), intent(in), optional :: seasonal_prior_rate(:) !! Seasonal prior rate.
      real(dp), intent(in), optional :: seasonal_state_normal_draws(:, :, :, :) !! Seasonal state normal simulation draws.
      real(dp), intent(in), optional :: seasonal_gamma_draws(:, :) !! Seasonal gamma simulation draws.
      logical, intent(in), optional :: enable_series_state !! Flag controlling enable series state.
      logical, intent(in), optional :: enable_trend_state !! Flag controlling enable trend state.
      logical, intent(in), optional :: enable_seasonal_state !! Flag controlling enable seasonal state.
      type(bsts_shared_local_level_t) :: out
      type(ssm_model_t) :: model
      real(dp), allocatable :: path(:, :), loadings(:, :)
      real(dp), allocatable :: observation_scale(:), factor_scale(:)
      real(dp), allocatable :: predictors(:, :), working_response(:)
      real(dp), allocatable :: precision(:, :), covariance(:, :), score(:)
      real(dp), allocatable :: posterior_mean(:), coefficient(:)
      real(dp), allocatable :: slab_covariance(:, :), active_draw(:)
      integer, allocatable :: observed_index(:)
      integer, allocatable :: active(:)
      logical, allocatable :: included(:, :), candidate(:)
      integer, allocatable :: regression_active(:)
      logical, allocatable :: regression_included(:, :)
      logical, allocatable :: regression_candidate(:)
      real(dp), allocatable :: regression_draw(:)
      real(dp), allocatable :: series_scale(:)
      real(dp), allocatable :: combined_state_draws(:, :)
      real(dp), allocatable :: trend_scale(:, :), seasonal_scale(:)
      real(dp) :: fixed_part, residual_sum, difference_sum
      real(dp) :: log_included, log_excluded, log_probability, collapsed_rate
      real(dp) :: probability
      integer :: observations, series_count, factors, iterations
      integer :: iteration, series, factor, time, free_count, observed, status
      integer :: flips, flip_limit, retained
      integer :: dimension, base_dimension, predictors_count
      integer :: coefficient_start, coefficient_end, series_state_start
      integer :: trend_state_start, seasonal_state_start, seasonal_dimension
      integer :: transition_slices, slice, component
      integer :: seasonal_duration_value
      integer :: regression_model_limit, regression_flip_limit
      logical :: use_selection, use_regression, use_regression_selection
      logical :: use_series_state
      logical :: use_trend_state, use_seasonal_state

      observations = size(response, 1)
      series_count = size(response, 2)
      factors = size(initial_mean)
      iterations = size(state_normal_draws, 3)
      use_selection = present(prior_inclusion_probability) .or. &
         present(initial_inclusion) .or. present(inclusion_uniform_draws)
      use_regression = present(regression_predictors) .or. &
         present(regression_prior_mean) .or. &
         present(regression_prior_covariance)
      use_regression_selection = &
         present(regression_prior_inclusion_probability) .or. &
         present(initial_regression_inclusion) .or. &
         present(regression_normal_draws) .or. &
         present(regression_uniform_draws)
      use_series_state = present(series_initial_mean) .or. &
         present(series_initial_variance) .or. present(series_variance) .or. &
         present(series_prior_shape) .or. present(series_prior_rate) .or. &
         present(series_state_normal_draws) .or. present(series_gamma_draws)
      use_trend_state = present(trend_initial_mean) .or. &
         present(trend_initial_covariance) .or. present(trend_variance) .or. &
         present(trend_prior_shape) .or. present(trend_prior_rate) .or. &
         present(trend_state_normal_draws) .or. present(trend_gamma_draws)
      use_seasonal_state = present(seasonal_nseasons) .or. &
         present(seasonal_duration) .or. &
         present(seasonal_initial_variance) .or. &
         present(seasonal_variance) .or. present(seasonal_prior_shape) .or. &
         present(seasonal_prior_rate) .or. &
         present(seasonal_state_normal_draws) .or. &
         present(seasonal_gamma_draws)
      if (present(enable_series_state)) use_series_state = enable_series_state
      if (present(enable_trend_state)) use_trend_state = enable_trend_state
      if (present(enable_seasonal_state)) &
         use_seasonal_state = enable_seasonal_state
      seasonal_duration_value = 1
      if (present(seasonal_duration)) &
         seasonal_duration_value = seasonal_duration
      predictors_count = 0
      if (present(regression_predictors)) &
         predictors_count = size(regression_predictors, 2)
      base_dimension = factors + series_count*predictors_count
      dimension = base_dimension
      if (use_series_state) dimension = dimension + series_count
      series_state_start = base_dimension
      trend_state_start = dimension
      if (use_trend_state) dimension = dimension + 2*series_count
      seasonal_state_start = dimension
      seasonal_dimension = 0
      if (use_seasonal_state .and. present(seasonal_nseasons)) &
         seasonal_dimension = seasonal_nseasons - 1
      if (use_seasonal_state) &
         dimension = dimension + seasonal_dimension*series_count
      transition_slices = 1
      if (use_seasonal_state) transition_slices = observations
      regression_model_limit = predictors_count
      if (present(regression_maximum_model_size)) &
         regression_model_limit = regression_maximum_model_size
      regression_flip_limit = -1
      if (present(regression_maximum_flips)) &
         regression_flip_limit = regression_maximum_flips
      flip_limit = -1
      if (present(maximum_flips)) flip_limit = maximum_flips
      if (observations < 2 .or. series_count < 1 .or. factors < 1 .or. &
         factors > series_count .or. &
         any([(count(ieee_is_finite(response(:, series))) < 1, &
         series=1, series_count)]) .or. &
         any(shape(initial_covariance) /= [factors, factors]) .or. &
         any(shape(initial_loadings) /= [series_count, factors]) .or. &
         any(shape(loading_prior_mean) /= [series_count, factors]) .or. &
         any(shape(loading_prior_variance) /= [series_count, factors]) .or. &
         size(observation_variance) /= series_count .or. &
         size(factor_variance) /= factors .or. &
         size(observation_prior_shape) /= series_count .or. &
         size(observation_prior_rate) /= series_count .or. &
         size(factor_prior_shape) /= factors .or. &
         size(factor_prior_rate) /= factors .or. &
         any(observation_variance <= 0.0_dp) .or. &
         any(factor_variance <= 0.0_dp) .or. &
         any(observation_prior_shape <= 0.0_dp) .or. &
         any(observation_prior_rate <= 0.0_dp) .or. &
         any(factor_prior_shape <= 0.0_dp) .or. &
         any(factor_prior_rate <= 0.0_dp) .or. iterations < 1 .or. &
         burn < 0 .or. burn >= iterations .or. &
         any(shape(state_normal_draws) /= &
         [base_dimension, observations, iterations]) .or. &
         any(shape(loading_normal_draws) /= &
         [factors, series_count, iterations]) .or. &
         any(shape(gamma_draws) /= [series_count + factors, iterations]) .or. &
         any(gamma_draws <= 0.0_dp) .or. &
         .not. all(ieee_is_finite(initial_mean)) .or. &
         .not. all(ieee_is_finite(initial_covariance)) .or. &
         .not. all(ieee_is_finite(initial_loadings)) .or. &
         .not. all(ieee_is_finite(loading_prior_mean)) .or. &
         .not. all(ieee_is_finite(loading_prior_variance)) .or. &
         .not. all(ieee_is_finite(observation_variance)) .or. &
         .not. all(ieee_is_finite(factor_variance)) .or. &
         .not. all(ieee_is_finite(observation_prior_shape)) .or. &
         .not. all(ieee_is_finite(observation_prior_rate)) .or. &
         .not. all(ieee_is_finite(factor_prior_shape)) .or. &
         .not. all(ieee_is_finite(factor_prior_rate)) .or. &
         .not. all(ieee_is_finite(state_normal_draws)) .or. &
         .not. all(ieee_is_finite(loading_normal_draws)) .or. &
         .not. all(ieee_is_finite(gamma_draws)) .or. flip_limit < -1) then
         out%info = 1
         return
      end if
      if (use_series_state) then
         if (.not. use_regression .or. &
            .not. present(series_initial_mean) .or. &
            .not. present(series_initial_variance) .or. &
            .not. present(series_variance) .or. &
            .not. present(series_prior_shape) .or. &
            .not. present(series_prior_rate) .or. &
            .not. present(series_state_normal_draws) .or. &
            .not. present(series_gamma_draws)) then
            out%info = 1
            return
         end if
         if (size(series_initial_mean) /= series_count .or. &
            size(series_initial_variance) /= series_count .or. &
            size(series_variance) /= series_count .or. &
            size(series_prior_shape) /= series_count .or. &
            size(series_prior_rate) /= series_count .or. &
            any(shape(series_state_normal_draws) /= &
            [series_count, observations, iterations]) .or. &
            any(shape(series_gamma_draws) /= [series_count, iterations]) .or. &
            any(series_initial_variance <= 0.0_dp) .or. &
            any(series_variance <= 0.0_dp) .or. &
            any(series_prior_shape <= 0.0_dp) .or. &
            any(series_prior_rate <= 0.0_dp) .or. &
            any(series_gamma_draws <= 0.0_dp) .or. &
            .not. all(ieee_is_finite(series_initial_mean)) .or. &
            .not. all(ieee_is_finite(series_initial_variance)) .or. &
            .not. all(ieee_is_finite(series_variance)) .or. &
            .not. all(ieee_is_finite(series_prior_shape)) .or. &
            .not. all(ieee_is_finite(series_prior_rate)) .or. &
            .not. all(ieee_is_finite(series_state_normal_draws)) .or. &
            .not. all(ieee_is_finite(series_gamma_draws))) then
            out%info = 1
            return
         end if
      end if
      if (use_trend_state) then
         if (.not. use_regression .or. .not. present(trend_initial_mean) .or. &
            .not. present(trend_initial_covariance) .or. &
            .not. present(trend_variance) .or. &
            .not. present(trend_prior_shape) .or. &
            .not. present(trend_prior_rate) .or. &
            .not. present(trend_state_normal_draws) .or. &
            .not. present(trend_gamma_draws)) then
            out%info = 1
            return
         end if
         if (any(shape(trend_initial_mean) /= [2, series_count]) .or. &
            any(shape(trend_initial_covariance) /= [2, 2, series_count]) .or. &
            any(shape(trend_variance) /= [2, series_count]) .or. &
            any(shape(trend_prior_shape) /= [2, series_count]) .or. &
            any(shape(trend_prior_rate) /= [2, series_count]) .or. &
            any(shape(trend_state_normal_draws) /= &
            [2, series_count, observations, iterations]) .or. &
            any(shape(trend_gamma_draws) /= [2, series_count, iterations]) .or. &
            any(trend_variance <= 0.0_dp) .or. &
            any(trend_prior_shape <= 0.0_dp) .or. &
            any(trend_prior_rate <= 0.0_dp) .or. &
            any(trend_gamma_draws <= 0.0_dp) .or. &
            .not. all(ieee_is_finite(trend_initial_mean)) .or. &
            .not. all(ieee_is_finite(trend_initial_covariance)) .or. &
            .not. all(ieee_is_finite(trend_variance)) .or. &
            .not. all(ieee_is_finite(trend_prior_shape)) .or. &
            .not. all(ieee_is_finite(trend_prior_rate)) .or. &
            .not. all(ieee_is_finite(trend_state_normal_draws)) .or. &
            .not. all(ieee_is_finite(trend_gamma_draws))) then
            out%info = 1
            return
         end if
         do series = 1, series_count
            call invert_matrix(trend_initial_covariance(:, :, series), &
               covariance, status)
            if (status /= 0) then
               out%info = 5
               return
            end if
         end do
      end if
      if (use_seasonal_state) then
         if (.not. use_regression .or. .not. present(seasonal_nseasons) .or. &
            .not. present(seasonal_duration) .or. &
            .not. present(seasonal_initial_variance) .or. &
            .not. present(seasonal_variance) .or. &
            .not. present(seasonal_prior_shape) .or. &
            .not. present(seasonal_prior_rate) .or. &
            .not. present(seasonal_state_normal_draws) .or. &
            .not. present(seasonal_gamma_draws)) then
            out%info = 1
            return
         end if
         if (seasonal_dimension < 1 .or. seasonal_duration < 1 .or. &
            size(seasonal_initial_variance) /= series_count .or. &
            size(seasonal_variance) /= series_count .or. &
            size(seasonal_prior_shape) /= series_count .or. &
            size(seasonal_prior_rate) /= series_count .or. &
            any(shape(seasonal_state_normal_draws) /= &
            [seasonal_dimension, series_count, observations, iterations]) .or. &
            any(shape(seasonal_gamma_draws) /= [series_count, iterations]) .or. &
            any(seasonal_initial_variance <= 0.0_dp) .or. &
            any(seasonal_variance <= 0.0_dp) .or. &
            any(seasonal_prior_shape <= 0.0_dp) .or. &
            any(seasonal_prior_rate <= 0.0_dp) .or. &
            any(seasonal_gamma_draws <= 0.0_dp) .or. &
            .not. all(ieee_is_finite(seasonal_state_normal_draws))) then
            out%info = 1
            return
         end if
      end if
      if (use_selection) then
         if (.not. present(prior_inclusion_probability) .or. &
            .not. present(initial_inclusion) .or. &
            .not. present(inclusion_uniform_draws)) then
            out%info = 1
            return
         end if
         if (any(shape(prior_inclusion_probability) /= &
            [series_count, factors]) .or. &
            any(shape(initial_inclusion) /= [series_count, factors]) .or. &
            any(shape(inclusion_uniform_draws) /= &
            [factors, series_count, iterations]) .or. &
            any(prior_inclusion_probability < 0.0_dp) .or. &
            any(prior_inclusion_probability > 1.0_dp) .or. &
            any(inclusion_uniform_draws <= 0.0_dp) .or. &
            any(inclusion_uniform_draws >= 1.0_dp) .or. &
            .not. all(ieee_is_finite(prior_inclusion_probability)) .or. &
            .not. all(ieee_is_finite(inclusion_uniform_draws))) then
            out%info = 1
            return
         end if
      end if
      if (use_regression) then
         if (.not. present(regression_predictors) .or. &
            .not. present(regression_prior_mean) .or. &
            .not. present(regression_prior_covariance)) then
            out%info = 1
            return
         end if
         if (predictors_count < 1 .or. &
            any(shape(regression_predictors) /= &
            [observations, predictors_count, series_count]) .or. &
            any(shape(regression_prior_mean) /= &
            [predictors_count, series_count]) .or. &
            any(shape(regression_prior_covariance) /= &
            [predictors_count, predictors_count, series_count]) .or. &
            .not. all(ieee_is_finite(regression_predictors)) .or. &
            .not. all(ieee_is_finite(regression_prior_mean)) .or. &
            .not. all(ieee_is_finite(regression_prior_covariance))) then
            out%info = 1
            return
         end if
         do series = 1, series_count
            call invert_matrix(regression_prior_covariance(:, :, series), &
               covariance, status)
            if (status /= 0) then
               out%info = 3
               return
            end if
         end do
      end if
      if (use_regression_selection) then
         if (.not. use_regression .or. &
            .not. present(regression_prior_inclusion_probability) .or. &
            .not. present(initial_regression_inclusion) .or. &
            .not. present(regression_normal_draws) .or. &
            .not. present(regression_uniform_draws)) then
            out%info = 1
            return
         end if
         if (any(shape(regression_prior_inclusion_probability) /= &
            [predictors_count, series_count]) .or. &
            any(shape(initial_regression_inclusion) /= &
            [predictors_count, series_count]) .or. &
            any(shape(regression_normal_draws) /= &
            [predictors_count, series_count, iterations]) .or. &
            any(shape(regression_uniform_draws) /= &
            [predictors_count, series_count, iterations]) .or. &
            any(regression_prior_inclusion_probability < 0.0_dp) .or. &
            any(regression_prior_inclusion_probability > 1.0_dp) .or. &
            any(regression_uniform_draws <= 0.0_dp) .or. &
            any(regression_uniform_draws >= 1.0_dp) .or. &
            .not. all(ieee_is_finite( &
            regression_prior_inclusion_probability)) .or. &
            .not. all(ieee_is_finite(regression_normal_draws)) .or. &
            .not. all(ieee_is_finite(regression_uniform_draws)) .or. &
            regression_model_limit < 0 .or. &
            regression_model_limit > predictors_count .or. &
            regression_flip_limit < -1) then
            out%info = 1
            return
         end if
      end if
      do series = 1, series_count
         free_count = min(factors, series - 1)
         if (free_count > 0) then
            if (any(loading_prior_variance(series, :free_count) <= 0.0_dp)) then
               out%info = 1
               return
            end if
         end if
      end do
      if (present(offset_draws)) then
         if (any(shape(offset_draws) /= &
            [observations, series_count, iterations]) .or. &
            .not. all(ieee_is_finite(offset_draws))) then
            out%info = 1
            return
         end if
      end if
      call invert_matrix(initial_covariance, covariance, status)
      if (status /= 0) then
         out%info = 2
         return
      end if
      allocate(model%y(observations, series_count), &
         model%z(series_count, dimension, observations), &
         model%h(series_count, series_count, 1), &
         model%transition(dimension, dimension, transition_slices), &
         model%r(dimension, dimension, transition_slices), &
         model%q(dimension, dimension, transition_slices), &
         model%a1(dimension), model%p1(dimension, dimension), &
         model%p1inf(dimension, dimension), &
         model%missing(observations, series_count), &
         loadings(series_count, factors), observation_scale(series_count), &
         factor_scale(factors), &
         out%state(factors, observations, iterations), &
         out%loadings(series_count, factors, iterations), &
         out%observation_variance(series_count, iterations), &
         out%factor_variance(factors, iterations), &
         out%contribution(series_count, observations, iterations), &
         out%included(series_count, factors, iterations), &
         out%inclusion_probability(series_count, factors), &
         included(series_count, factors))
      if (use_regression) then
         allocate(out%regression_coefficients(predictors_count, series_count, &
            iterations), &
            out%regression_contribution(series_count, observations, iterations), &
            out%fitted(series_count, observations, iterations), &
            out%residuals(series_count, observations, iterations))
      end if
      if (use_regression_selection) then
         allocate(regression_included(predictors_count, series_count), &
            out%regression_included(predictors_count, series_count, iterations), &
            out%regression_inclusion_probability(predictors_count, series_count))
         regression_included = initial_regression_inclusion
         where (regression_prior_inclusion_probability <= 0.0_dp) &
            regression_included = .false.
         where (regression_prior_inclusion_probability >= 1.0_dp) &
            regression_included = .true.
         if (any([(count(regression_included(:, series)) > &
            regression_model_limit, series=1, series_count)])) then
            out%info = 4
            return
         end if
      end if
      if (use_series_state) then
         allocate(series_scale(series_count), &
            out%series_state(series_count, observations, iterations), &
            out%series_variance(series_count, iterations), &
            out%series_contribution(series_count, observations, iterations))
         series_scale = series_variance
      end if
      if (use_trend_state) then
         allocate(trend_scale(2, series_count), &
            out%series_trend_state(2, series_count, observations, iterations), &
            out%series_trend_variance(2, series_count, iterations), &
            out%series_trend_contribution(series_count, observations, iterations))
         trend_scale = trend_variance
      end if
      if (use_seasonal_state) then
         allocate(seasonal_scale(series_count), &
            out%series_seasonal_state(seasonal_dimension, series_count, &
            observations, iterations), &
            out%series_seasonal_variance(series_count, iterations), &
            out%series_seasonal_contribution(series_count, observations, &
            iterations))
         seasonal_scale = seasonal_variance
      end if
      if (use_series_state .or. use_trend_state .or. use_seasonal_state) &
         allocate(combined_state_draws(dimension, observations))
      do slice = 1, transition_slices
         model%transition(:, :, slice) = identity_block(dimension, 1.0_dp)
         model%r(:, :, slice) = identity_block(dimension, 1.0_dp)
      end do
      if (use_trend_state) then
         do series = 1, series_count
            do slice = 1, transition_slices
               model%transition(trend_state_start + 2*series - 1, &
                  trend_state_start + 2*series, slice) = 1.0_dp
            end do
         end do
      end if
      if (use_seasonal_state) then
         do series = 1, series_count
            coefficient_start = seasonal_state_start + &
               (series - 1)*seasonal_dimension + 1
            coefficient_end = coefficient_start + seasonal_dimension - 1
            do slice = 1, observations
               if (modulo(slice, seasonal_duration_value) == 0) then
                  model%transition(coefficient_start:coefficient_end, &
                     coefficient_start:coefficient_end, slice) = 0.0_dp
                  model%transition(coefficient_start, &
                     coefficient_start:coefficient_end, slice) = -1.0_dp
                  do component = 2, seasonal_dimension
                     model%transition(coefficient_start + component - 1, &
                        coefficient_start + component - 2, slice) = 1.0_dp
                  end do
               end if
            end do
         end do
      end if
      model%q = 0.0_dp
      model%a1 = 0.0_dp
      model%a1(:factors) = initial_mean
      model%p1 = 0.0_dp
      model%p1(:factors, :factors) = initial_covariance
      if (use_regression) then
         do series = 1, series_count
            coefficient_start = factors + (series - 1)*predictors_count + 1
            coefficient_end = coefficient_start + predictors_count - 1
            model%a1(coefficient_start:coefficient_end) = &
               regression_prior_mean(:, series)
            model%p1(coefficient_start:coefficient_end, &
               coefficient_start:coefficient_end) = &
               regression_prior_covariance(:, :, series)
         end do
      end if
      if (use_series_state) then
         model%a1(series_state_start + 1: &
            series_state_start + series_count) = series_initial_mean
         do series = 1, series_count
            model%p1(series_state_start + series, &
               series_state_start + series) = series_initial_variance(series)
         end do
      end if
      if (use_trend_state) then
         do series = 1, series_count
            coefficient_start = trend_state_start + 2*series - 1
            coefficient_end = coefficient_start + 1
            model%a1(coefficient_start:coefficient_end) = &
               trend_initial_mean(:, series)
            model%p1(coefficient_start:coefficient_end, &
               coefficient_start:coefficient_end) = &
               trend_initial_covariance(:, :, series)
         end do
      end if
      if (use_seasonal_state) then
         do series = 1, series_count
            coefficient_start = seasonal_state_start + &
               (series - 1)*seasonal_dimension + 1
            coefficient_end = coefficient_start + seasonal_dimension - 1
            do component = coefficient_start, coefficient_end
               model%p1(component, component) = &
                  seasonal_initial_variance(series)
            end do
         end do
      end if
      model%p1inf = 0.0_dp
      model%missing = .not. ieee_is_finite(response)
      loadings = initial_loadings
      call impose_shared_loading_constraints(loadings)
      included = .true.
      do series = 1, series_count
         do factor = 1, factors
            if (factor > series) included(series, factor) = .false.
         end do
      end do
      if (use_selection) then
         included = initial_inclusion
         where (prior_inclusion_probability <= 0.0_dp) included = .false.
         where (prior_inclusion_probability >= 1.0_dp) included = .true.
         call impose_shared_inclusion_constraints(included)
         where (.not. included) loadings = 0.0_dp
         call impose_shared_loading_constraints(loadings)
      end if
      observation_scale = observation_variance
      factor_scale = factor_variance
      out%burn = burn
      out%spike_slab = use_selection
      out%is_mbsts = use_regression
      out%regression_spike_slab = use_regression_selection
      out%has_series_local_level = use_series_state
      out%has_series_local_linear_trend = use_trend_state
      out%has_series_seasonal = use_seasonal_state
      if (use_seasonal_state) then
         out%series_seasons = seasonal_nseasons
         out%series_season_duration = seasonal_duration_value
      end if
      do iteration = 1, iterations
         model%y = response
         if (present(offset_draws)) &
            model%y = model%y - offset_draws(:, :, iteration)
         model%z = 0.0_dp
         do time = 1, observations
            model%z(:, :factors, time) = loadings
            if (use_regression) then
               do series = 1, series_count
                  coefficient_start = factors + &
                     (series - 1)*predictors_count + 1
                  coefficient_end = coefficient_start + predictors_count - 1
                  model%z(series, coefficient_start:coefficient_end, time) = &
                     regression_predictors(time, :, series)
                  if (use_regression_selection) then
                     where (.not. regression_included(:, series))
                        model%z(series, coefficient_start:coefficient_end, &
                           time) = 0.0_dp
                     end where
                  end if
               end do
            end if
            if (use_series_state) then
               do series = 1, series_count
                  model%z(series, series_state_start + series, time) = 1.0_dp
               end do
            end if
            if (use_trend_state) then
               do series = 1, series_count
                  model%z(series, trend_state_start + 2*series - 1, time) = &
                     1.0_dp
               end do
            end if
            if (use_seasonal_state) then
               do series = 1, series_count
                  model%z(series, seasonal_state_start + &
                     (series - 1)*seasonal_dimension + 1, time) = 1.0_dp
               end do
            end if
         end do
         model%h = 0.0_dp
         do series = 1, series_count
            model%h(series, series, 1) = observation_scale(series)
         end do
         model%q = 0.0_dp
         do slice = 1, transition_slices
            do factor = 1, factors
               model%q(factor, factor, slice) = factor_scale(factor)
            end do
            if (use_series_state) then
               do series = 1, series_count
                  model%q(series_state_start + series, &
                     series_state_start + series, slice) = series_scale(series)
               end do
            end if
            if (use_trend_state) then
               do series = 1, series_count
                  do component = 1, 2
                     model%q(trend_state_start + 2*series - 2 + component, &
                        trend_state_start + 2*series - 2 + component, slice) = &
                        trend_scale(component, series)
                  end do
               end do
            end if
            if (use_seasonal_state .and. &
               modulo(slice, seasonal_duration_value) == 0) then
               do series = 1, series_count
                  component = seasonal_state_start + &
                     (series - 1)*seasonal_dimension + 1
                  model%q(component, component, slice) = seasonal_scale(series)
               end do
            end if
         end do
         if (use_series_state .or. use_trend_state .or. &
            use_seasonal_state) then
            combined_state_draws = 0.0_dp
            combined_state_draws(:base_dimension, :) = &
               state_normal_draws(:, :, iteration)
            if (use_series_state) &
            combined_state_draws(series_state_start + 1: &
               series_state_start + series_count, :) = &
               series_state_normal_draws(:, :, iteration)
            if (use_trend_state) then
               do series = 1, series_count
                  coefficient_start = trend_state_start + 2*series - 1
                  combined_state_draws(coefficient_start:coefficient_start + 1, &
                     :) = trend_state_normal_draws(:, series, :, iteration)
               end do
            end if
            if (use_seasonal_state) then
               do series = 1, series_count
                  coefficient_start = seasonal_state_start + &
                     (series - 1)*seasonal_dimension + 1
                  coefficient_end = coefficient_start + seasonal_dimension - 1
                  combined_state_draws(coefficient_start:coefficient_end, :) = &
                     seasonal_state_normal_draws(:, series, :, iteration)
               end do
            end if
            call state_path_from_draws(model, combined_state_draws, path, &
               status)
         else
            call state_path_from_draws(model, &
               state_normal_draws(:, :, iteration), path, status)
         end if
         if (status /= 0) then
            out%info = 10 + status
            return
         end if
         if (use_regression_selection) then
            do series = 1, series_count
               observed_index = pack([(time, time=1, observations)], &
                  ieee_is_finite(response(:, series)))
               observed = size(observed_index)
               allocate(predictors(observed, predictors_count), &
                  working_response(observed), &
                  regression_candidate(predictors_count), &
                  regression_draw(predictors_count))
               predictors = regression_predictors(observed_index, :, series)
               working_response = response(observed_index, series)
               if (present(offset_draws)) working_response = working_response - &
                  offset_draws(observed_index, series, iteration)
               do time = 1, observed
                  working_response(time) = working_response(time) - &
                     dot_product(loadings(series, :), &
                     path(:factors, observed_index(time)))
                  if (use_series_state) working_response(time) = &
                     working_response(time) - &
                     path(series_state_start + series, observed_index(time))
                  if (use_trend_state) working_response(time) = &
                     working_response(time) - &
                     path(trend_state_start + 2*series - 1, &
                     observed_index(time))
                  if (use_seasonal_state) working_response(time) = &
                     working_response(time) - path(seasonal_state_start + &
                     (series - 1)*seasonal_dimension + 1, observed_index(time))
               end do
               flips = 0
               do factor = 1, predictors_count
                  if (regression_prior_inclusion_probability( &
                     factor, series) <= 0.0_dp) then
                     regression_included(factor, series) = .false.
                     cycle
                  end if
                  if (regression_prior_inclusion_probability( &
                     factor, series) >= 1.0_dp) then
                     regression_included(factor, series) = .true.
                     cycle
                  end if
                  regression_candidate = regression_included(:, series)
                  regression_candidate(factor) = .false.
                  call spike_slab_model_moments(working_response, predictors, &
                     regression_candidate, regression_prior_mean(:, series), &
                     regression_prior_covariance(:, :, series), &
                     regression_prior_inclusion_probability(:, series), &
                     observation_scale(series), observation_prior_rate(series), &
                     log_excluded, collapsed_rate, posterior_mean, covariance, &
                     regression_active, status)
                  if (status /= 0) then
                     out%info = 80 + status
                     return
                  end if
                  regression_candidate(factor) = .true.
                  if (count(regression_candidate) > regression_model_limit) then
                     probability = 0.0_dp
                  else
                     call spike_slab_model_moments(working_response, predictors, &
                        regression_candidate, &
                        regression_prior_mean(:, series), &
                        regression_prior_covariance(:, :, series), &
                        regression_prior_inclusion_probability(:, series), &
                        observation_scale(series), &
                        observation_prior_rate(series), log_included, &
                        collapsed_rate, posterior_mean, covariance, &
                        regression_active, status)
                     if (status /= 0) then
                        out%info = 90 + status
                        return
                     end if
                     probability = logistic_log_odds( &
                        log_included - log_excluded)
                  end if
                  regression_candidate(factor) = regression_uniform_draws( &
                     factor, series, iteration) < probability
                  if (regression_candidate(factor) .neqv. &
                     regression_included(factor, series)) then
                     if (regression_flip_limit < 0 .or. &
                        flips < regression_flip_limit) then
                        regression_included(factor, series) = &
                           regression_candidate(factor)
                        flips = flips + 1
                     end if
                  end if
               end do
               call spike_slab_model_moments(working_response, predictors, &
                  regression_included(:, series), &
                  regression_prior_mean(:, series), &
                  regression_prior_covariance(:, :, series), &
                  regression_prior_inclusion_probability(:, series), &
                  observation_scale(series), observation_prior_rate(series), &
                  log_probability, collapsed_rate, posterior_mean, covariance, &
                  regression_active, status)
               if (status /= 0) then
                  out%info = 100 + status
                  return
               end if
               regression_draw = 0.0_dp
               if (size(regression_active) > 0) then
                  allocate(active_draw(size(regression_active)))
                  call multivariate_normal_from_psd(posterior_mean, &
                     observation_scale(series)*covariance, &
                     regression_normal_draws(regression_active, series, &
                     iteration), active_draw, status)
                  if (status /= 0) then
                     out%info = 110 + status
                     return
                  end if
                  regression_draw(regression_active) = active_draw
                  deallocate(active_draw)
               end if
               coefficient_start = factors + &
                  (series - 1)*predictors_count + 1
               coefficient_end = coefficient_start + predictors_count - 1
               path(coefficient_start:coefficient_end, :) = &
                  spread(regression_draw, 2, observations)
               deallocate(predictors, working_response, regression_candidate, &
                  regression_draw, posterior_mean)
            end do
         end if
         do series = 1, series_count
            free_count = min(factors, series - 1)
            loadings(series, :) = 0.0_dp
            if (series <= factors) loadings(series, series) = 1.0_dp
            if (free_count < 1) cycle
            observed_index = pack([(time, time=1, observations)], &
               ieee_is_finite(response(:, series)))
            observed = size(observed_index)
            allocate(predictors(observed, free_count), &
               working_response(observed), coefficient(free_count))
            predictors = transpose(path(:free_count, observed_index))
            working_response = response(observed_index, series)
            if (present(offset_draws)) working_response = working_response - &
               offset_draws(observed_index, series, iteration)
            if (use_regression) then
               coefficient_start = factors + &
                  (series - 1)*predictors_count + 1
               coefficient_end = coefficient_start + predictors_count - 1
               do time = 1, observed
                  working_response(time) = working_response(time) - &
                     dot_product(regression_predictors(observed_index(time), &
                     :, series), path(coefficient_start:coefficient_end, &
                     observed_index(time)))
               end do
            end if
            if (use_series_state) working_response = working_response - &
               path(series_state_start + series, observed_index)
            if (use_trend_state) working_response = working_response - &
               path(trend_state_start + 2*series - 1, observed_index)
            if (use_seasonal_state) working_response = working_response - &
               path(seasonal_state_start + &
               (series - 1)*seasonal_dimension + 1, observed_index)
            if (series <= factors) working_response = working_response - &
               path(series, observed_index)
            if (use_selection) then
               allocate(slab_covariance(free_count, free_count), &
                  candidate(free_count))
               slab_covariance = 0.0_dp
               do factor = 1, free_count
                  slab_covariance(factor, factor) = &
                     loading_prior_variance(series, factor)
               end do
               flips = 0
               do factor = 1, free_count
                  if (prior_inclusion_probability(series, factor) <= 0.0_dp) then
                     included(series, factor) = .false.
                     cycle
                  end if
                  if (prior_inclusion_probability(series, factor) >= 1.0_dp) then
                     included(series, factor) = .true.
                     cycle
                  end if
                  candidate = included(series, :free_count)
                  candidate(factor) = .false.
                  call spike_slab_model_moments(working_response, predictors, &
                     candidate, loading_prior_mean(series, :free_count), &
                     slab_covariance, &
                     prior_inclusion_probability(series, :free_count), &
                     observation_scale(series), &
                     observation_prior_rate(series), log_excluded, &
                     collapsed_rate, posterior_mean, covariance, active, status)
                  if (status /= 0) then
                     out%info = 20 + status
                     return
                  end if
                  candidate(factor) = .true.
                  call spike_slab_model_moments(working_response, predictors, &
                     candidate, loading_prior_mean(series, :free_count), &
                     slab_covariance, &
                     prior_inclusion_probability(series, :free_count), &
                     observation_scale(series), &
                     observation_prior_rate(series), log_included, &
                     collapsed_rate, posterior_mean, covariance, active, status)
                  if (status /= 0) then
                     out%info = 30 + status
                     return
                  end if
                  probability = logistic_log_odds(log_included - log_excluded)
                  candidate(factor) = inclusion_uniform_draws( &
                     factor, series, iteration) < probability
                  if (candidate(factor) .neqv. included(series, factor)) then
                     if (flip_limit < 0 .or. flips < flip_limit) then
                        included(series, factor) = candidate(factor)
                        flips = flips + 1
                     end if
                  end if
               end do
               call spike_slab_model_moments(working_response, predictors, &
                  included(series, :free_count), &
                  loading_prior_mean(series, :free_count), slab_covariance, &
                  prior_inclusion_probability(series, :free_count), &
                  observation_scale(series), observation_prior_rate(series), &
                  log_probability, collapsed_rate, posterior_mean, covariance, &
                  active, status)
               if (status /= 0) then
                  out%info = 40 + status
                  return
               end if
               coefficient = 0.0_dp
               if (size(active) > 0) then
                  allocate(active_draw(size(active)))
                  call multivariate_normal_from_psd(posterior_mean, &
                     observation_scale(series)*covariance, &
                     loading_normal_draws(active, series, iteration), &
                     active_draw, status)
                  if (status /= 0) then
                     out%info = 50 + status
                     return
                  end if
                  coefficient(active) = active_draw
                  deallocate(active_draw)
               end if
               deallocate(slab_covariance, candidate)
            else
               allocate(precision(free_count, free_count), score(free_count), &
                  posterior_mean(free_count))
               precision = matmul(transpose(predictors), predictors)/ &
                  observation_scale(series)
               score = matmul(transpose(predictors), working_response)/ &
                  observation_scale(series)
               do factor = 1, free_count
                  precision(factor, factor) = precision(factor, factor) + &
                     1.0_dp/loading_prior_variance(series, factor)
                  score(factor) = score(factor) + &
                     loading_prior_mean(series, factor)/ &
                     loading_prior_variance(series, factor)
               end do
               call invert_matrix(precision, covariance, status)
               if (status /= 0) then
                  out%info = 60 + status
                  return
               end if
               posterior_mean = matmul(covariance, score)
               call multivariate_normal_from_psd(posterior_mean, covariance, &
                  loading_normal_draws(:free_count, series, iteration), &
                  coefficient, status)
               if (status /= 0) then
                  out%info = 70 + status
                  return
               end if
               deallocate(precision, score)
            end if
            loadings(series, :free_count) = coefficient
            deallocate(predictors, working_response, posterior_mean, coefficient)
         end do
         do series = 1, series_count
            residual_sum = 0.0_dp
            observed = 0
            do time = 1, observations
               if (.not. ieee_is_finite(response(time, series))) cycle
               fixed_part = dot_product(loadings(series, :), &
                  path(:factors, time))
               if (use_regression) then
                  coefficient_start = factors + &
                     (series - 1)*predictors_count + 1
                  coefficient_end = coefficient_start + predictors_count - 1
                  fixed_part = fixed_part + dot_product( &
                     regression_predictors(time, :, series), &
                     path(coefficient_start:coefficient_end, time))
               end if
               if (use_series_state) fixed_part = fixed_part + &
                  path(series_state_start + series, time)
               if (use_trend_state) fixed_part = fixed_part + &
                  path(trend_state_start + 2*series - 1, time)
               if (use_seasonal_state) fixed_part = fixed_part + &
                  path(seasonal_state_start + &
                  (series - 1)*seasonal_dimension + 1, time)
               if (present(offset_draws)) fixed_part = fixed_part + &
                  offset_draws(time, series, iteration)
               residual_sum = residual_sum + &
                  (response(time, series) - fixed_part)**2
               observed = observed + 1
            end do
            observation_scale(series) = &
               (observation_prior_rate(series) + 0.5_dp*residual_sum)/ &
               gamma_draws(series, iteration)
         end do
         do factor = 1, factors
            difference_sum = sum((path(factor, 2:) - path(factor, &
               :observations - 1))**2)
            factor_scale(factor) = &
               (factor_prior_rate(factor) + 0.5_dp*difference_sum)/ &
               gamma_draws(series_count + factor, iteration)
         end do
         if (use_series_state) then
            do series = 1, series_count
               difference_sum = sum((path(series_state_start + series, 2:) - &
                  path(series_state_start + series, &
                  :observations - 1))**2)
               series_scale(series) = &
                  (series_prior_rate(series) + 0.5_dp*difference_sum)/ &
                  series_gamma_draws(series, iteration)
            end do
         end if
         if (use_trend_state) then
            do series = 1, series_count
               component = trend_state_start + 2*series - 1
               difference_sum = sum((path(component, 2:) - &
                  path(component, :observations - 1) - &
                  path(component + 1, :observations - 1))**2)
               trend_scale(1, series) = (trend_prior_rate(1, series) + &
                  0.5_dp*difference_sum)/trend_gamma_draws(1, series, iteration)
               difference_sum = sum((path(component + 1, 2:) - &
                  path(component + 1, :observations - 1))**2)
               trend_scale(2, series) = (trend_prior_rate(2, series) + &
                  0.5_dp*difference_sum)/trend_gamma_draws(2, series, iteration)
            end do
         end if
         if (use_seasonal_state) then
            do series = 1, series_count
               component = seasonal_state_start + &
                  (series - 1)*seasonal_dimension + 1
               difference_sum = 0.0_dp
               do time = 2, observations
                  if (modulo(time - 1, seasonal_duration_value) == 0) &
                     difference_sum = difference_sum + &
                     (path(component, time) + sum(path(component: &
                     component + seasonal_dimension - 1, time - 1)))**2
               end do
               seasonal_scale(series) = (seasonal_prior_rate(series) + &
                  0.5_dp*difference_sum)/seasonal_gamma_draws(series, iteration)
            end do
         end if
         out%state(:, :, iteration) = path(:factors, :)
         out%loadings(:, :, iteration) = loadings
         out%observation_variance(:, iteration) = observation_scale
         out%factor_variance(:, iteration) = factor_scale
         out%contribution(:, :, iteration) = matmul(loadings, path(:factors, :))
         if (use_series_state) then
            out%series_state(:, :, iteration) = &
               path(series_state_start + 1: &
               series_state_start + series_count, :)
            out%series_variance(:, iteration) = series_scale
            out%series_contribution(:, :, iteration) = &
               path(series_state_start + 1: &
               series_state_start + series_count, :)
         end if
         if (use_trend_state) then
            do series = 1, series_count
               component = trend_state_start + 2*series - 1
               out%series_trend_state(:, series, :, iteration) = &
                  path(component:component + 1, :)
               out%series_trend_contribution(series, :, iteration) = &
                  path(component, :)
            end do
            out%series_trend_variance(:, :, iteration) = trend_scale
         end if
         if (use_seasonal_state) then
            do series = 1, series_count
               component = seasonal_state_start + &
                  (series - 1)*seasonal_dimension + 1
               out%series_seasonal_state(:, series, :, iteration) = &
                  path(component:component + seasonal_dimension - 1, :)
               out%series_seasonal_contribution(series, :, iteration) = &
                  path(component, :)
            end do
            out%series_seasonal_variance(:, iteration) = seasonal_scale
         end if
         if (use_regression) then
            do series = 1, series_count
               coefficient_start = factors + &
                  (series - 1)*predictors_count + 1
               coefficient_end = coefficient_start + predictors_count - 1
               out%regression_coefficients(:, series, iteration) = &
                  path(coefficient_start:coefficient_end, observations)
               do time = 1, observations
                  out%regression_contribution(series, time, iteration) = &
                     dot_product(regression_predictors(time, :, series), &
                     out%regression_coefficients(:, series, iteration))
               end do
            end do
            out%fitted(:, :, iteration) = out%contribution(:, :, iteration) + &
               out%regression_contribution(:, :, iteration)
            if (use_series_state) out%fitted(:, :, iteration) = &
               out%fitted(:, :, iteration) + &
               out%series_contribution(:, :, iteration)
            if (use_trend_state) out%fitted(:, :, iteration) = &
               out%fitted(:, :, iteration) + &
               out%series_trend_contribution(:, :, iteration)
            if (use_seasonal_state) out%fitted(:, :, iteration) = &
               out%fitted(:, :, iteration) + &
               out%series_seasonal_contribution(:, :, iteration)
            if (present(offset_draws)) out%fitted(:, :, iteration) = &
               out%fitted(:, :, iteration) + &
               transpose(offset_draws(:, :, iteration))
            out%residuals(:, :, iteration) = &
               transpose(response) - out%fitted(:, :, iteration)
            if (use_regression_selection) &
               out%regression_included(:, :, iteration) = regression_included
         end if
         out%included(:, :, iteration) = included
      end do
      retained = iterations - burn
      do series = 1, series_count
         do factor = 1, factors
            out%inclusion_probability(series, factor) = real(count( &
               out%included(series, factor, burn + 1:)), dp)/real(retained, dp)
         end do
      end do
      if (use_regression_selection) then
         do series = 1, series_count
            do factor = 1, predictors_count
               out%regression_inclusion_probability(factor, series) = &
                  real(count(out%regression_included(factor, series, &
                  burn + 1:)), dp)/real(retained, dp)
            end do
         end do
      end if
   end function bsts_shared_local_level_draws

   function bsts_shared_local_level(response, nfactors, iterations, burn, &
      offset_draws, prior_inclusion_probability, initial_inclusion, &
      maximum_flips, select_loadings) result(out)
      !! Sample a shared multivariate local-level model using the random stream.
      real(dp), intent(in) :: response(:, :) !! Response observations.
      integer, intent(in) :: nfactors !! Nfactors.
      integer, intent(in) :: iterations !! Number of algorithm iterations.
      integer, intent(in), optional :: burn !! Number of initial simulation draws to discard.
      real(dp), intent(in), optional :: offset_draws(:, :, :) !! Offset simulation draws.
      real(dp), intent(in), optional :: prior_inclusion_probability(:, :) !! Prior inclusion probability.
      logical, intent(in), optional :: initial_inclusion(:, :) !! Initial inclusion.
      integer, intent(in), optional :: maximum_flips !! Maximum flips.
      logical, intent(in), optional :: select_loadings !! Flag controlling select loadings.
      type(bsts_shared_local_level_t) :: out
      real(dp), allocatable :: initial_mean(:), initial_covariance(:, :)
      real(dp), allocatable :: initial_loadings(:, :), loading_mean(:, :)
      real(dp), allocatable :: loading_variance(:, :), observation_scale(:)
      real(dp), allocatable :: factor_scale(:), observation_shape(:)
      real(dp), allocatable :: observation_rate(:), factor_shape(:)
      real(dp), allocatable :: factor_rate(:), state_normals(:, :, :)
      real(dp), allocatable :: loading_normals(:, :, :), gammas(:, :)
      real(dp), allocatable :: inclusion0(:, :), uniforms(:, :, :)
      logical, allocatable :: model0(:, :)
      real(dp) :: largest_variance
      integer :: observations, series_count, discarded, iteration, flip_limit
      integer :: series, factor, observed
      logical :: use_selection

      observations = size(response, 1)
      series_count = size(response, 2)
      if (observations < 2 .or. series_count < 1 .or. nfactors < 1 .or. &
         nfactors > series_count .or. iterations < 1 .or. &
         any([(count(ieee_is_finite(response(:, series))) < 1, &
         series=1, series_count)])) then
         out%info = 1
         return
      end if
      discarded = 0
      if (present(burn)) discarded = burn
      if (discarded < 0 .or. discarded >= iterations) then
         out%info = 1
         return
      end if
      use_selection = .true.
      if (present(select_loadings)) use_selection = select_loadings
      flip_limit = -1
      if (present(maximum_flips)) flip_limit = maximum_flips
      if (flip_limit < -1) then
         out%info = 1
         return
      end if
      allocate(initial_mean(nfactors), &
         initial_covariance(nfactors, nfactors), &
         initial_loadings(series_count, nfactors), &
         loading_mean(series_count, nfactors), &
         loading_variance(series_count, nfactors), &
         observation_scale(series_count), factor_scale(nfactors), &
         observation_shape(series_count), observation_rate(series_count), &
         factor_shape(nfactors), factor_rate(nfactors), &
         state_normals(nfactors, observations, iterations), &
         loading_normals(nfactors, series_count, iterations), &
         gammas(series_count + nfactors, iterations), &
         inclusion0(series_count, nfactors), model0(series_count, nfactors), &
         uniforms(nfactors, series_count, iterations))
      largest_variance = 1.0e-6_dp
      do series = 1, series_count
         observation_scale(series) = max(0.1_dp* &
            finite_sample_variance(response(:, series)), 1.0e-8_dp)
         largest_variance = max(largest_variance, &
            finite_sample_variance(response(:, series)))
      end do
      initial_mean = 0.0_dp
      initial_covariance = identity_block(nfactors, largest_variance)
      initial_loadings = 1.0_dp
      call impose_shared_loading_constraints(initial_loadings)
      loading_mean = 1.0_dp
      loading_variance = 1.0_dp
      inclusion0 = 0.5_dp
      if (present(prior_inclusion_probability)) then
         if (any(shape(prior_inclusion_probability) /= &
            [series_count, nfactors])) then
            out%info = 1
            return
         end if
         inclusion0 = prior_inclusion_probability
      end if
      model0 = inclusion0 > 0.5_dp
      if (present(initial_inclusion)) then
         if (any(shape(initial_inclusion) /= [series_count, nfactors])) then
            out%info = 1
            return
         end if
         model0 = initial_inclusion
      end if
      call impose_shared_inclusion_priors(inclusion0)
      where (inclusion0 <= 0.0_dp) model0 = .false.
      where (inclusion0 >= 1.0_dp) model0 = .true.
      call impose_shared_inclusion_constraints(model0)
      factor_scale = max(0.01_dp*largest_variance, 1.0e-10_dp)
      observation_shape = 2.0_dp
      observation_rate = observation_scale
      factor_shape = 2.0_dp
      factor_rate = factor_scale
      do iteration = 1, iterations
         do factor = 1, nfactors
            do series = 1, observations
               state_normals(factor, series, iteration) = &
                  random_standard_normal()
            end do
         end do
         do series = 1, series_count
            do factor = 1, nfactors
               loading_normals(factor, series, iteration) = &
                  random_standard_normal()
               uniforms(factor, series, iteration) = max(tiny(1.0_dp), &
                  min(1.0_dp - epsilon(1.0_dp), random_uniform()))
            end do
            observed = count(ieee_is_finite(response(:, series)))
            gammas(series, iteration) = random_gamma( &
               observation_shape(series) + 0.5_dp*real(observed, dp))
         end do
         do factor = 1, nfactors
            gammas(series_count + factor, iteration) = random_gamma( &
               factor_shape(factor) + 0.5_dp*real(observations - 1, dp))
         end do
      end do
      if (use_selection .and. present(offset_draws)) then
         out = bsts_shared_local_level_draws(response, initial_mean, &
            initial_covariance, initial_loadings, loading_mean, &
            loading_variance, observation_scale, factor_scale, &
            observation_shape, observation_rate, factor_shape, factor_rate, &
            discarded, state_normals, loading_normals, gammas, offset_draws, &
            inclusion0, model0, uniforms, flip_limit)
      else if (use_selection) then
         out = bsts_shared_local_level_draws(response, initial_mean, &
            initial_covariance, initial_loadings, loading_mean, &
            loading_variance, observation_scale, factor_scale, &
            observation_shape, observation_rate, factor_shape, factor_rate, &
            discarded, state_normals, loading_normals, gammas, &
            prior_inclusion_probability=inclusion0, initial_inclusion=model0, &
            inclusion_uniform_draws=uniforms, maximum_flips=flip_limit)
      else if (present(offset_draws)) then
         out = bsts_shared_local_level_draws(response, initial_mean, &
            initial_covariance, initial_loadings, loading_mean, &
            loading_variance, observation_scale, factor_scale, &
            observation_shape, observation_rate, factor_shape, factor_rate, &
            discarded, state_normals, loading_normals, gammas, offset_draws)
      else
         out = bsts_shared_local_level_draws(response, initial_mean, &
            initial_covariance, initial_loadings, loading_mean, &
            loading_variance, observation_scale, factor_scale, &
            observation_shape, observation_rate, factor_shape, factor_rate, &
            discarded, state_normals, loading_normals, gammas)
      end if
   end function bsts_shared_local_level

   pure function bsts_shared_local_level_predict_draws(fit, horizon, &
      state_normal_draws, observation_normal_draws) result(out)
      !! Forecast shared local levels from supplied independent normal draws.
      type(bsts_shared_local_level_t), intent(in) :: fit !! Previously fitted model.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      real(dp), intent(in) :: state_normal_draws(:, :, :) !! Independent standard-normal state draws.
      real(dp), intent(in) :: observation_normal_draws(:, :, :) !! Observation normal draws.
      type(bsts_multivariate_prediction_t) :: out
      real(dp), allocatable :: state(:), sorted(:)
      integer :: factors, series_count, retained, draw, source, step, series

      if (fit%info /= 0 .or. .not. allocated(fit%state) .or. &
         .not. allocated(fit%loadings) .or. &
         .not. allocated(fit%observation_variance) .or. &
         .not. allocated(fit%factor_variance)) then
         out%info = 1
         return
      end if
      factors = size(fit%state, 1)
      series_count = size(fit%loadings, 1)
      retained = size(fit%state, 3) - fit%burn
      if (horizon < 1 .or. retained < 1 .or. &
         any(shape(state_normal_draws) /= [factors, horizon, retained]) .or. &
         any(shape(observation_normal_draws) /= &
         [series_count, horizon, retained]) .or. &
         .not. all(ieee_is_finite(state_normal_draws)) .or. &
         .not. all(ieee_is_finite(observation_normal_draws))) then
         out%info = 1
         return
      end if
      allocate(out%draws(series_count, horizon, retained), &
         out%mean(series_count, horizon), &
         out%standard_deviation(series_count, horizon), &
         out%lower(series_count, horizon), out%upper(series_count, horizon), &
         state(factors))
      do draw = 1, retained
         source = fit%burn + draw
         state = fit%state(:, size(fit%state, 2), source)
         do step = 1, horizon
            state = state + sqrt(fit%factor_variance(:, source))* &
               state_normal_draws(:, step, draw)
            out%draws(:, step, draw) = &
               matmul(fit%loadings(:, :, source), state) + &
               sqrt(fit%observation_variance(:, source))* &
               observation_normal_draws(:, step, draw)
         end do
      end do
      do series = 1, series_count
         do step = 1, horizon
            out%mean(series, step) = sum(out%draws(series, step, :))/ &
               real(retained, dp)
            if (retained > 1) then
               out%standard_deviation(series, step) = &
                  sqrt(sum((out%draws(series, step, :) - &
                  out%mean(series, step))**2)/real(retained - 1, dp))
            else
               out%standard_deviation(series, step) = 0.0_dp
            end if
            sorted = out%draws(series, step, :)
            call insertion_sort(sorted)
            out%lower(series, step) = quantile(sorted, 0.025_dp)
            out%upper(series, step) = quantile(sorted, 0.975_dp)
         end do
      end do
   end function bsts_shared_local_level_predict_draws

   function bsts_shared_local_level_predict(fit, horizon) result(out)
      !! Forecast shared local levels using the shared random stream.
      type(bsts_shared_local_level_t), intent(in) :: fit !! Previously fitted model.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      type(bsts_multivariate_prediction_t) :: out
      real(dp), allocatable :: state_normals(:, :, :)
      real(dp), allocatable :: observation_normals(:, :, :)
      integer :: factors, series_count, retained, draw, step, component

      if (.not. allocated(fit%state) .or. .not. allocated(fit%loadings)) then
         out%info = 1
         return
      end if
      factors = size(fit%state, 1)
      series_count = size(fit%loadings, 1)
      retained = size(fit%state, 3) - fit%burn
      if (horizon < 1 .or. retained < 1) then
         out%info = 1
         return
      end if
      allocate(state_normals(factors, horizon, retained), &
         observation_normals(series_count, horizon, retained))
      do draw = 1, retained
         do step = 1, horizon
            do component = 1, factors
               state_normals(component, step, draw) = random_standard_normal()
            end do
            do component = 1, series_count
               observation_normals(component, step, draw) = &
                  random_standard_normal()
            end do
         end do
      end do
      out = bsts_shared_local_level_predict_draws(fit, horizon, state_normals, &
         observation_normals)
   end function bsts_shared_local_level_predict

   pure function bsts_mbsts_draws(response, regression_predictors, &
      regression_prior_mean, regression_prior_covariance, initial_mean, &
      initial_covariance, initial_loadings, loading_prior_mean, &
      loading_prior_variance, observation_variance, factor_variance, &
      observation_prior_shape, observation_prior_rate, factor_prior_shape, &
      factor_prior_rate, burn, state_normal_draws, loading_normal_draws, &
      gamma_draws, offset_draws, regression_prior_inclusion_probability, &
      initial_regression_inclusion, regression_normal_draws, &
      regression_uniform_draws, maximum_model_size, maximum_flips, &
      series_initial_mean, series_initial_variance, series_variance, &
      series_prior_shape, series_prior_rate, series_state_normal_draws, &
      series_gamma_draws, trend_initial_mean, trend_initial_covariance, &
      trend_variance, trend_prior_shape, trend_prior_rate, &
      trend_state_normal_draws, trend_gamma_draws, seasonal_nseasons, &
      seasonal_duration, seasonal_initial_variance, seasonal_variance, &
      seasonal_prior_shape, seasonal_prior_rate, &
      seasonal_state_normal_draws, seasonal_gamma_draws, &
      enable_series_state, enable_trend_state, &
      enable_seasonal_state) result(out)
      !! Sample joint factors, regressions, and series-specific states.
      real(dp), intent(in) :: response(:, :) !! Response observations.
      real(dp), intent(in) :: regression_predictors(:, :, :) !! Regression predictors.
      real(dp), intent(in) :: regression_prior_mean(:, :) !! Regression prior mean.
      real(dp), intent(in) :: regression_prior_covariance(:, :, :) !! Regression prior covariance matrix.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: initial_loadings(:, :) !! Initial loadings.
      real(dp), intent(in) :: loading_prior_mean(:, :) !! Loading prior mean.
      real(dp), intent(in) :: loading_prior_variance(:, :) !! Loading prior variance.
      real(dp), intent(in) :: observation_variance(:) !! Observation-error variance.
      real(dp), intent(in) :: factor_variance(:) !! Factor variance.
      real(dp), intent(in) :: observation_prior_shape(:) !! Observation prior shape.
      real(dp), intent(in) :: observation_prior_rate(:) !! Observation prior rate.
      real(dp), intent(in) :: factor_prior_shape(:) !! Factor prior shape.
      real(dp), intent(in) :: factor_prior_rate(:) !! Factor prior rate.
      integer, intent(in) :: burn !! Number of initial simulation draws to discard.
      real(dp), intent(in) :: state_normal_draws(:, :, :) !! Independent standard-normal state draws.
      real(dp), intent(in) :: loading_normal_draws(:, :, :) !! Loading normal simulation draws.
      real(dp), intent(in) :: gamma_draws(:, :) !! Gamma simulation draws.
      real(dp), intent(in), optional :: offset_draws(:, :, :) !! Offset simulation draws.
      real(dp), intent(in), optional :: regression_prior_inclusion_probability(:, :) !! Regression prior inclusion probability.
      logical, intent(in), optional :: initial_regression_inclusion(:, :) !! Initial regression inclusion.
      real(dp), intent(in), optional :: regression_normal_draws(:, :, :) !! Regression normal simulation draws.
      real(dp), intent(in), optional :: regression_uniform_draws(:, :, :) !! Regression uniform simulation draws.
      integer, intent(in), optional :: maximum_model_size !! Maximum model size.
      integer, intent(in), optional :: maximum_flips !! Maximum flips.
      real(dp), intent(in), optional :: series_initial_mean(:) !! Series initial mean.
      real(dp), intent(in), optional :: series_initial_variance(:) !! Series initial variance.
      real(dp), intent(in), optional :: series_variance(:) !! Series variance.
      real(dp), intent(in), optional :: series_prior_shape(:) !! Series prior shape.
      real(dp), intent(in), optional :: series_prior_rate(:) !! Series prior rate.
      real(dp), intent(in), optional :: series_state_normal_draws(:, :, :) !! Series state normal simulation draws.
      real(dp), intent(in), optional :: series_gamma_draws(:, :) !! Series gamma simulation draws.
      real(dp), intent(in), optional :: trend_initial_mean(:, :) !! Trend initial mean.
      real(dp), intent(in), optional :: trend_initial_covariance(:, :, :) !! Trend initial covariance matrix.
      real(dp), intent(in), optional :: trend_variance(:, :) !! Trend variance.
      real(dp), intent(in), optional :: trend_prior_shape(:, :) !! Trend prior shape.
      real(dp), intent(in), optional :: trend_prior_rate(:, :) !! Trend prior rate.
      real(dp), intent(in), optional :: trend_state_normal_draws(:, :, :, :) !! Trend state normal simulation draws.
      real(dp), intent(in), optional :: trend_gamma_draws(:, :, :) !! Trend gamma simulation draws.
      integer, intent(in), optional :: seasonal_nseasons !! Seasonal nseasons.
      integer, intent(in), optional :: seasonal_duration !! Seasonal duration.
      real(dp), intent(in), optional :: seasonal_initial_variance(:) !! Seasonal initial variance.
      real(dp), intent(in), optional :: seasonal_variance(:) !! Seasonal variance.
      real(dp), intent(in), optional :: seasonal_prior_shape(:) !! Seasonal prior shape.
      real(dp), intent(in), optional :: seasonal_prior_rate(:) !! Seasonal prior rate.
      real(dp), intent(in), optional :: seasonal_state_normal_draws(:, :, :, :) !! Seasonal state normal simulation draws.
      real(dp), intent(in), optional :: seasonal_gamma_draws(:, :) !! Seasonal gamma simulation draws.
      logical, intent(in), optional :: enable_series_state !! Flag controlling enable series state.
      logical, intent(in), optional :: enable_trend_state !! Flag controlling enable trend state.
      logical, intent(in), optional :: enable_seasonal_state !! Flag controlling enable seasonal state.
      type(bsts_shared_local_level_t) :: out

      if (present(offset_draws)) then
         out = bsts_shared_local_level_draws(response, initial_mean, &
            initial_covariance, initial_loadings, loading_prior_mean, &
            loading_prior_variance, observation_variance, factor_variance, &
            observation_prior_shape, observation_prior_rate, &
            factor_prior_shape, factor_prior_rate, burn, state_normal_draws, &
            loading_normal_draws, gamma_draws, offset_draws, &
            regression_predictors=regression_predictors, &
            regression_prior_mean=regression_prior_mean, &
            regression_prior_covariance=regression_prior_covariance, &
            regression_prior_inclusion_probability= &
            regression_prior_inclusion_probability, &
            initial_regression_inclusion=initial_regression_inclusion, &
            regression_normal_draws=regression_normal_draws, &
            regression_uniform_draws=regression_uniform_draws, &
            regression_maximum_model_size=maximum_model_size, &
            regression_maximum_flips=maximum_flips, &
            series_initial_mean=series_initial_mean, &
            series_initial_variance=series_initial_variance, &
            series_variance=series_variance, &
            series_prior_shape=series_prior_shape, &
            series_prior_rate=series_prior_rate, &
            series_state_normal_draws=series_state_normal_draws, &
            series_gamma_draws=series_gamma_draws, &
            trend_initial_mean=trend_initial_mean, &
            trend_initial_covariance=trend_initial_covariance, &
            trend_variance=trend_variance, &
            trend_prior_shape=trend_prior_shape, &
            trend_prior_rate=trend_prior_rate, &
            trend_state_normal_draws=trend_state_normal_draws, &
            trend_gamma_draws=trend_gamma_draws, &
            seasonal_nseasons=seasonal_nseasons, &
            seasonal_duration=seasonal_duration, &
            seasonal_initial_variance=seasonal_initial_variance, &
            seasonal_variance=seasonal_variance, &
            seasonal_prior_shape=seasonal_prior_shape, &
            seasonal_prior_rate=seasonal_prior_rate, &
            seasonal_state_normal_draws=seasonal_state_normal_draws, &
            seasonal_gamma_draws=seasonal_gamma_draws, &
            enable_series_state=enable_series_state, &
            enable_trend_state=enable_trend_state, &
            enable_seasonal_state=enable_seasonal_state)
      else
         out = bsts_shared_local_level_draws(response, initial_mean, &
            initial_covariance, initial_loadings, loading_prior_mean, &
            loading_prior_variance, observation_variance, factor_variance, &
            observation_prior_shape, observation_prior_rate, &
            factor_prior_shape, factor_prior_rate, burn, state_normal_draws, &
            loading_normal_draws, gamma_draws, &
            regression_predictors=regression_predictors, &
            regression_prior_mean=regression_prior_mean, &
            regression_prior_covariance=regression_prior_covariance, &
            regression_prior_inclusion_probability= &
            regression_prior_inclusion_probability, &
            initial_regression_inclusion=initial_regression_inclusion, &
            regression_normal_draws=regression_normal_draws, &
            regression_uniform_draws=regression_uniform_draws, &
            regression_maximum_model_size=maximum_model_size, &
            regression_maximum_flips=maximum_flips, &
            series_initial_mean=series_initial_mean, &
            series_initial_variance=series_initial_variance, &
            series_variance=series_variance, &
            series_prior_shape=series_prior_shape, &
            series_prior_rate=series_prior_rate, &
            series_state_normal_draws=series_state_normal_draws, &
            series_gamma_draws=series_gamma_draws, &
            trend_initial_mean=trend_initial_mean, &
            trend_initial_covariance=trend_initial_covariance, &
            trend_variance=trend_variance, &
            trend_prior_shape=trend_prior_shape, &
            trend_prior_rate=trend_prior_rate, &
            trend_state_normal_draws=trend_state_normal_draws, &
            trend_gamma_draws=trend_gamma_draws, &
            seasonal_nseasons=seasonal_nseasons, &
            seasonal_duration=seasonal_duration, &
            seasonal_initial_variance=seasonal_initial_variance, &
            seasonal_variance=seasonal_variance, &
            seasonal_prior_shape=seasonal_prior_shape, &
            seasonal_prior_rate=seasonal_prior_rate, &
            seasonal_state_normal_draws=seasonal_state_normal_draws, &
            seasonal_gamma_draws=seasonal_gamma_draws, &
            enable_series_state=enable_series_state, &
            enable_trend_state=enable_trend_state, &
            enable_seasonal_state=enable_seasonal_state)
      end if
   end function bsts_mbsts_draws

   function bsts_mbsts(response, regression_predictors, nfactors, iterations, &
      burn, offset_draws, regression_prior_inclusion_probability, &
      initial_regression_inclusion, maximum_model_size, maximum_flips, &
      select_regression, series_local_level, series_local_linear_trend, &
      series_nseasons, series_season_duration) result(out)
      !! Sample multivariate BSTS with optional series-specific states.
      real(dp), intent(in) :: response(:, :) !! Response observations.
      real(dp), intent(in) :: regression_predictors(:, :, :) !! Regression predictors.
      integer, intent(in) :: nfactors !! Nfactors.
      integer, intent(in) :: iterations !! Number of algorithm iterations.
      integer, intent(in), optional :: burn !! Number of initial simulation draws to discard.
      real(dp), intent(in), optional :: offset_draws(:, :, :) !! Offset simulation draws.
      real(dp), intent(in), optional :: regression_prior_inclusion_probability(:, :) !! Regression prior inclusion probability.
      logical, intent(in), optional :: initial_regression_inclusion(:, :) !! Initial regression inclusion.
      integer, intent(in), optional :: maximum_model_size !! Maximum model size.
      integer, intent(in), optional :: maximum_flips !! Maximum flips.
      logical, intent(in), optional :: select_regression !! Flag controlling select regression.
      logical, intent(in), optional :: series_local_level !! Flag controlling series local level.
      logical, intent(in), optional :: series_local_linear_trend !! Flag controlling series local linear trend.
      integer, intent(in), optional :: series_nseasons !! Series nseasons.
      integer, intent(in), optional :: series_season_duration !! Series season duration.
      type(bsts_shared_local_level_t) :: out
      real(dp), allocatable :: regression_mean(:, :)
      real(dp), allocatable :: regression_covariance(:, :, :)
      real(dp), allocatable :: initial_mean(:), initial_covariance(:, :)
      real(dp), allocatable :: initial_loadings(:, :), loading_mean(:, :)
      real(dp), allocatable :: loading_variance(:, :), observation_scale(:)
      real(dp), allocatable :: factor_scale(:), observation_shape(:)
      real(dp), allocatable :: observation_rate(:), factor_shape(:)
      real(dp), allocatable :: factor_rate(:), state_normals(:, :, :)
      real(dp), allocatable :: loading_normals(:, :, :), gammas(:, :)
      real(dp), allocatable :: regression_inclusion_probability(:, :)
      real(dp), allocatable :: regression_normals(:, :, :)
      real(dp), allocatable :: regression_uniforms(:, :, :)
      logical, allocatable :: regression_inclusion(:, :)
      real(dp), allocatable :: series_mean(:), series_initial_scale(:)
      real(dp), allocatable :: series_scale(:), series_shape(:), series_rate(:)
      real(dp), allocatable :: series_normals(:, :, :), series_gammas(:, :)
      real(dp), allocatable :: trend_mean(:, :), trend_initial_cov(:, :, :)
      real(dp), allocatable :: trend_scale(:, :), trend_shape(:, :)
      real(dp), allocatable :: trend_rate(:, :), trend_normals(:, :, :, :)
      real(dp), allocatable :: trend_gammas(:, :, :)
      real(dp), allocatable :: seasonal_initial_scale(:), seasonal_scale(:)
      real(dp), allocatable :: seasonal_shape(:), seasonal_rate(:)
      real(dp), allocatable :: seasonal_normals(:, :, :, :)
      real(dp), allocatable :: seasonal_gammas(:, :)
      real(dp) :: largest_variance, predictor_square, response_variance
      integer :: observations, series_count, predictors_count, dimension
      integer :: discarded, iteration, series, factor, predictor, observed
      integer :: model_limit, flip_limit, seasons, season_duration
      integer :: seasonal_dimension
      logical :: use_selection, use_series_state, use_trend_state
      logical :: use_seasonal_state

      observations = size(response, 1)
      series_count = size(response, 2)
      predictors_count = size(regression_predictors, 2)
      dimension = nfactors + series_count*predictors_count
      if (observations < 2 .or. series_count < 1 .or. predictors_count < 1 .or. &
         nfactors < 1 .or. nfactors > series_count .or. iterations < 1 .or. &
         any(shape(regression_predictors) /= &
         [observations, predictors_count, series_count]) .or. &
         .not. all(ieee_is_finite(regression_predictors)) .or. &
         any([(count(ieee_is_finite(response(:, series))) < 1, &
         series=1, series_count)])) then
         out%info = 1
         return
      end if
      discarded = 0
      if (present(burn)) discarded = burn
      if (discarded < 0 .or. discarded >= iterations) then
         out%info = 1
         return
      end if
      use_selection = .true.
      if (present(select_regression)) use_selection = select_regression
      use_series_state = .false.
      if (present(series_local_level)) use_series_state = series_local_level
      use_trend_state = .false.
      if (present(series_local_linear_trend)) &
         use_trend_state = series_local_linear_trend
      seasons = 0
      if (present(series_nseasons)) seasons = series_nseasons
      use_seasonal_state = seasons >= 2
      season_duration = 1
      if (present(series_season_duration)) &
         season_duration = series_season_duration
      if ((seasons /= 0 .and. seasons < 2) .or. season_duration < 1) then
         out%info = 1
         return
      end if
      seasonal_dimension = max(1, seasons - 1)
      model_limit = predictors_count
      if (present(maximum_model_size)) model_limit = maximum_model_size
      flip_limit = -1
      if (present(maximum_flips)) flip_limit = maximum_flips
      if (model_limit < 0 .or. model_limit > predictors_count .or. &
         flip_limit < -1) then
         out%info = 1
         return
      end if
      allocate(regression_mean(predictors_count, series_count), &
         regression_covariance(predictors_count, predictors_count, &
         series_count), initial_mean(nfactors), &
         initial_covariance(nfactors, nfactors), &
         initial_loadings(series_count, nfactors), &
         loading_mean(series_count, nfactors), &
         loading_variance(series_count, nfactors), &
         observation_scale(series_count), factor_scale(nfactors), &
         observation_shape(series_count), observation_rate(series_count), &
         factor_shape(nfactors), factor_rate(nfactors), &
         state_normals(dimension, observations, iterations), &
         loading_normals(nfactors, series_count, iterations), &
         gammas(series_count + nfactors, iterations))
      if (use_selection) then
         allocate(regression_inclusion_probability(predictors_count, &
            series_count), regression_inclusion(predictors_count, &
            series_count), regression_normals(predictors_count, &
            series_count, iterations), regression_uniforms(predictors_count, &
            series_count, iterations))
         regression_inclusion_probability = 0.5_dp
         if (present(regression_prior_inclusion_probability)) then
            if (any(shape(regression_prior_inclusion_probability) /= &
               [predictors_count, series_count]) .or. &
               .not. all(ieee_is_finite( &
               regression_prior_inclusion_probability)) .or. &
               any(regression_prior_inclusion_probability < 0.0_dp) .or. &
               any(regression_prior_inclusion_probability > 1.0_dp)) then
               out%info = 1
               return
            end if
            regression_inclusion_probability = &
               regression_prior_inclusion_probability
         end if
         regression_inclusion = regression_inclusion_probability > 0.5_dp
         if (present(initial_regression_inclusion)) then
            if (any(shape(initial_regression_inclusion) /= &
               [predictors_count, series_count])) then
               out%info = 1
               return
            end if
            regression_inclusion = initial_regression_inclusion
         end if
         where (regression_inclusion_probability <= 0.0_dp) &
            regression_inclusion = .false.
         where (regression_inclusion_probability >= 1.0_dp) &
            regression_inclusion = .true.
         if (any([(count(regression_inclusion(:, series)) > model_limit, &
            series=1, series_count)])) then
            out%info = 1
            return
         end if
      else if (present(regression_prior_inclusion_probability) .or. &
         present(initial_regression_inclusion) .or. &
         present(maximum_model_size) .or. present(maximum_flips)) then
         out%info = 1
         return
      end if
      allocate(series_mean(series_count), &
            series_initial_scale(series_count), series_scale(series_count), &
            series_shape(series_count), series_rate(series_count), &
            series_normals(series_count, observations, iterations), &
            series_gammas(series_count, iterations), &
            trend_mean(2, series_count), trend_initial_cov(2, 2, series_count), &
            trend_scale(2, series_count), trend_shape(2, series_count), &
            trend_rate(2, series_count), &
            trend_normals(2, series_count, observations, iterations), &
            trend_gammas(2, series_count, iterations), &
            seasonal_initial_scale(series_count), seasonal_scale(series_count), &
            seasonal_shape(series_count), seasonal_rate(series_count), &
            seasonal_normals(seasonal_dimension, series_count, observations, &
            iterations), seasonal_gammas(series_count, iterations))
      series_mean = 0.0_dp
      trend_mean = 0.0_dp
      trend_initial_cov = 0.0_dp
      largest_variance = 1.0e-6_dp
      regression_mean = 0.0_dp
      regression_covariance = 0.0_dp
      do series = 1, series_count
         response_variance = max(finite_sample_variance(response(:, series)), &
            1.0e-6_dp)
         observation_scale(series) = max(0.1_dp*response_variance, 1.0e-8_dp)
         series_initial_scale(series) = response_variance
         series_scale(series) = max(0.01_dp*response_variance, 1.0e-10_dp)
         series_shape(series) = 2.0_dp
         series_rate(series) = series_scale(series)
         trend_initial_cov(1, 1, series) = response_variance
         trend_initial_cov(2, 2, series) = response_variance
         trend_scale(:, series) = max(0.01_dp*response_variance, 1.0e-10_dp)
         trend_shape(:, series) = 2.0_dp
         trend_rate(:, series) = trend_scale(:, series)
         seasonal_initial_scale(series) = response_variance
         seasonal_scale(series) = max(0.0001_dp*response_variance, 1.0e-10_dp)
         seasonal_shape(series) = 2.0_dp
         seasonal_rate(series) = seasonal_scale(series)
         largest_variance = max(largest_variance, response_variance)
         do predictor = 1, predictors_count
            predictor_square = sum(regression_predictors(:, predictor, &
               series)**2)/real(observations, dp)
            regression_covariance(predictor, predictor, series) = &
               10.0_dp*response_variance/ &
               max(predictor_square, 1.0e-8_dp)
         end do
      end do
      initial_mean = 0.0_dp
      initial_covariance = identity_block(nfactors, largest_variance)
      initial_loadings = 1.0_dp
      call impose_shared_loading_constraints(initial_loadings)
      loading_mean = 1.0_dp
      loading_variance = 1.0_dp
      factor_scale = max(0.01_dp*largest_variance, 1.0e-10_dp)
      observation_shape = 2.0_dp
      observation_rate = observation_scale
      factor_shape = 2.0_dp
      factor_rate = factor_scale
      do iteration = 1, iterations
         do factor = 1, dimension
            do predictor = 1, observations
               state_normals(factor, predictor, iteration) = &
                  random_standard_normal()
            end do
         end do
         do series = 1, series_count
            do factor = 1, nfactors
               loading_normals(factor, series, iteration) = &
                  random_standard_normal()
            end do
            observed = count(ieee_is_finite(response(:, series)))
            gammas(series, iteration) = random_gamma( &
               observation_shape(series) + 0.5_dp*real(observed, dp))
            if (use_selection) then
               do predictor = 1, predictors_count
                  regression_normals(predictor, series, iteration) = &
                     random_standard_normal()
                  call random_number(regression_uniforms(predictor, series, &
                     iteration))
                  regression_uniforms(predictor, series, iteration) = max( &
                     tiny(1.0_dp), min(1.0_dp - epsilon(1.0_dp), &
                     regression_uniforms(predictor, series, iteration)))
               end do
            end if
            do predictor = 1, observations
               series_normals(series, predictor, iteration) = &
                  random_standard_normal()
               trend_normals(1, series, predictor, iteration) = &
                  random_standard_normal()
               trend_normals(2, series, predictor, iteration) = &
                  random_standard_normal()
               do factor = 1, seasonal_dimension
                  seasonal_normals(factor, series, predictor, iteration) = &
                     random_standard_normal()
               end do
            end do
            series_gammas(series, iteration) = random_gamma( &
               series_shape(series) + 0.5_dp*real(observations - 1, dp))
            do factor = 1, 2
               trend_gammas(factor, series, iteration) = random_gamma( &
                  trend_shape(factor, series) + &
                  0.5_dp*real(observations - 1, dp))
            end do
            seasonal_gammas(series, iteration) = random_gamma( &
               seasonal_shape(series) + 0.5_dp*real( &
               (observations - 1)/season_duration, dp))
         end do
         do factor = 1, nfactors
            gammas(series_count + factor, iteration) = random_gamma( &
               factor_shape(factor) + 0.5_dp*real(observations - 1, dp))
         end do
      end do
      if (use_selection) then
         out = bsts_mbsts_draws(response, regression_predictors, &
            regression_mean, regression_covariance, initial_mean, &
            initial_covariance, initial_loadings, loading_mean, &
            loading_variance, observation_scale, factor_scale, &
            observation_shape, observation_rate, factor_shape, factor_rate, &
            discarded, state_normals, loading_normals, gammas, offset_draws, &
            regression_inclusion_probability, regression_inclusion, &
            regression_normals, regression_uniforms, model_limit, flip_limit, &
            series_mean, series_initial_scale, series_scale, series_shape, &
            series_rate, series_normals, series_gammas, trend_mean, &
            trend_initial_cov, trend_scale, trend_shape, trend_rate, &
            trend_normals, trend_gammas, seasons, season_duration, &
            seasonal_initial_scale, seasonal_scale, seasonal_shape, &
            seasonal_rate, seasonal_normals, seasonal_gammas, &
            use_series_state, use_trend_state, use_seasonal_state)
      else
         out = bsts_mbsts_draws(response, regression_predictors, &
            regression_mean, regression_covariance, initial_mean, &
            initial_covariance, initial_loadings, loading_mean, &
            loading_variance, observation_scale, factor_scale, &
            observation_shape, observation_rate, factor_shape, factor_rate, &
            discarded, state_normals, loading_normals, gammas, offset_draws, &
            series_initial_mean=series_mean, &
            series_initial_variance=series_initial_scale, &
            series_variance=series_scale, series_prior_shape=series_shape, &
            series_prior_rate=series_rate, &
            series_state_normal_draws=series_normals, &
            series_gamma_draws=series_gammas, &
            trend_initial_mean=trend_mean, &
            trend_initial_covariance=trend_initial_cov, &
            trend_variance=trend_scale, trend_prior_shape=trend_shape, &
            trend_prior_rate=trend_rate, &
            trend_state_normal_draws=trend_normals, &
            trend_gamma_draws=trend_gammas, seasonal_nseasons=seasons, &
            seasonal_duration=season_duration, &
            seasonal_initial_variance=seasonal_initial_scale, &
            seasonal_variance=seasonal_scale, &
            seasonal_prior_shape=seasonal_shape, &
            seasonal_prior_rate=seasonal_rate, &
            seasonal_state_normal_draws=seasonal_normals, &
            seasonal_gamma_draws=seasonal_gammas, &
            enable_series_state=use_series_state, &
            enable_trend_state=use_trend_state, &
            enable_seasonal_state=use_seasonal_state)
      end if
   end function bsts_mbsts

   pure function bsts_mbsts_predict_draws(fit, future_predictors, &
      state_normal_draws, observation_normal_draws, &
      series_state_normal_draws, trend_state_normal_draws, &
      seasonal_state_normal_draws) result(out)
      !! Forecast shared, regression, and series-specific state contributions.
      type(bsts_shared_local_level_t), intent(in) :: fit !! Previously fitted model.
      real(dp), intent(in) :: future_predictors(:, :, :) !! Future predictors.
      real(dp), intent(in) :: state_normal_draws(:, :, :) !! Independent standard-normal state draws.
      real(dp), intent(in) :: observation_normal_draws(:, :, :) !! Observation normal draws.
      real(dp), intent(in), optional :: series_state_normal_draws(:, :, :) !! Series state normal simulation draws.
      real(dp), intent(in), optional :: trend_state_normal_draws(:, :, :, :) !! Trend state normal simulation draws.
      real(dp), intent(in), optional :: seasonal_state_normal_draws(:, :, :, :) !! Seasonal state normal simulation draws.
      type(bsts_multivariate_prediction_t) :: out
      real(dp), allocatable :: sorted(:)
      integer :: horizon, predictors_count, series_count, retained
      integer :: draw, source, step, series
      real(dp) :: series_level
      real(dp) :: trend_state(2)
      real(dp), allocatable :: seasonal_state(:), next_seasonal_state(:)

      if (.not. fit%is_mbsts .or. &
         .not. allocated(fit%regression_coefficients)) then
         out%info = 1
         return
      end if
      horizon = size(future_predictors, 1)
      predictors_count = size(fit%regression_coefficients, 1)
      series_count = size(fit%loadings, 1)
      retained = size(fit%state, 3) - fit%burn
      if (horizon < 1 .or. retained < 1 .or. &
         any(shape(future_predictors) /= &
         [horizon, predictors_count, series_count]) .or. &
         .not. all(ieee_is_finite(future_predictors)) .or. &
         (fit%has_series_local_level .neqv. &
         present(series_state_normal_draws)) .or. &
         (fit%has_series_local_linear_trend .neqv. &
         present(trend_state_normal_draws)) .or. &
         (fit%has_series_seasonal .neqv. &
         present(seasonal_state_normal_draws))) then
         out%info = 1
         return
      end if
      if (fit%has_series_local_linear_trend) then
         if (.not. allocated(fit%series_trend_state) .or. &
            .not. allocated(fit%series_trend_variance) .or. &
            any(shape(trend_state_normal_draws) /= &
            [2, series_count, horizon, retained]) .or. &
            .not. all(ieee_is_finite(trend_state_normal_draws))) then
            out%info = 1
            return
         end if
      end if
      if (fit%has_series_seasonal) then
         if (.not. allocated(fit%series_seasonal_state) .or. &
            .not. allocated(fit%series_seasonal_variance) .or. &
            any(shape(seasonal_state_normal_draws) /= &
            [fit%series_seasons - 1, series_count, horizon, retained]) .or. &
            .not. all(ieee_is_finite(seasonal_state_normal_draws))) then
            out%info = 1
            return
         end if
         allocate(seasonal_state(fit%series_seasons - 1), &
            next_seasonal_state(fit%series_seasons - 1))
      end if
      if (fit%has_series_local_level) then
         if (.not. allocated(fit%series_state) .or. &
            .not. allocated(fit%series_variance) .or. &
            any(shape(series_state_normal_draws) /= &
            [series_count, horizon, retained]) .or. &
            .not. all(ieee_is_finite(series_state_normal_draws))) then
            out%info = 1
            return
         end if
      end if
      out = bsts_shared_local_level_predict_draws(fit, horizon, &
         state_normal_draws, observation_normal_draws)
      if (out%info /= 0) return
      do draw = 1, retained
         source = fit%burn + draw
         do series = 1, series_count
            if (fit%has_series_local_level) &
               series_level = fit%series_state(series, &
               size(fit%series_state, 2), source)
            if (fit%has_series_local_linear_trend) &
               trend_state = fit%series_trend_state(:, series, &
               size(fit%series_trend_state, 3), source)
            if (fit%has_series_seasonal) seasonal_state = &
               fit%series_seasonal_state(:, series, &
               size(fit%series_seasonal_state, 3), source)
            do step = 1, horizon
               out%draws(series, step, draw) = &
                  out%draws(series, step, draw) + &
                  dot_product(future_predictors(step, :, series), &
                  fit%regression_coefficients(:, series, source))
               if (fit%has_series_local_level) then
                  series_level = series_level + &
                     sqrt(fit%series_variance(series, source))* &
                     series_state_normal_draws(series, step, draw)
                  out%draws(series, step, draw) = &
                     out%draws(series, step, draw) + series_level
               end if
               if (fit%has_series_local_linear_trend) then
                  trend_state(1) = trend_state(1) + trend_state(2) + &
                     sqrt(fit%series_trend_variance(1, series, source))* &
                     trend_state_normal_draws(1, series, step, draw)
                  trend_state(2) = trend_state(2) + &
                     sqrt(fit%series_trend_variance(2, series, source))* &
                     trend_state_normal_draws(2, series, step, draw)
                  out%draws(series, step, draw) = &
                     out%draws(series, step, draw) + trend_state(1)
               end if
               if (fit%has_series_seasonal) then
                  if (modulo(size(fit%series_seasonal_state, 3) + step - 1, &
                     fit%series_season_duration) == 0) then
                     next_seasonal_state = 0.0_dp
                     next_seasonal_state(1) = -sum(seasonal_state) + &
                        sqrt(fit%series_seasonal_variance(series, source))* &
                        seasonal_state_normal_draws(1, series, step, draw)
                     if (size(seasonal_state) > 1) &
                        next_seasonal_state(2:) = &
                        seasonal_state(:size(seasonal_state) - 1)
                     seasonal_state = next_seasonal_state
                  end if
                  out%draws(series, step, draw) = &
                     out%draws(series, step, draw) + seasonal_state(1)
               end if
            end do
         end do
      end do
      do series = 1, series_count
         do step = 1, horizon
            out%mean(series, step) = sum(out%draws(series, step, :))/ &
               real(retained, dp)
            if (retained > 1) then
               out%standard_deviation(series, step) = sqrt(sum( &
                  (out%draws(series, step, :) - out%mean(series, step))**2)/ &
                  real(retained - 1, dp))
            else
               out%standard_deviation(series, step) = 0.0_dp
            end if
            sorted = out%draws(series, step, :)
            call insertion_sort(sorted)
            out%lower(series, step) = quantile(sorted, 0.025_dp)
            out%upper(series, step) = quantile(sorted, 0.975_dp)
         end do
      end do
   end function bsts_mbsts_predict_draws

   function bsts_mbsts_predict(fit, future_predictors) result(out)
      !! Forecast multivariate BSTS regression using the shared random stream.
      type(bsts_shared_local_level_t), intent(in) :: fit !! Previously fitted model.
      real(dp), intent(in) :: future_predictors(:, :, :) !! Future predictors.
      type(bsts_multivariate_prediction_t) :: out
      real(dp), allocatable :: state_normals(:, :, :)
      real(dp), allocatable :: observation_normals(:, :, :)
      real(dp), allocatable :: series_normals(:, :, :)
      real(dp), allocatable :: trend_normals(:, :, :, :)
      real(dp), allocatable :: seasonal_normals(:, :, :, :)
      integer :: factors, series_count, horizon, retained
      integer :: draw, step, component, seasonal_component

      if (.not. allocated(fit%state) .or. .not. allocated(fit%loadings)) then
         out%info = 1
         return
      end if
      factors = size(fit%state, 1)
      series_count = size(fit%loadings, 1)
      horizon = size(future_predictors, 1)
      retained = size(fit%state, 3) - fit%burn
      if (horizon < 1 .or. retained < 1) then
         out%info = 1
         return
      end if
      allocate(state_normals(factors, horizon, retained), &
         observation_normals(series_count, horizon, retained))
      if (fit%has_series_local_level) &
         allocate(series_normals(series_count, horizon, retained))
      if (fit%has_series_local_linear_trend) &
         allocate(trend_normals(2, series_count, horizon, retained))
      if (fit%has_series_seasonal) &
         allocate(seasonal_normals(fit%series_seasons - 1, series_count, &
         horizon, retained))
      do draw = 1, retained
         do step = 1, horizon
            do component = 1, factors
               state_normals(component, step, draw) = random_standard_normal()
            end do
            do component = 1, series_count
               observation_normals(component, step, draw) = &
                  random_standard_normal()
               if (fit%has_series_local_level) &
                  series_normals(component, step, draw) = &
                  random_standard_normal()
               if (fit%has_series_local_linear_trend) then
                  trend_normals(1, component, step, draw) = &
                     random_standard_normal()
                  trend_normals(2, component, step, draw) = &
                     random_standard_normal()
               end if
               if (fit%has_series_seasonal) then
                  do seasonal_component = 1, fit%series_seasons - 1
                     seasonal_normals(seasonal_component, component, step, &
                        draw) = &
                        random_standard_normal()
                  end do
               end if
            end do
         end do
      end do
      if (fit%has_series_local_level .and. &
         fit%has_series_local_linear_trend .and. fit%has_series_seasonal) then
         out = bsts_mbsts_predict_draws(fit, future_predictors, state_normals, &
            observation_normals, series_normals, trend_normals, &
            seasonal_normals)
      else if (fit%has_series_local_level .and. &
         fit%has_series_local_linear_trend) then
         out = bsts_mbsts_predict_draws(fit, future_predictors, state_normals, &
            observation_normals, series_normals, trend_normals)
      else if (fit%has_series_local_level .and. fit%has_series_seasonal) then
         out = bsts_mbsts_predict_draws(fit, future_predictors, state_normals, &
            observation_normals, series_state_normal_draws=series_normals, &
            seasonal_state_normal_draws=seasonal_normals)
      else if (fit%has_series_local_linear_trend .and. &
         fit%has_series_seasonal) then
         out = bsts_mbsts_predict_draws(fit, future_predictors, state_normals, &
            observation_normals, trend_state_normal_draws=trend_normals, &
            seasonal_state_normal_draws=seasonal_normals)
      else if (fit%has_series_local_level) then
         out = bsts_mbsts_predict_draws(fit, future_predictors, state_normals, &
            observation_normals, series_normals)
      else if (fit%has_series_local_linear_trend) then
         out = bsts_mbsts_predict_draws(fit, future_predictors, state_normals, &
            observation_normals, trend_state_normal_draws=trend_normals)
      else if (fit%has_series_seasonal) then
         out = bsts_mbsts_predict_draws(fit, future_predictors, state_normals, &
            observation_normals, &
            seasonal_state_normal_draws=seasonal_normals)
      else
         out = bsts_mbsts_predict_draws(fit, future_predictors, state_normals, &
            observation_normals)
      end if
   end function bsts_mbsts_predict

   pure function bsts_trig_draws(y, period, frequencies, initial_mean, &
      initial_covariance, observation_variance, trig_variance, &
      observation_prior_shape, observation_prior_rate, trig_prior_shape, &
      trig_prior_rate, burn, state_normal_draws, gamma_draws) result(out)
      !! Draw a harmonic trigonometric posterior from supplied random variates.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: period !! Seasonal period.
      real(dp), intent(in) :: frequencies(:) !! Frequencies.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: observation_variance !! Observation-error variance.
      real(dp), intent(in) :: trig_variance !! Trig variance.
      real(dp), intent(in) :: observation_prior_shape !! Observation prior shape.
      real(dp), intent(in) :: observation_prior_rate !! Observation prior rate.
      real(dp), intent(in) :: trig_prior_shape !! Trig prior shape.
      real(dp), intent(in) :: trig_prior_rate !! Trig prior rate.
      integer, intent(in) :: burn !! Number of initial simulation draws to discard.
      real(dp), intent(in) :: state_normal_draws(:, :, :) !! Independent standard-normal state draws.
      real(dp), intent(in) :: gamma_draws(:, :) !! Gamma simulation draws.
      type(bsts_mcmc_t) :: out
      type(ssm_model_t) :: model
      real(dp), allocatable :: path(:, :), transition(:, :), state_scale(:)
      real(dp), allocatable :: previous(:), disturbance(:)
      real(dp) :: observation_scale, shared_scale, observation_sum, state_sum
      integer :: dimension, iterations, iteration, time, component, status

      dimension = 2*size(frequencies)
      iterations = size(state_normal_draws, 3)
      if (size(y) < 2 .or. period <= 0.0_dp .or. size(frequencies) < 1 .or. &
         any(frequencies <= 0.0_dp) .or. size(initial_mean) /= dimension .or. &
         any(shape(initial_covariance) /= [dimension, dimension]) .or. &
         observation_variance <= 0.0_dp .or. trig_variance <= 0.0_dp .or. &
         observation_prior_shape <= 0.0_dp .or. observation_prior_rate <= 0.0_dp .or. &
         trig_prior_shape <= 0.0_dp .or. trig_prior_rate <= 0.0_dp .or. &
         burn < 0 .or. burn >= iterations .or. &
         size(state_normal_draws, 1) /= dimension .or. &
         size(state_normal_draws, 2) /= size(y) .or. &
         any(shape(gamma_draws) /= [2, iterations]) .or. &
         any(gamma_draws <= 0.0_dp) .or. &
         .not. all(ieee_is_finite(state_normal_draws))) then
         out%info = 1
         return
      end if
      transition = trig_transition_matrix(period, frequencies)
      allocate(state_scale(dimension), previous(dimension), disturbance(dimension), &
         out%state(dimension, size(y), iterations), &
         out%observation_variance(iterations), &
         out%state_variance(dimension, iterations), &
         out%component_variance(1, iterations), &
         out%transition(dimension, dimension), out%observation(dimension), &
         out%transition_schedule(dimension, dimension, 1), &
         out%state_loading_schedule(dimension, dimension, 1))
      out%transition = transition
      out%transition_schedule(:, :, 1) = transition
      out%observation = 0.0_dp
      out%state_loading_schedule = 0.0_dp
      do component = 1, size(frequencies)
         out%observation(2*component - 1) = 1.0_dp
      end do
      do component = 1, dimension
         out%state_loading_schedule(component, component, 1) = 1.0_dp
      end do
      out%burn = burn
      observation_scale = observation_variance
      shared_scale = trig_variance
      do iteration = 1, iterations
         state_scale = shared_scale
         model = structural_model(y, initial_mean, initial_covariance, &
            transition, out%observation, observation_scale, state_scale)
         call state_path_from_draws(model, state_normal_draws(:, :, iteration), &
            path, status)
         if (status /= 0) then
            out%info = 10 + status
            return
         end if
         out%state(:, :, iteration) = path
         observation_sum = 0.0_dp
         do time = 1, size(y)
            if (ieee_is_finite(y(time))) observation_sum = observation_sum + &
               (y(time) - dot_product(out%observation, path(:, time)))**2
         end do
         observation_scale = (observation_prior_rate + &
            0.5_dp*observation_sum)/gamma_draws(1, iteration)
         state_sum = 0.0_dp
         do time = 2, size(y)
            previous = matmul(transition, path(:, time - 1))
            disturbance = path(:, time) - previous
            state_sum = state_sum + sum(disturbance**2)
         end do
         shared_scale = (trig_prior_rate + 0.5_dp*state_sum)/ &
            gamma_draws(2, iteration)
         out%observation_variance(iteration) = observation_scale
         out%state_variance(:, iteration) = shared_scale
         out%component_variance(1, iteration) = shared_scale
      end do
   end function bsts_trig_draws

   function bsts_trig(y, period, frequencies, iterations, burn, &
      initial_mean, initial_covariance, observation_variance, trig_variance, &
      observation_prior_shape, observation_prior_rate, trig_prior_shape, &
      trig_prior_rate) result(out)
      !! Sample a harmonic trigonometric posterior using the shared random stream.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: period !! Seasonal period.
      real(dp), intent(in) :: frequencies(:) !! Frequencies.
      integer, intent(in) :: iterations !! Number of algorithm iterations.
      integer, intent(in), optional :: burn !! Number of initial simulation draws to discard.
      real(dp), intent(in), optional :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in), optional :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in), optional :: observation_variance !! Observation-error variance.
      real(dp), intent(in), optional :: trig_variance !! Trig variance.
      real(dp), intent(in), optional :: observation_prior_shape !! Observation prior shape.
      real(dp), intent(in), optional :: observation_prior_rate !! Observation prior rate.
      real(dp), intent(in), optional :: trig_prior_shape !! Trig prior shape.
      real(dp), intent(in), optional :: trig_prior_rate !! Trig prior rate.
      type(bsts_mcmc_t) :: out
      real(dp), allocatable :: mean0(:), covariance0(:, :)
      real(dp), allocatable :: normal_draws(:, :, :), gamma_draws(:, :)
      real(dp) :: observation0, trig0, observation_shape, observation_rate
      real(dp) :: trig_shape, trig_rate, series_variance
      integer :: dimension, discarded, iteration, observed, component

      dimension = 2*size(frequencies)
      discarded = 0
      if (present(burn)) discarded = burn
      if (size(y) < 2 .or. period <= 0.0_dp .or. size(frequencies) < 1 .or. &
         any(frequencies <= 0.0_dp) .or. iterations < 1 .or. &
         discarded < 0 .or. discarded >= iterations) then
         out%info = 1
         return
      end if
      series_variance = finite_sample_variance(y)
      allocate(mean0(dimension), covariance0(dimension, dimension))
      mean0 = 0.0_dp
      if (present(initial_mean)) then
         if (size(initial_mean) /= dimension) then
            out%info = 1
            return
         end if
         mean0 = initial_mean
      end if
      covariance0 = 0.0_dp
      do component = 1, dimension
         covariance0(component, component) = max(series_variance, 1.0_dp)
      end do
      if (present(initial_covariance)) then
         if (any(shape(initial_covariance) /= [dimension, dimension])) then
            out%info = 1
            return
         end if
         covariance0 = initial_covariance
      end if
      observation0 = max(0.5_dp*series_variance, 1.0e-6_dp)
      if (present(observation_variance)) observation0 = observation_variance
      trig0 = max(0.0001_dp*series_variance, 1.0e-10_dp)
      if (present(trig_variance)) trig0 = trig_variance
      observation_shape = 0.01_dp
      if (present(observation_prior_shape)) &
         observation_shape = observation_prior_shape
      observation_rate = max(0.01_dp*series_variance, 1.0e-8_dp)
      if (present(observation_prior_rate)) observation_rate = &
         observation_prior_rate
      trig_shape = 0.01_dp
      if (present(trig_prior_shape)) trig_shape = trig_prior_shape
      trig_rate = max(0.0001_dp*series_variance, 1.0e-10_dp)
      if (present(trig_prior_rate)) trig_rate = trig_prior_rate
      if (observation0 <= 0.0_dp .or. trig0 <= 0.0_dp .or. &
         observation_shape <= 0.0_dp .or. observation_rate <= 0.0_dp .or. &
         trig_shape <= 0.0_dp .or. trig_rate <= 0.0_dp) then
         out%info = 1
         return
      end if
      allocate(normal_draws(dimension, size(y), iterations), &
         gamma_draws(2, iterations))
      observed = count(ieee_is_finite(y))
      do iteration = 1, iterations
         call fill_standard_normals(normal_draws(:, :, iteration))
         gamma_draws(1, iteration) = random_gamma( &
            observation_shape + 0.5_dp*real(observed, dp))
         gamma_draws(2, iteration) = random_gamma(trig_shape + &
            0.5_dp*real(dimension*(size(y) - 1), dp))
      end do
      out = bsts_trig_draws(y, period, frequencies, mean0, covariance0, &
         observation0, trig0, observation_shape, observation_rate, trig_shape, &
         trig_rate, discarded, normal_draws, gamma_draws)
   end function bsts_trig

   pure function bsts_spike_slab_draws(response, predictors, slab_mean, &
      slab_covariance, prior_inclusion_probability, variance_prior_shape, &
      variance_prior_rate, initial_variance, initial_inclusion, &
      maximum_model_size, maximum_flips, burn, coefficient_normal_draws, &
      inclusion_uniform_draws, gamma_draws, offset_draws) result(out)
      !! Sample static Gaussian spike-and-slab regression from supplied draws.
      real(dp), intent(in) :: response(:) !! Response observations.
      real(dp), intent(in) :: predictors(:, :) !! Predictor matrix.
      real(dp), intent(in) :: slab_mean(:) !! Slab mean.
      real(dp), intent(in) :: slab_covariance(:, :) !! Slab covariance matrix.
      real(dp), intent(in) :: prior_inclusion_probability(:) !! Prior inclusion probability.
      real(dp), intent(in) :: variance_prior_shape !! Variance prior shape.
      real(dp), intent(in) :: variance_prior_rate !! Variance prior rate.
      real(dp), intent(in) :: initial_variance !! Initial variance.
      logical, intent(in) :: initial_inclusion(:) !! Initial inclusion.
      integer, intent(in) :: maximum_model_size !! Maximum model size.
      integer, intent(in) :: maximum_flips !! Maximum flips.
      integer, intent(in) :: burn !! Number of initial simulation draws to discard.
      real(dp), intent(in) :: coefficient_normal_draws(:, :) !! Coefficient normal simulation draws.
      real(dp), intent(in) :: inclusion_uniform_draws(:, :) !! Inclusion uniform simulation draws.
      real(dp), intent(in) :: gamma_draws(:) !! Gamma simulation draws.
      real(dp), intent(in), optional :: offset_draws(:, :) !! Offset simulation draws.
      type(bsts_spike_slab_t) :: out
      real(dp), allocatable :: observed_x(:, :), observed_y(:), working_y(:)
      real(dp), allocatable :: posterior_mean(:), posterior_covariance(:, :)
      real(dp), allocatable :: coefficient_draw(:), active_draw(:)
      integer, allocatable :: observed_index(:), active(:)
      logical, allocatable :: included(:), candidate(:)
      real(dp) :: variance, log_included, log_excluded, probability
      real(dp) :: log_probability, collapsed_rate
      integer :: observations, variables, iterations, iteration, variable
      integer :: flips, status, retained

      observations = count(ieee_is_finite(response))
      variables = size(predictors, 2)
      iterations = size(coefficient_normal_draws, 2)
      if (size(response) < 1 .or. observations < 1 .or. variables < 1 .or. &
         size(predictors, 1) /= size(response) .or. &
         .not. all(ieee_is_finite(predictors)) .or. &
         size(slab_mean) /= variables .or. &
         any(shape(slab_covariance) /= [variables, variables]) .or. &
         size(prior_inclusion_probability) /= variables .or. &
         any(prior_inclusion_probability < 0.0_dp) .or. &
         any(prior_inclusion_probability > 1.0_dp) .or. &
         variance_prior_shape <= 0.0_dp .or. variance_prior_rate <= 0.0_dp .or. &
         initial_variance <= 0.0_dp .or. size(initial_inclusion) /= variables .or. &
         maximum_model_size < 0 .or. maximum_model_size > variables .or. &
         maximum_flips < -1 .or. burn < 0 .or. burn >= iterations .or. &
         size(coefficient_normal_draws, 1) /= variables .or. &
         any(shape(inclusion_uniform_draws) /= [variables, iterations]) .or. &
         size(gamma_draws) /= iterations .or. any(gamma_draws <= 0.0_dp) .or. &
         any(inclusion_uniform_draws <= 0.0_dp) .or. &
         any(inclusion_uniform_draws >= 1.0_dp) .or. &
         .not. all(ieee_is_finite(coefficient_normal_draws))) then
         out%info = 1
         return
      end if
      if (present(offset_draws)) then
         if (any(shape(offset_draws) /= [size(response), iterations]) .or. &
            .not. all(ieee_is_finite(offset_draws))) then
            out%info = 1
            return
         end if
      end if
      observed_index = pack([(variable, variable=1, size(response))], &
         ieee_is_finite(response))
      allocate(observed_x(observations, variables), observed_y(observations), &
         working_y(observations), included(variables), candidate(variables), &
         coefficient_draw(variables), out%coefficients(variables, iterations), &
         out%included(variables, iterations), out%residual_variance(iterations), &
         out%log_model_probability(iterations), &
         out%inclusion_probability(variables))
      observed_y = response(observed_index)
      observed_x = predictors(observed_index, :)
      included = initial_inclusion
      where (prior_inclusion_probability <= 0.0_dp) included = .false.
      where (prior_inclusion_probability >= 1.0_dp) included = .true.
      if (count(included) > maximum_model_size) then
         out%info = 2
         return
      end if
      variance = initial_variance
      out%burn = burn
      do iteration = 1, iterations
         working_y = observed_y
         if (present(offset_draws)) working_y = working_y - &
            offset_draws(observed_index, iteration)
         flips = 0
         do variable = 1, variables
            if (prior_inclusion_probability(variable) <= 0.0_dp) then
               included(variable) = .false.
               cycle
            end if
            if (prior_inclusion_probability(variable) >= 1.0_dp) then
               included(variable) = .true.
               cycle
            end if
            candidate = included
            candidate(variable) = .false.
            call spike_slab_model_moments(working_y, observed_x, candidate, &
               slab_mean, slab_covariance, prior_inclusion_probability, &
               variance, variance_prior_rate, log_excluded, collapsed_rate, &
               posterior_mean, posterior_covariance, active, status)
            if (status /= 0) then
               out%info = 10 + status
               return
            end if
            candidate(variable) = .true.
            if (count(candidate) > maximum_model_size) then
               probability = 0.0_dp
            else
               call spike_slab_model_moments(working_y, observed_x, candidate, &
                  slab_mean, slab_covariance, prior_inclusion_probability, &
                  variance, variance_prior_rate, log_included, collapsed_rate, &
                  posterior_mean, posterior_covariance, active, status)
               if (status /= 0) then
                  out%info = 20 + status
                  return
               end if
               probability = logistic_log_odds(log_included - log_excluded)
            end if
            candidate(variable) = inclusion_uniform_draws(variable, iteration) < &
               probability
            if (candidate(variable) .neqv. included(variable)) then
               if (maximum_flips < 0 .or. flips < maximum_flips) then
                  included(variable) = candidate(variable)
                  flips = flips + 1
               end if
            end if
         end do
         call spike_slab_model_moments(working_y, observed_x, included, &
            slab_mean, slab_covariance, prior_inclusion_probability, variance, &
            variance_prior_rate, log_probability, collapsed_rate, &
            posterior_mean, posterior_covariance, active, status)
         if (status /= 0) then
            out%info = 30 + status
            return
         end if
         variance = collapsed_rate/gamma_draws(iteration)
         coefficient_draw = 0.0_dp
         if (size(active) > 0) then
            allocate(active_draw(size(active)))
            call multivariate_normal_from_standard(posterior_mean, &
               variance*posterior_covariance, &
               coefficient_normal_draws(active, iteration), &
               active_draw, status)
            if (status /= 0) then
               out%info = 40 + status
               return
            end if
            coefficient_draw(active) = active_draw
            deallocate(active_draw)
         end if
         out%coefficients(:, iteration) = coefficient_draw
         out%included(:, iteration) = included
         out%residual_variance(iteration) = variance
         out%log_model_probability(iteration) = log_probability
      end do
      retained = iterations - burn
      do variable = 1, variables
         out%inclusion_probability(variable) = real(count( &
            out%included(variable, burn + 1:)), dp)/real(retained, dp)
      end do
   end function bsts_spike_slab_draws

   function bsts_spike_slab(response, predictors, iterations, burn, &
      slab_mean, slab_covariance, prior_inclusion_probability, &
      expected_model_size, variance_prior_shape, variance_prior_rate, &
      initial_variance, initial_inclusion, maximum_model_size, maximum_flips, &
      offset_draws) result(out)
      !! Sample static spike-and-slab regression using the shared random stream.
      real(dp), intent(in) :: response(:) !! Response observations.
      real(dp), intent(in) :: predictors(:, :) !! Predictor matrix.
      integer, intent(in) :: iterations !! Number of algorithm iterations.
      integer, intent(in), optional :: burn !! Number of initial simulation draws to discard.
      integer, intent(in), optional :: maximum_model_size !! Maximum model size.
      integer, intent(in), optional :: maximum_flips !! Maximum flips.
      real(dp), intent(in), optional :: slab_mean(:) !! Slab mean.
      real(dp), intent(in), optional :: slab_covariance(:, :) !! Slab covariance matrix.
      real(dp), intent(in), optional :: prior_inclusion_probability(:) !! Prior inclusion probability.
      real(dp), intent(in), optional :: expected_model_size !! Expected model size.
      real(dp), intent(in), optional :: variance_prior_shape !! Variance prior shape.
      real(dp), intent(in), optional :: variance_prior_rate !! Variance prior rate.
      real(dp), intent(in), optional :: initial_variance !! Initial variance.
      logical, intent(in), optional :: initial_inclusion(:) !! Initial inclusion.
      real(dp), intent(in), optional :: offset_draws(:, :) !! Offset simulation draws.
      type(bsts_spike_slab_t) :: out
      real(dp), allocatable :: mean0(:), covariance0(:, :), inclusion0(:)
      real(dp), allocatable :: normal_draws(:, :), uniform_draws(:, :)
      real(dp), allocatable :: gamma_draws(:), xtx(:, :), inverse(:, :)
      logical, allocatable :: model0(:)
      real(dp) :: expected, shape0, rate0, variance0, series_variance
      integer :: variables, discarded, model_limit, flip_limit
      integer :: iteration, variable, status, observations

      variables = size(predictors, 2)
      if (size(response) < 1 .or. size(predictors, 1) /= size(response) .or. &
         variables < 1 .or. iterations < 1 .or. &
         .not. all(ieee_is_finite(predictors))) then
         out%info = 1
         return
      end if
      discarded = 0
      if (present(burn)) discarded = burn
      model_limit = variables
      if (present(maximum_model_size)) model_limit = maximum_model_size
      flip_limit = -1
      if (present(maximum_flips)) flip_limit = maximum_flips
      series_variance = finite_sample_variance(response)
      allocate(mean0(variables), covariance0(variables, variables), &
         inclusion0(variables), model0(variables))
      mean0 = 0.0_dp
      if (present(slab_mean)) then
         if (size(slab_mean) /= variables) then
            out%info = 1
            return
         end if
         mean0 = slab_mean
      end if
      observations = count(ieee_is_finite(response))
      xtx = matmul(transpose(pack_predictors(response, predictors)), &
         pack_predictors(response, predictors))
      call invert_matrix(xtx, inverse, status)
      if (status == 0) then
         covariance0 = real(max(1, observations), dp)*inverse
      else
         covariance0 = 0.0_dp
         do variable = 1, variables
            covariance0(variable, variable) = &
               real(max(1, observations), dp)/max(xtx(variable, variable), &
               1.0_dp)
         end do
      end if
      if (present(slab_covariance)) then
         if (any(shape(slab_covariance) /= [variables, variables])) then
            out%info = 1
            return
         end if
         covariance0 = slab_covariance
      end if
      expected = real(min(3, variables), dp)
      if (present(expected_model_size)) expected = expected_model_size
      inclusion0 = max(0.0_dp, min(1.0_dp, expected/real(variables, dp)))
      if (present(prior_inclusion_probability)) then
         if (size(prior_inclusion_probability) /= variables) then
            out%info = 1
            return
         end if
         inclusion0 = prior_inclusion_probability
      end if
      model0 = inclusion0 >= 1.0_dp
      if (present(initial_inclusion)) then
         if (size(initial_inclusion) /= variables) then
            out%info = 1
            return
         end if
         model0 = initial_inclusion
      end if
      shape0 = 0.01_dp
      if (present(variance_prior_shape)) shape0 = variance_prior_shape
      rate0 = max(0.01_dp*series_variance, 1.0e-8_dp)
      if (present(variance_prior_rate)) rate0 = variance_prior_rate
      variance0 = max(series_variance, 1.0e-6_dp)
      if (present(initial_variance)) variance0 = initial_variance
      if (discarded < 0 .or. discarded >= iterations .or. expected < 0.0_dp .or. &
         model_limit < 0 .or. model_limit > variables .or. flip_limit < -1 .or. &
         shape0 <= 0.0_dp .or. rate0 <= 0.0_dp .or. variance0 <= 0.0_dp) then
         out%info = 1
         return
      end if
      allocate(normal_draws(variables, iterations), &
         uniform_draws(variables, iterations), gamma_draws(iterations))
      do iteration = 1, iterations
         do variable = 1, variables
            normal_draws(variable, iteration) = random_standard_normal()
            uniform_draws(variable, iteration) = max(tiny(1.0_dp), &
               min(1.0_dp - epsilon(1.0_dp), random_uniform()))
         end do
         gamma_draws(iteration) = random_gamma(shape0 + &
            0.5_dp*real(observations, dp))
      end do
      if (present(offset_draws)) then
         out = bsts_spike_slab_draws(response, predictors, mean0, covariance0, &
            inclusion0, shape0, rate0, variance0, model0, model_limit, &
            flip_limit, discarded, normal_draws, uniform_draws, gamma_draws, &
            offset_draws)
      else
         out = bsts_spike_slab_draws(response, predictors, mean0, covariance0, &
            inclusion0, shape0, rate0, variance0, model0, model_limit, &
            flip_limit, discarded, normal_draws, uniform_draws, gamma_draws)
      end if
   end function bsts_spike_slab

   pure function bsts_regression_predict_draws(fit, future_predictors, &
      observation_normal_draws) result(out)
      !! Form static-regression posterior predictions from supplied normals.
      type(bsts_spike_slab_t), intent(in) :: fit !! Previously fitted model.
      real(dp), intent(in) :: future_predictors(:, :) !! Future predictors.
      real(dp), intent(in) :: observation_normal_draws(:, :) !! Observation normal draws.
      type(bsts_prediction_t) :: out
      real(dp), allocatable :: sorted(:)
      integer :: horizon, retained, draw, source, step

      if (fit%info /= 0 .or. .not. allocated(fit%coefficients) .or. &
         .not. allocated(fit%residual_variance)) then
         out%info = 1
         return
      end if
      horizon = size(future_predictors, 1)
      retained = size(fit%coefficients, 2) - fit%burn
      if (horizon < 1 .or. retained < 1 .or. &
         size(future_predictors, 2) /= size(fit%coefficients, 1) .or. &
         any(shape(observation_normal_draws) /= [horizon, retained]) .or. &
         .not. all(ieee_is_finite(future_predictors)) .or. &
         .not. all(ieee_is_finite(observation_normal_draws))) then
         out%info = 1
         return
      end if
      allocate(out%draws(horizon, retained), out%mean(horizon), &
         out%standard_deviation(horizon), out%lower(horizon), &
         out%upper(horizon))
      do draw = 1, retained
         source = fit%burn + draw
         out%draws(:, draw) = matmul(future_predictors, &
            fit%coefficients(:, source)) + sqrt(fit%residual_variance(source))* &
            observation_normal_draws(:, draw)
      end do
      do step = 1, horizon
         out%mean(step) = sum(out%draws(step, :))/real(retained, dp)
         if (retained > 1) then
            out%standard_deviation(step) = sqrt(sum((out%draws(step, :) - &
               out%mean(step))**2)/real(retained - 1, dp))
         else
            out%standard_deviation(step) = 0.0_dp
         end if
         sorted = out%draws(step, :)
         call insertion_sort(sorted)
         out%lower(step) = quantile(sorted, 0.025_dp)
         out%upper(step) = quantile(sorted, 0.975_dp)
      end do
   end function bsts_regression_predict_draws

   function bsts_regression_predict(fit, future_predictors) result(out)
      !! Simulate static-regression predictions using the shared random stream.
      type(bsts_spike_slab_t), intent(in) :: fit !! Previously fitted model.
      real(dp), intent(in) :: future_predictors(:, :) !! Future predictors.
      type(bsts_prediction_t) :: out
      real(dp), allocatable :: normal_draws(:, :)
      integer :: retained, draw, step

      if (.not. allocated(fit%coefficients)) then
         out%info = 1
         return
      end if
      retained = size(fit%coefficients, 2) - fit%burn
      if (size(future_predictors, 1) < 1 .or. retained < 1) then
         out%info = 1
         return
      end if
      allocate(normal_draws(size(future_predictors, 1), retained))
      do draw = 1, retained
         do step = 1, size(future_predictors, 1)
            normal_draws(step, draw) = random_standard_normal()
         end do
      end do
      out = bsts_regression_predict_draws(fit, future_predictors, normal_draws)
   end function bsts_regression_predict

   pure function bsts_structural_prediction_errors(fit, response, &
      initial_mean, initial_covariance, standardize) result(out)
      !! Compute Gaussian filtering errors for retained structural draws.
      type(bsts_mcmc_t), intent(in) :: fit !! Previously fitted model.
      real(dp), intent(in) :: response(:) !! Response observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      logical, intent(in), optional :: standardize !! Flag controlling standardize.
      type(bsts_prediction_errors_t) :: out
      type(ssm_model_t) :: model
      type(kfs_filter_t) :: filtered
      logical :: scaled
      real(dp) :: missing_value
      integer :: dimension, retained, draw, source, time

      dimension = size(initial_mean)
      retained = 0
      if (allocated(fit%observation_variance)) retained = &
         size(fit%observation_variance) - fit%burn
      scaled = .false.
      if (present(standardize)) scaled = standardize
      if (fit%info /= 0 .or. size(response) < 1 .or. dimension < 1 .or. &
         retained < 1 .or. any(shape(initial_covariance) /= &
         [dimension, dimension]) .or. .not. allocated(fit%state_variance) .or. &
         fit%is_semilocal .or. fit%ar_order > 0 .or. &
         allocated(fit%state_weights)) then
         out%info = 1
         return
      end if
      allocate(out%draws(size(response), retained), &
         out%forecast_variance(size(response), retained), &
         out%mean(size(response)))
      missing_value = ieee_value(0.0_dp, ieee_quiet_nan)
      do draw = 1, retained
         source = fit%burn + draw
         model = structural_model_from_fit(response, fit, source, initial_mean, &
            initial_covariance)
         filtered = kfs_filter(model)
         if (filtered%info /= 0) then
            out%info = 10 + filtered%info
            return
         end if
         out%draws(:, draw) = filtered%innovation(:, 1)
         out%forecast_variance(:, draw) = filtered%innovation_cov(1, 1, :)
         do time = 1, size(response)
            if (.not. ieee_is_finite(response(time))) then
               out%draws(time, draw) = missing_value
               out%forecast_variance(time, draw) = missing_value
            else if (scaled) then
               if (out%forecast_variance(time, draw) <= 0.0_dp) then
                  out%info = 2
                  return
               end if
               out%draws(time, draw) = out%draws(time, draw)/ &
                  sqrt(out%forecast_variance(time, draw))
            end if
         end do
      end do
      call prediction_error_summaries(out)
      out%standardized = scaled
      out%training_size = size(response)
   end function bsts_structural_prediction_errors

   pure function bsts_regression_prediction_errors(fit, response, predictors, &
      standardize) result(out)
      !! Compute posterior Gaussian errors for static regression draws.
      type(bsts_spike_slab_t), intent(in) :: fit !! Previously fitted model.
      real(dp), intent(in) :: response(:) !! Response observations.
      real(dp), intent(in) :: predictors(:, :) !! Predictor matrix.
      logical, intent(in), optional :: standardize !! Flag controlling standardize.
      type(bsts_prediction_errors_t) :: out
      logical :: scaled
      real(dp) :: missing_value
      integer :: retained, draw, source, time

      retained = 0
      if (allocated(fit%coefficients)) retained = &
         size(fit%coefficients, 2) - fit%burn
      scaled = .false.
      if (present(standardize)) scaled = standardize
      if (fit%info /= 0 .or. size(response) < 1 .or. retained < 1 .or. &
         size(predictors, 1) /= size(response) .or. &
         size(predictors, 2) /= size(fit%coefficients, 1) .or. &
         .not. all(ieee_is_finite(predictors))) then
         out%info = 1
         return
      end if
      allocate(out%draws(size(response), retained), &
         out%forecast_variance(size(response), retained), &
         out%mean(size(response)))
      missing_value = ieee_value(0.0_dp, ieee_quiet_nan)
      do draw = 1, retained
         source = fit%burn + draw
         out%draws(:, draw) = response - &
            matmul(predictors, fit%coefficients(:, source))
         out%forecast_variance(:, draw) = fit%residual_variance(source)
         do time = 1, size(response)
            if (.not. ieee_is_finite(response(time))) then
               out%draws(time, draw) = missing_value
               out%forecast_variance(time, draw) = missing_value
            else if (scaled) then
               out%draws(time, draw) = out%draws(time, draw)/ &
                  sqrt(fit%residual_variance(source))
            end if
         end do
      end do
      call prediction_error_summaries(out)
      out%standardized = scaled
      out%training_size = size(response)
   end function bsts_regression_prediction_errors

   function bsts_local_level_holdout_errors(response, cutpoint, iterations, &
      burn, standardize) result(out)
      !! Refit a local-level model on a prefix and filter the holdout data.
      real(dp), intent(in) :: response(:) !! Response observations.
      integer, intent(in) :: cutpoint !! Cutpoint.
      integer, intent(in) :: iterations !! Number of algorithm iterations.
      integer, intent(in), optional :: burn !! Number of initial simulation draws to discard.
      logical, intent(in), optional :: standardize !! Flag controlling standardize.
      type(bsts_prediction_errors_t) :: out
      type(bsts_mcmc_t) :: fit
      real(dp) :: initial_mean, initial_variance
      integer :: discarded
      logical :: scaled

      discarded = 0
      if (present(burn)) discarded = burn
      scaled = .false.
      if (present(standardize)) scaled = standardize
      if (cutpoint < 2 .or. cutpoint >= size(response) .or. iterations < 1 .or. &
         discarded < 0 .or. discarded >= iterations) then
         out%info = 1
         return
      end if
      initial_mean = first_finite_value(response(:cutpoint))
      initial_variance = max(finite_sample_variance(response(:cutpoint)), 1.0_dp)
      fit = bsts_local_level(response(:cutpoint), iterations, burn=discarded, &
         initial_mean=initial_mean, initial_variance=initial_variance)
      if (fit%info /= 0) then
         out%info = 10 + fit%info
         return
      end if
      out = bsts_structural_prediction_errors(fit, response, [initial_mean], &
         reshape([initial_variance], [1, 1]), scaled)
      out%training_size = cutpoint
   end function bsts_local_level_holdout_errors

   pure function bsts_compare_prediction_errors(errors, start_index) result(out)
      !! Rank compatible models by posterior mean prediction-error losses.
      type(bsts_prediction_errors_t), intent(in) :: errors(:) !! Errors.
      integer, intent(in), optional :: start_index !! Index of start.
      type(bsts_model_comparison_t) :: out
      logical, allocatable :: finite(:)
      real(dp) :: cumulative
      integer :: models, time_points, first, model, time, destination
      integer :: count_finite

      models = size(errors)
      if (models < 2 .or. .not. allocated(errors(1)%mean)) then
         out%info = 1
         return
      end if
      time_points = size(errors(1)%mean)
      first = 1
      if (present(start_index)) first = start_index
      if (time_points < 1 .or. first < 1 .or. first > time_points) then
         out%info = 1
         return
      end if
      do model = 1, models
         if (errors(model)%info /= 0 .or. .not. allocated(errors(model)%mean) .or. &
            size(errors(model)%mean) /= time_points .or. &
            errors(model)%standardized .neqv. errors(1)%standardized) then
            out%info = 2
            return
         end if
      end do
      allocate(out%cumulative_absolute_error(time_points - first + 1, models), &
         out%rmse(models), out%mae(models), out%final_absolute_error(models), &
         out%rmse_rank(models), out%mae_rank(models))
      do model = 1, models
         finite = ieee_is_finite(errors(model)%mean(first:))
         count_finite = count(finite)
         if (count_finite < 1) then
            out%info = 3
            return
         end if
         out%rmse(model) = sqrt(sum(errors(model)%mean(first:)**2, mask=finite)/ &
            real(count_finite, dp))
         out%mae(model) = sum(abs(errors(model)%mean(first:)), mask=finite)/ &
            real(count_finite, dp)
         cumulative = 0.0_dp
         destination = 0
         do time = first, time_points
            destination = destination + 1
            if (ieee_is_finite(errors(model)%mean(time))) cumulative = &
               cumulative + abs(errors(model)%mean(time))
            out%cumulative_absolute_error(destination, model) = cumulative
         end do
         out%final_absolute_error(model) = cumulative
      end do
      call score_ranks(out%rmse, out%rmse_rank)
      call score_ranks(out%mae, out%mae_rank)
      out%best_rmse = minloc(out%rmse, dim=1)
      out%best_mae = minloc(out%mae, dim=1)
      out%start_index = first
      out%standardized = errors(1)%standardized
   end function bsts_compare_prediction_errors

   pure function bsts_static_intercept_draws(response, prior_mean, &
      prior_variance, observation_variance, observation_prior_shape, &
      observation_prior_rate, burn, intercept_normal_draws, gamma_draws, &
      offset_draws) result(out)
      !! Sample a constant Gaussian state from supplied random variates.
      real(dp), intent(in) :: response(:) !! Response observations.
      real(dp), intent(in) :: prior_mean !! Prior mean.
      real(dp), intent(in) :: prior_variance !! Prior variance.
      real(dp), intent(in) :: observation_variance !! Observation-error variance.
      real(dp), intent(in) :: observation_prior_shape !! Observation prior shape.
      real(dp), intent(in) :: observation_prior_rate !! Observation prior rate.
      integer, intent(in) :: burn !! Number of initial simulation draws to discard.
      real(dp), intent(in) :: intercept_normal_draws(:) !! Intercept normal simulation draws.
      real(dp), intent(in) :: gamma_draws(:) !! Gamma simulation draws.
      real(dp), intent(in), optional :: offset_draws(:, :) !! Offset simulation draws.
      type(bsts_static_intercept_t) :: out
      real(dp), allocatable :: working(:)
      real(dp) :: variance, posterior_variance, posterior_mean, residual_sum
      integer :: observations, observed, iterations, iteration, time

      observations = size(response)
      observed = count(ieee_is_finite(response))
      iterations = size(intercept_normal_draws)
      if (observations < 1 .or. observed < 1 .or. iterations < 1 .or. &
         prior_variance <= 0.0_dp .or. observation_variance <= 0.0_dp .or. &
         observation_prior_shape <= 0.0_dp .or. observation_prior_rate <= 0.0_dp .or. &
         burn < 0 .or. burn >= iterations .or. size(gamma_draws) /= iterations .or. &
         any(gamma_draws <= 0.0_dp) .or. &
         .not. all(ieee_is_finite(intercept_normal_draws))) then
         out%info = 1
         return
      end if
      if (present(offset_draws)) then
         if (any(shape(offset_draws) /= [observations, iterations]) .or. &
            .not. all(ieee_is_finite(offset_draws))) then
            out%info = 1
            return
         end if
      end if
      allocate(working(observations), out%intercept(iterations), &
         out%observation_variance(iterations), &
         out%contribution(observations, iterations), &
         out%fitted(observations, iterations), &
         out%residuals(observations, iterations))
      variance = observation_variance
      out%burn = burn
      do iteration = 1, iterations
         working = response
         if (present(offset_draws)) working = working - offset_draws(:, iteration)
         posterior_variance = 1.0_dp/(1.0_dp/prior_variance + &
            real(observed, dp)/variance)
         posterior_mean = posterior_variance*(prior_mean/prior_variance + &
            sum(working, mask=ieee_is_finite(response))/variance)
         out%intercept(iteration) = posterior_mean + sqrt(posterior_variance)* &
            intercept_normal_draws(iteration)
         residual_sum = 0.0_dp
         do time = 1, observations
            if (ieee_is_finite(response(time))) residual_sum = residual_sum + &
               (working(time) - out%intercept(iteration))**2
         end do
         variance = (observation_prior_rate + 0.5_dp*residual_sum)/ &
            gamma_draws(iteration)
         out%observation_variance(iteration) = variance
         out%contribution(:, iteration) = out%intercept(iteration)
         out%fitted(:, iteration) = out%intercept(iteration)
         if (present(offset_draws)) out%fitted(:, iteration) = &
            out%fitted(:, iteration) + offset_draws(:, iteration)
         out%residuals(:, iteration) = response - out%fitted(:, iteration)
      end do
   end function bsts_static_intercept_draws

   function bsts_static_intercept(response, iterations, burn, prior_mean, &
      prior_variance, observation_variance, observation_prior_shape, &
      observation_prior_rate, offset_draws) result(out)
      !! Sample a constant Gaussian state using the shared random stream.
      real(dp), intent(in) :: response(:) !! Response observations.
      integer, intent(in) :: iterations !! Number of algorithm iterations.
      integer, intent(in), optional :: burn !! Number of initial simulation draws to discard.
      real(dp), intent(in), optional :: prior_mean !! Prior mean.
      real(dp), intent(in), optional :: prior_variance !! Prior variance.
      real(dp), intent(in), optional :: observation_variance !! Observation-error variance.
      real(dp), intent(in), optional :: observation_prior_shape !! Observation prior shape.
      real(dp), intent(in), optional :: observation_prior_rate !! Observation prior rate.
      real(dp), intent(in), optional :: offset_draws(:, :) !! Offset simulation draws.
      type(bsts_static_intercept_t) :: out
      real(dp), allocatable :: normals(:), gammas(:)
      real(dp) :: mean0, prior0, observation0, shape0, rate0, series_variance
      integer :: discarded, observed, iteration

      if (size(response) < 1 .or. iterations < 1 .or. &
         count(ieee_is_finite(response)) < 1) then
         out%info = 1
         return
      end if
      discarded = 0
      if (present(burn)) discarded = burn
      series_variance = finite_sample_variance(response)
      mean0 = first_finite_value(response)
      if (present(prior_mean)) mean0 = prior_mean
      prior0 = max(series_variance, 1.0e-6_dp)
      if (present(prior_variance)) prior0 = prior_variance
      observation0 = max(series_variance, 1.0e-6_dp)
      if (present(observation_variance)) observation0 = observation_variance
      shape0 = 0.01_dp
      if (present(observation_prior_shape)) shape0 = observation_prior_shape
      rate0 = max(0.01_dp*series_variance, 1.0e-8_dp)
      if (present(observation_prior_rate)) rate0 = observation_prior_rate
      if (discarded < 0 .or. discarded >= iterations .or. prior0 <= 0.0_dp .or. &
         observation0 <= 0.0_dp .or. shape0 <= 0.0_dp .or. rate0 <= 0.0_dp) then
         out%info = 1
         return
      end if
      observed = count(ieee_is_finite(response))
      allocate(normals(iterations), gammas(iterations))
      do iteration = 1, iterations
         normals(iteration) = random_standard_normal()
         gammas(iteration) = random_gamma(shape0 + 0.5_dp*real(observed, dp))
      end do
      if (present(offset_draws)) then
         out = bsts_static_intercept_draws(response, mean0, prior0, &
            observation0, shape0, rate0, discarded, normals, gammas, offset_draws)
      else
         out = bsts_static_intercept_draws(response, mean0, prior0, &
            observation0, shape0, rate0, discarded, normals, gammas)
      end if
   end function bsts_static_intercept

   pure function bsts_static_intercept_predict_draws(fit, future_offset, &
      observation_normal_draws) result(out)
      !! Forecast a constant state from supplied observation-normal draws.
      type(bsts_static_intercept_t), intent(in) :: fit !! Previously fitted model.
      real(dp), intent(in) :: future_offset(:) !! Future offset.
      real(dp), intent(in) :: observation_normal_draws(:, :) !! Observation normal draws.
      type(bsts_prediction_t) :: out
      integer :: horizon, retained, draw, source

      horizon = size(future_offset)
      retained = 0
      if (allocated(fit%intercept)) retained = size(fit%intercept) - fit%burn
      if (fit%info /= 0 .or. horizon < 1 .or. retained < 1 .or. &
         any(shape(observation_normal_draws) /= [horizon, retained]) .or. &
         .not. all(ieee_is_finite(future_offset)) .or. &
         .not. all(ieee_is_finite(observation_normal_draws))) then
         out%info = 1
         return
      end if
      allocate(out%draws(horizon, retained))
      do draw = 1, retained
         source = fit%burn + draw
         out%draws(:, draw) = fit%intercept(source) + future_offset + &
            sqrt(fit%observation_variance(source))* &
            observation_normal_draws(:, draw)
      end do
      call summarize_prediction(out)
   end function bsts_static_intercept_predict_draws

   function bsts_static_intercept_predict(fit, future_offset) result(out)
      !! Forecast a constant state using the shared random stream.
      type(bsts_static_intercept_t), intent(in) :: fit !! Previously fitted model.
      real(dp), intent(in) :: future_offset(:) !! Future offset.
      type(bsts_prediction_t) :: out
      real(dp), allocatable :: normals(:, :)
      integer :: retained, draw, time

      if (.not. allocated(fit%intercept)) then
         out%info = 1
         return
      end if
      retained = size(fit%intercept) - fit%burn
      if (size(future_offset) < 1 .or. retained < 1) then
         out%info = 1
         return
      end if
      allocate(normals(size(future_offset), retained))
      do draw = 1, retained
         do time = 1, size(future_offset)
            normals(time, draw) = random_standard_normal()
         end do
      end do
      out = bsts_static_intercept_predict_draws(fit, future_offset, normals)
   end function bsts_static_intercept_predict

   pure function bsts_dirm_draws(response, predictors, time_index, &
      initial_level_mean, initial_level_variance, observation_variance, &
      level_variance, slab_mean, slab_covariance, &
      prior_inclusion_probability, observation_prior_shape, &
      observation_prior_rate, level_prior_shape, level_prior_rate, &
      initial_inclusion, maximum_model_size, maximum_flips, burn, &
      state_normal_draws, coefficient_normal_draws, inclusion_uniform_draws, &
      gamma_draws) result(out)
      !! Sample dynamic-intercept regression from supplied random variates.
      real(dp), intent(in) :: response(:) !! Response observations.
      real(dp), intent(in) :: predictors(:, :) !! Predictor matrix.
      integer, intent(in) :: time_index(:) !! Index of time.
      real(dp), intent(in) :: initial_level_mean !! Initial level mean.
      real(dp), intent(in) :: initial_level_variance !! Initial level variance.
      real(dp), intent(in) :: observation_variance !! Observation-error variance.
      real(dp), intent(in) :: level_variance !! Level variance.
      real(dp), intent(in) :: slab_mean(:) !! Slab mean.
      real(dp), intent(in) :: slab_covariance(:, :) !! Slab covariance matrix.
      real(dp), intent(in) :: prior_inclusion_probability(:) !! Prior inclusion probability.
      real(dp), intent(in) :: observation_prior_shape !! Observation prior shape.
      real(dp), intent(in) :: observation_prior_rate !! Observation prior rate.
      real(dp), intent(in) :: level_prior_shape !! Level prior shape.
      real(dp), intent(in) :: level_prior_rate !! Level prior rate.
      logical, intent(in) :: initial_inclusion(:) !! Initial inclusion.
      integer, intent(in) :: maximum_model_size !! Maximum model size.
      integer, intent(in) :: maximum_flips !! Maximum flips.
      integer, intent(in) :: burn !! Number of initial simulation draws to discard.
      real(dp), intent(in) :: state_normal_draws(:, :, :) !! Independent standard-normal state draws.
      real(dp), intent(in) :: coefficient_normal_draws(:, :) !! Coefficient normal simulation draws.
      real(dp), intent(in) :: inclusion_uniform_draws(:, :) !! Inclusion uniform simulation draws.
      real(dp), intent(in) :: gamma_draws(:, :) !! Gamma simulation draws.
      type(bsts_dirm_t) :: out
      type(ssm_model_t) :: model
      real(dp), allocatable :: observed_y(:), observed_x(:, :), working_y(:)
      real(dp), allocatable :: aggregated(:), posterior_mean(:)
      real(dp), allocatable :: posterior_covariance(:, :), active_draw(:)
      real(dp), allocatable :: coefficient(:), path(:, :)
      integer, allocatable :: observed_index(:), counts(:), active(:)
      logical, allocatable :: included(:), candidate(:)
      real(dp) :: observation_scale, level_scale, log_included, log_excluded
      real(dp) :: log_probability, collapsed_rate, probability, level_sum
      integer :: observations, observed, variables, time_points, iterations
      integer :: iteration, variable, row, time, flips, status, retained

      observations = size(response)
      variables = size(predictors, 2)
      iterations = size(coefficient_normal_draws, 2)
      if (.not. dirm_time_index_valid(time_index) .or. observations < 1 .or. &
         size(predictors, 1) /= observations .or. &
         size(time_index) /= observations .or. variables < 1 .or. &
         count(ieee_is_finite(response)) < 1 .or. &
         .not. all(ieee_is_finite(predictors)) .or. &
         size(slab_mean) /= variables .or. &
         any(shape(slab_covariance) /= [variables, variables]) .or. &
         size(prior_inclusion_probability) /= variables .or. &
         any(prior_inclusion_probability < 0.0_dp) .or. &
         any(prior_inclusion_probability > 1.0_dp) .or. &
         size(initial_inclusion) /= variables .or. &
         .not. ieee_is_finite(initial_level_mean) .or. &
         .not. all(ieee_is_finite(slab_mean)) .or. &
         .not. all(ieee_is_finite(slab_covariance)) .or. &
         initial_level_variance <= 0.0_dp .or. observation_variance <= 0.0_dp .or. &
         level_variance <= 0.0_dp .or. observation_prior_shape <= 0.0_dp .or. &
         observation_prior_rate <= 0.0_dp .or. level_prior_shape <= 0.0_dp .or. &
         level_prior_rate <= 0.0_dp .or. maximum_model_size < 0 .or. &
         maximum_model_size > variables .or. maximum_flips < -1 .or. &
         iterations < 1 .or. burn < 0 .or. burn >= iterations .or. &
         any(shape(coefficient_normal_draws) /= [variables, iterations]) .or. &
         any(shape(inclusion_uniform_draws) /= [variables, iterations]) .or. &
         any(shape(gamma_draws) /= [2, iterations]) .or. &
         size(state_normal_draws, 1) /= 1 .or. &
         size(state_normal_draws, 3) /= iterations .or. &
         any(inclusion_uniform_draws <= 0.0_dp) .or. &
         any(inclusion_uniform_draws >= 1.0_dp) .or. any(gamma_draws <= 0.0_dp) .or. &
         .not. all(ieee_is_finite(state_normal_draws)) .or. &
         .not. all(ieee_is_finite(coefficient_normal_draws))) then
         out%info = 1
         return
      end if
      time_points = maxval(time_index)
      if (size(state_normal_draws, 2) /= time_points) then
         out%info = 1
         return
      end if
      observed = count(ieee_is_finite(response))
      observed_index = pack([(row, row=1, observations)], &
         ieee_is_finite(response))
      allocate(observed_y(observed), observed_x(observed, variables), &
         working_y(observed), aggregated(time_points), counts(time_points), &
         coefficient(variables), included(variables), candidate(variables), &
         out%state(time_points, iterations), &
         out%coefficients(variables, iterations), &
         out%included(variables, iterations), &
         out%observation_variance(iterations), out%level_variance(iterations), &
         out%contribution(observations, iterations), &
         out%fitted(observations, iterations), &
         out%residuals(observations, iterations), &
         out%inclusion_probability(variables))
      observed_y = response(observed_index)
      observed_x = predictors(observed_index, :)
      coefficient = 0.0_dp
      included = initial_inclusion
      where (prior_inclusion_probability <= 0.0_dp) included = .false.
      where (prior_inclusion_probability >= 1.0_dp) included = .true.
      if (count(included) > maximum_model_size) then
         out%info = 2
         return
      end if
      observation_scale = observation_variance
      level_scale = level_variance
      out%time_points = time_points
      out%burn = burn
      out%state = initial_level_mean
      do iteration = 1, iterations
         do row = 1, observed
            working_y(row) = observed_y(row) - &
               out%state(time_index(observed_index(row)), max(1, iteration - 1))
         end do
         flips = 0
         do variable = 1, variables
            if (prior_inclusion_probability(variable) <= 0.0_dp) then
               included(variable) = .false.
               cycle
            end if
            if (prior_inclusion_probability(variable) >= 1.0_dp) then
               included(variable) = .true.
               cycle
            end if
            candidate = included
            candidate(variable) = .false.
            call spike_slab_model_moments(working_y, observed_x, candidate, &
               slab_mean, slab_covariance, prior_inclusion_probability, &
               observation_scale, observation_prior_rate, log_excluded, &
               collapsed_rate, posterior_mean, posterior_covariance, active, status)
            if (status /= 0) then
               out%info = 10 + status
               return
            end if
            candidate(variable) = .true.
            if (count(candidate) > maximum_model_size) then
               probability = 0.0_dp
            else
               call spike_slab_model_moments(working_y, observed_x, candidate, &
                  slab_mean, slab_covariance, prior_inclusion_probability, &
                  observation_scale, observation_prior_rate, log_included, &
                  collapsed_rate, posterior_mean, posterior_covariance, active, &
                  status)
               if (status /= 0) then
                  out%info = 20 + status
                  return
               end if
               probability = logistic_log_odds(log_included - log_excluded)
            end if
            candidate(variable) = inclusion_uniform_draws(variable, iteration) < &
               probability
            if (candidate(variable) .neqv. included(variable)) then
               if (maximum_flips < 0 .or. flips < maximum_flips) then
                  included(variable) = candidate(variable)
                  flips = flips + 1
               end if
            end if
         end do
         call spike_slab_model_moments(working_y, observed_x, included, &
            slab_mean, slab_covariance, prior_inclusion_probability, &
            observation_scale, observation_prior_rate, log_probability, &
            collapsed_rate, posterior_mean, posterior_covariance, active, status)
         if (status /= 0) then
            out%info = 30 + status
            return
         end if
         observation_scale = collapsed_rate/gamma_draws(1, iteration)
         coefficient = 0.0_dp
         if (size(active) > 0) then
            allocate(active_draw(size(active)))
            call multivariate_normal_from_standard(posterior_mean, &
               observation_scale*posterior_covariance, &
               coefficient_normal_draws(active, iteration), active_draw, status)
            if (status /= 0) then
               out%info = 40 + status
               return
            end if
            coefficient(active) = active_draw
            deallocate(active_draw)
         end if
         aggregated = 0.0_dp
         counts = 0
         do row = 1, observed
            time = time_index(observed_index(row))
            aggregated(time) = aggregated(time) + observed_y(row) - &
               dot_product(observed_x(row, :), coefficient)
            counts(time) = counts(time) + 1
         end do
         do time = 1, time_points
            if (counts(time) > 0) aggregated(time) = &
               aggregated(time)/real(counts(time), dp)
         end do
         model = dirm_local_level_model(aggregated, counts, initial_level_mean, &
            initial_level_variance, observation_scale, level_scale)
         call state_path_from_draws(model, state_normal_draws(:, :, iteration), &
            path, status)
         if (status /= 0) then
            out%info = 50 + status
            return
         end if
         level_sum = sum((path(1, 2:) - path(1, :time_points - 1))**2)
         level_scale = (level_prior_rate + 0.5_dp*level_sum)/ &
            gamma_draws(2, iteration)
         out%state(:, iteration) = path(1, :)
         out%coefficients(:, iteration) = coefficient
         out%included(:, iteration) = included
         out%observation_variance(iteration) = observation_scale
         out%level_variance(iteration) = level_scale
         do row = 1, observations
            out%contribution(row, iteration) = path(1, time_index(row))
            out%fitted(row, iteration) = out%contribution(row, iteration) + &
               dot_product(predictors(row, :), coefficient)
            out%residuals(row, iteration) = response(row) - &
               out%fitted(row, iteration)
         end do
      end do
      retained = iterations - burn
      do variable = 1, variables
         out%inclusion_probability(variable) = real(count( &
            out%included(variable, burn + 1:)), dp)/real(retained, dp)
      end do
   end function bsts_dirm_draws

   function bsts_dirm(response, predictors, time_index, iterations, burn, &
      initial_level_mean, initial_level_variance, observation_variance, &
      level_variance, slab_mean, slab_covariance, &
      prior_inclusion_probability, expected_model_size, &
      observation_prior_shape, observation_prior_rate, level_prior_shape, &
      level_prior_rate, initial_inclusion, maximum_model_size, maximum_flips) &
      result(out)
      !! Sample dynamic-intercept regression using the shared random stream.
      real(dp), intent(in) :: response(:) !! Response observations.
      real(dp), intent(in) :: predictors(:, :) !! Predictor matrix.
      integer, intent(in) :: time_index(:) !! Index of time.
      integer, intent(in) :: iterations !! Number of algorithm iterations.
      integer, intent(in), optional :: burn !! Number of initial simulation draws to discard.
      integer, intent(in), optional :: maximum_model_size !! Maximum model size.
      integer, intent(in), optional :: maximum_flips !! Maximum flips.
      real(dp), intent(in), optional :: initial_level_mean !! Initial level mean.
      real(dp), intent(in), optional :: initial_level_variance !! Initial level variance.
      real(dp), intent(in), optional :: observation_variance !! Observation-error variance.
      real(dp), intent(in), optional :: level_variance !! Level variance.
      real(dp), intent(in), optional :: slab_mean(:) !! Slab mean.
      real(dp), intent(in), optional :: slab_covariance(:, :) !! Slab covariance matrix.
      real(dp), intent(in), optional :: prior_inclusion_probability(:) !! Prior inclusion probability.
      real(dp), intent(in), optional :: expected_model_size !! Expected model size.
      real(dp), intent(in), optional :: observation_prior_shape !! Observation prior shape.
      real(dp), intent(in), optional :: observation_prior_rate !! Observation prior rate.
      real(dp), intent(in), optional :: level_prior_shape !! Level prior shape.
      real(dp), intent(in), optional :: level_prior_rate !! Level prior rate.
      logical, intent(in), optional :: initial_inclusion(:) !! Initial inclusion.
      type(bsts_dirm_t) :: out
      real(dp), allocatable :: mean0(:), covariance0(:, :), inclusion0(:)
      real(dp), allocatable :: state_normals(:, :, :), coefficient_normals(:, :)
      real(dp), allocatable :: uniforms(:, :), gammas(:, :), xtx(:, :), inverse(:, :)
      logical, allocatable :: model0(:)
      real(dp) :: level_mean0, level_initial0, observation0, level0
      real(dp) :: observation_shape0, observation_rate0, level_shape0, level_rate0
      real(dp) :: expected, series_variance
      integer :: variables, time_points, discarded, model_limit, flip_limit
      integer :: iteration, variable, observed, status

      variables = size(predictors, 2)
      if (iterations < 1 .or. variables < 1 .or. &
         size(response) /= size(time_index) .or. &
         size(predictors, 1) /= size(response) .or. &
         .not. dirm_time_index_valid(time_index) .or. &
         count(ieee_is_finite(response)) < 1 .or. &
         .not. all(ieee_is_finite(predictors))) then
         out%info = 1
         return
      end if
      time_points = maxval(time_index)
      observed = count(ieee_is_finite(response))
      discarded = 0
      if (present(burn)) discarded = burn
      model_limit = variables
      if (present(maximum_model_size)) model_limit = maximum_model_size
      flip_limit = -1
      if (present(maximum_flips)) flip_limit = maximum_flips
      series_variance = finite_sample_variance(response)
      level_mean0 = first_finite_value(response)
      if (present(initial_level_mean)) level_mean0 = initial_level_mean
      level_initial0 = max(series_variance, 1.0_dp)
      if (present(initial_level_variance)) level_initial0 = initial_level_variance
      observation0 = max(0.5_dp*series_variance, 1.0e-6_dp)
      if (present(observation_variance)) observation0 = observation_variance
      level0 = max(0.01_dp*series_variance, 1.0e-8_dp)
      if (present(level_variance)) level0 = level_variance
      observation_shape0 = 0.01_dp
      if (present(observation_prior_shape)) observation_shape0 = &
         observation_prior_shape
      observation_rate0 = max(0.01_dp*series_variance, 1.0e-8_dp)
      if (present(observation_prior_rate)) observation_rate0 = &
         observation_prior_rate
      level_shape0 = 0.01_dp
      if (present(level_prior_shape)) level_shape0 = level_prior_shape
      level_rate0 = max(0.0001_dp*series_variance, 1.0e-10_dp)
      if (present(level_prior_rate)) level_rate0 = level_prior_rate
      allocate(mean0(variables), covariance0(variables, variables), &
         inclusion0(variables), model0(variables))
      mean0 = 0.0_dp
      if (present(slab_mean)) then
         if (size(slab_mean) /= variables) then
            out%info = 1
            return
         end if
         mean0 = slab_mean
      end if
      xtx = matmul(transpose(pack_predictors(response, predictors)), &
         pack_predictors(response, predictors))
      call invert_matrix(xtx, inverse, status)
      if (status == 0) then
         covariance0 = real(max(1, observed), dp)*inverse
      else
         covariance0 = 0.0_dp
         do variable = 1, variables
            covariance0(variable, variable) = real(max(1, observed), dp)/ &
               max(xtx(variable, variable), 1.0_dp)
         end do
      end if
      if (present(slab_covariance)) then
         if (any(shape(slab_covariance) /= [variables, variables])) then
            out%info = 1
            return
         end if
         covariance0 = slab_covariance
      end if
      expected = real(min(3, variables), dp)
      if (present(expected_model_size)) expected = expected_model_size
      inclusion0 = max(0.0_dp, min(1.0_dp, expected/real(variables, dp)))
      if (present(prior_inclusion_probability)) then
         if (size(prior_inclusion_probability) /= variables) then
            out%info = 1
            return
         end if
         inclusion0 = prior_inclusion_probability
      end if
      do variable = 1, variables
         if (sum((predictors(:, variable) - sum(predictors(:, variable))/ &
            real(size(predictors, 1), dp))**2) <= 100.0_dp*epsilon(1.0_dp)) &
            inclusion0(variable) = 0.0_dp
      end do
      model0 = inclusion0 >= 1.0_dp
      if (present(initial_inclusion)) then
         if (size(initial_inclusion) /= variables) then
            out%info = 1
            return
         end if
         model0 = initial_inclusion
      end if
      if (discarded < 0 .or. discarded >= iterations .or. expected < 0.0_dp .or. &
         model_limit < 0 .or. model_limit > variables .or. flip_limit < -1 .or. &
         level_initial0 <= 0.0_dp .or. observation0 <= 0.0_dp .or. &
         level0 <= 0.0_dp .or. observation_shape0 <= 0.0_dp .or. &
         observation_rate0 <= 0.0_dp .or. level_shape0 <= 0.0_dp .or. &
         level_rate0 <= 0.0_dp) then
         out%info = 1
         return
      end if
      allocate(state_normals(1, time_points, iterations), &
         coefficient_normals(variables, iterations), &
         uniforms(variables, iterations), gammas(2, iterations))
      do iteration = 1, iterations
         do variable = 1, variables
            coefficient_normals(variable, iteration) = random_standard_normal()
            uniforms(variable, iteration) = max(tiny(1.0_dp), &
               min(1.0_dp - epsilon(1.0_dp), random_uniform()))
         end do
         do variable = 1, time_points
            state_normals(1, variable, iteration) = random_standard_normal()
         end do
         gammas(1, iteration) = random_gamma(observation_shape0 + &
            0.5_dp*real(observed, dp))
         gammas(2, iteration) = random_gamma(level_shape0 + &
            0.5_dp*real(max(0, time_points - 1), dp))
      end do
      out = bsts_dirm_draws(response, predictors, time_index, level_mean0, &
         level_initial0, observation0, level0, mean0, covariance0, inclusion0, &
         observation_shape0, observation_rate0, level_shape0, level_rate0, &
         model0, model_limit, flip_limit, discarded, state_normals, &
         coefficient_normals, uniforms, gammas)
   end function bsts_dirm

   pure function bsts_dirm_predict_draws(fit, future_predictors, &
      future_time_index, state_normal_draws, observation_normal_draws) result(out)
      !! Forecast grouped dynamic-intercept observations from supplied draws.
      type(bsts_dirm_t), intent(in) :: fit !! Previously fitted model.
      real(dp), intent(in) :: future_predictors(:, :) !! Future predictors.
      integer, intent(in) :: future_time_index(:) !! Future time index.
      real(dp), intent(in) :: state_normal_draws(:, :) !! Independent standard-normal state draws.
      real(dp), intent(in) :: observation_normal_draws(:, :) !! Observation normal draws.
      type(bsts_prediction_t) :: out
      real(dp), allocatable :: sorted(:)
      real(dp) :: level
      integer :: observations, time_points, retained, draw, source, row, time

      observations = size(future_predictors, 1)
      retained = 0
      if (allocated(fit%coefficients)) retained = &
         size(fit%coefficients, 2) - fit%burn
      if (fit%info /= 0 .or. observations < 1 .or. retained < 1 .or. &
         .not. dirm_time_index_valid(future_time_index) .or. &
         size(future_time_index) /= observations .or. &
         size(future_predictors, 2) /= size(fit%coefficients, 1) .or. &
         .not. all(ieee_is_finite(future_predictors))) then
         out%info = 1
         return
      end if
      time_points = maxval(future_time_index)
      if (any(shape(state_normal_draws) /= [time_points, retained]) .or. &
         any(shape(observation_normal_draws) /= [observations, retained]) .or. &
         .not. all(ieee_is_finite(state_normal_draws)) .or. &
         .not. all(ieee_is_finite(observation_normal_draws))) then
         out%info = 1
         return
      end if
      allocate(out%draws(observations, retained), out%mean(observations), &
         out%standard_deviation(observations), out%lower(observations), &
         out%upper(observations))
      do draw = 1, retained
         source = fit%burn + draw
         level = fit%state(fit%time_points, source)
         do time = 1, time_points
            level = level + sqrt(fit%level_variance(source))* &
               state_normal_draws(time, draw)
            do row = 1, observations
               if (future_time_index(row) /= time) cycle
               out%draws(row, draw) = level + &
                  dot_product(future_predictors(row, :), &
                  fit%coefficients(:, source)) + &
                  sqrt(fit%observation_variance(source))* &
                  observation_normal_draws(row, draw)
            end do
         end do
      end do
      do row = 1, observations
         out%mean(row) = sum(out%draws(row, :))/real(retained, dp)
         if (retained > 1) then
            out%standard_deviation(row) = sqrt(sum((out%draws(row, :) - &
               out%mean(row))**2)/real(retained - 1, dp))
         else
            out%standard_deviation(row) = 0.0_dp
         end if
         sorted = out%draws(row, :)
         call insertion_sort(sorted)
         out%lower(row) = quantile(sorted, 0.025_dp)
         out%upper(row) = quantile(sorted, 0.975_dp)
      end do
   end function bsts_dirm_predict_draws

   function bsts_dirm_predict(fit, future_predictors, future_time_index) &
      result(out)
      !! Forecast grouped dynamic-intercept observations using random draws.
      type(bsts_dirm_t), intent(in) :: fit !! Previously fitted model.
      real(dp), intent(in) :: future_predictors(:, :) !! Future predictors.
      integer, intent(in) :: future_time_index(:) !! Future time index.
      type(bsts_prediction_t) :: out
      real(dp), allocatable :: state_normals(:, :), observation_normals(:, :)
      integer :: retained, time_points, draw, row

      if (.not. allocated(fit%coefficients) .or. &
         .not. dirm_time_index_valid(future_time_index)) then
         out%info = 1
         return
      end if
      retained = size(fit%coefficients, 2) - fit%burn
      time_points = maxval(future_time_index)
      if (retained < 1 .or. size(future_predictors, 1) /= &
         size(future_time_index)) then
         out%info = 1
         return
      end if
      allocate(state_normals(time_points, retained), &
         observation_normals(size(future_time_index), retained))
      do draw = 1, retained
         do row = 1, time_points
            state_normals(row, draw) = random_standard_normal()
         end do
         do row = 1, size(future_time_index)
            observation_normals(row, draw) = random_standard_normal()
         end do
      end do
      out = bsts_dirm_predict_draws(fit, future_predictors, future_time_index, &
         state_normals, observation_normals)
   end function bsts_dirm_predict

   pure function bsts_mixed_draws(coarse_response, predictors, coarse_index, &
      membership_fraction, contains_end, initial_level_mean, &
      initial_level_variance, observation_variance, level_variance, slab_mean, &
      slab_covariance, prior_inclusion_probability, observation_prior_shape, &
      observation_prior_rate, level_prior_shape, level_prior_rate, &
      initial_inclusion, maximum_model_size, maximum_flips, burn, &
      latent_normal_draws, state_normal_draws, coefficient_normal_draws, &
      inclusion_uniform_draws, gamma_draws) result(out)
      !! Sample a local-level mixed-frequency model from supplied variates.
      real(dp), intent(in) :: coarse_response(:) !! Coarse response.
      real(dp), intent(in) :: predictors(:, :) !! Predictor matrix.
      integer, intent(in) :: coarse_index(:) !! Index of coarse.
      real(dp), intent(in) :: membership_fraction(:) !! Membership fraction.
      logical, intent(in) :: contains_end(:) !! Flag controlling contains end.
      real(dp), intent(in) :: initial_level_mean !! Initial level mean.
      real(dp), intent(in) :: initial_level_variance !! Initial level variance.
      real(dp), intent(in) :: observation_variance !! Observation-error variance.
      real(dp), intent(in) :: level_variance !! Level variance.
      real(dp), intent(in) :: slab_mean(:) !! Slab mean.
      real(dp), intent(in) :: slab_covariance(:, :) !! Slab covariance matrix.
      real(dp), intent(in) :: prior_inclusion_probability(:) !! Prior inclusion probability.
      real(dp), intent(in) :: observation_prior_shape !! Observation prior shape.
      real(dp), intent(in) :: observation_prior_rate !! Observation prior rate.
      real(dp), intent(in) :: level_prior_shape !! Level prior shape.
      real(dp), intent(in) :: level_prior_rate !! Level prior rate.
      logical, intent(in) :: initial_inclusion(:) !! Initial inclusion.
      integer, intent(in) :: maximum_model_size !! Maximum model size.
      integer, intent(in) :: maximum_flips !! Maximum flips.
      integer, intent(in) :: burn !! Number of initial simulation draws to discard.
      real(dp), intent(in) :: latent_normal_draws(:, :) !! Latent normal simulation draws.
      real(dp), intent(in) :: state_normal_draws(:, :, :) !! Independent standard-normal state draws.
      real(dp), intent(in) :: coefficient_normal_draws(:, :) !! Coefficient normal simulation draws.
      real(dp), intent(in) :: inclusion_uniform_draws(:, :) !! Inclusion uniform simulation draws.
      real(dp), intent(in) :: gamma_draws(:, :) !! Gamma simulation draws.
      type(bsts_mixed_t) :: out
      type(bsts_spike_slab_t) :: regression
      type(ssm_model_t) :: model
      real(dp), allocatable :: design(:, :), observed_design(:, :)
      real(dp), allocatable :: constraint_covariance(:, :), constraint_inverse(:, :)
      real(dp), allocatable :: target(:), latent(:), raw_latent(:), state(:)
      real(dp), allocatable :: coefficient(:), state_response(:), path(:, :)
      real(dp), allocatable :: offset(:, :)
      integer, allocatable :: observed_coarse(:)
      logical, allocatable :: included(:)
      real(dp) :: observation_scale, level_scale, level_sum
      integer :: fine_points, coarse_points, variables, iterations
      integer :: iteration, row, status, retained, variable

      fine_points = size(predictors, 1)
      coarse_points = size(coarse_response)
      variables = size(predictors, 2)
      iterations = size(latent_normal_draws, 2)
      if (fine_points < 2 .or. coarse_points < 1 .or. variables < 1 .or. &
         size(coarse_index) /= fine_points .or. &
         size(membership_fraction) /= fine_points .or. &
         size(contains_end) /= fine_points .or. &
         any(membership_fraction < 0.0_dp) .or. &
         any(membership_fraction > 1.0_dp) .or. &
         .not. all(ieee_is_finite(membership_fraction)) .or. &
         .not. all(ieee_is_finite(predictors)) .or. &
         size(slab_mean) /= variables .or. &
         any(shape(slab_covariance) /= [variables, variables]) .or. &
         size(prior_inclusion_probability) /= variables .or. &
         size(initial_inclusion) /= variables .or. &
         any(prior_inclusion_probability < 0.0_dp) .or. &
         any(prior_inclusion_probability > 1.0_dp) .or. &
         .not. ieee_is_finite(initial_level_mean) .or. &
         .not. all(ieee_is_finite(slab_mean)) .or. &
         .not. all(ieee_is_finite(slab_covariance)) .or. &
         initial_level_variance <= 0.0_dp .or. observation_variance <= 0.0_dp .or. &
         level_variance <= 0.0_dp .or. observation_prior_shape <= 0.0_dp .or. &
         observation_prior_rate <= 0.0_dp .or. level_prior_shape <= 0.0_dp .or. &
         level_prior_rate <= 0.0_dp .or. maximum_model_size < 0 .or. &
         maximum_model_size > variables .or. maximum_flips < -1 .or. &
         iterations < 1 .or. burn < 0 .or. burn >= iterations .or. &
         size(latent_normal_draws, 1) /= fine_points .or. &
         any(shape(state_normal_draws) /= [1, fine_points, iterations]) .or. &
         any(shape(coefficient_normal_draws) /= [variables, iterations]) .or. &
         any(shape(inclusion_uniform_draws) /= [variables, iterations]) .or. &
         any(shape(gamma_draws) /= [2, iterations]) .or. &
         any(inclusion_uniform_draws <= 0.0_dp) .or. &
         any(inclusion_uniform_draws >= 1.0_dp) .or. any(gamma_draws <= 0.0_dp) .or. &
         .not. all(ieee_is_finite(latent_normal_draws)) .or. &
         .not. all(ieee_is_finite(state_normal_draws)) .or. &
         .not. all(ieee_is_finite(coefficient_normal_draws))) then
         out%info = 1
         return
      end if
      design = mixed_frequency_design(coarse_index, membership_fraction, &
         contains_end, coarse_points)
      do row = 1, coarse_points
         if (ieee_is_finite(coarse_response(row)) .and. &
            sum(abs(design(row, :))) <= tiny(1.0_dp)) then
            out%info = 2
            return
         end if
      end do
      observed_coarse = pack([(row, row=1, coarse_points)], &
         ieee_is_finite(coarse_response))
      if (size(observed_coarse) > 0) then
         observed_design = design(observed_coarse, :)
         target = coarse_response(observed_coarse)
         constraint_covariance = matmul(observed_design, &
            transpose(observed_design))
         call invert_matrix(constraint_covariance, constraint_inverse, status)
         if (status /= 0) then
            out%info = 3
            return
         end if
      else
         allocate(observed_design(0, fine_points), target(0), &
            constraint_inverse(0, 0))
      end if
      allocate(latent(fine_points), raw_latent(fine_points), state(fine_points), &
         coefficient(variables), included(variables), state_response(fine_points), &
         offset(fine_points, 1), out%latent_fine(fine_points, iterations), &
         out%cumulator(fine_points, iterations), &
         out%state(fine_points, iterations), &
         out%coefficients(variables, iterations), &
         out%included(variables, iterations), &
         out%observation_variance(iterations), out%level_variance(iterations), &
         out%state_contribution(fine_points, iterations), &
         out%regression_contribution(fine_points, iterations), &
         out%coarse_fitted(coarse_points, iterations), &
         out%inclusion_probability(variables))
      state = initial_level_mean
      coefficient = 0.0_dp
      included = initial_inclusion
      where (prior_inclusion_probability <= 0.0_dp) included = .false.
      where (prior_inclusion_probability >= 1.0_dp) included = .true.
      if (count(included) > maximum_model_size) then
         out%info = 4
         return
      end if
      observation_scale = observation_variance
      level_scale = level_variance
      out%burn = burn
      do iteration = 1, iterations
         raw_latent = state + matmul(predictors, coefficient) + &
            sqrt(observation_scale)*latent_normal_draws(:, iteration)
         latent = raw_latent
         if (size(observed_coarse) > 0) latent = raw_latent + &
            matmul(transpose(observed_design), matmul(constraint_inverse, &
            target - matmul(observed_design, raw_latent)))
         offset(:, 1) = state
         regression = bsts_spike_slab_draws(latent, predictors, slab_mean, &
            slab_covariance, prior_inclusion_probability, &
            observation_prior_shape, observation_prior_rate, observation_scale, &
            included, maximum_model_size, maximum_flips, 0, &
            coefficient_normal_draws(:, iteration:iteration), &
            inclusion_uniform_draws(:, iteration:iteration), &
            gamma_draws(1, iteration:iteration), offset)
         if (regression%info /= 0) then
            out%info = 10 + regression%info
            return
         end if
         coefficient = regression%coefficients(:, 1)
         included = regression%included(:, 1)
         observation_scale = regression%residual_variance(1)
         state_response = latent - matmul(predictors, coefficient)
         model = structural_model(state_response, [initial_level_mean], &
            reshape([initial_level_variance], [1, 1]), reshape([1.0_dp], [1, 1]), &
            [1.0_dp], observation_scale, [level_scale])
         call state_path_from_draws(model, state_normal_draws(:, :, iteration), &
            path, status)
         if (status /= 0) then
            out%info = 20 + status
            return
         end if
         state = path(1, :)
         level_sum = sum((state(2:) - state(:fine_points - 1))**2)
         level_scale = (level_prior_rate + 0.5_dp*level_sum)/ &
            gamma_draws(2, iteration)
         out%latent_fine(:, iteration) = latent
         out%cumulator(:, iteration) = mixed_frequency_cumulator(latent, &
            membership_fraction, contains_end)
         out%state(:, iteration) = state
         out%coefficients(:, iteration) = coefficient
         out%included(:, iteration) = included
         out%observation_variance(iteration) = observation_scale
         out%level_variance(iteration) = level_scale
         out%state_contribution(:, iteration) = state
         out%regression_contribution(:, iteration) = &
            matmul(predictors, coefficient)
         out%coarse_fitted(:, iteration) = matmul(design, latent)
      end do
      retained = iterations - burn
      do variable = 1, variables
         out%inclusion_probability(variable) = real(count( &
            out%included(variable, burn + 1:)), dp)/real(retained, dp)
      end do
   end function bsts_mixed_draws

   function bsts_mixed(coarse_response, predictors, coarse_index, &
      membership_fraction, contains_end, iterations, burn, &
      initial_level_mean, initial_level_variance, observation_variance, &
      level_variance, slab_mean, slab_covariance, prior_inclusion_probability, &
      expected_model_size, observation_prior_shape, observation_prior_rate, &
      level_prior_shape, level_prior_rate, initial_inclusion, &
      maximum_model_size, maximum_flips) result(out)
      !! Sample a local-level mixed-frequency model using random variates.
      real(dp), intent(in) :: coarse_response(:) !! Coarse response.
      real(dp), intent(in) :: predictors(:, :) !! Predictor matrix.
      integer, intent(in) :: coarse_index(:) !! Index of coarse.
      integer, intent(in) :: iterations !! Number of algorithm iterations.
      real(dp), intent(in) :: membership_fraction(:) !! Membership fraction.
      logical, intent(in) :: contains_end(:) !! Flag controlling contains end.
      integer, intent(in), optional :: burn !! Number of initial simulation draws to discard.
      integer, intent(in), optional :: maximum_model_size !! Maximum model size.
      integer, intent(in), optional :: maximum_flips !! Maximum flips.
      real(dp), intent(in), optional :: initial_level_mean !! Initial level mean.
      real(dp), intent(in), optional :: initial_level_variance !! Initial level variance.
      real(dp), intent(in), optional :: observation_variance !! Observation-error variance.
      real(dp), intent(in), optional :: level_variance !! Level variance.
      real(dp), intent(in), optional :: slab_mean(:) !! Slab mean.
      real(dp), intent(in), optional :: slab_covariance(:, :) !! Slab covariance matrix.
      real(dp), intent(in), optional :: prior_inclusion_probability(:) !! Prior inclusion probability.
      real(dp), intent(in), optional :: expected_model_size !! Expected model size.
      real(dp), intent(in), optional :: observation_prior_shape !! Observation prior shape.
      real(dp), intent(in), optional :: observation_prior_rate !! Observation prior rate.
      real(dp), intent(in), optional :: level_prior_shape !! Level prior shape.
      real(dp), intent(in), optional :: level_prior_rate !! Level prior rate.
      logical, intent(in), optional :: initial_inclusion(:) !! Initial inclusion.
      type(bsts_mixed_t) :: out
      real(dp), allocatable :: mean0(:), covariance0(:, :), inclusion0(:)
      real(dp), allocatable :: latent_normals(:, :), state_normals(:, :, :)
      real(dp), allocatable :: coefficient_normals(:, :), uniforms(:, :)
      real(dp), allocatable :: gammas(:, :)
      logical, allocatable :: model0(:)
      real(dp) :: level_mean0, initial0, observation0, level0, series_variance
      real(dp) :: observation_shape0, observation_rate0, level_shape0, level_rate0
      real(dp) :: expected
      integer :: fine_points, variables, discarded, model_limit, flip_limit
      integer :: iteration, variable

      fine_points = size(predictors, 1)
      variables = size(predictors, 2)
      if (fine_points < 2 .or. variables < 1 .or. iterations < 1 .or. &
         size(coarse_index) /= fine_points .or. &
         size(membership_fraction) /= fine_points .or. &
         size(contains_end) /= fine_points) then
         out%info = 1
         return
      end if
      discarded = 0
      if (present(burn)) discarded = burn
      model_limit = variables
      if (present(maximum_model_size)) model_limit = maximum_model_size
      flip_limit = -1
      if (present(maximum_flips)) flip_limit = maximum_flips
      series_variance = finite_sample_variance(coarse_response)/ &
         real(max(1, fine_points/max(1, size(coarse_response))), dp)
      level_mean0 = first_finite_value(coarse_response)/ &
         real(max(1, fine_points/max(1, size(coarse_response))), dp)
      if (present(initial_level_mean)) level_mean0 = initial_level_mean
      initial0 = max(series_variance, 1.0_dp)
      if (present(initial_level_variance)) initial0 = initial_level_variance
      observation0 = max(0.5_dp*series_variance, 1.0e-6_dp)
      if (present(observation_variance)) observation0 = observation_variance
      level0 = max(0.01_dp*series_variance, 1.0e-8_dp)
      if (present(level_variance)) level0 = level_variance
      observation_shape0 = 0.01_dp
      if (present(observation_prior_shape)) observation_shape0 = &
         observation_prior_shape
      observation_rate0 = max(0.01_dp*series_variance, 1.0e-8_dp)
      if (present(observation_prior_rate)) observation_rate0 = &
         observation_prior_rate
      level_shape0 = 0.01_dp
      if (present(level_prior_shape)) level_shape0 = level_prior_shape
      level_rate0 = max(0.0001_dp*series_variance, 1.0e-10_dp)
      if (present(level_prior_rate)) level_rate0 = level_prior_rate
      allocate(mean0(variables), covariance0(variables, variables), &
         inclusion0(variables), model0(variables))
      mean0 = 0.0_dp
      if (present(slab_mean)) then
         if (size(slab_mean) /= variables) then
            out%info = 1
            return
         end if
         mean0 = slab_mean
      end if
      covariance0 = 0.0_dp
      do variable = 1, variables
         covariance0(variable, variable) = real(fine_points, dp)/ &
            max(dot_product(predictors(:, variable), predictors(:, variable)), &
            1.0_dp)
      end do
      if (present(slab_covariance)) then
         if (any(shape(slab_covariance) /= [variables, variables])) then
            out%info = 1
            return
         end if
         covariance0 = slab_covariance
      end if
      expected = real(min(3, variables), dp)
      if (present(expected_model_size)) expected = expected_model_size
      inclusion0 = max(0.0_dp, min(1.0_dp, expected/real(variables, dp)))
      if (present(prior_inclusion_probability)) then
         if (size(prior_inclusion_probability) /= variables) then
            out%info = 1
            return
         end if
         inclusion0 = prior_inclusion_probability
      end if
      do variable = 1, variables
         if (sum((predictors(:, variable) - sum(predictors(:, variable))/ &
            real(fine_points, dp))**2) <= 100.0_dp*epsilon(1.0_dp)) &
            inclusion0(variable) = 0.0_dp
      end do
      model0 = inclusion0 >= 1.0_dp
      if (present(initial_inclusion)) then
         if (size(initial_inclusion) /= variables) then
            out%info = 1
            return
         end if
         model0 = initial_inclusion
      end if
      if (discarded < 0 .or. discarded >= iterations .or. &
         model_limit < 0 .or. model_limit > variables .or. flip_limit < -1 .or. &
         expected < 0.0_dp .or. initial0 <= 0.0_dp .or. &
         observation0 <= 0.0_dp .or. level0 <= 0.0_dp .or. &
         observation_shape0 <= 0.0_dp .or. observation_rate0 <= 0.0_dp .or. &
         level_shape0 <= 0.0_dp .or. level_rate0 <= 0.0_dp) then
         out%info = 1
         return
      end if
      allocate(latent_normals(fine_points, iterations), &
         state_normals(1, fine_points, iterations), &
         coefficient_normals(variables, iterations), &
         uniforms(variables, iterations), gammas(2, iterations))
      do iteration = 1, iterations
         do variable = 1, fine_points
            latent_normals(variable, iteration) = random_standard_normal()
            state_normals(1, variable, iteration) = random_standard_normal()
         end do
         do variable = 1, variables
            coefficient_normals(variable, iteration) = random_standard_normal()
            uniforms(variable, iteration) = max(tiny(1.0_dp), &
               min(1.0_dp - epsilon(1.0_dp), random_uniform()))
         end do
         gammas(1, iteration) = random_gamma(observation_shape0 + &
            0.5_dp*real(fine_points, dp))
         gammas(2, iteration) = random_gamma(level_shape0 + &
            0.5_dp*real(fine_points - 1, dp))
      end do
      out = bsts_mixed_draws(coarse_response, predictors, coarse_index, &
         membership_fraction, contains_end, level_mean0, initial0, observation0, &
         level0, mean0, covariance0, inclusion0, observation_shape0, &
         observation_rate0, level_shape0, level_rate0, model0, model_limit, &
         flip_limit, discarded, latent_normals, state_normals, &
         coefficient_normals, uniforms, gammas)
   end function bsts_mixed

   pure function bsts_mixed_trend_seasonal_draws(coarse_response, predictors, &
      coarse_index, membership_fraction, contains_end, nseasons, duration, &
      initial_trend_mean, initial_trend_covariance, seasonal_initial_variance, &
      observation_variance, component_variance, slab_mean, slab_covariance, &
      prior_inclusion_probability, observation_prior_shape, &
      observation_prior_rate, component_prior_shape, component_prior_rate, &
      initial_inclusion, maximum_model_size, maximum_flips, burn, &
      latent_normal_draws, structural_normal_draws, coefficient_normal_draws, &
      inclusion_uniform_draws, gamma_draws) result(out)
      !! Sample a mixed-frequency trend-seasonal model from supplied variates.
      real(dp), intent(in) :: coarse_response(:) !! Coarse response.
      real(dp), intent(in) :: predictors(:, :) !! Predictor matrix.
      integer, intent(in) :: coarse_index(:) !! Index of coarse.
      integer, intent(in) :: nseasons !! Nseasons.
      integer, intent(in) :: duration !! Duration.
      real(dp), intent(in) :: membership_fraction(:) !! Membership fraction.
      logical, intent(in) :: contains_end(:) !! Flag controlling contains end.
      real(dp), intent(in) :: initial_trend_mean(2) !! Initial trend mean.
      real(dp), intent(in) :: initial_trend_covariance(2, 2) !! Initial trend covariance.
      real(dp), intent(in) :: seasonal_initial_variance !! Seasonal initial variance.
      real(dp), intent(in) :: observation_variance !! Observation-error variance.
      real(dp), intent(in) :: component_variance(3) !! Component variance.
      real(dp), intent(in) :: slab_mean(:) !! Slab mean.
      real(dp), intent(in) :: slab_covariance(:, :) !! Slab covariance matrix.
      real(dp), intent(in) :: prior_inclusion_probability(:) !! Prior inclusion probability.
      real(dp), intent(in) :: observation_prior_shape !! Observation prior shape.
      real(dp), intent(in) :: observation_prior_rate !! Observation prior rate.
      real(dp), intent(in) :: component_prior_shape(3) !! Component prior shape.
      real(dp), intent(in) :: component_prior_rate(3) !! Component prior rate.
      logical, intent(in) :: initial_inclusion(:) !! Initial inclusion.
      integer, intent(in) :: maximum_model_size !! Maximum model size.
      integer, intent(in) :: maximum_flips !! Maximum flips.
      integer, intent(in) :: burn !! Number of initial simulation draws to discard.
      real(dp), intent(in) :: latent_normal_draws(:, :) !! Latent normal simulation draws.
      real(dp), intent(in) :: structural_normal_draws(:, :, :) !! Structural normal simulation draws.
      real(dp), intent(in) :: coefficient_normal_draws(:, :) !! Coefficient normal simulation draws.
      real(dp), intent(in) :: inclusion_uniform_draws(:, :) !! Inclusion uniform simulation draws.
      real(dp), intent(in) :: gamma_draws(:, :) !! Gamma simulation draws.
      type(bsts_mixed_t) :: out
      type(bsts_spike_slab_t) :: regression
      type(ssm_model_t) :: model
      real(dp), allocatable :: design(:, :), observed_design(:, :)
      real(dp), allocatable :: constraint_covariance(:, :), constraint_inverse(:, :)
      real(dp), allocatable :: target(:), latent(:), raw_latent(:)
      real(dp), allocatable :: contribution(:), coefficient(:), state_response(:)
      real(dp), allocatable :: path(:, :), initial_mean(:), initial_covariance(:, :)
      real(dp), allocatable :: offset(:, :)
      integer, allocatable :: observed_coarse(:)
      logical, allocatable :: included(:)
      real(dp) :: observation_scale, variance(3), residual_sum
      integer :: fine_points, coarse_points, variables, dimension, iterations
      integer :: iteration, row, time, status, retained, variable, component

      fine_points = size(predictors, 1)
      coarse_points = size(coarse_response)
      variables = size(predictors, 2)
      dimension = nseasons + 1
      iterations = size(latent_normal_draws, 2)
      if (fine_points < 2 .or. coarse_points < 1 .or. variables < 1 .or. &
         nseasons < 2 .or. duration < 1 .or. &
         size(coarse_index) /= fine_points .or. &
         size(membership_fraction) /= fine_points .or. &
         size(contains_end) /= fine_points .or. &
         any(membership_fraction < 0.0_dp) .or. &
         any(membership_fraction > 1.0_dp) .or. &
         .not. all(ieee_is_finite(predictors)) .or. &
         size(slab_mean) /= variables .or. &
         any(shape(slab_covariance) /= [variables, variables]) .or. &
         size(prior_inclusion_probability) /= variables .or. &
         size(initial_inclusion) /= variables .or. &
         any(prior_inclusion_probability < 0.0_dp) .or. &
         any(prior_inclusion_probability > 1.0_dp) .or. &
         seasonal_initial_variance <= 0.0_dp .or. &
         observation_variance <= 0.0_dp .or. &
         any(component_variance <= 0.0_dp) .or. &
         observation_prior_shape <= 0.0_dp .or. observation_prior_rate <= 0.0_dp .or. &
         any(component_prior_shape <= 0.0_dp) .or. &
         any(component_prior_rate <= 0.0_dp) .or. maximum_model_size < 0 .or. &
         maximum_model_size > variables .or. maximum_flips < -1 .or. &
         iterations < 1 .or. burn < 0 .or. burn >= iterations .or. &
         size(latent_normal_draws, 1) /= fine_points .or. &
         any(shape(structural_normal_draws) /= &
         [dimension, fine_points, iterations]) .or. &
         any(shape(coefficient_normal_draws) /= [variables, iterations]) .or. &
         any(shape(inclusion_uniform_draws) /= [variables, iterations]) .or. &
         any(shape(gamma_draws) /= [4, iterations]) .or. &
         any(inclusion_uniform_draws <= 0.0_dp) .or. &
         any(inclusion_uniform_draws >= 1.0_dp) .or. any(gamma_draws <= 0.0_dp) .or. &
         .not. all(ieee_is_finite(latent_normal_draws)) .or. &
         .not. all(ieee_is_finite(structural_normal_draws)) .or. &
         .not. all(ieee_is_finite(coefficient_normal_draws))) then
         out%info = 1
         return
      end if
      design = mixed_frequency_design(coarse_index, membership_fraction, &
         contains_end, coarse_points)
      do row = 1, coarse_points
         if (ieee_is_finite(coarse_response(row)) .and. &
            sum(abs(design(row, :))) <= tiny(1.0_dp)) then
            out%info = 2
            return
         end if
      end do
      observed_coarse = pack([(row, row=1, coarse_points)], &
         ieee_is_finite(coarse_response))
      if (size(observed_coarse) > 0) then
         observed_design = design(observed_coarse, :)
         target = coarse_response(observed_coarse)
         constraint_covariance = matmul(observed_design, &
            transpose(observed_design))
         call invert_matrix(constraint_covariance, constraint_inverse, status)
         if (status /= 0) then
            out%info = 3
            return
         end if
      else
         allocate(observed_design(0, fine_points), target(0), &
            constraint_inverse(0, 0))
      end if
      allocate(latent(fine_points), raw_latent(fine_points), &
         contribution(fine_points), coefficient(variables), &
         included(variables), state_response(fine_points), offset(fine_points, 1), &
         initial_mean(dimension), initial_covariance(dimension, dimension), &
         out%latent_fine(fine_points, iterations), &
         out%cumulator(fine_points, iterations), out%state(fine_points, iterations), &
         out%coefficients(variables, iterations), &
         out%included(variables, iterations), &
         out%observation_variance(iterations), out%level_variance(iterations), &
         out%state_contribution(fine_points, iterations), &
         out%regression_contribution(fine_points, iterations), &
         out%coarse_fitted(coarse_points, iterations), &
         out%structural_state(dimension, fine_points, iterations), &
         out%component_variance(3, iterations), &
         out%inclusion_probability(variables))
      initial_mean = 0.0_dp
      initial_mean(:2) = initial_trend_mean
      initial_covariance = 0.0_dp
      initial_covariance(:2, :2) = initial_trend_covariance
      do component = 3, dimension
         initial_covariance(component, component) = seasonal_initial_variance
      end do
      contribution = initial_trend_mean(1)
      coefficient = 0.0_dp
      included = initial_inclusion
      where (prior_inclusion_probability <= 0.0_dp) included = .false.
      where (prior_inclusion_probability >= 1.0_dp) included = .true.
      if (count(included) > maximum_model_size) then
         out%info = 4
         return
      end if
      observation_scale = observation_variance
      variance = component_variance
      out%burn = burn
      out%nseasons = nseasons
      out%season_duration = duration
      do iteration = 1, iterations
         raw_latent = contribution + matmul(predictors, coefficient) + &
            sqrt(observation_scale)*latent_normal_draws(:, iteration)
         latent = raw_latent
         if (size(observed_coarse) > 0) latent = raw_latent + &
            matmul(transpose(observed_design), matmul(constraint_inverse, &
            target - matmul(observed_design, raw_latent)))
         offset(:, 1) = contribution
         regression = bsts_spike_slab_draws(latent, predictors, slab_mean, &
            slab_covariance, prior_inclusion_probability, &
            observation_prior_shape, observation_prior_rate, observation_scale, &
            included, maximum_model_size, maximum_flips, 0, &
            coefficient_normal_draws(:, iteration:iteration), &
            inclusion_uniform_draws(:, iteration:iteration), &
            gamma_draws(1, iteration:iteration), offset)
         if (regression%info /= 0) then
            out%info = 10 + regression%info
            return
         end if
         coefficient = regression%coefficients(:, 1)
         included = regression%included(:, 1)
         observation_scale = regression%residual_variance(1)
         state_response = latent - matmul(predictors, coefficient)
         model = mixed_trend_seasonal_model(state_response, nseasons, duration, &
            initial_mean, initial_covariance, observation_scale, variance)
         call state_path_from_draws(model, &
            structural_normal_draws(:, :, iteration), path, status)
         if (status /= 0) then
            out%info = 20 + status
            return
         end if
         contribution = path(1, :) + path(3, :)
         residual_sum = sum((path(1, 2:) - path(1, :fine_points - 1) - &
            path(2, :fine_points - 1))**2)
         variance(1) = (component_prior_rate(1) + 0.5_dp*residual_sum)/ &
            gamma_draws(2, iteration)
         residual_sum = sum((path(2, 2:) - path(2, :fine_points - 1))**2)
         variance(2) = (component_prior_rate(2) + 0.5_dp*residual_sum)/ &
            gamma_draws(3, iteration)
         residual_sum = 0.0_dp
         do time = 2, fine_points
            if (modulo(time - 1, duration) == 0) residual_sum = residual_sum + &
               (path(3, time) + sum(path(3:, time - 1)))**2
         end do
         variance(3) = (component_prior_rate(3) + 0.5_dp*residual_sum)/ &
            gamma_draws(4, iteration)
         out%latent_fine(:, iteration) = latent
         out%cumulator(:, iteration) = mixed_frequency_cumulator(latent, &
            membership_fraction, contains_end)
         out%state(:, iteration) = contribution
         out%structural_state(:, :, iteration) = path
         out%coefficients(:, iteration) = coefficient
         out%included(:, iteration) = included
         out%observation_variance(iteration) = observation_scale
         out%level_variance(iteration) = variance(1)
         out%component_variance(:, iteration) = variance
         out%state_contribution(:, iteration) = contribution
         out%regression_contribution(:, iteration) = matmul(predictors, coefficient)
         out%coarse_fitted(:, iteration) = matmul(design, latent)
      end do
      retained = iterations - burn
      do variable = 1, variables
         out%inclusion_probability(variable) = real(count( &
            out%included(variable, burn + 1:)), dp)/real(retained, dp)
      end do
   end function bsts_mixed_trend_seasonal_draws

   function bsts_mixed_trend_seasonal(coarse_response, predictors, &
      coarse_index, membership_fraction, contains_end, nseasons, duration, &
      iterations, burn, prior_inclusion_probability, expected_model_size, &
      maximum_model_size, maximum_flips) result(out)
      !! Sample a mixed-frequency trend-seasonal model using random variates.
      real(dp), intent(in) :: coarse_response(:) !! Coarse response.
      real(dp), intent(in) :: predictors(:, :) !! Predictor matrix.
      integer, intent(in) :: coarse_index(:) !! Index of coarse.
      integer, intent(in) :: nseasons !! Nseasons.
      integer, intent(in) :: duration !! Duration.
      integer, intent(in) :: iterations !! Number of algorithm iterations.
      real(dp), intent(in) :: membership_fraction(:) !! Membership fraction.
      logical, intent(in) :: contains_end(:) !! Flag controlling contains end.
      integer, intent(in), optional :: burn !! Number of initial simulation draws to discard.
      integer, intent(in), optional :: maximum_model_size !! Maximum model size.
      integer, intent(in), optional :: maximum_flips !! Maximum flips.
      real(dp), intent(in), optional :: prior_inclusion_probability(:) !! Prior inclusion probability.
      real(dp), intent(in), optional :: expected_model_size !! Expected model size.
      type(bsts_mixed_t) :: out
      real(dp), allocatable :: inclusion(:), slab_covariance(:, :)
      real(dp), allocatable :: latent_normals(:, :), structural_normals(:, :, :)
      real(dp), allocatable :: coefficient_normals(:, :), uniforms(:, :)
      real(dp), allocatable :: gammas(:, :)
      logical, allocatable :: initial_inclusion(:)
      real(dp) :: series_variance, expected
      integer :: fine_points, variables, dimension, discarded
      integer :: model_limit, flip_limit, iteration, variable, component

      fine_points = size(predictors, 1)
      variables = size(predictors, 2)
      dimension = nseasons + 1
      if (fine_points < 2 .or. variables < 1 .or. nseasons < 2 .or. &
         duration < 1 .or. iterations < 1) then
         out%info = 1
         return
      end if
      discarded = 0
      if (present(burn)) discarded = burn
      model_limit = variables
      if (present(maximum_model_size)) model_limit = maximum_model_size
      flip_limit = -1
      if (present(maximum_flips)) flip_limit = maximum_flips
      expected = real(min(3, variables), dp)
      if (present(expected_model_size)) expected = expected_model_size
      allocate(inclusion(variables), slab_covariance(variables, variables), &
         initial_inclusion(variables))
      inclusion = max(0.0_dp, min(1.0_dp, expected/real(variables, dp)))
      if (present(prior_inclusion_probability)) then
         if (size(prior_inclusion_probability) /= variables) then
            out%info = 1
            return
         end if
         inclusion = prior_inclusion_probability
      end if
      slab_covariance = 0.0_dp
      do variable = 1, variables
         slab_covariance(variable, variable) = real(fine_points, dp)/ &
            max(dot_product(predictors(:, variable), predictors(:, variable)), &
            1.0_dp)
         if (sum((predictors(:, variable) - sum(predictors(:, variable))/ &
            real(fine_points, dp))**2) <= 100.0_dp*epsilon(1.0_dp)) &
            inclusion(variable) = 0.0_dp
      end do
      initial_inclusion = inclusion >= 1.0_dp
      series_variance = finite_sample_variance(coarse_response)/ &
         real(max(1, fine_points/max(1, size(coarse_response))), dp)
      if (discarded < 0 .or. discarded >= iterations .or. expected < 0.0_dp .or. &
         model_limit < 0 .or. model_limit > variables .or. flip_limit < -1) then
         out%info = 1
         return
      end if
      allocate(latent_normals(fine_points, iterations), &
         structural_normals(dimension, fine_points, iterations), &
         coefficient_normals(variables, iterations), &
         uniforms(variables, iterations), gammas(4, iterations))
      do iteration = 1, iterations
         do variable = 1, fine_points
            latent_normals(variable, iteration) = random_standard_normal()
            do component = 1, dimension
               structural_normals(component, variable, iteration) = &
                  random_standard_normal()
            end do
         end do
         do variable = 1, variables
            coefficient_normals(variable, iteration) = random_standard_normal()
            uniforms(variable, iteration) = max(tiny(1.0_dp), &
               min(1.0_dp - epsilon(1.0_dp), random_uniform()))
         end do
         gammas(1, iteration) = random_gamma(2.0_dp + &
            0.5_dp*real(fine_points, dp))
         gammas(2, iteration) = random_gamma(2.0_dp + &
            0.5_dp*real(fine_points - 1, dp))
         gammas(3, iteration) = random_gamma(2.0_dp + &
            0.5_dp*real(fine_points - 1, dp))
         gammas(4, iteration) = random_gamma(2.0_dp + &
            0.5_dp*real((fine_points - 1)/duration, dp))
      end do
      out = bsts_mixed_trend_seasonal_draws(coarse_response, predictors, &
         coarse_index, membership_fraction, contains_end, nseasons, duration, &
         [first_finite_value(coarse_response)/real(max(1, fine_points/ &
         max(1, size(coarse_response))), dp), 0.0_dp], &
         reshape([max(series_variance, 1.0_dp), 0.0_dp, 0.0_dp, 1.0_dp], &
         [2, 2]), 1.0_dp, max(0.5_dp*series_variance, 1.0e-6_dp), &
         [max(0.01_dp*series_variance, 1.0e-8_dp), &
         max(0.001_dp*series_variance, 1.0e-10_dp), &
         max(0.001_dp*series_variance, 1.0e-10_dp)], &
         [(0.0_dp, variable=1, variables)], slab_covariance, inclusion, &
         2.0_dp, max(0.01_dp*series_variance, 1.0e-8_dp), &
         [2.0_dp, 2.0_dp, 2.0_dp], &
         [max(0.0001_dp*series_variance, 1.0e-10_dp), &
         max(0.00001_dp*series_variance, 1.0e-12_dp), &
         max(0.00001_dp*series_variance, 1.0e-12_dp)], &
         initial_inclusion, model_limit, flip_limit, discarded, latent_normals, &
         structural_normals, coefficient_normals, uniforms, gammas)
   end function bsts_mixed_trend_seasonal

   pure function bsts_mixed_trend_seasonal_predict_draws(fit, &
      future_predictors, future_coarse_index, membership_fraction, contains_end, &
      structural_normal_draws, observation_normal_draws) result(out)
      !! Forecast a mixed trend-seasonal model from supplied random variates.
      type(bsts_mixed_t), intent(in) :: fit !! Previously fitted model.
      real(dp), intent(in) :: future_predictors(:, :) !! Future predictors.
      integer, intent(in) :: future_coarse_index(:) !! Future coarse index.
      real(dp), intent(in) :: membership_fraction(:) !! Membership fraction.
      logical, intent(in) :: contains_end(:) !! Flag controlling contains end.
      real(dp), intent(in) :: structural_normal_draws(:, :, :) !! Structural normal simulation draws.
      real(dp), intent(in) :: observation_normal_draws(:, :) !! Observation normal draws.
      type(bsts_mixed_prediction_t) :: out
      real(dp), allocatable :: design(:, :), structural(:), next_seasonal(:)
      integer :: fine_points, coarse_points, retained, dimension
      integer :: draw, source, time, component

      fine_points = size(future_predictors, 1)
      retained = 0
      if (allocated(fit%coefficients)) retained = &
         size(fit%coefficients, 2) - fit%burn
      dimension = fit%nseasons + 1
      if (fit%info /= 0 .or. fit%nseasons < 2 .or. fine_points < 1 .or. &
         retained < 1 .or. size(future_coarse_index) /= fine_points .or. &
         size(membership_fraction) /= fine_points .or. &
         size(contains_end) /= fine_points .or. &
         size(future_predictors, 2) /= size(fit%coefficients, 1) .or. &
         minval(future_coarse_index) < 1 .or. &
         any(shape(structural_normal_draws) /= [3, fine_points, retained]) .or. &
         any(shape(observation_normal_draws) /= [fine_points, retained]) .or. &
         .not. all(ieee_is_finite(future_predictors)) .or. &
         .not. all(ieee_is_finite(structural_normal_draws)) .or. &
         .not. all(ieee_is_finite(observation_normal_draws))) then
         out%info = 1
         return
      end if
      coarse_points = maxval(future_coarse_index)
      design = mixed_frequency_design(future_coarse_index, membership_fraction, &
         contains_end, coarse_points)
      allocate(out%fine%draws(fine_points, retained), &
         out%coarse%draws(coarse_points, retained), structural(dimension), &
         next_seasonal(fit%nseasons - 1))
      do draw = 1, retained
         source = fit%burn + draw
         structural = fit%structural_state(:, &
            size(fit%structural_state, 2), source)
         do time = 1, fine_points
            structural(1) = structural(1) + structural(2) + &
               sqrt(fit%component_variance(1, source))* &
               structural_normal_draws(1, time, draw)
            structural(2) = structural(2) + &
               sqrt(fit%component_variance(2, source))* &
               structural_normal_draws(2, time, draw)
            if (modulo(size(fit%structural_state, 2) + time - 1, &
               fit%season_duration) == 0) then
               next_seasonal = 0.0_dp
               next_seasonal(1) = -sum(structural(3:)) + &
                  sqrt(fit%component_variance(3, source))* &
                  structural_normal_draws(3, time, draw)
               do component = 2, size(next_seasonal)
                  next_seasonal(component) = structural(component + 1)
               end do
               structural(3:) = next_seasonal
            end if
            out%fine%draws(time, draw) = structural(1) + structural(3) + &
               dot_product(future_predictors(time, :), &
               fit%coefficients(:, source)) + &
               sqrt(fit%observation_variance(source))* &
               observation_normal_draws(time, draw)
         end do
         out%coarse%draws(:, draw) = matmul(design, out%fine%draws(:, draw))
      end do
      call summarize_prediction(out%fine)
      call summarize_prediction(out%coarse)
   end function bsts_mixed_trend_seasonal_predict_draws

   function bsts_mixed_trend_seasonal_predict(fit, future_predictors, &
      future_coarse_index, membership_fraction, contains_end) result(out)
      !! Forecast a mixed trend-seasonal model using the shared random stream.
      type(bsts_mixed_t), intent(in) :: fit !! Previously fitted model.
      real(dp), intent(in) :: future_predictors(:, :) !! Future predictors.
      integer, intent(in) :: future_coarse_index(:) !! Future coarse index.
      real(dp), intent(in) :: membership_fraction(:) !! Membership fraction.
      logical, intent(in) :: contains_end(:) !! Flag controlling contains end.
      type(bsts_mixed_prediction_t) :: out
      real(dp), allocatable :: structural_normals(:, :, :)
      real(dp), allocatable :: observation_normals(:, :)
      integer :: retained, draw, time, component

      if (.not. allocated(fit%coefficients)) then
         out%info = 1
         return
      end if
      retained = size(fit%coefficients, 2) - fit%burn
      if (size(future_predictors, 1) < 1 .or. retained < 1) then
         out%info = 1
         return
      end if
      allocate(structural_normals(3, size(future_predictors, 1), retained), &
         observation_normals(size(future_predictors, 1), retained))
      do draw = 1, retained
         do time = 1, size(future_predictors, 1)
            do component = 1, 3
               structural_normals(component, time, draw) = &
                  random_standard_normal()
            end do
            observation_normals(time, draw) = random_standard_normal()
         end do
      end do
      out = bsts_mixed_trend_seasonal_predict_draws(fit, future_predictors, &
         future_coarse_index, membership_fraction, contains_end, &
         structural_normals, observation_normals)
   end function bsts_mixed_trend_seasonal_predict

   pure function bsts_mixed_predict_draws(fit, future_predictors, &
      future_coarse_index, membership_fraction, contains_end, &
      state_normal_draws, observation_normal_draws) result(out)
      !! Forecast fine and aggregated mixed-frequency observations.
      type(bsts_mixed_t), intent(in) :: fit !! Previously fitted model.
      real(dp), intent(in) :: future_predictors(:, :) !! Future predictors.
      integer, intent(in) :: future_coarse_index(:) !! Future coarse index.
      real(dp), intent(in) :: membership_fraction(:) !! Membership fraction.
      logical, intent(in) :: contains_end(:) !! Flag controlling contains end.
      real(dp), intent(in) :: state_normal_draws(:, :) !! Independent standard-normal state draws.
      real(dp), intent(in) :: observation_normal_draws(:, :) !! Observation normal draws.
      type(bsts_mixed_prediction_t) :: out
      real(dp), allocatable :: design(:, :)
      real(dp) :: level
      integer :: fine_points, coarse_points, retained, draw, source, time

      fine_points = size(future_predictors, 1)
      retained = 0
      if (allocated(fit%coefficients)) retained = &
         size(fit%coefficients, 2) - fit%burn
      if (fit%info /= 0 .or. fine_points < 1 .or. retained < 1 .or. &
         size(future_coarse_index) /= fine_points .or. &
         size(membership_fraction) /= fine_points .or. &
         size(contains_end) /= fine_points .or. &
         size(future_predictors, 2) /= size(fit%coefficients, 1) .or. &
         minval(future_coarse_index) < 1 .or. &
         any(membership_fraction < 0.0_dp) .or. &
         any(membership_fraction > 1.0_dp) .or. &
         any(shape(state_normal_draws) /= [fine_points, retained]) .or. &
         any(shape(observation_normal_draws) /= [fine_points, retained]) .or. &
         .not. all(ieee_is_finite(future_predictors)) .or. &
         .not. all(ieee_is_finite(state_normal_draws)) .or. &
         .not. all(ieee_is_finite(observation_normal_draws))) then
         out%info = 1
         return
      end if
      coarse_points = maxval(future_coarse_index)
      design = mixed_frequency_design(future_coarse_index, membership_fraction, &
         contains_end, coarse_points)
      allocate(out%fine%draws(fine_points, retained), &
         out%coarse%draws(coarse_points, retained))
      do draw = 1, retained
         source = fit%burn + draw
         level = fit%state(size(fit%state, 1), source)
         do time = 1, fine_points
            level = level + sqrt(fit%level_variance(source))* &
               state_normal_draws(time, draw)
            out%fine%draws(time, draw) = level + &
               dot_product(future_predictors(time, :), &
               fit%coefficients(:, source)) + &
               sqrt(fit%observation_variance(source))* &
               observation_normal_draws(time, draw)
         end do
         out%coarse%draws(:, draw) = matmul(design, out%fine%draws(:, draw))
      end do
      call summarize_prediction(out%fine)
      call summarize_prediction(out%coarse)
   end function bsts_mixed_predict_draws

   function bsts_mixed_predict(fit, future_predictors, future_coarse_index, &
      membership_fraction, contains_end) result(out)
      !! Forecast a mixed-frequency model using the shared random stream.
      type(bsts_mixed_t), intent(in) :: fit !! Previously fitted model.
      real(dp), intent(in) :: future_predictors(:, :) !! Future predictors.
      integer, intent(in) :: future_coarse_index(:) !! Future coarse index.
      real(dp), intent(in) :: membership_fraction(:) !! Membership fraction.
      logical, intent(in) :: contains_end(:) !! Flag controlling contains end.
      type(bsts_mixed_prediction_t) :: out
      real(dp), allocatable :: state_normals(:, :), observation_normals(:, :)
      integer :: retained, draw, time

      if (.not. allocated(fit%coefficients)) then
         out%info = 1
         return
      end if
      retained = size(fit%coefficients, 2) - fit%burn
      if (size(future_predictors, 1) < 1 .or. retained < 1) then
         out%info = 1
         return
      end if
      allocate(state_normals(size(future_predictors, 1), retained), &
         observation_normals(size(future_predictors, 1), retained))
      do draw = 1, retained
         do time = 1, size(future_predictors, 1)
            state_normals(time, draw) = random_standard_normal()
            observation_normals(time, draw) = random_standard_normal()
         end do
      end do
      out = bsts_mixed_predict_draws(fit, future_predictors, &
         future_coarse_index, membership_fraction, contains_end, &
         state_normals, observation_normals)
   end function bsts_mixed_predict

   pure function bsts_ar_draws(y, initial_mean, initial_covariance, &
      initial_ar_coefficients, observation_variance, innovation_variance, &
      observation_prior_shape, observation_prior_rate, innovation_prior_shape, &
      innovation_prior_rate, burn, state_normal_draws, ar_normal_draws, &
      gamma_draws, offset_draws) result(out)
      !! Sample a stationary AR(p) state component from supplied random draws.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: initial_ar_coefficients(:) !! Initial autoregressive coefficients.
      real(dp), intent(in) :: observation_variance !! Observation-error variance.
      real(dp), intent(in) :: innovation_variance !! Innovation variance.
      real(dp), intent(in) :: observation_prior_shape !! Observation prior shape.
      real(dp), intent(in) :: observation_prior_rate !! Observation prior rate.
      real(dp), intent(in) :: innovation_prior_shape !! Innovation prior shape.
      real(dp), intent(in) :: innovation_prior_rate !! Innovation prior rate.
      integer, intent(in) :: burn !! Number of initial simulation draws to discard.
      real(dp), intent(in) :: state_normal_draws(:, :, :) !! Independent standard-normal state draws.
      real(dp), intent(in) :: ar_normal_draws(:, :, :) !! Autoregressive normal simulation draws.
      real(dp), intent(in) :: gamma_draws(:, :) !! Gamma simulation draws.
      real(dp), intent(in), optional :: offset_draws(:, :) !! Offset simulation draws.
      type(bsts_mcmc_t) :: out
      type(ssm_model_t) :: model
      real(dp), allocatable :: path(:, :), working_y(:), design(:, :)
      real(dp), allocatable :: target(:), crossproduct(:, :), inverse(:, :)
      real(dp), allocatable :: posterior_mean(:), candidate(:)
      real(dp) :: ar(size(initial_ar_coefficients)), observation_scale
      real(dp) :: innovation_scale, observation_sum, innovation_sum
      integer :: order, observations, iterations, attempts, iteration, time
      integer :: proposal, status

      order = size(initial_ar_coefficients)
      observations = size(y)
      iterations = size(state_normal_draws, 3)
      attempts = size(ar_normal_draws, 2)
      if (observations < 2 .or. count(ieee_is_finite(y)) < 1 .or. order < 1 .or. &
         size(initial_mean) /= order .or. &
         any(shape(initial_covariance) /= [order, order]) .or. &
         .not. ar_coefficients_stationary(initial_ar_coefficients) .or. &
         observation_variance <= 0.0_dp .or. innovation_variance <= 0.0_dp .or. &
         observation_prior_shape <= 0.0_dp .or. observation_prior_rate <= 0.0_dp .or. &
         innovation_prior_shape <= 0.0_dp .or. innovation_prior_rate <= 0.0_dp .or. &
         iterations < 1 .or. attempts < 1 .or. burn < 0 .or. burn >= iterations .or. &
         any(shape(state_normal_draws) /= [order, observations, iterations]) .or. &
         any(shape(ar_normal_draws) /= [order, attempts, iterations]) .or. &
         any(shape(gamma_draws) /= [2, iterations]) .or. &
         any(gamma_draws <= 0.0_dp) .or. &
         .not. all(ieee_is_finite(state_normal_draws)) .or. &
         .not. all(ieee_is_finite(ar_normal_draws))) then
         out%info = 1
         return
      end if
      if (present(offset_draws)) then
         if (any(shape(offset_draws) /= [observations, iterations]) .or. &
            .not. all(ieee_is_finite(offset_draws))) then
            out%info = 1
            return
         end if
      end if
      allocate(model%y(observations, 1), model%z(1, order, 1), &
         model%h(1, 1, 1), model%transition(order, order, 1), &
         model%r(order, 1, 1), model%q(1, 1, 1), model%a1(order), &
         model%p1(order, order), model%p1inf(order, order), &
         model%missing(observations, 1), working_y(observations), &
         design(observations - 1, order), target(observations - 1), &
         crossproduct(order, order), inverse(order, order), &
         posterior_mean(order), candidate(order), &
         out%state(order, observations, iterations), &
         out%observation_variance(iterations), &
         out%state_variance(order, iterations), &
         out%component_variance(1, iterations), &
         out%ar_coefficients(order, iterations), &
         out%transition(order, order), out%observation(order), &
         out%transition_schedule(order, order, 1), &
         out%state_loading_schedule(order, 1, 1))
      model%z = 0.0_dp
      model%z(1, 1, 1) = 1.0_dp
      model%r = 0.0_dp
      model%r(1, 1, 1) = 1.0_dp
      model%a1 = initial_mean
      model%p1 = initial_covariance
      model%p1inf = 0.0_dp
      model%missing(:, 1) = .not. ieee_is_finite(y)
      out%observation = 0.0_dp
      out%observation(1) = 1.0_dp
      out%state_loading_schedule = 0.0_dp
      out%state_loading_schedule(1, 1, 1) = 1.0_dp
      out%state_variance = 0.0_dp
      out%burn = burn
      out%ar_order = order
      ar = initial_ar_coefficients
      observation_scale = observation_variance
      innovation_scale = innovation_variance
      do iteration = 1, iterations
         working_y = y
         if (present(offset_draws)) working_y = working_y - &
            offset_draws(:, iteration)
         model%y(:, 1) = working_y
         where (.not. ieee_is_finite(model%y)) model%y = 0.0_dp
         model%h(1, 1, 1) = observation_scale
         call set_ar_transition(model%transition(:, :, 1), ar)
         model%q(1, 1, 1) = innovation_scale
         call state_path_from_draws(model, state_normal_draws(:, :, iteration), &
            path, status)
         if (status /= 0) then
            out%info = 10 + status
            return
         end if
         out%state(:, :, iteration) = path
         observation_sum = 0.0_dp
         do time = 1, observations
            if (ieee_is_finite(y(time))) observation_sum = observation_sum + &
               (working_y(time) - path(1, time))**2
         end do
         observation_scale = (observation_prior_rate + &
            0.5_dp*observation_sum)/gamma_draws(1, iteration)
         target = path(1, 2:)
         design = transpose(path(:, :observations - 1))
         crossproduct = matmul(transpose(design), design)
         call invert_matrix(crossproduct, inverse, status)
         if (status /= 0) then
            out%info = 50 + status
            return
         end if
         posterior_mean = matmul(inverse, matmul(transpose(design), target))
         do proposal = 1, attempts
            call multivariate_normal_from_standard(posterior_mean, &
               innovation_scale*inverse, ar_normal_draws(:, proposal, iteration), &
               candidate, status)
            if (status /= 0) then
               out%info = 60 + status
               return
            end if
            if (ar_coefficients_stationary(candidate)) then
               ar = candidate
               exit
            end if
         end do
         innovation_sum = sum((target - matmul(design, ar))**2)
         innovation_scale = (innovation_prior_rate + 0.5_dp*innovation_sum)/ &
            gamma_draws(2, iteration)
         out%observation_variance(iteration) = observation_scale
         out%state_variance(1, iteration) = innovation_scale
         out%component_variance(1, iteration) = innovation_scale
         out%ar_coefficients(:, iteration) = ar
      end do
      call set_ar_transition(out%transition, ar)
      out%transition_schedule(:, :, 1) = out%transition
   end function bsts_ar_draws

   function bsts_ar(y, lags, iterations, burn, proposal_attempts, offset_draws) &
      result(out)
      !! Sample a stationary AR(p) state component using the random stream.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      integer, intent(in) :: lags !! Lags.
      integer, intent(in) :: iterations !! Number of algorithm iterations.
      integer, intent(in), optional :: burn !! Number of initial simulation draws to discard.
      integer, intent(in), optional :: proposal_attempts !! Proposal attempts.
      real(dp), intent(in), optional :: offset_draws(:, :) !! Offset simulation draws.
      type(bsts_mcmc_t) :: out
      real(dp), allocatable :: mean0(:), covariance0(:, :), ar0(:)
      real(dp), allocatable :: state_normals(:, :, :), ar_normals(:, :, :)
      real(dp), allocatable :: gammas(:, :)
      real(dp) :: series_variance, observation0, innovation0
      real(dp) :: observation_shape, observation_rate
      real(dp) :: innovation_shape, innovation_rate
      integer :: discarded, attempts, iteration, proposal, observed

      if (size(y) < 2 .or. count(ieee_is_finite(y)) < 1 .or. lags < 1 .or. &
         iterations < 1) then
         out%info = 1
         return
      end if
      discarded = 0
      if (present(burn)) discarded = burn
      attempts = 100
      if (present(proposal_attempts)) attempts = proposal_attempts
      if (discarded < 0 .or. discarded >= iterations .or. attempts < 1) then
         out%info = 1
         return
      end if
      series_variance = max(finite_sample_variance(y), 1.0e-6_dp)
      observation0 = max(0.1_dp*series_variance, 1.0e-8_dp)
      innovation0 = max(0.01_dp*series_variance, 1.0e-10_dp)
      observation_shape = 2.0_dp
      observation_rate = observation0
      innovation_shape = 2.0_dp
      innovation_rate = innovation0
      allocate(mean0(lags), covariance0(lags, lags), ar0(lags), &
         state_normals(lags, size(y), iterations), &
         ar_normals(lags, attempts, iterations), gammas(2, iterations))
      mean0 = 0.0_dp
      mean0(1) = first_finite_value(y)
      covariance0 = identity_block(lags, series_variance)
      ar0 = 0.0_dp
      if (lags == 1) ar0(1) = 0.5_dp
      observed = count(ieee_is_finite(y))
      do iteration = 1, iterations
         call fill_standard_normals(state_normals(:, :, iteration))
         do proposal = 1, attempts
            call fill_vector_normals(ar_normals(:, proposal, iteration))
         end do
         gammas(1, iteration) = random_gamma(observation_shape + &
            0.5_dp*real(observed, dp))
         gammas(2, iteration) = random_gamma(innovation_shape + &
            0.5_dp*real(size(y) - 1, dp))
      end do
      if (present(offset_draws)) then
         out = bsts_ar_draws(y, mean0, covariance0, ar0, observation0, &
            innovation0, observation_shape, observation_rate, innovation_shape, &
            innovation_rate, discarded, state_normals, ar_normals, gammas, &
            offset_draws)
      else
         out = bsts_ar_draws(y, mean0, covariance0, ar0, observation0, &
            innovation0, observation_shape, observation_rate, innovation_shape, &
            innovation_rate, discarded, state_normals, ar_normals, gammas)
      end if
   end function bsts_ar

   pure function bsts_auto_ar_draws(y, initial_mean, initial_covariance, &
      initial_ar_coefficients, initial_inclusion, observation_variance, &
      innovation_variance, slab_mean, slab_covariance, &
      prior_inclusion_probability, observation_prior_shape, &
      observation_prior_rate, innovation_prior_shape, innovation_prior_rate, &
      maximum_flips, truncate_stationary, burn, state_normal_draws, &
      coefficient_normal_draws, inclusion_uniform_draws, gamma_draws, &
      offset_draws) result(out)
      !! Sample a spike-and-slab AR state model from supplied random draws.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: initial_ar_coefficients(:) !! Initial autoregressive coefficients.
      logical, intent(in) :: initial_inclusion(:) !! Initial inclusion.
      real(dp), intent(in) :: observation_variance !! Observation-error variance.
      real(dp), intent(in) :: innovation_variance !! Innovation variance.
      real(dp), intent(in) :: slab_mean(:) !! Slab mean.
      real(dp), intent(in) :: slab_covariance(:, :) !! Slab covariance matrix.
      real(dp), intent(in) :: prior_inclusion_probability(:) !! Prior inclusion probability.
      real(dp), intent(in) :: observation_prior_shape !! Observation prior shape.
      real(dp), intent(in) :: observation_prior_rate !! Observation prior rate.
      real(dp), intent(in) :: innovation_prior_shape !! Innovation prior shape.
      real(dp), intent(in) :: innovation_prior_rate !! Innovation prior rate.
      integer, intent(in) :: maximum_flips !! Maximum flips.
      integer, intent(in) :: burn !! Number of initial simulation draws to discard.
      logical, intent(in) :: truncate_stationary !! Flag controlling truncate stationary.
      real(dp), intent(in) :: state_normal_draws(:, :, :) !! Independent standard-normal state draws.
      real(dp), intent(in) :: coefficient_normal_draws(:, :, :) !! Coefficient normal simulation draws.
      real(dp), intent(in) :: inclusion_uniform_draws(:, :) !! Inclusion uniform simulation draws.
      real(dp), intent(in) :: gamma_draws(:, :) !! Gamma simulation draws.
      real(dp), intent(in), optional :: offset_draws(:, :) !! Offset simulation draws.
      type(bsts_mcmc_t) :: out
      type(ssm_model_t) :: model
      real(dp), allocatable :: path(:, :), working_y(:), design(:, :), target(:)
      real(dp), allocatable :: posterior_mean(:), posterior_covariance(:, :)
      real(dp), allocatable :: active_draw(:)
      integer, allocatable :: active(:)
      real(dp) :: ar(size(initial_ar_coefficients)), proposed(size(initial_ar_coefficients))
      real(dp) :: observation_scale, innovation_scale, observation_sum
      real(dp) :: log_included, log_excluded, log_probability, collapsed_rate
      real(dp) :: probability
      logical, allocatable :: included(:), candidate_model(:)
      integer :: order, observations, iterations, attempts, iteration, time
      integer :: lag, proposal, flips, status, retained
      logical :: accepted

      order = size(initial_ar_coefficients)
      observations = size(y)
      iterations = size(state_normal_draws, 3)
      attempts = size(coefficient_normal_draws, 2)
      if (observations < 2 .or. count(ieee_is_finite(y)) < 1 .or. order < 1 .or. &
         size(initial_mean) /= order .or. &
         any(shape(initial_covariance) /= [order, order]) .or. &
         size(initial_inclusion) /= order .or. &
         size(slab_mean) /= order .or. &
         any(shape(slab_covariance) /= [order, order]) .or. &
         size(prior_inclusion_probability) /= order .or. &
         any(prior_inclusion_probability < 0.0_dp) .or. &
         any(prior_inclusion_probability > 1.0_dp) .or. &
         observation_variance <= 0.0_dp .or. innovation_variance <= 0.0_dp .or. &
         observation_prior_shape <= 0.0_dp .or. observation_prior_rate <= 0.0_dp .or. &
         innovation_prior_shape <= 0.0_dp .or. innovation_prior_rate <= 0.0_dp .or. &
         maximum_flips < -1 .or. iterations < 1 .or. attempts < 1 .or. &
         burn < 0 .or. burn >= iterations .or. &
         any(shape(state_normal_draws) /= [order, observations, iterations]) .or. &
         any(shape(coefficient_normal_draws) /= [order, attempts, iterations]) .or. &
         any(shape(inclusion_uniform_draws) /= [order, iterations]) .or. &
         any(shape(gamma_draws) /= [2, iterations]) .or. &
         any(inclusion_uniform_draws <= 0.0_dp) .or. &
         any(inclusion_uniform_draws >= 1.0_dp) .or. any(gamma_draws <= 0.0_dp) .or. &
         .not. all(ieee_is_finite(state_normal_draws)) .or. &
         .not. all(ieee_is_finite(coefficient_normal_draws))) then
         out%info = 1
         return
      end if
      if (truncate_stationary .and. &
         .not. ar_coefficients_stationary(initial_ar_coefficients)) then
         out%info = 2
         return
      end if
      if (present(offset_draws)) then
         if (any(shape(offset_draws) /= [observations, iterations]) .or. &
            .not. all(ieee_is_finite(offset_draws))) then
            out%info = 1
            return
         end if
      end if
      allocate(model%y(observations, 1), model%z(1, order, 1), &
         model%h(1, 1, 1), model%transition(order, order, 1), &
         model%r(order, 1, 1), model%q(1, 1, 1), model%a1(order), &
         model%p1(order, order), model%p1inf(order, order), &
         model%missing(observations, 1), working_y(observations), &
         design(observations - 1, order), target(observations - 1), &
         included(order), candidate_model(order), &
         out%state(order, observations, iterations), &
         out%observation_variance(iterations), &
         out%state_variance(order, iterations), &
         out%component_variance(1, iterations), &
         out%ar_coefficients(order, iterations), &
         out%ar_included(order, iterations), &
         out%ar_inclusion_probability(order), &
         out%transition(order, order), out%observation(order), &
         out%transition_schedule(order, order, 1), &
         out%state_loading_schedule(order, 1, 1))
      model%z = 0.0_dp
      model%z(1, 1, 1) = 1.0_dp
      model%r = 0.0_dp
      model%r(1, 1, 1) = 1.0_dp
      model%a1 = initial_mean
      model%p1 = initial_covariance
      model%p1inf = 0.0_dp
      model%missing(:, 1) = .not. ieee_is_finite(y)
      out%observation = 0.0_dp
      out%observation(1) = 1.0_dp
      out%state_loading_schedule = 0.0_dp
      out%state_loading_schedule(1, 1, 1) = 1.0_dp
      out%state_variance = 0.0_dp
      out%burn = burn
      out%ar_order = order
      ar = initial_ar_coefficients
      included = initial_inclusion
      where (.not. included) ar = 0.0_dp
      where (prior_inclusion_probability <= 0.0_dp) included = .false.
      where (prior_inclusion_probability >= 1.0_dp) included = .true.
      observation_scale = observation_variance
      innovation_scale = innovation_variance
      do iteration = 1, iterations
         working_y = y
         if (present(offset_draws)) working_y = working_y - &
            offset_draws(:, iteration)
         model%y(:, 1) = working_y
         where (.not. ieee_is_finite(model%y)) model%y = 0.0_dp
         model%h(1, 1, 1) = observation_scale
         call set_ar_transition(model%transition(:, :, 1), ar)
         model%q(1, 1, 1) = innovation_scale
         call state_path_from_draws(model, state_normal_draws(:, :, iteration), &
            path, status)
         if (status /= 0) then
            out%info = 10 + status
            return
         end if
         out%state(:, :, iteration) = path
         observation_sum = 0.0_dp
         do time = 1, observations
            if (ieee_is_finite(y(time))) observation_sum = observation_sum + &
               (working_y(time) - path(1, time))**2
         end do
         observation_scale = (observation_prior_rate + &
            0.5_dp*observation_sum)/gamma_draws(1, iteration)
         target = path(1, 2:)
         design = transpose(path(:, :observations - 1))
         flips = 0
         do lag = 1, order
            if (prior_inclusion_probability(lag) <= 0.0_dp) then
               included(lag) = .false.
               cycle
            end if
            if (prior_inclusion_probability(lag) >= 1.0_dp) then
               included(lag) = .true.
               cycle
            end if
            candidate_model = included
            candidate_model(lag) = .false.
            call spike_slab_model_moments(target, design, candidate_model, &
               slab_mean, slab_covariance, prior_inclusion_probability, &
               innovation_scale, innovation_prior_rate, log_excluded, &
               collapsed_rate, posterior_mean, posterior_covariance, active, &
               status)
            if (status /= 0) then
               out%info = 50 + status
               return
            end if
            candidate_model(lag) = .true.
            call spike_slab_model_moments(target, design, candidate_model, &
               slab_mean, slab_covariance, prior_inclusion_probability, &
               innovation_scale, innovation_prior_rate, log_included, &
               collapsed_rate, posterior_mean, posterior_covariance, active, &
               status)
            if (status /= 0) then
               out%info = 60 + status
               return
            end if
            probability = logistic_log_odds(log_included - log_excluded)
            candidate_model(lag) = inclusion_uniform_draws(lag, iteration) < &
               probability
            if (candidate_model(lag) .neqv. included(lag)) then
               if (maximum_flips < 0 .or. flips < maximum_flips) then
                  included(lag) = candidate_model(lag)
                  flips = flips + 1
               end if
            end if
         end do
         call spike_slab_model_moments(target, design, included, slab_mean, &
            slab_covariance, prior_inclusion_probability, innovation_scale, &
            innovation_prior_rate, log_probability, collapsed_rate, &
            posterior_mean, posterior_covariance, active, status)
         if (status /= 0) then
            out%info = 70 + status
            return
         end if
         proposed = 0.0_dp
         accepted = size(active) == 0
         if (size(active) > 0) then
            allocate(active_draw(size(active)))
            do proposal = 1, attempts
               call multivariate_normal_from_standard(posterior_mean, &
                  innovation_scale*posterior_covariance, &
                  coefficient_normal_draws(active, proposal, iteration), &
                  active_draw, status)
               if (status /= 0) then
                  out%info = 80 + status
                  return
               end if
               proposed = 0.0_dp
               proposed(active) = active_draw
               if (.not. truncate_stationary .or. &
                  ar_coefficients_stationary(proposed)) then
                  accepted = .true.
                  exit
               end if
            end do
            deallocate(active_draw)
         end if
         if (.not. accepted) then
            out%info = 3
            return
         end if
         ar = proposed
         innovation_scale = collapsed_rate/gamma_draws(2, iteration)
         out%observation_variance(iteration) = observation_scale
         out%state_variance(1, iteration) = innovation_scale
         out%component_variance(1, iteration) = innovation_scale
         out%ar_coefficients(:, iteration) = ar
         out%ar_included(:, iteration) = included
      end do
      retained = iterations - burn
      do lag = 1, order
         out%ar_inclusion_probability(lag) = real(count( &
            out%ar_included(lag, burn + 1:)), dp)/real(retained, dp)
      end do
      call set_ar_transition(out%transition, ar)
      out%transition_schedule(:, :, 1) = out%transition
   end function bsts_auto_ar_draws

   function bsts_auto_ar(y, lags, iterations, burn, maximum_flips, &
      truncate_stationary, proposal_attempts, prior_inclusion_probability, &
      slab_mean, slab_standard_deviation, offset_draws) result(out)
      !! Sample spike-and-slab AR lag selection using the random stream.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      integer, intent(in) :: lags !! Lags.
      integer, intent(in) :: iterations !! Number of algorithm iterations.
      integer, intent(in), optional :: burn !! Number of initial simulation draws to discard.
      integer, intent(in), optional :: maximum_flips !! Maximum flips.
      integer, intent(in), optional :: proposal_attempts !! Proposal attempts.
      logical, intent(in), optional :: truncate_stationary !! Flag controlling truncate stationary.
      real(dp), intent(in), optional :: prior_inclusion_probability(:) !! Prior inclusion probability.
      real(dp), intent(in), optional :: slab_mean(:) !! Slab mean.
      real(dp), intent(in), optional :: slab_standard_deviation(:) !! Slab standard deviation.
      real(dp), intent(in), optional :: offset_draws(:, :) !! Offset simulation draws.
      type(bsts_mcmc_t) :: out
      real(dp), allocatable :: mean0(:), covariance0(:, :), ar0(:)
      real(dp), allocatable :: inclusion0(:), slab0(:), slab_sd0(:)
      real(dp), allocatable :: slab_covariance(:, :)
      real(dp), allocatable :: state_normals(:, :, :), coefficient_normals(:, :, :)
      real(dp), allocatable :: uniforms(:, :), gammas(:, :)
      logical, allocatable :: model0(:)
      real(dp) :: series_variance, observation0, innovation0
      integer :: discarded, flip_limit, attempts, iteration, proposal, lag
      integer :: observed
      logical :: truncate

      if (size(y) < 2 .or. count(ieee_is_finite(y)) < 1 .or. lags < 1 .or. &
         iterations < 1) then
         out%info = 1
         return
      end if
      discarded = 0
      if (present(burn)) discarded = burn
      flip_limit = -1
      if (present(maximum_flips)) flip_limit = maximum_flips
      attempts = 100
      if (present(proposal_attempts)) attempts = proposal_attempts
      truncate = .true.
      if (present(truncate_stationary)) truncate = truncate_stationary
      if (discarded < 0 .or. discarded >= iterations .or. flip_limit < -1 .or. &
         attempts < 1) then
         out%info = 1
         return
      end if
      series_variance = max(finite_sample_variance(y), 1.0e-6_dp)
      observation0 = max(0.1_dp*series_variance, 1.0e-8_dp)
      innovation0 = max(0.01_dp*series_variance, 1.0e-10_dp)
      allocate(mean0(lags), covariance0(lags, lags), ar0(lags), &
         inclusion0(lags), slab0(lags), slab_sd0(lags), &
         slab_covariance(lags, lags), model0(lags), &
         state_normals(lags, size(y), iterations), &
         coefficient_normals(lags, attempts, iterations), &
         uniforms(lags, iterations), gammas(2, iterations))
      mean0 = 0.0_dp
      mean0(1) = first_finite_value(y)
      covariance0 = identity_block(lags, series_variance)
      ar0 = 0.0_dp
      model0 = .false.
      do lag = 1, lags
         inclusion0(lag) = 0.8_dp*0.8_dp**real(lag - 1, dp)
         slab_sd0(lag) = 0.5_dp*0.8_dp**real(lag - 1, dp)
      end do
      slab0 = 0.0_dp
      if (present(prior_inclusion_probability)) then
         if (size(prior_inclusion_probability) /= lags) then
            out%info = 1
            return
         end if
         inclusion0 = prior_inclusion_probability
      end if
      if (present(slab_mean)) then
         if (size(slab_mean) /= lags) then
            out%info = 1
            return
         end if
         slab0 = slab_mean
      end if
      if (present(slab_standard_deviation)) then
         if (size(slab_standard_deviation) /= lags) then
            out%info = 1
            return
         end if
         slab_sd0 = slab_standard_deviation
      end if
      if (any(inclusion0 < 0.0_dp) .or. any(inclusion0 > 1.0_dp) .or. &
         any(slab_sd0 <= 0.0_dp)) then
         out%info = 1
         return
      end if
      slab_covariance = 0.0_dp
      do lag = 1, lags
         slab_covariance(lag, lag) = slab_sd0(lag)**2
      end do
      observed = count(ieee_is_finite(y))
      do iteration = 1, iterations
         call fill_standard_normals(state_normals(:, :, iteration))
         do proposal = 1, attempts
            call fill_vector_normals(coefficient_normals(:, proposal, iteration))
         end do
         do lag = 1, lags
            uniforms(lag, iteration) = max(tiny(1.0_dp), &
               min(1.0_dp - epsilon(1.0_dp), random_uniform()))
         end do
         gammas(1, iteration) = random_gamma(2.0_dp + &
            0.5_dp*real(observed, dp))
         gammas(2, iteration) = random_gamma(2.0_dp + &
            0.5_dp*real(size(y) - 1, dp))
      end do
      if (present(offset_draws)) then
         out = bsts_auto_ar_draws(y, mean0, covariance0, ar0, model0, &
            observation0, innovation0, slab0, slab_covariance, inclusion0, &
            2.0_dp, observation0, 2.0_dp, innovation0, flip_limit, truncate, &
            discarded, state_normals, coefficient_normals, uniforms, gammas, &
            offset_draws)
      else
         out = bsts_auto_ar_draws(y, mean0, covariance0, ar0, model0, &
            observation0, innovation0, slab0, slab_covariance, inclusion0, &
            2.0_dp, observation0, 2.0_dp, innovation0, flip_limit, truncate, &
            discarded, state_normals, coefficient_normals, uniforms, gammas)
      end if
   end function bsts_auto_ar

   pure function bsts_ar_predict_draws(fit, horizon, state_normal_draws, &
      observation_normal_draws) result(out)
      !! Forecast a Bayesian AR state component from supplied normal draws.
      type(bsts_mcmc_t), intent(in) :: fit !! Previously fitted model.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      real(dp), intent(in) :: state_normal_draws(:, :) !! Independent standard-normal state draws.
      real(dp), intent(in) :: observation_normal_draws(:, :) !! Observation normal draws.
      type(bsts_prediction_t) :: out
      real(dp), allocatable :: state(:), sorted(:)
      real(dp) :: next_state
      integer :: order, retained, draw, source, step

      if (fit%info /= 0 .or. fit%ar_order < 1 .or. &
         .not. allocated(fit%state) .or. .not. allocated(fit%ar_coefficients) .or. &
         .not. allocated(fit%component_variance) .or. &
         .not. allocated(fit%observation_variance)) then
         out%info = 1
         return
      end if
      order = fit%ar_order
      retained = size(fit%state, 3) - fit%burn
      if (horizon < 1 .or. retained < 1 .or. &
         any(shape(state_normal_draws) /= [horizon, retained]) .or. &
         any(shape(observation_normal_draws) /= [horizon, retained]) .or. &
         .not. all(ieee_is_finite(state_normal_draws)) .or. &
         .not. all(ieee_is_finite(observation_normal_draws))) then
         out%info = 1
         return
      end if
      allocate(out%draws(horizon, retained), out%mean(horizon), &
         out%standard_deviation(horizon), out%lower(horizon), &
         out%upper(horizon), state(order))
      do draw = 1, retained
         source = fit%burn + draw
         state = fit%state(:, size(fit%state, 2), source)
         do step = 1, horizon
            next_state = dot_product(fit%ar_coefficients(:, source), state) + &
               sqrt(fit%component_variance(1, source))* &
               state_normal_draws(step, draw)
            if (order > 1) state(2:) = state(:order - 1)
            state(1) = next_state
            out%draws(step, draw) = state(1) + &
               sqrt(fit%observation_variance(source))* &
               observation_normal_draws(step, draw)
         end do
      end do
      do step = 1, horizon
         out%mean(step) = sum(out%draws(step, :))/real(retained, dp)
         if (retained > 1) then
            out%standard_deviation(step) = sqrt(sum((out%draws(step, :) - &
               out%mean(step))**2)/real(retained - 1, dp))
         else
            out%standard_deviation(step) = 0.0_dp
         end if
         sorted = out%draws(step, :)
         call insertion_sort(sorted)
         out%lower(step) = quantile(sorted, 0.025_dp)
         out%upper(step) = quantile(sorted, 0.975_dp)
      end do
   end function bsts_ar_predict_draws

   function bsts_ar_predict(fit, horizon) result(out)
      !! Forecast a Bayesian AR state component using the random stream.
      type(bsts_mcmc_t), intent(in) :: fit !! Previously fitted model.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      type(bsts_prediction_t) :: out
      real(dp), allocatable :: state_normals(:, :), observation_normals(:, :)
      integer :: retained, draw, step

      if (.not. allocated(fit%state)) then
         out%info = 1
         return
      end if
      retained = size(fit%state, 3) - fit%burn
      if (horizon < 1 .or. retained < 1) then
         out%info = 1
         return
      end if
      allocate(state_normals(horizon, retained), &
         observation_normals(horizon, retained))
      do draw = 1, retained
         do step = 1, horizon
            state_normals(step, draw) = random_standard_normal()
            observation_normals(step, draw) = random_standard_normal()
         end do
      end do
      out = bsts_ar_predict_draws(fit, horizon, state_normals, &
         observation_normals)
   end function bsts_ar_predict

   pure function bsts_dynamic_regression_draws(response, predictors, &
      initial_mean, initial_covariance, residual_variance, &
      innovation_variance, residual_prior_shape, residual_prior_rate, &
      innovation_prior_shape, innovation_prior_rate, burn, &
      state_normal_draws, gamma_draws, offset_draws) result(out)
      !! Sample independent random-walk coefficient variances from supplied draws.
      real(dp), intent(in) :: response(:) !! Response observations.
      real(dp), intent(in) :: predictors(:, :) !! Predictor matrix.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: residual_variance !! Residual variance.
      real(dp), intent(in) :: innovation_variance(:) !! Innovation variance.
      real(dp), intent(in) :: residual_prior_shape !! Residual prior shape.
      real(dp), intent(in) :: residual_prior_rate !! Residual prior rate.
      real(dp), intent(in) :: innovation_prior_shape(:) !! Innovation prior shape.
      real(dp), intent(in) :: innovation_prior_rate(:) !! Innovation prior rate.
      integer, intent(in) :: burn !! Number of initial simulation draws to discard.
      real(dp), intent(in) :: state_normal_draws(:, :, :) !! Independent standard-normal state draws.
      real(dp), intent(in) :: gamma_draws(:, :) !! Gamma simulation draws.
      real(dp), intent(in), optional :: offset_draws(:, :) !! Offset simulation draws.
      type(bsts_dynamic_regression_t) :: out
      real(dp) :: unused_variance(size(innovation_variance))

      unused_variance = 1.0_dp
      if (present(offset_draws)) then
         out = dynamic_regression_draws_core(response, predictors, initial_mean, &
            initial_covariance, residual_variance, innovation_variance, &
            residual_prior_shape, residual_prior_rate, innovation_prior_shape, &
            innovation_prior_rate, burn, state_normal_draws, gamma_draws, &
            .false., unused_variance, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, &
            offset_draws)
      else
         out = dynamic_regression_draws_core(response, predictors, initial_mean, &
            initial_covariance, residual_variance, innovation_variance, &
            residual_prior_shape, residual_prior_rate, innovation_prior_shape, &
            innovation_prior_rate, burn, state_normal_draws, gamma_draws, &
            .false., unused_variance, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp)
      end if
   end function bsts_dynamic_regression_draws

   pure function bsts_dynamic_regression_hierarchical_draws(response, &
      predictors, initial_mean, initial_covariance, residual_variance, &
      scaled_innovation_variance, predictor_variance, residual_prior_shape, &
      residual_prior_rate, shrinkage_shape, initial_hierarchy_rate, &
      hierarchy_rate_prior_shape, hierarchy_rate_prior_rate, burn, &
      state_normal_draws, gamma_draws, offset_draws) result(out)
      !! Sample a shared-rate hierarchy for random-walk coefficient variances.
      real(dp), intent(in) :: response(:) !! Response observations.
      real(dp), intent(in) :: predictors(:, :) !! Predictor matrix.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: residual_variance !! Residual variance.
      real(dp), intent(in) :: scaled_innovation_variance(:) !! Scaled innovation variance.
      real(dp), intent(in) :: predictor_variance(:) !! Predictor variance.
      real(dp), intent(in) :: residual_prior_shape !! Residual prior shape.
      real(dp), intent(in) :: residual_prior_rate !! Residual prior rate.
      real(dp), intent(in) :: shrinkage_shape !! Shrinkage shape.
      real(dp), intent(in) :: initial_hierarchy_rate !! Initial hierarchy rate.
      real(dp), intent(in) :: hierarchy_rate_prior_shape !! Hierarchy rate prior shape.
      real(dp), intent(in) :: hierarchy_rate_prior_rate !! Hierarchy rate prior rate.
      integer, intent(in) :: burn !! Number of initial simulation draws to discard.
      real(dp), intent(in) :: state_normal_draws(:, :, :) !! Independent standard-normal state draws.
      real(dp), intent(in) :: gamma_draws(:, :) !! Gamma simulation draws.
      real(dp), intent(in), optional :: offset_draws(:, :) !! Offset simulation draws.
      type(bsts_dynamic_regression_t) :: out
      real(dp) :: actual_variance(size(scaled_innovation_variance))
      real(dp) :: unused_shape(size(scaled_innovation_variance))
      real(dp) :: unused_rate(size(scaled_innovation_variance))

      actual_variance = scaled_innovation_variance/predictor_variance
      unused_shape = shrinkage_shape
      unused_rate = initial_hierarchy_rate
      if (present(offset_draws)) then
         out = dynamic_regression_draws_core(response, predictors, initial_mean, &
            initial_covariance, residual_variance, actual_variance, &
            residual_prior_shape, residual_prior_rate, unused_shape, &
            unused_rate, burn, state_normal_draws, gamma_draws, .true., &
            predictor_variance, shrinkage_shape, initial_hierarchy_rate, &
            hierarchy_rate_prior_shape, hierarchy_rate_prior_rate, offset_draws)
      else
         out = dynamic_regression_draws_core(response, predictors, initial_mean, &
            initial_covariance, residual_variance, actual_variance, &
            residual_prior_shape, residual_prior_rate, unused_shape, &
            unused_rate, burn, state_normal_draws, gamma_draws, .true., &
            predictor_variance, shrinkage_shape, initial_hierarchy_rate, &
            hierarchy_rate_prior_shape, hierarchy_rate_prior_rate)
      end if
   end function bsts_dynamic_regression_hierarchical_draws

   function bsts_dynamic_regression(response, predictors, iterations, burn, &
      initial_mean, initial_covariance, residual_variance, &
      innovation_variance, residual_prior_shape, residual_prior_rate, &
      innovation_prior_shape, innovation_prior_rate, offset_draws) result(out)
      !! Sample independent dynamic-regression variances using the random stream.
      real(dp), intent(in) :: response(:) !! Response observations.
      real(dp), intent(in) :: predictors(:, :) !! Predictor matrix.
      integer, intent(in) :: iterations !! Number of algorithm iterations.
      integer, intent(in), optional :: burn !! Number of initial simulation draws to discard.
      real(dp), intent(in), optional :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in), optional :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in), optional :: residual_variance !! Residual variance.
      real(dp), intent(in), optional :: innovation_variance(:) !! Innovation variance.
      real(dp), intent(in), optional :: residual_prior_shape !! Residual prior shape.
      real(dp), intent(in), optional :: residual_prior_rate !! Residual prior rate.
      real(dp), intent(in), optional :: innovation_prior_shape(:) !! Innovation prior shape.
      real(dp), intent(in), optional :: innovation_prior_rate(:) !! Innovation prior rate.
      real(dp), intent(in), optional :: offset_draws(:, :) !! Offset simulation draws.
      type(bsts_dynamic_regression_t) :: out
      real(dp), allocatable :: mean0(:), covariance0(:, :), variance0(:)
      real(dp), allocatable :: shape0(:), rate0(:), normals(:, :, :)
      real(dp), allocatable :: gammas(:, :), predictor_scale(:)
      real(dp) :: response_scale, observation0, observation_shape
      real(dp) :: observation_rate
      integer :: variables, discarded, iteration, component, observed

      variables = size(predictors, 2)
      if (.not. dynamic_regression_data_valid(response, predictors) .or. &
         iterations < 1) then
         out%info = 1
         return
      end if
      discarded = 0
      if (present(burn)) discarded = burn
      if (discarded < 0 .or. discarded >= iterations) then
         out%info = 1
         return
      end if
      response_scale = max(finite_sample_variance(response), 1.0e-6_dp)
      predictor_scale = predictor_variances(predictors)
      allocate(mean0(variables), covariance0(variables, variables), &
         variance0(variables), shape0(variables), rate0(variables))
      mean0 = 0.0_dp
      if (present(initial_mean)) then
         if (size(initial_mean) /= variables) then
            out%info = 1
            return
         end if
         mean0 = initial_mean
      end if
      covariance0 = 0.0_dp
      do component = 1, variables
         covariance0(component, component) = 10.0_dp*response_scale/ &
            predictor_scale(component)
      end do
      if (present(initial_covariance)) then
         if (any(shape(initial_covariance) /= [variables, variables])) then
            out%info = 1
            return
         end if
         covariance0 = initial_covariance
      end if
      observation0 = max(0.1_dp*response_scale, 1.0e-8_dp)
      if (present(residual_variance)) observation0 = residual_variance
      observation_shape = 2.0_dp
      if (present(residual_prior_shape)) observation_shape = residual_prior_shape
      observation_rate = observation0
      if (present(residual_prior_rate)) observation_rate = residual_prior_rate
      variance0 = 1.0e-4_dp*response_scale/predictor_scale
      if (present(innovation_variance)) then
         if (size(innovation_variance) /= variables) then
            out%info = 1
            return
         end if
         variance0 = innovation_variance
      end if
      shape0 = 2.0_dp
      if (present(innovation_prior_shape)) then
         if (size(innovation_prior_shape) /= variables) then
            out%info = 1
            return
         end if
         shape0 = innovation_prior_shape
      end if
      rate0 = variance0
      if (present(innovation_prior_rate)) then
         if (size(innovation_prior_rate) /= variables) then
            out%info = 1
            return
         end if
         rate0 = innovation_prior_rate
      end if
      if (observation0 <= 0.0_dp .or. observation_shape <= 0.0_dp .or. &
         observation_rate <= 0.0_dp .or. any(variance0 <= 0.0_dp) .or. &
         any(shape0 <= 0.0_dp) .or. any(rate0 <= 0.0_dp)) then
         out%info = 1
         return
      end if
      allocate(normals(variables, size(response), iterations), &
         gammas(variables + 1, iterations))
      observed = count(ieee_is_finite(response))
      do iteration = 1, iterations
         call fill_standard_normals(normals(:, :, iteration))
         gammas(1, iteration) = random_gamma(observation_shape + &
            0.5_dp*real(observed, dp))
         do component = 1, variables
            gammas(component + 1, iteration) = random_gamma(shape0(component) + &
               0.5_dp*real(size(response) - 1, dp))
         end do
      end do
      if (present(offset_draws)) then
         out = bsts_dynamic_regression_draws(response, predictors, mean0, &
            covariance0, observation0, variance0, observation_shape, &
            observation_rate, shape0, rate0, discarded, normals, gammas, &
            offset_draws)
      else
         out = bsts_dynamic_regression_draws(response, predictors, mean0, &
            covariance0, observation0, variance0, observation_shape, &
            observation_rate, shape0, rate0, discarded, normals, gammas)
      end if
   end function bsts_dynamic_regression

   function bsts_dynamic_regression_hierarchical(response, predictors, &
      iterations, burn, shrinkage_shape, hierarchy_rate_prior_shape, &
      hierarchy_rate_prior_rate, offset_draws) result(out)
      !! Sample scale-adjusted dynamic coefficients with a shared prior rate.
      real(dp), intent(in) :: response(:) !! Response observations.
      real(dp), intent(in) :: predictors(:, :) !! Predictor matrix.
      integer, intent(in) :: iterations !! Number of algorithm iterations.
      integer, intent(in), optional :: burn !! Number of initial simulation draws to discard.
      real(dp), intent(in), optional :: shrinkage_shape !! Shrinkage shape.
      real(dp), intent(in), optional :: hierarchy_rate_prior_shape !! Hierarchy rate prior shape.
      real(dp), intent(in), optional :: hierarchy_rate_prior_rate !! Hierarchy rate prior rate.
      real(dp), intent(in), optional :: offset_draws(:, :) !! Offset simulation draws.
      type(bsts_dynamic_regression_t) :: out
      real(dp), allocatable :: mean0(:), covariance0(:, :), scaled0(:)
      real(dp), allocatable :: predictor_scale(:), normals(:, :, :)
      real(dp), allocatable :: gammas(:, :)
      real(dp) :: response_scale, observation0, shape0, rate0
      real(dp) :: shrinkage, hierarchy0, hierarchy_shape, hierarchy_prior_rate
      integer :: variables, discarded, iteration, component, observed

      variables = size(predictors, 2)
      if (.not. dynamic_regression_data_valid(response, predictors) .or. &
         iterations < 1) then
         out%info = 1
         return
      end if
      discarded = 0
      if (present(burn)) discarded = burn
      shrinkage = 10.0_dp
      if (present(shrinkage_shape)) shrinkage = shrinkage_shape
      hierarchy_shape = 1.0_dp
      if (present(hierarchy_rate_prior_shape)) &
         hierarchy_shape = hierarchy_rate_prior_shape
      hierarchy_prior_rate = 1.0_dp
      if (present(hierarchy_rate_prior_rate)) &
         hierarchy_prior_rate = hierarchy_rate_prior_rate
      if (discarded < 0 .or. discarded >= iterations .or. shrinkage <= 0.0_dp .or. &
         hierarchy_shape <= 0.0_dp .or. hierarchy_prior_rate <= 0.0_dp) then
         out%info = 1
         return
      end if
      response_scale = max(finite_sample_variance(response), 1.0e-6_dp)
      predictor_scale = predictor_variances(predictors)
      observation0 = max(0.1_dp*response_scale, 1.0e-8_dp)
      shape0 = 2.0_dp
      rate0 = observation0
      hierarchy0 = max((shrinkage - 1.0_dp)*1.0e-4_dp*response_scale, &
         1.0e-10_dp)
      allocate(mean0(variables), covariance0(variables, variables), &
         scaled0(variables), normals(variables, size(response), iterations), &
         gammas(variables + 2, iterations))
      mean0 = 0.0_dp
      covariance0 = 0.0_dp
      scaled0 = max(hierarchy0/max(shrinkage, 1.0_dp), 1.0e-10_dp)
      do component = 1, variables
         covariance0(component, component) = 10.0_dp*response_scale/ &
            predictor_scale(component)
      end do
      observed = count(ieee_is_finite(response))
      do iteration = 1, iterations
         call fill_standard_normals(normals(:, :, iteration))
         gammas(1, iteration) = random_gamma(shape0 + &
            0.5_dp*real(observed, dp))
         do component = 1, variables
            gammas(component + 1, iteration) = random_gamma(shrinkage + &
               0.5_dp*real(size(response) - 1, dp))
         end do
         gammas(variables + 2, iteration) = random_gamma(hierarchy_shape + &
            real(variables, dp)*shrinkage)
      end do
      if (present(offset_draws)) then
         out = bsts_dynamic_regression_hierarchical_draws(response, predictors, &
            mean0, covariance0, observation0, scaled0, predictor_scale, shape0, &
            rate0, shrinkage, hierarchy0, hierarchy_shape, &
            hierarchy_prior_rate, discarded, normals, gammas, offset_draws)
      else
         out = bsts_dynamic_regression_hierarchical_draws(response, predictors, &
            mean0, covariance0, observation0, scaled0, predictor_scale, shape0, &
            rate0, shrinkage, hierarchy0, hierarchy_shape, &
            hierarchy_prior_rate, discarded, normals, gammas)
      end if
   end function bsts_dynamic_regression_hierarchical

   pure function bsts_dynamic_regression_predict_draws(fit, future_predictors, &
      state_normal_draws, observation_normal_draws) result(out)
      !! Forecast dynamic regression from supplied state and observation normals.
      type(bsts_dynamic_regression_t), intent(in) :: fit !! Previously fitted model.
      real(dp), intent(in) :: future_predictors(:, :) !! Future predictors.
      real(dp), intent(in) :: state_normal_draws(:, :, :) !! Independent standard-normal state draws.
      real(dp), intent(in) :: observation_normal_draws(:, :) !! Observation normal draws.
      type(bsts_prediction_t) :: out
      real(dp), allocatable :: coefficient(:), sorted(:)
      integer :: horizon, retained, variables, draw, source, step

      if (fit%info /= 0 .or. .not. allocated(fit%coefficients) .or. &
         .not. allocated(fit%innovation_variance) .or. &
         .not. allocated(fit%residual_variance)) then
         out%info = 1
         return
      end if
      variables = size(fit%coefficients, 1)
      horizon = size(future_predictors, 1)
      retained = size(fit%coefficients, 3) - fit%burn
      if (horizon < 1 .or. retained < 1 .or. &
         size(future_predictors, 2) /= variables .or. &
         any(shape(state_normal_draws) /= [variables, horizon, retained]) .or. &
         any(shape(observation_normal_draws) /= [horizon, retained]) .or. &
         .not. all(ieee_is_finite(future_predictors)) .or. &
         .not. all(ieee_is_finite(state_normal_draws)) .or. &
         .not. all(ieee_is_finite(observation_normal_draws))) then
         out%info = 1
         return
      end if
      allocate(out%draws(horizon, retained), out%mean(horizon), &
         out%standard_deviation(horizon), out%lower(horizon), &
         out%upper(horizon), coefficient(variables))
      do draw = 1, retained
         source = fit%burn + draw
         coefficient = fit%coefficients(:, size(fit%coefficients, 2), source)
         do step = 1, horizon
            coefficient = coefficient + sqrt(fit%innovation_variance(:, source))* &
               state_normal_draws(:, step, draw)
            out%draws(step, draw) = dot_product(future_predictors(step, :), &
               coefficient) + sqrt(fit%residual_variance(source))* &
               observation_normal_draws(step, draw)
         end do
      end do
      do step = 1, horizon
         out%mean(step) = sum(out%draws(step, :))/real(retained, dp)
         if (retained > 1) then
            out%standard_deviation(step) = sqrt(sum((out%draws(step, :) - &
               out%mean(step))**2)/real(retained - 1, dp))
         else
            out%standard_deviation(step) = 0.0_dp
         end if
         sorted = out%draws(step, :)
         call insertion_sort(sorted)
         out%lower(step) = quantile(sorted, 0.025_dp)
         out%upper(step) = quantile(sorted, 0.975_dp)
      end do
   end function bsts_dynamic_regression_predict_draws

   function bsts_dynamic_regression_predict(fit, future_predictors) result(out)
      !! Forecast dynamic regression using the shared random stream.
      type(bsts_dynamic_regression_t), intent(in) :: fit !! Previously fitted model.
      real(dp), intent(in) :: future_predictors(:, :) !! Future predictors.
      type(bsts_prediction_t) :: out
      real(dp), allocatable :: state_normals(:, :, :), observation_normals(:, :)
      integer :: variables, horizon, retained, draw, step, component

      if (.not. allocated(fit%coefficients)) then
         out%info = 1
         return
      end if
      variables = size(fit%coefficients, 1)
      horizon = size(future_predictors, 1)
      retained = size(fit%coefficients, 3) - fit%burn
      if (horizon < 1 .or. retained < 1) then
         out%info = 1
         return
      end if
      allocate(state_normals(variables, horizon, retained), &
         observation_normals(horizon, retained))
      do draw = 1, retained
         do step = 1, horizon
            do component = 1, variables
               state_normals(component, step, draw) = random_standard_normal()
            end do
            observation_normals(step, draw) = random_standard_normal()
         end do
      end do
      out = bsts_dynamic_regression_predict_draws(fit, future_predictors, &
         state_normals, observation_normals)
   end function bsts_dynamic_regression_predict

   pure function bsts_dynamic_regression_ar_draws(response, predictors, lags, &
      initial_state_mean, initial_state_covariance, initial_ar_coefficients, &
      residual_variance, scaled_innovation_variance, predictor_mean_square, &
      residual_prior_shape, residual_prior_rate, innovation_prior_shape, &
      innovation_prior_rate, burn, state_normal_draws, ar_normal_draws, &
      gamma_draws, offset_draws) result(out)
      !! Sample AR(p) dynamic-regression paths from supplied random variates.
      real(dp), intent(in) :: response(:) !! Response observations.
      real(dp), intent(in) :: predictors(:, :) !! Predictor matrix.
      integer, intent(in) :: lags !! Lags.
      real(dp), intent(in) :: initial_state_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_state_covariance(:, :) !! Initial state covariance.
      real(dp), intent(in) :: initial_ar_coefficients(:, :) !! Initial autoregressive coefficients.
      real(dp), intent(in) :: residual_variance !! Residual variance.
      real(dp), intent(in) :: scaled_innovation_variance(:) !! Scaled innovation variance.
      real(dp), intent(in) :: predictor_mean_square(:) !! Predictor mean square.
      real(dp), intent(in) :: residual_prior_shape !! Residual prior shape.
      real(dp), intent(in) :: residual_prior_rate !! Residual prior rate.
      real(dp), intent(in) :: innovation_prior_shape(:) !! Innovation prior shape.
      real(dp), intent(in) :: innovation_prior_rate(:) !! Innovation prior rate.
      integer, intent(in) :: burn !! Number of initial simulation draws to discard.
      real(dp), intent(in) :: state_normal_draws(:, :, :) !! Independent standard-normal state draws.
      real(dp), intent(in) :: ar_normal_draws(:, :, :, :) !! Autoregressive normal simulation draws.
      real(dp), intent(in) :: gamma_draws(:, :) !! Gamma simulation draws.
      real(dp), intent(in), optional :: offset_draws(:, :) !! Offset simulation draws.
      type(bsts_dynamic_regression_t) :: out
      type(ssm_model_t) :: model
      real(dp), allocatable :: path(:, :), working_response(:)
      real(dp), allocatable :: design(:, :), target(:), crossproduct(:, :)
      real(dp), allocatable :: inverse(:, :), posterior_mean(:), candidate(:)
      real(dp) :: ar(lags, size(predictors, 2))
      real(dp) :: scaled(size(predictors, 2)), state_variance(size(predictors, 2))
      real(dp) :: observation, residual_sum, innovation_sum
      integer :: observations, variables, dimension, iterations, attempts
      integer :: iteration, component, time, lag, start, proposal, status, observed

      observations = size(response)
      variables = size(predictors, 2)
      dimension = variables*lags
      iterations = size(state_normal_draws, 3)
      attempts = size(ar_normal_draws, 3)
      if (.not. dynamic_regression_data_valid(response, predictors) .or. &
         lags < 1 .or. size(initial_state_mean) /= dimension .or. &
         any(shape(initial_state_covariance) /= [dimension, dimension]) .or. &
         any(shape(initial_ar_coefficients) /= [lags, variables]) .or. &
         residual_variance <= 0.0_dp .or. &
         size(scaled_innovation_variance) /= variables .or. &
         any(scaled_innovation_variance <= 0.0_dp) .or. &
         size(predictor_mean_square) /= variables .or. &
         any(predictor_mean_square <= 0.0_dp) .or. &
         residual_prior_shape <= 0.0_dp .or. residual_prior_rate <= 0.0_dp .or. &
         size(innovation_prior_shape) /= variables .or. &
         size(innovation_prior_rate) /= variables .or. &
         any(innovation_prior_shape <= 0.0_dp) .or. &
         any(innovation_prior_rate <= 0.0_dp) .or. &
         iterations < 1 .or. attempts < 1 .or. burn < 0 .or. burn >= iterations .or. &
         any(shape(state_normal_draws) /= [dimension, observations, iterations]) .or. &
         any(shape(ar_normal_draws) /= [lags, variables, attempts, iterations]) .or. &
         any(shape(gamma_draws) /= [variables + 1, iterations]) .or. &
         any(gamma_draws <= 0.0_dp) .or. &
         .not. all(ieee_is_finite(state_normal_draws)) .or. &
         .not. all(ieee_is_finite(ar_normal_draws))) then
         out%info = 1
         return
      end if
      do component = 1, variables
         if (.not. ar_coefficients_stationary(initial_ar_coefficients(:, component))) then
            out%info = 2
            return
         end if
      end do
      if (present(offset_draws)) then
         if (any(shape(offset_draws) /= [observations, iterations]) .or. &
            .not. all(ieee_is_finite(offset_draws))) then
            out%info = 1
            return
         end if
      end if
      allocate(model%y(observations, 1), model%z(1, dimension, observations), &
         model%h(1, 1, 1), model%transition(dimension, dimension, 1), &
         model%r(dimension, variables, 1), model%q(variables, variables, 1), &
         model%a1(dimension), model%p1(dimension, dimension), &
         model%p1inf(dimension, dimension), model%missing(observations, 1), &
         working_response(observations), design(observations - 1, lags), &
         target(observations - 1), crossproduct(lags, lags), &
         inverse(lags, lags), posterior_mean(lags), candidate(lags), &
         out%coefficients(variables, observations, iterations), &
         out%ar_coefficients(lags, variables, iterations), &
         out%innovation_variance(variables, iterations), &
         out%scaled_innovation_variance(variables, iterations), &
         out%residual_variance(iterations), out%predictor_variance(variables))
      model%z = 0.0_dp
      model%r = 0.0_dp
      do component = 1, variables
         start = (component - 1)*lags + 1
         model%z(1, start, :) = predictors(:, component)
         model%r(start, component, 1) = 1.0_dp
      end do
      model%a1 = initial_state_mean
      model%p1 = initial_state_covariance
      model%p1inf = 0.0_dp
      model%missing(:, 1) = .not. ieee_is_finite(response)
      observation = residual_variance
      scaled = scaled_innovation_variance
      state_variance = scaled/predictor_mean_square
      ar = initial_ar_coefficients
      out%burn = burn
      out%ar_order = lags
      out%predictor_variance = predictor_mean_square
      observed = count(ieee_is_finite(response))
      do iteration = 1, iterations
         working_response = response
         if (present(offset_draws)) working_response = working_response - &
            offset_draws(:, iteration)
         model%y(:, 1) = working_response
         where (.not. ieee_is_finite(model%y)) model%y = 0.0_dp
         model%h(1, 1, 1) = observation
         call set_dynamic_ar_transition(model%transition(:, :, 1), ar, lags)
         model%q = 0.0_dp
         do component = 1, variables
            model%q(component, component, 1) = state_variance(component)
         end do
         call state_path_from_draws(model, state_normal_draws(:, :, iteration), &
            path, status)
         if (status /= 0) then
            out%info = 10 + status
            return
         end if
         do component = 1, variables
            start = (component - 1)*lags + 1
            out%coefficients(component, :, iteration) = path(start, :)
         end do
         residual_sum = 0.0_dp
         do time = 1, observations
            if (ieee_is_finite(response(time))) residual_sum = residual_sum + &
               (working_response(time) - dot_product(predictors(time, :), &
               out%coefficients(:, time, iteration)))**2
         end do
         observation = (residual_prior_rate + 0.5_dp*residual_sum)/ &
            gamma_draws(1, iteration)
         do component = 1, variables
            start = (component - 1)*lags + 1
            do time = 1, observations - 1
               target(time) = path(start, time + 1)
               do lag = 1, lags
                  design(time, lag) = path(start + lag - 1, time)
               end do
            end do
            crossproduct = matmul(transpose(design), design)
            call invert_matrix(crossproduct, inverse, status)
            if (status /= 0) then
               out%info = 50 + status
               return
            end if
            posterior_mean = matmul(inverse, matmul(transpose(design), target))
            do proposal = 1, attempts
               call multivariate_normal_from_standard(posterior_mean, &
                  state_variance(component)*inverse, &
                  ar_normal_draws(:, component, proposal, iteration), &
                  candidate, status)
               if (status /= 0) then
                  out%info = 60 + status
                  return
               end if
               if (ar_coefficients_stationary(candidate)) then
                  ar(:, component) = candidate
                  exit
               end if
            end do
            innovation_sum = sum((target - matmul(design, ar(:, component)))**2)
            scaled(component) = (innovation_prior_rate(component) + &
               0.5_dp*predictor_mean_square(component)*innovation_sum)/ &
               gamma_draws(component + 1, iteration)
            state_variance(component) = scaled(component)/ &
               predictor_mean_square(component)
         end do
         out%ar_coefficients(:, :, iteration) = ar
         out%residual_variance(iteration) = observation
         out%innovation_variance(:, iteration) = state_variance
         out%scaled_innovation_variance(:, iteration) = scaled
      end do
   end function bsts_dynamic_regression_ar_draws

   function bsts_dynamic_regression_ar(response, predictors, lags, iterations, &
      burn, proposal_attempts, offset_draws) result(out)
      !! Sample AR(p) dynamic regression using the shared random stream.
      real(dp), intent(in) :: response(:) !! Response observations.
      real(dp), intent(in) :: predictors(:, :) !! Predictor matrix.
      integer, intent(in) :: lags !! Lags.
      integer, intent(in) :: iterations !! Number of algorithm iterations.
      integer, intent(in), optional :: burn !! Number of initial simulation draws to discard.
      integer, intent(in), optional :: proposal_attempts !! Proposal attempts.
      real(dp), intent(in), optional :: offset_draws(:, :) !! Offset simulation draws.
      type(bsts_dynamic_regression_t) :: out
      real(dp), allocatable :: mean0(:), covariance0(:, :), ar0(:, :)
      real(dp), allocatable :: scaled0(:), mean_square(:), shape0(:), rate0(:)
      real(dp), allocatable :: state_normals(:, :, :), ar_normals(:, :, :, :)
      real(dp), allocatable :: gammas(:, :)
      real(dp) :: response_scale, observation0, observation_shape, observation_rate
      integer :: variables, dimension, discarded, attempts, iteration
      integer :: component, proposal, observed

      variables = size(predictors, 2)
      dimension = variables*lags
      if (.not. dynamic_regression_data_valid(response, predictors) .or. &
         lags < 1 .or. iterations < 1) then
         out%info = 1
         return
      end if
      discarded = 0
      if (present(burn)) discarded = burn
      attempts = 100
      if (present(proposal_attempts)) attempts = proposal_attempts
      if (discarded < 0 .or. discarded >= iterations .or. attempts < 1) then
         out%info = 1
         return
      end if
      response_scale = max(finite_sample_variance(response), 1.0e-6_dp)
      mean_square = predictor_mean_squares(predictors)
      observation0 = max(0.1_dp*response_scale, 1.0e-8_dp)
      observation_shape = 2.0_dp
      observation_rate = observation0
      allocate(mean0(dimension), covariance0(dimension, dimension), &
         ar0(lags, variables), scaled0(variables), shape0(variables), &
         rate0(variables), state_normals(dimension, size(response), iterations), &
         ar_normals(lags, variables, attempts, iterations), &
         gammas(variables + 1, iterations))
      mean0 = 0.0_dp
      covariance0 = 0.0_dp
      do component = 1, variables
         covariance0((component - 1)*lags + 1:component*lags, &
            (component - 1)*lags + 1:component*lags) = &
            identity_block(lags, 10.0_dp*response_scale/mean_square(component))
      end do
      ar0 = 0.0_dp
      if (lags == 1) ar0 = 0.5_dp
      scaled0 = max(1.0e-4_dp*response_scale, 1.0e-10_dp)
      shape0 = 2.0_dp
      rate0 = scaled0
      observed = count(ieee_is_finite(response))
      do iteration = 1, iterations
         call fill_standard_normals(state_normals(:, :, iteration))
         do proposal = 1, attempts
            do component = 1, variables
               call fill_vector_normals(ar_normals(:, component, proposal, iteration))
            end do
         end do
         gammas(1, iteration) = random_gamma(observation_shape + &
            0.5_dp*real(observed, dp))
         do component = 1, variables
            gammas(component + 1, iteration) = random_gamma(shape0(component) + &
               0.5_dp*real(size(response) - 1, dp))
         end do
      end do
      if (present(offset_draws)) then
         out = bsts_dynamic_regression_ar_draws(response, predictors, lags, &
            mean0, covariance0, ar0, observation0, scaled0, mean_square, &
            observation_shape, observation_rate, shape0, rate0, discarded, &
            state_normals, ar_normals, gammas, offset_draws)
      else
         out = bsts_dynamic_regression_ar_draws(response, predictors, lags, &
            mean0, covariance0, ar0, observation0, scaled0, mean_square, &
            observation_shape, observation_rate, shape0, rate0, discarded, &
            state_normals, ar_normals, gammas)
      end if
   end function bsts_dynamic_regression_ar

   pure function bsts_dynamic_regression_ar_predict_draws(fit, &
      future_predictors, state_normal_draws, observation_normal_draws) result(out)
      !! Forecast AR dynamic coefficients from supplied independent normals.
      type(bsts_dynamic_regression_t), intent(in) :: fit !! Previously fitted model.
      real(dp), intent(in) :: future_predictors(:, :) !! Future predictors.
      real(dp), intent(in) :: state_normal_draws(:, :, :) !! Independent standard-normal state draws.
      real(dp), intent(in) :: observation_normal_draws(:, :) !! Observation normal draws.
      type(bsts_prediction_t) :: out
      real(dp), allocatable :: history(:, :), sorted(:), coefficient(:)
      integer :: variables, lags, horizon, retained, draw, source, step, lag
      integer :: component

      if (fit%info /= 0 .or. fit%ar_order < 1 .or. &
         .not. allocated(fit%coefficients) .or. &
         .not. allocated(fit%ar_coefficients) .or. &
         .not. allocated(fit%innovation_variance) .or. &
         .not. allocated(fit%residual_variance)) then
         out%info = 1
         return
      end if
      variables = size(fit%coefficients, 1)
      lags = fit%ar_order
      horizon = size(future_predictors, 1)
      retained = size(fit%coefficients, 3) - fit%burn
      if (size(fit%coefficients, 2) < lags .or. horizon < 1 .or. retained < 1 .or. &
         size(future_predictors, 2) /= variables .or. &
         any(shape(state_normal_draws) /= [variables, horizon, retained]) .or. &
         any(shape(observation_normal_draws) /= [horizon, retained]) .or. &
         .not. all(ieee_is_finite(future_predictors)) .or. &
         .not. all(ieee_is_finite(state_normal_draws)) .or. &
         .not. all(ieee_is_finite(observation_normal_draws))) then
         out%info = 1
         return
      end if
      allocate(out%draws(horizon, retained), out%mean(horizon), &
         out%standard_deviation(horizon), out%lower(horizon), &
         out%upper(horizon), history(lags, variables), coefficient(variables))
      do draw = 1, retained
         source = fit%burn + draw
         do lag = 1, lags
            history(lag, :) = fit%coefficients(:, &
               size(fit%coefficients, 2) - lag + 1, source)
         end do
         do step = 1, horizon
            do component = 1, variables
               coefficient(component) = dot_product( &
                  fit%ar_coefficients(:, component, source), &
                  history(:, component)) + &
                  sqrt(fit%innovation_variance(component, source))* &
                  state_normal_draws(component, step, draw)
            end do
            if (lags > 1) history(2:, :) = history(:lags - 1, :)
            history(1, :) = coefficient
            out%draws(step, draw) = dot_product(future_predictors(step, :), &
               coefficient) + sqrt(fit%residual_variance(source))* &
               observation_normal_draws(step, draw)
         end do
      end do
      do step = 1, horizon
         out%mean(step) = sum(out%draws(step, :))/real(retained, dp)
         if (retained > 1) then
            out%standard_deviation(step) = sqrt(sum((out%draws(step, :) - &
               out%mean(step))**2)/real(retained - 1, dp))
         else
            out%standard_deviation(step) = 0.0_dp
         end if
         sorted = out%draws(step, :)
         call insertion_sort(sorted)
         out%lower(step) = quantile(sorted, 0.025_dp)
         out%upper(step) = quantile(sorted, 0.975_dp)
      end do
   end function bsts_dynamic_regression_ar_predict_draws

   function bsts_dynamic_regression_ar_predict(fit, future_predictors) result(out)
      !! Forecast AR dynamic regression using the shared random stream.
      type(bsts_dynamic_regression_t), intent(in) :: fit !! Previously fitted model.
      real(dp), intent(in) :: future_predictors(:, :) !! Future predictors.
      type(bsts_prediction_t) :: out
      real(dp), allocatable :: state_normals(:, :, :), observation_normals(:, :)
      integer :: variables, horizon, retained, draw, step

      if (.not. allocated(fit%coefficients)) then
         out%info = 1
         return
      end if
      variables = size(fit%coefficients, 1)
      horizon = size(future_predictors, 1)
      retained = size(fit%coefficients, 3) - fit%burn
      if (horizon < 1 .or. retained < 1) then
         out%info = 1
         return
      end if
      allocate(state_normals(variables, horizon, retained), &
         observation_normals(horizon, retained))
      do draw = 1, retained
         do step = 1, horizon
            call fill_vector_normals(state_normals(:, step, draw))
            observation_normals(step, draw) = random_standard_normal()
         end do
      end do
      out = bsts_dynamic_regression_ar_predict_draws(fit, future_predictors, &
         state_normals, observation_normals)
   end function bsts_dynamic_regression_ar_predict

   pure subroutine set_dynamic_ar_transition(transition, coefficients, lags)
      !! Fill independent companion blocks for AR dynamic coefficients.
      real(dp), intent(out) :: transition(:, :) !! State transition matrix.
      real(dp), intent(in) :: coefficients(:, :) !! Model coefficients.
      integer, intent(in) :: lags !! Lags.
      integer :: component, lag, start

      transition = 0.0_dp
      do component = 1, size(coefficients, 2)
         start = (component - 1)*lags + 1
         transition(start, start:start + lags - 1) = coefficients(:, component)
         do lag = 2, lags
            transition(start + lag - 1, start + lag - 2) = 1.0_dp
         end do
      end do
   end subroutine set_dynamic_ar_transition

   pure subroutine set_ar_transition(transition, coefficients)
      !! Fill a companion transition matrix for one AR process.
      real(dp), intent(out) :: transition(:, :) !! State transition matrix.
      real(dp), intent(in) :: coefficients(:) !! Model coefficients.
      integer :: lag

      transition = 0.0_dp
      transition(1, :) = coefficients
      do lag = 2, size(coefficients)
         transition(lag, lag - 1) = 1.0_dp
      end do
   end subroutine set_ar_transition

   pure function ar_coefficients_stationary(coefficients) result(stationary)
      !! Test AR stationarity through the inverse Levinson recursion.
      real(dp), intent(in) :: coefficients(:) !! Model coefficients.
      logical :: stationary
      real(dp) :: work(size(coefficients)), previous(size(coefficients))
      real(dp) :: reflection, denominator
      integer :: order, lag

      work = coefficients
      stationary = all(ieee_is_finite(coefficients))
      do order = size(coefficients), 1, -1
         reflection = work(order)
         if (abs(reflection) >= 1.0_dp) then
            stationary = .false.
            return
         end if
         if (order == 1) exit
         denominator = 1.0_dp - reflection**2
         previous = 0.0_dp
         do lag = 1, order - 1
            previous(lag) = (work(lag) + &
               reflection*work(order - lag))/denominator
         end do
         work(:order - 1) = previous(:order - 1)
      end do
   end function ar_coefficients_stationary

   pure function predictor_mean_squares(predictors) result(mean_square)
      !! Compute positive predictor second moments used by AR state scaling.
      real(dp), intent(in) :: predictors(:, :) !! Predictor matrix.
      real(dp) :: mean_square(size(predictors, 2))
      integer :: component

      do component = 1, size(predictors, 2)
         mean_square(component) = sum(predictors(:, component)**2)/ &
            real(size(predictors, 1), dp)
         mean_square(component) = max(mean_square(component), 1.0e-12_dp)
      end do
   end function predictor_mean_squares

   pure function identity_block(dimension, scale) result(matrix)
      !! Construct a scaled identity matrix.
      integer, intent(in) :: dimension !! Dimension.
      real(dp), intent(in) :: scale !! Scale.
      real(dp) :: matrix(dimension, dimension)
      integer :: component

      matrix = 0.0_dp
      do component = 1, dimension
         matrix(component, component) = scale
      end do
   end function identity_block

   subroutine fill_vector_normals(draws)
      !! Fill a vector with independent standard normal draws.
      real(dp), intent(out) :: draws(:) !! Draws.
      integer :: component

      do component = 1, size(draws)
         draws(component) = random_standard_normal()
      end do
   end subroutine fill_vector_normals

   pure subroutine gamma_from_proposals(shape_parameter, normal_draws, &
      uniform_draws, value, info)
      !! Generate a unit-scale gamma draw from supplied Marsaglia proposals.
      real(dp), intent(in) :: shape_parameter !! Shape parameter.
      real(dp), intent(in) :: normal_draws(:) !! Independent standard-normal draws.
      real(dp), intent(in) :: uniform_draws(:) !! Uniform simulation draws.
      real(dp), intent(out) :: value !! Input value.
      integer, intent(out) :: info !! Status code; zero indicates success.
      real(dp) :: d, c, candidate, normal, uniform
      integer :: proposal

      if (shape_parameter < 1.0_dp .or. &
         size(normal_draws) /= size(uniform_draws) .or. &
         size(normal_draws) < 1) then
         value = 0.0_dp
         info = 1
         return
      end if
      d = shape_parameter - 1.0_dp/3.0_dp
      c = 1.0_dp/sqrt(9.0_dp*d)
      do proposal = 1, size(normal_draws)
         normal = normal_draws(proposal)
         candidate = (1.0_dp + c*normal)**3
         if (candidate <= 0.0_dp) cycle
         uniform = uniform_draws(proposal)
         if (uniform < 1.0_dp - 0.0331_dp*normal**4 .or. &
            log(uniform) < 0.5_dp*normal**2 + &
            d*(1.0_dp - candidate + log(candidate))) then
            value = d*candidate
            info = 0
            return
         end if
      end do
      value = 0.0_dp
      info = 2
   end subroutine gamma_from_proposals

   pure elemental real(dp) function reflected_value(value, lower, upper) &
      result(reflected)
      !! Reflect a proposal into a finite interval without point masses.
      real(dp), intent(in) :: value !! Input value.
      real(dp), intent(in) :: lower !! Lower.
      real(dp), intent(in) :: upper !! Upper.
      real(dp) :: position, width

      width = upper - lower
      position = modulo(value - lower, 2.0_dp*width)
      if (position <= width) then
         reflected = lower + position
      else
         reflected = upper - (position - width)
      end if
   end function reflected_value

   pure function student_weight_log_density(weights, degrees) result(value)
      !! Evaluate latent Student precision weights up to a common constant.
      real(dp), intent(in) :: weights(:) !! Observation or objective weights.
      real(dp), intent(in) :: degrees !! Degrees.
      real(dp) :: value, shape_parameter, rate

      shape_parameter = 0.5_dp*degrees
      rate = 0.5_dp*degrees
      value = real(size(weights), dp)*(shape_parameter*log(rate) - &
         log_gamma(shape_parameter)) + (shape_parameter - 1.0_dp)* &
         sum(log(weights)) - rate*sum(weights)
   end function student_weight_log_density

   pure function dynamic_regression_draws_core(response, predictors, &
      initial_mean, initial_covariance, residual_variance, innovation_variance, &
      residual_prior_shape, residual_prior_rate, innovation_prior_shape, &
      innovation_prior_rate, burn, state_normal_draws, gamma_draws, &
      hierarchical, predictor_variance, shrinkage_shape, initial_hierarchy_rate, &
      hierarchy_rate_prior_shape, hierarchy_rate_prior_rate, &
      offset_draws) result(out)
      !! Run conjugate Gibbs updates for a random-walk dynamic regression.
      real(dp), intent(in) :: response(:) !! Response observations.
      real(dp), intent(in) :: predictors(:, :) !! Predictor matrix.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: residual_variance !! Residual variance.
      real(dp), intent(in) :: innovation_variance(:) !! Innovation variance.
      real(dp), intent(in) :: residual_prior_shape !! Residual prior shape.
      real(dp), intent(in) :: residual_prior_rate !! Residual prior rate.
      real(dp), intent(in) :: innovation_prior_shape(:) !! Innovation prior shape.
      real(dp), intent(in) :: innovation_prior_rate(:) !! Innovation prior rate.
      integer, intent(in) :: burn !! Number of initial simulation draws to discard.
      real(dp), intent(in) :: state_normal_draws(:, :, :) !! Independent standard-normal state draws.
      real(dp), intent(in) :: gamma_draws(:, :) !! Gamma simulation draws.
      logical, intent(in) :: hierarchical !! Flag controlling hierarchical.
      real(dp), intent(in) :: predictor_variance(:) !! Predictor variance.
      real(dp), intent(in) :: shrinkage_shape !! Shrinkage shape.
      real(dp), intent(in) :: initial_hierarchy_rate !! Initial hierarchy rate.
      real(dp), intent(in) :: hierarchy_rate_prior_shape !! Hierarchy rate prior shape.
      real(dp), intent(in) :: hierarchy_rate_prior_rate !! Hierarchy rate prior rate.
      real(dp), intent(in), optional :: offset_draws(:, :) !! Offset simulation draws.
      type(bsts_dynamic_regression_t) :: out
      type(ssm_model_t) :: model
      real(dp), allocatable :: path(:, :), working_response(:)
      real(dp) :: observation, state(size(innovation_variance))
      real(dp) :: scaled(size(innovation_variance)), hierarchy_rate
      real(dp) :: residual_sum, difference_sum, rate
      integer :: observations, variables, iterations, iteration, time
      integer :: component, status, observed, required_gamma

      observations = size(response)
      variables = size(predictors, 2)
      iterations = size(state_normal_draws, 3)
      required_gamma = variables + 1
      if (hierarchical) required_gamma = variables + 2
      if (.not. dynamic_regression_data_valid(response, predictors) .or. &
         size(initial_mean) /= variables .or. &
         any(shape(initial_covariance) /= [variables, variables]) .or. &
         residual_variance <= 0.0_dp .or. &
         size(innovation_variance) /= variables .or. &
         any(innovation_variance <= 0.0_dp) .or. residual_prior_shape <= 0.0_dp .or. &
         residual_prior_rate <= 0.0_dp .or. &
         size(innovation_prior_shape) /= variables .or. &
         size(innovation_prior_rate) /= variables .or. &
         any(innovation_prior_shape <= 0.0_dp) .or. &
         any(innovation_prior_rate <= 0.0_dp) .or. iterations < 1 .or. &
         burn < 0 .or. burn >= iterations .or. &
         any(shape(state_normal_draws) /= [variables, observations, iterations]) .or. &
         any(shape(gamma_draws) /= [required_gamma, iterations]) .or. &
         any(gamma_draws <= 0.0_dp) .or. &
         .not. all(ieee_is_finite(state_normal_draws))) then
         out%info = 1
         return
      end if
      if (hierarchical) then
         if (size(predictor_variance) /= variables .or. &
            any(predictor_variance <= 0.0_dp) .or. shrinkage_shape <= 0.0_dp .or. &
            initial_hierarchy_rate <= 0.0_dp .or. &
            hierarchy_rate_prior_shape <= 0.0_dp .or. &
            hierarchy_rate_prior_rate <= 0.0_dp) then
            out%info = 1
            return
         end if
      end if
      if (present(offset_draws)) then
         if (any(shape(offset_draws) /= [observations, iterations]) .or. &
            .not. all(ieee_is_finite(offset_draws))) then
            out%info = 1
            return
         end if
      end if
      allocate(model%y(observations, 1), model%z(1, variables, observations), &
         model%h(1, 1, 1), model%transition(variables, variables, 1), &
         model%r(variables, variables, 1), model%q(variables, variables, 1), &
         model%a1(variables), model%p1(variables, variables), &
         model%p1inf(variables, variables), model%missing(observations, 1), &
         working_response(observations), &
         out%coefficients(variables, observations, iterations), &
         out%innovation_variance(variables, iterations), &
         out%scaled_innovation_variance(variables, iterations), &
         out%residual_variance(iterations), out%predictor_variance(variables))
      if (hierarchical) allocate(out%hierarchy_rate(iterations))
      model%transition = 0.0_dp
      model%r = 0.0_dp
      do component = 1, variables
         model%transition(component, component, 1) = 1.0_dp
         model%r(component, component, 1) = 1.0_dp
      end do
      do time = 1, observations
         model%z(1, :, time) = predictors(time, :)
      end do
      model%a1 = initial_mean
      model%p1 = initial_covariance
      model%p1inf = 0.0_dp
      model%missing(:, 1) = .not. ieee_is_finite(response)
      observation = residual_variance
      state = innovation_variance
      scaled = state
      if (hierarchical) scaled = state*predictor_variance
      hierarchy_rate = initial_hierarchy_rate
      out%burn = burn
      out%hierarchical = hierarchical
      out%predictor_variance = predictor_variance
      observed = count(ieee_is_finite(response))
      do iteration = 1, iterations
         working_response = response
         if (present(offset_draws)) working_response = working_response - &
            offset_draws(:, iteration)
         model%y(:, 1) = working_response
         where (.not. ieee_is_finite(model%y)) model%y = 0.0_dp
         model%h(1, 1, 1) = observation
         model%q = 0.0_dp
         do component = 1, variables
            model%q(component, component, 1) = state(component)
         end do
         call state_path_from_draws(model, state_normal_draws(:, :, iteration), &
            path, status)
         if (status /= 0) then
            out%info = 10 + status
            return
         end if
         out%coefficients(:, :, iteration) = path
         residual_sum = 0.0_dp
         do time = 1, observations
            if (ieee_is_finite(response(time))) residual_sum = residual_sum + &
               (working_response(time) - dot_product(predictors(time, :), &
               path(:, time)))**2
         end do
         observation = (residual_prior_rate + 0.5_dp*residual_sum)/ &
            gamma_draws(1, iteration)
         do component = 1, variables
            difference_sum = sum((path(component, 2:) - &
               path(component, :observations - 1))**2)
            if (hierarchical) then
               rate = hierarchy_rate + &
                  0.5_dp*predictor_variance(component)*difference_sum
               scaled(component) = rate/gamma_draws(component + 1, iteration)
               state(component) = scaled(component)/predictor_variance(component)
            else
               rate = innovation_prior_rate(component) + 0.5_dp*difference_sum
               state(component) = rate/gamma_draws(component + 1, iteration)
               scaled(component) = state(component)
            end if
         end do
         if (hierarchical) then
            hierarchy_rate = gamma_draws(variables + 2, iteration)/ &
               (hierarchy_rate_prior_rate + sum(1.0_dp/scaled))
            out%hierarchy_rate(iteration) = hierarchy_rate
         end if
         out%residual_variance(iteration) = observation
         out%innovation_variance(:, iteration) = state
         out%scaled_innovation_variance(:, iteration) = scaled
      end do
   end function dynamic_regression_draws_core

   pure function dynamic_regression_data_valid(response, predictors) result(valid)
      !! Report whether dynamic-regression observations and predictors are usable.
      real(dp), intent(in) :: response(:) !! Response observations.
      real(dp), intent(in) :: predictors(:, :) !! Predictor matrix.
      logical :: valid

      valid = size(response) >= 2 .and. size(predictors, 1) == size(response) .and. &
         size(predictors, 2) >= 1 .and. count(ieee_is_finite(response)) >= 1 .and. &
         all(ieee_is_finite(predictors))
   end function dynamic_regression_data_valid

   pure function predictor_variances(predictors) result(variances)
      !! Compute positive sample variances used to scale dynamic coefficients.
      real(dp), intent(in) :: predictors(:, :) !! Predictor matrix.
      real(dp) :: variances(size(predictors, 2)), center
      integer :: component

      do component = 1, size(predictors, 2)
         center = sum(predictors(:, component))/real(size(predictors, 1), dp)
         variances(component) = sum((predictors(:, component) - center)**2)/ &
            real(max(1, size(predictors, 1) - 1), dp)
         variances(component) = max(variances(component), 1.0e-12_dp)
      end do
   end function predictor_variances

   pure function structural_gibbs_draws(y, initial_mean, initial_covariance, &
      transition, observation, observation_variance, state_variance, &
      observation_prior_shape, observation_prior_rate, state_prior_shape, &
      state_prior_rate, burn, state_normal_draws, gamma_draws) result(out)
      !! Run conjugate structural Gibbs updates for diagonal state disturbances.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: transition(:, :) !! State transition matrix.
      real(dp), intent(in) :: observation(:) !! Observed value or vector.
      real(dp), intent(in) :: observation_variance !! Observation-error variance.
      real(dp), intent(in) :: state_variance(:) !! State variance.
      real(dp), intent(in) :: observation_prior_shape !! Observation prior shape.
      real(dp), intent(in) :: observation_prior_rate !! Observation prior rate.
      real(dp), intent(in) :: state_prior_shape(:) !! State prior shape.
      real(dp), intent(in) :: state_prior_rate(:) !! State prior rate.
      integer, intent(in) :: burn !! Number of initial simulation draws to discard.
      real(dp), intent(in) :: state_normal_draws(:, :, :) !! Independent standard-normal state draws.
      real(dp), intent(in) :: gamma_draws(:, :) !! Gamma simulation draws.
      type(bsts_mcmc_t) :: out
      type(ssm_model_t) :: model
      real(dp), allocatable :: path(:, :), previous(:), disturbance(:)
      real(dp) :: observation_scale, state_scale(size(state_variance))
      real(dp) :: observation_sum, state_sum(size(state_variance))
      integer :: dimension, iterations, iteration, time, component, status

      dimension = size(initial_mean)
      iterations = size(state_normal_draws, 3)
      if (size(y) < 2 .or. dimension < 1 .or. iterations < 1 .or. &
         any(shape(initial_covariance) /= [dimension, dimension]) .or. &
         any(shape(transition) /= [dimension, dimension]) .or. &
         size(observation) /= dimension .or. size(state_variance) /= dimension .or. &
         size(state_prior_shape) /= dimension .or. &
         size(state_prior_rate) /= dimension .or. &
         size(state_normal_draws, 1) /= dimension .or. &
         size(state_normal_draws, 2) /= size(y) .or. &
         any(shape(gamma_draws) /= [dimension + 1, iterations]) .or. &
         burn < 0 .or. burn >= iterations .or. observation_variance <= 0.0_dp .or. &
         any(state_variance <= 0.0_dp) .or. observation_prior_shape <= 0.0_dp .or. &
         observation_prior_rate <= 0.0_dp .or. any(state_prior_shape <= 0.0_dp) .or. &
         any(state_prior_rate <= 0.0_dp) .or. any(gamma_draws <= 0.0_dp) .or. &
         .not. all(ieee_is_finite(state_normal_draws))) then
         out%info = 1
         return
      end if
      allocate(out%state(dimension, size(y), iterations), &
         out%observation_variance(iterations), &
         out%state_variance(dimension, iterations), &
         out%component_variance(dimension, iterations), &
         out%transition(dimension, dimension), out%observation(dimension), &
         out%transition_schedule(dimension, dimension, 1), &
         out%state_loading_schedule(dimension, dimension, 1), &
         previous(dimension), disturbance(dimension))
      observation_scale = observation_variance
      state_scale = state_variance
      out%transition = transition
      out%observation = observation
      out%transition_schedule(:, :, 1) = transition
      out%state_loading_schedule = 0.0_dp
      do component = 1, dimension
         out%state_loading_schedule(component, component, 1) = 1.0_dp
      end do
      out%burn = burn
      do iteration = 1, iterations
         model = structural_model(y, initial_mean, initial_covariance, &
            transition, observation, observation_scale, state_scale)
         call state_path_from_draws(model, state_normal_draws(:, :, iteration), &
            path, status)
         if (status /= 0) then
            out%info = 10 + status
            return
         end if
         out%state(:, :, iteration) = path
         observation_sum = 0.0_dp
         do time = 1, size(y)
            if (ieee_is_finite(y(time))) observation_sum = observation_sum + &
               (y(time) - dot_product(observation, path(:, time)))**2
         end do
         observation_scale = (observation_prior_rate + &
            0.5_dp*observation_sum)/gamma_draws(1, iteration)
         state_sum = 0.0_dp
         do time = 2, size(y)
            previous = matmul(transition, path(:, time - 1))
            disturbance = path(:, time) - previous
            state_sum = state_sum + disturbance**2
         end do
         do component = 1, dimension
            state_scale(component) = (state_prior_rate(component) + &
               0.5_dp*state_sum(component))/gamma_draws(component + 1, iteration)
         end do
         out%observation_variance(iteration) = observation_scale
         out%state_variance(:, iteration) = state_scale
         out%component_variance(:, iteration) = state_scale
      end do
   end function structural_gibbs_draws

   pure function seasonal_model(y, nseasons, season_duration, initial_mean, &
      initial_covariance, observation_variance, seasonal_variance) result(model)
      !! Construct bsts's duration-aware sum-to-zero seasonal state model.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      integer, intent(in) :: nseasons !! Nseasons.
      integer, intent(in) :: season_duration !! Season duration.
      real(dp), intent(in) :: observation_variance !! Observation-error variance.
      real(dp), intent(in) :: seasonal_variance !! Seasonal variance.
      type(ssm_model_t) :: model
      real(dp), allocatable :: seasonal_transition(:, :)
      integer :: dimension, time, component

      dimension = nseasons - 1
      allocate(seasonal_transition(dimension, dimension), &
         model%y(size(y), 1), model%z(1, dimension, 1), model%h(1, 1, 1), &
         model%transition(dimension, dimension, size(y)), &
         model%r(dimension, 1, size(y)), model%q(1, 1, size(y)), &
         model%a1(dimension), model%p1(dimension, dimension), &
         model%p1inf(dimension, dimension), model%missing(size(y), 1))
      seasonal_transition = 0.0_dp
      seasonal_transition(1, :) = -1.0_dp
      do component = 2, dimension
         seasonal_transition(component, component - 1) = 1.0_dp
      end do
      model%y(:, 1) = y
      where (.not. ieee_is_finite(model%y)) model%y = 0.0_dp
      model%z = 0.0_dp
      model%z(1, 1, 1) = 1.0_dp
      model%h(1, 1, 1) = observation_variance
      model%transition = 0.0_dp
      model%r = 0.0_dp
      model%q(1, 1, :) = seasonal_variance
      do time = 1, size(y)
         if (modulo(time, season_duration) == 0) then
            model%transition(:, :, time) = seasonal_transition
            model%r(1, 1, time) = 1.0_dp
         else
            do component = 1, dimension
               model%transition(component, component, time) = 1.0_dp
            end do
         end if
      end do
      model%a1 = initial_mean
      model%p1 = initial_covariance
      model%p1inf = 0.0_dp
      model%missing(:, 1) = .not. ieee_is_finite(y)
   end function seasonal_model

   pure function trig_transition_matrix(period, frequencies) result(transition)
      !! Build harmonic rotation blocks for bsts's trigonometric state model.
      real(dp), intent(in) :: period !! Seasonal period.
      real(dp), intent(in) :: frequencies(:) !! Frequencies.
      real(dp), allocatable :: transition(:, :)
      real(dp) :: angle
      integer :: component, first

      allocate(transition(2*size(frequencies), 2*size(frequencies)))
      transition = 0.0_dp
      do component = 1, size(frequencies)
         first = 2*component - 1
         angle = 2.0_dp*acos(-1.0_dp)*frequencies(component)/period
         transition(first, first) = cos(angle)
         transition(first, first + 1) = sin(angle)
         transition(first + 1, first) = -sin(angle)
         transition(first + 1, first + 1) = cos(angle)
      end do
   end function trig_transition_matrix

   pure function structural_model(y, initial_mean, initial_covariance, &
      transition, observation, observation_variance, state_variance) result(model)
      !! Construct a time-invariant Gaussian structural state-space model.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      real(dp), intent(in) :: transition(:, :) !! State transition matrix.
      real(dp), intent(in) :: observation(:) !! Observed value or vector.
      real(dp), intent(in) :: observation_variance !! Observation-error variance.
      real(dp), intent(in) :: state_variance(:) !! State variance.
      type(ssm_model_t) :: model
      integer :: dimension, component

      dimension = size(initial_mean)
      allocate(model%y(size(y), 1), model%z(1, dimension, 1), &
         model%h(1, 1, 1), model%transition(dimension, dimension, 1), &
         model%r(dimension, dimension, 1), model%q(dimension, dimension, 1), &
         model%a1(dimension), model%p1(dimension, dimension), &
         model%p1inf(dimension, dimension), model%missing(size(y), 1))
      model%y(:, 1) = y
      where (.not. ieee_is_finite(model%y)) model%y = 0.0_dp
      model%z(1, :, 1) = observation
      model%h(1, 1, 1) = observation_variance
      model%transition(:, :, 1) = transition
      model%r = 0.0_dp
      model%q = 0.0_dp
      do component = 1, dimension
         model%r(component, component, 1) = 1.0_dp
         model%q(component, component, 1) = state_variance(component)
      end do
      model%a1 = initial_mean
      model%p1 = initial_covariance
      model%p1inf = 0.0_dp
      model%missing(:, 1) = .not. ieee_is_finite(y)
   end function structural_model

   pure function structural_model_from_fit(y, fit, source, initial_mean, &
      initial_covariance) result(model)
      !! Reconstruct a Gaussian state-space model from one posterior draw.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      type(bsts_mcmc_t), intent(in) :: fit !! Previously fitted model.
      integer, intent(in) :: source !! Source.
      type(ssm_model_t) :: model
      real(dp), allocatable :: variance(:)
      integer :: dimension, disturbances, slices, z_slices, r_slices
      integer :: time, component

      dimension = size(initial_mean)
      disturbances = dimension
      if (allocated(fit%state_loading_schedule)) disturbances = &
         size(fit%state_loading_schedule, 2)
      slices = 1
      if (allocated(fit%transition_schedule)) slices = &
         size(fit%transition_schedule, 3)
      z_slices = 1
      if (allocated(fit%observation_schedule)) z_slices = &
         size(fit%observation_schedule, 2)
      r_slices = 1
      if (allocated(fit%state_loading_schedule)) r_slices = &
         size(fit%state_loading_schedule, 3)
      allocate(model%y(size(y), 1), model%z(1, dimension, z_slices), &
         model%h(1, 1, 1), &
         model%transition(dimension, dimension, slices), &
         model%r(dimension, disturbances, r_slices), &
         model%q(disturbances, disturbances, 1), model%a1(dimension), &
         model%p1(dimension, dimension), model%p1inf(dimension, dimension), &
         model%missing(size(y), 1))
      model%y(:, 1) = y
      where (.not. ieee_is_finite(model%y)) model%y = 0.0_dp
      if (allocated(fit%observation_schedule)) then
         do time = 1, size(fit%observation_schedule, 2)
            model%z(1, :, time) = fit%observation_schedule(:, time)
         end do
      else
         model%z(1, :, 1) = fit%observation
      end if
      model%h = fit%observation_variance(source)
      if (allocated(fit%transition_schedule)) then
         model%transition = fit%transition_schedule
      else
         model%transition(:, :, 1) = fit%transition
      end if
      model%r = 0.0_dp
      if (allocated(fit%state_loading_schedule)) then
         model%r = fit%state_loading_schedule
      else
         do component = 1, min(dimension, disturbances)
            model%r(component, component, 1) = 1.0_dp
         end do
      end if
      if (allocated(fit%component_variance)) then
         variance = fit%component_variance(:, source)
      else
         variance = fit%state_variance(:, source)
      end if
      model%q = 0.0_dp
      do component = 1, disturbances
         model%q(component, component, 1) = &
            variance(min(component, size(variance)))
      end do
      model%a1 = initial_mean
      model%p1 = initial_covariance
      model%p1inf = 0.0_dp
      model%missing(:, 1) = .not. ieee_is_finite(y)
   end function structural_model_from_fit

   pure function dirm_local_level_model(y, counts, initial_mean, &
      initial_variance, observation_variance, level_variance) result(model)
      !! Construct a local-level model for grouped Gaussian observations.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean !! Initial state mean.
      real(dp), intent(in) :: initial_variance !! Initial variance.
      integer, intent(in) :: counts(:) !! Counts.
      real(dp), intent(in) :: observation_variance !! Observation-error variance.
      real(dp), intent(in) :: level_variance !! Level variance.
      type(ssm_model_t) :: model
      integer :: time

      allocate(model%y(size(y), 1), model%z(1, 1, 1), &
         model%h(1, 1, size(y)), model%transition(1, 1, 1), &
         model%r(1, 1, 1), model%q(1, 1, 1), model%a1(1), &
         model%p1(1, 1), model%p1inf(1, 1), model%missing(size(y), 1))
      model%y(:, 1) = y
      model%z = 1.0_dp
      do time = 1, size(y)
         model%h(1, 1, time) = observation_variance/ &
            real(max(1, counts(time)), dp)
      end do
      model%transition = 1.0_dp
      model%r = 1.0_dp
      model%q = level_variance
      model%a1 = initial_mean
      model%p1 = initial_variance
      model%p1inf = 0.0_dp
      model%missing(:, 1) = counts <= 0
   end function dirm_local_level_model

   pure function mixed_trend_seasonal_model(y, nseasons, duration, &
      initial_mean, initial_covariance, observation_variance, &
      component_variance) result(model)
      !! Construct a local-linear-trend and dummy-seasonal state model.
      real(dp), intent(in) :: y(:) !! Response or time-series observations.
      real(dp), intent(in) :: initial_mean(:) !! Initial state mean.
      real(dp), intent(in) :: initial_covariance(:, :) !! Initial state covariance matrix.
      integer, intent(in) :: nseasons !! Nseasons.
      integer, intent(in) :: duration !! Duration.
      real(dp), intent(in) :: observation_variance !! Observation-error variance.
      real(dp), intent(in) :: component_variance(3) !! Component variance.
      type(ssm_model_t) :: model
      real(dp), allocatable :: seasonal_transition(:, :)
      integer :: dimension, seasonal_dimension, time, component

      dimension = nseasons + 1
      seasonal_dimension = nseasons - 1
      allocate(seasonal_transition(seasonal_dimension, seasonal_dimension), &
         model%y(size(y), 1), model%z(1, dimension, 1), &
         model%h(1, 1, 1), model%transition(dimension, dimension, size(y)), &
         model%r(dimension, 3, size(y)), model%q(3, 3, 1), &
         model%a1(dimension), model%p1(dimension, dimension), &
         model%p1inf(dimension, dimension), model%missing(size(y), 1))
      seasonal_transition = 0.0_dp
      seasonal_transition(1, :) = -1.0_dp
      do component = 2, seasonal_dimension
         seasonal_transition(component, component - 1) = 1.0_dp
      end do
      model%y(:, 1) = y
      model%z = 0.0_dp
      model%z(1, 1, 1) = 1.0_dp
      model%z(1, 3, 1) = 1.0_dp
      model%h = observation_variance
      model%transition = 0.0_dp
      model%r = 0.0_dp
      model%q = 0.0_dp
      model%q(1, 1, 1) = component_variance(1)
      model%q(2, 2, 1) = component_variance(2)
      model%q(3, 3, 1) = component_variance(3)
      do time = 1, size(y)
         model%transition(1, 1, time) = 1.0_dp
         model%transition(1, 2, time) = 1.0_dp
         model%transition(2, 2, time) = 1.0_dp
         model%r(1, 1, time) = 1.0_dp
         model%r(2, 2, time) = 1.0_dp
         if (modulo(time, duration) == 0) then
            model%transition(3:, 3:, time) = seasonal_transition
            model%r(3, 3, time) = 1.0_dp
         else
            do component = 3, dimension
               model%transition(component, component, time) = 1.0_dp
            end do
         end if
      end do
      model%a1 = initial_mean
      model%p1 = initial_covariance
      model%p1inf = 0.0_dp
      model%missing = .false.
   end function mixed_trend_seasonal_model

   pure logical function dirm_time_index_valid(time_index) result(valid)
      !! Check for ordered contiguous time groups beginning at one.
      integer, intent(in) :: time_index(:) !! Index of time.
      integer :: row

      valid = size(time_index) > 0
      if (.not. valid) return
      valid = time_index(1) == 1
      do row = 2, size(time_index)
         if (time_index(row) < time_index(row - 1) .or. &
            time_index(row) > time_index(row - 1) + 1) then
            valid = .false.
            return
         end if
      end do
   end function dirm_time_index_valid

   pure subroutine state_path_from_draws(model, normal_draws, path, info)
      !! Draw one state trajectory by backward sampling filtered KFAS moments.
      type(ssm_model_t), intent(in) :: model !! Model specification.
      real(dp), intent(in) :: normal_draws(:, :) !! Independent standard-normal draws.
      real(dp), allocatable, intent(out) :: path(:, :) !! Path.
      integer, intent(out) :: info !! Status code; zero indicates success.
      type(kfs_filter_t) :: filtered
      real(dp), allocatable :: inverse(:, :), gain(:, :), mean(:)
      real(dp), allocatable :: covariance(:, :), transition(:, :)
      integer :: dimension, observations, time, status, slice

      dimension = size(model%a1)
      observations = size(model%y, 1)
      if (any(shape(normal_draws) /= [dimension, observations])) then
         info = 1
         return
      end if
      filtered = kfs_filter(model)
      if (filtered%info /= 0) then
         info = 10 + filtered%info
         return
      end if
      allocate(path(dimension, observations), inverse(dimension, dimension), &
         gain(dimension, dimension), mean(dimension), &
         covariance(dimension, dimension), transition(dimension, dimension))
      call multivariate_normal_from_psd(filtered%a_filt(:, observations), &
         filtered%p_filt(:, :, observations), normal_draws(:, observations), &
         path(:, observations), status)
      if (status /= 0) then
         info = 20 + status
         return
      end if
      do time = observations - 1, 1, -1
         slice = min(time, size(model%transition, 3))
         transition = model%transition(:, :, slice)
         call invert_matrix(filtered%p_pred(:, :, time + 1), inverse, status)
         if (status /= 0) then
            info = 30 + status
            return
         end if
         gain = matmul(matmul(filtered%p_filt(:, :, time), &
            transpose(transition)), inverse)
         mean = filtered%a_filt(:, time) + matmul(gain, &
            path(:, time + 1) - filtered%a_pred(:, time + 1))
         covariance = symmetrize(filtered%p_filt(:, :, time) - &
            matmul(gain, matmul(filtered%p_pred(:, :, time + 1), &
            transpose(gain))))
         call multivariate_normal_from_psd(mean, covariance, &
            normal_draws(:, time), path(:, time), status)
         if (status /= 0) then
            info = 40 + status
            return
         end if
      end do
      info = 0
   end subroutine state_path_from_draws

   pure subroutine multivariate_normal_from_psd(mean, covariance, standard, &
      draw, info)
      !! Draw from a covariance after removing numerical negative eigenvalues.
      real(dp), intent(in) :: mean(:) !! Mean value or vector.
      real(dp), intent(in) :: covariance(:, :) !! Covariance matrix.
      real(dp), intent(in) :: standard(:) !! Standard.
      real(dp), intent(out) :: draw(:) !! Draw.
      integer, intent(out) :: info !! Status code; zero indicates success.
      real(dp), allocatable :: eigenvalues(:), eigenvectors(:, :)
      real(dp) :: a, b, c, root, first_value, second_value
      real(dp) :: first_vector(2), vector_norm

      if (size(covariance, 1) /= size(mean) .or. &
         size(covariance, 2) /= size(mean) .or. &
         size(standard) /= size(mean) .or. size(draw) /= size(mean)) then
         info = 1
         return
      end if
      if (.not. all(ieee_is_finite(mean)) .or. &
         .not. all(ieee_is_finite(covariance)) .or. &
         .not. all(ieee_is_finite(standard))) then
         info = 2
         return
      end if
      if (size(mean) == 1) then
         draw(1) = mean(1) + sqrt(max(covariance(1, 1), 0.0_dp))*standard(1)
         info = 0
         return
      end if
      if (size(mean) == 2) then
         a = covariance(1, 1)
         b = 0.5_dp*(covariance(1, 2) + covariance(2, 1))
         c = covariance(2, 2)
         root = hypot(a - c, 2.0_dp*b)
         first_value = 0.5_dp*(a + c + root)
         second_value = 0.5_dp*(a + c - root)
         if (abs(b) > epsilon(1.0_dp)*max(1.0_dp, abs(first_value))) then
            first_vector = [b, first_value - a]
            vector_norm = sqrt(sum(first_vector**2))
            first_vector = first_vector/vector_norm
         else if (a >= c) then
            first_vector = [1.0_dp, 0.0_dp]
         else
            first_vector = [0.0_dp, 1.0_dp]
         end if
         draw = mean + sqrt(max(first_value, 0.0_dp))*standard(1)* &
            first_vector + sqrt(max(second_value, 0.0_dp))*standard(2)* &
            [-first_vector(2), first_vector(1)]
         info = 0
         return
      end if
      call symmetric_eigen(symmetrize(covariance), eigenvalues, eigenvectors, info)
      if (info /= 0) return
      draw = mean + matmul(eigenvectors, &
         sqrt(max(eigenvalues, 0.0_dp))*standard)
   end subroutine multivariate_normal_from_psd

   pure function bsts_predict_draws(fit, horizon, state_normal_draws, &
      observation_normal_draws) result(out)
      !! Simulate posterior forecasts from supplied independent normal draws.
      type(bsts_mcmc_t), intent(in) :: fit !! Previously fitted model.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      real(dp), intent(in) :: state_normal_draws(:, :, :) !! Independent standard-normal state draws.
      real(dp), intent(in) :: observation_normal_draws(:, :) !! Observation normal draws.
      type(bsts_prediction_t) :: out
      real(dp), allocatable :: state(:), sorted(:)
      real(dp), allocatable :: loading(:, :), transition(:, :)
      integer :: retained, dimension, disturbances, draw, source, step, phase

      if (fit%ar_order > 0 .or. allocated(fit%degrees_of_freedom) .or. &
         fit%is_monthly .or. fit%is_holiday) then
         out%info = 2
         return
      end if

      if (.not. allocated(fit%state) .or. &
         .not. allocated(fit%state_variance) .or. &
         .not. allocated(fit%observation_variance) .or. &
         .not. allocated(fit%transition_schedule) .or. &
         .not. allocated(fit%state_loading_schedule) .or. &
         .not. allocated(fit%observation)) then
         out%info = 1
         return
      end if
      if (fit%is_semilocal) then
         if (.not. allocated(fit%slope_mean) .or. &
            .not. allocated(fit%slope_ar)) then
            out%info = 1
            return
         end if
      end if
      dimension = size(fit%state, 1)
      disturbances = size(fit%state_loading_schedule, 2)
      retained = size(fit%state, 3) - fit%burn
      if (fit%info /= 0 .or. horizon < 1 .or. retained < 1 .or. &
         size(fit%state_variance, 1) /= disturbances .or. &
         size(fit%state_variance, 2) /= size(fit%state, 3) .or. &
         size(fit%transition_schedule, 1) /= dimension .or. &
         size(fit%transition_schedule, 2) /= dimension .or. &
         size(fit%state_loading_schedule, 1) /= dimension .or. &
         size(fit%state_loading_schedule, 3) /= &
         size(fit%transition_schedule, 3) .or. fit%forecast_phase < 1 .or. &
         any(shape(state_normal_draws) /= [disturbances, horizon, retained]) .or. &
         any(shape(observation_normal_draws) /= [horizon, retained]) .or. &
         .not. all(ieee_is_finite(state_normal_draws)) .or. &
         .not. all(ieee_is_finite(observation_normal_draws))) then
         out%info = 1
         return
      end if
      allocate(out%draws(horizon, retained), out%mean(horizon), &
         out%standard_deviation(horizon), out%lower(horizon), &
         out%upper(horizon), state(dimension))
      allocate(loading(dimension, disturbances), &
         transition(dimension, dimension))
      do draw = 1, retained
         source = fit%burn + draw
         state = fit%state(:, size(fit%state, 2), source)
         do step = 1, horizon
            phase = 1 + modulo(fit%forecast_phase + step - 2, &
               size(fit%transition_schedule, 3))
            transition = fit%transition_schedule(:, :, phase)
            loading = fit%state_loading_schedule(:, :, phase)
            if (fit%is_semilocal) then
               state = [state(1) + state(2), fit%slope_mean(source) + &
                  fit%slope_ar(source)*(state(2) - fit%slope_mean(source))] + &
                  matmul(loading, sqrt(fit%state_variance(:, source))* &
                  state_normal_draws(:, step, draw))
            else
               state = matmul(transition, state) + matmul(loading, &
                  sqrt(fit%state_variance(:, source))* &
                  state_normal_draws(:, step, draw))
            end if
            out%draws(step, draw) = dot_product(fit%observation, state) + &
               sqrt(fit%observation_variance(source))* &
               observation_normal_draws(step, draw)
         end do
      end do
      do step = 1, horizon
         out%mean(step) = sum(out%draws(step, :))/real(retained, dp)
         if (retained > 1) then
            out%standard_deviation(step) = sqrt(sum((out%draws(step, :) - &
               out%mean(step))**2)/real(retained - 1, dp))
         else
            out%standard_deviation(step) = 0.0_dp
         end if
         sorted = out%draws(step, :)
         call insertion_sort(sorted)
         out%lower(step) = quantile(sorted, 0.025_dp)
         out%upper(step) = quantile(sorted, 0.975_dp)
      end do
   end function bsts_predict_draws

   function bsts_predict(fit, horizon) result(out)
      !! Simulate posterior forecasts using the shared random stream.
      type(bsts_mcmc_t), intent(in) :: fit !! Previously fitted model.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      type(bsts_prediction_t) :: out
      real(dp), allocatable :: state_draws(:, :, :), observation_draws(:, :)
      integer :: retained, draw, step

      if (.not. allocated(fit%state) .or. &
         .not. allocated(fit%state_loading_schedule)) then
         out%info = 1
         return
      end if
      retained = size(fit%state, 3) - fit%burn
      if (horizon < 1 .or. retained < 1) then
         out%info = 1
         return
      end if
      allocate(state_draws(size(fit%state_loading_schedule, 2), &
         horizon, retained), &
         observation_draws(horizon, retained))
      do draw = 1, retained
         do step = 1, horizon
            call fill_standard_normals(state_draws(:, step:step, draw))
            observation_draws(step, draw) = random_standard_normal()
         end do
      end do
      out = bsts_predict_draws(fit, horizon, state_draws, observation_draws)
   end function bsts_predict

   pure function bsts_logit_trend_seasonal_draws(successes, trials, nseasons, &
      season_duration, initial_trend_mean, initial_trend_covariance, &
      seasonal_initial_variance, component_variance, component_prior_shape, &
      component_prior_rate, burn, state_normal_draws, acceptance_uniform_draws, &
      gamma_draws) result(out)
      !! Sample binomial-logit local trend and seasonal structural states.
      real(dp), intent(in) :: successes(:) !! Successes.
      real(dp), intent(in) :: trials(:) !! Trials.
      integer, intent(in) :: nseasons !! Nseasons.
      integer, intent(in) :: season_duration !! Season duration.
      integer, intent(in) :: burn !! Number of initial simulation draws to discard.
      real(dp), intent(in) :: initial_trend_mean(2) !! Initial trend mean.
      real(dp), intent(in) :: initial_trend_covariance(2, 2) !! Initial trend covariance.
      real(dp), intent(in) :: seasonal_initial_variance !! Seasonal initial variance.
      real(dp), intent(in) :: component_variance(3) !! Component variance.
      real(dp), intent(in) :: component_prior_shape(3) !! Component prior shape.
      real(dp), intent(in) :: component_prior_rate(3) !! Component prior rate.
      real(dp), intent(in) :: state_normal_draws(:, :, :) !! Independent standard-normal state draws.
      real(dp), intent(in) :: acceptance_uniform_draws(:) !! Acceptance uniform simulation draws.
      real(dp), intent(in) :: gamma_draws(:, :) !! Gamma simulation draws.
      type(bsts_non_gaussian_t) :: out

      out = non_gaussian_trend_seasonal_draws(successes, trials, 1, nseasons, &
         season_duration, initial_trend_mean, initial_trend_covariance, &
         seasonal_initial_variance, component_variance, component_prior_shape, &
         component_prior_rate, burn, state_normal_draws, &
         acceptance_uniform_draws, gamma_draws)
   end function bsts_logit_trend_seasonal_draws

   pure function bsts_poisson_trend_seasonal_draws(counts, exposure, nseasons, &
      season_duration, initial_trend_mean, initial_trend_covariance, &
      seasonal_initial_variance, component_variance, component_prior_shape, &
      component_prior_rate, burn, state_normal_draws, acceptance_uniform_draws, &
      gamma_draws) result(out)
      !! Sample Poisson local trend and seasonal states with exposures.
      real(dp), intent(in) :: counts(:) !! Counts.
      real(dp), intent(in) :: exposure(:) !! Exposure.
      integer, intent(in) :: nseasons !! Nseasons.
      integer, intent(in) :: season_duration !! Season duration.
      integer, intent(in) :: burn !! Number of initial simulation draws to discard.
      real(dp), intent(in) :: initial_trend_mean(2) !! Initial trend mean.
      real(dp), intent(in) :: initial_trend_covariance(2, 2) !! Initial trend covariance.
      real(dp), intent(in) :: seasonal_initial_variance !! Seasonal initial variance.
      real(dp), intent(in) :: component_variance(3) !! Component variance.
      real(dp), intent(in) :: component_prior_shape(3) !! Component prior shape.
      real(dp), intent(in) :: component_prior_rate(3) !! Component prior rate.
      real(dp), intent(in) :: state_normal_draws(:, :, :) !! Independent standard-normal state draws.
      real(dp), intent(in) :: acceptance_uniform_draws(:) !! Acceptance uniform simulation draws.
      real(dp), intent(in) :: gamma_draws(:, :) !! Gamma simulation draws.
      type(bsts_non_gaussian_t) :: out

      out = non_gaussian_trend_seasonal_draws(counts, exposure, 2, nseasons, &
         season_duration, initial_trend_mean, initial_trend_covariance, &
         seasonal_initial_variance, component_variance, component_prior_shape, &
         component_prior_rate, burn, state_normal_draws, &
         acceptance_uniform_draws, gamma_draws)
   end function bsts_poisson_trend_seasonal_draws

   function bsts_logit_trend_seasonal(successes, trials, nseasons, iterations, &
      season_duration, burn) result(out)
      !! Sample binomial-logit trend and seasonal states using random draws.
      real(dp), intent(in) :: successes(:) !! Successes.
      real(dp), intent(in) :: trials(:) !! Trials.
      integer, intent(in) :: nseasons !! Nseasons.
      integer, intent(in) :: iterations !! Number of algorithm iterations.
      integer, intent(in), optional :: season_duration !! Season duration.
      integer, intent(in), optional :: burn !! Number of initial simulation draws to discard.
      type(bsts_non_gaussian_t) :: out

      out = random_non_gaussian_trend_seasonal(successes, trials, 1, &
         nseasons, iterations, season_duration, burn)
   end function bsts_logit_trend_seasonal

   function bsts_poisson_trend_seasonal(counts, exposure, nseasons, iterations, &
      season_duration, burn) result(out)
      !! Sample Poisson trend and seasonal states using random draws.
      real(dp), intent(in) :: counts(:) !! Counts.
      real(dp), intent(in) :: exposure(:) !! Exposure.
      integer, intent(in) :: nseasons !! Nseasons.
      integer, intent(in) :: iterations !! Number of algorithm iterations.
      integer, intent(in), optional :: season_duration !! Season duration.
      integer, intent(in), optional :: burn !! Number of initial simulation draws to discard.
      type(bsts_non_gaussian_t) :: out

      out = random_non_gaussian_trend_seasonal(counts, exposure, 2, &
         nseasons, iterations, season_duration, burn)
   end function bsts_poisson_trend_seasonal

   pure function bsts_logit_trend_seasonal_predict_draws(fit, future_trials, &
      state_normal_draws, observation_uniform_draws) result(out)
      !! Forecast binomial observations from trend and seasonal states.
      type(bsts_non_gaussian_t), intent(in) :: fit !! Previously fitted model.
      real(dp), intent(in) :: future_trials(:) !! Future trials.
      real(dp), intent(in) :: state_normal_draws(:, :, :) !! Independent standard-normal state draws.
      real(dp), intent(in) :: observation_uniform_draws(:, :) !! Observation uniform draws.
      type(bsts_prediction_t) :: out

      out = non_gaussian_trend_seasonal_predict_draws(fit, future_trials, &
         state_normal_draws, observation_uniform_draws, 1)
   end function bsts_logit_trend_seasonal_predict_draws

   pure function bsts_poisson_trend_seasonal_predict_draws(fit, future_exposure, &
      state_normal_draws, observation_uniform_draws) result(out)
      !! Forecast Poisson observations from trend and seasonal states.
      type(bsts_non_gaussian_t), intent(in) :: fit !! Previously fitted model.
      real(dp), intent(in) :: future_exposure(:) !! Future exposure.
      real(dp), intent(in) :: state_normal_draws(:, :, :) !! Independent standard-normal state draws.
      real(dp), intent(in) :: observation_uniform_draws(:, :) !! Observation uniform draws.
      type(bsts_prediction_t) :: out

      out = non_gaussian_trend_seasonal_predict_draws(fit, future_exposure, &
         state_normal_draws, observation_uniform_draws, 2)
   end function bsts_poisson_trend_seasonal_predict_draws

   pure function bsts_logit_regression_draws(successes, trials, predictors, &
      slab_mean, slab_variance, prior_inclusion, initial_inclusion, &
      maximum_model_size, maximum_flips, initial_mean, initial_variance, &
      state_variance, state_prior_shape, state_prior_rate, state_proposal_scale, &
      coefficient_proposal_scale, burn, state_normals, state_uniforms, &
      coefficient_normals, coefficient_uniforms, birth_normals, birth_uniforms, &
      gamma_draws) result(out)
      !! Sample sparse binomial-logit regression with a local-level state.
      real(dp), intent(in) :: successes(:) !! Successes.
      real(dp), intent(in) :: trials(:) !! Trials.
      real(dp), intent(in) :: predictors(:, :) !! Predictor matrix.
      real(dp), intent(in) :: slab_mean(:) !! Slab mean.
      real(dp), intent(in) :: slab_variance(:) !! Slab variance.
      real(dp), intent(in) :: prior_inclusion(:) !! Prior inclusion.
      logical, intent(in) :: initial_inclusion(:) !! Initial inclusion.
      integer, intent(in) :: maximum_model_size !! Maximum model size.
      integer, intent(in) :: maximum_flips !! Maximum flips.
      integer, intent(in) :: burn !! Number of initial simulation draws to discard.
      real(dp), intent(in) :: initial_mean !! Initial state mean.
      real(dp), intent(in) :: initial_variance !! Initial variance.
      real(dp), intent(in) :: state_variance !! State variance.
      real(dp), intent(in) :: state_prior_shape !! State prior shape.
      real(dp), intent(in) :: state_prior_rate !! State prior rate.
      real(dp), intent(in) :: state_proposal_scale !! State proposal scale.
      real(dp), intent(in) :: coefficient_proposal_scale !! Coefficient proposal scale.
      real(dp), intent(in) :: state_normals(:, :) !! State normals.
      real(dp), intent(in) :: state_uniforms(:, :) !! State uniforms.
      real(dp), intent(in) :: coefficient_normals(:, :) !! Coefficient normals.
      real(dp), intent(in) :: coefficient_uniforms(:, :) !! Coefficient uniforms.
      real(dp), intent(in) :: birth_normals(:, :) !! Birth normals.
      real(dp), intent(in) :: birth_uniforms(:, :) !! Birth uniforms.
      real(dp), intent(in) :: gamma_draws(:) !! Gamma simulation draws.
      type(bsts_non_gaussian_t) :: out

      out = non_gaussian_regression_draws(successes, trials, predictors, 1, &
         slab_mean, slab_variance, prior_inclusion, initial_inclusion, &
         maximum_model_size, maximum_flips, initial_mean, initial_variance, &
         state_variance, state_prior_shape, state_prior_rate, &
         state_proposal_scale, coefficient_proposal_scale, burn, state_normals, &
         state_uniforms, coefficient_normals, coefficient_uniforms, &
         birth_normals, birth_uniforms, gamma_draws)
   end function bsts_logit_regression_draws

   pure function bsts_poisson_regression_draws(counts, exposure, predictors, &
      slab_mean, slab_variance, prior_inclusion, initial_inclusion, &
      maximum_model_size, maximum_flips, initial_mean, initial_variance, &
      state_variance, state_prior_shape, state_prior_rate, state_proposal_scale, &
      coefficient_proposal_scale, burn, state_normals, state_uniforms, &
      coefficient_normals, coefficient_uniforms, birth_normals, birth_uniforms, &
      gamma_draws) result(out)
      !! Sample sparse Poisson regression with exposure and a local level.
      real(dp), intent(in) :: counts(:) !! Counts.
      real(dp), intent(in) :: exposure(:) !! Exposure.
      real(dp), intent(in) :: predictors(:, :) !! Predictor matrix.
      real(dp), intent(in) :: slab_mean(:) !! Slab mean.
      real(dp), intent(in) :: slab_variance(:) !! Slab variance.
      real(dp), intent(in) :: prior_inclusion(:) !! Prior inclusion.
      logical, intent(in) :: initial_inclusion(:) !! Initial inclusion.
      integer, intent(in) :: maximum_model_size !! Maximum model size.
      integer, intent(in) :: maximum_flips !! Maximum flips.
      integer, intent(in) :: burn !! Number of initial simulation draws to discard.
      real(dp), intent(in) :: initial_mean !! Initial state mean.
      real(dp), intent(in) :: initial_variance !! Initial variance.
      real(dp), intent(in) :: state_variance !! State variance.
      real(dp), intent(in) :: state_prior_shape !! State prior shape.
      real(dp), intent(in) :: state_prior_rate !! State prior rate.
      real(dp), intent(in) :: state_proposal_scale !! State proposal scale.
      real(dp), intent(in) :: coefficient_proposal_scale !! Coefficient proposal scale.
      real(dp), intent(in) :: state_normals(:, :) !! State normals.
      real(dp), intent(in) :: state_uniforms(:, :) !! State uniforms.
      real(dp), intent(in) :: coefficient_normals(:, :) !! Coefficient normals.
      real(dp), intent(in) :: coefficient_uniforms(:, :) !! Coefficient uniforms.
      real(dp), intent(in) :: birth_normals(:, :) !! Birth normals.
      real(dp), intent(in) :: birth_uniforms(:, :) !! Birth uniforms.
      real(dp), intent(in) :: gamma_draws(:) !! Gamma simulation draws.
      type(bsts_non_gaussian_t) :: out

      out = non_gaussian_regression_draws(counts, exposure, predictors, 2, &
         slab_mean, slab_variance, prior_inclusion, initial_inclusion, &
         maximum_model_size, maximum_flips, initial_mean, initial_variance, &
         state_variance, state_prior_shape, state_prior_rate, &
         state_proposal_scale, coefficient_proposal_scale, burn, state_normals, &
         state_uniforms, coefficient_normals, coefficient_uniforms, &
         birth_normals, birth_uniforms, gamma_draws)
   end function bsts_poisson_regression_draws

   function bsts_logit_regression(successes, trials, predictors, iterations, &
      burn, prior_inclusion_probability, initial_inclusion, &
      maximum_model_size, maximum_flips) result(out)
      !! Sample sparse binomial-logit regression using random draws.
      real(dp), intent(in) :: successes(:) !! Successes.
      real(dp), intent(in) :: trials(:) !! Trials.
      real(dp), intent(in) :: predictors(:, :) !! Predictor matrix.
      integer, intent(in) :: iterations !! Number of algorithm iterations.
      integer, intent(in), optional :: burn !! Number of initial simulation draws to discard.
      integer, intent(in), optional :: maximum_model_size !! Maximum model size.
      integer, intent(in), optional :: maximum_flips !! Maximum flips.
      real(dp), intent(in), optional :: prior_inclusion_probability(:) !! Prior inclusion probability.
      logical, intent(in), optional :: initial_inclusion(:) !! Initial inclusion.
      type(bsts_non_gaussian_t) :: out

      out = random_non_gaussian_regression(successes, trials, predictors, 1, &
         iterations, burn, prior_inclusion_probability, initial_inclusion, &
         maximum_model_size, maximum_flips)
   end function bsts_logit_regression

   function bsts_poisson_regression(counts, exposure, predictors, iterations, &
      burn, prior_inclusion_probability, initial_inclusion, &
      maximum_model_size, maximum_flips) result(out)
      !! Sample sparse Poisson regression with observation exposures.
      real(dp), intent(in) :: counts(:) !! Counts.
      real(dp), intent(in) :: exposure(:) !! Exposure.
      real(dp), intent(in) :: predictors(:, :) !! Predictor matrix.
      integer, intent(in) :: iterations !! Number of algorithm iterations.
      integer, intent(in), optional :: burn !! Number of initial simulation draws to discard.
      integer, intent(in), optional :: maximum_model_size !! Maximum model size.
      integer, intent(in), optional :: maximum_flips !! Maximum flips.
      real(dp), intent(in), optional :: prior_inclusion_probability(:) !! Prior inclusion probability.
      logical, intent(in), optional :: initial_inclusion(:) !! Initial inclusion.
      type(bsts_non_gaussian_t) :: out

      out = random_non_gaussian_regression(counts, exposure, predictors, 2, &
         iterations, burn, prior_inclusion_probability, initial_inclusion, &
         maximum_model_size, maximum_flips)
   end function bsts_poisson_regression

   pure function bsts_logit_regression_predict_draws(fit, future_trials, &
      future_predictors, state_normal_draws, observation_uniform_draws) &
      result(out)
      !! Forecast sparse binomial-logit regression from supplied draws.
      type(bsts_non_gaussian_t), intent(in) :: fit !! Previously fitted model.
      real(dp), intent(in) :: future_trials(:) !! Future trials.
      real(dp), intent(in) :: future_predictors(:, :) !! Future predictors.
      real(dp), intent(in) :: state_normal_draws(:, :) !! Independent standard-normal state draws.
      real(dp), intent(in) :: observation_uniform_draws(:, :) !! Observation uniform draws.
      type(bsts_prediction_t) :: out

      out = non_gaussian_regression_predict_draws(fit, future_trials, &
         future_predictors, state_normal_draws, observation_uniform_draws, 1)
   end function bsts_logit_regression_predict_draws

   pure function bsts_poisson_regression_predict_draws(fit, future_exposure, &
      future_predictors, state_normal_draws, observation_uniform_draws) &
      result(out)
      !! Forecast sparse Poisson regression with future exposures.
      type(bsts_non_gaussian_t), intent(in) :: fit !! Previously fitted model.
      real(dp), intent(in) :: future_exposure(:) !! Future exposure.
      real(dp), intent(in) :: future_predictors(:, :) !! Future predictors.
      real(dp), intent(in) :: state_normal_draws(:, :) !! Independent standard-normal state draws.
      real(dp), intent(in) :: observation_uniform_draws(:, :) !! Observation uniform draws.
      type(bsts_prediction_t) :: out

      out = non_gaussian_regression_predict_draws(fit, future_exposure, &
         future_predictors, state_normal_draws, observation_uniform_draws, 2)
   end function bsts_poisson_regression_predict_draws

   pure function bsts_logit_local_level_draws(successes, trials, &
      initial_mean, initial_variance, state_variance, state_prior_shape, &
      state_prior_rate, proposal_scale, burn, proposal_normal_draws, &
      acceptance_uniform_draws, gamma_draws) result(out)
      !! Sample a binomial-logit local-level posterior from supplied draws.
      real(dp), intent(in) :: successes(:) !! Successes.
      real(dp), intent(in) :: trials(:) !! Trials.
      real(dp), intent(in) :: initial_mean !! Initial state mean.
      real(dp), intent(in) :: initial_variance !! Initial variance.
      real(dp), intent(in) :: state_variance !! State variance.
      real(dp), intent(in) :: state_prior_shape !! State prior shape.
      real(dp), intent(in) :: state_prior_rate !! State prior rate.
      real(dp), intent(in) :: proposal_scale !! Proposal scale.
      integer, intent(in) :: burn !! Number of initial simulation draws to discard.
      real(dp), intent(in) :: proposal_normal_draws(:, :) !! Proposal normal draws.
      real(dp), intent(in) :: acceptance_uniform_draws(:, :) !! Acceptance uniform simulation draws.
      real(dp), intent(in) :: gamma_draws(:) !! Gamma simulation draws.
      type(bsts_non_gaussian_t) :: out

      out = non_gaussian_local_level_draws(successes, trials, 1, &
         initial_mean, initial_variance, state_variance, state_prior_shape, &
         state_prior_rate, proposal_scale, burn, proposal_normal_draws, &
         acceptance_uniform_draws, gamma_draws)
   end function bsts_logit_local_level_draws

   pure function bsts_poisson_local_level_draws(counts, exposure, &
      initial_mean, initial_variance, state_variance, state_prior_shape, &
      state_prior_rate, proposal_scale, burn, proposal_normal_draws, &
      acceptance_uniform_draws, gamma_draws) result(out)
      !! Sample a Poisson local-level posterior with exposure offsets.
      real(dp), intent(in) :: counts(:) !! Counts.
      real(dp), intent(in) :: exposure(:) !! Exposure.
      real(dp), intent(in) :: initial_mean !! Initial state mean.
      real(dp), intent(in) :: initial_variance !! Initial variance.
      real(dp), intent(in) :: state_variance !! State variance.
      real(dp), intent(in) :: state_prior_shape !! State prior shape.
      real(dp), intent(in) :: state_prior_rate !! State prior rate.
      real(dp), intent(in) :: proposal_scale !! Proposal scale.
      integer, intent(in) :: burn !! Number of initial simulation draws to discard.
      real(dp), intent(in) :: proposal_normal_draws(:, :) !! Proposal normal draws.
      real(dp), intent(in) :: acceptance_uniform_draws(:, :) !! Acceptance uniform simulation draws.
      real(dp), intent(in) :: gamma_draws(:) !! Gamma simulation draws.
      type(bsts_non_gaussian_t) :: out

      out = non_gaussian_local_level_draws(counts, exposure, 2, initial_mean, &
         initial_variance, state_variance, state_prior_shape, state_prior_rate, &
         proposal_scale, burn, proposal_normal_draws, &
         acceptance_uniform_draws, gamma_draws)
   end function bsts_poisson_local_level_draws

   function bsts_logit_local_level(successes, trials, iterations, burn, &
      proposal_scale) result(out)
      !! Sample a binomial-logit local-level posterior using random draws.
      real(dp), intent(in) :: successes(:) !! Successes.
      real(dp), intent(in) :: trials(:) !! Trials.
      integer, intent(in) :: iterations !! Number of algorithm iterations.
      integer, intent(in), optional :: burn !! Number of initial simulation draws to discard.
      real(dp), intent(in), optional :: proposal_scale !! Proposal scale.
      type(bsts_non_gaussian_t) :: out
      real(dp), allocatable :: normals(:, :), uniforms(:, :), gammas(:)
      real(dp) :: scale
      integer :: discarded, iteration, time

      discarded = 0
      if (present(burn)) discarded = burn
      scale = 0.35_dp
      if (present(proposal_scale)) scale = proposal_scale
      if (iterations < 1) then
         out%info = 1
         return
      end if
      allocate(normals(size(successes), iterations), &
         uniforms(size(successes), iterations), gammas(iterations))
      do iteration = 1, iterations
         do time = 1, size(successes)
            normals(time, iteration) = random_standard_normal()
            uniforms(time, iteration) = max(tiny(1.0_dp), &
               min(1.0_dp - epsilon(1.0_dp), random_uniform()))
         end do
         gammas(iteration) = random_gamma(2.0_dp + &
            0.5_dp*real(max(0, size(successes) - 1), dp))
      end do
      out = bsts_logit_local_level_draws(successes, trials, 0.0_dp, 10.0_dp, &
         0.05_dp, 2.0_dp, 0.05_dp, scale, discarded, normals, uniforms, gammas)
   end function bsts_logit_local_level

   function bsts_poisson_local_level(counts, exposure, iterations, burn, &
      proposal_scale) result(out)
      !! Sample a Poisson local-level posterior using random draws.
      real(dp), intent(in) :: counts(:) !! Counts.
      real(dp), intent(in) :: exposure(:) !! Exposure.
      integer, intent(in) :: iterations !! Number of algorithm iterations.
      integer, intent(in), optional :: burn !! Number of initial simulation draws to discard.
      real(dp), intent(in), optional :: proposal_scale !! Proposal scale.
      type(bsts_non_gaussian_t) :: out
      real(dp), allocatable :: normals(:, :), uniforms(:, :), gammas(:)
      real(dp) :: scale
      integer :: discarded, iteration, time

      discarded = 0
      if (present(burn)) discarded = burn
      scale = 0.25_dp
      if (present(proposal_scale)) scale = proposal_scale
      if (iterations < 1) then
         out%info = 1
         return
      end if
      allocate(normals(size(counts), iterations), &
         uniforms(size(counts), iterations), gammas(iterations))
      do iteration = 1, iterations
         do time = 1, size(counts)
            normals(time, iteration) = random_standard_normal()
            uniforms(time, iteration) = max(tiny(1.0_dp), &
               min(1.0_dp - epsilon(1.0_dp), random_uniform()))
         end do
         gammas(iteration) = random_gamma(2.0_dp + &
            0.5_dp*real(max(0, size(counts) - 1), dp))
      end do
      out = bsts_poisson_local_level_draws(counts, exposure, 0.0_dp, 10.0_dp, &
         0.05_dp, 2.0_dp, 0.05_dp, scale, discarded, normals, uniforms, gammas)
   end function bsts_poisson_local_level

   pure function bsts_logit_predict_draws(fit, future_trials, &
      state_normal_draws, observation_uniform_draws) result(out)
      !! Forecast binomial observations from supplied random draws.
      type(bsts_non_gaussian_t), intent(in) :: fit !! Previously fitted model.
      real(dp), intent(in) :: future_trials(:) !! Future trials.
      real(dp), intent(in) :: state_normal_draws(:, :) !! Independent standard-normal state draws.
      real(dp), intent(in) :: observation_uniform_draws(:, :) !! Observation uniform draws.
      type(bsts_prediction_t) :: out

      out = non_gaussian_predict_draws(fit, future_trials, &
         state_normal_draws, observation_uniform_draws, 1)
   end function bsts_logit_predict_draws

   pure function bsts_poisson_predict_draws(fit, future_exposure, &
      state_normal_draws, observation_uniform_draws) result(out)
      !! Forecast Poisson observations with supplied future exposures.
      type(bsts_non_gaussian_t), intent(in) :: fit !! Previously fitted model.
      real(dp), intent(in) :: future_exposure(:) !! Future exposure.
      real(dp), intent(in) :: state_normal_draws(:, :) !! Independent standard-normal state draws.
      real(dp), intent(in) :: observation_uniform_draws(:, :) !! Observation uniform draws.
      type(bsts_prediction_t) :: out

      out = non_gaussian_predict_draws(fit, future_exposure, &
         state_normal_draws, observation_uniform_draws, 2)
   end function bsts_poisson_predict_draws

   function bsts_logit_predict(fit, future_trials) result(out)
      !! Forecast binomial observations using the shared random stream.
      type(bsts_non_gaussian_t), intent(in) :: fit !! Previously fitted model.
      real(dp), intent(in) :: future_trials(:) !! Future trials.
      type(bsts_prediction_t) :: out
      real(dp), allocatable :: normals(:, :), uniforms(:, :)

      call fill_non_gaussian_forecast_draws(fit, size(future_trials), &
         normals, uniforms)
      if (.not. allocated(normals)) then
         out%info = 1
         return
      end if
      out = bsts_logit_predict_draws(fit, future_trials, normals, uniforms)
   end function bsts_logit_predict

   function bsts_poisson_predict(fit, future_exposure) result(out)
      !! Forecast Poisson observations using future exposure offsets.
      type(bsts_non_gaussian_t), intent(in) :: fit !! Previously fitted model.
      real(dp), intent(in) :: future_exposure(:) !! Future exposure.
      type(bsts_prediction_t) :: out
      real(dp), allocatable :: normals(:, :), uniforms(:, :)

      call fill_non_gaussian_forecast_draws(fit, size(future_exposure), &
         normals, uniforms)
      if (.not. allocated(normals)) then
         out%info = 1
         return
      end if
      out = bsts_poisson_predict_draws(fit, future_exposure, normals, uniforms)
   end function bsts_poisson_predict

   function random_non_gaussian_trend_seasonal(response, scale_data, family, &
      nseasons, iterations, season_duration, burn) result(out)
      !! Generate random streams for non-Gaussian structural states.
      real(dp), intent(in) :: response(:) !! Response observations.
      real(dp), intent(in) :: scale_data(:) !! Scale data.
      integer, intent(in) :: family !! Family.
      integer, intent(in) :: nseasons !! Nseasons.
      integer, intent(in) :: iterations !! Number of algorithm iterations.
      integer, intent(in), optional :: season_duration !! Season duration.
      integer, intent(in), optional :: burn !! Number of initial simulation draws to discard.
      type(bsts_non_gaussian_t) :: out
      real(dp), allocatable :: normals(:, :, :), uniforms(:), gammas(:, :)
      integer :: duration, discarded, dimension, iteration, time, component

      duration = 1
      if (present(season_duration)) duration = season_duration
      discarded = 0
      if (present(burn)) discarded = burn
      dimension = nseasons + 1
      if (dimension < 3 .or. iterations < 1) then
         out%info = 1
         return
      end if
      allocate(normals(dimension, size(response), iterations), &
         uniforms(iterations), gammas(3, iterations))
      do iteration = 1, iterations
         do time = 1, size(response)
            do component = 1, dimension
               normals(component, time, iteration) = random_standard_normal()
            end do
         end do
         uniforms(iteration) = max(tiny(1.0_dp), &
            min(1.0_dp - epsilon(1.0_dp), random_uniform()))
         do component = 1, 2
            gammas(component, iteration) = random_gamma(2.0_dp + &
               0.5_dp*real(max(0, size(response) - 1), dp))
         end do
         gammas(3, iteration) = random_gamma(2.0_dp + &
            0.5_dp*real(max(0, (size(response) - 1)/duration), dp))
      end do
      if (family == 1) then
         out = bsts_logit_trend_seasonal_draws(response, scale_data, nseasons, &
            duration, [0.0_dp, 0.0_dp], &
            reshape([10.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2]), 1.0_dp, &
            [0.05_dp, 0.01_dp, 0.01_dp], [2.0_dp, 2.0_dp, 2.0_dp], &
            [0.05_dp, 0.01_dp, 0.01_dp], discarded, normals, uniforms, gammas)
      else
         out = bsts_poisson_trend_seasonal_draws(response, scale_data, &
            nseasons, duration, [0.0_dp, 0.0_dp], &
            reshape([10.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2]), 1.0_dp, &
            [0.05_dp, 0.01_dp, 0.01_dp], [2.0_dp, 2.0_dp, 2.0_dp], &
            [0.05_dp, 0.01_dp, 0.01_dp], discarded, normals, uniforms, gammas)
      end if
   end function random_non_gaussian_trend_seasonal

   pure function non_gaussian_trend_seasonal_draws(response, scale_data, &
      family, nseasons, duration, initial_trend_mean, initial_trend_covariance, &
      seasonal_initial_variance, variance0, prior_shape, prior_rate, burn, &
      normals, uniforms, gammas) result(out)
      !! Sample an exact-likelihood trend-seasonal structural posterior.
      real(dp), intent(in) :: response(:) !! Response observations.
      real(dp), intent(in) :: scale_data(:) !! Scale data.
      integer, intent(in) :: family !! Family.
      integer, intent(in) :: nseasons !! Nseasons.
      integer, intent(in) :: duration !! Duration.
      integer, intent(in) :: burn !! Number of initial simulation draws to discard.
      real(dp), intent(in) :: initial_trend_mean(2) !! Initial trend mean.
      real(dp), intent(in) :: initial_trend_covariance(2, 2) !! Initial trend covariance.
      real(dp), intent(in) :: seasonal_initial_variance !! Seasonal initial variance.
      real(dp), intent(in) :: variance0(3) !! Variance0.
      real(dp), intent(in) :: prior_shape(3) !! Prior shape.
      real(dp), intent(in) :: prior_rate(3) !! Prior rate.
      real(dp), intent(in) :: normals(:, :, :) !! Independent standard-normal draws.
      real(dp), intent(in) :: uniforms(:) !! Uniforms.
      real(dp), intent(in) :: gammas(:, :) !! Gammas.
      type(bsts_non_gaussian_t) :: out
      real(dp), allocatable :: current(:, :), proposal(:, :)
      real(dp), allocatable :: contribution(:), proposed_contribution(:)
      real(dp), allocatable :: initial_draw(:)
      real(dp) :: variance(3), old_likelihood, new_likelihood, residual_sum
      integer :: dimension, observations, iterations, iteration, time, status
      integer :: component

      observations = size(response)
      dimension = nseasons + 1
      iterations = size(normals, 3)
      if (observations < 2 .or. nseasons < 2 .or. duration < 1 .or. &
         burn < 0 .or. burn >= iterations .or. &
         size(scale_data) /= observations .or. seasonal_initial_variance <= 0.0_dp .or. &
         any(variance0 <= 0.0_dp) .or. any(prior_shape <= 0.0_dp) .or. &
         any(prior_rate <= 0.0_dp) .or. &
         any(shape(normals) /= [dimension, observations, iterations]) .or. &
         size(uniforms) /= iterations .or. &
         any(shape(gammas) /= [3, iterations]) .or. any(gammas <= 0.0_dp) .or. &
         any(uniforms <= 0.0_dp) .or. any(uniforms >= 1.0_dp) .or. &
         .not. all(ieee_is_finite(scale_data)) .or. &
         .not. all(ieee_is_finite(normals)) .or. &
         (family /= 1 .and. family /= 2)) then
         out%info = 1
         return
      end if
      if ((family == 1 .and. (any(scale_data < 0.0_dp) .or. &
         any(response < 0.0_dp) .or. any(response > scale_data))) .or. &
         (family == 2 .and. (any(scale_data <= 0.0_dp) .or. &
         any(response < 0.0_dp)))) then
         out%info = 1
         return
      end if
      allocate(current(dimension, observations), proposal(dimension, observations), &
         contribution(observations), proposed_contribution(observations), &
         initial_draw(2), out%structural_state(dimension, observations, iterations), &
         out%state(observations, iterations), out%fitted_mean(observations, iterations), &
         out%component_variance(3, iterations), out%state_variance(iterations))
      current = 0.0_dp
      contribution = 0.0_dp
      variance = variance0
      out%burn = burn
      out%family = family
      out%nseasons = nseasons
      out%season_duration = duration
      do iteration = 1, iterations
         call multivariate_normal_from_psd(initial_trend_mean, &
            initial_trend_covariance, normals(:2, 1, iteration), initial_draw, status)
         if (status /= 0) then
            out%info = 2
            return
         end if
         proposal = 0.0_dp
         proposal(:2, 1) = initial_draw
         proposal(3:, 1) = sqrt(seasonal_initial_variance)* &
            normals(3:, 1, iteration)
         do time = 2, observations
            proposal(1, time) = proposal(1, time - 1) + &
               proposal(2, time - 1) + sqrt(variance(1))*normals(1, time, iteration)
            proposal(2, time) = proposal(2, time - 1) + &
               sqrt(variance(2))*normals(2, time, iteration)
            proposal(3:, time) = proposal(3:, time - 1)
            if (modulo(time - 1, duration) == 0) then
               proposal(3, time) = -sum(proposal(3:, time - 1)) + &
                  sqrt(variance(3))*normals(3, time, iteration)
               do component = 4, dimension
                  proposal(component, time) = proposal(component - 1, time - 1)
               end do
            end if
         end do
         proposed_contribution = proposal(1, :) + proposal(3, :)
         if (iteration == 1) then
            current = proposal
            contribution = proposed_contribution
         else
            old_likelihood = non_gaussian_log_likelihood(response, scale_data, &
               contribution, family)
            new_likelihood = non_gaussian_log_likelihood(response, scale_data, &
               proposed_contribution, family)
            if (log(uniforms(iteration)) < new_likelihood - old_likelihood) then
               current = proposal
               contribution = proposed_contribution
            end if
         end if
         residual_sum = sum((current(1, 2:) - current(1, :observations - 1) - &
            current(2, :observations - 1))**2)
         variance(1) = (prior_rate(1) + 0.5_dp*residual_sum)/gammas(1, iteration)
         residual_sum = sum((current(2, 2:) - current(2, :observations - 1))**2)
         variance(2) = (prior_rate(2) + 0.5_dp*residual_sum)/gammas(2, iteration)
         residual_sum = 0.0_dp
         do time = 2, observations
            if (modulo(time - 1, duration) == 0) residual_sum = residual_sum + &
               (current(3, time) + sum(current(3:, time - 1)))**2
         end do
         variance(3) = (prior_rate(3) + 0.5_dp*residual_sum)/gammas(3, iteration)
         out%structural_state(:, :, iteration) = current
         out%state(:, iteration) = contribution
         out%component_variance(:, iteration) = variance
         out%state_variance(iteration) = variance(1)
         if (family == 1) then
            out%fitted_mean(:, iteration) = scale_data*logistic_log_odds(contribution)
         else
            out%fitted_mean(:, iteration) = scale_data*exp(min(contribution, 700.0_dp))
         end if
      end do
   end function non_gaussian_trend_seasonal_draws

   function random_non_gaussian_regression(response, scale_data, predictors, &
      family, iterations, burn, prior_inclusion_probability, initial_inclusion, &
      maximum_model_size, maximum_flips) result(out)
      !! Generate random streams for sparse non-Gaussian regression.
      real(dp), intent(in) :: response(:) !! Response observations.
      real(dp), intent(in) :: scale_data(:) !! Scale data.
      real(dp), intent(in) :: predictors(:, :) !! Predictor matrix.
      integer, intent(in) :: family !! Family.
      integer, intent(in) :: iterations !! Number of algorithm iterations.
      integer, intent(in), optional :: burn !! Number of initial simulation draws to discard.
      integer, intent(in), optional :: maximum_model_size !! Maximum model size.
      integer, intent(in), optional :: maximum_flips !! Maximum flips.
      real(dp), intent(in), optional :: prior_inclusion_probability(:) !! Prior inclusion probability.
      logical, intent(in), optional :: initial_inclusion(:) !! Initial inclusion.
      type(bsts_non_gaussian_t) :: out
      real(dp), allocatable :: prior(:), slab_mean(:), slab_variance(:)
      real(dp), allocatable :: state_normals(:, :), state_uniforms(:, :)
      real(dp), allocatable :: coefficient_normals(:, :)
      real(dp), allocatable :: coefficient_uniforms(:, :), birth_normals(:, :)
      real(dp), allocatable :: birth_uniforms(:, :), gammas(:)
      logical, allocatable :: included(:)
      integer :: variables, discarded, model_limit, flip_limit
      integer :: iteration, time, variable

      variables = size(predictors, 2)
      discarded = 0
      if (present(burn)) discarded = burn
      model_limit = variables
      if (present(maximum_model_size)) model_limit = maximum_model_size
      flip_limit = -1
      if (present(maximum_flips)) flip_limit = maximum_flips
      if (iterations < 1 .or. variables < 1) then
         out%info = 1
         return
      end if
      allocate(prior(variables), slab_mean(variables), slab_variance(variables), &
         included(variables), state_normals(size(response), iterations), &
         state_uniforms(size(response), iterations), &
         coefficient_normals(variables, iterations), &
         coefficient_uniforms(variables, iterations), &
         birth_normals(variables, iterations), &
         birth_uniforms(variables, iterations), gammas(iterations))
      prior = 0.5_dp
      if (present(prior_inclusion_probability)) then
         if (size(prior_inclusion_probability) /= variables) then
            out%info = 1
            return
         end if
         prior = prior_inclusion_probability
      end if
      included = prior > 0.5_dp
      if (present(initial_inclusion)) then
         if (size(initial_inclusion) /= variables) then
            out%info = 1
            return
         end if
         included = initial_inclusion
      end if
      slab_mean = 0.0_dp
      slab_variance = 10.0_dp
      do iteration = 1, iterations
         do time = 1, size(response)
            state_normals(time, iteration) = random_standard_normal()
            state_uniforms(time, iteration) = max(tiny(1.0_dp), &
               min(1.0_dp - epsilon(1.0_dp), random_uniform()))
         end do
         do variable = 1, variables
            coefficient_normals(variable, iteration) = random_standard_normal()
            birth_normals(variable, iteration) = random_standard_normal()
            coefficient_uniforms(variable, iteration) = max(tiny(1.0_dp), &
               min(1.0_dp - epsilon(1.0_dp), random_uniform()))
            birth_uniforms(variable, iteration) = max(tiny(1.0_dp), &
               min(1.0_dp - epsilon(1.0_dp), random_uniform()))
         end do
         gammas(iteration) = random_gamma(2.0_dp + &
            0.5_dp*real(size(response) - 1, dp))
      end do
      if (family == 1) then
         out = bsts_logit_regression_draws(response, scale_data, predictors, &
            slab_mean, slab_variance, prior, included, model_limit, flip_limit, &
            0.0_dp, 10.0_dp, 0.05_dp, 2.0_dp, 0.05_dp, 0.3_dp, 0.2_dp, &
            discarded, state_normals, state_uniforms, coefficient_normals, &
            coefficient_uniforms, birth_normals, birth_uniforms, gammas)
      else
         out = bsts_poisson_regression_draws(response, scale_data, predictors, &
            slab_mean, slab_variance, prior, included, model_limit, flip_limit, &
            0.0_dp, 10.0_dp, 0.05_dp, 2.0_dp, 0.05_dp, 0.25_dp, 0.15_dp, &
            discarded, state_normals, state_uniforms, coefficient_normals, &
            coefficient_uniforms, birth_normals, birth_uniforms, gammas)
      end if
   end function random_non_gaussian_regression

   pure function non_gaussian_regression_draws(response, scale_data, predictors, &
      family, slab_mean, slab_variance, prior_inclusion, initial_inclusion, &
      maximum_model_size, maximum_flips, initial_mean, initial_variance, &
      variance0, prior_shape, prior_rate, state_proposal_scale, &
      coefficient_proposal_scale, burn, state_normals, state_uniforms, &
      coefficient_normals, coefficient_uniforms, birth_normals, birth_uniforms, &
      gammas) result(out)
      !! Run sparse non-Gaussian regression and local-level MCMC updates.
      real(dp), intent(in) :: response(:) !! Response observations.
      real(dp), intent(in) :: scale_data(:) !! Scale data.
      real(dp), intent(in) :: predictors(:, :) !! Predictor matrix.
      integer, intent(in) :: family !! Family.
      integer, intent(in) :: maximum_model_size !! Maximum model size.
      integer, intent(in) :: maximum_flips !! Maximum flips.
      integer, intent(in) :: burn !! Number of initial simulation draws to discard.
      real(dp), intent(in) :: slab_mean(:) !! Slab mean.
      real(dp), intent(in) :: slab_variance(:) !! Slab variance.
      real(dp), intent(in) :: prior_inclusion(:) !! Prior inclusion.
      logical, intent(in) :: initial_inclusion(:) !! Initial inclusion.
      real(dp), intent(in) :: initial_mean !! Initial state mean.
      real(dp), intent(in) :: initial_variance !! Initial variance.
      real(dp), intent(in) :: variance0 !! Variance0.
      real(dp), intent(in) :: prior_shape !! Prior shape.
      real(dp), intent(in) :: prior_rate !! Prior rate.
      real(dp), intent(in) :: state_proposal_scale !! State proposal scale.
      real(dp), intent(in) :: coefficient_proposal_scale !! Coefficient proposal scale.
      real(dp), intent(in) :: state_normals(:, :) !! State normals.
      real(dp), intent(in) :: state_uniforms(:, :) !! State uniforms.
      real(dp), intent(in) :: coefficient_normals(:, :) !! Coefficient normals.
      real(dp), intent(in) :: coefficient_uniforms(:, :) !! Coefficient uniforms.
      real(dp), intent(in) :: birth_normals(:, :) !! Birth normals.
      real(dp), intent(in) :: birth_uniforms(:, :) !! Birth uniforms.
      real(dp), intent(in) :: gammas(:) !! Gammas.
      type(bsts_non_gaussian_t) :: out
      real(dp), allocatable :: state(:), coefficient(:), regression(:)
      logical, allocatable :: included(:)
      real(dp) :: proposal, old_density, new_density, variance, difference_sum
      real(dp) :: log_ratio
      integer :: observations, variables, iterations, iteration, time, variable
      integer :: flips, retained

      observations = size(response)
      variables = size(predictors, 2)
      iterations = size(state_normals, 2)
      if (observations < 2 .or. variables < 1 .or. &
         any(shape(predictors) /= [observations, variables]) .or. &
         size(scale_data) /= observations .or. size(slab_mean) /= variables .or. &
         size(slab_variance) /= variables .or. &
         size(prior_inclusion) /= variables .or. &
         size(initial_inclusion) /= variables .or. &
         any(slab_variance <= 0.0_dp) .or. any(prior_inclusion < 0.0_dp) .or. &
         any(prior_inclusion > 1.0_dp) .or. maximum_model_size < 0 .or. &
         maximum_model_size > variables .or. maximum_flips < -1 .or. &
         burn < 0 .or. burn >= iterations .or. initial_variance <= 0.0_dp .or. &
         variance0 <= 0.0_dp .or. prior_shape <= 0.0_dp .or. &
         prior_rate <= 0.0_dp .or. state_proposal_scale <= 0.0_dp .or. &
         coefficient_proposal_scale <= 0.0_dp .or. &
         any(shape(state_normals) /= [observations, iterations]) .or. &
         any(shape(state_uniforms) /= [observations, iterations]) .or. &
         any(shape(coefficient_normals) /= [variables, iterations]) .or. &
         any(shape(coefficient_uniforms) /= [variables, iterations]) .or. &
         any(shape(birth_normals) /= [variables, iterations]) .or. &
         any(shape(birth_uniforms) /= [variables, iterations]) .or. &
         size(gammas) /= iterations .or. any(gammas <= 0.0_dp) .or. &
         any(state_uniforms <= 0.0_dp) .or. any(state_uniforms >= 1.0_dp) .or. &
         any(coefficient_uniforms <= 0.0_dp) .or. &
         any(coefficient_uniforms >= 1.0_dp) .or. &
         any(birth_uniforms <= 0.0_dp) .or. any(birth_uniforms >= 1.0_dp) .or. &
         .not. all(ieee_is_finite(predictors)) .or. &
         .not. all(ieee_is_finite(scale_data)) .or. &
         .not. all(ieee_is_finite(slab_mean)) .or. &
         .not. all(ieee_is_finite(slab_variance)) .or. &
         .not. all(ieee_is_finite(prior_inclusion)) .or. &
         .not. all(ieee_is_finite(state_normals)) .or. &
         .not. all(ieee_is_finite(coefficient_normals)) .or. &
         .not. all(ieee_is_finite(birth_normals)) .or. &
         (family /= 1 .and. family /= 2)) then
         out%info = 1
         return
      end if
      if ((family == 1 .and. (any(scale_data < 0.0_dp) .or. &
         any(response < 0.0_dp) .or. any(response > scale_data))) .or. &
         (family == 2 .and. (any(scale_data <= 0.0_dp) .or. &
         any(response < 0.0_dp)))) then
         out%info = 1
         return
      end if
      allocate(state(observations), coefficient(variables), &
         regression(observations), included(variables), &
         out%state(observations, iterations), &
         out%state_variance(iterations), out%fitted_mean(observations, iterations), &
         out%coefficients(variables, iterations), &
         out%regression_contribution(observations, iterations), &
         out%included(variables, iterations), out%inclusion_probability(variables))
      state = initial_mean
      coefficient = 0.0_dp
      included = initial_inclusion
      where (prior_inclusion <= 0.0_dp) included = .false.
      where (prior_inclusion >= 1.0_dp) included = .true.
      if (count(included) > maximum_model_size) then
         out%info = 1
         return
      end if
      where (included) coefficient = slab_mean
      regression = matmul(predictors, coefficient)
      variance = variance0
      out%burn = burn
      out%family = family
      do iteration = 1, iterations
         do time = 1, observations
            proposal = state(time) + &
               state_proposal_scale*state_normals(time, iteration)
            old_density = local_level_regression_log_density(state(time), time, &
               state, regression, response, scale_data, family, initial_mean, &
               initial_variance, variance)
            new_density = local_level_regression_log_density(proposal, time, &
               state, regression, response, scale_data, family, initial_mean, &
               initial_variance, variance)
            if (log(state_uniforms(time, iteration)) < &
               new_density - old_density) state(time) = proposal
         end do
         flips = 0
         do variable = 1, variables
            if (prior_inclusion(variable) <= 0.0_dp) then
               included(variable) = .false.
               coefficient(variable) = 0.0_dp
               cycle
            end if
            if (prior_inclusion(variable) >= 1.0_dp) then
               included(variable) = .true.
            else if (maximum_flips < 0 .or. flips < maximum_flips) then
               if (included(variable)) then
                  proposal = 0.0_dp
                  log_ratio = non_gaussian_log_likelihood(response, scale_data, &
                     state + regression + predictors(:, variable)* &
                     (proposal - coefficient(variable)), family) - &
                     non_gaussian_log_likelihood(response, scale_data, &
                     state + regression, family) + &
                     log((1.0_dp - prior_inclusion(variable))/ &
                     prior_inclusion(variable))
               else if (count(included) < maximum_model_size) then
                  proposal = slab_mean(variable) + sqrt(slab_variance(variable))* &
                     birth_normals(variable, iteration)
                  log_ratio = non_gaussian_log_likelihood(response, scale_data, &
                     state + regression + predictors(:, variable)*proposal, &
                     family) - non_gaussian_log_likelihood(response, scale_data, &
                     state + regression, family) + log(prior_inclusion(variable)/ &
                     (1.0_dp - prior_inclusion(variable)))
               else
                  cycle
               end if
               if (log(birth_uniforms(variable, iteration)) < log_ratio) then
                  regression = regression + predictors(:, variable)* &
                     (proposal - coefficient(variable))
                  coefficient(variable) = proposal
                  included(variable) = .not. included(variable)
                  flips = flips + 1
               end if
            end if
            if (included(variable)) then
               proposal = coefficient(variable) + coefficient_proposal_scale* &
                  coefficient_normals(variable, iteration)
               log_ratio = non_gaussian_log_likelihood(response, scale_data, &
                  state + regression + predictors(:, variable)* &
                  (proposal - coefficient(variable)), family) - &
                  non_gaussian_log_likelihood(response, scale_data, &
                  state + regression, family) - 0.5_dp* &
                  ((proposal - slab_mean(variable))**2 - &
                  (coefficient(variable) - slab_mean(variable))**2)/ &
                  slab_variance(variable)
               if (log(coefficient_uniforms(variable, iteration)) < log_ratio) then
                  regression = regression + predictors(:, variable)* &
                     (proposal - coefficient(variable))
                  coefficient(variable) = proposal
               end if
            end if
         end do
         difference_sum = sum((state(2:) - state(:observations - 1))**2)
         variance = (prior_rate + 0.5_dp*difference_sum)/gammas(iteration)
         out%state(:, iteration) = state
         out%state_variance(iteration) = variance
         out%coefficients(:, iteration) = coefficient
         out%regression_contribution(:, iteration) = regression
         out%included(:, iteration) = included
         if (family == 1) then
            out%fitted_mean(:, iteration) = scale_data* &
               logistic_log_odds(state + regression)
         else
            out%fitted_mean(:, iteration) = scale_data* &
               exp(min(state + regression, 700.0_dp))
         end if
      end do
      retained = iterations - burn
      do variable = 1, variables
         out%inclusion_probability(variable) = real(count( &
            out%included(variable, burn + 1:)), dp)/real(retained, dp)
      end do
   end function non_gaussian_regression_draws

   pure real(dp) function local_level_regression_log_density(value, time, state, &
      regression, response, scale_data, family, initial_mean, initial_variance, &
      variance) result(density)
      !! Evaluate a non-Gaussian local state conditional with regression.
      real(dp), intent(in) :: value !! Input value.
      real(dp), intent(in) :: state(:) !! State vector or state sequence.
      real(dp), intent(in) :: regression(:) !! Regression.
      real(dp), intent(in) :: response(:) !! Response observations.
      real(dp), intent(in) :: scale_data(:) !! Scale data.
      real(dp), intent(in) :: initial_mean !! Initial state mean.
      real(dp), intent(in) :: initial_variance !! Initial variance.
      real(dp), intent(in) :: variance !! Variance value or matrix.
      integer, intent(in) :: time !! Observation times.
      integer, intent(in) :: family !! Family.

      density = non_gaussian_observation_log_density(response(time), &
         scale_data(time), value + regression(time), family)
      if (time == 1) density = density - &
         0.5_dp*(value - initial_mean)**2/initial_variance
      if (time > 1) density = density - &
         0.5_dp*(value - state(time - 1))**2/variance
      if (time < size(state)) density = density - &
         0.5_dp*(state(time + 1) - value)**2/variance
   end function local_level_regression_log_density

   pure real(dp) function non_gaussian_log_likelihood(response, scale_data, &
      predictor, family) result(value)
      !! Sum binomial-logit or Poisson log likelihood contributions.
      real(dp), intent(in) :: response(:) !! Response observations.
      real(dp), intent(in) :: scale_data(:) !! Scale data.
      real(dp), intent(in) :: predictor(:) !! Predictor.
      integer, intent(in) :: family !! Family.
      integer :: time

      value = 0.0_dp
      do time = 1, size(response)
         value = value + non_gaussian_observation_log_density(response(time), &
            scale_data(time), predictor(time), family)
      end do
   end function non_gaussian_log_likelihood

   pure real(dp) function non_gaussian_observation_log_density(response, &
      scale_data, predictor, family) result(value)
      !! Evaluate one non-Gaussian observation log likelihood.
      real(dp), intent(in) :: response !! Response observations.
      real(dp), intent(in) :: scale_data !! Scale data.
      real(dp), intent(in) :: predictor !! Predictor.
      integer, intent(in) :: family !! Family.

      value = 0.0_dp
      if (.not. ieee_is_finite(response)) return
      if (family == 1) then
         value = response*predictor - scale_data* &
            stable_log_one_plus_exp(predictor)
      else
         value = response*predictor - scale_data*exp(min(predictor, 700.0_dp))
      end if
   end function non_gaussian_observation_log_density

   pure function non_gaussian_local_level_draws(response, scale_data, family, &
      initial_mean, initial_variance, variance0, prior_shape, prior_rate, &
      proposal_scale, burn, normals, uniforms, gammas) result(out)
      !! Run componentwise Metropolis updates for a non-Gaussian local level.
      real(dp), intent(in) :: response(:) !! Response observations.
      real(dp), intent(in) :: scale_data(:) !! Scale data.
      integer, intent(in) :: family !! Family.
      integer, intent(in) :: burn !! Number of initial simulation draws to discard.
      real(dp), intent(in) :: initial_mean !! Initial state mean.
      real(dp), intent(in) :: initial_variance !! Initial variance.
      real(dp), intent(in) :: variance0 !! Variance0.
      real(dp), intent(in) :: prior_shape !! Prior shape.
      real(dp), intent(in) :: prior_rate !! Prior rate.
      real(dp), intent(in) :: proposal_scale !! Proposal scale.
      real(dp), intent(in) :: normals(:, :) !! Independent standard-normal draws.
      real(dp), intent(in) :: uniforms(:, :) !! Uniforms.
      real(dp), intent(in) :: gammas(:) !! Gammas.
      type(bsts_non_gaussian_t) :: out
      real(dp), allocatable :: state(:)
      real(dp) :: proposal, old_density, new_density, variance, difference_sum
      integer :: observations, iterations, iteration, time

      observations = size(response)
      iterations = size(normals, 2)
      if (observations < 2 .or. size(scale_data) /= observations .or. &
         iterations < 1 .or. burn < 0 .or. burn >= iterations .or. &
         initial_variance <= 0.0_dp .or. variance0 <= 0.0_dp .or. &
         prior_shape <= 0.0_dp .or. prior_rate <= 0.0_dp .or. &
         proposal_scale <= 0.0_dp .or. &
         any(shape(normals) /= [observations, iterations]) .or. &
         any(shape(uniforms) /= [observations, iterations]) .or. &
         size(gammas) /= iterations .or. any(gammas <= 0.0_dp) .or. &
         any(uniforms <= 0.0_dp) .or. any(uniforms >= 1.0_dp) .or. &
         .not. all(ieee_is_finite(scale_data)) .or. &
         .not. all(ieee_is_finite(normals)) .or. &
         .not. all(ieee_is_finite(uniforms)) .or. &
         .not. all(ieee_is_finite(gammas))) then
         out%info = 1
         return
      end if
      if (family == 1) then
         if (any(scale_data < 0.0_dp) .or. any(response < 0.0_dp) .or. &
            any(response > scale_data)) then
            out%info = 1
            return
         end if
      else if (family == 2) then
         if (any(scale_data <= 0.0_dp) .or. any(response < 0.0_dp)) then
            out%info = 1
            return
         end if
      else
         out%info = 1
         return
      end if
      allocate(state(observations), out%state(observations, iterations), &
         out%fitted_mean(observations, iterations), &
         out%state_variance(iterations))
      state = initial_mean
      variance = variance0
      out%burn = burn
      out%family = family
      do iteration = 1, iterations
         do time = 1, observations
            proposal = state(time) + proposal_scale*normals(time, iteration)
            old_density = local_level_log_density(state(time), time, state, &
               response, scale_data, family, initial_mean, initial_variance, &
               variance)
            new_density = local_level_log_density(proposal, time, state, &
               response, scale_data, family, initial_mean, initial_variance, &
               variance)
            if (log(uniforms(time, iteration)) < new_density - old_density) &
               state(time) = proposal
         end do
         difference_sum = sum((state(2:) - state(:observations - 1))**2)
         variance = (prior_rate + 0.5_dp*difference_sum)/gammas(iteration)
         out%state(:, iteration) = state
         out%state_variance(iteration) = variance
         if (family == 1) then
            out%fitted_mean(:, iteration) = scale_data*logistic_log_odds(state)
         else
            out%fitted_mean(:, iteration) = scale_data*exp(min(state, 700.0_dp))
         end if
      end do
   end function non_gaussian_local_level_draws

   pure real(dp) function local_level_log_density(value, time, state, &
      response, scale_data, family, initial_mean, initial_variance, variance) &
      result(density)
      !! Evaluate one local-level full conditional up to a constant.
      real(dp), intent(in) :: value !! Input value.
      real(dp), intent(in) :: state(:) !! State vector or state sequence.
      real(dp), intent(in) :: response(:) !! Response observations.
      real(dp), intent(in) :: scale_data(:) !! Scale data.
      integer, intent(in) :: time !! Observation times.
      integer, intent(in) :: family !! Family.
      real(dp), intent(in) :: initial_mean !! Initial state mean.
      real(dp), intent(in) :: initial_variance !! Initial variance.
      real(dp), intent(in) :: variance !! Variance value or matrix.

      density = 0.0_dp
      if (ieee_is_finite(response(time))) then
         if (family == 1) then
            density = response(time)*value - scale_data(time)* &
               stable_log_one_plus_exp(value)
         else
            density = response(time)*value - &
               scale_data(time)*exp(min(value, 700.0_dp))
         end if
      end if
      if (time == 1) density = density - &
         0.5_dp*(value - initial_mean)**2/initial_variance
      if (time > 1) density = density - &
         0.5_dp*(value - state(time - 1))**2/variance
      if (time < size(state)) density = density - &
         0.5_dp*(state(time + 1) - value)**2/variance
   end function local_level_log_density

   pure elemental real(dp) function stable_log_one_plus_exp(value) result(ans)
      !! Evaluate log(1 + exp(value)) without overflow.
      real(dp), intent(in) :: value !! Input value.

      if (value > 0.0_dp) then
         ans = value + log(1.0_dp + exp(-value))
      else
         ans = log(1.0_dp + exp(value))
      end if
   end function stable_log_one_plus_exp

   pure function non_gaussian_trend_seasonal_predict_draws(fit, future_scale, &
      normals, uniforms, family) result(out)
      !! Forecast non-Gaussian local trend and seasonal structural states.
      type(bsts_non_gaussian_t), intent(in) :: fit !! Previously fitted model.
      real(dp), intent(in) :: future_scale(:) !! Future scale.
      real(dp), intent(in) :: normals(:, :, :) !! Independent standard-normal draws.
      real(dp), intent(in) :: uniforms(:, :) !! Uniforms.
      integer, intent(in) :: family !! Family.
      type(bsts_prediction_t) :: out
      real(dp), allocatable :: structural(:), next_seasonal(:), sorted(:)
      real(dp) :: predictor, mean_value
      integer :: horizon, retained, dimension, draw, step, source, component

      horizon = size(future_scale)
      retained = size(fit%structural_state, 3) - fit%burn
      dimension = size(fit%structural_state, 1)
      if (fit%family /= family .or. fit%nseasons < 2 .or. horizon < 1 .or. &
         any(shape(normals) /= [dimension, horizon, retained]) .or. &
         any(shape(uniforms) /= [horizon, retained])) then
         out%info = 1
         return
      end if
      allocate(structural(dimension), next_seasonal(dimension - 2), &
         out%draws(horizon, retained), out%mean(horizon), &
         out%standard_deviation(horizon), out%lower(horizon), out%upper(horizon))
      do draw = 1, retained
         source = fit%burn + draw
         structural = fit%structural_state(:, &
            size(fit%structural_state, 2), source)
         do step = 1, horizon
            structural(1) = structural(1) + structural(2) + &
               sqrt(fit%component_variance(1, source))*normals(1, step, draw)
            structural(2) = structural(2) + &
               sqrt(fit%component_variance(2, source))*normals(2, step, draw)
            if (modulo(size(fit%structural_state, 2) + step - 1, &
               fit%season_duration) == 0) then
               next_seasonal = 0.0_dp
               next_seasonal(1) = -sum(structural(3:)) + &
                  sqrt(fit%component_variance(3, source))*normals(3, step, draw)
               do component = 2, size(next_seasonal)
                  next_seasonal(component) = structural(component + 1)
               end do
               structural(3:) = next_seasonal
            end if
            predictor = structural(1) + structural(3)
            if (family == 1) then
               out%draws(step, draw) = real(binomial_quantile( &
                  nint(future_scale(step)), logistic_log_odds(predictor), &
                  uniforms(step, draw)), dp)
            else
               mean_value = future_scale(step)*exp(min(predictor, 700.0_dp))
               out%draws(step, draw) = real(poisson_quantile(mean_value, &
                  uniforms(step, draw)), dp)
            end if
         end do
      end do
      do step = 1, horizon
         out%mean(step) = sum(out%draws(step, :))/real(retained, dp)
         if (retained > 1) then
            out%standard_deviation(step) = sqrt(sum((out%draws(step, :) - &
               out%mean(step))**2)/real(retained - 1, dp))
         else
            out%standard_deviation(step) = 0.0_dp
         end if
         sorted = out%draws(step, :)
         call insertion_sort(sorted)
         out%lower(step) = quantile(sorted, 0.025_dp)
         out%upper(step) = quantile(sorted, 0.975_dp)
      end do
   end function non_gaussian_trend_seasonal_predict_draws

   pure function non_gaussian_regression_predict_draws(fit, future_scale, &
      future_predictors, normals, uniforms, family) result(out)
      !! Forecast sparse non-Gaussian regression from supplied draws.
      type(bsts_non_gaussian_t), intent(in) :: fit !! Previously fitted model.
      real(dp), intent(in) :: future_scale(:) !! Future scale.
      real(dp), intent(in) :: future_predictors(:, :) !! Future predictors.
      real(dp), intent(in) :: normals(:, :) !! Independent standard-normal draws.
      real(dp), intent(in) :: uniforms(:, :) !! Uniforms.
      integer, intent(in) :: family !! Family.
      type(bsts_prediction_t) :: out
      real(dp), allocatable :: sorted(:)
      real(dp) :: state, predictor, mean_value
      integer :: horizon, retained, variables, draw, step, source

      horizon = size(future_scale)
      retained = size(fit%state, 2) - fit%burn
      variables = size(fit%coefficients, 1)
      if (fit%family /= family .or. horizon < 1 .or. retained < 1 .or. &
         any(shape(future_predictors) /= [horizon, variables]) .or. &
         any(shape(normals) /= [horizon, retained]) .or. &
         any(shape(uniforms) /= [horizon, retained]) .or. &
         any(future_scale < 0.0_dp) .or. any(uniforms <= 0.0_dp) .or. &
         any(uniforms >= 1.0_dp)) then
         out%info = 1
         return
      end if
      allocate(out%draws(horizon, retained), out%mean(horizon), &
         out%standard_deviation(horizon), out%lower(horizon), out%upper(horizon))
      do draw = 1, retained
         source = fit%burn + draw
         state = fit%state(size(fit%state, 1), source)
         do step = 1, horizon
            state = state + sqrt(fit%state_variance(source))*normals(step, draw)
            predictor = state + dot_product(future_predictors(step, :), &
               fit%coefficients(:, source))
            if (family == 1) then
               out%draws(step, draw) = real(binomial_quantile( &
                  nint(future_scale(step)), logistic_log_odds(predictor), &
                  uniforms(step, draw)), dp)
            else
               mean_value = future_scale(step)*exp(min(predictor, 700.0_dp))
               out%draws(step, draw) = real(poisson_quantile(mean_value, &
                  uniforms(step, draw)), dp)
            end if
         end do
      end do
      do step = 1, horizon
         out%mean(step) = sum(out%draws(step, :))/real(retained, dp)
         if (retained > 1) then
            out%standard_deviation(step) = sqrt(sum((out%draws(step, :) - &
               out%mean(step))**2)/real(retained - 1, dp))
         else
            out%standard_deviation(step) = 0.0_dp
         end if
         sorted = out%draws(step, :)
         call insertion_sort(sorted)
         out%lower(step) = quantile(sorted, 0.025_dp)
         out%upper(step) = quantile(sorted, 0.975_dp)
      end do
   end function non_gaussian_regression_predict_draws

   pure function non_gaussian_predict_draws(fit, future_scale, normals, &
      uniforms, family) result(out)
      !! Forecast a non-Gaussian local-level model from supplied draws.
      type(bsts_non_gaussian_t), intent(in) :: fit !! Previously fitted model.
      real(dp), intent(in) :: future_scale(:) !! Future scale.
      real(dp), intent(in) :: normals(:, :) !! Independent standard-normal draws.
      real(dp), intent(in) :: uniforms(:, :) !! Uniforms.
      integer, intent(in) :: family !! Family.
      type(bsts_prediction_t) :: out
      real(dp), allocatable :: sorted(:)
      real(dp) :: state, mean_value
      integer :: horizon, retained, draw, step, source

      horizon = size(future_scale)
      retained = 0
      if (allocated(fit%state)) retained = size(fit%state, 2) - fit%burn
      if (fit%info /= 0 .or. allocated(fit%coefficients) .or. &
         fit%family /= family .or. horizon < 1 .or. &
         retained < 1 .or. any(shape(normals) /= [horizon, retained]) .or. &
         any(shape(uniforms) /= [horizon, retained]) .or. &
         any(future_scale < 0.0_dp) .or. any(uniforms <= 0.0_dp) .or. &
         any(uniforms >= 1.0_dp)) then
         out%info = 1
         return
      end if
      if (family == 2 .and. any(future_scale <= 0.0_dp)) then
         out%info = 1
         return
      end if
      allocate(out%draws(horizon, retained), out%mean(horizon), &
         out%standard_deviation(horizon), out%lower(horizon), &
         out%upper(horizon))
      do draw = 1, retained
         source = fit%burn + draw
         state = fit%state(size(fit%state, 1), source)
         do step = 1, horizon
            state = state + sqrt(fit%state_variance(source))*normals(step, draw)
            if (family == 1) then
               mean_value = future_scale(step)*logistic_log_odds(state)
               out%draws(step, draw) = real(binomial_quantile( &
                  nint(future_scale(step)), logistic_log_odds(state), &
                  uniforms(step, draw)), dp)
            else
               mean_value = future_scale(step)*exp(min(state, 700.0_dp))
               out%draws(step, draw) = real(poisson_quantile(mean_value, &
                  uniforms(step, draw)), dp)
            end if
         end do
      end do
      do step = 1, horizon
         out%mean(step) = sum(out%draws(step, :))/real(retained, dp)
         if (retained > 1) then
            out%standard_deviation(step) = sqrt(sum((out%draws(step, :) - &
               out%mean(step))**2)/real(retained - 1, dp))
         else
            out%standard_deviation(step) = 0.0_dp
         end if
         sorted = out%draws(step, :)
         call insertion_sort(sorted)
         out%lower(step) = quantile(sorted, 0.025_dp)
         out%upper(step) = quantile(sorted, 0.975_dp)
      end do
   end function non_gaussian_predict_draws

   pure integer function binomial_quantile(trials, probability, uniform) &
      result(value)
      !! Transform a uniform draw into a binomial count.
      integer, intent(in) :: trials !! Trials.
      real(dp), intent(in) :: probability !! Probability value.
      real(dp), intent(in) :: uniform !! Uniform.
      real(dp) :: mass, cumulative

      if (trials <= 0 .or. probability <= 0.0_dp) then
         value = 0
         return
      end if
      if (probability >= 1.0_dp) then
         value = trials
         return
      end if
      mass = (1.0_dp - probability)**trials
      cumulative = mass
      value = 0
      do while (uniform > cumulative .and. value < trials)
         mass = mass*real(trials - value, dp)/real(value + 1, dp)* &
            probability/(1.0_dp - probability)
         value = value + 1
         cumulative = cumulative + mass
      end do
   end function binomial_quantile

   pure integer function poisson_quantile(mean_value, uniform) result(value)
      !! Transform a uniform draw into a Poisson count.
      real(dp), intent(in) :: mean_value !! Mean value.
      real(dp), intent(in) :: uniform !! Uniform.
      real(dp) :: mass, cumulative

      if (mean_value <= 0.0_dp) then
         value = 0
         return
      end if
      if (mean_value > 700.0_dp) then
         value = max(0, nint(mean_value))
         return
      end if
      mass = exp(-mean_value)
      cumulative = mass
      value = 0
      do while (uniform > cumulative .and. value < 1000000)
         value = value + 1
         mass = mass*mean_value/real(value, dp)
         cumulative = cumulative + mass
      end do
   end function poisson_quantile

   subroutine fill_non_gaussian_forecast_draws(fit, horizon, normals, uniforms)
      !! Generate Gaussian state and uniform observation forecast draws.
      type(bsts_non_gaussian_t), intent(in) :: fit !! Previously fitted model.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      real(dp), allocatable, intent(out) :: normals(:, :) !! Independent standard-normal draws.
      real(dp), allocatable, intent(out) :: uniforms(:, :) !! Uniforms.
      integer :: retained, draw, step

      if (.not. allocated(fit%state) .or. horizon < 1) return
      retained = size(fit%state, 2) - fit%burn
      if (retained < 1) return
      allocate(normals(horizon, retained), uniforms(horizon, retained))
      do draw = 1, retained
         do step = 1, horizon
            normals(step, draw) = random_standard_normal()
            uniforms(step, draw) = max(tiny(1.0_dp), &
               min(1.0_dp - epsilon(1.0_dp), random_uniform()))
         end do
      end do
   end subroutine fill_non_gaussian_forecast_draws

   pure function bsts_geometric_sequence(length, initial_value, &
      discount_factor) result(sequence)
      !! Return bsts GeometricSequence values without graphics or R objects.
      integer, intent(in) :: length !! Length.
      real(dp), intent(in), optional :: initial_value !! Initial value.
      real(dp), intent(in), optional :: discount_factor !! Discount factor.
      real(dp), allocatable :: sequence(:)
      real(dp) :: initial, discount
      integer :: index

      initial = 1.0_dp
      if (present(initial_value)) initial = initial_value
      discount = 0.5_dp
      if (present(discount_factor)) discount = discount_factor
      if (length < 1 .or. initial == 0.0_dp .or. discount == 0.0_dp) then
         allocate(sequence(0))
         return
      end if
      allocate(sequence(length))
      do index = 1, length
         sequence(index) = initial*discount**(index - 1)
      end do
   end function bsts_geometric_sequence

   pure function bsts_harvey_cumulator(fine_series, contains_end, &
      membership_fraction) result(aggregate)
      !! Form Harvey's running partial aggregates across coarse boundaries.
      real(dp), intent(in) :: fine_series(:) !! Fine series.
      real(dp), intent(in) :: membership_fraction(:) !! Membership fraction.
      logical, intent(in) :: contains_end(:) !! Flag controlling contains end.
      real(dp), allocatable :: aggregate(:)
      real(dp) :: cumulator, fraction
      integer :: time

      if (size(fine_series) < 1 .or. &
         size(contains_end) /= size(fine_series) .or. &
         (size(membership_fraction) /= 1 .and. &
         size(membership_fraction) /= size(fine_series))) then
         allocate(aggregate(0))
         return
      end if
      if (size(fine_series) == 1) then
         aggregate = fine_series
         return
      end if
      allocate(aggregate(size(fine_series)))
      cumulator = 0.0_dp
      do time = 1, size(fine_series)
         fraction = membership_fraction(min(time, size(membership_fraction)))
         if (fraction < 0.0_dp .or. fraction > 1.0_dp) then
            deallocate(aggregate)
            allocate(aggregate(0))
            return
         end if
         if (contains_end(time)) then
            cumulator = cumulator + fraction*fine_series(time)
            aggregate(time) = cumulator
            cumulator = (1.0_dp - fraction)*fine_series(time)
         else
            cumulator = cumulator + fine_series(time)
            aggregate(time) = cumulator
         end if
      end do
   end function bsts_harvey_cumulator

   pure function mixed_frequency_design(coarse_index, membership_fraction, &
      contains_end, coarse_points) result(design)
      !! Build the Harvey flow-aggregation matrix for fine observations.
      integer, intent(in) :: coarse_index(:) !! Index of coarse.
      integer, intent(in) :: coarse_points !! Coarse points.
      real(dp), intent(in) :: membership_fraction(:) !! Membership fraction.
      logical, intent(in) :: contains_end(:) !! Flag controlling contains end.
      real(dp), allocatable :: design(:, :)
      real(dp), allocatable :: accumulated(:)
      integer :: fine_points, time, coarse

      fine_points = size(coarse_index)
      allocate(design(coarse_points, fine_points), accumulated(fine_points))
      design = 0.0_dp
      accumulated = 0.0_dp
      do time = 1, fine_points
         if (contains_end(time)) then
            coarse = coarse_index(time)
            if (coarse >= 1 .and. coarse <= coarse_points) then
               design(coarse, :) = accumulated
               design(coarse, time) = design(coarse, time) + &
                  membership_fraction(time)
            end if
            accumulated = 0.0_dp
            accumulated(time) = 1.0_dp - membership_fraction(time)
         else
            accumulated(time) = accumulated(time) + 1.0_dp
         end if
      end do
   end function mixed_frequency_design

   pure function mixed_frequency_cumulator(fine_series, membership_fraction, &
      contains_end) result(cumulator)
      !! Return the accumulated flow immediately before each fine observation.
      real(dp), intent(in) :: fine_series(:) !! Fine series.
      real(dp), intent(in) :: membership_fraction(:) !! Membership fraction.
      logical, intent(in) :: contains_end(:) !! Flag controlling contains end.
      real(dp), allocatable :: cumulator(:)
      real(dp) :: accumulated
      integer :: time

      allocate(cumulator(size(fine_series)))
      accumulated = 0.0_dp
      do time = 1, size(fine_series)
         cumulator(time) = accumulated
         if (contains_end(time)) then
            accumulated = (1.0_dp - membership_fraction(time))*fine_series(time)
         else
            accumulated = accumulated + fine_series(time)
         end if
      end do
   end function mixed_frequency_cumulator

   pure subroutine summarize_prediction(prediction)
      !! Fill pointwise posterior summaries for an allocated draw matrix.
      type(bsts_prediction_t), intent(inout) :: prediction !! Prediction, updated in place.
      real(dp), allocatable :: sorted(:)
      integer :: points, draws, point

      points = size(prediction%draws, 1)
      draws = size(prediction%draws, 2)
      allocate(prediction%mean(points), prediction%standard_deviation(points), &
         prediction%lower(points), prediction%upper(points))
      do point = 1, points
         prediction%mean(point) = sum(prediction%draws(point, :))/real(draws, dp)
         if (draws > 1) then
            prediction%standard_deviation(point) = sqrt(sum( &
               (prediction%draws(point, :) - prediction%mean(point))**2)/ &
               real(draws - 1, dp))
         else
            prediction%standard_deviation(point) = 0.0_dp
         end if
         sorted = prediction%draws(point, :)
         call insertion_sort(sorted)
         prediction%lower(point) = quantile(sorted, 0.025_dp)
         prediction%upper(point) = quantile(sorted, 0.975_dp)
      end do
   end subroutine summarize_prediction

   pure subroutine prediction_error_summaries(errors)
      !! Summarize posterior mean prediction errors and overall loss.
      type(bsts_prediction_errors_t), intent(inout) :: errors !! Errors, updated in place.
      logical, allocatable :: finite(:)
      integer :: time, count_finite

      do time = 1, size(errors%draws, 1)
         finite = ieee_is_finite(errors%draws(time, :))
         count_finite = count(finite)
         if (count_finite > 0) then
            errors%mean(time) = sum(errors%draws(time, :), mask=finite)/ &
               real(count_finite, dp)
         else
            errors%mean(time) = ieee_value(0.0_dp, ieee_quiet_nan)
         end if
      end do
      finite = ieee_is_finite(errors%mean)
      count_finite = count(finite)
      if (count_finite > 0) then
         errors%rmse = sqrt(sum(errors%mean**2, mask=finite)/ &
            real(count_finite, dp))
         errors%mae = sum(abs(errors%mean), mask=finite)/real(count_finite, dp)
      end if
   end subroutine prediction_error_summaries

   pure subroutine score_ranks(scores, ranks)
      !! Assign one-based ascending ranks, preserving tied ranks.
      real(dp), intent(in) :: scores(:) !! Scores.
      integer, intent(out) :: ranks(:) !! Ranks.
      integer :: model

      do model = 1, size(scores)
         ranks(model) = 1 + count(scores < scores(model))
      end do
   end subroutine score_ranks

   pure function sorted_unique_real(values) result(unique)
      !! Return sorted distinct real values.
      real(dp), intent(in) :: values(:) !! Input values.
      real(dp), allocatable :: sorted(:), unique(:)
      integer :: index, count_unique

      sorted = values
      call insertion_sort(sorted)
      if (size(sorted) == 0) then
         allocate(unique(0))
         return
      end if
      count_unique = 1
      do index = 2, size(sorted)
         if (sorted(index) /= sorted(index - 1)) count_unique = count_unique + 1
      end do
      allocate(unique(count_unique))
      unique(1) = sorted(1)
      count_unique = 1
      do index = 2, size(sorted)
         if (sorted(index) == sorted(index - 1)) cycle
         count_unique = count_unique + 1
         unique(count_unique) = sorted(index)
      end do
   end function sorted_unique_real

   pure function sorted_unique_integer(values) result(unique)
      !! Return sorted distinct integer values.
      integer, intent(in) :: values(:) !! Input values.
      integer, allocatable :: sorted(:), unique(:)
      integer :: index, position, candidate, count_unique

      sorted = values
      do index = 2, size(sorted)
         candidate = sorted(index)
         position = index - 1
         do while (position >= 1)
            if (sorted(position) <= candidate) exit
            sorted(position + 1) = sorted(position)
            position = position - 1
         end do
         sorted(position + 1) = candidate
      end do
      if (size(sorted) == 0) then
         allocate(unique(0))
         return
      end if
      count_unique = 1
      do index = 2, size(sorted)
         if (sorted(index) /= sorted(index - 1)) count_unique = count_unique + 1
      end do
      allocate(unique(count_unique))
      unique(1) = sorted(1)
      count_unique = 1
      do index = 2, size(sorted)
         if (sorted(index) == sorted(index - 1)) cycle
         count_unique = count_unique + 1
         unique(count_unique) = sorted(index)
      end do
   end function sorted_unique_integer

   pure integer function integer_position(values, target) result(position)
      !! Return the one-based position of an integer in a short vector.
      integer, intent(in) :: values(:) !! Input values.
      integer, intent(in) :: target !! Target.
      integer :: index

      position = 0
      do index = 1, size(values)
         if (values(index) == target) then
            position = index
            return
         end if
      end do
   end function integer_position

   pure elemental type(date_t) function add_calendar_months(value, months) &
      result(shifted)
      !! Add calendar months while clamping the day within the target month.
      type(date_t), intent(in) :: value !! Input value.
      integer, intent(in) :: months !! Months.
      integer :: total, year, month, day

      total = 12*value%year + value%month - 1 + months
      year = floor(real(total, dp)/12.0_dp)
      month = total - 12*year + 1
      day = min(value%day, date_days_in_month(year, month))
      shifted = date_t(year, month, day)
   end function add_calendar_months

   pure elemental type(date_t) function last_day_in_month(value) result(last)
      !! Return the final Gregorian day in a date's month.
      type(date_t), intent(in) :: value !! Input value.

      last = date_t(value%year, value%month, &
         date_days_in_month(value%year, value%month))
   end function last_day_in_month

   subroutine fill_standard_normals(draws)
      !! Fill an array with independent standard-normal draws.
      real(dp), intent(out) :: draws(:, :) !! Draws.
      integer :: row, column

      do column = 1, size(draws, 2)
         do row = 1, size(draws, 1)
            draws(row, column) = random_standard_normal()
         end do
      end do
   end subroutine fill_standard_normals

   pure subroutine impose_shared_loading_constraints(loadings)
      !! Impose the lower-triangular unit-diagonal factor identification.
      real(dp), intent(inout) :: loadings(:, :) !! Loadings, updated in place.
      integer :: series, factor

      do series = 1, size(loadings, 1)
         do factor = series + 1, size(loadings, 2)
            loadings(series, factor) = 0.0_dp
         end do
         if (series <= size(loadings, 2)) loadings(series, series) = 1.0_dp
      end do
   end subroutine impose_shared_loading_constraints

   pure subroutine impose_shared_inclusion_priors(probability)
      !! Force identifying loading probabilities to zero or one.
      real(dp), intent(inout) :: probability(:, :) !! Probability value, updated in place.
      integer :: series, factor

      do series = 1, size(probability, 1)
         do factor = series + 1, size(probability, 2)
            probability(series, factor) = 0.0_dp
         end do
         if (series <= size(probability, 2)) &
            probability(series, series) = 1.0_dp
      end do
   end subroutine impose_shared_inclusion_priors

   pure subroutine impose_shared_inclusion_constraints(included)
      !! Force identifying loading indicators to their fixed states.
      logical, intent(inout) :: included(:, :) !! Flag controlling included, updated in place.
      integer :: series, factor

      do series = 1, size(included, 1)
         do factor = series + 1, size(included, 2)
            included(series, factor) = .false.
         end do
         if (series <= size(included, 2)) included(series, series) = .true.
      end do
   end subroutine impose_shared_inclusion_constraints

   pure subroutine inverse_wishart_from_draws(degrees_of_freedom, scale, &
      normal_draws, gamma_draws, value, info)
      !! Transform Bartlett draws into an inverse-Wishart matrix.
      real(dp), intent(in) :: degrees_of_freedom !! Degrees of freedom.
      real(dp), intent(in) :: scale(:, :) !! Scale.
      real(dp), intent(in) :: normal_draws(:, :) !! Independent standard-normal draws.
      real(dp), intent(in) :: gamma_draws(:) !! Gamma simulation draws.
      real(dp), allocatable, intent(out) :: value(:, :) !! Input value.
      integer, intent(out) :: info !! Status code; zero indicates success.
      real(dp), allocatable :: inverse_scale(:, :), scale_factor(:, :)
      real(dp), allocatable :: bartlett(:, :), factor(:, :), precision(:, :)
      integer :: dimension, row, column, status

      dimension = size(scale, 1)
      if (dimension < 1 .or. size(scale, 2) /= dimension .or. &
         degrees_of_freedom <= real(dimension - 1, dp) .or. &
         any(shape(normal_draws) /= [dimension, dimension]) .or. &
         size(gamma_draws) /= dimension .or. any(gamma_draws <= 0.0_dp) .or. &
         .not. all(ieee_is_finite(normal_draws))) then
         allocate(value(0, 0))
         info = 1
         return
      end if
      call invert_matrix(scale, inverse_scale, status)
      if (status /= 0) then
         allocate(value(0, 0))
         info = 2
         return
      end if
      call cholesky_lower(inverse_scale, scale_factor, status)
      if (status /= 0) then
         allocate(value(0, 0))
         info = 3
         return
      end if
      allocate(bartlett(dimension, dimension))
      bartlett = 0.0_dp
      do row = 1, dimension
         bartlett(row, row) = sqrt(2.0_dp*gamma_draws(row))
         do column = 1, row - 1
            bartlett(row, column) = normal_draws(row, column)
         end do
      end do
      factor = matmul(scale_factor, bartlett)
      precision = matmul(factor, transpose(factor))
      call invert_matrix(precision, value, status)
      if (status /= 0) then
         info = 4
         return
      end if
      value = symmetrize(value)
      info = 0
   end subroutine inverse_wishart_from_draws

   pure real(dp) function finite_sample_variance(values) result(variance)
      !! Compute the sample variance after omitting nonfinite observations.
      real(dp), intent(in) :: values(:) !! Input values.
      real(dp) :: mean
      integer :: count_finite

      count_finite = count(ieee_is_finite(values))
      if (count_finite < 2) then
         variance = 1.0_dp
         return
      end if
      mean = sum(values, mask=ieee_is_finite(values))/real(count_finite, dp)
      variance = sum((values - mean)**2, mask=ieee_is_finite(values))/ &
         real(count_finite - 1, dp)
      variance = max(variance, tiny(1.0_dp))
   end function finite_sample_variance

   pure real(dp) function first_finite_value(values) result(value)
      !! Return the first finite observation or zero when none exists.
      real(dp), intent(in) :: values(:) !! Input values.
      integer :: index

      value = 0.0_dp
      do index = 1, size(values)
         if (ieee_is_finite(values(index))) then
            value = values(index)
            return
         end if
      end do
   end function first_finite_value

   pure real(dp) function finite_endpoint_slope(values) result(slope)
      !! Estimate the endpoint slope using the first and last finite values.
      real(dp), intent(in) :: values(:) !! Input values.
      integer :: first, last

      first = 1
      do while (first <= size(values))
         if (ieee_is_finite(values(first))) exit
         first = first + 1
      end do
      last = size(values)
      do while (last >= first)
         if (ieee_is_finite(values(last))) exit
         last = last - 1
      end do
      if (last > first) then
         slope = (values(last) - values(first))/real(last - first + 1, dp)
      else
         slope = 0.0_dp
      end if
   end function finite_endpoint_slope

   pure subroutine insertion_sort(values)
      !! Sort a short real vector in ascending order.
      real(dp), intent(inout) :: values(:) !! Input values, updated in place.
      real(dp) :: candidate
      integer :: i, j

      do i = 2, size(values)
         candidate = values(i)
         j = i - 1
         do while (j >= 1)
            if (values(j) <= candidate) exit
            values(j + 1) = values(j)
            j = j - 1
         end do
         values(j + 1) = candidate
      end do
   end subroutine insertion_sort

   pure function pack_predictors(response, predictors) result(observed)
      !! Retain predictor rows whose corresponding response is finite.
      real(dp), intent(in) :: response(:) !! Response observations.
      real(dp), intent(in) :: predictors(:, :) !! Predictor matrix.
      real(dp), allocatable :: observed(:, :)
      integer :: row, destination

      allocate(observed(count(ieee_is_finite(response)), size(predictors, 2)))
      destination = 0
      do row = 1, size(response)
         if (.not. ieee_is_finite(response(row))) cycle
         destination = destination + 1
         observed(destination, :) = predictors(row, :)
      end do
   end function pack_predictors

   pure subroutine spike_slab_model_moments(response, predictors, included, &
      slab_mean, slab_covariance, inclusion_probability, variance, prior_rate, &
      log_probability, collapsed_rate, posterior_mean, posterior_covariance, &
      active, info)
      !! Compute collapsed likelihood and active-coefficient posterior moments.
      real(dp), intent(in) :: response(:) !! Response observations.
      real(dp), intent(in) :: predictors(:, :) !! Predictor matrix.
      real(dp), intent(in) :: slab_mean(:) !! Slab mean.
      real(dp), intent(in) :: slab_covariance(:, :) !! Slab covariance matrix.
      real(dp), intent(in) :: inclusion_probability(:) !! Inclusion probability.
      real(dp), intent(in) :: variance !! Variance value or matrix.
      real(dp), intent(in) :: prior_rate !! Prior rate.
      logical, intent(in) :: included(:) !! Flag controlling included.
      real(dp), intent(out) :: log_probability !! Log probability.
      real(dp), intent(out) :: collapsed_rate !! Collapsed rate.
      real(dp), allocatable, intent(out) :: posterior_mean(:) !! Posterior mean.
      real(dp), allocatable, intent(out) :: posterior_covariance(:, :) !! Posterior covariance matrix.
      integer, allocatable, intent(out) :: active(:) !! Active.
      integer, intent(out) :: info !! Status code; zero indicates success.
      real(dp), allocatable :: design(:, :), covariance(:, :), precision(:, :)
      real(dp), allocatable :: prior_precision(:, :), prior_mean(:), score(:)
      real(dp) :: log_prior, logdet_prior, logdet_precision, sum_squares
      integer :: variable, status

      active = pack([(variable, variable=1, size(included))], included)
      log_prior = 0.0_dp
      do variable = 1, size(included)
         if (included(variable)) then
            if (inclusion_probability(variable) <= 0.0_dp) then
               log_probability = -huge(1.0_dp)
               collapsed_rate = huge(1.0_dp)
               allocate(posterior_mean(0), posterior_covariance(0, 0))
               info = 0
               return
            end if
            log_prior = log_prior + log(inclusion_probability(variable))
         else
            if (inclusion_probability(variable) >= 1.0_dp) then
               log_probability = -huge(1.0_dp)
               collapsed_rate = huge(1.0_dp)
               allocate(posterior_mean(0), posterior_covariance(0, 0))
               info = 0
               return
            end if
            log_prior = log_prior + log(1.0_dp - inclusion_probability(variable))
         end if
      end do
      if (size(active) == 0) then
         allocate(posterior_mean(0), posterior_covariance(0, 0))
         sum_squares = dot_product(response, response)
         collapsed_rate = prior_rate + 0.5_dp*sum_squares
         log_probability = log_prior - 0.5_dp*real(size(response), dp)* &
            log(variance) - 0.5_dp*sum_squares/variance
         info = 0
         return
      end if
      design = predictors(:, active)
      covariance = slab_covariance(active, active)
      prior_mean = slab_mean(active)
      allocate(prior_precision(size(active), size(active)), &
         posterior_covariance(size(active), size(active)))
      call inverse_logdet(covariance, prior_precision, logdet_prior, status, &
         100.0_dp*epsilon(1.0_dp))
      if (status /= 0) then
         info = 1
         return
      end if
      precision = matmul(transpose(design), design) + prior_precision
      call inverse_logdet(precision, posterior_covariance, logdet_precision, &
         status, 100.0_dp*epsilon(1.0_dp))
      if (status /= 0) then
         info = 2
         return
      end if
      score = matmul(transpose(design), response) + &
         matmul(prior_precision, prior_mean)
      posterior_mean = matmul(posterior_covariance, score)
      sum_squares = dot_product(response, response) + &
         dot_product(prior_mean, matmul(prior_precision, prior_mean)) - &
         dot_product(score, posterior_mean)
      sum_squares = max(0.0_dp, sum_squares)
      collapsed_rate = prior_rate + 0.5_dp*sum_squares
      log_probability = log_prior - 0.5_dp*(logdet_prior + &
         logdet_precision) - 0.5_dp*real(size(response), dp)*log(variance) - &
         0.5_dp*sum_squares/variance
      info = 0
   end subroutine spike_slab_model_moments

   pure elemental real(dp) function logistic_log_odds(log_odds) result(value)
      !! Convert log odds to probability without avoidable overflow.
      real(dp), intent(in) :: log_odds !! Log odds.

      if (log_odds >= 0.0_dp) then
         value = 1.0_dp/(1.0_dp + exp(-min(log_odds, 700.0_dp)))
      else
         value = exp(max(log_odds, -700.0_dp))/ &
            (1.0_dp + exp(max(log_odds, -700.0_dp)))
      end if
   end function logistic_log_odds

   pure real(dp) function truncated_normal_draw(mean, standard_deviation, &
      lower, upper, uniform) result(value)
      !! Transform a uniform variate into a bounded Gaussian draw.
      real(dp), intent(in) :: mean !! Mean value or vector.
      real(dp), intent(in) :: standard_deviation !! Standard deviation.
      real(dp), intent(in) :: lower !! Lower.
      real(dp), intent(in) :: upper !! Upper.
      real(dp), intent(in) :: uniform !! Uniform.
      real(dp) :: lower_probability, upper_probability, probability

      lower_probability = 0.0_dp
      if (lower > -huge(1.0_dp)/2.0_dp) lower_probability = &
         0.5_dp*erfc(-(lower - mean)/(sqrt(2.0_dp)*standard_deviation))
      upper_probability = 1.0_dp
      if (upper < huge(1.0_dp)/2.0_dp) upper_probability = &
         0.5_dp*erfc(-(upper - mean)/(sqrt(2.0_dp)*standard_deviation))
      probability = lower_probability + uniform* &
         max(upper_probability - lower_probability, epsilon(1.0_dp))
      probability = max(epsilon(1.0_dp), min(1.0_dp - epsilon(1.0_dp), &
         probability))
      value = mean + standard_deviation*normal_quantile(probability)
      value = max(lower, min(upper, value))
   end function truncated_normal_draw

end module bsts_mod
