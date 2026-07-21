! SPDX-License-Identifier: MIT
! SPDX-FileComment: Automatic univariate time-series profiling and model selection.
module automatic_modeling_mod
   !! Profile a univariate series, fit suitable candidates, and rank forecasts.
   use kind_mod, only: dp
   use forecast_mod, only: forecast_result_t, meanf, naive, rwf, snaive, &
      ses, holt, holt_winters, thetaf
   use time_series_stats_mod, only: yule_walker_result_t, yule_walker_fit, &
      acf_values, pacf_values
   use astsa_mod, only: astsa_sarima_fit_t, astsa_sarima_forecast_t, &
      astsa_sarima_likelihood_t, sarima_fit, sarima_forecast, &
      sarima_likelihood
   use stats_mod, only: ols_fit
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, &
      ieee_quiet_nan
   implicit none
   private

   integer, parameter :: model_mean = 1
   integer, parameter :: model_naive = 2
   integer, parameter :: model_drift = 3
   integer, parameter :: model_ses = 4
   integer, parameter :: model_holt = 5
   integer, parameter :: model_theta = 6
   integer, parameter :: model_seasonal_naive = 7
   integer, parameter :: model_holt_winters = 8
   integer, parameter :: model_ar = 9
   integer, parameter :: model_arma = 10
   character(len=*), parameter :: selection_validation = "validation"
   character(len=*), parameter :: selection_aicc = "aicc"
   character(len=*), parameter :: selection_bic = "bic"

   type, public :: automatic_model_options_t
      !! Controls profiling, validation, and forecast dimensions.
      integer :: frequency = 1
      integer :: horizon = 10
      integer :: max_lag = 24
      integer :: max_ar_order = 12
      integer :: validation_size = 0
      logical :: include_seasonal = .true.
      logical :: time_fits = .false.
      character(len=10) :: selection = selection_validation
   end type automatic_model_options_t

   type, public :: series_profile_t
      !! Statistical properties used to gate automatic-model candidates.
      integer :: observations = 0
      integer :: frequency = 1
      integer :: max_lag = 0
      real(dp) :: mean = 0.0_dp
      real(dp) :: variance = 0.0_dp
      real(dp) :: minimum = 0.0_dp
      real(dp) :: maximum = 0.0_dp
      real(dp) :: trend_strength = 0.0_dp
      real(dp) :: seasonal_strength = 0.0_dp
      real(dp) :: dependence_threshold = 0.0_dp
      logical :: trend_detected = .false.
      logical :: seasonality_detected = .false.
      logical :: autocorrelation_detected = .false.
      logical :: conditional_variance_detected = .false.
      logical :: differencing_suggested = .false.
      real(dp), allocatable :: acf(:)
      real(dp), allocatable :: pacf(:)
      real(dp), allocatable :: squared_acf(:)
      integer :: info = 0
   end type series_profile_t

   type, public :: order_search_result_t
      !! Fit diagnostics for one AR or ARMA order considered during selection.
      integer :: ar_order = 0
      integer :: ma_order = 0
      real(dp) :: mean = 0.0_dp
      real(dp) :: innovation_variance = 0.0_dp
      real(dp) :: search_score = huge(1.0_dp)
      real(dp) :: log_likelihood = -huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      real(dp) :: aicc = huge(1.0_dp)
      real(dp) :: bic = huge(1.0_dp)
      real(dp), allocatable :: ar_coefficients(:)
      real(dp), allocatable :: ma_coefficients(:)
      logical :: converged = .false.
      integer :: info = 0
   end type order_search_result_t

   type, public :: candidate_model_result_t
      !! Validation performance and fitted settings for one candidate family.
      character(len=:), allocatable :: name
      character(len=:), allocatable :: rationale
      integer :: model_code = 0
      integer :: order = 0
      integer :: ma_order = 0
      real(dp) :: alpha = 0.0_dp
      real(dp) :: beta = 0.0_dp
      real(dp) :: gamma = 0.0_dp
      real(dp) :: level = 0.0_dp
      real(dp) :: drift = 0.0_dp
      real(dp) :: innovation_variance = 0.0_dp
      real(dp), allocatable :: coefficients(:)
      real(dp), allocatable :: ma_coefficients(:)
      type(order_search_result_t), allocatable :: order_fits(:)
      real(dp) :: rmse = huge(1.0_dp)
      real(dp) :: mae = huge(1.0_dp)
      real(dp) :: log_likelihood = -huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      real(dp) :: aicc = huge(1.0_dp)
      real(dp) :: bic = huge(1.0_dp)
      integer :: parameter_count = 0
      integer :: effective_observations = 0
      real(dp) :: full_fit_seconds = 0.0_dp
      real(dp) :: validation_fit_seconds = 0.0_dp
      logical :: converged = .false.
      logical :: selected = .false.
      integer :: info = 0
   end type candidate_model_result_t

   type, public :: automatic_model_result_t
      !! Profile, ranked candidates, selected fit, and future forecasts.
      type(series_profile_t) :: profile
      type(candidate_model_result_t), allocatable :: candidates(:)
      real(dp), allocatable :: forecast(:)
      real(dp), allocatable :: fitted(:)
      real(dp), allocatable :: residuals(:)
      real(dp), allocatable :: validation_actual(:)
      integer :: validation_size = 0
      character(len=10) :: selection_criterion = selection_validation
      logical :: fit_times_recorded = .false.
      integer :: selected = 0
      integer :: info = 0
   end type automatic_model_result_t

   type :: candidate_fit_t
      real(dp), allocatable :: forecast(:)
      real(dp), allocatable :: fitted(:)
      real(dp), allocatable :: residuals(:)
      integer :: order = 0
      integer :: ma_order = 0
      real(dp) :: alpha = 0.0_dp
      real(dp) :: beta = 0.0_dp
      real(dp) :: gamma = 0.0_dp
      real(dp) :: level = 0.0_dp
      real(dp) :: drift = 0.0_dp
      real(dp), allocatable :: coefficients(:)
      real(dp), allocatable :: ma_coefficients(:)
      type(order_search_result_t), allocatable :: order_fits(:)
      integer :: parameter_count = 0
      integer :: residual_start = 1
      integer :: info = 0
   end type candidate_fit_t

   type :: information_criteria_t
      real(dp) :: log_likelihood = -huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      real(dp) :: aicc = huge(1.0_dp)
      real(dp) :: bic = huge(1.0_dp)
      real(dp) :: innovation_variance = 0.0_dp
      integer :: observations = 0
      integer :: info = 0
   end type information_criteria_t

   type :: arma_screen_t
      real(dp), allocatable :: ar(:)
      real(dp), allocatable :: ma(:)
      integer :: info = 0
   end type arma_screen_t

   interface display
      module procedure display_series_profile
      module procedure display_automatic_model_result
   end interface display

   public :: profile_series, automatic_model, display

contains

   pure function profile_series(series, frequency, max_lag) result(out)
      !! Compute dependence, trend, seasonality, and volatility indicators.
      real(dp), intent(in) :: series(:) !! Finite univariate observations.
      integer, intent(in), optional :: frequency !! Observations per seasonal cycle.
      integer, intent(in), optional :: max_lag !! Maximum reported correlation lag.
      type(series_profile_t) :: out
      real(dp), allocatable :: squared(:)
      real(dp) :: centered_time, centered_value, covariance_value
      real(dp) :: time_variance, value_variance
      integer :: i, lag_limit, selected_frequency

      out%observations = size(series)
      if (out%observations < 3) then
         out%info = 1
         return
      end if
      if (any(.not. ieee_is_finite(series))) then
         out%info = 2
         return
      end if
      selected_frequency = 1
      if (present(frequency)) selected_frequency = frequency
      if (selected_frequency < 1) then
         out%info = 3
         return
      end if
      lag_limit = min(24, out%observations - 1)
      if (present(max_lag)) lag_limit = min(max_lag, out%observations - 1)
      if (lag_limit < 1) then
         out%info = 4
         return
      end if

      out%frequency = selected_frequency
      out%max_lag = lag_limit
      out%mean = sum(series)/real(out%observations, dp)
      out%variance = sum((series - out%mean)**2)/ &
         real(out%observations - 1, dp)
      out%minimum = minval(series)
      out%maximum = maxval(series)
      out%acf = acf_values(series, lag_limit)
      out%pacf = pacf_values(series, lag_limit)
      squared = (series - out%mean)**2
      out%squared_acf = acf_values(squared, lag_limit)
      out%dependence_threshold = 1.96_dp/sqrt(real(out%observations, dp))

      covariance_value = 0.0_dp
      time_variance = 0.0_dp
      value_variance = sum((series - out%mean)**2)
      do i = 1, out%observations
         centered_time = real(i, dp) - 0.5_dp*real(out%observations + 1, dp)
         centered_value = series(i) - out%mean
         covariance_value = covariance_value + centered_time*centered_value
         time_variance = time_variance + centered_time**2
      end do
      if (time_variance > 0.0_dp .and. value_variance > 0.0_dp) then
         out%trend_strength = abs(covariance_value)/ &
            sqrt(time_variance*value_variance)
      end if
      out%trend_detected = out%trend_strength > 0.35_dp
      out%autocorrelation_detected = &
         maxval(abs(out%acf(2:))) > out%dependence_threshold
      out%conditional_variance_detected = &
         maxval(abs(out%squared_acf(2:))) > out%dependence_threshold
      out%differencing_suggested = abs(out%acf(2)) > 0.8_dp
      if (selected_frequency > 1 .and. selected_frequency <= lag_limit) then
         out%seasonal_strength = abs(out%acf(selected_frequency + 1))
         out%seasonality_detected = &
            out%seasonal_strength > out%dependence_threshold
      end if
   end function profile_series

   function automatic_model(series, options) result(out)
      !! Profile a series, validate suitable candidates, and refit the winner.
      real(dp), intent(in) :: series(:) !! Finite univariate observations.
      type(automatic_model_options_t), intent(in), optional :: options !! Search configuration.
      type(automatic_model_result_t) :: out
      type(automatic_model_options_t) :: selected_options
      type(candidate_fit_t), allocatable :: full_fits(:)
      type(candidate_fit_t) :: fit
      integer, allocatable :: codes(:)
      real(dp), allocatable :: full_fit_times(:), validation_fit_times(:)
      character(len=10) :: criterion
      integer :: candidate_count, clock_end, clock_max, clock_rate, clock_start
      integer :: common_start, i, training_size, validation_size

      selected_options = automatic_model_options_t()
      if (present(options)) selected_options = options
      criterion = lower_ascii(trim(selected_options%selection))
      if (size(series) < 8 .or. selected_options%horizon < 1 .or. &
         selected_options%max_lag < 1 .or. selected_options%max_ar_order < 1) then
         out%info = 1
         return
      end if
      if (criterion /= selection_validation .and. criterion /= selection_aicc .and. &
         criterion /= selection_bic) then
         out%info = 5
         return
      end if
      out%selection_criterion = criterion
      out%fit_times_recorded = selected_options%time_fits
      out%profile = profile_series(series, selected_options%frequency, &
         selected_options%max_lag)
      if (out%profile%info /= 0) then
         out%info = 10 + out%profile%info
         return
      end if

      validation_size = 0
      if (criterion == selection_validation) then
         validation_size = selected_options%validation_size
         if (validation_size <= 0) then
            validation_size = max(4, min(24, size(series)/5))
         end if
         validation_size = min(validation_size, size(series) - 6)
         if (selected_options%frequency > 1) then
            validation_size = min(validation_size, &
               max(1, size(series) - 2*selected_options%frequency))
         end if
         if (validation_size < 1) then
            out%info = 2
            return
         end if
      end if
      training_size = size(series) - validation_size
      out%validation_size = validation_size
      if (validation_size > 0) then
         out%validation_actual = series(training_size + 1:)
      else
         allocate(out%validation_actual(0))
      end if

      allocate(codes(10))
      candidate_count = 0
      call append_code(codes, candidate_count, model_mean)
      call append_code(codes, candidate_count, model_naive)
      call append_code(codes, candidate_count, model_drift)
      call append_code(codes, candidate_count, model_ses)
      call append_code(codes, candidate_count, model_theta)
      if (out%profile%trend_detected) then
         call append_code(codes, candidate_count, model_holt)
      end if
      if (out%profile%autocorrelation_detected) then
         call append_code(codes, candidate_count, model_ar)
         if (.not. out%profile%differencing_suggested) then
            call append_code(codes, candidate_count, model_arma)
         end if
      end if
      if (selected_options%include_seasonal .and. &
         out%profile%seasonality_detected .and. &
         training_size >= 2*selected_options%frequency) then
         call append_code(codes, candidate_count, model_seasonal_naive)
         call append_code(codes, candidate_count, model_holt_winters)
      end if

      allocate(out%candidates(candidate_count), full_fits(candidate_count), &
         full_fit_times(candidate_count), validation_fit_times(candidate_count))
      full_fit_times = 0.0_dp
      validation_fit_times = 0.0_dp
      do i = 1, candidate_count
         if (selected_options%time_fits) then
            call system_clock(clock_start, clock_rate, clock_max)
         end if
         full_fits(i) = fit_candidate(series, selected_options%horizon, codes(i), &
            selected_options%frequency, selected_options%max_ar_order, criterion)
         if (selected_options%time_fits) then
            call system_clock(clock_end)
            full_fit_times(i) = clock_seconds(clock_start, &
               clock_end, clock_rate, clock_max)
         end if
      end do
      common_start = maxval(full_fits%residual_start)
      do i = 1, candidate_count
         out%candidates(i) = candidate_summary(codes(i), full_fits(i), &
            out%profile, common_start)
         out%candidates(i)%full_fit_seconds = full_fit_times(i)
         if (selected_options%time_fits) then
            call system_clock(clock_start, clock_rate, clock_max)
         end if
         if (criterion == selection_validation) then
            fit = fit_candidate(series(:training_size), validation_size, codes(i), &
               selected_options%frequency, selected_options%max_ar_order, criterion)
            out%candidates(i) = add_validation_scores(out%candidates(i), &
               fit%forecast, out%validation_actual, fit%info)
         end if
         if (selected_options%time_fits .and. &
            criterion == selection_validation) then
            call system_clock(clock_end)
            validation_fit_times(i) = clock_seconds( &
               clock_start, clock_end, clock_rate, clock_max)
            out%candidates(i)%validation_fit_seconds = validation_fit_times(i)
         end if
      end do
      call sort_candidates(out%candidates, criterion)
      if (.not. out%candidates(1)%converged) then
         out%info = 3
         return
      end if
      out%candidates(1)%selected = .true.
      out%selected = 1
      do i = 1, candidate_count
         if (codes(i) == out%candidates(1)%model_code) exit
      end do
      fit = full_fits(i)
      if (fit%info /= 0) then
         out%info = 4
         return
      end if
      out%forecast = fit%forecast
      out%fitted = fit%fitted
      out%residuals = fit%residuals
   end function automatic_model

   pure subroutine append_code(codes, count, code)
      !! Append one internal model code to a preallocated candidate vector.
      integer, intent(inout) :: codes(:) !! Candidate model codes.
      integer, intent(inout) :: count !! Number of populated codes.
      integer, intent(in) :: code !! Model code to append.

      count = count + 1
      codes(count) = code
   end subroutine append_code

   pure function fit_candidate(series, horizon, code, frequency, max_ar_order, &
      criterion) &
      result(out)
      !! Dispatch one candidate family to its existing library implementation.
      real(dp), intent(in) :: series(:) !! Training observations.
      integer, intent(in) :: horizon !! Forecast horizon.
      integer, intent(in) :: code !! Internal candidate model code.
      integer, intent(in) :: frequency !! Observations per seasonal cycle.
      integer, intent(in) :: max_ar_order !! Maximum autoregressive order.
      character(len=*), intent(in) :: criterion !! Criterion used to select AR order.
      type(candidate_fit_t) :: out
      type(forecast_result_t) :: base

      select case (code)
      case (model_mean)
         base = meanf(series, horizon)
         out = from_forecast_result(base, 2, 1)
         out%level = sum(series)/real(size(series), dp)
      case (model_naive)
         base = naive(series, horizon)
         out = from_forecast_result(base, 1, 2)
         out%level = series(size(series))
      case (model_drift)
         base = rwf(series, horizon, .true.)
         out = from_forecast_result(base, 2, 2)
         out%level = series(size(series))
         out%drift = (series(size(series)) - series(1))/ &
            real(size(series) - 1, dp)
      case (model_ses)
         out = fit_ses(series, horizon)
      case (model_holt)
         out = fit_holt(series, horizon)
      case (model_theta)
         out = fit_theta(series, horizon)
      case (model_seasonal_naive)
         if (frequency < 2 .or. size(series) < frequency) then
            out%info = 1
         else
            base = snaive(series, frequency, horizon)
            out = from_forecast_result(base, 1, frequency + 1)
         end if
      case (model_holt_winters)
         out = fit_holt_winters(series, frequency, horizon)
      case (model_ar)
         out = fit_autoregression(series, horizon, max_ar_order, criterion)
      case (model_arma)
         out = fit_arma(series, horizon, criterion)
      case default
         out%info = 2
      end select
   end function fit_candidate

   pure function from_forecast_result(base, parameter_count, residual_start) result(out)
      !! Convert the shared forecast result into an automatic candidate fit.
      type(forecast_result_t), intent(in) :: base !! Existing forecast result.
      integer, intent(in) :: parameter_count !! Number of estimated parameters.
      integer, intent(in) :: residual_start !! First residual eligible for scoring.
      type(candidate_fit_t) :: out

      out%forecast = base%mean
      out%fitted = base%fitted
      out%residuals = base%residuals
      out%parameter_count = parameter_count
      out%residual_start = residual_start
      if (allocated(out%forecast) .and. all(ieee_is_finite(out%forecast))) then
         out%info = 0
      else
         out%info = 1
      end if
   end function from_forecast_result

   pure function fit_ses(series, horizon) result(out)
      !! Select simple exponential smoothing by in-sample squared error.
      real(dp), intent(in) :: series(:) !! Training observations.
      integer, intent(in) :: horizon !! Forecast horizon.
      type(candidate_fit_t) :: out
      type(forecast_result_t) :: current, best
      real(dp) :: alpha, score, best_score
      integer :: i

      best_score = huge(1.0_dp)
      do i = 1, 19
         alpha = 0.05_dp*real(i, dp)
         current = ses(series, horizon, alpha)
         score = residual_sse(current%residuals)
         if (score < best_score) then
            best_score = score
            best = current
            out%alpha = alpha
         end if
      end do
      out%forecast = best%mean
      out%fitted = best%fitted
      out%residuals = best%residuals
      out%parameter_count = 2
      out%residual_start = 2
      out%level = best%mean(1)
   end function fit_ses

   pure function fit_holt(series, horizon) result(out)
      !! Select Holt smoothing weights by in-sample squared error.
      real(dp), intent(in) :: series(:) !! Training observations.
      integer, intent(in) :: horizon !! Forecast horizon.
      type(candidate_fit_t) :: out
      type(forecast_result_t) :: current, best
      real(dp) :: alpha, beta, score, best_score
      integer :: i, j

      best_score = huge(1.0_dp)
      do i = 1, 9
         alpha = 0.1_dp*real(i, dp)
         do j = 1, 9
            beta = 0.1_dp*real(j, dp)
            current = holt(series, horizon, alpha, beta)
            score = residual_sse(current%residuals)
            if (score < best_score) then
               best_score = score
               best = current
               out%alpha = alpha
               out%beta = beta
            end if
         end do
      end do
      out%forecast = best%mean
      out%fitted = best%fitted
      out%residuals = best%residuals
      out%parameter_count = 4
      out%residual_start = 2
      out%level = best%mean(1)
   end function fit_holt

   pure function fit_theta(series, horizon) result(out)
      !! Select the Theta smoothing weight by in-sample squared error.
      real(dp), intent(in) :: series(:) !! Training observations.
      integer, intent(in) :: horizon !! Forecast horizon.
      type(candidate_fit_t) :: out
      type(forecast_result_t) :: current, best
      real(dp) :: alpha, score, best_score
      integer :: i

      best_score = huge(1.0_dp)
      do i = 1, 19
         alpha = 0.05_dp*real(i, dp)
         current = thetaf(series, horizon, alpha)
         score = residual_sse(current%residuals)
         if (score < best_score) then
            best_score = score
            best = current
            out%alpha = alpha
         end if
      end do
      out%forecast = best%mean
      out%fitted = best%fitted
      out%residuals = best%residuals
      out%parameter_count = 3
      out%residual_start = 2
      out%level = best%mean(1)
   end function fit_theta

   pure function fit_holt_winters(series, frequency, horizon) result(out)
      !! Select additive Holt-Winters weights on a compact parameter grid.
      real(dp), intent(in) :: series(:) !! Training observations.
      integer, intent(in) :: frequency !! Observations per seasonal cycle.
      integer, intent(in) :: horizon !! Forecast horizon.
      type(candidate_fit_t) :: out
      type(forecast_result_t) :: current, best
      real(dp), parameter :: grid(3) = [0.2_dp, 0.5_dp, 0.8_dp]
      real(dp) :: score, best_score
      integer :: i, j, k

      if (frequency < 2 .or. size(series) < 2*frequency) then
         out%info = 1
         return
      end if
      best_score = huge(1.0_dp)
      do i = 1, size(grid)
         do j = 1, size(grid)
            do k = 1, size(grid)
               current = holt_winters(series, frequency, horizon, grid(i), &
                  grid(j), grid(k))
               score = residual_sse(current%residuals)
               if (score < best_score) then
                  best_score = score
                  best = current
                  out%alpha = grid(i)
                  out%beta = grid(j)
                  out%gamma = grid(k)
               end if
            end do
         end do
      end do
      out%forecast = best%mean
      out%fitted = best%fitted
      out%residuals = best%residuals
      out%parameter_count = frequency + 4
      out%residual_start = frequency + 1
      out%level = best%mean(1)
   end function fit_holt_winters

   pure function fit_autoregression(series, horizon, max_order, criterion) result(out)
      !! Select a Yule-Walker AR order and recursively forecast about the mean.
      real(dp), intent(in) :: series(:) !! Training observations.
      integer, intent(in) :: horizon !! Forecast horizon.
      integer, intent(in) :: max_order !! Maximum autoregressive order.
      character(len=*), intent(in) :: criterion !! Criterion used to select the order.
      type(candidate_fit_t) :: out
      type(yule_walker_result_t) :: current, best
      type(information_criteria_t) :: order_criteria
      real(dp), allocatable :: history(:)
      real(dp) :: mean_value, best_criterion, order_score, prediction
      integer :: order, selected_max, i, lag

      selected_max = min(max_order, max(1, size(series)/3))
      allocate(out%order_fits(selected_max))
      best_criterion = huge(1.0_dp)
      do order = 1, selected_max
         out%order_fits(order)%ar_order = order
         current = yule_walker_fit(series, order)
         out%order_fits(order)%info = current%info
         order_score = huge(1.0_dp)
         if (current%info == 0) then
            order_criteria = autoregression_information_criteria(series, &
               current%coefficients, order + 2)
            out%order_fits(order)%info = order_criteria%info
            if (order_criteria%info == 0) then
               out%order_fits(order)%mean = sum(series)/real(size(series), dp)
               out%order_fits(order)%innovation_variance = &
                  order_criteria%innovation_variance
               out%order_fits(order)%log_likelihood = order_criteria%log_likelihood
               out%order_fits(order)%aic = order_criteria%aic
               out%order_fits(order)%aicc = order_criteria%aicc
               out%order_fits(order)%bic = order_criteria%bic
               out%order_fits(order)%ar_coefficients = current%coefficients
               out%order_fits(order)%converged = .true.
            end if
            select case (trim(criterion))
            case (selection_aicc, selection_bic)
               if (trim(criterion) == selection_aicc) then
                  order_score = order_criteria%aicc
               else
                  order_score = order_criteria%bic
               end if
            case default
               order_score = current%criterion
            end select
         end if
         out%order_fits(order)%search_score = order_score
         if (out%order_fits(order)%converged .and. &
            order_score < best_criterion) then
            best = current
            best_criterion = order_score
            out%order = order
         end if
      end do
      if (out%order == 0) then
         out%info = 1
         return
      end if
      mean_value = sum(series)/real(size(series), dp)
      allocate(out%fitted(size(series)), out%residuals(size(series)))
      out%fitted = ieee_value(0.0_dp, ieee_quiet_nan)
      do i = out%order + 1, size(series)
         prediction = mean_value
         do lag = 1, out%order
            prediction = prediction + best%coefficients(lag)* &
               (series(i - lag) - mean_value)
         end do
         out%fitted(i) = prediction
      end do
      out%residuals = series - out%fitted
      allocate(history(size(series) + horizon), out%forecast(horizon))
      history(:size(series)) = series
      do i = 1, horizon
         prediction = mean_value
         do lag = 1, out%order
            prediction = prediction + best%coefficients(lag)* &
               (history(size(series) + i - lag) - mean_value)
         end do
         history(size(series) + i) = prediction
         out%forecast(i) = prediction
      end do
      out%parameter_count = out%order + 2
      out%residual_start = out%order + 1
      out%level = mean_value
      out%coefficients = best%coefficients
   end function fit_autoregression

   pure function autoregression_information_criteria(series, coefficients, &
      parameter_count) result(out)
      !! Score one autoregressive order from its conditional residuals.
      real(dp), intent(in) :: series(:) !! Training observations.
      real(dp), intent(in) :: coefficients(:) !! Autoregressive coefficients.
      integer, intent(in) :: parameter_count !! Mean, AR, and variance parameters.
      type(information_criteria_t) :: out
      real(dp), allocatable :: residuals(:)
      real(dp) :: mean_value, prediction
      integer :: i, lag, order

      order = size(coefficients)
      if (order < 1 .or. size(series) <= order) then
         out%info = 1
         return
      end if
      mean_value = sum(series)/real(size(series), dp)
      allocate(residuals(size(series)))
      residuals = ieee_value(0.0_dp, ieee_quiet_nan)
      do i = order + 1, size(series)
         prediction = mean_value
         do lag = 1, order
            prediction = prediction + coefficients(lag)* &
               (series(i - lag) - mean_value)
         end do
         residuals(i) = series(i) - prediction
      end do
      out = gaussian_information_criteria(residuals, order + 1, parameter_count)
   end function autoregression_information_criteria

   pure function fit_arma(series, horizon, criterion) result(out)
      !! Screen an ARMA order grid and refine manageable samples by exact MLE.
      real(dp), intent(in) :: series(:) !! Training observations.
      integer, intent(in) :: horizon !! Forecast horizon.
      character(len=*), intent(in) :: criterion !! Criterion used to select the order pair.
      type(candidate_fit_t) :: out
      integer, parameter :: exact_refinement_limit = 2000
      type(astsa_sarima_fit_t) :: current, best
      type(astsa_sarima_forecast_t) :: prediction
      type(astsa_sarima_likelihood_t) :: likelihood, best_likelihood
      type(arma_screen_t) :: screened, best_screened
      real(dp), allocatable :: extended(:), extended_residuals(:)
      real(dp) :: best_score, mean_value, predicted, score
      integer :: ar_order, attempt, lag, ma_order, step, time, trial

      best_score = huge(1.0_dp)
      mean_value = sum(series)/real(size(series), dp)
      allocate(out%order_fits(9))
      trial = 0
      do ar_order = 0, 3
         do ma_order = 1, 3
            if (ar_order + ma_order > 4) cycle
            trial = trial + 1
            out%order_fits(trial)%ar_order = ar_order
            out%order_fits(trial)%ma_order = ma_order
            out%order_fits(trial)%mean = mean_value
            if (size(series) <= exact_refinement_limit) then
               current = sarima_fit(series, ar_order, 0, ma_order, &
                  0, 0, 0, 1, include_intercept=.true., &
                  max_iterations=100, tolerance=1.0e-5_dp, &
                  exact_likelihood=.true.)
               out%order_fits(trial)%info = current%likelihood%info
               if (current%likelihood%info /= 0) cycle
               likelihood = current%likelihood
            else
               screened = screen_arma_hannan_rissanen(series, ar_order, ma_order)
               out%order_fits(trial)%info = screened%info
               if (screened%info /= 0) cycle
               do attempt = 0, 12
                  likelihood = sarima_likelihood(series, ar=screened%ar, &
                     ma=screened%ma, intercept=mean_value)
                  if (likelihood%info == 0) exit
                  screened%ar = 0.8_dp*screened%ar
                  screened%ma = 0.8_dp*screened%ma
               end do
               out%order_fits(trial)%info = likelihood%info
               if (likelihood%info /= 0) cycle
            end if
            out%order_fits(trial)%innovation_variance = likelihood%sigma2
            out%order_fits(trial)%log_likelihood = likelihood%log_likelihood
            out%order_fits(trial)%aic = likelihood%aic
            out%order_fits(trial)%aicc = likelihood%aicc
            out%order_fits(trial)%bic = likelihood%bic
            if (size(series) <= exact_refinement_limit) then
               out%order_fits(trial)%ar_coefficients = &
                  current%coefficients(:ar_order)
               out%order_fits(trial)%ma_coefficients = &
                  current%coefficients(ar_order + 1:ar_order + ma_order)
               out%order_fits(trial)%mean = &
                  current%coefficients(ar_order + ma_order + 1)
            else
               out%order_fits(trial)%ar_coefficients = screened%ar
               out%order_fits(trial)%ma_coefficients = screened%ma
            end if
            if (trim(criterion) == selection_bic) then
               score = likelihood%bic
            else
               score = likelihood%aicc
            end if
            out%order_fits(trial)%search_score = score
            if (.not. ieee_is_finite(score)) then
               out%order_fits(trial)%info = 90
               cycle
            end if
            out%order_fits(trial)%converged = .true.
            if (score < best_score) then
               if (size(series) <= exact_refinement_limit) then
                  best = current
                  best%info = 0
               else
                  best_screened = screened
                  best_likelihood = likelihood
               end if
               best_score = score
            end if
         end do
      end do
      if (best_score >= huge(1.0_dp)) then
         out%info = 1
         return
      end if
      if (size(series) > exact_refinement_limit) then
         out%order = size(best_screened%ar)
         out%ma_order = size(best_screened%ma)
         out%parameter_count = out%order + out%ma_order + 2
         out%residual_start = max(out%order, out%ma_order) + 1
         out%coefficients = best_screened%ar
         out%ma_coefficients = best_screened%ma
         out%level = mean_value
         out%residuals = best_likelihood%residuals
         out%fitted = series - out%residuals
         allocate(extended(size(series) + horizon), &
            extended_residuals(size(series) + horizon), out%forecast(horizon))
         extended(:size(series)) = series - mean_value
         extended_residuals(:size(series)) = out%residuals
         extended_residuals(size(series) + 1:) = 0.0_dp
         do step = 1, horizon
            time = size(series) + step
            predicted = 0.0_dp
            do lag = 1, min(out%order, time - 1)
               predicted = predicted + out%coefficients(lag)* &
                  extended(time - lag)
            end do
            do lag = 1, min(out%ma_order, time - 1)
               predicted = predicted + out%ma_coefficients(lag)* &
                  extended_residuals(time - lag)
            end do
            extended(time) = predicted
            out%forecast(step) = mean_value + predicted
         end do
         return
      end if
      prediction = sarima_forecast(best, series, horizon)
      if (prediction%info /= 0 .or. .not. allocated(prediction%mean)) then
         out%info = 2
         return
      end if
      out%order = best%p
      out%ma_order = best%q
      out%parameter_count = best%p + best%q + 2
      out%residual_start = max(best%p, best%q) + 1
      out%coefficients = best%coefficients(:best%p)
      out%ma_coefficients = best%coefficients(best%p + 1:best%p + best%q)
      out%level = best%coefficients(best%p + best%q + 1)
      out%forecast = prediction%mean
      out%residuals = best%likelihood%residuals
      out%fitted = series - out%residuals
   end function fit_arma

   pure function screen_arma_hannan_rissanen(series, ar_order, ma_order) &
      result(out)
      !! Estimate ARMA coefficients by regression without exact diagnostics.
      real(dp), intent(in) :: series(:) !! Training observations.
      integer, intent(in) :: ar_order !! Autoregressive order.
      integer, intent(in) :: ma_order !! Moving-average order.
      type(arma_screen_t) :: out
      type(yule_walker_result_t) :: preliminary
      real(dp), allocatable :: centered(:), proxy(:), design(:, :)
      real(dp), allocatable :: response(:), coefficients(:)
      real(dp), allocatable :: standard_errors(:), residuals(:)
      real(dp) :: mean_value, rss
      integer :: columns, lag, offset, preliminary_order, row, rows, status, time

      preliminary_order = 20 + ar_order + ma_order
      offset = max(ar_order, ma_order)
      rows = size(series) - preliminary_order - offset
      columns = ar_order + ma_order
      if (ar_order < 0 .or. ma_order < 1 .or. rows <= columns .or. &
         preliminary_order >= size(series)) then
         out%info = 1
         return
      end if
      mean_value = sum(series)/real(size(series), dp)
      centered = series - mean_value
      preliminary = yule_walker_fit(centered, preliminary_order)
      if (preliminary%info /= 0) then
         out%info = 10 + preliminary%info
         return
      end if
      allocate(proxy(size(series)))
      proxy = 0.0_dp
      do time = preliminary_order + 1, size(series)
         proxy(time) = centered(time)
         do lag = 1, preliminary_order
            proxy(time) = proxy(time) - &
               preliminary%coefficients(lag)*centered(time - lag)
         end do
      end do
      allocate(design(rows, columns), response(rows))
      do row = 1, rows
         time = preliminary_order + offset + row
         response(row) = centered(time)
         do lag = 1, ar_order
            design(row, lag) = centered(time - lag)
         end do
         do lag = 1, ma_order
            design(row, ar_order + lag) = proxy(time - lag)
         end do
      end do
      call ols_fit(design, response, coefficients, standard_errors, &
         residuals, rss, status)
      if (status /= 0) then
         out%info = 20 + status
         return
      end if
      allocate(out%ar(ar_order), out%ma(ma_order))
      if (ar_order > 0) out%ar = coefficients(:ar_order)
      out%ma = coefficients(ar_order + 1:)
   end function screen_arma_hannan_rissanen

   pure function candidate_summary(code, fit, profile, scoring_start) result(out)
      !! Combine one full-data fit with information criteria and rationale.
      integer, intent(in) :: code !! Internal candidate model code.
      type(candidate_fit_t), intent(in) :: fit !! Candidate fit and forecasts.
      type(series_profile_t), intent(in) :: profile !! Series statistical profile.
      integer, intent(in) :: scoring_start !! Common first residual used by every model.
      type(candidate_model_result_t) :: out
      type(information_criteria_t) :: criteria

      out%model_code = code
      out%name = model_name(code)
      out%rationale = model_rationale(code, profile)
      out%order = fit%order
      out%ma_order = fit%ma_order
      out%alpha = fit%alpha
      out%beta = fit%beta
      out%gamma = fit%gamma
      out%level = fit%level
      out%drift = fit%drift
      if (allocated(fit%coefficients)) out%coefficients = fit%coefficients
      if (allocated(fit%ma_coefficients)) then
         out%ma_coefficients = fit%ma_coefficients
      end if
      if (allocated(fit%order_fits)) out%order_fits = fit%order_fits
      out%parameter_count = fit%parameter_count
      out%info = fit%info
      if (fit%info /= 0 .or. .not. allocated(fit%forecast) .or. &
         .not. allocated(fit%residuals)) return
      if (any(.not. ieee_is_finite(fit%forecast))) then
         out%info = 2
         return
      end if
      criteria = gaussian_information_criteria(fit%residuals, scoring_start, &
         fit%parameter_count)
      if (criteria%info /= 0) then
         out%info = 3
         return
      end if
      out%log_likelihood = criteria%log_likelihood
      out%aic = criteria%aic
      out%aicc = criteria%aicc
      out%bic = criteria%bic
      out%innovation_variance = criteria%innovation_variance
      out%effective_observations = criteria%observations
      out%converged = .true.
   end function candidate_summary

   pure function add_validation_scores(candidate, forecast, actual, fit_info) result(out)
      !! Add holdout errors to a full-data candidate summary.
      type(candidate_model_result_t), intent(in) :: candidate !! Full-data candidate summary.
      real(dp), intent(in), allocatable :: forecast(:) !! Forecasts for the held-out tail.
      real(dp), intent(in) :: actual(:) !! Held-out observations.
      integer, intent(in) :: fit_info !! Validation-fit status code.
      type(candidate_model_result_t) :: out
      real(dp), allocatable :: error(:)

      out = candidate
      if (.not. out%converged .or. fit_info /= 0 .or. .not. allocated(forecast)) then
         out%converged = .false.
         out%info = 4
         return
      end if
      if (size(forecast) /= size(actual) .or. &
         any(.not. ieee_is_finite(forecast))) then
         out%converged = .false.
         out%info = 5
         return
      end if
      error = actual - forecast
      out%rmse = sqrt(sum(error**2)/real(size(error), dp))
      out%mae = sum(abs(error))/real(size(error), dp)
   end function add_validation_scores

   pure function model_name(code) result(name)
      !! Return the display name for one internal model code.
      integer, intent(in) :: code !! Internal candidate model code.
      character(len=:), allocatable :: name

      select case (code)
      case (model_mean)
         name = "Mean"
      case (model_naive)
         name = "Naive"
      case (model_drift)
         name = "Random walk with drift"
      case (model_ses)
         name = "Simple exponential smoothing"
      case (model_holt)
         name = "Holt trend"
      case (model_theta)
         name = "Theta"
      case (model_seasonal_naive)
         name = "Seasonal naive"
      case (model_holt_winters)
         name = "Additive Holt-Winters"
      case (model_ar)
         name = "Yule-Walker autoregression"
      case (model_arma)
         name = "Gaussian ARMA"
      case default
         name = "Unknown"
      end select
   end function model_name

   pure function model_rationale(code, profile) result(rationale)
      !! Explain why one model family entered the candidate set.
      integer, intent(in) :: code !! Internal candidate model code.
      type(series_profile_t), intent(in) :: profile !! Series statistical profile.
      character(len=:), allocatable :: rationale

      select case (code)
      case (model_mean, model_naive, model_drift)
         rationale = "Benchmark model retained for forecast comparison."
      case (model_ses, model_theta)
         rationale = "General level forecast suitable for nonseasonal data."
      case (model_holt)
         rationale = "Trend correlation exceeded the profiling threshold."
      case (model_seasonal_naive, model_holt_winters)
         rationale = "Correlation at the supplied seasonal frequency was significant."
      case (model_ar)
         rationale = "At least one level autocorrelation exceeded the sampling threshold."
      case (model_arma)
         rationale = "Stationary level dependence supports joint AR and MA dynamics."
      case default
         rationale = "No rationale is available."
      end select
   end function model_rationale

   pure subroutine sort_candidates(candidates, criterion)
      !! Sort candidate results by the requested selection criterion.
      type(candidate_model_result_t), intent(inout) :: candidates(:) !! Candidate results to rank.
      character(len=*), intent(in) :: criterion !! Validation, AICc, or BIC ranking mode.
      type(candidate_model_result_t) :: held
      integer :: i, j

      do i = 2, size(candidates)
         held = candidates(i)
         j = i - 1
         do while (j >= 1)
            if (.not. candidate_precedes(held, candidates(j), criterion)) exit
            candidates(j + 1) = candidates(j)
            j = j - 1
         end do
         candidates(j + 1) = held
      end do
   end subroutine sort_candidates

   pure logical function candidate_precedes(first, second, criterion) result(precedes)
      !! Compare two candidates using the requested ranking statistic.
      type(candidate_model_result_t), intent(in) :: first !! Candidate on the left.
      type(candidate_model_result_t), intent(in) :: second !! Candidate on the right.
      character(len=*), intent(in) :: criterion !! Validation, AICc, or BIC ranking mode.

      if (first%converged .neqv. second%converged) then
         precedes = first%converged
      else
         select case (trim(criterion))
         case (selection_aicc)
            if (first%aicc /= second%aicc) then
               precedes = first%aicc < second%aicc
            else
               precedes = first%bic < second%bic
            end if
         case (selection_bic)
            if (first%bic /= second%bic) then
               precedes = first%bic < second%bic
            else
               precedes = first%aicc < second%aicc
            end if
         case default
            if (first%rmse /= second%rmse) then
               precedes = first%rmse < second%rmse
            else
               precedes = first%mae < second%mae
            end if
         end select
      end if
   end function candidate_precedes

   pure real(dp) function residual_sse(residuals) result(value)
      !! Sum squared finite residuals and reject an empty residual sample.
      real(dp), intent(in) :: residuals(:) !! Model residuals.
      integer :: i, count

      value = 0.0_dp
      count = 0
      do i = 1, size(residuals)
         if (ieee_is_finite(residuals(i))) then
            value = value + residuals(i)**2
            count = count + 1
         end if
      end do
      if (count == 0) value = huge(1.0_dp)
   end function residual_sse

   pure function gaussian_information_criteria(residuals, first, &
      parameter_count) result(out)
      !! Compute Gaussian log likelihood, AIC, AICc, and BIC.
      real(dp), intent(in) :: residuals(:) !! Model residuals.
      integer, intent(in) :: first !! First residual included in every model score.
      integer, intent(in) :: parameter_count !! Estimated parameter count.
      type(information_criteria_t) :: out
      real(dp) :: pi, sse, variance
      integer :: count

      if (first < 1 .or. first > size(residuals) .or. parameter_count < 1) then
         out%info = 1
         return
      end if
      if (any(.not. ieee_is_finite(residuals(first:)))) then
         out%info = 2
         return
      end if
      count = size(residuals) - first + 1
      sse = sum(residuals(first:)**2)
      pi = acos(-1.0_dp)
      variance = max(sse/real(count, dp), tiny(1.0_dp))
      out%innovation_variance = variance
      out%observations = count
      out%log_likelihood = -0.5_dp*real(count, dp)* &
         (log(2.0_dp*pi*variance) + 1.0_dp)
      out%aic = -2.0_dp*out%log_likelihood + &
         2.0_dp*real(parameter_count, dp)
      if (count > parameter_count + 1) then
         out%aicc = out%aic + &
            2.0_dp*real(parameter_count*(parameter_count + 1), dp)/ &
            real(count - parameter_count - 1, dp)
      end if
      out%bic = -2.0_dp*out%log_likelihood + &
         real(parameter_count, dp)*log(real(count, dp))
   end function gaussian_information_criteria

   pure function lower_ascii(value) result(lower)
      !! Convert ASCII letters to lowercase for option matching.
      character(len=*), intent(in) :: value !! Text to convert.
      character(len=len(value)) :: lower
      integer :: code, i

      lower = value
      do i = 1, len(value)
         code = iachar(value(i:i))
         if (code >= iachar("A") .and. code <= iachar("Z")) then
            lower(i:i) = achar(code + iachar("a") - iachar("A"))
         end if
      end do
   end function lower_ascii

   pure function clock_seconds(start_count, end_count, count_rate, &
      maximum_count) result(seconds)
      !! Convert system-clock counts to elapsed seconds, including one wrap.
      integer, intent(in) :: start_count !! System-clock count before an operation.
      integer, intent(in) :: end_count !! System-clock count after an operation.
      integer, intent(in) :: count_rate !! System-clock counts per second.
      integer, intent(in) :: maximum_count !! Largest count before the clock wraps.
      real(dp) :: seconds

      if (count_rate <= 0) then
         seconds = 0.0_dp
      else if (end_count >= start_count) then
         seconds = real(end_count - start_count, dp)/real(count_rate, dp)
      else
         seconds = (real(maximum_count - start_count, dp) + &
            real(end_count, dp) + 1.0_dp)/real(count_rate, dp)
      end if
   end function clock_seconds

   subroutine display_series_profile(profile, display_lags)
      !! Display the properties used by automatic candidate gating.
      type(series_profile_t), intent(in) :: profile !! Series statistical profile.
      integer, intent(in), optional :: display_lags !! Maximum positive lag to print.
      integer :: lag, lag_limit

      write(*, '(a)') "Univariate series profile"
      write(*, '(a,i0)') "  status: ", profile%info
      write(*, '(a,i0)') "  observations: ", profile%observations
      if (profile%info /= 0) return
      write(*, '(a,i0)') "  frequency: ", profile%frequency
      write(*, '(a,es14.6)') "  mean: ", profile%mean
      write(*, '(a,es14.6)') "  variance: ", profile%variance
      write(*, '(a,f8.4)') "  trend strength: ", profile%trend_strength
      write(*, '(a,f8.4)') "  seasonal strength: ", profile%seasonal_strength
      write(*, '(a,l1)') "  trend detected: ", profile%trend_detected
      write(*, '(a,l1)') "  seasonality detected: ", profile%seasonality_detected
      write(*, '(a,l1)') "  level dependence detected: ", &
         profile%autocorrelation_detected
      write(*, '(a,l1)') "  squared-value dependence detected: ", &
         profile%conditional_variance_detected
      write(*, '(a,l1)') "  differencing suggested: ", &
         profile%differencing_suggested
      lag_limit = min(5, profile%max_lag)
      if (present(display_lags)) lag_limit = min(display_lags, profile%max_lag)
      write(*, '(a)') ""
      write(*, '(a)') "  lag          ACF     squared ACF          PACF"
      do lag = 1, lag_limit
         write(*, '(2x,i3,3f16.6)') lag, profile%acf(lag + 1), &
            profile%squared_acf(lag + 1), profile%pacf(lag)
      end do
   end subroutine display_series_profile

   subroutine display_automatic_model_result(result, display_lags, &
      print_parameters, max_models, print_all_ar, print_all_arma)
      !! Display a profile, ranked candidates, and selected-model forecasts.
      type(automatic_model_result_t), intent(in) :: result !! Automatic modeling result.
      integer, intent(in), optional :: display_lags !! Maximum positive lag to print.
      logical, intent(in), optional :: print_parameters !! Print fitted model parameters.
      integer, intent(in), optional :: max_models !! Maximum candidate summaries, or zero for all.
      logical, intent(in), optional :: print_all_ar !! Print every tested AR order.
      logical, intent(in), optional :: print_all_arma !! Print every tested ARMA order.
      integer :: display_count, i
      logical :: show_all_ar, show_all_arma, show_parameters

      show_parameters = .false.
      if (present(print_parameters)) show_parameters = print_parameters
      show_all_ar = .false.
      if (present(print_all_ar)) show_all_ar = print_all_ar
      show_all_arma = .false.
      if (present(print_all_arma)) show_all_arma = print_all_arma

      if (result%profile%observations > 0) then
         call display_series_profile(result%profile, display_lags)
         write(*, '(a)') ""
      end if
      write(*, '(a)') "Automatic model comparison"
      write(*, '(a,i0)') "  status: ", result%info
      if (result%info /= 0) return
      write(*, '(a,a)') "  selection criterion: ", &
         trim(result%selection_criterion)
      if (result%selection_criterion == selection_validation) then
         write(*, '(a,i0)') "  validation observations: ", &
            result%validation_size
      else
         write(*, '(a)') "  models fitted to the full data set"
      end if
      display_count = size(result%candidates)
      if (present(max_models)) then
         if (max_models > 0) display_count = min(max_models, display_count)
      end if
      do i = 1, display_count
         call display_candidate_summary(result, i, show_parameters)
      end do
      if (result%selected > display_count) then
         write(*, '(a)') "Selected candidate outside the displayed ranking"
         call display_candidate_summary(result, result%selected, show_parameters)
      end if
      if (show_all_ar) then
         call display_order_search(result, model_ar, show_parameters)
      end if
      if (show_all_arma) then
         call display_order_search(result, model_arma, show_parameters)
      end if
      write(*, '(a,a)') "Selected model: ", result%candidates(result%selected)%name
      if (result%candidates(result%selected)%model_code == model_arma) then
         write(*, '(a,i0,a,i0,a)') "Selected ARMA order: (", &
            result%candidates(result%selected)%order, ",", &
            result%candidates(result%selected)%ma_order, ")"
      else if (result%candidates(result%selected)%order > 0) then
         write(*, '(a,i0)') "Selected AR order: ", &
            result%candidates(result%selected)%order
      end if
      if (result%profile%conditional_variance_detected) then
         write(*, '(a)') "Variance note: squared-value dependence suggests a conditional variance model."
      end if
      write(*, '(a)') "Forecasts"
      do i = 1, size(result%forecast)
         write(*, '(2x,i5,2x,es16.8)') i, result%forecast(i)
      end do
   end subroutine display_automatic_model_result

   subroutine display_order_search(result, model_code, print_parameters)
      !! Display all order fits considered within an AR or ARMA candidate.
      type(automatic_model_result_t), intent(in) :: result !! Automatic modeling result.
      integer, intent(in) :: model_code !! AR or ARMA internal model code.
      logical, intent(in) :: print_parameters !! Print coefficients for successful fits.
      integer :: candidate_index, i
      logical :: selected_order

      candidate_index = 0
      do i = 1, size(result%candidates)
         if (result%candidates(i)%model_code == model_code) then
            candidate_index = i
            exit
         end if
      end do
      write(*, '(a)') ""
      if (model_code == model_ar) then
         write(*, '(a)') "AR order search"
      else
         write(*, '(a)') "ARMA order search"
      end if
      if (candidate_index == 0 .or. &
         .not. allocated(result%candidates(candidate_index)%order_fits)) then
         write(*, '(2x,a)') "This model class was not fitted."
         return
      end if
      if (model_code == model_ar) then
         write(*, '(2x,a)') &
            " p  status          logLik           AIC          AICc           BIC"
      else
         write(*, '(2x,a)') &
            " p  q  status          logLik           AIC          AICc           BIC"
      end if
      do i = 1, size(result%candidates(candidate_index)%order_fits)
         selected_order = order_is_selected( &
            result%candidates(candidate_index), &
            result%candidates(candidate_index)%order_fits(i))
         call display_order_fit( &
            result%candidates(candidate_index)%order_fits(i), &
            model_code == model_arma, selected_order, print_parameters)
      end do
   end subroutine display_order_search

   subroutine display_order_fit(fit, is_arma, selected, print_parameters)
      !! Display diagnostics and optional coefficients for one order fit.
      type(order_search_result_t), intent(in) :: fit !! One tested model order.
      logical, intent(in) :: is_arma !! Include the moving-average order when true.
      logical, intent(in) :: selected !! Mark the order retained for its model class.
      logical, intent(in) :: print_parameters !! Print fitted coefficients when true.
      character(len=8) :: status_text
      integer :: i

      if (fit%converged) then
         status_text = "ok"
         if (selected) status_text = "selected"
         if (is_arma) then
            write(*, '(2x,i2,1x,i2,2x,a8,4(2x,es12.4))') fit%ar_order, &
               fit%ma_order, status_text, fit%log_likelihood, fit%aic, &
               fit%aicc, fit%bic
         else
            write(*, '(2x,i2,4x,a8,4(2x,es12.4))') fit%ar_order, &
               status_text, fit%log_likelihood, fit%aic, fit%aicc, fit%bic
         end if
      else if (is_arma) then
         write(*, '(2x,i2,1x,i2,2x,a,i0,a)') fit%ar_order, fit%ma_order, &
            "failed (status ", fit%info, ")"
      else
         write(*, '(2x,i2,4x,a,i0,a)') fit%ar_order, &
            "failed (status ", fit%info, ")"
      end if
      if (.not. print_parameters .or. .not. fit%converged) return
      write(*, '(9x,a,f16.8)') "innovation variance ", &
         fit%innovation_variance
      write(*, '(9x,a,f16.8)') "mean ", fit%mean
      if (allocated(fit%ar_coefficients)) then
         do i = 1, size(fit%ar_coefficients)
            write(*, '(9x,a,i0,a,f16.8)') "phi(", i, ") ", &
               fit%ar_coefficients(i)
         end do
      end if
      if (allocated(fit%ma_coefficients)) then
         do i = 1, size(fit%ma_coefficients)
            write(*, '(9x,a,i0,a,f16.8)') "theta(", i, ") ", &
               fit%ma_coefficients(i)
         end do
      end if
   end subroutine display_order_fit

   pure logical function order_is_selected(candidate, fit) result(selected)
      !! Test whether an order-search entry was retained for its model class.
      type(candidate_model_result_t), intent(in) :: candidate !! Selected class fit.
      type(order_search_result_t), intent(in) :: fit !! One order-search entry.

      selected = fit%converged .and. candidate%order == fit%ar_order .and. &
         candidate%ma_order == fit%ma_order
   end function order_is_selected

   subroutine display_candidate_summary(result, rank, print_parameters)
      !! Display one ranked candidate and its optional fitted parameters.
      type(automatic_model_result_t), intent(in) :: result !! Automatic modeling result.
      integer, intent(in) :: rank !! Candidate rank to display.
      logical, intent(in) :: print_parameters !! Print fitted parameter values when true.

      write(*, '(2x,i3,2x,a)') rank, result%candidates(rank)%name
      if (result%candidates(rank)%model_code == model_arma) then
         write(*, '(7x,a,i0,a,i0,a)') "order (", &
            result%candidates(rank)%order, ",", &
            result%candidates(rank)%ma_order, ")"
      end if
      if (result%selection_criterion == selection_validation) then
         write(*, '(7x,a,es12.4,2x,a,es12.4)') &
            "RMSE ", result%candidates(rank)%rmse, &
            "MAE ", result%candidates(rank)%mae
      end if
      write(*, '(7x,a,es12.4,2x,a,es12.4,2x,a,es12.4,2x,a,es12.4)') &
         "logLik ", result%candidates(rank)%log_likelihood, &
         "AIC ", result%candidates(rank)%aic, &
         "AICc ", result%candidates(rank)%aicc, &
         "BIC ", result%candidates(rank)%bic
      write(*, '(7x,a,i0,2x,a,i0)') "parameter count ", &
         result%candidates(rank)%parameter_count, "effective observations ", &
         result%candidates(rank)%effective_observations
      if (result%fit_times_recorded) then
         write(*, '(7x,a,f10.4,a)') "full-data fit time ", &
            result%candidates(rank)%full_fit_seconds, " seconds"
         if (result%selection_criterion == selection_validation) then
            write(*, '(7x,a,f10.4,a)') "validation fit time ", &
               result%candidates(rank)%validation_fit_seconds, " seconds"
         end if
      end if
      write(*, '(7x,a)') result%candidates(rank)%rationale
      if (print_parameters) then
         call display_candidate_parameters(result%candidates(rank))
      end if
      write(*, '(a)') ""
   end subroutine display_candidate_summary

   subroutine display_candidate_parameters(candidate)
      !! Display fitted parameters and states available for one candidate.
      type(candidate_model_result_t), intent(in) :: candidate !! Candidate to describe.
      integer :: i

      write(*, '(7x,a,f16.8)') "innovation variance ", &
         candidate%innovation_variance
      select case (candidate%model_code)
      case (model_mean)
         write(*, '(7x,a,f16.8)') "mean ", candidate%level
      case (model_naive)
         write(*, '(7x,a,f16.8)') "last level ", candidate%level
      case (model_drift)
         write(*, '(7x,a,f16.8,2x,a,f16.8)') &
            "last level ", candidate%level, "drift ", candidate%drift
      case (model_ses)
         write(*, '(7x,a,f16.8,2x,a,f16.8)') &
            "alpha ", candidate%alpha, "terminal level ", candidate%level
      case (model_holt)
         write(*, '(7x,a,f16.8,2x,a,f16.8)') &
            "alpha ", candidate%alpha, "beta ", candidate%beta
      case (model_theta)
         write(*, '(7x,a,f16.8)') "alpha ", candidate%alpha
      case (model_holt_winters)
         write(*, '(7x,a,f16.8,2x,a,f16.8,2x,a,f16.8)') &
            "alpha ", candidate%alpha, "beta ", candidate%beta, &
            "gamma ", candidate%gamma
      case (model_ar)
         write(*, '(7x,a,f16.8)') "mean ", candidate%level
         if (allocated(candidate%coefficients)) then
            do i = 1, size(candidate%coefficients)
               write(*, '(7x,a,i0,a,f16.8)') "phi(", i, ") ", &
                  candidate%coefficients(i)
            end do
         end if
      case (model_arma)
         write(*, '(7x,a,f16.8)') "mean ", candidate%level
         if (allocated(candidate%coefficients)) then
            do i = 1, size(candidate%coefficients)
               write(*, '(7x,a,i0,a,f16.8)') "phi(", i, ") ", &
                  candidate%coefficients(i)
            end do
         end if
         if (allocated(candidate%ma_coefficients)) then
            do i = 1, size(candidate%ma_coefficients)
               write(*, '(7x,a,i0,a,f16.8)') "theta(", i, ") ", &
                  candidate%ma_coefficients(i)
            end do
         end if
      end select
   end subroutine display_candidate_parameters

end module automatic_modeling_mod
