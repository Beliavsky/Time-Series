! SPDX-License-Identifier: GPL-3.0-only
! SPDX-FileComment: Algorithms translated from the R narfima package.
module narfima_mod
   !! Neural autoregression with lagged baseline-model innovations.
   use kind_mod, only: dp
   use arfima_mod, only: arfima_model_t, arfima_fit_t, &
      arfima_regression_fit_t, arfima_transfer_t, arfima_model, &
      arfima_fit, arfima_regression_fit, arfima_likelihood_css
   use arima2_mod, only: arima2_fit_t, arima2_fit
   use bsts_mod, only: bsts_mcmc_t, bsts_semilocal_trend, bsts_seasonal
   use neural_network_mod, only: neural_network_t, neural_network_fit, &
      neural_network_predict, neural_network_parameter_count
   use random_mod, only: random_uniform, random_standard_normal
   use stats_mod, only: ols_fit, quantile, sorted, standard_deviation
   use time_series_stats_mod, only: yule_walker_result_t, yule_walker_fit
   use utils_mod, only: quiet_nan
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use, intrinsic :: iso_fortran_env, only: output_unit
   implicit none
   private

   integer, parameter, public :: narfima_baseline_supplied = 0
   integer, parameter, public :: narfima_baseline_arfima = 1
   integer, parameter, public :: narfima_baseline_arima = 2
   integer, parameter, public :: narfima_baseline_bsts = 3
   integer, parameter, public :: narfima_baseline_naive = 4

   type, public :: narfima_scale_t
      !! Center and scale used to standardize one input variable.
      real(dp) :: center = 0.0_dp
      real(dp) :: scale = 1.0_dp
      logical :: active = .false.
   end type narfima_scale_t

   type, public :: narfima_member_t
      !! One neural-network ensemble member and optional direct connection.
      type(neural_network_t) :: network
      real(dp), allocatable :: direct_coefficients(:)
      integer :: info = 0
   end type narfima_member_t

   type, public :: narfima_model_t
      !! Fitted neural ARMA model driven by baseline-model innovations.
      type(narfima_member_t), allocatable :: members(:)
      type(narfima_scale_t) :: response_scale
      type(narfima_scale_t) :: error_scale
      type(narfima_scale_t), allocatable :: regressor_scales(:)
      real(dp), allocatable :: series(:)
      real(dp), allocatable :: errors(:)
      real(dp), allocatable :: regressors(:, :)
      real(dp), allocatable :: fitted(:)
      real(dp), allocatable :: residuals(:)
      integer, allocatable :: fitted_indices(:)
      integer, allocatable :: response_lags(:)
      integer, allocatable :: error_lags(:)
      real(dp) :: response_lambda = 0.5_dp
      real(dp) :: error_lambda = 0.5_dp
      real(dp) :: mse = huge(1.0_dp)
      integer :: period = 1
      integer :: baseline = narfima_baseline_supplied
      integer :: seasonal_order = 0
      integer :: hidden_count = 1
      integer :: repetitions = 0
      integer :: maximum_lag = 0
      integer :: info = 0
      logical :: direct = .true.
      logical :: transform_response = .true.
      logical :: transform_errors = .true.
      logical :: scale_inputs = .true.
   end type narfima_model_t

   type, public :: narfima_forecast_t
      !! Recursive point forecast and optional simulation-path summaries.
      real(dp), allocatable :: mean(:)
      real(dp), allocatable :: paths(:, :)
      real(dp), allocatable :: standard_deviation(:)
      real(dp), allocatable :: levels(:)
      real(dp), allocatable :: lower(:, :)
      real(dp), allocatable :: upper(:, :)
      integer :: info = 0
   end type narfima_forecast_t

   interface display
      module procedure display_narfima_model
      module procedure display_narfima_forecast
   end interface display

   public :: narfima_fit, narfima_forecast, narfima_forecast_from_innovations
   public :: narfima_auto_narfima, narfima_auto_narima
   public :: narfima_auto_nbsts, narfima_auto_nnaive
   public :: display, display_narfima_model, display_narfima_forecast

contains

   function narfima_fit(series, errors, p, q, period, seasonal_order, &
      hidden_count, repetitions, direct, response_lambda, error_lambda, &
      scale_inputs, regressors, max_iterations, tolerance, decay, &
      transform_response, transform_errors) result(model)
      !! Fit the common NARFIMA neural model from supplied baseline innovations.
      real(dp), intent(in) :: series(:) !! Univariate response observations.
      real(dp), intent(in) :: errors(:) !! Aligned baseline-model innovations.
      integer, intent(in), optional :: p !! Number of consecutive response lags.
      integer, intent(in), optional :: q !! Number of consecutive error lags.
      integer, intent(in), optional :: period !! Seasonal period; one is nonseasonal.
      integer, intent(in), optional :: seasonal_order !! Seasonal error-lag count.
      integer, intent(in), optional :: hidden_count !! Hidden units per network.
      integer, intent(in), optional :: repetitions !! Number of network fits.
      logical, intent(in), optional :: direct !! Include direct linear connections.
      real(dp), intent(in), optional :: response_lambda !! Response Box-Cox parameter.
      real(dp), intent(in), optional :: error_lambda !! Innovation Box-Cox parameter.
      logical, intent(in), optional :: scale_inputs !! Standardize transformed inputs.
      real(dp), intent(in), optional :: regressors(:, :) !! Observation-aligned regressors.
      integer, intent(in), optional :: max_iterations !! Maximum network BFGS iterations.
      real(dp), intent(in), optional :: tolerance !! Network gradient tolerance.
      real(dp), intent(in), optional :: decay !! Network weight decay.
      logical, intent(in), optional :: transform_response !! Apply response Box-Cox transformation.
      logical, intent(in), optional :: transform_errors !! Apply innovation Box-Cox transformation.
      type(narfima_model_t) :: model
      real(dp), allocatable :: transformed(:), transformed_errors(:)
      real(dp), allocatable :: scaled_regressors(:, :), predictors(:, :)
      real(dp), allocatable :: response(:), target(:, :), initial(:)
      real(dp), allocatable :: fitted_member(:), fitted_scaled(:)
      logical, allocatable :: valid(:)
      integer :: response_order, error_order, seasons, hidden, reps
      integer :: observations, variables, row, lag, column, member, parameters
      integer :: observation
      integer :: maximum_lag
      logical :: use_direct, use_scaling

      observations = size(series)
      model%period = 1
      if (present(period)) model%period = period
      response_order = select_ar_order(seasonally_adjusted( &
         interpolate_missing(series), model%period))
      if (present(p)) response_order = p
      error_order = 1
      if (present(q)) error_order = q
      seasons = merge(1, 0, model%period > 1)
      if (present(seasonal_order)) seasons = seasonal_order
      hidden = max(1, (response_order + error_order)/2)
      if (present(hidden_count)) hidden = hidden_count
      reps = 100
      if (present(repetitions)) reps = repetitions
      use_direct = .true.
      if (present(direct)) use_direct = direct
      use_scaling = .true.
      if (present(scale_inputs)) use_scaling = scale_inputs
      if (observations < 5 .or. size(errors) /= observations .or. &
         response_order < 1 .or. error_order < 1 .or. model%period < 1 .or. &
         seasons < 0 .or. hidden < 1 .or. reps < 1 .or. &
         count(ieee_is_finite(series)) < 5 .or. &
         count(ieee_is_finite(errors)) < 3) then
         model%info = 1
         return
      end if
      if (present(regressors)) then
         if (size(regressors, 1) /= observations) then
            model%info = 1
            return
         end if
      end if
      model%response_lambda = 0.5_dp
      if (present(response_lambda)) model%response_lambda = response_lambda
      model%error_lambda = 0.5_dp
      if (present(error_lambda)) model%error_lambda = error_lambda
      model%transform_response = .true.
      if (present(transform_response)) model%transform_response = transform_response
      model%transform_errors = .true.
      if (present(transform_errors)) model%transform_errors = transform_errors
      model%scale_inputs = use_scaling
      model%direct = use_direct
      model%seasonal_order = seasons
      model%hidden_count = hidden
      model%repetitions = reps
      model%series = series
      model%errors = errors
      model%response_lags = [(lag, lag=1, response_order)]
      model%error_lags = build_error_lags(error_order, model%period, seasons)
      maximum_lag = max(response_order, maxval(model%error_lags))
      model%maximum_lag = maximum_lag
      if (observations - maximum_lag < 3) then
         model%info = 2
         return
      end if
      if (model%transform_response) then
         transformed = signed_box_cox(series, model%response_lambda)
      else
         transformed = series
      end if
      if (model%transform_errors) then
         transformed_errors = signed_box_cox(errors, model%error_lambda)
      else
         transformed_errors = errors
      end if
      call standardize_vector(transformed, use_scaling, model%response_scale)
      call standardize_vector(transformed_errors, use_scaling, model%error_scale)
      variables = 0
      if (present(regressors)) variables = size(regressors, 2)
      allocate(model%regressor_scales(variables))
      allocate(scaled_regressors(observations, variables))
      if (variables > 0) then
         model%regressors = regressors
         scaled_regressors = regressors
         do column = 1, variables
            call standardize_vector(scaled_regressors(:, column), use_scaling, &
               model%regressor_scales(column))
         end do
      else
         allocate(model%regressors(observations, 0))
      end if
      allocate(valid(observations), source=.false.)
      do observation = maximum_lag + 1, observations
         valid(observation) = ieee_is_finite(transformed(observation))
         if (.not. valid(observation)) cycle
         do lag = 1, response_order
            valid(observation) = valid(observation) .and. &
               ieee_is_finite(transformed(observation - lag))
         end do
         do lag = 1, size(model%error_lags)
            valid(observation) = valid(observation) .and. &
               ieee_is_finite(transformed_errors(observation - &
               model%error_lags(lag)))
         end do
         if (variables > 0) valid(observation) = valid(observation) .and. &
            all(ieee_is_finite(scaled_regressors(observation, :)))
      end do
      model%fitted_indices = pack([(observation, observation=1, observations)], &
         valid)
      if (size(model%fitted_indices) < 3) then
         model%info = 3
         return
      end if
      allocate(response(size(model%fitted_indices)))
      allocate(predictors(size(response), response_order + &
         size(model%error_lags) + variables))
      do row = 1, size(response)
         observation = model%fitted_indices(row)
         response(row) = transformed(observation)
         column = 0
         do lag = 1, response_order
            column = column + 1
            predictors(row, column) = transformed(observation - lag)
         end do
         do lag = 1, size(model%error_lags)
            column = column + 1
            predictors(row, column) = transformed_errors(observation - &
               model%error_lags(lag))
         end do
         if (variables > 0) predictors(row, column + 1:) = &
            scaled_regressors(observation, :)
      end do
      allocate(model%members(reps))
      allocate(fitted_scaled(size(response)), source=0.0_dp)
      target = reshape(response, [size(response), 1])
      do member = 1, reps
         allocate(model%members(member)%direct_coefficients(0))
         parameters = neural_network_parameter_count(size(predictors, 2), &
            hidden, 1, use_direct)
         allocate(initial(parameters))
         do row = 1, parameters
            initial(row) = (2.0_dp*random_uniform() - 1.0_dp)/ &
               sqrt(real(size(predictors, 2), dp))
         end do
         model%members(member)%network = neural_network_fit(predictors, target, &
            hidden, max_iterations, tolerance, decay, initial, use_direct, &
            .false.)
         deallocate(initial)
         if (model%members(member)%network%info /= 0) then
            model%members(member)%info = model%members(member)%network%info
            model%info = 4
            return
         end if
         fitted_member = member_predict(model%members(member), predictors)
         fitted_scaled = fitted_scaled + fitted_member/real(reps, dp)
      end do
      fitted_member = undo_response_transform(model, fitted_scaled)
      allocate(model%fitted(observations), source=quiet_nan())
      allocate(model%residuals(observations), source=quiet_nan())
      model%fitted(model%fitted_indices) = fitted_member
      model%residuals(model%fitted_indices) = &
         series(model%fitted_indices) - fitted_member
      model%mse = sum(model%residuals(model%fitted_indices)**2)/ &
         real(size(model%fitted_indices), dp)
   end function narfima_fit

   function narfima_auto_nnaive(series, p, q, period, seasonal_order, &
      hidden_count, repetitions, direct, response_lambda, error_lambda, &
      scale_inputs, regressors, max_iterations, tolerance, decay, &
      transform_response, transform_errors) result(model)
      !! Fit NNaive using one-step naive innovations as neural error inputs.
      real(dp), intent(in) :: series(:) !! Univariate response observations.
      integer, intent(in), optional :: p !! Number of response lags.
      integer, intent(in), optional :: q !! Number of error lags.
      integer, intent(in), optional :: period !! Seasonal period.
      integer, intent(in), optional :: seasonal_order !! Seasonal error-lag count.
      integer, intent(in), optional :: hidden_count !! Hidden units per network.
      integer, intent(in), optional :: repetitions !! Number of network fits.
      logical, intent(in), optional :: direct !! Include direct linear connections.
      real(dp), intent(in), optional :: response_lambda !! Response Box-Cox parameter.
      real(dp), intent(in), optional :: error_lambda !! Innovation Box-Cox parameter.
      logical, intent(in), optional :: scale_inputs !! Standardize transformed inputs.
      real(dp), intent(in), optional :: regressors(:, :) !! Observation-aligned regressors.
      integer, intent(in), optional :: max_iterations !! Maximum network iterations.
      real(dp), intent(in), optional :: tolerance !! Network gradient tolerance.
      real(dp), intent(in), optional :: decay !! Network weight decay.
      logical, intent(in), optional :: transform_response !! Apply response Box-Cox transformation.
      logical, intent(in), optional :: transform_errors !! Apply innovation Box-Cox transformation.
      type(narfima_model_t) :: model
      real(dp), allocatable :: errors(:), filled(:)

      if (size(series) < 2) then
         model%info = 1
         return
      end if
      filled = interpolate_missing(series)
      allocate(errors(size(series)), source=0.0_dp)
      errors(2:) = filled(2:) - filled(:size(series) - 1)
      model = narfima_fit(series, errors, p, q, period, seasonal_order, &
         hidden_count, repetitions, direct, response_lambda, error_lambda, &
         scale_inputs, regressors, max_iterations, tolerance, decay, &
         transform_response, transform_errors)
      if (model%info == 0) model%baseline = narfima_baseline_naive
   end function narfima_auto_nnaive

   function narfima_auto_narfima(series, initial_model, p, q, period, &
      seasonal_order, hidden_count, repetitions, direct, response_lambda, &
      error_lambda, scale_inputs, regressors, max_iterations, tolerance, &
      decay, transform_response, transform_errors) result(model)
      !! Fit NARFIMA using innovations from the shared ARFIMA estimator.
      real(dp), intent(in) :: series(:) !! Univariate response observations.
      type(arfima_model_t), intent(in), optional :: initial_model !! Baseline ARFIMA start.
      integer, intent(in), optional :: p !! Number of neural response lags.
      integer, intent(in), optional :: q !! Number of neural error lags.
      integer, intent(in), optional :: period !! Seasonal period.
      integer, intent(in), optional :: seasonal_order !! Seasonal error-lag count.
      integer, intent(in), optional :: hidden_count !! Hidden units per network.
      integer, intent(in), optional :: repetitions !! Number of network fits.
      logical, intent(in), optional :: direct !! Include direct linear connections.
      real(dp), intent(in), optional :: response_lambda !! Response Box-Cox parameter.
      real(dp), intent(in), optional :: error_lambda !! Innovation Box-Cox parameter.
      logical, intent(in), optional :: scale_inputs !! Standardize transformed inputs.
      real(dp), intent(in), optional :: regressors(:, :) !! Observation-aligned regressors.
      integer, intent(in), optional :: max_iterations !! Maximum network iterations.
      real(dp), intent(in), optional :: tolerance !! Network gradient tolerance.
      real(dp), intent(in), optional :: decay !! Network weight decay.
      logical, intent(in), optional :: transform_response !! Apply response Box-Cox transformation.
      logical, intent(in), optional :: transform_errors !! Apply innovation Box-Cox transformation.
      type(narfima_model_t) :: model
      type(arfima_model_t) :: start
      type(arfima_fit_t) :: baseline
      type(arfima_regression_fit_t) :: regression_baseline
      type(arfima_transfer_t) :: no_transfers(0)
      real(dp), allocatable :: errors(:)
      real(dp), allocatable :: initial_coefficients(:)
      real(dp), allocatable :: no_transfer_regressors(:, :)
      real(dp), allocatable :: baseline_series(:), baseline_regressors(:, :)
      integer :: first, selected_q, variable

      if (present(initial_model)) then
         start = initial_model
      else
         start = arfima_model([0.1_dp], [0.1_dp], 0.1_dp)
      end if
      baseline_series = interpolate_missing(series)
      allocate(errors(size(series)), source=0.0_dp)
      if (present(regressors)) then
         allocate(baseline_regressors(size(regressors, 1), size(regressors, 2)))
         do variable = 1, size(regressors, 2)
            baseline_regressors(:, variable) = &
               interpolate_missing(regressors(:, variable))
         end do
         allocate(initial_coefficients(size(regressors, 2)), source=0.0_dp)
         allocate(no_transfer_regressors(size(series), 0))
         regression_baseline = arfima_regression_fit(baseline_series, start, &
            baseline_regressors, &
            initial_coefficients, no_transfers, no_transfer_regressors, &
            arfima_likelihood_css, .true.)
         if (regression_baseline%info /= 0 .or. &
            .not. allocated(regression_baseline%likelihood%innovation)) then
            model%info = 10 + regression_baseline%info
            return
         end if
         first = size(series) - &
            size(regression_baseline%likelihood%innovation) + 1
         errors(first:) = regression_baseline%likelihood%innovation
         selected_q = max(1, size(regression_baseline%model%theta))
      else
         baseline = arfima_fit(baseline_series, start, arfima_likelihood_css, &
            .true.)
         if (baseline%info /= 0 .or. &
            .not. allocated(baseline%likelihood%innovation)) then
            model%info = 10 + baseline%info
            return
         end if
         first = size(series) - size(baseline%likelihood%innovation) + 1
         errors(first:) = baseline%likelihood%innovation
         selected_q = max(1, size(baseline%model%theta))
      end if
      if (present(q)) selected_q = q
      model = narfima_fit(series, errors, p, selected_q, period, seasonal_order, &
         hidden_count, repetitions, direct, response_lambda, error_lambda, &
         scale_inputs, regressors, max_iterations, tolerance, decay, &
         transform_response, transform_errors)
      if (model%info == 0) model%baseline = narfima_baseline_arfima
   end function narfima_auto_narfima

   function narfima_auto_narima(series, ar_order, difference_order, ma_order, &
      p, q, period, seasonal_order, hidden_count, repetitions, direct, &
      response_lambda, error_lambda, scale_inputs, regressors, max_iterations, &
      tolerance, decay, transform_response, transform_errors) result(model)
      !! Fit NARIMA using innovations from the shared ARIMA estimator.
      real(dp), intent(in) :: series(:) !! Univariate response observations.
      integer, intent(in), optional :: ar_order !! Baseline AR order.
      integer, intent(in), optional :: difference_order !! Baseline difference order.
      integer, intent(in), optional :: ma_order !! Baseline MA order.
      integer, intent(in), optional :: p !! Number of neural response lags.
      integer, intent(in), optional :: q !! Number of neural error lags.
      integer, intent(in), optional :: period !! Seasonal period.
      integer, intent(in), optional :: seasonal_order !! Seasonal error-lag count.
      integer, intent(in), optional :: hidden_count !! Hidden units per network.
      integer, intent(in), optional :: repetitions !! Number of network fits.
      logical, intent(in), optional :: direct !! Include direct linear connections.
      real(dp), intent(in), optional :: response_lambda !! Response Box-Cox parameter.
      real(dp), intent(in), optional :: error_lambda !! Innovation Box-Cox parameter.
      logical, intent(in), optional :: scale_inputs !! Standardize transformed inputs.
      real(dp), intent(in), optional :: regressors(:, :) !! Observation-aligned regressors.
      integer, intent(in), optional :: max_iterations !! Maximum network iterations.
      real(dp), intent(in), optional :: tolerance !! Network gradient tolerance.
      real(dp), intent(in), optional :: decay !! Network weight decay.
      logical, intent(in), optional :: transform_response !! Apply response Box-Cox transformation.
      logical, intent(in), optional :: transform_errors !! Apply innovation Box-Cox transformation.
      type(narfima_model_t) :: model
      type(arima2_fit_t) :: baseline
      real(dp), allocatable :: errors(:), baseline_series(:)
      real(dp), allocatable :: baseline_regressors(:, :)
      integer :: ar, difference, ma, first, selected_q, variable

      ar = 1
      if (present(ar_order)) ar = ar_order
      difference = 0
      if (present(difference_order)) difference = difference_order
      ma = 1
      if (present(ma_order)) ma = ma_order
      baseline_series = interpolate_missing(series)
      if (present(regressors)) then
         allocate(baseline_regressors(size(regressors, 1), size(regressors, 2)))
         do variable = 1, size(regressors, 2)
            baseline_regressors(:, variable) = &
               interpolate_missing(regressors(:, variable))
         end do
         baseline = arima2_fit(baseline_series, ar, difference, ma, 0, 0, 0, &
            1, max_starts=8, max_repeats=2, max_iterations=200, &
            regressors=baseline_regressors)
      else
         baseline = arima2_fit(baseline_series, ar, difference, ma, 0, 0, 0, &
            1, max_starts=8, max_repeats=2, max_iterations=200)
      end if
      if (baseline%info /= 0 .or. &
         .not. allocated(baseline%fit%likelihood%residuals)) then
         model%info = 20 + baseline%info
         return
      end if
      allocate(errors(size(series)), source=0.0_dp)
      first = size(series) - size(baseline%fit%likelihood%residuals) + 1
      errors(first:) = baseline%fit%likelihood%residuals
      selected_q = max(1, ma)
      if (present(q)) selected_q = q
      model = narfima_fit(series, errors, p, selected_q, period, seasonal_order, &
         hidden_count, repetitions, direct, response_lambda, error_lambda, &
         scale_inputs, regressors, max_iterations, tolerance, decay, &
         transform_response, transform_errors)
      if (model%info == 0) model%baseline = narfima_baseline_arima
   end function narfima_auto_narima

   function narfima_auto_nbsts(series, iterations, burn, p, q, period, &
      seasonal_order, hidden_count, repetitions, direct, response_lambda, &
      error_lambda, scale_inputs, regressors, max_iterations, tolerance, &
      decay, transform_response, transform_errors) result(model)
      !! Fit NBSTS using semilocal-trend and seasonal posterior residuals.
      real(dp), intent(in) :: series(:) !! Univariate response observations.
      integer, intent(in), optional :: iterations !! Baseline MCMC iterations.
      integer, intent(in), optional :: burn !! Baseline MCMC burn-in.
      integer, intent(in), optional :: p !! Number of neural response lags.
      integer, intent(in), optional :: q !! Number of neural error lags.
      integer, intent(in), optional :: period !! Seasonal period.
      integer, intent(in), optional :: seasonal_order !! Seasonal error-lag count.
      integer, intent(in), optional :: hidden_count !! Hidden units per network.
      integer, intent(in), optional :: repetitions !! Number of network fits.
      logical, intent(in), optional :: direct !! Include direct linear connections.
      real(dp), intent(in), optional :: response_lambda !! Response Box-Cox parameter.
      real(dp), intent(in), optional :: error_lambda !! Innovation Box-Cox parameter.
      logical, intent(in), optional :: scale_inputs !! Standardize transformed inputs.
      real(dp), intent(in), optional :: regressors(:, :) !! Observation-aligned regressors.
      integer, intent(in), optional :: max_iterations !! Maximum network iterations.
      real(dp), intent(in), optional :: tolerance !! Network gradient tolerance.
      real(dp), intent(in), optional :: decay !! Network weight decay.
      logical, intent(in), optional :: transform_response !! Apply response Box-Cox transformation.
      logical, intent(in), optional :: transform_errors !! Apply innovation Box-Cox transformation.
      type(narfima_model_t) :: model
      type(bsts_mcmc_t) :: trend_fit, seasonal_fit
      real(dp), allocatable :: errors(:), structural_response(:)
      real(dp), allocatable :: trend_contribution(:), seasonal_contribution(:)
      real(dp), allocatable :: beta(:), standard_errors(:), regression_residuals(:)
      real(dp), allocatable :: baseline_series(:), baseline_regressors(:, :)
      real(dp) :: rss
      integer :: draws, discarded, seasonal_period, observation, cycle, status
      integer :: variable

      draws = 1000
      if (present(iterations)) draws = iterations
      discarded = max(1, draws/10)
      if (present(burn)) discarded = burn
      seasonal_period = 12
      if (present(period)) seasonal_period = period
      if (draws < 2 .or. discarded < 0 .or. discarded >= draws) then
         model%info = 30
         return
      end if
      baseline_series = interpolate_missing(series)
      structural_response = baseline_series
      if (present(regressors)) then
         allocate(baseline_regressors(size(regressors, 1), size(regressors, 2)))
         do variable = 1, size(regressors, 2)
            baseline_regressors(:, variable) = &
               interpolate_missing(regressors(:, variable))
         end do
         call ols_fit(baseline_regressors, baseline_series, beta, standard_errors, &
            regression_residuals, rss, status)
         if (status /= 0) then
            model%info = 31
            return
         end if
         structural_response = regression_residuals
      end if
      allocate(trend_contribution(size(series)), source=0.0_dp)
      allocate(seasonal_contribution(size(series)), source=0.0_dp)
      do cycle = 1, merge(2, 1, seasonal_period > 1)
         trend_fit = bsts_semilocal_trend( &
            structural_response - seasonal_contribution, draws, &
            burn=discarded, force_stationary=.true.)
         if (trend_fit%info /= 0) then
            model%info = 32 + trend_fit%info
            return
         end if
         do observation = 1, size(series)
            trend_contribution(observation) = sum( &
               trend_fit%state(:, observation, discarded + 1:))/ &
               real(size(trend_fit%state, 3) - discarded, dp)
         end do
         if (seasonal_period > 1) then
            seasonal_fit = bsts_seasonal( &
               structural_response - trend_contribution, seasonal_period, &
               draws, burn=discarded)
            if (seasonal_fit%info /= 0) then
               model%info = 40 + seasonal_fit%info
               return
            end if
            do observation = 1, size(series)
               seasonal_contribution(observation) = sum( &
                  seasonal_fit%state(:, observation, discarded + 1:))/ &
                  real(size(seasonal_fit%state, 3) - discarded, dp)
            end do
         end if
      end do
      errors = structural_response - trend_contribution - seasonal_contribution
      model = narfima_fit(series, errors, p, q, period, seasonal_order, &
         hidden_count, repetitions, direct, response_lambda, error_lambda, &
         scale_inputs, regressors, max_iterations, tolerance, decay, &
         transform_response, transform_errors)
      if (model%info == 0) model%baseline = narfima_baseline_bsts
   end function narfima_auto_nbsts

   pure function narfima_forecast_from_innovations(model, innovations, &
      future_regressors, levels) result(out)
      !! Recursively forecast from caller-supplied original-scale innovations.
      type(narfima_model_t), intent(in) :: model !! Fitted NARFIMA model.
      real(dp), intent(in) :: innovations(:, :) !! Horizon-by-path innovations.
      real(dp), intent(in), optional :: future_regressors(:, :) !! Future regressors.
      real(dp), intent(in), optional :: levels(:) !! Central interval percentages.
      type(narfima_forecast_t) :: out
      real(dp), allocatable :: response_history(:), error_history(:)
      real(dp), allocatable :: predictor(:, :), scaled_regressors(:, :)
      real(dp), allocatable :: member_values(:), prediction_values(:), ordered(:)
      real(dp) :: scaled_prediction, scaled_error, scaled_innovation, probability
      integer :: horizon, paths, variables, path, step, lag, column, member, level

      horizon = size(innovations, 1)
      paths = size(innovations, 2)
      variables = size(model%regressor_scales)
      if (model%info /= 0 .or. horizon < 1 .or. paths < 1 .or. &
         .not. all(ieee_is_finite(innovations))) then
         out%info = 1
         return
      end if
      if (variables > 0) then
         if (.not. present(future_regressors)) then
            out%info = 2
            return
         end if
         if (size(future_regressors, 1) /= horizon .or. &
            size(future_regressors, 2) /= variables .or. &
            .not. all(ieee_is_finite(future_regressors))) then
            out%info = 2
            return
         end if
         scaled_regressors = future_regressors
         do column = 1, variables
            scaled_regressors(:, column) = apply_scale( &
               scaled_regressors(:, column), model%regressor_scales(column))
         end do
      else
         allocate(scaled_regressors(horizon, 0))
      end if
      allocate(out%paths(horizon, paths))
      allocate(predictor(1, size(model%response_lags) + &
         size(model%error_lags) + variables))
      allocate(member_values(model%repetitions))
      do path = 1, paths
         if (model%transform_response) then
            response_history = signed_box_cox(model%series, &
               model%response_lambda)
         else
            response_history = model%series
         end if
         response_history = apply_scale(response_history, model%response_scale)
         if (model%transform_errors) then
            error_history = signed_box_cox(model%errors, model%error_lambda)
         else
            error_history = model%errors
         end if
         error_history = apply_scale(error_history, model%error_scale)
         do step = 1, horizon
            column = 0
            do lag = 1, size(model%response_lags)
               column = column + 1
               predictor(1, column) = response_history( &
                  size(response_history) + 1 - model%response_lags(lag))
            end do
            do lag = 1, size(model%error_lags)
               column = column + 1
               predictor(1, column) = error_history( &
                  size(error_history) + 1 - model%error_lags(lag))
            end do
            if (variables > 0) predictor(1, column + 1:) = &
               scaled_regressors(step, :)
            if (.not. all(ieee_is_finite(predictor))) then
               out%info = 4
               return
            end if
            do member = 1, model%repetitions
               prediction_values = member_predict(model%members(member), &
                  predictor)
               member_values(member) = prediction_values(1)
            end do
            scaled_innovation = innovations(step, path)
            if (model%response_scale%active) scaled_innovation = &
               scaled_innovation/model%response_scale%scale
            scaled_error = innovations(step, path)
            if (model%error_scale%active) scaled_error = &
               scaled_error/model%error_scale%scale
            scaled_prediction = sum(member_values)/real(model%repetitions, dp) + &
               scaled_innovation
            response_history = [response_history, scaled_prediction]
            error_history = [error_history, scaled_error]
            out%paths(step, path) = undo_response_transform_scalar(model, &
               scaled_prediction)
         end do
      end do
      allocate(out%mean(horizon), out%standard_deviation(horizon))
      out%mean = sum(out%paths, 2)/real(paths, dp)
      do step = 1, horizon
         out%standard_deviation(step) = standard_deviation(out%paths(step, :))
      end do
      if (present(levels)) then
         if (any(levels <= 0.0_dp) .or. any(levels >= 100.0_dp)) then
            out%info = 3
            return
         end if
         out%levels = levels
         allocate(out%lower(horizon, size(levels)))
         allocate(out%upper(horizon, size(levels)))
         do step = 1, horizon
            ordered = sorted(out%paths(step, :))
            do level = 1, size(levels)
               probability = 0.5_dp - levels(level)/200.0_dp
               out%lower(step, level) = quantile_type8(ordered, probability)
               out%upper(step, level) = quantile_type8(ordered, 1.0_dp - probability)
            end do
         end do
      else
         allocate(out%levels(0), out%lower(horizon, 0), out%upper(horizon, 0))
      end if
   end function narfima_forecast_from_innovations

   function narfima_forecast(model, horizon, simulations, future_regressors, &
      levels, bootstrap) result(out)
      !! Forecast with shared-RNG Gaussian or residual-bootstrap innovations.
      type(narfima_model_t), intent(in) :: model !! Fitted NARFIMA model.
      integer, intent(in) :: horizon !! Positive forecast horizon.
      integer, intent(in), optional :: simulations !! Number of forecast paths.
      real(dp), intent(in), optional :: future_regressors(:, :) !! Future regressors.
      real(dp), intent(in), optional :: levels(:) !! Central interval percentages.
      logical, intent(in), optional :: bootstrap !! Resample centered fit residuals.
      type(narfima_forecast_t) :: out
      real(dp), allocatable :: innovations(:, :), usable(:)
      real(dp) :: scale, draw
      integer :: paths, path, step, index
      logical :: use_bootstrap

      paths = 1
      if (present(simulations)) paths = simulations
      use_bootstrap = .false.
      if (present(bootstrap)) use_bootstrap = bootstrap
      if (model%info /= 0 .or. horizon < 1 .or. paths < 1) then
         out%info = 1
         return
      end if
      usable = pack(model%residuals, ieee_is_finite(model%residuals))
      usable = usable - sum(usable)/real(size(usable), dp)
      scale = standard_deviation(usable)
      allocate(innovations(horizon, paths))
      do path = 1, paths
         do step = 1, horizon
            if (use_bootstrap) then
               draw = random_uniform()
               index = min(size(usable), 1 + int(draw*real(size(usable), dp)))
               innovations(step, path) = usable(index)
            else if (paths == 1) then
               innovations(step, path) = 0.0_dp
            else
               innovations(step, path) = scale*random_standard_normal()
            end if
         end do
      end do
      out = narfima_forecast_from_innovations(model, innovations, &
         future_regressors, levels)
   end function narfima_forecast

   subroutine display_narfima_model(model, unit, print_obs)
      !! Display a fitted NARFIMA model and optionally its observations.
      type(narfima_model_t), intent(in) :: model !! Model to display.
      integer, intent(in), optional :: unit !! Output unit; defaults to standard output.
      logical, intent(in), optional :: print_obs !! Whether to print fitted observations.
      integer :: destination, observation
      logical :: show_observations

      destination = output_unit
      if (present(unit)) destination = unit
      show_observations = .false.
      if (present(print_obs)) show_observations = print_obs
      write(destination, '(a)') 'Neural autoregressive innovation model'
      write(destination, '(a, i0)') 'Status: ', model%info
      write(destination, '(a, i0)') 'Baseline type: ', model%baseline
      write(destination, '(a, *(i0, 1x))') 'Response lags: ', model%response_lags
      write(destination, '(a, *(i0, 1x))') 'Error lags: ', model%error_lags
      write(destination, '(a, i0)') 'Hidden nodes: ', model%hidden_count
      write(destination, '(a, i0)') 'Ensemble members: ', model%repetitions
      write(destination, '(a, l1)') 'Direct connections: ', model%direct
      write(destination, '(a, es14.6)') 'MSE: ', model%mse
      if (show_observations .and. allocated(model%fitted)) then
         write(destination, '(a)') 'Index, fitted, residual:'
         do observation = model%maximum_lag + 1, size(model%fitted)
            write(destination, '(i8, 2(1x, es14.6))') observation, &
               model%fitted(observation), model%residuals(observation)
         end do
      end if
   end subroutine display_narfima_model

   subroutine display_narfima_forecast(forecast, unit)
      !! Display NARFIMA forecast means and path standard deviations.
      type(narfima_forecast_t), intent(in) :: forecast !! Forecast to display.
      integer, intent(in), optional :: unit !! Output unit; defaults to standard output.
      integer :: destination, step

      destination = output_unit
      if (present(unit)) destination = unit
      write(destination, '(a)') 'NARFIMA forecast'
      write(destination, '(a, i0)') 'Status: ', forecast%info
      if (allocated(forecast%mean)) then
         write(destination, '(a)') 'Horizon, mean, standard deviation:'
         do step = 1, size(forecast%mean)
            write(destination, '(i8, 2(1x, es14.6))') step, &
               forecast%mean(step), forecast%standard_deviation(step)
         end do
      end if
   end subroutine display_narfima_forecast

   pure function build_error_lags(q, period, seasonal_order) result(lags)
      !! Combine ordinary and seasonal innovation lags without duplicates.
      integer, intent(in) :: q !! Consecutive ordinary lag count.
      integer, intent(in) :: period !! Seasonal period.
      integer, intent(in) :: seasonal_order !! Seasonal lag count.
      integer, allocatable :: lags(:)
      logical, allocatable :: included(:)
      integer :: maximum, lag, count_lags

      maximum = max(q, period*seasonal_order)
      allocate(included(maximum), source=.false.)
      included(:q) = .true.
      if (period > 1) then
         do lag = 1, seasonal_order
            included(period*lag) = .true.
         end do
      end if
      count_lags = count(included)
      allocate(lags(count_lags))
      count_lags = 0
      do lag = 1, maximum
         if (.not. included(lag)) cycle
         count_lags = count_lags + 1
         lags(count_lags) = lag
      end do
   end function build_error_lags

   pure integer function select_ar_order(series) result(order)
      !! Select the response-lag order by Yule-Walker AIC.
      real(dp), intent(in) :: series(:) !! Finite response observations.
      type(yule_walker_result_t) :: candidate
      real(dp) :: best_criterion
      integer :: trial, maximum

      order = 1
      best_criterion = huge(1.0_dp)
      maximum = max(1, min(size(series) - 1, &
         int(floor(10.0_dp*log10(real(size(series), dp))))))
      do trial = 1, maximum
         candidate = yule_walker_fit(series, trial)
         if (candidate%info == 0 .and. &
            candidate%criterion < best_criterion) then
            best_criterion = candidate%criterion
            order = trial
         end if
      end do
   end function select_ar_order

   pure function seasonally_adjusted(series, period) result(adjusted)
      !! Remove mean seasonal positions before automatic AR-order selection.
      real(dp), intent(in) :: series(:) !! Response observations.
      integer, intent(in) :: period !! Seasonal period.
      real(dp) :: adjusted(size(series))
      real(dp), allocatable :: seasonal_mean(:)
      real(dp) :: overall
      integer :: season, count_seasons

      adjusted = series
      if (period <= 1 .or. size(series) <= 2*period) return
      allocate(seasonal_mean(period))
      overall = sum(series)/real(size(series), dp)
      do season = 1, period
         count_seasons = 1 + (size(series) - season)/period
         seasonal_mean(season) = sum(series(season::period))/ &
            real(count_seasons, dp) - overall
         adjusted(season::period) = series(season::period) - &
            seasonal_mean(season)
      end do
   end function seasonally_adjusted

   pure function interpolate_missing(values) result(filled)
      !! Linearly interpolate internal nonfinite values and extend endpoints.
      real(dp), intent(in) :: values(:) !! Possibly incomplete observations.
      real(dp) :: filled(size(values))
      integer :: first, last, left, right, index

      filled = values
      first = 0
      do index = 1, size(values)
         if (ieee_is_finite(values(index))) then
            first = index
            exit
         end if
      end do
      if (first == 0) return
      last = first
      do index = size(values), first, -1
         if (ieee_is_finite(values(index))) then
            last = index
            exit
         end if
      end do
      filled(:first - 1) = values(first)
      filled(last + 1:) = values(last)
      left = first
      do while (left < last)
         right = left + 1
         do while (right <= last)
            if (ieee_is_finite(values(right))) exit
            right = right + 1
         end do
         do index = left + 1, right - 1
            filled(index) = values(left) + &
               real(index - left, dp)/real(right - left, dp)* &
               (values(right) - values(left))
         end do
         left = right
      end do
   end function interpolate_missing

   pure elemental real(dp) function signed_box_cox_scalar(value, lambda) &
      result(transformed)
      !! Apply forecast's signed Box-Cox transformation.
      real(dp), intent(in) :: value !! Value to transform.
      real(dp), intent(in) :: lambda !! Transformation parameter.

      if (abs(lambda) < epsilon(1.0_dp)) then
         transformed = log(value)
      else
         transformed = (sign(1.0_dp, value)*abs(value)**lambda - 1.0_dp)/lambda
      end if
   end function signed_box_cox_scalar

   pure function signed_box_cox(values, lambda) result(transformed)
      !! Apply the signed Box-Cox transformation to a vector.
      real(dp), intent(in) :: values(:) !! Values to transform.
      real(dp), intent(in) :: lambda !! Transformation parameter.
      real(dp) :: transformed(size(values))

      transformed = signed_box_cox_scalar(values, lambda)
   end function signed_box_cox

   pure elemental real(dp) function signed_inverse_box_cox(value, lambda) &
      result(transformed)
      !! Invert forecast's signed Box-Cox transformation.
      real(dp), intent(in) :: value !! Transformed value.
      real(dp), intent(in) :: lambda !! Transformation parameter.
      real(dp) :: shifted

      if (abs(lambda) < epsilon(1.0_dp)) then
         transformed = exp(value)
      else
         shifted = lambda*value + 1.0_dp
         transformed = sign(1.0_dp, shifted)*abs(shifted)**(1.0_dp/lambda)
      end if
   end function signed_inverse_box_cox

   pure subroutine standardize_vector(values, active, scale)
      !! Standardize a vector in place and retain its transformation.
      real(dp), intent(inout) :: values(:) !! Values replaced by standardized values.
      logical, intent(in) :: active !! Whether scaling is requested.
      type(narfima_scale_t), intent(out) :: scale !! Estimated scaling parameters.
      logical :: finite(size(values))
      integer :: observations

      scale%active = active
      if (.not. active) return
      finite = ieee_is_finite(values)
      observations = count(finite)
      if (observations < 1) then
         scale%active = .false.
         return
      end if
      scale%center = sum(values, mask=finite)/real(observations, dp)
      if (observations > 1) then
         scale%scale = sqrt(sum((values - scale%center)**2, mask=finite)/ &
            real(observations - 1, dp))
      end if
      scale%scale = max(scale%scale, sqrt(epsilon(1.0_dp)))
      where (finite) values = (values - scale%center)/scale%scale
   end subroutine standardize_vector

   pure function apply_scale(values, scale) result(transformed)
      !! Apply retained standardization to a vector.
      real(dp), intent(in) :: values(:) !! Values to standardize.
      type(narfima_scale_t), intent(in) :: scale !! Scaling parameters.
      real(dp) :: transformed(size(values))

      transformed = values
      if (scale%active) transformed = (values - scale%center)/scale%scale
   end function apply_scale

   pure real(dp) function apply_scale_scalar(value, scale) result(transformed)
      !! Apply retained standardization to one value.
      real(dp), intent(in) :: value !! Value to standardize.
      type(narfima_scale_t), intent(in) :: scale !! Scaling parameters.

      transformed = value
      if (scale%active) transformed = (value - scale%center)/scale%scale
   end function apply_scale_scalar

   pure function member_predict(member, predictors) result(prediction)
      !! Predict with one nonlinear member plus its direct connection.
      type(narfima_member_t), intent(in) :: member !! Fitted ensemble member.
      real(dp), intent(in) :: predictors(:, :) !! Predictor matrix.
      real(dp), allocatable :: prediction(:)
      real(dp), allocatable :: nonlinear(:, :), design(:, :)

      nonlinear = neural_network_predict(member%network, predictors)
      prediction = nonlinear(:, 1)
      if (size(member%direct_coefficients) > 0) then
         allocate(design(size(predictors, 1), size(predictors, 2) + 1))
         design(:, 1) = 1.0_dp
         design(:, 2:) = predictors
         prediction = prediction + matmul(design, member%direct_coefficients)
      end if
   end function member_predict

   pure function undo_response_transform(model, scaled) result(values)
      !! Undo response standardization and signed Box-Cox transformation.
      type(narfima_model_t), intent(in) :: model !! Fitted model.
      real(dp), intent(in) :: scaled(:) !! Scaled transformed values.
      real(dp) :: values(size(scaled))

      values = scaled
      if (model%response_scale%active) values = values* &
         model%response_scale%scale + model%response_scale%center
      if (model%transform_response) &
         values = signed_inverse_box_cox(values, model%response_lambda)
   end function undo_response_transform

   pure real(dp) function undo_response_transform_scalar(model, scaled) &
      result(value)
      !! Undo the fitted response transformations for one value.
      type(narfima_model_t), intent(in) :: model !! Fitted model.
      real(dp), intent(in) :: scaled !! Scaled transformed value.
      real(dp) :: transformed

      transformed = scaled
      if (model%response_scale%active) transformed = transformed* &
         model%response_scale%scale + model%response_scale%center
      if (model%transform_response) then
         value = signed_inverse_box_cox(transformed, model%response_lambda)
      else
         value = transformed
      end if
   end function undo_response_transform_scalar

   pure real(dp) function quantile_type8(ordered, probability) result(value)
      !! Interpolate an ordered sample using R quantile type 8.
      real(dp), intent(in) :: ordered(:) !! Ascending sample.
      real(dp), intent(in) :: probability !! Probability in [0,1].
      real(dp) :: position, fraction
      integer :: lower

      if (size(ordered) == 0) then
         value = 0.0_dp
         return
      end if
      position = (real(size(ordered), dp) + 1.0_dp/3.0_dp)* &
         max(0.0_dp, min(1.0_dp, probability)) + 1.0_dp/3.0_dp
      if (position <= 1.0_dp) then
         value = ordered(1)
      else if (position >= real(size(ordered), dp)) then
         value = ordered(size(ordered))
      else
         lower = int(floor(position))
         fraction = position - real(lower, dp)
         value = (1.0_dp - fraction)*ordered(lower) + &
            fraction*ordered(lower + 1)
      end if
   end function quantile_type8

end module narfima_mod
