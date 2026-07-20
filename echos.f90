! SPDX-License-Identifier: GPL-3.0-only
! SPDX-FileComment: Algorithms translated from the R echos package.
module echos_mod
   !! Echo State Networks for automatic univariate modeling and forecasting.
   use kind_mod, only: dp
   use linalg_mod, only: general_eigenvalues, invert_matrix, solve_matrix
   use random_mod, only: random_uniform, set_random_seed
   use stats_mod, only: quantile, sort, standard_deviation
   use utils_mod, only: quiet_nan
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use, intrinsic :: iso_fortran_env, only: output_unit
   implicit none
   private

   integer, parameter, public :: echos_ic_aic = 1
   integer, parameter, public :: echos_ic_aicc = 2
   integer, parameter, public :: echos_ic_bic = 3
   integer, parameter, public :: echos_ic_hqc = 4

   type, public :: echos_kpss_t
      !! KPSS stationarity statistic, critical value, and decision.
      real(dp) :: statistic = 0.0_dp
      real(dp) :: critical_value = 0.0_dp
      real(dp) :: alpha = 0.05_dp
      logical :: trend = .false.
      logical :: reject = .false.
      integer :: info = 0
   end type echos_kpss_t

   type, public :: echos_ridge_t
      !! Ridge readout coefficients, fitted values, and selection metrics.
      real(dp), allocatable :: coefficients(:)
      real(dp), allocatable :: fitted(:)
      real(dp), allocatable :: residuals(:)
      real(dp) :: log_likelihood = -huge(1.0_dp)
      real(dp) :: degrees_of_freedom = 0.0_dp
      real(dp) :: lambda = 0.0_dp
      real(dp) :: aic = huge(1.0_dp)
      real(dp) :: aicc = huge(1.0_dp)
      real(dp) :: bic = huge(1.0_dp)
      real(dp) :: hqc = huge(1.0_dp)
      real(dp) :: mse = huge(1.0_dp)
      real(dp) :: mae = huge(1.0_dp)
      integer :: info = 0
   end type echos_ridge_t

   type, public :: echos_model_t
      !! Trained ESN reservoir, ridge readout, preprocessing, and diagnostics.
      real(dp), allocatable :: input_weights(:, :)
      real(dp), allocatable :: reservoir_weights(:, :)
      real(dp), allocatable :: output_weights(:)
      real(dp), allocatable :: states(:, :)
      real(dp), allocatable :: series(:)
      real(dp), allocatable :: transformed(:)
      real(dp), allocatable :: actual(:)
      real(dp), allocatable :: fitted(:)
      real(dp), allocatable :: residuals(:)
      real(dp), allocatable :: scaled_residuals(:)
      real(dp), allocatable :: candidate_lambda(:)
      real(dp), allocatable :: candidate_criterion(:)
      real(dp), allocatable :: candidate_log_likelihood(:)
      real(dp), allocatable :: candidate_degrees_of_freedom(:)
      real(dp), allocatable :: candidate_aic(:)
      real(dp), allocatable :: candidate_aicc(:)
      real(dp), allocatable :: candidate_bic(:)
      real(dp), allocatable :: candidate_hqc(:)
      real(dp), allocatable :: candidate_mse(:)
      real(dp), allocatable :: candidate_mae(:)
      integer, allocatable :: lags(:)
      real(dp) :: original_minimum = 0.0_dp
      real(dp) :: original_maximum = 1.0_dp
      real(dp) :: scaled_minimum = -0.5_dp
      real(dp) :: scaled_maximum = 0.5_dp
      real(dp) :: alpha = 1.0_dp
      real(dp) :: rho = 1.0_dp
      real(dp) :: density = 0.5_dp
      real(dp) :: lambda = 0.0_dp
      real(dp) :: degrees_of_freedom = 0.0_dp
      integer :: difference_order = 0
      integer :: state_count = 0
      integer :: model_count = 0
      integer :: initial_count = 0
      integer :: information_criterion = echos_ic_bic
      integer :: info = 0
   end type echos_model_t

   type, public :: echos_forecast_t
      !! Recursive point forecasts and moving-block-bootstrap uncertainty.
      real(dp), allocatable :: point(:)
      real(dp), allocatable :: interval(:, :)
      real(dp), allocatable :: simulation(:, :)
      real(dp), allocatable :: standard_deviation(:)
      real(dp), allocatable :: levels(:)
      real(dp), allocatable :: actual(:)
      real(dp), allocatable :: fitted(:)
      integer :: horizon = 0
      integer :: simulations = 0
      integer :: info = 0
   end type echos_forecast_t

   type, public :: echos_tuning_t
      !! Expanding-window ESN hyperparameter grid results.
      real(dp), allocatable :: alpha(:)
      real(dp), allocatable :: rho(:)
      real(dp), allocatable :: tau(:)
      real(dp), allocatable :: mse(:)
      real(dp), allocatable :: mae(:)
      real(dp), allocatable :: forecasts(:, :)
      integer, allocatable :: split(:)
      integer, allocatable :: train_end(:)
      integer, allocatable :: test_start(:)
      integer, allocatable :: test_end(:)
      integer, allocatable :: configuration(:)
      real(dp), allocatable :: actual(:)
      integer :: best = 0
      integer :: horizon = 0
      integer :: info = 0
   end type echos_tuning_t

   interface display
      module procedure display_echos_model
      module procedure display_echos_forecast
      module procedure display_echos_tuning
   end interface display

   public :: echos_run_reservoir, echos_fit_ridge
   public :: echos_kpss, echos_estimate_differences
   public :: echos_fit, echos_forecast, echos_tune, display

contains

   pure function echos_run_reservoir(inputs, input_weights, &
      reservoir_weights, alpha, initial_state) result(states)
      !! Generate leaky-tanh reservoir states from fixed weights and inputs.
      real(dp), intent(in) :: inputs(:, :) !! Observation-by-input matrix.
      real(dp), intent(in) :: input_weights(:, :) !! State-by-input weights.
      real(dp), intent(in) :: reservoir_weights(:, :) !! Recurrent state weights.
      real(dp), intent(in) :: alpha !! Leakage rate in (0, 1].
      real(dp), intent(in), optional :: initial_state(:) !! State preceding the first input row.
      real(dp), allocatable :: states(:, :)
      real(dp), allocatable :: previous(:), update(:)
      integer :: observation, state_count, start_observation

      state_count = size(reservoir_weights, 1)
      if (size(inputs, 2) /= size(input_weights, 2) .or. &
         size(input_weights, 1) /= state_count .or. &
         size(reservoir_weights, 2) /= state_count .or. &
         alpha <= 0.0_dp .or. alpha > 1.0_dp) then
         allocate(states(0, 0))
         return
      end if
      allocate(states(size(inputs, 1), state_count), source=0.0_dp)
      allocate(previous(state_count), source=0.0_dp)
      if (present(initial_state)) then
         if (size(initial_state) /= state_count) then
            deallocate(states)
            allocate(states(0, 0))
            return
         end if
         previous = initial_state
      end if
      start_observation = 1
      if (.not. present(initial_state) .and. size(inputs, 1) > 0) then
         states(1, :) = previous
         start_observation = 2
      end if
      do observation = start_observation, size(inputs, 1)
         update = tanh(matmul(input_weights, inputs(observation, :)) + &
            matmul(reservoir_weights, previous))
         states(observation, :) = alpha*update + (1.0_dp - alpha)*previous
         previous = states(observation, :)
      end do
   end function echos_run_reservoir

   pure function echos_fit_ridge(design, response, lambda) result(fit)
      !! Fit an intercept-unpenalized ridge readout and calculate echos metrics.
      real(dp), intent(in) :: design(:, :) !! Design matrix with intercept first.
      real(dp), intent(in) :: response(:) !! Response observations.
      real(dp), intent(in) :: lambda !! Nonnegative ridge penalty.
      type(echos_ridge_t) :: fit
      real(dp), allocatable :: cross_product(:, :), right(:, :), solution(:, :)
      real(dp), allocatable :: inverse(:, :), penalized_df(:, :)
      real(dp) :: rss, observations
      integer :: variable, status

      if (size(design, 1) /= size(response) .or. size(response) < 2 .or. &
         size(design, 2) < 1 .or. lambda < 0.0_dp .or. &
         .not. all(ieee_is_finite(design)) .or. &
         .not. all(ieee_is_finite(response))) then
         fit%info = 1
         return
      end if
      cross_product = matmul(transpose(design), design)
      penalized_df = cross_product
      do variable = 1, size(design, 2)
         penalized_df(variable, variable) = penalized_df(variable, variable) + &
            lambda
      end do
      cross_product = penalized_df
      cross_product(1, 1) = cross_product(1, 1) - lambda
      allocate(right(size(design, 2), 1))
      right(:, 1) = matmul(transpose(design), response)
      allocate(solution(size(design, 2), 1))
      call solve_matrix(cross_product, right, solution, status)
      if (status /= 0) then
         fit%info = 2
         return
      end if
      fit%coefficients = solution(:, 1)
      fit%fitted = matmul(design, fit%coefficients)
      fit%residuals = response - fit%fitted
      fit%lambda = lambda
      call invert_matrix(penalized_df, inverse, status)
      if (status /= 0) then
         fit%info = 3
         return
      end if
      fit%degrees_of_freedom = trace_matrix(matmul(inverse, &
         matmul(transpose(design), design)))
      observations = real(size(response), dp)
      rss = max(sum(fit%residuals**2), tiny(1.0_dp))
      fit%log_likelihood = -0.5_dp*observations*(log(2.0_dp*acos(-1.0_dp)) + &
         1.0_dp + log(rss/observations))
      fit%aic = -2.0_dp*fit%log_likelihood + &
         2.0_dp*fit%degrees_of_freedom
      if (observations > fit%degrees_of_freedom + 1.0_dp) then
         fit%aicc = fit%aic + &
            (2.0_dp*fit%degrees_of_freedom**2 + &
            2.0_dp*fit%degrees_of_freedom)/ &
            (observations - fit%degrees_of_freedom - 1.0_dp)
      end if
      fit%bic = -2.0_dp*fit%log_likelihood + &
         log(observations)*fit%degrees_of_freedom
      fit%hqc = -2.0_dp*fit%log_likelihood + &
         2.0_dp*fit%degrees_of_freedom*log(log(observations))
      fit%mse = rss/observations
      fit%mae = sum(abs(fit%residuals))/observations
   end function echos_fit_ridge

   pure function echos_kpss(series, trend, alpha) result(test)
      !! Perform the echos KPSS stationarity test with Bartlett covariance.
      real(dp), intent(in) :: series(:) !! Finite univariate observations.
      logical, intent(in), optional :: trend !! Include a deterministic linear trend.
      real(dp), intent(in), optional :: alpha !! Significance level among package values.
      type(echos_kpss_t) :: test
      real(dp), allocatable :: residuals(:), cumulative(:)
      real(dp) :: mean_time, mean_series, slope, denominator, long_variance
      real(dp) :: covariance, weight
      integer :: observations, bandwidth, lag

      observations = size(series)
      test%trend = .false.
      if (present(trend)) test%trend = trend
      if (present(alpha)) test%alpha = alpha
      if (observations < 3 .or. .not. all(ieee_is_finite(series)) .or. &
         .not. valid_kpss_alpha(test%alpha)) then
         test%info = 1
         return
      end if
      mean_series = sum(series)/real(observations, dp)
      allocate(residuals(observations))
      if (test%trend) then
         mean_time = 0.5_dp*real(observations + 1, dp)
         denominator = 0.0_dp
         slope = 0.0_dp
         do lag = 1, observations
            denominator = denominator + (real(lag, dp) - mean_time)**2
            slope = slope + (real(lag, dp) - mean_time)* &
               (series(lag) - mean_series)
         end do
         slope = slope/denominator
         do lag = 1, observations
            residuals(lag) = series(lag) - mean_series - &
               slope*(real(lag, dp) - mean_time)
         end do
      else
         residuals = series - mean_series
      end if
      bandwidth = max(1, floor(3.0_dp*sqrt(real(observations, dp))/13.0_dp))
      long_variance = sum(residuals**2)/real(observations, dp)
      do lag = 1, bandwidth
         covariance = sum(residuals(lag + 1:)* &
            residuals(:observations - lag))/real(observations, dp)
         weight = 1.0_dp - real(lag, dp)/real(bandwidth + 1, dp)
         long_variance = long_variance + 2.0_dp*weight*covariance
      end do
      if (long_variance <= tiny(1.0_dp)) then
         test%info = 2
         return
      end if
      allocate(cumulative(observations))
      cumulative(1) = residuals(1)
      do lag = 2, observations
         cumulative(lag) = cumulative(lag - 1) + residuals(lag)
      end do
      test%statistic = sum(cumulative**2)/ &
         (real(observations, dp)**2*long_variance)
      test%critical_value = kpss_critical_value(test%trend, test%alpha)
      test%reject = test%statistic > test%critical_value
   end function echos_kpss

   pure function echos_estimate_differences(series, maximum, trend, alpha) &
      result(order)
      !! Select an ordinary difference order by repeated echos KPSS tests.
      real(dp), intent(in) :: series(:) !! Finite univariate observations.
      integer, intent(in), optional :: maximum !! Maximum order checked; defaults to one.
      logical, intent(in), optional :: trend !! Include a deterministic trend in KPSS tests.
      real(dp), intent(in), optional :: alpha !! KPSS significance level.
      integer :: order
      type(echos_kpss_t) :: test
      real(dp), allocatable :: transformed(:)
      integer :: limit, candidate

      limit = 1
      if (present(maximum)) limit = maximum
      if (limit < 0 .or. size(series) <= limit + 2) then
         order = -1
         return
      end if
      do candidate = 0, limit
         transformed = difference_series(series, candidate)
         test = echos_kpss(transformed, trend, alpha)
         if (test%info /= 0) then
            order = -1
            return
         end if
         if (.not. test%reject) then
            order = candidate
            return
         end if
      end do
      order = limit
   end function echos_estimate_differences

   function echos_fit(series, lags, information_criterion, difference_order, &
      model_count, state_count, initial_count, seed, alpha, rho, tau, density, &
      lambda_range, input_weight_scale, reservoir_weight_scale, scale_range) &
      result(model)
      !! Train an automatic echos ESN with random reservoir and ridge selection.
      real(dp), intent(in) :: series(:) !! Finite univariate observations.
      integer, intent(in), optional :: lags(:) !! Positive autoregressive input lags.
      integer, intent(in), optional :: information_criterion !! AIC, AICc, BIC, or HQC code.
      integer, intent(in), optional :: difference_order !! Ordinary differences; automatic when absent.
      integer, intent(in), optional :: model_count !! Number of random ridge penalties.
      integer, intent(in), optional :: state_count !! Number of reservoir states.
      integer, intent(in), optional :: initial_count !! Initial reservoir states discarded.
      integer, intent(in), optional :: seed !! Shared random-number seed.
      real(dp), intent(in), optional :: alpha !! Leakage rate in (0, 1].
      real(dp), intent(in), optional :: rho !! Target reservoir spectral radius.
      real(dp), intent(in), optional :: tau !! Default reservoir-size fraction.
      real(dp), intent(in), optional :: density !! Fraction of nonzero recurrent weights.
      real(dp), intent(in), optional :: lambda_range(:) !! Lower and upper ridge-penalty bounds.
      real(dp), intent(in), optional :: input_weight_scale !! Symmetric input-weight bound.
      real(dp), intent(in), optional :: reservoir_weight_scale !! Symmetric recurrent-weight bound.
      real(dp), intent(in), optional :: scale_range(:) !! Lower and upper transformed-data bounds.
      type(echos_model_t) :: model
      type(echos_ridge_t) :: ridge_fit, best_fit
      real(dp), allocatable :: differenced(:), scaled(:), inputs(:, :)
      real(dp), allocatable :: design(:, :), target(:)
      real(dp) :: selected_tau, lambda_lower, lambda_upper
      real(dp) :: input_scale, reservoir_scale, criterion, best_criterion
      integer :: observations, maximum_lag, training_rows, retained_rows
      integer :: candidate, selected_seed, target_start

      observations = size(series)
      model%lags = [1]
      if (present(lags)) model%lags = lags
      model%information_criterion = echos_ic_bic
      if (present(information_criterion)) &
         model%information_criterion = information_criterion
      model%difference_order = -1
      if (present(difference_order)) model%difference_order = difference_order
      selected_tau = 0.4_dp
      if (present(tau)) selected_tau = tau
      model%state_count = min(floor(real(observations, dp)*selected_tau), 200)
      if (present(state_count)) model%state_count = state_count
      model%model_count = 2*model%state_count
      if (present(model_count)) model%model_count = model_count
      model%initial_count = floor(0.05_dp*real(observations, dp))
      if (present(initial_count)) model%initial_count = initial_count
      model%alpha = 1.0_dp
      if (present(alpha)) model%alpha = alpha
      model%rho = 1.0_dp
      if (present(rho)) model%rho = rho
      model%density = 0.5_dp
      if (present(density)) model%density = density
      lambda_lower = 1.0e-4_dp
      lambda_upper = 2.0_dp
      if (present(lambda_range)) then
         if (size(lambda_range) == 2) then
            lambda_lower = lambda_range(1)
            lambda_upper = lambda_range(2)
         end if
      end if
      input_scale = 0.5_dp
      if (present(input_weight_scale)) input_scale = input_weight_scale
      reservoir_scale = 0.5_dp
      if (present(reservoir_weight_scale)) &
         reservoir_scale = reservoir_weight_scale
      model%scaled_minimum = -0.5_dp
      model%scaled_maximum = 0.5_dp
      if (present(scale_range)) then
         if (size(scale_range) == 2) then
            model%scaled_minimum = scale_range(1)
            model%scaled_maximum = scale_range(2)
         end if
      end if
      if (model%difference_order < 0) &
         model%difference_order = echos_estimate_differences(series)
      maximum_lag = 0
      if (size(model%lags) > 0) maximum_lag = maxval(model%lags)
      if (observations < 5 .or. size(model%lags) < 1 .or. &
         any(model%lags < 1) .or. model%difference_order < 0 .or. &
         observations - model%difference_order - maximum_lag - &
         model%initial_count < 2 .or. model%state_count < 1 .or. &
         model%model_count < 1 .or. model%initial_count < 0 .or. &
         model%alpha <= 0.0_dp .or. model%alpha > 1.0_dp .or. &
         model%rho <= 0.0_dp .or. selected_tau <= 0.0_dp .or. &
         selected_tau > 1.0_dp .or. model%density <= 0.0_dp .or. &
         model%density > 1.0_dp .or. lambda_lower <= 0.0_dp .or. &
         lambda_upper <= lambda_lower .or. input_scale <= 0.0_dp .or. &
         reservoir_scale <= 0.0_dp .or. &
         model%scaled_maximum <= model%scaled_minimum .or. &
         model%information_criterion < echos_ic_aic .or. &
         model%information_criterion > echos_ic_hqc .or. &
         .not. all(ieee_is_finite(series))) then
         model%info = 1
         return
      end if
      selected_seed = 42
      if (present(seed)) selected_seed = seed
      call set_random_seed(selected_seed)
      model%series = series
      differenced = difference_series(series, model%difference_order)
      model%original_minimum = minval(differenced)
      model%original_maximum = maxval(differenced)
      if (model%original_maximum <= model%original_minimum) then
         model%info = 2
         return
      end if
      scaled = scale_values(differenced, model%original_minimum, &
         model%original_maximum, model%scaled_minimum, model%scaled_maximum)
      model%transformed = scaled
      inputs = lag_matrix(scaled, model%lags)
      training_rows = size(inputs, 1)
      model%input_weights = random_uniform_matrix(model%state_count, &
         size(model%lags), input_scale)
      model%reservoir_weights = random_reservoir(model%state_count, model%rho, &
         model%density, reservoir_scale)
      if (size(model%reservoir_weights, 1) == 0) then
         model%info = 3
         return
      end if
      model%states = echos_run_reservoir(inputs, model%input_weights, &
         model%reservoir_weights, model%alpha)
      retained_rows = training_rows - model%initial_count
      allocate(design(retained_rows, model%state_count + 1))
      design(:, 1) = 1.0_dp
      design(:, 2:) = model%states(model%initial_count + 1:, :)
      target_start = maximum_lag + model%initial_count + 1
      target = scaled(target_start:)
      allocate(model%candidate_lambda(model%model_count))
      allocate(model%candidate_criterion(model%model_count))
      allocate(model%candidate_log_likelihood(model%model_count), &
         source=-huge(1.0_dp))
      allocate(model%candidate_degrees_of_freedom(model%model_count), &
         source=0.0_dp)
      allocate(model%candidate_aic(model%model_count), source=huge(1.0_dp))
      allocate(model%candidate_aicc(model%model_count), source=huge(1.0_dp))
      allocate(model%candidate_bic(model%model_count), source=huge(1.0_dp))
      allocate(model%candidate_hqc(model%model_count), source=huge(1.0_dp))
      allocate(model%candidate_mse(model%model_count), source=huge(1.0_dp))
      allocate(model%candidate_mae(model%model_count), source=huge(1.0_dp))
      best_criterion = huge(1.0_dp)
      do candidate = 1, model%model_count
         model%candidate_lambda(candidate) = lambda_lower + &
            (lambda_upper - lambda_lower)*random_uniform()
         ridge_fit = echos_fit_ridge(design, target, &
            model%candidate_lambda(candidate))
         if (ridge_fit%info /= 0) then
            model%candidate_criterion(candidate) = huge(1.0_dp)
            cycle
         end if
         criterion = ridge_criterion(ridge_fit, model%information_criterion)
         model%candidate_criterion(candidate) = criterion
         model%candidate_log_likelihood(candidate) = ridge_fit%log_likelihood
         model%candidate_degrees_of_freedom(candidate) = &
            ridge_fit%degrees_of_freedom
         model%candidate_aic(candidate) = ridge_fit%aic
         model%candidate_aicc(candidate) = ridge_fit%aicc
         model%candidate_bic(candidate) = ridge_fit%bic
         model%candidate_hqc(candidate) = ridge_fit%hqc
         model%candidate_mse(candidate) = ridge_fit%mse
         model%candidate_mae(candidate) = ridge_fit%mae
         if (criterion < best_criterion) then
            best_criterion = criterion
            best_fit = ridge_fit
         end if
      end do
      if (best_criterion >= 0.5_dp*huge(1.0_dp)) then
         model%info = 4
         return
      end if
      model%output_weights = best_fit%coefficients
      model%lambda = best_fit%lambda
      model%degrees_of_freedom = best_fit%degrees_of_freedom
      model%scaled_residuals = best_fit%residuals
      allocate(model%actual(observations), source=quiet_nan())
      allocate(model%fitted(observations), source=quiet_nan())
      allocate(model%residuals(observations), source=quiet_nan())
      model%actual(target_start + model%difference_order:) = &
         series(target_start + model%difference_order:)
      model%fitted(target_start + model%difference_order:) = &
         restore_one_step_fitted(series, unscale_values(best_fit%fitted, &
         model%original_minimum, model%original_maximum, &
         model%scaled_minimum, model%scaled_maximum), &
         model%difference_order, target_start + model%difference_order)
      model%residuals(target_start + model%difference_order:) = &
         model%actual(target_start + model%difference_order:) - &
         model%fitted(target_start + model%difference_order:)
   end function echos_fit

   function echos_forecast(model, horizon, levels, simulations, seed) result(out)
      !! Recursively forecast an ESN and form moving-block-bootstrap intervals.
      type(echos_model_t), intent(in) :: model !! Trained echos ESN.
      integer, intent(in) :: horizon !! Positive forecast horizon.
      real(dp), intent(in), optional :: levels(:) !! Central interval levels in percent.
      integer, intent(in), optional :: simulations !! Bootstrap path count; zero omits intervals.
      integer, intent(in), optional :: seed !! Shared random-number seed.
      type(echos_forecast_t) :: out
      real(dp), allocatable :: scaled_point(:), innovations(:), centered(:)
      real(dp), allocatable :: path(:), ordered(:)
      integer :: simulation, level, level_count, selected_seed, block_size

      out%horizon = horizon
      out%simulations = 100
      if (present(simulations)) out%simulations = simulations
      out%levels = [80.0_dp, 95.0_dp]
      if (present(levels)) out%levels = levels
      if (model%info /= 0 .or. horizon < 1 .or. out%simulations < 0 .or. &
         size(out%levels) < 1 .or. any(out%levels <= 0.0_dp) .or. &
         any(out%levels >= 100.0_dp)) then
         out%info = 1
         return
      end if
      selected_seed = 42
      if (present(seed)) selected_seed = seed
      call set_random_seed(selected_seed)
      call sort(out%levels)
      out%actual = model%actual
      out%fitted = model%fitted
      scaled_point = forecast_scaled(model, horizon)
      out%point = integrate_forecast(model%series, &
         unscale_values(scaled_point, model%original_minimum, &
         model%original_maximum, model%scaled_minimum, model%scaled_maximum), &
         model%difference_order)
      if (out%simulations == 0) then
         allocate(out%interval(horizon, 0))
         allocate(out%simulation(horizon, 0))
         allocate(out%standard_deviation(0))
         return
      end if
      centered = model%scaled_residuals - &
         sum(model%scaled_residuals)/real(size(model%scaled_residuals), dp)
      block_size = max(1, floor(real(size(centered), dp)**(1.0_dp/3.0_dp)))
      allocate(out%simulation(horizon, out%simulations))
      do simulation = 1, out%simulations
         innovations = moving_block_draw(centered, horizon, block_size)
         path = forecast_scaled(model, horizon, innovations)
         out%simulation(:, simulation) = integrate_forecast(model%series, &
            unscale_values(path, model%original_minimum, model%original_maximum, &
            model%scaled_minimum, model%scaled_maximum), model%difference_order)
      end do
      level_count = size(out%levels)
      allocate(out%interval(horizon, 2*level_count))
      allocate(out%standard_deviation(horizon))
      do simulation = 1, horizon
         ordered = out%simulation(simulation, :)
         call sort(ordered)
         do level = 1, level_count
            out%interval(simulation, level) = quantile(ordered, &
               0.5_dp - out%levels(level)/200.0_dp)
            out%interval(simulation, level_count + level) = quantile(ordered, &
               0.5_dp + out%levels(level)/200.0_dp)
         end do
         out%standard_deviation(simulation) = standard_deviation(ordered)
      end do
   end function echos_forecast

   function echos_tune(series, horizon, split_count, alphas, rhos, taus, &
      minimum_train, lags, information_criterion, difference_order, &
      model_count, state_count, initial_count, seed, density, lambda_range, &
      input_weight_scale, reservoir_weight_scale, scale_range) result(out)
      !! Tune leakage, spectral radius, and reservoir size by rolling forecasts.
      real(dp), intent(in) :: series(:) !! Finite univariate observations.
      integer, intent(in), optional :: horizon !! Forecast horizon; defaults to 12.
      integer, intent(in), optional :: split_count !! Requested expanding-window splits.
      real(dp), intent(in), optional :: alphas(:) !! Candidate leakage rates.
      real(dp), intent(in), optional :: rhos(:) !! Candidate spectral radii.
      real(dp), intent(in), optional :: taus(:) !! Candidate reservoir-size fractions.
      integer, intent(in), optional :: minimum_train !! Minimum first training size.
      integer, intent(in), optional :: lags(:) !! Positive autoregressive input lags.
      integer, intent(in), optional :: information_criterion !! Ridge selection criterion code.
      integer, intent(in), optional :: difference_order !! Ordinary difference order.
      integer, intent(in), optional :: model_count !! Ridge candidates per fit.
      integer, intent(in), optional :: state_count !! Fixed reservoir states, overriding tau.
      integer, intent(in), optional :: initial_count !! Initial states discarded.
      integer, intent(in), optional :: seed !! Shared random-number seed.
      real(dp), intent(in), optional :: density !! Fraction of nonzero recurrent weights.
      real(dp), intent(in), optional :: lambda_range(:) !! Ridge-penalty bounds.
      real(dp), intent(in), optional :: input_weight_scale !! Symmetric input-weight bound.
      real(dp), intent(in), optional :: reservoir_weight_scale !! Symmetric recurrent-weight bound.
      real(dp), intent(in), optional :: scale_range(:) !! Transformed-data scaling interval.
      type(echos_tuning_t) :: out
      type(echos_model_t) :: model
      type(echos_forecast_t) :: forecast
      real(dp), allocatable :: alpha_grid(:), rho_grid(:), tau_grid(:)
      real(dp), allocatable :: mean_mse(:)
      integer, allocatable :: train_ends(:)
      real(dp) :: error
      integer :: requested_splits, minimum, feasible, last_train, first
      integer :: configurations, rows, configuration, split, row
      integer :: alpha_index, rho_index, tau_index, selected_horizon

      selected_horizon = 12
      if (present(horizon)) selected_horizon = horizon
      requested_splits = 5
      if (present(split_count)) requested_splits = split_count
      minimum = max(30, 2*selected_horizon)
      if (present(minimum_train)) minimum = minimum_train
      alpha_grid = [(0.1_dp*real(row, dp), row = 1, 10)]
      if (present(alphas)) alpha_grid = alphas
      rho_grid = [(0.1_dp*real(row, dp), row = 1, 10)]
      if (present(rhos)) rho_grid = rhos
      tau_grid = [0.1_dp, 0.2_dp, 0.4_dp]
      if (present(taus)) tau_grid = taus
      out%horizon = selected_horizon
      out%actual = series
      if (selected_horizon < 1 .or. requested_splits < 1 .or. &
         minimum < 3 .or. size(series) < minimum + selected_horizon .or. &
         size(alpha_grid) < 1 .or. size(rho_grid) < 1 .or. &
         size(tau_grid) < 1 .or. any(alpha_grid <= 0.0_dp) .or. &
         any(alpha_grid > 1.0_dp) .or. any(rho_grid <= 0.0_dp) .or. &
         any(tau_grid <= 0.0_dp) .or. any(tau_grid > 1.0_dp) .or. &
         .not. all(ieee_is_finite(series))) then
         out%info = 1
         return
      end if
      last_train = size(series) - selected_horizon
      first = last_train - (requested_splits - 1)*selected_horizon
      feasible = 0
      do split = 1, requested_splits
         if (first + (split - 1)*selected_horizon >= minimum) &
            feasible = feasible + 1
      end do
      if (feasible < 1) then
         out%info = 2
         return
      end if
      allocate(train_ends(feasible))
      row = 0
      do split = 1, requested_splits
         last_train = first + (split - 1)*selected_horizon
         if (last_train < minimum) cycle
         row = row + 1
         train_ends(row) = last_train
      end do
      configurations = size(alpha_grid)*size(rho_grid)*size(tau_grid)
      rows = configurations*feasible
      allocate(out%alpha(rows), out%rho(rows), out%tau(rows))
      allocate(out%mse(rows), source=huge(1.0_dp))
      allocate(out%mae(rows), source=huge(1.0_dp))
      allocate(out%forecasts(rows, selected_horizon), source=quiet_nan())
      allocate(out%split(rows), out%train_end(rows), out%test_start(rows))
      allocate(out%test_end(rows), out%configuration(rows))
      configuration = 0
      row = 0
      do tau_index = 1, size(tau_grid)
         do rho_index = 1, size(rho_grid)
            do alpha_index = 1, size(alpha_grid)
               configuration = configuration + 1
               do split = 1, feasible
                  row = row + 1
                  out%configuration(row) = configuration
                  out%alpha(row) = alpha_grid(alpha_index)
                  out%rho(row) = rho_grid(rho_index)
                  out%tau(row) = tau_grid(tau_index)
                  out%split(row) = split
                  out%train_end(row) = train_ends(split)
                  out%test_start(row) = train_ends(split) + 1
                  out%test_end(row) = train_ends(split) + selected_horizon
                  model = echos_fit(series(:train_ends(split)), lags, &
                     information_criterion, difference_order, model_count, &
                     state_count, initial_count, seed, alpha_grid(alpha_index), &
                     rho_grid(rho_index), tau_grid(tau_index), density, &
                     lambda_range, input_weight_scale, reservoir_weight_scale, &
                     scale_range)
                  if (model%info /= 0) cycle
                  forecast = echos_forecast(model, selected_horizon, &
                     simulations=0, seed=seed)
                  if (forecast%info /= 0) cycle
                  out%forecasts(row, :) = forecast%point
                  error = sum((series(out%test_start(row):out%test_end(row)) - &
                     forecast%point)**2)/real(selected_horizon, dp)
                  out%mse(row) = error
                  out%mae(row) = sum(abs(series(out%test_start(row): &
                     out%test_end(row)) - forecast%point))/ &
                     real(selected_horizon, dp)
               end do
            end do
         end do
      end do
      allocate(mean_mse(configurations), source=huge(1.0_dp))
      do configuration = 1, configurations
         if (all(out%mse < 0.5_dp*huge(1.0_dp) .or. &
            out%configuration /= configuration)) then
            mean_mse(configuration) = sum(out%mse, &
               mask=out%configuration == configuration)/real(feasible, dp)
         end if
      end do
      out%best = minloc(mean_mse, dim=1)
      if (mean_mse(out%best) >= 0.5_dp*huge(1.0_dp)) out%info = 3
   end function echos_tune

   pure real(dp) function trace_matrix(matrix) result(value)
      !! Return the trace of a square matrix.
      real(dp), intent(in) :: matrix(:, :) !! Square matrix.
      integer :: index

      value = 0.0_dp
      do index = 1, min(size(matrix, 1), size(matrix, 2))
         value = value + matrix(index, index)
      end do
   end function trace_matrix

   pure logical function valid_kpss_alpha(alpha) result(valid)
      !! Test whether alpha is one of the echos KPSS significance levels.
      real(dp), intent(in) :: alpha !! Candidate significance level.

      valid = abs(alpha - 0.10_dp) < 1.0e-12_dp .or. &
         abs(alpha - 0.05_dp) < 1.0e-12_dp .or. &
         abs(alpha - 0.025_dp) < 1.0e-12_dp .or. &
         abs(alpha - 0.01_dp) < 1.0e-12_dp
   end function valid_kpss_alpha

   pure real(dp) function kpss_critical_value(trend, alpha) result(value)
      !! Return the package KPSS critical value for a supported test size.
      logical, intent(in) :: trend !! Include a deterministic trend.
      real(dp), intent(in) :: alpha !! Supported significance level.

      if (trend) then
         if (abs(alpha - 0.10_dp) < 1.0e-12_dp) value = 0.119_dp
         if (abs(alpha - 0.05_dp) < 1.0e-12_dp) value = 0.146_dp
         if (abs(alpha - 0.025_dp) < 1.0e-12_dp) value = 0.176_dp
         if (abs(alpha - 0.01_dp) < 1.0e-12_dp) value = 0.216_dp
      else
         if (abs(alpha - 0.10_dp) < 1.0e-12_dp) value = 0.347_dp
         if (abs(alpha - 0.05_dp) < 1.0e-12_dp) value = 0.463_dp
         if (abs(alpha - 0.025_dp) < 1.0e-12_dp) value = 0.574_dp
         if (abs(alpha - 0.01_dp) < 1.0e-12_dp) value = 0.739_dp
      end if
   end function kpss_critical_value

   pure function difference_series(series, order) result(transformed)
      !! Apply consecutive ordinary first differences.
      real(dp), intent(in) :: series(:) !! Input observations.
      integer, intent(in) :: order !! Nonnegative difference order.
      real(dp), allocatable :: transformed(:), next(:)
      integer :: stage

      transformed = series
      do stage = 1, order
         allocate(next(size(transformed) - 1))
         next = transformed(2:) - transformed(:size(transformed) - 1)
         call move_alloc(next, transformed)
      end do
   end function difference_series

   pure function scale_values(values, old_minimum, old_maximum, new_minimum, &
      new_maximum) result(scaled)
      !! Scale values linearly between two finite intervals.
      real(dp), intent(in) :: values(:) !! Values in the original interval.
      real(dp), intent(in) :: old_minimum !! Original lower bound.
      real(dp), intent(in) :: old_maximum !! Original upper bound.
      real(dp), intent(in) :: new_minimum !! New lower bound.
      real(dp), intent(in) :: new_maximum !! New upper bound.
      real(dp), allocatable :: scaled(:)

      scaled = new_minimum + (values - old_minimum)* &
         (new_maximum - new_minimum)/(old_maximum - old_minimum)
   end function scale_values

   pure function unscale_values(values, old_minimum, old_maximum, new_minimum, &
      new_maximum) result(original)
      !! Reverse a linear interval scaling transformation.
      real(dp), intent(in) :: values(:) !! Values in the scaled interval.
      real(dp), intent(in) :: old_minimum !! Original lower bound.
      real(dp), intent(in) :: old_maximum !! Original upper bound.
      real(dp), intent(in) :: new_minimum !! Scaled lower bound.
      real(dp), intent(in) :: new_maximum !! Scaled upper bound.
      real(dp), allocatable :: original(:)

      original = old_minimum + (values - new_minimum)* &
         (old_maximum - old_minimum)/(new_maximum - new_minimum)
   end function unscale_values

   pure function lag_matrix(series, lags) result(inputs)
      !! Construct complete lagged input rows from a transformed series.
      real(dp), intent(in) :: series(:) !! Transformed observations.
      integer, intent(in) :: lags(:) !! Positive lag indices.
      real(dp), allocatable :: inputs(:, :)
      integer :: maximum, row, lag

      maximum = maxval(lags)
      allocate(inputs(size(series) - maximum, size(lags)))
      do lag = 1, size(lags)
         do row = 1, size(inputs, 1)
            inputs(row, lag) = series(maximum + row - lags(lag))
         end do
      end do
   end function lag_matrix

   function random_uniform_matrix(rows, columns, bound) result(matrix)
      !! Draw a dense matrix from a symmetric uniform interval.
      integer, intent(in) :: rows !! Matrix rows.
      integer, intent(in) :: columns !! Matrix columns.
      real(dp), intent(in) :: bound !! Positive symmetric interval bound.
      real(dp), allocatable :: matrix(:, :)
      integer :: row, column

      allocate(matrix(rows, columns))
      do column = 1, columns
         do row = 1, rows
            matrix(row, column) = bound*(2.0_dp*random_uniform() - 1.0_dp)
         end do
      end do
   end function random_uniform_matrix

   function random_reservoir(state_count, rho, density, bound) result(matrix)
      !! Draw an exact-density reservoir and scale it to a spectral radius.
      integer, intent(in) :: state_count !! Number of reservoir states.
      real(dp), intent(in) :: rho !! Target spectral radius.
      real(dp), intent(in) :: density !! Fraction of retained weights.
      real(dp), intent(in) :: bound !! Symmetric uniform weight bound.
      real(dp), allocatable :: matrix(:, :)
      real(dp), allocatable :: flat(:)
      complex(dp), allocatable :: eigenvalues(:)
      integer, allocatable :: positions(:)
      real(dp) :: radius
      integer :: total, nonzero, index, other, held, status, attempt

      total = state_count*state_count
      nonzero = max(1, min(total, nint(real(total, dp)*density)))
      allocate(flat(total), positions(total))
      do attempt = 1, 20
         do index = 1, total
            flat(index) = bound*(2.0_dp*random_uniform() - 1.0_dp)
            positions(index) = index
         end do
         do index = total, 2, -1
            other = 1 + int(random_uniform()*real(index, dp))
            other = min(index, other)
            held = positions(index)
            positions(index) = positions(other)
            positions(other) = held
         end do
         if (nonzero < total) flat(positions(nonzero + 1:)) = 0.0_dp
         matrix = reshape(flat, [state_count, state_count])
         call general_eigenvalues(matrix, eigenvalues, status)
         if (status /= 0) cycle
         radius = maxval(abs(eigenvalues))
         if (radius > sqrt(tiny(1.0_dp))) then
            matrix = matrix*(rho/radius)
            return
         end if
      end do
      deallocate(matrix)
      allocate(matrix(0, 0))
   end function random_reservoir

   pure real(dp) function ridge_criterion(fit, criterion) result(value)
      !! Select one information criterion from ridge-fit metrics.
      type(echos_ridge_t), intent(in) :: fit !! Ridge readout fit.
      integer, intent(in) :: criterion !! Information-criterion code.

      select case (criterion)
      case (echos_ic_aic)
         value = fit%aic
      case (echos_ic_aicc)
         value = fit%aicc
      case (echos_ic_bic)
         value = fit%bic
      case default
         value = fit%hqc
      end select
   end function ridge_criterion

   pure function restore_one_step_fitted(series, differences, order, &
      first_index) result(fitted)
      !! Restore one-step fitted levels from fitted ordinary differences.
      real(dp), intent(in) :: series(:) !! Original observed levels.
      real(dp), intent(in) :: differences(:) !! Fitted order-d differences.
      integer, intent(in) :: order !! Difference order.
      integer, intent(in) :: first_index !! Original index of the first fit.
      real(dp), allocatable :: fitted(:)
      real(dp) :: value
      integer :: observation, lag

      allocate(fitted(size(differences)))
      do observation = 1, size(differences)
         value = differences(observation)
         do lag = 1, order
            value = value + (-1.0_dp)**(lag + 1)* &
               real(binomial_coefficient(order, lag), dp)* &
               series(first_index + observation - 1 - lag)
         end do
         fitted(observation) = value
      end do
   end function restore_one_step_fitted

   pure function forecast_scaled(model, horizon, innovations) result(forecast)
      !! Recursively forecast in the differenced and scaled model domain.
      type(echos_model_t), intent(in) :: model !! Trained ESN model.
      integer, intent(in) :: horizon !! Positive forecast horizon.
      real(dp), intent(in), optional :: innovations(:) !! Additive scaled innovations.
      real(dp), allocatable :: forecast(:)
      real(dp), allocatable :: history(:), state(:), input(:), update(:)
      real(dp) :: prediction
      integer :: observations, step, lag

      observations = size(model%transformed)
      allocate(history(observations + horizon))
      history(:observations) = model%transformed
      state = model%states(size(model%states, 1), :)
      allocate(input(size(model%lags)), forecast(horizon))
      do step = 1, horizon
         do lag = 1, size(model%lags)
            input(lag) = history(observations + step - model%lags(lag))
         end do
         update = tanh(matmul(model%input_weights, input) + &
            matmul(model%reservoir_weights, state))
         state = model%alpha*update + (1.0_dp - model%alpha)*state
         prediction = model%output_weights(1) + &
            dot_product(model%output_weights(2:), state)
         if (present(innovations)) prediction = prediction + innovations(step)
         forecast(step) = prediction
         history(observations + step) = prediction
      end do
   end function forecast_scaled

   pure function integrate_forecast(series, differences, order) result(levels)
      !! Recursively integrate future ordinary differences to the level scale.
      real(dp), intent(in) :: series(:) !! Historical observed levels.
      real(dp), intent(in) :: differences(:) !! Future order-d differences.
      integer, intent(in) :: order !! Difference order.
      real(dp), allocatable :: levels(:), history(:)
      real(dp) :: value
      integer :: observations, step, lag

      observations = size(series)
      allocate(history(observations + size(differences)))
      history(:observations) = series
      allocate(levels(size(differences)))
      do step = 1, size(differences)
         value = differences(step)
         do lag = 1, order
            value = value + (-1.0_dp)**(lag + 1)* &
               real(binomial_coefficient(order, lag), dp)* &
               history(observations + step - lag)
         end do
         levels(step) = value
         history(observations + step) = value
      end do
   end function integrate_forecast

   function moving_block_draw(values, count, block_size) result(draw)
      !! Draw one noncircular moving-block bootstrap sequence.
      real(dp), intent(in) :: values(:) !! Centered residual values.
      integer, intent(in) :: count !! Required output length.
      integer, intent(in) :: block_size !! Positive overlapping-block size.
      real(dp), allocatable :: draw(:)
      integer :: available, output, start, offset

      allocate(draw(count))
      available = size(values) - block_size + 1
      output = 0
      do while (output < count)
         start = 1 + int(random_uniform()*real(available, dp))
         start = min(available, start)
         do offset = 0, block_size - 1
            if (output == count) exit
            output = output + 1
            draw(output) = values(start + offset)
         end do
      end do
   end function moving_block_draw

   pure integer function binomial_coefficient(n, k) result(value)
      !! Return a small integer binomial coefficient.
      integer, intent(in) :: n !! Nonnegative upper argument.
      integer, intent(in) :: k !! Lower argument between zero and n.
      integer :: index, selected

      selected = min(k, n - k)
      value = 1
      do index = 1, selected
         value = value*(n - selected + index)/index
      end do
   end function binomial_coefficient

   subroutine display_echos_model(model, unit, print_obs)
      !! Display a trained ESN and optionally its aligned fitted observations.
      type(echos_model_t), intent(in) :: model !! Trained echos model.
      integer, intent(in), optional :: unit !! Output unit; defaults to standard output.
      logical, intent(in), optional :: print_obs !! Print actual, fitted, and residual values.
      integer :: destination, observation
      logical :: show_observations

      destination = output_unit
      if (present(unit)) destination = unit
      show_observations = .false.
      if (present(print_obs)) show_observations = print_obs
      write(destination, '(a)') 'Echo State Network fit'
      write(destination, '(a, i0)') '  status: ', model%info
      write(destination, '(a, i0)') '  reservoir states: ', model%state_count
      write(destination, '(a, *(i0, 1x))') '  input lags: ', model%lags
      write(destination, '(a, i0)') '  difference order: ', &
         model%difference_order
      write(destination, '(a, f8.4)') '  leakage rate: ', model%alpha
      write(destination, '(a, f8.4)') '  spectral radius: ', model%rho
      write(destination, '(a, f8.4)') '  reservoir density: ', model%density
      write(destination, '(a, es14.6)') '  selected ridge penalty: ', &
         model%lambda
      write(destination, '(a, f10.4)') '  effective degrees of freedom: ', &
         model%degrees_of_freedom
      if (.not. show_observations .or. model%info /= 0) return
      write(destination, '(a)') '  index          actual         fitted        residual'
      do observation = 1, size(model%series)
         write(destination, '(i7, 3(1x, es14.6))') observation, &
            model%actual(observation), model%fitted(observation), &
            model%residuals(observation)
      end do
   end subroutine display_echos_model

   subroutine display_echos_forecast(forecast, unit, print_obs)
      !! Display ESN forecast metadata and optionally point forecasts.
      type(echos_forecast_t), intent(in) :: forecast !! ESN forecast result.
      integer, intent(in), optional :: unit !! Output unit; defaults to standard output.
      logical, intent(in), optional :: print_obs !! Print horizon-specific forecasts.
      integer :: destination, step
      logical :: show_observations

      destination = output_unit
      if (present(unit)) destination = unit
      show_observations = .false.
      if (present(print_obs)) show_observations = print_obs
      write(destination, '(a)') 'Echo State Network forecast'
      write(destination, '(a, i0)') '  status: ', forecast%info
      write(destination, '(a, i0)') '  horizon: ', forecast%horizon
      write(destination, '(a, i0)') '  bootstrap simulations: ', &
         forecast%simulations
      if (.not. show_observations .or. forecast%info /= 0) return
      write(destination, '(a)') '  step           point'
      do step = 1, forecast%horizon
         write(destination, '(i7, 1x, es14.6)') step, forecast%point(step)
      end do
   end subroutine display_echos_forecast

   subroutine display_echos_tuning(tuning, unit, print_obs)
      !! Display ESN tuning dimensions and optionally all rolling errors.
      type(echos_tuning_t), intent(in) :: tuning !! Expanding-window tuning result.
      integer, intent(in), optional :: unit !! Output unit; defaults to standard output.
      logical, intent(in), optional :: print_obs !! Print every grid and split result.
      integer :: destination, row
      logical :: show_observations

      destination = output_unit
      if (present(unit)) destination = unit
      show_observations = .false.
      if (present(print_obs)) show_observations = print_obs
      write(destination, '(a)') 'Echo State Network tuning'
      write(destination, '(a, i0)') '  status: ', tuning%info
      write(destination, '(a, i0)') '  best configuration: ', tuning%best
      write(destination, '(a, i0)') '  evaluated rows: ', size(tuning%mse)
      if (.not. show_observations .or. tuning%info /= 0) return
      write(destination, '(a)') &
         ' config split      alpha        rho        tau          mse          mae'
      do row = 1, size(tuning%mse)
         write(destination, '(2i7, 5(1x, es12.4))') &
            tuning%configuration(row), tuning%split(row), tuning%alpha(row), &
            tuning%rho(row), tuning%tau(row), tuning%mse(row), tuning%mae(row)
      end do
   end subroutine display_echos_tuning

end module echos_mod
