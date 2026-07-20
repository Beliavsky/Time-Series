! SPDX-License-Identifier: BSD-3-Clause
! SPDX-FileComment: Algorithms translated from the R robustarima package.
! Filtered tau-estimation algorithms translated from robustarima.
module robustarima_mod
   use kind_mod, only: dp
   use tfarima_mod, only: tfarima_forecast_t, tfarima_arima_forecast, &
      tfarima_polynomial_multiply, tfarima_polynomial_power, &
      tfarima_polynomial_ratio, tfarima_difference, &
      tfarima_outlier_response
   use linalg_mod, only: invert_matrix
   use stats_mod, only: median
   use optimization_mod, only: optimization_result_t, &
      bfgs_minimize_fd, finite_difference_hessian
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   private

   integer, parameter, public :: robustarima_outlier_io = 1
   integer, parameter, public :: robustarima_outlier_ao = 2
   integer, parameter, public :: robustarima_outlier_ls = 3

   type, public :: robustarima_filter_t
      ! Bounded-innovation filter paths and robust scale diagnostics.
      real(dp), allocatable :: innovations(:), bounded_innovations(:)
      real(dp), allocatable :: cleaned(:), fitted(:), weights(:)
      real(dp), allocatable :: standardized_innovations(:)
      real(dp), allocatable :: prediction_scale(:)
      real(dp), allocatable :: state(:, :), state_covariance(:, :, :)
      logical, allocatable :: level_shift_candidate(:)
      real(dp) :: initial_scale = 0.0_dp
      real(dp) :: innovation_scale = 0.0_dp
      real(dp) :: tau_scale = 0.0_dp
      real(dp) :: objective = huge(1.0_dp)
      integer :: initial_observations = 0
      integer :: info = 0
   end type robustarima_filter_t

   type, public :: robustarima_fit_t
      ! Filtered tau estimates for regression with ARIMA errors.
      real(dp), allocatable :: ar(:), ma(:), regression_coefficients(:)
      real(dp), allocatable :: difference_polynomial(:)
      real(dp), allocatable :: parameter_covariance(:, :)
      real(dp), allocatable :: standard_error(:)
      real(dp), allocatable :: regression_covariance(:, :)
      real(dp), allocatable :: regression_standard_error(:)
      real(dp), allocatable :: regression_residuals(:)
      type(robustarima_filter_t) :: filter
      real(dp) :: seasonal_ma = 0.0_dp
      real(dp) :: tuning_constant = 1.0_dp
      real(dp) :: robust_aic = huge(1.0_dp)
      real(dp) :: tau_inverse_efficiency = 0.0_dp
      integer :: difference_order = 0
      integer :: seasonal_period = 1
      integer :: seasonal_difference_order = 0
      integer :: iterations = 0
      integer :: info = 0
      logical :: includes_seasonal_ma = .false.
      logical :: converged = .false.
   end type robustarima_fit_t

   type, public :: robustarima_tau_inference_t
      ! Tau sandwich covariance, filtered design, and efficiency diagnostics.
      real(dp), allocatable :: covariance(:, :), standard_error(:)
      real(dp), allocatable :: filtered_design(:, :), implicit_weight(:)
      real(dp) :: m_scale = 0.0_dp
      real(dp) :: inverse_efficiency = 0.0_dp
      real(dp) :: rho1_weight = 0.0_dp
      integer :: info = 0
   end type robustarima_tau_inference_t

   type, public :: robustarima_forecast_t
      ! Robust ARIMA forecasts and Gaussian innovation-scale intervals.
      real(dp), allocatable :: mean(:), standard_error(:)
      real(dp), allocatable :: lower(:, :), upper(:, :)
      integer :: info = 0
   end type robustarima_forecast_t

   type, public :: robustarima_outliers_t
      ! Robust innovation, additive, and level-shift detections.
      integer, allocatable :: outlier_type(:), position(:)
      real(dp), allocatable :: effect(:), t_statistic(:)
      real(dp), allocatable :: cleaned_series(:)
      real(dp) :: scale_before = 0.0_dp
      real(dp) :: scale_after = 0.0_dp
      integer :: info = 0
   end type robustarima_outliers_t

   type, public :: robustarima_order_selection_t
      ! Robust information criteria and selected autoregressive fit.
      type(robustarima_fit_t) :: fit
      real(dp), allocatable :: criterion(:)
      integer :: selected_order = 0
      integer :: info = 0
   end type robustarima_order_selection_t

   public :: robustarima_rho, robustarima_psi, robustarima_psi_derivative
   public :: robustarima_m_scale, robustarima_tau_scale
   public :: robustarima_correlation_series
   public :: robustarima_bounded_filter, robustarima_fit
   public :: robustarima_select_ar_order
   public :: robustarima_tau_inference
   public :: robustarima_forecast, robustarima_detect_outliers

contains

   pure elemental real(dp) function robustarima_rho(value) result(rho)
      !! Evaluate robustarima's bounded, redescending rho function.
      real(dp), intent(in) :: value !! Input value.
      real(dp) :: absolute_value, squared

      absolute_value = abs(value)
      squared = value*value
      if (absolute_value < 2.0_dp) then
         rho = 0.5_dp*squared
      else if (absolute_value > 3.0_dp) then
         rho = 3.25_dp
      else
         rho = -0.972_dp*squared + 0.432_dp*squared**2 - &
            0.052_dp*squared**3 + 0.002_dp*squared**4 + 1.792_dp
      end if
   end function robustarima_rho

   pure elemental real(dp) function robustarima_psi(value) result(psi)
      !! Evaluate the derivative of the bounded rho function.
      real(dp), intent(in) :: value !! Input value.
      real(dp) :: squared

      if (abs(value) > 3.0_dp) then
         psi = 0.0_dp
      else if (abs(value) <= 2.0_dp) then
         psi = value
      else
         squared = value*value
         psi = value*(0.016_dp*squared**3 - 0.312_dp*squared**2 + &
            1.728_dp*squared - 1.944_dp)
      end if
   end function robustarima_psi

   pure elemental real(dp) function robustarima_psi_derivative(value) &
      result(derivative)
      !! Evaluate the derivative of robustarima's psi function.
      real(dp), intent(in) :: value !! Input value.
      real(dp) :: squared

      if (abs(value) > 3.0_dp) then
         derivative = 0.0_dp
      else if (abs(value) <= 2.0_dp) then
         derivative = 1.0_dp
      else
         squared = value*value
         derivative = 0.112_dp*squared**3 - 1.56_dp*squared**2 + &
            5.184_dp*squared - 1.944_dp
      end if
   end function robustarima_psi_derivative

   pure real(dp) function robustarima_m_scale(residuals, first) result(scale)
      !! Compute the package's bisquare M-scale initialized by the MAD.
      real(dp), intent(in) :: residuals(:) !! Model residuals.
      integer, intent(in), optional :: first !! First operand.
      real(dp), allocatable :: absolute_residual(:)
      real(dp) :: mad, old_scale, new_scale, rho_sum
      integer :: start, iteration, observations

      start = 1
      if (present(first)) start = first
      if (start < 1 .or. start > size(residuals) .or. &
         .not. all(ieee_is_finite(residuals))) then
         scale = 0.0_dp
         return
      end if
      absolute_residual = abs(residuals(start:))
      observations = size(absolute_residual)
      mad = median(absolute_residual)/0.6745_dp
      mad = max(mad, 1.0e-20_dp)
      old_scale = 1.0_dp
      do iteration = 1, 10000
         rho_sum = sum(robustarima_rho(absolute_residual/ &
            (old_scale*mad*0.405_dp)))
         if (rho_sum <= tiny(1.0_dp)) then
            scale = mad
            return
         end if
         new_scale = old_scale*sqrt(rho_sum/(real(observations, dp)*1.625_dp))
         if (abs(new_scale - old_scale) <= &
            1.0e-9_dp*max(old_scale, tiny(1.0_dp))) exit
         old_scale = max(new_scale, tiny(1.0_dp))
      end do
      scale = max(new_scale*mad, tiny(1.0_dp))
   end function robustarima_m_scale

   pure real(dp) function robustarima_tau_scale(residuals, first) result(scale)
      !! Compute the tau-scale associated with the filtered residuals.
      real(dp), intent(in) :: residuals(:) !! Model residuals.
      integer, intent(in), optional :: first !! First operand.
      real(dp) :: mscale
      integer :: start, observations

      start = 1
      if (present(first)) start = first
      if (start < 1 .or. start > size(residuals)) then
         scale = 0.0_dp
         return
      end if
      mscale = robustarima_m_scale(residuals, start)
      observations = size(residuals) - start + 1
      if (mscale <= tiny(1.0_dp)) then
         scale = 0.0_dp
         return
      end if
      scale = mscale*sqrt(sum(robustarima_rho( &
         residuals(start:)/mscale))/(0.488_dp*real(observations, dp)))
   end function robustarima_tau_scale

   pure function robustarima_correlation_series(residuals, scales, first) &
      result(pseudo_observations)
      !! Form the clipped pseudo-series used for robust ACF and PACF estimates.
      real(dp), intent(in) :: residuals(:) !! Model residuals.
      real(dp), intent(in), optional :: scales(:) !! Scales.
      integer, intent(in), optional :: first !! First operand.
      real(dp), allocatable :: pseudo_observations(:), standardized(:)
      real(dp) :: scale
      integer :: start

      start = 1
      if (present(first)) start = first
      if (start < 1 .or. start > size(residuals) .or. &
         .not. all(ieee_is_finite(residuals))) then
         allocate(pseudo_observations(0))
         return
      end if
      standardized = residuals(start:)
      if (present(scales)) then
         if (size(scales) /= size(residuals) .or. &
            any(scales(start:) <= 0.0_dp) .or. &
            .not. all(ieee_is_finite(scales))) then
            allocate(pseudo_observations(0))
            return
         end if
         standardized = standardized/scales(start:)
      end if
      scale = robustarima_m_scale(standardized)
      if (scale <= tiny(1.0_dp)) then
         allocate(pseudo_observations(0))
         return
      end if
      pseudo_observations = max(-2.5_dp, min(2.5_dp, standardized/scale))
   end function robustarima_correlation_series

   pure function robustarima_bounded_filter(series, ar, ma, &
      tuning_constant, scale, seasonal_period, seasonal_ma) result(out)
      !! Filter an ARMA series with bounded propagation of large innovations.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      real(dp), intent(in) :: ar(:) !! Autoregressive coefficients.
      real(dp), intent(in) :: ma(:) !! Moving-average coefficients.
      real(dp), intent(in), optional :: tuning_constant !! Tuning constant.
      real(dp), intent(in), optional :: scale !! Scale.
      integer, intent(in), optional :: seasonal_period !! Seasonal period.
      real(dp), intent(in), optional :: seasonal_ma !! Seasonal moving-average.
      type(robustarima_filter_t) :: out
      real(dp), allocatable :: combined_ma(:)
      real(dp), allocatable :: transition(:, :), disturbance(:)
      real(dp), allocatable :: state(:), covariance(:, :)
      real(dp), allocatable :: predicted_state(:), predicted_covariance(:, :)
      real(dp), allocatable :: gain_numerator(:), saved_state(:)
      real(dp), allocatable :: saved_covariance(:, :)
      real(dp) :: tuning, working_scale, innovation, bounded
      real(dp) :: innovation_variance, standardized, weight
      real(dp) :: scale_correction, geometric_variance, rho_mean
      integer :: time, p, q, initial, dimension, i
      integer :: consecutive_extremes, first_extreme

      tuning = 1.0_dp
      if (present(tuning_constant)) tuning = tuning_constant
      if (size(series) < 1 .or. tuning <= 0.0_dp .or. &
         .not. all(ieee_is_finite(series)) .or. &
         .not. all(ieee_is_finite(ar)) .or. .not. all(ieee_is_finite(ma))) then
         out%info = 1
         return
      end if
      combined_ma = seasonal_ma_coefficients(ma, seasonal_period, seasonal_ma)
      p = size(ar)
      q = size(combined_ma)
      dimension = max(1, p, q + 1)
      initial = min(size(series) - 1, dimension)
      working_scale = robustarima_m_scale(series)
      if (present(scale)) working_scale = scale
      if (working_scale <= tiny(1.0_dp) .or. &
         .not. ieee_is_finite(working_scale)) then
         out%info = 2
         return
      end if
      allocate(out%innovations(size(series)), &
         out%bounded_innovations(size(series)), out%cleaned(size(series)), &
         out%fitted(size(series)), out%weights(size(series)), &
         out%standardized_innovations(size(series)), &
         out%prediction_scale(size(series)), &
         out%state(dimension, size(series)), &
         out%state_covariance(dimension, dimension, size(series)), &
         out%level_shift_candidate(size(series)))
      allocate(transition(dimension, dimension), disturbance(dimension), &
         state(dimension), covariance(dimension, dimension), &
         predicted_state(dimension), &
         predicted_covariance(dimension, dimension), &
         gain_numerator(dimension), saved_state(dimension), &
         saved_covariance(dimension, dimension))
      out%innovations = 0.0_dp
      out%bounded_innovations = 0.0_dp
      out%cleaned = series
      out%fitted = 0.0_dp
      out%weights = 1.0_dp
      out%standardized_innovations = 0.0_dp
      out%prediction_scale = working_scale
      out%state = 0.0_dp
      out%state_covariance = 0.0_dp
      out%level_shift_candidate = .false.
      transition = 0.0_dp
      transition(:p, 1) = ar
      do i = 1, dimension - 1
         transition(i, i + 1) = 1.0_dp
      end do
      disturbance = 0.0_dp
      disturbance(1) = 1.0_dp
      if (q > 0) disturbance(2:q + 1) = -combined_ma
      state = 0.0_dp
      covariance = 0.0_dp
      consecutive_extremes = 0
      first_extreme = 0
      time = 1
      do while (time <= size(series))
         predicted_state = matmul(transition, state)
         predicted_covariance = matmul(matmul(transition, covariance), &
            transpose(transition)) + working_scale**2* &
            spread(disturbance, 2, dimension)* &
            spread(disturbance, 1, dimension)
         predicted_covariance = 0.5_dp*(predicted_covariance + &
            transpose(predicted_covariance))
         gain_numerator = predicted_covariance(:, 1)
         innovation_variance = max(predicted_covariance(1, 1), &
            tiny(1.0_dp))
         out%prediction_scale(time) = sqrt(innovation_variance)
         innovation = series(time) - predicted_state(1)
         standardized = innovation/out%prediction_scale(time)
         bounded = tuning*robustarima_psi(standardized/tuning)* &
            out%prediction_scale(time)
         weight = 1.0_dp
         if (abs(innovation) > 1.0e-10_dp) weight = bounded/innovation
         if (abs(standardized) > 2.5_dp) then
            consecutive_extremes = consecutive_extremes + 1
            if (consecutive_extremes == 1) then
               first_extreme = time
               saved_state = predicted_state
               saved_covariance = predicted_covariance
            else if (consecutive_extremes > 2) then
               out%level_shift_candidate(first_extreme) = .true.
               innovation = series(first_extreme) - saved_state(1)
               gain_numerator = saved_covariance(:, 1)
               innovation_variance = max(saved_covariance(1, 1), tiny(1.0_dp))
               state = saved_state + gain_numerator*innovation/innovation_variance
               covariance = saved_covariance - spread(gain_numerator, 2, &
                  dimension)*spread(gain_numerator, 1, dimension)/ &
                  innovation_variance
               covariance = 0.5_dp*(covariance + transpose(covariance))
               out%innovations(first_extreme) = innovation
               out%bounded_innovations(first_extreme) = innovation
               out%weights(first_extreme) = 1.0_dp
               out%cleaned(first_extreme) = state(1)
               out%state(:, first_extreme) = state
               out%state_covariance(:, :, first_extreme) = covariance
               time = first_extreme + 1
               consecutive_extremes = 0
               first_extreme = 0
               cycle
            end if
         else
            consecutive_extremes = 0
            first_extreme = 0
         end if
         state = predicted_state + weight*gain_numerator* &
            innovation/innovation_variance
         covariance = predicted_covariance - weight* &
            spread(gain_numerator, 2, dimension)* &
            spread(gain_numerator, 1, dimension)/innovation_variance
         covariance = 0.5_dp*(covariance + transpose(covariance))
         out%innovations(time) = innovation
         out%bounded_innovations(time) = bounded
         out%fitted(time) = predicted_state(1)
         out%cleaned(time) = state(1)
         out%weights(time) = weight
         out%standardized_innovations(time) = standardized
         out%state(:, time) = state
         out%state_covariance(:, :, time) = covariance
         time = time + 1
      end do
      out%initial_scale = working_scale
      out%initial_observations = initial
      scale_correction = robustarima_m_scale( &
         out%standardized_innovations, initial + 1)
      if (scale_correction <= tiny(1.0_dp)) then
         out%info = 3
         return
      end if
      out%innovation_scale = working_scale*scale_correction
      geometric_variance = exp(sum(log(max(out%prediction_scale( &
         initial + 1:)**2, tiny(1.0_dp))))/ &
         real(size(series) - initial, dp))
      rho_mean = sum(robustarima_rho(out%standardized_innovations( &
         initial + 1:)/scale_correction))/ &
         real(size(series) - initial, dp)
      out%tau_scale = sqrt(max(tiny(1.0_dp), geometric_variance* &
         scale_correction**2*rho_mean/0.488_dp))
      out%standardized_innovations = out%innovations/ &
         max(out%prediction_scale*scale_correction, &
         tiny(1.0_dp))
      out%objective = out%tau_scale**2
      if (.not. ieee_is_finite(out%objective)) out%info = 4
   end function robustarima_bounded_filter

   pure function robustarima_fit(series, ar_order, ma_order, &
      difference_order, regressors, seasonal_period, &
      seasonal_difference_order, include_seasonal_ma, tuning_constant, &
      max_iterations, tolerance) result(out)
      !! Fit regression with ARIMA errors by filtered tau minimization.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      integer, intent(in) :: ar_order !! Autoregressive order.
      integer, intent(in) :: ma_order !! Moving-average order.
      integer, intent(in), optional :: difference_order !! Difference order.
      integer, intent(in), optional :: seasonal_period !! Seasonal period.
      integer, intent(in), optional :: seasonal_difference_order !! Seasonal difference order.
      real(dp), intent(in), optional :: regressors(:, :) !! Regression design matrix.
      real(dp), intent(in), optional :: tuning_constant !! Tuning constant.
      logical, intent(in), optional :: include_seasonal_ma !! Whether to include the seasonal moving-average.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      type(robustarima_fit_t) :: out
      type(optimization_result_t) :: optimized, candidate
      type(robustarima_tau_inference_t) :: inference
      real(dp), allocatable :: initial(:), tuning_values(:), hessian(:, :)
      real(dp), allocatable :: inverse(:, :), parameters(:), stationary(:)
      real(dp) :: tol
      integer :: d, period, sd, limit, regression_count, parameter_count
      integer :: regression_start
      integer :: tuning_index, status, i
      logical :: seasonal

      d = 0
      if (present(difference_order)) d = difference_order
      period = 1
      if (present(seasonal_period)) period = seasonal_period
      sd = 0
      if (present(seasonal_difference_order)) sd = seasonal_difference_order
      seasonal = .false.
      if (present(include_seasonal_ma)) seasonal = include_seasonal_ma
      limit = 300
      if (present(max_iterations)) limit = max_iterations
      tol = 1.0e-6_dp
      if (present(tolerance)) tol = tolerance
      regression_count = 0
      if (present(regressors)) regression_count = size(regressors, 2)
      if (size(series) < 4 .or. ar_order < 0 .or. ma_order < 0 .or. &
         d < 0 .or. d > 2 .or. sd < 0 .or. sd > 2 .or. period < 1 .or. &
         limit < 1 .or. tol <= 0.0_dp .or. &
         .not. all(ieee_is_finite(series))) then
         out%info = 1
         return
      end if
      if (present(regressors)) then
         if (size(regressors, 1) /= size(series) .or. &
            .not. all(ieee_is_finite(regressors))) then
            out%info = 1
            return
         end if
      end if
      out%difference_polynomial = arima_difference_polynomial(d, period, sd)
      if (size(series) <= size(out%difference_polynomial) - 1 + &
         max(ar_order, ma_order + merge(period, 0, seasonal) + 1)) then
         out%info = 2
         return
      end if
      parameter_count = ar_order + ma_order + regression_count + &
         merge(1, 0, seasonal)
      allocate(initial(parameter_count))
      initial = 0.0_dp
      if (regression_count > 0) call initial_regression_estimate( &
         series, regressors, out%difference_polynomial, &
         initial(ar_order + ma_order + merge(1, 0, seasonal) + 1:))
      if (present(tuning_constant)) then
         if (tuning_constant <= 0.0_dp) then
            out%info = 1
            return
         end if
         tuning_values = [tuning_constant]
      else
         tuning_values = [1.0_dp, 0.8_dp, 0.64_dp, 1000.0_dp]
      end if
      if (parameter_count == 0) then
         allocate(parameters(0))
         optimized%parameters = parameters
         optimized%objective = objective(parameters, tuning_values(1))
         optimized%converged = .true.
         out%tuning_constant = tuning_values(1)
      else
         optimized%objective = huge(1.0_dp)
         do tuning_index = 1, size(tuning_values)
            candidate = bfgs_minimize_fd(local_objective, initial, limit, tol)
            if (candidate%objective < optimized%objective) then
               optimized = candidate
               out%tuning_constant = tuning_values(tuning_index)
            end if
         end do
      end if
      if (.not. allocated(optimized%parameters)) then
         out%info = 3
         return
      end if
      parameters = optimized%parameters
      call unpack_parameters(parameters, ar_order, ma_order, regression_count, &
         seasonal, out%ar, out%ma, out%seasonal_ma, &
         out%regression_coefficients)
      out%regression_residuals = series
      if (regression_count > 0) out%regression_residuals = &
         out%regression_residuals - matmul(regressors, &
         out%regression_coefficients)
      stationary = tfarima_difference(out%regression_residuals, &
         out%difference_polynomial)
      out%filter = robustarima_bounded_filter(stationary, out%ar, out%ma, &
         out%tuning_constant, seasonal_period=period, &
         seasonal_ma=merge(out%seasonal_ma, 0.0_dp, seasonal))
      if (out%filter%info /= 0) then
         out%info = 4
         return
      end if
      out%difference_order = d
      out%seasonal_period = period
      out%seasonal_difference_order = sd
      out%includes_seasonal_ma = seasonal
      allocate(out%parameter_covariance(parameter_count, parameter_count), &
         out%standard_error(parameter_count))
      out%parameter_covariance = 0.0_dp
      out%standard_error = 0.0_dp
      if (parameter_count > 0) then
         hessian = finite_difference_hessian(final_objective, parameters)
         call invert_matrix(hessian, inverse, status)
         if (status == 0) then
            out%parameter_covariance = 0.5_dp*(inverse + transpose(inverse))
            do i = 1, parameter_count
               out%standard_error(i) = sqrt(max(0.0_dp, &
                  out%parameter_covariance(i, i)))
            end do
         end if
      end if
      allocate(out%regression_covariance(regression_count, regression_count), &
         out%regression_standard_error(regression_count))
      out%regression_covariance = 0.0_dp
      out%regression_standard_error = 0.0_dp
      if (regression_count > 0) then
         inference = robustarima_tau_inference(out, regressors)
         if (inference%info == 0) then
            out%regression_covariance = inference%covariance
            out%regression_standard_error = inference%standard_error
            out%tau_inverse_efficiency = inference%inverse_efficiency
            regression_start = ar_order + ma_order + &
               merge(1, 0, seasonal) + 1
            out%parameter_covariance(regression_start:, regression_start:) = &
               inference%covariance
            out%standard_error(regression_start:) = inference%standard_error
         end if
      end if
      out%iterations = optimized%iterations
      out%converged = optimized%converged
      out%robust_aic = real(size(stationary), dp)* &
         log(max(out%filter%objective, tiny(1.0_dp))) + &
         2.0_dp*real(parameter_count, dp)
      if (.not. optimized%converged) out%info = 5

   contains

      pure real(dp) function local_objective(trial) result(value)
         !! Evaluate one bandwidth candidate during optimization.
         real(dp), intent(in) :: trial(:) !! Trial.
         value = objective(trial, tuning_values(tuning_index))
      end function local_objective

      pure real(dp) function final_objective(trial) result(value)
         !! Evaluate the selected-bandwidth objective for inference.
         real(dp), intent(in) :: trial(:) !! Trial.
         value = objective(trial, out%tuning_constant)
      end function final_objective

      pure real(dp) function objective(trial, tuning) result(value)
         !! Return the filtered tau objective at transformed parameters.
         real(dp), intent(in) :: trial(:) !! Trial.
         real(dp), intent(in) :: tuning !! Tuning.
         real(dp), allocatable :: ar(:), ma(:), beta(:), residual(:), work(:)
         real(dp) :: seasonal_coefficient
         type(robustarima_filter_t) :: filtered

         call unpack_parameters(trial, ar_order, ma_order, regression_count, &
            seasonal, ar, ma, seasonal_coefficient, beta)
         residual = series
         if (regression_count > 0) residual = residual - matmul(regressors, beta)
         work = tfarima_difference(residual, out%difference_polynomial)
         filtered = robustarima_bounded_filter(work, ar, ma, tuning, &
            seasonal_period=period, seasonal_ma=merge(seasonal_coefficient, &
            0.0_dp, seasonal))
         if (filtered%info /= 0 .or. .not. ieee_is_finite(filtered%objective)) then
            value = huge(1.0_dp)
         else
            value = log(max(filtered%objective, tiny(1.0_dp)))
         end if
      end function objective

   end function robustarima_fit

   pure function robustarima_tau_inference(fit, regressors, rho1_constant) &
      result(out)
      !! Compute robustarima's tau sandwich covariance for regression terms.
      type(robustarima_fit_t), intent(in) :: fit !! Previously fitted model.
      real(dp), intent(in) :: regressors(:, :) !! Regression design matrix.
      real(dp), intent(in), optional :: rho1_constant !! Rho1 constant.
      type(robustarima_tau_inference_t) :: out
      real(dp), allocatable :: standardized(:), working(:), information(:, :)
      real(dp), allocatable :: inverse(:, :)
      real(dp) :: xk, sum_rho, sum_psi2, sum_psi1, sum_square, sum_derivative
      real(dp) :: value, combined_psi, mean_weight
      integer :: first, observations, time, i, j, status

      xk = 0.405_dp
      if (present(rho1_constant)) xk = rho1_constant
      if (fit%filter%info /= 0 .or. xk <= 0.0_dp .or. &
         size(regressors, 1) /= size(fit%regression_residuals) .or. &
         size(regressors, 2) /= size(fit%regression_coefficients) .or. &
         size(regressors, 2) < 1 .or. &
         .not. all(ieee_is_finite(regressors))) then
         out%info = 1
         return
      end if
      out%filtered_design = filtered_regression_design(fit, regressors)
      if (size(out%filtered_design, 1) /= size(fit%filter%innovations)) then
         out%info = 2
         return
      end if
      first = fit%filter%initial_observations + 1
      observations = size(fit%filter%innovations) - first + 1
      allocate(standardized(size(fit%filter%innovations)))
      standardized = fit%filter%innovations/ &
         max(fit%filter%prediction_scale, tiny(1.0_dp))
      out%m_scale = robustarima_m_scale(standardized, first)
      if (out%m_scale <= tiny(1.0_dp)) then
         out%info = 3
         return
      end if
      working = standardized(first:)/out%m_scale
      sum_rho = sum(robustarima_rho(working))
      sum_psi2 = sum(robustarima_psi(working)*working)
      sum_psi1 = sum(robustarima_psi(working/xk)*(working/xk))
      if (abs(sum_psi1) <= tiny(1.0_dp)) then
         out%info = 3
         return
      end if
      out%rho1_weight = (2.0_dp*sum_rho - sum_psi2)/sum_psi1
      allocate(out%implicit_weight(observations))
      sum_square = 0.0_dp
      sum_derivative = 0.0_dp
      do time = 1, observations
         value = working(time)
         combined_psi = out%rho1_weight*robustarima_psi(value/xk)/xk + &
            robustarima_psi(value)
         sum_square = sum_square + combined_psi**2
         sum_derivative = sum_derivative + out%rho1_weight* &
            robustarima_psi_derivative(value/xk)/xk**2 + &
            robustarima_psi_derivative(value)
         if (abs(value) > sqrt(epsilon(1.0_dp))) then
            out%implicit_weight(time) = combined_psi/value
         else
            out%implicit_weight(time) = out%rho1_weight/xk**2 + 1.0_dp
         end if
      end do
      mean_weight = sum(out%implicit_weight)/real(observations, dp)
      if (abs(sum_derivative) <= tiny(1.0_dp) .or. &
         abs(mean_weight) <= tiny(1.0_dp)) then
         out%info = 3
         return
      end if
      out%inverse_efficiency = real(observations, dp)*sum_square/ &
         sum_derivative**2
      allocate(information(size(regressors, 2), size(regressors, 2)))
      information = 0.0_dp
      do i = 1, size(regressors, 2)
         do j = 1, size(regressors, 2)
            information(i, j) = sum(out%filtered_design(first:, i)* &
               out%filtered_design(first:, j)*out%implicit_weight/ &
               fit%filter%prediction_scale(first:)**2)/mean_weight
         end do
      end do
      call invert_matrix(information, inverse, status)
      if (status /= 0) then
         out%info = 4
         return
      end if
      out%covariance = out%m_scale**2*out%inverse_efficiency*inverse
      out%covariance = 0.5_dp*(out%covariance + transpose(out%covariance))
      allocate(out%standard_error(size(regressors, 2)))
      do i = 1, size(regressors, 2)
         out%standard_error(i) = sqrt(max(0.0_dp, out%covariance(i, i)))
      end do
   end function robustarima_tau_inference

   pure function filtered_regression_design(fit, regressors) result(design)
      !! Filter differenced regressors with the fitted robust gains and weights.
      type(robustarima_fit_t), intent(in) :: fit !! Previously fitted model.
      real(dp), intent(in) :: regressors(:, :) !! Regression design matrix.
      real(dp), allocatable :: design(:, :)
      real(dp), allocatable :: combined_ma(:), transition(:, :)
      real(dp), allocatable :: disturbance(:), state(:), predicted_state(:)
      real(dp), allocatable :: predicted_covariance(:, :)
      real(dp), allocatable :: previous_covariance(:, :), gain_numerator(:)
      real(dp), allocatable :: filtered_regressor(:)
      real(dp) :: innovation, innovation_variance
      integer :: dimension, column, time, i, p, q

      if (size(regressors, 1) /= size(fit%regression_residuals) .or. &
         size(regressors, 2) < 1 .or. fit%filter%info /= 0) then
         allocate(design(0, 0))
         return
      end if
      dimension = size(fit%filter%state, 1)
      if (dimension < 1 .or. &
         size(fit%filter%state, 2) /= size(fit%filter%innovations) .or. &
         size(fit%filter%state_covariance, 1) /= dimension .or. &
         size(fit%filter%state_covariance, 2) /= dimension .or. &
         size(fit%filter%state_covariance, 3) /= &
         size(fit%filter%innovations)) then
         allocate(design(0, 0))
         return
      end if
      combined_ma = seasonal_ma_coefficients(fit%ma, &
         fit%seasonal_period, merge(fit%seasonal_ma, 0.0_dp, &
         fit%includes_seasonal_ma))
      p = size(fit%ar)
      q = size(combined_ma)
      allocate(design(size(fit%filter%innovations), size(regressors, 2)), &
         transition(dimension, dimension), disturbance(dimension), &
         state(dimension), predicted_state(dimension), &
         predicted_covariance(dimension, dimension), &
         previous_covariance(dimension, dimension), &
         gain_numerator(dimension))
      design = 0.0_dp
      transition = 0.0_dp
      if (p > 0) transition(:p, 1) = fit%ar
      do i = 1, dimension - 1
         transition(i, i + 1) = 1.0_dp
      end do
      disturbance = 0.0_dp
      disturbance(1) = 1.0_dp
      if (q > 0) disturbance(2:q + 1) = -combined_ma
      do column = 1, size(regressors, 2)
         filtered_regressor = tfarima_difference(regressors(:, column), &
            fit%difference_polynomial)
         if (size(filtered_regressor) /= size(design, 1)) then
            deallocate(design)
            allocate(design(0, 0))
            return
         end if
         state = 0.0_dp
         previous_covariance = 0.0_dp
         do time = 1, size(design, 1)
            if (time > 1) previous_covariance = &
               fit%filter%state_covariance(:, :, time - 1)
            predicted_state = matmul(transition, state)
            predicted_covariance = matmul(matmul(transition, &
               previous_covariance), transpose(transition)) + &
               fit%filter%initial_scale**2* &
               spread(disturbance, 2, dimension)* &
               spread(disturbance, 1, dimension)
            gain_numerator = predicted_covariance(:, 1)
            innovation_variance = max( &
               fit%filter%prediction_scale(time)**2, tiny(1.0_dp))
            innovation = filtered_regressor(time) - predicted_state(1)
            design(time, column) = innovation
            state = predicted_state + fit%filter%weights(time)* &
               gain_numerator*innovation/innovation_variance
         end do
      end do
   end function filtered_regression_design

   pure function robustarima_select_ar_order(series, maximum_order, &
      difference_order, regressors, seasonal_period, &
      seasonal_difference_order, include_seasonal_ma, tuning_constant, &
      max_iterations, tolerance) result(out)
      !! Select an AR order by the package-style robust Akaike criterion.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      integer, intent(in) :: maximum_order !! Maximum order.
      integer, intent(in), optional :: difference_order !! Difference order.
      integer, intent(in), optional :: seasonal_period !! Seasonal period.
      integer, intent(in), optional :: seasonal_difference_order !! Seasonal difference order.
      real(dp), intent(in), optional :: regressors(:, :) !! Regression design matrix.
      real(dp), intent(in), optional :: tuning_constant !! Tuning constant.
      logical, intent(in), optional :: include_seasonal_ma !! Whether to include the seasonal moving-average.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      type(robustarima_order_selection_t) :: out
      type(robustarima_fit_t) :: candidate
      integer :: order

      if (maximum_order < 0) then
         out%info = 1
         return
      end if
      allocate(out%criterion(0:maximum_order))
      out%criterion = huge(1.0_dp)
      do order = 0, maximum_order
         candidate = robustarima_fit(series, order, 0, difference_order, &
            regressors, seasonal_period, seasonal_difference_order, &
            include_seasonal_ma, tuning_constant, max_iterations, tolerance)
         if (candidate%info /= 0) cycle
         out%criterion(order) = candidate%robust_aic
         if (out%criterion(order) < out%criterion(out%selected_order)) then
            out%selected_order = order
            out%fit = candidate
         else if (order == 0) then
            out%fit = candidate
         end if
      end do
      if (.not. ieee_is_finite(out%criterion(out%selected_order)) .or. &
         out%criterion(out%selected_order) >= huge(1.0_dp)) out%info = 2
   end function robustarima_select_ar_order

   pure function robustarima_forecast(fit, series, horizon, regressors, &
      future_regressors, levels) result(out)
      !! Forecast from a robust fit using its cleaned innovation history.
      type(robustarima_fit_t), intent(in) :: fit !! Previously fitted model.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      real(dp), intent(in), optional :: regressors(:, :) !! Regression design matrix.
      real(dp), intent(in), optional :: future_regressors(:, :) !! Future regressors.
      real(dp), intent(in), optional :: levels(:) !! Levels.
      type(robustarima_forecast_t) :: out
      type(tfarima_forecast_t) :: base
      real(dp), allocatable :: residual(:), ar_polynomial(:), ma_polynomial(:)
      real(dp), allocatable :: combined_ma(:), future_signal(:), stationary(:)

      if (fit%info /= 0 .or. horizon < 1 .or. size(series) < 1) then
         out%info = 1
         return
      end if
      residual = series
      if (size(fit%regression_coefficients) > 0) then
         if (.not. present(regressors) .or. .not. present(future_regressors)) then
            out%info = 2
            return
         end if
         if (size(regressors, 1) /= size(series) .or. &
            size(regressors, 2) /= size(fit%regression_coefficients) .or. &
            size(future_regressors, 1) /= horizon .or. &
            size(future_regressors, 2) /= size(fit%regression_coefficients)) then
            out%info = 2
            return
         end if
         residual = residual - matmul(regressors, fit%regression_coefficients)
         future_signal = matmul(future_regressors, fit%regression_coefficients)
      else
         allocate(future_signal(horizon))
         future_signal = 0.0_dp
      end if
      stationary = tfarima_difference(residual, fit%difference_polynomial)
      if (size(stationary) == size(fit%filter%cleaned)) residual = &
         integrate_filtered_series(residual, fit%filter%cleaned, &
         fit%difference_polynomial)
      ar_polynomial = [1.0_dp, -fit%ar]
      ar_polynomial = tfarima_polynomial_multiply(ar_polynomial, &
         fit%difference_polynomial)
      combined_ma = seasonal_ma_coefficients(fit%ma, fit%seasonal_period, &
         merge(fit%seasonal_ma, 0.0_dp, fit%includes_seasonal_ma))
      ma_polynomial = [1.0_dp, -combined_ma]
      if (present(levels)) then
         base = tfarima_arima_forecast(residual, ar_polynomial, [1.0_dp], &
            ma_polynomial, fit%filter%innovation_scale**2, horizon, &
            levels=levels)
      else
         base = tfarima_arima_forecast(residual, ar_polynomial, [1.0_dp], &
            ma_polynomial, fit%filter%innovation_scale**2, horizon)
      end if
      if (base%info /= 0) then
         out%info = 3
         return
      end if
      out%mean = base%mean + future_signal
      out%standard_error = sqrt(max(0.0_dp, base%variance))
      out%lower = base%lower + spread(future_signal, 2, size(base%lower, 2))
      out%upper = base%upper + spread(future_signal, 2, size(base%upper, 2))
   end function robustarima_forecast

   pure function robustarima_detect_outliers(fit, series, critical_value, &
      include_innovation_outliers) result(out)
      !! Detect IO, AO, and LS effects from robust filtered innovations.
      type(robustarima_fit_t), intent(in) :: fit !! Previously fitted model.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      real(dp), intent(in), optional :: critical_value !! Critical value.
      logical, intent(in), optional :: include_innovation_outliers !! Whether to include the innovation outliers.
      type(robustarima_outliers_t) :: out
      real(dp), allocatable :: aligned(:), aligned_scale(:)
      real(dp), allocatable :: response(:, :), candidate(:)
      real(dp), allocatable :: selected_effect(:), selected_statistic(:)
      real(dp), allocatable :: impulse(:)
      real(dp) :: cutoff, effect, statistic, denominator, best_statistic
      real(dp) :: best_effect
      integer, allocatable :: selected_type(:)
      integer :: offset, position, outlier_type, best_type, count, i
      logical :: include_io

      cutoff = 3.0_dp
      if (size(series) > 200) cutoff = 3.5_dp
      if (size(series) > 500) cutoff = 4.0_dp
      if (present(critical_value)) cutoff = critical_value
      include_io = .false.
      if (present(include_innovation_outliers)) include_io = &
         include_innovation_outliers
      if (fit%info /= 0 .or. size(series) < 1 .or. cutoff <= 0.0_dp) then
         out%info = 1
         return
      end if
      offset = size(fit%difference_polynomial) - 1
      allocate(aligned(size(series)), aligned_scale(size(series)), &
         response(size(series), 3), &
         candidate(size(series)), selected_type(size(series)), &
         selected_effect(size(series)), selected_statistic(size(series)))
      aligned = 0.0_dp
      aligned(offset + 1:) = fit%filter%innovations
      aligned_scale = max(fit%filter%innovation_scale, tiny(1.0_dp))
      aligned_scale(offset + 1:) = fit%filter%prediction_scale* &
         fit%filter%innovation_scale/max(fit%filter%initial_scale, &
         tiny(1.0_dp))
      response(:, robustarima_outlier_io) = tfarima_outlier_response( &
         robustarima_outlier_io, size(series), [1.0_dp, -fit%ar], &
         fit%difference_polynomial, [1.0_dp, &
         -seasonal_ma_coefficients(fit%ma, fit%seasonal_period, &
         merge(fit%seasonal_ma, 0.0_dp, fit%includes_seasonal_ma))])
      response(:, robustarima_outlier_ao) = tfarima_outlier_response( &
         robustarima_outlier_ao, size(series), [1.0_dp, -fit%ar], &
         fit%difference_polynomial, [1.0_dp, &
         -seasonal_ma_coefficients(fit%ma, fit%seasonal_period, &
         merge(fit%seasonal_ma, 0.0_dp, fit%includes_seasonal_ma))])
      response(:, robustarima_outlier_ao) = 0.0_dp
      response(1, robustarima_outlier_ao) = 1.0_dp
      response(:, robustarima_outlier_ls) = tfarima_outlier_response( &
         robustarima_outlier_ls, size(series), [1.0_dp, -fit%ar], &
         fit%difference_polynomial, [1.0_dp, &
         -seasonal_ma_coefficients(fit%ma, fit%seasonal_period, &
         merge(fit%seasonal_ma, 0.0_dp, fit%includes_seasonal_ma))])
      count = 0
      candidate = 0.0_dp
      do position = offset + fit%filter%initial_observations + 1, size(series)
         best_statistic = 0.0_dp
         best_effect = 0.0_dp
         best_type = 0
         do outlier_type = merge(robustarima_outlier_io, &
            robustarima_outlier_ao, include_io), robustarima_outlier_ls
            denominator = sum((response(:size(series) - position + 1, &
               outlier_type)/aligned_scale(position:))**2)
            if (denominator <= tiny(1.0_dp)) cycle
            effect = sum(aligned(position:)* &
               response(:size(series) - position + 1, outlier_type)/ &
               aligned_scale(position:)**2)/denominator
            statistic = effect*sqrt(denominator)
            if (abs(statistic) > abs(best_statistic)) then
               best_statistic = statistic
               best_effect = effect
               best_type = outlier_type
            end if
         end do
         if (abs(best_statistic) >= cutoff) then
            count = count + 1
            candidate(count) = real(position, dp)
            aligned(position:) = aligned(position:) - &
               best_effect*response(:size(series) - position + 1, best_type)
            selected_type(count) = best_type
            selected_effect(count) = best_effect
            selected_statistic(count) = best_statistic
         end if
      end do
      allocate(out%position(count), out%outlier_type(count), &
         out%effect(count), out%t_statistic(count), out%cleaned_series(size(series)))
      out%cleaned_series = series
      do i = 1, count
         position = nint(candidate(i))
         out%position(i) = position
         out%outlier_type(i) = selected_type(i)
         out%effect(i) = selected_effect(i)
         out%t_statistic(i) = selected_statistic(i)
         select case (out%outlier_type(i))
         case (robustarima_outlier_ao)
            out%cleaned_series(position) = out%cleaned_series(position) - &
               out%effect(i)
         case (robustarima_outlier_ls)
            out%cleaned_series(position:) = out%cleaned_series(position:) - &
               out%effect(i)
         case (robustarima_outlier_io)
            impulse = tfarima_polynomial_ratio([1.0_dp, &
               -seasonal_ma_coefficients(fit%ma, fit%seasonal_period, &
               merge(fit%seasonal_ma, 0.0_dp, fit%includes_seasonal_ma))], &
               tfarima_polynomial_multiply([1.0_dp, -fit%ar], &
               fit%difference_polynomial), size(series) - position)
            out%cleaned_series(position:) = out%cleaned_series(position:) - &
               out%effect(i)*impulse
         end select
      end do
      out%scale_before = fit%filter%innovation_scale
      out%scale_after = fit%filter%initial_scale*robustarima_m_scale( &
         aligned/aligned_scale, offset + fit%filter%initial_observations + 1)
   end function robustarima_detect_outliers

   pure function arima_difference_polynomial(order, seasonal_period, &
      seasonal_order) result(polynomial)
      !! Construct regular and seasonal differencing operators.
      integer, intent(in) :: order !! Model or polynomial order.
      integer, intent(in) :: seasonal_period !! Seasonal period.
      integer, intent(in) :: seasonal_order !! Seasonal order.
      real(dp), allocatable :: polynomial(:), seasonal(:)

      polynomial = tfarima_polynomial_power([1.0_dp, -1.0_dp], order)
      allocate(seasonal(seasonal_period + 1))
      seasonal = 0.0_dp
      seasonal(1) = 1.0_dp
      seasonal(seasonal_period + 1) = -1.0_dp
      seasonal = tfarima_polynomial_power(seasonal, seasonal_order)
      polynomial = tfarima_polynomial_multiply(polynomial, seasonal)
   end function arima_difference_polynomial

   pure function seasonal_ma_coefficients(ma, seasonal_period, seasonal_ma) &
      result(combined)
      !! Expand regular and single seasonal MA factors without the leading one.
      real(dp), intent(in) :: ma(:) !! Moving-average coefficients.
      integer, intent(in), optional :: seasonal_period !! Seasonal period.
      real(dp), intent(in), optional :: seasonal_ma !! Seasonal moving-average.
      real(dp), allocatable :: combined(:), regular(:), seasonal(:), product(:)
      integer :: period

      period = 1
      if (present(seasonal_period)) period = seasonal_period
      if (period < 1) then
         allocate(combined(0))
         return
      end if
      regular = [1.0_dp, -ma]
      if (present(seasonal_ma)) then
         if (abs(seasonal_ma) > 0.0_dp) then
            allocate(seasonal(period + 1))
            seasonal = 0.0_dp
            seasonal(1) = 1.0_dp
            seasonal(period + 1) = -seasonal_ma
            product = tfarima_polynomial_multiply(regular, seasonal)
            combined = -product(2:)
            return
         end if
      end if
      combined = ma
   end function seasonal_ma_coefficients

   pure function integrate_filtered_series(original, stationary, &
      difference_polynomial) result(cleaned)
      !! Reconstruct a cleaned original-scale path from filtered differences.
      real(dp), intent(in) :: original(:) !! Original.
      real(dp), intent(in) :: stationary(:) !! Stationary.
      real(dp), intent(in) :: difference_polynomial(:) !! Difference polynomial coefficients.
      real(dp), allocatable :: cleaned(:)
      real(dp) :: value
      integer :: degree, time, lag

      degree = size(difference_polynomial) - 1
      allocate(cleaned(size(original)))
      cleaned = original
      if (size(stationary) /= size(original) - degree .or. &
         abs(difference_polynomial(1)) <= tiny(1.0_dp)) return
      do time = degree + 1, size(original)
         value = stationary(time - degree)
         do lag = 1, degree
            value = value - difference_polynomial(lag + 1)* &
               cleaned(time - lag)
         end do
         cleaned(time) = value/difference_polynomial(1)
      end do
   end function integrate_filtered_series

   pure subroutine unpack_parameters(parameters, ar_order, ma_order, &
      regression_count, seasonal, ar, ma, seasonal_ma, beta)
      !! Transform unconstrained coordinates into model coefficients.
      real(dp), intent(in) :: parameters(:) !! Model parameter values.
      integer, intent(in) :: ar_order !! Autoregressive order.
      integer, intent(in) :: ma_order !! Moving-average order.
      integer, intent(in) :: regression_count !! Number of regression.
      logical, intent(in) :: seasonal !! Flag controlling seasonal.
      real(dp), allocatable, intent(out) :: ar(:) !! Autoregressive coefficients.
      real(dp), allocatable, intent(out) :: ma(:) !! Moving-average coefficients.
      real(dp), allocatable, intent(out) :: beta(:) !! Regression or model coefficients.
      real(dp), intent(out) :: seasonal_ma !! Seasonal moving-average.
      real(dp), allocatable :: partial(:)
      integer :: first

      allocate(ar(ar_order), ma(ma_order), beta(regression_count))
      if (ar_order > 0) then
         partial = 2.0_dp*atan(parameters(:ar_order))/acos(-1.0_dp)
         ar = partial_to_coefficients(partial)
      end if
      if (ma_order > 0) then
         partial = 2.0_dp*atan(parameters(ar_order + 1: &
            ar_order + ma_order))/acos(-1.0_dp)
         ma = partial_to_coefficients(partial)
      end if
      first = ar_order + ma_order + 1
      seasonal_ma = 0.0_dp
      if (seasonal) then
         seasonal_ma = 2.0_dp*atan(parameters(first))/acos(-1.0_dp)
         first = first + 1
      end if
      if (regression_count > 0) beta = parameters(first:)
   end subroutine unpack_parameters

   pure function partial_to_coefficients(partial) result(coefficients)
      !! Apply the inverse Durbin transform to partial autocorrelations.
      real(dp), intent(in) :: partial(:) !! Partial.
      real(dp), allocatable :: coefficients(:), previous(:)
      integer :: order, lag

      allocate(coefficients(size(partial)))
      coefficients = 0.0_dp
      if (size(partial) < 1) return
      coefficients(1) = partial(1)
      do order = 2, size(partial)
         previous = coefficients
         coefficients(order) = partial(order)
         do lag = 1, order - 1
            coefficients(lag) = previous(lag) - &
               partial(order)*previous(order - lag)
         end do
      end do
   end function partial_to_coefficients

   pure subroutine initial_regression_estimate(series, regressors, &
      difference_polynomial, beta)
      !! Initialize regression coefficients by differenced least squares.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      real(dp), intent(in) :: regressors(:, :) !! Regression design matrix.
      real(dp), intent(in) :: difference_polynomial(:) !! Difference polynomial coefficients.
      real(dp), intent(out) :: beta(:) !! Regression or model coefficients.
      real(dp), allocatable :: y(:), x(:, :), inverse(:, :)
      integer :: column, status

      y = tfarima_difference(series, difference_polynomial)
      allocate(x(size(y), size(regressors, 2)))
      do column = 1, size(regressors, 2)
         x(:, column) = tfarima_difference(regressors(:, column), &
            difference_polynomial)
      end do
      call invert_matrix(matmul(transpose(x), x), inverse, status)
      if (status == 0) then
         beta = matmul(inverse, matmul(transpose(x), y))
      else
         beta = 0.0_dp
      end if
   end subroutine initial_regression_estimate

end module robustarima_mod
