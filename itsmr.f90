! SPDX-License-Identifier: BSD-2-Clause
! SPDX-FileComment: Algorithms translated from the FreeBSD-licensed R itsmr package.
! Numerical translations of distinct algorithms from the FreeBSD-licensed ITSMR package.
module itsmr_mod
   use kind_mod, only: dp
   use special_functions_mod, only: regularized_gamma_q
   use fourier_mod, only: real_dft, inverse_real_dft
   use polynomial_mod, only: polynomial_product
   use linalg_mod, only: invert_matrix
   use optimization_mod, only: optimization_result_t, bfgs_minimize_fd, &
      finite_difference_hessian
   use time_series_stats_mod, only: burg_result_t, burg_fit, harmonic_regression_result_t, &
      harmonic_regression, yule_walker_result_t, yule_walker_fit
   use stats_mod, only: ols_fit
   use utils_mod, only: inverse_standard_normal
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   private

   integer, parameter, public :: itsmr_transform_difference = 1
   integer, parameter, public :: itsmr_transform_harmonic = 2
   integer, parameter, public :: itsmr_transform_trend = 3
   integer, parameter, public :: itsmr_transform_seasonal = 4
   integer, parameter, public :: itsmr_transform_log = 5

   type, public :: itsmr_arma_model_t
      ! ARMA coefficients, innovation variance, corrected AIC, and standard errors.
      real(dp), allocatable :: ar(:), ma(:), ar_standard_error(:), ma_standard_error(:)
      real(dp) :: innovation_variance = 0.0_dp
      real(dp) :: aicc = 0.0_dp
      integer :: iterations = 0
      integer :: info = 0
      logical :: converged = .false.
   end type

   type, public :: itsmr_innovations_t
      ! Innovations-algorithm coefficients, prediction variances, and predictions.
      real(dp), allocatable :: theta(:, :), variance(:), prediction(:), innovations(:)
      integer :: info = 0
   end type

   type, public :: itsmr_forecast_t
      ! ARMA forecasts, uncertainty, recovered innovations, and MA-infinity weights.
      real(dp), allocatable :: mean(:), standard_error(:), lower(:), upper(:)
      real(dp), allocatable :: innovations(:), psi(:)
      real(dp) :: confidence_level = 0.95_dp
      integer :: info = 0
   end type

   type, public :: itsmr_arar_t
      ! ARAR forecasts, selected sparse model, filters, uncertainty, and status.
      real(dp), allocatable :: mean(:), standard_error(:), lower(:), upper(:)
      real(dp), allocatable :: coefficients(:), memory_polynomial(:), filter(:)
      integer, allocatable :: lags(:)
      real(dp) :: innovation_variance = 0.0_dp
      real(dp) :: shortened_mean = 0.0_dp
      real(dp) :: confidence_level = 0.95_dp
      integer :: info = 0
   end type

   type, public :: itsmr_residuals_t
      ! Fitted values, exact innovations, prediction variances, and standardized innovations.
      real(dp), allocatable :: fitted(:), innovations(:), variance(:), standardized(:)
      real(dp) :: mean = 0.0_dp
      integer :: info = 0
   end type

   type, public :: itsmr_transform_t
      ! One typed preprocessing step for transformed ARMA forecasting.
      integer :: kind = 0
      integer :: lag = 1
      integer :: order = 1
      real(dp), allocatable :: periods(:)
   end type

   type, public :: itsmr_transformed_forecast_t
      ! Forecasts and intervals after applying and reversing typed transforms.
      real(dp), allocatable :: mean(:), standard_error(:), lower(:), upper(:), transformed_series(:)
      real(dp) :: confidence_level = 0.95_dp
      integer :: info = 0
   end type

   type, public :: itsmr_rank_filter_t
      ! Rank-filtered reconstruction and retained positive-frequency diagnostics.
      real(dp), allocatable :: filtered(:), frequency(:), amplitude(:)
      complex(dp), allocatable :: coefficients(:)
      integer, allocatable :: retained_bins(:)
      integer :: info = 0
   end type

   type, public :: itsmr_randomness_tests_t
      ! ITSMR residual randomness statistics and upper-tail p-values.
      real(dp) :: ljung_box = 0.0_dp
      real(dp) :: ljung_box_p_value = 1.0_dp
      real(dp) :: mcleod_li = 0.0_dp
      real(dp) :: mcleod_li_p_value = 1.0_dp
      real(dp) :: turning_points = 0.0_dp
      real(dp) :: turning_point_z = 0.0_dp
      real(dp) :: turning_point_p_value = 1.0_dp
      real(dp) :: difference_signs = 0.0_dp
      real(dp) :: difference_sign_z = 0.0_dp
      real(dp) :: difference_sign_p_value = 1.0_dp
      real(dp) :: rank_statistic = 0.0_dp
      real(dp) :: rank_z = 0.0_dp
      real(dp) :: rank_p_value = 1.0_dp
      integer :: lag = 0
      integer :: info = 0
   end type

   type :: itsmr_transform_state_t
      ! Internal fitted state needed to reverse one transform.
      integer :: kind = 0
      integer :: lag = 0
      real(dp), allocatable :: history(:), restoration(:)
   end type

   public :: arma_acvf, innovations_algorithm, innovations_ma_fit, hannan_rissanen_fit
   public :: arma_mle_fit, arma_autofit
   public :: arma_infinite_ma, arma_forecast
   public :: arar_forecast
   public :: burg_ar_fit
   public :: arma_infinite_ar, arma_residuals
   public :: harmonic_regression_result_t, harmonic_regression
   public :: transformed_arma_forecast
   public :: spectral_rank_filter
   public :: residual_randomness_tests
   public :: regularized_gamma_q
   public :: itsmr_moving_average, itsmr_exponential_smooth, itsmr_fft_smooth
   public :: itsmr_seasonal_component

contains

   pure function arma_acvf(model, max_lag) result(covariance)
      !! Return theoretical autocovariances of a causal ARMA model.
      type(itsmr_arma_model_t), intent(in) :: model !! Model specification.
      integer, intent(in) :: max_lag !! Maximum lag to consider.
      real(dp), allocatable :: covariance(:)
      real(dp), allocatable :: transition(:, :), loading(:), state_covariance(:, :), next_covariance(:, :), vector(:)
      real(dp) :: convergence
      integer :: p, q, state_size, i, lag

      if (max_lag < 0 .or. model%innovation_variance < 0.0_dp) then
         allocate(covariance(0))
         return
      end if
      p = 0
      if (allocated(model%ar)) p = size(model%ar)
      q = 0
      if (allocated(model%ma)) q = size(model%ma)
      state_size = max(1, max(p, q + 1))
      allocate(transition(state_size, state_size), loading(state_size))
      allocate(state_covariance(state_size, state_size), next_covariance(state_size, state_size))
      transition = 0.0_dp
      loading = 0.0_dp
      if (p > 0) transition(:p, 1) = model%ar
      do i = 1, state_size - 1
         transition(i, i + 1) = 1.0_dp
      end do
      loading(1) = 1.0_dp
      if (q > 0) loading(2:q + 1) = model%ma
      state_covariance = 0.0_dp
      convergence = huge(1.0_dp)
      do i = 1, 100000
         next_covariance = matmul(matmul(transition, state_covariance), transpose(transition)) + &
            model%innovation_variance*spread(loading, 2, state_size)*spread(loading, 1, state_size)
         convergence = maxval(abs(next_covariance - state_covariance))
         state_covariance = next_covariance
         if (convergence <= 100.0_dp*epsilon(1.0_dp)) exit
      end do
      allocate(covariance(0:max_lag), vector(state_size))
      if (convergence > 1.0e-8_dp) then
         covariance = huge(1.0_dp)
         return
      end if
      vector = 0.0_dp
      vector(1) = 1.0_dp
      do lag = 0, max_lag
         covariance(lag) = dot_product(state_covariance(1, :), vector)
         vector = matmul(transpose(transition), vector)
      end do
   end function arma_acvf

   pure function innovations_algorithm(covariance, observations) result(out)
      !! Run the innovations recursion for a stationary covariance sequence.
      real(dp), intent(in) :: covariance(0:) !! Covariance matrix.
      real(dp), intent(in), optional :: observations(:) !! Observed time-series values.
      type(itsmr_innovations_t) :: out
      integer :: levels, n, k, j
      real(dp) :: accumulated

      levels = ubound(covariance, 1)
      if (levels < 0 .or. covariance(0) <= 0.0_dp) then
         out%info = 1
         return
      end if
      allocate(out%theta(levels, levels), out%variance(0:levels))
      out%theta = 0.0_dp
      out%variance(0) = covariance(0)
      do n = 1, levels
         do k = 0, n - 1
            accumulated = 0.0_dp
            do j = 0, k - 1
               accumulated = accumulated + out%theta(k, k - j)*out%theta(n, n - j)*out%variance(j)
            end do
            out%theta(n, n - k) = (covariance(n - k) - accumulated)/out%variance(k)
         end do
         accumulated = 0.0_dp
         do j = 0, n - 1
            accumulated = accumulated + out%theta(n, n - j)**2*out%variance(j)
         end do
         out%variance(n) = covariance(0) - accumulated
         if (out%variance(n) <= tiny(1.0_dp)) then
            out%info = 2
            return
         end if
      end do
      if (present(observations)) then
         if (size(observations) > levels + 1) then
            out%info = 3
            return
         end if
         allocate(out%prediction(size(observations)), out%innovations(size(observations)))
         out%prediction = 0.0_dp
         out%innovations(1) = observations(1)
         do n = 1, size(observations) - 1
            do k = 1, n
               out%prediction(n + 1) = out%prediction(n + 1) + &
                  out%theta(n, k)*out%innovations(n + 1 - k)
            end do
            out%innovations(n + 1) = observations(n + 1) - out%prediction(n + 1)
         end do
      end if
   end function innovations_algorithm

   pure function innovations_ma_fit(series, order, recursion_level) result(model)
      !! Estimate an MA model using ITSMR's innovations algorithm.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      integer, intent(in) :: order !! Model or polynomial order.
      integer, intent(in), optional :: recursion_level !! Recursion level.
      type(itsmr_arma_model_t) :: model
      type(itsmr_innovations_t) :: recursion, fitted
      real(dp), allocatable :: centered(:), sample_covariance(:), model_covariance(:)
      real(dp) :: mean_value, log_likelihood
      integer :: level, n, lag, parameter_count

      n = size(series)
      level = 17
      if (present(recursion_level)) level = recursion_level
      if (order < 1 .or. level < order .or. level + 1 >= n) then
         model%info = 1
         return
      end if
      mean_value = sum(series)/real(n, dp)
      centered = series - mean_value
      allocate(sample_covariance(0:level + 1))
      do lag = 0, level + 1
         sample_covariance(lag) = dot_product(centered(:n - lag), centered(lag + 1:))/real(n, dp)
      end do
      recursion = innovations_algorithm(sample_covariance)
      if (recursion%info /= 0) then
         model%info = 10 + recursion%info
         return
      end if
      allocate(model%ar(0), model%ma(order), model%ar_standard_error(0), model%ma_standard_error(order))
      model%ma = recursion%theta(level, 1:order)
      model%ma_standard_error(1) = sqrt(1.0_dp/real(n, dp))
      do lag = 2, order
         model%ma_standard_error(lag) = sqrt((1.0_dp + sum(recursion%theta(level, 1:lag - 1)**2))/real(n, dp))
      end do
      model%innovation_variance = 1.0_dp
      model_covariance = arma_acvf(model, n - 1)
      fitted = innovations_algorithm(model_covariance, centered)
      if (fitted%info /= 0) then
         model%info = 20 + fitted%info
         return
      end if
      model%innovation_variance = sum(fitted%innovations**2/fitted%variance(:n - 1))/real(n, dp)
      parameter_count = order + 1
      log_likelihood = -0.5_dp*real(n, dp)*log(2.0_dp*acos(-1.0_dp)*model%innovation_variance) - &
         0.5_dp*sum(log(fitted%variance(:n - 1))) - 0.5_dp*real(n, dp)
      model%aicc = -2.0_dp*log_likelihood + &
         2.0_dp*real(parameter_count*n, dp)/real(n - order - 2, dp)
   end function innovations_ma_fit

   pure function hannan_rissanen_fit(series, ar_order, ma_order) result(model)
      !! Estimate an ARMA model by ITSMR's Hannan-Rissanen regression.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      integer, intent(in) :: ar_order !! Autoregressive order.
      integer, intent(in) :: ma_order !! Moving-average order.
      type(itsmr_arma_model_t) :: model
      type(yule_walker_result_t) :: preliminary
      real(dp), allocatable :: centered(:), proxy(:), design(:, :), response(:)
      real(dp), allocatable :: coefficients(:), standard_errors(:), residuals(:)
      real(dp) :: mean_value, rss, scale
      integer :: n, preliminary_order, offset, rows, columns, row, t, lag, status

      n = size(series)
      preliminary_order = 20 + ar_order + ma_order
      offset = max(ar_order, ma_order)
      rows = n - preliminary_order - offset
      columns = ar_order + ma_order
      if (ar_order < 0 .or. ma_order < 1 .or. rows <= columns .or. preliminary_order >= n) then
         model%info = 1
         return
      end if
      mean_value = sum(series)/real(n, dp)
      centered = series - mean_value
      preliminary = yule_walker_fit(centered, preliminary_order)
      if (preliminary%info /= 0) then
         model%info = 10 + preliminary%info
         return
      end if
      allocate(proxy(n))
      proxy = 0.0_dp
      do t = preliminary_order + 1, n
         proxy(t) = centered(t)
         do lag = 1, preliminary_order
            proxy(t) = proxy(t) - preliminary%coefficients(lag)*centered(t - lag)
         end do
      end do
      allocate(design(rows, columns), response(rows))
      do row = 1, rows
         t = preliminary_order + offset + row
         response(row) = centered(t)
         do lag = 1, ar_order
            design(row, lag) = centered(t - lag)
         end do
         do lag = 1, ma_order
            design(row, ar_order + lag) = proxy(t - lag)
         end do
      end do
      call ols_fit(design, response, coefficients, standard_errors, residuals, rss, status)
      if (status /= 0) then
         model%info = 20 + status
         return
      end if
      scale = sqrt(real(rows - columns, dp)/real(rows, dp))
      standard_errors = scale*standard_errors
      allocate(model%ar(ar_order), model%ma(ma_order))
      allocate(model%ar_standard_error(ar_order), model%ma_standard_error(ma_order))
      if (ar_order > 0) then
         model%ar = coefficients(:ar_order)
         model%ar_standard_error = standard_errors(:ar_order)
      end if
      model%ma = coefficients(ar_order + 1:)
      model%ma_standard_error = standard_errors(ar_order + 1:)
      call update_innovation_statistics(centered, model)
   end function hannan_rissanen_fit

   pure function arma_mle_fit(series, ar_order, ma_order, max_iterations, tolerance) result(model)
      !! Fit a causal and invertible ARMA model by exact innovations likelihood.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      integer, intent(in) :: ar_order !! Autoregressive order.
      integer, intent(in) :: ma_order !! Moving-average order.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      type(itsmr_arma_model_t) :: model
      type(itsmr_arma_model_t) :: initial_model
      type(yule_walker_result_t) :: ar_initial
      type(optimization_result_t) :: optimization
      real(dp), allocatable :: centered(:), initial(:), hessian(:, :), inverse(:, :), jacobian(:, :)
      real(dp) :: mean_value, gradient_tolerance, step
      integer :: parameter_count, limit, status, i, j

      parameter_count = ar_order + ma_order
      if (ar_order < 0 .or. ma_order < 0 .or. size(series) <= parameter_count + 2) then
         model%info = 1
         return
      end if
      mean_value = sum(series)/real(size(series), dp)
      centered = series - mean_value
      allocate(initial(parameter_count))
      initial = 0.0_dp
      if (ma_order > 0) then
         initial_model = hannan_rissanen_fit(centered, ar_order, ma_order)
         if (initial_model%info == 0) then
            if (ar_order > 0) initial(:ar_order) = ar_to_unconstrained(initial_model%ar)
            initial(ar_order + 1:) = ar_to_unconstrained(-initial_model%ma)
         end if
      else if (ar_order > 0) then
         ar_initial = yule_walker_fit(centered, ar_order)
         if (ar_initial%info == 0) initial = ar_to_unconstrained(ar_initial%coefficients)
      end if
      if (parameter_count == 0) then
         allocate(model%ar(0), model%ma(0), model%ar_standard_error(0), model%ma_standard_error(0))
         call update_innovation_statistics(centered, model)
         model%converged = model%info == 0
         return
      end if
      limit = 200
      if (present(max_iterations)) limit = max_iterations
      gradient_tolerance = 1.0e-6_dp
      if (present(tolerance)) gradient_tolerance = tolerance
      optimization = bfgs_minimize_fd(objective, initial, limit, gradient_tolerance)
      model = model_from_coordinates(optimization%parameters)
      call update_innovation_statistics(centered, model)
      model%iterations = optimization%iterations
      model%converged = optimization%converged
      if (optimization%info /= 0 .and. model%info == 0) model%info = 100 + optimization%info
      hessian = finite_difference_hessian(objective, optimization%parameters)
      call invert_matrix(hessian, inverse, status)
      if (status == 0) then
         allocate(jacobian(parameter_count, parameter_count))
         step = epsilon(1.0_dp)**0.25_dp
         do j = 1, parameter_count
            jacobian(:, j) = transformed_difference(optimization%parameters, j, step)
         end do
         inverse = matmul(matmul(jacobian, inverse), transpose(jacobian))
         do i = 1, ar_order
            model%ar_standard_error(i) = sqrt(max(0.0_dp, inverse(i, i)))
         end do
         do i = 1, ma_order
            model%ma_standard_error(i) = sqrt(max(0.0_dp, inverse(ar_order + i, ar_order + i)))
         end do
      else if (model%info == 0) then
         model%info = 200 + status
      end if

   contains

      pure function objective(coordinates) result(value)
         !! Return the profiled negative innovations log likelihood.
         real(dp), intent(in) :: coordinates(:) !! Coordinates.
         real(dp) :: value
         type(itsmr_arma_model_t) :: candidate

         candidate = model_from_coordinates(coordinates)
         call update_innovation_statistics(centered, candidate)
         if (candidate%info == 0 .and. ieee_is_finite(candidate%aicc)) then
            value = 0.5_dp*(candidate%aicc - &
               2.0_dp*real((parameter_count + 1)*size(centered), dp)/ &
               real(size(centered) - parameter_count - 2, dp))
         else
            value = 1.0e30_dp + dot_product(coordinates, coordinates)
         end if
      end function objective

      pure function model_from_coordinates(coordinates) result(candidate)
         !! Map unconstrained coordinates to causal AR and invertible MA coefficients.
         real(dp), intent(in) :: coordinates(:) !! Coordinates.
         type(itsmr_arma_model_t) :: candidate

         allocate(candidate%ar(ar_order), candidate%ma(ma_order))
         allocate(candidate%ar_standard_error(ar_order), candidate%ma_standard_error(ma_order))
         candidate%ar_standard_error = 0.0_dp
         candidate%ma_standard_error = 0.0_dp
         if (ar_order > 0) candidate%ar = partial_to_ar(tanh(coordinates(:ar_order)))
         if (ma_order > 0) candidate%ma = -partial_to_ar(tanh(coordinates(ar_order + 1:)))
      end function model_from_coordinates

      pure function transformed_difference(coordinates, column, difference_step) result(derivative)
         !! Differentiate reported coefficients with respect to optimizer coordinates.
         real(dp), intent(in) :: coordinates(:) !! Coordinates.
         real(dp), intent(in) :: difference_step !! Difference step.
         integer, intent(in) :: column !! Column.
         real(dp) :: derivative(size(coordinates)), plus(size(coordinates)), minus(size(coordinates)), h
         type(itsmr_arma_model_t) :: plus_model, minus_model

         h = difference_step*max(1.0_dp, abs(coordinates(column)))
         plus = coordinates
         minus = coordinates
         plus(column) = plus(column) + h
         minus(column) = minus(column) - h
         plus_model = model_from_coordinates(plus)
         minus_model = model_from_coordinates(minus)
         derivative = ([plus_model%ar, plus_model%ma] - [minus_model%ar, minus_model%ma])/(2.0_dp*h)
      end function transformed_difference
   end function arma_mle_fit

   pure function arma_autofit(series, ar_orders, ma_orders, max_iterations, tolerance) result(best)
      !! Select the requested ARMA order pair by minimum corrected AIC.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      integer, intent(in) :: ar_orders(:) !! Autoregressive orders.
      integer, intent(in) :: ma_orders(:) !! Moving-average orders.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      type(itsmr_arma_model_t) :: best
      type(itsmr_arma_model_t) :: candidate
      integer :: i, j

      best%aicc = huge(1.0_dp)
      best%info = 1
      do i = 1, size(ar_orders)
         do j = 1, size(ma_orders)
            candidate = arma_mle_fit(series, ar_orders(i), ma_orders(j), max_iterations, tolerance)
            if (candidate%info == 0 .and. candidate%aicc < best%aicc) best = candidate
         end do
      end do
   end function arma_autofit

   pure function partial_to_ar(partial) result(coefficients)
      !! Convert partial autocorrelations to stable AR coefficients.
      real(dp), intent(in) :: partial(:) !! Partial.
      real(dp) :: coefficients(size(partial)), previous(size(partial))
      integer :: order, j

      coefficients = 0.0_dp
      do order = 1, size(partial)
         previous = coefficients
         coefficients(order) = partial(order)
         do j = 1, order - 1
            coefficients(j) = previous(j) - partial(order)*previous(order - j)
         end do
      end do
   end function partial_to_ar

   pure function ar_to_unconstrained(coefficients) result(coordinates)
      !! Convert stable AR coefficients to unconstrained optimizer coordinates.
      real(dp), intent(in) :: coefficients(:) !! Model coefficients.
      real(dp) :: coordinates(size(coefficients)), current(size(coefficients)), previous(size(coefficients)), reflection
      integer :: order, j

      current = coefficients
      do order = size(coefficients), 1, -1
         reflection = max(-0.9999_dp, min(0.9999_dp, current(order)))
         coordinates(order) = atanh(reflection)
         previous = current
         do j = 1, order - 1
            current(j) = (previous(j) + reflection*previous(order - j))/(1.0_dp - reflection**2)
         end do
      end do
   end function ar_to_unconstrained

   pure function arma_infinite_ma(model, max_lag) result(psi)
      !! Return MA-infinity coefficients from lag zero through max_lag.
      type(itsmr_arma_model_t), intent(in) :: model !! Model specification.
      integer, intent(in) :: max_lag !! Maximum lag to consider.
      real(dp), allocatable :: psi(:)
      real(dp) :: ma_coefficient
      integer :: p, q, lag, ar_lag

      if (max_lag < 0 .or. .not. allocated(model%ar) .or. .not. allocated(model%ma)) then
         allocate(psi(0))
         return
      end if
      p = size(model%ar)
      q = size(model%ma)
      allocate(psi(0:max_lag))
      psi = 0.0_dp
      psi(0) = 1.0_dp
      do lag = 1, max_lag
         ma_coefficient = 0.0_dp
         if (lag <= q) ma_coefficient = model%ma(lag)
         psi(lag) = ma_coefficient
         do ar_lag = 1, min(p, lag)
            psi(lag) = psi(lag) + model%ar(ar_lag)*psi(lag - ar_lag)
         end do
      end do
   end function arma_infinite_ma

   pure function arma_infinite_ar(model, max_lag) result(pi_coefficients)
      !! Return AR-infinity coefficients from lag zero through max_lag.
      type(itsmr_arma_model_t), intent(in) :: model !! Model specification.
      integer, intent(in) :: max_lag !! Maximum lag to consider.
      real(dp), allocatable :: pi_coefficients(:)
      real(dp) :: ar_coefficient
      integer :: p, q, lag, ma_lag

      if (max_lag < 0 .or. .not. allocated(model%ar) .or. .not. allocated(model%ma)) then
         allocate(pi_coefficients(0))
         return
      end if
      p = size(model%ar)
      q = size(model%ma)
      allocate(pi_coefficients(0:max_lag))
      pi_coefficients = 0.0_dp
      pi_coefficients(0) = 1.0_dp
      do lag = 1, max_lag
         ar_coefficient = 0.0_dp
         if (lag <= p) ar_coefficient = model%ar(lag)
         pi_coefficients(lag) = -ar_coefficient
         do ma_lag = 1, min(q, lag)
            pi_coefficients(lag) = pi_coefficients(lag) - &
               model%ma(ma_lag)*pi_coefficients(lag - ma_lag)
         end do
      end do
   end function arma_infinite_ar

   pure function arma_residuals(series, model) result(out)
      !! Return exact stationary innovations and their standardized values.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      type(itsmr_arma_model_t), intent(in) :: model !! Model specification.
      type(itsmr_residuals_t) :: out
      type(itsmr_innovations_t) :: recursion
      real(dp), allocatable :: covariance(:), centered(:)
      integer :: n

      n = size(series)
      if (n < 1 .or. model%innovation_variance <= 0.0_dp .or. &
         .not. allocated(model%ar) .or. .not. allocated(model%ma)) then
         out%info = 1
         return
      end if
      out%mean = sum(series)/real(n, dp)
      centered = series - out%mean
      covariance = arma_acvf(model, n - 1)
      if (size(covariance) /= n .or. any(.not. ieee_is_finite(covariance))) then
         out%info = 2
         return
      end if
      recursion = innovations_algorithm(covariance, centered)
      if (recursion%info /= 0) then
         out%info = 10 + recursion%info
         return
      end if
      allocate(out%fitted(n), out%innovations(n), out%variance(n), out%standardized(n))
      out%fitted = recursion%prediction + out%mean
      out%innovations = recursion%innovations
      out%variance = recursion%variance(:n - 1)
      out%standardized = out%innovations/sqrt(out%variance)
   end function arma_residuals

   pure function arma_forecast(series, model, horizon, confidence_level) result(out)
      !! Forecast a fitted ARMA model with normal prediction intervals.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      type(itsmr_arma_model_t), intent(in) :: model !! Model specification.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      real(dp), intent(in), optional :: confidence_level !! Confidence level.
      type(itsmr_forecast_t) :: out
      type(itsmr_innovations_t) :: history
      real(dp), allocatable :: covariance(:), extended(:), extended_innovations(:)
      real(dp) :: mean_value, critical_value, accumulated_variance, ar_part, ma_part
      integer :: n, p, q, step, t, lag

      n = size(series)
      out%confidence_level = 0.95_dp
      if (present(confidence_level)) out%confidence_level = confidence_level
      if (n < 1 .or. horizon < 1 .or. model%innovation_variance <= 0.0_dp .or. &
         out%confidence_level <= 0.0_dp .or. out%confidence_level >= 1.0_dp .or. &
         .not. allocated(model%ar) .or. .not. allocated(model%ma)) then
         out%info = 1
         return
      end if
      p = size(model%ar)
      q = size(model%ma)
      mean_value = sum(series)/real(n, dp)
      covariance = arma_acvf(model, n - 1)
      if (size(covariance) /= n .or. any(.not. ieee_is_finite(covariance))) then
         out%info = 2
         return
      end if
      history = innovations_algorithm(covariance, series - mean_value)
      if (history%info /= 0) then
         out%info = 10 + history%info
         return
      end if
      allocate(extended(n + horizon), extended_innovations(n + horizon))
      extended(:n) = series - mean_value
      extended(n + 1:) = 0.0_dp
      extended_innovations = 0.0_dp
      extended_innovations(:n) = history%innovations
      do step = 1, horizon
         t = n + step
         ar_part = 0.0_dp
         do lag = 1, min(p, t - 1)
            ar_part = ar_part + model%ar(lag)*extended(t - lag)
         end do
         ma_part = 0.0_dp
         do lag = 1, min(q, t - 1)
            ma_part = ma_part + model%ma(lag)*extended_innovations(t - lag)
         end do
         extended(t) = ar_part + ma_part
      end do
      allocate(out%mean(horizon), out%standard_error(horizon), out%lower(horizon), out%upper(horizon))
      allocate(out%innovations(n))
      out%mean = extended(n + 1:) + mean_value
      out%innovations = history%innovations
      allocate(out%psi(0:horizon - 1))
      out%psi = arma_infinite_ma(model, horizon - 1)
      accumulated_variance = 0.0_dp
      do step = 1, horizon
         accumulated_variance = accumulated_variance + out%psi(step - 1)**2
         out%standard_error(step) = sqrt(model%innovation_variance*accumulated_variance)
      end do
      critical_value = inverse_standard_normal(0.5_dp*(1.0_dp + out%confidence_level))
      out%lower = out%mean - critical_value*out%standard_error
      out%upper = out%mean + critical_value*out%standard_error
   end function arma_forecast

   pure function arar_forecast(series, horizon, confidence_level) result(out)
      !! Forecast by ITSMR's automatic memory-shortened sparse AR algorithm.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      real(dp), intent(in), optional :: confidence_level !! Confidence level.
      type(itsmr_arar_t) :: out
      real(dp), allocatable :: shortened(:), next_shortened(:), memory(:), next_memory(:)
      real(dp), allocatable :: centered(:), covariance(:), matrix(:, :), inverse(:, :), right(:)
      real(dp), allocatable :: candidate(:), best_coefficients(:), filter(:), extended(:), impulse(:)
      real(dp) :: phi(15), error_ratio(15), denominator, numerator, best_variance
      real(dp) :: variance, intercept, critical_value
      integer :: n, iteration, lag, selected_lag, i, j, k, row, column, status, filter_lag, step, t
      integer :: candidate_lags(4), best_lags(4)

      out%confidence_level = 0.95_dp
      if (present(confidence_level)) out%confidence_level = confidence_level
      if (size(series) <= 26 .or. horizon < 1 .or. out%confidence_level <= 0.0_dp .or. &
         out%confidence_level >= 1.0_dp) then
         out%info = 1
         return
      end if
      shortened = series
      allocate(memory(0:0))
      memory = 1.0_dp
      do iteration = 1, 3
         n = size(shortened)
         if (n <= 15) exit
         do lag = 1, 15
            denominator = sum(shortened(:n - lag)**2)
            if (denominator <= tiny(1.0_dp)) then
               phi(lag) = 0.0_dp
               error_ratio(lag) = huge(1.0_dp)
            else
               phi(lag) = dot_product(shortened(lag + 1:), shortened(:n - lag))/denominator
               numerator = sum((shortened(lag + 1:) - phi(lag)*shortened(:n - lag))**2)
               denominator = sum(shortened(lag + 1:)**2)
               error_ratio(lag) = numerator/max(denominator, tiny(1.0_dp))
            end if
         end do
         selected_lag = minloc(error_ratio, dim=1)
         if (error_ratio(selected_lag) <= 8.0_dp/real(n, dp) .or. &
            (phi(selected_lag) >= 0.93_dp .and. selected_lag > 2)) then
            if (n - selected_lag <= 26) exit
            next_shortened = shortened(selected_lag + 1:) - &
               phi(selected_lag)*shortened(:n - selected_lag)
            allocate(next_memory(0:ubound(memory, 1) + selected_lag))
            next_memory = 0.0_dp
            next_memory(:ubound(memory, 1)) = memory
            next_memory(selected_lag:) = next_memory(selected_lag:) - phi(selected_lag)*memory
         else if (phi(selected_lag) >= 0.93_dp) then
            if (n - 2 <= 26) exit
            allocate(matrix(2, 2), right(2))
            matrix(1, 1) = sum(shortened(2:n - 1)**2)
            matrix(1, 2) = dot_product(shortened(:n - 2), shortened(2:n - 1))
            matrix(2, 1) = matrix(1, 2)
            matrix(2, 2) = sum(shortened(:n - 2)**2)
            right(1) = dot_product(shortened(3:), shortened(2:n - 1))
            right(2) = dot_product(shortened(3:), shortened(:n - 2))
            call invert_matrix(matrix, inverse, status)
            if (status /= 0) then
               out%info = 10 + status
               return
            end if
            candidate = matmul(inverse, right)
            next_shortened = shortened(3:) - candidate(1)*shortened(2:n - 1) - &
               candidate(2)*shortened(:n - 2)
            allocate(next_memory(0:ubound(memory, 1) + 2))
            next_memory = 0.0_dp
            next_memory(:ubound(memory, 1)) = memory
            next_memory(1:ubound(memory, 1) + 1) = &
               next_memory(1:ubound(memory, 1) + 1) - candidate(1)*memory
            next_memory(2:) = next_memory(2:) - candidate(2)*memory
            deallocate(matrix, right, inverse)
         else
            exit
         end if
         shortened = next_shortened
         memory = next_memory
         deallocate(next_shortened, next_memory)
      end do
      n = size(shortened)
      if (n <= 26) then
         out%info = 2
         return
      end if
      out%shortened_mean = sum(shortened)/real(n, dp)
      centered = shortened - out%shortened_mean
      allocate(covariance(0:26))
      do lag = 0, 26
         covariance(lag) = dot_product(centered(:n - lag), centered(lag + 1:))/real(n, dp)
      end do
      allocate(matrix(4, 4), right(4), best_coefficients(4))
      best_variance = huge(1.0_dp)
      do i = 2, 24
         do j = i + 1, 25
            do k = j + 1, 26
               candidate_lags = [1, i, j, k]
               do row = 1, 4
                  right(row) = covariance(candidate_lags(row))
                  do column = 1, 4
                     matrix(row, column) = covariance(abs(candidate_lags(row) - candidate_lags(column)))
                  end do
               end do
               call invert_matrix(matrix, inverse, status)
               if (status /= 0) cycle
               candidate = matmul(inverse, right)
               variance = covariance(0) - dot_product(candidate, right)
               if (variance > tiny(1.0_dp) .and. variance < best_variance) then
                  best_variance = variance
                  best_coefficients = candidate
                  best_lags = candidate_lags
               end if
            end do
         end do
      end do
      if (.not. ieee_is_finite(best_variance) .or. best_variance == huge(1.0_dp)) then
         out%info = 3
         return
      end if
      filter_lag = ubound(memory, 1) + best_lags(4)
      allocate(filter(0:filter_lag))
      filter = 0.0_dp
      filter(:ubound(memory, 1)) = memory
      do row = 1, 4
         lag = best_lags(row)
         filter(lag:lag + ubound(memory, 1)) = &
            filter(lag:lag + ubound(memory, 1)) - best_coefficients(row)*memory
      end do
      if (size(series) < filter_lag) then
         out%info = 4
         return
      end if
      allocate(extended(size(series) + horizon))
      extended(:size(series)) = series
      extended(size(series) + 1:) = 0.0_dp
      intercept = (1.0_dp - sum(best_coefficients))*out%shortened_mean
      do step = 1, horizon
         t = size(series) + step
         extended(t) = intercept
         do lag = 1, filter_lag
            extended(t) = extended(t) - filter(lag)*extended(t - lag)
         end do
      end do
      allocate(impulse(0:horizon - 1))
      impulse = 0.0_dp
      impulse(0) = 1.0_dp
      do step = 1, horizon - 1
         do lag = 1, min(step, filter_lag)
            impulse(step) = impulse(step) - filter(lag)*impulse(step - lag)
         end do
      end do
      allocate(out%mean(horizon), out%standard_error(horizon), out%lower(horizon), out%upper(horizon))
      allocate(out%lags(4), out%coefficients(4), out%memory_polynomial(0:ubound(memory, 1)))
      allocate(out%filter(0:filter_lag))
      out%mean = extended(size(series) + 1:)
      do step = 1, horizon
         out%standard_error(step) = sqrt(best_variance*sum(impulse(:step - 1)**2))
      end do
      critical_value = inverse_standard_normal(0.5_dp*(1.0_dp + out%confidence_level))
      out%lower = out%mean - critical_value*out%standard_error
      out%upper = out%mean + critical_value*out%standard_error
      out%lags = best_lags
      out%coefficients = best_coefficients
      out%memory_polynomial = memory
      out%filter = filter
      out%innovation_variance = best_variance
   end function arar_forecast

   pure function burg_ar_fit(series, order) result(model)
      !! Estimate an ITSMR AR model using the shared Burg recursion.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      integer, intent(in) :: order !! Model or polynomial order.
      type(itsmr_arma_model_t) :: model
      type(burg_result_t) :: fitted
      real(dp), allocatable :: centered(:)
      integer :: i

      fitted = burg_fit(series, order)
      if (fitted%info /= 0) then
         model%info = fitted%info
         return
      end if
      allocate(model%ar(order), model%ma(0), model%ar_standard_error(order))
      allocate(model%ma_standard_error(0))
      model%ar = fitted%coefficients
      do i = 1, order
         model%ar_standard_error(i) = sqrt(max(0.0_dp, fitted%covariance(i, i)))
      end do
      centered = series - sum(series)/real(size(series), dp)
      call update_innovation_statistics(centered, model)
   end function burg_ar_fit

   pure function transformed_arma_forecast(series, model, transforms, horizon, &
      confidence_level) result(out)
      !! Forecast an ARMA residual model through a typed reversible transform pipeline.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      type(itsmr_arma_model_t), intent(in) :: model !! Model specification.
      type(itsmr_transform_t), intent(in) :: transforms(:) !! Transformation specifications.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      real(dp), intent(in), optional :: confidence_level !! Confidence level.
      type(itsmr_transformed_forecast_t) :: out
      type(itsmr_transform_state_t), allocatable :: states(:)
      type(harmonic_regression_result_t) :: regression
      type(itsmr_forecast_t) :: base
      real(dp), allocatable :: work(:), next_work(:), component(:), polynomial(:), factor(:)
      real(dp), allocatable :: effective_ar(:), psi(:), restored(:)
      real(dp) :: critical_value, accumulated_variance
      integer :: count, index, n, lag, step
      logical :: has_log

      count = size(transforms)
      out%confidence_level = 0.95_dp
      if (present(confidence_level)) out%confidence_level = confidence_level
      if (size(series) < 2 .or. horizon < 1 .or. out%confidence_level <= 0.0_dp .or. &
         out%confidence_level >= 1.0_dp) then
         out%info = 1
         return
      end if
      has_log = .false.
      do index = 1, count
         if (transforms(index)%kind == itsmr_transform_log) then
            if (index /= 1 .or. has_log) then
               out%info = 2
               return
            end if
            has_log = .true.
         end if
      end do
      allocate(states(count))
      work = series
      do index = 1, count
         states(index)%kind = transforms(index)%kind
         select case (transforms(index)%kind)
         case (itsmr_transform_log)
            if (any(work <= 0.0_dp)) then
               out%info = 10 + index
               return
            end if
            work = log(work)
         case (itsmr_transform_difference)
            lag = transforms(index)%lag
            n = size(work)
            if (lag < 1 .or. lag >= n) then
               out%info = 20 + index
               return
            end if
            states(index)%lag = lag
            states(index)%history = work(n - lag + 1:n)
            next_work = work(lag + 1:) - work(:n - lag)
            work = next_work
         case (itsmr_transform_harmonic)
            if (.not. allocated(transforms(index)%periods)) then
               out%info = 30 + index
               return
            end if
            regression = harmonic_regression(work, transforms(index)%periods, horizon)
            if (regression%info /= 0) then
               out%info = 40 + index
               return
            end if
            states(index)%restoration = regression%forecast
            work = regression%residuals
         case (itsmr_transform_trend)
            regression = harmonic_regression(work, [real(dp) ::], horizon, transforms(index)%order)
            if (regression%info /= 0) then
               out%info = 50 + index
               return
            end if
            states(index)%restoration = regression%forecast
            work = regression%residuals
         case (itsmr_transform_seasonal)
            component = seasonal_component(work, transforms(index)%lag, horizon)
            if (size(component) /= size(work) + horizon) then
               out%info = 60 + index
               return
            end if
            states(index)%restoration = component(size(work) + 1:)
            work = work - component(:size(work))
         case default
            out%info = 70 + index
            return
         end select
      end do
      out%transformed_series = work
      base = arma_forecast(work, model, horizon, out%confidence_level)
      if (base%info /= 0) then
         out%info = 100 + base%info
         return
      end if
      restored = base%mean
      if (allocated(next_work)) deallocate(next_work)
      do index = count, 1, -1
         select case (states(index)%kind)
         case (itsmr_transform_difference)
            lag = states(index)%lag
            allocate(next_work(horizon))
            do step = 1, horizon
               if (step <= lag) then
                  next_work(step) = restored(step) + states(index)%history(step)
               else
                  next_work(step) = restored(step) + next_work(step - lag)
               end if
            end do
            restored = next_work
            deallocate(next_work)
         case (itsmr_transform_harmonic, itsmr_transform_trend, itsmr_transform_seasonal)
            restored = restored + states(index)%restoration
         end select
      end do
      polynomial = [1.0_dp, -model%ar]
      do index = 1, count
         if (transforms(index)%kind /= itsmr_transform_difference) cycle
         lag = transforms(index)%lag
         allocate(factor(0:lag))
         factor = 0.0_dp
         factor(0) = 1.0_dp
         factor(lag) = -1.0_dp
         polynomial = polynomial_product(polynomial, factor)
         deallocate(factor)
      end do
      effective_ar = -polynomial(2:)
      allocate(psi(0:horizon - 1))
      psi = arma_weights(effective_ar, model%ma, horizon - 1)
      allocate(out%mean(horizon), out%standard_error(horizon), out%lower(horizon), out%upper(horizon))
      out%mean = restored
      accumulated_variance = 0.0_dp
      do step = 1, horizon
         accumulated_variance = accumulated_variance + psi(step - 1)**2
         out%standard_error(step) = sqrt(model%innovation_variance*accumulated_variance)
      end do
      critical_value = inverse_standard_normal(0.5_dp*(1.0_dp + out%confidence_level))
      out%lower = out%mean - critical_value*out%standard_error
      out%upper = out%mean + critical_value*out%standard_error
      if (has_log) then
         out%mean = exp(out%mean)
         out%lower = exp(out%lower)
         out%upper = exp(out%upper)
         out%standard_error = out%mean*out%standard_error
      end if
   end function transformed_arma_forecast

   pure function arma_weights(ar, ma, max_lag) result(weights)
      !! Return zero-based MA-infinity weights for supplied ARMA coefficients.
      real(dp), intent(in) :: ar(:) !! Autoregressive coefficients.
      real(dp), intent(in) :: ma(:) !! Moving-average coefficients.
      integer, intent(in) :: max_lag !! Maximum lag to consider.
      real(dp) :: weights(0:max_lag), ma_coefficient
      integer :: lag, ar_lag

      weights = 0.0_dp
      weights(0) = 1.0_dp
      do lag = 1, max_lag
         ma_coefficient = 0.0_dp
         if (lag <= size(ma)) ma_coefficient = ma(lag)
         weights(lag) = ma_coefficient
         do ar_lag = 1, min(size(ar), lag)
            weights(lag) = weights(lag) + ar(ar_lag)*weights(lag - ar_lag)
         end do
      end do
   end function arma_weights

   pure function seasonal_component(series, period, horizon) result(component)
      !! Estimate and extend ITSMR's classical seasonal component.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      integer, intent(in) :: period !! Seasonal period.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      real(dp), allocatable :: component(:)
      real(dp), allocatable :: moving_average(:), deviations(:), seasonal_means(:)
      real(dp) :: total
      integer :: n, q, t, phase, first, last, count

      n = size(series)
      q = period/2
      if (period < 2 .or. horizon < 0 .or. n <= 2*q) then
         allocate(component(0))
         return
      end if
      allocate(moving_average(n), deviations(n), seasonal_means(period))
      if (mod(period, 2) == 0) then
         moving_average = 0.0_dp
         do t = q + 1, n - q
            moving_average(t) = 0.5_dp*series(t - q) + 0.5_dp*series(t + q)
            if (q > 1) moving_average(t) = moving_average(t) + sum(series(t - q + 1:t + q - 1))
            moving_average(t) = moving_average(t)/real(period, dp)
         end do
      else
         do t = 1, n
            total = 0.0_dp
            do phase = -q, q
               total = total + series(max(1, min(n, t + phase)))
            end do
            moving_average(t) = total/real(period, dp)
         end do
      end if
      deviations = series - moving_average
      do phase = 1, period
         first = phase + q
         last = n - q
         if (first > last) then
            allocate(component(0))
            return
         end if
         total = 0.0_dp
         count = 0
         do t = first, last, period
            total = total + deviations(t)
            count = count + 1
         end do
         seasonal_means(phase) = total/real(count, dp)
      end do
      seasonal_means = seasonal_means - sum(seasonal_means)/real(period, dp)
      allocate(component(n + horizon))
      do t = 1, n + horizon
         phase = modulo(q + t - 1, period) + 1
         component(t) = seasonal_means(phase)
      end do
   end function seasonal_component

   pure function spectral_rank_filter(series, retained_count) result(out)
      !! Reconstruct a series from its strongest positive-frequency DFT bins.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      integer, intent(in) :: retained_count !! Number of retained.
      type(itsmr_rank_filter_t) :: out
      complex(dp), allocatable :: coefficients(:)
      real(dp), allocatable :: magnitude(:)
      integer, allocatable :: order(:), retained(:)
      real(dp) :: value
      integer :: n, positive_count, i, j, selected, temporary

      n = size(series)
      positive_count = n/2
      if (n < 2 .or. retained_count < 0 .or. retained_count > positive_count) then
         out%info = 1
         return
      end if
      coefficients = real_dft(series)
      allocate(magnitude(positive_count), order(positive_count))
      do i = 1, positive_count
         magnitude(i) = abs(coefficients(i + 1))
         order(i) = i
      end do
      do i = 1, positive_count - 1
         selected = i
         do j = i + 1, positive_count
            if (magnitude(order(j)) < magnitude(order(selected))) selected = j
         end do
         if (selected /= i) then
            temporary = order(i)
            order(i) = order(selected)
            order(selected) = temporary
         end if
      end do
      do i = 1, positive_count - retained_count
         coefficients(order(i) + 1) = cmplx(0.0_dp, 0.0_dp, dp)
         coefficients(n - order(i) + 1) = cmplx(0.0_dp, 0.0_dp, dp)
      end do
      allocate(retained(retained_count))
      if (retained_count > 0) retained = order(positive_count - retained_count + 1:)
      do i = 1, retained_count - 1
         selected = i
         do j = i + 1, retained_count
            if (retained(j) < retained(selected)) selected = j
         end do
         if (selected /= i) then
            temporary = retained(i)
            retained(i) = retained(selected)
            retained(selected) = temporary
         end if
      end do
      allocate(out%retained_bins(retained_count), out%frequency(retained_count))
      allocate(out%amplitude(retained_count))
      out%retained_bins = retained
      do i = 1, retained_count
         out%frequency(i) = real(retained(i), dp)/real(n, dp)
         value = 2.0_dp*abs(coefficients(retained(i) + 1))
         if (mod(n, 2) == 0 .and. retained(i) == n/2) value = 0.5_dp*value
         out%amplitude(i) = value
      end do
      out%coefficients = coefficients
      out%filtered = inverse_real_dft(coefficients)
   end function spectral_rank_filter

   pure function itsmr_moving_average(series, half_window) result(smoothed)
      !! Apply ITSMR's endpoint-replicated centered moving average.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      integer, intent(in) :: half_window !! Half window.
      real(dp), allocatable :: smoothed(:)
      real(dp) :: total
      integer :: n, t, offset

      n = size(series)
      if (n < 1 .or. half_window < 0) then
         allocate(smoothed(0))
         return
      end if
      allocate(smoothed(n))
      do t = 1, n
         total = 0.0_dp
         do offset = -half_window, half_window
            total = total + series(max(1, min(n, t + offset)))
         end do
         smoothed(t) = total/real(2*half_window + 1, dp)
      end do
   end function itsmr_moving_average

   pure function itsmr_exponential_smooth(series, alpha) result(smoothed)
      !! Apply ITSMR's recursively initialized exponential smoother.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      real(dp), intent(in) :: alpha !! Significance, smoothing, or model coefficient.
      real(dp), allocatable :: smoothed(:)
      integer :: n, t

      n = size(series)
      if (n < 1 .or. alpha < 0.0_dp .or. alpha > 1.0_dp) then
         allocate(smoothed(0))
         return
      end if
      allocate(smoothed(n))
      smoothed(1) = series(1)
      do t = 2, n
         smoothed(t) = alpha*series(t) + (1.0_dp - alpha)*smoothed(t - 1)
      end do
   end function itsmr_exponential_smooth

   pure function itsmr_fft_smooth(series, pass_fraction) result(smoothed)
      !! Apply ITSMR's symmetric-index low-pass Fourier filter.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      real(dp), intent(in) :: pass_fraction !! Pass fraction.
      real(dp), allocatable :: smoothed(:)
      complex(dp), allocatable :: coefficients(:)
      complex(dp) :: value
      real(dp) :: angle, root_n
      integer :: n, width, retained, frequency, t, index

      n = size(series)
      if (n < 1 .or. pass_fraction < 0.0_dp .or. pass_fraction > 1.0_dp) then
         allocate(smoothed(0))
         return
      end if
      width = n/2
      retained = floor(pass_fraction*real(width, dp))
      root_n = sqrt(real(n, dp))
      allocate(coefficients(-width:width), smoothed(n))
      coefficients = cmplx(0.0_dp, 0.0_dp, dp)
      do frequency = -retained, retained
         do t = 1, n
            angle = -2.0_dp*acos(-1.0_dp)*real(frequency*t, dp)/real(n, dp)
            coefficients(frequency) = coefficients(frequency) + &
               series(t)*cmplx(cos(angle), sin(angle), dp)
         end do
         coefficients(frequency) = coefficients(frequency)/root_n
      end do
      do t = 1, n
         value = cmplx(0.0_dp, 0.0_dp, dp)
         do index = -width, width
            angle = 2.0_dp*acos(-1.0_dp)*real(index*t, dp)/real(n, dp)
            value = value + coefficients(index)*cmplx(cos(angle), sin(angle), dp)
         end do
         smoothed(t) = real(value, dp)/root_n
      end do
   end function itsmr_fft_smooth

   pure function itsmr_seasonal_component(series, period) result(component)
      !! Return ITSMR's classical seasonal component over the observed sample.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      integer, intent(in) :: period !! Seasonal period.
      real(dp), allocatable :: component(:)
      component = seasonal_component(series, period, 0)
   end function itsmr_seasonal_component

   pure function residual_randomness_tests(residuals, max_lag) result(out)
      !! Compute ITSMR's five numerical residual randomness tests.
      real(dp), intent(in) :: residuals(:) !! Model residuals.
      integer, intent(in), optional :: max_lag !! Maximum lag to consider.
      type(itsmr_randomness_tests_t) :: out
      real(dp), allocatable :: correlation(:), squared(:)
      real(dp) :: mean_value, expected, standard_deviation
      integer :: n, lag, i, j, count

      n = size(residuals)
      lag = min(20, n - 1)
      if (present(max_lag)) lag = max_lag
      if (n < 3 .or. lag < 1 .or. lag >= n) then
         out%info = 1
         return
      end if
      out%lag = lag
      allocate(correlation(0:lag))
      correlation = biased_correlations(residuals, lag)
      do i = 1, lag
         out%ljung_box = out%ljung_box + &
            real(n*(n + 2), dp)*correlation(i)**2/real(n - i, dp)
      end do
      out%ljung_box_p_value = regularized_gamma_q(0.5_dp*real(lag, dp), 0.5_dp*out%ljung_box)
      squared = residuals**2
      correlation = biased_correlations(squared, lag)
      do i = 1, lag
         out%mcleod_li = out%mcleod_li + &
            real(n*(n + 2), dp)*correlation(i)**2/real(n - i, dp)
      end do
      out%mcleod_li_p_value = regularized_gamma_q(0.5_dp*real(lag, dp), 0.5_dp*out%mcleod_li)
      count = 0
      do i = 2, n - 1
         if ((residuals(i) > residuals(i - 1) .and. residuals(i) > residuals(i + 1)) .or. &
            (residuals(i) < residuals(i - 1) .and. residuals(i) < residuals(i + 1))) count = count + 1
      end do
      out%turning_points = real(count, dp)
      expected = 2.0_dp*real(n - 2, dp)/3.0_dp
      standard_deviation = sqrt(real(16*n - 29, dp)/90.0_dp)
      out%turning_point_z = (out%turning_points - expected)/standard_deviation
      out%turning_point_p_value = two_sided_normal_p(out%turning_point_z)
      count = 0
      do i = 2, n
         if (residuals(i) > residuals(i - 1)) count = count + 1
      end do
      out%difference_signs = real(count, dp)
      expected = 0.5_dp*real(n - 1, dp)
      standard_deviation = sqrt(real(n + 1, dp)/12.0_dp)
      out%difference_sign_z = (out%difference_signs - expected)/standard_deviation
      out%difference_sign_p_value = two_sided_normal_p(out%difference_sign_z)
      count = 0
      do i = 1, n - 1
         do j = i + 1, n
            if (residuals(j) > residuals(i)) count = count + 1
         end do
      end do
      out%rank_statistic = real(count, dp)
      expected = real(n*(n - 1), dp)/4.0_dp
      standard_deviation = sqrt(real(n*(n - 1)*(2*n + 5), dp)/72.0_dp)
      out%rank_z = (out%rank_statistic - expected)/standard_deviation
      out%rank_p_value = two_sided_normal_p(out%rank_z)
   end function residual_randomness_tests

   pure function biased_correlations(series, max_lag) result(correlation)
      !! Return biased sample correlations from lag zero through max_lag.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      integer, intent(in) :: max_lag !! Maximum lag to consider.
      real(dp) :: correlation(0:max_lag), covariance, mean_value
      integer :: n, lag

      n = size(series)
      mean_value = sum(series)/real(n, dp)
      covariance = sum((series - mean_value)**2)/real(n, dp)
      if (covariance <= tiny(1.0_dp)) then
         correlation = 0.0_dp
         correlation(0) = 1.0_dp
         return
      end if
      correlation(0) = 1.0_dp
      do lag = 1, max_lag
         correlation(lag) = dot_product(series(:n - lag) - mean_value, &
            series(lag + 1:) - mean_value)/(real(n, dp)*covariance)
      end do
   end function biased_correlations

   pure elemental real(dp) function two_sided_normal_p(value) result(probability)
      !! Return a two-sided standard-normal p-value.
      real(dp), intent(in) :: value !! Input value.
      probability = erfc(abs(value)/sqrt(2.0_dp))
   end function two_sided_normal_p


   pure subroutine update_innovation_statistics(centered, model)
      !! Update an ARMA model's innovation variance and corrected AIC.
      real(dp), intent(in) :: centered(:) !! Centered.
      type(itsmr_arma_model_t), intent(inout) :: model !! Model specification, updated in place.
      type(itsmr_innovations_t) :: fitted
      real(dp), allocatable :: covariance(:)
      real(dp) :: log_likelihood
      integer :: n, parameter_count

      n = size(centered)
      parameter_count = size(model%ar) + size(model%ma) + 1
      if (n <= parameter_count + 1) then
         model%info = 30
         return
      end if
      model%innovation_variance = 1.0_dp
      covariance = arma_acvf(model, n - 1)
      if (any(covariance >= huge(1.0_dp))) then
         model%info = 31
         return
      end if
      fitted = innovations_algorithm(covariance, centered)
      if (fitted%info /= 0) then
         model%info = 40 + fitted%info
         return
      end if
      model%innovation_variance = sum(fitted%innovations**2/fitted%variance(:n - 1))/real(n, dp)
      if (model%innovation_variance <= tiny(1.0_dp)) then
         model%info = 32
         return
      end if
      log_likelihood = -0.5_dp*real(n, dp)*log(2.0_dp*acos(-1.0_dp)*model%innovation_variance) - &
         0.5_dp*sum(log(fitted%variance(:n - 1))) - 0.5_dp*real(n, dp)
      model%aicc = -2.0_dp*log_likelihood + &
         2.0_dp*real(parameter_count*n, dp)/real(n - parameter_count - 1, dp)
   end subroutine update_innovation_statistics
end module itsmr_mod
