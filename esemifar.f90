! SPDX-License-Identifier: GPL-3.0-only
! SPDX-FileComment: Algorithms translated from the R esemifar package.
! Extended SEMIFAR smoothing and forecasting algorithms translated from esemifar.
module esemifar_mod
   use kind_mod, only: dp
   use arfima_mod, only: arfima_fractional_weights
   use fracdiff_mod, only: fracdiff_fit_t, fracdiff_simulation_t, fracdiff_fit, &
      fracdiff_simulate_from_innovations
   use linalg_mod, only: invert_matrix
   use random_mod, only: set_random_seed, random_uniform
   use stats_mod, only: normal_quantile, sort, quantile
   use polynomial_mod, only: polynomial_product => polynomial_product_truncated
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   private

   integer, parameter, public :: esemifar_inflation_optimal = 1
   integer, parameter, public :: esemifar_inflation_naive = 2
   integer, parameter, public :: esemifar_inflation_stable = 3

   type, public :: esemifar_smooth_t
      ! Local-polynomial estimates, residuals, weights, and smoothing settings.
      real(dp), allocatable :: estimate(:), residuals(:), weights(:, :)
      real(dp) :: bandwidth = 0.0_dp
      integer :: derivative_order = 0
      integer :: polynomial_order = 1
      integer :: kernel_smoothness = 1
      integer :: info = 0
   end type esemifar_smooth_t

   type, public :: esemifar_order_selection_t
      ! FARIMA information-criterion grid and selected model.
      real(dp), allocatable :: criterion(:, :)
      type(fracdiff_fit_t) :: fit
      integer :: ar_order = 0
      integer :: ma_order = 0
      integer :: info = 0
   end type esemifar_order_selection_t

   type, public :: esemifar_model_t
      ! IPI trend estimate, selected FARIMA errors, and convergence history.
      type(esemifar_smooth_t) :: smoother
      type(fracdiff_fit_t) :: farima
      real(dp), allocatable :: bandwidth_history(:)
      real(dp) :: variance_factor = 0.0_dp
      real(dp) :: curvature_integral = 0.0_dp
      integer :: iterations = 0
      integer :: info = 0
      logical :: converged = .false.
   end type esemifar_model_t

   type, public :: esemifar_forecast_t
      ! Trend-plus-FARIMA forecasts and interval bounds.
      real(dp), allocatable :: mean(:), lower(:, :), upper(:, :)
      integer :: info = 0
   end type esemifar_forecast_t

   public :: esemifar_arma_to_ar, esemifar_arma_to_ma, esemifar_d_coefficients
   public :: esemifar_farima_to_ar, esemifar_farima_to_ma
   public :: esemifar_kdf, esemifar_smooth, esemifar_order_selection
   public :: esemifar_trend_fit, esemifar_derivative_fit
   public :: esemifar_derivative_ipi
   public :: esemifar_forecast_normal, esemifar_forecast_bootstrap
   public :: esemifar_forecast_bootstrap_advanced

contains

   pure function esemifar_arma_to_ma(ar, ma, max_lag) result(coefficients)
      !! Expand the positive-sign ARMA model into its infinite MA representation.
      real(dp), intent(in) :: ar(:) !! Autoregressive coefficients.
      real(dp), intent(in) :: ma(:) !! Moving-average coefficients.
      integer, intent(in) :: max_lag !! Maximum lag to consider.
      real(dp), allocatable :: coefficients(:)
      integer :: lag, j

      if (max_lag < 0) then
         allocate(coefficients(0))
         return
      end if
      allocate(coefficients(0:max_lag))
      coefficients = 0.0_dp
      coefficients(0) = 1.0_dp
      do lag = 1, max_lag
         if (lag <= size(ma)) coefficients(lag) = ma(lag)
         do j = 1, min(lag, size(ar))
            coefficients(lag) = coefficients(lag) + ar(j)*coefficients(lag - j)
         end do
      end do
   end function esemifar_arma_to_ma

   pure function esemifar_arma_to_ar(ar, ma, max_lag) result(coefficients)
      !! Expand minus the ARMA innovation filter into its infinite AR representation.
      real(dp), intent(in) :: ar(:) !! Autoregressive coefficients.
      real(dp), intent(in) :: ma(:) !! Moving-average coefficients.
      integer, intent(in) :: max_lag !! Maximum lag to consider.
      real(dp), allocatable :: coefficients(:), inverse(:)
      integer :: lag, j

      if (max_lag < 0) then
         allocate(coefficients(0))
         return
      end if
      allocate(inverse(0:max_lag), coefficients(0:max_lag))
      inverse = 0.0_dp
      inverse(0) = 1.0_dp
      do lag = 1, max_lag
         if (lag <= size(ar)) inverse(lag) = -ar(lag)
         do j = 1, min(lag, size(ma))
            inverse(lag) = inverse(lag) - ma(j)*inverse(lag - j)
         end do
      end do
      coefficients = -inverse
   end function esemifar_arma_to_ar

   pure function esemifar_d_coefficients(d, max_lag) result(coefficients)
      !! Expand the fractional differencing operator.
      real(dp), intent(in) :: d !! Fractional-differencing parameter or differencing order.
      integer, intent(in) :: max_lag !! Maximum lag to consider.
      real(dp), allocatable :: coefficients(:)

      coefficients = arfima_fractional_weights(d, max_lag)
   end function esemifar_d_coefficients

   pure function esemifar_farima_to_ar(ar, ma, d, max_lag) result(coefficients)
      !! Expand minus the complete FARIMA innovation filter.
      real(dp), intent(in) :: ar(:) !! Autoregressive coefficients.
      real(dp), intent(in) :: ma(:) !! Moving-average coefficients.
      real(dp), intent(in) :: d !! Fractional-differencing parameter or differencing order.
      integer, intent(in) :: max_lag !! Maximum lag to consider.
      real(dp), allocatable :: coefficients(:), short_filter(:), fractional(:)

      short_filter = -esemifar_arma_to_ar(ar, ma, max_lag)
      fractional = arfima_fractional_weights(d, max_lag)
      coefficients = -polynomial_product(short_filter, fractional, max_lag)
   end function esemifar_farima_to_ar

   pure function esemifar_farima_to_ma(ar, ma, d, max_lag) result(coefficients)
      !! Expand the complete FARIMA impulse-response filter.
      real(dp), intent(in) :: ar(:) !! Autoregressive coefficients.
      real(dp), intent(in) :: ma(:) !! Moving-average coefficients.
      real(dp), intent(in) :: d !! Fractional-differencing parameter or differencing order.
      integer, intent(in) :: max_lag !! Maximum lag to consider.
      real(dp), allocatable :: coefficients(:), short_memory(:), fractional(:)

      short_memory = esemifar_arma_to_ma(ar, ma, max_lag)
      fractional = arfima_fractional_weights(-d, max_lag)
      coefficients = polynomial_product(short_memory, fractional, max_lag)
   end function esemifar_farima_to_ma

   pure real(dp) function esemifar_kdf(l, m, d) result(value)
      !! Evaluate esemifar's long-memory kernel double-integral constant.
      integer, intent(in) :: l !! L.
      integer, intent(in) :: m !! M.
      real(dp), intent(in) :: d !! Fractional-differencing parameter or differencing order.
      real(dp) :: inner
      integer :: i, j

      value = 0.0_dp
      do i = 0, m
         inner = 0.0_dp
         do j = i, l + m
            inner = inner + (-1.0_dp)**(j - i)*binomial(l + m - i, j - i)* &
               2.0_dp**(2.0_dp*d + real(j + 1, dp))/ &
               (2.0_dp*d + real(j + 1, dp))
         end do
         value = value + 2.0_dp*binomial(m, i)*inner/(2.0_dp*d + real(i, dp))
      end do
   end function esemifar_kdf

   pure function esemifar_smooth(series, derivative_order, polynomial_order, &
      kernel_smoothness, bandwidth, boundary_nearest) result(out)
      !! Estimate a trend derivative by boundary-aware local polynomial regression.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      real(dp), intent(in) :: bandwidth !! Smoothing or spectral bandwidth.
      integer, intent(in), optional :: derivative_order !! Derivative order.
      integer, intent(in), optional :: polynomial_order !! Polynomial order.
      integer, intent(in), optional :: kernel_smoothness !! Kernel smoothness.
      logical, intent(in), optional :: boundary_nearest !! Flag controlling boundary nearest.
      type(esemifar_smooth_t) :: out
      real(dp), allocatable :: design(:, :), weighted_cross(:, :), inverse(:, :)
      real(dp), allocatable :: response(:), weight(:), offset(:)
      integer :: n, derivative, polynomial, mu, half_width, target, left, right
      integer :: count, i, j, status
      logical :: nearest

      n = size(series)
      derivative = 0
      if (present(derivative_order)) derivative = derivative_order
      polynomial = derivative + 1
      if (present(polynomial_order)) polynomial = polynomial_order
      mu = 1
      if (present(kernel_smoothness)) mu = kernel_smoothness
      nearest = .true.
      if (present(boundary_nearest)) nearest = boundary_nearest
      out%bandwidth = bandwidth
      out%derivative_order = derivative
      out%polynomial_order = polynomial
      out%kernel_smoothness = mu
      if (n < 3 .or. derivative < 0 .or. polynomial < derivative .or. &
         mod(polynomial - derivative, 2) /= 1 .or. mu < 0 .or. &
         bandwidth <= 0.0_dp .or. bandwidth >= 0.5_dp .or. &
         .not. all(ieee_is_finite(series))) then
         out%info = 1
         return
      end if
      half_width = max(polynomial + 1, nint(real(n, dp)*bandwidth))
      half_width = min(half_width, (n - 1)/2)
      allocate(out%estimate(n), out%weights(n, n))
      out%weights = 0.0_dp
      do target = 1, n
         if (nearest) then
            left = max(1, min(target - half_width, n - 2*half_width))
            right = min(n, left + 2*half_width)
         else
            left = max(1, target - half_width)
            right = min(n, target + half_width)
         end if
         count = right - left + 1
         allocate(design(count, polynomial + 1), response(count), weight(count), &
            offset(count))
         do i = 1, count
            offset(i) = real(left + i - 1 - target, dp)/real(n, dp)
            weight(i) = max(0.0_dp, 1.0_dp - &
               (offset(i)/bandwidth)**2)**mu
            design(i, 1) = 1.0_dp
            do j = 1, polynomial
               design(i, j + 1) = offset(i)**j
            end do
         end do
         response = series(left:right)
         weighted_cross = matmul(transpose(design), &
            design*spread(weight, 2, polynomial + 1))
         call invert_matrix(weighted_cross, inverse, status)
         if (status /= 0) then
            out%info = 2
            return
         end if
         out%weights(target, left:right) = factorial_real(derivative)* &
            matmul(inverse(derivative + 1, :), transpose(design))*weight
         out%estimate(target) = dot_product(out%weights(target, left:right), response)
         deallocate(design, response, weight, offset, weighted_cross, inverse)
      end do
      if (derivative == 0) then
         out%residuals = series - out%estimate
      else
         allocate(out%residuals(0))
      end if
   end function esemifar_smooth

   pure function esemifar_order_selection(series, p_max, q_max, use_bic, &
      initial_d, truncation, max_iterations) result(out)
      !! Select FARIMA orders by a grid of Haslett-Raftery fits.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      integer, intent(in) :: p_max !! P max.
      integer, intent(in) :: q_max !! Q max.
      logical, intent(in), optional :: use_bic !! Whether to use the bic.
      real(dp), intent(in), optional :: initial_d !! Initial d.
      integer, intent(in), optional :: truncation !! Truncation.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      type(esemifar_order_selection_t) :: out
      type(fracdiff_fit_t) :: candidate
      real(dp), allocatable :: ar(:), ma(:)
      real(dp) :: d_start, value, best
      integer :: p, q
      logical :: bic

      bic = .true.
      if (present(use_bic)) bic = use_bic
      d_start = 0.1_dp
      if (present(initial_d)) d_start = initial_d
      if (size(series) < 4 .or. p_max < 0 .or. q_max < 0) then
         out%info = 1
         return
      end if
      allocate(out%criterion(0:p_max, 0:q_max))
      best = huge(1.0_dp)
      do p = 0, p_max
         allocate(ar(p))
         ar = 0.05_dp
         do q = 0, q_max
            allocate(ma(q))
            ma = 0.05_dp
            candidate = fracdiff_fit(series, ar, ma, d_start, truncation=truncation, &
               max_iterations=max_iterations)
            value = candidate%aic
            if (bic) value = candidate%bic
            out%criterion(p, q) = value
            if (candidate%info == 0 .and. value < best) then
               best = value
               out%ar_order = p
               out%ma_order = q
               out%fit = candidate
            end if
            deallocate(ma)
         end do
         deallocate(ar)
      end do
      if (best == huge(1.0_dp)) out%info = 2
   end function esemifar_order_selection

   pure function esemifar_trend_fit(series, polynomial_order, kernel_smoothness, &
      initial_bandwidth, inflation, p_max, q_max, boundary_nearest, &
      max_iterations) result(out)
      !! Select a trend bandwidth by the ESEMIFAR iterative plug-in algorithm.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      integer, intent(in), optional :: polynomial_order !! Polynomial order.
      integer, intent(in), optional :: kernel_smoothness !! Kernel smoothness.
      integer, intent(in), optional :: inflation !! Inflation.
      integer, intent(in), optional :: p_max !! P max.
      integer, intent(in), optional :: q_max !! Q max.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      real(dp), intent(in), optional :: initial_bandwidth !! Initial bandwidth.
      logical, intent(in), optional :: boundary_nearest !! Flag controlling boundary nearest.
      type(esemifar_model_t) :: out
      type(esemifar_smooth_t) :: smooth, derivative_smooth
      type(esemifar_order_selection_t) :: selection
      real(dp), allocatable :: history(:)
      real(dp) :: bandwidth, old_bandwidth, older_bandwidth, inflated, d, kernel_d
      real(dp) :: c1, c2, c3, rp, moment, i2, variance_factor, candidate
      integer :: polynomial, mu, inflation_method, max_p, max_q, limit, iteration
      integer :: trim, k
      logical :: nearest

      polynomial = 1
      if (present(polynomial_order)) polynomial = polynomial_order
      mu = 1
      if (present(kernel_smoothness)) mu = kernel_smoothness
      inflation_method = esemifar_inflation_optimal
      if (present(inflation)) inflation_method = inflation
      max_p = 0
      if (present(p_max)) max_p = p_max
      max_q = 0
      if (present(q_max)) max_q = q_max
      bandwidth = 0.15_dp
      if (present(initial_bandwidth)) bandwidth = initial_bandwidth
      nearest = .true.
      if (present(boundary_nearest)) nearest = boundary_nearest
      limit = 40
      if (present(max_iterations)) limit = max_iterations
      if ((polynomial /= 1 .and. polynomial /= 3) .or. mu < 0 .or. mu > 3 .or. &
         bandwidth <= 0.0_dp .or. bandwidth >= 0.5_dp .or. limit < 1) then
         out%info = 1
         return
      end if
      allocate(history(limit))
      history = 0.0_dp
      older_bandwidth = bandwidth
      do iteration = 1, limit
         old_bandwidth = bandwidth
         smooth = esemifar_smooth(series, 0, polynomial, mu, bandwidth, nearest)
         selection = esemifar_order_selection(smooth%residuals, max_p, max_q, &
            truncation=min(100, size(series)), max_iterations=120)
         if (selection%info /= 0) then
            out%info = 2
            return
         end if
         d = selection%fit%model%long_memory_parameter
         kernel_d = max(d, sqrt(epsilon(1.0_dp)))
         variance_factor = farima_variance_factor(selection%fit)
         inflated = inflated_bandwidth(bandwidth, d, polynomial, inflation_method)
         inflated = min(0.49_dp, inflated)
         k = polynomial + 1
         derivative_smooth = esemifar_smooth(series, k, polynomial + 2, mu, &
            inflated, nearest)
         trim = int(0.05_dp*real(size(series), dp))
         i2 = sum(derivative_smooth%estimate(max(1, trim): &
            max(1, size(series) - trim))**2)/real(max(1, size(series) - 2*trim), dp)
         rp = 2.0_dp*gamma(1.0_dp - 2.0_dp*kernel_d)* &
            sin(acos(-1.0_dp)*kernel_d)*kernel_rp(polynomial, mu, kernel_d)
         moment = kernel_moment(polynomial, mu, k)
         c1 = factorial_real(k)**2/(2.0_dp*real(k, dp))
         c2 = 0.9_dp*(1.0_dp - 2.0_dp*d)*rp/(moment**2)
         c3 = variance_factor/max(i2, tiny(1.0_dp))
         candidate = (max(c1*c2*c3, tiny(1.0_dp)))** &
            (1.0_dp/(real(2*polynomial + 3, dp) - 2.0_dp*d))* &
            real(size(series), dp)**((2.0_dp*d - 1.0_dp)/ &
            (real(2*polynomial + 3, dp) - 2.0_dp*d))
         candidate = max(candidate, real(size(series), dp)** &
            (-real(2*polynomial + 3, dp)/real(2*polynomial + 5, dp)))
         bandwidth = min(0.49_dp, candidate)
         history(iteration) = bandwidth
         if (iteration > 2) then
            if (abs(old_bandwidth - bandwidth)/bandwidth < &
               1.0_dp/real(size(series), dp)) then
               out%converged = .true.
               exit
            end if
         end if
         if (iteration > 3) then
            if (abs(older_bandwidth - bandwidth)/bandwidth < &
               1.0_dp/real(size(series), dp)) then
               bandwidth = 0.5_dp*(old_bandwidth + bandwidth)
               out%converged = .true.
               exit
            end if
         end if
         older_bandwidth = old_bandwidth
      end do
      out%iterations = min(iteration, limit)
      out%bandwidth_history = history(:out%iterations)
      out%smoother = esemifar_smooth(series, 0, polynomial, mu, bandwidth, nearest)
      selection = esemifar_order_selection(out%smoother%residuals, max_p, max_q, &
         truncation=min(100, size(series)), max_iterations=150)
      out%farima = selection%fit
      out%variance_factor = farima_variance_factor(out%farima)
      out%curvature_integral = i2
      if (selection%info /= 0 .or. out%smoother%info /= 0) out%info = 3
   end function esemifar_trend_fit

   pure function esemifar_derivative_fit(series, derivative_order, bandwidth, &
      polynomial_order, kernel_smoothness, boundary_nearest) result(out)
      !! Estimate a requested trend derivative at a fixed bandwidth.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      real(dp), intent(in) :: bandwidth !! Smoothing or spectral bandwidth.
      integer, intent(in) :: derivative_order !! Derivative order.
      integer, intent(in), optional :: polynomial_order !! Polynomial order.
      integer, intent(in), optional :: kernel_smoothness !! Kernel smoothness.
      logical, intent(in), optional :: boundary_nearest !! Flag controlling boundary nearest.
      type(esemifar_smooth_t) :: out
      integer :: polynomial

      polynomial = derivative_order + 1
      if (present(polynomial_order)) polynomial = polynomial_order
      out = esemifar_smooth(series, derivative_order, polynomial, kernel_smoothness, &
         bandwidth, boundary_nearest)
   end function esemifar_derivative_fit

   pure function esemifar_derivative_ipi(series, derivative_order, pilot_order, &
      kernel_smoothness, pilot_kernel_smoothness, initial_bandwidth, inflation, &
      p_max, q_max, max_iterations) result(out)
      !! Select a derivative bandwidth after an IPI pilot trend fit.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      integer, intent(in) :: derivative_order !! Derivative order.
      integer, intent(in), optional :: pilot_order !! Pilot order.
      integer, intent(in), optional :: kernel_smoothness !! Kernel smoothness.
      integer, intent(in), optional :: pilot_kernel_smoothness !! Pilot kernel smoothness.
      integer, intent(in), optional :: inflation !! Inflation.
      integer, intent(in), optional :: p_max !! P max.
      integer, intent(in), optional :: q_max !! Q max.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      real(dp), intent(in), optional :: initial_bandwidth !! Initial bandwidth.
      type(esemifar_model_t) :: out
      type(esemifar_model_t) :: pilot
      type(esemifar_smooth_t) :: high_derivative
      real(dp), allocatable :: history(:)
      real(dp) :: bandwidth, old_bandwidth, older_bandwidth, inflated, d, kernel_d
      real(dp) :: c1, c2, c3, rp, moment, i2, candidate
      integer :: pilot_p, mu, pilot_mu, inflation_method, max_p, max_q
      integer :: limit, iteration, trim, k, denominator_power

      pilot_p = 1
      if (present(pilot_order)) pilot_p = pilot_order
      mu = 1
      if (present(kernel_smoothness)) mu = kernel_smoothness
      pilot_mu = 1
      if (present(pilot_kernel_smoothness)) pilot_mu = pilot_kernel_smoothness
      inflation_method = esemifar_inflation_optimal
      if (present(inflation)) inflation_method = inflation
      max_p = 0
      if (present(p_max)) max_p = p_max
      max_q = 0
      if (present(q_max)) max_q = q_max
      limit = 40
      if (present(max_iterations)) limit = max_iterations
      if (derivative_order < 1 .or. derivative_order > 2 .or. &
         (pilot_p /= 1 .and. pilot_p /= 3) .or. mu < 0 .or. mu > 3) then
         out%info = 1
         return
      end if
      pilot = esemifar_trend_fit(series, pilot_p, pilot_mu, initial_bandwidth, &
         inflation_method, max_p, max_q, max_iterations=limit)
      if (pilot%info /= 0) then
         out%info = 2
         return
      end if
      bandwidth = pilot%smoother%bandwidth
      d = pilot%farima%model%long_memory_parameter
      kernel_d = max(d, sqrt(epsilon(1.0_dp)))
      allocate(history(limit))
      history = 0.0_dp
      older_bandwidth = bandwidth
      k = derivative_order + 2
      denominator_power = 2*derivative_order + 5
      do iteration = 1, limit
         old_bandwidth = bandwidth
         inflated = derivative_inflated_bandwidth(bandwidth, d, derivative_order, &
            inflation_method)
         high_derivative = esemifar_smooth(series, k, derivative_order + 3, mu, &
            min(0.49_dp, inflated))
         trim = int(0.05_dp*real(size(series), dp))
         i2 = sum(high_derivative%estimate(max(1, trim): &
            max(1, size(series) - trim))**2)/real(max(1, size(series) - 2*trim), dp)
         rp = 2.0_dp*gamma(1.0_dp - 2.0_dp*kernel_d)* &
            sin(acos(-1.0_dp)*kernel_d)* &
            derivative_kernel_rp(derivative_order, mu, kernel_d)
         moment = derivative_kernel_moment(derivative_order, mu, k)
         c1 = factorial_real(k)**2*real(2*derivative_order + 1, dp)/4.0_dp
         c2 = 0.9_dp*(1.0_dp - 2.0_dp*d)*rp/(moment**2)
         c3 = pilot%variance_factor/max(i2, tiny(1.0_dp))
         candidate = max(c1*c2*c3, tiny(1.0_dp))** &
            (1.0_dp/(real(denominator_power, dp) - 2.0_dp*d))* &
            real(size(series), dp)**((2.0_dp*d - 1.0_dp)/ &
            (real(denominator_power, dp) - 2.0_dp*d))
         candidate = max(candidate, real(size(series), dp)** &
            (-real(denominator_power, dp)/real(denominator_power + 2, dp)))
         bandwidth = min(0.49_dp, candidate)
         history(iteration) = bandwidth
         if (iteration > 2 .and. abs(old_bandwidth - bandwidth)/bandwidth < &
            1.0_dp/real(size(series), dp)) then
            out%converged = .true.
            exit
         end if
         if (iteration > 3 .and. abs(older_bandwidth - bandwidth)/bandwidth < &
            1.0_dp/real(size(series), dp)) then
            bandwidth = 0.5_dp*(old_bandwidth + bandwidth)
            out%converged = .true.
            exit
         end if
         older_bandwidth = old_bandwidth
      end do
      out%iterations = min(iteration, limit)
      out%bandwidth_history = history(:out%iterations)
      out%smoother = esemifar_smooth(series, derivative_order, derivative_order + 1, &
         mu, bandwidth)
      out%farima = pilot%farima
      out%variance_factor = pilot%variance_factor
      out%curvature_integral = i2
      if (out%smoother%info /= 0) out%info = 3
   end function esemifar_derivative_ipi

   pure function esemifar_forecast_normal(model, horizon, levels, linear_trend, &
      exponentiate) result(out)
      !! Forecast trend plus FARIMA errors with Gaussian analytic intervals.
      type(esemifar_model_t), intent(in) :: model !! Model specification.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      real(dp), intent(in) :: levels(:) !! Levels.
      logical, intent(in), optional :: linear_trend !! Flag controlling linear trend.
      logical, intent(in), optional :: exponentiate !! Flag controlling exponentiate.
      type(esemifar_forecast_t) :: out
      real(dp), allocatable :: trend(:), error_forecast(:), psi(:), sd(:)
      real(dp) :: residual_mean, step, quantile
      integer :: h, level
      logical :: linear, expo

      linear = .true.
      if (present(linear_trend)) linear = linear_trend
      expo = .false.
      if (present(exponentiate)) expo = exponentiate
      if (model%info /= 0 .or. horizon < 1 .or. size(levels) < 1 .or. &
         any(levels <= 0.0_dp) .or. any(levels >= 1.0_dp)) then
         out%info = 1
         return
      end if
      allocate(trend(horizon))
      step = 0.0_dp
      if (linear .and. size(model%smoother%estimate) > 1) step = &
         model%smoother%estimate(size(model%smoother%estimate)) - &
         model%smoother%estimate(size(model%smoother%estimate) - 1)
      do h = 1, horizon
         trend(h) = model%smoother%estimate(size(model%smoother%estimate)) + &
            real(h, dp)*step
      end do
      residual_mean = sum(model%smoother%residuals)/ &
         real(size(model%smoother%residuals), dp)
      error_forecast = farima_point_forecast(model%smoother%residuals, model%farima, &
         horizon, residual_mean)
      out%mean = trend + error_forecast
      psi = esemifar_farima_to_ma(model%farima%model%ar, -model%farima%model%theta, &
         model%farima%model%long_memory_parameter, horizon - 1)
      allocate(sd(horizon), out%lower(horizon, size(levels)), &
         out%upper(horizon, size(levels)))
      do h = 1, horizon
         sd(h) = sqrt(model%farima%innovation_variance*sum(psi(:h)**2))
      end do
      do level = 1, size(levels)
         quantile = normal_quantile(0.5_dp + 0.5_dp*levels(level))
         out%lower(:, level) = out%mean - quantile*sd
         out%upper(:, level) = out%mean + quantile*sd
      end do
      if (expo) then
         out%mean = exp(out%mean)
         out%lower = exp(out%lower)
         out%upper = exp(out%upper)
      end if
   end function esemifar_forecast_normal

   function esemifar_forecast_bootstrap(model, horizon, levels, paths, seed, &
      linear_trend, exponentiate) result(out)
      !! Forecast with residual-resampling intervals and fixed FARIMA parameters.
      type(esemifar_model_t), intent(in) :: model !! Model specification.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      integer, intent(in) :: paths !! Paths.
      real(dp), intent(in) :: levels(:) !! Levels.
      integer, intent(in), optional :: seed !! Random-number seed.
      logical, intent(in), optional :: linear_trend !! Flag controlling linear trend.
      logical, intent(in), optional :: exponentiate !! Flag controlling exponentiate.
      type(esemifar_forecast_t) :: out
      real(dp), allocatable :: simulations(:, :), innovations(:), trend(:), centered(:)
      real(dp) :: step
      integer :: path, h, index, level
      logical :: linear, expo

      linear = .true.
      if (present(linear_trend)) linear = linear_trend
      expo = .false.
      if (present(exponentiate)) expo = exponentiate
      if (model%info /= 0 .or. horizon < 1 .or. paths < 2 .or. &
         size(levels) < 1 .or. any(levels <= 0.0_dp) .or. &
         any(levels >= 1.0_dp)) then
         out%info = 1
         return
      end if
      if (present(seed)) call set_random_seed(seed)
      centered = model%farima%residuals - sum(model%farima%residuals)/ &
         real(size(model%farima%residuals), dp)
      allocate(simulations(horizon, paths), innovations(horizon), trend(horizon))
      step = 0.0_dp
      if (linear .and. size(model%smoother%estimate) > 1) step = &
         model%smoother%estimate(size(model%smoother%estimate)) - &
         model%smoother%estimate(size(model%smoother%estimate) - 1)
      do h = 1, horizon
         trend(h) = model%smoother%estimate(size(model%smoother%estimate)) + &
            real(h, dp)*step
      end do
      do path = 1, paths
         do h = 1, horizon
            index = min(size(centered), 1 + int(random_uniform()*real(size(centered), dp)))
            innovations(h) = centered(index)
         end do
         simulations(:, path) = trend + farima_future_path(model%smoother%residuals, &
            model%farima, innovations)
      end do
      out%mean = trend + farima_point_forecast(model%smoother%residuals, model%farima, &
         horizon, sum(model%smoother%residuals)/real(size(model%smoother%residuals), dp))
      allocate(out%lower(horizon, size(levels)), out%upper(horizon, size(levels)))
      do h = 1, horizon
         call sort(simulations(h, :))
         do level = 1, size(levels)
            out%lower(h, level) = quantile(simulations(h, :), &
               0.5_dp - 0.5_dp*levels(level))
            out%upper(h, level) = quantile(simulations(h, :), &
               0.5_dp + 0.5_dp*levels(level))
         end do
      end do
      if (expo) then
         out%mean = exp(out%mean)
         out%lower = exp(out%lower)
         out%upper = exp(out%upper)
      end if
   end function esemifar_forecast_bootstrap

   function esemifar_forecast_bootstrap_advanced(model, horizon, levels, paths, &
      seed, linear_trend, exponentiate, burn_in, max_iterations) result(out)
      !! Forecast with FARIMA-refitted predictive-root bootstrap intervals.
      type(esemifar_model_t), intent(in) :: model !! Model specification.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      integer, intent(in) :: paths !! Paths.
      real(dp), intent(in) :: levels(:) !! Levels.
      integer, intent(in), optional :: seed !! Random-number seed.
      integer, intent(in), optional :: burn_in !! Burn in.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      logical, intent(in), optional :: linear_trend !! Flag controlling linear trend.
      logical, intent(in), optional :: exponentiate !! Flag controlling exponentiate.
      type(esemifar_forecast_t) :: out
      type(fracdiff_simulation_t) :: simulation
      type(fracdiff_fit_t) :: refit
      real(dp), allocatable :: roots(:, :), innovations(:), future_innovations(:)
      real(dp), allocatable :: centered(:), trend(:), base(:), future(:), refitted(:)
      real(dp) :: step, simulated_mean
      integer :: path, h, level, selected_burn_in, selected_iterations, required
      logical :: linear, expo

      linear = .true.
      if (present(linear_trend)) linear = linear_trend
      expo = .false.
      if (present(exponentiate)) expo = exponentiate
      selected_burn_in = 5000
      if (present(burn_in)) selected_burn_in = burn_in
      selected_iterations = 120
      if (present(max_iterations)) selected_iterations = max_iterations
      if (model%info /= 0 .or. horizon < 1 .or. paths < 2 .or. &
         size(levels) < 1 .or. any(levels <= 0.0_dp) .or. &
         any(levels >= 1.0_dp) .or. selected_burn_in < 0 .or. &
         selected_iterations < 1) then
         out%info = 1
         return
      end if
      if (present(seed)) call set_random_seed(seed)
      centered = model%farima%residuals - sum(model%farima%residuals)/ &
         real(size(model%farima%residuals), dp)
      required = size(model%smoother%residuals) + selected_burn_in + &
         size(model%farima%model%theta)
      allocate(roots(horizon, paths), innovations(required), &
         future_innovations(horizon), trend(horizon))
      step = 0.0_dp
      if (linear .and. size(model%smoother%estimate) > 1) step = &
         model%smoother%estimate(size(model%smoother%estimate)) - &
         model%smoother%estimate(size(model%smoother%estimate) - 1)
      do h = 1, horizon
         trend(h) = model%smoother%estimate(size(model%smoother%estimate)) + &
            real(h, dp)*step
      end do
      base = farima_point_forecast(model%smoother%residuals, model%farima, &
         horizon, sum(model%smoother%residuals)/ &
         real(size(model%smoother%residuals), dp))
      do path = 1, paths
         call resample_centered(centered, innovations)
         simulation = fracdiff_simulate_from_innovations( &
            size(model%smoother%residuals), model%farima%model%ar, &
            model%farima%model%theta, model%farima%model%long_memory_parameter, &
            innovations, burn_in=selected_burn_in)
         if (simulation%info /= 0) then
            out%info = 2
            return
         end if
         simulated_mean = sum(simulation%series)/real(size(simulation%series), dp)
         refit = fracdiff_fit(simulation%series, model%farima%model%ar, &
            model%farima%model%theta, model%farima%model%long_memory_parameter, &
            truncation=min(model%farima%truncation, size(simulation%series)), &
            max_iterations=selected_iterations)
         if (refit%info /= 0 .and. .not. allocated(refit%residuals)) then
            out%info = 3
            return
         end if
         refitted = farima_point_forecast(model%smoother%residuals, refit, horizon, &
            simulated_mean)
         call resample_centered(centered, future_innovations)
         future = farima_future_path(model%smoother%residuals, model%farima, &
            future_innovations)
         roots(:, path) = future - refitted
      end do
      out%mean = trend + base
      allocate(out%lower(horizon, size(levels)), out%upper(horizon, size(levels)))
      do h = 1, horizon
         call sort(roots(h, :))
         do level = 1, size(levels)
            out%lower(h, level) = trend(h) + base(h) + &
               quantile(roots(h, :), 0.5_dp - 0.5_dp*levels(level))
            out%upper(h, level) = trend(h) + base(h) + &
               quantile(roots(h, :), 0.5_dp + 0.5_dp*levels(level))
         end do
      end do
      if (expo) then
         out%mean = exp(out%mean)
         out%lower = exp(out%lower)
         out%upper = exp(out%upper)
      end if
   end function esemifar_forecast_bootstrap_advanced

   subroutine resample_centered(sample, result)
      !! Draw a replacement sample from a centered innovation vector.
      real(dp), intent(in) :: sample(:) !! Sample.
      real(dp), intent(out) :: result(:) !! Result.
      integer :: i, index

      do i = 1, size(result)
         index = min(size(sample), 1 + int(random_uniform()*real(size(sample), dp)))
         result(i) = sample(index)
      end do
   end subroutine resample_centered

   pure function farima_point_forecast(observations, fit, horizon, mean_value) result(forecast)
      !! Forecast from esemifar's truncated infinite-AR representation.
      real(dp), intent(in) :: observations(:) !! Observed time-series values.
      real(dp), intent(in) :: mean_value !! Mean value.
      type(fracdiff_fit_t), intent(in) :: fit !! Previously fitted model.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      real(dp), allocatable :: forecast(:), history(:), ar_weights(:)
      integer :: n, h, lag

      n = size(observations)
      ar_weights = esemifar_farima_to_ar(fit%model%ar, -fit%model%theta, &
         fit%model%long_memory_parameter, n + horizon - 1)
      allocate(history(n + horizon), forecast(horizon))
      history(:n) = observations - mean_value
      do h = 1, horizon
         history(n + h) = 0.0_dp
         do lag = 1, n + h - 1
            history(n + h) = history(n + h) + &
               ar_weights(lag + 1)*history(n + h - lag)
         end do
         forecast(h) = history(n + h) + mean_value
      end do
   end function farima_point_forecast

   pure function farima_future_path(observations, fit, innovations) result(path)
      !! Generate future observations from supplied FARIMA innovations.
      real(dp), intent(in) :: observations(:) !! Observed time-series values.
      real(dp), intent(in) :: innovations(:) !! Model innovations.
      type(fracdiff_fit_t), intent(in) :: fit !! Previously fitted model.
      real(dp), allocatable :: path(:), history(:), ar_weights(:)
      real(dp) :: mean_value
      integer :: n, h, lag

      n = size(observations)
      mean_value = sum(observations)/real(n, dp)
      ar_weights = esemifar_farima_to_ar(fit%model%ar, -fit%model%theta, &
         fit%model%long_memory_parameter, n + size(innovations) - 1)
      allocate(history(n + size(innovations)), path(size(innovations)))
      history(:n) = observations - mean_value
      do h = 1, size(innovations)
         history(n + h) = innovations(h)
         do lag = 1, n + h - 1
            history(n + h) = history(n + h) + &
               ar_weights(lag + 1)*history(n + h - lag)
         end do
         path(h) = history(n + h) + mean_value
      end do
   end function farima_future_path

   pure real(dp) function farima_variance_factor(fit) result(value)
      !! Compute the low-frequency FARIMA variance factor.
      type(fracdiff_fit_t), intent(in) :: fit !! Previously fitted model.
      real(dp) :: ar_sum, ma_sum

      ar_sum = sum(fit%model%ar)
      ma_sum = sum(fit%model%theta)
      value = ((1.0_dp - ma_sum)/(1.0_dp - ar_sum))**2* &
         fit%innovation_variance/(2.0_dp*acos(-1.0_dp))
   end function farima_variance_factor

   pure real(dp) function inflated_bandwidth(bandwidth, d, polynomial, method) result(value)
      !! Apply one of esemifar's optimal, naive, or stable inflation rates.
      real(dp), intent(in) :: bandwidth !! Smoothing or spectral bandwidth.
      real(dp), intent(in) :: d !! Fractional-differencing parameter or differencing order.
      integer, intent(in) :: polynomial !! Polynomial.
      integer, intent(in) :: method !! Algorithm or estimation method.
      real(dp) :: numerator, denominator

      if (method == esemifar_inflation_stable) then
         value = sqrt(bandwidth)
         return
      end if
      numerator = real(2*polynomial + 3, dp) - 2.0_dp*d
      denominator = real(2*polynomial + 5, dp) - 2.0_dp*d
      if (method == esemifar_inflation_naive) denominator = &
         real(2*polynomial + 7, dp) - 2.0_dp*d
      value = bandwidth**(numerator/denominator)
   end function inflated_bandwidth

   pure real(dp) function kernel_rp(polynomial, mu, d) result(value)
      !! Evaluate the squared equivalent-kernel long-memory constant.
      integer, intent(in) :: polynomial !! Polynomial.
      integer, intent(in) :: mu !! Mu.
      real(dp), intent(in) :: d !! Fractional-differencing parameter or differencing order.
      real(dp) :: coefficient(0:4)
      integer :: i, j, degree

      coefficient = 0.0_dp
      if (polynomial == 1) then
         degree = mu
         do i = 0, mu
            coefficient(i) = (-1.0_dp)**i*binomial(mu, i)
         end do
      else
         degree = mu + 2
         select case (mu)
         case (0)
            coefficient(:1) = [3.0_dp, -5.0_dp]
            coefficient = 0.375_dp*coefficient
            degree = 1
         case (1)
            coefficient(:2) = [3.0_dp, -10.0_dp, 7.0_dp]
         case (2)
            coefficient(:3) = [1.0_dp, -5.0_dp, 7.0_dp, -3.0_dp]
         case default
            coefficient = [3.0_dp, -20.0_dp, 42.0_dp, -36.0_dp, 11.0_dp]
         end select
      end if
      value = 0.0_dp
      do i = 0, degree
         do j = 0, degree
            value = value + coefficient(i)*coefficient(j)* &
               esemifar_kdf(2*i, 2*j, d)
         end do
      end do
   end function kernel_rp

   pure real(dp) function kernel_moment(polynomial, mu, power) result(value)
      !! Evaluate the equivalent kernel's requested even moment.
      integer, intent(in) :: polynomial !! Polynomial.
      integer, intent(in) :: mu !! Mu.
      integer, intent(in) :: power !! Power.
      real(dp) :: coefficient(0:4)
      integer :: i, degree

      coefficient = 0.0_dp
      if (polynomial == 1) then
         degree = mu
         do i = 0, mu
            coefficient(i) = (-1.0_dp)**i*binomial(mu, i)
         end do
      else
         select case (mu)
         case (0)
            coefficient(:1) = 0.375_dp*[3.0_dp, -5.0_dp]
            degree = 1
         case (1)
            coefficient(:2) = [3.0_dp, -10.0_dp, 7.0_dp]
            degree = 2
         case (2)
            coefficient(:3) = [1.0_dp, -5.0_dp, 7.0_dp, -3.0_dp]
            degree = 3
         case default
            coefficient = [3.0_dp, -20.0_dp, 42.0_dp, -36.0_dp, 11.0_dp]
            degree = 4
         end select
      end if
      value = 0.0_dp
      do i = 0, degree
         if (mod(power, 2) == 0) value = value + &
            2.0_dp*coefficient(i)/real(power + 2*i + 1, dp)
      end do
   end function kernel_moment

   pure real(dp) function derivative_inflated_bandwidth(bandwidth, d, &
      derivative, method) result(value)
      !! Apply esemifar's derivative pilot-bandwidth inflation rate.
      real(dp), intent(in) :: bandwidth !! Smoothing or spectral bandwidth.
      real(dp), intent(in) :: d !! Fractional-differencing parameter or differencing order.
      integer, intent(in) :: derivative !! Derivative.
      integer, intent(in) :: method !! Algorithm or estimation method.
      real(dp) :: numerator, denominator

      if (method == esemifar_inflation_stable) then
         value = sqrt(bandwidth)
         return
      end if
      numerator = real(2*derivative + 5, dp) - 2.0_dp*d
      denominator = real(2*derivative + 7, dp) - 2.0_dp*d
      if (method == esemifar_inflation_naive) denominator = &
         real(2*derivative + 9, dp) - 2.0_dp*d
      value = bandwidth**(numerator/denominator)
   end function derivative_inflated_bandwidth

   pure real(dp) function derivative_kernel_rp(derivative, mu, d) result(value)
      !! Evaluate a derivative equivalent-kernel long-memory constant.
      integer, intent(in) :: derivative !! Derivative.
      integer, intent(in) :: mu !! Mu.
      real(dp), intent(in) :: d !! Fractional-differencing parameter or differencing order.
      real(dp) :: coefficient(0:8)
      integer :: i, j, degree

      call derivative_kernel_coefficients(derivative, mu, coefficient, degree)
      value = 0.0_dp
      do i = 0, degree
         do j = 0, degree
            value = value + coefficient(i)*coefficient(j)*esemifar_kdf(i, j, d)
         end do
      end do
   end function derivative_kernel_rp

   pure real(dp) function derivative_kernel_moment(derivative, mu, power) result(value)
      !! Evaluate a derivative equivalent kernel's requested moment.
      integer, intent(in) :: derivative !! Derivative.
      integer, intent(in) :: mu !! Mu.
      integer, intent(in) :: power !! Power.
      real(dp) :: coefficient(0:8)
      integer :: i, degree

      call derivative_kernel_coefficients(derivative, mu, coefficient, degree)
      value = 0.0_dp
      do i = 0, degree
         if (mod(power + i, 2) == 0) value = value + &
            2.0_dp*coefficient(i)/real(power + i + 1, dp)
      end do
   end function derivative_kernel_moment

   pure subroutine derivative_kernel_coefficients(derivative, mu, coefficient, degree)
      !! Return esemifar's first- or second-derivative equivalent kernel.
      integer, intent(in) :: derivative !! Derivative.
      integer, intent(in) :: mu !! Mu.
      real(dp), intent(out) :: coefficient(0:8) !! Coefficient.
      integer, intent(out) :: degree !! Degree.

      coefficient = 0.0_dp
      if (derivative == 1) then
         select case (mu)
         case (0)
            coefficient(1) = 1.0_dp
            degree = 1
         case (1)
            coefficient(1) = 1.0_dp
            coefficient(3) = -1.0_dp
            degree = 3
         case (2)
            coefficient(1) = -1.0_dp
            coefficient(3) = 2.0_dp
            coefficient(5) = -1.0_dp
            degree = 5
         case default
            coefficient(1) = -1.0_dp
            coefficient(3) = 3.0_dp
            coefficient(5) = -3.0_dp
            coefficient(7) = 1.0_dp
            degree = 7
         end select
      else
         select case (mu)
         case (0)
            coefficient(0) = 1.0_dp
            coefficient(2) = -3.0_dp
            degree = 2
         case (1)
            coefficient(0) = -1.0_dp
            coefficient(2) = 6.0_dp
            coefficient(4) = -5.0_dp
            degree = 4
         case (2)
            coefficient(0) = -1.0_dp
            coefficient(2) = 9.0_dp
            coefficient(4) = -15.0_dp
            coefficient(6) = 7.0_dp
            degree = 6
         case default
            coefficient(0) = -1.0_dp
            coefficient(2) = 12.0_dp
            coefficient(4) = -30.0_dp
            coefficient(6) = 28.0_dp
            coefficient(8) = -9.0_dp
            degree = 8
         end select
      end if
   end subroutine derivative_kernel_coefficients

   pure elemental real(dp) function binomial(n, k) result(value)
      !! Return a binomial coefficient in double precision.
      integer, intent(in) :: n !! Number of observations or elements.
      integer, intent(in) :: k !! K.
      integer :: i

      if (k < 0 .or. k > n) then
         value = 0.0_dp
         return
      end if
      value = 1.0_dp
      do i = 1, min(k, n - k)
         value = value*real(n - i + 1, dp)/real(i, dp)
      end do
   end function binomial

   pure elemental real(dp) function factorial_real(n) result(value)
      !! Return a factorial in double precision.
      integer, intent(in) :: n !! Number of observations or elements.
      integer :: i

      value = 1.0_dp
      do i = 2, n
         value = value*real(i, dp)
      end do
   end function factorial_real

end module esemifar_mod
