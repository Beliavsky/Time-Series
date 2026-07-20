! SPDX-License-Identifier: GPL-3.0-only
! SPDX-FileComment: Algorithms translated from the R garma package.
! Gegenbauer ARMA algorithms translated from garma.
module garma_mod
   use kind_mod, only: dp
   use fourier_mod, only: fft_transform
   use linalg_mod, only: invert_matrix
   use optimization_mod, only: optimization_result_t, &
      nelder_mead_minimize, finite_difference_hessian
   use stats_mod, only: ols_fit, sorted, variance_of => variance
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   private

   integer, parameter, public :: garma_semiparametric_gsp = 1
   integer, parameter, public :: garma_semiparametric_lpr = 2
   integer, parameter, public :: garma_method_css = 1
   integer, parameter, public :: garma_method_whittle = 2
   integer, parameter, public :: garma_method_wll = 3

   type, public :: garma_periodogram_t
      ! Fourier frequencies and raw detrended periodogram ordinates.
      real(dp), allocatable :: frequency(:), spectrum(:)
      integer :: info = 0
   end type garma_periodogram_t

   type, public :: garma_factors_t
      ! Gegenbauer pole locations, periods, exponents, and status.
      real(dp), allocatable :: u(:), frequency(:), period(:), d(:)
      integer :: bandwidth = 0
      integer :: info = 0
   end type garma_factors_t

   type, public :: garma_fit_t
      ! Fitted GARMA coefficients, inference, residuals, and fit criteria.
      real(dp), allocatable :: u(:), d(:), ar(:), ma(:)
      real(dp), allocatable :: covariance(:, :), standard_error(:)
      real(dp), allocatable :: fitted(:), residuals(:)
      real(dp) :: innovation_variance = 0.0_dp
      real(dp) :: objective = huge(1.0_dp)
      real(dp) :: log_likelihood = -huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      integer :: method = garma_method_whittle
      integer :: iterations = 0
      integer :: info = 0
      logical :: estimated_frequencies = .false.
      logical :: converged = .false.
   end type garma_fit_t

   type, public :: garma_forecast_t
      ! GARMA point forecasts and status.
      real(dp), allocatable :: mean(:)
      integer :: info = 0
   end type garma_forecast_t

   type, public :: garma_gof_t
      ! Sequential Bartlett periodogram test statistics and probabilities.
      real(dp), allocatable :: statistic(:), p_value(:)
      real(dp) :: minimum_p_value = 1.0_dp
      real(dp) :: minimum_frequency = 0.0_dp
      real(dp) :: minimum_period = 0.0_dp
      integer :: info = 0
   end type garma_gof_t

   type, public :: garma_regression_fit_t
      ! GARMA dynamics combined with deterministic and external regression terms.
      type(garma_fit_t) :: dynamic
      real(dp), allocatable :: regression_coefficients(:), regression_standard_error(:)
      real(dp), allocatable :: regression_covariance(:, :), regression_residuals(:)
      real(dp), allocatable :: differenced_series(:), fitted(:), residuals(:)
      integer :: difference_order = 0
      integer :: external_regressors = 0
      logical :: includes_mean = .false.
      logical :: includes_drift = .false.
      integer :: info = 0
   end type garma_regression_fit_t

   type, public :: garma_accuracy_t
      ! Standard point-forecast accuracy measures.
      real(dp) :: mean_error = 0.0_dp
      real(dp) :: root_mean_squared_error = 0.0_dp
      real(dp) :: mean_absolute_error = 0.0_dp
      real(dp) :: mean_percentage_error = 0.0_dp
      real(dp) :: mean_absolute_percentage_error = 0.0_dp
      real(dp) :: mean_absolute_scaled_error = 0.0_dp
      real(dp) :: residual_acf1 = 0.0_dp
      integer :: observations = 0
      integer :: info = 0
   end type garma_accuracy_t

   public :: garma_gegenbauer_coefficients, garma_periodogram
   public :: garma_semiparametric, garma_spectral_inverse
   public :: garma_css_residuals, garma_extract_arma
   public :: garma_whittle_objective, garma_wll_objective, garma_fit
   public :: garma_fitted, garma_forecast, garma_goodness_of_fit
   public :: garma_regression_fit, garma_regression_forecast, garma_accuracy

contains

   pure function garma_gegenbauer_coefficients(n, d, u) result(coefficients)
      !! Expand the Gegenbauer filter through n coefficients.
      integer, intent(in) :: n !! Number of observations or elements.
      real(dp), intent(in) :: d !! Fractional-differencing parameter or differencing order.
      real(dp), intent(in) :: u !! Input vector or random variate.
      real(dp), allocatable :: coefficients(:)
      integer :: j

      if (n < 1 .or. abs(u) > 1.0_dp .or. .not. ieee_is_finite(d)) then
         allocate(coefficients(0))
         return
      end if
      allocate(coefficients(n))
      coefficients(1) = 1.0_dp
      if (n >= 2) coefficients(2) = 2.0_dp*d*u
      if (n >= 3) coefficients(3) = 2.0_dp*d*(d + 1.0_dp)*u*u - d
      do j = 4, n
         coefficients(j) = (2.0_dp*u*(real(j - 2, dp) + d)*coefficients(j - 1) - &
            (real(j - 3, dp) + 2.0_dp*d)*coefficients(j - 2))/real(j - 1, dp)
      end do
   end function garma_gegenbauer_coefficients

   pure function garma_periodogram(series) result(out)
      !! Compute garma's linearly detrended, unsmoothed raw periodogram.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      type(garma_periodogram_t) :: out
      complex(dp), allocatable :: transformed(:), padded(:)
      real(dp), allocatable :: detrended(:)
      real(dp) :: mean_x, mean_t, slope, denominator
      integer :: n, count, i

      n = size(series)
      if (n < 4 .or. .not. all(ieee_is_finite(series))) then
         out%info = 1
         return
      end if
      mean_x = sum(series)/real(n, dp)
      mean_t = 0.5_dp*real(n + 1, dp)
      denominator = real(n*(n*n - 1), dp)/12.0_dp
      slope = sum(([(real(i, dp), i=1,n)] - mean_t)*(series - mean_x))/denominator
      detrended = series - mean_x - slope*([(real(i, dp), i=1,n)] - mean_t)
      allocate(padded(n))
      padded = cmplx(detrended, 0.0_dp, dp)
      transformed = discrete_transform(padded)
      count = n/2
      allocate(out%frequency(count), out%spectrum(count))
      do i = 1, count
         out%frequency(i) = real(i, dp)/real(n, dp)
         out%spectrum(i) = abs(transformed(i + 1))**2/real(n, dp)
      end do
   end function garma_periodogram

   pure function garma_semiparametric(series, periods, factor_count, alpha, &
      method, peak_bandwidth) result(out)
      !! Estimate Gegenbauer poles and exponents by GSP or log-periodogram regression.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      real(dp), intent(in), optional :: periods(:) !! Periods.
      real(dp), intent(in), optional :: alpha !! Significance, smoothing, or model coefficient.
      integer, intent(in), optional :: factor_count !! Number of factor.
      integer, intent(in), optional :: method !! Algorithm or estimation method.
      integer, intent(in), optional :: peak_bandwidth !! Peak bandwidth.
      type(garma_factors_t) :: out
      type(garma_periodogram_t) :: pgram
      real(dp), allocatable :: work(:)
      real(dp) :: selected_alpha
      integer :: count, selected_method, bandwidth, factor, index, left, right

      selected_alpha = 0.8_dp
      if (present(alpha)) selected_alpha = alpha
      selected_method = garma_semiparametric_gsp
      if (present(method)) selected_method = method
      bandwidth = 1
      if (present(peak_bandwidth)) bandwidth = peak_bandwidth
      count = 1
      if (present(factor_count)) count = factor_count
      if (present(periods)) count = size(periods)
      if (size(series) < 4 .or. count < 1 .or. selected_alpha <= 0.0_dp .or. &
         selected_alpha >= 1.0_dp .or. bandwidth < 0 .or. &
         (selected_method /= garma_semiparametric_gsp .and. &
         selected_method /= garma_semiparametric_lpr)) then
         out%info = 1
         return
      end if
      pgram = garma_periodogram(series)
      if (pgram%info /= 0) then
         out%info = 2
         return
      end if
      allocate(out%u(count), out%frequency(count), out%period(count), out%d(count))
      out%bandwidth = floor((real(size(series), dp)/2.0_dp)**selected_alpha)
      if (out%bandwidth < 2 .or. out%bandwidth > size(pgram%spectrum)) then
         out%info = 3
         return
      end if
      if (present(periods)) then
         if (any(periods <= 2.0_dp)) then
            out%info = 1
            return
         end if
         out%period = periods
         out%frequency = 1.0_dp/periods
      else
         work = pgram%spectrum
         do factor = 1, count
            index = maxloc(work, dim=1)
            out%frequency(factor) = pgram%frequency(index)
            out%period(factor) = 1.0_dp/out%frequency(factor)
            left = max(1, index - bandwidth)
            right = min(size(work), index + bandwidth)
            work(left:right) = 0.0_dp
         end do
      end if
      out%u = cos(2.0_dp*acos(-1.0_dp)*out%frequency)
      do factor = 1, count
         index = minloc(abs(pgram%frequency - out%frequency(factor)), dim=1)
         if (selected_method == garma_semiparametric_gsp) then
            out%d(factor) = gsp_estimate(pgram, index, out%bandwidth)
         else
            out%d(factor) = lpr_estimate(pgram, index, out%bandwidth)
         end if
      end do
   end function garma_semiparametric

   pure function garma_spectral_inverse(frequency, u, d, ar, ma) result(inverse)
      !! Evaluate the GARMA spectrum denominator without the innovation scale.
      real(dp), intent(in) :: frequency(:) !! Number of observations per seasonal cycle.
      real(dp), intent(in) :: u(:) !! Input vector or random variate.
      real(dp), intent(in) :: d(:) !! Fractional-differencing parameter or differencing order.
      real(dp), intent(in) :: ar(:) !! Autoregressive coefficients.
      real(dp), intent(in) :: ma(:) !! Moving-average coefficients.
      real(dp), allocatable :: inverse(:)
      real(dp) :: omega, ar_real, ar_imag, ma_real, ma_imag
      integer :: i, j

      if (size(u) /= size(d) .or. any(abs(u) > 1.0_dp)) then
         allocate(inverse(0))
         return
      end if
      allocate(inverse(size(frequency)))
      do i = 1, size(frequency)
         omega = 2.0_dp*acos(-1.0_dp)*frequency(i)
         ar_real = 1.0_dp
         ar_imag = 0.0_dp
         do j = 1, size(ar)
            ar_real = ar_real - ar(j)*cos(real(j, dp)*omega)
            ar_imag = ar_imag - ar(j)*sin(real(j, dp)*omega)
         end do
         ma_real = 1.0_dp
         ma_imag = 0.0_dp
         do j = 1, size(ma)
            ma_real = ma_real - ma(j)*cos(real(j, dp)*omega)
            ma_imag = ma_imag - ma(j)*sin(real(j, dp)*omega)
         end do
         inverse(i) = (ar_real**2 + ar_imag**2)/(ma_real**2 + ma_imag**2)
         do j = 1, size(u)
            inverse(i) = inverse(i)*(4.0_dp*(cos(omega) - u(j))**2)**d(j)
         end do
      end do
   end function garma_spectral_inverse

   pure function garma_css_residuals(series, u, d, ar, ma) result(residuals)
      !! Filter observations through the short-memory and Gegenbauer denominators.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      real(dp), intent(in) :: u(:) !! Input vector or random variate.
      real(dp), intent(in) :: d(:) !! Fractional-differencing parameter or differencing order.
      real(dp), intent(in) :: ar(:) !! Autoregressive coefficients.
      real(dp), intent(in) :: ma(:) !! Moving-average coefficients.
      real(dp), allocatable :: residuals(:)
      real(dp), allocatable :: numerator(:), denominator(:), work(:)
      integer :: factor

      if (size(u) /= size(d) .or. size(series) < 1) then
         allocate(residuals(0))
         return
      end if
      allocate(numerator(0:size(ar)), denominator(0:size(ma)))
      numerator = 0.0_dp
      denominator = 0.0_dp
      numerator(0) = 1.0_dp
      denominator(0) = 1.0_dp
      if (size(ar) > 0) numerator(1:) = -ar
      if (size(ma) > 0) denominator(1:) = ma
      work = rational_filter(series, numerator, denominator)
      do factor = 1, size(u)
         denominator = garma_gegenbauer_coefficients(size(series), d(factor), u(factor))
         work = rational_filter(work, [1.0_dp], denominator)
      end do
      residuals = work
   end function garma_css_residuals

   pure function garma_extract_arma(series, u, d) result(short_memory)
      !! Remove all Gegenbauer long-memory factors from a series.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      real(dp), intent(in) :: u(:) !! Input vector or random variate.
      real(dp), intent(in) :: d(:) !! Fractional-differencing parameter or differencing order.
      real(dp), allocatable :: short_memory(:), denominator(:)
      integer :: factor

      if (size(u) /= size(d)) then
         allocate(short_memory(0))
         return
      end if
      short_memory = series
      do factor = 1, size(u)
         denominator = garma_gegenbauer_coefficients(size(series), d(factor), u(factor))
         short_memory = rational_filter(short_memory, [1.0_dp], denominator)
      end do
   end function garma_extract_arma

   pure real(dp) function garma_whittle_objective(series, u, d, ar, ma, &
      include_log_term) result(value)
      !! Evaluate garma's short or full Whittle spectral objective.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      real(dp), intent(in) :: u(:) !! Input vector or random variate.
      real(dp), intent(in) :: d(:) !! Fractional-differencing parameter or differencing order.
      real(dp), intent(in) :: ar(:) !! Autoregressive coefficients.
      real(dp), intent(in) :: ma(:) !! Moving-average coefficients.
      logical, intent(in), optional :: include_log_term !! Whether to include the log term.
      type(garma_periodogram_t) :: pgram
      real(dp), allocatable :: inverse(:)
      logical :: full

      full = .false.
      if (present(include_log_term)) full = include_log_term
      pgram = garma_periodogram(series)
      inverse = garma_spectral_inverse(pgram%frequency, u, d, ar, ma)
      if (pgram%info /= 0 .or. size(inverse) == 0 .or. any(inverse < 0.0_dp) .or. &
         .not. all(ieee_is_finite(inverse))) then
         value = huge(1.0_dp)/100.0_dp
         return
      end if
      value = sum(pgram%spectrum*inverse)/(2.0_dp*acos(-1.0_dp)* &
         real(size(inverse), dp))
      if (full) value = value - sum(log(pack(inverse, inverse > 0.0_dp)))
   end function garma_whittle_objective

   pure real(dp) function garma_wll_objective(series, u, d, ar, ma, variance) result(value)
      !! Evaluate Hunt-Peiris-Weber's Whittle-like-log objective.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      real(dp), intent(in) :: u(:) !! Input vector or random variate.
      real(dp), intent(in) :: d(:) !! Fractional-differencing parameter or differencing order.
      real(dp), intent(in) :: ar(:) !! Autoregressive coefficients.
      real(dp), intent(in) :: ma(:) !! Moving-average coefficients.
      real(dp), intent(in) :: variance !! Variance value or matrix.
      type(garma_periodogram_t) :: pgram
      real(dp), allocatable :: inverse(:), ratio(:)

      pgram = garma_periodogram(series)
      inverse = garma_spectral_inverse(pgram%frequency, u, d, ar, ma)
      if (pgram%info /= 0 .or. variance <= 0.0_dp .or. size(inverse) == 0) then
         value = huge(1.0_dp)/100.0_dp
         return
      end if
      ratio = pgram%spectrum*(2.0_dp*acos(-1.0_dp)/variance)*inverse
      if (count(ratio > 0.0_dp .and. ieee_is_finite(ratio)) == 0) then
         value = huge(1.0_dp)/100.0_dp
      else
         value = sum(log(pack(ratio, ratio > 0.0_dp .and. ieee_is_finite(ratio)))**2)
      end if
   end function garma_wll_objective

   pure function garma_fit(series, initial_u, initial_d, initial_ar, initial_ma, &
      method, estimate_frequencies, d_limits, max_iterations, tolerance) result(out)
      !! Fit a GARMA model by CSS, Whittle, or WLL with derivative-free optimization.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      real(dp), intent(in) :: initial_u(:) !! Initial u.
      real(dp), intent(in) :: initial_d(:) !! Initial d.
      real(dp), intent(in) :: initial_ar(:) !! Initial autoregressive.
      real(dp), intent(in) :: initial_ma(:) !! Initial moving-average.
      integer, intent(in), optional :: method !! Algorithm or estimation method.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      logical, intent(in), optional :: estimate_frequencies !! Whether to estimate the frequencies.
      real(dp), intent(in), optional :: d_limits(:) !! D limits.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      type(garma_fit_t) :: out
      type(optimization_result_t) :: optimized
      real(dp), allocatable :: initial(:), covariance(:, :)
      real(dp) :: lower_d, upper_d, selected_tolerance, variance
      integer :: selected_method, limit, count, status, i
      logical :: estimate_u

      selected_method = garma_method_whittle
      if (present(method)) selected_method = method
      estimate_u = .false.
      if (present(estimate_frequencies)) estimate_u = estimate_frequencies
      lower_d = 0.0_dp
      upper_d = 0.5_dp
      if (present(d_limits)) then
         if (size(d_limits) == 2) then
            lower_d = minval(d_limits)
            upper_d = maxval(d_limits)
         end if
      end if
      limit = 2000
      if (present(max_iterations)) limit = max_iterations
      selected_tolerance = 1.0e-7_dp
      if (present(tolerance)) selected_tolerance = tolerance
      if (size(series) < 16 .or. size(initial_u) /= size(initial_d) .or. &
         any(abs(initial_u) >= 1.0_dp) .or. lower_d >= upper_d .or. &
         any(initial_d < lower_d) .or. any(initial_d > upper_d) .or. &
         selected_method < garma_method_css .or. selected_method > garma_method_wll) then
         out%info = 1
         return
      end if
      variance = variance_of(series)
      initial = encode_fit(initial_u, initial_d, initial_ar, initial_ma, estimate_u, &
         lower_d, upper_d, variance, selected_method == garma_method_wll)
      optimized = nelder_mead_minimize(objective, initial, limit, selected_tolerance)
      call decode_fit(optimized%parameters, initial_u, estimate_u, lower_d, upper_d, &
         size(initial_ar), size(initial_ma), selected_method == garma_method_wll, &
         out%u, out%d, out%ar, out%ma, variance)
      out%estimated_frequencies = estimate_u
      out%method = selected_method
      out%objective = optimized%objective
      out%iterations = optimized%iterations
      out%converged = optimized%converged
      out%info = optimized%info
      if (out%info == 4 .and. ieee_is_finite(out%objective)) out%info = 0
      if (selected_method == garma_method_css) then
         out%residuals = garma_css_residuals(series, out%u, out%d, out%ar, out%ma)
         out%innovation_variance = sum(out%residuals**2)/real(size(series), dp)
      else if (selected_method == garma_method_whittle) then
         out%innovation_variance = garma_whittle_objective(series, out%u, out%d, &
            out%ar, out%ma)
         out%residuals = garma_css_residuals(series, out%u, out%d, out%ar, out%ma)
      else
         out%innovation_variance = variance*exp(0.5772156649015329_dp)/ &
            (2.0_dp*acos(-1.0_dp))
         out%residuals = garma_css_residuals(series, out%u, out%d, out%ar, out%ma)
      end if
      out%fitted = series - out%residuals
      out%log_likelihood = -0.5_dp*real(size(series), dp)* &
         log(2.0_dp*acos(-1.0_dp)*max(out%innovation_variance, tiny(1.0_dp))) - &
         0.5_dp*sum(out%residuals**2)/max(out%innovation_variance, tiny(1.0_dp))
      count = size(initial_d) + size(initial_ar) + size(initial_ma) + &
         merge(size(initial_u), 0, estimate_u)
      out%aic = -2.0_dp*out%log_likelihood + 2.0_dp*real(count + 1, dp)
      if (selected_method == garma_method_css) then
         call css_covariance(out, series, estimate_u, covariance, status)
      else if (selected_method == garma_method_whittle) then
         call spectral_covariance(out, series, estimate_u, covariance, status)
      else
         allocate(covariance(count, count))
         covariance = 0.0_dp
         status = 0
      end if
      out%covariance = covariance
      allocate(out%standard_error(size(covariance, 1)))
      do i = 1, size(out%standard_error)
         out%standard_error(i) = sqrt(max(0.0_dp, covariance(i, i)))
      end do
      if (selected_method == garma_method_wll) then
         call set_wll_standard_errors(out, series, estimate_u)
      end if
      if (status /= 0 .and. selected_method /= garma_method_wll) out%info = max(out%info, 2)

   contains

      pure function objective(parameters) result(value)
         !! Decode and evaluate the selected GARMA fit objective.
         real(dp), intent(in) :: parameters(:) !! Model parameter values.
         real(dp) :: value, candidate_variance
         real(dp), allocatable :: u(:), d(:), ar(:), ma(:), residual(:)

         call decode_fit(parameters, initial_u, estimate_u, lower_d, upper_d, &
            size(initial_ar), size(initial_ma), selected_method == garma_method_wll, &
            u, d, ar, ma, candidate_variance)
         if (selected_method == garma_method_css) then
            residual = garma_css_residuals(series, u, d, ar, ma)
            value = sum(residual**2)
         else if (selected_method == garma_method_whittle) then
            value = garma_whittle_objective(series, u, d, ar, ma)
         else
            value = garma_wll_objective(series, u, d, ar, ma, candidate_variance)
         end if
      end function objective

   end function garma_fit

   pure function garma_fitted(series, u, d, ar, ma) result(fitted)
      !! Return fitted values from GARMA conditional residuals.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      real(dp), intent(in) :: u(:) !! Input vector or random variate.
      real(dp), intent(in) :: d(:) !! Fractional-differencing parameter or differencing order.
      real(dp), intent(in) :: ar(:) !! Autoregressive coefficients.
      real(dp), intent(in) :: ma(:) !! Moving-average coefficients.
      real(dp), allocatable :: fitted(:), residual(:)

      residual = garma_css_residuals(series, u, d, ar, ma)
      fitted = series - residual
   end function garma_fitted

   pure function garma_forecast(series, u, d, ar, ma, horizon) result(out)
      !! Forecast a GARMA process from its truncated infinite-AR representation.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      real(dp), intent(in) :: u(:) !! Input vector or random variate.
      real(dp), intent(in) :: d(:) !! Fractional-differencing parameter or differencing order.
      real(dp), intent(in) :: ar(:) !! Autoregressive coefficients.
      real(dp), intent(in) :: ma(:) !! Moving-average coefficients.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      type(garma_forecast_t) :: out
      real(dp), allocatable :: operator(:), factor(:), denominator(:), history(:)
      integer :: n, j, h

      n = size(series)
      if (n < 1 .or. horizon < 1 .or. size(u) /= size(d)) then
         out%info = 1
         return
      end if
      allocate(operator(0:n + horizon - 1))
      operator = 0.0_dp
      operator(0) = 1.0_dp
      if (size(ar) > 0) operator(1:size(ar)) = -ar
      do j = 1, size(u)
         factor = garma_gegenbauer_coefficients(n + horizon, -d(j), u(j))
         operator = polynomial_product(operator, factor, n + horizon - 1)
      end do
      allocate(denominator(0:size(ma)))
      denominator = 0.0_dp
      denominator(0) = 1.0_dp
      if (size(ma) > 0) denominator(1:) = ma
      operator = polynomial_division(operator, denominator, n + horizon - 1)
      allocate(history(n + horizon), out%mean(horizon))
      history(:n) = series
      do h = 1, horizon
         history(n + h) = 0.0_dp
         do j = 1, min(n + h - 1, ubound(operator, 1))
            history(n + h) = history(n + h) - operator(j)*history(n + h - j)
         end do
         out%mean(h) = history(n + h)
      end do
   end function garma_forecast

   pure function garma_regression_fit(series, initial_u, initial_d, initial_ar, &
      initial_ma, regressors, include_mean, include_drift, difference_order, method, &
      estimate_frequencies, d_limits, max_iterations, tolerance) result(out)
      !! Fit deterministic regression first, then GARMA dynamics to differenced residuals.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      real(dp), intent(in) :: initial_u(:) !! Initial u.
      real(dp), intent(in) :: initial_d(:) !! Initial d.
      real(dp), intent(in) :: initial_ar(:) !! Initial autoregressive.
      real(dp), intent(in) :: initial_ma(:) !! Initial moving-average.
      real(dp), intent(in), optional :: regressors(:, :) !! Regression design matrix.
      real(dp), intent(in), optional :: d_limits(:) !! D limits.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      logical, intent(in), optional :: include_mean !! Whether to include a mean term.
      logical, intent(in), optional :: include_drift !! Whether to include the drift.
      logical, intent(in), optional :: estimate_frequencies !! Whether to estimate the frequencies.
      integer, intent(in), optional :: difference_order !! Difference order.
      integer, intent(in), optional :: method !! Algorithm or estimation method.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      type(garma_regression_fit_t) :: out
      real(dp), allocatable :: design(:, :), regression_fitted(:), beta(:), se(:)
      real(dp), allocatable :: regression_residual(:), dynamic_fitted(:), covariance(:, :)
      real(dp) :: rss
      integer :: order, columns, status, inversion_status
      logical :: mean_term, drift_term

      mean_term = .false.
      if (present(include_mean)) mean_term = include_mean
      drift_term = .false.
      if (present(include_drift)) drift_term = include_drift
      order = 0
      if (present(difference_order)) order = difference_order
      out%difference_order = order
      out%includes_mean = mean_term .and. order == 0
      out%includes_drift = drift_term
      if (present(regressors)) out%external_regressors = size(regressors, 2)
      if (size(series) < 16 .or. order < 0 .or. order >= size(series) .or. &
         (present(regressors) .and. size(regressors, 1) /= size(series))) then
         out%info = 1
         return
      end if
      design = regression_design(size(series), regressors, out%includes_mean, drift_term, 0)
      columns = size(design, 2)
      if (columns > 0) then
         call ols_fit(design, series, beta, se, regression_residual, rss, status)
         if (status /= 0) then
            out%info = 2
            return
         end if
         out%regression_coefficients = beta
         out%regression_standard_error = se
         call invert_matrix(matmul(transpose(design), design), covariance, inversion_status)
         if (inversion_status == 0 .and. size(series) > columns) then
            out%regression_covariance = covariance*rss/real(size(series) - columns, dp)
         else
            allocate(out%regression_covariance(columns, columns))
            out%regression_covariance = 0.0_dp
         end if
         regression_fitted = matmul(design, beta)
      else
         allocate(out%regression_coefficients(0), out%regression_standard_error(0))
         allocate(out%regression_covariance(0, 0), regression_fitted(size(series)))
         regression_fitted = 0.0_dp
         regression_residual = series
      end if
      out%regression_residuals = regression_residual
      out%differenced_series = repeated_difference(regression_residual, order)
      if (present(d_limits)) then
         out%dynamic = garma_fit(out%differenced_series, initial_u, initial_d, &
            initial_ar, initial_ma, method, estimate_frequencies, d_limits, &
            max_iterations, tolerance)
      else
         out%dynamic = garma_fit(out%differenced_series, initial_u, initial_d, &
            initial_ar, initial_ma, method, estimate_frequencies, &
            max_iterations=max_iterations, tolerance=tolerance)
      end if
      if (out%dynamic%info /= 0) then
         out%info = 10 + out%dynamic%info
         return
      end if
      dynamic_fitted = reintegrate_fitted(out%dynamic%fitted, regression_residual, order)
      out%fitted = regression_fitted + dynamic_fitted
      out%residuals = series - out%fitted
   end function garma_regression_fit

   pure function garma_regression_forecast(fit, future_regressors, horizon) result(out)
      !! Forecast a regression GARMA fit and restore integration and deterministic terms.
      type(garma_regression_fit_t), intent(in) :: fit !! Previously fitted model.
      real(dp), intent(in), optional :: future_regressors(:, :) !! Future regressors.
      integer, intent(in) :: horizon !! Number of periods to forecast.
      type(garma_forecast_t) :: out
      type(garma_forecast_t) :: dynamic_forecast
      real(dp), allocatable :: design(:, :), restored(:)

      if (fit%info /= 0 .or. horizon < 1 .or. &
         (fit%external_regressors > 0 .and. .not. present(future_regressors))) then
         out%info = 1
         return
      end if
      if (present(future_regressors)) then
         if (size(future_regressors, 1) /= horizon .or. &
            size(future_regressors, 2) /= fit%external_regressors) then
            out%info = 1
            return
         end if
      end if
      dynamic_forecast = garma_forecast(fit%differenced_series, fit%dynamic%u, &
         fit%dynamic%d, fit%dynamic%ar, fit%dynamic%ma, horizon)
      if (dynamic_forecast%info /= 0) then
         out%info = 2
         return
      end if
      restored = reintegrate_forecast(dynamic_forecast%mean, fit%regression_residuals, &
         fit%difference_order)
      design = regression_design(horizon, future_regressors, fit%includes_mean, &
         fit%includes_drift, size(fit%regression_residuals))
      out%mean = restored
      if (size(design, 2) > 0) then
         out%mean = out%mean + matmul(design, fit%regression_coefficients)
      end if
   end function garma_regression_forecast

   pure function garma_accuracy(actual, predicted, training, seasonal_lag) result(out)
      !! Compute forecast package point-accuracy measures.
      real(dp), intent(in) :: actual(:) !! Observed values used for evaluation.
      real(dp), intent(in) :: predicted(:) !! Predicted values.
      real(dp), intent(in) :: training(:) !! Training observations.
      integer, intent(in), optional :: seasonal_lag !! Seasonal lag.
      type(garma_accuracy_t) :: out
      real(dp), allocatable :: error(:), scaled_difference(:)
      real(dp) :: center_error, denominator
      integer :: lag

      lag = 1
      if (present(seasonal_lag)) lag = seasonal_lag
      if (size(actual) < 1 .or. size(predicted) /= size(actual) .or. &
         size(training) <= lag .or. lag < 1) then
         out%info = 1
         return
      end if
      error = actual - predicted
      out%observations = size(error)
      out%mean_error = sum(error)/real(size(error), dp)
      out%root_mean_squared_error = sqrt(sum(error**2)/real(size(error), dp))
      out%mean_absolute_error = sum(abs(error))/real(size(error), dp)
      if (all(abs(actual) > tiny(1.0_dp))) then
         out%mean_percentage_error = 100.0_dp*sum(error/actual)/real(size(error), dp)
         out%mean_absolute_percentage_error = 100.0_dp*sum(abs(error/actual))/ &
            real(size(error), dp)
      end if
      scaled_difference = abs(training(lag + 1:) - training(:size(training) - lag))
      denominator = sum(scaled_difference)/real(size(scaled_difference), dp)
      if (denominator > tiny(1.0_dp)) then
         out%mean_absolute_scaled_error = out%mean_absolute_error/denominator
      end if
      if (size(error) > 1) then
         center_error = out%mean_error
         denominator = sum((error - center_error)**2)
         if (denominator > tiny(1.0_dp)) out%residual_acf1 = &
            dot_product(error(:size(error) - 1) - center_error, &
            error(2:) - center_error)/denominator
      end if
   end function garma_accuracy

   pure function garma_goodness_of_fit(residuals) result(out)
      !! Apply garma's sequential Bartlett periodogram white-noise diagnostic.
      real(dp), intent(in) :: residuals(:) !! Model residuals.
      type(garma_gof_t) :: out
      type(garma_periodogram_t) :: pgram
      real(dp), allocatable :: sorted_full(:), sorted_prefix(:)
      real(dp) :: empirical, d_plus, d_minus, root_n, probability, term
      integer :: frequencies, prefix, i, j, minimum_index

      pgram = undetrended_periodogram(residuals)
      frequencies = size(pgram%spectrum)
      if (pgram%info /= 0 .or. frequencies < 2) then
         out%info = 1
         return
      end if
      sorted_full = sorted(pgram%spectrum)
      allocate(out%statistic(frequencies - 1), out%p_value(frequencies - 1))
      do prefix = 2, frequencies
         sorted_prefix = sorted(pgram%spectrum(:prefix))
         d_plus = 0.0_dp
         d_minus = 0.0_dp
         do i = 1, prefix
            empirical = real(count(sorted_full <= sorted_prefix(i)), dp)/ &
               real(frequencies, dp)
            d_plus = max(d_plus, real(i, dp)/real(prefix, dp) - empirical)
            d_minus = max(d_minus, empirical - real(i - 1, dp)/real(prefix, dp))
         end do
         out%statistic(prefix - 1) = max(d_plus, d_minus)
         root_n = sqrt(real(prefix, dp))
         probability = 0.0_dp
         do j = 1, 100
            term = 2.0_dp*(-1.0_dp)**(j - 1)*exp(-2.0_dp*real(j*j, dp)* &
               (root_n*out%statistic(prefix - 1))**2)
            probability = probability + term
            if (abs(term) < 1.0e-14_dp) exit
         end do
         out%p_value(prefix - 1) = max(0.0_dp, min(1.0_dp, probability))
      end do
      minimum_index = minloc(out%p_value, dim=1)
      out%minimum_p_value = out%p_value(minimum_index)
      out%minimum_frequency = 2.0_dp*acos(-1.0_dp)*real(minimum_index + 1, dp)/ &
         real(size(residuals), dp)
      out%minimum_period = real(size(residuals), dp)/real(minimum_index + 1, dp)
   end function garma_goodness_of_fit

   pure function undetrended_periodogram(series) result(out)
      !! Compute the raw periodogram used by the residual Bartlett diagnostic.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      type(garma_periodogram_t) :: out
      complex(dp), allocatable :: transformed(:)
      integer :: n, i

      n = size(series)
      if (n < 4 .or. .not. all(ieee_is_finite(series))) then
         out%info = 1
         return
      end if
      transformed = discrete_transform(cmplx(series, 0.0_dp, dp))
      allocate(out%frequency(n/2), out%spectrum(n/2))
      do i = 1, n/2
         out%frequency(i) = real(i, dp)/real(n, dp)
         out%spectrum(i) = abs(transformed(i + 1))**2/real(n, dp)
      end do
   end function undetrended_periodogram

   pure function gsp_estimate(pgram, pole_index, bandwidth) result(estimate)
      !! Minimize Arteche's Gaussian semiparametric criterion by golden search.
      type(garma_periodogram_t), intent(in) :: pgram !! Periodogram values and frequencies.
      integer, intent(in) :: pole_index !! Index of pole.
      integer, intent(in) :: bandwidth !! Smoothing or spectral bandwidth.
      real(dp) :: estimate, lower, upper, first, second, f_first, f_second
      real(dp), parameter :: ratio = 0.6180339887498948482_dp
      integer :: iteration

      lower = -10.0_dp
      upper = 10.0_dp
      first = upper - ratio*(upper - lower)
      second = lower + ratio*(upper - lower)
      f_first = gsp_criterion(first, pgram, pole_index, bandwidth)
      f_second = gsp_criterion(second, pgram, pole_index, bandwidth)
      do iteration = 1, 120
         if (f_first < f_second) then
            upper = second
            second = first
            f_second = f_first
            first = upper - ratio*(upper - lower)
            f_first = gsp_criterion(first, pgram, pole_index, bandwidth)
         else
            lower = first
            first = second
            f_first = f_second
            second = lower + ratio*(upper - lower)
            f_second = gsp_criterion(second, pgram, pole_index, bandwidth)
         end if
      end do
      estimate = 0.5_dp*(lower + upper)
   end function gsp_estimate

   pure real(dp) function gsp_criterion(d, pgram, pole_index, bandwidth) result(value)
      !! Evaluate the local Gaussian semiparametric criterion.
      real(dp), intent(in) :: d !! Fractional-differencing parameter or differencing order.
      type(garma_periodogram_t), intent(in) :: pgram !! Periodogram values and frequencies.
      integer, intent(in) :: pole_index !! Index of pole.
      integer, intent(in) :: bandwidth !! Smoothing or spectral bandwidth.
      real(dp) :: omega(bandwidth), spectrum(bandwidth)
      integer :: i, index, total

      total = size(pgram%spectrum)
      do i = 1, bandwidth
         omega(i) = 2.0_dp*acos(-1.0_dp)*pgram%frequency(i)
         index = pole_index + i - 1
         if (index > total) index = 2*total - index + 1
         index = max(1, min(total, index))
         spectrum(i) = pgram%spectrum(index)
      end do
      value = log(sum(omega**(2.0_dp*d)*spectrum)/real(bandwidth, dp)) - &
         2.0_dp*d*sum(log(omega))/real(bandwidth, dp)
   end function gsp_criterion

   pure real(dp) function lpr_estimate(pgram, pole_index, bandwidth) result(estimate)
      !! Compute the local symmetric log-periodogram regression estimate.
      type(garma_periodogram_t), intent(in) :: pgram !! Periodogram values and frequencies.
      integer, intent(in) :: pole_index !! Index of pole.
      integer, intent(in) :: bandwidth !! Smoothing or spectral bandwidth.
      real(dp), allocatable :: centered_log(:)
      real(dp) :: denominator, numerator
      integer :: i, plus_index, minus_index, total

      if (bandwidth < 2) then
         estimate = 0.0_dp
         return
      end if
      allocate(centered_log(bandwidth - 1))
      do i = 1, bandwidth - 1
         centered_log(i) = log(real(i, dp))
      end do
      centered_log = centered_log - sum(centered_log)/real(bandwidth - 1, dp)
      denominator = 4.0_dp*sum(centered_log**2)
      numerator = 0.0_dp
      total = size(pgram%spectrum)
      do i = 1, bandwidth - 1
         plus_index = reflected_spectrum_index(pole_index + i, total)
         minus_index = reflected_spectrum_index(pole_index - i + 2, total)
         numerator = numerator + centered_log(i)*(log(pgram%spectrum(plus_index)) + &
            log(pgram%spectrum(minus_index)))
      end do
      estimate = -0.5_dp*numerator/denominator
   end function lpr_estimate

   pure elemental integer function reflected_spectrum_index(index, total) result(reflected)
      !! Reflect a periodogram index across zero and the Nyquist frequency.
      integer, intent(in) :: index !! Element or observation index.
      integer, intent(in) :: total !! Total.

      reflected = index
      if (reflected < 1) reflected = 2 - reflected
      if (reflected > total) reflected = 2*total - reflected + 1
      reflected = max(1, min(total, reflected))
   end function reflected_spectrum_index

   pure function rational_filter(series, numerator, denominator) result(filtered)
      !! Apply a causal rational lag filter with unit-leading denominator.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      real(dp), intent(in) :: numerator(0:) !! Numerator polynomial coefficients.
      real(dp), intent(in) :: denominator(0:) !! Denominator polynomial coefficients.
      real(dp), allocatable :: filtered(:)
      integer :: t, lag

      allocate(filtered(size(series)))
      filtered = 0.0_dp
      do t = 1, size(series)
         do lag = 0, min(ubound(numerator, 1), t - 1)
            filtered(t) = filtered(t) + numerator(lag)*series(t - lag)
         end do
         do lag = 1, min(ubound(denominator, 1), t - 1)
            filtered(t) = filtered(t) - denominator(lag)*filtered(t - lag)
         end do
         filtered(t) = filtered(t)/denominator(0)
      end do
   end function rational_filter

   pure function discrete_transform(values) result(transformed)
      !! Use the radix-2 FFT when possible and an exact DFT otherwise.
      complex(dp), intent(in) :: values(:) !! Input values.
      complex(dp), allocatable :: transformed(:)
      real(dp) :: angle
      integer :: n, k, j

      n = size(values)
      if (n > 0 .and. iand(n, n - 1) == 0) then
         transformed = fft_transform(values)
         return
      end if
      allocate(transformed(n))
      transformed = cmplx(0.0_dp, 0.0_dp, dp)
      do k = 0, n - 1
         do j = 0, n - 1
            angle = -2.0_dp*acos(-1.0_dp)*real(j*k, dp)/real(n, dp)
            transformed(k + 1) = transformed(k + 1) + values(j + 1)* &
               cmplx(cos(angle), sin(angle), dp)
         end do
      end do
   end function discrete_transform

   pure function encode_fit(u, d, ar, ma, estimate_u, lower_d, upper_d, variance, &
      include_variance) result(parameters)
      !! Map bounded GARMA parameters into unconstrained coordinates.
      real(dp), intent(in) :: u(:) !! Input vector or random variate.
      real(dp), intent(in) :: d(:) !! Fractional-differencing parameter or differencing order.
      real(dp), intent(in) :: ar(:) !! Autoregressive coefficients.
      real(dp), intent(in) :: ma(:) !! Moving-average coefficients.
      real(dp), intent(in) :: lower_d !! Lower bound for d.
      real(dp), intent(in) :: upper_d !! Upper bound for d.
      real(dp), intent(in) :: variance !! Variance value or matrix.
      logical, intent(in) :: estimate_u !! Whether to estimate the u.
      logical, intent(in) :: include_variance !! Whether to include the variance.
      real(dp), allocatable :: parameters(:)
      real(dp) :: scaled
      integer :: factor, offset

      allocate(parameters(size(d) + size(ar) + size(ma) + &
         merge(size(u), 0, estimate_u) + merge(1, 0, include_variance)))
      offset = 0
      do factor = 1, size(d)
         if (estimate_u) then
            scaled = max(1.0e-8_dp, min(1.0_dp - 1.0e-8_dp, u(factor)))
            parameters(offset + 1) = log(scaled/(1.0_dp - scaled))
            offset = offset + 1
         end if
         scaled = max(1.0e-8_dp, min(1.0_dp - 1.0e-8_dp, &
            (d(factor) - lower_d)/(upper_d - lower_d)))
         parameters(offset + 1) = log(scaled/(1.0_dp - scaled))
         offset = offset + 1
      end do
      parameters(offset + 1:offset + size(ar)) = atanh(max(-0.999999_dp, &
         min(0.999999_dp, ar/0.999_dp)))
      offset = offset + size(ar)
      parameters(offset + 1:offset + size(ma)) = atanh(max(-0.999999_dp, &
         min(0.999999_dp, ma/0.999_dp)))
      offset = offset + size(ma)
      if (include_variance) parameters(offset + 1) = log(max(variance, tiny(1.0_dp)))
   end function encode_fit

   pure subroutine decode_fit(parameters, fixed_u, estimate_u, lower_d, upper_d, p, q, &
      include_variance, u, d, ar, ma, variance)
      !! Recover bounded GARMA parameters from optimizer coordinates.
      real(dp), intent(in) :: parameters(:) !! Model parameter values.
      real(dp), intent(in) :: fixed_u(:) !! Fixed u.
      real(dp), intent(in) :: lower_d !! Lower bound for d.
      real(dp), intent(in) :: upper_d !! Upper bound for d.
      logical, intent(in) :: estimate_u !! Whether to estimate the u.
      logical, intent(in) :: include_variance !! Whether to include the variance.
      integer, intent(in) :: p !! Autoregressive order or model dimension.
      integer, intent(in) :: q !! Model order, dimension, or parameter.
      real(dp), allocatable, intent(out) :: u(:) !! Input vector or random variate.
      real(dp), allocatable, intent(out) :: d(:) !! Fractional-differencing parameter or differencing order.
      real(dp), allocatable, intent(out) :: ar(:) !! Autoregressive coefficients.
      real(dp), allocatable, intent(out) :: ma(:) !! Moving-average coefficients.
      real(dp), intent(out) :: variance !! Variance value or matrix.
      integer :: factor, offset

      allocate(u(size(fixed_u)), d(size(fixed_u)), ar(p), ma(q))
      u = fixed_u
      offset = 0
      do factor = 1, size(u)
         if (estimate_u) then
            u(factor) = 0.5_dp*(1.0_dp + tanh(0.5_dp*parameters(offset + 1)))
            offset = offset + 1
         end if
         d(factor) = lower_d + (upper_d - lower_d)* &
            0.5_dp*(1.0_dp + tanh(0.5_dp*parameters(offset + 1)))
         offset = offset + 1
      end do
      if (estimate_u .and. size(u) > 1) call sort_pairs(u, d)
      ar = 0.999_dp*tanh(parameters(offset + 1:offset + p))
      offset = offset + p
      ma = 0.999_dp*tanh(parameters(offset + 1:offset + q))
      offset = offset + q
      variance = 1.0_dp
      if (include_variance) variance = exp(min(300.0_dp, parameters(offset + 1)))
   end subroutine decode_fit

   pure subroutine sort_pairs(u, d)
      !! Sort pole locations while retaining their paired exponents.
      real(dp), intent(inout) :: u(:) !! Input vector or random variate, updated in place.
      real(dp), intent(inout) :: d(:) !! Fractional-differencing parameter or differencing order, updated in place.
      real(dp) :: saved_u, saved_d
      integer :: i, j

      do i = 2, size(u)
         saved_u = u(i)
         saved_d = d(i)
         j = i - 1
         do while (j >= 1)
            if (u(j) <= saved_u) exit
            u(j + 1) = u(j)
            d(j + 1) = d(j)
            j = j - 1
         end do
         u(j + 1) = saved_u
         d(j + 1) = saved_d
      end do
   end subroutine sort_pairs

   pure subroutine spectral_covariance(fit, series, estimate_u, covariance, info)
      !! Estimate Whittle covariance from frequency-domain score products.
      type(garma_fit_t), intent(in) :: fit !! Previously fitted model.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      logical, intent(in) :: estimate_u !! Whether to estimate the u.
      real(dp), allocatable, intent(out) :: covariance(:, :) !! Covariance matrix.
      integer, intent(out) :: info !! Status code; zero indicates success.
      type(garma_periodogram_t) :: pgram
      real(dp), allocatable :: physical(:), plus(:), minus(:), score(:, :), omega(:, :)
      real(dp) :: step
      integer :: count, frequency, parameter

      count = size(fit%d) + size(fit%ar) + size(fit%ma) + &
         merge(size(fit%u), 0, estimate_u)
      allocate(covariance(count, count))
      covariance = 0.0_dp
      info = 0
      if (fit%method == garma_method_wll .or. count == 0) return
      physical = pack_physical(fit, estimate_u)
      pgram = garma_periodogram(series)
      allocate(score(size(pgram%frequency), count), plus(count), minus(count))
      do parameter = 1, count
         step = 1.0e-5_dp*max(1.0_dp, abs(physical(parameter)))
         plus = physical
         minus = physical
         plus(parameter) = plus(parameter) + step
         minus(parameter) = minus(parameter) - step
         score(:, parameter) = (log(max(spectral_from_physical(pgram%frequency, plus, &
            fit, estimate_u), tiny(1.0_dp))) - &
            log(max(spectral_from_physical(pgram%frequency, minus, fit, estimate_u), &
            tiny(1.0_dp))))/(2.0_dp*step)
      end do
      allocate(omega(count, count))
      omega = 0.0_dp
      do frequency = 1, size(pgram%frequency)
         omega = omega + spread(score(frequency, :), 2, count)* &
            spread(score(frequency, :), 1, count)
      end do
      omega = omega/real(size(pgram%frequency), dp)
      call invert_matrix(real(size(series), dp)*omega, covariance, info)
   end subroutine spectral_covariance

   pure subroutine css_covariance(fit, series, estimate_u, covariance, info)
      !! Estimate CSS covariance from the physical-parameter objective Hessian.
      type(garma_fit_t), intent(in) :: fit !! Previously fitted model.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      logical, intent(in) :: estimate_u !! Whether to estimate the u.
      real(dp), allocatable, intent(out) :: covariance(:, :) !! Covariance matrix.
      integer, intent(out) :: info !! Status code; zero indicates success.
      real(dp), allocatable :: parameters(:), hessian(:, :), inverse(:, :)

      parameters = pack_physical(fit, estimate_u)
      if (size(parameters) == 0) then
         allocate(covariance(0, 0))
         info = 0
         return
      end if
      hessian = finite_difference_hessian(objective, parameters, 1.0e-4_dp)
      call invert_matrix(hessian, inverse, info)
      if (info == 0) then
         covariance = 2.0_dp*fit%innovation_variance*inverse
      else
         allocate(covariance(size(parameters), size(parameters)))
         covariance = 0.0_dp
      end if

   contains

      pure function objective(parameters) result(value)
         !! Evaluate CSS at a packed physical parameter vector.
         real(dp), intent(in) :: parameters(:) !! Model parameter values.
         real(dp) :: value
         real(dp), allocatable :: u(:), d(:), ar(:), ma(:), residual(:)

         call unpack_physical(parameters, fit, estimate_u, u, d, ar, ma)
         residual = garma_css_residuals(series, u, d, ar, ma)
         value = sum(residual**2)
      end function objective

   end subroutine css_covariance

   pure subroutine set_wll_standard_errors(fit, series, estimate_u)
      !! Set the specialized WLL standard errors available for Gegenbauer exponents.
      type(garma_fit_t), intent(inout) :: fit !! Previously fitted model, updated in place.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      logical, intent(in) :: estimate_u !! Whether to estimate the u.
      type(garma_periodogram_t) :: pgram
      real(dp), allocatable :: score(:)
      integer :: factor, index

      pgram = garma_periodogram(series)
      do factor = 1, size(fit%d)
         score = log(max(4.0_dp*(cos(2.0_dp*acos(-1.0_dp)*pgram%frequency) - &
            fit%u(factor))**2, tiny(1.0_dp)))
         index = factor
         if (estimate_u) index = 2*factor
         fit%standard_error(index) = acos(-1.0_dp)**2/ &
            (6.0_dp*sqrt(sum(score**2)))
      end do
   end subroutine set_wll_standard_errors

   pure subroutine unpack_physical(parameters, fit, estimate_u, u, d, ar, ma)
      !! Unpack package-ordered physical parameters using a fit as the shape template.
      real(dp), intent(in) :: parameters(:) !! Model parameter values.
      type(garma_fit_t), intent(in) :: fit !! Previously fitted model.
      logical, intent(in) :: estimate_u !! Whether to estimate the u.
      real(dp), allocatable, intent(out) :: u(:) !! Input vector or random variate.
      real(dp), allocatable, intent(out) :: d(:) !! Fractional-differencing parameter or differencing order.
      real(dp), allocatable, intent(out) :: ar(:) !! Autoregressive coefficients.
      real(dp), allocatable, intent(out) :: ma(:) !! Moving-average coefficients.
      integer :: factor, offset

      allocate(u(size(fit%u)), d(size(fit%d)), ar(size(fit%ar)), ma(size(fit%ma)))
      u = fit%u
      offset = 0
      do factor = 1, size(d)
         if (estimate_u) then
            u(factor) = parameters(offset + 1)
            offset = offset + 1
         end if
         d(factor) = parameters(offset + 1)
         offset = offset + 1
      end do
      ar = parameters(offset + 1:offset + size(ar))
      offset = offset + size(ar)
      ma = parameters(offset + 1:)
   end subroutine unpack_physical

   pure function pack_physical(fit, estimate_u) result(parameters)
      !! Pack GARMA parameters in package coefficient order.
      type(garma_fit_t), intent(in) :: fit !! Previously fitted model.
      logical, intent(in) :: estimate_u !! Whether to estimate the u.
      real(dp), allocatable :: parameters(:)
      integer :: factor, offset

      allocate(parameters(size(fit%d) + size(fit%ar) + size(fit%ma) + &
         merge(size(fit%u), 0, estimate_u)))
      offset = 0
      do factor = 1, size(fit%d)
         if (estimate_u) then
            parameters(offset + 1) = fit%u(factor)
            offset = offset + 1
         end if
         parameters(offset + 1) = fit%d(factor)
         offset = offset + 1
      end do
      parameters(offset + 1:offset + size(fit%ar)) = fit%ar
      offset = offset + size(fit%ar)
      parameters(offset + 1:) = fit%ma
   end function pack_physical

   pure function spectral_from_physical(frequency, parameters, fit, estimate_u) result(values)
      !! Evaluate spectral inverse from a packed physical parameter vector.
      real(dp), intent(in) :: frequency(:) !! Number of observations per seasonal cycle.
      real(dp), intent(in) :: parameters(:) !! Model parameter values.
      type(garma_fit_t), intent(in) :: fit !! Previously fitted model.
      logical, intent(in) :: estimate_u !! Whether to estimate the u.
      real(dp), allocatable :: values(:), u(:), d(:), ar(:), ma(:)
      integer :: factor, offset

      call unpack_physical(parameters, fit, estimate_u, u, d, ar, ma)
      values = garma_spectral_inverse(frequency, u, d, ar, ma)
   end function spectral_from_physical

   pure function polynomial_product(first, second, max_degree) result(product)
      !! Multiply two zero-origin polynomials with truncation.
      real(dp), intent(in) :: first(0:) !! First operand.
      real(dp), intent(in) :: second(:) !! Second operand.
      integer, intent(in) :: max_degree !! Maximum degree.
      real(dp), allocatable :: product(:)
      integer :: i, j

      allocate(product(0:max_degree))
      product = 0.0_dp
      do i = 0, min(ubound(first, 1), max_degree)
         do j = 0, min(size(second) - 1, max_degree - i)
            product(i + j) = product(i + j) + first(i)*second(j + 1)
         end do
      end do
   end function polynomial_product

   pure function polynomial_division(numerator, denominator, max_degree) result(quotient)
      !! Divide zero-origin lag polynomials through a requested degree.
      real(dp), intent(in) :: numerator(0:) !! Numerator polynomial coefficients.
      real(dp), intent(in) :: denominator(0:) !! Denominator polynomial coefficients.
      integer, intent(in) :: max_degree !! Maximum degree.
      real(dp), allocatable :: quotient(:)
      integer :: degree, lag

      allocate(quotient(0:max_degree))
      quotient = 0.0_dp
      do degree = 0, max_degree
         if (degree <= ubound(numerator, 1)) quotient(degree) = numerator(degree)
         do lag = 1, min(degree, ubound(denominator, 1))
            quotient(degree) = quotient(degree) - denominator(lag)*quotient(degree - lag)
         end do
         quotient(degree) = quotient(degree)/denominator(0)
      end do
   end function polynomial_division

   pure function regression_design(observations, regressors, include_mean, &
      include_drift, time_offset) result(design)
      !! Build external, intercept, and drift columns in package regression order.
      integer, intent(in) :: observations !! Observed time-series values.
      integer, intent(in) :: time_offset !! Time offset.
      real(dp), intent(in), optional :: regressors(:, :) !! Regression design matrix.
      logical, intent(in) :: include_mean !! Whether to include a mean term.
      logical, intent(in) :: include_drift !! Whether to include the drift.
      real(dp), allocatable :: design(:, :)
      integer :: columns, offset, i

      columns = merge(1, 0, include_mean) + merge(1, 0, include_drift)
      if (present(regressors)) columns = columns + size(regressors, 2)
      allocate(design(observations, columns))
      offset = 0
      if (include_mean) then
         design(:, 1) = 1.0_dp
         offset = 1
      end if
      if (include_drift) then
         do i = 1, observations
            design(i, offset + 1) = real(time_offset + i, dp)
         end do
         offset = offset + 1
      end if
      if (present(regressors)) design(:, offset + 1:) = regressors
   end function regression_design

   pure function repeated_difference(series, order) result(differenced)
      !! Apply ordinary differencing repeatedly.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      integer, intent(in) :: order !! Model or polynomial order.
      real(dp), allocatable :: differenced(:), next(:)
      integer :: iteration

      differenced = series
      do iteration = 1, order
         next = differenced(2:) - differenced(:size(differenced) - 1)
         call move_alloc(next, differenced)
      end do
   end function repeated_difference

   pure function reintegrate_fitted(fitted_difference, original, order) result(fitted)
      !! Restore fitted differences using leading observed initial conditions.
      real(dp), intent(in) :: fitted_difference(:) !! Fitted difference.
      real(dp), intent(in) :: original(:) !! Original.
      integer, intent(in) :: order !! Model or polynomial order.
      real(dp), allocatable :: fitted(:)
      integer :: t, lag

      if (order == 0) then
         fitted = fitted_difference
         return
      end if
      allocate(fitted(size(fitted_difference) + order))
      fitted(:order) = original(:order)
      do t = order + 1, size(fitted)
         fitted(t) = fitted_difference(t - order)
         do lag = 1, order
            fitted(t) = fitted(t) + (-1.0_dp)**(lag + 1)* &
               binomial_coefficient(order, lag)*fitted(t - lag)
         end do
      end do
   end function reintegrate_fitted

   pure function reintegrate_forecast(forecast_difference, history, order) result(forecast)
      !! Restore future differences from the latest observed integration states.
      real(dp), intent(in) :: forecast_difference(:) !! Forecast difference.
      real(dp), intent(in) :: history(:) !! History.
      integer, intent(in) :: order !! Model or polynomial order.
      real(dp), allocatable :: forecast(:), states(:), differences(:)
      integer :: level, h, n

      if (order == 0) then
         forecast = forecast_difference
         return
      end if
      n = size(history)
      allocate(states(order), differences(n))
      differences = history
      do level = 1, order
         states(level) = differences(size(differences))
         if (level < order) differences = differences(2:) - &
            differences(:size(differences) - 1)
      end do
      allocate(forecast(size(forecast_difference)))
      do h = 1, size(forecast)
         states(order) = states(order) + forecast_difference(h)
         do level = order - 1, 1, -1
            states(level) = states(level) + states(level + 1)
         end do
         forecast(h) = states(1)
      end do
   end function reintegrate_forecast

   pure elemental real(dp) function binomial_coefficient(n, k) result(value)
      !! Return a small integer binomial coefficient in double precision.
      integer, intent(in) :: n !! Number of observations or elements.
      integer, intent(in) :: k !! K.
      integer :: i

      value = 1.0_dp
      do i = 1, min(k, n - k)
         value = value*real(n - i + 1, dp)/real(i, dp)
      end do
   end function binomial_coefficient

end module garma_mod
