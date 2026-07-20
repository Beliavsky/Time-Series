! SPDX-License-Identifier: GPL-2.0-or-later
! SPDX-FileComment: Algorithms translated from the R fracdiff package.
! Distinct fractional-differencing algorithms translated from fracdiff.
module fracdiff_mod
   use kind_mod, only: dp
   use arfima_mod, only: arfima_model_t, arfima_model, arfima_acvf_t, &
      arfima_simulation_t, arfima_fdwn_acvf, arfima_durbin_levinson_simulate
   use arima2_mod, only: arima2_roots_t, arma_polynomial_roots, &
      durbin_levinson_coefficients
   use fourier_mod, only: fft_transform
   use linalg_mod, only: invert_matrix
   use optimization_mod, only: optimization_result_t, &
      bfgs_minimize_fd, finite_difference_hessian
   use random_mod, only: set_random_seed, random_standard_normal
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   private

   type, public :: fracdiff_semiparametric_t
      ! Semiparametric memory estimate and asymptotic and regression errors.
      real(dp) :: d = 0.0_dp
      real(dp) :: asymptotic_standard_error = 0.0_dp
      real(dp) :: regression_standard_error = 0.0_dp
      integer :: bandwidth = 0
      integer :: retained_frequencies = 0
      integer :: info = 0
   end type fracdiff_semiparametric_t

   type, public :: fracdiff_filter_t
      ! Haslett-Raftery filtered series, profiled mean, and variance determinant.
      real(dp), allocatable :: filtered(:), prediction_variance(:)
      real(dp) :: mean = 0.0_dp
      real(dp) :: log_variance = 0.0_dp
      integer :: truncation = 0
      integer :: info = 0
   end type fracdiff_filter_t

   type, public :: fracdiff_fit_t
      ! Haslett-Raftery ARFIMA fit and observed parameter inference.
      type(arfima_model_t) :: model
      real(dp), allocatable :: covariance(:, :), standard_error(:)
      real(dp), allocatable :: unconstrained_parameters(:)
      real(dp), allocatable :: residuals(:), fitted(:)
      real(dp) :: log_likelihood = -huge(1.0_dp)
      real(dp) :: innovation_variance = 0.0_dp
      real(dp) :: aic = huge(1.0_dp)
      real(dp) :: bic = huge(1.0_dp)
      integer :: truncation = 100
      integer :: iterations = 0
      integer :: info = 0
      logical :: converged = .false.
   end type fracdiff_fit_t

   type, public :: fracdiff_simulation_t
      ! Simulated fracdiff series and the burn-in length discarded.
      real(dp), allocatable :: series(:)
      integer :: burn_in = 0
      integer :: info = 0
   end type fracdiff_simulation_t

   public :: fracdiff_difference
   public :: fracdiff_gph, fracdiff_sperio
   public :: fracdiff_hr_filter, fracdiff_fit
   public :: fracdiff_simulate_from_innovations, fracdiff_simulate

contains

   pure function fracdiff_difference(series, d) result(differenced)
      !! Fractionally difference a demeaned series by FFT convolution.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      real(dp), intent(in) :: d !! Fractional-differencing parameter or differencing order.
      real(dp), allocatable :: differenced(:)
      complex(dp), allocatable :: series_transform(:), weight_transform(:)
      complex(dp), allocatable :: padded_series(:), weights(:), product(:)
      integer :: n, nfft, k

      n = size(series)
      if (n < 2 .or. .not. ieee_is_finite(d) .or. &
         .not. all(ieee_is_finite(series))) then
         allocate(differenced(0))
         return
      end if
      nfft = 1
      do while (nfft < 2*n - 1)
         nfft = 2*nfft
      end do
      allocate(padded_series(nfft), weights(nfft))
      padded_series = cmplx(0.0_dp, 0.0_dp, dp)
      weights = cmplx(0.0_dp, 0.0_dp, dp)
      padded_series(:n) = cmplx(series - sum(series)/real(n, dp), 0.0_dp, dp)
      weights(1) = cmplx(1.0_dp, 0.0_dp, dp)
      do k = 1, n - 1
         weights(k + 1) = weights(k)*real(k - 1, dp)/real(k, dp) - &
            weights(k)*d/real(k, dp)
      end do
      series_transform = fft_transform(padded_series)
      weight_transform = fft_transform(weights)
      product = series_transform*weight_transform
      product = fft_transform(product, inverse=.true.)
      allocate(differenced(n))
      differenced = real(product(:n), dp)
   end function fracdiff_difference

   pure function fracdiff_gph(series, bandwidth_exponent) result(out)
      !! Estimate fractional memory by Geweke-Porter-Hudak regression.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      real(dp), intent(in), optional :: bandwidth_exponent !! Bandwidth exponent.
      type(fracdiff_semiparametric_t) :: out
      real(dp) :: exponent

      exponent = 0.5_dp
      if (present(bandwidth_exponent)) exponent = bandwidth_exponent
      out = semiparametric_fit(series, exponent, .false., 0.9_dp)
   end function fracdiff_gph

   pure function fracdiff_sperio(series, bandwidth_exponent, taper_exponent) result(out)
      !! Estimate fractional memory from the tapered Sperio periodogram.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      real(dp), intent(in), optional :: bandwidth_exponent !! Bandwidth exponent.
      real(dp), intent(in), optional :: taper_exponent !! Taper exponent.
      type(fracdiff_semiparametric_t) :: out
      real(dp) :: exponent, beta

      exponent = 0.5_dp
      if (present(bandwidth_exponent)) exponent = bandwidth_exponent
      beta = 0.9_dp
      if (present(taper_exponent)) beta = taper_exponent
      out = semiparametric_fit(series, exponent, .true., beta)
   end function fracdiff_sperio

   pure function fracdiff_hr_filter(series, d, truncation) result(out)
      !! Apply the Haslett-Raftery truncated fractional likelihood filter.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      real(dp), intent(in) :: d !! Fractional-differencing parameter or differencing order.
      integer, intent(in), optional :: truncation !! Truncation.
      type(fracdiff_filter_t) :: out
      real(dp), allocatable :: prediction(:), mean_loading(:), reflection(:)
      real(dp), allocatable :: pi_weight(:)
      real(dp) :: g0, value, denominator, numerator, sum_pi, accumulated
      integer :: n, selected_truncation, active, k, j, previous, offset

      n = size(series)
      selected_truncation = 100
      if (present(truncation)) selected_truncation = truncation
      out%truncation = selected_truncation
      if (n < 2 .or. selected_truncation < 1 .or. d < -0.5_dp .or. d >= 0.5_dp .or. &
         .not. all(ieee_is_finite(series))) then
         out%info = 1
         return
      end if
      active = min(selected_truncation, n)
      allocate(prediction(n), mean_loading(n), out%prediction_variance(n))
      allocate(reflection(max(1, active - 1)), pi_weight(max(1, active)))
      prediction = 0.0_dp
      mean_loading = 1.0_dp
      out%prediction_variance = 1.0_dp
      reflection = 0.0_dp
      pi_weight = 0.0_dp
      g0 = exp(log_gamma(1.0_dp - 2.0_dp*d) - &
         2.0_dp*log_gamma(1.0_dp - d))
      out%prediction_variance(1) = g0
      if (active >= 2) then
         value = d/(1.0_dp - d)
         reflection(1) = value
         prediction(2) = value*series(1)
         mean_loading(2) = 1.0_dp - value
         out%prediction_variance(2) = g0*(1.0_dp - value**2)
      end if
      do k = 3, active
         previous = k - 1
         do j = 1, previous - 1
            value = real(previous - j, dp)
            reflection(j) = reflection(j)*real(previous, dp)*(value - d)/ &
               ((real(previous, dp) - d)*value)
         end do
         reflection(previous) = d/(real(previous, dp) - d)
         out%prediction_variance(k) = out%prediction_variance(k - 1)* &
            (1.0_dp - reflection(previous)**2)
         prediction(k) = 0.0_dp
         mean_loading(k) = 1.0_dp
         do j = 1, previous
            prediction(k) = prediction(k) + reflection(j)*series(k - j)
            mean_loading(k) = mean_loading(k) - reflection(j)
         end do
      end do
      if (active < n) then
         pi_weight(1) = d
         sum_pi = d
         do j = 2, active
            pi_weight(j) = pi_weight(j - 1)*(real(j - 1, dp) - d)/real(j, dp)
            sum_pi = sum_pi + pi_weight(j)
         end do
         accumulated = 0.0_dp
         do k = active + 1, n
            offset = k - active
            prediction(k) = 0.0_dp
            do j = 1, active
               prediction(k) = prediction(k) + pi_weight(j)*series(k - j)
            end do
            if (offset > 1 .and. abs(d) > sqrt(epsilon(1.0_dp))) then
               value = real(active, dp)*pi_weight(active)* &
                  (1.0_dp - (real(active, dp)/real(k, dp))**d)/d
               prediction(k) = prediction(k) + &
                  value*accumulated/real(offset - 1, dp)
               mean_loading(k) = 1.0_dp - sum_pi - value
            else
               mean_loading(k) = 1.0_dp - sum_pi
            end if
            accumulated = accumulated + series(offset)
            out%prediction_variance(k) = 1.0_dp
         end do
      end if
      numerator = 0.0_dp
      denominator = 0.0_dp
      do k = 1, n
         value = mean_loading(k)
         if (k <= active) then
            numerator = numerator + (series(k) - prediction(k))*value/ &
               out%prediction_variance(k)
            denominator = denominator + value**2/out%prediction_variance(k)
         else
            numerator = numerator + (series(k) - prediction(k))*value
            denominator = denominator + value**2
         end if
      end do
      if (denominator <= tiny(1.0_dp) .or. &
         any(out%prediction_variance(:active) <= tiny(1.0_dp))) then
         out%info = 2
         return
      end if
      out%mean = numerator/denominator
      allocate(out%filtered(n))
      out%filtered = series - prediction - out%mean*mean_loading
      out%filtered(:active) = out%filtered(:active)/ &
         sqrt(out%prediction_variance(:active))
      out%log_variance = sum(log(out%prediction_variance(:active)))
   end function fracdiff_hr_filter

   pure function fracdiff_fit(series, initial_ar, initial_ma, initial_d, d_range, &
      truncation, max_iterations, tolerance) result(out)
      !! Fit Haslett-Raftery ARFIMA parameters by stable-coordinate BFGS.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      real(dp), intent(in) :: initial_ar(:) !! Initial autoregressive.
      real(dp), intent(in) :: initial_ma(:) !! Initial moving-average.
      real(dp), intent(in) :: initial_d !! Initial d.
      real(dp), intent(in), optional :: d_range(2) !! D range.
      real(dp), intent(in), optional :: tolerance !! Numerical convergence tolerance.
      integer, intent(in), optional :: truncation !! Truncation.
      integer, intent(in), optional :: max_iterations !! Maximum number of algorithm iterations.
      type(fracdiff_fit_t) :: out
      type(optimization_result_t) :: optimization
      type(arfima_model_t) :: initial_model
      real(dp), allocatable :: initial(:), hessian(:, :), inverse(:, :)
      real(dp), allocatable :: jacobian(:, :), operating_covariance(:, :)
      real(dp), allocatable :: ar(:), ma(:), residuals(:), differenced(:)
      real(dp) :: selected_range(2), variance, likelihood, step
      integer :: selected_truncation, selected_iterations, status, i, j, npar
      real(dp) :: selected_tolerance

      selected_range = [0.0_dp, 0.5_dp]
      if (present(d_range)) selected_range = d_range
      selected_truncation = 100
      if (present(truncation)) selected_truncation = truncation
      selected_iterations = 200
      if (present(max_iterations)) selected_iterations = max_iterations
      selected_tolerance = 1.0e-6_dp
      if (present(tolerance)) selected_tolerance = tolerance
      initial_model = arfima_model(initial_ar, initial_ma, initial_d)
      if (size(series) < max(3, max(size(initial_ar), size(initial_ma)) + 2) .or. &
         initial_model%info /= 0 .or. selected_range(1) < -0.5_dp .or. &
         selected_range(2) > 0.5_dp .or. selected_range(1) >= selected_range(2) .or. &
         selected_truncation < 1) then
         out%info = 1
         return
      end if
      initial = encode_parameters(initial_ar, initial_ma, initial_d, selected_range)
      optimization = bfgs_minimize_fd(objective, initial, &
         max_iterations=selected_iterations, gradient_tolerance=selected_tolerance)
      out%iterations = optimization%iterations
      out%converged = optimization%converged
      out%unconstrained_parameters = optimization%parameters
      call decode_parameters(optimization%parameters, size(initial_ar), size(initial_ma), &
         selected_range, ar, ma, out%model%long_memory_parameter)
      out%model = arfima_model(ar, ma, out%model%long_memory_parameter)
      call hr_likelihood(series, ar, ma, out%model%long_memory_parameter, &
         selected_truncation, likelihood, variance, residuals, status)
      if (status /= 0) then
         out%info = status
         return
      end if
      out%log_likelihood = likelihood
      out%innovation_variance = variance
      out%model%innovation_variance = variance
      out%truncation = selected_truncation
      differenced = fracdiff_difference(series, out%model%long_memory_parameter)
      out%residuals = arma_residuals(differenced, ar, ma)
      out%fitted = series - out%residuals
      npar = size(initial)
      hessian = finite_difference_hessian(objective, optimization%parameters)
      call invert_matrix(hessian, inverse, status)
      allocate(out%covariance(npar, npar), out%standard_error(npar))
      out%covariance = 0.0_dp
      out%standard_error = 0.0_dp
      if (status == 0) then
         allocate(jacobian(npar, npar))
         do i = 1, npar
            step = 1.0e-5_dp*max(1.0_dp, abs(optimization%parameters(i)))
            jacobian(:, i) = operating_parameters( &
               optimization%parameters + unit_step(npar, i, step), &
               size(initial_ar), size(initial_ma), selected_range) - &
               operating_parameters(optimization%parameters - &
               unit_step(npar, i, step), size(initial_ar), size(initial_ma), &
               selected_range)
            jacobian(:, i) = jacobian(:, i)/(2.0_dp*step)
         end do
         operating_covariance = matmul(matmul(jacobian, inverse), transpose(jacobian))
         if (all([(operating_covariance(i, i) > 0.0_dp, i=1, npar)])) then
            out%covariance = operating_covariance
            do i = 1, npar
               out%standard_error(i) = sqrt(out%covariance(i, i))
            end do
         end if
      end if
      out%aic = -2.0_dp*out%log_likelihood + 2.0_dp*real(npar + 1, dp)
      out%bic = -2.0_dp*out%log_likelihood + &
         log(real(size(series), dp))*real(npar + 1, dp)
      if (.not. optimization%converged) out%info = optimization%info

   contains

      pure function objective(parameters) result(value)
         !! Return the profiled Haslett-Raftery negative log likelihood.
         real(dp), intent(in) :: parameters(:) !! Model parameter values.
         real(dp) :: value, candidate_d, candidate_variance, candidate_likelihood
         real(dp), allocatable :: candidate_ar(:), candidate_ma(:)
         real(dp), allocatable :: candidate_residuals(:)
         integer :: candidate_status

         call decode_parameters(parameters, size(initial_ar), size(initial_ma), &
            selected_range, candidate_ar, candidate_ma, candidate_d)
         call hr_likelihood(series, candidate_ar, candidate_ma, candidate_d, &
            selected_truncation, candidate_likelihood, candidate_variance, &
            candidate_residuals, candidate_status)
         if (candidate_status /= 0 .or. .not. ieee_is_finite(candidate_likelihood)) then
            value = huge(1.0_dp)/100.0_dp
         else
            value = -candidate_likelihood
         end if
      end function objective

   end function fracdiff_fit

   pure function fracdiff_simulate_from_innovations(observations, ar, ma, d, &
      innovations, mean, burn_in) result(out)
      !! Simulate fracdiff's FDWN process and opposite-sign MA recursion.
      integer, intent(in) :: observations !! Observed time-series values.
      real(dp), intent(in) :: ar(:) !! Autoregressive coefficients.
      real(dp), intent(in) :: ma(:) !! Moving-average coefficients.
      real(dp), intent(in) :: d !! Fractional-differencing parameter or differencing order.
      real(dp), intent(in) :: innovations(:) !! Model innovations.
      real(dp), intent(in), optional :: mean !! Mean value or vector.
      integer, intent(in), optional :: burn_in !! Burn in.
      type(fracdiff_simulation_t) :: out
      type(arfima_acvf_t) :: covariance
      type(arfima_simulation_t) :: fractional
      real(dp), allocatable :: filtered(:)
      real(dp) :: selected_mean
      integer :: selected_burn_in, total, required, k, lag

      selected_mean = 0.0_dp
      if (present(mean)) selected_mean = mean
      selected_burn_in = default_burn_in(ar, ma)
      if (present(burn_in)) selected_burn_in = burn_in
      out%burn_in = selected_burn_in
      total = observations + selected_burn_in
      required = total + size(ma)
      if (observations < 1 .or. selected_burn_in < 0 .or. &
         size(innovations) < required .or. d < -0.5_dp .or. d > 0.5_dp .or. &
         .not. all(ieee_is_finite(innovations(:required)))) then
         out%info = 1
         return
      end if
      covariance = arfima_fdwn_acvf(d, required - 1)
      if (covariance%info /= 0) then
         out%info = covariance%info
         return
      end if
      fractional = arfima_durbin_levinson_simulate(innovations(:required), &
         covariance%covariance)
      if (fractional%info /= 0) then
         out%info = fractional%info
         return
      end if
      allocate(filtered(total))
      filtered = 0.0_dp
      do k = 1, total
         do lag = 1, min(size(ar), k - 1)
            filtered(k) = filtered(k) + ar(lag)*filtered(k - lag)
         end do
         do lag = 1, size(ma)
            filtered(k) = filtered(k) - &
               ma(lag)*fractional%series(k + size(ma) - lag)
         end do
         filtered(k) = filtered(k) + fractional%series(k + size(ma))
      end do
      out%series = filtered(selected_burn_in + 1:) + selected_mean
   end function fracdiff_simulate_from_innovations

   function fracdiff_simulate(observations, ar, ma, d, mean, burn_in, seed) result(out)
      !! Simulate a Gaussian fracdiff process using the shared random stream.
      integer, intent(in) :: observations !! Observed time-series values.
      real(dp), intent(in) :: ar(:) !! Autoregressive coefficients.
      real(dp), intent(in) :: ma(:) !! Moving-average coefficients.
      real(dp), intent(in) :: d !! Fractional-differencing parameter or differencing order.
      real(dp), intent(in), optional :: mean !! Mean value or vector.
      integer, intent(in), optional :: burn_in !! Burn in.
      integer, intent(in), optional :: seed !! Random-number seed.
      type(fracdiff_simulation_t) :: out
      real(dp), allocatable :: innovations(:)
      real(dp) :: selected_mean
      integer :: selected_burn_in, i, required

      selected_mean = 0.0_dp
      if (present(mean)) selected_mean = mean
      selected_burn_in = default_burn_in(ar, ma)
      if (present(burn_in)) selected_burn_in = burn_in
      required = observations + selected_burn_in + size(ma)
      if (required < 1) then
         out%info = 1
         return
      end if
      if (present(seed)) call set_random_seed(seed)
      allocate(innovations(required))
      do i = 1, required
         innovations(i) = random_standard_normal()
      end do
      out = fracdiff_simulate_from_innovations(observations, ar, ma, d, &
         innovations, selected_mean, selected_burn_in)
   end function fracdiff_simulate

   pure function semiparametric_fit(series, exponent, tapered, beta) result(out)
      !! Fit the common log-periodogram regression used by GPH and Sperio.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      real(dp), intent(in) :: exponent !! Exponent.
      real(dp), intent(in) :: beta !! Regression or model coefficients.
      logical, intent(in) :: tapered !! Flag controlling tapered.
      type(fracdiff_semiparametric_t) :: out
      real(dp), allocatable :: centered(:), covariance(:), periodogram(:)
      real(dp), allocatable :: predictor(:), response(:), taper(:), residual(:)
      real(dp) :: frequency, predictor_mean, response_mean, slope, xss, rss, ratio
      integer :: n, bandwidth, maximum_taper, half_taper, lag, j, retained

      n = size(series)
      if (n < 4 .or. exponent <= 0.0_dp .or. exponent >= 1.0_dp .or. &
         beta <= 0.0_dp .or. beta > 1.0_dp .or. &
         .not. all(ieee_is_finite(series))) then
         out%info = 1
         return
      end if
      bandwidth = int(real(n, dp)**exponent)
      if (bandwidth < 2 .or. bandwidth >= n) then
         out%info = 1
         return
      end if
      out%bandwidth = bandwidth
      centered = series - sum(series)/real(n, dp)
      allocate(covariance(n - 1), taper(n - 1), periodogram(bandwidth))
      do lag = 1, n - 1
         covariance(lag) = dot_product(centered(:n - lag), centered(lag + 1:))/ &
            real(n, dp)
      end do
      taper = 1.0_dp
      maximum_taper = int(real(n, dp)**beta)
      half_taper = maximum_taper/2
      if (tapered) then
         do lag = 1, n - 1
            ratio = real(lag, dp)/real(maximum_taper, dp)
            if (lag <= half_taper) then
               taper(lag) = 1.0_dp - 6.0_dp*ratio**2*(1.0_dp - ratio)
            else if (lag <= maximum_taper) then
               taper(lag) = 2.0_dp*(1.0_dp - ratio)**3
            else
               taper(lag) = 0.0_dp
            end if
         end do
      end if
      do j = 1, bandwidth
         frequency = 2.0_dp*acos(-1.0_dp)*real(j, dp)/real(n, dp)
         periodogram(j) = sum(centered**2)/real(n, dp)
         do lag = 1, n - 1
            periodogram(j) = periodogram(j) + &
               2.0_dp*covariance(lag)*taper(lag)*cos(frequency*real(lag, dp))
         end do
      end do
      retained = count(periodogram > 0.0_dp)
      out%retained_frequencies = retained
      if (retained < 2) then
         out%info = 2
         return
      end if
      allocate(predictor(retained), response(retained), residual(retained))
      retained = 0
      do j = 1, bandwidth
         if (periodogram(j) <= 0.0_dp) cycle
         retained = retained + 1
         frequency = 2.0_dp*acos(-1.0_dp)*real(j, dp)/real(n, dp)
         predictor(retained) = 2.0_dp*log(2.0_dp*sin(0.5_dp*frequency))
         response(retained) = log(periodogram(j)/(2.0_dp*acos(-1.0_dp)))
      end do
      predictor_mean = sum(predictor)/real(retained, dp)
      response_mean = sum(response)/real(retained, dp)
      xss = sum((predictor - predictor_mean)**2)
      if (xss <= tiny(1.0_dp)) then
         out%info = 2
         return
      end if
      slope = dot_product(predictor - predictor_mean, &
         response - response_mean)/xss
      residual = response - response_mean - slope*(predictor - predictor_mean)
      rss = sum(residual**2)
      out%d = -slope
      if (tapered) then
         out%asymptotic_standard_error = &
            sqrt((0.539285_dp*real(maximum_taper, dp)/real(n, dp))/xss)
      else
         out%asymptotic_standard_error = sqrt(acos(-1.0_dp)**2/(6.0_dp*xss))
      end if
      out%regression_standard_error = &
         sqrt(rss/(real(bandwidth - 1, dp)*xss))
   end function semiparametric_fit

   pure subroutine hr_likelihood(series, ar, ma, d, truncation, likelihood, &
      variance, residuals, info)
      !! Evaluate the profiled Haslett-Raftery likelihood for one parameter set.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      real(dp), intent(in) :: ar(:) !! Autoregressive coefficients.
      real(dp), intent(in) :: ma(:) !! Moving-average coefficients.
      real(dp), intent(in) :: d !! Fractional-differencing parameter or differencing order.
      integer, intent(in) :: truncation !! Truncation.
      real(dp), intent(out) :: likelihood !! Likelihood.
      real(dp), intent(out) :: variance !! Variance value or matrix.
      real(dp), allocatable, intent(out) :: residuals(:) !! Model residuals.
      integer, intent(out) :: info !! Status code; zero indicates success.
      type(fracdiff_filter_t) :: filtered
      integer :: n, start, t, lag

      filtered = fracdiff_hr_filter(series, d, truncation)
      if (filtered%info /= 0) then
         info = filtered%info
         likelihood = -huge(1.0_dp)
         variance = huge(1.0_dp)
         allocate(residuals(0))
         return
      end if
      n = size(series)
      start = max(size(ar), size(ma)) + 1
      if (start > n) then
         info = 2
         likelihood = -huge(1.0_dp)
         variance = huge(1.0_dp)
         return
      end if
      if (size(ar) + size(ma) == 0) then
         residuals = filtered%filtered
         variance = sum(residuals**2)/real(n, dp)
      else
         if (n - start < 1) then
            info = 2
            likelihood = -huge(1.0_dp)
            variance = huge(1.0_dp)
            return
         end if
         allocate(residuals(n))
         residuals = 0.0_dp
         do t = start, n
            residuals(t) = filtered%filtered(t)
            do lag = 1, size(ar)
               residuals(t) = residuals(t) - ar(lag)*filtered%filtered(t - lag)
            end do
            do lag = 1, min(size(ma), t - start)
               residuals(t) = residuals(t) + ma(lag)*residuals(t - lag)
            end do
         end do
         variance = sum(residuals(start:)**2)/real(n - start, dp)
      end if
      if (variance <= tiny(1.0_dp) .or. .not. ieee_is_finite(variance)) then
         info = 2
         likelihood = -huge(1.0_dp)
         return
      end if
      likelihood = -0.5_dp*(real(n, dp)*(log(variance) + 2.8378_dp) + &
         filtered%log_variance)
      info = 0
   end subroutine hr_likelihood

   pure function arma_residuals(series, ar, ma) result(residuals)
      !! Apply fracdiff's AR and opposite-sign MA recursion.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      real(dp), intent(in) :: ar(:) !! Autoregressive coefficients.
      real(dp), intent(in) :: ma(:) !! Moving-average coefficients.
      real(dp), allocatable :: residuals(:)
      integer :: t, lag

      allocate(residuals(size(series)))
      residuals = 0.0_dp
      do t = 1, size(series)
         residuals(t) = series(t)
         do lag = 1, min(size(ar), t - 1)
            residuals(t) = residuals(t) - ar(lag)*series(t - lag)
         end do
         do lag = 1, min(size(ma), t - 1)
            residuals(t) = residuals(t) + ma(lag)*residuals(t - lag)
         end do
      end do
   end function arma_residuals

   pure function encode_parameters(ar, ma, d, d_range) result(parameters)
      !! Map initial operating parameters into stable optimizer coordinates.
      real(dp), intent(in) :: ar(:) !! Autoregressive coefficients.
      real(dp), intent(in) :: ma(:) !! Moving-average coefficients.
      real(dp), intent(in) :: d !! Fractional-differencing parameter or differencing order.
      real(dp), intent(in) :: d_range(2) !! D range.
      real(dp), allocatable :: parameters(:), partial(:)
      real(dp) :: midpoint, half_width, scaled
      integer :: offset

      allocate(parameters(1 + size(ar) + size(ma)))
      midpoint = 0.5_dp*(d_range(1) + d_range(2))
      half_width = 0.5_dp*(d_range(2) - d_range(1))
      scaled = max(-1.0_dp + 1.0e-10_dp, &
         min(1.0_dp - 1.0e-10_dp, (d - midpoint)/half_width))
      parameters(1) = atanh(scaled)
      offset = 1
      if (size(ar) > 0) then
         partial = ar_to_partial(ar)
         parameters(offset + 1:offset + size(ar)) = guarded_atanh(partial)
         offset = offset + size(ar)
      end if
      if (size(ma) > 0) then
         partial = ar_to_partial(ma)
         parameters(offset + 1:offset + size(ma)) = guarded_atanh(partial)
      end if
   end function encode_parameters

   pure subroutine decode_parameters(parameters, ar_order, ma_order, d_range, &
      ar, ma, d)
      !! Decode stable optimizer coordinates into fracdiff operating parameters.
      real(dp), intent(in) :: parameters(:) !! Model parameter values.
      real(dp), intent(in) :: d_range(2) !! D range.
      integer, intent(in) :: ar_order !! Autoregressive order.
      integer, intent(in) :: ma_order !! Moving-average order.
      real(dp), allocatable, intent(out) :: ar(:) !! Autoregressive coefficients.
      real(dp), allocatable, intent(out) :: ma(:) !! Moving-average coefficients.
      real(dp), intent(out) :: d !! Fractional-differencing parameter or differencing order.
      real(dp) :: midpoint, half_width
      integer :: offset

      midpoint = 0.5_dp*(d_range(1) + d_range(2))
      half_width = 0.5_dp*(d_range(2) - d_range(1))
      d = midpoint + half_width*tanh(parameters(1))
      allocate(ar(ar_order), ma(ma_order))
      offset = 1
      if (ar_order > 0) then
         ar = durbin_levinson_coefficients(tanh(parameters(offset + 1: &
            offset + ar_order)))
         offset = offset + ar_order
      end if
      if (ma_order > 0) then
         ma = durbin_levinson_coefficients(tanh(parameters(offset + 1: &
            offset + ma_order)))
      end if
   end subroutine decode_parameters

   pure function operating_parameters(parameters, ar_order, ma_order, d_range) &
      result(operating)
      !! Return flattened d, AR, and MA operating parameters.
      real(dp), intent(in) :: parameters(:) !! Model parameter values.
      real(dp), intent(in) :: d_range(2) !! D range.
      integer, intent(in) :: ar_order !! Autoregressive order.
      integer, intent(in) :: ma_order !! Moving-average order.
      real(dp), allocatable :: operating(:), ar(:), ma(:)
      real(dp) :: d

      call decode_parameters(parameters, ar_order, ma_order, d_range, ar, ma, d)
      operating = [d, ar, ma]
   end function operating_parameters

   pure function unit_step(n, index, step) result(vector)
      !! Return a coordinate finite-difference step vector.
      integer, intent(in) :: n !! Number of observations or elements.
      integer, intent(in) :: index !! Element or observation index.
      real(dp), intent(in) :: step !! Step.
      real(dp) :: vector(n)

      vector = 0.0_dp
      vector(index) = step
   end function unit_step

   pure function guarded_atanh(values) result(transformed)
      !! Apply inverse hyperbolic tangent inside the open unit interval.
      real(dp), intent(in) :: values(:) !! Input values.
      real(dp) :: transformed(size(values))

      transformed = atanh(max(-1.0_dp + 1.0e-10_dp, &
         min(1.0_dp - 1.0e-10_dp, values)))
   end function guarded_atanh

   pure function ar_to_partial(coefficients) result(partial)
      !! Convert stable AR coefficients to Durbin-Levinson partial coefficients.
      real(dp), intent(in) :: coefficients(:) !! Model coefficients.
      real(dp), allocatable :: partial(:)
      real(dp), allocatable :: current(:), previous(:)
      real(dp) :: denominator
      integer :: order, stage, j

      order = size(coefficients)
      allocate(partial(order), current(order), previous(order))
      current = coefficients
      do stage = order, 1, -1
         partial(stage) = current(stage)
         if (stage == 1) exit
         denominator = 1.0_dp - partial(stage)**2
         if (denominator <= tiny(1.0_dp)) then
            partial = 0.0_dp
            return
         end if
         previous = current
         do j = 1, stage - 1
            current(j) = (previous(j) + partial(stage)* &
               previous(stage - j))/denominator
         end do
      end do
   end function ar_to_partial

   pure integer function default_burn_in(ar, ma) result(burn_in)
      !! Reproduce fracdiff.sim's AR-root burn-in rule.
      real(dp), intent(in) :: ar(:) !! Autoregressive coefficients.
      real(dp), intent(in) :: ma(:) !! Moving-average coefficients.
      type(arima2_roots_t) :: roots
      real(dp) :: minimum_root

      burn_in = size(ar) + size(ma)
      if (size(ar) == 0) return
      roots = arma_polynomial_roots(ar)
      if (roots%info /= 0 .or. size(roots%roots) == 0) then
         minimum_root = 1.01_dp
      else
         minimum_root = minval(abs(roots%roots))
         if (minimum_root <= 1.0_dp) minimum_root = 1.01_dp
      end if
      burn_in = burn_in + ceiling(6.0_dp/log(minimum_root))
   end function default_burn_in

end module fracdiff_mod
