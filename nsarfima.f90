! SPDX-License-Identifier: GPL-3.0-or-later
! SPDX-FileComment: Algorithms translated from the R nsarfima package.
! Nonstationary ARFIMA algorithms translated from nsarfima.
module nsarfima_mod
   use kind_mod, only: dp
   use arfima_mod, only: arfima_model_t, arfima_model, arfima_information_t, &
      arfima_long_memory_none, arfima_fractional_weights, arfima_fisher_information
   use time_series_fourier_mod, only: fft_transform
   use time_series_linalg_mod, only: invert_matrix
   use time_series_optimization_mod, only: optimization_result_t, bfgs_minimize_fd
   use time_series_random_mod, only: set_random_seed, random_standard_normal
   use itsmr_mod, only: regularized_gamma_q
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   private

   type, public :: nsarfima_filter_t
      ! Filtered innovations, residual autocorrelations, and status.
      real(dp), allocatable :: residuals(:), autocorrelation(:)
      integer :: info = 0
   end type nsarfima_filter_t

   type, public :: nsarfima_fit_t
      ! ARFIMA parameter estimates, inference, diagnostics, and optimizer status.
      real(dp), allocatable :: ar(:), ma(:), covariance(:, :), standard_error(:)
      real(dp), allocatable :: residuals(:), residual_autocorrelation(:)
      real(dp) :: d = 0.0_dp
      real(dp) :: mean = 0.0_dp
      real(dp) :: innovation_variance = 0.0_dp
      real(dp) :: objective = huge(1.0_dp)
      real(dp) :: portmanteau_statistic = 0.0_dp
      real(dp) :: portmanteau_p_value = 0.0_dp
      integer :: lag_max = 0
      integer :: iterations = 0
      integer :: info = 0
      logical :: converged = .false.
   end type nsarfima_fit_t

   type, public :: nsarfima_simulation_t
      ! Simulated series and the effective burn-in length.
      real(dp), allocatable :: series(:)
      integer :: burn_in = 0
      integer :: info = 0
   end type nsarfima_simulation_t

   public :: nsarfima_convolve, nsarfima_residuals, nsarfima_residual_acf
   public :: nsarfima_mde, nsarfima_pml
   public :: nsarfima_simulate_from_innovations, nsarfima_simulate

contains

   pure function nsarfima_convolve(first, second) result(convolution)
      ! Return the leading causal convolution using zero-padded FFTs.
      real(dp), intent(in) :: first(:), second(:)
      real(dp), allocatable :: convolution(:)
      complex(dp), allocatable :: a(:), b(:), transformed(:)
      integer :: n, nfft

      n = size(first)
      if (n < 1 .or. size(second) /= n .or. .not. all(ieee_is_finite(first)) .or. &
         .not. all(ieee_is_finite(second))) then
         allocate(convolution(0))
         return
      end if
      nfft = 1
      do while (nfft < 2*n - 1)
         nfft = 2*nfft
      end do
      allocate(a(nfft), b(nfft))
      a = cmplx(0.0_dp, 0.0_dp, dp)
      b = cmplx(0.0_dp, 0.0_dp, dp)
      a(:n) = cmplx(first, 0.0_dp, dp)
      b(:n) = cmplx(second, 0.0_dp, dp)
      transformed = fft_transform(a)*fft_transform(b)
      transformed = fft_transform(transformed, inverse=.true.)
      allocate(convolution(n))
      convolution = real(transformed(:n), dp)
   end function nsarfima_convolve

   pure function nsarfima_residuals(series, d, ar, ma, subtract_mean) result(residuals)
      ! Filter a possibly nonstationary ARFIMA series into innovations.
      real(dp), intent(in) :: series(:), d, ar(:), ma(:)
      logical, intent(in), optional :: subtract_mean
      real(dp), allocatable :: residuals(:)
      real(dp), allocatable :: work(:), weights(:)
      logical :: centered
      integer :: n, integer_order

      n = size(series)
      if (n < 1 .or. d < -1.0_dp .or. .not. ieee_is_finite(d) .or. &
         .not. all(ieee_is_finite(series)) .or. .not. all(ieee_is_finite(ar)) .or. &
         .not. all(ieee_is_finite(ma))) then
         allocate(residuals(0))
         return
      end if
      centered = .true.
      if (present(subtract_mean)) centered = subtract_mean
      allocate(work(n))
      work = series
      if (centered) work = work - sum(work)/real(n, dp)
      if (d /= 0.0_dp) then
         integer_order = floor(d + 0.5_dp)
         if (integer_order > 0) then
            work = leading_integer_difference(work, integer_order)
            weights = arfima_fractional_weights(d - real(integer_order, dp), n - 1)
         else
            weights = arfima_fractional_weights(d, n - 1)
         end if
         work = nsarfima_convolve(weights, work)
      end if
      if (size(ar) + size(ma) > 0) then
         weights = inverse_arma_weights(ar, ma, n - 1)
         work = nsarfima_convolve(weights, work)
      end if
      residuals = work
   end function nsarfima_residuals

   pure function nsarfima_residual_acf(series, d, ar, ma, lag_max, &
      subtract_mean) result(out)
      ! Compute Mayoral residual autocorrelations through the selected lag.
      real(dp), intent(in) :: series(:), d, ar(:), ma(:)
      integer, intent(in), optional :: lag_max
      logical, intent(in), optional :: subtract_mean
      type(nsarfima_filter_t) :: out
      real(dp) :: denominator
      integer :: selected_lag, n, lag
      logical :: centered

      n = size(series)
      selected_lag = floor(sqrt(real(max(0, n), dp)))
      if (present(lag_max)) selected_lag = lag_max
      centered = .true.
      if (present(subtract_mean)) centered = subtract_mean
      if (n < 2 .or. selected_lag < 0 .or. selected_lag >= n) then
         out%info = 1
         return
      end if
      out%residuals = nsarfima_residuals(series, d, ar, ma, centered)
      if (size(out%residuals) /= n) then
         out%info = 1
         return
      end if
      out%residuals = out%residuals - sum(out%residuals)/real(n, dp)
      denominator = sum(out%residuals**2)
      if (denominator <= tiny(1.0_dp)) then
         out%info = 2
         return
      end if
      allocate(out%autocorrelation(0:selected_lag))
      out%autocorrelation(0) = 1.0_dp
      do lag = 1, selected_lag
         out%autocorrelation(lag) = dot_product(out%residuals(:n - lag), &
            out%residuals(lag + 1:))/denominator
      end do
   end function nsarfima_residual_acf

   pure function nsarfima_mde(series, initial_ar, initial_ma, initial_d, d_range, &
      lag_max, estimate_mean, max_iterations, tolerance) result(out)
      ! Fit Mayoral's minimum-distance estimator to an ARFIMA model.
      real(dp), intent(in) :: series(:), initial_ar(:), initial_ma(:), initial_d
      real(dp), intent(in) :: d_range(:)
      integer, intent(in), optional :: lag_max, max_iterations
      logical, intent(in), optional :: estimate_mean
      real(dp), intent(in), optional :: tolerance
      type(nsarfima_fit_t) :: out
      type(optimization_result_t) :: optimized
      type(nsarfima_filter_t) :: diagnostic
      real(dp), allocatable :: initial(:), ar(:), ma(:), j_matrix(:, :), xi(:, :)
      real(dp), allocatable :: dynamic_covariance(:, :), dynamic_error(:), weights(:)
      real(dp) :: d, lower, upper, selected_tolerance
      integer :: p, q, count, selected_lag, selected_iterations, status, column, lag
      logical :: include_d, include_mean

      p = size(initial_ar)
      q = size(initial_ma)
      include_d = size(d_range) == 2
      include_mean = .true.
      if (present(estimate_mean)) include_mean = estimate_mean
      selected_lag = floor(sqrt(real(size(series), dp)))
      if (present(lag_max)) selected_lag = lag_max
      selected_iterations = 200
      if (present(max_iterations)) selected_iterations = max_iterations
      selected_tolerance = 1.0e-6_dp
      if (present(tolerance)) selected_tolerance = tolerance
      if (.not. valid_fit_input(series, initial_ar, initial_ma, initial_d, d_range, &
         selected_lag) .or. (.not. include_d .and. p + q == 0)) then
         out%info = 1
         return
      end if
      lower = 0.0_dp
      upper = 0.0_dp
      if (include_d) then
         lower = minval(d_range)
         upper = maxval(d_range)
      end if
      initial = encode_parameters(initial_d, initial_ar, initial_ma, include_d, lower, upper)
      optimized = bfgs_minimize_fd(objective, initial, selected_iterations, selected_tolerance)
      call decode_parameters(optimized%parameters, include_d, initial_d, lower, upper, &
         p, q, d, ar, ma)
      out%d = d
      out%ar = ar
      out%ma = ma
      out%objective = optimized%objective
      out%iterations = optimized%iterations
      out%converged = optimized%converged
      out%info = optimized%info
      if (optimized%info == 4 .and. ieee_is_finite(optimized%objective)) out%info = 0
      out%mean = 0.0_dp
      if (include_mean) out%mean = sum(series)/real(size(series), dp)
      out%residuals = nsarfima_residuals(series, d, ar, ma, include_mean)
      diagnostic = nsarfima_residual_acf(series, d, ar, ma, &
         floor(real(size(series), dp)**0.25_dp), include_mean)
      if (diagnostic%info == 0) then
         out%residual_autocorrelation = diagnostic%autocorrelation
         call set_portmanteau(out, size(series), p, q)
      end if
      count = p + q + merge(1, 0, include_d)
      allocate(j_matrix(selected_lag, count))
      j_matrix = 0.0_dp
      column = 0
      if (include_d) then
         column = 1
         do lag = 1, selected_lag
            j_matrix(lag, column) = -1.0_dp/real(lag, dp)
         end do
      end if
      do column = 1, p
         weights = arma_impulse_weights(ar, [real(dp) ::], selected_lag - column)
         j_matrix(column:, merge(1, 0, include_d) + column) = weights
      end do
      do column = 1, q
         weights = inverse_arma_weights([real(dp) ::], ma, selected_lag - column)
         j_matrix(column:, merge(1, 0, include_d) + p + column) = weights
      end do
      xi = matmul(transpose(j_matrix), j_matrix)
      call invert_matrix(xi, dynamic_covariance, status)
      if (status == 0) then
         dynamic_covariance = dynamic_covariance/real(size(series), dp)
      else
         dynamic_covariance = 0.0_dp
         out%info = max(out%info, 2)
      end if
      out%covariance = dynamic_covariance
      allocate(dynamic_error(count))
      do column = 1, count
         dynamic_error(column) = sqrt(max(0.0_dp, dynamic_covariance(column, column)))
      end do
      if (include_mean) then
         allocate(out%standard_error(count + 1))
         out%standard_error(1) = sample_standard_deviation(series)/sqrt(real(size(series), dp))
         out%standard_error(2:) = dynamic_error
      else
         out%standard_error = dynamic_error
      end if

   contains

      pure function objective(parameters) result(value)
         ! Return the sum of squared residual autocorrelations.
         real(dp), intent(in) :: parameters(:)
         real(dp) :: value, candidate_d
         real(dp), allocatable :: candidate_ar(:), candidate_ma(:)
         type(nsarfima_filter_t) :: filtered

         call decode_parameters(parameters, include_d, initial_d, lower, upper, p, q, &
            candidate_d, candidate_ar, candidate_ma)
         filtered = nsarfima_residual_acf(series, candidate_d, candidate_ar, &
            candidate_ma, selected_lag, include_mean)
         if (filtered%info /= 0) then
            value = huge(1.0_dp)/100.0_dp
         else
            value = sum(filtered%autocorrelation**2)
         end if
      end function objective

   end function nsarfima_mde

   pure function nsarfima_pml(series, initial_ar, initial_ma, initial_d, d_range, &
      estimate_mean, max_iterations, tolerance, information_resolution) result(out)
      ! Fit Beran's residual pseudo-likelihood ARFIMA estimator.
      real(dp), intent(in) :: series(:), initial_ar(:), initial_ma(:), initial_d
      real(dp), intent(in) :: d_range(:)
      logical, intent(in), optional :: estimate_mean
      integer, intent(in), optional :: max_iterations, information_resolution
      real(dp), intent(in), optional :: tolerance
      type(nsarfima_fit_t) :: out
      type(optimization_result_t) :: optimized
      type(nsarfima_filter_t) :: diagnostic
      real(dp), allocatable :: initial(:), ar(:), ma(:)
      real(dp) :: d, lower, upper, selected_tolerance
      integer :: p, q, selected_iterations, resolution, denominator
      logical :: include_d, include_mean

      p = size(initial_ar)
      q = size(initial_ma)
      include_d = size(d_range) == 2
      include_mean = .true.
      if (present(estimate_mean)) include_mean = estimate_mean
      selected_iterations = 200
      if (present(max_iterations)) selected_iterations = max_iterations
      selected_tolerance = 1.0e-6_dp
      if (present(tolerance)) selected_tolerance = tolerance
      resolution = 4096
      if (present(information_resolution)) resolution = information_resolution
      if (.not. valid_fit_input(series, initial_ar, initial_ma, initial_d, d_range, 1) .or. &
         (.not. include_d .and. p + q == 0) .or. resolution < 16) then
         out%info = 1
         return
      end if
      lower = 0.0_dp
      upper = 0.0_dp
      if (include_d) then
         lower = minval(d_range)
         upper = maxval(d_range)
      end if
      initial = encode_parameters(initial_d, initial_ar, initial_ma, include_d, lower, upper)
      optimized = bfgs_minimize_fd(objective, initial, selected_iterations, selected_tolerance)
      call decode_parameters(optimized%parameters, include_d, initial_d, lower, upper, &
         p, q, d, ar, ma)
      out%d = d
      out%ar = ar
      out%ma = ma
      out%objective = optimized%objective
      out%iterations = optimized%iterations
      out%converged = optimized%converged
      out%info = optimized%info
      if (optimized%info == 4 .and. ieee_is_finite(optimized%objective)) out%info = 0
      out%mean = 0.0_dp
      if (include_mean) out%mean = sum(series)/real(size(series), dp)
      out%residuals = nsarfima_residuals(series, d, ar, ma, include_mean)
      denominator = size(series) - 1 - floor(d + 0.5_dp)
      if (denominator > 0) out%innovation_variance = out%objective/real(denominator, dp)
      diagnostic = nsarfima_residual_acf(series, d, ar, ma, &
         floor(real(size(series), dp)**0.25_dp), include_mean)
      if (diagnostic%info == 0) then
         out%residual_autocorrelation = diagnostic%autocorrelation
         call set_portmanteau(out, size(series), p, q)
      end if
      call set_pml_covariance(out, include_d, include_mean, resolution, size(series), &
         sample_standard_deviation(series)/sqrt(real(size(series), dp)))

   contains

      pure function objective(parameters) result(value)
         ! Return the pseudo-likelihood residual sum of squares.
         real(dp), intent(in) :: parameters(:)
         real(dp) :: value, candidate_d
         real(dp), allocatable :: candidate_ar(:), candidate_ma(:), residual(:)

         call decode_parameters(parameters, include_d, initial_d, lower, upper, p, q, &
            candidate_d, candidate_ar, candidate_ma)
         residual = nsarfima_residuals(series, candidate_d, candidate_ar, &
            candidate_ma, include_mean)
         if (size(residual) /= size(series)) then
            value = huge(1.0_dp)/100.0_dp
         else
            value = sum(residual(2:)**2)
         end if
      end function objective

   end function nsarfima_pml

   pure function nsarfima_simulate_from_innovations(observations, d, ar, ma, &
      innovations, mean, burn_in, stationary_integration) result(out)
      ! Simulate nsarfima's truncated causal filter from supplied innovations.
      integer, intent(in) :: observations
      real(dp), intent(in) :: d, ar(:), ma(:), innovations(:)
      real(dp), intent(in), optional :: mean
      integer, intent(in), optional :: burn_in
      logical, intent(in), optional :: stationary_integration
      type(nsarfima_simulation_t) :: out
      real(dp), allocatable :: work(:), weights(:)
      real(dp) :: location, fractional_d
      integer :: selected_burn, integer_order, iteration, n
      logical :: stationary

      selected_burn = merge(100, 10, d >= 0.0_dp)
      if (present(burn_in)) selected_burn = burn_in
      location = 0.0_dp
      if (present(mean)) location = mean
      stationary = .false.
      if (present(stationary_integration)) stationary = stationary_integration
      if (observations < 1 .or. selected_burn < 0 .or. size(innovations) < observations .or. &
         d < -1.0_dp .or. .not. all(ieee_is_finite(innovations))) then
         out%info = 1
         return
      end if
      if (size(innovations) < observations + selected_burn) selected_burn = 0
      out%burn_in = selected_burn
      n = size(innovations)
      work = innovations
      integer_order = floor(d + 0.5_dp)
      fractional_d = d - real(integer_order, dp)
      if (stationary .and. d >= 0.5_dp) then
         weights = arfima_fractional_weights(-fractional_d, n - 1)
         work = nsarfima_convolve(weights, work)
         weights = arma_impulse_weights(ar, ma, n - 1)
         work = nsarfima_convolve(weights, work) + location
         do iteration = 1, integer_order
            work = cumulative_sum(work)
         end do
      else
         if (integer_order > 0) then
            do iteration = 1, integer_order
               work = cumulative_sum(work)
            end do
         else if (integer_order < 0) then
            do iteration = 1, -integer_order
               work = ordinary_difference(work)
            end do
         end if
         n = size(work)
         weights = arfima_fractional_weights(-fractional_d, n - 1)
         work = nsarfima_convolve(weights, work)
         weights = arma_impulse_weights(ar, ma, n - 1)
         work = nsarfima_convolve(weights, work) + location
      end if
      if (size(work) < observations) then
         out%info = 2
         return
      end if
      out%series = work(size(work) - observations + 1:)
   end function nsarfima_simulate_from_innovations

   function nsarfima_simulate(observations, d, ar, ma, mean, innovation_variance, &
      burn_in, stationary_integration, seed) result(out)
      ! Simulate an nsarfima model with Gaussian innovations.
      integer, intent(in) :: observations
      real(dp), intent(in) :: d, ar(:), ma(:)
      real(dp), intent(in), optional :: mean, innovation_variance
      integer, intent(in), optional :: burn_in, seed
      logical, intent(in), optional :: stationary_integration
      type(nsarfima_simulation_t) :: out
      real(dp), allocatable :: innovations(:)
      real(dp) :: location, variance
      integer :: selected_burn, i
      logical :: stationary

      selected_burn = merge(100, 10, d >= 0.0_dp)
      if (present(burn_in)) selected_burn = burn_in
      location = 0.0_dp
      if (present(mean)) location = mean
      variance = 1.0_dp
      if (present(innovation_variance)) variance = innovation_variance
      stationary = .false.
      if (present(stationary_integration)) stationary = stationary_integration
      if (observations < 1 .or. selected_burn < 0 .or. variance < 0.0_dp) then
         out%info = 1
         return
      end if
      if (present(seed)) call set_random_seed(seed)
      allocate(innovations(observations + selected_burn))
      do i = 1, size(innovations)
         innovations(i) = sqrt(variance)*random_standard_normal()
      end do
      out = nsarfima_simulate_from_innovations(observations, d, ar, ma, innovations, &
         location, selected_burn, stationary)
   end function nsarfima_simulate

   pure function leading_integer_difference(series, order) result(differenced)
      ! Difference a zero-prefixed series without shortening the result.
      real(dp), intent(in) :: series(:)
      integer, intent(in) :: order
      real(dp), allocatable :: differenced(:), work(:), next(:)
      integer :: iteration, n

      n = size(series)
      allocate(work(n + order))
      work = 0.0_dp
      work(order + 1:) = series
      do iteration = 1, order
         allocate(next(size(work) - 1))
         next = work(2:) - work(:size(work) - 1)
         call move_alloc(next, work)
      end do
      differenced = work
   end function leading_integer_difference

   pure function inverse_arma_weights(ar, ma, max_lag) result(weights)
      ! Expand the inverse ARFIMA short-memory filter.
      real(dp), intent(in) :: ar(:), ma(:)
      integer, intent(in) :: max_lag
      real(dp), allocatable :: weights(:)
      integer :: lag, j

      if (max_lag < 0) then
         allocate(weights(0))
         return
      end if
      allocate(weights(0:max_lag))
      weights = 0.0_dp
      weights(0) = 1.0_dp
      do lag = 1, max_lag
         if (lag <= size(ar)) weights(lag) = -ar(lag)
         do j = 1, min(lag, size(ma))
            weights(lag) = weights(lag) - ma(j)*weights(lag - j)
         end do
      end do
   end function inverse_arma_weights

   pure function arma_impulse_weights(ar, ma, max_lag) result(weights)
      ! Expand the R-sign ARMA impulse-response filter.
      real(dp), intent(in) :: ar(:), ma(:)
      integer, intent(in) :: max_lag
      real(dp), allocatable :: weights(:)
      integer :: lag, j

      if (max_lag < 0) then
         allocate(weights(0))
         return
      end if
      allocate(weights(0:max_lag))
      weights = 0.0_dp
      weights(0) = 1.0_dp
      do lag = 1, max_lag
         if (lag <= size(ma)) weights(lag) = ma(lag)
         do j = 1, min(lag, size(ar))
            weights(lag) = weights(lag) + ar(j)*weights(lag - j)
         end do
      end do
   end function arma_impulse_weights

   pure function encode_parameters(d, ar, ma, include_d, lower, upper) result(parameters)
      ! Map bounded package parameters to unconstrained optimizer coordinates.
      real(dp), intent(in) :: d, ar(:), ma(:), lower, upper
      logical, intent(in) :: include_d
      real(dp), allocatable :: parameters(:)
      real(dp) :: scaled
      integer :: offset

      allocate(parameters(size(ar) + size(ma) + merge(1, 0, include_d)))
      offset = 0
      if (include_d) then
         scaled = max(1.0e-8_dp, min(1.0_dp - 1.0e-8_dp, (d - lower)/(upper - lower)))
         parameters(1) = log(scaled/(1.0_dp - scaled))
         offset = 1
      end if
      parameters(offset + 1:offset + size(ar)) = atanh(max(-0.999999_dp, &
         min(0.999999_dp, ar/0.999_dp)))
      parameters(offset + size(ar) + 1:) = atanh(max(-0.999999_dp, &
         min(0.999999_dp, ma/0.999_dp)))
   end function encode_parameters

   pure subroutine decode_parameters(parameters, include_d, fixed_d, lower, upper, &
      p, q, d, ar, ma)
      ! Recover bounded package parameters from optimizer coordinates.
      real(dp), intent(in) :: parameters(:), fixed_d, lower, upper
      logical, intent(in) :: include_d
      integer, intent(in) :: p, q
      real(dp), intent(out) :: d
      real(dp), allocatable, intent(out) :: ar(:), ma(:)
      real(dp) :: logistic
      integer :: offset

      allocate(ar(p), ma(q))
      offset = 0
      d = fixed_d
      if (include_d) then
         logistic = 0.5_dp*(1.0_dp + tanh(0.5_dp*parameters(1)))
         d = lower + (upper - lower)*logistic
         offset = 1
      end if
      ar = 0.999_dp*tanh(parameters(offset + 1:offset + p))
      ma = 0.999_dp*tanh(parameters(offset + p + 1:offset + p + q))
   end subroutine decode_parameters

   pure logical function valid_fit_input(series, ar, ma, d, d_range, lag_max) result(valid)
      ! Check common MDE and pseudo-likelihood fit arguments.
      real(dp), intent(in) :: series(:), ar(:), ma(:), d, d_range(:)
      integer, intent(in) :: lag_max

      valid = size(series) >= 3 .and. lag_max >= 1 .and. lag_max < size(series) .and. &
         (size(d_range) == 1 .or. size(d_range) == 2) .and. &
         all(d_range >= -1.0_dp) .and. all(ieee_is_finite(series)) .and. &
         all(ieee_is_finite(ar)) .and. all(ieee_is_finite(ma)) .and. ieee_is_finite(d)
      if (size(d_range) == 1) valid = valid .and. d_range(1) == 0.0_dp
      if (size(d_range) == 2) valid = valid .and. &
         maxval(d_range) > minval(d_range) .and. d >= minval(d_range) .and. &
         d <= maxval(d_range)
   end function valid_fit_input

   pure subroutine set_portmanteau(out, observations, p, q)
      ! Compute the package's residual portmanteau statistic and probability.
      type(nsarfima_fit_t), intent(inout) :: out
      integer, intent(in) :: observations, p, q
      integer :: lag, degrees

      out%lag_max = ubound(out%residual_autocorrelation, 1)
      do lag = 1, out%lag_max
         out%portmanteau_statistic = out%portmanteau_statistic + &
            out%residual_autocorrelation(lag)**2/real(observations - lag, dp)
      end do
      out%portmanteau_statistic = real(observations*(observations + 2), dp)* &
         out%portmanteau_statistic
      degrees = out%lag_max - p - q
      if (degrees > 0) out%portmanteau_p_value = 1.0_dp - &
         regularized_gamma_q(0.5_dp*real(degrees, dp), &
         0.5_dp*out%portmanteau_statistic)
   end subroutine set_portmanteau

   pure subroutine set_pml_covariance(out, include_d, include_mean, resolution, &
      observations, mean_standard_error)
      ! Reorder shared spectral information into nsarfima's parameter convention.
      type(nsarfima_fit_t), intent(inout) :: out
      logical, intent(in) :: include_d, include_mean
      integer, intent(in) :: resolution, observations
      real(dp), intent(in) :: mean_standard_error
      type(arfima_model_t) :: model
      type(arfima_information_t) :: information
      real(dp), allocatable :: transform(:, :), dynamic_covariance(:, :), errors(:)
      integer :: p, q, count, source_count, i, offset

      p = size(out%ar)
      q = size(out%ma)
      count = 1 + p + q + merge(1, 0, include_d)
      if (include_d) then
         model = arfima_model(out%ar, -out%ma, out%d - &
            real(floor(out%d + 0.5_dp), dp))
      else
         model = arfima_model(out%ar, -out%ma, 0.0_dp, &
            long_memory_type=arfima_long_memory_none)
      end if
      information = arfima_fisher_information(model, exact=.true., &
         resolution=resolution, observations=observations)
      source_count = p + q + merge(1, 0, include_d)
      allocate(out%covariance(count, count))
      out%covariance = 0.0_dp
      out%covariance(1, 1) = 2.0_dp*out%innovation_variance**2/real(observations, dp)
      if (information%info == 0 .and. size(information%covariance, 1) == source_count) then
         allocate(transform(source_count, source_count))
         transform = 0.0_dp
         offset = merge(1, 0, include_d)
         if (include_d) transform(1, p + q + 1) = 1.0_dp
         do i = 1, p
            transform(offset + i, i) = 1.0_dp
         end do
         do i = 1, q
            transform(offset + p + i, p + i) = -1.0_dp
         end do
         dynamic_covariance = matmul(matmul(transform, information%covariance), &
            transpose(transform))
         out%covariance(2:, 2:) = dynamic_covariance
      else
         out%info = max(out%info, 3)
      end if
      allocate(errors(count))
      do i = 1, count
         errors(i) = sqrt(max(0.0_dp, out%covariance(i, i)))
      end do
      if (include_mean) then
         allocate(out%standard_error(count + 1))
         out%standard_error(1) = mean_standard_error
         out%standard_error(2:) = errors
      else
         out%standard_error = errors
      end if
   end subroutine set_pml_covariance

   pure real(dp) function sample_standard_deviation(values) result(value)
      ! Return the ordinary sample standard deviation.
      real(dp), intent(in) :: values(:)
      real(dp) :: center

      if (size(values) < 2) then
         value = 0.0_dp
         return
      end if
      center = sum(values)/real(size(values), dp)
      value = sqrt(sum((values - center)**2)/real(size(values) - 1, dp))
   end function sample_standard_deviation

   pure function cumulative_sum(values) result(sums)
      ! Return cumulative sums without modifying the input.
      real(dp), intent(in) :: values(:)
      real(dp) :: sums(size(values))
      integer :: i

      if (size(values) < 1) return
      sums(1) = values(1)
      do i = 2, size(values)
         sums(i) = sums(i - 1) + values(i)
      end do
   end function cumulative_sum

   pure function ordinary_difference(values) result(difference)
      ! Return one ordinary first difference.
      real(dp), intent(in) :: values(:)
      real(dp), allocatable :: difference(:)

      if (size(values) < 2) then
         allocate(difference(0))
      else
         difference = values(2:) - values(:size(values) - 1)
      end if
   end function ordinary_difference

end module nsarfima_mod
