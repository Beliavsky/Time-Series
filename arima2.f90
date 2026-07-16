! SPDX-License-Identifier: GPL-3.0-or-later
! SPDX-FileComment: Algorithms translated from the R arima2 package.
! Likelihood inference algorithms translated from the GPL-3 arima2 package.
module arima2_mod
   use kind_mod, only: dp
   use astsa_mod, only: astsa_sarima_fit_t, sarima_fit
   use time_series_random_mod, only: random_uniform, random_standard_normal
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   private

   integer, parameter, public :: arima2_method_css = 1
   integer, parameter, public :: arima2_method_ml = 2
   integer, parameter, public :: arima2_method_css_ml = 3
   integer, parameter, public :: arima2_sampling_dl = 1
   integer, parameter, public :: arima2_sampling_uniform_roots = 2

   type, public :: arima2_roots_t
      ! Polynomial roots and numerical convergence status.
      complex(dp), allocatable :: roots(:)
      integer :: info = 0
   end type

   type, public :: arima2_fit_t
      ! Best random-restart SARIMA fit and likelihood search history.
      type(astsa_sarima_fit_t) :: fit
      real(dp), allocatable :: log_likelihoods(:)
      integer :: starts_attempted = 0
      integer :: improvements = 0
      integer :: info = 0
   end type

   type, public :: arima2_aic_table_t
      ! AIC or AICc values indexed by AR and MA order.
      real(dp), allocatable :: values(:, :)
      integer :: differencing = 0
      integer :: info = 0
   end type

   type, public :: arima2_profile_t
      ! Fixed parameter grid, fitted coefficients, and profile log likelihoods.
      real(dp), allocatable :: parameter_values(:), coefficients(:, :), log_likelihood(:)
      integer :: parameter_index = 0
      integer :: info = 0
   end type

   type, public :: arima2_coefficient_samples_t
      ! Sampled ordinary, seasonal, and optional intercept coefficients.
      real(dp), allocatable :: values(:, :)
      integer :: attempts = 0
      integer :: info = 0
   end type

   public :: arma_polynomial_roots, inverse_roots_to_coefficients
   public :: durbin_levinson_coefficients, arima2_fit_from_starts
   public :: arima2_aic_table, arima2_profile
   public :: sample_arma_coefficients, arima2_fit
   public :: sample_inverse_roots, invert_ma_coefficients

contains

   pure function arma_polynomial_roots(coefficients, moving_average) result(out)
      ! Return roots of an AR or MA lag polynomial.
      real(dp), intent(in) :: coefficients(:)
      logical, intent(in), optional :: moving_average
      type(arima2_roots_t) :: out
      real(dp), allocatable :: polynomial(:)
      logical :: is_ma

      is_ma = .false.
      if (present(moving_average)) is_ma = moving_average
      if (size(coefficients) < 1) then
         out%info = 1
         return
      end if
      allocate(polynomial(0:size(coefficients)))
      polynomial = 0.0_dp
      polynomial(0) = 1.0_dp
      if (is_ma) then
         polynomial(1:) = coefficients
      else
         polynomial(1:) = -coefficients
      end if
      out = polynomial_roots(polynomial)
   end function arma_polynomial_roots

   pure function inverse_roots_to_coefficients(inverse_roots, moving_average) result(coefficients)
      ! Convert inverse polynomial roots to real AR or MA coefficients.
      complex(dp), intent(in) :: inverse_roots(:)
      logical, intent(in), optional :: moving_average
      real(dp), allocatable :: coefficients(:)
      complex(dp), allocatable :: polynomial(:), next_polynomial(:)
      logical :: is_ma
      integer :: order, i

      is_ma = .false.
      if (present(moving_average)) is_ma = moving_average
      allocate(polynomial(0:0))
      polynomial = cmplx(1.0_dp, 0.0_dp, dp)
      do order = 1, size(inverse_roots)
         allocate(next_polynomial(0:order))
         next_polynomial = cmplx(0.0_dp, 0.0_dp, dp)
         do i = 0, order - 1
            next_polynomial(i) = next_polynomial(i) + polynomial(i)
            next_polynomial(i + 1) = next_polynomial(i + 1) - inverse_roots(order)*polynomial(i)
         end do
         polynomial = next_polynomial
         deallocate(next_polynomial)
      end do
      allocate(coefficients(size(inverse_roots)))
      if (is_ma) then
         coefficients = real(polynomial(1:), dp)
      else
         coefficients = -real(polynomial(1:), dp)
      end if
   end function inverse_roots_to_coefficients

   pure function durbin_levinson_coefficients(partial, moving_average) result(coefficients)
      ! Convert sampled partial autocorrelations to stable AR or MA coefficients.
      real(dp), intent(in) :: partial(:)
      logical, intent(in), optional :: moving_average
      real(dp) :: coefficients(size(partial)), previous(size(partial))
      logical :: is_ma
      integer :: order, j

      is_ma = .false.
      if (present(moving_average)) is_ma = moving_average
      coefficients = 0.0_dp
      do order = 1, size(partial)
         previous = coefficients
         coefficients(order) = partial(order)
         do j = 1, order - 1
            coefficients(j) = previous(j) - partial(order)*previous(order - j)
         end do
      end do
      if (is_ma) coefficients = -coefficients
   end function durbin_levinson_coefficients

   pure function arima2_fit_from_starts(series, p, d, q, seasonal_p, seasonal_difference, &
      seasonal_q, season, starts, include_mean, max_repeats, epsilon_tolerance, &
      max_inverse_root, min_inverse_root_distance, max_iterations, tolerance, regressors, &
      estimated, exact_likelihood) result(out)
      ! Select the best exact SARIMA likelihood from caller-supplied restart values.
      real(dp), intent(in) :: series(:), starts(:, :)
      integer, intent(in) :: p, d, q, seasonal_p, seasonal_difference, seasonal_q, season
      logical, intent(in), optional :: include_mean
      integer, intent(in), optional :: max_repeats, max_iterations
      real(dp), intent(in), optional :: epsilon_tolerance, max_inverse_root
      real(dp), intent(in), optional :: min_inverse_root_distance, tolerance
      real(dp), intent(in), optional :: regressors(:, :)
      logical, intent(in), optional :: estimated(:), exact_likelihood
      type(arima2_fit_t) :: out
      type(astsa_sarima_fit_t) :: candidate
      real(dp) :: best_log_likelihood, improvement_tolerance, root_limit, distance_limit
      integer :: parameter_count, repeat_limit, repeats, start_index
      integer :: regression_count
      logical :: use_mean, valid, use_exact

      use_mean = d + seasonal_difference == 0
      if (present(include_mean)) use_mean = include_mean
      regression_count = 0
      if (present(regressors)) regression_count = size(regressors, 2)
      parameter_count = p + q + seasonal_p + seasonal_q + regression_count + merge(1, 0, use_mean)
      use_exact = .true.
      if (present(exact_likelihood)) use_exact = exact_likelihood
      repeat_limit = 10
      if (present(max_repeats)) repeat_limit = max_repeats
      improvement_tolerance = 1.0e-4_dp
      if (present(epsilon_tolerance)) improvement_tolerance = epsilon_tolerance
      root_limit = 1.0_dp
      if (present(max_inverse_root)) root_limit = max_inverse_root
      distance_limit = 0.0_dp
      if (present(min_inverse_root_distance)) distance_limit = min_inverse_root_distance
      if (size(starts, 1) /= parameter_count .or. size(starts, 2) < 1 .or. &
         repeat_limit < 1 .or. root_limit <= 0.0_dp .or. root_limit > 1.0_dp .or. &
         distance_limit < 0.0_dp .or. distance_limit >= 1.0_dp) then
         out%info = 1
         return
      end if
      allocate(out%log_likelihoods(size(starts, 2)))
      out%log_likelihoods = -huge(1.0_dp)
      best_log_likelihood = -huge(1.0_dp)
      repeats = 0
      do start_index = 1, size(starts, 2)
         candidate = sarima_fit(series, p, d, q, seasonal_p, seasonal_difference, seasonal_q, &
            season, initial=starts(:, start_index), include_intercept=use_mean, &
            max_iterations=max_iterations, tolerance=tolerance, regressors=regressors, &
            estimated=estimated, exact_likelihood=use_exact)
         out%starts_attempted = start_index
         if (candidate%likelihood%info == 0 .and. &
            ieee_is_finite(candidate%likelihood%log_likelihood)) then
            out%log_likelihoods(start_index) = candidate%likelihood%log_likelihood
         end if
         valid = candidate%likelihood%info == 0 .and. &
            ieee_is_finite(candidate%likelihood%log_likelihood)
         if (valid .and. start_index > 1) valid = acceptable_roots(candidate%coefficients, &
            p, q, seasonal_p, seasonal_q, root_limit, distance_limit)
         if (valid .and. candidate%likelihood%log_likelihood > &
            best_log_likelihood + improvement_tolerance) then
            out%fit = candidate
            best_log_likelihood = candidate%likelihood%log_likelihood
            out%improvements = out%improvements + 1
            repeats = 0
         else
            repeats = repeats + 1
         end if
         if (repeats >= repeat_limit) exit
      end do
      out%log_likelihoods = out%log_likelihoods(:out%starts_attempted)
      if (best_log_likelihood == -huge(1.0_dp)) out%info = 2
   end function arima2_fit_from_starts

   pure function arima2_aic_table(series, max_ar, differencing, max_ma, corrected, &
      max_iterations, tolerance) result(out)
      ! Fit every nonseasonal ARIMA order and return its AIC or AICc.
      real(dp), intent(in) :: series(:)
      integer, intent(in) :: max_ar, differencing, max_ma
      logical, intent(in), optional :: corrected
      integer, intent(in), optional :: max_iterations
      real(dp), intent(in), optional :: tolerance
      type(arima2_aic_table_t) :: out
      type(astsa_sarima_fit_t) :: fitted
      real(dp), allocatable :: initial(:)
      integer :: p, q, count
      logical :: use_corrected, use_mean

      use_corrected = .false.
      if (present(corrected)) use_corrected = corrected
      if (max_ar < 0 .or. max_ma < 0 .or. differencing < 0) then
         out%info = 1
         return
      end if
      out%differencing = differencing
      allocate(out%values(0:max_ar, 0:max_ma))
      out%values = huge(1.0_dp)
      use_mean = differencing == 0
      do p = 0, max_ar
         do q = 0, max_ma
            count = p + q + merge(1, 0, use_mean)
            allocate(initial(count))
            initial = 0.0_dp
            if (use_mean) initial(count) = sum(series)/real(size(series), dp)
            fitted = sarima_fit(series, p, differencing, q, 0, 0, 0, 1, initial=initial, &
               include_intercept=use_mean, max_iterations=max_iterations, tolerance=tolerance, &
               exact_likelihood=.true.)
            if (fitted%likelihood%info == 0) then
               if (use_corrected) then
                  out%values(p, q) = fitted%likelihood%aicc
               else
                  out%values(p, q) = fitted%likelihood%aic
               end if
            end if
            deallocate(initial)
         end do
      end do
      if (all(out%values == huge(1.0_dp))) out%info = 2
   end function arima2_aic_table

   pure function arima2_profile(series, fitted, parameter_index, parameter_values, &
      max_iterations, tolerance) result(out)
      ! Refit an ARIMA model over a grid with one reported parameter fixed.
      real(dp), intent(in) :: series(:), parameter_values(:)
      type(astsa_sarima_fit_t), intent(in) :: fitted
      integer, intent(in) :: parameter_index
      integer, intent(in), optional :: max_iterations
      real(dp), intent(in), optional :: tolerance
      type(arima2_profile_t) :: out
      type(astsa_sarima_fit_t) :: candidate
      real(dp), allocatable :: initial(:)
      logical, allocatable :: estimated(:)
      integer :: count, grid_index

      count = size(fitted%coefficients)
      if (parameter_index < 1 .or. parameter_index > count .or. size(parameter_values) < 1) then
         out%info = 1
         return
      end if
      out%parameter_index = parameter_index
      out%parameter_values = parameter_values
      allocate(out%coefficients(count, size(parameter_values)))
      allocate(out%log_likelihood(size(parameter_values)), initial(count), estimated(count))
      out%coefficients = 0.0_dp
      out%log_likelihood = -huge(1.0_dp)
      estimated = .true.
      estimated(parameter_index) = .false.
      do grid_index = 1, size(parameter_values)
         initial = fitted%coefficients
         initial(parameter_index) = parameter_values(grid_index)
         candidate = sarima_fit(series, fitted%p, fitted%d, fitted%q, fitted%seasonal_p, &
            fitted%seasonal_difference, fitted%seasonal_q, fitted%season, initial=initial, &
            include_intercept=fitted%includes_intercept, include_drift=fitted%includes_drift, &
            estimated=estimated, transform_parameters=.false., max_iterations=max_iterations, &
            tolerance=tolerance, exact_likelihood=.true.)
         if (candidate%likelihood%info == 0) then
            out%coefficients(:, grid_index) = candidate%coefficients
            out%log_likelihood(grid_index) = candidate%likelihood%log_likelihood
         end if
      end do
      if (all(out%log_likelihood == -huge(1.0_dp))) out%info = 2
   end function arima2_profile

   function sample_arma_coefficients(ar_order, ma_order, seasonal_ar_order, &
      seasonal_ma_order, sample_count, min_inverse_root_distance, include_intercept, &
      intercept_mean, intercept_sd, max_attempts, sampling_method, modulus_bounds) result(out)
      ! Sample stable ARMA coefficients by the arima2 Durbin-Levinson method.
      integer, intent(in) :: ar_order, ma_order, seasonal_ar_order, seasonal_ma_order
      integer, intent(in) :: sample_count
      real(dp), intent(in), optional :: min_inverse_root_distance, intercept_mean, intercept_sd
      logical, intent(in), optional :: include_intercept
      integer, intent(in), optional :: max_attempts
      integer, intent(in), optional :: sampling_method
      real(dp), intent(in), optional :: modulus_bounds(2)
      type(arima2_coefficient_samples_t) :: out
      real(dp), allocatable :: partial(:), coefficients(:)
      real(dp) :: distance_limit, mean_value, sd_value
      integer :: dynamic_count, total_count, limit, sample, offset, order, attempt
      integer :: sampler
      logical :: use_intercept, valid
      real(dp) :: bounds(2)

      distance_limit = 0.0_dp
      if (present(min_inverse_root_distance)) distance_limit = min_inverse_root_distance
      use_intercept = .false.
      if (present(include_intercept)) use_intercept = include_intercept
      mean_value = 0.0_dp
      if (present(intercept_mean)) mean_value = intercept_mean
      sd_value = 0.05_dp
      if (present(intercept_sd)) sd_value = intercept_sd
      limit = 10000
      if (present(max_attempts)) limit = max_attempts
      sampler = arima2_sampling_dl
      if (present(sampling_method)) sampler = sampling_method
      bounds = [0.0_dp, 1.0_dp]
      if (present(modulus_bounds)) bounds = modulus_bounds
      dynamic_count = ar_order + ma_order + seasonal_ar_order + seasonal_ma_order
      total_count = dynamic_count + merge(1, 0, use_intercept)
      if (min(ar_order, ma_order, seasonal_ar_order, seasonal_ma_order) < 0 .or. &
         sample_count < 1 .or. distance_limit < 0.0_dp .or. distance_limit >= 1.0_dp .or. &
         sd_value < 0.0_dp .or. limit < sample_count .or. &
         sampler < arima2_sampling_dl .or. sampler > arima2_sampling_uniform_roots .or. &
         bounds(1) < 0.0_dp .or. bounds(2) > 1.0_dp .or. bounds(1) >= bounds(2)) then
         out%info = 1
         return
      end if
      allocate(out%values(total_count, sample_count), coefficients(dynamic_count))
      out%values = 0.0_dp
      sample = 0
      do attempt = 1, limit
         out%attempts = attempt
         offset = 0
         do order = 1, 4
            select case (order)
            case (1)
               call sampled_block(ar_order, .false., coefficients, offset)
            case (2)
               call sampled_block(ma_order, .true., coefficients, offset)
            case (3)
               call sampled_block(seasonal_ar_order, .false., coefficients, offset)
            case (4)
               call sampled_block(seasonal_ma_order, .true., coefficients, offset)
            end select
         end do
         valid = acceptable_roots(coefficients, ar_order, ma_order, seasonal_ar_order, &
            seasonal_ma_order, 1.0_dp, distance_limit)
         if (.not. valid) cycle
         sample = sample + 1
         if (dynamic_count > 0) out%values(:dynamic_count, sample) = coefficients
         if (use_intercept) out%values(total_count, sample) = mean_value + sd_value*random_standard_normal()
         if (sample == sample_count) return
      end do
      out%info = 2

   contains

      subroutine sampled_block(block_order, moving_average, values, start)
         ! Draw one stable coefficient block from bounded partial autocorrelations.
         integer, intent(in) :: block_order
         logical, intent(in) :: moving_average
         real(dp), intent(inout) :: values(:)
         integer, intent(inout) :: start
         integer :: j
         complex(dp), allocatable :: inverse_roots(:)

         if (block_order == 0) return
         if (sampler == arima2_sampling_dl) then
            allocate(partial(block_order))
            do j = 1, block_order
               partial(j) = -0.99_dp + 1.98_dp*random_uniform()
            end do
            values(start + 1:start + block_order) = &
               durbin_levinson_coefficients(partial, moving_average)
            deallocate(partial)
         else
            inverse_roots = sample_inverse_roots(block_order, bounds(1), bounds(2))
            values(start + 1:start + block_order) = &
               inverse_roots_to_coefficients(inverse_roots, moving_average)
         end if
         start = start + block_order
      end subroutine sampled_block
   end function sample_arma_coefficients

   function sample_inverse_roots(order, minimum_modulus, maximum_modulus) result(inverse_roots)
      ! Sample real roots and conjugate pairs within an inverse-root annulus.
      integer, intent(in) :: order
      real(dp), intent(in) :: minimum_modulus, maximum_modulus
      complex(dp), allocatable :: inverse_roots(:)
      real(dp) :: radius, angle, sign_value
      integer :: index

      if (order < 0 .or. minimum_modulus < 0.0_dp .or. maximum_modulus > 1.0_dp .or. &
         minimum_modulus >= maximum_modulus) then
         allocate(inverse_roots(0))
         return
      end if
      allocate(inverse_roots(order))
      index = 1
      do while (index <= order)
         if (index < order .and. random_uniform() < 1.0_dp - sqrt(0.5_dp)) then
            radius = minimum_modulus + (maximum_modulus - minimum_modulus)*random_uniform()
            angle = acos(-1.0_dp)*random_uniform()
            inverse_roots(index) = radius*cmplx(cos(angle), sin(angle), dp)
            inverse_roots(index + 1) = conjg(inverse_roots(index))
            index = index + 2
         else
            radius = minimum_modulus + (maximum_modulus - minimum_modulus)*random_uniform()
            sign_value = merge(1.0_dp, -1.0_dp, random_uniform() >= 0.5_dp)
            inverse_roots(index) = cmplx(sign_value*radius, 0.0_dp, dp)
            index = index + 1
         end if
      end do
   end function sample_inverse_roots

   pure function invert_ma_coefficients(coefficients) result(inverted)
      ! Reflect MA polynomial roots inside the unit circle to obtain invertibility.
      real(dp), intent(in) :: coefficients(:)
      real(dp), allocatable :: inverted(:)
      type(arima2_roots_t) :: root_result
      complex(dp), allocatable :: roots(:)
      integer :: i

      root_result = arma_polynomial_roots(coefficients, .true.)
      if (root_result%info /= 0) then
         allocate(inverted(0))
         return
      end if
      roots = root_result%roots
      do i = 1, size(roots)
         if (abs(roots(i)) < 1.0_dp) roots(i) = 1.0_dp/roots(i)
      end do
      inverted = inverse_roots_to_coefficients(1.0_dp/roots, .true.)
   end function invert_ma_coefficients

   function arima2_fit(series, p, d, q, seasonal_p, seasonal_difference, seasonal_q, &
      season, max_starts, max_repeats, epsilon_tolerance, max_inverse_root, &
      min_inverse_root_distance, include_mean, max_iterations, tolerance, method, &
      initial, estimated, regressors) result(out)
      ! Fit SARIMA by a conditional baseline followed by random exact-likelihood restarts.
      real(dp), intent(in) :: series(:)
      integer, intent(in) :: p, d, q, seasonal_p, seasonal_difference, seasonal_q, season
      integer, intent(in), optional :: max_starts, max_repeats, max_iterations
      real(dp), intent(in), optional :: epsilon_tolerance, max_inverse_root
      real(dp), intent(in), optional :: min_inverse_root_distance, tolerance
      logical, intent(in), optional :: include_mean
      integer, intent(in), optional :: method
      real(dp), intent(in), optional :: initial(:), regressors(:, :)
      logical, intent(in), optional :: estimated(:)
      type(arima2_fit_t) :: out
      type(astsa_sarima_fit_t) :: baseline
      type(arima2_coefficient_samples_t) :: samples
      real(dp), allocatable :: starts(:, :), initial_values(:)
      real(dp) :: distance_limit
      logical, allocatable :: estimate_mask(:)
      integer :: start_count, dynamic_count, parameter_count, random_count, first_random
      integer :: fitting_method, regression_count, column
      logical :: use_mean

      start_count = 100
      if (present(max_starts)) start_count = max_starts
      use_mean = d + seasonal_difference == 0
      if (present(include_mean)) use_mean = include_mean
      dynamic_count = p + q + seasonal_p + seasonal_q
      regression_count = 0
      if (present(regressors)) regression_count = size(regressors, 2)
      parameter_count = dynamic_count + regression_count + merge(1, 0, use_mean)
      fitting_method = arima2_method_css_ml
      if (present(method)) fitting_method = method
      distance_limit = 0.0_dp
      if (present(min_inverse_root_distance)) distance_limit = min_inverse_root_distance
      if (start_count < 1 .or. fitting_method < arima2_method_css .or. &
         fitting_method > arima2_method_css_ml) then
         out%info = 1
         return
      end if
      allocate(initial_values(parameter_count), starts(parameter_count, start_count))
      initial_values = 0.0_dp
      if (use_mean) initial_values(dynamic_count + 1) = sum(series)/real(size(series), dp)
      if (present(initial)) then
         if (size(initial) /= parameter_count) then
            out%info = 2
            return
         end if
         initial_values = initial
      end if
      allocate(estimate_mask(parameter_count))
      estimate_mask = .true.
      if (present(estimated)) then
         if (size(estimated) /= parameter_count) then
            out%info = 2
            return
         end if
         estimate_mask = estimated
      end if
      baseline = sarima_fit(series, p, d, q, seasonal_p, seasonal_difference, seasonal_q, &
         season, initial=initial_values, include_intercept=use_mean, max_iterations=max_iterations, &
         tolerance=tolerance, regressors=regressors, estimated=estimate_mask, &
         exact_likelihood=.false.)
      if (fitting_method == arima2_method_css) then
         out%fit = baseline
         out%starts_attempted = 1
         allocate(out%log_likelihoods(1))
         out%log_likelihoods(1) = baseline%likelihood%log_likelihood
         if (baseline%likelihood%info /= 0) out%info = 3
         return
      end if
      starts(:, 1) = initial_values
      first_random = 2
      if (fitting_method == arima2_method_css_ml .and. start_count > 1 .and. &
         baseline%likelihood%info == 0) then
         starts(:, 2) = baseline%coefficients
         first_random = 3
      end if
      random_count = start_count - first_random + 1
      if (random_count > 0) then
         samples = sample_arma_coefficients(p, q, seasonal_p, seasonal_q, random_count, &
            distance_limit, .false., max_attempts= &
            max(10000, 1000*random_count))
         if (samples%info /= 0) then
            out%info = 10 + samples%info
            return
         end if
         do column = first_random, start_count
            starts(:, column) = initial_values
            if (dynamic_count > 0) starts(:dynamic_count, column) = &
               samples%values(:, column - first_random + 1)
            where (.not. estimate_mask) starts(:, column) = initial_values
         end do
      end if
      out = arima2_fit_from_starts(series, p, d, q, seasonal_p, seasonal_difference, &
         seasonal_q, season, starts, use_mean, max_repeats, epsilon_tolerance, &
         max_inverse_root, min_inverse_root_distance, max_iterations, tolerance, &
         regressors, estimate_mask, .true.)
   end function arima2_fit

   pure logical function acceptable_roots(coefficients, p, q, seasonal_p, seasonal_q, &
      root_limit, distance_limit) result(acceptable)
      ! Test inverse-root magnitude and AR-MA separation constraints.
      real(dp), intent(in) :: coefficients(:), root_limit, distance_limit
      integer, intent(in) :: p, q, seasonal_p, seasonal_q
      type(arima2_roots_t) :: ar_roots, ma_roots, sar_roots, sma_roots
      integer :: offset

      acceptable = .true.
      offset = 0
      if (p > 0) ar_roots = arma_polynomial_roots(coefficients(1:p))
      offset = offset + p
      if (q > 0) ma_roots = arma_polynomial_roots(coefficients(offset + 1:offset + q), .true.)
      offset = offset + q
      if (seasonal_p > 0) sar_roots = arma_polynomial_roots(coefficients(offset + 1:offset + seasonal_p))
      offset = offset + seasonal_p
      if (seasonal_q > 0) sma_roots = arma_polynomial_roots(coefficients(offset + 1:offset + seasonal_q), .true.)
      if (p > 0) acceptable = acceptable .and. roots_within_limit(ar_roots, root_limit)
      if (q > 0) acceptable = acceptable .and. roots_within_limit(ma_roots, root_limit)
      if (seasonal_p > 0) acceptable = acceptable .and. roots_within_limit(sar_roots, root_limit)
      if (seasonal_q > 0) acceptable = acceptable .and. roots_within_limit(sma_roots, root_limit)
      if (distance_limit > 0.0_dp .and. p > 0 .and. q > 0) &
         acceptable = acceptable .and. inverse_root_distance(ar_roots, ma_roots) > distance_limit
      if (distance_limit > 0.0_dp .and. seasonal_p > 0 .and. seasonal_q > 0) &
         acceptable = acceptable .and. inverse_root_distance(sar_roots, sma_roots) > distance_limit
   end function acceptable_roots

   pure logical function roots_within_limit(root_result, limit) result(acceptable)
      ! Test that converged inverse roots do not exceed a limit.
      type(arima2_roots_t), intent(in) :: root_result
      real(dp), intent(in) :: limit
      acceptable = root_result%info == 0
      if (acceptable) acceptable = maxval(1.0_dp/abs(root_result%roots)) <= limit
   end function roots_within_limit

   pure real(dp) function inverse_root_distance(first, second) result(distance)
      ! Return the minimum distance between two sets of inverse roots.
      type(arima2_roots_t), intent(in) :: first, second
      integer :: i, j

      distance = huge(1.0_dp)
      if (first%info /= 0 .or. second%info /= 0) then
         distance = 0.0_dp
         return
      end if
      do i = 1, size(first%roots)
         do j = 1, size(second%roots)
            distance = min(distance, abs(1.0_dp/first%roots(i) - 1.0_dp/second%roots(j)))
         end do
      end do
   end function inverse_root_distance

   pure function polynomial_roots(polynomial) result(out)
      ! Compute complex polynomial roots by Durand-Kerner iteration.
      real(dp), intent(in) :: polynomial(0:)
      type(arima2_roots_t) :: out
      complex(dp), allocatable :: previous(:)
      complex(dp) :: numerator, denominator
      real(dp) :: radius, change, angle
      integer :: degree, iteration, i, j, power

      degree = ubound(polynomial, 1)
      if (degree < 1 .or. abs(polynomial(degree)) <= tiny(1.0_dp)) then
         out%info = 1
         return
      end if
      allocate(out%roots(degree), previous(degree))
      if (degree == 1) then
         out%roots(1) = cmplx(-polynomial(0)/polynomial(1), 0.0_dp, dp)
         return
      end if
      radius = 1.0_dp + maxval(abs(polynomial(:degree - 1)/polynomial(degree)))
      do i = 1, degree
         angle = 2.0_dp*acos(-1.0_dp)*real(i - 1, dp)/real(degree, dp) + 0.1_dp
         out%roots(i) = radius*cmplx(cos(angle), sin(angle), dp)
      end do
      do iteration = 1, 2000
         previous = out%roots
         change = 0.0_dp
         do i = 1, degree
            numerator = cmplx(polynomial(degree), 0.0_dp, dp)
            do power = degree - 1, 0, -1
               numerator = numerator*previous(i) + polynomial(power)
            end do
            denominator = cmplx(1.0_dp, 0.0_dp, dp)
            do j = 1, degree
               if (j /= i) denominator = denominator*(previous(i) - previous(j))
            end do
            if (abs(denominator) <= tiny(1.0_dp)) then
               out%info = 2
               return
            end if
            out%roots(i) = previous(i) - numerator/denominator
            change = max(change, abs(out%roots(i) - previous(i)))
         end do
         if (change <= 100.0_dp*epsilon(1.0_dp)) return
      end do
      out%info = 3
   end function polynomial_roots
end module arima2_mod
