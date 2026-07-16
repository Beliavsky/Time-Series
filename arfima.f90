! SPDX-License-Identifier: MIT
! SPDX-FileComment: Algorithms translated from the R arfima package.
! Long-memory covariance and simulation algorithms translated from arfima.
module arfima_mod
   use kind_mod, only: dp
   use itsmr_mod, only: itsmr_arma_model_t, arma_acvf, regularized_gamma_q
   use arima2_mod, only: arima2_roots_t, arma_polynomial_roots, &
      durbin_levinson_coefficients
   use time_series_optimization_mod, only: optimization_result_t, bfgs_minimize_fd, &
      finite_difference_hessian
   use time_series_linalg_mod, only: invert_matrix, cholesky_lower
   use time_series_random_mod, only: set_random_seed, random_uniform, &
      random_standard_normal, &
      multivariate_normal_from_standard
   use time_series_stats_mod, only: normal_quantile
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   private

   integer, parameter, public :: arfima_long_memory_none = 0
   integer, parameter, public :: arfima_long_memory_fdwn = 1
   integer, parameter, public :: arfima_long_memory_fgn = 2
   integer, parameter, public :: arfima_long_memory_hd = 3
   integer, parameter, public :: arfima_likelihood_exact = 1
   integer, parameter, public :: arfima_likelihood_css = 2
   integer, parameter, public :: arfima_weed_none = 0
   integer, parameter, public :: arfima_weed_operator = 1
   integer, parameter, public :: arfima_weed_pacf = 2
   integer, parameter, public :: arfima_weed_both = 3

   type, public :: arfima_model_t
      ! Nonseasonal and seasonal ARFIMA parameters in Box-Jenkins convention.
      real(dp), allocatable :: ar(:), theta(:), seasonal_ar(:), seasonal_theta(:)
      real(dp) :: long_memory_parameter = 0.0_dp
      real(dp) :: seasonal_long_memory_parameter = 0.0_dp
      real(dp) :: innovation_variance = 1.0_dp
      real(dp) :: mean = 0.0_dp
      integer :: long_memory_type = arfima_long_memory_fdwn
      integer :: seasonal_long_memory_type = arfima_long_memory_fdwn
      integer :: difference_order = 0
      integer :: seasonal_difference_order = 0
      integer :: period = 0
      integer :: info = 0
   end type arfima_model_t

   type, public :: arfima_acvf_t
      ! Theoretical autocovariances at lags zero through max_lag.
      real(dp), allocatable :: covariance(:)
      integer :: max_lag = 0
      integer :: info = 0
   end type arfima_acvf_t

   type, public :: arfima_simulation_t
      ! Covariance-driven simulation and prediction-error variances.
      real(dp), allocatable :: series(:), stationary_series(:), covariance(:)
      real(dp), allocatable :: innovation_variance(:)
      integer :: info = 0
   end type arfima_simulation_t

   type, public :: arfima_likelihood_t
      ! Profiled Gaussian likelihood, innovations, and conditional variances.
      real(dp), allocatable :: innovation(:), prediction_variance(:)
      real(dp) :: mean = 0.0_dp
      real(dp) :: innovation_variance = 0.0_dp
      real(dp) :: sum_squares = 0.0_dp
      real(dp) :: log_likelihood = -huge(1.0_dp)
      integer :: observations = 0
      integer :: method = arfima_likelihood_exact
      integer :: info = 0
   end type arfima_likelihood_t

   type, public :: arfima_fit_t
      ! Estimated model, likelihood, inference, and information criteria.
      type(arfima_model_t) :: model
      type(arfima_likelihood_t) :: likelihood
      real(dp), allocatable :: unconstrained_parameters(:)
      real(dp), allocatable :: parameter_covariance(:, :), standard_error(:)
      logical, allocatable :: fixed_parameters(:)
      real(dp) :: aic = huge(1.0_dp)
      real(dp) :: aicc = huge(1.0_dp)
      real(dp) :: bic = huge(1.0_dp)
      integer :: iterations = 0
      integer :: parameter_count = 0
      integer :: estimated_parameter_count = 0
      integer :: info = 0
      logical :: converged = .false.
   end type arfima_fit_t

   type, public :: arfima_multifit_t
      ! Multiple local likelihood modes and their originating starts.
      type(arfima_fit_t), allocatable :: modes(:)
      real(dp), allocatable :: starting_parameters(:, :), weight(:)
      integer, allocatable :: starting_index(:), long_memory_type(:)
      integer :: attempted = 0
      integer :: converged = 0
      integer :: info = 0
      logical :: weeded = .false.
   end type arfima_multifit_t

   type, public :: arfima_forecast_t
      ! Conditional forecasts, uncertainty, intervals, and optional paths.
      real(dp), allocatable :: mean(:), standard_error(:), covariance(:, :)
      real(dp), allocatable :: lower(:), upper(:), paths(:, :)
      real(dp), allocatable :: stationary_mean(:), stationary_covariance(:, :)
      real(dp), allocatable :: regression_mean(:)
      real(dp) :: level = 0.95_dp
      integer :: horizon = 0
      integer :: simulations = 0
      integer :: info = 0
   end type arfima_forecast_t

   type, public :: arfima_transfer_t
      ! One dynamic-regression transfer function with an integer delay.
      real(dp), allocatable :: denominator(:), numerator(:)
      integer :: delay = 0
      integer :: info = 0
   end type arfima_transfer_t

   type, public :: arfima_regression_fit_t
      ! Joint ARFIMA, regression, and transfer-function estimates.
      type(arfima_model_t) :: model
      type(arfima_likelihood_t) :: likelihood
      type(arfima_transfer_t), allocatable :: transfers(:)
      real(dp), allocatable :: coefficients(:), unconstrained_parameters(:)
      real(dp), allocatable :: parameter_covariance(:, :), standard_error(:)
      logical, allocatable :: fixed_parameters(:)
      real(dp) :: aic = huge(1.0_dp)
      real(dp) :: aicc = huge(1.0_dp)
      real(dp) :: bic = huge(1.0_dp)
      integer :: iterations = 0
      integer :: parameter_count = 0
      integer :: estimated_parameter_count = 0
      integer :: info = 0
      logical :: converged = .false.
   end type arfima_regression_fit_t

   type, public :: arfima_identifiability_t
      ! Stability, common-factor distance, and overall identifiability.
      real(dp) :: minimum_root_distance = huge(1.0_dp)
      logical :: causal = .false.
      logical :: invertible = .false.
      logical :: seasonal_causal = .false.
      logical :: seasonal_invertible = .false.
      logical :: transfer_stable = .false.
      logical :: identifiable = .false.
      integer :: info = 0
   end type arfima_identifiability_t

   type, public :: arfima_information_t
      ! Per-observation Fisher information and optional sample covariance.
      real(dp), allocatable :: information(:, :), covariance(:, :)
      integer :: parameters = 0
      integer :: resolution = 0
      integer :: observations = 0
      integer :: info = 0
      logical :: exact = .true.
      logical :: positive_definite = .false.
   end type arfima_information_t

   type, public :: arfima_diagnostics_t
      ! Fitted values, innovations, portmanteau tests, and normal Q-Q data.
      real(dp), allocatable :: fitted(:), residuals(:), standardized_residuals(:)
      real(dp), allocatable :: prediction_variance(:), regression_residuals(:)
      real(dp), allocatable :: residual_acf(:), squared_residual_acf(:)
      real(dp), allocatable :: ljung_box(:), ljung_box_p_value(:)
      real(dp), allocatable :: squared_ljung_box(:), squared_ljung_box_p_value(:)
      real(dp), allocatable :: qq_sample(:), qq_theoretical(:)
      integer, allocatable :: lags(:), degrees_of_freedom(:)
      integer :: observations = 0
      integer :: info = 0
   end type arfima_diagnostics_t

   type, public :: arfima_covariance_t
      ! Observed and expected covariance and correlation matrices.
      real(dp), allocatable :: observed(:, :), expected(:, :)
      real(dp), allocatable :: observed_correlation(:, :)
      real(dp), allocatable :: expected_correlation(:, :)
      integer :: info = 0
      logical :: expected_available = .false.
   end type arfima_covariance_t

   public :: arfima_model
   public :: arfima_fractional_weights, arfima_psi_weights, arfima_seasonal_psi_weights
   public :: arfima_pi_weights
   public :: arfima_fdwn_acvf, arfima_fgn_acvf, arfima_hd_acvf
   public :: arfima_farma_acvf, arfima_model_acvf
   public :: arfima_integrate, arfima_durbin_levinson_simulate
   public :: arfima_simulate_from_innovations, arfima_simulate
   public :: arfima_likelihood, arfima_fit, arfima_fit_fixed
   public :: arfima_fit_modes, arfima_multistart_fit, arfima_select_long_memory
   public :: arfima_fit_modes_fixed, arfima_multistart_fit_fixed
   public :: arfima_mode_distance, arfima_mode_distances, arfima_weed_modes
   public :: arfima_best_modes, arfima_remove_mode
   public :: arfima_mode_forecast, arfima_mode_information
   public :: arfima_mode_identifiability
   public :: arfima_diagnostics, arfima_mode_diagnostics
   public :: arfima_fit_covariance, arfima_average_forecast
   public :: arfima_forecast, arfima_forecast_paths
   public :: arfima_transfer, arfima_transfer_response
   public :: arfima_regression_likelihood, arfima_regression_fit
   public :: arfima_regression_fit_fixed
   public :: arfima_regression_forecast
   public :: arfima_regression_simulate_from_innovations, arfima_regression_simulate
   public :: arfima_identifiability
   public :: arfima_fisher_information, arfima_fit_information

   interface arfima_forecast
      module procedure arfima_forecast_model
      module procedure arfima_forecast_fit
   end interface arfima_forecast

   interface arfima_forecast_paths
      module procedure arfima_forecast_model_paths
      module procedure arfima_forecast_fit_paths
   end interface arfima_forecast_paths

   interface arfima_diagnostics
      module procedure arfima_fit_diagnostics
      module procedure arfima_regression_diagnostics
   end interface arfima_diagnostics

contains

   pure function arfima_model(ar, theta, d, seasonal_ar, seasonal_theta, seasonal_d, &
      period, innovation_variance, mean, difference_order, seasonal_difference_order, &
      long_memory_type, seasonal_long_memory_type) result(model)
      ! Construct an ARFIMA model with allocated coefficient arrays.
      real(dp), intent(in) :: ar(:), theta(:), d
      real(dp), intent(in), optional :: seasonal_ar(:), seasonal_theta(:), seasonal_d
      real(dp), intent(in), optional :: innovation_variance, mean
      integer, intent(in), optional :: period, difference_order, seasonal_difference_order
      integer, intent(in), optional :: long_memory_type, seasonal_long_memory_type
      type(arfima_model_t) :: model

      model%ar = ar
      model%theta = theta
      allocate(model%seasonal_ar(0), model%seasonal_theta(0))
      if (present(seasonal_ar)) model%seasonal_ar = seasonal_ar
      if (present(seasonal_theta)) model%seasonal_theta = seasonal_theta
      model%long_memory_parameter = d
      if (present(seasonal_d)) model%seasonal_long_memory_parameter = seasonal_d
      if (present(period)) model%period = period
      if (present(innovation_variance)) model%innovation_variance = innovation_variance
      if (present(mean)) model%mean = mean
      if (present(difference_order)) model%difference_order = difference_order
      if (present(seasonal_difference_order)) then
         model%seasonal_difference_order = seasonal_difference_order
      end if
      if (present(long_memory_type)) model%long_memory_type = long_memory_type
      if (present(seasonal_long_memory_type)) then
         model%seasonal_long_memory_type = seasonal_long_memory_type
      end if
      if (.not. valid_model(model)) model%info = 1
   end function arfima_model

   pure function arfima_transfer(denominator, numerator, delay) result(transfer)
      ! Construct and validate one dynamic-regression transfer function.
      real(dp), intent(in) :: denominator(:), numerator(:)
      integer, intent(in), optional :: delay
      type(arfima_transfer_t) :: transfer

      transfer%denominator = denominator
      transfer%numerator = numerator
      if (present(delay)) transfer%delay = delay
      if (.not. valid_transfer(transfer)) transfer%info = 1
   end function arfima_transfer

   pure function arfima_transfer_response(regressor, transfer) result(response)
      ! Filter one regressor through a delayed rational transfer function.
      real(dp), intent(in) :: regressor(:)
      type(arfima_transfer_t), intent(in) :: transfer
      real(dp), allocatable :: response(:)
      integer :: i, j, start

      if (.not. valid_transfer(transfer) .or. &
         .not. all(ieee_is_finite(regressor))) then
         allocate(response(0))
         return
      end if
      allocate(response(size(regressor)))
      response = 0.0_dp
      start = max(size(transfer%denominator) + 1, &
         size(transfer%numerator) + transfer%delay)
      do i = start, size(regressor)
         do j = 1, size(transfer%denominator)
            response(i) = response(i) + transfer%denominator(j)*response(i - j)
         end do
         response(i) = response(i) + &
            transfer%numerator(1)*regressor(i - transfer%delay)
         do j = 2, size(transfer%numerator)
            response(i) = response(i) - &
               transfer%numerator(j)*regressor(i - transfer%delay - j + 1)
         end do
      end do
   end function arfima_transfer_response

   pure function arfima_regression_likelihood(series, model, regressors, coefficients, &
      transfers, transfer_regressors, method, estimate_mean) result(out)
      ! Evaluate ARFIMA likelihood after removing static and dynamic regressors.
      real(dp), intent(in) :: series(:), regressors(:, :), coefficients(:)
      type(arfima_model_t), intent(in) :: model
      type(arfima_transfer_t), intent(in) :: transfers(:)
      real(dp), intent(in) :: transfer_regressors(:, :)
      integer, intent(in), optional :: method
      logical, intent(in), optional :: estimate_mean
      type(arfima_likelihood_t) :: out
      real(dp), allocatable :: effect(:)
      integer :: selected_method
      logical :: profile_mean

      selected_method = arfima_likelihood_exact
      if (present(method)) selected_method = method
      profile_mean = .true.
      if (present(estimate_mean)) profile_mean = estimate_mean
      if (.not. valid_regression_inputs(series, regressors, coefficients, transfers, &
         transfer_regressors)) then
         out%info = 1
         return
      end if
      effect = regression_effect(regressors, coefficients, transfers, transfer_regressors)
      out = arfima_likelihood(series - effect, model, selected_method, profile_mean)
   end function arfima_regression_likelihood

   pure function arfima_identifiability(model, transfers, tolerance) result(out)
      ! Diagnose stability and AR/MA common polynomial factors.
      type(arfima_model_t), intent(in) :: model
      type(arfima_transfer_t), intent(in), optional :: transfers(:)
      real(dp), intent(in), optional :: tolerance
      type(arfima_identifiability_t) :: out
      type(arima2_roots_t) :: ar_roots, ma_roots
      real(dp), allocatable :: combined_ar(:), combined_ma(:)
      real(dp) :: selected_tolerance
      integer :: i, j

      selected_tolerance = sqrt(epsilon(1.0_dp))
      if (present(tolerance)) selected_tolerance = tolerance
      if (selected_tolerance < 0.0_dp .or. .not. allocated(model%ar) .or. &
         .not. allocated(model%theta) .or. .not. allocated(model%seasonal_ar) .or. &
         .not. allocated(model%seasonal_theta)) then
         out%info = 1
         return
      end if
      out%causal = stable_coefficients(model%ar)
      out%invertible = stable_coefficients(model%theta)
      out%seasonal_causal = stable_coefficients(model%seasonal_ar)
      out%seasonal_invertible = stable_coefficients(model%seasonal_theta)
      out%transfer_stable = .true.
      if (present(transfers)) then
         do i = 1, size(transfers)
            if (.not. valid_transfer(transfers(i))) out%transfer_stable = .false.
         end do
      end if
      combined_ar = combined_operator_coefficients(model%ar, model%seasonal_ar, model%period)
      combined_ma = combined_operator_coefficients(model%theta, model%seasonal_theta, &
         model%period)
      if (size(combined_ar) > 0 .and. size(combined_ma) > 0) then
         ar_roots = arma_polynomial_roots(combined_ar)
         ma_roots = arma_polynomial_roots(combined_ma)
         if (ar_roots%info /= 0 .or. ma_roots%info /= 0) then
            out%info = 2
            return
         end if
         do i = 1, size(ar_roots%roots)
            do j = 1, size(ma_roots%roots)
               out%minimum_root_distance = min(out%minimum_root_distance, &
                  abs(ar_roots%roots(i) - ma_roots%roots(j)))
            end do
         end do
      end if
      out%identifiable = out%causal .and. out%invertible .and. &
         out%seasonal_causal .and. out%seasonal_invertible .and. &
         out%transfer_stable .and. valid_long_memory(model%long_memory_type, &
         model%long_memory_parameter) .and. &
         (model%period == 0 .or. valid_long_memory(model%seasonal_long_memory_type, &
         model%seasonal_long_memory_parameter)) .and. &
         out%minimum_root_distance > selected_tolerance
   end function arfima_identifiability

   pure function arfima_fractional_weights(d, max_lag, period) result(weights)
      ! Expand (1-B**period)**d through the requested ordinary lag.
      real(dp), intent(in) :: d
      integer, intent(in) :: max_lag
      integer, intent(in), optional :: period
      real(dp), allocatable :: weights(:)
      integer :: k, selected_period

      if (max_lag < 0 .or. .not. ieee_is_finite(d)) then
         allocate(weights(0))
         return
      end if
      selected_period = 1
      if (present(period)) selected_period = period
      if (selected_period < 1) then
         allocate(weights(0))
         return
      end if
      allocate(weights(0:max_lag))
      weights = 0.0_dp
      weights(0) = 1.0_dp
      do k = 1, max_lag/selected_period
         weights(k*selected_period) = weights((k - 1)*selected_period)* &
            (real(k - 1, dp) - d)/real(k, dp)
      end do
   end function arfima_fractional_weights

   pure function arfima_psi_weights(ar, theta, d, max_lag) result(weights)
      ! Compute FARMA impulse weights using arfima's Box-Jenkins MA sign.
      real(dp), intent(in) :: ar(:), theta(:), d
      integer, intent(in) :: max_lag
      real(dp), allocatable :: weights(:)
      real(dp), allocatable :: fractional(:), short_memory(:)
      integer :: i, j

      if (max_lag < 0 .or. .not. ieee_is_finite(d) .or. &
         .not. all(ieee_is_finite(ar)) .or. .not. all(ieee_is_finite(theta))) then
         allocate(weights(0))
         return
      end if
      allocate(short_memory(0:max_lag))
      short_memory = 0.0_dp
      short_memory(0) = 1.0_dp
      do i = 1, max_lag
         do j = 1, min(i, size(ar))
            short_memory(i) = short_memory(i) + ar(j)*short_memory(i - j)
         end do
         if (i <= size(theta)) short_memory(i) = short_memory(i) - theta(i)
      end do
      fractional = arfima_fractional_weights(-d, max_lag)
      allocate(weights(0:max_lag))
      weights = polynomial_product_truncated(short_memory, fractional, max_lag)
   end function arfima_psi_weights

   pure function arfima_seasonal_psi_weights(ar, theta, d, seasonal_ar, &
      seasonal_theta, seasonal_d, period, max_lag) result(weights)
      ! Combine ordinary and seasonal FARMA impulse weights.
      real(dp), intent(in) :: ar(:), theta(:), d
      real(dp), intent(in) :: seasonal_ar(:), seasonal_theta(:), seasonal_d
      integer, intent(in) :: period, max_lag
      real(dp), allocatable :: weights(:)
      real(dp), allocatable :: ordinary(:), seasonal_base(:), seasonal(:)
      integer :: lag

      if (period < 2 .or. max_lag < 0) then
         allocate(weights(0))
         return
      end if
      ordinary = arfima_psi_weights(ar, theta, d, max_lag)
      seasonal_base = arfima_psi_weights(seasonal_ar, seasonal_theta, &
         seasonal_d, max_lag/period)
      if (size(ordinary) == 0 .or. size(seasonal_base) == 0) then
         allocate(weights(0))
         return
      end if
      allocate(seasonal(0:max_lag))
      seasonal = 0.0_dp
      do lag = 0, max_lag/period
         seasonal(lag*period) = seasonal_base(lag + 1)
      end do
      allocate(weights(0:max_lag))
      weights = polynomial_product_truncated(ordinary, seasonal, max_lag)
   end function arfima_seasonal_psi_weights

   pure function arfima_pi_weights(model, max_lag) result(weights)
      ! Expand the inverse FD-ARFIMA filter through the requested lag.
      type(arfima_model_t), intent(in) :: model
      integer, intent(in) :: max_lag
      real(dp), allocatable :: weights(:)
      real(dp), allocatable :: numerator(:), factor(:), ma(:)
      integer :: lag, j

      if (max_lag < 0 .or. .not. valid_model(model) .or. &
         (model%long_memory_type /= arfima_long_memory_none .and. &
         model%long_memory_type /= arfima_long_memory_fdwn) .or. &
         (model%period > 0 .and. &
         model%seasonal_long_memory_type /= arfima_long_memory_none .and. &
         model%seasonal_long_memory_type /= arfima_long_memory_fdwn)) then
         allocate(weights(0))
         return
      end if
      allocate(numerator(0:max_lag))
      numerator = 0.0_dp
      numerator(0) = 1.0_dp
      factor = coefficient_lag_polynomial(model%ar, 1, max_lag)
      numerator = polynomial_product_truncated(numerator, factor, max_lag)
      if (model%period > 0) then
         factor = coefficient_lag_polynomial(model%seasonal_ar, model%period, max_lag)
         numerator = polynomial_product_truncated(numerator, factor, max_lag)
      end if
      if (model%long_memory_type == arfima_long_memory_fdwn) then
         factor = arfima_fractional_weights(model%long_memory_parameter, max_lag)
         numerator = polynomial_product_truncated(numerator, factor, max_lag)
      end if
      if (model%period > 0 .and. &
         model%seasonal_long_memory_type == arfima_long_memory_fdwn) then
         factor = arfima_fractional_weights(model%seasonal_long_memory_parameter, &
            max_lag, model%period)
         numerator = polynomial_product_truncated(numerator, factor, max_lag)
      end if
      if (model%difference_order > 0) then
         factor = arfima_fractional_weights(real(model%difference_order, dp), max_lag)
         numerator = polynomial_product_truncated(numerator, factor, max_lag)
      end if
      if (model%seasonal_difference_order > 0) then
         factor = arfima_fractional_weights(real(model%seasonal_difference_order, dp), &
            max_lag, model%period)
         numerator = polynomial_product_truncated(numerator, factor, max_lag)
      end if
      allocate(ma(0:max_lag))
      ma = coefficient_lag_polynomial(model%theta, 1, max_lag)
      if (model%period > 0) then
         factor = coefficient_lag_polynomial(model%seasonal_theta, model%period, max_lag)
         ma = polynomial_product_truncated(ma, factor, max_lag)
      end if
      allocate(weights(0:max_lag))
      weights = 0.0_dp
      do lag = 0, max_lag
         weights(lag) = numerator(lag)
         do j = 1, lag
            weights(lag) = weights(lag) - ma(j)*weights(lag - j)
         end do
      end do
   end function arfima_pi_weights

   pure function arfima_fdwn_acvf(d, max_lag) result(out)
      ! Compute fractional-differenced white-noise autocovariances.
      real(dp), intent(in) :: d
      integer, intent(in) :: max_lag
      type(arfima_acvf_t) :: out
      integer :: lag

      out%max_lag = max_lag
      if (max_lag < 0 .or. d <= -1.0_dp .or. d >= 0.5_dp .or. &
         .not. ieee_is_finite(d)) then
         out%info = 1
         return
      end if
      allocate(out%covariance(0:max_lag))
      out%covariance(0) = exp(log_gamma(1.0_dp - 2.0_dp*d) - &
         2.0_dp*log_gamma(1.0_dp - d))
      do lag = 1, max_lag
         out%covariance(lag) = ((real(lag - 1, dp) + d)/ &
            (real(lag, dp) - d))*out%covariance(lag - 1)
      end do
   end function arfima_fdwn_acvf

   pure function arfima_fgn_acvf(hurst, max_lag) result(out)
      ! Compute unit-variance fractional Gaussian-noise autocovariances.
      real(dp), intent(in) :: hurst
      integer, intent(in) :: max_lag
      type(arfima_acvf_t) :: out
      real(dp) :: exponent
      integer :: lag

      out%max_lag = max_lag
      if (max_lag < 0 .or. hurst <= 0.0_dp .or. hurst >= 1.0_dp .or. &
         .not. ieee_is_finite(hurst)) then
         out%info = 1
         return
      end if
      exponent = 2.0_dp*hurst
      allocate(out%covariance(0:max_lag))
      out%covariance(0) = 1.0_dp
      do lag = 1, max_lag
         out%covariance(lag) = 0.5_dp*(real(lag + 1, dp)**exponent - &
            2.0_dp*real(lag, dp)**exponent + real(lag - 1, dp)**exponent)
      end do
   end function arfima_fgn_acvf

   pure function arfima_hd_acvf(alpha, max_lag, zeta_terms) result(out)
      ! Compute normalized hyperbolic-decay autocovariances.
      real(dp), intent(in) :: alpha
      integer, intent(in) :: max_lag
      integer, intent(in), optional :: zeta_terms
      type(arfima_acvf_t) :: out
      real(dp) :: scale, zeta_value
      integer :: lag, terms

      out%max_lag = max_lag
      terms = 20
      if (present(zeta_terms)) terms = zeta_terms
      if (max_lag < 0 .or. alpha <= 0.0_dp .or. alpha >= 3.0_dp .or. &
         terms < 2 .or. .not. ieee_is_finite(alpha)) then
         out%info = 1
         return
      end if
      if (abs(alpha - 1.0_dp) <= epsilon(1.0_dp)) then
         scale = 0.0_dp
      else
         zeta_value = borwein_zeta(alpha, terms)
         if (.not. ieee_is_finite(zeta_value) .or. zeta_value == 0.0_dp) then
            out%info = 2
            return
         end if
         scale = -0.5_dp/zeta_value
      end if
      allocate(out%covariance(0:max_lag))
      out%covariance(0) = 1.0_dp
      do lag = 1, max_lag
         out%covariance(lag) = scale*real(lag, dp)**(-alpha)
      end do
   end function arfima_hd_acvf

   pure function arfima_farma_acvf(ar, theta, d, max_lag, innovation_variance, &
      truncation_lag) result(out)
      ! Combine FDWN and ARMA covariances by truncated spectral convolution.
      real(dp), intent(in) :: ar(:), theta(:), d, innovation_variance
      integer, intent(in) :: max_lag
      integer, intent(in), optional :: truncation_lag
      type(arfima_acvf_t) :: out
      type(arfima_acvf_t) :: fractional
      type(itsmr_arma_model_t) :: arma_model
      real(dp), allocatable :: arma_covariance(:)
      integer :: truncation

      out%max_lag = max_lag
      truncation = max(max_lag, 256)
      if (present(truncation_lag)) truncation = truncation_lag
      if (max_lag < 0 .or. truncation < max_lag .or. innovation_variance < 0.0_dp .or. &
         .not. ieee_is_finite(innovation_variance) .or. &
         .not. all(ieee_is_finite(ar)) .or. .not. all(ieee_is_finite(theta))) then
         out%info = 1
         return
      end if
      fractional = arfima_fdwn_acvf(d, truncation)
      if (fractional%info /= 0) then
         out%info = fractional%info
         return
      end if
      allocate(arma_model%ar(size(ar)), arma_model%ma(size(theta)))
      arma_model%ar = ar
      arma_model%ma = -theta
      arma_model%innovation_variance = 1.0_dp
      arma_covariance = arma_acvf(arma_model, truncation)
      if (.not. all(ieee_is_finite(arma_covariance))) then
         out%info = 2
         return
      end if
      allocate(out%covariance(0:max_lag))
      out%covariance = innovation_variance* &
         symmetric_covariance_product(fractional%covariance, arma_covariance, max_lag)
   end function arfima_farma_acvf

   pure function arfima_model_acvf(model, max_lag, truncation_lag) result(out)
      ! Compute a seasonal ARFIMA covariance by multiplying component spectra.
      type(arfima_model_t), intent(in) :: model
      integer, intent(in) :: max_lag
      integer, intent(in), optional :: truncation_lag
      type(arfima_acvf_t) :: out
      type(arfima_acvf_t) :: ordinary, seasonal_base
      real(dp), allocatable :: seasonal(:)
      integer :: lag, truncation

      out%max_lag = max_lag
      if (model%period >= 2) then
         truncation = 2*max(128, next_power_of_two(max_lag))* &
            next_power_of_two(model%period)
      else
         truncation = max(max_lag, 256)
      end if
      if (present(truncation_lag)) truncation = truncation_lag
      if (.not. valid_model(model) .or. max_lag < 0 .or. truncation < max_lag) then
         out%info = 1
         return
      end if
      ordinary = component_acvf(model%ar, model%theta, model%long_memory_type, &
         model%long_memory_parameter, truncation)
      if (ordinary%info /= 0) then
         out%info = ordinary%info
         return
      end if
      if (model%period == 0) then
         allocate(out%covariance(0:max_lag))
         out%covariance = model%innovation_variance*ordinary%covariance(:max_lag)
         return
      end if

      seasonal_base = component_acvf(model%seasonal_ar, model%seasonal_theta, &
         model%seasonal_long_memory_type, model%seasonal_long_memory_parameter, &
         truncation/model%period)
      if (seasonal_base%info /= 0) then
         out%info = seasonal_base%info
         return
      end if
      allocate(seasonal(0:truncation))
      seasonal = 0.0_dp
      do lag = 0, truncation/model%period
         seasonal(lag*model%period) = seasonal_base%covariance(lag)
      end do
      allocate(out%covariance(0:max_lag))
      out%covariance = model%innovation_variance* &
         symmetric_covariance_product(ordinary%covariance, seasonal, max_lag)
   end function arfima_model_acvf

   pure function arfima_integrate(series, initial_values, difference_order, &
      seasonal_difference_order, period) result(integrated)
      ! Reverse ordinary and seasonal integer differencing using initial values.
      real(dp), intent(in) :: series(:), initial_values(:)
      integer, intent(in) :: difference_order, seasonal_difference_order, period
      real(dp), allocatable :: integrated(:)
      real(dp), allocatable :: ordinary(:), seasonal(:), coefficient(:)
      integer :: active, i, lag

      active = difference_order + period*seasonal_difference_order
      if (difference_order < 0 .or. seasonal_difference_order < 0 .or. &
         (seasonal_difference_order > 0 .and. period < 2) .or. &
         size(initial_values) /= active .or. .not. all(ieee_is_finite(series)) .or. &
         .not. all(ieee_is_finite(initial_values))) then
         allocate(integrated(0))
         return
      end if
      allocate(integrated(size(series) + active))
      if (active == 0) then
         integrated = series
         return
      end if
      integrated(:active) = initial_values
      allocate(ordinary(0:active), seasonal(0:active), coefficient(0:active))
      ordinary = arfima_fractional_weights(real(difference_order, dp), active)
      seasonal = 0.0_dp
      seasonal(0) = 1.0_dp
      if (seasonal_difference_order > 0) then
         seasonal = arfima_fractional_weights(real(seasonal_difference_order, dp), &
            active, period)
      end if
      coefficient = polynomial_product_truncated(ordinary, seasonal, active)
      do i = active + 1, size(integrated)
         integrated(i) = series(i - active)
         do lag = 1, active
            integrated(i) = integrated(i) - coefficient(lag)*integrated(i - lag)
         end do
      end do
   end function arfima_integrate

   pure function arfima_durbin_levinson_simulate(innovations, covariance) result(out)
      ! Simulate a stationary process from its autocovariances and innovations.
      real(dp), intent(in) :: innovations(:)
      real(dp), intent(in) :: covariance(0:)
      type(arfima_simulation_t) :: out
      real(dp), allocatable :: previous(:), current(:)
      real(dp) :: reflection, accumulated
      integer :: i, k, n

      n = size(innovations)
      if (n < 1 .or. ubound(covariance, 1) < n - 1 .or. covariance(0) <= 0.0_dp .or. &
         .not. all(ieee_is_finite(innovations)) .or. &
         .not. all(ieee_is_finite(covariance(:n - 1)))) then
         out%info = 1
         return
      end if
      allocate(out%series(n), out%innovation_variance(n), previous(max(1, n - 1)))
      out%innovation_variance(1) = covariance(0)
      out%series(1) = sqrt(covariance(0))*innovations(1)
      if (n == 1) return
      previous = 0.0_dp
      do k = 1, n - 1
         accumulated = 0.0_dp
         do i = 1, k - 1
            accumulated = accumulated + previous(i)*covariance(k - i)
         end do
         reflection = (covariance(k) - accumulated)/out%innovation_variance(k)
         allocate(current(k))
         do i = 1, k - 1
            current(i) = previous(i) - reflection*previous(k - i)
         end do
         current(k) = reflection
         out%innovation_variance(k + 1) = out%innovation_variance(k)* &
            (1.0_dp - reflection**2)
         if (out%innovation_variance(k + 1) <= epsilon(1.0_dp)*covariance(0)) then
            out%info = 2
            return
         end if
         accumulated = 0.0_dp
         do i = 1, k
            accumulated = accumulated + current(i)*out%series(k + 1 - i)
         end do
         out%series(k + 1) = accumulated + &
            sqrt(out%innovation_variance(k + 1))*innovations(k + 1)
         previous(:k) = current
         deallocate(current)
      end do
   end function arfima_durbin_levinson_simulate

   pure function arfima_simulate_from_innovations(model, innovations, &
      initial_values) result(out)
      ! Simulate, center, and optionally inverse-difference an ARFIMA process.
      type(arfima_model_t), intent(in) :: model
      real(dp), intent(in) :: innovations(:)
      real(dp), intent(in), optional :: initial_values(:)
      type(arfima_simulation_t) :: out
      type(arfima_acvf_t) :: covariance_result
      real(dp), allocatable :: initial(:), integrated(:)
      integer :: initial_count, n

      n = size(innovations)
      if (n < 1 .or. .not. valid_model(model)) then
         out%info = 1
         return
      end if
      covariance_result = arfima_model_acvf(model, n - 1)
      if (covariance_result%info /= 0) then
         out%info = covariance_result%info
         return
      end if
      out = arfima_durbin_levinson_simulate(innovations, &
         covariance_result%covariance)
      if (out%info /= 0) return
      out%covariance = covariance_result%covariance
      out%stationary_series = out%series - sum(out%series)/real(n, dp) + model%mean
      out%series = out%stationary_series

      initial_count = model%difference_order + &
         model%period*model%seasonal_difference_order
      if (initial_count == 0) return
      allocate(initial(initial_count))
      initial = 0.0_dp
      if (present(initial_values)) then
         if (size(initial_values) /= initial_count) then
            out%info = 1
            return
         end if
         initial = initial_values
      end if
      integrated = arfima_integrate(out%stationary_series, initial, &
         model%difference_order, model%seasonal_difference_order, model%period)
      if (size(integrated) /= n + initial_count) then
         out%info = 2
         return
      end if
      out%series = integrated(initial_count + 1:)
   end function arfima_simulate_from_innovations

   function arfima_simulate(model, observations, seed, initial_values) result(out)
      ! Simulate an ARFIMA process using standard-normal innovations.
      type(arfima_model_t), intent(in) :: model
      integer, intent(in) :: observations
      integer, intent(in), optional :: seed
      real(dp), intent(in), optional :: initial_values(:)
      type(arfima_simulation_t) :: out
      real(dp), allocatable :: innovations(:)
      integer :: i

      if (observations < 1) then
         out%info = 1
         return
      end if
      if (present(seed)) call set_random_seed(seed)
      allocate(innovations(observations))
      do i = 1, observations
         innovations(i) = random_standard_normal()
      end do
      if (present(initial_values)) then
         out = arfima_simulate_from_innovations(model, innovations, initial_values)
      else
         out = arfima_simulate_from_innovations(model, innovations)
      end if
   end function arfima_simulate

   pure function arfima_likelihood(series, model, method, estimate_mean) result(out)
      ! Evaluate a profiled exact or conditional Gaussian likelihood.
      real(dp), intent(in) :: series(:)
      type(arfima_model_t), intent(in) :: model
      integer, intent(in), optional :: method
      logical, intent(in), optional :: estimate_mean
      type(arfima_likelihood_t) :: out
      integer :: selected_method
      logical :: profile_mean

      selected_method = arfima_likelihood_exact
      if (present(method)) selected_method = method
      profile_mean = .true.
      if (present(estimate_mean)) profile_mean = estimate_mean
      if (selected_method == arfima_likelihood_exact) then
         out = exact_likelihood(series, model, profile_mean)
      else if (selected_method == arfima_likelihood_css) then
         out = css_likelihood(series, model, profile_mean)
      else
         out%info = 1
      end if
   end function arfima_likelihood

   pure function exact_likelihood(series, model, estimate_mean) result(out)
      ! Evaluate the exact likelihood using the Durbin-Levinson recursion.
      real(dp), intent(in) :: series(:)
      type(arfima_model_t), intent(in) :: model
      logical, intent(in) :: estimate_mean
      type(arfima_likelihood_t) :: out
      type(arfima_model_t) :: unit_model
      type(arfima_acvf_t) :: covariance_result
      real(dp), allocatable :: work(:), base_innovation(:), mean_loading(:)
      real(dp), allocatable :: previous(:), current(:)
      real(dp) :: accumulated, denominator, log_determinant, reflection
      integer :: i, k, n

      out%method = arfima_likelihood_exact
      work = difference_series(series, model%difference_order, &
         model%seasonal_difference_order, model%period)
      n = size(work)
      out%observations = n
      if (n < 2 .or. .not. valid_model(model) .or. .not. all(ieee_is_finite(work))) then
         out%info = 1
         return
      end if
      unit_model = model
      unit_model%innovation_variance = 1.0_dp
      covariance_result = arfima_model_acvf(unit_model, n - 1)
      if (covariance_result%info /= 0) then
         out%info = covariance_result%info
         return
      end if

      allocate(base_innovation(n), mean_loading(n), out%prediction_variance(n))
      allocate(previous(n - 1), current(n - 1))
      previous = 0.0_dp
      current = 0.0_dp
      base_innovation(1) = work(1)
      mean_loading(1) = 1.0_dp
      out%prediction_variance(1) = covariance_result%covariance(0)
      do k = 1, n - 1
         accumulated = 0.0_dp
         do i = 1, k - 1
            accumulated = accumulated + previous(i)* &
               covariance_result%covariance(k - i)
         end do
         reflection = (covariance_result%covariance(k) - accumulated)/ &
            out%prediction_variance(k)
         do i = 1, k - 1
            current(i) = previous(i) - reflection*previous(k - i)
         end do
         current(k) = reflection
         out%prediction_variance(k + 1) = out%prediction_variance(k)* &
            (1.0_dp - reflection**2)
         if (out%prediction_variance(k + 1) <= &
            epsilon(1.0_dp)*covariance_result%covariance(0)) then
            out%info = 2
            return
         end if
         base_innovation(k + 1) = work(k + 1)
         do i = 1, k
            base_innovation(k + 1) = base_innovation(k + 1) - &
               current(i)*work(k + 1 - i)
         end do
         mean_loading(k + 1) = 1.0_dp - sum(current(:k))
         previous(:k) = current(:k)
      end do

      if (estimate_mean) then
         denominator = sum(mean_loading**2/out%prediction_variance)
         if (denominator <= tiny(1.0_dp)) then
            out%info = 2
            return
         end if
         out%mean = sum(base_innovation*mean_loading/out%prediction_variance)/denominator
      else
         out%mean = model%mean
      end if
      out%innovation = base_innovation - out%mean*mean_loading
      out%sum_squares = sum(out%innovation**2/out%prediction_variance)
      out%innovation_variance = out%sum_squares/real(n, dp)
      if (out%innovation_variance <= tiny(1.0_dp)) then
         out%info = 2
         return
      end if
      log_determinant = sum(log(out%prediction_variance))
      out%log_likelihood = -0.5_dp*(real(n, dp)*(log(2.0_dp*acos(-1.0_dp)) + &
         1.0_dp + log(out%innovation_variance)) + log_determinant)
      out%prediction_variance = out%innovation_variance*out%prediction_variance
   end function exact_likelihood

   pure function css_likelihood(series, model, estimate_mean) result(out)
      ! Evaluate a conditional likelihood from truncated impulse weights.
      real(dp), intent(in) :: series(:)
      type(arfima_model_t), intent(in) :: model
      logical, intent(in) :: estimate_mean
      type(arfima_likelihood_t) :: out
      real(dp), allocatable :: work(:), weights(:), base_innovation(:), mean_loading(:)
      real(dp) :: denominator
      integer :: i, lag, n

      out%method = arfima_likelihood_css
      work = difference_series(series, model%difference_order, &
         model%seasonal_difference_order, model%period)
      n = size(work)
      out%observations = n
      if (n < 2 .or. .not. valid_model(model) .or. &
         (model%long_memory_type /= arfima_long_memory_fdwn .and. &
         model%long_memory_type /= arfima_long_memory_none) .or. &
         (model%period > 0 .and. &
         model%seasonal_long_memory_type /= arfima_long_memory_fdwn .and. &
         model%seasonal_long_memory_type /= arfima_long_memory_none)) then
         out%info = 1
         return
      end if
      if (model%period > 0) then
         weights = arfima_seasonal_psi_weights(model%ar, model%theta, &
            model%long_memory_parameter, model%seasonal_ar, model%seasonal_theta, &
            model%seasonal_long_memory_parameter, model%period, n - 1)
      else
         weights = arfima_psi_weights(model%ar, model%theta, &
            model%long_memory_parameter, n - 1)
      end if
      allocate(base_innovation(n), mean_loading(n), out%prediction_variance(n))
      base_innovation = 0.0_dp
      mean_loading = 0.0_dp
      do i = 1, n
         base_innovation(i) = work(i)
         mean_loading(i) = 1.0_dp
         do lag = 1, i - 1
            base_innovation(i) = base_innovation(i) - weights(lag + 1)* &
               base_innovation(i - lag)
            mean_loading(i) = mean_loading(i) - weights(lag + 1)* &
               mean_loading(i - lag)
         end do
      end do
      if (estimate_mean) then
         denominator = sum(mean_loading**2)
         if (denominator <= tiny(1.0_dp)) then
            out%info = 2
            return
         end if
         out%mean = sum(base_innovation*mean_loading)/denominator
      else
         out%mean = model%mean
      end if
      out%innovation = base_innovation - out%mean*mean_loading
      out%sum_squares = sum(out%innovation**2)
      out%innovation_variance = out%sum_squares/real(n, dp)
      if (out%innovation_variance <= tiny(1.0_dp)) then
         out%info = 2
         return
      end if
      out%prediction_variance = out%innovation_variance
      out%log_likelihood = -0.5_dp*real(n, dp)*(log(2.0_dp*acos(-1.0_dp)) + &
         1.0_dp + log(out%innovation_variance))
   end function css_likelihood

   pure function arfima_fit(series, initial_model, method, estimate_mean, &
      max_iterations, tolerance) result(out)
      ! Estimate ARFIMA parameters by finite-difference BFGS.
      real(dp), intent(in) :: series(:)
      type(arfima_model_t), intent(in) :: initial_model
      integer, intent(in), optional :: method, max_iterations
      logical, intent(in), optional :: estimate_mean
      real(dp), intent(in), optional :: tolerance
      type(arfima_fit_t) :: out
      type(optimization_result_t) :: optimization
      real(dp), allocatable :: initial(:), hessian(:, :)
      integer :: hessian_info, information_parameters, selected_method
      logical :: profile_mean

      selected_method = arfima_likelihood_exact
      if (present(method)) selected_method = method
      profile_mean = .true.
      if (present(estimate_mean)) profile_mean = estimate_mean
      if (.not. valid_model(initial_model) .or. size(series) < 3 .or. &
         (selected_method /= arfima_likelihood_exact .and. &
         selected_method /= arfima_likelihood_css)) then
         out%info = 1
         return
      end if
      initial = encode_model_parameters(initial_model)
      out%parameter_count = size(initial)
      out%estimated_parameter_count = size(initial)
      allocate(out%fixed_parameters(size(initial)))
      out%fixed_parameters = .false.
      if (size(initial) == 0) then
         out%model = initial_model
         out%likelihood = arfima_likelihood(series, out%model, selected_method, profile_mean)
         if (out%likelihood%info /= 0) then
            out%info = out%likelihood%info
            return
         end if
         allocate(out%unconstrained_parameters(0), out%parameter_covariance(0, 0))
         allocate(out%standard_error(0))
         out%converged = .true.
      else
         optimization = bfgs_minimize_fd(objective, initial, &
            max_iterations=max_iterations, gradient_tolerance=tolerance)
         out%iterations = optimization%iterations
         out%converged = optimization%converged
         out%unconstrained_parameters = optimization%parameters
         out%model = decode_model_parameters(initial_model, optimization%parameters)
         out%likelihood = arfima_likelihood(series, out%model, selected_method, profile_mean)
         if (out%likelihood%info /= 0) then
            out%info = out%likelihood%info
            return
         end if
         hessian = finite_difference_hessian(objective, optimization%parameters)
         call invert_matrix(hessian, out%parameter_covariance, hessian_info)
         allocate(out%standard_error(size(initial)))
         if (hessian_info == 0 .and. all([(out%parameter_covariance( &
            information_parameters, information_parameters) > 0.0_dp, &
            information_parameters=1, size(initial))])) then
            do information_parameters = 1, size(initial)
               out%standard_error(information_parameters) = sqrt(out%parameter_covariance( &
                  information_parameters, information_parameters))
            end do
         else
            out%parameter_covariance = 0.0_dp
            out%standard_error = 0.0_dp
         end if
         if (.not. optimization%converged) out%info = optimization%info
      end if

      out%model%mean = out%likelihood%mean
      out%model%innovation_variance = out%likelihood%innovation_variance
      information_parameters = out%parameter_count + merge(1, 0, profile_mean) + 1
      out%aic = -2.0_dp*out%likelihood%log_likelihood + &
         2.0_dp*real(information_parameters, dp)
      out%bic = -2.0_dp*out%likelihood%log_likelihood + &
         log(real(out%likelihood%observations, dp))*real(information_parameters, dp)
      if (out%likelihood%observations > information_parameters + 1) then
         out%aicc = out%aic + 2.0_dp*real(information_parameters* &
            (information_parameters + 1), dp)/ &
            real(out%likelihood%observations - information_parameters - 1, dp)
      end if

   contains

      pure function objective(parameters) result(value)
         ! Return the profiled negative log-likelihood for BFGS.
         real(dp), intent(in) :: parameters(:)
         real(dp) :: value
         type(arfima_model_t) :: candidate
         type(arfima_likelihood_t) :: likelihood

         candidate = decode_model_parameters(initial_model, parameters)
         likelihood = arfima_likelihood(series, candidate, selected_method, profile_mean)
         if (likelihood%info /= 0 .or. .not. ieee_is_finite(likelihood%log_likelihood)) then
            value = huge(1.0_dp)/100.0_dp
         else
            value = -likelihood%log_likelihood
         end if
      end function objective

   end function arfima_fit

   pure function arfima_fit_fixed(series, initial_model, fixed_parameters, method, &
      estimate_mean, max_iterations, tolerance) result(out)
      ! Estimate free ARFIMA coordinates while retaining selected initial values.
      real(dp), intent(in) :: series(:)
      type(arfima_model_t), intent(in) :: initial_model
      logical, intent(in) :: fixed_parameters(:)
      integer, intent(in), optional :: method, max_iterations
      logical, intent(in), optional :: estimate_mean
      real(dp), intent(in), optional :: tolerance
      type(arfima_fit_t) :: out
      type(optimization_result_t) :: optimization
      real(dp), allocatable :: initial(:), free_initial(:), full(:)
      real(dp), allocatable :: hessian(:, :), free_covariance(:, :)
      integer, allocatable :: free_index(:)
      integer :: hessian_info, i, j, information_parameters, selected_method
      logical :: profile_mean, valid_covariance

      selected_method = arfima_likelihood_exact
      if (present(method)) selected_method = method
      profile_mean = .true.
      if (present(estimate_mean)) profile_mean = estimate_mean
      initial = encode_model_parameters(initial_model)
      out%parameter_count = size(initial)
      out%estimated_parameter_count = count(.not. fixed_parameters)
      if (.not. valid_model(initial_model) .or. size(series) < 3 .or. &
         size(fixed_parameters) /= size(initial) .or. &
         (selected_method /= arfima_likelihood_exact .and. &
         selected_method /= arfima_likelihood_css)) then
         out%info = 1
         return
      end if
      out%fixed_parameters = fixed_parameters
      free_index = pack([(i, i=1, size(initial))], .not. fixed_parameters)
      free_initial = initial(free_index)
      full = initial
      allocate(out%parameter_covariance(size(initial), size(initial)))
      allocate(out%standard_error(size(initial)))
      out%parameter_covariance = 0.0_dp
      out%standard_error = 0.0_dp
      if (size(free_initial) == 0) then
         out%unconstrained_parameters = full
         out%model = initial_model
         out%likelihood = arfima_likelihood(series, out%model, selected_method, profile_mean)
         out%converged = out%likelihood%info == 0
      else
         optimization = bfgs_minimize_fd(objective, free_initial, &
            max_iterations=max_iterations, gradient_tolerance=tolerance)
         out%iterations = optimization%iterations
         out%converged = optimization%converged
         full(free_index) = optimization%parameters
         out%unconstrained_parameters = full
         out%model = decode_model_parameters(initial_model, full)
         call impose_fixed_model(initial_model, fixed_parameters, out%model)
         out%likelihood = arfima_likelihood(series, out%model, selected_method, profile_mean)
         if (out%likelihood%info /= 0) then
            out%info = out%likelihood%info
            return
         end if
         hessian = finite_difference_hessian(objective, optimization%parameters)
         call invert_matrix(hessian, free_covariance, hessian_info)
         valid_covariance = hessian_info == 0
         if (valid_covariance) then
            do i = 1, size(free_initial)
               if (free_covariance(i, i) <= 0.0_dp) valid_covariance = .false.
            end do
         end if
         if (valid_covariance) then
            do j = 1, size(free_initial)
               do i = 1, size(free_initial)
                  out%parameter_covariance(free_index(i), free_index(j)) = &
                     free_covariance(i, j)
               end do
               out%standard_error(free_index(j)) = sqrt(free_covariance(j, j))
            end do
         end if
         if (.not. optimization%converged) out%info = optimization%info
      end if
      if (out%likelihood%info /= 0) then
         out%info = out%likelihood%info
         return
      end if
      out%model%mean = out%likelihood%mean
      out%model%innovation_variance = out%likelihood%innovation_variance
      information_parameters = out%estimated_parameter_count + &
         merge(1, 0, profile_mean) + 1
      out%aic = -2.0_dp*out%likelihood%log_likelihood + &
         2.0_dp*real(information_parameters, dp)
      out%bic = -2.0_dp*out%likelihood%log_likelihood + &
         log(real(out%likelihood%observations, dp))*real(information_parameters, dp)
      if (out%likelihood%observations > information_parameters + 1) then
         out%aicc = out%aic + 2.0_dp*real(information_parameters* &
            (information_parameters + 1), dp)/ &
            real(out%likelihood%observations - information_parameters - 1, dp)
      end if

   contains

      pure function objective(free_parameters) result(value)
         ! Return the negative likelihood after inserting fixed coordinates.
         real(dp), intent(in) :: free_parameters(:)
         real(dp) :: value
         real(dp) :: candidate_parameters(size(initial))
         type(arfima_model_t) :: candidate
         type(arfima_likelihood_t) :: likelihood

         candidate_parameters = initial
         candidate_parameters(free_index) = free_parameters
         candidate = decode_model_parameters(initial_model, candidate_parameters)
         call impose_fixed_model(initial_model, fixed_parameters, candidate)
         likelihood = arfima_likelihood(series, candidate, selected_method, profile_mean)
         if (likelihood%info /= 0 .or. .not. ieee_is_finite(likelihood%log_likelihood)) then
            value = huge(1.0_dp)/100.0_dp
         else
            value = -likelihood%log_likelihood
         end if
      end function objective

   end function arfima_fit_fixed

   pure function arfima_fit_modes(series, initial_model, starts, method, &
      estimate_mean, max_iterations, tolerance) result(out)
      ! Fit independent local modes from caller-supplied optimizer coordinates.
      real(dp), intent(in) :: series(:), starts(:, :)
      type(arfima_model_t), intent(in) :: initial_model
      integer, intent(in), optional :: method, max_iterations
      logical, intent(in), optional :: estimate_mean
      real(dp), intent(in), optional :: tolerance
      type(arfima_multifit_t) :: out
      type(arfima_fit_t), allocatable :: candidates(:)
      type(arfima_model_t) :: start_model
      logical, allocatable :: keep(:)
      integer :: i, retained_count, selected_method, iteration_limit
      logical :: profile_mean
      real(dp) :: selected_tolerance

      selected_method = arfima_likelihood_exact
      if (present(method)) selected_method = method
      iteration_limit = 200
      if (present(max_iterations)) iteration_limit = max_iterations
      selected_tolerance = 1.0e-6_dp
      if (present(tolerance)) selected_tolerance = tolerance
      profile_mean = .true.
      if (present(estimate_mean)) profile_mean = estimate_mean
      out%attempted = size(starts, 2)
      out%starting_parameters = starts
      if (size(starts, 1) /= size(encode_model_parameters(initial_model)) .or. &
         size(starts, 2) < 1 .or. iteration_limit < 1 .or. &
         selected_tolerance <= 0.0_dp) then
         out%info = 1
         call allocate_empty_multifit(out)
         return
      end if
      allocate(candidates(size(starts, 2)), keep(size(starts, 2)))
      keep = .false.
      do i = 1, size(starts, 2)
         start_model = decode_model_parameters(initial_model, starts(:, i))
         candidates(i) = arfima_fit(series, start_model, selected_method, &
            profile_mean, iteration_limit, selected_tolerance)
         keep(i) = candidates(i)%converged .and. &
            candidates(i)%likelihood%info == 0 .and. &
            ieee_is_finite(candidates(i)%likelihood%log_likelihood)
      end do
      retained_count = count(keep)
      out%converged = retained_count
      allocate(out%modes(retained_count), out%starting_index(retained_count))
      allocate(out%long_memory_type(retained_count), out%weight(retained_count))
      retained_count = 0
      do i = 1, size(candidates)
         if (.not. keep(i)) cycle
         retained_count = retained_count + 1
         out%modes(retained_count) = candidates(i)
         out%starting_index(retained_count) = i
         out%long_memory_type(retained_count) = &
            candidates(i)%model%long_memory_type
      end do
      out%weight = mode_weights(out%modes)
      if (size(out%modes) == 0) out%info = 2
   end function arfima_fit_modes

   pure function arfima_fit_modes_fixed(series, initial_model, fixed_parameters, &
      starts, method, estimate_mean, max_iterations, tolerance) result(out)
      ! Fit local modes while retaining selected coordinates at initial values.
      real(dp), intent(in) :: series(:), starts(:, :)
      type(arfima_model_t), intent(in) :: initial_model
      logical, intent(in) :: fixed_parameters(:)
      integer, intent(in), optional :: method, max_iterations
      logical, intent(in), optional :: estimate_mean
      real(dp), intent(in), optional :: tolerance
      type(arfima_multifit_t) :: out
      type(arfima_fit_t), allocatable :: candidates(:)
      type(arfima_model_t) :: start_model
      real(dp), allocatable :: base(:), effective_starts(:, :)
      logical, allocatable :: keep(:)
      integer :: i, retained_count, selected_method, iteration_limit
      logical :: profile_mean
      real(dp) :: selected_tolerance

      base = encode_model_parameters(initial_model)
      out%attempted = size(starts, 2)
      if (size(fixed_parameters) /= size(base) .or. &
         size(starts, 1) /= size(base) .or. size(starts, 2) < 1) then
         out%info = 1
         call allocate_empty_multifit(out)
         return
      end if
      effective_starts = starts
      do i = 1, size(starts, 2)
         where (fixed_parameters) effective_starts(:, i) = base
      end do
      out%starting_parameters = effective_starts
      selected_method = arfima_likelihood_exact
      if (present(method)) selected_method = method
      iteration_limit = 200
      if (present(max_iterations)) iteration_limit = max_iterations
      selected_tolerance = 1.0e-6_dp
      if (present(tolerance)) selected_tolerance = tolerance
      profile_mean = .true.
      if (present(estimate_mean)) profile_mean = estimate_mean
      allocate(candidates(size(starts, 2)), keep(size(starts, 2)))
      keep = .false.
      do i = 1, size(starts, 2)
         start_model = decode_model_parameters(initial_model, effective_starts(:, i))
         call impose_fixed_model(initial_model, fixed_parameters, start_model)
         candidates(i) = arfima_fit_fixed(series, start_model, fixed_parameters, &
            selected_method, profile_mean, iteration_limit, selected_tolerance)
         keep(i) = candidates(i)%converged .and. &
            candidates(i)%likelihood%info == 0 .and. &
            ieee_is_finite(candidates(i)%likelihood%log_likelihood)
      end do
      retained_count = count(keep)
      out%converged = retained_count
      allocate(out%modes(retained_count), out%starting_index(retained_count))
      allocate(out%long_memory_type(retained_count), out%weight(retained_count))
      retained_count = 0
      do i = 1, size(candidates)
         if (.not. keep(i)) cycle
         retained_count = retained_count + 1
         out%modes(retained_count) = candidates(i)
         out%starting_index(retained_count) = i
         out%long_memory_type(retained_count) = &
            candidates(i)%model%long_memory_type
      end do
      out%weight = mode_weights(out%modes)
      if (size(out%modes) == 0) out%info = 2
   end function arfima_fit_modes_fixed

   function arfima_multistart_fit_fixed(series, initial_model, fixed_parameters, &
      number_starts, seed, method, estimate_mean, max_iterations, tolerance) result(out)
      ! Fit fixed-parameter modes from seeded random stable starts.
      real(dp), intent(in) :: series(:)
      type(arfima_model_t), intent(in) :: initial_model
      logical, intent(in) :: fixed_parameters(:)
      integer, intent(in) :: number_starts
      integer, intent(in), optional :: seed, method, max_iterations
      logical, intent(in), optional :: estimate_mean
      real(dp), intent(in), optional :: tolerance
      type(arfima_multifit_t) :: out
      real(dp), allocatable :: starts(:, :)
      integer :: selected_method, iteration_limit
      logical :: profile_mean
      real(dp) :: selected_tolerance

      if (present(seed)) call set_random_seed(seed)
      starts = random_model_starts(initial_model, number_starts)
      selected_method = arfima_likelihood_exact
      if (present(method)) selected_method = method
      iteration_limit = 200
      if (present(max_iterations)) iteration_limit = max_iterations
      selected_tolerance = 1.0e-6_dp
      if (present(tolerance)) selected_tolerance = tolerance
      profile_mean = .true.
      if (present(estimate_mean)) profile_mean = estimate_mean
      out = arfima_fit_modes_fixed(series, initial_model, fixed_parameters, starts, &
         selected_method, profile_mean, iteration_limit, selected_tolerance)
   end function arfima_multistart_fit_fixed

   function arfima_multistart_fit(series, initial_model, number_starts, seed, &
      method, estimate_mean, max_iterations, tolerance) result(out)
      ! Fit local modes from seeded random stable PACF and memory starts.
      real(dp), intent(in) :: series(:)
      type(arfima_model_t), intent(in) :: initial_model
      integer, intent(in) :: number_starts
      integer, intent(in), optional :: seed, method, max_iterations
      logical, intent(in), optional :: estimate_mean
      real(dp), intent(in), optional :: tolerance
      type(arfima_multifit_t) :: out
      real(dp), allocatable :: starts(:, :)
      integer :: selected_method, iteration_limit
      logical :: profile_mean
      real(dp) :: selected_tolerance

      if (present(seed)) call set_random_seed(seed)
      starts = random_model_starts(initial_model, number_starts)
      selected_method = arfima_likelihood_exact
      if (present(method)) selected_method = method
      iteration_limit = 200
      if (present(max_iterations)) iteration_limit = max_iterations
      selected_tolerance = 1.0e-6_dp
      if (present(tolerance)) selected_tolerance = tolerance
      profile_mean = .true.
      if (present(estimate_mean)) profile_mean = estimate_mean
      out = arfima_fit_modes(series, initial_model, starts, selected_method, &
         profile_mean, iteration_limit, selected_tolerance)
   end function arfima_multistart_fit

   function arfima_select_long_memory(series, initial_model, starts_per_family, &
      seed, method, estimate_mean, max_iterations, tolerance, weed_tolerance) result(out)
      ! Compare none, FDWN, FGN, and hyperbolic-decay likelihood modes.
      real(dp), intent(in) :: series(:)
      type(arfima_model_t), intent(in) :: initial_model
      integer, intent(in) :: starts_per_family
      integer, intent(in), optional :: seed, method, max_iterations
      logical, intent(in), optional :: estimate_mean
      real(dp), intent(in), optional :: tolerance, weed_tolerance
      type(arfima_multifit_t) :: out
      type(arfima_multifit_t) :: groups(4)
      type(arfima_model_t) :: candidate
      integer :: family, selected_method, iteration_limit
      logical :: profile_mean
      real(dp) :: selected_tolerance, duplicate_tolerance

      if (present(seed)) call set_random_seed(seed)
      selected_method = arfima_likelihood_exact
      if (present(method)) selected_method = method
      iteration_limit = 200
      if (present(max_iterations)) iteration_limit = max_iterations
      selected_tolerance = 1.0e-6_dp
      if (present(tolerance)) selected_tolerance = tolerance
      duplicate_tolerance = 0.025_dp
      if (present(weed_tolerance)) duplicate_tolerance = weed_tolerance
      profile_mean = .true.
      if (present(estimate_mean)) profile_mean = estimate_mean
      do family = arfima_long_memory_none, arfima_long_memory_hd
         candidate = initial_model
         candidate%long_memory_type = family
         select case (family)
         case (arfima_long_memory_none)
            candidate%long_memory_parameter = 0.0_dp
         case (arfima_long_memory_fdwn)
            candidate%long_memory_parameter = 0.1_dp
         case (arfima_long_memory_fgn)
            candidate%long_memory_parameter = 0.7_dp
         case (arfima_long_memory_hd)
            candidate%long_memory_parameter = 1.0_dp
         end select
         groups(family + 1) = arfima_fit_modes(series, candidate, &
            random_model_starts(candidate, starts_per_family), selected_method, &
            profile_mean, iteration_limit, selected_tolerance)
      end do
      out = combine_multifits(groups)
      if (size(out%modes) > 0) then
         out = arfima_weed_modes(out, arfima_weed_both, duplicate_tolerance, .true.)
      end if
   end function arfima_select_long_memory

   pure real(dp) function arfima_mode_distance(first, second, p, transformed) &
      result(distance)
      ! Return a p-norm distance in coefficient or PACF parameter space.
      type(arfima_fit_t), intent(in) :: first, second
      real(dp), intent(in), optional :: p
      logical, intent(in), optional :: transformed
      real(dp), allocatable :: first_parameters(:), second_parameters(:)
      real(dp) :: norm_order
      logical :: use_transformed

      norm_order = 2.0_dp
      if (present(p)) norm_order = p
      use_transformed = .false.
      if (present(transformed)) use_transformed = transformed
      first_parameters = mode_parameters(first%model, use_transformed)
      second_parameters = mode_parameters(second%model, use_transformed)
      if (norm_order <= 0.0_dp .or. &
         .not. compatible_mode_models(first%model, second%model) .or. &
         size(first_parameters) /= size(second_parameters)) then
         distance = huge(1.0_dp)
      else if (size(first_parameters) == 0) then
         distance = 0.0_dp
      else
         distance = sum(abs(first_parameters - second_parameters)**norm_order)** &
            (1.0_dp/norm_order)
      end if
   end function arfima_mode_distance

   pure function arfima_mode_distances(multifit, p, transformed) result(distances)
      ! Return the symmetric pairwise distance matrix for retained modes.
      type(arfima_multifit_t), intent(in) :: multifit
      real(dp), intent(in), optional :: p
      logical, intent(in), optional :: transformed
      real(dp), allocatable :: distances(:, :)
      real(dp) :: norm_order
      logical :: use_transformed
      integer :: i, j

      norm_order = 2.0_dp
      if (present(p)) norm_order = p
      use_transformed = .false.
      if (present(transformed)) use_transformed = transformed
      allocate(distances(size(multifit%modes), size(multifit%modes)))
      do j = 1, size(multifit%modes)
         do i = 1, j
            distances(i, j) = arfima_mode_distance(multifit%modes(i), &
               multifit%modes(j), norm_order, use_transformed)
            distances(j, i) = distances(i, j)
         end do
      end do
   end function arfima_mode_distances

   pure function arfima_weed_modes(multifit, space, tolerance, adaptive, p, walls, &
      wall_tolerance) result(out)
      ! Remove nearby modes and secondary modes sharing a parameter boundary.
      type(arfima_multifit_t), intent(in) :: multifit
      integer, intent(in), optional :: space
      real(dp), intent(in), optional :: tolerance, p, wall_tolerance
      logical, intent(in), optional :: adaptive, walls
      type(arfima_multifit_t) :: out
      integer, allocatable :: order(:), retained(:)
      integer :: selected_space, i, j, kept, dimension
      real(dp) :: radius, norm_order, selected_wall_tolerance
      logical :: adapt, distinct, use_walls

      selected_space = arfima_weed_operator
      if (present(space)) selected_space = space
      radius = 0.025_dp
      if (present(tolerance)) radius = tolerance
      norm_order = 2.0_dp
      if (present(p)) norm_order = p
      adapt = .true.
      if (present(adaptive)) adapt = adaptive
      use_walls = .false.
      if (present(walls)) use_walls = walls
      selected_wall_tolerance = 0.01_dp
      if (present(wall_tolerance)) selected_wall_tolerance = wall_tolerance
      if (selected_space < arfima_weed_none .or. selected_space > arfima_weed_both .or. &
         radius < 0.0_dp .or. norm_order <= 0.0_dp .or. &
         selected_wall_tolerance <= 0.0_dp) then
         out = multifit
         out%info = 1
         return
      end if
      if (size(multifit%modes) == 0 .or. &
         (selected_space == arfima_weed_none .and. .not. use_walls)) then
         out = multifit
         out%weeded = selected_space /= arfima_weed_none
         return
      end if
      dimension = size(mode_parameters(multifit%modes(1)%model, .false.))
      if (adapt) radius = (1.0_dp + radius)**dimension - 1.0_dp
      order = likelihood_order(multifit%modes)
      allocate(retained(size(order)))
      kept = 0
      do i = 1, size(order)
         distinct = .true.
         do j = 1, kept
            if (selected_space == arfima_weed_operator .or. &
               selected_space == arfima_weed_both) then
               distinct = arfima_mode_distance(multifit%modes(order(i)), &
                  multifit%modes(retained(j)), norm_order, .false.) > radius
            end if
            if (distinct .and. (selected_space == arfima_weed_pacf .or. &
               selected_space == arfima_weed_both)) then
               distinct = arfima_mode_distance(multifit%modes(order(i)), &
                  multifit%modes(retained(j)), norm_order, .true.) > radius
            end if
            if (distinct .and. use_walls) then
               distinct = .not. modes_share_wall(multifit%modes(order(i)), &
                  multifit%modes(retained(j)), selected_wall_tolerance)
            end if
            if (.not. distinct) exit
         end do
         if (.not. distinct) cycle
         kept = kept + 1
         retained(kept) = order(i)
      end do
      out = subset_multifit(multifit, retained(:kept))
      out%weeded = .true.
   end function arfima_weed_modes

   pure function arfima_best_modes(multifit, number) result(out)
      ! Retain the requested number of highest-likelihood modes.
      type(arfima_multifit_t), intent(in) :: multifit
      integer, intent(in) :: number
      type(arfima_multifit_t) :: out
      integer, allocatable :: order(:)

      if (number < 0 .or. number > size(multifit%modes)) then
         out = multifit
         out%info = 1
         return
      end if
      order = likelihood_order(multifit%modes)
      out = subset_multifit(multifit, order(:number))
   end function arfima_best_modes

   pure function arfima_remove_mode(multifit, index) result(out)
      ! Return a multimode fit with one selected mode removed.
      type(arfima_multifit_t), intent(in) :: multifit
      integer, intent(in) :: index
      type(arfima_multifit_t) :: out
      integer, allocatable :: retained(:)
      integer :: i, offset

      if (index < 1 .or. index > size(multifit%modes)) then
         out = multifit
         out%info = 1
         return
      end if
      allocate(retained(size(multifit%modes) - 1))
      offset = 0
      do i = 1, size(multifit%modes)
         if (i == index) cycle
         offset = offset + 1
         retained(offset) = i
      end do
      out = subset_multifit(multifit, retained)
   end function arfima_remove_mode

   pure function arfima_mode_forecast(series, multifit, index, horizon, level) result(out)
      ! Forecast from one selected local likelihood mode.
      real(dp), intent(in) :: series(:)
      type(arfima_multifit_t), intent(in) :: multifit
      integer, intent(in) :: index, horizon
      real(dp), intent(in), optional :: level
      type(arfima_forecast_t) :: out

      if (index < 1 .or. index > size(multifit%modes)) then
         out%info = 1
      else if (present(level)) then
         out = arfima_forecast(series, multifit%modes(index), horizon, level)
      else
         out = arfima_forecast(series, multifit%modes(index), horizon)
      end if
   end function arfima_mode_forecast

   pure function arfima_mode_information(multifit, index, exact, resolution) result(out)
      ! Compute Fisher information for one selected local mode.
      type(arfima_multifit_t), intent(in) :: multifit
      integer, intent(in) :: index
      logical, intent(in), optional :: exact
      integer, intent(in), optional :: resolution
      type(arfima_information_t) :: out

      if (index < 1 .or. index > size(multifit%modes)) then
         out%info = 1
      else if (present(exact) .and. present(resolution)) then
         out = arfima_fit_information(multifit%modes(index), exact, resolution)
      else if (present(exact)) then
         out = arfima_fit_information(multifit%modes(index), exact=exact)
      else if (present(resolution)) then
         out = arfima_fit_information(multifit%modes(index), resolution=resolution)
      else
         out = arfima_fit_information(multifit%modes(index))
      end if
   end function arfima_mode_information

   pure function arfima_mode_identifiability(multifit, index) result(out)
      ! Diagnose stability and common factors for one selected local mode.
      type(arfima_multifit_t), intent(in) :: multifit
      integer, intent(in) :: index
      type(arfima_identifiability_t) :: out

      if (index < 1 .or. index > size(multifit%modes)) then
         out%info = 1
      else
         out = arfima_identifiability(multifit%modes(index)%model)
      end if
   end function arfima_mode_identifiability

   pure function arfima_regression_fit(series, initial_model, regressors, &
      initial_coefficients, initial_transfers, transfer_regressors, method, &
      estimate_mean, max_iterations, tolerance) result(out)
      ! Jointly estimate ARFIMA, static regression, and transfer parameters.
      real(dp), intent(in) :: series(:), regressors(:, :), initial_coefficients(:)
      type(arfima_model_t), intent(in) :: initial_model
      type(arfima_transfer_t), intent(in) :: initial_transfers(:)
      real(dp), intent(in) :: transfer_regressors(:, :)
      integer, intent(in), optional :: method, max_iterations
      logical, intent(in), optional :: estimate_mean
      real(dp), intent(in), optional :: tolerance
      type(arfima_regression_fit_t) :: out
      type(optimization_result_t) :: optimization
      real(dp), allocatable :: initial(:), hessian(:, :)
      integer :: hessian_info, i, information_parameters, selected_method
      logical :: profile_mean, valid_covariance

      selected_method = arfima_likelihood_exact
      if (present(method)) selected_method = method
      profile_mean = .true.
      if (present(estimate_mean)) profile_mean = estimate_mean
      if (.not. valid_model(initial_model) .or. &
         .not. valid_regression_inputs(series, regressors, initial_coefficients, &
         initial_transfers, transfer_regressors)) then
         out%info = 1
         return
      end if
      initial = encode_regression_parameters(initial_model, initial_coefficients, &
         initial_transfers)
      out%parameter_count = size(initial)
      out%estimated_parameter_count = size(initial)
      allocate(out%fixed_parameters(size(initial)))
      out%fixed_parameters = .false.
      if (size(initial) == 0) then
         out%model = initial_model
         out%coefficients = initial_coefficients
         out%transfers = initial_transfers
         out%likelihood = arfima_regression_likelihood(series, out%model, regressors, &
            out%coefficients, out%transfers, transfer_regressors, selected_method, &
            profile_mean)
         allocate(out%unconstrained_parameters(0), out%parameter_covariance(0, 0))
         allocate(out%standard_error(0))
         out%converged = out%likelihood%info == 0
      else
         optimization = bfgs_minimize_fd(objective, initial, &
            max_iterations=max_iterations, gradient_tolerance=tolerance)
         out%iterations = optimization%iterations
         out%converged = optimization%converged
         out%unconstrained_parameters = optimization%parameters
         call decode_regression_parameters(initial_model, initial_coefficients, &
            initial_transfers, optimization%parameters, out%model, out%coefficients, &
            out%transfers)
         out%likelihood = arfima_regression_likelihood(series, out%model, regressors, &
            out%coefficients, out%transfers, transfer_regressors, selected_method, &
            profile_mean)
         if (out%likelihood%info /= 0) then
            out%info = out%likelihood%info
            return
         end if
         hessian = finite_difference_hessian(objective, optimization%parameters)
         call invert_matrix(hessian, out%parameter_covariance, hessian_info)
         allocate(out%standard_error(size(initial)))
         valid_covariance = hessian_info == 0
         if (valid_covariance) then
            do i = 1, size(initial)
               if (out%parameter_covariance(i, i) <= 0.0_dp) valid_covariance = .false.
            end do
         end if
         if (valid_covariance) then
            do i = 1, size(initial)
               out%standard_error(i) = sqrt(out%parameter_covariance(i, i))
            end do
         else
            out%parameter_covariance = 0.0_dp
            out%standard_error = 0.0_dp
         end if
         if (.not. optimization%converged) out%info = optimization%info
      end if
      if (out%likelihood%info /= 0) then
         out%info = out%likelihood%info
         return
      end if
      out%model%mean = out%likelihood%mean
      out%model%innovation_variance = out%likelihood%innovation_variance
      information_parameters = out%parameter_count + merge(1, 0, profile_mean) + 1
      out%aic = -2.0_dp*out%likelihood%log_likelihood + &
         2.0_dp*real(information_parameters, dp)
      out%bic = -2.0_dp*out%likelihood%log_likelihood + &
         log(real(out%likelihood%observations, dp))*real(information_parameters, dp)
      if (out%likelihood%observations > information_parameters + 1) then
         out%aicc = out%aic + 2.0_dp*real(information_parameters* &
            (information_parameters + 1), dp)/ &
            real(out%likelihood%observations - information_parameters - 1, dp)
      end if

   contains

      pure function objective(parameters) result(value)
         ! Return the joint profiled negative log-likelihood.
         real(dp), intent(in) :: parameters(:)
         real(dp) :: value
         type(arfima_model_t) :: candidate_model
         type(arfima_transfer_t), allocatable :: candidate_transfers(:)
         real(dp), allocatable :: candidate_coefficients(:)
         type(arfima_likelihood_t) :: candidate_likelihood

         call decode_regression_parameters(initial_model, initial_coefficients, &
            initial_transfers, parameters, candidate_model, candidate_coefficients, &
            candidate_transfers)
         candidate_likelihood = arfima_regression_likelihood(series, candidate_model, &
            regressors, candidate_coefficients, candidate_transfers, transfer_regressors, &
            selected_method, profile_mean)
         if (candidate_likelihood%info /= 0 .or. &
            .not. ieee_is_finite(candidate_likelihood%log_likelihood)) then
            value = huge(1.0_dp)/100.0_dp
         else
            value = -candidate_likelihood%log_likelihood
         end if
      end function objective

   end function arfima_regression_fit

   pure function arfima_regression_fit_fixed(series, initial_model, regressors, &
      initial_coefficients, initial_transfers, transfer_regressors, fixed_parameters, &
      method, estimate_mean, max_iterations, tolerance) result(out)
      ! Estimate a joint regression model with selected coordinates held fixed.
      real(dp), intent(in) :: series(:), regressors(:, :), initial_coefficients(:)
      type(arfima_model_t), intent(in) :: initial_model
      type(arfima_transfer_t), intent(in) :: initial_transfers(:)
      real(dp), intent(in) :: transfer_regressors(:, :)
      logical, intent(in) :: fixed_parameters(:)
      integer, intent(in), optional :: method, max_iterations
      logical, intent(in), optional :: estimate_mean
      real(dp), intent(in), optional :: tolerance
      type(arfima_regression_fit_t) :: out
      type(optimization_result_t) :: optimization
      real(dp), allocatable :: initial(:), free_initial(:), full(:)
      real(dp), allocatable :: hessian(:, :), free_covariance(:, :)
      integer, allocatable :: free_index(:)
      integer :: hessian_info, i, j, information_parameters, selected_method
      logical :: profile_mean, valid_covariance

      selected_method = arfima_likelihood_exact
      if (present(method)) selected_method = method
      profile_mean = .true.
      if (present(estimate_mean)) profile_mean = estimate_mean
      initial = encode_regression_parameters(initial_model, initial_coefficients, &
         initial_transfers)
      out%parameter_count = size(initial)
      out%estimated_parameter_count = count(.not. fixed_parameters)
      if (.not. valid_model(initial_model) .or. &
         size(fixed_parameters) /= size(initial) .or. &
         .not. valid_regression_inputs(series, regressors, initial_coefficients, &
         initial_transfers, transfer_regressors)) then
         out%info = 1
         return
      end if
      out%fixed_parameters = fixed_parameters
      free_index = pack([(i, i=1, size(initial))], .not. fixed_parameters)
      free_initial = initial(free_index)
      full = initial
      allocate(out%parameter_covariance(size(initial), size(initial)))
      allocate(out%standard_error(size(initial)))
      out%parameter_covariance = 0.0_dp
      out%standard_error = 0.0_dp
      if (size(free_initial) == 0) then
         out%unconstrained_parameters = full
         out%model = initial_model
         out%coefficients = initial_coefficients
         out%transfers = initial_transfers
         out%likelihood = arfima_regression_likelihood(series, out%model, regressors, &
            out%coefficients, out%transfers, transfer_regressors, selected_method, &
            profile_mean)
         out%converged = out%likelihood%info == 0
      else
         optimization = bfgs_minimize_fd(objective, free_initial, &
            max_iterations=max_iterations, gradient_tolerance=tolerance)
         out%iterations = optimization%iterations
         out%converged = optimization%converged
         full(free_index) = optimization%parameters
         out%unconstrained_parameters = full
         call decode_regression_parameters(initial_model, initial_coefficients, &
            initial_transfers, full, out%model, out%coefficients, out%transfers)
         call impose_fixed_regression(initial_model, initial_coefficients, &
            initial_transfers, fixed_parameters, out%model, out%coefficients, &
            out%transfers)
         out%likelihood = arfima_regression_likelihood(series, out%model, regressors, &
            out%coefficients, out%transfers, transfer_regressors, selected_method, &
            profile_mean)
         if (out%likelihood%info /= 0) then
            out%info = out%likelihood%info
            return
         end if
         hessian = finite_difference_hessian(objective, optimization%parameters)
         call invert_matrix(hessian, free_covariance, hessian_info)
         valid_covariance = hessian_info == 0
         if (valid_covariance) then
            do i = 1, size(free_initial)
               if (free_covariance(i, i) <= 0.0_dp) valid_covariance = .false.
            end do
         end if
         if (valid_covariance) then
            do j = 1, size(free_initial)
               do i = 1, size(free_initial)
                  out%parameter_covariance(free_index(i), free_index(j)) = &
                     free_covariance(i, j)
               end do
               out%standard_error(free_index(j)) = sqrt(free_covariance(j, j))
            end do
         end if
         if (.not. optimization%converged) out%info = optimization%info
      end if
      if (out%likelihood%info /= 0) then
         out%info = out%likelihood%info
         return
      end if
      out%model%mean = out%likelihood%mean
      out%model%innovation_variance = out%likelihood%innovation_variance
      information_parameters = out%estimated_parameter_count + &
         merge(1, 0, profile_mean) + 1
      out%aic = -2.0_dp*out%likelihood%log_likelihood + &
         2.0_dp*real(information_parameters, dp)
      out%bic = -2.0_dp*out%likelihood%log_likelihood + &
         log(real(out%likelihood%observations, dp))*real(information_parameters, dp)
      if (out%likelihood%observations > information_parameters + 1) then
         out%aicc = out%aic + 2.0_dp*real(information_parameters* &
            (information_parameters + 1), dp)/ &
            real(out%likelihood%observations - information_parameters - 1, dp)
      end if

   contains

      pure function objective(free_parameters) result(value)
         ! Return the joint negative likelihood after inserting fixed values.
         real(dp), intent(in) :: free_parameters(:)
         real(dp) :: value
         real(dp) :: candidate_parameters(size(initial))
         type(arfima_model_t) :: candidate_model
         type(arfima_transfer_t), allocatable :: candidate_transfers(:)
         real(dp), allocatable :: candidate_coefficients(:)
         type(arfima_likelihood_t) :: candidate_likelihood

         candidate_parameters = initial
         candidate_parameters(free_index) = free_parameters
         call decode_regression_parameters(initial_model, initial_coefficients, &
            initial_transfers, candidate_parameters, candidate_model, &
            candidate_coefficients, candidate_transfers)
         call impose_fixed_regression(initial_model, initial_coefficients, &
            initial_transfers, fixed_parameters, candidate_model, &
            candidate_coefficients, candidate_transfers)
         candidate_likelihood = arfima_regression_likelihood(series, candidate_model, &
            regressors, candidate_coefficients, candidate_transfers, transfer_regressors, &
            selected_method, profile_mean)
         if (candidate_likelihood%info /= 0 .or. &
            .not. ieee_is_finite(candidate_likelihood%log_likelihood)) then
            value = huge(1.0_dp)/100.0_dp
         else
            value = -candidate_likelihood%log_likelihood
         end if
      end function objective

   end function arfima_regression_fit_fixed

   pure function arfima_fisher_information(model, exact, resolution, &
      observations) result(out)
      ! Compute spectral or truncated-score Fisher information.
      type(arfima_model_t), intent(in) :: model
      logical, intent(in), optional :: exact
      integer, intent(in), optional :: resolution, observations
      type(arfima_information_t) :: out
      real(dp), allocatable :: lower(:, :)
      integer :: inversion_info, status

      out%exact = .true.
      if (present(exact)) out%exact = exact
      out%resolution = merge(4096, 2048, out%exact)
      if (present(resolution)) out%resolution = resolution
      if (present(observations)) out%observations = observations
      out%parameters = information_parameter_count(model)
      if (.not. valid_model(model) .or. out%resolution < 16 .or. &
         out%observations < 0 .or. &
         (model%long_memory_type /= arfima_long_memory_none .and. &
         model%long_memory_type /= arfima_long_memory_fdwn) .or. &
         (model%period > 0 .and. &
         model%seasonal_long_memory_type /= arfima_long_memory_none .and. &
         model%seasonal_long_memory_type /= arfima_long_memory_fdwn)) then
         out%info = 1
         return
      end if
      if (out%parameters == 0) then
         allocate(out%information(0, 0), out%covariance(0, 0))
         out%positive_definite = .true.
         return
      end if
      if (out%exact) then
         out%information = spectral_information(model, out%resolution)
      else
         out%information = truncated_score_information(model, out%resolution)
      end if
      out%information = 0.5_dp*(out%information + transpose(out%information))
      call cholesky_lower(out%information, lower, status)
      out%positive_definite = status == 0
      if (out%observations > 0 .and. out%positive_definite) then
         call invert_matrix(real(out%observations, dp)*out%information, &
            out%covariance, inversion_info)
         if (inversion_info /= 0) out%info = 2
      else
         allocate(out%covariance(out%parameters, out%parameters))
         out%covariance = 0.0_dp
      end if
   end function arfima_fisher_information

   pure function arfima_fit_information(fit, exact, resolution) result(out)
      ! Compute information and covariance for a fitted ARFIMA model.
      type(arfima_fit_t), intent(in) :: fit
      logical, intent(in), optional :: exact
      integer, intent(in), optional :: resolution
      type(arfima_information_t) :: out
      logical :: selected_exact
      integer :: selected_resolution

      selected_exact = .true.
      if (present(exact)) selected_exact = exact
      selected_resolution = merge(4096, 2048, selected_exact)
      if (present(resolution)) selected_resolution = resolution
      if (fit%likelihood%info /= 0) then
         out%info = 1
         return
      end if
      out = arfima_fisher_information(fit%model, selected_exact, selected_resolution, &
         fit%likelihood%observations)
   end function arfima_fit_information

   pure function arfima_fit_diagnostics(series, fit, max_lag) result(out)
      ! Diagnose one fitted ARFIMA model from its one-step innovations.
      real(dp), intent(in) :: series(:)
      type(arfima_fit_t), intent(in) :: fit
      integer, intent(in), optional :: max_lag
      type(arfima_diagnostics_t) :: out
      type(arfima_likelihood_t) :: likelihood
      integer :: selected_lag

      if (fit%likelihood%info /= 0) then
         out%info = 1
         return
      end if
      likelihood = arfima_likelihood(series, fit%model, fit%likelihood%method, .false.)
      selected_lag = min(max(20, fit%parameter_count + 8), &
         max(0, likelihood%observations - 1))
      if (present(max_lag)) selected_lag = max_lag
      out = diagnostics_from_likelihood(series, likelihood, &
         fit%estimated_parameter_count, selected_lag)
   end function arfima_fit_diagnostics

   pure function arfima_regression_diagnostics(series, fit, regressors, &
      transfer_regressors, max_lag) result(out)
      ! Diagnose fitted regression and transfer-function ARFIMA innovations.
      real(dp), intent(in) :: series(:), regressors(:, :), transfer_regressors(:, :)
      type(arfima_regression_fit_t), intent(in) :: fit
      integer, intent(in), optional :: max_lag
      type(arfima_diagnostics_t) :: out
      type(arfima_likelihood_t) :: likelihood
      real(dp), allocatable :: effect(:), adjusted(:)
      integer :: selected_lag

      if (fit%likelihood%info /= 0 .or. &
         .not. valid_regression_inputs(series, regressors, fit%coefficients, &
         fit%transfers, transfer_regressors)) then
         out%info = 1
         return
      end if
      effect = regression_effect(regressors, fit%coefficients, fit%transfers, &
         transfer_regressors)
      adjusted = series - effect
      likelihood = arfima_likelihood(adjusted, fit%model, fit%likelihood%method, .false.)
      selected_lag = min(max(20, size(encode_model_parameters(fit%model)) + 8), &
         max(0, likelihood%observations - 1))
      if (present(max_lag)) selected_lag = max_lag
      out = diagnostics_from_likelihood(series, likelihood, &
         size(encode_model_parameters(fit%model)), selected_lag, adjusted)
   end function arfima_regression_diagnostics

   pure function arfima_mode_diagnostics(series, multifit, index, max_lag) result(out)
      ! Diagnose one selected mode from a multimode ARFIMA fit.
      real(dp), intent(in) :: series(:)
      type(arfima_multifit_t), intent(in) :: multifit
      integer, intent(in) :: index
      integer, intent(in), optional :: max_lag
      type(arfima_diagnostics_t) :: out

      if (index < 1 .or. index > size(multifit%modes)) then
         out%info = 1
      else if (present(max_lag)) then
         out = arfima_fit_diagnostics(series, multifit%modes(index), max_lag)
      else
         out = arfima_fit_diagnostics(series, multifit%modes(index))
      end if
   end function arfima_mode_diagnostics

   pure function arfima_fit_covariance(fit, exact, resolution) result(out)
      ! Extract observed and, when defined, expected covariance correlations.
      type(arfima_fit_t), intent(in) :: fit
      logical, intent(in), optional :: exact
      integer, intent(in), optional :: resolution
      type(arfima_covariance_t) :: out
      type(arfima_information_t) :: expected_information
      real(dp), allocatable :: jacobian(:, :), free_information(:, :)
      real(dp), allocatable :: free_covariance(:, :)
      integer, allocatable :: free_index(:)
      logical :: selected_exact
      integer :: i, j, inversion_info, selected_resolution

      if (fit%info /= 0 .or. .not. allocated(fit%parameter_covariance) .or. &
         .not. allocated(fit%unconstrained_parameters) .or. &
         .not. allocated(fit%fixed_parameters)) then
         out%info = 1
         return
      end if
      jacobian = operating_parameter_jacobian(fit)
      out%observed = matmul(matmul(jacobian, fit%parameter_covariance), &
         transpose(jacobian))
      out%observed_correlation = covariance_to_correlation(out%observed)
      allocate(out%expected(0, 0), out%expected_correlation(0, 0))
      selected_exact = .true.
      if (present(exact)) selected_exact = exact
      selected_resolution = merge(4096, 2048, selected_exact)
      if (present(resolution)) selected_resolution = resolution
      expected_information = arfima_fit_information(fit, selected_exact, &
         selected_resolution)
      if (expected_information%info == 0) then
         free_index = pack([(i, i=1, fit%parameter_count)], &
            .not. fit%fixed_parameters)
         deallocate(out%expected)
         allocate(out%expected(fit%parameter_count, fit%parameter_count))
         out%expected = 0.0_dp
         if (size(free_index) > 0) then
            allocate(free_information(size(free_index), size(free_index)))
            do j = 1, size(free_index)
               do i = 1, size(free_index)
                  free_information(i, j) = expected_information%information( &
                     free_index(i), free_index(j))
               end do
            end do
            call invert_matrix(real(fit%likelihood%observations, dp)*free_information, &
               free_covariance, inversion_info)
            if (inversion_info /= 0) return
            do j = 1, size(free_index)
               do i = 1, size(free_index)
                  out%expected(free_index(i), free_index(j)) = free_covariance(i, j)
               end do
            end do
         end if
         out%expected_correlation = covariance_to_correlation(out%expected)
         out%expected_available = .true.
      end if
   end function arfima_fit_covariance

   pure function arfima_average_forecast(series, multifit, horizon, level) result(out)
      ! Average mode forecasts and include between-mode mean uncertainty.
      real(dp), intent(in) :: series(:)
      type(arfima_multifit_t), intent(in) :: multifit
      integer, intent(in) :: horizon
      real(dp), intent(in), optional :: level
      type(arfima_forecast_t) :: out
      type(arfima_forecast_t), allocatable :: forecasts(:)
      real(dp), allocatable :: weights(:), difference(:)
      real(dp) :: selected_level, cutoff
      integer :: i, j

      selected_level = 0.95_dp
      if (present(level)) selected_level = level
      if (horizon < 1 .or. selected_level <= 0.0_dp .or. selected_level >= 1.0_dp .or. &
         size(multifit%modes) < 1) then
         out%info = 1
         return
      end if
      weights = mode_weights(multifit%modes)
      allocate(forecasts(size(multifit%modes)))
      do i = 1, size(forecasts)
         forecasts(i) = arfima_forecast(series, multifit%modes(i), horizon, &
            selected_level)
         if (forecasts(i)%info /= 0) then
            out%info = 2
            return
         end if
      end do
      out%horizon = horizon
      out%level = selected_level
      allocate(out%mean(horizon), out%covariance(horizon, horizon))
      allocate(out%stationary_mean(horizon))
      allocate(out%stationary_covariance(horizon, horizon))
      out%mean = 0.0_dp
      out%stationary_mean = 0.0_dp
      do i = 1, size(forecasts)
         out%mean = out%mean + weights(i)*forecasts(i)%mean
         out%stationary_mean = out%stationary_mean + &
            weights(i)*forecasts(i)%stationary_mean
      end do
      out%covariance = 0.0_dp
      out%stationary_covariance = 0.0_dp
      do i = 1, size(forecasts)
         difference = forecasts(i)%mean - out%mean
         out%covariance = out%covariance + weights(i)*forecasts(i)%covariance
         do j = 1, horizon
            out%covariance(:, j) = out%covariance(:, j) + &
               weights(i)*difference*difference(j)
         end do
         difference = forecasts(i)%stationary_mean - out%stationary_mean
         out%stationary_covariance = out%stationary_covariance + &
            weights(i)*forecasts(i)%stationary_covariance
         do j = 1, horizon
            out%stationary_covariance(:, j) = out%stationary_covariance(:, j) + &
               weights(i)*difference*difference(j)
         end do
      end do
      allocate(out%standard_error(horizon), out%lower(horizon), out%upper(horizon))
      do i = 1, horizon
         out%standard_error(i) = sqrt(max(0.0_dp, out%covariance(i, i)))
      end do
      cutoff = normal_quantile(0.5_dp + 0.5_dp*selected_level)
      out%lower = out%mean - cutoff*out%standard_error
      out%upper = out%mean + cutoff*out%standard_error
   end function arfima_average_forecast

   pure function arfima_forecast_model(series, model, horizon, level) result(out)
      ! Compute exact conditional Gaussian forecasts from an ARFIMA model.
      real(dp), intent(in) :: series(:)
      type(arfima_model_t), intent(in) :: model
      integer, intent(in) :: horizon
      real(dp), intent(in), optional :: level
      type(arfima_forecast_t) :: out
      type(arfima_acvf_t) :: covariance_result
      real(dp), allocatable :: work(:), past_covariance(:, :), inverse(:, :)
      real(dp), allocatable :: cross_covariance(:, :), future_covariance(:, :)
      real(dp), allocatable :: conditioning(:, :), initial(:), integrated(:)
      real(dp), allocatable :: integration_weight(:), transform(:, :)
      real(dp) :: cutoff
      integer :: active, i, j, n, status

      out%horizon = horizon
      if (present(level)) out%level = level
      if (horizon < 1 .or. out%level <= 0.0_dp .or. out%level >= 1.0_dp .or. &
         .not. valid_model(model) .or. .not. all(ieee_is_finite(series))) then
         out%info = 1
         return
      end if
      work = difference_series(series, model%difference_order, &
         model%seasonal_difference_order, model%period)
      n = size(work)
      if (n < 1) then
         out%info = 1
         return
      end if
      covariance_result = arfima_model_acvf(model, n + horizon - 1)
      if (covariance_result%info /= 0) then
         out%info = covariance_result%info
         return
      end if
      allocate(past_covariance(n, n), cross_covariance(horizon, n))
      allocate(future_covariance(horizon, horizon))
      do i = 1, n
         do j = 1, n
            past_covariance(i, j) = covariance_result%covariance(abs(i - j))
         end do
      end do
      do i = 1, horizon
         do j = 1, n
            cross_covariance(i, j) = covariance_result%covariance(n + i - j)
         end do
         do j = 1, horizon
            future_covariance(i, j) = covariance_result%covariance(abs(i - j))
         end do
      end do
      call invert_matrix(past_covariance, inverse, status)
      if (status /= 0) then
         out%info = 2
         return
      end if
      conditioning = matmul(cross_covariance, inverse)
      out%stationary_mean = model%mean + matmul(conditioning, work - model%mean)
      out%stationary_covariance = future_covariance - &
         matmul(conditioning, transpose(cross_covariance))
      out%stationary_covariance = 0.5_dp*(out%stationary_covariance + &
         transpose(out%stationary_covariance))

      active = model%difference_order + model%period*model%seasonal_difference_order
      if (active == 0) then
         out%mean = out%stationary_mean
         out%covariance = out%stationary_covariance
      else
         initial = series(size(series) - active + 1:)
         integrated = arfima_integrate(out%stationary_mean, initial, &
            model%difference_order, model%seasonal_difference_order, model%period)
         out%mean = integrated(active + 1:)
         allocate(integration_weight(0:horizon - 1), transform(horizon, horizon))
         integration_weight = inverse_difference_weights(model%difference_order, &
            model%seasonal_difference_order, model%period, horizon - 1)
         transform = 0.0_dp
         do i = 1, horizon
            do j = 1, i
               transform(i, j) = integration_weight(i - j)
            end do
         end do
         out%covariance = matmul(matmul(transform, out%stationary_covariance), &
            transpose(transform))
      end if
      allocate(out%standard_error(horizon), out%lower(horizon), out%upper(horizon))
      do i = 1, horizon
         out%standard_error(i) = sqrt(max(0.0_dp, out%covariance(i, i)))
      end do
      cutoff = normal_quantile(0.5_dp + 0.5_dp*out%level)
      out%lower = out%mean - cutoff*out%standard_error
      out%upper = out%mean + cutoff*out%standard_error
   end function arfima_forecast_model

   pure function arfima_forecast_fit(series, fit, horizon, level) result(out)
      ! Compute exact forecasts directly from a fitted ARFIMA object.
      real(dp), intent(in) :: series(:)
      type(arfima_fit_t), intent(in) :: fit
      integer, intent(in) :: horizon
      real(dp), intent(in), optional :: level
      type(arfima_forecast_t) :: out

      if (present(level)) then
         out = arfima_forecast_model(series, fit%model, horizon, level)
      else
         out = arfima_forecast_model(series, fit%model, horizon)
      end if
      if (fit%likelihood%info /= 0) out%info = 1
   end function arfima_forecast_fit

   function arfima_forecast_model_paths(series, model, horizon, simulations, &
      seed, level) result(out)
      ! Draw conditional Gaussian forecast paths from an ARFIMA model.
      real(dp), intent(in) :: series(:)
      type(arfima_model_t), intent(in) :: model
      integer, intent(in) :: horizon, simulations
      integer, intent(in), optional :: seed
      real(dp), intent(in), optional :: level
      type(arfima_forecast_t) :: out
      real(dp), allocatable :: standard(:), draw(:)
      integer :: i, simulation, status

      if (present(level)) then
         out = arfima_forecast_model(series, model, horizon, level)
      else
         out = arfima_forecast_model(series, model, horizon)
      end if
      if (out%info /= 0 .or. simulations < 1) then
         out%info = 1
         return
      end if
      if (present(seed)) call set_random_seed(seed)
      allocate(out%paths(horizon, simulations), standard(horizon), draw(horizon))
      do simulation = 1, simulations
         do i = 1, horizon
            standard(i) = random_standard_normal()
         end do
         call multivariate_normal_from_standard(out%mean, out%covariance, standard, &
            draw, status)
         if (status /= 0) then
            out%info = 2
            return
         end if
         out%paths(:, simulation) = draw
      end do
      out%simulations = simulations
   end function arfima_forecast_model_paths

   function arfima_forecast_fit_paths(series, fit, horizon, simulations, seed, &
      level) result(out)
      ! Draw conditional forecast paths directly from an ARFIMA fit.
      real(dp), intent(in) :: series(:)
      type(arfima_fit_t), intent(in) :: fit
      integer, intent(in) :: horizon, simulations
      integer, intent(in), optional :: seed
      real(dp), intent(in), optional :: level
      type(arfima_forecast_t) :: out

      if (present(seed) .and. present(level)) then
         out = arfima_forecast_model_paths(series, fit%model, horizon, simulations, &
            seed, level)
      else if (present(seed)) then
         out = arfima_forecast_model_paths(series, fit%model, horizon, simulations, &
            seed=seed)
      else if (present(level)) then
         out = arfima_forecast_model_paths(series, fit%model, horizon, simulations, &
            level=level)
      else
         out = arfima_forecast_model_paths(series, fit%model, horizon, simulations)
      end if
      if (fit%likelihood%info /= 0) out%info = 1
   end function arfima_forecast_fit_paths

   pure function arfima_regression_forecast(series, fit, regressors, future_regressors, &
      transfer_regressors, horizon, level) result(out)
      ! Forecast a fitted regression ARFIMA model using future regressors.
      real(dp), intent(in) :: series(:), regressors(:, :), future_regressors(:, :)
      type(arfima_regression_fit_t), intent(in) :: fit
      real(dp), intent(in) :: transfer_regressors(:, :)
      integer, intent(in) :: horizon
      real(dp), intent(in), optional :: level
      type(arfima_forecast_t) :: out
      real(dp), allocatable :: historical_effect(:), full_static(:, :), full_effect(:)
      integer :: n

      n = size(series)
      if (fit%likelihood%info /= 0 .or. horizon < 1 .or. &
         size(regressors, 1) /= n .or. size(regressors, 2) /= size(fit%coefficients) .or. &
         size(future_regressors, 1) /= horizon .or. &
         size(future_regressors, 2) /= size(fit%coefficients) .or. &
         size(transfer_regressors, 1) /= n + horizon .or. &
         size(transfer_regressors, 2) /= size(fit%transfers)) then
         out%info = 1
         return
      end if
      allocate(full_static(n + horizon, size(fit%coefficients)))
      if (size(fit%coefficients) > 0) then
         full_static(:n, :) = regressors
         full_static(n + 1:, :) = future_regressors
      end if
      full_effect = regression_effect(full_static, fit%coefficients, fit%transfers, &
         transfer_regressors)
      historical_effect = full_effect(:n)
      if (present(level)) then
         out = arfima_forecast_model(series - historical_effect, fit%model, horizon, level)
      else
         out = arfima_forecast_model(series - historical_effect, fit%model, horizon)
      end if
      if (out%info /= 0) return
      out%regression_mean = full_effect(n + 1:)
      out%mean = out%mean + out%regression_mean
      out%lower = out%lower + out%regression_mean
      out%upper = out%upper + out%regression_mean
   end function arfima_regression_forecast

   pure function arfima_regression_simulate_from_innovations(fit, innovations, &
      regressors, transfer_regressors, initial_values) result(out)
      ! Simulate a fitted regression ARFIMA model from supplied innovations.
      type(arfima_regression_fit_t), intent(in) :: fit
      real(dp), intent(in) :: innovations(:), regressors(:, :), transfer_regressors(:, :)
      real(dp), intent(in), optional :: initial_values(:)
      type(arfima_simulation_t) :: out
      real(dp), allocatable :: effect(:)

      if (fit%likelihood%info /= 0 .or. &
         .not. valid_regression_inputs(innovations, regressors, fit%coefficients, &
         fit%transfers, transfer_regressors)) then
         out%info = 1
         return
      end if
      if (present(initial_values)) then
         out = arfima_simulate_from_innovations(fit%model, innovations, initial_values)
      else
         out = arfima_simulate_from_innovations(fit%model, innovations)
      end if
      if (out%info /= 0) return
      effect = regression_effect(regressors, fit%coefficients, fit%transfers, &
         transfer_regressors)
      out%series = out%series + effect
   end function arfima_regression_simulate_from_innovations

   function arfima_regression_simulate(fit, regressors, transfer_regressors, seed, &
      initial_values) result(out)
      ! Simulate a fitted regression ARFIMA model with Gaussian innovations.
      type(arfima_regression_fit_t), intent(in) :: fit
      real(dp), intent(in) :: regressors(:, :), transfer_regressors(:, :)
      integer, intent(in), optional :: seed
      real(dp), intent(in), optional :: initial_values(:)
      type(arfima_simulation_t) :: out
      real(dp), allocatable :: innovations(:)
      integer :: i

      if (size(regressors, 1) < 1) then
         out%info = 1
         return
      end if
      if (present(seed)) call set_random_seed(seed)
      allocate(innovations(size(regressors, 1)))
      do i = 1, size(innovations)
         innovations(i) = random_standard_normal()
      end do
      if (present(initial_values)) then
         out = arfima_regression_simulate_from_innovations(fit, innovations, regressors, &
            transfer_regressors, initial_values)
      else
         out = arfima_regression_simulate_from_innovations(fit, innovations, regressors, &
            transfer_regressors)
      end if
   end function arfima_regression_simulate

   function random_model_starts(model, number_starts) result(starts)
      ! Generate stable starts in unconstrained PACF and memory coordinates.
      type(arfima_model_t), intent(in) :: model
      integer, intent(in) :: number_starts
      real(dp), allocatable :: starts(:, :)
      real(dp), allocatable :: initial(:)
      integer :: i, j

      initial = encode_model_parameters(model)
      allocate(starts(size(initial), max(0, number_starts)))
      if (number_starts < 1) return
      starts(:, 1) = initial
      do j = 2, number_starts
         do i = 1, size(initial)
            starts(i, j) = atanh(0.95_dp*(2.0_dp*random_uniform() - 1.0_dp))
         end do
      end do
   end function random_model_starts

   pure function combine_multifits(groups) result(out)
      ! Concatenate mode searches that may have different parameter counts.
      type(arfima_multifit_t), intent(in) :: groups(:)
      type(arfima_multifit_t) :: out
      integer :: i, j, mode_offset, start_offset, total_modes, total_starts
      integer :: maximum_parameters

      total_modes = 0
      total_starts = 0
      maximum_parameters = 0
      do i = 1, size(groups)
         total_modes = total_modes + size(groups(i)%modes)
         total_starts = total_starts + groups(i)%attempted
         maximum_parameters = max(maximum_parameters, &
            size(groups(i)%starting_parameters, 1))
         out%converged = out%converged + groups(i)%converged
      end do
      out%attempted = total_starts
      allocate(out%modes(total_modes), out%starting_index(total_modes))
      allocate(out%long_memory_type(total_modes), out%weight(total_modes))
      allocate(out%starting_parameters(maximum_parameters, total_starts))
      out%starting_parameters = 0.0_dp
      mode_offset = 0
      start_offset = 0
      do i = 1, size(groups)
         if (groups(i)%attempted > 0 .and. maximum_parameters > 0) then
            out%starting_parameters(:size(groups(i)%starting_parameters, 1), &
               start_offset + 1:start_offset + groups(i)%attempted) = &
               groups(i)%starting_parameters
         end if
         do j = 1, size(groups(i)%modes)
            out%modes(mode_offset + j) = groups(i)%modes(j)
            out%starting_index(mode_offset + j) = &
               start_offset + groups(i)%starting_index(j)
            out%long_memory_type(mode_offset + j) = &
               groups(i)%long_memory_type(j)
         end do
         mode_offset = mode_offset + size(groups(i)%modes)
         start_offset = start_offset + groups(i)%attempted
      end do
      out%weight = mode_weights(out%modes)
      if (total_modes == 0) out%info = 2
   end function combine_multifits

   pure function subset_multifit(multifit, indices) result(out)
      ! Copy selected modes while retaining search provenance.
      type(arfima_multifit_t), intent(in) :: multifit
      integer, intent(in) :: indices(:)
      type(arfima_multifit_t) :: out
      integer :: i

      out%attempted = multifit%attempted
      out%converged = multifit%converged
      out%info = multifit%info
      out%weeded = multifit%weeded
      out%starting_parameters = multifit%starting_parameters
      allocate(out%modes(size(indices)), out%starting_index(size(indices)))
      allocate(out%long_memory_type(size(indices)), out%weight(size(indices)))
      do i = 1, size(indices)
         out%modes(i) = multifit%modes(indices(i))
         out%starting_index(i) = multifit%starting_index(indices(i))
         out%long_memory_type(i) = multifit%long_memory_type(indices(i))
      end do
      out%weight = mode_weights(out%modes)
   end function subset_multifit

   pure function mode_weights(modes) result(weights)
      ! Normalize relative likelihoods into numerically stable mode weights.
      type(arfima_fit_t), intent(in) :: modes(:)
      real(dp), allocatable :: weights(:)
      real(dp) :: maximum_log_likelihood, total
      integer :: i

      allocate(weights(size(modes)))
      if (size(modes) == 0) return
      maximum_log_likelihood = maxval([(modes(i)%likelihood%log_likelihood, &
         i=1, size(modes))])
      do i = 1, size(modes)
         weights(i) = exp(modes(i)%likelihood%log_likelihood - maximum_log_likelihood)
      end do
      total = sum(weights)
      if (total > 0.0_dp) weights = weights/total
   end function mode_weights

   pure function likelihood_order(modes) result(order)
      ! Return mode indices sorted by decreasing log likelihood.
      type(arfima_fit_t), intent(in) :: modes(:)
      integer, allocatable :: order(:)
      integer :: i, j, best, temporary

      order = [(i, i=1, size(modes))]
      do i = 1, size(order) - 1
         best = i
         do j = i + 1, size(order)
            if (modes(order(j))%likelihood%log_likelihood > &
               modes(order(best))%likelihood%log_likelihood) best = j
         end do
         temporary = order(i)
         order(i) = order(best)
         order(best) = temporary
      end do
   end function likelihood_order

   pure function mode_parameters(model, transformed) result(parameters)
      ! Flatten operating coefficients or PACFs and long-memory parameters.
      type(arfima_model_t), intent(in) :: model
      logical, intent(in) :: transformed
      real(dp), allocatable :: parameters(:), block(:)
      integer :: count, offset

      count = size(model%ar) + size(model%theta) + size(model%seasonal_ar) + &
         size(model%seasonal_theta) + &
         merge(1, 0, model%long_memory_type /= arfima_long_memory_none) + &
         merge(1, 0, model%period > 0 .and. &
         model%seasonal_long_memory_type /= arfima_long_memory_none)
      allocate(parameters(count))
      offset = 0
      if (size(model%ar) > 0) then
         block = model%ar
         if (transformed) block = ar_to_partial(model%ar)
         parameters(offset + 1:offset + size(block)) = block
         offset = offset + size(block)
      end if
      if (size(model%theta) > 0) then
         block = model%theta
         if (transformed) block = ar_to_partial(model%theta)
         parameters(offset + 1:offset + size(block)) = block
         offset = offset + size(block)
      end if
      if (size(model%seasonal_ar) > 0) then
         block = model%seasonal_ar
         if (transformed) block = ar_to_partial(model%seasonal_ar)
         parameters(offset + 1:offset + size(block)) = block
         offset = offset + size(block)
      end if
      if (size(model%seasonal_theta) > 0) then
         block = model%seasonal_theta
         if (transformed) block = ar_to_partial(model%seasonal_theta)
         parameters(offset + 1:offset + size(block)) = block
         offset = offset + size(block)
      end if
      if (model%long_memory_type /= arfima_long_memory_none) then
         offset = offset + 1
         parameters(offset) = model%long_memory_parameter
      end if
      if (model%period > 0 .and. &
         model%seasonal_long_memory_type /= arfima_long_memory_none) then
         offset = offset + 1
         parameters(offset) = model%seasonal_long_memory_parameter
      end if

   end function mode_parameters

   pure logical function compatible_mode_models(first, second) result(compatible)
      ! Test whether two fitted modes share one parameterization.
      type(arfima_model_t), intent(in) :: first, second

      compatible = size(first%ar) == size(second%ar) .and. &
         size(first%theta) == size(second%theta) .and. &
         size(first%seasonal_ar) == size(second%seasonal_ar) .and. &
         size(first%seasonal_theta) == size(second%seasonal_theta) .and. &
         first%long_memory_type == second%long_memory_type .and. &
         first%seasonal_long_memory_type == second%seasonal_long_memory_type .and. &
         first%period == second%period
   end function compatible_mode_models

   pure logical function modes_share_wall(first, second, tolerance) result(shared)
      ! Test whether two compatible modes occupy the same parameter boundary.
      type(arfima_fit_t), intent(in) :: first, second
      real(dp), intent(in) :: tolerance
      real(dp), allocatable :: first_parameters(:), second_parameters(:)
      real(dp), allocatable :: lower(:), upper(:)
      integer :: i

      shared = .false.
      if (.not. compatible_mode_models(first%model, second%model)) return
      first_parameters = mode_parameters(first%model, .true.)
      second_parameters = mode_parameters(second%model, .true.)
      call mode_parameter_bounds(first%model, lower, upper)
      do i = 1, size(first_parameters)
         if (abs(first_parameters(i) - lower(i)) < tolerance .and. &
            abs(second_parameters(i) - lower(i)) < tolerance) then
            shared = .true.
            return
         end if
         if (abs(first_parameters(i) - upper(i)) < tolerance .and. &
            abs(second_parameters(i) - upper(i)) < tolerance) then
            shared = .true.
            return
         end if
      end do
   end function modes_share_wall

   pure subroutine mode_parameter_bounds(model, lower, upper)
      ! Return PACF and long-memory bounds in flattened mode order.
      type(arfima_model_t), intent(in) :: model
      real(dp), allocatable, intent(out) :: lower(:), upper(:)
      real(dp) :: midpoint, half_width
      integer :: arma_count, count, offset

      arma_count = size(model%ar) + size(model%theta) + size(model%seasonal_ar) + &
         size(model%seasonal_theta)
      count = size(mode_parameters(model, .true.))
      allocate(lower(count), upper(count))
      lower = -1.0_dp
      upper = 1.0_dp
      offset = arma_count
      if (model%long_memory_type /= arfima_long_memory_none) then
         offset = offset + 1
         call long_memory_bounds(model%long_memory_type, midpoint, half_width)
         lower(offset) = midpoint - half_width
         upper(offset) = midpoint + half_width
      end if
      if (model%period > 0 .and. &
         model%seasonal_long_memory_type /= arfima_long_memory_none) then
         offset = offset + 1
         call long_memory_bounds(model%seasonal_long_memory_type, midpoint, half_width)
         lower(offset) = midpoint - half_width
         upper(offset) = midpoint + half_width
      end if
   end subroutine mode_parameter_bounds

   pure subroutine allocate_empty_multifit(out)
      ! Allocate empty mode arrays for an invalid or unsuccessful search.
      type(arfima_multifit_t), intent(inout) :: out

      if (.not. allocated(out%modes)) allocate(out%modes(0))
      if (.not. allocated(out%starting_index)) allocate(out%starting_index(0))
      if (.not. allocated(out%long_memory_type)) allocate(out%long_memory_type(0))
      if (.not. allocated(out%weight)) allocate(out%weight(0))
      if (.not. allocated(out%starting_parameters)) then
         allocate(out%starting_parameters(0, 0))
      end if
   end subroutine allocate_empty_multifit

   pure function encode_model_parameters(model) result(parameters)
      ! Map stable model parameters to unconstrained optimizer coordinates.
      type(arfima_model_t), intent(in) :: model
      real(dp), allocatable :: parameters(:)
      real(dp), allocatable :: partial(:)
      integer :: count, offset

      count = size(model%ar) + size(model%theta) + size(model%seasonal_ar) + &
         size(model%seasonal_theta) + &
         merge(1, 0, model%long_memory_type /= arfima_long_memory_none) + &
         merge(1, 0, model%period > 0 .and. &
         model%seasonal_long_memory_type /= arfima_long_memory_none)
      allocate(parameters(count))
      offset = 0
      if (size(model%ar) > 0) then
         partial = ar_to_partial(model%ar)
         parameters(offset + 1:offset + size(partial)) = inverse_hyperbolic(partial)
         offset = offset + size(partial)
      end if
      if (size(model%theta) > 0) then
         partial = ar_to_partial(model%theta)
         parameters(offset + 1:offset + size(partial)) = inverse_hyperbolic(partial)
         offset = offset + size(partial)
      end if
      if (size(model%seasonal_ar) > 0) then
         partial = ar_to_partial(model%seasonal_ar)
         parameters(offset + 1:offset + size(partial)) = inverse_hyperbolic(partial)
         offset = offset + size(partial)
      end if
      if (size(model%seasonal_theta) > 0) then
         partial = ar_to_partial(model%seasonal_theta)
         parameters(offset + 1:offset + size(partial)) = inverse_hyperbolic(partial)
         offset = offset + size(partial)
      end if
      if (model%long_memory_type /= arfima_long_memory_none) then
         offset = offset + 1
         parameters(offset) = encode_long_memory(model%long_memory_type, &
            model%long_memory_parameter)
      end if
      if (model%period > 0 .and. &
         model%seasonal_long_memory_type /= arfima_long_memory_none) then
         offset = offset + 1
         parameters(offset) = encode_long_memory(model%seasonal_long_memory_type, &
            model%seasonal_long_memory_parameter)
      end if
   end function encode_model_parameters

   pure function decode_model_parameters(template, parameters) result(model)
      ! Convert unconstrained optimizer coordinates into a stable model.
      type(arfima_model_t), intent(in) :: template
      real(dp), intent(in) :: parameters(:)
      type(arfima_model_t) :: model
      integer :: offset, width

      model = template
      if (size(parameters) /= size(encode_model_parameters(template))) then
         model%info = 1
         return
      end if
      offset = 0
      width = size(model%ar)
      if (width > 0) then
         model%ar = durbin_levinson_coefficients(tanh(parameters(offset + 1:offset + width)))
         offset = offset + width
      end if
      width = size(model%theta)
      if (width > 0) then
         model%theta = durbin_levinson_coefficients(tanh(parameters(offset + 1:offset + width)))
         offset = offset + width
      end if
      width = size(model%seasonal_ar)
      if (width > 0) then
         model%seasonal_ar = &
            durbin_levinson_coefficients(tanh(parameters(offset + 1:offset + width)))
         offset = offset + width
      end if
      width = size(model%seasonal_theta)
      if (width > 0) then
         model%seasonal_theta = &
            durbin_levinson_coefficients(tanh(parameters(offset + 1:offset + width)))
         offset = offset + width
      end if
      if (model%long_memory_type /= arfima_long_memory_none) then
         offset = offset + 1
         model%long_memory_parameter = decode_long_memory(model%long_memory_type, &
            parameters(offset))
      end if
      if (model%period > 0 .and. &
         model%seasonal_long_memory_type /= arfima_long_memory_none) then
         offset = offset + 1
         model%seasonal_long_memory_parameter = decode_long_memory( &
            model%seasonal_long_memory_type, parameters(offset))
      end if
      model%info = 0
   end function decode_model_parameters

   pure function encode_regression_parameters(model, coefficients, transfers) &
      result(parameters)
      ! Flatten joint model parameters into unconstrained optimizer coordinates.
      type(arfima_model_t), intent(in) :: model
      real(dp), intent(in) :: coefficients(:)
      type(arfima_transfer_t), intent(in) :: transfers(:)
      real(dp), allocatable :: parameters(:)
      real(dp), allocatable :: model_parameters(:), partial(:)
      integer :: count, i, offset

      model_parameters = encode_model_parameters(model)
      count = size(model_parameters) + size(coefficients)
      do i = 1, size(transfers)
         count = count + size(transfers(i)%denominator) + size(transfers(i)%numerator)
      end do
      allocate(parameters(count))
      offset = size(model_parameters)
      if (offset > 0) parameters(:offset) = model_parameters
      if (size(coefficients) > 0) then
         parameters(offset + 1:offset + size(coefficients)) = coefficients
         offset = offset + size(coefficients)
      end if
      do i = 1, size(transfers)
         if (size(transfers(i)%denominator) > 0) then
            partial = ar_to_partial(transfers(i)%denominator)
            parameters(offset + 1:offset + size(partial)) = inverse_hyperbolic(partial)
            offset = offset + size(partial)
         end if
         parameters(offset + 1:offset + size(transfers(i)%numerator)) = &
            transfers(i)%numerator
         offset = offset + size(transfers(i)%numerator)
      end do
   end function encode_regression_parameters

   pure subroutine decode_regression_parameters(model_template, coefficient_template, &
      transfer_template, parameters, model, coefficients, transfers)
      ! Reconstruct joint model objects from optimizer coordinates.
      type(arfima_model_t), intent(in) :: model_template
      real(dp), intent(in) :: coefficient_template(:), parameters(:)
      type(arfima_transfer_t), intent(in) :: transfer_template(:)
      type(arfima_model_t), intent(out) :: model
      real(dp), allocatable, intent(out) :: coefficients(:)
      type(arfima_transfer_t), allocatable, intent(out) :: transfers(:)
      integer :: i, model_count, offset, width

      model_count = size(encode_model_parameters(model_template))
      model = decode_model_parameters(model_template, parameters(:model_count))
      allocate(coefficients(size(coefficient_template)), transfers(size(transfer_template)))
      offset = model_count
      if (size(coefficients) > 0) then
         coefficients = parameters(offset + 1:offset + size(coefficients))
         offset = offset + size(coefficients)
      end if
      do i = 1, size(transfers)
         transfers(i) = transfer_template(i)
         width = size(transfers(i)%denominator)
         if (width > 0) then
            transfers(i)%denominator = &
               durbin_levinson_coefficients(tanh(parameters(offset + 1:offset + width)))
            offset = offset + width
         end if
         width = size(transfers(i)%numerator)
         transfers(i)%numerator = parameters(offset + 1:offset + width)
         offset = offset + width
         transfers(i)%info = 0
      end do
   end subroutine decode_regression_parameters

   pure subroutine impose_fixed_model(template, fixed_parameters, model)
      ! Restore fixed operating coefficients after a stable-coordinate decode.
      type(arfima_model_t), intent(in) :: template
      logical, intent(in) :: fixed_parameters(:)
      type(arfima_model_t), intent(inout) :: model
      integer :: offset, width

      offset = 0
      width = size(model%ar)
      if (width > 0) then
         where (fixed_parameters(offset + 1:offset + width)) model%ar = template%ar
         offset = offset + width
      end if
      width = size(model%theta)
      if (width > 0) then
         where (fixed_parameters(offset + 1:offset + width)) model%theta = template%theta
         offset = offset + width
      end if
      width = size(model%seasonal_ar)
      if (width > 0) then
         where (fixed_parameters(offset + 1:offset + width))
            model%seasonal_ar = template%seasonal_ar
         end where
         offset = offset + width
      end if
      width = size(model%seasonal_theta)
      if (width > 0) then
         where (fixed_parameters(offset + 1:offset + width))
            model%seasonal_theta = template%seasonal_theta
         end where
         offset = offset + width
      end if
      if (model%long_memory_type /= arfima_long_memory_none) then
         offset = offset + 1
         if (fixed_parameters(offset)) then
            model%long_memory_parameter = template%long_memory_parameter
         end if
      end if
      if (model%period > 0 .and. &
         model%seasonal_long_memory_type /= arfima_long_memory_none) then
         offset = offset + 1
         if (fixed_parameters(offset)) then
            model%seasonal_long_memory_parameter = &
               template%seasonal_long_memory_parameter
         end if
      end if
   end subroutine impose_fixed_model

   pure subroutine impose_fixed_regression(model_template, coefficient_template, &
      transfer_template, fixed_parameters, model, coefficients, transfers)
      ! Restore fixed model, regression, and transfer operating coefficients.
      type(arfima_model_t), intent(in) :: model_template
      real(dp), intent(in) :: coefficient_template(:)
      type(arfima_transfer_t), intent(in) :: transfer_template(:)
      logical, intent(in) :: fixed_parameters(:)
      type(arfima_model_t), intent(inout) :: model
      real(dp), intent(inout) :: coefficients(:)
      type(arfima_transfer_t), intent(inout) :: transfers(:)
      integer :: i, model_count, offset, width

      model_count = size(encode_model_parameters(model_template))
      call impose_fixed_model(model_template, fixed_parameters(:model_count), model)
      offset = model_count
      width = size(coefficients)
      if (width > 0) then
         where (fixed_parameters(offset + 1:offset + width))
            coefficients = coefficient_template
         end where
         offset = offset + width
      end if
      do i = 1, size(transfers)
         width = size(transfers(i)%denominator)
         if (width > 0) then
            where (fixed_parameters(offset + 1:offset + width))
               transfers(i)%denominator = transfer_template(i)%denominator
            end where
            offset = offset + width
         end if
         width = size(transfers(i)%numerator)
         if (width > 0) then
            where (fixed_parameters(offset + 1:offset + width))
               transfers(i)%numerator = transfer_template(i)%numerator
            end where
            offset = offset + width
         end if
      end do
   end subroutine impose_fixed_regression

   pure integer function information_parameter_count(model) result(count)
      ! Count dynamic coefficients represented in ARFIMA information matrices.
      type(arfima_model_t), intent(in) :: model

      count = size(model%ar) + size(model%theta) + size(model%seasonal_ar) + &
         size(model%seasonal_theta)
      if (model%long_memory_type == arfima_long_memory_fdwn) count = count + 1
      if (model%period > 0 .and. &
         model%seasonal_long_memory_type == arfima_long_memory_fdwn) count = count + 1
   end function information_parameter_count

   pure function spectral_information(model, grid_size) result(information)
      ! Integrate log-spectrum score products over positive frequencies.
      type(arfima_model_t), intent(in) :: model
      integer, intent(in) :: grid_size
      real(dp), allocatable :: information(:, :)
      real(dp), allocatable :: score(:)
      complex(dp) :: z, ar_value, theta_value, seasonal_ar_value, seasonal_theta_value
      real(dp) :: frequency
      integer :: count, frequency_index, j, offset

      count = information_parameter_count(model)
      allocate(information(count, count), score(count))
      information = 0.0_dp
      do frequency_index = 1, grid_size
         frequency = (real(frequency_index, dp) - 0.5_dp)*acos(-1.0_dp)/ &
            real(grid_size, dp)
         z = cmplx(cos(frequency), -sin(frequency), dp)
         ar_value = lag_polynomial(model%ar, z, 1)
         theta_value = lag_polynomial(model%theta, z, 1)
         seasonal_ar_value = lag_polynomial(model%seasonal_ar, z, model%period)
         seasonal_theta_value = lag_polynomial(model%seasonal_theta, z, model%period)
         offset = 0
         do j = 1, size(model%ar)
            offset = offset + 1
            score(offset) = 2.0_dp*real(z**j/ar_value, dp)
         end do
         do j = 1, size(model%theta)
            offset = offset + 1
            score(offset) = -2.0_dp*real(z**j/theta_value, dp)
         end do
         do j = 1, size(model%seasonal_ar)
            offset = offset + 1
            score(offset) = 2.0_dp*real(z**(j*model%period)/seasonal_ar_value, dp)
         end do
         do j = 1, size(model%seasonal_theta)
            offset = offset + 1
            score(offset) = -2.0_dp*real(z**(j*model%period)/ &
               seasonal_theta_value, dp)
         end do
         if (model%long_memory_type == arfima_long_memory_fdwn) then
            offset = offset + 1
            score(offset) = -2.0_dp*log(abs(1.0_dp - z))
         end if
         if (model%period > 0 .and. &
            model%seasonal_long_memory_type == arfima_long_memory_fdwn) then
            offset = offset + 1
            score(offset) = -2.0_dp*log(abs(1.0_dp - z**model%period))
         end if
         information = information + spread(score, 2, count)*spread(score, 1, count)
      end do
      information = 0.5_dp*information/real(grid_size, dp)
   end function spectral_information

   pure function truncated_score_information(model, max_lag) result(information)
      ! Approximate information using finite causal score-filter weights.
      type(arfima_model_t), intent(in) :: model
      integer, intent(in) :: max_lag
      real(dp), allocatable :: information(:, :)
      real(dp), allocatable :: score(:, :), inverse_weights(:)
      integer :: count, j, lag, offset, shift

      count = information_parameter_count(model)
      allocate(score(count, max_lag))
      score = 0.0_dp
      offset = 0
      inverse_weights = inverse_operator_weights(model%ar, 1, max_lag)
      do j = 1, size(model%ar)
         offset = offset + 1
         do lag = j, max_lag
            score(offset, lag) = inverse_weights(lag - j + 1)
         end do
      end do
      inverse_weights = inverse_operator_weights(model%theta, 1, max_lag)
      do j = 1, size(model%theta)
         offset = offset + 1
         do lag = j, max_lag
            score(offset, lag) = -inverse_weights(lag - j + 1)
         end do
      end do
      inverse_weights = inverse_operator_weights(model%seasonal_ar, model%period, max_lag)
      do j = 1, size(model%seasonal_ar)
         offset = offset + 1
         shift = j*model%period
         do lag = shift, max_lag
            score(offset, lag) = inverse_weights(lag - shift + 1)
         end do
      end do
      inverse_weights = inverse_operator_weights(model%seasonal_theta, model%period, max_lag)
      do j = 1, size(model%seasonal_theta)
         offset = offset + 1
         shift = j*model%period
         do lag = shift, max_lag
            score(offset, lag) = -inverse_weights(lag - shift + 1)
         end do
      end do
      if (model%long_memory_type == arfima_long_memory_fdwn) then
         offset = offset + 1
         do lag = 1, max_lag
            score(offset, lag) = 1.0_dp/real(lag, dp)
         end do
      end if
      if (model%period > 0 .and. &
         model%seasonal_long_memory_type == arfima_long_memory_fdwn) then
         offset = offset + 1
         do lag = model%period, max_lag, model%period
            score(offset, lag) = 1.0_dp/real(lag/model%period, dp)
         end do
      end if
      information = matmul(score, transpose(score))
   end function truncated_score_information

   pure complex(dp) function lag_polynomial(coefficients, z, period) result(value)
      ! Evaluate 1-sum(coefficients(j)*z**(j*period)).
      real(dp), intent(in) :: coefficients(:)
      complex(dp), intent(in) :: z
      integer, intent(in) :: period
      integer :: j

      value = cmplx(1.0_dp, 0.0_dp, dp)
      do j = 1, size(coefficients)
         value = value - coefficients(j)*z**(j*period)
      end do
   end function lag_polynomial

   pure function inverse_operator_weights(coefficients, period, max_lag) result(weights)
      ! Expand the inverse of one ordinary or seasonal AR-form polynomial.
      real(dp), intent(in) :: coefficients(:)
      integer, intent(in) :: period, max_lag
      real(dp), allocatable :: weights(:), base(:)
      integer :: lag

      allocate(weights(0:max_lag))
      weights = 0.0_dp
      if (size(coefficients) == 0) then
         weights(0) = 1.0_dp
         return
      end if
      base = arfima_psi_weights(coefficients, [real(dp) ::], 0.0_dp, max_lag/period)
      do lag = 0, max_lag/period
         weights(lag*period) = base(lag + 1)
      end do
   end function inverse_operator_weights

   pure function combined_operator_coefficients(ordinary, seasonal, period) result(coefficients)
      ! Multiply ordinary and expanded seasonal AR-form polynomials.
      real(dp), intent(in) :: ordinary(:), seasonal(:)
      integer, intent(in) :: period
      real(dp), allocatable :: coefficients(:)
      real(dp), allocatable :: ordinary_polynomial(:), seasonal_polynomial(:), product(:)
      integer :: degree, j, seasonal_degree

      if (size(ordinary) == 0 .and. size(seasonal) == 0) then
         allocate(coefficients(0))
         return
      end if
      seasonal_degree = 0
      if (period > 0) seasonal_degree = size(seasonal)*period
      degree = size(ordinary) + seasonal_degree
      allocate(ordinary_polynomial(0:degree), seasonal_polynomial(0:degree))
      allocate(product(0:degree), coefficients(degree))
      ordinary_polynomial = 0.0_dp
      seasonal_polynomial = 0.0_dp
      ordinary_polynomial(0) = 1.0_dp
      seasonal_polynomial(0) = 1.0_dp
      if (size(ordinary) > 0) ordinary_polynomial(1:size(ordinary)) = -ordinary
      if (period > 0) then
         do j = 1, size(seasonal)
            seasonal_polynomial(j*period) = -seasonal(j)
         end do
      end if
      product = polynomial_product_truncated(ordinary_polynomial, &
         seasonal_polynomial, degree)
      coefficients = -product(1:)
   end function combined_operator_coefficients

   pure function ar_to_partial(coefficients) result(partial)
      ! Convert stable AR coefficients to partial autocorrelations.
      real(dp), intent(in) :: coefficients(:)
      real(dp) :: partial(size(coefficients)), work(size(coefficients))
      real(dp) :: previous(size(coefficients)), reflection, denominator
      integer :: j, order

      work = coefficients
      partial = 0.0_dp
      do order = size(coefficients), 1, -1
         reflection = work(order)
         partial(order) = reflection
         if (order == 1) exit
         denominator = 1.0_dp - reflection**2
         previous = 0.0_dp
         do j = 1, order - 1
            previous(j) = (work(j) + reflection*work(order - j))/denominator
         end do
         work(:order - 1) = previous(:order - 1)
      end do
   end function ar_to_partial

   pure function inverse_hyperbolic(values) result(transformed)
      ! Apply a guarded inverse hyperbolic tangent elementally.
      real(dp), intent(in) :: values(:)
      real(dp) :: transformed(size(values)), bounded
      integer :: i

      do i = 1, size(values)
         bounded = max(-1.0_dp + 1.0e-10_dp, min(1.0_dp - 1.0e-10_dp, values(i)))
         transformed(i) = atanh(bounded)
      end do
   end function inverse_hyperbolic

   pure real(dp) function encode_long_memory(long_memory_type, parameter) result(value)
      ! Map one bounded long-memory parameter to the real line.
      integer, intent(in) :: long_memory_type
      real(dp), intent(in) :: parameter
      real(dp) :: midpoint, half_width, scaled

      call long_memory_bounds(long_memory_type, midpoint, half_width)
      scaled = (parameter - midpoint)/half_width
      scaled = max(-1.0_dp + 1.0e-10_dp, min(1.0_dp - 1.0e-10_dp, scaled))
      value = atanh(scaled)
   end function encode_long_memory

   pure real(dp) function decode_long_memory(long_memory_type, value) result(parameter)
      ! Map one optimizer coordinate into its long-memory domain.
      integer, intent(in) :: long_memory_type
      real(dp), intent(in) :: value
      real(dp) :: midpoint, half_width

      call long_memory_bounds(long_memory_type, midpoint, half_width)
      parameter = midpoint + half_width*tanh(value)
   end function decode_long_memory

   pure subroutine long_memory_bounds(long_memory_type, midpoint, half_width)
      ! Return the midpoint and half-width of a long-memory domain.
      integer, intent(in) :: long_memory_type
      real(dp), intent(out) :: midpoint, half_width

      select case (long_memory_type)
      case (arfima_long_memory_fdwn)
         midpoint = -0.25_dp
         half_width = 0.75_dp
      case (arfima_long_memory_fgn)
         midpoint = 0.5_dp
         half_width = 0.5_dp
      case (arfima_long_memory_hd)
         midpoint = 1.5_dp
         half_width = 1.5_dp
      case default
         midpoint = 0.0_dp
         half_width = 1.0_dp
      end select
   end subroutine long_memory_bounds

   pure function difference_series(series, difference_order, &
      seasonal_difference_order, period) result(differenced)
      ! Apply ordinary and seasonal integer differencing in one polynomial.
      real(dp), intent(in) :: series(:)
      integer, intent(in) :: difference_order, seasonal_difference_order, period
      real(dp), allocatable :: differenced(:)
      real(dp), allocatable :: ordinary(:), seasonal(:), coefficient(:)
      integer :: active, i, lag

      active = difference_order + period*seasonal_difference_order
      if (difference_order < 0 .or. seasonal_difference_order < 0 .or. &
         (seasonal_difference_order > 0 .and. period < 2) .or. &
         size(series) <= active) then
         allocate(differenced(0))
         return
      end if
      if (active == 0) then
         differenced = series
         return
      end if
      allocate(ordinary(0:active), seasonal(0:active), coefficient(0:active))
      ordinary = arfima_fractional_weights(real(difference_order, dp), active)
      seasonal = 0.0_dp
      seasonal(0) = 1.0_dp
      if (seasonal_difference_order > 0) then
         seasonal = arfima_fractional_weights(real(seasonal_difference_order, dp), &
            active, period)
      end if
      coefficient = polynomial_product_truncated(ordinary, seasonal, active)
      allocate(differenced(size(series) - active))
      do i = active + 1, size(series)
         differenced(i - active) = 0.0_dp
         do lag = 0, active
            differenced(i - active) = differenced(i - active) + &
               coefficient(lag)*series(i - lag)
         end do
      end do
   end function difference_series

   pure function inverse_difference_weights(difference_order, &
      seasonal_difference_order, period, max_lag) result(weights)
      ! Expand the inverse ordinary and seasonal differencing polynomial.
      integer, intent(in) :: difference_order, seasonal_difference_order, period, max_lag
      real(dp), allocatable :: weights(:)
      real(dp), allocatable :: ordinary(:), seasonal(:)

      if (max_lag < 0 .or. difference_order < 0 .or. seasonal_difference_order < 0 .or. &
         (seasonal_difference_order > 0 .and. period < 2)) then
         allocate(weights(0))
         return
      end if
      allocate(ordinary(0:max_lag), seasonal(0:max_lag), weights(0:max_lag))
      ordinary = arfima_fractional_weights(-real(difference_order, dp), max_lag)
      seasonal = 0.0_dp
      seasonal(0) = 1.0_dp
      if (seasonal_difference_order > 0) then
         seasonal = arfima_fractional_weights(-real(seasonal_difference_order, dp), &
            max_lag, period)
      end if
      weights = polynomial_product_truncated(ordinary, seasonal, max_lag)
   end function inverse_difference_weights

   pure logical function valid_transfer(transfer) result(valid)
      ! Check one transfer function for finite stable coefficients.
      type(arfima_transfer_t), intent(in) :: transfer

      valid = allocated(transfer%denominator) .and. allocated(transfer%numerator)
      if (.not. valid) return
      valid = size(transfer%numerator) > 0 .and. transfer%delay >= 0 .and. &
         all(ieee_is_finite(transfer%denominator)) .and. &
         all(ieee_is_finite(transfer%numerator)) .and. &
         stable_coefficients(transfer%denominator)
   end function valid_transfer

   pure logical function valid_regression_inputs(series, regressors, coefficients, &
      transfers, transfer_regressors) result(valid)
      ! Check static and dynamic regression dimensions and values.
      real(dp), intent(in) :: series(:), regressors(:, :), coefficients(:)
      type(arfima_transfer_t), intent(in) :: transfers(:)
      real(dp), intent(in) :: transfer_regressors(:, :)
      integer :: i

      valid = size(series) > 1 .and. size(regressors, 1) == size(series) .and. &
         size(regressors, 2) == size(coefficients) .and. &
         size(transfer_regressors, 1) == size(series) .and. &
         size(transfer_regressors, 2) == size(transfers) .and. &
         all(ieee_is_finite(series)) .and. all(ieee_is_finite(regressors)) .and. &
         all(ieee_is_finite(coefficients)) .and. &
         all(ieee_is_finite(transfer_regressors))
      if (.not. valid) return
      do i = 1, size(transfers)
         if (.not. valid_transfer(transfers(i))) then
            valid = .false.
            return
         end if
      end do
   end function valid_regression_inputs

   pure function regression_effect(regressors, coefficients, transfers, &
      transfer_regressors) result(effect)
      ! Combine static regression and dynamic transfer-function effects.
      real(dp), intent(in) :: regressors(:, :), coefficients(:)
      type(arfima_transfer_t), intent(in) :: transfers(:)
      real(dp), intent(in) :: transfer_regressors(:, :)
      real(dp), allocatable :: effect(:), response(:)
      integer :: i

      allocate(effect(size(regressors, 1)))
      effect = 0.0_dp
      if (size(coefficients) > 0) effect = matmul(regressors, coefficients)
      do i = 1, size(transfers)
         response = arfima_transfer_response(transfer_regressors(:, i), transfers(i))
         effect = effect + response
      end do
   end function regression_effect

   pure function diagnostics_from_likelihood(series, likelihood, parameter_count, &
      max_lag, regression_residuals) result(out)
      ! Assemble fitted values and numerical innovation diagnostics.
      real(dp), intent(in) :: series(:)
      type(arfima_likelihood_t), intent(in) :: likelihood
      integer, intent(in) :: parameter_count, max_lag
      real(dp), intent(in), optional :: regression_residuals(:)
      type(arfima_diagnostics_t) :: out
      real(dp), allocatable :: squared(:)
      real(dp) :: statistic, squared_statistic, probability
      integer :: n, first_lag, count_lags, h, i

      n = likelihood%observations
      out%observations = n
      first_lag = parameter_count + 1
      if (likelihood%info /= 0 .or. n < 3 .or. size(series) < n .or. &
         max_lag < first_lag .or. max_lag >= n .or. &
         .not. allocated(likelihood%innovation) .or. &
         .not. allocated(likelihood%prediction_variance)) then
         out%info = 1
         return
      end if
      out%residuals = likelihood%innovation
      out%prediction_variance = likelihood%prediction_variance
      out%standardized_residuals = out%residuals/sqrt(out%prediction_variance)
      out%fitted = series(size(series) - n + 1:) - out%residuals
      if (present(regression_residuals)) then
         out%regression_residuals = regression_residuals
      else
         allocate(out%regression_residuals(0))
      end if
      out%residual_acf = centered_acf(out%residuals, max_lag)
      squared = out%residuals**2
      out%squared_residual_acf = centered_acf(squared, max_lag)
      count_lags = max_lag - first_lag + 1
      allocate(out%lags(count_lags), out%degrees_of_freedom(count_lags))
      allocate(out%ljung_box(count_lags), out%ljung_box_p_value(count_lags))
      allocate(out%squared_ljung_box(count_lags))
      allocate(out%squared_ljung_box_p_value(count_lags))
      statistic = 0.0_dp
      squared_statistic = 0.0_dp
      i = 0
      do h = 1, max_lag
         statistic = statistic + real(n*(n + 2), dp)* &
            out%residual_acf(h + 1)**2/real(n - h, dp)
         squared_statistic = squared_statistic + real(n*(n + 2), dp)* &
            out%squared_residual_acf(h + 1)**2/real(n - h, dp)
         if (h < first_lag) cycle
         i = i + 1
         out%lags(i) = h
         out%degrees_of_freedom(i) = h - parameter_count
         out%ljung_box(i) = statistic
         probability = regularized_gamma_q( &
            0.5_dp*real(out%degrees_of_freedom(i), dp), 0.5_dp*statistic)
         out%ljung_box_p_value(i) = max(0.0_dp, min(1.0_dp, probability))
         out%squared_ljung_box(i) = squared_statistic
         probability = regularized_gamma_q(0.5_dp*real(h, dp), &
            0.5_dp*squared_statistic)
         out%squared_ljung_box_p_value(i) = max(0.0_dp, min(1.0_dp, probability))
      end do
      out%qq_sample = sorted_real(out%standardized_residuals)
      allocate(out%qq_theoretical(n))
      do i = 1, n
         if (n <= 10) then
            probability = (real(i, dp) - 0.375_dp)/(real(n, dp) + 0.25_dp)
         else
            probability = (real(i, dp) - 0.5_dp)/real(n, dp)
         end if
         out%qq_theoretical(i) = normal_quantile(probability)
      end do
   end function diagnostics_from_likelihood

   pure function centered_acf(values, max_lag) result(correlation)
      ! Compute biased sample autocorrelations after centering.
      real(dp), intent(in) :: values(:)
      integer, intent(in) :: max_lag
      real(dp), allocatable :: correlation(:)
      real(dp), allocatable :: centered(:)
      real(dp) :: denominator
      integer :: lag

      allocate(correlation(max_lag + 1))
      centered = values - sum(values)/real(size(values), dp)
      denominator = sum(centered**2)
      if (denominator <= tiny(1.0_dp)) then
         correlation = 0.0_dp
         correlation(1) = 1.0_dp
         return
      end if
      correlation(1) = 1.0_dp
      do lag = 1, max_lag
         correlation(lag + 1) = dot_product(centered(:size(values) - lag), &
            centered(lag + 1:))/denominator
      end do
   end function centered_acf

   pure function sorted_real(values) result(sorted)
      ! Return real values in ascending order using insertion sort.
      real(dp), intent(in) :: values(:)
      real(dp), allocatable :: sorted(:)
      real(dp) :: value
      integer :: i, j

      sorted = values
      do i = 2, size(sorted)
         value = sorted(i)
         j = i - 1
         do while (j >= 1)
            if (sorted(j) <= value) exit
            sorted(j + 1) = sorted(j)
            j = j - 1
         end do
         sorted(j + 1) = value
      end do
   end function sorted_real

   pure function covariance_to_correlation(covariance) result(correlation)
      ! Standardize a covariance matrix while guarding zero variances.
      real(dp), intent(in) :: covariance(:, :)
      real(dp), allocatable :: correlation(:, :)
      real(dp) :: scale
      integer :: i, j, n

      n = size(covariance, 1)
      if (size(covariance, 2) /= n) then
         allocate(correlation(0, 0))
         return
      end if
      allocate(correlation(n, n))
      correlation = 0.0_dp
      do j = 1, n
         do i = 1, n
            scale = covariance(i, i)*covariance(j, j)
            if (scale > 0.0_dp) then
               correlation(i, j) = covariance(i, j)/sqrt(scale)
            end if
         end do
      end do
   end function covariance_to_correlation

   pure function operating_parameter_jacobian(fit) result(jacobian)
      ! Differentiate operating parameters with respect to optimizer coordinates.
      type(arfima_fit_t), intent(in) :: fit
      real(dp), allocatable :: jacobian(:, :)
      type(arfima_model_t) :: lower_model, upper_model
      real(dp), allocatable :: lower(:), upper(:), point(:)
      real(dp) :: step
      integer :: i, n

      n = size(fit%unconstrained_parameters)
      allocate(jacobian(n, n))
      if (n == 0) return
      point = fit%unconstrained_parameters
      do i = 1, n
         step = 1.0e-5_dp*max(1.0_dp, abs(point(i)))
         lower = point
         upper = point
         lower(i) = lower(i) - step
         upper(i) = upper(i) + step
         lower_model = decode_model_parameters(fit%model, lower)
         upper_model = decode_model_parameters(fit%model, upper)
         call impose_fixed_model(fit%model, fit%fixed_parameters, lower_model)
         call impose_fixed_model(fit%model, fit%fixed_parameters, upper_model)
         jacobian(:, i) = (mode_parameters(upper_model, .false.) - &
            mode_parameters(lower_model, .false.))/(2.0_dp*step)
      end do
   end function operating_parameter_jacobian

   pure function coefficient_lag_polynomial(coefficients, period, max_lag) &
      result(polynomial)
      ! Form one minus a coefficient polynomial at a selected lag spacing.
      real(dp), intent(in) :: coefficients(:)
      integer, intent(in) :: period, max_lag
      real(dp), allocatable :: polynomial(:)
      integer :: i

      allocate(polynomial(0:max_lag))
      polynomial = 0.0_dp
      polynomial(0) = 1.0_dp
      do i = 1, min(size(coefficients), max_lag/period)
         polynomial(i*period) = -coefficients(i)
      end do
   end function coefficient_lag_polynomial

   pure function component_acvf(ar, theta, long_memory_type, parameter, &
      max_lag) result(out)
      ! Combine one long-memory covariance with one ARMA covariance.
      real(dp), intent(in) :: ar(:), theta(:), parameter
      integer, intent(in) :: long_memory_type, max_lag
      type(arfima_acvf_t) :: out
      type(arfima_acvf_t) :: long_memory
      type(itsmr_arma_model_t) :: arma_model
      real(dp), allocatable :: arma_covariance(:)

      out%max_lag = max_lag
      if (max_lag < 0 .or. .not. valid_long_memory(long_memory_type, parameter)) then
         out%info = 1
         return
      end if
      select case (long_memory_type)
      case (arfima_long_memory_none)
         allocate(long_memory%covariance(0:max_lag))
         long_memory%covariance = 0.0_dp
         long_memory%covariance(0) = 1.0_dp
      case (arfima_long_memory_fdwn)
         long_memory = arfima_fdwn_acvf(parameter, max_lag)
      case (arfima_long_memory_fgn)
         long_memory = arfima_fgn_acvf(parameter, max_lag)
      case (arfima_long_memory_hd)
         long_memory = arfima_hd_acvf(parameter, max_lag)
      end select
      if (long_memory%info /= 0) then
         out%info = long_memory%info
         return
      end if
      if (size(ar) == 0 .and. size(theta) == 0) then
         out%covariance = long_memory%covariance
         return
      end if

      allocate(arma_model%ar(size(ar)), arma_model%ma(size(theta)))
      arma_model%ar = ar
      arma_model%ma = -theta
      arma_model%innovation_variance = 1.0_dp
      arma_covariance = arma_acvf(arma_model, max_lag)
      if (.not. all(ieee_is_finite(arma_covariance)) .or. &
         maxval(abs(arma_covariance)) >= sqrt(huge(1.0_dp))) then
         out%info = 2
         return
      end if
      allocate(out%covariance(0:max_lag))
      out%covariance = symmetric_covariance_product(long_memory%covariance, &
         arma_covariance, max_lag)
   end function component_acvf

   pure logical function valid_model(model) result(valid)
      ! Check dimensions and scalar domains of an ARFIMA model.
      type(arfima_model_t), intent(in) :: model

      valid = allocated(model%ar) .and. allocated(model%theta) .and. &
         allocated(model%seasonal_ar) .and. allocated(model%seasonal_theta)
      if (.not. valid) return
      valid = all(ieee_is_finite(model%ar)) .and. &
         all(ieee_is_finite(model%theta)) .and. &
         all(ieee_is_finite(model%seasonal_ar)) .and. &
         all(ieee_is_finite(model%seasonal_theta)) .and. &
         stable_coefficients(model%ar) .and. stable_coefficients(model%theta) .and. &
         stable_coefficients(model%seasonal_ar) .and. &
         stable_coefficients(model%seasonal_theta) .and. &
         ieee_is_finite(model%innovation_variance) .and. &
         ieee_is_finite(model%mean) .and. model%innovation_variance > 0.0_dp .and. &
         model%difference_order >= 0 .and. model%seasonal_difference_order >= 0 .and. &
         valid_long_memory(model%long_memory_type, model%long_memory_parameter)
      if (.not. valid) return
      if (model%period == 0) then
         valid = size(model%seasonal_ar) == 0 .and. size(model%seasonal_theta) == 0 .and. &
            model%seasonal_difference_order == 0 .and. &
            model%seasonal_long_memory_parameter == 0.0_dp
      else
         valid = model%period >= 2 .and. &
            valid_long_memory(model%seasonal_long_memory_type, &
            model%seasonal_long_memory_parameter)
      end if
   end function valid_model

   pure logical function stable_coefficients(coefficients) result(stable)
      ! Test an AR-form polynomial through its partial autocorrelations.
      real(dp), intent(in) :: coefficients(:)
      real(dp) :: partial(size(coefficients))

      if (size(coefficients) == 0) then
         stable = .true.
         return
      end if
      partial = ar_to_partial(coefficients)
      stable = all(ieee_is_finite(partial)) .and. all(abs(partial) < 1.0_dp)
   end function stable_coefficients

   pure logical function valid_long_memory(long_memory_type, parameter) result(valid)
      ! Check one generalized long-memory parameter.
      integer, intent(in) :: long_memory_type
      real(dp), intent(in) :: parameter

      valid = ieee_is_finite(parameter)
      if (.not. valid) return
      select case (long_memory_type)
      case (arfima_long_memory_none)
         valid = parameter == 0.0_dp
      case (arfima_long_memory_fdwn)
         valid = parameter > -1.0_dp .and. parameter < 0.5_dp
      case (arfima_long_memory_fgn)
         valid = parameter > 0.0_dp .and. parameter < 1.0_dp
      case (arfima_long_memory_hd)
         valid = parameter > 0.0_dp .and. parameter < 3.0_dp
      case default
         valid = .false.
      end select
   end function valid_long_memory

   pure integer function next_power_of_two(value) result(power)
      ! Return the smallest power of two not below a nonnegative integer.
      integer, intent(in) :: value

      power = 1
      do while (power < value)
         power = 2*power
      end do
   end function next_power_of_two

   pure function polynomial_product_truncated(first, second, max_lag) result(product)
      ! Multiply two lag polynomials through max_lag.
      real(dp), intent(in) :: first(0:), second(0:)
      integer, intent(in) :: max_lag
      real(dp) :: product(0:max_lag)
      integer :: i, j

      product = 0.0_dp
      do i = 0, min(max_lag, ubound(first, 1))
         do j = 0, min(max_lag - i, ubound(second, 1))
            product(i + j) = product(i + j) + first(i)*second(j)
         end do
      end do
   end function polynomial_product_truncated

   pure function symmetric_covariance_product(first, second, max_lag) result(product)
      ! Reproduce arfima's circular convolution of symmetric covariances.
      real(dp), intent(in) :: first(0:), second(0:)
      integer, intent(in) :: max_lag
      real(dp) :: product(0:max_lag)
      integer :: i, lag, limit, sequence_length, target, second_index

      product = 0.0_dp
      limit = min(ubound(first, 1), ubound(second, 1))
      if (limit == 0) then
         product(0) = first(0)*second(0)
         return
      end if
      sequence_length = 2*limit
      do lag = 0, max_lag
         target = 2*limit - 2 - lag
         do i = 0, sequence_length - 1
            second_index = modulo(target - i, sequence_length)
            product(lag) = product(lag) + &
               first(abs(i - limit + 1))*second(abs(second_index - limit + 1))
         end do
      end do
   end function symmetric_covariance_product

   pure real(dp) function borwein_zeta(value, terms) result(zeta_value)
      ! Approximate the Riemann zeta function using Borwein's finite sum.
      real(dp), intent(in) :: value
      integer, intent(in) :: terms
      real(dp), allocatable :: coefficient(:)
      real(dp) :: term
      integer :: k

      allocate(coefficient(0:terms))
      coefficient(0) = 1.0_dp
      do k = 1, terms
         term = exp(log_gamma(real(terms + k, dp)) - &
            log_gamma(real(terms - k + 1, dp)) - log_gamma(real(2*k + 1, dp)))
         coefficient(k) = coefficient(k - 1) + real(terms, dp)*term*4.0_dp**k
      end do
      zeta_value = 0.0_dp
      do k = 0, terms - 1
         if (mod(k, 2) == 0) then
            zeta_value = zeta_value + (coefficient(k) - coefficient(terms))/ &
               real(k + 1, dp)**value
         else
            zeta_value = zeta_value - (coefficient(k) - coefficient(terms))/ &
               real(k + 1, dp)**value
         end if
      end do
      zeta_value = -zeta_value/(coefficient(terms)*(1.0_dp - 2.0_dp**(1.0_dp - value)))
   end function borwein_zeta

end module arfima_mod
