! SPDX-License-Identifier: GPL-2.0-or-later
! SPDX-FileComment: Algorithms translated from the R tfarima package.
! Transfer-function and lag-polynomial algorithms translated from tfarima.
module tfarima_mod
   use kind_mod, only: dp
   use arima2_mod, only: arima2_roots_t, arma_polynomial_roots
   use forecast_mod, only: acf_values, pacf_values, ccf_values
   use itsmr_mod, only: itsmr_arma_model_t, arma_acvf
   use kfas_mod, only: ssm_model_t, kfs_filter_t, kfs_smoother_t, &
      kfs_filter, kfs_filter_diffuse, kfs_smooth
   use time_series_calendar_mod, only: date_t, date_valid, date_days_in_month, &
      date_day_of_week, date_day_number, date_easter, operator(+), operator(-)
   use time_series_linalg_mod, only: invert_matrix, cholesky_lower
   use time_series_optimization_mod, only: optimization_result_t, &
      bfgs_minimize_fd, nelder_mead_minimize, finite_difference_hessian
   use time_series_stats_mod, only: normal_quantile
   use time_series_diagnostics_mod, only: weighted_box_test_t, &
      weighted_box_test, box_test_ljung_box, residual_raw
   use time_series_random_mod, only: set_random_seed, random_standard_normal_matrix
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   private

   integer, parameter, public :: tfarima_intervention_pulse = 1
   integer, parameter, public :: tfarima_intervention_step = 2
   integer, parameter, public :: tfarima_intervention_ramp = 3
   integer, parameter, public :: tfarima_outlier_io = 1
   integer, parameter, public :: tfarima_outlier_ao = 2
   integer, parameter, public :: tfarima_outlier_ls = 3
   integer, parameter, public :: tfarima_outlier_tc = 4
   integer, parameter, public :: tfarima_calendar_dif = 1
   integer, parameter, public :: tfarima_calendar_td = 2
   integer, parameter, public :: tfarima_calendar_td7 = 3
   integer, parameter, public :: tfarima_calendar_td6 = 4
   integer, parameter, public :: tfarima_calendar_wd = 5
   integer, parameter, public :: tfarima_calendar_null = 6
   integer, parameter, public :: tfarima_selection_regressor = 1
   integer, parameter, public :: tfarima_selection_input = 2
   integer, parameter, public :: tfarima_factor_newton = 1
   integer, parameter, public :: tfarima_factor_roots = 2
   integer, parameter, public :: tfarima_factor_bauer = 3
   integer, parameter, public :: tfarima_factor_laurie = 4
   integer, parameter, public :: tfarima_factor_wilson = 5
   integer, parameter, public :: tfarima_factor_best = 6
   integer, parameter, public :: tfarima_decomposition_forecast = 1
   integer, parameter, public :: tfarima_decomposition_backcast = 2
   integer, parameter, public :: tfarima_decomposition_mixed = 3
   integer, parameter, public :: tfarima_component_trend = 1
   integer, parameter, public :: tfarima_component_seasonal = 2
   integer, parameter, public :: tfarima_component_exponential = 3
   integer, parameter, public :: tfarima_component_cycle = 4

   type, public :: tfarima_polynomial_division_t
      ! Quotient and remainder from lag-polynomial long division.
      real(dp), allocatable :: quotient(:), remainder(:)
      integer :: info = 0
   end type tfarima_polynomial_division_t

   type, public :: tfarima_transfer_t
      ! Rational transfer function omega(B) B^delay / delta(B).
      real(dp), allocatable :: numerator(:), denominator(:)
      integer :: delay = 0
      integer :: info = 0
   end type tfarima_transfer_t

   type, public :: tfarima_forecast_t
      ! ARIMA point forecasts, variances, and normal interval bounds.
      real(dp), allocatable :: mean(:), variance(:), lower(:, :), upper(:, :)
      integer :: info = 0
   end type tfarima_forecast_t

   type, public :: tfarima_ma_factor_t
      ! Finite MA coefficients whose products reproduce supplied autocovariances.
      real(dp), allocatable :: coefficients(:)
      real(dp) :: residual_norm = huge(1.0_dp)
      integer :: iterations = 0
      integer :: method = tfarima_factor_newton
      integer :: info = 0
      logical :: converged = .false.
   end type tfarima_ma_factor_t

   type, public :: tfarima_lag_polynomial_t
      ! Restricted sparse lag polynomial and its expanded coefficients.
      integer, allocatable :: lags(:)
      real(dp), allocatable :: offset(:), loading(:, :), parameters(:)
      real(dp), allocatable :: base(:), expanded(:)
      integer :: exponent = 1
      integer :: info = 0
   end type tfarima_lag_polynomial_t

   type, public :: tfarima_transfer_fit_t
      ! Conditional Gaussian transfer-function estimates and inference.
      type(tfarima_transfer_t) :: transfer
      real(dp), allocatable :: parameters(:), covariance(:, :), standard_error(:)
      real(dp), allocatable :: fitted(:), residuals(:)
      real(dp) :: mean = 0.0_dp
      real(dp) :: innovation_variance = 0.0_dp
      real(dp) :: log_likelihood = -huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      real(dp) :: bic = huge(1.0_dp)
      integer :: iterations = 0
      integer :: info = 0
      logical :: converged = .false.
   end type tfarima_transfer_fit_t

   type, public :: tfarima_transfer_spec_t
      ! Rational dynamic input and masks selecting its estimated coefficients.
      type(tfarima_transfer_t) :: transfer
      logical, allocatable :: estimate_numerator(:), estimate_denominator(:)
      integer :: info = 0
   end type tfarima_transfer_spec_t

   type, public :: tfarima_prewhitened_ccf_t
      ! Prewhitened input-output correlations and scaled impulse weights.
      real(dp), allocatable :: correlation(:), impulse_response(:)
      real(dp), allocatable :: input_residuals(:), output_residuals(:)
      integer, allocatable :: lag(:)
      logical, allocatable :: significant(:)
      real(dp) :: critical_value = 0.0_dp
      integer :: info = 0
   end type tfarima_prewhitened_ccf_t

   type, public :: tfarima_transfer_identification_t
      ! Identified delay and rational transfer starting specification.
      type(tfarima_transfer_spec_t) :: specification
      type(tfarima_prewhitened_ccf_t) :: diagnostic
      integer :: selected_delay = 0
      integer :: info = 0
   end type tfarima_transfer_identification_t

   type, public :: tfarima_outlier_result_t
      ! Detected outlier positions, types, effects, t-ratios, and TC decay.
      integer, allocatable :: position(:), outlier_type(:)
      real(dp), allocatable :: effect(:), t_ratio(:), decay(:)
      integer :: passes = 0
      integer :: info = 0
   end type tfarima_outlier_result_t

   type, public :: tfarima_calendar_result_t
      ! Monthly calendar regressors and selected coding information.
      real(dp), allocatable :: values(:, :)
      integer :: form = tfarima_calendar_null
      integer :: reference_weekday = 0
      integer :: info = 0
   end type tfarima_calendar_result_t

   type, public :: tfarima_ucarima_component_t
      ! Independent ARIMA component in normalized lag-polynomial form.
      real(dp), allocatable :: ar_polynomial(:), difference_polynomial(:)
      real(dp), allocatable :: ma_polynomial(:)
      real(dp) :: innovation_variance = 0.0_dp
      integer :: info = 0
   end type tfarima_ucarima_component_t

   type, public :: tfarima_exact_transfer_fit_t
      ! Exact multi-input transfer estimates, noise model, and inference.
      type(tfarima_transfer_spec_t), allocatable :: inputs(:)
      type(tfarima_ucarima_component_t) :: noise
      real(dp), allocatable :: parameters(:), covariance(:, :), standard_error(:)
      real(dp), allocatable :: regression_coefficients(:)
      real(dp), allocatable :: signal(:, :), regression_signal(:)
      real(dp), allocatable :: fitted(:), residuals(:), standardized_residuals(:)
      real(dp) :: log_likelihood = -huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      real(dp) :: bic = huge(1.0_dp)
      integer :: iterations = 0
      integer :: info = 0
      logical :: converged = .false.
   end type tfarima_exact_transfer_fit_t

   type, public :: tfarima_exact_transfer_forecast_t
      ! Known-input transfer forecasts and Gaussian noise uncertainty.
      real(dp), allocatable :: mean(:), variance(:), lower(:, :), upper(:, :)
      integer :: info = 0
   end type tfarima_exact_transfer_forecast_t

   type, public :: tfarima_transfer_selection_t
      ! Backward-selected exact transfer model and original-variable history.
      type(tfarima_exact_transfer_fit_t) :: fit
      logical, allocatable :: retained_regressor(:), retained_input(:)
      real(dp), allocatable :: regressor_p_value(:), input_p_value(:)
      integer, allocatable :: removed_kind(:), removed_index(:)
      integer :: steps = 0
      integer :: info = 0
      logical :: converged = .false.
   end type tfarima_transfer_selection_t

   type, public :: tfarima_transfer_diagnostics_t
      ! Numerical residual, portmanteau, spectrum, and input diagnostics.
      real(dp), allocatable :: residual_acf(:), residual_pacf(:)
      real(dp), allocatable :: cumulative_periodogram(:)
      type(weighted_box_test_t) :: ljung_box, weighted_ljung_box
      type(tfarima_prewhitened_ccf_t), allocatable :: input_ccf(:)
      real(dp) :: residual_mean = 0.0_dp
      real(dp) :: residual_variance = 0.0_dp
      real(dp) :: skewness = 0.0_dp
      real(dp) :: excess_kurtosis = 0.0_dp
      real(dp) :: normality_statistic = 0.0_dp
      real(dp) :: normality_p_value = 1.0_dp
      real(dp) :: cumulative_periodogram_statistic = 0.0_dp
      real(dp) :: cumulative_periodogram_critical = 0.0_dp
      real(dp) :: cumulative_periodogram_p_value = 1.0_dp
      integer :: info = 0
   end type tfarima_transfer_diagnostics_t

   type, public :: tfarima_transfer_signal_t
      ! Per-input, regression, and aggregate deterministic transfer signals.
      real(dp), allocatable :: input_signal(:, :), regression_signal(:), total(:)
      integer :: info = 0
   end type tfarima_transfer_signal_t

   type, public :: tfarima_transfer_simulation_t
      ! Simulated outputs and their retained input, signal, and noise paths.
      real(dp), allocatable :: output(:, :), noise(:, :), innovation(:, :)
      real(dp), allocatable :: input(:, :, :), input_signal(:, :, :)
      real(dp), allocatable :: input_innovation(:, :, :)
      real(dp), allocatable :: regression_signal(:, :)
      integer :: observations = 0
      integer :: simulations = 0
      integer :: burn_in = 0
      integer :: info = 0
   end type tfarima_transfer_simulation_t

   type, public :: tfarima_arima_algebra_t
      ! Combined ARIMA model, spectral factor, and reconstruction diagnostics.
      type(tfarima_ucarima_component_t) :: model
      type(tfarima_ma_factor_t) :: factor
      real(dp), allocatable :: target_autocovariance(:)
      real(dp), allocatable :: reconstructed_autocovariance(:)
      real(dp) :: reconstruction_error = huge(1.0_dp)
      integer :: info = 0
   end type tfarima_arima_algebra_t

   type, public :: tfarima_weights_t
      ! Finite response weights, cumulative weights, and delta-method inference.
      real(dp), allocatable :: weight(:), cumulative(:)
      real(dp), allocatable :: covariance(:, :), standard_error(:)
      real(dp), allocatable :: cumulative_covariance(:, :)
      real(dp), allocatable :: cumulative_standard_error(:)
      integer :: info = 0
   end type tfarima_weights_t

   type, public :: tfarima_leverrier_t
      ! Characteristic polynomial and observation-weighted adjoint coefficients.
      real(dp), allocatable :: characteristic(:), adjoint_loading(:, :)
      integer :: info = 0
   end type tfarima_leverrier_t

   type, public :: tfarima_ssm_form_t
      ! Time-invariant univariate SSM with joint observation and state noise.
      real(dp), allocatable :: observation_loading(:), transition(:, :)
      real(dp), allocatable :: disturbance_covariance(:, :)
      logical :: state_noise_contemporaneous = .false.
      integer :: info = 0
   end type tfarima_ssm_form_t

   type, public :: tfarima_ssm_reduction_t
      ! Reduced ARIMA spectrum and state-disturbance polynomial diagnostics.
      type(tfarima_ucarima_component_t) :: model
      type(tfarima_ma_factor_t) :: factor
      real(dp), allocatable :: characteristic(:)
      real(dp), allocatable :: disturbance_numerator(:, :)
      real(dp), allocatable :: disturbance_covariance(:, :)
      real(dp), allocatable :: target_autocovariance(:)
      real(dp), allocatable :: reconstructed_autocovariance(:)
      real(dp) :: reconstruction_error = huge(1.0_dp)
      integer :: info = 0
   end type tfarima_ssm_reduction_t

   type, public :: tfarima_decomposition_basis_t
      ! Root table, component classification, and deterministic root basis.
      real(dp), allocatable :: root_table(:, :), classification(:, :)
      real(dp), allocatable :: basis(:, :)
      integer :: info = 0
   end type tfarima_decomposition_basis_t

   type, public :: tfarima_root_decomposition_t
      ! Root-based component effects and seasonal-adjustment diagnostics.
      type(tfarima_decomposition_basis_t) :: decomposition_basis
      real(dp), allocatable :: state_coefficients(:, :), innovations(:)
      real(dp), allocatable :: component(:, :), irregular(:), reconstruction(:)
      real(dp), allocatable :: seasonally_adjusted(:)
      integer :: method = tfarima_decomposition_mixed
      integer :: info = 0
   end type tfarima_root_decomposition_t

   type, public :: tfarima_structural_ssm_t
      ! Direct ARIMA structural form and disturbance-variance matching diagnostics.
      type(tfarima_ssm_form_t) :: form
      type(tfarima_decomposition_basis_t) :: decomposition_basis
      type(tfarima_leverrier_t) :: leverrier
      real(dp), allocatable :: target_autocovariance(:), variance_design(:, :)
      real(dp), allocatable :: disturbance_variance(:), fitted_autocovariance(:)
      real(dp) :: reconstruction_error = huge(1.0_dp)
      integer :: info = 0
      logical :: multiple_sources = .true.
      logical :: admissible = .false.
      logical :: used_nonnegative_fit = .false.
   end type tfarima_structural_ssm_t

   type, public :: tfarima_structural_initialization_t
      ! GLS initial state and covariance for a joint-disturbance structural model.
      real(dp), allocatable :: state(:), covariance(:, :)
      integer :: observations = 0
      integer :: info = 0
   end type tfarima_structural_initialization_t

   type, public :: tfarima_structural_filter_t
      ! Joint-disturbance innovations and observation/next-state filter moments.
      real(dp), allocatable :: observation_state(:, :)
      real(dp), allocatable :: observation_covariance(:, :, :)
      real(dp), allocatable :: filtered_state(:, :), filtered_covariance(:, :, :)
      real(dp), allocatable :: gain(:, :), innovation(:), innovation_variance(:)
      real(dp), allocatable :: standardized_innovation(:)
      real(dp) :: log_likelihood = -huge(1.0_dp)
      integer :: info = 0
   end type tfarima_structural_filter_t

   type, public :: tfarima_structural_smoother_t
      ! Smoothed observation-aligned states and conditional covariance matrices.
      type(tfarima_structural_filter_t) :: filter
      real(dp), allocatable :: state(:, :), covariance(:, :, :)
      integer :: info = 0
   end type tfarima_structural_smoother_t

   type, public :: tfarima_structural_forecast_t
      ! Structural-model forecasts on latent and optional original scales.
      real(dp), allocatable :: latent_mean(:), latent_variance(:)
      real(dp), allocatable :: mean(:), variance(:)
      real(dp), allocatable :: state(:, :), state_covariance(:, :, :)
      integer :: info = 0
   end type tfarima_structural_forecast_t

   type, public :: tfarima_band_cholesky_t
      ! Lower band Cholesky factor and its logarithmic determinant.
      real(dp), allocatable :: factor(:, :)
      real(dp) :: log_determinant = 0.0_dp
      integer :: info = 0
   end type tfarima_band_cholesky_t

   type, public :: tfarima_reduced_likelihood_t
      ! Profiled exact likelihood and reduced-form residual diagnostics.
      real(dp), allocatable :: residuals(:), standardized_residuals(:)
      real(dp), allocatable :: covariance_band(:, :), cholesky_band(:, :)
      real(dp) :: innovation_variance = 0.0_dp
      real(dp) :: log_likelihood = -huge(1.0_dp)
      real(dp) :: log_determinant = 0.0_dp
      integer :: bandwidth = 0
      integer :: info = 0
   end type tfarima_reduced_likelihood_t

   type, public :: tfarima_ucarima_model_t
      ! UCARIMA components and their Cramer-Wold aggregate representation.
      type(tfarima_ucarima_component_t), allocatable :: components(:)
      real(dp), allocatable :: denominator(:), ma_polynomial(:)
      real(dp), allocatable :: numerator_autocovariance(:)
      real(dp), allocatable :: component_numerator(:, :)
      real(dp) :: innovation_variance = 0.0_dp
      real(dp) :: factor_residual_norm = huge(1.0_dp)
      integer :: factor_iterations = 0
      integer :: factor_method = tfarima_factor_best
      integer :: info = 0
      logical :: factor_converged = .false.
   end type tfarima_ucarima_model_t

   type, public :: tfarima_wk_filter_t
      ! Rational Wiener-Kolmogorov filter and symmetric finite weights.
      real(dp), allocatable :: symmetric_numerator(:), denominator(:), weights(:)
      real(dp) :: expansion_tail_energy = 0.0_dp
      integer :: component = 0
      integer :: max_lag = 0
      integer :: info = 0
   end type tfarima_wk_filter_t

   type, public :: tfarima_ucarima_decomposition_t
      ! Wiener-Kolmogorov component estimates and reconstruction residual.
      real(dp), allocatable :: component(:, :), reconstruction(:), remainder(:)
      integer :: info = 0
   end type tfarima_ucarima_decomposition_t

   type, public :: tfarima_bezout_t
      ! Polynomial GCD and coefficients satisfying first*u + second*v = gcd.
      real(dp), allocatable :: u(:), v(:), gcd(:)
      integer :: info = 0
   end type tfarima_bezout_t

   type, public :: tfarima_partial_fraction_t
      ! Numerators over supplied coprime factors and reconstruction diagnostics.
      real(dp), allocatable :: numerator(:, :), reconstruction(:)
      integer, allocatable :: numerator_length(:)
      real(dp) :: reconstruction_error = 0.0_dp
      integer :: info = 0
   end type tfarima_partial_fraction_t

   type, public :: tfarima_ucarima_conversion_t
      ! Aggregate-ARIMA to UCARIMA conversion and spectral diagnostics.
      type(tfarima_ucarima_model_t) :: model
      type(tfarima_partial_fraction_t) :: partial_fraction
      real(dp), allocatable :: wold_numerator(:), wold_denominator(:)
      real(dp), allocatable :: quotient(:), remainder(:)
      real(dp) :: reconstruction_error = huge(1.0_dp)
      integer :: info = 0
      logical :: admissible = .false.
   end type tfarima_ucarima_conversion_t

   type, public :: tfarima_ucarima_ssm_t
      ! Block innovations state-space form and component observation loadings.
      type(ssm_model_t) :: model
      real(dp), allocatable :: component_loading(:, :)
      integer, allocatable :: state_start(:), state_size(:)
      integer :: info = 0
      logical :: diffuse = .false.
   end type tfarima_ucarima_ssm_t

   type, public :: tfarima_ucarima_smoothing_t
      ! Smoothed components, uncertainty, likelihood, and component forecasts.
      real(dp), allocatable :: component(:, :), variance(:, :)
      real(dp), allocatable :: forecast(:, :), forecast_variance(:, :)
      real(dp), allocatable :: lower(:, :, :), upper(:, :, :)
      real(dp) :: log_likelihood = 0.0_dp
      integer :: info = 0
   end type tfarima_ucarima_smoothing_t

   type, public :: tfarima_ucarima_fit_t
      ! Exact KFAS estimates, inference, residuals, and fitted UCARIMA model.
      type(tfarima_ucarima_model_t) :: model
      real(dp), allocatable :: parameters(:), covariance(:, :), standard_error(:)
      real(dp), allocatable :: regression_coefficients(:)
      real(dp), allocatable :: residuals(:), standardized_residuals(:)
      real(dp) :: log_likelihood = -huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      real(dp) :: bic = huge(1.0_dp)
      integer :: iterations = 0
      integer :: info = 0
      logical :: converged = .false.
   end type tfarima_ucarima_fit_t

   public :: tfarima_polynomial_multiply, tfarima_polynomial_divide
   public :: tfarima_polynomial_gcd, tfarima_polynomial_power
   public :: tfarima_polynomial_ratio, tfarima_polynomial_derivative
   public :: tfarima_transfer, tfarima_filter, tfarima_impulse_response
   public :: tfarima_difference, tfarima_arima_forecast, tfarima_arima_backcast
   public :: tfarima_intervention, tfarima_seasonal_dummies
   public :: tfarima_harmonic_regressors, tfarima_standardize
   public :: tfarima_autocovariance, tfarima_partial_autocorrelation
   public :: tfarima_autocovariance_to_ma
   public :: tfarima_cramer_wold_factor
   public :: tfarima_palindromic_to_wold, tfarima_wold_to_palindromic
   public :: tfarima_lag_polynomial, tfarima_update_lag_polynomial
   public :: tfarima_polynomial_admissible, tfarima_transfer_fit
   public :: tfarima_restricted_transfer_fit
   public :: tfarima_transfer_spec, tfarima_exact_transfer_fit
   public :: tfarima_exact_transfer_forecast
   public :: tfarima_select_transfer
   public :: tfarima_prewhitened_ccf, tfarima_identify_transfer
   public :: tfarima_exact_transfer_ccf
   public :: tfarima_diagnose_transfer
   public :: tfarima_simulate_transfer, tfarima_simulate_exact_model
   public :: tfarima_combine_arima, tfarima_psi_weights, tfarima_pi_weights
   public :: tfarima_leverrier_faddeev, tfarima_ssm_form
   public :: tfarima_switch_ssm_form, tfarima_reduce_ssm, tfarima_ssm_to_arima
   public :: tfarima_root_decomposition_basis, tfarima_root_decompose
   public :: tfarima_arima_to_structural_ssm
   public :: tfarima_structural_initialize, tfarima_structural_filter
   public :: tfarima_structural_smooth, tfarima_structural_forecast
   public :: tfarima_band_cholesky, tfarima_band_forward_solve
   public :: tfarima_reduced_likelihood
   public :: tfarima_ucarima_reduced_likelihood
   public :: tfarima_structural_reduced_likelihood
   public :: tfarima_outlier_response, tfarima_detect_outliers
   public :: tfarima_month_lengths, tfarima_weekday_counts
   public :: tfarima_leap_year_regressor, tfarima_easter_regressor
   public :: tfarima_calendar_regressors
   public :: tfarima_ucarima_component, tfarima_build_ucarima
   public :: tfarima_wiener_kolmogorov_filter, tfarima_ucarima_decompose
   public :: tfarima_extended_polynomial_gcd, tfarima_partial_fractions
   public :: tfarima_arima_to_ucarima
   public :: tfarima_ucarima_state_space, tfarima_ucarima_smooth
   public :: tfarima_ucarima_fit

contains

   pure function tfarima_ucarima_component(ar_polynomial, &
      difference_polynomial, ma_polynomial, innovation_variance) result(out)
      ! Construct one normalized independent UCARIMA component.
      real(dp), intent(in) :: ar_polynomial(:), difference_polynomial(:)
      real(dp), intent(in) :: ma_polynomial(:), innovation_variance
      type(tfarima_ucarima_component_t) :: out
      real(dp) :: scale

      if (size(ar_polynomial) < 1 .or. size(difference_polynomial) < 1 .or. &
         size(ma_polynomial) < 1 .or. innovation_variance < 0.0_dp .or. &
         .not. ieee_is_finite(innovation_variance) .or. &
         .not. all(ieee_is_finite(ar_polynomial)) .or. &
         .not. all(ieee_is_finite(difference_polynomial)) .or. &
         .not. all(ieee_is_finite(ma_polynomial)) .or. &
         abs(ar_polynomial(1)) <= tiny(1.0_dp) .or. &
         abs(difference_polynomial(1)) <= tiny(1.0_dp) .or. &
         abs(ma_polynomial(1)) <= tiny(1.0_dp)) then
         out%info = 1
         return
      end if
      if (.not. tfarima_polynomial_admissible(ar_polynomial)) then
         out%info = 2
         return
      end if
      scale = ma_polynomial(1)/(ar_polynomial(1)*difference_polynomial(1))
      out%ar_polynomial = ar_polynomial/ar_polynomial(1)
      out%difference_polynomial = &
         difference_polynomial/difference_polynomial(1)
      out%ma_polynomial = ma_polynomial/ma_polynomial(1)
      out%innovation_variance = innovation_variance*scale**2
   end function tfarima_ucarima_component

   pure function tfarima_build_ucarima(components, tolerance, max_iterations) &
      result(out)
      ! Aggregate independent components by Cramer-Wold spectral factorization.
      type(tfarima_ucarima_component_t), intent(in) :: components(:)
      real(dp), intent(in), optional :: tolerance
      integer, intent(in), optional :: max_iterations
      type(tfarima_ucarima_model_t) :: out
      type(tfarima_ma_factor_t) :: factor
      real(dp), allocatable :: component_denominator(:), lifted(:), gcd(:)
      real(dp), allocatable :: covariance(:)
      real(dp) :: tol
      integer :: limit, component, other, lag, max_length

      tol = 1.0e-10_dp
      if (present(tolerance)) tol = tolerance
      limit = 500
      if (present(max_iterations)) limit = max_iterations
      if (size(components) < 1 .or. tol <= 0.0_dp .or. limit < 1) then
         out%info = 1
         return
      end if
      allocate(out%components(size(components)))
      do component = 1, size(components)
         if (.not. allocated(components(component)%ar_polynomial) .or. &
            .not. allocated(components(component)%difference_polynomial) .or. &
            .not. allocated(components(component)%ma_polynomial)) then
            out%info = 1
            return
         end if
         out%components(component) = tfarima_ucarima_component( &
            components(component)%ar_polynomial, &
            components(component)%difference_polynomial, &
            components(component)%ma_polynomial, &
            components(component)%innovation_variance)
         if (components(component)%info /= 0 .or. &
            out%components(component)%info /= 0) then
            out%info = 1
            return
         end if
      end do

      do component = 1, size(components) - 1
         component_denominator = tfarima_polynomial_multiply( &
            out%components(component)%ar_polynomial, &
            out%components(component)%difference_polynomial)
         do other = component + 1, size(components)
            lifted = tfarima_polynomial_multiply( &
               out%components(other)%ar_polynomial, &
               out%components(other)%difference_polynomial)
            gcd = tfarima_polynomial_gcd(component_denominator, lifted, tol)
            if (size(gcd) < 1) then
               out%info = 1
               return
            end if
            if (size(gcd) > 1) then
               out%info = 2
               return
            end if
         end do
      end do

      out%denominator = [1.0_dp]
      do component = 1, size(components)
         component_denominator = tfarima_polynomial_multiply( &
            out%components(component)%ar_polynomial, &
            out%components(component)%difference_polynomial)
         out%denominator = tfarima_polynomial_multiply(out%denominator, &
            component_denominator)
      end do
      max_length = 1
      do component = 1, size(components)
         lifted = out%components(component)%ma_polynomial
         do other = 1, size(components)
            if (other == component) cycle
            component_denominator = tfarima_polynomial_multiply( &
               out%components(other)%ar_polynomial, &
               out%components(other)%difference_polynomial)
            lifted = tfarima_polynomial_multiply(lifted, component_denominator)
         end do
         max_length = max(max_length, size(lifted))
      end do
      allocate(out%component_numerator(max_length, size(components)))
      out%component_numerator = 0.0_dp
      do component = 1, size(components)
         lifted = out%components(component)%ma_polynomial
         do other = 1, size(components)
            if (other == component) cycle
            component_denominator = tfarima_polynomial_multiply( &
               out%components(other)%ar_polynomial, &
               out%components(other)%difference_polynomial)
            lifted = tfarima_polynomial_multiply(lifted, component_denominator)
         end do
         out%component_numerator(:size(lifted), component) = lifted
      end do

      allocate(covariance(max_length))
      covariance = 0.0_dp
      do lag = 0, max_length - 1
         do component = 1, size(components)
            covariance(lag + 1) = covariance(lag + 1) + &
               out%components(component)%innovation_variance*dot_product( &
               out%component_numerator(:max_length - lag, component), &
               out%component_numerator(lag + 1:, component))
         end do
      end do
      covariance = trim_polynomial(covariance, tol)
      out%numerator_autocovariance = covariance
      factor = tfarima_cramer_wold_factor(covariance, tfarima_factor_best, &
         tolerance=tol, max_iterations=limit)
      out%factor_iterations = factor%iterations
      out%factor_method = factor%method
      out%factor_residual_norm = factor%residual_norm
      out%factor_converged = factor%converged
      if (factor%info /= 0 .or. .not. allocated(factor%coefficients) .or. &
         abs(factor%coefficients(1)) <= tiny(1.0_dp)) then
         out%info = 3
         return
      end if
      out%innovation_variance = factor%coefficients(1)**2
      out%ma_polynomial = factor%coefficients/factor%coefficients(1)
      out%info = 0
   end function tfarima_build_ucarima

   pure function tfarima_wiener_kolmogorov_filter(model, component, max_lag, &
      tolerance) result(out)
      ! Construct one component's symmetric Wiener-Kolmogorov filter.
      type(tfarima_ucarima_model_t), intent(in) :: model
      integer, intent(in) :: component, max_lag
      real(dp), intent(in), optional :: tolerance
      type(tfarima_wk_filter_t) :: out
      real(dp), allocatable :: lifted(:), impulse(:)
      real(dp) :: scale, tol
      integer :: expansion_order, lag, length, tail_start

      tol = 1.0e-10_dp
      if (present(tolerance)) tol = tolerance
      if (model%info /= 0 .or. component < 1 .or. &
         .not. allocated(model%components) .or. &
         component > size(model%components) .or. max_lag < 0 .or. &
         tol <= 0.0_dp .or. model%innovation_variance <= 0.0_dp .or. &
         .not. allocated(model%ma_polynomial) .or. &
         .not. allocated(model%component_numerator)) then
         out%info = 1
         return
      end if
      out%component = component
      out%max_lag = max_lag
      lifted = trim_polynomial(model%component_numerator(:, component), tol)
      scale = model%components(component)%innovation_variance/ &
         model%innovation_variance
      allocate(out%symmetric_numerator(size(lifted)))
      do lag = 0, size(lifted) - 1
         out%symmetric_numerator(lag + 1) = scale*dot_product( &
            lifted(:size(lifted) - lag), lifted(lag + 1:))
      end do
      out%denominator = model%ma_polynomial

      expansion_order = max_lag + 4096
      impulse = tfarima_polynomial_ratio(lifted, model%ma_polynomial, &
         expansion_order)
      if (size(impulse) /= expansion_order + 1 .or. &
         .not. all(ieee_is_finite(impulse))) then
         out%info = 2
         return
      end if
      length = size(impulse)
      allocate(out%weights(2*max_lag + 1))
      do lag = 0, max_lag
         out%weights(max_lag + 1 + lag) = scale*dot_product( &
            impulse(:length - lag), impulse(lag + 1:))
         out%weights(max_lag + 1 - lag) = &
            out%weights(max_lag + 1 + lag)
      end do
      tail_start = max(1, length - 63)
      out%expansion_tail_energy = scale*sum(impulse(tail_start:)**2)
      if (out%expansion_tail_energy > sqrt(tol)) out%info = 3
   end function tfarima_wiener_kolmogorov_filter

   pure function tfarima_ucarima_decompose(series, model, max_lag, tolerance) &
      result(out)
      ! Extract all UCARIMA components with endpoint-extended WK filters.
      real(dp), intent(in) :: series(:)
      type(tfarima_ucarima_model_t), intent(in) :: model
      integer, intent(in), optional :: max_lag
      real(dp), intent(in), optional :: tolerance
      type(tfarima_ucarima_decomposition_t) :: out
      type(tfarima_wk_filter_t) :: filter
      type(tfarima_forecast_t) :: forecast
      real(dp), allocatable :: ar_polynomial(:), difference_polynomial(:)
      real(dp), allocatable :: backcast(:), extended(:)
      real(dp) :: tol
      integer :: selected_lag, component, observation, offset, index, n

      selected_lag = 50
      if (present(max_lag)) selected_lag = max_lag
      tol = 1.0e-10_dp
      if (present(tolerance)) tol = tolerance
      if (size(series) < 2 .or. model%info /= 0 .or. &
         .not. all(ieee_is_finite(series)) .or. &
         .not. allocated(model%components) .or. selected_lag < 0 .or. &
         tol <= 0.0_dp) then
         out%info = 1
         return
      end if
      ar_polynomial = [1.0_dp]
      difference_polynomial = [1.0_dp]
      do component = 1, size(model%components)
         ar_polynomial = tfarima_polynomial_multiply(ar_polynomial, &
            model%components(component)%ar_polynomial)
         difference_polynomial = tfarima_polynomial_multiply( &
            difference_polynomial, &
            model%components(component)%difference_polynomial)
      end do
      n = size(series)
      if (selected_lag > 0) then
         forecast = tfarima_arima_forecast(series, ar_polynomial, &
            difference_polynomial, model%ma_polynomial, &
            model%innovation_variance, selected_lag)
         backcast = tfarima_arima_backcast(series, ar_polynomial, &
            difference_polynomial, model%ma_polynomial, &
            model%innovation_variance, selected_lag)
         if (forecast%info /= 0 .or. size(backcast) /= selected_lag) then
            out%info = 2
            return
         end if
         extended = [backcast, series, forecast%mean]
      else
         extended = series
      end if
      allocate(out%component(n, size(model%components)))
      out%component = 0.0_dp
      do component = 1, size(model%components)
         filter = tfarima_wiener_kolmogorov_filter(model, component, &
            selected_lag, tol)
         if (filter%info /= 0) then
            out%info = 3
            return
         end if
         do observation = 1, n
            do offset = -selected_lag, selected_lag
               index = selected_lag + observation + offset
               out%component(observation, component) = &
                  out%component(observation, component) + &
                  filter%weights(selected_lag + 1 + offset)*extended(index)
            end do
         end do
      end do
      allocate(out%reconstruction(n), out%remainder(n))
      out%reconstruction = sum(out%component, dim=2)
      out%remainder = series - out%reconstruction
   end function tfarima_ucarima_decompose

   pure function tfarima_extended_polynomial_gcd(first, second, tolerance) &
      result(out)
      ! Compute a polynomial GCD together with normalized Bezout coefficients.
      real(dp), intent(in) :: first(:), second(:)
      real(dp), intent(in), optional :: tolerance
      type(tfarima_bezout_t) :: out
      type(tfarima_polynomial_division_t) :: division
      real(dp), allocatable :: r0(:), r1(:), r2(:), s0(:), s1(:), s2(:)
      real(dp), allocatable :: t0(:), t1(:), t2(:)
      real(dp) :: scale, tol
      integer :: iteration, limit

      tol = 1.0e-10_dp
      if (present(tolerance)) tol = tolerance
      if (size(first) < 1 .or. size(second) < 1 .or. tol <= 0.0_dp .or. &
         .not. all(ieee_is_finite(first)) .or. &
         .not. all(ieee_is_finite(second))) then
         out%info = 1
         return
      end if
      r0 = trim_polynomial(first, tol)
      r1 = trim_polynomial(second, tol)
      s0 = [1.0_dp]
      s1 = [0.0_dp]
      t0 = [0.0_dp]
      t1 = [1.0_dp]
      limit = 2*(size(first) + size(second))
      do iteration = 1, limit
         if (polynomial_is_zero(r1, tol)) exit
         division = tfarima_polynomial_divide(r0, r1, tol)
         if (division%info /= 0) then
            out%info = 2
            return
         end if
         r2 = division%remainder
         s2 = polynomial_subtract(s0, &
            tfarima_polynomial_multiply(division%quotient, s1), tol)
         t2 = polynomial_subtract(t0, &
            tfarima_polynomial_multiply(division%quotient, t1), tol)
         r0 = r1
         r1 = r2
         s0 = s1
         s1 = s2
         t0 = t1
         t1 = t2
      end do
      if (.not. polynomial_is_zero(r1, tol) .or. polynomial_is_zero(r0, tol)) then
         out%info = 3
         return
      end if
      scale = r0(size(r0))
      if (abs(scale) <= tol) then
         out%info = 3
         return
      end if
      out%gcd = trim_polynomial(r0/scale, tol)
      out%u = trim_polynomial(s0/scale, tol)
      out%v = trim_polynomial(t0/scale, tol)
   end function tfarima_extended_polynomial_gcd

   pure function tfarima_partial_fractions(numerator, factors, factor_lengths, &
      tolerance) result(out)
      ! Split a proper rational numerator over pairwise-coprime factors.
      real(dp), intent(in) :: numerator(:), factors(:, :)
      integer, intent(in) :: factor_lengths(:)
      real(dp), intent(in), optional :: tolerance
      type(tfarima_partial_fraction_t) :: out
      type(tfarima_bezout_t) :: bezout
      real(dp), allocatable :: matrix(:, :), inverse(:, :), rhs(:), solution(:)
      real(dp), allocatable :: product(:), term(:), factor_i(:), factor_j(:)
      real(dp) :: tol
      integer :: count, total_degree, max_length, column, component, other
      integer :: lag, degree, first, last, status

      tol = 1.0e-10_dp
      if (present(tolerance)) tol = tolerance
      count = size(factor_lengths)
      if (count < 1 .or. size(factors, 2) /= count .or. tol <= 0.0_dp .or. &
         size(numerator) < 1 .or. any(factor_lengths < 2) .or. &
         any(factor_lengths > size(factors, 1)) .or. &
         .not. all(ieee_is_finite(numerator)) .or. &
         .not. all(ieee_is_finite(factors))) then
         out%info = 1
         return
      end if
      total_degree = sum(factor_lengths - 1)
      if (size(trim_polynomial(numerator, tol)) > total_degree) then
         out%info = 1
         return
      end if
      do component = 1, count - 1
         factor_i = factors(:factor_lengths(component), component)
         do other = component + 1, count
            factor_j = factors(:factor_lengths(other), other)
            bezout = tfarima_extended_polynomial_gcd(factor_i, factor_j, tol)
            if (bezout%info /= 0 .or. size(bezout%gcd) /= 1) then
               out%info = 2
               return
            end if
         end do
      end do
      allocate(matrix(total_degree, total_degree), rhs(total_degree))
      matrix = 0.0_dp
      rhs = 0.0_dp
      rhs(:size(trim_polynomial(numerator, tol))) = &
         trim_polynomial(numerator, tol)
      column = 0
      do component = 1, count
         product = [1.0_dp]
         do other = 1, count
            if (other == component) cycle
            product = tfarima_polynomial_multiply(product, &
               factors(:factor_lengths(other), other))
         end do
         degree = factor_lengths(component) - 1
         do lag = 0, degree - 1
            column = column + 1
            matrix(lag + 1:lag + size(product), column) = product
         end do
      end do
      call invert_matrix(matrix, inverse, status)
      if (status /= 0) then
         out%info = 3
         return
      end if
      solution = matmul(inverse, rhs)
      max_length = maxval(factor_lengths)
      allocate(out%numerator(max_length, count), &
         out%numerator_length(count))
      out%numerator = 0.0_dp
      first = 1
      do component = 1, count
         degree = factor_lengths(component) - 1
         last = first + degree - 1
         term = trim_polynomial(solution(first:last), tol)
         out%numerator(:size(term), component) = term
         out%numerator_length(component) = size(term)
         first = last + 1
      end do
      allocate(out%reconstruction(total_degree))
      out%reconstruction = 0.0_dp
      do component = 1, count
         product = [1.0_dp]
         do other = 1, count
            if (other == component) cycle
            product = tfarima_polynomial_multiply(product, &
               factors(:factor_lengths(other), other))
         end do
         term = tfarima_polynomial_multiply( &
            out%numerator(:out%numerator_length(component), component), product)
         out%reconstruction(:size(term)) = &
            out%reconstruction(:size(term)) + term
      end do
      out%reconstruction_error = maxval(abs(out%reconstruction - rhs))
      if (out%reconstruction_error > 100.0_dp*tol) out%info = 4
   end function tfarima_partial_fractions

   pure function tfarima_arima_to_ucarima(ar_polynomial, &
      difference_polynomial, ma_polynomial, innovation_variance, factors, &
      factor_lengths, factor_is_ar, canonical, tolerance, max_iterations) &
      result(out)
      ! Convert an aggregate ARIMA spectrum to independent UCARIMA components.
      real(dp), intent(in) :: ar_polynomial(:), difference_polynomial(:)
      real(dp), intent(in) :: ma_polynomial(:), innovation_variance
      real(dp), intent(in) :: factors(:, :)
      integer, intent(in) :: factor_lengths(:)
      logical, intent(in) :: factor_is_ar(:)
      logical, intent(in), optional :: canonical
      real(dp), intent(in), optional :: tolerance
      integer, intent(in), optional :: max_iterations
      type(tfarima_ucarima_conversion_t) :: out
      type(tfarima_polynomial_division_t) :: division
      type(tfarima_ma_factor_t) :: ma_factor
      type(tfarima_ucarima_component_t), allocatable :: components(:)
      real(dp), allocatable :: normalized(:, :), wold_factors(:, :)
      real(dp), allocatable :: ar_product(:), difference_product(:), gamma(:)
      real(dp), allocatable :: irregular_gamma(:)
      real(dp), allocatable :: numerator(:), source_covariance(:), adjusted(:)
      real(dp), allocatable :: product(:), term(:), target(:)
      real(dp) :: canonical_shift, minimum, tol
      integer :: count, component, length, total_components, limit, other
      integer :: total_degree
      logical :: make_canonical, has_irregular

      tol = 1.0e-8_dp
      if (present(tolerance)) tol = tolerance
      limit = 500
      if (present(max_iterations)) limit = max_iterations
      make_canonical = .false.
      if (present(canonical)) make_canonical = canonical
      count = size(factor_lengths)
      if (size(ar_polynomial) < 1 .or. size(difference_polynomial) < 1 .or. &
         size(ma_polynomial) < 1 .or. &
         abs(ar_polynomial(1)) <= tiny(1.0_dp) .or. &
         abs(difference_polynomial(1)) <= tiny(1.0_dp) .or. &
         abs(ma_polynomial(1)) <= tiny(1.0_dp) .or. &
         .not. all(ieee_is_finite(ar_polynomial)) .or. &
         .not. all(ieee_is_finite(difference_polynomial)) .or. &
         .not. all(ieee_is_finite(ma_polynomial)) .or. &
         .not. ieee_is_finite(innovation_variance) .or. &
         count < 1 .or. size(factor_is_ar) /= count .or. &
         size(factors, 2) /= count .or. any(factor_lengths < 2) .or. &
         any(factor_lengths > size(factors, 1)) .or. &
         innovation_variance <= 0.0_dp .or. tol <= 0.0_dp .or. limit < 1 .or. &
         .not. all(ieee_is_finite(factors))) then
         out%info = 1
         return
      end if
      allocate(normalized(size(factors, 1), count))
      normalized = 0.0_dp
      do component = 1, count
         if (abs(factors(1, component)) <= tiny(1.0_dp)) then
            out%info = 1
            return
         end if
         length = factor_lengths(component)
         normalized(:length, component) = &
            factors(:length, component)/factors(1, component)
      end do
      ar_product = [1.0_dp]
      difference_product = [1.0_dp]
      do component = 1, count
         if (factor_is_ar(component)) then
            ar_product = tfarima_polynomial_multiply(ar_product, &
               normalized(:factor_lengths(component), component))
         else
            difference_product = tfarima_polynomial_multiply( &
               difference_product, &
               normalized(:factor_lengths(component), component))
         end if
      end do
      if (polynomial_max_difference(ar_product, &
         ar_polynomial/ar_polynomial(1)) > tol .or. &
         polynomial_max_difference(difference_product, &
         difference_polynomial/difference_polynomial(1)) > tol) then
         out%info = 2
         return
      end if

      out%wold_numerator = wold_from_palindromic( &
         polynomial_autocovariance(ma_polynomial/ma_polynomial(1)))
      out%wold_denominator = wold_from_palindromic( &
         polynomial_autocovariance(tfarima_polynomial_multiply( &
         ar_product, difference_product)))
      division = tfarima_polynomial_divide(out%wold_numerator, &
         out%wold_denominator, tol)
      if (division%info /= 0) then
         out%info = 3
         return
      end if
      out%quotient = division%quotient
      out%remainder = division%remainder
      allocate(wold_factors(size(factors, 1), count))
      wold_factors = 0.0_dp
      do component = 1, count
         length = factor_lengths(component)
         adjusted = wold_from_palindromic(polynomial_autocovariance( &
            normalized(:length, component)))
         wold_factors(:size(adjusted), component) = adjusted
      end do
      out%partial_fraction = tfarima_partial_fractions(out%remainder, &
         wold_factors, factor_lengths, tol)
      if (out%partial_fraction%info /= 0) then
         out%info = 4
         return
      end if
      out%admissible = .true.
      canonical_shift = 0.0_dp
      do component = 1, count
         length = factor_lengths(component)
         minimum = rational_spectral_minimum( &
            out%partial_fraction%numerator(:, component), &
            wold_factors(:length, component), tol)
         if (.not. ieee_is_finite(minimum)) then
            out%info = 5
            out%admissible = .false.
            return
         end if
         if (minimum < -tol .and. .not. make_canonical) &
            out%admissible = .false.
         if (make_canonical) then
            out%partial_fraction%numerator(:length, component) = &
               out%partial_fraction%numerator(:length, component) - &
               minimum*wold_factors(:length, component)
            out%partial_fraction%numerator_length(component) = length
            out%quotient(1) = out%quotient(1) + minimum
            canonical_shift = canonical_shift + minimum
         end if
      end do
      if (.not. out%admissible) then
         out%info = 6
         return
      end if
      if (make_canonical) then
         total_degree = sum(factor_lengths - 1)
         allocate(target(total_degree + 1))
         target = 0.0_dp
         target(:size(out%remainder)) = out%remainder
         target(:size(out%wold_denominator)) = &
            target(:size(out%wold_denominator)) - &
            canonical_shift*out%wold_denominator
         out%remainder = trim_polynomial(target, tol)
         deallocate(out%partial_fraction%reconstruction)
         allocate(out%partial_fraction%reconstruction(total_degree + 1))
         out%partial_fraction%reconstruction = 0.0_dp
         do component = 1, count
            product = [1.0_dp]
            do other = 1, count
               if (other == component) cycle
               product = tfarima_polynomial_multiply(product, &
                  wold_factors(:factor_lengths(other), other))
            end do
            term = tfarima_polynomial_multiply( &
               out%partial_fraction%numerator( &
               :out%partial_fraction%numerator_length(component), component), &
               product)
            out%partial_fraction%reconstruction(:size(term)) = &
               out%partial_fraction%reconstruction(:size(term)) + term
         end do
         out%partial_fraction%reconstruction_error = &
            polynomial_max_difference( &
            out%partial_fraction%reconstruction, out%remainder)
      end if

      irregular_gamma = palindromic_from_wold(out%quotient)
      has_irregular = .not. polynomial_is_zero(irregular_gamma, tol)
      total_components = count + merge(1, 0, has_irregular)
      allocate(components(total_components))
      do component = 1, count
         numerator = trim_polynomial( &
            out%partial_fraction%numerator(:, component), tol)
         gamma = trim_polynomial(palindromic_from_wold(numerator), tol)
         ma_factor = tfarima_cramer_wold_factor(gamma, &
            tfarima_factor_best, tolerance=tol, max_iterations=limit)
         if (ma_factor%info /= 0) then
            out%info = 7
            out%admissible = .false.
            return
         end if
         if (factor_is_ar(component)) then
            components(component) = tfarima_ucarima_component( &
               normalized(:factor_lengths(component), component), [1.0_dp], &
               ma_factor%coefficients, innovation_variance)
         else
            components(component) = tfarima_ucarima_component([1.0_dp], &
               normalized(:factor_lengths(component), component), &
               ma_factor%coefficients, innovation_variance)
         end if
      end do
      if (has_irregular) then
         ma_factor = tfarima_cramer_wold_factor( &
            trim_polynomial(irregular_gamma, tol), tfarima_factor_best, &
            tolerance=tol, max_iterations=limit)
         if (ma_factor%info /= 0) then
            out%info = 8
            out%admissible = .false.
            return
         end if
         components(total_components) = tfarima_ucarima_component( &
            [1.0_dp], [1.0_dp], ma_factor%coefficients, innovation_variance)
      end if
      out%model = tfarima_build_ucarima(components, tol, limit)
      if (out%model%info /= 0) then
         out%info = 9
         out%admissible = .false.
         return
      end if
      source_covariance = innovation_variance*polynomial_autocovariance( &
         ma_polynomial/ma_polynomial(1))
      out%reconstruction_error = polynomial_max_difference( &
         out%model%numerator_autocovariance, source_covariance)
      if (out%reconstruction_error > 100.0_dp*tol) then
         out%info = 10
         return
      end if
   end function tfarima_arima_to_ucarima

   pure function tfarima_ucarima_state_space(series, ucarima, tolerance) &
      result(out)
      ! Assemble independent UCARIMA components into one KFAS-compatible model.
      real(dp), intent(in) :: series(:)
      type(tfarima_ucarima_model_t), intent(in) :: ucarima
      real(dp), intent(in), optional :: tolerance
      type(tfarima_ucarima_ssm_t) :: out
      real(dp), allocatable :: denominator(:), ma_polynomial(:)
      real(dp), allocatable :: transition(:, :), loading(:)
      real(dp), allocatable :: covariance(:, :)
      real(dp) :: tol
      integer :: component, count, degree, ma_degree, states, first, last, i
      logical :: converged, integrated

      tol = 1.0e-10_dp
      if (present(tolerance)) tol = tolerance
      if (size(series) < 1 .or. ucarima%info /= 0 .or. &
         .not. allocated(ucarima%components) .or. &
         tol <= 0.0_dp) then
         out%info = 1
         return
      end if
      count = size(ucarima%components)
      allocate(out%state_start(count), out%state_size(count))
      states = 0
      do component = 1, count
         denominator = trim_polynomial(tfarima_polynomial_multiply( &
            ucarima%components(component)%ar_polynomial, &
            ucarima%components(component)%difference_polynomial), tol)
         ma_polynomial = trim_polynomial( &
            ucarima%components(component)%ma_polynomial, tol)
         degree = size(denominator) - 1
         ma_degree = size(ma_polynomial) - 1
         out%state_start(component) = states + 1
         out%state_size(component) = max(1, degree, ma_degree + 1)
         states = states + out%state_size(component)
      end do
      allocate(out%component_loading(count, states))
      out%component_loading = 0.0_dp
      allocate(out%model%y(size(series), 1), out%model%z(1, states, 1))
      allocate(out%model%h(1, 1, 1), out%model%transition(states, states, 1))
      allocate(out%model%r(states, count, 1), out%model%q(count, count, 1))
      allocate(out%model%a1(states), out%model%p1(states, states))
      allocate(out%model%p1inf(states, states), &
         out%model%missing(size(series), 1))
      out%model%y(:, 1) = series
      out%model%z = 0.0_dp
      out%model%h = 0.0_dp
      out%model%transition = 0.0_dp
      out%model%r = 0.0_dp
      out%model%q = 0.0_dp
      out%model%a1 = 0.0_dp
      out%model%p1 = 0.0_dp
      out%model%p1inf = 0.0_dp
      out%model%missing(:, 1) = .not. ieee_is_finite(series)
      do component = 1, count
         first = out%state_start(component)
         last = first + out%state_size(component) - 1
         denominator = tfarima_polynomial_multiply( &
            ucarima%components(component)%ar_polynomial, &
            ucarima%components(component)%difference_polynomial)
         denominator = trim_polynomial(denominator, tol)
         ma_polynomial = trim_polynomial( &
            ucarima%components(component)%ma_polynomial, tol)
         allocate(transition(out%state_size(component), &
            out%state_size(component)), loading(out%state_size(component)))
         transition = 0.0_dp
         loading = 0.0_dp
         if (size(denominator) > 1) &
            transition(:size(denominator) - 1, 1) = -denominator(2:)
         do i = 1, out%state_size(component) - 1
            transition(i, i + 1) = 1.0_dp
         end do
         loading(1) = 1.0_dp
         if (size(ma_polynomial) > 1) then
            loading(2:size(ma_polynomial)) = ma_polynomial(2:)
         end if
         out%model%transition(first:last, first:last, 1) = transition
         out%model%r(first:last, component, 1) = loading
         out%model%q(component, component, 1) = &
            ucarima%components(component)%innovation_variance
         out%component_loading(component, first) = 1.0_dp
         out%model%z(1, first, 1) = 1.0_dp
         integrated = size( &
            ucarima%components(component)%difference_polynomial) > 1
         if (integrated) then
            out%diffuse = .true.
            do i = first, last
               out%model%p1inf(i, i) = 1.0_dp
               out%model%p1(i, i) = tol
            end do
            out%model%p1(first:last, first:last) = &
               out%model%p1(first:last, first:last) + &
               ucarima%components(component)%innovation_variance* &
               spread(loading, 2, size(loading))* &
               spread(loading, 1, size(loading))
         else
            call stationary_state_covariance(transition, loading, &
               ucarima%components(component)%innovation_variance, tol, &
               covariance, converged)
            if (.not. converged) then
               out%info = 2
               return
            end if
            out%model%p1(first:last, first:last) = covariance
         end if
         deallocate(transition, loading)
      end do
   end function tfarima_ucarima_state_space

   pure function tfarima_ucarima_smooth(series, ucarima, forecast_horizon, &
      levels, tolerance) result(out)
      ! Smooth UCARIMA components and propagate component forecast uncertainty.
      real(dp), intent(in) :: series(:)
      type(tfarima_ucarima_model_t), intent(in) :: ucarima
      integer, intent(in), optional :: forecast_horizon
      real(dp), intent(in), optional :: levels(:), tolerance
      type(tfarima_ucarima_smoothing_t) :: out
      type(tfarima_ucarima_ssm_t) :: state_space
      type(kfs_filter_t) :: filtered
      type(kfs_smoother_t) :: smoothed
      real(dp), allocatable :: state(:), covariance(:, :), process(:, :)
      real(dp), allocatable :: next_state(:), next_covariance(:, :)
      real(dp) :: tol, quantile
      integer :: horizon, component, observation, level, count, h

      tol = 1.0e-10_dp
      if (present(tolerance)) tol = tolerance
      horizon = 0
      if (present(forecast_horizon)) horizon = forecast_horizon
      if (horizon < 0 .or. tol <= 0.0_dp) then
         out%info = 1
         return
      end if
      if (present(levels)) then
         if (any(levels <= 0.0_dp) .or. any(levels >= 1.0_dp)) then
            out%info = 1
            return
         end if
      end if
      state_space = tfarima_ucarima_state_space(series, ucarima, tol)
      if (state_space%info /= 0) then
         out%info = 2
         return
      end if
      if (state_space%diffuse) then
         filtered = kfs_filter_diffuse(state_space%model, tol)
      else
         filtered = kfs_filter(state_space%model, tol)
      end if
      if (filtered%info /= 0) then
         out%info = 3
         return
      end if
      smoothed = kfs_smooth(state_space%model, filtered, tol)
      if (smoothed%info /= 0) then
         out%info = 4
         return
      end if
      count = size(ucarima%components)
      allocate(out%component(size(series), count), &
         out%variance(size(series), count))
      do component = 1, count
         do observation = 1, size(series)
            out%component(observation, component) = dot_product( &
               state_space%component_loading(component, :), &
               smoothed%state(:, observation))
            out%variance(observation, component) = max(0.0_dp, dot_product( &
               state_space%component_loading(component, :), matmul( &
               smoothed%covariance(:, :, observation), &
               state_space%component_loading(component, :))))
         end do
      end do
      out%log_likelihood = filtered%log_likelihood
      allocate(out%forecast(horizon, count), &
         out%forecast_variance(horizon, count))
      if (horizon < 1) return
      state = filtered%a_filt(:, size(series))
      covariance = filtered%p_filt(:, :, size(series))
      process = matmul(matmul(state_space%model%r(:, :, 1), &
         state_space%model%q(:, :, 1)), &
         transpose(state_space%model%r(:, :, 1)))
      do h = 1, horizon
         next_state = matmul(state_space%model%transition(:, :, 1), state)
         next_covariance = matmul(matmul( &
            state_space%model%transition(:, :, 1), covariance), &
            transpose(state_space%model%transition(:, :, 1))) + process
         do component = 1, count
            out%forecast(h, component) = dot_product( &
               state_space%component_loading(component, :), next_state)
            out%forecast_variance(h, component) = max(0.0_dp, dot_product( &
               state_space%component_loading(component, :), matmul( &
               next_covariance, &
               state_space%component_loading(component, :))))
         end do
         state = next_state
         covariance = next_covariance
      end do
      if (present(levels)) then
         allocate(out%lower(horizon, count, size(levels)), &
            out%upper(horizon, count, size(levels)))
         do level = 1, size(levels)
            quantile = normal_quantile(0.5_dp + 0.5_dp*levels(level))
            out%lower(:, :, level) = out%forecast - &
               quantile*sqrt(out%forecast_variance)
            out%upper(:, :, level) = out%forecast + &
               quantile*sqrt(out%forecast_variance)
         end do
      end if
   end function tfarima_ucarima_smooth

   pure function tfarima_ucarima_fit(series, initial_model, estimate_ar, &
      estimate_ma, estimate_variance, regressors, initial_regression, &
      max_iterations, tolerance) result(out)
      ! Estimate UCARIMA dynamics and variances by exact state-space likelihood.
      real(dp), intent(in) :: series(:)
      type(tfarima_ucarima_model_t), intent(in) :: initial_model
      logical, intent(in), optional :: estimate_ar(:, :), estimate_ma(:, :)
      logical, intent(in), optional :: estimate_variance(:)
      real(dp), intent(in), optional :: regressors(:, :), initial_regression(:)
      integer, intent(in), optional :: max_iterations
      real(dp), intent(in), optional :: tolerance
      type(tfarima_ucarima_fit_t) :: out
      type(optimization_result_t) :: warm_start, optimization
      type(kfs_filter_t) :: filtered
      type(tfarima_ucarima_component_t), allocatable :: fitted_components(:)
      logical, allocatable :: ar_mask(:, :), ma_mask(:, :), variance_mask(:)
      logical, allocatable :: variance_parameter(:)
      real(dp), allocatable :: initial(:), coordinates(:), beta(:), adjusted(:)
      real(dp), allocatable :: hessian(:, :), inverse(:, :), jacobian(:, :)
      real(dp) :: tol
      integer :: component_count, max_ar, max_ma, nreg, npar
      integer :: position, component, lag
      integer :: limit, status, observations, i, bfgs_iterations

      tol = 1.0e-6_dp
      if (present(tolerance)) tol = tolerance
      limit = 250
      if (present(max_iterations)) limit = max_iterations
      if (size(series) < 3 .or. initial_model%info /= 0 .or. &
         .not. allocated(initial_model%components) .or. tol <= 0.0_dp .or. &
         limit < 1 .or. .not. all(ieee_is_finite(series))) then
         out%info = 1
         return
      end if
      component_count = size(initial_model%components)
      max_ar = 0
      max_ma = 0
      do component = 1, component_count
         max_ar = max(max_ar, &
            size(initial_model%components(component)%ar_polynomial) - 1)
         max_ma = max(max_ma, &
            size(initial_model%components(component)%ma_polynomial) - 1)
      end do
      allocate(ar_mask(max_ar, component_count), &
         ma_mask(max_ma, component_count), variance_mask(component_count))
      ar_mask = .false.
      ma_mask = .false.
      variance_mask = .true.
      if (present(estimate_ar)) then
         if (any(shape(estimate_ar) /= [max_ar, component_count])) then
            out%info = 2
            return
         end if
         ar_mask = estimate_ar
      end if
      if (present(estimate_ma)) then
         if (any(shape(estimate_ma) /= [max_ma, component_count])) then
            out%info = 2
            return
         end if
         ma_mask = estimate_ma
      end if
      if (present(estimate_variance)) then
         if (size(estimate_variance) /= component_count) then
            out%info = 2
            return
         end if
         variance_mask = estimate_variance
      end if
      do component = 1, component_count
         lag = size(initial_model%components(component)%ar_polynomial) - 1
         if (lag < max_ar) then
            if (any(ar_mask(lag + 1:, component))) then
               out%info = 2
               return
            end if
         end if
         lag = size(initial_model%components(component)%ma_polynomial) - 1
         if (lag < max_ma) then
            if (any(ma_mask(lag + 1:, component))) then
               out%info = 2
               return
            end if
         end if
      end do
      nreg = 0
      if (present(regressors)) then
         if (size(regressors, 1) /= size(series) .or. &
            .not. all(ieee_is_finite(regressors))) then
            out%info = 3
            return
         end if
         nreg = size(regressors, 2)
      end if
      if (present(initial_regression)) then
         if (.not. present(regressors) .or. size(initial_regression) /= nreg) then
            out%info = 3
            return
         end if
      end if
      npar = nreg + count(ar_mask) + count(ma_mask) + count(variance_mask)
      allocate(initial(npar), variance_parameter(npar))
      variance_parameter = .false.
      position = 0
      if (nreg > 0) then
         initial(1:nreg) = 0.0_dp
         if (present(initial_regression)) initial(1:nreg) = initial_regression
         position = nreg
      end if
      do component = 1, component_count
         do lag = 1, size(initial_model%components(component)%ar_polynomial) - 1
            if (.not. ar_mask(lag, component)) cycle
            position = position + 1
            initial(position) = &
               initial_model%components(component)%ar_polynomial(lag + 1)
         end do
      end do
      do component = 1, component_count
         do lag = 1, size(initial_model%components(component)%ma_polynomial) - 1
            if (.not. ma_mask(lag, component)) cycle
            position = position + 1
            initial(position) = &
               initial_model%components(component)%ma_polynomial(lag + 1)
         end do
      end do
      do component = 1, component_count
         if (.not. variance_mask(component)) cycle
         position = position + 1
         initial(position) = 0.5_dp*log(max( &
            initial_model%components(component)%innovation_variance, tol))
         variance_parameter(position) = .true.
      end do

      if (npar > 0) then
         warm_start = nelder_mead_minimize(objective, initial, &
            max_iterations=limit, tolerance=tol, initial_step=0.05_dp)
         optimization = bfgs_minimize_fd(objective, warm_start%parameters, &
            max_iterations=limit, gradient_tolerance=tol)
         bfgs_iterations = optimization%iterations
         if (.not. optimization%converged .and. warm_start%converged) then
            optimization = warm_start
         end if
         coordinates = optimization%parameters
         out%iterations = warm_start%iterations + bfgs_iterations
         out%converged = optimization%converged
      else
         allocate(coordinates(0))
         out%converged = .true.
      end if
      call evaluate(coordinates, fitted_components, beta, adjusted, filtered, status)
      if (status /= 0) then
         out%info = 10 + status
         return
      end if
      out%model = tfarima_build_ucarima(fitted_components, tol, 500)
      if (out%model%info /= 0) then
         out%info = 20 + out%model%info
         return
      end if
      out%regression_coefficients = beta
      out%log_likelihood = filtered%log_likelihood
      observations = count(ieee_is_finite(filtered%innovation(:, 1)))
      out%aic = -2.0_dp*out%log_likelihood + 2.0_dp*real(npar, dp)
      out%bic = -2.0_dp*out%log_likelihood + &
         log(real(max(1, observations), dp))*real(npar, dp)
      allocate(out%residuals(size(series)), &
         out%standardized_residuals(size(series)))
      out%residuals = filtered%innovation(:, 1)
      out%standardized_residuals = 0.0_dp
      do i = 1, size(series)
         if (.not. ieee_is_finite(out%residuals(i))) cycle
         if (filtered%innovation_cov(1, 1, i) <= 0.0_dp) cycle
         out%standardized_residuals(i) = out%residuals(i)/ &
            sqrt(filtered%innovation_cov(1, 1, i))
      end do
      allocate(out%parameters(npar), out%covariance(npar, npar), &
         out%standard_error(npar))
      call natural_parameters(coordinates, out%parameters)
      out%covariance = 0.0_dp
      out%standard_error = 0.0_dp
      if (npar > 0) then
         hessian = finite_difference_hessian(objective, coordinates)
         call invert_matrix(hessian, inverse, status)
         if (status == 0) then
            allocate(jacobian(npar, npar))
            jacobian = 0.0_dp
            do i = 1, npar
               jacobian(i, i) = 1.0_dp
               if (variance_parameter(i)) &
                  jacobian(i, i) = 2.0_dp*out%parameters(i)
            end do
            out%covariance = matmul(matmul(jacobian, inverse), &
               transpose(jacobian))
            do i = 1, npar
               out%standard_error(i) = &
                  sqrt(max(0.0_dp, out%covariance(i, i)))
            end do
         end if
      end if
      if (.not. out%converged .and. npar > 0) out%info = 30

   contains

      pure function objective(values) result(value)
         ! Return the negative exact Gaussian likelihood for optimizer values.
         real(dp), intent(in) :: values(:)
         real(dp) :: value
         type(tfarima_ucarima_component_t), allocatable :: trial_components(:)
         type(kfs_filter_t) :: trial_filter
         real(dp), allocatable :: trial_beta(:), trial_series(:)
         integer :: evaluation_status

         call evaluate(values, trial_components, trial_beta, trial_series, &
            trial_filter, evaluation_status)
         if (evaluation_status == 0 .and. &
            ieee_is_finite(trial_filter%log_likelihood)) then
            value = -trial_filter%log_likelihood
         else
            value = 1.0e30_dp + dot_product(values, values)
         end if
      end function objective

      pure subroutine evaluate(values, components, regression, transformed, &
         filter_result, evaluation_status)
         ! Decode parameters and run the shared ordinary or diffuse KFAS filter.
         real(dp), intent(in) :: values(:)
         type(tfarima_ucarima_component_t), allocatable, intent(out) :: components(:)
         real(dp), allocatable, intent(out) :: regression(:), transformed(:)
         type(kfs_filter_t), intent(out) :: filter_result
         integer, intent(out) :: evaluation_status
         type(tfarima_ucarima_model_t) :: trial_model
         type(tfarima_ucarima_ssm_t) :: state_space
         real(dp), allocatable :: ar(:), ma(:)
         integer :: trial_position, trial_component, trial_lag

         evaluation_status = 1
         allocate(components(component_count), regression(nreg))
         trial_position = nreg
         if (nreg > 0) regression = values(:nreg)
         do trial_component = 1, component_count
            ar = initial_model%components(trial_component)%ar_polynomial
            do trial_lag = 1, size(ar) - 1
               if (.not. ar_mask(trial_lag, trial_component)) cycle
               trial_position = trial_position + 1
               ar(trial_lag + 1) = values(trial_position)
            end do
            if (.not. tfarima_polynomial_admissible(ar)) return
            components(trial_component)%ar_polynomial = ar
         end do
         do trial_component = 1, component_count
            ma = initial_model%components(trial_component)%ma_polynomial
            do trial_lag = 1, size(ma) - 1
               if (.not. ma_mask(trial_lag, trial_component)) cycle
               trial_position = trial_position + 1
               ma(trial_lag + 1) = values(trial_position)
            end do
            if (.not. tfarima_polynomial_admissible(ma)) return
            components(trial_component)%ma_polynomial = ma
            components(trial_component)%difference_polynomial = &
               initial_model%components(trial_component)%difference_polynomial
            components(trial_component)%innovation_variance = &
               initial_model%components(trial_component)%innovation_variance
         end do
         do trial_component = 1, component_count
            if (.not. variance_mask(trial_component)) cycle
            trial_position = trial_position + 1
            if (abs(values(trial_position)) > 350.0_dp) return
            components(trial_component)%innovation_variance = &
               exp(2.0_dp*values(trial_position))
         end do
         transformed = series
         if (nreg > 0) transformed = transformed - matmul(regressors, regression)
         allocate(trial_model%components(component_count))
         trial_model%components = components
         trial_model%info = 0
         state_space = tfarima_ucarima_state_space(transformed, trial_model, tol)
         if (state_space%info /= 0) then
            evaluation_status = 2
            return
         end if
         if (state_space%diffuse) then
            filter_result = kfs_filter_diffuse(state_space%model, tol)
         else
            filter_result = kfs_filter(state_space%model, tol)
         end if
         if (filter_result%info /= 0) then
            evaluation_status = 3
            return
         end if
         evaluation_status = 0
      end subroutine evaluate

      pure subroutine natural_parameters(values, natural)
         ! Transform log-standard-deviation coordinates to component variances.
         real(dp), intent(in) :: values(:)
         real(dp), intent(out) :: natural(:)
         integer :: parameter

         natural = values
         do parameter = 1, size(values)
            if (variance_parameter(parameter)) &
               natural(parameter) = exp(2.0_dp*values(parameter))
         end do
      end subroutine natural_parameters

   end function tfarima_ucarima_fit

   pure function tfarima_polynomial_multiply(first, second) result(product)
      ! Multiply two increasing-lag coefficient vectors.
      real(dp), intent(in) :: first(:), second(:)
      real(dp), allocatable :: product(:)
      integer :: i, j

      if (size(first) < 1 .or. size(second) < 1) then
         allocate(product(0))
         return
      end if
      allocate(product(size(first) + size(second) - 1))
      product = 0.0_dp
      do i = 1, size(first)
         do j = 1, size(second)
            product(i + j - 1) = product(i + j - 1) + first(i)*second(j)
         end do
      end do
   end function tfarima_polynomial_multiply

   pure function tfarima_polynomial_divide(numerator, denominator, tolerance) &
      result(out)
      ! Divide lag polynomials and return trimmed quotient and remainder.
      real(dp), intent(in) :: numerator(:), denominator(:)
      real(dp), intent(in), optional :: tolerance
      type(tfarima_polynomial_division_t) :: out
      real(dp), allocatable :: work(:), divisor(:), quotient(:)
      real(dp) :: tol, scale
      integer :: numerator_degree, denominator_degree, quotient_degree, i, j

      tol = 1.0e-10_dp
      if (present(tolerance)) tol = tolerance
      work = trim_polynomial(numerator, tol)
      divisor = trim_polynomial(denominator, tol)
      if (size(work) < 1 .or. size(divisor) < 1 .or. &
         maxval(abs(divisor)) <= tol) then
         out%info = 1
         allocate(out%quotient(0), out%remainder(0))
         return
      end if
      numerator_degree = size(work) - 1
      denominator_degree = size(divisor) - 1
      if (numerator_degree < denominator_degree) then
         out%quotient = [0.0_dp]
         out%remainder = work
         return
      end if
      quotient_degree = numerator_degree - denominator_degree
      allocate(quotient(quotient_degree + 1))
      quotient = 0.0_dp
      do i = quotient_degree, 0, -1
         scale = work(denominator_degree + i + 1)/divisor(denominator_degree + 1)
         quotient(i + 1) = scale
         do j = 0, denominator_degree
            work(i + j + 1) = work(i + j + 1) - scale*divisor(j + 1)
         end do
      end do
      out%quotient = trim_polynomial(quotient, tol)
      if (denominator_degree == 0) then
         out%remainder = [0.0_dp]
      else
         out%remainder = trim_polynomial(work(:denominator_degree), tol)
      end if
   end function tfarima_polynomial_divide

   pure function tfarima_polynomial_gcd(first, second, tolerance) result(gcd)
      ! Compute the monic-at-lag-zero Euclidean polynomial GCD.
      real(dp), intent(in) :: first(:), second(:)
      real(dp), intent(in), optional :: tolerance
      real(dp), allocatable :: gcd(:), left(:), right(:), swap(:)
      type(tfarima_polynomial_division_t) :: division
      real(dp) :: tol
      integer :: iteration, limit

      tol = 1.0e-10_dp
      if (present(tolerance)) tol = tolerance
      left = trim_polynomial(first, tol)
      right = trim_polynomial(second, tol)
      if (size(left) < size(right)) then
         swap = left
         left = right
         right = swap
      end if
      limit = max(1, size(left) + size(right))
      do iteration = 1, limit
         if (size(right) == 1 .and. abs(right(1)) <= tol) exit
         division = tfarima_polynomial_divide(left, right, tol)
         if (division%info /= 0) then
            allocate(gcd(0))
            return
         end if
         left = right
         right = division%remainder
      end do
      if (size(left) < 1 .or. abs(left(1)) <= tol) then
         allocate(gcd(0))
      else
         gcd = left/left(1)
      end if
   end function tfarima_polynomial_gcd

   pure function tfarima_polynomial_power(polynomial, exponent) result(power)
      ! Raise a lag polynomial to a nonnegative integer power.
      real(dp), intent(in) :: polynomial(:)
      integer, intent(in) :: exponent
      real(dp), allocatable :: power(:), factor(:)
      integer :: i

      if (exponent < 0 .or. size(polynomial) < 1) then
         allocate(power(0))
         return
      end if
      power = [1.0_dp]
      factor = polynomial
      do i = 1, exponent
         power = tfarima_polynomial_multiply(power, factor)
      end do
   end function tfarima_polynomial_power

   pure function tfarima_polynomial_ratio(numerator, denominator, degree) result(ratio)
      ! Expand a rational lag polynomial through a requested degree.
      real(dp), intent(in) :: numerator(:), denominator(:)
      integer, intent(in) :: degree
      real(dp), allocatable :: ratio(:)
      real(dp) :: value
      integer :: lag, j

      if (degree < 0 .or. size(numerator) < 1 .or. size(denominator) < 1 .or. &
         abs(denominator(1)) <= tiny(1.0_dp)) then
         allocate(ratio(0))
         return
      end if
      allocate(ratio(degree + 1))
      ratio = 0.0_dp
      do lag = 0, degree
         value = 0.0_dp
         do j = 1, min(lag, size(denominator) - 1)
            value = value + denominator(j + 1)*ratio(lag - j + 1)
         end do
         if (lag < size(numerator)) then
            ratio(lag + 1) = (numerator(lag + 1) - value)/denominator(1)
         else
            ratio(lag + 1) = -value/denominator(1)
         end if
      end do
   end function tfarima_polynomial_ratio

   pure function tfarima_combine_arima(first, second, subtract_second, &
      tolerance, max_iterations) result(out)
      ! Add or subtract two ARIMA spectra and cancel common polynomial factors.
      type(tfarima_ucarima_component_t), intent(in) :: first, second
      logical, intent(in), optional :: subtract_second
      real(dp), intent(in), optional :: tolerance
      integer, intent(in), optional :: max_iterations
      type(tfarima_arima_algebra_t) :: out
      type(tfarima_polynomial_division_t) :: division
      real(dp), allocatable :: ar(:), difference(:), denominator(:)
      real(dp), allocatable :: first_denominator(:), second_denominator(:)
      real(dp), allocatable :: first_ma(:), second_ma(:), gcd(:)
      real(dp), allocatable :: ma(:), reconstructed(:)
      real(dp) :: tol, second_sign, scale
      integer :: limit, lag, length

      tol = 1.0e-10_dp
      if (present(tolerance)) tol = tolerance
      limit = 500
      if (present(max_iterations)) limit = max_iterations
      second_sign = 1.0_dp
      if (present(subtract_second)) then
         if (subtract_second) second_sign = -1.0_dp
      end if
      if (.not. valid_component(first) .or. .not. valid_component(second) .or. &
         tol <= 0.0_dp .or. limit < 1) then
         out%info = 1
         return
      end if

      ar = polynomial_lcm(first%ar_polynomial, second%ar_polynomial, tol)
      difference = polynomial_lcm(first%difference_polynomial, &
         second%difference_polynomial, tol)
      if (size(ar) < 1 .or. size(difference) < 1) then
         out%info = 2
         return
      end if
      denominator = tfarima_polynomial_multiply(ar, difference)
      first_denominator = tfarima_polynomial_multiply(first%ar_polynomial, &
         first%difference_polynomial)
      second_denominator = tfarima_polynomial_multiply(second%ar_polynomial, &
         second%difference_polynomial)
      division = tfarima_polynomial_divide(denominator, first_denominator, tol)
      if (.not. exact_division(division, tol)) then
         out%info = 2
         return
      end if
      first_ma = tfarima_polynomial_multiply(first%ma_polynomial, &
         division%quotient)
      division = tfarima_polynomial_divide(denominator, second_denominator, tol)
      if (.not. exact_division(division, tol)) then
         out%info = 2
         return
      end if
      second_ma = tfarima_polynomial_multiply(second%ma_polynomial, &
         division%quotient)
      length = max(size(first_ma), size(second_ma))
      allocate(out%target_autocovariance(length))
      out%target_autocovariance = 0.0_dp
      do lag = 0, length - 1
         if (lag < size(first_ma)) out%target_autocovariance(lag + 1) = &
            first%innovation_variance*dot_product( &
            first_ma(:size(first_ma) - lag), first_ma(lag + 1:))
         if (lag < size(second_ma)) out%target_autocovariance(lag + 1) = &
            out%target_autocovariance(lag + 1) + second_sign* &
            second%innovation_variance*dot_product( &
            second_ma(:size(second_ma) - lag), second_ma(lag + 1:))
      end do
      out%factor = tfarima_cramer_wold_factor(out%target_autocovariance, &
         tfarima_factor_best, tolerance=tol, max_iterations=limit)
      if (out%factor%info /= 0 .or. .not. allocated(out%factor%coefficients) .or. &
         abs(out%factor%coefficients(1)) <= tol) then
         out%info = 3
         return
      end if
      allocate(reconstructed(size(out%factor%coefficients)))
      call factor_moments(out%factor%coefficients, reconstructed)
      out%reconstructed_autocovariance = reconstructed
      out%reconstruction_error = sqrt(sum((reconstructed - &
         out%target_autocovariance(:size(reconstructed)))**2))
      scale = out%factor%coefficients(1)
      ma = out%factor%coefficients/scale

      call cancel_factor(ma, difference, tol)
      call cancel_factor(ma, ar, tol)
      out%model = tfarima_ucarima_component(ar, difference, ma, scale**2)
      if (out%model%info /= 0) then
         out%info = 4
         return
      end if

   contains

      pure logical function valid_component(component) result(valid)
         ! Report whether a component has finite normalized polynomials.
         type(tfarima_ucarima_component_t), intent(in) :: component

         valid = component%info == 0 .and. &
            allocated(component%ar_polynomial) .and. &
            allocated(component%difference_polynomial) .and. &
            allocated(component%ma_polynomial)
         if (.not. valid) return
         valid = component%innovation_variance >= 0.0_dp .and. &
            ieee_is_finite(component%innovation_variance) .and. &
            all(ieee_is_finite(component%ar_polynomial)) .and. &
            all(ieee_is_finite(component%difference_polynomial)) .and. &
            all(ieee_is_finite(component%ma_polynomial))
      end function valid_component

      pure function polynomial_lcm(left, right, local_tolerance) result(lcm)
         ! Form a normalized least-common multiple of two lag polynomials.
         real(dp), intent(in) :: left(:), right(:), local_tolerance
         real(dp), allocatable :: lcm(:), common(:)
         type(tfarima_polynomial_division_t) :: quotient

         common = tfarima_polynomial_gcd(left, right, local_tolerance)
         if (size(common) < 1) then
            allocate(lcm(0))
            return
         end if
         quotient = tfarima_polynomial_divide(right, common, local_tolerance)
         if (.not. exact_division(quotient, local_tolerance)) then
            allocate(lcm(0))
            return
         end if
         lcm = trim_polynomial(tfarima_polynomial_multiply(left, &
            quotient%quotient), local_tolerance)
         if (size(lcm) > 0 .and. abs(lcm(1)) > local_tolerance) lcm = lcm/lcm(1)
      end function polynomial_lcm

      pure logical function exact_division(quotient, local_tolerance) &
         result(exact)
         ! Report whether polynomial division has a negligible remainder.
         type(tfarima_polynomial_division_t), intent(in) :: quotient
         real(dp), intent(in) :: local_tolerance

         exact = quotient%info == 0 .and. allocated(quotient%remainder)
         if (exact) exact = maxval(abs(quotient%remainder)) <= local_tolerance
      end function exact_division

      pure subroutine cancel_factor(numerator, denominator_polynomial, &
         local_tolerance)
         ! Cancel a nonconstant common factor from numerator and denominator.
         real(dp), allocatable, intent(inout) :: numerator(:)
         real(dp), allocatable, intent(inout) :: denominator_polynomial(:)
         real(dp), intent(in) :: local_tolerance
         real(dp), allocatable :: common(:)
         type(tfarima_polynomial_division_t) :: numerator_division
         type(tfarima_polynomial_division_t) :: denominator_division

         common = tfarima_polynomial_gcd(numerator, denominator_polynomial, &
            local_tolerance)
         if (size(common) <= 1) return
         numerator_division = tfarima_polynomial_divide(numerator, common, &
            local_tolerance)
         denominator_division = tfarima_polynomial_divide( &
            denominator_polynomial, common, local_tolerance)
         if (.not. exact_division(numerator_division, local_tolerance) .or. &
            .not. exact_division(denominator_division, local_tolerance)) return
         numerator = numerator_division%quotient
         denominator_polynomial = denominator_division%quotient
      end subroutine cancel_factor

   end function tfarima_combine_arima

   pure function tfarima_psi_weights(model, max_lag, include_difference, &
      coefficient_covariance) result(out)
      ! Expand finite ARIMA PSI weights with optional delta-method inference.
      type(tfarima_ucarima_component_t), intent(in) :: model
      integer, intent(in) :: max_lag
      logical, intent(in), optional :: include_difference
      real(dp), intent(in), optional :: coefficient_covariance(:, :)
      type(tfarima_weights_t) :: out

      if (present(coefficient_covariance)) then
         out = arima_weights(model, max_lag, include_difference, .false., &
            coefficient_covariance)
      else
         out = arima_weights(model, max_lag, include_difference, .false.)
      end if
   end function tfarima_psi_weights

   pure function tfarima_pi_weights(model, max_lag, include_difference, &
      coefficient_covariance) result(out)
      ! Expand finite ARIMA PI weights with optional delta-method inference.
      type(tfarima_ucarima_component_t), intent(in) :: model
      integer, intent(in) :: max_lag
      logical, intent(in), optional :: include_difference
      real(dp), intent(in), optional :: coefficient_covariance(:, :)
      type(tfarima_weights_t) :: out

      if (present(coefficient_covariance)) then
         out = arima_weights(model, max_lag, include_difference, .true., &
            coefficient_covariance)
      else
         out = arima_weights(model, max_lag, include_difference, .true.)
      end if
   end function tfarima_pi_weights

   pure function arima_weights(model, max_lag, include_difference, inverse, &
      coefficient_covariance) result(out)
      ! Compute rational response weights and finite-difference delta inference.
      type(tfarima_ucarima_component_t), intent(in) :: model
      integer, intent(in) :: max_lag
      logical, intent(in), optional :: include_difference
      logical, intent(in) :: inverse
      real(dp), intent(in), optional :: coefficient_covariance(:, :)
      type(tfarima_weights_t) :: out
      real(dp), allocatable :: ar(:), ma(:), numerator(:), denominator(:)
      real(dp), allocatable :: trial(:), jacobian(:, :), cumulative_jacobian(:, :)
      real(dp) :: step
      integer :: p, q, parameters, column, lag
      logical :: integrated

      integrated = .false.
      if (present(include_difference)) integrated = include_difference
      if (model%info /= 0 .or. max_lag < 0 .or. &
         .not. allocated(model%ar_polynomial) .or. &
         .not. allocated(model%difference_polynomial) .or. &
         .not. allocated(model%ma_polynomial)) then
         out%info = 1
         return
      end if
      ar = model%ar_polynomial
      if (integrated) ar = tfarima_polynomial_multiply(ar, &
         model%difference_polynomial)
      ma = model%ma_polynomial
      if (inverse) then
         numerator = ar
         denominator = ma
      else
         numerator = ma
         denominator = ar
      end if
      out%weight = tfarima_polynomial_ratio(numerator, denominator, max_lag)
      if (size(out%weight) /= max_lag + 1) then
         out%info = 1
         return
      end if
      allocate(out%cumulative(max_lag + 1))
      out%cumulative(1) = out%weight(1)
      do lag = 2, max_lag + 1
         out%cumulative(lag) = out%cumulative(lag - 1) + out%weight(lag)
      end do

      p = size(model%ar_polynomial) - 1
      q = size(model%ma_polynomial) - 1
      parameters = p + q
      allocate(jacobian(max_lag + 1, parameters))
      jacobian = 0.0_dp
      do column = 1, parameters
         ar = model%ar_polynomial
         ma = model%ma_polynomial
         if (column <= p) then
            step = sqrt(epsilon(1.0_dp))*(1.0_dp + abs(ar(column + 1)))
            ar(column + 1) = ar(column + 1) + step
         else
            step = sqrt(epsilon(1.0_dp))*(1.0_dp + &
               abs(ma(column - p + 1)))
            ma(column - p + 1) = ma(column - p + 1) + step
         end if
         if (integrated) ar = tfarima_polynomial_multiply(ar, &
            model%difference_polynomial)
         if (inverse) then
            trial = tfarima_polynomial_ratio(ar, ma, max_lag)
         else
            trial = tfarima_polynomial_ratio(ma, ar, max_lag)
         end if
         jacobian(:, column) = (trial - out%weight)/step
      end do
      allocate(out%covariance(max_lag + 1, max_lag + 1), &
         out%standard_error(max_lag + 1), &
         out%cumulative_covariance(max_lag + 1, max_lag + 1), &
         out%cumulative_standard_error(max_lag + 1))
      out%covariance = 0.0_dp
      out%cumulative_covariance = 0.0_dp
      if (present(coefficient_covariance)) then
         if (any(shape(coefficient_covariance) /= [parameters, parameters]) .or. &
            .not. all(ieee_is_finite(coefficient_covariance))) then
            out%info = 2
            return
         end if
         out%covariance = matmul(matmul(jacobian, coefficient_covariance), &
            transpose(jacobian))
         out%covariance = 0.5_dp*(out%covariance + transpose(out%covariance))
         cumulative_jacobian = jacobian
         do lag = 2, max_lag + 1
            cumulative_jacobian(lag, :) = cumulative_jacobian(lag - 1, :) + &
               jacobian(lag, :)
         end do
         out%cumulative_covariance = matmul(matmul(cumulative_jacobian, &
            coefficient_covariance), transpose(cumulative_jacobian))
         out%cumulative_covariance = 0.5_dp*(out%cumulative_covariance + &
            transpose(out%cumulative_covariance))
      end if
      do lag = 1, max_lag + 1
         out%standard_error(lag) = sqrt(max(0.0_dp, out%covariance(lag, lag)))
         out%cumulative_standard_error(lag) = sqrt(max(0.0_dp, &
            out%cumulative_covariance(lag, lag)))
      end do
   end function arima_weights

   pure function tfarima_leverrier_faddeev(observation_loading, transition) &
      result(out)
      ! Compute characteristic and observation-weighted adjoint polynomials.
      real(dp), intent(in) :: observation_loading(:), transition(:, :)
      type(tfarima_leverrier_t) :: out
      real(dp), allocatable :: identity(:, :), recursion(:, :), product(:, :)
      integer :: n, order, i

      n = size(transition, 1)
      if (n < 1 .or. size(transition, 2) /= n .or. &
         size(observation_loading) /= n .or. &
         .not. all(ieee_is_finite(transition)) .or. &
         .not. all(ieee_is_finite(observation_loading))) then
         out%info = 1
         return
      end if
      allocate(identity(n, n), recursion(n, n), &
         out%characteristic(n + 1), out%adjoint_loading(n, n))
      identity = 0.0_dp
      do i = 1, n
         identity(i, i) = 1.0_dp
      end do
      recursion = identity
      out%characteristic = 0.0_dp
      out%characteristic(1) = 1.0_dp
      out%adjoint_loading = 0.0_dp
      out%adjoint_loading(1, :) = observation_loading
      do order = 2, n + 1
         product = matmul(transition, recursion)
         out%characteristic(order) = -sum([(product(i, i), i=1,n)])/ &
            real(order - 1, dp)
         recursion = product + &
            out%characteristic(order)*identity
         if (order <= n) out%adjoint_loading(order, :) = &
            matmul(observation_loading, recursion)
      end do
   end function tfarima_leverrier_faddeev

   pure function tfarima_ssm_form(observation_loading, transition, &
      disturbance_covariance, state_noise_contemporaneous) result(out)
      ! Construct a covariance-aware time-invariant univariate SSM form.
      real(dp), intent(in) :: observation_loading(:), transition(:, :)
      real(dp), intent(in) :: disturbance_covariance(:, :)
      logical, intent(in), optional :: state_noise_contemporaneous
      type(tfarima_ssm_form_t) :: out
      integer :: n

      n = size(observation_loading)
      if (n < 1 .or. any(shape(transition) /= [n, n]) .or. &
         any(shape(disturbance_covariance) /= [n + 1, n + 1]) .or. &
         .not. all(ieee_is_finite(observation_loading)) .or. &
         .not. all(ieee_is_finite(transition)) .or. &
         .not. all(ieee_is_finite(disturbance_covariance))) then
         out%info = 1
         return
      end if
      if (maxval(abs(disturbance_covariance - &
         transpose(disturbance_covariance))) > 100.0_dp*epsilon(1.0_dp)) then
         out%info = 2
         return
      end if
      out%observation_loading = observation_loading
      out%transition = transition
      out%disturbance_covariance = 0.5_dp*(disturbance_covariance + &
         transpose(disturbance_covariance))
      if (present(state_noise_contemporaneous)) &
         out%state_noise_contemporaneous = state_noise_contemporaneous
   end function tfarima_ssm_form

   pure function tfarima_switch_ssm_form(form) result(out)
      ! Switch state disturbances between contemporaneous and one-lag forms.
      type(tfarima_ssm_form_t), intent(in) :: form
      type(tfarima_ssm_form_t) :: out
      real(dp), allocatable :: loading(:), covariance(:, :), inverse(:, :)
      real(dp), allocatable :: state_covariance(:, :), cross_covariance(:)
      integer :: n, status

      if (.not. valid_ssm_form(form)) then
         out%info = 1
         return
      end if
      n = size(form%observation_loading)
      loading = form%observation_loading
      covariance = form%disturbance_covariance
      state_covariance = covariance(2:n + 1, 2:n + 1)
      cross_covariance = covariance(1, 2:n + 1)
      if (form%state_noise_contemporaneous) then
         covariance(1, 1) = covariance(1, 1) + &
            dot_product(loading, matmul(state_covariance, loading)) + &
            2.0_dp*dot_product(loading, cross_covariance)
         covariance(1, 2:n + 1) = cross_covariance + &
            matmul(loading, state_covariance)
         covariance(2:n + 1, 1) = covariance(1, 2:n + 1)
         loading = matmul(loading, form%transition)
         out = tfarima_ssm_form(loading, form%transition, covariance, .false.)
      else
         call invert_matrix(form%transition, inverse, status)
         if (status /= 0) then
            out%info = 2
            return
         end if
         loading = matmul(loading, inverse)
         covariance(1, 1) = covariance(1, 1) + &
            dot_product(loading, matmul(state_covariance, loading)) - &
            2.0_dp*dot_product(loading, cross_covariance)
         covariance(1, 2:n + 1) = cross_covariance - &
            matmul(loading, state_covariance)
         covariance(2:n + 1, 1) = covariance(1, 2:n + 1)
         out = tfarima_ssm_form(loading, form%transition, covariance, .true.)
      end if
   end function tfarima_switch_ssm_form

   pure function tfarima_ssm_to_arima(model, tolerance, max_iterations) &
      result(out)
      ! Reduce a time-invariant univariate KFAS model to one ARIMA spectrum.
      type(ssm_model_t), intent(in) :: model
      real(dp), intent(in), optional :: tolerance
      integer, intent(in), optional :: max_iterations
      type(tfarima_ssm_reduction_t) :: out
      type(tfarima_ssm_form_t) :: form
      real(dp), allocatable :: covariance(:, :), process(:, :)
      integer :: n

      if (.not. allocated(model%z) .or. .not. allocated(model%h) .or. &
         .not. allocated(model%transition) .or. .not. allocated(model%r) .or. &
         .not. allocated(model%q)) then
         out%info = 1
         return
      end if
      n = size(model%transition, 1)
      if (n < 1 .or. size(model%z, 1) /= 1 .or. size(model%z, 2) /= n .or. &
         size(model%z, 3) /= 1 .or. any(shape(model%h) /= [1, 1, 1]) .or. &
         any(shape(model%transition) /= [n, n, 1]) .or. &
         size(model%r, 1) /= n .or. size(model%r, 3) /= 1 .or. &
         size(model%q, 1) /= size(model%q, 2) .or. &
         size(model%q, 1) /= size(model%r, 2) .or. size(model%q, 3) /= 1) then
         out%info = 1
         return
      end if
      process = matmul(matmul(model%r(:, :, 1), model%q(:, :, 1)), &
         transpose(model%r(:, :, 1)))
      allocate(covariance(n + 1, n + 1))
      covariance = 0.0_dp
      covariance(1, 1) = model%h(1, 1, 1)
      covariance(2:, 2:) = process
      form = tfarima_ssm_form(model%z(1, :, 1), model%transition(:, :, 1), &
         covariance, .false.)
      if (form%info /= 0) then
         out%info = 2
         return
      end if
      out = tfarima_reduce_ssm(form, tolerance, max_iterations)
   end function tfarima_ssm_to_arima

   pure function tfarima_reduce_ssm(form, tolerance, max_iterations) result(out)
      ! Reduce a covariance-aware state-space form to an ARIMA spectrum.
      type(tfarima_ssm_form_t), intent(in) :: form
      real(dp), intent(in), optional :: tolerance
      integer, intent(in), optional :: max_iterations
      type(tfarima_ssm_reduction_t) :: out
      type(tfarima_leverrier_t) :: leverrier
      type(arima2_roots_t) :: roots
      real(dp), allocatable :: ar(:), difference(:), ma(:)
      real(dp), allocatable :: reconstructed(:)
      logical, allocatable :: unit_root(:)
      real(dp) :: tol, scale
      integer :: limit, n, lag, i, j, coefficient, factor_length

      tol = 1.0e-5_dp
      if (present(tolerance)) tol = tolerance
      limit = 500
      if (present(max_iterations)) limit = max_iterations
      if (.not. valid_ssm_form(form) .or. tol <= 0.0_dp .or. limit < 1) then
         out%info = 1
         return
      end if
      n = size(form%observation_loading)
      leverrier = tfarima_leverrier_faddeev(form%observation_loading, &
         form%transition)
      if (leverrier%info /= 0) then
         out%info = 2
         return
      end if
      out%characteristic = leverrier%characteristic
      roots = arma_polynomial_roots(-leverrier%characteristic(2:))
      if (roots%info /= 0) then
         out%info = 3
         return
      end if
      allocate(unit_root(n))
      unit_root = abs(abs(roots%roots) - 1.0_dp) <= tol
      ar = polynomial_from_selected_roots(roots%roots, .not. unit_root, tol)
      difference = polynomial_from_selected_roots(roots%roots, unit_root, tol)
      if (size(ar) < 1 .or. size(difference) < 1) then
         out%info = 3
         return
      end if

      allocate(out%disturbance_numerator(n + 1, n + 1))
      out%disturbance_numerator = 0.0_dp
      out%disturbance_numerator(1, :) = leverrier%characteristic
      if (form%state_noise_contemporaneous) then
         out%disturbance_numerator(2:, :n) = &
            transpose(leverrier%adjoint_loading)
      else
         out%disturbance_numerator(2:, 2:) = &
            transpose(leverrier%adjoint_loading)
      end if
      out%disturbance_covariance = form%disturbance_covariance
      allocate(out%target_autocovariance(n + 1))
      out%target_autocovariance = 0.0_dp
      do lag = 0, n
         do i = 1, n + 1
            do j = 1, n + 1
               do coefficient = 1, n + 1 - lag
                  out%target_autocovariance(lag + 1) = &
                     out%target_autocovariance(lag + 1) + &
                     form%disturbance_covariance(i, j)* &
                     out%disturbance_numerator(i, coefficient + lag)* &
                     out%disturbance_numerator(j, coefficient)
               end do
            end do
         end do
      end do
      out%factor = tfarima_cramer_wold_factor(out%target_autocovariance, &
         tfarima_factor_best, tolerance=tol, max_iterations=limit)
      if (out%factor%info /= 0 .or. .not. allocated(out%factor%coefficients) .or. &
         abs(out%factor%coefficients(1)) <= tol) then
         out%info = 4
         return
      end if
      factor_length = size(out%factor%coefficients)
      allocate(reconstructed(factor_length))
      call factor_moments(out%factor%coefficients, reconstructed)
      out%reconstructed_autocovariance = reconstructed
      out%reconstruction_error = sum((reconstructed - &
         out%target_autocovariance(:factor_length))**2)
      if (factor_length < size(out%target_autocovariance)) &
         out%reconstruction_error = out%reconstruction_error + &
         sum(out%target_autocovariance(factor_length + 1:)**2)
      out%reconstruction_error = sqrt(out%reconstruction_error)
      scale = out%factor%coefficients(1)
      ma = out%factor%coefficients/scale
      call cancel_reduction_factor(ma, difference, tol)
      call cancel_reduction_factor(ma, ar, tol)
      out%model = tfarima_ucarima_component(ar, difference, ma, scale**2)
      if (out%model%info /= 0) out%info = 5

   contains

      pure function polynomial_from_selected_roots(all_roots, selected, &
         local_tolerance) result(polynomial)
         ! Construct a real normalized polynomial from selected ordinary roots.
         complex(dp), intent(in) :: all_roots(:)
         logical, intent(in) :: selected(:)
         real(dp), intent(in) :: local_tolerance
         real(dp), allocatable :: polynomial(:)
         complex(dp), allocatable :: work(:), next(:)
         integer :: root_index, order, k

         allocate(work(1))
         work = cmplx(1.0_dp, 0.0_dp, dp)
         order = 0
         do root_index = 1, size(all_roots)
            if (.not. selected(root_index)) cycle
            order = order + 1
            allocate(next(order + 1))
            next = cmplx(0.0_dp, 0.0_dp, dp)
            do k = 1, order
               next(k) = next(k) + work(k)
               next(k + 1) = next(k + 1) - work(k)/all_roots(root_index)
            end do
            call move_alloc(next, work)
         end do
         if (maxval(abs(aimag(work))) > 100.0_dp*local_tolerance) then
            allocate(polynomial(0))
         else
            polynomial = real(work, dp)
            polynomial = polynomial/polynomial(1)
         end if
      end function polynomial_from_selected_roots

      pure subroutine cancel_reduction_factor(numerator, denominator, &
         local_tolerance)
         ! Cancel one common polynomial factor in a reduced ARIMA spectrum.
         real(dp), allocatable, intent(inout) :: numerator(:), denominator(:)
         real(dp), intent(in) :: local_tolerance
         real(dp), allocatable :: common(:)
         type(tfarima_polynomial_division_t) :: numerator_division
         type(tfarima_polynomial_division_t) :: denominator_division

         common = tfarima_polynomial_gcd(numerator, denominator, local_tolerance)
         if (size(common) <= 1) return
         numerator_division = tfarima_polynomial_divide(numerator, common, &
            local_tolerance)
         denominator_division = tfarima_polynomial_divide(denominator, common, &
            local_tolerance)
         if (numerator_division%info /= 0 .or. denominator_division%info /= 0) return
         if (maxval(abs(numerator_division%remainder)) > local_tolerance .or. &
            maxval(abs(denominator_division%remainder)) > local_tolerance) return
         numerator = numerator_division%quotient
         denominator = denominator_division%quotient
      end subroutine cancel_reduction_factor

   end function tfarima_reduce_ssm

   pure logical function valid_ssm_form(form) result(valid)
      ! Report whether a covariance-aware state-space form is complete.
      type(tfarima_ssm_form_t), intent(in) :: form
      integer :: n

      valid = form%info == 0 .and. allocated(form%observation_loading) .and. &
         allocated(form%transition) .and. &
         allocated(form%disturbance_covariance)
      if (.not. valid) return
      n = size(form%observation_loading)
      valid = n > 0 .and. all(shape(form%transition) == [n, n]) .and. &
         all(shape(form%disturbance_covariance) == [n + 1, n + 1]) .and. &
         all(ieee_is_finite(form%observation_loading)) .and. &
         all(ieee_is_finite(form%transition)) .and. &
         all(ieee_is_finite(form%disturbance_covariance))
      if (valid) valid = maxval(abs(form%disturbance_covariance - &
         transpose(form%disturbance_covariance))) <= &
         100.0_dp*epsilon(1.0_dp)
   end function valid_ssm_form

   pure function tfarima_root_decomposition_basis(ar_polynomial, &
      difference_polynomial, mean_value, tolerance) result(out)
      ! Build TFARIMA root classifications and deterministic decomposition basis.
      real(dp), intent(in) :: ar_polynomial(:), difference_polynomial(:)
      real(dp), intent(in), optional :: mean_value, tolerance
      type(tfarima_decomposition_basis_t) :: out
      type(arima2_roots_t) :: roots
      complex(dp), allocatable :: inverse_roots(:)
      real(dp), allocatable :: table(:, :), grouped(:, :), temporary(:)
      logical, allocatable :: retained(:)
      real(dp) :: mu, tol, modulus, frequency, coefficient, argument
      integer :: p, d, total, i, j, h, order, columns, column, group_count

      mu = 0.0_dp
      if (present(mean_value)) mu = mean_value
      tol = 1.0e-5_dp
      if (present(tolerance)) tol = tolerance
      p = size(ar_polynomial) - 1
      d = size(difference_polynomial) - 1
      total = p + d
      if (p < 0 .or. d < 0 .or. total < 1 .or. tol <= 0.0_dp .or. &
         abs(ar_polynomial(1)) <= tiny(1.0_dp) .or. &
         abs(difference_polynomial(1)) <= tiny(1.0_dp)) then
         out%info = 1
         return
      end if
      allocate(inverse_roots(total))
      if (p > 0) then
         roots = decomposition_polynomial_roots(ar_polynomial/ar_polynomial(1))
         if (roots%info /= 0) then
            out%info = 2
            return
         end if
         inverse_roots(:p) = 1.0_dp/roots%roots
      end if
      if (d > 0) then
         roots = decomposition_polynomial_roots( &
            difference_polynomial/difference_polynomial(1))
         if (roots%info /= 0) then
            out%info = 2
            return
         end if
         inverse_roots(p + 1:) = 1.0_dp/roots%roots
      end if
      allocate(table(total, 6))
      do i = 1, total
         table(i, 1) = real(inverse_roots(i), dp)
         table(i, 2) = aimag(inverse_roots(i))
         modulus = abs(inverse_roots(i))
         table(i, 3) = modulus
         if (modulus <= tol) then
            frequency = 0.0_dp
         else
            frequency = acos(max(-1.0_dp, min(1.0_dp, &
               table(i, 1)/modulus)))/(2.0_dp*acos(-1.0_dp))
         end if
         if (abs(table(i, 2)**2) <= tol) then
            table(i, 2) = 0.0_dp
            if (table(i, 1) >= 0.0_dp) then
               frequency = 0.0_dp
            else
               frequency = 0.5_dp
            end if
         end if
         table(i, 4) = frequency
         if (frequency <= tol) then
            table(i, 5) = huge(1.0_dp)
         else
            table(i, 5) = 1.0_dp/frequency
         end if
         table(i, 6) = 1.0_dp
      end do
      do i = 1, total - 1
         do j = i + 1, total
            if (table(j, 3) < table(i, 3) - tol .or. &
               (abs(table(j, 3) - table(i, 3)) <= tol .and. &
               table(j, 4) < table(i, 4))) then
               temporary = table(i, :)
               table(i, :) = table(j, :)
               table(j, :) = temporary
            end if
         end do
      end do
      allocate(retained(total))
      retained = .true.
      do i = 1, total
         if (.not. retained(i)) cycle
         do j = i + 1, total
            if (abs(table(j, 4) - table(i, 4)) > tol) cycle
            if (abs(table(j, 1) - table(i, 1)) <= tol .and. &
               abs(table(j, 2) - table(i, 2)) <= tol) then
               table(i, 6) = table(i, 6) + 1.0_dp
               retained(j) = .false.
            end if
         end do
      end do
      group_count = count(retained)
      allocate(grouped(group_count, 6))
      j = 0
      do i = 1, total
         if (.not. retained(i)) cycle
         j = j + 1
         grouped(j, :) = table(i, :)
      end do
      out%root_table = grouped
      columns = total
      if (mu /= 0.0_dp) columns = columns + 1
      allocate(out%classification(4, columns), out%basis(columns + 1, columns))
      out%classification = 0.0_dp
      out%basis = 0.0_dp
      column = 1
      do j = 1, group_count
         frequency = grouped(j, 4)
         order = nint(grouped(j, 6))
         if (abs(frequency) <= tol) then
            coefficient = grouped(j, 1)
            if (abs(coefficient - 1.0_dp) <= tol) then
               if (mu /= 0.0_dp) order = order + 1
               out%classification(tfarima_component_trend, &
                  column:column + order - 1) = 1.0_dp
            else
               out%classification(tfarima_component_exponential, &
                  column:column + order - 1) = 1.0_dp
            end if
            do h = 0, order - 1
               do i = 1, columns + 1
                  out%basis(i, column) = coefficient**i*real(i, dp)**h
               end do
               column = column + 1
            end do
         else
            coefficient = grouped(j, 3)
            if (abs(frequency - 0.5_dp) <= tol) then
               if (abs(coefficient - 1.0_dp) <= tol) then
                  out%classification(tfarima_component_seasonal, &
                     column:column + order - 1) = 1.0_dp
               else
                  out%classification(tfarima_component_cycle, &
                     column:column + order - 1) = 1.0_dp
               end if
               coefficient = grouped(j, 1)
               do h = 0, order - 1
                  do i = 1, columns + 1
                     out%basis(i, column) = coefficient**i*real(i, dp)**h
                  end do
                  column = column + 1
               end do
            else if (grouped(j, 2) > 0.0_dp) then
               if (abs(coefficient - 1.0_dp) <= tol) then
                  out%classification(tfarima_component_seasonal, &
                     column:column + 2*order - 1) = 1.0_dp
               else
                  out%classification(tfarima_component_cycle, &
                     column:column + 2*order - 1) = 1.0_dp
               end if
               argument = 2.0_dp*acos(-1.0_dp)*frequency
               do h = 0, order - 1
                  do i = 1, columns + 1
                     out%basis(i, column) = coefficient**i*cos(argument*i)* &
                        real(i, dp)**h
                  end do
                  column = column + 1
               end do
               do h = 0, order - 1
                  do i = 1, columns + 1
                     out%basis(i, column) = coefficient**i*sin(argument*i)* &
                        real(i, dp)**h
                  end do
                  column = column + 1
               end do
            end if
         end if
      end do
      if (column /= columns + 1) out%info = 3
   end function tfarima_root_decomposition_basis

   pure function tfarima_root_decompose(series, ar_polynomial, &
      difference_polynomial, ma_polynomial, innovation_variance, mean_value, &
      method, log_transform, tolerance) result(out)
      ! Decompose an ARIMA series by forecast, backcast, or mixed root effects.
      real(dp), intent(in) :: series(:), ar_polynomial(:)
      real(dp), intent(in) :: difference_polynomial(:), ma_polynomial(:)
      real(dp), intent(in) :: innovation_variance
      real(dp), intent(in), optional :: mean_value, tolerance
      integer, intent(in), optional :: method
      logical, intent(in), optional :: log_transform
      type(tfarima_root_decomposition_t) :: out
      real(dp), allocatable :: transformed(:), effects(:, :), loading(:)
      real(dp) :: mu, tol
      integer :: selected_method, component
      logical :: logarithm

      mu = 0.0_dp
      if (present(mean_value)) mu = mean_value
      tol = 1.0e-5_dp
      if (present(tolerance)) tol = tolerance
      selected_method = tfarima_decomposition_mixed
      if (present(method)) selected_method = method
      logarithm = .false.
      if (present(log_transform)) logarithm = log_transform
      if (size(series) < 2 .or. innovation_variance < 0.0_dp .or. &
         .not. all(ieee_is_finite(series)) .or. &
         (logarithm .and. any(series <= 0.0_dp)) .or. &
         selected_method < tfarima_decomposition_forecast .or. &
         selected_method > tfarima_decomposition_mixed) then
         out%info = 1
         return
      end if
      out%method = selected_method
      out%decomposition_basis = tfarima_root_decomposition_basis( &
         ar_polynomial, difference_polynomial, mu, tol)
      if (out%decomposition_basis%info /= 0) then
         out%info = 2
         return
      end if
      transformed = series
      if (logarithm) transformed = log(series)
      effects = root_decomposition_effects(transformed, mu, ar_polynomial, &
         difference_polynomial, ma_polynomial, out%decomposition_basis%basis, &
         selected_method)
      if (size(effects, 1) /= size(series)) then
         out%info = 3
         return
      end if
      out%state_coefficients = effects(:, :size(effects, 2) - 1)
      out%innovations = effects(:, size(effects, 2))
      allocate(out%component(size(series), 4))
      do component = 1, 4
         loading = out%decomposition_basis%classification(component, :)* &
            out%decomposition_basis%basis(1, :)
         out%component(:, component) = matmul(out%state_coefficients, loading)
      end do
      out%irregular = out%innovations
      out%reconstruction = sum(out%component, dim=2) + out%irregular
      if (logarithm) then
         out%seasonally_adjusted = exp(transformed - &
            out%component(:, tfarima_component_seasonal))
      else
         out%seasonally_adjusted = series - &
            out%component(:, tfarima_component_seasonal)
      end if
   end function tfarima_root_decompose

   pure function tfarima_arima_to_structural_ssm(model, mean_value, &
      multiple_sources, contemporaneous, grouping, tolerance, max_iterations) &
      result(out)
      ! Convert an ARIMA model to TFARIMA's eventual-forecast structural form.
      type(tfarima_ucarima_component_t), intent(in) :: model
      real(dp), intent(in), optional :: mean_value, grouping(:, :), tolerance
      logical, intent(in), optional :: multiple_sources, contemporaneous
      integer, intent(in), optional :: max_iterations
      type(tfarima_structural_ssm_t) :: out
      real(dp), allocatable :: first_basis(:, :), next_basis(:, :), inverse(:, :)
      real(dp), allocatable :: transition(:, :), loading(:), forcing(:), psi(:)
      real(dp), allocatable :: covariance(:, :), design(:, :), mapping(:, :)
      real(dp), allocatable :: effective_design(:, :), parameters(:), unrestricted(:)
      real(dp) :: mu, tol, kappa
      integer :: r, i, status, limit
      logical :: use_multiple, use_contemporaneous

      mu = 0.0_dp
      if (present(mean_value)) mu = mean_value
      tol = 1.0e-8_dp
      if (present(tolerance)) tol = tolerance
      limit = 10000
      if (present(max_iterations)) limit = max_iterations
      use_multiple = .true.
      if (present(multiple_sources)) use_multiple = multiple_sources
      use_contemporaneous = .true.
      if (present(contemporaneous)) use_contemporaneous = contemporaneous
      out%multiple_sources = use_multiple
      if (model%info /= 0 .or. .not. allocated(model%ar_polynomial) .or. &
         .not. allocated(model%difference_polynomial) .or. &
         .not. allocated(model%ma_polynomial) .or. &
         model%innovation_variance < 0.0_dp .or. tol <= 0.0_dp .or. limit < 1) then
         out%info = 1
         return
      end if
      if (.not. all(ieee_is_finite(model%ar_polynomial)) .or. &
         .not. all(ieee_is_finite(model%difference_polynomial)) .or. &
         .not. all(ieee_is_finite(model%ma_polynomial)) .or. &
         .not. ieee_is_finite(model%innovation_variance)) then
         out%info = 1
         return
      end if
      out%decomposition_basis = tfarima_root_decomposition_basis( &
         model%ar_polynomial, model%difference_polynomial, mu, tol)
      if (out%decomposition_basis%info /= 0) then
         out%info = 2
         return
      end if
      r = size(out%decomposition_basis%basis, 2)
      first_basis = out%decomposition_basis%basis(:r, :)
      next_basis = out%decomposition_basis%basis(2:r + 1, :)
      call invert_matrix(first_basis, inverse, status)
      if (status /= 0) then
         out%info = 3
         return
      end if
      transition = matmul(inverse, next_basis)
      psi = tfarima_polynomial_ratio(model%ma_polynomial, &
         tfarima_polynomial_multiply(model%ar_polynomial, &
         model%difference_polynomial), r)
      if (size(psi) /= r + 1) then
         out%info = 3
         return
      end if
      forcing = matmul(inverse, psi(2:r + 1))
      if (use_contemporaneous) then
         call invert_matrix(transition, inverse, status)
         if (status /= 0) then
            out%info = 3
            return
         end if
         loading = matmul(out%decomposition_basis%basis(1, :), inverse)
         kappa = 1.0_dp - dot_product(loading, forcing)
      else
         loading = out%decomposition_basis%basis(1, :)
         kappa = 1.0_dp
      end if
      out%target_autocovariance = finite_ma_autocovariance( &
         model%ma_polynomial, model%innovation_variance, r)
      allocate(covariance(r + 1, r + 1))
      covariance = 0.0_dp
      if (.not. use_multiple) then
         forcing = [kappa, forcing]
         covariance = model%innovation_variance* &
            spread(forcing, 2, r + 1)*spread(forcing, 1, r + 1)
         out%disturbance_variance = [(covariance(i, i), i=1,r + 1)]
         allocate(out%variance_design(r + 1, 0))
         out%admissible = .true.
      else
         out%leverrier = tfarima_leverrier_faddeev(loading, transition)
         if (out%leverrier%info /= 0) then
            out%info = 4
            return
         end if
         allocate(design(r + 1, r + 1))
         design(:, 1) = finite_ma_autocovariance( &
            out%leverrier%characteristic, 1.0_dp, r)
         do i = 1, r
            design(:, i + 1) = finite_ma_autocovariance( &
               out%leverrier%adjoint_loading(:, i), 1.0_dp, r)
         end do
         out%variance_design = design
         if (present(grouping)) then
            if (size(grouping, 1) /= r + 1 .or. size(grouping, 2) < 1 .or. &
               .not. all(ieee_is_finite(grouping)) .or. any(grouping < 0.0_dp)) then
               out%info = 5
               return
            end if
            mapping = grouping
         else
            allocate(mapping(r + 1, r + 1))
            mapping = 0.0_dp
            do i = 1, r + 1
               mapping(i, i) = 1.0_dp
            end do
         end if
         effective_design = matmul(design, mapping)
         unrestricted = least_squares_coefficients(effective_design, &
            out%target_autocovariance, tol)
         out%admissible = size(unrestricted) == size(mapping, 2)
         if (out%admissible) out%admissible = all(ieee_is_finite(unrestricted))
         if (out%admissible) out%admissible = &
            all(matmul(mapping, unrestricted) >= -tol)
         if (out%admissible) then
            parameters = unrestricted
         else
            parameters = nonnegative_least_squares(effective_design, &
               out%target_autocovariance, limit, tol)
            out%used_nonnegative_fit = .true.
         end if
         if (size(parameters) /= size(mapping, 2)) then
            out%info = 6
            return
         end if
         out%disturbance_variance = max(0.0_dp, matmul(mapping, parameters))
         do i = 1, r + 1
            covariance(i, i) = out%disturbance_variance(i)
         end do
      end if
      out%form = tfarima_ssm_form(loading, transition, covariance, &
         use_contemporaneous)
      if (out%form%info /= 0) then
         out%info = 7
         return
      end if
      out%fitted_autocovariance = structural_ssm_autocovariance(out%form)
      if (size(out%fitted_autocovariance) /= r + 1) then
         out%info = 8
         return
      end if
      out%reconstruction_error = sqrt(sum((out%fitted_autocovariance - &
         out%target_autocovariance)**2))
   end function tfarima_arima_to_structural_ssm

   pure function finite_ma_autocovariance(polynomial, variance, max_lag) &
      result(covariance)
      ! Return finite-MA autocovariances through a requested lag.
      real(dp), intent(in) :: polynomial(:), variance
      integer, intent(in) :: max_lag
      real(dp), allocatable :: covariance(:)
      integer :: lag, n

      if (size(polynomial) < 1 .or. max_lag < 0 .or. variance < 0.0_dp) then
         allocate(covariance(0))
         return
      end if
      n = size(polynomial)
      allocate(covariance(max_lag + 1))
      covariance = 0.0_dp
      do lag = 0, min(max_lag, n - 1)
         covariance(lag + 1) = variance*dot_product( &
            polynomial(:n - lag), polynomial(lag + 1:))
      end do
   end function finite_ma_autocovariance

   pure function least_squares_coefficients(matrix, rhs, tolerance) &
      result(coefficients)
      ! Solve a full-rank least-squares system through normal equations.
      real(dp), intent(in) :: matrix(:, :), rhs(:), tolerance
      real(dp), allocatable :: coefficients(:), normal(:, :), inverse(:, :)
      integer :: status

      if (size(matrix, 1) /= size(rhs) .or. size(matrix, 2) < 1) then
         allocate(coefficients(0))
         return
      end if
      normal = matmul(transpose(matrix), matrix)
      call invert_matrix(normal, inverse, status)
      if (status /= 0) then
         normal = normal + tolerance*identity_array(size(normal, 1))
         call invert_matrix(normal, inverse, status)
      end if
      if (status /= 0) then
         allocate(coefficients(0))
         return
      end if
      coefficients = matmul(inverse, matmul(transpose(matrix), rhs))
   end function least_squares_coefficients

   pure function nonnegative_least_squares(matrix, rhs, max_iterations, &
      tolerance) result(coefficients)
      ! Minimize a linear least-squares residual under nonnegative coordinates.
      real(dp), intent(in) :: matrix(:, :), rhs(:), tolerance
      integer, intent(in) :: max_iterations
      real(dp), allocatable :: coefficients(:), residual(:)
      real(dp) :: denominator, updated, change
      integer :: iteration, column

      allocate(coefficients(size(matrix, 2)), residual(size(rhs)))
      coefficients = 0.0_dp
      residual = rhs
      do iteration = 1, max_iterations
         change = 0.0_dp
         do column = 1, size(matrix, 2)
            residual = residual + matrix(:, column)*coefficients(column)
            denominator = dot_product(matrix(:, column), matrix(:, column))
            if (denominator > tiny(1.0_dp)) then
               updated = max(0.0_dp, dot_product(matrix(:, column), residual)/ &
                  denominator)
            else
               updated = 0.0_dp
            end if
            change = max(change, abs(updated - coefficients(column)))
            coefficients(column) = updated
            residual = residual - matrix(:, column)*coefficients(column)
         end do
         if (change <= tolerance*(1.0_dp + maxval(coefficients))) exit
      end do
   end function nonnegative_least_squares

   pure function structural_ssm_autocovariance(form) result(covariance)
      ! Reconstruct reduced numerator autocovariances from an SSM disturbance form.
      type(tfarima_ssm_form_t), intent(in) :: form
      real(dp), allocatable :: covariance(:), numerator(:, :)
      type(tfarima_leverrier_t) :: leverrier
      integer :: r, lag, i, j, coefficient

      if (.not. valid_ssm_form(form)) then
         allocate(covariance(0))
         return
      end if
      r = size(form%observation_loading)
      leverrier = tfarima_leverrier_faddeev(form%observation_loading, &
         form%transition)
      if (leverrier%info /= 0) then
         allocate(covariance(0))
         return
      end if
      allocate(numerator(r + 1, r + 1), covariance(r + 1))
      numerator = 0.0_dp
      numerator(1, :) = leverrier%characteristic
      if (form%state_noise_contemporaneous) then
         numerator(2:, :r) = transpose(leverrier%adjoint_loading)
      else
         numerator(2:, 2:) = transpose(leverrier%adjoint_loading)
      end if
      covariance = 0.0_dp
      do lag = 0, r
         do i = 1, r + 1
            do j = 1, r + 1
               do coefficient = 1, r + 1 - lag
                  covariance(lag + 1) = covariance(lag + 1) + &
                     form%disturbance_covariance(i, j)* &
                     numerator(i, coefficient + lag)*numerator(j, coefficient)
               end do
            end do
         end do
      end do
   end function structural_ssm_autocovariance

   pure function identity_array(order) result(identity)
      ! Return a square identity matrix of requested order.
      integer, intent(in) :: order
      real(dp) :: identity(order, order)
      integer :: i

      identity = 0.0_dp
      do i = 1, order
         identity(i, i) = 1.0_dp
      end do
   end function identity_array

   pure function tfarima_structural_initialize(series, form, observations, &
      log_transform, regressors, regression_coefficients) result(out)
      ! Estimate structural initial moments by filtering a GLS design matrix.
      real(dp), intent(in) :: series(:)
      type(tfarima_ssm_form_t), intent(in) :: form
      integer, intent(in), optional :: observations
      logical, intent(in), optional :: log_transform
      real(dp), intent(in), optional :: regressors(:, :)
      real(dp), intent(in), optional :: regression_coefficients(:)
      type(tfarima_structural_initialization_t) :: out
      type(tfarima_ssm_form_t) :: lagged
      real(dp), allocatable :: adjusted(:)
      integer :: count

      adjusted = prepare_structural_series(series, log_transform, regressors, &
         regression_coefficients)
      if (size(adjusted) /= size(series) .or. .not. valid_ssm_form(form)) then
         out%info = 1
         return
      end if
      lagged = form
      if (form%state_noise_contemporaneous) lagged = tfarima_switch_ssm_form(form)
      if (lagged%info /= 0) then
         out%info = 2
         return
      end if
      count = size(series)
      if (present(observations)) count = observations
      out = initialize_structural_series(adjusted, lagged, count)
   end function tfarima_structural_initialize

   pure function tfarima_structural_filter(series, form, initial_state, &
      initial_covariance, log_transform, regressors, regression_coefficients) &
      result(out)
      ! Filter a structural model with correlated observation and state noise.
      real(dp), intent(in) :: series(:)
      type(tfarima_ssm_form_t), intent(in) :: form
      real(dp), intent(in), optional :: initial_state(:)
      real(dp), intent(in), optional :: initial_covariance(:, :)
      logical, intent(in), optional :: log_transform
      real(dp), intent(in), optional :: regressors(:, :)
      real(dp), intent(in), optional :: regression_coefficients(:)
      type(tfarima_structural_filter_t) :: out
      type(tfarima_structural_initialization_t) :: initialization
      type(tfarima_ssm_form_t) :: lagged
      real(dp), allocatable :: adjusted(:), state(:), covariance(:, :)
      integer :: dimension

      adjusted = prepare_structural_series(series, log_transform, regressors, &
         regression_coefficients)
      if (size(adjusted) /= size(series) .or. .not. valid_ssm_form(form)) then
         out%info = 1
         return
      end if
      lagged = form
      if (form%state_noise_contemporaneous) lagged = tfarima_switch_ssm_form(form)
      if (lagged%info /= 0) then
         out%info = 2
         return
      end if
      dimension = size(lagged%observation_loading)
      if (present(initial_state) .neqv. present(initial_covariance)) then
         out%info = 3
         return
      end if
      if (present(initial_state)) then
         if (size(initial_state) /= dimension .or. &
            any(shape(initial_covariance) /= [dimension, dimension])) then
            out%info = 3
            return
         end if
         state = initial_state
         covariance = initial_covariance
      else
         initialization = initialize_structural_series(adjusted, lagged, &
            size(adjusted))
         if (initialization%info /= 0) then
            out%info = 4
            return
         end if
         state = initialization%state
         covariance = initialization%covariance
      end if
      out = run_structural_filter(adjusted, lagged, state, covariance)
   end function tfarima_structural_filter

   pure function tfarima_structural_smooth(series, form, initial_state, &
      initial_covariance, log_transform, regressors, regression_coefficients) &
      result(out)
      ! Smooth observation-aligned states in a joint-disturbance structural model.
      real(dp), intent(in) :: series(:)
      type(tfarima_ssm_form_t), intent(in) :: form
      real(dp), intent(in), optional :: initial_state(:)
      real(dp), intent(in), optional :: initial_covariance(:, :)
      logical, intent(in), optional :: log_transform
      real(dp), intent(in), optional :: regressors(:, :)
      real(dp), intent(in), optional :: regression_coefficients(:)
      type(tfarima_structural_smoother_t) :: out
      type(tfarima_structural_initialization_t) :: initialization
      type(tfarima_ssm_form_t) :: lagged
      real(dp), allocatable :: adjusted(:), state(:), covariance(:, :)
      real(dp), allocatable :: rvector(:), information(:, :), propagation(:, :)
      real(dp), allocatable :: next_rvector(:), next_information(:, :), prior(:, :)
      integer :: dimension, time

      adjusted = prepare_structural_series(series, log_transform, regressors, &
         regression_coefficients)
      if (size(adjusted) /= size(series) .or. .not. valid_ssm_form(form)) then
         out%info = 1
         return
      end if
      lagged = form
      if (form%state_noise_contemporaneous) lagged = tfarima_switch_ssm_form(form)
      if (lagged%info /= 0) then
         out%info = 2
         return
      end if
      dimension = size(lagged%observation_loading)
      if (present(initial_state) .neqv. present(initial_covariance)) then
         out%info = 3
         return
      end if
      if (present(initial_state)) then
         if (size(initial_state) /= dimension .or. &
            any(shape(initial_covariance) /= [dimension, dimension])) then
            out%info = 3
            return
         end if
         state = initial_state
         covariance = initial_covariance
      else
         initialization = initialize_structural_series(adjusted, lagged, &
            size(adjusted))
         if (initialization%info /= 0) then
            out%info = 4
            return
         end if
         state = initialization%state
         covariance = initialization%covariance
      end if
      out%filter = run_structural_filter(adjusted, lagged, state, covariance)
      if (out%filter%info /= 0) then
         out%info = 5
         return
      end if
      allocate(out%state(dimension, size(series)), &
         out%covariance(dimension, dimension, size(series)))
      allocate(rvector(dimension), information(dimension, dimension))
      rvector = 0.0_dp
      information = 0.0_dp
      do time = size(series), 1, -1
         propagation = lagged%transition - spread(out%filter%gain(:, time), &
            2, dimension)*spread(lagged%observation_loading, 1, dimension)
         next_rvector = lagged%observation_loading* &
            out%filter%innovation(time)/out%filter%innovation_variance(time) + &
            matmul(transpose(propagation), rvector)
         next_information = spread(lagged%observation_loading, 2, dimension)* &
            spread(lagged%observation_loading, 1, dimension)/ &
            out%filter%innovation_variance(time) + &
            matmul(matmul(transpose(propagation), information), propagation)
         prior = out%filter%observation_covariance(:, :, time)
         out%state(:, time) = out%filter%observation_state(:, time) + &
            matmul(prior, next_rvector)
         out%covariance(:, :, time) = prior - &
            matmul(matmul(prior, next_information), prior)
         out%covariance(:, :, time) = 0.5_dp*(out%covariance(:, :, time) + &
            transpose(out%covariance(:, :, time)))
         rvector = next_rvector
         information = next_information
      end do
   end function tfarima_structural_smooth

   pure function tfarima_structural_forecast(filtered, form, horizon, &
      future_regressors, regression_coefficients, log_transform) result(out)
      ! Forecast a filtered structural model with optional regression and log scale.
      type(tfarima_structural_filter_t), intent(in) :: filtered
      type(tfarima_ssm_form_t), intent(in) :: form
      integer, intent(in) :: horizon
      real(dp), intent(in), optional :: future_regressors(:, :)
      real(dp), intent(in), optional :: regression_coefficients(:)
      logical, intent(in), optional :: log_transform
      type(tfarima_structural_forecast_t) :: out
      type(tfarima_ssm_form_t) :: lagged
      real(dp), allocatable :: state(:), covariance(:, :), process(:, :)
      real(dp) :: regression, latent_variance
      integer :: dimension, step
      logical :: logarithm

      logarithm = .false.
      if (present(log_transform)) logarithm = log_transform
      if (filtered%info /= 0 .or. horizon < 1 .or. &
         .not. allocated(filtered%filtered_state) .or. &
         .not. valid_ssm_form(form)) then
         out%info = 1
         return
      end if
      lagged = form
      if (form%state_noise_contemporaneous) lagged = tfarima_switch_ssm_form(form)
      if (lagged%info /= 0) then
         out%info = 2
         return
      end if
      if (present(future_regressors) .neqv. present(regression_coefficients)) then
         out%info = 3
         return
      end if
      if (present(future_regressors)) then
         if (size(future_regressors, 1) /= horizon .or. &
            size(future_regressors, 2) /= size(regression_coefficients)) then
            out%info = 3
            return
         end if
      end if
      dimension = size(lagged%observation_loading)
      state = filtered%filtered_state(:, size(filtered%filtered_state, 2))
      covariance = filtered%filtered_covariance(:, :, &
         size(filtered%filtered_covariance, 3))
      process = lagged%disturbance_covariance(2:, 2:)
      allocate(out%latent_mean(horizon), out%latent_variance(horizon), &
         out%mean(horizon), out%variance(horizon), &
         out%state(dimension, horizon), &
         out%state_covariance(dimension, dimension, horizon))
      do step = 1, horizon
         regression = 0.0_dp
         if (present(future_regressors)) regression = &
            dot_product(future_regressors(step, :), regression_coefficients)
         out%latent_mean(step) = dot_product(lagged%observation_loading, state) + &
            regression
         latent_variance = dot_product(lagged%observation_loading, &
            matmul(covariance, lagged%observation_loading)) + &
            lagged%disturbance_covariance(1, 1)
         out%latent_variance(step) = max(0.0_dp, latent_variance)
         if (logarithm) then
            out%mean(step) = exp(out%latent_mean(step) + &
               0.5_dp*out%latent_variance(step))
            out%variance(step) = (exp(out%latent_variance(step)) - 1.0_dp)* &
               exp(2.0_dp*out%latent_mean(step) + out%latent_variance(step))
         else
            out%mean(step) = out%latent_mean(step)
            out%variance(step) = out%latent_variance(step)
         end if
         out%state(:, step) = state
         out%state_covariance(:, :, step) = covariance
         state = matmul(lagged%transition, state)
         covariance = matmul(matmul(lagged%transition, covariance), &
            transpose(lagged%transition)) + process
         covariance = 0.5_dp*(covariance + transpose(covariance))
      end do
   end function tfarima_structural_forecast

   pure function tfarima_band_cholesky(covariance_band) result(out)
      ! Factor a symmetric positive-definite matrix in lower-band storage.
      real(dp), intent(in) :: covariance_band(:, :)
      type(tfarima_band_cholesky_t) :: out
      real(dp) :: value
      integer :: row, lag, k, observations, width

      observations = size(covariance_band, 1)
      width = size(covariance_band, 2)
      if (observations < 1 .or. width < 1 .or. width > observations .or. &
         .not. all(ieee_is_finite(covariance_band))) then
         out%info = 1
         return
      end if
      allocate(out%factor(observations, width))
      out%factor = 0.0_dp
      do row = 1, observations
         value = covariance_band(row, 1)
         if (width > 1) value = value - &
            dot_product(out%factor(row, 2:), out%factor(row, 2:))
         if (value <= 100.0_dp*epsilon(1.0_dp)) then
            out%info = row
            return
         end if
         out%factor(row, 1) = sqrt(value)
         out%log_determinant = out%log_determinant + &
            log(out%factor(row, 1))
         do lag = 1, width - 1
            if (row + lag > observations) exit
            value = covariance_band(row + lag, lag + 1)
            do k = 1, width - lag - 1
               value = value - out%factor(row + lag, lag + k + 1)* &
                  out%factor(row, k + 1)
            end do
            out%factor(row + lag, lag + 1) = value/out%factor(row, 1)
         end do
      end do
   end function tfarima_band_cholesky

   pure function tfarima_band_forward_solve(cholesky_band, right_hand_side) &
      result(solution)
      ! Solve L*x=b for a lower triangular matrix in band storage.
      real(dp), intent(in) :: cholesky_band(:, :), right_hand_side(:)
      real(dp), allocatable :: solution(:)
      real(dp) :: value
      integer :: row, column, observations, width

      observations = size(cholesky_band, 1)
      width = size(cholesky_band, 2)
      if (observations < 1 .or. width < 1 .or. width > observations .or. &
         size(right_hand_side) /= observations .or. &
         any(cholesky_band(:, 1) <= 0.0_dp) .or. &
         .not. all(ieee_is_finite(cholesky_band)) .or. &
         .not. all(ieee_is_finite(right_hand_side))) then
         allocate(solution(0))
         return
      end if
      allocate(solution(observations))
      do row = 1, observations
         value = right_hand_side(row)
         do column = max(1, row - width + 1), row - 1
            value = value - cholesky_band(row, row - column + 1)* &
               solution(column)
         end do
         solution(row) = value/cholesky_band(row, 1)
      end do
   end function tfarima_band_forward_solve

   pure function tfarima_reduced_likelihood(series, ar_polynomial, &
      disturbance_numerator, disturbance_covariance) result(out)
      ! Evaluate TFARIMA's exact profiled reduced-form Gaussian likelihood.
      real(dp), intent(in) :: series(:), ar_polynomial(:)
      real(dp), intent(in) :: disturbance_numerator(:, :)
      real(dp), intent(in) :: disturbance_covariance(:, :)
      type(tfarima_reduced_likelihood_t) :: out
      type(tfarima_band_cholesky_t) :: factorization
      real(dp), allocatable :: ar(:), numerator(:, :), psi(:, :)
      real(dp), allocatable :: phi0(:, :), g0_matrix(:, :), g0_inverse(:, :)
      real(dp), allocatable :: gamma_rhs(:), gamma_initial(:), toeplitz(:, :)
      real(dp), allocatable :: theta0(:, :), psi0(:, :), block_covariance(:, :)
      real(dp), allocatable :: cross(:, :), correction(:, :), ma_covariance(:)
      real(dp), allocatable :: conditional(:), whitened(:)
      real(dp) :: scale, value, sum_squares
      integer :: p, q, r, disturbances, observations
      integer :: i, j, h, lag, time, status

      observations = size(series)
      p = size(ar_polynomial) - 1
      q = size(disturbance_numerator, 2) - 1
      disturbances = size(disturbance_numerator, 1)
      r = max(p, q + 1)
      if (observations < max(1, r) .or. p < 0 .or. q < 0 .or. &
         disturbances < 1 .or. &
         any(shape(disturbance_covariance) /= [disturbances, disturbances]) .or. &
         abs(ar_polynomial(1)) <= tiny(1.0_dp) .or. &
         .not. all(ieee_is_finite(series)) .or. &
         .not. all(ieee_is_finite(ar_polynomial)) .or. &
         .not. all(ieee_is_finite(disturbance_numerator)) .or. &
         .not. all(ieee_is_finite(disturbance_covariance)) .or. &
         maxval(abs(disturbance_covariance - &
         transpose(disturbance_covariance))) > 100.0_dp*epsilon(1.0_dp)) then
         out%info = 1
         return
      end if
      ar = ar_polynomial/ar_polynomial(1)
      numerator = disturbance_numerator/ar_polynomial(1)
      conditional = conditional_residuals(series, ar, [1.0_dp])
      allocate(out%covariance_band(observations, r))
      out%covariance_band = 0.0_dp

      if (p > 0) then
         psi = numerator
         do i = 2, q
            do j = 1, min(p, i - 1)
               psi(:, i) = psi(:, i) - ar(j + 1)*psi(:, i - j)
            end do
         end do
         allocate(phi0(p, p), g0_matrix(p + 1, p + 1), &
            gamma_rhs(p + 1))
         phi0 = 0.0_dp
         g0_matrix = 0.0_dp
         gamma_rhs = 0.0_dp
         do i = 0, p - 1
            do j = i, p - 1
               phi0(i + 1, j - i + 1) = -ar(j + 2)
            end do
         end do
         do i = 0, p
            do j = 0, i
               g0_matrix(i + 1, j + 1) = ar(i - j + 1)
            end do
            do j = i + 1, p
               g0_matrix(i + 1, j - i + 1) = &
                  g0_matrix(i + 1, j - i + 1) + ar(j + 1)
            end do
            do j = i, q
               gamma_rhs(i + 1) = gamma_rhs(i + 1) + &
                  dot_product(matmul(disturbance_covariance, &
                  numerator(:, j - i + 1)), numerator(:, j + 1))
            end do
         end do
         call invert_matrix(g0_matrix, g0_inverse, status)
         if (status /= 0) then
            out%info = 2
            return
         end if
         gamma_initial = matmul(g0_inverse, gamma_rhs)
         allocate(toeplitz(p, p))
         do i = 1, p
            do j = 1, p
               toeplitz(i, j) = gamma_initial(abs(i - j) + 1)
            end do
         end do
         toeplitz = matmul(matmul(phi0, toeplitz), transpose(phi0))
         do i = 0, p - 1
            do j = 0, i
               out%covariance_band(i + 1, i - j + 1) = toeplitz(i + 1, j + 1)
            end do
         end do

         allocate(theta0(q, q*disturbances), psi0(p, q*disturbances))
         theta0 = 0.0_dp
         psi0 = 0.0_dp
         do i = 0, q - 1
            do j = i, q - 1
               do h = 0, disturbances - 1
                  theta0(i + 1, (j - i)*disturbances + h + 1) = &
                     numerator(h + 1, j + 2)
               end do
            end do
         end do
         do i = 0, p - 1
            do j = i, q - 1
               do h = 0, disturbances - 1
                  psi0(i + 1, j*disturbances + h + 1) = &
                     psi(h + 1, j - i + 1)
               end do
            end do
         end do
         block_covariance = repeated_block_diagonal(disturbance_covariance, q)
         cross = matmul(matmul(matmul(phi0, psi0), block_covariance), &
            transpose(theta0))
         allocate(correction(r, r))
         correction = 0.0_dp
         do i = 1, p
            do j = 1, q
               correction(i, j) = correction(i, j) + cross(i, j)
               correction(j, i) = correction(j, i) + cross(i, j)
            end do
         end do
         do i = 0, r - 1
            do j = 0, i
               out%covariance_band(i + 1, i - j + 1) = &
                  out%covariance_band(i + 1, i - j + 1) + correction(i + 1, j + 1)
            end do
         end do
      end if

      allocate(ma_covariance(q + 1))
      ma_covariance = 0.0_dp
      do lag = 0, q
         do j = lag, q
            ma_covariance(lag + 1) = ma_covariance(lag + 1) + &
               dot_product(matmul(disturbance_covariance, &
               numerator(:, j - lag + 1)), numerator(:, j + 1))
         end do
         do time = lag + 1, observations
            out%covariance_band(time, lag + 1) = &
               out%covariance_band(time, lag + 1) + ma_covariance(lag + 1)
         end do
      end do
      scale = maxval(out%covariance_band(:, 1))
      if (.not. ieee_is_finite(scale) .or. scale <= tiny(1.0_dp)) then
         out%info = 3
         return
      end if
      if (scale < 1.0_dp) scale = 1.0_dp
      out%covariance_band = out%covariance_band/scale
      factorization = tfarima_band_cholesky(out%covariance_band)
      if (factorization%info /= 0) then
         out%info = 4
         return
      end if
      whitened = tfarima_band_forward_solve(factorization%factor, conditional)
      if (size(whitened) /= observations) then
         out%info = 5
         return
      end if
      sum_squares = dot_product(whitened, whitened)
      if (sum_squares <= tiny(1.0_dp)) then
         out%info = 6
         return
      end if
      out%innovation_variance = sum_squares/(real(observations, dp)*scale)
      out%log_determinant = factorization%log_determinant + &
         0.5_dp*real(observations, dp)*log(scale)
      value = sum_squares/real(observations, dp)
      out%log_likelihood = -0.5_dp*real(observations, dp)*(1.0_dp + &
         log(2.0_dp*acos(-1.0_dp)) + log(value)) - &
         factorization%log_determinant
      out%residuals = whitened/sqrt(scale)
      out%standardized_residuals = out%residuals* &
         exp(factorization%log_determinant/real(observations, dp))
      out%cholesky_band = factorization%factor*sqrt(scale)
      out%covariance_band = out%covariance_band*scale
      out%bandwidth = r - 1
   end function tfarima_reduced_likelihood

   pure function tfarima_ucarima_reduced_likelihood(series, model) result(out)
      ! Evaluate a built UCARIMA model by exact reduced-form factorization.
      real(dp), intent(in) :: series(:)
      type(tfarima_ucarima_model_t), intent(in) :: model
      type(tfarima_reduced_likelihood_t) :: out
      real(dp), allocatable :: ar(:), difference(:), adjusted(:)
      real(dp), allocatable :: numerator(:, :), covariance(:, :)
      integer :: component

      if (model%info /= 0 .or. .not. allocated(model%components) .or. &
         .not. allocated(model%component_numerator)) then
         out%info = 1
         return
      end if
      ar = [1.0_dp]
      difference = [1.0_dp]
      do component = 1, size(model%components)
         ar = tfarima_polynomial_multiply(ar, &
            model%components(component)%ar_polynomial)
         difference = tfarima_polynomial_multiply(difference, &
            model%components(component)%difference_polynomial)
      end do
      adjusted = tfarima_difference(series, difference)
      numerator = transpose(model%component_numerator)
      allocate(covariance(size(model%components), size(model%components)))
      covariance = 0.0_dp
      do component = 1, size(model%components)
         covariance(component, component) = &
            model%components(component)%innovation_variance
      end do
      out = tfarima_reduced_likelihood(adjusted, ar, numerator, covariance)
   end function tfarima_ucarima_reduced_likelihood

   pure function tfarima_structural_reduced_likelihood(series, form, &
      difference_polynomial) result(out)
      ! Evaluate a structural SSM through its exact reduced representation.
      real(dp), intent(in) :: series(:), difference_polynomial(:)
      type(tfarima_ssm_form_t), intent(in) :: form
      type(tfarima_reduced_likelihood_t) :: out
      type(tfarima_ssm_reduction_t) :: reduction
      type(tfarima_polynomial_division_t) :: division
      real(dp), allocatable :: adjusted(:)

      if (size(difference_polynomial) < 1 .or. &
         abs(difference_polynomial(1)) <= tiny(1.0_dp)) then
         out%info = 1
         return
      end if
      reduction = tfarima_reduce_ssm(form)
      if (reduction%info /= 0) then
         out%info = 2
         return
      end if
      division = tfarima_polynomial_divide(reduction%characteristic, &
         difference_polynomial, 1.0e-8_dp)
      if (division%info /= 0 .or. &
         maxval(abs(division%remainder)) > 1.0e-7_dp) then
         out%info = 3
         return
      end if
      adjusted = tfarima_difference(series, difference_polynomial)
      out = tfarima_reduced_likelihood(adjusted, division%quotient, &
         reduction%disturbance_numerator, reduction%disturbance_covariance)
   end function tfarima_structural_reduced_likelihood

   pure function repeated_block_diagonal(block, repetitions) result(matrix)
      ! Repeat one square matrix along a block diagonal.
      real(dp), intent(in) :: block(:, :)
      integer, intent(in) :: repetitions
      real(dp), allocatable :: matrix(:, :)
      integer :: first, last, repetition, n

      n = size(block, 1)
      allocate(matrix(n*repetitions, n*repetitions))
      matrix = 0.0_dp
      do repetition = 1, repetitions
         first = (repetition - 1)*n + 1
         last = repetition*n
         matrix(first:last, first:last) = block
      end do
   end function repeated_block_diagonal

   pure function prepare_structural_series(series, log_transform, regressors, &
      regression_coefficients) result(adjusted)
      ! Apply an optional log transform and known deterministic regression signal.
      real(dp), intent(in) :: series(:)
      logical, intent(in), optional :: log_transform
      real(dp), intent(in), optional :: regressors(:, :)
      real(dp), intent(in), optional :: regression_coefficients(:)
      real(dp), allocatable :: adjusted(:)
      logical :: logarithm

      logarithm = .false.
      if (present(log_transform)) logarithm = log_transform
      if (size(series) < 1 .or. .not. all(ieee_is_finite(series)) .or. &
         (logarithm .and. any(series <= 0.0_dp)) .or. &
         (present(regressors) .neqv. present(regression_coefficients))) then
         allocate(adjusted(0))
         return
      end if
      if (present(regressors)) then
         if (size(regressors, 1) /= size(series) .or. &
            size(regressors, 2) /= size(regression_coefficients) .or. &
            .not. all(ieee_is_finite(regressors)) .or. &
            .not. all(ieee_is_finite(regression_coefficients))) then
            allocate(adjusted(0))
            return
         end if
      end if
      adjusted = series
      if (logarithm) adjusted = log(adjusted)
      if (present(regressors)) adjusted = adjusted - &
         matmul(regressors, regression_coefficients)
   end function prepare_structural_series

   pure function initialize_structural_series(series, form, observations) &
      result(out)
      ! Estimate initial state moments from filtered GLS observation equations.
      real(dp), intent(in) :: series(:)
      type(tfarima_ssm_form_t), intent(in) :: form
      integer, intent(in) :: observations
      type(tfarima_structural_initialization_t) :: out
      type(tfarima_structural_filter_t) :: filtered
      real(dp), allocatable :: zero(:), zero_covariance(:, :), design(:, :)
      real(dp), allocatable :: whitened(:, :), response(:), normal(:, :), inverse(:, :)
      real(dp), allocatable :: coefficients(:), residual(:), loading(:)
      real(dp) :: scale
      integer :: n, dimension, time, column, status

      n = observations
      dimension = size(form%observation_loading)
      if (n < dimension .or. n > size(series)) n = size(series)
      if (n < 1) then
         out%info = 1
         return
      end if
      allocate(zero(dimension), zero_covariance(dimension, dimension))
      zero = 0.0_dp
      zero_covariance = 0.0_dp
      filtered = run_structural_filter(series(:n), form, zero, zero_covariance)
      if (filtered%info /= 0) then
         out%info = 2
         return
      end if
      response = filtered%standardized_innovation
      allocate(design(n, dimension), whitened(n, dimension))
      loading = form%observation_loading
      do time = 1, n
         design(time, :) = loading
         loading = matmul(loading, form%transition)
      end do
      do column = 1, dimension
         filtered = run_structural_filter(design(:, column), form, zero, &
            zero_covariance)
         if (filtered%info /= 0) then
            out%info = 2
            return
         end if
         whitened(:, column) = filtered%standardized_innovation
      end do
      normal = matmul(transpose(whitened), whitened)
      call invert_matrix(normal, inverse, status)
      if (status /= 0) then
         normal = normal + 1.0e-8_dp*identity_array(dimension)
         call invert_matrix(normal, inverse, status)
      end if
      if (status /= 0) then
         out%info = 3
         return
      end if
      coefficients = matmul(inverse, matmul(transpose(whitened), response))
      residual = response - matmul(whitened, coefficients)
      scale = sum(residual**2)/real(max(1, n - dimension), dp)
      out%state = coefficients
      if (n <= dimension) then
         out%covariance = 1.0e4_dp*identity_array(dimension)
      else
         out%covariance = scale*inverse
      end if
      out%observations = n
   end function initialize_structural_series

   pure function run_structural_filter(series, form, initial_state, &
      initial_covariance) result(out)
      ! Run the lagged joint-disturbance Kalman innovation recursion.
      real(dp), intent(in) :: series(:), initial_state(:)
      type(tfarima_ssm_form_t), intent(in) :: form
      real(dp), intent(in) :: initial_covariance(:, :)
      type(tfarima_structural_filter_t) :: out
      real(dp), allocatable :: state(:), covariance(:, :), next_state(:)
      real(dp), allocatable :: next_covariance(:, :), cross(:), process(:, :)
      real(dp) :: variance
      integer :: dimension, n, time

      dimension = size(form%observation_loading)
      n = size(series)
      if (n < 1 .or. size(initial_state) /= dimension .or. &
         any(shape(initial_covariance) /= [dimension, dimension]) .or. &
         form%state_noise_contemporaneous) then
         out%info = 1
         return
      end if
      state = initial_state
      covariance = 0.5_dp*(initial_covariance + transpose(initial_covariance))
      cross = form%disturbance_covariance(2:, 1)
      process = form%disturbance_covariance(2:, 2:)
      allocate(out%observation_state(dimension, n), &
         out%observation_covariance(dimension, dimension, n), &
         out%filtered_state(dimension, n), &
         out%filtered_covariance(dimension, dimension, n), &
         out%gain(dimension, n), out%innovation(n), &
         out%innovation_variance(n), out%standardized_innovation(n))
      out%log_likelihood = 0.0_dp
      do time = 1, n
         out%observation_state(:, time) = state
         out%observation_covariance(:, :, time) = covariance
         out%innovation(time) = series(time) - &
            dot_product(form%observation_loading, state)
         variance = dot_product(form%observation_loading, &
            matmul(covariance, form%observation_loading)) + &
            form%disturbance_covariance(1, 1)
         if (.not. ieee_is_finite(variance) .or. variance <= 1.0e-12_dp) then
            out%info = 2
            return
         end if
         out%innovation_variance(time) = variance
         out%standardized_innovation(time) = out%innovation(time)/sqrt(variance)
         out%gain(:, time) = (matmul(form%transition, &
            matmul(covariance, form%observation_loading)) + cross)/variance
         next_state = matmul(form%transition, state)
         next_covariance = matmul(matmul(form%transition, covariance), &
            transpose(form%transition)) + process
         state = next_state + out%gain(:, time)*out%innovation(time)
         covariance = next_covariance - variance* &
            spread(out%gain(:, time), 2, dimension)* &
            spread(out%gain(:, time), 1, dimension)
         covariance = 0.5_dp*(covariance + transpose(covariance))
         out%filtered_state(:, time) = state
         out%filtered_covariance(:, :, time) = covariance
         out%log_likelihood = out%log_likelihood - 0.5_dp*( &
            log(2.0_dp*acos(-1.0_dp)) + log(variance) + &
            out%innovation(time)**2/variance)
      end do
   end function run_structural_filter

   pure recursive function root_decomposition_effects(series, mean_value, ar_polynomial, &
      difference_polynomial, ma_polynomial, basis, method) result(effects)
      ! Propagate root-effect coefficients from an exact ARMA innovation sequence.
      real(dp), intent(in) :: series(:), mean_value, ar_polynomial(:)
      real(dp), intent(in) :: difference_polynomial(:), ma_polynomial(:)
      real(dp), intent(in) :: basis(:, :)
      integer, intent(in) :: method
      real(dp), allocatable :: effects(:, :), forward(:, :), backward(:, :)
      real(dp), allocatable :: stationary(:), raw_innovations(:), innovations(:)
      real(dp), allocatable :: psi(:), initial_values(:), coefficient(:)
      real(dp), allocatable :: first_basis(:, :), next_basis(:, :), inverse(:, :)
      integer :: n, r, i, j, status

      if (method == tfarima_decomposition_backcast) then
         backward = root_decomposition_effects(series(size(series):1:-1), &
            mean_value, ar_polynomial, difference_polynomial, ma_polynomial, &
            basis, tfarima_decomposition_forecast)
         if (size(backward, 1) == size(series)) then
            effects = backward(size(series):1:-1, :)
         else
            allocate(effects(0, 0))
         end if
         return
      else if (method == tfarima_decomposition_mixed) then
         forward = root_decomposition_effects(series, mean_value, ar_polynomial, &
            difference_polynomial, ma_polynomial, basis, &
            tfarima_decomposition_forecast)
         backward = root_decomposition_effects(series, mean_value, ar_polynomial, &
            difference_polynomial, ma_polynomial, basis, &
            tfarima_decomposition_backcast)
         if (size(forward, 1) == size(series) .and. &
            size(backward, 1) == size(series)) then
            effects = 0.5_dp*(forward + backward)
         else
            allocate(effects(0, 0))
         end if
         return
      end if
      n = size(series)
      r = size(basis, 2)
      if (r < 1 .or. size(basis, 1) /= r + 1 .or. n < r) then
         allocate(effects(0, 0))
         return
      end if
      stationary = tfarima_difference(series, difference_polynomial) - mean_value
      raw_innovations = exact_arma_residuals(stationary, ar_polynomial, &
         ma_polynomial)
      if (size(raw_innovations) < 1) then
         allocate(effects(0, 0))
         return
      end if
      allocate(innovations(n))
      innovations = 0.0_dp
      if (size(raw_innovations) >= n) then
         innovations = raw_innovations(size(raw_innovations) - n + 1:)
      else
         innovations(n - size(raw_innovations) + 1:) = raw_innovations
      end if
      psi = tfarima_polynomial_ratio(ma_polynomial, &
         tfarima_polynomial_multiply(ar_polynomial, difference_polynomial), r)
      allocate(initial_values(r))
      do j = 1, r
         initial_values(j) = series(j) - innovations(j)
         do i = 1, j - 1
            initial_values(j) = initial_values(j) - &
               psi(i + 1)*innovations(j - i)
         end do
      end do
      first_basis = basis(:r, :)
      next_basis = basis(2:r + 1, :)
      call invert_matrix(first_basis, inverse, status)
      if (status /= 0) then
         allocate(effects(0, 0))
         return
      end if
      coefficient = matmul(inverse, initial_values)
      next_basis = matmul(inverse, next_basis)
      psi = matmul(inverse, psi(2:r + 1))
      allocate(effects(n, r + 1))
      do i = 1, n
         effects(i, :r) = coefficient
         effects(i, r + 1) = innovations(i)
         coefficient = matmul(next_basis, coefficient) + psi*innovations(i)
      end do
   end function root_decomposition_effects

   pure function exact_arma_residuals(series, ar_polynomial, ma_polynomial) &
      result(residuals)
      ! Compute exact ARMA residuals using estimated Gaussian presample values.
      real(dp), intent(in) :: series(:), ar_polynomial(:), ma_polynomial(:)
      real(dp), allocatable :: residuals(:), initial(:)
      real(dp) :: value
      integer :: n, p, q, time, lag, source

      n = size(series)
      p = size(ar_polynomial) - 1
      q = size(ma_polynomial) - 1
      if (n < 1 .or. p < 0 .or. q < 0) then
         allocate(residuals(0))
         return
      end if
      if (p + q == 0) then
         residuals = series
         return
      end if
      initial = arma_initial_conditions(series, ar_polynomial, ma_polynomial)
      if (size(initial) /= p + q) then
         allocate(residuals(0))
         return
      end if
      allocate(residuals(n + q))
      residuals = 0.0_dp
      if (q > 0) residuals(:q) = initial(p + 1:p + q)
      do time = 1, n
         value = 0.0_dp
         do lag = 0, p
            source = time - lag
            if (source >= 1) then
               value = value + ar_polynomial(lag + 1)*series(source)
            else
               value = value + ar_polynomial(lag + 1)*initial(p + source)
            end if
         end do
         do lag = 1, q
            value = value - ma_polynomial(lag + 1)* &
               residuals(time - lag + q)
         end do
         residuals(time + q) = value
      end do
   end function exact_arma_residuals

   pure function arma_initial_conditions(series, ar_polynomial, ma_polynomial) &
      result(initial)
      ! Estimate stationary ARMA presample series and innovation conditions.
      real(dp), intent(in) :: series(:), ar_polynomial(:), ma_polynomial(:)
      real(dp), allocatable :: initial(:)
      type(itsmr_arma_model_t) :: model
      real(dp), allocatable :: conditional(:), covariance(:), psi(:), pu(:, :)
      real(dp), allocatable :: lower(:, :), fpu(:, :), state(:), xx(:, :)
      real(dp), allocatable :: xy(:), information(:, :), inverse(:, :)
      real(dp) :: value
      integer :: n, p, q, r, s, i, j, time, status

      n = size(series)
      p = size(ar_polynomial) - 1
      q = size(ma_polynomial) - 1
      s = p + q
      r = max(p, q)
      if (s < 1 .or. n < r) then
         allocate(initial(0))
         return
      end if
      conditional = conditional_residuals(series, ar_polynomial, ma_polynomial)
      allocate(pu(s, s), fpu(r, s), state(r), xx(r, r), xy(r))
      pu = 0.0_dp
      fpu = 0.0_dp
      state = 0.0_dp
      xx = 0.0_dp
      xy = 0.0_dp
      model%ar = -ar_polynomial(2:)
      model%ma = ma_polynomial(2:)
      model%innovation_variance = 1.0_dp
      if (p > 0) then
         covariance = arma_acvf(model, p - 1)
         if (size(covariance) /= p) then
            allocate(initial(0))
            return
         end if
         do i = 1, p
            do j = 1, p
               pu(i, j) = covariance(abs(i - j))
            end do
         end do
      end if
      do i = 1, q
         pu(p + i, p + i) = 1.0_dp
      end do
      psi = tfarima_polynomial_ratio(ma_polynomial, ar_polynomial, q)
      do i = 0, min(p, q) - 1
         do j = i, min(p, q) - 1
            pu(p + q - j, p - j + i) = psi(i + 1)
            pu(p - j + i, p + q - j) = psi(i + 1)
         end do
      end do
      do i = 0, p - 1
         do j = i, p - 1
            fpu(i + 1, j + 1) = -ar_polynomial(p - j + i + 1)
         end do
      end do
      do i = 0, q - 1
         do j = i, q - 1
            fpu(i + 1, p + j + 1) = ma_polynomial(q - j + i + 1)
         end do
      end do
      call cholesky_lower(pu, lower, status)
      if (status /= 0) then
         allocate(initial(0))
         return
      end if
      fpu = matmul(fpu, lower)
      if (q > 0) then
         state(1) = 1.0_dp
         do time = 1, n
            xx = xx + spread(state, 2, r)*spread(state, 1, r)
            xy = xy + state*conditional(time)
            value = -dot_product(ma_polynomial(2:q + 1), state(:q))
            do j = r, 2, -1
               state(j) = state(j - 1)
            end do
            state(1) = value
         end do
      else
         do time = 1, p
            xy(time) = xy(time) + conditional(time)
            xx(time, time) = xx(time, time) + 1.0_dp
         end do
      end if
      information = matmul(matmul(transpose(fpu), xx), fpu)
      do i = 1, s
         information(i, i) = information(i, i) + 1.0_dp
      end do
      call invert_matrix(information, inverse, status)
      if (status /= 0) then
         allocate(initial(0))
         return
      end if
      xy = matmul(transpose(fpu), xy)
      initial = matmul(lower, matmul(inverse, xy))
   end function arma_initial_conditions

   pure function decomposition_polynomial_roots(polynomial) result(out)
      ! Find validated polynomial roots with an Aberth fallback for sparse cases.
      real(dp), intent(in) :: polynomial(:)
      type(arima2_roots_t) :: out
      complex(dp), allocatable :: previous(:), next(:)
      complex(dp) :: value, derivative, correction, interaction, denominator
      real(dp) :: radius, angle, change, residual, scale
      integer :: degree, iteration, i, j, power

      degree = size(polynomial) - 1
      if (degree < 1 .or. abs(polynomial(size(polynomial))) <= tiny(1.0_dp)) then
         out%info = 1
         return
      end if
      if (degree == 1) then
         allocate(out%roots(1))
         out%roots(1) = cmplx(-polynomial(1)/polynomial(2), 0.0_dp, dp)
         return
      end if

      allocate(out%roots(degree), previous(degree), next(degree))
      radius = abs(polynomial(1)/polynomial(size(polynomial)))** &
         (1.0_dp/real(degree, dp))
      radius = max(radius, 0.1_dp)
      do i = 1, degree
         angle = 2.0_dp*acos(-1.0_dp)*(real(i - 1, dp) + 0.173_dp)/ &
            real(degree, dp)
         out%roots(i) = radius*cmplx(cos(angle), sin(angle), dp)
      end do
      scale = max(1.0_dp, maxval(abs(polynomial)))
      do iteration = 1, 5000
         previous = out%roots
         change = 0.0_dp
         do i = 1, degree
            value = cmplx(polynomial(size(polynomial)), 0.0_dp, dp)
            derivative = cmplx(0.0_dp, 0.0_dp, dp)
            do power = degree - 1, 0, -1
               derivative = derivative*previous(i) + value
               value = value*previous(i) + polynomial(power + 1)
            end do
            if (abs(derivative) <= tiny(1.0_dp)) then
               next(i) = previous(i) + cmplx(1.0e-8_dp, 1.0e-8_dp, dp)
               cycle
            end if
            correction = value/derivative
            interaction = cmplx(0.0_dp, 0.0_dp, dp)
            do j = 1, degree
               if (j == i) cycle
               denominator = previous(i) - previous(j)
               if (abs(denominator) > tiny(1.0_dp)) &
                  interaction = interaction + 1.0_dp/denominator
            end do
            denominator = 1.0_dp - correction*interaction
            if (abs(denominator) > tiny(1.0_dp)) correction = correction/denominator
            next(i) = previous(i) - correction
            change = max(change, abs(correction))
         end do
         out%roots = next
         if (.not. all(ieee_is_finite(real(out%roots, dp))) .or. &
            .not. all(ieee_is_finite(aimag(out%roots)))) exit
         if (change <= 1.0e-12_dp*(1.0_dp + maxval(abs(out%roots)))) exit
      end do
      residual = polynomial_root_residual(polynomial, out%roots)
      if (.not. ieee_is_finite(residual) .or. residual > 1.0e-7_dp*scale) then
         out%info = 3
      else
         out%info = 0
      end if
   end function decomposition_polynomial_roots

   pure real(dp) function polynomial_root_residual(polynomial, roots) &
      result(residual)
      ! Return the largest absolute polynomial residual over candidate roots.
      real(dp), intent(in) :: polynomial(:)
      complex(dp), intent(in) :: roots(:)
      complex(dp) :: value
      integer :: i, power

      residual = 0.0_dp
      do i = 1, size(roots)
         value = cmplx(polynomial(size(polynomial)), 0.0_dp, dp)
         do power = size(polynomial) - 2, 0, -1
            value = value*roots(i) + polynomial(power + 1)
         end do
         residual = max(residual, abs(value))
      end do
   end function polynomial_root_residual

   pure complex(dp) function tfarima_polynomial_derivative(polynomial, z, order) &
      result(value)
      ! Evaluate a requested derivative of a lag polynomial at a complex point.
      real(dp), intent(in) :: polynomial(:)
      complex(dp), intent(in) :: z
      integer, intent(in) :: order
      integer :: degree

      value = cmplx(0.0_dp, 0.0_dp, dp)
      if (order < 0) return
      do degree = order, size(polynomial) - 1
         value = value + polynomial(degree + 1)*falling_factorial(degree, order)* &
            z**(degree - order)
      end do
   end function tfarima_polynomial_derivative

   pure function tfarima_transfer_spec(numerator, denominator, delay, &
      estimate_numerator, estimate_denominator) result(out)
      ! Construct a rational input specification with optional free masks.
      real(dp), intent(in) :: numerator(:), denominator(:)
      integer, intent(in), optional :: delay
      logical, intent(in), optional :: estimate_numerator(:)
      logical, intent(in), optional :: estimate_denominator(:)
      type(tfarima_transfer_spec_t) :: out
      integer :: selected_delay

      selected_delay = 0
      if (present(delay)) selected_delay = delay
      out%transfer = tfarima_transfer(numerator, denominator, selected_delay)
      if (out%transfer%info /= 0) then
         out%info = 1
         return
      end if
      allocate(out%estimate_numerator(size(numerator)), &
         out%estimate_denominator(size(denominator)))
      out%estimate_numerator = .true.
      out%estimate_denominator = .true.
      out%estimate_denominator(1) = .false.
      if (present(estimate_numerator)) then
         if (size(estimate_numerator) /= size(numerator)) then
            out%info = 2
            return
         end if
         out%estimate_numerator = estimate_numerator
      end if
      if (present(estimate_denominator)) then
         if (size(estimate_denominator) /= size(denominator) .or. &
            estimate_denominator(1)) then
            out%info = 2
            return
         end if
         out%estimate_denominator = estimate_denominator
      end if
   end function tfarima_transfer_spec

   pure function tfarima_prewhitened_ccf(input, output, input_model, max_lag, &
      output_model, confidence_level) result(out)
      ! Correlate input and output after conditional ARIMA prewhitening.
      real(dp), intent(in) :: input(:), output(:)
      type(tfarima_ucarima_component_t), intent(in) :: input_model
      integer, intent(in) :: max_lag
      type(tfarima_ucarima_component_t), intent(in), optional :: output_model
      real(dp), intent(in), optional :: confidence_level
      type(tfarima_prewhitened_ccf_t) :: out
      real(dp), allocatable :: input_residuals(:), output_residuals(:)
      real(dp) :: confidence, input_mean, output_mean, scale
      integer :: n, start_input, start_output, lag

      confidence = 0.95_dp
      if (present(confidence_level)) confidence = confidence_level
      if (size(input) < 2 .or. size(output) < 2 .or. max_lag < 0 .or. &
         confidence <= 0.0_dp .or. confidence >= 1.0_dp .or. &
         .not. valid_prewhitening_model(input_model)) then
         out%info = 1
         return
      end if
      if (present(output_model)) then
         if (.not. valid_prewhitening_model(output_model)) then
            out%info = 1
            return
         end if
         output_residuals = prewhiten_series(output, output_model)
      else
         output_residuals = prewhiten_series(output, input_model)
      end if
      input_residuals = prewhiten_series(input, input_model)
      n = min(size(input_residuals), size(output_residuals))
      if (n <= max_lag + 1) then
         out%info = 2
         return
      end if
      start_input = size(input_residuals) - n + 1
      start_output = size(output_residuals) - n + 1
      out%input_residuals = input_residuals(start_input:)
      out%output_residuals = output_residuals(start_output:)
      input_mean = sum(out%input_residuals)/real(n, dp)
      output_mean = sum(out%output_residuals)/real(n, dp)
      if (.not. all(ieee_is_finite(out%input_residuals)) .or. &
         .not. all(ieee_is_finite(out%output_residuals)) .or. &
         sum((out%input_residuals - input_mean)**2) <= tiny(1.0_dp) .or. &
         sum((out%output_residuals - output_mean)**2) <= tiny(1.0_dp)) then
         out%info = 3
         return
      end if
      out%correlation = ccf_values(out%input_residuals, out%output_residuals, &
         max_lag)
      allocate(out%lag(2*max_lag + 1), out%significant(2*max_lag + 1))
      out%lag = [(lag, lag=-max_lag,max_lag)]
      scale = sqrt(sum((out%output_residuals - output_mean)**2)/ &
         sum((out%input_residuals - input_mean)**2))
      out%impulse_response = out%correlation*scale
      out%critical_value = normal_quantile(0.5_dp + 0.5_dp*confidence)/sqrt(real(n, dp))
      out%significant = abs(out%correlation) > out%critical_value
   end function tfarima_prewhitened_ccf

   pure function tfarima_identify_transfer(output, input, input_model, &
      numerator_order, denominator_order, max_lag, output_model, delay, &
      confidence_level) result(out)
      ! Identify a rational transfer delay and starting coefficients.
      real(dp), intent(in) :: output(:), input(:)
      type(tfarima_ucarima_component_t), intent(in) :: input_model
      integer, intent(in) :: numerator_order, denominator_order, max_lag
      type(tfarima_ucarima_component_t), intent(in), optional :: output_model
      integer, intent(in), optional :: delay
      real(dp), intent(in), optional :: confidence_level
      type(tfarima_transfer_identification_t) :: out
      real(dp), allocatable :: impulse(:), numerator(:), denominator(:)
      real(dp), allocatable :: design(:, :), rhs(:), normal_matrix(:, :)
      real(dp), allocatable :: inverse(:, :)
      integer :: selected_delay, tail_lag, first_equation, equation_count
      integer :: equation, coefficient, status, shrink

      if (present(output_model)) then
         if (present(confidence_level)) then
            out%diagnostic = tfarima_prewhitened_ccf(input, output, input_model, &
               max_lag, output_model, confidence_level)
         else
            out%diagnostic = tfarima_prewhitened_ccf(input, output, input_model, &
               max_lag, output_model)
         end if
      else if (present(confidence_level)) then
         out%diagnostic = tfarima_prewhitened_ccf(input, output, input_model, &
            max_lag, confidence_level=confidence_level)
      else
         out%diagnostic = tfarima_prewhitened_ccf(input, output, input_model, max_lag)
      end if
      if (out%diagnostic%info /= 0 .or. numerator_order < 0 .or. &
         denominator_order < 0) then
         out%info = 1
         return
      end if
      if (present(delay)) then
         selected_delay = delay
      else
         selected_delay = first_significant_delay(out%diagnostic)
      end if
      if (selected_delay < 0 .or. selected_delay > max_lag) then
         out%info = 2
         return
      end if
      tail_lag = max_lag - selected_delay
      first_equation = max(numerator_order + 1, denominator_order)
      if (tail_lag < numerator_order .or. &
         (denominator_order > 0 .and. tail_lag - first_equation + 1 < &
         denominator_order)) then
         out%info = 3
         return
      end if
      allocate(impulse(0:tail_lag))
      impulse = out%diagnostic%impulse_response( &
         selected_delay + max_lag + 1:2*max_lag + 1)
      allocate(denominator(denominator_order + 1))
      denominator = 0.0_dp
      denominator(1) = 1.0_dp
      if (denominator_order > 0) then
         equation_count = tail_lag - first_equation + 1
         allocate(design(equation_count, denominator_order), rhs(equation_count))
         do equation = 1, equation_count
            tail_lag = first_equation + equation - 1
            rhs(equation) = -impulse(tail_lag)
            do coefficient = 1, denominator_order
               design(equation, coefficient) = impulse(tail_lag - coefficient)
            end do
         end do
         normal_matrix = matmul(transpose(design), design)
         call invert_matrix(normal_matrix, inverse, status)
         if (status /= 0) then
            out%info = 4
            return
         end if
         denominator(2:) = matmul(inverse, matmul(transpose(design), rhs))
         do shrink = 1, 100
            if (tfarima_polynomial_admissible(denominator, strict=.false.)) exit
            denominator(2:) = 0.95_dp*denominator(2:)
         end do
         if (.not. tfarima_polynomial_admissible(denominator, strict=.false.)) then
            out%info = 5
            return
         end if
      end if
      allocate(numerator(numerator_order + 1))
      do coefficient = 0, numerator_order
         numerator(coefficient + 1) = impulse(coefficient)
         do equation = 1, min(denominator_order, coefficient)
            numerator(coefficient + 1) = numerator(coefficient + 1) + &
               denominator(equation + 1)*impulse(coefficient - equation)
         end do
      end do
      out%selected_delay = selected_delay
      out%specification = tfarima_transfer_spec(numerator, denominator, &
         delay=selected_delay)
      if (out%specification%info /= 0) out%info = 6
   end function tfarima_identify_transfer

   pure function tfarima_exact_transfer_ccf(fit, input, input_models, &
      max_lag, confidence_level) result(out)
      ! Check fitted residuals against each separately prewhitened input.
      type(tfarima_exact_transfer_fit_t), intent(in) :: fit
      real(dp), intent(in) :: input(:, :)
      type(tfarima_ucarima_component_t), intent(in) :: input_models(:)
      integer, intent(in) :: max_lag
      real(dp), intent(in), optional :: confidence_level
      type(tfarima_prewhitened_ccf_t), allocatable :: out(:)
      type(tfarima_ucarima_component_t) :: identity_model
      integer :: dynamic_input

      allocate(out(size(input_models)))
      if (.not. allocated(fit%residuals) .or. size(input, 1) /= size(fit%residuals) .or. &
         size(input, 2) /= size(input_models)) then
         out%info = 1
         return
      end if
      identity_model = tfarima_ucarima_component([1.0_dp], [1.0_dp], &
         [1.0_dp], 1.0_dp)
      do dynamic_input = 1, size(input_models)
         if (present(confidence_level)) then
            out(dynamic_input) = tfarima_prewhitened_ccf( &
               input(:, dynamic_input), fit%residuals, input_models(dynamic_input), &
               max_lag, identity_model, confidence_level)
         else
            out(dynamic_input) = tfarima_prewhitened_ccf( &
               input(:, dynamic_input), fit%residuals, input_models(dynamic_input), &
               max_lag, identity_model)
         end if
      end do
   end function tfarima_exact_transfer_ccf

   pure function tfarima_diagnose_transfer(fit, max_lag, input, input_models, &
      confidence_level, fitted_parameter_count) result(out)
      ! Summarize fitted transfer residual dependence, shape, and spectrum.
      type(tfarima_exact_transfer_fit_t), intent(in) :: fit
      integer, intent(in) :: max_lag
      real(dp), intent(in), optional :: input(:, :)
      type(tfarima_ucarima_component_t), intent(in), optional :: input_models(:)
      real(dp), intent(in), optional :: confidence_level
      integer, intent(in), optional :: fitted_parameter_count
      type(tfarima_transfer_diagnostics_t) :: out
      real(dp), allocatable :: residuals(:), centered(:), periodogram(:)
      real(dp) :: confidence, alpha, scale, cosine_sum, sine_sum
      real(dp) :: cumulative, expected, pi
      integer :: n, parameter_count, frequencies, frequency, time

      confidence = 0.95_dp
      if (present(confidence_level)) confidence = confidence_level
      parameter_count = 0
      if (allocated(fit%parameters)) parameter_count = size(fit%parameters)
      if (present(fitted_parameter_count)) &
         parameter_count = fitted_parameter_count
      if (.not. allocated(fit%standardized_residuals) .or. max_lag < 1 .or. &
         confidence <= 0.0_dp .or. confidence >= 1.0_dp .or. &
         parameter_count < 0 .or. &
         (present(input) .neqv. present(input_models))) then
         out%info = 1
         return
      end if
      residuals = fit%standardized_residuals
      n = size(residuals)
      if (n <= max_lag .or. n < 4 .or. &
         .not. all(ieee_is_finite(residuals))) then
         out%info = 1
         return
      end if
      out%residual_mean = sum(residuals)/real(n, dp)
      centered = residuals - out%residual_mean
      scale = sum(centered**2)/real(n, dp)
      if (scale <= tiny(1.0_dp)) then
         out%info = 2
         return
      end if
      out%residual_variance = scale
      out%skewness = sum(centered**3)/real(n, dp)/scale**1.5_dp
      out%excess_kurtosis = sum(centered**4)/real(n, dp)/scale**2 - 3.0_dp
      out%normality_statistic = real(n, dp)*(out%skewness**2 + &
         0.25_dp*out%excess_kurtosis**2)/6.0_dp
      out%normality_p_value = exp(-0.5_dp*out%normality_statistic)
      out%residual_acf = acf_values(residuals, max_lag)
      out%residual_pacf = pacf_values(residuals, max_lag)
      out%ljung_box = weighted_box_test(residuals, max_lag, &
         box_test_ljung_box, parameter_count, .false., residual_raw)
      out%weighted_ljung_box = weighted_box_test(residuals, max_lag, &
         box_test_ljung_box, parameter_count, .true., residual_raw)
      if (out%ljung_box%info /= 0 .or. out%weighted_ljung_box%info /= 0) then
         out%info = 3
         return
      end if
      frequencies = n/2
      allocate(periodogram(frequencies), out%cumulative_periodogram(frequencies))
      pi = acos(-1.0_dp)
      do frequency = 1, frequencies
         cosine_sum = 0.0_dp
         sine_sum = 0.0_dp
         do time = 1, n
            cosine_sum = cosine_sum + centered(time)* &
               cos(2.0_dp*pi*real(frequency*time, dp)/real(n, dp))
            sine_sum = sine_sum - centered(time)* &
               sin(2.0_dp*pi*real(frequency*time, dp)/real(n, dp))
         end do
         periodogram(frequency) = cosine_sum**2 + sine_sum**2
      end do
      if (sum(periodogram) <= tiny(1.0_dp)) then
         out%info = 4
         return
      end if
      periodogram = periodogram/sum(periodogram)
      cumulative = 0.0_dp
      do frequency = 1, frequencies
         cumulative = cumulative + periodogram(frequency)
         out%cumulative_periodogram(frequency) = cumulative
         expected = real(frequency, dp)/real(frequencies, dp)
         out%cumulative_periodogram_statistic = max( &
            out%cumulative_periodogram_statistic, abs(cumulative - expected))
      end do
      alpha = 1.0_dp - confidence
      out%cumulative_periodogram_critical = &
         sqrt(-0.5_dp*log(0.5_dp*alpha)/real(frequencies, dp))
      out%cumulative_periodogram_p_value = min(1.0_dp, &
         2.0_dp*exp(-2.0_dp*real(frequencies, dp)* &
         out%cumulative_periodogram_statistic**2))
      if (present(input)) then
         out%input_ccf = tfarima_exact_transfer_ccf(fit, input, input_models, &
            max_lag, confidence)
         if (any(out%input_ccf%info /= 0)) out%info = 5
      else
         allocate(out%input_ccf(0))
      end if
   end function tfarima_diagnose_transfer

   pure function tfarima_simulate_transfer(input, specifications, regressors, &
      regression_coefficients) result(out)
      ! Evaluate multiple delayed rational inputs and deterministic regressors.
      real(dp), intent(in) :: input(:, :)
      type(tfarima_transfer_spec_t), intent(in) :: specifications(:)
      real(dp), intent(in), optional :: regressors(:, :), regression_coefficients(:)
      type(tfarima_transfer_signal_t) :: out
      integer :: n, input_count, dynamic_input

      n = size(input, 1)
      input_count = size(specifications)
      if (size(input, 2) /= input_count .or. n < 1 .or. &
         .not. all(ieee_is_finite(input)) .or. &
         (present(regressors) .neqv. present(regression_coefficients))) then
         out%info = 1
         return
      end if
      if (present(regressors)) then
         if (size(regressors, 1) /= n .or. &
            size(regressors, 2) /= size(regression_coefficients) .or. &
            .not. all(ieee_is_finite(regressors)) .or. &
            .not. all(ieee_is_finite(regression_coefficients))) then
            out%info = 1
            return
         end if
      end if
      allocate(out%input_signal(n, input_count), out%regression_signal(n), &
         out%total(n))
      out%input_signal = 0.0_dp
      do dynamic_input = 1, input_count
         if (specifications(dynamic_input)%info /= 0 .or. &
            specifications(dynamic_input)%transfer%info /= 0) then
            out%info = 2
            return
         end if
         out%input_signal(:, dynamic_input) = tfarima_filter( &
            input(:, dynamic_input), specifications(dynamic_input)%transfer)
      end do
      out%regression_signal = 0.0_dp
      if (present(regressors)) out%regression_signal = &
         matmul(regressors, regression_coefficients)
      out%total = out%regression_signal
      do dynamic_input = 1, input_count
         out%total = out%total + out%input_signal(:, dynamic_input)
      end do
   end function tfarima_simulate_transfer

   function tfarima_simulate_exact_model(fit, observations, simulations, input, &
      regressors, input_models, input_history, noise_history, &
      noise_innovations, input_innovations, noise_innovation_history, &
      input_innovation_history, burn_in, seed) result(out)
      ! Simulate exact transfer models with fixed or ARIMA-generated inputs.
      type(tfarima_exact_transfer_fit_t), intent(in) :: fit
      integer, intent(in) :: observations, simulations
      real(dp), intent(in), optional :: input(:, :), regressors(:, :)
      type(tfarima_ucarima_component_t), intent(in), optional :: input_models(:)
      real(dp), intent(in), optional :: input_history(:, :), noise_history(:)
      real(dp), intent(in), optional :: noise_innovations(:, :)
      real(dp), intent(in), optional :: input_innovations(:, :, :)
      real(dp), intent(in), optional :: noise_innovation_history(:)
      real(dp), intent(in), optional :: input_innovation_history(:, :)
      integer, intent(in), optional :: burn_in, seed
      type(tfarima_transfer_simulation_t) :: out
      type(tfarima_transfer_signal_t) :: signal
      real(dp), allocatable :: full_input(:, :, :), full_noise(:, :)
      real(dp), allocatable :: full_noise_innovation(:, :), full_regressors(:, :)
      real(dp), allocatable :: full_input_innovation(:, :, :), draws(:, :)
      real(dp), allocatable :: series_history(:), innovation_history(:)
      integer :: total, discarded, input_count, regression_count
      integer :: path, dynamic_input, first

      discarded = 0
      if (present(burn_in)) discarded = burn_in
      input_count = 0
      if (allocated(fit%inputs)) input_count = size(fit%inputs)
      regression_count = 0
      if (allocated(fit%regression_coefficients)) &
         regression_count = size(fit%regression_coefficients)
      total = observations + discarded
      if (observations < 1 .or. simulations < 1 .or. discarded < 0 .or. &
         .not. allocated(fit%inputs) .or. &
         .not. allocated(fit%regression_coefficients) .or. &
         .not. valid_prewhitening_model(fit%noise) .or. &
         fit%noise%innovation_variance < 0.0_dp .or. &
         (present(input) .and. present(input_models)) .or. &
         (input_count > 0 .and. .not. present(input) .and. &
         .not. present(input_models))) then
         out%info = 1
         return
      end if
      if (present(input)) then
         if (size(input, 1) /= total .or. size(input, 2) /= input_count .or. &
            .not. all(ieee_is_finite(input)) .or. present(input_innovations)) then
            out%info = 1
            return
         end if
      end if
      if (present(input_models)) then
         if (size(input_models) /= input_count) then
            out%info = 1
            return
         end if
         do dynamic_input = 1, input_count
            if (.not. valid_prewhitening_model(input_models(dynamic_input)) .or. &
               input_models(dynamic_input)%innovation_variance < 0.0_dp) then
               out%info = 1
               return
            end if
         end do
      end if
      if (present(regressors)) then
         if (size(regressors, 1) /= total .or. &
            size(regressors, 2) /= regression_count .or. &
            .not. all(ieee_is_finite(regressors))) then
            out%info = 1
            return
         end if
      else if (regression_count > 0) then
         out%info = 1
         return
      end if
      if (present(input_history)) then
         if (.not. present(input_models) .or. &
            size(input_history, 2) /= input_count .or. &
            .not. all(ieee_is_finite(input_history))) then
            out%info = 1
            return
         end if
      end if
      if (present(input_innovation_history)) then
         if (.not. present(input_models) .or. &
            size(input_innovation_history, 2) /= input_count .or. &
            .not. all(ieee_is_finite(input_innovation_history))) then
            out%info = 1
            return
         end if
      end if
      if (present(noise_history)) then
         if (.not. all(ieee_is_finite(noise_history))) then
            out%info = 1
            return
         end if
      end if
      if (present(noise_innovation_history)) then
         if (.not. all(ieee_is_finite(noise_innovation_history))) then
            out%info = 1
            return
         end if
      end if
      if (present(noise_innovations)) then
         if (any(shape(noise_innovations) /= [total, simulations]) .or. &
            .not. all(ieee_is_finite(noise_innovations))) then
            out%info = 1
            return
         end if
      end if
      if (present(input_innovations)) then
         if (.not. present(input_models) .or. &
            any(shape(input_innovations) /= [total, input_count, simulations]) .or. &
            .not. all(ieee_is_finite(input_innovations))) then
            out%info = 1
            return
         end if
      end if
      if (present(seed)) call set_random_seed(seed)
      allocate(full_noise_innovation(total, simulations))
      if (present(noise_innovations)) then
         full_noise_innovation = noise_innovations
      else
         call random_standard_normal_matrix(full_noise_innovation)
         full_noise_innovation = sqrt(fit%noise%innovation_variance)* &
            full_noise_innovation
      end if
      allocate(series_history(0), innovation_history(0))
      if (present(noise_history)) series_history = noise_history
      if (present(noise_innovation_history)) &
         innovation_history = noise_innovation_history
      allocate(full_noise(total, simulations))
      do path = 1, simulations
         full_noise(:, path) = simulate_arima_innovations(fit%noise, &
            full_noise_innovation(:, path), series_history, innovation_history)
      end do
      allocate(full_input(total, input_count, simulations), &
         full_input_innovation(total, input_count, simulations))
      full_input = 0.0_dp
      full_input_innovation = 0.0_dp
      if (present(input)) then
         do path = 1, simulations
            full_input(:, :, path) = input
         end do
      else if (input_count > 0) then
         do dynamic_input = 1, input_count
            allocate(draws(total, simulations))
            if (present(input_innovations)) then
               draws = input_innovations(:, dynamic_input, :)
            else
               call random_standard_normal_matrix(draws)
               draws = sqrt(input_models(dynamic_input)%innovation_variance)*draws
            end if
            full_input_innovation(:, dynamic_input, :) = draws
            deallocate(series_history, innovation_history)
            allocate(series_history(0), innovation_history(0))
            if (present(input_history)) &
               series_history = input_history(:, dynamic_input)
            if (present(input_innovation_history)) &
               innovation_history = input_innovation_history(:, dynamic_input)
            do path = 1, simulations
               full_input(:, dynamic_input, path) = simulate_arima_innovations( &
                  input_models(dynamic_input), draws(:, path), series_history, &
                  innovation_history)
            end do
            deallocate(draws)
         end do
      end if
      allocate(full_regressors(total, regression_count))
      if (regression_count > 0) full_regressors = regressors
      first = discarded + 1
      allocate(out%output(observations, simulations), &
         out%noise(observations, simulations), &
         out%innovation(observations, simulations), &
         out%input(observations, input_count, simulations), &
         out%input_innovation(observations, input_count, simulations), &
         out%input_signal(observations, input_count, simulations), &
         out%regression_signal(observations, simulations))
      do path = 1, simulations
         signal = tfarima_simulate_transfer(full_input(:, :, path), fit%inputs, &
            full_regressors, fit%regression_coefficients)
         if (signal%info /= 0) then
            out%info = 10 + signal%info
            return
         end if
         out%noise(:, path) = full_noise(first:, path)
         out%innovation(:, path) = full_noise_innovation(first:, path)
         out%input(:, :, path) = full_input(first:, :, path)
         out%input_innovation(:, :, path) = &
            full_input_innovation(first:, :, path)
         out%input_signal(:, :, path) = signal%input_signal(first:, :)
         out%regression_signal(:, path) = signal%regression_signal(first:)
         out%output(:, path) = signal%total(first:) + out%noise(:, path)
      end do
      out%observations = observations
      out%simulations = simulations
      out%burn_in = discarded
   end function tfarima_simulate_exact_model

   pure function tfarima_transfer(numerator, denominator, delay) result(out)
      ! Construct and validate a rational transfer function.
      real(dp), intent(in) :: numerator(:), denominator(:)
      integer, intent(in), optional :: delay
      type(tfarima_transfer_t) :: out

      out%delay = 0
      if (present(delay)) out%delay = abs(delay)
      if (size(numerator) < 1 .or. size(denominator) < 1 .or. &
         abs(denominator(1)) <= tiny(1.0_dp) .or. &
         .not. all(ieee_is_finite(numerator)) .or. &
         .not. all(ieee_is_finite(denominator))) then
         out%info = 1
         return
      end if
      out%numerator = numerator/denominator(1)
      out%denominator = denominator/denominator(1)
   end function tfarima_transfer

   pure function tfarima_filter(input, transfer) result(output)
      ! Filter an input through a delayed rational lag polynomial.
      real(dp), intent(in) :: input(:)
      type(tfarima_transfer_t), intent(in) :: transfer
      real(dp), allocatable :: output(:)
      real(dp) :: value
      integer :: time, lag, source

      allocate(output(size(input)))
      output = 0.0_dp
      if (transfer%info /= 0 .or. .not. allocated(transfer%numerator) .or. &
         .not. allocated(transfer%denominator)) return
      do time = 1, size(input)
         value = 0.0_dp
         do lag = 0, size(transfer%numerator) - 1
            source = time - transfer%delay - lag
            if (source >= 1) value = value + transfer%numerator(lag + 1)*input(source)
         end do
         do lag = 1, size(transfer%denominator) - 1
            if (time - lag >= 1) value = value - &
               transfer%denominator(lag + 1)*output(time - lag)
         end do
         output(time) = value
      end do
   end function tfarima_filter

   pure function tfarima_impulse_response(transfer, max_lag, cumulative) result(response)
      ! Compute a transfer function's impulse or cumulative step response.
      type(tfarima_transfer_t), intent(in) :: transfer
      integer, intent(in) :: max_lag
      logical, intent(in), optional :: cumulative
      real(dp), allocatable :: response(:), ratio(:)
      logical :: step
      integer :: lag

      step = .false.
      if (present(cumulative)) step = cumulative
      if (transfer%info /= 0 .or. max_lag < 0 .or. transfer%delay > max_lag) then
         allocate(response(0))
         return
      end if
      allocate(response(max_lag + 1))
      response = 0.0_dp
      ratio = tfarima_polynomial_ratio(transfer%numerator, transfer%denominator, &
         max_lag - transfer%delay)
      response(transfer%delay + 1:) = ratio
      if (step) then
         do lag = 2, size(response)
            response(lag) = response(lag) + response(lag - 1)
         end do
      end if
   end function tfarima_impulse_response

   pure function tfarima_difference(series, difference_polynomial, log_transform) &
      result(differenced)
      ! Apply an optional log transform and lag-difference polynomial.
      real(dp), intent(in) :: series(:), difference_polynomial(:)
      logical, intent(in), optional :: log_transform
      real(dp), allocatable :: differenced(:), transformed(:)
      logical :: logarithm
      integer :: degree, time, lag

      logarithm = .false.
      if (present(log_transform)) logarithm = log_transform
      degree = size(difference_polynomial) - 1
      if (degree < 0 .or. size(series) <= degree .or. &
         (logarithm .and. any(series <= 0.0_dp))) then
         allocate(differenced(0))
         return
      end if
      transformed = series
      if (logarithm) transformed = log(series)
      allocate(differenced(size(series) - degree))
      differenced = 0.0_dp
      do time = degree + 1, size(series)
         do lag = 0, degree
            differenced(time - degree) = differenced(time - degree) + &
               difference_polynomial(lag + 1)*transformed(time - lag)
         end do
      end do
   end function tfarima_difference

   pure function tfarima_arima_forecast(series, ar_polynomial, &
      difference_polynomial, ma_polynomial, innovation_variance, horizon, &
      mean_value, levels, log_transform) result(out)
      ! Forecast a tfarima-form ARIMA model by its ARMA recursions.
      real(dp), intent(in) :: series(:), ar_polynomial(:), difference_polynomial(:)
      real(dp), intent(in) :: ma_polynomial(:), innovation_variance
      integer, intent(in) :: horizon
      real(dp), intent(in), optional :: mean_value, levels(:)
      logical, intent(in), optional :: log_transform
      type(tfarima_forecast_t) :: out
      real(dp), allocatable :: transformed(:), stationary(:), innovations(:)
      real(dp), allocatable :: full_stationary(:), full_innovations(:), psi(:)
      real(dp) :: mu, value, quantile
      integer :: n, degree, p, q, h, lag, level
      logical :: logarithm

      mu = 0.0_dp
      if (present(mean_value)) mu = mean_value
      logarithm = .false.
      if (present(log_transform)) logarithm = log_transform
      degree = size(difference_polynomial) - 1
      p = size(ar_polynomial) - 1
      q = size(ma_polynomial) - 1
      if (horizon < 1 .or. degree < 0 .or. p < 0 .or. q < 0 .or. &
         size(series) <= degree .or. innovation_variance < 0.0_dp .or. &
         abs(ar_polynomial(1)) <= tiny(1.0_dp) .or. &
         abs(ma_polynomial(1)) <= tiny(1.0_dp) .or. &
         (logarithm .and. any(series <= 0.0_dp))) then
         out%info = 1
         return
      end if
      transformed = series
      if (logarithm) transformed = log(series)
      stationary = tfarima_difference(series, difference_polynomial, logarithm) - mu
      innovations = conditional_residuals(stationary, ar_polynomial, ma_polynomial)
      n = size(series)
      allocate(full_stationary(n + horizon), full_innovations(n + horizon))
      full_stationary = 0.0_dp
      full_innovations = 0.0_dp
      full_stationary(degree + 1:n) = stationary
      full_innovations(degree + 1:n) = innovations
      allocate(out%mean(horizon), out%variance(horizon))
      do h = 1, horizon
         value = 0.0_dp
         do lag = 1, p
            if (n + h - lag >= 1) value = value - &
               ar_polynomial(lag + 1)*full_stationary(n + h - lag)
         end do
         do lag = 1, q
            if (n + h - lag >= 1) value = value + &
               ma_polynomial(lag + 1)*full_innovations(n + h - lag)
         end do
         full_stationary(n + h) = value
         value = full_stationary(n + h) + mu
         do lag = 1, degree
            value = value - difference_polynomial(lag + 1)*transformed(n + h - lag)
         end do
         transformed = [transformed, value]
         if (logarithm) then
            out%mean(h) = exp(value)
         else
            out%mean(h) = value
         end if
      end do
      psi = tfarima_polynomial_ratio(ma_polynomial, &
         tfarima_polynomial_multiply(ar_polynomial, difference_polynomial), horizon - 1)
      value = 0.0_dp
      do h = 1, horizon
         value = value + psi(h)**2
         out%variance(h) = innovation_variance*value
      end do
      if (present(levels)) then
         if (any(levels <= 0.0_dp) .or. any(levels >= 1.0_dp)) then
            out%info = 2
            return
         end if
         allocate(out%lower(horizon, size(levels)), out%upper(horizon, size(levels)))
         do level = 1, size(levels)
            quantile = normal_quantile(0.5_dp + 0.5_dp*levels(level))
            if (logarithm) then
               out%lower(:, level) = out%mean*exp(-quantile*sqrt(out%variance))
               out%upper(:, level) = out%mean*exp(quantile*sqrt(out%variance))
            else
               out%lower(:, level) = out%mean - quantile*sqrt(out%variance)
               out%upper(:, level) = out%mean + quantile*sqrt(out%variance)
            end if
         end do
      end if
   end function tfarima_arima_forecast

   pure function tfarima_arima_backcast(series, ar_polynomial, &
      difference_polynomial, ma_polynomial, innovation_variance, horizon, &
      mean_value, log_transform) result(backcast)
      ! Backcast an ARIMA model by forecasting the reversed series.
      real(dp), intent(in) :: series(:), ar_polynomial(:), difference_polynomial(:)
      real(dp), intent(in) :: ma_polynomial(:), innovation_variance
      integer, intent(in) :: horizon
      real(dp), intent(in), optional :: mean_value
      logical, intent(in), optional :: log_transform
      real(dp), allocatable :: backcast(:), reversed(:)
      type(tfarima_forecast_t) :: forecast
      integer :: i

      allocate(reversed(size(series)))
      do i = 1, size(series)
         reversed(i) = series(size(series) - i + 1)
      end do
      forecast = tfarima_arima_forecast(reversed, ar_polynomial, &
         difference_polynomial, ma_polynomial, innovation_variance, horizon, &
         mean_value, log_transform=log_transform)
      if (forecast%info /= 0) then
         allocate(backcast(0))
         return
      end if
      allocate(backcast(horizon))
      do i = 1, horizon
         backcast(i) = forecast%mean(horizon - i + 1)
      end do
   end function tfarima_arima_backcast

   pure function tfarima_intervention(length, position, intervention_type) result(variable)
      ! Create pulse, step, or ramp intervention values.
      integer, intent(in) :: length, position, intervention_type
      real(dp), allocatable :: variable(:)
      integer :: i

      if (length < 1 .or. position < 1 .or. position > length .or. &
         intervention_type < tfarima_intervention_pulse .or. &
         intervention_type > tfarima_intervention_ramp) then
         allocate(variable(0))
         return
      end if
      allocate(variable(length))
      variable = 0.0_dp
      variable(position) = 1.0_dp
      if (intervention_type >= tfarima_intervention_step) then
         do i = position + 1, length
            variable(i) = variable(i - 1)
         end do
      end if
      if (intervention_type == tfarima_intervention_ramp) then
         do i = position + 1, length
            variable(i) = variable(i) + variable(i - 1)
         end do
      end if
   end function tfarima_intervention

   pure function tfarima_seasonal_dummies(length, period, start_season, &
      reference, constant) result(dummy)
      ! Create tfarima reference-coded seasonal dummy regressors.
      integer, intent(in) :: length, period, start_season, reference
      logical, intent(in), optional :: constant
      real(dp), allocatable :: dummy(:, :)
      integer :: columns, offset, time, season, target, column
      logical :: include_constant

      include_constant = .false.
      if (present(constant)) include_constant = constant
      if (length < 1 .or. period < 2 .or. start_season < 1 .or. &
         start_season > period .or. reference < 1 .or. reference > period) then
         allocate(dummy(0, 0))
         return
      end if
      offset = merge(1, 0, include_constant)
      columns = period - 1 + offset
      allocate(dummy(length, columns))
      dummy = 0.0_dp
      if (include_constant) dummy(:, 1) = 1.0_dp
      do time = 1, length
         season = modulo(start_season + time - 2, period) + 1
         column = offset
         do target = 1, period
            if (target == reference) cycle
            column = column + 1
            if (season == target) dummy(time, column) = 1.0_dp
            if (season == reference) dummy(time, column) = -1.0_dp
         end do
      end do
   end function tfarima_seasonal_dummies

   pure function tfarima_harmonic_regressors(length, period, start_season, &
      constant) result(regressors)
      ! Create tfarima cosine-sine seasonal regressors.
      integer, intent(in) :: length, period, start_season
      logical, intent(in), optional :: constant
      real(dp), allocatable :: regressors(:, :)
      real(dp) :: angle
      integer :: harmonics, columns, offset, time, harmonic, column
      logical :: include_constant

      include_constant = .false.
      if (present(constant)) include_constant = constant
      if (length < 1 .or. period < 2 .or. start_season < 1 .or. &
         start_season > period) then
         allocate(regressors(0, 0))
         return
      end if
      harmonics = (period - 1)/2
      offset = merge(1, 0, include_constant)
      columns = 2*harmonics + merge(1, 0, mod(period, 2) == 0) + offset
      allocate(regressors(length, columns))
      regressors = 0.0_dp
      if (include_constant) regressors(:, 1) = 1.0_dp
      do time = 1, length
         angle = 2.0_dp*acos(-1.0_dp)*real(start_season + time - 1, dp)/ &
            real(period, dp)
         column = offset
         do harmonic = 1, harmonics
            column = column + 1
            regressors(time, column) = cos(real(harmonic, dp)*angle)
            column = column + 1
            regressors(time, column) = sin(real(harmonic, dp)*angle)
         end do
         if (mod(period, 2) == 0) then
            regressors(time, columns) = cos(0.5_dp*real(period, dp)*angle)
         end if
      end do
   end function tfarima_harmonic_regressors

   pure function tfarima_standardize(values) result(standardized)
      ! Center and scale a vector using R's sample standard deviation.
      real(dp), intent(in) :: values(:)
      real(dp), allocatable :: standardized(:)
      real(dp) :: mean_value, standard_deviation

      if (size(values) < 2) then
         allocate(standardized(0))
         return
      end if
      mean_value = sum(values)/real(size(values), dp)
      standard_deviation = sqrt(sum((values - mean_value)**2)/ &
         real(size(values) - 1, dp))
      allocate(standardized(size(values)))
      if (standard_deviation <= tiny(1.0_dp)) then
         standardized = 0.0_dp
      else
         standardized = (values - mean_value)/standard_deviation
      end if
   end function tfarima_standardize

   pure function tfarima_autocovariance(ar_polynomial, ma_polynomial, &
      innovation_variance, max_lag) result(covariance)
      ! Return ARMA autocovariances using tfarima polynomial signs.
      real(dp), intent(in) :: ar_polynomial(:), ma_polynomial(:)
      real(dp), intent(in) :: innovation_variance
      integer, intent(in) :: max_lag
      real(dp), allocatable :: covariance(:)
      type(itsmr_arma_model_t) :: model

      if (size(ar_polynomial) < 1 .or. size(ma_polynomial) < 1 .or. &
         abs(ar_polynomial(1)) <= tiny(1.0_dp) .or. &
         abs(ma_polynomial(1)) <= tiny(1.0_dp) .or. &
         innovation_variance < 0.0_dp .or. max_lag < 0) then
         allocate(covariance(0))
         return
      end if
      allocate(model%ar(max(0, size(ar_polynomial) - 1)))
      allocate(model%ma(max(0, size(ma_polynomial) - 1)))
      if (size(model%ar) > 0) model%ar = -ar_polynomial(2:)/ar_polynomial(1)
      if (size(model%ma) > 0) model%ma = ma_polynomial(2:)/ma_polynomial(1)
      model%innovation_variance = innovation_variance* &
         (ma_polynomial(1)/ar_polynomial(1))**2
      covariance = arma_acvf(model, max_lag)
   end function tfarima_autocovariance

   pure function tfarima_partial_autocorrelation(ar_polynomial, ma_polynomial, &
      max_lag) result(partial)
      ! Compute theoretical ARMA partial autocorrelations by Levinson recursion.
      real(dp), intent(in) :: ar_polynomial(:), ma_polynomial(:)
      integer, intent(in) :: max_lag
      real(dp), allocatable :: partial(:), correlation(:), previous(:), current(:)
      real(dp) :: denominator, numerator
      integer :: order, j

      if (max_lag < 1) then
         allocate(partial(0))
         return
      end if
      correlation = tfarima_autocovariance(ar_polynomial, ma_polynomial, &
         1.0_dp, max_lag)
      if (size(correlation) /= max_lag + 1 .or. &
         correlation(1) <= tiny(1.0_dp)) then
         allocate(partial(0))
         return
      end if
      correlation = correlation/correlation(1)
      allocate(partial(max_lag), previous(max_lag), current(max_lag))
      partial = 0.0_dp
      previous = 0.0_dp
      current = 0.0_dp
      do order = 1, max_lag
         numerator = correlation(order + 1)
         denominator = 1.0_dp
         do j = 1, order - 1
            numerator = numerator - previous(j)*correlation(order - j + 1)
            denominator = denominator - previous(j)*correlation(j + 1)
         end do
         if (abs(denominator) <= tiny(1.0_dp)) then
            partial = huge(1.0_dp)
            return
         end if
         partial(order) = numerator/denominator
         current = 0.0_dp
         current(order) = partial(order)
         do j = 1, order - 1
            current(j) = previous(j) - partial(order)*previous(order - j)
         end do
         previous = current
      end do
   end function tfarima_partial_autocorrelation

   pure function tfarima_autocovariance_to_ma(covariance, tolerance, &
      max_iterations) result(out)
      ! Recover finite MA coefficients with tfarima's Newton factorization.
      real(dp), intent(in) :: covariance(:)
      real(dp), intent(in), optional :: tolerance
      integer, intent(in), optional :: max_iterations
      type(tfarima_ma_factor_t) :: out
      real(dp), allocatable :: coefficient(:), next_coefficient(:), &
         residual(:), jacobian(:, :), inverse(:, :)
      real(dp) :: tol
      integer :: limit, iteration, i, j, q, status

      tol = 1.0e-8_dp
      if (present(tolerance)) tol = tolerance
      limit = 500
      if (present(max_iterations)) limit = max_iterations
      out%method = tfarima_factor_newton
      q = size(covariance)
      if (q < 1 .or. covariance(1) <= 0.0_dp .or. tol <= 0.0_dp .or. &
         limit < 1 .or. .not. all(ieee_is_finite(covariance))) then
         out%info = 1
         return
      end if
      allocate(coefficient(q), next_coefficient(q), residual(q), jacobian(q, q))
      coefficient = 0.0_dp
      coefficient(1) = sqrt(covariance(1))
      do iteration = 1, limit
         jacobian = 0.0_dp
         do i = 1, q
            do j = 1, q - i + 1
               jacobian(i, j) = coefficient(i + j - 1)
            end do
            do j = i, q
               jacobian(i, j) = jacobian(i, j) + coefficient(j - i + 1)
            end do
            residual(i) = -covariance(i) + &
               dot_product(coefficient(:q - i + 1), coefficient(i:q))
         end do
         if (all(abs(residual) < tol)) then
            out%converged = .true.
            exit
         end if
         call invert_matrix(jacobian, inverse, status)
         if (status /= 0) then
            out%info = 2
            exit
         end if
         next_coefficient = coefficient - matmul(inverse, residual)
         coefficient = next_coefficient
      end do
      out%iterations = min(iteration, limit)
      out%coefficients = coefficient
      out%residual_norm = factor_residual_norm(covariance, coefficient)
      if (.not. out%converged .and. out%info == 0) out%info = 3
   end function tfarima_autocovariance_to_ma

   pure function tfarima_cramer_wold_factor(covariance, method, initial, &
      tolerance, max_iterations) result(out)
      ! Factor finite autocovariances with selectable robust Cramer-Wold methods.
      real(dp), intent(in) :: covariance(:)
      integer, intent(in), optional :: method, max_iterations
      real(dp), intent(in), optional :: initial(:), tolerance
      type(tfarima_ma_factor_t) :: out
      type(tfarima_ma_factor_t) :: candidate
      real(dp), allocatable :: g(:), start(:)
      real(dp) :: tol
      integer :: selected_method, limit

      selected_method = tfarima_factor_best
      if (present(method)) selected_method = method
      tol = 1.0e-8_dp
      if (present(tolerance)) tol = tolerance
      limit = 500
      if (present(max_iterations)) limit = max_iterations
      if (size(covariance) < 1 .or. covariance(1) <= 0.0_dp .or. &
         .not. all(ieee_is_finite(covariance)) .or. tol <= 0.0_dp .or. &
         limit < 1 .or. selected_method < tfarima_factor_newton .or. &
         selected_method > tfarima_factor_best) then
         out%info = 1
         return
      end if
      g = trim_polynomial(covariance, tol)
      allocate(start(0))
      if (present(initial)) then
         if (size(initial) /= size(g) .or. .not. all(ieee_is_finite(initial))) then
            out%info = 1
            return
         end if
         start = initial
      end if
      select case (selected_method)
      case (tfarima_factor_newton)
         out = tfarima_autocovariance_to_ma(g, tol, limit)
      case (tfarima_factor_roots)
         out = cramer_wold_roots(g, tol)
      case (tfarima_factor_bauer)
         out = cramer_wold_bauer(g, tol, limit)
      case (tfarima_factor_laurie)
         out = cramer_wold_laurie(g, start, tol, limit)
      case (tfarima_factor_wilson)
         out = cramer_wold_wilson(g, start, tol, limit)
      case (tfarima_factor_best)
         out = tfarima_autocovariance_to_ma(g, tol, limit)
         candidate = cramer_wold_roots(g, tol)
         call retain_better_factor(out, candidate)
         candidate = cramer_wold_bauer(g, tol, limit)
         call retain_better_factor(out, candidate)
         candidate = cramer_wold_laurie(g, start, tol, limit)
         call retain_better_factor(out, candidate)
         candidate = cramer_wold_wilson(g, start, tol, limit)
         call retain_better_factor(out, candidate)
      end select
   end function tfarima_cramer_wold_factor

   pure function cramer_wold_roots(covariance, tolerance) result(out)
      ! Select the minimum-phase half of the palindromic covariance roots.
      real(dp), intent(in) :: covariance(:), tolerance
      type(tfarima_ma_factor_t) :: out
      type(arima2_roots_t) :: roots
      real(dp), allocatable :: palindromic(:)
      complex(dp), allocatable :: ordered(:), polynomial(:), next(:)
      complex(dp) :: temporary
      real(dp) :: scale
      integer :: n, order, i, j, selected

      out%method = tfarima_factor_roots
      n = size(covariance)
      if (n == 1) then
         out%coefficients = [sqrt(covariance(1))]
         out%residual_norm = 0.0_dp
         out%converged = .true.
         return
      end if
      palindromic = [covariance(n:1:-1), covariance(2:n)]
      if (abs(palindromic(1)) <= tiny(1.0_dp)) then
         out%info = 1
         return
      end if
      roots = arma_polynomial_roots(-palindromic(2:)/palindromic(1))
      if (roots%info /= 0 .or. size(roots%roots) /= 2*(n - 1)) then
         out%info = 2
         return
      end if
      ordered = roots%roots
      do i = 1, size(ordered) - 1
         selected = i
         do j = i + 1, size(ordered)
            if (abs(ordered(j)) < abs(ordered(selected))) selected = j
         end do
         if (selected /= i) then
            temporary = ordered(i)
            ordered(i) = ordered(selected)
            ordered(selected) = temporary
         end if
      end do
      allocate(polynomial(1))
      polynomial(1) = cmplx(1.0_dp, 0.0_dp, dp)
      do order = 1, n - 1
         allocate(next(order + 1))
         next = cmplx(0.0_dp, 0.0_dp, dp)
         next(:order) = next(:order) + polynomial
         next(2:) = next(2:) - ordered(order)*polynomial
         call move_alloc(next, polynomial)
      end do
      if (abs(real(polynomial(n), dp)) <= tiny(1.0_dp) .or. &
         maxval(abs(aimag(polynomial))) > 100.0_dp*tolerance) then
         out%info = 3
         return
      end if
      scale = covariance(n)/real(polynomial(n), dp)
      out%coefficients = real(polynomial, dp)*sign(sqrt(abs(scale)), scale)
      out%residual_norm = factor_residual_norm(covariance, out%coefficients)
      out%converged = ieee_is_finite(out%residual_norm) .and. &
         out%residual_norm <= tolerance*max(1.0_dp, sqrt(sum(covariance**2)))
      if (.not. out%converged) out%info = 4
   end function cramer_wold_roots

   pure function cramer_wold_bauer(covariance, tolerance, max_iterations) &
      result(out)
      ! Initialize a Cramer-Wold factor with Laurie's Bauer recursion.
      real(dp), intent(in) :: covariance(:), tolerance
      integer, intent(in) :: max_iterations
      type(tfarima_ma_factor_t) :: out
      real(dp), allocatable :: g(:), theta(:), previous(:)
      real(dp) :: ratio, g0
      integer :: n, iteration, i

      out%method = tfarima_factor_bauer
      n = size(covariance)
      g0 = covariance(1)
      g = covariance/g0
      allocate(theta(n), previous(n))
      theta = [g(2:), 0.0_dp]
      do iteration = 1, max_iterations
         previous = theta
         if (g(1) <= abs(g(n))) exit
         ratio = -theta(1)/g(1)
         if (abs(ratio) > 1.0_dp) exit
         do i = 1, n - 1
            g(i) = g(i) + ratio*theta(i)
            theta(i) = theta(i + 1) + ratio*g(i + 1)
         end do
         if (maxval(abs(theta - previous)) < tolerance) then
            out%converged = .true.
            exit
         end if
      end do
      if (g(1) > 0.0_dp) theta = sqrt(g0)*g/sqrt(g(1))
      out%coefficients = theta
      out%iterations = min(iteration, max_iterations)
      out%residual_norm = factor_residual_norm(covariance, theta)
      if (.not. out%converged .and. out%residual_norm <= tolerance) &
         out%converged = .true.
      if (.not. out%converged) out%info = 3
   end function cramer_wold_bauer

   pure function cramer_wold_laurie(covariance, initial, tolerance, &
      max_iterations) result(out)
      ! Refine a Cramer-Wold factor with Laurie algorithm AS 175.
      real(dp), intent(in) :: covariance(:), initial(:), tolerance
      integer, intent(in) :: max_iterations
      type(tfarima_ma_factor_t) :: out
      type(tfarima_ma_factor_t) :: bauer
      real(dp), allocatable :: theta(:), next(:)
      integer :: n, iteration

      out%method = tfarima_factor_laurie
      n = size(covariance)
      if (size(initial) == n) then
         theta = initial
      else
         bauer = cramer_wold_bauer(covariance, tolerance, max_iterations)
         if (.not. allocated(bauer%coefficients)) then
            out%info = 2
            return
         end if
         theta = bauer%coefficients
         theta(1) = abs(theta(1))
      end if
      do iteration = 1, max_iterations
         next = laurie_iteration(covariance, theta)
         if (next(1) <= 0.0_dp) exit
         next = next + 0.5_dp*theta
         if (abs(next(n)/next(1)) > 1.0_dp .and. &
            abs(next(n)) > tiny(1.0_dp)) next(n) = 1.0_dp/next(n)
         if (maxval(abs(next - theta))/max(abs(theta(1)), tiny(1.0_dp)) < &
            tolerance) then
            theta = next
            out%converged = .true.
            exit
         end if
         theta = next
      end do
      out%coefficients = theta
      out%iterations = min(iteration, max_iterations)
      out%residual_norm = factor_residual_norm(covariance, theta)
      if (.not. out%converged .and. out%residual_norm <= tolerance) &
         out%converged = .true.
      if (.not. out%converged) out%info = 3
   end function cramer_wold_laurie

   pure function laurie_iteration(covariance, theta) result(next)
      ! Perform one Laurie AS 175 factor recursion.
      real(dp), intent(in) :: covariance(:), theta(:)
      real(dp) :: next(size(theta)), g(size(theta))
      real(dp) :: ratio, saved
      integer :: order, m, i

      order = size(theta) - 1
      m = order + 1
      g = covariance
      g(1) = 0.5_dp*g(1)
      next = theta
      do i = 1, order
         if (next(1) <= 0.0_dp) return
         ratio = g(m)/next(1)
         g(1:m) = g(1:m) - ratio*next(m:1:-1)
         g(m) = ratio
         ratio = next(m)/next(1)
         next(1:m) = next(1:m) - ratio*next(m:1:-1)
         next(m) = ratio
         m = m - 1
      end do
      if (next(1) <= 0.0_dp) return
      next(1) = g(1)/next(1)
      do i = 2, order + 1
         saved = next(i)
         next(i) = 0.0_dp
         next(1:i) = next(1:i) - saved*next(i:1:-1)
         next(i) = next(i) + g(i)
      end do
   end function laurie_iteration

   pure function cramer_wold_wilson(covariance, initial, tolerance, &
      max_iterations) result(out)
      ! Solve covariance moment equations by Wilson's matrix iteration.
      real(dp), intent(in) :: covariance(:), initial(:), tolerance
      integer, intent(in) :: max_iterations
      type(tfarima_ma_factor_t) :: out
      real(dp), allocatable :: g(:), theta(:), moments(:), residual(:)
      real(dp), allocatable :: matrix(:, :), inverse(:, :)
      real(dp) :: g0, length
      integer :: n, iteration, column, status

      out%method = tfarima_factor_wilson
      n = size(covariance)
      g0 = covariance(1)
      g = covariance/g0
      allocate(theta(n), moments(n), residual(n), matrix(n, n))
      theta = 0.0_dp
      theta(1) = 1.0_dp
      if (size(initial) == n) theta = initial/sqrt(g0)
      do iteration = 1, max_iterations
         call factor_moments(theta, moments)
         residual = moments - g
         if (maxval(abs(residual)) < tolerance) then
            out%converged = .true.
            exit
         end if
         matrix = 0.0_dp
         do column = 1, n
            matrix(:n - column + 1, column) = &
               theta(column:n)
            matrix(:column, column) = matrix(:column, column) + &
               theta(column:1:-1)
         end do
         call invert_matrix(matrix, inverse, status)
         if (status /= 0) then
            do column = 1, n
               matrix(column, column) = matrix(column, column) + 0.001_dp
            end do
            call invert_matrix(matrix, inverse, status)
            if (status /= 0) then
               out%info = 2
               exit
            end if
         end if
         theta = theta - matmul(inverse, residual)
         length = sqrt(sum(theta**2))
         if (length <= tiny(1.0_dp)) then
            out%info = 2
            exit
         end if
         theta = theta/length
      end do
      out%coefficients = theta*sqrt(g0)
      out%iterations = min(iteration, max_iterations)
      out%residual_norm = factor_residual_norm(covariance, out%coefficients)
      if (.not. out%converged .and. out%residual_norm <= tolerance) &
         out%converged = .true.
      if (.not. out%converged .and. out%info == 0) out%info = 3
   end function cramer_wold_wilson

   pure subroutine factor_moments(coefficients, moments)
      ! Reconstruct one-sided autocovariances from finite MA coefficients.
      real(dp), intent(in) :: coefficients(:)
      real(dp), intent(out) :: moments(:)
      integer :: lag, n

      n = size(coefficients)
      do lag = 0, n - 1
         moments(lag + 1) = dot_product(coefficients(:n - lag), &
            coefficients(lag + 1:))
      end do
   end subroutine factor_moments

   pure real(dp) function factor_residual_norm(covariance, coefficients) &
      result(value)
      ! Return the Euclidean autocovariance reconstruction error.
      real(dp), intent(in) :: covariance(:), coefficients(:)
      real(dp), allocatable :: moments(:)

      if (size(covariance) /= size(coefficients)) then
         value = huge(1.0_dp)
         return
      end if
      allocate(moments(size(coefficients)))
      call factor_moments(coefficients, moments)
      value = sqrt(sum((moments - covariance)**2))
   end function factor_residual_norm

   pure subroutine retain_better_factor(best, candidate)
      ! Retain a valid factor when it improves the reconstruction residual.
      type(tfarima_ma_factor_t), intent(inout) :: best
      type(tfarima_ma_factor_t), intent(in) :: candidate

      if (candidate%info /= 0 .or. .not. allocated(candidate%coefficients)) return
      if (best%info /= 0 .or. .not. allocated(best%coefficients) .or. &
         candidate%residual_norm < best%residual_norm) best = candidate
   end subroutine retain_better_factor

   pure function tfarima_palindromic_to_wold(polynomial) result(wold)
      ! Convert symmetric palindromic coordinates to a Wold polynomial.
      real(dp), intent(in) :: polynomial(:)
      real(dp), allocatable :: wold(:)

      wold = wold_from_palindromic(polynomial)
   end function tfarima_palindromic_to_wold

   pure function tfarima_wold_to_palindromic(wold) result(polynomial)
      ! Convert a Wold polynomial back to symmetric palindromic coordinates.
      real(dp), intent(in) :: wold(:)
      real(dp), allocatable :: polynomial(:)

      polynomial = palindromic_from_wold(wold)
   end function tfarima_wold_to_palindromic

   pure function tfarima_lag_polynomial(lags, offset, loading, parameters, &
      exponent) result(out)
      ! Construct a sparse lag polynomial with linear parameter restrictions.
      integer, intent(in) :: lags(:)
      real(dp), intent(in) :: offset(:), loading(:, :), parameters(:)
      integer, intent(in), optional :: exponent
      type(tfarima_lag_polynomial_t) :: out
      real(dp), allocatable :: coefficients(:)

      out%exponent = 1
      if (present(exponent)) out%exponent = exponent
      if (size(lags) < 1 .or. size(offset) /= size(lags) .or. &
         size(loading, 1) /= size(lags) .or. &
         size(loading, 2) /= size(parameters) .or. any(lags < 1) .or. &
         any(lags(2:) <= lags(:size(lags) - 1)) .or. out%exponent < 1 .or. &
         .not. all(ieee_is_finite(offset)) .or. &
         .not. all(ieee_is_finite(loading)) .or. &
         .not. all(ieee_is_finite(parameters))) then
         out%info = 1
         return
      end if
      out%lags = lags
      out%offset = offset
      out%loading = loading
      out%parameters = parameters
      coefficients = offset + matmul(loading, parameters)
      allocate(out%base(maxval(lags) + 1))
      out%base = 0.0_dp
      out%base(1) = 1.0_dp
      out%base(lags + 1) = -coefficients
      out%expanded = tfarima_polynomial_power(out%base, out%exponent)
   end function tfarima_lag_polynomial

   pure function tfarima_update_lag_polynomial(model, parameters) result(out)
      ! Evaluate a restricted lag polynomial at new parameter values.
      type(tfarima_lag_polynomial_t), intent(in) :: model
      real(dp), intent(in) :: parameters(:)
      type(tfarima_lag_polynomial_t) :: out

      if (model%info /= 0 .or. .not. allocated(model%lags) .or. &
         .not. allocated(model%offset) .or. .not. allocated(model%loading)) then
         out%info = 1
         return
      end if
      out = tfarima_lag_polynomial(model%lags, model%offset, model%loading, &
         parameters, model%exponent)
   end function tfarima_update_lag_polynomial

   pure logical function tfarima_polynomial_admissible(polynomial, strict, &
      tolerance) result(admissible)
      ! Test whether every lag-polynomial root is outside the unit circle.
      real(dp), intent(in) :: polynomial(:)
      logical, intent(in), optional :: strict
      real(dp), intent(in), optional :: tolerance
      type(arima2_roots_t) :: root_result
      real(dp) :: tol
      logical :: outside

      outside = .true.
      if (present(strict)) outside = strict
      tol = 1.0e-8_dp
      if (present(tolerance)) tol = tolerance
      if (size(polynomial) < 1 .or. abs(polynomial(1)) <= tiny(1.0_dp) .or. &
         .not. all(ieee_is_finite(polynomial))) then
         admissible = .false.
         return
      end if
      if (size(polynomial) == 1) then
         admissible = .true.
         return
      end if
      root_result = arma_polynomial_roots(-polynomial(2:)/polynomial(1))
      if (root_result%info /= 0) then
         admissible = .false.
      else if (outside) then
         admissible = all(abs(root_result%roots) > 1.0_dp + tol)
      else
         admissible = all(abs(root_result%roots) >= 1.0_dp - tol)
      end if
   end function tfarima_polynomial_admissible

   pure function tfarima_exact_transfer_fit(output, input, specifications, &
      noise_model, regressors, initial_regression, estimate_noise_ar, &
      estimate_noise_ma, estimate_noise_variance, max_iterations, tolerance) &
      result(out)
      ! Fit multiple rational inputs and ARIMA noise by exact KFAS likelihood.
      real(dp), intent(in) :: output(:), input(:, :)
      type(tfarima_transfer_spec_t), intent(in) :: specifications(:)
      type(tfarima_ucarima_component_t), intent(in) :: noise_model
      real(dp), intent(in), optional :: regressors(:, :), initial_regression(:)
      logical, intent(in), optional :: estimate_noise_ar(:), estimate_noise_ma(:)
      logical, intent(in), optional :: estimate_noise_variance
      integer, intent(in), optional :: max_iterations
      real(dp), intent(in), optional :: tolerance
      type(tfarima_exact_transfer_fit_t) :: out
      type(optimization_result_t) :: warm_start, optimization
      type(kfs_filter_t) :: filtered
      type(tfarima_transfer_spec_t), allocatable :: fitted_specs(:)
      type(tfarima_ucarima_component_t) :: fitted_noise
      logical, allocatable :: noise_ar_mask(:), noise_ma_mask(:)
      logical, allocatable :: variance_parameter(:)
      real(dp), allocatable :: initial(:), coordinates(:), beta(:), adjusted(:)
      real(dp), allocatable :: hessian(:, :), inverse(:, :), jacobian(:, :)
      real(dp) :: tol
      integer :: input_count, nreg, npar, position, dynamic_input, coefficient
      integer :: limit, status, observations, i, bfgs_iterations
      logical :: fit_variance

      tol = 1.0e-6_dp
      if (present(tolerance)) tol = tolerance
      limit = 250
      if (present(max_iterations)) limit = max_iterations
      input_count = size(specifications)
      if (size(output) < 4 .or. size(input, 1) /= size(output) .or. &
         size(input, 2) /= input_count .or. noise_model%info /= 0 .or. &
         .not. allocated(noise_model%ar_polynomial) .or. &
         .not. allocated(noise_model%difference_polynomial) .or. &
         .not. allocated(noise_model%ma_polynomial) .or. &
         noise_model%innovation_variance <= 0.0_dp .or. &
         .not. all(ieee_is_finite(output)) .or. &
         .not. all(ieee_is_finite(input)) .or. tol <= 0.0_dp .or. limit < 1) then
         out%info = 1
         return
      end if
      if (.not. tfarima_polynomial_admissible(noise_model%ar_polynomial) .or. &
         .not. tfarima_polynomial_admissible(noise_model%ma_polynomial)) then
         out%info = 1
         return
      end if
      do dynamic_input = 1, input_count
         if (specifications(dynamic_input)%info /= 0 .or. &
            specifications(dynamic_input)%transfer%info /= 0 .or. &
            .not. allocated(specifications(dynamic_input)%transfer%numerator) .or. &
            .not. allocated(specifications(dynamic_input)%transfer%denominator) .or. &
            .not. allocated(specifications(dynamic_input)%estimate_numerator) .or. &
            .not. allocated(specifications(dynamic_input)%estimate_denominator)) then
            out%info = 2
            return
         end if
         if (size(specifications(dynamic_input)%estimate_numerator) /= &
            size(specifications(dynamic_input)%transfer%numerator) .or. &
            size(specifications(dynamic_input)%estimate_denominator) /= &
            size(specifications(dynamic_input)%transfer%denominator) .or. &
            specifications(dynamic_input)%estimate_denominator(1) .or. &
            .not. tfarima_polynomial_admissible( &
            specifications(dynamic_input)%transfer%denominator)) then
            out%info = 2
            return
         end if
      end do
      allocate(noise_ar_mask(size(noise_model%ar_polynomial) - 1), &
         noise_ma_mask(size(noise_model%ma_polynomial) - 1))
      noise_ar_mask = .false.
      noise_ma_mask = .false.
      if (present(estimate_noise_ar)) then
         if (size(estimate_noise_ar) /= size(noise_ar_mask)) then
            out%info = 3
            return
         end if
         noise_ar_mask = estimate_noise_ar
      end if
      if (present(estimate_noise_ma)) then
         if (size(estimate_noise_ma) /= size(noise_ma_mask)) then
            out%info = 3
            return
         end if
         noise_ma_mask = estimate_noise_ma
      end if
      fit_variance = .true.
      if (present(estimate_noise_variance)) &
         fit_variance = estimate_noise_variance
      nreg = 0
      if (present(regressors)) then
         if (size(regressors, 1) /= size(output) .or. &
            .not. all(ieee_is_finite(regressors))) then
            out%info = 4
            return
         end if
         nreg = size(regressors, 2)
      end if
      if (present(initial_regression)) then
         if (.not. present(regressors) .or. size(initial_regression) /= nreg) then
            out%info = 4
            return
         end if
      end if
      npar = nreg + count(noise_ar_mask) + count(noise_ma_mask) + &
         merge(1, 0, fit_variance)
      do dynamic_input = 1, input_count
         npar = npar + count(specifications(dynamic_input)%estimate_numerator)
         npar = npar + count(specifications(dynamic_input)%estimate_denominator)
      end do
      allocate(initial(npar), variance_parameter(npar))
      variance_parameter = .false.
      position = 0
      if (nreg > 0) then
         initial(:nreg) = 0.0_dp
         if (present(initial_regression)) initial(:nreg) = initial_regression
         position = nreg
      end if
      do dynamic_input = 1, input_count
         do coefficient = 1, &
            size(specifications(dynamic_input)%transfer%numerator)
            if (.not. specifications(dynamic_input)% &
               estimate_numerator(coefficient)) cycle
            position = position + 1
            initial(position) = &
               specifications(dynamic_input)%transfer%numerator(coefficient)
         end do
         do coefficient = 1, &
            size(specifications(dynamic_input)%transfer%denominator)
            if (.not. specifications(dynamic_input)% &
               estimate_denominator(coefficient)) cycle
            position = position + 1
            initial(position) = &
               specifications(dynamic_input)%transfer%denominator(coefficient)
         end do
      end do
      do coefficient = 1, size(noise_ar_mask)
         if (.not. noise_ar_mask(coefficient)) cycle
         position = position + 1
         initial(position) = noise_model%ar_polynomial(coefficient + 1)
      end do
      do coefficient = 1, size(noise_ma_mask)
         if (.not. noise_ma_mask(coefficient)) cycle
         position = position + 1
         initial(position) = noise_model%ma_polynomial(coefficient + 1)
      end do
      if (fit_variance) then
         position = position + 1
         initial(position) = 0.5_dp*log(max(noise_model%innovation_variance, tol))
         variance_parameter(position) = .true.
      end if
      if (npar > 0) then
         warm_start = nelder_mead_minimize(objective, initial, &
            max_iterations=limit, tolerance=tol, initial_step=0.05_dp)
         optimization = bfgs_minimize_fd(objective, warm_start%parameters, &
            max_iterations=limit, gradient_tolerance=tol)
         bfgs_iterations = optimization%iterations
         if (.not. optimization%converged .and. warm_start%converged) &
            optimization = warm_start
         coordinates = optimization%parameters
         out%iterations = warm_start%iterations + bfgs_iterations
         out%converged = optimization%converged
      else
         allocate(coordinates(0))
         out%converged = .true.
      end if
      call evaluate(coordinates, fitted_specs, fitted_noise, beta, adjusted, &
         filtered, status)
      if (status /= 0) then
         out%info = 10 + status
         return
      end if
      out%inputs = fitted_specs
      out%noise = fitted_noise
      out%regression_coefficients = beta
      allocate(out%signal(size(output), input_count), &
         out%regression_signal(size(output)))
      out%signal = 0.0_dp
      do dynamic_input = 1, input_count
         out%signal(:, dynamic_input) = tfarima_filter( &
            input(:, dynamic_input), fitted_specs(dynamic_input)%transfer)
      end do
      out%regression_signal = 0.0_dp
      if (nreg > 0) out%regression_signal = matmul(regressors, beta)
      allocate(out%residuals(size(output)), out%standardized_residuals(size(output)), &
         out%fitted(size(output)))
      out%residuals = filtered%innovation(:, 1)
      out%fitted = output - out%residuals
      out%standardized_residuals = 0.0_dp
      do i = 1, size(output)
         if (filtered%innovation_cov(1, 1, i) <= 0.0_dp) cycle
         out%standardized_residuals(i) = out%residuals(i)/ &
            sqrt(filtered%innovation_cov(1, 1, i))
      end do
      out%log_likelihood = filtered%log_likelihood
      observations = count(ieee_is_finite(filtered%innovation(:, 1)))
      out%aic = -2.0_dp*out%log_likelihood + 2.0_dp*real(npar, dp)
      out%bic = -2.0_dp*out%log_likelihood + &
         log(real(max(1, observations), dp))*real(npar, dp)
      allocate(out%parameters(npar), out%covariance(npar, npar), &
         out%standard_error(npar))
      call natural_parameters(coordinates, out%parameters)
      out%covariance = 0.0_dp
      out%standard_error = 0.0_dp
      if (npar > 0) then
         hessian = finite_difference_hessian(objective, coordinates)
         call invert_matrix(hessian, inverse, status)
         if (status == 0) then
            allocate(jacobian(npar, npar))
            jacobian = 0.0_dp
            do i = 1, npar
               jacobian(i, i) = 1.0_dp
               if (variance_parameter(i)) &
                  jacobian(i, i) = 2.0_dp*out%parameters(i)
            end do
            out%covariance = matmul(matmul(jacobian, inverse), &
               transpose(jacobian))
            do i = 1, npar
               out%standard_error(i) = &
                  sqrt(max(0.0_dp, out%covariance(i, i)))
            end do
         end if
      end if
      if (.not. out%converged .and. npar > 0) out%info = 30

   contains

      pure function objective(values) result(value)
         ! Return the negative exact likelihood for transfer optimizer values.
         real(dp), intent(in) :: values(:)
         real(dp) :: value
         type(tfarima_transfer_spec_t), allocatable :: trial_specs(:)
         type(tfarima_ucarima_component_t) :: trial_noise
         type(kfs_filter_t) :: trial_filter
         real(dp), allocatable :: trial_beta(:), trial_series(:)
         integer :: evaluation_status

         call evaluate(values, trial_specs, trial_noise, trial_beta, trial_series, &
            trial_filter, evaluation_status)
         if (evaluation_status == 0 .and. &
            ieee_is_finite(trial_filter%log_likelihood)) then
            value = -trial_filter%log_likelihood
         else
            value = 1.0e30_dp + dot_product(values, values)
         end if
      end function objective

      pure subroutine evaluate(values, specs, noise, regression, transformed, &
         filter_result, evaluation_status)
         ! Decode all signals and evaluate their exact ARIMA noise likelihood.
         real(dp), intent(in) :: values(:)
         type(tfarima_transfer_spec_t), allocatable, intent(out) :: specs(:)
         type(tfarima_ucarima_component_t), intent(out) :: noise
         real(dp), allocatable, intent(out) :: regression(:), transformed(:)
         type(kfs_filter_t), intent(out) :: filter_result
         integer, intent(out) :: evaluation_status
         type(tfarima_ucarima_model_t) :: trial_model
         type(tfarima_ucarima_ssm_t) :: state_space
         real(dp), allocatable :: numerator(:), denominator(:), ar(:), ma(:)
         integer :: trial_position, trial_input, trial_coefficient

         evaluation_status = 1
         allocate(specs(input_count), regression(nreg))
         trial_position = nreg
         if (nreg > 0) regression = values(:nreg)
         transformed = output
         if (nreg > 0) transformed = transformed - matmul(regressors, regression)
         do trial_input = 1, input_count
            numerator = specifications(trial_input)%transfer%numerator
            denominator = specifications(trial_input)%transfer%denominator
            do trial_coefficient = 1, size(numerator)
               if (.not. specifications(trial_input)% &
                  estimate_numerator(trial_coefficient)) cycle
               trial_position = trial_position + 1
               numerator(trial_coefficient) = values(trial_position)
            end do
            do trial_coefficient = 1, size(denominator)
               if (.not. specifications(trial_input)% &
                  estimate_denominator(trial_coefficient)) cycle
               trial_position = trial_position + 1
               denominator(trial_coefficient) = values(trial_position)
            end do
            if (.not. tfarima_polynomial_admissible(denominator)) return
            specs(trial_input) = tfarima_transfer_spec(numerator, denominator, &
               specifications(trial_input)%transfer%delay, &
               specifications(trial_input)%estimate_numerator, &
               specifications(trial_input)%estimate_denominator)
            transformed = transformed - tfarima_filter( &
               input(:, trial_input), specs(trial_input)%transfer)
         end do
         ar = noise_model%ar_polynomial
         ma = noise_model%ma_polynomial
         do trial_coefficient = 1, size(noise_ar_mask)
            if (.not. noise_ar_mask(trial_coefficient)) cycle
            trial_position = trial_position + 1
            ar(trial_coefficient + 1) = values(trial_position)
         end do
         do trial_coefficient = 1, size(noise_ma_mask)
            if (.not. noise_ma_mask(trial_coefficient)) cycle
            trial_position = trial_position + 1
            ma(trial_coefficient + 1) = values(trial_position)
         end do
         if (.not. tfarima_polynomial_admissible(ar) .or. &
            .not. tfarima_polynomial_admissible(ma)) return
         noise%ar_polynomial = ar
         noise%difference_polynomial = noise_model%difference_polynomial
         noise%ma_polynomial = ma
         noise%innovation_variance = noise_model%innovation_variance
         if (fit_variance) then
            trial_position = trial_position + 1
            if (abs(values(trial_position)) > 350.0_dp) return
            noise%innovation_variance = exp(2.0_dp*values(trial_position))
         end if
         allocate(trial_model%components(1))
         trial_model%components(1) = noise
         trial_model%info = 0
         state_space = tfarima_ucarima_state_space(transformed, trial_model, tol)
         if (state_space%info /= 0) then
            evaluation_status = 2
            return
         end if
         if (state_space%diffuse) then
            filter_result = kfs_filter_diffuse(state_space%model, tol)
         else
            filter_result = kfs_filter(state_space%model, tol)
         end if
         if (filter_result%info /= 0) then
            evaluation_status = 3
            return
         end if
         evaluation_status = 0
      end subroutine evaluate

      pure subroutine natural_parameters(values, natural)
         ! Transform the final noise variance coordinate to its natural scale.
         real(dp), intent(in) :: values(:)
         real(dp), intent(out) :: natural(:)
         integer :: parameter

         natural = values
         do parameter = 1, size(values)
            if (variance_parameter(parameter)) &
               natural(parameter) = exp(2.0_dp*values(parameter))
         end do
      end subroutine natural_parameters

   end function tfarima_exact_transfer_fit

   pure function tfarima_select_transfer(output, input, specifications, &
      noise_model, regressors, initial_regression, estimate_noise_ar, &
      estimate_noise_ma, estimate_noise_variance, significance_level, &
      keep_regressors, keep_inputs, max_steps, max_iterations, tolerance) &
      result(out)
      ! Backward-select regressors and complete inputs using exact-fit p-values.
      real(dp), intent(in) :: output(:), input(:, :)
      type(tfarima_transfer_spec_t), intent(in) :: specifications(:)
      type(tfarima_ucarima_component_t), intent(in) :: noise_model
      real(dp), intent(in), optional :: regressors(:, :), initial_regression(:)
      logical, intent(in), optional :: estimate_noise_ar(:), estimate_noise_ma(:)
      logical, intent(in), optional :: estimate_noise_variance
      real(dp), intent(in), optional :: significance_level, tolerance
      logical, intent(in), optional :: keep_regressors(:), keep_inputs(:)
      integer, intent(in), optional :: max_steps, max_iterations
      type(tfarima_transfer_selection_t) :: out
      type(tfarima_transfer_spec_t), allocatable :: current_specs(:)
      type(tfarima_ucarima_component_t) :: current_noise
      logical, allocatable :: noise_ar_mask(:), noise_ma_mask(:)
      logical, allocatable :: protected_regressor(:), protected_input(:)
      real(dp), allocatable :: current_input(:, :), current_regressors(:, :)
      real(dp), allocatable :: current_beta(:)
      integer, allocatable :: regressor_index(:), input_index(:)
      integer, allocatable :: history_kind(:), history_index(:)
      real(dp) :: alpha, tol, p_value, worst_p
      integer :: nreg, ninput, step_limit, iteration_limit, step
      integer :: candidate_kind, candidate_position, original_index
      integer :: parameter_position, dynamic_input, status
      logical :: fit_variance

      nreg = 0
      if (present(regressors)) nreg = size(regressors, 2)
      ninput = size(specifications)
      alpha = 0.10_dp
      if (present(significance_level)) alpha = significance_level
      tol = 1.0e-6_dp
      if (present(tolerance)) tol = tolerance
      iteration_limit = 250
      if (present(max_iterations)) iteration_limit = max_iterations
      step_limit = nreg + ninput
      if (present(max_steps)) step_limit = max_steps
      if (size(input, 1) /= size(output) .or. size(input, 2) /= ninput .or. &
         alpha <= 0.0_dp .or. alpha >= 1.0_dp .or. tol <= 0.0_dp .or. &
         iteration_limit < 1 .or. step_limit < 0) then
         out%info = 1
         return
      end if
      if (.not. allocated(noise_model%ar_polynomial) .or. &
         .not. allocated(noise_model%ma_polynomial)) then
         out%info = 1
         return
      end if
      if (present(regressors)) then
         if (size(regressors, 1) /= size(output)) then
            out%info = 1
            return
         end if
      end if
      if (present(initial_regression)) then
         if (.not. present(regressors) .or. size(initial_regression) /= nreg) then
            out%info = 1
            return
         end if
      end if
      allocate(protected_regressor(nreg), protected_input(ninput))
      protected_regressor = .false.
      protected_input = .false.
      if (present(keep_regressors)) then
         if (size(keep_regressors) /= nreg) then
            out%info = 1
            return
         end if
         protected_regressor = keep_regressors
      end if
      if (present(keep_inputs)) then
         if (size(keep_inputs) /= ninput) then
            out%info = 1
            return
         end if
         protected_input = keep_inputs
      end if
      allocate(noise_ar_mask(size(noise_model%ar_polynomial) - 1), &
         noise_ma_mask(size(noise_model%ma_polynomial) - 1))
      noise_ar_mask = .false.
      noise_ma_mask = .false.
      if (present(estimate_noise_ar)) then
         if (size(estimate_noise_ar) /= size(noise_ar_mask)) then
            out%info = 1
            return
         end if
         noise_ar_mask = estimate_noise_ar
      end if
      if (present(estimate_noise_ma)) then
         if (size(estimate_noise_ma) /= size(noise_ma_mask)) then
            out%info = 1
            return
         end if
         noise_ma_mask = estimate_noise_ma
      end if
      fit_variance = .true.
      if (present(estimate_noise_variance)) fit_variance = estimate_noise_variance
      allocate(current_regressors(size(output), nreg), current_beta(nreg))
      if (nreg > 0) current_regressors = regressors
      current_beta = 0.0_dp
      if (present(initial_regression)) current_beta = initial_regression
      current_input = input
      current_specs = specifications
      current_noise = noise_model
      regressor_index = [(status, status=1,nreg)]
      input_index = [(status, status=1,ninput)]
      allocate(out%retained_regressor(nreg), out%retained_input(ninput))
      allocate(out%regressor_p_value(nreg), out%input_p_value(ninput))
      allocate(history_kind(nreg + ninput), history_index(nreg + ninput))
      out%retained_regressor = .true.
      out%retained_input = .true.
      out%regressor_p_value = 0.0_dp
      out%input_p_value = 0.0_dp
      do step = 0, step_limit
         out%fit = tfarima_exact_transfer_fit(output, current_input, &
            current_specs, current_noise, regressors=current_regressors, &
            initial_regression=current_beta, estimate_noise_ar=noise_ar_mask, &
            estimate_noise_ma=noise_ma_mask, &
            estimate_noise_variance=fit_variance, &
            max_iterations=iteration_limit, tolerance=tol)
         if (out%fit%info /= 0) then
            out%info = 10 + out%fit%info
            exit
         end if
         current_specs = out%fit%inputs
         current_noise = out%fit%noise
         current_beta = out%fit%regression_coefficients
         worst_p = alpha
         candidate_kind = 0
         candidate_position = 0
         do status = 1, size(regressor_index)
            original_index = regressor_index(status)
            p_value = coefficient_p_value(current_beta(status), &
               out%fit%standard_error(status))
            if (p_value < 0.0_dp) then
               out%info = 3
               exit
            end if
            out%regressor_p_value(original_index) = p_value
            if (protected_regressor(original_index) .or. p_value <= worst_p) cycle
            worst_p = p_value
            candidate_kind = tfarima_selection_regressor
            candidate_position = status
         end do
         if (out%info /= 0) exit
         parameter_position = size(regressor_index)
         do dynamic_input = 1, size(input_index)
            original_index = input_index(dynamic_input)
            if (current_specs(dynamic_input)%estimate_numerator(1)) then
               p_value = coefficient_p_value( &
                  out%fit%parameters(parameter_position + 1), &
                  out%fit%standard_error(parameter_position + 1))
               if (p_value < 0.0_dp) then
                  out%info = 3
                  exit
               end if
               out%input_p_value(original_index) = p_value
               if (.not. protected_input(original_index) .and. p_value > worst_p) then
                  worst_p = p_value
                  candidate_kind = tfarima_selection_input
                  candidate_position = dynamic_input
               end if
            end if
            parameter_position = parameter_position + &
               count(current_specs(dynamic_input)%estimate_numerator) + &
               count(current_specs(dynamic_input)%estimate_denominator)
         end do
         if (out%info /= 0) exit
         if (candidate_kind == 0) then
            out%converged = .true.
            exit
         end if
         if (step == step_limit) then
            out%info = 2
            exit
         end if
         out%steps = out%steps + 1
         history_kind(out%steps) = candidate_kind
         if (candidate_kind == tfarima_selection_regressor) then
            original_index = regressor_index(candidate_position)
            history_index(out%steps) = original_index
            out%retained_regressor(original_index) = .false.
            call omit_regressor(current_regressors, current_beta, &
               regressor_index, candidate_position)
         else
            original_index = input_index(candidate_position)
            history_index(out%steps) = original_index
            out%retained_input(original_index) = .false.
            call omit_transfer_input(current_input, current_specs, input_index, &
               candidate_position)
         end if
      end do
      allocate(out%removed_kind(out%steps), out%removed_index(out%steps))
      if (out%steps > 0) then
         out%removed_kind = history_kind(:out%steps)
         out%removed_index = history_index(:out%steps)
      end if

   contains

      pure real(dp) function coefficient_p_value(coefficient, standard_error) &
         result(value)
         ! Return a two-sided normal p-value or -1 when inference is unavailable.
         real(dp), intent(in) :: coefficient, standard_error

         if (.not. ieee_is_finite(coefficient) .or. &
            .not. ieee_is_finite(standard_error) .or. &
            standard_error <= 0.0_dp) then
            value = -1.0_dp
         else
            value = erfc(abs(coefficient/standard_error)/sqrt(2.0_dp))
         end if
      end function coefficient_p_value

      pure subroutine omit_regressor(matrix, beta, index, position)
         ! Remove one regressor column and its current coefficient.
         real(dp), allocatable, intent(inout) :: matrix(:, :), beta(:)
         integer, allocatable, intent(inout) :: index(:)
         integer, intent(in) :: position
         integer, allocatable :: retained(:)
         integer :: item

         retained = pack([(item, item=1,size(index))], &
            [(item /= position, item=1,size(index))])
         matrix = matrix(:, retained)
         beta = beta(retained)
         index = index(retained)
      end subroutine omit_regressor

      pure subroutine omit_transfer_input(matrix, specs, index, position)
         ! Remove one dynamic input, its specification, and original index.
         real(dp), allocatable, intent(inout) :: matrix(:, :)
         type(tfarima_transfer_spec_t), allocatable, intent(inout) :: specs(:)
         integer, allocatable, intent(inout) :: index(:)
         integer, intent(in) :: position
         integer, allocatable :: retained(:)
         integer :: item

         retained = pack([(item, item=1,size(index))], &
            [(item /= position, item=1,size(index))])
         matrix = matrix(:, retained)
         specs = specs(retained)
         index = index(retained)
      end subroutine omit_transfer_input

   end function tfarima_select_transfer

   pure function tfarima_exact_transfer_forecast(fit, output_history, &
      input_history, future_input, future_regressors, levels) result(out)
      ! Forecast a fitted transfer model for known future dynamic inputs.
      type(tfarima_exact_transfer_fit_t), intent(in) :: fit
      real(dp), intent(in) :: output_history(:), input_history(:, :)
      real(dp), intent(in) :: future_input(:, :)
      real(dp), intent(in), optional :: future_regressors(:, :), levels(:)
      type(tfarima_exact_transfer_forecast_t) :: out
      type(tfarima_forecast_t) :: noise_forecast
      real(dp), allocatable :: noise_history(:), complete_input(:), signal(:)
      real(dp) :: quantile
      integer :: horizon, dynamic_input, level, nreg

      horizon = size(future_input, 1)
      if (fit%info /= 0 .or. .not. allocated(fit%inputs) .or. horizon < 1 .or. &
         size(input_history, 1) /= size(output_history) .or. &
         size(input_history, 2) /= size(fit%inputs) .or. &
         size(future_input, 2) /= size(fit%inputs) .or. &
         .not. all(ieee_is_finite(output_history)) .or. &
         .not. all(ieee_is_finite(input_history)) .or. &
         .not. all(ieee_is_finite(future_input))) then
         out%info = 1
         return
      end if
      nreg = size(fit%regression_coefficients)
      if (nreg > 0) then
         if (.not. present(future_regressors)) then
            out%info = 2
            return
         end if
         if (any(shape(future_regressors) /= [horizon, nreg])) then
            out%info = 2
            return
         end if
      end if
      noise_history = output_history
      if (nreg > 0) then
         if (size(fit%regression_signal) /= size(output_history)) then
            out%info = 2
            return
         end if
         noise_history = noise_history - fit%regression_signal
      end if
      do dynamic_input = 1, size(fit%inputs)
         signal = tfarima_filter(input_history(:, dynamic_input), &
            fit%inputs(dynamic_input)%transfer)
         noise_history = noise_history - signal
      end do
      noise_forecast = tfarima_arima_forecast(noise_history, &
         fit%noise%ar_polynomial, fit%noise%difference_polynomial, &
         fit%noise%ma_polynomial, fit%noise%innovation_variance, horizon)
      if (noise_forecast%info /= 0) then
         out%info = 3
         return
      end if
      out%mean = noise_forecast%mean
      out%variance = noise_forecast%variance
      if (nreg > 0) out%mean = out%mean + &
         matmul(future_regressors, fit%regression_coefficients)
      do dynamic_input = 1, size(fit%inputs)
         complete_input = [input_history(:, dynamic_input), &
            future_input(:, dynamic_input)]
         signal = tfarima_filter(complete_input, &
            fit%inputs(dynamic_input)%transfer)
         out%mean = out%mean + signal(size(output_history) + 1:)
      end do
      if (present(levels)) then
         if (any(levels <= 0.0_dp) .or. any(levels >= 1.0_dp)) then
            out%info = 4
            return
         end if
         allocate(out%lower(horizon, size(levels)), &
            out%upper(horizon, size(levels)))
         do level = 1, size(levels)
            quantile = normal_quantile(0.5_dp + 0.5_dp*levels(level))
            out%lower(:, level) = out%mean - quantile*sqrt(out%variance)
            out%upper(:, level) = out%mean + quantile*sqrt(out%variance)
         end do
      end if
   end function tfarima_exact_transfer_forecast

   pure function tfarima_transfer_fit(output, input, delay, denominator_order, &
      numerator_order, noise_ar_polynomial, noise_ma_polynomial, &
      initial_parameters, include_mean, max_iterations, tolerance) result(out)
      ! Fit a rational transfer function by conditional Gaussian likelihood.
      real(dp), intent(in) :: output(:), input(:)
      integer, intent(in) :: delay, denominator_order, numerator_order
      real(dp), intent(in) :: noise_ar_polynomial(:), noise_ma_polynomial(:)
      real(dp), intent(in), optional :: initial_parameters(:), tolerance
      logical, intent(in), optional :: include_mean
      integer, intent(in), optional :: max_iterations
      type(tfarima_transfer_fit_t) :: out
      type(optimization_result_t) :: optimization, warm_start
      real(dp), allocatable :: initial(:), denominator(:), numerator(:), signal(:)
      real(dp), allocatable :: hessian(:, :), inverse(:, :)
      real(dp) :: tol, variance
      integer :: npar, selected_iterations, status, i
      logical :: with_mean

      with_mean = .false.
      if (present(include_mean)) with_mean = include_mean
      selected_iterations = 200
      if (present(max_iterations)) selected_iterations = max_iterations
      tol = 1.0e-6_dp
      if (present(tolerance)) tol = tolerance
      npar = numerator_order + 1 + denominator_order + merge(1, 0, with_mean)
      if (size(output) /= size(input) .or. size(output) < 4 .or. delay < 0 .or. &
         denominator_order < 0 .or. numerator_order < 0 .or. &
         size(noise_ar_polynomial) < 1 .or. size(noise_ma_polynomial) < 1 .or. &
         selected_iterations < 1 .or. tol <= 0.0_dp .or. &
         .not. all(ieee_is_finite(output)) .or. &
         .not. all(ieee_is_finite(input))) then
         out%info = 1
         return
      end if
      allocate(initial(npar))
      initial = 0.0_dp
      initial(1) = dot_product(output, input)/max(dot_product(input, input), &
         tiny(1.0_dp))
      if (with_mean) initial(npar) = sum(output)/real(size(output), dp)
      if (present(initial_parameters)) then
         if (size(initial_parameters) /= npar) then
            out%info = 1
            return
         end if
         initial = initial_parameters
      end if
      warm_start = nelder_mead_minimize(objective, initial, &
         max_iterations=selected_iterations, tolerance=tol, initial_step=0.05_dp)
      optimization = bfgs_minimize_fd(objective, warm_start%parameters, &
         max_iterations=selected_iterations, gradient_tolerance=tol)
      if (.not. optimization%converged .and. warm_start%converged) then
         optimization%parameters = warm_start%parameters
         optimization%objective = warm_start%objective
         optimization%converged = .true.
         optimization%info = 0
      end if
      out%parameters = optimization%parameters
      out%iterations = warm_start%iterations + optimization%iterations
      out%converged = optimization%converged
      call unpack_parameters(out%parameters, numerator_order, denominator_order, &
         with_mean, numerator, denominator, out%mean)
      out%transfer = tfarima_transfer(numerator, denominator, delay)
      signal = tfarima_filter(input, out%transfer)
      out%fitted = signal + out%mean
      out%residuals = conditional_residuals(output - out%fitted, &
         noise_ar_polynomial, noise_ma_polynomial)
      variance = dot_product(out%residuals, out%residuals)/real(size(output), dp)
      out%innovation_variance = variance
      out%log_likelihood = -0.5_dp*real(size(output), dp)* &
         (log(2.0_dp*acos(-1.0_dp)*variance) + 1.0_dp)
      out%aic = -2.0_dp*out%log_likelihood + 2.0_dp*real(npar + 1, dp)
      out%bic = -2.0_dp*out%log_likelihood + &
         log(real(size(output), dp))*real(npar + 1, dp)
      hessian = finite_difference_hessian(objective, out%parameters)
      call invert_matrix(hessian, inverse, status)
      allocate(out%covariance(npar, npar), out%standard_error(npar))
      out%covariance = 0.0_dp
      out%standard_error = 0.0_dp
      if (status == 0) then
         out%covariance = inverse
         do i = 1, npar
            if (inverse(i, i) > 0.0_dp) out%standard_error(i) = sqrt(inverse(i, i))
         end do
      end if
      if (.not. optimization%converged) out%info = optimization%info

   contains

      pure function objective(parameters) result(value)
         ! Return the profiled conditional Gaussian negative log likelihood.
         real(dp), intent(in) :: parameters(:)
         real(dp) :: value, candidate_mean, candidate_variance
         real(dp), allocatable :: candidate_numerator(:), candidate_denominator(:)
         real(dp), allocatable :: candidate_signal(:), candidate_residuals(:)
         type(tfarima_transfer_t) :: candidate_transfer

         call unpack_parameters(parameters, numerator_order, denominator_order, &
            with_mean, candidate_numerator, candidate_denominator, candidate_mean)
         if (.not. tfarima_polynomial_admissible(candidate_denominator)) then
            value = huge(1.0_dp)/100.0_dp
            return
         end if
         candidate_transfer = tfarima_transfer(candidate_numerator, &
            candidate_denominator, delay)
         candidate_signal = tfarima_filter(input, candidate_transfer)
         candidate_residuals = conditional_residuals( &
            output - candidate_signal - candidate_mean, noise_ar_polynomial, &
            noise_ma_polynomial)
         candidate_variance = dot_product(candidate_residuals, candidate_residuals)/ &
            real(size(output), dp)
         value = 0.5_dp*real(size(output), dp)* &
            (log(2.0_dp*acos(-1.0_dp)*max(candidate_variance, tiny(1.0_dp))) + 1.0_dp)
      end function objective

   end function tfarima_transfer_fit

   pure function tfarima_restricted_transfer_fit(output, input, delay, &
      numerator_model, denominator_model, noise_ar_polynomial, &
      noise_ma_polynomial, initial_gain, include_mean, max_iterations, &
      tolerance) result(out)
      ! Fit a transfer function whose lag polynomials have linear restrictions.
      real(dp), intent(in) :: output(:), input(:), initial_gain
      integer, intent(in) :: delay
      type(tfarima_lag_polynomial_t), intent(in) :: numerator_model
      type(tfarima_lag_polynomial_t), intent(in) :: denominator_model
      real(dp), intent(in) :: noise_ar_polynomial(:), noise_ma_polynomial(:)
      logical, intent(in), optional :: include_mean
      integer, intent(in), optional :: max_iterations
      real(dp), intent(in), optional :: tolerance
      type(tfarima_transfer_fit_t) :: out
      type(optimization_result_t) :: optimization, warm_start
      type(tfarima_lag_polynomial_t) :: fitted_numerator, fitted_denominator
      real(dp), allocatable :: initial(:), signal(:), hessian(:, :), inverse(:, :)
      real(dp) :: tol, gain, variance
      integer :: numerator_parameters, denominator_parameters, npar
      integer :: selected_iterations, status, i
      logical :: with_mean

      with_mean = .false.
      if (present(include_mean)) with_mean = include_mean
      selected_iterations = 200
      if (present(max_iterations)) selected_iterations = max_iterations
      tol = 1.0e-6_dp
      if (present(tolerance)) tol = tolerance
      if (numerator_model%info /= 0 .or. denominator_model%info /= 0 .or. &
         .not. allocated(numerator_model%parameters) .or. &
         .not. allocated(denominator_model%parameters)) then
         out%info = 1
         return
      end if
      numerator_parameters = size(numerator_model%parameters)
      denominator_parameters = size(denominator_model%parameters)
      npar = 1 + numerator_parameters + denominator_parameters + &
         merge(1, 0, with_mean)
      if (size(output) /= size(input) .or. size(output) < 4 .or. delay < 0 .or. &
         selected_iterations < 1 .or. tol <= 0.0_dp) then
         out%info = 1
         return
      end if
      allocate(initial(npar))
      initial(1) = initial_gain
      initial(2:1 + numerator_parameters) = numerator_model%parameters
      initial(2 + numerator_parameters:1 + numerator_parameters + &
         denominator_parameters) = denominator_model%parameters
      if (with_mean) initial(npar) = sum(output)/real(size(output), dp)
      warm_start = nelder_mead_minimize(objective, initial, &
         max_iterations=selected_iterations, tolerance=tol, initial_step=0.05_dp)
      optimization = bfgs_minimize_fd(objective, warm_start%parameters, &
         max_iterations=selected_iterations, gradient_tolerance=tol)
      if (.not. optimization%converged .and. warm_start%converged) then
         optimization%parameters = warm_start%parameters
         optimization%objective = warm_start%objective
         optimization%converged = .true.
         optimization%info = 0
      end if
      out%parameters = optimization%parameters
      out%iterations = warm_start%iterations + optimization%iterations
      out%converged = optimization%converged
      call decode_restricted(out%parameters, gain, fitted_numerator, &
         fitted_denominator, out%mean)
      out%transfer = tfarima_transfer(gain*fitted_numerator%expanded, &
         fitted_denominator%expanded, delay)
      signal = tfarima_filter(input, out%transfer)
      out%fitted = signal + out%mean
      out%residuals = conditional_residuals(output - out%fitted, &
         noise_ar_polynomial, noise_ma_polynomial)
      variance = dot_product(out%residuals, out%residuals)/real(size(output), dp)
      out%innovation_variance = variance
      out%log_likelihood = -0.5_dp*real(size(output), dp)* &
         (log(2.0_dp*acos(-1.0_dp)*variance) + 1.0_dp)
      out%aic = -2.0_dp*out%log_likelihood + 2.0_dp*real(npar + 1, dp)
      out%bic = -2.0_dp*out%log_likelihood + &
         log(real(size(output), dp))*real(npar + 1, dp)
      hessian = finite_difference_hessian(objective, out%parameters)
      call invert_matrix(hessian, inverse, status)
      allocate(out%covariance(npar, npar), out%standard_error(npar))
      out%covariance = 0.0_dp
      out%standard_error = 0.0_dp
      if (status == 0) then
         out%covariance = inverse
         do i = 1, npar
            if (inverse(i, i) > 0.0_dp) out%standard_error(i) = sqrt(inverse(i, i))
         end do
      end if
      if (.not. optimization%converged) out%info = optimization%info

   contains

      pure subroutine decode_restricted(parameters, candidate_gain, &
         candidate_numerator, candidate_denominator, candidate_mean)
         ! Decode gain, restricted operators, and optional mean.
         real(dp), intent(in) :: parameters(:)
         real(dp), intent(out) :: candidate_gain, candidate_mean
         type(tfarima_lag_polynomial_t), intent(out) :: candidate_numerator
         type(tfarima_lag_polynomial_t), intent(out) :: candidate_denominator
         integer :: first_denominator

         candidate_gain = parameters(1)
         candidate_numerator = tfarima_update_lag_polynomial(numerator_model, &
            parameters(2:1 + numerator_parameters))
         first_denominator = 2 + numerator_parameters
         candidate_denominator = tfarima_update_lag_polynomial(denominator_model, &
            parameters(first_denominator:first_denominator + &
            denominator_parameters - 1))
         candidate_mean = 0.0_dp
         if (with_mean) candidate_mean = parameters(size(parameters))
      end subroutine decode_restricted

      pure function objective(parameters) result(value)
         ! Return the restricted profiled conditional negative log likelihood.
         real(dp), intent(in) :: parameters(:)
         real(dp) :: value, candidate_gain, candidate_mean, candidate_variance
         real(dp), allocatable :: candidate_signal(:), candidate_residuals(:)
         type(tfarima_lag_polynomial_t) :: candidate_numerator
         type(tfarima_lag_polynomial_t) :: candidate_denominator
         type(tfarima_transfer_t) :: candidate_transfer

         call decode_restricted(parameters, candidate_gain, candidate_numerator, &
            candidate_denominator, candidate_mean)
         if (.not. tfarima_polynomial_admissible(candidate_denominator%expanded)) then
            value = huge(1.0_dp)/100.0_dp
            return
         end if
         candidate_transfer = tfarima_transfer( &
            candidate_gain*candidate_numerator%expanded, &
            candidate_denominator%expanded, delay)
         candidate_signal = tfarima_filter(input, candidate_transfer)
         candidate_residuals = conditional_residuals( &
            output - candidate_signal - candidate_mean, noise_ar_polynomial, &
            noise_ma_polynomial)
         candidate_variance = dot_product(candidate_residuals, candidate_residuals)/ &
            real(size(output), dp)
         value = 0.5_dp*real(size(output), dp)* &
            (log(2.0_dp*acos(-1.0_dp)*max(candidate_variance, tiny(1.0_dp))) + 1.0_dp)
      end function objective

   end function tfarima_restricted_transfer_fit

   pure function tfarima_outlier_response(outlier_type, observations, &
      ar_polynomial, difference_polynomial, ma_polynomial, decay) result(response)
      ! Return an outlier's unit response in the ARIMA residual domain.
      integer, intent(in) :: outlier_type, observations
      real(dp), intent(in) :: ar_polynomial(:), difference_polynomial(:)
      real(dp), intent(in) :: ma_polynomial(:)
      real(dp), intent(in), optional :: decay
      real(dp), allocatable :: response(:), pi_weight(:)
      real(dp) :: selected_decay
      integer :: lag

      selected_decay = 0.7_dp
      if (present(decay)) selected_decay = decay
      if (observations < 1 .or. outlier_type < tfarima_outlier_io .or. &
         outlier_type > tfarima_outlier_tc .or. selected_decay < 0.0_dp .or. &
         selected_decay >= 1.0_dp) then
         allocate(response(0))
         return
      end if
      allocate(response(observations))
      response = 0.0_dp
      if (outlier_type == tfarima_outlier_io) then
         response(1) = 1.0_dp
         return
      end if
      pi_weight = tfarima_polynomial_ratio( &
         tfarima_polynomial_multiply(ar_polynomial, difference_polynomial), &
         ma_polynomial, observations - 1)
      select case (outlier_type)
      case (tfarima_outlier_ao)
         response = pi_weight
      case (tfarima_outlier_ls)
         response(1) = pi_weight(1)
         do lag = 2, observations
            response(lag) = response(lag - 1) + pi_weight(lag)
         end do
      case (tfarima_outlier_tc)
         response(1) = 1.0_dp
         do lag = 2, observations
            response(lag) = selected_decay*response(lag - 1) + pi_weight(lag)
         end do
      end select
   end function tfarima_outlier_response

   pure function tfarima_detect_outliers(series, ar_polynomial, &
      difference_polynomial, ma_polynomial, enabled_types, candidate_positions, &
      mean_value, cutoff, retention_cutoff, tc_decay, max_passes, &
      log_transform) result(out)
      ! Detect and jointly refit tfarima AO, LS, TC, and IO effects.
      real(dp), intent(in) :: series(:), ar_polynomial(:), difference_polynomial(:)
      real(dp), intent(in) :: ma_polynomial(:)
      integer, intent(in) :: enabled_types(:)
      integer, intent(in), optional :: candidate_positions(:), max_passes
      real(dp), intent(in), optional :: mean_value, cutoff, retention_cutoff
      real(dp), intent(in), optional :: tc_decay
      logical, intent(in), optional :: log_transform
      type(tfarima_outlier_result_t) :: out
      real(dp), allocatable :: stationary(:), residual(:), original_residual(:)
      real(dp), allocatable :: aligned(:), response(:, :), selected_response(:, :)
      real(dp), allocatable :: effect(:), t_ratio(:), decay(:), xtx(:, :), inverse(:, :)
      real(dp), allocatable :: xty(:), beta(:), covariance(:, :)
      integer, allocatable :: candidates(:), selected_position(:), selected_type(:)
      logical, allocatable :: handled(:), retained(:)
      real(dp) :: mu, threshold, retain_threshold, selected_decay, scale
      real(dp) :: candidate_effect, candidate_t, best_effect, best_t, variance
      integer :: limit, pass, item, type_index, position, n, degree
      integer :: selected, status, i, k
      logical :: logarithm

      mu = 0.0_dp
      if (present(mean_value)) mu = mean_value
      threshold = 3.0_dp
      if (present(cutoff)) threshold = abs(cutoff)
      retain_threshold = 1.64_dp
      if (present(retention_cutoff)) retain_threshold = abs(retention_cutoff)
      selected_decay = 0.7_dp
      if (present(tc_decay)) selected_decay = tc_decay
      limit = 5
      if (present(max_passes)) limit = max_passes
      logarithm = .false.
      if (present(log_transform)) logarithm = log_transform
      n = size(series)
      degree = size(difference_polynomial) - 1
      if (n < 3 .or. size(enabled_types) < 1 .or. &
         any(enabled_types < tfarima_outlier_io) .or. &
         any(enabled_types > tfarima_outlier_tc) .or. threshold <= 0.0_dp .or. &
         retain_threshold <= 0.0_dp .or. selected_decay < 0.0_dp .or. &
         selected_decay >= 1.0_dp .or. limit < 1 .or. &
         has_duplicate_integer(enabled_types)) then
         out%info = 1
         return
      end if
      if (present(candidate_positions)) then
         if (size(candidate_positions) < 1 .or. any(candidate_positions < 1) .or. &
            any(candidate_positions > n) .or. &
            has_duplicate_integer(candidate_positions)) then
            out%info = 1
            return
         end if
         candidates = candidate_positions
      else
         allocate(candidates(n))
         candidates = [(i, i=1,n)]
      end if
      stationary = tfarima_difference(series, difference_polynomial, logarithm) - mu
      if (size(stationary) < 1) then
         out%info = 1
         return
      end if
      residual = conditional_residuals(stationary, ar_polynomial, ma_polynomial)
      allocate(aligned(n))
      aligned = 0.0_dp
      aligned(degree + 1:) = residual
      original_residual = aligned
      allocate(response(n, 4))
      do type_index = 1, 4
         response(:, type_index) = tfarima_outlier_response(type_index, n, &
            ar_polynomial, difference_polynomial, ma_polynomial, selected_decay)
      end do
      allocate(handled(size(candidates)), selected_position(size(candidates)), &
         selected_type(size(candidates)), effect(size(candidates)), &
         t_ratio(size(candidates)), decay(size(candidates)))
      handled = .false.
      selected_position = 0
      selected_type = 0
      effect = 0.0_dp
      t_ratio = 0.0_dp
      decay = 0.0_dp
      selected = 0
      do pass = 1, limit
         scale = sample_standard_deviation(aligned)
         if (scale <= tiny(1.0_dp)) exit
         k = 0
         do item = 1, size(candidates)
            if (handled(item)) cycle
            position = candidates(item)
            if (.not. present(candidate_positions)) then
               if (abs(aligned(position)/scale) <= threshold) cycle
            end if
            best_t = 0.0_dp
            best_effect = 0.0_dp
            type_index = 0
            do i = 1, size(enabled_types)
               call outlier_t_ratio(position, response(:, enabled_types(i)), &
                  aligned, candidate_effect, candidate_t)
               if (abs(candidate_t) > abs(best_t)) then
                  best_t = candidate_t
                  best_effect = candidate_effect
                  type_index = enabled_types(i)
               end if
            end do
            if (type_index == 0) cycle
            handled(item) = .true.
            selected = selected + 1
            selected_position(selected) = position
            selected_type(selected) = type_index
            effect(selected) = best_effect
            t_ratio(selected) = best_t
            if (type_index == tfarima_outlier_tc) decay(selected) = selected_decay
            if (abs(best_t) > retain_threshold) then
               aligned(position:) = aligned(position:) - &
                  best_effect*response(:n - position + 1, type_index)
               k = k + 1
            end if
         end do
         out%passes = pass
         if (k == 0) exit
      end do
      if (selected == 0) then
         allocate(out%position(0), out%outlier_type(0), out%effect(0), &
            out%t_ratio(0), out%decay(0))
         return
      end if
      allocate(selected_response(n, selected))
      selected_response = 0.0_dp
      do item = 1, selected
         position = selected_position(item)
         selected_response(position:, item) = &
            response(:n - position + 1, selected_type(item))
      end do
      xtx = matmul(transpose(selected_response), selected_response)
      xty = matmul(transpose(selected_response), original_residual)
      call invert_matrix(xtx, inverse, status)
      if (status /= 0) then
         out%info = 2
         return
      end if
      beta = matmul(inverse, xty)
      variance = max(dot_product(original_residual - matmul(selected_response, beta), &
         original_residual - matmul(selected_response, beta))/ &
         real(max(1, n - selected), dp), tiny(1.0_dp))
      covariance = variance*inverse
      allocate(retained(selected))
      do item = 1, selected
         t_ratio(item) = beta(item)/sqrt(max(covariance(item, item), tiny(1.0_dp)))
         effect(item) = beta(item)
         retained(item) = abs(t_ratio(item)) >= retain_threshold
      end do
      k = count(retained)
      allocate(out%position(k), out%outlier_type(k), out%effect(k), &
         out%t_ratio(k), out%decay(k))
      out%position = pack(selected_position(:selected), retained)
      out%outlier_type = pack(selected_type(:selected), retained)
      out%effect = pack(effect(:selected), retained)
      out%t_ratio = pack(t_ratio(:selected), retained)
      out%decay = pack(decay(:selected), retained)
   end function tfarima_detect_outliers

   pure function tfarima_month_lengths(start_year, start_month, observations) &
      result(lengths)
      ! Return calendar lengths for a sequence of monthly observations.
      integer, intent(in) :: start_year, start_month, observations
      real(dp), allocatable :: lengths(:)
      integer :: time, year, month

      if (start_month < 1 .or. start_month > 12 .or. observations < 1) then
         allocate(lengths(0))
         return
      end if
      allocate(lengths(observations))
      do time = 1, observations
         call monthly_period(start_year, start_month, time, year, month)
         lengths(time) = real(date_days_in_month(year, month), dp)
      end do
   end function tfarima_month_lengths

   pure function tfarima_weekday_counts(start_year, start_month, observations) &
      result(counts)
      ! Count Sunday-through-Saturday occurrences in monthly periods.
      integer, intent(in) :: start_year, start_month, observations
      real(dp), allocatable :: counts(:, :)
      type(date_t) :: current
      integer :: time, year, month, day, weekday

      if (start_month < 1 .or. start_month > 12 .or. observations < 1) then
         allocate(counts(0, 0))
         return
      end if
      allocate(counts(observations, 7))
      counts = 0.0_dp
      do time = 1, observations
         call monthly_period(start_year, start_month, time, year, month)
         do day = 1, date_days_in_month(year, month)
            current = date_t(year, month, day)
            weekday = modulo(date_day_of_week(current), 7) + 1
            counts(time, weekday) = counts(time, weekday) + 1.0_dp
         end do
      end do
   end function tfarima_weekday_counts

   pure function tfarima_leap_year_regressor(start_year, start_month, &
      observations, working_day_coding) result(regressor)
      ! Create tfarima's February leap-year or working-day contrast.
      integer, intent(in) :: start_year, start_month, observations
      logical, intent(in), optional :: working_day_coding
      real(dp), allocatable :: regressor(:)
      integer :: time, year, month
      logical :: working_day

      working_day = .false.
      if (present(working_day_coding)) working_day = working_day_coding
      if (start_month < 1 .or. start_month > 12 .or. observations < 1) then
         allocate(regressor(0))
         return
      end if
      allocate(regressor(observations))
      regressor = 0.0_dp
      do time = 1, observations
         call monthly_period(start_year, start_month, time, year, month)
         if (month /= 2) cycle
         if (working_day) then
            if (date_days_in_month(year, 2) == 29) then
               regressor(time) = 0.75_dp
            else
               regressor(time) = -0.25_dp
            end if
         else if (date_days_in_month(year, 2) == 29) then
            regressor(time) = 1.0_dp
         end if
      end do
   end function tfarima_leap_year_regressor

   pure function tfarima_easter_regressor(start_year, start_month, observations, &
      window_length, easter_monday) result(regressor)
      ! Aggregate the final pre-Easter days into monthly regressors.
      integer, intent(in) :: start_year, start_month, observations, window_length
      logical, intent(in), optional :: easter_monday
      real(dp), allocatable :: regressor(:)
      type(date_t) :: easter_date, window_start, month_start, month_end
      integer :: time, year, month, first_day, last_day, overlap
      logical :: monday

      monday = .false.
      if (present(easter_monday)) monday = easter_monday
      if (start_month < 1 .or. start_month > 12 .or. observations < 1 .or. &
         window_length < 1 .or. window_length > 22) then
         allocate(regressor(0))
         return
      end if
      allocate(regressor(observations))
      regressor = 0.0_dp
      do time = 1, observations
         call monthly_period(start_year, start_month, time, year, month)
         easter_date = date_easter(year)
         if (monday) easter_date = easter_date + 1
         window_start = easter_date - (window_length - 1)
         month_start = date_t(year, month, 1)
         month_end = date_t(year, month, date_days_in_month(year, month))
         if (.not. date_valid(easter_date)) cycle
         first_day = max(date_day_number(window_start), date_day_number(month_start))
         last_day = min(date_day_number(easter_date), date_day_number(month_end))
         overlap = max(0, last_day - first_day + 1)
         regressor(time) = real(overlap, dp)/real(window_length, dp)
      end do
   end function tfarima_easter_regressor

   pure function tfarima_calendar_regressors(start_year, start_month, &
      observations, form, reference_weekday, include_month_length, &
      include_leap_year, include_easter, easter_window, easter_monday) result(out)
      ! Construct tfarima monthly trading-day and holiday regressors.
      integer, intent(in) :: start_year, start_month, observations, form
      integer, intent(in), optional :: reference_weekday, easter_window
      logical, intent(in), optional :: include_month_length, include_leap_year
      logical, intent(in), optional :: include_easter, easter_monday
      type(tfarima_calendar_result_t) :: out
      real(dp), allocatable :: counts(:, :), excess(:, :), lengths(:)
      real(dp), allocatable :: leap(:), easter(:)
      integer :: reference, window, base_columns, columns, column, weekday
      logical :: month_length, leap_year, easter_effect, monday

      reference = 0
      if (present(reference_weekday)) reference = reference_weekday
      window = 4
      if (present(easter_window)) window = easter_window
      month_length = .true.
      if (present(include_month_length)) month_length = include_month_length
      leap_year = .true.
      if (present(include_leap_year)) leap_year = include_leap_year
      easter_effect = .false.
      if (present(include_easter)) easter_effect = include_easter
      monday = .false.
      if (present(easter_monday)) monday = easter_monday
      out%form = form
      out%reference_weekday = reference
      if (start_month < 1 .or. start_month > 12 .or. observations < 1 .or. &
         form < tfarima_calendar_dif .or. form > tfarima_calendar_null .or. &
         reference < 0 .or. reference > 6 .or. window < 1 .or. window > 22) then
         out%info = 1
         return
      end if
      select case (form)
      case (tfarima_calendar_dif, tfarima_calendar_td)
         base_columns = 6
         if (month_length .and. leap_year) then
            if (form == tfarima_calendar_td) then
               month_length = .false.
            else
               leap_year = .false.
            end if
         end if
      case (tfarima_calendar_td7)
         base_columns = 7
         month_length = .false.
         leap_year = .false.
      case (tfarima_calendar_td6)
         base_columns = 6
      case (tfarima_calendar_wd)
         base_columns = 1
         if (leap_year) month_length = .false.
      case default
         base_columns = 0
      end select
      columns = base_columns + merge(1, 0, month_length) + &
         merge(1, 0, leap_year) + merge(1, 0, easter_effect)
      allocate(out%values(observations, columns))
      out%values = 0.0_dp
      counts = tfarima_weekday_counts(start_year, start_month, observations)
      excess = counts - 4.0_dp
      column = 0
      select case (form)
      case (tfarima_calendar_dif, tfarima_calendar_td)
         do weekday = 1, 7
            if (weekday == reference + 1) cycle
            column = column + 1
            out%values(:, column) = excess(:, weekday) - excess(:, reference + 1)
         end do
      case (tfarima_calendar_td7)
         out%values(:, 1:7) = excess
         column = 7
      case (tfarima_calendar_td6)
         do weekday = 1, 7
            if (weekday == reference + 1) cycle
            column = column + 1
            out%values(:, column) = excess(:, weekday)
         end do
      case (tfarima_calendar_wd)
         column = 1
         out%values(:, 1) = sum(excess(:, 2:6), dim=2) - &
            2.5_dp*(excess(:, 1) + excess(:, 7))
      end select
      if (month_length) then
         column = column + 1
         lengths = tfarima_month_lengths(start_year, start_month, observations)
         out%values(:, column) = lengths - 28.0_dp
      end if
      if (leap_year) then
         column = column + 1
         leap = tfarima_leap_year_regressor(start_year, start_month, observations, &
            working_day_coding=.true.)
         out%values(:, column) = leap
      end if
      if (easter_effect) then
         column = column + 1
         easter = tfarima_easter_regressor(start_year, start_month, observations, &
            window, monday)
         out%values(:, column) = easter
      end if
   end function tfarima_calendar_regressors

   pure subroutine monthly_period(start_year, start_month, offset, year, month)
      ! Map a one-origin monthly offset to Gregorian year and month.
      integer, intent(in) :: start_year, start_month, offset
      integer, intent(out) :: year, month
      integer :: zero_origin

      zero_origin = start_month - 1 + offset - 1
      year = start_year + zero_origin/12
      month = modulo(zero_origin, 12) + 1
   end subroutine monthly_period

   pure subroutine outlier_t_ratio(position, response, residual, effect, t_ratio)
      ! Estimate one shifted outlier response and its conditional t-ratio.
      integer, intent(in) :: position
      real(dp), intent(in) :: response(:), residual(:)
      real(dp), intent(out) :: effect, t_ratio
      real(dp) :: sxx, rss
      integer :: length

      length = size(residual) - position + 1
      sxx = dot_product(response(:length), response(:length))
      if (sxx <= tiny(1.0_dp)) then
         effect = 0.0_dp
         t_ratio = 0.0_dp
         return
      end if
      effect = dot_product(residual(position:), response(:length))/sxx
      rss = sum(residual(:position - 1)**2) + &
         sum((residual(position:) - effect*response(:length))**2)
      rss = rss/real(max(1, size(residual) - 1), dp)
      t_ratio = effect/sqrt(max(rss/sxx, tiny(1.0_dp)))
   end subroutine outlier_t_ratio

   pure subroutine unpack_parameters(parameters, numerator_order, denominator_order, &
      include_mean, numerator, denominator, mean_value)
      ! Decode unrestricted transfer-function optimizer coordinates.
      real(dp), intent(in) :: parameters(:)
      integer, intent(in) :: numerator_order, denominator_order
      logical, intent(in) :: include_mean
      real(dp), allocatable, intent(out) :: numerator(:), denominator(:)
      real(dp), intent(out) :: mean_value
      integer :: first_denominator

      first_denominator = numerator_order + 2
      numerator = parameters(:numerator_order + 1)
      allocate(denominator(denominator_order + 1))
      denominator(1) = 1.0_dp
      if (denominator_order > 0) denominator(2:) = &
         parameters(first_denominator:first_denominator + denominator_order - 1)
      mean_value = 0.0_dp
      if (include_mean) mean_value = parameters(size(parameters))
   end subroutine unpack_parameters

   pure function simulate_arima_innovations(model, innovations, series_history, &
      innovation_history) result(series)
      ! Generate one ARIMA path from supplied scaled innovations and histories.
      type(tfarima_ucarima_component_t), intent(in) :: model
      real(dp), intent(in) :: innovations(:), series_history(:)
      real(dp), intent(in) :: innovation_history(:)
      real(dp), allocatable :: series(:), stationary(:), stationary_history(:)
      real(dp) :: value
      integer :: time, lag, index

      allocate(series(size(innovations)), stationary(size(innovations)))
      stationary_history = tfarima_difference(series_history, &
         model%difference_polynomial)
      do time = 1, size(innovations)
         value = model%ma_polynomial(1)*innovations(time)
         do lag = 1, size(model%ar_polynomial) - 1
            index = time - lag
            if (index >= 1) then
               value = value - model%ar_polynomial(lag + 1)*stationary(index)
            else if (size(stationary_history) + index >= 1) then
               value = value - model%ar_polynomial(lag + 1)* &
                  stationary_history(size(stationary_history) + index)
            end if
         end do
         do lag = 1, size(model%ma_polynomial) - 1
            index = time - lag
            if (index >= 1) then
               value = value + model%ma_polynomial(lag + 1)*innovations(index)
            else if (size(innovation_history) + index >= 1) then
               value = value + model%ma_polynomial(lag + 1)* &
                  innovation_history(size(innovation_history) + index)
            end if
         end do
         stationary(time) = value/model%ar_polynomial(1)
         value = stationary(time)
         do lag = 1, size(model%difference_polynomial) - 1
            index = time - lag
            if (index >= 1) then
               value = value - model%difference_polynomial(lag + 1)*series(index)
            else if (size(series_history) + index >= 1) then
               value = value - model%difference_polynomial(lag + 1)* &
                  series_history(size(series_history) + index)
            end if
         end do
         series(time) = value/model%difference_polynomial(1)
      end do
   end function simulate_arima_innovations

   pure logical function valid_prewhitening_model(model) result(valid)
      ! Check an ARIMA component before using it as a prewhitening filter.
      type(tfarima_ucarima_component_t), intent(in) :: model

      valid = model%info == 0 .and. allocated(model%ar_polynomial) .and. &
         allocated(model%difference_polynomial) .and. &
         allocated(model%ma_polynomial)
      if (.not. valid) return
      valid = size(model%ar_polynomial) > 0 .and. &
         size(model%difference_polynomial) > 0 .and. &
         size(model%ma_polynomial) > 0 .and. &
         abs(model%ar_polynomial(1)) > tiny(1.0_dp) .and. &
         abs(model%ma_polynomial(1)) > tiny(1.0_dp) .and. &
         all(ieee_is_finite(model%ar_polynomial)) .and. &
         all(ieee_is_finite(model%difference_polynomial)) .and. &
         all(ieee_is_finite(model%ma_polynomial))
   end function valid_prewhitening_model

   pure function prewhiten_series(series, model) result(residuals)
      ! Difference a series and apply a zero-presample ARMA inverse filter.
      real(dp), intent(in) :: series(:)
      type(tfarima_ucarima_component_t), intent(in) :: model
      real(dp), allocatable :: residuals(:), stationary(:)

      stationary = tfarima_difference(series, model%difference_polynomial)
      residuals = conditional_residuals(stationary, model%ar_polynomial, &
         model%ma_polynomial)
   end function prewhiten_series

   pure integer function first_significant_delay(diagnostic) result(delay)
      ! Return the first significant nonnegative prewhitened correlation lag.
      type(tfarima_prewhitened_ccf_t), intent(in) :: diagnostic
      integer :: position

      delay = -1
      do position = 1, size(diagnostic%lag)
         if (diagnostic%lag(position) < 0) cycle
         if (.not. diagnostic%significant(position)) cycle
         delay = diagnostic%lag(position)
         return
      end do
   end function first_significant_delay

   pure function conditional_residuals(series, ar_polynomial, ma_polynomial) &
      result(residuals)
      ! Compute zero-presample conditional ARMA residuals.
      real(dp), intent(in) :: series(:), ar_polynomial(:), ma_polynomial(:)
      real(dp), allocatable :: residuals(:)
      real(dp) :: value
      integer :: time, lag

      allocate(residuals(size(series)))
      residuals = 0.0_dp
      do time = 1, size(series)
         value = series(time)
         do lag = 1, min(size(ar_polynomial) - 1, time - 1)
            value = value + ar_polynomial(lag + 1)*series(time - lag)
         end do
         do lag = 1, min(size(ma_polynomial) - 1, time - 1)
            value = value - ma_polynomial(lag + 1)*residuals(time - lag)
         end do
         residuals(time) = value/ma_polynomial(1)
      end do
   end function conditional_residuals

   pure function trim_polynomial(polynomial, tolerance) result(trimmed)
      ! Remove negligible trailing high-lag coefficients.
      real(dp), intent(in) :: polynomial(:), tolerance
      real(dp), allocatable :: trimmed(:)
      integer :: last

      if (size(polynomial) < 1) then
         allocate(trimmed(0))
         return
      end if
      last = size(polynomial)
      do while (last > 1 .and. abs(polynomial(last)) <= tolerance)
         last = last - 1
      end do
      trimmed = polynomial(:last)
   end function trim_polynomial

   pure function polynomial_subtract(first, second, tolerance) result(difference)
      ! Subtract coefficient vectors and remove negligible trailing terms.
      real(dp), intent(in) :: first(:), second(:), tolerance
      real(dp), allocatable :: difference(:), work(:)
      integer :: length

      length = max(size(first), size(second))
      allocate(work(length))
      work = 0.0_dp
      work(:size(first)) = first
      work(:size(second)) = work(:size(second)) - second
      difference = trim_polynomial(work, tolerance)
   end function polynomial_subtract

   pure logical function polynomial_is_zero(polynomial, tolerance) result(is_zero)
      ! Test whether every polynomial coefficient is negligible.
      real(dp), intent(in) :: polynomial(:), tolerance

      is_zero = size(polynomial) < 1
      if (size(polynomial) > 0) is_zero = all(abs(polynomial) <= tolerance)
   end function polynomial_is_zero

   pure function polynomial_autocovariance(polynomial) result(covariance)
      ! Form nonnegative-lag coefficients of p(B)*p(F).
      real(dp), intent(in) :: polynomial(:)
      real(dp), allocatable :: covariance(:)
      integer :: lag

      allocate(covariance(size(polynomial)))
      do lag = 0, size(polynomial) - 1
         covariance(lag + 1) = dot_product( &
            polynomial(:size(polynomial) - lag), polynomial(lag + 1:))
      end do
   end function polynomial_autocovariance

   pure function wold_from_palindromic(palindromic) result(wold)
      ! Express a symmetric Laurent polynomial in x = B + B**(-1).
      real(dp), intent(in) :: palindromic(:)
      real(dp), allocatable :: wold(:), p0(:), p1(:), p2(:)
      integer :: degree, n

      n = size(palindromic)
      allocate(wold(n))
      wold = 0.0_dp
      if (n < 1) return
      wold(1) = palindromic(1)
      if (n == 1) return
      wold(2) = palindromic(2)
      allocate(p0(n), p1(n), p2(n))
      p0 = 0.0_dp
      p1 = 0.0_dp
      p0(1) = 2.0_dp
      p1(2) = 1.0_dp
      do degree = 2, n - 1
         p2 = 0.0_dp
         p2(2:degree + 1) = p1(:degree)
         p2 = p2 - p0
         wold = wold + palindromic(degree + 1)*p2
         p0 = p1
         p1 = p2
      end do
   end function wold_from_palindromic

   pure function palindromic_from_wold(wold) result(palindromic)
      ! Convert a polynomial in B + B**(-1) to nonnegative Laurent lags.
      real(dp), intent(in) :: wold(:)
      real(dp), allocatable :: palindromic(:)
      integer :: degree, selection, index, n

      n = size(wold)
      allocate(palindromic(n))
      palindromic = 0.0_dp
      if (n < 1) return
      palindromic(1) = wold(1)
      if (n == 1) return
      palindromic(2) = wold(2)
      do degree = 2, n - 1
         do selection = 0, degree/2
            index = degree - 2*selection + 1
            palindromic(index) = palindromic(index) + &
               wold(degree + 1)*binomial_coefficient(degree, selection)
         end do
      end do
   end function palindromic_from_wold

   pure real(dp) function polynomial_max_difference(first, second) result(error)
      ! Return the largest coefficient difference after zero padding.
      real(dp), intent(in) :: first(:), second(:)
      real(dp), allocatable :: left(:), right(:)
      integer :: length

      length = max(size(first), size(second))
      allocate(left(length), right(length))
      left = 0.0_dp
      right = 0.0_dp
      left(:size(first)) = first
      right(:size(second)) = second
      error = maxval(abs(left - right))
   end function polynomial_max_difference

   pure real(dp) function rational_spectral_minimum(numerator, denominator, &
      tolerance) result(minimum)
      ! Approximate a rational spectrum minimum over x in [-2, 2].
      real(dp), intent(in) :: numerator(:), denominator(:), tolerance
      real(dp) :: den, value, x
      integer :: point

      minimum = huge(1.0_dp)
      do point = 0, 10000
         x = -2.0_dp + 4.0_dp*real(point, dp)/10000.0_dp
         den = abs(polynomial_value(denominator, x))
         if (den <= tolerance) cycle
         value = polynomial_value(numerator, x)/den
         minimum = min(minimum, value)
      end do
   end function rational_spectral_minimum

   pure real(dp) function polynomial_value(polynomial, x) result(value)
      ! Evaluate an increasing-power real polynomial by Horner's rule.
      real(dp), intent(in) :: polynomial(:), x
      integer :: coefficient

      value = 0.0_dp
      do coefficient = size(polynomial), 1, -1
         value = value*x + polynomial(coefficient)
      end do
   end function polynomial_value

   pure real(dp) function binomial_coefficient(n, k) result(value)
      ! Return an integer binomial coefficient in double precision.
      integer, intent(in) :: n, k
      integer :: i, selected

      selected = min(k, n - k)
      value = 1.0_dp
      do i = 1, selected
         value = value*real(n - selected + i, dp)/real(i, dp)
      end do
   end function binomial_coefficient

   pure subroutine stationary_state_covariance(transition, loading, variance, &
      tolerance, covariance, converged)
      ! Solve a stable discrete Lyapunov equation by covariance iteration.
      real(dp), intent(in) :: transition(:, :), loading(:), variance, tolerance
      real(dp), allocatable, intent(out) :: covariance(:, :)
      logical, intent(out) :: converged
      real(dp), allocatable :: next_covariance(:, :), noise(:, :)
      real(dp) :: scale
      integer :: iteration, states

      states = size(loading)
      allocate(covariance(states, states), &
         next_covariance(states, states), noise(states, states))
      noise = variance*spread(loading, 2, states)*spread(loading, 1, states)
      covariance = 0.0_dp
      converged = .false.
      do iteration = 1, 100000
         next_covariance = matmul(matmul(transition, covariance), &
            transpose(transition)) + noise
         next_covariance = 0.5_dp*(next_covariance + transpose(next_covariance))
         scale = max(1.0_dp, maxval(abs(next_covariance)))
         if (maxval(abs(next_covariance - covariance)) <= tolerance*scale) then
            converged = .true.
            covariance = next_covariance
            return
         end if
         covariance = next_covariance
      end do
   end subroutine stationary_state_covariance

   pure real(dp) function sample_standard_deviation(values) result(standard_deviation)
      ! Return the sample standard deviation of a real vector.
      real(dp), intent(in) :: values(:)
      real(dp) :: mean_value

      if (size(values) < 2) then
         standard_deviation = 0.0_dp
         return
      end if
      mean_value = sum(values)/real(size(values), dp)
      standard_deviation = sqrt(sum((values - mean_value)**2)/ &
         real(size(values) - 1, dp))
   end function sample_standard_deviation

   pure logical function has_duplicate_integer(values) result(duplicate)
      ! Report whether an integer vector contains repeated values.
      integer, intent(in) :: values(:)
      integer :: i, j

      duplicate = .false.
      do i = 1, size(values) - 1
         do j = i + 1, size(values)
            if (values(i) == values(j)) then
               duplicate = .true.
               return
            end if
         end do
      end do
   end function has_duplicate_integer

   pure elemental real(dp) function falling_factorial(value, order) result(product)
      ! Return value times its requested descending predecessors.
      integer, intent(in) :: value, order
      integer :: i

      product = 1.0_dp
      do i = 0, order - 1
         product = product*real(value - i, dp)
      end do
   end function falling_factorial

end module tfarima_mod
