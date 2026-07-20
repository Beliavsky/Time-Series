! SPDX-License-Identifier: GPL-3.0-or-later
! SPDX-FileComment: Algorithms translated from the R tseriesTARMA package.
module tseriestarma_mod
   !! Two-regime threshold autoregressive moving-average models.
   use kind_mod, only: dp
   use linalg_mod, only: symmetric_pseudoinverse
   use optimization_mod, only: optimization_result_t, bfgs_minimize, &
      bfgs_minimize_fd, finite_difference_hessian
   use astsa_mod, only: astsa_sarima_fit_t, astsa_arma_diagnostic_t, &
      sarima_fit, arma_check
   use random_mod, only: random_uniform, random_standard_normal, &
      random_standard_normal_matrix
   use stats_mod, only: sorted, quantile
   use tseriestarma_tables_mod, only: &
      tseriestarma_unit_root_critical_values, &
      tseriestarma_unit_root_critical_probabilities, &
      tseriestarma_garch_critical_values, &
      tseriestarma_garch_critical_probabilities
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, &
      ieee_quiet_nan
   use, intrinsic :: iso_fortran_env, only: output_unit
   implicit none
   private

   integer, parameter, public :: tseriestarma_innovation_normal = 1
   integer, parameter, public :: tseriestarma_innovation_student = 2
   integer, parameter, public :: tseriestarma_bootstrap_iid = 1
   integer, parameter, public :: tseriestarma_bootstrap_hansen = 2
   integer, parameter, public :: tseriestarma_bootstrap_rademacher = 3
   integer, parameter, public :: tseriestarma_bootstrap_normal = 4

   type, public :: tseriestarma_model_t
      !! Parameters of a sparse two-regime TARMA model.
      real(dp), allocatable :: phi_lower(:)
      real(dp), allocatable :: phi_upper(:)
      real(dp), allocatable :: theta_lower(:)
      real(dp), allocatable :: theta_upper(:)
      integer, allocatable :: ar_lags_lower(:)
      integer, allocatable :: ar_lags_upper(:)
      integer, allocatable :: ma_lags_lower(:)
      integer, allocatable :: ma_lags_upper(:)
      real(dp) :: threshold = 0.0_dp
      real(dp) :: innovation_sd_lower = 1.0_dp
      real(dp) :: innovation_sd_upper = 1.0_dp
      integer :: delay = 1
      integer :: info = 0
   end type tseriestarma_model_t

   type, public :: tseriestarma_evaluation_t
      !! Conditional fitted values, innovations, regimes, and objective value.
      real(dp), allocatable :: fitted(:)
      real(dp), allocatable :: residuals(:)
      integer, allocatable :: regime(:)
      real(dp) :: rss = huge(1.0_dp)
      integer :: start_index = 0
      integer :: effective_observations = 0
      integer :: info = 0
   end type tseriestarma_evaluation_t

   type, public :: tseriestarma_simulation_t
      !! Simulated observations, realized innovations, and regime path.
      type(tseriestarma_model_t) :: model
      real(dp), allocatable :: observations(:)
      real(dp), allocatable :: innovations(:)
      integer, allocatable :: regime(:)
      integer :: burnin = 0
      integer :: info = 0
   end type tseriestarma_simulation_t

   type, public :: tseriestarma_fit_t
      !! Conditional least-squares TARMA fit and threshold profile.
      type(tseriestarma_model_t) :: model
      type(tseriestarma_evaluation_t) :: evaluation
      real(dp), allocatable :: parameters(:)
      real(dp), allocatable :: threshold_values(:)
      real(dp), allocatable :: profile_rss(:)
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: standard_errors(:)
      real(dp), allocatable :: weights(:)
      real(dp) :: sigma2 = huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      real(dp) :: bic = huge(1.0_dp)
      integer :: iterations = 0
      integer :: optimizer_info = 0
      integer :: irls_iterations = 0
      integer :: innovation_family = tseriestarma_innovation_normal
      real(dp) :: robustness_alpha = 0.0_dp
      real(dp) :: student_degrees = 0.0_dp
      integer :: info = 0
      logical :: converged = .false.
      logical :: robust = .false.
   end type tseriestarma_fit_t

   type, public :: tseriestarma_derivatives_t
      !! Conditional innovations and their parameter derivatives.
      real(dp), allocatable :: residuals(:)
      real(dp), allocatable :: derivative(:, :)
      integer :: start_index = 0
      integer :: info = 0
   end type tseriestarma_derivatives_t

   type, public :: tseriestarma_forecast_t
      !! Point path and simulation-based forecast distribution.
      real(dp), allocatable :: point(:)
      real(dp), allocatable :: mean(:)
      real(dp), allocatable :: standard_deviation(:)
      real(dp), allocatable :: probabilities(:)
      real(dp), allocatable :: quantiles(:, :)
      real(dp), allocatable :: paths(:, :)
      integer, allocatable :: point_regime(:)
      integer, allocatable :: regimes(:, :)
      integer :: info = 0
   end type tseriestarma_forecast_t

   type, public :: tseriestarma_fit2_t
      !! Exact common-MA threshold-regression fit corresponding to TARMA.fit2.
      type(astsa_sarima_fit_t) :: noise_fit
      type(tseriestarma_model_t) :: forecast_model
      real(dp), allocatable :: common_ar(:)
      real(dp), allocatable :: lower_ar(:)
      real(dp), allocatable :: upper_ar(:)
      real(dp), allocatable :: regular_ma(:)
      real(dp), allocatable :: seasonal_ma(:)
      real(dp), allocatable :: regression_coefficients(:)
      integer, allocatable :: common_ar_lags(:)
      integer, allocatable :: lower_ar_lags(:)
      integer, allocatable :: upper_ar_lags(:)
      real(dp), allocatable :: threshold_values(:)
      real(dp), allocatable :: profile_aic(:)
      real(dp), allocatable :: fitted(:)
      real(dp), allocatable :: residuals(:)
      integer, allocatable :: regime(:)
      real(dp) :: intercept_lower = 0.0_dp
      real(dp) :: intercept_upper = 0.0_dp
      real(dp) :: threshold = 0.0_dp
      real(dp) :: rss = huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      real(dp) :: bic = huge(1.0_dp)
      integer :: delay = 1
      integer :: period = 0
      integer :: start_index = 0
      integer :: info = 0
      logical :: regime_intercepts = .true.
      logical :: external_threshold = .false.
   end type tseriestarma_fit2_t

   type, public :: tseriestarma_tar_test_t
      !! Supremum LM test results for an AR null against a TAR alternative.
      real(dp), allocatable :: test_values(:, :)
      real(dp), allocatable :: threshold_values(:)
      real(dp), allocatable :: coefficients(:)
      real(dp), allocatable :: residuals(:)
      real(dp) :: statistic(2) = 0.0_dp
      real(dp) :: maximizing_threshold(2) = 0.0_dp
      real(dp) :: lower_regime_proportion(2) = 0.0_dp
      real(dp) :: sigma2 = 0.0_dp
      real(dp) :: lower_probability = 0.25_dp
      real(dp) :: upper_probability = 0.75_dp
      integer :: ar_order = 0
      integer :: delay = 1
      integer :: degrees_of_freedom = 0
      integer :: effective_observations = 0
      integer :: info = 0
   end type tseriestarma_tar_test_t

   type, extends(tseriestarma_tar_test_t), public :: tseriestarma_tarma_test_t
      !! Supremum LM test results for an ARMA null against a TARMA alternative.
      type(astsa_sarima_fit_t) :: noise_fit
      integer :: ma_order = 0
      logical :: ma_fixed = .true.
   end type tseriestarma_tarma_test_t

   type, extends(tseriestarma_tar_test_t), public :: tseriestarma_tar_bootstrap_t
      !! Bootstrap inference for the AR-versus-TAR supremum LM statistics.
      real(dp), allocatable :: bootstrap_statistics(:, :)
      real(dp) :: p_values(2) = 1.0_dp
      real(dp) :: p_value = 1.0_dp
      integer :: bootstrap_type = tseriestarma_bootstrap_iid
      integer :: replications = 0
      integer :: burnin = 0
   end type tseriestarma_tar_bootstrap_t

   type, public :: tseriestarma_unit_root_test_t
      !! Supremum LM unit-root test of an IMA(1) null against TARMA(1,1).
      type(astsa_sarima_fit_t) :: noise_fit
      real(dp), allocatable :: test_values(:)
      real(dp), allocatable :: threshold_values(:)
      real(dp), allocatable :: residuals(:)
      real(dp) :: critical_probabilities(4) = 0.0_dp
      real(dp) :: critical_values(4) = 0.0_dp
      real(dp) :: statistic = 0.0_dp
      real(dp) :: maximizing_threshold = 0.0_dp
      real(dp) :: sigma2 = 0.0_dp
      real(dp) :: ma_coefficient = 0.0_dp
      real(dp) :: lower_probability = 0.25_dp
      real(dp) :: upper_probability = 0.75_dp
      integer :: effective_observations = 0
      integer :: info = 0
   end type tseriestarma_unit_root_test_t

   type, extends(tseriestarma_unit_root_test_t), public :: &
      tseriestarma_unit_root_bootstrap_t
      !! Bootstrap inference for the IMA-versus-TARMA unit-root test.
      real(dp), allocatable :: bootstrap_statistics(:)
      real(dp) :: p_value = 1.0_dp
      integer :: bootstrap_type = tseriestarma_bootstrap_rademacher
      integer :: replications = 0
   end type tseriestarma_unit_root_bootstrap_t

   type, public :: tseriestarma_arma_garch_fit_t
      !! Joint Gaussian ARMA-GARCH null-model fit used by TARMAGARCH.test.
      real(dp), allocatable :: ar(:)
      real(dp), allocatable :: ma(:)
      real(dp), allocatable :: arch(:)
      real(dp), allocatable :: garch(:)
      real(dp), allocatable :: residuals(:)
      real(dp), allocatable :: conditional_variance(:)
      real(dp) :: mean = 0.0_dp
      real(dp) :: variance_intercept = 0.0_dp
      real(dp) :: log_likelihood = -huge(1.0_dp)
      integer :: iterations = 0
      integer :: optimizer_info = 0
      integer :: info = 0
      logical :: converged = .false.
      logical :: causal = .false.
      logical :: invertible = .false.
      logical :: variance_stationary = .false.
   end type tseriestarma_arma_garch_fit_t

   type, public :: tseriestarma_garch_test_t
      !! GARCH-aware supremum LM test of ARMA against TARMA nonlinearity.
      type(tseriestarma_arma_garch_fit_t) :: null_fit
      real(dp), allocatable :: test_values(:)
      real(dp), allocatable :: threshold_values(:)
      real(dp) :: critical_probabilities(3) = 0.0_dp
      real(dp) :: critical_values(3) = 0.0_dp
      real(dp) :: statistic = 0.0_dp
      real(dp) :: maximizing_threshold = 0.0_dp
      real(dp) :: lower_regime_proportion = 0.0_dp
      real(dp) :: standardized_variance = 0.0_dp
      real(dp) :: lower_probability = 0.25_dp
      real(dp) :: upper_probability = 0.75_dp
      integer :: ar_order = 1
      integer :: ma_order = 1
      integer :: arch_order = 1
      integer :: garch_order = 1
      integer :: delay = 1
      integer :: degrees_of_freedom = 0
      integer :: effective_observations = 0
      integer :: info = 0
   end type tseriestarma_garch_test_t

   interface display
      module procedure display_tseriestarma_model
      module procedure display_tseriestarma_evaluation
      module procedure display_tseriestarma_simulation
      module procedure display_tseriestarma_fit
      module procedure display_tseriestarma_forecast
      module procedure display_tseriestarma_fit2
      module procedure display_tseriestarma_tar_test
      module procedure display_tseriestarma_tarma_test
      module procedure display_tseriestarma_tar_bootstrap
      module procedure display_tseriestarma_unit_root_test
      module procedure display_tseriestarma_unit_root_bootstrap
      module procedure display_tseriestarma_arma_garch_fit
      module procedure display_tseriestarma_garch_test
   end interface display

   interface tseriestarma_coefficients
      module procedure tseriestarma_model_coefficients
      module procedure tseriestarma_fit_coefficients
      module procedure tseriestarma_fit2_coefficients
      module procedure tseriestarma_arma_garch_coefficients
   end interface tseriestarma_coefficients

   interface tseriestarma_covariance
      module procedure tseriestarma_fit_covariance
      module procedure tseriestarma_fit2_covariance
   end interface tseriestarma_covariance

   interface tseriestarma_residuals
      module procedure tseriestarma_fit_residuals
      module procedure tseriestarma_fit2_residuals
      module procedure tseriestarma_arma_garch_residuals
   end interface tseriestarma_residuals

   interface tseriestarma_fitted_values
      module procedure tseriestarma_fit_fitted_values
      module procedure tseriestarma_fit2_fitted_values
   end interface tseriestarma_fitted_values

   public :: tseriestarma_model, tseriestarma_evaluate
   public :: tseriestarma_simulate_from_standard, tseriestarma_simulate
   public :: tseriestarma_initial_parameters, tseriestarma_fit
   public :: tseriestarma_residual_derivatives
   public :: tseriestarma_least_squares_gradient
   public :: tseriestarma_robust_weights, tseriestarma_robust_fit
   public :: tseriestarma_fit2
   public :: tseriestarma_forecast_from_standard, tseriestarma_forecast
   public :: tseriestarma_tar_test
   public :: tseriestarma_tarma_test
   public :: tseriestarma_tar_bootstrap_from_random
   public :: tseriestarma_tar_bootstrap
   public :: tseriestarma_unit_root_test
   public :: tseriestarma_unit_root_bootstrap_from_random
   public :: tseriestarma_unit_root_bootstrap
   public :: tseriestarma_unit_root_critical_values
   public :: tseriestarma_unit_root_critical_probabilities
   public :: tseriestarma_arma_garch_fit
   public :: tseriestarma_garch_test
   public :: tseriestarma_garch_critical_values
   public :: tseriestarma_garch_critical_probabilities
   public :: display, tseriestarma_coefficients, tseriestarma_covariance
   public :: tseriestarma_residuals, tseriestarma_fitted_values
   public :: display_tseriestarma_model, display_tseriestarma_evaluation
   public :: display_tseriestarma_simulation, display_tseriestarma_fit
   public :: display_tseriestarma_forecast, display_tseriestarma_fit2
   public :: display_tseriestarma_tar_test, display_tseriestarma_tarma_test
   public :: display_tseriestarma_tar_bootstrap
   public :: display_tseriestarma_unit_root_test
   public :: display_tseriestarma_unit_root_bootstrap
   public :: display_tseriestarma_arma_garch_fit
   public :: display_tseriestarma_garch_test
   public :: tseriestarma_model_coefficients, tseriestarma_fit_coefficients
   public :: tseriestarma_fit2_coefficients
   public :: tseriestarma_arma_garch_coefficients
   public :: tseriestarma_fit_covariance, tseriestarma_fit2_covariance
   public :: tseriestarma_fit_residuals, tseriestarma_fit2_residuals
   public :: tseriestarma_arma_garch_residuals
   public :: tseriestarma_fit_fitted_values
   public :: tseriestarma_fit2_fitted_values

contains

   subroutine display_tseriestarma_model(model, unit)
      !! Display a concise TARMA model specification.
      type(tseriestarma_model_t), intent(in) :: model !! TARMA model to display.
      integer, intent(in), optional :: unit !! Output unit; defaults to standard output.
      integer :: destination

      destination = output_unit
      if (present(unit)) destination = unit
      write(destination, '(a)') 'Two-regime TARMA model'
      write(destination, '(a, i0)') 'Status: ', model%info
      write(destination, '(a, es14.6)') 'Threshold: ', model%threshold
      write(destination, '(a, i0)') 'Threshold delay: ', model%delay
      call display_real_vector(destination, 'Lower-regime AR coefficients:', &
         model%phi_lower)
      call display_real_vector(destination, 'Upper-regime AR coefficients:', &
         model%phi_upper)
      call display_real_vector(destination, 'Lower-regime MA coefficients:', &
         model%theta_lower)
      call display_real_vector(destination, 'Upper-regime MA coefficients:', &
         model%theta_upper)
      write(destination, '(a, es14.6)') 'Lower innovation SD: ', &
         model%innovation_sd_lower
      write(destination, '(a, es14.6)') 'Upper innovation SD: ', &
         model%innovation_sd_upper
   end subroutine display_tseriestarma_model

   subroutine display_tseriestarma_evaluation(evaluation, unit, print_obs)
      !! Display a conditional TARMA evaluation and optionally its series arrays.
      type(tseriestarma_evaluation_t), intent(in) :: evaluation !! Evaluation to display.
      integer, intent(in), optional :: unit !! Output unit; defaults to standard output.
      logical, intent(in), optional :: print_obs !! Whether to print fitted values and residuals.
      integer :: destination, index
      logical :: show_observations

      destination = output_unit
      if (present(unit)) destination = unit
      show_observations = .false.
      if (present(print_obs)) show_observations = print_obs
      write(destination, '(a)') 'TARMA conditional evaluation'
      write(destination, '(a, i0)') 'Status: ', evaluation%info
      write(destination, '(a, i0)') 'First usable index: ', evaluation%start_index
      write(destination, '(a, i0)') 'Effective observations: ', &
         evaluation%effective_observations
      write(destination, '(a, es14.6)') 'Residual sum of squares: ', evaluation%rss
      if (show_observations .and. allocated(evaluation%fitted) .and. &
         allocated(evaluation%residuals)) then
         write(destination, '(a)') 'Index, fitted value, residual:'
         do index = 1, size(evaluation%fitted)
            write(destination, '(i8, 2(1x, es14.6))') index, &
               evaluation%fitted(index), evaluation%residuals(index)
         end do
      end if
   end subroutine display_tseriestarma_evaluation

   subroutine display_tseriestarma_simulation(simulation, unit, print_obs)
      !! Display a TARMA simulation specification and optionally its observations.
      type(tseriestarma_simulation_t), intent(in) :: simulation !! Simulation to display.
      integer, intent(in), optional :: unit !! Output unit; defaults to standard output.
      logical, intent(in), optional :: print_obs !! Whether to print simulated observations.
      integer :: destination, index
      logical :: show_observations

      destination = output_unit
      if (present(unit)) destination = unit
      show_observations = .false.
      if (present(print_obs)) show_observations = print_obs
      write(destination, '(a)') 'TARMA simulation'
      write(destination, '(a, i0)') 'Status: ', simulation%info
      write(destination, '(a, i0)') 'Burn-in observations: ', simulation%burnin
      if (allocated(simulation%observations)) &
         write(destination, '(a, i0)') 'Observations: ', size(simulation%observations)
      call display_tseriestarma_model(simulation%model, destination)
      if (show_observations .and. allocated(simulation%observations)) then
         write(destination, '(a)') 'Index, observation, innovation, regime:'
         do index = 1, size(simulation%observations)
            write(destination, '(i8, 2(1x, es14.6), 1x, i0)') index, &
               simulation%observations(index), simulation%innovations(index), &
               simulation%regime(index)
         end do
      end if
   end subroutine display_tseriestarma_simulation

   subroutine display_tseriestarma_fit(fit, unit, print_obs)
      !! Display a conditional least-squares TARMA fit.
      type(tseriestarma_fit_t), intent(in) :: fit !! Fitted TARMA model to display.
      integer, intent(in), optional :: unit !! Output unit; defaults to standard output.
      logical, intent(in), optional :: print_obs !! Whether to print fitted values and residuals.
      integer :: destination

      destination = output_unit
      if (present(unit)) destination = unit
      write(destination, '(a)') 'Conditional least-squares TARMA fit'
      write(destination, '(a, i0)') 'Status: ', fit%info
      write(destination, '(a, l1)') 'Converged: ', fit%converged
      write(destination, '(a, es14.6)') 'Innovation variance: ', fit%sigma2
      write(destination, '(a, es14.6)') 'AIC: ', fit%aic
      write(destination, '(a, es14.6)') 'BIC: ', fit%bic
      call display_tseriestarma_model(fit%model, destination)
      call display_real_vector(destination, 'Coefficients:', fit%parameters)
      if (present(print_obs)) then
         call display_tseriestarma_evaluation(fit%evaluation, destination, print_obs)
      end if
   end subroutine display_tseriestarma_fit

   subroutine display_tseriestarma_forecast(forecast, unit, print_obs)
      !! Display a TARMA forecast summary and optionally its simulated paths.
      type(tseriestarma_forecast_t), intent(in) :: forecast !! Forecast to display.
      integer, intent(in), optional :: unit !! Output unit; defaults to standard output.
      logical, intent(in), optional :: print_obs !! Whether to print simulated forecast paths.
      integer :: destination, row
      logical :: show_observations

      destination = output_unit
      if (present(unit)) destination = unit
      show_observations = .false.
      if (present(print_obs)) show_observations = print_obs
      write(destination, '(a)') 'TARMA forecast'
      write(destination, '(a, i0)') 'Status: ', forecast%info
      call display_real_vector(destination, 'Point forecast:', forecast%point)
      call display_real_vector(destination, 'Simulation mean:', forecast%mean)
      call display_real_vector(destination, 'Simulation standard deviation:', &
         forecast%standard_deviation)
      if (show_observations .and. allocated(forecast%paths)) then
         write(destination, '(a)') 'Simulated forecast paths:'
         do row = 1, size(forecast%paths, 1)
            write(destination, '(*(es14.6, 1x))') forecast%paths(row, :)
         end do
      end if
   end subroutine display_tseriestarma_forecast

   subroutine display_tseriestarma_fit2(fit, unit, print_obs)
      !! Display an exact common-MA threshold-regression fit.
      type(tseriestarma_fit2_t), intent(in) :: fit !! Exact TARMA fit to display.
      integer, intent(in), optional :: unit !! Output unit; defaults to standard output.
      logical, intent(in), optional :: print_obs !! Whether to print fitted values and residuals.
      integer :: destination, index
      logical :: show_observations

      destination = output_unit
      if (present(unit)) destination = unit
      show_observations = .false.
      if (present(print_obs)) show_observations = print_obs
      write(destination, '(a)') 'Exact common-MA threshold-regression fit'
      write(destination, '(a, i0)') 'Status: ', fit%info
      write(destination, '(a, es14.6)') 'Threshold: ', fit%threshold
      write(destination, '(a, i0)') 'Threshold delay: ', fit%delay
      write(destination, '(a, es14.6)') 'Residual sum of squares: ', fit%rss
      write(destination, '(a, es14.6)') 'AIC: ', fit%aic
      write(destination, '(a, es14.6)') 'BIC: ', fit%bic
      call display_real_vector(destination, 'Coefficients:', &
         fit%noise_fit%coefficients)
      if (show_observations .and. allocated(fit%fitted) .and. &
         allocated(fit%residuals)) then
         write(destination, '(a)') 'Index, fitted value, residual:'
         do index = 1, size(fit%fitted)
            write(destination, '(i8, 2(1x, es14.6))') index, &
               fit%fitted(index), fit%residuals(index)
         end do
      end if
   end subroutine display_tseriestarma_fit2

   subroutine display_tseriestarma_tar_test(test, unit)
      !! Display AR-versus-TAR supremum LM test results.
      type(tseriestarma_tar_test_t), intent(in) :: test !! TAR test result to display.
      integer, intent(in), optional :: unit !! Output unit; defaults to standard output.
      integer :: destination

      destination = output_unit
      if (present(unit)) destination = unit
      write(destination, '(a)') 'AR versus TAR supremum LM test'
      call display_tar_test_summary(test, destination)
   end subroutine display_tseriestarma_tar_test

   subroutine display_tseriestarma_tarma_test(test, unit)
      !! Display ARMA-versus-TARMA supremum LM test results.
      type(tseriestarma_tarma_test_t), intent(in) :: test !! TARMA test result to display.
      integer, intent(in), optional :: unit !! Output unit; defaults to standard output.
      integer :: destination

      destination = output_unit
      if (present(unit)) destination = unit
      write(destination, '(a)') 'ARMA versus TARMA supremum LM test'
      write(destination, '(a, i0)') 'MA order: ', test%ma_order
      write(destination, '(a, l1)') 'MA coefficients fixed under alternative: ', &
         test%ma_fixed
      call display_tar_test_summary(test%tseriestarma_tar_test_t, destination)
   end subroutine display_tseriestarma_tarma_test

   subroutine display_tseriestarma_tar_bootstrap(test, unit)
      !! Display bootstrap AR-versus-TAR test results.
      type(tseriestarma_tar_bootstrap_t), intent(in) :: test !! Bootstrap result to display.
      integer, intent(in), optional :: unit !! Output unit; defaults to standard output.
      integer :: destination

      destination = output_unit
      if (present(unit)) destination = unit
      write(destination, '(a)') 'Bootstrap AR versus TAR test'
      call display_tar_test_summary(test%tseriestarma_tar_test_t, destination)
      write(destination, '(a, i0)') 'Bootstrap method: ', test%bootstrap_type
      write(destination, '(a, i0)') 'Replications: ', test%replications
      write(destination, '(a, 2(es14.6, 1x))') 'Component p-values: ', test%p_values
      write(destination, '(a, es14.6)') 'Overall p-value: ', test%p_value
   end subroutine display_tseriestarma_tar_bootstrap

   subroutine display_tseriestarma_unit_root_test(test, unit)
      !! Display IMA-versus-TARMA unit-root test results.
      type(tseriestarma_unit_root_test_t), intent(in) :: test !! Unit-root test result to display.
      integer, intent(in), optional :: unit !! Output unit; defaults to standard output.
      integer :: destination

      destination = output_unit
      if (present(unit)) destination = unit
      write(destination, '(a)') 'IMA(1) versus TARMA(1,1) unit-root test'
      call display_unit_root_summary(test, destination)
   end subroutine display_tseriestarma_unit_root_test

   subroutine display_tseriestarma_unit_root_bootstrap(test, unit)
      !! Display bootstrap IMA-versus-TARMA unit-root test results.
      type(tseriestarma_unit_root_bootstrap_t), intent(in) :: test !! Bootstrap result to display.
      integer, intent(in), optional :: unit !! Output unit; defaults to standard output.
      integer :: destination

      destination = output_unit
      if (present(unit)) destination = unit
      write(destination, '(a)') 'Bootstrap IMA versus TARMA unit-root test'
      call display_unit_root_summary(test%tseriestarma_unit_root_test_t, &
         destination)
      write(destination, '(a, i0)') 'Bootstrap method: ', test%bootstrap_type
      write(destination, '(a, i0)') 'Replications: ', test%replications
      write(destination, '(a, es14.6)') 'Bootstrap p-value: ', test%p_value
   end subroutine display_tseriestarma_unit_root_bootstrap

   subroutine display_tseriestarma_arma_garch_fit(fit, unit, print_obs)
      !! Display a joint Gaussian ARMA-GARCH fit.
      type(tseriestarma_arma_garch_fit_t), intent(in) :: fit !! ARMA-GARCH fit to display.
      integer, intent(in), optional :: unit !! Output unit; defaults to standard output.
      logical, intent(in), optional :: print_obs !! Whether to print residuals and variances.
      integer :: destination, index
      logical :: show_observations

      destination = output_unit
      if (present(unit)) destination = unit
      show_observations = .false.
      if (present(print_obs)) show_observations = print_obs
      write(destination, '(a)') 'Gaussian ARMA-GARCH fit'
      write(destination, '(a, i0)') 'Status: ', fit%info
      write(destination, '(a, l1)') 'Converged: ', fit%converged
      write(destination, '(a, es14.6)') 'Mean: ', fit%mean
      call display_real_vector(destination, 'AR coefficients:', fit%ar)
      call display_real_vector(destination, 'MA coefficients:', fit%ma)
      write(destination, '(a, es14.6)') 'Variance intercept: ', &
         fit%variance_intercept
      call display_real_vector(destination, 'ARCH coefficients:', fit%arch)
      call display_real_vector(destination, 'GARCH coefficients:', fit%garch)
      write(destination, '(a, es14.6)') 'Log likelihood: ', fit%log_likelihood
      if (show_observations .and. allocated(fit%residuals) .and. &
         allocated(fit%conditional_variance)) then
         write(destination, '(a)') 'Index, residual, conditional variance:'
         do index = 1, size(fit%residuals)
            write(destination, '(i8, 2(1x, es14.6))') index, &
               fit%residuals(index), fit%conditional_variance(index)
         end do
      end if
   end subroutine display_tseriestarma_arma_garch_fit

   subroutine display_tseriestarma_garch_test(test, unit)
      !! Display GARCH-aware ARMA-versus-TARMA test results.
      type(tseriestarma_garch_test_t), intent(in) :: test !! GARCH test result to display.
      integer, intent(in), optional :: unit !! Output unit; defaults to standard output.
      integer :: destination

      destination = output_unit
      if (present(unit)) destination = unit
      write(destination, '(a)') 'GARCH-aware ARMA versus TARMA test'
      write(destination, '(a, i0)') 'Status: ', test%info
      write(destination, '(a, es14.6)') 'Supremum LM statistic: ', test%statistic
      write(destination, '(a, es14.6)') 'Maximizing threshold: ', &
         test%maximizing_threshold
      write(destination, '(a, es14.6)') 'Lower-regime proportion: ', &
         test%lower_regime_proportion
      write(destination, '(a, i0)') 'Degrees of freedom: ', test%degrees_of_freedom
      write(destination, '(a, 3(es14.6, 1x))') 'Critical probabilities: ', &
         test%critical_probabilities
      write(destination, '(a, 3(es14.6, 1x))') 'Critical values: ', &
         test%critical_values
   end subroutine display_tseriestarma_garch_test

   pure function tseriestarma_model_coefficients(model) result(values)
      !! Return packed lower and upper AR and MA model coefficients.
      type(tseriestarma_model_t), intent(in) :: model !! TARMA model.
      real(dp), allocatable :: values(:)

      if (allocated(model%phi_lower) .and. allocated(model%phi_upper) .and. &
         allocated(model%theta_lower) .and. allocated(model%theta_upper)) then
         values = [model%phi_lower, model%phi_upper, model%theta_lower, &
            model%theta_upper]
      else
         allocate(values(0))
      end if
   end function tseriestarma_model_coefficients

   pure function tseriestarma_fit_coefficients(fit) result(values)
      !! Return the packed coefficients of a conditional TARMA fit.
      type(tseriestarma_fit_t), intent(in) :: fit !! Conditional TARMA fit.
      real(dp), allocatable :: values(:)

      if (allocated(fit%parameters)) then
         values = fit%parameters
      else
         allocate(values(0))
      end if
   end function tseriestarma_fit_coefficients

   pure function tseriestarma_fit2_coefficients(fit) result(values)
      !! Return the coefficients of an exact common-MA TARMA fit.
      type(tseriestarma_fit2_t), intent(in) :: fit !! Exact TARMA fit.
      real(dp), allocatable :: values(:)

      if (allocated(fit%noise_fit%coefficients)) then
         values = fit%noise_fit%coefficients
      else
         allocate(values(0))
      end if
   end function tseriestarma_fit2_coefficients

   pure function tseriestarma_arma_garch_coefficients(fit) result(values)
      !! Return packed mean, AR, MA, variance, ARCH, and GARCH coefficients.
      type(tseriestarma_arma_garch_fit_t), intent(in) :: fit !! ARMA-GARCH fit.
      real(dp), allocatable :: values(:)

      if (allocated(fit%ar) .and. allocated(fit%ma) .and. &
         allocated(fit%arch) .and. allocated(fit%garch)) then
         values = [fit%mean, fit%ar, fit%ma, fit%variance_intercept, &
            fit%arch, fit%garch]
      else
         allocate(values(0))
      end if
   end function tseriestarma_arma_garch_coefficients

   pure function tseriestarma_fit_covariance(fit) result(values)
      !! Return the coefficient covariance matrix of a conditional TARMA fit.
      type(tseriestarma_fit_t), intent(in) :: fit !! Conditional TARMA fit.
      real(dp), allocatable :: values(:, :)

      if (allocated(fit%covariance)) then
         values = fit%covariance
      else
         allocate(values(0, 0))
      end if
   end function tseriestarma_fit_covariance

   pure function tseriestarma_fit2_covariance(fit) result(values)
      !! Return the coefficient covariance matrix of an exact TARMA fit.
      type(tseriestarma_fit2_t), intent(in) :: fit !! Exact TARMA fit.
      real(dp), allocatable :: values(:, :)

      if (allocated(fit%noise_fit%covariance)) then
         values = fit%noise_fit%covariance
      else
         allocate(values(0, 0))
      end if
   end function tseriestarma_fit2_covariance

   pure function tseriestarma_fit_residuals(fit) result(values)
      !! Return conditional innovations from a TARMA fit.
      type(tseriestarma_fit_t), intent(in) :: fit !! Conditional TARMA fit.
      real(dp), allocatable :: values(:)

      if (allocated(fit%evaluation%residuals)) then
         values = fit%evaluation%residuals
      else
         allocate(values(0))
      end if
   end function tseriestarma_fit_residuals

   pure function tseriestarma_fit2_residuals(fit) result(values)
      !! Return innovations from an exact common-MA TARMA fit.
      type(tseriestarma_fit2_t), intent(in) :: fit !! Exact TARMA fit.
      real(dp), allocatable :: values(:)

      if (allocated(fit%residuals)) then
         values = fit%residuals
      else
         allocate(values(0))
      end if
   end function tseriestarma_fit2_residuals

   pure function tseriestarma_arma_garch_residuals(fit) result(values)
      !! Return innovations from a joint ARMA-GARCH fit.
      type(tseriestarma_arma_garch_fit_t), intent(in) :: fit !! ARMA-GARCH fit.
      real(dp), allocatable :: values(:)

      if (allocated(fit%residuals)) then
         values = fit%residuals
      else
         allocate(values(0))
      end if
   end function tseriestarma_arma_garch_residuals

   pure function tseriestarma_fit_fitted_values(fit) result(values)
      !! Return conditional fitted values from a TARMA fit.
      type(tseriestarma_fit_t), intent(in) :: fit !! Conditional TARMA fit.
      real(dp), allocatable :: values(:)

      if (allocated(fit%evaluation%fitted)) then
         values = fit%evaluation%fitted
      else
         allocate(values(0))
      end if
   end function tseriestarma_fit_fitted_values

   pure function tseriestarma_fit2_fitted_values(fit) result(values)
      !! Return fitted values from an exact common-MA TARMA fit.
      type(tseriestarma_fit2_t), intent(in) :: fit !! Exact TARMA fit.
      real(dp), allocatable :: values(:)

      if (allocated(fit%fitted)) then
         values = fit%fitted
      else
         allocate(values(0))
      end if
   end function tseriestarma_fit2_fitted_values

   subroutine display_real_vector(unit, heading, values)
      !! Write an allocated real vector with a heading.
      integer, intent(in) :: unit !! Output unit.
      character(len=*), intent(in) :: heading !! Heading written before the values.
      real(dp), allocatable, intent(in) :: values(:) !! Vector to write when allocated.

      write(unit, '(a)') heading
      if (allocated(values)) then
         if (size(values) > 0) write(unit, '(*(es14.6, 1x))') values
      else
         write(unit, '(a)') 'not allocated'
      end if
   end subroutine display_real_vector

   subroutine display_tar_test_summary(test, unit)
      !! Write fields shared by TAR and TARMA supremum LM tests.
      type(tseriestarma_tar_test_t), intent(in) :: test !! Test result to display.
      integer, intent(in) :: unit !! Output unit.

      write(unit, '(a, i0)') 'Status: ', test%info
      write(unit, '(a, i0)') 'AR order: ', test%ar_order
      write(unit, '(a, i0)') 'Threshold delay: ', test%delay
      write(unit, '(a, 2(es14.6, 1x))') 'LM statistics: ', test%statistic
      write(unit, '(a, 2(es14.6, 1x))') 'Maximizing thresholds: ', &
         test%maximizing_threshold
      write(unit, '(a, 2(es14.6, 1x))') 'Lower-regime proportions: ', &
         test%lower_regime_proportion
   end subroutine display_tar_test_summary

   subroutine display_unit_root_summary(test, unit)
      !! Write fields shared by asymptotic and bootstrap unit-root tests.
      type(tseriestarma_unit_root_test_t), intent(in) :: test !! Test result to display.
      integer, intent(in) :: unit !! Output unit.

      write(unit, '(a, i0)') 'Status: ', test%info
      write(unit, '(a, es14.6)') 'Supremum LM statistic: ', test%statistic
      write(unit, '(a, es14.6)') 'Maximizing threshold: ', &
         test%maximizing_threshold
      write(unit, '(a, es14.6)') 'Null MA coefficient: ', test%ma_coefficient
      write(unit, '(a, 4(es14.6, 1x))') 'Critical probabilities: ', &
         test%critical_probabilities
      write(unit, '(a, 4(es14.6, 1x))') 'Critical values: ', &
         test%critical_values
   end subroutine display_unit_root_summary

   pure function tseriestarma_tar_test(series, ar_order, delay, &
      lower_probability, upper_probability) result(out)
      !! Compute classic and heteroskedasticity-robust AR-versus-TAR supLM tests.
      real(dp), intent(in) :: series(:) !! Univariate time-series observations.
      integer, intent(in) :: ar_order !! Positive order of the null autoregression.
      integer, intent(in), optional :: delay !! Positive delay of the threshold variable.
      real(dp), intent(in), optional :: lower_probability !! Lower threshold-search quantile.
      real(dp), intent(in), optional :: upper_probability !! Upper threshold-search quantile.
      type(tseriestarma_tar_test_t) :: out
      real(dp), allocatable :: design(:, :), response(:), threshold_variable(:)
      real(dp), allocatable :: ordered_thresholds(:), cross(:, :), cross_inverse(:, :)
      real(dp), allocatable :: dphi(:, :), dpsi(:, :), weighted_dpsi(:, :)
      real(dp), allocatable :: null_score_map(:, :), score_terms(:, :)
      real(dp), allocatable :: m21(:, :), information(:, :), information_inverse(:, :)
      real(dp), allocatable :: robust_covariance(:, :), robust_inverse(:, :)
      real(dp), allocatable :: score(:), robust_score(:), indicator(:)
      real(dp) :: rss, threshold
      integer :: first, last, grid_index, info, k, n, neff, q, row, time
      integer :: classic_location, robust_location

      out%ar_order = ar_order
      if (present(delay)) out%delay = delay
      if (present(lower_probability)) out%lower_probability = lower_probability
      if (present(upper_probability)) out%upper_probability = upper_probability
      n = size(series)
      k = max(ar_order, out%delay)
      neff = n - k
      q = ar_order + 1
      out%effective_observations = max(0, neff)
      out%degrees_of_freedom = q
      if (ar_order < 1 .or. out%delay < 1 .or. neff <= q) then
         out%info = 1
         return
      end if
      if (out%lower_probability <= 0.0_dp .or. &
         out%upper_probability >= 1.0_dp .or. &
         out%lower_probability >= out%upper_probability) then
         out%info = 2
         return
      end if
      first = ceiling(real(neff, dp)*out%lower_probability)
      last = floor(real(neff, dp)*out%upper_probability)
      first = max(1, first)
      last = min(neff, last)
      if (first > last) then
         out%info = 2
         return
      end if

      allocate(design(neff, q), response(neff), threshold_variable(neff))
      design(:, 1) = 1.0_dp
      do row = 1, neff
         time = k + row
         response(row) = series(time)
         threshold_variable(row) = series(time - out%delay)
         do grid_index = 1, ar_order
            design(row, grid_index + 1) = series(time - grid_index)
         end do
      end do
      allocate(cross(q, q), cross_inverse(q, q))
      cross = matmul(transpose(design), design)
      call symmetric_pseudoinverse(cross, cross_inverse, info)
      if (info /= 0) then
         out%info = 3
         return
      end if
      allocate(out%coefficients(q), out%residuals(neff))
      out%coefficients = matmul(cross_inverse, matmul(transpose(design), response))
      out%residuals = response - matmul(design, out%coefficients)
      rss = sum(out%residuals**2)
      out%sigma2 = rss/real(neff, dp)
      if (out%sigma2 <= tiny(1.0_dp)) then
         out%info = 3
         return
      end if

      ordered_thresholds = sorted(threshold_variable)
      out%threshold_values = ordered_thresholds(first:last)
      allocate(out%test_values(size(out%threshold_values), 2))
      allocate(dphi(neff, q), dpsi(neff, q), weighted_dpsi(neff, q))
      allocate(null_score_map(q, neff), score_terms(q, neff))
      allocate(m21(q, q), information(q, q), information_inverse(q, q))
      allocate(robust_covariance(q, q), robust_inverse(q, q))
      allocate(score(q), robust_score(q), indicator(neff))
      dphi = -design
      null_score_map = matmul(cross_inverse, transpose(dphi* &
         spread(out%residuals, 2, q)))
      do grid_index = 1, size(out%threshold_values)
         threshold = out%threshold_values(grid_index)
         indicator = merge(1.0_dp, 0.0_dp, threshold_variable <= threshold)
         dpsi = -design*spread(indicator, 2, q)
         weighted_dpsi = dpsi*spread(out%residuals, 2, q)
         m21 = matmul(transpose(dpsi), dphi)
         score = -sum(weighted_dpsi, dim=1)
         score_terms = -transpose(weighted_dpsi) + matmul(m21, null_score_map)
         robust_score = sum(score_terms, dim=2)

         information = m21 - matmul(m21, matmul(cross_inverse, m21))
         robust_covariance = matmul(score_terms, transpose(score_terms))
         call symmetric_pseudoinverse(information, information_inverse, info)
         if (info /= 0) then
            out%info = 4
            return
         end if
         call symmetric_pseudoinverse(robust_covariance, robust_inverse, info)
         if (info /= 0) then
            out%info = 4
            return
         end if
         out%test_values(grid_index, 1) = max(0.0_dp, &
            dot_product(score, matmul(information_inverse, score))/out%sigma2)
         out%test_values(grid_index, 2) = max(0.0_dp, &
            dot_product(robust_score, matmul(robust_inverse, robust_score)))
      end do

      out%statistic = maxval(out%test_values, dim=1)
      classic_location = maxloc(out%test_values(:, 1), dim=1)
      robust_location = maxloc(out%test_values(:, 2), dim=1)
      out%maximizing_threshold = [out%threshold_values(classic_location), &
         out%threshold_values(robust_location)]
      out%lower_regime_proportion(1) = real(count(threshold_variable <= &
         out%maximizing_threshold(1)), dp)/real(neff, dp)
      out%lower_regime_proportion(2) = real(count(threshold_variable <= &
         out%maximizing_threshold(2)), dp)/real(neff, dp)
   end function tseriestarma_tar_test

   pure function tseriestarma_tar_bootstrap_from_random(series, ar_order, &
      bootstrap_type, random_values, delay, lower_probability, &
      upper_probability, burnin) result(out)
      !! Bootstrap an AR-versus-TAR test using supplied resampling variates.
      real(dp), intent(in) :: series(:) !! Univariate time-series observations.
      integer, intent(in) :: ar_order !! Positive order of the null autoregression.
      integer, intent(in) :: bootstrap_type !! Bootstrap method constant.
      real(dp), intent(in) :: random_values(:, :) !! Uniforms or bootstrap multipliers by replication.
      integer, intent(in), optional :: delay !! Positive delay of the threshold variable.
      real(dp), intent(in), optional :: lower_probability !! Lower threshold-search quantile.
      real(dp), intent(in), optional :: upper_probability !! Upper threshold-search quantile.
      integer, intent(in), optional :: burnin !! Transient length for simulated-null bootstraps.
      type(tseriestarma_tar_bootstrap_t) :: out
      type(tseriestarma_tar_test_t) :: observed, replicate_test
      real(dp), allocatable :: centered(:), bootstrap_residuals(:), bootstrap_series(:)
      real(dp), allocatable :: design(:, :), threshold_variable(:), ordered_thresholds(:)
      real(dp), allocatable :: thresholds(:), dphi(:, :), dpsi(:, :)
      real(dp), allocatable :: cross(:, :), cross_inverse(:, :), null_score_map(:, :)
      real(dp), allocatable :: weighted_dpsi(:, :), score_terms(:, :)
      real(dp), allocatable :: m21(:, :), information(:, :), information_inverse(:, :)
      real(dp), allocatable :: robust_covariance(:, :), robust_inverse(:, :)
      real(dp), allocatable :: draw_scores(:, :), indicator(:)
      real(dp) :: pa, pb, threshold
      integer :: selected_delay, selected_burnin, n, k, neff, np, q
      integer :: first, last, grid_index, replication, row, time, lag, index, info

      selected_delay = 1
      if (present(delay)) selected_delay = delay
      pa = 0.25_dp
      if (present(lower_probability)) pa = lower_probability
      pb = 0.75_dp
      if (present(upper_probability)) pb = upper_probability
      observed = tseriestarma_tar_test(series, ar_order, selected_delay, pa, pb)
      out%tseriestarma_tar_test_t = observed
      out%bootstrap_type = bootstrap_type
      out%replications = size(random_values, 2)
      if (observed%info /= 0) return
      if (bootstrap_type < tseriestarma_bootstrap_iid .or. &
         bootstrap_type > tseriestarma_bootstrap_normal .or. &
         size(random_values, 2) < 1 .or. &
         .not. all(ieee_is_finite(random_values))) then
         out%info = 5
         return
      end if
      n = size(series)
      k = max(ar_order, selected_delay)
      neff = n - k
      q = ar_order + 1
      selected_burnin = max(k, n/3)
      if (present(burnin)) selected_burnin = burnin
      out%burnin = merge(0, selected_burnin, &
         bootstrap_type == tseriestarma_bootstrap_hansen)
      if (bootstrap_type == tseriestarma_bootstrap_hansen) then
         if (size(random_values, 1) /= neff) then
            out%info = 5
            return
         end if
      else
         if (selected_burnin < k .or. selected_burnin > neff) then
            out%info = 5
            return
         end if
         np = neff + selected_burnin
         if (size(random_values, 1) /= np) then
            out%info = 5
            return
         end if
         if (bootstrap_type == tseriestarma_bootstrap_iid .and. &
            (any(random_values < 0.0_dp) .or. any(random_values >= 1.0_dp))) then
            out%info = 5
            return
         end if
         if (bootstrap_type == tseriestarma_bootstrap_rademacher .and. &
            any(abs(abs(random_values) - 1.0_dp) > 10.0_dp*epsilon(1.0_dp))) then
            out%info = 5
            return
         end if
      end if

      allocate(out%bootstrap_statistics(out%replications, 2), source=0.0_dp)
      if (bootstrap_type /= tseriestarma_bootstrap_hansen) then
         centered = observed%residuals - &
            sum(observed%residuals)/real(neff, dp)
         allocate(bootstrap_residuals(np), bootstrap_series(np))
         do replication = 1, out%replications
            if (bootstrap_type == tseriestarma_bootstrap_iid) then
               do time = 1, np
                  index = min(neff, 1 + int(random_values(time, replication)* &
                     real(neff, dp)))
                  bootstrap_residuals(time) = centered(index)
               end do
            else
               do time = 1, np
                  if (time <= selected_burnin) then
                     bootstrap_residuals(time) = centered(time)* &
                        random_values(time, replication)
                  else
                     bootstrap_residuals(time) = centered(time - selected_burnin)* &
                        random_values(time, replication)
                  end if
               end do
            end if
            bootstrap_series = bootstrap_residuals
            do time = ar_order + 1, np
               bootstrap_series(time) = observed%coefficients(1) + &
                  bootstrap_residuals(time)
               do lag = 1, ar_order
                  bootstrap_series(time) = bootstrap_series(time) + &
                     observed%coefficients(lag + 1)*bootstrap_series(time - lag)
               end do
            end do
            replicate_test = tseriestarma_tar_test(bootstrap_series( &
               selected_burnin - k + 1:), ar_order, selected_delay, pa, pb)
            if (replicate_test%info /= 0) then
               out%info = 6
               return
            end if
            out%bootstrap_statistics(replication, :) = replicate_test%statistic
         end do
      else
         first = max(1, ceiling(real(neff - 1, dp)*pa))
         last = min(neff, floor(real(neff - 1, dp)*pb))
         if (first > last) then
            out%info = 5
            return
         end if
         allocate(design(neff, q), threshold_variable(neff))
         design(:, 1) = 1.0_dp
         do row = 1, neff
            time = k + row
            threshold_variable(row) = series(time - selected_delay)
            do lag = 1, ar_order
               design(row, lag + 1) = series(time - lag)
            end do
         end do
         ordered_thresholds = sorted(threshold_variable)
         thresholds = ordered_thresholds(first:last)
         allocate(dphi(neff, q), cross(q, q), cross_inverse(q, q))
         allocate(null_score_map(q, neff), dpsi(neff, q))
         allocate(weighted_dpsi(neff, q), score_terms(q, neff))
         allocate(m21(q, q), information(q, q), information_inverse(q, q))
         allocate(robust_covariance(q, q), robust_inverse(q, q))
         allocate(draw_scores(q, out%replications), indicator(neff))
         dphi = -design
         cross = matmul(transpose(dphi), dphi)
         call symmetric_pseudoinverse(cross, cross_inverse, info)
         if (info /= 0) then
            out%info = 4
            return
         end if
         null_score_map = matmul(cross_inverse, transpose(dphi* &
            spread(observed%residuals, 2, q)))
         do grid_index = 1, size(thresholds)
            threshold = thresholds(grid_index)
            indicator = merge(1.0_dp, 0.0_dp, threshold_variable <= threshold)
            dpsi = -design*spread(indicator, 2, q)
            weighted_dpsi = dpsi*spread(observed%residuals, 2, q)
            m21 = matmul(transpose(dpsi), dphi)
            information = m21 - matmul(m21, matmul(cross_inverse, m21))
            score_terms = -transpose(weighted_dpsi) + &
               matmul(m21, null_score_map)
            robust_covariance = matmul(score_terms, transpose(score_terms))
            call symmetric_pseudoinverse(information, information_inverse, info)
            if (info /= 0) then
               out%info = 4
               return
            end if
            call symmetric_pseudoinverse(robust_covariance, robust_inverse, info)
            if (info /= 0) then
               out%info = 4
               return
            end if
            draw_scores = matmul(score_terms, random_values)
            do replication = 1, out%replications
               out%bootstrap_statistics(replication, 1) = max( &
                  out%bootstrap_statistics(replication, 1), &
                  dot_product(draw_scores(:, replication), &
                  matmul(information_inverse, draw_scores(:, replication)))/ &
                  observed%sigma2)
               out%bootstrap_statistics(replication, 2) = max( &
                  out%bootstrap_statistics(replication, 2), &
                  dot_product(draw_scores(:, replication), &
                  matmul(robust_inverse, draw_scores(:, replication))))
            end do
         end do
      end if
      do index = 1, 2
         out%p_values(index) = real(count(out%bootstrap_statistics(:, index) > &
            observed%statistic(index)), dp)/real(out%replications, dp)
      end do
      out%p_value = out%p_values(1)
   end function tseriestarma_tar_bootstrap_from_random

   function tseriestarma_tar_bootstrap(series, ar_order, replications, &
      bootstrap_type, delay, lower_probability, upper_probability, &
      burnin) result(out)
      !! Bootstrap an AR-versus-TAR test using the shared random stream.
      real(dp), intent(in) :: series(:) !! Univariate time-series observations.
      integer, intent(in) :: ar_order !! Positive order of the null autoregression.
      integer, intent(in), optional :: replications !! Number of bootstrap replications.
      integer, intent(in), optional :: bootstrap_type !! Bootstrap method constant.
      integer, intent(in), optional :: delay !! Positive delay of the threshold variable.
      real(dp), intent(in), optional :: lower_probability !! Lower threshold-search quantile.
      real(dp), intent(in), optional :: upper_probability !! Upper threshold-search quantile.
      integer, intent(in), optional :: burnin !! Transient length for simulated-null bootstraps.
      type(tseriestarma_tar_bootstrap_t) :: out
      real(dp), allocatable :: random_values(:, :)
      real(dp) :: pa, pb
      integer :: selected_replications, selected_type, selected_delay
      integer :: selected_burnin, k, neff, rows, row, replication

      selected_replications = 1000
      if (present(replications)) selected_replications = replications
      selected_type = tseriestarma_bootstrap_iid
      if (present(bootstrap_type)) selected_type = bootstrap_type
      selected_delay = 1
      if (present(delay)) selected_delay = delay
      pa = 0.25_dp
      if (present(lower_probability)) pa = lower_probability
      pb = 0.75_dp
      if (present(upper_probability)) pb = upper_probability
      k = max(ar_order, selected_delay)
      neff = size(series) - k
      selected_burnin = max(k, size(series)/3)
      if (present(burnin)) selected_burnin = burnin
      rows = neff + selected_burnin
      if (selected_type == tseriestarma_bootstrap_hansen) rows = neff
      if (selected_replications < 1 .or. rows < 1) then
         out%info = 5
         return
      end if
      allocate(random_values(rows, selected_replications))
      select case (selected_type)
      case (tseriestarma_bootstrap_iid)
         do replication = 1, selected_replications
            do row = 1, rows
               random_values(row, replication) = random_uniform()
            end do
         end do
      case (tseriestarma_bootstrap_rademacher)
         do replication = 1, selected_replications
            do row = 1, rows
               random_values(row, replication) = &
                  merge(-1.0_dp, 1.0_dp, random_uniform() < 0.5_dp)
            end do
         end do
      case (tseriestarma_bootstrap_hansen, tseriestarma_bootstrap_normal)
         call random_standard_normal_matrix(random_values)
      case default
         out%info = 5
         return
      end select
      out = tseriestarma_tar_bootstrap_from_random(series, ar_order, &
         selected_type, random_values, selected_delay, pa, pb, &
         selected_burnin)
   end function tseriestarma_tar_bootstrap

   pure function tseriestarma_unit_root_test(series, threshold_range, &
      lower_probability, upper_probability, initial, max_iterations, &
      tolerance) result(out)
      !! Test an integrated MA(1) null against a stationary TARMA(1,1).
      real(dp), intent(in) :: series(:) !! Univariate time-series observations.
      real(dp), intent(in), optional :: threshold_range(:) !! Threshold values to test.
      real(dp), intent(in), optional :: lower_probability !! Lower threshold-search quantile.
      real(dp), intent(in), optional :: upper_probability !! Upper threshold-search quantile.
      real(dp), intent(in), optional :: initial(:) !! Initial MA and drift estimates.
      integer, intent(in), optional :: max_iterations !! Maximum null-fit iterations.
      real(dp), intent(in), optional :: tolerance !! Null-fit convergence tolerance.
      type(tseriestarma_unit_root_test_t) :: out
      real(dp), allocatable :: starting(:), ordered(:), derivative_input(:, :)
      real(dp), allocatable :: nuisance_derivative(:, :), alternative_input(:, :)
      real(dp), allocatable :: alternative_derivative(:, :), residuals(:)
      real(dp), allocatable :: indicator(:), score(:), scale(:)
      real(dp), allocatable :: alternative_information(:, :), cross_information(:, :)
      real(dp), allocatable :: efficient_information(:, :), efficient_inverse(:, :)
      real(dp) :: selected_tolerance, nuisance_information, threshold
      integer :: first, last, grid_index, info, iteration_limit, n, neff
      integer :: location

      n = size(series)
      neff = n - 1
      out%effective_observations = max(0, neff)
      if (n < 6 .or. .not. all(ieee_is_finite(series))) then
         out%info = 1
         return
      end if
      if (present(lower_probability)) out%lower_probability = lower_probability
      if (present(upper_probability)) out%upper_probability = upper_probability
      if (out%lower_probability <= 0.0_dp .or. &
         out%upper_probability >= 1.0_dp .or. &
         out%lower_probability >= out%upper_probability) then
         out%info = 2
         return
      end if
      allocate(starting(2), source=0.0_dp)
      if (present(initial)) then
         if (size(initial) /= 2 .or. .not. all(ieee_is_finite(initial))) then
            out%info = 1
            return
         end if
         starting = initial
      end if
      iteration_limit = 200
      if (present(max_iterations)) iteration_limit = max_iterations
      selected_tolerance = 1.0e-6_dp
      if (present(tolerance)) selected_tolerance = tolerance
      if (iteration_limit < 1 .or. selected_tolerance <= 0.0_dp) then
         out%info = 1
         return
      end if
      out%noise_fit = sarima_fit(series, 0, 1, 1, 0, 0, 0, 1, &
         initial=starting, include_intercept=.false., include_drift=.true., &
         max_iterations=iteration_limit, tolerance=selected_tolerance, &
         transform_parameters=.true., exact_likelihood=.true.)
      if (out%noise_fit%likelihood%info /= 0 .or. &
         .not. allocated(out%noise_fit%coefficients) .or. &
         size(out%noise_fit%coefficients) /= 2 .or. &
         .not. allocated(out%noise_fit%likelihood%residuals)) then
         out%info = 3
         return
      end if
      out%ma_coefficient = out%noise_fit%coefficients(1)
      out%sigma2 = out%noise_fit%likelihood%sigma2
      if (out%sigma2 <= tiny(1.0_dp)) then
         out%info = 3
         return
      end if
      residuals = out%noise_fit%likelihood%residuals(2:n)
      out%residuals = residuals
      if (present(threshold_range)) then
         if (size(threshold_range) < 1 .or. &
            .not. all(ieee_is_finite(threshold_range))) then
            out%info = 2
            return
         end if
         out%threshold_values = threshold_range
      else
         ordered = sorted(series)
         first = max(1, ceiling(real(neff, dp)*out%lower_probability))
         last = min(n, floor(real(neff, dp)*out%upper_probability))
         if (first > last) then
            out%info = 2
            return
         end if
         out%threshold_values = ordered(first:last)
      end if

      allocate(derivative_input(neff, 1), alternative_input(neff, 2))
      allocate(indicator(neff), score(2), scale(2))
      derivative_input(:, 1) = -1.0_dp
      nuisance_derivative = recursive_filter_matrix(derivative_input, &
         [-out%ma_coefficient])
      nuisance_information = sum(nuisance_derivative(:, 1)**2)/ &
         (out%sigma2*real(neff, dp))
      if (nuisance_information <= tiny(1.0_dp)) then
         out%info = 4
         return
      end if
      scale = [sqrt(real(neff, dp)), real(neff, dp)]
      allocate(out%test_values(size(out%threshold_values)))
      allocate(alternative_information(2, 2), cross_information(2, 1))
      allocate(efficient_information(2, 2), efficient_inverse(2, 2))
      do grid_index = 1, size(out%threshold_values)
         threshold = out%threshold_values(grid_index)
         indicator = merge(1.0_dp, 0.0_dp, series(:neff) <= threshold)
         alternative_input(:, 1) = -indicator
         alternative_input(:, 2) = -series(:neff)*indicator
         alternative_derivative = recursive_filter_matrix(alternative_input, &
            [-out%ma_coefficient])
         score = -sum(alternative_derivative*spread(residuals, 2, 2), dim=1)/ &
            (out%sigma2*scale)
         alternative_information = matmul(transpose(alternative_derivative), &
            alternative_derivative)/out%sigma2
         alternative_information = alternative_information/ &
            (spread(scale, 2, 2)*spread(scale, 1, 2))
         cross_information = matmul(transpose(alternative_derivative), &
            nuisance_derivative)/out%sigma2
         cross_information(:, 1) = cross_information(:, 1)/ &
            (scale*sqrt(real(neff, dp)))
         efficient_information = alternative_information - &
            matmul(cross_information, transpose(cross_information))/ &
            nuisance_information
         call symmetric_pseudoinverse(efficient_information, &
            efficient_inverse, info)
         if (info /= 0) then
            out%info = 4
            return
         end if
         out%test_values(grid_index) = max(0.0_dp, &
            dot_product(score, matmul(efficient_inverse, score)))
      end do
      out%statistic = maxval(out%test_values)
      location = maxloc(out%test_values, dim=1)
      out%maximizing_threshold = out%threshold_values(location)
      out%critical_probabilities = &
         tseriestarma_unit_root_critical_probabilities()
      out%critical_values = tseriestarma_unit_root_critical_values(n, &
         out%lower_probability, out%ma_coefficient)
   end function tseriestarma_unit_root_test

   pure function tseriestarma_unit_root_bootstrap_from_random(series, &
      bootstrap_type, random_values, threshold_range, lower_probability, &
      upper_probability, max_iterations, tolerance) result(out)
      !! Bootstrap the IMA-versus-TARMA test using supplied random variates.
      real(dp), intent(in) :: series(:) !! Univariate time-series observations.
      integer, intent(in) :: bootstrap_type !! IID or wild-bootstrap method constant.
      real(dp), intent(in) :: random_values(:, :) !! IID uniforms, or an initial normal seed followed by wild multipliers.
      real(dp), intent(in), optional :: threshold_range(:) !! Threshold values to test.
      real(dp), intent(in), optional :: lower_probability !! Lower threshold-search quantile.
      real(dp), intent(in), optional :: upper_probability !! Upper threshold-search quantile.
      integer, intent(in), optional :: max_iterations !! Maximum null-fit iterations.
      real(dp), intent(in), optional :: tolerance !! Null-fit convergence tolerance.
      type(tseriestarma_unit_root_bootstrap_t) :: out
      type(tseriestarma_unit_root_test_t) :: observed, replicate_test
      real(dp), allocatable :: centered(:), bootstrap_residuals(:)
      real(dp), allocatable :: bootstrap_series(:), starting(:)
      real(dp) :: pa, pb, level
      integer :: n, np, replication, time, index, iteration_limit
      real(dp) :: selected_tolerance

      pa = 0.25_dp
      if (present(lower_probability)) pa = lower_probability
      pb = 0.75_dp
      if (present(upper_probability)) pb = upper_probability
      iteration_limit = 200
      if (present(max_iterations)) iteration_limit = max_iterations
      selected_tolerance = 1.0e-6_dp
      if (present(tolerance)) selected_tolerance = tolerance
      if (present(threshold_range)) then
         observed = tseriestarma_unit_root_test(series, threshold_range, pa, &
            pb, max_iterations=iteration_limit, tolerance=selected_tolerance)
      else
         observed = tseriestarma_unit_root_test(series, &
            lower_probability=pa, upper_probability=pb, &
            max_iterations=iteration_limit, tolerance=selected_tolerance)
      end if
      out%tseriestarma_unit_root_test_t = observed
      out%bootstrap_type = bootstrap_type
      out%replications = size(random_values, 2)
      if (observed%info /= 0) return
      n = size(series)
      np = n + 1
      if ((bootstrap_type /= tseriestarma_bootstrap_iid .and. &
         bootstrap_type /= tseriestarma_bootstrap_rademacher .and. &
         bootstrap_type /= tseriestarma_bootstrap_normal) .or. &
         size(random_values, 2) < 1 .or. &
         .not. all(ieee_is_finite(random_values))) then
         out%info = 5
         return
      end if
      if (bootstrap_type == tseriestarma_bootstrap_iid) then
         if (size(random_values, 1) /= np .or. any(random_values < 0.0_dp) .or. &
            any(random_values >= 1.0_dp)) then
            out%info = 5
            return
         end if
      else
         if (size(random_values, 1) /= np + 1) then
            out%info = 5
            return
         end if
         if (bootstrap_type == tseriestarma_bootstrap_rademacher .and. &
            any(abs(abs(random_values(2:, :)) - 1.0_dp) > &
            10.0_dp*epsilon(1.0_dp))) then
            out%info = 5
            return
         end if
      end if
      centered = observed%noise_fit%likelihood%residuals - &
         sum(observed%noise_fit%likelihood%residuals)/real(n, dp)
      allocate(out%bootstrap_statistics(out%replications))
      allocate(bootstrap_residuals(np), bootstrap_series(n))
      starting = observed%noise_fit%coefficients
      do replication = 1, out%replications
         if (bootstrap_type == tseriestarma_bootstrap_iid) then
            do time = 1, np
               index = min(n, 1 + int(random_values(time, replication)* &
                  real(n, dp)))
               bootstrap_residuals(time) = centered(index)
            end do
         else
            bootstrap_residuals(1) = random_values(1, replication)* &
               random_values(2, replication)
            do time = 2, np
               bootstrap_residuals(time) = centered(time - 1)* &
                  random_values(time + 1, replication)
            end do
         end if
         level = series(1)
         do time = 1, n
            level = level + bootstrap_residuals(time + 1) + &
               observed%ma_coefficient*bootstrap_residuals(time)
            bootstrap_series(time) = level
         end do
         if (present(threshold_range)) then
            replicate_test = tseriestarma_unit_root_test(bootstrap_series, &
               threshold_range, pa, pb, initial=starting, &
               max_iterations=iteration_limit, tolerance=selected_tolerance)
         else
            replicate_test = tseriestarma_unit_root_test(bootstrap_series, &
               lower_probability=pa, upper_probability=pb, initial=starting, &
               max_iterations=iteration_limit, tolerance=selected_tolerance)
         end if
         if (replicate_test%info /= 0) then
            out%info = 6
            return
         end if
         out%bootstrap_statistics(replication) = replicate_test%statistic
      end do
      out%p_value = real(count(out%bootstrap_statistics > observed%statistic), &
         dp)/real(out%replications, dp)
   end function tseriestarma_unit_root_bootstrap_from_random

   function tseriestarma_unit_root_bootstrap(series, replications, &
      bootstrap_type, threshold_range, lower_probability, upper_probability, &
      max_iterations, tolerance) result(out)
      !! Bootstrap the IMA-versus-TARMA test using the shared random stream.
      real(dp), intent(in) :: series(:) !! Univariate time-series observations.
      integer, intent(in), optional :: replications !! Number of bootstrap replications.
      integer, intent(in), optional :: bootstrap_type !! IID or wild-bootstrap method constant.
      real(dp), intent(in), optional :: threshold_range(:) !! Threshold values to test.
      real(dp), intent(in), optional :: lower_probability !! Lower threshold-search quantile.
      real(dp), intent(in), optional :: upper_probability !! Upper threshold-search quantile.
      integer, intent(in), optional :: max_iterations !! Maximum null-fit iterations.
      real(dp), intent(in), optional :: tolerance !! Null-fit convergence tolerance.
      type(tseriestarma_unit_root_bootstrap_t) :: out
      real(dp), allocatable :: random_values(:, :)
      real(dp) :: pa, pb
      integer :: selected_replications, selected_type, rows, row, replication
      integer :: iteration_limit
      real(dp) :: selected_tolerance

      selected_replications = 1000
      if (present(replications)) selected_replications = replications
      selected_type = tseriestarma_bootstrap_rademacher
      if (present(bootstrap_type)) selected_type = bootstrap_type
      pa = 0.25_dp
      if (present(lower_probability)) pa = lower_probability
      pb = 0.75_dp
      if (present(upper_probability)) pb = upper_probability
      iteration_limit = 200
      if (present(max_iterations)) iteration_limit = max_iterations
      selected_tolerance = 1.0e-6_dp
      if (present(tolerance)) selected_tolerance = tolerance
      rows = size(series) + 1
      if (selected_type /= tseriestarma_bootstrap_iid) rows = rows + 1
      if (selected_replications < 1 .or. rows < 1) then
         out%info = 5
         return
      end if
      allocate(random_values(rows, selected_replications))
      select case (selected_type)
      case (tseriestarma_bootstrap_iid)
         do replication = 1, selected_replications
            do row = 1, rows
               random_values(row, replication) = random_uniform()
            end do
         end do
      case (tseriestarma_bootstrap_rademacher)
         do replication = 1, selected_replications
            random_values(1, replication) = random_standard_normal()
            do row = 2, rows
               random_values(row, replication) = &
                  merge(-1.0_dp, 1.0_dp, random_uniform() < 0.5_dp)
            end do
         end do
      case (tseriestarma_bootstrap_normal)
         call random_standard_normal_matrix(random_values)
      case default
         out%info = 5
         return
      end select
      if (present(threshold_range)) then
         out = tseriestarma_unit_root_bootstrap_from_random(series, &
            selected_type, random_values, threshold_range, pa, pb, &
            iteration_limit, selected_tolerance)
      else
         out = tseriestarma_unit_root_bootstrap_from_random(series, &
            selected_type, random_values, lower_probability=pa, &
            upper_probability=pb, max_iterations=iteration_limit, &
            tolerance=selected_tolerance)
      end if
   end function tseriestarma_unit_root_bootstrap

   pure function tseriestarma_arma_garch_fit(series, ar_order, ma_order, &
      arch_order, garch_order, initial, max_iterations, tolerance) result(out)
      !! Jointly fit a Gaussian ARMA model with ARCH or GARCH innovations.
      real(dp), intent(in) :: series(:) !! Univariate time-series observations.
      integer, intent(in) :: ar_order !! Positive autoregressive order.
      integer, intent(in) :: ma_order !! Positive moving-average order.
      integer, intent(in) :: arch_order !! Positive ARCH order.
      integer, intent(in) :: garch_order !! Nonnegative GARCH order.
      real(dp), intent(in), optional :: initial(:) !! Initial mean, AR, MA, variance, ARCH, and GARCH values.
      integer, intent(in), optional :: max_iterations !! Maximum optimizer iterations.
      real(dp), intent(in), optional :: tolerance !! Optimizer convergence tolerance.
      type(tseriestarma_arma_garch_fit_t) :: out
      type(astsa_sarima_fit_t) :: arma_fit
      type(astsa_arma_diagnostic_t) :: diagnostic
      type(optimization_result_t) :: optimized
      real(dp), allocatable :: starting(:), selected(:), residuals(:), variance_path(:)
      real(dp) :: center, sample_variance, selected_tolerance, objective_value
      integer :: count, info, iteration_limit, offset

      count = 2 + ar_order + ma_order + arch_order + garch_order
      if (ar_order < 1 .or. ma_order < 1 .or. arch_order < 1 .or. &
         garch_order < 0 .or. size(series) <= 3*count .or. &
         .not. all(ieee_is_finite(series))) then
         out%info = 1
         return
      end if
      iteration_limit = 300
      if (present(max_iterations)) iteration_limit = max_iterations
      selected_tolerance = 1.0e-6_dp
      if (present(tolerance)) selected_tolerance = tolerance
      if (iteration_limit < 1 .or. selected_tolerance <= 0.0_dp) then
         out%info = 1
         return
      end if
      center = sum(series)/real(size(series), dp)
      sample_variance = sum((series - center)**2)/ &
         real(max(1, size(series) - 1), dp)
      sample_variance = max(sample_variance, 1.0e-6_dp)
      allocate(starting(count), source=0.0_dp)
      if (.not. present(initial)) then
         arma_fit = sarima_fit(series, ar_order, 0, ma_order, 0, 0, 0, 1, &
            include_intercept=.true., max_iterations=min(150, iteration_limit), &
            tolerance=selected_tolerance, transform_parameters=.true., &
            exact_likelihood=.true.)
         if (allocated(arma_fit%coefficients) .and. &
            size(arma_fit%coefficients) == ar_order + ma_order + 1) then
            starting(2:ar_order + 1) = arma_fit%coefficients(:ar_order)
            starting(ar_order + 2:ar_order + ma_order + 1) = &
               arma_fit%coefficients(ar_order + 1:ar_order + ma_order)
            starting(1) = arma_fit%coefficients(ar_order + ma_order + 1)* &
               (1.0_dp - sum(starting(2:ar_order + 1)))
         end if
      end if
      offset = ar_order + ma_order + 1
      starting(offset + 1) = 0.07_dp*sample_variance
      starting(offset + 2:offset + 1 + arch_order) = 0.08_dp/ &
         real(arch_order, dp)
      if (garch_order > 0) starting(offset + 2 + arch_order:) = &
         0.85_dp/real(garch_order, dp)
      if (present(initial)) then
         if (size(initial) /= count .or. .not. all(ieee_is_finite(initial))) then
            out%info = 1
            return
         end if
         starting = initial
      end if
      call evaluate_parameters(starting, residuals, variance_path, &
         objective_value, info)
      if (info /= 0) then
         out%info = 2
         return
      end if
      optimized = bfgs_minimize_fd(objective, starting, iteration_limit, &
         selected_tolerance)
      selected = starting
      if (allocated(optimized%parameters)) then
         call evaluate_parameters(optimized%parameters, residuals, variance_path, &
            objective_value, info)
         if (info == 0 .and. objective_value <= objective(starting)) &
            selected = optimized%parameters
      end if
      call evaluate_parameters(selected, residuals, variance_path, &
         objective_value, info)
      if (info /= 0) then
         out%info = 3
         return
      end if
      out%mean = selected(1)
      out%ar = selected(2:ar_order + 1)
      out%ma = selected(ar_order + 2:ar_order + ma_order + 1)
      offset = ar_order + ma_order + 1
      out%variance_intercept = selected(offset + 1)
      out%arch = selected(offset + 2:offset + 1 + arch_order)
      out%garch = selected(offset + 2 + arch_order:)
      out%residuals = residuals
      out%conditional_variance = variance_path
      out%log_likelihood = -objective_value
      out%iterations = optimized%iterations
      out%optimizer_info = optimized%info
      out%converged = optimized%converged
      diagnostic = arma_check(out%ar, out%ma)
      out%causal = diagnostic%info == 0 .and. diagnostic%causal
      out%invertible = diagnostic%info == 0 .and. diagnostic%invertible
      out%variance_stationary = sum(out%arch) + sum(out%garch) < 1.0_dp

   contains

      pure real(dp) function objective(parameters) result(value)
         !! Return the joint Gaussian ARMA-GARCH negative log likelihood.
         real(dp), intent(in) :: parameters(:) !! Physical null-model parameters.
         real(dp), allocatable :: trial_residuals(:), trial_variance(:)
         integer :: status

         call evaluate_parameters(parameters, trial_residuals, trial_variance, &
            value, status)
         if (status /= 0) value = 1.0e20_dp + dot_product(parameters, parameters)
      end function objective

      pure subroutine evaluate_parameters(parameters, innovations, variance_values, &
         value, status)
         !! Evaluate one physical ARMA-GARCH parameter vector.
         real(dp), intent(in) :: parameters(:) !! Mean, AR, MA, variance, ARCH, and GARCH values.
         real(dp), allocatable, intent(out) :: innovations(:) !! Conditional innovations.
         real(dp), allocatable, intent(out) :: variance_values(:) !! Conditional variances.
         real(dp), intent(out) :: value !! Gaussian negative log likelihood.
         integer, intent(out) :: status !! Status code; zero indicates success.
         type(astsa_arma_diagnostic_t) :: checked
         real(dp), allocatable :: ar(:), ma(:), arch(:), garch(:)
         real(dp) :: intercept, persistence, prediction, unconditional
         integer :: i, j, variance_offset

         status = 0
         value = huge(1.0_dp)
         if (size(parameters) /= count .or. &
            .not. all(ieee_is_finite(parameters))) then
            status = 1
            return
         end if
         ar = parameters(2:ar_order + 1)
         ma = parameters(ar_order + 2:ar_order + ma_order + 1)
         variance_offset = ar_order + ma_order + 1
         intercept = parameters(variance_offset + 1)
         arch = parameters(variance_offset + 2:variance_offset + 1 + arch_order)
         garch = parameters(variance_offset + 2 + arch_order:)
         persistence = sum(arch) + sum(garch)
         checked = arma_check(ar, ma)
         if (checked%info /= 0 .or. .not. checked%causal .or. &
            .not. checked%invertible .or. intercept <= tiny(1.0_dp) .or. &
            any(arch < 0.0_dp) .or. any(garch < 0.0_dp) .or. &
            persistence >= 0.999_dp) then
            status = 2
            return
         end if
         unconditional = intercept/max(1.0e-6_dp, 1.0_dp - persistence)
         allocate(innovations(size(series)), variance_values(size(series)))
         innovations = 0.0_dp
         variance_values = unconditional
         value = 0.0_dp
         do i = 1, size(series)
            prediction = parameters(1)
            do j = 1, min(ar_order, i - 1)
               prediction = prediction + ar(j)*series(i - j)
            end do
            do j = 1, min(ma_order, i - 1)
               prediction = prediction + ma(j)*innovations(i - j)
            end do
            innovations(i) = series(i) - prediction
            variance_values(i) = intercept
            do j = 1, arch_order
               if (i > j) then
                  variance_values(i) = variance_values(i) + &
                     arch(j)*innovations(i - j)**2
               else
                  variance_values(i) = variance_values(i) + arch(j)*unconditional
               end if
            end do
            do j = 1, garch_order
               if (i > j) then
                  variance_values(i) = variance_values(i) + &
                     garch(j)*variance_values(i - j)
               else
                  variance_values(i) = variance_values(i) + garch(j)*unconditional
               end if
            end do
            if (variance_values(i) <= tiny(1.0_dp) .or. &
               .not. ieee_is_finite(variance_values(i))) then
               status = 3
               return
            end if
            value = value + 0.5_dp*(log(2.0_dp*acos(-1.0_dp)* &
               variance_values(i)) + innovations(i)**2/variance_values(i))
         end do
         if (.not. ieee_is_finite(value)) status = 3
      end subroutine evaluate_parameters

   end function tseriestarma_arma_garch_fit

   pure function tseriestarma_garch_test(series, ar_order, ma_order, &
      arch_order, garch_order, delay, threshold_range, lower_probability, &
      upper_probability, initial, max_iterations, tolerance) result(out)
      !! Test an ARMA-GARCH null against a TARMA-GARCH alternative.
      real(dp), intent(in) :: series(:) !! Univariate time-series observations.
      integer, intent(in) :: ar_order !! Positive autoregressive order.
      integer, intent(in) :: ma_order !! Positive moving-average order.
      integer, intent(in) :: arch_order !! Positive ARCH order.
      integer, intent(in) :: garch_order !! Nonnegative GARCH order.
      integer, intent(in), optional :: delay !! Positive threshold delay.
      real(dp), intent(in), optional :: threshold_range(:) !! Threshold values to test.
      real(dp), intent(in), optional :: lower_probability !! Lower threshold-search quantile.
      real(dp), intent(in), optional :: upper_probability !! Upper threshold-search quantile.
      real(dp), intent(in), optional :: initial(:) !! Initial null-model parameters.
      integer, intent(in), optional :: max_iterations !! Maximum null-fit iterations.
      real(dp), intent(in), optional :: tolerance !! Null-fit convergence tolerance.
      type(tseriestarma_garch_test_t) :: out
      real(dp), allocatable :: ordered(:), design(:, :), residual_lags(:, :)
      real(dp), allocatable :: threshold_variable(:), residuals(:), variance_values(:)
      real(dp), allocatable :: null_derivative(:, :), alternative_derivative(:, :)
      real(dp), allocatable :: null_variance_derivative(:, :)
      real(dp), allocatable :: alternative_variance_derivative(:, :)
      real(dp), allocatable :: null_standardized(:, :), null_variance_scaled(:, :)
      real(dp), allocatable :: alternative_standardized(:, :)
      real(dp), allocatable :: alternative_variance_scaled(:, :)
      real(dp), allocatable :: null_information(:, :), null_inverse(:, :)
      real(dp), allocatable :: alternative_information(:, :), cross_information(:, :)
      real(dp), allocatable :: efficient_information(:, :), efficient_inverse(:, :)
      real(dp), allocatable :: null_score(:), alternative_score(:), efficient_score(:)
      real(dp), allocatable :: indicator(:), ones(:)
      real(dp) :: pa, pb, selected_tolerance, threshold
      integer :: first, last, grid_index, info, iteration_limit, k, n, neff
      integer :: parameter_count, row, time, lag, location

      out%ar_order = ar_order
      out%ma_order = ma_order
      out%arch_order = arch_order
      out%garch_order = garch_order
      if (present(delay)) out%delay = delay
      pa = 0.25_dp
      if (present(lower_probability)) pa = lower_probability
      pb = 0.75_dp
      if (present(upper_probability)) pb = upper_probability
      out%lower_probability = pa
      out%upper_probability = pb
      parameter_count = ar_order + ma_order + 1
      out%degrees_of_freedom = parameter_count
      n = size(series)
      k = max(ar_order, max(ma_order, out%delay))
      neff = n - k - arch_order
      out%effective_observations = max(0, neff)
      if (ar_order < 1 .or. ma_order < 1 .or. arch_order < 1 .or. &
         garch_order < 0 .or. out%delay < 1 .or. neff <= parameter_count .or. &
         pa <= 0.0_dp .or. pb >= 1.0_dp .or. pa >= pb .or. &
         .not. all(ieee_is_finite(series))) then
         out%info = 1
         return
      end if
      iteration_limit = 300
      if (present(max_iterations)) iteration_limit = max_iterations
      selected_tolerance = 1.0e-6_dp
      if (present(tolerance)) selected_tolerance = tolerance
      if (present(initial)) then
         out%null_fit = tseriestarma_arma_garch_fit(series, ar_order, &
            ma_order, arch_order, garch_order, initial, iteration_limit, &
            selected_tolerance)
      else
         out%null_fit = tseriestarma_arma_garch_fit(series, ar_order, &
            ma_order, arch_order, garch_order, max_iterations=iteration_limit, &
            tolerance=selected_tolerance)
      end if
      if (out%null_fit%info /= 0 .or. .not. out%null_fit%causal .or. &
         .not. out%null_fit%invertible .or. &
         .not. out%null_fit%variance_stationary) then
         out%info = 2
         return
      end if
      if (present(threshold_range)) then
         if (size(threshold_range) < 1 .or. &
            .not. all(ieee_is_finite(threshold_range))) then
            out%info = 1
            return
         end if
         out%threshold_values = threshold_range
      else
         ordered = sorted(series(out%delay + 1:n))
         first = max(1, ceiling(real(n - out%delay, dp)*pa))
         last = min(n - out%delay, floor(real(n - out%delay, dp)*pb))
         if (first > last) then
            out%info = 1
            return
         end if
         out%threshold_values = ordered(first:last)
      end if

      allocate(design(neff, ar_order + 1), residual_lags(neff, ma_order))
      allocate(threshold_variable(neff), residuals(neff), variance_values(neff))
      design(:, 1) = 1.0_dp
      do row = 1, neff
         time = k + arch_order + row
         threshold_variable(row) = series(time - out%delay)
         residuals(row) = out%null_fit%residuals(time)
         variance_values(row) = out%null_fit%conditional_variance(time)
         do lag = 1, ar_order
            design(row, lag + 1) = series(time - lag)
         end do
         do lag = 1, ma_order
            residual_lags(row, lag) = out%null_fit%residuals(time - lag)
         end do
      end do
      allocate(null_derivative(neff, parameter_count))
      null_derivative(:, :ar_order + 1) = &
         recursive_filter_matrix(-design, -out%null_fit%ma)
      null_derivative(:, ar_order + 2:) = &
         recursive_filter_matrix(residual_lags, -out%null_fit%ma)
      allocate(ones(neff), source=1.0_dp)
      null_variance_derivative = variance_score_derivative(ones)
      null_standardized = null_derivative/ &
         spread(sqrt(variance_values), 2, parameter_count)
      null_variance_scaled = null_variance_derivative/ &
         spread(variance_values, 2, parameter_count)
      allocate(null_information(parameter_count, parameter_count))
      allocate(null_inverse(parameter_count, parameter_count))
      null_information = (matmul(transpose(null_standardized), &
         null_standardized) + 0.5_dp*matmul(transpose(null_variance_scaled), &
         null_variance_scaled))/real(neff, dp)
      call symmetric_pseudoinverse(null_information, null_inverse, info)
      if (info /= 0) then
         out%info = 3
         return
      end if
      allocate(null_score(parameter_count))
      null_score = sum(-null_derivative*spread(residuals/variance_values, 2, &
         parameter_count) + 0.5_dp*null_variance_scaled*spread( &
         residuals**2/variance_values - 1.0_dp, 2, parameter_count), dim=1)/ &
         sqrt(real(neff, dp))

      allocate(out%test_values(size(out%threshold_values)))
      allocate(alternative_derivative(neff, parameter_count))
      allocate(alternative_information(parameter_count, parameter_count))
      allocate(cross_information(parameter_count, parameter_count))
      allocate(efficient_information(parameter_count, parameter_count))
      allocate(efficient_inverse(parameter_count, parameter_count))
      allocate(alternative_score(parameter_count), efficient_score(parameter_count))
      allocate(indicator(neff))
      do grid_index = 1, size(out%threshold_values)
         threshold = out%threshold_values(grid_index)
         indicator = merge(1.0_dp, 0.0_dp, threshold_variable <= threshold)
         alternative_derivative(:, :ar_order + 1) = &
            recursive_filter_matrix(-design*spread(indicator, 2, &
            ar_order + 1), -out%null_fit%ma)
         alternative_derivative(:, ar_order + 2:) = &
            recursive_filter_matrix(residual_lags*spread(indicator, 2, &
            ma_order), -out%null_fit%ma)
         alternative_variance_derivative = variance_score_derivative(indicator)
         alternative_standardized = alternative_derivative/ &
            spread(sqrt(variance_values), 2, parameter_count)
         alternative_variance_scaled = alternative_variance_derivative/ &
            spread(variance_values, 2, parameter_count)
         alternative_information = (matmul(transpose(alternative_standardized), &
            alternative_standardized) + 0.5_dp*matmul( &
            transpose(alternative_variance_scaled), &
            alternative_variance_scaled))/real(neff, dp)
         cross_information = (matmul(transpose(alternative_standardized), &
            null_standardized) + 0.5_dp*matmul( &
            transpose(alternative_variance_scaled), &
            null_variance_scaled))/real(neff, dp)
         alternative_score = sum(-alternative_derivative*spread( &
            residuals/variance_values, 2, parameter_count) + &
            0.5_dp*alternative_variance_scaled*spread( &
            residuals**2/variance_values - 1.0_dp, 2, parameter_count), &
            dim=1)/sqrt(real(neff, dp))
         efficient_score = alternative_score - &
            matmul(cross_information, matmul(null_inverse, null_score))
         efficient_information = alternative_information - &
            matmul(cross_information, matmul(null_inverse, &
            transpose(cross_information)))
         call symmetric_pseudoinverse(efficient_information, &
            efficient_inverse, info)
         if (info /= 0) then
            out%info = 3
            return
         end if
         out%test_values(grid_index) = max(0.0_dp, &
            dot_product(efficient_score, matmul(efficient_inverse, &
            efficient_score)))
      end do
      out%statistic = maxval(out%test_values)
      location = maxloc(out%test_values, dim=1)
      out%maximizing_threshold = out%threshold_values(location)
      out%lower_regime_proportion = real(count(series(out%delay + 1:n) <= &
         out%maximizing_threshold), dp)/real(n - out%delay, dp)
      out%standardized_variance = sum(residuals**2/variance_values)/ &
         real(neff, dp)
      out%critical_probabilities = tseriestarma_garch_critical_probabilities()
      out%critical_values = tseriestarma_garch_critical_values(parameter_count, pa)

   contains

      pure function variance_score_derivative(weights) result(derivative)
         !! Propagate mean-parameter effects through ARCH and GARCH variance recursions.
         real(dp), intent(in) :: weights(:) !! Threshold indicators or unit weights.
         real(dp) :: derivative(neff, parameter_count)
         real(dp), allocatable :: raw(:, :), lag_design(:, :), lag_residuals(:, :)
         real(dp), allocatable :: ar_part(:, :), ma_part(:, :)
         integer :: arch_lag, column, local_row, local_time

         allocate(raw(neff, parameter_count), source=0.0_dp)
         allocate(lag_design(neff, ar_order + 1))
         allocate(lag_residuals(neff, ma_order))
         do arch_lag = 1, arch_order
            lag_design(:, 1) = 1.0_dp
            do local_row = 1, neff
               local_time = k + arch_order + local_row
               do column = 1, ar_order
                  lag_design(local_row, column + 1) = &
                     series(local_time - arch_lag - column)
               end do
               do column = 1, ma_order
                  lag_residuals(local_row, column) = &
                     out%null_fit%residuals(local_time - arch_lag - column)
               end do
            end do
            ar_part = recursive_filter_matrix(-lag_design*spread(weights, 2, &
               ar_order + 1), -out%null_fit%ma)
            ma_part = recursive_filter_matrix(lag_residuals*spread(weights, 2, &
               ma_order), -out%null_fit%ma)
            do local_row = 1, neff
               local_time = k + arch_order + local_row
               raw(local_row, :) = raw(local_row, :) + 2.0_dp* &
                  out%null_fit%arch(arch_lag)* &
                  out%null_fit%residuals(local_time - arch_lag)* &
                  [ar_part(local_row, :), ma_part(local_row, :)]
            end do
         end do
         derivative = recursive_filter_matrix(raw, out%null_fit%garch)
      end function variance_score_derivative

   end function tseriestarma_garch_test

   pure function tseriestarma_tarma_test(series, ar_order, ma_order, delay, &
      ma_fixed, threshold_range, lower_probability, upper_probability, &
      initial, max_iterations, tolerance) result(out)
      !! Compute classic and robust ARMA-versus-TARMA supremum LM tests.
      real(dp), intent(in) :: series(:) !! Univariate time-series observations.
      integer, intent(in) :: ar_order !! Nonnegative order of the null autoregression.
      integer, intent(in) :: ma_order !! Nonnegative order of the null moving average.
      integer, intent(in), optional :: delay !! Positive delay of the threshold variable.
      logical, intent(in), optional :: ma_fixed !! Test only AR threshold effects when true.
      real(dp), intent(in), optional :: threshold_range(:) !! Threshold values to test.
      real(dp), intent(in), optional :: lower_probability !! Lower threshold-search quantile.
      real(dp), intent(in), optional :: upper_probability !! Upper threshold-search quantile.
      real(dp), intent(in), optional :: initial(:) !! Initial AR, MA, and intercept estimates.
      integer, intent(in), optional :: max_iterations !! Maximum null-fit iterations.
      real(dp), intent(in), optional :: tolerance !! Null-fit convergence tolerance.
      type(tseriestarma_tarma_test_t) :: out
      real(dp), allocatable :: design(:, :), threshold_variable(:), ordered_thresholds(:)
      real(dp), allocatable :: residuals(:), residual_lags(:, :), ma(:)
      real(dp), allocatable :: ar_derivative(:, :), ma_derivative(:, :)
      real(dp), allocatable :: threshold_ar(:, :), threshold_ma(:, :)
      real(dp), allocatable :: null_derivative(:, :), alternative_derivative(:, :)
      real(dp), allocatable :: null_information(:, :), null_inverse(:, :)
      real(dp), allocatable :: cross_information(:, :), alternative_information(:, :)
      real(dp), allocatable :: efficient_information(:, :), efficient_inverse(:, :)
      real(dp), allocatable :: score_terms(:, :), robust_covariance(:, :), robust_inverse(:, :)
      real(dp), allocatable :: score(:), robust_score(:), indicator(:)
      real(dp) :: selected_tolerance, threshold
      integer :: first, last, grid_index, info, k, n, neff, null_count
      integer :: tested_count, row, time, lag, iteration_limit
      integer :: classic_location, robust_location

      out%ar_order = ar_order
      out%ma_order = ma_order
      if (present(delay)) out%delay = delay
      if (present(ma_fixed)) out%ma_fixed = ma_fixed
      if (present(lower_probability)) out%lower_probability = lower_probability
      if (present(upper_probability)) out%upper_probability = upper_probability
      n = size(series)
      if (ar_order < 0 .or. ma_order < 0 .or. ar_order + ma_order < 1 .or. &
         out%delay < 1 .or. n <= max(ar_order, max(ma_order, out%delay)) + &
         ar_order + ma_order + 1 .or. .not. all(ieee_is_finite(series))) then
         out%info = 1
         return
      end if
      if (out%lower_probability <= 0.0_dp .or. &
         out%upper_probability >= 1.0_dp .or. &
         out%lower_probability >= out%upper_probability) then
         out%info = 2
         return
      end if
      if (present(initial)) then
         if (size(initial) /= ar_order + ma_order + 1 .or. &
            .not. all(ieee_is_finite(initial))) then
            out%info = 1
            return
         end if
      end if
      iteration_limit = 200
      if (present(max_iterations)) iteration_limit = max_iterations
      selected_tolerance = 1.0e-6_dp
      if (present(tolerance)) selected_tolerance = tolerance
      if (iteration_limit < 1 .or. selected_tolerance <= 0.0_dp) then
         out%info = 1
         return
      end if

      if (present(initial)) then
         out%noise_fit = sarima_fit(series, ar_order, 0, ma_order, 0, 0, 0, 1, &
            initial=initial, include_intercept=.true., &
            max_iterations=iteration_limit, tolerance=selected_tolerance, &
            transform_parameters=.true., exact_likelihood=.true.)
      else
         out%noise_fit = sarima_fit(series, ar_order, 0, ma_order, 0, 0, 0, 1, &
            include_intercept=.true., max_iterations=iteration_limit, &
            tolerance=selected_tolerance, transform_parameters=.true., &
            exact_likelihood=.true.)
      end if
      if (out%noise_fit%likelihood%info /= 0 .or. &
         .not. allocated(out%noise_fit%coefficients) .or. &
         .not. allocated(out%noise_fit%likelihood%residuals)) then
         out%info = 3
         return
      end if
      out%coefficients = out%noise_fit%coefficients
      out%sigma2 = out%noise_fit%likelihood%sigma2
      if (out%sigma2 <= tiny(1.0_dp)) then
         out%info = 3
         return
      end if

      if (present(threshold_range)) then
         if (size(threshold_range) < 1 .or. &
            .not. all(ieee_is_finite(threshold_range))) then
            out%info = 2
            return
         end if
         out%threshold_values = threshold_range
      else
         ordered_thresholds = sorted(series(out%delay + 1:n))
         first = max(1, ceiling(real(n - out%delay, dp)*out%lower_probability))
         last = min(n - out%delay, &
            floor(real(n - out%delay, dp)*out%upper_probability))
         if (first > last) then
            out%info = 2
            return
         end if
         out%threshold_values = ordered_thresholds(first:last)
      end if

      if (out%ma_fixed) then
         k = max(ar_order, out%delay)
         tested_count = ar_order + 1
         null_count = tested_count
      else
         k = max(ar_order, max(ma_order, out%delay))
         tested_count = ar_order + ma_order + 1
         null_count = tested_count
      end if
      neff = n - k
      out%effective_observations = neff
      out%degrees_of_freedom = tested_count
      allocate(design(neff, ar_order + 1), threshold_variable(neff), residuals(neff))
      design(:, 1) = 1.0_dp
      do row = 1, neff
         time = k + row
         threshold_variable(row) = series(time - out%delay)
         residuals(row) = out%noise_fit%likelihood%residuals(time)
         do lag = 1, ar_order
            design(row, lag + 1) = series(time - lag)
         end do
      end do
      out%residuals = residuals
      allocate(ma(ma_order))
      if (ma_order > 0) ma = out%noise_fit%coefficients(ar_order + 1: &
         ar_order + ma_order)
      ar_derivative = recursive_filter_matrix(-design, -ma)
      allocate(null_derivative(neff, null_count))
      null_derivative(:, :ar_order + 1) = ar_derivative
      if (.not. out%ma_fixed .and. ma_order > 0) then
         allocate(residual_lags(neff, ma_order))
         do row = 1, neff
            time = k + row
            do lag = 1, ma_order
               residual_lags(row, lag) = &
                  out%noise_fit%likelihood%residuals(time - lag)
            end do
         end do
         ma_derivative = recursive_filter_matrix(residual_lags, -ma)
         null_derivative(:, ar_order + 2:) = ma_derivative
      end if
      allocate(null_information(null_count, null_count), &
         null_inverse(null_count, null_count))
      null_information = matmul(transpose(null_derivative), null_derivative)
      call symmetric_pseudoinverse(null_information, null_inverse, info)
      if (info /= 0) then
         out%info = 4
         return
      end if

      allocate(out%test_values(size(out%threshold_values), 2))
      allocate(threshold_ar(neff, ar_order + 1), indicator(neff))
      allocate(alternative_derivative(neff, tested_count))
      allocate(cross_information(tested_count, null_count))
      allocate(alternative_information(tested_count, tested_count))
      allocate(efficient_information(tested_count, tested_count))
      allocate(efficient_inverse(tested_count, tested_count))
      allocate(score_terms(tested_count, neff))
      allocate(robust_covariance(tested_count, tested_count))
      allocate(robust_inverse(tested_count, tested_count))
      allocate(score(tested_count), robust_score(tested_count))
      if (.not. out%ma_fixed .and. ma_order > 0) &
         allocate(threshold_ma(neff, ma_order))
      do grid_index = 1, size(out%threshold_values)
         threshold = out%threshold_values(grid_index)
         indicator = merge(1.0_dp, 0.0_dp, threshold_variable <= threshold)
         threshold_ar = recursive_filter_matrix(-design* &
            spread(indicator, 2, ar_order + 1), -ma)
         alternative_derivative(:, :ar_order + 1) = threshold_ar
         if (.not. out%ma_fixed .and. ma_order > 0) then
            threshold_ma = recursive_filter_matrix(residual_lags* &
               spread(indicator, 2, ma_order), -ma)
            alternative_derivative(:, ar_order + 2:) = threshold_ma
         end if
         score = -sum(alternative_derivative* &
            spread(residuals, 2, tested_count), dim=1)
         alternative_information = matmul(transpose(alternative_derivative), &
            alternative_derivative)
         cross_information = matmul(transpose(alternative_derivative), &
            null_derivative)
         efficient_information = alternative_information - &
            matmul(cross_information, matmul(null_inverse, &
            transpose(cross_information)))
         call symmetric_pseudoinverse(efficient_information, &
            efficient_inverse, info)
         if (info /= 0) then
            out%info = 4
            return
         end if
         score_terms = -transpose(alternative_derivative* &
            spread(residuals, 2, tested_count)) + &
            matmul(matmul(cross_information, null_inverse), &
            transpose(null_derivative*spread(residuals, 2, null_count)))
         robust_score = sum(score_terms, dim=2)
         robust_covariance = matmul(score_terms, transpose(score_terms))
         call symmetric_pseudoinverse(robust_covariance, robust_inverse, info)
         if (info /= 0) then
            out%info = 4
            return
         end if
         out%test_values(grid_index, 1) = max(0.0_dp, &
            dot_product(score, matmul(efficient_inverse, score))/out%sigma2)
         out%test_values(grid_index, 2) = max(0.0_dp, &
            dot_product(robust_score, matmul(robust_inverse, robust_score)))
      end do

      out%statistic = maxval(out%test_values, dim=1)
      classic_location = maxloc(out%test_values(:, 1), dim=1)
      robust_location = maxloc(out%test_values(:, 2), dim=1)
      out%maximizing_threshold = [out%threshold_values(classic_location), &
         out%threshold_values(robust_location)]
      out%lower_regime_proportion(1) = real(count(series(out%delay + 1:n) <= &
         out%maximizing_threshold(1)), dp)/real(n - out%delay, dp)
      out%lower_regime_proportion(2) = real(count(series(out%delay + 1:n) <= &
         out%maximizing_threshold(2)), dp)/real(n - out%delay, dp)
   end function tseriestarma_tarma_test

   pure function tseriestarma_model(phi_lower, phi_upper, theta_lower, &
      theta_upper, ar_lags_lower, ar_lags_upper, ma_lags_lower, ma_lags_upper, &
      delay, threshold, innovation_sd_lower, innovation_sd_upper) result(out)
      !! Construct and validate a sparse two-regime TARMA model.
      real(dp), intent(in) :: phi_lower(:) !! Lower intercept followed by sparse AR coefficients.
      real(dp), intent(in) :: phi_upper(:) !! Upper intercept followed by sparse AR coefficients.
      real(dp), intent(in) :: theta_lower(:) !! Lower sparse MA coefficients.
      real(dp), intent(in) :: theta_upper(:) !! Upper sparse MA coefficients.
      integer, intent(in) :: ar_lags_lower(:) !! Lower positive AR lag indices.
      integer, intent(in) :: ar_lags_upper(:) !! Upper positive AR lag indices.
      integer, intent(in) :: ma_lags_lower(:) !! Lower positive MA lag indices.
      integer, intent(in) :: ma_lags_upper(:) !! Upper positive MA lag indices.
      integer, intent(in), optional :: delay !! Positive threshold delay.
      real(dp), intent(in), optional :: threshold !! Regime threshold.
      real(dp), intent(in), optional :: innovation_sd_lower !! Lower innovation standard deviation.
      real(dp), intent(in), optional :: innovation_sd_upper !! Upper innovation standard deviation.
      type(tseriestarma_model_t) :: out

      out%phi_lower = phi_lower
      out%phi_upper = phi_upper
      out%theta_lower = theta_lower
      out%theta_upper = theta_upper
      out%ar_lags_lower = ar_lags_lower
      out%ar_lags_upper = ar_lags_upper
      out%ma_lags_lower = ma_lags_lower
      out%ma_lags_upper = ma_lags_upper
      if (present(delay)) out%delay = delay
      if (present(threshold)) out%threshold = threshold
      if (present(innovation_sd_lower)) &
         out%innovation_sd_lower = innovation_sd_lower
      if (present(innovation_sd_upper)) &
         out%innovation_sd_upper = innovation_sd_upper
      if (.not. valid_model(out)) out%info = 1
   end function tseriestarma_model

   pure function tseriestarma_evaluate(model, series) result(out)
      !! Evaluate conditional TARMA residuals with zero presample innovations.
      type(tseriestarma_model_t), intent(in) :: model !! TARMA model parameters.
      real(dp), intent(in) :: series(:) !! Observed time series.
      type(tseriestarma_evaluation_t) :: out
      real(dp) :: prediction, nan_value
      integer :: time, lag, start

      start = model_start(model) + 1
      if (.not. valid_model(model) .or. size(series) < start .or. &
         .not. all(ieee_is_finite(series))) then
         out%info = 1
         return
      end if
      allocate(out%fitted(size(series)), out%residuals(size(series)))
      allocate(out%regime(size(series)), source=0)
      nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
      out%fitted = nan_value
      out%residuals = 0.0_dp
      do time = start, size(series)
         if (series(time - model%delay) <= model%threshold) then
            out%regime(time) = 1
            prediction = model%phi_lower(1)
            do lag = 1, size(model%ar_lags_lower)
               prediction = prediction + model%phi_lower(lag + 1)* &
                  series(time - model%ar_lags_lower(lag))
            end do
            do lag = 1, size(model%ma_lags_lower)
               prediction = prediction + model%theta_lower(lag)* &
                  out%residuals(time - model%ma_lags_lower(lag))
            end do
         else
            out%regime(time) = 2
            prediction = model%phi_upper(1)
            do lag = 1, size(model%ar_lags_upper)
               prediction = prediction + model%phi_upper(lag + 1)* &
                  series(time - model%ar_lags_upper(lag))
            end do
            do lag = 1, size(model%ma_lags_upper)
               prediction = prediction + model%theta_upper(lag)* &
                  out%residuals(time - model%ma_lags_upper(lag))
            end do
         end if
         out%fitted(time) = prediction
         out%residuals(time) = series(time) - prediction
      end do
      out%start_index = start
      out%effective_observations = size(series) - start + 1
      out%rss = sum(out%residuals(start:)**2)
   end function tseriestarma_evaluate

   pure function tseriestarma_simulate_from_standard(model, standard, burnin) &
      result(out)
      !! Simulate from caller-supplied standard innovations and discard burn-in.
      type(tseriestarma_model_t), intent(in) :: model !! TARMA model parameters.
      real(dp), intent(in) :: standard(:) !! Standard innovations including burn-in.
      integer, intent(in), optional :: burnin !! Leading simulated values to discard.
      type(tseriestarma_simulation_t) :: out
      real(dp), allocatable :: series(:), innovations(:), draws(:)
      integer, allocatable :: regimes(:)
      real(dp) :: value, scale
      integer :: discarded, start, total, time, lag

      discarded = 0
      if (present(burnin)) discarded = burnin
      start = model_start(model)
      total = start + size(standard)
      if (.not. valid_model(model) .or. discarded < 0 .or. &
         discarded >= size(standard) .or. &
         .not. all(ieee_is_finite(standard))) then
         out%info = 1
         return
      end if
      allocate(series(total), source=0.0_dp)
      allocate(innovations(total), source=0.0_dp)
      allocate(draws(total), source=0.0_dp)
      allocate(regimes(total), source=0)
      draws(start + 1:) = standard
      do time = start + 1, total
         if (series(time - model%delay) <= model%threshold) then
            regimes(time) = 1
            scale = model%innovation_sd_lower
            value = model%phi_lower(1)
            do lag = 1, size(model%ar_lags_lower)
               value = value + model%phi_lower(lag + 1)* &
                  series(time - model%ar_lags_lower(lag))
            end do
            do lag = 1, size(model%ma_lags_lower)
               value = value + scale*model%theta_lower(lag)* &
                  draws(time - model%ma_lags_lower(lag))
            end do
         else
            regimes(time) = 2
            scale = model%innovation_sd_upper
            value = model%phi_upper(1)
            do lag = 1, size(model%ar_lags_upper)
               value = value + model%phi_upper(lag + 1)* &
                  series(time - model%ar_lags_upper(lag))
            end do
            do lag = 1, size(model%ma_lags_upper)
               value = value + scale*model%theta_upper(lag)* &
                  draws(time - model%ma_lags_upper(lag))
            end do
         end if
         innovations(time) = scale*draws(time)
         series(time) = value + innovations(time)
      end do
      out%model = model
      out%burnin = discarded
      out%observations = series(start + discarded + 1:)
      out%innovations = innovations(start + discarded + 1:)
      out%regime = regimes(start + discarded + 1:)
   end function tseriestarma_simulate_from_standard

   function tseriestarma_simulate(model, observations, burnin) result(out)
      !! Simulate TARMA observations using the shared standard-normal generator.
      type(tseriestarma_model_t), intent(in) :: model !! TARMA model parameters.
      integer, intent(in) :: observations !! Number of retained observations.
      integer, intent(in), optional :: burnin !! Number of discarded observations.
      type(tseriestarma_simulation_t) :: out
      real(dp), allocatable :: standard(:)
      integer :: discarded, draw

      discarded = 500
      if (present(burnin)) discarded = burnin
      if (observations < 1 .or. discarded < 0) then
         out%info = 1
         return
      end if
      allocate(standard(observations + discarded))
      do draw = 1, size(standard)
         standard(draw) = random_standard_normal()
      end do
      out = tseriestarma_simulate_from_standard(model, standard, discarded)
   end function tseriestarma_simulate

   pure function tseriestarma_initial_parameters(series, ar_lags_lower, &
      ar_lags_upper, ma_lags_lower, ma_lags_upper, threshold, delay) &
      result(parameters)
      !! Obtain regimewise AR least-squares starts with zero MA coefficients.
      real(dp), intent(in) :: series(:) !! Observed time series.
      integer, intent(in) :: ar_lags_lower(:) !! Lower positive AR lag indices.
      integer, intent(in) :: ar_lags_upper(:) !! Upper positive AR lag indices.
      integer, intent(in) :: ma_lags_lower(:) !! Lower positive MA lag indices.
      integer, intent(in) :: ma_lags_upper(:) !! Upper positive MA lag indices.
      real(dp), intent(in) :: threshold !! Fixed threshold used to classify rows.
      integer, intent(in) :: delay !! Positive threshold delay.
      real(dp), allocatable :: parameters(:)
      real(dp), allocatable :: beta_lower(:), beta_upper(:)
      integer :: offset

      call regime_ar_start(series, ar_lags_lower, threshold, delay, 1, beta_lower)
      call regime_ar_start(series, ar_lags_upper, threshold, delay, 2, beta_upper)
      allocate(parameters(2 + size(ar_lags_lower) + size(ar_lags_upper) + &
         size(ma_lags_lower) + size(ma_lags_upper)), source=0.0_dp)
      parameters(:size(beta_lower)) = beta_lower
      offset = size(beta_lower)
      parameters(offset + 1:offset + size(beta_upper)) = beta_upper
   end function tseriestarma_initial_parameters

   pure function tseriestarma_fit(series, ar_lags_lower, ar_lags_upper, &
      ma_lags_lower, ma_lags_upper, threshold, delay, lower_probability, &
      upper_probability, initial, max_iterations, tolerance) result(out)
      !! Fit sparse TARMA parameters by conditional least squares.
      real(dp), intent(in) :: series(:) !! Observed time series.
      integer, intent(in) :: ar_lags_lower(:) !! Lower positive AR lag indices.
      integer, intent(in) :: ar_lags_upper(:) !! Upper positive AR lag indices.
      integer, intent(in) :: ma_lags_lower(:) !! Lower positive MA lag indices.
      integer, intent(in) :: ma_lags_upper(:) !! Upper positive MA lag indices.
      real(dp), intent(in), optional :: threshold !! Fixed threshold; absent profiles the central sample.
      integer, intent(in), optional :: delay !! Positive threshold delay.
      real(dp), intent(in), optional :: lower_probability !! Lower threshold-search fraction.
      real(dp), intent(in), optional :: upper_probability !! Upper threshold-search fraction.
      real(dp), intent(in), optional :: initial(:) !! Initial packed parameter vector.
      integer, intent(in), optional :: max_iterations !! Maximum BFGS iterations per threshold.
      real(dp), intent(in), optional :: tolerance !! BFGS gradient tolerance.
      type(tseriestarma_fit_t) :: out
      type(optimization_result_t) :: optimized
      type(tseriestarma_derivatives_t) :: final_derivatives
      real(dp), allocatable :: ordered(:), starting(:), best_parameters(:)
      real(dp), allocatable :: cross(:, :), inverse(:, :)
      real(dp) :: pa, pb, best_rss, current_threshold
      integer :: selected_delay, start, effective, lower, upper, candidate
      integer :: parameter_count, best_index, inference_info, parameter

      selected_delay = 1
      if (present(delay)) selected_delay = delay
      pa = 0.25_dp
      if (present(lower_probability)) pa = lower_probability
      pb = 0.75_dp
      if (present(upper_probability)) pb = upper_probability
      start = maximum_lag(ar_lags_lower, ar_lags_upper, ma_lags_lower, &
         ma_lags_upper, selected_delay)
      effective = size(series) - start
      parameter_count = 2 + size(ar_lags_lower) + size(ar_lags_upper) + &
         size(ma_lags_lower) + size(ma_lags_upper)
      if (selected_delay < 1 .or. pa < 0.0_dp .or. pb > 1.0_dp .or. &
         pa >= pb .or. effective <= parameter_count .or. &
         any(ar_lags_lower < 1) .or. any(ar_lags_upper < 1) .or. &
         any(ma_lags_lower < 1) .or. any(ma_lags_upper < 1) .or. &
         .not. all(ieee_is_finite(series))) then
         out%info = 1
         return
      end if
      if (present(threshold)) then
         out%threshold_values = [threshold]
      else
         ordered = sorted(series(1:size(series) - selected_delay))
         lower = max(1, ceiling(real(effective - 1, dp)*pa))
         upper = min(size(ordered), floor(real(effective - 1, dp)*pb))
         if (upper < lower) then
            out%info = 1
            return
         end if
         out%threshold_values = ordered(lower:upper)
      end if
      allocate(out%profile_rss(size(out%threshold_values)), &
         source=huge(1.0_dp))
      best_rss = huge(1.0_dp)
      best_index = 0
      do candidate = 1, size(out%threshold_values)
         current_threshold = out%threshold_values(candidate)
         if (present(initial)) then
            if (size(initial) /= parameter_count) then
               out%info = 1
               return
            end if
            starting = initial
         else
            starting = tseriestarma_initial_parameters(series, ar_lags_lower, &
               ar_lags_upper, ma_lags_lower, ma_lags_upper, &
               current_threshold, selected_delay)
         end if
         optimized = bfgs_minimize(objective, objective_gradient, starting, &
            max_iterations, tolerance)
         if (.not. allocated(optimized%parameters) .or. &
            .not. ieee_is_finite(optimized%objective)) cycle
         out%profile_rss(candidate) = optimized%objective
         if (optimized%objective < best_rss) then
            best_rss = optimized%objective
            best_index = candidate
            best_parameters = optimized%parameters
            out%iterations = optimized%iterations
            out%optimizer_info = optimized%info
            out%converged = optimized%converged
         end if
      end do
      if (best_index == 0) then
         out%info = 2
         return
      end if
      out%parameters = best_parameters
      out%model = unpack_model(best_parameters, ar_lags_lower, ar_lags_upper, &
         ma_lags_lower, ma_lags_upper, selected_delay, &
         out%threshold_values(best_index))
      out%evaluation = tseriestarma_evaluate(out%model, series)
      out%sigma2 = best_rss/real(effective, dp)
      out%model%innovation_sd_lower = sqrt(out%sigma2)
      out%model%innovation_sd_upper = sqrt(out%sigma2)
      out%aic = real(effective, dp)*log(out%sigma2) + &
         2.0_dp*real(parameter_count + 1, dp)
      out%bic = real(effective, dp)*log(out%sigma2) + &
         log(real(effective, dp))*real(parameter_count + 1, dp)
      final_derivatives = tseriestarma_residual_derivatives(out%model, series)
      cross = matmul(transpose(final_derivatives%derivative(start + 1:, :)), &
         final_derivatives%derivative(start + 1:, :))
      allocate(inverse(parameter_count, parameter_count))
      call symmetric_pseudoinverse(cross, inverse, inference_info)
      if (inference_info == 0) then
         out%covariance = out%sigma2*inverse
         allocate(out%standard_errors(parameter_count))
         do parameter = 1, parameter_count
            out%standard_errors(parameter) = &
               sqrt(max(0.0_dp, out%covariance(parameter, parameter)))
         end do
      end if

   contains

      pure function objective(parameters) result(value)
         !! Evaluate conditional RSS at the current threshold candidate.
         real(dp), intent(in) :: parameters(:) !! Packed TARMA parameters.
         real(dp) :: value
         type(tseriestarma_model_t) :: objective_model
         type(tseriestarma_evaluation_t) :: objective_evaluation

         objective_model = unpack_model(parameters, ar_lags_lower, &
            ar_lags_upper, ma_lags_lower, ma_lags_upper, selected_delay, &
            current_threshold)
         objective_evaluation = tseriestarma_evaluate(objective_model, series)
         value = objective_evaluation%rss
         if (objective_evaluation%info /= 0) value = huge(1.0_dp)
      end function objective

      pure function objective_gradient(parameters) result(gradient)
         !! Evaluate the analytic conditional least-squares gradient.
         real(dp), intent(in) :: parameters(:) !! Packed TARMA parameters.
         real(dp) :: gradient(size(parameters))
         type(tseriestarma_model_t) :: gradient_model

         gradient_model = unpack_model(parameters, ar_lags_lower, &
            ar_lags_upper, ma_lags_lower, ma_lags_upper, selected_delay, &
            current_threshold)
         gradient = tseriestarma_least_squares_gradient(gradient_model, series)
      end function objective_gradient

   end function tseriestarma_fit

   pure function tseriestarma_residual_derivatives(model, series) result(out)
      !! Differentiate the conditional innovation recursion analytically.
      type(tseriestarma_model_t), intent(in) :: model !! TARMA model parameters.
      real(dp), intent(in) :: series(:) !! Observed time series.
      type(tseriestarma_derivatives_t) :: out
      real(dp) :: prediction
      integer :: parameter_count, p1, p2, q1, start, time, lag, first

      start = model_start(model) + 1
      p1 = size(model%phi_lower)
      p2 = size(model%phi_upper)
      q1 = size(model%theta_lower)
      parameter_count = p1 + p2 + q1 + size(model%theta_upper)
      if (.not. valid_model(model) .or. size(series) < start .or. &
         .not. all(ieee_is_finite(series))) then
         out%info = 1
         return
      end if
      allocate(out%residuals(size(series)), source=0.0_dp)
      allocate(out%derivative(size(series), parameter_count), source=0.0_dp)
      do time = start, size(series)
         if (series(time - model%delay) <= model%threshold) then
            prediction = model%phi_lower(1)
            do lag = 1, q1
               out%derivative(time, :) = out%derivative(time, :) - &
                  model%theta_lower(lag)*out%derivative( &
                  time - model%ma_lags_lower(lag), :)
            end do
            out%derivative(time, 1) = out%derivative(time, 1) - 1.0_dp
            do lag = 1, size(model%ar_lags_lower)
               prediction = prediction + model%phi_lower(lag + 1)* &
                  series(time - model%ar_lags_lower(lag))
               out%derivative(time, lag + 1) = &
                  out%derivative(time, lag + 1) - &
                  series(time - model%ar_lags_lower(lag))
            end do
            first = p1 + p2
            do lag = 1, q1
               prediction = prediction + model%theta_lower(lag)* &
                  out%residuals(time - model%ma_lags_lower(lag))
               out%derivative(time, first + lag) = &
                  out%derivative(time, first + lag) - &
                  out%residuals(time - model%ma_lags_lower(lag))
            end do
         else
            prediction = model%phi_upper(1)
            do lag = 1, size(model%theta_upper)
               out%derivative(time, :) = out%derivative(time, :) - &
                  model%theta_upper(lag)*out%derivative( &
                  time - model%ma_lags_upper(lag), :)
            end do
            first = p1 + 1
            out%derivative(time, first) = &
               out%derivative(time, first) - 1.0_dp
            do lag = 1, size(model%ar_lags_upper)
               prediction = prediction + model%phi_upper(lag + 1)* &
                  series(time - model%ar_lags_upper(lag))
               out%derivative(time, first + lag) = &
                  out%derivative(time, first + lag) - &
                  series(time - model%ar_lags_upper(lag))
            end do
            first = p1 + p2 + q1
            do lag = 1, size(model%theta_upper)
               prediction = prediction + model%theta_upper(lag)* &
                  out%residuals(time - model%ma_lags_upper(lag))
               out%derivative(time, first + lag) = &
                  out%derivative(time, first + lag) - &
                  out%residuals(time - model%ma_lags_upper(lag))
            end do
         end if
         out%residuals(time) = series(time) - prediction
      end do
      out%start_index = start
   end function tseriestarma_residual_derivatives

   pure function tseriestarma_least_squares_gradient(model, series, weights) &
      result(gradient)
      !! Return the analytic conditional weighted least-squares gradient.
      type(tseriestarma_model_t), intent(in) :: model !! TARMA model parameters.
      real(dp), intent(in) :: series(:) !! Observed time series.
      real(dp), intent(in), optional :: weights(:) !! Nonnegative observation weights.
      real(dp), allocatable :: gradient(:)
      type(tseriestarma_derivatives_t) :: derivatives
      real(dp), allocatable :: selected_weights(:)
      integer :: start

      derivatives = tseriestarma_residual_derivatives(model, series)
      if (derivatives%info /= 0) then
         allocate(gradient(0))
         return
      end if
      allocate(selected_weights(size(series)), source=1.0_dp)
      if (present(weights)) then
         if (size(weights) /= size(series) .or. any(weights < 0.0_dp)) then
            allocate(gradient(0))
            return
         end if
         selected_weights = weights
      end if
      start = derivatives%start_index
      gradient = 2.0_dp*matmul(transpose(derivatives%derivative(start:, :)), &
         selected_weights(start:)*derivatives%residuals(start:))
   end function tseriestarma_least_squares_gradient

   pure function tseriestarma_robust_weights(residuals, alpha, &
      innovation_family, student_degrees) result(weights)
      !! Compute normalized density-power IRLS weights for TARMA innovations.
      real(dp), intent(in) :: residuals(:) !! Conditional innovations.
      real(dp), intent(in) :: alpha !! Nonnegative density-power tuning parameter.
      integer, intent(in), optional :: innovation_family !! Normal or Student innovation code.
      real(dp), intent(in), optional :: student_degrees !! Student degrees of freedom.
      real(dp), allocatable :: weights(:)
      real(dp) :: scale2, degrees, exponent, maximum_weight
      integer :: family

      family = tseriestarma_innovation_normal
      if (present(innovation_family)) family = innovation_family
      degrees = 5.0_dp
      if (present(student_degrees)) degrees = student_degrees
      allocate(weights(size(residuals)), source=0.0_dp)
      if (size(residuals) < 1 .or. alpha < 0.0_dp .or. &
         (family /= tseriestarma_innovation_normal .and. &
         family /= tseriestarma_innovation_student) .or. &
         (family == tseriestarma_innovation_student .and. degrees <= 2.0_dp)) return
      scale2 = max(tiny(1.0_dp), &
         sum(residuals**2)/real(size(residuals), dp))
      select case (family)
      case (tseriestarma_innovation_normal)
         weights = exp(-alpha*residuals**2/(2.0_dp*scale2))/scale2
      case (tseriestarma_innovation_student)
         exponent = 0.5_dp*alpha*(degrees + 1.0_dp)
         weights = (degrees + 1.0_dp)/(degrees + residuals**2)* &
            (1.0_dp + residuals**2/degrees)**(-exponent)
      end select
      maximum_weight = maxval(weights)
      if (maximum_weight > 0.0_dp) weights = weights/maximum_weight
   end function tseriestarma_robust_weights

   pure function tseriestarma_robust_fit(series, initial_fit, alpha, &
      innovation_family, student_degrees, trim_probabilities, &
      max_irls_iterations, irls_tolerance, max_iterations, tolerance) &
      result(out)
      !! Refine an ordinary TARMA fit by density-power IRLS and sandwich inference.
      real(dp), intent(in) :: series(:) !! Observed time series.
      type(tseriestarma_fit_t), intent(in) :: initial_fit !! Ordinary fit supplying model and threshold.
      real(dp), intent(in) :: alpha !! Nonnegative density-power tuning parameter.
      integer, intent(in), optional :: innovation_family !! Normal or Student innovation code.
      real(dp), intent(in), optional :: student_degrees !! Fixed Student degrees of freedom.
      real(dp), intent(in), optional :: trim_probabilities(:) !! Two initial response-trimming probabilities.
      integer, intent(in), optional :: max_irls_iterations !! Maximum IRLS updates.
      real(dp), intent(in), optional :: irls_tolerance !! Relative parameter-change tolerance.
      integer, intent(in), optional :: max_iterations !! Maximum BFGS iterations per IRLS update.
      real(dp), intent(in), optional :: tolerance !! BFGS gradient tolerance.
      type(tseriestarma_fit_t) :: out
      type(tseriestarma_model_t) :: current_model
      type(tseriestarma_derivatives_t) :: derivatives
      type(optimization_result_t) :: optimized
      real(dp), allocatable :: parameters(:), new_parameters(:), effective_weights(:)
      real(dp), allocatable :: final_weights(:)
      real(dp), allocatable :: ordered(:), scores(:, :), hessian(:, :)
      real(dp), allocatable :: inverse(:, :), meat(:, :)
      real(dp) :: selected_alpha, degrees, change, scale2, lower_value, upper_value
      real(dp) :: irls_tol
      integer :: family, limit, iteration, start, effective, info, parameter_count

      family = tseriestarma_innovation_normal
      if (present(innovation_family)) family = innovation_family
      degrees = 5.0_dp
      if (present(student_degrees)) degrees = student_degrees
      selected_alpha = alpha
      limit = 100
      if (present(max_irls_iterations)) limit = max_irls_iterations
      irls_tol = 1.0e-4_dp
      if (present(irls_tolerance)) irls_tol = irls_tolerance
      if (initial_fit%info /= 0 .or. .not. valid_model(initial_fit%model) .or. &
         selected_alpha < 0.0_dp .or. limit < 1 .or. irls_tol <= 0.0_dp .or. &
         (family /= tseriestarma_innovation_normal .and. &
         family /= tseriestarma_innovation_student) .or. &
         (family == tseriestarma_innovation_student .and. degrees <= 2.0_dp)) then
         out%info = 1
         return
      end if
      current_model = initial_fit%model
      parameters = pack_model(current_model)
      parameter_count = size(parameters)
      start = model_start(current_model) + 1
      effective = size(series) - start + 1
      if (effective <= parameter_count) then
         out%info = 1
         return
      end if
      allocate(out%weights(size(series)), source=0.0_dp)
      out%weights(start:) = 1.0_dp
      if (present(trim_probabilities)) then
         if (size(trim_probabilities) /= 2 .or. &
            trim_probabilities(1) < 0.0_dp .or. &
            trim_probabilities(2) > 1.0_dp .or. &
            trim_probabilities(1) >= trim_probabilities(2)) then
            out%info = 1
            return
         end if
         ordered = sorted(series(start:))
         lower_value = quantile(ordered, trim_probabilities(1))
         upper_value = quantile(ordered, trim_probabilities(2))
         where (series(start:) <= lower_value .or. series(start:) >= upper_value)
            out%weights(start:) = 0.0_dp
         end where
         optimized = bfgs_minimize(weighted_objective, weighted_gradient, &
            parameters, max_iterations, tolerance)
         if (allocated(optimized%parameters)) parameters = optimized%parameters
      end if
      do iteration = 1, limit
         current_model = unpack_like(parameters, initial_fit%model)
         derivatives = tseriestarma_residual_derivatives(current_model, series)
         effective_weights = tseriestarma_robust_weights( &
            derivatives%residuals(start:), selected_alpha, family, degrees)
         out%weights = 0.0_dp
         out%weights(start:) = effective_weights
         optimized = bfgs_minimize(weighted_objective, weighted_gradient, &
            parameters, max_iterations, tolerance)
         if (.not. allocated(optimized%parameters)) then
            out%info = 2
            return
         end if
         new_parameters = optimized%parameters
         change = norm2((new_parameters - parameters)/ &
            max(1.0_dp, abs(parameters)))
         parameters = new_parameters
         if (change < irls_tol) exit
      end do
      final_weights = out%weights
      out = initial_fit
      out%weights = final_weights
      out%parameters = parameters
      out%model = unpack_like(parameters, initial_fit%model)
      out%evaluation = tseriestarma_evaluate(out%model, series)
      out%weights = 0.0_dp
      out%weights(start:) = tseriestarma_robust_weights( &
         out%evaluation%residuals(start:), selected_alpha, family, degrees)
      scale2 = sum(out%evaluation%residuals(start:)**2)/real(effective, dp)
      out%sigma2 = scale2/max(tiny(1.0_dp), 1.0_dp - 0.5_dp*selected_alpha)
      out%model%innovation_sd_lower = sqrt(out%sigma2)
      out%model%innovation_sd_upper = sqrt(out%sigma2)
      out%evaluation%rss = sum(out%weights(start:)* &
         out%evaluation%residuals(start:)**2)
      out%robust = .true.
      out%robustness_alpha = selected_alpha
      out%innovation_family = family
      if (family == tseriestarma_innovation_student) out%student_degrees = degrees
      out%irls_iterations = min(iteration, limit)
      out%iterations = optimized%iterations
      out%optimizer_info = optimized%info
      out%converged = change < irls_tol
      if (allocated(out%covariance)) deallocate(out%covariance)
      if (allocated(out%standard_errors)) deallocate(out%standard_errors)
      derivatives = tseriestarma_residual_derivatives(out%model, series)
      scores = spread(2.0_dp*out%weights(start:)* &
         derivatives%residuals(start:), 2, parameter_count)* &
         derivatives%derivative(start:, :)
      meat = matmul(transpose(scores), scores)/real(effective, dp)
      hessian = finite_difference_hessian(weighted_objective, parameters)/ &
         real(effective, dp)
      allocate(inverse(parameter_count, parameter_count))
      call symmetric_pseudoinverse(hessian, inverse, info)
      if (info /= 0) then
         hessian = 2.0_dp*matmul(transpose( &
            derivatives%derivative(start:, :)), &
            spread(out%weights(start:), 2, parameter_count)* &
            derivatives%derivative(start:, :))/real(effective, dp)
         call symmetric_pseudoinverse(hessian, inverse, info)
      end if
      if (info == 0) then
         out%covariance = matmul(matmul(inverse, meat), inverse)/ &
            real(effective, dp)
         allocate(out%standard_errors(parameter_count))
         do iteration = 1, parameter_count
            out%standard_errors(iteration) = &
               sqrt(max(0.0_dp, out%covariance(iteration, iteration)))
         end do
      end if
      out%info = 0

   contains

      pure function weighted_objective(candidate) result(value)
         !! Evaluate the current frozen-weight IRLS objective.
         real(dp), intent(in) :: candidate(:) !! Packed TARMA parameters.
         real(dp) :: value
         type(tseriestarma_model_t) :: objective_model
         type(tseriestarma_evaluation_t) :: objective_evaluation

         objective_model = unpack_like(candidate, initial_fit%model)
         objective_evaluation = tseriestarma_evaluate(objective_model, series)
         value = sum(out%weights(start:)* &
            objective_evaluation%residuals(start:)**2)
      end function weighted_objective

      pure function weighted_gradient(candidate) result(gradient)
         !! Evaluate the analytic frozen-weight IRLS gradient.
         real(dp), intent(in) :: candidate(:) !! Packed TARMA parameters.
         real(dp) :: gradient(size(candidate))
         type(tseriestarma_model_t) :: gradient_model
         type(tseriestarma_derivatives_t) :: gradient_derivatives

         gradient_model = unpack_like(candidate, initial_fit%model)
         gradient_derivatives = &
            tseriestarma_residual_derivatives(gradient_model, series)
         gradient = 2.0_dp*matmul(transpose( &
            gradient_derivatives%derivative(start:, :)), &
            out%weights(start:)*gradient_derivatives%residuals(start:))
      end function weighted_gradient

   end function tseriestarma_robust_fit

   pure function tseriestarma_fit2(series, common_ar_lags, lower_ar_lags, &
      upper_ar_lags, ma_order, seasonal_ma_order, period, threshold, delay, &
      lower_probability, upper_probability, threshold_variable, &
      regime_intercepts, regressors, initial, max_iterations, tolerance) &
      result(out)
      !! Fit common AR/MA errors and regime-specific threshold regressions by exact ML.
      real(dp), intent(in) :: series(:) !! Observed time series.
      integer, intent(in) :: common_ar_lags(:) !! Positive AR lags shared by both regimes.
      integer, intent(in) :: lower_ar_lags(:) !! Positive lower-regime threshold AR lags.
      integer, intent(in) :: upper_ar_lags(:) !! Positive upper-regime threshold AR lags.
      integer, intent(in) :: ma_order !! Nonseasonal common MA order.
      integer, intent(in), optional :: seasonal_ma_order !! Seasonal common MA order.
      integer, intent(in), optional :: period !! Seasonal period when seasonal MA terms are used.
      real(dp), intent(in), optional :: threshold !! Fixed threshold; absent profiles the central sample.
      integer, intent(in), optional :: delay !! Positive internal threshold delay.
      real(dp), intent(in), optional :: lower_probability !! Lower threshold-search fraction.
      real(dp), intent(in), optional :: upper_probability !! Upper threshold-search fraction.
      real(dp), intent(in), optional :: threshold_variable(:) !! Externally supplied aligned threshold variable.
      logical, intent(in), optional :: regime_intercepts !! Use separate rather than common intercepts.
      real(dp), intent(in), optional :: regressors(:, :) !! Additional contemporaneous regressors.
      real(dp), intent(in), optional :: initial(:) !! Initial MA, seasonal MA, and regression parameters.
      integer, intent(in), optional :: max_iterations !! Maximum exact-likelihood BFGS iterations.
      real(dp), intent(in), optional :: tolerance !! Optimizer gradient tolerance.
      type(tseriestarma_fit2_t) :: out
      type(astsa_sarima_fit_t) :: fitted_noise, best_noise
      real(dp), allocatable :: threshold_data(:), ordered(:), design(:, :)
      real(dp), allocatable :: starting(:), coefficients(:), ma_values(:)
      real(dp), allocatable :: lower_phi(:), upper_phi(:)
      integer, allocatable :: ma_lags(:)
      real(dp) :: pa, pb, best_aic, candidate_threshold, nan_value
      real(dp) :: variance_lower, variance_upper
      integer :: selected_delay, seasonal_order, selected_period, first
      integer :: effective, lower, upper, candidate, best_index, column, offset
      integer :: design_columns, regression_count, parameter_count, time
      integer :: lower_count, upper_count
      logical :: separate_intercepts

      selected_delay = 1
      if (present(delay)) selected_delay = delay
      seasonal_order = 0
      if (present(seasonal_ma_order)) seasonal_order = seasonal_ma_order
      selected_period = 1
      if (present(period)) selected_period = period
      separate_intercepts = .true.
      if (present(regime_intercepts)) separate_intercepts = regime_intercepts
      pa = 0.25_dp
      if (present(lower_probability)) pa = lower_probability
      pb = 0.75_dp
      if (present(upper_probability)) pb = upper_probability
      regression_count = 0
      if (present(regressors)) regression_count = size(regressors, 2)
      first = max(selected_delay, max_or_zero(common_ar_lags), &
         max_or_zero(lower_ar_lags), max_or_zero(upper_ar_lags)) + 1
      effective = size(series) - first + 1
      design_columns = size(common_ar_lags) + size(lower_ar_lags) + &
         size(upper_ar_lags) + regression_count + merge(2, 1, separate_intercepts)
      parameter_count = ma_order + seasonal_order + design_columns
      if (size(series) < 2 .or. selected_delay < 1 .or. ma_order < 0 .or. &
         seasonal_order < 0 .or. &
         (seasonal_order > 0 .and. selected_period < 2) .or. &
         pa < 0.0_dp .or. pb > 1.0_dp .or. pa >= pb .or. &
         effective <= parameter_count .or. any(common_ar_lags < 1) .or. &
         any(lower_ar_lags < 1) .or. any(upper_ar_lags < 1) .or. &
         sets_overlap(common_ar_lags, lower_ar_lags) .or. &
         sets_overlap(common_ar_lags, upper_ar_lags) .or. &
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
      allocate(threshold_data(size(series)), source=0.0_dp)
      if (present(threshold_variable)) then
         if (size(threshold_variable) /= size(series) .or. &
            .not. all(ieee_is_finite(threshold_variable))) then
            out%info = 1
            return
         end if
         threshold_data = threshold_variable
         out%external_threshold = .true.
      else
         threshold_data(first:) = series(first - selected_delay: &
            size(series) - selected_delay)
      end if
      if (present(threshold)) then
         out%threshold_values = [threshold]
      else
         ordered = sorted(threshold_data(first:))
         lower = max(1, ceiling(real(effective - 1, dp)*pa))
         upper = min(effective, floor(real(effective - 1, dp)*pb))
         if (upper < lower) then
            out%info = 1
            return
         end if
         out%threshold_values = ordered(lower:upper)
      end if
      allocate(out%profile_aic(size(out%threshold_values)), &
         source=huge(1.0_dp))
      if (present(initial)) then
         if (size(initial) /= parameter_count) then
            out%info = 1
            return
         end if
         starting = initial
      end if
      best_aic = huge(1.0_dp)
      best_index = 0
      do candidate = 1, size(out%threshold_values)
         candidate_threshold = out%threshold_values(candidate)
         design = fit2_design(candidate_threshold)
         if (allocated(starting)) then
            fitted_noise = sarima_fit(series(first:), 0, 0, ma_order, 0, 0, &
               seasonal_order, selected_period, initial=starting, &
               include_intercept=.false., regressors=design, &
               max_iterations=max_iterations, tolerance=tolerance, &
               transform_parameters=.true., exact_likelihood=.true.)
         else
            fitted_noise = sarima_fit(series(first:), 0, 0, ma_order, 0, 0, &
               seasonal_order, selected_period, &
               include_intercept=.false., regressors=design, &
               max_iterations=max_iterations, tolerance=tolerance, &
               transform_parameters=.true., exact_likelihood=.true.)
         end if
         if (fitted_noise%likelihood%info /= 0) cycle
         out%profile_aic(candidate) = fitted_noise%likelihood%aic
         starting = fitted_noise%coefficients
         if (out%profile_aic(candidate) < best_aic) then
            best_aic = out%profile_aic(candidate)
            best_index = candidate
            best_noise = fitted_noise
         end if
      end do
      if (best_index == 0) then
         out%info = 2
         return
      end if
      out%threshold = out%threshold_values(best_index)
      out%noise_fit = best_noise
      out%delay = selected_delay
      out%period = merge(selected_period, 0, seasonal_order > 0)
      out%start_index = first
      out%regime_intercepts = separate_intercepts
      out%common_ar_lags = common_ar_lags
      out%lower_ar_lags = lower_ar_lags
      out%upper_ar_lags = upper_ar_lags
      coefficients = best_noise%coefficients
      out%regular_ma = coefficients(:ma_order)
      offset = ma_order
      out%seasonal_ma = coefficients(offset + 1:offset + seasonal_order)
      offset = offset + seasonal_order
      if (separate_intercepts) then
         out%common_ar = coefficients(offset + 1: &
            offset + size(common_ar_lags))
         offset = offset + size(common_ar_lags)
         out%intercept_lower = coefficients(offset + 1)
         offset = offset + 1
         out%lower_ar = coefficients(offset + 1: &
            offset + size(lower_ar_lags))
         offset = offset + size(lower_ar_lags)
         out%intercept_upper = coefficients(offset + 1)
         offset = offset + 1
         out%upper_ar = coefficients(offset + 1: &
            offset + size(upper_ar_lags))
         offset = offset + size(upper_ar_lags)
      else
         out%intercept_lower = coefficients(offset + 1)
         out%intercept_upper = out%intercept_lower
         offset = offset + 1
         out%common_ar = coefficients(offset + 1: &
            offset + size(common_ar_lags))
         offset = offset + size(common_ar_lags)
         out%lower_ar = coefficients(offset + 1: &
            offset + size(lower_ar_lags))
         offset = offset + size(lower_ar_lags)
         out%upper_ar = coefficients(offset + 1: &
            offset + size(upper_ar_lags))
         offset = offset + size(upper_ar_lags)
      end if
      out%regression_coefficients = coefficients(offset + 1: &
         offset + regression_count)
      nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
      allocate(out%fitted(size(series)), out%residuals(size(series)))
      allocate(out%regime(size(series)), source=0)
      out%fitted = nan_value
      out%residuals = nan_value
      out%fitted(first:) = best_noise%likelihood%fitted
      out%residuals(first:) = best_noise%likelihood%residuals
      where (threshold_data(first:) <= out%threshold)
         out%regime(first:) = 1
      elsewhere
         out%regime(first:) = 2
      end where
      lower_count = count(out%regime(first:) == 1)
      upper_count = effective - lower_count
      if (lower_count < 1 .or. upper_count < 1) then
         out%info = 3
         return
      end if
      variance_lower = sum(pack(out%residuals(first:)**2, &
         out%regime(first:) == 1))/real(lower_count, dp)
      variance_upper = sum(pack(out%residuals(first:)**2, &
         out%regime(first:) == 2))/real(upper_count, dp)
      out%rss = sum(out%residuals(first:)**2)
      out%aic = real(lower_count, dp)*log(max(tiny(1.0_dp), variance_lower)) + &
         real(upper_count, dp)*log(max(tiny(1.0_dp), variance_upper)) + &
         real(effective, dp)*(1.0_dp + log(2.0_dp*acos(-1.0_dp))) + &
         2.0_dp*real(parameter_count + 1, dp)
      out%bic = real(lower_count, dp)*log(max(tiny(1.0_dp), variance_lower)) + &
         real(upper_count, dp)*log(max(tiny(1.0_dp), variance_upper)) + &
         real(effective, dp)*(1.0_dp + log(2.0_dp*acos(-1.0_dp))) + &
         log(real(effective, dp))*real(parameter_count + 1, dp)
      ma_lags = [(column, column=1, ma_order), &
         (column*selected_period, column=1, seasonal_order)]
      ma_values = [out%regular_ma, out%seasonal_ma]
      lower_phi = [out%intercept_lower, out%common_ar, out%lower_ar]
      upper_phi = [out%intercept_upper, out%common_ar, out%upper_ar]
      out%forecast_model = tseriestarma_model(lower_phi, upper_phi, ma_values, &
         ma_values, [common_ar_lags, lower_ar_lags], &
         [common_ar_lags, upper_ar_lags], ma_lags, ma_lags, selected_delay, &
         out%threshold, sqrt(best_noise%likelihood%sigma2), &
         sqrt(best_noise%likelihood%sigma2))
      if (out%external_threshold) out%forecast_model%info = 2

   contains

      pure function fit2_design(candidate_value) result(matrix)
         !! Construct common and regime-specific threshold regressors.
         real(dp), intent(in) :: candidate_value !! Candidate threshold.
         real(dp), allocatable :: matrix(:, :)
         integer :: row, lag_index, source_time, design_column
         logical :: lower_regime

         allocate(matrix(effective, design_columns), source=0.0_dp)
         do row = 1, effective
            source_time = first + row - 1
            lower_regime = threshold_data(source_time) <= candidate_value
            design_column = 0
            if (.not. separate_intercepts) then
               design_column = 1
               matrix(row, design_column) = 1.0_dp
            end if
            do lag_index = 1, size(common_ar_lags)
               design_column = design_column + 1
               matrix(row, design_column) = &
                  series(source_time - common_ar_lags(lag_index))
            end do
            if (separate_intercepts) then
               design_column = design_column + 1
               matrix(row, design_column) = merge(1.0_dp, 0.0_dp, lower_regime)
            end if
            do lag_index = 1, size(lower_ar_lags)
               design_column = design_column + 1
               if (lower_regime) matrix(row, design_column) = &
                  series(source_time - lower_ar_lags(lag_index))
            end do
            if (separate_intercepts) then
               design_column = design_column + 1
               matrix(row, design_column) = merge(0.0_dp, 1.0_dp, lower_regime)
            end if
            do lag_index = 1, size(upper_ar_lags)
               design_column = design_column + 1
               if (.not. lower_regime) matrix(row, design_column) = &
                  series(source_time - upper_ar_lags(lag_index))
            end do
            do lag_index = 1, regression_count
               design_column = design_column + 1
               matrix(row, design_column) = regressors(source_time, lag_index)
            end do
         end do
      end function fit2_design

   end function tseriestarma_fit2

   pure function tseriestarma_forecast_from_standard(model, series, &
      standard_draws, probabilities) result(out)
      !! Simulate TARMA forecast paths from supplied future standard innovations.
      type(tseriestarma_model_t), intent(in) :: model !! Fitted TARMA model.
      real(dp), intent(in) :: series(:) !! Observed history through the forecast origin.
      real(dp), intent(in) :: standard_draws(:, :) !! Horizon-by-simulation standard innovations.
      real(dp), intent(in), optional :: probabilities(:) !! Requested marginal quantile probabilities.
      type(tseriestarma_forecast_t) :: out
      type(tseriestarma_evaluation_t) :: evaluated
      real(dp), allocatable :: ordered(:), zero_draws(:)
      real(dp) :: center
      integer :: horizon, simulations, draw, step, probability_index, info

      horizon = size(standard_draws, 1)
      simulations = size(standard_draws, 2)
      if (.not. valid_model(model) .or. horizon < 1 .or. simulations < 1 .or. &
         size(series) < model_start(model) .or. &
         .not. all(ieee_is_finite(series)) .or. &
         .not. all(ieee_is_finite(standard_draws))) then
         out%info = 1
         return
      end if
      if (present(probabilities)) then
         if (size(probabilities) < 1 .or. any(probabilities < 0.0_dp) .or. &
            any(probabilities > 1.0_dp)) then
            out%info = 2
            return
         end if
         out%probabilities = probabilities
      else
         out%probabilities = [0.05_dp, 0.95_dp]
      end if
      evaluated = tseriestarma_evaluate(model, series)
      if (evaluated%info /= 0) then
         out%info = 3
         return
      end if
      allocate(out%paths(horizon, simulations))
      allocate(out%regimes(horizon, simulations))
      do draw = 1, simulations
         call tarma_forecast_path(model, series, evaluated%residuals, &
            standard_draws(:, draw), out%paths(:, draw), &
            out%regimes(:, draw), info)
         if (info /= 0) then
            out%info = 100 + draw
            return
         end if
      end do
      allocate(zero_draws(horizon), source=0.0_dp)
      allocate(out%point(horizon))
      allocate(out%point_regime(horizon))
      call tarma_forecast_path(model, series, evaluated%residuals, zero_draws, &
         out%point, out%point_regime, info)
      allocate(out%mean(horizon), out%standard_deviation(horizon))
      allocate(out%quantiles(horizon, size(out%probabilities)))
      do step = 1, horizon
         out%mean(step) = sum(out%paths(step, :))/real(simulations, dp)
         center = out%mean(step)
         if (simulations > 1) then
            out%standard_deviation(step) = sqrt(sum( &
               (out%paths(step, :) - center)**2)/real(simulations - 1, dp))
         else
            out%standard_deviation(step) = 0.0_dp
         end if
         ordered = sorted(out%paths(step, :))
         do probability_index = 1, size(out%probabilities)
            out%quantiles(step, probability_index) = quantile(ordered, &
               out%probabilities(probability_index))
         end do
      end do
   end function tseriestarma_forecast_from_standard

   function tseriestarma_forecast(model, series, horizon, simulations, &
      probabilities) result(out)
      !! Simulate a TARMA forecast distribution with the shared normal generator.
      type(tseriestarma_model_t), intent(in) :: model !! Fitted TARMA model.
      real(dp), intent(in) :: series(:) !! Observed history through the forecast origin.
      integer, intent(in) :: horizon !! Positive forecast horizon.
      integer, intent(in), optional :: simulations !! Monte Carlo path count; default one thousand.
      real(dp), intent(in), optional :: probabilities(:) !! Requested marginal quantile probabilities.
      type(tseriestarma_forecast_t) :: out
      real(dp), allocatable :: standard_draws(:, :)
      integer :: count, draw, step

      count = 1000
      if (present(simulations)) count = simulations
      if (horizon < 1 .or. count < 1) then
         out%info = 1
         return
      end if
      allocate(standard_draws(horizon, count))
      do draw = 1, count
         do step = 1, horizon
            standard_draws(step, draw) = random_standard_normal()
         end do
      end do
      out = tseriestarma_forecast_from_standard(model, series, &
         standard_draws, probabilities)
   end function tseriestarma_forecast

   pure subroutine tarma_forecast_path(model, history, residual_history, &
      standard, path, regimes, info)
      !! Advance one TARMA path from observed values and reconstructed innovations.
      type(tseriestarma_model_t), intent(in) :: model !! Fitted TARMA model.
      real(dp), intent(in) :: history(:) !! Observed values through the origin.
      real(dp), intent(in) :: residual_history(:) !! Reconstructed historical innovations.
      real(dp), intent(in) :: standard(:) !! Future standard innovations.
      real(dp), intent(out) :: path(:) !! Simulated future observations.
      integer, intent(out) :: regimes(:) !! Simulated future regime labels.
      integer, intent(out) :: info !! Status code; zero indicates success.
      real(dp), allocatable :: values(:), residuals(:)
      real(dp) :: prediction, scale
      integer :: observations, time, step, lag

      observations = size(history)
      info = 0
      if (size(residual_history) /= observations .or. &
         size(path) /= size(standard) .or. size(regimes) /= size(standard)) then
         info = 1
         return
      end if
      allocate(values(observations + size(standard)))
      allocate(residuals(observations + size(standard)))
      values(:observations) = history
      residuals(:observations) = residual_history
      do step = 1, size(standard)
         time = observations + step
         if (values(time - model%delay) <= model%threshold) then
            regimes(step) = 1
            scale = model%innovation_sd_lower
            prediction = model%phi_lower(1)
            do lag = 1, size(model%ar_lags_lower)
               prediction = prediction + model%phi_lower(lag + 1)* &
                  values(time - model%ar_lags_lower(lag))
            end do
            do lag = 1, size(model%ma_lags_lower)
               prediction = prediction + model%theta_lower(lag)* &
                  residuals(time - model%ma_lags_lower(lag))
            end do
         else
            regimes(step) = 2
            scale = model%innovation_sd_upper
            prediction = model%phi_upper(1)
            do lag = 1, size(model%ar_lags_upper)
               prediction = prediction + model%phi_upper(lag + 1)* &
                  values(time - model%ar_lags_upper(lag))
            end do
            do lag = 1, size(model%ma_lags_upper)
               prediction = prediction + model%theta_upper(lag)* &
                  residuals(time - model%ma_lags_upper(lag))
            end do
         end if
         residuals(time) = scale*standard(step)
         values(time) = prediction + residuals(time)
         path(step) = values(time)
      end do
   end subroutine tarma_forecast_path

   pure subroutine regime_ar_start(series, lags, threshold, delay, regime, beta)
      !! Fit one regime's intercept and sparse AR coefficients for initialization.
      real(dp), intent(in) :: series(:) !! Observed time series.
      integer, intent(in) :: lags(:) !! Positive AR lag indices.
      real(dp), intent(in) :: threshold !! Regime threshold.
      integer, intent(in) :: delay !! Threshold delay.
      integer, intent(in) :: regime !! One for lower or two for upper regime.
      real(dp), allocatable, intent(out) :: beta(:) !! Initial regression coefficients.
      real(dp), allocatable :: design(:, :), response(:), cross(:, :), inverse(:, :)
      integer, allocatable :: indices(:)
      integer :: start, time, row, lag, info

      start = max(delay, max_or_zero(lags)) + 1
      indices = pack([(time, time=start, size(series))], &
         (series(start - delay:size(series) - delay) <= threshold) .eqv. &
         (regime == 1))
      allocate(beta(size(lags) + 1), source=0.0_dp)
      if (size(indices) <= size(beta)) return
      allocate(design(size(indices), size(beta)), response(size(indices)))
      design(:, 1) = 1.0_dp
      do row = 1, size(indices)
         response(row) = series(indices(row))
         do lag = 1, size(lags)
            design(row, lag + 1) = series(indices(row) - lags(lag))
         end do
      end do
      cross = matmul(transpose(design), design)
      allocate(inverse(size(beta), size(beta)))
      call symmetric_pseudoinverse(cross, inverse, info)
      if (info == 0) beta = matmul(inverse, matmul(transpose(design), response))
   end subroutine regime_ar_start

   pure function pack_model(model) result(parameters)
      !! Pack a TARMA model in lower AR, upper AR, lower MA, upper MA order.
      type(tseriestarma_model_t), intent(in) :: model !! TARMA model.
      real(dp), allocatable :: parameters(:)
      integer :: first, last, count

      count = size(model%phi_lower) + size(model%phi_upper) + &
         size(model%theta_lower) + size(model%theta_upper)
      allocate(parameters(count))
      first = 1
      last = size(model%phi_lower)
      parameters(first:last) = model%phi_lower
      first = last + 1
      last = first + size(model%phi_upper) - 1
      parameters(first:last) = model%phi_upper
      first = last + 1
      last = first + size(model%theta_lower) - 1
      parameters(first:last) = model%theta_lower
      first = last + 1
      last = first + size(model%theta_upper) - 1
      parameters(first:last) = model%theta_upper
   end function pack_model

   pure function unpack_like(parameters, template) result(model)
      !! Unpack parameters using a template model's sparse lag specification.
      real(dp), intent(in) :: parameters(:) !! Packed TARMA parameters.
      type(tseriestarma_model_t), intent(in) :: template !! Model supplying lag and threshold metadata.
      type(tseriestarma_model_t) :: model

      model = unpack_model(parameters, template%ar_lags_lower, &
         template%ar_lags_upper, template%ma_lags_lower, &
         template%ma_lags_upper, template%delay, template%threshold)
      model%innovation_sd_lower = template%innovation_sd_lower
      model%innovation_sd_upper = template%innovation_sd_upper
   end function unpack_like

   pure function unpack_model(parameters, ar_lags_lower, ar_lags_upper, &
      ma_lags_lower, ma_lags_upper, delay, threshold) result(model)
      !! Convert the package's packed parameter order into a TARMA model.
      real(dp), intent(in) :: parameters(:) !! Packed lower AR, upper AR, lower MA, upper MA values.
      integer, intent(in) :: ar_lags_lower(:) !! Lower AR lag indices.
      integer, intent(in) :: ar_lags_upper(:) !! Upper AR lag indices.
      integer, intent(in) :: ma_lags_lower(:) !! Lower MA lag indices.
      integer, intent(in) :: ma_lags_upper(:) !! Upper MA lag indices.
      integer, intent(in) :: delay !! Threshold delay.
      real(dp), intent(in) :: threshold !! Threshold value.
      type(tseriestarma_model_t) :: model
      integer :: first, last

      first = 1
      last = 1 + size(ar_lags_lower)
      model%phi_lower = parameters(first:last)
      first = last + 1
      last = first + size(ar_lags_upper)
      model%phi_upper = parameters(first:last)
      first = last + 1
      last = first + size(ma_lags_lower) - 1
      model%theta_lower = parameters(first:last)
      first = last + 1
      last = first + size(ma_lags_upper) - 1
      model%theta_upper = parameters(first:last)
      model%ar_lags_lower = ar_lags_lower
      model%ar_lags_upper = ar_lags_upper
      model%ma_lags_lower = ma_lags_lower
      model%ma_lags_upper = ma_lags_upper
      model%delay = delay
      model%threshold = threshold
      if (.not. valid_model(model)) model%info = 1
   end function unpack_model

   pure logical function valid_model(model) result(valid)
      !! Check TARMA dimensions, lags, scales, and finite parameters.
      type(tseriestarma_model_t), intent(in) :: model !! TARMA model.

      valid = model%info == 0 .and. allocated(model%phi_lower) .and. &
         allocated(model%phi_upper) .and. allocated(model%theta_lower) .and. &
         allocated(model%theta_upper) .and. allocated(model%ar_lags_lower) .and. &
         allocated(model%ar_lags_upper) .and. allocated(model%ma_lags_lower) .and. &
         allocated(model%ma_lags_upper)
      if (.not. valid) return
      valid = size(model%phi_lower) == size(model%ar_lags_lower) + 1 .and. &
         size(model%phi_upper) == size(model%ar_lags_upper) + 1 .and. &
         size(model%theta_lower) == size(model%ma_lags_lower) .and. &
         size(model%theta_upper) == size(model%ma_lags_upper) .and. &
         model%delay >= 1 .and. model%innovation_sd_lower > 0.0_dp .and. &
         model%innovation_sd_upper > 0.0_dp .and. &
         all(model%ar_lags_lower >= 1) .and. all(model%ar_lags_upper >= 1) .and. &
         all(model%ma_lags_lower >= 1) .and. all(model%ma_lags_upper >= 1) .and. &
         all(ieee_is_finite(model%phi_lower)) .and. &
         all(ieee_is_finite(model%phi_upper)) .and. &
         all(ieee_is_finite(model%theta_lower)) .and. &
         all(ieee_is_finite(model%theta_upper)) .and. &
         ieee_is_finite(model%threshold)
   end function valid_model

   pure integer function model_start(model) result(start)
      !! Return the largest AR, MA, or threshold lag required by a model.
      type(tseriestarma_model_t), intent(in) :: model !! TARMA model.

      start = maximum_lag(model%ar_lags_lower, model%ar_lags_upper, &
         model%ma_lags_lower, model%ma_lags_upper, model%delay)
   end function model_start

   pure integer function maximum_lag(ar_lower, ar_upper, ma_lower, ma_upper, &
      delay) result(value)
      !! Return the largest lag across four sparse lag sets and the delay.
      integer, intent(in) :: ar_lower(:) !! Lower AR lags.
      integer, intent(in) :: ar_upper(:) !! Upper AR lags.
      integer, intent(in) :: ma_lower(:) !! Lower MA lags.
      integer, intent(in) :: ma_upper(:) !! Upper MA lags.
      integer, intent(in) :: delay !! Threshold delay.

      value = max(delay, max_or_zero(ar_lower), max_or_zero(ar_upper), &
         max_or_zero(ma_lower), max_or_zero(ma_upper))
   end function maximum_lag

   pure integer function max_or_zero(values) result(value)
      !! Return an integer-vector maximum or zero for an empty vector.
      integer, intent(in) :: values(:) !! Integer values.

      value = 0
      if (size(values) > 0) value = maxval(values)
   end function max_or_zero

   pure logical function sets_overlap(first, second) result(overlap)
      !! Report whether two integer lag sets share at least one value.
      integer, intent(in) :: first(:) !! First lag set.
      integer, intent(in) :: second(:) !! Second lag set.
      integer :: index

      overlap = .false.
      do index = 1, size(first)
         if (any(second == first(index))) then
            overlap = .true.
            return
         end if
      end do
   end function sets_overlap

   pure function recursive_filter_matrix(input, feedback) result(filtered)
      !! Apply one recursive filter independently to every input column.
      real(dp), intent(in) :: input(:, :) !! Matrix of unfiltered derivative inputs.
      real(dp), intent(in) :: feedback(:) !! Recursive feedback coefficients.
      real(dp) :: filtered(size(input, 1), size(input, 2))
      integer :: lag, time

      filtered = input
      do time = 1, size(input, 1)
         do lag = 1, min(size(feedback), time - 1)
            filtered(time, :) = filtered(time, :) + &
               feedback(lag)*filtered(time - lag, :)
         end do
      end do
   end function recursive_filter_matrix

end module tseriestarma_mod
