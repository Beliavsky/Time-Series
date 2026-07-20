! SPDX-License-Identifier: GPL-2.0-or-later
! SPDX-FileComment: Algorithms translated from the R mixAR package.
module mixar_mod
   !! Gaussian mixture autoregression likelihood, simulation, and EM estimation.
   use kind_mod, only: dp
   use linalg_mod, only: solve_matrix, general_eigenvalues, kronecker_product, &
      symmetric_pseudoinverse, inverse_logdet, symmetric_eigen
   use random_mod, only: random_uniform, random_standard_normal, &
      random_standard_student, random_standard_normal_matrix, &
      multivariate_normal_from_standard
   use special_functions_mod, only: regularized_beta, &
      multivariate_normal_log_density, regularized_gamma_q
   use stats_mod, only: sorted, quantile, variance
   use stats_mod, only: normal_quantile
   use time_series_stats_mod, only: acf_values
   use time_series_diagnostics_mod, only: weighted_box_test_t, &
      weighted_box_test, box_test_ljung_box, residual_raw, residual_squared, &
      multivariate_white_noise_test_t, multivariate_white_noise_test
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   private

   real(dp), parameter :: log_two_pi = 1.83787706640934548356_dp
   integer, parameter, public :: mixar_order_prior_flat = 0
   integer, parameter, public :: mixar_order_prior_ratio = 1
   integer, parameter, public :: mixar_order_prior_poisson = 2

   type, public :: mixar_model_t
      !! Gaussian MAR parameters with padded columns for ragged AR coefficients.
      real(dp), allocatable :: probability(:)
      integer, allocatable :: order(:)
      real(dp), allocatable :: shift(:)
      real(dp), allocatable :: scale(:)
      real(dp), allocatable :: ar(:, :)
      real(dp), allocatable :: degrees_of_freedom(:)
      integer :: info = 0
   end type mixar_model_t

   type, public :: mixar_filter_t
      !! Conditional component locations, probabilities, and mixture diagnostics.
      real(dp), allocatable :: location(:, :)
      real(dp), allocatable :: responsibility(:, :)
      real(dp), allocatable :: fitted(:)
      real(dp), allocatable :: variance(:)
      real(dp), allocatable :: residual(:)
      real(dp) :: log_likelihood = -huge(1.0_dp)
      integer :: start = 0
      integer :: info = 0
   end type mixar_filter_t

   type, public :: mixar_fit_t
      !! Gaussian MAR EM estimate and convergence diagnostics.
      type(mixar_model_t) :: model
      type(mixar_filter_t) :: filter
      integer :: iterations = 0
      logical :: converged = .false.
      real(dp) :: log_likelihood = -huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      real(dp) :: bic = huge(1.0_dp)
      integer :: info = 0
   end type mixar_fit_t

   type, public :: mixar_simulation_t
      !! Simulated Gaussian MAR observations and component regimes.
      real(dp), allocatable :: series(:)
      integer, allocatable :: regime(:)
      type(mixar_model_t) :: model
      integer :: burnin = 0
      integer :: info = 0
   end type mixar_simulation_t

   type, public :: mixar_predictive_distribution_t
      !! Exact Gaussian mixture distribution at one forecast horizon.
      real(dp), allocatable :: probability(:)
      real(dp), allocatable :: location(:)
      real(dp), allocatable :: scale(:)
      integer, allocatable :: regime_path(:, :)
      integer :: horizon = 0
      integer :: info = 0
   end type mixar_predictive_distribution_t

   type, public :: mixar_exact_forecast_t
      !! Exact path-mixture forecasts and pointwise distribution summaries.
      type(mixar_predictive_distribution_t), allocatable :: distribution(:)
      real(dp), allocatable :: mean(:)
      real(dp), allocatable :: variance(:)
      real(dp), allocatable :: lower(:)
      real(dp), allocatable :: upper(:)
      real(dp) :: level = 0.95_dp
      integer :: info = 0
   end type mixar_exact_forecast_t

   type, public :: mixar_forecast_sample_t
      !! Simulated forecast paths and pointwise empirical summaries.
      real(dp), allocatable :: paths(:, :)
      real(dp), allocatable :: mean(:)
      real(dp), allocatable :: variance(:)
      real(dp), allocatable :: lower(:)
      real(dp), allocatable :: upper(:)
      real(dp) :: level = 0.95_dp
      integer :: info = 0
   end type mixar_forecast_sample_t

   type, public :: mixar_regression_model_t
      !! Mixture autoregression with shared or component-specific regression effects.
      type(mixar_model_t) :: mar
      real(dp), allocatable :: coefficient(:, :)
      logical :: component_specific = .false.
      integer :: info = 0
   end type mixar_regression_model_t

   type, public :: mixar_regression_fit_t
      !! Fitted mixture autoregression with regression covariates.
      type(mixar_regression_model_t) :: model
      type(mixar_filter_t) :: filter
      integer :: iterations = 0
      logical :: converged = .false.
      real(dp) :: log_likelihood = -huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      real(dp) :: bic = huge(1.0_dp)
      integer :: info = 0
   end type mixar_regression_fit_t

   type, public :: mixar_regression_simulation_t
      !! Simulated regression-MAR observations, regimes, and innovations.
      real(dp), allocatable :: series(:)
      integer, allocatable :: regime(:)
      real(dp), allocatable :: innovation(:)
      type(mixar_regression_model_t) :: model
      integer :: info = 0
   end type mixar_regression_simulation_t

   type, public :: mixar_seasonal_model_t
      !! Additive ordinary- and seasonal-lag mixture autoregression.
      type(mixar_model_t) :: mar
      integer :: seasonal_period = 0
      integer, allocatable :: seasonal_order(:)
      real(dp), allocatable :: seasonal_ar(:, :)
      integer :: info = 0
   end type mixar_seasonal_model_t

   type, public :: mixar_seasonal_fit_t
      !! Fitted seasonal MAR parameters and convergence diagnostics.
      type(mixar_seasonal_model_t) :: model
      type(mixar_filter_t) :: filter
      integer :: iterations = 0
      logical :: converged = .false.
      real(dp) :: log_likelihood = -huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      real(dp) :: bic = huge(1.0_dp)
      integer :: info = 0
   end type mixar_seasonal_fit_t

   type, public :: mixar_seasonal_simulation_t
      !! Simulated seasonal MAR observations and component regimes.
      real(dp), allocatable :: series(:)
      integer, allocatable :: regime(:)
      type(mixar_seasonal_model_t) :: model
      integer :: burnin = 0
      integer :: info = 0
   end type mixar_seasonal_simulation_t

   type, public :: mixar_inference_t
      !! Observed information, covariance, standard errors, and confidence limits.
      character(len=48), allocatable :: parameter(:)
      real(dp), allocatable :: estimate(:)
      real(dp), allocatable :: standard_error(:)
      real(dp), allocatable :: lower(:)
      real(dp), allocatable :: upper(:)
      real(dp), allocatable :: observed_information(:, :)
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: probability_standard_error(:)
      real(dp) :: level = 0.95_dp
      integer :: info = 0
   end type mixar_inference_t

   type, public :: mixar_diagnostics_t
      !! Residual sequences, serial-correlation tests, and classification diagnostics.
      real(dp), allocatable :: residual(:)
      real(dp), allocatable :: standardized_residual(:)
      real(dp), allocatable :: uniform_residual(:)
      real(dp), allocatable :: quantile_residual(:)
      real(dp), allocatable :: classified_residual(:)
      real(dp), allocatable :: residual_acf(:)
      real(dp), allocatable :: squared_residual_acf(:)
      real(dp), allocatable :: posterior_probability(:, :)
      real(dp), allocatable :: maximum_posterior(:)
      integer, allocatable :: classification(:)
      type(weighted_box_test_t) :: residual_test
      type(weighted_box_test_t) :: squared_residual_test
      real(dp) :: mean_entropy = 0.0_dp
      real(dp) :: bic = huge(1.0_dp)
      integer :: info = 0
   end type mixar_diagnostics_t

   type, public :: mixar_bic_selection_t
      !! BIC scores and selected candidate model index.
      real(dp), allocatable :: bic(:)
      integer :: selected = 0
      integer :: info = 0
   end type mixar_bic_selection_t

   type, public :: mixar_var_model_t
      !! Gaussian mixture VAR parameters with padded ragged lag matrices.
      real(dp), allocatable :: probability(:)
      integer, allocatable :: order(:)
      real(dp), allocatable :: intercept(:, :)
      real(dp), allocatable :: ar(:, :, :, :)
      real(dp), allocatable :: innovation_covariance(:, :, :)
      integer :: info = 0
   end type mixar_var_model_t

   type, public :: mixar_var_filter_t
      !! Mixture VAR conditional moments and posterior classifications.
      real(dp), allocatable :: location(:, :, :)
      real(dp), allocatable :: responsibility(:, :)
      real(dp), allocatable :: fitted(:, :)
      real(dp), allocatable :: residual(:, :)
      real(dp), allocatable :: conditional_covariance(:, :, :)
      real(dp) :: log_likelihood = -huge(1.0_dp)
      integer :: start = 0
      integer :: info = 0
   end type mixar_var_filter_t

   type, public :: mixar_var_fit_t
      !! Gaussian mixture VAR EM estimate and convergence diagnostics.
      type(mixar_var_model_t) :: model
      type(mixar_var_filter_t) :: filter
      integer :: iterations = 0
      logical :: converged = .false.
      real(dp) :: log_likelihood = -huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      real(dp) :: bic = huge(1.0_dp)
      integer :: info = 0
   end type mixar_var_fit_t

   type, public :: mixar_var_simulation_t
      !! Simulated mixture VAR observations, regimes, and parameters.
      real(dp), allocatable :: series(:, :)
      integer, allocatable :: regime(:)
      type(mixar_var_model_t) :: model
      integer :: burnin = 0
      integer :: info = 0
   end type mixar_var_simulation_t

   type, public :: mixar_var_forecast_t
      !! Simulated mixture VAR forecast paths and pointwise summaries.
      real(dp), allocatable :: paths(:, :, :)
      real(dp), allocatable :: mean(:, :)
      real(dp), allocatable :: covariance(:, :, :)
      real(dp), allocatable :: lower(:, :)
      real(dp), allocatable :: upper(:, :)
      real(dp) :: level = 0.95_dp
      integer :: info = 0
   end type mixar_var_forecast_t

   type, public :: mixar_var_diagnostics_t
      !! Multivariate residual, whiteness, and classification diagnostics.
      real(dp), allocatable :: residual(:, :)
      real(dp), allocatable :: standardized_residual(:, :)
      real(dp), allocatable :: posterior_probability(:, :)
      real(dp), allocatable :: maximum_posterior(:)
      integer, allocatable :: classification(:)
      type(multivariate_white_noise_test_t) :: white_noise_test
      real(dp) :: mean_entropy = 0.0_dp
      integer :: info = 0
   end type mixar_var_diagnostics_t

   type, public :: mixar_bayesian_random_t
      !! Independent variates used by the reproducible Bayesian MAR sampler.
      real(dp), allocatable :: allocation_uniform(:, :)
      real(dp), allocatable :: probability_uniform(:, :)
      real(dp), allocatable :: mean_normal(:, :)
      real(dp), allocatable :: lambda_uniform(:)
      real(dp), allocatable :: precision_uniform(:, :)
      real(dp), allocatable :: ar_normal(:, :, :)
      real(dp), allocatable :: acceptance_uniform(:, :)
   end type mixar_bayesian_random_t

   type, public :: mixar_bayesian_draws_t
      !! Retained Bayesian MAR parameter draws and posterior summaries.
      real(dp), allocatable :: probability(:, :)
      real(dp), allocatable :: scale(:, :)
      real(dp), allocatable :: precision(:, :)
      real(dp), allocatable :: shift(:, :)
      real(dp), allocatable :: component_mean(:, :)
      real(dp), allocatable :: ar(:, :, :)
      real(dp), allocatable :: lambda(:)
      integer, allocatable :: last_allocation(:)
      real(dp), allocatable :: acceptance_rate(:)
      integer, allocatable :: order(:)
      type(mixar_model_t) :: posterior_mean_model
      integer :: burnin = 0
      logical :: fixed_shift = .false.
      integer :: info = 0
   end type mixar_bayesian_draws_t

   type, public :: mixar_marginal_likelihood_t
      !! Chib-style marginal-likelihood components at a representative posterior point.
      real(dp) :: log_marginal_likelihood = -huge(1.0_dp)
      real(dp) :: log_likelihood = -huge(1.0_dp)
      real(dp) :: log_prior_density = -huge(1.0_dp)
      real(dp) :: log_posterior_ordinate = -huge(1.0_dp)
      type(mixar_model_t) :: representative_model
      integer :: info = 0
   end type mixar_marginal_likelihood_t

   type, public :: mixar_order_random_t
      !! Independent variates used by reproducible mixture-order RJMCMC.
      type(mixar_bayesian_random_t), allocatable :: parameter_random(:)
      real(dp), allocatable :: component_uniform(:)
      real(dp), allocatable :: direction_uniform(:)
      real(dp), allocatable :: coefficient_uniform(:)
      real(dp), allocatable :: acceptance_uniform(:)
   end type mixar_order_random_t

   type, public :: mixar_order_selection_t
      !! Mixture component-order RJMCMC trace and posterior model frequencies.
      integer, allocatable :: order_trace(:, :)
      real(dp), allocatable :: probability_trace(:, :)
      real(dp), allocatable :: shift_trace(:, :)
      real(dp), allocatable :: scale_trace(:, :)
      real(dp), allocatable :: ar_trace(:, :, :)
      integer, allocatable :: proposed_component(:)
      integer, allocatable :: proposed_order(:)
      logical, allocatable :: accepted(:)
      integer, allocatable :: unique_order(:, :)
      real(dp), allocatable :: posterior_probability(:)
      integer, allocatable :: modal_order(:)
      type(mixar_model_t) :: modal_model
      real(dp) :: acceptance_rate = 0.0_dp
      integer :: maximum_order = 0
      integer :: prior_method = mixar_order_prior_flat
      real(dp) :: prior_parameter = 1.0_dp
      logical :: fixed_shift = .false.
      integer :: info = 0
   end type mixar_order_selection_t

   type, public :: mixar_initialization_t
      !! Data-driven Gaussian MAR starting model and component residuals.
      type(mixar_model_t) :: model
      real(dp), allocatable :: residual(:, :)
      real(dp), allocatable :: responsibility(:, :)
      integer, allocatable :: subsample_index(:, :)
      integer :: stability_contractions = 0
      integer :: info = 0
   end type mixar_initialization_t

   type, public :: mixar_multistart_t
      !! Gaussian MAR initializations, EM fits, and best-likelihood start.
      type(mixar_initialization_t), allocatable :: initialization(:)
      type(mixar_fit_t), allocatable :: fit(:)
      real(dp), allocatable :: initial_log_likelihood(:)
      integer :: best = 0
      integer :: successful = 0
      integer :: info = 0
   end type mixar_multistart_t

   type, public :: mixar_moment_t
      !! Scalar analytical moment with existence and numerical status.
      real(dp) :: value = 0.0_dp
      integer :: order = 0
      logical :: exists = .false.
      integer :: info = 0
   end type mixar_moment_t

   type, public :: mixar_stationary_moments_t
      !! Stable MAR companion-state mean, covariance, and marginal moments.
      real(dp), allocatable :: state_mean(:)
      real(dp), allocatable :: state_covariance(:, :)
      real(dp) :: mean = 0.0_dp
      real(dp) :: variance = 0.0_dp
      integer :: info = 0
   end type mixar_stationary_moments_t

   public :: mixar_model, mixar_component_locations, mixar_conditional_density
   public :: mixar_conditional_cdf, mixar_filter, mixar_log_likelihood
   public :: mixar_simulate_from_draws, mixar_simulate, mixar_is_stable
   public :: mixar_fit, mixar_general_fit, mixar_standard_density
   public :: mixar_standard_cdf, mixar_exact_forecast
   public :: mixar_predictive_density, mixar_predictive_cdf
   public :: mixar_predictive_quantile, mixar_forecast_from_draws
   public :: mixar_forecast, display
   public :: mixar_regression_model, mixar_regression_filter
   public :: mixar_regression_log_likelihood, mixar_regression_fit
   public :: mixar_regression_simulate_from_draws, mixar_regression_simulate
   public :: mixar_regression_exact_forecast
   public :: mixar_regression_forecast_from_draws, mixar_regression_forecast
   public :: mixar_seasonal_model, mixar_seasonal_expanded_model
   public :: mixar_seasonal_filter, mixar_seasonal_log_likelihood
   public :: mixar_seasonal_fit, mixar_seasonal_is_stable
   public :: mixar_seasonal_simulate_from_draws, mixar_seasonal_simulate
   public :: mixar_seasonal_exact_forecast
   public :: mixar_seasonal_forecast_from_draws, mixar_seasonal_forecast
   public :: mixar_observed_inference, mixar_seasonal_observed_inference
   public :: mixar_regression_observed_inference
   public :: mixar_diagnose, mixar_bic, mixar_select_bic
   public :: mixar_var_model, mixar_var_component_locations, mixar_var_filter
   public :: mixar_var_log_likelihood, mixar_var_fit, mixar_var_is_stable
   public :: mixar_var_simulate_from_draws, mixar_var_simulate
   public :: mixar_var_forecast_from_draws, mixar_var_forecast
   public :: mixar_var_diagnose
   public :: mixar_bayesian_random, mixar_bayesian_sample_from_random
   public :: mixar_bayesian_sample, mixar_bayesian_relabel
   public :: mixar_marginal_likelihood
   public :: mixar_order_birth_death_probability, mixar_order_random
   public :: mixar_order_select_from_random, mixar_order_select
   public :: mixar_initialize_from_indices, mixar_random_subsample_indices
   public :: mixar_random_initialize, mixar_multistart_fit
   public :: mixar_random_multistart_fit
   public :: mixar_standard_moment, mixar_standard_absolute_moment
   public :: mixar_innovation_moment, mixar_conditional_moment
   public :: mixar_conditional_central_moment, mixar_conditional_kurtosis
   public :: mixar_conditional_excess_kurtosis, mixar_stationary_moments

   interface display
      module procedure display_mixar_model
      module procedure display_mixar_fit
      module procedure display_mixar_simulation
      module procedure display_mixar_exact_forecast
      module procedure display_mixar_forecast_sample
      module procedure display_mixar_regression_model
      module procedure display_mixar_regression_fit
      module procedure display_mixar_regression_simulation
      module procedure display_mixar_seasonal_model
      module procedure display_mixar_seasonal_fit
      module procedure display_mixar_seasonal_simulation
      module procedure display_mixar_inference
      module procedure display_mixar_diagnostics
      module procedure display_mixar_var_model
      module procedure display_mixar_var_fit
      module procedure display_mixar_var_simulation
      module procedure display_mixar_var_forecast
      module procedure display_mixar_var_diagnostics
      module procedure display_mixar_bayesian_draws
      module procedure display_mixar_marginal_likelihood
      module procedure display_mixar_order_selection
      module procedure display_mixar_initialization
      module procedure display_mixar_multistart
      module procedure display_mixar_moment
      module procedure display_mixar_stationary_moments
   end interface display

   abstract interface
      pure function scalar_objective_t(parameters) result(value)
         !! Evaluate a scalar likelihood objective for numerical differentiation.
         import dp
         real(dp), intent(in) :: parameters(:) !! Candidate parameter vector.
         real(dp) :: value
      end function scalar_objective_t
   end interface

contains

   pure function mixar_model(probability, order, shift, scale, ar, &
      degrees_of_freedom) result(model)
      !! Construct and validate a Gaussian mixture autoregression model.
      real(dp), intent(in) :: probability(:) !! Component mixing probabilities.
      integer, intent(in) :: order(:) !! Autoregressive order of each component.
      real(dp), intent(in) :: shift(:) !! Component regression intercepts.
      real(dp), intent(in) :: scale(:) !! Component innovation standard deviations.
      real(dp), intent(in) :: ar(:, :) !! Padded lag coefficients by lag and component.
      real(dp), intent(in), optional :: degrees_of_freedom(:) !! Zero for Gaussian or Student-t degrees above two.
      type(mixar_model_t) :: model

      model%probability = probability
      model%order = order
      model%shift = shift
      model%scale = scale
      model%ar = ar
      allocate(model%degrees_of_freedom(size(probability)), source=0.0_dp)
      if (present(degrees_of_freedom)) model%degrees_of_freedom = degrees_of_freedom
      if (.not. valid_model(model)) model%info = 1
   end function mixar_model

   pure elemental real(dp) function mixar_standard_density(value, &
      degrees_of_freedom) result(density)
      !! Evaluate a standard normal or unit-variance Student-t density.
      real(dp), intent(in) :: value !! Standardized innovation value.
      real(dp), intent(in) :: degrees_of_freedom !! Zero for Gaussian or Student-t degrees above two.

      if (degrees_of_freedom == 0.0_dp) then
         density = exp(-0.5_dp*(log_two_pi + value**2))
      else if (degrees_of_freedom > 2.0_dp) then
         density = exp(log_gamma(0.5_dp*(degrees_of_freedom + 1.0_dp)) - &
            log_gamma(0.5_dp*degrees_of_freedom) - &
            0.5_dp*log(acos(-1.0_dp)*(degrees_of_freedom - 2.0_dp)) - &
            0.5_dp*(degrees_of_freedom + 1.0_dp)* &
            log(1.0_dp + value**2/(degrees_of_freedom - 2.0_dp)))
      else
         density = 0.0_dp
      end if
   end function mixar_standard_density

   pure elemental real(dp) function mixar_standard_cdf(value, &
      degrees_of_freedom) result(probability)
      !! Evaluate a standard normal or unit-variance Student-t distribution function.
      real(dp), intent(in) :: value !! Standardized innovation value.
      real(dp), intent(in) :: degrees_of_freedom !! Zero for Gaussian or Student-t degrees above two.
      real(dp) :: beta_value

      if (degrees_of_freedom == 0.0_dp) then
         probability = 0.5_dp*erfc(-value/sqrt(2.0_dp))
      else if (degrees_of_freedom > 2.0_dp) then
         beta_value = regularized_beta((degrees_of_freedom - 2.0_dp)/ &
            (degrees_of_freedom - 2.0_dp + value**2), &
            0.5_dp*degrees_of_freedom, 0.5_dp)
         if (value >= 0.0_dp) then
            probability = 1.0_dp - 0.5_dp*beta_value
         else
            probability = 0.5_dp*beta_value
         end if
      else
         probability = 0.0_dp
      end if
   end function mixar_standard_cdf

   pure function mixar_component_locations(model, series, start) result(location)
      !! Evaluate every component's conditional location over a series.
      type(mixar_model_t), intent(in) :: model !! Gaussian MAR model.
      real(dp), intent(in) :: series(:) !! Observed series including required lags.
      integer, intent(in), optional :: start !! First observation to evaluate.
      real(dp), allocatable :: location(:, :)
      integer :: first, time, component, lag

      first = maxval(model%order) + 1
      if (present(start)) first = start
      if (.not. valid_model(model) .or. first < maxval(model%order) + 1 .or. &
         first > size(series) + 1) then
         allocate(location(0, 0))
         return
      end if
      allocate(location(max(0, size(series) - first + 1), size(model%order)))
      do component = 1, size(model%order)
         location(:, component) = model%shift(component)
         do lag = 1, model%order(component)
            do time = first, size(series)
               location(time - first + 1, component) = &
                  location(time - first + 1, component) + &
                  model%ar(lag, component)*series(time - lag)
            end do
         end do
      end do
   end function mixar_component_locations

   pure function mixar_conditional_density(model, value, history) result(density)
      !! Evaluate the one-step Gaussian mixture conditional density.
      type(mixar_model_t), intent(in) :: model !! Gaussian MAR model.
      real(dp), intent(in) :: value !! Candidate next observation.
      real(dp), intent(in) :: history(:) !! History ending with the most recent observation.
      real(dp) :: density
      real(dp) :: location, standardized
      integer :: component, lag

      density = 0.0_dp
      if (.not. valid_model(model) .or. size(history) < maxval(model%order)) return
      do component = 1, size(model%order)
         location = model%shift(component)
         do lag = 1, model%order(component)
            location = location + model%ar(lag, component)* &
               history(size(history) - lag + 1)
         end do
         standardized = (value - location)/model%scale(component)
         density = density + model%probability(component)* &
            mixar_standard_density(standardized, &
            model%degrees_of_freedom(component))/model%scale(component)
      end do
   end function mixar_conditional_density

   pure function mixar_conditional_cdf(model, value, history) result(probability)
      !! Evaluate the one-step Gaussian mixture conditional distribution function.
      type(mixar_model_t), intent(in) :: model !! Gaussian MAR model.
      real(dp), intent(in) :: value !! Candidate next observation.
      real(dp), intent(in) :: history(:) !! History ending with the most recent observation.
      real(dp) :: probability
      real(dp) :: location
      integer :: component, lag

      probability = 0.0_dp
      if (.not. valid_model(model) .or. size(history) < maxval(model%order)) return
      do component = 1, size(model%order)
         location = model%shift(component)
         do lag = 1, model%order(component)
            location = location + model%ar(lag, component)* &
               history(size(history) - lag + 1)
         end do
         probability = probability + model%probability(component)* &
            mixar_standard_cdf((value - location)/model%scale(component), &
            model%degrees_of_freedom(component))
      end do
   end function mixar_conditional_cdf

   pure function mixar_filter(model, series) result(out)
      !! Filter a series and compute posterior component responsibilities.
      type(mixar_model_t), intent(in) :: model !! Gaussian MAR model.
      real(dp), intent(in) :: series(:) !! Observed time series.
      type(mixar_filter_t) :: out
      real(dp), allocatable :: log_kernel(:)
      real(dp) :: maximum, denominator
      integer :: observations, components, row, component

      if (.not. valid_model(model) .or. size(series) <= maxval(model%order) .or. &
         .not. all(ieee_is_finite(series))) then
         out%info = 1
         return
      end if
      out%start = maxval(model%order) + 1
      observations = size(series) - out%start + 1
      components = size(model%order)
      out%location = mixar_component_locations(model, series, out%start)
      allocate(out%responsibility(observations, components), out%fitted(observations))
      allocate(out%variance(observations), out%residual(observations), log_kernel(components))
      out%log_likelihood = 0.0_dp
      do row = 1, observations
         do component = 1, components
            log_kernel(component) = log(max(model%probability(component), tiny(1.0_dp))) - &
               log(model%scale(component)) + log(max(tiny(1.0_dp), &
               mixar_standard_density((series(out%start + row - 1) - &
               out%location(row, component))/model%scale(component), &
               model%degrees_of_freedom(component))))
         end do
         maximum = maxval(log_kernel)
         denominator = sum(exp(log_kernel - maximum))
         out%responsibility(row, :) = exp(log_kernel - maximum)/denominator
         out%log_likelihood = out%log_likelihood + maximum + log(denominator)
         out%fitted(row) = dot_product(model%probability, out%location(row, :))
         out%variance(row) = dot_product(model%probability, model%scale**2 + &
            out%location(row, :)**2) - out%fitted(row)**2
         out%residual(row) = series(out%start + row - 1) - out%fitted(row)
      end do
   end function mixar_filter

   pure real(dp) function mixar_log_likelihood(model, series) result(log_likelihood)
      !! Evaluate the conditional Gaussian mixture log likelihood.
      type(mixar_model_t), intent(in) :: model !! Gaussian MAR model.
      real(dp), intent(in) :: series(:) !! Observed time series.
      type(mixar_filter_t) :: filtered

      filtered = mixar_filter(model, series)
      log_likelihood = filtered%log_likelihood
   end function mixar_log_likelihood

   pure function mixar_simulate_from_draws(model, observations, initial, &
      uniforms, normals, burnin) result(out)
      !! Simulate a Gaussian MAR path from supplied regime and normal draws.
      type(mixar_model_t), intent(in) :: model !! Gaussian MAR model.
      integer, intent(in) :: observations !! Number of retained observations.
      real(dp), intent(in) :: initial(:) !! Initial values ending at the simulation origin.
      real(dp), intent(in) :: uniforms(:) !! Uniform component-selection draws.
      real(dp), intent(in) :: normals(:) !! Standard-normal innovation draws.
      integer, intent(in), optional :: burnin !! Number of discarded observations.
      type(mixar_simulation_t) :: out
      real(dp), allocatable :: work(:)
      real(dp) :: cumulative
      integer :: discarded, total, p, time, component, lag

      discarded = 100
      if (present(burnin)) discarded = burnin
      p = maxval(model%order)
      total = observations + discarded
      if (.not. valid_model(model) .or. observations < 1 .or. discarded < 0 .or. &
         size(initial) < p .or. size(uniforms) < total .or. size(normals) < total .or. &
         any(uniforms(1:total) < 0.0_dp) .or. any(uniforms(1:total) > 1.0_dp)) then
         out%info = 1
         return
      end if
      allocate(work(p + total), out%regime(observations), out%series(observations))
      work(1:p) = initial(size(initial) - p + 1:size(initial))
      do time = 1, total
         cumulative = 0.0_dp
         component = size(model%probability)
         do lag = 1, size(model%probability)
            cumulative = cumulative + model%probability(lag)
            if (uniforms(time) <= cumulative) then
               component = lag
               exit
            end if
         end do
         work(p + time) = model%shift(component) + model%scale(component)*normals(time)
         do lag = 1, model%order(component)
            work(p + time) = work(p + time) + model%ar(lag, component)*work(p + time - lag)
         end do
         if (time > discarded) out%regime(time - discarded) = component
      end do
      out%series = work(p + discarded + 1:p + total)
      out%model = model
      out%burnin = discarded
   end function mixar_simulate_from_draws

   function mixar_simulate(model, observations, initial, burnin) result(out)
      !! Simulate a Gaussian MAR path using the library's shared random generator.
      type(mixar_model_t), intent(in) :: model !! Gaussian MAR model.
      integer, intent(in) :: observations !! Number of retained observations.
      real(dp), intent(in) :: initial(:) !! Initial values ending at the simulation origin.
      integer, intent(in), optional :: burnin !! Number of discarded observations.
      type(mixar_simulation_t) :: out
      real(dp), allocatable :: uniforms(:), normals(:)
      integer :: discarded, total, i

      discarded = 100
      if (present(burnin)) discarded = burnin
      total = max(0, observations + discarded)
      allocate(uniforms(total), normals(total))
      do i = 1, total
         uniforms(i) = random_uniform()
         normals(i) = random_standard_normal()
      end do
      if (valid_model(model)) then
         do i = 1, total
            normals(i) = component_random_innovation(model, uniforms(i), normals(i))
         end do
      end if
      out = mixar_simulate_from_draws(model, observations, initial, uniforms, normals, discarded)
   end function mixar_simulate

   function component_random_innovation(model, uniform, normal_draw) result(innovation)
      !! Generate a standardized innovation for the component selected by a uniform draw.
      type(mixar_model_t), intent(in) :: model !! Gaussian or Student-t MAR model.
      real(dp), intent(in) :: uniform !! Regime-selection uniform draw.
      real(dp), intent(in) :: normal_draw !! Already generated standard-normal draw.
      real(dp) :: innovation
      real(dp) :: cumulative
      integer :: component, candidate

      cumulative = 0.0_dp
      component = size(model%probability)
      do candidate = 1, size(model%probability)
         cumulative = cumulative + model%probability(candidate)
         component = candidate
         if (uniform <= cumulative) exit
      end do
      if (model%degrees_of_freedom(component) > 2.0_dp) then
         innovation = random_standard_student(model%degrees_of_freedom(component))
      else
         innovation = normal_draw
      end if
   end function component_random_innovation

   pure logical function mixar_is_stable(model) result(stable)
      !! Test mixAR second-order stability by its Kronecker companion criterion.
      type(mixar_model_t), intent(in) :: model !! Gaussian MAR model.
      real(dp), allocatable :: companion(:, :), aggregate(:, :)
      complex(dp), allocatable :: eigenvalues(:)
      integer :: p, component, i, info

      stable = .false.
      if (.not. valid_model(model)) return
      p = maxval(model%order)
      if (p == 0) then
         stable = .true.
         return
      end if
      allocate(aggregate(p*p, p*p), companion(p, p))
      aggregate = 0.0_dp
      do component = 1, size(model%order)
         companion = 0.0_dp
         if (model%order(component) > 0) &
            companion(1, 1:model%order(component)) = &
            model%ar(1:model%order(component), component)
         do i = 2, p
            companion(i, i - 1) = 1.0_dp
         end do
         aggregate = aggregate + model%probability(component)* &
            kronecker_product(companion, companion)
      end do
      call general_eigenvalues(aggregate, eigenvalues, info)
      if (info == 0) stable = maxval(abs(eigenvalues)) < 1.0_dp
   end function mixar_is_stable

   pure function mixar_fit(series, initial_model, estimate_shift, &
      max_iterations, tolerance, minimum_scale) result(out)
      !! Fit a Gaussian MAR model by mixAR's fixed-point EM algorithm.
      real(dp), intent(in) :: series(:) !! Observed time series.
      type(mixar_model_t), intent(in) :: initial_model !! Initial Gaussian MAR parameters.
      logical, intent(in), optional :: estimate_shift !! Estimate component intercepts.
      integer, intent(in), optional :: max_iterations !! Maximum EM iterations.
      real(dp), intent(in), optional :: tolerance !! Relative likelihood tolerance.
      real(dp), intent(in), optional :: minimum_scale !! Innovation scale floor.
      type(mixar_fit_t) :: out
      type(mixar_filter_t) :: filtered
      type(mixar_model_t) :: current, updated
      real(dp), allocatable :: a(:, :), b(:, :), solution(:, :), residual(:)
      real(dp) :: previous, current_ll, tol, scale_floor, weight_sum
      integer :: limit, iteration, component, row, lag, time, q, info, parameters
      logical :: fit_shift

      limit = 200
      tol = 1.0e-10_dp
      scale_floor = 1.0e-7_dp
      fit_shift = .true.
      if (present(max_iterations)) limit = max_iterations
      if (present(tolerance)) tol = tolerance
      if (present(minimum_scale)) scale_floor = minimum_scale
      if (present(estimate_shift)) fit_shift = estimate_shift
      if (.not. valid_model(initial_model) .or. size(series) <= maxval(initial_model%order) .or. &
         any(initial_model%degrees_of_freedom /= 0.0_dp) .or. &
         limit < 1 .or. tol <= 0.0_dp .or. scale_floor <= 0.0_dp) then
         out%info = 1
         return
      end if
      current = initial_model
      previous = mixar_log_likelihood(current, series)
      do iteration = 1, limit
         filtered = mixar_filter(current, series)
         if (filtered%info /= 0) then
            out%info = 2
            return
         end if
         updated = current
         do component = 1, size(current%order)
            weight_sum = sum(filtered%responsibility(:, component))
            updated%probability(component) = weight_sum/real(size(filtered%responsibility, 1), dp)
            q = current%order(component) + merge(1, 0, fit_shift)
            if (q > 0 .and. weight_sum > tiny(1.0_dp)) then
               allocate(a(q, q), b(q, 1), solution(q, 1))
               a = 0.0_dp
               b = 0.0_dp
               do row = 1, size(filtered%responsibility, 1)
                  time = filtered%start + row - 1
                  call accumulate_normal_equations(series, time, &
                     filtered%responsibility(row, component), current%order(component), &
                     fit_shift, a, b(:, 1))
               end do
               call solve_matrix(a, b, solution, info)
               if (info /= 0) then
                  out%info = 3
                  return
               end if
               if (fit_shift) then
                  updated%shift(component) = solution(1, 1)
                  if (current%order(component) > 0) &
                     updated%ar(1:current%order(component), component) = solution(2:q, 1)
               else if (current%order(component) > 0) then
                  updated%ar(1:current%order(component), component) = solution(:, 1)
               end if
               deallocate(a, b, solution)
            end if
         end do
         updated%probability = updated%probability/sum(updated%probability)
         updated%scale = 0.0_dp
         allocate(residual(size(filtered%responsibility, 1)))
         do component = 1, size(current%order)
            do row = 1, size(residual)
               time = filtered%start + row - 1
               residual(row) = series(time) - updated%shift(component)
               do lag = 1, updated%order(component)
                  residual(row) = residual(row) - updated%ar(lag, component)*series(time - lag)
               end do
            end do
            weight_sum = sum(filtered%responsibility(:, component))
            if (weight_sum > tiny(1.0_dp)) updated%scale(component) = max(scale_floor, &
               sqrt(dot_product(filtered%responsibility(:, component), residual**2)/weight_sum))
         end do
         deallocate(residual)
         current_ll = mixar_log_likelihood(updated, series)
         if (.not. ieee_is_finite(current_ll)) then
            out%info = 4
            return
         end if
         current = updated
         out%iterations = iteration
         if (abs(current_ll - previous) <= tol*(1.0_dp + abs(previous))) then
            out%converged = .true.
            exit
         end if
         previous = current_ll
      end do
      out%model = current
      out%filter = mixar_filter(current, series)
      out%log_likelihood = out%filter%log_likelihood
      parameters = size(current%probability) - 1 + size(current%scale) + &
         sum(current%order) + merge(size(current%shift), 0, fit_shift)
      out%aic = -2.0_dp*out%log_likelihood + 2.0_dp*real(parameters, dp)
      out%bic = -2.0_dp*out%log_likelihood + log(real(size(out%filter%residual), dp))* &
         real(parameters, dp)
   end function mixar_fit

   pure function mixar_general_fit(series, initial_model, estimate_shift, &
      max_iterations, tolerance, minimum_scale, inner_iterations, &
      estimate_degrees) result(out)
      !! Fit mixed Gaussian and Student-t MAR components by generalized EM.
      real(dp), intent(in) :: series(:) !! Observed time series.
      type(mixar_model_t), intent(in) :: initial_model !! Initial mixed-distribution MAR parameters.
      logical, intent(in), optional :: estimate_shift !! Estimate component intercepts.
      integer, intent(in), optional :: max_iterations !! Maximum generalized EM iterations.
      real(dp), intent(in), optional :: tolerance !! Relative likelihood tolerance.
      real(dp), intent(in), optional :: minimum_scale !! Innovation scale floor.
      integer, intent(in), optional :: inner_iterations !! Location and scale IRLS iterations per M-step.
      logical, intent(in), optional :: estimate_degrees(:) !! Select Student-t degrees of freedom to estimate.
      type(mixar_fit_t) :: out
      type(mixar_filter_t) :: filtered
      type(mixar_model_t) :: current, updated, candidate
      real(dp), allocatable :: cross(:, :), rhs(:, :), solution(:, :)
      real(dp), allocatable :: residual(:), robust_weight(:)
      logical, allocatable :: fit_degrees(:)
      real(dp) :: previous, candidate_likelihood, tolerance_value, scale_floor
      real(dp) :: effective_count, degrees, standardized, old_scale
      integer :: limit, inner_limit, iteration, inner, component, row, lag
      integer :: time, q, solve_info, parameters
      logical :: fit_shift

      limit = 200
      inner_limit = 5
      tolerance_value = 1.0e-9_dp
      scale_floor = 1.0e-7_dp
      fit_shift = .true.
      if (present(max_iterations)) limit = max_iterations
      if (present(inner_iterations)) inner_limit = inner_iterations
      if (present(tolerance)) tolerance_value = tolerance
      if (present(minimum_scale)) scale_floor = minimum_scale
      if (present(estimate_shift)) fit_shift = estimate_shift
      if (.not. valid_model(initial_model)) then
         out%info = 1
         return
      end if
      allocate(fit_degrees(size(initial_model%probability)), source=.false.)
      if (present(estimate_degrees)) then
         if (size(estimate_degrees) /= size(fit_degrees)) then
            out%info = 1
            return
         end if
         fit_degrees = estimate_degrees
      end if
      if (size(series) <= maxval(initial_model%order) .or. &
         limit < 1 .or. inner_limit < 1 .or. tolerance_value <= 0.0_dp .or. &
         scale_floor <= 0.0_dp .or. any(fit_degrees .and. &
         initial_model%degrees_of_freedom <= 2.0_dp)) then
         out%info = 1
         return
      end if
      current = initial_model
      previous = mixar_log_likelihood(current, series)
      do iteration = 1, limit
         filtered = mixar_filter(current, series)
         if (filtered%info /= 0) then
            out%info = 2
            return
         end if
         updated = current
         do component = 1, size(current%order)
            effective_count = sum(filtered%responsibility(:, component))
            updated%probability(component) = effective_count/ &
               real(size(filtered%responsibility, 1), dp)
            q = current%order(component) + merge(1, 0, fit_shift)
            allocate(residual(size(filtered%responsibility, 1)))
            allocate(robust_weight(size(filtered%responsibility, 1)))
            do inner = 1, inner_limit
               call component_residuals(series, filtered%start, updated, component, residual)
               degrees = current%degrees_of_freedom(component)
               if (degrees > 2.0_dp) then
                  do row = 1, size(residual)
                     standardized = residual(row)/updated%scale(component)
                     robust_weight(row) = (degrees + 1.0_dp)/ &
                        (degrees - 2.0_dp + standardized**2)
                  end do
               else
                  robust_weight = 1.0_dp
               end if
               if (q > 0 .and. effective_count > tiny(1.0_dp)) then
                  allocate(cross(q, q), rhs(q, 1), solution(q, 1))
                  cross = 0.0_dp
                  rhs = 0.0_dp
                  do row = 1, size(residual)
                     time = filtered%start + row - 1
                     call accumulate_normal_equations(series, time, &
                        filtered%responsibility(row, component)*robust_weight(row), &
                        current%order(component), fit_shift, cross, rhs(:, 1))
                  end do
                  call solve_matrix(cross, rhs, solution, solve_info)
                  if (solve_info /= 0) then
                     out%info = 3
                     return
                  end if
                  if (fit_shift) then
                     updated%shift(component) = solution(1, 1)
                     if (current%order(component) > 0) &
                        updated%ar(1:current%order(component), component) = solution(2:q, 1)
                  else if (current%order(component) > 0) then
                     updated%ar(1:current%order(component), component) = solution(:, 1)
                  end if
                  deallocate(cross, rhs, solution)
               end if
               call component_residuals(series, filtered%start, updated, component, residual)
               old_scale = updated%scale(component)
               if (effective_count > tiny(1.0_dp)) updated%scale(component) = &
                  max(scale_floor, sqrt(dot_product(filtered%responsibility(:, component)* &
                  robust_weight, residual**2)/effective_count))
               if (abs(updated%scale(component) - old_scale) <= &
                  1.0e-8_dp*(1.0_dp + old_scale)) exit
            end do
            if (fit_degrees(component)) updated%degrees_of_freedom(component) = &
               estimate_student_degrees(residual/updated%scale(component), &
               filtered%responsibility(:, component))
            deallocate(residual, robust_weight)
         end do
         updated%probability = updated%probability/sum(updated%probability)
         candidate = updated
         candidate_likelihood = mixar_log_likelihood(candidate, series)
         if (.not. ieee_is_finite(candidate_likelihood)) then
            out%info = 4
            return
         end if
         if (candidate_likelihood + 1.0e-10_dp < previous) then
            call damp_general_step(current, updated, series, previous, candidate, &
               candidate_likelihood)
         end if
         current = candidate
         out%iterations = iteration
         if (abs(candidate_likelihood - previous) <= &
            tolerance_value*(1.0_dp + abs(previous))) then
            out%converged = .true.
            exit
         end if
         previous = candidate_likelihood
      end do
      out%model = current
      out%filter = mixar_filter(current, series)
      out%log_likelihood = out%filter%log_likelihood
      parameters = size(current%probability) - 1 + size(current%scale) + &
         sum(current%order) + merge(size(current%shift), 0, fit_shift) + &
         count(fit_degrees)
      out%aic = -2.0_dp*out%log_likelihood + 2.0_dp*real(parameters, dp)
      out%bic = -2.0_dp*out%log_likelihood + &
         log(real(size(out%filter%residual), dp))*real(parameters, dp)
   end function mixar_general_fit

   pure real(dp) function estimate_student_degrees(standardized_residual, &
      weight) result(degrees)
      !! Maximize the weighted standardized Student-t likelihood over its degrees.
      real(dp), intent(in) :: standardized_residual(:) !! Residuals divided by component scale.
      real(dp), intent(in) :: weight(:) !! Posterior component responsibilities.
      real(dp), parameter :: lower = 2.05_dp, upper = 200.0_dp
      real(dp), parameter :: ratio = 0.6180339887498948482_dp
      real(dp) :: left, right, first, second, first_value, second_value
      integer :: iteration

      left = lower
      right = upper
      first = right - ratio*(right - left)
      second = left + ratio*(right - left)
      first_value = weighted_student_log_likelihood(standardized_residual, &
         weight, first)
      second_value = weighted_student_log_likelihood(standardized_residual, &
         weight, second)
      do iteration = 1, 100
         if (first_value < second_value) then
            left = first
            first = second
            first_value = second_value
            second = left + ratio*(right - left)
            second_value = weighted_student_log_likelihood(standardized_residual, &
               weight, second)
         else
            right = second
            second = first
            second_value = first_value
            first = right - ratio*(right - left)
            first_value = weighted_student_log_likelihood(standardized_residual, &
               weight, first)
         end if
         if (right - left <= 1.0e-7_dp*(1.0_dp + left)) exit
      end do
      degrees = 0.5_dp*(left + right)
   end function estimate_student_degrees

   pure real(dp) function weighted_student_log_likelihood(standardized_residual, &
      weight, degrees) result(value)
      !! Evaluate a responsibility-weighted standardized Student-t log likelihood.
      real(dp), intent(in) :: standardized_residual(:) !! Standardized component residuals.
      real(dp), intent(in) :: weight(:) !! Posterior component responsibilities.
      real(dp), intent(in) :: degrees !! Candidate degrees of freedom above two.
      real(dp) :: constant

      constant = log_gamma(0.5_dp*(degrees + 1.0_dp)) - &
         log_gamma(0.5_dp*degrees) - &
         0.5_dp*log(acos(-1.0_dp)*(degrees - 2.0_dp))
      value = sum(weight*(constant - 0.5_dp*(degrees + 1.0_dp)* &
         log(1.0_dp + standardized_residual**2/(degrees - 2.0_dp))))
   end function weighted_student_log_likelihood

   pure subroutine component_residuals(series, start, model, component, residual)
      !! Compute one component's conditional residuals.
      real(dp), intent(in) :: series(:) !! Observed time series.
      integer, intent(in) :: start !! First fitted observation index.
      type(mixar_model_t), intent(in) :: model !! Current mixed-distribution MAR model.
      integer, intent(in) :: component !! Component index.
      real(dp), intent(out) :: residual(:) !! Conditional residuals.
      integer :: row, time, lag

      do row = 1, size(residual)
         time = start + row - 1
         residual(row) = series(time) - model%shift(component)
         do lag = 1, model%order(component)
            residual(row) = residual(row) - model%ar(lag, component)*series(time - lag)
         end do
      end do
   end subroutine component_residuals

   pure subroutine damp_general_step(old_model, proposed_model, series, &
      old_likelihood, accepted_model, accepted_likelihood)
      !! Backtrack a generalized EM step until it does not reduce likelihood.
      type(mixar_model_t), intent(in) :: old_model !! Model before the M-step.
      type(mixar_model_t), intent(in) :: proposed_model !! Full proposed M-step model.
      real(dp), intent(in) :: series(:) !! Observed time series.
      real(dp), intent(in) :: old_likelihood !! Likelihood before the M-step.
      type(mixar_model_t), intent(out) :: accepted_model !! Accepted damped model.
      real(dp), intent(out) :: accepted_likelihood !! Accepted model likelihood.
      real(dp) :: fraction
      integer :: attempt

      fraction = 0.5_dp
      do attempt = 1, 30
         accepted_model = old_model
         accepted_model%probability = old_model%probability + fraction* &
            (proposed_model%probability - old_model%probability)
         accepted_model%shift = old_model%shift + fraction* &
            (proposed_model%shift - old_model%shift)
         accepted_model%scale = old_model%scale + fraction* &
            (proposed_model%scale - old_model%scale)
         accepted_model%ar = old_model%ar + fraction*(proposed_model%ar - old_model%ar)
         accepted_likelihood = mixar_log_likelihood(accepted_model, series)
         if (accepted_likelihood >= old_likelihood - 1.0e-10_dp) return
         fraction = 0.5_dp*fraction
      end do
      accepted_model = old_model
      accepted_likelihood = old_likelihood
   end subroutine damp_general_step

   pure function mixar_regression_model(mar, coefficient, &
      component_specific) result(model)
      !! Construct a MAR regression model from an innovation model and coefficients.
      type(mixar_model_t), intent(in) :: mar !! MAR model for regression errors.
      real(dp), intent(in) :: coefficient(:, :) !! Regressor-by-one or regressor-by-component coefficients.
      logical, intent(in), optional :: component_specific !! Require one coefficient column per component.
      type(mixar_regression_model_t) :: model

      model%mar = mar
      model%coefficient = coefficient
      model%component_specific = size(coefficient, 2) > 1
      if (present(component_specific)) model%component_specific = component_specific
      if (.not. valid_regression_model(model)) model%info = 1
   end function mixar_regression_model

   pure function mixar_regression_filter(model, response, regressors) result(out)
      !! Filter a MAR regression and compute component responsibilities.
      type(mixar_regression_model_t), intent(in) :: model !! MAR regression model.
      real(dp), intent(in) :: response(:) !! Response time series.
      real(dp), intent(in) :: regressors(:, :) !! Observation-by-regressor design matrix.
      type(mixar_filter_t) :: out
      real(dp), allocatable :: log_kernel(:)
      real(dp) :: maximum, denominator
      integer :: observations, components, row, component, time

      if (.not. valid_regression_model(model) .or. &
         size(regressors, 1) /= size(response) .or. &
         size(regressors, 2) /= size(model%coefficient, 1) .or. &
         size(response) <= maxval(model%mar%order) .or. &
         .not. all(ieee_is_finite(response)) .or. &
         .not. all(ieee_is_finite(regressors))) then
         out%info = 1
         return
      end if
      out%start = maxval(model%mar%order) + 1
      observations = size(response) - out%start + 1
      components = size(model%mar%probability)
      allocate(out%location(observations, components))
      allocate(out%responsibility(observations, components), out%fitted(observations))
      allocate(out%variance(observations), out%residual(observations), log_kernel(components))
      out%log_likelihood = 0.0_dp
      do row = 1, observations
         time = out%start + row - 1
         do component = 1, components
            out%location(row, component) = regression_component_location(model, &
               response, regressors, time, component)
            log_kernel(component) = log(max(model%mar%probability(component), &
               tiny(1.0_dp))) - log(model%mar%scale(component)) + &
               log(max(tiny(1.0_dp), mixar_standard_density((response(time) - &
               out%location(row, component))/model%mar%scale(component), &
               model%mar%degrees_of_freedom(component))))
         end do
         maximum = maxval(log_kernel)
         denominator = sum(exp(log_kernel - maximum))
         out%responsibility(row, :) = exp(log_kernel - maximum)/denominator
         out%log_likelihood = out%log_likelihood + maximum + log(denominator)
         out%fitted(row) = dot_product(model%mar%probability, out%location(row, :))
         out%variance(row) = dot_product(model%mar%probability, &
            model%mar%scale**2 + out%location(row, :)**2) - out%fitted(row)**2
         out%residual(row) = response(time) - out%fitted(row)
      end do
   end function mixar_regression_filter

   pure real(dp) function mixar_regression_log_likelihood(model, response, &
      regressors) result(log_likelihood)
      !! Evaluate a MAR regression conditional log likelihood.
      type(mixar_regression_model_t), intent(in) :: model !! MAR regression model.
      real(dp), intent(in) :: response(:) !! Response time series.
      real(dp), intent(in) :: regressors(:, :) !! Observation-by-regressor design matrix.
      type(mixar_filter_t) :: filtered

      filtered = mixar_regression_filter(model, response, regressors)
      log_likelihood = filtered%log_likelihood
   end function mixar_regression_log_likelihood

   pure real(dp) function regression_component_location(model, response, &
      regressors, time, component) result(location)
      !! Evaluate one regression-MAR component's conditional location.
      type(mixar_regression_model_t), intent(in) :: model !! MAR regression model.
      real(dp), intent(in) :: response(:) !! Response time series.
      real(dp), intent(in) :: regressors(:, :) !! Observation-by-regressor design matrix.
      integer, intent(in) :: time !! Observation index.
      integer, intent(in) :: component !! Mixture component index.
      real(dp) :: lag_residual
      integer :: column, lag

      column = regression_column(model, component)
      location = model%mar%shift(component) + &
         dot_product(regressors(time, :), model%coefficient(:, column))
      do lag = 1, model%mar%order(component)
         lag_residual = response(time - lag) - &
            dot_product(regressors(time - lag, :), model%coefficient(:, column))
         location = location + model%mar%ar(lag, component)*lag_residual
      end do
   end function regression_component_location

   pure integer function regression_column(model, component) result(column)
      !! Select the shared or component-specific regression coefficient column.
      type(mixar_regression_model_t), intent(in) :: model !! MAR regression model.
      integer, intent(in) :: component !! Mixture component index.

      if (model%component_specific) then
         column = component
      else
         column = 1
      end if
   end function regression_column

   pure function mixar_regression_fit(response, regressors, initial_model, &
      estimate_shift, max_iterations, tolerance, minimum_scale, &
      estimate_degrees) result(out)
      !! Fit shared or component-specific MAR regression effects by generalized EM.
      real(dp), intent(in) :: response(:) !! Response time series.
      real(dp), intent(in) :: regressors(:, :) !! Observation-by-regressor design matrix.
      type(mixar_regression_model_t), intent(in) :: initial_model !! Initial MAR regression model.
      logical, intent(in), optional :: estimate_shift !! Estimate MAR component shifts.
      integer, intent(in), optional :: max_iterations !! Maximum generalized EM iterations.
      real(dp), intent(in), optional :: tolerance !! Relative likelihood tolerance.
      real(dp), intent(in), optional :: minimum_scale !! Innovation scale floor.
      logical, intent(in), optional :: estimate_degrees(:) !! Select Student-t degrees to estimate.
      type(mixar_regression_fit_t) :: out
      type(mixar_regression_model_t) :: current, updated, candidate
      type(mixar_filter_t) :: filtered
      real(dp), allocatable :: cross(:, :), rhs(:, :), solution(:, :)
      real(dp), allocatable :: error_series(:), robust_weight(:, :), residual(:)
      real(dp), allocatable :: transformed_regressor(:)
      logical, allocatable :: fit_degrees(:)
      real(dp) :: previous, candidate_likelihood, tolerance_value, scale_floor
      real(dp) :: effective_count, degrees, standardized, weight
      integer :: limit, iteration, component, row, time, lag, q, column
      integer :: solve_info, components, observations, parameters
      logical :: fit_shift

      limit = 200
      tolerance_value = 1.0e-7_dp
      scale_floor = 1.0e-7_dp
      fit_shift = .false.
      if (present(max_iterations)) limit = max_iterations
      if (present(tolerance)) tolerance_value = tolerance
      if (present(minimum_scale)) scale_floor = minimum_scale
      if (present(estimate_shift)) fit_shift = estimate_shift
      if (.not. valid_regression_model(initial_model) .or. &
         size(regressors, 1) /= size(response) .or. &
         size(regressors, 2) /= size(initial_model%coefficient, 1) .or. &
         size(response) <= maxval(initial_model%mar%order) .or. limit < 1 .or. &
         tolerance_value <= 0.0_dp .or. scale_floor <= 0.0_dp) then
         out%info = 1
         return
      end if
      components = size(initial_model%mar%probability)
      allocate(fit_degrees(components), source=.false.)
      if (present(estimate_degrees)) then
         if (size(estimate_degrees) /= components) then
            out%info = 1
            return
         end if
         fit_degrees = estimate_degrees
      end if
      if (any(fit_degrees .and. initial_model%mar%degrees_of_freedom <= 2.0_dp)) then
         out%info = 1
         return
      end if
      current = initial_model
      previous = mixar_regression_log_likelihood(current, response, regressors)
      do iteration = 1, limit
         filtered = mixar_regression_filter(current, response, regressors)
         if (filtered%info /= 0) then
            out%info = 2
            return
         end if
         observations = size(filtered%responsibility, 1)
         allocate(robust_weight(observations, components))
         do component = 1, components
            degrees = current%mar%degrees_of_freedom(component)
            do row = 1, observations
               standardized = (response(filtered%start + row - 1) - &
                  filtered%location(row, component))/current%mar%scale(component)
               if (degrees > 2.0_dp) then
                  robust_weight(row, component) = (degrees + 1.0_dp)/ &
                     (degrees - 2.0_dp + standardized**2)
               else
                  robust_weight(row, component) = 1.0_dp
               end if
            end do
         end do
         updated = current
         do component = 1, components
            effective_count = sum(filtered%responsibility(:, component))
            updated%mar%probability(component) = effective_count/real(observations, dp)
            q = updated%mar%order(component) + merge(1, 0, fit_shift)
            if (q > 0 .and. effective_count > tiny(1.0_dp)) then
               column = regression_column(updated, component)
               error_series = response - matmul(regressors, updated%coefficient(:, column))
               allocate(cross(q, q), rhs(q, 1), solution(q, 1))
               cross = 0.0_dp
               rhs = 0.0_dp
               do row = 1, observations
                  time = filtered%start + row - 1
                  weight = filtered%responsibility(row, component)* &
                     robust_weight(row, component)
                  call accumulate_normal_equations(error_series, time, weight, &
                     updated%mar%order(component), fit_shift, cross, rhs(:, 1))
               end do
               call solve_matrix(cross, rhs, solution, solve_info)
               if (solve_info /= 0) then
                  out%info = 3
                  return
               end if
               if (fit_shift) then
                  updated%mar%shift(component) = solution(1, 1)
                  if (updated%mar%order(component) > 0) &
                     updated%mar%ar(1:updated%mar%order(component), component) = &
                     solution(2:q, 1)
               else if (updated%mar%order(component) > 0) then
                  updated%mar%ar(1:updated%mar%order(component), component) = &
                     solution(:, 1)
               end if
               deallocate(cross, rhs, solution)
            end if
         end do
         updated%mar%probability = updated%mar%probability/ &
            sum(updated%mar%probability)
         do column = 1, size(updated%coefficient, 2)
            q = size(updated%coefficient, 1)
            allocate(cross(q, q), rhs(q, 1), solution(q, 1))
            allocate(transformed_regressor(q))
            cross = 0.0_dp
            rhs = 0.0_dp
            do component = 1, components
               if (updated%component_specific .and. component /= column) cycle
               do row = 1, observations
                  time = filtered%start + row - 1
                  transformed_regressor = regressors(time, :)
                  standardized = response(time) - updated%mar%shift(component)
                  do lag = 1, updated%mar%order(component)
                     transformed_regressor = transformed_regressor - &
                        updated%mar%ar(lag, component)*regressors(time - lag, :)
                     standardized = standardized - &
                        updated%mar%ar(lag, component)*response(time - lag)
                  end do
                  weight = filtered%responsibility(row, component)* &
                     robust_weight(row, component)
                  cross = cross + weight*spread(transformed_regressor, 2, q)* &
                     spread(transformed_regressor, 1, q)
                  rhs(:, 1) = rhs(:, 1) + weight*transformed_regressor*standardized
               end do
            end do
            call solve_matrix(cross, rhs, solution, solve_info)
            if (solve_info /= 0) then
               out%info = 4
               return
            end if
            updated%coefficient(:, column) = solution(:, 1)
            deallocate(cross, rhs, solution, transformed_regressor)
         end do
         allocate(residual(observations))
         do component = 1, components
            do row = 1, observations
               time = filtered%start + row - 1
               residual(row) = response(time) - regression_component_location( &
                  updated, response, regressors, time, component)
            end do
            effective_count = sum(filtered%responsibility(:, component))
            if (effective_count > tiny(1.0_dp)) updated%mar%scale(component) = &
               max(scale_floor, sqrt(dot_product(filtered%responsibility(:, component)* &
               robust_weight(:, component), residual**2)/effective_count))
            if (fit_degrees(component)) updated%mar%degrees_of_freedom(component) = &
               estimate_student_degrees(residual/updated%mar%scale(component), &
               filtered%responsibility(:, component))
         end do
         deallocate(residual, robust_weight)
         candidate = updated
         candidate_likelihood = mixar_regression_log_likelihood(candidate, &
            response, regressors)
         if (.not. ieee_is_finite(candidate_likelihood)) then
            out%info = 5
            return
         end if
         if (candidate_likelihood + 1.0e-10_dp < previous) &
            call damp_regression_step(current, updated, response, regressors, &
            previous, candidate, candidate_likelihood)
         current = candidate
         out%iterations = iteration
         if (abs(candidate_likelihood - previous) <= &
            tolerance_value*(1.0_dp + abs(previous))) then
            out%converged = .true.
            exit
         end if
         previous = candidate_likelihood
      end do
      out%model = current
      out%filter = mixar_regression_filter(current, response, regressors)
      out%log_likelihood = out%filter%log_likelihood
      parameters = components - 1 + components + sum(current%mar%order) + &
         merge(components, 0, fit_shift) + size(current%coefficient) + &
         count(fit_degrees)
      out%aic = -2.0_dp*out%log_likelihood + 2.0_dp*real(parameters, dp)
      out%bic = -2.0_dp*out%log_likelihood + &
         log(real(size(out%filter%residual), dp))*real(parameters, dp)
   end function mixar_regression_fit

   pure subroutine damp_regression_step(old_model, proposed_model, response, &
      regressors, old_likelihood, accepted_model, accepted_likelihood)
      !! Backtrack a MAR regression update until likelihood does not decrease.
      type(mixar_regression_model_t), intent(in) :: old_model !! Model before the update.
      type(mixar_regression_model_t), intent(in) :: proposed_model !! Full proposed update.
      real(dp), intent(in) :: response(:) !! Response time series.
      real(dp), intent(in) :: regressors(:, :) !! Observation-by-regressor design matrix.
      real(dp), intent(in) :: old_likelihood !! Likelihood before the update.
      type(mixar_regression_model_t), intent(out) :: accepted_model !! Accepted model.
      real(dp), intent(out) :: accepted_likelihood !! Accepted model likelihood.
      real(dp) :: fraction
      integer :: attempt

      fraction = 0.5_dp
      do attempt = 1, 30
         accepted_model = old_model
         accepted_model%coefficient = old_model%coefficient + fraction* &
            (proposed_model%coefficient - old_model%coefficient)
         accepted_model%mar%probability = old_model%mar%probability + fraction* &
            (proposed_model%mar%probability - old_model%mar%probability)
         accepted_model%mar%shift = old_model%mar%shift + fraction* &
            (proposed_model%mar%shift - old_model%mar%shift)
         accepted_model%mar%scale = old_model%mar%scale + fraction* &
            (proposed_model%mar%scale - old_model%mar%scale)
         accepted_model%mar%ar = old_model%mar%ar + fraction* &
            (proposed_model%mar%ar - old_model%mar%ar)
         accepted_likelihood = mixar_regression_log_likelihood(accepted_model, &
            response, regressors)
         if (accepted_likelihood >= old_likelihood - 1.0e-10_dp) return
         fraction = 0.5_dp*fraction
      end do
      accepted_model = old_model
      accepted_likelihood = old_likelihood
   end subroutine damp_regression_step

   pure function mixar_regression_simulate_from_draws(model, initial_response, &
      initial_regressors, future_regressors, uniforms, standard_innovations) &
      result(out)
      !! Simulate a MAR regression from supplied regimes and standardized innovations.
      type(mixar_regression_model_t), intent(in) :: model !! MAR regression model.
      real(dp), intent(in) :: initial_response(:) !! Response history ending at the origin.
      real(dp), intent(in) :: initial_regressors(:, :) !! Regressor history aligned with responses.
      real(dp), intent(in) :: future_regressors(:, :) !! Future observation-by-regressor matrix.
      real(dp), intent(in) :: uniforms(:) !! Future regime-selection uniform draws.
      real(dp), intent(in) :: standard_innovations(:) !! Future standardized innovations.
      type(mixar_regression_simulation_t) :: out
      real(dp), allocatable :: response_work(:), regressor_work(:, :)
      real(dp) :: location, lag_residual
      integer :: p, observations, time, regime, column, lag

      p = maxval(model%mar%order)
      observations = size(future_regressors, 1)
      if (.not. valid_regression_model(model) .or. size(initial_response) < p .or. &
         size(initial_regressors, 1) < p .or. &
         size(initial_regressors, 2) /= size(model%coefficient, 1) .or. &
         size(future_regressors, 2) /= size(model%coefficient, 1) .or. &
         observations < 1 .or. size(uniforms) < observations .or. &
         size(standard_innovations) < observations .or. &
         any(uniforms(:observations) < 0.0_dp) .or. &
         any(uniforms(:observations) > 1.0_dp)) then
         out%info = 1
         return
      end if
      allocate(response_work(p + observations))
      allocate(regressor_work(p + observations, size(model%coefficient, 1)))
      allocate(out%series(observations), out%regime(observations))
      allocate(out%innovation(observations))
      if (p > 0) then
         response_work(:p) = initial_response(size(initial_response) - p + 1:)
         regressor_work(:p, :) = initial_regressors( &
            size(initial_regressors, 1) - p + 1:, :)
      end if
      regressor_work(p + 1:, :) = future_regressors
      do time = 1, observations
         regime = selected_component(model%mar%probability, uniforms(time))
         column = regression_column(model, regime)
         location = model%mar%shift(regime) + &
            dot_product(regressor_work(p + time, :), model%coefficient(:, column))
         do lag = 1, model%mar%order(regime)
            lag_residual = response_work(p + time - lag) - &
               dot_product(regressor_work(p + time - lag, :), &
               model%coefficient(:, column))
            location = location + model%mar%ar(lag, regime)*lag_residual
         end do
         out%innovation(time) = model%mar%scale(regime)*standard_innovations(time)
         response_work(p + time) = location + out%innovation(time)
         out%series(time) = response_work(p + time)
         out%regime(time) = regime
      end do
      out%model = model
   end function mixar_regression_simulate_from_draws

   function mixar_regression_simulate(model, initial_response, &
      initial_regressors, future_regressors) result(out)
      !! Simulate a MAR regression using the shared random generator.
      type(mixar_regression_model_t), intent(in) :: model !! MAR regression model.
      real(dp), intent(in) :: initial_response(:) !! Response history ending at the origin.
      real(dp), intent(in) :: initial_regressors(:, :) !! Regressor history aligned with responses.
      real(dp), intent(in) :: future_regressors(:, :) !! Future observation-by-regressor matrix.
      type(mixar_regression_simulation_t) :: out
      real(dp), allocatable :: uniforms(:), innovations(:)
      integer :: observations, time, regime

      observations = size(future_regressors, 1)
      allocate(uniforms(observations), innovations(observations))
      do time = 1, observations
         uniforms(time) = random_uniform()
         regime = selected_component(model%mar%probability, uniforms(time))
         if (model%mar%degrees_of_freedom(regime) > 2.0_dp) then
            innovations(time) = random_standard_student( &
               model%mar%degrees_of_freedom(regime))
         else
            innovations(time) = random_standard_normal()
         end if
      end do
      out = mixar_regression_simulate_from_draws(model, initial_response, &
         initial_regressors, future_regressors, uniforms, innovations)
   end function mixar_regression_simulate

   pure integer function selected_component(probability, uniform) result(component)
      !! Select a mixture component from one unit-interval draw.
      real(dp), intent(in) :: probability(:) !! Component mixing probabilities.
      real(dp), intent(in) :: uniform !! Unit-interval draw.
      real(dp) :: cumulative

      cumulative = 0.0_dp
      component = size(probability)
      do component = 1, size(probability)
         cumulative = cumulative + probability(component)
         if (uniform <= cumulative) return
      end do
      component = size(probability)
   end function selected_component

   pure function mixar_regression_exact_forecast(model, initial_response, &
      initial_regressors, future_regressors, level, maximum_paths) result(out)
      !! Construct exact Gaussian MAR regression forecasts for known future regressors.
      type(mixar_regression_model_t), intent(in) :: model !! Gaussian MAR regression model.
      real(dp), intent(in) :: initial_response(:) !! Response history ending at the origin.
      real(dp), intent(in) :: initial_regressors(:, :) !! Regressor history aligned with responses.
      real(dp), intent(in) :: future_regressors(:, :) !! Future observation-by-regressor matrix.
      real(dp), intent(in), optional :: level !! Central prediction interval probability.
      integer, intent(in), optional :: maximum_paths !! Maximum regime paths at any horizon.
      type(mixar_exact_forecast_t) :: out
      real(dp), allocatable :: response_history(:), regressor_history(:, :)
      real(dp) :: interval_level, tail_probability
      integer :: horizons, horizon, paths, path_limit, p

      interval_level = 0.95_dp
      path_limit = 1000000
      if (present(level)) interval_level = level
      if (present(maximum_paths)) path_limit = maximum_paths
      horizons = size(future_regressors, 1)
      p = maxval(model%mar%order)
      if (.not. valid_regression_model(model) .or. &
         any(model%mar%degrees_of_freedom /= 0.0_dp) .or. &
         size(initial_response) < p .or. size(initial_regressors, 1) < p .or. &
         size(initial_regressors, 2) /= size(model%coefficient, 1) .or. &
         size(future_regressors, 2) /= size(model%coefficient, 1) .or. &
         horizons < 1 .or. interval_level <= 0.0_dp .or. interval_level >= 1.0_dp) then
         out%info = 1
         return
      end if
      paths = 1
      do horizon = 1, horizons
         if (paths > path_limit/size(model%mar%probability)) then
            out%info = 2
            return
         end if
         paths = paths*size(model%mar%probability)
      end do
      allocate(response_history(p), regressor_history(p + horizons, &
         size(model%coefficient, 1)))
      if (p > 0) then
         response_history = initial_response(size(initial_response) - p + 1:)
         regressor_history(:p, :) = initial_regressors( &
            size(initial_regressors, 1) - p + 1:, :)
      end if
      regressor_history(p + 1:, :) = future_regressors
      allocate(out%distribution(horizons), out%mean(horizons))
      allocate(out%variance(horizons), out%lower(horizons), out%upper(horizons))
      out%level = interval_level
      tail_probability = 0.5_dp*(1.0_dp - interval_level)
      do horizon = 1, horizons
         out%distribution(horizon) = exact_regression_horizon_distribution( &
            model, response_history, regressor_history, p, horizon)
         out%mean(horizon) = dot_product(out%distribution(horizon)%probability, &
            out%distribution(horizon)%location)
         out%variance(horizon) = dot_product(out%distribution(horizon)%probability, &
            out%distribution(horizon)%scale**2 + &
            out%distribution(horizon)%location**2) - out%mean(horizon)**2
         out%lower(horizon) = mixar_predictive_quantile(out%distribution(horizon), &
            tail_probability)
         out%upper(horizon) = mixar_predictive_quantile(out%distribution(horizon), &
            1.0_dp - tail_probability)
      end do
   end function mixar_regression_exact_forecast

   pure function exact_regression_horizon_distribution(model, initial_response, &
      regressors, history_rows, horizon) result(out)
      !! Propagate each regression-MAR regime path into a Gaussian component.
      type(mixar_regression_model_t), intent(in) :: model !! Gaussian MAR regression model.
      real(dp), intent(in) :: initial_response(:) !! Response history used as state.
      real(dp), intent(in) :: regressors(:, :) !! Historical and future regressors.
      integer, intent(in) :: history_rows !! Number of historical regressor rows.
      integer, intent(in) :: horizon !! Forecast horizon.
      type(mixar_predictive_distribution_t) :: out
      real(dp), allocatable :: state(:), next_state(:), noise(:, :), next_noise(:, :)
      real(dp), allocatable :: coefficient(:)
      integer, allocatable :: regimes(:)
      real(dp) :: intercept
      integer :: p, paths, path, code, step, regime, row, lag, column

      p = max(1, maxval(model%mar%order))
      paths = size(model%mar%probability)**horizon
      allocate(out%probability(paths), out%location(paths), out%scale(paths))
      allocate(out%regime_path(horizon, paths), state(p), next_state(p))
      allocate(noise(p, horizon), next_noise(p, horizon), coefficient(p), regimes(horizon))
      out%horizon = horizon
      do path = 1, paths
         code = path - 1
         do step = horizon, 1, -1
            regimes(step) = modulo(code, size(model%mar%probability)) + 1
            code = code/size(model%mar%probability)
         end do
         out%regime_path(:, path) = regimes
         state = 0.0_dp
         do row = 1, size(initial_response)
            state(row) = initial_response(size(initial_response) - row + 1)
         end do
         noise = 0.0_dp
         out%probability(path) = 1.0_dp
         do step = 1, horizon
            regime = regimes(step)
            column = regression_column(model, regime)
            coefficient = 0.0_dp
            if (model%mar%order(regime) > 0) coefficient(1:model%mar%order(regime)) = &
               model%mar%ar(1:model%mar%order(regime), regime)
            intercept = model%mar%shift(regime) + &
               dot_product(regressors(history_rows + step, :), &
               model%coefficient(:, column))
            do lag = 1, model%mar%order(regime)
               intercept = intercept - model%mar%ar(lag, regime)* &
                  dot_product(regressors(history_rows + step - lag, :), &
                  model%coefficient(:, column))
            end do
            next_state = 0.0_dp
            next_state(1) = intercept + dot_product(coefficient, state)
            if (p > 1) next_state(2:p) = state(1:p - 1)
            next_noise = 0.0_dp
            next_noise(1, :) = matmul(coefficient, noise)
            next_noise(1, step) = next_noise(1, step) + model%mar%scale(regime)
            if (p > 1) next_noise(2:p, :) = noise(1:p - 1, :)
            state = next_state
            noise = next_noise
            out%probability(path) = out%probability(path)*model%mar%probability(regime)
         end do
         out%location(path) = state(1)
         out%scale(path) = sqrt(sum(noise(1, :)**2))
      end do
      out%probability = out%probability/sum(out%probability)
   end function exact_regression_horizon_distribution

   pure function mixar_regression_forecast_from_draws(model, initial_response, &
      initial_regressors, future_regressors, uniforms, standard_innovations, &
      level) result(out)
      !! Forecast a MAR regression from supplied path draws and future regressors.
      type(mixar_regression_model_t), intent(in) :: model !! MAR regression model.
      real(dp), intent(in) :: initial_response(:) !! Response history ending at the origin.
      real(dp), intent(in) :: initial_regressors(:, :) !! Regressor history aligned with responses.
      real(dp), intent(in) :: future_regressors(:, :) !! Future observation-by-regressor matrix.
      real(dp), intent(in) :: uniforms(:, :) !! Horizon-by-simulation regime uniform draws.
      real(dp), intent(in) :: standard_innovations(:, :) !! Horizon-by-simulation standardized innovations.
      real(dp), intent(in), optional :: level !! Central empirical interval probability.
      type(mixar_forecast_sample_t) :: out
      type(mixar_regression_simulation_t) :: simulation
      real(dp), allocatable :: ordered(:)
      real(dp) :: interval_level, tail_probability
      integer :: horizons, simulations, sample, horizon

      interval_level = 0.95_dp
      if (present(level)) interval_level = level
      horizons = size(future_regressors, 1)
      if (.not. valid_regression_model(model) .or. &
         any(shape(uniforms) /= shape(standard_innovations)) .or. &
         size(uniforms, 1) /= horizons .or. size(uniforms, 2) < 1 .or. &
         interval_level <= 0.0_dp .or. interval_level >= 1.0_dp) then
         out%info = 1
         return
      end if
      simulations = size(uniforms, 2)
      allocate(out%paths(horizons, simulations), out%mean(horizons))
      allocate(out%variance(horizons), out%lower(horizons), out%upper(horizons))
      do sample = 1, simulations
         simulation = mixar_regression_simulate_from_draws(model, initial_response, &
            initial_regressors, future_regressors, uniforms(:, sample), &
            standard_innovations(:, sample))
         if (simulation%info /= 0) then
            out%info = 2
            return
         end if
         out%paths(:, sample) = simulation%series
      end do
      out%level = interval_level
      tail_probability = 0.5_dp*(1.0_dp - interval_level)
      do horizon = 1, horizons
         ordered = sorted(out%paths(horizon, :))
         out%mean(horizon) = sum(ordered)/real(simulations, dp)
         out%variance(horizon) = variance(ordered)
         out%lower(horizon) = quantile(ordered, tail_probability)
         out%upper(horizon) = quantile(ordered, 1.0_dp - tail_probability)
      end do
   end function mixar_regression_forecast_from_draws

   function mixar_regression_forecast(model, initial_response, &
      initial_regressors, future_regressors, simulations, level) result(out)
      !! Forecast a MAR regression with future regressors and the shared RNG.
      type(mixar_regression_model_t), intent(in) :: model !! MAR regression model.
      real(dp), intent(in) :: initial_response(:) !! Response history ending at the origin.
      real(dp), intent(in) :: initial_regressors(:, :) !! Regressor history aligned with responses.
      real(dp), intent(in) :: future_regressors(:, :) !! Future observation-by-regressor matrix.
      integer, intent(in) :: simulations !! Positive number of simulated forecast paths.
      real(dp), intent(in), optional :: level !! Central empirical interval probability.
      type(mixar_forecast_sample_t) :: out
      type(mixar_regression_simulation_t) :: simulation
      real(dp), allocatable :: ordered(:)
      real(dp) :: interval_level, tail_probability
      integer :: horizons, sample, horizon

      interval_level = 0.95_dp
      if (present(level)) interval_level = level
      horizons = size(future_regressors, 1)
      if (.not. valid_regression_model(model) .or. horizons < 1 .or. &
         simulations < 1 .or. interval_level <= 0.0_dp .or. &
         interval_level >= 1.0_dp) then
         out%info = 1
         return
      end if
      allocate(out%paths(horizons, simulations), out%mean(horizons))
      allocate(out%variance(horizons), out%lower(horizons), out%upper(horizons))
      do sample = 1, simulations
         simulation = mixar_regression_simulate(model, initial_response, &
            initial_regressors, future_regressors)
         if (simulation%info /= 0) then
            out%info = 2
            return
         end if
         out%paths(:, sample) = simulation%series
      end do
      out%level = interval_level
      tail_probability = 0.5_dp*(1.0_dp - interval_level)
      do horizon = 1, horizons
         ordered = sorted(out%paths(horizon, :))
         out%mean(horizon) = sum(ordered)/real(simulations, dp)
         out%variance(horizon) = variance(ordered)
         out%lower(horizon) = quantile(ordered, tail_probability)
         out%upper(horizon) = quantile(ordered, 1.0_dp - tail_probability)
      end do
   end function mixar_regression_forecast

   pure function mixar_seasonal_model(mar, seasonal_period, seasonal_order, &
      seasonal_ar) result(model)
      !! Construct an additive seasonal mixture autoregression.
      type(mixar_model_t), intent(in) :: mar !! Ordinary-lag MAR parameters.
      integer, intent(in) :: seasonal_period !! Positive seasonal lag spacing.
      integer, intent(in) :: seasonal_order(:) !! Seasonal order of each component.
      real(dp), intent(in) :: seasonal_ar(:, :) !! Padded seasonal coefficients by order and component.
      type(mixar_seasonal_model_t) :: model

      model%mar = mar
      model%seasonal_period = seasonal_period
      model%seasonal_order = seasonal_order
      model%seasonal_ar = seasonal_ar
      if (.not. valid_seasonal_model(model)) model%info = 1
   end function mixar_seasonal_model

   pure function mixar_seasonal_expanded_model(model) result(expanded)
      !! Expand ordinary and seasonal coefficients onto their actual lag positions.
      type(mixar_seasonal_model_t), intent(in) :: model !! Seasonal MAR model.
      type(mixar_model_t) :: expanded
      real(dp), allocatable :: ar(:, :)
      integer, allocatable :: order(:)
      integer :: maximum_lag, component, lag, actual_lag

      if (.not. valid_seasonal_model(model)) then
         expanded%info = 1
         return
      end if
      maximum_lag = max(maxval(model%mar%order), &
         model%seasonal_period*maxval(model%seasonal_order))
      allocate(ar(maximum_lag, size(model%mar%probability)), source=0.0_dp)
      allocate(order(size(model%mar%order)))
      do component = 1, size(order)
         order(component) = max(model%mar%order(component), &
            model%seasonal_period*model%seasonal_order(component))
         do lag = 1, model%mar%order(component)
            ar(lag, component) = ar(lag, component) + model%mar%ar(lag, component)
         end do
         do lag = 1, model%seasonal_order(component)
            actual_lag = model%seasonal_period*lag
            ar(actual_lag, component) = ar(actual_lag, component) + &
               model%seasonal_ar(lag, component)
         end do
      end do
      expanded = mixar_model(model%mar%probability, order, model%mar%shift, &
         model%mar%scale, ar, model%mar%degrees_of_freedom)
   end function mixar_seasonal_expanded_model

   pure function mixar_seasonal_filter(model, series) result(out)
      !! Filter a seasonal MAR through its sparse expanded lag representation.
      type(mixar_seasonal_model_t), intent(in) :: model !! Seasonal MAR model.
      real(dp), intent(in) :: series(:) !! Observed time series.
      type(mixar_filter_t) :: out

      out = mixar_filter(mixar_seasonal_expanded_model(model), series)
   end function mixar_seasonal_filter

   pure real(dp) function mixar_seasonal_log_likelihood(model, series) &
      result(log_likelihood)
      !! Evaluate the seasonal MAR conditional log likelihood.
      type(mixar_seasonal_model_t), intent(in) :: model !! Seasonal MAR model.
      real(dp), intent(in) :: series(:) !! Observed time series.
      type(mixar_filter_t) :: filtered

      filtered = mixar_seasonal_filter(model, series)
      log_likelihood = filtered%log_likelihood
   end function mixar_seasonal_log_likelihood

   pure function mixar_seasonal_fit(series, initial_model, estimate_shift, &
      max_iterations, tolerance, minimum_scale, estimate_degrees) result(out)
      !! Fit additive ordinary and seasonal MAR coefficients by generalized EM.
      real(dp), intent(in) :: series(:) !! Observed time series.
      type(mixar_seasonal_model_t), intent(in) :: initial_model !! Initial seasonal MAR model.
      logical, intent(in), optional :: estimate_shift !! Estimate component shifts.
      integer, intent(in), optional :: max_iterations !! Maximum generalized EM iterations.
      real(dp), intent(in), optional :: tolerance !! Relative likelihood tolerance.
      real(dp), intent(in), optional :: minimum_scale !! Innovation scale floor.
      logical, intent(in), optional :: estimate_degrees(:) !! Select Student-t degrees to estimate.
      type(mixar_seasonal_fit_t) :: out
      type(mixar_seasonal_model_t) :: current, updated, candidate
      type(mixar_filter_t) :: filtered
      real(dp), allocatable :: cross(:, :), inverse(:, :), rhs(:), solution(:)
      real(dp), allocatable :: regressor(:), residual(:), robust_weight(:)
      logical, allocatable :: fit_degrees(:)
      real(dp) :: previous, candidate_likelihood, tolerance_value, scale_floor
      real(dp) :: effective_count, degrees, standardized, weight
      integer :: limit, iteration, component, row, time, lag, q, offset
      integer :: info, components, observations, parameters
      logical :: fit_shift

      limit = 200
      tolerance_value = 1.0e-9_dp
      scale_floor = 1.0e-7_dp
      fit_shift = .false.
      if (present(max_iterations)) limit = max_iterations
      if (present(tolerance)) tolerance_value = tolerance
      if (present(minimum_scale)) scale_floor = minimum_scale
      if (present(estimate_shift)) fit_shift = estimate_shift
      if (.not. valid_seasonal_model(initial_model) .or. &
         size(series) <= max(maxval(initial_model%mar%order), &
         initial_model%seasonal_period*maxval(initial_model%seasonal_order)) .or. &
         limit < 1 .or. tolerance_value <= 0.0_dp .or. scale_floor <= 0.0_dp) then
         out%info = 1
         return
      end if
      components = size(initial_model%mar%probability)
      allocate(fit_degrees(components), source=.false.)
      if (present(estimate_degrees)) then
         if (size(estimate_degrees) /= components) then
            out%info = 1
            return
         end if
         fit_degrees = estimate_degrees
      end if
      if (any(fit_degrees .and. initial_model%mar%degrees_of_freedom <= 2.0_dp)) then
         out%info = 1
         return
      end if
      current = initial_model
      previous = mixar_seasonal_log_likelihood(current, series)
      do iteration = 1, limit
         filtered = mixar_seasonal_filter(current, series)
         if (filtered%info /= 0) then
            out%info = 2
            return
         end if
         observations = size(filtered%responsibility, 1)
         updated = current
         do component = 1, components
            effective_count = sum(filtered%responsibility(:, component))
            updated%mar%probability(component) = effective_count/real(observations, dp)
            q = updated%mar%order(component) + updated%seasonal_order(component) + &
               merge(1, 0, fit_shift)
            allocate(cross(q, q), inverse(q, q), rhs(q), solution(q), regressor(q))
            allocate(robust_weight(observations), residual(observations))
            cross = 0.0_dp
            rhs = 0.0_dp
            degrees = current%mar%degrees_of_freedom(component)
            do row = 1, observations
               time = filtered%start + row - 1
               standardized = (series(time) - filtered%location(row, component))/ &
                  current%mar%scale(component)
               if (degrees > 2.0_dp) then
                  robust_weight(row) = (degrees + 1.0_dp)/ &
                     (degrees - 2.0_dp + standardized**2)
               else
                  robust_weight(row) = 1.0_dp
               end if
               call seasonal_regressor(series, time, updated, component, &
                  fit_shift, regressor)
               weight = filtered%responsibility(row, component)*robust_weight(row)
               cross = cross + weight*spread(regressor, 2, q)*spread(regressor, 1, q)
               standardized = series(time)
               if (.not. fit_shift) standardized = standardized - &
                  updated%mar%shift(component)
               rhs = rhs + weight*regressor*standardized
            end do
            call symmetric_pseudoinverse(cross, inverse, info)
            if (info /= 0) then
               out%info = 3
               return
            end if
            solution = matmul(inverse, rhs)
            offset = merge(1, 0, fit_shift)
            if (fit_shift) updated%mar%shift(component) = solution(1)
            if (updated%mar%order(component) > 0) &
               updated%mar%ar(1:updated%mar%order(component), component) = &
               solution(offset + 1:offset + updated%mar%order(component))
            if (updated%seasonal_order(component) > 0) &
               updated%seasonal_ar(1:updated%seasonal_order(component), component) = &
               solution(offset + updated%mar%order(component) + 1:q)
            do row = 1, observations
               time = filtered%start + row - 1
               residual(row) = series(time) - seasonal_component_location( &
                  updated, series, time, component)
            end do
            if (effective_count > tiny(1.0_dp)) updated%mar%scale(component) = &
               max(scale_floor, sqrt(dot_product(filtered%responsibility(:, component)* &
               robust_weight, residual**2)/effective_count))
            if (fit_degrees(component)) updated%mar%degrees_of_freedom(component) = &
               estimate_student_degrees(residual/updated%mar%scale(component), &
               filtered%responsibility(:, component))
            deallocate(cross, inverse, rhs, solution, regressor, robust_weight, residual)
         end do
         updated%mar%probability = updated%mar%probability/ &
            sum(updated%mar%probability)
         candidate = updated
         candidate_likelihood = mixar_seasonal_log_likelihood(candidate, series)
         if (.not. ieee_is_finite(candidate_likelihood)) then
            out%info = 4
            return
         end if
         if (candidate_likelihood + 1.0e-10_dp < previous) &
            call damp_seasonal_step(current, updated, series, previous, candidate, &
            candidate_likelihood)
         current = candidate
         out%iterations = iteration
         if (abs(candidate_likelihood - previous) <= &
            tolerance_value*(1.0_dp + abs(previous))) then
            out%converged = .true.
            exit
         end if
         previous = candidate_likelihood
      end do
      out%model = current
      out%filter = mixar_seasonal_filter(current, series)
      out%log_likelihood = out%filter%log_likelihood
      parameters = components - 1 + components + sum(current%mar%order) + &
         sum(current%seasonal_order) + merge(components, 0, fit_shift) + &
         count(fit_degrees)
      out%aic = -2.0_dp*out%log_likelihood + 2.0_dp*real(parameters, dp)
      out%bic = -2.0_dp*out%log_likelihood + &
         log(real(size(out%filter%residual), dp))*real(parameters, dp)
   end function mixar_seasonal_fit

   pure subroutine seasonal_regressor(series, time, model, component, &
      include_shift, regressor)
      !! Build one component's ordinary and seasonal lag regression row.
      real(dp), intent(in) :: series(:) !! Observed time series.
      integer, intent(in) :: time !! Response observation index.
      type(mixar_seasonal_model_t), intent(in) :: model !! Seasonal MAR model.
      integer, intent(in) :: component !! Mixture component index.
      logical, intent(in) :: include_shift !! Include an intercept column.
      real(dp), intent(out) :: regressor(:) !! Ordinary and seasonal lag regressors.
      integer :: offset, lag

      offset = merge(1, 0, include_shift)
      if (include_shift) regressor(1) = 1.0_dp
      do lag = 1, model%mar%order(component)
         regressor(offset + lag) = series(time - lag)
      end do
      offset = offset + model%mar%order(component)
      do lag = 1, model%seasonal_order(component)
         regressor(offset + lag) = series(time - model%seasonal_period*lag)
      end do
   end subroutine seasonal_regressor

   pure real(dp) function seasonal_component_location(model, series, time, &
      component) result(location)
      !! Evaluate one seasonal MAR component's conditional location.
      type(mixar_seasonal_model_t), intent(in) :: model !! Seasonal MAR model.
      real(dp), intent(in) :: series(:) !! Observed time series.
      integer, intent(in) :: time !! Observation index.
      integer, intent(in) :: component !! Mixture component index.
      integer :: lag

      location = model%mar%shift(component)
      do lag = 1, model%mar%order(component)
         location = location + model%mar%ar(lag, component)*series(time - lag)
      end do
      do lag = 1, model%seasonal_order(component)
         location = location + model%seasonal_ar(lag, component)* &
            series(time - model%seasonal_period*lag)
      end do
   end function seasonal_component_location

   pure subroutine damp_seasonal_step(old_model, proposed_model, series, &
      old_likelihood, accepted_model, accepted_likelihood)
      !! Backtrack a seasonal MAR update until likelihood does not decrease.
      type(mixar_seasonal_model_t), intent(in) :: old_model !! Model before the update.
      type(mixar_seasonal_model_t), intent(in) :: proposed_model !! Full proposed update.
      real(dp), intent(in) :: series(:) !! Observed time series.
      real(dp), intent(in) :: old_likelihood !! Likelihood before the update.
      type(mixar_seasonal_model_t), intent(out) :: accepted_model !! Accepted model.
      real(dp), intent(out) :: accepted_likelihood !! Accepted model likelihood.
      real(dp) :: fraction
      integer :: attempt

      fraction = 0.5_dp
      do attempt = 1, 30
         accepted_model = old_model
         accepted_model%mar%probability = old_model%mar%probability + fraction* &
            (proposed_model%mar%probability - old_model%mar%probability)
         accepted_model%mar%shift = old_model%mar%shift + fraction* &
            (proposed_model%mar%shift - old_model%mar%shift)
         accepted_model%mar%scale = old_model%mar%scale + fraction* &
            (proposed_model%mar%scale - old_model%mar%scale)
         accepted_model%mar%ar = old_model%mar%ar + fraction* &
            (proposed_model%mar%ar - old_model%mar%ar)
         accepted_model%seasonal_ar = old_model%seasonal_ar + fraction* &
            (proposed_model%seasonal_ar - old_model%seasonal_ar)
         accepted_likelihood = mixar_seasonal_log_likelihood(accepted_model, series)
         if (accepted_likelihood >= old_likelihood - 1.0e-10_dp) return
         fraction = 0.5_dp*fraction
      end do
      accepted_model = old_model
      accepted_likelihood = old_likelihood
   end subroutine damp_seasonal_step

   pure logical function mixar_seasonal_is_stable(model) result(stable)
      !! Apply mixAR's second-order stability test to combined seasonal lags.
      type(mixar_seasonal_model_t), intent(in) :: model !! Seasonal MAR model.

      stable = mixar_is_stable(mixar_seasonal_expanded_model(model))
   end function mixar_seasonal_is_stable

   pure function mixar_seasonal_simulate_from_draws(model, observations, &
      initial, uniforms, standard_innovations, burnin) result(out)
      !! Simulate a seasonal MAR from supplied regime and innovation draws.
      type(mixar_seasonal_model_t), intent(in) :: model !! Seasonal MAR model.
      integer, intent(in) :: observations !! Number of retained observations.
      real(dp), intent(in) :: initial(:) !! Initial values ending at the simulation origin.
      real(dp), intent(in) :: uniforms(:) !! Regime-selection uniform draws.
      real(dp), intent(in) :: standard_innovations(:) !! Standardized innovation draws.
      integer, intent(in), optional :: burnin !! Number of discarded observations.
      type(mixar_seasonal_simulation_t) :: out
      type(mixar_simulation_t) :: simulation
      integer :: discarded

      discarded = 100
      if (present(burnin)) discarded = burnin
      simulation = mixar_simulate_from_draws(mixar_seasonal_expanded_model(model), &
         observations, initial, uniforms, standard_innovations, discarded)
      out%info = simulation%info
      if (simulation%info /= 0) return
      out%series = simulation%series
      out%regime = simulation%regime
      out%model = model
      out%burnin = discarded
   end function mixar_seasonal_simulate_from_draws

   function mixar_seasonal_simulate(model, observations, initial, burnin) &
      result(out)
      !! Simulate a seasonal MAR with the shared random generator.
      type(mixar_seasonal_model_t), intent(in) :: model !! Seasonal MAR model.
      integer, intent(in) :: observations !! Number of retained observations.
      real(dp), intent(in) :: initial(:) !! Initial values ending at the simulation origin.
      integer, intent(in), optional :: burnin !! Number of discarded observations.
      type(mixar_seasonal_simulation_t) :: out
      type(mixar_simulation_t) :: simulation
      integer :: discarded

      discarded = 100
      if (present(burnin)) discarded = burnin
      simulation = mixar_simulate(mixar_seasonal_expanded_model(model), &
         observations, initial, discarded)
      out%info = simulation%info
      if (simulation%info /= 0) return
      out%series = simulation%series
      out%regime = simulation%regime
      out%model = model
      out%burnin = discarded
   end function mixar_seasonal_simulate

   pure function mixar_seasonal_exact_forecast(model, history, horizons, &
      level, maximum_paths) result(out)
      !! Construct exact Gaussian forecasts for an additive seasonal MAR.
      type(mixar_seasonal_model_t), intent(in) :: model !! Seasonal MAR model.
      real(dp), intent(in) :: history(:) !! Conditioning history ending with the latest value.
      integer, intent(in) :: horizons !! Positive maximum forecast horizon.
      real(dp), intent(in), optional :: level !! Central prediction interval probability.
      integer, intent(in), optional :: maximum_paths !! Maximum regime paths at any horizon.
      type(mixar_exact_forecast_t) :: out
      real(dp) :: interval_level
      integer :: path_limit

      interval_level = 0.95_dp
      path_limit = 1000000
      if (present(level)) interval_level = level
      if (present(maximum_paths)) path_limit = maximum_paths
      out = mixar_exact_forecast(mixar_seasonal_expanded_model(model), history, &
         horizons, interval_level, path_limit)
   end function mixar_seasonal_exact_forecast

   pure function mixar_seasonal_forecast_from_draws(model, history, uniforms, &
      standard_innovations, level) result(out)
      !! Forecast a seasonal MAR from supplied regime and innovation draws.
      type(mixar_seasonal_model_t), intent(in) :: model !! Seasonal MAR model.
      real(dp), intent(in) :: history(:) !! Conditioning history ending with the latest value.
      real(dp), intent(in) :: uniforms(:, :) !! Horizon-by-simulation regime uniform draws.
      real(dp), intent(in) :: standard_innovations(:, :) !! Horizon-by-simulation standardized innovations.
      real(dp), intent(in), optional :: level !! Central empirical interval probability.
      type(mixar_forecast_sample_t) :: out
      real(dp) :: interval_level

      interval_level = 0.95_dp
      if (present(level)) interval_level = level
      out = mixar_forecast_from_draws(mixar_seasonal_expanded_model(model), &
         history, uniforms, standard_innovations, interval_level)
   end function mixar_seasonal_forecast_from_draws

   function mixar_seasonal_forecast(model, history, horizons, simulations, &
      level) result(out)
      !! Forecast a seasonal MAR with the shared random generator.
      type(mixar_seasonal_model_t), intent(in) :: model !! Seasonal MAR model.
      real(dp), intent(in) :: history(:) !! Conditioning history ending with the latest value.
      integer, intent(in) :: horizons !! Positive forecast horizon count.
      integer, intent(in) :: simulations !! Positive number of simulated paths.
      real(dp), intent(in), optional :: level !! Central empirical interval probability.
      type(mixar_forecast_sample_t) :: out
      real(dp) :: interval_level

      interval_level = 0.95_dp
      if (present(level)) interval_level = level
      out = mixar_forecast(mixar_seasonal_expanded_model(model), history, &
         horizons, simulations, interval_level)
   end function mixar_seasonal_forecast

   pure function mixar_diagnose(model, series, lag) result(out)
      !! Compute residual, quantile, portmanteau, and classification diagnostics.
      type(mixar_model_t), intent(in) :: model !! Gaussian or Student-t MAR model.
      real(dp), intent(in) :: series(:) !! Observed time series.
      integer, intent(in), optional :: lag !! Maximum ACF and portmanteau lag.
      type(mixar_diagnostics_t) :: out
      type(mixar_filter_t) :: filtered
      real(dp) :: probability, entropy_term
      integer :: selected_lag, row, component, time, selected

      filtered = mixar_filter(model, series)
      if (filtered%info /= 0) then
         out%info = 1
         return
      end if
      selected_lag = min(20, size(filtered%residual) - 1)
      if (present(lag)) selected_lag = lag
      if (selected_lag < 1 .or. selected_lag >= size(filtered%residual)) then
         out%info = 1
         return
      end if
      out%residual = filtered%residual
      allocate(out%standardized_residual(size(filtered%residual)))
      allocate(out%uniform_residual(size(filtered%residual)))
      allocate(out%quantile_residual(size(filtered%residual)))
      allocate(out%classified_residual(size(filtered%residual)))
      allocate(out%classification(size(filtered%residual)))
      allocate(out%maximum_posterior(size(filtered%residual)))
      out%posterior_probability = filtered%responsibility
      out%standardized_residual = filtered%residual/sqrt(filtered%variance)
      out%mean_entropy = 0.0_dp
      do row = 1, size(filtered%residual)
         time = filtered%start + row - 1
         probability = 0.0_dp
         do component = 1, size(model%probability)
            probability = probability + model%probability(component)* &
               mixar_standard_cdf((series(time) - filtered%location(row, component))/ &
               model%scale(component), model%degrees_of_freedom(component))
            if (filtered%responsibility(row, component) > 0.0_dp) then
               entropy_term = -filtered%responsibility(row, component)* &
                  log(filtered%responsibility(row, component))
               out%mean_entropy = out%mean_entropy + entropy_term
            end if
         end do
         probability = max(epsilon(1.0_dp), min(1.0_dp - epsilon(1.0_dp), probability))
         out%uniform_residual(row) = probability
         out%quantile_residual(row) = normal_quantile(probability)
         selected = maxloc(filtered%responsibility(row, :), dim=1)
         out%classification(row) = selected
         out%maximum_posterior(row) = filtered%responsibility(row, selected)
         out%classified_residual(row) = (series(time) - &
            filtered%location(row, selected))/model%scale(selected)
      end do
      out%mean_entropy = out%mean_entropy/real(size(filtered%residual), dp)
      out%residual_acf = acf_values(out%standardized_residual, selected_lag)
      out%squared_residual_acf = acf_values(out%standardized_residual**2, selected_lag)
      out%residual_test = weighted_box_test(out%standardized_residual, selected_lag, &
         box_test_ljung_box, 0, .false., residual_raw)
      out%squared_residual_test = weighted_box_test(out%standardized_residual, &
         selected_lag, box_test_ljung_box, 0, .false., residual_squared)
      out%bic = mixar_bic(model, series)
   end function mixar_diagnose

   pure real(dp) function mixar_bic(model, series, estimate_shift, &
      estimate_degrees) result(value)
      !! Compute mixAR's conditional BIC with explicit fitted-parameter flags.
      type(mixar_model_t), intent(in) :: model !! Gaussian or Student-t MAR model.
      real(dp), intent(in) :: series(:) !! Observed time series.
      logical, intent(in), optional :: estimate_shift !! Count component shifts as fitted.
      logical, intent(in), optional :: estimate_degrees(:) !! Count selected Student-t degrees.
      type(mixar_filter_t) :: filtered
      integer :: parameters
      logical :: count_shift

      value = huge(1.0_dp)
      count_shift = .true.
      if (present(estimate_shift)) count_shift = estimate_shift
      filtered = mixar_filter(model, series)
      if (filtered%info /= 0) return
      parameters = size(model%probability) - 1 + size(model%scale) + &
         sum(model%order) + merge(size(model%shift), 0, count_shift)
      if (present(estimate_degrees)) then
         if (size(estimate_degrees) /= size(model%probability)) return
         parameters = parameters + count(estimate_degrees)
      end if
      value = -2.0_dp*filtered%log_likelihood + &
         real(parameters, dp)*log(real(size(filtered%residual), dp))
   end function mixar_bic

   pure function mixar_select_bic(models, series, estimate_shift) result(out)
      !! Compare candidate MAR models by conditional BIC.
      type(mixar_model_t), intent(in) :: models(:) !! Candidate MAR models.
      real(dp), intent(in) :: series(:) !! Common observed time series.
      logical, intent(in), optional :: estimate_shift !! Count component shifts as fitted.
      type(mixar_bic_selection_t) :: out
      logical :: count_shift
      integer :: candidate

      count_shift = .true.
      if (present(estimate_shift)) count_shift = estimate_shift
      if (size(models) < 1) then
         out%info = 1
         return
      end if
      allocate(out%bic(size(models)))
      do candidate = 1, size(models)
         out%bic(candidate) = mixar_bic(models(candidate), series, count_shift)
      end do
      out%selected = minloc(out%bic, dim=1)
      if (.not. ieee_is_finite(out%bic(out%selected))) out%info = 1
   end function mixar_select_bic

   pure function mixar_observed_inference(model, series, estimate_shift, &
      estimate_degrees, level) result(out)
      !! Estimate MAR covariance and confidence intervals from the observed Hessian.
      type(mixar_model_t), intent(in) :: model !! Fitted Gaussian or Student-t MAR model.
      real(dp), intent(in) :: series(:) !! Observed time series.
      logical, intent(in), optional :: estimate_shift !! Include component shifts in inference.
      logical, intent(in), optional :: estimate_degrees(:) !! Include selected Student-t degrees.
      real(dp), intent(in), optional :: level !! Central confidence level.
      type(mixar_inference_t) :: out
      real(dp), allocatable :: parameters(:), steps(:)
      logical, allocatable :: fit_degrees(:)
      logical :: fit_shift
      integer :: components

      fit_shift = .true.
      if (present(estimate_shift)) fit_shift = estimate_shift
      components = size(model%probability)
      allocate(fit_degrees(components), source=.false.)
      if (present(estimate_degrees)) then
         if (size(estimate_degrees) /= components) then
            out%info = 1
            return
         end if
         fit_degrees = estimate_degrees
      end if
      if (.not. valid_model(model) .or. any(fit_degrees .and. &
         model%degrees_of_freedom <= 2.0_dp)) then
         out%info = 1
         return
      end if
      call pack_base_parameters(model, fit_shift, fit_degrees, parameters, &
         steps, out%parameter)
      call finite_difference_inference(parameters, steps, objective, out)
      if (present(level)) out%level = level
      call finish_inference(out)
      call probability_errors(out, components)

   contains

      pure real(dp) function objective(candidate_parameters) result(value)
         !! Evaluate the MAR likelihood at one numerical-Hessian parameter vector.
         real(dp), intent(in) :: candidate_parameters(:) !! Candidate independent parameters.
         type(mixar_model_t) :: candidate

         candidate = unpack_base_parameters(model, candidate_parameters, &
            fit_shift, fit_degrees)
         value = mixar_log_likelihood(candidate, series)
      end function objective

   end function mixar_observed_inference

   pure subroutine pack_base_parameters(model, fit_shift, fit_degrees, &
      parameters, steps, labels)
      !! Pack independent MAR parameters and stable finite-difference steps.
      type(mixar_model_t), intent(in) :: model !! Fitted MAR model.
      logical, intent(in) :: fit_shift !! Include component shifts.
      logical, intent(in) :: fit_degrees(:) !! Include selected Student-t degrees.
      real(dp), allocatable, intent(out) :: parameters(:) !! Packed natural parameters.
      real(dp), allocatable, intent(out) :: steps(:) !! Constraint-aware differentiation steps.
      character(len=48), allocatable, intent(out) :: labels(:) !! Parameter labels.
      real(dp) :: probability_step
      integer :: count_value, position, component, lag

      count_value = size(model%probability) - 1 + sum(model%order) + &
         size(model%scale) + merge(size(model%shift), 0, fit_shift) + &
         count(fit_degrees)
      allocate(parameters(count_value), steps(count_value), labels(count_value))
      position = 0
      probability_step = max(1.0e-7_dp, 0.05_dp*minval(model%probability))
      do component = 1, size(model%probability) - 1
         position = position + 1
         parameters(position) = model%probability(component)
         steps(position) = probability_step
         write(labels(position), '(a,i0)') 'probability_', component
      end do
      do component = 1, size(model%probability)
         if (fit_shift) then
            position = position + 1
            parameters(position) = model%shift(component)
            steps(position) = 5.0e-4_dp*(1.0_dp + abs(parameters(position)))
            write(labels(position), '(a,i0)') 'shift_', component
         end if
         do lag = 1, model%order(component)
            position = position + 1
            parameters(position) = model%ar(lag, component)
            steps(position) = 5.0e-4_dp*(1.0_dp + abs(parameters(position)))
            write(labels(position), '(a,i0,a,i0)') 'ar_', component, '_', lag
         end do
         position = position + 1
         parameters(position) = model%scale(component)
         steps(position) = min(5.0e-4_dp*(1.0_dp + parameters(position)), &
            0.1_dp*parameters(position))
         write(labels(position), '(a,i0)') 'scale_', component
         if (fit_degrees(component)) then
            position = position + 1
            parameters(position) = model%degrees_of_freedom(component)
            steps(position) = min(5.0e-4_dp*(1.0_dp + parameters(position)), &
               0.1_dp*(parameters(position) - 2.0_dp))
            write(labels(position), '(a,i0)') 'degrees_', component
         end if
      end do
   end subroutine pack_base_parameters

   pure function unpack_base_parameters(template, parameters, fit_shift, &
      fit_degrees) result(model)
      !! Restore a MAR model from its independent natural parameter vector.
      type(mixar_model_t), intent(in) :: template !! Model supplying fixed values and dimensions.
      real(dp), intent(in) :: parameters(:) !! Packed independent parameters.
      logical, intent(in) :: fit_shift !! Include component shifts.
      logical, intent(in) :: fit_degrees(:) !! Include selected Student-t degrees.
      type(mixar_model_t) :: model
      integer :: position, component, lag, components

      model = template
      components = size(model%probability)
      position = components - 1
      if (components > 1) model%probability(:components - 1) = &
         parameters(:components - 1)
      model%probability(components) = 1.0_dp - &
         sum(model%probability(:components - 1))
      do component = 1, components
         if (fit_shift) then
            position = position + 1
            model%shift(component) = parameters(position)
         end if
         do lag = 1, model%order(component)
            position = position + 1
            model%ar(lag, component) = parameters(position)
         end do
         position = position + 1
         model%scale(component) = parameters(position)
         if (fit_degrees(component)) then
            position = position + 1
            model%degrees_of_freedom(component) = parameters(position)
         end if
      end do
      if (.not. valid_model(model)) model%info = 1
   end function unpack_base_parameters

   pure subroutine finite_difference_inference(parameters, steps, objective, out)
      !! Numerically differentiate a scalar likelihood and invert observed information.
      real(dp), intent(in) :: parameters(:) !! Parameter vector at the fitted model.
      real(dp), intent(in) :: steps(:) !! Positive central-difference steps.
      procedure(scalar_objective_t) :: objective
      type(mixar_inference_t), intent(inout) :: out !! Inference result receiving matrices.
      real(dp), allocatable :: plus(:), minus(:), first_second(:), first_minus(:)
      real(dp), allocatable :: minus_second(:), both_minus(:)
      real(dp) :: center, f_plus, f_minus
      integer :: first, second, n, info

      n = size(parameters)
      out%estimate = parameters
      allocate(out%observed_information(n, n), out%covariance(n, n))
      out%observed_information = 0.0_dp
      center = objective(parameters)
      if (.not. ieee_is_finite(center)) then
         out%info = 2
         return
      end if
      do first = 1, n
         plus = parameters
         minus = parameters
         plus(first) = plus(first) + steps(first)
         minus(first) = minus(first) - steps(first)
         f_plus = objective(plus)
         f_minus = objective(minus)
         out%observed_information(first, first) = &
            -(f_plus - 2.0_dp*center + f_minus)/steps(first)**2
         do second = first + 1, n
            first_second = parameters
            first_minus = parameters
            minus_second = parameters
            both_minus = parameters
            first_second(first) = first_second(first) + steps(first)
            first_second(second) = first_second(second) + steps(second)
            first_minus(first) = first_minus(first) + steps(first)
            first_minus(second) = first_minus(second) - steps(second)
            minus_second(first) = minus_second(first) - steps(first)
            minus_second(second) = minus_second(second) + steps(second)
            both_minus(first) = both_minus(first) - steps(first)
            both_minus(second) = both_minus(second) - steps(second)
            out%observed_information(first, second) = -(objective(first_second) - &
               objective(first_minus) - objective(minus_second) + &
               objective(both_minus))/(4.0_dp*steps(first)*steps(second))
            out%observed_information(second, first) = &
               out%observed_information(first, second)
         end do
      end do
      call symmetric_pseudoinverse(out%observed_information, out%covariance, info)
      if (info /= 0) out%info = 3
   end subroutine finite_difference_inference

   pure subroutine finish_inference(out)
      !! Derive standard errors and central normal confidence intervals.
      type(mixar_inference_t), intent(inout) :: out !! Inference result with covariance.
      real(dp) :: critical
      integer :: parameter

      if (out%info /= 0 .or. out%level <= 0.0_dp .or. out%level >= 1.0_dp) then
         if (out%info == 0) out%info = 1
         return
      end if
      allocate(out%standard_error(size(out%estimate)))
      allocate(out%lower(size(out%estimate)), out%upper(size(out%estimate)))
      do parameter = 1, size(out%estimate)
         out%standard_error(parameter) = sqrt(max(0.0_dp, &
            out%covariance(parameter, parameter)))
      end do
      critical = normal_quantile(0.5_dp*(1.0_dp + out%level))
      out%lower = out%estimate - critical*out%standard_error
      out%upper = out%estimate + critical*out%standard_error
   end subroutine finish_inference

   pure subroutine probability_errors(out, components)
      !! Recover standard errors for all weights including the dependent final weight.
      type(mixar_inference_t), intent(inout) :: out !! Inference result with covariance.
      integer, intent(in) :: components !! Number of mixture components.
      integer :: independent

      allocate(out%probability_standard_error(components), source=0.0_dp)
      independent = components - 1
      if (out%info /= 0 .or. independent == 0) return
      out%probability_standard_error(:independent) = &
         out%standard_error(:independent)
      out%probability_standard_error(components) = sqrt(max(0.0_dp, &
         sum(out%covariance(:independent, :independent))))
   end subroutine probability_errors

   pure function mixar_seasonal_observed_inference(model, series, &
      estimate_shift, estimate_degrees, level) result(out)
      !! Estimate seasonal MAR covariance and intervals from the observed Hessian.
      type(mixar_seasonal_model_t), intent(in) :: model !! Fitted seasonal MAR model.
      real(dp), intent(in) :: series(:) !! Observed time series.
      logical, intent(in), optional :: estimate_shift !! Include component shifts.
      logical, intent(in), optional :: estimate_degrees(:) !! Include selected Student-t degrees.
      real(dp), intent(in), optional :: level !! Central confidence level.
      type(mixar_inference_t) :: out
      real(dp), allocatable :: parameters(:), steps(:)
      logical, allocatable :: fit_degrees(:)
      logical :: fit_shift
      integer :: components

      fit_shift = .true.
      if (present(estimate_shift)) fit_shift = estimate_shift
      components = size(model%mar%probability)
      allocate(fit_degrees(components), source=.false.)
      if (present(estimate_degrees)) then
         if (size(estimate_degrees) /= components) then
            out%info = 1
            return
         end if
         fit_degrees = estimate_degrees
      end if
      if (.not. valid_seasonal_model(model) .or. any(fit_degrees .and. &
         model%mar%degrees_of_freedom <= 2.0_dp)) then
         out%info = 1
         return
      end if
      call pack_seasonal_parameters(model, fit_shift, fit_degrees, parameters, &
         steps, out%parameter)
      call finite_difference_inference(parameters, steps, objective, out)
      if (present(level)) out%level = level
      call finish_inference(out)
      call probability_errors(out, components)

   contains

      pure real(dp) function objective(candidate_parameters) result(value)
         !! Evaluate seasonal MAR likelihood for numerical differentiation.
         real(dp), intent(in) :: candidate_parameters(:) !! Candidate independent parameters.
         type(mixar_seasonal_model_t) :: candidate

         candidate = unpack_seasonal_parameters(model, candidate_parameters, &
            fit_shift, fit_degrees)
         value = mixar_seasonal_log_likelihood(candidate, series)
      end function objective

   end function mixar_seasonal_observed_inference

   pure subroutine pack_seasonal_parameters(model, fit_shift, fit_degrees, &
      parameters, steps, labels)
      !! Pack seasonal MAR parameters and constraint-aware difference steps.
      type(mixar_seasonal_model_t), intent(in) :: model !! Seasonal MAR model.
      logical, intent(in) :: fit_shift !! Include component shifts.
      logical, intent(in) :: fit_degrees(:) !! Include selected Student-t degrees.
      real(dp), allocatable, intent(out) :: parameters(:) !! Packed parameters.
      real(dp), allocatable, intent(out) :: steps(:) !! Differentiation steps.
      character(len=48), allocatable, intent(out) :: labels(:) !! Parameter labels.
      real(dp) :: probability_step
      integer :: count_value, position, component, lag

      count_value = size(model%mar%probability) - 1 + sum(model%mar%order) + &
         sum(model%seasonal_order) + size(model%mar%scale) + &
         merge(size(model%mar%shift), 0, fit_shift) + count(fit_degrees)
      allocate(parameters(count_value), steps(count_value), labels(count_value))
      position = 0
      probability_step = max(1.0e-7_dp, 0.05_dp*minval(model%mar%probability))
      do component = 1, size(model%mar%probability) - 1
         position = position + 1
         parameters(position) = model%mar%probability(component)
         steps(position) = probability_step
         write(labels(position), '(a,i0)') 'probability_', component
      end do
      do component = 1, size(model%mar%probability)
         if (fit_shift) then
            position = position + 1
            parameters(position) = model%mar%shift(component)
            steps(position) = 5.0e-4_dp*(1.0_dp + abs(parameters(position)))
            write(labels(position), '(a,i0)') 'shift_', component
         end if
         do lag = 1, model%mar%order(component)
            position = position + 1
            parameters(position) = model%mar%ar(lag, component)
            steps(position) = 5.0e-4_dp*(1.0_dp + abs(parameters(position)))
            write(labels(position), '(a,i0,a,i0)') 'ar_', component, '_', lag
         end do
         do lag = 1, model%seasonal_order(component)
            position = position + 1
            parameters(position) = model%seasonal_ar(lag, component)
            steps(position) = 5.0e-4_dp*(1.0_dp + abs(parameters(position)))
            write(labels(position), '(a,i0,a,i0)') 'seasonal_ar_', component, '_', lag
         end do
         position = position + 1
         parameters(position) = model%mar%scale(component)
         steps(position) = min(5.0e-4_dp*(1.0_dp + parameters(position)), &
            0.1_dp*parameters(position))
         write(labels(position), '(a,i0)') 'scale_', component
         if (fit_degrees(component)) then
            position = position + 1
            parameters(position) = model%mar%degrees_of_freedom(component)
            steps(position) = min(5.0e-4_dp*(1.0_dp + parameters(position)), &
               0.1_dp*(parameters(position) - 2.0_dp))
            write(labels(position), '(a,i0)') 'degrees_', component
         end if
      end do
   end subroutine pack_seasonal_parameters

   pure function unpack_seasonal_parameters(template, parameters, fit_shift, &
      fit_degrees) result(model)
      !! Restore a seasonal MAR model from independent parameters.
      type(mixar_seasonal_model_t), intent(in) :: template !! Template seasonal model.
      real(dp), intent(in) :: parameters(:) !! Packed independent parameters.
      logical, intent(in) :: fit_shift !! Include component shifts.
      logical, intent(in) :: fit_degrees(:) !! Include selected Student-t degrees.
      type(mixar_seasonal_model_t) :: model
      integer :: position, component, lag, components

      model = template
      components = size(model%mar%probability)
      position = components - 1
      if (components > 1) model%mar%probability(:components - 1) = &
         parameters(:components - 1)
      model%mar%probability(components) = 1.0_dp - &
         sum(model%mar%probability(:components - 1))
      do component = 1, components
         if (fit_shift) then
            position = position + 1
            model%mar%shift(component) = parameters(position)
         end if
         do lag = 1, model%mar%order(component)
            position = position + 1
            model%mar%ar(lag, component) = parameters(position)
         end do
         do lag = 1, model%seasonal_order(component)
            position = position + 1
            model%seasonal_ar(lag, component) = parameters(position)
         end do
         position = position + 1
         model%mar%scale(component) = parameters(position)
         if (fit_degrees(component)) then
            position = position + 1
            model%mar%degrees_of_freedom(component) = parameters(position)
         end if
      end do
      if (.not. valid_seasonal_model(model)) model%info = 1
   end function unpack_seasonal_parameters

   pure function mixar_regression_observed_inference(model, response, regressors, &
      estimate_shift, estimate_degrees, level) result(out)
      !! Estimate MAR regression covariance and intervals from the observed Hessian.
      type(mixar_regression_model_t), intent(in) :: model !! Fitted MAR regression model.
      real(dp), intent(in) :: response(:) !! Response time series.
      real(dp), intent(in) :: regressors(:, :) !! Observation-by-regressor design matrix.
      logical, intent(in), optional :: estimate_shift !! Include component shifts.
      logical, intent(in), optional :: estimate_degrees(:) !! Include selected Student-t degrees.
      real(dp), intent(in), optional :: level !! Central confidence level.
      type(mixar_inference_t) :: out
      real(dp), allocatable :: parameters(:), steps(:)
      logical, allocatable :: fit_degrees(:)
      logical :: fit_shift
      integer :: components

      fit_shift = .true.
      if (present(estimate_shift)) fit_shift = estimate_shift
      components = size(model%mar%probability)
      allocate(fit_degrees(components), source=.false.)
      if (present(estimate_degrees)) then
         if (size(estimate_degrees) /= components) then
            out%info = 1
            return
         end if
         fit_degrees = estimate_degrees
      end if
      if (.not. valid_regression_model(model) .or. &
         size(regressors, 1) /= size(response) .or. &
         size(regressors, 2) /= size(model%coefficient, 1) .or. &
         any(fit_degrees .and. model%mar%degrees_of_freedom <= 2.0_dp)) then
         out%info = 1
         return
      end if
      call pack_regression_parameters(model, fit_shift, fit_degrees, parameters, &
         steps, out%parameter)
      call finite_difference_inference(parameters, steps, objective, out)
      if (present(level)) out%level = level
      call finish_inference(out)
      call probability_errors(out, components)

   contains

      pure real(dp) function objective(candidate_parameters) result(value)
         !! Evaluate MAR regression likelihood for numerical differentiation.
         real(dp), intent(in) :: candidate_parameters(:) !! Candidate independent parameters.
         type(mixar_regression_model_t) :: candidate

         candidate = unpack_regression_parameters(model, candidate_parameters, &
            fit_shift, fit_degrees)
         value = mixar_regression_log_likelihood(candidate, response, regressors)
      end function objective

   end function mixar_regression_observed_inference

   pure subroutine pack_regression_parameters(model, fit_shift, fit_degrees, &
      parameters, steps, labels)
      !! Pack MAR regression and coefficient parameters for differentiation.
      type(mixar_regression_model_t), intent(in) :: model !! MAR regression model.
      logical, intent(in) :: fit_shift !! Include component shifts.
      logical, intent(in) :: fit_degrees(:) !! Include selected Student-t degrees.
      real(dp), allocatable, intent(out) :: parameters(:) !! Packed parameters.
      real(dp), allocatable, intent(out) :: steps(:) !! Differentiation steps.
      character(len=48), allocatable, intent(out) :: labels(:) !! Parameter labels.
      real(dp), allocatable :: base_parameters(:), base_steps(:)
      character(len=48), allocatable :: base_labels(:)
      integer :: base_count, position, column, regressor

      call pack_base_parameters(model%mar, fit_shift, fit_degrees, base_parameters, &
         base_steps, base_labels)
      base_count = size(base_parameters)
      allocate(parameters(base_count + size(model%coefficient)))
      allocate(steps(size(parameters)), labels(size(parameters)))
      parameters(:base_count) = base_parameters
      steps(:base_count) = base_steps
      labels(:base_count) = base_labels
      position = base_count
      do column = 1, size(model%coefficient, 2)
         do regressor = 1, size(model%coefficient, 1)
            position = position + 1
            parameters(position) = model%coefficient(regressor, column)
            steps(position) = 5.0e-4_dp*(1.0_dp + abs(parameters(position)))
            write(labels(position), '(a,i0,a,i0)') 'regression_', column, '_', regressor
         end do
      end do
   end subroutine pack_regression_parameters

   pure function unpack_regression_parameters(template, parameters, fit_shift, &
      fit_degrees) result(model)
      !! Restore a MAR regression model from independent parameters.
      type(mixar_regression_model_t), intent(in) :: template !! Template regression model.
      real(dp), intent(in) :: parameters(:) !! Packed independent parameters.
      logical, intent(in) :: fit_shift !! Include component shifts.
      logical, intent(in) :: fit_degrees(:) !! Include selected Student-t degrees.
      type(mixar_regression_model_t) :: model
      integer :: base_count, position, column, regressor

      model = template
      base_count = size(model%mar%probability) - 1 + sum(model%mar%order) + &
         size(model%mar%scale) + merge(size(model%mar%shift), 0, fit_shift) + &
         count(fit_degrees)
      model%mar = unpack_base_parameters(template%mar, parameters(:base_count), &
         fit_shift, fit_degrees)
      position = base_count
      do column = 1, size(model%coefficient, 2)
         do regressor = 1, size(model%coefficient, 1)
            position = position + 1
            model%coefficient(regressor, column) = parameters(position)
         end do
      end do
      if (.not. valid_regression_model(model)) model%info = 1
   end function unpack_regression_parameters

   pure function mixar_exact_forecast(model, history, horizons, level, &
      maximum_paths) result(out)
      !! Construct exact multi-step Gaussian path-mixture forecast distributions.
      type(mixar_model_t), intent(in) :: model !! Gaussian MAR model.
      real(dp), intent(in) :: history(:) !! Conditioning history ending with the latest value.
      integer, intent(in) :: horizons !! Positive maximum forecast horizon.
      real(dp), intent(in), optional :: level !! Central prediction interval probability.
      integer, intent(in), optional :: maximum_paths !! Maximum regime paths at any horizon.
      type(mixar_exact_forecast_t) :: out
      integer :: horizon, paths, path_limit
      real(dp) :: interval_level, tail_probability

      interval_level = 0.95_dp
      path_limit = 1000000
      if (present(level)) interval_level = level
      if (present(maximum_paths)) path_limit = maximum_paths
      if (.not. valid_model(model) .or. any(model%degrees_of_freedom /= 0.0_dp) .or. &
         horizons < 1 .or. size(history) < maxval(model%order) .or. &
         interval_level <= 0.0_dp .or. interval_level >= 1.0_dp .or. &
         path_limit < 1) then
         out%info = 1
         return
      end if
      paths = 1
      do horizon = 1, horizons
         if (paths > path_limit/size(model%probability)) then
            out%info = 2
            return
         end if
         paths = paths*size(model%probability)
      end do
      allocate(out%distribution(horizons), out%mean(horizons))
      allocate(out%variance(horizons), out%lower(horizons), out%upper(horizons))
      out%level = interval_level
      tail_probability = 0.5_dp*(1.0_dp - interval_level)
      do horizon = 1, horizons
         out%distribution(horizon) = exact_horizon_distribution(model, history, horizon)
         out%mean(horizon) = dot_product(out%distribution(horizon)%probability, &
            out%distribution(horizon)%location)
         out%variance(horizon) = dot_product(out%distribution(horizon)%probability, &
            out%distribution(horizon)%scale**2 + &
            out%distribution(horizon)%location**2) - out%mean(horizon)**2
         out%lower(horizon) = mixar_predictive_quantile( &
            out%distribution(horizon), tail_probability)
         out%upper(horizon) = mixar_predictive_quantile( &
            out%distribution(horizon), 1.0_dp - tail_probability)
      end do
   end function mixar_exact_forecast

   pure function exact_horizon_distribution(model, history, horizon) result(out)
      !! Propagate every regime sequence into a Gaussian forecast component.
      type(mixar_model_t), intent(in) :: model !! Gaussian MAR model.
      real(dp), intent(in) :: history(:) !! Conditioning history.
      integer, intent(in) :: horizon !! Forecast horizon.
      type(mixar_predictive_distribution_t) :: out
      real(dp), allocatable :: state(:), next_state(:)
      real(dp), allocatable :: noise(:, :), next_noise(:, :), coefficient(:)
      integer, allocatable :: regimes(:)
      integer :: p, paths, path, code, step, regime, row

      p = max(1, maxval(model%order))
      paths = size(model%probability)**horizon
      allocate(out%probability(paths), out%location(paths), out%scale(paths))
      allocate(out%regime_path(horizon, paths), state(p), next_state(p))
      allocate(noise(p, horizon), next_noise(p, horizon), coefficient(p))
      allocate(regimes(horizon))
      out%horizon = horizon
      do path = 1, paths
         code = path - 1
         do step = horizon, 1, -1
            regimes(step) = modulo(code, size(model%probability)) + 1
            code = code/size(model%probability)
         end do
         out%regime_path(:, path) = regimes
         state = 0.0_dp
         do row = 1, maxval(model%order)
            state(row) = history(size(history) - row + 1)
         end do
         noise = 0.0_dp
         out%probability(path) = 1.0_dp
         do step = 1, horizon
            regime = regimes(step)
            coefficient = 0.0_dp
            if (model%order(regime) > 0) coefficient(1:model%order(regime)) = &
               model%ar(1:model%order(regime), regime)
            next_state = 0.0_dp
            next_state(1) = model%shift(regime) + dot_product(coefficient, state)
            if (p > 1) next_state(2:p) = state(1:p - 1)
            next_noise = 0.0_dp
            next_noise(1, :) = matmul(coefficient, noise)
            next_noise(1, step) = next_noise(1, step) + model%scale(regime)
            if (p > 1) next_noise(2:p, :) = noise(1:p - 1, :)
            state = next_state
            noise = next_noise
            out%probability(path) = out%probability(path)*model%probability(regime)
         end do
         out%location(path) = state(1)
         out%scale(path) = sqrt(sum(noise(1, :)**2))
      end do
      out%probability = out%probability/sum(out%probability)
   end function exact_horizon_distribution

   pure real(dp) function mixar_predictive_density(distribution, value) &
      result(density)
      !! Evaluate an exact Gaussian path-mixture predictive density.
      type(mixar_predictive_distribution_t), intent(in) :: distribution !! Exact horizon distribution.
      real(dp), intent(in) :: value !! Forecast value.
      real(dp), allocatable :: standardized(:)

      density = 0.0_dp
      if (.not. allocated(distribution%probability) .or. &
         .not. allocated(distribution%location) .or. &
         .not. allocated(distribution%scale)) return
      standardized = (value - distribution%location)/distribution%scale
      density = sum(distribution%probability* &
         exp(-0.5_dp*(log_two_pi + standardized**2))/distribution%scale)
   end function mixar_predictive_density

   pure real(dp) function mixar_predictive_cdf(distribution, value) &
      result(probability)
      !! Evaluate an exact Gaussian path-mixture predictive distribution function.
      type(mixar_predictive_distribution_t), intent(in) :: distribution !! Exact horizon distribution.
      real(dp), intent(in) :: value !! Forecast value.

      probability = 0.0_dp
      if (.not. allocated(distribution%probability) .or. &
         .not. allocated(distribution%location) .or. &
         .not. allocated(distribution%scale)) return
      probability = sum(0.5_dp*distribution%probability* &
         erfc(-(value - distribution%location)/(sqrt(2.0_dp)*distribution%scale)))
   end function mixar_predictive_cdf

   pure real(dp) function mixar_predictive_quantile(distribution, probability) &
      result(value)
      !! Invert an exact Gaussian path-mixture predictive distribution.
      type(mixar_predictive_distribution_t), intent(in) :: distribution !! Exact horizon distribution.
      real(dp), intent(in) :: probability !! Cumulative probability in the unit interval.
      real(dp) :: left, right, middle, center, width
      integer :: iteration

      value = 0.0_dp
      if (.not. allocated(distribution%probability) .or. probability < 0.0_dp .or. &
         probability > 1.0_dp) return
      center = dot_product(distribution%probability, distribution%location)
      width = max(1.0_dp, maxval(abs(distribution%location - center) + &
         8.0_dp*distribution%scale))
      left = center - width
      right = center + width
      if (probability == 0.0_dp) then
         value = -huge(1.0_dp)
         return
      else if (probability == 1.0_dp) then
         value = huge(1.0_dp)
         return
      end if
      do iteration = 1, 120
         middle = 0.5_dp*(left + right)
         if (mixar_predictive_cdf(distribution, middle) < probability) then
            left = middle
         else
            right = middle
         end if
      end do
      value = 0.5_dp*(left + right)
   end function mixar_predictive_quantile

   pure function mixar_forecast_from_draws(model, history, uniforms, &
      standard_innovations, level) result(out)
      !! Simulate forecast paths from supplied regime uniforms and standardized innovations.
      type(mixar_model_t), intent(in) :: model !! Gaussian or Student-t MAR model.
      real(dp), intent(in) :: history(:) !! Conditioning history ending with the latest value.
      real(dp), intent(in) :: uniforms(:, :) !! Horizon-by-simulation regime uniform draws.
      real(dp), intent(in) :: standard_innovations(:, :) !! Horizon-by-simulation standardized innovations.
      real(dp), intent(in), optional :: level !! Central empirical interval probability.
      type(mixar_forecast_sample_t) :: out
      type(mixar_simulation_t) :: simulation
      real(dp), allocatable :: ordered(:)
      real(dp) :: interval_level, tail_probability
      integer :: horizons, simulations, sample, horizon

      interval_level = 0.95_dp
      if (present(level)) interval_level = level
      if (.not. valid_model(model) .or. any(shape(uniforms) /= shape(standard_innovations)) .or. &
         size(uniforms, 1) < 1 .or. size(uniforms, 2) < 1 .or. &
         size(history) < maxval(model%order) .or. interval_level <= 0.0_dp .or. &
         interval_level >= 1.0_dp) then
         out%info = 1
         return
      end if
      horizons = size(uniforms, 1)
      simulations = size(uniforms, 2)
      allocate(out%paths(horizons, simulations), out%mean(horizons))
      allocate(out%variance(horizons), out%lower(horizons), out%upper(horizons))
      do sample = 1, simulations
         simulation = mixar_simulate_from_draws(model, horizons, history, &
            uniforms(:, sample), standard_innovations(:, sample), 0)
         if (simulation%info /= 0) then
            out%info = 2
            return
         end if
         out%paths(:, sample) = simulation%series
      end do
      out%level = interval_level
      tail_probability = 0.5_dp*(1.0_dp - interval_level)
      do horizon = 1, horizons
         ordered = sorted(out%paths(horizon, :))
         out%mean(horizon) = sum(ordered)/real(simulations, dp)
         out%variance(horizon) = variance(ordered)
         out%lower(horizon) = quantile(ordered, tail_probability)
         out%upper(horizon) = quantile(ordered, 1.0_dp - tail_probability)
      end do
   end function mixar_forecast_from_draws

   function mixar_forecast(model, history, horizons, simulations, level) &
      result(out)
      !! Simulate Gaussian or Student-t MAR forecast paths with the shared RNG.
      type(mixar_model_t), intent(in) :: model !! Gaussian or Student-t MAR model.
      real(dp), intent(in) :: history(:) !! Conditioning history ending with the latest value.
      integer, intent(in) :: horizons !! Positive forecast horizon count.
      integer, intent(in) :: simulations !! Positive number of simulated paths.
      real(dp), intent(in), optional :: level !! Central empirical interval probability.
      type(mixar_forecast_sample_t) :: out
      type(mixar_simulation_t) :: simulation
      real(dp), allocatable :: ordered(:)
      real(dp) :: interval_level, tail_probability
      integer :: sample, horizon

      interval_level = 0.95_dp
      if (present(level)) interval_level = level
      if (.not. valid_model(model) .or. size(history) < maxval(model%order) .or. &
         horizons < 1 .or. simulations < 1 .or. interval_level <= 0.0_dp .or. &
         interval_level >= 1.0_dp) then
         out%info = 1
         return
      end if
      allocate(out%paths(horizons, simulations), out%mean(horizons))
      allocate(out%variance(horizons), out%lower(horizons), out%upper(horizons))
      do sample = 1, simulations
         simulation = mixar_simulate(model, horizons, history, 0)
         if (simulation%info /= 0) then
            out%info = 2
            return
         end if
         out%paths(:, sample) = simulation%series
      end do
      out%level = interval_level
      tail_probability = 0.5_dp*(1.0_dp - interval_level)
      do horizon = 1, horizons
         ordered = sorted(out%paths(horizon, :))
         out%mean(horizon) = sum(ordered)/real(simulations, dp)
         out%variance(horizon) = variance(ordered)
         out%lower(horizon) = quantile(ordered, tail_probability)
         out%upper(horizon) = quantile(ordered, 1.0_dp - tail_probability)
      end do
   end function mixar_forecast

   pure subroutine accumulate_normal_equations(series, time, weight, order, &
      include_shift, cross, rhs)
      !! Add one weighted AR regression row to normal equations.
      real(dp), intent(in) :: series(:) !! Observed time series.
      integer, intent(in) :: time !! Response observation index.
      real(dp), intent(in) :: weight !! Posterior component responsibility.
      integer, intent(in) :: order !! Component autoregressive order.
      logical, intent(in) :: include_shift !! Include an intercept column.
      real(dp), intent(inout) :: cross(:, :) !! Weighted regressor cross-product.
      real(dp), intent(inout) :: rhs(:) !! Weighted regressor-response product.
      real(dp) :: regressor(size(rhs))
      integer :: offset, lag

      offset = merge(1, 0, include_shift)
      if (include_shift) regressor(1) = 1.0_dp
      do lag = 1, order
         regressor(offset + lag) = series(time - lag)
      end do
      cross = cross + weight*spread(regressor, 2, size(rhs))* &
         spread(regressor, 1, size(rhs))
      rhs = rhs + weight*regressor*series(time)
   end subroutine accumulate_normal_equations

   pure function mixar_var_model(probability, order, intercept, ar, &
      innovation_covariance) result(model)
      !! Construct and validate a Gaussian mixture vector autoregression.
      real(dp), intent(in) :: probability(:) !! Component mixing probabilities.
      integer, intent(in) :: order(:) !! Autoregressive order of each component.
      real(dp), intent(in) :: intercept(:, :) !! Intercepts by variable and component.
      real(dp), intent(in) :: ar(:, :, :, :) !! Padded AR matrices by response, predictor, lag, and component.
      real(dp), intent(in) :: innovation_covariance(:, :, :) !! Innovation covariance by component.
      type(mixar_var_model_t) :: model

      model%probability = probability
      model%order = order
      model%intercept = intercept
      model%ar = ar
      model%innovation_covariance = innovation_covariance
      if (.not. valid_var_model(model)) model%info = 1
   end function mixar_var_model

   pure function mixar_var_component_locations(model, series) result(location)
      !! Compute component conditional means for every usable observation.
      type(mixar_var_model_t), intent(in) :: model !! Mixture VAR parameters.
      real(dp), intent(in) :: series(:, :) !! Observations by time and variable.
      real(dp), allocatable :: location(:, :, :)
      integer :: component, lag, maximum_order, observation, row

      maximum_order = maxval(model%order)
      allocate(location(max(0, size(series, 1) - maximum_order), &
         size(series, 2), size(model%probability)), source=0.0_dp)
      if (.not. valid_var_model(model) .or. size(series, 2) /= &
         size(model%intercept, 1)) return
      do observation = maximum_order + 1, size(series, 1)
         row = observation - maximum_order
         do component = 1, size(model%probability)
            location(row, :, component) = model%intercept(:, component)
            do lag = 1, model%order(component)
               location(row, :, component) = location(row, :, component) + &
                  matmul(model%ar(:, :, lag, component), &
                  series(observation - lag, :))
            end do
         end do
      end do
   end function mixar_var_component_locations

   pure function mixar_var_filter(model, series) result(filtered)
      !! Evaluate the Gaussian mixture VAR likelihood and posterior weights.
      type(mixar_var_model_t), intent(in) :: model !! Mixture VAR parameters.
      real(dp), intent(in) :: series(:, :) !! Observations by time and variable.
      type(mixar_var_filter_t) :: filtered
      real(dp), allocatable :: log_weight(:), difference(:)
      real(dp) :: maximum_log, row_log
      integer :: component, maximum_order, observation, row, variables

      maximum_order = maxval(model%order)
      variables = size(series, 2)
      filtered%start = maximum_order + 1
      if (.not. valid_var_model(model) .or. variables /= &
         size(model%intercept, 1) .or. size(series, 1) <= maximum_order) then
         filtered%info = 1
         return
      end if
      filtered%location = mixar_var_component_locations(model, series)
      allocate(filtered%responsibility(size(filtered%location, 1), &
         size(model%probability)))
      allocate(filtered%fitted(size(filtered%location, 1), variables))
      allocate(filtered%residual(size(filtered%location, 1), variables))
      allocate(filtered%conditional_covariance(variables, variables, &
         size(filtered%location, 1)))
      allocate(log_weight(size(model%probability)), difference(variables))
      filtered%log_likelihood = 0.0_dp
      do observation = maximum_order + 1, size(series, 1)
         row = observation - maximum_order
         do component = 1, size(model%probability)
            log_weight(component) = log(model%probability(component)) + &
               multivariate_normal_log_density(series(observation, :), &
               filtered%location(row, :, component), &
               model%innovation_covariance(:, :, component))
         end do
         maximum_log = maxval(log_weight)
         if (.not. ieee_is_finite(maximum_log)) then
            filtered%info = 2
            return
         end if
         row_log = maximum_log + log(sum(exp(log_weight - maximum_log)))
         filtered%responsibility(row, :) = exp(log_weight - row_log)
         filtered%log_likelihood = filtered%log_likelihood + row_log
         filtered%fitted(row, :) = matmul(filtered%location(row, :, :), &
            filtered%responsibility(row, :))
         filtered%residual(row, :) = series(observation, :) - &
            filtered%fitted(row, :)
         filtered%conditional_covariance(:, :, row) = 0.0_dp
         do component = 1, size(model%probability)
            difference = filtered%location(row, :, component) - &
               filtered%fitted(row, :)
            filtered%conditional_covariance(:, :, row) = &
               filtered%conditional_covariance(:, :, row) + &
               filtered%responsibility(row, component)*( &
               model%innovation_covariance(:, :, component) + &
               spread(difference, 2, variables)*spread(difference, 1, variables))
         end do
      end do
   end function mixar_var_filter

   pure real(dp) function mixar_var_log_likelihood(model, series) result(value)
      !! Return the conditional Gaussian mixture VAR log likelihood.
      type(mixar_var_model_t), intent(in) :: model !! Mixture VAR parameters.
      real(dp), intent(in) :: series(:, :) !! Observations by time and variable.
      type(mixar_var_filter_t) :: filtered

      filtered = mixar_var_filter(model, series)
      value = filtered%log_likelihood
   end function mixar_var_log_likelihood

   pure function mixar_var_fit(series, initial_model, max_iterations, &
      tolerance, covariance_floor) result(fit)
      !! Fit a ragged Gaussian mixture VAR by expectation maximization.
      real(dp), intent(in) :: series(:, :) !! Observations by time and variable.
      type(mixar_var_model_t), intent(in) :: initial_model !! Starting parameters and fixed component orders.
      integer, intent(in), optional :: max_iterations !! Maximum EM iterations.
      real(dp), intent(in), optional :: tolerance !! Relative log-likelihood convergence tolerance.
      real(dp), intent(in), optional :: covariance_floor !! Diagonal covariance regularization floor.
      type(mixar_var_fit_t) :: fit
      type(mixar_var_filter_t) :: current_filter, candidate_filter
      type(mixar_var_model_t) :: candidate
      real(dp), allocatable :: design(:, :), cross(:, :), gram(:, :)
      real(dp), allocatable :: gram_inverse(:, :), coefficient(:, :), difference(:)
      real(dp) :: effective, floor_value, previous, requested_tolerance
      integer :: component, iteration, lag, maximum, observations, parameter_count
      integer :: predictor, row, status, variables, requested_iterations

      requested_iterations = 200
      if (present(max_iterations)) requested_iterations = max_iterations
      requested_tolerance = 1.0e-8_dp
      if (present(tolerance)) requested_tolerance = tolerance
      floor_value = 1.0e-8_dp
      if (present(covariance_floor)) floor_value = covariance_floor
      if (.not. valid_var_model(initial_model) .or. requested_iterations < 1 .or. &
         requested_tolerance <= 0.0_dp .or. floor_value <= 0.0_dp) then
         fit%info = 1
         return
      end if
      variables = size(series, 2)
      maximum = maxval(initial_model%order)
      observations = size(series, 1) - maximum
      if (variables /= size(initial_model%intercept, 1) .or. observations < 2) then
         fit%info = 1
         return
      end if
      fit%model = initial_model
      current_filter = mixar_var_filter(fit%model, series)
      if (current_filter%info /= 0) then
         fit%info = 2
         return
      end if
      previous = current_filter%log_likelihood
      do iteration = 1, requested_iterations
         candidate = fit%model
         candidate%probability = sum(current_filter%responsibility, dim=1)/ &
            real(observations, dp)
         do component = 1, size(candidate%probability)
            allocate(design(observations, 1 + variables*candidate%order(component)))
            design(:, 1) = 1.0_dp
            do row = 1, observations
               do lag = 1, candidate%order(component)
                  design(row, 2 + (lag - 1)*variables:1 + lag*variables) = &
                     series(maximum + row - lag, :)
               end do
            end do
            allocate(gram(size(design, 2), size(design, 2)))
            allocate(cross(size(design, 2), variables))
            gram = matmul(transpose(design), design*spread( &
               current_filter%responsibility(:, component), 2, size(design, 2)))
            cross = matmul(transpose(design), &
               series(maximum + 1:, :)*spread( &
               current_filter%responsibility(:, component), 2, variables))
            allocate(gram_inverse(size(gram, 1), size(gram, 2)))
            call symmetric_pseudoinverse(gram, gram_inverse, status)
            if (status /= 0) then
               fit%info = 3
               return
            end if
            coefficient = matmul(gram_inverse, cross)
            candidate%intercept(:, component) = coefficient(1, :)
            candidate%ar(:, :, :, component) = 0.0_dp
            do lag = 1, candidate%order(component)
               do predictor = 1, variables
                  candidate%ar(:, predictor, lag, component) = &
                     coefficient(1 + (lag - 1)*variables + predictor, :)
               end do
            end do
            effective = sum(current_filter%responsibility(:, component))
            candidate%innovation_covariance(:, :, component) = 0.0_dp
            allocate(difference(variables))
            do row = 1, observations
               difference = series(maximum + row, :) - &
                  matmul(design(row, :), coefficient)
               candidate%innovation_covariance(:, :, component) = &
                  candidate%innovation_covariance(:, :, component) + &
                  current_filter%responsibility(row, component)* &
                  spread(difference, 2, variables)*spread(difference, 1, variables)
            end do
            candidate%innovation_covariance(:, :, component) = &
               candidate%innovation_covariance(:, :, component)/max(effective, tiny(1.0_dp))
            do predictor = 1, variables
               candidate%innovation_covariance(predictor, predictor, component) = &
                  candidate%innovation_covariance(predictor, predictor, component) + floor_value
            end do
            deallocate(design, gram, cross, gram_inverse, coefficient, difference)
         end do
         candidate_filter = mixar_var_filter(candidate, series)
         if (candidate_filter%info /= 0 .or. candidate_filter%log_likelihood < &
            previous - 100.0_dp*epsilon(1.0_dp)*max(1.0_dp, abs(previous))) then
            fit%info = 4
            return
         end if
         fit%model = candidate
         current_filter = candidate_filter
         fit%iterations = iteration
         if (abs(current_filter%log_likelihood - previous) <= &
            requested_tolerance*(1.0_dp + abs(previous))) then
            fit%converged = .true.
            exit
         end if
         previous = current_filter%log_likelihood
      end do
      fit%filter = current_filter
      fit%log_likelihood = current_filter%log_likelihood
      parameter_count = size(fit%model%probability) - 1
      do component = 1, size(fit%model%probability)
         parameter_count = parameter_count + variables + &
            variables*variables*fit%model%order(component) + &
            variables*(variables + 1)/2
      end do
      fit%aic = -2.0_dp*fit%log_likelihood + 2.0_dp*real(parameter_count, dp)
      fit%bic = -2.0_dp*fit%log_likelihood + &
         log(real(observations, dp))*real(parameter_count, dp)
   end function mixar_var_fit

   pure logical function mixar_var_is_stable(model, tolerance) result(stable)
      !! Test mixture VAR second-order stability by the Kronecker companion matrix.
      type(mixar_var_model_t), intent(in) :: model !! Mixture VAR parameters.
      real(dp), intent(in), optional :: tolerance !! Margin below unit spectral radius.
      complex(dp), allocatable :: eigenvalues(:)
      real(dp), allocatable :: companion(:, :), transition(:, :)
      real(dp) :: margin
      integer :: block, component, maximum, status, variables

      stable = .false.
      if (.not. valid_var_model(model)) return
      margin = 1.0e-10_dp
      if (present(tolerance)) margin = tolerance
      maximum = maxval(model%order)
      variables = size(model%intercept, 1)
      if (maximum == 0) then
         stable = .true.
         return
      end if
      allocate(companion(variables*maximum, variables*maximum), source=0.0_dp)
      allocate(transition((variables*maximum)**2, (variables*maximum)**2), &
         source=0.0_dp)
      do component = 1, size(model%probability)
         companion = 0.0_dp
         do block = 1, model%order(component)
            companion(1:variables, (block - 1)*variables + 1:block*variables) = &
               model%ar(:, :, block, component)
         end do
         do block = 2, maximum
            companion((block - 1)*variables + 1:block*variables, &
               (block - 2)*variables + 1:(block - 1)*variables) = &
               identity_block(variables)
         end do
         transition = transition + model%probability(component)* &
            kronecker_product(companion, companion)
      end do
      call general_eigenvalues(transition, eigenvalues, status)
      if (status == 0) stable = maxval(abs(eigenvalues)) < 1.0_dp - margin
   end function mixar_var_is_stable

   pure function mixar_var_simulate_from_draws(model, observations, initial, &
      uniforms, standard_draws, burnin) result(simulation)
      !! Simulate a mixture VAR from supplied uniforms and standard-normal draws.
      type(mixar_var_model_t), intent(in) :: model !! Mixture VAR parameters.
      integer, intent(in) :: observations !! Number of retained observations.
      real(dp), intent(in) :: initial(:, :) !! Initial history by time and variable.
      real(dp), intent(in) :: uniforms(:) !! Component-selection uniforms.
      real(dp), intent(in) :: standard_draws(:, :) !! Independent standard normals by time and variable.
      integer, intent(in), optional :: burnin !! Number of discarded observations.
      type(mixar_var_simulation_t) :: simulation
      real(dp), allocatable :: all_series(:, :), innovation(:), mean(:)
      real(dp) :: cumulative
      integer :: component, discard, i, lag, maximum, status, total, variables

      discard = 0
      if (present(burnin)) discard = burnin
      maximum = maxval(model%order)
      variables = size(model%intercept, 1)
      total = observations + discard
      if (.not. valid_var_model(model) .or. observations < 1 .or. discard < 0 .or. &
         size(initial, 1) < maximum .or. size(initial, 2) /= variables .or. &
         size(uniforms) < total .or. size(standard_draws, 1) < total .or. &
         size(standard_draws, 2) /= variables) then
         simulation%info = 1
         return
      end if
      allocate(all_series(maximum + total, variables), source=0.0_dp)
      if (maximum > 0) all_series(:maximum, :) = initial( &
         size(initial, 1) - maximum + 1:, :)
      allocate(simulation%regime(observations), simulation%series(observations, variables))
      allocate(innovation(variables), mean(variables))
      do i = 1, total
         component = size(model%probability)
         cumulative = 0.0_dp
         do lag = 1, size(model%probability)
            cumulative = cumulative + model%probability(lag)
            if (uniforms(i) <= cumulative) then
               component = lag
               exit
            end if
         end do
         mean = model%intercept(:, component)
         do lag = 1, model%order(component)
            mean = mean + matmul(model%ar(:, :, lag, component), &
               all_series(maximum + i - lag, :))
         end do
         call multivariate_normal_from_standard(mean, &
            model%innovation_covariance(:, :, component), &
            standard_draws(i, :), innovation, status)
         if (status /= 0) then
            simulation%info = 2
            return
         end if
         all_series(maximum + i, :) = innovation
         if (i > discard) simulation%regime(i - discard) = component
      end do
      simulation%series = all_series(maximum + discard + 1:, :)
      simulation%model = model
      simulation%burnin = discard
   end function mixar_var_simulate_from_draws

   function mixar_var_simulate(model, observations, initial, burnin) result(simulation)
      !! Simulate a mixture VAR using the shared random-number stream.
      type(mixar_var_model_t), intent(in) :: model !! Mixture VAR parameters.
      integer, intent(in) :: observations !! Number of retained observations.
      real(dp), intent(in) :: initial(:, :) !! Initial history by time and variable.
      integer, intent(in), optional :: burnin !! Number of discarded observations.
      type(mixar_var_simulation_t) :: simulation
      real(dp), allocatable :: uniforms(:), standard_draws(:, :)
      integer :: discard, i

      discard = 0
      if (present(burnin)) discard = burnin
      allocate(uniforms(max(0, observations + discard)))
      do i = 1, size(uniforms)
         uniforms(i) = random_uniform()
      end do
      allocate(standard_draws(max(0, observations + discard), &
         size(model%intercept, 1)))
      call random_standard_normal_matrix(standard_draws)
      simulation = mixar_var_simulate_from_draws(model, observations, initial, &
         uniforms, standard_draws, discard)
   end function mixar_var_simulate

   pure function mixar_var_forecast_from_draws(model, history, uniforms, &
      standard_draws, level) result(forecast)
      !! Generate mixture VAR forecast paths from supplied random draws.
      type(mixar_var_model_t), intent(in) :: model !! Mixture VAR parameters.
      real(dp), intent(in) :: history(:, :) !! Observed history by time and variable.
      real(dp), intent(in) :: uniforms(:, :) !! Selection uniforms by horizon and simulation.
      real(dp), intent(in) :: standard_draws(:, :, :) !! Standard normals by horizon, variable, and simulation.
      real(dp), intent(in), optional :: level !! Central prediction interval coverage.
      type(mixar_var_forecast_t) :: forecast
      type(mixar_var_simulation_t) :: simulation
      real(dp), allocatable :: ordered(:), difference(:)
      integer :: horizon, i, j, simulations, variables
      real(dp) :: alpha

      horizon = size(uniforms, 1)
      simulations = size(uniforms, 2)
      variables = size(history, 2)
      forecast%level = 0.95_dp
      if (present(level)) forecast%level = level
      if (.not. valid_var_model(model) .or. horizon < 1 .or. simulations < 2 .or. &
         size(history, 1) < maxval(model%order) .or. &
         any(shape(standard_draws) /= [horizon, variables, simulations]) .or. &
         forecast%level <= 0.0_dp .or. forecast%level >= 1.0_dp) then
         forecast%info = 1
         return
      end if
      allocate(forecast%paths(horizon, variables, simulations))
      do i = 1, simulations
         simulation = mixar_var_simulate_from_draws(model, horizon, history, &
            uniforms(:, i), standard_draws(:, :, i))
         if (simulation%info /= 0) then
            forecast%info = 2
            return
         end if
         forecast%paths(:, :, i) = simulation%series
      end do
      allocate(forecast%mean(horizon, variables), forecast%lower(horizon, variables))
      allocate(forecast%upper(horizon, variables))
      allocate(forecast%covariance(variables, variables, horizon), source=0.0_dp)
      alpha = 0.5_dp*(1.0_dp - forecast%level)
      allocate(difference(variables))
      do i = 1, horizon
         forecast%mean(i, :) = sum(forecast%paths(i, :, :), dim=2)/ &
            real(simulations, dp)
         do j = 1, simulations
            difference = forecast%paths(i, :, j) - forecast%mean(i, :)
            forecast%covariance(:, :, i) = forecast%covariance(:, :, i) + &
               spread(difference, 2, variables)*spread(difference, 1, variables)
         end do
         forecast%covariance(:, :, i) = forecast%covariance(:, :, i)/ &
            real(simulations - 1, dp)
         do j = 1, variables
            ordered = sorted(forecast%paths(i, j, :))
            forecast%lower(i, j) = quantile(ordered, alpha)
            forecast%upper(i, j) = quantile(ordered, 1.0_dp - alpha)
         end do
      end do
   end function mixar_var_forecast_from_draws

   function mixar_var_forecast(model, history, horizon, simulations, level) &
      result(forecast)
      !! Simulate mixture VAR point, interval, covariance, and path forecasts.
      type(mixar_var_model_t), intent(in) :: model !! Mixture VAR parameters.
      real(dp), intent(in) :: history(:, :) !! Observed history by time and variable.
      integer, intent(in) :: horizon !! Forecast horizon.
      integer, intent(in) :: simulations !! Number of simulated paths.
      real(dp), intent(in), optional :: level !! Central prediction interval coverage.
      type(mixar_var_forecast_t) :: forecast
      real(dp), allocatable :: uniforms(:, :), standard_draws(:, :, :)
      integer :: i, j

      if (horizon < 1 .or. simulations < 2) then
         forecast%info = 1
         return
      end if
      allocate(uniforms(horizon, simulations))
      allocate(standard_draws(horizon, size(history, 2), simulations))
      do j = 1, simulations
         do i = 1, horizon
            uniforms(i, j) = random_uniform()
         end do
         call random_standard_normal_matrix(standard_draws(:, :, j))
      end do
      if (present(level)) then
         forecast = mixar_var_forecast_from_draws(model, history, uniforms, &
            standard_draws, level)
      else
         forecast = mixar_var_forecast_from_draws(model, history, uniforms, &
            standard_draws)
      end if
   end function mixar_var_forecast

   pure function mixar_var_diagnose(model, series, lag) result(diagnostics)
      !! Compute standardized residual, whiteness, and classification diagnostics.
      type(mixar_var_model_t), intent(in) :: model !! Mixture VAR parameters.
      real(dp), intent(in) :: series(:, :) !! Observations by time and variable.
      integer, intent(in) :: lag !! Maximum lag for white-noise testing.
      type(mixar_var_diagnostics_t) :: diagnostics
      type(mixar_var_filter_t) :: filtered
      real(dp), allocatable :: eigenvalues(:), eigenvectors(:, :), residual(:)
      integer :: component, i, status, variables

      filtered = mixar_var_filter(model, series)
      if (filtered%info /= 0 .or. lag < 1 .or. lag >= size(filtered%residual, 1)) then
         diagnostics%info = 1
         return
      end if
      variables = size(series, 2)
      diagnostics%residual = filtered%residual
      diagnostics%posterior_probability = filtered%responsibility
      allocate(diagnostics%standardized_residual(size(filtered%residual, 1), variables))
      allocate(diagnostics%classification(size(filtered%residual, 1)))
      allocate(diagnostics%maximum_posterior(size(filtered%residual, 1)))
      allocate(residual(variables))
      do i = 1, size(filtered%residual, 1)
         component = maxloc(filtered%responsibility(i, :), dim=1)
         diagnostics%classification(i) = component
         diagnostics%maximum_posterior(i) = filtered%responsibility(i, component)
         call symmetric_eigen(filtered%conditional_covariance(:, :, i), &
            eigenvalues, eigenvectors, status)
         if (status /= 0 .or. minval(eigenvalues) <= tiny(1.0_dp)) then
            diagnostics%info = 2
            return
         end if
         residual = matmul(transpose(eigenvectors), filtered%residual(i, :))
         diagnostics%standardized_residual(i, :) = matmul(eigenvectors, &
            residual/sqrt(eigenvalues))
         do component = 1, size(model%probability)
            if (filtered%responsibility(i, component) > 0.0_dp) &
               diagnostics%mean_entropy = diagnostics%mean_entropy - &
               filtered%responsibility(i, component)* &
               log(filtered%responsibility(i, component))
         end do
      end do
      diagnostics%mean_entropy = diagnostics%mean_entropy/ &
         real(size(filtered%residual, 1), dp)
      diagnostics%white_noise_test = multivariate_white_noise_test( &
         diagnostics%standardized_residual, lag)
      if (diagnostics%white_noise_test%info /= 0) diagnostics%info = 3
   end function mixar_var_diagnose

   pure logical function valid_var_model(model) result(valid)
      !! Check dimensions, probabilities, orders, and covariance definiteness.
      type(mixar_var_model_t), intent(in) :: model !! Mixture VAR parameters.
      real(dp), allocatable :: inverse(:, :)
      real(dp) :: log_determinant
      integer :: component, components, status, variables

      valid = .false.
      if (.not. allocated(model%probability) .or. .not. allocated(model%order) .or. &
         .not. allocated(model%intercept) .or. .not. allocated(model%ar) .or. &
         .not. allocated(model%innovation_covariance)) return
      components = size(model%probability)
      variables = size(model%intercept, 1)
      if (components < 1 .or. variables < 1 .or. size(model%order) /= components .or. &
         size(model%intercept, 2) /= components .or. &
         any(shape(model%ar) /= [variables, variables, size(model%ar, 3), components]) .or. &
         any(shape(model%innovation_covariance) /= [variables, variables, components]) .or. &
         any(model%order < 0) .or. any(model%order > size(model%ar, 3)) .or. &
         any(model%probability <= 0.0_dp) .or. &
         abs(sum(model%probability) - 1.0_dp) > 1.0e-10_dp .or. &
         .not. all(ieee_is_finite(model%intercept)) .or. &
         .not. all(ieee_is_finite(model%ar))) return
      allocate(inverse(variables, variables))
      do component = 1, components
         call inverse_logdet(model%innovation_covariance(:, :, component), &
            inverse, log_determinant, status, 1.0e-12_dp)
         if (status /= 0 .or. .not. ieee_is_finite(log_determinant)) return
      end do
      valid = .true.
   end function valid_var_model

   pure function identity_block(size_) result(identity)
      !! Construct an identity matrix for a companion-matrix block.
      integer, intent(in) :: size_ !! Matrix dimension.
      real(dp) :: identity(size_, size_)
      integer :: i

      identity = 0.0_dp
      do i = 1, size_
         identity(i, i) = 1.0_dp
      end do
   end function identity_block

   function mixar_bayesian_random(observations, components, maximum_order, &
      iterations) result(random_values)
      !! Generate independent variates for reproducible Bayesian MAR sampling.
      integer, intent(in) :: observations !! Number of conditional observations.
      integer, intent(in) :: components !! Number of mixture components.
      integer, intent(in) :: maximum_order !! Largest component AR order.
      integer, intent(in) :: iterations !! Total burn-in and retained iterations.
      type(mixar_bayesian_random_t) :: random_values
      integer :: component, i, lag

      if (observations < 1 .or. components < 1 .or. maximum_order < 0 .or. &
         iterations < 1) return
      allocate(random_values%allocation_uniform(observations, iterations))
      allocate(random_values%probability_uniform(components, iterations))
      allocate(random_values%mean_normal(components, iterations))
      allocate(random_values%lambda_uniform(iterations))
      allocate(random_values%precision_uniform(components, iterations))
      allocate(random_values%ar_normal(maximum_order, components, iterations))
      allocate(random_values%acceptance_uniform(components, iterations))
      do i = 1, iterations
         random_values%lambda_uniform(i) = random_uniform()
         do component = 1, components
            random_values%probability_uniform(component, i) = random_uniform()
            random_values%mean_normal(component, i) = random_standard_normal()
            random_values%precision_uniform(component, i) = random_uniform()
            random_values%acceptance_uniform(component, i) = random_uniform()
            do lag = 1, maximum_order
               random_values%ar_normal(lag, component, i) = &
                  random_standard_normal()
            end do
         end do
         do component = 1, observations
            random_values%allocation_uniform(component, i) = random_uniform()
         end do
      end do
   end function mixar_bayesian_random

   pure function mixar_bayesian_sample_from_random(series, initial_model, &
      retained_draws, burnin, random_values, fix_shift, precision_shape, &
      precision_hyper_shape, proposal_scale, dirichlet_prior) result(draws)
      !! Sample the Bayesian Gaussian MAR posterior from supplied independent variates.
      real(dp), intent(in) :: series(:) !! Observed time series.
      type(mixar_model_t), intent(in) :: initial_model !! Initial Gaussian MAR parameters.
      integer, intent(in) :: retained_draws !! Number of posterior draws retained after burn-in.
      integer, intent(in) :: burnin !! Number of discarded initial draws.
      type(mixar_bayesian_random_t), intent(in) :: random_values !! Supplied independent random variates.
      logical, intent(in), optional :: fix_shift !! Hold component intercepts fixed when true.
      real(dp), intent(in), optional :: precision_shape !! Component-precision Gamma prior shape.
      real(dp), intent(in), optional :: precision_hyper_shape !! Common precision-rate Gamma prior shape.
      real(dp), intent(in), optional :: proposal_scale(:) !! AR random-walk standard deviations by component.
      real(dp), intent(in), optional :: dirichlet_prior(:) !! Dirichlet concentration parameters.
      type(mixar_bayesian_draws_t) :: draws
      type(mixar_model_t) :: current, proposal
      real(dp), allocatable :: allocation_probability(:), alpha(:), candidate_ar(:)
      real(dp), allocatable :: component_mean(:), precision(:), tau(:)
      real(dp) :: a, b, b_coefficient, c, cumulative, effective_precision
      real(dp) :: lambda, log_ratio, maximum_log, posterior_mean, proposal_ssr
      real(dp) :: range_value, reference_mean, scale_value, sum_error, current_ssr
      integer, allocatable :: allocation(:), count_value(:)
      integer :: component, i, iteration, lag, maximum_order, observation
      integer :: retained, total, usable
      logical :: hold_shift, proposal_allowed

      hold_shift = .false.
      if (present(fix_shift)) hold_shift = fix_shift
      a = 0.2_dp
      if (present(precision_hyper_shape)) a = precision_hyper_shape
      c = 2.0_dp
      if (present(precision_shape)) c = precision_shape
      maximum_order = maxval(initial_model%order)
      usable = size(series) - maximum_order
      total = retained_draws + burnin
      if (.not. valid_model(initial_model) .or. &
         any(initial_model%degrees_of_freedom /= 0.0_dp) .or. &
         retained_draws < 1 .or. burnin < 0 .or. usable < 2 .or. &
         a <= 0.0_dp .or. c <= 0.0_dp .or. &
         .not. valid_bayesian_random(random_values, usable, &
         size(initial_model%probability), maximum_order, total)) then
         draws%info = 1
         return
      end if
      range_value = maxval(series) - minval(series)
      if (range_value <= sqrt(epsilon(1.0_dp))) then
         draws%info = 1
         return
      end if
      allocate(alpha(size(initial_model%probability)), source=1.0_dp)
      if (present(dirichlet_prior)) then
         if (size(dirichlet_prior) /= size(alpha) .or. &
            any(dirichlet_prior <= 0.0_dp)) then
            draws%info = 1
            return
         end if
         alpha = dirichlet_prior
      end if
      allocate(tau(size(initial_model%probability)), source=0.05_dp)
      if (present(proposal_scale)) then
         if (size(proposal_scale) /= size(tau) .or. any(proposal_scale <= 0.0_dp)) then
            draws%info = 1
            return
         end if
         tau = proposal_scale
      end if
      current = initial_model
      allocate(component_mean(size(alpha)), precision(size(alpha)))
      do component = 1, size(alpha)
         b_coefficient = 1.0_dp - sum(current%ar( &
            :current%order(component), component))
         if (abs(b_coefficient) <= sqrt(epsilon(1.0_dp))) then
            draws%info = 2
            return
         end if
         component_mean(component) = current%shift(component)/b_coefficient
      end do
      precision = 1.0_dp/current%scale**2
      reference_mean = 0.5_dp*(minval(series) + maxval(series))
      effective_precision = 1.0_dp/range_value
      b = 100.0_dp*a/(c*range_value**2)
      allocate(allocation(usable), count_value(size(alpha)), &
         allocation_probability(size(alpha)), candidate_ar(maximum_order))
      allocate(draws%probability(retained_draws, size(alpha)))
      allocate(draws%scale(retained_draws, size(alpha)))
      allocate(draws%precision(retained_draws, size(alpha)))
      allocate(draws%shift(retained_draws, size(alpha)))
      allocate(draws%component_mean(retained_draws, size(alpha)))
      allocate(draws%ar(retained_draws, maximum_order, size(alpha)), source=0.0_dp)
      allocate(draws%lambda(retained_draws))
      allocate(draws%acceptance_rate(size(alpha)), source=0.0_dp)
      draws%order = current%order
      draws%burnin = burnin
      draws%fixed_shift = hold_shift
      do iteration = 1, total
         count_value = 0
         do i = 1, usable
            observation = maximum_order + i
            do component = 1, size(alpha)
               posterior_mean = current%shift(component)
               do lag = 1, current%order(component)
                  posterior_mean = posterior_mean + current%ar(lag, component)* &
                     series(observation - lag)
               end do
               allocation_probability(component) = log(current%probability(component)) - &
                  log(current%scale(component)) - 0.5_dp*log_two_pi - &
                  0.5_dp*((series(observation) - posterior_mean)/ &
                  current%scale(component))**2
            end do
            maximum_log = maxval(allocation_probability)
            allocation_probability = exp(allocation_probability - maximum_log)
            allocation_probability = allocation_probability/sum(allocation_probability)
            cumulative = 0.0_dp
            allocation(i) = size(alpha)
            do component = 1, size(alpha)
               cumulative = cumulative + allocation_probability(component)
               if (random_values%allocation_uniform(i, iteration) <= cumulative) then
                  allocation(i) = component
                  exit
               end if
            end do
            count_value(allocation(i)) = count_value(allocation(i)) + 1
         end do
         do component = 1, size(alpha)
            current%probability(component) = gamma_quantile( &
               random_values%probability_uniform(component, iteration), &
               alpha(component) + real(count_value(component), dp))
         end do
         current%probability = current%probability/sum(current%probability)
         if (.not. hold_shift) then
            do component = 1, size(alpha)
               b_coefficient = 1.0_dp - sum(current%ar( &
                  :current%order(component), component))
               sum_error = bayesian_unshifted_error_sum(series, allocation, &
                  component, maximum_order, current%order(component), &
                  current%ar(:, component))
               scale_value = precision(component)*real(count_value(component), dp)* &
                  b_coefficient**2 + effective_precision
               posterior_mean = (precision(component)*b_coefficient*sum_error + &
                  effective_precision*reference_mean)/scale_value
               component_mean(component) = posterior_mean + &
                  random_values%mean_normal(component, iteration)/sqrt(scale_value)
               current%shift(component) = component_mean(component)*b_coefficient
            end do
         end if
         lambda = gamma_quantile(random_values%lambda_uniform(iteration), &
            a + real(size(alpha), dp)*c)/(b + sum(precision))
         do component = 1, size(alpha)
            current_ssr = bayesian_component_ssr(series, allocation, component, &
               maximum_order, current%shift(component), current%order(component), &
               current%ar(:, component))
            precision(component) = gamma_quantile( &
               random_values%precision_uniform(component, iteration), &
               c + 0.5_dp*real(count_value(component), dp))/ &
               (lambda + 0.5_dp*current_ssr)
            current%scale(component) = 1.0_dp/sqrt(precision(component))
         end do
         do component = 1, size(alpha)
            candidate_ar = current%ar(:, component)
            candidate_ar(:current%order(component)) = &
               candidate_ar(:current%order(component)) + tau(component)* &
               random_values%ar_normal(:current%order(component), component, iteration)
            current_ssr = bayesian_component_ssr(series, allocation, component, &
               maximum_order, current%shift(component), current%order(component), &
               current%ar(:, component))
            proposal = current
            proposal%ar(:, component) = candidate_ar
            proposal_allowed = .true.
            if (hold_shift) then
               b_coefficient = 1.0_dp - sum(candidate_ar(:current%order(component)))
               proposal_allowed = abs(b_coefficient) > sqrt(epsilon(1.0_dp))
               if (proposal_allowed) &
                  component_mean(component) = current%shift(component)/b_coefficient
            else
               proposal%shift(component) = component_mean(component)*(1.0_dp - &
                  sum(candidate_ar(:current%order(component))))
            end if
            proposal_ssr = bayesian_component_ssr(series, allocation, component, &
               maximum_order, proposal%shift(component), current%order(component), &
               candidate_ar)
            log_ratio = -0.5_dp*precision(component)*(proposal_ssr - current_ssr)
            if (log(max(random_values%acceptance_uniform(component, iteration), &
               tiny(1.0_dp))) <= min(0.0_dp, log_ratio) .and. &
               proposal_allowed .and. mixar_is_stable(proposal)) then
               current = proposal
               if (iteration > burnin) &
                  draws%acceptance_rate(component) = &
                  draws%acceptance_rate(component) + 1.0_dp
            end if
         end do
         if (hold_shift) then
            do component = 1, size(alpha)
               b_coefficient = 1.0_dp - sum(current%ar( &
                  :current%order(component), component))
               component_mean(component) = current%shift(component)/b_coefficient
            end do
         end if
         if (iteration > burnin) then
            retained = iteration - burnin
            draws%probability(retained, :) = current%probability
            draws%scale(retained, :) = current%scale
            draws%precision(retained, :) = precision
            draws%shift(retained, :) = current%shift
            draws%component_mean(retained, :) = component_mean
            draws%ar(retained, :, :) = current%ar
            draws%lambda(retained) = lambda
         end if
      end do
      draws%last_allocation = allocation
      draws%acceptance_rate = draws%acceptance_rate/real(retained_draws, dp)
      draws%posterior_mean_model = bayesian_posterior_mean_model(draws)
   end function mixar_bayesian_sample_from_random

   function mixar_bayesian_sample(series, initial_model, retained_draws, &
      burnin, fix_shift, precision_shape, precision_hyper_shape, &
      proposal_scale, dirichlet_prior) result(draws)
      !! Sample the Bayesian Gaussian MAR posterior using the shared random stream.
      real(dp), intent(in) :: series(:) !! Observed time series.
      type(mixar_model_t), intent(in) :: initial_model !! Initial Gaussian MAR parameters.
      integer, intent(in) :: retained_draws !! Number of retained posterior draws.
      integer, intent(in), optional :: burnin !! Number of discarded initial draws.
      logical, intent(in), optional :: fix_shift !! Hold component intercepts fixed when true.
      real(dp), intent(in), optional :: precision_shape !! Component-precision Gamma prior shape.
      real(dp), intent(in), optional :: precision_hyper_shape !! Common precision-rate Gamma prior shape.
      real(dp), intent(in), optional :: proposal_scale(:) !! AR random-walk standard deviations by component.
      real(dp), intent(in), optional :: dirichlet_prior(:) !! Dirichlet concentration parameters.
      type(mixar_bayesian_draws_t) :: draws
      type(mixar_bayesian_random_t) :: random_values
      real(dp), allocatable :: alpha(:), tau(:)
      real(dp) :: a, c
      integer :: discarded, maximum_order, total
      logical :: held

      discarded = 0
      if (present(burnin)) discarded = burnin
      held = .false.
      if (present(fix_shift)) held = fix_shift
      a = 0.2_dp
      if (present(precision_hyper_shape)) a = precision_hyper_shape
      c = 2.0_dp
      if (present(precision_shape)) c = precision_shape
      allocate(tau(size(initial_model%probability)), source=0.05_dp)
      if (present(proposal_scale)) tau = proposal_scale
      allocate(alpha(size(initial_model%probability)), source=1.0_dp)
      if (present(dirichlet_prior)) alpha = dirichlet_prior
      total = retained_draws + discarded
      maximum_order = maxval(initial_model%order)
      random_values = mixar_bayesian_random(size(series) - maximum_order, &
         size(initial_model%probability), maximum_order, total)
      draws = mixar_bayesian_sample_from_random(series, initial_model, &
         retained_draws, discarded, random_values, held, c, a, tau, alpha)
   end function mixar_bayesian_sample

   pure function mixar_bayesian_relabel(input_draws) result(draws)
      !! Resolve equal-order component labels by sorting each draw on component means.
      type(mixar_bayesian_draws_t), intent(in) :: input_draws !! Bayesian MAR draws to relabel.
      type(mixar_bayesian_draws_t) :: draws
      integer, allocatable :: permutation(:)
      integer :: component, i, j, held

      draws = input_draws
      if (draws%info /= 0 .or. .not. allocated(draws%order) .or. &
         any(draws%order /= draws%order(1))) then
         draws%info = 1
         return
      end if
      allocate(permutation(size(draws%order)))
      do i = 1, size(draws%probability, 1)
         permutation = [(component, component=1, size(draws%order))]
         do component = 2, size(permutation)
            held = permutation(component)
            j = component - 1
            do while (j >= 1)
               if (draws%component_mean(i, permutation(j)) <= &
                  draws%component_mean(i, held)) exit
               permutation(j + 1) = permutation(j)
               j = j - 1
            end do
            permutation(j + 1) = held
         end do
         draws%probability(i, :) = input_draws%probability(i, permutation)
         draws%scale(i, :) = input_draws%scale(i, permutation)
         draws%precision(i, :) = input_draws%precision(i, permutation)
         draws%shift(i, :) = input_draws%shift(i, permutation)
         draws%component_mean(i, :) = input_draws%component_mean(i, permutation)
         draws%ar(i, :, :) = input_draws%ar(i, :, permutation)
      end do
      if (allocated(draws%last_allocation)) then
         do component = 1, size(draws%last_allocation)
            draws%last_allocation(component) = find_integer( &
               permutation, input_draws%last_allocation(component))
         end do
      end if
      draws%posterior_mean_model = bayesian_posterior_mean_model(draws)
   end function mixar_bayesian_relabel

   pure function mixar_marginal_likelihood(series, draws, precision_shape, &
      precision_hyper_shape, dirichlet_prior) result(out)
      !! Estimate a Chib-style marginal likelihood using posterior ordinate KDEs.
      real(dp), intent(in) :: series(:) !! Observed time series.
      type(mixar_bayesian_draws_t), intent(in) :: draws !! Retained Bayesian MAR draws.
      real(dp), intent(in), optional :: precision_shape !! Component-precision Gamma prior shape.
      real(dp), intent(in), optional :: precision_hyper_shape !! Common precision-rate Gamma prior shape.
      real(dp), intent(in), optional :: dirichlet_prior(:) !! Dirichlet concentration parameters.
      type(mixar_marginal_likelihood_t) :: out
      real(dp), allocatable :: alpha(:), mean_value(:), precision_value(:)
      real(dp) :: a, b, c, ka, lambda_value, range_value, zeta
      integer :: component, lag

      if (draws%info /= 0 .or. size(draws%probability, 1) < 10) then
         out%info = 1
         return
      end if
      a = 0.2_dp
      if (present(precision_hyper_shape)) a = precision_hyper_shape
      c = 2.0_dp
      if (present(precision_shape)) c = precision_shape
      allocate(alpha(size(draws%probability, 2)), source=1.0_dp)
      if (present(dirichlet_prior)) then
         if (size(dirichlet_prior) /= size(alpha) .or. any(dirichlet_prior <= 0.0_dp)) then
            out%info = 1
            return
         end if
         alpha = dirichlet_prior
      end if
      range_value = maxval(series) - minval(series)
      if (range_value <= sqrt(epsilon(1.0_dp))) then
         out%info = 1
         return
      end if
      zeta = 0.5_dp*(minval(series) + maxval(series))
      ka = 1.0_dp/range_value
      b = 100.0_dp*a/(c*range_value**2)
      out%representative_model = draws%posterior_mean_model
      out%log_likelihood = mixar_log_likelihood(out%representative_model, series)
      mean_value = sum(draws%component_mean, dim=1)/ &
         real(size(draws%component_mean, 1), dp)
      precision_value = sum(draws%precision, dim=1)/ &
         real(size(draws%precision, 1), dp)
      lambda_value = sum(draws%lambda)/real(size(draws%lambda), dp)
      out%log_prior_density = gamma_log_density(lambda_value, a, b) + &
         dirichlet_log_density(out%representative_model%probability, alpha)
      out%log_posterior_ordinate = kde_log_at(draws%lambda, lambda_value)
      do component = 1, size(alpha)
         out%log_prior_density = out%log_prior_density + &
            gamma_log_density(precision_value(component), c, lambda_value)
         out%log_posterior_ordinate = out%log_posterior_ordinate + &
            kde_log_at(draws%precision(:, component), precision_value(component))
         if (.not. draws%fixed_shift) then
            out%log_prior_density = out%log_prior_density + &
               0.5_dp*(log(ka) - log_two_pi) - &
               0.5_dp*ka*(mean_value(component) - zeta)**2
            out%log_posterior_ordinate = out%log_posterior_ordinate + &
               kde_log_at(draws%component_mean(:, component), mean_value(component))
         end if
         do lag = 1, draws%order(component)
            out%log_prior_density = out%log_prior_density - 0.5_dp*( &
               log_two_pi + out%representative_model%ar(lag, component)**2)
            out%log_posterior_ordinate = out%log_posterior_ordinate + &
               kde_log_at(draws%ar(:, lag, component), &
               out%representative_model%ar(lag, component))
         end do
      end do
      do component = 1, size(alpha) - 1
         out%log_posterior_ordinate = out%log_posterior_ordinate + &
            kde_log_at(draws%probability(:, component), &
            out%representative_model%probability(component))
      end do
      if (.not. ieee_is_finite(out%log_likelihood) .or. &
         .not. ieee_is_finite(out%log_prior_density) .or. &
         .not. ieee_is_finite(out%log_posterior_ordinate)) then
         out%info = 2
         return
      end if
      out%log_marginal_likelihood = out%log_likelihood + &
         out%log_prior_density - out%log_posterior_ordinate
   end function mixar_marginal_likelihood

   pure subroutine mixar_order_birth_death_probability(prior_method, &
      prior_parameter, order, birth_probability, death_probability, info)
      !! Compute mixAR birth and death proposal probabilities at one order.
      integer, intent(in) :: prior_method !! Flat, ratio, or Poisson proposal method.
      real(dp), intent(in) :: prior_parameter !! Positive ratio or Poisson parameter.
      integer, intent(in) :: order !! Current nonnegative AR order.
      real(dp), intent(out) :: birth_probability !! Probability of proposing a birth.
      real(dp), intent(out) :: death_probability !! Probability of proposing a death.
      integer, intent(out) :: info !! Status code; zero indicates success.

      info = 0
      if (order < 0 .or. (prior_method /= mixar_order_prior_flat .and. &
         prior_method /= mixar_order_prior_ratio .and. &
         prior_method /= mixar_order_prior_poisson) .or. &
         (prior_method /= mixar_order_prior_flat .and. prior_parameter <= 0.0_dp)) then
         birth_probability = 0.0_dp
         death_probability = 0.0_dp
         info = 1
         return
      end if
      select case (prior_method)
      case (mixar_order_prior_ratio)
         birth_probability = prior_parameter/(prior_parameter + real(order, dp))
      case (mixar_order_prior_poisson)
         birth_probability = prior_parameter/real(order + 1, dp)
      case default
         birth_probability = 0.5_dp
      end select
      birth_probability = max(0.0_dp, min(1.0_dp, birth_probability))
      death_probability = 1.0_dp - birth_probability
   end subroutine mixar_order_birth_death_probability

   function mixar_order_random(observations, components, maximum_order, &
      iterations) result(random_values)
      !! Generate independent variates for reproducible mixture-order RJMCMC.
      integer, intent(in) :: observations !! Total observed series length.
      integer, intent(in) :: components !! Number of mixture components.
      integer, intent(in) :: maximum_order !! Largest permitted component order.
      integer, intent(in) :: iterations !! Number of RJMCMC iterations.
      type(mixar_order_random_t) :: random_values
      integer :: i

      if (observations < 2 .or. components < 1 .or. maximum_order < 0 .or. &
         iterations < 1) return
      allocate(random_values%parameter_random(iterations))
      allocate(random_values%component_uniform(iterations))
      allocate(random_values%direction_uniform(iterations))
      allocate(random_values%coefficient_uniform(iterations))
      allocate(random_values%acceptance_uniform(iterations))
      do i = 1, iterations
         random_values%parameter_random(i) = mixar_bayesian_random( &
            observations, components, maximum_order, 2)
         random_values%component_uniform(i) = random_uniform()
         random_values%direction_uniform(i) = random_uniform()
         random_values%coefficient_uniform(i) = random_uniform()
         random_values%acceptance_uniform(i) = random_uniform()
      end do
   end function mixar_order_random

   pure function mixar_order_select_from_random(series, initial_model, &
      maximum_order, iterations, random_values, proposal_scale, prior_method, &
      prior_parameter, fix_shift) result(selection)
      !! Select component AR orders by supplied-random reversible-jump MCMC.
      real(dp), intent(in) :: series(:) !! Observed time series.
      type(mixar_model_t), intent(in) :: initial_model !! Initial Gaussian MAR model.
      integer, intent(in) :: maximum_order !! Largest permitted component order.
      integer, intent(in) :: iterations !! Number of RJMCMC iterations.
      type(mixar_order_random_t), intent(in) :: random_values !! Supplied independent variates.
      real(dp), intent(in) :: proposal_scale(:) !! AR random-walk scales by component.
      integer, intent(in), optional :: prior_method !! Flat, ratio, or Poisson proposal method.
      real(dp), intent(in), optional :: prior_parameter !! Positive ratio or Poisson parameter.
      logical, intent(in), optional :: fix_shift !! Hold component intercepts fixed when true.
      type(mixar_order_selection_t) :: selection
      type(mixar_bayesian_draws_t) :: parameter_draw
      type(mixar_model_t) :: current, proposal
      real(dp), allocatable :: component_mean(:), precision(:)
      real(dp) :: birth_current, birth_proposed, death_current, death_proposed
      real(dp) :: current_ssr, log_ratio, proposed_ssr, removed_density
      integer, allocatable :: common_allocation(:)
      integer :: component, current_maximum, current_order, direction_status
      integer :: i, method_value, proposed_order, start_index
      real(dp) :: parameter_value
      logical :: held

      method_value = mixar_order_prior_flat
      if (present(prior_method)) method_value = prior_method
      parameter_value = 1.0_dp
      if (present(prior_parameter)) parameter_value = prior_parameter
      held = .false.
      if (present(fix_shift)) held = fix_shift
      selection%maximum_order = maximum_order
      selection%prior_method = method_value
      selection%prior_parameter = parameter_value
      selection%fixed_shift = held
      call mixar_order_birth_death_probability(method_value, parameter_value, &
         0, birth_current, death_current, direction_status)
      if (.not. valid_model(initial_model) .or. &
         any(initial_model%degrees_of_freedom /= 0.0_dp) .or. &
         maximum_order < maxval(initial_model%order) .or. &
         maximum_order >= size(series) - 1 .or. iterations < 1 .or. &
         size(proposal_scale) /= size(initial_model%probability) .or. &
         any(proposal_scale <= 0.0_dp) .or. direction_status /= 0 .or. &
         .not. valid_order_random(random_values, size(series), &
         size(initial_model%probability), maximum_order, iterations)) then
         selection%info = 1
         return
      end if
      current = padded_mixar_model(initial_model, maximum_order)
      allocate(selection%order_trace(iterations, size(current%probability)))
      allocate(selection%probability_trace(iterations, size(current%probability)))
      allocate(selection%shift_trace(iterations, size(current%probability)))
      allocate(selection%scale_trace(iterations, size(current%probability)))
      allocate(selection%ar_trace(iterations, maximum_order, &
         size(current%probability)), source=0.0_dp)
      allocate(selection%proposed_component(iterations))
      allocate(selection%proposed_order(iterations))
      allocate(selection%accepted(iterations), source=.false.)
      allocate(component_mean(size(current%probability)))
      allocate(precision(size(current%probability)))
      do i = 1, iterations
         parameter_draw = mixar_bayesian_sample_from_random(series, current, &
            1, 1, random_values%parameter_random(i), held, &
            proposal_scale=proposal_scale)
         if (parameter_draw%info /= 0) then
            selection%info = 2
            return
         end if
         current = padded_mixar_model(parameter_draw%posterior_mean_model, &
            maximum_order)
         component_mean = parameter_draw%component_mean(1, :)
         precision = parameter_draw%precision(1, :)
         current_maximum = maxval(current%order)
         start_index = maximum_order - current_maximum + 1
         common_allocation = parameter_draw%last_allocation(start_index:)
         component = min(size(current%probability), 1 + int( &
            random_values%component_uniform(i)*real(size(current%probability), dp)))
         current_order = current%order(component)
         if (maximum_order == 0) then
            proposed_order = 0
         else if (current_order == 0) then
            proposed_order = 1
         else if (current_order == maximum_order) then
            proposed_order = maximum_order - 1
         else
            call mixar_order_birth_death_probability(method_value, parameter_value, &
               current_order, birth_current, death_current, direction_status)
            if (random_values%direction_uniform(i) < birth_current) then
               proposed_order = current_order + 1
            else
               proposed_order = current_order - 1
            end if
         end if
         selection%proposed_component(i) = component
         selection%proposed_order(i) = proposed_order
         if (proposed_order /= current_order) then
            proposal = current
            proposal%order(component) = proposed_order
            if (proposed_order > current_order) then
               proposal%ar(proposed_order, component) = &
                  3.0_dp*random_values%coefficient_uniform(i) - 1.5_dp
            else
               proposal%ar(current_order, component) = 0.0_dp
            end if
            if (.not. held) proposal%shift(component) = component_mean(component)* &
               (1.0_dp - sum(proposal%ar(:proposed_order, component)))
            current_ssr = bayesian_component_ssr(series, common_allocation, &
               component, maximum_order, &
               current%shift(component), current_order, current%ar(:, component))
            proposed_ssr = bayesian_component_ssr(series, common_allocation, &
               component, maximum_order, proposal%shift(component), &
               proposed_order, proposal%ar(:, component))
            log_ratio = -0.5_dp*precision(component)*(proposed_ssr - current_ssr)
            call mixar_order_birth_death_probability(method_value, parameter_value, &
               current_order, birth_current, death_current, direction_status)
            call mixar_order_birth_death_probability(method_value, parameter_value, &
               proposed_order, birth_proposed, death_proposed, direction_status)
            if (proposed_order > current_order) then
               log_ratio = log_ratio + log(3.0_dp) + &
                  log(max(death_proposed, tiny(1.0_dp))) - &
                  log(max(birth_current, tiny(1.0_dp)))
            else
               removed_density = exp(-0.5_dp*log_two_pi)/ &
                  proposal_scale(component)
               log_ratio = log_ratio + log(removed_density) + &
                  log(max(birth_proposed, tiny(1.0_dp))) - &
                  log(max(death_current, tiny(1.0_dp)))
            end if
            if (log(max(random_values%acceptance_uniform(i), tiny(1.0_dp))) <= &
               min(0.0_dp, log_ratio) .and. mixar_is_stable(proposal)) then
               current = proposal
               selection%accepted(i) = .true.
            end if
         end if
         selection%order_trace(i, :) = current%order
         selection%probability_trace(i, :) = current%probability
         selection%shift_trace(i, :) = current%shift
         selection%scale_trace(i, :) = current%scale
         selection%ar_trace(i, :, :) = current%ar
      end do
      selection%acceptance_rate = real(count(selection%accepted), dp)/ &
         real(iterations, dp)
      call summarize_order_selection(selection)
   end function mixar_order_select_from_random

   function mixar_order_select(series, initial_model, maximum_order, &
      iterations, proposal_scale, prior_method, prior_parameter, fix_shift) &
      result(selection)
      !! Select component AR orders using shared-RNG reversible-jump MCMC.
      real(dp), intent(in) :: series(:) !! Observed time series.
      type(mixar_model_t), intent(in) :: initial_model !! Initial Gaussian MAR model.
      integer, intent(in) :: maximum_order !! Largest permitted component order.
      integer, intent(in) :: iterations !! Number of RJMCMC iterations.
      real(dp), intent(in) :: proposal_scale(:) !! AR random-walk scales by component.
      integer, intent(in), optional :: prior_method !! Flat, ratio, or Poisson proposal method.
      real(dp), intent(in), optional :: prior_parameter !! Positive ratio or Poisson parameter.
      logical, intent(in), optional :: fix_shift !! Hold component intercepts fixed when true.
      type(mixar_order_selection_t) :: selection
      type(mixar_order_random_t) :: random_values
      integer :: method_value
      real(dp) :: parameter_value
      logical :: held

      method_value = mixar_order_prior_flat
      if (present(prior_method)) method_value = prior_method
      parameter_value = 1.0_dp
      if (present(prior_parameter)) parameter_value = prior_parameter
      held = .false.
      if (present(fix_shift)) held = fix_shift
      random_values = mixar_order_random(size(series), &
         size(initial_model%probability), maximum_order, iterations)
      selection = mixar_order_select_from_random(series, initial_model, &
         maximum_order, iterations, random_values, proposal_scale, &
         method_value, parameter_value, held)
   end function mixar_order_select

   pure subroutine summarize_order_selection(selection)
      !! Tabulate posterior order frequencies and construct the modal model.
      type(mixar_order_selection_t), intent(inout) :: selection !! RJMCMC result updated with summaries.
      integer, allocatable :: unique_work(:, :), modal(:)
      real(dp), allocatable :: frequency_work(:), probability(:), shift(:), scale(:)
      real(dp), allocatable :: ar(:, :)
      integer :: components, i, j, modal_index, states
      real(dp) :: count_modal

      components = size(selection%order_trace, 2)
      allocate(unique_work(size(selection%order_trace, 1), components), source=0)
      allocate(frequency_work(size(selection%order_trace, 1)), source=0.0_dp)
      states = 0
      do i = 1, size(selection%order_trace, 1)
         j = 0
         if (states > 0) then
            do modal_index = 1, states
               if (all(unique_work(modal_index, :) == selection%order_trace(i, :))) then
                  j = modal_index
                  exit
               end if
            end do
         end if
         if (j == 0) then
            states = states + 1
            unique_work(states, :) = selection%order_trace(i, :)
            j = states
         end if
         frequency_work(j) = frequency_work(j) + 1.0_dp
      end do
      selection%unique_order = unique_work(:states, :)
      selection%posterior_probability = frequency_work(:states)/ &
         real(size(selection%order_trace, 1), dp)
      modal_index = maxloc(selection%posterior_probability, dim=1)
      modal = selection%unique_order(modal_index, :)
      selection%modal_order = modal
      allocate(probability(components), shift(components), scale(components))
      allocate(ar(selection%maximum_order, components), source=0.0_dp)
      probability = 0.0_dp
      shift = 0.0_dp
      scale = 0.0_dp
      count_modal = 0.0_dp
      do i = 1, size(selection%order_trace, 1)
         if (.not. all(selection%order_trace(i, :) == modal)) cycle
         count_modal = count_modal + 1.0_dp
         probability = probability + selection%probability_trace(i, :)
         shift = shift + selection%shift_trace(i, :)
         scale = scale + selection%scale_trace(i, :)
         ar = ar + selection%ar_trace(i, :, :)
      end do
      probability = probability/count_modal
      shift = shift/count_modal
      scale = scale/count_modal
      ar = ar/count_modal
      selection%modal_model = mixar_model(probability/sum(probability), modal, &
         shift, scale, ar)
      if (selection%modal_model%info /= 0) selection%info = 3
   end subroutine summarize_order_selection

   pure function padded_mixar_model(input_model, maximum_order) result(model)
      !! Pad a Gaussian MAR coefficient matrix without changing its parameters.
      type(mixar_model_t), intent(in) :: input_model !! Gaussian MAR model to pad.
      integer, intent(in) :: maximum_order !! Requested coefficient row count.
      type(mixar_model_t) :: model
      real(dp), allocatable :: ar(:, :)
      integer :: copied_order

      allocate(ar(maximum_order, size(input_model%probability)), source=0.0_dp)
      copied_order = min(maximum_order, size(input_model%ar, 1))
      if (copied_order > 0) ar(:copied_order, :) = input_model%ar(:copied_order, :)
      model = mixar_model(input_model%probability, input_model%order, &
         input_model%shift, input_model%scale, ar, &
         input_model%degrees_of_freedom)
   end function padded_mixar_model

   pure logical function valid_order_random(random_values, observations, &
      components, maximum_order, iterations) result(valid)
      !! Validate dimensions and support of supplied order-selection variates.
      type(mixar_order_random_t), intent(in) :: random_values !! Supplied RJMCMC variates.
      integer, intent(in) :: observations !! Total observed series length.
      integer, intent(in) :: components !! Number of mixture components.
      integer, intent(in) :: maximum_order !! Largest permitted component order.
      integer, intent(in) :: iterations !! Number of RJMCMC iterations.
      integer :: i

      valid = allocated(random_values%parameter_random) .and. &
         allocated(random_values%component_uniform) .and. &
         allocated(random_values%direction_uniform) .and. &
         allocated(random_values%coefficient_uniform) .and. &
         allocated(random_values%acceptance_uniform)
      if (.not. valid) return
      valid = size(random_values%parameter_random) == iterations .and. &
         size(random_values%component_uniform) == iterations .and. &
         size(random_values%direction_uniform) == iterations .and. &
         size(random_values%coefficient_uniform) == iterations .and. &
         size(random_values%acceptance_uniform) == iterations .and. &
         all(random_values%component_uniform > 0.0_dp) .and. &
         all(random_values%component_uniform < 1.0_dp) .and. &
         all(random_values%direction_uniform > 0.0_dp) .and. &
         all(random_values%direction_uniform < 1.0_dp) .and. &
         all(random_values%coefficient_uniform > 0.0_dp) .and. &
         all(random_values%coefficient_uniform < 1.0_dp) .and. &
         all(random_values%acceptance_uniform > 0.0_dp) .and. &
         all(random_values%acceptance_uniform < 1.0_dp)
      if (.not. valid) return
      do i = 1, iterations
         if (.not. valid_bayesian_random(random_values%parameter_random(i), &
            observations - maximum_order, components, maximum_order, 2)) then
            valid = .false.
            return
         end if
      end do
   end function valid_order_random

   pure function mixar_initialize_from_indices(series, order, subsample_index, &
      estimate_shift, fixed_shift, minimum_scale) result(initialization)
      !! Build a stable Gaussian MAR starting model from component subsamples.
      real(dp), intent(in) :: series(:) !! Observed time series.
      integer, intent(in) :: order(:) !! Requested AR order by component.
      integer, intent(in) :: subsample_index(:, :) !! Response-time indices by subsample row and component.
      logical, intent(in), optional :: estimate_shift !! Estimate component intercepts when true.
      real(dp), intent(in), optional :: fixed_shift(:) !! Fixed component intercepts when shifts are not estimated.
      real(dp), intent(in), optional :: minimum_scale !! Innovation scale floor.
      type(mixar_initialization_t) :: initialization
      real(dp), allocatable :: ar(:, :), shift(:), scale(:), probability(:)
      real(dp), allocatable :: design(:, :), gram(:, :), inverse(:, :)
      real(dp), allocatable :: right(:, :), coefficient(:, :), feature(:)
      real(dp) :: denominator, floor_value, residual_value, weight_sum
      integer :: component, components, i, info, lag, maximum_order, observation
      integer :: q, row, usable
      logical :: fit_shift

      fit_shift = .true.
      if (present(estimate_shift)) fit_shift = estimate_shift
      floor_value = 1.0e-7_dp
      if (present(minimum_scale)) floor_value = minimum_scale
      components = size(order)
      if (components < 1 .or. any(order < 0) .or. &
         size(subsample_index, 2) /= components .or. &
         size(subsample_index, 1) < max(2, maxval(order) + merge(1, 0, fit_shift)) .or. &
         floor_value <= 0.0_dp) then
         initialization%info = 1
         return
      end if
      maximum_order = maxval(order)
      usable = size(series) - maximum_order
      if (usable < 2 .or. any(subsample_index < maximum_order + 1) .or. &
         any(subsample_index > size(series))) then
         initialization%info = 1
         return
      end if
      if (.not. fit_shift .and. present(fixed_shift)) then
         if (size(fixed_shift) /= components) then
            initialization%info = 1
            return
         end if
      end if
      allocate(ar(maximum_order, components), source=0.0_dp)
      allocate(shift(components), source=0.0_dp)
      if (.not. fit_shift .and. present(fixed_shift)) shift = fixed_shift
      allocate(scale(components), probability(components))
      allocate(initialization%residual(usable, components))
      do component = 1, components
         q = order(component) + merge(1, 0, fit_shift)
         if (q > 0) then
            allocate(design(size(subsample_index, 1), q))
            allocate(right(size(subsample_index, 1), 1))
            do row = 1, size(subsample_index, 1)
               observation = subsample_index(row, component)
               if (fit_shift) design(row, 1) = 1.0_dp
               do lag = 1, order(component)
                  design(row, lag + merge(1, 0, fit_shift)) = &
                     series(observation - lag)
               end do
               right(row, 1) = series(observation)
               if (.not. fit_shift) right(row, 1) = right(row, 1) - shift(component)
            end do
            gram = matmul(transpose(design), design)
            allocate(inverse(q, q))
            call symmetric_pseudoinverse(gram, inverse, info)
            if (info /= 0) then
               initialization%info = 2
               return
            end if
            coefficient = matmul(inverse, matmul(transpose(design), right))
            if (fit_shift) shift(component) = coefficient(1, 1)
            if (order(component) > 0) ar(:order(component), component) = &
               coefficient(1 + merge(1, 0, fit_shift):q, 1)
            deallocate(design, right, gram, inverse, coefficient)
         end if
         do row = 1, usable
            observation = maximum_order + row
            residual_value = series(observation) - shift(component)
            do lag = 1, order(component)
               residual_value = residual_value - ar(lag, component)* &
                  series(observation - lag)
            end do
            initialization%residual(row, component) = residual_value
         end do
      end do
      allocate(initialization%responsibility(usable, components))
      if (components == 1) then
         initialization%responsibility = 1.0_dp
      else
         do row = 1, usable
            denominator = sum(abs(initialization%residual(row, :)))
            if (denominator <= tiny(1.0_dp)) then
               initialization%responsibility(row, :) = 1.0_dp/real(components, dp)
            else
               initialization%responsibility(row, :) = &
                  (1.0_dp - abs(initialization%residual(row, :))/denominator)/ &
                  real(components - 1, dp)
            end if
         end do
      end if
      do component = 1, components
         weight_sum = sum(initialization%responsibility(:, component))
         if (weight_sum <= tiny(1.0_dp)) then
            initialization%info = 3
            return
         end if
         probability(component) = weight_sum/real(usable, dp)
         q = order(component) + merge(1, 0, fit_shift)
         if (q > 0) then
            allocate(gram(q, q), right(q, 1), feature(q))
            gram = 0.0_dp
            right = 0.0_dp
            do row = 1, usable
               observation = maximum_order + row
               if (fit_shift) feature(1) = 1.0_dp
               do lag = 1, order(component)
                  feature(lag + merge(1, 0, fit_shift)) = series(observation - lag)
               end do
               residual_value = series(observation)
               if (.not. fit_shift) residual_value = residual_value - shift(component)
               gram = gram + initialization%responsibility(row, component)* &
                  spread(feature, 2, q)*spread(feature, 1, q)
               right(:, 1) = right(:, 1) + &
                  initialization%responsibility(row, component)*feature*residual_value
            end do
            allocate(inverse(q, q))
            call symmetric_pseudoinverse(gram, inverse, info)
            if (info /= 0) then
               initialization%info = 4
               return
            end if
            coefficient = matmul(inverse, right)
            if (fit_shift) shift(component) = coefficient(1, 1)
            if (order(component) > 0) ar(:order(component), component) = &
               coefficient(1 + merge(1, 0, fit_shift):q, 1)
            deallocate(gram, right, feature, inverse, coefficient)
         end if
         do row = 1, usable
            observation = maximum_order + row
            residual_value = series(observation) - shift(component)
            do lag = 1, order(component)
               residual_value = residual_value - ar(lag, component)* &
                  series(observation - lag)
            end do
            initialization%residual(row, component) = residual_value
         end do
         scale(component) = max(floor_value, sqrt(dot_product( &
            initialization%responsibility(:, component), &
            initialization%residual(:, component)**2)/weight_sum))
      end do
      initialization%model = mixar_model(probability/sum(probability), order, &
         shift, scale, ar)
      do while (.not. mixar_is_stable(initialization%model) .and. &
         initialization%stability_contractions < 200)
         initialization%model%ar = 0.9_dp*initialization%model%ar
         initialization%stability_contractions = &
            initialization%stability_contractions + 1
      end do
      if (.not. mixar_is_stable(initialization%model)) then
         initialization%info = 5
         return
      end if
      do component = 1, components
         do row = 1, usable
            observation = maximum_order + row
            residual_value = series(observation) - &
               initialization%model%shift(component)
            do lag = 1, order(component)
               residual_value = residual_value - &
                  initialization%model%ar(lag, component)* &
                  series(observation - lag)
            end do
            initialization%residual(row, component) = residual_value
         end do
         weight_sum = sum(initialization%responsibility(:, component))
         initialization%model%scale(component) = max(floor_value, &
            sqrt(dot_product(initialization%responsibility(:, component), &
            initialization%residual(:, component)**2)/weight_sum))
      end do
      initialization%subsample_index = subsample_index
   end function mixar_initialize_from_indices

   function mixar_random_subsample_indices(series_length, maximum_order, &
      components, starts, subsample_size) result(indices)
      !! Draw component subsamples without replacement using the shared RNG.
      integer, intent(in) :: series_length !! Total observed series length.
      integer, intent(in) :: maximum_order !! Largest requested component order.
      integer, intent(in) :: components !! Number of mixture components.
      integer, intent(in) :: starts !! Number of independent starts.
      integer, intent(in) :: subsample_size !! Observations in each component subsample.
      integer, allocatable :: indices(:, :, :)
      integer, allocatable :: pool(:)
      integer :: available, component, held, i, selected, start

      available = series_length - maximum_order
      if (available < 1 .or. components < 1 .or. starts < 1 .or. &
         subsample_size < 1 .or. subsample_size > available) then
         allocate(indices(0, 0, 0))
         return
      end if
      allocate(indices(subsample_size, components, starts))
      allocate(pool(available))
      do start = 1, starts
         do component = 1, components
            pool = [(maximum_order + i, i=1, available)]
            do i = 1, subsample_size
               selected = i + int(random_uniform()*real(available - i + 1, dp))
               selected = min(available, selected)
               held = pool(i)
               pool(i) = pool(selected)
               pool(selected) = held
               indices(i, component, start) = pool(i)
            end do
         end do
      end do
   end function mixar_random_subsample_indices

   function mixar_random_initialize(series, order, subsample_size, &
      estimate_shift, fixed_shift, minimum_scale) result(initialization)
      !! Build one stable Gaussian MAR start from shared-RNG subsamples.
      real(dp), intent(in) :: series(:) !! Observed time series.
      integer, intent(in) :: order(:) !! Requested AR order by component.
      integer, intent(in), optional :: subsample_size !! Observations per component subsample.
      logical, intent(in), optional :: estimate_shift !! Estimate component intercepts when true.
      real(dp), intent(in), optional :: fixed_shift(:) !! Fixed component intercepts.
      real(dp), intent(in), optional :: minimum_scale !! Innovation scale floor.
      type(mixar_initialization_t) :: initialization
      integer, allocatable :: indices(:, :, :)
      integer :: requested_size
      logical :: fit_shift
      real(dp) :: floor_value
      real(dp), allocatable :: held_shift(:)

      if (size(order) < 1 .or. any(order < 0)) then
         initialization%info = 1
         return
      end if
      fit_shift = .true.
      if (present(estimate_shift)) fit_shift = estimate_shift
      requested_size = max(10, maxval(order) + merge(1, 0, fit_shift))
      if (present(subsample_size)) requested_size = subsample_size
      floor_value = 1.0e-7_dp
      if (present(minimum_scale)) floor_value = minimum_scale
      allocate(held_shift(size(order)), source=0.0_dp)
      if (present(fixed_shift)) then
         if (size(fixed_shift) /= size(order)) then
            initialization%info = 1
            return
         end if
         held_shift = fixed_shift
      end if
      indices = mixar_random_subsample_indices(size(series), maxval(order), &
         size(order), 1, requested_size)
      if (size(indices) == 0) then
         initialization%info = 1
         return
      end if
      initialization = mixar_initialize_from_indices(series, order, &
         indices(:, :, 1), fit_shift, held_shift, floor_value)
   end function mixar_random_initialize

   pure function mixar_multistart_fit(series, order, subsample_index, &
      estimate_shift, fixed_shift, max_iterations, tolerance, minimum_scale) &
      result(out)
      !! Fit supplied subsample starts and retain the best converged Gaussian MAR model.
      real(dp), intent(in) :: series(:) !! Observed time series.
      integer, intent(in) :: order(:) !! Requested AR order by component.
      integer, intent(in) :: subsample_index(:, :, :) !! Subsample indices by row, component, and start.
      logical, intent(in), optional :: estimate_shift !! Estimate component intercepts when true.
      real(dp), intent(in), optional :: fixed_shift(:) !! Fixed component intercepts.
      integer, intent(in), optional :: max_iterations !! Maximum EM iterations per start.
      real(dp), intent(in), optional :: tolerance !! Relative likelihood tolerance.
      real(dp), intent(in), optional :: minimum_scale !! Innovation scale floor.
      type(mixar_multistart_t) :: out
      real(dp), allocatable :: held_shift(:)
      real(dp) :: floor_value, tolerance_value
      integer :: iteration_limit, start, starts
      logical :: fit_shift

      starts = size(subsample_index, 3)
      fit_shift = .true.
      if (present(estimate_shift)) fit_shift = estimate_shift
      iteration_limit = 200
      if (present(max_iterations)) iteration_limit = max_iterations
      tolerance_value = 1.0e-10_dp
      if (present(tolerance)) tolerance_value = tolerance
      floor_value = 1.0e-7_dp
      if (present(minimum_scale)) floor_value = minimum_scale
      allocate(held_shift(size(order)), source=0.0_dp)
      if (present(fixed_shift)) then
         if (size(fixed_shift) /= size(order)) then
            out%info = 1
            return
         end if
         held_shift = fixed_shift
      end if
      if (starts < 1 .or. size(subsample_index, 2) /= size(order) .or. &
         iteration_limit < 1 .or. tolerance_value <= 0.0_dp .or. &
         floor_value <= 0.0_dp) then
         out%info = 1
         return
      end if
      allocate(out%initialization(starts), out%fit(starts))
      allocate(out%initial_log_likelihood(starts), source=-huge(1.0_dp))
      do start = 1, starts
         out%initialization(start) = mixar_initialize_from_indices(series, order, &
            subsample_index(:, :, start), fit_shift, held_shift, floor_value)
         if (out%initialization(start)%info /= 0) cycle
         out%initial_log_likelihood(start) = mixar_log_likelihood( &
            out%initialization(start)%model, series)
         out%fit(start) = mixar_fit(series, out%initialization(start)%model, &
            fit_shift, iteration_limit, tolerance_value, floor_value)
         if (out%fit(start)%info /= 0 .or. .not. out%fit(start)%converged) cycle
         out%successful = out%successful + 1
         if (out%best == 0) then
            out%best = start
         else if (out%fit(start)%log_likelihood > &
            out%fit(out%best)%log_likelihood) then
            out%best = start
         end if
      end do
      if (out%best == 0) out%info = 2
   end function mixar_multistart_fit

   function mixar_random_multistart_fit(series, order, starts, subsample_size, &
      estimate_shift, fixed_shift, max_iterations, tolerance, minimum_scale) &
      result(out)
      !! Generate shared-RNG subsample starts and select the best converged EM fit.
      real(dp), intent(in) :: series(:) !! Observed time series.
      integer, intent(in) :: order(:) !! Requested AR order by component.
      integer, intent(in) :: starts !! Number of random EM starts.
      integer, intent(in), optional :: subsample_size !! Observations per component subsample.
      logical, intent(in), optional :: estimate_shift !! Estimate component intercepts when true.
      real(dp), intent(in), optional :: fixed_shift(:) !! Fixed component intercepts.
      integer, intent(in), optional :: max_iterations !! Maximum EM iterations per start.
      real(dp), intent(in), optional :: tolerance !! Relative likelihood tolerance.
      real(dp), intent(in), optional :: minimum_scale !! Innovation scale floor.
      type(mixar_multistart_t) :: out
      integer, allocatable :: indices(:, :, :)
      integer :: iteration_limit, requested_size
      logical :: fit_shift
      real(dp) :: floor_value, tolerance_value
      real(dp), allocatable :: held_shift(:)

      if (size(order) < 1 .or. any(order < 0) .or. starts < 1) then
         out%info = 1
         return
      end if
      fit_shift = .true.
      if (present(estimate_shift)) fit_shift = estimate_shift
      requested_size = max(10, maxval(order) + merge(1, 0, fit_shift))
      if (present(subsample_size)) requested_size = subsample_size
      iteration_limit = 200
      if (present(max_iterations)) iteration_limit = max_iterations
      tolerance_value = 1.0e-10_dp
      if (present(tolerance)) tolerance_value = tolerance
      floor_value = 1.0e-7_dp
      if (present(minimum_scale)) floor_value = minimum_scale
      allocate(held_shift(size(order)), source=0.0_dp)
      if (present(fixed_shift)) then
         if (size(fixed_shift) /= size(order)) then
            out%info = 1
            return
         end if
         held_shift = fixed_shift
      end if
      indices = mixar_random_subsample_indices(size(series), maxval(order), &
         size(order), starts, requested_size)
      if (size(indices) == 0) then
         out%info = 1
         return
      end if
      out = mixar_multistart_fit(series, order, indices, fit_shift, held_shift, &
         iteration_limit, tolerance_value, floor_value)
   end function mixar_random_multistart_fit

   pure function mixar_standard_absolute_moment(degrees_of_freedom, order) &
      result(moment)
      !! Return an absolute moment of a standard normal or unit-variance Student-t variable.
      real(dp), intent(in) :: degrees_of_freedom !! Zero for Gaussian or Student-t degrees above two.
      integer, intent(in) :: order !! Nonnegative integer moment order.
      type(mixar_moment_t) :: moment
      real(dp) :: log_value

      moment%order = order
      if (order < 0 .or. (degrees_of_freedom /= 0.0_dp .and. &
         degrees_of_freedom <= 2.0_dp)) then
         moment%info = 1
         return
      end if
      if (degrees_of_freedom > 0.0_dp .and. &
         real(order, dp) >= degrees_of_freedom) then
         moment%value = huge(1.0_dp)
         return
      end if
      if (degrees_of_freedom == 0.0_dp) then
         log_value = 0.5_dp*real(order, dp)*log(2.0_dp) + &
            log_gamma(0.5_dp*real(order + 1, dp)) - &
            0.5_dp*log(acos(-1.0_dp))
      else
         log_value = 0.5_dp*real(order, dp)*log(degrees_of_freedom - 2.0_dp) + &
            log_gamma(0.5_dp*real(order + 1, dp)) + &
            log_gamma(0.5_dp*(degrees_of_freedom - real(order, dp))) - &
            0.5_dp*log(acos(-1.0_dp)) - &
            log_gamma(0.5_dp*degrees_of_freedom)
      end if
      moment%exists = .true.
      if (log_value >= log(huge(1.0_dp))) then
         moment%value = huge(1.0_dp)
         moment%info = 2
      else
         moment%value = exp(log_value)
      end if
   end function mixar_standard_absolute_moment

   pure function mixar_standard_moment(degrees_of_freedom, order) result(moment)
      !! Return an ordinary moment of a standard normal or unit-variance Student-t variable.
      real(dp), intent(in) :: degrees_of_freedom !! Zero for Gaussian or Student-t degrees above two.
      integer, intent(in) :: order !! Nonnegative integer moment order.
      type(mixar_moment_t) :: moment

      moment = mixar_standard_absolute_moment(degrees_of_freedom, order)
      if (moment%exists .and. mod(order, 2) == 1) moment%value = 0.0_dp
   end function mixar_standard_moment

   pure function mixar_innovation_moment(model, order) result(moment)
      !! Return a moment of the model's unconditional innovation mixture.
      type(mixar_model_t), intent(in) :: model !! Mixture autoregression parameters.
      integer, intent(in) :: order !! Nonnegative integer moment order.
      type(mixar_moment_t) :: moment
      type(mixar_moment_t) :: component_moment
      integer :: component

      moment%order = order
      if (.not. valid_model(model) .or. order < 0) then
         moment%info = 1
         return
      end if
      moment%exists = .true.
      do component = 1, size(model%probability)
         component_moment = mixar_standard_moment( &
            model%degrees_of_freedom(component), order)
         if (.not. component_moment%exists) then
            moment%exists = .false.
            moment%value = huge(1.0_dp)
            return
         end if
         if (component_moment%info /= 0) then
            moment%info = component_moment%info
            return
         end if
         moment%value = moment%value + model%probability(component)* &
            model%scale(component)**order*component_moment%value
      end do
   end function mixar_innovation_moment

   pure function mixar_conditional_moment(model, history, order) result(moment)
      !! Return a raw one-step conditional moment after the supplied history.
      type(mixar_model_t), intent(in) :: model !! Mixture autoregression parameters.
      real(dp), intent(in) :: history(:) !! Observed history ending at the forecast origin.
      integer, intent(in) :: order !! Nonnegative integer moment order.
      type(mixar_moment_t) :: moment

      moment = conditional_moment_kernel(model, history, order, .false.)
   end function mixar_conditional_moment

   pure function mixar_conditional_central_moment(model, history, order) &
      result(moment)
      !! Return a central one-step conditional moment after the supplied history.
      type(mixar_model_t), intent(in) :: model !! Mixture autoregression parameters.
      real(dp), intent(in) :: history(:) !! Observed history ending at the forecast origin.
      integer, intent(in) :: order !! Nonnegative integer moment order.
      type(mixar_moment_t) :: moment

      moment = conditional_moment_kernel(model, history, order, .true.)
   end function mixar_conditional_central_moment

   pure function mixar_conditional_kurtosis(model, history) result(moment)
      !! Return one-step conditional kurtosis after the supplied history.
      type(mixar_model_t), intent(in) :: model !! Mixture autoregression parameters.
      real(dp), intent(in) :: history(:) !! Observed history ending at the forecast origin.
      type(mixar_moment_t) :: moment
      type(mixar_moment_t) :: second, fourth

      moment%order = 4
      second = mixar_conditional_central_moment(model, history, 2)
      fourth = mixar_conditional_central_moment(model, history, 4)
      if (second%info /= 0 .or. fourth%info /= 0) then
         moment%info = max(second%info, fourth%info)
         return
      end if
      if (.not. second%exists .or. .not. fourth%exists) return
      if (second%value <= tiny(1.0_dp)) then
         moment%info = 2
         return
      end if
      moment%value = fourth%value/second%value**2
      moment%exists = .true.
   end function mixar_conditional_kurtosis

   pure function mixar_conditional_excess_kurtosis(model, history) result(moment)
      !! Return one-step conditional excess kurtosis after the supplied history.
      type(mixar_model_t), intent(in) :: model !! Mixture autoregression parameters.
      real(dp), intent(in) :: history(:) !! Observed history ending at the forecast origin.
      type(mixar_moment_t) :: moment

      moment = mixar_conditional_kurtosis(model, history)
      if (moment%exists) moment%value = moment%value - 3.0_dp
   end function mixar_conditional_excess_kurtosis

   pure function mixar_stationary_moments(model) result(moments)
      !! Solve stable companion-state equations for stationary mean and variance.
      type(mixar_model_t), intent(in) :: model !! Stable mixture autoregression parameters.
      type(mixar_stationary_moments_t) :: moments
      real(dp), allocatable :: aggregate(:, :), companion(:, :), intercept(:)
      real(dp), allocatable :: mean_matrix(:, :), mean_right(:, :), solution(:, :)
      real(dp), allocatable :: second_matrix(:, :), second_right(:, :)
      real(dp), allocatable :: raw_second(:, :), innovation(:, :), transformed_mean(:)
      real(dp), allocatable :: contribution(:, :)
      integer :: component, dimension, i, info, maximum_order, squared_dimension

      if (.not. valid_model(model) .or. .not. mixar_is_stable(model)) then
         moments%info = 1
         return
      end if
      maximum_order = maxval(model%order)
      dimension = max(1, maximum_order)
      allocate(aggregate(dimension, dimension), source=0.0_dp)
      allocate(mean_right(dimension, 1), source=0.0_dp)
      allocate(companion(dimension, dimension), intercept(dimension))
      do component = 1, size(model%probability)
         companion = 0.0_dp
         if (model%order(component) > 0) companion(1, &
            :model%order(component)) = model%ar(:model%order(component), component)
         do i = 2, dimension
            companion(i, i - 1) = 1.0_dp
         end do
         intercept = 0.0_dp
         intercept(1) = model%shift(component)
         aggregate = aggregate + model%probability(component)*companion
         mean_right(:, 1) = mean_right(:, 1) + &
            model%probability(component)*intercept
      end do
      mean_matrix = identity_block(dimension) - aggregate
      allocate(solution(dimension, 1))
      call solve_matrix(mean_matrix, mean_right, solution, info)
      if (info /= 0) then
         moments%info = 2
         return
      end if
      allocate(moments%state_mean(dimension))
      moments%state_mean = solution(:, 1)
      squared_dimension = dimension*dimension
      deallocate(aggregate)
      allocate(aggregate(squared_dimension, squared_dimension), source=0.0_dp)
      allocate(second_right(squared_dimension, 1), source=0.0_dp)
      allocate(innovation(dimension, dimension), contribution(dimension, dimension))
      allocate(transformed_mean(dimension))
      do component = 1, size(model%probability)
         companion = 0.0_dp
         if (model%order(component) > 0) companion(1, &
            :model%order(component)) = model%ar(:model%order(component), component)
         do i = 2, dimension
            companion(i, i - 1) = 1.0_dp
         end do
         intercept = 0.0_dp
         intercept(1) = model%shift(component)
         innovation = 0.0_dp
         innovation(1, 1) = model%scale(component)**2
         transformed_mean = matmul(companion, moments%state_mean)
         contribution = spread(intercept, 2, dimension)* &
            spread(intercept, 1, dimension) + &
            spread(transformed_mean, 2, dimension)* &
            spread(intercept, 1, dimension) + &
            spread(intercept, 2, dimension)* &
            spread(transformed_mean, 1, dimension) + innovation
         aggregate = aggregate + model%probability(component)* &
            kronecker_product(companion, companion)
         second_right(:, 1) = second_right(:, 1) + &
            model%probability(component)*reshape(contribution, [squared_dimension])
      end do
      second_matrix = identity_block(squared_dimension) - aggregate
      deallocate(solution)
      allocate(solution(squared_dimension, 1))
      call solve_matrix(second_matrix, second_right, solution, info)
      if (info /= 0) then
         moments%info = 3
         return
      end if
      raw_second = reshape(solution(:, 1), [dimension, dimension])
      allocate(moments%state_covariance(dimension, dimension))
      moments%state_covariance = raw_second - &
         spread(moments%state_mean, 2, dimension)* &
         spread(moments%state_mean, 1, dimension)
      moments%state_covariance = 0.5_dp*(moments%state_covariance + &
         transpose(moments%state_covariance))
      moments%mean = moments%state_mean(1)
      moments%variance = moments%state_covariance(1, 1)
      if (moments%variance < -1.0e-10_dp) then
         moments%info = 4
      else
         moments%variance = max(0.0_dp, moments%variance)
      end if
   end function mixar_stationary_moments

   pure function conditional_moment_kernel(model, history, order, central) &
      result(moment)
      !! Evaluate a raw or central moment of the one-step mixture distribution.
      type(mixar_model_t), intent(in) :: model !! Mixture autoregression parameters.
      real(dp), intent(in) :: history(:) !! Observed history ending at the forecast origin.
      integer, intent(in) :: order !! Nonnegative integer moment order.
      logical, intent(in) :: central !! Center component locations on the mixture mean.
      type(mixar_moment_t) :: moment
      type(mixar_moment_t) :: innovation_moment
      real(dp), allocatable :: location(:)
      real(dp) :: center, coefficient
      integer :: component, j, lag

      moment%order = order
      if (.not. valid_model(model) .or. order < 0 .or. &
         size(history) < maxval(model%order)) then
         moment%info = 1
         return
      end if
      allocate(location(size(model%probability)))
      do component = 1, size(model%probability)
         location(component) = model%shift(component)
         do lag = 1, model%order(component)
            location(component) = location(component) + model%ar(lag, component)* &
               history(size(history) - lag + 1)
         end do
      end do
      center = 0.0_dp
      if (central) center = dot_product(model%probability, location)
      moment%exists = .true.
      do component = 1, size(model%probability)
         coefficient = 1.0_dp
         do j = 0, order
            innovation_moment = mixar_standard_moment( &
               model%degrees_of_freedom(component), j)
            if (.not. innovation_moment%exists) then
               moment%exists = .false.
               moment%value = huge(1.0_dp)
               return
            end if
            if (innovation_moment%info /= 0) then
               moment%info = innovation_moment%info
               return
            end if
            if (j > 0) coefficient = coefficient* &
               real(order - j + 1, dp)/real(j, dp)
            moment%value = moment%value + model%probability(component)* &
               coefficient*(location(component) - center)**(order - j)* &
               model%scale(component)**j*innovation_moment%value
         end do
      end do
   end function conditional_moment_kernel

   pure function bayesian_posterior_mean_model(draws) result(model)
      !! Construct a Gaussian MAR model from posterior parameter means.
      type(mixar_bayesian_draws_t), intent(in) :: draws !! Retained Bayesian MAR draws.
      type(mixar_model_t) :: model
      real(dp), allocatable :: ar(:, :), probability(:), scale(:), shift(:)
      real(dp) :: divisor

      divisor = real(size(draws%probability, 1), dp)
      probability = sum(draws%probability, dim=1)/divisor
      scale = sum(draws%scale, dim=1)/divisor
      shift = sum(draws%shift, dim=1)/divisor
      ar = sum(draws%ar, dim=1)/divisor
      model = mixar_model(probability/sum(probability), draws%order, shift, scale, ar)
   end function bayesian_posterior_mean_model

   pure logical function valid_bayesian_random(random_values, observations, &
      components, maximum_order, iterations) result(valid)
      !! Validate dimensions and support of supplied Bayesian sampler variates.
      type(mixar_bayesian_random_t), intent(in) :: random_values !! Supplied independent variates.
      integer, intent(in) :: observations !! Number of conditional observations.
      integer, intent(in) :: components !! Number of mixture components.
      integer, intent(in) :: maximum_order !! Largest component AR order.
      integer, intent(in) :: iterations !! Total sampler iterations.

      valid = allocated(random_values%allocation_uniform) .and. &
         allocated(random_values%probability_uniform) .and. &
         allocated(random_values%mean_normal) .and. &
         allocated(random_values%lambda_uniform) .and. &
         allocated(random_values%precision_uniform) .and. &
         allocated(random_values%ar_normal) .and. &
         allocated(random_values%acceptance_uniform)
      if (.not. valid) return
      valid = size(random_values%allocation_uniform, 1) >= observations .and. &
         size(random_values%allocation_uniform, 2) >= iterations .and. &
         size(random_values%probability_uniform, 1) >= components .and. &
         size(random_values%probability_uniform, 2) >= iterations .and. &
         size(random_values%mean_normal, 1) >= components .and. &
         size(random_values%mean_normal, 2) >= iterations .and. &
         size(random_values%lambda_uniform) >= iterations .and. &
         size(random_values%precision_uniform, 1) >= components .and. &
         size(random_values%precision_uniform, 2) >= iterations .and. &
         size(random_values%ar_normal, 1) >= maximum_order .and. &
         size(random_values%ar_normal, 2) >= components .and. &
         size(random_values%ar_normal, 3) >= iterations .and. &
         size(random_values%acceptance_uniform, 1) >= components .and. &
         size(random_values%acceptance_uniform, 2) >= iterations .and. &
         all(random_values%allocation_uniform > 0.0_dp) .and. &
         all(random_values%allocation_uniform < 1.0_dp) .and. &
         all(random_values%probability_uniform > 0.0_dp) .and. &
         all(random_values%probability_uniform < 1.0_dp) .and. &
         all(random_values%lambda_uniform > 0.0_dp) .and. &
         all(random_values%lambda_uniform < 1.0_dp) .and. &
         all(random_values%precision_uniform > 0.0_dp) .and. &
         all(random_values%precision_uniform < 1.0_dp) .and. &
         all(random_values%acceptance_uniform > 0.0_dp) .and. &
         all(random_values%acceptance_uniform < 1.0_dp)
   end function valid_bayesian_random

   pure real(dp) function bayesian_component_ssr(series, allocation, component, &
      maximum_order, shift, order, ar) result(ssr)
      !! Return one allocated component's autoregressive residual sum of squares.
      real(dp), intent(in) :: series(:) !! Observed time series.
      integer, intent(in) :: allocation(:) !! Component allocation by usable observation.
      integer, intent(in) :: component !! Selected component index.
      integer, intent(in) :: maximum_order !! Largest model order.
      real(dp), intent(in) :: shift !! Component intercept.
      integer, intent(in) :: order !! Selected component order.
      real(dp), intent(in) :: ar(:) !! Padded component AR coefficients.
      real(dp) :: residual
      integer :: i, lag, observation

      ssr = 0.0_dp
      do i = 1, size(allocation)
         if (allocation(i) /= component) cycle
         observation = maximum_order + i
         residual = series(observation) - shift
         do lag = 1, order
            residual = residual - ar(lag)*series(observation - lag)
         end do
         ssr = ssr + residual**2
      end do
   end function bayesian_component_ssr

   pure real(dp) function bayesian_unshifted_error_sum(series, allocation, &
      component, maximum_order, order, ar) result(error_sum)
      !! Sum allocated responses after removing autoregressive lag terms.
      real(dp), intent(in) :: series(:) !! Observed time series.
      integer, intent(in) :: allocation(:) !! Component allocation by usable observation.
      integer, intent(in) :: component !! Selected component index.
      integer, intent(in) :: maximum_order !! Largest model order.
      integer, intent(in) :: order !! Selected component order.
      real(dp), intent(in) :: ar(:) !! Padded component AR coefficients.
      real(dp) :: residual
      integer :: i, lag, observation

      error_sum = 0.0_dp
      do i = 1, size(allocation)
         if (allocation(i) /= component) cycle
         observation = maximum_order + i
         residual = series(observation)
         do lag = 1, order
            residual = residual - ar(lag)*series(observation - lag)
         end do
         error_sum = error_sum + residual
      end do
   end function bayesian_unshifted_error_sum

   pure real(dp) function gamma_quantile(probability, shape) result(value)
      !! Invert a unit-rate Gamma distribution by bracketed bisection.
      real(dp), intent(in) :: probability !! Probability strictly between zero and one.
      real(dp), intent(in) :: shape !! Gamma shape parameter.
      real(dp) :: lower, upper, midpoint
      integer :: iteration

      lower = 0.0_dp
      upper = max(1.0_dp, shape)
      do while (1.0_dp - regularized_gamma_q(shape, upper) < probability)
         upper = 2.0_dp*upper
         if (upper >= huge(1.0_dp)/4.0_dp) exit
      end do
      do iteration = 1, 80
         midpoint = 0.5_dp*(lower + upper)
         if (1.0_dp - regularized_gamma_q(shape, midpoint) < probability) then
            lower = midpoint
         else
            upper = midpoint
         end if
      end do
      value = 0.5_dp*(lower + upper)
   end function gamma_quantile

   pure real(dp) function gamma_log_density(value, shape, rate) result(log_density)
      !! Evaluate a shape-rate Gamma log density.
      real(dp), intent(in) :: value !! Positive evaluation point.
      real(dp), intent(in) :: shape !! Gamma shape parameter.
      real(dp), intent(in) :: rate !! Gamma rate parameter.

      if (value <= 0.0_dp .or. shape <= 0.0_dp .or. rate <= 0.0_dp) then
         log_density = -huge(1.0_dp)
      else
         log_density = shape*log(rate) - log_gamma(shape) + &
            (shape - 1.0_dp)*log(value) - rate*value
      end if
   end function gamma_log_density

   pure real(dp) function dirichlet_log_density(probability, alpha) result(log_density)
      !! Evaluate a Dirichlet log density on the probability simplex.
      real(dp), intent(in) :: probability(:) !! Probability vector.
      real(dp), intent(in) :: alpha(:) !! Positive concentration parameters.
      integer :: component

      if (size(probability) /= size(alpha) .or. any(probability <= 0.0_dp) .or. &
         any(alpha <= 0.0_dp)) then
         log_density = -huge(1.0_dp)
         return
      end if
      log_density = log_gamma(sum(alpha))
      do component = 1, size(alpha)
         log_density = log_density - log_gamma(alpha(component)) + &
            (alpha(component) - 1.0_dp)*log(probability(component))
      end do
   end function dirichlet_log_density

   pure real(dp) function kde_log_at(values, point) result(log_density)
      !! Estimate a univariate log density with a Gaussian kernel.
      real(dp), intent(in) :: values(:) !! Posterior sample values.
      real(dp), intent(in) :: point !! Evaluation point.
      real(dp) :: bandwidth, sample_standard_deviation, maximum_log
      real(dp), allocatable :: log_kernel(:)

      sample_standard_deviation = sqrt(max(variance(values), 0.0_dp))
      bandwidth = max(1.06_dp*sample_standard_deviation* &
         real(size(values), dp)**(-0.2_dp), 1.0e-8_dp*max(1.0_dp, abs(point)))
      log_kernel = -0.5_dp*((values - point)/bandwidth)**2
      maximum_log = maxval(log_kernel)
      log_density = maximum_log + log(sum(exp(log_kernel - maximum_log))) - &
         log(real(size(values), dp)*bandwidth) - 0.5_dp*log_two_pi
   end function kde_log_at

   pure integer function find_integer(values, target) result(position)
      !! Return the position of an integer in a permutation vector.
      integer, intent(in) :: values(:) !! Permutation vector.
      integer, intent(in) :: target !! Integer value to locate.
      integer :: i

      position = 0
      do i = 1, size(values)
         if (values(i) == target) then
            position = i
            return
         end if
      end do
   end function find_integer

   pure logical function valid_model(model) result(valid)
      !! Check dimensions and constraints of a Gaussian MAR model.
      type(mixar_model_t), intent(in) :: model !! Gaussian MAR model.
      integer :: components

      valid = .false.
      if (.not. allocated(model%probability) .or. .not. allocated(model%order) .or. &
         .not. allocated(model%shift) .or. .not. allocated(model%scale) .or. &
         .not. allocated(model%ar) .or. &
         .not. allocated(model%degrees_of_freedom)) return
      components = size(model%probability)
      if (components < 1 .or. size(model%order) /= components .or. &
         size(model%shift) /= components .or. size(model%scale) /= components .or. &
         size(model%ar, 2) /= components .or. &
         size(model%degrees_of_freedom) /= components .or. any(model%order < 0)) return
      if (maxval(model%order) > size(model%ar, 1) .or. &
         any(model%probability < 0.0_dp) .or. any(model%scale <= 0.0_dp) .or. &
         abs(sum(model%probability) - 1.0_dp) > 1.0e-10_dp) return
      if (any(model%degrees_of_freedom /= 0.0_dp .and. &
         model%degrees_of_freedom <= 2.0_dp)) return
      valid = all(ieee_is_finite(model%probability)) .and. &
         all(ieee_is_finite(model%shift)) .and. all(ieee_is_finite(model%scale)) .and. &
         all(ieee_is_finite(model%ar)) .and. &
         all(ieee_is_finite(model%degrees_of_freedom))
   end function valid_model

   pure logical function valid_regression_model(model) result(valid)
      !! Check dimensions and constraints of a MAR regression model.
      type(mixar_regression_model_t), intent(in) :: model !! MAR regression model.
      integer :: columns, components

      valid = .false.
      if (.not. valid_model(model%mar) .or. &
         .not. allocated(model%coefficient)) return
      columns = size(model%coefficient, 2)
      components = size(model%mar%probability)
      if (size(model%coefficient, 1) < 1 .or. &
         .not. all(ieee_is_finite(model%coefficient))) return
      if (model%component_specific) then
         if (columns /= components) return
      else
         if (columns /= 1) return
      end if
      valid = .true.
   end function valid_regression_model

   pure logical function valid_seasonal_model(model) result(valid)
      !! Check dimensions and constraints of an additive seasonal MAR model.
      type(mixar_seasonal_model_t), intent(in) :: model !! Seasonal MAR model.
      integer :: components

      valid = .false.
      if (.not. valid_model(model%mar) .or. model%seasonal_period < 1 .or. &
         .not. allocated(model%seasonal_order) .or. &
         .not. allocated(model%seasonal_ar)) return
      components = size(model%mar%probability)
      if (size(model%seasonal_order) /= components .or. &
         size(model%seasonal_ar, 2) /= components .or. &
         any(model%seasonal_order < 0) .or. &
         maxval(model%seasonal_order) > size(model%seasonal_ar, 1) .or. &
         .not. all(ieee_is_finite(model%seasonal_ar))) return
      valid = .true.
   end function valid_seasonal_model

   subroutine display_mixar_model(model)
      !! Display Gaussian MAR parameters without package-internal storage details.
      type(mixar_model_t), intent(in) :: model !! Gaussian MAR model.
      integer :: component

      print '(a,i0)', 'Gaussian MAR components: ', size(model%probability)
      do component = 1, size(model%probability)
         print '(a,i0,a,f10.5,a,i0,a,f10.5,a,f10.5)', 'component ', component, &
            ': probability=', model%probability(component), ', order=', model%order(component), &
            ', shift=', model%shift(component), ', scale=', model%scale(component)
         if (model%degrees_of_freedom(component) > 2.0_dp) &
            print '(a,f10.4)', '  standardized Student-t degrees of freedom: ', &
            model%degrees_of_freedom(component)
         if (model%order(component) > 0) print '(a,*(f10.5,1x))', &
            '  AR: ', model%ar(1:model%order(component), component)
      end do
   end subroutine display_mixar_model

   subroutine display_mixar_fit(fit)
      !! Display Gaussian MAR fit convergence and parameter estimates.
      type(mixar_fit_t), intent(in) :: fit !! Fitted Gaussian MAR result.

      print '(a,l1)', 'Converged: ', fit%converged
      print '(a,i0)', 'Iterations: ', fit%iterations
      print '(a,f14.5)', 'Log likelihood: ', fit%log_likelihood
      print '(a,f14.5)', 'AIC: ', fit%aic
      print '(a,f14.5)', 'BIC: ', fit%bic
      call display_mixar_model(fit%model)
   end subroutine display_mixar_fit

   subroutine display_mixar_simulation(simulation, print_obs)
      !! Display simulation parameters and optionally the generated observations.
      type(mixar_simulation_t), intent(in) :: simulation !! Gaussian MAR simulation result.
      logical, intent(in), optional :: print_obs !! Print simulated observations when true.
      logical :: show_observations

      show_observations = .false.
      if (present(print_obs)) show_observations = print_obs
      print '(a,i0)', 'Simulated observations: ', size(simulation%series)
      print '(a,i0)', 'Burn-in: ', simulation%burnin
      call display_mixar_model(simulation%model)
      if (show_observations) print '(*(f12.6,1x))', simulation%series
   end subroutine display_mixar_simulation

   subroutine display_mixar_exact_forecast(forecast)
      !! Display exact forecast means, standard deviations, and intervals.
      type(mixar_exact_forecast_t), intent(in) :: forecast !! Exact Gaussian MAR forecast.
      integer :: horizon

      print '(a,f7.3)', 'Prediction interval level: ', forecast%level
      print '(a)', ' horizon          mean            sd         lower         upper      components'
      do horizon = 1, size(forecast%mean)
         print '(i8,4f14.6,i16)', horizon, forecast%mean(horizon), &
            sqrt(max(0.0_dp, forecast%variance(horizon))), &
            forecast%lower(horizon), forecast%upper(horizon), &
            size(forecast%distribution(horizon)%probability)
      end do
   end subroutine display_mixar_exact_forecast

   subroutine display_mixar_forecast_sample(forecast)
      !! Display simulation forecast means, standard deviations, and intervals.
      type(mixar_forecast_sample_t), intent(in) :: forecast !! Simulated MAR forecast.
      integer :: horizon

      print '(a,i0)', 'Simulated paths: ', size(forecast%paths, 2)
      print '(a,f7.3)', 'Prediction interval level: ', forecast%level
      print '(a)', ' horizon          mean            sd         lower         upper'
      do horizon = 1, size(forecast%mean)
         print '(i8,4f14.6)', horizon, forecast%mean(horizon), &
            sqrt(max(0.0_dp, forecast%variance(horizon))), &
            forecast%lower(horizon), forecast%upper(horizon)
      end do
   end subroutine display_mixar_forecast_sample

   subroutine display_mixar_regression_model(model)
      !! Display MAR regression coefficients and innovation-model parameters.
      type(mixar_regression_model_t), intent(in) :: model !! MAR regression model.
      integer :: column

      if (model%component_specific) then
         print '(a)', 'Component-specific regression coefficients'
      else
         print '(a)', 'Shared regression coefficients'
      end if
      do column = 1, size(model%coefficient, 2)
         print '(a,i0,a,*(f12.6,1x))', 'coefficient column ', column, ': ', &
            model%coefficient(:, column)
      end do
      call display_mixar_model(model%mar)
   end subroutine display_mixar_regression_model

   subroutine display_mixar_regression_fit(fit)
      !! Display MAR regression convergence, criteria, and model parameters.
      type(mixar_regression_fit_t), intent(in) :: fit !! Fitted MAR regression.

      print '(a,l1)', 'Converged: ', fit%converged
      print '(a,i0)', 'Iterations: ', fit%iterations
      print '(a,f14.5)', 'Log likelihood: ', fit%log_likelihood
      print '(a,f14.5)', 'AIC: ', fit%aic
      print '(a,f14.5)', 'BIC: ', fit%bic
      call display_mixar_regression_model(fit%model)
   end subroutine display_mixar_regression_fit

   subroutine display_mixar_regression_simulation(simulation, print_obs)
      !! Display MAR regression simulation parameters and optionally observations.
      type(mixar_regression_simulation_t), intent(in) :: simulation !! MAR regression simulation.
      logical, intent(in), optional :: print_obs !! Print simulated observations when true.
      logical :: show_observations

      show_observations = .false.
      if (present(print_obs)) show_observations = print_obs
      print '(a,i0)', 'Simulated observations: ', size(simulation%series)
      call display_mixar_regression_model(simulation%model)
      if (show_observations) print '(*(f12.6,1x))', simulation%series
   end subroutine display_mixar_regression_simulation

   subroutine display_mixar_seasonal_model(model)
      !! Display seasonal MAR lag structure and component parameters.
      type(mixar_seasonal_model_t), intent(in) :: model !! Seasonal MAR model.
      integer :: component

      print '(a,i0)', 'Seasonal period: ', model%seasonal_period
      call display_mixar_model(model%mar)
      do component = 1, size(model%seasonal_order)
         print '(a,i0,a,i0)', 'component ', component, &
            ' seasonal order: ', model%seasonal_order(component)
         if (model%seasonal_order(component) > 0) print '(a,*(f10.5,1x))', &
            '  seasonal AR: ', &
            model%seasonal_ar(1:model%seasonal_order(component), component)
      end do
   end subroutine display_mixar_seasonal_model

   subroutine display_mixar_seasonal_fit(fit)
      !! Display seasonal MAR convergence, criteria, and parameters.
      type(mixar_seasonal_fit_t), intent(in) :: fit !! Fitted seasonal MAR result.

      print '(a,l1)', 'Converged: ', fit%converged
      print '(a,i0)', 'Iterations: ', fit%iterations
      print '(a,f14.5)', 'Log likelihood: ', fit%log_likelihood
      print '(a,f14.5)', 'AIC: ', fit%aic
      print '(a,f14.5)', 'BIC: ', fit%bic
      call display_mixar_seasonal_model(fit%model)
   end subroutine display_mixar_seasonal_fit

   subroutine display_mixar_seasonal_simulation(simulation, print_obs)
      !! Display seasonal MAR simulation parameters and optionally observations.
      type(mixar_seasonal_simulation_t), intent(in) :: simulation !! Seasonal MAR simulation.
      logical, intent(in), optional :: print_obs !! Print simulated observations when true.
      logical :: show_observations

      show_observations = .false.
      if (present(print_obs)) show_observations = print_obs
      print '(a,i0)', 'Simulated observations: ', size(simulation%series)
      print '(a,i0)', 'Burn-in: ', simulation%burnin
      call display_mixar_seasonal_model(simulation%model)
      if (show_observations) print '(*(f12.6,1x))', simulation%series
   end subroutine display_mixar_seasonal_simulation

   subroutine display_mixar_inference(inference)
      !! Display parameter estimates, standard errors, and confidence intervals.
      type(mixar_inference_t), intent(in) :: inference !! Observed-information inference result.
      integer :: parameter

      print '(a,f7.3)', 'Confidence level: ', inference%level
      print '(a)', 'parameter                                      estimate      std. error'// &
         '         lower         upper'
      do parameter = 1, size(inference%estimate)
         print '(a48,4f14.6)', trim(inference%parameter(parameter)), &
            inference%estimate(parameter), inference%standard_error(parameter), &
            inference%lower(parameter), inference%upper(parameter)
      end do
   end subroutine display_mixar_inference

   subroutine display_mixar_diagnostics(diagnostics)
      !! Display principal residual and classification diagnostic statistics.
      type(mixar_diagnostics_t), intent(in) :: diagnostics !! MAR diagnostic result.

      print '(a,f12.5,a,f10.6)', 'Residual Ljung-Box: ', &
         diagnostics%residual_test%statistic, ', p=', &
         diagnostics%residual_test%p_value
      print '(a,f12.5,a,f10.6)', 'Squared-residual Ljung-Box: ', &
         diagnostics%squared_residual_test%statistic, ', p=', &
         diagnostics%squared_residual_test%p_value
      print '(a,f12.6)', 'Mean classification entropy: ', diagnostics%mean_entropy
      print '(a,f14.5)', 'BIC: ', diagnostics%bic
   end subroutine display_mixar_diagnostics

   subroutine display_mixar_var_model(model)
      !! Display mixture VAR probabilities, orders, coefficients, and covariances.
      type(mixar_var_model_t), intent(in) :: model !! Mixture VAR parameters.
      integer :: component, lag

      print '(a,i0)', 'Mixture VAR components: ', size(model%probability)
      print '(a,i0)', 'Variables: ', size(model%intercept, 1)
      do component = 1, size(model%probability)
         print '(a,i0,a,f10.6,a,i0)', 'component ', component, &
            ': probability=', model%probability(component), &
            ', order=', model%order(component)
         print '(a,*(f12.6,1x))', '  intercept: ', model%intercept(:, component)
         do lag = 1, model%order(component)
            print '(a,i0,a)', '  AR lag ', lag, ':'
            print '(*(f12.6,1x))', model%ar(:, :, lag, component)
         end do
         print '(a)', '  innovation covariance:'
         print '(*(f12.6,1x))', model%innovation_covariance(:, :, component)
      end do
   end subroutine display_mixar_var_model

   subroutine display_mixar_var_fit(fit)
      !! Display mixture VAR EM convergence, criteria, and parameters.
      type(mixar_var_fit_t), intent(in) :: fit !! Fitted mixture VAR result.

      print '(a,l1)', 'Converged: ', fit%converged
      print '(a,i0)', 'Iterations: ', fit%iterations
      print '(a,f14.5)', 'Log likelihood: ', fit%log_likelihood
      print '(a,f14.5)', 'AIC: ', fit%aic
      print '(a,f14.5)', 'BIC: ', fit%bic
      call display_mixar_var_model(fit%model)
   end subroutine display_mixar_var_fit

   subroutine display_mixar_var_simulation(simulation, print_obs)
      !! Display mixture VAR simulation parameters and optionally observations.
      type(mixar_var_simulation_t), intent(in) :: simulation !! Mixture VAR simulation.
      logical, intent(in), optional :: print_obs !! Print simulated observations when true.
      logical :: show_observations

      show_observations = .false.
      if (present(print_obs)) show_observations = print_obs
      print '(a,i0)', 'Simulated observations: ', size(simulation%series, 1)
      print '(a,i0)', 'Burn-in: ', simulation%burnin
      call display_mixar_var_model(simulation%model)
      if (show_observations) print '(*(f12.6,1x))', simulation%series
   end subroutine display_mixar_var_simulation

   subroutine display_mixar_var_forecast(forecast)
      !! Display mixture VAR forecast means and marginal intervals.
      type(mixar_var_forecast_t), intent(in) :: forecast !! Mixture VAR simulated forecast.
      integer :: horizon

      print '(a,i0)', 'Simulated paths: ', size(forecast%paths, 3)
      print '(a,f7.3)', 'Prediction interval level: ', forecast%level
      do horizon = 1, size(forecast%mean, 1)
         print '(a,i0)', 'horizon ', horizon
         print '(a,*(f12.6,1x))', '  mean:  ', forecast%mean(horizon, :)
         print '(a,*(f12.6,1x))', '  lower: ', forecast%lower(horizon, :)
         print '(a,*(f12.6,1x))', '  upper: ', forecast%upper(horizon, :)
      end do
   end subroutine display_mixar_var_forecast

   subroutine display_mixar_var_diagnostics(diagnostics)
      !! Display mixture VAR whiteness and classification diagnostics.
      type(mixar_var_diagnostics_t), intent(in) :: diagnostics !! Mixture VAR diagnostics.

      print '(a,f12.6)', 'Mean classification entropy: ', diagnostics%mean_entropy
      print '(a,f12.5,a,f10.6)', 'Multivariate Q: ', &
         diagnostics%white_noise_test%multivariate_q_statistic, ', p=', &
         diagnostics%white_noise_test%multivariate_q_p_value
   end subroutine display_mixar_var_diagnostics

   subroutine display_mixar_bayesian_draws(draws)
      !! Display Bayesian MAR sample size, acceptance rates, and posterior means.
      type(mixar_bayesian_draws_t), intent(in) :: draws !! Bayesian MAR posterior draws.

      print '(a,i0)', 'Retained posterior draws: ', size(draws%probability, 1)
      print '(a,i0)', 'Burn-in: ', draws%burnin
      print '(a,*(f10.5,1x))', 'AR acceptance rates: ', draws%acceptance_rate
      call display_mixar_model(draws%posterior_mean_model)
   end subroutine display_mixar_bayesian_draws

   subroutine display_mixar_marginal_likelihood(result)
      !! Display the marginal-likelihood decomposition.
      type(mixar_marginal_likelihood_t), intent(in) :: result !! Marginal-likelihood estimate.

      print '(a,f16.6)', 'Log likelihood: ', result%log_likelihood
      print '(a,f16.6)', 'Log prior density: ', result%log_prior_density
      print '(a,f16.6)', 'Log posterior ordinate: ', result%log_posterior_ordinate
      print '(a,f16.6)', 'Log marginal likelihood: ', &
         result%log_marginal_likelihood
   end subroutine display_mixar_marginal_likelihood

   subroutine display_mixar_order_selection(selection)
      !! Display RJMCMC order frequencies, modal orders, and acceptance rate.
      type(mixar_order_selection_t), intent(in) :: selection !! Mixture-order selection result.
      integer :: state

      print '(a,i0)', 'Maximum component order: ', selection%maximum_order
      print '(a,f10.6)', 'Order-move acceptance rate: ', selection%acceptance_rate
      print '(a,*(i0,1x))', 'Modal component orders: ', selection%modal_order
      print '(a)', 'Posterior order probabilities:'
      do state = 1, size(selection%posterior_probability)
         print '(*(i0,1x),f12.6)', selection%unique_order(state, :), &
            selection%posterior_probability(state)
      end do
   end subroutine display_mixar_order_selection

   subroutine display_mixar_initialization(initialization)
      !! Display data-driven MAR initialization and stability contractions.
      type(mixar_initialization_t), intent(in) :: initialization !! Gaussian MAR initialization.

      print '(a,i0)', 'Stability contractions: ', &
         initialization%stability_contractions
      call display_mixar_model(initialization%model)
   end subroutine display_mixar_initialization

   subroutine display_mixar_multistart(result)
      !! Display multistart convergence, likelihoods, and selected fit.
      type(mixar_multistart_t), intent(in) :: result !! Gaussian MAR multistart result.
      integer :: start

      print '(a,i0)', 'Successful converged starts: ', result%successful
      print '(a,i0)', 'Best start: ', result%best
      do start = 1, size(result%fit)
         print '(a,i0,a,l1,a,i0,a,f14.5)', 'start ', start, &
            ': converged=', result%fit(start)%converged, &
            ', info=', result%fit(start)%info, &
            ', log likelihood=', result%fit(start)%log_likelihood
      end do
      if (result%best > 0) call display_mixar_fit(result%fit(result%best))
   end subroutine display_mixar_multistart

   subroutine display_mixar_moment(moment)
      !! Display an analytical moment and whether it exists.
      type(mixar_moment_t), intent(in) :: moment !! Analytical scalar moment.

      print '(a,i0)', 'Moment order: ', moment%order
      print '(a,l1)', 'Exists: ', moment%exists
      if (moment%exists) print '(a,es18.8)', 'Value: ', moment%value
   end subroutine display_mixar_moment

   subroutine display_mixar_stationary_moments(moments)
      !! Display stable MAR stationary marginal mean and variance.
      type(mixar_stationary_moments_t), intent(in) :: moments !! Stationary MAR moments.

      print '(a,f14.6)', 'Stationary mean: ', moments%mean
      print '(a,f14.6)', 'Stationary variance: ', moments%variance
   end subroutine display_mixar_stationary_moments

end module mixar_mod
