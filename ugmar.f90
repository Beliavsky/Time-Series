! SPDX-License-Identifier: GPL-3.0-only
! SPDX-FileComment: Algorithms translated from the R uGMAR package.
module ugmar_mod
   !! Univariate Gaussian and Student-t mixture autoregressive algorithms.
   use kind_mod, only: dp
   use random_mod, only: random_uniform, random_standard_normal
   use gmvarkit_mod, only: gmvarkit_model_t, gmvarkit_evaluation_t, &
      gmvarkit_simulation_t, gmvarkit_forecast_t, &
      gmvarkit_unconditional_moments_t, gmvarkit_fit_t, &
      gmvarkit_genetic_fit_t, gmvarkit_multistart_fit_t, &
      gmvarkit_constraints_t, gmvarkit_inference_t, &
      gmvarkit_likelihood_profile_t, gmvarkit_regime_conversion_t, &
      gmvarkit_hypothesis_test_t, &
      gmvarkit_companion_eigen_t, gmvarkit_quantile_residuals_t, &
      gmvarkit_residual_tests_t, gmvarkit_evaluate, gmvarkit_simulate, &
      gmvarkit_forecast, gmvarkit_unconditional_moments, &
      gmvarkit_estimate, gmvarkit_companion_eigenvalues, &
      gmvarkit_genetic_estimate, gmvarkit_multistart_estimate, &
      gmvarkit_estimate_constrained, gmvarkit_inference, &
      gmvarkit_wald_test, gmvarkit_likelihood_ratio, &
      gmvarkit_profile_likelihood, gmvarkit_convert_student_regimes, &
      gmvarkit_quantile_residuals, gmvarkit_quantile_residual_tests
   implicit none
   private

   type, public :: ugmar_model_t
      !! Parameters of a GMAR, StMAR, or G-StMAR model.
      real(dp), allocatable :: intercept(:)
      real(dp), allocatable :: ar(:, :)
      real(dp), allocatable :: innovation_variance(:)
      real(dp), allocatable :: weight(:)
      real(dp), allocatable :: degrees_of_freedom(:)
      integer :: gaussian_regimes = 0
      integer :: info = 0
   end type ugmar_model_t

   type, public :: ugmar_evaluation_t
      !! Conditional regime quantities and log likelihood.
      real(dp), allocatable :: mixing_weight(:, :)
      real(dp), allocatable :: regime_mean(:, :)
      real(dp), allocatable :: regime_variance(:, :)
      real(dp), allocatable :: conditional_mean(:)
      real(dp), allocatable :: conditional_variance(:)
      real(dp), allocatable :: log_likelihood_term(:)
      real(dp) :: log_likelihood = -huge(1.0_dp)
      integer :: info = 0
   end type ugmar_evaluation_t

   type, public :: ugmar_stationary_moments_t
      !! Regime-specific and mixture stationary moments through the AR order.
      real(dp), allocatable :: regime_mean(:)
      real(dp), allocatable :: regime_variance(:)
      real(dp), allocatable :: regime_autocovariance(:, :)
      real(dp) :: mean = 0.0_dp
      real(dp) :: variance = 0.0_dp
      real(dp), allocatable :: autocovariance(:)
      real(dp), allocatable :: autocorrelation(:)
      integer :: info = 0
   end type ugmar_stationary_moments_t

   type, public :: ugmar_simulation_t
      !! Simulated observations, regimes, and endogenous mixing weights.
      real(dp), allocatable :: series(:, :)
      integer, allocatable :: regime(:, :)
      real(dp), allocatable :: mixing_weight(:, :, :)
      integer :: info = 0
   end type ugmar_simulation_t

   type, public :: ugmar_forecast_t
      !! Monte Carlo forecasts and future mixing-weight summaries.
      real(dp), allocatable :: mean(:)
      real(dp), allocatable :: median(:)
      real(dp), allocatable :: quantile(:, :)
      real(dp), allocatable :: probability(:)
      real(dp), allocatable :: mixing_weight_mean(:, :)
      integer :: simulations = 0
      integer :: info = 0
   end type ugmar_forecast_t

   type, public :: ugmar_fit_t
      !! Maximum-likelihood model estimate and optimization diagnostics.
      type(ugmar_model_t) :: model
      real(dp) :: log_likelihood = -huge(1.0_dp)
      integer :: iterations = 0
      integer :: parameter_count = 0
      logical :: converged = .false.
      integer :: info = 0
   end type ugmar_fit_t

   type, public :: ugmar_genetic_fit_t
      !! Genetic-search path and locally refined best model.
      type(ugmar_fit_t) :: fit
      real(dp), allocatable :: best_objective(:)
      integer :: population_size = 0
      integer :: generations = 0
      integer :: evaluations = 0
      integer :: info = 0
   end type ugmar_genetic_fit_t

   type, public :: ugmar_multistart_fit_t
      !! Ranked local fits with duplicate-solution classifications.
      type(ugmar_fit_t), allocatable :: fit(:)
      integer, allocatable :: order(:)
      logical, allocatable :: successful(:)
      logical, allocatable :: distinct(:)
      integer, allocatable :: duplicate_of(:)
      integer :: best_index = 0
      integer :: successful_count = 0
      integer :: distinct_count = 0
      integer :: info = 0
   end type ugmar_multistart_fit_t

   type, public :: ugmar_constraints_t
      !! Linear AR, shared-mean, and fixed-weight restrictions.
      real(dp), allocatable :: ar_mapping(:, :)
      integer, allocatable :: mean_group(:)
      real(dp), allocatable :: fixed_weight(:)
      integer :: info = 0
   end type ugmar_constraints_t

   type, public :: ugmar_inference_t
      !! Numerical likelihood derivatives, covariance, and standard errors.
      real(dp), allocatable :: parameter(:)
      real(dp), allocatable :: gradient(:)
      real(dp), allocatable :: hessian(:, :)
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: standard_error(:)
      integer :: info = 0
   end type ugmar_inference_t

   type, public :: ugmar_likelihood_profile_t
      !! Fixed-coordinate likelihood profiles on transformed parameters.
      integer, allocatable :: parameter(:)
      real(dp), allocatable :: center(:)
      real(dp), allocatable :: value(:, :)
      real(dp), allocatable :: log_likelihood(:, :)
      logical, allocatable :: valid(:, :)
      integer :: info = 0
   end type ugmar_likelihood_profile_t

   type, public :: ugmar_regime_conversion_t
      !! Student-to-Gaussian regime conversion and permutation metadata.
      type(ugmar_fit_t) :: fit
      integer, allocatable :: old_to_new(:)
      integer, allocatable :: new_to_old(:)
      logical, allocatable :: converted(:)
      logical :: reestimated = .false.
      integer :: info = 0
   end type ugmar_regime_conversion_t

   type, public :: ugmar_ar_roots_t
      !! Characteristic-polynomial root moduli by regime.
      real(dp), allocatable :: modulus(:, :)
      integer :: info = 0
   end type ugmar_ar_roots_t

   type, public :: ugmar_quantile_residuals_t
      !! Probability-integral-transform residuals and conditional CDF values.
      real(dp), allocatable :: residual(:)
      real(dp), allocatable :: cdf(:)
      integer :: info = 0
   end type ugmar_quantile_residuals_t

   type, public :: ugmar_hypothesis_test_t
      !! Chi-square residual-test result.
      real(dp) :: statistic = 0.0_dp
      integer :: degrees_of_freedom = 0
      real(dp) :: p_value = 1.0_dp
      integer :: info = 0
   end type ugmar_hypothesis_test_t

   type, public :: ugmar_residual_tests_t
      !! Quantile-residual normality, correlation, and variance tests.
      type(ugmar_hypothesis_test_t) :: normality
      type(ugmar_hypothesis_test_t), allocatable :: autocorrelation(:)
      type(ugmar_hypothesis_test_t), allocatable :: heteroskedasticity(:)
      integer, allocatable :: autocorrelation_lag(:)
      integer, allocatable :: heteroskedasticity_lag(:)
      logical :: parameter_corrected = .false.
      integer :: info = 0
   end type ugmar_residual_tests_t

   public :: ugmar_model, ugmar_is_stationary
   public :: ugmar_random_model, ugmar_smart_model, ugmar_ar_roots
   public :: ugmar_regime_means, ugmar_model_from_regime_means
   public :: ugmar_evaluate, ugmar_log_likelihood, ugmar_mixing_weights
   public :: ugmar_stationary_moments, ugmar_estimate, ugmar_genetic_estimate
   public :: ugmar_multistart_estimate, ugmar_estimate_constrained
   public :: ugmar_inference, ugmar_wald_test, ugmar_likelihood_ratio
   public :: ugmar_profile_likelihood, ugmar_convert_student_regimes
   public :: ugmar_simulate, ugmar_forecast
   public :: ugmar_quantile_residuals, ugmar_quantile_residual_tests

contains

   pure function ugmar_model(intercept, ar, innovation_variance, weight, &
      gaussian_regimes, degrees_of_freedom) result(out)
      !! Construct and validate a stationary GMAR, StMAR, or G-StMAR model.
      real(dp), intent(in) :: intercept(:) !! Regime intercepts.
      real(dp), intent(in) :: ar(:, :) !! AR coefficients by lag and regime.
      real(dp), intent(in) :: innovation_variance(:) !! Positive regime innovation variances.
      real(dp), intent(in) :: weight(:) !! Positive regime mixing proportions summing to one.
      integer, intent(in), optional :: gaussian_regimes !! Number of leading Gaussian regimes.
      real(dp), intent(in), optional :: degrees_of_freedom(:) !! Student-t degrees of freedom by regime.
      type(ugmar_model_t) :: out
      integer :: regimes, gaussian_count

      regimes = size(intercept)
      gaussian_count = regimes
      if (present(gaussian_regimes)) gaussian_count = gaussian_regimes
      if (regimes < 1 .or. size(ar, 1) < 1 .or. size(ar, 2) /= regimes .or. &
         size(innovation_variance) /= regimes .or. size(weight) /= regimes) then
         out%info = 1
         return
      end if
      if (any(innovation_variance <= 0.0_dp) .or. any(weight <= 0.0_dp) .or. &
         abs(sum(weight) - 1.0_dp) > 1.0e-10_dp) then
         out%info = 2
         return
      end if
      if (gaussian_count < 0 .or. gaussian_count > regimes) then
         out%info = 3
         return
      end if
      allocate(out%degrees_of_freedom(regimes), source=0.0_dp)
      if (present(degrees_of_freedom)) then
         if (size(degrees_of_freedom) /= regimes) then
            out%info = 3
            return
         end if
         out%degrees_of_freedom = degrees_of_freedom
      end if
      if (gaussian_count < regimes) then
         if (.not. present(degrees_of_freedom) .or. &
            any(out%degrees_of_freedom(gaussian_count + 1:) <= 2.0_dp)) then
            out%info = 3
            return
         end if
      end if
      out%intercept = intercept
      out%ar = ar
      out%innovation_variance = innovation_variance
      out%weight = weight
      out%gaussian_regimes = gaussian_count
      if (.not. ugmar_is_stationary(out)) out%info = 4
   end function ugmar_model

   function ugmar_random_model(order, regimes, gaussian_regimes, &
      location_mean, location_standard_deviation, variance_scale, &
      minimum_degrees_of_freedom, maximum_degrees_of_freedom) result(out)
      !! Generate a stationary random model for global-search initialization.
      integer, intent(in) :: order !! Positive common autoregressive order.
      integer, intent(in) :: regimes !! Positive number of mixture regimes.
      integer, intent(in), optional :: gaussian_regimes !! Number of leading Gaussian regimes.
      real(dp), intent(in), optional :: location_mean !! Center of stationary regime means.
      real(dp), intent(in), optional :: location_standard_deviation !! Positive scale of regime means.
      real(dp), intent(in), optional :: variance_scale !! Positive scale of innovation variances.
      real(dp), intent(in), optional :: minimum_degrees_of_freedom !! Lower Student-t degrees bound.
      real(dp), intent(in), optional :: maximum_degrees_of_freedom !! Upper Student-t degrees bound.
      type(ugmar_model_t) :: out
      real(dp), allocatable :: intercept(:), ar(:, :), variance(:), weight(:)
      real(dp), allocatable :: degrees(:), mean(:)
      real(dp) :: mean_center, mean_scale, innovation_scale, minimum_df, maximum_df
      integer :: gaussian_count, regime, lag, attempt

      gaussian_count = regimes
      if (present(gaussian_regimes)) gaussian_count = gaussian_regimes
      mean_center = 0.0_dp
      if (present(location_mean)) mean_center = location_mean
      mean_scale = 1.0_dp
      if (present(location_standard_deviation)) &
         mean_scale = location_standard_deviation
      innovation_scale = 1.0_dp
      if (present(variance_scale)) innovation_scale = variance_scale
      minimum_df = 2.1_dp
      if (present(minimum_degrees_of_freedom)) &
         minimum_df = minimum_degrees_of_freedom
      maximum_df = 30.0_dp
      if (present(maximum_degrees_of_freedom)) &
         maximum_df = maximum_degrees_of_freedom
      if (order < 1 .or. regimes < 1 .or. gaussian_count < 0 .or. &
         gaussian_count > regimes .or. mean_scale <= 0.0_dp .or. &
         innovation_scale <= 0.0_dp .or. minimum_df <= 2.0_dp .or. &
         maximum_df <= minimum_df) then
         out%info = 1
         return
      end if
      allocate(intercept(regimes), ar(order, regimes), variance(regimes))
      allocate(weight(regimes), degrees(regimes), mean(regimes))
      do attempt = 1, 42
         do regime = 1, regimes
            mean(regime) = mean_center + mean_scale*random_standard_normal()
            variance(regime) = max(epsilon(1.0_dp), &
               abs(innovation_scale*random_standard_normal()))
            weight(regime) = max(epsilon(1.0_dp), random_uniform())
            do lag = 1, order
               ar(lag, regime) = (2.0_dp*random_uniform() - 1.0_dp)/ &
                  real(order, dp)
            end do
            intercept(regime) = mean(regime)*(1.0_dp - sum(ar(:, regime)))
         end do
         weight = weight/sum(weight)
         degrees = 0.0_dp
         do regime = gaussian_count + 1, regimes
            degrees(regime) = minimum_df + &
               (maximum_df - minimum_df)*random_uniform()
         end do
         out = ugmar_model(intercept, ar, variance, weight, gaussian_count, &
            degrees)
         if (out%info == 0) return
      end do
      out%info = 2
   end function ugmar_random_model

   function ugmar_smart_model(model, accuracy, random_regimes, location_mean, &
      location_standard_deviation, variance_scale) result(out)
      !! Perturb a fitted model locally, optionally replacing selected regimes.
      type(ugmar_model_t), intent(in) :: model !! Stationary model around which to mutate.
      real(dp), intent(in) :: accuracy !! Positive inverse mutation scale.
      integer, intent(in), optional :: random_regimes(:) !! Regimes replaced by random draws.
      real(dp), intent(in), optional :: location_mean !! Center for randomly replaced regime means.
      real(dp), intent(in), optional :: location_standard_deviation !! Scale for random regime means.
      real(dp), intent(in), optional :: variance_scale !! Scale for random innovation variances.
      type(ugmar_model_t) :: out
      type(ugmar_model_t) :: random_model
      real(dp), allocatable :: mean(:), weight(:), random_mean(:)
      integer :: regimes, regime, lag
      logical :: replace

      out = model
      if (model%info /= 0 .or. accuracy <= 0.0_dp) then
         out%info = 1
         return
      end if
      regimes = size(model%intercept)
      if (present(random_regimes)) then
         if (any(random_regimes < 1) .or. any(random_regimes > regimes)) then
            out%info = 1
            return
         end if
         random_model = ugmar_random_model(size(model%ar, 1), regimes, &
            model%gaussian_regimes, location_mean, &
            location_standard_deviation, variance_scale)
         if (random_model%info /= 0) then
            out%info = 2
            return
         end if
         random_mean = ugmar_regime_means(random_model)
      end if
      mean = ugmar_regime_means(model)
      weight = model%weight
      do regime = 1, regimes
         replace = .false.
         if (present(random_regimes)) replace = any(random_regimes == regime)
         if (replace) then
            mean(regime) = random_mean(regime)
            out%ar(:, regime) = random_model%ar(:, regime)
            out%innovation_variance(regime) = &
               random_model%innovation_variance(regime)
            out%degrees_of_freedom(regime) = &
               random_model%degrees_of_freedom(regime)
         else
            mean(regime) = mean(regime) + &
               abs(mean(regime)/accuracy)*random_standard_normal()
            do lag = 1, size(model%ar, 1)
               out%ar(lag, regime) = model%ar(lag, regime) + &
                  abs(model%ar(lag, regime)/accuracy)*random_standard_normal()
            end do
            out%innovation_variance(regime) = max(epsilon(1.0_dp), &
               abs(model%innovation_variance(regime) + &
               model%innovation_variance(regime)* &
               random_standard_normal()/accuracy))
            if (regime > model%gaussian_regimes) then
               out%degrees_of_freedom(regime) = max(2.1_dp, &
                  model%degrees_of_freedom(regime) + &
                  model%degrees_of_freedom(regime)* &
                  random_standard_normal()/accuracy)
            end if
         end if
         weight(regime) = abs(model%weight(regime) + &
            0.2_dp*random_standard_normal())
      end do
      weight = max(weight, epsilon(1.0_dp))
      out%weight = weight/sum(weight)
      out%intercept = mean*(1.0_dp - sum(out%ar, dim=1))
      out%info = 0
      if (.not. ugmar_is_stationary(out)) out = model
   end function ugmar_smart_model

   pure function ugmar_ar_roots(model) result(out)
      !! Return absolute characteristic-polynomial roots for every regime.
      type(ugmar_model_t), intent(in) :: model !! Model whose AR roots are calculated.
      type(ugmar_ar_roots_t) :: out
      type(gmvarkit_companion_eigen_t) :: eigen

      eigen = gmvarkit_companion_eigenvalues(to_gmvarkit(model))
      out%info = eigen%info
      if (eigen%info /= 0) return
      allocate(out%modulus(size(eigen%modulus, 1), size(eigen%modulus, 2)))
      where (eigen%modulus > tiny(1.0_dp))
         out%modulus = 1.0_dp/eigen%modulus
      elsewhere
         out%modulus = huge(1.0_dp)
      end where
   end function ugmar_ar_roots

   pure logical function ugmar_is_stationary(model, tolerance) result(stationary)
      !! Test whether every regime AR polynomial is stationary.
      type(ugmar_model_t), intent(in) :: model !! Model whose regime roots are tested.
      real(dp), intent(in), optional :: tolerance !! Distance below one treated as the stability boundary.
      type(gmvarkit_companion_eigen_t) :: eigen
      real(dp) :: selected_tolerance

      stationary = .false.
      selected_tolerance = 0.0_dp
      if (present(tolerance)) selected_tolerance = tolerance
      if (model%info /= 0 .and. model%info /= 4) return
      eigen = gmvarkit_companion_eigenvalues(to_gmvarkit(model), tolerance)
      if (eigen%info == 0) stationary = all(eigen%spectral_radius < &
         1.0_dp - selected_tolerance)
   end function ugmar_is_stationary

   pure function ugmar_regime_means(model) result(mean)
      !! Convert the stored intercepts to stationary regime means.
      type(ugmar_model_t), intent(in) :: model !! Intercept-parameterized model.
      real(dp), allocatable :: mean(:)

      if (model%info /= 0 .or. .not. allocated(model%intercept) .or. &
         .not. allocated(model%ar)) then
         allocate(mean(0))
         return
      end if
      mean = model%intercept/(1.0_dp - sum(model%ar, dim=1))
   end function ugmar_regime_means

   pure function ugmar_model_from_regime_means(template, mean) result(out)
      !! Replace a model's location parameters using stationary regime means.
      type(ugmar_model_t), intent(in) :: template !! Model supplying dynamics and distributions.
      real(dp), intent(in) :: mean(:) !! Stationary mean for each regime.
      type(ugmar_model_t) :: out

      out = template
      if (template%info /= 0 .or. .not. allocated(template%intercept) .or. &
         .not. allocated(template%ar) .or. size(mean) /= size(template%intercept)) then
         out%info = 1
         return
      end if
      out%intercept = mean*(1.0_dp - sum(template%ar, dim=1))
   end function ugmar_model_from_regime_means

   pure function ugmar_evaluate(series, model, conditional) result(out)
      !! Evaluate endogenous weights, conditional moments, and log likelihood.
      real(dp), intent(in) :: series(:) !! Observed univariate time series.
      type(ugmar_model_t), intent(in) :: model !! GMAR, StMAR, or G-StMAR model.
      logical, intent(in), optional :: conditional !! Omit the initial stationary-density term.
      type(ugmar_evaluation_t) :: out
      type(gmvarkit_evaluation_t) :: base
      real(dp), allocatable :: matrix(:, :)
      integer :: regime

      matrix = reshape(series, [size(series), 1])
      base = gmvarkit_evaluate(matrix, to_gmvarkit(model), conditional)
      out%info = base%info
      if (base%info /= 0) return
      out%mixing_weight = base%mixing_weight
      out%regime_mean = base%regime_mean(1, :, :)
      out%regime_variance = base%arch_scalar
      do regime = 1, size(model%intercept)
         out%regime_variance(:, regime) = out%regime_variance(:, regime)* &
            model%innovation_variance(regime)
      end do
      out%conditional_mean = base%conditional_mean(1, :)
      out%conditional_variance = base%conditional_covariance(1, 1, :)
      out%log_likelihood_term = base%log_likelihood_term
      out%log_likelihood = base%log_likelihood
   end function ugmar_evaluate

   pure real(dp) function ugmar_log_likelihood(series, model, conditional) &
      result(log_likelihood)
      !! Return the conditional or exact log likelihood.
      real(dp), intent(in) :: series(:) !! Observed univariate time series.
      type(ugmar_model_t), intent(in) :: model !! GMAR, StMAR, or G-StMAR model.
      logical, intent(in), optional :: conditional !! Omit the initial stationary-density term.
      type(ugmar_evaluation_t) :: evaluation

      evaluation = ugmar_evaluate(series, model, conditional)
      log_likelihood = evaluation%log_likelihood
   end function ugmar_log_likelihood

   pure function ugmar_mixing_weights(series, model) result(weight)
      !! Return observation-dependent regime probabilities for each usable time.
      real(dp), intent(in) :: series(:) !! Observed univariate time series.
      type(ugmar_model_t), intent(in) :: model !! GMAR, StMAR, or G-StMAR model.
      real(dp), allocatable :: weight(:, :)
      type(ugmar_evaluation_t) :: evaluation

      evaluation = ugmar_evaluate(series, model)
      if (evaluation%info == 0) then
         weight = evaluation%mixing_weight
      else
         allocate(weight(0, 0))
      end if
   end function ugmar_mixing_weights

   pure function ugmar_stationary_moments(model) result(out)
      !! Calculate regime and mixture stationary moments through lag p.
      type(ugmar_model_t), intent(in) :: model !! GMAR, StMAR, or G-StMAR model.
      type(ugmar_stationary_moments_t) :: out
      type(gmvarkit_unconditional_moments_t) :: base

      base = gmvarkit_unconditional_moments(to_gmvarkit(model))
      out%info = base%info
      if (base%info /= 0) return
      out%regime_mean = model%intercept/(1.0_dp - sum(model%ar, dim=1))
      out%regime_autocovariance = base%regime_autocovariance(1, 1, :, :)
      out%regime_variance = out%regime_autocovariance(1, :)
      out%mean = base%mean(1)
      out%variance = base%autocovariance(1, 1, 1)
      out%autocovariance = base%autocovariance(1, 1, :)
      out%autocorrelation = base%autocorrelation(1, 1, :)
   end function ugmar_stationary_moments

   pure function ugmar_estimate(series, initial_model, conditional, &
      max_iterations, tolerance) result(out)
      !! Refine a model by finite-difference BFGS maximum likelihood.
      real(dp), intent(in) :: series(:) !! Observed univariate time series.
      type(ugmar_model_t), intent(in) :: initial_model !! Valid starting model.
      logical, intent(in), optional :: conditional !! Optimize the conditional likelihood.
      integer, intent(in), optional :: max_iterations !! Maximum BFGS iterations.
      real(dp), intent(in), optional :: tolerance !! Gradient convergence tolerance.
      type(ugmar_fit_t) :: out
      type(gmvarkit_fit_t) :: base
      real(dp), allocatable :: matrix(:, :)

      matrix = reshape(series, [size(series), 1])
      base = gmvarkit_estimate(matrix, to_gmvarkit(initial_model), conditional, &
         max_iterations, tolerance)
      out%model = from_gmvarkit(base%model)
      out%log_likelihood = base%log_likelihood
      out%iterations = base%iterations
      out%parameter_count = base%parameter_count
      out%converged = base%converged
      out%info = base%info
   end function ugmar_estimate

   function ugmar_genetic_estimate(series, initial_model, conditional, &
      population_size, generations, mutation_scale, local_iterations) result(out)
      !! Search globally with an elitist genetic algorithm and local refinement.
      real(dp), intent(in) :: series(:) !! Observed univariate time series.
      type(ugmar_model_t), intent(in) :: initial_model !! Center of the initial population.
      logical, intent(in), optional :: conditional !! Optimize the conditional likelihood.
      integer, intent(in), optional :: population_size !! Candidates retained per generation.
      integer, intent(in), optional :: generations !! Number of genetic generations.
      real(dp), intent(in), optional :: mutation_scale !! Initial transformed-space mutation scale.
      integer, intent(in), optional :: local_iterations !! Final BFGS iteration limit.
      type(ugmar_genetic_fit_t) :: out
      type(gmvarkit_genetic_fit_t) :: base
      real(dp), allocatable :: matrix(:, :)

      matrix = reshape(series, [size(series), 1])
      base = gmvarkit_genetic_estimate(matrix, to_gmvarkit(initial_model), &
         conditional, population_size, generations, mutation_scale, &
         local_iterations)
      out%fit = copy_fit(base%fit)
      if (allocated(base%best_objective)) out%best_objective = base%best_objective
      out%population_size = base%population_size
      out%generations = base%generations
      out%evaluations = base%evaluations
      out%info = base%info
   end function ugmar_genetic_estimate

   pure function ugmar_multistart_estimate(series, initial_model, conditional, &
      max_iterations, tolerance, likelihood_tolerance, parameter_tolerance) &
      result(out)
      !! Fit, rank, and deduplicate multiple local starting models.
      real(dp), intent(in) :: series(:) !! Observed univariate time series.
      type(ugmar_model_t), intent(in) :: initial_model(:) !! Models sharing one specification.
      logical, intent(in), optional :: conditional !! Optimize the conditional likelihood.
      integer, intent(in), optional :: max_iterations !! Maximum BFGS iterations per start.
      real(dp), intent(in), optional :: tolerance !! Gradient convergence tolerance.
      real(dp), intent(in), optional :: likelihood_tolerance !! Relative likelihood duplicate tolerance.
      real(dp), intent(in), optional :: parameter_tolerance !! Relative parameter duplicate tolerance.
      type(ugmar_multistart_fit_t) :: out
      type(gmvarkit_multistart_fit_t) :: base
      type(gmvarkit_model_t), allocatable :: starts(:)
      real(dp), allocatable :: matrix(:, :)
      integer :: start

      allocate(starts(size(initial_model)))
      do start = 1, size(initial_model)
         starts(start) = to_gmvarkit(initial_model(start))
      end do
      matrix = reshape(series, [size(series), 1])
      base = gmvarkit_multistart_estimate(matrix, starts, conditional, &
         max_iterations, tolerance, likelihood_tolerance, parameter_tolerance)
      out%info = base%info
      if (.not. allocated(base%fit)) return
      allocate(out%fit(size(base%fit)))
      do start = 1, size(base%fit)
         out%fit(start) = copy_fit(base%fit(start))
      end do
      out%order = base%order
      out%successful = base%successful
      out%distinct = base%distinct
      out%duplicate_of = base%duplicate_of
      out%best_index = base%best_index
      out%successful_count = base%successful_count
      out%distinct_count = base%distinct_count
   end function ugmar_multistart_estimate

   pure function ugmar_estimate_constrained(series, initial_model, constraints, &
      conditional, max_iterations, tolerance) result(out)
      !! Estimate under linear AR, shared-mean, and fixed-weight restrictions.
      real(dp), intent(in) :: series(:) !! Observed univariate time series.
      type(ugmar_model_t), intent(in) :: initial_model !! Valid starting model.
      type(ugmar_constraints_t), intent(in) :: constraints !! Parameter restrictions.
      logical, intent(in), optional :: conditional !! Optimize the conditional likelihood.
      integer, intent(in), optional :: max_iterations !! Maximum BFGS iterations.
      real(dp), intent(in), optional :: tolerance !! Gradient convergence tolerance.
      type(ugmar_fit_t) :: out
      type(gmvarkit_fit_t) :: base
      real(dp), allocatable :: matrix(:, :)

      matrix = reshape(series, [size(series), 1])
      base = gmvarkit_estimate_constrained(matrix, to_gmvarkit(initial_model), &
         to_gmvarkit_constraints(constraints), conditional, max_iterations, &
         tolerance)
      out = copy_fit(base)
   end function ugmar_estimate_constrained

   pure function ugmar_inference(series, fit, conditional, difference_step) &
      result(out)
      !! Calculate numerical gradient, Hessian covariance, and standard errors.
      real(dp), intent(in) :: series(:) !! Observed univariate time series.
      type(ugmar_fit_t), intent(in) :: fit !! Fitted unrestricted model.
      logical, intent(in), optional :: conditional !! Use the conditional likelihood.
      real(dp), intent(in), optional :: difference_step !! Relative finite-difference step.
      type(ugmar_inference_t) :: out
      type(gmvarkit_inference_t) :: base
      real(dp), allocatable :: matrix(:, :)

      matrix = reshape(series, [size(series), 1])
      base = gmvarkit_inference(matrix, to_gmvarkit_fit(fit), conditional, &
         difference_step)
      out = copy_inference(base)
   end function ugmar_inference

   pure function ugmar_wald_test(inference, restriction, null_value) result(out)
      !! Test linear restrictions on transformed model parameters.
      type(ugmar_inference_t), intent(in) :: inference !! Numerical parameter inference.
      real(dp), intent(in) :: restriction(:, :) !! Full-row-rank restriction matrix.
      real(dp), intent(in) :: null_value(:) !! Values imposed under the null hypothesis.
      type(ugmar_hypothesis_test_t) :: out

      out = copy_test(gmvarkit_wald_test(to_gmvarkit_inference(inference), &
         restriction, null_value))
   end function ugmar_wald_test

   pure function ugmar_likelihood_ratio(unrestricted, restricted) result(out)
      !! Compare nested models with a likelihood-ratio test.
      type(ugmar_fit_t), intent(in) :: unrestricted !! Freely estimated model.
      type(ugmar_fit_t), intent(in) :: restricted !! Nested restricted model.
      type(ugmar_hypothesis_test_t) :: out

      out = copy_test(gmvarkit_likelihood_ratio( &
         to_gmvarkit_fit(unrestricted), to_gmvarkit_fit(restricted)))
   end function ugmar_likelihood_ratio

   pure function ugmar_profile_likelihood(series, model, parameter, scale, &
      points, conditional) result(out)
      !! Evaluate fixed-coordinate likelihood profiles around a model.
      real(dp), intent(in) :: series(:) !! Observed univariate time series.
      type(ugmar_model_t), intent(in) :: model !! Model at the profile center.
      integer, intent(in), optional :: parameter(:) !! One-based transformed coordinates.
      real(dp), intent(in), optional :: scale !! Relative profile half-width.
      integer, intent(in), optional :: points !! Equally spaced values per profile.
      logical, intent(in), optional :: conditional !! Evaluate conditional likelihoods.
      type(ugmar_likelihood_profile_t) :: out
      type(gmvarkit_likelihood_profile_t) :: base
      real(dp), allocatable :: matrix(:, :)

      matrix = reshape(series, [size(series), 1])
      base = gmvarkit_profile_likelihood(matrix, to_gmvarkit(model), &
         parameter, scale, points, conditional)
      out%info = base%info
      if (allocated(base%parameter)) out%parameter = base%parameter
      if (allocated(base%center)) out%center = base%center
      if (allocated(base%value)) out%value = base%value
      if (allocated(base%log_likelihood)) &
         out%log_likelihood = base%log_likelihood
      if (allocated(base%valid)) out%valid = base%valid
   end function ugmar_profile_likelihood

   pure function ugmar_convert_student_regimes(model, &
      maximum_degrees_of_freedom, series, estimate, conditional, &
      max_iterations, tolerance) result(out)
      !! Replace high-degree Student regimes by Gaussian regimes and reorder.
      type(ugmar_model_t), intent(in) :: model !! Source StMAR or G-StMAR model.
      real(dp), intent(in), optional :: maximum_degrees_of_freedom !! Strict conversion threshold.
      real(dp), intent(in), optional :: series(:) !! Data used for evaluation or re-estimation.
      logical, intent(in), optional :: estimate !! Re-estimate the converted model.
      logical, intent(in), optional :: conditional !! Use conditional likelihood when data are supplied.
      integer, intent(in), optional :: max_iterations !! Maximum optional BFGS iterations.
      real(dp), intent(in), optional :: tolerance !! Optional gradient tolerance.
      type(ugmar_regime_conversion_t) :: out
      type(gmvarkit_regime_conversion_t) :: base
      real(dp), allocatable :: matrix(:, :)

      if (present(series)) then
         matrix = reshape(series, [size(series), 1])
         base = gmvarkit_convert_student_regimes(to_gmvarkit(model), &
            maximum_degrees_of_freedom=maximum_degrees_of_freedom, &
            series=matrix, estimate=estimate, conditional=conditional, &
            max_iterations=max_iterations, tolerance=tolerance)
      else
         base = gmvarkit_convert_student_regimes(to_gmvarkit(model), &
            maximum_degrees_of_freedom=maximum_degrees_of_freedom, &
            estimate=estimate, conditional=conditional, &
            max_iterations=max_iterations, tolerance=tolerance)
      end if
      out%fit = copy_fit(base%fit)
      if (allocated(base%old_to_new)) out%old_to_new = base%old_to_new
      if (allocated(base%new_to_old)) out%new_to_old = base%new_to_old
      if (allocated(base%converted)) out%converted = base%converted
      out%reestimated = base%reestimated
      out%info = base%info
   end function ugmar_convert_student_regimes

   function ugmar_simulate(model, initial_values, observations, paths) result(out)
      !! Simulate paths with endogenous regime probabilities.
      type(ugmar_model_t), intent(in) :: model !! GMAR, StMAR, or G-StMAR model.
      real(dp), intent(in) :: initial_values(:) !! Initial history in chronological order.
      integer, intent(in) :: observations !! Number of observations generated per path.
      integer, intent(in), optional :: paths !! Number of independent paths.
      type(ugmar_simulation_t) :: out
      type(gmvarkit_simulation_t) :: base
      real(dp), allocatable :: history(:, :)

      history = reshape(initial_values, [size(initial_values), 1])
      base = gmvarkit_simulate(to_gmvarkit(model), history, observations, paths)
      out%info = base%info
      if (base%info /= 0) return
      out%series = base%series(:, 1, :)
      out%regime = base%regime
      out%mixing_weight = base%mixing_weight
   end function ugmar_simulate

   function ugmar_forecast(model, initial_values, horizon, simulations, &
      probabilities) result(out)
      !! Summarize recursive simulation draws as predictive forecasts.
      type(ugmar_model_t), intent(in) :: model !! GMAR, StMAR, or G-StMAR model.
      real(dp), intent(in) :: initial_values(:) !! Initial history in chronological order.
      integer, intent(in) :: horizon !! Positive forecast horizon.
      integer, intent(in), optional :: simulations !! Number of Monte Carlo paths.
      real(dp), intent(in), optional :: probabilities(:) !! Empirical quantile probabilities.
      type(ugmar_forecast_t) :: out
      type(gmvarkit_forecast_t) :: base
      real(dp), allocatable :: history(:, :)

      history = reshape(initial_values, [size(initial_values), 1])
      base = gmvarkit_forecast(to_gmvarkit(model), history, horizon, &
         simulations, probabilities)
      out%info = base%info
      if (base%info /= 0) return
      out%mean = base%mean(:, 1)
      out%median = base%median(:, 1)
      out%quantile = base%quantile(:, :, 1)
      out%probability = base%probability
      out%mixing_weight_mean = base%mixing_weight_mean
      out%simulations = base%simulations
   end function ugmar_forecast

   pure function ugmar_quantile_residuals(series, model) result(out)
      !! Calculate Gaussian quantile residuals from conditional mixture CDFs.
      real(dp), intent(in) :: series(:) !! Observed univariate time series.
      type(ugmar_model_t), intent(in) :: model !! Fitted GMAR, StMAR, or G-StMAR model.
      type(ugmar_quantile_residuals_t) :: out
      type(gmvarkit_quantile_residuals_t) :: base
      real(dp), allocatable :: matrix(:, :)

      matrix = reshape(series, [size(series), 1])
      base = gmvarkit_quantile_residuals(matrix, to_gmvarkit(model))
      out%info = base%info
      if (base%info /= 0) return
      out%residual = base%residual(:, 1)
      out%cdf = base%cdf(:, 1)
   end function ugmar_quantile_residuals

   pure function ugmar_quantile_residual_tests(residuals, &
      autocorrelation_lags, heteroskedasticity_lags, series, model, &
      difference_step, conditional) result(out)
      !! Test quantile-residual normality, correlation, and heteroskedasticity.
      real(dp), intent(in) :: residuals(:) !! Univariate quantile residuals.
      integer, intent(in) :: autocorrelation_lags(:) !! Maximum lags for correlation tests.
      integer, intent(in) :: heteroskedasticity_lags(:) !! Maximum lags for squared-residual tests.
      real(dp), intent(in), optional :: series(:) !! Data used for parameter-estimation correction.
      type(ugmar_model_t), intent(in), optional :: model !! Model used for parameter-estimation correction.
      real(dp), intent(in), optional :: difference_step !! Relative numerical derivative step.
      logical, intent(in), optional :: conditional !! Use conditional likelihood scores in correction.
      type(ugmar_residual_tests_t) :: out
      type(gmvarkit_residual_tests_t) :: base
      real(dp), allocatable :: residual_matrix(:, :), series_matrix(:, :)
      type(gmvarkit_model_t) :: base_model
      integer :: test

      residual_matrix = reshape(residuals, [size(residuals), 1])
      if (present(series) .neqv. present(model)) then
         out%info = 1
         return
      end if
      if (present(series) .and. present(model)) then
         series_matrix = reshape(series, [size(series), 1])
         base_model = to_gmvarkit(model)
         base = gmvarkit_quantile_residual_tests(residual_matrix, &
            autocorrelation_lags, heteroskedasticity_lags, series_matrix, &
            base_model, difference_step, conditional)
      else
         base = gmvarkit_quantile_residual_tests(residual_matrix, &
            autocorrelation_lags, heteroskedasticity_lags)
      end if
      out%info = base%info
      if (base%info /= 0) return
      out%normality = copy_test(base%normality)
      out%autocorrelation_lag = base%autocorrelation_lag
      out%heteroskedasticity_lag = base%heteroskedasticity_lag
      allocate(out%autocorrelation(size(base%autocorrelation)))
      allocate(out%heteroskedasticity(size(base%heteroskedasticity)))
      do test = 1, size(base%autocorrelation)
         out%autocorrelation(test) = copy_test(base%autocorrelation(test))
      end do
      do test = 1, size(base%heteroskedasticity)
         out%heteroskedasticity(test) = &
            copy_test(base%heteroskedasticity(test))
      end do
      out%parameter_corrected = base%parameter_corrected
   end function ugmar_quantile_residual_tests

   pure function to_gmvarkit(model) result(out)
      !! Expand a univariate uGMAR model into singleton VAR dimensions.
      type(ugmar_model_t), intent(in) :: model !! Univariate model parameters.
      type(gmvarkit_model_t) :: out
      integer :: lags, regimes, regime, lag

      if (.not. allocated(model%intercept) .or. .not. allocated(model%ar) .or. &
         .not. allocated(model%innovation_variance) .or. &
         .not. allocated(model%weight) .or. &
         .not. allocated(model%degrees_of_freedom)) then
         out%info = 1
         return
      end if
      lags = size(model%ar, 1)
      regimes = size(model%intercept)
      allocate(out%intercept(1, regimes), out%ar(1, 1, lags, regimes))
      allocate(out%covariance(1, 1, regimes))
      out%intercept(1, :) = model%intercept
      do regime = 1, regimes
         do lag = 1, lags
            out%ar(1, 1, lag, regime) = model%ar(lag, regime)
         end do
         out%covariance(1, 1, regime) = model%innovation_variance(regime)
      end do
      out%weight = model%weight
      out%degrees_of_freedom = model%degrees_of_freedom
      out%gaussian_regimes = model%gaussian_regimes
      out%info = merge(model%info, 0, model%info /= 4)
   end function to_gmvarkit

   pure function from_gmvarkit(model) result(out)
      !! Collapse a singleton mixture VAR model to uGMAR dimensions.
      type(gmvarkit_model_t), intent(in) :: model !! Singleton-variable mixture VAR model.
      type(ugmar_model_t) :: out
      integer :: lags, regimes, regime, lag

      if (.not. allocated(model%intercept) .or. size(model%intercept, 1) /= 1) then
         out%info = 1
         return
      end if
      lags = size(model%ar, 3)
      regimes = size(model%intercept, 2)
      allocate(out%intercept(regimes), out%ar(lags, regimes))
      allocate(out%innovation_variance(regimes))
      out%intercept = model%intercept(1, :)
      do regime = 1, regimes
         do lag = 1, lags
            out%ar(lag, regime) = model%ar(1, 1, lag, regime)
         end do
         out%innovation_variance(regime) = model%covariance(1, 1, regime)
      end do
      out%weight = model%weight
      out%degrees_of_freedom = model%degrees_of_freedom
      out%gaussian_regimes = model%gaussian_regimes
      out%info = model%info
   end function from_gmvarkit

   pure function copy_fit(source) result(out)
      !! Convert a shared-engine fit to the uGMAR fit type.
      type(gmvarkit_fit_t), intent(in) :: source !! Shared-engine model fit.
      type(ugmar_fit_t) :: out

      out%model = from_gmvarkit(source%model)
      out%log_likelihood = source%log_likelihood
      out%iterations = source%iterations
      out%parameter_count = source%parameter_count
      out%converged = source%converged
      out%info = source%info
   end function copy_fit

   pure function to_gmvarkit_fit(source) result(out)
      !! Convert a uGMAR fit to the shared-engine fit type.
      type(ugmar_fit_t), intent(in) :: source !! Univariate model fit.
      type(gmvarkit_fit_t) :: out

      out%model = to_gmvarkit(source%model)
      out%log_likelihood = source%log_likelihood
      out%iterations = source%iterations
      out%parameter_count = source%parameter_count
      out%converged = source%converged
      out%info = source%info
   end function to_gmvarkit_fit

   pure function to_gmvarkit_constraints(source) result(out)
      !! Convert univariate restrictions to the shared-engine representation.
      type(ugmar_constraints_t), intent(in) :: source !! Univariate model restrictions.
      type(gmvarkit_constraints_t) :: out

      if (allocated(source%ar_mapping)) out%ar_mapping = source%ar_mapping
      if (allocated(source%mean_group)) out%mean_group = source%mean_group
      if (allocated(source%fixed_weight)) out%fixed_weight = source%fixed_weight
      out%info = source%info
   end function to_gmvarkit_constraints

   pure function copy_inference(source) result(out)
      !! Convert shared-engine numerical inference to the uGMAR type.
      type(gmvarkit_inference_t), intent(in) :: source !! Shared-engine inference result.
      type(ugmar_inference_t) :: out

      if (allocated(source%parameter)) out%parameter = source%parameter
      if (allocated(source%gradient)) out%gradient = source%gradient
      if (allocated(source%hessian)) out%hessian = source%hessian
      if (allocated(source%covariance)) out%covariance = source%covariance
      if (allocated(source%standard_error)) &
         out%standard_error = source%standard_error
      out%info = source%info
   end function copy_inference

   pure function to_gmvarkit_inference(source) result(out)
      !! Convert uGMAR numerical inference to the shared-engine type.
      type(ugmar_inference_t), intent(in) :: source !! Univariate inference result.
      type(gmvarkit_inference_t) :: out

      if (allocated(source%parameter)) out%parameter = source%parameter
      if (allocated(source%gradient)) out%gradient = source%gradient
      if (allocated(source%hessian)) out%hessian = source%hessian
      if (allocated(source%covariance)) out%covariance = source%covariance
      if (allocated(source%standard_error)) &
         out%standard_error = source%standard_error
      out%info = source%info
   end function to_gmvarkit_inference

   pure elemental function copy_test(source) result(out)
      !! Convert a shared-engine hypothesis test to the uGMAR result type.
      type(gmvarkit_hypothesis_test_t), intent(in) :: source !! Shared-engine test result.
      type(ugmar_hypothesis_test_t) :: out

      out%statistic = source%statistic
      out%degrees_of_freedom = source%degrees_of_freedom
      out%p_value = source%p_value
      out%info = source%info
   end function copy_test

end module ugmar_mod
